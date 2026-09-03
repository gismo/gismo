/** @file composed_integration_test.cpp

    @brief Quadrature study for composed maps: finest-mesh vs super-mesh, and refinement rate

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gismo.h>
#include <iomanip>

using namespace gismo;
typedef real_t T;

/*
    Experiment. The physical domain of every geometry below is the unit
    square: S is built from Greville-point control points displaced by
    A*(sin(2*pi*x)*sin(pi*y), sin(pi*x)*sin(2*pi*y)), a perturbation that
    vanishes at x in {0,1} and y in {0,1}, so the boundary is unchanged
    while the p-th derivative of S jumps across every interior knot line.
    Composing S with sigma therefore reproduces, on the unit square, the
    mechanism that makes quadrature over a composed map G = S o sigma
    harder than quadrature over S alone: sigma pulls the preimages of S's
    knot lines into curves that a rule built only from sigma's own mesh
    does not resolve.

    sigma is one of three deformations of the identity on [0,1]^2:
      mode 0 (identity)      : sigma(u,v) = (u,v)
      mode 1 (direction-wise) : sigma(u,v) = (u + delta*sin(pi*u), v + delta/2*sin(2*pi*v))
      mode 2 (coupling)      : sigma(u,v) = (u + delta*sin(pi*u)*sin(pi*v), v + delta*sin(2*pi*u)*sin(pi*v))
    all built on a tensor B-spline basis of degree degreeMap with
    numKnotsMap interior knots per direction; every component vanishes on
    the boundary, so sigma maps [0,1]^2 onto itself.

    The integrand is phi(x,y) = (1-x)^2 + (y-x)^2; its exact integral over
    the unit square is 1/2, independent of S and sigma, so every relative
    error below is |value - 1/2| / (1/2) computed against this closed form,
    never against a numerically converged surrogate.

    Two quadrature strategies are compared, both with the library default
    node count (quA=1, quB=1, i.e. p*d+1 Gauss nodes per direction, where
    p*d is the degree of the integration basis):
      finest mesh : the basis (S's or sigma's) with more elements per
                    direction, degree-raised to p*d -- the rule a user gets
                    by picking whichever mesh looks finer, without merging
                    the other map's knot lines in;
      super mesh  : S's basis with sigma's interior knots inserted and the
                    degree raised to p*d (gsAdaptiveParametrization::
                    makeIntegrationBasis), optionally dyadically refined.

    Panel (a) (p = --degree, coupling sigma, delta = --delta) sweeps the
    number of analysis elements per direction (numEL = 1..6) and reports
    the relative error and quadrature-point count of both strategies;
    columns numEL, err_fine, err_super, nqp_fine, nqp_super.

    Panel (b) (numKnots = 3, delta = --delta, p in {2,3}) sweeps the
    dyadic refinement level r = 0..maxRefine of the super-mesh for the
    three sigma modes; columns level, h, err_id_p2, err_dir_p2,
    err_coup_p2, nqp_p2, err_id_p3, err_dir_p3, err_coup_p3, nqp_p3, with
    h the smallest super-mesh cell width 1/(12*2^r) (knots of S at
    multiples of 1/4, of sigma at multiples of 1/3).

    Theory and the derivation of the O(h^p) behaviour:
    INTEGRATION_DERIVATIONS.md, Sections 2-7 (repository root).
*/

// sigma: degree d, nM interior knots per direction.
//  mode 0: identity (Greville coefficients)
//  mode 1: direction-wise, sigma(u,v) = (a(u), b(v))
//  mode 2: coupling, each component depends on both u and v
// The perturbations vanish on the boundary of [0,1]^2 so sigma maps the
// boundary onto itself, and the coefficients stay inside [0,1]^2.
gsTensorBSpline<2,T> makeSigma(int d, int nM, int mode, T delta)
{
    gsKnotVector<T> kv(0,1,nM,d+1);
    gsTensorBSplineBasis<2,T> basis(kv,kv);
    gsMatrix<T> g = basis.anchors();
    gsMatrix<T> coefs(g.cols(),2);
    for (index_t i=0;i<g.cols();++i)
    {
        const T x = g(0,i), y = g(1,i);
        T sx = x, sy = y;
        if (mode==1)
        {
            sx = x + delta*std::sin(M_PI*x);
            sy = y + 0.5*delta*std::sin(2*M_PI*y);
        }
        else if (mode==2)
        {
            sx = x + delta*std::sin(M_PI*x)*std::sin(M_PI*y);
            sy = y + delta*std::sin(2*M_PI*x)*std::sin(M_PI*y);
        }
        coefs(i,0) = sx; coefs(i,1) = sy;
    }
    return gsTensorBSpline<2,T>(basis,coefs);
}

// Geometry S: degree p, nQ interior knots per direction, control points =
// Greville points displaced by A*(sin(2*pi*x)*sin(pi*y), sin(pi*x)*sin(2*pi*y)).
// The displacement vanishes at x in {0,1} and y in {0,1}, so the boundary
// controls are unchanged and the physical domain is the unit square.
gsTensorBSpline<2,T> makeS(int p, int nQ, T A)
{
    gsKnotVector<T> kv(0,1,nQ,p+1);
    gsTensorBSplineBasis<2,T> basis(kv,kv);
    gsMatrix<T> coefs = basis.anchors().transpose();
    for (index_t i=0;i<coefs.rows();++i)
    {
        const T x = coefs(i,0), y = coefs(i,1);
        coefs(i,0) += A*std::sin(2*M_PI*x)*std::sin(M_PI*y);
        coefs(i,1) += A*std::sin(M_PI*x)*std::sin(2*M_PI*y);
    }
    return gsTensorBSpline<2,T>(basis,coefs);
}

// Super-mesh basis exactly as gsAdaptiveParametrization::makeIntegrationBasis
// (src/gsAssembler/gsAdaptiveParametrization.hpp:2145-2210): S's basis,
// sigma's interior knots inserted, degree raised to p*d. Then `refine`
// dyadic uniform refinements of the result.
gsTensorBSplineBasis<2,T> superMesh(const gsTensorBSplineBasis<2,T> & bS,
                                    const gsTensorBSplineBasis<2,T> & bM,
                                    int refine)
{
    gsTensorBSplineBasis<2,T> ib(bS);
    for (int dir=0;dir<2;++dir)
    {
        for (auto it = std::next(bM.knots(dir).ubegin());
                  it != std::prev(bM.knots(dir).uend()); ++it)
            if (!ib.knots(dir).has(*it)) ib.insertKnot(*it,dir);
        const index_t target = bS.degree(dir)*bM.degree(dir);
        ib.degreeIncrease(target-ib.degree(dir),dir);
    }
    for (int r=0;r<refine;++r) ib.uniformRefine();
    return ib;
}

// Finest-mesh basis: whichever of bS, bM has more elements (in total, both
// directions), copied and degree-raised to p*d per direction, WITHOUT
// inserting the other map's interior knots. This is the "pick the finer
// mesh" strategy that leaves the other map's preimage curves unresolved.
gsTensorBSplineBasis<2,T> finestMesh(const gsTensorBSplineBasis<2,T> & bS,
                                     const gsTensorBSplineBasis<2,T> & bM)
{
    gsTensorBSplineBasis<2,T> ib = (bS.numElements() >= bM.numElements()) ? bS : bM;
    for (int dir=0;dir<2;++dir)
    {
        const index_t target = bS.degree(dir)*bM.degree(dir);
        ib.degreeIncrease(target-ib.degree(dir),dir);
    }
    return ib;
}

T minDetJ(const gsGeometry<T> & sig)
{
    gsVector<T> lo = gsVector<T>::Zero(2), hi = gsVector<T>::Ones(2);
    gsMatrix<T> pts = uniformPointGrid<T>(lo,hi,101*101);
    gsMatrix<T> der;
    sig.deriv_into(pts,der);
    T m = std::numeric_limits<T>::max();
    for (index_t c=0;c<der.cols();++c)
        m = std::min(m, der(0,c)*der(3,c)-der(1,c)*der(2,c));
    return m;
}

// Integrates int_[0,1]^2 f(map(.)) |J_map| over the integration basis ib,
// with SameElement=false since quadrature points of ib generally cross
// element boundaries of map's own basis (gsExprEvaluator::defaultOptions(),
// src/gsAssembler/gsExprEvaluator.h:79-89, does not register SameElement,
// so it must be added before any other option is touched).
T integrate(const gsFunctionSet<T> & map, const gsFunction<T> & f,
            const gsBasis<T> & ib, T quA, index_t quB, index_t & nQP)
{
    gsExprEvaluator<T> ev;
    gsOptionList o;
    o.addReal("quA","",quA);
    o.addInt ("quB","",quB);
    o.addSwitch("SameElement","",false);
    ev.options().update(o, gsOptionList::addIfUnknown);
    auto G = ev.getMap(map);
    auto F = ev.getVariable(f,G);
    gsMultiBasis<T> mb(ib);
    ev.setIntegrationElements(mb);
    nQP = gsQuadrature::getAllNodes(ib, ev.options()).cols();
    return ev.integral(F*meas(G));
}

int main(int argc, char *argv[])
{
    index_t degree      = 3;
    index_t degreeMap    = 2;
    index_t numKnotsMap  = 2;
    T amplitude          = 0.05;
    T delta              = 0.3;
    index_t maxRefine    = 4;
    std::string output   = "composed_integration";
    bool writeGeom        = false;
    bool plot             = false;

    gsCmdLine cmd("Quadrature study for composed maps: finest-mesh vs super-mesh, "
                  "and refinement rate of the super-mesh.");
    cmd.addInt   ("p","degree",     "Degree of the analysis geometry S (panel a)", degree);
    cmd.addInt   ("q","degreeMap",  "Degree of sigma", degreeMap);
    cmd.addInt   ("m","numKnotsMap","Interior knots of sigma per direction", numKnotsMap);
    cmd.addReal  ("A","amplitude",  "Amplitude of the sinusoidal displacement of S", amplitude);
    cmd.addReal  ("d","delta",      "Amplitude of the deformation of sigma", delta);
    cmd.addInt   ("r","maxRefine",  "Maximum dyadic refinement level of the super-mesh (panel b)", maxRefine);
    cmd.addString("o","output",     "Output file prefix", output);
    cmd.addSwitch("writeGeom",      "Write the p=3, numKnots=3 geometry to <prefix>_geometry.xml", writeGeom);
    cmd.addSwitch("plot",           "Export ParaView files", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFunctionExpr<T> phi("(1-x)^2 + (y-x)^2",2);
    const T exact = 0.5;

    gsInfo << std::scientific << std::setprecision(4);

    T minDetJ_S_all     = std::numeric_limits<T>::max();
    T minDetJ_sigma_all = std::numeric_limits<T>::max();

    // ---------------------------------------------------------------
    // Panel (b): refinement rate, p = 2 and p = 3, nQ = 3, three sigma modes
    // ---------------------------------------------------------------
    gsInfo << "\n### Panel (b): super-mesh error vs refinement level\n";
    const index_t nQb = 3;
    gsMatrix<T> refineTable(maxRefine+1, 10);
    std::vector<std::string> refineHeaders = {"level","h",
        "err_id_p2","err_dir_p2","err_coup_p2","nqp_p2",
        "err_id_p3","err_dir_p3","err_coup_p3","nqp_p3"};

    gsTensorBSpline<2,T> S_p3 = makeS(3, nQb, amplitude);
    T errIdP2_r0 = 0.0, errIdP3_r0 = 0.0;

    for (int pIdx=0; pIdx<2; ++pIdx)
    {
        const index_t p = (pIdx==0) ? 2 : 3;
        gsTensorBSpline<2,T> S = (p==3) ? S_p3 : makeS(2, nQb, amplitude);
        const T mS = minDetJ(S);
        minDetJ_S_all = std::min(minDetJ_S_all, mS);
        gsInfo << "# S (panel b): p=" << p << " nQ=" << nQb
               << "  minDetJ_S=" << mS << "\n";

        for (int mode=0; mode<3; ++mode)
        {
            gsTensorBSpline<2,T> sigma = makeSigma(degreeMap, numKnotsMap, mode, delta);
            const T mSig = minDetJ(sigma);
            minDetJ_sigma_all = std::min(minDetJ_sigma_all, mSig);
            gsInfo << "# sigma (panel b): mode=" << mode << " p=" << p
                   << "  minDetJ_sigma=" << mSig << "\n";

            gsComposedGeometry<T> G(sigma,S);
            for (index_t r=0; r<=maxRefine; ++r)
            {
                gsTensorBSplineBasis<2,T> ib = superMesh(S.basis(), sigma.basis(), (int)r);
                index_t nQP;
                const T val = integrate(G, phi, ib, 1, 1, nQP);
                const T err = math::abs(val-exact)/exact;

                refineTable(r,0) = static_cast<T>(r);
                refineTable(r,1) = 1.0/(12.0*std::pow(2.0,(double)r));
                const index_t errCol = (p==2) ? 2+mode : 6+mode;
                refineTable(r,errCol) = err;
                refineTable(r, (p==2) ? 5 : 9) = static_cast<T>(nQP);

                if (mode==0 && r==0)
                {
                    if (p==2) errIdP2_r0 = err; else errIdP3_r0 = err;
                }

                gsInfo << "b, p=" << p << ", mode=" << mode << ", r=" << r
                       << ", h=" << refineTable(r,1) << ", nEl=" << ib.numElements()
                       << ", nQP=" << nQP << ", relerr=" << err << "\n";
            }
        }
    }
    gsWriteCsv(output+"_refine.csv", refineTable, refineHeaders);

    // ---------------------------------------------------------------
    // Panel (a): finest mesh vs super mesh, p = degree, coupling sigma
    // ---------------------------------------------------------------
    gsInfo << "\n### Panel (a): finest mesh vs super mesh, p=" << degree << ", coupling sigma\n";
    gsMatrix<T> supermeshTable(6,5);
    std::vector<std::string> supermeshHeaders = {"numEL","err_fine","err_super","nqp_fine","nqp_super"};

    gsTensorBSpline<2,T> sigma_a = makeSigma(degreeMap, numKnotsMap, 2, delta);
    {
        const T mSig = minDetJ(sigma_a);
        minDetJ_sigma_all = std::min(minDetJ_sigma_all, mSig);
        gsInfo << "# sigma (panel a): coupling, minDetJ_sigma=" << mSig << "\n";
    }

    for (index_t numEL=1; numEL<=6; ++numEL)
    {
        const index_t nQ = numEL-1;
        gsTensorBSpline<2,T> S = makeS((int)degree,(int)nQ,amplitude);
        const T mS = minDetJ(S);
        minDetJ_S_all = std::min(minDetJ_S_all, mS);
        gsInfo << "# S (panel a): numEL=" << numEL << " minDetJ_S=" << mS << "\n";

        gsComposedGeometry<T> G(sigma_a,S);

        gsTensorBSplineBasis<2,T> superB = superMesh(S.basis(), sigma_a.basis(), 0);
        gsTensorBSplineBasis<2,T> fineB  = finestMesh(S.basis(), sigma_a.basis());
        index_t nqpSuper;
        const T valSuper = integrate(G, phi, superB, 1, 1, nqpSuper);
        const T errSuper = math::abs(valSuper-exact)/exact;

        T errFine; index_t nqpFine;
        if (fineB.knots(0)==superB.knots(0) && fineB.knots(1)==superB.knots(1))
        {
            // The two constructions produced the same knot vectors (either
            // S's and sigma's element counts tie, or sigma's knots are
            // already nested in S's): reuse the super-mesh result instead
            // of a second quadrature call, so the two columns are
            // bit-identical rather than merely equal up to floating-point
            // summation order.
            errFine = errSuper;
            nqpFine = nqpSuper;
        }
        else
        {
            const T valFine = integrate(G, phi, fineB, 1, 1, nqpFine);
            errFine = math::abs(valFine-exact)/exact;
        }

        supermeshTable(numEL-1,0) = static_cast<T>(numEL);
        supermeshTable(numEL-1,1) = errFine;
        supermeshTable(numEL-1,2) = errSuper;
        supermeshTable(numEL-1,3) = static_cast<T>(nqpFine);
        supermeshTable(numEL-1,4) = static_cast<T>(nqpSuper);

        gsInfo << "a, numEL=" << numEL << ", err_fine=" << errFine
               << ", err_super=" << errSuper << ", nqp_fine=" << nqpFine
               << ", nqp_super=" << nqpSuper << "\n";
    }
    gsWriteCsv(output+"_supermesh_p"+util::to_string(degree)+".csv", supermeshTable, supermeshHeaders);

    // ---------------------------------------------------------------
    // Self-check
    // ---------------------------------------------------------------
    gsInfo << "\n### self-check\n";
    index_t nQPref;
    const T ref = integrate(S_p3, phi, S_p3.basis(), 3, 5, nQPref);
    const T refErr = math::abs(ref-exact);
    const bool refPass  = refErr <= 1e-12;
    const bool idP2Pass = errIdP2_r0 <= 1e-12;
    const bool idP3Pass = errIdP3_r0 <= 1e-12;

    gsInfo << "|ref-0.5| = " << refErr << "  " << (refPass ? "PASS" : "FAIL") << "\n";
    gsInfo << "identity-sigma super-mesh error, p=2, r=0: " << errIdP2_r0
           << "  " << (idP2Pass ? "PASS" : "FAIL") << "\n";
    gsInfo << "identity-sigma super-mesh error, p=3, r=0: " << errIdP3_r0
           << "  " << (idP3Pass ? "PASS" : "FAIL") << "\n";
    gsInfo << "minDetJ_S (min over all cases)     = " << minDetJ_S_all << "\n";
    gsInfo << "minDetJ_sigma (min over all cases) = " << minDetJ_sigma_all << "\n";

    if (writeGeom)
    {
        gsTensorBSpline<2,T> Sgeom = makeS(3,3,amplitude);
        gsWrite(Sgeom, output+"_geometry");
        gsInfo << "Wrote " << output << "_geometry.xml\n";
    }

    if (plot)
    {
        gsWriteParaview(S_p3, output+"_S", 1000, true, true);
        gsTensorBSpline<2,T> sigmaPlot = makeSigma(degreeMap, numKnotsMap, 2, delta);
        gsWriteParaview(sigmaPlot, output+"_sigma", 1000, true, true);
        gsComposedGeometry<T> Gplot(sigmaPlot, S_p3);
        gsWriteParaview(Gplot, output+"_composed", 1000, true, true);

        gsTensorBSplineBasis<2,T> ibPlot = superMesh(S_p3.basis(), sigmaPlot.basis(), 0);
        gsOptionList plotOpt;
        plotOpt.addReal("quA","",1);
        plotOpt.addInt ("quB","",1);
        gsMatrix<T> nodesI = gsQuadrature::getAllNodes(ibPlot, plotOpt);
        gsMatrix<T> cnodesI;
        sigmaPlot.eval_into(nodesI,cnodesI);
        gsWriteParaviewPoints(nodesI, output+"_nodesI");
        gsWriteParaviewPoints(cnodesI, output+"_cnodesI");
    }

    const bool ok = refPass && idP2Pass && idP3Pass
                  && (minDetJ_S_all >= 0.1) && (minDetJ_sigma_all > 0.0);
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
