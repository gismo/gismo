/** @file reproduce_poisson_convergence.cpp

    @brief Reproduces the surface-Poisson convergence study of the
           composite-spline-relocation paper (the L2/H1 table in Section 4.3).

    It solves, on a spherical patch, the Poisson problem with the manufactured
    solution  u_ex = tanh( kappa*(1/2 - v) ) + 1  (v = polar angle, kappa = 20),
    for three parameterisations:
        - Original : identity composition sigma;
        - d2 NxN   : sigma re-parameterised with the gradient-based monitor
                     (paper Eq. 12) using the exact solution as indicator.
    For each it reports the L2 and H1 errors under uniform h-refinement of the
    analysis space, so that the three columns can be compared with the paper.

    Setup (r-adaptivity): the geometry is kept exact through the composed map
    G = S o sigma, while the analysis/solution space is an independent
    tensor-product B-spline of bi-degree (p,p) on the fixed computational square
    [0,1]^2, refined uniformly.  The mesh concentration enters solely through G.
    Integration uses the product-degree super-mesh (paper Section 4.4).

    Usage:
        reproduce_poisson_convergence -i <sphere_patch.xml> [-p 2] [-r 3]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji (y.ji-1@tudelft.nl)
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsOptim/gsOptim.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

using namespace gismo;

// Product-degree super-mesh integration basis (paper Section 4.4): parent
// analysis basis with the composition knot-lines inserted and degree raised to
// the product of the two degrees.
static gsTensorBSplineBasis<2,real_t>
integrationBasis(const gsTensorBSplineBasis<2,real_t> & analysis,
                 const gsTensorBSplineBasis<2,real_t> & composition)
{
    gsTensorBSplineBasis<2,real_t> ib(analysis);
    for (short_t d = 0; d != 2; ++d)
    {
        for (gsKnotVector<real_t>::uiterator it = std::next(composition.knots(d).ubegin());
             it != std::prev(composition.knots(d).uend()); ++it)
            if (!ib.knots(d).has(*it)) ib.insertKnot(*it,d);
        const index_t target = ib.degree(d) * composition.degree(d);
        ib.degreeIncrease(target - ib.degree(d), d);
    }
    return ib;
}

// Solve the surface Poisson problem for a given composition sigma, returning
// L2/H1 errors over (numRefine+1) uniform refinements of the deg-p analysis space.
static void solvePoisson(const gsGeometry<real_t> & S,       // exact surface (2->3)
                         const gsGeometry<real_t> & sigma,   // composition (2->2)
                         index_t p, index_t numRefine,
                         gsVector<real_t> & l2, gsVector<real_t> & h1,
                         gsVector<real_t> & dofs)
{
    // Composed geometry map G = S o sigma on [0,1]^2.
    gsMultiPatch<> cmp;  cmp.addPatch(gsComposedGeometry<real_t>(sigma, S));

    // Independent analysis space: deg-p tensor B-spline on [0,1]^2 (start 1 elem).
    gsKnotVector<> kv(0,1,0,p+1);
    gsTensorBSplineBasis<2> analysis(kv,kv);
    gsMultiBasis<> mb; mb.addBasis(analysis.clone().release());

    const gsTensorBSplineBasis<2> & compBasis =
        static_cast<const gsTensorBSplineBasis<2>&>(sigma.basis());

    // Manufactured solution and source in physical coordinates (kappa = 1/t = 20).
    const index_t dim = S.targetDim();
    std::string vstr = "atan2(sqrt(x^2+y^2),z)";
    std::string ms = "tanh((1/2-("+vstr+"))/0.05) + 1";
    std::string ff = "u:=atan2(y,x); v:="+vstr+"; w:=sqrt(x^2+y^2+z^2); t:=0.05;"
                     "-(-cos(v)*(1 - tanh((1/2 - v)/t)^2)/t - 2*sin(v)*tanh((1/2 - v)/t)*(1 - tanh((1/2 - v)/t)^2)/t^2)/(w^2*sin(v))";
    gsFunctionExpr<> msFun(ms,dim), fFun(ff,dim);
    gsComposedFunction<real_t> cms(cmp.patch(0), msFun), cf(cmp.patch(0), fFun);

    gsBoundaryConditions<> bc;
    for (short_t s = 1; s <= 4; ++s)
        bc.addCondition(s, condition_type::dirichlet, &cms, 0, true);
    bc.setGeoMap(cmp);

    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsMatrix<> solVector;
    gsExprEvaluator<> ev(A);
    A.options().setReal("quA",1.0);  A.options().setInt("quB",1);
    ev.options().setReal("quA",1.0); ev.options().setInt("quB",1);
    gsOptionList seOff; seOff.addSwitch("SameElement","",false);
    A.options().update(seOff, gsOptionList::addIfUnknown);
    ev.options().update(seOff, gsOptionList::addIfUnknown);

    geometryMap G = A.getMap(cmp);
    space u = A.getSpace(mb);
    auto fsrc = A.getCoeff(cf);
    auto uex  = ev.getVariable(cms);
    solution usol = A.getSolution(u, solVector);
    typename gsSparseSolver<>::QR solver;

    l2.setZero(numRefine+1); h1.setZero(numRefine+1); dofs.setZero(numRefine+1);
    for (index_t i = 0; i <= numRefine; ++i)
    {
        mb.uniformRefine();
        const gsTensorBSplineBasis<2> & ab =
            static_cast<const gsTensorBSplineBasis<2>&>(mb.basis(0));
        gsMultiBasis<> ib; ib.addBasis(integrationBasis(ab, compBasis).clone().release());
        A.setIntegrationElements(ib);
        ev.setIntegrationElements(ib);

        u.setup(bc, dirichlet::l2Projection, 0);
        A.initSystem();
        A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G),
                    u * fsrc * meas(G) );
        solver.compute(A.matrix());
        solVector = solver.solve(A.rhs());

        l2[i] = math::sqrt( ev.integral( (uex - usol).sqNorm() * meas(G) ) );
        h1[i] = l2[i] + math::sqrt( ev.integral(
                    ( igrad(uex,G) - igrad(usol,G) ).sqNorm() * meas(G) ) );
        dofs[i] = A.numDofs();
    }
}

// Build an identity or gradient-based re-parameterised composition of bi-degree
// (compDeg,compDeg) with NxN elements on [0,1]^2.
static gsGeometry<real_t>::uPtr makeSigma(const gsGeometry<real_t> & S,
                                          index_t compDeg, index_t N, bool reparam,
                                          real_t smoothing, real_t penalty)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(compDeg);
    const index_t nref = (index_t)std::round(std::log2((double)N));
    for (index_t r = 0; r < nref; ++r) comp->uniformRefine();
    gsSquareDomain<real_t> domain(*comp);
    domain.options().addSwitch("Slide","",false);   // fixed boundary (paper Sec. 4.3)
    domain.applyOptions();

    if (reparam && domain.nControls() > 0)
    {
        // Indicator = exact solution; gradient-based monitor (paper Eq. 12).
        std::string ms = "tanh((1/2-(atan2(sqrt(x^2+y^2),z)))/0.05) + 1";
        gsFunctionExpr<> indicator(ms, S.targetDim());
        gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",250);
        gsAdaptiveParametrization<real_t,MonitorMode::GradientBased>
            relo(domain, S, indicator, opt, /*parametric*/false);
        relo.options().setReal("Smoothing",smoothing);
        relo.options().setReal("Penalty",penalty);

        // Guard against L-BFGS silently returning the initial (identity) design after a
        // failed first line search -- that would tabulate the ORIGINAL geometry's errors
        // in the "re-parameterised" columns. See reproduce_surface_metrics.cpp.
        const gsVector<real_t> before = domain.getControls();
        relo.solve();
        const bool moved = (domain.getControls()-before).norm() > 1e-12;

        const real_t minJs = relo.computeMinJacobian();
        GISMO_ENSURE(math::isfinite(minJs) && minJs > 0,
                     "Re-parameterisation folded/degenerate (min det Jsigma = "<<minJs<<").");
        if (!moved)
            gsWarn<<"  L-BFGS did not move; the re-parameterised columns will equal the "
                    "Original ones. Try a larger --penalty.\n";
        gsInfo << "  [reparam d="<<compDeg<<" "<<N<<"x"<<N
               << ", min det Jsigma = "<<minJs<<"]\n";
    }
    return domain.domain().clone();
}

int main(int argc, char** argv)
{
    std::string input = "parametrization/monitor_example_sphere_patch_with_gradient_function.xml";
    index_t p = 2, numRefine = 3, compDeg = 2;
    real_t smoothing = 1e-2, penalty = 1e-8;
    gsCmdLine cmd("Reproduce the Section-4.3 Poisson convergence table.");
    cmd.addString("i","input","Sphere-patch xml (uses 'geometry' label or first geometry)",input);
    cmd.addInt("p","analysisDegree","Degree of the analysis space",p);
    cmd.addInt("r","refine","Number of uniform analysis refinements",numRefine);
    cmd.addInt("d","compDegree","Degree of the composition sigma",compDeg);
    cmd.addReal("S","smoothing","Monitor smoothing theta",smoothing);
    cmd.addReal("P","penalty","Fold-barrier penalty",penalty);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Fail loudly on a missing input: the default path is relative to the CWD, and a
    // silent empty run is far worse than an error message.
    GISMO_ENSURE(gsFileManager::fileExists(input),
                 "Input file not found: '"<<input<<"'.\n"
                 "  The default path is relative to the current directory; pass an "
                 "absolute path with -i <sphere_patch.xml>.");

    gsFileData<> fd(input);
    gsMultiPatch<> geomMP;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",geomMP);
    else                         geomMP.addPatch(*fd.getFirst<gsGeometry<>>());
    const gsGeometry<>& S = geomMP.patch(0);
    gsInfo << "Analysis degree p="<<p<<", composition degree d="<<compDeg
           << ", "<<numRefine+1<<" refinement levels.\n\n";

    struct Col { std::string name; bool reparam; index_t N; };
    std::vector<Col> cols = { {"Original",false,4}, {"d"+std::to_string(compDeg)+" 4x4",true,4},
                              {"d"+std::to_string(compDeg)+" 8x8",true,8} };

    std::vector<gsVector<real_t>> L2, H1, DOF;
    for (const Col & c : cols)
    {
        gsInfo << "=== "<<c.name<<" ===\n";
        gsGeometry<>::uPtr sigma = makeSigma(S, compDeg, c.N, c.reparam, smoothing, penalty);
        gsVector<real_t> l2,h1,dofs;
        solvePoisson(S, *sigma, p, numRefine, l2, h1, dofs);
        L2.push_back(l2); H1.push_back(h1); DOF.push_back(dofs);
    }

    gsInfo << std::scientific << std::setprecision(4);
    gsInfo << "\n---- L2 error ----\n" << std::left << std::setw(8) << "level";
    for (const Col & c : cols) gsInfo << std::setw(14) << c.name;
    gsInfo << "\n";
    for (index_t i = 0; i <= numRefine; ++i)
    {
        gsInfo << std::left << std::setw(8) << ((int)DOF[0][i]);
        for (size_t c = 0; c < cols.size(); ++c) gsInfo << std::setw(14) << L2[c][i];
        gsInfo << "\n";
    }
    gsInfo << "\n---- H1 error ----\n" << std::left << std::setw(8) << "level";
    for (const Col & c : cols) gsInfo << std::setw(14) << c.name;
    gsInfo << "\n";
    for (index_t i = 0; i <= numRefine; ++i)
    {
        gsInfo << std::left << std::setw(8) << ((int)DOF[0][i]);
        for (size_t c = 0; c < cols.size(); ++c) gsInfo << std::setw(14) << H1[c][i];
        gsInfo << "\n";
    }
    gsInfo << "\n(level = analysis DoF; compare trends with the paper's Section-4.3 table)\n";
    return 0;
}
