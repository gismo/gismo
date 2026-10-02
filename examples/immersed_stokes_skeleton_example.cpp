/** @file immersed_stokes_skeleton_example.cpp

    @brief Skeleton/ghost-stabilized immersed Stokes on a quarter annulus.

    Solves the incompressible Stokes system

        -2*mu*div(grad^s u) + grad p = f,   div u = 0   in Omega = { 1 < r < 4, x > 0, y > 0 }

    immersed in a [0,L]^2 background mesh, with an equal-order velocity/
    pressure pair (one vector velocity space and one scalar pressure space
    on the same basis and degree, in a two-block gsExprAssembler<>(2,2)).
    The volume form is the symmetric-gradient one,
    mu*(2*grad^s u : grad^s w) - (p, div w) - (q, div u), matching the
    traction sigma(u,p)*n = (2*mu*grad^s u - p*I)*n used by the Nitsche
    terms below.

    u = 0 is imposed on the level set r = 1, r = 4 by Nitsche's method with
    a traction consistent with the volume form (both consistency and
    adjoint-consistency terms, plus a penalty and the two pressure
    couplings); u = 0 on the conforming box edges x = 0, y = 0 is imposed
    strongly by DOF elimination. Ghost-penalty stabilization on the
    velocity jump and skeleton (pressure) stabilization on cut faces cure
    the small-cut-cell conditioning and the equal-order inf-sup deficiency;
    exterior-DOF elimination removes basis functions with no physical
    support.

    Geometry and level set are expressed in ONE coordinate system: the
    background patch is the IDENTITY map (BSplineRectangleWithPara), so a
    gsImplicitTrimmedDomain, which classifies in parameter space, and
    getCoeff(f, G), which evaluates f in physical space via the geometry
    map, see the same coordinates. The box side is L = 4.2, not 4: on
    [0,4]^2 the outer circle r = 4 is exactly tangent to the box edges at
    (4,0) and (0,4), and the Algoim surface rule hangs on exact tangency;
    4.2 keeps every power-of-two element count away from a mesh line
    landing on x = 4 (that would require 21 | N). x = 0 and y = 0 remain
    conforming box edges.

    The level set is written in the smooth polynomial form
    (x^2+y^2-1)*(x^2+y^2-16) rather than max(1-r^2, r^2-16): both have the
    same sign field, but the polynomial form is C^inf, which is what the
    Algoim rule's Bernstein approximation needs, whereas max() has an
    interior kink (at r^2 = 8.5) inside the physical domain.

    Manufactured solution: u carries the factors (r^2-1)(r^2-16), so it
    vanishes on the ENTIRE immersed boundary (both circles), and g = 0 on
    the conforming edges too. Every Nitsche right-hand-side term is
    therefore identically zero and is omitted; only the matrix terms and
    the volume load (f, w) survive. selfCheck() below is what licenses
    this: it certifies the pasted solution strings before any assembly
    runs.

    The Nitsche traction pairs a differentiated vector trial function with
    an undifferentiated vector test function -- a "vertical" (na*dim) x dim
    local-matrix shape that has no other live call site in this tree. Its
    two consistency halves are certified once, numerically, by
    nitscheTractionCheck() before any convergence table is trusted (see
    that function and the derivation in the comments around tracT/tracN).

    Reference: the paper this driver reproduces calibrates gamma (skeleton)
    and gammaT (ghost) on a [0,4]^2, 11x11 background mesh; the scaling
    exponents h^(2k+1) and h^(2k-1) are applied explicitly here so the same
    constants transfer to other refinement levels.

    Example command lines:
      ./immersed_stokes_skeleton_example --study rates -p 1 -r 3
      ./immersed_stokes_skeleton_example --study rates -p 2 -r 3
      ./immersed_stokes_skeleton_example --study inertia
      ./immersed_stokes_skeleton_example --study gammaSweep
      ./immersed_stokes_skeleton_example --study gammaTSweep

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsCore/gsDofMapper.h>
#include <gsAssembler/gsDofMapperCreator.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

namespace {

std::string fmtSci(real_t v, int prec = 6)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

const std::string faceShiftLegend =
    "Face integrals are evaluated at x -/+ faceShift*h (faceShift = 1e-6): quantities of "
    "order (faceShift*h)^2 are the evaluation-offset floor, not convergence to zero.";

// Integer grid coordinate / flat level-0 element index of a point on the
// uniform N-element grid of [0,L]^2, used only by the --sliverDelta guard.
index_t gridCoord(const real_t x, const real_t L, const index_t n)
{ return static_cast<index_t>(math::round(x / L * (real_t)n)); }

size_t flatIndex(const gsVector<real_t> & p, const real_t L, const index_t n)
{
    size_t flat = 0, stride = 1;
    for (index_t j = 0; j != p.rows(); ++j)
    { flat += (size_t)gridCoord(p[j], L, n) * stride; stride *= (size_t)n; }
    return flat;
}

short_t sliverSign(const gsImplicitTrimmedDomain<2,real_t> & dom, size_t targetFlat,
                    real_t L, index_t N)
{
    gsDomain<real_t>::iterator eInt = dom.end<InteriorSign>();
    for (gsDomain<real_t>::iterator it = dom.beginInterior(); it < eInt; ++it)
        if (flatIndex(it.lowerCorner(), L, N) == targetFlat) return -1;
    gsDomain<real_t>::iterator eCut = dom.endBdr(boundary::none);
    for (gsDomain<real_t>::iterator it = dom.beginBdr(boundary::none); it < eCut; ++it)
        if (flatIndex(it.lowerCorner(), L, N) == targetFlat) return 0;
    return +1;
}

/// Four checks on the pasted manufactured-solution strings (pointwise
/// values at six reference points, divergence, boundary vanishing, and a
/// finite-difference PDE residual), certified once at startup rather than
/// trusted by inspection: mistyped exponents or a wrong sign in an
/// 800-character string would otherwise silently degrade the observed
/// convergence rate instead of failing visibly.
bool selfCheck(const gsFunctionExpr<real_t> & u_exact,
               const gsFunctionExpr<real_t> & p_exact,
               const gsFunctionExpr<real_t> & f_rhs)
{
    static const real_t ref[6][7] = {
        {1.3000, 0.7000,  7.03481753539671e-04, -4.20661016860827e-04, -3.81049573059987e-01, -5.55056348093649e-01,  2.63439754440807e-01},
        {2.0000, 2.0000,  1.06086400000000e+00, -4.58752000000000e-01,  0.00000000000000e+00,  8.79019361667012e-01, -1.70845936166701e+00},
        {3.0000, 1.1000,  1.75721407720466e-01,  4.24296068004257e-03, -5.84460995653630e-01,  2.96976982485217e-01, -6.39177591370341e-01},
        {0.5000, 1.9000,  3.47962790638383e-02, -3.38385876270254e-02,  4.78556608689874e-01, -1.30962207684048e-01,  8.26171412237808e-01},
        {1.1000, 0.4500,  1.24109285143761e-05, -9.84962470619402e-06, -2.35783617346219e-01, -1.01317095003407e+00, -4.26348348507552e-01},
        {2.7000, 2.6000, -1.46580892963441e+00,  2.22684374513687e+00, -1.00923053392994e-02, -3.46304322855461e+01,  4.58167530796724e+01}
    };

    const real_t atol = 1e-12, rtol = 1e-10;
    auto mixedOk = [&](real_t a, real_t b)
    { return math::abs(a - b) <= atol + rtol * math::abs(b); };

    // Item 1: pointwise values against the reference table.
    real_t worst1 = 0;
    gsMatrix<real_t> pt(2,1), val;
    for (index_t i = 0; i != 6; ++i)
    {
        pt(0,0) = ref[i][0]; pt(1,0) = ref[i][1];
        u_exact.eval_into(pt, val);
        const real_t u1 = val(0,0), u2 = val(1,0);
        p_exact.eval_into(pt, val);
        const real_t p = val(0,0);
        f_rhs.eval_into(pt, val);
        const real_t f1 = val(0,0), f2 = val(1,0);

        const real_t got[5]  = { u1, u2, p, f1, f2 };
        const real_t want[5] = { ref[i][2], ref[i][3], ref[i][4], ref[i][5], ref[i][6] };
        for (index_t c = 0; c != 5; ++c)
        {
            if (!mixedOk(got[c], want[c]))
            {
                gsWarn << "selfCheck item 1 FAILED at point (" << ref[i][0] << "," << ref[i][1]
                       << "), component " << c << ": got " << got[c] << " want " << want[c] << "\n";
                return false;
            }
            worst1 = math::max(worst1, math::abs(got[c] - want[c]));
        }
    }

    // Item 2: divergence of u by central differences, step h = 1e-4.
    const real_t h = 1e-4;
    real_t worst2 = 0;
    for (index_t i = 0; i != 6; ++i)
    {
        const real_t x = ref[i][0], y = ref[i][1];
        gsMatrix<real_t> pxp(2,1), pxm(2,1), pyp(2,1), pym(2,1), v;
        pxp << x+h, y; pxm << x-h, y; pyp << x, y+h; pym << x, y-h;
        u_exact.eval_into(pxp, v); const real_t u1xp = v(0,0);
        u_exact.eval_into(pxm, v); const real_t u1xm = v(0,0);
        u_exact.eval_into(pyp, v); const real_t u2yp = v(1,0);
        u_exact.eval_into(pym, v); const real_t u2ym = v(1,0);
        const real_t div = (u1xp - u1xm) / (2*h) + (u2yp - u2ym) / (2*h);
        worst2 = math::max(worst2, math::abs(div));
    }
    if (worst2 > 1e-6)
    { gsWarn << "selfCheck item 2 FAILED: worst |div u| = " << worst2 << " > 1e-6\n"; return false; }

    // Item 3: u vanishes on the immersed boundary and on the conforming
    // edges, within the physical annulus.
    real_t worst3 = 0;
    for (index_t a = 0; a != 17; ++a)
    {
        const real_t theta = (real_t)a / 16.0 * (M_PI / 2.0);
        for (real_t r : { 1.0, 4.0 })
        {
            pt(0,0) = r * math::cos(theta); pt(1,0) = r * math::sin(theta);
            u_exact.eval_into(pt, val);
            worst3 = math::max(worst3, math::max(math::abs(val(0,0)), math::abs(val(1,0))));
        }
    }
    for (index_t a = 0; a != 17; ++a)
    {
        const real_t r = 1.0 + (real_t)a / 16.0 * 3.0;
        pt(0,0) = 0.0; pt(1,0) = r; u_exact.eval_into(pt, val);
        worst3 = math::max(worst3, math::max(math::abs(val(0,0)), math::abs(val(1,0))));
        pt(0,0) = r; pt(1,0) = 0.0; u_exact.eval_into(pt, val);
        worst3 = math::max(worst3, math::max(math::abs(val(0,0)), math::abs(val(1,0))));
    }
    if (worst3 > 1e-10)
    { gsWarn << "selfCheck item 3 FAILED: worst |u| on boundary = " << worst3 << " > 1e-10\n"; return false; }

    // Item 4: pasted f against a finite-difference -mu*Lap(u) + grad(p).
    const real_t mu = 1.0;
    real_t worst4 = 0;
    for (index_t i = 0; i != 6; ++i)
    {
        const real_t x = ref[i][0], y = ref[i][1];
        gsMatrix<real_t> pc(2,1), v;
        auto U = [&](real_t xx, real_t yy, index_t comp)->real_t
        { pc << xx, yy; u_exact.eval_into(pc, v); return v(comp,0); };
        auto P = [&](real_t xx, real_t yy)->real_t
        { pc << xx, yy; p_exact.eval_into(pc, v); return v(0,0); };

        const real_t lap_u1 = (U(x+h,y,0) - 2*U(x,y,0) + U(x-h,y,0)) / (h*h)
                             + (U(x,y+h,0) - 2*U(x,y,0) + U(x,y-h,0)) / (h*h);
        const real_t lap_u2 = (U(x+h,y,1) - 2*U(x,y,1) + U(x-h,y,1)) / (h*h)
                             + (U(x,y+h,1) - 2*U(x,y,1) + U(x,y-h,1)) / (h*h);
        const real_t px = (P(x+h,y) - P(x-h,y)) / (2*h);
        const real_t py = (P(x,y+h) - P(x,y-h)) / (2*h);

        const real_t fd1 = -mu*lap_u1 + px;
        const real_t fd2 = -mu*lap_u2 + py;

        pc << x, y; f_rhs.eval_into(pc, v);
        const real_t f1 = v(0,0), f2 = v(1,0);

        const real_t rel1 = math::abs(fd1 - f1) / math::max(math::abs(f1), 1.0);
        const real_t rel2 = math::abs(fd2 - f2) / math::max(math::abs(f2), 1.0);
        worst4 = math::max(worst4, math::max(rel1, rel2));
    }
    if (worst4 > 1e-3)
    { gsWarn << "selfCheck item 4 FAILED: worst relative f mismatch = " << worst4 << " > 1e-3\n"; return false; }

    gsInfo << "selfCheck: worst pointwise deviation = " << fmtSci(worst1)
           << ", worst |div u| = " << fmtSci(worst2)
           << ", worst |u| on boundary = " << fmtSci(worst3)
           << ", worst relative f mismatch = " << fmtSci(worst4) << "\n";
    gsInfo << "selfCheck: OK\n";
    return true;
}

/// Numerically certifies the two consistency halves of the Nitsche
/// traction (tracT, tracN below) against a plain Space==0 evaluation of
/// the same physical quantity on a random discrete field, since no live
/// in-tree call site pairs a differentiated vector trial function with an
/// undifferentiated vector test function. See the file header.
bool nitscheTractionCheck(const gsFunctionExpr<real_t> & phi,
                           const gsFunctionExpr<real_t> & nImm)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    const real_t L = 4.2, mu = 1.0;
    const index_t N0 = 8, degree = 1;

    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineRectangleWithPara(0.0, 0.0, L, L));
    gsMultiBasis<real_t> dbasis(mp, true);
    dbasis.setDegree(degree);
    dbasis.uniformRefine(N0 - 1);

    gsTensorBSplineBasis<2,real_t> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");
    memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > tr =
        memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phi, *tbsPtr));
    GISMO_ENSURE(1 == tr->numLevels(), "gsTrimmedDomain requires a single kd-tree level.");
    tr->numGhostFaces();   // warm the lazily-built face sign grid, single-threaded

    std::vector<patchSide> bdr(1);
    bdr[0] = patchSide(0, boundary::none);

    gsExprAssembler<> A(1, 1);
    A.options().setInt("quRule", gsQuadrature::AlgoimRule);

    geometryMap G = A.getMap(mp);
    space w = A.getSpace(dbasis, 2, 0);
    gsDofMapper mapper = createMapper(dbasis, 2, false); mapper.finalize(); w.setupMapper(mapper);
    const_cast<expr::gsFeSpace<real_t>&>(w).fixedPart().setZero(mapper.boundarySize(), 1);

    A.setIntegrationElements(dbasis);
    A.initSystem();
    A.setIntegrationDomain(tr);

    gsMatrix<real_t> a = gsMatrix<real_t>::Random(A.numDofs(), 1);
    solution ua = A.getSolution(w, a);

    auto n_imm    = A.getCoeff(nImm, G);
    auto surfMeas = meas(G) * (jac(G).inv().tr() * n_imm).norm();

    // (grad w)^T . n and (grad w) . n, each paired with the undifferentiated
    // trial function w: tracT contracts the Jacobian's second index against
    // n after a block-wise transpose (keeps Space == 2, the trial side),
    // tracN contracts the first index and re-transposes with cwisetr() so
    // Space stays 2 as well -- see solveOne() for the identical spelling
    // used in the assembled system.
    auto tracT = ijac(w, G).tr() * n_imm;
    auto tracN = (n_imm.tr() * ijac(w, G).tr()).cwisetr();

    A.options().addInt("quDim", "Surface (phi==0) quadrature selector", 2);
    A.assembleBdr(bdr, -mu*surfMeas*(w*tracT), -mu*surfMeas*(w*tracN));
    A.assembleBdr(bdr, -mu*surfMeas*(w * ((igrad(ua,G).cwisetr() + igrad(ua,G)) * n_imm)));
    A.options().setInt("quDim", -1);

    const real_t rel = (A.matrix()*a - A.rhs()).norm() / math::max(A.rhs().norm(), 1e-30);
    gsInfo << "nitscheTractionCheck: relative residual = " << fmtSci(rel) << "\n";
    return rel <= 1e-10;
}

} // anonymous namespace

/// Plain-scalar summary of one (degree, refinement level, parameters)
/// Stokes solve.
struct RunResult
{
    index_t nV = 0, nP = 0, ndofs = 0, nCut = 0, nGhost = 0, nSkel = 0;
    real_t  l2u = std::numeric_limits<real_t>::quiet_NaN();
    real_t  h1u = std::numeric_limits<real_t>::quiet_NaN();
    real_t  l2p = std::numeric_limits<real_t>::quiet_NaN();
    real_t  l2pShift = std::numeric_limits<real_t>::quiet_NaN();
    real_t  minRowNorm = 0;
    real_t  symRel = 0;
    index_t nPos = 0, nNeg = 0, nZero = 0;
    real_t  lambdaMinAbsOverMax = 0;
    real_t  smallestRatios[5] = {0,0,0,0,0};
    bool    haveSpectrum = false;
    bool    finite = true;
    real_t  L = 0, h = 0;
};

/// Assembles and solves one configuration. Builds a fresh gsExprAssembler
/// and a fresh trimmed domain (both from a per-call \a phi), never reused
/// across sweep points: assemble*() accumulates into the matrix and never
/// clears it.
RunResult solveOne(index_t degree, index_t N, real_t L,
                    const gsFunctionExpr<real_t> & phi,
                    const gsFunctionExpr<real_t> & nImm,
                    const gsFunctionExpr<real_t> & u_exact,
                    const gsFunctionExpr<real_t> & p_exact,
                    const gsFunctionExpr<real_t> & f_rhs,
                    real_t beta, real_t gamma, real_t gammaT,
                    index_t quRule, index_t quA,
                    bool noGhost, bool noSkeleton,
                    const std::string & pressureMode,
                    bool computeSpectrum, bool doPlot, const std::string & outPath)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    const real_t mu = 1.0;
    const index_t k = degree;

    RunResult result;
    result.L = L; result.h = L / (real_t)N;
    const real_t hmax = result.h;

    // 1. Assembler + quadrature options; quDim added once, volume default -1.
    gsExprAssembler<> A(2, 2);
    A.options().setInt("quRule", quRule);
    if (quA > 0) A.options().setReal("quA", (real_t)quA);
    A.options().addInt("quDim", "Surface (phi==0) quadrature selector", -1);

    // 2. Trimmed domain; single kd-tree level; warm the face sign grid.
    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineRectangleWithPara(0.0, 0.0, L, L));
    gsMultiBasis<real_t> dbasis(mp, true);
    dbasis.setDegree(degree);       // elevate BEFORE refining, or continuity drops to C^0
    dbasis.uniformRefine(7);        // level 0: 8 elements per direction
    {
        index_t n0 = 8;
        while (n0 < N) { n0 <<= 1; dbasis.uniformRefine(); }
    }

    gsTensorBSplineBasis<2,real_t> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");
    memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > tr_domain =
        memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phi, *tbsPtr));
    GISMO_ENSURE(1 == tr_domain->numLevels(), "gsTrimmedDomain requires a single kd-tree level.");
    result.nGhost = static_cast<index_t>(tr_domain->numGhostFaces());
    result.nCut   = static_cast<index_t>(tr_domain->numElementsBdr(boundary::none));
    result.nSkel  = static_cast<index_t>(tr_domain->numSkeletonFaces());

    // 3. Evaluator with its own quRule/quA options.
    gsExprEvaluator<> ev(A);
    ev.options().setInt("quRule", quRule);
    if (quA > 0) ev.options().setReal("quA", (real_t)quA);

    std::vector<patchSide> bdr_immersed(1);
    bdr_immersed[0] = patchSide(0, boundary::none);

    // 4. Map, spaces, coefficients, named solution vector.
    geometryMap G = A.getMap(mp);
    space v = A.getSpace(dbasis, 2, 0);   // vector velocity, block 0
    space p = A.getSpace(dbasis, 1, 1);   // scalar pressure, block 1

    auto ff       = A.getCoeff(f_rhs, G);
    auto n_imm    = A.getCoeff(nImm, G);
    auto surfMeas = meas(G) * (jac(G).inv().tr() * n_imm).norm();

    gsMatrix<real_t> solVector;   // NAMED: getSolution keeps a pointer to it
    solution v_sol = A.getSolution(v, solVector);
    solution p_sol = A.getSolution(p, solVector);

    // 5. Mappers: exterior elimination (both spaces), strong Dirichlet on
    //    x=0/y=0 (velocity), optional pressure pin.
    gsDofMapper mapperV = createMapper(dbasis, 2, false);
    gsDofMapper mapperP = createMapper(dbasis, 1, false);
    {
        std::vector<bool> keepV(dbasis.basis(0).size(), false);
        std::vector<bool> keepP(dbasis.basis(0).size(), false);
        std::vector<bool> interiorP(dbasis.basis(0).size(), false);
        gsMatrix<real_t> centre(2, 1);
        gsMatrix<index_t> act;

        gsDomain<real_t>::iterator eInt = tr_domain->end<InteriorSign>();
        for (gsDomain<real_t>::iterator it = tr_domain->beginInterior(); it < eInt; ++it)
        {
            centre.col(0) = 0.5 * (it.lowerCorner() + it.upperCorner());
            dbasis.basis(0).active_into(centre, act);
            for (index_t i = 0; i != act.rows(); ++i)
            { keepV[act(i,0)] = true; keepP[act(i,0)] = true; interiorP[act(i,0)] = true; }
        }
        gsDomain<real_t>::iterator eCut = tr_domain->endBdr(boundary::none);
        for (gsDomain<real_t>::iterator it = tr_domain->beginBdr(boundary::none); it < eCut; ++it)
        {
            centre.col(0) = 0.5 * (it.lowerCorner() + it.upperCorner());
            dbasis.basis(0).active_into(centre, act);
            for (index_t i = 0; i != act.rows(); ++i)
            { keepV[act(i,0)] = true; keepP[act(i,0)] = true; }
        }

        // Strong u = 0 on the conforming edges x = 0 (west), y = 0 (south).
        gsMatrix<index_t> bWest = dbasis.basis(0).boundary(boundary::west);
        gsMatrix<index_t> bSouth = dbasis.basis(0).boundary(boundary::south);
        for (index_t i = 0; i != bWest.rows(); ++i) keepV[bWest(i,0)] = false;
        for (index_t i = 0; i != bSouth.rows(); ++i) keepV[bSouth(i,0)] = false;

        for (index_t i = 0; i != (index_t)keepV.size(); ++i)
        {
            for (index_t c = 0; c != 2; ++c)
                if (!keepV[i]) mapperV.eliminateDof(i, 0, c);
        }
        for (index_t i = 0; i != (index_t)keepP.size(); ++i)
            if (!keepP[i]) mapperP.eliminateDof(i, 0);

        if ("pin" == pressureMode)
        {
            for (index_t i = 0; i != (index_t)keepP.size(); ++i)
            {
                if (keepP[i] && interiorP[i]) { mapperP.eliminateDof(i, 0); break; }
            }
        }

        mapperV.finalize();
        v.setupMapper(mapperV);
        const_cast<expr::gsFeSpace<real_t>&>(v).fixedPart().setZero(mapperV.boundarySize(), 1);

        mapperP.finalize();
        p.setupMapper(mapperP);
        const_cast<expr::gsFeSpace<real_t>&>(p).fixedPart().setZero(mapperP.boundarySize(), 1);
    }

    // 6. Integration elements, initSystem, dof-layout asserts.
    A.setIntegrationElements(dbasis);
    A.initSystem();

    const index_t nV = v.mapper().freeSize();
    const index_t nP = p.mapper().freeSize();
    GISMO_ENSURE(A.numDofs() == nV + nP,
                 "Multi-block dof layout mismatch: numDofs=" << A.numDofs()
                 << " nV+nP=" << nV+nP);
    GISMO_ENSURE(p.mapper().firstIndex() == nV,
                 "Pressure block does not start right after the velocity block: firstIndex="
                 << p.mapper().firstIndex() << " nV=" << nV);
    result.nV = nV; result.nP = nP; result.ndofs = A.numDofs();

    // 7. Sparsity pattern from the three volume matrix terms.
    auto Jv = ijac(v, G);
    A.computePattern(
          mu * ((Jv.cwisetr() + Jv) % Jv.tr()) * meas(G)
        , -idiv(v, G) * p.tr() * meas(G)
        , -p * idiv(v, G).tr() * meas(G)
    );

    // 8. Install the trimmed integration domain (after setIntegrationElements).
    A.setIntegrationDomain(tr_domain);

    // 9. Volume assembly: symmetric-gradient stiffness, the two pressure
    //    couplings, and the load.
    A.assemble(
          mu * ((Jv.cwisetr() + Jv) % Jv.tr()) * meas(G)
        , -idiv(v, G) * p.tr() * meas(G)
        , -p * idiv(v, G).tr() * meas(G)
        , v * ff * meas(G)
    );

    // 10. N3/N4/N5 -- already a symmetric set: N3 is the (0,0) penalty mass
    //     matrix v.v, and N4/N5 are transposes of each other by
    //     construction ((v.n)*p^T versus p*(v.n)^T), so this block alone
    //     assembles to a symmetric matrix before N1/N2 are added below.
    A.options().setInt("quDim", 2);
    A.assembleBdr(bdr_immersed,
          beta * mu / hmax * v * v.tr() * surfMeas          // N3: penalty
        ,  (v * n_imm) * p.tr()         * surfMeas          // N4: pressure consistency
        ,  p * (v * n_imm).tr()         * surfMeas          // N5: pressure adjoint, N4^T
    );
    A.options().setInt("quDim", -1);

    // 11. Ghost (velocity) and skeleton (pressure) stabilization -- BEFORE
    //     the M0 snapshot below. Ghost/skeleton are symmetric on their own,
    //     so assembling them after the snapshot would still pass the
    //     symmetry check while the N1 -> C -> K transposition trick in
    //     steps 12-14 doubled their contribution silently.
    auto dJv = dnk(v.jump(), G.left(), k);
    auto dJp = dnk(p.jump(), G.left(), k);
    if (!noGhost)
        A.assembleGhost(gammaT * math::pow(hmax, 2*k - 1) * dJv * dJv.tr());
    if (!noSkeleton)
        // The minus sign is load-bearing: it enters the (1,1) block
        // negatively, giving the saddle-point signature [[A,B^T],[B,-C]]
        // that the inertia diagnostic in --study inertia assumes.
        A.assembleSkeleton(-gamma * math::pow(hmax, 2*k + 1) / mu * dJp * dJp.tr());

    // 12. Snapshot, value copy, before the last matrix assembly (N1).
    const gsSparseMatrix<real_t> M0 = A.matrix();

    // 13. N1: the two halves of -mu*(2*grad^s u . n).w, paired as a trial
    //     (differentiated) x test (undifferentiated) block; certified by
    //     nitscheTractionCheck() at startup. N2 (the transposed pairing,
    //     test differentiated / trial undifferentiated) cannot be spelled
    //     directly in the expression language (see file header), so it is
    //     recovered below by transposing N1 instead.
    auto tracT = ijac(v, G).tr() * n_imm;
    auto tracN = (n_imm.tr() * ijac(v, G).tr()).cwisetr();
    A.options().setInt("quDim", 2);
    A.assembleBdr(bdr_immersed, -mu*surfMeas*(v*tracT), -mu*surfMeas*(v*tracN));
    A.options().setInt("quDim", -1);

    // 14. K = M0 + N1 + N1^T. The zero right-hand-side counterpart to N1^T
    //     relies on every eliminated DOF being fixed to zero (exterior
    //     elimination, strong Dirichlet on x=0/y=0, and the optional
    //     pressure pin): the eliminated-column term in the scatter loop is
    //     then localMat(...)*0 == 0 and needs no explicit RHS assembly.
    const gsSparseMatrix<real_t> C = A.matrix() - M0;
    gsSparseMatrix<real_t> K = A.matrix() + gsSparseMatrix<real_t>(C.transpose());

    // 15a. Symmetry check.
    result.symRel = (K - gsSparseMatrix<real_t>(K.transpose())).norm() / math::max(K.norm(), 1e-30);
    GISMO_ENSURE(result.symRel <= 1e-10, "K is not symmetric to 1e-10 relative: " << result.symRel);

    // 15b. Zero-row check.
    {
        gsVector<real_t> rowSq = gsVector<real_t>::Zero(K.rows());
        for (index_t c = 0; c != K.outerSize(); ++c)
            for (gsSparseMatrix<real_t>::InnerIterator it(K, c); it; ++it)
                rowSq[it.row()] += it.value() * it.value();
        index_t nZeroRows = 0;
        result.minRowNorm = std::numeric_limits<real_t>::max();
        for (index_t i = 0; i != rowSq.rows(); ++i)
        {
            result.minRowNorm = math::min(result.minRowNorm, math::sqrt(rowSq[i]));
            if (0.0 == rowSq[i]) ++nZeroRows;
        }
        GISMO_ENSURE(0 == nZeroRows, "K has " << nZeroRows << " identically-zero rows out of "
                     << rowSq.rows() << " (DOF elimination did not remove every unsupported DOF).");
    }

    gsInfo << "  nV=" << result.nV << " nP=" << result.nP << " numDofs=" << result.ndofs
           << " nCut=" << result.nCut << " nGhost=" << result.nGhost << " nSkel=" << result.nSkel
           << " minRowNorm=" << fmtSci(result.minRowNorm) << "\n";

    // 15c. Optional dense spectrum / inertia signature.
    if (computeSpectrum)
    {
        const gsMatrix<real_t> dense = K.toDense();
        gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> es(dense);
        GISMO_ENSURE(gsEigen::Success == es.info(), "Eigen solver failed.");
        const gsVector<real_t> ev_ = es.eigenvalues();
        const real_t lmaxAbs = ev_.cwiseAbs().maxCoeff();
        const real_t zeroTol = 1e-10;
        index_t nPos = 0, nNeg = 0, nZero = 0;
        std::vector<real_t> ratios(ev_.rows());
        for (index_t i = 0; i != ev_.rows(); ++i)
        {
            ratios[i] = math::abs(ev_[i]) / lmaxAbs;
            if (math::abs(ev_[i]) < zeroTol * lmaxAbs) ++nZero;
            else if (ev_[i] >= zeroTol * lmaxAbs) ++nPos;
            else ++nNeg;
        }
        std::sort(ratios.begin(), ratios.end());
        result.nPos = nPos; result.nNeg = nNeg; result.nZero = nZero;
        for (index_t i = 0; i != 5 && i != (index_t)ratios.size(); ++i)
            result.smallestRatios[i] = ratios[i];
        result.haveSpectrum = true;
    }

    // 15d. Solve on K (LU: the system is indefinite, so neither CG nor
    //      pivotless LDLT applies).
    gsSparseSolver<real_t>::LU solver;
    solver.compute(K);
    solVector = solver.solve(A.rhs());
    result.finite = solVector.allFinite();

    if (result.finite)
    {
        auto u_ex = ev.getVariable(u_exact, G);
        auto p_ex = ev.getVariable(p_exact, G);
        result.l2u = math::sqrt(ev.integral((u_ex - v_sol).sqNorm() * meas(G)));
        result.h1u = result.l2u + math::sqrt(ev.integral((igrad(u_ex) - igrad(v_sol, G)).sqNorm() * meas(G)));

        const real_t area = ev.integral(meas(G));
        const real_t I1 = ev.integral((p_ex - p_sol) * meas(G));
        const real_t I2 = ev.integral((p_ex - p_sol).sqNorm() * meas(G));
        result.l2p = math::sqrt(I2);
        const real_t cshift = I1 / area;
        result.l2pShift = math::sqrt(math::max((real_t)0.0, I2 - cshift*cshift*area));

        if (doPlot)
        {
            ev.options().setSwitch("plot.elements", true);
            ev.options().setInt("plot.npts", 100000);
            ev.writeParaview(v_sol, G, outPath + "/immersed_stokes_skeleton_velocity");
            ev.writeParaview(p_sol, G, outPath + "/immersed_stokes_skeleton_pressure");
            gsMesh<> mesh(dbasis.basis(0));
            gsWriteParaview(mesh, outPath + "/immersed_stokes_skeleton_background_mesh");
        }
    }

    return result;
}

namespace {

void runRates(index_t degree, index_t numRefine, real_t beta, real_t gamma, real_t gammaT,
              index_t quRule, index_t quA, bool noGhost, bool noSkeleton,
              const std::string & pressureMode, real_t L, bool plot, const std::string & outPath)
{
    gsFunctionExpr<real_t> u_exact(
        "x^2*y^4*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 18*x^2*y^2 - 85*x^2 + 13*y^4 - 153*y^2 + 80)/1000000",
        "-x*y^5*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 6*x^2*y^2 - 51*x^2 + y^4 - 17*y^2 + 16)/500000", 2);
    gsFunctionExpr<real_t> p_exact(
        "-x*y*(x - y)*(x + y)*(x^2 + y^2 - 16)^2*(x^2 + y^2 - 1)^2/10000000*exp(14/sqrt(x^2 + y^2))", 2);
    gsFunctionExpr<real_t> f_rhs(
        "(-y^2*(30*x^10 + 645*x^8*y^2 - 1020*x^8 + 2296*x^6*y^4 - 15470*x^6*y^2 + 9630*x^6 + 2790*x^4*y^6 - 36414*x^4*y^4 + 91485*x^4*y^2 - 16320*x^4 + 1122*x^2*y^8 - 22338*x^2*y^6 + 107856*x^2*y^4 - 73440*x^2*y^2 + 7680*x^2 + 13*y^10 - 374*y^8 + 2889*y^6 - 3808*y^4 + 1280*y^2)/500000) + exp(14/sqrt(x^2 + y^2))*(-y*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(11*x^8*sqrt(x^2 + y^2) - 14*x^8 + 16*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 119*x^6*sqrt(x^2 + y^2) + 238*x^6 - 2*x^4*y^4*sqrt(x^2 + y^2) + 14*x^4*y^4 - 85*x^4*y^2*sqrt(x^2 + y^2) + 48*x^4*sqrt(x^2 + y^2) - 224*x^4 - 8*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 51*x^2*y^4*sqrt(x^2 + y^2) - 238*x^2*y^4 + 32*x^2*y^2*sqrt(x^2 + y^2) + 224*x^2*y^2 - y^8*sqrt(x^2 + y^2) + 17*y^6*sqrt(x^2 + y^2) - 16*y^4*sqrt(x^2 + y^2))/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))",
        "(x*y^3*(25*x^8 + 258*x^6*y^2 - 680*x^6 + 492*x^4*y^4 - 4641*x^4*y^2 + 4815*x^4 + 310*x^2*y^6 - 5202*x^2*y^4 + 18297*x^2*y^2 - 5440*x^2 + 51*y^8 - 1241*y^6 + 7704*y^4 - 7344*y^2 + 1280)/125000) + exp(14/sqrt(x^2 + y^2))*(-x*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(x^8*sqrt(x^2 + y^2) + 8*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 17*x^6*sqrt(x^2 + y^2) + 2*x^4*y^4*sqrt(x^2 + y^2) - 14*x^4*y^4 - 51*x^4*y^2*sqrt(x^2 + y^2) + 238*x^4*y^2 + 16*x^4*sqrt(x^2 + y^2) - 16*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 85*x^2*y^4*sqrt(x^2 + y^2) - 32*x^2*y^2*sqrt(x^2 + y^2) - 224*x^2*y^2 - 11*y^8*sqrt(x^2 + y^2) + 14*y^8 + 119*y^6*sqrt(x^2 + y^2) - 238*y^6 - 48*y^4*sqrt(x^2 + y^2) + 224*y^4)/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))", 2);

    gsFunctionExpr<real_t> phi("(x^2 + y^2 - 1)*(x^2 + y^2 - 16)", 2);
    // Unit outward normal to the level set; only valid on the zero set {phi = 0}.
    gsFunctionExpr<real_t> nImm("x*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))",
                                "y*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))", 2);

    gsVector<real_t> l2u(numRefine+1), h1u(numRefine+1), l2p(numRefine+1), l2pS(numRefine+1);

    gsInfo << "\n=== rates  degree=" << degree << "  beta=" << beta
           << "  gamma=" << gamma << "  gammaT=" << gammaT
           << "  pressureMode=" << pressureMode << " ===\n";
    gsInfo << std::right << std::setw(4) << "r" << std::setw(6) << "N" << std::setw(9) << "h"
           << std::setw(8) << "dofs" << std::setw(13) << "L2-u" << std::setw(13) << "H1-u"
           << std::setw(13) << "L2-p" << std::setw(13) << "L2-p*" << "\n";

    for (index_t r = 0; r <= numRefine; ++r)
    {
        const index_t N = 8 << r;
        const bool doPlot = (plot && r == numRefine);
        RunResult res = solveOne(degree, N, L, phi, nImm, u_exact, p_exact, f_rhs,
                                  beta, gamma, gammaT, quRule, quA, noGhost, noSkeleton,
                                  pressureMode, /*computeSpectrum*/ false, doPlot, outPath);
        l2u[r] = res.l2u; h1u[r] = res.h1u;
        l2p[r] = res.l2p; l2pS[r] = res.l2pShift;

        gsInfo << std::setw(4) << r << std::setw(6) << N << std::setw(9) << fmtSci(res.h,2)
               << std::setw(8) << res.ndofs
               << std::setw(13) << fmtSci(res.l2u,3) << std::setw(13) << fmtSci(res.h1u,3)
               << std::setw(13) << fmtSci(res.l2p,3) << std::setw(13) << fmtSci(res.l2pShift,3) << "\n";
    }

    gsInfo << "Pressure rate criterion uses the shifted L2-p* in both pressure modes: under "
              "pressureMode=none the constant pressure mode is an exact discrete kernel vector "
              "(Dirichlet velocity on the whole boundary makes B^T*1 = 0 and dnk(const) = 0), so "
              "the constant is not recoverable and raw L2-p is reported for information only.\n";
    gsInfo << faceShiftLegend << "\n";

    if (numRefine > 0)
    {
        gsInfo << "\nEoC (L2-u): " << std::fixed << std::setprecision(2)
               << (l2u.head(numRefine).array() / l2u.tail(numRefine).array()).log().transpose() / std::log(2.0) << "\n";
        gsInfo << "EoC (H1-u): "
               << (h1u.head(numRefine).array() / h1u.tail(numRefine).array()).log().transpose() / std::log(2.0) << "\n";
        gsInfo << "EoC (L2-p): "
               << (l2p.head(numRefine).array() / l2p.tail(numRefine).array()).log().transpose() / std::log(2.0) << "\n";
        gsInfo << "EoC (L2-p*): "
               << (l2pS.head(numRefine).array() / l2pS.tail(numRefine).array()).log().transpose() / std::log(2.0) << "\n";
    }
}

void runGammaSweep(bool tildeSweep, index_t degree, real_t beta, real_t gammaBase, real_t gammaTBase,
                    index_t quRule, index_t quA, bool noGhost, bool noSkeleton,
                    const std::string & pressureMode, real_t L)
{
    gsFunctionExpr<real_t> u_exact(
        "x^2*y^4*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 18*x^2*y^2 - 85*x^2 + 13*y^4 - 153*y^2 + 80)/1000000",
        "-x*y^5*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 6*x^2*y^2 - 51*x^2 + y^4 - 17*y^2 + 16)/500000", 2);
    gsFunctionExpr<real_t> p_exact(
        "-x*y*(x - y)*(x + y)*(x^2 + y^2 - 16)^2*(x^2 + y^2 - 1)^2/10000000*exp(14/sqrt(x^2 + y^2))", 2);
    gsFunctionExpr<real_t> f_rhs(
        "(-y^2*(30*x^10 + 645*x^8*y^2 - 1020*x^8 + 2296*x^6*y^4 - 15470*x^6*y^2 + 9630*x^6 + 2790*x^4*y^6 - 36414*x^4*y^4 + 91485*x^4*y^2 - 16320*x^4 + 1122*x^2*y^8 - 22338*x^2*y^6 + 107856*x^2*y^4 - 73440*x^2*y^2 + 7680*x^2 + 13*y^10 - 374*y^8 + 2889*y^6 - 3808*y^4 + 1280*y^2)/500000) + exp(14/sqrt(x^2 + y^2))*(-y*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(11*x^8*sqrt(x^2 + y^2) - 14*x^8 + 16*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 119*x^6*sqrt(x^2 + y^2) + 238*x^6 - 2*x^4*y^4*sqrt(x^2 + y^2) + 14*x^4*y^4 - 85*x^4*y^2*sqrt(x^2 + y^2) + 48*x^4*sqrt(x^2 + y^2) - 224*x^4 - 8*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 51*x^2*y^4*sqrt(x^2 + y^2) - 238*x^2*y^4 + 32*x^2*y^2*sqrt(x^2 + y^2) + 224*x^2*y^2 - y^8*sqrt(x^2 + y^2) + 17*y^6*sqrt(x^2 + y^2) - 16*y^4*sqrt(x^2 + y^2))/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))",
        "(x*y^3*(25*x^8 + 258*x^6*y^2 - 680*x^6 + 492*x^4*y^4 - 4641*x^4*y^2 + 4815*x^4 + 310*x^2*y^6 - 5202*x^2*y^4 + 18297*x^2*y^2 - 5440*x^2 + 51*y^8 - 1241*y^6 + 7704*y^4 - 7344*y^2 + 1280)/125000) + exp(14/sqrt(x^2 + y^2))*(-x*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(x^8*sqrt(x^2 + y^2) + 8*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 17*x^6*sqrt(x^2 + y^2) + 2*x^4*y^4*sqrt(x^2 + y^2) - 14*x^4*y^4 - 51*x^4*y^2*sqrt(x^2 + y^2) + 238*x^4*y^2 + 16*x^4*sqrt(x^2 + y^2) - 16*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 85*x^2*y^4*sqrt(x^2 + y^2) - 32*x^2*y^2*sqrt(x^2 + y^2) - 224*x^2*y^2 - 11*y^8*sqrt(x^2 + y^2) + 14*y^8 + 119*y^6*sqrt(x^2 + y^2) - 238*y^6 - 48*y^4*sqrt(x^2 + y^2) + 224*y^4)/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))", 2);

    gsFunctionExpr<real_t> phi("(x^2 + y^2 - 1)*(x^2 + y^2 - 16)", 2);
    gsFunctionExpr<real_t> nImm("x*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))",
                                "y*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))", 2);

    const real_t base = tildeSweep ? gammaTBase : gammaBase;
    const real_t values[5] = { base*1e-2, base*1e-1, base, base*1e1, base*1e2 };
    const index_t N = 8;

    gsInfo << "\n=== " << (tildeSweep ? "gammaTSweep" : "gammaSweep") << "  degree=" << degree
           << "  N=" << N << " ===\n";
    gsInfo << std::right << std::setw(13) << (tildeSweep ? "gammaT" : "gamma")
           << std::setw(8) << "dofs" << std::setw(13) << "L2-u" << std::setw(13) << "H1-u"
           << std::setw(13) << "L2-p" << std::setw(8) << "nZero" << std::setw(15) << "lminAbs/lmax" << "\n";

    for (index_t j = 0; j != 5; ++j)
    {
        const real_t g  = tildeSweep ? gammaBase  : values[j];
        const real_t gt = tildeSweep ? values[j]  : gammaTBase;
        RunResult res = solveOne(degree, N, L, phi, nImm, u_exact, p_exact, f_rhs,
                                  beta, g, gt, quRule, quA, noGhost, noSkeleton,
                                  pressureMode, /*computeSpectrum*/ true, false, "");
        gsInfo << std::setw(13) << fmtSci(values[j],2) << std::setw(8) << res.ndofs
               << std::setw(13) << fmtSci(res.l2u,3) << std::setw(13) << fmtSci(res.h1u,3)
               << std::setw(13) << fmtSci(res.l2p,3) << std::setw(8) << res.nZero
               << std::setw(15) << fmtSci(res.smallestRatios[0],3) << "\n";
    }
}

void runInertia(index_t degree, real_t beta, real_t gamma, real_t gammaT,
                 index_t quRule, index_t quA, bool noGhost, const std::string & pressureMode, real_t L)
{
    gsFunctionExpr<real_t> u_exact(
        "x^2*y^4*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 18*x^2*y^2 - 85*x^2 + 13*y^4 - 153*y^2 + 80)/1000000",
        "-x*y^5*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 6*x^2*y^2 - 51*x^2 + y^4 - 17*y^2 + 16)/500000", 2);
    gsFunctionExpr<real_t> p_exact(
        "-x*y*(x - y)*(x + y)*(x^2 + y^2 - 16)^2*(x^2 + y^2 - 1)^2/10000000*exp(14/sqrt(x^2 + y^2))", 2);
    gsFunctionExpr<real_t> f_rhs(
        "(-y^2*(30*x^10 + 645*x^8*y^2 - 1020*x^8 + 2296*x^6*y^4 - 15470*x^6*y^2 + 9630*x^6 + 2790*x^4*y^6 - 36414*x^4*y^4 + 91485*x^4*y^2 - 16320*x^4 + 1122*x^2*y^8 - 22338*x^2*y^6 + 107856*x^2*y^4 - 73440*x^2*y^2 + 7680*x^2 + 13*y^10 - 374*y^8 + 2889*y^6 - 3808*y^4 + 1280*y^2)/500000) + exp(14/sqrt(x^2 + y^2))*(-y*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(11*x^8*sqrt(x^2 + y^2) - 14*x^8 + 16*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 119*x^6*sqrt(x^2 + y^2) + 238*x^6 - 2*x^4*y^4*sqrt(x^2 + y^2) + 14*x^4*y^4 - 85*x^4*y^2*sqrt(x^2 + y^2) + 48*x^4*sqrt(x^2 + y^2) - 224*x^4 - 8*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 51*x^2*y^4*sqrt(x^2 + y^2) - 238*x^2*y^4 + 32*x^2*y^2*sqrt(x^2 + y^2) + 224*x^2*y^2 - y^8*sqrt(x^2 + y^2) + 17*y^6*sqrt(x^2 + y^2) - 16*y^4*sqrt(x^2 + y^2))/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))",
        "(x*y^3*(25*x^8 + 258*x^6*y^2 - 680*x^6 + 492*x^4*y^4 - 4641*x^4*y^2 + 4815*x^4 + 310*x^2*y^6 - 5202*x^2*y^4 + 18297*x^2*y^2 - 5440*x^2 + 51*y^8 - 1241*y^6 + 7704*y^4 - 7344*y^2 + 1280)/125000) + exp(14/sqrt(x^2 + y^2))*(-x*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(x^8*sqrt(x^2 + y^2) + 8*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 17*x^6*sqrt(x^2 + y^2) + 2*x^4*y^4*sqrt(x^2 + y^2) - 14*x^4*y^4 - 51*x^4*y^2*sqrt(x^2 + y^2) + 238*x^4*y^2 + 16*x^4*sqrt(x^2 + y^2) - 16*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 85*x^2*y^4*sqrt(x^2 + y^2) - 32*x^2*y^2*sqrt(x^2 + y^2) - 224*x^2*y^2 - 11*y^8*sqrt(x^2 + y^2) + 14*y^8 + 119*y^6*sqrt(x^2 + y^2) - 238*y^6 - 48*y^4*sqrt(x^2 + y^2) + 224*y^4)/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))", 2);

    gsFunctionExpr<real_t> phi("(x^2 + y^2 - 1)*(x^2 + y^2 - 16)", 2);
    gsFunctionExpr<real_t> nImm("x*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))",
                                "y*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))", 2);

    const index_t N = 8;
    gsInfo << "\n=== inertia  degree=" << degree << "  N=" << N
           << "  (--noSkeleton is ignored here; --noGhost is honoured) ===\n";
    gsInfo << std::right << std::setw(14) << "config" << std::setw(6) << "nPos" << std::setw(6) << "nNeg"
           << std::setw(6) << "nZero" << std::setw(10) << "exp(nV)" << std::setw(8) << "exp(nP)"
           << "   five smallest |lambda|/lmax\n";

    for (int j = 0; j != 2; ++j)
    {
        const bool skelOff = (0 == j);
        RunResult res = solveOne(degree, N, L, phi, nImm, u_exact, p_exact, f_rhs,
                                  beta, gamma, gammaT, quRule, quA, noGhost, skelOff,
                                  pressureMode, /*computeSpectrum*/ true, false, "");
        gsInfo << std::setw(14) << (skelOff ? "skeleton off" : "skeleton on")
               << std::setw(6) << res.nPos << std::setw(6) << res.nNeg << std::setw(6) << res.nZero
               << std::setw(10) << res.nV << std::setw(8) << res.nP << "   ";
        for (index_t i = 0; i != 5; ++i) gsInfo << fmtSci(res.smallestRatios[i],3) << " ";
        gsInfo << "\n";
    }
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    index_t degree      = 1;
    index_t numRefine   = 2;
    real_t  gamma        = -1;
    real_t  gammaT       = -1;
    real_t  beta         = -1;
    index_t quRule       = gsQuadrature::AlgoimRule;
    index_t quA          = -1;
    bool    noGhost       = false;
    bool    noSkeleton    = false;
    real_t  sliverDelta  = 0;
    std::string pressureMode = "pin";
    std::string study        = "rates";
    bool    plot          = false;
    std::string outFolder = "output_immersed_stokes_skeleton";

    gsCmdLine cmd("Skeleton/ghost-stabilized immersed Stokes on the quarter annulus "
                  "1 < r < 4, x > 0, y > 0.");
    cmd.addInt   ("p", "degree",       "Spline degree of velocity/pressure (1 or 2)", degree);
    cmd.addInt   ("r", "uniformRefine","Refinement levels of the rates study (r = 0..r)", numRefine);
    cmd.addReal  ("",  "gamma",       "Skeleton (pressure) stabilization coefficient (<0: auto)", gamma);
    cmd.addReal  ("",  "gammaT",      "Ghost (velocity) stabilization coefficient (<0: auto)", gammaT);
    cmd.addReal  ("",  "beta",        "Nitsche penalty coefficient (<=0: auto 6*(k+1)^2)", beta);
    cmd.addInt   ("q", "quRule",      "Immersed volume rule (11=CutCell, 12=Algoim; 14=MomentFitting "
                                       "is rejected, see below)", quRule);
    cmd.addInt   ("",  "quA",         "Quadrature order parameter quA (<=0: library default)", quA);
    cmd.addSwitch("",  "noGhost",     "Disable ghost-penalty (velocity) stabilization", noGhost);
    cmd.addSwitch("",  "noSkeleton",  "Disable skeleton (pressure) stabilization", noSkeleton);
    cmd.addReal  ("",  "sliverDelta", "Sliver depth d in [0,1): places a mesh line at x=4-d*h", sliverDelta);
    cmd.addString("",  "pressureMode","Pressure treatment: none | pin (lagrange is out of scope)", pressureMode);
    cmd.addString("",  "study",       "Study: rates | gammaSweep | gammaTSweep | inertia", study);
    cmd.addSwitch("plot", "Create ParaView output at the finest rates level", plot);
    cmd.addString("o", "output",      "Output folder", outFolder);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if ("rates" != study && "gammaSweep" != study && "gammaTSweep" != study && "inertia" != study)
    { gsWarn << "--study must be rates, gammaSweep, gammaTSweep or inertia\n"; return EXIT_FAILURE; }

    if ("none" != pressureMode && "pin" != pressureMode)
    { gsWarn << "--pressureMode must be none or pin (lagrange is deliberately out of scope)\n"; return EXIT_FAILURE; }

    if (1 != degree && 2 != degree)
    { gsWarn << "-p/--degree must be 1 or 2\n"; return EXIT_FAILURE; }

    if (numRefine < 0)
    { gsWarn << "-r/--uniformRefine must be >= 0\n"; return EXIT_FAILURE; }
    if (numRefine > 4)
        gsWarn << "-r/--uniformRefine = " << numRefine << " is large; the dense inertia diagnostic "
                  "is only ever run at one pinned small level, never inside this loop.\n";

    if (quRule != gsQuadrature::AlgoimRule && quRule != gsQuadrature::CutCellRule
        && quRule != gsQuadrature::MomentFittingRule)
    { gsWarn << "-q/--quRule must be an immersed volume rule (11, 12 or 14)\n"; return EXIT_FAILURE; }
    if (gsQuadrature::MomentFittingRule == quRule)
    {
        // gsQuadrature::makeMomentFittingPtr GISMO_ENSUREs quDim < 0 (it
        // compresses points onto a volume tensor grid, which is meaningless
        // on the lower-dimensional zero level set). This driver's Nitsche
        // terms unconditionally assemble with quDim == 2, so quRule ==
        // MomentFittingRule can never complete a run here.
        gsWarn << "-q/--quRule 14 (MomentFittingRule) cannot be used: it refuses surface "
                  "quadrature (quDim >= 0), and every Nitsche term in this driver assembles "
                  "with quDim == 2. Use 11 (CutCellRule) or 12 (AlgoimRule).\n";
        return EXIT_FAILURE;
    }
    index_t quAi = -1;
    if (quA > 0) quAi = static_cast<index_t>(quA);

    if (sliverDelta < 0 || sliverDelta >= 1)
    { gsWarn << "--sliverDelta must lie in [0,1)\n"; return EXIT_FAILURE; }
    if (sliverDelta > 0 && "rates" == study)
    {
        gsWarn << "--sliverDelta > 0 changes h (see file: L = 4*N/(j+d)), so it is not compatible "
                  "with --study rates, which needs a fixed geometry across refinement levels.\n";
        return EXIT_FAILURE;
    }

    const real_t betaEff   = (beta <= 0) ? 6.0 * (degree + 1) * (degree + 1) : beta;
    const real_t gammaEff  = (gamma < 0) ? ((1 == degree) ? 10.0 : 0.1) : gamma;
    const real_t gammaTEff = (gammaT < 0) ? math::pow(10.0, -(degree + 1)) : gammaT;

    real_t L = 4.2;
    if (sliverDelta > 0)
    {
        const index_t N = 8;
        const index_t j = (index_t)std::floor(4.0 * N / 4.2);
        L = 4.0 * N / ((real_t)j + sliverDelta);
        const real_t h = L / N;

        gsFunctionExpr<real_t> phiCheck("(x^2 + y^2 - 1)*(x^2 + y^2 - 16)", 2);
        gsMultiPatch<real_t> mpCheck(*gsNurbsCreator<real_t>::BSplineRectangleWithPara(0.0, 0.0, L, L));
        gsMultiBasis<real_t> dbasisCheck(mpCheck, true);
        dbasisCheck.setDegree(degree);
        dbasisCheck.uniformRefine(N - 1);
        gsTensorBSplineBasis<2,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasisCheck.basis(0));
        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > trCheck =
            memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phiCheck, *tbsPtr));
        const size_t sliverFlat = (size_t)j;   // cell (j, 0)
        const short_t sign = sliverSign(*trCheck, sliverFlat, L, N);
        if (0 != sign)
        {
            gsWarn << "--sliverDelta " << sliverDelta << ": the engineered sliver cell (" << j
                   << ",0) was not classified cut by the trimmed domain (sign=" << (int)sign
                   << "); reporting numbers for a sliver the classifier never saw would be "
                      "meaningless. L=" << L << " h=" << h << "\n";
            return EXIT_FAILURE;
        }
        gsInfo << "sliverDelta=" << sliverDelta << ": effective L=" << L << " h=" << h << "\n";
    }

    // Manufactured solution for the self-check and traction oracle.
    gsFunctionExpr<real_t> u_exact(
        "x^2*y^4*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 18*x^2*y^2 - 85*x^2 + 13*y^4 - 153*y^2 + 80)/1000000",
        "-x*y^5*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(5*x^4 + 6*x^2*y^2 - 51*x^2 + y^4 - 17*y^2 + 16)/500000", 2);
    gsFunctionExpr<real_t> p_exact(
        "-x*y*(x - y)*(x + y)*(x^2 + y^2 - 16)^2*(x^2 + y^2 - 1)^2/10000000*exp(14/sqrt(x^2 + y^2))", 2);
    gsFunctionExpr<real_t> f_rhs(
        "(-y^2*(30*x^10 + 645*x^8*y^2 - 1020*x^8 + 2296*x^6*y^4 - 15470*x^6*y^2 + 9630*x^6 + 2790*x^4*y^6 - 36414*x^4*y^4 + 91485*x^4*y^2 - 16320*x^4 + 1122*x^2*y^8 - 22338*x^2*y^6 + 107856*x^2*y^4 - 73440*x^2*y^2 + 7680*x^2 + 13*y^10 - 374*y^8 + 2889*y^6 - 3808*y^4 + 1280*y^2)/500000) + exp(14/sqrt(x^2 + y^2))*(-y*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(11*x^8*sqrt(x^2 + y^2) - 14*x^8 + 16*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 119*x^6*sqrt(x^2 + y^2) + 238*x^6 - 2*x^4*y^4*sqrt(x^2 + y^2) + 14*x^4*y^4 - 85*x^4*y^2*sqrt(x^2 + y^2) + 48*x^4*sqrt(x^2 + y^2) - 224*x^4 - 8*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 51*x^2*y^4*sqrt(x^2 + y^2) - 238*x^2*y^4 + 32*x^2*y^2*sqrt(x^2 + y^2) + 224*x^2*y^2 - y^8*sqrt(x^2 + y^2) + 17*y^6*sqrt(x^2 + y^2) - 16*y^4*sqrt(x^2 + y^2))/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))",
        "(x*y^3*(25*x^8 + 258*x^6*y^2 - 680*x^6 + 492*x^4*y^4 - 4641*x^4*y^2 + 4815*x^4 + 310*x^2*y^6 - 5202*x^2*y^4 + 18297*x^2*y^2 - 5440*x^2 + 51*y^8 - 1241*y^6 + 7704*y^4 - 7344*y^2 + 1280)/125000) + exp(14/sqrt(x^2 + y^2))*(-x*(x^2 + y^2 - 16)*(x^2 + y^2 - 1)*(x^8*sqrt(x^2 + y^2) + 8*x^6*y^2*sqrt(x^2 + y^2) - 14*x^6*y^2 - 17*x^6*sqrt(x^2 + y^2) + 2*x^4*y^4*sqrt(x^2 + y^2) - 14*x^4*y^4 - 51*x^4*y^2*sqrt(x^2 + y^2) + 238*x^4*y^2 + 16*x^4*sqrt(x^2 + y^2) - 16*x^2*y^6*sqrt(x^2 + y^2) + 14*x^2*y^6 + 85*x^2*y^4*sqrt(x^2 + y^2) - 32*x^2*y^2*sqrt(x^2 + y^2) - 224*x^2*y^2 - 11*y^8*sqrt(x^2 + y^2) + 14*y^8 + 119*y^6*sqrt(x^2 + y^2) - 238*y^6 - 48*y^4*sqrt(x^2 + y^2) + 224*y^4)/10000000)/((x^2 + y^2)*sqrt(x^2 + y^2))", 2);
    gsFunctionExpr<real_t> phi("(x^2 + y^2 - 1)*(x^2 + y^2 - 16)", 2);
    gsFunctionExpr<real_t> nImm("x*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))",
                                "y*(2*(x^2 + y^2) - 17)/(15*sqrt(x^2 + y^2))", 2);

    if (!selfCheck(u_exact, p_exact, f_rhs))
    { gsWarn << "selfCheck FAILED: the manufactured-solution strings do not match the reference "
                "table; refusing to assemble.\n"; return EXIT_FAILURE; }

    if (!nitscheTractionCheck(phi, nImm))
    { gsWarn << "nitscheTractionCheck FAILED: the Nitsche traction matrix terms do not reproduce "
                "a numerical evaluation of (2*mu*grad^s u).n to 1e-10 relative.\n"; return EXIT_FAILURE; }

    const index_t N0 = 8;
    gsInfo << "k=" << degree << " N0=" << N0 << " h0=" << fmtSci(L/N0,4)
           << " beta=" << betaEff << " gamma=" << gammaEff << " gammaT=" << gammaTEff
           << " quRule=" << quRule << " quA=" << (quAi>0?quAi:0)
           << " pressureMode=" << pressureMode << " noGhost=" << noGhost << " noSkeleton=" << noSkeleton
           << "\n" << faceShiftLegend << "\n";

    std::string outPath = outFolder;
    if (outPath.empty()) outPath = "output_immersed_stokes_skeleton";
    const bool isAbsolutePath = (!outPath.empty() && outPath[0] == '/');
    if (!isAbsolutePath) outPath = gsFileManager::getCurrentPath() + "/" + outPath;
    const std::string out = gsFileManager::getCanonicRepresentation(outPath);
    if (plot) gsFileManager::mkdir(out);

    if ("rates" == study)
        runRates(degree, numRefine, betaEff, gammaEff, gammaTEff, quRule, quAi,
                 noGhost, noSkeleton, pressureMode, L, plot, out);
    else if ("gammaSweep" == study)
        runGammaSweep(false, degree, betaEff, gammaEff, gammaTEff, quRule, quAi,
                      noGhost, noSkeleton, pressureMode, L);
    else if ("gammaTSweep" == study)
        runGammaSweep(true, degree, betaEff, gammaEff, gammaTEff, quRule, quAi,
                      noGhost, noSkeleton, pressureMode, L);
    else if ("inertia" == study)
        runInertia(degree, betaEff, gammaEff, gammaTEff, quRule, quAi, noGhost, pressureMode, L);

    return EXIT_SUCCESS;
}
