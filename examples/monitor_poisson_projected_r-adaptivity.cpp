/** @file monitor_poisson_projected_r-adaptivity.cpp

    @brief Quasi-optimal r-adaptivity by projection: r-adapt a composed map
    G = geom∘σ (monitor-based relocation of σ), then L2-project G onto a
    plain spline basis and solve the Poisson problem on the projected
    geometry with STANDARD quadrature.

    Rationale (see RH_FITTING_MATH.md and the moving-mesh literature): the
    r-adapted map G is smooth — only its parameterization is graded — so
    projecting it onto a plain basis incurs a full-order-small geometric
    consistency error, while assembly on the projected geometry avoids the
    composed-map quadrature machinery entirely (union integration mesh,
    SameElement=false). The price is a "quasi"-optimal parameterization and
    a possible loss of injectivity (checked per element here).

    Three Poisson solves are compared per refinement level:
      U  uniform     : identity geometry, plain basis, standard quadrature;
      C  composed    : map G = geom∘σ, plain basis, union integration basis
                       + SameElement=false (the expensive reference);
      P  projected   : map ΠG on a plain basis (analysis basis + `-k` extra
                       refinements), standard quadrature.
    All paths use the SAME discretization space dimension at k=0; σ is
    relocated ONCE (it does not depend on the analysis basis) and reused on
    every level.

    FINDINGS (2026-07-24, tanh-front test, theta=10..100, 3 levels):
    - P is not a lossy compromise: it beats the composed reference C by
      3-4x in L2 error at matched DoFs (and C can even fall behind the
      uniform U at strong grading). Mechanism: G = geom∘σ carries
      second-derivative kinks along σ's knot lines, which cut through the
      INTERIOR of analysis elements; the pulled-back solution u∘G then has
      reduced smoothness inside elements and the approximation constant
      degrades. Projection onto the analysis basis relocates the map's
      kinks to element boundaries and repairs the approximation property.
    - Nested corollary, PRECISE FORM: exactness of the projection needs
      G itself to live in a spline space, which holds only for
      identity/affine templates (then G = σ; confirmed to machine
      precision with -R 7 on dyadic meshes: C = P = standard assembly to
      all digits, no L2 solve needed — pure knot insertion). For a CURVED
      template, G = geom∘σ has degree p_G*p_σ on curved pullback pieces
      and is NEVER exactly representable (--geomWave 0.1 measures
      J-distance ~1e-2 with nested knots) — BUT the loss is negligible:
      P ≈ C to within 8% (3.02e-5 vs 3.27e-5 at 4096 DoFs), both 4.5x
      better than U. Nesting σ's knots also helps the COMPOSED path
      itself (C: 3.3e-5 nested vs 1.15e-4 non-nested at equal DoFs),
      because the strong (grading-carrying) kinks of the map belong to σ
      and nesting aligns them with element boundaries; geom's curved
      pullback kinks are C^{p_G-1}-weak and ungraded. Two design rules:
      (a) always nest σ's knots in the analysis mesh; (b) project and
      discard the composition — costs <=8% accuracy, removes all composed
      machinery downstream.
    - Diagnostics S (σ as a plain spline map = C to all digits) verified
      the gsComposedGeometry evaluation chain; S2 (ignoring σ's knots in
      the integration mesh) loses 1-2 orders of magnitude — the union
      integration basis is load-bearing for composed assembly.
    - The projErr returned by gsL2Projection is a SQUARED quantity; the
      driver prints true sup-norm map/Jacobian distances separately.

    Notes:
    - The template geometry is the identity on [0,1]^2, hence G = σ exactly
      (a plain spline on σ's knot mesh). The composed path still exercises
      the full composed-assembly machinery. For curved template geometries
      G is genuinely not a spline; nesting σ's knots then still aligns σ's
      kinks with element boundaries, but the geometry-knot pullbacks stay
      curved and a (small, smooth-map) projection error returns.
    - The ValueBased monitor ω = 1/sqrt(1+θ f) CONTRACTS σ (small cells)
      where f is large; f = |∇u_ex| (normalized to [0,1]) therefore
      concentrates cells at the solution front, which is what a Poisson
      solve wants (opposite of the fitting driver, which inverts f).
    - θ (-T/--theta) and --steepness do DIFFERENT things. ω = 1/sqrt(1+θ f)
      is bounded to [1/sqrt(1+θ), 1]: θ ALONE sets how small a cell can get,
      i.e. the reachable min det J_σ (measured: θ=1 → 0.726, θ=10 → 0.416,
      θ=100 → folds at the default penalty). --steepness applies a
      sharpening exponent p to the monitor, f ← max(f,0)^p, which NARROWS
      the high-weight band around the front but leaves f's endpoints (0 and
      1) — and therefore ω's bounds — untouched. So: θ moves min det J,
      steepness narrows the band. Reach for θ first.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsUtils/gsL2Projection.h>
#include <gsUtils/gsQuasiInterpolate.h>

#include <fstream>
#include <iomanip>
#include <sstream>

using namespace gismo;

namespace {

struct PoissonRun
{
    index_t dofs;
    real_t  tAssemble, tSolve, l2err, h1err;
    gsMultiPatch<real_t> sol; // parametric solution spline (plain basis)
};

// Per-level diagnostics of the projected map; reported on the P row.
struct ProjDiag
{
    real_t projErr    = -1; // SQUARED, weighted -- see the header comment
    real_t minDetJ    =  0; // of PiG, on its own basis
    real_t mapDist    = -1; // max|G - PiG| on a 201^2 grid
    real_t jacDist    = -1; // max|J_G - J_PiG| on the same grid
    real_t bdrMapDist = -1; // max|G - PiG| on the four edges
    real_t minDetJG   = -1; // of G, on the UNION mesh
};

// L2 error and H1 error of one solution field. The H1 entry follows the
// convention of examples/poisson2_example.cpp:180-181 -- h1 = l2 + the H1
// SEMINORM of the error, NOT sqrt(l2^2 + seminorm^2). That is the established
// convention in this codebase; keep it so the numbers are comparable with the
// library's own Poisson tutorial.
template <class T>
struct ErrorNorms
{
    T l2, h1;
};

// L2 and H1 errors of a parametric solution spline against the physical exact
// solution ms, integrated through the map mp on the elements of ib. Both norms
// share ONE evaluator, one integration mesh and one quadrature rule.
template <class T>
ErrorNorms<T> evalErrorNorms(const gsMultiPatch<T> & mp,
                             const gsMultiBasis<T> & ib,
                             const gsMultiPatch<T> & sol,
                             const gsFunctionExpr<T> & ms,
                             bool sameElement,
                             T errQuA = 4.0, index_t errQuB = 2)
{
    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(ib);
    ev.options().setReal("quA", errQuA);
    ev.options().setInt ("quB", errQuB);
    // a standalone evaluator has no SameElement switch in its list yet
    ev.options().addSwitch("SameElement", "", sameElement);
    auto G = ev.getMap(mp);
    auto u_sol = ev.getVariable(sol);   // PARAMETRIC: no map attached, hence
                                        // the explicit chain rule igrad(.,G)
    auto u_ex  = ev.getVariable(ms, G); // PHYSICAL: composed with G already
    ErrorNorms<T> e;
    e.l2 = math::sqrt(ev.integral((u_sol - u_ex).sqNorm() * meas(G)));
    // h1 = l2 + H1-SEMINORM of the error: the convention of
    // examples/poisson2_example.cpp:180-181, kept deliberately so the numbers
    // are comparable with the library's own Poisson tutorial.
    e.h1 = e.l2 + math::sqrt(ev.integral(
               (igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));
    return e;
}

// Solve -Δu = f (physical) with Dirichlet data u = ms (physical) on the
// geometry mp, discretization basis mb, integration elements ib. The L2
// error is computed against ms with elevated quadrature on the same mesh.
template <class T>
PoissonRun solvePoisson(const gsMultiPatch<T> & mp,
                        const gsMultiBasis<T> & mb,
                        const gsMultiBasis<T> & ib,
                        const gsFunctionExpr<T> & f,
                        const gsFunctionExpr<T> & ms,
                        bool sameElement,
                        T quA = 2.0, index_t quB = 2,
                        T errQuA = 4.0, index_t errQuB = 2)
{
    gsExprAssembler<T> A(1, 1);
    A.setIntegrationElements(ib);
    A.options().setReal("quA", quA);
    A.options().setInt ("quB", quB);
    A.options().setSwitch("SameElement", sameElement);

    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    auto ff = A.getCoeff(f, G); // physical source composed with the map

    gsBoundaryConditions<T> bc;
    for (index_t s = 1; s <= 4; s++)
        bc.addCondition(0, boundary::side(s), condition_type::dirichlet,
                        const_cast<gsFunctionExpr<T>*>(&ms), 0, false);
    bc.setGeoMap(mp);

    gsMatrix<T> solVector;
    auto u_sol = A.getSolution(u, solVector);

    gsStopwatch timer;
    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G),
               u * ff * meas(G));
    const T tAssemble = timer.stop();

    timer.restart();
    typename gsSparseSolver<T>::SimplicialLDLT solver; // direct: removes the
    solver.compute(A.matrix());                        // iterative-solver
    solVector = solver.solve(A.rhs());                 // variable on graded maps
    const T tSolve = timer.stop();

    PoissonRun res;
    res.dofs = A.numDofs();
    res.tAssemble = tAssemble;
    res.tSolve = tSolve;
    u_sol.extract(res.sol);
    const ErrorNorms<T> e =
        evalErrorNorms(mp, ib, res.sol, ms, sameElement, errQuA, errQuB);
    res.l2err = e.l2;
    res.h1err = e.h1;
    return res;
}

// Unconstrained L2 best approximation of the physical function ms in the
// discrete space V_h(M) = span{ N_i o M^{-1} } carried by the map mp:
//
//      M c = b,   M_ij = int N_i N_j |det J_M|,   b_i = int N_i (ms o M) |det J_M|
//
// There are NO boundary conditions (u.setup(-1) leaves every dof free), no
// stiffness matrix and no Poisson solve, so the Dirichlet treatment, the
// consistency error and the linear solver are all absent by construction.
// What remains is a pure property of the SPACE, which is exactly what
// separates "the projected space approximates better" from "the composed
// Galerkin machinery loses accuracy". Returns both the L2 and H1 best-
// approximation errors (see ErrorNorms / evalErrorNorms above).
template <class T>
ErrorNorms<T> bestApproxError(const gsMultiPatch<T> & mp,
                  const gsMultiBasis<T> & mb,
                  const gsMultiBasis<T> & ib,
                  const gsFunctionExpr<T> & ms,
                  bool sameElement,
                  T quA, index_t quB,
                  T errQuA, index_t errQuB,
                  index_t & dofs)
{
    gsExprAssembler<T> A(1, 1);
    A.setIntegrationElements(ib);
    A.options().setReal("quA", quA);
    A.options().setInt ("quB", quB);
    A.options().setSwitch("SameElement", sameElement);

    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    auto uex = A.getCoeff(ms, G); // physical exact solution composed with the map

    gsMatrix<T> solVector;
    auto u_sol = A.getSolution(u, solVector);

    u.setup(-1); // no boundary conditions: the full space
    A.initSystem();
    A.assemble(u * u.tr() * meas(G), u * uex * meas(G));

    typename gsSparseSolver<T>::SimplicialLDLT solver;
    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());

    dofs = A.numDofs();
    gsMultiPatch<T> sol;
    u_sol.extract(sol);
    return evalErrorNorms(mp, ib, sol, ms, sameElement, errQuA, errQuB);
}

// Minimum det(J) of a plain geometry, sampled on a 7x7 grid per element of
// its own basis (same policy as the rh-adaptive fitting driver).
template <class T>
T minDetJ(const gsGeometry<T> & geo)
{
    const gsBasis<T> & b = geo.basis();
    gsVector<unsigned> np(2);
    np << 7, 7;
    gsMatrix<T> pts, der;
    T mn = std::numeric_limits<T>::max();
    for (auto & elem : b.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        geo.deriv_into(pts, der);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, der(0, q) * der(3, q) - der(1, q) * der(2, q));
    }
    return mn;
}

// Minimum det(J) of ANY map (in particular a gsComposedGeometry, whose own
// basis reports the analysis elements and therefore misses sigma's knots)
// sampled on a 7x7 grid per element of the basis \a elems -- pass the union
// basis to resolve sigma's knot mesh.
template <class T>
T minDetJOn(const gsFunction<T> & geo, const gsBasis<T> & elems)
{
    gsVector<unsigned> np(2);
    np << 7, 7;
    gsMatrix<T> pts, der;
    T mn = std::numeric_limits<T>::max();
    for (auto & elem : elems.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        geo.deriv_into(pts, der);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, der(0, q) * der(3, q) - der(1, q) * der(2, q));
    }
    return mn;
}

// max |G - PiG| restricted to the four edges of the parameter domain: the
// projection moves the PHYSICAL boundary, which the interior sup-norm
// distance does not expose. Matters for curved templates.
template <class T>
T boundaryMapDist(const gsFunction<T> & g1, const gsFunction<T> & g2, index_t n = 401)
{
    gsMatrix<T> t = gsMatrix<T>::Zero(1, n);
    for (index_t i = 0; i != n; i++) t(0, i) = static_cast<T>(i) / (n - 1);
    gsMatrix<T> pts(2, 4 * n);
    pts.block(0, 0 * n, 1, n) = t; pts.block(1, 0 * n, 1, n).setZero();  // south
    pts.block(0, 1 * n, 1, n) = t; pts.block(1, 1 * n, 1, n).setOnes();  // north
    pts.block(0, 2 * n, 1, n).setZero(); pts.block(1, 2 * n, 1, n) = t;  // west
    pts.block(0, 3 * n, 1, n).setOnes(); pts.block(1, 3 * n, 1, n) = t;  // east
    gsMatrix<T> v1, v2;
    g1.eval_into(pts, v1);
    g2.eval_into(pts, v2);
    return (v1 - v2).cwiseAbs().maxCoeff();
}

// sigma's knot vector: uniform interior knots k/(R+1), optionally displaced
// by \a shift measured in LEVEL-0 ANALYSIS ELEMENTS (h0). shift = 0 with
// R+1 a power of two nests sigma's knots in a dyadic analysis mesh; any
// other value breaks the nesting by a controlled amount. NOTE the sweep must
// be run at a FIXED analysis level -- one sigma reused across levels changes
// its alignment with every refinement.
gsKnotVector<real_t> sigmaKnotVector(index_t nInterior, index_t deg,
                                     real_t shift, real_t h0)
{
    if (0 == shift || 0 == nInterior)
        return gsKnotVector<real_t>(0, 1, nInterior, deg + 1);

    std::vector<real_t> kn;
    kn.reserve(2 * (deg + 1) + nInterior);
    for (index_t i = 0; i <= deg; i++) kn.push_back(0.0);
    for (index_t k = 1; k <= nInterior; k++)
    {
        const real_t v = static_cast<real_t>(k) / (nInterior + 1) + shift * h0;
        GISMO_ENSURE(v > 0 && v < 1, "shifted sigma knot " << v << " left (0,1)");
        kn.push_back(v);
    }
    for (index_t i = 0; i <= deg; i++) kn.push_back(1.0);
    return gsKnotVector<real_t>(kn, deg);
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    index_t numElevate = 1;   // analysis degree = 1 + numElevate
    index_t numRefine  = 4;   // analysis refinements of the coarsest level
    index_t numLoops   = 1;   // additional uniform refinement levels
    index_t sigmaDeg   = 2;
    index_t sigmaKnots = 5;   // -R 7 makes sigma's knots nested in dyadic meshes
    index_t projRefine = 0;   // extra refinements of the projection basis
    real_t  theta      = 1.0; // monitor smoothing
    real_t  steepness  = 1.0; // monitor sharpening exponent (1: unchanged)
    real_t  penalty    = 1e-2;// fold penalty of the relocation energy
    index_t optIter    = 200;
    index_t optVerbose = 0;
    real_t  optTol     = 1e-6;
    index_t testCase   = 0;
    bool    noSlide    = false;
    bool    plot       = false;
    real_t  cQuA       = 2.0;  // quadrature multiplier of the composed path
    index_t cQuB       = 2;    // extra quadrature nodes of the composed path
    bool    cSameElem  = false;// SameElement flag of the composed path
    bool    bestApprox = false;// also report the L2 BEST APPROXIMATION per path
    real_t  baQuA      = 4.0;  // quadrature of the best-approximation solve
    index_t baQuB      = 2;
    real_t  errQuA     = 4.0;  // quadrature of every L2 error integral
    index_t errQuB     = 2;
    real_t  sigmaShift = 0.0;  // sigma knot displacement, in level-0 elements
    index_t sigmaMode  = 0;    // 0: relocated spline sigma (default)
                               // 1: ANALYTIC sigma (gsFunctionExpr, no knots)
                               // 2: spline interpolant of that analytic sigma
    real_t  sigmaAmp   = 0.5;  // grading amplitude of the analytic sigma
    index_t sigmaFreq  = 2;    // grading frequency n: cells contract where
                               // sin(2 pi n t) has its descending zero
    bool    annulus    = false;// curved PHYSICAL domain: quarter annulus 1<=r<=2
    real_t  frontR0    = -1;   // front radius (<0: 0.25 square / 1.5 annulus)
    real_t  frontW     = 0.05; // front width
    real_t  geomWave   = 0.0;  // curved template: interior control points
                               // displaced by A*(sin(pi x)sin(2pi y),
                               // sin(2pi x)sin(pi y)) -- boundary and
                               // physical domain stay the unit square, but
                               // G = geom o sigma is then degree p_G*p_sigma
                               // on curved pullback pieces (never exactly
                               // representable; tests the nested corollary's
                               // limits)
    std::string output = "proj_r_adaptivity_output";

    gsCmdLine cmd("Quasi-optimal r-adaptivity by projection: relocate sigma "
                  "against a solution monitor, project G = geom o sigma onto "
                  "a plain basis, and compare Poisson solves (uniform / "
                  "composed / projected) in accuracy and assembly cost.");
    cmd.addInt ("e", "elevAnalysis", "Degree elevations of the analysis basis", numElevate);
    cmd.addInt ("r", "refAnalysis", "Uniform refinements of the analysis basis", numRefine);
    cmd.addInt ("l", "loops", "Number of refinement levels to test", numLoops);
    cmd.addInt ("E", "sigmaDeg", "Degree of the sigma map", sigmaDeg);
    cmd.addInt ("R", "sigmaKnots", "Interior knots of the sigma basis (each direction)", sigmaKnots);
    cmd.addInt ("k", "projRefine", "Extra uniform refinements of the projection basis", projRefine);
    cmd.addReal("T", "theta", "Monitor smoothing parameter theta", theta);
    cmd.addReal("", "steepness", "Sharpening exponent p of the monitor: "
                "mon = max(1 - tanh(...)^2, 0)^p. p > 1 narrows the band "
                "around the front; it does NOT change the monitor's range, "
                "so the reachable min det J_sigma is still set by -T/--theta "
                "(floor 1/sqrt(1+theta))", steepness);
    cmd.addReal("P", "penalty", "Fold penalty of the relocation energy", penalty);
    cmd.addInt ("", "optIter", "Optimizer iterations of the relocation", optIter);
    cmd.addInt ("", "optVerbose", "Optimizer verbosity", optVerbose);
    cmd.addReal("", "optTol", "Optimizer minimal gradient length", optTol);
    cmd.addInt ("t", "testCase", "0: circular tanh front at the center; 1: at the bottom-left corner", testCase);
    cmd.addSwitch("noslide", "Fix the boundary control points of sigma", noSlide);
    cmd.addReal("", "cQuA", "Quadrature multiplier of the composed path", cQuA);
    cmd.addInt ("", "cQuB", "Extra quadrature nodes of the composed path", cQuB);
    cmd.addSwitch("cSameElem", "Assemble the composed path with SameElement=true "
                  "(valid here: the space is plain)", cSameElem);
    cmd.addReal("", "geomWave", "Amplitude of a curved (wavy) template geometry "
                "(0: identity)", geomWave);
    cmd.addSwitch("bestApprox", "Also report the unconstrained L2 best approximation "
                  "per path (rows UB/CB/PB): a pure property of the discrete space, "
                  "with no BCs, no stiffness matrix and no Poisson solve", bestApprox);
    cmd.addReal("", "baQuA", "Quadrature multiplier of the best-approximation solve", baQuA);
    cmd.addInt ("", "baQuB", "Extra quadrature nodes of the best-approximation solve", baQuB);
    cmd.addReal("", "errQuA", "Quadrature multiplier of every L2 error integral", errQuA);
    cmd.addInt ("", "errQuB", "Extra quadrature nodes of every L2 error integral", errQuB);
    cmd.addReal("", "sigmaShift", "Displace sigma's interior knots by this many "
                "LEVEL-0 analysis elements (0 = uniform; sweep at a fixed level)", sigmaShift);
    cmd.addInt ("", "sigmaMode", "0: relocated spline sigma; 1: analytic sigma "
                "(no knots at all); 2: spline interpolant of that analytic sigma", sigmaMode);
    cmd.addReal("", "sigmaAmp", "Grading amplitude a of the analytic sigma "
                "s(t) = t + a*sin(2 pi n t)/(2 pi n), |a| < 1", sigmaAmp);
    cmd.addInt ("", "sigmaFreq", "Grading frequency n of the analytic sigma: "
                "cells contract at t = 1/4, 3/4 for n = 2 (square front) and "
                "at t = 1/2 for n = 1 (annulus front)", sigmaFreq);
    cmd.addSwitch("annulus", "Curved PHYSICAL domain: quarter annulus 1<=r<=2 "
                  "(instead of the unit square); PiG is then never exact", annulus);
    cmd.addReal("", "frontR0", "Radius of the tanh front (<0: 0.25 square / 1.5 annulus)", frontR0);
    cmd.addReal("", "frontW", "Width of the tanh front", frontW);
    cmd.addSwitch("plot", "Create ParaView visualization files", plot);
    cmd.addString("o", "output", "Output directory", output);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    GISMO_ENSURE(steepness > 0, "--steepness must be positive");

    if (!output.empty() && output.back() != gsFileManager::getNativePathSeparator())
        output += gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);

    // ------------------------------------------------------------------
    // Manufactured solution (tanh front), source and monitor
    // ------------------------------------------------------------------
    // ms = tanh((0.25 - r)/0.05) + 1,  f = -Δms  (same fronts as
    // monitor_poisson_composed_r-adaptivity_single.cpp); the monitor is the
    // normalized front indicator |∇ms|·0.05 = 1 - tanh(s)^2 in [0,1].
    // ms = tanh((r0 - r)/w) + 1 with r the distance to the front's centre.
    // With T = tanh((r0-r)/w) the source follows analytically from
    // Delta = d_rr + (1/r) d_r:
    //      f = -Delta ms = (1 - T^2) * ( 2T/w^2 + 1/(w r) ).
    // Generating it (instead of hard-coding one expansion) is what lets the
    // front and the domain be varied; the legacy string is kept below purely
    // as a self-check.
    if (frontR0 < 0) frontR0 = annulus ? 1.5 : 0.25; // per-domain default
    auto num = [](real_t v)
    { std::ostringstream o; o << std::setprecision(17) << v; return o.str(); };

    const std::string ctr = (testCase == 1 && !annulus) ? "" : "-0.5";
    const std::string R = annulus ? "sqrt(x^2 + y^2)"
                                  : "sqrt((x" + ctr + ")^2 + (y" + ctr + ")^2)";
    const std::string T = "tanh((" + num(frontR0) + "-" + R + ")/" + num(frontW) + ")";

    const std::string msstring  = T + "+1";
    const std::string fstring   = "(1-(" + T + ")^2)*(2*(" + T + ")/" + num(frontW * frontW)
                                + " + 1/(" + num(frontW) + "*" + R + "))";
    // Sharpening exponent p (--steepness): mon = max(1 - T^2, 0)^p. The
    // max() guards the tiny negative round-off of 1 - tanh^2 near the
    // plateaus, which a fractional power would turn into NaN. At p == 1 the
    // ORIGINAL string is emitted verbatim, so the default run is unchanged.
    const std::string monstring = (1.0 == steepness)
        ? "1 - (" + T + ")^2"
        : "max(1-(" + T + ")^2,0)^" + num(steepness);

    gsFunctionExpr<> f (fstring , 2);
    gsFunctionExpr<> ms(msstring, 2);
    gsFunctionExpr<> monitor(monstring, 2);

    if (!annulus && testCase == 0 && frontR0 == 0.25 && frontW == 0.05)
    {   // self-check: the generated source must equal the hand-expanded string
        // this driver used before, on which every stored result was produced
        gsFunctionExpr<> fRef("((tanh(20*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 5)^2 - 1)*(20*(x-0.5)^2 + 20*(y-0.5)^2)*(40*tanh(20*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 5)*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 1))/((x-0.5)^2 + (y-0.5)^2)^(3/2)", 2);
        gsVector<real_t> a(2), b(2);
        a.setConstant(0.01); b.setConstant(0.99);
        gsVector<unsigned> np(2); np << 51, 51;
        gsMatrix<> pts = gsPointGrid(a, b, np), v1, v2;
        f.eval_into(pts, v1);
        fRef.eval_into(pts, v2);
        const real_t rel = (v1 - v2).cwiseAbs().maxCoeff() / v2.cwiseAbs().maxCoeff();
        gsInfo << "Source self-check       : max rel |f_generated - f_legacy| = "
               << rel << "\n";
        GISMO_ENSURE(rel < 1e-12, "generated source disagrees with the legacy string");
    }

    // ------------------------------------------------------------------
    // Sigma map and relocation (once: sigma does not depend on the
    // analysis basis)
    // ------------------------------------------------------------------
    const real_t h0 = 1.0 / static_cast<real_t>(1 << numRefine); // level-0 element size
    // In analytic mode sigma carries NO knots; the (knot-free) basis is still
    // built so that makeIntegrationBasis raises the integration degree to
    // p*p_sigma exactly as in the spline case -- only the knots differ.
    gsKnotVector<> ks = sigmaKnotVector(1 == sigmaMode ? 0 : sigmaKnots,
                                        sigmaDeg, sigmaShift, h0);
    gsTensorBSplineBasis<2> sbasis(ks, ks);

    // Analytic grading s(t) = t + a*sin(2 pi n t)/(2 pi n): a bijection of
    // [0,1] onto itself (s(0)=0, s(1)=1, s' = 1 + a cos(2 pi n t) > 0 for
    // |a|<1) contracting cells at t = 1/4, 3/4 (n=2, where the square's
    // circular front crosses the axes) or at t = 1/2 (n=1, the annulus
    // front). Being ANALYTIC it has no knot lines at all, so u_ex o G is
    // smooth inside every analysis element no matter how the meshes align.
    const std::string sa = util::to_string(sigmaAmp);
    const std::string sf = util::to_string(2 * sigmaFreq); // 2n, as in 2n*pi
    gsFunctionExpr<real_t> sigmaExpr(
        "x + " + sa + "*sin(" + sf + "*pi*x)/(" + sf + "*pi)",
        "y + " + sa + "*sin(" + sf + "*pi*y)/(" + sf + "*pi)", 2);

    gsSquareDomain<real_t> domain(sbasis);
    if (2 == sigmaMode)
    {   // spline interpolant of the SAME analytic map: identical grading,
        // different representation smoothness -- the controlled pair that
        // isolates smoothness from grading strength
        gsMatrix<> scoefs;
        gsQuasiInterpolate<real_t>::localIntpl(sbasis, sigmaExpr, scoefs);
        gsGeometry<>::uPtr sgeo = sbasis.makeGeometry(give(scoefs));
        domain = gsSquareDomain<real_t>(*sgeo);
    }
    domain.options().addSwitch("Slide", "Boundary controls slide along the boundary", !noSlide);
    domain.applyOptions();

    // sigma actually used as the composition
    const gsFunction<real_t> & sigmaFun = (1 == sigmaMode)
        ? static_cast<const gsFunction<real_t> &>(sigmaExpr)
        : static_cast<const gsFunction<real_t> &>(domain);

    // Coarsest analysis basis; identity geometry on it
    gsKnotVector<> kv({0, 0, 1, 1}, 1);
    gsTensorBSplineBasis<2> tbasis0(kv, kv);
    tbasis0.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbasis0.uniformRefine();

    gsInfo << "Analysis basis (level 0): " << tbasis0 << "\n";
    gsInfo << "Sigma basis             : deg " << sigmaDeg << ", "
           << domain.nControls() << " free controls\n";
    gsInfo << "Monitor shape           : theta = " << theta
           << " (monitor floor 1/sqrt(1+theta) = "
           << 1.0 / math::sqrt(1.0 + theta)
           << ", which bounds how far min det J_sigma can drop), steepness = "
           << steepness << " (band width only)\n";

    // Template geometry: identity, or a wavy reparameterization of the unit
    // square (displacement vanishes on the boundary, so the physical domain
    // and the manufactured solution are unchanged for any amplitude).
    gsTensorBSpline<2> geo0 = []( const gsTensorBSplineBasis<2> & b, real_t A)
    {
        gsMatrix<> c = b.anchors().transpose();
        for (index_t i = 0; i != c.rows(); i++)
        {
            const real_t x = c(i, 0), y = c(i, 1);
            c(i, 0) += A * math::sin(EIGEN_PI * x) * math::sin(2 * EIGEN_PI * y);
            c(i, 1) += A * math::sin(2 * EIGEN_PI * x) * math::sin(EIGEN_PI * y);
        }
        return gsTensorBSpline<2>(b, give(c));
    }(tbasis0, geomWave);
    if (geomWave != 0)
        gsInfo << "Template geometry       : wavy, A = " << geomWave
               << ", min det J = " << minDetJ<real_t>(geo0) << "\n";

    if (annulus)
    {   // Quarter annulus 1 <= r <= 2, (u,v) -> ((1+u) cos(pi v/2), (1+u) sin(pi v/2)):
        // a genuinely CURVED physical domain, so G = geom o sigma has degree
        // p_G*p_sigma on curved pullback pieces and PiG is never exact -- the
        // hard case for the projection. Orientation chosen so det J > 0.
        // The geometry is built ONCE in the level-0 analysis basis and then
        // refined by knot insertion (exact), so the physical domain is
        // identical on every level.
        gsFunctionExpr<real_t> amap("(1+x)*cos(pi*y/2)", "(1+x)*sin(pi*y/2)", 2);
        gsMatrix<> acoefs;
        gsQuasiInterpolate<real_t>::localIntpl(tbasis0, amap, acoefs);
        geo0 = gsTensorBSpline<2>(tbasis0, give(acoefs));

        // how far the spline template is from the EXACT annulus (a fixed
        // constant, independent of the analysis level -- must be reported
        // because it floors every error below)
        gsVector<real_t> a(2), b(2);
        a.setZero(); b.setOnes();
        gsVector<unsigned> np(2); np << 201, 201;
        gsMatrix<> pts = gsPointGrid(a, b, np), v1, v2;
        geo0.eval_into(pts, v1);
        amap.eval_into(pts, v2);
        gsInfo << "Template geometry       : quarter annulus (curved), "
               << "max|geo - exact| = " << (v1 - v2).cwiseAbs().maxCoeff()
               << ", min det J = " << minDetJ<real_t>(geo0) << "\n";
    }

    gsStopwatch timer;
    if (0 != sigmaMode)
    {   // prescribed sigma: no relocation, no optimizer, no monitor -- the
        // grading is fixed by hand so that every dial below varies ONE thing
        gsInfo << "Sigma                   : "
               << (1 == sigmaMode ? "ANALYTIC (no knots), a = "
                                  : "spline interpolant of the analytic map, a = ")
               << sigmaAmp << ", min det J_sigma = "
               << minDetJOn<real_t>(sigmaFun, tbasis0) << "\n";
    }
    else
    {
        gsHLBFGS<real_t> optimizer;
        optimizer.options().setInt ("MaxIterations", optIter);
        optimizer.options().setInt ("Verbose", optVerbose);
        optimizer.options().setReal("MinGradLen", optTol);

        // Monitor in physical coordinates for a curved template
        // (parametric == physical when the template is the identity)
        gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>
            reloc(domain, geo0, monitor, tbasis0, optimizer,
                  geomWave == 0 && !annulus);
        reloc.options().setReal("Penalty", penalty);
        reloc.options().setReal("Smoothing", theta);
        reloc.solve();
        gsInfo << "Relocation: " << timer.stop() << " s, min det J_sigma = "
               << reloc.computeMinJacobian() << "\n";
    }
    const real_t tReloc = timer.stop();

    // ------------------------------------------------------------------
    // Per-level comparison
    // ------------------------------------------------------------------
    std::ofstream csv(output + "comparison.csv");
    csv << "level,path,dofs,tAssemble,tSolve,l2err,h1err,projErr,minDetJ,"
           "mapDist,jacDist,bdrMapDist,minDetJ_G,p,psigma,sigmaKnots,theta,"
           "steepness\n";
    // config echo, appended to every row so sweep CSVs are self-describing
    const std::string cfg = "," + util::to_string(numElevate + 1)
                          + "," + util::to_string(sigmaDeg)
                          + "," + util::to_string(sigmaKnots)
                          + "," + util::to_string(theta)
                          + "," + util::to_string(steepness);

    gsInfo << "\n"
           << std::setw(6)  << "level" << std::setw(6) << "path"
           << std::setw(8)  << "dofs"
           << std::setw(12) << "assemble[s]" << std::setw(10) << "solve[s]"
           << std::setw(13) << "L2 error" << std::setw(13) << "H1 error"
           << std::setw(13) << "proj error" << std::setw(12) << "min detJ"
           << "\n";

    gsTensorBSpline<2> geo(geo0);
    for (index_t lvl = 0; lvl != numLoops; lvl++, geo.uniformRefine())
    {
        const gsTensorBSplineBasis<2> & tbasis =
            static_cast<const gsTensorBSplineBasis<2> &>(geo.basis());

        auto report = [&](const char * path, const PoissonRun & r,
                          const ProjDiag & d = ProjDiag())
        {
            gsInfo << std::setw(6) << lvl << std::setw(6) << path
                   << std::setw(8) << r.dofs
                   << std::fixed << std::setprecision(3)
                   << std::setw(12) << r.tAssemble << std::setw(10) << r.tSolve
                   << std::scientific << std::setprecision(3)
                   << std::setw(13) << r.l2err << std::setw(13) << r.h1err;
            if (d.projErr >= 0)
                gsInfo << std::setw(13) << d.projErr << std::setw(12) << d.minDetJ;
            gsInfo << std::defaultfloat << "\n";
            auto opt = [](real_t v) { return v >= 0 ? util::to_string(v) : std::string(); };
            csv << lvl << "," << path << "," << r.dofs << "," << r.tAssemble
                << "," << r.tSolve << "," << r.l2err << "," << r.h1err << ","
                << opt(d.projErr) << ","
                << (d.projErr >= 0 ? util::to_string(d.minDetJ) : "") << ","
                << opt(d.mapDist) << "," << opt(d.jacDist) << ","
                << opt(d.bdrMapDist) << "," << opt(d.minDetJG)
                << cfg << "\n";
        };

        // best-approximation row: no timings, no BCs, no solve -- see
        // bestApproxError(). Reported as UB/CB/PB next to U/C/P.
        auto reportBA = [&](const char * path, index_t dofs,
                            const ErrorNorms<real_t> & e)
        {
            gsInfo << std::setw(6) << lvl << std::setw(6) << path
                   << std::setw(8) << dofs
                   << std::setw(12) << "-" << std::setw(10) << "-"
                   << std::scientific << std::setprecision(3)
                   << std::setw(13) << e.l2 << std::setw(13) << e.h1
                   << std::defaultfloat << "\n";
            // 7 commas after h1err: the six columns projErr..minDetJ_G are
            // empty on a best-approximation row, and cfg.substr(1) then starts
            // exactly at the `p` column. (The pre-edit line had six commas and
            // was therefore one field SHORT -- `p` landed under minDetJ_G.
            // Fixed here because this is the very line the task rewrites.)
            csv << lvl << "," << path << "," << dofs << ",,,"
                << e.l2 << "," << e.h1 << ",,,,,,," << cfg.substr(1) << "\n";
        };

        // Union integration basis (analysis knots + sigma's knots, degree
        // p*p_sigma): built once per level, used by C, S and by every
        // best-approximation solve so that all paths share ONE integration
        // mesh and ONE quadrature rule.
        gsTensorBSplineBasis<2> ibasis =
            gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                makeIntegrationBasis<2>(tbasis, sbasis);
        gsMultiBasis<> ib(ibasis);

        // --- U: uniform reference -------------------------------------
        {
            gsMultiPatch<> mp;
            mp.addPatch(geo);
            gsMultiBasis<> mb(mp);
            report("U", solvePoisson(mp, mb, mb, f, ms, true, 2.0, 2, errQuA, errQuB));
            if (bestApprox)
            {
                index_t nd;
                const ErrorNorms<real_t> e = bestApproxError(mp, mb, ib, ms, true,
                                                 baQuA, baQuB, errQuA, errQuB, nd);
                reportBA("UB", nd, e);
            }
        }

        // --- projection of the composed map (needed by C's diagnostic
        //     and by the P path) ------------------------------------------
        gsComposedGeometry<real_t> cgeo(sigmaFun, geo);
        gsTensorBSplineBasis<2> pbasis(tbasis);
        for (index_t i = 0; i < projRefine; i++)
            pbasis.uniformRefine();

        // One-time projection, integrated accurately on the union mesh
        gsTensorBSplineBasis<2> pibasis =
            gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                makeIntegrationBasis<2>(pbasis, sbasis);
        timer.restart();
        gsMatrix<> pcoefs;
        const real_t projErr =
            gsL2Projection<real_t>::project(pbasis, pibasis, cgeo, pcoefs);
        const real_t tProj = timer.stop();

        pcoefs.resize(pbasis.size(), 2);
        gsGeometry<>::uPtr pgeo = pbasis.makeGeometry(give(pcoefs));
        ProjDiag diag;
        diag.projErr = projErr;
        {   // diagnostic: sup-norm distance between G and PiG (values and
            // Jacobians) on a fine grid -- the projection error alone does
            // not bound the Jacobian (hence space) perturbation
            gsVector<real_t> a(2), b(2);
            a.setZero(); b.setOnes();
            gsVector<unsigned> np(2);
            np << 201, 201;
            gsMatrix<> pts = gsPointGrid(a, b, np), vC, vP, dC, dP;
            cgeo.eval_into(pts, vC);
            pgeo->eval_into(pts, vP);
            cgeo.deriv_into(pts, dC);
            pgeo->deriv_into(pts, dP);
            diag.mapDist = (vC - vP).cwiseAbs().maxCoeff();
            diag.jacDist = (dC - dP).cwiseAbs().maxCoeff();
            diag.bdrMapDist = boundaryMapDist<real_t>(cgeo, *pgeo);
            diag.minDetJG = minDetJOn<real_t>(cgeo, ibasis);
            gsInfo << "         (map dist: max|G-PiG| = "
                   << std::scientific << std::setprecision(3)
                   << diag.mapDist
                   << ", max|J_G-J_PiG| = " << diag.jacDist
                   << ", bdr = " << diag.bdrMapDist
                   << ", min detJ_G = " << diag.minDetJG
                   << std::defaultfloat << ")\n";
        }
        const real_t mdj = minDetJ(*pgeo);
        diag.minDetJ = mdj;
        if (mdj <= 0)
            gsWarn << "Projected geometry is folded (min det J = " << mdj
                   << ") -- its solve is unreliable; refine the "
                      "projection basis (-k) or reduce theta.\n";
        gsMultiPatch<> pmp;
        pmp.addPatch(*pgeo);
        gsMultiBasis<> pmb(pmp);

        // --- C: composed map, plain space, union quadrature ------------
        {
            gsMultiPatch<> cmp;
            cmp.addPatch(cgeo);
            gsMultiBasis<> mb(tbasis); // plain space: B∘G^{-1} is r-adapted
            PoissonRun r = solvePoisson(cmp, mb, ib, f, ms, cSameElem, cQuA, cQuB,
                                        errQuA, errQuB);
            report("C", r);
            if (bestApprox)
            {
                index_t nd;
                const ErrorNorms<real_t> e = bestApproxError(cmp, mb, ib, ms, cSameElem,
                                                 baQuA, baQuB, errQuA, errQuB, nd);
                reportBA("CB", nd, e);
            }
            // Diagnostic: the same solution spline evaluated through the
            // projected map (valid: ||G - PiG|| ~ projErr). A discrepancy
            // against r.l2err isolates composed-map evaluation/quadrature
            // artifacts in the error integral itself.
            if (projRefine == 0)
            {
                const ErrorNorms<real_t> ce =
                    evalErrorNorms(pmp, pmb, r.sol, ms, true, errQuA, errQuB);
                gsInfo << "         (C error via projected map: L2 "
                       << std::scientific << std::setprecision(3) << ce.l2
                       << ", H1 " << ce.h1
                       << std::defaultfloat << ")\n";
            }
        }

        // --- S: diagnostic — identity template means G = sigma EXACTLY,
        //     so assembling with sigma as a plain spline map must agree
        //     with C up to quadrature. A discrepancy isolates a defect in
        //     the gsComposedGeometry evaluation chain. (identity only)
        if (geomWave == 0 && !annulus && 1 != sigmaMode) // G = sigma only then
        {
            gsMultiPatch<> smp;
            smp.addPatch(domain.domain().clone());
            gsMultiBasis<> mb(tbasis);
            report("S", solvePoisson(smp, mb, ib, f, ms, true, 2.0, 2, errQuA, errQuB));
            // S2: same, but integrating on the ANALYSIS mesh (ignoring
            // sigma's knots). Isolates the union-integration-basis machinery.
            report("S2", solvePoisson(smp, mb, mb, f, ms, true, 2.0, 2, errQuA, errQuB));
        }

        // --- P: projected map, plain space, standard quadrature --------
        {
            PoissonRun r = solvePoisson(pmp, pmb, pmb, f, ms, true, 2.0, 2,
                                        errQuA, errQuB);
            r.tAssemble += tProj; // projection is part of the P pipeline cost
            report("P", r, diag);
            if (bestApprox)
            {   // NB the union mesh ib and the same (baQuA,baQuB) as CB: the
                // two best-approximation numbers differ ONLY in the map.
                index_t nd;
                const ErrorNorms<real_t> e = bestApproxError(pmp, pmb, ib, ms, true,
                                                 baQuA, baQuB, errQuA, errQuB, nd);
                reportBA("PB", nd, e);
            }

            if (plot && lvl == numLoops - 1)
            {
                gsWriteParaview(*pgeo, output + "projected_geometry", 1000, true, true);
                if (1 != sigmaMode)
                    gsWriteParaview(domain.domain(), output + "sigma", 1000, true, true);
                gsWriteParaview(geo0, monitor, output + "monitor", 1000);
            }
        }
    }
    gsInfo << "\nRelocation time (once): " << tReloc << " s\n";
    gsInfo << "Results written to " << output << "comparison.csv\n";

    return 0;
}
