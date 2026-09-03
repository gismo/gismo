/** @file rh-adaptive_poisson_projected_example.cpp

    @brief Schedule-driven ({S,R,U}) restructuring of
    monitor_poisson_projected_r-adaptivity.cpp: quasi-optimal r-adaptivity by
    projection, but with the control flow expressed as a repeatable
    --schedule string instead of a fixed "relocate once, then loop over
    uniform refinement levels" sequence. The NUMERICS are unchanged; see
    that file's header for the full method description. Three Poisson solves
    are compared per S step:
      U  uniform     : identity geometry, plain basis, standard quadrature;
      C  composed    : map G = geom∘σ, plain basis, union integration basis
                       + SameElement=false (the expensive reference);
      P  projected   : map ΠG on a plain basis (analysis basis + `-k` extra
                       refinements), standard quadrature.

    Schedule step letters:
      S  (solve)     the U/C/[S/S2]/P Poisson-solve comparison at the CURRENT
                     state of geo/domain; writes one group of rows to
                     comparison.csv and increments the solve counter `lvl`
                     (the CSV's `level` column) by one.
      R  (relocate)  relocate sigma against the monitor via
                     gsAdaptiveParametrization + HLBFGS. Re-invokable: geo and
                     its basis are read at CALL time (not baked in), so an R
                     step after a U step sees the refined analysis basis.
                     With --sigmaMode != 0 (prescribed sigma) this step does
                     no optimization; it only prints the fixed sigma's min
                     det J on the CURRENT analysis basis.
      U  (refine)    geo.uniformRefine(): a plain uniform h-refinement of the
                     analysis geometry/basis; sigma (the gsSquareDomain) is
                     untouched.

      H  (THB refine) THB h-refinement of the analysis basis/geometry,
                     marking cells from the pointwise PHYSICAL error of the
                     most recent P-path solve (sampled 7x7 per element,
                     see hStep()); controlled by --threshold/--refPercent/
                     --extension. --extension is now a boolean switch on the
                     marker's own element-extension (>0 enables it) rather
                     than a numeric cell count -- see the Notes below.
                     Bootstraps an S step if none has run yet.
                     Only reachable when --schedule contains 'H'; in that
                     case the analysis basis/geometry become a
                     gsTHBSplineBasis<2>/THB geometry (see "THB mode" in the
                     Notes below).

    An empty --schedule (the default) is auto-built from -l/--loops as
      R S (U S)^(loops-1)
    which reproduces the original driver's fixed sequence exactly: relocate
    sigma once, solve, then (uniformly refine, solve) loops-1 more times. Any
    non-empty --schedule is used verbatim (case-insensitively) and repeated
    --iter times, e.g.
      --schedule RSUSUS --iter 1   (equivalent to -l 3, no --schedule)
      --schedule RSU     --iter 3   (relocate/solve/refine, three cycles)
      --schedule RS      --iter 1   (relocate once, solve once, no refinement)

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
    - THB mode (schedule contains 'H'): the analysis basis becomes a
      gsTHBSplineBasis<2> and the analysis geometry its THB counterpart.
      gsAdaptiveParametrization only accepts a plain tensor integration
      basis, so the R step and every union-integration-mesh construction
      substitute thb.tensorLevel(thb.maxLevel()) for the plain tbasis (see
      rh-adaptive_fitting_example.cpp:796-798 for the same workaround). In H
      (THB) mode the projection basis is the FINEST THB level, taken
      uniformly over the whole domain, so the P path's DoF count no longer
      reflects THB adaptivity (it is that of a uniform mesh at the finest
      level). U and C keep true adaptive DoF counts: their discretization
      space is `gsMultiBasis<>(thb)` directly and is unaffected by the
      `pbasis` choice. Reason: `gsL2Projection::project` onto a THB target
      basis is unverified here; keeping the projection plain is deliberate
      scope control.
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
    - The H step marks through the library's gsHElementMarker<2,real_t>
      (src/gsHSplines/gsHElementMarker.h) instead of the ad-hoc point-cloud
      the ad-hoc point-cloud getBoxes/refineThreshold helpers: it is fed one
      MAX-over-samples error per hierarchical element (GARU/relative-
      threshold rule) and marks with "Admissible" in force. The library
      pipeline works on true hierarchical elements and enforces an
      admissible mesh, which the point-cloud helper knows nothing about;
      a file-local refineThreshold (percentile over the sample points) is
      still used to turn --threshold/--refPercent into the scalar threshold
      t, now mapped onto the marker's relative RefineParam.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsHSplines/gsHElementMarker.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsUtils/gsL2Projection.h>
#include <gsUtils/gsQuasiInterpolate.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <vector>

using namespace gismo;

namespace {

// Error value such that refPercent of the points lie above it (adapted from
// gsHFitting<2,T>::setRefineThreshold; formerly examples/gsTHBBoxMarking.h).
real_t refineThreshold(const std::vector<real_t> & errors, real_t refPercent)
{
    GISMO_ENSURE(refPercent > 0 && refPercent <= 1,
                 "refPercent must lie in (0,1], got " << refPercent);
    std::vector<real_t> errorsCopy = errors;
    const size_t i = std::min(cast<real_t, size_t>(errorsCopy.size() * (1.0 - refPercent)),
                              errorsCopy.size() - 1);
    typename std::vector<real_t>::iterator pos = errorsCopy.begin() + i;
    std::nth_element(errorsCopy.begin(), pos, errorsCopy.end());
    return *pos;
}

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

// Per-point Jacobian scaling factor of a 2D-parametric map, from a column of
// deriv_into output (rows grouped by target component: 2*i, 2*i+1 = d/du,
// d/dv of component i). targetDim==2: the SIGNED 2x2 determinant (used
// everywhere in this file for fold detection, sign < 0 <=> folded) -- this
// is the historical formula, unchanged. targetDim>2 (a genuinely curved
// template loaded via --file, e.g. embedded in R^3): there is no
// orientation to "fold" the same way, so this returns the UNSIGNED area
// element sqrt(det(J^T J)) = sqrt(E*G - F^2) (first fundamental form) --
// it detects degeneracy (-> 0) but not sign flips.
template <class T>
T jacobianScale(const gsMatrix<T> & der, index_t col, index_t targetDim)
{
    if (2 == targetDim)
        return der(0, col) * der(3, col) - der(1, col) * der(2, col);
    T E = 0, F = 0, Gm = 0;
    for (index_t i = 0; i != targetDim; i++)
    {
        const T du = der(2 * i,     col);
        const T dv = der(2 * i + 1, col);
        E += du * du; F += du * dv; Gm += dv * dv;
    }
    return math::sqrt(std::max(T(0), E * Gm - F * F));
}

// Minimum Jacobian scale of a plain geometry, sampled on a 7x7 grid per
// element of its own basis (same policy as the rh-adaptive fitting driver).
template <class T>
T minDetJ(const gsGeometry<T> & geo)
{
    const gsBasis<T> & b = geo.basis();
    const index_t td = geo.targetDim();
    gsVector<unsigned> np(2);
    np << 7, 7;
    gsMatrix<T> pts, der;
    T mn = std::numeric_limits<T>::max();
    for (auto & elem : b.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        geo.deriv_into(pts, der);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, jacobianScale<T>(der, q, td));
    }
    return mn;
}

// Minimum Jacobian scale of ANY map (in particular a gsComposedGeometry,
// whose own basis reports the analysis elements and therefore misses
// sigma's knots) sampled on a 7x7 grid per element of the basis \a elems --
// pass the union basis to resolve sigma's knot mesh.
template <class T>
T minDetJOn(const gsFunction<T> & geo, const gsBasis<T> & elems)
{
    const index_t td = geo.targetDim();
    gsVector<unsigned> np(2);
    np << 7, 7;
    gsMatrix<T> pts, der;
    T mn = std::numeric_limits<T>::max();
    for (auto & elem : elems.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        geo.deriv_into(pts, der);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, jacobianScale<T>(der, q, td));
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
    index_t numLoops   = 1;   // additional uniform refinement levels (auto-build only)
    index_t sigmaDeg   = 2;
    index_t sigmaRef   = 3;   // sigma mesh = level-0 analysis basis refined this often
    index_t projRefine = 0;   // extra refinements of the projection basis
    real_t  theta      = 1.0; // monitor smoothing
    real_t  penalty    = 1e-2;// fold penalty of the relocation energy
    index_t optIter    = 200;
    index_t optVerbose = 0;
    real_t  optTol     = 1e-6;
    index_t testCase   = 0;
    std::string geomFile;      // optional --file template geometry (see cmd.addString below)
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
    std::string schedule = "";   // empty: auto-built from -l/--loops
    index_t     maxIter  = 1;    // repetitions of the schedule string
    real_t  threshold  = -1.0;   // H step: marking threshold (<0: use refPercent percentile)
    real_t  refPercent = 0.1;    // H step: fraction of sample points marked for refinement
    index_t extension  = 2;      // H step: >0 enables the marker's element-extension
    bool    noAdmissible = false; // H step: disable the marker's admissible closure (diagnostic)

    gsCmdLine cmd("Schedule-driven quasi-optimal r-adaptivity by projection: "
                  "relocate sigma against a solution monitor, project "
                  "G = geom o sigma onto a plain basis, and compare Poisson "
                  "solves (uniform / composed / projected) in accuracy and "
                  "assembly cost, over a repeatable --schedule of {S,R,U} steps.");
    cmd.addInt ("e", "elevAnalysis", "Degree elevations of the analysis basis", numElevate);
    cmd.addInt ("r", "refAnalysis", "Uniform refinements of the analysis basis", numRefine);
    cmd.addInt ("l", "loops", "Number of refinement levels to test (only used "
                "to auto-build --schedule when it is empty)", numLoops);
    cmd.addString("", "schedule", "Cycle string over {S,R,U} (empty: built from "
                  "-l/--loops as R S (U S)^(l-1))", schedule);
    cmd.addInt ("i", "iter", "Number of schedule cycles (repetitions of --schedule)", maxIter);
    cmd.addInt ("E", "sigmaDeg", "Degree of the sigma map", sigmaDeg);
    cmd.addInt ("R", "sigmaRef", "Uniform refinements of sigma's mesh w.r.t. the "
                "LEVEL-0 analysis basis (2^sigmaRef elements per direction, so sigma "
                "is dyadically nested with every analysis level -- unless --sigmaShift "
                "deliberately breaks that)", sigmaRef);
    cmd.addInt ("k", "projRefine", "Extra uniform refinements of the projection basis", projRefine);
    cmd.addReal("T", "theta", "Monitor smoothing parameter theta", theta);
    cmd.addReal("P", "penalty", "Fold penalty of the relocation energy", penalty);
    cmd.addInt ("", "optIter", "Optimizer iterations of the relocation", optIter);
    cmd.addInt ("", "optVerbose", "Optimizer verbosity", optVerbose);
    cmd.addReal("", "optTol", "Optimizer minimal gradient length", optTol);
    cmd.addInt ("t", "testCase", "0: circular tanh front at the center; 1: at the bottom-left "
                "corner; 2: anisotropic curved ridge u=tanh((y-m(x))/frontW)+1, "
                "m(x)=1/2+0.3 sin(pi x) (from reproduce_rh_anisotropic.cpp); ignored "
                "when --file is given AND the file supplies its own source/solution", testCase);
    cmd.addString("f", "file", "Optional template geometry (single patch, 2D-parametric, any "
                  "target dim), replacing the built-in identity/wavy/annulus template -- "
                  "e.g. filedata/parametrization/monitor_example_surface.xml (geoDim 3) or "
                  "filedata/pde/surfacepoisson_sphere_patch_bvp.xml (geoDim 3, also supplies "
                  "the manufactured problem, see below). If the file also has a "
                  "gsFunctionExpr labelled 'geometry' it is preferred over the first "
                  "MultiPatch found; a MultiPatch at id=0 is tried next. If the file has "
                  "id=1 (source) and id=3 (exact solution) FunctionExpr entries (the "
                  "poisson2_example.cpp / filedata/pde/surfacepoisson_*.xml convention), "
                  "those replace --testCase's built-in planar formulas -- REQUIRED when "
                  "the geometry's target dim != 2, since the flat-Laplacian --testCase "
                  "formulas are not valid on a curved surface. A FunctionExpr labelled "
                  "'function' (monitor_example_surface.xml's own convention) is used as "
                  "the R step's monitor if present.", geomFile);
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
    cmd.addReal("", "threshold",  "H step: marking threshold (-1: use the refPercent percentile)", threshold);
    cmd.addReal("", "refPercent", "H step: fraction of sample points marked for refinement", refPercent);
    cmd.addInt ("", "extension",  "H step: >0 enables the marker's element-extension when "
                "converting marked elements to refinement boxes", extension);
    cmd.addSwitch("", "noAdmissible", "H step: disable the marker's admissible "
                  "closure (diagnostic; the default marking IS admissible)", noAdmissible);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(numLoops >= 1, "--loops must be at least 1");
    GISMO_ENSURE(maxIter >= 1, "--iter must be at least 1");
    GISMO_ENSURE(extension >= 0, "--extension must be non-negative");
    GISMO_ENSURE(testCase != 2 || !annulus,
                 "--testCase 2 (anisotropic ridge) has no annulus variant");
    GISMO_ENSURE(geomFile.empty() || !annulus,
                 "--file replaces the template geometry entirely; --annulus is unused with it");
    GISMO_ENSURE(geomFile.empty() || geomWave == 0,
                 "--file replaces the template geometry entirely; --geomWave is unused with it");

    // Auto-build first, validate afterwards: one code path, no special cases.
    if (schedule.empty())
    {
        // reproduce the base driver's fixed sequence: relocate sigma once,
        // solve, then (uniform refine, solve) numLoops-1 more times
        schedule = "RS";
        for (index_t i = 1; i < numLoops; i++) schedule += "US";
        if (maxIter != 1)
        {
            gsInfo << "Note: an auto-built schedule is run exactly once; "
                      "--iter " << maxIter << " ignored.\n";
            maxIter = 1;
        }
    }
    else if (numLoops != 1)
        gsInfo << "Note: --schedule was given explicitly; -l/--loops ("
               << numLoops << ") is ignored for schedule construction.\n";

    std::string schedUp;   // upper-cased copy: the schedule is dispatched
                           // case-insensitively, so validate the same way
    for (char c : schedule)
    {
        const char C = std::toupper(c);
        GISMO_ENSURE(C == 'S' || C == 'R' || C == 'U' || C == 'H',
                     "Invalid schedule character '" << c << "' (use S, R, U, H)");
        schedUp.push_back(C);
    }
    GISMO_ENSURE(schedUp.find('S') != std::string::npos,
                "A schedule without an 'S' step never solves anything");

    // THB mode switch: an 'H' anywhere in the (normalized) schedule makes the
    // analysis basis/geometry hierarchical for the whole run (see the "THB
    // mode" header Notes bullet for what this substitutes and where).
    const bool useTHB = (schedUp.find('H') != std::string::npos);

    gsInfo << "Schedule                : " << schedule << " x " << maxIter << "\n";
    gsInfo << "THB mode                : " << (useTHB ? "ON (schedule contains 'H')" : "OFF") << "\n";

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

    // haveFileMS: the manufactured problem (source f, exact solution ms) came
    // from --file's id=1/id=3 entries rather than --testCase's built-in
    // planar formulas -- REQUIRED for a geoDim!=2 template (checked once the
    // template is loaded below, since geoDim isn't known yet here).
    // haveMonitor: an R step has something to relocate against: the
    // built-in testCase monstring, or --file's label="function" monitor.
    bool haveFileMS  = false;
    bool haveMonitor = true;
    gsFunctionExpr<> f, ms, monitor;

    if (!geomFile.empty())
    {
        gsFileData<> gfd(geomFile);
        if (gfd.hasId(1) && gfd.hasId(3))
        {
            gfd.getId(1, f);
            gfd.getId(3, ms);
            haveFileMS = true;
            gsInfo << "Manufactured solution   : loaded from --file (id=1 source, "
                      "id=3 exact solution)\n";
            haveMonitor = gfd.hasLabel("function");
            if (haveMonitor)
            {
                gfd.getLabel("function", monitor);
                gsInfo << "Monitor                 : loaded from --file (label=\"function\")\n";
            }
        }
    }

    if (!haveFileMS)
    {
    std::string msstring, fstring, monstring;
    if (2 == testCase)
    {   // Anisotropic curved ridge (reproduce_rh_anisotropic.cpp): a sharp,
        // CURVED front u = tanh((y-m(x))/t)+1, m(x) = 1/2 + 0.3 sin(pi x),
        // unlike the radial fronts above this is not a distance-to-a-point
        // front, so it gets its own T/f derivation. --frontW plays the role
        // of the layer width t; --frontR0 is unused here.
        const std::string T   = "tanh((y-(1/2+0.3*sin(pi*x)))/" + num(frontW) + ")";
        const std::string mp  = "(0.3*pi*cos(pi*x))";   // m'(x)
        const std::string mpp = "(-0.3*pi^2*sin(pi*x))"; // m''(x)
        // f = -Delta u = (2/t^2) sech^2(s) tanh(s) (1+m'^2) + (m''/t) sech^2(s),
        // s=(y-m)/t; sech^2(s) = 1-tanh(s)^2 = 1-T^2, tanh(s) = T.
        msstring  = T + "+1";
        fstring   = "(1-(" + T + ")^2)*((2/" + num(frontW * frontW) + ")*(" + T
                  + ")*(1+(" + mp + ")^2) + (" + mpp + ")/" + num(frontW) + ")";
        monstring = "1 - (" + T + ")^2";
    }
    else
    {
        const std::string ctr = (testCase == 1 && !annulus) ? "" : "-0.5";
        const std::string R = annulus ? "sqrt(x^2 + y^2)"
                                      : "sqrt((x" + ctr + ")^2 + (y" + ctr + ")^2)";
        const std::string T = "tanh((" + num(frontR0) + "-" + R + ")/" + num(frontW) + ")";

        msstring  = T + "+1";
        fstring   = "(1-(" + T + ")^2)*(2*(" + T + ")/" + num(frontW * frontW)
                  + " + 1/(" + num(frontW) + "*" + R + "))";
        monstring = "1 - (" + T + ")^2";
    }

    f = gsFunctionExpr<>(fstring , 2);
    ms = gsFunctionExpr<>(msstring, 2);
    monitor = gsFunctionExpr<>(monstring, 2);   // haveMonitor is already true here

    if (!annulus && geomFile.empty() && testCase == 0 && frontR0 == 0.25 && frontW == 0.05)
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
    } // !haveFileMS

    // ------------------------------------------------------------------
    // Sigma map (relocated by R steps against the monitor; does not depend
    // on the analysis basis except through the R step's tbasis argument)
    // ------------------------------------------------------------------
    const real_t h0 = 1.0 / static_cast<real_t>(1 << numRefine); // level-0 element size
    // In analytic mode sigma carries NO knots; the (knot-free) basis is still
    // built so that makeIntegrationBasis raises the integration degree to
    // p*p_sigma exactly as in the spline case -- only the knots differ.
    // Sigma's mesh is a refinement LEVEL of the level-0 analysis basis, not a
    // free knot count.  sigmaKnotVector unions into the paper's super mesh via
    // makeIntegrationBasis; a non-nested sigma/analysis pair costs roughly their
    // SUM of elements, a nested pair only the finer of the two -- cf. Sec.
    // "Numerical integration scheme".  Every analysis level is dyadic, so tying
    // sigma to the same ladder makes an ACCIDENTALLY non-nested sigma
    // inexpressible.  --sigmaShift still breaks nesting DELIBERATELY, which is
    // the point of that study (see the note on sigmaKnotVector).
    GISMO_ENSURE(sigmaRef >= 0, "sigmaRef must be non-negative");
    const index_t sigmaKnots = (1 << sigmaRef) - 1;
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

    // sigma actually used as the composition. MUST stay a reference: domain
    // mutates in place when an R step runs, and a repeated R step is only
    // visible to later S steps because of this aliasing.
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
           << ", which bounds how far min det J_sigma can drop)\n";

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

    // geoDim: target dimension of the template geometry. 2 for the built-in
    // identity/wavy/annulus templates; may be >2 (a true curved surface,
    // e.g. embedded in R^3) when --file supplies its own geometry.
    index_t geoDim = 2;
    if (!geomFile.empty())
    {   // Replace geo0 wholesale with a quasi-interpolant of the loaded
        // geometry on tbasis0 -- same technique as the annulus template
        // above, generalized to an arbitrary (and possibly non-planar)
        // target dimension.
        gsFileData<> gfd(geomFile);
        gsMultiPatch<real_t> gmp;
        bool haveGeo = false;
        if (gfd.hasLabel("geometry"))
        {   gfd.getLabel("geometry", gmp); haveGeo = true; }
        else if (gfd.hasId(0))
        {   gfd.getId(0, gmp); haveGeo = true; }
        else
            haveGeo = gfd.getFirst(gmp);
        GISMO_ENSURE(haveGeo && gmp.nPatches() >= 1,
                    "--file " << geomFile << ": no MultiPatch geometry found "
                    "(expected label=\"geometry\", id=0, or the first MultiPatch "
                    "in the file)");
        GISMO_ENSURE(gmp.patch(0).parDim() == 2,
                    "--file geometry must be a 2D-parametric patch (planar or a "
                    "surface), got parDim=" << gmp.patch(0).parDim());

        geoDim = gmp.patch(0).targetDim();
        GISMO_ENSURE(haveFileMS || geoDim == 2,
                    "--file geometry has geoDim=" << geoDim << " (a true curved "
                    "surface) but no id=1 source / id=3 exact solution was found "
                    "-- --testCase's flat-Laplacian formulas are not valid there. "
                    "Add id=1/id=3 FunctionExpr entries to the file (see "
                    "filedata/pde/surfacepoisson_sphere_patch_bvp.xml for the "
                    "convention), or use a geoDim=2 file.");

        gsMatrix<> gcoefs;
        gsQuasiInterpolate<real_t>::localIntpl(tbasis0, gmp.patch(0), gcoefs);
        geo0 = gsTensorBSpline<2>(tbasis0, give(gcoefs));
        gsInfo << "Template geometry       : loaded from " << geomFile
               << " (patch 0), geoDim = " << geoDim
               << ", min/max detJ-scale = " << minDetJ<real_t>(geo0) << "\n";
    }

    GISMO_ENSURE(haveMonitor || schedUp.find('R') == std::string::npos || sigmaMode != 0,
                "--schedule has an 'R' step but no monitor is available (--file has "
                "no label=\"function\" entry and --testCase's built-in formulas do "
                "not apply) -- add a monitor to the file, drop 'R' from --schedule, "
                "or use --sigmaMode 1/2 (prescribed sigma, no relocation needed)");

    // THB analysis basis/geometry (built unconditionally; only touched when
    // useTHB). A single-level THB basis built from tbasis0 has the same
    // functions in the same order as tbasis0 itself, so geo0's coefficients
    // carry over unchanged -- this is exactly what Gate 1 (no-op superset)
    // checks. gsTHBSplineBasis::makeGeometry mirrors
    // rh-adaptive_fitting_example.cpp:731.
    gsTHBSplineBasis<2> thb(tbasis0);
    gsGeometry<>::uPtr  hgeo = thb.makeGeometry(geo0.coefs());

    // ------------------------------------------------------------------
    // Schedule state (mutated by the R/S/U/H closures below)
    // ------------------------------------------------------------------
    gsTensorBSpline<2> geo(geo0);   // current analysis geometry, r/h-adapted
                                    // by the R and U steps (non-THB path)
    index_t lvl    = 0;             // solve counter: incremented once per
                                    // completed S step (the CSV's `level`)
    real_t  tReloc = 0;             // accumulated wall-clock time over all R steps

    // Accessors so the S/R/U bodies stay single-copy across both modes: in
    // non-THB mode they route to the plain geo/geo.basis(); in THB mode to
    // hgeo/thb. gsAdaptiveParametrization and the union-integration-basis
    // helper only accept a plain tensor basis, so curTensor() substitutes
    // the finest THB level (see the header Notes bullet on THB mode).
    auto curGeo    = [&]() -> const gsGeometry<real_t> &
    { return useTHB ? *hgeo : static_cast<const gsGeometry<real_t> &>(geo); };
    auto curSpace  = [&]() -> const gsBasis<real_t> &
    { return useTHB ? static_cast<const gsBasis<real_t> &>(thb) : geo.basis(); };
    auto curTensor = [&]() -> const gsTensorBSplineBasis<2> &
    {
        return useTHB ? thb.tensorLevel(thb.maxLevel())
                       : static_cast<const gsTensorBSplineBasis<2> &>(geo.basis());
    };

    // H-step state: the pointwise-error mark is read off the most recent
    // P-path solve (parametric solution + projected geometry), set at the
    // end of every S step.
    gsMultiPatch<real_t>  lastPSol;
    gsGeometry<>::uPtr    lastPGeo;
    bool                  haveSolve = false;

    // ------------------------------------------------------------------
    // CSV/console header
    // ------------------------------------------------------------------
    std::ofstream csv(output + "comparison.csv");
    csv << "level,path,dofs,tAssemble,tSolve,l2err,h1err,projErr,minDetJ,"
           "mapDist,jacDist,bdrMapDist,minDetJ_G,p,psigma,sigmaKnots,theta\n";
    // config echo, appended to every row so sweep CSVs are self-describing
    const std::string cfg = "," + util::to_string(numElevate + 1)
                          + "," + util::to_string(sigmaDeg)
                          + "," + util::to_string(sigmaKnots)
                          + "," + util::to_string(theta);

    gsInfo << "\n"
           << std::setw(6)  << "level" << std::setw(6) << "path"
           << std::setw(8)  << "dofs"
           << std::setw(12) << "assemble[s]" << std::setw(10) << "solve[s]"
           << std::setw(13) << "L2 error" << std::setw(13) << "H1 error"
           << std::setw(13) << "proj error" << std::setw(12) << "min detJ"
           << "\n";

    // Total number of S steps across the whole (repeated) schedule: used
    // only to decide which S step is the LAST one, for --plot.
    const index_t nSolves = static_cast<index_t>(
        std::count(schedUp.begin(), schedUp.end(), 'S')) * maxIter;

    // ------------------------------------------------------------------
    // Step closures
    // ------------------------------------------------------------------

    // R -- relocate sigma against the monitor on the CURRENT analysis basis.
    // Re-invokable: geo/tbasis are read at call time, so a repeated R step
    // after a U step sees the refined geometry. In the auto-built schedule R
    // is the FIRST character, so at that moment geo == geo0 and
    // geo.basis() == tbasis0 -- the arguments are literally the old ones,
    // which is why this reproduces the original driver's single relocation.
    auto rStep = [&]()
    {
        // Site 1 (task 04 spec): the integration basis passed to
        // gsAdaptiveParametrization must stay a plain tensor basis even in
        // THB mode -- curTensor() substitutes the finest THB level.
        const gsTensorBSplineBasis<2> & tbasis = curTensor();
        gsStopwatch timer;
        if (0 != sigmaMode)
        {   // prescribed sigma: no relocation, no optimizer, no monitor -- the
            // grading is fixed by hand so that every dial below varies ONE thing
            gsInfo << "Sigma                   : "
                   << (1 == sigmaMode ? "ANALYTIC (no knots), a = "
                                      : "spline interpolant of the analytic map, a = ")
                   << sigmaAmp << ", min det J_sigma = "
                   << minDetJOn<real_t>(sigmaFun, tbasis) << "\n";
            return;
        }
        gsHLBFGS<real_t> optimizer;
        optimizer.options().setInt ("MaxIterations", optIter);
        optimizer.options().setInt ("Verbose", optVerbose);
        optimizer.options().setReal("MinGradLen", optTol);

        // Monitor in physical coordinates for a curved template
        // (parametric == physical when the template is the identity)
        gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>
            reloc(domain, curGeo(), monitor, tbasis, optimizer,
                  geomWave == 0 && !annulus && geomFile.empty());
        reloc.options().setReal("Penalty", penalty);
        reloc.options().setReal("Smoothing", theta);
        reloc.solve();
        const real_t tR = timer.stop();      // read the watch ONCE: the old
                                             // file called stop() twice,
                                             // harmless there, ill-defined here
        gsInfo << "Relocation: " << tR
               << " s, min det J_sigma = "
               << reloc.computeMinJacobian() << "\n";
        tReloc += tR;
    };

    // U -- plain uniform h-refinement of the analysis geometry/basis; sigma
    // is untouched. In THB mode both thb and hgeo must be refined (defaults
    // only: gsHTensorBasis::uniformRefine(dir) other than -1 is a silent
    // wrong answer under NDEBUG, per rh-adaptive_fitting_example.cpp:973-977).
    auto uStep = [&]()
    {
        if (useTHB) { thb.uniformRefine(); hgeo->uniformRefine(); }
        else geo.uniformRefine();
    };

    // S -- the U/C/[S/S2]/P Poisson-solve comparison at the current geo/sigma
    // state. Rebuilds every per-level quantity fresh (ibasis/ib, cgeo, the
    // projected geometry, ...) so it always picks up the current geo and the
    // current (possibly relocated) sigma.
    auto sStep = [&]()
    {
        // Site 2 (task 04 spec): the union integration mesh below needs a
        // plain tensor basis -- curTensor() substitutes the finest THB level.
        const gsTensorBSplineBasis<2> & tbasis = curTensor();

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

        // --- U: uniform reference ---------------------------------------
        // Discretization space is curSpace() (site (e) of the task 04 spec):
        // gsMultiBasis<>(thb) in THB mode, so U keeps a true adaptive DoF
        // count.
        {
            gsMultiPatch<> mp;
            mp.addPatch(curGeo());
            gsMultiBasis<> mb(curSpace());
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
        gsComposedGeometry<real_t> cgeo(sigmaFun, curGeo());
        // Scoping decision (task 04 spec, section 1(d)): the projection basis
        // stays a PLAIN tensor basis even in THB mode -- built from the
        // finest THB level (curTensor()/tbasis above), never from thb
        // itself. gsL2Projection::project onto a THB target basis is
        // unverified in this codebase; see the header Notes bullet on THB
        // mode for the DoF-count consequence (P's DoF count is then that of
        // a uniform mesh at the finest THB level, not the adaptive one).
        gsTensorBSplineBasis<2> pbasis(tbasis);
        for (index_t i = 0; i < projRefine; i++)
            pbasis.uniformRefine();

        // One-time projection, integrated accurately on the union mesh
        gsTensorBSplineBasis<2> pibasis =
            gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                makeIntegrationBasis<2>(pbasis, sbasis);
        gsStopwatch timer;
        gsMatrix<> pcoefs;
        const real_t projErr =
            gsL2Projection<real_t>::project(pbasis, pibasis, cgeo, pcoefs);
        const real_t tProj = timer.stop();

        pcoefs.resize(pbasis.size(), geoDim);
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
            gsMultiBasis<> mb(curSpace()); // B∘G^{-1} is r-adapted (THB-aware
                                           // discretization space, site (e))
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
        if (geomWave == 0 && !annulus && geomFile.empty() && 1 != sigmaMode) // G = sigma only then
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

            // H-step state: the most recent P solve is what H marks from.
            lastPSol = r.sol;
            lastPGeo = pgeo->clone();
            haveSolve = true;
            if (bestApprox)
            {   // NB the union mesh ib and the same (baQuA,baQuB) as CB: the
                // two best-approximation numbers differ ONLY in the map.
                index_t nd;
                const ErrorNorms<real_t> e = bestApproxError(pmp, pmb, ib, ms, true,
                                                 baQuA, baQuB, errQuA, errQuB, nd);
                reportBA("PB", nd, e);
            }

            if (plot && lvl == nSolves - 1)
            {
                gsWriteParaview(*pgeo, output + "projected_geometry", 1000, true, true);
                if (1 != sigmaMode)
                    gsWriteParaview(domain.domain(), output + "sigma", 1000, true, true);
                gsWriteParaview(geo0, monitor, output + "monitor", 1000);
            }
        }

        lvl++; // AFTER the body: report()/reportBA() must see the pre-increment value
    };

    // H -- THB h-refinement of the analysis basis/geometry, marking cells
    // from the pointwise PHYSICAL error of the last P-path solve. Only
    // reachable when useTHB (an 'H' in the schedule implies it); bootstraps
    // an S step if none has run yet (dispatch switch below), mirroring
    // rh-adaptive_fitting_example.cpp:1029 (`if (!geom) fitStep();`).
    auto hStep = [&]() -> bool
    {
        GISMO_ENSURE(useTHB, "hStep() requires useTHB (schedule must contain 'H')");

        // 1. sample a 7x7-per-element grid of the analysis parameter domain
        //    (the same density/iteration pattern as minDetJ/minDetJOn), AND
        //    record each element's id() so the per-element error below can
        //    be indexed exactly as gsHElementMarker::setErrors requires.
        //    beginAll()/endAll() (not allElements()) is used explicitly:
        //    under _OPENMP, allElements() returns a per-thread chunk, and
        //    id() is only guaranteed to match visit order on the full range.
        gsVector<unsigned> np(2);
        np << 7, 7;
        std::vector<real_t> pts;
        std::vector<index_t> eIds; // one entry per element, in visit order
        {
            auto dom = curSpace().domain(); // keep the gsHDomain alive (see
                                            // the lifetime caution, task 04 spec)
            for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
            {
                const gsMatrix<real_t> grid =
                    gsPointGrid(it.lowerCorner(), it.upperCorner(), np);
                for (index_t c = 0; c != grid.cols(); c++)
                {
                    pts.push_back(grid(0, c));
                    pts.push_back(grid(1, c));
                }
                eIds.push_back(it.id());
            }
        }
        const index_t M = static_cast<index_t>(pts.size() / 2);
        gsMatrix<real_t> xi(2, M);
        for (index_t i = 0; i != M; i++)
        {
            xi(0, i) = pts[2 * i];
            xi(1, i) = pts[2 * i + 1];
        }

        // 2. pointwise PHYSICAL error of the last P solve at those parameters
        gsMatrix<real_t> uvals, xvals, exvals;
        lastPSol.patch(0).eval_into(xi, uvals); // parametric solution at xi
        lastPGeo->eval_into(xi, xvals);         // physical points  pgeo(xi)
        ms.eval_into(xvals, exvals);            // exact solution there
        std::vector<real_t> errors(M);
        for (index_t i = 0; i != M; i++)
            errors[i] = math::abs(uvals(0, i) - exvals(0, i));

        // 3. per-element error aggregation: MAX over each element's 49
        //    samples. This reproduces the old getBoxes semantics ("does any
        //    sample point in this cell exceed the threshold") as closely as
        //    possible while feeding gsHElementMarker one error per
        //    hierarchical element, indexed by the domain iterator's id().
        const index_t nElem = static_cast<index_t>(eIds.size());
        GISMO_ENSURE(M == 49 * nElem, "H step: sample count does not match 7x7 per element");
        GISMO_ENSURE(nElem == (index_t)thb.numElements(),
                     "H step: element count mismatch between the sampling loop and thb");
        std::vector<real_t> elementErrors(nElem, 0.0);
        for (index_t k = 0; k != nElem; k++)
        {
            real_t emax = 0.0;
            for (index_t j = 0; j != 49; j++) emax = math::max(emax, errors[49 * k + j]);
            elementErrors[eIds[k]] = emax; // NOTE: indexed by id(), not by k
        }

        // 4. threshold + marking through the library's admissible pipeline.
        //    t's computation stays byte-identical to before (still driven by
        //    the POINTWISE sample list, so --refPercent keeps its meaning);
        //    GARU has no absolute threshold, only a relative RefineParam, so
        //    t is mapped onto RefineParam = t / max(elementErrors), which
        //    marks exactly elemErr >= t (gsHElementMarker.hpp:205,209).
        const real_t t = (threshold >= 0) ? threshold : refineThreshold(errors, refPercent);
        const index_t dofsBefore = curSpace().size();

        const real_t maxElemErr = *std::max_element(elementErrors.begin(), elementErrors.end());
        if (maxElemErr <= 0.0) // degenerate: RefineParam would be a 0 threshold and mark EVERYTHING
        {
            gsInfo << "    [H] marked 0 boxes | analysis DoFs " << dofsBefore
                   << " -> " << dofsBefore << " | threshold " << t << "\n";
            gsInfo << "    (no cells marked for refinement)\n";
            return false;
        }

        gsHElementMarker<2,real_t> marker(thb);
        marker.options().setInt   ("RefineRule", 1);          // GARU: relative threshold rule
        marker.options().setInt   ("MaxLevel",   10);         // library default (3) would silently
                                                               // cap refinement; 10 matches
                                                               // gsHElement_marking_example.cpp:82
        marker.options().setSwitch("Admissible", !noAdmissible);
        marker.options().setSwitch("Extension",  extension > 0);
        marker.options().setReal  ("RefineParam", t / maxElemErr); // => marks exactly elemErr >= t

        marker.setErrors(elementErrors);
        auto marked = marker.markRef();
        std::vector<index_t> boxes = marker.toRefBoxes(marked);
        if (boxes.empty())
        {
            gsInfo << "    [H] marked 0 boxes | analysis DoFs " << dofsBefore
                   << " -> " << dofsBefore << " | threshold " << t << "\n";
            gsInfo << "    (no cells marked for refinement)\n";
            return false;
        }
        thb.refineElements(boxes);
        hgeo->refineElements(boxes);
        const index_t dofsAfter = curSpace().size();
        gsInfo << "    [H] marked " << marked.size() << " boxes | analysis DoFs "
               << dofsBefore << " -> " << dofsAfter << " | threshold " << t << "\n";
        return true;
    };

    // ------------------------------------------------------------------
    // Schedule loop
    // ------------------------------------------------------------------
    for (index_t it = 0; it != maxIter; ++it)
    {
        gsInfo << "----------------\n";
        gsInfo << "Cycle " << it << " [" << schedule << "]\n";
        for (char c : schedUp)
        {
            switch (c)
            {
            case 'R': rStep(); break;
            case 'S': sStep(); break;
            case 'U': uStep(); break;
            case 'H':
                if (!haveSolve) sStep(); // H marks from the last P solve: bootstrap
                hStep();
                break;
            }
        }
    }

    gsInfo << "\nRelocation time (total over R steps): " << tReloc << " s\n";
    gsInfo << "Results written to " << output << "comparison.csv\n";

    return 0;
}
