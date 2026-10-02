/** @file poisson2_ghost_penalty_example.cpp

    @brief Immersed Poisson on a circle inside [0,1]^2, comparing full-domain
    FCM regularization against exterior-DOF elimination plus ghost-penalty
    face stabilization.

    Forked from poisson2_nitsche_immersed_example.cpp (kept byte-identical),
    this driver isolates the two jobs the fork's single FCM term
    (alpha * igrad(u,G)*igrad(u,G).tr() over the FULL background mesh) does
    at once:

      1. it gives DOFs with ZERO physical support a non-singular row (pure
         regularization of rows trimmed integration never touches), and
      2. it masks the SLIVER conditioning of DOFs whose physical support is
         tiny but non-zero (well-posed, but numerically hostile).

    Exterior-DOF elimination (a gsDofMapper that keeps only the DOFs active
    on an interior-or-cut element and removes the rest via eliminateDof())
    replaces job 1 only. Ghost-penalty face stabilization (a jump penalty on
    the highest-order normal derivative across skeleton faces touching a cut
    element) replaces job 2 only. Neither replaces the other, so this driver
    exposes three arms behind --arm:

      baseline    alpha=1e-10 FCM term, dirichlet::none setup, no
                  elimination, no ghost. Exactly the fork's formulation:
                  both jobs done by one term. kappa stays flat as delta/h
                  shrinks -- the mask is doing its job, at the cost of an
                  O(alpha) perturbation of the physical solution.
      control     alpha=0 (no FCM term at all), elimination on, ghost off.
                  Job 1 is solved, job 2 is not: kappa is EXPECTED to grow
                  by orders of magnitude as delta/h -> 0, and the matrix may
                  become indefinite or the direct solve may fail to produce
                  a finite vector at the smallest classified ratio. That is
                  the measurement, not a bug -- it isolates what ghost adds
                  by showing what elimination alone does not fix. Do not
                  add a fallback regularization to this arm.
      stabilized  alpha=0, elimination on, ghost penalty on. Both jobs
                  solved by the mechanism that actually targets each: kappa
                  should stay roughly flat across the sweep, like baseline,
                  but without alpha's O(alpha) solution perturbation.

    Ghost term and continuity. The background basis is built by
    dbasis.setDegree(degree) BEFORE uniformRefine(), which inserts single
    interior knots, so the space is C^{degree-1}: every derivative order
    j < degree has an identically-zero jump [[d^j u]] in exact arithmetic.
    Only the j = degree term is assembled -- adding lower orders would not
    add stabilization, it would inject the faceShift evaluation-offset floor
    (see below) into the stiffness matrix as signal.

    Face integrals never evaluate exactly on the face: gsExprAssembler and
    gsExprEvaluator shift face-loop points by -+faceShift*h (faceShift =
    1e-6 by default). Any derivative order below the jump order therefore
    reports a deterministic (faceShift*h)^2 evaluation offset, not zero and
    not roundoff -- about 20 orders above the true roundoff floor. The
    driver prints this statement next to the one face-integral diagnostic it
    reports (the ghost jump energy of the computed solution, Eghost).

    Two independent studies, --study sweep|convergence|both:

      sweep        One PINNED refinement level (--sweepRefine), five
                   delta/h ratios poking a sliver of depth delta = ratio*h
                   into one fixed background cell. A dense
                   SelfAdjointEigenSolver (gsSpectra is not enabled in this
                   checkout) gives exact condition numbers -- affordable
                   only because the level is pinned (O(n^3), never run
                   inside a refinement loop). Rows are only ever printed for
                   ratios whose sliver cell is independently verified cut
                   (sign == 0) by the same kd-tree classification the
                   assembler uses; kd-tree classification uses Lobatto
                   sampling and stops resolving the sliver below some
                   delta/h, called the classification cliff below which the
                   sliver silently vanishes and a good-looking kappa would
                   be measuring nothing.
      convergence  A LEVEL-INDEPENDENT radius (--radius, default 0.4, same
                   centre) run through a genuine refinement loop, no
                   condition numbers. A radius grazing the finest level (the
                   sweep's construction) would put delta/h BELOW the
                   classification cliff at every coarser level, putting
                   unverified geometry at the coarse end of a loop whose
                   whole purpose is a convergence rate; a fixed radius keeps
                   cut cells and ghost faces present at every level instead.

    The level-set centre is offset from 0.5 by 1e-3 in both studies: an
    exactly tangent level set (centre on a mesh line) hangs the Algoim
    surface rule.

    Example command lines:
      ./poisson2_ghost_penalty_example --arm baseline
      ./poisson2_ghost_penalty_example --arm control
      ./poisson2_ghost_penalty_example --arm stabilized
      ./poisson2_ghost_penalty_example --arm stabilized --study sweep
      ./poisson2_ghost_penalty_example --arm stabilized --study convergence --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsCore/gsDofMapper.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

namespace {

// Integer grid coordinate / flat level-0 element index of a point on the
// uniform N-element grid of [0,1]^2. Duplicated (not shared) from
// examples/immersed_face_iterators_example.cpp: kept independent on
// purpose, see that file's header for the rationale.
index_t gridCoord(const real_t x, const index_t n)
{ return static_cast<index_t>(math::round(x * (real_t)n)); }

size_t flatIndex(const gsVector<real_t> & p, const index_t n)
{
    size_t flat = 0, stride = 1;
    for (index_t j = 0; j != p.rows(); ++j)
    { flat += (size_t)gridCoord(p[j], n) * stride; stride *= (size_t)n; }
    return flat;
}

// Sign of the level-0 background element at flat index \a targetFlat, read
// from the trimmed domain's own interior/cut classification (sign < 0:
// interior, sign == 0: cut, sign > 0: exterior/unvisited). A single-cell
// lookup, not a full sign grid: the sweep only ever needs the sign of the
// one cell the sliver is engineered into.
short_t sliverSign(const gsImplicitTrimmedDomain<2,real_t> & dom, size_t targetFlat, index_t N)
{
    gsDomain<real_t>::iterator eInt = dom.end<InteriorSign>();
    for (gsDomain<real_t>::iterator it = dom.beginInterior(); it < eInt; ++it)
        if (flatIndex(it.lowerCorner(), N) == targetFlat) return -1;
    gsDomain<real_t>::iterator eCut = dom.endBdr(boundary::none);
    for (gsDomain<real_t>::iterator it = dom.beginBdr(boundary::none); it < eCut; ++it)
        if (flatIndex(it.lowerCorner(), N) == targetFlat) return 0;
    return +1;
}

// Pretty-print helpers, idiom from examples/momfit_operator_probe_example.cpp.
std::string fmtSci(real_t v, int prec = 6)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

std::string fmtCond(real_t lmin, real_t lmax)
{
    if (lmin <= 0) return "(indef)";
    return fmtSci(lmax / lmin, 2);
}

const std::string faceShiftLegend =
    "Face integrals are evaluated at x -/+ faceShift*h (faceShift = 1e-6): quantities of "
    "order (faceShift*h)^2 ~ 1e-15 are the evaluation-offset floor, not convergence to zero.";

} // anonymous namespace

/// Plain-scalar summary of one (arm, geometry, refinement level) solve.
/// Never holds the trimmed domain or anything referencing the level set:
/// gsImplicitTrimmedDomain holds phi non-owned, and phi is a per-call local
/// in solveOne() below.
struct RunResult
{
    index_t ndofs        = 0;
    real_t  l2            = std::numeric_limits<real_t>::quiet_NaN();
    real_t  h1            = std::numeric_limits<real_t>::quiet_NaN();
    real_t  lmin          = 0;
    real_t  lmax          = 0;
    real_t  minRowNorm    = 0;
    real_t  Eghost        = std::numeric_limits<real_t>::quiet_NaN();
    index_t nCut          = 0;
    index_t nGhost        = 0;
    bool    haveSpectrum  = false;
    bool    finite        = true;
};

/// Assembles and solves one configuration. Builds a fresh gsExprAssembler
/// (so A.options().addInt("quDim",...) runs at most once per options
/// object) and a fresh trimmed domain from \a phi. \a ratio is the sweep's
/// delta/h (negative selects "n/a" for the convergence loop, which has no
/// sliver ratio); it is only used to label the zero-row GISMO_ENSURE.
RunResult solveOne(const std::string & armName,
                    const gsMultiPatch<real_t> & mp,
                    gsMultiBasis<real_t> & dbasis,
                    const gsFunctionExpr<real_t> & phi,
                    const gsFunctionExpr<real_t> & nImm,
                    const gsFunctionExpr<real_t> & u_exact,
                    const gsFunctionExpr<real_t> & f_rhs,
                    const gsBoundaryConditions<real_t> & bc,
                    index_t degree,
                    index_t N,
                    real_t ratio,
                    real_t hmax,
                    real_t gammaEff,
                    real_t gtEff,
                    index_t quRule,
                    bool computeSpectrum,
                    bool doPlot,
                    const std::string & outPath)
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    const bool isBaseline   = ("baseline"   == armName);
    const bool isStabilized = ("stabilized" == armName);
    const index_t k = degree;

    RunResult result;

    gsExprAssembler<> A(1, 1);
    A.options().setInt("quRule", quRule);

    // Trimmed (implicit) integration domain -- not installed on A yet.
    gsTensorBSplineBasis<2,real_t> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");
    memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > tr_domain =
        memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phi, *tbsPtr));
    GISMO_ENSURE(1 == tr_domain->numLevels(), "gsTrimmedDomain requires a single kd-tree level.");

    // Warm the (lazily built, non-thread-safe) face sign grid once, single-
    // threaded, before any OpenMP assembly region touches it.
    result.nGhost = static_cast<index_t>(tr_domain->numGhostFaces());
    result.nCut   = static_cast<index_t>(tr_domain->numElementsBdr(boundary::none));

    gsExprEvaluator<> ev(A);
    ev.options().setInt("quRule", quRule);

    std::vector<patchSide> bdr_immersed(1);
    bdr_immersed[0] = patchSide(0, boundary::none);

    geometryMap G = A.getMap(mp);
    space u = A.getSpace(dbasis);
    auto ff    = A.getCoeff(f_rhs, G);
    auto n_imm = A.getCoeff(nImm, G);
    auto u_ex  = ev.getVariable(u_exact, G);
    gsMatrix<real_t> solVector;                    // NAMED: getSolution keeps a pointer to it
    solution u_sol = A.getSolution(u, solVector);

    // Lazily-tagged jump expression: constructing it has no evaluation-time
    // side effect (symbol_expr::jump() only sets a side tag), so it is safe
    // to declare unconditionally and use it only on the stabilized arm.
    auto dJ = dnk(u.jump(), G.left(), k);

    // Arm-dependent dof setup.
    if (isBaseline)
    {
        u.setup(bc, dirichlet::none, 0);
    }
    else
    {
        // Exterior-DOF elimination: keep exactly the DOFs active on an
        // interior-or-cut element; eliminate (fix to zero) every other DOF,
        // whose entire support lies where trimmed integration never
        // touches. Zero is both the correct extension of the discrete
        // solution outside Omega and what gsFeSpace::getCoeffs needs to
        // reconstruct plotted/evaluated fields correctly.
        gsDofMapper mapper(dbasis, 1);
        std::vector<bool> keep(dbasis.basis(0).size(), false);
        gsMatrix<real_t> centre(2, 1);
        gsMatrix<index_t> act;

        gsDomain<real_t>::iterator eInt = tr_domain->end<InteriorSign>();
        for (gsDomain<real_t>::iterator it = tr_domain->beginInterior(); it < eInt; ++it)
        {
            centre.col(0) = 0.5 * (it.lowerCorner() + it.upperCorner());
            dbasis.basis(0).active_into(centre, act);
            for (index_t i = 0; i != act.rows(); ++i) keep[act(i,0)] = true;
        }
        gsDomain<real_t>::iterator eCut = tr_domain->endBdr(boundary::none);
        for (gsDomain<real_t>::iterator it = tr_domain->beginBdr(boundary::none); it < eCut; ++it)
        {
            centre.col(0) = 0.5 * (it.lowerCorner() + it.upperCorner());
            dbasis.basis(0).active_into(centre, act);
            for (index_t i = 0; i != act.rows(); ++i) keep[act(i,0)] = true;
        }

        for (index_t i = 0; i != (index_t)keep.size(); ++i)
            if (!keep[i]) mapper.eliminateDof(i, 0);
        mapper.finalize();
        u.setupMapper(mapper);
        // gsExprAssembler::space is `const expr::gsFeSpace<T>`; the non-const
        // fixedPart() overload is reached through the same const_cast idiom
        // used throughout the library (e.g. gsAssembler/gsDirichletValues.h).
        const_cast<expr::gsFeSpace<real_t>&>(u).fixedPart().setZero(mapper.boundarySize(), 1);
    }

    A.setIntegrationElements(dbasis);
    A.initSystem();

    A.computePattern(igrad(u) * igrad(u).tr());

    // FCM stabilization, baseline arm only, over the FULL background
    // mesh with standard Gauss quadrature (see file header, job 1 and 2).
    if (isBaseline)
    {
        const real_t alpha = 1e-10;
        A.options().setInt("quRule", gsQuadrature::GaussLegendre);
        A.assemble(alpha * igrad(u, G) * igrad(u, G).tr() * meas(G));
    }

    // Restrict to the trimmed physical domain for the remaining terms.
    A.setIntegrationDomain(tr_domain);
    A.options().setInt("quRule", quRule);
    if (isStabilized)
        A.computePatternGhost(dJ * dJ.tr());

    A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G),
        u * ff * meas(G)
    );

    // Nitsche imposition of the immersed Dirichlet datum; surfMeas is the
    // Nanson factor (identically 1 for the identity map).
    auto g_D = A.getCoeff(u_exact, G);
    auto surfMeas = meas(G) * (jac(G).inv().tr() * n_imm).norm();
    A.options().addInt("quDim", "Surface (phi==0) quadrature selector", 2);
    A.assembleBdr(
        bdr_immersed,
        - (igrad(u, G) * n_imm) * u.tr()       * surfMeas
        - u * (igrad(u, G) * n_imm).tr()       * surfMeas
        + gammaEff / hmax * u * u.tr()         * surfMeas
    );
    A.assembleBdr(
        bdr_immersed,
        - (igrad(u, G) * n_imm) * g_D          * surfMeas
        + gammaEff / hmax * u * g_D            * surfMeas
    );
    A.options().setInt("quDim", -1);

    // Ghost penalty, stabilized arm only, LAST assembly call.
    if (isStabilized)
        A.assembleGhost(gtEff * math::pow(hmax, 2 * k - 1) * dJ * dJ.tr());

    // Exterior-DOF elimination (or the alpha mask on
    // baseline) must have left no identically-zero matrix row.
    {
        const gsSparseMatrix<real_t> & M = A.matrix();
        gsVector<real_t> rowSq = gsVector<real_t>::Zero(M.rows());
        for (index_t c = 0; c != M.outerSize(); ++c)
            for (gsSparseMatrix<real_t>::iterator it(M, c); it; ++it)
                rowSq[it.row()] += it.value() * it.value();
        index_t nZeroRows = 0;
        result.minRowNorm = std::numeric_limits<real_t>::max();
        for (index_t i = 0; i != rowSq.rows(); ++i)
        {
            result.minRowNorm = math::min(result.minRowNorm, math::sqrt(rowSq[i]));
            if (0.0 == rowSq[i]) ++nZeroRows;
        }
        const std::string ratioLabel = (ratio < 0) ? std::string("n/a") : fmtSci(ratio, 2);
        GISMO_ENSURE(0 == nZeroRows, "arm=" << armName << " N=" << N << " degree=" << degree
                     << " delta/h=" << ratioLabel << ": " << nZeroRows
                     << " of " << rowSq.rows() << " matrix rows are identically zero"
                     << " (exterior-DOF elimination did not remove every unsupported DOF).");
    }

    result.ndofs = A.numDofs();

    // Condition number, sweep-only: dense copy at the pinned level, plain
    // Eigen (gsSpectra is not enabled in this checkout).
    if (computeSpectrum)
    {
        const gsMatrix<real_t> dense = A.matrix().toDense();
        gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> es(dense);
        GISMO_ENSURE(gsEigen::Success == es.info(), "Eigen solver failed.");
        result.lmin = es.eigenvalues().minCoeff();
        result.lmax = es.eigenvalues().maxCoeff();
        result.haveSpectrum = true;
    }

    // Solve with a direct sparse solver: the control arm's kappa can pass
    // 1e15 by design, where an iterative solver hits its cap and returns
    // garbage rather than a diagnosable non-finite vector.
    gsSparseSolver<real_t>::LU solver;
    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());
    result.finite = solVector.allFinite();

    if (result.finite)
    {
        result.l2 = math::sqrt(ev.integral((u_ex - u_sol).sqNorm() * meas(G)));
        result.h1 = result.l2 + math::sqrt(ev.integral((igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));

        // The one face-integral diagnostic (see file header for the
        // faceShift floor it is subject to), gated on a finite solve and a
        // non-empty ghost set.
        if (result.nGhost > 0)
            result.Eghost = ev.integralGhost(dnk(u_sol.jump(), G.left(), k) * dnk(u_sol.jump(), G.left(), k));
    }

    if (doPlot && result.finite)
    {
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 100000);
        ev.writeParaview(u_sol, G, outPath + "/poisson2_ghost_penalty_solution");
        ev.writeParaview(u_ex, G, outPath + "/poisson2_ghost_penalty_exact");
        gsMesh<> mesh(dbasis.basis(0));
        gsWriteParaview(mesh, outPath + "/poisson2_ghost_penalty_background_mesh");
    }

    return result;
}

/// Sliver conditioning sweep at one pinned refinement level: a
/// circle radius is engineered so its rightmost point pokes a sliver of
/// depth delta = ratio*h into one fixed background cell, for a hard-coded
/// ratio list. Every printed row is gated on that cell being independently
/// classified cut by the trimmed domain's own kd-tree sign, exactly like
/// the assembler sees it; ratios that fail this check are SKIPPED, not
/// silently reported with a meaningless kappa.
void runSweep(const std::string & armName, index_t degree, index_t sweepRefine,
              real_t gammaEff, real_t gtEff, index_t quRule,
              const gsMultiPatch<real_t> & mp, const gsBoundaryConditions<real_t> & bc,
              const gsFunctionExpr<real_t> & u_exact, const gsFunctionExpr<real_t> & f_rhs,
              index_t N0, real_t c0)
{
    gsMultiBasis<real_t> dbasis(mp, true);
    dbasis.setDegree(degree);            // elevate BEFORE refining, or continuity drops to C^0
    dbasis.uniformRefine(N0 - 1);        // level 0: N0 elements per direction
    for (index_t s = 0; s != sweepRefine; ++s) dbasis.uniformRefine();

    const index_t N = N0 << sweepRefine;
    real_t hmax = 0;
    for (size_t p = 0; p != dbasis.nBases(); ++p)
        hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());

    const real_t h0    = 1.0 / (real_t)N0;
    const index_t m    = (index_t)std::floor((c0 + 0.3) / h0);
    const real_t xStar = (real_t)m * h0;

    const index_t sliverIdx0 = (index_t)math::round(xStar / hmax);
    const index_t sliverIdx1 = (index_t)std::floor(c0 / hmax);
    const size_t  sliverFlat = (size_t)sliverIdx0 + (size_t)sliverIdx1 * (size_t)N;

    const real_t ratios[] = { 0.5, 1e-1, 1e-2, 1e-3, 1e-4 };
    const index_t nRatios = 5;

    gsInfo << "\n=== Sliver sweep  arm=" << armName << "  N=" << N
           << "  h=" << fmtSci(hmax, 4) << "  degree=" << degree
           << "  gamma=" << gammaEff << "  ghost=" << fmtSci(gtEff, 0) << " ===\n";
    gsInfo << std::right
           << std::setw(10) << "delta/h" << std::setw(13) << "delta" << std::setw(13) << "radius"
           << std::setw(6) << "#cut" << std::setw(7) << "#ghost" << std::setw(7) << "dofs"
           << std::setw(15) << "lambdaMin" << std::setw(15) << "lambdaMax" << std::setw(11) << "kappa"
           << std::setw(13) << "L2err" << std::setw(13) << "H1err"
           << std::setw(13) << "minRowNorm" << std::setw(11) << "Eghost" << "\n";

    bool cliffFound = false;
    real_t cliffRatio = 0;

    for (index_t j = 0; j != nRatios; ++j)
    {
        const real_t ratio  = ratios[j];
        const real_t delta  = ratio * hmax;
        const real_t radius = xStar - c0 + delta;

        std::ostringstream phiStr; phiStr << std::setprecision(17)
            << "sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)-" << radius;
        std::ostringstream nxStr; nxStr << std::setprecision(17)
            << "(x-" << c0 << ")/sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)";
        std::ostringstream nyStr; nyStr << std::setprecision(17)
            << "(y-" << c0 << ")/sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)";
        gsFunctionExpr<real_t> phi(phiStr.str(), 2);      // outlives every trim built below it
        gsFunctionExpr<real_t> nImm(nxStr.str(), nyStr.str(), 2);

        gsTensorBSplineBasis<2,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");
        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > classifyDomain =
            memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phi, *tbsPtr));
        const short_t sign = sliverSign(*classifyDomain, sliverFlat, N);

        // Non-vacuity: the coarsest sweep point MUST classify cut, or the
        // whole geometric construction is broken and every row is meaningless.
        // Checked BEFORE the skip branch below, or a coarsest-ratio failure
        // would `continue` past it and the ENSURE would never fire.
        GISMO_ENSURE(0 != j || 0 == sign, "Sliver sweep: the coarsest ratio (delta/h=0.5) did "
                     "not classify the sliver cell as cut -- the sweep geometry is broken.");

        if (0 != sign)
        {
            gsInfo << std::setw(10) << fmtSci(ratio, 2)
                   << "  SKIPPED: sliver cell (" << sliverIdx0 << "," << sliverIdx1
                   << ") sign=" << (int)sign << " -- not classified cut\n";
            if (!cliffFound) { cliffFound = true; cliffRatio = ratio; }
            continue;
        }

        RunResult r = solveOne(armName, mp, dbasis, phi, nImm, u_exact, f_rhs, bc,
                                degree, N, ratio, hmax, gammaEff, gtEff, quRule,
                                /*computeSpectrum*/ true, /*doPlot*/ false, "");

        gsInfo << std::setw(10) << fmtSci(ratio, 2)
               << std::setw(13) << fmtSci(delta, 4) << std::setw(13) << fmtSci(radius, 4)
               << std::setw(6) << r.nCut << std::setw(7) << r.nGhost << std::setw(7) << r.ndofs
               << std::setw(15) << fmtSci(r.lmin, 3) << std::setw(15) << fmtSci(r.lmax, 3)
               << std::setw(11) << fmtCond(r.lmin, r.lmax)
               << std::setw(13) << (r.finite ? fmtSci(r.l2, 3) : std::string("diverged"))
               << std::setw(13) << (r.finite ? fmtSci(r.h1, 3) : std::string("n/a"))
               << std::setw(13) << fmtSci(r.minRowNorm, 3)
               << std::setw(11) << (r.finite && r.nGhost > 0 ? fmtSci(r.Eghost, 3) : std::string("n/a"))
               << "\n";
    }

    if (cliffFound)
        gsInfo << "Classification cliff: sliver detection lost at delta/h = " << fmtSci(cliffRatio, 2)
               << " (first ratio whose sliver cell was not classified cut; below it the sliver "
               << "silently vanishes and kappa looks good for the wrong reason).\n";
    else
        gsInfo << "Classification cliff: not reached within the swept range.\n";

    gsInfo << faceShiftLegend << "\n";
}

/// Refinement-rate study at a level-independent geometry: no
/// condition numbers are ever computed inside this loop (O(n^3) dense
/// eigensolves are affordable only at the sweep's single pinned level).
void runConvergence(const std::string & armName, index_t degree, index_t numRefine,
                     real_t radius, real_t gammaEff, real_t gtEff, index_t quRule,
                     const gsMultiPatch<real_t> & mp, const gsBoundaryConditions<real_t> & bc,
                     const gsFunctionExpr<real_t> & u_exact, const gsFunctionExpr<real_t> & f_rhs,
                     index_t N0, real_t c0, bool plot, const std::string & outPath)
{
    std::ostringstream phiStr; phiStr << std::setprecision(17)
        << "sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)-" << radius;
    std::ostringstream nxStr; nxStr << std::setprecision(17)
        << "(x-" << c0 << ")/sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)";
    std::ostringstream nyStr; nyStr << std::setprecision(17)
        << "(y-" << c0 << ")/sqrt((x-" << c0 << ")^2+(y-" << c0 << ")^2)";
    gsFunctionExpr<real_t> phi(phiStr.str(), 2);          // outlives every trim built below it
    gsFunctionExpr<real_t> nImm(nxStr.str(), nyStr.str(), 2);

    gsMultiBasis<real_t> dbasis(mp, true);
    dbasis.setDegree(degree);
    dbasis.uniformRefine(N0 - 1);

    gsVector<real_t> l2err(numRefine + 1), h1err(numRefine + 1);

    gsInfo << "\n=== Convergence  arm=" << armName << "  radius=" << radius
           << "  degree=" << degree << " ===\n";
    gsInfo << std::right << std::setw(4) << "r" << std::setw(6) << "N" << std::setw(8) << "dofs"
           << std::setw(13) << "L2err" << std::setw(13) << "H1err" << "\n";

    for (index_t r = 0; r <= numRefine; ++r)
    {
        const index_t N = N0 << r;
        real_t hmax = 0;
        for (size_t p = 0; p != dbasis.nBases(); ++p)
            hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());

        const bool doPlot = (plot && r == numRefine);
        RunResult res = solveOne(armName, mp, dbasis, phi, nImm, u_exact, f_rhs, bc,
                                  degree, N, /*ratio (n/a)*/ -1, hmax, gammaEff, gtEff, quRule,
                                  /*computeSpectrum*/ false, doPlot, outPath);

        l2err[r] = res.l2;
        h1err[r] = res.h1;

        gsInfo << std::setw(4) << r << std::setw(6) << N << std::setw(8) << res.ndofs
               << std::setw(13) << fmtSci(res.l2, 3) << std::setw(13) << fmtSci(res.h1, 3) << "\n";

        if (r != numRefine) dbasis.uniformRefine();   // refine at the END of the iteration
    }

    gsInfo << faceShiftLegend << "\n";

    gsInfo << "\nL2 error: " << std::scientific << std::setprecision(3)
           << l2err.transpose() << "\n";
    gsInfo << "H1 error: " << std::scientific << h1err.transpose() << "\n";
    if (numRefine > 0)
    {
        gsInfo << "\nEoC (L2): " << std::fixed << std::setprecision(2)
               << (l2err.head(numRefine).array() / l2err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
        gsInfo << "EoC (H1): " << std::fixed << std::setprecision(2)
               << (h1err.head(numRefine).array() / h1err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
    }
}

int main(int argc, char *argv[])
{
    std::string arm   = "stabilized";
    std::string study = "both";
    index_t numRefine   = 4;
    index_t sweepRefine = 3;
    index_t degree      = 2;
    real_t  radius       = 0.4;
    real_t  gamma         = -1;
    real_t  ghostCoef     = -1;
    index_t quRule       = gsQuadrature::AlgoimRule;
    bool    plot          = false;
    std::string outFolder = "output_ghost_penalty";

    gsCmdLine cmd("Immersed Poisson with three stabilization arms: a sliver conditioning "
                  "sweep and a convergence study.");
    cmd.addString("",     "arm",           "Arm: baseline | control | stabilized", arm);
    cmd.addString("",     "study",         "Study: sweep | convergence | both", study);
    cmd.addInt   ("r",    "uniformRefine", "Refinement levels of the convergence study (r = 0..r)", numRefine);
    cmd.addInt   ("",     "sweepRefine",   "Pinned refinement level of the sliver sweep (N = 4<<level)", sweepRefine);
    cmd.addInt   ("",     "degree",        "Spline degree of the background basis", degree);
    cmd.addReal  ("",     "radius",        "Level-set radius of the convergence study", radius);
    cmd.addReal  ("g",    "gamma",         "Nitsche penalty coefficient (<=0: auto 6*(k+1)^2)", gamma);
    cmd.addReal  ("",     "ghost",         "Ghost-penalty coefficient gamma-tilde (<0: auto 10^(-k-1))", ghostCoef);
    cmd.addInt   ("q",    "quRule",        "Immersed volume rule: 11=CutCell, 12=Algoim", quRule);
    cmd.addSwitch("plot", "Create ParaView output at the finest convergence level", plot);
    cmd.addString("o",    "output",        "Output folder", outFolder);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if ("baseline" != arm && "control" != arm && "stabilized" != arm)
    { gsWarn << "--arm must be baseline, control or stabilized\n"; return EXIT_FAILURE; }

    if ("sweep" != study && "convergence" != study && "both" != study)
    { gsWarn << "--study must be sweep, convergence or both\n"; return EXIT_FAILURE; }

    // A negative degree faults inside gsBasis::setDegree rather than being
    // rejected there, so it must not reach it.
    if (degree < 1)
    { gsWarn << "--degree must be >= 1\n"; return EXIT_FAILURE; }

    if (quRule != gsQuadrature::AlgoimRule && quRule != gsQuadrature::CutCellRule)
    {
        gsWarn << "Unsupported quRule=" << quRule << ". Falling back to AlgoimRule (12).\n";
        quRule = gsQuadrature::AlgoimRule;
    }

    // Fixed constants (not CLI options -- see file header on why alpha in
    // particular must stay hard-coded).
    const index_t N0 = 4;
    const real_t  c0 = 0.5 + 1e-3;

    // The convergence study needs a genuine immersed boundary. A non-positive
    // radius puts the level set outside [0,1]^2 entirely, so every DOF is
    // eliminated and the empty system divides by zero in the error norms; a
    // radius reaching the farthest corner encloses the whole square, leaving no
    // Dirichlet data and a singular pure-Neumann system that solves to garbage.
    const real_t rMax = math::sqrt(2.0) * math::max(c0, 1.0 - c0);
    if (radius <= 0 || radius >= rMax)
    {
        gsWarn << "--radius must lie in (0, " << rMax << "): outside that interval the level "
                  "set either misses [0,1]^2 or encloses it, leaving no immersed boundary.\n";
        return EXIT_FAILURE;
    }

    const real_t gammaEff = (gamma <= 0) ? 6.0 * (degree + 1) * (degree + 1) : gamma;
    const real_t gtEff    = (ghostCoef < 0) ? math::pow(10.0, -(degree + 1)) : ghostCoef;

    // Geometry-independent manufactured solution, valid for every radius.
    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare());
    gsFunctionExpr<real_t> u_exact("sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<real_t> f_rhs("2*pi^2*sin(pi*x)*sin(pi*y)", 2);

    gsBoundaryConditions<real_t> bc;
    bc.setGeoMap(mp);

    std::string outPath = outFolder;
    if (outPath.empty()) outPath = "output_ghost_penalty";
    const bool isAbsolutePath = (!outPath.empty() && outPath[0] == '/');
    if (!isAbsolutePath) outPath = gsFileManager::getCurrentPath() + "/" + outPath;
    const std::string out = gsFileManager::getCanonicRepresentation(outPath);
    if (plot) gsFileManager::mkdir(out);

    if ("sweep" == study || "both" == study)
        runSweep(arm, degree, sweepRefine, gammaEff, gtEff, quRule, mp, bc, u_exact, f_rhs, N0, c0);

    if ("convergence" == study || "both" == study)
        runConvergence(arm, degree, numRefine, radius, gammaEff, gtEff, quRule, mp, bc,
                        u_exact, f_rhs, N0, c0, plot, out);

    return EXIT_SUCCESS;
}
