/** @file l2projection_rh_schedule_example.cpp

    @brief Schedule-driven time-dependent L2 projection with rh-adaptivity.

    One of four "rh-driver-unification" drivers (poisson, l2projection,
    fitting, lrfitting) that share a CLI core and a schedule alphabet, but
    implement only the letters their problem supports.

    | letter | meaning                                    | poisson | l2 | fitting | lrfitting |
    |--------|---------------------------------------------|---------|----|---------|-----------|
    | S      | solve / project / fit (the primal step)      | yes     | yes| yes     | yes       |
    | R      | monitor-driven relocation of sigma (Winslow) | yes     | yes| --      | --        |
    | D      | error-driven relocation of sigma             | --      | yes| yes     | yes       |
    | H      | local refinement (THB / LR)                  | yes     | yes| yes     | yes       |
    | U      | uniform refinement                           | yes     | yes| yes     | --        |

    This driver implements S, R, D, H, U. Its primal-step letter is 'S' (not
    the older 'P'), and 'P'/'F' are NOT accepted as aliases -- any letter
    outside S,R,D,H,U is rejected by name.

    D19 -- SURFACE by default. The geometry (-f) is read from a template XML
    in the filedata/parametrization convention (shipped default: a spherical
    patch), the target function f is generated in code (RingMonitor, below)
    in AMBIENT coordinates, and the R step therefore runs with
    parametric=false. -f "" selects the planar unit square [-0.5,0.5]^2
    (dim=2) -- the historical merging-rings case, and the only configuration
    comparable with the archived pre-D19 baseline.

    CLI flags:
      -f/--file       template geometry XML (filedata/parametrization
                       convention); "" selects the planar fallback. Default
                       "parametrization/monitor_example_sphere_patch_with_function.xml".
      --options       method option list XML (gsOptionList), layered over the
                       built-in defaultDriverOptions() (partial lists stay
                       valid). "-f" carries the PROBLEM only, never options
                       (ORCHESTRATOR RULING 1).
      -o/--output     output prefix/directory
      --plot          write the ParaView output set (D12)
      -p/--degree     analysis polynomial degree
      -r/--refine     initial uniform refinements
      -E/--sigmaDeg   degree of the sigma map
      -R/--sigmaRef   refinement LEVEL of sigma's mesh (2^sigmaRef elements
                       per direction, dyadically nested with every analysis
                       level -- see the nesting comment at sigma's construction)
      --schedule      letters S,R,D,H,U (default "SRS")
      -i/--iter       schedule repetitions PER TIME STEP (hard cap)
      --project       also project onto the L2-projected composed geometry
                       at every S step (D8)
      -n/--nsteps     number of time steps (this driver only, D11)

    -d/--data is deliberately NOT declared (plan D18e): this driver has no
    external data set, only the analytic RingMonitor.

    XML options (defaultDriverOptions(), layered with --options). This is one
    of four rh drivers sharing ONE unified 27-key option table
    (rh-driver-unification task 18); DirSkip and Lambda are declared but
    IGNORED here purely for key-set parity across the four reference XMLs:
      RefineRule, RefineParam, CoarsenRule, CoarsenParam, Coarsen, MaxLevel,
      Admissible, Extension          -- gsHElementMarker (H step)
      Smoothing, Penalty             -- R step (Winslow smoothing / fold penalty)
      Target, Band                   -- H step error band (D9)
      MaxRefIt                       -- inner adapt iterations per H step
      Optimizer, MaxIterations, OptTol, Verbose -- R/D-step optimizer backend
      MonitorMode                    -- R step monitor: "value" | "gradient"
      BarrierMu, BarrierEps          -- D step fold-barrier weight/floor (the
                                         values formerly hardcoded as
                                         dMu=1000.0/dEps=5e-2, renamed to keys)
      DirSkip                        -- IGNORED, not implemented in this
                                         driver's D step
      Lambda                         -- IGNORED, no gsFitting in this driver
      Slide                          -- sigma boundary controls may slide
      quA, quB, SameElement          -- IGA error/assembly quadrature, mirrored
                                         onto every element-wise error sweep
                                         (ORCHESTRATOR ADDENDUM, task 04's
                                         review: the default quA=1,quB=1
                                         element sweep under-integrates by
                                         ~39%). Now defaults to 2.0/2, the
                                         converged setting (task 18, defect 2)
                                         -- a DIFFERENT quantity from the
                                         fitting drivers' quB (D-step barrier
                                         Gauss order on sigma's knot mesh)
      RingRadius, RingWidth, RingCx0, RingAxis -- the moving test function's
                                         shape (D19: the ring geometry knobs
                                         belong to the options XML, not the CLI)

    D9 -- error-band rule for the H step, driven by the GLOBAL L2
    error of the last S step:

      err >  Target*Band                   -> refine only, keep going
      Target/Band < err <= Target*Band     -> ONE refine+coarsen sweep, STOP
      err <= Target/Band                   -> coarsen only, keep going

    The band is a target zone the H step drives the error INTO, from either
    side. Reaching it stops the CURRENT TIME STEP (not the outer ts loop --
    the RingMonitor moves between time steps, D19.3, so convergence at ts
    says nothing about ts+1; see the `stepConverged` declaration at the top
    of the ts loop).

    2026-08-18 CHANGE: this header always documented an early stop, but the
    implementation computed `inBandStop` and then discarded it
    (`(void)inBandStop`), so the driver never stopped early at all. Fixed to
    match the documented rule.

    2026-08-20 CHANGE: the dead zone used to refine as well, and only a
    fourth branch below it ("hold", at or under NoCoarsenBelow =
    Target/Band^2) stopped the step -- so it had to overshoot its requested
    Target by a factor Band^2 before it was allowed to finish, and on a case
    whose error floor lies inside the band it never did. NoCoarsenBelow no
    longer takes part in the decision, and coarsening follows the branch
    rather than running on every H step: above the band there is nothing to
    give back, below it refinement is not what is called for, and inside it
    both act exactly once.

    Target < 0 disables the band entirely (always refine, coarsen iff
    Coarsen). -i/--iter is
    a hard cap on schedule repetitions regardless of the band; MaxRefIt caps
    the inner adapt loop of a single H letter.

    D13 -- <output>/convergence.csv, one row per EXECUTED schedule letter,
    frozen header (identical in all four drivers):

      cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time

    cycle = timeStep * iterations + rep (so timeStep = cycle / iterations,
    rep = cycle % iterations -- the frozen schema has no separate time-step
    column). path is 'C' (composed) or 'P' (projected) on S rows only;
    non-primal rows (R/D/H/U) are not tied to a solve path and carry the
    literal '-'. err = global L2 error of
    the row's own S step (or the last known one for R/D/H/U rows);
    minErr/maxErr = min/max per-ELEMENT L2 error from
    gsExprEvaluator::integralElWise (D2); pctBelowTol (ORCHESTRATOR ADDENDUM,
    task 07's review) = percentage of elements whose per-element error is
    below Target/sqrt(nElements) -- the per-element share of a global target
    under equidistribution, 0 when Target<0; minDetJsigma =
    sigma.minJacobian(), RECOMPUTED on every row.

    Column vintage (task 20 -- state coherence). minErr/maxErr/err/
    pctBelowTol on EVERY row describe the state AFTER that row's own step,
    never a mix: they are all written by ONE update path (updateState() in
    the driver source), called at the end of every step that changes the
    carried solution/error/mesh:
      - an S step (C, or under --project the following P) -- a fresh
        projection, the source of truth for the others;
      - R's post-relocation re-measurement of the frozen solution against the
        new sigma (review fix 1: R moves the composed geometry S o sigma
        exactly as D does, so it needs the identical re-measurement);
      - D's post-relocation re-measurement of the frozen solution (task 19);
      - U's re-measurement of the frozen solution against the refined basis
        (review fix 3: uniformRefine() changes the mesh the carried state was
        measured on, same rationale as R/D);
      - the H band's inner re-projection between refinement iterations, AND
        (review fix 2) a terminal re-measurement of the frozen solution
        against the FINAL mesh once the H letter's inner loop exits, on any
        exit path that changed the mesh (normal exhaustion, the `!changed`
        break, or the `inBand` break) -- MaxRefIt defaults to 1, so in the
        default configuration this terminal refresh is the ONLY update the H
        branch performs.
    In particular a D (or R, or U) row's four error columns are the
    POST-step values (not the pre-step ones the step started from), and
    lastErr is refreshed for real -- not merely substituted into that one CSV
    row -- so the NEXT row (typically an H band decision, D9) also reads the
    post-step error, not a stale pre-step one.

    --project geometry convention on a D row. A D step's re-measurement
    always evaluates the frozen lastSol against the COMPOSED geometry
    S o sigma (the only geometry sigma's own relocation can possibly
    affect). Once an 'S,P' row has run under --project, lastSol is the P
    path's solution living on a FROZEN, ordinary (non-composed) geometry
    that a later D step's sigma relocation cannot touch at all -- so in that
    specific case (--project, D immediately after an S,P row) a D row's err
    is not a like-for-like before/after number: it mixes a P-path solution
    with a C-path (composed) measure. This is a genuine convention clash,
    not a bookkeeping bug this driver can clean up on its own; a harvester
    comparing a --project run's D row against its preceding P row (e.g.
    rh_harvest.py's `drops`) must account for it -- fixing rh_harvest.py
    itself is out of scope here.

    Under --plot the driver writes, in <output>: `geometry` (the final
    composed map), `sigma` (the sigma map itself), `solution` (one
    gsParaviewCollection per path -- `solution` alone without --project,
    `solution_C`/`solution_P` with it) and `mesh` (the final analysis mesh).
    <output>/options.xml and <output>/convergence.csv are written
    UNCONDITIONALLY.

    MaxIterations DEFAULT (task-08 measurement, D-step of THIS driver):
    10000 (the historical schedule-driver default) makes the D step a
    guaranteed REJECT-AND-RESTORE NO-OP costing 170-290s, on both surface and
    planar geometries -- the long run wanders into a fold, the reject-and-
    restore policy then discards the whole step, and minutes are burned to
    change nothing. 50 converges in ~1.7s with a genuine ~34% error drop.
    More iterations make the result STRICTLY WORSE here, so the shipped
    default is 100 (top of the sanctioned 50-100 range) and ONE shared key
    drives both the R and the D step (no MaxIterationsD). Measured on THIS
    driver (not the poisson baseline the D-step numbers above come from): the
    R step typically REACHES the 100 cap rather than converging under it, so
    MaxIterations also truncates the R step here -- raise it on the command
    line when a fully converged R step is wanted. This is an ASYMMETRY with
    the poisson driver, which keeps MaxIterations
    at its historical 10000: poisson has no D step (the reject-and-restore
    no-op measured here does not apply to it), and its archived R-step
    baseline was produced with 10000, so changing it there would break the
    bit-identity comparison against that baseline.

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
#include <gsOptim/gsOptim.h>
#include <gsUtils/gsProjection.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <memory>
#include <set>
#include <sstream>
#include <type_traits>

using namespace gismo;

namespace {

// ---------------------------------------------------------------------------
// Two-merging-rings test function, generalized from domainDim 2 (D19.3) to
// domainDim `dim` (2 or 3): f lives in AMBIENT coordinates, matching what the
// filedata/parametrization XMLs already do. This is the single definition of
// f used by all three steps that need it (projection source, R-step monitor,
// error/exact field) -- see the class-level warning at deriv2_into.
//
// Per ring i (centre c_i, radius R, front width w), with local offset
// e = x - c_i (size n = dim):
//   d = max(||e||, 1e-14)  -- floor: a, and the 1/d terms in the Hessian, are
//                             singular at the ring centre. Quadrature points
//                             never land exactly on it, but the floor turns a
//                             would-be NaN into a large finite number that an
//                             isfinite() guard downstream can act on.
//   s = (d-R)/w,  t = tanh(s),  r = 1 - t^2   (annulus bump, r=1 on d=R)
//   a_i = e_i/d,  k = -2 r t / w,  dk = (-6 r^2 + 4 r) / w^2
//   dr/dx_i      = k * a_i
//   d2r/dx_i dx_j = dk * a_j * a_i + k * (delta_ij - a_i*a_j) / d
// For n=2 these reduce EXACTLY to the historical 2D formulas
// (rxx = dkdx*a + k*(1-a*a)/d, rxy = dkdy*a - k*a*b/d).
//
// Soft union of the two rings (stays in [0,1] even where the rings overlap):
//   f    = r1 + r2 - r1*r2
//   f_i  = r1_i (1-r2) + r2_i (1-r1)
//   f_ij = r1_ij(1-r2) + r2_ij(1-r1) - r1_i r2_j - r1_j r2_i
// (for i==j the last two terms collapse to -2 r1_i r2_i).
// ---------------------------------------------------------------------------
class RingMonitor : public gsFunction<real_t>
{
public:
    // \a c1, \a c2 are the two ring centres at the current time (ambient,
    // size dim).
    RingMonitor(const gsVector<real_t> & c1, const gsVector<real_t> & c2,
                real_t R, real_t w)
    : m_c1(c1), m_c2(c2), m_R(R), m_w(w)
    { GISMO_ENSURE(c1.size()==c2.size(), "ring centres must have the same size"); }

    GISMO_CLONE_FUNCTION(RingMonitor)

    short_t domainDim() const override { return static_cast<short_t>(m_c1.size()); }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        Ring g1(domainDim()), g2(domainDim());
        for (index_t p = 0; p != u.cols(); ++p)
        {
            ring(u.col(p) - m_c1, g1, false);
            ring(u.col(p) - m_c2, g2, false);
            result(0,p) = g1.r + g2.r - g1.r * g2.r;
        }
    }

    void deriv_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        const short_t n = domainDim();
        result.resize(n, u.cols());
        Ring g1(n), g2(n);
        for (index_t p = 0; p != u.cols(); ++p)
        {
            ring(u.col(p) - m_c1, g1, true);
            ring(u.col(p) - m_c2, g2, true);
            for (short_t i = 0; i != n; ++i)
                result(i,p) = g1.grad[i] * (1 - g2.r) + g2.grad[i] * (1 - g1.r);
        }
    }

    // Row order follows gsFunctionSet's documented convention (pure
    // derivatives first, then the upper triangle lexicographically):
    // n=2 -> [f_xx,f_yy,f_xy] (3 rows), n=3 -> [f_xx,f_yy,f_zz,f_xy,f_xz,f_yz]
    // (6 rows). Built via an index loop (not a hand-written list) so the
    // n=2/n=3 cases cannot drift.
    void deriv2_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        const short_t n = domainDim();
        std::vector<std::pair<short_t,short_t> > idx;
        for (short_t i = 0; i != n; ++i) idx.push_back(std::make_pair(i,i));
        for (short_t i = 0; i != n; ++i)
            for (short_t j = i+1; j != n; ++j)
                idx.push_back(std::make_pair(i,j));

        result.resize((index_t)idx.size(), u.cols());
        Ring g1(n), g2(n);
        for (index_t p = 0; p != u.cols(); ++p)
        {
            ring(u.col(p) - m_c1, g1, true);
            ring(u.col(p) - m_c2, g2, true);
            for (size_t k = 0; k != idx.size(); ++k)
            {
                const short_t i = idx[k].first, j = idx[k].second;
                real_t val = g1.hess(i,j) * (1 - g2.r) + g2.hess(i,j) * (1 - g1.r)
                           - g1.grad[i]*g2.grad[j] - g1.grad[j]*g2.grad[i];
                if (i==j) val = g1.hess(i,j) * (1 - g2.r) + g2.hess(i,j) * (1 - g1.r)
                               - 2 * g1.grad[i] * g2.grad[i];
                result((index_t)k, p) = val;
            }
        }
    }

private:
    struct Ring
    {
        explicit Ring(short_t n) : grad(n), hess(n,n), r(0) { grad.setZero(); hess.setZero(); }
        gsVector<real_t> grad;
        gsMatrix<real_t> hess;
        real_t r;
    };

    // One ring bump at local offset e = x - c. Fills the value always and
    // the derivatives only when \a wantDerivs.
    void ring(const gsVector<real_t> & e, Ring & g, bool wantDerivs) const
    {
        const short_t n = static_cast<short_t>(e.size());
        const real_t d = std::max(e.norm(), real_t(1e-14));
        const real_t s = (d - m_R) / m_w;
        const real_t t = math::tanh(s);
        g.r = 1 - t*t;
        if (!wantDerivs) return;

        gsVector<real_t> a = e / d;
        const real_t k  = -2 * g.r * t / m_w;
        const real_t dk = (-6 * g.r * g.r + 4 * g.r) / (m_w * m_w);
        for (short_t i = 0; i != n; ++i)
            g.grad[i] = k * a[i];
        for (short_t i = 0; i != n; ++i)
            for (short_t j = 0; j != n; ++j)
                g.hess(i,j) = dk * a[j] * a[i]
                            + k * ((i==j ? real_t(1) : real_t(0)) - a[i]*a[j]) / d;
    }

    gsVector<real_t> m_c1, m_c2;
    real_t m_R, m_w;
};

// ---------------------------------------------------------------------------
// Driver default options (D6). Names/defaults are frozen by the task spec.
// ---------------------------------------------------------------------------
gsOptionList defaultDriverOptions()
{
    gsOptionList opt;
    opt.addInt   ("RefineRule",     "gsHElementMarker refinement rule (1=GARU,2=PUCA,3=BULK)", 3);
    opt.addReal  ("RefineParam",    "gsHElementMarker refinement parameter", 0.5);
    opt.addInt   ("CoarsenRule",    "gsHElementMarker coarsening rule", 3);
    opt.addReal  ("CoarsenParam",   "gsHElementMarker coarsening parameter", 0.1);
    opt.addSwitch("Coarsen",        "enable coarsening in the H step", false);
    opt.addInt   ("MaxLevel",       "gsHElementMarker level cap (library default 3 is far too low here)", 10);
    opt.addSwitch("Admissible",     "gsHElementMarker admissible closure", true);
    opt.addSwitch("Extension",      "gsHElementMarker box extension", true);
    opt.addReal  ("Smoothing",      "R step: Winslow smoothing theta", 1.0);
    opt.addReal  ("Penalty",        "R step: fold penalty", 1e-3);
    opt.addReal  ("Target",         "H band: target L2 error (<0 disables the band)", -1);
    opt.addReal  ("Band",           "H band: dead-zone width factor (>=1)", 2.0);
    opt.addReal  ("NoCoarsenBelow", "H band: UNUSED -- the band decision no "
                                    "longer reads it (2026-08-20)", 0.0);
    opt.addInt   ("MaxRefIt",       "H band: inner adapt iterations per H step", 1);
    opt.addString("Optimizer",      "R/D-step optimizer backend: gsOptim | HLBFGS", "gsOptim");
    // task-08 measurement (l2 driver's D step): MaxIterations=10000 (the
    // historical schedule-driver default) makes the D step a guaranteed
    // reject-and-restore NO-OP costing 170-290s on both surface and planar
    // geometries -- the long run wanders into a fold, then the whole step is
    // discarded; MaxIterations=50 converges in ~1.7s with a genuine ~34%
    // error drop. More iterations make the result WORSE here, so the shipped
    // default is 100 (top of the sanctioned 50-100 range, matching the two
    // fitting drivers' historical default). This key drives BOTH the R and
    // the D step (one shared key, see makeOptimizer() below). Measured on
    // THIS driver: the R step typically REACHES the 100 cap rather than
    // converging under it (the poisson driver's "~20 iterations" baseline is
    // a different geometry/monitor/objective and does not transfer), so this
    // key also truncates the R step here -- raise it on the command line
    // when a fully converged R step is wanted.
    //
    // NOTE the resulting poisson-vs-l2 asymmetry: the poisson driver's own
    // MaxIterations default stays 10000 (its archived R-step baseline was
    // produced with it, and the bit-identity comparison against that
    // baseline depends on it) -- poisson has no D step, so the reject-and-
    // restore no-op measured here does not apply there.
    opt.addInt   ("MaxIterations",  "R/D-step optimizer iterations (see the file "
                                     "header for why a large value is "
                                     "counter-productive on the D step)", 100);
    opt.addReal  ("OptTol",         "R/D-step optimizer gradient tolerance (unified name; "
                                     "NEVER GradTol -- gsHLBFGS already owns a GradTol that is "
                                     "the line-search curvature constant, default 0.9, a "
                                     "different quantity)", 1e-4);
    opt.addInt   ("Verbose",        "R/D-step optimizer verbosity", 0);
    opt.addString("MonitorMode",    "R step monitor: value | gradient. LEGITIMATELY differs "
                                     "from poisson's \"value\" default: this driver supplies an "
                                     "exact analytic monitor (RingMonitor, with value/gradient/"
                                     "Hessian), which is what GradientBased requires; poisson "
                                     "has only a discrete solution", "gradient");
    opt.addSwitch("Slide",          "sigma boundary control points may slide", true);
    // quA=2/quB=2 is the CONVERGED error/assembly quadrature (task 08's review
    // measured quA=1/quB=1 to under-integrate the err column by ~4.5%, task
    // 04's review by ~39% on the default element sweep) -- rh-driver-
    // unification task 18, defect 2.
    opt.addReal  ("quA",            "quadrature points-per-direction factor "
                                     "(error evaluators + gsAdaptiveParametrization)", 2.0);
    opt.addInt   ("quB",            "quadrature extra points per direction -- IGA error/assembly "
                                     "quadrature, a DIFFERENT quantity from the fitting drivers' "
                                     "quB (the D-step fold-barrier Gauss order on sigma's own "
                                     "knot mesh)", 2);
    opt.addSwitch("SameElement",    "quadrature: integration/analysis elements assumed to "
                                     "coincide -- must stay false under the sigma composition "
                                     "(a quadrature element does not map into a single analysis "
                                     "element)", false);
    opt.addReal  ("BarrierMu",      "D-step fold-barrier weight (the value formerly hardcoded "
                                     "as dMu, unrenamed except for the key -- do not confuse "
                                     "with the R-step 'Penalty')", 1000.0);
    opt.addReal  ("BarrierEps",     "D-step fold-barrier floor (the value formerly hardcoded "
                                     "as dEps)", 5e-2);
    opt.addInt   ("BarrierMode",    "D-step fold barrier: 0 = sampled Gauss rule on sigma's "
                                     "knot mesh (quB), 1 = Bezier-coefficient barrier (no quadrature)", 0);
    opt.addInt   ("DirSkip",        "D-step direction skip -- IGNORED, not implemented in this "
                                     "driver's D step", 0);
    opt.addReal  ("Lambda",         "gsFitting smoothing weight -- IGNORED, no gsFitting in "
                                     "this driver", 1e-6);
    opt.addReal  ("RingRadius",     "ring radius R", 0.25);
    opt.addReal  ("RingWidth",      "ring front width w", 0.05);
    opt.addReal  ("RingCx0",        "initial ring centre offset from the patch midpoint", 0.30);
    opt.addInt   ("RingAxis",       "ambient axis the two centres separate along: 0=x,1=y,2=z", 0);
    return opt;
}

// Both backends are reachable; the unified names (MaxIterations/OptTol/
// Verbose) are mapped onto each backend's own so the mismatch is visible in
// the log instead of hidden in a header. The unified tolerance key is
// deliberately called "OptTol", not "GradTol" -- gsHLBFGS already OWNS an
// option literally named "GradTol" (the line-search curvature constant,
// default 0.9), a completely different quantity (ORCHESTRATOR RULING 2).
std::unique_ptr<gsOptimizer<real_t> > makeOptimizer(const gsOptionList & o)
{
    const std::string which = o.askString("Optimizer","gsOptim");
    std::unique_ptr<gsOptimizer<real_t> > optimizer;
    if (which=="HLBFGS")
    {
        optimizer.reset(new gsHLBFGS<real_t>());
        optimizer->options().setReal("MinGradLen", o.getReal("OptTol"));
    }
    else if (which=="gsOptim")
    {
        optimizer.reset(new gsOptim<real_t>::LBFGS());
        optimizer->options().setReal("GradErrTol", o.getReal("OptTol"));
    }
    else
        GISMO_ERROR("Unknown Optimizer '"<<which<<"' (gsOptim | HLBFGS)");
    optimizer->options().setInt("MaxIterations", o.getInt("MaxIterations"));
    optimizer->options().setInt("Verbose",       o.getInt("Verbose"));
    return optimizer;
}

// One S-step result: the extracted (scalar) projection, its global L2 error
// (from the SAME evaluator/quadrature as the per-element sweep, D2), and the
// per-element L2 error contributions on the ANALYSIS basis (needed by both
// the CSV min/max columns and the H-step marker).
struct ProjResult
{
    gsMultiPatch<> sol;
    real_t err = 0.0;
    std::vector<real_t> elErr; // per element, on active's elements, sqrt'd
};

// Per-element and global L2 error of an ALREADY-COMPUTED solution \a sol
// against \a monitor over the composed geometry \a cmp, using integration
// elements \a mb and analysis basis \a active -- the D2 sweep
// (gsExprEvaluator::integralElWise), factored out of doProject so both a
// freshly projected solution (doProject below) and a re-measurement of a
// FROZEN solution (task 19/20's D step: no re-projection, plan D5) go
// through the identical quadrature/convention instead of two hand-written
// integrals that can silently drift apart.
ProjResult measureError(
    const gsMultiPatch<>& cmp, const gsMultiBasis<>& mb, const gsBasis<>& active,
    const gsMultiPatch<>& sol, const gsFunction<real_t>& monitor, const gsOptionList& opt)
{
    ProjResult res;
    res.sol = sol;

    gsExprEvaluator<> eev;
    eev.options().setReal ("quA", opt.getReal("quA"));
    eev.options().setInt  ("quB", opt.getInt ("quB"));
    // SameElement=false: a quadrature element of the integration basis does
    // not map into a single element of the (composed) analysis space.
    eev.options().addSwitch("SameElement","",opt.getSwitch("SameElement"));
    eev.setIntegrationElements(mb);
    auto G  = eev.getMap(cmp);
    auto uh = eev.getVariable(res.sol);
    auto fe = eev.getVariable(monitor, G);
    const real_t err2 = eev.integralElWise((fe - uh).sqNorm() * meas(G));
    const std::vector<real_t> & elw = eev.elementwise();
    GISMO_ENSURE((index_t)elw.size() == active.numElements(),
                 "element error count does not match the analysis basis");
    res.elErr.assign(active.numElements(), 0.0);
    index_t k = 0;
    for (auto it = active.domain()->beginAll(); it != active.domain()->endAll(); ++it, ++k)
        res.elErr[it.id()] = math::sqrt(elw[k]);
    res.err = math::sqrt(err2);
    return res;
}

// Project \a monitor onto \a mb using the composed geometry \a cmp with
// integration \a ib, then compute its L2 error via integralElWise on the
// ANALYSIS basis \a active (D2: gsExprEvaluator::integralElWise, not a
// hand-rolled sample grid). \a pev is a persistent, sharing evaluator; when
// \a collection is given the current time step is appended to it (its
// integration elements are (re)set to \a mb here so plotElements shows the
// analysis mesh, not the union integration mesh).
ProjResult doProject(
    const gsMultiPatch<>& cmp, const gsMultiBasis<>& mb, const gsMultiBasis<>& ib,
    const gsBasis<>& active, const gsFunction<real_t>& monitor,
    const gsOptionList& opt, real_t t, index_t step,
    gsExprEvaluator<>* pev, gsParaviewCollection* collection, const char* tag)
{
    // gsL2Projection registers the source with the PARAMETRIC getCoeff
    // overload ("auto f = A.getCoeff(sourceFunction)", gsProjection.hpp), so
    // it expects a function of the reference variable v -- NOT of the
    // physical point. Handing it `monitor` directly would project the rings
    // as read on the reference square, which for the surface default is not
    // even the same space as monitor's ambient R^3 domain. Compose with the
    // geometry so the source really is f(G(v)).
    const gsComposedFunction<real_t> funcOnParam(cmp.patch(0), monitor);
    gsMatrix<> coefs;
    const real_t projErr = gsL2Projection<real_t>::project(mb, ib, cmp, funcOnParam, coefs);
    coefs.resize(mb.size(), 1);

    gsGeometry<>::uPtr solGeom = mb.basis(0).makeGeometry(give(coefs));
    gsMultiPatch<> solmp; solmp.addPatch(*solGeom);

    // D2: per-element L2 error via integralElWise on the ANALYSIS basis
    // (not the union integration basis -- gsHElementMarker::setErrors
    // indexes by the id() of the basis it was constructed on). ORCHESTRATOR
    // ADDENDUM (task 04's review): quA/quB set EXPLICITLY here, mirroring
    // the resolved driver options -- the default quA=1,quB=1 element sweep
    // was measured to under-integrate by ~39%. Routed through measureError()
    // (task 20) so this sweep and the D step's post-relocation re-measurement
    // are the exact same code.
    ProjResult res = measureError(cmp, mb, active, solmp, monitor, opt);

    gsInfo << "  " << std::left << std::setw(2) << tag << std::right << "| step ";
    if (step < 0) gsInfo << std::setw(3) << "-";
    else          gsInfo << std::setw(3) << step;
    gsInfo << " | t=" << std::fixed << std::setprecision(3) << t
           << " | dofs " << std::setw(6) << mb.totalSize()
           << " | L2 err " << std::scientific << std::setprecision(3) << res.err
           << " (gsL2Projection reference: " << projErr << ")"
           << std::defaultfloat << "\n";

    if (collection && pev)
    {
        pev->setIntegrationElements(mb);
        auto pG  = pev->getMap(cmp);
        auto pf  = pev->getVariable(monitor, pG);
        auto puh = pev->getVariable(res.sol);
        collection->newTimeStep(const_cast<gsMultiPatch<>*>(&cmp), t);
        collection->addField(puh,                 "projected solution");
        collection->addField(pf,                  "exact function");
        collection->addField((pf - puh).sqNorm(), "error");
        collection->saveTimeStep();
    }
    return res;
}

} // namespace

int main(int argc, char** argv)
{
    index_t degree     = 2;
    index_t initialRef  = 2;
    index_t nSteps      = 6;
    index_t sigmaDeg    = 2;
    index_t sigmaRef    = 3;
    index_t iterations  = 1;
    bool    project      = false;
    bool    plot         = false;
    bool    coarsenFlag  = false;
    real_t  targetCli    = std::numeric_limits<real_t>::quiet_NaN();
    std::string schedule = "SRS";
    std::string file     = "parametrization/monitor_example_sphere_patch_with_function.xml";
    std::string options  = "";
    std::string output   = "l2projection_rh_output";

    gsCmdLine cmd("Schedule-driven time-dependent L2 projection with rh-adaptivity "
                  "(default: a moving pair of rings projected onto a spherical patch).");
    cmd.addString("f","file","template geometry XML in the filedata/parametrization "
                             "convention (labels 'geometry' [required], 'composition' "
                             "[optional]); \"\" selects the built-in PLANAR unit square "
                             "[-0.5,0.5]^2 (the historical merging-rings case)",file);
    cmd.addString("","options","method option list XML (layered with addIfUnknown over "
                               "the built-in defaults; -f carries the PROBLEM only, "
                               "never options -- ORCHESTRATOR RULING 1)",options);
    cmd.addString("o","output","output prefix/directory",output);
    cmd.addSwitch("plot","write the ParaView output set",plot);
    cmd.addInt("p","degree","analysis degree",degree);
    cmd.addInt("r","refine","initial uniform refinements",initialRef);
    cmd.addInt("E","sigmaDeg","degree of the sigma map",sigmaDeg);
    cmd.addInt("R","sigmaRef","uniform refinement LEVEL of sigma's mesh w.r.t. the "
               "LEVEL-0 analysis basis (2^sigmaRef elements per direction, so sigma "
               "is always dyadically nested with every analysis level)",sigmaRef);
    cmd.addString("","schedule","letters S,R,D,H,U",schedule);
    cmd.addInt("i","iter","schedule repetitions PER TIME STEP (hard cap)",iterations);
    cmd.addSwitch("project","also project onto the L2-projected composed geometry, "
                            "at every S step",project);
    cmd.addSwitch("coarsen","H step: force Coarsen=true regardless of the XML/--options "
                            "value (omit to use whatever the option list resolves to)",
                            coarsenFlag);
    cmd.addReal("","target","H step: override the error-band Target (L2 units) "
                            "regardless of the XML/--options value; NaN (default, omit "
                            "the flag) means 'no override, use the option list'",
                            targetCli);
    cmd.addInt("n","nsteps","number of time steps",nSteps);
    // -d/--data is deliberately NOT declared (plan D18e): this driver has no
    // external data set, only the analytic RingMonitor.

    // ORCHESTRATOR RULING 2: -E/-R "given explicitly on the command line"
    // cannot be read back from gsCmdLine (no isSet()/wasSet() API), so scan
    // argv directly BEFORE cmd.getValues() consumes it. This is what lets a
    // sweep over -E/-R always win over a file's composition (see sigma
    // construction below), regardless of which value was passed.
    bool sigmaGivenOnCLI = false;
    for (int ai = 1; ai < argc; ++ai)
    {
        const std::string a = argv[ai];
        if (a=="-E" || a=="--sigmaDeg" || a=="-R" || a=="--sigmaRef")
            sigmaGivenOnCLI = true;
    }

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // -r and -R are EXPONENTS (2^r elements per direction). 1<<n is UB for
    // n>=31 and silently ALIASES for n>=32 (measured: -r 34 behaves as -r 2),
    // so a typo'd sweep exponent would otherwise produce plausible numbers for
    // the wrong mesh. GISMO_ENSURE, not GISMO_ASSERT: asserts are inert in
    // build_rel. 20 is already 2^20 elements per direction -- far past any
    // runnable size.
    GISMO_ENSURE(initialRef >= 0 && initialRef <= 20,
                 "-r/--refine must be in [0,20] (got " << initialRef << ")");
    GISMO_ENSURE(sigmaRef   >= 0 && sigmaRef   <= 20,
                 "-R/--sigmaRef must be in [0,20] (got " << sigmaRef << ")");
    GISMO_ENSURE(sigmaDeg >= 1, "-E/--sigmaDeg must be at least 1 (got " << sigmaDeg << ")");

    // Validate schedule BEFORE creating the output directory, so an invalid
    // schedule fails fast without leaving an empty directory behind
    // (D4: this driver implements S,R,D,H,U -- 'P' and 'F' are NOT accepted
    // as aliases for 'S').
    std::string sched;
    for (char c : schedule)
    {
        const char C = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        GISMO_ENSURE(C=='S'||C=='R'||C=='D'||C=='H'||C=='U',
            "l2projection_rh_schedule_example: unsupported schedule letter '" << c
            << "'. This driver implements S (project), R (monitor relocation), "
            "D (error-driven relocation), H (local THB refinement), U (uniform "
            "refinement). 'P' and 'F' are NOT accepted as aliases for 'S'.");
        sched += C;
    }
    const bool useH = sched.find('H') != std::string::npos;
    GISMO_ENSURE(nSteps > 0, "nSteps must be positive");
    GISMO_ENSURE(iterations > 0, "iterations must be positive");
    GISMO_ENSURE(sigmaRef >= 0, "sigmaRef must be non-negative");

    if (output.empty() || output.back() != gsFileManager::getNativePathSeparator())
        output += gsFileManager::getNativePathSeparator();
    // Written unconditionally (not only under --plot), so a run directory is
    // reproducible from its own output (options.xml, convergence.csv).
    gsFileManager::mkdir(output);

    // -------------------------------------------------------------------
    // D19.1/D19.2: load the template geometry (or build the planar
    // fallback), and derive the analysis space from it.
    // -------------------------------------------------------------------
    gsMultiPatch<> physical, composition;
    // id=4 -- the driver's OWN method options, added 2026-08-19 alongside
    // poisson/fitting/lrfitting's identical convention -- is read into a
    // plain gsOptionList here (NOT by keeping the gsFileData<> itself
    // alive: it has a custom destructor over pImpl-style internals, so
    // copying/reassigning a gsFileData<> after construction double-frees;
    // read out everything needed while `fd` is in scope instead) so it can
    // be layered into `opt` in the options section further down even
    // though -f's gsFileData<> goes out of scope at the end of this block.
    gsOptionList fileOpts;
    bool haveFileOpts = false;
    if (!file.empty())
    {
        gsFileData<> fd(file);
        gsInfo << "Loaded geometry file " << fd.lastPath() << "\n";
        if (fd.hasId(4)) { fd.getId(4, fileOpts); haveFileOpts = true; }
        fd.getLabel("geometry", physical);
        GISMO_ENSURE(physical.nPatches()==1,
            "l2projection_rh_schedule_example: the 'geometry' label must carry a "
            "single patch");
        if (fd.hasLabel("composition"))
            fd.getLabel("composition", composition);
        if (fd.hasLabel("function"))
            gsInfo << "Note: the file's 'function' label is IGNORED -- this driver's "
                      "target function is always the time-dependent RingMonitor "
                      "(D19), never a static function read from file.\n";
        if (fd.hasLabel("optimizer_options") || fd.hasLabel("parametrization_options"))
            gsInfo << "Note: the file's 'optimizer_options'/'parametrization_options' "
                      "are IGNORED -- method options come from --options only "
                      "(ORCHESTRATOR RULING 1).\n";
    }
    else
    {
        // Planar fallback (D19): today's unit square of degree p, shifted so
        // the domain is [-0.5,0.5]^2 -- the historical merging-rings case,
        // and the only configuration comparable with the archived baseline.
        physical.addPatch(*gsNurbsCreator<>::BSplineSquareDeg(degree));
        physical.patch(0).coefs().array() -= 0.5;
        gsInfo << "-f \"\": using the built-in PLANAR square [-0.5,0.5]^2 (dim=2) "
                  "instead of a template geometry file.\n";
    }
    GISMO_ENSURE(physical.nPatches()==1, "single-patch only");
    const short_t dim = physical.patch(0).targetDim();
    gsInfo << "Ambient dimension dim = " << dim << " (geoDim = " << dim << ")\n";

    // The analysis space stays a plain (non-rational) tensor B-spline space
    // even when the geometry itself is rational (the sphere patch carries a
    // gsTensorNurbsBasis<2>): take the source basis of a rational one, or
    // the basis itself if it is already a plain tensor B-spline basis. Both
    // cases put f in AMBIENT coordinates (below), which is what makes
    // parametric=false necessary for the R step (D19.4).
    gsTensorBSplineBasis<2> basis;
    {
        const gsBasis<real_t>& gb = physical.patch(0).basis();
        if (const gsTensorNurbsBasis<2,real_t>* nb =
                dynamic_cast<const gsTensorNurbsBasis<2,real_t>*>(&gb))
            basis = nb->source();
        else if (const gsTensorBSplineBasis<2>* bb =
                dynamic_cast<const gsTensorBSplineBasis<2>*>(&gb))
            basis = *bb;
        else
            GISMO_ERROR("expected a tensor-product (B-spline or NURBS) geometry basis");
    }
    basis.setDegree(degree);
    for (index_t i = 0; i < initialRef; ++i)
        basis.uniformRefine();

    // THB basis wrapping the tensor basis (used even when not refining, for
    // consistency with the H branch).
    gsTHBSplineBasis<2> thb(basis);

    // -------------------------------------------------------------------
    // Options: built-in defaults layered with --options (partial lists stay
    // valid, ORCHESTRATOR RULING 1 pattern).
    // -------------------------------------------------------------------
    gsOptionList opt = defaultDriverOptions();
    // Every key the driver actually knows about, captured BEFORE the file is
    // layered in, so a typo'd key in the XML can be named instead of silently
    // dropped (ignoreIfUnknown) or silently appended (addIfUnknown).
    std::set<std::string> knownKeys;
    for (const gsOptionList::OptionListEntry & e : opt.getAllEntries())
        knownKeys.insert(e.label);
    if (haveFileOpts)
    {
        for (const gsOptionList::OptionListEntry & e : fileOpts.getAllEntries())
            if (!knownKeys.count(e.label))
                gsInfo << "WARNING: option '" << e.label << "' from " << file
                       << " (id=4) is not used by l2projection_rh_schedule_example (typo?)\n";
        opt.update(fileOpts, gsOptionList::addIfUnknown);
        gsInfo << "Method options from -f id=4: " << file << "\n";
    }
    if (!options.empty())
    {
        gsFileData<> fdo(options);
        gsOptionList fo;
        if (fdo.hasId(4)) fdo.getId(4, fo); else fdo.getFirst(fo);
        for (const gsOptionList::OptionListEntry & e : fo.getAllEntries())
            if (!knownKeys.count(e.label))
                gsInfo << "WARNING: option '" << e.label << "' from " << fdo.lastPath()
                       << " is not used by l2projection_rh_schedule_example (typo?)\n";
        opt.update(fo, gsOptionList::addIfUnknown);
        gsInfo << "Method options from --options: " << fdo.lastPath() << "\n";
    }
    // CLI overrides, applied LAST so they win over both the built-in default
    // and any XML: --coarsen only ever turns Coarsen ON (there is no
    // --no-coarsen; omit the flag and rely on the XML/default for "off");
    // --target replaces Target outright (any sign), letting a sweep script
    // switch the band on/off and tune it per run without touching the XML.
    // Target is the error the time step finishes at: the H step stops as
    // soon as the error is inside [Target/Band, Target*Band], so --target
    // alone fixes the endpoint and no second threshold has to be kept in
    // step with it.
    if (coarsenFlag) opt.setSwitch("Coarsen", true);
    if (!math::isnan(targetCli))
    {
        opt.setReal("Target", targetCli);
    }
    GISMO_ENSURE(opt.getInt("BarrierMode") == 0 || opt.getInt("BarrierMode") == 1,
                 "BarrierMode must be 0 (sampled) or 1 (coefficient), got "
                 << opt.getInt("BarrierMode"));

    const real_t target = opt.getReal("Target"), band = opt.getReal("Band");
    const bool useBand = target >= 0;
    if (useBand) GISMO_ENSURE(band >= 1, "Band must be >= 1");
    const bool coarsenSwitch = opt.getSwitch("Coarsen");
    const index_t maxRefIt = opt.getInt("MaxRefIt");
    GISMO_ENSURE(maxRefIt > 0, "MaxRefIt must be positive");

    // -------------------------------------------------------------------
    // Sigma (r-adaptivity map). ORCHESTRATOR RULING 2 precedence, highest
    // first:
    //   1. -E/-R given explicitly on the CLI -> always build from the
    //      task-01 factory (the paper's sweep axes; a file that silently
    //      overrode them would make a sweep produce identical runs).
    //   2. Otherwise, if the file carries a 'composition' that is actually
    //      usable (degree>=2 or more than one element) -> pin sigma from it.
    //   3. Otherwise -> factory from the -E/-R defaults.
    // The shipped default file's composition is a degree-1 single-element
    // identity square (4 corner controls): pinning to it would leave the R
    // and D steps with ZERO design freedom in the headline surface run --
    // exactly the trap this precedence prevents.
    //
    // Sigma's mesh is a REFINEMENT LEVEL of the level-0 analysis basis (not
    // a free knot count): makeIntegrationBasis (the paper's super-mesh)
    // unions the sigma and analysis knot lines, and a non-nested pair costs
    // roughly their SUM of elements while a nested pair costs only the finer
    // of the two. Every analysis level is dyadic, so tying sigma to the same
    // ladder makes a non-nested sigma inexpressible; sigma's DEGREE
    // (--sigmaDeg) stays independent of this.
    // -------------------------------------------------------------------
    std::unique_ptr<gsSquareDomain<real_t> > sigmaPtr;
    bool pinned = false;
    if (!sigmaGivenOnCLI && composition.nPatches() > 0)
    {
        const gsGeometry<real_t>& sgeo = composition.patch(0);
        const gsTensorBSplineBasis<2>* sb =
            dynamic_cast<const gsTensorBSplineBasis<2>*>(&sgeo.basis());
        GISMO_ENSURE(sb, "the 'composition' entry must carry a tensor-product "
                         "B-spline basis");
        const bool pinnable = (sb->degree(0) >= 2 || sb->degree(1) >= 2 ||
                               sb->numElements() > 1);
        if (pinnable)
        {
            sigmaPtr.reset(new gsSquareDomain<real_t>(sgeo));
            sigmaPtr->options().addSwitch("Slide",
                "Boundary controls slide along the boundary", opt.getSwitch("Slide"));
            sigmaPtr->applyOptions();
            pinned = true;
            gsInfo << "sigma source: pinned from the file's 'composition' (degree "
                   << sb->degree(0) << "," << sb->degree(1) << ", "
                   << sb->numElements() << " elements) -- overrides -E/-R\n";
        }
        else
            gsInfo << "composition in file is a trivial bilinear identity; "
                      "building sigma from -E/-R instead\n";
    }
    if (!pinned)
    {
        gsKnotVector<real_t> ks(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
        gsTensorBSplineBasis<2,real_t> b(ks, ks);
        sigmaPtr.reset(new gsSquareDomain<real_t>(b, opt.getSwitch("Slide")));
        gsInfo << "sigma source: factory ("
               << (sigmaGivenOnCLI ? "-E/-R given explicitly on the CLI"
                                    : "-E/-R defaults; no usable file composition")
               << ") -- degree " << sigmaDeg << ", " << ((index_t(1)<<sigmaRef)-1)
               << " interior knots per direction\n";
    }
    gsSquareDomain<real_t>& sigma = *sigmaPtr;
    const gsTensorBSplineBasis<2> * sb =
        dynamic_cast<const gsTensorBSplineBasis<2>*>(&sigma.domain().basis());
    GISMO_ENSURE(sb, "sigma must carry a tensor-product B-spline basis");
    const gsTensorBSplineBasis<2> & sbasis = *sb;

    // -------------------------------------------------------------------
    // Ring monitor parameters (D19: the ring geometry knobs live in the
    // options XML, not the CLI).
    // -------------------------------------------------------------------
    const real_t ringR    = opt.getReal("RingRadius");
    const real_t ringW    = opt.getReal("RingWidth");
    const real_t ringCx0  = opt.getReal("RingCx0");
    const index_t ringAxis = opt.getInt("RingAxis");
    GISMO_ENSURE(ringAxis >= 0 && ringAxis < dim,
        "RingAxis must satisfy 0 <= RingAxis < dim (dim=" << dim << ")");

    // c0 = the physical midpoint of the patch, so the construction reduces
    // EXACTLY to the historical planar case (c0=(0,0) on [-0.5,0.5]^2) and,
    // on the sphere default, is the north pole (0,0,1).
    gsMatrix<real_t> mid(2,1); mid << 0.5, 0.5;
    gsMatrix<real_t> c0m;
    physical.patch(0).eval_into(mid, c0m);
    const gsVector<real_t> c0 = c0m.col(0);

    gsInfo << "schedule=" << sched << " nSteps=" << nSteps << " iterations=" << iterations
           << " project=" << project << " useH=" << useH << "\n";
    gsInfo << "Quadrature: quA=" << opt.getReal("quA") << " quB=" << opt.getInt("quB")
           << " SameElement=" << opt.getSwitch("SameElement") << "\n";
    if (opt.getSwitch("SameElement"))
        gsWarn << "WARNING: SameElement=1 is set, but this driver always integrates "
                  "on a composed map: quadrature points will cross element boundaries "
                  "and results will be WRONG. Set SameElement=0 in --options.\n";

    gsInfo << opt << "\n";
    gsFileData<> fdOut;
    fdOut << opt;
    fdOut.dump(output + "options");

    // -------------------------------------------------------------------
    // Persistent evaluator (shared with the ParaView collection so that
    // addField uses the same expression tree the collection was built
    // with).
    // -------------------------------------------------------------------
    gsExprEvaluator<> pev;

    std::unique_ptr<gsParaviewCollection> solcolC, solcolP;
    if (plot)
    {
        const std::string cname = project ? (output+"solution_C") : (output+"solution");
        solcolC.reset(new gsParaviewCollection(cname, &pev));
        solcolC->options().setSwitch("plotElements", true);
        solcolC->options().setInt("plotElements.resolution", 4);
        solcolC->options().setInt("numPoints", 1000);
        solcolC->options().setInt("precision", 12);
        if (project)
        {
            solcolP.reset(new gsParaviewCollection(output+"solution_P", &pev));
            solcolP->options().setSwitch("plotElements", true);
            solcolP->options().setInt("plotElements.resolution", 4);
            solcolP->options().setInt("numPoints", 1000);
            solcolP->options().setInt("precision", 12);
        }
    }

    std::ofstream csv(output + "convergence.csv");
    // Unbuffered: a sweep kills a run that overruns its time budget, and a
    // buffered stream loses every row it has not flushed -- a 2 h run then
    // yields a 0-byte file instead of the convergence history it did earn.
    csv << std::unitbuf;
    csv << "cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time\n";
    csv << std::scientific << std::setprecision(8);

    gsMultiPatch<> lastSol;
    bool haveProj = false;
    real_t lastErr = 0, lastMinErr = 0, lastMaxErr = 0, lastPctBelowTol = 0;
    // The per-element L2 error vector driving the H-step marker: follows
    // WHICHEVER path last set lastSol/lastErr -- the COMPOSED (C) path's
    // elErr after a plain S, the PROJECTED (P) path's after an 'S,P' row
    // under --project (same convention as the poisson driver's lastElErr:
    // R/D/H see whatever the driver just computed, not always C).
    std::vector<real_t> lastElErr;
    index_t projStep = 0;

    // task 20: ONE update path for all carried state. Every step that
    // changes lastSol/lastErr/lastElErr -- an S row (C or, under --project,
    // P), the H band's inner re-projection, and D's post-relocation
    // re-measurement of the frozen solution -- calls this exactly once, so
    // the invariant "after step X, the carried state describes the state
    // AFTER X" is enforced by construction instead of by convention (task
    // 19's review: a scattered per-branch update left the H band decision
    // reading a stale pre-D error, and a D-after-H measuring a stale
    // pre-refinement solution).
    auto updateState = [&](const ProjResult& r)
    {
        lastSol = r.sol;
        lastErr = r.err;
        lastElErr = r.elErr;
        // r.elErr is empty only if the analysis basis has zero elements,
        // which never happens in practice; guarded anyway (task 20 review
        // notes) rather than dereferencing *std::min/max_element on an empty
        // range.
        if (!r.elErr.empty())
        {
            lastMinErr = *std::min_element(r.elErr.begin(), r.elErr.end());
            lastMaxErr = *std::max_element(r.elErr.begin(), r.elErr.end());
        }
        else
        {
            lastMinErr = lastMaxErr = 0.0;
        }
        // ORCHESTRATOR ADDENDUM (task 07's review): pctBelowTol is the
        // percentage of ELEMENTS below the EQUIDISTRIBUTION threshold
        // Target/sqrt(nElements), not the raw global Target.
        index_t nBelow = 0;
        if (useBand)
        {
            const real_t thrEq = target / math::sqrt((real_t)r.elErr.size());
            for (real_t e : r.elErr) if (e <= thrEq) ++nBelow;
        }
        lastPctBelowTol = useBand ? (100.0*nBelow/(real_t)r.elErr.size()) : 0.0;
        haveProj = true;
    };

    // D13: before the first S step the four error columns are emitted as
    // empty fields, not a literal 0.0 (indistinguishable from a perfect
    // projection to any downstream harvester).
    auto errFields = [&]() -> std::string
    {
        if (!haveProj) return ",,,";
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(8)
            << lastMinErr << "," << lastMaxErr << "," << lastErr << ","
            << lastPctBelowTol;
        return oss.str();
    };

    for (index_t ts = 0; ts < nSteps; ++ts)
    {
        const real_t t = (nSteps > 1) ? static_cast<real_t>(ts) / (nSteps - 1) : 0.0;

        // Ring centres at time t (D19.3): converge from c0 +- RingCx0 along
        // RingAxis at t=0 to c0 (merged) at t=1.
        gsVector<real_t> c1 = c0, c2 = c0;
        c1[ringAxis] += ringCx0 * (1.0 - t);
        c2[ringAxis] -= ringCx0 * (1.0 - t);
        const RingMonitor monitor(c1, c2, ringR, ringW);

        gsInfo << "Time step " << ts << " / " << nSteps - 1
               << "  t=" << std::fixed << std::setprecision(3) << t
               << "  centres c1=[" << c1.transpose() << "]  c2=[" << c2.transpose()
               << "]" << std::defaultfloat << "\n";

        // D9, 2026-08-18: "hold" (deep convergence, err <= NoCoarsenBelow)
        // stops working on THIS time step -- not the whole run. The
        // RingMonitor moves between time steps (D19.3), so a mesh converged
        // for time step ts is not necessarily converged for ts+1; only the
        // rep/op loops below break early, the outer ts loop always runs to
        // nSteps. See the D9 header comment for the "band" vs "hold"
        // distinction this flag is built on.
        bool stepConverged = false;

        for (index_t rep = 0; rep < iterations && !stepConverged; ++rep)
        {
            const index_t cycle = ts * iterations + rep;
            gsInfo << "  Cycle " << cycle << " (t-step " << ts << ", rep " << rep
                   << ") [" << sched << "]\n";

            for (char op : sched)
            {
                if (stepConverged) break;
                const index_t dofs = useH ? thb.size() : basis.size();
                gsStopwatch stepTimer;

                if (op == 'S')
                {
                    const gsTensorBSplineBasis<2>& tb = useH ? thb.tensorLevel(thb.maxLevel()) : basis;
                    const gsBasis<>& active = useH ? static_cast<const gsBasis<>&>(thb)
                                                    : static_cast<const gsBasis<>&>(basis);
                    gsMultiBasis<> mb(active);
                    gsTensorBSplineBasis<2> ibasis =
                        gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                        makeIntegrationBasis(tb, sbasis);
                    gsMultiBasis<> ib(ibasis);
                    gsComposedGeometry<> composed(sigma, physical.patch(0));
                    gsMultiPatch<> cmp; cmp.addPatch(composed);

                    const index_t step = projStep++;
                    ProjResult cres = doProject(cmp, mb, ib, active, monitor, opt, t, step,
                                                 plot ? &pev : nullptr,
                                                 plot ? solcolC.get() : nullptr, "C");
                    // task 20: single update path -- refreshes lastSol,
                    // lastErr, lastElErr, lastMinErr/lastMaxErr/lastPctBelowTol
                    // together, from this C step's own result.
                    updateState(cres);
                    csv << cycle << ",S,C," << active.size() << "," << sigma.nControls() << ","
                        << lastMinErr << "," << lastMaxErr << "," << lastErr << "," << lastPctBelowTol << ","
                        << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";

                    if (project)
                    {
                        // D8: project the COMPOSED geometry S o sigma onto the
                        // analysis basis, then run the function projection
                        // AGAIN on that projected (ordinary, non-composed)
                        // geometry -- no super-mesh is needed once the map is
                        // an ordinary spline, so pmb serves as both
                        // projection and integration basis.
                        gsStopwatch pTimer;
                        gsMatrix<> pcoefs;
                        gsL2Projection<real_t>::project(active, ibasis, composed, pcoefs);
                        pcoefs.resize(active.size(), physical.patch(0).targetDim());
                        gsGeometry<>::uPtr pg = active.makeGeometry(give(pcoefs));
                        gsMultiPatch<> pmp; pmp.addPatch(*pg);
                        gsMultiBasis<> pmb(active);
                        ProjResult pres = doProject(pmp, pmb, pmb, active, monitor, opt, t, step,
                                                     plot ? &pev : nullptr,
                                                     plot ? solcolP.get() : nullptr, "P");
                        // The P path's solution/error are what the following
                        // R/D/H steps see (same choice as the poisson driver);
                        // task 20: routed through the single update path, so
                        // this overwrite of the C-step state above is total
                        // and coherent (lastSol/lastErr/lastElErr/lastMinErr/
                        // lastMaxErr/lastPctBelowTol all become the P step's).
                        updateState(pres);
                        csv << cycle << ",S,P," << active.size() << "," << sigma.nControls() << ","
                            << lastMinErr << "," << lastMaxErr << "," << lastErr << "," << lastPctBelowTol << ","
                            << sigma.minJacobian(7) << "," << pTimer.stop() << "\n";
                    }
                }
                else if (op == 'R')
                {
                    if (!haveProj)
                    {
                        gsInfo << "  R | no previous projection, skipping\n";
                        csv << cycle << ",R,-," << dofs << "," << sigma.nControls() << ","
                            << errFields() << "," << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";
                        continue;
                    }

                    // Generic lambda instantiated at compile time for each
                    // MonitorMode (a non-type template parameter of
                    // gsAdaptiveParametrization).
                    auto relocate = [&](auto mode)
                    {
                        gsStopwatch timer;
                        std::unique_ptr<gsOptimizer<real_t> > optimizer = makeOptimizer(opt);

                        // NOTE: no integration basis is handed to the R step.
                        // The R step's integrand is
                        //   omega(||grad_x f||) * ||J_g J_sigma||^2 / chi,
                        // whose only ingredients are sigma, the COARSE
                        // level-0 geometry physical.patch(0) (no H step ever
                        // touches it), and the monitor -- the ANALYSIS basis
                        // appears nowhere in it. Constructing without an
                        // integration basis defaults it to the geometry's own
                        // (level-0) basis, into which makeIntegrationBasis
                        // then inserts sigma's knots: the result is exactly
                        // sigma's own mesh, degree-raised, whose size is
                        // fixed by --sigmaRef instead of growing with h.
                        //
                        // parametric=false: f (RingMonitor) is defined in
                        // AMBIENT coordinates ([-0.5,0.5]^2 planar, R^3 on
                        // the surface default), not on the parametric square
                        // [0,1]^2 sigma maps INTO. With parametric=true the
                        // monitor is sampled at sigma(xi) in [0,1]^2 instead
                        // of at the physical point G(sigma(xi)): sigma then
                        // adapts to a shifted/wrong domain and drifts into a
                        // corner (bitten once already for the planar case,
                        // and worse on a surface, where [0,1]^2 is not even
                        // in the same space as f's ambient R^3 domain).
                        //
                        // Nothing else changes for a surface: the td>dd
                        // branch of the objective
                        // (gsAdaptiveParametrization.hpp:448-463) already
                        // carries the surface area element
                        // g_S=sqrt(det(J_g^T J_g)), and the GradientBased
                        // weight already contains the tangential projection
                        // via Cg^{-1} -- adding either of those again in this
                        // driver would be a known trap (D19.5) and is
                        // FORBIDDEN.
                        gsAdaptiveParametrization<real_t, mode.value> rel(
                            sigma, physical.patch(0), monitor, *optimizer, false);
                        rel.options().setReal("Smoothing", opt.getReal("Smoothing"));
                        rel.options().setReal("Penalty",   opt.getReal("Penalty"));
                        // quA/quB mirrored from the resolved driver options
                        // (same rule as the element-wise error sweeps: the
                        // assembler/optimizer and the norms that measure it
                        // must use the SAME quadrature).
                        rel.options().setReal("quA", opt.getReal("quA"));
                        rel.options().setInt ("quB", opt.getInt ("quB"));
                        rel.solve();
                        const real_t tR = timer.stop();

                        gsInfo << "  R | mode "
                               << (mode.value==MonitorMode::GradientBased?"gradient":"value")
                               << " | parametric=false (ambient-coordinate monitor, D19.4)"
                               << " | dofs " << std::setw(6) << (useH ? thb.size() : basis.size())
                               << " | sigma controls " << std::setw(4) << sigma.nControls()
                               << " | min det J_sigma " << std::scientific << std::setprecision(3)
                               << rel.computeMinJacobian()
                               << " | min det J_sigma (certificate) " << std::scientific
                               << std::setprecision(3) << sigma.minDetJCoefficient();
                        gsInfo << " | maxIt " << opt.getInt("MaxIterations")
                               << " | " << std::fixed << std::setprecision(2) << tR << " s\n"
                               << std::defaultfloat;
                    };

                    const std::string mm = opt.askString("MonitorMode","gradient");
                    if (mm=="gradient")
                        relocate(std::integral_constant<enum MonitorMode, MonitorMode::GradientBased>());
                    else if (mm=="value")
                        relocate(std::integral_constant<enum MonitorMode, MonitorMode::ValueBased>());
                    else
                        GISMO_ERROR("Unknown MonitorMode '"<<mm<<"' (value | gradient)");

                    // task 20 (review fix 1): the R step relocates sigma, i.e.
                    // it changes the composed geometry S o sigma the frozen
                    // lastSol is measured against -- exactly the mutation that
                    // made the D step's stale-lastErr bug real. Re-measure the
                    // FROZEN lastSol (no re-projection) against the POST-R
                    // sigma, through the same single update path as D and H,
                    // so the next row (typically an H band decision) reads the
                    // post-R error instead of the pre-R one.
                    {
                        const gsBasis<>& activeR = useH ? static_cast<const gsBasis<>&>(thb)
                                                         : static_cast<const gsBasis<>&>(basis);
                        gsMultiBasis<> mbR(activeR);
                        gsComposedGeometry<> composedR(sigma, physical.patch(0));
                        gsMultiPatch<> cmpR; cmpR.addPatch(composedR);
                        ProjResult rres = measureError(cmpR, mbR, activeR, lastSol, monitor, opt);
                        updateState(rres);
                    }

                    csv << cycle << ",R,-," << dofs << "," << sigma.nControls() << ","
                        << errFields() << "," << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";
                }
                else if (op == 'D')
                {
                    if (!haveProj)
                    {
                        gsInfo << "  D | no previous projection, skipping\n";
                    }
                    else
                    {
                        // gsOptL2 now supports the surface case (targetDim()>2)
                        // via the area element g_S = sqrt(det(J_S^T J_S)),
                        // landed by task 04/14 (see gsAdaptiveParametrization.h
                        // :495-535, the @note above the gsOptL2 declaration).
                        // No dimension guard is needed here any more.

                        gsStopwatch timer;
                        const gsTensorBSplineBasis<2>& tb = useH ? thb.tensorLevel(thb.maxLevel()) : basis;
                        gsTensorBSplineBasis<2> ibasis =
                            gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                            makeIntegrationBasis(tb, sbasis);

                        // Fold-barrier params for the D step: BarrierMu/
                        // BarrierEps (task 18) are the D-specific keys --
                        // Penalty (default 1e-3) is the R step's Winslow fold
                        // weight, a DIFFERENT objective, never reused here;
                        // the fitting driver's own D step needs mu>=1e3 to
                        // dominate the data term (mu=1 "provably lets sigma
                        // fold and the step gets rejected", see
                        // rh-adaptive_fitting_example.cpp:354-357).
                        // BarrierMode=0 (the default) selects Gauss quadrature
                        // on sigma's own knot mesh (order dQuB). dQuB stays a
                        // local hardcoded 8 -- task 18 explicitly does NOT
                        // wire it to the driver's own quB key: after task 18
                        // quB=2 (the IGA error quadrature), and 2 barrier-
                        // Gauss points per element per direction lets HLBFGS
                        // fold sigma between the nodes (task 09), where >=10
                        // are required.
                        const real_t  dMu   = opt.getReal("BarrierMu");
                        const real_t  dEps  = opt.getReal("BarrierEps");
                        const gsFoldBarrierMode dMode = opt.getInt("BarrierMode") == 0
                            ? gsFoldBarrierMode::Sampled : gsFoldBarrierMode::Coefficient;
                        const index_t dQuB  = 8;

                        std::unique_ptr<gsOptimizer<real_t> > optimizer = makeOptimizer(opt);
                        gsOptL2<real_t> obj(sigma, physical.patch(0), lastSol.patch(0), monitor,
                                            &ibasis, dMu, dEps, dMode, dQuB, /*parametric=*/false);
                        // Element-sweep quA/quB (distinct from the barrier's
                        // own dMode/dQuB above), mirrored from the resolved
                        // driver options -- same rule as the R step and the
                        // error evaluators.
                        obj.options().setReal("quA", opt.getReal("quA"));
                        obj.options().setInt ("quB", opt.getInt ("quB"));
                        optimizer->setProblem(&obj);

                        const gsVector<real_t> uPrev = sigma.getControls();
                        optimizer->solve(uPrev);
                        gsVector<real_t> uNew = optimizer->currentDesign();
                        sigma.setControls(uNew);
                        real_t minDetJ = sigma.minJacobian(7);
                        if (!uNew.allFinite() || minDetJ <= 0)
                        {
                            gsWarn << "D step rejected (min det J_sigma = " << minDetJ
                                   << "); sigma restored.\n";
                            sigma.setControls(uPrev);
                            minDetJ = sigma.minJacobian(7);
                        }
                        const real_t tD = timer.stop();
                        gsInfo << "  D | dofs " << std::setw(6) << (useH ? thb.size() : basis.size())
                               << " | sigma controls " << std::setw(4) << sigma.nControls()
                               << " | min det J_sigma " << std::scientific << std::setprecision(3)
                               << minDetJ
                               << " | maxIt " << opt.getInt("MaxIterations")
                               << " | " << std::fixed << std::setprecision(2) << tD << " s\n"
                               << std::defaultfloat;

                        // task 19/20: re-measure the FROZEN lastSol (no
                        // re-projection -- that would destroy the
                        // variable-projection structure the D step depends
                        // on, plan D5) against the POST-step sigma (either
                        // relocated, or restored to uPrev above on rejection,
                        // in which case this reproduces the pre-D state
                        // exactly). Mirrors rh-adaptive_fitting_example.cpp's
                        // D step (updateParams() then computeErrors()): here
                        // sigma IS the updated parametrization already, so
                        // only the error evaluation needs repeating, via the
                        // same D2 element-wise sweep (measureError, shared
                        // with doProject) as every other L2 error in this
                        // file -- not just the global integral, so
                        // minErr/maxErr/pctBelowTol come out post-D too.
                        // Pushed through the single update path
                        // (updateState): this is what makes the NEXT step
                        // (typically an H band decision, D9) see the POST-D
                        // error instead of a stale pre-D one -- the
                        // behavioural bug task 19's review found.
                        const gsBasis<>& activeD = useH ? static_cast<const gsBasis<>&>(thb)
                                                         : static_cast<const gsBasis<>&>(basis);
                        gsMultiBasis<> mbD(activeD);
                        gsComposedGeometry<> composedD(sigma, physical.patch(0));
                        gsMultiPatch<> cmpD; cmpD.addPatch(composedD);
                        ProjResult dres = measureError(cmpD, mbD, activeD, lastSol, monitor, opt);
                        updateState(dres);
                    }
                    csv << cycle << ",D,-," << dofs << "," << sigma.nControls() << ","
                        << errFields() << "," << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";
                }
                else if (op == 'U')
                {
                    if (useH) thb.uniformRefine();
                    else      basis.uniformRefine();

                    // task 20 (review fix 3): uniformRefine() changes the
                    // analysis basis/mesh the carried lastSol/lastElErr were
                    // measured on -- the U row must not emit post-U dofs
                    // beside pre-U error columns (and a stale lastElErr trips
                    // the H-step "mesh changed since the last projection"
                    // guard). Re-measure the FROZEN lastSol (no
                    // re-projection) against the new mesh, same D-step
                    // pattern as R and the H terminal refresh.
                    if (haveProj)
                    {
                        const gsBasis<>& activeU = useH ? static_cast<const gsBasis<>&>(thb)
                                                         : static_cast<const gsBasis<>&>(basis);
                        gsMultiBasis<> mbU(activeU);
                        gsComposedGeometry<> composedU(sigma, physical.patch(0));
                        gsMultiPatch<> cmpU; cmpU.addPatch(composedU);
                        ProjResult ures = measureError(cmpU, mbU, activeU, lastSol, monitor, opt);
                        updateState(ures);
                    }

                    gsInfo << "  U | dofs " << std::setw(6) << (useH ? thb.size() : basis.size())
                           << " | " << std::fixed << std::setprecision(2) << stepTimer.stop() << " s\n"
                           << std::defaultfloat;
                    csv << cycle << ",U,-," << (useH ? thb.size() : basis.size()) << ","
                        << sigma.nControls() << "," << errFields() << ","
                        << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";
                }
                else if (op == 'H')
                {
                    if (!(useH && haveProj))
                    {
                        // ORCHESTRATOR ADDENDUM (task 07's review): a
                        // schedule step with an unmet precondition SKIPS
                        // with a warning; it never aborts -- uniform at
                        // every position in the schedule.
                        gsInfo << "  H | " << (useH ? "no previous projection"
                                                    : "'H' not in --schedule")
                               << ", skipping\n";
                    }
                    else
                    {
                        bool inBandStop = false;
                        // D9, 2026-08-18: set ONLY on the true "hold" branch
                        // (err <= NoCoarsenBelow), never on "band" (the dead
                        // zone) or "below"-without-coarsening -- see the
                        // holdReached computation below for why a plain
                        // "!doRef && !doCrs" check is NOT sufficient (it is
                        // also true for "below" when Coarsen is off).
                        bool letterHoldReached = false;
                        // task 20 (review fix 2): tracks whether ANY iteration
                        // of this H letter actually changed the mesh, so the
                        // terminal refresh below (which re-measures the
                        // frozen lastSol on whatever mesh the loop ends on)
                        // only fires when it can matter.
                        bool anyChanged = false;
                        for (index_t hit = 0; hit < maxRefIt; ++hit)
                        {
                            const gsBasis<>& active = static_cast<const gsBasis<>&>(thb);
                            // lastElErr is cached from the last S step's
                            // elements; if the mesh changed since then with
                            // no intervening S (e.g. this is a re-entry
                            // after a bare U), it no longer lines up.
                            if (lastElErr.size() != (size_t)active.numElements())
                            {
                                gsInfo << "  H | mesh changed since the last projection, "
                                          "skipping (insert an S before this H)\n";
                                break;
                            }

                            // D9 band decision on the GLOBAL L2 error of the
                            // last projection. 2026-08-20: the band is a
                            // TARGET ZONE the H step drives the error INTO,
                            // from whichever side it starts on, and reaching
                            // it ends the current time step:
                            //
                            //   err >  Target*Band  "above" refine only
                            //   in the band         "band"  ONE combined
                            //                               refine+coarsen
                            //                               sweep, then stop
                            //   err <= Target/Band  "below" coarsen only
                            //
                            // The previous rule refined inside the band too
                            // and only stopped at Target/Band^2, so a run had
                            // to overshoot its requested Target by Band^2
                            // before it could finish -- and on a case whose
                            // error floor lies inside the band it never did.
                            // Coarsening follows the branch now rather than
                            // running on every H step.
                            bool doRef=true, doCrs=coarsenSwitch; bool inBand=false;
                            const char* branch = "-";
                            if (useBand)
                            {
                                if      (lastErr >  target*band)
                                { doRef=true;  doCrs=false;         branch="above"; }
                                else if (lastErr >  target/band)
                                { doRef=true;  doCrs=coarsenSwitch; branch="band"; inBand=true; }
                                else
                                { doRef=false; doCrs=coarsenSwitch; branch="below";
                                  // Below the band with coarsening off, no
                                  // action can raise the error back into it,
                                  // so the loop would spin. Stop instead.
                                  if (!coarsenSwitch) inBand=true; }
                            }
                            // Reaching the band is the stop signal: the run
                            // has arrived at the error it was asked for.
                            const bool holdReached = inBand;
                            if (!doRef && !doCrs)
                            {
                                gsInfo << "  H | it " << (hit+1) << "/" << maxRefIt
                                       << " | err " << std::scientific << std::setprecision(3)
                                       << lastErr << " (" << branch << ") | no action"
                                       << " | elements " << thb.numElements() << "\n"
                                       << std::defaultfloat;
                                inBandStop = true;
                                letterHoldReached = holdReached;
                                break;
                            }

                            typedef gsHElementMarker<2,real_t> marker_t;
                            marker_t marker(thb);
                            marker.options().setInt   ("RefineRule",   opt.getInt ("RefineRule"));
                            marker.options().setReal  ("RefineParam",  opt.getReal("RefineParam"));
                            marker.options().setInt   ("CoarsenRule",  opt.getInt ("CoarsenRule"));
                            marker.options().setReal  ("CoarsenParam", opt.getReal("CoarsenParam"));
                            marker.options().setInt   ("MaxLevel",     opt.getInt ("MaxLevel"));
                            marker.options().setSwitch("Admissible",   opt.getSwitch("Admissible"));
                            marker.options().setSwitch("Extension",    opt.getSwitch("Extension"));
                            marker.setErrors(lastElErr);

                            const marker_t::HElementContainer markedRef =
                                doRef ? marker.markRef() : marker_t::HElementContainer();
                            std::vector<index_t> boxes = doRef ? marker.toRefBoxes(markedRef)
                                                                : std::vector<index_t>();
                            if (doRef && markedRef.empty())
                            {
                                const index_t nEl = active.numElements();
                                const index_t rule = opt.getInt("RefineRule");
                                const real_t param = opt.getReal("RefineParam");
                                gsInfo << "  H | WARNING: RefineRule " << rule << " marked 0 of " << nEl
                                       << " elements (RefineParam = " << param << ")";
                                if (rule == 2) gsInfo << "; PUCA marks floor(RefineParam*nElements), which is "
                                                         "0 below " << (param>0 ? 1.0/param : 0.0) << " elements";
                                gsInfo << "\n";
                            }
                            std::vector<index_t> crsBoxes; size_t nCrs = 0;
                            if (doCrs)
                            {
                                const marker_t::HElementContainer markedCrs = marker.markCrs(markedRef);
                                nCrs = markedCrs.size();
                                crsBoxes = marker.toCrsBoxes(markedCrs);
                            }

                            // Both box lists refer to the SAME (pre-update)
                            // mesh; the refined and coarsened regions are
                            // disjoint by construction.
                            const bool changed = !boxes.empty() || !crsBoxes.empty();
                            anyChanged = anyChanged || changed;
                            if (!boxes.empty())    thb.refineElements(boxes);
                            if (!crsBoxes.empty()) thb.unrefineElements(crsBoxes);

                            gsInfo << "  H | it " << (hit+1) << "/" << maxRefIt
                                   << " | err " << std::scientific << std::setprecision(3)
                                   << lastErr << " (" << branch << ")"
                                   << " | refined " << std::setw(4) << markedRef.size();
                            if (doCrs) gsInfo << " | coarsened " << std::setw(4) << nCrs;
                            gsInfo << " | elements " << thb.numElements() << std::defaultfloat << "\n";

                            // No boxes marked at all -> no further progress
                            // possible, stop instead of spinning to maxRefIt.
                            if (!changed) break;
                            // Error is inside the band: done for THIS letter
                            // (inner loop only) -- "hold" never reaches this
                            // line (it always breaks above, doRef=doCrs=false),
                            // so holdReached is always false here; kept for
                            // symmetry with the other break point above.
                            if (inBand) { inBandStop = true; letterHoldReached = holdReached; break; }

                            // Re-project on the new mesh so the next
                            // iteration's band decision sees the CURRENT
                            // state; skipped on the last iteration (the
                            // schedule's next S does it anyway). This is an
                            // INNER re-projection: no CSV row, no Paraview
                            // frame.
                            if (hit + 1 < maxRefIt)
                            {
                                const gsTensorBSplineBasis<2>& tb2 = thb.tensorLevel(thb.maxLevel());
                                gsMultiBasis<> mb2(active);
                                gsTensorBSplineBasis<2> ibasis2 =
                                    gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                                    makeIntegrationBasis(tb2, sbasis);
                                gsMultiBasis<> ib2(ibasis2);
                                gsComposedGeometry<> composed2(sigma, physical.patch(0));
                                gsMultiPatch<> cmp2; cmp2.addPatch(composed2);
                                ProjResult r2 = doProject(cmp2, mb2, ib2, active, monitor, opt,
                                                          t, -1, nullptr, nullptr, "S*");
                                // task 20: single update path, including
                                // lastSol -- the old code refreshed only
                                // lastErr/lastElErr here and left lastSol at
                                // the PRE-refinement solution, so a D step
                                // run immediately after this H measured the
                                // wrong solution on the new mesh (task 19's
                                // review, ~3.6% shift).
                                updateState(r2);
                            }
                        }
                        (void)inBandStop;
                        // Propagate genuine convergence out to the op/rep
                        // loops (see the `stepConverged` declaration at the
                        // top of the ts loop for why this does NOT also
                        // break ts: the moving RingMonitor means convergence
                        // at time step ts says nothing about ts+1).
                        if (letterHoldReached) stepConverged = true;

                        // task 20 (review fix 2): every exit path above --
                        // normal exhaustion, the `!changed` break, and the
                        // `inBand` break -- can leave lastSol/lastErr/lastElErr
                        // describing a PRE-refinement mesh (the inner
                        // re-projection only runs strictly between iterations,
                        // never on the terminal one, and MaxRefIt defaults to
                        // 1, so in the default configuration the H branch
                        // never called updateState at all). Chosen fix (see
                        // task 20 review, required fix 2): the cheap,
                        // spec-conforming D-step pattern -- re-measure the
                        // FROZEN lastSol (no re-projection) against the FINAL
                        // mesh once, rather than promoting the inner
                        // re-projection to run unconditionally (which would
                        // cost an extra full projection per H letter).
                        if (anyChanged)
                        {
                            const gsBasis<>& activeHend = static_cast<const gsBasis<>&>(thb);
                            gsMultiBasis<> mbHend(activeHend);
                            gsComposedGeometry<> composedHend(sigma, physical.patch(0));
                            gsMultiPatch<> cmpHend; cmpHend.addPatch(composedHend);
                            ProjResult hres = measureError(cmpHend, mbHend, activeHend, lastSol, monitor, opt);
                            updateState(hres);
                        }
                    }
                    csv << cycle << ",H,-," << thb.size() << "," << sigma.nControls() << ","
                        << errFields() << "," << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";
                }
            } // end schedule loop
        } // end rep loop
    } // end time-step loop

    // Sigma sanity summary (D19): reports whether sigma has collapsed to a
    // corner of [0,1]^2 -- the parametric=true failure mode of D19.4 -- so a
    // run's log alone is enough to check it, without --plot.
    {
        const gsVector<real_t> ctrl = sigma.getControls();
        gsInfo << "Final sigma controls: min=" << ctrl.minCoeff()
               << " max=" << ctrl.maxCoeff() << " (n=" << ctrl.size() << ")\n";
    }

    if (plot)
    {
        solcolC->save();
        if (project) solcolP->save();
        // The final composed geometry G = S o sigma, freshly built from the
        // final sigma/analysis state.
        gsComposedGeometry<> finalComposed(sigma, physical.patch(0));
        gsWriteParaview(finalComposed, output+"geometry", 1000, true, true);
        gsWriteParaview(sigma.domain(), output+"sigma", 1000, true, true);
        const gsBasis<>& active = useH ? static_cast<const gsBasis<>&>(thb)
                                        : static_cast<const gsBasis<>&>(basis);
        gsWriteParaview(active, output+"mesh", 1000, true);
    }
    return 0;
}
