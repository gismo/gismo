/** @file rh-adaptive_lrfitting_example.cpp

    @brief Schedule-driven adaptive fitting of a parametrized point cloud
    with a composed LR-spline geometry G = S o sigma (case 4 of the
    hr-adaptivity study).

    One of four "rh-driver-unification" drivers (poisson, l2projection,
    fitting, lrfitting) that share a 10-flag CLI core and a schedule
    alphabet, but implement only the letters their problem supports (D4):

        letter  meaning                                poisson  l2  fitting  lrfitting
        S       solve / project / fit (the primal step)   x      x     x        x
        R       monitor-driven relocation of sigma        x      x     -        -
        D       error-driven relocation of sigma          -      x     x        x
        H       local refinement (THB / LR)                x      x     x        x
        U       uniform refinement                         x      x     x        -

    This driver is the LR-spline twin of rh-adaptive_fitting_example.cpp
    (THB): it implements ONLY S (fit), D (error-driven relocation of sigma)
    and H (LR refinement). It does NOT implement R (for point-cloud fitting
    the r-adaptivity IS the D step -- there is no separate monitor-driven
    relocation) and does NOT implement U: gsLRSplineBasis::uniformRefine is
    GISMO_NO_IMPLEMENTATION (optional/gsLRSplines/src/gsLRSplineBasis.hpp:1497),
    so there is no library call to route -U through; use -r to set the
    initial refinement level instead. There are no F/P aliases either (the
    shared alphabet has none; the old "F" of this driver is now "S"). Every
    letter outside S, D, H is rejected with a message naming this driver and
    the reason.

    CLI (11 flags, D6-amendment: --options replaces the original spec's -f,
    same rationale as the THB driver -- this driver has no separate PROBLEM
    file beyond the point cloud):

        -d, --data      Input point-cloud XML        (default "fitting/ogre.xml")
        -o, --output    Output directory              (default "rh_lrfitting_output")
            --plot      Write ParaView files
        -p, --degree    Degree of the LR basis, both directions (default 2)
        -r, --refine    Initial refinement LEVEL of the S basis (default 2):
                        nKnots = (1<<r)-1 interior knots in EACH direction of
                        the starting tensor knot vector, applied AT
                        CONSTRUCTION (uniformRefine throws on an LR basis).
        -E, --sigmaDeg  Degree of the sigma map        (default 2)
        -R, --sigmaRef  Refinement LEVEL of sigma's mesh (default 3)
            --schedule  Cycle string over {S,D,H}       (default "SH")
        -i, --iter      Maximum number of schedule cycles (default 3)
            --aniso     H-step refinement: 0 = isotropic (baseline, default),
                        1 = per-element anisotropic (element-gradient
                        criterion; see "ANISOTROPIC REFINEMENT" below). The
                        one driver-specific flag: it selects the arm a
                        figure compares.
            --options   Method option list XML (layered over the driver
                        defaults below via gsOptionList::update with
                        ignoreIfUnknown)

    --project is NOT offered here (D8): the projected path L2-projects
    S o sigma onto an ordinary spline space over the super-mesh built by
    gsAdaptiveParametrization<...>::makeIntegrationBasis(tb, sbasis)
    (src/gsAssembler/gsAdaptiveParametrization.h:435), which takes TWO
    gsTensorBSplineBasis. An LR basis has no tensor level, so the super-mesh
    -- and therefore the projected path -- does not exist for this driver.
    The CSV `path` column is therefore always 'C' on S rows here (D/H rows
    carry '-'; see the D13 block below). A user copying a
    command line from the THB driver gets a stub switch that explains this
    (see main(), below the schedule validation) instead of "unknown option".

    Driver keys carried on the resolved gsOptionList (D6; declared here even
    when NOT honoured, so a single reference XML serves every driver --
    "honoured" column states this driver's behaviour). This is one of four
    rh drivers sharing ONE unified 27-key option table (rh-driver-
    unification task 18); the D-step fold-barrier weight, formerly named
    "Penalty" here, is now "BarrierMu" -- "Penalty" is reserved, across all
    four drivers, for the R-step Winslow fold penalty (which this driver
    does not implement, hence IGNORED below):

        RefineRule      3        gsHElementMarker-style refinement rule: 1=GARU,
                                  2=PUCA, 3=BULK (D7, inline reimplementation --
                                  honoured). Unified on BULK (task 18, defect
                                  1): with RefineParam=0.1, PUCA marked
                                  floor(0.1*nElements) = 0 elements below 10
                                  elements, refining nothing in early cycles.
        RefineParam     0.5      refinement rule parameter (honoured; unified
                                  value, task 18)
        CoarsenRule     3        coarsening rule -- IGNORED, no LR coarsening (D9)
        CoarsenParam    0.1      coarsening rule parameter -- IGNORED, idem
        Coarsen         false    enable coarsening in the H step -- IGNORED, idem
        MaxLevel        10       marker level cap -- IGNORED, LR elements carry
                                  no level
        Admissible      true     admissible closure -- IGNORED, no hierarchical
                                  admissibility on LR
        Extension       true     marker box extension -- IGNORED, LR marks
                                  without a halo (documented asymmetry vs THB)
        Smoothing       1.0      R-step Winslow theta -- IGNORED, no R step here
        Penalty         1e-3     IGNORED: R-step Winslow fold penalty, no R
                                  step in this driver (kept for key-set parity)
        BarrierMu       1000.0   D-step fold-barrier weight mu (honoured; old
                                  --dirMu / old "Penalty", renamed task 18)
        BarrierEps      5e-2     D-step fold barrier: det J_sigma floor
                                  (honoured; promoted from the named constant
                                  barrierEps, same value -- task 18)
        Target          1e-2     error band: absolute rmse target (honoured)
        Band            2.0      error band: dead-zone factor >= 1 (honoured)
        NoCoarsenBelow  0.0      error band: below this rmse, hold (honoured,
                                  only to reach the `hold` branch)
        MaxRefIt        1        inner adapt iterations per H letter (honoured)
        Lambda          1e-6     gsFitting smoothing weight (honoured)
        Optimizer       "HLBFGS" D-step optimizer backend: gsOptim | HLBFGS
                                  (honoured)
        MaxIterations   100      optimizer iteration cap -- honoured; see the
                                  D-step comment: LARGE values are
                                  COUNTER-PRODUCTIVE here
        OptTol          1e-6     optimizer gradient tolerance (unified name;
                                  mapped onto gsOptim's GradErrTol or gsHLBFGS's
                                  MinGradLen -- NOT gsHLBFGS's own "GradTol",
                                  the line-search curvature constant) (honoured)
        Verbose         0        optimizer/driver verbosity (honoured)
        MonitorMode     "value"  R-step monitor mode -- IGNORED, no R step
                                  (spelling unified with l2's vocabulary,
                                  task 18; was "ValueBased")
        Slide           true     sigma boundary controls may slide (honoured)
        quA             1.0      IGNORED -- no projection path in this driver
                                  (--project rejected, D8); the D-step barrier
                                  hard-codes its own quA=1 inside
                                  gsFoldBarrier's Sampled mode, so this driver
                                  could never honour a quA of its own
        quB             8        quadrature extra nodes, D-step barrier
                                  (honoured; old --dirQuB, 8 is the tuned
                                  value). The D-step barrier quadrature is a
                                  Gauss rule on sigma's own knot mesh, a
                                  DIFFERENT quantity from the IGA drivers'
                                  quB (the error/assembly quadrature).
        SameElement     false    IGNORED -- no gsExprAssembler path in this
                                  driver (LR fitting uses gsLRFitting's own
                                  BiCGSTABILUT solve, not an assembled system)

    Driver-specific option entries (beyond D6's shared list), carrying the
    surviving aniso* tuning knobs from the old CLI:

        AnisoTol        0.1      direction tie band (relative unless AnisoAbsTol)
        AnisoComp       -2       target component driving the indicator (-2:
                                  last row = height, -1: norm over all, >=0:
                                  that component)
        AnisoAbs        false    integrate |dS/dxi| instead of |int dS/dxi|
        AnisoAbsTol     false    paper-faithful unnormalised ABSOLUTE tie test
        AnisoQuB        0        extra Gauss nodes per direction for the
                                  direction indicator
        AnisoDump       false    dump a per-element direction record to
                                  <output>direction_dump.csv
        DirSkip         0        skip D steps in the first DirSkip cycles
        CheckGradient   false    finite-difference check of the D-step gradient
                                  (first D step only)
        nSamples        1000     ParaView sampling resolution

    The D-step barrier floor eps is the option "BarrierEps" (task 18;
    formerly a named constant, barrierEps, in this file -- promoted to an
    option, same value 5e-2, so this driver's key set matches the other
    three): it is a tuned floor -- below it the barrier stops preventing
    folds, above it it distorts sigma. The legacy uniform barrier grid is
    gone; the barrier's representation is now selected by "BarrierMode"
    (0 = sampled Gauss rule on sigma's own knot mesh, the shipped default;
    1 = Bezier-coefficient barrier, no quadrature).

    D7 -- MARKING (RELAXED FOR THIS DRIVER, on user instruction: the LR
    driver is experimental). gsHElementMarker itself cannot be used here:
    its gsBasis constructor does
    dynamic_cast<const gsHTensorBasis<d,T>&>(basis) (gsHElementMarker.hpp:31),
    which throws std::bad_cast on an LR basis, and its element type
    gsHElement<d,T> carries a hierarchical level() LR elements do not have.
    The three rules (GARU/PUCA/BULK) are therefore reproduced INLINE on a
    std::vector<std::pair<index_t,real_t>> of (element index, error), on a
    best-effort basis; no equivalence test against gsHElementMarker is
    required or attempted, and a numerical mismatch against the THB marker is
    NOT a driver defect. The MaxLevel/Admissible guards inside the library
    rules have no LR analogue and are dropped (see the ignored table above).
    Per-element error for a point cloud: the MAXIMUM pointwise error over the
    data points falling in that element (elements with no data point get
    error 0) -- the direct analogue of gsHFitting's "refine the cells
    containing an over-tolerance point", independent of how many points land
    in a cell.

    D9 -- error-band rule for the H step, driven by the rmse of the
    last S step (mirrors the THB driver's band):

        rmse >  Target*Band                  -> refine                    branch "above"
        Target/Band < rmse <= Target*Band     -> refine (+coarsen)         branch "band"
        NoCoarsenBelow < rmse <= Target/Band  -> coarsen only              branch "below"
        rmse <= NoCoarsenBelow                -> hold (no action), STOP   branch "hold"

    2026-08-18 CHANGE: the dead zone ("band") no longer stops the run --
    same fix as the other three drivers, same reasoning: see the poisson
    driver's D9 comment. Only "hold" (deep convergence) stops the outer
    schedule; NoCoarsenBelow must sit meaningfully below Target/Band for
    "hold" to ever be reached (--target derives it as Target/Band^2, see
    its own comment above).

    -i/--iter remains a hard cap regardless of the band. COARSENING IS
    UNAVAILABLE: gsLRSplineBasis has no unrefine (unrefine / unrefineElements
    are GISMO_NO_IMPLEMENTATION, optional/gsLRSplines/src/gsLRSplineBasis.hpp:
    1263, 1485). Whenever the "below" branch (or the coarsen half of "band")
    is reached, this driver prints an explicit message (see hStep()) and
    takes NO basis action -- silently doing nothing is explicitly not
    acceptable. Max pointwise error is a diagnostic in the log and the
    `maxErr` CSV column only; it does not stop the run.

    A schedule step whose precondition is unmet (D or H with no previous S)
    SKIPS with a warning and continues; it never aborts -- uniformly, at
    every position in the schedule, including the first letter.

    D13 -- <output>/convergence.csv, one row per EXECUTED schedule step,
    FROZEN header (byte-identical to the other three drivers):

        cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time

    `step` is the schedule letter actually executed (S, D or H); `S` rows are
    always 'C' here (D8, no projected path); `D`/`H` rows carry the literal
    '-' (not tied to a solve path). `dofs` = lr.size(), `dofs_sigma`
    = domain.nControls(). `err` = rmse (the band quantity); minErr/maxErr are
    the pointwise extremes. `pctBelowTol` (fitting-driver convention, unlike
    the IGA drivers): percentage of DATA POINTS with error below Target.
    `minDetJsigma` = domain.minJacobian(), RECOMPUTED on every row (never
    carried forward). The per-direction counters nRefU/nRefV/nRefIso are NOT
    appended to this file (the header must stay byte-identical across
    drivers); they go to the run log and, when AnisoDump is on, to
    <output>/direction_dump.csv (header: cycle,elem,cu,cv,hu,hv,gu,gv,dir).
    A run that ends on D or H appends one more S row (step=S) so the S-row
    curve reaches the finest basis (see the "terminal S row" comment below).

    Under --plot the driver writes, in <output>: `geometry` (the composed map
    G = S o sigma), `sigma` (the sigma map), `mesh` (the LR element mesh),
    `points` (the input point cloud, fitting-specific), plus `options.xml`
    and `convergence.csv` (written UNCONDITIONALLY, so a run directory is
    reproducible on its own). There is NO `solution` ParaView collection in
    this driver: a fitting driver has no PDE field and no exact solution to
    animate across cycles, unlike the poisson/l2 drivers.

    ANISOTROPIC REFINEMENT (--aniso 1)
    ----------------------------------
    Direction criterion, ported from the LHB-spline paper (LHBSplines.tex,
    eq. (elementGradient)) and its driver
    gsRemappedBasis/examples/LHB_fitting_example.cpp:346-352,469-532:
    for a marked element e, integrate the fitted field's derivative over e,

        g_j = (1/|e|) * | int_e  dS_c / dxi_j |,      j in {u,v}

    and refine every direction whose g_j lies within a tie band of max_j g_j;
    ties in both directions fall back to isotropic. Three adaptations of the
    published form, each switchable so the paper form stays reproducible:

      (a) --anisoComp / AnisoComp. The paper fits a SCALAR height field, so
          grad(u_h) is unambiguous. Here S is a SURFACE whose first two
          components are close to the identity (dx/du ~ 1, dy/dv ~ 1); those
          O(1) constants would swamp the height signal and reduce the
          criterion to noise. The indicator is therefore driven by the
          height component (last row of xyz) by default. AnisoComp = -1 uses
          the norm over all components, for genuinely 3D point clouds.
      (b) --anisoAbsTol / AnisoAbsTol. gsClose(a,b,tol) is an ABSOLUTE test
          (gsMath.h:439), applied in the LHB driver to an unnormalised
          integral that scales like h^2. As elements shrink every g_j -> 0,
          all directions fall inside the tolerance, and the criterion
          silently degenerates to isotropic refinement on exactly the fine
          meshes where anisotropy should pay. Default here: divide by |e|
          (which is what the paper's text says -- "element-averaged" -- even
          though its formula omits it) and use the scale-invariant relative
          band g_j >= (1-eps)*max_j g_j.
      (c) --anisoAbs / AnisoAbs. The signed integrand vanishes for a feature
          symmetric within an element (a ridge crest centred in the cell),
          degenerating the choice; int_e |dS_c/dxi_j| is robust. Default is
          the signed (paper-faithful) form; the run log reports how often the
          two disagree.

    MESH-REPAIR ASYMMETRY BETWEEN THE TWO ARMS -- a confound, not a bug, and
    it must be quoted whenever the iso and aniso curves are compared. The
    isotropic path (gsLRSplineBasis::refineElements ->
    LRSplineSurface::refineBasisFunction, which starts at
    LRSplineSurface.cpp:814) ends with aPosterioriFixes() (called at
    LRSplineSurface.cpp:837, defined at LRSplineSurface.cpp:1130), which
    enforces maxTjoints_, closeGaps_ and maxAspectRatio_. The anisotropic path
    (refine_direction) calls m_lr->insert_line() directly and NEVER calls it.
    The two arms therefore differ in mesh repair as well as in direction. The
    per-H-step aspect-ratio summary printed by this driver is the measurement
    of that difference.

    COST OF AN H STEP (kept from the pre-unification header, still true).
    make_N2S2() ends with a compute() (gsLRFitting.h:190) before the explicit
    compute() that follows an H step's refine call, so an H step runs TWO
    least-squares solves and BOTH are timed into the CSV `time` column: LR
    refinement has no coefficient transfer here, so every refinement re-solves
    the fit from scratch. The ONLY difference between the isotropic and
    anisotropic arms is the refine call itself (refineElements vs
    refine_direction); everything else -- marking, N2S2, solve count -- is
    shared code.

    H ROWS CARRY POST-REFIT ERRORS: an H row pairs post-refinement DoFs with
    POST-REFIT errors (because the H step re-fits internally), unlike the THB
    driver whose H row pairs post-refinement DoFs with pre-refinement errors.

    TERMINAL S ROW: exactly as in the THB driver, a run that ends on a D or H
    step appends one more S row so that the S-row curve reaches the finest
    basis; because this driver's H already re-fits, that terminal S re-solves
    the same system at the same basis and xi and reports (near-)identical
    numbers -- a deliberate duplication so both drivers' curves are directly
    comparable at equal -i.

    Solver note: the S step here is gsFitting::compute, i.e. a BiCGSTABILUT
    ITERATIVE solve of the normal equations (LR bases have no
    gsExprAssembler path). Parameter correction is HARD-WIRED OFF
    (maxPcIter = 0, never touched by this driver): LR's per-point parameter
    correction would move xi_i away from sigma(uv_i) and fight the D step,
    confounding this driver's case with the baseline it is compared against.

    Expected input (same as fitting_example.cpp / the THB driver): an XML
    file containing
        Matrix id 0 : 2 x N parameter values (rescaled here to [0,1]^2)
        Matrix id 1 : 3 x N (or 2 x N) point coordinates

    Reference datasets (filedata/fitting/): ogre.xml (default), face.xml,
    nefertiti_lpsp.xml (user-provided point clouds), circle_band.xml (the
    non-axis-aligned case; lshape-style data is axis-aligned and understates
    the method).

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>   // gsOptFit (task 04)
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsOptim/gsOptim.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>

// !!! KEEP THIS INCLUDE LAST !!!
// gsLRFitting.h:20 does `#define d 2` and never #undef's it.  Any header
// included afterwards that declares `template<short_t d, ...>` is corrupted
// into `template<short_t 2, ...>`, and any identifier named `d` in the code
// below becomes the literal 2.  No single-letter `d` is used past this point.
// <gsAssembler/gsAdaptiveParametrization.h> (gsOptFit), <gsNurbs/gsSquareDomain.h>,
// <gsOptim/gsOptim.h> and <gsHLBFGS/gsHLBFGS.h> all declare or pull in
// template<short_t d, ...> and must therefore come BEFORE this include.
#include <gsLRSplines/gsLRFitting.h>

using namespace gismo;

namespace {

// D-step barrier floor: below it the barrier stops preventing folds, above
// it it distorts sigma. Promoted from a named constant to the option
// "BarrierEps" (task 18, same value 5e-2) so this driver's key set matches
// the other three rh drivers.

// Both optimizer backends are reachable through the unified option names
// (MaxIterations/OptTol/Verbose), mapped onto each backend's OWN option.
// NOTE: the unified tolerance key is deliberately called "OptTol", not
// "GradTol" -- gsHLBFGS already OWNS an option literally named "GradTol"
// (the line-search curvature constant, default 0.9), a completely different
// quantity; reusing that name here would silently clobber it
// (ORCHESTRATOR RULING 2).
std::unique_ptr<gsOptimizer<real_t> > makeOptimizer(const gsOptionList & o)
{
    const std::string which = o.askString("Optimizer", "HLBFGS");
    std::unique_ptr<gsOptimizer<real_t> > optimizer;
    if (which == "HLBFGS")
    {
        optimizer.reset(new gsHLBFGS<real_t>());
        optimizer->options().setReal("MinGradLen", o.getReal("OptTol"));
    }
    else if (which == "gsOptim")
    {
        optimizer.reset(new gsOptim<real_t>::LBFGS());
        optimizer->options().setReal("GradErrTol", o.getReal("OptTol"));
    }
    else
        GISMO_ERROR("Unknown Optimizer '" << which << "' (gsOptim | HLBFGS)");
    optimizer->options().setInt("MaxIterations", o.getInt("MaxIterations"));
    optimizer->options().setInt("Verbose",       o.getInt("Verbose"));
    return optimizer;
}

// Refinement direction codes of gsLRSplineBasis::refine_direction (taken
// from the CODE at gsLRSplineBasis.hpp:1329ff, NOT from the doc comment at
// gsLRSplineBasis.h:455, which contradicts it):
//   0 : isotropic (both directions)
//   1 : "vertical", i.e. const-u mesh lines splitting the largest u-span
//       -> adds resolution in u
//   2 : "horizontal", i.e. const-v mesh lines -> adds resolution in v
enum RefDir { DIR_ISO = 0, DIR_U = 1, DIR_V = 2 };

// Element aspect ratios h_max/h_min over the whole LR mesh; the diagnostic
// the LHB paper reads its own meshes with (LHBSplines.tex:2672,3106) and the
// measurement of the aPosterioriFixes asymmetry documented in the file header.
void aspectRatioStats(const gsLRSplineBasis<2, real_t> & lr,
                       real_t & mn, real_t & med, real_t & mx)
{
    const index_t nEl = static_cast<index_t>(lr.numElements());
    std::vector<real_t> ar;
    ar.reserve(nEl);
    gsVector<real_t> b;
    for (index_t e = 0; e != nEl; e++)
    {
        lr.getElementBounds(static_cast<std::size_t>(e), b);
        const real_t hu = b(2) - b(0), hv = b(3) - b(1);
        if (hu <= 0 || hv <= 0) continue;
        ar.push_back(std::max(hu, hv) / std::min(hu, hv));
    }
    if (ar.empty()) { mn = med = mx = 1; return; }
    std::sort(ar.begin(), ar.end());
    mn  = ar.front();
    mx  = ar.back();
    med = ar[ar.size() / 2];
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    // ------------------------------------------------------------------
    // CLI (D6-amendment: --options replaces the original spec's -f; this
    // driver has no separate PROBLEM file beyond the point cloud)
    // ------------------------------------------------------------------
    std::string dataFile = "fitting/ogre.xml";
    std::string output   = "rh_lrfitting_output";
    bool plot = false;
    index_t degree = 2, initialRef = 2;
    index_t sigmaDeg = 2, sigmaRef = 3;
    std::string schedule = "SH";
    index_t maxIter = 3;
    index_t aniso = 0;         // 0: isotropic (baseline), 1: per-element anisotropic
    bool project = false;      // rejection stub only, see D8
    bool coarsenFlag = false;  // rejection stub only, no LR unrefine (see below)
    real_t targetCli = std::numeric_limits<real_t>::quiet_NaN();
    std::string optFile = "";

    gsCmdLine cmd("Schedule-driven adaptive fitting of a parametrized point "
                  "cloud with a composed LR-spline geometry G = S(sigma(u)). "
                  "The --schedule string over {S,D,H} is repeated up to -i "
                  "times. LR parameter correction is always off (it would "
                  "fight sigma).");
    cmd.addString("d", "data",     "Input point-cloud XML", dataFile);
    cmd.addString("o", "output",   "Output directory", output);
    cmd.addSwitch("",  "plot",     "Write ParaView files", plot);
    cmd.addInt   ("p", "degree",   "Degree of the LR basis in both directions", degree);
    cmd.addInt   ("r", "refine",   "Initial refinement LEVEL of the S basis", initialRef);
    cmd.addInt   ("E", "sigmaDeg", "Degree of the sigma map", sigmaDeg);
    cmd.addInt   ("R", "sigmaRef", "Refinement LEVEL of sigma's mesh", sigmaRef);
    cmd.addString("",  "schedule", "Cycle string over {S,D,H}", schedule);
    cmd.addInt   ("i", "iter",     "Maximum number of schedule cycles", maxIter);
    cmd.addInt   ("",  "aniso",    "H-step refinement: 0 = isotropic (baseline), "
                                    "1 = per-element anisotropic", aniso);
    // D8: --project is a REJECTION STUB only. See the file header. A user
    // copying a THB command line gets an explanation instead of "unknown
    // option".
    cmd.addSwitch("",  "project",  "NOT AVAILABLE for LR bases -- rejected, "
                                    "see the file header (D8)", project);
    cmd.addSwitch("",  "coarsen",  "NOT AVAILABLE for LR bases -- rejected: "
                                    "gsLRSplineBasis has no unrefine (see the "
                                    "'coarsening requested but...' message in hStep)",
                                    coarsenFlag);
    cmd.addReal  ("",  "target",   "H step: override the error-band Target (rmse "
                                    "units) regardless of the XML/--options value; "
                                    "NaN (default, omit the flag) means 'no override, "
                                    "use the option list'", targetCli);
    cmd.addString("",  "options",  "Method option list XML (layered over the "
                                    "driver defaults, unknown keys ignored)", optFile);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(!project,
        "rh-adaptive_lrfitting_example: --project is not available for LR "
        "bases: makeIntegrationBasis (src/gsAssembler/gsAdaptiveParametrization.h:435) "
        "takes two gsTensorBSplineBasis and an LR basis has no tensor level.");
    GISMO_ENSURE(!coarsenFlag,
        "rh-adaptive_lrfitting_example: --coarsen is not available for LR bases: "
        "gsLRSplineBasis::unrefine/unrefineElements are GISMO_NO_IMPLEMENTATION "
        "(optional/gsLRSplines/src/gsLRSplineBasis.hpp:1263,1485). The Coarsen "
        "option is declared but IGNORED for the same reason (see the XML header).");

    // -r and -R are EXPONENTS (2^r elements per direction, via the (1<<r)-1
    // knot-count expressions below). 1<<n is UB for n>=31 and silently
    // ALIASES for n>=32 (measured: -r 34 behaves as -r 2), so a typo'd sweep
    // exponent would otherwise produce plausible numbers for the wrong mesh.
    // GISMO_ENSURE, not GISMO_ASSERT: asserts are inert in build_rel. 20 is
    // already 2^20 elements per direction -- far past any runnable size.
    GISMO_ENSURE(initialRef >= 0 && initialRef <= 20,
                 "-r/--refine must be in [0,20] (got " << initialRef << ")");
    GISMO_ENSURE(sigmaRef   >= 0 && sigmaRef   <= 20,
                 "-R/--sigmaRef must be in [0,20] (got " << sigmaRef << ")");
    GISMO_ENSURE(sigmaDeg >= 1, "-E/--sigmaDeg must be at least 1 (got " << sigmaDeg << ")");

    GISMO_ENSURE(degree > 0, "Degree must be positive");
    GISMO_ENSURE(maxIter > 0, "iter must be positive");
    GISMO_ENSURE(aniso == 0 || aniso == 1, "--aniso must be 0 or 1");

    // Schedule alphabet (D4): the full shared alphabet is S R D H U; this
    // driver implements ONLY S, D, H.
    std::string schedUp;
    for (char c : schedule)
    {
        const char C = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (C == 'R')
            GISMO_ERROR("rh-adaptive_lrfitting_example: no R step. Fitting "
                "r-adaptivity is the error-driven D step; the monitor-driven "
                "Winslow relocation exists only in the PDE drivers.");
        if (C == 'U')
            GISMO_ERROR("rh-adaptive_lrfitting_example: no U step. "
                "gsLRSplineBasis::uniformRefine is GISMO_NO_IMPLEMENTATION "
                "(optional/gsLRSplines/src/gsLRSplineBasis.hpp:1497), so uniform "
                "refinement of the LR basis is not available; use -r to set "
                "the initial refinement level instead.");
        if (C == 'P' || C == 'F')
            GISMO_ERROR("rh-adaptive_lrfitting_example: '" << c << "' is not "
                "an alias in the unified alphabet; the primal step is S.");
        GISMO_ENSURE(C == 'S' || C == 'D' || C == 'H',
            "rh-adaptive_lrfitting_example: invalid schedule letter '" << c
            << "' (this driver implements S, D, H)");
        schedUp.push_back(C);
    }

    if (output.empty() || output.back() != gsFileManager::getNativePathSeparator())
        output += gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);

    // Opened early (before the options section) so id=4 of -d can be
    // layered in at the same point as --options, giving every per-case
    // point-cloud file (2026-08-19) the same self-sufficiency poisson's -f
    // PDE files already have: -d ALONE, no --options, fully determines the
    // run's method options too.
    gsFileData<> fd_in(dataFile);

    // ------------------------------------------------------------------
    // Method options (D6): built-in defaults, layered with id=4 of -d, then
    // --options (highest precedence -- an explicit --options always wins)
    // ------------------------------------------------------------------
    gsOptionList opt;
    opt.addInt   ("RefineRule",     "Refinement rule (D7, inline): 1=GARU, 2=PUCA, 3=BULK", 3);
    opt.addReal  ("RefineParam",    "Refinement rule parameter", 0.5);
    opt.addInt   ("CoarsenRule",    "Coarsening rule -- IGNORED, no LR coarsening", 3);
    opt.addReal  ("CoarsenParam",   "Coarsening rule parameter -- IGNORED, idem", 0.1);
    opt.addSwitch("Coarsen",        "Enable coarsening -- IGNORED, idem", false);
    opt.addInt   ("MaxLevel",       "Level cap -- IGNORED, LR elements carry no level", 10);
    opt.addSwitch("Admissible",     "Admissible closure -- IGNORED, no LR admissibility", true);
    opt.addSwitch("Extension",      "Marker box extension -- IGNORED, LR marks without a halo", true);
    opt.addReal  ("Smoothing",      "R-step Winslow theta -- IGNORED, no R step here", 1.0);
    opt.addReal  ("Penalty",        "R-step Winslow fold penalty -- IGNORED, no R step here", 1e-3);
    opt.addReal  ("BarrierMu",      "D-step fold-barrier weight mu (renamed from Penalty, task 18)", 1000.0);
    opt.addInt   ("BarrierMode",    "D-step fold barrier: 0 = sampled Gauss rule on sigma's "
                                     "knot mesh (quB), 1 = Bezier-coefficient barrier (no quadrature)", 0);
    opt.addReal  ("BarrierEps",     "D-step fold barrier: det J_sigma floor (promoted from the "
                                     "named constant barrierEps, same value -- task 18)", 5e-2);
    opt.addReal  ("Target",         "Error band: absolute rmse target", 1e-2);
    opt.addReal  ("Band",           "Error band: dead-zone factor >= 1", 2.0);
    opt.addReal  ("NoCoarsenBelow", "Error band: below this rmse, hold (no LR coarsening happens)", 0.0);
    opt.addInt   ("MaxRefIt",       "Inner adapt iterations per H step", 1);
    opt.addReal  ("Lambda",         "gsFitting smoothing weight", 1e-6);
    opt.addString("Optimizer",      "D-step optimizer: gsOptim | HLBFGS", "HLBFGS");
    // task-08 measurement (shared across the fitting drivers): MaxIterations
    // = 10000 (the historical schedule-driver default) makes the D step a
    // guaranteed reject-and-restore NO-OP costing 170-290s (the long run
    // wanders into a fold, then gets discarded); MaxIterations = 50
    // converges in ~1.7s with a genuine ~34% error drop. More iterations
    // make the result WORSE here, so the shipped default stays in the
    // historical fitting-driver range (100).
    opt.addInt   ("MaxIterations",  "Optimizer iteration cap", 100);
    opt.addReal  ("OptTol",         "Optimizer gradient tolerance (unified name; NOT "
                                     "gsHLBFGS's own line-search GradTol)", 1e-6);
    opt.addInt   ("Verbose",        "Optimizer/driver verbosity", 0);
    opt.addString("MonitorMode",    "R-step monitor -- IGNORED, no R step here", "value");
    opt.addSwitch("Slide",          "Boundary controls of sigma slide along the boundary", true);
    // quA is DECLARED (for key-set parity, task 18) but not honoured: the
    // D-step barrier quadrature (gsFoldBarrier's Sampled mode, src/gsAssembler/
    // gsFoldBarrier.hpp) hard-codes
    // quadOptions.addReal("quA", "", 1.0) before calling gsQuadrature::get --
    // there is no consumer this driver could plumb a quA value into, so
    // claiming to honour one would be false advertising.
    opt.addReal  ("quA",            "Quadrature degree factor -- IGNORED, no projection path "
                                     "in this driver (--project rejected, D8); the D-step "
                                     "barrier hard-codes its own quA=1", 1.0);
    opt.addInt   ("quB",            "Quadrature: extra nodes, D-step barrier (tuned: 8; a "
                                     "DIFFERENT quantity from the IGA drivers' quB, the "
                                     "error/assembly quadrature)", 8);
    opt.addSwitch("SameElement",    "IGNORED -- no gsExprAssembler path in this driver", false);
    // Driver-specific (beyond D6's shared list): surviving aniso* knobs.
    opt.addReal  ("AnisoTol",       "Direction tie band (relative unless AnisoAbsTol)", 0.1);
    opt.addInt   ("AnisoComp",      "Target component for the direction indicator "
                                     "(-2: height, -1: norm, >=0: that component)", -2);
    opt.addSwitch("AnisoAbs",       "Integrate |dS/dxi| instead of |int dS/dxi|", false);
    opt.addSwitch("AnisoAbsTol",    "Paper-faithful unnormalised ABSOLUTE tie test", false);
    opt.addInt   ("AnisoQuB",       "Extra Gauss nodes per direction, direction indicator", 0);
    opt.addSwitch("AnisoDump",      "Dump a per-element direction record to "
                                     "<output>direction_dump.csv", false);
    opt.addInt   ("DirSkip",        "Skip D steps in the first DirSkip cycles", 0);
    opt.addSwitch("CheckGradient",  "FD-validate the D-step gradient once", false);
    opt.addInt   ("nSamples",       "ParaView sampling resolution", 1000);

    // Every key the driver actually knows about, captured BEFORE the file is
    // layered in, so a typo'd key in the XML can be named instead of silently
    // dropped (ignoreIfUnknown) or silently appended (addIfUnknown).
    std::set<std::string> knownKeys;
    for (const gsOptionList::OptionListEntry & e : opt.getAllEntries())
        knownKeys.insert(e.label);

    if (fd_in.hasId(4))
    {
        gsOptionList fo; fd_in.getId(4, fo);
        for (const gsOptionList::OptionListEntry & e : fo.getAllEntries())
            if (!knownKeys.count(e.label))
                gsInfo << "WARNING: option '" << e.label << "' from " << fd_in.lastPath()
                       << " (id=4) is not used by rh-adaptive_lrfitting_example (typo?)\n";
        opt.update(fo, gsOptionList::ignoreIfUnknown);
        gsInfo << "Method options from -d id=4: " << fd_in.lastPath() << "\n";
    }

    if (!optFile.empty())
    {
        gsFileData<> fdo(optFile);
        gsOptionList fo;
        GISMO_ENSURE(fdo.getFirst(fo), "No gsOptionList found in " << optFile);
        for (const gsOptionList::OptionListEntry & e : fo.getAllEntries())
            if (!knownKeys.count(e.label))
                gsInfo << "WARNING: option '" << e.label << "' from " << fdo.lastPath()
                       << " is not used by rh-adaptive_lrfitting_example (typo?)\n";
        opt.update(fo, gsOptionList::ignoreIfUnknown);
        gsInfo << "Method options from --options: " << fdo.lastPath() << "\n";
    }
    // --target replaces Target outright (any sign), applied LAST so it wins
    // over both the built-in default and any XML -- lets a sweep script
    // switch the band on/off and tune it per run without touching the XML.
    // It ALSO derives NoCoarsenBelow = Target/Band^2 (D9, 2026-08-18): the
    // dead zone no longer stops the run (only "hold" does), so
    // NoCoarsenBelow must sit meaningfully below Target/Band for "hold" to
    // ever be reachable -- there is deliberately no separate
    // --noCoarsenBelow flag (D6, "not too large a number of flags").
    // (--coarsen has no equivalent line here: it is rejected above, LR has
    // no unrefine.)
    if (!math::isnan(targetCli))
    {
        opt.setReal("Target", targetCli);
        const real_t bandFactor = opt.getReal("Band");
        opt.setReal("NoCoarsenBelow", targetCli / (bandFactor * bandFactor));
    }
    {
        const std::string checkWhich = opt.askString("Optimizer", "HLBFGS");
        GISMO_ENSURE(checkWhich == "HLBFGS" || checkWhich == "gsOptim",
                     "Unknown Optimizer '" << checkWhich << "' (gsOptim | HLBFGS)");
    }
    GISMO_ENSURE(opt.getInt("BarrierMode") == 0 || opt.getInt("BarrierMode") == 1,
                 "BarrierMode must be 0 (sampled) or 1 (coefficient), got "
                 << opt.getInt("BarrierMode"));

    gsInfo << opt;
    gsFileData<> fdout; fdout << opt; fdout.save(output + "options");

    gsInfo << "Ignored (this driver): CoarsenRule/CoarsenParam/Coarsen (no LR "
              "unrefine), MaxLevel/Admissible/Extension (no LR level or "
              "hierarchy), Smoothing/MonitorMode (no R step), SameElement "
              "(no gsExprAssembler path).\n";

    const std::string which = opt.askString("Optimizer", "HLBFGS");
    gsInfo << "Optimizer: " << which << " | MaxIterations = " << opt.getInt("MaxIterations")
           << " | OptTol = " << opt.getReal("OptTol")
           << " | quB = " << opt.getInt("quB") << "\n";

    //! [Read data]
    // id 0: 2 x N parameter values; id 1: dim x N point coordinates (fd_in
    // was already opened above, before the options section)
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv);
    fd_in.getId<gsMatrix<> >(1, xyz);
    GISMO_ENSURE(uv.cols() == xyz.cols() && uv.rows() == 2, "Wrong input");
    const index_t N = uv.cols();

    // Rescale the parameters to the sigma domain [0,1]^2
    for (index_t row = 0; row != 2; row++)
    {
        const real_t mn = uv.row(row).minCoeff(), mx = uv.row(row).maxCoeff();
        GISMO_ENSURE(mx > mn, "Degenerate parameter range in direction " << row);
        uv.row(row) = (uv.row(row).array() - mn) / (mx - mn);
    }
    //! [Read data]

    //! [Create bases and sigma]
    // LR basis for S on [0,1]^2. Held as the CONCRETE type on purpose: only
    // gsLRSplineBasis<2,real_t>::refineElements(gsVector<index_t>) exists --
    // it does NOT override gsBasis::refineElements(std::vector<index_t>), so
    // through a gsBasis& one lands on GISMO_NO_IMPLEMENTATION (a real
    // throw). -r is applied AT CONSTRUCTION, as a LEVEL (like -R), never
    // through uniformRefine (which throws on this class). The ctor takes
    // NON-const knot-vector references: named lvalues required.
    const index_t nKnots = (1 << initialRef) - 1;
    gsKnotVector<> u_knots(0, 1, nKnots, degree + 1);
    gsKnotVector<> v_knots(0, 1, nKnots, degree + 1);
    gsLRSplineBasis<2> lr(u_knots, v_knots);

    // Sigma map (identity initially). Sigma's mesh is a refinement LEVEL of
    // the level-0 analysis basis, not a free knot count: a non-nested pair
    // costs roughly the SUM of their elements in any union mesh, a nested
    // pair only the finer of the two.
    gsKnotVector<real_t> sigmaKv(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
    gsTensorBSplineBasis<2,real_t> sigmaBasis(sigmaKv, sigmaKv);
    gsSquareDomain<real_t> domain(sigmaBasis, opt.getSwitch("Slide"));

    // Upper clamp of xi. The LR element lookup used for marking
    // (gsLRSplineBasis::getSingleElement -> LR::LRSplineSurface::getElementContaining)
    // is half-open, and the sanity check that follows it (gsLRFitting.h:123)
    // is a plain assert(), i.e. gone under NDEBUG -- a point exactly on u=1
    // or v=1 would feed an out-of-range element index straight into
    // refineElements. gsTHBSplineBasis is safe there because
    // gsKnotVector::uFind handles the right endpoint; LR is not. The offset
    // is 1e-12, twelve orders below any fitting error, so the composed
    // geometry and the THB comparison are unaffected by the clamp.
    const real_t xiMax = 1.0 - 1e-12;
    //! [Create bases and sigma]

    gsInfo << "----------------\n";
    gsInfo << "Fitting " << N << " samples from " << fd_in.lastPath() << "\n";
    gsInfo << "Schedule           : " << schedUp << " x " << maxIter << "\n";
    gsInfo << "S basis (LR)       : deg " << degree << ", " << lr.size() << " DoFs\n";
    gsInfo << "Sigma basis        : deg " << sigmaDeg << ", " << domain.nControls()
           << " free controls\n";
    gsInfo << "H refinement       : " << (aniso == 0 ? "isotropic" : "ANISOTROPIC") << "\n";
    if (aniso != 0)
        gsInfo << "NOTE: the anisotropic refine path (refine_direction) does NOT\n"
                  "      call aPosterioriFixes(), which the isotropic path does.\n"
                  "      The two arms differ in mesh repair as well as direction;\n"
                  "      quote the aspect-ratio summary when comparing curves.\n";
    gsInfo << "----------------\n";

    // ------------------------------------------------------------------
    // State shared by the S/D/H steps
    // ------------------------------------------------------------------
    gsMatrix<> xi;                        // xi_i = sigma(uv_i), 2 x N
    std::vector<real_t> errors(N, 0.0);   // pointwise errors ||G(uv_i)-xyz_i||
    real_t minErr = 0, maxErr = 0, rmse = 0;

    // ONE persistent fitting object for the whole run (sigma is the identity
    // at construction, so uv is the correct initial parameter matrix).
    // make_N2S2() flips the object member `vertical_expansion`
    // (gsLRFitting.h:196), so rebuilding this object per S step would
    // silently pin every H step to one expansion direction; updateParams()
    // below instead writes xi into the persistent object every S/D step.
    gsLRFitting<real_t> fit(uv, xyz, lr, 0.1, opt.getReal("Lambda"));

    auto updateParams = [&]()
    {
        domain.eval_into(uv, xi);
        xi = xi.cwiseMax(0.0).cwiseMin(xiMax);
        // gsFitting caches the parameters: never let them drift from sigma(uv)
        fit.returnParamValues() = xi;
    };

    // Mirror the fitting object's own point errors into the local statistics
    auto syncErrors = [&]()
    {
        const std::vector<real_t> & pe = fit.pointWiseErrors();
        GISMO_ENSURE(static_cast<index_t>(pe.size()) == N,
                     "Fitting object reported " << pe.size() << " errors, expected " << N);
        errors.assign(pe.begin(), pe.end());
        minErr = fit.minPointError();
        maxErr = fit.maxPointError();
        rmse = 0;
        for (index_t i = 0; i != N; i++)
            rmse += errors[i] * errors[i];
        rmse = math::sqrt(rmse / N);
    };

    // pctBelowTol (fitting drivers): percentage of DATA POINTS below Target.
    const real_t target = opt.getReal("Target");
    auto pctBelow = [&]() -> real_t
    {
        index_t n = 0;
        for (real_t e : errors) if (e < target) n++;
        return 100.0 * n / static_cast<real_t>(errors.size());
    };

    // S: penalized least-squares fit of S at xi = sigma(uv), fixed LR basis.
    auto sStep = [&](index_t cycle)
    {
        gsStopwatch timer;
        updateParams();
        fit.compute(opt.getReal("Lambda"));
        // BAIL-OUT SEMANTICS OF gsFitting::compute -- read this before
        // trusting the checks below. compute() `delete m_result` FIRST
        // (gsFitting.hpp:83-85) and only then assembles and solves; if the
        // ILUT preconditioner fails it prints "The preconditioner failed.
        // Aborting." and RETURNS WITHOUT NULLING m_result (gsFitting.hpp:
        // 125-130): the object is left holding a DANGLING pointer, not the
        // previous geometry. The size check below only discriminates when
        // lr.size() changed since the last successful fit; on a fixed-basis
        // S (S after D) it would silently pass on stale data. THE RELIABLE
        // DETECTOR IS THE RUN LOG: grep -c "preconditioner failed" <log> and
        // discard any run that trips it.
        GISMO_ENSURE(fit.result() != nullptr, "LR least-squares fit produced no result");
        GISMO_ENSURE(fit.result()->basis().size() == lr.size(),
                     "LR fit result has " << fit.result()->basis().size()
                     << " coefficients but the basis has " << lr.size()
                     << ": the least-squares solve bailed out (preconditioner?)");
        GISMO_ENSURE(fit.result()->coefs().allFinite(),
                     "LR least-squares fit produced non-finite coefficients");
        fit.computeErrors();
        syncErrors();
        const real_t time = timer.stop();

        gsInfo << "  S | S DoFs: " << std::setw(6) << lr.size()
               << " | sigma DoFs: " << std::setw(4) << domain.nControls()
               << " | err (min/max): " << std::scientific << std::setprecision(3)
               << minErr << " / " << maxErr
               << " | rmse: " << rmse
               << " | min det Js: " << std::scientific << std::setprecision(2)
               << domain.minJacobian(7)
               << " | " << std::fixed << std::setprecision(2) << time << "s\n"
               << std::defaultfloat;
    };

    // Direction criterion of the H step (--aniso 1). Which target component
    // drives the indicator: see adaptation (a) in the file header.
    const index_t targetDim = xyz.rows();
    const index_t anisoComp = opt.getInt("AnisoComp");
    const index_t drvComp = (anisoComp == -2) ? targetDim - 1
                          : (anisoComp == -1) ? -1
                                              : anisoComp;
    GISMO_ENSURE(drvComp < targetDim,
                 "AnisoComp " << anisoComp << " exceeds the data dimension " << targetDim);

    index_t nRefU = 0, nRefV = 0, nRefIso = 0;
    index_t nVariantDisagree = 0, nDirTotal = 0;
    index_t curCycle = 0;
    std::ofstream dumpCsv;
    const bool anisoDump = opt.getSwitch("AnisoDump");
    const real_t anisoTol = opt.getReal("AnisoTol");
    const bool anisoAbs = opt.getSwitch("AnisoAbs");
    const bool anisoAbsTol = opt.getSwitch("AnisoAbsTol");
    const index_t anisoQuB = opt.getInt("AnisoQuB");

    auto elementDirection = [&](const gsVector<real_t> & b,
                                 real_t & guOut, real_t & gvOut) -> index_t
    {
        gsVector<index_t> nnodes(2);
        nnodes.setConstant(degree + 1 + anisoQuB);
        gsGaussRule<real_t> qr(nnodes);
        gsVector<real_t> lo(2), up(2), wts;
        lo << b(0), b(1);
        up << b(2), b(3);
        gsMatrix<real_t> pts;
        qr.mapTo(lo, up, pts, wts);

        gsMatrix<real_t> D;
        fit.result()->deriv_into(pts, D);

        gsVector<real_t> gSig(2), gAbs(2);
        for (index_t j = 0; j != 2; j++)
        {
            real_t accSig = 0, accAbs = 0;
            for (index_t q = 0; q != pts.cols(); q++)
            {
                real_t v;
                if (drvComp >= 0)
                    v = D(2 * drvComp + j, q);
                else
                {
                    real_t s = 0;
                    for (index_t c = 0; c != targetDim; c++)
                        s += D(2 * c + j, q) * D(2 * c + j, q);
                    v = math::sqrt(s);
                }
                accSig += wts[q] * v;
                accAbs += wts[q] * math::abs(v);
            }
            gSig[j] = math::abs(accSig);
            gAbs[j] = accAbs;
        }

        if (!anisoAbsTol)
        {
            const real_t area = (b(2) - b(0)) * (b(3) - b(1));
            if (area > 0) { gSig /= area; gAbs /= area; }
        }

        auto decide = [&](const gsVector<real_t> & g) -> index_t
        {
            const real_t gmax = g.maxCoeff();
            bool keepU, keepV;
            if (anisoAbsTol)
            {
                keepU = gsClose(g[0], gmax, anisoTol);
                keepV = gsClose(g[1], gmax, anisoTol);
            }
            else
            {
                if (gmax <= 0) return DIR_ISO;
                keepU = g[0] >= (1 - anisoTol) * gmax;
                keepV = g[1] >= (1 - anisoTol) * gmax;
            }
            if (keepU && keepV) return DIR_ISO;
            return keepU ? DIR_U : DIR_V;
        };

        const index_t dSig = decide(gSig), dAbs = decide(gAbs);
        nDirTotal++;
        if (dSig != dAbs) nVariantDisagree++;
        const gsVector<real_t> & gUsed = anisoAbs ? gAbs : gSig;
        guOut = gUsed[0];
        gvOut = gUsed[1];
        return anisoAbs ? dAbs : dSig;
    };

    // D: direct minimization of the true LS fitting error over the sigma
    // controls, via gsOptFit (task 04; variable projection with S: this
    // cannot increase the LS error beyond line-search tolerance, since the
    // optimizer descends on the reported objective itself). The fold
    // barrier's representation is selected by BarrierMode (0 = sampled Gauss
    // rule on sigma's own knot mesh, gone is only the legacy uniform grid;
    // 1 = Bezier-coefficient barrier).
    bool dBarrierPrinted = false;
    real_t minDetJ = 1.0;
    auto dStep = [&](index_t cycle)
    {
        gsStopwatch timer;
        const gsFoldBarrierMode dirMode = opt.getInt("BarrierMode") == 0
            ? gsFoldBarrierMode::Sampled : gsFoldBarrierMode::Coefficient;
        gsOptFit<real_t> obj(domain, *fit.result(), uv, xyz,
                              opt.getReal("BarrierMu"), opt.getReal("BarrierEps"),
                              dirMode, opt.getInt("quB"));
        if (!dBarrierPrinted)
        {
            gsInfo << "    [D] barrier quadrature: " << obj.nBarrierPts()
                   << " points ("
                   << (opt.getInt("BarrierMode") == 1 ? "Bezier-coefficient barrier"
                                                       : "Gauss on the sigma knot mesh")
                   << ")\n";
            dBarrierPrinted = true;
        }
        if (opt.getSwitch("CheckGradient"))
        {
            gsInfo << "    [D diag] max rel FD gradient error = "
                   << gsCheckSigmaGradient(obj, domain) << "\n";
            opt.setSwitch("CheckGradient", false); // first D step only
        }
        std::unique_ptr<gsOptimizer<real_t> > solver = makeOptimizer(opt);
        solver->setProblem(&obj);
        const gsVector<real_t> uPrev = domain.getControls();
        const real_t rmsePrev = rmse;
        solver->solve(uPrev);
        gsVector<real_t> uNew = solver->currentDesign();
        domain.setControls(uNew);
        minDetJ = domain.minJacobian(7);
        if (!uNew.allFinite() || minDetJ <= 0)
        {
            gsWarn << "Direct step rejected (min det Js = " << minDetJ
                   << "); sigma restored.\n";
            domain.setControls(uPrev);
            minDetJ = domain.minJacobian(7);
        }
        updateParams();
        // sigma moved but S did not: report the errors of the CURRENT S at
        // the new xi (the next S does the re-fit) while keeping the fitting
        // object's own error state in sync, so a following H marks against
        // the same numbers.
        fit.updateGeometry(fit.result()->coefs(), xi);
        syncErrors();
        if (rmse > rmsePrev && rmsePrev > 0)
            gsInfo << "    [D diag] rmse increased (" << rmsePrev << " -> " << rmse
                   << "): fold/box penalty active or line-search tolerance.\n";
        const real_t time = timer.stop();

        gsInfo << "  D | S DoFs: " << std::setw(6) << lr.size()
               << " | rmse: " << std::scientific << std::setprecision(3) << rmse
               << " | min det Js: " << std::scientific << std::setprecision(2) << minDetJ
               << " | " << std::fixed << std::setprecision(2) << time << "s\n"
               << std::defaultfloat;
        (void)cycle;
    };

    // H: LR refinement of the elements holding over-threshold points, then a
    // re-fit (mark -> [direction] -> refine -> make_N2S2() -> compute()).
    // MaxRefIt wraps the inner mark+refine iteration (band logic, D9).
    auto hStep = [&](index_t cycle) -> bool
    {
        bool convergedLocal = false;
        const real_t band = opt.getReal("Band"), nocrs = opt.getReal("NoCoarsenBelow");
        for (index_t inner = 0; inner < opt.getInt("MaxRefIt"); inner++)
        {
            gsStopwatch timer;
            bool doRef = true, doCrs = false, inBand = false, holdReached = false;
            const char * branch = "-";
            const real_t refTol = target * band, crsTol = target / band;
            GISMO_ENSURE(band >= 1, "Band must be >= 1");
            if      (rmse > refTol) { doRef = true;  doCrs = false; branch = "above"; }
            else if (rmse > crsTol) { doRef = true;  doCrs = true;  branch = "band";  inBand = true; }
            else if (rmse > nocrs)  { doRef = false; doCrs = true;  branch = "below"; }
            else                    { doRef = false; doCrs = false; branch = "hold";  inBand = true; holdReached = true; }

            if (doCrs)
                gsInfo << "    [H] coarsening requested (err " << branch << ") but "
                          "gsLRSplineBasis has no unrefine:\n"
                          "        unrefine / unrefineElements are GISMO_NO_IMPLEMENTATION\n"
                          "        (optional/gsLRSplines/src/gsLRSplineBasis.hpp:1263, :1485).\n"
                          "        The basis is held; only the D step can still reduce the error.\n";

            bool refined = false;
            if (doRef)
            {
                // Per-element error for a point cloud: MAXIMUM pointwise
                // error over the data points falling in that element (D7);
                // 0 for elements with no data point.
                const index_t nEl = static_cast<index_t>(lr.numElements());
                std::vector<real_t> elErr(nEl, 0.0);
                gsMatrix<real_t> pt(2, 1);
                for (index_t i = 0; i != N; i++)
                {
                    pt = xi.col(i);
                    const index_t el = lr.getSingleElement(pt);
                    GISMO_ENSURE(el >= 0 && el < nEl,
                                 "Point " << i << " at (" << pt(0,0) << "," << pt(1,0)
                                 << ") maps to LR element " << el << ", outside [0," << nEl << ")");
                    elErr[el] = std::max(elErr[el], errors[i]);
                }

                // D7 -- inline GARU/PUCA/BULK on a sorted (index,error) list.
                std::vector<std::pair<index_t,real_t> > sorted;
                sorted.reserve(nEl);
                for (index_t e = 0; e != nEl; e++)
                    if (elErr[e] > 0) sorted.push_back(std::make_pair(e, elErr[e]));
                std::sort(sorted.begin(), sorted.end(),
                          [](const std::pair<index_t,real_t> & a,
                             const std::pair<index_t,real_t> & b) { return a.second < b.second; });

                std::vector<index_t> elems;
                const index_t rule = opt.getInt("RefineRule");
                const real_t param = opt.getReal("RefineParam");
                if (!sorted.empty())
                {
                    if (rule == 1) // GARU: err >= param * maxError, from the top
                    {
                        const real_t maxErrEl = sorted.back().second;
                        for (auto it = sorted.rbegin(); it != sorted.rend(); ++it)
                        {
                            if (it->second < param * maxErrEl) break;
                            elems.push_back(it->first);
                        }
                    }
                    else if (rule == 2) // PUCA: floor(param * numElements) elements
                    {
                        index_t nMark = static_cast<index_t>(
                            std::floor(param * static_cast<real_t>(sorted.size())));
                        nMark = std::min(nMark, static_cast<index_t>(sorted.size()));
                        for (index_t k = 0; k != nMark; k++)
                            elems.push_back(sorted[sorted.size() - 1 - k].first);
                    }
                    else // BULK (default rule==3): cumulate from the top until
                         // param * sum(all errors)
                    {
                        real_t total = 0;
                        for (const auto & p : sorted) total += p.second;
                        real_t acc = 0;
                        for (auto it = sorted.rbegin(); it != sorted.rend(); ++it)
                        {
                            elems.push_back(it->first);
                            acc += it->second;
                            if (acc >= param * total) break;
                        }
                    }
                }

                if (elems.empty())
                {
                    // nEl here is sorted.size() -- the number of elements
                    // carrying at least one data point (0-error elements are
                    // never marking candidates), not lr.numElements().
                    const index_t nEl = static_cast<index_t>(sorted.size());
                    gsInfo << "  H | WARNING: RefineRule " << rule << " marked 0 of " << nEl
                           << " elements (RefineParam = " << param << ")";
                    if (rule == 2) gsInfo << "; PUCA marks floor(RefineParam*nElements), which is "
                                             "0 below " << (param>0 ? 1.0/param : 0.0) << " elements";
                    gsInfo << "\n";
                }

                if (!elems.empty())
                {
                    refined = true;
                    nRefU = nRefV = nRefIso = 0;
                    if (aniso == 0)
                    {
                        gsVector<index_t> els(elems.size());
                        for (size_t k = 0; k != elems.size(); k++) els[k] = elems[k];
                        gsInfo << "Refining " << elems.size() << " elements.\n";
                        lr.refineElements(els);
                    }
                    else
                    {
                        // Per-element direction, then a vote onto the
                        // support functions (getElementSupport, task 05):
                        // refine_direction addresses BASIS FUNCTIONS, whose
                        // supports straddle several elements. Conflicting
                        // votes fall back to isotropic (conservative: it can
                        // only add resolution, never withhold it).
                        std::map<index_t, index_t> vote;
                        gsVector<real_t> b;
                        gsVector<index_t> fns;
                        for (size_t k = 0; k != elems.size(); k++)
                        {
                            lr.getElementBounds(static_cast<std::size_t>(elems[k]), b);
                            real_t gu = 0, gv = 0;
                            const index_t dir = elementDirection(b, gu, gv);
                            if      (dir == DIR_U) nRefU++;
                            else if (dir == DIR_V) nRefV++;
                            else                   nRefIso++;

                            if (anisoDump && dumpCsv.is_open())
                            {
                                const real_t cu = 0.5 * (b(0) + b(2));
                                const real_t cv = 0.5 * (b(1) + b(3));
                                const real_t hu = b(2) - b(0);
                                const real_t hv = b(3) - b(1);
                                dumpCsv << curCycle << "," << elems[k] << ","
                                        << cu << "," << cv << "," << hu << "," << hv << ","
                                        << gu << "," << gv << "," << dir << "\n";
                            }

                            // ORCHESTRATOR ADDENDUM (task 05 review): only
                            // call getElementSupport AFTER refinement/ID
                            // generation has run at least once; here the
                            // first H step is preceded by at least one S, so
                            // this is safe.
                            lr.getElementSupport(static_cast<std::size_t>(elems[k]), fns);
                            for (index_t f = 0; f != fns.size(); f++)
                            {
                                GISMO_ASSERT(fns[f] >= 0, "getElementSupport returned a "
                                             "negative id -- basis IDs not generated yet?");
                                std::map<index_t, index_t>::iterator vit = vote.find(fns[f]);
                                if (vit == vote.end())       vote[fns[f]] = dir;
                                else if (vit->second != dir) vit->second = DIR_ISO;
                            }
                        }
                        gsInfo << "Refining " << elems.size() << " elements (u: " << nRefU
                               << ", v: " << nRefV << ", iso: " << nRefIso << ") through "
                               << vote.size() << " basis functions.\n";

                        gsVector<index_t> fnv(vote.size()), dirs(vote.size());
                        index_t k = 0;
                        for (std::map<index_t, index_t>::const_iterator it = vote.begin();
                             it != vote.end(); ++it, ++k)
                        {
                            fnv[k]  = it->first;
                            dirs[k] = it->second;
                        }
                        lr.refine_direction(fnv, dirs);
                    }

                    // Exactly nextIteration's tail: make_N2S2() ends with a
                    // compute() (gsLRFitting.h:190) and one more compute()
                    // follows -- TWO solves, do not collapse them (see the
                    // file header cost-model note).
                    fit.make_N2S2();
                    fit.compute(opt.getReal("Lambda"));
                    GISMO_ENSURE(fit.result() != nullptr &&
                                 fit.result()->basis().size() == lr.size(),
                                 "LR re-fit after refinement bailed out (result has "
                                 << (fit.result() ? fit.result()->basis().size() : 0)
                                 << " coefficients, basis has " << lr.size() << ")");
                    fit.computeErrors();
                    syncErrors();
                }
            }
            const real_t time = timer.stop();

            gsInfo << "  H | branch " << branch << " | marked " << (refined ? "some" : "none")
                   << " elements | LR elements " << lr.numElements()
                   << " | " << std::fixed << std::setprecision(2) << time << "s\n"
                   << std::defaultfloat;
            if (refined && aniso != 0)
            {
                real_t arMin, arMed, arMax;
                aspectRatioStats(lr, arMin, arMed, arMax);
                gsInfo << "    [H] split u/v/iso: " << nRefU << "/" << nRefV << "/" << nRefIso
                       << " | aspect ratio (min/med/max): " << std::fixed << std::setprecision(2)
                       << arMin << " / " << arMed << " / " << arMax << "\n" << std::defaultfloat;
            }
            (void)cycle;

            // inBand breaks the INNER loop for both "band" and "hold"; only
            // "hold" also stops the OUTER schedule (see the D9 comment).
            if (inBand) { convergedLocal = holdReached; break; }
            if (!refined) break; // nothing left to do
        }
        return convergedLocal;
    };

    // ------------------------------------------------------------------
    // Adaptive loop
    // ------------------------------------------------------------------
    std::ofstream csv(output + "convergence.csv");
    // D13: FROZEN header, byte-identical across the four drivers.
    csv << "cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time\n";
    csv << std::scientific << std::setprecision(8);

    if (opt.getSwitch("AnisoDump"))
    {
        dumpCsv.open(output + "direction_dump.csv");
        dumpCsv << "cycle,elem,cu,cv,hu,hv,gu,gv,dir\n";
    }

    bool converged = false;
    bool dirty = false;   // true when sigma or the basis changed since the last S
    for (index_t cyc = 0; cyc != maxIter && !converged; ++cyc)
    {
        curCycle = cyc;
        gsInfo << "----------------\n";
        gsInfo << "Cycle " << cyc << " [" << schedUp << "]\n";
        for (char C : schedUp)
        {
            if (converged) break;
            gsStopwatch stepTimer;
            if (C == 'S')
            {
                sStep(cyc);
                dirty = false;
                csv << cyc << ",S,C," << lr.size() << "," << domain.nControls() << ","
                    << minErr << "," << maxErr << "," << rmse << "," << pctBelow() << ","
                    << domain.minJacobian(7) << "," << stepTimer.stop() << "\n";
            }
            else if (C == 'D')
            {
                if (cyc < opt.getInt("DirSkip"))
                {
                    gsInfo << "  D | DirSkip active (cycle " << cyc << " < "
                           << opt.getInt("DirSkip") << "), skipping\n";
                    continue;
                }
                if (!fit.result())
                {
                    gsInfo << "  D | no previous solve, skipping\n";
                    continue;
                }
                dStep(cyc);
                dirty = true;
                csv << cyc << ",D,-," << lr.size() << "," << domain.nControls() << ","
                    << minErr << "," << maxErr << "," << rmse << "," << pctBelow() << ","
                    << minDetJ << "," << stepTimer.stop() << "\n";
            }
            else if (C == 'H')
            {
                if (!fit.result())
                {
                    gsInfo << "  H | no previous solve, skipping\n";
                    continue;
                }
                if (hStep(cyc)) converged = true;
                dirty = true;
                csv << cyc << ",H,-," << lr.size() << "," << domain.nControls() << ","
                    << minErr << "," << maxErr << "," << rmse << "," << pctBelow() << ","
                    << domain.minJacobian(7) << "," << stepTimer.stop() << "\n";
            }
        }
    }

    // TERMINAL S ROW (D13/D14): the schedule ended dirty (last executed
    // letter was D or H): conclude with one more S so the reported result --
    // and the S-row curve a harvester would build from it -- reflects the
    // final sigma AND the final basis, and the terminal DoF count and rmse
    // stay mutually consistent (an H "band" branch both refines and sets
    // convergence in the SAME inner iteration, so its post-refine mesh's
    // error has not been re-measured against the pre-refine numbers yet).
    if (dirty)
    {
        gsStopwatch timer;
        sStep(maxIter);
        csv << maxIter << ",S,C," << lr.size() << "," << domain.nControls() << ","
            << minErr << "," << maxErr << "," << rmse << "," << pctBelow() << ","
            << domain.minJacobian(7) << "," << timer.stop() << "\n";
    }
    csv.close();
    if (opt.getSwitch("AnisoDump")) dumpCsv.close();

    gsInfo << "----------------\n";
    gsInfo << "Finished: " << (converged ? "error band reached" : "iteration budget spent")
           << ". Final S DoFs: " << lr.size()
           << ", max error: " << maxErr << ", rmse: " << rmse << "\n";
    if (aniso != 0 && nDirTotal > 0)
    {
        real_t arMin, arMed, arMax;
        aspectRatioStats(lr, arMin, arMed, arMax);
        gsInfo << "Direction decisions: " << nDirTotal
               << " | signed vs |.|-integrated disagree on " << nVariantDisagree
               << " (" << (100.0 * nVariantDisagree / nDirTotal) << "%)\n";
        gsInfo << "Final element aspect ratio (min/med/max): "
               << arMin << " / " << arMed << " / " << arMax << "\n";
    }

    // ------------------------------------------------------------------
    // Output (D12: no `solution` collection -- no PDE field, no exact
    // solution to animate across cycles in a fitting driver)
    // ------------------------------------------------------------------
    if (plot && fit.result())
    {
        const index_t nSamples = opt.getInt("nSamples");
        gsComposedGeometry<real_t> geometry(domain.domain(), *fit.result());
        gsMultiPatch<> mp;
        mp.addPatch(geometry);
        gsWriteParaview(mp, output + "geometry", nSamples, true);
        gsWriteParaview(domain.domain(), output + "sigma", nSamples, true, true);
        gsMesh<real_t> lrMesh(lr, 0);
        gsWriteParaview(lrMesh, output + "mesh");
        gsWriteParaviewPoints(xyz, output + "points");
        gsInfo << "ParaView files written to " << output << "\n";
    }
    else if (!plot)
        gsInfo << "Done. Re-run with --plot to export S and sigma to ParaView.\n";

    return 0;
}
