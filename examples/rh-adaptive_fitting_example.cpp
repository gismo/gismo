/** @file rh-adaptive_fitting_example.cpp

    @brief Schedule-driven adaptive fitting of a parametrized point cloud
    with a composed geometry G = S o sigma.

    One of four "rh-driver-unification" drivers (poisson, l2projection,
    fitting, lrfitting) that share a 10-flag CLI core and a schedule
    alphabet, but implement only the letters their problem supports:

        letter  meaning                                poisson  l2  fitting  lrfitting
        S       solve / project / fit (the primal step)   x      x     x        x
        R       monitor-driven relocation of sigma        x      x     -        -
        D       error-driven relocation of sigma          -      x     x        x
        H       local refinement (THB / LR)                x      x     x        x
        U       uniform refinement                         x      x     x        -

    This driver implements S (fit), D (error-driven relocation of sigma) and
    H (THB refinement) and U (uniform refinement of the THB basis). It does
    NOT implement R: for point-cloud fitting the r-adaptivity IS the D step
    (there is no separate monitor-driven relocation), so R is rejected by
    name, and so are T (directional tensor refine, removed: THB
    uniformRefine(dir) is unsafe, see the U-step comment) and the legacy F/P
    aliases (this shared alphabet has no aliases). Every other letter is
    rejected too, with a message naming this driver and the reason.

    CLI (11 flags, identical meaning across all four drivers except D6-amendment
    below):

        -d, --data      Input point-cloud XML        (default "fitting/face.xml")
        -o, --output    Output directory              (default "rh_fitting_output")
            --plot      Write ParaView files
        -p, --degree    Degree of S in both directions (default 2)
        -r, --refine    Initial uniform refinement LEVEL of S (default 2):
                        S starts from gsKnotVector<>(0,1,0,p+1) in both
                        directions, then tbasis.uniformRefine((1<<r)-1), i.e.
                        2^r elements per direction (mirrors sigmaKnots=2^R-1).
        -E, --sigmaDeg  Degree of the sigma map        (default 2)
        -R, --sigmaRef  Refinement LEVEL of sigma's mesh (default 3)
            --schedule  Cycle string over {S,D,H,U}     (default "SDH")
        -i, --iter      Maximum number of schedule cycles (default 3)
            --project   Also measure the projected (P) path (D8)
            --options   Method option list XML (layered over the driver
                        defaults below via gsOptionList::update with
                        ignoreIfUnknown)

    ORCHESTRATOR RULING (D6-amendment, binding across all four drivers):
    -f/--file is reserved for a PROBLEM file; this driver has no separate
    problem definition beyond the point cloud (-d), so -f is dropped and
    --options takes its slot in the 11-flag budget (the original task spec
    called this flag "-f, --file"; the amendment supersedes that name).

    Driver keys carried on the resolved gsOptionList (name / default / use).
    This is one of four rh drivers sharing ONE unified 27-key option table
    (rh-driver-unification task 18); keys this driver does not implement are
    still DECLARED, with a "-- IGNORED, <reason>" desc, purely so the four
    reference XMLs carry the same key set and the "not used by ... (typo?)"
    warning stays silent across all four. NOTE (task 18): the D-step
    fold-barrier weight, formerly named "Penalty" here, is now "BarrierMu" --
    "Penalty" is reserved, across all four drivers, for the R-step Winslow
    fold penalty (which this driver does not implement, hence IGNORED below).

        RefineRule      3        gsHElementMarker refinement rule (1=GARU,2=PUCA,3=BULK)
        RefineParam     0.5      refinement parameter
        CoarsenRule     3        coarsening rule
        CoarsenParam    0.1      coarsening parameter
        Coarsen         false    enable coarsening in the H step
        MaxLevel        10       gsHElementMarker level cap
        Admissible      true     admissible closure
        Extension       true     marker box extension
        Lambda          1e-6     fitting smoothing weight (gsFitting::compute)
        Penalty         1e-3     IGNORED: R-step Winslow fold penalty, no R
                                  step in this driver (kept for key-set parity)
        BarrierMu       1000     D-step fold-barrier weight mu (renamed from
                                  "Penalty", same value -- task 18)
        Target          -1       error band: absolute rmse target (<0 disables the band)
        Band            2.0      error band: dead-zone factor >= 1
        NoCoarsenBelow  0.0      UNUSED since 2026-08-20 (see the band table)
        MaxRefIt        1        inner adapt iterations per H letter
        Optimizer       "HLBFGS" D-step optimizer backend: gsOptim | HLBFGS
        MaxIterations   100      optimizer iteration cap (see the D-step comment:
                                  large values are COUNTER-PRODUCTIVE here)
        OptTol          1e-6     optimizer gradient tolerance (unified name;
                                  mapped onto gsOptim's GradErrTol or gsHLBFGS's
                                  MinGradLen -- NOT gsHLBFGS's own "GradTol",
                                  which is the line-search curvature constant)
        Verbose         0        optimizer/driver verbosity
        MonitorMode     "value"  IGNORED: no R step in this driver
        Smoothing       1.0      IGNORED: no R step in this driver
        Slide           true     sigma boundary controls may slide
        quA             1.0      quadrature degree factor (projection, --project
                                  path only)
        quB             8        quadrature extra nodes (projection AND the
                                  D-step barrier -- ~10 nodes per element per
                                  direction are needed or HLBFGS folds sigma
                                  between the barrier's quadrature nodes; a
                                  DIFFERENT quantity from the IGA drivers'
                                  quB, the error/assembly quadrature)
        SameElement     false    quadrature: MUST be false under a composition
                                  (quadrature points cross element boundaries
                                  under sigma)
        BarrierEps      5e-2     D-step fold barrier: det J_sigma floor
        DirSkip         0        skip D steps in the first DirSkip cycles:
                                  adapting sigma against an S that resolves
                                  none of the features distorts the map and
                                  poisons later cycles

    D9 -- error-band rule for the H step, driven by the rmse of the
    last S step (mirrors l2projection_rh_schedule_example's band):

        rmse >  Target*Band                -> refine only, keep going
        Target/Band < rmse <= Target*Band  -> ONE refine+coarsen sweep, STOP
        rmse <= Target/Band                -> coarsen only, keep going

    The band is a target zone the H step drives the rmse INTO, from either
    side, and reaching it ends the run: pass the rmse you want the run to
    finish at as Target and it finishes there.

    2026-08-20 CHANGE: the dead zone used to refine as well, and only a
    fourth branch below it ("hold", at or under NoCoarsenBelow =
    Target/Band^2) stopped the run -- so a run had to overshoot its
    requested Target by a factor Band^2 before it was allowed to finish,
    and on a case whose error floor lies inside the band it never finished
    at all (measured here: rmse pinned at 5.87e-05 against a 1e-4 Target
    while the DoFs ran 4779 -> 27716, a 5.8x mesh for a 0.3% error gain).
    NoCoarsenBelow no longer takes part in the decision, and coarsening
    follows the branch instead of running on every H step: above the band
    there is nothing to give back, below it refinement is not what is
    called for, and inside it both act exactly once. Per-element marking is
    unchanged (gsHElementMarker::markCrs, relative to the CURRENT max
    element error).

    Target < 0 disables the band entirely: H always refines, and coarsens
    iff Coarsen is set.
    -i/--iter remains a hard cap regardless of the band.
    Band >= 1: LARGER means a WIDER dead zone (not the "bandwidth" of
    benchmark_Wrinkling_DWR.cpp).

    A schedule step whose precondition is unmet (D or H with no previous S)
    SKIPS with a warning and continues; it never aborts -- uniformly, at
    every position in the schedule, including the first letter (aborting on
    one position and skipping on another would be exactly the inconsistency
    a teaching example must not have).

    D13 -- <output>/convergence.csv, one row per EXECUTED schedule step (so
    an S with --project writes two rows, C then P), frozen header:

        cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time

    err is the rmse (the band quantity); minErr/maxErr are the pointwise
    min/max of ||G(uv_i)-x_i||; pctBelowTol (FITTING drivers only, unlike the
    IGA drivers) is the percentage of DATA POINTS with error below Target
    when Target >= 0, else the documented 1e-2 fallback; minDetJsigma is
    domain.minJacobian(), RECOMPUTED on every row (never carried forward).
    The closing "final S step when the schedule ended dirty" row (dirty =
    the schedule's last letter was D/H/U, not S) is kept, with step = S.
    path is 'C' (composed) or 'P' (projected) on S rows only; non-primal rows
    (D/H/U) are not tied to a solve path and carry the literal '-'.

    --project (D8) additionally, at every S step, L2-projects the composed
    map G = S o sigma onto the SAME THB basis as S (same DoF count -- the
    entire point of the comparison is "does the composed representation
    survive export to an ordinary CAD spline?") and re-measures the data
    error of the projected map at the ORIGINAL parameters uv; writes a
    second CSV row with path = 'P' for the same cycle/step (the composed row
    keeps path = 'C'). SameElement=false is mandatory for the projection
    quadrature too, for the same reason as above. gsProjection<L2,T>::project
    returns the INTEGRAL OF THE SQUARED error, not an L2 norm; the driver
    reports its square root as "projL2Err" so it is directly comparable to
    rmse (a norm). A schedule with no D step leaves sigma at the identity, so
    every P row of such a run is projErr ~ 0 / rmse(P) == rmse(C) by
    construction and is not a measurement of anything.

    Under --plot the driver writes, in <output>: `geometry` (the final
    composed map G = S o sigma), `sigma` (the sigma map, with mesh and
    control net), `mesh` (S with its THB element mesh), `solution` (one
    gsParaviewCollection, one timestep per S step), `points` (the input
    point cloud) plus `options.xml` and `convergence.csv`, which are written
    UNCONDITIONALLY so a run directory is reproducible on its own.

    Expected input (same as fitting_example.cpp): an XML file containing
        Matrix id 0 : 2 x N parameter values (rescaled here to [0,1]^2)
        Matrix id 1 : 3 x N (or 2 x N) point coordinates

    Reference datasets (filedata/fitting/): ogre.xml, face.xml (default),
    nefertiti_lpsp.xml, circle_band.xml (the non-axis-aligned case;
    lshape-style data is axis-aligned and understates the method).

    Cost note: a D step is one integrand sweep per objective/gradient call
    (no inner solve), because S is frozen while sigma is optimized --
    "variable projection" alternating with S. Sigma's mesh is kept
    dyadically nested with every THB analysis level (see the sigma
    construction comment below) so the paper's union integration basis stays
    cheap: a non-nested pair costs roughly the SUM of their elements, a
    nested pair only the finer of the two.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsHSplines/gsHElementMarker.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsOptim/gsOptim.h>
#include <gsUtils/gsProjection.h>   // gsL2Projection<T> alias, NOT gsUtils/gsL2Projection.h
                                    // (that file is a 0-byte legacy shim)

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>

using namespace gismo;

namespace {

// Both optimizer backends are reachable through the unified option names
// (MaxIterations/OptTol/Verbose), mapped onto each backend's OWN option --
// gsHLBFGS already owns an option literally named "GradTol" (the
// line-search curvature constant, default 0.9, a different quantity), so
// the unified tolerance key is deliberately NOT called GradTol.
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

} // namespace

int main(int argc, char *argv[])
{
    // ------------------------------------------------------------------
    // CLI (11 flags, D6-amendment: --options replaces the spec's -f)
    // ------------------------------------------------------------------
    std::string dataFile = "fitting/face.xml";
    std::string output   = "rh_fitting_output";
    bool plot = false;
    index_t degree = 2, initialRef = 2;
    index_t sigmaDeg = 2, sigmaRef = 3;
    std::string schedule = "SDH";
    index_t maxIter = 3;
    bool project = false;
    bool coarsenFlag = false;
    real_t targetCli = std::numeric_limits<real_t>::quiet_NaN();
    std::string optFile = "";

    gsCmdLine cmd("Schedule-driven adaptive fitting of a parametrized point cloud "
                  "with a composed THB geometry G = S(sigma(u)). The --schedule "
                  "string over {S,D,H,U} is repeated up to -i times.");
    cmd.addString("d", "data",     "Input point-cloud XML", dataFile);
    cmd.addString("o", "output",   "Output directory", output);
    cmd.addSwitch("",  "plot",     "Write ParaView files", plot);
    cmd.addInt   ("p", "degree",   "Degree of S in both directions", degree);
    cmd.addInt   ("r", "refine",   "Initial uniform refinement LEVEL of S", initialRef);
    cmd.addInt   ("E", "sigmaDeg", "Degree of the sigma map", sigmaDeg);
    cmd.addInt   ("R", "sigmaRef", "Refinement LEVEL of sigma's mesh", sigmaRef);
    cmd.addString("",  "schedule", "Cycle string over {S,D,H,U}", schedule);
    cmd.addInt   ("i", "iter",     "Maximum number of schedule cycles", maxIter);
    cmd.addSwitch("",  "project",  "Also measure the projected (P) path", project);
    cmd.addSwitch("",  "coarsen",  "H step: force Coarsen=true regardless of the "
                                    "XML/--options value (omit to use whatever the "
                                    "option list resolves to)", coarsenFlag);
    cmd.addReal  ("",  "target",   "H step: override the error-band Target (rmse "
                                    "units) regardless of the XML/--options value; "
                                    "NaN (default, omit the flag) means 'no override, "
                                    "use the option list'", targetCli);
    cmd.addString("",  "options",  "Method option list XML (layered over the "
                                    "driver defaults, unknown keys ignored)", optFile);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

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

    // Schedule alphabet (D4): the full shared alphabet is S R D H U; this
    // driver implements S, D, H, U. R is refused because for point-cloud
    // fitting the r-adaptivity IS the D step (no separate monitor-driven
    // relocation exists); T (directional tensor refine) is gone entirely
    // (see the U-step comment); F/P are not aliases in the shared alphabet.
    std::string schedUp;
    for (char c : schedule)
    {
        const char C = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        GISMO_ENSURE(C=='S' || C=='D' || C=='H' || C=='U',
            "rh-adaptive_fitting_example does not implement schedule step '" << c
            << "'. This driver implements S (fit), D (error-driven sigma relocation), "
               "H (THB refinement) and U (uniform refinement). "
            << (C=='R' ? "R (monitor-driven Winslow relocation) is not available for "
                         "point-cloud fitting: the fitting r-adaptivity IS the D step. "
               : C=='T' ? "T (directional tensor refinement) has been removed; use U. "
               : C=='F' ? "F is not an alias for S; the shared alphabet has no aliases. "
               : C=='P' ? "P is not an alias for S; the shared alphabet has no aliases. "
               : "")
            << "Alphabet: S R D H U.");
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
    opt.addInt   ("RefineRule",     "Refinement rule: 1=GARU, 2=PUCA, 3=BULK", 3);
    opt.addReal  ("RefineParam",    "Refinement rule parameter", 0.5);
    opt.addInt   ("CoarsenRule",    "Coarsening rule: 1=GARU, 2=PUCA, 3=BULK", 3);
    opt.addReal  ("CoarsenParam",   "Coarsening rule parameter", 0.1);
    opt.addSwitch("Coarsen",        "Enable coarsening in the H step", false);
    opt.addInt   ("MaxLevel",       "Maximum THB refinement level", 10);
    opt.addSwitch("Admissible",     "Keep the hierarchy admissible", true);
    opt.addSwitch("Extension",      "Extend marked element regions", true);
    opt.addReal  ("Lambda",         "Fitting smoothing weight", 1e-6);
    opt.addReal  ("Penalty",        "R-step Winslow fold penalty -- IGNORED, no R step in this driver", 1e-3);
    opt.addReal  ("BarrierMu",      "D-step fold-barrier weight mu (renamed from Penalty, task 18)", 1000.0);
    opt.addInt   ("BarrierMode",    "D-step fold barrier: 0 = sampled Gauss rule on sigma's "
                                     "knot mesh (quB), 1 = Bezier-coefficient barrier (no quadrature)", 0);
    opt.addReal  ("Target",         "Error band: absolute rmse target (<0 disables the band)", -1.0);
    opt.addReal  ("Band",           "Error band: dead-zone factor >= 1", 2.0);
    opt.addReal  ("NoCoarsenBelow", "Error band: UNUSED -- the band decision no "
                                    "longer reads it (2026-08-20)", 0.0);
    opt.addInt   ("MaxRefIt",       "Inner adapt iterations per H step", 1);
    opt.addString("Optimizer",      "D-step optimizer: gsOptim | HLBFGS", "HLBFGS");
    // task-08 measurement: MaxIterations=10000 (the historical schedule-driver
    // default) makes the D step a guaranteed reject-and-restore NO-OP costing
    // 170-290s (the long run wanders into a fold, then gets discarded);
    // MaxIterations=50 converges in ~1.7s with a genuine ~34% error drop. More
    // iterations make the result WORSE here, so the shipped default stays in
    // the historical fitting-driver range (100), not the schedule-driver one.
    opt.addInt   ("MaxIterations",  "Optimizer iteration cap", 100);
    opt.addReal  ("OptTol",         "Optimizer gradient tolerance (unified name)", 1e-6);
    opt.addInt   ("Verbose",        "Optimizer/driver verbosity", 0);
    opt.addSwitch("Slide",          "Boundary controls of sigma slide along the boundary", true);
    opt.addReal  ("quA",            "Quadrature: degree factor (projection)", 1.0);
    opt.addInt   ("quB",            "Quadrature: extra nodes (projection and D-step barrier)", 8);
    opt.addSwitch("SameElement",    "Quadrature: assume points share an element (must be "
                                     "false under a composition)", false);
    opt.addReal  ("BarrierEps",     "D-step fold barrier: det J_sigma floor", 5e-2);
    opt.addInt   ("DirSkip",        "Skip D steps in the first DirSkip cycles", 0);
    opt.addString("MonitorMode",    "R step monitor -- IGNORED, no R step in this driver", "value");
    opt.addReal  ("Smoothing",      "R step: Winslow smoothing theta -- IGNORED, no R step in this driver", 1.0);

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
                       << " (id=4) is not used by rh-adaptive_fitting_example (typo?)\n";
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
                       << " is not used by rh-adaptive_fitting_example (typo?)\n";
        // ignoreIfUnknown: a key from a shared XML that this driver does not
        // declare at all is silently dropped instead of being carried into
        // the dump as if it did something. Keys the driver does not
        // IMPLEMENT but is still expected to carry for key-set parity
        // (MonitorMode, Smoothing, Penalty -- R-step only) are declared
        // above with an IGNORED desc, so they survive this layering.
        opt.update(fo, gsOptionList::ignoreIfUnknown);
        gsInfo << "Method options from --options: " << fdo.lastPath() << "\n";
    }
    // CLI overrides, applied LAST so they win over both the built-in default
    // and any XML: --coarsen only ever turns Coarsen ON (there is no
    // --no-coarsen; omit the flag and rely on the XML/default for "off");
    // --target replaces Target outright (any sign), letting a sweep script
    // switch the band on/off and tune it per run without touching the XML.
    // Target is the rmse the run finishes at: the H step stops as soon as
    // the rmse is inside [Target/Band, Target*Band], so --target alone
    // fixes the endpoint and no second threshold has to be kept in step
    // with it.
    if (coarsenFlag) opt.setSwitch("Coarsen", true);
    if (!math::isnan(targetCli))
    {
        opt.setReal("Target", targetCli);
    }
    // Validate the Optimizer string HERE, at startup, not lazily at the
    // first D step: makeOptimizer() is only reached inside dStep, so an
    // unrecognised backend would otherwise complete a D-free run (e.g.
    // --schedule S) successfully, get printed as if resolved, and get
    // written into options.xml -- and only fail on an unrelated run that
    // happens to contain a D.
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

    if (opt.getSwitch("SameElement"))
        gsWarn << "WARNING: SameElement=1 is set, but this driver always evaluates "
                  "on a composed map: quadrature points will cross element boundaries "
                  "and results will be WRONG. Set SameElement=0.\n";

    const std::string which = opt.askString("Optimizer", "HLBFGS");
    gsInfo << "Optimizer: " << which << " | MaxIterations = " << opt.getInt("MaxIterations")
           << " | quA = " << opt.getReal("quA") << " | quB = " << opt.getInt("quB")
           << " | SameElement = " << opt.getSwitch("SameElement") << "\n";

    //! [Read data]
    // id 0: 2 x N parameter values; id 1: d x N point coordinates (fd_in
    // was already opened above, before the options section)
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv);
    fd_in.getId<gsMatrix<> >(1, xyz);
    GISMO_ENSURE(uv.cols() == xyz.cols() && uv.rows() == 2, "Wrong input");
    const index_t N = uv.cols();

    // Rescale the parameters to the sigma domain [0,1]^2
    for (index_t r = 0; r != 2; r++)
    {
        const real_t mn = uv.row(r).minCoeff(), mx = uv.row(r).maxCoeff();
        GISMO_ENSURE(mx > mn, "Degenerate parameter range in direction " << r);
        uv.row(r) = (uv.row(r).array() - mn) / (mx - mn);
    }
    //! [Read data]

    //! [Create bases and sigma]
    // THB basis for S on [0,1]^2, starting from a single element per
    // direction and doubling 2^r times (mirrors sigmaKnots = 2^R - 1 below).
    gsKnotVector<> ku(0, 1, 0, degree + 1);
    gsKnotVector<> kv(0, 1, 0, degree + 1);
    gsTensorBSplineBasis<2> tbasis(ku, kv);
    tbasis.uniformRefine((1 << initialRef) - 1);
    gsTHBSplineBasis<2> thb(tbasis);

    // Sigma map (identity initially). Sigma's mesh is a refinement LEVEL of
    // the level-0 analysis basis, not a free knot count. makeIntegrationBasis
    // (the paper's super mesh) unions the sigma and analysis knot lines; a
    // non-nested pair costs roughly the SUM of their elements in the super
    // mesh, a nested pair only the finer of the two -- cf. Sec. "Numerical
    // integration scheme": "if either the space associated with sigma or
    // the analysis space is nested within the other, no additional knot
    // lines are introduced". Every analysis level is dyadic, so tying sigma
    // to the same ladder makes a non-nested sigma inexpressible.
    // sbasis is also makeIntegrationBasis's second argument below.
    gsKnotVector<> ks(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
    gsTensorBSplineBasis<2> sbasis(ks, ks);
    gsSquareDomain<real_t> domain(sbasis, opt.getSwitch("Slide"));
    //! [Create bases and sigma]

    gsInfo << "----------------\n";
    gsInfo << "Fitting " << N << " samples from " << fd_in.lastPath() << "\n";
    gsInfo << "Schedule           : " << schedUp << " x " << maxIter << "\n";
    gsInfo << "S basis            : deg " << degree << ", " << thb.size() << " DoFs\n";
    gsInfo << "Sigma basis        : deg " << sigmaDeg << ", " << domain.nControls()
           << " free controls\n";
    gsInfo << "Project            : " << (project ? "yes" : "no") << "\n";
    gsInfo << "----------------\n";

    // ------------------------------------------------------------------
    // State shared by the S/D/H/U steps
    // ------------------------------------------------------------------
    gsGeometry<>::uPtr geom;              // current S
    gsMatrix<> xi;                        // xi_i = sigma(uv_i), 2 x N
    std::vector<real_t> errors(N, 0.0);   // pointwise errors ||G(uv_i)-xyz_i||
    real_t minErr = 0, maxErr = 0, rmse = 0;

    auto updateParams = [&]()
    {
        domain.eval_into(uv, xi);
        // Guard against numerical over/undershoot of sigma
        xi = xi.cwiseMax(0.0).cwiseMin(1.0);
    };

    auto computeErrors = [&]()
    {
        gsMatrix<> vals;
        geom->eval_into(xi, vals);
        minErr = std::numeric_limits<real_t>::max();
        maxErr = rmse = 0;
        for (index_t i = 0; i != N; i++)
        {
            errors[i] = (vals.col(i) - xyz.col(i)).norm();
            minErr = std::min(minErr, errors[i]);
            maxErr = std::max(maxErr, errors[i]);
            rmse += errors[i] * errors[i];
        }
        rmse = math::sqrt(rmse / N);
    };

    // pctBelowTol (fitting drivers): percentage of DATA POINTS below Target
    // when Target >= 0, else the documented 1e-2 fallback.
    const real_t target = opt.getReal("Target");
    const real_t pctThreshold = (target >= 0) ? target : 1e-2;
    gsInfo << "pctBelowTol threshold: " << pctThreshold
           << (target >= 0 ? " (Target)" : " (1e-2 fallback, Target < 0)") << "\n";
    auto pctBelow = [&](const std::vector<real_t> & errs) -> real_t
    {
        index_t n = 0;
        for (real_t e : errs) if (e < pctThreshold) n++;
        return 100.0 * n / static_cast<real_t>(errs.size());
    };

    std::ofstream csv(output + "convergence.csv");
    // Unbuffered: a sweep kills a run that overruns its time budget, and a
    // buffered stream loses every row it has not flushed -- a 2 h run then
    // yields a 0-byte file instead of the convergence history it did earn.
    csv << std::unitbuf;
    csv << "cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time\n";
    csv << std::scientific << std::setprecision(8);

    const index_t nSamples = 1000;
    std::unique_ptr<gsParaviewCollection> solcol;
    index_t solveStep = 0;
    if (plot)
        solcol.reset(new gsParaviewCollection(output + "solution"));

    // S: gsFitting penalized least-squares fit of S at xi = sigma(uv).
    // NUMERICAL CHANGE from the deleted hand-rolled leastSquaresFit: that
    // helper solved with the DIRECT SimplicialLDLT and a 1e-12 Tikhonov
    // diagonal shift; gsFitting<T>::compute solves with the ITERATIVE
    // BiCGSTABILUT and applies NO diagonal shift, and on preconditioner
    // failure leaves result() == nullptr (hence the GISMO_ENSURE below).
    auto sStep = [&](index_t cycle)
    {
        gsStopwatch timer;
        updateParams();
        gsFitting<real_t> fitter(xi, xyz, thb);
        fitter.compute(opt.getReal("Lambda"));
        GISMO_ENSURE(fitter.result() != nullptr,
                     "gsFitting::compute failed (BiCGSTABILUT preconditioner) at "
                     << thb.size() << " DoFs");
        geom = fitter.result()->clone();   // result() dies with fitter
        computeErrors();
        const real_t time = timer.stop();

        gsInfo << "  S | S DoFs: " << std::setw(6) << thb.size()
               << " | sigma DoFs: " << std::setw(4) << domain.nControls()
               << " | err (min/max): " << std::scientific << std::setprecision(3)
               << minErr << " / " << maxErr
               << " | rmse: " << rmse
               << " | min det Js: " << std::scientific << std::setprecision(2)
               << domain.minJacobian(7)
               << " | " << std::fixed << std::setprecision(2) << time << "s\n"
               << std::defaultfloat;
        csv << cycle << ",S,C," << thb.size() << "," << domain.nControls() << ","
            << minErr << "," << maxErr << "," << rmse << "," << pctBelow(errors) << ","
            << domain.minJacobian(7) << "," << time << "\n";

        if (plot)
        {
            // The single-gsGeometry overload of gsWriteParaview writes
            // exactly "solutionK.vts" (no per-patch "_0" suffix), matching
            // the filename the collection below references; the
            // gsMultiPatch overload used for the final "geometry" output
            // below names the file "..._0.vts" instead, which is fine there
            // because nothing else points at it.
            gsComposedGeometry<real_t> cgeom(domain.domain(), *geom);
            gsWriteParaview(cgeom, output + "solution" + std::to_string(solveStep), nSamples);
            solcol->addPart("solution" + std::to_string(solveStep) + ".vts", solveStep, "Solution");
            solveStep++;
        }

        // --project (D8): L2-project the composed map G onto the SAME THB
        // basis as S (same DoF count -- otherwise the "does the composed
        // representation survive CAD export" comparison is meaningless) and
        // re-measure the data error at the ORIGINAL parameters uv.
        // SameElement=false is mandatory here too: a quadrature element of
        // the integration basis does not map into a single element of the
        // analysis space under the composition.
        if (project)
        {
            gsStopwatch pTimer;
            gsTensorBSplineBasis<2> ibasis =
                gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>::
                    makeIntegrationBasis<2>(thb.tensorLevel(thb.maxLevel()), sbasis);
            gsComposedGeometry<real_t> cgeom(domain.domain(), *geom);
            gsMatrix<> pcoefs;
            gsOptionList projOpt = opt; // carries quA/quB/SameElement
            // gsProjection<L2,T>::project returns the INTEGRAL OF THE
            // SQUARED error (integral(sqNorm(s-f)*meas(G))), not an L2 norm
            // -- take the square root before comparing it against rmse
            // (which IS a norm), or the two look many orders apart when they
            // are actually the same order.
            const real_t projErrSq = gsL2Projection<real_t>::project(
                thb, ibasis, cgeom, pcoefs, gsBoundaryConditions<real_t>(), projOpt);
            // A negative value here can only mean broken quadrature (the integrand is
            // a squared norm); do NOT hide it behind math::abs(). Round-off can still
            // produce a tiny negative when the projection is near-exact, so clamp only
            // within round-off and fail loudly outside it.
            GISMO_ENSURE(projErrSq > -1e-12,
                         "gsL2Projection returned a negative squared error ("
                         << projErrSq << "): the projection quadrature is broken");
            const real_t projL2Err = math::sqrt(math::max(projErrSq, (real_t)0));
            pcoefs.resize(thb.size(), xyz.rows());
            gsGeometry<>::uPtr pgeo = thb.makeGeometry(give(pcoefs));

            gsMatrix<> pvals;
            pgeo->eval_into(uv, pvals);
            real_t pMinErr = std::numeric_limits<real_t>::max(), pMaxErr = 0, pRmse = 0;
            std::vector<real_t> pErrors(N);
            for (index_t i = 0; i != N; i++)
            {
                pErrors[i] = (pvals.col(i) - xyz.col(i)).norm();
                pMinErr = std::min(pMinErr, pErrors[i]);
                pMaxErr = std::max(pMaxErr, pErrors[i]);
                pRmse += pErrors[i] * pErrors[i];
            }
            pRmse = math::sqrt(pRmse / N);
            const real_t pTime = pTimer.stop();

            gsInfo << "  P | projL2Err: " << std::scientific << std::setprecision(3) << projL2Err
                   << " | rmse(P): " << pRmse
                   << " | " << std::fixed << std::setprecision(2) << pTime << "s\n"
                   << std::defaultfloat;
            csv << cycle << ",S,P," << thb.size() << "," << domain.nControls() << ","
                << pMinErr << "," << pMaxErr << "," << pRmse << "," << pctBelow(pErrors) << ","
                << domain.minJacobian(7) << "," << pTime << "\n";
        }
    };

    // D: direct minimization of the true LS fitting error over the sigma
    // controls (gsOptFit -- variable projection with S: this cannot
    // increase the LS error beyond line-search tolerance, unlike a
    // monitor-driven relocation, since the optimizer descends on the
    // reported objective itself). The fold barrier's representation is
    // selected by BarrierMode (0 = sampled Gauss rule; 1 = Bezier-coefficient
    // barrier).
    bool dBarrierPrinted = false;
    auto dStep = [&](index_t cycle)
    {
        gsStopwatch timer;
        const gsFoldBarrierMode dirMode = opt.getInt("BarrierMode") == 0
            ? gsFoldBarrierMode::Sampled : gsFoldBarrierMode::Coefficient;
        gsOptFit<real_t> obj(domain, *geom, uv, xyz,
                              opt.getReal("BarrierMu"), opt.getReal("BarrierEps"),
                              dirMode, opt.getInt("quB"));
        if (!dBarrierPrinted && opt.getInt("Verbose") > 0)
        {
            gsInfo << "    [D] barrier quadrature: " << obj.nBarrierPts()
                   << " points ("
                   << (opt.getInt("BarrierMode") == 1 ? "Bezier-coefficient barrier"
                                                       : "Gauss on the sigma knot mesh")
                   << ")\n";
            dBarrierPrinted = true;
        }
        std::unique_ptr<gsOptimizer<real_t> > solver = makeOptimizer(opt);
        solver->setProblem(&obj);
        const gsVector<real_t> uPrev = domain.getControls();
        const real_t rmsePrev = rmse;
        solver->solve(uPrev);
        gsVector<real_t> uNew = solver->currentDesign();
        domain.setControls(uNew);
        real_t minDetJ = domain.minJacobian(7);
        if (!uNew.allFinite() || minDetJ <= 0)
        {
            gsWarn << "Direct step rejected (min det Js = " << minDetJ
                   << "); sigma restored.\n";
            domain.setControls(uPrev);
            minDetJ = domain.minJacobian(7);
        }
        updateParams();
        computeErrors();
        if (rmse > rmsePrev && rmsePrev > 0)
            gsInfo << "    [D diag] rmse increased (" << rmsePrev << " -> " << rmse
                   << "): fold/box penalty active or line-search tolerance.\n";
        const real_t time = timer.stop();

        gsInfo << "  D | S DoFs: " << std::setw(6) << thb.size()
               << " | rmse: " << std::scientific << std::setprecision(3) << rmse
               << " | min det Js: " << std::scientific << std::setprecision(2) << minDetJ
               << " | " << std::fixed << std::setprecision(2) << time << "s\n"
               << std::defaultfloat;
        csv << cycle << ",D,-," << thb.size() << "," << domain.nControls() << ","
            << minErr << "," << maxErr << "," << rmse << "," << pctBelow(errors) << ","
            << minDetJ << "," << time << "\n";
    };

    // Maps a point in the S-domain to the (level, cell0, cell1) key of the
    // THB element that currently contains it (finest-level cell index,
    // dropped down to the containing element's own level via query3).
    auto elementKey = [&](real_t x, real_t y) -> std::array<index_t,3>
    {
        const int maxLvl = thb.maxLevel();
        const gsTensorBSplineBasis<2> & tB = thb.tensorLevel(maxLvl);
        gsVector<index_t,2> c;
        c(0) = tB.component(0).knots().uFind(x).uIndex();
        c(1) = tB.component(1).knots().uFind(y).uIndex();
        const int lvl = thb.tree().query3(c, c + gsVector<index_t,2>::Ones(), maxLvl);
        return { static_cast<index_t>(lvl), c(0) >> (maxLvl - lvl), c(1) >> (maxLvl - lvl) };
    };

    // Per-element errors from per-point errors (H step marking). Built ONCE
    // per call: elemErr[e] = sqrt(sum over points i in element e of
    // ||G(uv_i)-x_i||^2), 0 for elements with no data points, indexed by
    // it.id() (the order gsHElementMarker::setErrors expects). The
    // consistency GISMO_ENSURE (sum of squares reproduces N*rmse^2) is the
    // cheap check that the point-to-element mapping is right.
    auto elementErrors = [&]() -> std::vector<real_t>
    {
        std::map<std::array<index_t,3>, index_t> elemOf;
        for (auto it = thb.domain()->beginAll(); it != thb.domain()->endAll(); ++it)
        {
            const int lvl = static_cast<const gsHDomainIterator<real_t,2>*>(it.get())->getLevel();
            const gsTensorBSplineBasis<2> & tb = thb.tensorLevel(lvl);
            const index_t i0 = tb.component(0).knots().uFind(it.lowerCorner()(0)).uIndex();
            const index_t i1 = tb.component(1).knots().uFind(it.lowerCorner()(1)).uIndex();
            elemOf[{ static_cast<index_t>(lvl), i0, i1 }] = it.id();
        }
        std::vector<real_t> elemErrSq(thb.numElements(), 0.0);
        for (index_t i = 0; i != N; i++)
        {
            const std::array<index_t,3> key = elementKey(xi(0,i), xi(1,i));
            auto found = elemOf.find(key);
            GISMO_ENSURE(found != elemOf.end(),
                         "H step: data point " << i << " could not be located in any THB element");
            elemErrSq[found->second] += errors[i] * errors[i];
        }
        real_t sumSq = 0;
        for (real_t e : elemErrSq) sumSq += e;
        GISMO_ENSURE(math::abs(sumSq - N * rmse * rmse) <=
                     1e-6 * std::max(static_cast<real_t>(1), N * rmse * rmse),
                     "H step: per-element error aggregation does not reproduce N*rmse^2 ("
                     << sumSq << " vs " << N * rmse * rmse << "): point-to-element mapping is wrong.");
        std::vector<real_t> elemErr(elemErrSq.size());
        for (size_t k = 0; k != elemErrSq.size(); k++)
            elemErr[k] = math::sqrt(elemErrSq[k]);
        return elemErr;
    };

    const real_t band = opt.getReal("Band");
    const bool useBand = (target >= 0);
    GISMO_ENSURE(!useBand || band >= 1, "Band must be >= 1");
    const bool coarsenOpt = opt.getSwitch("Coarsen");

    // H: gsHElementMarker (task 02: markCrs dispatches on CoarsenRule
    // directly, no RefineRule save/restore workaround needed). Refining the
    // basis AND the geometry together keeps G = S o sigma unchanged (and
    // the cached errors valid) until the next S step. Both box lists are
    // computed on the SAME pre-update mesh and the two regions are disjoint
    // by construction.
    auto hStep = [&](index_t cycle) -> bool
    {
        bool convergedLocal = false;
        for (index_t inner = 0; inner < opt.getInt("MaxRefIt"); inner++)
        {
            gsStopwatch timer;
            // 2026-08-18: doCrs is UNCONDITIONAL on coarsenOpt -- see the
            // D9 header comment for why (coarsening is now a per-element,
            // every-loop criterion, decoupled from the global band).
            // 2026-08-20: the band is a TARGET ZONE the H step drives the
            // rmse INTO, from whichever side it starts on, and reaching it
            // ends the run:
            //
            //   rmse >  Target*Band   "above"  refine only, and keep going
            //   in the band           "band"   ONE combined refine+coarsen
            //                                  sweep, then stop
            //   rmse <= Target/Band   "below"  coarsen only, and keep going
            //
            // The previous rule refined inside the band too and only stopped
            // at Target/Band^2. On a case whose error floor lies inside the
            // band that never terminates -- measured here: rmse pinned at
            // 5.87e-05 against a 1e-4 Target while the DoFs ran 4779 ->
            // 27716, a 5.8x mesh for a 0.3% error gain. Coarsening now
            // follows the branch rather than running on every H step.
            bool doRef = true, doCrs = coarsenOpt, inBand = false, holdReached = false;
            const char * branch = "-";
            if (useBand)
            {
                const real_t refTol = target * band, crsTol = target / band;
                if      (rmse > refTol) { doRef = true;  doCrs = false;     branch = "above"; }
                else if (rmse > crsTol) { doRef = true;  doCrs = coarsenOpt; branch = "band";
                                          inBand = true; holdReached = true; }
                else                    { doRef = false; doCrs = coarsenOpt; branch = "below";
                    // Below the band with coarsening off, no action can
                    // raise the rmse back into it, so the loop would spin.
                    if (!coarsenOpt) { inBand = true; holdReached = true; } }
            }

            std::vector<index_t> boxes, crsBoxes;
            size_t nCrs = 0;
            if (doRef || doCrs)
            {
                std::vector<real_t> elemErr = elementErrors();
                gsHElementMarker<2,real_t> marker(thb);
                marker.options().setInt   ("RefineRule",   opt.getInt ("RefineRule"));
                marker.options().setReal  ("RefineParam",  opt.getReal("RefineParam"));
                marker.options().setInt   ("CoarsenRule",  opt.getInt ("CoarsenRule"));
                marker.options().setReal  ("CoarsenParam", opt.getReal("CoarsenParam"));
                marker.options().setInt   ("MaxLevel",     opt.getInt ("MaxLevel"));
                marker.options().setSwitch("Admissible",   opt.getSwitch("Admissible"));
                marker.options().setSwitch("Extension",    opt.getSwitch("Extension"));
                marker.setErrors(elemErr);
                typedef gsHElementMarker<2,real_t>::HElementContainer HElementContainer;
                const HElementContainer markedRef = doRef ? marker.markRef() : HElementContainer();
                if (doRef) boxes = marker.toRefBoxes(markedRef);
                if (doRef && markedRef.empty())
                {
                    const index_t nEl = thb.numElements();
                    const index_t rule = opt.getInt("RefineRule");
                    const real_t param = opt.getReal("RefineParam");
                    gsInfo << "  H | WARNING: RefineRule " << rule << " marked 0 of " << nEl
                           << " elements (RefineParam = " << param << ")";
                    if (rule == 2) gsInfo << "; PUCA marks floor(RefineParam*nElements), which is "
                                             "0 below " << (param>0 ? 1.0/param : 0.0) << " elements";
                    gsInfo << "\n";
                }
                if (doCrs)
                {
                    const HElementContainer markedCrs = marker.markCrs(markedRef);
                    nCrs = markedCrs.size();
                    crsBoxes = marker.toCrsBoxes(markedCrs);
                }
            }

            if (!boxes.empty())    { thb.refineElements(boxes);    geom->refineElements(boxes); }
            if (!crsBoxes.empty()) { thb.unrefineElements(crsBoxes); geom->unrefineElements(crsBoxes); }
            const real_t time = timer.stop();

            gsInfo << "  H | branch " << branch << " | marked " << boxes.size()
                   << " | coarsened " << nCrs << " | elements " << thb.numElements()
                   << " | " << std::fixed << std::setprecision(2) << time << "s\n"
                   << std::defaultfloat;
            if (boxes.empty() && crsBoxes.empty() && (doRef || doCrs))
                gsInfo << "    (no cells marked)\n";
            csv << cycle << ",H,-," << thb.size() << "," << domain.nControls() << ","
                << minErr << "," << maxErr << "," << rmse << "," << pctBelow(errors) << ","
                << domain.minJacobian(7) << "," << time << "\n";

            // Reaching the band breaks the INNER loop (there is no point
            // re-marking within the same H letter after its one
            // refine+coarsen sweep) and stops the OUTER schedule: the run
            // has arrived at the error it was asked for.
            if (inBand) { convergedLocal = holdReached; break; }
            if (boxes.empty() && crsBoxes.empty()) break; // nothing left to do
        }
        return convergedLocal;
    };

    // U: uniform refinement of the THB basis and the geometry together (so
    // G stays unchanged, exactly as the H step does). No direction argument
    // and no T step: gsHTensorBasis::uniformRefine(dir) routes through a
    // direction-agnostic m_tree.multiplyByTwo(), and its dir==-1 assert is
    // compiled out under NDEBUG (silent UB in a release build) -- that
    // hazard is the reason directional refinement is gone from this driver.
    auto uStep = [&](index_t cycle)
    {
        gsStopwatch timer;
        thb.uniformRefine();
        if (geom) geom->uniformRefine();
        const real_t time = timer.stop();
        gsInfo << "  U | S DoFs: " << std::setw(6) << thb.size()
               << " | " << std::fixed << std::setprecision(2) << time << "s\n"
               << std::defaultfloat;
        csv << cycle << ",U,-," << thb.size() << "," << domain.nControls() << ","
            << minErr << "," << maxErr << "," << rmse << "," << pctBelow(errors) << ","
            << domain.minJacobian(7) << "," << time << "\n";
    };

    // ------------------------------------------------------------------
    // Adaptive loop
    // ------------------------------------------------------------------
    bool converged = false;
    bool dirty = false;   // true when sigma or the basis changed since the last S
    for (index_t cyc = 0; cyc != maxIter && !converged; ++cyc)
    {
        gsInfo << "----------------\n";
        gsInfo << "Cycle " << cyc << " [" << schedUp << "]\n";
        for (char C : schedUp)
        {
            if (converged) break;
            if (C == 'S')
            {
                sStep(cyc);
                dirty = false;
            }
            else if (C == 'D')
            {
                if (cyc < opt.getInt("DirSkip"))
                {
                    gsInfo << "  D | DirSkip active (cycle " << cyc << " < "
                           << opt.getInt("DirSkip") << "), skipping\n";
                    continue;
                }
                if (!geom)
                {
                    gsInfo << "  D | no previous solve, skipping\n";
                    continue;
                }
                dStep(cyc);
                dirty = true;
            }
            else if (C == 'H')
            {
                if (!geom)
                {
                    gsInfo << "  H | no previous solve, skipping\n";
                    continue;
                }
                if (hStep(cyc)) converged = true;
                dirty = true;
            }
            else if (C == 'U')
            {
                uStep(cyc);
                dirty = true;
            }
        }
    }

    // The schedule ended dirty (D13): conclude with a fit so that the
    // reported result and the final CSV row reflect the final basis and
    // sigma. This includes the band's "converged" case: the "band" branch
    // of hStep both refines AND sets inBand/converged in the SAME inner
    // iteration, so the mesh that triggered convergence is a POST-refinement
    // mesh whose error has not been re-measured yet -- without this refit,
    // the terminal DoF count and the reported rmse would be mismatched (the
    // rmse would be the PRE-refinement value). A "hold" convergence (no
    // refine/coarsen at all) merely re-fits an unchanged basis, which is
    // harmless.
    if (dirty)
        sStep(maxIter);

    gsInfo << "----------------\n";
    gsInfo << "Finished: " << (converged ? "error band reached" : "iteration budget spent")
           << ". Final S DoFs: " << thb.size()
           << ", max error: " << maxErr << ", rmse: " << rmse << "\n";

    // ------------------------------------------------------------------
    // Output
    // ------------------------------------------------------------------
    if (plot)
    {
        solcol->save();
        if (geom)
        {
            gsComposedGeometry<real_t> cgeomFinal(domain.domain(), *geom);
            gsMultiPatch<> mpFinal; mpFinal.addPatch(cgeomFinal);
            gsWriteParaview(mpFinal, output + "geometry", nSamples, true);
            gsWriteParaview(*geom, output + "mesh", nSamples, true, true);
        }
        gsWriteParaview(domain.domain(), output + "sigma", nSamples, true, true);
        gsWriteParaviewPoints(xyz, output + "points");
        gsInfo << "ParaView files written to " << output << "\n";
    }

    return 0;
}
