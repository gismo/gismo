/** @file poisson_rh_schedule_example.cpp
    @brief Schedule-driven manufactured-solution Poisson rh-adaptivity driver.

    One of four "rh-driver-unification" drivers (poisson, l2projection,
    fitting, lrfitting) that share a 10-flag CLI core and a schedule alphabet,
    but implement only the letters their problem supports:

    | letter | meaning                                    | poisson |
    |--------|---------------------------------------------|---------|
    | S      | solve / project / fit (the primal step)      | yes     |
    | R      | monitor-driven relocation of sigma (Winslow) | yes     |
    | D      | error-driven relocation of sigma             | --      |
    | H      | local refinement (THB / LR)                  | yes     |
    | U      | uniform refinement                           | yes     |

    This driver does NOT implement D (an error-driven relocation of sigma for
    the Poisson problem needs an error estimator, which does not exist in
    src/ yet), nor the older P/F aliases used by the other three drivers --
    any letter outside S,R,H,U is rejected by name.

    The problem is defined ENTIRELY by the XML file (-f), using the same ID
    scheme as poisson2_example.cpp:

      id=0  MultiPatch   physical domain (one patch, tensor-product basis)
      id=1  FunctionExpr source term
      id=2  boundaryConditions
      id=3  FunctionExpr exact solution
      id=4  OptionList   ALL options (assembler + driver), layered over the
                          built-in defaults below (a partial list is valid)
      id=5  Geometry     composition geometry for sigma        [OPTIONAL,
                          overrides -E/-R when present]

    Option resolution: the built-in defaults below are layered with id=4 of
    -f (if present), and that result is layered AGAIN with --options <file>
    (if given) -- gsOptionList::update overwrites matching keys and appends
    unknown ones, so precedence is PER KEY, not per source: --options wins
    over id=4 wins over the built-in default for any key it sets, but a key
    that only id=4 sets (e.g. the file's own quA/quB/SameElement) survives
    even when --options is given and does not mention it. Either layering
    uses gsOptionList::addIfUnknown, so a PARTIAL list stays valid and
    unknown extra keys are simply carried along.

    Driver keys carried on the resolved gsOptionList (name / default / use).
    This is one of four rh drivers (poisson, l2, fitting, lrfitting) that
    share ONE unified 27-key option table (rh-driver-unification task 18);
    keys this driver does not implement are still DECLARED, with a
    "-- IGNORED, <reason>" desc, purely so the four reference XMLs carry the
    same key set and the "not used by ... (typo?)" warning stays silent
    across all four:

      RefineRule      3        gsHElementMarker refinement rule (1=GARU,2=PUCA,3=BULK)
      RefineParam      0.5     refinement parameter
      CoarsenRule     3        coarsening rule
      CoarsenParam     0.1     coarsening parameter
      Coarsen         false    enable coarsening in the H step
      MaxLevel        10       marker level cap (library default 3 is far too low here)
      Admissible      true     admissible closure
      Extension       true     marker box extension
      MaxRefIt        1        IGNORED: this driver does one mark/refine per H
                                letter; put more H letters in the schedule
                                instead
      Smoothing        1.0     R step, Winslow smoothing theta
      Penalty          1e-3    R step, fold penalty
      Target          -1       band target L2 error; < 0 disables the band (D9)
      Band             2.0     band width factor, must be >= 1
      NoCoarsenBelow  -1       UNUSED since 2026-08-20 (see the band table below)
      Optimizer       "gsOptim" R-step backend: "gsOptim" or "HLBFGS"
      MaxIterations   10000    optimizer iterations
      OptTol           1e-12   optimizer gradient tolerance (unified name; mapped
                                onto each backend's OWN tolerance option, see
                                makeOptimizer() -- gsHLBFGS already has a
                                different "GradTol", the line-search curvature
                                constant, so the unified key is deliberately
                                NOT called GradTol)
      Verbose         0        optimizer verbosity
      MonitorMode     "value"  R step monitor: "value" | "gradient" (string,
                                unified with the l2/lrfitting drivers)
      DirSkip         0        IGNORED: no D step in this driver (D-step
                                direction skip)
      Lambda          1e-6     IGNORED: no gsFitting in this driver
      Slide           true     sigma boundary control points may slide

    plus quA/quB/SameElement (and every other gsExprAssembler<>::defaultOptions()
    key), which are mirrored onto the error evaluators AND onto the R step's
    gsAdaptiveParametrization so the system, the norms that measure it and the
    R-step objective are all integrated identically -- quA/quB are pinned to
    2/2 as the driver's BUILT-IN default (the converged error quadrature; the
    old assembler default of 1/1 under-integrated the marking error, see the
    comment where opt is constructed). MaxRefIt is not honoured: one
    mark/refine happens per H letter; put more H letters in the schedule for
    more refinement rounds. PROVIDED an S separates them: the per-element
    error vector driving the marker is cached from the last S step (see the H
    branch), so an H with no S since the last mesh change (e.g. the second H
    of "SHH", or any H after a bare U) has nothing valid to mark against and
    is skipped with a message instead of refining stale/mismatched data;
    insert an S between H letters to get another marking/refinement round.

    D9 -- error-band rule for the H step (mirrors
    l2projection_rh_schedule_example's band, driven by the GLOBAL L2 error of
    the last S step):

      err >  Target*Band                -> refine only, keep going
      Target/Band < err <= Target*Band  -> ONE refine+coarsen sweep, then STOP
      err <= Target/Band                -> coarsen only, keep going

    The band is a target zone the H step drives the error INTO, from either
    side, and reaching it ends the run. Pass the error you want the run to
    finish at as Target and it finishes there.

    2026-08-20 CHANGE: the dead zone used to refine as well, and only a
    fourth branch below it ("hold", at or under NoCoarsenBelow = Target/Band^2)
    stopped the run -- so a run had to overshoot its requested Target by a
    factor Band^2 before it was allowed to finish, and on a case whose error
    floor lies inside the band it never finished at all (measured on the
    fitting driver: rmse pinned at 5.87e-05 against a 1e-4 Target while the
    DoFs ran 4779 -> 27716). NoCoarsenBelow no longer takes part in the band
    decision. Coarsening likewise follows the branch now rather than running
    on every H step.

    (b) Coarsening is now UNCONDITIONAL on the Coarsen switch -- it no
        longer depends on which branch the GLOBAL error falls into, and runs
        on EVERY H step (while doRef alone still follows the four branches
        above). It marks on the PER-ELEMENT error via
        gsHElementMarker::markCrs, relative to the CURRENT max element
        error: CoarsenRule=1 (GARU) with CoarsenParam=0.001, for instance,
        coarsens every element whose error is below 0.001x that step's max
        -- a genuinely local, every-loop criterion, independent of the
        global Target entirely. (The shipped default is CoarsenRule=3/BULK,
        CoarsenParam=0.1 -- change both via --options if GARU's
        max-relative threshold is what you want.)

    Target < 0 disables the band entirely: the H step always refines, and
    coarsens iff Coarsen is set (unchanged by (b) -- this was already
    unconditional in that case).
    -i/--iter remains a hard cap on the cycle count regardless of the band.

    D13 -- <output>/convergence.csv, one row per EXECUTED schedule step
    (so an S with --project writes two rows, C then P), frozen header:

      cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time

    path is 'C' (composed) or 'P' (projected) on S rows only; non-primal rows
    (R/D/H/U) are not tied to a solve path and carry the literal '-'.

    err = global L2 error of the row's own solve (S rows) / last known solve
    (other rows); minErr/maxErr = min/max per-ELEMENT L2 error contribution
    from gsExprEvaluator::integralElWise; pctBelowTol (matching the l2
    driver's formula) = percentage of elements whose per-element error is
    below Target/sqrt(nElements) -- the per-element share of a global target
    under equidistribution, 0 when Target < 0; minDetJsigma =
    sigma.minJacobian(), RECOMPUTED on every row (not carried forward).

    --project (D8) additionally solves, at every S step, on the L2 projection
    of the composed geometry onto the active analysis basis, and writes a
    second CSV row with path=P.

    Under --plot the driver writes, in <output>: `geometry` (the final
    composed map), `sigma` (the sigma map itself), `solution` (one
    gsParaviewCollection, numerical/exact/error fields, one timestep per S
    row of the CSV, in the same order) and `mesh` (the final analysis mesh).
    <output>/options.xml and <output>/convergence.csv are written
    UNCONDITIONALLY, so a run directory is reproducible on its own.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

// One S-step result: the extracted solution, its L2/H1 errors, and the
// per-element L2 error contributions on the ANALYSIS basis (needed by both
// the CSV min/max columns and the H-step marker).
struct SolveResult
{
    gsMultiPatch<> sol;
    real_t l2err = 0.0;
    real_t h1err = 0.0;
    std::vector<real_t> elErr; // per element, on mb's elements, SQUARE-ROOTED already
};

// Assemble and solve -Delta u = f with Dirichlet data u = exact on the map
// mp, on the analysis basis mb with integration elements ib. The solution is
// extracted, its L2/H1 errors and per-element L2 error contributions are
// computed, and (when \a collection is given) the current step is appended
// to a gsParaviewCollection as a new timestep with the numerical solution,
// the exact solution and the pointwise error fields.
// \a A and \a ev are a persistent assembler and a sharing solver evaluator
// (see gsExprEvaluator(const gsExprAssembler<>&)) created once in main().
// \a pev is a sharing PLOT/ERROR evaluator: its integration elements are set
// to the ANALYSIS basis \a mb here (not the merged integration basis \a ib),
// both so the collection element mesh (plotElements) shows the space's own
// elements, and -- the reason this matters beyond plotting -- because
// gsHElementMarker::setErrors() indexes its error vector by the id() of the
// elements of the basis the marker was constructed on (mb's basis), so the
// per-element error vector MUST have exactly that length and indexing; the
// union integration basis \a ib is a different (finer, super-meshed) basis
// and would silently produce a vector of the wrong size.
SolveResult solve(gsExprAssembler<>& A, gsExprEvaluator<>& ev,
                   gsExprEvaluator<>& pev,
                   const gsMultiPatch<>& mp, const gsMultiBasis<>& mb,
                   const gsMultiBasis<>& ib, const gsFunctionExpr<>& f,
                   const gsFunctionExpr<>& exact, gsBoundaryConditions<>& bc,
                   const char path, index_t step,
                   gsParaviewCollection* collection)
{
    // Quadrature/SameElement come from the file's id=4 (or --options) and are
    // set once on A and ev in main(), so the system and the norms that
    // measure it are integrated identically.
    A.setIntegrationElements(ib);
    ev.setIntegrationElements(ib);
    auto G=A.getMap(mp);
    auto u=A.getSpace(mb);
    auto ff=A.getCoeff(f,G);
    auto ue=ev.getVariable(exact, G); // exact composed with the map: physical values
    bc.setGeoMap(mp);
    gsMatrix<> x;
    auto us=A.getSolution(u,x);
    u.setup(bc,dirichlet::l2Projection,0);

    A.initSystem();

    gsStopwatch timer;
    A.assemble(igrad(u,G)*igrad(u,G).tr()*meas(G),u*ff*meas(G));
    const real_t tAssemble=timer.stop();

    timer.restart();
    typename gsSparseSolver<>::SimplicialLDLT solver(A.matrix());
    x=solver.solve(A.rhs());
    const real_t tSolve=timer.stop();

    SolveResult result;
    us.extract(result.sol);

    timer.restart();
    result.l2err=math::sqrt(ev.integral((us-ue).sqNorm()*meas(G)));
    // ue is bound to G, so igrad(ue) is the AMBIENT gradient in R^geoDim; on an
    // embedded geometry (geoDim != parDim) it must be projected onto the tangent
    // space before it can be compared with the tangential igrad(us,G).
    // P = jac(G)*jac(G).ginv() is that (symmetric) projector and reduces to the
    // identity whenever geoDim == parDim, so the planar cases are unaffected.
    result.h1err=result.l2err+math::sqrt(ev.integral(
        (igrad(ue)*jac(G)*jac(G).ginv()-igrad(us,G)).sqNorm()*meas(G)));
    const real_t tError=timer.stop();

    // Per-element L2 error contributions, on the ELEMENTS OF THE ANALYSIS
    // BASIS (mb) -- the marker indexes its error vector by the id() of
    // those elements. This is the SAME integrand as the reported L2 norm
    // above, so the quantity that steers refinement and the quantity that
    // is reported are one and the same (this replaces the old hand-rolled
    // 5x5-point-per-element max-error sample loop).
    pev.setIntegrationElements(mb);
    auto pG=pev.getMap(mp);
    auto pus=pev.getVariable(result.sol);
    auto pue=pev.getVariable(exact, pG);
    pev.integralElWise((pue-pus).sqNorm()*meas(pG));
    result.elErr = pev.elementwise();          // squared, per element
    for (real_t & e : result.elErr) e = math::sqrt(e); // element L2 error
    GISMO_ENSURE(result.elErr.size()==(size_t)mb.basis(0).numElements(),
                 "element error vector does not match the analysis basis");

    gsInfo << "  " << path << " | step " << std::setw(3) << step
           << " | dofs " << std::setw(6) << A.numDofs()
           << " | assemble " << std::fixed << std::setprecision(3) << tAssemble << " s"
           << " | solve " << tSolve << " s"
           << " | error " << tError << " s"
           << " | L2 " << std::scientific << std::setprecision(3) << result.l2err
           << " | H1 " << result.h1err
           << std::defaultfloat << "\n";

    if (collection)
    {
        collection->newTimeStep(const_cast<gsMultiPatch<>*>(&mp),
                                static_cast<real_t>(step));
        collection->addField(pus, "numerical solution");
        collection->addField(pue, "exact solution");
        collection->addField((pue-pus).sqNorm(), "error");
        collection->saveTimeStep();
    }
    return result;
}

// Both backends are reachable; the unified names (MaxIterations/OptTol/
// Verbose) are mapped onto each backend's own so that the 100-vs-10000
// default mismatch is visible in the log instead of hidden in a header.
// NOTE: the unified tolerance key is deliberately called "OptTol", not
// "GradTol" -- gsHLBFGS already OWNS an option literally named "GradTol"
// (the line-search curvature constant, default 0.9), a completely different
// quantity; reusing that name here would silently clobber it.
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

} // namespace

int main(int argc, char** argv)
{
    index_t degree=2, initialRef=2, iterations=3;
    index_t sigmaDeg=2, sigmaRef=3;
    bool project=false, plot=false, coarsenFlag=false;
    real_t targetCli=std::numeric_limits<real_t>::quiet_NaN();
    std::string schedule="SRHS";
    std::string file="pde/poisson2d_rh_center_bvp.xml";
    std::string options="";
    std::string output="poisson_rh_output";
    gsCmdLine cmd("Schedule-driven manufactured-solution Poisson rh-adaptivity example.");
    cmd.addString("f","file","PDE XML file (ids 0-5 as in poisson2_example, "
                             "optional id=5 = sigma composition geometry)",file);
    cmd.addString("o","output","output prefix/directory",output);
    cmd.addSwitch("plot","write the ParaView output set",plot);
    cmd.addInt("p","degree","analysis degree",degree);
    cmd.addInt("r","refine","initial uniform refinements",initialRef);
    cmd.addInt("E","sigmaDeg","degree of the sigma map",sigmaDeg);
    cmd.addInt("R","sigmaRef","uniform refinement LEVEL of sigma's mesh w.r.t. the "
               "LEVEL-0 analysis basis (2^sigmaRef elements per direction, so sigma "
               "is always dyadically nested with every analysis level)",sigmaRef);
    cmd.addString("","schedule","letters S,R,H,U",schedule);
    cmd.addInt("i","iter","schedule repetitions (hard cap)",iterations);
    cmd.addSwitch("project","also solve on the L2-projected composed geometry",project);
    cmd.addSwitch("coarsen","H step: force Coarsen=true regardless of the XML/--options "
                            "value (omit to use whatever the option list resolves to)",
                            coarsenFlag);
    cmd.addReal("","target","H step: override the error-band Target (rmse/L2 units) "
                            "regardless of the XML/--options value; NaN (default, omit "
                            "the flag) means 'no override, use the option list'",
                            targetCli);
    // ORCHESTRATOR RULING (D6-amendment): --options is the shared "method
    // option list" flag in all four drivers; -f means the PROBLEM only.
    // Poisson resolves options from the built-in defaults, layered with
    // id=4 of -f, layered AGAIN with --options when given -- per KEY, not
    // per source (see the header comment).
    cmd.addString("","options","method option list XML (layered PER KEY over "
                               "id=4 of -f, taking precedence for any key it "
                               "sets)",options);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

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

    if (output.empty() || output.back()!=gsFileManager::getNativePathSeparator())
        output += gsFileManager::getNativePathSeparator();

    GISMO_ENSURE(iterations>0, "invalid iteration count");

    // Validated BEFORE the output directory is created, so an invalid
    // schedule (unknown letter) fails fast without leaving an empty
    // directory behind.
    std::string sched;
    for (char c:schedule)
    {
        const char C=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        GISMO_ENSURE(C=='S'||C=='R'||C=='H'||C=='U',
            "poisson_rh_schedule_example does not implement schedule letter '"<<c<<"'. "
            "This driver implements S (solve), R (monitor relocation of sigma), "
            "H (local THB refinement) and U (uniform refinement). In particular there "
            "is no D step here: an error-driven relocation of sigma for the Poisson "
            "problem needs an error estimator, which does not exist in src/ yet.");
        sched+=C;
    }
    const bool useH=sched.find('H')!=std::string::npos;

    // Written unconditionally (not only under --plot), so a run directory is
    // reproducible from its own output.
    gsFileManager::mkdir(output);

    //! [Read input file] -- ID scheme shared with poisson2_example.cpp
    gsFileData<> fd(file);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";

    gsMultiPatch<> physical;
    fd.getId(0, physical);                  // id=0: physical domain
    GISMO_ENSURE(physical.nPatches()==1,"this example expects one patch");

    gsFunctionExpr<> f, ms;
    fd.getId(1, f);                         // id=1: source term
    gsInfo << "Source function "  << f  << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc);                        // id=2: boundary conditions

    fd.getId(3, ms);                        // id=3: exact solution
    gsInfo << "Exact solution "   << ms << "\n";

    // id=4 (or --options): ALL options -- assembler AND driver -- layered
    // over the built-in defaults, so a partial list stays valid. Pushed to
    // the error evaluators too (below), so the system and the norms that
    // measure it use the same quadrature.
    gsOptionList opt = gsExprAssembler<>::defaultOptions();
    // Pin the converged error quadrature as the driver's OWN built-in
    // default (rh-driver-unification task 18, defect 2): the
    // gsExprAssembler default above is quA=1/quB=1, which under-integrates
    // the marking error; shipped bvp files already override this at id=4,
    // but a run WITHOUT such a file must not silently fall back to 1/1.
    opt.setReal("quA", 2.0);
    opt.setInt ("quB", 2);
    opt.addInt   ("RefineRule",     "gsHElementMarker refinement rule (1=GARU,2=PUCA,3=BULK)", 3);
    opt.addReal  ("RefineParam",    "gsHElementMarker refinement parameter", 0.5);
    opt.addInt   ("CoarsenRule",    "gsHElementMarker coarsening rule", 3);
    opt.addReal  ("CoarsenParam",   "gsHElementMarker coarsening parameter", 0.1);
    opt.addSwitch("Coarsen",        "enable coarsening in the H step", false);
    opt.addInt   ("MaxLevel",       "gsHElementMarker level cap", 10);
    opt.addSwitch("Admissible",     "gsHElementMarker admissible closure", true);
    opt.addSwitch("Extension",      "gsHElementMarker box extension", true);
    opt.addInt   ("MaxRefIt",       "H step refinement iterations -- IGNORED, this driver does one mark/refine per H letter", 1);
    opt.addReal  ("Smoothing",      "R step: Winslow smoothing theta", 1.0);
    opt.addReal  ("Penalty",        "R step: fold penalty", 1e-3);
    opt.addReal  ("Target",         "H band: target L2 error (<0 disables the band)", -1);
    opt.addReal  ("Band",           "H band: dead-zone width factor (>=1)", 2.0);
    opt.addReal  ("NoCoarsenBelow", "H band: UNUSED -- the band decision no "
                                    "longer reads it (2026-08-20)", -1);
    opt.addString("Optimizer",      "R-step optimizer backend: gsOptim | HLBFGS", "gsOptim");
    opt.addInt   ("MaxIterations",  "R-step optimizer iterations", 10000);
    opt.addReal  ("OptTol",         "R-step optimizer gradient tolerance (unified name)", 1e-12);
    opt.addInt   ("Verbose",        "R-step optimizer verbosity", 0);
    opt.addString("MonitorMode",    "R step monitor: value | gradient", "value");
    opt.addInt   ("DirSkip",        "D-step direction skip -- IGNORED, no D step in this driver", 0);
    opt.addReal  ("Lambda",         "gsFitting smoothing weight -- IGNORED, no gsFitting in this driver", 1e-6);
    opt.addSwitch("Slide",          "sigma boundary control points may slide", true);
    // Every key the driver actually knows about, captured BEFORE the file is
    // layered in, so a typo'd key in the XML can be named instead of silently
    // dropped (ignoreIfUnknown) or silently appended (addIfUnknown).
    std::set<std::string> knownKeys;
    for (const gsOptionList::OptionListEntry & e : opt.getAllEntries())
        knownKeys.insert(e.label);
    // Layer id=4 of -f first, then --options on top: gsOptionList::update
    // overwrites keys it finds and appends the rest (addIfUnknown), so this
    // gives PER-KEY precedence (--options > id=4 > built-in default) instead
    // of discarding one source wholesale when the other is present -- a file
    // like id=4's SameElement=0/quA=2/quB=2 (mandatory under the sigma
    // composition, see the Quadrature comment below) must survive even when
    // --options only sets an unrelated key such as Target.
    if (fd.hasId(4))
    {
        gsOptionList fo; fd.getId(4, fo);
        for (const gsOptionList::OptionListEntry & e : fo.getAllEntries())
            if (!knownKeys.count(e.label))
                gsInfo << "WARNING: option '" << e.label << "' from " << fd.lastPath()
                       << " is not used by poisson_rh_schedule_example (typo?)\n";
        opt.update(fo, gsOptionList::addIfUnknown);
    }
    if (!options.empty())
    {
        gsFileData<> fdOpt(options);
        gsOptionList fo;
        GISMO_ENSURE(fdOpt.getFirst<gsOptionList>(fo),
            "--options file "<<options<<" contains no OptionList");
        for (const gsOptionList::OptionListEntry & e : fo.getAllEntries())
            if (!knownKeys.count(e.label))
                gsInfo << "WARNING: option '" << e.label << "' from " << fdOpt.lastPath()
                       << " is not used by poisson_rh_schedule_example (typo?)\n";
        opt.update(fo, gsOptionList::addIfUnknown);
        gsInfo << "Method options from --options: " << fdOpt.lastPath() << "\n";
    }
    // CLI overrides, applied LAST so they win over both the built-in default
    // and any XML: --coarsen only ever turns Coarsen ON (there is no
    // --no-coarsen; omit the flag and rely on the XML/default for "off");
    // --target replaces Target outright (any sign), letting a sweep script
    // switch the band on/off and tune it per run without touching the XML.
    // Target is the error the run finishes at: the H step stops as soon as
    // the error is inside [Target/Band, Target*Band], so --target alone
    // fixes the endpoint and no second threshold has to be kept in step
    // with it.
    if (coarsenFlag) opt.setSwitch("Coarsen", true);
    if (!math::isnan(targetCli))
    {
        opt.setReal("Target", targetCli);
    }
    //! [Read input file]

    const real_t target=opt.getReal("Target"), band=opt.getReal("Band");
    const bool useBand=target>=0;
    if (useBand) GISMO_ENSURE(band>=1, "Band must be >= 1");
    const bool coarsen=opt.getSwitch("Coarsen");

    // The analysis space stays a plain (non-rational) tensor B-spline space
    // even when the geometry itself is rational (e.g. a NURBS sphere patch):
    // take the source basis of a rational one, or the basis itself if it is
    // already a plain tensor B-spline basis (2026-08-19, mirrors
    // l2projection_rh_schedule_example.cpp's identical extraction). This
    // does NOT approximate the geometry -- physical stays the exact
    // rational multipatch everywhere else (assembly, R-step monitor,
    // error evaluation); only the DISCRETE SOLUTION SPACE is non-rational.
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
            GISMO_ERROR("This example expects a tensor-product (B-spline or NURBS) basis");
    }
    basis.setDegree(degree);
    for (index_t i=0;i<initialRef;++i)
        basis.uniformRefine();

    gsTHBSplineBasis<2> thb(basis);

    // id=5 (optional): pin the sigma space from the file.  Absent -> build it
    // from -E/-R via the gsSquareDomain(basis, slide) constructor so the sigma
    // discretisation stays sweepable.
    // Sigma's mesh is a refinement LEVEL of the level-0 analysis basis, not a
    // free knot count.  makeIntegrationBasis (the paper's super mesh) unions the
    // sigma and analysis knot lines; a non-nested pair costs roughly their SUM of
    // elements, a nested pair only the finer of the two -- cf. Sec. "Numerical
    // integration scheme": "if either the space associated with sigma or the
    // analysis space is nested within the other, no additional knot lines are
    // introduced, leading to a more efficient quadrature process".  Every
    // analysis level is dyadic, so tying sigma to the same ladder makes a
    // non-nested sigma inexpressible.
    GISMO_ENSURE(sigmaRef>=0,"sigmaRef must be non-negative");
    std::unique_ptr<gsSquareDomain<real_t> > sigmaPtr;
    if (fd.hasId(5))
    {
        gsGeometry<real_t>::uPtr sgeo = fd.getId<gsGeometry<real_t> >(5);
        sigmaPtr.reset(new gsSquareDomain<real_t>(*sgeo));
        gsInfo << "Composition geometry from file (id=5): overrides -E/-R\n";
    }
    else
    {
        gsKnotVector<real_t> ks(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
        gsTensorBSplineBasis<2,real_t> b(ks, ks);
        sigmaPtr.reset(new gsSquareDomain<real_t>(b, opt.getSwitch("Slide")));
    }
    gsSquareDomain<real_t>& sigma = *sigmaPtr;
    if (fd.hasId(5))
    {
        // Both gsSquareDomain constructors already add "Slide" (defaulted
        // false). On the id=5 path we still call addSwitch explicitly to
        // overwrite that ctor default with the CLI value --
        // gsOptionList::addSwitch assigns into m_switches unconditionally.
        sigma.options().addSwitch("Slide","Boundary controls slide along the boundary",
                                   opt.getSwitch("Slide"));
        sigma.applyOptions();
    }
    const gsTensorBSplineBasis<2> * sb =
        dynamic_cast<const gsTensorBSplineBasis<2>*>(&sigma.domain().basis());
    GISMO_ENSURE(sb, "sigma must carry a tensor-product B-spline basis");
    const gsTensorBSplineBasis<2> & sbasis = *sb;
    gsInfo << "schedule=" << sched << " project=" << project << " H=" << useH << "\n";

    // Persistent assembler and two sharing evaluators (constructed FROM the
    // assembler, so all see the same expression tree) used by every S step.
    // ev is the solver evaluator (integration elements = union basis); pev is
    // the plot/error evaluator, whose integration elements are set to the
    // ANALYSIS basis in solve() (see the long comment there).
    gsExprAssembler<> A(1,1);
    A.setOptions(opt);
    gsExprEvaluator<> ev(A);
    gsExprEvaluator<> pev(A);
    // Mirror the resolved quadrature onto the evaluators -- and onto every
    // element-wise error sweep in particular: the default element-sweep
    // quadrature (quA=1,quB=1) was measured to UNDER-INTEGRATE the marking
    // error by ~39%, which silently mis-ranks elements for refinement. The
    // norms must be integrated exactly like the system they measure.
    // SameElement=0 is mandatory under the sigma composition -- a quadrature
    // element of the integration basis does not map into a single element of
    // the space.
    for (gsExprEvaluator<>* e : {&ev,&pev})
    {
        e->options().setReal("quA", opt.getReal("quA"));
        e->options().setInt ("quB", opt.getInt ("quB"));
        e->options().addSwitch("SameElement","",opt.askSwitch("SameElement",false));
    }
    gsInfo << "Quadrature: quA=" << opt.getReal("quA") << " quB=" << opt.getInt("quB")
           << " SameElement=" << opt.askSwitch("SameElement",false) << "\n";
    // Recommended safety net (task-07 review, note 1): under the sigma
    // composition a quadrature element of the integration basis does not map
    // into a single element of the analysis space, so SameElement=true is
    // always wrong here regardless of where it was set from (built-in
    // default / id=4 / --options).
    if (opt.askSwitch("SameElement",false))
        gsWarn << "WARNING: SameElement=1 is set, but this driver always "
                  "assembles/integrates on a composed map: quadrature points "
                  "will cross element boundaries and results will be WRONG. "
                  "Set SameElement=0 in id=4 of -f or in --options.\n";

    gsInfo << opt << "\n";
    gsFileData<> fdout;
    fdout.add(opt);
    fdout.save(output+"options.xml");

    std::unique_ptr<gsParaviewCollection> solcol;
    if (plot)
    {
        solcol.reset(new gsParaviewCollection(output+"solution",&pev));
        solcol->options().setSwitch("plotElements", true);
        solcol->options().setInt("plotElements.resolution", 4);
        solcol->options().setInt("numPoints", 1000);
        solcol->options().setInt("precision", 12);
    }

    std::ofstream csv(output+"convergence.csv");
    // Unbuffered: a sweep kills a run that overruns its time budget, and a
    // buffered stream loses every row it has not flushed -- a 2 h run then
    // yields a 0-byte file instead of the convergence history it did earn.
    csv << std::unitbuf;
    csv << "cycle,step,path,dofs,dofs_sigma,minErr,maxErr,err,pctBelowTol,minDetJsigma,time\n";
    csv << std::scientific << std::setprecision(8);

    gsMultiPatch<> lastSol; bool haveSolve=false;
    real_t lastMinErr=0, lastMaxErr=0, lastErr=0, lastPctBelowTol=0;
    // The per-element L2 error vector driving the H-step marker: always the
    // COMPOSED (C) path's, even when --project also computed a P solve (D9:
    // "err is the L2 error of the last S step, the composed path when both
    // are computed").
    std::vector<real_t> lastElErr;
    index_t solveStep=0;
    bool done=false;
    // D13: "Non-S rows repeat the last known error values (empty fields
    // before the first solve)". A literal 0.0 there is indistinguishable
    // from a perfect solve for any downstream harvester, so the four error
    // columns (minErr,maxErr,err,pctBelowTol) are ONLY ever the last known
    // values once haveSolve is true; before that they are emitted as three
    // empty fields.
    auto errFields = [&]() -> std::string
    {
        if (!haveSolve) return ",,,";
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(8)
            << lastMinErr << "," << lastMaxErr << "," << lastErr << ","
            << lastPctBelowTol;
        return oss.str();
    };
    for (index_t cycle=0;cycle<iterations && !done;++cycle)
    {
        gsInfo << "Cycle " << cycle << " [" << sched << "]\n";
        for (char op:sched)
        {
        if (done) break;
        // dofs/minDetJsigma for the CSV are read at EACH branch's own write
        // point below (not cached here), so minDetJsigma in particular is
        // RECOMPUTED on every row (D18c), never carried forward.
        const index_t dofs = useH?thb.size():basis.size();
        gsStopwatch stepTimer;

        if (op=='R')
        {
            // Generic lambda instantiated at compile time for each MonitorMode:
            // the mode is a non-type template parameter of gsAdaptiveParametrization.
            auto relocate=[&](auto mode) {
                gsStopwatch timer;
                std::unique_ptr<gsOptimizer<real_t> > optimizer = makeOptimizer(opt);
                const gsTensorBSplineBasis<2>& tb=useH?thb.tensorLevel(thb.maxLevel()):basis;
                // Choose the integration mesh: the hierarchical one built on
                // the analysis THB's actual element partition is cheaper than
                // the finest tensor grid EXCEPT on anisotropic level-0 meshes,
                // where L_sigma (a single isotropic level) makes it more
                // expensive -- probe both candidates' element counts and take
                // the smaller (useHier = nH < nT is a strict inequality, so a
                // tie keeps the tensor arm). Both arms are pure knot/box
                // arithmetic (no quadrature, no assembly), so probing both is
                // negligible against the R step's own MaxIterations
                // objective+gradient sweeps -- measured below.
                // Choose the integration mesh, and KEEP the winning basis: it is
                // handed to the constructor via the pass-through
                // (integrationBasisIsFinal) overload below instead of the RAW
                // thb/tb, which would otherwise make the constructor rebuild
                // (knot union + degree raise) the very basis just probed here.
                bool useHier = false;
                typename gsBasis<real_t>::uPtr ib;
                if (useH)
                {
                    gsStopwatch probeTimer;
                    typename gsBasis<real_t>::uPtr hIb =
                        gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::
                        makeIntegrationBasis<2>(thb,sbasis);
                    typename gsBasis<real_t>::uPtr tIb =
                        memory::make_unique(new gsTensorBSplineBasis<2>(
                            gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::
                            makeIntegrationBasis(tb,sbasis)));
                    const index_t nH = (index_t)hIb->numElements();
                    const index_t nT = (index_t)tIb->numElements();
                    useHier = (nH < nT);
                    gsInfo << "  R | integration mesh: " << (useHier?"hierarchical":"tensor")
                           << " | hierarchical elements " << nH
                           << " | tensor elements " << nT
                           << " | probe " << probeTimer.stop() << " s\n";
                    ib = useHier ? give(hIb) : give(tIb);
                }
                else
                {
                    // Tensor-only path: no hierarchical candidate to probe
                    // against, just build the union basis once.
                    ib = memory::make_unique(new gsTensorBSplineBasis<2>(
                        gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::
                        makeIntegrationBasis(tb,sbasis)));
                }
                // parametric=true: the monitor is the discrete solution in the
                // PARAMETRIC domain (the analysis basis lives on [0,1]^2, and
                // that is where lastSol is expressed), matching the ctor's
                // documented "composition function defined in the parametric
                // domain" semantics.
                gsAdaptiveParametrization<real_t,mode.value> rel(
                    sigma, physical.patch(0), &lastSol.patch(0), *ib, *optimizer,
                    /*parametric=*/true, integrationBasisIsFinal);
                rel.options().setReal("Smoothing", opt.getReal("Smoothing"));
                rel.options().setReal("Penalty",   opt.getReal("Penalty"));
                rel.options().setReal("quA", opt.getReal("quA"));
                rel.options().setInt ("quB", opt.getInt ("quB"));
                rel.solve();
                const real_t tR=timer.stop();

                gsInfo << "  R | mode "
                       << (mode.value==MonitorMode::GradientBased?"gradient":"value")
                       << " | dofs " << std::setw(6) << (useH?thb.size():basis.size())
                       << " | sigma controls " << std::setw(4) << sigma.nControls()
                       << " | min det J_sigma " << std::scientific << std::setprecision(3)
                       << rel.computeMinJacobian()
                       << " | min det J_sigma (certificate) " << std::scientific
                       << std::setprecision(3) << sigma.minDetJCoefficient();
                gsInfo << " | " << std::fixed << std::setprecision(2) << tR << " s\n"
                       << std::defaultfloat;
            };

            if (lastSol.nPatches()==0)
                gsInfo<<"  R | no previous solution, skipping relocation\n";
            else
            {
                const std::string mm = opt.askString("MonitorMode","value");
                if (mm=="gradient")
                    relocate(std::integral_constant<enum MonitorMode,
                                                    MonitorMode::GradientBased>());
                else if (mm=="value")
                    relocate(std::integral_constant<enum MonitorMode,
                                                    MonitorMode::ValueBased>());
                else
                    GISMO_ERROR("Unknown MonitorMode '"<<mm<<"' (value | gradient)");
            }
            const real_t tStep=stepTimer.stop();
            csv << cycle << ",R,-," << dofs << "," << sigma.nControls() << ","
                << errFields() << "," << sigma.minJacobian(7) << "," << tStep << "\n";
        }
        else if (op=='U')
        {
            gsStopwatch timer;
            if (useH) { thb.uniformRefine(); }
            else { basis.uniformRefine(); }
            const real_t tU=timer.stop();
            gsInfo << "  U | dofs " << std::setw(6) << (useH?thb.size():basis.size())
                   << " | " << std::fixed << std::setprecision(2) << tU << " s\n"
                   << std::defaultfloat;
            const real_t tStep=stepTimer.stop();
            csv << cycle << ",U,-," << (useH?thb.size():basis.size()) << ","
                << sigma.nControls() << "," << errFields() << ","
                << sigma.minJacobian(7) << "," << tStep << "\n";
        }
        else if (op=='S')
        {
            const gsTensorBSplineBasis<2>& tb=useH?thb.tensorLevel(thb.maxLevel()):basis;
            const gsBasis<>& active=useH?static_cast<const gsBasis<>&>(thb):
                                             static_cast<const gsBasis<>&>(basis);
            gsMultiBasis<> mb(active);
            gsTensorBSplineBasis<2> ibasis=
                gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>::
                makeIntegrationBasis(tb,sbasis);
            gsMultiBasis<> ib(ibasis);
            gsComposedGeometry<> composed(sigma,physical.patch(0));
            gsMultiPatch<> cmp; cmp.addPatch(composed);
            const index_t step=solveStep++;
            SolveResult cres = solve(A,ev,pev,cmp,mb,ib,f,ms,bc,'C',step,solcol.get());
            lastSol=cres.sol;
            haveSolve=true;

            const real_t minC=*std::min_element(cres.elErr.begin(),cres.elErr.end());
            const real_t maxC=*std::max_element(cres.elErr.begin(),cres.elErr.end());
            index_t nBelowC=0;
            // ORCHESTRATOR ADDENDUM (task 07's review): pctBelowTol is the
            // percentage of ELEMENTS below the EQUIDISTRIBUTION threshold
            // Target/sqrt(nElements), not the raw global Target (which would
            // make the column read ~100% at any realistic target).
            if (useBand)
            {
                const real_t thrEq = target / math::sqrt((real_t)cres.elErr.size());
                for (real_t e : cres.elErr) if (e <= thrEq) ++nBelowC;
            }
            const real_t pctC = useBand ? (100.0*nBelowC/(real_t)cres.elErr.size()) : 0.0;
            // This is the state the R/H/U branches and the H-step band rule
            // read: always the COMPOSED path's (see the member comment).
            lastMinErr=minC; lastMaxErr=maxC; lastErr=cres.l2err; lastPctBelowTol=pctC;
            lastElErr=cres.elErr;
            csv << cycle << ",S,C," << active.size() << "," << sigma.nControls() << ","
                << minC << "," << maxC << "," << cres.l2err << "," << pctC << ","
                << sigma.minJacobian(7) << "," << stepTimer.stop() << "\n";

            if (project)
            {
                gsStopwatch pTimer;
                gsMatrix<> coefs;
                gsL2Projection<real_t>::project(active,ibasis,composed,coefs);
                coefs.resize(active.size(), physical.patch(0).targetDim());
                gsGeometry<>::uPtr pg=active.makeGeometry(give(coefs));
                gsMultiPatch<> pmp; pmp.addPatch(*pg); gsMultiBasis<> pmb(active);
                SolveResult pres = solve(A,ev,pev,pmp,pmb,pmb,f,ms,bc,'P',step,solcol.get());
                // The P solve is reported in its own CSV row and updates
                // lastSol for the NEXT R step (matching the pre-existing
                // composed-vs-projected mechanism), but does NOT become the
                // driving state for the H-step band rule (see the member
                // comment on lastElErr): only its OWN row uses
                // minP/maxP/pres.l2err.
                lastSol=pres.sol;

                const real_t minP=*std::min_element(pres.elErr.begin(),pres.elErr.end());
                const real_t maxP=*std::max_element(pres.elErr.begin(),pres.elErr.end());
                index_t nBelowP=0;
                if (useBand)
                {
                    const real_t thrEq = target / math::sqrt((real_t)pres.elErr.size());
                    for (real_t e : pres.elErr) if (e <= thrEq) ++nBelowP;
                }
                const real_t pctP = useBand ? (100.0*nBelowP/(real_t)pres.elErr.size()) : 0.0;
                csv << cycle << ",S,P," << active.size() << "," << sigma.nControls() << ","
                    << minP << "," << maxP << "," << pres.l2err << "," << pctP << ","
                    << sigma.minJacobian(7) << "," << pTimer.stop() << "\n";
            }
        }
        else if (op=='H')
        {
            gsStopwatch timer;
            // ORCHESTRATOR ADDENDUM (task 07's review, point 2): a schedule
            // step with an unmet precondition SKIPS with a warning; it never
            // aborts -- uniformly at every position in the schedule,
            // including the first letter.
            if (!(useH && haveSolve))
            {
                gsInfo << "  H | " << (useH ? "no previous solve"
                                            : "'H' not in --schedule")
                       << ", skipping\n";
                const real_t tStep=stepTimer.stop();
                csv << cycle << ",H,-," << (useH?thb.size():basis.size()) << ","
                    << sigma.nControls() << "," << errFields() << ","
                    << sigma.minJacobian(7) << "," << tStep << "\n";
                continue;
            }
            const gsBasis<>& active=useH?static_cast<const gsBasis<>&>(thb):
                                             static_cast<const gsBasis<>&>(basis);

            // Band decision on the GLOBAL L2 error of the last S step (the
            // composed path when --project made both C and P available).
            // Without a band (Target < 0) this reduces to the historical
            // "always refine, coarsen iff Coarsen" behaviour, so the default
            // run is comparable with the pre-band baseline.
            //
            // 2026-08-20: the band is a TARGET ZONE the H step drives the
            // error INTO, from whichever side it starts on, and reaching it
            // ends the run:
            //
            //   err >  Target*Band   "above"  refine only, and keep going
            //   in the band          "band"   ONE combined refine+coarsen
            //                                 sweep, then stop
            //   err <= Target/Band   "below"  coarsen only, and keep going
            //
            // The previous rule refined inside the band as well and only
            // stopped at Target/Band^2, i.e. it had to overshoot the
            // requested Target by a factor Band^2 before it was allowed to
            // finish. On a case with an error floor inside the band that
            // never terminates: measured on the fitting driver, the rmse sat
            // at 5.87e-05 against a 1e-4 Target while the DoFs ran from 4779
            // to 27716 -- a 5.8x mesh for a 0.3% error gain. It also meant a
            // sweep arm could not land on its reference line's last point,
            // which is the whole purpose of passing that point as --target.
            //
            // Coarsening now follows the branch instead of running on every
            // H step (reversing the 2026-08-18 "unconditional" rule): above
            // the band there is nothing to give back, below it refinement is
            // not what is called for, and in the band both act once. The
            // per-element marking is unchanged (gsHElementMarker::markCrs,
            // relative to the CURRENT max element error).
            bool doRef=true, doCrs=coarsen;
            const char* branch="-";
            if (useBand)
            {
                if      (lastErr > target*band) { doRef=true;  doCrs=false;   branch="above"; }
                else if (lastErr > target/band) { doRef=true;  doCrs=coarsen; branch="band";
                                                  done=true; }
                else                            { doRef=false; doCrs=coarsen; branch="below";
                    // Below the band with coarsening switched off there is
                    // no action left that could raise the error back into
                    // it, so the loop would spin. Stop instead.
                    if (!coarsen) done=true; }
            }

            std::vector<index_t> boxes, crsBoxes;
            size_t nCrs=0;
            // lastElErr is cached from the last S step's per-element error
            // vector, indexed by the id() of THAT S step's analysis-basis
            // elements. If the mesh changed since then with no intervening S
            // (e.g. the second H of "SHH", or an H right after a bare U),
            // that vector no longer lines up with active's current elements
            // (gsHElementMarker::setErrors would reject the size mismatch);
            // rather than aborting mid-run, skip this H's marking and say so
            // -- the fix documented next to the "more H letters" sentence in
            // the file header.
            const bool haveValidElErr =
                lastElErr.size()==(size_t)active.numElements();
            if ((doRef || doCrs) && !haveValidElErr)
                gsInfo << "  H | mesh changed since the last solve, skipping "
                          "(insert an S before this H)\n";
            if ((doRef || doCrs) && haveValidElErr)
            {
                gsHElementMarker<2,real_t> marker(active);
                marker.options().setInt   ("RefineRule",   opt.getInt ("RefineRule"));
                marker.options().setReal  ("RefineParam",  opt.getReal("RefineParam"));
                marker.options().setInt   ("CoarsenRule",  opt.getInt ("CoarsenRule"));
                marker.options().setReal  ("CoarsenParam", opt.getReal("CoarsenParam"));
                marker.options().setInt   ("MaxLevel",     opt.getInt ("MaxLevel"));
                marker.options().setSwitch("Admissible",   opt.getSwitch("Admissible"));
                marker.options().setSwitch("Extension",    opt.getSwitch("Extension"));
                marker.setErrors(lastElErr);
                const gsHElementMarker<2,real_t>::HElementContainer markedRef =
                    doRef ? marker.markRef() : gsHElementMarker<2,real_t>::HElementContainer();
                if (doRef) boxes = marker.toRefBoxes(markedRef);
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

                // Coarsening: elements with the SMALLEST errors are un-refined,
                // which lets the mesh recover from h-refinement that a later R
                // step (or a better solution) makes obsolete. markCrs() gets
                // the CLOSED refined set (elements about to be refined, and
                // their siblings, must not be coarsened) -- task 02 changed
                // gsHElementMarker::markCrs to dispatch on "CoarsenRule"
                // directly, so no save/restore of "RefineRule" is needed here.
                if (doCrs)
                {
                    const gsHElementMarker<2,real_t>::HElementContainer markedCrs =
                        marker.markCrs(markedRef);
                    nCrs=markedCrs.size();
                    crsBoxes=marker.toCrsBoxes(markedCrs);
                }
            }

            // Both box lists are computed on the SAME (pre-update) mesh; the
            // refined and coarsened regions are disjoint by construction, so
            // applying them in sequence is safe.
            if (!boxes.empty())    thb.refineElements(boxes);
            if (!crsBoxes.empty()) thb.unrefineElements(crsBoxes);
            const real_t tH=timer.stop();
            if (done)
                gsInfo << "  H | err " << std::scientific << std::setprecision(3)
                       << lastErr << " (" << branch << ") | band reached, stopping\n"
                       << std::defaultfloat;
            gsInfo << "  H | marked boxes " << std::setw(4) << boxes.size();
            if (doCrs) gsInfo << " | coarsened " << std::setw(4) << nCrs;
            gsInfo  << " | elements " << thb.numElements()
                   << " | " << std::fixed << std::setprecision(2) << tH << " s\n"
                   << std::defaultfloat;
            const real_t tStep=stepTimer.stop();
            csv << cycle << ",H,-," << thb.size() << "," << sigma.nControls() << ","
                << errFields() << "," << sigma.minJacobian(7) << "," << tStep << "\n";
        }
        }
    }
    if (plot)
    {
        solcol->save();
        // The final composed geometry G = S o sigma, freshly built from the
        // final sigma/analysis state (D12) -- independent of whether the
        // last S step also solved on a --project'ed geometry.
        gsComposedGeometry<> finalComposed(sigma,physical.patch(0));
        gsWriteParaview(finalComposed, output+"geometry", 1000, true, true);
        gsWriteParaview(sigma.domain(), output+"sigma", 1000, true, true);
        const gsBasis<>& active=useH?static_cast<const gsBasis<>&>(thb):
                                         static_cast<const gsBasis<>&>(basis);
        gsMesh<real_t> msh(active);
        gsWriteParaview(msh, output+"mesh", false);
    }
    return 0;
}
