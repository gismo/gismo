/** @file monitor_parameterization_newton_example.cpp

    @brief Composite spline relocation via Newton / Picard / HLBFGS.

    Same input format as monitor_parameterization_example.cpp; adds
    --solver newton|picard|hlbfgs|all to compare the exact-Hessian Newton
    solver and the SPD Picard surrogate (gsAdaptiveParametrizationNewton)
    against the quasi-Newton HLBFGS baseline (gsAdaptiveParametrization).

    See doc/derivation_hessian.md for the assembled operators.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveParametrizationNewton.h>

#include <gsHLBFGS/gsHLBFGS.h>

#include <gsUtils/gsStopwatch.h>
#include <iomanip>

using namespace gismo;

struct RunResult
{
    std::string name;
    index_t iters;
    real_t energy;
    real_t minJ;
    real_t time;
};

template<enum MonitorMode MODE>
RunResult runNewtonOrPicard(gsSquareDomain<real_t> & domain,
                            const gsGeometry<real_t> & geom,
                            const gsFunctionExpr<> & function,
                            const gsBasis<real_t> & ibasis,
                            bool parametric, bool picard,
                            real_t penalty, real_t smoothing,
                            const gsVector<real_t> & u0,
                            const std::string & logFile)
{
    domain.setControls(u0);

    typedef gsAdaptiveParametrizationNewton<real_t,MODE> Solver;
    memory::unique_ptr<Solver> solver;
    if (function.domainDim() == 0)
        solver.reset(new Solver(domain, geom, ibasis));
    else
        solver.reset(new Solver(domain, geom, function, ibasis, parametric));

    solver->options().setReal("Penalty",   penalty);
    solver->options().setReal("Smoothing", smoothing);
    solver->options().setString("LogFile", logFile);
    solver->options().setInt("MaxIter",1000);

    gsStopwatch watch;
    const index_t iters = picard ? solver->solvePicard() : solver->solveNewton();
    const real_t elapsed = watch.stop();

    RunResult res;
    res.name   = picard ? "Picard" : "Newton";
    res.iters  = iters;
    res.energy = solver->evalObj();
    res.minJ   = solver->computeMinJacobian();
    res.time   = elapsed;
    return res;
}

template<enum MonitorMode MODE>
RunResult runHLBFGS(gsSquareDomain<real_t> & domain,
                    const gsGeometry<real_t> & geom,
                    const gsFunctionExpr<> & function,
                    const gsBasis<real_t> & gbasis,      // raw geometry basis (relocator builds its own union)
                    const gsBasis<real_t> & ibasis,      // union basis (for the energy evaluation)
                    bool parametric,
                    real_t penalty, real_t smoothing,
                    const gsVector<real_t> & u0,
                    const gsOptionList & OPToptions)
{
    domain.setControls(u0);

    gsHLBFGS<real_t> optimizer;
    optimizer.options().update(OPToptions);

    typedef gsAdaptiveParametrization<real_t,MODE> Reloc;
    memory::unique_ptr<Reloc> relocator;
    if (function.domainDim() == 0)
        relocator.reset(new Reloc(domain, geom, gbasis, optimizer, parametric));
    else
        relocator.reset(new Reloc(domain, geom, function, gbasis, optimizer, parametric));

    relocator->options().setReal("Penalty",   penalty);
    relocator->options().setReal("Smoothing", smoothing);

    gsStopwatch watch;
    relocator->solve();
    const real_t elapsed = watch.stop();

    RunResult res;
    res.name   = "HLBFGS";
    res.iters  = optimizer.iterations();
    // Energy/minJ via a throwaway Newton object (shares evalObj/minJacobian)
    {
        typedef gsAdaptiveParametrizationNewton<real_t,MODE> Eval;
        memory::unique_ptr<Eval> ev;
        if (function.domainDim() == 0)
            ev.reset(new Eval(domain, geom, ibasis));
        else
            ev.reset(new Eval(domain, geom, function, ibasis, parametric));
        ev->options().setReal("Penalty",   penalty);
        ev->options().setReal("Smoothing", smoothing);
        res.energy = ev->evalObj();
        res.minJ   = ev->computeMinJacobian();
    }
    res.time = elapsed;
    return res;
}

int main(int arg, char *argv[])
{
    index_t numRefineG=0, numRefineD=1, numElevateG=0, numElevateD=1;
    index_t mode = MonitorMode::ValueBased;
    real_t penalty = -1;
    real_t smoothing = -1;
    std::string input = "parametrization/monitor_example_planar.xml";
    std::string output = "output";
    std::string solverName = "all";

    bool plot = false;
    index_t nSamples = 1000;

    gsCmdLine cmd("Composite spline relocation: Newton vs Picard vs HLBFGS.");
    cmd.addReal("P","penalty","Override penalty parameter",penalty);
    cmd.addReal("S","smoothing","Override smoothing parameter",smoothing);
    cmd.addInt("r","numRefG","Uniform h-refinements of the geometry",numRefineG);
    cmd.addInt("e","numElevG","Degree elevations of the geometry",numElevateG);
    cmd.addInt("R","numRefD","Uniform h-refinements of the domain",numRefineD);
    cmd.addInt("E","numElevD","Degree elevations of the domain",numElevateD);
    cmd.addInt("m","mode","Monitor mode: 0 ValueBased, 1 GradientBased",mode);
    cmd.addInt("n","nSamples","Number of samples for ParaView output",nSamples);
    cmd.addString("i","input","Input file",input);
    cmd.addString("o","output","Output directory (CSV iteration logs)",output);
    cmd.addString("","solver","Solver: newton|picard|hlbfgs|all",solverName);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution",plot);
    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }

    gsFileManager::mkdir(output);
    output += gsFileManager::getNativePathSeparator();

    // Read input (same labels as monitor_parameterization_example)
    gsFileData<> fd(input);
    gsMultiPatch<>   geom;
    gsMultiPatch<>   composition;
    gsFunctionExpr<> function;
    gsOptionList     OPToptions;
    gsOptionList     PARoptions;

    if (fd.hasLabel("geometry"))
        fd.getLabel("geometry", geom);
    else
    {
        GISMO_ENSURE(fd.getFirst(geom),
                     "Input file '" << input << "' has neither a 'geometry'-labelled "
                     "nor an unlabelled MultiPatch to use as the geometry.");
        gsInfo << "No 'geometry' label found in " << input
               << "; using the first MultiPatch in the file.\n";
    }
    if (fd.hasLabel("composition"))
        fd.getLabel("composition", composition);
    else
    {
        gsInfo << "No 'composition' label found in " << input
               << "; using the geometry basis as the sigma domain.\n";
        composition = geom;
    }
    if (fd.hasLabel("function"))
        fd.getLabel("function", function);
    if (fd.hasLabel("optimizer_options"))
        fd.getLabel("optimizer_options", OPToptions);
    if (fd.hasLabel("parametrization_options"))
        fd.getLabel("parametrization_options", PARoptions);

    GISMO_ENSURE(geom.nPatches()==1, "Only one patch is supported");

    for (size_t p=0; p!=geom.nPatches(); p++)
    {
        geom.patch(p).degreeElevate(numElevateG);
        for (index_t i = 0; i < numRefineG; i++)
            geom.patch(p).uniformRefine();
        composition.patch(p).degreeElevate(numElevateD);
        for (index_t i = 0; i < numRefineD; i++)
            composition.patch(p).uniformRefine();
    }

    const gsBasis<>    & cbasis0 = composition.basis(0);
    const gsGeometry<> & gpatch0 = geom.patch(0);
    const gsBasis<>    & gbasis0 = gpatch0.basis();

    // Union integration basis (geometry knots + composition knots, degree
    // p_G * p_sigma) — same construction as gsAdaptiveParametrization.
    // Required: the integration elements must not be coarser than the
    // sigma-basis elements (assembly assumes one active set per element).
    const gsTensorBSplineBasis<2> * comp_tb =
        dynamic_cast<const gsTensorBSplineBasis<2> *>(&cbasis0);
    GISMO_ENSURE(comp_tb, "Composition basis must be a 2D tensor B-spline basis");
    const gsTensorBSplineBasis<2> * geom_tb =
        dynamic_cast<const gsTensorBSplineBasis<2> *>(&gbasis0);
    if (!geom_tb)
    {
        const gsTensorNurbsBasis<2,real_t> * geom_nb =
            dynamic_cast<const gsTensorNurbsBasis<2,real_t> *>(&gbasis0);
        GISMO_ENSURE(geom_nb, "Geometry basis must be a 2D tensor B-spline/NURBS basis");
        geom_tb = &geom_nb->source();
    }
    gsTensorBSplineBasis<2> ibasis =
        gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>
            ::makeIntegrationBasis<2>(*geom_tb, *comp_tb);

    gsSquareDomain<real_t> domain(cbasis0);
    domain.options().addSwitch("Slide","",PARoptions.askSwitch("Slide",false));
    domain.applyOptions();

    const bool parametric = PARoptions.askSwitch("Parametric",false);
    const real_t pen   = (penalty==-1)   ? PARoptions.askReal("Penalty",  1e-2) : penalty;
    const real_t smth  = (smoothing==-1) ? PARoptions.askReal("Smoothing",1e-2) : smoothing;

    gsInfo << "Optimizer DoFs: " << domain.nControls() << "\n";
    if (domain.nControls() == 0)
    {
        gsWarn << "The composition has zero free controls: the sigma domain has "
                  "no interior control points, so there is nothing to relocate. "
                  "Increase -R/--numRefD or -E/--numElevD to give the "
                  "composition interior control points.\n";
        return 0;
    }
    gsVector<real_t> u0 = domain.getControls();

    std::vector<RunResult> results;
    const bool doNewton = (solverName=="newton" || solverName=="all");
    // const bool doPicard = (solverName=="picard" || solverName=="all");
    const bool doPicard = false;
    const bool doHLBFGS = (solverName=="hlbfgs" || solverName=="all");

#   define RUN_MODE(MODEVAL)                                                          \
    {                                                                                 \
        if (doNewton)                                                                 \
            results.push_back(runNewtonOrPicard<MODEVAL>(domain, gpatch0, function,   \
                ibasis, parametric, false, pen, smth, u0, output+"newton.csv"));      \
        if (doPicard)                                                                 \
            results.push_back(runNewtonOrPicard<MODEVAL>(domain, gpatch0, function,   \
                ibasis, parametric, true, pen, smth, u0, output+"picard.csv"));       \
        if (doHLBFGS)                                                                 \
            results.push_back(runHLBFGS<MODEVAL>(domain, gpatch0, function, gbasis0,  \
                ibasis, parametric, pen, smth, u0, OPToptions));                      \
    }

    if (mode == MonitorMode::ValueBased)
        RUN_MODE(MonitorMode::ValueBased)
    else
        RUN_MODE(MonitorMode::GradientBased)
#   undef RUN_MODE

    // Comparison table
    gsInfo << "\n" << std::string(70,'=') << "\n";
    gsInfo << std::setw(10) << "Solver" << std::setw(8) << "iters"
           << std::setw(18) << "energy" << std::setw(14) << "min det Js"
           << std::setw(12) << "time [s]" << "\n";
    gsInfo << std::string(70,'-') << "\n";
    for (const auto & r : results)
        gsInfo << std::setw(10) << r.name << std::setw(8) << r.iters
               << std::setw(18) << r.energy << std::setw(14) << r.minJ
               << std::setw(12) << r.time << "\n";
    gsInfo << std::string(70,'=') << "\n";

    // //////////////////////////////////////////////////
    // // PLOTTING
    // //////////////////////////////////////////////////
    if (plot)
    {
        gsComposedGeometry<real_t> cspline(domain.domain(),gpatch0);
        gsMultiPatch<> mp;
        mp.addPatch(cspline);
        // mp.embed(3);
        gsMultiBasis<> mb(mp);

        gsExprEvaluator<> ev;
        ev.options().setInt("plot.npts",nSamples);
        // ev.options().setSwitch("SameElement",false);
        ev.setIntegrationElements(mb);
        auto Gnew = ev.getMap(mp);
        auto Gold = ev.getMap(gpatch0);

        gsWriteParaview(mp,output+"cgeom",nSamples,true,true);
        gsWriteParaview(domain.domain(),output+"domain",nSamples,true,true);
        
        if (function.domainDim()!=0)
        {
            gsComposedFunction<real_t> cfun = (PARoptions.askSwitch("Parametric",true)) ? gsComposedFunction<real_t>(domain,function) : gsComposedFunction<real_t>(cspline,function);
            gsWriteParaview(mp,cfun,output+"cfun",nSamples);
            gsField<> fun(mp,function);
            gsWriteParaview(fun,output+"fun",nSamples);
        }

        // Export jacobian determinants
        if (mp.domainDim()==mp.targetDim())
        {
            auto detGnew = jac(Gnew).det();
            auto detGold = jac(Gold).det();
            ev.writeParaview(detGnew,Gnew,output+"jacobian_determinant_new");
            ev.writeParaview(detGold,Gold,output+"jacobian_determinant_old");
        }
        else if (mp.domainDim()<mp.targetDim())
        {
            auto fform_new = jac(Gnew).tr()*jac(Gnew);
            auto detGnew = pow(fform_new.det().val(),0.5); //jacobian determinant for a surface, i.e. the measure
            auto fform_old = jac(Gold).tr()*jac(Gold);
            auto detGold = pow(fform_old.det().val(),0.5); //jacobian determinant for a surface, i.e. the measure

            ev.writeParaview(detGnew,Gnew,output+"jacobian_determinant_new");
            ev.writeParaview(detGold,Gold,output+"jacobian_determinant_old");
        }
    }

    return 0;
}
