/** @file monitor_template_composed_r-adaptivity.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

#include <gsHLBFGS/gsHLBFGS.h>
#include <gsOptimizer/gsGradientDescent.h>
#include <gsOptim/gsOptim.h>

int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool gradient=false;
    index_t numRefineG=0, numRefineD=0, numElevateG=0, numElevateD=0;
    index_t nIsolines = 20;
    index_t nSamples = 10000;
    std::string input = "domain2d/lake.xml";
    std::string output = "output";

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt("r","numRefG","Number of Uniform h-refinement loops for the geometry",numRefineG);
    cmd.addInt("e","numElevG","Number of degree elevation steps to perform for the geometry",numElevateG);
    cmd.addInt("R","numRefD","Number of Uniform h-refinement loops for the domain",numRefineD);
    cmd.addInt("E","numElevD","Number of degree elevation steps to perform for the domain",numElevateD);
    cmd.addString("i", "input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("gradient", "Use gradient-based monitor", gradient);
    cmd.addInt("n","plotLines","number of isolines to export",nIsolines);
    cmd.addInt("s","plotPoints","number of sampling points for plot",nSamples);
    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // Make output directory
    gsFileManager::mkdir(output);
    output += gsFileManager::getNativePathSeparator();

    // Read input
    gsFileData<> fd(input);
    gsMultiPatch<>      geom;
    gsMultiPatch<>      composition;
    gsFunctionExpr<>    function;
    gsOptionList        OPToptions;
    gsOptionList        PARoptions;

    // Get the geometry
    fd.getLabel("geometry",geom);
    if (geom.nInterfaces()==0)
        geom.computeTopology();
    gsDebugVar(geom.nInterfaces());
    // Get the composition
    fd.getLabel("composition",composition);
    // Get the function
    if (fd.hasLabel("function"))
        fd.getLabel("function",function);
    gsInfo<<"Function:\n"<<function<<"\n";
    // Get the optimizer options
    fd.getLabel("optimizer_options",OPToptions);
    gsInfo<<OPToptions<<"\n";
    // Get the parametric options
    fd.getLabel("parametrization_options",PARoptions);
    gsInfo<<PARoptions<<"\n";

    GISMO_ENSURE(geom.nPatches()==composition.nPatches(),"The number of patches in the geometry and the composition must be the same");
    GISMO_ENSURE((index_t)geom.nPatches()==function.nPieces(),"The number of patches in the geometry and the function must be the same");
    GISMO_ENSURE(geom.nPatches()==1,"Only one patch is supported");

    // Refine and elevate the geometry and composition if needed
    for (size_t p=0; p!=geom.nPatches(); p++)
    {
        geom.patch(p).degreeElevate(numElevateG);
        for (index_t i = 0; i < numRefineG; i++)
            geom.patch(p).uniformRefine();

        composition.patch(p).degreeElevate(numElevateD);
        for (index_t i = 0; i < numRefineD; i++)
            composition.patch(p).uniformRefine();
    }

    const gsBasis<> & cbasis0 = composition.basis(0);
    const gsGeometry<> & gpatch0 = geom.patch(0);
    const gsBasis<> & gbasis0 = gpatch0.basis();

    GISMO_ENSURE(cbasis0.domainDim()==2,"The composition basis must be 2D (for now)");
    const gsTensorBSplineBasis<2,real_t> * cbasis0_ptr = dynamic_cast<const gsTensorBSplineBasis<2,real_t> *>(&cbasis0);
    GISMO_ENSURE(cbasis0_ptr,"The composition basis must be a tensor B-spline basis");

    gsInfo<<"Mapper basis:\n"<<*cbasis0_ptr<<"\n";
    gsSquareDomain<real_t> domain(*cbasis0_ptr);
    for (auto ifc = geom.iBegin(); ifc != geom.iEnd(); ++ifc)
        domain.addInterface(*ifc);

    domain.options().addSwitch("Slide","",PARoptions.askSwitch("Slide",false));
    domain.applyOptions();

    // domain.perturb(1e-1);


/*
    PERFORM R-ADAPTIVITY
 */

    gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";
    gsOptimizer<real_t> * optimizer;
    std::string OPTstring = OPToptions.askString("Optimizer","LBFGS");
    if (OPTstring=="HLBFGS")
        optimizer = new gsHLBFGS<real_t>();
    else
        optimizer = gsOptim<real_t>::get(OPToptions.askString("Optimizer","LBFGS")).release();

    optimizer->options().update(OPToptions);

    gsInfo<<domain.domain().coefs()<<"\n";
    // gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> * relocator;
    gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> * relocator;

    if (function.domainDim()==0)
        relocator = new gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>(domain,gpatch0,gbasis0,*optimizer,PARoptions.askSwitch("Parametric",false));
    else
        relocator = new gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>(domain,gpatch0,function,gbasis0,*optimizer,PARoptions.askSwitch("Parametric",false));

    gsDebugVar(PARoptions.askReal("Penalty",1e-2));
    relocator->options().setReal("Penalty",PARoptions.askReal("Penalty",1e-2));
    relocator->options().setReal("Smoothing",PARoptions.askReal("Smoothing",1e-2));
    gsDebugVar(relocator->options().getReal("Penalty"));
    relocator->solve();

    gsInfo<<domain.domain().coefs()<<"\n";

    // //////////////////////////////////////////////////
    // // PLOTTING
    // //////////////////////////////////////////////////
    gsComposedGeometry<real_t> cspline(domain.domain(),gpatch0);
    gsMultiPatch<> mp;
    mp.addPatch(cspline);
    // mp.embed(3);
    gsMultiBasis<> mb(mp);

    gsExprEvaluator<> ev;
    ev.options().setInt("plot.npts",nSamples);
    ev.options().setSwitch("SameElement",false);
    ev.setIntegrationElements(mb);
    auto Gnew = ev.getMap(mp);
    auto Gold = ev.getMap(gpatch0);

    gsWriteParaview(mp,output+"cgeom",1000,true,true);
    gsWriteParaview(domain.domain(),output+"domain",1000,true,true);

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

    // Plot iso lines
    gsKnotVector<> kv(0,1,nIsolines-2,1,1,0);
    gsTensorBSplineBasis<2> meshBasis(kv,kv);
    gsMesh<> meshnew(meshBasis,2);
    gsMesh<> meshold(meshnew);
    cspline.evaluateMesh(meshnew);
    gsWriteParaview(meshnew,output+"isolines_new");
    gpatch0.evaluateMesh(meshold);
    gsWriteParaview(meshold,output+"isolines_old");

    // Export geometry to file
    gsFileData<> fdout;
    fdout.add(cspline);
    fdout.save(output+"composedGeometry.xml");


    // delete fun;
    delete optimizer;

    delete relocator;
    return 0;
}// end main
