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

template <short_t comp>
struct slice
{
    slice(const std::vector<std::pair<index_t,index_t>> & indices)
    :
    m_indices(indices)
    {}

    std::vector<std::pair<index_t,index_t>> m_indices;

    index_t size() const { return m_indices.size(); }
    index_t operator[](index_t i) const { return std::get<comp>(m_indices[i]); }
};

int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefineG=0, numRefineD=0, numElevateG=0, numElevateD=0;
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
    gsSquareDomain<2,real_t> domain(*cbasis0_ptr);
    domain.options().addSwitch("Slide","",PARoptions.askSwitch("Slide",false));
    domain.applyOptions();

    // domain.perturb(1e-1);


/*
    PERFORM R-ADAPTIVITY
 */

    gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";
    gsOptim<real_t>::uPtr optimizer = gsOptim<real_t>::get(OPToptions.askString("Optimizer","LBFGS"));
    optimizer->options().update(OPToptions);

    gsInfo<<domain.domain().coefs()<<"\n";
    gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> * relocator;
    if (function.domainDim()==0)
        relocator = new gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>(domain,gpatch0,gbasis0,*optimizer,PARoptions.askSwitch("Parametric",false));
    else
        relocator = new gsAdaptiveParametrization<real_t,MonitorMode::ValueBased>(domain,gpatch0,function,gbasis0,*optimizer,PARoptions.askSwitch("Parametric",false));

    relocator->options().setReal("Penalty",PARoptions.askReal("Penalty",1e-2));
    relocator->options().setReal("Smoothing",PARoptions.askReal("Smoothing",1e-2));
    relocator->solve();

    gsInfo<<domain.domain().coefs()<<"\n";

    // //////////////////////////////////////////////////
    // // PLOTTING
    // //////////////////////////////////////////////////
    gsComposedGeometry<real_t> cspline(domain,gpatch0);
    gsMultiPatch<> mp;
    mp.addPatch(cspline);
    // mp.embed(3);
    // gsMultiBasis<> mb(mp);

    // gsExprEvaluator<> ev;
    // ev.setIntegrationElements(mb);
    // auto G = ev.getMap(mp);
    // // auto f = ev.getVariable(*fun);

    gsWriteParaview(mp,output+"cgeom",100000,true,true);
    gsWriteParaview(domain.domain(),output+"domain",100000,true,true);

    if (function.domainDim()!=0)
    {
        gsComposedFunction<real_t> cfun = (PARoptions.askSwitch("Parametric",true)) ? gsComposedFunction<real_t>(domain,function) : gsComposedFunction<real_t>(cspline,function);
        gsWriteParaview(mp,cfun,output+"cfun",100000);
    }



    // ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant");

    // delete fun;
    // delete optimizer;

    delete relocator;
    return 0;
}// end main
