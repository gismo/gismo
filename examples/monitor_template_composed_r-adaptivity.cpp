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
    index_t numRefine  = 2;
    index_t numElevate = 1;
    index_t maxIt = 100;
    real_t tol_g = 5e-5;
    real_t eps = 1e-2;
    bool slide = true;
    index_t testCase = 0;
    index_t opt = 2;
    index_t nSamplingPoints = 50;
    bool parametric = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal("g", "tolG", "relative tol", tol_g);
    cmd.addInt( "i", "maxIt", "max num iterations",  maxIt );
    cmd.addReal("", "eps", "eps",  eps );
    cmd.addSwitch("noslide", "Do not slide the boundaries",  slide );
    cmd.addInt( "t", "testCase", "Function to be used: 0: cosine waves, 1: spiral.",  testCase );
    cmd.addInt( "o", "opt", "Optimizer: 0: gsGradientDescent, 1: gsHLBFGS, 2: gsOptim::LBFGS.",  opt );
    cmd.addInt( "S", "nSamplingPoints", "Number of sampling points in each parametric direction",  nSamplingPoints );
    cmd.addSwitch("n", "parametric", "Define the function parametrically",  parametric );
    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_S"+util::to_string(nSamplingPoints)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();

    gsMultiPatch<> geom;
    // geom.addPatch(*gsNurbsCreator<>::BSplineSquare());
    geom.addPatch(*gsNurbsCreator<>::BSplineRectangle(0,0,2,1));
    geom.patch(0).uniformRefine(1,1,0);

    // Basis for the square domain
    gsKnotVector<> kv({0,0,1,1},1);
    gsTensorBSplineBasis<2> tbbasis(kv,kv);
    tbbasis.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbbasis.uniformRefine();

    tbbasis.uniformRefine(1,1,0);


    gsInfo<<"Mapper basis:\n"<<tbbasis<<"\n";

    gsSquareDomain<2,real_t> domain(tbbasis);
    domain.options().addSwitch("Slide","",slide);
    domain.applyOptions();
    // domain.perturb(1e-1);

    gsComposedGeometry<real_t> cspline(domain,geom.patch(0));
    gsFunction<> * fun;
    if      (testCase==0)
    {
        fun = new gsFunctionExpr<>("1 + 5 * exp( -50 * abs( (x-0.5)^2 + (y-0.5)^2 - 0.09 ) ) ",2);
    }
    else if (testCase==1)
    {
        std::string R     = "sqrt( (x-0.5)^2 + (y-0.5)^2 )";
        fun = new gsFunctionExpr<>("1 /(2 + cos( 8 * pi * "+R+"))",2);
    }
    else if (testCase==2)
    {
        std::string R     = "sqrt( (x-0.7)^2 + (y-0.5)^2 )";
        std::string Theta = "atan2((y-0.5),(x-0.7))";
        fun = new gsFunctionExpr<>("1 + 9/(1 + ( 10 * "+R+" * cos(" + Theta +" - 20 * "+R+"^2 ) )^2)",2);
    }
    else
    {
        GISMO_ERROR("Unknown test case");
    }


/*
    PERFORM R-ADAPTIVITY
 */

    gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";

    gsOptimizer<real_t> * optimizer;
    if      (opt==0) // gsGradientDescent
    {
        optimizer = new gsGradientDescent<real_t>;
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",2);
        optimizer->options().setReal("MinGradientLength",tol_g);

    }
    else if (opt==1) // gsHLBFGS
    {
        optimizer = new gsHLBFGS<real_t>;
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",2);
        optimizer->options().setReal("tolRelG",tol_g);
    }
    else if (opt==2) //gsOptim::LBFGS
    {
        optimizer = new gsOptim<real_t>::LBFGS;
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",1);
        optimizer->options().setReal("GradErrTol",tol_g);
    }
    else
    {
        GISMO_ERROR("Unknown optimizer");
    }

    gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,geom.patch(0),*fun,geom.basis(0),*optimizer,parametric);
    relocator.solve();

    //////////////////////////////////////////////////
    // PLOTTING
    //////////////////////////////////////////////////
    gsComposedFunction<real_t> cfun = (parametric) ? gsComposedFunction<real_t>(domain,*fun) : gsComposedFunction<real_t>(cspline,*fun);

    gsMultiPatch<> mp;
    mp.addPatch(cspline);
    // mp.embed(3);
    gsMultiBasis<> mb(mp);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(mp);
    // auto f = ev.getVariable(*fun);

    gsWriteParaview(mp,cfun,dirname+"cfun");
    gsWriteParaview(domain.domain(),dirname+"domain",1000,true,true);
    ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant");

    delete fun;
    delete optimizer;
    return 0;
}// end main
