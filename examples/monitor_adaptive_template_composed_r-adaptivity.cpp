/** @file monitor_adaptive_template_composed_r-adaptivity.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsCore/gsComposedFunction.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>


int main(int arg, char *argv[])
{
    //! [Parse command line]
    
    index_t numURef  = 3;
    index_t numRefine  = 2;
    index_t numElevate = 1;
    index_t numRefineI = 5;
    index_t numElevateI = 3;
    index_t maxIt = 1000;
    real_t tol_g = 5e-5;
    real_t eps = 1e-4;
    bool slide = true;
    index_t testCase = 0;
    index_t opt = 2;
    index_t numAdaptiveLoops = 1;
    index_t rule = 1;
    real_t perc_ref = 0.1;
    bool parametric = false;
    
    

    gsCmdLine cmd("hr adaptivity with monitor function and composite maps, for a given error field.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis space and the sigma mapper", numElevate );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis space and the sigma mapper",  numRefine );
    
    cmd.addInt( "E", "elevIntegral","Number of degree elevation steps to perform for the integration", numElevateI );
    cmd.addInt( "R", "refIntegral", "Number of Uniform h-refinement loops for the integration",  numRefineI );
    
    
    cmd.addReal("g", "tolG", "relative tol", tol_g);
    cmd.addInt( "i", "maxIt", "max num iterations",  maxIt );
    cmd.addReal("", "eps", "eps",  eps );
    cmd.addSwitch("noslide", "Do not slide the boundaries",  slide );
    cmd.addInt( "t", "testCase", "Function to be used: 0: cosine waves, 1: spiral.",  testCase );
    cmd.addInt( "o", "opt", "Optimizer: 0: gsGradientDescent, 1: gsHLBFGS, 2: gsOptim::LBFGS.",  opt );
    cmd.addInt( "j", "ell", "Maximum number of iterations for the adaptive loop.",  numAdaptiveLoops );
    cmd.addInt( "w", "rule", "Refinement rule: 1 - GARU; 2 - PUKA; 3 - BULK.",  rule );
    cmd.addReal("p", "refp", "Refine percentage",  perc_ref );
    cmd.addSwitch("n", "parametric", "Define the function parametrically",  parametric );
    

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_S"+util::to_string(nSamplingPoints)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateI)+"_R"+util::to_string(numRefineI)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();


    // Geometry (square domain, represented by thb.)
    gsTensorBSpline<2,real_t> nurbs = *gsNurbsCreator<>::BSplineSquare(1,0,0);
    nurbs.uniformRefine( (1<<numURef)-1 );
    gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&nurbs);
    gsTHBSpline<2,real_t> thb = gsTHBSpline<2,real_t>(*geo);
    
    gsMultiPatch<> geom;
    geom.addPatch(thb);
    //geom.addPatch(*geo);
    gsWriteParaview(thb, dirname+"init_geom", 10000, true, true);



    // Basis for the sigma mapper
    gsKnotVector<> kv({0,0,1,1},1);
    gsTensorBSplineBasis<2> tbbasis(kv,kv);
    tbbasis.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbbasis.uniformRefine();
    
    gsInfo<<"Mapper basis:\n"<<tbbasis<<"\n";

    // // sigma mapper
    // gsSquareDomain<2,real_t> domain(tbbasis);
    // domain.options().addSwitch("Slide","",slide);
    // domain.applyOptions();


    // error field to test
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
        PERFORM H-ADAPTIVITY
    */
    //! [beginRefLoop]
    for( int refLoop = 0; refLoop < numAdaptiveLoops; refLoop++)
    {
            gsInfo << "#loop = "<< refLoop + 1<< "\n";
        /* 
            PERFORM R-ADAPTIVITY
        */

        // sigma mapper
        gsSquareDomain<2,real_t> domain(tbbasis);
        domain.options().addSwitch("Slide","",slide);
        domain.applyOptions();


        gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";

        gsWriteParaview(geom,dirname+"geom_in"+std::to_string(refLoop),1,true);
        
        gsMultiPatch<> mp;
        gsComposedGeometry<real_t> cspline(domain,geom.patch(0));
        mp.addPatch(cspline);
        gsWriteParaview(geom,dirname+"cspline_in"+std::to_string(refLoop),1,true);


        //gsOptMesh<> optMesh(domain,geom.patch(0),*fun,eps);
        //gsVector<> controls(domain.nControls());
        
        //optMesh.options().setInt("nRefine",numRefineI);
        //optMesh.options().setInt("nElevate",numElevateI);
        //for (size_t k=0; k!=domain.nControls(); k++)
        //    controls[k] = domain.control(k);


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

        gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> relocator(domain,
                                                                            geom.patch(0),
                                                                            *fun,
                                                                            geom.basis(0),
                                                                            *optimizer,
                                                                            parametric);
        relocator.solve();

        //optimizer->solve(controls);
        //gsVector<> optSol = optimizer->currentDesign();
        
        //for (size_t k=0; k!=optSol.rows(); k++)
        //    domain.control(k) = optSol[k];

        // gsInfo << "Domain after optimization:\n" << domain.domain() << "\n";
        // gsInfo << "Domain #coefficients after optimization:\n" << domain.domain().coefs().rows() << "\n";
        // gsInfo << "Domain #controls after optimization:\n" << domain.nControls() << "\n";

        
        // mp.embed(3);


        //////////////////////////////////////////////////
        // PLOTTING
        //////////////////////////////////////////////////
        gsComposedFunction<real_t> cfun = (parametric) ? gsComposedFunction<real_t>(domain,*fun) : gsComposedFunction<real_t>(cspline,*fun);

        gsMultiPatch<> mmmp;
        mmmp.addPatch(cspline);
        // mp.embed(3);
        gsMultiBasis<> mmmb(mmmp);

        gsExprEvaluator<> mmmev;
        mmmev.setIntegrationElements(mmmb);
        auto mmmG = mmmev.getMap(mp);
        // auto f = ev.getVariable(*fun);

        gsWriteParaview(mmmp,cfun,dirname+"mmmcfun");
        

        // MIO
        gsWriteParaview(cspline,dirname+"cspline_r"+std::to_string(refLoop),1000,true,true); // plots in xi, eta because it uses the support of the basis of cspline (thus the support of the cbasis, thus the support of the original basis)
        //gsWriteParaview(cspline.basis(),dirname+"cbasis_init"+std::to_string(refLoop),1000);

        gsMultiBasis<> mb(geom);

        gsInfo << "mmmb gsMultiBasis:\n" << mmmb << "\n";
        gsInfo << "  mb gsMultiBasis:\n" <<   mb << "\n";

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(mb);
        auto G = ev.getMap(mp); // geometry map
        auto f  = ev.getVariable(*fun,G);

        ev.writeParaview(f,G,dirname+"fun");
        gsWriteParaview(domain.domain(),dirname+"domain"+std::to_string(refLoop),1000,true,true);
        ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant"+std::to_string(refLoop));
        

        // Get the element-wise errors.
        //ev.integralElWise( f*meas(G) );
        ev.maxElWise(f);

        
        const std::vector<real_t> & elErrs = ev.elementwise();
        gsInfo << elErrs.size() << " element-wise errors\n";

        // to check
        gsElementErrorPlotter<real_t> err_eh(thb.basis(),elErrs);
        // gsComposedFunction<real_t> cerr_eh(domain,err_eh);
        // gsWriteParaview<>(err_eh, cspline.support(), dirname+"error_elem_ref", 10000); 
        gsWriteParaview<>(cspline,err_eh, dirname+"error_elem_ref"+std::to_string(refLoop), 1000); 

        gsAdaptiveMeshing<real_t> mesher(geom);
        mesher.options().setInt("RefineRule",rule);
        mesher.options().setInt("CoarsenRule",rule);
        mesher.options().setSwitch("Admissible",true);
        mesher.options().setReal("RefineParam",perc_ref);
        
        
        mesher.getOptions();
        gsHBoxContainer<2,real_t> refine;
        mesher.markRef_into(elErrs,refine);

        // gsInfo<<"Cells marked for refinement:\n";
        // gsInfo<<refine<<"\n";
        gsWriteParaview(refine,dirname+"marked4ref_"+std::to_string(refLoop));


        mesher.refine(refine);
        gsWriteParaview(geom,dirname+"geom_ref"+std::to_string(refLoop),1,true);

        gsComposedGeometry<real_t> refcspline(domain,geom.patch(0));
        gsWriteParaview(refcspline,dirname+"cspline_ref"+std::to_string(refLoop),1000,true,true);

        delete optimizer;

    }//! [end H Loop]

    delete fun;
    return 0;
}// end main
