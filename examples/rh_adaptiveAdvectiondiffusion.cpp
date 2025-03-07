/** @file rh_adaptiveAdvectiondiffusion.cpp

    @brief Tutorial on how to use expression assembler to solve a advection-diffusion equation in rh-adaptive mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>  // Include the new class of r_refinement

using namespace std;
using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    bool plot             = false;
    index_t numRefine     = 3;// for local refinement:  0 means no local h-refinement
    index_t UnifRefine    = 4;// initial refinement: for MAE resolution take at least >=3 for Bejictive mapping 
    index_t DegElevate    = 2; // degree Elevation
    index_t NumArMarEl    = 1; // Number of ring of cells around marked elements
    index_t maxIter       = 30;
    double IntensityMAE   = 9.;
    real_t adaptRefParam  = 0.;     // ... adapt parameter.
    index_t FactRefPar    = 0;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    bool ErrorPrint       = true, export_b64 =false;
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-APsiere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "DegElevate",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", DegElevate);
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  UnifRefine );
    cmd.addInt( "l", "numRefine", "Number of local h-refinement loops",  numRefine );
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addInt( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addInt( "c", "NumArMarEl", "augement NumArMarEl with such quantity in local h-refinement loops",  NumArMarEl );
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addSwitch( "ErrorPrint", "print Error", ErrorPrint);
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    // .... one single patch
    gsMultiPatch<> Psi = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    Psi.degreeElevate(DegElevate);
    Psi.computeTopology();

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //..... Test 1 : POISSON EQUATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // // Define Stabilization method
    // auto Stabilizationtype = stabilizerCDR::none;
    // // convection coefficient
    // gsFunctionExpr<> coeff_conv("0.","0.",2);
    // // diffusion coefficient:
    // gsFunctionExpr<> coeff_diff("1.","0","0","1.",2);
    // // For a posterior error estimate
    // gsFunctionExpr<> coeff_diffMax("1.",2);
    // // reaction coefficient:
    // gsFunctionExpr<> coeff_reac("0",2);
    /* *********************** Test 1 *********************** */
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("1./(1.+exp((y - x  - 0.2)/0.01))", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1./(1.+exp((y - x  - 0.2)/0.01))",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("4.12230724487712e-5*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2 - 1.69934170211664e-13*exp(-200.0*x + 200.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**3",2);
    /* *********************** Test 2 *********************** */
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("1./(1.+exp(100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y)) ))", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1./(1.+exp(100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y)) ))",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("40000.0*x**2*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 - 80000.0*x**2*exp(200*x**2 + 200*(y - 0.5)**2 - 150.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**3 + 1.0*(75.0*pi**2*sin(pi*y) + 200)*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 + 1.0*(40000*(y - 0.375*pi*cos(pi*y) - 0.5)**2)*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 + 200.0*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 - 80000.0*(y - 0.375*pi*cos(pi*y) - 0.5)**2*exp(200*x**2 + 200*(y - 0.5)**2 - 150.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**3",2);


    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // //..... Test 2 ADVECTION DUFFFUSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // Define Stabilization method
    auto Stabilizationtype = stabilizerCDR::SUPG;
    // Define  Dirichlet boundary conditions
    gsFunctionExpr<> Dg("if( y <= 0.2*(1.-x), 1,0)", 2);
    // Manufactured solition
    gsFunctionExpr<> s("if( y <1., 1./(1.+exp((y - x  - 0.2)/0.0001))-1/(1.+exp((1.-x)/0.0001)),0)",2);
    // convection coefficient:
    gsFunctionExpr<> coeff_conv("cos(pi/4)","sin(pi/4)",2);
    // diffusion coefficient:
    gsFunctionExpr<> coeff_diff("0.000001","0","0","0.000001",2);
    // For a posterior error estimate
    gsFunctionExpr<> coeff_diffMax("0.000001",2);
    // reaction coefficient:
    gsFunctionExpr<> coeff_reac("0",2);
    // // Right-hand side function
    gsFunctionExpr<> SourceFunc("0.",2);
    //Manufactured density function 1./cosh(100. * ( -x - 0.2 + y ))
    gsFunctionExpr<> f("( 1./cosh( 10.*( -x+y -0.2 ) )**2 + 1/(1.+exp((0.95-x)/0.01)) )",2);
    
    gsInfo<<"The Initial domain is "<< Psi.detail() << "\n";


    //gsOptionList Aopt;

    //! [Refinement]
    gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    //dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    //gsInfo << dbasis.degree(0) << " degree  \n";

    gsInfo << "Patches: "<< Psi.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< oPsi_get_max_threads() <<"\n";
#endif

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<>  h1err(numRefine+1), l2err(numRefine+1);
    gsVector<int>  DoFPDE(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;

     //! [Refinement]
    Psi.uniformRefine();
    for (int r=0; r<=UnifRefine; ++r)
    {
        dbasis.uniformRefine();
        //Psi.uniformRefine();
        Psi.uniformRefine();
    }
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    
    // Set Dirichlet boundary conditions
    gsBoundaryConditions<> bc;
    bc.setGeoMap(Psi);
    // given by exact solution Dg on all boundaries:
    for ( gsMultiPatch<>::const_biterator
                bit = Psi.bBegin(); bit != Psi.bEnd(); ++bit)
    {
        bc.addCondition( *bit, condition_type::dirichlet, &Dg );
    }
    // --------------- define Pde ---------------
    timer.restart();
    //! [definePde]
    gsConvDiffRePde<real_t> cdrPde(Psi, bc, & coeff_diff,& coeff_conv, & coeff_reac, & SourceFunc);
    //! [definePde]
    //! [constructAssembler]
    // Construct assembler
    gsCDRAssembler<real_t> cdrAss( cdrPde, dbasis);
    // Set stabilization flag to 1 = SUPG
    cdrAss.options().setInt("Stabilization", Stabilizationtype);
    // CoPsiute Dirichlet values by L2-projection
    // Caution: Interpolation does not work for locally refined (T)HB-splines!
    cdrAss.options().setInt("DirichletValues",dirichlet::l2Projection);
    //! [constructAssembler]
    gsInfo<< "." <<std::flush;// Assemblying done
    ma_time += timer.stop();

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(cdrAss.multiBasis());

    // // Set the geometry optimal map
    geometryMap PP    = ev.getMap(Psi);  

    // Recover rhs for Poisson equation
    auto SFunc        = ev.getVariable(SourceFunc, PP);
    auto u_ex         = ev.getVariable(s, PP);
    auto frho         = ev.getVariable(f, PP);

    // Coeffs for advection-reaction diffusion equation
    auto coeff_convPP = ev.getVariable(coeff_conv, PP);
    auto coeff_diffPP = ev.getVariable(coeff_diffMax, PP);
    auto coeff_reacPP = ev.getVariable(coeff_reac, PP);
    // numerical solutionas Vector
    gsMatrix<> rsolVector;

    // Generate system matrix and load vector
    cdrAss.assemble();
    // Solution vector and solution variable
    timer.restart();
    // Solve the system
    rsolVector = gsSparseSolver<>::BiCGSTABILUT( cdrAss.matrix() ).solve( cdrAss.rhs() );

    slv_time += timer.stop();
    gsInfo<< "DoFs in PDEs " << cdrAss.numDofs() <<std::flush;
    gsField<> solField;        // Construct the solution as a scalar field
    solField = cdrAss.constructSolution(rsolVector);
    ev.setIntegrationElements(cdrAss.multiBasis());
    gsExprEvaluator<>::variable ru_sol = ev.getVariable(solField.fields());

    // Get the element-wise norms.
    //ev.integralElWise( ( coeff_diffPP * ilapl(ru_sol,PP) - igrad(ru_sol, PP)[0]*0.7071067811865476- igrad(ru_sol, PP)[1]*0.7071067811865476 - coeff_reacPP * ru_sol + SFunc).sqNorm() );
    ev.integralElWise( frho);
    auto elwise = ev.elementwise();
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 1-2 : CoPsiutes the density function
    ###         and the multipatch adaptove mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, Psi, DegElevate, maxIter, IntensityMAE);
    auto density = MAE.buildDensity(elwise, UnifRefine, 0);
    auto Psitp   = MAE.buildMultiPatch(density);
    if (true){
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsMultiPatch<> Psi;
    for(size_t i =0; i<Psitp.nPatches(); ++i)
        Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 5: local h-refinement
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsInfo << "Patches: "<< Psi.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";

    for (int r=0; r<=numRefine; ++r)
    {
        // Generate system matrix and load vector
        cdrAss.assemble();
        // Solution vector and solution variable
        timer.restart();
        // Solve the system
        rsolVector = gsSparseSolver<>::BiCGSTABILUT( cdrAss.matrix() ).solve( cdrAss.rhs() );

        slv_time += timer.stop();
        gsInfo<< "DoFs in PDEs " << cdrAss.numDofs() <<std::flush;
        DoFPDE[r] = cdrAss.numDofs();
        gsField<> solField;        // Construct the solution as a scalar field
        solField = cdrAss.constructSolution(rsolVector);

        ev.setIntegrationElements(cdrAss.multiBasis());

        gsExprEvaluator<>::variable ru_sol = ev.getVariable(solField.fields());

        // oPsi_set_dynamic(0);     // Explicitly disable dynamic teams
        // oPsi_set_num_threads(1); // Use these threads for later parallel regions

        //! [errorCoPsiutation]
        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(PP) ) );
        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(ru_sol,PP) ).sqNorm() * meas(PP) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error coPsiutations done

        if(r < numRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
        // --------------- error estimation/coPsiutation ---------------
        // for test :ev.integralElWise( ( ilapl(ru_sol, PP) + SFunc ).sqNorm()*meas(PP) );
        // Get the element-wise norms.
        ev.integralElWise( ( coeff_diffPP * ilapl(ru_sol,PP) - igrad(ru_sol, PP)[0]*0.7071067811865476- igrad(ru_sol, PP)[1]*0.7071067811865476 - coeff_reacPP * ru_sol + SFunc).sqNorm() );

        const std::vector<real_t> eltErrs  = ev.elementwise();
        //! [errorCoPsiutation]
        //gsInfo<< "\n ..size of "<< eltErrs.size() << " "<< rsolVector.size() <<std::flush<< "\n"; // Linear solving done

        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the coPsiuted local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked( eltErrs.size() );
        gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( cdrAss.multiBasis(), elMarked, NumArMarEl);
        gsRefineMarkedElements( dbasis, elMarked, NumArMarEl);
        gsRefineMarkedElements( Psi, elMarked, NumArMarEl);

       // Call repair interfaces to make sure that the new meshes
       // match along patch interfaces.
       cdrAss.multiBasis().repairInterfaces( Psi.interfaces() );
       //! [refreshAssembler]
       cdrAss.refresh();
       if(r%2==0 && r>0)
       NumArMarEl = NumArMarEl + FactRefPar;
        }
    if(plot && r == numRefine){

        // gsInfo<<"Storing paraview...\n";
        // // Write the coPsiuted solution to paraview files
        // gsWriteParaview<>( solField, "ParaviewOutput/solution", 1000, true);
        gsInfo<<"Making in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ru_sol,"numerical solution");
        // collection.addField(igrad(ru_sol,PP),"gradient_numerical solution");
        collection.addField(jac(PP).det(), "Jacobian function");
        // collection.addField(u_ex, "exact solution");
        collection.saveTimeStep();
        collection.save();
         }
    }
    //! [Solver loop]    

    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";
    

    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1_error= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";

    if (ErrorPrint && numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]
       //! [Export to Paraview]
       // Export the final solution
    if(plot){
        //------------------------------------
        gsInfo<<"Plotting in Paraview...\n";
        // Run paraview
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main