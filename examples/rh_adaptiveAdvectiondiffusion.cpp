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
    //! [Parse command line]
    bool plot             = false;
    index_t numRefine     = 4;
    index_t numLRefine    = 3;
    index_t numElevate    = 0;
    index_t maxIter       = 50;
    index_t NumArMarEl    = 0; // Number of ring of cells around marked elements
    double IntensityMAE   = 6.;
    bool errorsave        = false;
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy... 
    index_t adaptRefCrit  = 2;  // 1: GARU, 2: PUCA, 3: BULK, 4: PBULK
    real_t  adaptRefParam = 0.; // ... adapt parameter.
    index_t FactRefPar    = 0;  // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    // Specify the file path
    std::string fn("pde/infinit_plate.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addInt( "c", "NumArMarEl", "augement NumArMarEl with such quantity in local h-refinement loops",  NumArMarEl );
    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addInt( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addInt( "r", "adaptRefCrit", "Adaptive refinement criterion [1:GARU,2:PUCA,3:BULK,4:PBULK]",  adaptRefCrit );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);

    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("errorsave", "Create a file in ... and save errors", errorsave);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //..... Test 1 : POISSON EQUATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    //     // Load the file
    // gsFileData<> fd(fn);
    // gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // // Create a gsMultipatch and add the loaded geometry
    // gsMultiPatch<> Psi;
    // fd.getId(1,Psi);
    // // Elevate and p-refine the basis to order p + numElevate
    // // where p is the highest degree in the bases
    // Psi.degreeElevate(numElevate);
    // Psi.computeTopology();
    // //...
    // // Define Stabilization method
    // auto Stabilizationtype = stabilizerCDR::SUPG;
    // // convection coefficient
    // gsFunctionExpr<> coeff_conv("1/(1.+exp((x +y  - 0.)/(0.001*2)))","1/(1.+exp((x +y  - 0.)/(0.001*2)))",2);
    // // diffusion coefficient:
    // gsFunctionExpr<> coeff_diff("0.001","0","0","0.001",2);
    // // For a posterior error estimate
    // gsFunctionExpr<> coeff_diffMax("0.001",2);
    // // reaction coefficient:
    // gsFunctionExpr<> coeff_reac("1./(0.001*2)*( 1.- 1/(1.+exp((x +y  - 0.)/(0.001*2))) )",2);
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("1/(1.+exp((x +y  - 0.)/(0.001*2)))", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1/(1.+exp((x +y  - 0.)/(0.001*2)))",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("0.",2);
    // // analytic density function
    // gsFunctionExpr<> f("1./cosh( 25.*( x+y -0. ) )",2);

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // //..... Test 2 ADVECTION DUFFFUSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // .... one single patch
    gsMultiPatch<> Psi = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    Psi.degreeElevate(numElevate);
    Psi.computeTopology();
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

    gsVector<>  h1err(numLRefine+1), l2err(numLRefine+1);
    gsVector<int>  DoFPDE(numLRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;

     //! [Refinement]
    Psi.uniformRefine();
    for (int r=0; r<= numRefine; ++r)
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
    //ev.integralElWise( ( coeff_diffPP * ilapl(ru_sol,PP) - igrad(ru_sol, PP) * coeff_convPP- coeff_reacPP * ru_sol + SFunc).sqNorm() );
    ev.integralElWise( (igrad(ru_sol, PP)).sqNorm() );
    // ev.integralElWise( frho);
    auto elwise = ev.elementwise();
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 1-2 : Computes the density function
    ###         and the multipatch adaptove mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(dbasis, Psi, numElevate, maxIter, IntensityMAE);
    // auto density = MAE.buildDensity(elwise, 0.005, 0);
    auto density = MAE.buildStrategyDensity(elwise, 0.95);

    // auto density = MAE.buildAnalyticDensity(f);
    auto Psitp   = MAE.buildMultiPatch(density, false); // false means we work on the computational domain
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


   // --------------- add bonudary conditions ---------------
   //! [Boundary conditions]
   gsBoundaryConditions<> bcInfo;

   // For simplicity, set Dirichlet boundary conditions
   // given by exact solution g on all boundaries:
   for ( gsMultiPatch<>::const_biterator
            bit = Psi.bBegin(); bit != Psi.bEnd(); ++bit)
   {
       bcInfo.addCondition( *bit, condition_type::dirichlet, &Dg ,0, false);
   }
   //! [Boundary conditions]

   // --------------- set up basis ---------------

   //! [GetBasisFromTHB]
   // Copy basis from the geometry
   gsMultiBasis<> bases( Psi, true );//true: poly-splines (not NURBS)
   //! [GetBasisFromTHB]

   // --------------- define Pde ---------------
   //! [definePde]
   gsConvDiffRePde<real_t> cdrPde(Psi, bcInfo, & coeff_diff,& coeff_conv, & coeff_reac, & SourceFunc);
   //! [definePde]

   // --------------- set up adaptive refinement loop ---------------

   //! [constructAssembler]
   // Construct assembler
   gsCDRAssembler<real_t> cdrAss( cdrPde, bases);
   // Set stabilization flag to 1 = SUPG
   cdrAss.options().setInt("Stabilization", Stabilizationtype);
   // Compute Dirichlet values by L2-projection
   // Caution: Interpolation does not work for locally refined (T)HB-splines!
   cdrAss.options().setInt("DirichletValues",dirichlet::l2Projection);
   //! [constructAssembler]

   // --------------- adaptive refinement loop ---------------
   //! [beginRefLoop]
   for( int refLoop = 0; refLoop <= numLRefine; refLoop++)
   {
   //! [beginRefLoop]
       gsInfo << "====== Loop " << refLoop << " of "
              <<numLRefine<< " ======" << "\n";
       // --------------- solving ---------------

       //! [solverPart]
       // Generate system matrix and load vector
       cdrAss.assemble();

       // Solve the system
       gsMatrix<real_t> solVector =
           gsSparseSolver<>::BiCGSTABILUT( cdrAss.matrix() ).solve( cdrAss.rhs() );

       // Construct the solution as a scalar field
       gsField<> solField;
       solField = cdrAss.constructSolution(solVector);
       //! [solverPart]

       // --------------- error estimation/computation ---------------

       //! [errorComputation]
       // Compute the H1-seminorm of the computed solution
       // ( which is, at least, equivalent to the energy norm in this example )
       // using the known exact solution.
       gsExprEvaluator<> ev;
       ev.setIntegrationElements(cdrAss.multiBasis());
       gsExprEvaluator<>::geometryMap Gm = ev.getMap(Psi);
       gsExprEvaluator<>::variable is = ev.getVariable(solField.fields());
       auto ms = ev.getVariable(Dg, Gm);

       // Get the element-wise norms.
        // Recover rhs for Poisson equation
        auto SFunc        = ev.getVariable(SourceFunc, Gm);
        // Coeffs for advection-reaction diffusion equation
        auto coeff_convGm = ev.getVariable(coeff_conv, Gm);
        auto coeff_diffGm = ev.getVariable(coeff_diffMax, Gm);
        auto coeff_reacGm = ev.getVariable(coeff_reac, Gm);

        // Recover manufactured solution for Poisson equation
        auto u_ex         = ev.getVariable(s, Gm);
        l2err[refLoop]    = math::sqrt(ev.integralElWise( (u_ex - is).sqNorm()*meas(Gm) ) );
        h1err[refLoop]    = math::sqrt(ev.integralElWise( (igrad(u_ex,Gm) - igrad(is,Gm)).sqNorm()*meas(Gm) ) );
        DoFPDE[refLoop]   = cdrAss.numDofs();

       // --------------- adaptive refinement ---------------
       ev.integralElWise( ( coeff_diffGm * ilapl(is,Gm) - igrad(is, Gm)*coeff_convGm - coeff_reacGm * is + SFunc).sqNorm()*meas(Gm) );
       //ev.integralElWise( ( igrad(is,Gm)).sqNorm() );
       const std::vector<real_t> eltErrs  = ev.elementwise();
       //! [errorComputation]

       //! [adaptRefinementPart]
       // Mark elements for refinement, based on the computed local errors and
       // the refinement-criterion and -parameter.
       std::vector<bool> elMarked( eltErrs.size() );
       gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
        if (IntensityMAE >1.){
            double Minvalue     = *std::max_element(eltErrs.begin(), eltErrs.end());                
            for(size_t i=0; i<eltErrs.size(); ++i)
            {
                if (elMarked[i] == true and Minvalue > eltErrs[i]) // Avoid numerical issues
                {
                    Minvalue = eltErrs[i];
                }
            }
            for(size_t i=0; i<eltErrs.size(); ++i)
            {
                if (Minvalue/pow(DoFPDE[refLoop],adaptRefParam) < eltErrs[i]) // Avoid numerical issues
                {
                    elMarked[i] = true;
                }
            }
        }
       gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";

       // Refine the marked elements with a 1-ring of cells around marked elements
       gsRefineMarkedElements( cdrAss.multiBasis(), elMarked, NumArMarEl );
       //! [adaptRefinementPart]


       //! [repairInterfaces]
       // Call repair interfaces to make sure that the new meshes
       // match along patch interfaces.
       cdrAss.multiBasis().repairInterfaces( Psi.interfaces() );

       //! [repairInterfaces]

       //! [refreshAssembler]
       cdrAss.refresh();
       //! [refreshAssembler]
        NumArMarEl = NumArMarEl + FactRefPar;
        // if (r%2==0){
        FactRefPar = 1*FactRefPar;
        //}

       //! [Export to Paraview]
       // Export the final solution
       if( plot && refLoop == numLRefine)
       {
           // Write the computed solution to paraview files
           gsWriteParaview<>( solField, "adaptRef", 1000, true);
       }
       //! [Export to Paraview]

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

    if (errorsave)
    {
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("error_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: " << adaptRefParam <<" "<< NumArMarEl <<" " << IntensityMAE << " \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        outFile << "#H1_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }
    }
    else
    {
        gsInfo << "Errors are not saved. To save them, try with --errorsave.\n";
    }

    //! [Error and convergence rates]
    if (numLRefine>0 && false)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numLRefine).array()  /
                   l2err.tail(numLRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numLRefine).array() /
                  h1err.tail(numLRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]
       //! [Export to Paraview]
       // Export the final solution
    if(plot){
        //------------------------------------
        gsInfo<<"Plotting in Paraview...\n";
        // Run paraview
        gsFileManager::open("adaptRef.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main