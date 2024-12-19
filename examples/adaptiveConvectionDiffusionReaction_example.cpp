/** @file adaptiveConvectionDiffusionReaction_example.cGm

    @brief Tutorial on how to use G+Smo to solve a convection-diffusion-reaction problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

//! [Include namespace]
# include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace std;
using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
   //! [Parse command line]
   bool plot                 = false;
   bool notAssemble          = false;
   // Number of initial uniform refinement steps:
   int numInitUniformRefine  = 4;
   //! [adaptRefSettings]
   // Number of refinement loops to be done
   int numRefinementLoops    = 5;

   gsCmdLine cmd("Example for solving a convection-diffusion problem.");
   cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
   cmd.addSwitch("notAssemble", "Use a solver based on predefined PDEs or assemble one manually", notAssemble);
   cmd.addInt("l","numRefinementLoops", "Number of local h-refinement loops", numRefinementLoops);
   cmd.addInt("u","numInitUniformRefine", "Number of uniform h-refinement loops", numInitUniformRefine);


   try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
   //! [Parse command line]

   // --------------- specify exact solution and right-hand-side ---------------

   //! [Function data]
    /* *********************** Test 1 ************************************* */
//     /// Stabilisation method
//     auto Stabilizationtype = stabilizerCDR::SUPG;
//    // Define exact solution (will be used for specifying Dirichlet boundary conditions
//    //gsFunctionExpr<> g("if( y<-x/2-1/2, 1, 0 )", 2);
//    gsFunctionExpr<> g("if( y>=0, if( x <=-1, 1, 0 ), 0 )", 2);
//    // Define source function
//    gsFunctionExpr<> rhs("0",2);

//    // diffusion coefficient:
//    gsFunctionExpr<> coeff_diff("0.000001","0","0","0.000001",2);
//    // convection coefficient:
//    gsFunctionExpr<> coeff_conv("3/sqrt(13)","-2/sqrt(13)",2);
//    // reaction coefficient:
//    gsFunctionExpr<> coeff_reac("0",2);
//    //! [Function data]

//    // Print out source function and solution
//    gsInfo<<"Source function " << rhs << "\n";
//    gsInfo<<"Dirichlet boundary conditions "  << g << "\n\n";
    /* *********************** End Assemble Test 1 *********************** */

   //! [Function data]
    /* *********************** Test 2 :  Poisson equation ************************************* */
    // Define Stabilization method 
    auto Stabilizationtype = stabilizerCDR::none;
    // Dirichlet boundary condition
    gsFunctionExpr<> g("1./(1.+exp((y - x  - 0.2)/0.01))", 2);
    // Define source function
    gsFunctionExpr<> rhs("4.12230724487712e-5*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2 - 1.69934170211664e-13*exp(-200.0*x + 200.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**3",2);
    // convection coefficient
    gsFunctionExpr<> coeff_conv("0.","0.",2);
    // diffusion coefficient:
    gsFunctionExpr<> coeff_diff("1.","0","0","1.",2);
    gsFunctionExpr<> coeff_diffMax("1.",2);// For a posteriori error estimate
    // reaction coefficient:
    gsFunctionExpr<> coeff_reac("0",2);
    // Manufactured solition
    gsFunctionExpr<> Sol_ex("1./(1.+exp((y - x  - 0.2)/0.01))",2);
    /* *********************** End Assemble Test 2 *********************** */

   // Read geometry from file
   //! [GetGeometryData]
   // Read xml and create gsMultiPatch
   // --------------- read geometry from file ---------------
//    string fileSrc( "planar/lshape2d_3patches_thb.xml" );
//    gsMultiPatch<real_t> patches;
//    gsReadFile<real_t>( fileSrc, patches);

    // .... one single patch
    gsMultiPatch<> patches, mptp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // ... need regularity to be at least C^1
    for(size_t i =0; i<mptp.nPatches(); ++i)
        patches.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mptp.patch(i)) ));
   //! [GetGeometryData]
   gsInfo << "The domain is a "<< patches <<"\n";

   //! [computeTopology]
   // Get all interfaces and boundaries:
   patches.computeTopology();
   //! [computeTopology]

   // --------------- add bonudary conditions ---------------
   //! [Boundary conditions]
   gsBoundaryConditions<> bcInfo;

   // For simplicity, set Dirichlet boundary conditions
   // given by exact solution g on all boundaries:
   for ( gsMultiPatch<>::const_biterator
            bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
   {
       bcInfo.addCondition( *bit, condition_type::dirichlet, &g ,0, false);
   }
   //! [Boundary conditions]

   // --------------- set up basis ---------------

   //! [GetBasisFromTHB]
   // Copy basis from the geometry
   gsMultiBasis<> bases( patches, true );//true: poly-splines (not NURBS)
   //! [GetBasisFromTHB]


   //! [initialRefinements]
   for (int i = 0; i < numInitUniformRefine; ++i){
     bases.uniformRefine();
     patches.uniformRefine();
    }

   // Specify cell-marking strategy...
   MarkingStrategy adaptRefCrit = PUCA;
   //MarkingStrategy adaptRefCrit = GARU;
   //MarkingStrategy adaptRefCrit = errorFraction;

   // ... and parameter.
   const real_t adaptRefParam = 0.7;
   //! [adaptRefSettings]
    // error and DoFs
    gsVector<>  l2err(numRefinementLoops+1);
    gsVector<int>  DoFPDE(numRefinementLoops+1);

   //! [initialRefinements]
   if (notAssemble){
   // --------------- define Pde ---------------
   //! [definePde]
   gsConvDiffRePde<real_t> cdrPde(patches, bcInfo, & coeff_diff,& coeff_conv, & coeff_reac, & rhs);
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
   for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
   {
   //! [beginRefLoop]
       gsInfo << "====== Loop " << refLoop << " of "
              <<numRefinementLoops<< " ======" << "\n";
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
       gsExprEvaluator<>::geometryMap Gm = ev.getMap(patches);
       gsExprEvaluator<>::variable is = ev.getVariable(solField.fields());
       auto ms = ev.getVariable(g, Gm);

       // Get the element-wise norms.
        // Recover rhs for Poisson equation
        auto SFunc        = ev.getVariable(rhs, Gm);
        // Coeffs for advection-reaction diffusion equation
        auto coeff_convGm = ev.getVariable(coeff_conv, Gm);
        auto coeff_diffGm = ev.getVariable(coeff_diffMax, Gm);
        auto coeff_reacGm = ev.getVariable(coeff_reac, Gm);
        //ev.integralElWise( ( coeff_diffGm * ilapl(is,Gm) - igrad(is, Gm)*coeff_convGm - coeff_reacGm * is + SFunc).sqNorm()*meas(Gm) );
        ev.integralElWise( ( ilapl(is,Gm) + SFunc).sqNorm()*meas(Gm) );
       const std::vector<real_t> & eltErrs  = ev.elementwise();
       //! [errorComputation]
        // Recover manufactured solution for Poisson equation
        auto u_ex         = ev.getVariable(Sol_ex, Gm);
        l2err[refLoop]    = math::sqrt(ev.integralElWise( (u_ex - is).sqNorm()*meas(Gm) ) );
        DoFPDE[refLoop]   = cdrAss.numDofs();
       // --------------- adaptive refinement ---------------

       //! [adaptRefinementPart]
       // Mark elements for refinement, based on the computed local errors and
       // the refinement-criterion and -parameter.
       std::vector<bool> elMarked( eltErrs.size() );
       gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
       gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";

       // Refine the marked elements with a 1-ring of cells around marked elements
       gsRefineMarkedElements( cdrAss.multiBasis(), elMarked, 1 );
       //! [adaptRefinementPart]


       //! [repairInterfaces]
       // Call repair interfaces to make sure that the new meshes
       // match along patch interfaces.
       cdrAss.multiBasis().repairInterfaces( patches.interfaces() );

       //! [repairInterfaces]

       //! [refreshAssembler]
       cdrAss.refresh();
       //! [refreshAssembler]

       //! [Export to Paraview]
       // Export the final solution
       if( plot && refLoop == numRefinementLoops )
       {
           // Write the computed solution to paraview files
           gsWriteParaview<>( solField, "adaptRef", 1000, true);
       }
       //! [Export to Paraview]

   }
   }
   else{
    gsSparseSolver<>::CGDiagonal solver;
    patches.addAutoBoundaries();
    gsBoundaryConditions<> bc;
    bc.setGeoMap(patches);
    // For simplicity, set Dirichlet boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
       bc.addCondition( *bit, condition_type::dirichlet, &g,0, false);
    }
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprEvaluator<> ev(A);

   for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
    {    
        // Elements used for numerical integration
        A.setIntegrationElements(bases);
        gsExprEvaluator<> ev(A);

        // Set the geometry map
        geometryMap Gm   = A.getMap(patches);

        // Set the discretization space
        space ru = A.getSpace(bases);
        // Solution vector and solution variable
        gsMatrix<> rsolVector;
        solution ru_sol = A.getSolution(ru, rsolVector);

        // Set the source term for Poisson equation
        auto SFunc      = A.getCoeff(rhs, Gm);  
        auto u_ex       = A.getCoeff(Sol_ex, Gm);    

       //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::

        ru.setup(bc, dirichlet::l2Projection, 0);

        // Compute the system matrix and right-hand side
        // Initialize the system
        A.initSystem();

        gsInfo<< "Solving PDEs " <<std::flush;
        gsInfo<< A.numDofs() <<std::flush;
        
        //auto h_Tau =  m_h/(2.*coeff_conv.squaredNorm()+m_h);

        A.assemble(
        igrad(ru, Gm) * igrad(ru, Gm).tr() * meas(Gm) //matrix
        // + ru * (coeff_conv * igrad(ru, Gm).tr()) * meas(Gm) //matrix
        // +coeff_reac * ru * ru.tr() * meas(Gm) //matrix
        ,
        ru * SFunc * meas(Gm) //rhs vector
        );

        gsInfo<< "." <<std::flush;// Assemblying done
        solver.compute( A.matrix() );
        rsolVector = solver.solve(A.rhs());
        gsInfo<< "." <<std::flush; // Linear solving done

        // ..
        DoFPDE[refLoop] = A.numDofs();
        l2err[refLoop]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(Gm) ) );

        gsInfo<< ". " <<std::flush; // Error computations done

        //! [beginRefLoop]
            gsInfo << "====== Loop " << refLoop << " of "
                    <<numRefinementLoops << " ======" << "\n";
        // --------------- error estimation/computation ---------------
        // Get the element-wise norms.
        ev.integralElWise( ( ilapl(ru_sol, Gm)+ SFunc ).sqNorm()*meas(Gm) );
        const std::vector<real_t> & eltErrs  = ev.elementwise();
        //! [errorComputation]

        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the computed local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked( eltErrs.size() );
        gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( bases, elMarked, 1);
        gsRefineMarkedElements( patches, elMarked, 1);
        }
   }
    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
   //! [Plot in Paraview]
   if( plot )
   {
       // Run paraview
       gsFileManager::open("adaptRef.pvd");
   }
   //! [Plot in Paraview]
   else
   {
       gsInfo<<"Done. No output created, re-run with --plot to get a ParaView "
               "file containing Plotting image data.\n";
   }
   return EXIT_SUCCESS;

}// end main
