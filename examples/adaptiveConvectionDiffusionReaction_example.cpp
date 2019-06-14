/** @file adaptiveConvectionDiffusionReaction_example.cpp

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
   bool plot = false;

   gsCmdLine cmd("Example for solving a convection-diffusion problem.");
   cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
   try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
   //! [Parse command line]

   // --------------- specify exact solution and right-hand-side ---------------

   //! [Function data]
   // Define exact solution (will be used for specifying Dirichlet boundary conditions
   //gsFunctionExpr<> g("if( y<-x/2-1/2, 1, 0 )", 2);
   gsFunctionExpr<> g("if( y>=0, if( x <=-1, 1, 0 ), 0 )", 2);
   // Define source function
   gsFunctionExpr<> rhs("0",2);

   // diffusion coefficient:
   gsFunctionExpr<> coeff_diff("0.000001","0","0","0.000001",2);
   // convection coefficient:
   gsFunctionExpr<> coeff_conv("3/sqrt(13)","-2/sqrt(13)",2);
   // reaction coefficient:
   gsFunctionExpr<> coeff_reac("0",2);
   //! [Function data]

   // Print out source function and solution
   gsInfo<<"Source function " << rhs << "\n";
   gsInfo<<"Dirichlet boundary conditions "  << g << "\n\n";


   // --------------- read geometry from file ---------------

   // Read geometry from file
   //! [GetGeometryData]
   // Read xml and create gsMultiPatch
   string fileSrc( "planar/lshape2d_3patches_thb.xml" );
   gsMultiPatch<real_t> patches;
   gsReadFile<real_t>( fileSrc, patches);
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
       bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
   }
   //! [Boundary conditions]

   // --------------- define Pde ---------------
   //! [definePde]
   gsConvDiffRePde<real_t> cdrPde(patches, bcInfo, & coeff_diff,& coeff_conv, & coeff_reac, & rhs);
   //! [definePde]


   // --------------- set up basis ---------------

   //! [GetBasisFromTHB]
   // Copy basis from the geometry
   gsMultiBasis<> bases( patches );
   //! [GetBasisFromTHB]


   //! [initialRefinements]
   // Number of initial uniform refinement steps:
   int numInitUniformRefine  = 2;

   for (int i = 0; i < numInitUniformRefine; ++i)
     bases.uniformRefine();
   //! [initialRefinements]


   // --------------- set up adaptive refinement loop ---------------

   //! [adaptRefSettings]
   // Number of refinement loops to be done
   int numRefinementLoops = 3;

   // Specify cell-marking strategy...
   MarkingStrategy adaptRefCrit = PUCA;
   //MarkingStrategy adaptRefCrit = GARU;
   //MarkingStrategy adaptRefCrit = errorFraction;

   // ... and parameter.
   const real_t adaptRefParam = 0.7;

   //! [adaptRefSettings]


   //! [constructAssembler]
   // Construct assembler
   gsCDRAssembler<real_t> cdrAss( cdrPde, bases);
   // Set stabilization flag to 1 = SUPG
   cdrAss.options().setInt("Stabilization", stabilizerCDR::SUPG);
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
       gsExprEvaluator<>::variable ms = ev.getVariable(g, Gm);

       // Get the element-wise norms.
       ev.integralElWise( ( igrad(is,Gm) - igrad(ms)).sqNorm()*meas(Gm) );
       const std::vector<real_t> & eltErrs  = ev.elementwise();
       //! [errorComputation]

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
