/** @file adaptRefinementThb_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

//! [Include namespace]
# include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
   //! [Parse command line]
  bool plot = false, thb = false;
   // Number of refinement loops to be done
   int numElevate = 1, numRefinementLoops = 4;
   // Number of initial uniform refinement steps:
   int numInitUniformRefine  = 1;
   
   real_t adaptRefParam = 0.9;

   gsCmdLine cmd("Tutorial on solving a Poisson problem.");
   cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
   cmd.addSwitch("tbnp", "Use truncated HB splines", thb );
    cmd.addReal( "a", "adaptiveParam",
                "Refinement parameter)", adaptRefParam );

    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps)", numElevate );
    cmd.addInt( "u", "urefine",
                "Number of a priori uniform refinement steps)", numInitUniformRefine );
   cmd.addInt( "r", "refine",
                "Number of adaptive refinement loops)", numRefinementLoops );
   try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
   //! [Parse command line]

   // --------------- specify exact solution and right-hand-side ---------------

   //! [Function data]
   // Define exact solution (will be used for specifying Dirichlet boundary conditions

   //   gsFunctionExpr<> g("if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x)+3*pi)/3.0 ) )", 2);

   gsFunctionExpr<> g("(x+1)^2",2);
   gsFunctionExpr<> dist("x-.01","y-.01",2);

   // Define source function
   gsFunctionExpr<> f("0",2);
   //! [Function data]

   // Print out source function and solution
   gsInfo<<"Source function " << f << "\n";
   gsInfo<<"Exact solution "  << g << "\n\n";

   // --------------- read geometry from file ---------------

   //! [GetGeometryData]
   
   /*
   // Read xml and create gsMultiPatch
   std::string fileSrc( "planar/lshape2d_3patches_thb.xml" );
   gsMultiPatch<real_t> patches;
   gsReadFile<real_t>( fileSrc, patches);
   //! [GetGeometryData]
   gsInfo << "The domain is a "<< patches <<"\n";

   //! [computeTopology]
   // Get all interfaces and boundaries:
   patches.computeTopology();
   //! [computeTopology]
   */

   //! [GetGeometryDataTens]
   std::string fileSrcTens( "planar/lshape2d_3patches_tens.xml" );
   gsMultiPatch<real_t> patchesTens;
   gsReadFile<real_t>( fileSrcTens, patchesTens);
   patchesTens.computeTopology();
   //! [GetGeometryDataTens]

   // --------------- add bonudary conditions ---------------
   //! [Boundary conditions]
   gsBoundaryConditions<> bcInfo;

   // For simplicity, set Dirichlet boundary conditions
   // given by exact solution g on all boundaries:
   for ( gsMultiPatch<>::const_biterator
         bit = patchesTens.bBegin(); bit != patchesTens.bEnd(); ++bit)
   {
       bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
   }
   //! [Boundary conditions]


   gsVector<> pt(2); 
   pt.setConstant(.99);
   /*
   gsTensorBSpline<2> & p1 = 
     const_cast<gsTensorBSpline<2>&>(
     dynamic_cast<const gsTensorBSpline<2>&>
     (patchesTens.piece(1)) );
   gsTensorBSpline<2> & p0 = 
     const_cast<gsTensorBSpline<2>&>(
     dynamic_cast<const gsTensorBSpline<2>&>
     (patchesTens.piece(0)) );
   gsTensorBSpline<2> & p2 = 
     const_cast<gsTensorBSpline<2>&>(
     dynamic_cast<const gsTensorBSpline<2>&>
     (patchesTens.piece(2)) );

   gsWriteParaview(patchesTens,"patches", 500, true);
   */
   // --------------- set up basis ---------------

   //! [GetBasisFromTens]
   // Copy tensor basis
   gsMultiBasis<real_t> basesTens( patchesTens );

   // Create a "basisContainer"
   std::vector< gsBasis<real_t>* > basisContainer;

   basesTens.degreeElevate(numElevate);
   for (int i = 0; i < numInitUniformRefine; ++i)
     basesTens.uniformRefine();

   gsTensorBSplineBasis<2> & p0 = 
     const_cast<gsTensorBSplineBasis<2>&>(
     dynamic_cast<const gsTensorBSplineBasis<2>&>(basesTens.piece(0)) );
   gsTensorBSplineBasis<2> & p1 = 
     const_cast<gsTensorBSplineBasis<2>&>(
     dynamic_cast<const gsTensorBSplineBasis<2>&>(basesTens.piece(1)) );
   gsTensorBSplineBasis<2> & p2 = 
     const_cast<gsTensorBSplineBasis<2>&>(
     dynamic_cast<const gsTensorBSplineBasis<2>&>(basesTens.piece(2)) );
   p0.insertKnot(pt[0],1,numElevate+1);
   p1.insertKnot(pt[0],0,numElevate+1); //knot, dir, mult
   p1.insertKnot(pt[1],1,numElevate+1);
   p2.insertKnot(pt[1],0,numElevate+1);

   // fill the "basisContainer" with patch-wise...
   for ( size_t i = 0; i < basesTens.nBases(); i++)
     if (thb)
       basisContainer.push_back(new gsTHBSplineBasis<2,real_t>(basesTens.basis(i)));
     else
       basisContainer.push_back(new gsHBSplineBasis<2,real_t>(basesTens.basis(i)));

   // finally, create the gsMultiBasis containing gsTHBSpline ...
   gsMultiBasis<real_t> bases( basisContainer, patchesTens );
   //! [GetBasisFromTens]

   // --------------- set up adaptive refinement loop ---------------

   //! [adaptRefSettings]

   // Specify cell-marking strategy...
   MarkingStrategy adaptRefCrit = PUCA;
   //MarkingStrategy adaptRefCrit = GARU;
   //MarkingStrategy adaptRefCrit = errorFraction;

   //! [adaptRefSettings]


   // --------------- adaptive refinement loop ---------------

   //! [beginRefLoop]
   for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
   {
   //! [beginRefLoop]

       // --------------- solving ---------------

       //! [solverPart]
       // Construct assembler
       gsPoissonAssembler<real_t> PoissonAssembler(patchesTens,bases,bcInfo,f);
       PoissonAssembler.options().setInt("DirichletValues", dirichlet::l2Projection);

       // Generate system matrix and load vector
       PoissonAssembler.assemble();

       // Initialize the conjugate gradient solver
       gsSparseSolver<>::CGDiagonal solver( PoissonAssembler.matrix() );

       // Solve the linear system
       gsMatrix<> solVector = solver.solve( PoissonAssembler.rhs() );

       // Construct the isogeometric solution
       gsMultiPatch<> sol;
       PoissonAssembler.constructSolution(solVector, sol);
       // Associate the solution to the patches (isogeometric field)
       gsField<> solField(patchesTens, sol);
       gsWriteParaview<>(solField, "adaptRef", 500, true);
       //! [solverPart]

       // --------------- error estimation/computation ---------------

       //! [errorComputation]
       // Compute the error in the H1-seminorm ( = energy norm in this example )
       // using the known exact solution.
       gsExprEvaluator<> ev;
       ev.setIntegrationElements(bases);
       gsExprEvaluator<>::geometryMap Gm = ev.getMap(patchesTens);
       gsExprEvaluator<>::variable is = ev.getVariable(sol);
       gsExprEvaluator<>::variable xy = ev.getVariable(dist, Gm);

       real_t val = ev.eval(is, pt  ,1).value();
       gsInfo << std::fixed << "-Value: "<< std::setprecision(16) 
	      << val        <<"   reference: 1.02679192610";
       gsInfo<< "( dofs: "<< bases.totalSize() <<").\n";

       //gsExprEvaluator<>::element el = ev.getElement();
       ev.options().setReal("quA", .0);
       ev.options().setInt ("quB", 1 );
       ev.minElWise( 1.0/xy.sqNorm() ); // distance from singularity (0,0)
       const std::vector<real_t> & eltErrs  = ev.elementwise();
       //! [errorComputation]
       
       if (refLoop == numRefinementLoops)
	 {
	   gsInfo<<std::setprecision(5) <<" Basis on 1 patch: "<< bases.piece(1) <<"\n";

	   //! [Export to Paraview]
	   // Export the final solution
	   if( plot )
	   {
	     // Write the computed solution to paraview files
	     gsWriteParaview<>(solField, "adaptRef", 1000, true);
	   }
	   //! [Export to Paraview]

	 break;
	 }
       // --------------- adaptive refinement ---------------

       //! [adaptRefinementPart]
       // Mark elements for refinement, based on the computed local errors and
       // the refinement-criterion and -parameter.
       std::vector<bool> elMarked( eltErrs.size() );

       gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);

       // Refine the marked elements with a 1-ring of cells around marked elements
       gsRefineMarkedElements( bases, elMarked, 1 );
       //! [adaptRefinementPart]


       //! [repairInterfaces]
       // Call repair interfaces to make sure that the new meshes
       // match along patch interfaces.
       bases.repairInterfaces( patchesTens.interfaces() );
       //! [repairInterfaces]
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
