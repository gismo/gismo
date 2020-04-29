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
    bool plot = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    bool last = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsMultiPatch<> mp0;
    mp0 = *gsNurbsCreator<>::BSplineRectangle(0,0,1,1);

    //! [Read input file]

    //! [Refinement]

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mp0.degreeElevate(numElevate);
    // h-refine each basis
    for (int r =0; r < numRefine-1; ++r)
        mp0.uniformRefine();

    // gsWriteParaview<>(mp0, "mp", 1000, true);

    // Cast all patches of the mp object to THB splines
    gsMultiPatch<> mp;
    gsTHBSpline<2,real_t> thb;
    for (index_t k=0; k!=mp0.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp0.patch(k));

        thb = gsTHBSpline<2,real_t>(*geo);
    // gsDebugVar(*geo);
    // gsDebugVar(*basis0);
        mp.addPatch(thb);
    }

    gsWriteParaview<>(mp, "mp", 1000, true);

    numRefine = 0;

    gsMultiBasis<> basis(mp);
    gsInfo<<"Basis Primal: "<<basis.basis(0)<<"\n";
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< basis.minCwiseDegree() <<"\n";


    gsTensorBSplineBasis<2,real_t> *basis0 = dynamic_cast< gsTensorBSplineBasis<2,real_t> * > (&mp0.basis(0));


    gsTHBSplineBasis<2,real_t> THBbasis(*basis0);
    THBbasis.refineBassFunction(0);

    gsWriteParaview<>(*basis0, "bases", 1000, true);









   // // --------------- set up adaptive refinement loop ---------------

   // //! [adaptRefSettings]
   // // Number of refinement loops to be done
   // int numRefinementLoops = 4;

   // // Specify cell-marking strategy...
   // MarkingStrategy adaptRefCrit = PUCA;
   // //MarkingStrategy adaptRefCrit = GARU;
   // //MarkingStrategy adaptRefCrit = errorFraction;

   // // ... and parameter.
   // const real_t adaptRefParam = 0.9;
   // //! [adaptRefSettings]


   // // --------------- adaptive refinement loop ---------------

   // //! [beginRefLoop]
   // for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
   // {
   // //! [beginRefLoop]

   //     // --------------- solving ---------------

   //     //! [solverPart]
   //     // Construct assembler
   //     gsPoissonAssembler<real_t> PoissonAssembler(patches,bases,bcInfo,f);
   //     PoissonAssembler.options().setInt("DirichletValues", dirichlet::l2Projection);

   //     // Generate system matrix and load vector
   //     PoissonAssembler.assemble();

   //     // Initialize the conjugate gradient solver
   //     gsSparseSolver<>::CGDiagonal solver( PoissonAssembler.matrix() );

   //     // Solve the linear system
   //     gsMatrix<> solVector = solver.solve( PoissonAssembler.rhs() );

   //     // Construct the isogeometric solution
   //     gsMultiPatch<> sol;
   //     PoissonAssembler.constructSolution(solVector, sol);
   //     // Associate the solution to the patches (isogeometric field)
   //     gsField<> solField(patches, sol);
   //     //! [solverPart]

   //     // --------------- error estimation/computation ---------------

   //     //! [errorComputation]
   //     // Compute the error in the H1-seminorm ( = energy norm in this example )
   //     // using the known exact solution.
   //     gsExprEvaluator<> ev;
   //     ev.setIntegrationElements(PoissonAssembler.multiBasis());
   //     gsExprEvaluator<>::geometryMap Gm = ev.getMap(patches);
   //     gsExprEvaluator<>::variable is = ev.getVariable(sol);
   //     gsExprEvaluator<>::variable ms = ev.getVariable(g, Gm);

   //     // Get the element-wise norms.
   //     ev.integralElWise( ( igrad(is,Gm) - igrad(ms)).sqNorm()*meas(Gm) );
   //     const std::vector<real_t> & eltErrs  = ev.elementwise();
   //     //! [errorComputation]

   //     // --------------- adaptive refinement ---------------

   //     //! [adaptRefinementPart]
   //     // Mark elements for refinement, based on the computed local errors and
   //     // the refinement-criterion and -parameter.
   //     std::vector<bool> elMarked( eltErrs.size() );
   //     gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);

   //     // Refine the marked elements with a 1-ring of cells around marked elements
   //     gsRefineMarkedElements( bases, elMarked, 1 );
   //     //! [adaptRefinementPart]


   //     //! [repairInterfaces]
   //     // Call repair interfaces to make sure that the new meshes
   //     // match along patch interfaces.
   //     bases.repairInterfaces( patches.interfaces() );
   //     //! [repairInterfaces]


   //     //! [Export to Paraview]
   //     // Export the final solution
   //     if( plot && refLoop == numRefinementLoops )
   //     {
   //         // Write the computed solution to paraview files
   //         gsWriteParaview<>(solField, "adaptRef", 1000, true);
   //     }
   //     //! [Export to Paraview]

   // }

   // //! [Plot in Paraview]
   // if( plot )
   // {
   //     // Run paraview
   //     gsFileManager::open("adaptRef.pvd");
   // }
   // //! [Plot in Paraview]
   // else
   // {
   //     gsInfo<<"Done. No output created, re-run with --plot to get a ParaView "
   //             "file containing Plotting image data.\n";
   // }

   return EXIT_SUCCESS;

}// end main
