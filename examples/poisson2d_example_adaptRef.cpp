/** @file poisson2d_example_adaptRef.cpp

    @brief Example for using the gsPoissonSolver with adaptive refinement with THB-splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/




#include <iostream>

#include <gismo.h>

#include <gsAssembler/gsAdaptiveRefUtils.h>

#include <gsAssembler/gsErrEstPoissonResidual.h>

using namespace std;
using namespace gismo;

//S.Kleiss
//
//This is a test example for a illustrating the adaptive
//refinement procedure implemented in the gsPoissonAssembler
//
//Flags, parameters, geometry and prescribed exact solution
//are specified within the main() function

int main()
{
  int result = 0;

  // Number of initial uniform mesh refinements
  int initUnifRef = 2;
  // Number of adaptive refinement loops
  int RefineLoopMax; // ...specified below with the examples

  // Flag for refinemet criterion
  // (see doxygen documentation of the free function
  // gsMarkElementsForRef explanation)
  const int refCriterion = 2;
  // Parameter for computing adaptive refinement threshold
  // (see doxygen documentation of the free function
  // gsMarkElementsForRef explanation)
  real_t refParameter;  // ...specified below with the examples

  // Flag whether final mesh should be plotted in ParaView
  const bool plot = true;

  // ****** Prepared test examples ******
  //
  // f       ... source term
  // g       ... exact solution
  // patches ... the computational domain given as object of gsMultiPatch
  //


  /*
  // ------ Example 1 ------

  // --- Unit square, with a spike of the source function at (0.25, 0.6)
  gsMFunctionExpr<>  f("if( (x-0.25)^2 + (y-0.6)^2 < 0.2^2, 1, 0 )", 2);
  gsMFunctionExpr<>  g("0", 2);
  gsMultiPatch<> patches( *safe(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,2.0,1.0) ));

  RefineLoopMax = 6;
  refParameter = 0.6;

  // ^^^^^^ Example 1 ^^^^^^
  //*/

//  /*
  // ------ Example 2 ------

  // The classical example associated with the L-Shaped domain.
  // The exact solution has a singularity at the re-entrant corner at the origin.

  gsMFunctionExpr<>  f("0", 2);
  gsMFunctionExpr<>  g("if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )", 2);

  gsMultiPatch<> patches( *safe(gsNurbsCreator<>::BSplineLShape_p2C1()) );

  RefineLoopMax = 4;
  refParameter = 0.85;

  // ^^^^^^ Example 2 ^^^^^^
  //*/

  cout<<"Source function "<< f << endl;
  cout<<"Exact solution "<< g <<".\n" << endl;

  // Define Boundary conditions
  gsBoundaryConditions<> bcInfo;
  // Dirichlet BCs
  bcInfo.addCondition( boundary::west,  boundary::dirichlet, &g );
  bcInfo.addCondition( boundary::east,  boundary::dirichlet, &g );
  bcInfo.addCondition( boundary::north, boundary::dirichlet, &g );
  bcInfo.addCondition( boundary::south, boundary::dirichlet, &g );

  gsTensorBSpline<2,real_t> * geo = dynamic_cast< gsTensorBSpline<2,real_t> * >( & patches.patch(0) );
  cout << " --- Geometry:\n" << *geo << endl;
  cout << "Number of patches: " << patches.nPatches() << endl;

  // With this gsTensorBSplineBasis, it's possible to call the THB-Spline constructor
  gsTHBSplineBasis<2,real_t> THB( geo->basis() );

  // Finally, create a vector (of length one) of this gsTHBSplineBasis
  gsMultiBasis<real_t> bases(THB);
  
  for (int i = 0; i < initUnifRef; ++i)
      bases.uniformRefine();

  // So, ready to start the adaptive refinement loop:
  for( int RefineLoop = 0; RefineLoop <= RefineLoopMax ; RefineLoop++ )
  {
      cout << "\n\n ====== Loop " << RefineLoop << " of " << RefineLoopMax << " ======" << endl << endl;

      // Create solver... maybe not the smartest thing to set up a new solver
      // in each iteration loop, but good enough for now.
      gsPoissonAssembler<real_t> pa(patches,bases,bcInfo,f,
                                    dirichlet::elimination,iFace::glue);

      // Assemble matrix and rhs
      pa.assemble();

      // Solve system
      gsMatrix<> solVector = Eigen::ConjugateGradient<gsSparseMatrix<> >(pa.matrix() ).solve( pa.rhs() );
      
      // Construct the solution for plotting the mesh later
      gsField<> * sol = pa.constructSolution(solVector);

      // Set up and compute the L2-error to the known exact solution
      gsNormL2<real_t> norm(*sol,g);
      gsErrEstPoissonResidual<real_t> errEst(*sol,g);

      norm.compute(1);
      errEst.compute(1);

      // Get the vector with element-wise local errors...
      const std::vector<real_t> & elError = norm.elementNorms();
      // ...or the vector with element-wise local error estimates.
      const std::vector<real_t> & elErrEst = errEst.elementNorms();

      // Mark elements for refinement, based on the computed local errors and
      // refCriterion and refParameter.
      std::vector<bool> elMarked( elError.size() );
      // Use the (in this case known) exact error...
      //gsMarkElementsForRef( elError, refCriterion, refParameter, elMarked);
      // ...or the error estimate.
      gsMarkElementsForRef( elErrEst, refCriterion, refParameter, elMarked);

      // Refine the elements of the mesh, based on elMarked.
      gsRefineMarkedElements( bases, elMarked);


      if ( (RefineLoop == RefineLoopMax) && plot)
      {
          // Write approximate solution to paraview files
          std::cout<<"Plotting in Paraview...\n";
          gsWriteParaview<>(*sol, "p2d_adaRef_sol", 1001, true);
          // Run paraview and plot the last mesh
          result = system("paraview p2d_adaRef_sol0_mesh.vtp &");
      }

      delete sol;
  }

  cout << "Test is done: Exiting" << endl;
  return  result;
}
