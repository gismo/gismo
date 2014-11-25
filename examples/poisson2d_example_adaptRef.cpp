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

#include <gsIO/gsIOUtils.h>


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
  // Number of initial uniform mesh refinements
  int initUnifRef = 1;
  // Number of adaptive refinement loops
  const int RefineLoopMax = 6;

  // Flag for refinemet criterion
  // (see doxygen documentation for explanation)
  const int Ref_Crit = 2;
  // Parameter for computing adaptive refinement threshold
  // (see doxygen documentation for explanation)
  const real_t Ref_Alpha = 0.85;

  // Flag whether final mesh should be plotted in ParaView
  const bool plot = false;

  int result = 0;

  // *** Prepared test examples

  // --- Example with a spike at (0.25, -0.5)
//  gsMFunctionExpr<>  f("if( (x-0.25)^2 + (y+0.5)^2 < 0.05^2, 1, 0 )", 2);
//  gsMFunctionExpr<>  g("0", 2);

  // --- The classical example associated with the L-Shaped domain
  // Singularity at the re-entrant corner at the origin.
  gsMFunctionExpr<>  f("0", 2);
  gsMFunctionExpr<>  g("if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) +3*pi)/3.0 ) )", 2);

  cout<<"Source function "<< f << endl;
  cout<<"Exact solution "<< g <<".\n" << endl;



  // *** Create geometry
  //gsMultiPatch<> * patches;

  // L-shaped domain with C1-continuous discretization
  //patches = new gsMultiPatch<>(gsNurbsCreator<>::BSplineLShape_p2C1());
  gsMultiPatch<> patches( *safe(gsNurbsCreator<>::BSplineLShape_p2C1()) );
  // L-shaped domain with C0-continuous discretization (C0 at diagonal)
  // patches = new gsMultiPatch<>(gsNurbsCreator<>::BSplineLShape_p2C0());

  gsTensorBSpline<2,real_t> * geo = dynamic_cast< gsTensorBSpline<2,real_t> * >( & patches.patch(0) );
  cout << " --- Geometry:\n" << *geo << endl;
  cout << "Number of patches: " << patches.nPatches() << endl;

  // Define Boundary conditions
  gsBoundaryConditions<> bcInfo;
  // Dirichlet BCs
  bcInfo.addCondition( boundary::west,  boundary::dirichlet, &g );
  bcInfo.addCondition( boundary::east,  boundary::dirichlet, &g );
  bcInfo.addCondition( boundary::north, boundary::dirichlet, &g );
  bcInfo.addCondition( boundary::south, boundary::dirichlet, &g );

  // With this gsTensorBSplineBasis, it's possible to call the THB-Spline constructor
  gsTHBSplineBasis<2,real_t> THB( geo->basis() );

  // Finally, create a vector (of length one) of this gsTHBSplineBasis
  gsMultiBasis<> bases(THB);
  for (int i = 0; i < initUnifRef; ++i)
      bases.uniformRefine();
  cout << endl << " --- Initial basis of discretization space is: " << endl << bases << endl << endl;

  // So, ready to start the adaptive refinement loop:
  for( int RefineLoop = 0; RefineLoop <= RefineLoopMax ; RefineLoop++ )
  {
      cout << "\n\n ====== Loop " << RefineLoop << " of " << RefineLoopMax << " ======" << endl << endl;

      // Create solver... maybe not the smartest thing to set up a new solver
      // in each iteration loop, but good enough for now.
      gsPoissonAssembler<real_t> pa(patches,bases,bcInfo,f,
                                    dirichlet::nitsche,iFace::glue);

      // Assemble matrix and rhs
      pa.assemble();
      // Solve system
      gsMatrix<> solVector = Eigen::ConjugateGradient<gsSparseMatrix<> >(pa.matrix() ).solve( pa.rhs() );
      
      // Construct the solution for plotting the mesh later
      gsField<> * sol = pa.constructSolution(solVector);

      
      gsNormL2<real_t> norm(*sol,g);

      norm.compute(1);
      const std::vector<real_t> & elError = norm.elementNorms();
      std::cout<<"Error per element: "<< gsAsConstMatrix<>(elError) <<"\n";

      // Compute error indicator (at the moment (06.May 2014) only computation via bubble functions)
      //gsVector< gsMatrix<> > EIv = PoissonSolver.errorIndicator();

      // Refine cells based on the error estimation
      //PoissonSolver.adaptiveRefine( EIv, Ref_Crit, Ref_Alpha );


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
