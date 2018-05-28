/** @file poisson_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

//! [Include namespace]
# include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Function data]
    // Define source function
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
                              "((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)",2);
    // For homogeneous term, we can use this (last argument is the dimension of the domain)
    //gsConstantFunction<> f(0.0, 0.0, 2);

    // Define exact solution (optional)
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",
                              "sin(pi*x*3)*sin(pi*y*4)-pi/10",2);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n\n";
    //! [Function data]

    //! [Geometry data]
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;
    // Create 4 (2 x 2) patches of squares:
    //
    // Square/patch 0 is in lower left  corner
    // Square/patch 1 is in upper left  corner
    // Square/patch 2 is in lower right corner
    // Square/patch 3 is in upper right corner
    //
    // The last argument scale the squares such that we
    // get the unit square as domain.
    patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);
    gsInfo << "The domain is a "<< patches <<"\n";
    //! [Geometry data]

    // For single patch unit square of quadratic elements use (Note:
    // you need to update the bounadry conditions section for this to
    // work properly!) :
    // patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    // Geometry can also be read from file (if gsMultiPatch):
    // gsReadFile<>("planar/lshape_p2.xml", patches);

    // Define Boundary conditions. Note that if one boundary is
    // "free", eg. if no condition is defined, then it is a natural
    // boundary (zero Neumann condition)
    // Also, remember that a pure Neumann problem has no unique
    // solution, thereforer implies a singular matrix. In this case
    // a corner DoF can be fixed to a given value to obtain a unique solution.
    // (example: bcInfo.addCornerValue(boundary::southwest, value, patch);)

    //! [Boundary conditions]
    gsBoundaryConditions<> bcInfo;
    // Every patch with a boundary need to be specified. In this
    // there are in total 8 sides (two for each patch)

    // Dirichlet Boundary conditions
    // First argument is the patch number
    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g);

    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);

    // Neumann Boundary conditions
    gsFunctionExpr<> hEast ("1*pi*cos(pi*1)*sin(pi*2*y)", "3*pi*cos(pi*3)*sin(pi*4*y)",2);
    gsFunctionExpr<> hSouth("-pi*2*sin(pi*x*1)","-pi*4*sin(pi*x*3)",2);

    bcInfo.addCondition(3, boundary::east,  condition_type::neumann, &hEast);
    bcInfo.addCondition(2, boundary::east,  condition_type::neumann, &hEast);

    bcInfo.addCondition(0, boundary::south, condition_type::neumann, &hSouth);
    bcInfo.addCondition(2, boundary::south, condition_type::neumann, &hSouth);
    //! [Boundary conditions]

    /*
      //Alternatively: You can automatically create Dirichlet boundary
      //conditions using one function (the exact solution) for all
      //boundaries like this:

    for (gsMultiPatch<>::const_biterator
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }
    */

    //! [Refinement]
    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    // Number for h-refinement of the computational (trial/test) basis.
    const int numRefine  = 2;

    // Number for p-refinement of the computational (trial/test) basis.
    const int degree     = 2;

    // h-refine each basis (4, one for each patch)
    for ( int i = 0; i < numRefine; ++i)
      refine_bases.uniformRefine();

    // k-refinement (set degree)
    for ( size_t i = 0; i < refine_bases.nBases(); ++ i )
        refine_bases[i].setDegreePreservingMultiplicity(degree);

    //! [Refinement]

    ////////////// Setup solver and solve //////////////
    // Initialize Solver
    // Setup method for handling Dirichlet boundaries, options:
    //
    // * elimination: Eliminate the Dirichlet DoFs from the linear system.
    //
    // * nitsche: Keep the Dirichlet DoFs and enforce the boundary
    //
    // condition weakly by a penalty term.
    // Setup method for handling patch interfaces, options:
    //
    // * glue:Glue patches together by merging DoFs across an interface into one.
    //   This only works for conforming interfaces
    //
    // * dg: Use discontinuous Galerkin-like coupling between adjacent patches.
    //       (This option might not be available yet)
    //! [Assemble]
    gsPoissonAssembler<real_t> assembler(patches,refine_bases,bcInfo,f,
                                       //dirichlet::elimination, iFace::glue);
                                         dirichlet::nitsche    , iFace::glue);

    // Generate system matrix and load vector
    gsInfo<< "Assembling...\n";
    assembler.assemble();
    gsInfo << "Have assembled a system (matrix and load vector) with "
           << assembler.numDofs() << " dofs.\n";
    //! [Assemble]

    //! [Solve]
    // Initialize the conjugate gradient solver
    gsInfo << "Solving...\n";
    gsSparseSolver<>::CGDiagonal solver( assembler.matrix() );
    gsMatrix<> solVector = solver.solve( assembler.rhs() );
    gsInfo << "Solved the system with CG solver.\n";
    //! [Solve]

    //! [Construct solution]
    // Construct the solution as a scalar field
    gsMultiPatch<> mpsol;
    assembler.constructSolution(solVector, mpsol);
    gsField<> sol( assembler.patches(), mpsol);
    //! [Construct solution]

    if (plot)
    {
        //! [Plot in Paraview]
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "poisson2d", 1000);
        const gsField<> exact( assembler.patches(), g, false );
        gsWriteParaview<>( exact, "poisson2d_exact", 1000);

        // Run paraview
        gsFileManager::open("poisson2d.pvd");
        gsFileManager::open("poisson2d_exact.pvd");
        //! [Plot in Paraview]
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    }
    return EXIT_SUCCESS;

}// end main
