/** @file iterativeSolvers_example.cpp

    @brief Tutorial on solving a Poisson problem with iterative solvers and preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/

//! [Include namespace]
# include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    // TODO: choose geometry
    index_t numRefine  = 2;
    // TODO: choose degree
    // TODO: choose smoothness
    bool useNitsche = false;
    std::string preconder("hyb");
    real_t tol = 1.e-8;
    index_t maxIter = 200;
    bool plot = false;
    // TODO: update docs

    gsCmdLine cmd("Tutorial on solving a Poisson problem with iterative solvers and preconditioners." );
    // TODO: choose geometry
    cmd.addInt   ("r", "refine",  "Number of refinement levels",                            numRefine );
    // TODO: choose degree
    // TODO: choose smoothness
    cmd.addSwitch(     "nitsche", "Use Nitsche approach to realize boundary conditions",    useNitsche);
    cmd.addReal  ("t", "tol",     "Tolerance for iterative solver",                         tol       );
    cmd.addString("",  "prec",    "Preconditioner to be used",                              preconder );
    cmd.addInt   ("",  "maxIter", "Maximum number of iterations",                           maxIter   );
    cmd.addSwitch(     "plot",    "Create a ParaView visualization file with the solution", plot      );
    cmd.getValues(argc,argv);

    if ( preconder != "none" && preconder != "j" && preconder != "gs" && preconder != "fd" && preconder != "hyb" )
    {
        gsInfo << "Unknwon preconditioner chosen. Known are only:\n"
                    "  \"none\" ... No preconditioner, i.e., identity matrix.\n"
                    "  \"j\"    ... Jacobi preconditioner.\n"
                    "  \"gs\"   ... Symmetric Gauss Seidel preconditioner.\n"
                    "  \"fd\"   ... Fast diagonalization method.\n"
                    "  \"hyb\"  ... Hybric preconditioner: GS, Fast diagonalization, reverse GS.\n"
                    "\n";
        return EXIT_FAILURE;
    }

    //! [Parse command line]

    //! [Function data]
    // Define source function
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",2);
    // For homogeneous term, we can use this (last argument is the dimension of the domain)
    //gsConstantFunction<> f(0.0, 0.0, 2);

    // We also define exact solution and use it for the Dirichlet
    // boundary conditions
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",2);

    // Print out source function and solution
    gsInfo << "Source function :" << f << "\n";
    gsInfo << "Exact solution  :" << g << "\n\n";
    //! [Function data]

    //! [Geometry data]
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;

    // For single patch unit square of quadratic elements, use:
    if (true ) patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineQuarterAnnulus(2));

    // Geometry can also be read from file (if gsMultiPatch):
    if (false) gsReadFile<>("planar/lshape_p2.xml", patches);

    // Create 4 (2 x 2) patches of squares:
    //
    // Square/patch 0 is in lower left  corner
    // Square/patch 1 is in upper left  corner
    // Square/patch 2 is in lower right corner
    // Square/patch 3 is in upper right corner
    //
    // The last argument scale the squares such that we
    // get the unit square as domain.
    if (false) patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);

    gsInfo << "The domain is a "<< patches <<"\n";
    //! [Geometry data]


    //! [Boundary conditions]
    // Define Boundary conditions. Note that if one boundary is
    // "free", eg. if no condition is defined, then it is a natural
    // boundary (zero Neumann condition)
    gsBoundaryConditions<> bcInfo;

    // Here, we just define Dirichlet boundary conditions everywhere,
    // based on the made up solution
    for (gsMultiPatch<>::const_biterator
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }

    // We can also set the boundary conditions directly.
    //
    // Also, remember that a pure Neumann problem has no unique
    // solution, thereforer implies a singular matrix. In this case
    // a corner DoF can be fixed to a given value to obtain a unique solution.
    // (example: bcInfo.addCornerValue(boundary::southwest, value, patch);)
    //
    // For the case of the BSplineSquareGrid, we could define for example
    // mixed Dirichlet and Neumann boundary conditions
    //
    // bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);
    // bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g);
    //
    // bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
    // bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);
    //
    // gsFunctionExpr<> hEast ("1*pi*cos(pi*1)*sin(pi*2*y)",2);
    // gsFunctionExpr<> hSouth("-pi*2*sin(pi*x*1)",2);
    //
    // bcInfo.addCondition(3, boundary::east,  condition_type::neumann, &hEast);
    // bcInfo.addCondition(2, boundary::east,  condition_type::neumann, &hEast);
    //
    // bcInfo.addCondition(0, boundary::south, condition_type::neumann, &hSouth);
    // bcInfo.addCondition(2, boundary::south, condition_type::neumann, &hSouth);
    //! [Boundary conditions]

    //! [Refinement]
    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    // h-refine each basis
    for (index_t i = 0; i < numRefine; ++i)
      refine_bases.uniformRefine();

    // Number for p-refinement of the computational (trial/test) basis.
    index_t numElevate = 2; //TODO

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        index_t max_tmp = refine_bases.minCwiseDegree();

        // Elevate all degrees uniformly
        max_tmp += numElevate;
        refine_bases.setDegree(max_tmp);
    }
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
    const dirichlet::strategy dir = useNitsche ? dirichlet::nitsche : dirichlet::elimination;

    gsPoissonAssembler<real_t> assembler(
        patches,
        refine_bases,
        bcInfo,
        f,
        dir,
        iFace::glue
    );

    // Generate system matrix and load vector
    gsInfo<< "Assembling...\n";
    assembler.assemble();
    gsInfo << "Assembled a system (matrix and load vector) with "
           << assembler.numDofs() << " dofs.\n";
    //! [Assemble]

    //! [Solve]
    // Initialize the conjugate gradient solver
    gsInfo << "Solving...\n";
    gsLinearOperator<>::Ptr preconditioner;

    if (refine_bases.nBases() > 1 && ( preconder=="fd" || preconder=="hyb" ) )
    {
        gsInfo << "The chosen preconditioner only works for single-patch geometries.\n";
        return EXIT_FAILURE;
    }

    if (preconder=="j")
        preconditioner = makeJacobiOp( assembler.matrix() );
    else if (preconder=="gs")
        preconditioner = makeSymmetricGaussSeidelOp( assembler.matrix() );
    else if (preconder=="fd")
        preconditioner = gsParameterDomainPreconditioners<>(refine_bases[0],bcInfo,dir).getFastDiagonalizationOp();
    else if (preconder=="hyb")
        preconditioner = gsCompositionOfPreconditionersOp<>::make(
            makeGaussSeidelOp( assembler.matrix() ),
            gsPreconditionerFromOp<>::make(
                makeMatrixOp( assembler.matrix() ),
                gsParameterDomainPreconditioners<>(refine_bases[0],bcInfo,dir).getFastDiagonalizationOp()
            ),
            makeReverseGaussSeidelOp( assembler.matrix() )
        );

    gsConjugateGradient<> solver( assembler.matrix(), preconditioner );
    solver.setTolerance(tol);
    solver.setMaxIterations(maxIter);
    gsMatrix<> solVector, errorHistory;
    solver.solveDetailed( assembler.rhs(), solVector, errorHistory );

    // Checking for success and printing corresponding messages:
    bool success = solver.error() <= solver.tolerance();
    if ( success )
        gsInfo << "Solved the system with CG using ";
    else
        gsInfo << "CG did not reach the desired error goal after ";

    gsInfo << ( errorHistory.rows() - 1 ) << " iterations:\n";
    if (errorHistory.rows() < 20)
        gsInfo << errorHistory.transpose() << "\n\n";
    else
        gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";
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
        gsInfo << "Plotting in Paraview.\n";
        gsWriteParaview<>(sol, "poisson2d", 1000);
        const gsField<> exact( assembler.patches(), g, false );
        gsWriteParaview<>( exact, "poisson2d_exact", 1000);

        system("paraview poisson2d.pvd &");
        system("paraview poisson2d_exact.pvd &");
        //! [Plot in Paraview]
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;

}
