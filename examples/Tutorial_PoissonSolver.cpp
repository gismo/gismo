/** @file Tutorial_PoissonSolver.cpp

    @brief Tutorial example on how the use one of the poisson assemblers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

// Tutorial on how to use gismo to solve the Poisson equation
// using the classes gsPoissonAssembler and gsPdeAssembler.
//
// Try also running ./bin/poissonvector.cpp -h to see
// command line options and how to load problem from
// file.
//
// We consider the vector valued Poisson problem

# include <gismo.h>


using std::cout;
using std::endl;
using namespace gismo;

int main(int argc, char *argv[]) 
{


    //////// Right-hand side and analytical solution ////////
    // Define source function
    gsMFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
                              "((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)");
    // For homogeneous term, we can use this (last argument dimension of the domain)
    //gsConstantFunction<> f(0.0, 0.0, 2);


    // Define exact solution (optional)
    gsMFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",
                              "sin(pi*x*3)*sin(pi*y*4)-pi/10");

    // Print out source function and solution
    cout<<"Source function "<< f << endl;
    cout<<"Exact solution "<< g <<"\n" << endl;
  

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> * patches;
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

    // For single patch unit square of quadratic elements use
    //patches = new gsMultiPatch<>(gsNurbsCreator<>::BSplineSquare(2));

    // Geometry can also be read from file
    //patch= gsReadFile<>( "filenameofgeometry.xml" )



    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
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
    gsMFunctionExpr<> hEast ("1*pi*cos(pi*1)*sin(pi*2*y)", "3*pi*cos(pi*3)*sin(pi*4*y)");
    gsMFunctionExpr<> hSouth("-pi*2*sin(pi*x*1)","-pi*4*sin(pi*x*3)");

    bcInfo.addCondition(3, boundary::east,  condition_type::neumann, &hEast);
    bcInfo.addCondition(2, boundary::east,  condition_type::neumann, &hEast);

    bcInfo.addCondition(0, boundary::south, condition_type::neumann, &hSouth);
    bcInfo.addCondition(2, boundary::south, condition_type::neumann, &hSouth);




    ////////////////////// Refinement h and p //////////////////////
    // Refinement

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( *patches );
    //std::vector<gsBasis<> *> refine_bases = patches->basesCopy();

    // Number for h-refinement of the computational (trail/test) basis.
    int numRefine  = 2;
    // Number for p-refinement of the computational (trail/test) basis.
    int numElevate = 2;

    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < numRefine; ++i)
      refine_bases.uniformRefine();

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.maxDegree(0);
        for (int j = 1; j < patches->parDim(); ++j )
            if ( max_tmp < refine_bases.maxDegree(j) )
                max_tmp = refine_bases.maxDegree(j);                
        
        // Elevate all degrees uniformly
        max_tmp += numElevate;
        refine_bases.setDegree(max_tmp);
    }

    ////////////// Setup solver and solve //////////////
    // Initialize Solver
    // Setup method for handling Dirichlet boundaries, options:
    //
    // elimination: Eliminate the Dirichlet DoFs from the linear system.
    //
    // nitsche: Keep the Dirichlet DoFs and enforce the boundary
    // condition weakly by a penalty term.
    // Setup method for handling patch interfaces, options:
    //
    // glue:Glue patches together by merging DoFs across an interface into one.
    // This only works for conforming interfaces
    //
    // dg: Use discontinuous Galerkin-like coupling between adjacent patches.
    // NB! This only works in combination with the Nitsche handling of Dirichlet BC
    gsPoissonAssembler<real_t> PoissonAssembler(*patches,refine_bases,bcInfo,f,
                                                dirichlet::nitsche, iFace::dg);

    // Generate system matrix and load vector
    std::cout<<"Assembling...\n";
    PoissonAssembler.assemble();

    // Initialize the congugate gradient solver
    std::cout<<"Solving...\n";
    Eigen::ConjugateGradient< gsSparseMatrix<> > solver( PoissonAssembler.matrix() );
    gsMatrix<> solVector = solver.solve( PoissonAssembler.rhs() );

    /////////////////// Post processing ///////////////////

    // Construct the solution as a scalar field
    gsField<>::uPtr sol = safe(PoissonAssembler.constructSolution(solVector));

    // Plot solution in paraview
    bool plot = false;
    int result = 0;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        std::cout<<"Plotting in Paraview...\n";
        gsWriteParaview<>(*sol, "poisson2d", 1000);
        const gsField<> exact( PoissonAssembler.patches(), g, false );
        gsWriteParaview<>( exact, "poisson2d_exact", 1000);

        // Run paraview
        result = system("paraview poisson2d.pvd &");
    }


    // Clean up and exit
    delete patches;

    return  result;
}

