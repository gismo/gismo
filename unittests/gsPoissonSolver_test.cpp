/** @file gsPoissonSolver_test.cpp

    @brief Tests for Poisson solvers.

    Verification test for the vector valued Poisson solver
    with the gsPdeAssembler structure.
    The test uses the Method of Manufactured Solutions (MMS)
    and measure the convergence rate of the L2 error. If the
    convergence rate is not sufficient, the test fails.

    Test inclueds:

    Inhomogeneous source term (right hand side)
    Inhomogeneous Dirichlet boundary conditions
    Inhomogeneous Neumann boundary conditions
    Nitche handling of Dirichlet BC
    DG handling of multipatch

    x,y \in (0,1)

    Source function:

        f(x,y)_x = ((pi*k0)^2 + (pi*k1)^2)*sin(pi*x*k0)*sin(pi*y*k1)
        f(x,y)_y = ((pi*k2)^2 + (pi*k3)^2)*sin(pi*x*k2)*sin(pi*y*k3)

    Solution:

        u(x,y)_x = sin(pi*x*k0)*sin(pi*y*k1)+pi/10
        u(x,y)_y = sin(pi*x*k2)*sin(pi*y*k3)-pi/10

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/

#include "gismo_unittest.h"
#include <gsSolver/gsSolverUtils.h>

using namespace gismo;


void runPoissonSolverTest( dirichlet::strategy Dstrategy, iFace::strategy Istrategy )
{
    int numRefine = 2;
    int maxIterations = 3;
    real_t convratetol = 1.85; // Convergence rate should be around 2

    // List of element sizes
    std::vector <real_t> h_list;

    // List of L2-errors
    std::vector <real_t> error_list;

    // Source function
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
                              "((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)",2);

    // Exact solution
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",
                              "sin(pi*x*3)*sin(pi*y*4)-pi/10",2);

    // Define Geometry (Unit square with 4 patches)
    gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);

    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    //Dirichlet BCs
    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);

    // Neumann BCs
    gsFunctionExpr<> gEast ("1*pi*cos(pi*1)*sin(pi*2*y)", "3*pi*cos(pi*3)*sin(pi*4*y)",2);
    gsFunctionExpr<> gSouth("-pi*2*sin(pi*x*1)","-pi*4*sin(pi*x*3)",2);

    bcInfo.addCondition(3, boundary::east,  condition_type::neumann, &gEast);
    bcInfo.addCondition(2, boundary::east,  condition_type::neumann, &gEast);
    bcInfo.addCondition(0, boundary::south, condition_type::neumann, &gSouth);
    bcInfo.addCondition(2, boundary::south, condition_type::neumann, &gSouth);

    // Copy bases for refinement
    gsMultiBasis<> refine_bases( patches );

    // Define discretization space by initial refining the basis of the geometry
    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();

    
    // linear solver
    gsSparseSolver<>::CGDiagonal solver;
    gsMatrix<> solVector;

    // Start Loop
    // For each iteration we h-refine the basis, then set up and solve the
    // Poisson problem. Then find the L2 error and element size.
    for (int i = 0; i < maxIterations; ++i)
    {
        refine_bases.uniformRefine();

        // Initilize Assembler
        gsPoissonAssembler<real_t> poisson(patches,refine_bases,bcInfo,f,Dstrategy,Istrategy);
        //gsPoissonAssembler<> poisson(*patches, bcInfo, refine_bases, f);
        
        // Assemble and solve
        poisson.assemble();
        solVector = solver.compute( poisson.matrix() ).solve( poisson.rhs() );

        // Find the element size
        real_t h = math::pow( (real_t) refine_bases.size(0), -1.0 / refine_bases.dim() );
        h_list.push_back(h);

        // Access the solution
        const gsField<> sol = poisson.constructSolution(solVector);

        // Find the l2 error
        real_t l2error = sol.distanceL2(g);
        error_list.push_back(l2error);

    }

    // Finding convergence rate in to differnt ways

    // Convergence rate found by last to errors are elemet sizes
    real_t convratelast = math::log(error_list[error_list.size()-2]/error_list[error_list.size()-1])/
            math::log(h_list[h_list.size()-2]/h_list[h_list.size()-1]);

    // Convergence rate found by a least square method
    // PS: This method of finding convergence rate mights require some initialt
    // refinement such that convergence have startet for the largest element
    // size value.
    real_t convrateavg =  gsSolverUtils<>::convergenceRateLS(error_list,h_list);

    CHECK( convrateavg > convratetol );
    CHECK( convratelast > convratetol );

    //Clean up before next test
    error_list.clear();
    h_list.clear();
}


SUITE(gsPoissonSolver_test)
{

    TEST(Galerkin_test)
    {
        runPoissonSolverTest(dirichlet::elimination, iFace::glue);
    }

    TEST(dG_test)
    {
        runPoissonSolverTest(dirichlet::elimination, iFace::dg);
    }

    TEST(Nitsche_dG_test)
    {
        runPoissonSolverTest(dirichlet::nitsche, iFace::dg);
    }
    
}

