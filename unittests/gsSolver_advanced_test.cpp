/** @file gsSolver_advanced_test.cpp

    @brief Additional comprehensive tests for gsSolver module (linear algebra solvers)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Copilot
*/

#include "gismo_unittest.h"

SUITE(gsSolver_advanced_test)
{
    TEST(gsConjugateGradient_IdentityMatrix)
    {
        int n = 10;
        gsSparseMatrix<> A(n, n);
        A.setIdentity();
        
        gsMatrix<> b(n, 1);
        b.setOnes();
        
        gsMatrix<> x(n, 1);
        x.setZero();
        
        gsConjugateGradient<> solver(A);
        solver.solve(b, x);
        
        // Solution should be b since A is identity
        for (int i = 0; i < n; i++)
            CHECK_CLOSE(x(i, 0), 1.0, 1e-8);
    }

    TEST(gsConjugateGradient_DiagonalMatrix)
    {
        int n = 5;
        gsSparseMatrix<> A(n, n);
        for (int i = 0; i < n; i++)
            A.insert(i, i) = i + 2.0;
        
        gsMatrix<> b(n, 1);
        for (int i = 0; i < n; i++)
            b(i, 0) = (i + 2.0) * (i + 1.0);
        
        gsMatrix<> x(n, 1);
        x.setZero();
        
        gsConjugateGradient<> solver(A);
        solver.solve(b, x);
        
        // Check solution
        for (int i = 0; i < n; i++)
            CHECK_CLOSE(x(i, 0), i + 1.0, 1e-8);
    }

    TEST(gsConjugateGradient_Convergence)
    {
        int n = 20;
        gsSparseMatrix<> A(n, n);
        
        // Create symmetric positive definite matrix
        for (int i = 0; i < n; i++)
        {
            A.insert(i, i) = 4.0;
            if (i > 0)
                A.insert(i, i-1) = -1.0;
            if (i < n-1)
                A.insert(i, i+1) = -1.0;
        }
        
        gsMatrix<> b(n, 1);
        b.setOnes();
        
        gsMatrix<> x(n, 1);
        x.setZero();
        
        gsConjugateGradient<> solver(A);
        solver.setTolerance(1e-10);
        solver.setMaxIterations(100);
        solver.solve(b, x);
        
        // Check that solver converged
        CHECK(solver.iterations() < 100);
    }

    TEST(gsMinRes_SimpleSystem)
    {
        int n = 8;
        gsSparseMatrix<> A(n, n);
        A.setIdentity();
        
        gsMatrix<> b(n, 1);
        b.setConstant(2.0);
        
        gsMatrix<> x(n, 1);
        x.setZero();
        
        gsMinRes<> solver(A);
        solver.solve(b, x);
        
        for (int i = 0; i < n; i++)
            CHECK_CLOSE(x(i, 0), 2.0, 1e-8);
    }

    TEST(gsGMRes_Restart)
    {
        int n = 10;
        gsSparseMatrix<> A(n, n);
        A.setIdentity();
        
        gsMatrix<> b(n, 1);
        b.setOnes();
        
        gsMatrix<> x(n, 1);
        x.setZero();
        
        gsGMRes<> solver(A);
        solver.setRestart(5);
        solver.solve(b, x);
        
        for (int i = 0; i < n; i++)
            CHECK_CLOSE(x(i, 0), 1.0, 1e-8);
    }

    TEST(gsSparseQR_Solve)
    {
        int n = 6;
        gsSparseMatrix<> A(n, n);
        
        // Create a simple SPD system
        for (int i = 0; i < n; i++)
        {
            A.insert(i, i) = 3.0;
            if (i > 0)
                A.insert(i, i-1) = 1.0;
        }
        
        gsMatrix<> b(n, 1);
        b.setConstant(1.0);
        
        gsSparseSolver<>::QR solver;
        solver.compute(A);
        gsMatrix<> x = solver.solve(b);
        
        CHECK(x.rows() == n);
    }

    TEST(gsSparseLU_Solve)
    {
        int n = 5;
        gsSparseMatrix<> A(n, n);
        
        for (int i = 0; i < n; i++)
        {
            A.insert(i, i) = 2.0;
            if (i > 0)
                A.insert(i, i-1) = -1.0;
            if (i < n-1)
                A.insert(i, i+1) = -1.0;
        }
        
        gsMatrix<> b(n, 1);
        b.setOnes();
        
        gsSparseSolver<>::LU solver;
        solver.compute(A);
        gsMatrix<> x = solver.solve(b);
        
        // Verify solution
        gsMatrix<> residual = A * x - b;
        CHECK(residual.norm() < 1e-8);
    }

    TEST(gsIterativeSolver_Tolerance)
    {
        int n = 10;
        gsSparseMatrix<> A(n, n);
        A.setIdentity();
        
        gsMatrix<> b(n, 1);
        b.setOnes();
        
        gsMatrix<> x(n, 1);
        x.setZero();
        
        gsConjugateGradient<> solver(A);
        
        double tol = 1e-12;
        solver.setTolerance(tol);
        CHECK_CLOSE(solver.tolerance(), tol, 1e-15);
    }

    TEST(gsIterativeSolver_MaxIterations)
    {
        gsSparseMatrix<> A(5, 5);
        A.setIdentity();
        
        gsConjugateGradient<> solver(A);
        
        int maxIter = 50;
        solver.setMaxIterations(maxIter);
        CHECK_EQUAL(solver.maxIterations(), maxIter);
    }

    /* 
     * Step-by-step instructions for additional complex gsSolver tests:
     * 
     * 1. Preconditioner tests:
     *    - Test Jacobi preconditioner with iterative solvers
     *    - Test incomplete LU (ILU) preconditioner
     *    - Test incomplete Cholesky preconditioner
     *    - Test algebraic multigrid preconditioner
     *    - Compare convergence with/without preconditioning
     * 
     * 2. gsSimpleOps and gsProductOp tests:
     *    - Create matrix-free linear operators
     *    - Test operator composition
     *    - Test with iterative solvers
     *    - Test scaling operators
     * 
     * 3. gsKroneckerOp tests:
     *    - Create Kronecker product operators for tensor problems
     *    - Test efficient tensor-structured solves
     *    - Test with Sylvester equations
     * 
     * 4. Advanced iterative solver tests:
     *    - Test BiCGSTAB for non-symmetric systems
     *    - Test different restart strategies for GMRES
     *    - Test flexible GMRES (FGMRES)
     *    - Test IDR(s) solver
     * 
     * 5. Direct solver tests with larger systems:
     *    - Test Pardiso solver (if available)
     *    - Test SuperLU solver
     *    - Test Umfpack solver
     *    - Compare performance of different direct solvers
     * 
     * 6. Block solver tests:
     *    - Create block-structured matrices
     *    - Test block preconditioners
     *    - Test Schur complement solvers
     *    - Test saddle point problems (Stokes, mixed formulations)
     * 
     * 7. Eigenvalue solver tests:
     *    - Compute smallest/largest eigenvalues
     *    - Test Arnoldi/Lanczos methods
     *    - Test shift-invert mode
     *    - Verify orthogonality of eigenvectors
     * 
     * 8. gsLinearOperator advanced tests:
     *    - Implement custom linear operators
     *    - Test apply() and apply Adjoint()
     *    - Test operator norms
     *    - Test with matrix-free methods
     */
}