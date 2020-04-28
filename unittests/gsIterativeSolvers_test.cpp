/** @file gsIterativeSolvers_test.cpp

    @brief Tests for iterative solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/

#include "gismo_unittest.h"


//Create a tri-diagonal matrix with -1 of the off diagonals and 2 in the diagonal.
//This matrix is equivalent to discretizing the 1D Poisson equation with homogenius
//Dirichlet boundary condition using a finite difference method. It is a SPD matrix.
//The solution is sin(pi*x);
void poissonDiscretization(gsSparseMatrix<> &mat, gsMatrix<> &rhs, index_t N)
{
    rhs.setZero(N,1);

    mat.resize(N,N);
    mat.setZero();
    real_t meshSize = 1./(N+1);

    //Reserving space in the sparse matrix (Speeds up the assemble time of the matrix)
    mat.reservePerColumn( 3 ); //Reserve 3 non-zero entry per column

    mat(0,0) = 2;
    mat(0,1) = -1;
    mat(N-1, N-1) = 2;
    mat(N-1, N-2) = -1;
    for (index_t k = 1; k < N-1; ++k)
    {
        mat(k,k) = 2;
        mat(k,k-1) = -1;
        mat(k,k+1) = -1;
    }

    for (index_t k = 0; k < N; ++k)
        rhs(k,0) = EIGEN_PI*EIGEN_PI*meshSize*meshSize*math::cos(meshSize*(1+k)*EIGEN_PI);

    //Compress the matrix
    mat.makeCompressed();
}

SUITE(gsIterativeSolvers_test)
{

    TEST(Gradient_test)
    {
        index_t          N = 10;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.25);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsGradientMethod<>::defaultOptions();
        opt.setInt ("MaxIterations", 200 );
        opt.setReal("Tolerance"    , tol );

        gsGradientMethod<> solver(mat);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(Gradient_fixed_test)
    {
        index_t          N = 10;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.25);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsGradientMethod<>::defaultOptions();
        opt.setInt   ("MaxIterations"   , 200         );
        opt.setReal  ("Tolerance"       , tol         );
        opt.setSwitch("AdaptiveStepSize", false       );
        opt.setReal  ("StepSize"        , (real_t)1/2 );

        gsGradientMethod<> solver(mat);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(CG_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsConjugateGradient<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        // This test checks also that we can use our iterative solver
        // as a linear operator
        gsIterativeSolverOp< gsConjugateGradient<> >::uPtr solverOp
            = gsIterativeSolverOp< gsConjugateGradient<> >::make(mat);

        solverOp->solver().setOptions(opt);
        solverOp->apply(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(MinRes_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsMinimalResidual<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsMinimalResidual<> solver(mat);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs, x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(MinRes_InexactResidual_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsMinimalResidual<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsMinimalResidual<> solver(mat);
        solver.setOptions(opt);
        solver.setInexactResidual(true);

        x.setZero(N,1);
        solver.solve(rhs, x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(GMRes_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsGMRes<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsGMRes<> solver(mat);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(CG_Jacobi_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsConjugateGradient<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsLinearOperator<>::Ptr preConMat = makeJacobiOp(mat);
        gsConjugateGradient<> solver(mat,preConMat);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(CG_SGS_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsConjugateGradient<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsLinearOperator<>::Ptr precon = makeSymmetricGaussSeidelOp(mat);
        gsConjugateGradient<> solver(mat,precon);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(GMRES_GS_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsGMRes<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsLinearOperator<>::Ptr precon = makeGaussSeidelOp(mat);
        gsGMRes<> solver(mat,precon);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs, x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(GMRES_RGS_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsGMRes<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsLinearOperator<>::Ptr precon = makeReverseGaussSeidelOp(mat);
        gsGMRes<> solver(mat,precon);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs, x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

    TEST(MinRes_Rich_test)
    {
        index_t          N = 100;
        real_t           tol = std::pow(10.0, - REAL_DIG * 0.75);

        gsSparseMatrix<> mat;
        gsMatrix<>       rhs;
        gsMatrix<>       x;

        poissonDiscretization(mat, rhs, N);

        gsOptionList opt = gsMinimalResidual<>::defaultOptions();
        opt.setInt ("MaxIterations", N  );
        opt.setReal("Tolerance"    , tol);

        gsRichardsonOp< gsSparseMatrix<> >::Ptr precon = gsRichardsonOp< gsSparseMatrix<> >::make(mat);
        precon->setNumOfSweeps(3);
        precon->setDamping((real_t)1/5);

        gsMinimalResidual<> solver(mat,precon);
        solver.setOptions(opt);

        x.setZero(N,1);
        solver.solve(rhs,x);

        CHECK( (mat*x-rhs).norm()/rhs.norm() <= tol );
    }

}
