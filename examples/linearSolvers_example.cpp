/** @file linearSolvers_example.cpp

    @brief Example on how the solve a system of linear equations with
    Gismo's own iterative solvers and Eigen's direct and iterative
    solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>
#include <gismo.h>


using namespace gismo;

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
    real_t pi = EIGEN_PI;

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
        rhs(k,0) = pi*pi*meshSize*meshSize*math::cos(meshSize*(1+k)*pi);

    //Compress the matrix
    mat.makeCompressed();
}

//Print out information of the iterative solver
template<typename SolverType>
void gsIterativeSolverInfo(const SolverType &method,
                           real_t error, double time, bool& succeeded )
{
    gsInfo << method.detail();
    gsInfo << " Computed res. error  : " << error << "\n";
    gsInfo << " Time to solve:       : " << time << "\n";
    if ( method.error() <= method.tolerance() && error <= method.tolerance() )
    {
        gsInfo <<" Test passed.\n";
    }
    else
    {
        gsInfo <<" TEST FAILED!\n";
        succeeded = false;
    }
}

int main(int argc, char *argv[])
{

    bool succeeded = true;

    index_t N = 100;
    //Tolerance
    real_t tol = std::pow(10.0, - REAL_DIG * 0.75);

    gsCmdLine cmd("Solves a 1D PDE with a Courant discretization with several solvers.");
    cmd.addInt ("n", "number", "Number of unknowns",                  N  );
    cmd.addReal("",  "tol",    "Tolerance for the iterative solvers", tol);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsSparseMatrix<> mat;
    gsMatrix<>       rhs;

    //Assemble the 1D Poisson equation
    poissonDiscretization(mat, rhs, N);

    //The minimal residual implementation requires a preconditioner.
    //We initialize an identity preconditioner (does nothing).
    gsLinearOperator<>::Ptr preConMat = gsIdentityOp<>::make(N);

    gsStopwatch clock;

    //initial guess
    gsMatrix<> x0;

#ifndef GISMO_WITH_GMP

    //Maximum number of iterations
    index_t maxIters = 3*N;

    gsOptionList opt = gsIterativeSolver<>::defaultOptions();
    opt.setInt ("MaxIterations", 3*N);
    opt.setReal("Tolerance"    , tol);

    ///----------------------GISMO-SOLVERS----------------------///
    gsInfo << "Testing G+Smo's solvers:\n";



    //Initialize the MinRes solver
    gsMinimalResidual<> MinRes(mat,preConMat);
    MinRes.setOptions(opt);

    //Set the initial guess to zero
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nMinRes: Started solving... ";
    clock.restart();
    MinRes.solve(rhs, x0);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(MinRes, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);


    //Initialize the MinRes solver with inexact residual error norm
    gsMinimalResidual<> MinResIR(mat,preConMat);
    MinResIR.setOptions(opt);
    MinResIR.setInexactResidual(true);

    //Set the initial guess to zero
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nMinResIR: Started solving... ";
    clock.restart();
    MinResIR.solve(rhs, x0);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(MinResIR, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);

    //Initialize the GMResSolver solver
    gsGMRes<> GMResSolver(mat,preConMat);
    GMResSolver.setOptions(opt);

    //Set the initial guess to zero
    x0.setZero(N,1);

    if (N < 200)
    {
        //Solve system with given preconditioner (solution is stored in x0)
        gsInfo << "\nGMRes: Started solving... ";
        clock.restart();
        GMResSolver.solve(rhs,x0);
        gsInfo << "done.\n";
        gsIterativeSolverInfo(GMResSolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);
    }
    else
        gsInfo << "\nSkipping GMRes due to high number of iterations...\n";


    //Initialize the CG solver
    gsConjugateGradient<> CGSolver(mat,preConMat);
    CGSolver.setOptions(opt);

    //Set the initial guess to zero
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nCG: Started solving... ";
    clock.restart();
    CGSolver.solve(rhs,x0);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(CGSolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);


    ///----------------------EIGEN-ITERATIVE-SOLVERS----------------------///
    gsInfo << "\nTesting Eigen's interative solvers:\n";

    gsSparseSolver<>::CGIdentity EigenCGIsolver;
    EigenCGIsolver.setMaxIterations(maxIters);
    EigenCGIsolver.setTolerance(tol);
    gsInfo << "\nEigen's CG + identity prec.: Started solving... ";
    clock.restart();
    EigenCGIsolver.compute(mat);
    x0 = EigenCGIsolver.solve(rhs);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(EigenCGIsolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);

    gsSparseSolver<>::CGDiagonal EigenCGDsolver;
    EigenCGDsolver.setMaxIterations(maxIters);
    EigenCGDsolver.setTolerance(tol);
    gsInfo << "\nEigen's CG + diagonal prec.: Started solving... ";
    clock.restart();
    EigenCGDsolver.compute(mat);
    x0 = EigenCGDsolver.solve(rhs);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(EigenCGDsolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);

    gsSparseSolver<>::BiCGSTABIdentity EigenBCGIsolver;
    EigenBCGIsolver.setMaxIterations(maxIters);
    EigenBCGIsolver.setTolerance(tol);
    gsInfo << "\nEigen's bi conjugate gradient stabilized + identity prec.: Started solving... ";
    clock.restart();
    EigenBCGIsolver.compute(mat);
    x0 = EigenBCGIsolver.solve(rhs);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(EigenBCGIsolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);

    gsSparseSolver<>::BiCGSTABDiagonal EigenBCGDsolver;
    EigenBCGDsolver.setMaxIterations(maxIters);
    EigenBCGDsolver.setTolerance(tol);
    gsInfo << "\nEigen's bi conjugate gradient stabilized solver + diagonal prec.: Started solving... ";
    clock.restart();
    EigenBCGDsolver.compute(mat);
    x0 = EigenBCGDsolver.solve(rhs);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(EigenBCGDsolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);

    gsSparseSolver<>::BiCGSTABILUT EigenBCGILUsolver;
    //EigenBCGILUsolver.preconditioner().setFillfactor(1);
    EigenBCGILUsolver.setMaxIterations(maxIters);
    EigenBCGILUsolver.setTolerance(tol);
    gsInfo << "\nEigen's bi conjugate gradient stabilized solver + ILU prec.: Started solving... ";
    clock.restart();
    EigenBCGILUsolver.compute(mat);
    x0 = EigenBCGILUsolver.solve(rhs);
    gsInfo << "done.\n";
    gsIterativeSolverInfo(EigenBCGILUsolver, (mat*x0-rhs).norm()/rhs.norm(), clock.stop(), succeeded);

    ///----------------------EIGEN-DIRECT-SOLVERS----------------------///
    gsSparseSolver<>::SimplicialLDLT EigenSLDLTsolver;
    gsInfo << "\nEigen's Simplicial LDLT: Started solving... ";
    clock.restart();
    EigenSLDLTsolver.compute(mat);
    x0 = EigenSLDLTsolver.solve(rhs);
    gsInfo << "done.\n";
    gsInfo << "Eigen's Simplicial LDLT: Time to solve       : " << clock.stop() << "\n";

    gsSparseSolver<>::QR solverQR;
    gsInfo << "\nEigen's QR: Started solving... ";
    clock.restart();
    solverQR.compute(mat);
    x0 = solverQR.solve(rhs);
    gsInfo << "done.\n";
    gsInfo << "Eigen's QR: Time to solve       : " << clock.stop() << "\n";

#endif

    gsSparseSolver<>::LU solverLU;
    gsInfo << "\nEigen's LU: Started solving... ";
    clock.restart();
    solverLU.compute(mat);
    x0 = solverLU.solve(rhs);
    gsInfo << "done.\n";
    gsInfo << "Eigen's LU: Time to solve       : " << clock.stop() << "\n";


    return succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
}
