/** @file iterativeSolvers.cpp

    @brief Example on how the solve a system of linear equation with the MINRES, GMRes and CG method.

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
    real_t pi = M_PI;

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
    {
        rhs(k,0) = pi*pi*meshSize*meshSize*math::cos(meshSize*(1+k)*pi);
    }

    //Compress the matrix
    mat.makeCompressed();
}

//Print out information of the iterative solver
void gsIterativeSolverInfo(const gsIterativeSolver &method, std::string methodName, double time)
{
    gsInfo << methodName +": System size         : " << method.size() << "\n";
    gsInfo << methodName +": Tolerance           : " << method.tolerance() << "\n";
    gsInfo << methodName +": Residual error      : " << method.error() << "\n";
    gsInfo << methodName +": Number of iterations: " << method.iterations() << "\n";
    gsInfo << methodName +": Time to solve:      : " << time << "\n";
}

int main(int argc, char *argv[])
{
    //Size of linear system
    index_t N = 100;
    if (argc >= 2)
        N = atoi(argv[1]);

    gsSparseMatrix<> mat;
    gsMatrix<>       rhs;

    //Assemble the 1D Poisson equation
    poissonDiscretization(mat, rhs, N);

    //The minimal residual implementation requires a preconditioner.
    //We initialize an identity preconditioner (does nothing).
    gsIdentityOp preConMat(N);

    //Tolerance
    real_t tol = std::pow(10.0, - REAL_DIG * 0.75);
    gsStopwatch clock;

    //initial guess
    gsMatrix<> x0;
    x0.setZero(N,1);

#ifndef GISMO_WITH_MPQ 

    //Maximum number of iterations
    index_t maxIters = 3*N;

    ///----------------------GISMO-SOLVERS----------------------///
    gsInfo << "Testing G+Smo's solvers:\n";

    //Initialize the MinRes solver
    gsMinimalResidual MinRes(mat,maxIters,tol);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nMinRes: Started solving..."  << "\n";
    clock.restart();
    MinRes.solve(rhs,x0,preConMat);
    gsIterativeSolverInfo(MinRes, "MinRes", clock.stop());

    //Initialize the CG solver
    gsGMRes GMResSolver(mat,maxIters,tol);

    //Set the initial guess to zero
    x0.setZero(N,1);

    if (N < 200)
    {
        //Solve system with given preconditioner (solution is stored in x0)
        gsInfo << "\nGMRes: Started solving..."  << "\n";
        clock.restart();
        GMResSolver.solve(rhs,x0,preConMat);
        gsIterativeSolverInfo(GMResSolver, "GMRes", clock.stop());
    }
    else
        gsInfo << "\nSkipping GMRes due to high number of iterations...\n";


    //Initialize the CG solver
    gsConjugateGradient CGSolver(mat,maxIters,tol);

    //Set the initial guess to zero
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nCG: Started solving..."  << "\n";
    clock.restart();
    CGSolver.solve(rhs,x0,preConMat);
    gsIterativeSolverInfo(CGSolver, "CG", clock.stop());


    ///----------------------EIGEN-ITERATIVE-SOLVERS----------------------///
    gsInfo << "Testing Eigen's interative solvers:\n";

    gsSparseSolver<>::CGIdentity EigenCGIsolver;
    EigenCGIsolver.setMaxIterations(maxIters);
    EigenCGIsolver.setTolerance(tol);
    gsInfo << "\nEigen's CG identity preconditioner: Started solving..."  << "\n";
    clock.restart();
    EigenCGIsolver.compute(mat);
    x0 = EigenCGIsolver.solve(rhs);
    gsInfo << "Eigen's CG: Tolerance           : " << EigenCGIsolver.tolerance() << "\n";
    gsInfo << "Eigen's CG: Residual error      : " << EigenCGIsolver.error() << "\n";
    gsInfo << "Eigen's CG: Number of iterations: " << EigenCGIsolver.iterations() << "\n";
    gsInfo << "Eigen's CG: Time to solve       : " << clock.stop() << "\n";


    gsSparseSolver<>::CGDiagonal EigenCGDsolver;
    EigenCGDsolver.setMaxIterations(maxIters);
    EigenCGDsolver.setTolerance(tol);
    gsInfo << "\nEigen's CG diagonal preconditioner: Started solving..."  << "\n";
    clock.restart();
    EigenCGDsolver.compute(mat);
    x0 = EigenCGDsolver.solve(rhs);
    gsInfo << "Eigen's CG: Tolerance           : " << EigenCGDsolver.tolerance() << "\n";
    gsInfo << "Eigen's CG: Residual error      : " << EigenCGDsolver.error() << "\n";
    gsInfo << "Eigen's CG: Number of iterations: " << EigenCGDsolver.iterations() << "\n";
    gsInfo << "Eigen's CG: Time to solve       : " << clock.stop() << "\n";

    gsSparseSolver<>::BiCGSTABIdentity EigenBCGIsolver;
    EigenBCGIsolver.setMaxIterations(maxIters);
    EigenBCGIsolver.setTolerance(tol);
    gsInfo << "\nEigen's bi conjugate gradient stabilized solver identity preconditioner: Started solving..."  << "\n";
    clock.restart();
    EigenBCGIsolver.compute(mat);
    x0 = EigenBCGIsolver.solve(rhs);
    gsInfo << "Eigen's BiCGSTAB: Tolerance           : " << EigenBCGIsolver.tolerance() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Residual error      : " << EigenBCGIsolver.error() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Number of iterations: " << EigenBCGIsolver.iterations() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Time to solve       : " << clock.stop() << "\n";

    gsSparseSolver<>::BiCGSTABDiagonal EigenBCGDsolver;
    EigenBCGDsolver.setMaxIterations(maxIters);
    EigenBCGDsolver.setTolerance(tol);
    gsInfo << "\nEigen's bi conjugate gradient stabilized solver diagonal preconditioner: Started solving..."  << "\n";
    clock.restart();
    EigenBCGDsolver.compute(mat);
    x0 = EigenBCGDsolver.solve(rhs);
    gsInfo << "Eigen's BiCGSTAB: Tolerance           : " << EigenBCGDsolver.tolerance() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Residual error      : " << EigenBCGDsolver.error() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Number of iterations: " << EigenBCGDsolver.iterations() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Time to solve       : " << clock.stop() << "\n";

    gsSparseSolver<>::BiCGSTABILUT EigenBCGILUsolver;
    //EigenBCGILUsolver.preconditioner().setFillfactor(1);
    EigenBCGILUsolver.setMaxIterations(maxIters);
    EigenBCGILUsolver.setTolerance(tol);
    gsInfo << "\nEigen's bi conjugate gradient stabilized solver ILU preconditioner: Started solving..."  << "\n";
    clock.restart();
    EigenBCGILUsolver.compute(mat);
    x0 = EigenBCGILUsolver.solve(rhs);
    gsInfo << "Eigen's BiCGSTAB: Tolerance           : " << EigenBCGILUsolver.tolerance() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Residual error      : " << EigenBCGILUsolver.error() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Number of iterations: " << EigenBCGILUsolver.iterations() << "\n";
    gsInfo << "Eigen's BiCGSTAB: Time to solve       : " << clock.stop() << "\n";


    ///----------------------EIGEN-DIRECT-SOLVERS----------------------///
    gsSparseSolver<>::SimplicialLDLT EigenSLDLTsolver;
    gsInfo << "\nEigen's Simplicial LDLT: Started solving..."  << "\n";
    clock.restart();
    EigenSLDLTsolver.compute(mat);
    x0 = EigenSLDLTsolver.solve(rhs);
    gsInfo << "Eigen's Simplicial LDLT: Time to solve       : " << clock.stop() << "\n";

    gsSparseSolver<>::QR solverQR;
    gsInfo << "\nEigen's QR: Started solving..."  << "\n";
    clock.restart();
    solverQR.compute(mat);
    x0 = solverQR.solve(rhs);
    gsInfo << "Eigen's QR: Time to solve       : " << clock.stop() << "\n";

#endif

    gsSparseSolver<>::LU solverLU;
    gsInfo << "\nEigen's LU: Started solving..."  << "\n";
    clock.restart();
    solverLU.compute(mat);
    x0 = solverLU.solve(rhs);
    gsInfo << "Eigen's LU: Time to solve       : " << clock.stop() << "\n";


    return 0;
}
