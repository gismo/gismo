/** @file iterativeSolvers.cpp

    @brief Example on how the solve a system of linear equation with the MINRES and CG method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>
#include <gismo.h>

using std::cout;
using std::endl;
using namespace gismo;

int main(int argc, char *argv[])
{
    //Size of linear system
    index_t N = 200;
    if (argc >= 2)
        N = atoi(argv[1]);

    gsMatrix<> mat;
    gsMatrix<> rhs;

    mat.setZero(N,N);
    rhs.setRandom(N,1);

    //Create a tri-diagonal matrix with -1 of the off diagonals and 2 in the diagonal.
    //This matrix is equivalent to discretizing the 1D Poisson equation with homogenius
    //Dirichlet boundary condition using a finit difference schema. It is a SPD matrix.
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

    //The minimal residual implementation requires a preconditioner.
    //We initialize an Identity preconditioner (does nothing).
    gsIdentityPreconditioner preConMat(N);

    //Set maximum number of iterations
    index_t maxIters = 1000;
    //Set tolerance 
    real_t tol = math::pow(10.0, - REAL_DIG / 3);

    //Initialize the MinRes solver
    gsMinimalResidual MinRes(mat,maxIters,tol);

    //Create the initial guess
    gsMatrix<> x0;
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nMinRes: Before solve"  << std::endl;
    MinRes.solve(rhs,x0,preConMat);
    gsInfo << "MinRes: After solve"  << std::endl;

    gsInfo << "MinRes: Solved a system of size " << N << "\n";
    gsInfo << "MinRes: Tolerance: " << tol << "\n";
    gsInfo << "MinRes: Residual error: " << MinRes.error() << "\n";
    gsInfo << "MinRes: Number of iterations: " << MinRes.iterations() << "\n";

    //Initialize the CG solver
    gsGMRes GMResSolver(mat,maxIters,tol);

    //Set the initial guess to zero
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nGMRes: Before solve"  << std::endl;
    GMResSolver.solve(rhs,x0,preConMat);
    gsInfo << "GMRes: After solve"  << std::endl;

    gsInfo << "GMRes: Solved a system of size " << N << "\n";
    gsInfo << "GMRes: Tolerance: " << tol << "\n";
    gsInfo << "GMRes: Residual error: " << GMResSolver.error() << "\n";
    gsInfo << "GMRes: Number of iterations: " << GMResSolver.iterations() << "\n";

    //Initialize the CG solver
    gsConjugateGradient CGSolver(mat,maxIters,tol);

    //Set the initial guess to zero
    x0.setZero(N,1);

    //Solve system with given preconditioner (solution is stored in x0)
    gsInfo << "\nCG: Before solve"  << std::endl;
    CGSolver.solve(rhs,x0,preConMat);
    gsInfo << "CG: After solve"  << std::endl;

    gsInfo << "CG: Solved a system of size " << N << "\n";
    gsInfo << "CG: Tolerance: " << tol << "\n";
    gsInfo << "CG: Residual error: " << CGSolver.error() << "\n";
    gsInfo << "CG: Number of iterations: " << CGSolver.iterations() << "\n";

    int result = (MinRes.error()<tol)?0:1;
    result += (CGSolver.error()<tol)?0:1;
    return result;
}
