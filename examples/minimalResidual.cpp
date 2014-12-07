/** @file minimalResidual.cpp

    @brief Example on how the solve a system og linear equation with the min. res. method.

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

    //Initialize an Identity preconditioner
    gsIdentityPreconditioner preConMat(N);

    //Set maximun number of iterations
    index_t maxIters = 1000;
    //Set tolerance
    real_t tol = 1e-08;

    //Initialize the solver
    gsMinimalResidual<gsMatrix<> > MinRes(mat,maxIters,tol);

    //Create the initial guess
    gsVector<> x0;
    x0.setZero(N);

    //Solve system with given preconditioner (solution is stored in x0)
    MinRes.solve(rhs,x0,preConMat);

    gsInfo << "Solved a system of size " << N << "\n";
    gsInfo << "Tolerance: " << tol << "\n";
    gsInfo << "Residual error: " << MinRes.error() << "\n";
    gsInfo << "Number of iterations: " << MinRes.iterations() << "\n";

    return (MinRes.error()<tol)?0:1;
}
