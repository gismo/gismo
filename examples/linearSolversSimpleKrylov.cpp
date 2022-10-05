/** @file linearSolversDirect.cpp

    @brief Example on how the solve a system of linear equations with
    unpreconditioned iterative solvers.

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
void poisson1D(gsSparseMatrix<> &mat, index_t N)
{
    mat.resize(N,N);
    mat.setZero();

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

    //Compress the matrix
    mat.makeCompressed();
}

void poisson2D(gsSparseMatrix<> &mat, index_t n)
{
    gsSparseMatrix<> K;
    gsSparseMatrix<> M;
    poisson1D(K, n);
    M.resize(n,n);
    M.setIdentity();
    mat = K.kron(M) + M.kron(K);    
}

void poisson3D(gsSparseMatrix<> &mat, index_t n)
{
    gsSparseMatrix<> K;
    gsSparseMatrix<> M;
    poisson1D(K, n);
    M.resize(n,n);
    M.setIdentity();
    mat = K.kron(M).kron(M) + M.kron(K).kron(M) + M.kron(M).kron(K);
}

//Print out information of the iterative solver
template<typename SolverType>
void gsIterativeSolverInfo(const SolverType &method,
                           real_t error, double time, bool& succeeded )
{
        
    gsInfo << method.detail();
    gsInfo << " Computed residual error  : " << error << "\n";
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
    gsStopwatch clock;
    index_t n   = 50;
    index_t dim = 2;
    bool useCG     = false;
    bool useMinRes = false;
    bool useGMRes  = false;
    real_t tol = std::pow(10.0, - REAL_DIG * 0.75);

    
    gsCmdLine cmd ("Solves a PDE with a Courant discretization with several unpreconditioned iterative solvers.");
    cmd.addInt    ("n", "number" , "Number of unknowns", n        );
    cmd.addInt    ("d", "dim"    , "Spartial dimension", dim      );
    cmd.addSwitch ("" , "cg"    , "Use     CG solver",   useCG    );
    cmd.addSwitch ("" , "minres", "Use MinRes solver",   useMinRes);
    cmd.addSwitch ("" , "gmres" , "Use  GMRes solver",   useGMRes );
    cmd.addReal   ("",  "tol"    , "Tolerance for the iterative solvers", tol);
    
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    gsInfo << cmd <<"\n";
    gsSparseMatrix<> mat;
    gsMatrix<> rhs;
    gsMatrix<> x; 
    
    //Assemble the Poisson equation
    if (dim==1)
        poisson1D(mat, n);
    else if (dim==2)
        poisson2D(mat, n);
    else if (dim==3)        
        poisson3D(mat, n);
    else
        GISMO_ERROR("Chosen dimension in invalid, possible choices are d=1,2,3.");

    const index_t N = mat.rows();
    rhs.setRandom(N,1);

    gsInfo << "Solving a "<<dim<<"D Poisson problem. Matrix size is "<<mat.rows()<<"\n";


    //Initialize a preconditioner (Linear operator) 
    gsLinearOperator<>::Ptr IdentityOp = gsIdentityOp<>::make(mat.rows());
    
    if (useCG)
    {
        gsConjugateGradient<> CGSolver(mat,IdentityOp);
        CGSolver.setTolerance(tol);
        CGSolver.setMaxIterations(N*3);
        x.setZero(N,1);
        gsInfo << "\nCG: Started solving... ";
        clock.restart();
        CGSolver.solve(rhs,x);
        gsInfo << "done.\n";
        gsIterativeSolverInfo(CGSolver, (mat*x-rhs).norm()/rhs.norm(), clock.stop(), succeeded);
    }
    if (useMinRes)
    {
        gsMinimalResidual<> MinRes(mat,IdentityOp);
        MinRes.setTolerance(tol);
        MinRes.setMaxIterations(N*3);
        x.setZero(N,1);
        gsInfo << "\nMinRes: Started solving... ";
        clock.restart();
        MinRes.solve(rhs,x);
        gsInfo << "done.\n";
        gsIterativeSolverInfo(MinRes, (mat*x-rhs).norm()/rhs.norm(), clock.stop(), succeeded);
    }
    if(useGMRes)
    {
        gsGMRes<> GMRes(mat,IdentityOp);
        GMRes.setTolerance(tol);
        GMRes.setMaxIterations(N*3);
        x.setZero(N,1);
        gsInfo << "\nGMRes: Started solving... ";
        clock.restart();
        GMRes.solve(rhs,x);
        gsInfo << "done.\n";
        gsIterativeSolverInfo(GMRes, (mat*x-rhs).norm()/rhs.norm(), clock.stop(), succeeded);
    }
    
    return succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
}
