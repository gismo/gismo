/** @file linearSolversDirect.cpp

    @brief Example on how the solve a system of linear equations with
    Eigen's direct solvers.

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

int main(int argc, char *argv[])
{
    
    bool succeeded = true;

    gsStopwatch clock;
    index_t n   = 50;
    index_t dim = 2;
    bool useQR      = false;
    bool usePardiso = false;
    bool noEigen    = false;
    index_t maxThreads = -1;


    gsCmdLine cmd("Solves a PDE with a Courant discretization with several solvers.");
    cmd.addInt    ("n", "number"     , "Number of unknowns",      n         );
    cmd.addInt    ("d", "dim"        , "Spartial dimension",      dim       );
    cmd.addSwitch ("" , "useQR"      , "Use QR solver as well",   useQR     );
    cmd.addSwitch ("" , "usePardiso" , "Use Pardiso solvers",     usePardiso);
    cmd.addSwitch ("" , "noEigen"    , "Don't use Eigen solvers", noEigen   );
    cmd.addInt    ("p", "maxThreads" , "max threads for pardiso", maxThreads);

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
    
    rhs.setZero(mat.rows(),1);    

    gsInfo << "Solving a "<<dim<<"D Poisson problem. Matrix size is "<<mat.rows()<<"\n";

    
    ///----------------------EIGEN-DIRECT-SOLVERS----------------------///
    if (!noEigen)
    {
        gsInfo << "\nEigen's Simplicial LDLT: Started solving... ";
        clock.restart();
        {
            gsSparseSolver<>::SimplicialLDLT EigenSLDLTsolver;    
            EigenSLDLTsolver.compute(mat);
            x = EigenSLDLTsolver.solve(rhs);
        }
        gsInfo << "done.\n";
        gsInfo << "Eigen's Simplicial LDLT: Time to solve       : " << clock.stop() << "\n";

        if (useQR)
        {
            gsInfo << "\nEigen's QR: Started solving... ";
            clock.restart();
            gsSparseSolver<>::QR solverQR;
            solverQR.compute(mat);
            gsInfo << "done.\n";
            gsInfo << "Eigen's QR: Time to solve       : " << clock.stop() << "\n";
            x = solverQR.solve(rhs);
        }
    
        gsInfo << "\nEigen's LU: Started solving... ";
        clock.restart();
        {
            gsSparseSolver<>::LU solverLU;
            solverLU.compute(mat);
            x = solverLU.solve(rhs);
        }
        gsInfo << "done.\n";
        gsInfo << "Eigen's LU: Time to solve       : " << clock.stop() << "\n";
    }

    ///----------------------PARDISO-DIRECT-SOLVERS----------------------///
    
#   ifdef GISMO_WITH_PARDISO
    if (maxThreads>0)    
        omp_set_num_threads(maxThreads);
    
    if (usePardiso)
    {
            gsInfo << "\nPardiso's LDLT: Started solving... ";
            clock.restart();
            {
                gsSparseSolver<>::PardisoLDLT solverPLDLT;
                solverPLDLT.compute(mat);
                x = solverPLDLT.solve(rhs);
            }
            gsInfo << "done.\n";
            gsInfo << "Pardiso's LDLT: Time to solve       : " << clock.stop() << "\n";
            
            gsInfo << "\nPardiso's LU: Started solving... ";
            clock.restart();
            {
                gsSparseSolver<>::PardisoLU solverPLU;
                solverPLU.compute(mat);
                x = solverPLU.solve(rhs);
            }
            gsInfo << "done.\n";
            gsInfo << "Pardiso's LU: Time to solve       : " << clock.stop() << "\n";
    }
    //PardisoLDLT
    //PardisoLLT
    //PardisoLU
#   endif
    
    return succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
}
