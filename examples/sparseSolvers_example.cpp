/** @file sparseSolvers_example.cpp

    @brief Testing the use of sparse linear solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

void report( const gsVector<>& computedSolution, const gsVector<>& exactSolution, bool& succeeded )
{
    gsInfo << "  Computed solution: " << computedSolution.transpose() << "\n";
    if ( (computedSolution-exactSolution).norm() <= 1.e-10 )
    {
        gsInfo << "  Test passed.\n";
    }
    else
    {
        gsInfo << "  Test faild.\n";
        succeeded = false;
    }
    gsInfo << "\n";
}

int main(int argc, char** argv)
{
#ifdef EIGEN_USE_MKL
        gsInfo << "EIGEN_USE_MKL=true.\n";
#endif

    index_t mat_size = 10;

    gsCmdLine cmd("Testing the use of sparse linear solvers.");

    cmd.addInt("n", "size", "Size of the matrices", mat_size);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsSparseMatrix<>  Q(mat_size,mat_size);
    gsVector<>        b(mat_size), x(mat_size), x0(mat_size);
    x0.setOnes();

    bool succeeded = true;

    Q.reserve( gsVector<int>::Constant(mat_size,1) ); // Reserve memory for 1 non-zero entry per column
    for (index_t i = 0; i!=mat_size; ++i)
        Q(i,i) = b[i] = i+1;

    Q.makeCompressed(); // always call makeCompressed after sparse matrix has been filled

//    /*
    gsSparseSolver<>::CGIdentity solverCGI;
    solverCGI.compute(Q);
    x = solverCGI.solve(b);
    gsInfo << "Solve Ax = b with Eigen's CG identity preconditioner.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::CGDiagonal solverCGD;
    solverCGD.compute(Q);
    x = solverCGD.solve(b);
    gsInfo << "Solve Ax = b with Eigen's CG diagonal preconditioner.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::BiCGSTABILUT solverBCGILU;
    solverBCGILU.compute(Q);
    x = solverBCGILU.solve(b);
    gsInfo << "Solve Ax = b with Eigen's BiCG with ILU preconditioner.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::BiCGSTABDiagonal solverBCGD;
    solverBCGD.compute(Q);
    x = solverBCGD.solve(b);
    gsInfo << "Solve Ax = b with Eigen's BiCG with diagonal preconditioner.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::BiCGSTABIdentity solverBCDI;
    solverBCDI.compute(Q);
    x = solverBCDI.solve(b);
    gsInfo << "Solve Ax = b with Eigen's BiCG without preconditioner.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::SimplicialLDLT solverSLDLT;
    solverSLDLT.compute(Q);
    x = solverSLDLT.solve(b);
    gsInfo << "Solve Ax = b with Eigen's Simplicial LDLT.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::QR solverQR;
    solverQR.compute(Q);
    x = solverQR.solve(b);
    gsInfo << "Solve Ax = b with Eigen's QR factorization.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::LU solverLU;
    solverLU.compute(Q);
    x = solverLU.solve(b);
    gsInfo << "Solve Ax = b with Eigen's LU factorization.\n";
    report( x, x0, succeeded );

//*/

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLU solverpLU;
    solverpLU.compute(Q);
    x = solverpLU.solve(b);
    gsInfo << "Error code of pardiso "<< solverpLU.info() <<"\n";
    gsInfo << "Solve Ax = b with PardisoLU.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::PardisoLDLT solverLDLT;
    solverLDLT.compute(Q);
    x = solverLDLT.solve(b);
    gsInfo << "Error code of pardiso "<< solverLDLT.info() <<"\n";
    gsInfo << "Solve Ax = b with PardisoLDLT.\n";
    report( x, x0, succeeded );

    gsSparseSolver<>::PardisoLLT solverLLT;
    solverLLT.compute(Q);
    x = solverLLT.solve(b);
    gsInfo << "Error code of pardiso "<< solverLLT.info() <<"\n";
    gsInfo << "Solve Ax = b with PardisoLLT.\n";
    report( x, x0, succeeded );
#   else
    gsInfo << "PARDISO is not available.\n";
#   endif

#ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solverSLU;
    solverSLU.compute(Q);
    x = solverSLU.solve(b);
    gsInfo << "Solve Ax = b with Super.\n";
    report( x, x0, succeeded );

#   else
    gsInfo << "SuperLU is not available.\n";
#   endif

    #ifdef GISMO_WITH_PASTIX
    gsInfo << "PastiX is not available.\n";
#   else
    gsInfo << "PastiX is not available.\n";
#   endif

    return succeeded ? 0 : 1;
}
