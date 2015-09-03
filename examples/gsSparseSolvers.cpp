/** @file gsPardiso_test.cpp

    @brief Testing the use of sparse linear solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

int main()
{
    index_t mat_size = 10;

    gsSparseMatrix<>  Q(mat_size,mat_size);
    gsVector<>        b(mat_size), x(mat_size);

    Q.reserve( gsVector<int>::Constant(mat_size,1) ); // Reserve memory for 1 non-zero entry per column
    Q(0,0) = 1;
    Q(1,1) = 2;
    Q(2,2) = 3;
    Q(3,3) = 4;
    Q(4,4) = 5;
    Q(5,5) = 6;
    Q(6,6) = 7;
    Q(7,7) = 8;
    Q(8,8) = 9;
    Q(9,9) = 10;
    Q.makeCompressed(); // always call makeCompressed after sparse matrix has been filled

    // fill-in right-hand side
    b << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10; 
   
    gsSparseSolver<>::CGIdentity solverCGI;
    solverCGI.compute(Q);
    x = solverCGI.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's CG identity preconditioner. Solution: "<< x.transpose() <<"\n";

    gsSparseSolver<>::CGDiagonal solverCGD;
    solverCGD.compute(Q);
    x = solverCGD.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's CG diagonal preconditioner. Solution: "<< x.transpose() <<"\n";
        
    gsSparseSolver<>::BiCGSTABILUT solverBCGILU;
    solverBCGILU.compute(Q);
    x = solverBCGILU.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's BiCG with ILU preconditioner. Solution: "<< x.transpose() <<"\n";
        
    gsSparseSolver<>::BiCGSTABDiagonal solverBCGD;
    solverBCGD.compute(Q);
    x = solverBCGD.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's BiCG with diagonal preconditioner. Solution: "<< x.transpose() <<"\n";
        
    gsSparseSolver<>::BiCGSTABIdentity solverBCDI;
    solverBCDI.compute(Q);
    x = solverBCDI.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's BiCG without preconditioner. Solution: "<< x.transpose() <<"\n";

    gsSparseSolver<>::SimplicialLDLT solverSLDLT;
    solverSLDLT.compute(Q);
    x = solverSLDLT.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's Simplicial LDLT. Solution: "<< x.transpose() <<"\n";

    gsSparseSolver<>::QR solverQR;
    solverQR.compute(Q);
    x = solverQR.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's QR factorization. Solution: "<< x.transpose() <<"\n";

    gsSparseSolver<>::LU solverLU;
    solverLU.compute(Q);
    x = solverLU.solve(b);
    gsInfo <<"Solved Ax = b with Eigen's LU factorization. Solution: "<< x.transpose() <<"\n";

    #ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solverLDLT;
    solverLDLT.compute(Q);
    x = solverLDLT.solve(b);
    gsInfo <<"Solved Ax = b with PardisoLDLT. Solution: "<< x.transpose() <<"\n";

    gsSparseSolver<>::PardisoLLT solverLLT;
    solverLLT.compute(Q);
    x = solverLLT.solve(b);
    gsInfo <<"Solved Ax = b with PardisoLLT. Solution: "<< x.transpose() <<"\n";    

    gsSparseSolver<>::PardisoLU solverpLU;
    solverpLU.compute(Q);
    gsInfo <<"Solved Ax = b with PardisoLU. Solution: "<< x.transpose() <<"\n";
    x = solverpLU.solve(b);
#   else
    gsInfo <<"PARDISO is not available.\n";
#   endif

    #ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solverSLU;
    solverSLU.compute(Q);
    gsInfo <<"Solved Ax = b with Super. Solution: "<< x.transpose() <<"\n";
    x = solverSLU.solve(b);
#   else
    gsInfo <<"SuperLU is not available.\n";
#   endif

    #ifdef GISMO_WITH_PASTIX
    gsInfo <<"PastiX is not available.\n";
#   else
    gsInfo <<"PastiX is not available.\n";
#   endif

    return 0;
}
