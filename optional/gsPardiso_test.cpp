/** @file gsPardiso_test.cpp

    @brief Testing the use of PARDISO paraller solver.

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

    gsSparseMatrix<double>  Q(mat_size,mat_size);
    gsVector<> b(mat_size), x(mat_size);
    Q.reserve(gsVector<int>::Constant(mat_size,6));
    
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
    b << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    Eigen::PardisoLDLT<gsSparseMatrix<> > solverLDLT;
    solverLDLT.compute(Q);
    x = solverLDLT.solve(b);
    gsInfo <<"Solved Ax = b with LDLT. Solution: "<< x.transpose() <<"\n";

    Eigen::PardisoLLT<gsSparseMatrix<> > solverLLT;
    solverLLT.compute(Q);
    x = solverLLT.solve(b);
    gsInfo <<"Solved Ax = b with LLT. Solution: "<< x.transpose() <<"\n";    

    Eigen::PardisoLU<gsSparseMatrix<> > solverLU;
    solverLU.compute(Q);
    gsInfo <<"Solved Ax = b with LU. Solution: "<< x.transpose() <<"\n";
    x = solverLU.solve(b);


}
