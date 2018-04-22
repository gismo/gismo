/** @file gsKronecker_test.cpp

    @brief Tests the Kronecker operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris, S. Takacs
*/

#include "gismo_unittest.h"

gsMatrix<> A = (gsMatrix<>(3,3) <<
                3, 7, 3,
                2, -2, 1,
                9, 5, 2
    ).finished();

gsMatrix<> B = (gsMatrix<>(3,3) <<
                6, 1, 4,
                -3, 4, 2,
                1, 2, 5
    ).finished();

gsMatrix<> C = (gsMatrix<>(9,9) << // kron(A,B) computed in Matlab
    18,    3,   12,   42,    7,   28,   18,    3,   12,
    -9,    12,   6,  -21,   28,   14,   -9,   12,    6,
    3,     6,   15,    7,   14,   35,    3,    6,   15,
    12,    2,    8,  -12,   -2,   -8,    6,    1,    4,
    -6,    8,    4,    6,   -8,   -4,   -3,    4,    2,
    2,     4,   10,   -2,   -4,  -10,    1,    2,    5,
    54,    9,   36,   30,    5,   20,   12,    2,    8,
    -27,   36,  18,  -15,   20,   10,   -6,    8,    4,
    9,     18,  45,    5,   10,   25,    2,    4,   10
    ).finished();

gsMatrix<> x = (gsMatrix<>(9,1) << 1,2,3,4,5,6,7,8,9).finished();

SUITE(gsKronecker_test)
{
    TEST(testKroneckerOp)
    {
        gsKroneckerOp<> kron( makeMatrixOp(A), makeMatrixOp(B) );
        gsMatrix<>  y;
        kron.apply(x, y);
        // compute Kronecker product directly and compare
        CHECK_EQUAL ( y, C * x );
    }

    TEST(DenseKronecker)
    {        
        gsMatrix<> D = A.kron(B);
        CHECK_EQUAL(C, D);
    }

    TEST(DenseKhatriRao)
    {        
        gsMatrix<> D = A.kathriRao(B);
        CHECK_EQUAL(C, D);
    }

    TEST(SparseKronecker)
    {
        gsSparseMatrix<> sA = A.sparseView();
        gsSparseMatrix<> sB = B.sparseView();
        gsSparseMatrix<> sD = sA.kron(sB);
        CHECK_EQUAL ( C,  sD.toDense() );
    }
}
