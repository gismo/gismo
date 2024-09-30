/** @file gsKronecker_test.cpp

    @brief Tests the Kronecker operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris, S. Takacs
*/

#include "gismo_unittest.h"

SUITE(gsKronecker_test)
{

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

    gsMatrix<> KP = (gsMatrix<>(9,9) << // kron(A,B) computed in Matlab
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

    gsMatrix<> KRP = (gsMatrix<>(9,3) <<
                      18,   7,  12,
                      -9,  28,   6,
                      3 ,  14,  15,
                      12,  -2,   4,
                      -6,  -8,   2,
                      2 ,  -4,   5,
                      54,   5,   8,
                      -27, 20,   4,
                      9 ,  10,  10
        ).finished();

    TEST(gsKroneckerOp)
    {
        gsKroneckerOp<> kron( makeMatrixOp(A), makeMatrixOp(B) );
        gsMatrix<> y, x = (gsMatrix<>(9,1) << 1,2,3,4,5,6,7,8,9).finished();
        kron.apply(x, y);
        // compute Kronecker product directly and compare
        CHECK_EQUAL ( y, KP * x );
    }

    TEST(DenseKronecker)
    {        
        gsMatrix<> C = A.kron(B);
        CHECK_EQUAL(KP, C);
    }

    TEST(DenseKhatriRao)
    {
        gsMatrix<> C = A.khatriRao(B);
        CHECK_EQUAL(KRP, C);
    }

    TEST(SparseKronecker)
    {
        gsSparseMatrix<> sA = A.sparseView();
        gsSparseMatrix<> sB = B.sparseView();
        gsSparseMatrix<> sC = sA.kron(sB);
        CHECK_EQUAL ( KP, sC.toDense() );
    }

    TEST(SparseKroneckerRowMajor)
    {
        gsSparseMatrix<real_t,RowMajor> sA = A.sparseView();
        gsSparseMatrix<real_t,RowMajor> sB = B.sparseView();
        gsSparseMatrix<real_t,RowMajor> sC = sA.kron(sB);
        CHECK_EQUAL ( KP, sC.toDense() );
    }

    TEST(SparseKroneckerEmpty)
    {
        gsSparseMatrix<> sA = A.sparseView();
        gsSparseMatrix<> sB;
        gsSparseMatrix<> sC = sA.kron(sB);
        CHECK_EQUAL ( sB.toDense(), sC.toDense() );
    }

}
