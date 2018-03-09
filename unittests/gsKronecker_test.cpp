/** @file gsKronecker_test.cpp

    @brief Tests the Kronecker operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include "gismo_unittest.h"

SUITE(gsKronecker_test)
{

    TEST(testKroneckerOp)
    {

        gsMatrix<index_t> A(3,3), B(3,3);

        A << 3, 7, 3,
            2, -2, 1,
            9, 5, 2;

        B << 6, 1, 4,
            -3, 4, 2,
            1, 2, 5;

        gsMatrix<index_t> C(9,9);
        C <<    // kron(A,B) computed in Matlab
            18,    3,   12,   42,    7,   28,   18,    3,   12,
            -9,    12,   6,  -21,   28,   14,   -9,   12,    6,
            3,     6,   15,    7,   14,   35,    3,    6,   15,
            12,    2,    8,  -12,   -2,   -8,    6,    1,    4,
            -6,    8,    4,    6,   -8,   -4,   -3,    4,    2,
            2,     4,   10,   -2,   -4,  -10,    1,    2,    5,
            54,    9,   36,   30,    5,   20,   12,    2,    8,
            -27,   36,  18,  -15,   20,   10,   -6,    8,    4,
            9,     18,  45,    5,   10,   25,    2,    4,   10;

        gsMatrix<index_t> x(9,1), y, y2;
        x << 1,2,3,4,5,6,7,8,9;

        gsKroneckerOp<> kron( makeMatrixOp(A), makeMatrixOp(B) );
        kron.apply(x, y);

        // compute Kronecker product directly and compare
        y2 = C * x;

        CHECK_EQUAL ( y2,  y );
    }

}
