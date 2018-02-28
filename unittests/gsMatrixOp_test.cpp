/** @file gsMatrixOp_test.cpp

    @brief Tests for gsMatrixOp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include "gismo_unittest.h"

SUITE(gsMatrixOp_test)
{

    TEST(DenseMatrix)
    {
        gsMatrix<> A (3,3);
        A << 2,2,3,  4,5,6,  7,8,10;

        gsLinearOperator<>::Ptr Aop = makeMatrixOp(A);

        A(0,0) = 1; // check that gsMatrixOp holds no copy

        gsMatrix<> C;
        Aop->toMatrix(C);

        CHECK( ( A - C ).norm() <= 1.e-10 );
    }

    TEST(DenseMatrixTransposed)
    {
        gsMatrix<> A (3,3);
        A << 2,2,3,  4,5,6,  7,8,10;

        gsLinearOperator<>::Ptr Aop = makeMatrixOp(A.transpose());

        A(0,0) = 1; // check that gsMatrixOp holds no copy

        gsMatrix<> C;
        Aop->toMatrix(C);

        CHECK( ( A.transpose() - C ).norm() <= 1.e-10 );
    }

    TEST(DenseMatrixSymm2)
    {
        gsMatrix<> A (3,3);
        A << 2,2,3,  4,5,6,  7,8,10;
        gsMatrix<> B (3,3);
        B << 1,4,7,  4,5,8,  7,8,10;

        gsMatrix<>::Nested C = A.selfadjointView<Lower>().derived()
        
        A(0,0) = 1; // check that C holds no copy

        CHECK( ( B - C ).norm() <= 1.e-10 );
    }

    TEST(DenseMatrixSymm)
    {
        gsMatrix<> A (3,3);
        A << 2,2,3,  4,5,6,  7,8,10;
        gsMatrix<> B (3,3);
        B << 1,4,7,  4,5,8,  7,8,10;

        gsLinearOperator<>::Ptr Aop = makeMatrixOp(A.selfadjointView<Lower>());

        A(0,0) = 1; // check that gsMatrixOp holds no copy

        gsMatrix<> C;
        Aop->toMatrix(C);

        CHECK( ( B - C ).norm() <= 1.e-10 );
    }

}
