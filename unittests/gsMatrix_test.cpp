/** @file gsMatrix_test.cpp

    @brief Tests for gsMatrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Linus Mussmaecher
**/

#include "gismo_unittest.h"

SUITE(gsMatrix_test)
{

// -------------------------------------------------------------------------
// rotate_cw
// -------------------------------------------------------------------------

TEST(rotate_cw_square)
{
    // 3x3 input:
    gsMatrix<real_t> A(3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    // Expected 3x3 result:
    gsMatrix<real_t> expected(3, 3);
    expected << 3, 6, 9,
                2, 5, 8,
                1, 4, 7;

    gsMatrix<real_t> result = A.rotate_cw();

    CHECK_EQUAL(3, result.rows());
    CHECK_EQUAL(3, result.cols());
    CHECK(result.isApprox(expected));
}

TEST(rotate_cw_non_square)
{
    // 2x3 input:
    gsMatrix<real_t> A(2, 3);
    A << 1, 2, 3,
         4, 5, 6;

    // Expected 3x2 result:
    gsMatrix<real_t> expected(3, 2);
    expected << 3, 6,
                2, 5,
                1, 4;

    gsMatrix<real_t> result = A.rotate_cw();

    CHECK_EQUAL(3, result.rows());
    CHECK_EQUAL(2, result.cols());
    CHECK(result.isApprox(expected));
}

TEST(rotate_cw_four_times_is_identity)
{
    // 2x3 input:
    gsMatrix<real_t> A(2, 3);
    A << 1, 2, 3,
         4, 5, 6;

    gsMatrix<real_t> result = A.rotate_cw().rotate_cw().rotate_cw().rotate_cw();

    CHECK_EQUAL(A.rows(), result.rows());
    CHECK_EQUAL(A.cols(), result.cols());
    CHECK(result.isApprox(A));
}

// -------------------------------------------------------------------------
// rotate_ccw
// -------------------------------------------------------------------------

TEST(rotate_ccw_square)
{
    // 3x3 input:
    gsMatrix<real_t> A(3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;

    // Expected 3x3 result:
    gsMatrix<real_t> expected(3, 3);
    expected << 7, 4, 1,
                8, 5, 2,
                9, 6, 3;

    gsMatrix<real_t> result = A.rotate_ccw();

    CHECK_EQUAL(3, result.rows());
    CHECK_EQUAL(3, result.cols());
    CHECK(result.isApprox(expected));
}

TEST(rotate_ccw_non_square)
{
    // 2x3 input:
    gsMatrix<real_t> A(2, 3);
    A << 1, 2, 3,
         4, 5, 6;

    // Expected 3x2 result:
    gsMatrix<real_t> expected(3, 2);
    expected << 4, 1,
                5, 2,
                6, 3;

    gsMatrix<real_t> result = A.rotate_ccw();

    CHECK_EQUAL(3, result.rows());
    CHECK_EQUAL(2, result.cols());
    CHECK(result.isApprox(expected));
}

TEST(rotate_ccw_four_times_is_identity_non_square)
{
    // 2x3 input:
    gsMatrix<real_t> A(2, 3);
    A << 1, 2, 3,
         4, 5, 6;

    gsMatrix<real_t> result = A.rotate_ccw().rotate_ccw().rotate_ccw().rotate_ccw();

    CHECK_EQUAL(A.rows(), result.rows());
    CHECK_EQUAL(A.cols(), result.cols());
    CHECK(result.isApprox(A));
}

// -------------------------------------------------------------------------
// rotate_cw combination
// -------------------------------------------------------------------------

TEST(rotate_inverses)
{
    gsMatrix<real_t> A(2, 3);
    A << 1, 2, 3,
         4, 5, 6;

    CHECK(A.rotate_cw().rotate_ccw().isApprox(A));
    CHECK(A.rotate_ccw().rotate_cw().isApprox(A));
}

} // SUITE
