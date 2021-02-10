/** @file gsAbs_Diff_test.cpp

    @brief Tests math::abs_diff, which is part of CXX11 and
    for CXX98 defined in gsTemplateTools.h.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl
**/

#include "gismo_unittest.h"

#define __STDC_LIMIT_MACROS
#include <stdint.h>

// This Unittest file tests math::abs_diff, and therefore
// util::make_unsigned, which is native C++11 and has an
// C++98 conform custom implementation in gsCore/gsTemplateTools.h
SUITE(gsAbs_Diff_test)
{
    TEST(int32_m1_1)
    {
        int i0 = -1;
        int i1 = 1;
        unsigned int exp = 2;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_min_max)
    {
        int i0 = INT32_MIN;
        int i1 = INT32_MAX;
        unsigned int exp = UINT32_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_min_m1)
    {
        int i0 = INT32_MIN;
        int i1 = -1;
        unsigned int exp = INT32_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_min_0)
    {
        int i0 = INT32_MIN;
        int i1 = 0;
        unsigned int exp = INT32_MAX + 1U;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_min_1)
    {
        int i0 = INT32_MIN;
        int i1 = 1;
        unsigned int exp = INT32_MAX + 2U;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_m1_max)
    {
        int i0 = -1;
        int i1 = INT32_MAX;
        unsigned int exp = INT32_MAX + 1U;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_0_max)
    {
        int i0 = 0;
        int i1 = INT32_MAX;
        unsigned int exp = INT32_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_1_max)
    {
        int i0 = 1;
        int i1 = INT32_MAX;
        unsigned int exp = INT32_MAX - 1;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(int32_max_max)
    {
        int i0 = INT32_MAX;
        int i1 = INT32_MAX;
        unsigned int exp = 0;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(uint32_0_imax)
    {
        unsigned int i0 = 0;
        unsigned int i1 = INT32_MAX;
        unsigned int exp = INT32_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(uint32_imax_imax)
    {
        unsigned int i0 = INT32_MAX;
        unsigned int i1 = INT32_MAX;
        unsigned int exp = 0;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_m1_1)
    {
        long long i0 = -1;
        long long i1 = 1;
        unsigned long long exp = 2;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_min_max)
    {
        long long i0 = INT64_MIN;
        long long i1 = INT64_MAX;
        unsigned long long exp = UINT64_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_min_m1)
    {
        long long i0 = INT64_MIN;
        long long i1 = -1;
        unsigned long long exp = INT64_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

/*
    TEST(INT64_min_0)
    {
        long long i0 = INT64_MIN;
        long long i1 = 0;
        unsigned long long exp = INT64_MAX + 1UL; //msvc warning C4307, C4245

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_min_1)
    {
        long long i0 = INT64_MIN;
        long long i1 = 1;
        unsigned long long exp = INT64_MAX + 2UL; //msvc warning C4307, C4245

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_m1_max)
    {
        long long i0 = -1;
        long long i1 = INT64_MAX;
        unsigned long long exp = INT64_MAX + 1UL; //msvc warning C4307, C4245

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }
*/

    TEST(INT64_0_max)
    {
        long long i0 = 0;
        long long i1 = INT64_MAX;
        unsigned long long exp = INT64_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_1_max)
    {
        long long i0 = 1;
        long long i1 = INT64_MAX;
        unsigned long long exp = INT64_MAX - 1;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(INT64_max_max)
    {
        long long i0 = INT64_MAX;
        long long i1 = INT64_MAX;
        unsigned long long exp = 0;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(uINT64_0_imax)
    {
        unsigned long long i0 = 0;
        unsigned long long i1 = INT64_MAX;
        unsigned long long exp = INT64_MAX;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }

    TEST(uINT64_imax_imax)
    {
        unsigned long long i0 = INT64_MAX;
        unsigned long long i1 = INT64_MAX;
        unsigned long long exp = 0;

        CHECK(math::abs_diff(i0, i1) == exp);
        CHECK(math::abs_diff(i1, i0) == exp);
    }
}
