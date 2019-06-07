/** @file gsUtils_test.cpp

    @brief Tests helper functions in util namespace.
    Most of them exist native in CXX11 and are handmade for CXX98.
    This Unittests help us to prove that the work as expected.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl
**/

#include "gismo_unittest.h"

SUITE(gsUtils)
{
TEST(to_string)
{
    CHECK_EQUAL("1", util::to_string(int(1)));
    CHECK_EQUAL("-1", util::to_string(int(-1)));
    CHECK_EQUAL("1", util::to_string(long(1)));
    CHECK_EQUAL("abcd", util::to_string("abcd"));
    CHECK_EQUAL("abcd", util::to_string("abcd\0"));
    CHECK_EQUAL("ab", util::to_string("ab\0cd"));
    char array[] = { 'a', 'b', 'c' };
    CHECK_EQUAL("abc", util::to_string(array));
}

TEST(starts_with)
{
    CHECK(util::starts_with("abcd", "a"));
    CHECK(util::starts_with("abcd", "ab"));
    CHECK(util::starts_with("abcd", "abc"));
    CHECK(util::starts_with("abcd", "abcd"));

    CHECK(!util::starts_with("abcd", "bcd"));
    CHECK(!util::starts_with("abcd", "acd"));
    CHECK(!util::starts_with("abcd", "abd"));

    CHECK(!util::starts_with("abc", "abcd"));

    CHECK(util::starts_with("a\0bcd", "a"));
    CHECK(!util::starts_with("a\0bcd", "ab"));
    CHECK(util::starts_with("acd", "a\0b"));
}
TEST(ends_with)
{
    CHECK(util::ends_with("abcd", "d"));
    CHECK(util::ends_with("abcd", "cd"));
    CHECK(util::ends_with("abcd", "bcd"));
    CHECK(util::ends_with("abcd", "abcd"));

    CHECK(!util::ends_with("abcd", "acd"));
    CHECK(!util::ends_with("abcd", "abd"));
    CHECK(!util::ends_with("abcd", "abc"));

    CHECK(!util::ends_with("abc", "abcd"));

    CHECK(util::ends_with("abc\0d", "bc"));
    CHECK(util::ends_with("a\0bcd", "a"));
    CHECK(util::ends_with("abc", "bc\0d"));
}
TEST(iota)
{
    {
        std::vector<index_t> v(10);
        util::iota(v.begin(), v.end(), -4);

        for (size_t i = 0; i < v.size(); ++i)
        {
            CHECK_EQUAL((index_t)i - 4, v[i]);
        }
    }

    {
        std::list<index_t> l(10);
        util::iota(l.begin(), l.end(), -4);

        for (size_t i = 0; i < l.size(); ++i)
        {
            CHECK_EQUAL((index_t)i - 4, l.front());
            l.pop_front();
        }

    }
}

TEST(stoi)
{
    CHECK_EQUAL(1, util::stoi("0001"));
    CHECK_EQUAL(1234567890, util::stoi("1234567890"));
    CHECK_EQUAL(-1, util::stoi("-0001"));
    CHECK_EQUAL(-1234567890, util::stoi("-1234567890"));
    CHECK_EQUAL(42, util::stoi("42"));
    CHECK_EQUAL(42, util::stoi(" 42"));
    CHECK_EQUAL(42, util::stoi("42 "));
    CHECK_EQUAL(42, util::stoi(" 42 "));
    CHECK_EQUAL(4, util::stoi(" 4 2 "));
}

}