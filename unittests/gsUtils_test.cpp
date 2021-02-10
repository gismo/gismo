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

SUITE(gsUtils_test)
{
TEST(to_string)
{
    CHECK_EQUAL("1", util::to_string(int(1)));
    CHECK_EQUAL("-1", util::to_string(int(-1)));
    CHECK_EQUAL("1", util::to_string(long(1)));
    CHECK_EQUAL("abcd", util::to_string("abcd"));
    CHECK_EQUAL("abcd", util::to_string("abcd\0"));
    CHECK_EQUAL("ab", util::to_string("ab\0cd"));
    //char array[] = { 'a', 'b', 'c' };
    //CHECK_EQUAL("abc", util::to_string(array)); // should fail to build
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
    // based on CXX11+ std::stoi(str) <=> std::stoi(str, 0, 10);
    CHECK_EQUAL(1, util::stoi("0001"));
    CHECK_EQUAL(1234567890, util::stoi("1234567890"));
    CHECK_EQUAL(-1, util::stoi("-0001"));
    CHECK_EQUAL(-1234567890, util::stoi("-1234567890"));
    CHECK_EQUAL(42, util::stoi("42"));
    CHECK_EQUAL(42, util::stoi(" 42"));
    CHECK_EQUAL(42, util::stoi("42 "));
    CHECK_EQUAL(42, util::stoi(" 42 "));
    CHECK_EQUAL(4, util::stoi(" 4 2 "));

    CHECK_THROW(util::stoi("a1"), std::invalid_argument);
    CHECK_EQUAL(1, util::stoi("1a"));
    CHECK_EQUAL(0, util::stoi("0x42"));
}

TEST(stod)
{
    CHECK_EQUAL(0.1, util::stod(".1"));
    CHECK_EQUAL(1.1, util::stod("0001.1000"));
    CHECK_EQUAL(1234567890., util::stod("1234567890"));
    CHECK_EQUAL(1234567890.123456789, util::stod("1234567890.123456789"));
    CHECK_EQUAL(-1234567890.123456789, util::stod(" -1234567890.123456789 "));
    CHECK_EQUAL(4.2, util::stod("4.2"));
    CHECK_EQUAL(4.2, util::stod(" 4.2"));
    CHECK_EQUAL(4.2, util::stod("4.2 "));
    CHECK_EQUAL(4.2, util::stod(" 4.2 "));
    CHECK_EQUAL(4., util::stod(" 4 2 "));

    CHECK_THROW(util::stod("a0.5"), std::invalid_argument);
    CHECK_EQUAL(0.4, util::stod("0.4")); //Note: util::stod("0.4a") fails on C98
#if !defined(_MSC_VER) || _MSC_VER >= 1900
    CHECK_EQUAL(-255.99609375, util::stod("-0xFF.FF"));
    CHECK_EQUAL(-255.99609375, util::stod("-0XFF.FF"));
    CHECK_EQUAL(-255.99609375, util::stod("-0xff.ff"));
    CHECK_EQUAL(-255.99609375, util::stod("-0Xff.ff"));
#endif
}

TEST(string_replace)
{
    std::string s0 = "abcdefghijk";
    util::string_replace(s0, "abc", "xxx");
    CHECK_EQUAL("xxxdefghijk", s0);

    util::string_replace(s0, "ijk", "yyy");
    CHECK_EQUAL("xxxdefghyyy", s0);

    util::string_replace(s0, "efg", "zzz");
    CHECK_EQUAL("xxxdzzzhyyy", s0);

    util::string_replace(s0, "xxx", "x");
    util::string_replace(s0, "yyy", "y");
    util::string_replace(s0, "zzz", "z");

    CHECK_EQUAL("xdzhy", s0);

    util::string_replace(s0, "x", "abc");
    util::string_replace(s0, "y", "ijk");
    util::string_replace(s0, "z", "efg");

    CHECK_EQUAL("abcdefghijk", s0);
}

TEST(tokenize)
{
    std::string s0("ab,ba,aa,bb");
    CHECK_EQUAL("ab", util::tokenize(s0, ",", 0));
    CHECK_EQUAL("ba", util::tokenize(s0, ",", 1));
    CHECK_EQUAL("aa", util::tokenize(s0, ",", 2));
    CHECK_EQUAL("bb", util::tokenize(s0, ",", 3));
    CHECK_THROW(util::tokenize(s0, ",", 4), std::runtime_error);
    CHECK_EQUAL("ab,ba,aa,bb", s0);

    std::string s1("ab..ba");
    CHECK_EQUAL("ab", util::tokenize(s1, ".", 0));
    CHECK_EQUAL("ba", util::tokenize(s1, ".", 1));
    CHECK_THROW(util::tokenize(s1, ".", 2), std::runtime_error);

    // ab ... ba
    //    ^^^
    std::string s2("ab...ba");
    CHECK_EQUAL("ab", util::tokenize(s2, "..", 0));
    CHECK_EQUAL("ba", util::tokenize(s2, "..", 1));
    CHECK_THROW(util::tokenize(s2, "..", 2), std::runtime_error);

    // xxx ab c ba zzz
    //     ^^   ^^
    std::string s3("xxxabcbazzz");
    CHECK_EQUAL("xxx", util::tokenize(s3, "ab", 0));
    CHECK_EQUAL("c", util::tokenize(s3, "ab", 1));
    CHECK_EQUAL("zzz", util::tokenize(s3, "ab", 2));
    CHECK_THROW(util::tokenize(s3, "ab", 3), std::runtime_error);

    /// xxx abbab zzz
    ///     ^^^^^
    std::string s4("xxxabbabzzz");
    CHECK_EQUAL("xxx", util::tokenize(s4, "ab", 0));
    CHECK_EQUAL("zzz", util::tokenize(s4, "ab", 1));
    CHECK_THROW(util::tokenize(s4, "ab", 3), std::runtime_error);

    /// xxx aba b aba zzz
    std::string s5("xxxabababazzz");
    CHECK_EQUAL("xxx", util::tokenize(s5, "aba", 0));
    CHECK_EQUAL("zzz", util::tokenize(s5, "aba", 1));
    CHECK_THROW(util::tokenize(s5, "aba", 2), std::runtime_error);
}

TEST(capitalize)
{
    std::string s0("abc");
    util::capitalize(s0);
    CHECK_EQUAL("Abc", s0);

    s0 = "Abc";
    util::capitalize(s0);
    CHECK_EQUAL("Abc", s0);

    s0 = "aBc";
    util::capitalize(s0);
    CHECK_EQUAL("ABc", s0);

    s0 = "abC";
    util::capitalize(s0);
    CHECK_EQUAL("AbC", s0);

    s0 = "AbC";
    util::capitalize(s0);
    CHECK_EQUAL("AbC", s0);

    s0 = "ABc";
    util::capitalize(s0);
    CHECK_EQUAL("ABc", s0);

    s0 = "aBC";
    util::capitalize(s0);
    CHECK_EQUAL("ABC", s0);

    s0 = "ABC";
    util::capitalize(s0);
    CHECK_EQUAL("ABC", s0);

    s0 = "def";
    util::capitalize(s0);
    CHECK_EQUAL("Def", s0);
}

TEST(returnCapitalized)
{
    std::string s0("abc");
    CHECK_EQUAL("Abc", util::returnCapitalized(s0));

    s0 = "Abc";
    CHECK_EQUAL("Abc", util::returnCapitalized(s0));

    s0 = "aBc";
    CHECK_EQUAL("ABc", util::returnCapitalized(s0));

    s0 = "abC";
    CHECK_EQUAL("AbC", util::returnCapitalized(s0));

    s0 = "AbC";
    CHECK_EQUAL("AbC", util::returnCapitalized(s0));

    s0 = "ABc";
    CHECK_EQUAL("ABc", util::returnCapitalized(s0));

    s0 = "aBC";
    CHECK_EQUAL("ABC", util::returnCapitalized(s0));

    s0 = "ABC";
    CHECK_EQUAL("ABC", util::returnCapitalized(s0));

    s0 = "def";
    CHECK_EQUAL("Def", util::returnCapitalized(s0));
}

TEST(name)
{
    // Fundamental types
    CHECK_EQUAL("bool", util::type<bool>::name());
    CHECK_EQUAL("signed char", util::type<signed char>::name());
    CHECK_EQUAL("unsigned char", util::type<unsigned char>::name());
    CHECK_EQUAL("char", util::type<char>::name());
    CHECK_EQUAL("short", util::type<signed short>::name());
    CHECK_EQUAL("short", util::type<signed short>::name());
    CHECK_EQUAL("unsigned short", util::type<unsigned short int>::name());
    CHECK_EQUAL("int", util::type<signed int>::name());
    CHECK_EQUAL("unsigned int", util::type<unsigned>::name());
    CHECK_EQUAL("long", util::type<long>::name());
    CHECK_EQUAL("unsigned long", util::type<unsigned long>::name());
#ifndef _WIN32
    CHECK_EQUAL("long long", util::type<signed long long int>::name());
    CHECK_EQUAL("unsigned long long", util::type<unsigned long long int>::name());
    //CHECK_EQUAL("gismo::gsGenericGeometry<(" + util::type<short_t>::name() + ")2, " + util::type<real_t>::name() + ">", util::type<gsGenericGeometry<2> >::name());
#else
    CHECK_EQUAL("__int64", util::type<signed long long int>::name());
    CHECK_EQUAL("unsigned __int64", util::type<unsigned long long int>::name());
    CHECK_EQUAL("class gismo::gsGenericGeometry<2," + util::type<real_t>::name() + ">",
        util::type<gsGenericGeometry<2> >::name());
#endif
}
}
