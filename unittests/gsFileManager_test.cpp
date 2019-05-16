/** @file gsFileManager_test.cpp

    @brief Tests for gsFileManager

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, J. Vogl
*/

#include "gismo_unittest.h"

SUITE(gsFileManager_test)
{
TEST(PathSeperators)
{
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), gsFileManager::getValidPathSeparators()[0]);

#if defined _WIN32
    CHECK_EQUAL(gsFileManager::getValidPathSeparators()[0], '\\');
    CHECK_EQUAL(gsFileManager::getValidPathSeparators()[1], '/');
#else
    CHECK_EQUAL(gsFileManager::getValidPathSeparators()[0], '/');
#endif
}

TEST(isFullyQualified)
{
    // for any OS
    std::string verum0("/");
    std::string verum1("/foo");
    std::string verum2("/Foo/Bar");
    std::string verum3("/foo bar/bar foo");

    CHECK(gsFileManager::isFullyQualified(verum0));
    CHECK(gsFileManager::isFullyQualified(verum1));
    CHECK(gsFileManager::isFullyQualified(verum2));
    CHECK(gsFileManager::isFullyQualified(verum3));

    std::string falsum0("");
    std::string falsum1("foo");
    std::string falsum2("foo\\bar");
    std::string falsum3("\"" + verum2 + "\"");

    CHECK(!gsFileManager::isFullyQualified(falsum0));
    CHECK(!gsFileManager::isFullyQualified(falsum1));
    CHECK(!gsFileManager::isFullyQualified(falsum2));
    CHECK(!gsFileManager::isFullyQualified(falsum3));

    // OS specific
#if defined _WIN32
    std::string verum4("\\");
    std::string verum5("E:\\Foo\\Bar");
    std::string verum6("\\foo bar\\bar foo");
    std::string verum7("f:/foo bar/bar foo");

    CHECK(gsFileManager::isFullyQualified(verum6));
    CHECK(gsFileManager::isFullyQualified(verum7));

    std::string falsum4("\\c:\\Foo\\Bar");
    std::string falsum5("/c:/Foo/Bar");
#else
    std::string verum4("\\c:\\Foo\\Bar");
    std::string verum5("/c:/Foo/Bar");

    std::string falsum4("\\");
    std::string falsum5("E:\\Foo\\Bar");
    std::string falsum6("\\foo bar\\bar foo");
    std::string falsum7("f:/foo bar/bar foo");

    CHECK(!gsFileManager::isFullyQualified(falsum6));
    CHECK(!gsFileManager::isFullyQualified(falsum7));
#endif

    CHECK(gsFileManager::isFullyQualified(verum4));
    CHECK(gsFileManager::isFullyQualified(verum5));

    CHECK(!gsFileManager::isFullyQualified(falsum4));
    CHECK(!gsFileManager::isFullyQualified(falsum5));
}

TEST(isExplicitlyRelative) {
    
}

TEST(getExePath)
{
#if defined _WIN32 // || defined __CYGWIN__
    std::string own_fn("unittests.exe");
#else
    std::string own_fn("unittests");
#endif
    CHECK(gsFileManager::fileExists(gsFileManager::getExePath() + own_fn));
}

TEST(SearchPaths)
{
    std::string verum0 = gsFileManager::getExePath();
    std::string verum1 = gsFileManager::getTempPath();
    CHECK_ASSERT(verum0 != verum1);
    std::string falsum("/fuubar");

    gsFileManager::setSearchPaths(falsum);
    CHECK_EQUAL(gsFileManager::getSearchPaths(), "");

    gsFileManager::setSearchPaths(verum1);
    CHECK_EQUAL(gsFileManager::getSearchPaths(), verum1);

    gsFileManager::addSearchPaths(verum1);
    CHECK_EQUAL(gsFileManager::getSearchPaths(), verum0 + ";" + verum1);
}

TEST(find)
{
#if defined _WIN32
    std::string path("C:\\Windows");
    std::string file("notepad.exe");
#else
    std::string path("/bin");
    std::string file("mkdir");
#endif
    std::string falsum("fuubar");

    gsFileManager::setSearchPaths(path);
    CHECK(gsFileManager::find(file) == (gsFileManager::getCanonicRepresentation(path)
        + gsFileManager::getNativePathSeparator() + file));

    CHECK(gsFileManager::find(falsum) == "");
}

}
