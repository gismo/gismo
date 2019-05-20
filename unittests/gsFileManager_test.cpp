/** @file gsFileManager_test.cpp

    @brief Tests for gsFileManager

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, J. Vogl
*/

#include "gismo_unittest.h"

#if defined _WIN32 || defined __CYGWIN__
#define own_fn "unittests.exe"
#else
#define own_fn "unittests"
#endif

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

TEST(Paths_absolut_relative)
{
    // for any OS
    std::string s00("/");
    std::string s01("/foo");
    std::string s02("/Foo/Bar");
    std::string s03("/foo bar/bar foo");
    std::string s04("");
    std::string s05("foo");
    std::string s06("foo\\bar");
    std::string s07("\"" + s02 + "\"");
    std::string s08("./");
    std::string s09("./foo");
    std::string s10("../");
    std::string s11("../foo");

    // isFullyQualified verum
    CHECK(gsFileManager::isFullyQualified(s00));
    CHECK(gsFileManager::isFullyQualified(s01));
    CHECK(gsFileManager::isFullyQualified(s02));
    CHECK(gsFileManager::isFullyQualified(s03));
    // isFullyQualified falsum
    CHECK(!gsFileManager::isFullyQualified(s04));
    CHECK(!gsFileManager::isFullyQualified(s05));
    CHECK(!gsFileManager::isFullyQualified(s06));
    CHECK(!gsFileManager::isFullyQualified(s07));
    CHECK(!gsFileManager::isFullyQualified(s08));
    CHECK(!gsFileManager::isFullyQualified(s09));
    CHECK(!gsFileManager::isFullyQualified(s10));
    CHECK(!gsFileManager::isFullyQualified(s11));

    // isExplicitlyRelative verum
    CHECK(gsFileManager::isExplicitlyRelative(s08));
    CHECK(gsFileManager::isExplicitlyRelative(s09));
    CHECK(gsFileManager::isExplicitlyRelative(s10));
    CHECK(gsFileManager::isExplicitlyRelative(s11));
    // isExplicitlyRelative falsum
    CHECK(!gsFileManager::isExplicitlyRelative(s00));
    CHECK(!gsFileManager::isExplicitlyRelative(s01));
    CHECK(!gsFileManager::isExplicitlyRelative(s02));
    CHECK(!gsFileManager::isExplicitlyRelative(s03));
    CHECK(!gsFileManager::isExplicitlyRelative(s04));
    CHECK(!gsFileManager::isExplicitlyRelative(s05));
    CHECK(!gsFileManager::isExplicitlyRelative(s06));
    CHECK(!gsFileManager::isExplicitlyRelative(s07));

    // OS specific
    std::string c00("\\");
    std::string c01("E:\\Foo\\Bar");
    std::string c02("\\foo bar\\bar foo");
    std::string c03("f:/foo bar/bar foo");
    std::string c04("\\c:\\Foo\\Bar");
    std::string c05("/c:/Foo/Bar");
    std::string c06(".\\");
    std::string c07(".\\foo");
    std::string c08("..\\");
    std::string c09("..\\foo");

#if defined _WIN32
    // isFullyQualified verum
    CHECK(gsFileManager::isFullyQualified(c00));
    CHECK(gsFileManager::isFullyQualified(c01));
    CHECK(gsFileManager::isFullyQualified(c02));
    CHECK(gsFileManager::isFullyQualified(c03));
    // isFullyQualified falsum
    CHECK(!gsFileManager::isFullyQualified(c04));
    CHECK(!gsFileManager::isFullyQualified(c05));
    CHECK(!gsFileManager::isFullyQualified(c06));
    CHECK(!gsFileManager::isFullyQualified(c07));
    CHECK(!gsFileManager::isFullyQualified(c08));
    CHECK(!gsFileManager::isFullyQualified(c09));

    // isExplicitlyRelative verum
    CHECK(gsFileManager::isExplicitlyRelative(c06));
    CHECK(gsFileManager::isExplicitlyRelative(c07));
    CHECK(gsFileManager::isExplicitlyRelative(c08));
    CHECK(gsFileManager::isExplicitlyRelative(c09));
    // isExplicitlyRelative falsum
    CHECK(!gsFileManager::isExplicitlyRelative(c00));
    CHECK(!gsFileManager::isExplicitlyRelative(c01));
    CHECK(!gsFileManager::isExplicitlyRelative(c02));
    CHECK(!gsFileManager::isExplicitlyRelative(c03));
    CHECK(!gsFileManager::isExplicitlyRelative(c04));
    CHECK(!gsFileManager::isExplicitlyRelative(c05));
#else

    // isFullyQualified verum
    CHECK(gsFileManager::isFullyQualified(c05));
    // isFullyQualified falsum
    CHECK(!gsFileManager::isFullyQualified(c00));
    CHECK(!gsFileManager::isFullyQualified(c01));
    CHECK(!gsFileManager::isFullyQualified(c02));
    CHECK(!gsFileManager::isFullyQualified(c03));
    CHECK(!gsFileManager::isFullyQualified(c04));
    CHECK(!gsFileManager::isFullyQualified(c06));
    CHECK(!gsFileManager::isFullyQualified(c07));
    CHECK(!gsFileManager::isFullyQualified(c08));
    CHECK(!gsFileManager::isFullyQualified(c09));

    // isExplicitlyRelative verum
    // isExplicitlyRelative falsum
    CHECK(!gsFileManager::isExplicitlyRelative(c00));
    CHECK(!gsFileManager::isExplicitlyRelative(c01));
    CHECK(!gsFileManager::isExplicitlyRelative(c02));
    CHECK(!gsFileManager::isExplicitlyRelative(c03));
    CHECK(!gsFileManager::isExplicitlyRelative(c04));
    CHECK(!gsFileManager::isExplicitlyRelative(c05));
    CHECK(!gsFileManager::isExplicitlyRelative(c06));
    CHECK(!gsFileManager::isExplicitlyRelative(c07));
    CHECK(!gsFileManager::isExplicitlyRelative(c08));
    CHECK(!gsFileManager::isExplicitlyRelative(c09));
#endif
}

TEST(SearchPaths)
{
    std::string defaultPath = gsFileManager::getSearchPaths();
    gsFileManager::setSearchPaths("");
    CHECK_ASSERT(gsFileManager::getSearchPaths() == "");

    std::string verum0 = gsFileManager::getExePath();
    std::string verum1 = gsFileManager::getTempPath();
    CHECK_ASSERT(verum0 != verum1);
    std::string falsum("/fuubar");

    std::string result;

    CHECK(!gsFileManager::setSearchPaths(falsum));
    CHECK_EQUAL(gsFileManager::getSearchPaths(), "");

    CHECK(gsFileManager::setSearchPaths(verum0));
    result = gsFileManager::getSearchPaths();
    CHECK_EQUAL(result, verum0);
    CHECK_EQUAL(result[result.length()-1], gsFileManager::getNativePathSeparator());

    CHECK(gsFileManager::addSearchPaths(verum1));
    result = gsFileManager::getSearchPaths();
    CHECK_EQUAL(result, verum0 + ";" + verum1);
    CHECK_EQUAL(result[result.length()-1], gsFileManager::getNativePathSeparator());

    // clear SearchPaths
    CHECK(!gsFileManager::setSearchPaths(""));
    CHECK_EQUAL(gsFileManager::getSearchPaths(), "");

    // set more values at once
    CHECK(gsFileManager::setSearchPaths(verum0 + ";" + verum1));
    result = gsFileManager::getSearchPaths();
    CHECK_EQUAL(result, verum0 + ";" + verum1);
    CHECK_EQUAL(result[result.length()-1], gsFileManager::getNativePathSeparator());

    gsFileManager::setSearchPaths(defaultPath);
    CHECK_ASSERT(gsFileManager::getSearchPaths() == defaultPath);
}

TEST(find)
{
    std::string defaultPath = gsFileManager::getSearchPaths();
    gsFileManager::setSearchPaths("");
    CHECK_ASSERT(gsFileManager::getSearchPaths() == "");

    std::string relative("./");                         // relative
    std::string absolute = gsFileManager::getExePath(); // absolute
    std::string falsum("fuubar");                       // fails

    CHECK_EQUAL(gsFileManager::find(relative + own_fn), relative + own_fn);
    CHECK_EQUAL(gsFileManager::find(absolute + own_fn), absolute + own_fn);

    CHECK_EQUAL(gsFileManager::find(own_fn), "");
    CHECK_EQUAL(gsFileManager::find(falsum), "");

    gsFileManager::setSearchPaths(relative);
    CHECK_EQUAL(gsFileManager::find(own_fn), absolute + own_fn);
    CHECK_EQUAL(gsFileManager::find(falsum), "");

    gsFileManager::setSearchPaths("");
    CHECK_ASSERT(gsFileManager::getSearchPaths() == "");

    gsFileManager::setSearchPaths(absolute);
    CHECK_EQUAL(gsFileManager::find(own_fn), absolute + own_fn);
    CHECK_EQUAL(gsFileManager::find(falsum), "");

    gsFileManager::setSearchPaths(defaultPath);
    CHECK_ASSERT(gsFileManager::getSearchPaths() == defaultPath);
}

TEST(fileExists)
{
    // FILE
    std::string defaultPath = gsFileManager::getSearchPaths();
    gsFileManager::setSearchPaths("");
    CHECK_ASSERT(gsFileManager::getSearchPaths() == "");

    std::string relative("./");                         // relative
    std::string absolute = gsFileManager::getExePath(); // absolute
    std::string falsum("fuubar");                       // fails

    CHECK(gsFileManager::fileExists(relative + own_fn));
    CHECK(gsFileManager::fileExists(absolute + own_fn));

    CHECK(!gsFileManager::fileExists(own_fn));
    CHECK(!gsFileManager::fileExists(falsum));

    gsFileManager::setSearchPaths(relative);
    CHECK(gsFileManager::fileExists(own_fn));
    CHECK(!gsFileManager::fileExists(falsum));

    gsFileManager::setSearchPaths("");
    CHECK_ASSERT(gsFileManager::getSearchPaths() == "");

    gsFileManager::setSearchPaths(absolute);
    CHECK(gsFileManager::fileExists(own_fn));
    CHECK(!gsFileManager::fileExists(falsum));

    gsFileManager::setSearchPaths(defaultPath);
    CHECK_ASSERT(gsFileManager::getSearchPaths() == defaultPath);
}

TEST(findInDataDir)
{
    std::string verum("options/assembler_options.xml");
    std::string falsum("fuu/bar");

    CHECK_EQUAL(gsFileManager::findInDataDir(verum), GISMO_DATA_DIR + verum);
    CHECK(gsFileManager::isFullyQualified(gsFileManager::findInDataDir(verum)));
    CHECK(!gsFileManager::isExplicitlyRelative(gsFileManager::findInDataDir(verum)));

    CHECK_EQUAL(gsFileManager::findInDataDir(falsum), "");
}

TEST(fileExistsInDataDir)
{
    std::string verum("options/assembler_options.xml");
    std::string falsum("fuu/bar");

    CHECK(gsFileManager::fileExistsInDataDir(verum));

    CHECK(!gsFileManager::fileExistsInDataDir(falsum));
}

TEST(getTempPath)
{
    std::string testString = gsFileManager::getTempPath();
    CHECK(testString != "");
    CHECK(gsFileManager::isFullyQualified(testString));
    CHECK(gsFileManager::fileExists(testString));
    CHECK_EQUAL(testString[testString.length()-1], gsFileManager::getNativePathSeparator());
}

TEST(getCurrentPath)
{
    std::string testString = gsFileManager::getTempPath();
    CHECK(testString != "");
    CHECK(gsFileManager::isFullyQualified(testString));
    CHECK(gsFileManager::fileExists(testString));
    CHECK_EQUAL(testString[testString.length()-1], gsFileManager::getNativePathSeparator());
}

TEST(getExePath)
{
    std::string testString = gsFileManager::getTempPath();
    CHECK(testString != "");
    CHECK(gsFileManager::isFullyQualified(testString));
    CHECK(gsFileManager::fileExists(testString));
    CHECK_EQUAL(testString[testString.length()-1], gsFileManager::getNativePathSeparator());
    CHECK(gsFileManager::fileExists(testString + own_fn));
}

TEST(mkdir)
{
    std::string temp = gsFileManager::getTempPath();
    if(temp != "" && gsFileManager::fileExists(temp)
       && temp[temp.length() - 1] == gsFileManager::getNativePathSeparator())
    {
        std::stringstream stream;
        for (int i = 0; i < 0xFFFF; ++i)
        {
            stream << temp << "gsMkDir" << std::hex << i;
            if (!gsFileManager::fileExists(stream.str()))
            {
                gsFileManager::mkdir(stream.str());
                CHECK(gsFileManager::fileExists(stream.str()));
                CHECK_EQUAL(stream.str() + gsFileManager::getNativePathSeparator(),
                    gsFileManager::find(stream.str()));
                break;
            }
            stream.str("");
        }
    }
    else
        CHECK(false);
}



}
