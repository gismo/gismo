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
    CHECK_EQUAL('\\', gsFileManager::getValidPathSeparators()[0]);
    CHECK_EQUAL('/', gsFileManager::getValidPathSeparators()[1]);
#else
    CHECK_EQUAL('/', gsFileManager::getValidPathSeparators()[0]);
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
    GISMO_ASSERT(gsFileManager::getSearchPaths() == "", "gsFileManager::getSearchPaths() not empty");

    std::string verum0 = gsFileManager::getExePath();
    std::string verum1 = gsFileManager::getTempPath();
    GISMO_ASSERT(verum0 != verum1, "gsFileManager::getExePath == gsFileManager::getTempPath()");
    std::string falsum("/fuubar");

    std::string result;

    CHECK(!gsFileManager::setSearchPaths(falsum));
    CHECK_EQUAL("", gsFileManager::getSearchPaths());

    CHECK(gsFileManager::setSearchPaths(verum0));
    result = gsFileManager::getSearchPaths();
    CHECK_EQUAL(verum0 + ";", result);
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), result[result.length() - 2]);

    CHECK(gsFileManager::addSearchPaths(verum1));
    result = gsFileManager::getSearchPaths();
    CHECK_EQUAL(verum0 + ";" + verum1 + ";", result);
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), result[result.length() - 2]);

    // clear SearchPaths
    CHECK(gsFileManager::setSearchPaths(""));
    CHECK_EQUAL("", gsFileManager::getSearchPaths());

    // set more values at once
    CHECK(gsFileManager::setSearchPaths(verum0 + ";" + verum1));
    result = gsFileManager::getSearchPaths();
    CHECK_EQUAL(verum0 + ";" + verum1 + ";", result);
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), result[result.length() - 2]);

    gsFileManager::setSearchPaths(defaultPath);
    GISMO_ASSERT(gsFileManager::getSearchPaths() == defaultPath, "Can't set back to default getSearchPaths!");
}

TEST(find)
{
    std::string defaultPath = gsFileManager::getSearchPaths();
    gsFileManager::setSearchPaths("");
    GISMO_ASSERT(gsFileManager::getSearchPaths() == "", "gsFileManager::getSearchPaths() not empty");

    // calculate relative path
    std::string absolute = GISMO_DATA_DIR;              // absolute
    std::string current = gsFileManager::getCurrentPath();
    std::string relative = gsFileManager::makeRelative(current, absolute);
    GISMO_ASSERT(gsFileManager::isExplicitlyRelative(relative),
        "variable relative isn't gsFilemanager::isExplicitlyRealtive");

    std::string verum("options/assembler_options.xml"); // success, if path known
    std::string falsum("fuubar");                       // fails

    // check without SearchPaths
    CHECK_EQUAL(relative + verum, gsFileManager::find(relative + verum));
    CHECK_EQUAL(absolute + verum, gsFileManager::find(absolute + verum));
    CHECK_EQUAL("", gsFileManager::find(verum));
    CHECK_EQUAL("", gsFileManager::find(falsum));

    // check with SearchPath set to relative
    gsFileManager::setSearchPaths(relative);
    CHECK_EQUAL(relative + verum, gsFileManager::find(relative + verum));
    CHECK_EQUAL(absolute + verum, gsFileManager::find(absolute + verum));
    CHECK_EQUAL(absolute + verum, gsFileManager::find(verum));
    CHECK_EQUAL("", gsFileManager::find(falsum));

    // clear SearchPaths
    gsFileManager::setSearchPaths("");
    GISMO_ASSERT(gsFileManager::getSearchPaths() == "", "gsFileManager::getSearchPaths() not empty");

    // check with SearchPath set to absolute
    gsFileManager::setSearchPaths(absolute);
    CHECK_EQUAL(relative + verum, gsFileManager::find(relative + verum));
    CHECK_EQUAL(absolute + verum, gsFileManager::find(absolute + verum));
    CHECK_EQUAL(absolute + verum, gsFileManager::find(verum));
    CHECK_EQUAL("", gsFileManager::find(falsum));

    gsFileManager::setSearchPaths(defaultPath);
    GISMO_ASSERT(gsFileManager::getSearchPaths() == defaultPath, "Can't set back to default getSearchPaths!");
}

TEST(fileExists)
{
    // FILE
    std::string defaultPath = gsFileManager::getSearchPaths();
    gsFileManager::setSearchPaths("");
    GISMO_ASSERT(gsFileManager::getSearchPaths() == "", "gsFileManager::getSearchPaths() not empty");

    std::string relative("./");                         // relative
    std::string absolute = gsFileManager::getExePath(); // absolute
    std::string falsum("fuubar");                       // fails

    CHECK(gsFileManager::fileExists(relative + own_fn));
    CHECK(gsFileManager::fileExists(absolute + own_fn));

    CHECK(!gsFileManager::fileExists(falsum));

    gsFileManager::setSearchPaths(relative);
    CHECK(gsFileManager::fileExists(own_fn));
    CHECK(!gsFileManager::fileExists(falsum));

    gsFileManager::setSearchPaths("");
    GISMO_ASSERT(gsFileManager::getSearchPaths() == "", "gsFileManager::getSearchPaths() not empty");

    gsFileManager::setSearchPaths(absolute);
    CHECK(gsFileManager::fileExists(own_fn));
    CHECK(!gsFileManager::fileExists(falsum));

    gsFileManager::setSearchPaths(defaultPath);
    GISMO_ASSERT(gsFileManager::getSearchPaths() == defaultPath, "Can't set back to default getSearchPaths!");
}

TEST(findInDataDir)
{
    std::string verum("options/assembler_options.xml");
    std::string falsum("fuu/bar");

    CHECK_EQUAL(GISMO_DATA_DIR + verum, gsFileManager::findInDataDir(verum));
    CHECK(gsFileManager::isFullyQualified(gsFileManager::findInDataDir(verum)));
    CHECK(!gsFileManager::isExplicitlyRelative(gsFileManager::findInDataDir(verum)));

    CHECK_EQUAL("", gsFileManager::findInDataDir(falsum));
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
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), testString[testString.length() - 1]);
}

TEST(getCurrentPath)
{
    std::string testString = gsFileManager::getCurrentPath();
    CHECK(testString != "");
    CHECK(gsFileManager::isFullyQualified(testString));
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), testString[testString.length() - 1]);
}

TEST(getExePath)
{
    std::string testString = gsFileManager::getExePath();
    CHECK(testString != "");
    CHECK(gsFileManager::isFullyQualified(testString));
    CHECK_EQUAL(gsFileManager::getNativePathSeparator(), testString[testString.length() - 1]);
    CHECK(gsFileManager::fileExists(testString + own_fn));
}

TEST(mkdir)
{
    std::string temp = gsFileManager::getTempPath();
    if (temp != "")
    {
        std::stringstream stream;
        for (int i = 0; i < 0xFFFF; ++i)
        {
            stream << temp << "gsMkDir" << std::hex << i << gsFileManager::getNativePathSeparator();
            if (!gsFileManager::fileExists(stream.str() + gsFileManager::getNativePathSeparator() + to_string(i) + ".xml"))
            {
                CHECK(gsFileManager::mkdir(stream.str()));  // create new directory
                CHECK(gsFileManager::mkdir(stream.str()));  // already existing directory
                stream << i;
                {
                    gsBSplineBasis<> geo;
                    gsWrite(geo, stream.str());
                }
                stream << ".xml";
                CHECK(gsFileManager::fileExists(stream.str()));
                CHECK(!gsFileManager::mkdir(stream.str())); // failing create of directory
                break;
            }
            stream.str("");
        }
    }
    else
        CHECK(false);
}

TEST(pathEqual)
{
    // absolute
    CHECK(gsFileManager::pathEqual("/", "/"));
    CHECK(gsFileManager::pathEqual("/foo", "/foo"));
    CHECK(gsFileManager::pathEqual("/foo/bar", "/foo/bar"));
    CHECK(!gsFileManager::pathEqual("/foo", "/bar"));
    CHECK(!gsFileManager::pathEqual("/foo/bar", "/bar/foo"));
    // canonical
    CHECK(gsFileManager::pathEqual("/foo/./bar", "/foo/bar"));
    CHECK(gsFileManager::pathEqual("/foo/baz/../bar", "/foo/bar"));
    // path separator - 2nd param
    CHECK(gsFileManager::pathEqual("/foo", "/foo/"));
    CHECK(gsFileManager::pathEqual("/foo/bar", "/foo/bar/"));
    CHECK(!gsFileManager::pathEqual("/foo", "/bar/"));
    CHECK(!gsFileManager::pathEqual("/foo/bar", "/bar/foo/"));
    CHECK(gsFileManager::pathEqual("/foo/./bar", "/foo/bar/"));
    CHECK(gsFileManager::pathEqual("/foo/baz/../bar", "/foo/bar/"));
    // path separator - 1st param
    CHECK(gsFileManager::pathEqual("/foo/", "/foo"));
    CHECK(gsFileManager::pathEqual("/foo/bar/", "/foo/bar"));
    CHECK(!gsFileManager::pathEqual("/foo/", "/bar"));
    CHECK(!gsFileManager::pathEqual("/foo/bar/", "/bar/foo"));
    CHECK(gsFileManager::pathEqual("/foo/./bar/", "/foo/bar"));
    CHECK(gsFileManager::pathEqual("/foo/baz/../bar/", "/foo/bar"));

    // relative
    CHECK(gsFileManager::pathEqual("", ""));
    CHECK(gsFileManager::pathEqual("", "./"));
    CHECK(gsFileManager::pathEqual("", gsFileManager::getCurrentPath()));
    CHECK(gsFileManager::pathEqual("foo", "foo"));
    CHECK(gsFileManager::pathEqual("foo/bar", "foo/bar"));
    CHECK(!gsFileManager::pathEqual("foo", "bar"));
    CHECK(!gsFileManager::pathEqual("foo/bar", "bar/foo"));
    // canonical
    CHECK(gsFileManager::pathEqual("foo/./bar", "foo/bar"));
    CHECK(gsFileManager::pathEqual("foo/baz/../bar", "foo/bar"));
    CHECK(gsFileManager::pathEqual("../bar", "./buz/../../bar"));
    // path separator - 2nd param
    CHECK(gsFileManager::pathEqual("foo", "foo/"));
    CHECK(gsFileManager::pathEqual("foo/bar", "foo/bar/"));
    CHECK(!gsFileManager::pathEqual("foo", "bar/"));
    CHECK(!gsFileManager::pathEqual("foo/bar", "bar/foo/"));
    CHECK(gsFileManager::pathEqual("foo/./bar", "foo/bar/"));
    CHECK(gsFileManager::pathEqual("foo/baz/../bar", "foo/bar/"));
    // path separator - 1st param
    CHECK(gsFileManager::pathEqual("foo/", "foo"));
    CHECK(gsFileManager::pathEqual("foo/bar/", "foo/bar"));
    CHECK(!gsFileManager::pathEqual("foo/", "bar"));
    CHECK(!gsFileManager::pathEqual("foo/bar/", "bar/foo"));
    CHECK(gsFileManager::pathEqual("foo/./bar/", "foo/bar"));
    CHECK(gsFileManager::pathEqual("foo/baz/../bar/", "foo/bar"));

#if defined _WIN32
    CHECK(gsFileManager::pathEqual("/foo/bar", "\\foo\\bar"));
    CHECK(gsFileManager::pathEqual("foo/bar", "foo\\bar"));
    CHECK(gsFileManager::pathEqual("../bar", ".\\buz\\..\\..\\bar"));
#endif
}

TEST(getExtension)
{
    CHECK_EQUAL("bar", gsFileManager::getExtension("foo.bar"));
    CHECK_EQUAL("bar", gsFileManager::getExtension("/foo.bar"));
    CHECK_EQUAL("bar", gsFileManager::getExtension("./foo.bar"));
    CHECK_EQUAL("bar", gsFileManager::getExtension("../foo.bar"));

    CHECK_EQUAL("bar", gsFileManager::getExtension("foo.baz.bar"));
    CHECK_EQUAL("bar", gsFileManager::getExtension("/foo.baz.bar"));
    CHECK_EQUAL("bar", gsFileManager::getExtension("./foo.baz.bar"));
    CHECK_EQUAL("bar", gsFileManager::getExtension("../foo.baz.bar"));

    CHECK_EQUAL("bar", gsFileManager::getExtension("/some/../baz/other/../foo.baz.bar"));

    CHECK_EQUAL("foo", gsFileManager::getExtension("bar.foo"));
}

TEST(getBasename)
{
    CHECK_EQUAL("foo", gsFileManager::getBasename("foo.bar"));
    CHECK_EQUAL("foo", gsFileManager::getBasename("/foo.bar"));
    CHECK_EQUAL("foo", gsFileManager::getBasename("./foo.bar"));
    CHECK_EQUAL("foo", gsFileManager::getBasename("../foo.bar"));

    CHECK_EQUAL("foo.baz", gsFileManager::getBasename("foo.baz.bar"));
    CHECK_EQUAL("foo.baz", gsFileManager::getBasename("/foo.baz.bar"));
    CHECK_EQUAL("foo.baz", gsFileManager::getBasename("./foo.baz.bar"));
    CHECK_EQUAL("foo.baz", gsFileManager::getBasename("../foo.baz.bar"));

    CHECK_EQUAL("foo.baz", gsFileManager::getBasename("/some/../bax/other/../foo.baz.bar"));

    CHECK_EQUAL("bar", gsFileManager::getBasename("bar.foo"));
}

TEST(getFilename)
{
    CHECK_EQUAL("foo", gsFileManager::getFilename("/foo"));
    CHECK_EQUAL("foo.bar", gsFileManager::getFilename("/foo.bar"));
    CHECK_EQUAL("foo.baz.bar", gsFileManager::getFilename("/foo.baz.bar"));
    CHECK_EQUAL("foo", gsFileManager::getFilename("../foo"));
    CHECK_EQUAL("foo.bar", gsFileManager::getFilename("../foo.bar"));
    CHECK_EQUAL("foo.baz.bar", gsFileManager::getFilename("../foo.baz.bar"));
    CHECK_EQUAL("foo", gsFileManager::getFilename("some/other/foo"));
    CHECK_EQUAL("foo.bar", gsFileManager::getFilename("some/other/foo.bar"));
    CHECK_EQUAL("foo.baz.bar", gsFileManager::getFilename("some/other/foo.baz.bar"));
    CHECK_EQUAL("foo", gsFileManager::getFilename("/some/other/foo"));
    CHECK_EQUAL("foo.bar", gsFileManager::getFilename("/some/other/foo.bar"));
    CHECK_EQUAL("foo.baz.bar", gsFileManager::getFilename("/some/other/foo.baz.bar"));
    CHECK_EQUAL("foo", gsFileManager::getFilename("./some/other/foo"));
    CHECK_EQUAL("foo.bar", gsFileManager::getFilename("./some/other/foo.bar"));
    CHECK_EQUAL("foo.baz.bar", gsFileManager::getFilename("./some/other/foo.baz.bar"));
}

TEST(getPath)
{
    CHECK_EQUAL("", gsFileManager::getPath("foo"));
    CHECK_EQUAL("", gsFileManager::getPath("foo.bar"));
    CHECK_EQUAL("", gsFileManager::getPath("foo.baz.bar"));
    CHECK_EQUAL("/", gsFileManager::getPath("/foo"));
    CHECK_EQUAL("/", gsFileManager::getPath("/foo.bar"));
    CHECK_EQUAL("/", gsFileManager::getPath("/foo.baz.bar"));
    CHECK_EQUAL("./", gsFileManager::getPath("./foo"));
    CHECK_EQUAL("./", gsFileManager::getPath("./foo.bar"));
    CHECK_EQUAL("./", gsFileManager::getPath("./foo.baz.bar"));
    CHECK_EQUAL("../", gsFileManager::getPath("../foo"));
    CHECK_EQUAL("../", gsFileManager::getPath("../foo.bar"));
    CHECK_EQUAL("../", gsFileManager::getPath("../foo.baz.bar"));
    CHECK_EQUAL("some/other/", gsFileManager::getPath("some/other/foo"));
    CHECK_EQUAL("some/other/", gsFileManager::getPath("some/other/foo.bar"));
    CHECK_EQUAL("some/other/", gsFileManager::getPath("some/other/foo.baz.bar"));
    CHECK_EQUAL("/some/other/", gsFileManager::getPath("/some/other/foo"));
    CHECK_EQUAL("/some/other/", gsFileManager::getPath("/some/other/foo.bar"));
    CHECK_EQUAL("/some/other/", gsFileManager::getPath("/some/other/foo.baz.bar"));
    CHECK_EQUAL("./some/other/", gsFileManager::getPath("./some/other/foo"));
    CHECK_EQUAL("./some/other/", gsFileManager::getPath("./some/other/foo.bar"));
    CHECK_EQUAL("./some/other/", gsFileManager::getPath("./some/other/foo.baz.bar"));
    CHECK_EQUAL("../some/other/", gsFileManager::getPath("../some/other/foo"));
    CHECK_EQUAL("../some/other/", gsFileManager::getPath("../some/other/foo.bar"));
    CHECK_EQUAL("../some/other/", gsFileManager::getPath("../some/other/foo.baz.bar"));
    // some canonical check
    CHECK_EQUAL("../some/other/", gsFileManager::getPath("../some/buz/../other/foo"));
    CHECK_EQUAL("../some/other/", gsFileManager::getPath("../some/buz/../other/foo.bar"));
    CHECK_EQUAL("../some/other/", gsFileManager::getPath("../some/buz/../other/foo.baz.bar"));
}

TEST(getCanonicRepresentation)
{

    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("foo/./bar"));
    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("foo/baz/../bar"));
    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("foo/baz/baz/../../bar"));
    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("foo/baz/../baz/../bar"));
    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("foo/baz/baz/.././././../bar"));
    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("foo/baz/./.././baz/./.././bar"));
    CHECK_EQUAL("foo/bar", gsFileManager::getCanonicRepresentation("baz/baz/../../foo/bar"));

    CHECK_EQUAL("./foo", gsFileManager::getCanonicRepresentation("./././foo"));
    CHECK_EQUAL("../foo", gsFileManager::getCanonicRepresentation("./.././foo"));

    CHECK_EQUAL("", gsFileManager::getCanonicRepresentation(""));
    CHECK_EQUAL("./", gsFileManager::getCanonicRepresentation("", true));
    CHECK_EQUAL("./foo/", gsFileManager::getCanonicRepresentation("./././foo", true));
    CHECK_EQUAL("../foo/", gsFileManager::getCanonicRepresentation("./.././foo", true));

    CHECK_EQUAL("./foo/", gsFileManager::getCanonicRepresentation("./././foo/"));
    CHECK_EQUAL("../foo/", gsFileManager::getCanonicRepresentation("./.././foo/"));
}

}
