/** @file gsFileManager_test.cpp

    @brief Tests for gsFileManager

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, J. Vogl
*/

#include "gismo_unittest.h"
#include "../../../../../../usr/lib/gcc/x86_64-pc-linux-gnu/9.1.0/include/g++-v9/bits/basic_string.h"

SUITE(gsFileManager_test)
{
    TEST(getExePath)
    {
#if defined _WIN32 // || defined __CYGWIN__
        std::string own_fn("unittests.exe");
#else
        std::string own_fn("unittests");
#endif
        CHECK( gsFileManager::fileExists( gsFileManager::getExePath() + own_fn ) );
    }

    TEST(isFullyQualified)
    {
#if defined _WIN32
        std::string verum0("C:\\Fuu\\Bar");
        std::string verum1("c:\\fuu\\bar");
        std::string verum2("\"c:\\fuu bar\\bar fuu\"");
        std::string verum3("c:/fuu bar/bar fuu");

        std::string falsum0("\\c:\\Fuu\\Bar");
        std::string falsum1("");
        std::string falsum2("fuu\\bar");
        std::string falsum3("\\fuu\\bar");
#else
        std::string verum0("/Fuu/Bar");
        std::string verum1("/fuu bar/bar fu");
        std::string verum2("/c:/fuu bar/bar fuu");
        std::string verum3("/");

        std::string falsum0("C:\\Fuu\\Bar");
        std::string falsum1("");
        std::string falsum2("\"c:\\fuu bar\\bar fuu\"");
        std::string falsum3("c:\\fuu bar\\bar fuu");
#endif
        CHECK(gsFileManager::isFullyQualified(verum0));
        CHECK(gsFileManager::isFullyQualified(verum1));
        CHECK(gsFileManager::isFullyQualified(verum2));
        CHECK(gsFileManager::isFullyQualified(verum3));

        CHECK(!gsFileManager::isFullyQualified(falsum0));
        CHECK(!gsFileManager::isFullyQualified(falsum1));
        CHECK(!gsFileManager::isFullyQualified(falsum2));
        CHECK(!gsFileManager::isFullyQualified(falsum3));
    }

    TEST(SearchPaths)
    {
#if defined _WIN32
        std::string verum("C:\\Windows"); // Assert, this exists
        std::string falsum("C:\\fuubar"); // Assert, this doesn't exists
#else
        std::string verum("/boot");
        std::string falsum("/fuubar");
#endif
        gsFileManager::setSearchPaths(falsum);
        CHECK(gsFileManager::getSearchPaths() == "");
        gsFileManager::setSearchPaths(verum);
        CHECK(gsFileManager::getSearchPaths().find(verum) != std::string::npos);
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
