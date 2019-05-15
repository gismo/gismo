/** @file gsFileManager_test.cpp

    @brief Tests for gsFileManager

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include "gismo_unittest.h"

SUITE(gsFileManager_test)
{
    TEST(getExePath)
    {
#if defined _WIN32
        std::string own_fn("unittests.exe");
#else
        std::string own_fn("unittests");
#endif
        CHECK( gsFileManager::fileExists( gsFileManager::getExePath() + own_fn ) )
    }

}
