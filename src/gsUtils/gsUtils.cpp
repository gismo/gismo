/** @file gsUtils.cpp

 @brief Several utility functions for miscellaneous tasks

 This file is part of the G+Smo library.

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Author(s): A. Mantzaflaris, Harald Weiner
 */

#include <gsCore/gsForwardDeclarations.h>
#include <gsUtils/gsUtils.h>

#include <cstdlib>

#if defined(_WIN32)
#include <windows.h>
#include <direct.h>
#define getcwd _getcwd
#else
#include <dlfcn.h>
#include <unistd.h>
#endif

namespace gismo
{

namespace util
{

std::string getTempPath()
{
#       if   defined(_WIN32)
    TCHAR _temp[MAX_PATH];
    (void)GetTempPath(MAX_PATH, // length of the buffer
            _temp);// buffer for path
    return std::string(_temp);
#       else
    char * _temp;
#       if defined(__APPLE__)
    _temp = getenv ("TMPDIR");
#       elif defined(__unix)
    _temp = getenv("TEMP");
#       endif

    std::string path;
    if (_temp != NULL)
    {
        path = _temp;
        std::free(_temp);
        return path;
    }

    _temp = getcwd(NULL, 0);
    path = _temp;
    std::free(_temp);
    return path;
#       endif
}

} // end namespace util

} // end namespace gismo
