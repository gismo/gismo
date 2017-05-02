/** @file gsUtils.cpp

 @brief Several utility functions for miscellaneous tasks

 This file is part of the G+Smo library.

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Author(s): A. Mantzaflaris, H. Weiner, S. Takacs
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

    // Typically, we should consider TMPDIR
    //   http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap08.html#tag_08_03
    //   https://en.wikipedia.org/wiki/TMPDIR&oldid=728654758
    char * _temp = getenv ("TMPDIR");
    // getenv returns NULL ptr if the variable is unknown (http://en.cppreference.com/w/cpp/utility/program/getenv).
    // If it is an empty string, we should also exclude it.
    if (_temp != NULL && _temp[0] != '\0')
    {
        // note: env variable needs no free
        return std::string(_temp);
    }

    // Okey, if first choice did not work, try this:
    _temp = getenv("TEMP");
    if (_temp != NULL && _temp[0] != '\0')
    {
        // note: env variable needs no free
        return std::string(_temp);
    }
    
    // And as third choice, use just current directory
    // http://man7.org/linux/man-pages/man2/getcwd.2.html
    _temp = getcwd(NULL, 0);
    GISMO_ASSERT(NULL!=_temp, "getcwd returned NULL.");
    std::string path(_temp);
    // The string is allocated using malloc, see the reference above
    std::free(_temp);
    return path;
#       endif
}

} // end namespace util

} // end namespace gismo
