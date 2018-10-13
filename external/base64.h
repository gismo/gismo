/** @file base64.h

    @brief base-64 encoding for file I/O

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <string>
#include <cstring>//for memcpy

#include <gsCore/gsExport.h>
#include <gsCore/gsDebug.h>

std::string base64_encode(unsigned char const*, unsigned int);

std::string base64_decode(const char*);

//template<class C> struct podType
//{ static const std::string name = ""; };

template<class C> inline
std::string encode64_array(const C* arr, unsigned int sz)
{
    return base64_encode(reinterpret_cast<const unsigned char*>(arr.data()),sz*sizeof(C));
}

template<class C> inline
void decode64_array(const char* arr, C* out)
{
    // POD types: type="double"
    // integer, float (mpq?)
    //if (sizeof(C)!=dsize) read+cast POD to C
    // todo.
    //else
    // check POD
    std::string dc = base64_decode(arr); // todo: decode in-place in out
    GISMO_ENSURE(dc.size()%sizeof(C)==0, // ==dsize
                 "Decoding seems to have failed, mod="<< dc.size()%sizeof(C));
    std::memcpy(reinterpret_cast<char*>(out),dc.c_str(),dc.size());
}
