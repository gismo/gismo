/** @file gsUtils.h

    @brief Several utility functions for miscellaneous tasks

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Author(s): A. Mantzaflaris, Harald Weiner
*/

#pragma once

#include <sstream>

#include <gsCore/gsExport.h>

#define STRINGIGY(x) #x

namespace gismo
{

/** @namespace util
    
    @brief This namespace gathers several utility functions for
    miscellaneous tasks
*/
namespace util
{

#if __cplusplus > 199711L
using std::to_string;
using std::iota;
using std::stod;
using std::stoi;
#else
/// \brief Converts value to string
/// \ingroup Utils
template<typename C>
std::string to_string(const C & value)
{
    std::ostringstream convert;
    convert << value;
    return convert.str();
}

// Fills the range [first, last) with sequentially increasing values,
// starting with value and rep//etitively evaluating ++value.             
template<class ForwardIterator, class T>
void iota(ForwardIterator first, ForwardIterator last, T value)
{
    while(first != last) {
        *first++ = value;
        ++value;
    }
}

inline int stoi(const std::string& str)
{
    std::istringstream ss(str);
    int i;
    if (!(ss >> std::noskipws >> i))
        //Extracting an int failed
        return 0;

    char c;
    if (ss >> c)
        //There was something after the number
        return 0;
    
    return i;
}

inline double stod(const std::string& str)
{
    std::istringstream ss(str);
    double i;
    if (!(ss >> std::noskipws >> i))
        //Extracting double failed
        return 0;

    char c;
    if (ss >> c)
        //There was something after the number
        return 0;

    return i;
}

#endif

/// \brief Replaces appeareances of \a oldStr with \a newStr inside the
/// string \a str
/// \ingroup Utils
inline void string_replace(std::string& str, 
                           const std::string& oldStr, 
                           const std::string& newStr)
{
    size_t pos = 0;
    while((pos = str.find(oldStr, pos)) != std::string::npos)
    {
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
}

/// \brief Auto-detect temp directory
/// \ingroup Utils
GISMO_EXPORT
std::string getTempPath();

} // end namespace util

} // end namespace gismo

