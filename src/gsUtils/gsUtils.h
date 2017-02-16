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
#include <numeric>

#include <gsCore/gsExport.h>
#include <gsCore/gsDebug.h>

#ifdef __GNUC__ 
#include <cxxabi.h>
#include <cstdlib>
#endif

#define STRINGIGY(x) #x

namespace gismo
{

/** @namespace util
    
    @brief This namespace gathers several utility functions for
    miscellaneous tasks
*/
namespace util
{

/// \brief Converts value to string, assuming "operator<<" defined on C
/// \ingroup Utils
template<typename C>
std::string to_string(const C & value)
{
    std::ostringstream convert;
    convert << value;
    return convert.str();
}

#if __cplusplus > 199711L
using std::to_string;
using std::iota;
using std::stod;
using std::stoi;

#else

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

/// \brief Replaces appearance of \a oldStr with \a newStr inside the
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

/// \brief Returns the \a i-th token of the string \a str using delimiter \a delim
/// \ingroup Utils
inline std::string tokenize(const std::string& str,
                            const std::string& delim,
                            const std::size_t token)
{
    std::size_t token_begin = 0;
    std::size_t token_end   = str.find_first_of(delim);
    
    for (std::size_t i=0; i<token; i++) {
        
        GISMO_ENSURE(token_end < std::string::npos,
                     "Requested token exceeds the number of tokens");
        
        token_begin += ++token_end;
        token_end    = str.substr(token_begin).find_first_of(delim);
    }
    
    return str.substr(token_begin,token_end);
}

/// \brief Capitalize string in situ
/// \ingroup Utils
inline void capitalize(std::string& str)
{
    str[0] = static_cast<char>(toupper(str[0]));
}

/// \brief Capitalize string
/// \ingroup Utils
inline std::string capitalize(const std::string& str)
{
    std::string newStr = str;
    capitalize(newStr);
    return newStr;
}

/// \brief Remove pointer from type
/// \ingroup Utils
template<typename T> struct remove_pointer {typedef T type;};
template<typename T> struct remove_pointer<T*> {typedef typename remove_pointer<T>::type type;};

/// \brief Print name of template type as a string
/// \ingroup Utils
template<typename T>
struct type
{
public:
    static std::string name()
    {
#ifdef __GNUC__ 
        int status = 0;
        // Note: C++11 style:
        //std::unique_ptr<char,decltype(std::free)*> dm(__cxxabiv1::__cxa_demangle( typeid(T).name(), NULL, NULL, &status ), std::free);
        char * dm = __cxxabiv1::__cxa_demangle( typeid(T).name(), NULL, NULL, &status );
        GISMO_ASSERT(0==status, "Demangling failed");
        std::string res(dm);
        free(dm);
        return res;
#else
        return typeid(T).name();
#endif
    }
};


} // end namespace util

} // end namespace gismo

