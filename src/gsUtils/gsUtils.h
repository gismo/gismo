/** @file gsUtils.h

    @brief Several utility functions for miscellaneous tasks

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Author(s): A. Mantzaflaris, Harald Weiner, J. Vogl
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

/// \brief Checks if a string \a haystack begins with the string \a needle
/// \ingroup Utils
inline bool starts_with( const std::string & haystack, const std::string & needle )
{
    std::string::const_iterator it1 = haystack.begin();
    std::string::const_iterator it2 = needle.begin();
    while ( it2!=needle.end() )
    {
        if ( it1 == haystack.end() || *it1 != *it2) return false;
        it1++; it2++;
    }
    return true;
}

#if __cplusplus > 199711L || (defined(_MSC_VER) && _MSC_VER >= 1600)
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
inline std::string returnCapitalized(const std::string& str)
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
#ifdef  __GNUC__ 
        int status = 0;
#if     __cplusplus > 199711L
        memory::unique_ptr<char,decltype(std::free)*>
            dm(__cxxabiv1::__cxa_demangle( typeid(T).name(), NULL, NULL, &status ), std::free);
        return (status==0) ? dm.get() : typeid(T).name();
#else
        char * dm = __cxxabiv1::__cxa_demangle( typeid(T).name(), NULL, NULL, &status );
        if (status!=0)
        {
            std::free(dm);
            return typeid(T).name();
        }
        std::string res(dm);
        std::free(dm);
        return res;
#endif
#else
        return typeid(T).name();
#endif
    }
};

/// \brief Create hash key for a rangle of (integral) numbers
template<typename T>
std::size_t hash_range(T const * start, const T * const end)
{
    std::size_t seed = end - start;
    for(; start!=end; ++start) 
        seed ^= *start + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

/// \brief Systemspecific path-separator symbol
/// \ingroup Utils
const char SEPARATOR =
#if defined _WIN32 || defined __CYGWIN__
    '\\';
#else
    '/';
#endif

} // end namespace util

} // end namespace gismo

