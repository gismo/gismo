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
#include <gsCore/gsMemory.h>

#ifdef __GNUC__ 
#include <cxxabi.h>
#include <cstdlib>
#endif

#define STRINGIFY(x) #x

namespace gismo
{

/** @namespace util
    
    @brief This namespace gathers several utility functions for
    miscellaneous tasks
*/
namespace util
{

#if __cplusplus >= 201103L || _MSC_VER >= 1600
template <class C, size_t N> // we catch up char arrays
std::string to_string(C (& value)[N])
{
    static_assert(!std::is_same<C[N], char[N]>::value, "Character arrays are not allowed");
    std::ostringstream convert;
    convert << value;
    return convert.str();
}
#endif

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

/// \brief Checks if a string \a haystack ends with the string \a needle
/// \ingroup Utils
inline bool ends_with( const std::string & haystack, const std::string & needle )
{
    if (needle.size() > haystack.size()) return false;
    //std::transform(haystack.begin(), value.end(), tmp.begin(), ::tolower); 
    return std::equal(needle.rbegin(), needle.rend(), haystack.rbegin());
}

#if __cplusplus > 199711L || _MSC_VER >= 1600
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

/// \brief equivalent to std::stoi(str), and therefore std::stoi(str, 0, 10)
inline int stoi(const std::string& str)
{
    std::istringstream ss(str);
    int i;
    if (!(ss >> std::skipws >> i)) // leading whitespaces are ignored by std::stoi
        //Extracting an int failed
        throw std::invalid_argument("stoi");    // if CXX11 code throws, CXX98 should do too, or?

    // std::stoi ignores all after a valid number
    //char c;
    //if (ss >> c)
    //    //There was something after the number
    //    return 0;
    
    return i;
}

/// \brief equivalent to std::stod(str)
inline double stod(const std::string& str)
{
    std::istringstream ss(str);
    double i;
    if (!(ss >> i))
        //Extracting double failed
        throw std::invalid_argument("stod");

    size_t pos;

    // hex is valid for std::stod - not a good implementation yet
    // ssi and ssd correct, better move to double need be done => convert to decimal string and let it parse at usual way.
    if(i == 0 && (((pos = str.find("0x")) != std::string::npos) || ((pos = str.find("0X")) != std::string::npos)))
    {
        bool negative = false;
        if (pos > 0)
            if (str[pos - 1] == '-')
                negative = true;

        size_t comma = str.find(".", pos+2);

        size_t integer, decimal;
        std::istringstream ssi(str.substr(pos+2, comma - pos - 2));
        std::istringstream ssd(str.substr(++comma)); // we need always comma+1

        if (!(ssi >> std::hex >> integer))
            throw std::invalid_argument("stod");

        if (!(ssd >> std::hex >> decimal))
            throw std::invalid_argument("stod");

        size_t lenght = str.find_first_not_of("0123456789abcdefABCDEF", comma);
        if (lenght == std::string::npos)
            lenght = str.length() - comma;
        else
            lenght -= comma;

        i = (integer + (decimal/pow(16, lenght))) * (negative ? -1. : 1.);
    }

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

/// \brief Returns the \a i-th token of the string \a str using any character in \a delim as delimiter without counting
/// empty sequences.
/// \example util::tokenize("abcdbca", "bd", ...) => {a, c, ca}
/// \ingroup Utils
inline std::string tokenize(const std::string& str,
                            const std::string& delim,
                            const size_t token)
{
    size_t token_end = std::string::npos;
    size_t token_begin = 0;
    size_t token_count = 0;
    bool catched = false;

    do
    {
        GISMO_ENSURE(!catched,
                     "Requested token exceeds the number of tokens");

        token_begin = token_end + 1;
        token_end = str.find_first_of(delim, token_begin);

        if(token_end == std::string::npos)  // catch in next iteration
            catched = true;

        if (token_end != token_begin) // ignore empty sequences
            ++token_count;
    }
    while (token_count <= token);

    return str.substr(token_begin, token_end - token_begin);
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
#else   // not __GNUC__
        return typeid(T).name();
#endif  // __GNUC__
    }
};

/// \brief Create hash key for a rangle of (integral) numbers
template<typename T>
size_t hash_range(T const * start, const T * const end)
{
    size_t seed = end - start;
    for(; start!=end; ++start) 
        seed ^= *start + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

#if __cplusplus >= 201703L || _MSVC_LANG >= 201703L
using std::size;
#else
template <class T, size_t N>
size_t size(const T (&)[N])
{
    return N;
}
template <class T>
size_t size(const T& t)
{
    return t.size();
}
#endif

} // end namespace util

// This macro assumes the operators == and < to be present and
// defines other four operators !=, >, <= and >=
#define GISMO_DELEGATING_COMPARISON_OPERATORS( T )                  \
inline bool operator!= (const T& a, const T& b) { return !(a==b); } \
inline bool operator>  (const T& a, const T& b) { return b<a;     } \
inline bool operator<= (const T& a, const T& b) { return !(b<a);  } \
inline bool operator>= (const T& a, const T& b) { return !(a<b);  }

// This macro deletes the operators ==, !=, <, >, <= and >=
// for operations that involve the types S and T (in either
// order)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
#define GISMO_DELETE_COMPARISON_OPERATORS( S, T )         \
inline bool operator== (const S& a, const T& b) = delete; \
inline bool operator!= (const S& a, const T& b) = delete; \
inline bool operator<  (const S& a, const T& b) = delete; \
inline bool operator>  (const S& a, const T& b) = delete; \
inline bool operator<= (const S& a, const T& b) = delete; \
inline bool operator>= (const S& a, const T& b) = delete; \
inline bool operator== (const T& a, const S& b) = delete; \
inline bool operator!= (const T& a, const S& b) = delete; \
inline bool operator<  (const T& a, const S& b) = delete; \
inline bool operator>  (const T& a, const S& b) = delete; \
inline bool operator<= (const T& a, const S& b) = delete; \
inline bool operator>= (const T& a, const S& b) = delete;
#else
#define GISMO_DELETE_COMPARISON_OPERATORS( S, T ) \
inline bool operator== (const S& a, const T& b); \
inline bool operator!= (const S& a, const T& b); \
inline bool operator<  (const S& a, const T& b); \
inline bool operator>  (const S& a, const T& b); \
inline bool operator<= (const S& a, const T& b); \
inline bool operator>= (const S& a, const T& b); \
inline bool operator== (const T& a, const S& b); \
inline bool operator!= (const T& a, const S& b); \
inline bool operator<  (const T& a, const S& b); \
inline bool operator>  (const T& a, const S& b); \
inline bool operator<= (const T& a, const S& b); \
inline bool operator>= (const T& a, const S& b);

#endif

} // end namespace gismo

