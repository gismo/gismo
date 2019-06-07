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

/// converts a \a std::string to \a T
template<typename T> inline T convert_to             (const std::string&);
template<> inline         int convert_to<int>        (const std::string& s) { return stoi(s); }
template<> inline      double convert_to<double>     (const std::string& s) { return stod(s); }
template<> inline std::string convert_to<std::string>(const std::string& s) { return s; }

/// \brief Creates a container \a T, splits and adds the std::string \a str
/// into parts of type \a T::value_type with char \a ch as separator.
/// T can be any container from std, that implements insert(iterator pos, const T::value_type& val);
/// Conversions from std::String to T::value_type must be defined as
/// TYPE convert_to<TYPE>(const std::string& s){return ...}
/// \example
/// std::set<index_t> s split<std::set<index_t> >("1,3,2,4,3,2,1", ','); s => {1,2,3,4}
/// std::vector<std::string> v split<std::vector<std::string> >("ab.ba.aa", '.'); v => {"ab","ba","aa"}
/// \tparam T STL container type, like std::vector<std::string> >, but not std::map
/// \param str input string
/// \param ch separator
/// \param empty if true, adds also empty parts. default is false. Makes only sense with std::string.
/// \return a container of T<T::value_type>
template <class T>
inline T split(const std::string& str, char ch, bool empty = false)
{
    T result;
    size_t start = 0;
    size_t found;
    while ((found = str.find(ch, start)) != std::string::npos)
    {
        if (!empty && (found == start)) // we ignore if some separators lined up together
        {
            ++start;
            continue;
        }
        result.insert(result.end(), convert_to<typename T::value_type>(str.substr(start, found - start)));
        start = found + 1;
    }
    if(start < str.length())
        result.insert(result.end(), convert_to<typename T::value_type>(str.substr(start, str.length())));
    return result;
}

/// \brief Returns the \a i-th token of the string \a str using delimiter \a delim
/// \ingroup Utils
inline std::string tokenize(const std::string& str,
                            const std::string& delim,
                            const size_t token)
{
    size_t token_begin = 0;
    size_t token_end   = str.find_first_of(delim);
    
    for (size_t i=0; i<token; i++) {
        
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
size_t hash_range(T const * start, const T * const end)
{
    size_t seed = end - start;
    for(; start!=end; ++start) 
        seed ^= *start + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

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
#if __cplusplus >= 201103L
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

