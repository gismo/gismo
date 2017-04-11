/** @file gsTemplateTools.h

    @brief Utilities related to template programming.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsExport.h>
#include <utility>

namespace gismo
{

namespace util {

#if __cplusplus >= 201103 || (defined(_MSC_VER) && _MSC_VER > 1700)
//see also http://lists.boost.org/Archives/boost/2009/04/151209.php
template <typename T> struct has_move_constructor
{
    typedef char yes[1];
    typedef char no[2];
    
    struct AmbiguousConverter
    {
        operator T&& ();
        operator const T& ();
    };
    template <typename C> static no& test(decltype( new C( AmbiguousConverter{} )));
    template <typename> static yes& test(...);
    enum { value = (sizeof(test<T>(0)) == sizeof(yes)) };
};

using std::conditional;
using std::enable_if;
using std::is_same;
using std::is_base_of;

#else

// template <typename T> struct has_move_constructor { enum { value = 0 }; };

template<bool B, class T, class F> struct conditional { typedef T type; };
template<class T, class F> struct conditional<false, T, F> { typedef F type; };

template<bool B, class T = void> struct enable_if {};
template<class T> struct enable_if<true, T> { typedef T type;};

template<class T, class U> struct is_same { enum { value = 0 }; };
template<class T>          struct is_same<T, T> { enum { value = 1 }; };

template <typename B, typename D> struct Host
{ operator B*() const; operator D*(); };
template <typename B, typename D>
struct is_base_of
{
    typedef char (&yes)[1];
    typedef char (&no)[2];
    template <typename T> static yes check(D*, T);
    static no check(B*, int);
    static const bool value = sizeof(check(Host<B,D>(), int())) == sizeof(yes);
};
#endif

} // end namespace util

} // end namespace gismo
