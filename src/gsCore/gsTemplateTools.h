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

namespace gismo
{

/*
  Compile time type switching.
  Example: The type

  gismo::conditional<k==1,double, gsVector<double> >::type

  evaluates (at compile time) to "double" if k is equal to 1, 
  and to "gsVector<double>" otherwise
*/

template <bool flag, class IsTrue, class IsFalse>
struct conditional;

template <class IsTrue, class IsFalse>
struct conditional<true, IsTrue, IsFalse> { typedef IsTrue type; };

template <class IsTrue, class IsFalse>
struct conditional<false, IsTrue, IsFalse> {
   typedef IsFalse type;
};

namespace util {

template<bool B, class T, class F> struct conditional { typedef T type; };
template<class T, class F> struct conditional<false, T, F> { typedef F type; };

template<bool B, class T = void> struct enable_if {};
template<class T> struct enable_if<true, T> { typedef T type;};

template<class T, class U> struct is_same { enum { value = 0 }; };
template<class T>          struct is_same<T, T> { enum { value = 1 }; };

template <typename B, typename D> struct Host
{
    operator B*() const;
    operator D*();
};

template <typename B, typename D>
struct is_base_of
{
    typedef char (&yes)[1];
    typedef char (&no)[2];
    template <typename T> 
    static yes check(D*, T);
    static no check(B*, int);
    static const bool value = sizeof(check(Host<B,D>(), int())) == sizeof(yes);
};

} // end namespace internal

} // end namespace gismo
