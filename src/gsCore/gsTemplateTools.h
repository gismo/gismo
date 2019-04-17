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
#include <complex>
#include <stdint.h>

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
using std::false_type;
using std::integral_constant;
using std::is_base_of;
using std::is_integral;
using std::is_same;
using std::reference_wrapper;
using std::remove_const;
using std::remove_cv;
using std::remove_volatile;
using std::true_type;
using std::make_unsigned;

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

template <class T>
class reference_wrapper
{
public:
    typedef T type;
    
    reference_wrapper(const T& ref) : _ptr(&const_cast<T&>(ref)) { }
    reference_wrapper(const reference_wrapper&o) : _ptr(o._ptr)  { }
    
    reference_wrapper& operator=(const reference_wrapper& o)
    { _ptr = o._ptr; return *this; }
    
    operator T& () const { return *_ptr; }
    T& get() const { return *_ptr; }

private:
    T* _ptr;
};

template<class T, T v>
    struct integral_constant {
    static const T value = v;
    typedef T value_type;
    typedef integral_constant type;
    value_type operator()() const { return value; }
};

typedef integral_constant<bool, true>  true_type;
typedef integral_constant<bool, false> false_type;

template<typename> struct is_integral_base             : false_type {};

template<> struct is_integral_base<bool>               : true_type {};
template<> struct is_integral_base<char>               : true_type {};
template<> struct is_integral_base<signed char>        : true_type {};
template<> struct is_integral_base<unsigned char>      : true_type {};
template<> struct is_integral_base<wchar_t>            : true_type {};
template<> struct is_integral_base<short>              : true_type {};
template<> struct is_integral_base<int>                : true_type {};
template<> struct is_integral_base<long>               : true_type {};
template<> struct is_integral_base<long long>          : true_type {};
template<> struct is_integral_base<unsigned short>     : true_type {};
template<> struct is_integral_base<unsigned int>       : true_type {};
template<> struct is_integral_base<unsigned long>      : true_type {};
template<> struct is_integral_base<unsigned long long> : true_type {};

template< class T > struct remove_const                  { typedef T type; };
template< class T > struct remove_const<const T>         { typedef T type; };

template< class T > struct remove_volatile               { typedef T type; };
template< class T > struct remove_volatile<volatile T>   { typedef T type; };

template< class T >
struct remove_cv { typedef typename remove_volatile<typename remove_const<T>::type>::type type; };

template<typename T> struct is_integral: is_integral_base<typename remove_cv<T>::type> {};

template<class T>
struct make_unsigned;

#define GISMO_MAKE_UNSIGNED(signed_type) \
template<>                               \
struct make_unsigned<signed_type> {      \
    typedef u##signed_type type;         \
};                                       \
template<>                               \
struct make_unsigned<u##signed_type> {   \
    typedef u##signed_type type;         \
};

GISMO_MAKE_UNSIGNED(int8_t)
GISMO_MAKE_UNSIGNED(int16_t)
GISMO_MAKE_UNSIGNED(int32_t)
GISMO_MAKE_UNSIGNED(int64_t)
#undef GISMO_MAKE_UNSIGNED

#endif

/// \brief Remove pointer from type
/// \ingroup Utils
template<typename T> struct remove_pointer {typedef T type;};
template<typename T> struct remove_pointer<T*> {typedef typename remove_pointer<T>::type type;};

/// \brief Type trait is_complex<T> checks if type T is of type
/// std::complex<...>
template <class T> struct is_complex : public false_type {};
template <class T> struct is_complex<const T > : public is_complex<T>{};
template <class T> struct is_complex<volatile const T > : public is_complex<T>{};
template <class T> struct is_complex<volatile T > : public is_complex<T>{};
template <class T> struct is_complex<std::complex<T> > : public true_type{};

/*
template<typename T>
struct is_complex : integral_constant<bool,
                    is_same<T,std::complex<short int    > >::value ||
                    is_same<T,std::complex<int          > >::value ||
                    is_same<T,std::complex<long int     > >::value ||
                    is_same<T,std::complex<long long int> >::value ||
                    is_same<T,std::complex<float        > >::value ||
                    is_same<T,std::complex<double       > >::value ||
                    is_same<T,std::complex<long double  > >::value ||
#ifdef GISMO_WITH_ADIFF
                    is_same<T,std::complex<gismo::ad::DScalar1<real_t, -1> > >::value ||
                    is_same<T,std::complex<gismo::ad::DScalar1<real_t,  2> > >::value ||
                    is_same<T,std::complex<gismo::ad::DScalar1<real_t,  3> > >::value ||
                    is_same<T,std::complex<gismo::ad::DScalar2<real_t, -1> > >::value ||
                    is_same<T,std::complex<gismo::ad::DScalar2<real_t,  2> > >::value ||
                    is_same<T,std::complex<gismo::ad::DScalar2<real_t,  3> > >::value ||
#endif
#ifdef GISMO_WITH_CODIPACK
                    is_same<T,std::complex<codi::RealForward> >::value ||
                    is_same<T,std::complex<codi::RealReverse> >::value ||
#endif
#ifdef GISMO_WITH_MPFR
                    is_same<T,std::complex<mpfr::mpreal> >::value      ||
#endif
#ifdef GISMO_WITH_MPQ
                    is_same<T,std::complex<mpq_class> >::value         ||
#endif
#ifdef GISMO_WITH_UNUM
                    is_same<T,std::complex<posit_32_2> >::value        ||
#endif
                    is_same<T,std::complex<real_t> >::value
                    > {};
*/
} // end namespace util

} // end namespace gismo
