/** @file gsMath.h

    @brief Mathematical functions for use in G+Smo.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan, A. Mantzaflaris
*/

#pragma once

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <limits>

namespace gismo {

/** @namespace gismo::math

    @brief
    This namespace contains common mathematical functions.
    
    \ingroup Core
*/
namespace math {

using std::numeric_limits;

// Math functions
using std::abs;
using std::pow;
using std::fabs;
using std::sqrt;
using std::ceil;
using std::floor;
using std::cos;
using std::cosh;
using std::sin;
using std::sinh;
using std::tan;
using std::tanh;
using std::acos;
using std::log;
using std::log10;

template <typename T>
T min(T a, T b) {return  (a < b ? a : b); }

template <typename T>
T max(T a, T b) {return  (a < b ? b : a); }

// functions to check for floating point errors
// Note: Some duplication in gsUtils/gsDebug.h
// Get isnan/isinf working on different compilers
#ifdef _MSC_VER
#include <float.h>
template <typename T>
int isnan (T a) {return _isnan(a); }
template <typename T>
int isfinite(T a) {return _finite(a);}
template <typename T>
int isinf(T a) {return (_FPCLASS_PINF|_FPCLASS_NINF) & _fpclass(a);}

#ifndef NAN
// MSVC doesn't have the NAN constant in cmath, so we use the C++ standard definition
#define NAN (std::numeric_limits<real_t>::quiet_NaN())
#endif

#else
#ifdef _INTEL_COMPILER
#include <mathimf.h>
#else
using std::isnan;
using std::isfinite;
using std::isnormal;
using std::isinf;
#endif
#endif

#ifdef __MPREAL_H__
// Math functions for MPFR
using mpfr::abs;
using mpfr::pow;
using mpfr::fabs;
using mpfr::sqrt;
using mpfr::ceil;
using mpfr::floor;
using mpfr::cos;
using mpfr::sin;
using mpfr::acos;
using mpfr::log;
using std::log10;
#endif



/// Returns the sign of \a val
template <typename T> int getSign(T val)
{
    return (T(0) < val) - (val < T(0));
}

/// integer power
inline int ipow(int x, unsigned p)
{
    if (p == 0) return 1;
    if (p == 1) return x;

    int tmp = ipow(x, p/2);
    if (p % 2 == 0)
        return tmp * tmp;
    else
        return x * tmp * tmp;
}

/// Returns convex combination of \a a and \a b with weight \a t
template<class T>
T mix(T const & a, T const & b, T const & t)
{
    return (1 - t) * a + t * b;
}

// numerical comparison
template<class T>
bool lessthan (T const& a, T const& b)
{ return ( b-a > 1e-6); }

// numerical equality
template<class T>
bool almostEqual (T const& a, T const& b)
{  return ( math::abs(b-a) < 1e-6); }

// numerical equality with \a prec decimal digits
template<int prec, class T>
bool almostEqual (T const& a, T const& b)
{  return ( math::log10(math::abs(b-a)) < -prec); }

// static const double pi      =  3.141592653589793238462;
// static const double e       =  2.718281828459045235360;
// static const double pi_2    =  1.570796326794896619231;
// static const double pi_4    =  0.785398163397448309616;
// static const double pi_180  =  0.017453292519943295769;
// static const double _1_pi   =  0.318309886183790671538;
// static const double _2_pi   =  0.636619772367581343076;
// static const double _180_pi = 57.295779513082320876798;

}

/**
    \brief tests if the difference between two numbers is below tolerance
    **/
template <typename T>
inline bool gsClose (const T &a, const T &b, const T &tol)
{
    return math::abs(a-b) <= tol;
}

/**
    \brief tests if the difference between two matrices is bounded by tol in L^\infty norm

    The tolerance is relative to maximum absolute values of the entries of the matrices.
    **/
template <typename matrix_t1, typename matrix_t2>
inline bool gsAllCloseRelativeToMax (const matrix_t1 &a, const matrix_t2 &b, const typename matrix_t1::Scalar &tol )
{
    GISMO_ASSERT( a.cols()==b.cols() && a.rows()==b.rows(), "Only matrices of the same size can be compared" );
    return (a-b).array().abs().maxCoeff() <= tol*math::max(a.array().abs().maxCoeff(), b.array().abs().maxCoeff());
}

/**
    \brief tests if the difference between two matrices is bounded by tol in L^\infty norm

    The tolerance is absolute, therefore independent of the matrix entries.
    **/
template <typename matrix_t1, typename matrix_t2>
inline bool gsAllCloseAbsolute (const matrix_t1 &a, const matrix_t2 &b, const typename matrix_t1::Scalar &tol )
{
    GISMO_ASSERT( a.cols()==b.cols() && a.rows()==b.rows(), "Only matrices of the same size can be compared" );
    return (a-b).array().abs().maxCoeff() <= tol;
}

/**
    \brief tests if the difference between two matrices is bounded by tol in L^\infty norm

    The tolerance is absolute below the reference, but relative for bigger numbers.
    **/
template <typename matrix_t1, typename matrix_t2>
inline bool gsAllCloseRelAndAbsWithRef (const matrix_t1 &a, const matrix_t2 &b, const typename matrix_t1::Scalar &tol, const typename matrix_t1::Scalar &ref )
{
    GISMO_ASSERT( a.cols()==b.cols() && a.rows()==b.rows(), "Only matrices of the same size can be compared" );
    return (a-b).array().abs().maxCoeff() <= tol*math::max(ref,a.array().abs().maxCoeff(), b.array().abs().maxCoeff());
}


/*
  CompileTimeLog2 computes the logarithm base 2 of the argument
  using template recursion.  
*/
template <unsigned arg>
struct CompileTimeLog2
{
    enum { result = CompileTimeLog2< arg / 2 >::result + 1 };
};
template <>
struct CompileTimeLog2<1>
{
    enum { result = 0 };
};
template <>
struct CompileTimeLog2<0>
{
    // empty because the logarithm of 0 is -infinity
};
template <unsigned arg>
unsigned CTLog2 (void) {return CompileTimeLog2<arg>::result;}

/*
//Template to get the minimum of two numbers at compile time
template <unsigned a, unsigned b>
struct min
{
    enum { result = a<b ? a:b};
};
*/

}
