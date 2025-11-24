/**
   Automatic differentiation data type for C++, depends on the Eigen
   linear algebra library.

   Copyright (c) 2012 by Wenzel Jakob. Based on code by Jon Kaldor
   and Eitan Grinspun.

   Modifications for G+Smo, Angelos Mantzaflaris, 2015

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
*/

#pragma once

#define Eigen gsEigen
#include <autodiff/forward/real.hpp>
// Note: Do not include eigen.hpp files here as they pull in Eigen/Core too early
// They will be included from gsAutoDiffEigen.h after Eigen plugins are set up
// #include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual.hpp>
// #include <autodiff/forward/dual/eigen.hpp>
#include <cmath>
#include <sstream>

// Short-hand type alias for autodiff dual numbers
using autodiff::detail::Dual;

namespace autodiff {
   // Helper to extract the base scalar type from autodiff types
   template<class T>
   struct autodiff_real_helper {
      typedef T type;
   };

   // Specialization for Dual numbers - recursively extract the value type
   template<typename T, typename G>
   struct autodiff_real_helper<autodiff::detail::Dual<T, G>> {
      typedef typename autodiff_real_helper<T>::type type;
   };

   template<class T>
   using autodiff_real = typename autodiff_real_helper<T>::type;
}

namespace gismo {
namespace math {

// Bring in std math functions by default
using std::floor;
using std::ceil;
using std::round;
using std::trunc;
using std::abs;
using std::log10;
using std::isinf;
using std::isnan;
using std::isfinite;

// Extend autodiff with missing math functions
// Add log10 for Dual types
template<typename T, typename G>
inline Dual<T, G> log10(const Dual<T, G>& x)
{
   using autodiff::detail::log10;
   return log10(x);
}

template<typename T, typename G>
inline Dual<T, G> floor(const Dual<T, G>& x)
{
   return Dual<T, G>(std::floor(x.val));
}

template<typename T, typename G>
inline Dual<T, G> ceil(const Dual<T, G>& x)
{
   return Dual<T, G>(std::ceil(x.val));
}

template<typename T, typename G>
inline Dual<T, G> round(const Dual<T, G>& x)
{
   return Dual<T, G>(std::round(x.val));
}

template<typename T, typename G>
inline Dual<T, G> trunc(const Dual<T, G>& x)
{
   return Dual<T, G>(std::trunc(x.val));
}


// Add log10 for expression templates
template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type log10(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   using autodiff::detail::log10;
   return log10(evaluated);
}

// For expression templates, evaluate to Dual first then apply
template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type floor(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return DualT(std::floor(evaluated.val));
}

template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type ceil(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return DualT(std::ceil(evaluated.val));
}

template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type round(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return DualT(std::round(evaluated.val));
}

template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type trunc(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return DualT(std::trunc(evaluated.val));
}

// Max/min for Dual types
template<typename T, typename G>
inline Dual<T,G> max(const Dual<T,G>& a, const Dual<T,G>& b)
{
    return a.val > b.val ? a : b;
}

template<typename T, typename G>
inline Dual<T,G> min(const Dual<T,G>& a, const Dual<T,G>& b)
{
    return a.val < b.val ? a : b;
}

// Max with expression template on left
template<typename Expr, typename T, typename G>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    Dual<T,G>
>::type max(const Expr& a, const Dual<T,G>& b)
{
    Dual<T,G> a_eval(a);
    return a_eval.val > b.val ? a_eval : b;
}

// Max with expression template on right
template<typename T, typename G, typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    Dual<T,G>
>::type max(const Dual<T,G>& a, const Expr& b)
{
    Dual<T,G> b_eval(b);
    return a.val > b_eval.val ? a : b_eval;
}

// Max with both expression templates
template<typename Expr1, typename Expr2>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr1> && autodiff::detail::isExpr<Expr2>,
    typename autodiff::detail::DualType<Expr1>
>::type max(const Expr1& a, const Expr2& b)
{
    typedef typename autodiff::detail::DualType<Expr1> DualT;
    DualT a_eval(a);
    DualT b_eval(b);
    return a_eval.val > b_eval.val ? a_eval : b_eval;
}

// Min with expression template on left
template<typename Expr, typename T, typename G>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    Dual<T,G>
>::type min(const Expr& a, const Dual<T,G>& b)
{
    Dual<T,G> a_eval(a);
    return a_eval.val < b.val ? a_eval : b;
}

// Min with expression template on right
template<typename T, typename G, typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    Dual<T,G>
>::type min(const Dual<T,G>& a, const Expr& b)
{
    Dual<T,G> b_eval(b);
    return a.val < b_eval.val ? a : b_eval;
}

// Min with both expression templates
template<typename Expr1, typename Expr2>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr1> && autodiff::detail::isExpr<Expr2>,
    typename autodiff::detail::DualType<Expr1>
>::type min(const Expr1& a, const Expr2& b)
{
    typedef typename autodiff::detail::DualType<Expr1> DualT;
    DualT a_eval(a);
    DualT b_eval(b);
    return a_eval.val < b_eval.val ? a_eval : b_eval;
}

// Boolean functions for Dual types
template<typename T, typename G>
inline bool isinf(const Dual<T,G>& x)
{
   return std::isinf(x.val);
}

template<typename T, typename G>
inline bool isnan(const Dual<T,G>& x)
{
   return std::isnan(x.val);
}

template<typename T, typename G>
inline bool isfinite(const Dual<T,G>& x)
{
   return std::isfinite(x.val);
}

// Boolean functions for expression templates
template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    bool
>::type isinf(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return std::isinf(evaluated.val);
}

template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    bool
>::type isnan(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return std::isnan(evaluated.val);
}

template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    bool
>::type isfinite(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return std::isfinite(evaluated.val);
}

} // namespace math

// Also add the math functions to global namespace for Eigen compatibility
using gismo::math::ceil;
using gismo::math::floor;
using gismo::math::round;
using gismo::math::trunc;
using gismo::math::isinf;
using gismo::math::isnan;
using gismo::math::isfinite;

// Add stream output operators for autodiff expression types
template<typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    std::ostream&
>::type operator<<(std::ostream& os, const Expr& expr)
{
    typedef typename autodiff::detail::DualType<Expr> DualT;
    DualT evaluated(expr);
    return os << evaluated.val;
}

} // namespace gismo

namespace std
{
   // Add std::to_string specialization for autodiff dual numbers
   template<typename T, typename G>
   string to_string(const Dual<T,G>& dual)
   {
      std::ostringstream oss;
      oss << dual.val;
      return oss.str();
   }

   // Make the other functions available in std namespace for compatibility
   using ::gismo::math::ceil;
   using ::gismo::math::floor;
   using ::gismo::math::round;
   using ::gismo::math::trunc;
   using ::gismo::math::isinf;
   using ::gismo::math::isnan;
   using ::gismo::math::isfinite;
}

// // Stream input operator for Dual numbers
// template<typename T, typename G>
// inline std::istream& operator>>(std::istream& is, autodiff::detail::Dual<T,G>& dual)
// {
//    T val;
//    is >> val;
//    dual = autodiff::detail::Dual<T,G>(val); // Create dual with zero gradient
//    return is;
// }

// // Stream output operator for Dual numbers (if not already present)
// template<typename T, typename G>
// inline std::ostream& operator<<(std::ostream& os, const autodiff::detail::Dual<T,G>& dual)
// {
//    return os << dual.val;
// }

// // Stream input operator for autodiff expression templates
// template<typename Expr>
// inline typename std::enable_if<
//     autodiff::detail::isExpr<Expr>,
//     std::istream&
// >::type operator>>(std::istream& is, Expr& expr)
// {
//    // Cannot directly read into an expression, so this is not supported
//    // Expressions are temporaries and should not be lvalues for >>
//    // If this is hit, the code should use a Dual variable instead
//    static_assert(!autodiff::detail::isExpr<Expr>,
//                  "Cannot read into autodiff expression template. Use a Dual variable.");
//    return is;
// }

// Provide only operator>> for Dual, in autodiff::detail (ADL)
namespace autodiff { namespace detail {

template<typename T, typename G>
inline std::istream& operator>>(std::istream& is, Dual<T,G>& dual)
{
    T v;
    if (is >> v)
        dual = Dual<T,G>(v); // gradient zero by default
    return is;
}
}};
#   undef Eigen

