/**
   Core automatic differentiation functionality for a specific AUTODIFF_TYPE

   This file is included multiple times from gsAutoDiff2.h with different
   AUTODIFF_TYPE definitions to provide support for dual, real, and var types.

   Copyright (c) 2012 by Wenzel Jakob. Based on code by Jon Kaldor
   and Eitan Grinspun.

   Modifications for G+Smo, Angelos Mantzaflaris, 2015
   Additional modifications for multi-mode support, H.M. Verhelst, 2025

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
*/

// This file should only be included from gsAutoDiff2.h
#ifndef AUTODIFF_TYPE
#error "AUTODIFF_TYPE must be defined before including this file"
#endif

#ifndef GSAUTODIFF2_DETAIL_FUNCTIONS_DEFINED
#define GSAUTODIFF2_DETAIL_FUNCTIONS_DEFINED
// Define functions in autodiff::detail namespace first for ADL
// This allows Eigen and other code to find them via ADL when working with Dual types
namespace autodiff {
namespace detail {

// Rounding functions for Dual types
template<typename T, typename G>
inline Dual<T, G> ceil(const Dual<T, G>& x)
{
   return Dual<T, G>(std::ceil(x.val));
}

template<typename T, typename G>
inline Dual<T, G> floor(const Dual<T, G>& x)
{
   return Dual<T, G>(std::floor(x.val));
}

// Rounding functions for expression templates
template<typename Expr>
inline typename std::enable_if<
    isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type ceil(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return DualT(std::ceil(evaluated.val));
}

template<typename Expr>
inline typename std::enable_if<
    isExpr<Expr>,
    typename autodiff::detail::DualType<Expr>
>::type floor(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return DualT(std::floor(evaluated.val));
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
    isExpr<Expr>,
    bool
>::type isinf(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return std::isinf(evaluated.val);
}

template<typename Expr>
inline typename std::enable_if<
    isExpr<Expr>,
    bool
>::type isnan(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return std::isnan(evaluated.val);
}

template<typename Expr>
inline typename std::enable_if<
    isExpr<Expr>,
    bool
>::type isfinite(const Expr& expr)
{
   typedef typename autodiff::detail::DualType<Expr> DualT;
   DualT evaluated(expr);
   return std::isfinite(evaluated.val);
}

} // namespace detail
} // namespace autodiff
#endif // GSAUTODIFF2_DETAIL_FUNCTIONS_DEFINED

namespace gismo {
namespace math {

#ifndef GSAUTODIFF2_DUAL_MATH_FUNCTIONS_DEFINED
#define GSAUTODIFF2_DUAL_MATH_FUNCTIONS_DEFINED
// Math function extensions for Dual types
template<typename T, typename G>
inline autodiff::detail::Dual<T, G> log10(const autodiff::detail::Dual<T, G>& x)
{
   using autodiff::detail::log10;
   return log10(x);
}

// Bring ceil and floor from autodiff::detail into gismo::math for uniform access
using autodiff::detail::ceil;
using autodiff::detail::floor;

template<typename T, typename G>
inline autodiff::detail::Dual<T, G> round(const autodiff::detail::Dual<T, G>& x)
{
   return autodiff::detail::Dual<T, G>(std::round(x.val));
}

template<typename T, typename G>
inline autodiff::detail::Dual<T, G> trunc(const autodiff::detail::Dual<T, G>& x)
{
   return autodiff::detail::Dual<T, G>(std::trunc(x.val));
}

// Expression template versions
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

// ceil and floor for expression templates are now in autodiff::detail namespace (see below)

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
inline autodiff::detail::Dual<T,G> max(const autodiff::detail::Dual<T,G>& a, const autodiff::detail::Dual<T,G>& b)
{
    return a.val > b.val ? a : b;
}

template<typename T, typename G>
inline autodiff::detail::Dual<T,G> min(const autodiff::detail::Dual<T,G>& a, const autodiff::detail::Dual<T,G>& b)
{
    return a.val < b.val ? a : b;
}

// Max with expression template on left
template<typename Expr, typename T, typename G>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    autodiff::detail::Dual<T,G>
>::type max(const Expr& a, const autodiff::detail::Dual<T,G>& b)
{
    autodiff::detail::Dual<T,G> a_eval(a);
    return a_eval.val > b.val ? a_eval : b;
}

// Max with expression template on right
template<typename T, typename G, typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    autodiff::detail::Dual<T,G>
>::type max(const autodiff::detail::Dual<T,G>& a, const Expr& b)
{
    autodiff::detail::Dual<T,G> b_eval(b);
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
    autodiff::detail::Dual<T,G>
>::type min(const Expr& a, const autodiff::detail::Dual<T,G>& b)
{
    autodiff::detail::Dual<T,G> a_eval(a);
    return a_eval.val < b.val ? a_eval : b;
}

// Min with expression template on right
template<typename T, typename G, typename Expr>
inline typename std::enable_if<
    autodiff::detail::isExpr<Expr>,
    autodiff::detail::Dual<T,G>
>::type min(const autodiff::detail::Dual<T,G>& a, const Expr& b)
{
    autodiff::detail::Dual<T,G> b_eval(b);
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

// Bring isinf, isnan, isfinite from autodiff::detail into gismo::math for uniform access
using autodiff::detail::isinf;
using autodiff::detail::isnan;
using autodiff::detail::isfinite;
#endif // GSAUTODIFF2_DUAL_MATH_FUNCTIONS_DEFINED

} // namespace math
} // namespace gismo

#ifndef GSAUTODIFF2_COMMON_OPERATORS_DEFINED
#define GSAUTODIFF2_COMMON_OPERATORS_DEFINED
// Add stream output operators for autodiff expression types
namespace gismo {
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

namespace std {
   // Add std::to_string specialization for autodiff dual numbers
   template<typename T, typename G>
   string to_string(const autodiff::detail::Dual<T,G>& dual)
   {
      std::ostringstream oss;
      oss << dual.val;
      return oss.str();
   }
} // namespace std

// Provide operator>> for Dual in autodiff::detail (ADL)
namespace autodiff { namespace detail {
template<typename T, typename G>
inline std::istream& operator>>(std::istream& is, Dual<T,G>& dual)
{
    T v;
    if (is >> v)
        dual = Dual<T,G>(v); // gradient zero by default
    return is;
}
}} // namespace autodiff::detail
#endif // GSAUTODIFF2_COMMON_OPERATORS_DEFINED

// ============================================================================
// Support for reverse-mode autodiff::var type
// ============================================================================
#ifdef AUTODIFF_VAR_TYPE

namespace gismo {
namespace math {

// Math functions for var type
inline autodiff::var log10(const autodiff::var& x)
{
   return log(x) / log(10.0);
}

inline autodiff::var abs(const autodiff::var& x)
{
   // Use autodiff's built-in abs function which handles the derivatives correctly
   return autodiff::abs(x);
}

inline autodiff::var floor(const autodiff::var& x)
{
   return autodiff::var(std::floor(val(x)));
}

inline autodiff::var ceil(const autodiff::var& x)
{
   return autodiff::var(std::ceil(val(x)));
}

inline autodiff::var round(const autodiff::var& x)
{
   return autodiff::var(std::round(val(x)));
}

inline autodiff::var trunc(const autodiff::var& x)
{
   return autodiff::var(std::trunc(val(x)));
}

inline autodiff::var max(const autodiff::var& a, const autodiff::var& b)
{
    return val(a) > val(b) ? a : b;
}

inline autodiff::var min(const autodiff::var& a, const autodiff::var& b)
{
    return val(a) < val(b) ? a : b;
}

inline bool isinf(const autodiff::var& x)
{
   return std::isinf(val(x));
}

inline bool isnan(const autodiff::var& x)
{
   return std::isnan(val(x));
}

inline bool isfinite(const autodiff::var& x)
{
   return std::isfinite(val(x));
}

} // namespace math
} // namespace gismo

// Add isnan, isinf, isfinite to autodiff namespace for ADL
namespace autodiff {
using gismo::math::isinf;
using gismo::math::isnan;
using gismo::math::isfinite;
} // namespace autodiff

namespace gismo {
// Stream output for var
inline std::ostream& operator<<(std::ostream& os, const autodiff::var& v)
{
    return os << val(v);
}
} // namespace gismo

namespace std {
// Add std::to_string specialization for var
inline string to_string(const autodiff::var& v)
{
   std::ostringstream oss;
   oss << val(v);
   return oss.str();
}
} // namespace std

namespace autodiff {
// Provide operator>> for var
inline std::istream& operator>>(std::istream& is, var& v)
{
    double val_v;
    if (is >> val_v)
        v = var(val_v);
    return is;
}
} // namespace autodiff

#endif // AUTODIFF_VAR_TYPE
