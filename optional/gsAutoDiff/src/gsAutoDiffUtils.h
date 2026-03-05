/**
   Automatic differentiation data type for C++, depends on the Eigen
   linear algebra library.

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

#pragma once

#include "gsAutoDiffTraits.h"
#include <type_traits>
#include <utility>
#include <cmath>
#include <gsCore/gsForwardDeclarations.h>


namespace gismo {
namespace math {
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> abs(const autodiff::reverse::detail::Variable<T>& x) { using std::abs; return abs(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sqrt(const autodiff::reverse::detail::Variable<T>& x) { using std::sqrt; return sqrt(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::Variable<T>& x, const autodiff::reverse::detail::Variable<T>& y) { using std::pow; return pow(x, y); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::Variable<T>& x, T y) { using std::pow; return pow(x, y); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::Variable<T>& x, int y) { using std::pow; return pow(x, static_cast<T>(y)); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> exp(const autodiff::reverse::detail::Variable<T>& x) { using std::exp; return exp(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> log(const autodiff::reverse::detail::Variable<T>& x) { using std::log; return log(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sin(const autodiff::reverse::detail::Variable<T>& x) { using std::sin; return sin(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> cos(const autodiff::reverse::detail::Variable<T>& x) { using std::cos; return cos(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> tan(const autodiff::reverse::detail::Variable<T>& x) { using std::tan; return tan(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> asin(const autodiff::reverse::detail::Variable<T>& x) { using std::asin; return asin(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> acos(const autodiff::reverse::detail::Variable<T>& x) { using std::acos; return acos(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> atan(const autodiff::reverse::detail::Variable<T>& x) { using std::atan; return atan(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sinh(const autodiff::reverse::detail::Variable<T>& x) { using std::sinh; return sinh(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> cosh(const autodiff::reverse::detail::Variable<T>& x) { using std::cosh; return cosh(x); }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> tanh(const autodiff::reverse::detail::Variable<T>& x) { using std::tanh; return tanh(x); }
    
    // Math functions for ExprPtr<T> (needed when autodiff operations return ExprPtr directly)
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> abs(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::abs;
        return autodiff::reverse::detail::Variable<T>(abs(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sqrt(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::sqrt;
        return autodiff::reverse::detail::Variable<T>(sqrt(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::ExprPtr<T>& x, T y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x, y));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::ExprPtr<T>& x, int y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x, static_cast<T>(y)));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> exp(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::exp;
        return autodiff::reverse::detail::Variable<T>(exp(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> log(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::log;
        return autodiff::reverse::detail::Variable<T>(log(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sin(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::sin;
        return autodiff::reverse::detail::Variable<T>(sin(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> cos(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::cos;
        return autodiff::reverse::detail::Variable<T>(cos(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> tan(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::tan;
        return autodiff::reverse::detail::Variable<T>(tan(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> asin(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::asin;
        return autodiff::reverse::detail::Variable<T>(asin(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> acos(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::acos;
        return autodiff::reverse::detail::Variable<T>(acos(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> atan(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::atan;
        return autodiff::reverse::detail::Variable<T>(atan(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sinh(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::sinh;
        return autodiff::reverse::detail::Variable<T>(sinh(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> cosh(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::cosh;
        return autodiff::reverse::detail::Variable<T>(cosh(x));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> tanh(const autodiff::reverse::detail::ExprPtr<T>& x) {
        using autodiff::reverse::detail::tanh;
        return autodiff::reverse::detail::Variable<T>(tanh(x));
    }
    
    // Note: isinf, isnan, isfinite are now defined in gsAutoDiffTraits.h in std:: namespace
    
    // Helper to recursively extract innermost real value for floor/ceil
    template<typename T>
    inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
    floor_val(const T& v) { using std::floor; return floor(v); }
    
    template<typename T, typename G>
    inline GISMO_COEFF_TYPE floor_val(const autodiff::detail::Dual<T, G>& v)
    { return floor_val(v.val); }
    
    template<typename T>
    inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
    ceil_val(const T& v) { using std::ceil; return ceil(v); }
    
    template<typename T, typename G>
    inline GISMO_COEFF_TYPE ceil_val(const autodiff::detail::Dual<T, G>& v)
    { return ceil_val(v.val); }
    
    // ceil/floor/round for autodiff::detail::Dual directly (base case: arithmetic T)
    template<typename G>
    inline autodiff::detail::Dual<GISMO_COEFF_TYPE, G> ceil(const autodiff::detail::Dual<GISMO_COEFF_TYPE, G>& v) {
        using std::ceil;
        // ceil is non-differentiable, return value with zero derivative
        return autodiff::detail::Dual<GISMO_COEFF_TYPE, G>(ceil(v.val));
    }
    
    // ceil for nested dual (dual2nd_t)
    template<typename T, typename G1, typename G2>
    inline autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2> ceil(const autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2>& v) {
        // ceil is non-differentiable, return value with zero derivative
        return autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2>(ceil_val(v));
    }
    
    // ceil for var types (reverse mode)
    template <typename T>
    inline autodiff::reverse::detail::Variable<T> ceil(const autodiff::reverse::detail::Variable<T>& v) {
        using std::ceil;
        // ceil is non-differentiable, return value only with zero derivative
        return autodiff::reverse::detail::Variable<T>(ceil(autodiff::val(v)));
    }

    // floor for dual (base case: arithmetic T)
    template<typename G>
    inline autodiff::detail::Dual<GISMO_COEFF_TYPE, G> floor(const autodiff::detail::Dual<GISMO_COEFF_TYPE, G>& v) {
        using std::floor;
        // floor is non-differentiable, return value with zero derivative
        return autodiff::detail::Dual<GISMO_COEFF_TYPE, G>(floor(v.val));
    }
    
    // floor for nested dual (dual2nd_t)
    template<typename T, typename G1, typename G2>
    inline autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2> floor(const autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2>& v) {
        // floor is non-differentiable, return value with zero derivative
        return autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2>(floor_val(v));
    }
    
    // floor for var types (reverse mode)
    template <typename T>
    inline autodiff::reverse::detail::Variable<T> floor(const autodiff::reverse::detail::Variable<T>& v) {
        using std::floor;
        // floor is non-differentiable, return value only with zero derivative
        return autodiff::reverse::detail::Variable<T>(floor(autodiff::val(v)));
    }
    
    // floor for ExprPtr (reverse mode expression)
    template <typename T>
    inline T floor(const std::shared_ptr<autodiff::reverse::detail::Expr<T>>& expr) {
        using std::floor;
        // Extract the value directly from the expression pointer
        return floor(expr->val);
    }
    
    // ceil for ExprPtr (reverse mode expression)
    template <typename T>
    inline T ceil(const std::shared_ptr<autodiff::reverse::detail::Expr<T>>& expr) {
        using std::ceil;
        // Extract the value directly from the expression pointer
        return ceil(expr->val);
    }

} // namespace math
} // namespace gismo

// Also add ceil/floor to autodiff::reverse::detail namespace for ADL
namespace autodiff {
namespace reverse {
namespace detail {
    template<typename T>
    inline T ceil(const ExprPtr<T>& expr) {
        using std::ceil;
        return ceil(expr->val);
    }
    
    template<typename T>
    inline T floor(const ExprPtr<T>& expr) {
        using std::floor;
        return floor(expr->val);
    }
} // namespace detail
} // namespace reverse
} // namespace autodiff

// ============================================================================
// Variable-returning operator overloads for reverse mode
// ============================================================================
// autodiff's operators return ExprPtr<T>, which causes ternary operator
// ambiguity with Variable<T> (both types convert to each other). These
// overloads shadow autodiff's and return Variable<T> directly.
namespace autodiff {
namespace reverse {
namespace detail {

    inline Variable<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr + r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr - r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr * r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator/(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr / r.expr);
    }

    inline Variable<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& v) {
        return Variable<GISMO_COEFF_TYPE>(-v.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& v) {
        return Variable<GISMO_COEFF_TYPE>(v.expr);
    }

    inline Variable<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr + r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator+(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l + r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr - r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator-(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l - r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr * r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator*(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l * r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator/(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr / r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator/(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l / r.expr);
    }

    inline Variable<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, int r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr * static_cast<GISMO_COEFF_TYPE>(r));
    }
    inline Variable<GISMO_COEFF_TYPE> operator*(int l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(static_cast<GISMO_COEFF_TYPE>(l) * r.expr);
    }

    inline Variable<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr + r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator+(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l + r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr - r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator-(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l - r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr * r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator*(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l * r.expr);
    }
    inline Variable<GISMO_COEFF_TYPE> operator/(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l.expr / r);
    }
    inline Variable<GISMO_COEFF_TYPE> operator/(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return Variable<GISMO_COEFF_TYPE>(l / r.expr);
    }

    inline Variable<GISMO_COEFF_TYPE> conj(const Variable<GISMO_COEFF_TYPE>& x) {
        return x;
    }

} // namespace detail
} // namespace reverse
} // namespace autodiff

namespace gismo {
namespace math {

    template<typename T, typename G>
    inline T round(const autodiff::detail::Dual<T, G>& v) {
        using std::round;
        return round(v.val);
    }

    // ceil for autodiff expression types (UnaryExpr, BinaryExpr, etc.)
    // These expressions arise from operations like -log10(...)
    template<typename Op, typename ExprType>
    inline auto ceil(const autodiff::detail::UnaryExpr<Op, ExprType>& expr)
        -> decltype(std::ceil(autodiff::val(expr)))
    {
        using std::ceil;
        return ceil(autodiff::val(expr));
    }

    template<typename Op, typename L, typename R>
    inline auto ceil(const autodiff::detail::BinaryExpr<Op, L, R>& expr)
        -> decltype(std::ceil(autodiff::val(expr)))
    {
        using std::ceil;
        return ceil(autodiff::val(expr));
    }

    template<typename Op, typename ExprType>
    inline auto floor(const autodiff::detail::UnaryExpr<Op, ExprType>& expr)
        -> decltype(std::floor(autodiff::val(expr)))
    {
        using std::floor;
        return floor(autodiff::val(expr));
    }

    template<typename Op, typename L, typename R>
    inline auto floor(const autodiff::detail::BinaryExpr<Op, L, R>& expr)
        -> decltype(std::floor(autodiff::val(expr)))
    {
        using std::floor;
        return floor(autodiff::val(expr));
    }

    template<typename Op, typename ExprType>
    inline auto round(const autodiff::detail::UnaryExpr<Op, ExprType>& expr)
        -> decltype(std::round(autodiff::val(expr)))
    {
        using std::round;
        return round(autodiff::val(expr));
    }

    template<typename Op, typename L, typename R>
    inline auto round(const autodiff::detail::BinaryExpr<Op, L, R>& expr)
        -> decltype(std::round(autodiff::val(expr)))
    {
        using std::round;
        return round(autodiff::val(expr));
    }
}

} // namespace gismo

// ============================================================================
// Expression template support for autodiff scalars
// ============================================================================
// The expression system uses gismo::util::is_arithmetic (which is std::is_arithmetic
// in C++11 mode) to determine if a type should be treated as a scalar constant 
// (_expr<T, true>). We specialize std::is_arithmetic for autodiff types so they work
// with gismo's expression templates.
//
// IMPORTANT: We do NOT specialize autodiff::detail::ArithmeticTraits because that would
// cause the autodiff library to treat Dual<double, double> as a primitive arithmetic type,
// leading to circular recursion in operators. The autodiff library already knows how to
// handle Dual types via its expression template system.

// Mark autodiff types as arithmetic for gismo's expression system only.
// We use gismo::util::is_arithmetic, NOT std::is_arithmetic, to avoid
// interfering with autodiff's own traits that rely on std::is_arithmetic.
// Note: gismo::util::is_arithmetic specializations for autodiff types
// are defined centrally in src/gsCore/gsTemplateTools.h to avoid duplication.



// ============================================================================
// Comparison operators for autodiff types
// ============================================================================
// These operators are needed for mesh operations (vertex comparison, sorting)
// and other algorithms that use relational operators. They compare the values
// only (derivatives are not compared).

namespace autodiff {
namespace reverse {
namespace detail {
    // Comparison operators for Variable<T> (var_t)
    // These compare only the scalar values, ignoring derivatives
    
    template<typename T>
    inline bool operator==(const Variable<T>& l, const Variable<T>& r) {
        return l.expr->val == r.expr->val;
    }
    
    template<typename T>
    inline bool operator!=(const Variable<T>& l, const Variable<T>& r) {
        return l.expr->val != r.expr->val;
    }
    
    template<typename T>
    inline bool operator<(const Variable<T>& l, const Variable<T>& r) {
        return l.expr->val < r.expr->val;
    }
    
    template<typename T>
    inline bool operator>(const Variable<T>& l, const Variable<T>& r) {
        return l.expr->val > r.expr->val;
    }
    
    template<typename T>
    inline bool operator<=(const Variable<T>& l, const Variable<T>& r) {
        return l.expr->val <= r.expr->val;
    }
    
    template<typename T>
    inline bool operator>=(const Variable<T>& l, const Variable<T>& r) {
        return l.expr->val >= r.expr->val;
    }
    
    // Mixed comparisons: Variable<T> with arithmetic types
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator==(const Variable<T>& l, U r) {
        return l.expr->val == r;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator==(U l, const Variable<T>& r) {
        return l == r.expr->val;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator!=(const Variable<T>& l, U r) {
        return l.expr->val != r;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator!=(U l, const Variable<T>& r) {
        return l != r.expr->val;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator<(const Variable<T>& l, U r) {
        return l.expr->val < r;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator<(U l, const Variable<T>& r) {
        return l < r.expr->val;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator>(const Variable<T>& l, U r) {
        return l.expr->val > r;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator>(U l, const Variable<T>& r) {
        return l > r.expr->val;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator<=(const Variable<T>& l, U r) {
        return l.expr->val <= r;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator<=(U l, const Variable<T>& r) {
        return l <= r.expr->val;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator>=(const Variable<T>& l, U r) {
        return l.expr->val >= r;
    }
    
    template<typename T, typename U>
    inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
    operator>=(U l, const Variable<T>& r) {
        return l >= r.expr->val;
    }

} // namespace detail
} // namespace reverse
} // namespace autodiff


