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

// Forward declarations for reverse mode
namespace autodiff { namespace reverse { namespace detail {
    template <typename T> class Variable;
    template <typename T> class Expr;
}}}
#include <memory>  // for std::shared_ptr


namespace gismo
{
    // Unified val() accessor that works for both forward and reverse autodiff types
    // For forward mode (dual_t): returns .val member (recursively for nested duals)
    // For reverse mode (var_t): returns expr->val
    // For arithmetic types: returns the value itself
    
    // Primary template for arithmetic types
    template<typename T>
    inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
    gismo_val(const T& v) { return v; }
    
    // Forward mode dual with arithmetic inner type (base case)
    template<typename G>
    inline GISMO_COEFF_TYPE gismo_val(const autodiff::detail::Dual<GISMO_COEFF_TYPE, G>& v)
    { return v.val; }
    
    // Forward mode dual with dual inner type (recursive case - dual2nd_t)
    template<typename T, typename G1, typename G2>
    inline GISMO_COEFF_TYPE gismo_val(const autodiff::detail::Dual<autodiff::detail::Dual<T, G1>, G2>& v)
    { return gismo_val(v.val); }
    
    // Reverse mode variable
    template<typename T>
    inline T gismo_val(const autodiff::reverse::detail::Variable<T>& v) { return v.expr->val; }

    template<typename T, typename Solver>
    struct SolveExpr : autodiff::reverse::detail::BinaryExpr<T> {
        // Using declarations for base class
        using Base = autodiff::reverse::detail::BinaryExpr<T>;
        using Base::l;  // A (matrix input, as ExprPtr<T> or extended)
        using Base::r;  // b (vector input)
        using Base::val;

        Solver solver;  // Template parameter for solver type

        SolveExpr(const T& v, const autodiff::reverse::detail::ExprPtr<T>& a, const autodiff::reverse::detail::ExprPtr<T>& bb, Solver s)
            : Base(v, a, bb), solver(std::move(s)) {}

        void propagate(const T& wprime) override {
            // Assuming solver.solve_transpose(wprime) computes A^T v = wprime
            T adjoint_b = solver.solve_transpose(wprime);
            r->propagate(adjoint_b);
            // Handle A if variable (e.g., propagate to l if needed)
        }

        void update() override {
            l->update();
            r->update();
            val = solver.solve(r->val);  // Forward compute
        }
    };

    // Convenience function
    template<typename T, typename Solver>
    autodiff::reverse::detail::ExprPtr<T> solve(const autodiff::reverse::detail::ExprPtr<T>& a, const autodiff::reverse::detail::ExprPtr<T>& b, Solver solver) 
    {
        T x_val = solver.solve(b->val);
        return std::make_shared<SolveExpr<T, Solver>>(x_val, a, b, std::move(solver));
    }
}

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
    
    // Math functions for ExprPtrWrapper (needed when expressions produce ExprPtrWrapper)
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> abs(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::abs;
        return autodiff::reverse::detail::Variable<T>(abs(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sqrt(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::sqrt;
        return autodiff::reverse::detail::Variable<T>(sqrt(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::ExprPtrWrapper<T>& x, T y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x.ptr, y));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::ExprPtrWrapper<T>& x, int y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x.ptr, static_cast<T>(y)));
    }
    // pow(Variable, ExprPtrWrapper) and pow(ExprPtrWrapper, Variable)
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::Variable<T>& x, const autodiff::reverse::detail::ExprPtrWrapper<T>& y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x.expr, y.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::ExprPtrWrapper<T>& x, const autodiff::reverse::detail::Variable<T>& y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x.ptr, y.expr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> pow(const autodiff::reverse::detail::ExprPtrWrapper<T>& x, const autodiff::reverse::detail::ExprPtrWrapper<T>& y) {
        using autodiff::reverse::detail::pow;
        return autodiff::reverse::detail::Variable<T>(pow(x.ptr, y.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> exp(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::exp;
        return autodiff::reverse::detail::Variable<T>(exp(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> log(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::log;
        return autodiff::reverse::detail::Variable<T>(log(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sin(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::sin;
        return autodiff::reverse::detail::Variable<T>(sin(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> cos(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::cos;
        return autodiff::reverse::detail::Variable<T>(cos(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> tan(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::tan;
        return autodiff::reverse::detail::Variable<T>(tan(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> asin(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::asin;
        return autodiff::reverse::detail::Variable<T>(asin(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> acos(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::acos;
        return autodiff::reverse::detail::Variable<T>(acos(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> atan(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::atan;
        return autodiff::reverse::detail::Variable<T>(atan(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> sinh(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::sinh;
        return autodiff::reverse::detail::Variable<T>(sinh(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> cosh(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::cosh;
        return autodiff::reverse::detail::Variable<T>(cosh(x.ptr));
    }
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> tanh(const autodiff::reverse::detail::ExprPtrWrapper<T>& x) {
        using autodiff::reverse::detail::tanh;
        return autodiff::reverse::detail::Variable<T>(tanh(x.ptr));
    }
    
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
// Ternary operator fix for autodiff::var (reverse mode)
// ============================================================================
// Problem: Variable<T> has both a constructor from ExprPtr<T> AND an implicit
// conversion to ExprPtr<T>. When binary operators return ExprPtr<T> and are
// used in ternary expressions, C++ cannot decide which conversion to use.
// Solution: Wrap ExprPtr to provide controlled implicit conversions.
// ExprPtrWrapper is defined in gsAutoDiffTraits.h (needed by isnan/isinf/isfinite)

// ============================================================================
// Operator overloads for ExprPtrWrapper
// ============================================================================
namespace autodiff {
namespace reverse {
namespace detail {
    
    // ========================================================================
    // Binary operator overloads: Variable op Variable -> ExprPtrWrapper
    // Note: These MUST be explicit (not templated) to avoid ambiguity with
    // autodiff's own templated operators that return ExprPtr<T>
    // Using GISMO_COEFF_TYPE to support both double and long double
    // ========================================================================
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr + r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr - r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr * r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator/(const Variable<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr / r.expr);
    }
    
    // ========================================================================
    // Unary operator overloads - must also be explicit for same reason
    // ========================================================================
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& v) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(-v.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& v) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(+v.expr);
    }
    
    // ========================================================================
    // Mixed operations: Variable op ExprPtr -> ExprPtrWrapper
    // Must be explicit to avoid ambiguity with autodiff's operators
    // ========================================================================
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr + r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr - r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr * r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator/(const Variable<GISMO_COEFF_TYPE>& l, const ExprPtr<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr / r);
    }
    
    // ========================================================================
    // Mixed operations: ExprPtr op Variable -> ExprPtrWrapper
    // ========================================================================
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator+(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l + r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator-(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l - r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l * r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator/(const ExprPtr<GISMO_COEFF_TYPE>& l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l / r.expr);
    }
    
    // ========================================================================
    // Operations with scalars: Variable op GISMO_COEFF_TYPE -> ExprPtrWrapper
    // Must be explicit to avoid conflicts with autodiff
    // ========================================================================
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator+(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr + r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator+(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l + r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator-(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr - r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator-(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l - r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr * r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l * r.expr);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator/(const Variable<GISMO_COEFF_TYPE>& l, GISMO_COEFF_TYPE r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr / r);
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator/(GISMO_COEFF_TYPE l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l / r.expr);
    }
    
    // With int (cast to GISMO_COEFF_TYPE)
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(const Variable<GISMO_COEFF_TYPE>& l, int r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(l.expr * static_cast<GISMO_COEFF_TYPE>(r));
    }
    
    inline ExprPtrWrapper<GISMO_COEFF_TYPE> operator*(int l, const Variable<GISMO_COEFF_TYPE>& r) {
        return ExprPtrWrapper<GISMO_COEFF_TYPE>(static_cast<GISMO_COEFF_TYPE>(l) * r.expr);
    }
    
    // ========================================================================
    // ExprPtrWrapper op scalar
    // ========================================================================
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(const ExprPtrWrapper<T>& l, T r) {
        return ExprPtrWrapper<T>(l.ptr + r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(T l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l + r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(const ExprPtrWrapper<T>& l, T r) {
        return ExprPtrWrapper<T>(l.ptr - r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(T l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l - r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(const ExprPtrWrapper<T>& l, T r) {
        return ExprPtrWrapper<T>(l.ptr * r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(T l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l * r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(const ExprPtrWrapper<T>& l, T r) {
        return ExprPtrWrapper<T>(l.ptr / r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(T l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l / r.ptr);
    }
    
    // ========================================================================
    // Mixed operations: Variable op ExprPtrWrapper -> ExprPtrWrapper
    // ========================================================================
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.expr + r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.expr - r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.expr * r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.expr / r.ptr);
    }
    
    // ========================================================================
    // Mixed operations: ExprPtrWrapper op Variable -> ExprPtrWrapper
    // ========================================================================
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return ExprPtrWrapper<T>(l.ptr + r.expr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return ExprPtrWrapper<T>(l.ptr - r.expr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return ExprPtrWrapper<T>(l.ptr * r.expr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return ExprPtrWrapper<T>(l.ptr / r.expr);
    }
    
    // ========================================================================
    // Mixed operations: ExprPtrWrapper op ExprPtrWrapper -> ExprPtrWrapper
    // ========================================================================
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.ptr + r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.ptr - r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.ptr * r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l.ptr / r.ptr);
    }
    
    // ========================================================================
    // Mixed operations: ExprPtrWrapper op ExprPtr -> ExprPtrWrapper
    // ========================================================================
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return ExprPtrWrapper<T>(l.ptr + r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return ExprPtrWrapper<T>(l.ptr - r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return ExprPtrWrapper<T>(l.ptr * r);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return ExprPtrWrapper<T>(l.ptr / r);
    }
    
    // ========================================================================
    // Mixed operations: ExprPtr op ExprPtrWrapper -> ExprPtrWrapper
    // ========================================================================
    
    template<typename T>
    inline ExprPtrWrapper<T> operator+(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l + r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator-(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l - r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator*(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l * r.ptr);
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> operator/(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return ExprPtrWrapper<T>(l / r.ptr);
    }
    
    // ========================================================================
    // Utility functions for ExprPtrWrapper
    // ========================================================================
    
    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, const ExprPtrWrapper<T>& w) {
        return os << w.ptr->val;
    }
    
    // ========================================================================
    // Compound assignment operators: Variable op= ExprPtrWrapper
    // These are needed because Variable::operator+= takes ExprPtr, not Wrapper
    // ========================================================================
    
    template<typename T>
    inline Variable<T>& operator+=(Variable<T>& l, const ExprPtrWrapper<T>& r) {
        l += r.ptr;
        return l;
    }
    
    template<typename T>
    inline Variable<T>& operator-=(Variable<T>& l, const ExprPtrWrapper<T>& r) {
        l -= r.ptr;
        return l;
    }
    
    template<typename T>
    inline Variable<T>& operator*=(Variable<T>& l, const ExprPtrWrapper<T>& r) {
        l *= r.ptr;
        return l;
    }
    
    template<typename T>
    inline Variable<T>& operator/=(Variable<T>& l, const ExprPtrWrapper<T>& r) {
        l /= r.ptr;
        return l;
    }
    
} // namespace detail
} // namespace reverse
} // namespace autodiff

// floor/ceil for ExprPtrWrapper
namespace gismo {
namespace math {
    template<typename T>
    inline T floor(const autodiff::reverse::detail::ExprPtrWrapper<T>& w) {
        using std::floor;
        return floor(w.ptr->val);
    }
    
    template<typename T>
    inline T ceil(const autodiff::reverse::detail::ExprPtrWrapper<T>& w) {
        using std::ceil;
        return ceil(w.ptr->val);
    }
    
    // min/max with ExprPtrWrapper (for gismo::math:: qualified calls)
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> min(
        const autodiff::reverse::detail::Variable<T>& a, 
        const autodiff::reverse::detail::ExprPtrWrapper<T>& b) {
        return a.expr->val < b.ptr->val ? a : autodiff::reverse::detail::Variable<T>(b.ptr);
    }
    
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> min(
        const autodiff::reverse::detail::ExprPtrWrapper<T>& a,
        const autodiff::reverse::detail::Variable<T>& b) {
        return a.ptr->val < b.expr->val ? autodiff::reverse::detail::Variable<T>(a.ptr) : b;
    }
    
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> max(
        const autodiff::reverse::detail::Variable<T>& a, 
        const autodiff::reverse::detail::ExprPtrWrapper<T>& b) {
        return a.expr->val > b.ptr->val ? a : autodiff::reverse::detail::Variable<T>(b.ptr);
    }
    
    template<typename T>
    inline autodiff::reverse::detail::Variable<T> max(
        const autodiff::reverse::detail::ExprPtrWrapper<T>& a,
        const autodiff::reverse::detail::Variable<T>& b) {
        return a.ptr->val > b.expr->val ? autodiff::reverse::detail::Variable<T>(a.ptr) : b;
    }
}
}

// Note: isfinite/isinf/isnan and val() for ExprPtrWrapper are defined in gsAutoDiffTraits.h

// Math functions for ExprPtrWrapper (needed by Eigen)
namespace autodiff {
namespace reverse {
namespace detail {
    // Math functions - keep templated as they don't conflict with autodiff
    template<typename T>
    inline ExprPtrWrapper<T> sqrt(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::sqrt;
        return ExprPtrWrapper<T>(sqrt(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> abs(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::abs;
        return ExprPtrWrapper<T>(abs(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> conj(const ExprPtrWrapper<T>& w) {
        // conj for real types is identity
        return w;
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> exp(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::exp;
        return ExprPtrWrapper<T>(exp(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> log(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::log;
        return ExprPtrWrapper<T>(log(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> log10(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::log10;
        return ExprPtrWrapper<T>(log10(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> sin(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::sin;
        return ExprPtrWrapper<T>(sin(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> cos(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::cos;
        return ExprPtrWrapper<T>(cos(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> tan(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::tan;
        return ExprPtrWrapper<T>(tan(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> asin(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::asin;
        return ExprPtrWrapper<T>(asin(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> acos(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::acos;
        return ExprPtrWrapper<T>(acos(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> atan(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::atan;
        return ExprPtrWrapper<T>(atan(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> pow(const ExprPtrWrapper<T>& base, const ExprPtrWrapper<T>& exp) {
        using autodiff::reverse::detail::pow;
        return ExprPtrWrapper<T>(pow(base.ptr, exp.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> pow(const ExprPtrWrapper<T>& base, T exp) {
        using autodiff::reverse::detail::pow;
        return ExprPtrWrapper<T>(pow(base.ptr, exp));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> erf(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::erf;
        return ExprPtrWrapper<T>(erf(w.ptr));
    }
    
    template<typename T>
    inline ExprPtrWrapper<T> erfc(const ExprPtrWrapper<T>& w) {
        using autodiff::reverse::detail::erfc;
        return ExprPtrWrapper<T>(erfc(w.ptr));
    }
    
    // val() extractor for ExprPtrWrapper (needed by autodiff utilities)
    template<typename T>
    inline T val(const ExprPtrWrapper<T>& w) {
        return w.ptr->val;
    }
    
    // min/max functions for Variable and ExprPtrWrapper combinations
    template<typename T>
    inline Variable<T> min(const Variable<T>& a, const ExprPtrWrapper<T>& b) {
        return a.expr->val < b.ptr->val ? a : Variable<T>(b.ptr);
    }
    
    template<typename T>
    inline Variable<T> min(const ExprPtrWrapper<T>& a, const Variable<T>& b) {
        return a.ptr->val < b.expr->val ? Variable<T>(a.ptr) : b;
    }
    
    template<typename T>
    inline Variable<T> max(const Variable<T>& a, const ExprPtrWrapper<T>& b) {
        return a.expr->val > b.ptr->val ? a : Variable<T>(b.ptr);
    }
    
    template<typename T>
    inline Variable<T> max(const ExprPtrWrapper<T>& a, const Variable<T>& b) {
        return a.ptr->val > b.expr->val ? Variable<T>(a.ptr) : b;
    }
    
    // Comparison operators for ExprPtrWrapper
    template<typename T>
    inline bool operator<(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return l.ptr->val < r.ptr->val;
    }
    
    template<typename T>
    inline bool operator<=(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return l.ptr->val <= r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return l.ptr->val > r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>=(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return l.ptr->val >= r.ptr->val;
    }
    
    template<typename T>
    inline bool operator==(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return l.ptr->val == r.ptr->val;
    }
    
    template<typename T>
    inline bool operator!=(const ExprPtrWrapper<T>& l, const ExprPtrWrapper<T>& r) {
        return l.ptr->val != r.ptr->val;
    }
    
    // Mixed comparisons: ExprPtrWrapper vs Variable
    template<typename T>
    inline bool operator<(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return l.ptr->val < r.expr->val;
    }
    
    template<typename T>
    inline bool operator<(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return l.expr->val < r.ptr->val;
    }
    
    template<typename T>
    inline bool operator<=(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return l.ptr->val <= r.expr->val;
    }
    
    template<typename T>
    inline bool operator<=(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return l.expr->val <= r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return l.ptr->val > r.expr->val;
    }
    
    template<typename T>
    inline bool operator>(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return l.expr->val > r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>=(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return l.ptr->val >= r.expr->val;
    }
    
    template<typename T>
    inline bool operator>=(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return l.expr->val >= r.ptr->val;
    }
    
    template<typename T>
    inline bool operator==(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return l.ptr->val == r.expr->val;
    }
    
    template<typename T>
    inline bool operator==(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return l.expr->val == r.ptr->val;
    }
    
    template<typename T>
    inline bool operator!=(const ExprPtrWrapper<T>& l, const Variable<T>& r) {
        return l.ptr->val != r.expr->val;
    }
    
    template<typename T>
    inline bool operator!=(const Variable<T>& l, const ExprPtrWrapper<T>& r) {
        return l.expr->val != r.ptr->val;
    }
    
    // Comparison with scalars - templated
    template<typename T>
    inline bool operator<(const ExprPtrWrapper<T>& l, T r) {
        return l.ptr->val < r;
    }
    
    template<typename T>
    inline bool operator<(T l, const ExprPtrWrapper<T>& r) {
        return l < r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>(const ExprPtrWrapper<T>& l, T r) {
        return l.ptr->val > r;
    }
    
    template<typename T>
    inline bool operator>(T l, const ExprPtrWrapper<T>& r) {
        return l > r.ptr->val;
    }
    
    template<typename T>
    inline bool operator<=(const ExprPtrWrapper<T>& l, int r) {
        return l.ptr->val <= static_cast<T>(r);
    }
    
    template<typename T>
    inline bool operator>=(const ExprPtrWrapper<T>& l, int r) {
        return l.ptr->val >= static_cast<T>(r);
    }
    
    template<typename T>
    inline bool operator<(const ExprPtrWrapper<T>& l, int r) {
        return l.ptr->val < static_cast<T>(r);
    }
    
     template<typename T>
     inline bool operator>(const ExprPtrWrapper<T>& l, int r) {
         return l.ptr->val > static_cast<T>(r);
     }
     
     // Comparison: ExprPtrWrapper vs scalar types - templated for double and other scalars
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator>=(const ExprPtrWrapper<T>& l, U r) {
         return l.ptr->val >= static_cast<T>(r);
     }
     
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator<=(const ExprPtrWrapper<T>& l, U r) {
         return l.ptr->val <= static_cast<T>(r);
     }
     
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator<(const ExprPtrWrapper<T>& l, U r) {
         return l.ptr->val < static_cast<T>(r);
     }
     
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator>(const ExprPtrWrapper<T>& l, U r) {
         return l.ptr->val > static_cast<T>(r);
     }

     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator>=(U l, const ExprPtrWrapper<T>& r) {
         return static_cast<T>(l) >= r.ptr->val;
     }
     
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator<=(U l, const ExprPtrWrapper<T>& r) {
         return static_cast<T>(l) <= r.ptr->val;
     }
     
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator<(U l, const ExprPtrWrapper<T>& r) {
         return static_cast<T>(l) < r.ptr->val;
     }
     
     template<typename T, typename U>
     inline typename std::enable_if<std::is_arithmetic<U>::value, bool>::type
     operator>(U l, const ExprPtrWrapper<T>& r) {
         return static_cast<T>(l) > r.ptr->val;
     }
     
     // Comparison: ExprPtr vs ExprPtrWrapper - templated
    template<typename T>
    inline bool operator<(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return l->val < r.ptr->val;
    }
    
    template<typename T>
    inline bool operator<(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return l.ptr->val < r->val;
    }
    
    template<typename T>
    inline bool operator>(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return l->val > r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return l.ptr->val > r->val;
    }
    
    template<typename T>
    inline bool operator<=(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return l->val <= r.ptr->val;
    }
    
    template<typename T>
    inline bool operator<=(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return l.ptr->val <= r->val;
    }
    
    template<typename T>
    inline bool operator>=(const ExprPtr<T>& l, const ExprPtrWrapper<T>& r) {
        return l->val >= r.ptr->val;
    }
    
    template<typename T>
    inline bool operator>=(const ExprPtrWrapper<T>& l, const ExprPtr<T>& r) {
        return l.ptr->val >= r->val;
    }
}
}
}

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

// ============================================================================
// Matrix multiplication for reverse mode autodiff (var_t)
// ============================================================================
// Eigen's optimized matrix multiplication doesn't preserve reverse AD
// derivatives. This function performs matrix multiplication that properly
// tracks derivatives through the computation graph.
// ============================================================================

namespace gismo {

/**
 * @brief Matrix multiplication for gsMatrix<var_t> that preserves reverse AD derivatives
 *
 * This function performs C = A * B where A and B contain var_t (reverse mode AD) types.
 * Unlike Eigen's optimized multiplication, this properly tracks dependencies in the
 * autodiff computation graph.
 *
 * @param A Left matrix (m x n)
 * @param B Right matrix (n x p)
 * @return Result matrix C (m x p) with proper derivative tracking
 */
template<typename T>
inline gsMatrix<T> autodiffMatMul(const gsMatrix<T>& A, const gsMatrix<T>& B) {
    static_assert(std::is_same<T, var_t>::value,
                  "autodiffMatMul is only needed for var_t types");

    GISMO_ASSERT(A.cols() == B.rows(), "Matrix dimensions mismatch in autodiffMatMul");

    gsMatrix<T> C(A.rows(), B.cols());

    // Explicit loop-based multiplication that uses var_t operations
    // This ensures the autodiff graph is properly built
    for (index_t i = 0; i < A.rows(); ++i) {
        for (index_t j = 0; j < B.cols(); ++j) {
            T sum = T(0.0);
            for (index_t k = 0; k < A.cols(); ++k) {
                sum = sum + A(i, k) * B(k, j);
            }
            C(i, j) = sum;
        }
    }

    return C;
}

/**
 * @brief Matrix-vector multiplication for gsMatrix<var_t>
 */
template<typename T>
inline gsMatrix<T> autodiffMatMul(const gsMatrix<T>& A, const gsVector<T>& B) {
    static_assert(std::is_same<T, var_t>::value,
                  "autodiffMatMul is only needed for var_t types");

    GISMO_ASSERT(A.cols() == B.rows(), "Matrix dimensions mismatch in autodiffMatMul");

    gsMatrix<T> C(A.rows(), 1);

    for (index_t i = 0; i < A.rows(); ++i) {
        T sum = T(0.0);
        for (index_t k = 0; k < A.cols(); ++k) {
            sum = sum + A(i, k) * B(k);
        }
        C(i, 0) = sum;
    }

    return C;
}

/**
 * @brief Vector-matrix multiplication for gsMatrix<var_t>
 */
template<typename T>
inline gsMatrix<T> autodiffMatMul(const gsVector<T>& A, const gsMatrix<T>& B) {
    static_assert(std::is_same<T, var_t>::value,
                  "autodiffMatMul is only needed for var_t types");

    GISMO_ASSERT(A.cols() == B.rows(), "Matrix dimensions mismatch in autodiffMatMul");

    gsMatrix<T> C(1, B.cols());

    for (index_t j = 0; j < B.cols(); ++j) {
        T sum = T(0.0);
        for (index_t k = 0; k < A.cols(); ++k) {
            sum = sum + A(k) * B(k, j);
        }
        C(0, j) = sum;
    }

    return C;
}

/**
 * @brief Matrix transpose-multiplication for gsMatrix<var_t>: C = A^T * B
 */
template<typename T>
inline gsMatrix<T> autodiffMatMulT(const gsMatrix<T>& A, const gsMatrix<T>& B) {
    static_assert(std::is_same<T, var_t>::value,
                  "autodiffMatMulT is only needed for var_t types");

    GISMO_ASSERT(A.rows() == B.rows(), "Matrix dimensions mismatch in autodiffMatMulT");

    gsMatrix<T> C(A.cols(), B.cols());

    for (index_t i = 0; i < A.cols(); ++i) {
        for (index_t j = 0; j < B.cols(); ++j) {
            T sum = T(0.0);
            for (index_t k = 0; k < A.rows(); ++k) {
                sum = sum + A(k, i) * B(k, j);
            }
            C(i, j) = sum;
        }
    }

    return C;
}

/**
 * @brief Matrix multiplication-transpose for gsMatrix<var_t>: C = A * B^T
 */
template<typename T>
inline gsMatrix<T> autodiffMatMulTT(const gsMatrix<T>& A, const gsMatrix<T>& B) {
    static_assert(std::is_same<T, var_t>::value,
                  "autodiffMatMulTT is only needed for var_t types");

    GISMO_ASSERT(A.cols() == B.cols(), "Matrix dimensions mismatch in autodiffMatMulTT");

    gsMatrix<T> C(A.rows(), B.rows());

    for (index_t i = 0; i < A.rows(); ++i) {
        for (index_t j = 0; j < B.rows(); ++j) {
            T sum = T(0.0);
            for (index_t k = 0; k < A.cols(); ++k) {
                sum = sum + A(i, k) * B(j, k);
            }
            C(i, j) = sum;
        }
    }

    return C;
}

} // namespace gismo
