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
    // For forward mode (dual_t): returns .val member
    // For reverse mode (var_t): returns expr->val
    // For arithmetic types: returns the value itself
    template<typename T, typename G>
    inline T gismo_val(const autodiff::detail::Dual<T, G>& v) { return v.val; }
    
    template<typename T>
    inline T gismo_val(const autodiff::reverse::detail::Variable<T>& v) { return v.expr->val; }
    
    template<typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
    inline T gismo_val(const T& v) { return v; }

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

        void propagatex(const autodiff::reverse::detail::ExprPtr<T>& wprime) override {
            // Similar for higher-order, if implemented
            auto adjoint_b = solver.solve_transpose(wprime);  // Assuming overloaded for ExprPtr
            r->propagatex(adjoint_b);
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
    
    // ceil/floor/round for autodiff::detail::Dual directly
    template<typename T, typename G>
    inline autodiff::detail::Dual<T, G> ceil(const autodiff::detail::Dual<T, G>& v) {
        using std::ceil;
        // ceil is non-differentiable, return value with zero derivative
        return autodiff::detail::Dual<T, G>(ceil(v.val));
    }
    
    // ceil for var types (reverse mode)
    template <typename T>
    inline autodiff::reverse::detail::Variable<T> ceil(const autodiff::reverse::detail::Variable<T>& v) {
        using std::ceil;
        // ceil is non-differentiable, return value only with zero derivative
        return autodiff::reverse::detail::Variable<T>(ceil(autodiff::val(v)));
    }

    template<typename T, typename G>
    inline autodiff::detail::Dual<T, G> floor(const autodiff::detail::Dual<T, G>& v) {
        using std::floor;
        // floor is non-differentiable, return value with zero derivative
        return autodiff::detail::Dual<T, G>(floor(v.val));
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
}
