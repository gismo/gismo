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

namespace gismo
{
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
    using autodiff::var;
    inline var abs(const var& x) { using std::abs; return abs(x); }
    inline var sqrt(const var& x) { using std::sqrt; return sqrt(x); }
    inline var pow(const var& x, const var& y) { using std::pow; return pow(x, y); }
    inline var pow(const var& x, double y) { using std::pow; return pow(x, y); }
    inline var exp(const var& x) { using std::exp; return exp(x); }
    inline var log(const var& x) { using std::log; return log(x); }
    inline var sin(const var& x) { using std::sin; return sin(x); }
    inline var cos(const var& x) { using std::cos; return cos(x); }
    inline var tan(const var& x) { using std::tan; return tan(x); }
    inline var asin(const var& x) { using std::asin; return asin(x); }
    inline var acos(const var& x) { using std::acos; return acos(x); }
    inline var atan(const var& x) { using std::atan; return atan(x); }
    inline var sinh(const var& x) { using std::sinh; return sinh(x); }
    inline var cosh(const var& x) { using std::cosh; return cosh(x); }
    inline var tanh(const var& x) { using std::tanh; return tanh(x); }
    
    inline bool isinf(const var& x) { return std::isinf(autodiff::val(x)); }
    inline bool isnan(const var& x) { return std::isnan(autodiff::val(x)); }
    inline bool isfinite(const var& x) { return std::isfinite(autodiff::val(x)); }
}
}
