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

using namespace autodiff;

using autodiff::reverse::detail::BinaryExpr;
using autodiff::reverse::detail::ExprPtr;
#include <type_traits>
#include <utility>
#include <gsCore/gsForwardDeclarations.h>

// If the helpers are not yet defined, provide them here too
#ifndef GISMO_AUTODIFF_TO_STRING_OVERLOADS
#define GISMO_AUTODIFF_TO_STRING_OVERLOADS
namespace gismo {
namespace util {

template <class C, class = void>
struct has_val : std::false_type {};
template <class C>
struct has_val<C, std::void_t<decltype(std::declval<const C&>().val)>> : std::true_type {};

template <typename E, std::enable_if_t<has_val<E>::value && !std::is_same<std::decay_t<E>, real_t>::value, int> = 0>
inline std::string to_string(const E &value)
{
    return util::to_string(static_cast<real_t>(value.val));
}

template <typename T>
inline std::string to_string(const autodiff::reverse::detail::ExprPtr<T> &expr)
{
    if (expr) return util::to_string(static_cast<real_t>(expr->val));
    return std::string("null");
}

} // namespace util
} // namespace gismo
#endif

namespace gismo
{
    template<typename T, typename Solver>
    struct SolveExpr : BinaryExpr<T> {
        // Using declarations for base class
        using BinaryExpr<T>::l;  // A (matrix input, as ExprPtr<T> or extended)
        using BinaryExpr<T>::r;  // b (vector input)
        using BinaryExpr<T>::val;

        Solver solver;  // Template parameter for solver type

        SolveExpr(const T& v, const ExprPtr<T>& a, const ExprPtr<T>& bb, Solver s)
            : BinaryExpr<T>(v, a, bb), solver(std::move(s)) {}

        void propagate(const T& wprime) override {
            // Assuming solver.solve_transpose(wprime) computes A^T v = wprime
            T adjoint_b = solver.solve_transpose(wprime);
            r->propagate(adjoint_b);
            // Handle A if variable (e.g., propagate to l if needed)
        }

        void propagatex(const ExprPtr<T>& wprime) override {
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
    ExprPtr<T> solve(const ExprPtr<T>& a, const ExprPtr<T>& b, Solver solver) 
    {
        T x_val = solver.solve(b->val);
        return std::make_shared<SolveExpr<T, Solver>>(x_val, a, b, std::move(solver));
    }
}
