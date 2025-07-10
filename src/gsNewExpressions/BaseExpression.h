/** @file BaseExpression.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsMatrix/gsMatrix.h>

namespace gismo
{
namespace Expr
{

    template <typename E>
    struct ExpressionTraits
    {
    public:
        typedef real_t Scalar;//todo
        typedef const E Nested_t;
    };

    // IsConstant: Flag that indicates if the expression is constant, e.g., its derivatives are zero
    template <typename E, typename S, short_t _order, bool _isConstant = true>
    class BaseExpression
    {
    public:
        typedef S Scalar; // Publicly expose Scalar type

        // The order of the expression: 0: scalar, 1: vector, 2: matrix, etc.
        static constexpr int order = _order;
        // Flag that indicates if the expression is constant, e.g., its derivatives are zero
        static constexpr bool isConstant = _isConstant;

        // The size of the expression, which is an array of sizes for each dimension
        const std::array<size_t, _order> sizes_;
        // Access to sizes is direct
        const std::array<size_t, _order>& sizes() const
        {
            return sizes_;
        }

    public:

        gsMatrix<Scalar> eval(const index_t k) const { return static_cast<E const&>(*this).eval(k); }

    protected:

        explicit BaseExpression(const std::array<size_t, order>& input_sizes)
        :
        sizes_(input_sizes)
        {}

        BaseExpression(const BaseExpression&) = default;
        BaseExpression& operator=(const BaseExpression&) = default;
        ~BaseExpression() = default;

        // Helper to calculate total elements from sizes (used internally by derived classes)
        static size_t tensorSize(const std::array<size_t, order>& dims)
        {
            if (order == 0) return 1; // Scalar has 1 element
            size_t total = 1;
            for (size_t dim_size : dims)
            {
                total *= dim_size;
            }
            return total;
        }

    };
}//namespace Expr
}//namespace gismo