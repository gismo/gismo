/** @file NullExpression.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

namespace gismo
{
namespace Expr
{

    template <typename S, short_t _order>
    class NullExpression : public BaseExpression<NullExpression<S,_order>, T, _order, true>
    {
        using Base = BaseExpression<NullExpression<S,_order>, T, _order, true>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr int order = Base::order;
        static constexpr bool isConstant = Base::isConstant;

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

        gsMatrix<Scalar> eval(const index_t k) const
        {
            GISMO_ERROR("The null expression should not be evaluated");
        }

    protected:

        explicit NullExpression(const std::array<size_t, order>& input_sizes)
        :
        sizes_(input_sizes)
        {}

        NullExpression(const NullExpression&) = default;
        NullExpression& operator=(const NullExpression&) = default;
        ~NullExpression() = default;

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