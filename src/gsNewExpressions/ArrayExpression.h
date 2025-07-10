/** @file BaseExpression.h

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
    // Expression to wrap an existing expression (e.g., a gsMatrix, gsVector, or another expression) in an array
    template <typename WrappedExpr>
    class ArrayExpression : public BaseExpression<ArrayExpression<WrappedExpr>, typename WrappedExpr::Scalar, WrappedExpr::order, WrappedExpr::isConstant>
    {
    public:
        typedef typename WrappedExpr::Scalar Scalar;
        static constexpr int order = WrappedExpr::order;
        static constexpr bool isConstant = WrappedExpr::isConstant; // ArrayExpression inherits constant-ness

        explicit ArrayExpression(const WrappedExpr& expr)
            : BaseExpression<ArrayExpression<WrappedExpr>, Scalar, order>(expr.sizes()), // Keep original sizes
            wrapped_expr_(expr)
        {}

        gsMatrix<Scalar> eval(const index_t k) const
        {
            return wrapped_expr_.eval(k); // Simply evaluate the wrapped expression
        }
        // value() method is not strictly necessary for ArrayExpression unless it wraps a scalar
        // For consistency, if WrappedExpr has value(), ArrayExpression can expose it
        Scalar value() const
        {
            return wrapped_expr_.value(); // Forward value() call to the wrapped expression
        }

    private:
        const WrappedExpr& wrapped_expr_;
    };

// Helper function to create ArrayExpression
template <typename Expr>
ArrayExpression<Expr> array(const Expr& expr)
{
    return ArrayExpression<Expr>(expr);
}

}//namespace Expr
}//namespace gismo