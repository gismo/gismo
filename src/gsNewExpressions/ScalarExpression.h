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

    template <class T, bool _isConstant = true>
    class ScalarExpression : public BaseExpression<ScalarExpression<T>, T, 0, _isConstant>
    {
        using Base = BaseExpression<ScalarExpression<T>, T, 0, _isConstant>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr int order = Base::order;
        static constexpr bool isConstant = Base::isConstant;

    private:
        Scalar m_value; // Store scalar value directly
    public:
        // Explicit constructor for a scalar value
        explicit ScalarExpression(Scalar value = 0.0)
        :
        BaseExpression<ScalarExpression<Scalar>,Scalar, 0>({}),    // Pass empty array for sizes_ for order 0 (scalar)
        m_value(value)
        {}

        gsMatrix<Scalar> eval(const index_t k) const
        {
            // For a scalar, k might be ignored or used as a parameter if it's a field.
            // Returns a 1x1 gsMatrix containing the scalar value.
            gsMatrix<Scalar> result(1, 1);
            result(0, 0) = m_value;
            return result;
        }

    };
}//namespace Expr
}//namespace gismo