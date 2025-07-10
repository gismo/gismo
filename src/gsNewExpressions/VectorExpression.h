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
    class VectorExpression : public BaseExpression<VectorExpression<T>, T, 1, _isConstant>
    {
        using Base = BaseExpression<VectorExpression<T>, T, 1, _isConstant>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr int order = Base::order;
        static constexpr bool isConstant = Base::isConstant;

    private:
        gsVector<Scalar> m_value;
    public:

        // Constructor for a matrix filled with the same value
        VectorExpression(size_t rows, T value = 0.0)
        :
        BaseExpression<VectorExpression<Scalar>, Scalar, 1>({rows}),
        m_value(gsVector<T>::Constant(rows,value))
        {
        }

        // Explicit constructor for a gsVector
        explicit VectorExpression(const gsVector<Scalar> & value)
        :
        BaseExpression<VectorExpression<Scalar>, Scalar, 1>({static_cast<size_t>(value.size())}),// Pass the size of the vector to the base class's sizes_ array
        m_value(value)
        {
        }

        gsMatrix<Scalar> eval(const index_t k) const // override keyword for clarity
        {
            return m_value;
        }
    };
}//namespace Expr
}//namespace gismo