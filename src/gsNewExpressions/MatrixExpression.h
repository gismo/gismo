/** @file BaseExpression.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{
namespace Expr
{
template <class T, bool _isConstant = true>
class MatrixExpression : public BaseExpression<MatrixExpression<T>, T, 2, _isConstant>
{
    using Base = BaseExpression<MatrixExpression<T>, T, 2, _isConstant>;
public:
    typedef T Scalar; // Define Scalar type for this expression
    static constexpr int order = Base::order;
    static constexpr bool isConstant = Base::isConstant;
private:
    gsMatrix<Scalar> m_value;
public:

    // Constructor for a matrix filled with the same value
    MatrixExpression(size_t rows, size_t cols, T value = 0.0)
    :
    BaseExpression<MatrixExpression<Scalar>, Scalar, 2>({rows,cols}),
    m_value(gsMatrix<T>::Constant(rows,cols,value))
    {
    }

    // Explicit constructor for a gsMatrix
    explicit MatrixExpression(const gsMatrix<Scalar> & value)
    :
    BaseExpression<MatrixExpression<Scalar>, Scalar, 2>({static_cast<size_t>(value.rows()),static_cast<size_t>(value.cols())}),// Pass the size of the vector to the base class's sizes_ array
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