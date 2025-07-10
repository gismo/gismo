/** @file BaseExpression.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gsCore/gsFunctionSet.h>

#pragma once

namespace gismo
{
namespace Expr
{

    template <class T, bool _isConstant = true>
    class VectorVariableExpression : public BaseExpression<VectorVariableExpression<T>, T, 0, _isConstant>
    {
        using Base = BaseExpression<VectorVariableExpression<T>, T, 0, _isConstant>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr int order = Base::order;
        static constexpr bool isConstant = Base::isConstant;

    private:
        const gsFunctionSet<Scalar> & m_function;
    public:
        // Explicit constructor for a scalar value
        explicit VectorVariableExpression(const gsFunctionSet<Scalar> & function)
        :
        BaseExpression<VectorVariableExpression<Scalar>,Scalar, 1>({funciton.targetDim()}),    // Pass empty array for sizes_ for order 0 (scalar)
        m_function(function)
        {
        }

        gsMatrix<Scalar> eval(const index_t k) const
        {
            // TEMPORARY
            gsMatrix<Scalar> pts(function.domainDim(),1);
            pts.setZero();
            return m_function.eval(pts);
        }

    };
}//namespace Expr
}//namespace gismo