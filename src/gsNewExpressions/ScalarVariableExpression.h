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
    class ScalarVariableExpression : public BaseExpression<ScalarVariableExpression<T>, T, 0, _isConstant>
    {
        using Base = BaseExpression<ScalarVariableExpression<T>, T, 0, _isConstant>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr int order = Base::order;
        static constexpr bool isConstant = Base::isConstant;

    private:
        const gsFunctionSet<Scalar> & m_function;
    public:
        // Explicit constructor for a scalar value
        explicit ScalarVariableExpression(const gsFunctionSet<Scalar> & function)
        :
        BaseExpression<ScalarVariableExpression<Scalar>,Scalar, 0>({}),    // Pass empty array for sizes_ for order 0 (scalar)
        m_function(function)
        {
            GISMO_ASSERT(m_function.targetDim()==1,"The function is not a scalar function, but maps from R^"<<m_function.domainDim()<<" to R^"<<m_function.targetDim());
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