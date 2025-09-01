/** @file TransposeExpression.h

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

template <typename E>
struct ExpressionTraits<TransposeExpression<E>>
{
    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order;
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv;
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

template<typename E>
class TransposeExpression : public UnaryOperator< TransposeExpression<E> >
{
    static_assert(ExpressionTraits<E>::order == 1 || ExpressionTraits<E>::order == 2,
                    "TransposeExpression: Only vector (order 1) or matrix (order 2) expressions can be transposed.");

    using Base = UnaryOperator<TransposeExpression<E>>;

public:
    TransposeExpression(const E& expr)
    :
    Base(expr)
    {
    }

    const std::array<size_t, Base::order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        return expr_.eval(k).transpose();
    }

    void parse(ExpressionHelper<typename Base::Scalar> & helper) const override
    {
        expr_.parse(helper);
    }

    void print(std::ostream & os) const override
    {
        os<<"("<<expr_<<")\u1D40";
    }

protected:
    using Base::expr_;

};

// Factory function for easier creation
template <typename E>
TransposeExpression<E> transpose(const E& expr)
{
    return TransposeExpression<E>(expr);
}

}//namespace Expr
}//namespace gismo