/** @file ArrayExpression.h

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
struct ExpressionTraits<ArrayExpression<E>>
{
    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order;
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv;
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;

};

template<typename E>
class ArrayExpression : public UnaryOperator< ArrayExpression<E> >
{
    using Base = UnaryOperator<ArrayExpression<E>>;


public:
    ArrayExpression(const E& expr)
    :
    Base(expr)
    {
    }

    std::array<size_t, Base::order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        return expr_.eval(k);
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        expr_.parse(helper);
    }

    void print(std::ostream & os) const override
    {
        os<<"array("<<expr_<<")";
    }

protected:
    using Base::expr_;

};

}//namespace Expr
}//namespace gismo