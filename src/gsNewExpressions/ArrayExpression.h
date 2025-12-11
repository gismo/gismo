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
    static constexpr size_t Order = ExpressionTraits<E>::Order;
    static constexpr size_t Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;

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

    std::array<size_t, Base::Order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    ExpressionValue<Scalar> eval(const index_t k) const
    {
        return expr_.eval(k);
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        expr_.parse(helper);
    }

    void print(std::ostream & os) const
    {
        os<<"array("<<expr_<<")";
    }

protected:
    using Base::expr_;

};

}//namespace Expr
}//namespace gismo