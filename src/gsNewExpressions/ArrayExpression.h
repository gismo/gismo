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
    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order;
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv;
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;

};

template<typename E>
class ArrayExpression : public BaseExpression< ArrayExpression<E> >
{
    using Base = BaseExpression<ArrayExpression<E>>;

public:
    typedef typename ExpressionTraits<ArrayExpression<E>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<ArrayExpression<E>>::order;
    static constexpr size_t space = ExpressionTraits<ArrayExpression<E>>::space;
    static constexpr size_t deriv = ExpressionTraits<ArrayExpression<E>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<ArrayExpression<E>>::isConstant;

    // sizes
    const std::array<size_t, order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

private:
    const E& expr_;

public:
    ArrayExpression(const E& expr)
    :
    BaseExpression<ArrayExpression<E>>(),
    expr_(expr)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        return expr_.eval(k);
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        expr_.parse(helper);
    }

    const SpaceObject<Scalar,space,order> & rowVar() const {return expr_.rowVar();}
    const SpaceObject<Scalar,space,order> & colVar() const {return expr_.colVar();}

    void print(std::ostream & os) const
    {
        os<<"array("<<expr_<<")";
    }

};

}//namespace Expr
}//namespace gismo