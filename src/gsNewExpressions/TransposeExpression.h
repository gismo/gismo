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
    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order;
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv;
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;

};

template<typename E>
class TransposeExpression : public BaseExpression< TransposeExpression<E> >
{
    static_assert(ExpressionTraits<E>::order == 1 || ExpressionTraits<E>::order == 2,
                    "TransposeExpression: Only vector (order 1) or matrix (order 2) expressions can be transposed.");

    using Base = BaseExpression<TransposeExpression<E>>;

public:
    typedef typename ExpressionTraits<TransposeExpression<E>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<TransposeExpression<E>>::order;
    static constexpr size_t space = ExpressionTraits<TransposeExpression<E>>::space;
    static constexpr size_t deriv = ExpressionTraits<TransposeExpression<E>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<TransposeExpression<E>>::isConstant;

    // sizes
    const std::array<size_t, order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

public:
    // typedef typename E::Scalar Scalar; // Define Scalar type for this expression
    // static constexpr size_t order = Base::order;
    // static constexpr bool isConstant = Base::isConstant;
    // static constexpr size_t space = Base::space;

private:
    const E& expr_;

public:
    TransposeExpression(const E& expr)
    :
    BaseExpression<TransposeExpression<E>>(),
    expr_(expr)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        return expr_.eval(k).transpose();
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        expr_.parse(helper);
    }

    const SpaceObject<Scalar,space,order> & rowVar() const {return expr_.rowVar();}
    const SpaceObject<Scalar,space,order> & colVar() const {return expr_.colVar();}

    void print(std::ostream & os) const
    {
        os<<"("<<expr_<<")\u1D40";
    }

};

// Factory function for easier creation
template <typename E>
TransposeExpression<E> transpose(const E& expr)
{
    return TransposeExpression<E>(expr);
}

}//namespace Expr
}//namespace gismo