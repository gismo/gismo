/** @file BinaryOperator.h

    @brief Base class for binary operators

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

// ExpressionTraits for BinaryOperator - forwards template parameters
template <typename Operator>
struct ExpressionTraits<BinaryOperator<Operator>>
{
    typedef typename ExpressionTraits<Operator>::LhsType LhsType;
    typedef typename ExpressionTraits<Operator>::RhsType RhsType;

    typedef typename ExpressionTraits<Operator>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<Operator>::order;
    static constexpr size_t space = ExpressionTraits<Operator>::space;
    static constexpr size_t deriv = ExpressionTraits<Operator>::deriv;
    static constexpr bool isConstant = ExpressionTraits<Operator>::isConstant;
};

/**
 * \brief Base class for Binary Operators - provides common functionality
 * \tparam Operator The operator expression type (e.g., AddExpression)
 *
 */
template <typename Operator>
class BinaryOperator : public BaseExpression<BinaryOperator<Operator>>
{
public:
    typedef typename ExpressionTraits<Operator>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<Operator>::order;
    static constexpr size_t space = ExpressionTraits<Operator>::space;
    static constexpr size_t deriv = ExpressionTraits<Operator>::deriv;
    static constexpr bool isConstant = ExpressionTraits<Operator>::isConstant;

    typedef typename ExpressionTraits<Operator>::LhsType LhsType;
    typedef typename ExpressionTraits<Operator>::RhsType RhsType;

    const LhsType& lhs() const {return lhs_expr_;}
    const RhsType& rhs() const {return rhs_expr_;}

    // const SpaceObject<Scalar, LhsExpr::space, LhsExpr::order> & rowVar() const
    // {
    //     return lhs_expr_.rowVar(); // Use left operand's row variable
    // }

    // const SpaceObject<Scalar, RhsExpr::space, RhsExpr::order> & colVar() const
    // {
    //     return rhs_expr_.colVar(); // Use right operand's column variable
    // }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const override
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

protected:
    BinaryOperator(const LhsType& lhs, const RhsType& rhs)
    : lhs_expr_(lhs), rhs_expr_(rhs)
    {
    }

    const LhsType lhs_expr_;
    const RhsType rhs_expr_;
};

}//namespace Expr
}//namespace gismo