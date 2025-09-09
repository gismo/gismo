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
    using Base = BaseExpression<BinaryOperator<Operator>>;
protected:
    typedef typename Base::Scalar T;
    using Base::order;
public:
    typedef typename ExpressionTraits<Operator>::LhsType LhsType;
    typedef typename ExpressionTraits<Operator>::RhsType RhsType;

    const LhsType& lhs() const {return lhs_expr_;}
    const RhsType& rhs() const {return rhs_expr_;}

    // const SpaceObject<T, LhsExpr::space, LhsExpr::order> & rowVar() const
    // {
    //     return lhs_expr_.rowVar(); // Use left operand's row variable
    // }

    // const SpaceObject<T, RhsExpr::space, RhsExpr::order> & colVar() const
    // {
    //     return rhs_expr_.colVar(); // Use right operand's column variable
    // }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    void print(std::ostream & os) const
    {
        static_cast<const Operator&>(*this).print(os);
    }

    const std::array<size_t, order>& sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return static_cast<const Operator&>(*this).domainDim();
    }

    // Copy constructor
    BinaryOperator(const BinaryOperator& other)
        : BaseExpression<BinaryOperator<Operator>>(), lhs_expr_(other.lhs_expr_), rhs_expr_(other.rhs_expr_), sizes_(other.sizes_) {
    }

    // Move constructor
    BinaryOperator(BinaryOperator&& other) noexcept
    : BaseExpression<BinaryOperator<Operator>>(), lhs_expr_(give(other.lhs_expr_)), rhs_expr_(give(other.rhs_expr_)), sizes_(other.sizes_)
    {
    }

    // Copy assignment operator
    BinaryOperator& operator=(const BinaryOperator& other)
    {
        if (this != &other) {
            const_cast<LhsType&>(lhs_expr_) = other.lhs_expr_;
            const_cast<RhsType&>(rhs_expr_) = other.rhs_expr_;
            sizes_ = other.sizes_;
        }
        return *this;
    }

    // Move assignment operator
    BinaryOperator& operator=(BinaryOperator&& other) noexcept
    {
        if (this != &other) {
            const_cast<LhsType&>(lhs_expr_) = give(other.lhs_expr_);
            const_cast<RhsType&>(rhs_expr_) = give(other.rhs_expr_);
            sizes_ = other.sizes_;
        }
        return *this;
    }

protected:
    BinaryOperator(const LhsType& lhs, const RhsType& rhs)
    : lhs_expr_(lhs), rhs_expr_(rhs), sizes_{}
    {
    }

    const LhsType lhs_expr_;
    const RhsType rhs_expr_;

    // Common sizes array for all binary operators
    std::array<size_t, order> sizes_;
};

}//namespace Expr
}//namespace gismo