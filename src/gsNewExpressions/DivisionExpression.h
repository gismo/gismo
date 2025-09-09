/** @file DivisionExpression.h

    @brief Division expression class

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

// --- DivisionExpression using Partial Specialization (Redesigned) ---

// --- Partial Specialization: Division of an expression by another expression of order 0 (scalar)
template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<DivisionExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
    //               "DivisionExpression requires same space operands");
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = LhsOrder;
    static constexpr size_t space = ExpressionTraits<LhsExpr>::space;
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr, size_t LhsOrder, size_t LhsSpace, size_t RhsSpace>
class DivisionExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>
 : public BinaryOperator<DivisionExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<DivisionExpression<LhsExpr, RhsExpr, LhsOrder, 0, LhsSpace, RhsSpace>>;

public:
    DivisionExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        // Division: result has same size as left operand (we're dividing by scalar)
        for (size_t d = 0; d < Base::order; ++d) {
            this->sizes_[d] = lhs.sizes()[d];
        }
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        gsMatrix<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        lhs_val.array() /= rhs_val.value(); // Element-wise division
        return lhs_val; // Return the modified lhs_val
    }

    size_t domainDim() const
    {
        return this->lhs_expr_.domainDim();
    }

    void print(std::ostream & os) const
    {
        os<<this->lhs_expr_<<"/"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// Generic operator/ to create DivisionExpression instances using SFINAE
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<ExpressionTraits<RhsExpr>::order == 0,
                        DivisionExpression<LhsExpr, RhsExpr, ExpressionTraits<LhsExpr>::order, 0, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>>::type
operator/(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return DivisionExpression<LhsExpr, RhsExpr, ExpressionTraits<LhsExpr>::order, 0, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

}//namespace Expr
}//namespace gismo