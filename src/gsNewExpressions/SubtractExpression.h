/** @file SubtractExpression.h

    @brief Subtraction expression class

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

// --- SubtractExpression using Partial Specialization (Redesigned) ---

// --- ExpressionTraits specializations ---

// ExpressionTraits for same-order subtraction (N,N) -> N
template <typename LhsExpr, typename RhsExpr, size_t Order, size_t LhsSpace, size_t RhsSpace>
struct ExpressionTraits<SubtractExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>>
{
    // Static assertions to ensure compatibility
    static_assert(ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order,
                  "SubtractExpression requires same order operands");
    static_assert(ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
                  "SubtractExpression requires same space operands");

    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = Order;
    static constexpr size_t space = ExpressionTraits<LhsExpr>::space;
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

// --- Partial Specialization: Subtraction of two expressions of the SAME ORDER (X - X) ---
template <typename LhsExpr, typename RhsExpr, size_t Order, size_t LhsSpace, size_t RhsSpace>
class SubtractExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>
 : public BinaryOperator<SubtractExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<SubtractExpression<LhsExpr, RhsExpr, Order, Order, LhsSpace, RhsSpace>>;
public:
    SubtractExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        for (size_t d=0; d!=Base::order; ++d) {
            GISMO_ENSURE(this->lhs_expr_.sizes()[d] == this->rhs_expr_.sizes()[d],"SubtractExpression requires same sizes in each dimension");
            this->sizes_[d] = this->lhs_expr_.sizes()[d];  // Fill Base::sizes_
        }
    }

    size_t domainDim() const
    {
        return this->lhs_expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        gsMatrix<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        lhs_val.array() -= rhs_val.array(); // Element-wise subtraction
        return lhs_val; // Return the modified lhs_val
    }

    void print(std::ostream & os) const
    {
        os << this->lhs_expr_ << "-" << this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// Generic operator- to create SubtractExpression instances using SFINAE
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == ExpressionTraits<RhsExpr>::order &&
    ExpressionTraits<LhsExpr>::space == ExpressionTraits<RhsExpr>::space,
    SubtractExpression<LhsExpr, RhsExpr, ExpressionTraits<LhsExpr>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
operator-(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return SubtractExpression<LhsExpr, RhsExpr, ExpressionTraits<LhsExpr>::order, ExpressionTraits<RhsExpr>::order, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

}//namespace Expr
}//namespace gismo