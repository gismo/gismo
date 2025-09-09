/** @file OuterProductExpression.h

    @brief Outer product expression class

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
// --- OuterProductExpression using Partial Specialization (Redesigned) ---

// --- ExpressionTraits specializations ---

// ExpressionTraits for vector outer product (1,1) -> matrix (2)
template <typename LhsExpr, typename RhsExpr>
struct ExpressionTraits<OuterProductExpression<LhsExpr, RhsExpr, 1, 1, Space::None, Space::None>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = 2;
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

// --- Partial specialization: Vector outer product
template <typename LhsExpr, typename RhsExpr>
class OuterProductExpression<LhsExpr, RhsExpr, 1, 1, Space::None, Space::None>
 : public BinaryOperator<OuterProductExpression<LhsExpr, RhsExpr, 1, 1, Space::None, Space::None>>
{
    using Base = BinaryOperator<OuterProductExpression<LhsExpr, RhsExpr, 1, 1, Space::None, Space::None>>;
    using Scalar = typename Base::Scalar;
    using Base::order;
protected:

public:
    OuterProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        this->sizes_[0] = lhs.sizes()[0];
        this->sizes_[1] = rhs.sizes()[0];
    }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return this->lhs_expr_.domainDim();
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = this->rhs_expr_.eval(k);
        GISMO_ASSERT(lhs_val.rows()==1 || lhs_val.cols()==1,"lhs should be a vector (ncols!=1)");
        GISMO_ASSERT(rhs_val.rows()==1 || rhs_val.cols()==1,"rhs should be a vector (ncols!=1)");
        lhs_val.resize(lhs_val.size(),1);
        rhs_val.resize(rhs_val.size(),1);
        return lhs_val.operator*(rhs_val.transpose());
    }

    void print(std::ostream & os) const
    {
        os<<this->lhs_expr_<<"\u2297"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// --- Function overloads

template<typename LhsExpr, typename RhsExpr>
typename std::enable_if<ExpressionTraits<LhsExpr>::order == 1 &&
                        ExpressionTraits<RhsExpr>::order == 1,
                        OuterProductExpression<LhsExpr, RhsExpr, 1, 1, Space::None, Space::None>>::type
outer(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return OuterProductExpression<LhsExpr, RhsExpr, 1, 1, Space::None, Space::None>(lhs, rhs);
}

}//namespace Expr
}//namespace gismo