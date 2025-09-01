/** @file CrossProductExpression.h

    @brief Cross product expression class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

// Forward declarations
namespace gismo { namespace Expr {
template <typename E> class TransposeExpression;
} }

namespace gismo
{
namespace Expr
{

// --- ExpressionTraits specializations ---
// ExpressionTraits for vector cross product (1,1) -> vector (1)
template <typename LhsExpr, typename RhsExpr>
struct ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr, 1, 1, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = 1;
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv; //TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

// --- Partial Specialization: Cross product of 3D vectors ---
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
class CrossProductExpression<LhsExpr, RhsExpr, 1, 1, LhsSpace, RhsSpace>
 : public BinaryOperator<CrossProductExpression<LhsExpr, RhsExpr, 1, 1, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<CrossProductExpression<LhsExpr, RhsExpr, 1, 1, LhsSpace, RhsSpace>>;

protected:

    std::array<size_t,Base::order> sizes_;

public:
    CrossProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs),
          sizes_({3})
    {
        GISMO_ENSURE(lhs.sizes()[0]==3,"lhs must be a vector of size 3");
        GISMO_ENSURE(rhs.sizes()[0]==3,"rhs must be a vector of size 3");
    }

    const std::array<size_t, Base::order> & sizes() const { return sizes_; }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return this->lhs_expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        gsMatrix<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        GISMO_ASSERT((lhs_val.rows()==3 && lhs_val.cols()==1) || (lhs_val.cols()==3 && lhs_val.rows()==1),"lhs should be a vector (ncols!=1)");
        GISMO_ASSERT((rhs_val.rows()==3 && rhs_val.cols()==1) || (rhs_val.cols()==3 && rhs_val.rows()==1),"rhs should be a vector (ncols!=1)");
        gsVector<typename Base::Scalar,3> lhs_vec = gsVector<typename Base::Scalar,3>::vec(lhs_val.at(0),lhs_val.at(1),lhs_val.at(2));
        gsVector<typename Base::Scalar,3> rhs_vec = gsVector<typename Base::Scalar,3>::vec(rhs_val.at(0),rhs_val.at(1),rhs_val.at(2));
        return lhs_vec.cross(rhs_vec);
    }

    void print(std::ostream & os) const override
    {
        os<<this->lhs_expr_<<"\u00D7"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// Generic cross operator to create CrossProductExpression instances using SFINAE
template <typename LeftExpr, typename RightExpr>
typename std::enable_if<
    ExpressionTraits<LeftExpr>::order == 1 && ExpressionTraits<LeftExpr>::space == Space::None &&
    ExpressionTraits<RightExpr>::order == 1 && ExpressionTraits<RightExpr>::space == Space::None,
    CrossProductExpression<LeftExpr, RightExpr, 1, 1, ExpressionTraits<LeftExpr>::space, ExpressionTraits<RightExpr>::space>
>::type
cross(const LeftExpr& lhs, const RightExpr& rhs)
{
    return CrossProductExpression<LeftExpr, RightExpr, 1, 1, ExpressionTraits<LeftExpr>::space, ExpressionTraits<RightExpr>::space>(lhs, rhs);
}

// Specializations for transpose expressions (these should trigger compile errors)
template <typename LeftExpr, typename RightExpr>
void cross(const TransposeExpression<LeftExpr>& lhs, const RightExpr& rhs)
{
    GISMO_ERROR("Cross product is not defined for vectors with different transposition");
}

template <typename LeftExpr, typename RightExpr>
void cross(const LeftExpr& lhs, const TransposeExpression<RightExpr>& rhs)
{
    GISMO_ERROR("Cross product is not defined for vectors with different transposition");
}

template <typename LeftExpr, typename RightExpr>
auto cross(const TransposeExpression<LeftExpr>& lhs, const TransposeExpression<RightExpr>& rhs)
    -> decltype(cross(lhs.expr(), rhs.expr()).transpose())
{
    return cross(lhs.expr(), rhs.expr()).transpose();
}

}//namespace Expr
}//namespace gismo