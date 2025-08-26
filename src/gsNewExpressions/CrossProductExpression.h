/** @file BaseObject.h

    @brief

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

// --- CrossProductExpression using Partial Specialization (Redesigned) ---

// Primary template: Catches all unsupported combinations with a compile-time error
template <typename LhsExpr, typename RhsExpr, typename Enable = void>
class CrossProductExpression
{
    static_assert(std::is_same<LhsExpr, void>::value,
                  "CrossProductExpression: Unsupported tensor order combination for addition.");
};

// --- Partial specialization 1: Vector cross product
template <typename LhsExpr, typename RhsExpr, typename Enable>
struct ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr, Enable>>
{
    using Scalar = typename ExpressionTraits<LhsExpr>::Scalar;
    static constexpr size_t order = 1;
    static constexpr size_t space = Space::None;
    // CORRECT??
    static constexpr size_t deriv = LhsExpr::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};


template <typename LhsExpr, typename RhsExpr>
class CrossProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==1 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::order==1 &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type
> : public BaseExpression<CrossProductExpression<LhsExpr, RhsExpr>>
{
public:
    typedef typename ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr>>::order;
    static constexpr size_t space = ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr>>::space;
    static constexpr size_t deriv = ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<CrossProductExpression<LhsExpr, RhsExpr>>::isConstant;

    const std::array<size_t, order> & sizes() const { return sizes_; }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return lhs_expr_.domainDim();
    }

    const LhsExpr& lhs() const {return lhs_expr_;}
    const RhsExpr& rhs() const {return rhs_expr_;}

protected:

    std::array<size_t,order> sizes_;

public:
    CrossProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs),
          sizes_({3})
    {
        GISMO_ENSURE(lhs.sizes()[0]==3,"lhs must be a vector of size 3");
        GISMO_ENSURE(rhs.sizes()[0]==3,"rhs must be a vector of size 3");
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        GISMO_ASSERT((lhs_val.rows()==3 && lhs_val.cols()==1) || (lhs_val.cols()==3 && lhs_val.rows()==1),"lhs should be a vector (ncols!=1)");
        GISMO_ASSERT((rhs_val.rows()==3 && rhs_val.cols()==1) || (rhs_val.cols()==3 && rhs_val.rows()==1),"rhs should be a vector (ncols!=1)");
        gsVector<Scalar,3> lhs_vec;
        lhs_vec<<(lhs_val.at(0),lhs_val.at(1),lhs_val.at(2));
        gsVector<Scalar,3> rhs_vec;
        rhs_vec<<(rhs_val.at(0),rhs_val.at(1),rhs_val.at(2));
        return lhs_vec.cross(rhs_vec);
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    const SpaceObject<Scalar, space, order> & rowVar() const
    {
        return lhs_expr_.rowVar(); // Use left operand's row variable
    }

    const SpaceObject<Scalar, space, order> & colVar() const
    {
        return lhs_expr_.colVar(); // Use left operand's column variable
    }

    void print(std::ostream & os) const
    {
        os<<lhs_expr_<<"\u2022"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};

// Generic dot operator to create CrossProductExpression instances
template <typename LhsExpr, typename RhsExpr>
CrossProductExpression<LhsExpr, RhsExpr>
cross(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return CrossProductExpression<LhsExpr, RhsExpr>(lhs, rhs);
}

template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order==1 && ExpressionTraits<LhsExpr>::space == Space::None &&
    ExpressionTraits<RhsExpr>::order==1 && ExpressionTraits<RhsExpr>::space == Space::None,
    CrossProductExpression<LhsExpr,RhsExpr>
>::type
cross(const TransposeExpression<LhsExpr>& lhs,
      const RhsExpr& rhs)
{
    GISMO_ERROR("Cross product is not defined for vectors with different transposition");
}

template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order==1 && ExpressionTraits<LhsExpr>::space == Space::None &&
    ExpressionTraits<RhsExpr>::order==1 && ExpressionTraits<RhsExpr>::space == Space::None,
    CrossProductExpression<LhsExpr,RhsExpr>
>::type
cross(const LhsExpr& lhs,
      const TransposeExpression<RhsExpr>& rhs)
{
    GISMO_ERROR("Cross product is not defined for vectors with different transposition");
}

template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order==1 && ExpressionTraits<LhsExpr>::space == Space::None &&
    ExpressionTraits<RhsExpr>::order==1 && ExpressionTraits<RhsExpr>::space == Space::None,
    TransposeExpression<
        CrossProductExpression<
            TransposeExpression<LhsExpr>,
            TransposeExpression<RhsExpr>
    >>
>::type
cross(const TransposeExpression<LhsExpr>& lhs,
      const TransposeExpression<RhsExpr>& rhs)
{
    return TransposeExpression<
                CrossProductExpression<
                    TransposeExpression<LhsExpr>,
                    TransposeExpression<RhsExpr>
    >>(CrossProductExpression<
                    TransposeExpression<LhsExpr>,
                    TransposeExpression<RhsExpr>
                >(lhs,rhs));
}

}//namespace Expr
}//namespace gismo