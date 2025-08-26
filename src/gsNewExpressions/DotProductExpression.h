/** @file BaseObject.h

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
// --- DotProductExpression using Partial Specialization (Redesigned) ---

// Primary template: Catches all unsupported combinations with a compile-time error
template <typename LhsExpr, typename RhsExpr, typename Enable = void>
class DotProductExpression
{
    static_assert(std::is_same<LhsExpr, void>::value,
                  "DotProductExpression: Unsupported tensor order combination for addition.");
};

// --- Partial specialization 1: Vector contraction


template <typename LhsExpr, typename RhsExpr, typename Enable>
struct ExpressionTraits<DotProductExpression<LhsExpr, RhsExpr, Enable>>
{
    using Scalar = typename ExpressionTraits<LhsExpr>::Scalar;
    static constexpr size_t order = 0;
    static constexpr size_t space = Space::None;
    // CORRECT??
    static constexpr size_t deriv = LhsExpr::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

template <typename LhsExpr, typename RhsExpr>
class DotProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==1 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::order==1 &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type
> : public BaseExpression<DotProductExpression<LhsExpr, RhsExpr>>
{
public:
    typedef typename ExpressionTraits<DotProductExpression<LhsExpr, RhsExpr>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<DotProductExpression<LhsExpr, RhsExpr>>::order;
    static constexpr size_t space = ExpressionTraits<DotProductExpression<LhsExpr, RhsExpr>>::space;
    static constexpr size_t deriv = ExpressionTraits<DotProductExpression<LhsExpr, RhsExpr>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<DotProductExpression<LhsExpr, RhsExpr>>::isConstant;

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
    DotProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> res(1,1);
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        GISMO_ASSERT(lhs_val.rows()==1 || lhs_val.cols()==1,"lhs should be a vector (ncols!=1)");
        GISMO_ASSERT(rhs_val.rows()==1 || rhs_val.cols()==1,"rhs should be a vector (ncols!=1)");
        lhs_val.resize(lhs_val.size(),1);
        rhs_val.resize(rhs_val.size(),1);
        res(0,0) = lhs_val.col(0).dot(rhs_val.col(0));
        return res;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    // const SpaceObject<Scalar, LhsExpr::space, LhsExpr::order> & rowVar() const
    // {
    //     return lhs_expr_.rowVar(); // Use left operand's row variable
    // }

    // const SpaceObject<Scalar, RhsExpr::space, RhsExpr::order> & colVar() const
    // {
    //     return rhs_expr_.colVar(); // Use right operand's column variable
    // }

    void print(std::ostream & os) const
    {
        os<<lhs_expr_<<"\u2022"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};

// --- Partial specialization 2: Matrix contraction
template <typename LhsExpr, typename RhsExpr>
class DotProductExpression<LhsExpr, RhsExpr,
    typename std::enable_if<(ExpressionTraits<LhsExpr>::order==2 &&
                             ExpressionTraits<LhsExpr>::space==Space::None &&
                             ExpressionTraits<RhsExpr>::order==2 &&
                             ExpressionTraits<RhsExpr>::space==Space::None)>::type
> : public BaseExpression<DotProductExpression<LhsExpr, RhsExpr>>
{
public:
// Scalar and Order are from ExpressionTraits
    typedef typename ExpressionTraits<RhsExpr>::Scalar Scalar;
    static constexpr size_t order = 0;
    static constexpr size_t space = ExpressionTraits<RhsExpr>::space;
    // TODO:
    // static constexpr size_t deriv = ExpressionTraits<RhsExpr>::deriv;
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;

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
    DotProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : lhs_expr_(lhs),
          rhs_expr_(rhs)
    {
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
        gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
        return (lhs_val.array()*rhs_val.array()).sum().matrix();
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        lhs_expr_.parse(helper);
        rhs_expr_.parse(helper);
    }

    const SpaceObject<Scalar, LhsExpr::space, LhsExpr::order> & rowVar() const
    {
        return lhs_expr_.rowVar(); // Use left operand's row variable
    }

    const SpaceObject<Scalar, RhsExpr::space, RhsExpr::order> & colVar() const
    {
        return rhs_expr_.colVar(); // Use right operand's column variable
    }

    void print(std::ostream & os) const
    {
        os<<lhs_expr_<<":"<<rhs_expr_;
    }

private:
    const LhsExpr& lhs_expr_;
    const RhsExpr& rhs_expr_;
};

// Generic dot operator to create DotProductExpression instances
template <typename LhsExpr, typename RhsExpr>
DotProductExpression<LhsExpr, RhsExpr>
dot(const LhsExpr& lhs, const RhsExpr& rhs)
{
    return DotProductExpression<LhsExpr, RhsExpr>(lhs, rhs);
}

}//namespace Expr
}//namespace gismo