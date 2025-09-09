/** @file InnerProductExpression.h

    @brief Inner product expression class

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
// --- InnerProductExpression using Partial Specialization (Redesigned) ---

// --- ExpressionTraits specializations ---

// ExpressionTraits for vector inner product (1,1) -> scalar (0)
template <typename LhsExpr, typename RhsExpr>
struct ExpressionTraits<InnerProductExpression<LhsExpr, RhsExpr, 1, 1, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = 0;
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = 0;//TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

// ExpressionTraits for matrix inner product (2,2) -> scalar (0)
template <typename LhsExpr, typename RhsExpr>
struct ExpressionTraits<InnerProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>>
{
    typedef LhsExpr LhsType;
    typedef RhsExpr RhsType;

    typedef typename ExpressionTraits<LhsExpr>::Scalar Scalar;
    static constexpr size_t order = 0;
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = 0;//TODO
    static constexpr bool isConstant = ExpressionTraits<LhsExpr>::isConstant && ExpressionTraits<RhsExpr>::isConstant;
};

// --- Partial Specialization: Inner product of vectors (dot product) ---
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
class InnerProductExpression<LhsExpr, RhsExpr, 1, 1, LhsSpace, RhsSpace>
 : public BinaryOperator<InnerProductExpression<LhsExpr, RhsExpr, 1, 1, LhsSpace, RhsSpace>>
{
    using Base = BinaryOperator<InnerProductExpression<LhsExpr, RhsExpr, 1, 1, LhsSpace, RhsSpace>>;

public:
    InnerProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        // Inner product of vectors results in a scalar (order 0 - no dimensions)
        // sizes_ is an empty array for order 0, so nothing to initialize
    }

    size_t domainDim() const
    {
        return this->lhs_expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        // For vectors (order 1), inner product is scalar (order 0)
        gsMatrix<typename Base::Scalar> lhs_eval = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_eval = this->rhs_expr_.eval(k);

        // Compute inner product
        typename Base::Scalar result = 0;
        for (index_t i = 0; i < lhs_eval.rows(); ++i)
        {
            result += lhs_eval(i, 0) * rhs_eval(i, 0);
        }

        gsMatrix<typename Base::Scalar> res(1, 1);
        res(0, 0) = result;
        return res;
    }

    void print(std::ostream & os) const
    {
        os<<this->lhs_expr_<<"\u2027"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// --- Partial specialization 2: Matrix inner product ---
template <typename LhsExpr, typename RhsExpr>
class InnerProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
 : public BinaryOperator<InnerProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>>
{
    using Base = BinaryOperator<InnerProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>>;

protected:
public:
    InnerProductExpression(const LhsExpr& lhs, const RhsExpr& rhs)
        : Base(lhs, rhs)
    {
        // Inner product results in a scalar (order 0), so sizes_ is empty
    }

    size_t domainDim() const
    {
        gsWarn<<"Correct?\n";
        return this->lhs_expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        gsMatrix<typename Base::Scalar> lhs_val = this->lhs_expr_.eval(k);
        gsMatrix<typename Base::Scalar> rhs_val = this->rhs_expr_.eval(k);
        return (lhs_val.array()*rhs_val.array()).sum().matrix();
    }

    void print(std::ostream & os) const
    {
        os<<this->lhs_expr_<<":"<<this->rhs_expr_;
    }

private:
    using Base::lhs_expr_;
    using Base::rhs_expr_;
};

// Generic dot operator to create InnerProductExpression instances using SFINAE
// Vector-vector inner product
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 1 &&
    ExpressionTraits<RhsExpr>::order == 1,
    InnerProductExpression<LhsExpr, RhsExpr, 1, 1, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
dot(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return InnerProductExpression<LhsExpr, RhsExpr, 1, 1, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

// `inner` alias for `dot`
template <typename LhsExpr, typename RhsExpr>
auto inner(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
-> decltype(dot(lhs, rhs))  // Use decltype to deduce return type
{
    return dot(lhs, rhs);    // Call the existing dot function
}

// Matrix-matrix inner product
template <typename LhsExpr, typename RhsExpr>
typename std::enable_if<
    ExpressionTraits<LhsExpr>::order == 2 &&
    ExpressionTraits<RhsExpr>::order == 2,
    InnerProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>
>::type
ddot(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
{
    return InnerProductExpression<LhsExpr, RhsExpr, 2, 2, ExpressionTraits<LhsExpr>::space, ExpressionTraits<RhsExpr>::space>(lhs, rhs);
}

// `inner` alias for `ddot`
template <typename LhsExpr, typename RhsExpr>
auto inner(const BaseExpression<LhsExpr>& lhs, const BaseExpression<RhsExpr>& rhs)
-> decltype(ddot(lhs, rhs))  // Use decltype to deduce return type
{
    return ddot(lhs, rhs);    // Call the existing ddot function
}

// Specialized dot product for vector · gradient (directional derivative)
// This computes (A·∇)B where A is a vector and ∇B is the gradient matrix
// Mathematically: A^T * ∇B = ∇B^T*A, but we need the result as a vector (not transpose)
template <typename LhsExpr, typename RhsExpr>
auto dot(const BaseExpression<LhsExpr>& lhs, const GradExpression<RhsExpr, ExpressionTraits<RhsExpr>::order, ExpressionTraits<RhsExpr>::space, ExpressionTraits<RhsExpr>::isConstant>& rhs)
-> decltype(transpose(rhs) * lhs)
{
    return transpose(rhs) * lhs;
}

}//namespace Expr
}//namespace gismo