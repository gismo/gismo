/** @file GradExpression.h

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

template <typename E, typename Enable>
struct ExpressionTraits<GradExpression<E, Enable>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<E>::order == 0,
    //               "GradExpression requires a scalar (order 0) operand");
    // static_assert(ExpressionTraits<E>::space != Space::None,
    //               "GradExpression requires the operand to be defined in a space");

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order+1; // Gradient results in a vector (order 1)
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv + 1; // Increment derivative order
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

// --- GradExpression using Partial Specialization (Redesigned) ---
// Primary template: Catches all unsupported combinations with a compile-time error
template <typename E, typename Enable = void>
class GradExpression : public BaseExpression<GradExpression<E>>
{
    static_assert(std::is_same<Enable, void>::value,
                  "TODO");
};

// --- Partial Specialization 1: Gradient of a constant expression ---
template <typename E>
class GradExpression<E,
    typename std::enable_if<(ExpressionTraits<E>::isConstant)>::type
> : public BaseExpression<GradExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant)>::type>>
{
    using Base = BaseExpression<GradExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant)>::type>>;

public:
    typedef typename ExpressionTraits<GradExpression<E>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<GradExpression<E>>::order;
    static constexpr size_t space = ExpressionTraits<GradExpression<E>>::space;
    static constexpr size_t deriv = ExpressionTraits<GradExpression<E>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<GradExpression<E>>::isConstant;

    const std::array<size_t, order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

private:
    std::array<size_t,order> sizes_;
    mutable gsMatrix<Scalar> tmp;
    const E& expr_;

public:
    GradExpression(const E& expr)
    :
    BaseExpression<GradExpression<E>>(),
    expr_(expr)
    {
        sizes_[0] = expr_.domainDim();
        for (short_t d=0; d!=E::order; d++)
            sizes_[d+1] = expr_.sizes()[d];
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.source().domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
    }

};

// --- Partial Specialization 2: gradient of a variable object ---
template <typename E>
class GradExpression<E,
    typename std::enable_if<(ExpressionTraits<E>::isConstant==false)>::type
> : public BaseExpression<GradExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant==false)>::type>>
{
    using Base = BaseExpression<GradExpression<E, typename std::enable_if<(ExpressionTraits<E>::isConstant==false)>::type>>;
public:
    typedef typename ExpressionTraits<GradExpression<E>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<GradExpression<E>>::order;
    static constexpr size_t space = ExpressionTraits<GradExpression<E>>::space;
    static constexpr size_t deriv = ExpressionTraits<GradExpression<E>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<GradExpression<E>>::isConstant;

    const std::array<size_t, order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

private:
    std::array<size_t,order> sizes_;
    mutable gsMatrix<Scalar> tmp;
    E & expr_;

public:
    GradExpression(E & expr)
    :
    BaseExpression<GradExpression<E>>(),
    expr_(expr)
    {
        sizes_[0] = expr_.domainDim();
        for (short_t d=0; d!=E::order; d++)
            sizes_[d+1] = expr_.sizes()[d];
    }

    gsMatrix<Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<Scalar> & helper) const
    {
        expr_.parse(helper);
        expr_.setDerivative(deriv);
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207("<<expr_<<")";
    }

};

// Generic factory function for easy creation
template <typename E>
GradExpression<E> grad(E& expr)
{
    return GradExpression<E>(expr);
}

// Partial specialization for addition
template <typename LhsExpr, typename RhsExpr>
auto grad(const AddExpression<LhsExpr,RhsExpr> expr)
-> AddExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>
{
    return AddExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>(grad(expr.lhs()),grad(expr.rhs()));
}

// // Partial specialization for multiplication
// template <typename LhsExpr, typename RhsExpr>
// auto grad(const ProductExpression<LhsExpr,RhsExpr> expr)
// -> ProductExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>
// {
//     return ProductExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>(grad(expr.lhs()),grad(expr.rhs()));
// }

// // Partial specialization for dot product
// template <typename LhsExpr, typename RhsExpr>
// auto grad(const DotProductExpression<LhsExpr,RhsExpr> expr)
// -> DotProductExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>
// {
//     return DotProductExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>(grad(expr.lhs()),grad(expr.rhs()));
// }

// // Partial specialization for cross product
// template <typename LhsExpr, typename RhsExpr>
// auto grad(const CrossProductExpression<LhsExpr,RhsExpr> expr)
// -> CrossProductExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>
// {
//     return CrossProductExpression<GradExpression<LhsExpr>, GradExpression<RhsExpr>>(grad(expr.lhs()),grad(expr.rhs()));
// }

}//namespace Expr
}//namespace gismo