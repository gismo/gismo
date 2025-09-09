/** @file LaplExpression.h

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

// --- LaplExpression ---
// Primary template: Catches all unsupported combinations with a compile-time error
template <typename E, size_t Order, size_t Space, size_t IsConstant>
struct ExpressionTraits<LaplExpression<E, Order, Space, IsConstant>>
{
    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order; // Laplacian decreases order by 2
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv + 2; // Increment derivative order
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

template <typename E, size_t Order, size_t Space, size_t IsConstant>
class LaplExpression : public UnaryOperator<LaplExpression<E, Order, Space, IsConstant>>
{
    static_assert(Order == 0,
                  "LaplExpression: Unsupported tensor order. Only scalars (order 0) are supported for Laplacian operations.");
    static_assert(IsConstant == 0 || IsConstant == 1,
                  "LaplExpression: IsConstant must be 0 (variable) or 1 (constant).");
};

// --- Partial Specialization 1: Laplacian of a constant expression ---
template <typename E, size_t Order, size_t Space>
class LaplExpression<E, Order, Space, 1> : public UnaryOperator<LaplExpression<E, Order, Space, 1>>
{
    using Base = UnaryOperator<LaplExpression<E, Order, Space, 1>>;
    typedef typename Base::Scalar T;
    using Base::order;
private:
    mutable gsMatrix<T> tmp;
    using Base::expr_;
public:
    LaplExpression(const E& expr)
    :
    Base(expr)
    {
    }

    const std::array<size_t, order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    gsMatrix<T> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        GISMO_UNUSED(helper);
    }

    void print(std::ostream & os) const
    {
        os<<"\u0394("<<expr_<<")";
    }
};

// --- Partial Specialization 2: Laplacian of a variable object ---
template <typename E, size_t Order, size_t Space>
class LaplExpression<E, Order, Space, 0> : public UnaryOperator<LaplExpression<E, Order, Space, 0>>
{
    using Base = UnaryOperator<LaplExpression<E, Order, Space, 0>>;
    typedef typename Base::Scalar T;
    using Base::order;

private:
    mutable gsMatrix<T> tmp;
    using Base::expr_;

public:
    LaplExpression(const E & expr)
    :
    Base(expr)
    {
    }

    const std::array<size_t, Base::order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    gsMatrix<T> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void print(std::ostream & os) const
    {
        os<<"\u0394("<<expr_<<")";
    }

};

// Generic factory function for easy creation
template <typename E>
LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> lapl(const E& expr)
{
    return LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>(expr);
}

// Partial specialization for addition
template <typename LhsExpr, typename RhsExpr>
auto lapl(const AddExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(lapl(expr.lhs()) + lapl(expr.rhs()))
{
    return lapl(expr.lhs()) + lapl(expr.rhs()); // ∇²(ψ + φ) = ∇²ψ + ∇²φ
}

// Partial specialization for subtraction
// ∇²(ψ - φ) = ∇²ψ - ∇²φ
template <typename LhsExpr, typename RhsExpr>
auto lapl(const SubtractExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(lapl(expr.lhs()) - lapl(expr.rhs()))
{
    return lapl(expr.lhs()) - lapl(expr.rhs());
}

// Partial specialization for product (second derivative identity)
// ∇²(ψφ) = ψ∇²φ + 2∇ψ•∇φ + φ∇²ψ
// For scalar functions ψ,φ (order 0), result is scalar (order 0)
template <typename LhsExpr, typename RhsExpr>
auto lapl(const ProductExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(expr.lhs() * lapl(expr.rhs()) + expr.rhs() * lapl(expr.lhs()) + ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs())))
{
    static_assert(ExpressionTraits<LhsExpr>::order == 0 && ExpressionTraits<RhsExpr>::order == 0,
                  "Laplacian of product only implemented for scalar functions (order 0)");
    return expr.lhs() * lapl(expr.rhs()) + expr.rhs() * lapl(expr.lhs()) + ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs()));
}

// Partial specialization for division of scalar by scalar
// ∇²(ψ/φ) = (φ∇²ψ - ψ∇²φ - 2∇ψ·∇φ)/φ²
// For scalar functions ψ,φ (order 0), result is scalar (order 0)
template <typename LhsExpr, typename RhsExpr>
auto lapl(const DivisionExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype((expr.rhs() * lapl(expr.lhs()) - expr.lhs() * lapl(expr.rhs()) - ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs()))) / (expr.rhs() * expr.rhs()))
{
    static_assert(ExpressionTraits<LhsExpr>::order == 0 && ExpressionTraits<RhsExpr>::order == 0,
                  "Laplacian of division only implemented for scalar functions (order 0)");
    return (expr.rhs() * lapl(expr.lhs()) - expr.lhs() * lapl(expr.rhs()) - ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs()))) / (expr.rhs() * expr.rhs());
}

// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Laplacian of gradient produces third derivatives
// ∇²(∇ψ) = ∇(∇²ψ) is a valid operation (third derivative)
// This could be implemented as a third derivative expression if needed
template <typename E>
auto lapl(const GradExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> decltype(grad(lapl(expr.expr())))
{
    return grad(lapl(expr.expr()));
}

// Laplacian of divergence is undefined for vector fields
// ∇²(∇•A) is undefined because divergence of a vector produces a scalar,
// and we already have ∇²(scalar), so this would be valid if the input expression is appropriate
template <typename E>
auto lapl(const DivExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> void
{
    static_assert(false, "∇²(∇•A) is undefined for vector fields.");
}

// Laplacian of curl is defined: ∇²(∇×A) = ∇×(∇²A)
// This is a valid vector calculus identity
template <typename E>
auto lapl(const CurlExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> decltype(curl(lapl(expr.expr())))
{
    return curl(lapl(expr.expr()));
}

// Laplacian of Laplacian (fourth derivative)
// ∇²(∇²ψ) = ∇⁴ψ is a valid operation (biharmonic operator)
// However, it is not implemented as a separate expression type
template <typename E>
auto lapl(const LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> void
{
    static_assert(false,"∇²(∇²ψ) = ∇⁴ψ is not implemented as a separate expression type.");
    // return lapl(lapl(expr.expr()));
}

}//namespace Expr
}//namespace gismo