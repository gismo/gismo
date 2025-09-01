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

template <typename E, size_t Order, size_t Space, size_t IsConstant>
struct ExpressionTraits<GradExpression<E, Order, Space, IsConstant>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<E>::order == 0,
    //               "GradExpression requires a scalar (order 0) operand");
    // static_assert(ExpressionTraits<E>::space != Space::None,
    //               "GradExpression requires the operand to be defined in a space");

    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = Order + 1; // Gradient increases order by 1
    static constexpr size_t space = Space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv + 1; // Increment derivative order
    static constexpr bool isConstant = IsConstant;
};

// --- GradExpression using Partial Specialization (Redesigned) ---
// Primary template: Catches all unsupported combinations with a compile-time error
template <typename E, size_t Order, size_t Space, size_t IsConstant>
class GradExpression : public UnaryOperator<GradExpression<E, Order, Space, IsConstant>>
{
    static_assert(std::is_same<E, void>::value,
                  "TODO");
};

// --- Partial Specialization 1: Gradient of a constant expression ---
template <typename E, size_t Order>
class GradExpression<E, Order, Space::None, true> : public UnaryOperator<GradExpression<E, Order, Space::None, true>>
{
    using Base = UnaryOperator<GradExpression<E, Order, Space::None, true>>;
private:
    std::array<size_t,Base::order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;

public:
    GradExpression(const E& expr)
    : Base(expr)
    {
        sizes_[0] = expr.domainDim();
        for (short_t d=0; d!=ExpressionTraits<E>::order; d++)
            sizes_[d+1] = expr.sizes()[d];
    }

    const std::array<size_t, Base::order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return this->expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(this->expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<typename Base::Scalar> & helper) const
    {
        GISMO_UNUSED(helper);
    }

    void print(std::ostream & os) const override
    {
        os<<"\u2207("<<this->expr_<<")";
    }

};

// --- Partial Specialization 2: gradient of a variable object ---
template <typename E, size_t Order>
class GradExpression<E, Order, Space::None, false> : public UnaryOperator<GradExpression<E, Order, Space::None, false>>
{
    using Base = UnaryOperator<GradExpression<E, Order, Space::None, false>>;

private:
    std::array<size_t,Base::order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;

public:
    GradExpression(const E& expr)
    : Base(expr)
    {
        sizes_[0] = this->expr_.domainDim();
        for (short_t d=0; d!=ExpressionTraits<E>::order; d++)
            sizes_[d+1] = this->expr_.sizes()[d];
    }

    const std::array<size_t, Base::order> & sizes() const
    {
    return sizes_;
    }

    size_t domainDim() const
    {
        return this->expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(this->expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void print(std::ostream & os) const override
    {
        os<<"\u2207("<<this->expr_<<")";
    }

};

// Partial specialization for dot product with space=0
// ∇(A·B) = (A·∇)B + (B·∇)A + A×(∇×B) + B×(∇×A)
// For vectors A,B (order 1), result is a vector (order 1)
template <typename LhsExpr, typename RhsExpr>
auto grad(const InnerProductExpression<LhsExpr,RhsExpr,1,1,0,0>& expr)
-> decltype(dot(expr.lhs(), grad(expr.rhs())) + dot(expr.rhs(), grad(expr.lhs())) + cross(expr.lhs(), curl(expr.rhs())) + cross(expr.rhs(), curl(expr.lhs())))
{
    return (dot(expr.lhs(), grad(expr.rhs())) + dot(expr.rhs(), grad(expr.lhs())) + cross(expr.lhs(), curl(expr.rhs())) + cross(expr.rhs(), curl(expr.lhs())));
}

// Generic factory function for easy creation
template <typename E>
GradExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> grad(const E& expr)
{
    return GradExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>(const_cast<E&>(expr));
}

// Partial specialization for addition
// ∇(X+Y) = ∇X + ∇Y
template <typename LhsExpr, typename RhsExpr>
auto grad(const AddExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(grad(expr.lhs()) + grad(expr.rhs()))
{
    return grad(expr.lhs()) + grad(expr.rhs());
}

// Partial specialization for multiplication
// ∇(X*Y) = (∇X)*Y + X*(∇Y)
template <typename LhsExpr, typename RhsExpr>
auto grad(const ProductExpression<LhsExpr,RhsExpr,
                                    ExpressionTraits<LhsExpr>::order,
                                    ExpressionTraits<RhsExpr>::order,
                                    ExpressionTraits<LhsExpr>::space,
                                    ExpressionTraits<RhsExpr>::space> expr)
-> decltype(grad(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs()))
{
    return grad(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs());
}

// Partial specialization for multiplication of scalars
// ∇(f*g) = (∇f)*g + f*(∇g)
// where f,g are scalar functions (order 0)
// and ∇f,∇g are vector functions (order 1)
template <typename LhsExpr, typename RhsExpr>
auto grad(const ProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(expr.lhs() * grad(expr.rhs()) + expr.rhs() * grad(expr.lhs()))
{
    return expr.lhs() * grad(expr.rhs()) + expr.rhs() * grad(expr.lhs());
}

// Partial specialization for multiplication of scalar and vector
// ∇(f*V) = (∇f)⊗V + f*(∇V)
// where f is a scalar function (order 0)
// V is a vector function (order 1)
// ∇f is a vector function (order 1)
// ∇V is a matrix function (order 2)
// ⊗ denotes the outer product
template <typename LhsExpr, typename RhsExpr>
auto grad(const ProductExpression<LhsExpr,RhsExpr,0,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(expr.lhs() * grad(expr.rhs()) + outer(grad(expr.lhs()), expr.rhs()))
{
    return expr.lhs() * grad(expr.rhs()) + outer(grad(expr.lhs()), expr.rhs());
}

// Partial specialization for multiplication of vector and scalar
// ∇(V*f) = (∇V)⊗f + V⊗(∇f)
// where V is a vector function (order 1)
// f is a scalar function (order 0)
// ∇V is a matrix function (order 2)
// ∇f is a vector function (order 1)
// ⊗ denotes the outer product
template <typename LhsExpr, typename RhsExpr>
auto grad(const ProductExpression<LhsExpr,RhsExpr,1,0,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(outer(grad(expr.lhs()),expr.rhs()) + expr.lhs() * grad(expr.rhs()))
{
    return outer(grad(expr.lhs()), expr.rhs()) + expr.lhs() * grad(expr.rhs());
}

// Partial specialization for cross product
// ∇(A×B) = (∇A)×B - (∇B)×A
// For vectors A,B (order 1), result is a matrix (order 2)
template <typename LhsExpr, typename RhsExpr>
auto grad(const CrossProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(cross(grad(expr.lhs()), expr.rhs()) - cross(grad(expr.rhs()), expr.lhs()))
{
    return (cross(grad(expr.lhs()), expr.rhs()) - cross(grad(expr.rhs()), expr.lhs()));
}

// Partial specialization for outer product
// ∇(A⊗B) = (∇A)⊗B + A⊗(∇B)
// For vectors A,B (order 1), result is a tensor (order 3)
template <typename LhsExpr, typename RhsExpr>
auto grad(const OuterProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(outer(grad(expr.lhs()), expr.rhs()) + outer(expr.lhs(), grad(expr.rhs())))
{
    return (outer(grad(expr.lhs()), expr.rhs()) + outer(expr.lhs(), grad(expr.rhs())));
}

// Quotient rule for division by a scalar
// ∇(ψ/φ) = (φ∇ψ - ψ∇φ)/φ²
// For scalar functions ψ,φ (order 0), result is a vector (order 1)
template <typename LhsExpr, typename RhsExpr>
auto grad(const DivisionExpression<LhsExpr,RhsExpr,0,0,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype((expr.rhs() * grad(expr.lhs()) - expr.lhs() * grad(expr.rhs())) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * grad(expr.lhs()) - expr.lhs() * grad(expr.rhs())) / (expr.rhs() * expr.rhs());
}

// Quotient rule for vector divided by scalar
// ∇(A/φ) = (φ∇A - A⊗∇φ)/φ²
// For vector A (order 1) and scalar φ (order 0), result is a matrix (order 2)
template <typename LhsExpr, typename RhsExpr>
auto grad(const DivisionExpression<LhsExpr,RhsExpr,1,0,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype((expr.rhs() * grad(expr.lhs()) - outer(expr.lhs(), grad(expr.rhs()))) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * grad(expr.lhs()) - outer(expr.lhs(), grad(expr.rhs()))) / (expr.rhs() * expr.rhs());
}

// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Gradient of divergence is undefined
// ∇(∇•A) is undefined because divergence of a vector produces a scalar,
// but we cannot take gradient of a divergence operation directly
template <typename E>
auto grad(const DivExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> void
{
    GISMO_ERROR("∇(∇•): gradient of divergence is not implemented");
}

// Gradient of curl for vectors is undefined in 3D
// ∇(∇×A) is undefined because curl of a vector produces a vector,
// but taking gradient of curl is not a standard vector calculus operation
template <typename E>
auto grad(const CurlExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> void
{
    static_assert(false,"∇(∇×) is undefined: gradient of curl is not a valid vector calculus operation");
}

// Gradient of Laplacian is undefined
// ∇(∇²ψ) is undefined as a direct operation
template <typename E>
auto grad(const LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> void
{
    static_assert(false,"∇(∇²) = ∇³ is not implemented as a separate expression type.");
}

} // namespace Expr
} // namespace gismo