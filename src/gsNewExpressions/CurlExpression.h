/** @file CurlExpression.h

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
struct ExpressionTraits<CurlExpression<E, Order, Space, IsConstant>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<E>::order == 0,
    //               "CurlExpression requires a scalar (order 0) operand");
    // static_assert(ExpressionTraits<E>::space != Space::None,
    //               "CurlExpression requires the operand to be defined in a space");

    static_assert(ExpressionTraits<E>::order == 1,
                  "CurlExpression requires a vector (order 1) operand");

    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order; // Divergence results in a vector (order 1)
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv + 1; // Increment derivative order
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

// --- Partial Specialization 1: Curl of a constant expression ---
template <typename E, size_t Order, size_t Space>
class CurlExpression<E, Order, Space, 1> : public UnaryOperator<CurlExpression<E, Order, Space, 1>>
{
    static_assert(Order == 1,
                  "CurlExpression: Unsupported tensor order. Only vectors (order 1) are supported for curl operations.");


    using Base = UnaryOperator<CurlExpression<E, Order, Space, 1>>;

private:
    mutable gsMatrix<typename Base::Scalar> tmp;
    using Base::expr_;

public:
    CurlExpression(const E& expr)
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

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(gismo::ExpressionHelper<typename Base::Scalar> & helper) const
    {
        GISMO_UNUSED(helper);
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207\u00D7("<<expr_<<")";
    }

};

// --- Partial Specialization 2: Curl of a variable object ---
template <typename E, size_t Order, size_t Space>
class CurlExpression<E, Order, Space, 0> : public UnaryOperator<CurlExpression<E, Order, Space, 0>>
{
    using Base = UnaryOperator<CurlExpression<E, Order, Space, 0>>;

private:
    mutable gsMatrix<typename Base::Scalar> tmp;
    using Base::expr_;

public:
    CurlExpression(const E& expr)
    :
    Base(expr)
    {
        GISMO_ENSURE(expr_.domainDim() == 3, "The domain dimension must be equal to 3 for the curl operator");
        for (short_t d = 0; d != E::order; d++)
            GISMO_ENSURE(expr_.sizes()[d] == 3, "All sizes must equal 3 for the curl operator");
    }

    const std::array<size_t, Base::order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207\u00D7("<<expr_<<")";
    }

};

// Generic factory function for easy creation
template <typename E>
CurlExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> curl(const E& expr)
{
    return CurlExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>(expr);
}

// Curl of curl identity: ∇×(∇×A) = ∇(∇•A) - ∇²A
// Specialized for space=0
template <typename E>
auto curl(const CurlExpression<E, 1, 0, true>& expr)
-> decltype(/* grad(div(expr.expr())) -  */lapl(expr.expr()))
{
    gsWarn<<"Warning: Curl of curl identity is not fully implemented (missing grad(div))!\n";
    return /* grad(div(expr.expr())) -  */lapl(expr.expr());
}

// Partial specialization for addition
template <typename LhsExpr, typename RhsExpr>
auto curl(const AddExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(curl(expr.lhs()) + curl(expr.rhs()))
{
    return curl(expr.lhs()) + curl(expr.rhs());
}

// Partial specialization for subtraction
// ∇×(A - B) = ∇×A - ∇×B
template <typename LhsExpr, typename RhsExpr>
auto curl(const SubtractExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(curl(expr.lhs()) - curl(expr.rhs()))
{
    return curl(expr.lhs()) - curl(expr.rhs());
}

// Partial specialization for multiplication by a scalar
// ∇×(ψA) = ψ(∇×A) + (∇ψ)×A
// For scalar ψ (order 0) and vector A (order 1)
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
auto curl(const ProductExpression<LhsExpr,RhsExpr,0,1,LhsSpace,RhsSpace>& expr)
-> decltype(expr.lhs() * curl(expr.rhs()) + cross(grad(expr.lhs()), expr.rhs()))
{
    return expr.lhs() * curl(expr.rhs()) + cross(grad(expr.lhs()), expr.rhs());
}

// Partial specialization for vector × scalar
// ∇×(Aψ) = ψ(∇×A) + A×(∇ψ)
// For vector A (order 1) and scalar ψ (order 0)
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
auto curl(const ProductExpression<LhsExpr,RhsExpr,1,0,LhsSpace,RhsSpace>& expr)
-> decltype(expr.rhs() * curl(expr.lhs()) + cross(expr.lhs(), grad(expr.rhs())))
{
    return expr.rhs() * curl(expr.lhs()) + cross(expr.lhs(), grad(expr.rhs()));
}

// Partial specialization for cross product
// ∇×(A×B) = A(∇·B) - B(∇·A) + (B·∇)A - (A·∇)B
// For vectors A,B (order 1), result is a vector (order 1)
template <typename LhsExpr, typename RhsExpr>
auto curl(const CrossProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(expr.lhs() * div(expr.rhs()) - expr.rhs() * div(expr.lhs()) + transpose(grad(expr.rhs())) * expr.lhs() - transpose(grad(expr.lhs())) * expr.rhs())
{
    return expr.lhs() * div(expr.rhs()) - expr.rhs() * div(expr.lhs()) + transpose(grad(expr.rhs())) * expr.lhs() - transpose(grad(expr.lhs())) * expr.rhs();
}

// Partial specialization for gradient (second derivative identity)
// ∇×(∇ψ) = 0 (curl of gradient is zero)
// For scalar field ψ (order 0), result is zero vector (order 1)
template <typename E>
auto curl(const GradExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>& expr)
-> ConstantObject<typename ExpressionTraits<E>::Scalar, 1>
{
    GISMO_UNUSED(expr);
    return ConstantObject<typename ExpressionTraits<E>::Scalar, 1>(std::array<size_t, 1>{3},"0");  // Curl of gradient is always zero vector of size 3
}

// Partial specialization for outer product
// ∇×(A⊗B) = (∇×A)⊗B - A×(∇B)
// For vectors A,B (order 1), result is a tensor (order 2)
template <typename LhsExpr, typename RhsExpr>
auto curl(const OuterProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(curl(expr.lhs()) * expr.rhs() - expr.lhs() * grad(expr.rhs()))
{
    return curl(expr.lhs()) * expr.rhs() - expr.lhs() * grad(expr.rhs());
}

// Partial specialization for division of a vector by a scalar
// ∇×(A/φ) = (φ∇×A - ∇φ×A)/φ²
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
auto curl(const DivisionExpression<LhsExpr,RhsExpr,1,0,LhsSpace,RhsSpace>& expr)
-> decltype((expr.rhs() * curl(expr.lhs()) - cross(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * curl(expr.lhs()) - cross(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs());
}

// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Curl of gradient is always zero (this is a valid identity, already implemented above)
// ∇×(∇ψ) = 0 is already implemented

// Curl of divergence is undefined
// ∇×(∇•A) is undefined because divergence of a vector produces a scalar,
// and curl of a scalar is not defined
template <typename E>
auto curl(const DivExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>& expr) -> void
{
    GISMO_ERROR("∇×(∇•) is undefined: curl of divergence is not defined (scalar has no curl)");
}

// Curl of Laplacian is undefined
// ∇×(∇²ψ) is undefined because Laplacian of a scalar produces a scalar,
// and curl of a scalar is not defined
template <typename E>
auto curl(const LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>& expr) -> void
{
    GISMO_ERROR("∇×(∇²) is undefined: curl of Laplacian is not defined (scalar has no curl)");
}

}//namespace Expr
}//namespace gismo
