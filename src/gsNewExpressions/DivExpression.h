/** @file DivExpression.h

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
struct ExpressionTraits<DivExpression<E, Order, Space, IsConstant>>
{
    static_assert(Order == 1,
                "DivExpression: Unsupported tensor order. Only vectors (order 1) are supported for divergence operations.");

    typedef E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<E>::order-1; // Divergence results in a vector (order 1)
    static constexpr size_t space = ExpressionTraits<E>::space;
    static constexpr size_t deriv = ExpressionTraits<E>::deriv + 1; // Increment derivative order
    static constexpr bool isConstant = ExpressionTraits<E>::isConstant;
};

// --- Partial Specialization 1: Divergence of a constant expression ---
template <typename E, size_t Order, size_t Space>
class DivExpression<E, Order, Space, 1> : public UnaryOperator<DivExpression<E, Order, Space, 1>>
{
    static_assert(Order == 1,
                  "DivExpression: Unsupported tensor order. Only vectors (order 1) are supported for divergence operations.");


    using Base = UnaryOperator<DivExpression<E, Order, Space, 1>>;

private:
    std::array<size_t,Base::order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;
    using Base::expr_;

public:
    DivExpression(const E& expr)
    :
    Base(expr)
    {
        // Assumes symmetry!
        for (short_t d=0; d!=ExpressionTraits<E>::order; d++)
        {
            GISMO_ENSURE(expr_.sizes()[d] == expr_.sizes()[0], "All sizes must be equal for the div operator");
            // GISMO_ENSURE(expr_.sizes()[d] == expr_.domainDim(), "All sizes must be equal to the domain dimension for the div operator");
        }

        for (short_t d=1; d!=ExpressionTraits<E>::order; d++)
            sizes_[d-1] = expr_.sizes()[d];
    }


    const std::array<size_t, Base::order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        GISMO_ASSERT(expr_.domainDim() == 0,"A constant expression should have domain dimension 0");
        return 0;
    }

    gsMatrix<typename Base::Scalar> eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        tmp.resize(expr_.domainDim(),1);
        tmp.setZero();
        return tmp;
    }

    void parse(ExpressionHelper<typename Base::Scalar> & helper) const
    {
        GISMO_UNUSED(helper);
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207\u2027("<<expr_<<")";
    }

};

// --- Partial Specialization 2: Divergence of a variable object ---
template <typename E, size_t Order, size_t Space>
class DivExpression<E, Order, Space, 0> : public UnaryOperator<DivExpression<E, Order, Space, 0>>
{
    using Base = UnaryOperator<DivExpression<E, Order, Space, 0>>;
private:
    std::array<size_t,Base::order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;
    using Base::expr_;

public:
    DivExpression(const E & expr)
    :
    Base(expr)
    {
        // Assumes symmetry!
        for (short_t d=0; d!=ExpressionTraits<E>::order; d++)
        {
            GISMO_ENSURE(expr_.sizes()[d] == expr_.sizes()[0], "All sizes must be equal for the div operator");
            GISMO_ENSURE(expr_.sizes()[d] == expr_.domainDim(), "All sizes must be equal to the domain dimension for the div operator");
        }

        for (short_t d=1; d!=ExpressionTraits<E>::order; d++)
            sizes_[d-1] = expr_.sizes()[d];
    }

    const std::array<size_t, Base::order> & sizes() const
    {
        return sizes_;
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
        os<<"\u2207\u2027("<<expr_<<")";
    }

};

// Generic factory function for easy creation
template <typename E>
DivExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>
div(const E& expr)
{
    return DivExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>(expr);
}

// Partial specialization for addition
// ∇•(X + Y) = ∇•X + ∇•Y
template <typename LhsExpr, typename RhsExpr>
auto div(const AddExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(div(expr.lhs()) + div(expr.rhs()))
{
    return (div(expr.lhs()) + div(expr.rhs()));
}

// Partial specialization for subtraction
// ∇•(X - Y) = ∇•X - ∇•Y
template <typename LhsExpr, typename RhsExpr>
auto div(const SubtractExpression<LhsExpr,RhsExpr,ExpressionTraits<LhsExpr>::order,ExpressionTraits<RhsExpr>::order,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space> expr)
-> decltype(div(expr.lhs()) - div(expr.rhs()))
{
    return (div(expr.lhs()) - div(expr.rhs()));
}

// Partial specialization for multiplication of a scalar and a vector
// ∇•(fV) = ∇f • V + f ∇•V
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
auto div(const ProductExpression<LhsExpr,RhsExpr,0,1,LhsSpace,RhsSpace>& expr)
-> decltype(dot(grad(expr.lhs()), expr.rhs()) + expr.lhs() * div(expr.rhs()))
{
    return dot(grad(expr.lhs()), expr.rhs()) + expr.lhs() * div(expr.rhs());
}

// Partial specialization for cross product
// ∇•(A×B) = (∇×A)•B - A•(∇×B) = B•(∇×A) - A•(∇×B)
// For vectors A,B (order 1), result is a scalar (order 0)
template <typename LhsExpr, typename RhsExpr>
auto div(const CrossProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(dot(expr.rhs(), curl(expr.lhs())) - dot(expr.lhs(), curl(expr.rhs())))
{
    return dot(expr.rhs(), curl(expr.lhs())) - dot(expr.lhs(), curl(expr.rhs()));
}

// Partial specialization for outer product
// ∇•(A⊗B) = (∇•A)B + (A•∇)B
// For vectors A,B (order 1), result is a vector (order 1)
template <typename LhsExpr, typename RhsExpr>
auto div(const OuterProductExpression<LhsExpr,RhsExpr,1,1,ExpressionTraits<LhsExpr>::space,ExpressionTraits<RhsExpr>::space>& expr)
-> decltype(div(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs()))
{
    GISMO_UNUSED(expr);
    return div(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs());
}

// Partial specialization for gradient (second derivative identity)
// ∇•(∇ψ) = ∇²ψ (Laplacian)
// For scalar field ψ (order 0), result is a scalar (order 0)
template <typename E>
auto div(const GradExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>& expr)
-> LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>
{
    return LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant>(expr.expr());
}

// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Divergence of gradient produces Laplacian, not undefined, but divergence of divergence is undefined
// ∇•(∇•A) is undefined because divergence of a vector produces a scalar,
// and divergence of a scalar is not defined
template <typename E>
auto div(const DivExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr) -> void
{
    GISMO_ERROR("∇•(∇•) is undefined: divergence of divergence is not defined (scalar has no divergence)");
}

// Divergence of curl is always zero (this is a valid identity)
// ∇•(∇×A) = 0, but for completeness, we could add this as a zero expression
// This is actually defined and equals zero, so we should implement it properly
template <typename E>
auto div(const CurlExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr)
-> ConstantObject<typename ExpressionTraits<E>::Scalar, 0>
{
    GISMO_UNUSED(expr);
    return ConstantObject<typename ExpressionTraits<E>::Scalar, 0>(std::array<size_t, 0>{},"0");  // Divergence of curl is always zero scalar
}

// Divergence of Laplacian is undefined
// ∇•(∇²ψ) is undefined because Laplacian of a scalar produces a scalar,
// and divergence of a scalar is not defined
template <typename E>
auto div(const LaplExpression<E, ExpressionTraits<E>::order, ExpressionTraits<E>::space, ExpressionTraits<E>::isConstant> expr) -> void
{
    GISMO_ERROR("∇•(∇²) is undefined: divergence of Laplacian is not defined (scalar has no divergence)");
}

// Partial specialization for division of a vector by a scalar
// ∇•(A/φ) = (φ∇•A - ∇φ•A)/φ²
template <typename LhsExpr, typename RhsExpr, size_t LhsSpace, size_t RhsSpace>
auto div(const DivisionExpression<LhsExpr,RhsExpr,1,0,LhsSpace,RhsSpace>& expr)
-> decltype((expr.rhs() * div(expr.lhs()) - dot(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * div(expr.lhs()) - dot(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs());
}

// --- Partial Specialization 2: Divergence of a VariableObject ---



// // --- Partial Specialization 1: Addition of two expressions of the SAME ORDER (X + X) ---
// template <typename LhsExpr, typename RhsExpr>
// class AddExpression<ArrayExpression<LhsExpr>, RhsExpr,
//     typename std::enable_if<0 == (RhsExpr::order)>::type // Simplified condition
// > : public BaseObject<LhsExpr,
//                           typename LhsExpr::Scalar,
//                           LhsExpr::order,
//                           LhsExpr::isConstant && RhsExpr::isConstant,
//                           0> // Use LhsExpr's Scalar and Order directly
// {
// public:
// // Scalar and Order are directly from LhsExpr/RhsExpr
//     typedef typename LhsExpr::Scalar Scalar;
//     static constexpr int order = LhsExpr::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<RhsExpr, typename LhsExpr::Scalar, order, LhsExpr::isConstant && RhsExpr::isConstant, 0>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
//           lhs_expr_(lhs),
//           rhs_expr_(rhs)
//     {
//     }

//     gsMatrix<Scalar> eval(const index_t k) const
//     {
//         gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
//         gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
//         lhs_val.array() += rhs_val; // Element-wise addition
//         return lhs_val; // Return the modified lhs_val
//     }

// private:
//     const LhsExpr& lhs_expr_;
//     const RhsExpr& rhs_expr_;
// };



// // --- Partial Specialization 2: Scalar (Order 0) + Higher Order (Order N > 0) ---
// template <typename LhsExpr, typename RhsExpr>
// class AddExpression<LhsExpr, RhsExpr,
//     typename std::enable_if<(LhsExpr::order == 0) && (RhsExpr::order > 0)>::type // Simplified condition
// > : public BaseObject<RhsExpr, typename LhsExpr::Scalar, RhsExpr::order> // Base on RhsExpr for Scalar and Order
// {
// public:
//     typedef typename RhsExpr::Scalar Scalar;
//     static constexpr int order = RhsExpr::order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<RhsExpr, typename LhsExpr::Scalar, order>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
//           lhs_expr_(lhs),
//           rhs_expr_(rhs)
//     {}

//     gsMatrix<Scalar> eval(const index_t k) const
//     {
//         gsMatrix<Scalar> lhs_val = lhs_expr_.eval(k);
//         gsMatrix<Scalar> rhs_val = rhs_expr_.eval(k);
//         rhs_val.array() += lhs_val.value(); // Add scalar to each element of lhs_val
//         return rhs_val;
//     }

// private:
//     const LhsExpr& lhs_expr_;
//     const RhsExpr& rhs_expr_;
// };


// // Scalar + Scalar
// AddExpression<ScalarExpression<real_t>, ScalarExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, ScalarExpression<real_t>>(lhs, rhs);
// }

// // Scalar + Vector
// AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const VectorExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>>(lhs, rhs);
// }

// // Vector + Scalar
// AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>> operator+(const VectorExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, VectorExpression<real_t>>(rhs, lhs);
// }

// // Vector + Vector
// AddExpression<VectorExpression<real_t>, VectorExpression<real_t>> operator+(const VectorExpression<real_t>& lhs, const VectorExpression<real_t>& rhs)
// {
//     return AddExpression<VectorExpression<real_t>, VectorExpression<real_t>>(lhs, rhs);
// }

// // Scalar + Matrix
// AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>> operator+(const ScalarExpression<real_t>& lhs, const MatrixExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>>(lhs, rhs);
// }

// // Matrix + Scalar
// AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>> operator+(const MatrixExpression<real_t>& lhs, const ScalarExpression<real_t>& rhs)
// {
//     return AddExpression<ScalarExpression<real_t>, MatrixExpression<real_t>>(rhs, lhs);
// }

}//namespace Expr
}//namespace gismo