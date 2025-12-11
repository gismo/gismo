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
template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
struct ExpressionTraits<LaplExpression<_E, _Order, _Space, _IsConstant>>
{
    typedef _E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<_E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_E>::Order; // Laplacian decreases order by 2
    static constexpr size_t Space = ExpressionTraits<_E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<_E>::Deriv + 2; // Increment derivative order
    static constexpr bool IsConstant = ExpressionTraits<_E>::IsConstant;
};

// --- Unified LaplExpression using enable_if for Space-aware eval ---
template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
class LaplExpression : public UnaryOperator<LaplExpression<_E, _Order, _Space, _IsConstant>>
{
    using Base = UnaryOperator<LaplExpression<_E, _Order, _Space, _IsConstant>>;
    typedef typename Base::Scalar T;
    using Base::Order;
private:
    mutable gsMatrix<T> tmp;
    using Base::expr_;
public:
    LaplExpression(const _E& expr)
    :
    Base(expr)
    {
    }

    const std::array<size_t, Order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    // Eval for Space = None (constant)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 1, ExpressionValue<T>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        ExpressionValue<T> result(1, 1);
        tmp.resize(1, 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = None (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 0, ExpressionValue<T>>::type
    eval(const index_t k) const
    {
        // Laplacian of variable expression: sum of second derivatives
        // Return: laplacians.col(k)
        ExpressionValue<T> result(1, 1);
        result(0, 0) = expr_.data().laplacians.col(k);
        return result;
    }

    // Eval for Space = Test or Trial (constant)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 1, ExpressionValue<T>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        ExpressionValue<T> result(1, 1);
        tmp.resize(1, 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 0, ExpressionValue<T>>::type
    eval(const index_t k) const
    {
        // Laplacian for basis functions: returns laplacians per basis function
        const index_t numActive = expr_.data().laplacians.rows();
        
        ExpressionValue<T> result(
            S == SpaceType::Test ? numActive : 1,
            S == SpaceType::Trial ? numActive : 1
        );
        
        // Fill result based on space type
        for (index_t i = 0; i < numActive; ++i)
        {
            if (S == SpaceType::Test)
                result(i, 0) = expr_.data().laplacians.col(k).row(i);
            else // Trial
                result(0, i) = expr_.data().laplacians.col(k).row(i);
        }
        
        return result;
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        helper.add(expr_);
        expr_.data().flags |= NEED_LAPLACIAN;
    }

    void print(std::ostream & os) const
    {
        os<<"\u0394("<<expr_<<")";
    }
};

// Old specializations removed - now using unified class with enable_if

// Generic factory function for easy creation
template <typename _E>
LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> lapl(const _E& expr)
{
    return LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr);
}

// Partial specialization for addition
template <typename _LhsExpr, typename _RhsExpr>
auto lapl(const AddExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(lapl(expr.lhs()) + lapl(expr.rhs()))
{
    return lapl(expr.lhs()) + lapl(expr.rhs()); // ∇²(ψ + φ) = ∇²ψ + ∇²φ
}

// Partial specialization for subtraction
// ∇²(ψ - φ) = ∇²ψ - ∇²φ
template <typename _LhsExpr, typename _RhsExpr>
auto lapl(const SubtractExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(lapl(expr.lhs()) - lapl(expr.rhs()))
{
    return lapl(expr.lhs()) - lapl(expr.rhs());
}

// Partial specialization for product (second derivative identity)
// ∇²(ψφ) = ψ∇²φ + 2∇ψ•∇φ + φ∇²ψ
// For scalar functions ψ,φ (order 0), result is scalar (order 0)
template <typename _LhsExpr, typename _RhsExpr>
auto lapl(const ProductExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(expr.lhs() * lapl(expr.rhs()) + expr.rhs() * lapl(expr.lhs()) + ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs())))
{
    static_assert(_LhsExpr::Order == 0 && _RhsExpr::Order == 0,
                  "Laplacian of product only implemented for scalar functions (order 0)");
    return expr.lhs() * lapl(expr.rhs()) + expr.rhs() * lapl(expr.lhs()) + ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs()));
}

// Partial specialization for division of scalar by scalar
// ∇²(ψ/φ) = (φ∇²ψ - ψ∇²φ - 2∇ψ·∇φ)/φ²
// For scalar functions ψ,φ (order 0), result is scalar (order 0)
template <typename _LhsExpr, typename _RhsExpr>
auto lapl(const DivisionExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype((expr.rhs() * lapl(expr.lhs()) - expr.lhs() * lapl(expr.rhs()) - ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs()))) / (expr.rhs() * expr.rhs()))
{
    static_assert(_LhsExpr::Order == 0 && _RhsExpr::Order == 0,
                  "Laplacian of division only implemented for scalar functions (order 0)");
    return (expr.rhs() * lapl(expr.lhs()) - expr.lhs() * lapl(expr.rhs()) - ConstantObject<double,0>(2.0) * dot(grad(expr.lhs()), grad(expr.rhs()))) / (expr.rhs() * expr.rhs());
}


// Laplacian of gradient produces third derivatives
// ∇²(∇ψ) = ∇(∇²ψ) is a valid operation (third derivative)
// This could be implemented as a third derivative expression if needed
template <typename _E>
auto lapl(const GradExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> decltype(grad(lapl(expr.expr())))
{
    return grad(lapl(expr.expr()));
}


// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Laplacian of divergence is undefined for vector fields
// ∇²(∇•A) is undefined because divergence of a vector produces a scalar,
// and we already have ∇²(scalar), so this would be valid if the input expression is appropriate
template <typename _E>
auto lapl(const DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> void
{
    static_assert(false, "∇²(∇•A) is undefined for vector fields.");
}

// Laplacian of curl is defined: ∇²(∇×A) = ∇×(∇²A)
// This is a valid vector calculus identity
template <typename _E>
auto lapl(const CurlExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> decltype(curl(lapl(expr.expr())))
{
    return curl(lapl(expr.expr()));
}

// Laplacian of Laplacian (fourth derivative)
// ∇²(∇²ψ) = ∇⁴ψ is a valid operation (biharmonic operator)
// However, it is not implemented as a separate expression type
template <typename _E>
auto lapl(const LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> void
{
    // TODO: Take care of 4th order derivative in gsFuncData
    static_assert(false,"∇²(∇²ψ) = ∇⁴ψ is not implemented as a separate expression type.");
    // return lapl(lapl(expr.expr()));
}

template <class T, size_t _Space, size_t _Order, typename _SpaceObject>
auto variation(const LaplExpression<SolutionObject<T,_Space,_Order>,_Order,SpaceType::None,false> & expr,
          const _SpaceObject & space)
-> decltype(lapl(variation(expr.expr(), space)))
{
    return lapl(variation(expr.expr(), space));
}

}//namespace Expr
}//namespace gismo