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

template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
struct ExpressionTraits<CurlExpression<_E, _Order, _Space, _IsConstant>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<E>::Order == 0,
    //               "CurlExpression requires a scalar (order 0) operand");
    // static_assert(ExpressionTraits<E>::Space != SpaceType::None,
    //               "CurlExpression requires the operand to be defined in a space");

    static_assert(ExpressionTraits<_E>::Order == 1,
                  "CurlExpression requires a vector (order 1) operand");

    typedef _E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<_E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_E>::Order; // Divergence results in a vector (order 1)
    static constexpr size_t Space = ExpressionTraits<_E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<_E>::Deriv + 1; // Increment derivative order
    static constexpr bool IsConstant = ExpressionTraits<_E>::IsConstant;
};

// --- Unified CurlExpression using enable_if for Space-aware eval ---
template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
class CurlExpression : public UnaryOperator<CurlExpression<_E, _Order, _Space, _IsConstant>>
{
    static_assert(_Order == 1,
                  "CurlExpression: Unsupported tensor order. Only vectors (order 1) are supported for curl operations.");

    using Base = UnaryOperator<CurlExpression<_E, _Order, _Space, _IsConstant>>;

private:
    mutable gsMatrix<typename Base::Scalar> tmp;
    using Base::expr_;

public:
    CurlExpression(const _E& expr)
    :
    Base(expr)
    {
        GISMO_ENSURE(expr_.domainDim() == 3, "The domain dimension must be equal to 3 for the curl operator");
        for (short_t d = 0; d != _Order; d++)
            GISMO_ENSURE(expr_.sizes()[d] == 3, "All sizes must equal 3 for the curl operator");
    }

    const std::array<size_t, Base::Order> & sizes() const
    {
        return expr_.sizes();
    }

    size_t domainDim() const
    {
        return expr_.domainDim();
    }

    // Eval for Space = None (constant)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 1, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        ExpressionValue<typename Base::Scalar> result(1, 1);
        tmp.resize(3, 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = None (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 0, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Curl computation for 3D vector field: curl = ∇ × V
        // Result is a 3D vector
        const index_t numActive = expr_.data().values[0].rows();
        gsAsConstMatrix<typename Base::Scalar, Dynamic, Dynamic> pd =
            expr_.data().values[1].reshapeCol(k, 3, numActive);
        
        tmp.resize(3, 1);
        tmp(0) = pd.row(2).sum();  // ∂V_z/∂y - ∂V_y/∂z
        tmp(1) = pd.row(0).sum();  // ∂V_x/∂z - ∂V_z/∂x
        tmp(2) = pd.row(1).sum();  // ∂V_y/∂x - ∂V_x/∂y
        
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (constant)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 1, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        ExpressionValue<typename Base::Scalar> result(1, 1);
        tmp.resize(3, 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 0, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Curl for basis functions: compute curl for each basis function
        const index_t numActive = expr_.data().values[0].rows();
        gsAsConstMatrix<typename Base::Scalar, Dynamic, Dynamic> pd =
            expr_.data().values[1].reshapeCol(k, 3, numActive);
        
        // Compute curl components for each basis function
        // curl(V) = (∂V_z/∂y - ∂V_y/∂z, ∂V_x/∂z - ∂V_z/∂x, ∂V_y/∂x - ∂V_x/∂y)
        ExpressionValue<typename Base::Scalar> result(
            S == SpaceType::Test ? numActive : 1,
            S == SpaceType::Trial ? numActive : 1
        );
        
        for (index_t i = 0; i < numActive; ++i)
        {
            tmp.resize(3, 1);
            tmp(0) = -pd(2, i);
            tmp(1) = pd(2, i);
            tmp(2) = -pd(0, i);
            tmp(0) += pd(1, i);
            tmp(1) += -pd(0, i);
            tmp(2) += pd(0, i);
            
            if (S == SpaceType::Test)
                result(i, 0) = tmp;
            else // Trial
                result(0, i) = tmp;
        }
        
        return result;
    }

    void parse(gismo::ExpressionHelper<typename Base::Scalar> & helper) const
    {
        helper.add(expr_);
        expr_.data().flags |= NEED_GRAD;
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207\u00D7("<<expr_<<")";
    }

};

// Old specializations removed - now using unified class with enable_if

// Generic factory function for easy creation
template <typename _E>
CurlExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> curl(const _E& expr)
{
    return CurlExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr);
}

// Curl of curl identity: ∇×(∇×A) = ∇(∇•A) - ∇²A
// Specialized for Space=0
template <typename _E>
auto curl(const CurlExpression<_E, 1, 0, true>& expr)
-> decltype(/* grad(div(expr.expr())) -  */lapl(expr.expr()))
{
    gsWarn<<"Warning: Curl of curl identity is not fully implemented (missing grad(div))!\n";
    return /* grad(div(expr.expr())) -  */lapl(expr.expr());
}

// Partial specialization for addition
template <typename _LhsExpr, typename _RhsExpr>
auto curl(const AddExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(curl(expr.lhs()) + curl(expr.rhs()))
{
    return curl(expr.lhs()) + curl(expr.rhs());
}

// Partial specialization for subtraction
// ∇×(A - B) = ∇×A - ∇×B
template <typename _LhsExpr, typename _RhsExpr>
auto curl(const SubtractExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(curl(expr.lhs()) - curl(expr.rhs()))
{
    return curl(expr.lhs()) - curl(expr.rhs());
}

// Partial specialization for multiplication by a scalar
// ∇×(ψA) = ψ(∇×A) + (∇ψ)×A
// For scalar ψ (order 0) and vector A (order 1)
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsSpace, size_t _RhsSpace>
auto curl(const ProductExpression<_LhsExpr,_RhsExpr,0,1,_LhsSpace,_RhsSpace>& expr)
-> decltype(expr.lhs() * curl(expr.rhs()) + cross(grad(expr.lhs()), expr.rhs()))
{
    return expr.lhs() * curl(expr.rhs()) + cross(grad(expr.lhs()), expr.rhs());
}

// Partial specialization for vector × scalar
// ∇×(Aψ) = ψ(∇×A) + A×(∇ψ)
// For vector A (order 1) and scalar ψ (order 0)
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsSpace, size_t _RhsSpace>
auto curl(const ProductExpression<_LhsExpr,_RhsExpr,1,0,_LhsSpace,_RhsSpace>& expr)
-> decltype(expr.rhs() * curl(expr.lhs()) + cross(expr.lhs(), grad(expr.rhs())))
{
    return expr.rhs() * curl(expr.lhs()) + cross(expr.lhs(), grad(expr.rhs()));
}

// Partial specialization for cross product
// ∇×(A×B) = A(∇·B) - B(∇·A) + (B·∇)A - (A·∇)B
// For vectors A,B (order 1), result is a vector (order 1)
template <typename _LhsExpr, typename _RhsExpr>
auto curl(const CrossProductExpression<_LhsExpr,_RhsExpr,1,1,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(expr.lhs() * div(expr.rhs()) - expr.rhs() * div(expr.lhs()) + transpose(grad(expr.rhs())) * expr.lhs() - transpose(grad(expr.lhs())) * expr.rhs())
{
    return expr.lhs() * div(expr.rhs()) - expr.rhs() * div(expr.lhs()) + transpose(grad(expr.rhs())) * expr.lhs() - transpose(grad(expr.lhs())) * expr.rhs();
}

// Partial specialization for gradient (second derivative identity)
// ∇×(∇ψ) = 0 (curl of gradient is zero)
// For scalar field ψ (order 0), result is zero vector (order 1)
template <typename _E>
auto curl(const GradExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>& expr)
-> ConstantObject<typename ExpressionTraits<_E>::Scalar, 1>
{
    GISMO_UNUSED(expr);
    return ConstantObject<typename ExpressionTraits<_E>::Scalar, 1>(std::array<size_t, 1>{3},"0");  // Curl of gradient is always zero vector of size 3
}

// Partial specialization for outer product
// ∇×(A⊗B) = (∇×A)⊗B - A×(∇B)
// For vectors A,B (order 1), result is a tensor (order 2)
template <typename _LhsExpr, typename _RhsExpr>
auto curl(const OuterProductExpression<_LhsExpr,_RhsExpr,1,1,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(curl(expr.lhs()) * expr.rhs() - expr.lhs() * grad(expr.rhs()))
{
    return curl(expr.lhs()) * expr.rhs() - expr.lhs() * grad(expr.rhs());
}

// Partial specialization for division of a vector by a scalar
// ∇×(A/φ) = (φ∇×A - ∇φ×A)/φ²
template <typename _LhsExpr, typename _RhsExpr, size_t _LhsSpace, size_t _RhsSpace>
auto curl(const DivisionExpression<_LhsExpr,_RhsExpr,1,0,_LhsSpace,_RhsSpace>& expr)
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
template <typename _E>
auto curl(const DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>& expr) -> void
{
    GISMO_ERROR("∇×(∇•) is undefined: curl of divergence is not defined (scalar has no curl)");
}

// Curl of Laplacian is undefined
// ∇×(∇²ψ) is undefined because Laplacian of a scalar produces a scalar,
// and curl of a scalar is not defined
template <typename _E>
auto curl(const LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>& expr) -> void
{
    GISMO_ERROR("∇×(∇²) is undefined: curl of Laplacian is not defined (scalar has no curl)");
}

template <class T, size_t _Space, size_t _Order, typename _SpaceObject>
auto variation(const CurlExpression<SolutionObject<T,_Space,_Order>,_Order,SpaceType::None,false> & expr,
          const _SpaceObject & space)
-> decltype(curl(variation(expr.expr(), space)))
{
    return curl(variation(expr.expr(), space));
}

}//namespace Expr
}//namespace gismo
