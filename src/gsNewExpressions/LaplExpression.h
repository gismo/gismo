/** @file LaplExpression.h

    @brief Laplacian expression class

    Computes the Laplacian (trace of Hessian) of an expression. For a scalar
    field f(x), lapl(f) = ∇²f = ∂²f/∂x² + ∂²f/∂y² + ... 
    For a vector field u(x), lapl(u) applies Laplacian component-wise.

    Key features:
    - Preserves tensor order (scalar → scalar, vector → vector)
    - Requires second derivatives (Deriv + 2)
    - Preserves Space type from operand
    - For vectors, requires ComponentLaplExpression to handle each component

    Why ComponentLaplExpression is needed:
    The Laplacian of a vector field [u1, u2, ...] must be computed as
    [lapl(u1), lapl(u2), ...]. We cannot directly compute the trace of the
    Jacobian's Hessian. ComponentLaplExpression extracts each component,
    applies Laplacian, then recombines them.

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
template <typename _E, size_t _Order, enum SpaceType _Space, size_t _IsConstant>
struct ExpressionTraits<LaplExpression<_E, _Order, _Space, _IsConstant>>
{
    typedef _E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<_E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_E>::Order; // Laplacian preserves order (scalar->scalar, handled by specializations for vectors)
    static constexpr SpaceType Space = ExpressionTraits<_E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<_E>::Deriv + 2; // Increment derivative order
    static constexpr bool IsConstant = ExpressionTraits<_E>::IsConstant;
};

// --- Unified LaplExpression using enable_if for Space-aware eval ---
template <typename _E, size_t _Order, enum SpaceType _Space, size_t _IsConstant>
class LaplExpression : public UnaryOperator<LaplExpression<_E, _Order, _Space, _IsConstant>>
{
    using Base = UnaryOperator<LaplExpression<_E, _Order, _Space, _IsConstant>>;
private:
    std::array<size_t, _Order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;

public:
    LaplExpression(const _E& expr)
    :
    Base(expr)
    {
        // Copy sizes from the underlying expression (Laplacian preserves tensor structure)
        for (size_t d = 0; d < _Order; ++d)
            sizes_[d] = this->expr_.sizes()[d];
    }

    const std::array<size_t, _Order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return this->expr_.domainDim();
    }

    // Eval for Space = None (constant)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 1, ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        ExpressionResult<typename Base::Scalar> result(1, 1);
        tmp.resize(1, 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = None (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 0, ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Laplacian of variable expression: sum of second derivatives
        // Return: laplacians.col(k)
        ExpressionResult<typename Base::Scalar> result(1, 1);
        result(0, 0) = this->expr_.data().laplacians.col(k);
        return result;
    }

    // Eval for Space = Test or Trial (constant)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 1, ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        ExpressionResult<typename Base::Scalar> result(1, 1);
        tmp.resize(1, 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 0, ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Laplacian for basis functions: returns laplacians per basis function
        const index_t numActive = this->expr_.data().laplacians.rows();
        
        ExpressionResult<typename Base::Scalar> result(
            S == SpaceType::Test ? numActive : 1,
            S == SpaceType::Trial ? numActive : 1
        );
        
        // Fill result based on space type
        for (index_t i = 0; i < numActive; ++i)
        {
            if (S == SpaceType::Test)
                result(i, 0) = this->expr_.data().laplacians.col(k).row(i);
            else // Trial
                result(0, i) = this->expr_.data().laplacians.col(k).row(i);
        }
        
        return result;
    }

    // Helper to detect whether ExprType has expr() method
    template <typename U>
    static auto has_expr_test(int) -> decltype(std::declval<const U>().expr(), char(0));
    template <typename U>
    static int has_expr_test(...);
    static constexpr bool has_expr = std::is_same<decltype(has_expr_test<_E>(0)), char>::value;

    // parse implementation for expressions that expose .expr()
    void parse_impl(gismo::ExpressionHelper<typename Base::Scalar> & helper, TrueTag) const
    {
        // Set derivative order first
        this->expr_.expr().setDerivative(Base::Deriv);
        // Parse the underlying expression
        this->expr_.expr().parse(helper);
        // Ensure laplacian flag is set
        this->expr_.expr().data().flags |= NEED_LAPLACIAN;
    }

    // parse implementation for expressions that are variable objects directly
    void parse_impl(gismo::ExpressionHelper<typename Base::Scalar> & helper, FalseTag) const
    {
        // Set derivative order first
        this->expr_.setDerivative(Base::Deriv);
        // Parse the underlying expression
        this->expr_.parse(helper);
        // Ensure laplacian flag is set
        this->expr_.data().flags |= NEED_LAPLACIAN;
    }

    void parse(gismo::ExpressionHelper<typename Base::Scalar> & helper) const
    {
        parse_impl(helper, typename BoolToTag<has_expr>::type());
    }

    void print(std::ostream & os) const
    {
        os<<"\u0394("<<this->expr_<<")";
    }
};

// Old specializations removed - now using unified class with enable_if

// Generic factory function for easy creation
// template <typename _E>
// LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> lapl(const _E& expr)
// {
//     return LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr);
// }

// Partial specialization for scalars
template <typename _E>
typename std::enable_if<ExpressionTraits<_E>::Order == 0, LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> >::type
lapl(const _E& expr)
{
    gsDebug<<"Using scalar Laplacian\n";
    return LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr);
}

template <typename _E>
typename std::enable_if<ExpressionTraits<_E>::Order == 2, LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> >::type
lapl(const _E& expr)
{
    gsDebug<<"Using MATRIX Laplacian\n";
    return LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr);
}

// Partial specialization for vectors
// Safer implementation for VariableObject vector fields: sum component-wise Laplacians
// This avoids complex curl/grad identities and return-type SFINAE issues
template <typename T, bool _IsConstant>
auto lapl(const VariableObject<T, 1, _IsConstant> & expr)
{
    gsDebug<<"Using vector Laplacian (component-wise sum)\n";
    // Build component expressions and sum their scalar laplacians
    auto c0 = ComponentExpression<VariableObject<T,1,_IsConstant>,1>(expr,0);
    auto c1 = ComponentExpression<VariableObject<T,1,_IsConstant>,1>(expr,1);
    return lapl(c0) + lapl(c1);
}

// Forward declaration for ComponentLaplExpression
template <typename E, size_t _Order> class ComponentLaplExpression;

// ExpressionTraits for ComponentLaplExpression
template <typename E, size_t _Order>
struct ExpressionTraits<ComponentLaplExpression<E, _Order>>
{
    typedef E ExprType;
    typedef typename ExpressionTraits<E>::Scalar Scalar;
    static constexpr size_t Order = 0; // Laplacian of component is scalar
    static constexpr SpaceType Space = ExpressionTraits<E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<E>::Deriv + 2;
    static constexpr bool IsConstant = ExpressionTraits<E>::IsConstant;
};

/**
 * \brief Laplacian of a component of a vector expression
 * 
 * This specialized expression is needed because computing lapl(u[i]) where u
 * is a vector field requires extracting the i-th component's laplacian from
 * the parent vector field's second derivative data (values[2]).
 * 
 * Why we need a separate class:
 * - The generic LaplExpression works on the expression's laplacians member
 * - For ComponentExpression<E,1>, there is no direct laplacians data
 * - Instead, we must access the parent's values[2] and extract the relevant
 *   second derivatives for the specific component
 * 
 * Memory layout in values[2] for a 2D vector field in 2D domain:
 * [d²f₀/dx², d²f₀/dy², d²f₀/dxdy, d²f₁/dx², d²f₁/dy², d²f₁/dxdy]
 * 
 * The laplacian of component i is: Σⱼ d²fᵢ/dxⱼ² (sum of diagonal second derivs)
 */
template <typename E, size_t _Order>
class ComponentLaplExpression : public BaseExpression<ComponentLaplExpression<E, _Order>>
{
    using Base = BaseExpression<ComponentLaplExpression<E, _Order>>;
    typedef typename Base::Scalar T;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;
    
    typedef typename ExpressionTraits<ComponentLaplExpression<E, _Order>>::ExprType ExprType;
    
    const ExprType& expr() const { return expr_; }
    
    ComponentLaplExpression(const ComponentExpression<E, _Order>& compExpr)
    : BaseExpression<ComponentLaplExpression<E, _Order>>(),
      expr_(compExpr.expr()),
      component_(compExpr.component())
    {
        static_assert(_Order == 1, "ComponentLaplExpression only for vector components");
    }
    
    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        // Set derivative order first (Deriv = parent's Deriv + 2)
        expr_.setDerivative(Base::Deriv);
        // Parse the underlying expression
        expr_.parse(helper);
        // Request second derivatives - we compute component laplacian from values[2]
        expr_.data().flags |= NEED_DERIV2;
    }
    
    void print(std::ostream & os) const
    {
        os << "\u0394(" << expr_ << "[" << component_ << "])";
    }
    
    std::array<size_t, 0> sizes() const { return std::array<size_t, 0>(); }
    
    size_t domainDim() const { return expr_.domainDim(); }
    
    ExpressionResult<T> eval(const index_t k) const
    {
        // Compute component-wise laplacian from second derivatives
        // For a 2D domain with 2D target (vector field), values[2] has 6 rows per point:
        // [d²f₀/dx², d²f₀/dy², d²f₀/dxdy, d²f₁/dx², d²f₁/dy², d²f₁/dxdy]
        // The laplacian of component i is: d²fᵢ/dx² + d²fᵢ/dy² (sum of first dim.first diagonal entries)
        
        const auto& data = expr_.data();
        const index_t d = data.dim.first;  // domain dimension
        
        // deriv2Size per component (number of unique second derivatives)
        // For 2D: 3 entries (d²/dx², d²/dy², d²/dxdy)
        const index_t d2sz = d * (d + 1) / 2;
        
        // Starting row for component in values[2]
        const index_t start = component_ * d2sz;
        
        // Sum the diagonal second derivatives (first d entries for this component)
        T lapl_val = 0;
        for (index_t i = 0; i < d; ++i)
            lapl_val += data.values[2](start + i, k);
        
        ExpressionResult<T> result(1, 1);
        gsMatrix<T> tmp(1, 1);
        tmp(0, 0) = lapl_val;
        result(0, 0) = tmp;
        return result;
    }
    
    // Forward test/trial to parent expression
    auto test() const -> decltype(std::declval<const ComponentLaplExpression&>().expr().test())
    {
        return expr_.test();
    }
    
    auto trial() const -> decltype(std::declval<const ComponentLaplExpression&>().expr().trial())
    {
        return expr_.trial();
    }
    
protected:
    const ExprType expr_;
    const index_t component_;
};

// Specialized lapl() for ComponentExpression - extracts the component's laplacian from parent
template <typename E>
ComponentLaplExpression<E, 1> lapl(const ComponentExpression<E, 1>& expr)
{
    gsDebug << "Using ComponentExpression Laplacian\n";
    return ComponentLaplExpression<E, 1>(expr);
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

// Helper to make static_assert dependent so overload resolution can proceed
template <typename T>
struct dependent_false { static constexpr bool value = false; };

// Laplacian of divergence is undefined for vector fields
// ∇²(∇•A) is undefined because divergence of a vector produces a scalar,
// and we already have ∇²(scalar), so this would be valid if the input expression is appropriate
template <typename _E>
auto lapl(const DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> void
{
    static_assert(dependent_false<_E>::value, "∇²(∇•A) is undefined for vector fields.");
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
    static_assert(dependent_false<_E>::value,"∇²(∇²ψ) = ∇⁴ψ is not implemented as a separate expression type.");
    // return lapl(lapl(expr.expr()));
}

template <class T, enum SpaceType _Space, size_t _Order, typename _SpaceObject>
auto variation(const LaplExpression<SolutionObject<T,_Space,_Order>,_Order,SpaceType::None,false> & expr,
          const _SpaceObject & space)
-> decltype(lapl(variation(expr.expr(), space)))
{
    return lapl(variation(expr.expr(), space));
}

}//namespace Expr
}//namespace gismo