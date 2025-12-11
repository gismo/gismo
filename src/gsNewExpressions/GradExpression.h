/** @file GradExpression.h

    @brief Gradient expression class

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
struct ExpressionTraits<GradExpression<_E, _Order, _Space, _IsConstant>>
{
    // Static assertions to ensure compatibility
    // static_assert(ExpressionTraits<_E>::Order == 0,
    //               "GradExpression requires a scalar (order 0) operand");
    // static_assert(ExpressionTraits<_E>::Space != SpaceType::None,
    //               "GradExpression requires the operand to be defined in a space");

    typedef _E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<_E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_E>::Order + 1; // Gradient increases order by 1
    static constexpr size_t Space = ExpressionTraits<_E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<_E>::Deriv + 1; // Increment derivative order
    static constexpr bool IsConstant = _IsConstant;
};

// --- Generic GradExpression using SFINAE-based eval functions ---
template <typename _E, size_t _Order, size_t _Space, size_t _IsConstant>
class GradExpression : public UnaryOperator<GradExpression<_E, _Order, _Space, _IsConstant>>
{
    using Base = UnaryOperator<GradExpression<_E, _Order, _Space, _IsConstant>>;
private:
    std::array<size_t,Base::Order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;

public:
    GradExpression(const _E& expr)
    : Base(expr)
    {
        sizes_[0] = this->expr_.domainDim();
        for (short_t d=0; d!=ExpressionTraits<_E>::Order; d++)
            sizes_[d+1] = this->expr_.sizes()[d];
    }

    const std::array<size_t, Base::Order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return this->expr_.domainDim();
    }

    // Eval for Space = None
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 1, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        // Constant expression with no space: return zero gradient
        ExpressionValue<typename Base::Scalar> result(1, 1);
        tmp.resize(this->expr_.domainDim(), 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = None (variable)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 0, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Variable expression with no space: compute gradient from derivatives
        // Assumes: derivatives are in expr_.data().values[1]
        const index_t numActive = this->expr_.data().values[0].rows();
        tmp = this->expr_.data().values[1].reshapeCol(k, this->domainDim(), numActive).transpose();
        
        ExpressionValue<typename Base::Scalar> result(1, 1);
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (constant - rare case)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 1, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        GISMO_UNUSED(k);
        // Constant expression with space: return zero gradient
        ExpressionValue<typename Base::Scalar> result(1, 1);
        tmp.resize(this->expr_.domainDim(), 1);
        tmp.setZero();
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (variable - main case for basis functions)
    template <size_t S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 0, ExpressionValue<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Variable expression with space: compute gradient for each basis function
        // Assumes: derivatives are in expr_.data().values[1]
        const index_t numActive = this->expr_.data().values[0].rows();
        tmp = this->expr_.data().values[1].reshapeCol(k, this->domainDim(), numActive).transpose();
        
        // For Space = Test/Trial, cardinality affects result organization
        ExpressionValue<typename Base::Scalar> result(
            S == SpaceType::Test ? numActive : 1,
            S == SpaceType::Trial ? numActive : 1
        );
        
        // Fill result based on space type
        for (index_t i = 0; i < numActive; ++i)
        {
            if (S == SpaceType::Test)
                result(i, 0) = tmp.row(i);
            else // Trial
                result(0, i) = tmp.row(i);
        }
        
        return result;
    }

    void parse(gismo::ExpressionHelper<typename Base::Scalar> & helper) const
    {
        helper.add(this->expr_);
        this->expr_.data().flags |= NEED_GRAD;
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207("<<this->expr_<<")";
    }
};

// Partial specialization for dot product with space=0
// ∇(A·B) = (A·∇)B + (B·∇)A + A×(∇×B) + B×(∇×A)
// For vectors A,B (order 1), result is a vector (order 1)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const InnerProductExpression<_LhsExpr,_RhsExpr,1,1,0,0>& expr)
-> decltype(dot(expr.lhs(), grad(expr.rhs())) + dot(expr.rhs(), grad(expr.lhs())) + cross(expr.lhs(), curl(expr.rhs())) + cross(expr.rhs(), curl(expr.lhs())))
{
    return (dot(expr.lhs(), grad(expr.rhs())) + dot(expr.rhs(), grad(expr.lhs())) + cross(expr.lhs(), curl(expr.rhs())) + cross(expr.rhs(), curl(expr.lhs())));
}

// Generic factory function for easy creation
template <typename E>
GradExpression<E, ExpressionTraits<E>::Order, ExpressionTraits<E>::Space, ExpressionTraits<E>::IsConstant> grad(const E& expr)
{
    return GradExpression<E, ExpressionTraits<E>::Order, ExpressionTraits<E>::Space, ExpressionTraits<E>::IsConstant>(const_cast<E&>(expr));
}

// Partial specialization for addition
// ∇(X+Y) = ∇X + ∇Y
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const AddExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,ExpressionTraits<_LhsExpr>::Space,_RhsExpr::Space>& expr)
-> decltype(grad(expr.lhs()) + grad(expr.rhs()))
{
    return grad(expr.lhs()) + grad(expr.rhs());
}

// Partial specialization for multiplication
// ∇(X*Y) = (∇X)*Y + X*(∇Y)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const ProductExpression<_LhsExpr,_RhsExpr,
                                    _LhsExpr::Order,
                                    _RhsExpr::Order,
                                    ExpressionTraits<_LhsExpr>::Space,
                                    _RhsExpr::Space> expr)
-> decltype(grad(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs()))
{
    return grad(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs());
}

// Partial specialization for multiplication of scalars
// ∇(f*g) = (∇f)*g + f*(∇g)
// where f,g are scalar functions (order 0)
// and ∇f,∇g are vector functions (order 1)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const ProductExpression<_LhsExpr,_RhsExpr,1,1,ExpressionTraits<_LhsExpr>::Space,_RhsExpr::Space> expr)
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
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const ProductExpression<_LhsExpr,_RhsExpr,0,1,ExpressionTraits<_LhsExpr>::Space,_RhsExpr::Space> expr)
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
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const ProductExpression<_LhsExpr,_RhsExpr,1,0,ExpressionTraits<_LhsExpr>::Space,_RhsExpr::Space> expr)
-> decltype(outer(grad(expr.lhs()),expr.rhs()) + expr.lhs() * grad(expr.rhs()))
{
    return outer(grad(expr.lhs()), expr.rhs()) + expr.lhs() * grad(expr.rhs());
}

// Partial specialization for cross product
// ∇(A×B) = (∇A)×B - (∇B)×A
// For vectors A,B (order 1), result is a matrix (order 2)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const CrossProductExpression<_LhsExpr,_RhsExpr,1,1,ExpressionTraits<_LhsExpr>::Space,_RhsExpr::Space> expr)
-> decltype(cross(grad(expr.lhs()), expr.rhs()) - cross(grad(expr.rhs()), expr.lhs()))
{
    return (cross(grad(expr.lhs()), expr.rhs()) - cross(grad(expr.rhs()), expr.lhs()));
}

// Partial specialization for outer product
// ∇(A⊗B) = (∇A)⊗B + A⊗(∇B)
// For vectors A,B (order 1), result is a tensor (order 3)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const OuterProductExpression<_LhsExpr,_RhsExpr,1,1,ExpressionTraits<_LhsExpr>::Space,_RhsExpr::Space> expr)
-> decltype(outer(grad(expr.lhs()), expr.rhs()) + outer(expr.lhs(), grad(expr.rhs())))
{
    return (outer(grad(expr.lhs()), expr.rhs()) + outer(expr.lhs(), grad(expr.rhs())));
}

// Quotient rule for division by a scalar
// ∇(ψ/φ) = (φ∇ψ - ψ∇φ)/φ²
// For scalar functions ψ,φ (order 0), result is a vector (order 1)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const DivisionExpression<_LhsExpr,_RhsExpr,0,0,ExpressionTraits<_LhsExpr>::Space, ExpressionTraits<_RhsExpr>::Space> expr)
-> decltype((expr.rhs() * grad(expr.lhs()) - expr.lhs() * grad(expr.rhs())) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * grad(expr.lhs()) - expr.lhs() * grad(expr.rhs())) / (expr.rhs() * expr.rhs());
}

// Quotient rule for vector divided by scalar
// ∇(A/φ) = (φ∇A - A⊗∇φ)/φ²
// For vector A (order 1) and scalar φ (order 0), result is a matrix (order 2)
template <typename _LhsExpr, typename _RhsExpr>
auto grad(const DivisionExpression<_LhsExpr,_RhsExpr,1,0,ExpressionTraits<_LhsExpr>::Space, ExpressionTraits<_RhsExpr>::Space> expr)
-> decltype((expr.rhs() * grad(expr.lhs()) - outer(expr.lhs(), grad(expr.rhs()))) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * grad(expr.lhs()) - outer(expr.lhs(), grad(expr.rhs()))) / (expr.rhs() * expr.rhs());
}

// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Gradient of divergence is undefined
// ∇(∇•A) is undefined because divergence of a vector produces a scalar,
// but we cannot take gradient of a divergence operation directly
template <typename _E>
auto grad(const DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> void
{
    GISMO_ERROR("∇(∇•): gradient of divergence is not implemented");
}

// Gradient of curl for vectors is undefined in 3D
// ∇(∇×A) is undefined because curl of a vector produces a vector,
// but taking gradient of curl is not a standard vector calculus operation
template <typename _E>
auto grad(const CurlExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> void
{
    static_assert(false,"∇(∇×) is undefined: gradient of curl is not a valid vector calculus operation");
}

// Gradient of Laplacian is undefined
// ∇(∇²ψ) is undefined as a direct operation
template <typename _E>
auto grad(const LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> void
{
    static_assert(false,"∇(∇²) = ∇³ is not implemented as a separate expression type.");
}

// Variation of GradExpression where the inner expression is a SolutionObject
template <class T, size_t _Space, size_t _Order, typename _SpaceObject>
auto variation(const GradExpression<SolutionObject<T,_Space,_Order>,_Order,SpaceType::None,false> & expr,
          const _SpaceObject & space)
-> decltype(grad(variation(expr.expr(), space)))
{
    return grad(variation(expr.expr(), space));
}

// Partial specialization for addition
template <typename _LhsExpr, typename _RhsExpr, size_t _Order, size_t _Space, size_t _IsConstant, typename _SpaceObject>
auto variation(const GradExpression<AddExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>,_Order,_Space,_IsConstant> expr,
            const _SpaceObject & space)
-> decltype(grad(variation(expr.expr().lhs(), space)) + grad(variation(expr.expr().rhs(), space)))
{
    return grad(variation(expr.expr().lhs(), space)) + grad(variation(expr.expr().rhs(), space));
}

// Partial specialization for subtraction
template <typename _LhsExpr, typename _RhsExpr, size_t _Order, size_t _Space, size_t _IsConstant, typename _SpaceObject>
auto variation(const GradExpression<SubtractExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space>,_Order,_Space,_IsConstant> expr,
            const _SpaceObject & space)
-> decltype(grad(variation(expr.expr().lhs(), space)) - grad(variation(expr.expr().rhs(), space)))
{
    return grad(variation(expr.expr().lhs(), space)) - grad(variation(expr.expr().rhs(), space));
}

// Partial specialization for multiplication
template <typename _LhsExpr, typename _RhsExpr, size_t _Order, size_t _Space, size_t _IsConstant, typename _SpaceObject>
auto variation(const GradExpression<ProductExpression<_LhsExpr,_RhsExpr,
                                        _LhsExpr::Order,
                                        _RhsExpr::Order,
                                        ExpressionTraits<_LhsExpr>::Space,
                                        _RhsExpr::Space>,_Order,_Space,_IsConstant> expr,
            const _SpaceObject & space)
-> decltype(grad(variation(expr.expr().lhs(), space)) * expr.expr().rhs() + expr.expr().lhs() * grad(variation(expr.expr().rhs(), space)))
{
    return grad(variation(expr.expr().lhs(), space)) * expr.expr().rhs() + expr.expr().lhs() * grad(variation(expr.expr().rhs(), space));
}

// Partial specialization for division
template <typename _LhsExpr, typename _RhsExpr, size_t _Order, size_t _Space, size_t _IsConstant, typename _SpaceObject>
auto variation(const GradExpression<DivisionExpression<_LhsExpr,_RhsExpr,
                                        _LhsExpr::Order,
                                        _RhsExpr::Order,
                                        ExpressionTraits<_LhsExpr>::Space,
                                        ExpressionTraits<_RhsExpr>::Space>,_Order,_Space,_IsConstant> expr,
            const _SpaceObject & space) 
-> decltype(grad(variation(expr.expr().lhs(), space)) * expr.expr().rhs() - expr.expr().lhs() * grad(variation(expr.expr().rhs(), space)) / (expr.expr().rhs() * expr.expr().rhs()))
{
    return (grad(variation(expr.expr().lhs(), space)) * expr.expr().rhs() - expr.expr().lhs() * grad(variation(expr.expr().rhs(), space))) / (expr.expr().rhs() * expr.expr().rhs());
}

} // namespace Expr
} // namespace gismo