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

template <typename _E, size_t _Order, enum SpaceType _Space, size_t _IsConstant>
struct ExpressionTraits<DivExpression<_E, _Order, _Space, _IsConstant>>
{
    static_assert(_Order == 1,
                "DivExpression: Unsupported tensor order. Only vectors (order 1) are supported for divergence operations.");

    typedef _E ExprType; // Needed for UnaryOperator

    typedef typename ExpressionTraits<_E>::Scalar Scalar;
    static constexpr size_t Order = ExpressionTraits<_E>::Order-1; // Divergence results in a vector (order 1)
    static constexpr SpaceType Space = ExpressionTraits<_E>::Space;
    static constexpr size_t Deriv = ExpressionTraits<_E>::Deriv + 1; // Increment derivative order
    static constexpr bool IsConstant = ExpressionTraits<_E>::IsConstant;
};

// --- Unified DivExpression using enable_if for Space-aware eval ---
template <typename _E, size_t _Order, enum SpaceType _Space, size_t _IsConstant>
class DivExpression : public UnaryOperator<DivExpression<_E, _Order, _Space, _IsConstant>>
{
    static_assert(_Order == 1,
                  "DivExpression: Unsupported tensor order. Only vectors (order 1) are supported for divergence operations.");

    using Base = UnaryOperator<DivExpression<_E, _Order, _Space, _IsConstant>>;

private:
    std::array<size_t,Base::Order> sizes_;
    mutable gsMatrix<typename Base::Scalar> tmp;

public:
    typedef typename Base::Scalar Scalar;

    explicit DivExpression(const _E& e) : Base(e), tmp(1, 1)
    {
        // Assumes symmetry!
        for (short_t d=0; d!=ExpressionTraits<_E>::Order; d++)
        {
            GISMO_ENSURE(this->expr_.sizes()[d] == this->expr_.sizes()[0], "All sizes must be equal for the div operator");
        }

        for (short_t d=1; d!=ExpressionTraits<_E>::Order; d++)
            sizes_[d-1] = this->expr_.sizes()[d];
    }

    const std::array<size_t, Base::Order> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const
    {
        return this->expr_.domainDim();
    }

    // Eval for Space = None (constant)
    template <enum SpaceType S = _Space, size_t C = _IsConstant>
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
    template <enum SpaceType S = _Space, size_t C = _IsConstant>
    typename std::enable_if<S == SpaceType::None && C == 0, ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Divergence computation: sum of diagonal derivatives
        // ∇‧V = ∂V_x/∂x + ∂V_y/∂y + ∂V_z/∂z
        const index_t dim = this->expr_.domainDim();
        gsAsConstMatrix<Scalar, Dynamic, Dynamic> deriv = 
            this->expr_.data().values[1].reshapeCol(k, dim, this->expr_.data().values[0].rows());
        
        Scalar divVal = 0;
        for (index_t d = 0; d < dim; ++d)
            divVal += deriv(d, d);
        
        tmp.resize(1, 1);
        tmp(0, 0) = divVal;
        
        ExpressionResult<typename Base::Scalar> result(1, 1);
        result(0, 0) = tmp;
        return result;
    }

    // Eval for Space = Test or Trial (constant)
    template <enum SpaceType S = _Space, size_t C = _IsConstant>
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
    template <enum SpaceType S = _Space, size_t C = _IsConstant>
    typename std::enable_if<(S == SpaceType::Test || S == SpaceType::Trial) && C == 0, ExpressionResult<typename Base::Scalar>>::type
    eval(const index_t k) const
    {
        // Divergence for basis functions (vector space)
        // For a vector space with `dim` components, each component uses the same scalar basis.
        // The derivative data `values[1]` contains derivatives of the scalar basis functions.
        // Shape: (domainDim * numActive, numPts) where numActive = number of scalar basis functions
        //
        // For a vector-valued basis function [0, ..., N_i, ..., 0]^T with N_i in component c:
        //   div([0, ..., N_i, ..., 0]) = ∂N_i/∂x_c
        //
        // The divergence is scalar-valued, so we return one value per scalar basis function.
        // The value is the sum of all diagonal derivatives: div(N_i * e_c) for c = 0..dim-1
        // which equals ∑_c ∂N_i/∂x_c = ∑_d ∂N_i/∂x_d (trace of the gradient)
        
        const index_t numActive = this->expr_.data().values[0].rows();
        const index_t dim = this->expr_.domainDim();
        
        // Reshape: values[1] has shape (domainDim * numActive, numPts)
        // We want to access it as (domainDim, numActive) for a single point
        gsAsConstMatrix<Scalar, Dynamic, Dynamic> deriv = 
            this->expr_.data().values[1].reshapeCol(k, dim, numActive);
        
        ExpressionResult<typename Base::Scalar> result(
            S == SpaceType::Test ? numActive : 1,
            S == SpaceType::Trial ? numActive : 1
        );
        
        // Compute divergence for each scalar basis function
        // div(N_i) = ∂N_i/∂x_0 + ∂N_i/∂x_1 + ... (sum over all spatial dimensions)
        for (index_t i = 0; i < numActive; ++i)
        {
            Scalar divVal = 0;
            for (index_t d = 0; d < dim; ++d)
                divVal += deriv(d, i);  // deriv(d, i) = ∂N_i/∂x_d
            
            tmp.resize(1, 1);
            tmp(0, 0) = divVal;
            
            if (S == SpaceType::Test)
                result(i, 0) = tmp;
            else // Trial
                result(0, i) = tmp;
        }
        
        return result;
    }

    void parse(ExpressionHelper<typename Base::Scalar> & helper) const
    {
        // Set derivative order first so child expression knows what to request
        this->expr_.setDerivative(Base::Deriv);
        // Parse the underlying expression
        this->expr_.parse(helper);
        // Additionally ensure derivative flag is set
        this->expr_.data().flags |= NEED_DERIV;
    }

    void print(std::ostream & os) const
    {
        os<<"\u2207\u2027("<<this->expr_<<")";
    }

};

// Generic factory function for easy creation
template <typename _E>
DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>
div(const _E& expr)
{
    return DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr);
}

// Partial specialization for addition
// ∇•(X + Y) = ∇•X + ∇•Y
template <typename _LhsExpr, typename _RhsExpr>
auto div(const AddExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> expr)
-> decltype(div(expr.lhs()) + div(expr.rhs()))
{
    return (div(expr.lhs()) + div(expr.rhs()));
}

// Partial specialization for subtraction
// ∇•(X - Y) = ∇•X - ∇•Y
template <typename _LhsExpr, typename _RhsExpr>
auto div(const SubtractExpression<_LhsExpr,_RhsExpr,_LhsExpr::Order,_RhsExpr::Order,_LhsExpr::Space,_RhsExpr::Space> expr)
-> decltype(div(expr.lhs()) - div(expr.rhs()))
{
    return (div(expr.lhs()) - div(expr.rhs()));
}

// Partial specialization for multiplication of a scalar and a vector
// ∇•(fV) = ∇f • V + f ∇•V
template <typename _LhsExpr, typename _RhsExpr, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
auto div(const ProductExpression<_LhsExpr,_RhsExpr,0,1,_LhsSpace,_RhsSpace>& expr)
-> decltype(dot(grad(expr.lhs()), expr.rhs()) + expr.lhs() * div(expr.rhs()))
{
    return dot(grad(expr.lhs()), expr.rhs()) + expr.lhs() * div(expr.rhs());
}

// Partial specialization for cross product
// ∇•(A×B) = (∇×A)•B - A•(∇×B) = B•(∇×A) - A•(∇×B)
// For vectors A,B (order 1), result is a scalar (order 0)
template <typename _LhsExpr, typename _RhsExpr>
auto div(const CrossProductExpression<_LhsExpr,_RhsExpr,1,1,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(dot(expr.rhs(), curl(expr.lhs())) - dot(expr.lhs(), curl(expr.rhs())))
{
    return dot(expr.rhs(), curl(expr.lhs())) - dot(expr.lhs(), curl(expr.rhs()));
}

// Partial specialization for outer product
// ∇•(A⊗B) = (∇•A)B + (A•∇)B
// For vectors A,B (order 1), result is a vector (order 1)
template <typename _LhsExpr, typename _RhsExpr>
auto div(const OuterProductExpression<_LhsExpr,_RhsExpr,1,1,_LhsExpr::Space,_RhsExpr::Space>& expr)
-> decltype(div(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs()))
{
    GISMO_UNUSED(expr);
    return div(expr.lhs()) * expr.rhs() + expr.lhs() * grad(expr.rhs());
}

// Partial specialization for gradient (second derivative identity)
// ∇•(∇ψ) = ∇²ψ (Laplacian)
// For scalar field ψ (order 0), result is a scalar (order 0)
template <typename _E>
auto div(const GradExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>& expr)
-> LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>
{
    return LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant>(expr.expr());
}

// === UNDEFINED OPERATIONS ===
// These operations are mathematically undefined and will produce compile-time errors

// Divergence of gradient produces Laplacian, not undefined, but divergence of divergence is undefined
// ∇•(∇•A) is undefined because divergence of a vector produces a scalar,
// and divergence of a scalar is not defined
template <typename _E>
auto div(const DivExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr) -> void
{
    GISMO_ERROR("∇•(∇•) is undefined: divergence of divergence is not defined (scalar has no divergence)");
}

// Divergence of curl is always zero (this is a valid identity)
// ∇•(∇×A) = 0, but for completeness, we could add this as a zero expression
// This is actually defined and equals zero, so we should implement it properly
template <typename _E>
auto div(const CurlExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr)
-> ConstantObject<typename ExpressionTraits<_E>::Scalar, 0>
{
    GISMO_UNUSED(expr);
    return ConstantObject<typename ExpressionTraits<_E>::Scalar, 0>(std::array<size_t, 0>{},"0");  // Divergence of curl is always zero scalar
}

// Divergence of Laplacian is undefined
// ∇•(∇²ψ) is undefined because Laplacian of a scalar produces a scalar,
// and divergence of a scalar is not defined
template <typename _E>
auto div(const LaplExpression<_E, ExpressionTraits<_E>::Order, ExpressionTraits<_E>::Space, ExpressionTraits<_E>::IsConstant> expr) -> void
{
    GISMO_ERROR("∇•(∇²) is undefined: divergence of Laplacian is not defined (scalar has no divergence)");
}

// Partial specialization for division of a vector by a scalar
// ∇•(A/φ) = (φ∇•A - ∇φ•A)/φ²
template <typename _LhsExpr, typename _RhsExpr, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace>
auto div(const DivisionExpression<_LhsExpr,_RhsExpr,1,0,_LhsSpace,_RhsSpace>& expr)
-> decltype((expr.rhs() * div(expr.lhs()) - dot(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs()))
{
    return (expr.rhs() * div(expr.lhs()) - dot(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs());
}

// --- Partial Specialization 2: Divergence of a VariableObject ---



// // --- Partial Specialization 1: Addition of two expressions of the SAME ORDER (X + X) ---
// template <typename _LhsExpr, typename _RhsExpr>
// class AddExpression<ArrayExpression<LhsExpr>, RhsExpr,
//     typename std::enable_if<0 == (RhsExpr::Order)>::type // Simplified condition
// > : public BaseObject<LhsExpr,
//                           typename LhsExpr::Scalar,
//                           LhsExpr::Order,
//                           LhsExpr::IsConstant && RhsExpr::IsConstant,
//                           0> // Use LhsExpr's Scalar and Order directly
// {
// public:
// // Scalar and Order are directly from LhsExpr/RhsExpr
//     typedef typename LhsExpr::Scalar Scalar;
//     static constexpr int Order = LhsExpr::Order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<RhsExpr, typename LhsExpr::Scalar, Order, LhsExpr::IsConstant && RhsExpr::IsConstant, 0>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
//           lhs_this->expr_(lhs),
//           rhs_this->expr_(rhs)
//     {
//     }

//     gsMatrix<Scalar> eval(const index_t k) const
//     {
//         gsMatrix<Scalar> lhs_val = lhs_this->expr_.eval(k);
//         gsMatrix<Scalar> rhs_val = rhs_this->expr_.eval(k);
//         lhs_val.array() += rhs_val; // Element-wise addition
//         return lhs_val; // Return the modified lhs_val
//     }

// private:
//     const LhsExpr& lhs_this->expr_;
//     const RhsExpr& rhs_this->expr_;
// };



// // --- Partial Specialization 2: Scalar (Order 0) + Higher Order (Order N > 0) ---
// template <typename _LhsExpr, typename _RhsExpr>
// class AddExpression<LhsExpr, RhsExpr,
//     typename std::enable_if<(LhsExpr::Order == 0) && (RhsExpr::Order > 0)>::type // Simplified condition
// > : public BaseObject<RhsExpr, typename LhsExpr::Scalar, RhsExpr::Order> // Base on RhsExpr for Scalar and Order
// {
// public:
//     typedef typename RhsExpr::Scalar Scalar;
//     static constexpr int Order = RhsExpr::Order;

// public:
//     AddExpression(const LhsExpr& lhs, const RhsExpr& rhs)
//         : BaseObject<RhsExpr, typename LhsExpr::Scalar, Order>(rhs.sizes()), // Pass RhsExpr as Derived to BaseObject
//           lhs_this->expr_(lhs),
//           rhs_this->expr_(rhs)
//     {}

//     gsMatrix<Scalar> eval(const index_t k) const
//     {
//         gsMatrix<Scalar> lhs_val = lhs_this->expr_.eval(k);
//         gsMatrix<Scalar> rhs_val = rhs_this->expr_.eval(k);
//         rhs_val.array() += lhs_val.value(); // Add scalar to each element of lhs_val
//         return rhs_val;
//     }

// private:
//     const LhsExpr& lhs_this->expr_;
//     const RhsExpr& rhs_this->expr_;
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

// Variation of DivExpression where the inner expression is a SolutionObject
template <class T, enum SpaceType _Space, size_t _Order, typename _SpaceObject>
auto variation(const DivExpression<SolutionObject<T,_Space,_Order>,_Order,SpaceType::None,false> & expr,
          const _SpaceObject & space)
-> decltype(div(variation(expr.expr(), space)))
{
    return div(variation(expr.expr(), space));
}

// Partial specialization for addition
template <typename _LhsExpr, typename _RhsExpr, size_t _Order, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace, typename _SpaceObject>
auto variation(const DivExpression<AddExpression<_LhsExpr,_RhsExpr,_Order,_Order,_LhsSpace,_RhsSpace>,_Order, static_cast<SpaceType>(_LhsSpace | _RhsSpace), (_LhsExpr::IsConstant && _RhsExpr::IsConstant)> & expr,
          const _SpaceObject & space)
-> decltype(variation(div(expr.lhs()), space) + variation(div(expr.rhs()), space))
{
    return variation(div(expr.lhs()), space) + variation(div(expr.rhs()), space);   
}

// Partial specialization for subtraction
template <typename _LhsExpr, typename _RhsExpr, size_t _Order, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace, typename _SpaceObject>
auto variation(const DivExpression<SubtractExpression<_LhsExpr,_RhsExpr,_Order,_Order,_LhsSpace,_RhsSpace>,_Order, static_cast<SpaceType>(_LhsSpace | _RhsSpace), (_LhsExpr::IsConstant && _RhsExpr::IsConstant)> & expr,
          const _SpaceObject & space)
-> decltype(variation(div(expr.lhs()), space) - variation(div(expr.rhs()), space))
{
    return variation(div(expr.lhs()), space) - variation(div(expr.rhs()), space);
}

// Partial specialization for multiplication 
template <typename _LhsExpr, typename _RhsExpr, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace, typename _SpaceObject>
auto variation(const DivExpression<ProductExpression<_LhsExpr,_RhsExpr,0,1,_LhsSpace,_RhsSpace>,1, static_cast<SpaceType>(_LhsSpace | _RhsSpace), (_LhsExpr::IsConstant && _RhsExpr::IsConstant)> & expr,
          const _SpaceObject & space)
-> decltype(variation(dot(grad(expr.lhs()), expr.rhs()), space) + variation(expr.lhs() * div(expr.rhs()), space))
{
    return variation(dot(grad(expr.lhs()), expr.rhs()), space) + variation(expr.lhs() * div(expr.rhs()), space);    
}

// Partial specialization for cross product
template <typename _LhsExpr, typename _RhsExpr, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace, typename _SpaceObject>
auto variation(const DivExpression<CrossProductExpression<_LhsExpr,_RhsExpr,1,1,_LhsSpace,_RhsSpace>,0, static_cast<SpaceType>(_LhsSpace | _RhsSpace), (_LhsExpr::IsConstant && _RhsExpr::IsConstant)> & expr,
          const _SpaceObject & space)
-> decltype(variation(dot(expr.rhs(), curl(expr.lhs())), space) - variation(dot(expr.lhs(), curl(expr.rhs())), space))
{
    return variation(dot(expr.rhs(), curl(expr.lhs())), space) - variation(dot(expr.lhs(), curl(expr.rhs())), space);    
}

// Partial specialization for outer product
template <typename _LhsExpr, typename _RhsExpr, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace, typename _SpaceObject>
auto variation(const DivExpression<OuterProductExpression<_LhsExpr,_RhsExpr,1,1,_LhsSpace,_RhsSpace>,1, static_cast<SpaceType>(_LhsSpace | _RhsSpace), (_LhsExpr::IsConstant && _RhsExpr::IsConstant)> & expr,
          const _SpaceObject & space)
-> decltype(variation(div(expr.lhs()) * expr.rhs(), space) + variation(expr.lhs() * grad(expr.rhs()), space))
{
    return variation(div(expr.lhs()) * expr.rhs(), space) + variation(expr.lhs() * grad(expr.rhs()), space);    
}

// Partial specialization for gradient
template <typename _E, enum SpaceType _Space, size_t _IsConstant, typename _SpaceObject>
auto variation(const DivExpression<GradExpression<_E, ExpressionTraits<_E>::Order, _Space, _IsConstant>, ExpressionTraits<_E>::Order - 1, _Space, _IsConstant> & expr,
          const _SpaceObject & space)
-> decltype(variation(LaplExpression<_E, ExpressionTraits<_E>::Order, _Space, _IsConstant>(expr.expr().expr()), space))
{
    return variation(LaplExpression<_E, ExpressionTraits<_E>::Order, _Space, _IsConstant>(expr.expr().expr()), space);
}

// Partial specialization for division
template <typename _LhsExpr, typename _RhsExpr, enum SpaceType _LhsSpace, enum SpaceType _RhsSpace, typename _SpaceObject>
auto variation(const DivExpression<DivisionExpression<_LhsExpr,_RhsExpr,1,0,_LhsSpace,_RhsSpace>,1, static_cast<SpaceType>(_LhsSpace | _RhsSpace), (_LhsExpr::IsConstant && _RhsExpr::IsConstant)> & expr,
          const _SpaceObject & space)
-> decltype(variation((expr.rhs() * div(expr.lhs()) - dot(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs()), space))
{
    return variation((expr.rhs() * div(expr.lhs()) - dot(grad(expr.rhs()), expr.lhs())) / (expr.rhs() * expr.rhs()), space);
}

}//namespace Expr
}//namespace gismo