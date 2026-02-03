/** @file JacobianExpressions.h

    @brief Jacobian and inverse Jacobian expressions for geometry maps

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

// ============================================================
// Forward declarations
// ============================================================

template <class T> class JacobianExpression;
template <class T> class InverseJacobianExpression;

// ============================================================
// JacobianExpression - Jacobian matrix of geometry map
// ============================================================

template <class T>
struct ExpressionTraits<JacobianExpression<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 2; // Matrix-valued
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true;
};

/**
 * @brief Expression for the Jacobian matrix of a geometry map
 * 
 * The Jacobian DG is the matrix of partial derivatives:
 * DG[i,j] = ∂G_i/∂ξ_j
 * 
 * where G: Ω̂ → Ω maps from parameter domain to physical domain.
 * 
 * Size: targetDim × domainDim
 * 
 * @tparam T The scalar type
 */
template<class T>
class JacobianExpression : public UnaryOperator<JacobianExpression<T>>
{
    using Base = UnaryOperator<JacobianExpression<T>>;
    
    mutable gsMatrix<T> tmp;
    std::array<size_t, 2> sizes_;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<JacobianExpression<T>>::ExprType ExprType;

    JacobianExpression(const GeometryMap<T> & geomMap)
    : Base(geomMap)
    {
        sizes_[0] = geomMap.targetDim();
        sizes_[1] = geomMap.domainDim();
    }

    const GeometryMap<T> & geomMap() const 
    { 
        return static_cast<const GeometryMap<T>&>(this->expr());
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        // Set derivative order first
        this->expr().setDerivative(Base::Deriv);
        this->expr().parse(helper);
        // Jacobian is already computed in gsMapData
    }

    void print(std::ostream & os) const
    {
        os << "jac(";
        this->expr().print(os);
        os << ")";
    }

    const std::array<size_t, 2> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const { return this->expr().domainDim(); }

    /// Evaluate the Jacobian at quadrature point k
    /// Returns targetDim × domainDim matrix
    ExpressionResult<T> eval(const index_t k) const
    {
        return ExpressionResult<T>(geomMap().jacobian(k));
    }

    std::string label() const
    {
        return "jac(G)";
    }
};

// ============================================================
// InverseJacobianExpression - Inverse/Pseudo-inverse Jacobian
// ============================================================

template <class T>
struct ExpressionTraits<InverseJacobianExpression<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 2; // Matrix-valued
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true;
};

/**
 * @brief Expression for the inverse (or pseudo-inverse) Jacobian
 * 
 * For square Jacobians: (DG)^{-1}
 * For non-square: generalized inverse (DG^T DG)^{-1} DG^T
 * 
 * This is used for pullbacks: grad_physical(u) = grad_param(u) * (DG)^{-1}
 * 
 * Size: domainDim × targetDim
 * 
 * @tparam T The scalar type
 */
template<class T>
class InverseJacobianExpression : public UnaryOperator<InverseJacobianExpression<T>>
{
    using Base = UnaryOperator<InverseJacobianExpression<T>>;
    
    mutable gsMatrix<T> tmp;
    std::array<size_t, 2> sizes_;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<InverseJacobianExpression<T>>::ExprType ExprType;

    InverseJacobianExpression(const GeometryMap<T> & geomMap)
    : Base(geomMap)
    {
        sizes_[0] = geomMap.domainDim();
        sizes_[1] = geomMap.targetDim();
    }

    const GeometryMap<T> & geomMap() const 
    { 
        return static_cast<const GeometryMap<T>&>(this->expr());
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        // Set derivative order first
        this->expr().setDerivative(Base::Deriv);
        this->expr().parse(helper);
        // TODO: Set flag for computing inverse Jacobian
        // helper.setFlag(NEED_GRAD_TRANSFORM);
    }

    void print(std::ostream & os) const
    {
        os << "jacInv(";
        this->expr().print(os);
        os << ")";
    }

    const std::array<size_t, 2> & sizes() const
    {
        return sizes_;
    }

    size_t domainDim() const { return this->expr().domainDim(); }

    /// Evaluate the inverse Jacobian at quadrature point k
    /// Returns domainDim × targetDim matrix
    ExpressionResult<T> eval(const index_t k) const
    {
        // gsMapData stores jacInvTr (transpose of inverse Jacobian)
        const auto & jacInvTr = geomMap().data().jacInvTr;
        
        tmp.resize(sizes_[0], sizes_[1]);
        tmp = jacInvTr.reshapeCol(k, sizes_[1], sizes_[0]).transpose();
        
        return ExpressionResult<T>(tmp);
    }

    /**
     * @brief Generalized inverse for non-square Jacobians
     * Returns (J^T J)^{-1} J^T
     */
    gsMatrix<T> generalizedInverse(const index_t k) const
    {
        auto jac = geomMap().jacobian(k);
        
        if (jac.rows() == jac.cols())
        {
            // Square: regular inverse
            return jac.inverse();
        }
        else
        {
            // Non-square: pseudo-inverse
            return (jac.transpose() * jac).inverse() * jac.transpose();
        }
    }

    std::string label() const
    {
        return "jacInv(G)";
    }
};

// ============================================================
// Helper functions
// ============================================================

/**
 * @brief Create a Jacobian expression for a GeometryMap
 * @param G The geometry map
 * @return JacobianExpression
 */
template<class T>
JacobianExpression<T> jac(const GeometryMap<T> & G)
{
    return JacobianExpression<T>(G);
}

/**
 * @brief Create an inverse Jacobian expression for a GeometryMap
 * @param G The geometry map
 * @return InverseJacobianExpression
 */
template<class T>
InverseJacobianExpression<T> jacInv(const GeometryMap<T> & G)
{
    return InverseJacobianExpression<T>(G);
}

// ============================================================
// Pullback operations: grad(u, G) transforms to physical domain
// ============================================================

/**
 * @brief Physical gradient - gradient with respect to physical coordinates
 * 
 * For a function u defined on the parameter domain:
 * grad_physical(u) = grad_param(u) * (DG)^{-1}
 * 
 * This is the "igrad" operation in the old expressions.
 * 
 * @param u The expression to take gradient of
 * @param G The geometry map for pullback
 * @return Gradient in physical coordinates
 */
template<typename E>
auto igrad(const E & u, const GeometryMap<typename E::Scalar> & G)
    -> decltype(grad(u) * jacInv(G))
{
    return grad(u) * jacInv(G);
}

/**
 * @brief Physical divergence - divergence with respect to physical coordinates
 * 
 * For a vector field u:
 * div_physical(u) = trace(ijac(u, G))
 * 
 * where ijac is the physical Jacobian.
 * 
 * @param u The expression to take divergence of
 * @param G The geometry map for pullback
 * @return Divergence in physical coordinates
 */
template<typename E>
auto idiv(const E & u, const GeometryMap<typename E::Scalar> & G)
    -> decltype(ijac(u, G))  // TODO: Add .trace() when implemented
{
    // TODO: Implement trace operation
    // return ijac(u, G).trace();
    return ijac(u, G);
}

/**
 * @brief Physical Jacobian of u with respect to physical coordinates
 * 
 * For a function u: Ω̂ → ℝ^m and geometry G: Ω̂ → Ω:
 * ijac(u) = jac(u) * (DG)^{-1}
 * 
 * This transforms derivatives from parameter space to physical space.
 * 
 * @param u The expression
 * @param G The geometry map for pullback
 * @return Physical Jacobian
 */
template<typename E>
auto ijac(const E & u, const GeometryMap<typename E::Scalar> & G)
    -> decltype(grad(u) * jacInv(G))
{
    // TODO: This should work for any expression with a Jacobian
    // For now, return grad * jacInv (works for scalar functions)
    // Full implementation needs jac(u) to work for all expression types
    return grad(u) * jacInv(G);
}

// Note: For grad(G), we want it to return jac(G) as an alias
// This is handled by specializing GradExpression for GeometryMap

} // namespace Expr
} // namespace gismo
