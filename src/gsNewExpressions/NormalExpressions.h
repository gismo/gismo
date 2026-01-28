/** @file NormalExpressions.h

    @brief Expressions for surface normals, boundary normals, and tangent vectors

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

template <class T> class SurfaceNormalExpression;
template <class T> class BoundaryNormalExpression;
template <class T> class BoundaryTangentExpression;

// ============================================================
// SurfaceNormalExpression - Out-of-plane surface normal
// ============================================================

template <class T>
struct ExpressionTraits<SurfaceNormalExpression<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 1; // Vector-valued
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true;
};

/**
 * @brief Expression for the out-of-plane surface normal of a geometry map
 * 
 * This is valid for codimension-1 geometries (e.g., surfaces in 3D space).
 * The surface normal is perpendicular to the tangent plane.
 * 
 * @tparam T The scalar type
 */
template<class T>
class SurfaceNormalExpression : public UnaryOperator<SurfaceNormalExpression<T>>
{
    using Base = UnaryOperator<SurfaceNormalExpression<T>>;
    
    mutable gsMatrix<T> tmp;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<SurfaceNormalExpression<T>>::ExprType ExprType;

    SurfaceNormalExpression(const GeometryMap<T> & geomMap)
    : Base(geomMap)
    {
        GISMO_ENSURE(geomMap.source().domainDim() + 1 == geomMap.source().targetDim(),
                     "Surface normal requires codimension 1 (e.g., 2D surface in 3D)");
    }

    const GeometryMap<T> & geomMap() const 
    { 
        return static_cast<const GeometryMap<T>&>(this->expr());
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        this->expr().parse(helper);
        // TODO: Set flag for computing normals
        // helper.setFlag(NEED_NORMAL);
    }

    void print(std::ostream & os) const
    {
        os << "sn(";
        this->expr().print(os);
        os << ")";
    }

    std::array<size_t, 1> sizes() const
    {
        return {static_cast<size_t>(geomMap().targetDim())};
    }

    size_t domainDim() const { return this->expr().domainDim(); }

    /// Evaluate the surface normal at quadrature point k
    ExpressionResult<T> eval(const index_t k) const
    {
        // Get the normals from gsMapData
        // gsMapData::normals contains the surface normals
        const auto & normals = geomMap().data().normals;
        
        tmp.resize(geomMap().targetDim(), 1);
        tmp = normals.col(k);
        
        return tmp;
    }

    std::string label() const
    {
        return "sn(G)";
    }
};

// ============================================================
// BoundaryNormalExpression - Outward boundary normal
// ============================================================

template <class T>
struct ExpressionTraits<BoundaryNormalExpression<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 1; // Vector-valued
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true;
};

/**
 * @brief Expression for the outward-pointing boundary normal
 * 
 * This is valid only at the boundaries of a geometric patch.
 * The boundary normal points outward from the domain.
 * 
 * @tparam T The scalar type
 */
template<class T>
class BoundaryNormalExpression : public UnaryOperator<BoundaryNormalExpression<T>>
{
    using Base = UnaryOperator<BoundaryNormalExpression<T>>;
    
    mutable gsMatrix<T> tmp;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<BoundaryNormalExpression<T>>::ExprType ExprType;

    explicit BoundaryNormalExpression(const GeometryMap<T> & geomMap)
    : Base(geomMap)
    {
    }

    const GeometryMap<T> & geomMap() const 
    { 
        return static_cast<const GeometryMap<T>&>(this->expr());
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        this->expr().parse(helper);
        // TODO: Set flag for computing outer normals
        // helper.setFlag(NEED_OUTER_NORMAL);
    }

    void print(std::ostream & os) const
    {
        os << "nv(";
        this->expr().print(os);
        os << ")";
    }

    std::array<size_t, 1> sizes() const
    {
        return {static_cast<size_t>(geomMap().targetDim())};
    }

    size_t domainDim() const { return this->expr().domainDim(); }

    /// Evaluate the boundary normal at quadrature point k
    ExpressionResult<T> eval(const index_t k) const
    {
        // Get the outer normals from gsMapData
        // gsMapData::outNormals contains the boundary normals
        const auto & outNormals = geomMap().data().outNormals;
        
        tmp.resize(geomMap().targetDim(), 1);
        tmp = outNormals.col(k);
        
        return tmp;
    }

    std::string label() const
    {
        return "nv(G)";
    }
};

// ============================================================
// BoundaryTangentExpression - Boundary tangent vector
// ============================================================

template <class T>
struct ExpressionTraits<BoundaryTangentExpression<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 1; // Vector-valued
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true;
};

/**
 * @brief Expression for the boundary tangent vector
 * 
 * This is valid only at the boundaries of a geometric patch.
 * For 2D: tangent is perpendicular to the boundary normal
 * For 3D: tangent is computed as sn × nv (surface normal × boundary normal)
 * 
 * @tparam T The scalar type
 */
template<class T>
class BoundaryTangentExpression : public UnaryOperator<BoundaryTangentExpression<T>>
{
    using Base = UnaryOperator<BoundaryTangentExpression<T>>;
    
    mutable gsMatrix<T> tmp;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<BoundaryTangentExpression<T>>::ExprType ExprType;

    BoundaryTangentExpression(const GeometryMap<T> & geomMap)
    : Base(geomMap)
    {
    }

    const GeometryMap<T> & geomMap() const 
    { 
        return static_cast<const GeometryMap<T>&>(this->expr());
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        this->expr().parse(helper);
        // TODO: Set flags for computing both normals
        // helper.setFlag(NEED_NORMAL | NEED_OUTER_NORMAL);
    }

    void print(std::ostream & os) const
    {
        os << "tv(";
        this->expr().print(os);
        os << ")";
    }

    std::array<size_t, 1> sizes() const
    {
        return {static_cast<size_t>(geomMap().targetDim())};
    }

    size_t domainDim() const { return this->expr().domainDim(); }

    /// Evaluate the boundary tangent at quadrature point k
    ExpressionResult<T> eval(const index_t k) const
    {
        const auto & outNormals = geomMap().data().outNormals;
        
        tmp.resize(geomMap().targetDim(), 1);
        
        if (geomMap().targetDim() == 2)
        {
            // 2D: Rotate the outer normal by 90 degrees
            // tangent = (-ny, nx) where normal = (nx, ny)
            tmp(0, 0) = -outNormals(1, k);
            tmp(1, 0) = outNormals(0, k);
        }
        else if (geomMap().targetDim() == 3)
        {
            // 3D: Cross product of surface normal and boundary normal
            const auto & normals = geomMap().data().normals;
            
            // tangent = sn × nv
            tmp(0, 0) = normals(1, k) * outNormals(2, k) - normals(2, k) * outNormals(1, k);
            tmp(1, 0) = normals(2, k) * outNormals(0, k) - normals(0, k) * outNormals(2, k);
            tmp(2, 0) = normals(0, k) * outNormals(1, k) - normals(1, k) * outNormals(0, k);
        }
        else
        {
            GISMO_ERROR("Boundary tangent not implemented for dimension " << geomMap().targetDim());
        }
        
        return tmp;
    }

    std::string label() const
    {
        return "tv(G)";
    }
};

// ============================================================
// Helper functions to create normal/tangent expressions
// ============================================================

/**
 * @brief Create a surface normal expression
 * @param G The geometry map
 * @return SurfaceNormalExpression
 */
template<class T>
SurfaceNormalExpression<T> sn(const GeometryMap<T> & G)
{
    return SurfaceNormalExpression<T>(G);
}

/**
 * @brief Create a boundary normal expression
 * @param G The geometry map
 * @return BoundaryNormalExpression
 */
template<class T>
BoundaryNormalExpression<T> nv(const GeometryMap<T> & G)
{
    return BoundaryNormalExpression<T>(G);
}

/**
 * @brief Create a boundary tangent expression
 * @param G The geometry map
 * @return BoundaryTangentExpression
 */
template<class T>
BoundaryTangentExpression<T> tv(const GeometryMap<T> & G)
{
    return BoundaryTangentExpression<T>(G);
}

// Normalized versions (commonly used)

/**
 * @brief Create a normalized boundary normal expression
 * @param G The geometry map
 * @return Normalized boundary normal
 */
template<class T>
auto unv(const GeometryMap<T> & G) -> decltype(nv(G))
{
    // TODO: Add normalization when NormExpression is implemented
    return nv(G);
}

/**
 * @brief Create a normalized surface normal expression
 * @param G The geometry map
 * @return Normalized surface normal
 */
template<class T>
auto usn(const GeometryMap<T> & G) -> decltype(sn(G))
{
    // TODO: Add normalization when NormExpression is implemented
    return sn(G);
}

} // namespace Expr
} // namespace gismo
