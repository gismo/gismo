/** @file GeometryMap.h

    @brief Expression for a geometry map

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{
namespace Expr
{

// Forward declaration
template <class T> class GeometryMap;

// ExpressionTraits for GeometryMap
template <class T>
struct ExpressionTraits<GeometryMap<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 1; // Geometry map is vector-valued (returns coordinates)
    static constexpr SpaceType Space = SpaceType::None; // Geometry is not a space function
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true; // Geometry doesn't change during assembly
};

/**
 * @brief Expression for a geometry map
 * 
 * This represents the parametrization G: Omega_hat -> Omega
 * where Omega_hat is the parameter domain and Omega is the physical domain.
 * 
 * @tparam T The scalar type
 */
template<class T>
class GeometryMap : public BaseObject<GeometryMap<T>>
{
    using Base = BaseObject<GeometryMap<T>>;
    
protected:
    const gsFunctionSet<T> * m_source;  ///< The geometry (typically gsMultiPatch)
    const gsMapData<T> * m_data;        ///< Evaluation data (Jacobians, etc.)
    
    mutable gsMatrix<T> tmp;  ///< Temporary storage for evaluation
    
    bool m_isAcross;  ///< true when evaluating across an interface

public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    /**
     * @brief Constructor
     * @param source The geometry function set (e.g., gsMultiPatch)
     */
    GeometryMap(const gsFunctionSet<T> & source)
    : Base(source.domainDim(), std::array<size_t,1>{static_cast<size_t>(source.targetDim())}, "G"),
      m_source(&source), m_data(nullptr), m_isAcross(false)
    {
    }

    /**
     * @brief Sets the evaluation data
     * @param data Pointer to gsMapData containing Jacobians, etc.
     */
    void setData(const gsMapData<T> * data)
    {
        m_data = data;
    }

    /// Returns the geometry function source
    const gsFunctionSet<T> & source() const 
    { 
        GISMO_ASSERT(m_source != nullptr, "GeometryMap: source is null");
        return *m_source;
    }

    /// Returns the evaluation data
    const gsMapData<T> & data() const
    {
        GISMO_ASSERT(m_data != nullptr, "GeometryMap: data is null");
        return *m_data;
    }

    /// Target dimension (physical space dimension)
    index_t targetDim() const { return m_source->targetDim(); }
    
    /// Domain dimension (parameter space dimension)
    size_t domainDim() const { return m_source->domainDim(); }

    std::array<size_t, 1> sizes() const
    {
        return {static_cast<size_t>(targetDim())};
    }

    bool isAcross() const { return m_isAcross; }

    /// Create a copy for evaluation on the right side of an interface
    GeometryMap right() const
    {
        GeometryMap gm(*m_source);
        gm.m_isAcross = true;
        gm.m_data = m_data;
        return gm;
    }

    /// Create a copy for evaluation on the left side of an interface
    GeometryMap left() const
    {
        GeometryMap gm(*m_source);
        gm.m_isAcross = false;
        gm.m_data = m_data;
        return gm;
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        // Register the geometry map with the helper
        helper.parse(*this);
    }

    void print(std::ostream & os) const
    {
        os << "G";
        if (m_isAcross) os << "_across";
    }

    /// Evaluate the geometry map at quadrature point k
    /// Returns the physical coordinates as a column vector
    ExpressionResult<T> eval(const index_t k) const
    {
        GISMO_ASSERT(m_data != nullptr, "GeometryMap: evaluation data not set");
        
        // Get the evaluated values from gsMapData
        // values(0) contains the physical coordinates
        const auto & values = m_data->values[0];
        
        tmp.resize(targetDim(), 1);
        tmp = values.col(k);
        
        return tmp;
    }

    /// Return the Jacobian matrix at quadrature point k
    gsMatrix<T> jacobian(const index_t k) const
    {
        GISMO_ASSERT(m_data != nullptr, "GeometryMap: evaluation data not set");
        
        const auto & jacs = m_data->values[1]; // Jacobian data
        
        gsMatrix<T> jac(targetDim(), domainDim());
        jac = jacs.reshapeCol(k, targetDim(), domainDim());
        
        return jac;
    }

    /// Return the Jacobian determinant at quadrature point k
    T jacobianDet(const index_t k) const
    {
        GISMO_ASSERT(m_data != nullptr, "GeometryMap: evaluation data not set");
        
        // gsMapData stores measure (Jacobian determinant * quadrature weight)
        // We need to extract just the determinant
        return m_data->measures(k);
    }

    std::string label() const
    {
        return "G";
    }
};

} // namespace Expr
} // namespace gismo
