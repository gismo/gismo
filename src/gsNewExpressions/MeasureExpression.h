/** @file MeasureExpression.h

    @brief Expression for the measure (Jacobian determinant) of a geometry map

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

// Forward declaration
template <class T> class MeasureExpression;

// ExpressionTraits for MeasureExpression
template <class T>
struct ExpressionTraits<MeasureExpression<T>>
{
    typedef GeometryMap<T> ExprType;
    
    typedef T Scalar;
    static constexpr size_t Order = 0; // Measure is scalar-valued
    static constexpr SpaceType Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true; // Geometry doesn't change during assembly
};

/**
 * @brief Expression for the measure (Jacobian determinant) of a geometry map
 * 
 * This represents |det(DG)| where DG is the Jacobian of the geometry map G.
 * It is used in integration: ∫_Ω f dx = ∫_Ω̂ f(G(ξ)) |det(DG(ξ))| dξ
 * 
 * @tparam T The scalar type
 */
template<class T>
class MeasureExpression : public UnaryOperator<MeasureExpression<T>>
{
    using Base = UnaryOperator<MeasureExpression<T>>;
    
public:
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    typedef typename ExpressionTraits<MeasureExpression<T>>::ExprType ExprType;

    /**
     * @brief Constructor
     * @param geomMap The geometry map expression
     */
    MeasureExpression(const GeometryMap<T> & geomMap)
    : Base(geomMap)
    {
    }

    const GeometryMap<T> & geomMap() const 
    { 
        return static_cast<const GeometryMap<T>&>(this->expr());
    }

    void parse(gismo::ExpressionHelper<T> & helper) const
    {
        // Parse the geometry map
        this->expr().parse(helper);
        
        // TODO: Set flags to compute measure/Jacobian determinant
        // helper.setFlag(NEED_MEASURE);
    }

    void print(std::ostream & os) const
    {
        os << "meas(";
        this->expr().print(os);
        os << ")";
    }

    std::array<size_t, 0> sizes() const
    {
        return std::array<size_t, 0>(); // Scalar
    }

    size_t domainDim() const { return this->expr().domainDim(); }

    /// Evaluate the measure at quadrature point k
    /// Returns the Jacobian determinant |det(DG)|
    ExpressionResult<T> eval(const index_t k) const
    {
        // Get the Jacobian determinant from the geometry map data
        return geomMap().jacobianDet(k);
    }

    std::string label() const
    {
        return "meas(G)";
    }
};

/**
 * @brief Helper function to create a measure expression
 * @param G The geometry map
 * @return MeasureExpression wrapping the geometry map
 */
template<class T>
MeasureExpression<T> meas(const GeometryMap<T> & G)
{
    return MeasureExpression<T>(G);
}

} // namespace Expr
} // namespace gismo
