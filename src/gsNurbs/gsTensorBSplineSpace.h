/** @file gsTensorBSplineSpace.h

    @brief Tensor B-spline space descriptions and exact transfer maps.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <vector>

namespace gismo
{

/// @brief Selects the knot multiplicities used by derived tensor spaces.
///
/// Minimal preserves the continuity implied by the operation whenever possible.
/// Bezier inserts all interior breakpoints to multiplicity degree+1 in the
/// resulting space.
enum class gsTensorBSplineSpacePolicy
{
    Minimal,
    Bezier
};

/// @brief Lightweight description of a tensor-product B-spline space.
///
/// The spec stores exactly one knot vector per tensor direction and can be
/// converted back to a gsTensorBSplineBasis. It is intended for construction
/// planning: derived spaces, common superspaces, and exact transfer maps.
///
/// Current scope is deliberately narrow: non-periodic gsTensorBSplineBasis
/// only, no multi-patch coupling, no geometry map, and no rational weights.
///
/// \tparam d parameter dimension of the tensor-product space
/// \tparam T scalar type of the knot values
template <short_t d, class T>
class gsTensorBSplineSpaceSpec
{
public:
    typedef gsKnotVector<T> KnotVectorType;
    typedef gsTensorBSplineBasis<d,T> Basis;

public:
    gsTensorBSplineSpaceSpec();

    /// @brief Construct from one knot vector per direction.
    ///
    /// Throws if \a knots does not contain exactly d entries.
    explicit gsTensorBSplineSpaceSpec(const std::vector<KnotVectorType>& knots);

    /// @brief Copy the tensor knot vectors from a non-periodic basis.
    ///
    /// Throws for periodic bases, because the current space operations work on
    /// explicit clamped knot vectors and coefficient transfer by knot insertion.
    static gsTensorBSplineSpaceSpec fromBasis(const Basis& basis);

    /// @brief Build a tensor B-spline basis with the stored knot vectors.
    Basis toBasis() const;

    /// @brief Return the knot vector in direction \a dir.
    const KnotVectorType& knots(short_t dir) const;

    /// @brief Return the polynomial degree in direction \a dir.
    short_t degree(short_t dir) const;

    /// @brief Return the number of basis functions in the tensor space.
    index_t size() const;

    bool isValid() const { return m_knots.size() == static_cast<size_t>(d); }

private:
    std::vector<KnotVectorType> m_knots;
};

/// @brief Return the full Bezier refinement of \a spec.
///
/// Every interior breakpoint receives multiplicity degree+1 in its direction.
/// Boundary multiplicities are preserved.
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> bezierSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec);

/// @brief Return a tensor space containing pointwise products from \a a and \a b.
///
/// The input spaces must have identical parameter domains in each direction.
/// Breakpoints are merged. The resulting degree is p_a+p_b per direction.
/// With Minimal policy, interior multiplicities are chosen to preserve the
/// lower continuity of both factors. With Bezier policy, all interior
/// multiplicities are degree+1.
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> productSpace(
    const gsTensorBSplineSpaceSpec<d,T>& a,
    const gsTensorBSplineSpaceSpec<d,T>& b,
    gsTensorBSplineSpacePolicy policy = gsTensorBSplineSpacePolicy::Minimal);

/// @brief Return the tensor space containing integer powers from \a spec.
///
/// The power must be at least one. The operation is formed by repeated
/// productSpace() calls, so degrees and multiplicities follow product rules.
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> powerSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    index_t power,
    gsTensorBSplineSpacePolicy policy = gsTensorBSplineSpacePolicy::Minimal);

/// @brief Return the tensor space for mixed derivatives of \a spec.
///
/// \a orders[k] is the derivative order in direction k. Each order must not
/// exceed the corresponding degree. Boundary knots are dropped according to the
/// derivative order, matching the exact derivative basis produced by
/// gsTensorBSpline::grad().
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> derivativeSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    const gsVector<index_t,d>& orders);

/// @brief Return a common scalar space for the divergence of a vector field.
///
/// The input spec describes the component basis of a d-vector field. The
/// returned space contains sum_k d f_k / d x_k after degree-elevating each term
/// back in its differentiated direction.
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> divergenceSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec);

/// @brief Return a common scalar space for the Laplacian of \a spec.
///
/// The returned Minimal space contains sum_k d^2 f / d x_k^2 after
/// degree-elevating each second derivative term back in its differentiated
/// direction. Bezier policy additionally returns the full Bezier refinement.
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> laplacianSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    gsTensorBSplineSpacePolicy policy = gsTensorBSplineSpacePolicy::Minimal);

/// @brief Return a common tensor superspace of all input specs.
///
/// All input domains must match. Degrees are elevated to the maximum degree in
/// each direction and knot multiplicities are the maximum multiplicity required
/// by any input after that degree elevation.
template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> commonSpace(
    const std::vector< gsTensorBSplineSpaceSpec<d,T> >& specs);

/// @brief Build the exact coefficient transfer from \a source to \a target.
///
/// The target must be a tensor superspace of the source: no lower degrees, no
/// missing source breakpoints, and no lower knot multiplicities after degree
/// elevation. The returned sparse matrix has size target.size() x source.size()
/// and maps source coefficients to target coefficients.
///
/// Complexity note: the implementation applies exact degree elevation and knot
/// insertion to identity coefficients, so intermediate storage can be dense in
/// target.size() x source.size(). Treat this as a construction-time map, not as
/// a repeated hot-path operation.
template <short_t d, class T>
gsSparseMatrix<T,RowMajor> transferMatrix(
    const gsTensorBSplineSpaceSpec<d,T>& source,
    const gsTensorBSplineSpaceSpec<d,T>& target);

} // namespace gismo

#include GISMO_HPP_HEADER(gsTensorBSplineSpace.hpp)
