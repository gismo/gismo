/** @file gsSquareDomain.h

    @brief Provides declaration of the gsSquareDomain class, a free-form
    reparametrization of the [0,1]^2 parameter domain used as the sigma map
    in composed/relocated-mesh (r-adaptive) discretizations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsDofMapper.h>
#include <gsIO/gsOptionList.h>
#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

template <class T>
class gsSquareDomain : public gsFunction<T>
{
    using Base = gsFunction<T> ;

public:

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  domain  The domain
     * @param[in]  slide   Value of the "Slide" option (boundary controls
     *                     slide along the boundary), applied before the dof
     *                     mapper is built. Default false.
     */
    gsSquareDomain(const gsGeometry<T> & domain, bool slide = false);

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  basis  The basis
     * @param[in]  slide  Value of the "Slide" option (boundary controls
     *                    slide along the boundary), applied before the dof
     *                    mapper is built. Default false.
     */
    gsSquareDomain(const gsBasis<T> & basis, bool slide = false);

    // Copy constructor
    gsSquareDomain(const gsSquareDomain<T> & other);

    // Copy assignment
    gsSquareDomain<T> & operator=(const gsSquareDomain<T> & other);

    GISMO_CLONE_FUNCTION(gsSquareDomain)

    /// @brief Couples the DoFs along the interface \a bi
    /// @param bi the interface
    void addInterface(const boundaryInterface & bi);

    gsOptionList & options();

    void applyOptions();

    const gsGeometry<T> & domain() const;

    const gsDofMapper & mapper() const;

    gsMatrix<T> support() const override;

    short_t maxDegree() const;

    short_t domainDim() const override;

    short_t targetDim() const override;

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override;

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override;

    // Exact second derivative (delegates to the underlying B-spline geometry).
    // Without this override gsFunction would fall back to finite differences,
    // which would spoil the second-derivative chain rule of gsComposedBasis.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override;

    /** @brief Evaluates values and derivatives in ONE pass.
     *
     * Without this override, gsFunctionSet::compute() falls back to calling
     * eval_into() / deriv_into() / deriv2_into() in turn, and each of those
     * delegates to the underlying B-spline geometry SEPARATELY -- so asking for
     * value+derivative evaluates the sigma basis twice, and
     * value+derivative+second derivative three times.  Delegating to
     * gsGeometry::compute() instead runs the basis' evalAllDers_into() once and
     * contracts each requested order with the control points.
     *
     * This matters because gsOptMesh's objective and gradient sweeps call
     * compute() on this object for every element on every optimizer iteration,
     * making it a hot-path routine on which redundant basis evaluation is
     * directly visible in wall time.
     */
    void compute(const gsMatrix<T> & u, gsFuncData<T> & out) const override;

    /// Sets the controls of the domain
    void setControls(const gsVector<T> & controls);

    /// Gets the controls of the domain
    /// NOTE: This makes a copy
    gsVector<T> getControls() const;

    /// Returns the number of controls of the function
    size_t nControls() const;

    /// Returns the control derivative
    virtual void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const;

    /// @brief det(J) at each column of \a jacobians, which must be in the
    ///        component-major (row c*dd+j) layout deriv_into produces.
    /// @param[in]  jacobians  (dd*dd) x nPoints. Precondition: dd >= 1 and
    ///                        jacobians.rows() == dd*dd.
    /// @param[in]  dd         the square Jacobian's dimension
    /// @param[out] result     1 x nPoints
    static void detFromJacobian_into(const gsMatrix<T> & jacobians, short_t dd,
                                     gsMatrix<T> & result);

    /// @brief det(J_sigma) at \a points. Result is 1 x points.cols().
    void detJacobian_into(const gsMatrix<T> & points, gsMatrix<T> & result) const;

    /**
     * @brief Computes the derivative of det(J_sigma) w.r.t. the free controls
     *
     * Given sigma(xi,eta) = sum_k alpha_k * phi_k(xi,eta), the Jacobian is
     * J_sigma = d(sigma)/d(xi,eta). This method computes d(det(J_sigma))/d(alpha_i)
     * for each free control alpha_i.
     *
     * @param[in]  points    The evaluation points in the parameter domain
     * @param[out] result    A vector of size nControls x nPoints, where
     *                       result(i, p) = d(det(J_sigma))/d(alpha_i) at point p
     * @param[out] jacobian  Optional: if non-null, also filled with J_sigma at
     *                       \a points in the same layout as deriv_into's output
     *                       (component-major: row c*domainDim+j). Computing
     *                       J_sigma is a byproduct of this method's own basis
     *                       evaluation, so a caller that needs both J_sigma and
     *                       this derivative (the common case: a fold-barrier
     *                       gradient) gets both for the cost of one basis pass
     *                       instead of two.
     */
    void detJacobianDeriv_into(const gsMatrix<T> & points, gsMatrix<T> & result,
                               gsMatrix<T> * jacobian = nullptr) const;

    /// Perturb the control points by a \a factor
    void perturb(T factor = 1e-3);

    /**
     * @brief Minimum of det(J_sigma), sampled on an \a nPerElement^d point
     *        grid PER ELEMENT of sigma's knot mesh. NECESSARY but never
     *        SUFFICIENT fold check (a fold between sample points is
     *        invisible to it). Derivation and the certified-bound
     *        comparison against minDetJCoefficient(): \ref adaptparam_foldcert.
     *
     * @param[in] nPerElement  Number of sample points per direction, per
     *                         element of sigma's knot mesh. Must be >= 2.
     * @return The minimum sampled value of det(J_sigma) over the domain.
     */
    T minJacobian(index_t nPerElement) const;

    /**
     * @brief Represents det(J_sigma) EXACTLY (to round-off) as a scalar
     *        tensor B-spline, in the one spline space that actually
     *        contains it. Bernstein/Leibniz construction: \ref adaptparam_foldcert.
     *
     * @param[in] keepBezier  If true, skip the interior-knot compression
     *                        inside \c gsTensorBSpline::multiply and return
     *                        the raw Bezier-form product -- the tightest
     *                        possible coefficient bound for
     *                        minDetJCoefficient(), at the cost of a larger basis.
     * @return A scalar (targetDim() == 1) tensor B-spline geometry equal to
     *         det(J_sigma) to round-off, owned by the caller.
     */
    typename gsGeometry<T>::uPtr detJacobianSpline(bool keepBezier = false) const;

    /**
     * @brief Guaranteed LOWER BOUND on det(J_sigma) over the whole domain,
     *        computed WITHOUT sampling anywhere: the partition-of-unity
     *        proof and the certified-vs-sampled ordering
     *        \code minDetJCoefficient() <= true min_x det(J_sigma) <= minJacobian(nPerElement) \endcode
     *        are derived in \ref adaptparam_foldcert.
     *
     * @return The minimum B-spline coefficient of det(J_sigma)'s exact
     *         representation -- a certified lower bound, never an
     *         approximation.
     */
    T minDetJCoefficient() const;

private:
    /// @brief Dimension-generic implementation of detJacobianSpline(),
    ///        dispatched on the runtime domainDim() by detJacobianSpline().
    template<short_t d>
    typename gsGeometry<T>::uPtr _detJacobianSplineImpl(bool keepBezier) const;

    /// @brief Initialize the dof mapper
    /// @param domain The domain as a tensor B-spline
    /// @param mapper The dof mapper (output)
    void _initMapper(const gsGeometry<T> & domain, gsDofMapper & mapper) const;

    /// @brief Initialize the indices of the free control points
    /// @param domain   The domain as a tensor B-spline
    /// @param mapper   The dof mapper
    /// @param indices  The indices of the free control points (output)
    void _initIndices(const gsGeometry<T> & domain, const gsDofMapper & mapper, std::vector<std::pair<index_t,index_t>> & indices) const;

protected:
    typename gsGeometry<T>::uPtr m_domain;
    gsDofMapper m_mapper;
    std::vector<boundaryInterface> m_interfaces;
    std::vector<std::pair<index_t,index_t>> m_indices;
    gsOptionList m_options;
};

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSquareDomain.hpp)
#endif