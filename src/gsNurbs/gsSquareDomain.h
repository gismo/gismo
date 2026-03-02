/** @file gsMultiBasis.h

    @brief Provides declaration of MultiBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
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
    // default constructor
    // gsSquareDomain()

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  domain  The domain
     */
    gsSquareDomain(const gsGeometry<T> & domain);

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  basis  The basis
     */
    gsSquareDomain(const gsBasis<T> & basis);


    // /**
    //  * @brief      Constructs a new instance.
    //  *
    //  * @param[in]  numElevation  The number elevation
    //  * @param[in]  numRefine     The number refine
    //  */
    // gsSquareDomain(index_t numElevation = 0, index_t numRefine = 0);

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

    // /// Sets the controls of the domain (does not override)
    // void setControls(const gsMatrix<T> & coefs)
    // {
    //     GISMO_ASSERT(coefs.rows()==m_domain.coefs().rows() && coefs.cols()==m_domain.coefs().cols(),"Wrong size of controls matrix");
    //     for (index_t i = 0; i!=m_indices.size(); i++)
    //         m_domain.coefs()(m_indices[i].first,m_indices[i].second) = coefs(m_indices[i].first,m_indices[i].second);
    // }

    /// Sets the controls of the domain
    void setControls(const gsVector<T> & controls) override;

    /// Gets the controls of the domain
    /// NOTE: This makes a copy
    gsVector<T> getControls() const override;

    /// Returns the \a i th control of the function
    // const T & control(index_t i) const override;
    //       T & control(index_t i)       override;


    // const gsVector<T> & parameters() const { return m_parameters; };
    //       gsVector<T> & parameters()       { return m_parameters; };

    /// Returns the number of controls of the function
    size_t nControls() const override;

    /// Returns the control derivative
    virtual void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override;

    /**
     * @brief Computes the derivative of det(J_sigma) w.r.t. the free controls
     *
     * Given sigma(xi,eta) = sum_k alpha_k * phi_k(xi,eta), the Jacobian is
     * J_sigma = d(sigma)/d(xi,eta). This method computes d(det(J_sigma))/d(alpha_i)
     * for each free control alpha_i.
     *
     * @param[in]  points  The evaluation points in the parameter domain
     * @param[out] result  A vector of size nControls x nPoints, where
     *                     result(i, p) = d(det(J_sigma))/d(alpha_i) at point p
     */
    void control_jacobian_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const;

    /// Perturb the control points by a \a factor
    void perturb(T factor = 1e-3);

private:
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