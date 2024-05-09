/** @file gsMultiBasis.h

    @brief Provides declaration of MultiBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsDofMapper.h>

namespace gismo
{


/**
 * @brief      This class describes a domain via the Möbius transformation.
 * NB: the class is only for 2d
 *
 * @tparam     T     { description }
 */
template <short_t DIM, class T>
class gsMobiusDomain : public gsFunction<T>
{
    using Base = gsFunction<T>;

public:
    // default constructor
    gsMobiusDomain()
    {
        GISMO_ASSERT(DIM==2,"The Möbius domain is only implemented in 2D");
        m_alpha.setOnes(2,DIM);
    }

    explicit gsMobiusDomain(const gsMatrix<T, 2, DIM> alpha)
    :
    m_alpha(alpha)
    {
        GISMO_ASSERT(DIM==2,"The Möbius domain is only implemented in 2D");
    }

    short_t domainDim() const override
    {
        return DIM;
    }

    short_t targetDim() const override
    {
        return 2;
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(DIM, u.cols());

//    T alpha = alpha_1 * t + alpha_2 * (1-t);
//    T beta  = beta_1 * s + beta_2 * (1-s);

        gsMatrix<T> alpha = m_alpha(0,0) * u.row(1).array() + m_alpha(1,0) * (1.0 - u.row(1).array());
        gsMatrix<T> beta  = m_alpha(0,1) * u.row(0).array() + m_alpha(1,1) * (1.0-u.row(0).array());

        gsMatrix<T> xi_denominator = 2 * alpha.array() * u.row(0).array() - u.row(0).array() - alpha.array();
        GISMO_ASSERT((xi_denominator.array()!=0).any(),"xi_denominator is zero!\n xi_denominator = "<<xi_denominator);
        result.row(0) = (alpha.array()-1)*u.row(0).array() / xi_denominator.array();
        gsMatrix<T> eta_denominator = 2 * beta.array() * u.row(1).array() - u.row(1).array() - beta.array();
        GISMO_ASSERT((eta_denominator.array()!=0).any(),"eta_denominator is zero!\n eta_denominator = "<<eta_denominator);
        result.row(1) = (beta.array()-1)*u.row(1).array() / eta_denominator.array();
    }

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(4, u.cols());

        gsMatrix<T> alpha = m_alpha(0,0) * u.row(1).array() + m_alpha(1,0) * (1.0-u.row(1).array());
        gsMatrix<T> beta  = m_alpha(0,1) * u.row(0).array() + m_alpha(1,1) * (1.0-u.row(0).array());

        gsMatrix<T> xi_denominator  = 2 * alpha.array() * u.row(0).array() - u.row(0).array() - alpha.array();
        gsMatrix<T> eta_denominator = 2 * beta.array() * u.row(1).array() - u.row(1).array() - beta.array();

        gsMatrix<T> xi  = (alpha.array()-1)*u.row(0).array() / xi_denominator.array();
        gsMatrix<T> eta = (beta.array()-1)*u.row(1).array() / eta_denominator.array();
        T deltaAlpha = m_alpha(0,0) - m_alpha(1,0);
        T deltaBeta  = m_alpha(0,1) - m_alpha(1,1);

        result.row(0) = (alpha.array()-1)*(xi_denominator.array() - (2*alpha.array()-1)*u.row(0).array())/(xi_denominator.array()*xi_denominator.array());
        result.row(1) = deltaAlpha*u.row(0).array()*(xi_denominator.array()-(alpha.array()-1)*(2*u.row(0).array()-1))/(xi_denominator.array()*xi_denominator.array());
        result.row(2) = deltaBeta*u.row(1).array()*(eta_denominator.array()-(beta.array()-1)*(2*u.row(1).array()-1))/(eta_denominator.array()*eta_denominator.array());
        result.row(3) = (beta.array()-1)*(eta_denominator.array()-(2*beta.array()-1)*u.row(1).array())/(eta_denominator.array()*eta_denominator.array());
    }

    /// Returns the controls of the function
    gsAsConstVector<T> controls() const override { return m_alpha.asVector(); };
    gsAsVector<T>      controls()       override { return m_alpha.asVector(); };

    /// Returns the \a i th control of the function
    const T & control(index_t i) const override { return m_alpha.coeff(i % DIM, math::floor(i / DIM)); }
          T & control(index_t i)       override { return m_alpha.coeffRef(i % DIM, math::floor(i / DIM)); }

    // const gsAsConstVector<T> controls() const { return m_alpha.asVector(); }
    //       gsAsVector<T>      controls()       { return m_alpha.asVector(); }


    /// Returns the number of controls of the function
    size_t nControls() const override
    {
        return m_alpha.size();
    }

  /// Returns the control derivative
//  virtual void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override
//  {
//    gsMatrix<T> tmp;
//
//    result.resize(targetDim()*nControls(), points.cols());
//    result.setZero();
//    for (index_t p = 0; p!=points.cols(); p++)
//    {
//      gsAsMatrix<T> res = result.reshapeCol(p,nControls(),targetDim());
//      for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
//        for (index_t d = 0; d!=m_domain.targetDim(); d++)
//          if (m_mapper.is_free(k,0,d))
//          {
//            m_domain.basis().evalSingle_into(k,points.col(p),tmp); // evaluate basis function k
//            res(m_mapper.index(k,0,d),d) = tmp(0,0); // tmp is a single value (1 point, 1 basis function)
//          }
//    }
//  }

protected:
//  T m_alpha1, m_alpha2, m_beta1, m_beta2;
    gsMatrix<T> m_alpha;
};

} // namespace gismo
