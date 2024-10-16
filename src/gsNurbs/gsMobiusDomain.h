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

#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

template <short_t DIM, class T>
class gsMobiusDomain : public gsFunction<T>
{
  using Base = gsFunction<T> ;

 public:
  // default constructor
  gsMobiusDomain() { m_alpha.setOnes(2,DIM); }

  explicit gsMobiusDomain(const gsMatrix<T, 2, DIM> alpha) : m_alpha(alpha) {}

//  gsMobiusDomain(index_t numElevation = 0, index_t numRefine = 0)
//  {
//    m_domain = *gsNurbsCreator<T>::BSplineSquare();
//    m_domain.degreeElevate(numElevation);
//    index_t numKts = pow(2, numRefine) - 1;
//    m_domain.uniformRefine(numKts);
//    // m_domain.uniformRefine(15);
//    gsInfo << " m_domain bi-degree = (" << m_domain.degree(0) <<", " << m_domain.degree(1) << ")\n";
//    gsInfo << " m_domain.coefsSize() = " << m_domain.coefsSize() << "\n";
////      gsDebugVar(m_domain.coefsSize());
//    // m_domain.uniformRefine();
//    // m_domain.uniformRefine();
//    // Mapper storing control points
//    m_mapper = gsDofMapper(m_domain.basis(),m_domain.targetDim());
//
//    gsMatrix<index_t> boundary = m_domain.basis().allBoundary();
//    for (index_t a = 0; a!=boundary.rows(); a++)
//      for (index_t d = 0; d!=m_domain.targetDim(); d++)
//        m_mapper.eliminateDof(boundary(a,0),0,d);
//    m_mapper.finalize();
//
//    m_parameters.resize(m_mapper.freeSize());
//    // std::vector<index_t> i(m_mapper.freeSize());
//    // std::vector<index_t> j(m_mapper.freeSize());
//    for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
//      for (index_t d = 0; d!=m_domain.targetDim(); d++)
//        if (m_mapper.is_free(k,0,d))
//        {
//          m_parameters[m_mapper.index(k,0,d)] = m_domain.coefs()(k,d);
//          // i[m_mapper.index(k,0,d)] = k; // i index of free entries
//          // j[m_mapper.index(k,0,d)] = d; // j index of free entries
//        }
//
//    // This is a way to cast only the free coefficients to a vector, and change an entry of that vector.
//    // However, it cannot be used in ''gsVector<T> & controls() override { return m_parameters; };''
//    //
//    // gsDebugVar(m_domain.coefs()(i,j).diagonal()(0));
//    // m_domain.coefs()(i,j).diagonal()(0) = 0.5;
//    // gsDebugVar(m_domain.coefs()(i,j).diagonal()(0));
//
//  }

//  const gsTensorBSpline<DIM,T> & domain() const
//  {
//    return m_domain;
//  }
//
//  gsMatrix<T> support() const override
//  {
//    return m_domain.support();
//  }
//

  short_t domainDim() const override
  {
    return m_domain.domainDim();
  }

  short_t targetDim() const override
  {
    return m_domain.domainDim();
  }

  void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
  {
    result.resize(2, u.cols());

//    T alpha = alpha_1 * t + alpha_2 * (1-t);
//    T beta  = beta_1 * s + beta_2 * (1-s);
    gsMatrix<T> alpha = m_alpha(0,0) * u.row(1).array() + m_alpha(1,0) * (1.0 - u.row(1).array());
    gsMatrix<T> beta  = m_alpha(0,1) * u.row(0).array() + m_alpha(1,1) * (1.0-u.row(0).array());

    gsMatrix<T> xi_denominator = 2 * alpha.array() * u.row(0).array() - u.row(0).array() - alpha.array();
    result.row(0) = (alpha.array()-1)*u.row(0).array() / xi_denominator.array();
    gsMatrix<T> eta_denominator = 2 * beta.array() * u.row(1).array() - u.row(1).array() - beta.array();
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

  void updateGeom(const gsMatrix<T>& alpha) {
    m_alpha = alpha;
  }

  void updateGeom(const gsAsConstVector<T>& alpha) {
    m_alpha(0,0) = alpha(0);
    m_alpha(1,0) = alpha(1);
    m_alpha(0,1) = alpha(2);
    m_alpha(1,1) = alpha(3);
  }

  /// Returns the controls of the function
  const gsVector<T> & controls() const { return m_parameters; };
  gsVector<T> & controls() { return m_parameters; };

  /// Returns the number of controls of the function
  size_t nControls() const
  {
    return m_mapper.freeSize();
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
  gsTensorBSpline<DIM,T> m_domain;
  gsDofMapper m_mapper;
  gsVector<T> m_parameters;
};

} // namespace gismo
