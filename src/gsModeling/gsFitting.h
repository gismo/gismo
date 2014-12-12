/** @file gsFitting.h

    @brief Provides declaration of data fitting algorithms by least
    squares approximation.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl, G. Kiss, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <vector>


namespace gismo
{

/**
    Class for performing a least squares fit to get a open/closed
    B-Spline curve for some given data
*/
template<class T>
class gsFitting
{
public:
  /// default constructor
  gsFitting()
  {
      m_basis = NULL;
      m_result= NULL ;
  }

  /// constructor
  gsFitting(gsMatrix<T> const & param_values, 
            gsMatrix<T> const & points, 
            gsBasis<T>  & basis);

  /// Destructor
  virtual ~gsFitting();

public:

  /// Computes the least squares fit for a gsBasis
  void compute(T lambda = 0);

  /// Computes the approximation error of the fitted curve to the original point cloud
  void computeApproxError(T & error, int type = 0) const;

  ///return the errors for each point
  void get_Error(std::vector<T>& errors, int type = 0) const;


  /// Computes the least squares fit for a gsBasis
  void iterativeCompute( T const & tolerance, unsigned const & num_iters = 10);


public:

  /// gives back the computed approximation
  gsGeometry<T> * result() const { return m_result; }

  /// Returns the basis of the approximation
  const gsBasis<T> & getBasis() const {return *m_basis;}

  void setBasis(gsBasis<T> & basis) {m_basis=&basis;}

  /// returns the parameter values
  gsMatrix<T> & getreturnParamValues() {return m_param_values;}
  gsMatrix<T> & returnParamValues() {return m_param_values;}

  /// returns the points
  gsMatrix<T> returnPoints() const {return m_points;}


protected:

  /// the parameter values of the point cloud
  gsMatrix<T> m_param_values;
  /// the points of the point cloud
  gsMatrix<T> m_points;
  /// Pointer keeping the basis
  gsBasis<T> * m_basis;
  /// Pointer keeping the resulting geometry
  gsGeometry<T> * m_result;

private:
  void applySmoothing(T lambda, gsSparseMatrix<T> & A_mat);
    //void applySmoothing(T lambda, gsMatrix<T> & A_mat);
}; // class gsFitting


}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFitting.hpp)
#endif
