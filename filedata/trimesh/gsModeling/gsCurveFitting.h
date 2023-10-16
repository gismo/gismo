/** @file gsCurveFitting.h

    @brief Fits a (closed) B-spline curve to some given data

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Kapl
*/

#pragma once

#include <iostream>
#include <stdexcept>

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsBSpline.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{
/**
    @brief Class for performing a least squares fit to get a open/closed B-Spline curve for some given data
    
    \ingroup Modeling
*/

template<class T>
class gsCurveFitting
{
public:
  /// default constructor
  gsCurveFitting(){
       //m_curve = NULL;
       m_closed=false;
    };

  /// constructor
  gsCurveFitting(gsMatrix<T> const & param_values, gsMatrix<T> const & points, gsKnotVector<T> const & knots, const bool & closed = false )
  {
      //m_curve = NULL;
      m_param_values=param_values;
      m_points=points;
      m_knots=knots;
      m_closed=closed;
   };

  /// Destructor
  ~gsCurveFitting() {
      //delete m_curve; // deleting the B-spline curve
  };


  /// computes the least squares fit for a (closed) B-spline curve
  void compute();

  // computes the least squares fit for a periodic B-spline curve (experimental)
  void compute_periodic();

  /// computes the approximation error of the fitted curve to the original point cloud
  void computeApproxError(T & error);

  /// set m_closed to true/or false since we want to have a closed/open curve
  void setClosedCurve(bool closed){
      m_closed=closed;
  };

  /// gives back the computed B-spline curve
  const gsBSpline<T>& curve() const { return m_curve; }

  /// returns the knot vector
  gsKnotVector<T> returnKnotVector() const {return m_knots;}

  /// returns the parameter values
  gsMatrix<T> returnParamValues() const {return m_param_values;}

  /// returns the points
  gsMatrix<T> returnPoints() const {return m_points;}

protected:

  /// pointer to a B-spline curve
  gsBSpline<T> m_curve;
  /// the parameter values of the point cloud
  gsMatrix<T> m_param_values;
  /// the points of the point cloud
  gsMatrix<T> m_points;
  /// the knot vector of the desired B-spline curve
  gsKnotVector<T> m_knots;
  /// closed or not closed curve
  bool m_closed;


}; // class gsCurveFitting

template<class T>
void gsCurveFitting<T>::compute()
{
    //degree of knot vector i.e. degree of wanted B-spline curve
    int m_degree=m_knots.degree();
    // number of points
    int num_points=m_points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=m_points.cols();
    // basis function definition
    gsBSplineBasis<T> *curveBasis = new gsBSplineBasis<T>(m_knots);

    //number of basis functions
    int num_basis=curveBasis->size();


    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    //computing the values of the basis functions at some position
    curveBasis->eval_into(m_param_values.transpose(),values);
    // which functions have been computed i.e. which are active
    curveBasis->active_into(m_param_values.transpose(),actives);

    //how many rows and columns has the A matrix and how many rows has the b vector
    int num_rows=num_basis;
    if(m_closed==true){
        num_rows=num_rows-m_degree;
    }

    //left side matrix
    gsMatrix<T> m_A(num_rows,num_rows);
    m_A.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_rows,m_dimension);
    m_B.setZero(); // enusure that all entris are zero in the beginning


    // building the matrix A and the vector b of the system of linear equations A*x==b(uses automatically the information of closed or not)
    for(index_t k=0;k<num_points;k++){
        for(index_t i=0;i<actives.rows();i++){
            m_B.row(actives(i,k)%num_rows) += values(i,k)*m_points.row(k);
            for(index_t j=0;j<actives.rows();j++){
                m_A(actives(i,k)%num_rows,actives(j,k)%num_rows) += values(i,k)*values(j,k);
            }
        }
    }

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);

    // making the coeffiecients ready if it was a closed or not closed curve
    gsMatrix<T> coefs=x;
    if(m_closed==true){
        coefs.conservativeResize(num_rows+m_degree,m_dimension);
        for(index_t i=0;i<m_degree;i++){
            for(index_t j=0;j<m_dimension;j++){
                coefs(num_rows+i,j)=coefs(i,j);
            }
        }

    }
    // finally generate the B-spline curve
    m_curve = gsBSpline<T >(*curveBasis, give(coefs));
    delete curveBasis;
}

template<class T>
void gsCurveFitting<T>::compute_periodic()
{
    //degree of knot vector i.e. degree of wanted B-spline curve
    int m_degree=m_knots.degree();
    // number of points
    int num_points=m_points.rows();
    // dimension of points will be also later dimension of coefficients
    int m_dimension=m_points.cols();
    // basis function definition
    gsBSplineBasis<T> curveBasis(m_knots, m_closed);

    //number of basis functions
    int num_basis=curveBasis.size();


    //for computing the value of the basis function
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    //computing the values of the basis functions at some position
    curveBasis.eval_into(m_param_values.transpose(),values);
    // which functions have been computed i.e. which are active
    curveBasis.active_into(m_param_values.transpose(),actives);

    //how many rows and columns has the A matrix and how many rows has the b vector
    int num_rows=num_basis;

    // No longer necessary:
    //if(m_closed==true){
    //    num_rows=num_rows-m_degree;
    //}

    //left side matrix
    gsMatrix<T> m_A(num_rows,num_rows);
    m_A.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_rows,m_dimension);
    m_B.setZero(); // enusure that all entris are zero in the beginning


    // building the matrix A and the vector b of the system of linear equations A*x==b(uses automatically the information of closed or not)
    for(index_t k=0;k<num_points;k++){
        for(index_t i=0;i<actives.rows();i++){
            m_B.row(actives(i,k)) += values(i,k)*m_points.row(k);
            for(index_t j=0;j<actives.rows();j++){
                m_A(actives(i,k),actives(j,k)) += values(i,k)*values(j,k);
            }
        }
    }

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)
    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);

    // making the coeffiecients ready if it was a closed or not closed curve
    gsMatrix<T> coefs=x;
    if(m_closed==true){
        coefs.conservativeResize(num_rows+m_degree,m_dimension);
        for(index_t i=0;i<m_degree;i++){
            for(index_t j=0;j<m_dimension;j++){
                coefs(num_rows+i,j)=coefs(i,j);
            }
        }

    }
    curveBasis.setPeriodic(false);
    //coefs = curveBasis->perCoefs(coefs);
    // finally generate the B-spline curve
    m_curve= gsBSpline<T>(curveBasis, give(coefs));
    m_curve.setPeriodic(m_closed);
}

template<class T>
void gsCurveFitting<T>::computeApproxError(T & error)
{
    gsMatrix<T> results;
    // the points of the curve for the corresponding parameter values
    m_curve.eval_into(m_param_values.transpose(),results);
    results.transposeInPlace();

    // computing the approximation error = sum_i ||x(u_i)-p_i||^2
    error = (m_points - results).squaredNorm();
}


template<class T>
std::ostream & operator <<( std::ostream &os, gsCurveFitting<T> const & cf)
{
    os <<"with "<< (cf.returnPoints()).rows() << " points of dimension " << (cf.returnPoints()).cols() <<
         ". The desired B-spline curve should have degree "<< (cf.returnKnotVector()).degree() <<" with the following knot vector: " << cf.returnKnotVector() << "\n";  // cf.points
    return os;
}

} // namespace gismo
