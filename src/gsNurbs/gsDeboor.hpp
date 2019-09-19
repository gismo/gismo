/** @file gsDeboor.hpp

    @brief Implementation of deBoor and tensor deBoor algorithm

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

//#include <gsNurbs/gsBSplineAlgorithms.h>
#include <gsNurbs/gsBoehm.h>

#include <gsUtils/gsCombinatorics.h>

namespace gismo {


/// Executes deBoor's algorithm on the absissae (row vector) u, knot vector \a knots,
/// degree \a deg and coefficients matrix \a coefs
// to do: make local version with knot iterators as arguments
template<class T, class KnotVectorType> inline
void gsDeboor( 
    const gsMatrix<T> &u,
    const KnotVectorType & knots,
    int deg,
    const gsMatrix<T> & coefs,
    gsMatrix<T>& result
    )
{
  GISMO_ASSERT( u.rows() == 1, "Waiting for 1D values" ) ;
  GISMO_ASSERT( coefs.rows() == index_t(knots.size() - deg-1) ,
                "coefs.rows(): " << coefs.rows() << ", knots.size(): " << knots.size() << ", deg: " << deg ) ;
  
  result.resize( coefs.cols(), u.cols() ) ;
  int ind, k;
  gsMatrix<T> points;
  T tmp;
  
  for ( index_t j=0; j< u.cols(); j++ ) // for all points (entries of u)
  {
    //De Boor's algorithm for parameter u(0,j)
    GISMO_ASSERT( (u(0,j)>knots[deg]-T(1e-4) ) && (u(0,j) < *(knots.end()-deg-1)+T(1e-4) ), 
                  "Parametric point "<< u(0,j) <<" outside knot domain ["
                  << knots[deg]<<","<<*(knots.end()-deg-1) <<"]."); 

    ind = (knots.iFind( u(0,j) ) - knots.begin()) - deg;
    
    //int s= knots.multiplicity( u(0,j) ) ; // TO DO: improve using multiplicity s
    points = coefs.middleRows( ind, deg+1 );
    
    for ( int r=0; r< deg; r++ ) // triangle step
      //for ( int r=s; r< deg; r++ ) // TO DO: improve using multiplicity s
      for ( int i=deg; i>r; i-- ) // recursive computation
        {
          k= ind + i;
          tmp= ( u(0,j) -  knots[k] ) / (knots[k+deg-r]-knots[k]);
          points.row(i) = (T(1)-tmp)*points.row(i-1) + tmp* points.row(i) ;
        }        
    result.col(j)= points.row(deg);
  }
}

/// Executes deBoor's algorithm on the absissae (row vector) u, for
/// the derivative of the B-spline defined by knot vector \a knots,
/// degree \a deg and coefficients matrix \a coefs
template<class T, class KnotVectorType> inline
void gsDeboorDeriv( 
    const gsMatrix<T> &u,
    const KnotVectorType & knots,
    int deg,
    const gsMatrix<T> & coefs,
    gsMatrix<T>& result
    )
{
  GISMO_ASSERT( u.rows() == 1, "Waiting for 1D values" ) ;
  GISMO_ASSERT( coefs.rows() == index_t(knots.size() - deg-1) ,
                "coefs.rows(): " << coefs.rows() << ", knots.size(): " << knots.size() << ", deg: " << deg ) ;

  deg--;// Derivative has degree reduced by one

  result.resize( coefs.cols(), u.cols() ) ;
  int ind, k;
  gsMatrix<T> points(deg+1, coefs.cols() );
  T tmp;

  for ( index_t j=0; j< u.cols(); j++ ) // for all points (entries of u)
  {
      //De Boor's algorithm for parameter u(0,j)
      GISMO_ASSERT( (u(0,j)>knots[deg]-T(1e-4) ) && (u(0,j) < *(knots.end()-deg-2)+T(1e-4) ), 
                    "Parametric point "<< u(0,j) <<" outside knot domain ["
                    << knots[deg]<<","<<*(knots.end()-deg-2) <<"]."); 
      
      ind = knots.findspan( u(0,j) ) - deg - 1;
      
      //int s= knots.multiplicity( u(0,j) ) ; // TO DO: improve using multiplicity s

      // Get the coefficients of the first derivative
      for ( int r=ind; r< deg+ind+1; r++ )
          points.row(r-ind) = ( T(deg+1) / (knots[r+deg+2]-knots[r+1]) )
                               * ( coefs.row(r+1) - coefs.row(r) );
      
      for ( int r=0; r< deg; r++ ) // triangle step
          //for ( int r=s; r< deg; r++ ) // TO DO: improve using multiplicity s
          for ( int i=deg; i>r; i-- ) // recursive computation
          {
              k= ind + i;
              tmp= ( u(0,j) -  knots[k] ) / (knots[k+deg+1-r]-knots[k]);
              points.row(i) = (T(1)-tmp)*points.row(i-1) + tmp* points.row(i) ;
          }        
      result.col(j)= points.row(deg);
  }
}
    

// =============================================================================
// ===== temporal version of gsTensorDeboor
// =============================================================================


template<short_t d, typename T, typename KnotVectorType, typename Mat>
inline
void gsTensorDeboor( //LC
        const gsMatrix<T>& u,
        const gsTensorBSplineBasis<d, T>& base,
        const Mat& coefs, // works just for 1D coefficients
        gsMatrix<T>& result
        )
{

    // Default version - linear combination of basis functions
    result.resize( coefs.cols(), u.cols() ) ;
    gsMatrix<T> B ;
    gsMatrix<index_t> ind;

    // "eval" of gsTensorBasis
    base.eval_into(u, B);   // col j = nonzero basis functions at column point u(..,j)
    base.active_into(u, ind);// col j = indices of active functions at column point u(..,j)

//    gsDebug << "u: \n" << u << std::endl;
//    gsDebug << "B: \n" << B << std::endl;
//    gsDebug << "ind: \n" << ind << std::endl;

    for ( index_t j=0; j< u.cols() ; j++ ) // for all points (columns of u)
    {
        result(0, j) = coefs(ind(0, j)) * B(0, j);
        for ( index_t i = 1; i < ind.rows(); ++i ) // for all non-zero basis functions
            result(0, j) += coefs(ind(i, j)) * B(i, j);
    }
}


// =============================================================================
// ===== derivatives for BSpline
// =============================================================================

template<short_t d, typename T, typename KnotVectorType, typename Mat>
inline //LC
void gsTensorDeriv_into(const gsMatrix<T>& u,
                        const gsTensorBSplineBasis<d, T>& base,
                        const Mat& coefs,
                        gsMatrix<T>& result)
{
    const unsigned nPts = u.cols(); // number of points
    const unsigned nDer = d; // number of derivatives

    result.setZero(d, nPts);
    gsMatrix<T> deriv; // derivatives
    gsMatrix<index_t> ind;

    base.deriv_into(u, deriv);   // col j = nonzero basis functions at column point u(..,j)
    base.active_into(u, ind);// col j = indices of active functions at column point u(..,j)

//    gsDebug << "jaka n = d: " << d << std::endl;
//    gsDebug << "ind: " << ind.rows() << " x " << ind.cols() << std::endl;
//    gsDebug << "B: " << B.rows() << " x " << B.cols() << std::endl;


    for (unsigned j = 0; j < nPts; ++j)
        for (unsigned k = 0; k < nDer; ++k) // for all rows of the jacobian
            for (index_t i = 0; i < ind.rows(); i++) // for all nonzero basis functions)
                    result(k, j) += coefs(ind(i, j)) * deriv(k + i * nDer, j);
}


// =============================================================================
// ===== second derivatives for BSpline
// =============================================================================

template<short_t d, typename T, typename KnotVectorType, typename Mat>
inline
void gsTensorDeriv2_into(const gsMatrix<T>& u,
                        const gsTensorBSplineBasis<d, T>& base,
                        const Mat& coefs,
                        gsMatrix<T>& result)
{//LC
    const unsigned nPts = u.cols(); // number points
    const unsigned n2der = (d * (d + 1)) / 2; // number of second derivatives

    result.setZero(n2der, nPts);
    gsMatrix<T> deriv2; // second derivatives
    gsMatrix<index_t> ind;


    base.deriv2_into(u, deriv2);
    base.active_into(u, ind);

    for (unsigned j = 0; j < nPts; j++)
        for (unsigned k = 0; k < n2der; k++)
            for (index_t i = 0; i < ind.rows(); i++)
                result(k, j) += coefs(ind(i, j)) * deriv2(k + i * n2der, j);
}

// =============================================================================
// other version of evaluation via knot insertion
// =============================================================================

template<short_t d, typename T, typename KnotVectorType, typename Mat>
inline
void gsTensorDeboor_v2(
        const gsMatrix<T>& u,
        const gsTensorBSplineBasis<d, T>& base,
        const Mat& coefs,
        gsMatrix<T>& result
        )
{

//    int log = -3;

//    if (0 < log)
//    {
//        gsDebug << "coefs: " << coefs.rows() << " x " << coefs.cols()
//                  << std::endl;
//    }

    // **************************************************************
    // ** find out how many coefficients one needs for calculation **

    result.resize(coefs.cols(), u.cols());
    gsVector<index_t, d> size_of_tmp_coefs(d);

    unsigned nmb_of_tmp_coefs = 1;
    for (unsigned dim = 0; dim < d; dim++)
    {
        const KnotVectorType& kv = base.knots(dim);
        size_of_tmp_coefs(dim) =  kv.degree() + 1;
        nmb_of_tmp_coefs *= size_of_tmp_coefs(dim);
    }

    Mat tmp_coefs(nmb_of_tmp_coefs, coefs.cols());
    tmp_coefs.fill(0);

    // for all points
    for (index_t i = 0; i < u.cols(); i++)
    {
//        if (0 < log)
//        {
//            gsDebug << "\n---------------------------------------------"
//                      << std::endl;

//            gsDebug << "u( " << i << " ): ";
//            for (unsigned dim = 0; dim < d; dim++)
//                gsDebug << u(dim, i) << " ";
//            gsDebug << std::endl;
//        }

        std::vector<bool> is_last(d, false);

        for (short_t dim = 0; dim < d; dim++)
        {
            const KnotVectorType& kv = base.knots(dim);
//            if (u(dim, i) == *(--kv.end()))
            if ( u(dim, i) == kv.last() )
            {
                is_last[dim] = true;
            }
        }



        gsVector<index_t, d> low, upp;
        base.active_cwise(u.col(i), low, upp);

//        if (1 < log)
//        {
//            gsDebug << "low: " << low.transpose() << "\nupp: "
//                      << upp.transpose() << std::endl;
//        }
//        if (0 < log)
//        {
//            gsDebug << "size_of_coefs: " << size_of_coefs.transpose() << "\n"
//                      << "number_of_coefs: " << nmb_of_coefs
//                      << std::endl;
//        }


        //*******************************************************
        //** copy appropriate coefficients to temporary matrix **

        // position of the coefficient
        gsVector<index_t, d> coefs_position(low);

        gsVector<index_t, d> zero(d);
        zero.fill(0);

        // position in tmp_coefs
        gsVector<index_t, d> tmp_coefs_position(zero);

        // strides in tmp_coefs
        gsVector<index_t, d> tmp_coefs_str(d);
        bspline::buildCoeffsStrides<d>(size_of_tmp_coefs, tmp_coefs_str);

        gsVector<index_t, d> tmp_coefs_last(d);
        bspline::getLastIndexLocal<d>(size_of_tmp_coefs, tmp_coefs_last);

        do
        {
//            if (10 < log)
//            {
//                gsDebug << "coefs_position: "
//                          << coefs_position.transpose() << "\n"
//                          << "tmp_coefs_position: "
//                          << tmp_coefs_position.transpose() << "\n"
//                          << std::endl;
//            }

            unsigned flat_index = base.index(coefs_position);
            unsigned tmp_coefs_index = bspline::getIndex<d>(tmp_coefs_str,
                                                            tmp_coefs_position);

            tmp_coefs.row(tmp_coefs_index) = coefs.row(flat_index);


            nextCubePoint<gsVector<index_t, d> >(tmp_coefs_position, zero,
                                                  tmp_coefs_last);
        } while(nextCubePoint<gsVector<index_t, d> >(coefs_position, low, upp));


//        if (1 < log)
//        {
//            gsDebug << "tmp_coefs: \n" << tmp_coefs << std::endl;
//        }


        // ********************************************************
        // ** perform knot insertion appropriate number of times **

        gsVector<index_t, d> start(d);
        start.fill(0);
        gsVector<index_t, d> end(size_of_tmp_coefs);
        for (short_t dim = 0; dim < d; dim++)
            end(dim)--;

        for (short_t dim = 0; dim < d; dim++) //in each dimension insert a knot
        {
            if (is_last[dim])
            {
                start(dim) = end(dim);
            }
            else
            {
                 const KnotVectorType& kv = base.knots(dim);

                 gismo::gsTensorInsertKnotDegreeTimes<d, T, KnotVectorType, Mat>
                         (kv, tmp_coefs, size_of_tmp_coefs, u(dim, i), dim,
                          start, end);

                 start(dim) = 0;
                 end(dim) = 0;
            }

        }

        // index of a point at a parameter u.col(i)
        index_t flat_index = bspline::getIndex<d>(tmp_coefs_str,
                                                   start);

//        if (0 < log)
//        {
//            gsDebug << "Point is: \n"
//                      << tmp_coefs.row(flat_index)
//                      << std::endl;
//        }

        result.col(i) = tmp_coefs.row(flat_index).transpose();

    } // end for each point
}



}; //end namespace gismo
