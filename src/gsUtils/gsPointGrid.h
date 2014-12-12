
#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo {

/// Construct a grid of points in the rectangle [a,b] using np[i] points in direction i
/// \param a
/// \param b
/// \param np
template<class T>
typename gsMatrix<T>::uPtr gsPointGrid( gsVector<T> const & a, gsVector<T> const & b, 
                                        gsVector<unsigned> const & np );


/// Construct a grid of points by coordinate vectors in the container cwise
template<class T>
inline typename gsMatrix<T>::uPtr gsPointGrid( std::vector< gsVector<T>* > const & cwise)
{
  gsMatrix<T> * res = new gsMatrix<T>;
  gsPointGrid(cwise, *res);
  return typename gsMatrix<T>::uPtr( res );
}

/// Construct a grid of points by coordinate vectors in the container
/// cwise, use out argument res
template<class T>
void gsPointGrid( std::vector< gsVector<T>* > const & cwise, gsMatrix<T>& res);

/// Construct a grid of points by coordinate vectors in the container cwise
template<class T>
inline typename gsMatrix<T>::uPtr gsPointGrid( std::vector< gsVector<T> > const & cwise)
{
  gsMatrix<T> * res = new gsMatrix<T>;
  gsPointGrid(cwise, *res);
  return typename gsMatrix<T>::uPtr( res );
}

/// Construct a grid of points by coordinate vectors in the container
/// cwise, use out argument res
template<class T>
void gsPointGrid( std::vector< gsVector<T> > const & cwise, gsMatrix<T>& res);

/// Compute the tensor product of the vectors in \a cwise and store it
/// in lexicographic order into the vector \a res
/// \todo move to gsTensorTools
template<class T>
void tensorProduct( std::vector< gsVector<T>* > const & cwise, gsVector<T>& res);

/// Specialization of the arguments for the 1D case
template<class T> inline
typename gsMatrix<T>::uPtr gsPointGrid( T const & t1, T const & t2, unsigned const & n = 100)
{
    gsVector<T> a(1) ; 
    gsVector<T> b(1) ;
    gsVector<unsigned> np(1);
    a<< t1;
    b<< t2;
    np<< n;
    return gsPointGrid(a,b,np);
}

/// Approximately uniformly spaced grid in every direction, with
/// approximately numPoints total points
template<typename T>
typename gsMatrix<T>::uPtr uniformPointGrid(const gsVector<T>& lower, const gsVector<T>& upper, int numPoints = 1000);


template<typename T>
gsVector<unsigned> uniformSampleCount (const gsVector<T>& lower, 
                                       const gsVector<T>& upper, 
                                       int numPoints = 1000);

template<typename T>
void uniformIntervals(const gsVector<T>& lower, const gsVector<T>& upper, 
                      std::vector< std::vector<T> >& intervals, int numIntervals = 1000);


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPointGrid.hpp)
#endif
