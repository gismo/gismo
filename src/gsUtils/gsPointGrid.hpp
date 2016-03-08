/** @file gsPointGrid.hpp

    @brief Provides implementation of functions which generate and iterate over structured point data

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/
#pragma once

#include <gsUtils/gsCombinatorics.h>

namespace gismo 
{

template<class T>
typename gsMatrix<T>::uPtr gsPointGrid( gsVector<T> const & a, 
                                        gsVector<T> const & b, 
                                        gsVector<unsigned> const & np )
{
  const index_t d = a.size();
  assert( d == b.size() );
  assert( d == np.size() );

  gsVector<T> span = b - a;

  gsMatrix<T> * res = new gsMatrix<T>( d, np.prod() );

  gsVector<unsigned> v(d);
  v.setZero();

  gsVector<T> segments ( (np.array() - 1).matrix().cast<T>() );
  for (index_t i = 0; i < d; ++i)
      if (segments[i] == 0)
          segments[i] = 1;     // avoid division by zero

  unsigned r = 0;
  do {
      // Make tensor product of sample points
      //res->col(r) = a.array() + v.array().cast<T>() * span.array() / segments.array();
      for (index_t i = 0; i < d; ++i)
          if ( v[i] == 0 )
              res->coeffRef(i,r) = a[i]; // avoid numerical error in the start of the interval
          else if ( v[i] == np[i]-1 )
              res->coeffRef(i,r) = b[i]; // avoid numerical error in the end of the interval
          else
              res->coeffRef(i,r) = a[i] + v[i] * span[i] / segments[i];
      ++r ;
  } while ( nextLexicographic(v, np) );

  return typename gsMatrix<T>::uPtr( res );
}


template<class T>
void gsPointGrid( std::vector< gsVector<T>* > const & cwise, gsMatrix<T>& res)
{
  unsigned d= cwise.size();

  gsVector<unsigned> np(d);
  for (unsigned i=0; i<d; ++i )
    np[i] = cwise[i]->size();
  
  gsVector<unsigned> v(d);
  v.setZero();

  res.resize( d, np.prod() );

  unsigned r = 0;
  do {
    // Make tensor product of sample points
    for (unsigned i=0; i<d; ++i )
      res(i,r) =  (*cwise[i])( v[i] )  ;
    ++r ;
  } while (nextLexicographic(v, np));
}

template<class T>
void gsPointGrid( std::vector< gsVector<T> > const & cwise, gsMatrix<T>& res)
{
  unsigned d= cwise.size();

  gsVector<unsigned> np(d);
  for (unsigned i=0; i<d; ++i )
    np[i] = cwise[i].size();
  
  gsVector<unsigned> v(d);
  v.setZero();

  res.resize( d, np.prod() );

  unsigned r = 0;
  do {
    // Make tensor product of sample points
    for (unsigned i=0; i<d; ++i )
        res(i,r) =  (cwise[i])( v[i] )  ;
    ++r ;
  } while (nextLexicographic(v, np));
}


template<class T>
void tensorProduct( std::vector< gsVector<T>* > const & cwise, gsVector<T>& res)
{
  unsigned d = cwise.size();

  gsVector<unsigned> np(d);
  for (unsigned i=0; i<d; ++i )
    np[i] = cwise[i]->size();
  
  gsVector<unsigned> v(d);
  v.setZero();

  res.resize( np.prod() );

  unsigned r = 0;
  do {
    // Compute tensor product of vectors
    res[r] = (*cwise[0])[v[0]];
    for (unsigned i=1; i<d; ++i)
      res[r] *= (*cwise[i])[v[i]];
    ++r ;
  } while (nextLexicographic(v, np));
}


template<typename T>
gsVector<unsigned> uniformSampleCount (const gsVector<T>& lower, 
                                       const gsVector<T>& upper, 
                                       int numPoints)
{
    const index_t d = lower.rows();
    assert( d == upper.rows() );

    // TO do : phys. volume criterion
    gsVector<T> span = upper - lower;
    const T volume = span.prod();
    const T h = math::pow(volume / numPoints, 1.0 / d);

    gsVector<unsigned> np(d);

    for (index_t i = 0; i < d; ++i)
    {
        np[i] = cast<T,unsigned>(math::ceil( span[i] / h ) );
        GISMO_ASSERT( np[i] > 0, "Something went wrong, number of points is zero..");
    }

    return np;
}

template<typename T>
typename gsMatrix<T>::uPtr uniformPointGrid(const gsVector<T>& lower, 
                                            const gsVector<T>& upper, 
                                            int numPoints)
{
    const gsVector<unsigned> cwisePoints = uniformSampleCount(lower, upper, numPoints);
    
    return gsPointGrid(lower, upper, cwisePoints);
}

template<class T>
void uniformPointGrid( gsMatrix<T> const        & box,
                       index_t                    np,
                       gsMatrix<T>              & result
                        )
{
    const gsVector<T> low = box.col(0);
    const gsVector<T> upp = box.col(1);
    const gsVector<unsigned> npts = uniformSampleCount(low, upp, np);
    result = gsPointGrid(low, upp, npts);
}

template<typename T>
void uniformIntervals(const gsVector<T>& lower, 
                      const gsVector<T>& upper, 
                      std::vector< std::vector<T> >& intervals, 
                      int numIntervals)
{
    const int d = lower.rows();
    assert( d == upper.rows() );

    gsVector<unsigned> np = uniformSampleCount( lower, upper, numIntervals );

    // resize to dimension d without copying old contents, if any
    intervals.clear();
    intervals.resize(d);

    for (int i = 0; i < d; ++i)
    {
        int numInt = np[i] - 1;
        if (numInt <= 1)
            numInt = 1;

        const T h = T(1) / numInt;

        intervals[i].resize(numInt + 1);
        for (int j = 0; j <= numInt; ++j)
            intervals[i][j] = j * h;
    }
}


};// namespace gismo
