/** @file gsBSplineBasis.hpp

    @brief Implementation of 1D B-spline basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mokris
*/

#pragma once

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsBSplineAlgorithms.h>

#include <gsNurbs/gsDeboor.hpp>
#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsUtils/gsMesh/gsMesh.h>

#include <gsMatrix/gsSparseRows.hpp>

#include <gsIO/gsXml.h>

namespace gismo
{


template <class T>
typename gsTensorBSplineBasis<1,T>::Self_t *
gsTensorBSplineBasis<1,T>::New(std::vector<gsBasis<T>*> & bb )
{
    GISMO_ASSERT( bb.size() == 1, "Expecting one component");
    Self_t * c = dynamic_cast<Self_t*>(bb.front());
    if ( NULL != c )
    {
        bb.clear();
        return c;
    }
    else
    {
        gsWarn<<"Something went wrong during BSpline construction.\n";
        return NULL;
    }
}

template <class T>
size_t gsTensorBSplineBasis<1,T>::elementIndex(const gsVector<T> & u ) const
{
    return m_knots.uFind(u(0,0)).uIndex();
}

template <class T>
size_t gsTensorBSplineBasis<1,T>::elementIndex(T u ) const
{
    return m_knots.uFind( u ).uIndex();
}

template <class T>
void gsTensorBSplineBasis<1,T>::connectivity(const gsMatrix<T> & nodes,
                                             gsMesh<T> & mesh) const
{
    const index_t sz  = size();
    GISMO_ASSERT( nodes.rows() == sz, "Invalid input.");

    // Add first vertex
    if ( sz != 0 )
        mesh.addVertex( nodes.row(0).transpose() );

    // Add adges and vertices
    for(index_t i = 1; i< sz; ++i )
    {
        mesh.addVertex( nodes.row(i).transpose() );
        mesh.addEdge( i-1, i );
    }

    if ( m_periodic )
        mesh.addEdge( sz-1, 0);
}

template <class T>
void gsTensorBSplineBasis<1,T>::matchWith(const boundaryInterface & bi,
                                          const gsBasis<T> & other,
                                          gsMatrix<index_t> & bndThis,
                                          gsMatrix<index_t> & bndOther) const
{
    if ( const TensorSelf_t * _other = dynamic_cast<const TensorSelf_t*>(&other) )
    {
        bndThis .resize(1,1);
        bndOther.resize(1,1);
        bndThis (0,0) =  bi.first() .side() == 1 ? 0 :         m_knots.size()        -m_p-2;
        bndOther(0,0) =  bi.second().side() == 1 ? 0 : _other->m_knots.size()-_other->m_p-2;
        return;
    }

    gsWarn<< "Cannot match with "<< other <<"\n";
}

template <class T>
void gsTensorBSplineBasis<1,T>::active_into(const gsMatrix<T>& u,
                                            gsMatrix<index_t>& result ) const
{
    result.resize(m_p+1, u.cols());

    if ( m_periodic )
    {
        // We want to keep the non-periodic case unaffected wrt
        // complexity, therefore we keep the modulo operation of the
        // periodic case separate
        const index_t s = size();
        for (index_t j = 0; j < u.cols(); ++j)
        {
            unsigned first = firstActive(u(0,j));
            for (int i = 0; i != m_p+1; ++i)
                result(i,j) = (first++) % s;
        }
    }
    else
    {
        for (index_t j = 0; j < u.cols(); ++j)
        {
            unsigned first = firstActive(u(0,j));
            for (int i = 0; i != m_p+1; ++i)
                result(i,j) = first++;
        }
    }
}

template <class T>
bool gsTensorBSplineBasis<1,T>::isActive(const index_t i, const gsVector<T>& u) const
{
    GISMO_ASSERT( u.rows() == 1, "Invalid input.");
    // Note: right end of the support will be considered active
    return( (u.value() >= m_knots[i]) && (u.value() <= m_knots[i+m_p+1]) );
}

template <class T>
gsMatrix<index_t> gsTensorBSplineBasis<1,T>::allBoundary() const
{
    if( m_periodic ) // Periodic basis does not have such things as boundaries.
    {
        gsWarn << "Periodic basis does not have such things as boundaries.\n";
        // return NULL;
        gsMatrix<index_t> matrix;
        return matrix;
    }
    else
    {
        gsMatrix<index_t> res(2,1);
        res(0,0) = 0;
        res(1,0) = m_knots.size()-m_p-2;
        return res;
    }
}


template <class T>
gsMatrix<index_t> gsTensorBSplineBasis<1,T>::boundaryOffset(boxSide const & s,
                                                               index_t offset ) const
{
    if( m_periodic )
        gsWarn << "Periodic basis does not have such things as boundaries.\n";

    gsMatrix<index_t> res(1,1);
    GISMO_ASSERT(offset+m_p+1 < static_cast<index_t>(m_knots.size()),
                 "Offset cannot be bigger than the amount of basis functions orthogonal to Boxside s!");
    switch (s) {
    case boundary::left : // left
        res(0,0)= offset;
        break;
    case boundary::right : // right
        res(0,0)= m_knots.size()-m_p-2-offset;
        break;
    default:
        GISMO_ERROR("gsBSplineBasis: valid sides is left(west) and right(east).");
    };
    return res;
}

/// @cond
template <class T>
gsConstantBasis<T> * gsTensorBSplineBasis<1,T>::boundaryBasis_impl(boxSide const &) const
{
    return new gsConstantBasis<T>(1.0);
}
/// @endcond

template <class T>
gsMatrix<T> gsTensorBSplineBasis<1,T>::support() const
{
    // The support of a the whole B-spline basis is the interval
    // between m_knots[m_p] and m_knots[i+m_p+1].
    // First m_p and last m_p knots are the "ghost" knots and are not
    // part of the domain
    gsMatrix<T> res(1,2);
    res << domainStart() , domainEnd() ;
    return res ;
}

template <class T>
gsMatrix<T> gsTensorBSplineBasis<1,T>::support(const index_t & i) const
{
    // Note: in the periodic case last index is
    // m_knots.size() - m_p - 1 - m_periodic
    // but we may still accept the original index of the basis function.
    // One issue is that in the periodic case the basis function
    // support has two connected components, so probably one should
    // call this function with twice for the two twins

    GISMO_ASSERT( static_cast<size_t>(i) < m_knots.size()-m_p-1,
                  "Invalid index of basis function." );
    gsMatrix<T> res(1,2);
    res << ( i > m_p ? m_knots[i] : m_knots[m_p] ),
        ( static_cast<size_t>(i) < (m_knots.size()-2*m_p-2) ? m_knots[i+m_p+1] :
          m_knots[m_knots.size()-m_p-1] );
    return res ;
}

template <class T>
index_t gsTensorBSplineBasis<1,T>::twin(index_t i) const
{
    if( m_periodic == 0 )
        return i;
    const index_t s = size();
    if ( i < static_cast<index_t>(m_periodic) )
        i += s;
    else if ( i > s )
        i -= s;
    return i;
}

template <class T>
void gsTensorBSplineBasis<1,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    result.resize(m_p+1, u.cols() );

#if (FALSE)

    typename KnotVectorType::const_iterator kspan;

    // Low degree specializations
    switch (m_p)
    {
    case 0: // constant B-splines
        result.setOnes();
        return;

    case 1:
        for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
            if ( ! inDomain( u(0,v) ) )
                result.col(v).setZero();
            else
            {
                kspan = m_knots.findspanIter ( u(0,v) );
                bspline::evalDeg1Basis( u(0,v), kspan, m_p, result.col(v) );
            }
        return;
    case 2:
        for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
            if ( ! inDomain( u(0,v) ) )
                result.col(v).setZero();
            else
            {
                kspan = m_knots.findspanIter ( u(0,v) );
                bspline::evalDeg2Basis( u(0,v), kspan, m_p, result.col(v) );
            }
        return;
    case 3:
        for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
            if ( ! inDomain( u(0,v) ) )
                result.col(v).setZero();
            else
            {
                kspan = m_knots.findspanIter ( u(0,v) );
                bspline::evalDeg3Basis( u(0,v), kspan, m_p, result.col(v) );
            }
        return;
    };

    STACK_ARRAY(T, left, m_p + 1);
    STACK_ARRAY(T, right, m_p + 1);

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
    {
        // Check if the point is in the domain
        if ( ! inDomain( u(0,v) ) )
        {
            // gsWarn<< "Point "<< u(0,v) <<" not in the BSpline domain.\n";
            result.col(v).setZero();
        }
        else
        {   // Locate the point in the knot-vector
            kspan = m_knots.findspanIter ( u(0,v) );
            // Run evaluation algorithm
            bspline::evalBasis( u(0,v), kspan, m_p, left, right, result.col(v) );
        }
    }

#endif
//#if (FALSE)
    STACK_ARRAY(T, left, m_p + 1);
    STACK_ARRAY(T, right, m_p + 1);

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
    {
        // Check if the point is in the domain
        if ( ! inDomain( u(0,v) ) )
        {
            // gsWarn<< "Point "<< u(0,v) <<" not in the BSpline domain.\n";
            result.col(v).setZero();
            continue;
        }

        // Run evaluation algorithm

        // Get span of absissae
        unsigned span = m_knots.iFind( u(0,v) ) - m_knots.begin() ;

        //ndu[0]   = T(1);  // 0-th degree function value
        result(0,v)= T(1);  // 0-th degree function value

        for(int j=1; j<= m_p; ++j) // For all degrees ( ndu column)
        {
            left[j]  = u(0,v) - m_knots[span+1-j];
            right[j] = m_knots[span+j] - u(0,v);
            T saved = T(0) ;

            for(int r=0; r!=j ; ++r) // For all (except the last)  basis functions of degree j ( ndu row)
            {
                // Strictly lower triangular part: Knot differences of distance j
                //ndu[j*p1 + r] = right[r+1]+left[j-r] ;

                //const T temp = ndu[r*p1 + j-1] / ndu[j*p1 + r] ;
                const T temp = result(r,v) / ( right[r+1] + left[j-r] );
                // Upper triangular part: Basis functions of degree j
                //ndu[r*p1 + j] = saved + right[r+1] * temp ;// r-th function value of degree j
                result(r,v)     = saved + right[r+1] * temp ;// r-th function value of degree j
                saved = left[j-r] * temp ;
            }
            //ndu[j*p1 + j] = saved ;// Diagonal: j-th (last) function value of degree j
            result(j,v)     = saved;
        }

    }// end for all columns v
//#endif
}


template <class T>
void gsTensorBSplineBasis<1,T>::evalSingle_into(index_t i,
                                                const gsMatrix<T> & u,
                                                gsMatrix<T>& result) const
{
    GISMO_ASSERT( static_cast<size_t>(i) < m_knots.size()-m_p-1,"Invalid index of basis function." );

    result.resize(1, u.cols() );
    STACK_ARRAY(T, N, m_p + 1);

    for (index_t s=0;s<u.cols(); ++s)
    {
        //Special cases

        // Periodic basis ?
        // Note: for periodic, i is understood modulo the basis size
        if( m_periodic )
        {
            // Basis is periodic and i could be among the crossing functions.
            if( static_cast<int>(i) <= m_periodic )
            {
                gsMatrix<T> supp = support(i);
                if( u(0,s) < supp(0) || u(0,s) > supp(1)) // u is outside the support of B_i
                    i += size(); // we use the basis function ignored ny the periodic construction
            }
        }

        // Special case of C^{-1} on right end of support
        if ( (static_cast<size_t>(i) == m_knots.size()-m_p-2) &&
             (u(0,s) == m_knots.last()) &&  (u(0,s)== m_knots[m_knots.size()-m_p-1]) )
        {
            result(0,s)= T(1.0);
            continue;
        }

        // Locality property
        if ( (u(0,s) < m_knots[i]) || (u(0,s) >= m_knots[i+m_p+1]) )
        {
            result(0,s)= T(0.0);
            continue;
        }

        //bspline::deBoorTriangle( u(0,s), m_knots.begin() + i, m_p, N );
        //result(0,s) = N[ m_p ];

        // Initialize zeroth degree functions
        for (int j=0;j<=m_p; ++j)
            if ( u(0,s) >= m_knots[i+j] && u(0,s) < m_knots[i+j+1] )
                N[j] = T(1.0);
            else
                N[j] = T(0.0);
        // Compute according to the trangular table
        for (int k=1;k<=m_p; ++k)
        {
            T saved;
            if (N[0]==0)
                saved= 0.0;
            else
                saved= ((u(0,s) - m_knots[i] )* N[0]) / (m_knots[k+i] - m_knots[i]);
            for (int j=0;j<m_p-k+1; ++j)
            {
                const T kleft  = m_knots[i+j+1];
                const T kright = m_knots[i+j+k+1];
                if ( N[j+1] == 0.0 )
                {
                    N[j] = saved;
                    saved = 0.0 ;
                }
                else
                {
                    const T temp = N[j+1] / ( kright - kleft );
                    N[j]= saved + ( kright-u(0,s) )*temp;
                    saved= (u(0,s) -kleft ) * temp;
                }
            }
        }
        result(0,s)= N[0];
    }
}

template <class T>
void gsTensorBSplineBasis<1,T>::evalDerSingle_into(index_t i,
                                                   const gsMatrix<T> & u,
                                                   int n,
                                                   gsMatrix<T>& result) const
{
    GISMO_ASSERT( u.rows() == 1 , "gsBSplineBasis accepts points with one coordinate.");
    // This is a copy of the interior of evalAllDerSingle_into() except for the last for cycle with k.

    // Notation from the NURBS book:
    // p - degree
    // U = {u[0],...,u[m]} the knot vector
    // i - id of the basis function (hopefully the same ;)
    // u - parameter (where to evaluate)

    // Note that N[i][k] = N[ i*table_size + k ].

    T saved;
    T Uleft;
    T Uright;
    T temp;
    const unsigned int table_size = m_p + 1;
    STACK_ARRAY(T, N, table_size*table_size);
    STACK_ARRAY(T, ND, table_size);

    result.resize(1, u.cols()); // (1 derivative) times number of points of evaluation.

    // From here on, we follow the book.
    for( index_t s = 0; s < u.cols(); s++ )
    {
        // Special cases
        // Periodic basis
        if( (int)(i) <= m_periodic ) // Basis is periodic and i could be among the crossing functions.
        {
            gsMatrix<T> supp = support(i);
            if( u(0,s) < supp(0) || u(0,s) > supp(1)) // u is outside the support of B_i
                i += size(); // we use the basis function ignored ny the periodic construction
        }

        if( u(0,s) == m_knots.last() )
            gsWarn << "evalDerSingle_into sometimes gives strange results for the last knot.\n";

        if( u(0,s) < m_knots[i] || u(0,s) >= m_knots[i+m_p+1]) // Local property
        {
            // If we are outside the support, the value and all the
            // derivatives there are equal to zero.
            result(0,s ) = 0;
            continue;// to the next point
        }

        bspline::deBoorTriangle( u(0,s), m_knots.begin() + i, m_p, N );

        /*
          for( int j = 0; j <= m_p; j++ ) // Initialize zeroth-degree functions
          if( u(0,s) >= m_knots[i+j] && u(0,s) < m_knots[i+j+1] )
          N[ j*table_size ] = 1;
          else
          N[ j*table_size ] = 0;
          for( int k = 1; k <= m_p; k++ ) // Compute full triangular table
          {
          if( N[(k-1)] == 0 )
          saved = 0;
          else
          saved = ((u(0,s)-m_knots[i])*N[ k-1 ])/(m_knots[i+k]-m_knots[i]);
          for( int j = 0; j < m_p-k+1; j++)
          {
          Uleft = m_knots[i+j+1];
          Uright = m_knots[i+j+k+1];
          if( N[ (j+1)*table_size + (k-1) ] == 0 )
          {
          N[ j*table_size + k ] = saved;
          saved = 0;
          }
          else
          {
          temp = N[ (j+1)*table_size + (k-1) ]/(Uright-Uleft);
          N[ j*table_size + k ] = saved + (Uright-u(0,s))*temp;
          saved = (u(0,s)-Uleft)*temp;
          }
          }
          }*/

        // Not tested few lines:
        if( n == 0 )
        {
            result(0,s) = N[ m_p ]; // The function value
            continue;
        }

        for( int j = 0; j <= n; j++ ) // Load appropriate column
            ND[j] = N[ j*table_size + (m_p-n) ];
        for( int jj = 1; jj <= n; jj++ ) // Compute table of width n
        {
            if( ND[0] == 0 )
                saved = 0;
            else
                saved = ND[0]/(m_knots[ i+m_p-n+jj ] - m_knots[i]);
            for( int j = 0; j < n-jj+1; j++ )
            {
                Uleft = m_knots[i+j+1];
                Uright = m_knots[i+j+m_p-n+jj+1];
                if( ND[j+1] == 0 )
                {
                    ND[j] = (m_p-n+jj)*saved;
                    saved = 0;
                }
                else
                {
                    temp = ND[j+1]/(Uright - Uleft);
                    ND[j] = (m_p-n+jj)*(saved - temp);
                    saved = temp;
                }
            }
        }
        result(0,s) = ND[0]; // n-th derivative
    }
    // Changes from the book: notation change, stack arrays instead of normal C arrays, changed all 0.0s and 1.0s into 0 and 1. No return in the first check
    // Warning: Segmentation fault if n > m_p.

}

template <class T> inline
void gsTensorBSplineBasis<1,T>::evalFunc_into(const gsMatrix<T> &u,
                                              const gsMatrix<T> & coefs,
                                              gsMatrix<T>& result) const
{
    GISMO_ASSERT( u.rows() == 1 , "gsBSplineBasis accepts points with one coordinate (got "
                  <<u.rows()<<").");
    if( m_periodic == 0 )
        gsDeboor(u, m_knots, m_p, coefs, result); // De Boor's algorithm
    else
    {
        gsDeboor(u, m_knots, m_p, perCoefs(coefs), result); // Make sure that the ghost coefficients agree with the beginning of the curve.
    }
    //  return gsBasis<T>::eval(u,coefs); // defalult gsBasis implementation
}


template <class T> inline
void gsTensorBSplineBasis<1,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    GISMO_ASSERT( u.rows() == 1 , "gsBSplineBasis accepts points with one coordinate.");

    const int pk = m_p-1 ;
    const int p1 = m_p + 1;       // degree plus one
    STACK_ARRAY(T, ndu  , m_p);
    STACK_ARRAY(T, left , p1 );
    STACK_ARRAY(T, right, p1 );

    result.resize( m_p + 1, u.cols() ) ;

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
    {
        // Check if the point is in the domain
        if ( ! inDomain( u(0,v) ) )
        {
            // gsWarn<< "Point "<< u(0,v) <<" not in the BSpline domain.\n";
            result.col(v).setZero();
            continue;
        }

        // Run evaluation algorithm and keep first derivative

        // Get span of absissae
        typename KnotVectorType::iterator span = m_knots.iFind( u(0,v) );

        ndu[0]  = T(1); // 0-th degree function value
        left[0] = 0;

        for(int j=1; j<m_p ; j++) // For all degrees
        {
            // Compute knot splits
            left[j]  = u(0,v) - *(span+1-j);
            right[j] = *(span+j) - u(0,v);

            // Compute Basis functions of degree m_p-1 ( ndu[] )
            T saved = T(0) ;
            for(int r=0; r<j ; r++) // For all (except the last) basis functions of degree 1..j
            {
                const T temp = ndu[r] / ( right[r+1] + left[j-r] ) ;
                ndu[r] = saved + right[r+1] * temp ;// r-th function value of degree 1..j
                saved = left[j-r] * temp ;
            }
            ndu[j] = saved ;// last function value of degree 1..j
        }

        // Compute last knot split
        left[m_p]  = u(0,v) - *(span+1-m_p);
        right[m_p] = *(span+m_p) - u(0,v);

        // Compute the first derivatives (using ndu[] and left+right)
        right[0] = right[1]+left[m_p] ;
        result(0  , v) = - m_p * ndu[0]  /  right[0] ;
        for(int r = 1; r < m_p; r++)
        {
            // Compute knot difference r of distance m_p (overwrite right[])
            right[r] = right[r+1]+left[m_p-r] ;
            result(r, v)  = m_p * (  ndu[r-1] / right[r-1] - ndu[r]  /  right[r] );
        }
        result(m_p, v) =   m_p * ndu[pk] / right[pk];

    }// end for all columns v
}

template <class T> inline
void gsTensorBSplineBasis<1,T>::deriv2_into(const gsMatrix<T> & u,
                                            gsMatrix<T>& result ) const
{
    std::vector<gsMatrix<T> > ev;
    this->evalAllDers_into(u, 2, ev);
    result.swap(ev[2]);
}

template <class T>  inline
void gsTensorBSplineBasis<1,T>::derivSingle_into(index_t i,
                                                 const gsMatrix<T> & u,
                                                 gsMatrix<T>& result ) const
{
    // \todo Redo an efficient implementation p. 76, Alg. A2.5 Nurbs book
    result.resize(1, u.cols() );
    gsMatrix<T> tmp;
    gsTensorBSplineBasis<1,T>::deriv_into(u, tmp);

    for (index_t j = 0; j < u.cols(); ++j)
    {
        const index_t first = firstActive(u(0,j));
        if ( (i>= first) && (i<= first + m_p) )
            result(0,j) = tmp(i-first,j);
        else
            result(0,j) = T(0.0);
    }
}

template <class T>  inline
void
gsTensorBSplineBasis<1,T>::evalAllDersSingle_into(index_t i,
                                                  const gsMatrix<T> & u,
                                                  int n,
                                                  gsMatrix<T>& result) const
{
    GISMO_ASSERT( u.rows() == 1 , "gsBSplineBasis accepts points with one coordinate.");
    //gsWarn << "You're about to use evalAllDersSingle_into(...) that has not been tested at all.\n";

    // Notation from the NURBS book:
    // p - degree
    // U = {u[0],...,u[m]} the knot vector
    // i - id of the basis function (hopefully the same ;)
    // u - parameter (where to evaluate)

    // Note that N[i][k] = N[ i*table_size + k ].

    T saved;
    T Uleft;
    T Uright;
    T temp;
    const unsigned int table_size = m_p + 1;
    STACK_ARRAY(T, N, table_size*table_size);
    STACK_ARRAY(T, ND, table_size);

    result.resize(( n + 1 ), u.cols()); // (1 value + n derivatives) times number of points of evaluation.

    gsWarn << "evalAllDersSingle_into(...) has a bug in the second derivatives.\n";

    // Commented parts were attempts to solve the bug in the second derivatives at the right endpoint.



    for( index_t s = 0; s < u.cols(); s++ )
    {
        // Periodic basis
        if( (int)(i) <= m_periodic ) // Basis is periodic and i could be among the crossing functions.
        {
            gsMatrix<T> supp = support(i);
            if( u(0,s) < supp(0) || u(0,s) > supp(1)) // u is outside the support of B_i
                i += size(); // we use the basis function ignored ny the periodic construction
        }
        if( u(0,s) < m_knots[i] || u(0,s) > m_knots[i+m_p+1]) // Local property
        {
            if( u( 0,s ) != m_knots[ m_knots.size() - 1 ] )
            {
                for( int j = 0; j <= n; j++ )
                    result( j,s ) = 0.0;
                continue;
            }
        }

        bspline::deBoorTriangle( u(0,s), m_knots.begin() + i, m_p, N );

        // The following commented area has been replaced by the deBoorTriangle(...); line.
/*        for( int j = 0; j <= m_p; j++ ) // Initialize zeroth-degree functions
          if( u(0,s) >= m_knots[i+j] && u(0,s) < m_knots[i+j+1] )
          N[ j*table_size ] = 1;
          else
          N[ j*table_size ] = 0;
          for( int k = 1; k <= m_p; k++ ) // Compute full triangular table
          {
          if( N[(k-1)] == 0 )
          saved = 0;
          else
          saved = ((u(0,s)-m_knots[i])*N[ k-1 ])/(m_knots[i+k]-m_knots[i]);
          for( int j = 0; j < m_p-k+1; j++)
          {
          Uleft = m_knots[i+j+1];
          Uright = m_knots[i+j+k+1];
          if( N[ (j+1)*table_size + (k-1) ] == 0 )
          {
          N[ j*table_size + k ] = saved;
          saved = 0;
          }
          else
          {
          temp = N[ (j+1)*table_size + (k-1) ]/(Uright-Uleft);
          N[ j*table_size + k ] = saved + (Uright-u(0,s))*temp;
          saved = (u(0,s)-Uleft)*temp;
          }
          }
          }
*/
        result(0,s) = N[ m_p ]; // The function value
        for( int k = 1; k <= n; k++ ) // Compute the derivatives
        {
            for( int j = 0; j <= k; j++ ) // Load appropriate column
                ND[j] = N[ j*table_size + (m_p-k) ];
            for( int jj = 1; jj <= k; jj++ ) // Compute table of width k
            {
                if( ND[0] == 0 )
                    saved = 0;
                else
                    saved = ND[0]/(m_knots[ i+m_p-k+jj ] - m_knots[i]);
                for( int j = 0; j < k-jj+1; j++ )
                {
                    Uleft = m_knots[i+j+1];
                    Uright = m_knots[i+j+m_p-k+jj+1];
                    if( ND[j+1] == 0 )
                    {
                        ND[j] = (m_p-k+jj)*saved;
                        saved = 0;
                    }
                    else
                    {
                        temp = ND[j+1]/(Uright - Uleft);
                        ND[j] = (m_p-k+jj)*(saved - temp);
                        saved = temp;
                    }
                }
            }
            result(k,s) = ND[0]; // k-th derivative
        }
    }

    /*
      if( clamped_end )
      {
      m_knots[ m_knots.size() - 1 ] -= 0.001;
      m_knots[ m_knots.size() - 2 ] -= 0.001;
      }
      // Changes from the book: notation change, stack arrays instead of normal C arrays, changed all 0.0s and 1.0s into 0 and 1. No return in the first check
      // Warning: Segmentation fault if n > m_p.
      gsDebug << "After:" << std::endl;
      for( int j = 0; j < m_knots.size(); j++ )
      {
      gsDebug << "m_knots[" << j << "] =" << m_knots[j] << std::endl;
      }*/
}


template <class T> inline
void gsTensorBSplineBasis<1,T>::deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const
{
    // TO DO specialized computation for gsBSplineBasis
    if( m_periodic == 0 )
        gsBasis<T>::derivFunc_into(u, coefs, result);
    else
        gsBasis<T>::derivFunc_into(u, perCoefs(coefs), result);
}

template <class T> inline
void gsTensorBSplineBasis<1,T>::deriv2_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const
{
    // TO DO specialized computation for gsBSplineBasis
    if( m_periodic == 0 )
        gsBasis<T>::deriv2Func_into(u, coefs, result);
    else
        gsBasis<T>::deriv2Func_into(u, perCoefs(coefs), result);
}

template <class T>  inline
void gsTensorBSplineBasis<1,T>::deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    // \todo Redo an efficient implementation p. 76, Alg. A2.5 Nurbs book
    result.resize(1, u.cols() );
    gsMatrix<T> tmp;
    gsTensorBSplineBasis<1,T>::deriv2_into(u, tmp);

    for (index_t j = 0; j < u.cols(); ++j)
    {
        const index_t first = firstActive(u(0,j));
        if ( (i>= first) && (i<= first + m_p) )
            result(0,j) = tmp(i-first,j);
        else
            result(0,j) = T(0.0);
    }
}

template <class T> inline
gsMatrix<T> gsTensorBSplineBasis<1,T>::laplacian(const gsMatrix<T> & u ) const
{
    std::vector<gsMatrix<T> > ev;
    this->evalAllDers_into(u, 2, ev);
    return ev[2].colwise().sum();
}

template <class T>
typename gsBasis<T>::uPtr
gsTensorBSplineBasis<1,T>::tensorize(const gsBasis<T> & other) const
{
    typename Self_t::uPtr ptr1 = memory::convert_ptr<Self_t>(other.clone());

    if ( ptr1 )
        return typename gsBasis<T>::uPtr(new gsTensorBSplineBasis<2,T>( ptr1.release(), this->clone().release() ));
    else
    {
        gsWarn<<"gsTensorBSplineBasis::tensorize: Invalid basis "<< other <<"\n";
        return typename gsBasis<T>::uPtr();
    }
}

template <class T>
memory::unique_ptr<gsGeometry<T> > gsBSplineBasis<T>::makeGeometry( gsMatrix<T> coefs ) const
{
    return typename gsGeometry<T>::uPtr(new GeometryType(*this, give(coefs)));
}

template <class T>
void gsTensorBSplineBasis<1,T>::
evalAllDers_into(const gsMatrix<T> & u, int n,
                 std::vector<gsMatrix<T> >& result) const
{
    // TO DO : Use less memory proportionally to n
    // Only last n+1 columns and last n rows of ndu are needed
    // Also a's size is proportional to n
    GISMO_ASSERT( u.rows() == 1 , "gsBSplineBasis accepts points with one coordinate.");

    const int p1 = m_p + 1;       // degree plus one

    STACK_ARRAY(T, ndu,  p1 * p1 );
    STACK_ARRAY(T, left, p1);
    STACK_ARRAY(T, right, p1);
    STACK_ARRAY(T, a, 2 * p1);

    result.resize(n+1);
    for(int k=0; k<=n; k++)
        result[k].resize(m_p + 1, u.cols());

#if FALSE

    const int pn = m_p - n;

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
    {
        // Check if the point is in the domain
        if ( ! inDomain( u(0,v) ) )
        {
            // gsWarn<< "Point "<< u(0,v) <<" not in the BSpline domain.\n";
            for(int k=0; k<=n; k++)
                result[k].col(v).setZero();
            continue;
        }

        // Run evaluation algorithm and keep the function values triangle & the knot differences
        unsigned span = m_knots.findspan( u(0,v) ) ;     // Get span of absissae

        for(int j=1; j<= m_p; j++) // For all degrees ( ndu column)
        {
            // Compute knot splits
            left[j]  = u(0,v) - m_knots[span+1-j];
            right[j] = m_knots[span+j] - u(0,v);
        }

        ndu[0] = T(1) ; // 0-th degree function value
        T saved = T(0) ;
        // Compute Basis functions of degree m_p-n ( ndu[] )
        for(int r=0; r<pn ; r++) // For all (except the last) basis functions of degree m_p-n
        {
            const T temp = ndu[r*p1 + pn -1] / ( right[r+1] + left[pn-r] ) ;
            ndu[r*p1 + pn] = saved + right[r+1] * temp ;// r-th function value of degree j
            saved = left[pn-r] * temp ;
        }
        ndu[pn*p1 + pn] = saved ;

        // Compute n-th derivative and continue to n-1 ... 0

        for(int r=0; r<pn ; r++) // For all (except the last)  basis functions of degree j ( ndu row)
        {
            // Strictly lower triangular part: Knot differences of distance j
            ndu[j*p1 + r] = right[r+1]+left[j-r] ;
            const T temp = ndu[r*p1 + j-1] / ndu[j*p1 + r] ;
            // Upper triangular part: Basis functions of degree j
            ndu[r*p1 + j] = saved + right[r+1] * temp ;// r-th function value of degree j
            saved = left[j-r] * temp ;
        }
        // Diagonal: j-th (last) function value of degree j
        ndu[j*p1 + j] = saved ;
    }

#endif

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
    {
        // Check if the point is in the domain
        if ( ! inDomain( u(0,v) ) )
        {
            //gsDebug<< "Point "<< u(0,v) <<" not in the BSpline domain ["
            //      << *(m_knots.begin()+m_p)<< ", "<<*(m_knots.end()-m_p-1)<<"].\n";
            for(int k=0; k<=n; k++)
                result[k].col(v).setZero();
            continue;
        }

        // Run evaluation algorithm and keep the function values triangle & the knot differences
        typename KnotVectorType::iterator span = m_knots.iFind( u(0,v) );

        ndu[0] = T(1) ; // 0-th degree function value
        for(int j=1; j<= m_p; j++) // For all degrees ( ndu column)
        {
            // Compute knot splits
            left[j] = u(0,v) - *(span+1-j);
            right[j] = *(span+j) - u(0,v);

            T saved = T(0) ;

            for(int r=0; r<j ; r++) // For all (except the last)  basis functions of degree j ( ndu row)
            {
                // Strictly lower triangular part: Knot differences of distance j
                ndu[j*p1 + r] = right[r+1]+left[j-r] ;
                const T temp = ndu[r*p1 + j-1] / ndu[j*p1 + r] ;
                // Upper triangular part: Basis functions of degree j
                ndu[r*p1 + j] = saved + right[r+1] * temp ;// r-th function value of degree j
                saved = left[j-r] * temp ;
            }
            // Diagonal: j-th (last) function value of degree j
            ndu[j*p1 + j] = saved ;
        }

        // Assign 0-derivative equal to function values
        //result.front().block(0,v, p1,1) = ndu.col(m_p);
        for (int j=0; j <= m_p ; ++j )
            result.front()(j,v) = ndu[j*p1 + m_p];

        // Compute the derivatives
        for(int r = 0; r <= m_p; r++)
        {
            // alternate rows in array a
            T* a1 = &a[0];
            T* a2 = &a[p1];

            a1[0] = T(1) ;

            // Compute the k-th derivative of the r-th basis function
            for(int k=1; k<=n; k++)
            {
                int rk,pk,j1,j2 ;
                T d(0) ;
                rk = r-k ; pk = m_p-k ;

                if(r >= k)
                {
                    a2[0] = a1[0] / ndu[ (pk+1)*p1 + rk] ;
                    d = a2[0] * ndu[rk*p1 + pk] ;
                }

                j1 = ( rk >= -1  ? 1   : -rk     );
                j2 = ( r-1 <= pk ? k-1 : m_p - r );

                for(int j = j1; j <= j2; j++)
                {
                    a2[j] = (a1[j] - a1[j-1]) / ndu[(pk+1)*p1 + rk+j] ;
                    d += a2[j] * ndu[(rk+j)*p1 + pk] ;
                }

                if(r <= pk)
                {
                    a2[k] = -a1[k-1] / ndu[(pk+1)*p1 + r] ;
                    d += a2[k] * ndu[r*p1 + pk] ;
                }

                result[k](r, v) = d;

                std::swap(a1, a2);              // Switch rows
            }
        }
    }// end for all columns v

    // Multiply through by the factor factorial(m_p)/factorial(m_p-k)
    int r = m_p ;
    for(int k=1; k<=n; k++)
    {
        result[k].array() *= T(r) ;
        r *= m_p - k ;
    }
}


template <class T>
void gsTensorBSplineBasis<1,T>::refine_withCoefs(gsMatrix<T>& coefs, const std::vector<T>& knots)
{
    // Probably does not work with periodic basis. We would need to shift the knots in the ghost areas by the active length,
    // sort the knot vector, swap the coefs accordingly, call this with perCoefs (or coefs?) and stretch outer knots.
    gsBoehmRefine(this->knots(), coefs, m_p, knots.begin(), knots.end());
}


template <class T>
void gsTensorBSplineBasis<1,T>::refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector<T>& knots)
{
    // See remark about periodic basis in refine_withCoefs, please.
    gsSparseRows<T> trans;
    trans.setIdentity( this->size() );
    gsBoehmRefine(this->knots(), trans, m_p, knots.begin(), knots.end());
    trans.toSparseMatrix( transfer );
}


template <class T>
void gsTensorBSplineBasis<1,T>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots, int mul)
{
    // See remark about periodic basis in refine_withCoefs, please.
    std::vector<T> newKnots;
    this->knots().getUniformRefinementKnots(numKnots, newKnots,mul);
    this->refine_withCoefs(coefs, newKnots);
}


template <class T>
void gsTensorBSplineBasis<1,T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots, int mul)
{
    // See remark about periodic basis in refine_withCoefs, please.
    std::vector<T> newKnots;
    this->knots().getUniformRefinementKnots(numKnots, newKnots,mul);
    this->refine_withTransfer(transfer, newKnots);
}

template <class T>
void gsTensorBSplineBasis<1,T>::uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots)
{
    // Simple implementation: coarsen and refine again.
    // Could be done more efficiently if needed.
    std::vector<T> removedKnots = m_knots.coarsen(numKnots);
    this->clone()->refine_withTransfer( transfer, removedKnots );
}

template <class T>
index_t gsTensorBSplineBasis<1,T>::functionAtCorner(boxCorner const & c) const
{
    GISMO_ASSERT(c<3,"Invalid corner for 1D basis.");
    return ( c == 1 ? 0 : this->size()-1);
}

template <class T>
std::ostream & gsTensorBSplineBasis<1,T>::print(std::ostream &os) const
{
    os << "BSplineBasis: deg=" <<degree()
       << ", size="<< size() << ", knot vector:\n";
    os << this->knots();
    if( m_periodic > 0 )
        os << ",\n m_periodic= " << m_periodic;
    return os;
}
template <class T>
std::string gsTensorBSplineBasis<1,T>::detail() const
{
    std::stringstream os;
    os << "BSplineBasis: deg=" <<degree()
       << ", size="<< size() << ", knot vector:\n";
    os << this->knots().detail();
    if( m_periodic > 0 )
        os << ",\n m_periodic= " << m_periodic;
    return os.str();
}


template<class T>
void gsTensorBSplineBasis<1,T>::enforceOuterKnotsPeriodic()
{
    int borderKnotMult = this->borderKnotMult(); // Otherwise complains that error: 'borderKnotMult' cannot be used as a function
    std::vector<T> newKnots(m_knots.begin(), m_knots.end());
    for( int i = 0; i <= m_p - borderKnotMult; i++ )
    {
        newKnots[i] = newKnots[ newKnots.size() - 2 * m_p - 2 + i + borderKnotMult ] - _activeLength();
        newKnots[ newKnots.size() - i - 1 ] = newKnots[ 2 * m_p - borderKnotMult + 1 - i ] + _activeLength();
    }
    m_knots=KnotVectorType(m_p, newKnots.begin(), newKnots.end());
}

template <class T>
void gsTensorBSplineBasis<1,T>::_convertToPeriodic()
{
    gsWarn << "gsBSplineBasis: Converting basis to periodic"<< *this<< "\n";

    int borderKnotMult = this->borderKnotMult();
    if( m_knots.size() < static_cast<size_t>(2 * m_p + 2) ) // We need at least one internal knot span.
    {
        gsWarn << "Your basis cannot be changed into periodic:\n Not enough internal control points for our construction of periodic basis.\n";
        m_periodic = 0;
    }

    else if( isClamped() )
    {
        // Clamped knots are OK, if BOTH of them are clamped.
        _stretchEndKnots();
        m_periodic = m_p + 2 - borderKnotMult; // +2 and not +1 because borderKnotMult has been computed before stretching
    }

    else
    {
        m_periodic = m_p + 1 - borderKnotMult; // We suppose it will be possible to convert the basis to periodic.

        T i1, i2;
        for( int i = 2; i <= 2 * m_p + 1 - borderKnotMult; i ++ )
        {
            // Compare the knot intervals in the beginning with the corresponding knot intervals at the end.
            i1 = m_knots[i] - m_knots[i-1];
            i2 = m_knots[m_knots.size() - (2*m_p) + i - 2 + borderKnotMult] - m_knots[m_knots.size() - (2*m_p) + i - 3 + borderKnotMult];
            if( math::abs( i1 - i2 ) > 1e-8 )
            {
                gsWarn << "Your basis cannot be changed into periodic:\n Trouble stretching interior knots.\n";
                //std::cerr << "i: " << i << ", i1: " << i1 << ", i2: " << i2 << std::endl;
                m_periodic = 0;
                return;
            }
        }
    }
}

template <class T>
int gsTensorBSplineBasis<1,T>::borderKnotMult() const
{
    // Terminology: Blue knots are the m_p + 1st from the beginning and from the end of knot vector.
    if( isClamped() )
        return m_p + 1;

    int multiFirst = 0;
    int multiLast = 0;
    int lastBlueId = m_knots.size() - m_p - 1; // Index of the m_p + 1st knot from the end.

    for( int i = 0; i < m_p; i++ )
    {
        if( m_knots[m_p - i] == m_knots[m_p])
            multiFirst++;
        else
            break;
    }

    for( int i = 0; i < m_p; i++ )
    {
        if( m_knots[lastBlueId + i] == m_knots[ lastBlueId ])
            multiLast++;
        else
            break;
    }

    if( multiFirst == multiLast )
        return multiLast;
    else
    {
        gsWarn << "Different multiplicity of the blue knots.\n";
        return 0;
    }

}

template <class T>
void gsTensorBSplineBasis<1,T>::_stretchEndKnots()
{
    GISMO_ASSERT( isClamped(), "_stretchEndKnots() is intended for use only to knot vectors with clamped end knots.");
    T curFirst=m_knots[0];
    T curLast=m_knots[m_knots.size() - 1];
    m_knots.remove(curFirst); m_knots.insert(m_knots[m_knots.size() - m_p - 2] - _activeLength());
    m_knots.remove(curLast);  m_knots.insert(m_knots[m_p + 1] + _activeLength());
}

/* ********************************************** */

template <class T>
gsBSplineBasis<T> & gsBSplineBasis<T>::component(short_t i)
{
    GISMO_UNUSED(i);
    GISMO_ASSERT(i==0,"gsBSplineBasis has only one component");
    return const_cast<gsBSplineBasis&>(*this);
}

template <class T>
const gsBSplineBasis<T> & gsBSplineBasis<T>::component(short_t i) const
{
    GISMO_UNUSED(i);
    GISMO_ASSERT(i==0,"gsBSplineBasis has only one component");
    return const_cast<gsBSplineBasis&>(*this);
}


namespace internal
{

/// Get a BSplineBasis from XML data
template<class T>
class gsXml< gsBSplineBasis<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBSplineBasis<T>);
    static std::string tag  () { return "Basis"; }
    static std::string type () { return "BSplineBasis"; }

    static gsBSplineBasis<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( !strcmp( node->name(),"Basis"), "Wrong tag." );

        if (!strcmp(node->first_attribute("type")->value(),"TensorBSplineBasis1"))
            node = node->first_node("BSplineBasis");

        GISMO_ASSERT( !strcmp(node->first_attribute("type")->value(),"BSplineBasis"),
                      "Wrong XML type, expected BSplineBasis." );

        gsXmlNode * tmp = node->first_node("KnotVector");
        // if type: == Plain, == Compact ..
        GISMO_ASSERT(tmp, "Did not find a KnotVector tag in the Xml file.");
        gsKnotVector<T> kv;
        gsXml<gsKnotVector<T> >::get_into(tmp, kv);

        return new gsBSplineBasis<T>( kv );
    }

    static gsXmlNode * put (const gsBSplineBasis<T> & obj,
                            gsXmlTree & data)
    {
        // Add a new node for obj (without data)
        gsXmlNode* bs_node = internal::makeNode("Basis" , data);
        bs_node->append_attribute( makeAttribute("type", "BSplineBasis", data) );

        // Write the knot vector
        gsXmlNode* tmp =
            internal::gsXml<gsKnotVector<T> >::put(obj.knots(), data );
        bs_node->append_node(tmp);

        // All set, return the BSPlineBasis node
        return bs_node;
    }

};


}// namespace internal

} // namespace gismo
