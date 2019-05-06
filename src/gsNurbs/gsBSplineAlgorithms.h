/** @file gsBSplineAlgorithms.h

    @brief Implementation of common algorithms for B-splines

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/** @namespace gismo::bspline

    @brief
    This namespace contains implementation of B-spline related algorithms.
    
    \ingroup Nurbs
*/
namespace bspline
{
    /// Input: an iterator pointing to the the biggest knot less than \a u,
    /// evaluation point \a u, output the value of all basis functions
    /// which are active at \a u.
    /// Uses B-spline recursion for evaluation
    template <class T, typename KnotIterator, typename Derived>
    void evalBasis( T u,
                    KnotIterator knot,
                    int deg,
                    Eigen::MatrixBase<Derived> const & result )
    {       
        STACK_ARRAY(T, left, deg  + 1);
        STACK_ARRAY(T, right, deg + 1);
        evalBasis( u, knot, deg, left, right, result);
    }

    template <class T, typename KnotIterator, typename Derived>
    void evalBasis( T u,
                    KnotIterator knot,
                    int deg,
                    T left[], T right[], 
                    Eigen::MatrixBase<Derived> const & result )
    {
        Eigen::MatrixBase<Derived>& res = const_cast<Eigen::MatrixBase<Derived>&>(result);
        
        res(0,0)= T(1.0);  // 0-th degree function value
        
        for(int j=1; j<= deg; j++) // For all degrees ( ndu column)
        {
            left[j]  = u - *(knot+1-j);
            right[j] = *(knot+j) - u;
            T saved = T(0) ;
            for(int r=0; r<j ; r++) // For all (except the last)  basis functions of degree j ( ndu row)
            {
                const T temp = res(r,0) / ( right[r+1] + left[j-r] );
                res(r,0)     = saved + right[r+1] * temp ;// r-th function value of degree j
                saved = left[j-r] * temp ;
            }
            res(j,0)     = saved;
        }
    }

    /// Evaluation for degree 1 B-spline basis
    template <class T, typename KnotIterator, typename Derived>
    void evalDeg1Basis(  const T & u,
                         const KnotIterator knot,
                         Eigen::MatrixBase<Derived> const & result )
    {
        result(0,0) =  u       - (*knot+1);
        result(1,0) =  (*knot) - u        ;
        const T dn  = (*knot)  - (*knot+1) ;
        if( dn !=  0.0 )
            result.array /= dn;
    }
    
    /// Evaluation for degree 2 B-spline basis
    template <class T, typename KnotIterator, typename Derived>
    void evalDeg2Basis(  const T & u,
                         const KnotIterator knot,
                         Eigen::MatrixBase<Derived> const & result )
    {
        
    }
    
    /// Evaluation for degree 3 B-spline basis
    template <class T, typename KnotIterator, typename MatrixType>
    void evalDeg3Basis(  const T & u,
                         const KnotIterator knot,
                         MatrixType & result )
    {
        
    }
    

    /// Input: parameter position \a u, KnotIterator \a knot identifying the active interval,
    /// degree \a deg, Output: table \a N.
    /// Computes deBoor table used in Algorithm A2.5 in NURBS book.
    /// Used by evalDerSingle_into(...), evalAllDerSingle_into(...) and evalSingle_into(...).
    template <class T, typename KnotIterator>
    void deBoorTriangle( T u,
                         KnotIterator knot,
                         int deg,
                         T N[] )
    {
        const int table_size = deg+1;
        T saved;

        for( int j = 0; j <= deg; j++ ) // Initialize zeroth-degree functions
            if( u >= *(knot+j) && u < *(knot+j+1) )
                N[ j*table_size ] = 1;
            else
                N[ j*table_size ] = 0;
        for( int k = 1; k <= deg; k++ ) // Compute full triangular table
        {
            if( N[(k-1)] == 0 )
                saved = 0;
            else
                saved = ( (u-*(knot) )*N[ k-1 ])/( *(knot+k)- *(knot));
            for( int j = 0; j < deg-k+1; j++)
            {
                if( N[ (j+1)*table_size + (k-1) ] == 0 )
                {
                    N[ j*table_size + k ] = saved;
                    saved = 0;
                }
                else
                {
                    const T Uleft = *(knot+j+1);
                    const T Uright = *(knot+j+k+1);
                    const T temp = N[ (j+1)*table_size + (k-1) ]/(Uright-Uleft);
                    N[ j*table_size + k ] = saved + (Uright-u)*temp;
                    saved = (u-Uleft)*temp;
                }
            }
        }
    }


    /// Input: a sequence of 2*p+1 knots starting at iterator position
    /// \a knot, evaluation point \a u, output the value of the
    /// B-spline defined by \a coefs
    /// Uses DeBoor's algorithm for evaluation
    /// \todo add an overload with stride parameter for \a coefs
    template <class T, typename KnotIterator, typename MatrixType>
    void evalGeo(  const T & u,
                   const KnotIterator & knot,
                   int deg,
                   MatrixType coefs, // makes a local copy of coefs
                   MatrixType & result )
    { 
        result.resize(deg+1, 1 ) ;
        T tmp;    
        
        for ( int r=0; r!= deg; r++ ) // triangle step
            //for ( int r=s; r< deg; r++ ) // TO DO: improve using multiplicity s
            for ( int i=deg; i>r; i-- ) // recursive computation
            {
                tmp= ( u -  *(knot+i) ) / ( *(knot+i+deg-r) - *(knot+i) );
                coefs.row(i) = (T(1)-tmp)*coefs.row(i-1) + tmp* coefs.row(i) ;
            }
        result.noalias() = coefs.row(deg);
    }
    
    /// Input: a sequence of 2*p+1 knots starting at iterator position \a knot,
    /// evaluation point \a u, output the value of all basis functions
    /// which are active on \a u.
    /// Uses DeBoor's algorithm for evaluation
    template <class T, typename KnotIterator, typename MatrixType>
    void evalBasis2(  const T & u,
                      const KnotIterator & knot,
                      int deg,
                      MatrixType & result )
    {
        MatrixType pseudo_coefs;
        pseudo_coefs.setIdentity(deg+1,deg+1);
        evalGeoAlg(u, knot, deg, pseudo_coefs, result);
    }
    
    
    /// Input: a sequence of p+2 knots, evaluation point. Output: the
    /// value of the basis function supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void evalBasisSingle( ) 
    { 


    }

    /// Input: a sequence of 2*p+1 knots, evaluation point, output
    /// the derivatives evaluated at all functions supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void derivBasis( ) 
    { }

    /// Input: a sequence of p+2 knots, evaluation point. Output: the
    /// value of the derivative of basis function supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void derivBasisSingle( ) 
    { }

    /// Input: a sequence of 2*p+1 knots, evaluation point, output
    /// the value of all derivatives up to order k supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void allDersBasis( ) 
    { }

    /// Input: a sequence of p+2 knots, evaluation point, output
    /// the value of all derivatives up to order k supported on those knots
    template <class T, typename KnotIterator, typename MatrixType>
    void allDersSingle( ) 
    { }



/// Increase the degree of a 1D B-spline from degree p to degree p + m.
template<class Basis_t>
void degreeElevateBSpline(Basis_t &basis, 
                          gsMatrix<typename Basis_t::Scalar_t> & coefs,
                          short_t m)
{
    typedef typename Basis_t::Scalar_t T;

    GISMO_ASSERT(m >= 0 && m<512, "Can only elevate degree by a positive (not enormous) amount.");
    GISMO_ASSERT(basis.size() == coefs.rows(), "Invalid coefficients");

    if (m==0) return;

    const short_t p      = basis.degree();
    const index_t ncoefs = coefs.rows();
    const index_t n      = coefs.cols();
    const gsKnotVector<T> & knots = basis.knots();

    // compute original derivative coefficients P (recurrence)
    // original derivative coefficients
// #   if defined(__GNUC__)
//     STACK_ARRAY(gsMatrix<T>,P,p+1);
// #   else // Note: No non-POD VLAs in clang/MSVC compiler
    std::vector<gsMatrix<T> > P(p+1);
// #   endif
    for(short_t i=0;i<p+1;i++)
        P[i].setZero(ncoefs - i, n);

    // insert first row
    P[0].swap(coefs);
    // fill table of derivative coefficients
    for(short_t j=1; j<=p;j++)
        for(index_t i=0; i<ncoefs-j; i++)
        {
            if(knots[i+p+1]>knots[i+j])
                P[j].row(i).noalias() = 
                    ( P[j-1].row(i+1) - P[j-1].row(i) ) / ( knots[i+p+1] - knots[i+j] );
        }

    // loop over unique knot values (exept last one)
    const typename gsKnotVector<T>::multContainer &
        mult = knots.multiplicities(); // vector of mulitplicities

    // degree elevate basis
    basis.degreeElevate(m);
    const index_t ncoefs_new = basis.size();
    const short_t p_new      = basis.degree();

    // new (elevated) derivative coefficients
// #   if defined(__GNUC__)
//     STACK_ARRAY(gsMatrix<T>,Q,p_new+1);
// #   else // Note: No non-POD VLAs in clang/MSVC compiler
    std::vector<gsMatrix<T> > Q(p_new+1);
// #   endif
    for(short_t i=0; i<p_new+1; i++)
        Q[i].setZero(ncoefs_new - i, n);

    // loop over knot intervals (with positive measure):
    // prescibe known coefficients (see Huang Theorem 2) and fill
    // corresponding triangular table using backwards recurrence

    // precompute factors
    gsVector<T> factor = gsVector<T>::Ones(p+1);
    for(short_t j=0; j<=p; j++)
        for(short_t l=1; l<=j; l++)
            factor[j] *= static_cast<T>((p+1-l)) / (p_new+1-l);

    //set known coefficients
    // for(int j=0; j<=p; ++j)
    //     Q[j].row(0)=factor[j]*P[j].row(0);

    int betak = 0; // sum of interior mulitplicities
    for(size_t k=0; k<mult.size()-1; k++)
    {
        //set known coefficients
        for(short_t j=p+1-mult[k]; j<=p; ++j)
            Q[j].row(betak+k*m) = factor[j] * P[j].row(betak);

        for(short_t j=p_new-1; j>=0;j--)
        {
            // fill triangular table
            for(short_t i=1; i<=p_new-j; i++)
            {
                const int ik= i+betak+k*m; // update index i for the considered knot value
                if(knots[ik+p_new]>knots[ik+j])
                    Q[j].row(ik).noalias() = 
                        Q[j].row(ik-1) + Q[j+1].row(ik-1) * (knots[ik-1+p_new+1]-knots[ik-1+j+1]);

            }
        }
        
        betak += mult[k+1];
    }

    // Insert new coefficients into coefs
    coefs.swap( Q[0] );
}


} // namespace bspline

}// namespace gismo
