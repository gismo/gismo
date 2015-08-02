/** @file gsFunction.hpp

    @brief Provides implementation of of Function common operations.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsLinearAlgebra.h>

#pragma once

namespace gismo
{

template <class T>
typename gsFunction<T>::uMatrixPtr
gsFunction<T>::eval(const gsMatrix<T>& u) const
{
    gsMatrix<T>* result = new gsMatrix<T>;
    this->eval_into( u, *result );
    return uMatrixPtr(result);
}

template <class T>
typename gsFunction<T>::uMatrixPtr
gsFunction<T>::deriv(const gsMatrix<T>& u) const
{
    gsMatrix<T>* result = new gsMatrix<T>;
    this->deriv_into( u, *result );
    return uMatrixPtr(result);
}

template <class T>
typename gsFunction<T>::uMatrixPtr
gsFunction<T>::deriv2(const gsMatrix<T>& u) const
{
    gsMatrix<T>* result = new gsMatrix<T>;
    this->deriv2_into( u, *result );
    return uMatrixPtr(result);
}

template <class T>
typename gsFunction<T>::uMatrixPtr
gsFunction<T>::jacobian(const gsMatrix<T>& u) const
{
    gsMatrix<T>* result = new gsMatrix<T>;
    this->jacobian_into( u, *result );
    return uMatrixPtr(result);
}

template <class T>
void gsFunction<T>::jacobian_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    const index_t tarDim = targetDim(); // dimension of codomain

    // Compute component gradients as columns of result
    deriv_into(u,result);

    // Reshape the matrix to get one Jacobian block per evaluation point
    result.resize(tarDim, result.size()/tarDim);
    result.blockTransposeInPlace(domainDim());
}

template <class T>
void gsFunction<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    gsDebug<< "Using finite differences (gsFunction::deriv_into) for derivatives.\n";
    const index_t parDim = u.rows();                // dimension of domain
    const index_t tarDim = targetDim();             // dimension of codomain
    const index_t numPts = u.cols();                // number of points to compute at

    gsVector<T> delta(parDim);
    gsVector<T> tmp(tarDim);

    gsMatrix<T> ev, uc(parDim,4);
    result.resize( parDim *tarDim,  numPts );

    for ( index_t p = 0; p < numPts; p++ ) // for all evaluation points
    {
        for ( index_t j = 0; j<parDim; j++ ) // for all variables
        {
            delta.setZero();
            delta(j)  = T(0.00001);
            uc.col(0) = u.col(p)+delta;
            uc.col(1) = u.col(p)-delta;
            delta(j)  = T(0.00002);
            uc.col(2) = u.col(p)+delta;
            uc.col(3) = u.col(p)-delta;
            this->eval_into(uc, ev );
            tmp=(8*( ev.col(0)- ev.col(1)) + ev.col(3) - ev.col(2) ) / T(0.00012);

            for (index_t c=0; c<tarDim; ++c)  // for all components
                result(c*parDim+j,p)=tmp(c);
        }
    }
}


template <class T>
void gsFunction<T>::deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const
{
    gsDebug << "Using finite differences (gsFunction::deriv2_into) for second derivatives.\n";
    const int d = u.rows();                     // dimension of domain
    const int n = targetDim();                  // dimension of codomain
    const int numPts = u.cols();                // number of points to compute at
    gsVector<T> tmp( d );
    gsMatrix<T> ev, uc( d, 3 ), ucm( d, 4 );
    const int stride = d + d * ( d - 1 ) / 2;   // number of second derivatives per component
    result.resize( n * stride , numPts );
    for ( int thisPt = 0; thisPt < numPts; thisPt++ ) {
        int r = d;
        for ( int j = 0; j < d; j++ )
        {   
            // pure 2nd derivs
            tmp.setZero();
            tmp( j )    = T( 0.00001 );
            uc.col( 0 ) = u.col( thisPt ) + tmp;
            uc.col( 1 ) = u.col( thisPt )    ;
            uc.col( 2 ) = u.col( thisPt ) - tmp;
            this->eval_into( uc, ev );
            for ( int k = 0; k < n; k++ ) // for all coordinates
                result( k * stride + j, thisPt ) =
                    ( ev( k, 0 ) - 2 * ev( k, 1 ) + ev( k, 2 ) ) / T( 0.0000000001 ) ;
            // mixed 2nd derivs
            for ( int l = j + 1; l < d; l++ )
            {
                tmp( l )     = T( 0.00001 );
                ucm.col( 0 ) = u.col( thisPt ) + tmp;
                ucm.col( 3 ) = u.col( thisPt ) - tmp;
                tmp( l )     = -T( 0.00001 );
                ucm.col( 1 ) = u.col( thisPt ) + tmp;
                ucm.col( 2 ) = u.col( thisPt ) - tmp;
                tmp( l ) = T( 0 );
                this->eval_into( ucm, ev );
                for ( int k = 0; k < n; k++ ) // for all coordinates
                    result( k * stride + r, thisPt ) =
                        ( ev( k, 0 ) - ev( k, 1 ) - ev( k, 2 ) + ev( k, 3 ) ) / T( 0.0000000004 ) ;
                r++;
            }
        }
    }
}


template <class T>
gsMatrix<T>* gsFunction<T>::laplacian( const gsMatrix<T>& u ) const
{
    gsDebug << "Using finite differences (gsFunction::laplacian) for computing Laplacian.\n";
    int d = u.rows();
    gsVector<T> tmp( d );
    gsMatrix<T>* res = new gsMatrix<T>( d, u.cols() );
    for ( int j = 0; j < d; j++ ) 
    {
        tmp.setZero();
        tmp( j, 0 )    = T( 0.0000000001 );
        res->row( j )  = 16 * ( *this->eval(u.colwise() + tmp ) + 
                                *this->eval(u.colwise() - tmp ) ) - 
                         30 * ( *this->eval( u ) );
        tmp( j, 0 )    = T( 0.0000000002 );
        res->row( j ) -=  ( *this->eval( u.colwise() - tmp ) + 
                            *this->eval( u.colwise() + tmp ) ) ;
        res->row( j ) /= T( 0.0000000012 ) ;
    }
    return res;
}

template <class T>
int gsFunction<T>::newtonRaphson(const gsVector<T> & value,
                                  gsVector<T> & arg,
                                  bool withSupport, 
                                  const T accuracy,
                                  int max_loop) const
{
    const index_t n = targetDim();
    GISMO_ASSERT( value.rows() == n, "Invalid input values");
    const bool squareJac = (n == domainDim());

    gsMatrix<T> delta, jac, supp;
    if (withSupport)
        supp = support();

    int iter = 0;
    do {
        // compute residual: value - f(arg)
        eval_into (arg, delta);
        delta = value - delta;

        // compute Jacobian 
        jacobian_into(arg, jac);

        // Solve for next update
        if (squareJac)
            delta = jac.partialPivLu().solve( delta );
        else// use pseudo-inverse
            delta = jac.colPivHouseholderQr().solve(
                gsMatrix<T>::Identity(n,n)) * delta;

        // update arg
        arg += delta;

        // clamp x to the support of the function
        if (withSupport)
            arg = arg.cwiseMax( supp.col(0) ).cwiseMin( supp.col(1) );

        if (delta.norm() <= accuracy)
            return iter;
    } while (++iter < max_loop);

    // no solution found within max_loop iterations
    return -1;
}

template <class T>
gsFunction<T> * gsFunction<T>::clone() const 
{ GISMO_NO_IMPLEMENTATION }

template <class T>
int gsFunction<T>::domainDim() const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
int gsFunction<T>::targetDim() const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
gsMatrix<T> gsFunction<T>::support() const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
void gsFunction<T>::eval_component_into(const gsMatrix<T>& u, 
                                     const index_t comp, 
                                     gsMatrix<T>& result) const 
{ GISMO_NO_IMPLEMENTATION }

template <class T>
typename gsFunction<T>::uMatrixPtr
gsFunction<T>::hess(const gsMatrix<T>& u, unsigned coord) const    
{ GISMO_NO_IMPLEMENTATION }


} // namespace gismo
