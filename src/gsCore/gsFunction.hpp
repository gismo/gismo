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
void gsFunction<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    gsDebug<< "Using finite differences (gsFunction::deriv_into) for derivatives.\n";
    const int d = u.rows();                     // dimension of domain
    const int n = eval(u.col(0))->rows();       // dimension of codomain
    const int numPts = u.cols();                // number of points to compute at

    gsVector<T> tmp(d);
    gsMatrix<T> ev, uc(d,4);
    result.resize( n, d * numPts );

    for ( int thisPt = 0; thisPt < numPts; thisPt++ )
    {
        for ( int j = 0; j<d; j++ )
        {
            int outputCol = thisPt * d + j;
            tmp.setZero();
            tmp(j) = T(0.00001);
            uc.col(0).noalias() = u.col(thisPt)+tmp;
            uc.col(1).noalias() = u.col(thisPt)-tmp;
            tmp(j) = T(0.00002);
            uc.col(2).noalias() = u.col(thisPt)+tmp;
            uc.col(3).noalias() = u.col(thisPt)-tmp;
            this->eval_into(uc, ev );

            result.col(outputCol) = (8*( ev.col(0)- ev.col(1)) + ev.col(3) - ev.col(2) ) / T(0.00012);
        }
    }
}


template <class T>
void gsFunction<T>::deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const
{
    gsDebug << "Using finite differences (gsFunction::deriv2_into) for second derivatives.\n";
    const int d = u.rows();                     // dimension of domain
    const int n = eval( u.col( 0 ) )->rows();   // dimension of codomain
    const int numPts = u.cols();                // number of points to compute at
    gsVector<T> tmp( d );
    gsMatrix<T> ev, uc( d, 3 ), ucm( d, 4 );
    const int stride = d + d * ( d - 1 ) / 2;
    result.resize( n * stride , numPts );
    for ( int thisPt = 0; thisPt < numPts; thisPt++ ) {
        int r = d;
        for ( int j = 0; j < d; j++ ) { // pure 2nd derivs
            tmp.setZero();
            tmp( j ) = T( 0.00001 );
            uc.col( 0 ).noalias() = u.col( thisPt ) + tmp;
            uc.col( 1 ).noalias() = u.col( thisPt )    ;
            uc.col( 2 ).noalias() = u.col( thisPt ) - tmp;
            this->eval_into( uc, ev );
            for ( int k = 0; k < n; k++ ) // for all coordinates
                result( k * stride + j, thisPt ) =
                    ( ev( k, 0 ) - 2 * ev( k, 1 ) + ev( k, 2 ) ) / T( 0.0000000001 ) ;
            for ( int l = j + 1; l < d; l++ ) { // pure 2nd derivs
                tmp( l ) = T( 0.00001 );
                ucm.col( 0 ).noalias() = u.col( thisPt ) + tmp;
                ucm.col( 3 ).noalias() = u.col( thisPt ) - tmp;
                tmp( l ) = -T( 0.00001 );
                ucm.col( 1 ).noalias() = u.col( thisPt ) + tmp;
                ucm.col( 2 ).noalias() = u.col( thisPt ) - tmp;
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
    //int n = eval(u.col(0))->rows();
    gsVector<T> tmp( d );
    gsMatrix<T>* res = new gsMatrix<T>( d, u.cols() );
    for ( int j = 0; j < d; j++ ) {
        tmp.setZero();
        tmp( j, 0 ) = T( 0.0000000001 );
        res->row( j ) = 16 * ( *this->eval( u.colwise() + tmp ) + *this->eval(
                    u.colwise() - tmp ) ) - 30 * ( *this->eval( u ) );
        tmp( j, 0 ) = T( 0.0000000002 );
        res->row( j ) += - ( *this->eval( u.colwise() - tmp ) + *this->eval( u.colwise() + tmp ) ) ;
        res->row( j ) /= T( 0.0000000012 ) ;
    }
    return res;
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
