/** @file gsFunction.hpp

    @brief Provides implementation of of Function common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFuncData.h>

#pragma once

namespace gismo
{

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
    // Compute component gradients as columns of result
    deriv_into(u,result);

    // Reshape the matrix to get one Jacobian block per evaluation point
    const index_t d = domainDim();     // dimension of domain
    result.resize(d, result.size()/d); //transposed Jacobians
    result.blockTransposeInPlace( targetDim() );
}

template <class T>
void gsFunction<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    //gsDebug<< "Using finite differences (gsFunction::deriv_into) for derivatives.\n";
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
    //gsDebug << "Using finite differences (gsFunction::deriv2_into) for second derivatives.\n";
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
    //gsDebug << "Using finite differences (gsFunction::laplacian) for computing Laplacian.\n";
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

/*
template <class T>
int gsFunction<T>::domainDim() const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
int gsFunction<T>::targetDim() const
{ GISMO_NO_IMPLEMENTATION }
*/

template <class T>
void gsFunction<T>::eval_component_into(const gsMatrix<T>& u, 
                                        const index_t comp,
                                        gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
typename gsFunction<T>::uMatrixPtr
gsFunction<T>::hess(const gsMatrix<T>& u, unsigned coord) const    
{ GISMO_NO_IMPLEMENTATION }


template <typename T, int domDim, int tarDim>
void computeAuxiliaryData (gsMapData<T> & InOut)
{
    const index_t  numPts = InOut.points.cols();

    if (InOut.flags & NEED_GRAD_TRANSFORM)
    {
        InOut.fundForms.resize(domDim*tarDim,numPts);
        for (index_t p=0; p<numPts; ++p)
        {
            gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(),domDim,tarDim);
            gsAsMatrix<T,tarDim,domDim>(InOut.fundForms.col(p).data(),tarDim,domDim) = jacT.transpose()*(jacT*jacT.transpose()).inverse().eval();
        }
    }

    if (InOut.flags & NEED_MEASURE)
    {
        InOut.measures.resize(1,numPts);
        for (index_t p = 0; p < numPts; ++p) // for all points
        {
            Eigen::Transpose< typename gsAsConstMatrix<T,domDim,tarDim>::Base > jac=gsAsConstMatrix<T,domDim,tarDim>(InOut.values[1].col(p).data(),domDim,tarDim).transpose();
            InOut.measures(0,p) = math::sqrt( ( jac.transpose()*jac  ).determinant() );
        }
    }

    if (InOut.flags & NEED_NORMAL)
    {
        GISMO_ASSERT( tarDim - domDim == 1, "Codimension should be equal to one");

        gsMatrix<T,domDim,tarDim-1> minor;
        InOut.normals.resize(tarDim, numPts);
        for (index_t p = 0; p != numPts; ++p) // for all points
        {
            const gsAsMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(),domDim,tarDim);
            real_t alt_sgn(1);
            for (int i = 0; i !=tarDim; ++i) // for all components of the normal vector
            {
                jacT.colMinor(i, minor);
                InOut.normals(i,p) = alt_sgn * minor.determinant();
                alt_sgn = -alt_sgn;
            }
        }
    }

    if (InOut.flags & NEED_OUTER_NORMAL )
    {
        const T   sgn = sideOrientation(InOut.side);
        const int dir = InOut.side.direction();
        InOut.outNormals.resize(tarDim,numPts);
        typename gsMatrix<T,domDim,tarDim>::FirstMinorMatrixType   minor;
        gsMatrix<T,domDim,tarDim>   jacT;
        for (index_t p=0; p<numPts ;++p)
        {
            jacT=gsAsMatrix<T,domDim,tarDim>(InOut.values[1].col(p).data(),domDim,tarDim);
            T alt_sgn = sgn * (jacT.rows()==jacT.cols() && jacT.determinant()<0 ? -1 : 1);
            for (int i = 0; i != tarDim; ++i) // for all components of the normal
            {
                jacT.firstMinor(dir,i, minor);
                InOut.outNormals(i,p) = alt_sgn * minor.determinant();
                alt_sgn  *= -1;
            }
        }
    }

}


// Computes map data out of this map
template <class T>
void gsFunction<T>::computeMap(gsMapData<T> & InOut) const
{
    // Fill function data
    if (InOut.flags & NEED_GRAD_TRANSFORM ||InOut.flags & NEED_MEASURE || InOut.flags & NEED_NORMAL || InOut.flags & NEED_OUTER_NORMAL)
        InOut.flags = InOut.flags | NEED_GRAD;

    this->compute(InOut.points, InOut);

    // Fill extra data
    const gsFuncInfo info = this->info();

    if (info.domainDim<=4 && info.targetDim <=4)
        switch ((info.domainDim-1)*4+info.targetDim-1)
        {
        // curves
        case  0: computeAuxiliaryData<T,1,1>(InOut); break;
        case  1: computeAuxiliaryData<T,1,2>(InOut); break;
        case  2: computeAuxiliaryData<T,1,3>(InOut); break;
        case  3: computeAuxiliaryData<T,1,4>(InOut); break;
            // surfaces
        case  4: computeAuxiliaryData<T,2,1>(InOut); break;
        case  5: computeAuxiliaryData<T,2,2>(InOut); break;
        case  6: computeAuxiliaryData<T,2,3>(InOut); break;
        case  7: computeAuxiliaryData<T,2,4>(InOut); break;
            // volumes
        case  8: computeAuxiliaryData<T,3,1>(InOut); break;
        case  9: computeAuxiliaryData<T,3,2>(InOut); break;
        case 10: computeAuxiliaryData<T,3,3>(InOut); break;
        case 11: computeAuxiliaryData<T,3,4>(InOut); break;
            // 4 dimensional volumes
        case 12: computeAuxiliaryData<T,4,1>(InOut); break;
        case 13: computeAuxiliaryData<T,4,2>(InOut); break;
        case 14: computeAuxiliaryData<T,4,3>(InOut); break;
        case 15: computeAuxiliaryData<T,4,4>(InOut); break;
        default:
            break;
        }
    else
    {

        const index_t  numPts = InOut.points.cols();



        if (InOut.flags & NEED_GRAD_TRANSFORM)
        {
            InOut.fundForms.resize(info.derivSize(),numPts);
            for (index_t p=0; p<numPts; ++p)
            {
                InOut.fundForms.reshapeCol(p, info.targetDim, info.domainDim) = InOut.jacobian(p)*(InOut.jacobian(p).transpose()*InOut.jacobian(p)).inverse().eval();
            }
        }

        if (InOut.flags & NEED_MEASURE)
        {
            InOut.measures.resize(1,numPts);
            for (index_t p = 0; p < numPts; ++p) // for all points
            {
                InOut.measures(0,p) = math::sqrt( ( InOut.jacobian(p).transpose()*InOut.jacobian(p)  ).determinant() );
            }
        }

        if (InOut.flags & NEED_NORMAL)
        {
            GISMO_ASSERT( info.targetDim - info.domainDim == 1, "Codimension should be equal to one");

            gsMatrix<T> minor;
            InOut.normals.resize(info.targetDim, numPts);
            for (index_t p = 0; p != numPts; ++p) // for all points
            {
                const gsAsMatrix<T> jac = InOut.values[1].reshapeCol(p, info.domainDim, info.targetDim);
                real_t alt_sgn(1);
                for (int i = 0; i != info.targetDim; ++i) // for all components of the normal vector
                {
                    jac.colMinor(i, minor);
                    InOut.normals(i,p) = alt_sgn * minor.determinant();
                    alt_sgn = -alt_sgn;
                }
            }
        }

        if (InOut.flags & NEED_OUTER_NORMAL )
        {
            const T   sgn = sideOrientation(InOut.side); // (!) * m_jacSign;
            const int dir = InOut.side.direction();
            InOut.outNormals.resize(info.targetDim,numPts);
            gsMatrix<T> jac;
            gsMatrix<T> minor;
            for (index_t p=0; p<numPts ;++p)
            {
                jac=InOut.jacobian(p);
                T alt_sgn = sgn * (jac.rows()==jac.cols() && jac.determinant()<0 ? -1 : 1);
                for (int i = 0; i != info.targetDim; ++i) // for all components of the normal
                {
                    jac.firstMinor(i, dir, minor);
                    InOut.outNormals(i,p) = alt_sgn * minor.determinant();
                    alt_sgn  *= -1;
                }
            }
        }
    }
}


} // namespace gismo
