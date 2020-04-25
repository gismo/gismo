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
gsMatrix<T>
gsFunction<T>::jacobian(const gsMatrix<T>& u) const
{
    gsMatrix<T> result;
    this->jacobian_into( u, result );
    return result;
}

template <class T>
void gsFunction<T>::jacobian_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    // Compute component gradients as columns of result
    deriv_into(u, result);

    // Reshape the matrix to get one Jacobian block per evaluation point
    const short_t d = domainDim();     // dimension of domain
    result.resize(d, result.size()/d); //transposed Jacobians
    result.blockTransposeInPlace( targetDim() );
}

template <class T>
void gsFunction<T>::div_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{

    // Compute component gradients as columns of result
    deriv_into(u, result);
    //gsDebug << "deriv_into result:\n"   << result << "\n";


    const index_t numPts = u.cols();    // number of points to compute at
    const index_t d = domainDim();      // dimension of domain

    gsVector<T> tmp_div(numPts);        // tmp. divergence storage
    tmp_div.setZero();                  // initialize by zeros
    gsVector<T> resCol(d * d);

    for ( index_t p = 0; p < numPts; p++ ) { // for all evaluation points
        resCol = result.col(p);
        for ( index_t i = 0; i < d; i++ ) { tmp_div(p) += resCol(i * d + i); }
    }
    // Resizing the result to store the
    result.resize(1, numPts);
    result = tmp_div.transpose();

    //gsDebug << "div_into:\n"   << result << "\n";

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
    result.resize( parDim * tarDim,  numPts );

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
            //m_geo.eval_into(u, tmp);
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
gsMatrix<T> gsFunction<T>::laplacian( const gsMatrix<T>& u ) const
{
    //gsDebug << "Using finite differences (gsFunction::laplacian) for computing Laplacian.\n";
    int d = u.rows();
    gsVector<T> tmp( d );
    gsMatrix<T> res( d, u.cols());
    for ( int j = 0; j < d; j++ )
    {
        tmp.setZero();
        tmp( j, 0 )    = T( 0.0000000001 );
        res.row( j )  = 16 * ( this->eval(u.colwise() + tmp ) +
                                this->eval(u.colwise() - tmp ) ) -
                30 * ( this->eval( u ) );
        tmp( j, 0 )    = T( 0.0000000002 );
        res.row( j ) -=  ( this->eval( u.colwise() - tmp ) +
                            this->eval( u.colwise() + tmp ) ) ;
        res.row( j ) /= T( 0.0000000012 ) ;
    }
    return res;
}

template <class T>
int gsFunction<T>::newtonRaphson(const gsVector<T> & value,
                                 gsVector<T> & arg,
                                 bool withSupport,
                                 const T accuracy,
                                 int max_loop,
                                 double damping_factor) const
{
    const index_t n = targetDim();
    GISMO_ASSERT( value.rows() == n, "Invalid input values");
    const bool squareJac = (n == domainDim());

    gsMatrix<T> delta, jac, supp;
    if (withSupport)
        supp = support();

    int cnt = 0, iter = 0;

    do {
        //gsDebugVar(arg.transpose());
        if ( withSupport &&
             ((arg.array()<supp.col(0).array()).any() ||
              (arg.array()>supp.col(1).array()).any() ) )
        {
                if ( ++cnt > 3 ) 
                    return -1;
                else
                    arg = arg.cwiseMax( supp.col(0) ).cwiseMin( supp.col(1) );
        }
        else
            cnt = 0;

        // compute residual: value - f(arg)
        eval_into (arg, delta);
        delta = value - delta;

        if(delta.norm() <= accuracy)
            return iter;

        // compute Jacobian
        jacobian_into(arg, jac);

        // Solve for next update
        if (squareJac)
            delta = jac.partialPivLu().solve( delta );
        else// use pseudo-inverse
            delta = jac.colPivHouseholderQr().solve(
                        gsMatrix<T>::Identity(n,n)) * delta;

        // update arg with damping factor
        arg += damping_factor * delta;

    } while (++iter <= max_loop);

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
void gsFunction<T>::eval_component_into(const gsMatrix<T>&,
                                        const index_t,
                                        gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template <class T> gsMatrix<T>
gsFunction<T>::hess(const gsMatrix<T>& u, unsigned coord) const
{
    gsMatrix<T> hessian, secDers;
    this->deriv2_into(u, secDers);

    const index_t dim = this->domainDim();
    const index_t sz  = dim*(dim+1)/2;
    typename gsMatrix<T>::Rows ders = secDers.middleRows(coord*sz, sz);
    //const gsAsConstMatrix<T> ders(secDers.data(), sz, secDers.size() / sz );
    hessian.resize(dim*dim, ders.cols() );

    switch ( dim )
    {
    case 1:
        hessian = secDers; // ders
        break;
    case 2:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(1)=//1,0
        hessian.row(2)=ders.row(2);//0,1
        hessian.row(3)=ders.row(1);//1,1
        break;
    case 3:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(3)=//0,1
        hessian.row(1)=ders.row(3);//1,0
        hessian.row(6)=//0,2
        hessian.row(2)=ders.row(4);//2,0
        hessian.row(4)=ders.row(1);//1,1
        hessian.row(7)=//1,2
        hessian.row(5)=ders.row(5);//2,1
        hessian.row(8)=ders.row(2);//2,2
        break;
    default:
        break;
    }
    return hessian;
}

template <typename T, short_t domDim, short_t tarDim>
inline void computeAuxiliaryData (gsMapData<T> & InOut, int d, int n)
{
    //GISMO_ASSERT( domDim*tarDim == 1, "Both domDim and tarDim must have the same sign");
    const index_t numPts = InOut.points.cols();

    // Gradient transformation
    if (InOut.flags & NEED_GRAD_TRANSFORM)
    {
        // domDim<=tarDim makes sense

        InOut.fundForms.resize(domDim*tarDim, numPts);
        for (index_t p=0; p!=numPts; ++p)
        {
            const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);

            if ( tarDim == domDim && tarDim!=-1 )
                gsAsMatrix<T,tarDim,domDim>(InOut.fundForms.col(p).data(), n, d)
                        = jacT.inverse();
            else
                gsAsMatrix<T,tarDim,domDim>(InOut.fundForms.col(p).data(), n, d)
                        = jacT.transpose()*(jacT*jacT.transpose()).inverse();
        }
    }

    // Measure
    if (InOut.flags & NEED_MEASURE)
    {
        InOut.measures.resize(1,numPts);
        for (index_t p = 0; p < numPts; ++p) // for all points
        {
            typename gsAsConstMatrix<T,domDim,tarDim>::Tr jac =
                    gsAsConstMatrix<T,domDim,tarDim>(InOut.values[1].col(p).data(),d, n).transpose();
            if (tarDim == domDim && tarDim!=-1)
				InOut.measures(0,p) = math::abs(jac.determinant());
			else
				InOut.measures(0,p) = math::sqrt( ( jac.transpose()*jac  ).determinant() );


        }
    }

    // Normal vector of hypersurface
    if (tarDim!=-1 && tarDim==domDim+1 && InOut.flags & NEED_NORMAL)
    {
        GISMO_ASSERT( n == d + 1, "Codimension should be equal to one");

        typename gsMatrix<T,domDim,tarDim>::ColMinorMatrixType   minor;
        InOut.normals.resize(tarDim, numPts);

        for (index_t p = 0; p != numPts; ++p) // for all points
        {
            const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
            T alt_sgn(1);
            for (int i = 0; i != tarDim; ++i) //for all components of the normal
            {
                jacT.colMinor(i, minor);
                InOut.normals(i,p) = alt_sgn * minor.determinant();
                alt_sgn = -alt_sgn;
            }
        }
    }

    // Outer normal vector
    if ( InOut.flags & NEED_OUTER_NORMAL)
    {
        const T   sgn = sideOrientation(InOut.side);
        const int dir = InOut.side.direction();
        InOut.outNormals.resize(n,numPts);

        if (tarDim!=-1 && tarDim==domDim)
        {
            if ( 1==tarDim ) { InOut.outNormals.setConstant(sgn); return; } // 1D case

            typename gsMatrix<T,domDim,tarDim>::FirstMinorMatrixType minor;
            for (index_t p=0;  p!=numPts; ++p)
            {
                const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
                T alt_sgn = sgn * ( //jacT.rows()==jacT.cols() &&
                                    jacT.determinant()<0 ? -1 : 1);
                for (int i = 0; i != tarDim; ++i) //for all components of the normal
                {
                    jacT.firstMinor(dir, i, minor);
                    InOut.outNormals(i,p) = alt_sgn * minor.determinant();
                    alt_sgn  *= -1;
                }
            }
        }
        else
        {
            gsMatrix<T,domDim,domDim> metric(d,d);
            gsVector<T,domDim>      param(d);
            typename gsMatrix<T,domDim,domDim>::FirstMinorMatrixType minor;
            for (index_t p=0;  p!=numPts; ++p)
            {
                const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
                metric= (jacT*jacT.transpose());
                T alt_sgn = sgn;
                for (int i = 0; i != (domDim!=-1?domDim:d); ++i) //for all components of the normal
                {
                    metric.firstMinor(dir, i, minor);
                    param(i) = alt_sgn * minor.determinant();
                    alt_sgn  *= -1;
                }
                //note: metric.determinant() == InOut.measures.at(p)
                InOut.outNormals.col(p)=jacT.transpose()*param/metric.determinant();
            }
        }

    }

}


// Computes map data out of this map
template <class T>
void gsFunction<T>::computeMap(gsMapData<T> & InOut) const
{
    // Fill function data
    if (InOut.flags & NEED_GRAD_TRANSFORM || InOut.flags & NEED_MEASURE    ||
            InOut.flags & NEED_NORMAL         || InOut.flags & NEED_OUTER_NORMAL)
        InOut.flags = InOut.flags | NEED_GRAD;

    this->compute(InOut.points, InOut);

    // Fill extra data
    std::pair<short_t, short_t> Dim = this->dimensions();

    GISMO_ASSERT(Dim.first<10,             "Domain dimension is too big");
    GISMO_ASSERT(Dim.first<=Dim.second, "Singular map: target dimension is lower then the domain dimension");
    switch (10 * Dim.second + Dim.first)
    {
    // curves
    case 11: computeAuxiliaryData<T,1,1>(InOut, Dim.first, Dim.second); break;
    case 21: computeAuxiliaryData<T,1,2>(InOut, Dim.first, Dim.second); break;
//    case 31: computeAuxiliaryData<T,1,3>(InOut, Dim.first, Dim.second); break;
//    case 41: computeAuxiliaryData<T,1,4>(InOut, Dim.first, Dim.second); break;
    // surfaces
//  case 12: computeAuxiliaryData<T,2,1>(InOut, Dim.first, Dim.second); break;
    case 22: computeAuxiliaryData<T,2,2>(InOut, Dim.first, Dim.second); break;
    case 32: computeAuxiliaryData<T,2,3>(InOut, Dim.first, Dim.second); break;
//    case 42: computeAuxiliaryData<T,2,4>(InOut, Dim.first, Dim.second); break;
// volumes
//  case 13: computeAuxiliaryData<T,3,1>(InOut, Dim.first, Dim.second); break;
//  case 23: computeAuxiliaryData<T,3,2>(InOut, Dim.first, Dim.second); break;
    case 33: computeAuxiliaryData<T,3,3>(InOut, Dim.first, Dim.second); break;
//    case 43: computeAuxiliaryData<T,3,4>(InOut, Dim.first, Dim.second); break;
// 4D bulks
//  case 14: computeAuxiliaryData<T,4,1>(InOut, Dim.first, Dim.second); break;
//  case 24: computeAuxiliaryData<T,4,2>(InOut, Dim.first, Dim.second); break;
//  case 34: computeAuxiliaryData<T,4,3>(InOut, Dim.first, Dim.second); break;
    case 44: computeAuxiliaryData<T,4,4>(InOut, Dim.first, Dim.second); break;
    default: computeAuxiliaryData<T,-1,-1>(InOut, Dim.first, Dim.second); break;
    }

}


} // namespace gismo
