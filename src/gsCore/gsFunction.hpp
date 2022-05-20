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
#include <gsCore/gsFuncCoordinate.h>
#include <gsTensor/gsGridIterator.h>

#pragma once


namespace gismo
{

template <class T>
gsFuncCoordinate<T> gsFunction<T>::coord(const index_t c) const
{
    return gsFuncCoordinate<T>(*this,c);
}

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
            delta(j)  = (T)(0.00001);
            uc.col(0) = u.col(p)+delta;
            uc.col(1) = u.col(p)-delta;
            delta(j)  = (T)(0.00002);
            uc.col(2) = u.col(p)+delta;
            uc.col(3) = u.col(p)-delta;
            //m_geo.eval_into(u, tmp);
            this->eval_into(uc, ev );
            tmp=(8*( ev.col(0)- ev.col(1)) + ev.col(3) - ev.col(2) ) / (T)(0.00012);

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
    const int stride = d + d * ( d - 1 ) / 2;   // number of second derivatives per component [= d(d+1)/2]
    result.resize( n * stride , numPts );
    for ( int thisPt = 0; thisPt < numPts; thisPt++ )
    {
        int r = d;
        for ( int j = 0; j < d; j++ )
        {
            // pure 2nd derivs
            tmp.setZero();
            tmp( j )    = (T)( 0.00001 );
            uc.col( 0 ) = u.col( thisPt ) + tmp;
            uc.col( 1 ) = u.col( thisPt )    ;
            uc.col( 2 ) = u.col( thisPt ) - tmp;
            this->eval_into( uc, ev );
            for ( int k = 0; k < n; k++ ) // for all coordinates
                result( k * stride + j, thisPt ) =
                        ( ev( k, 0 ) - (T)(2) * ev( k, 1 ) + ev( k, 2 ) ) / (T)( 0.0000000001 ) ;
            // mixed 2nd derivs
            for ( int l = j + 1; l < d; l++ )
            {
                tmp( l )     = (T)( 0.00001 );
                ucm.col( 0 ) = u.col( thisPt ) + tmp;
                ucm.col( 3 ) = u.col( thisPt ) - tmp;
                tmp( l )     = (T)( -0.00001 );
                ucm.col( 1 ) = u.col( thisPt ) + tmp;
                ucm.col( 2 ) = u.col( thisPt ) - tmp;
                tmp( l ) = (T)( 0 );
                this->eval_into( ucm, ev );
                for ( int k = 0; k < n; k++ ) // for all coordinates
                    result( k * stride + r, thisPt ) =
                            ( ev( k, 0 ) - ev( k, 1 ) - ev( k, 2 ) + ev( k, 3 ) ) / (T)( 0.0000000004 ) ;
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
        tmp( j, 0 )    = (T)( 0.0000000001 );
        res.row( j )  = 16 * ( this->eval(u.colwise() + tmp ) +
                                this->eval(u.colwise() - tmp ) ) -
                30 * ( this->eval( u ) );
        tmp( j, 0 )    = (T)( 0.0000000002 );
        res.row( j ) -=  ( this->eval( u.colwise() - tmp ) +
                            this->eval( u.colwise() + tmp ) ) ;
        res.row( j ) /= (T)( 0.0000000012 ) ;
    }
    return res;
}

template <class T>
template <int mode,int _Dim>
int gsFunction<T>::newtonRaphson_impl(
    const gsVector<T> & value,
    gsVector<T> & arg,
    bool withSupport,
    const T accuracy,
    int max_loop,
    double damping_factor, T scale) const
{
    const index_t n = value.rows();
    const T norm = value.norm() == 0 ? 1 : value.norm();
    const bool squareJac = (n == domainDim());//assumed for _Dim!=-1

    GISMO_ASSERT( arg.size() == domainDim(),
                  "Input argument has wrong dimensions: "<< arg.transpose() );

    gsMatrix<T> supp;
    gsMatrix<T,_Dim,(_Dim==-1?-1:1)> delta , residual;
    gsMatrix<T,_Dim,_Dim> jac;

    if (withSupport)
    {
        supp = support();
        GISMO_ASSERT( (arg.array()>=supp.col(0).array()).all() &&
                      (arg.array()<=supp.col(1).array()).all(),
                      "Initial point is outside the domain.");
    }
    int iter = 0;
    T rnorm[2]; rnorm[1]=1;
    //T alpha=.5, beta=.5;
    gsFuncData<> fd(0==mode?(NEED_VALUE|NEED_DERIV):(NEED_DERIV|NEED_HESSIAN));

    do {
        //gsInfo <<"Newton it: "<< arg.transpose()<<"\n";
        this->compute(arg,fd);
        residual = (0==mode?fd.values[0]:fd.values[1]);

        residual.noalias() = (value - scale*residual) / norm;// -->>>>>> NORMALIZE THIS??? / (residual.norm();
        rnorm[iter%2] = residual.norm();

        if(rnorm[iter%2] <= accuracy) // residual below threshold
        {
            // gsInfo <<"--- OK: Accuracy "<<rnorm[iter%2]<<" reached.\n";
            return iter;
        }

        if( iter>4 && (rnorm[(iter-1)%2]/rnorm[iter%2]) <1.1)
        {
            // gsInfo <<"--- OK: Converged to residual "<<rnorm[iter%2]<<" ("<<rnorm[(iter-1)%2]/rnorm[iter%2]<<"), niter = "<<iter<<".\n";
            return iter; //std::pair<iter,rnorm>
        }

        // compute Jacobian
        if (0==mode)
            jac = scale * fd.jacobian(0);
        else // (1==mode)
            jac = scale * fd.hessian(0);

        // Solve for next update
        if (squareJac)
        {
            //const T ddet = jac.determinant();
            //if (math::abs(ddet)<1e-4)
            //    gsWarn<< "Singular Jacobian: "<< ddet <<"\n";

            if (-1==_Dim)
                delta.noalias() = jac.partialPivLu().solve( residual );
            else
                delta.noalias() = jac.inverse() * residual;
        }
        else// use pseudo-inverse
            delta.noalias() = jac.colPivHouseholderQr().solve(
                        gsMatrix<T>::Identity(n,n)) * residual;

        // update arg
        arg += damping_factor * delta;

        if ( withSupport )
            arg = arg.cwiseMax( supp.col(0) ).cwiseMin( supp.col(1) );

    } while (++iter <= max_loop);

    gsWarn <<"--- Newton method did not converge after "<< max_loop
           <<" iterations. Residual norm: "<< rnorm[max_loop%2]<<".\n";
    return -1;
}

template <class T>
int gsFunction<T>::newtonRaphson(const gsVector<T> & value,
                                 gsVector<T> & arg,
                                 bool withSupport,
                                 const T accuracy,
                                 int max_loop,
                                 double damping_factor) const
{
    GISMO_ASSERT( value.rows() == targetDim(),
                  "Invalid input values:"<< value.rows()<<"!="<<targetDim());
    return newtonRaphson_impl<0>(value,arg,withSupport,accuracy,
                                 max_loop,damping_factor);
}

template <class T>
gsMatrix<T> gsFunction<T>::argMin(const T accuracy,
                                  int max_loop,
                                  gsMatrix<T> init,
                                  double damping_factor) const
{
    GISMO_ASSERT(1==targetDim(), "Currently argMin works for scalar functions");
    const index_t dd = domainDim();
    gsVector<T> result;

    // Initial point
    if ( 0 != init.size() )
        result = give(init);
    else
    {
        gsMatrix<T,2> supp = support();
        if (0!=supp.size())
        {
            gsGridIterator<T,CUBE> pt(supp, 5);//per direction
            T val, mval = std::numeric_limits<T>::max();
            for(;pt; ++pt)
            {
                if ( (val = this->eval(*pt).value())<mval )
                {
                    mval   = val;
                    result = *pt;
                }
            }
            //result = 0.5 * ( supp.col(0) + supp.col(1) );
        }
        else
            result.setZero( dd );
    }

    switch (dd)
    {
    case 2:
        newtonRaphson_impl<1,2>(gsVector<T>::Zero(dd), result, true,
                                accuracy,max_loop,damping_factor,(T)1);//argMax: (T)(-1)
        break;
    default:
        newtonRaphson_impl<1>(gsVector<T>::Zero(dd), result, true,
                              accuracy,max_loop,damping_factor,(T)1);//argMax: (T)(-1)
    }

    return result;
}

//argMax

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

template <class T> void
gsFunction<T>::hessian_into(const gsMatrix<T>& u, gsMatrix<T> & result,
                            index_t coord) const
{
    gsMatrix<T> secDers;
    this->deriv2_into(u, secDers);
    const index_t dim = this->domainDim();
    result = util::secDerToHessian(secDers, dim);
}

template <typename T, short_t domDim, short_t tarDim>
inline void computeAuxiliaryData(gsMapData<T> & InOut, int d, int n)
{
    //GISMO_ASSERT( domDim*tarDim == 1, "Both domDim and tarDim must have the same sign");
    const index_t numPts = InOut.points.cols();

    // Measure
    if (InOut.flags & NEED_MEASURE)
    {
        InOut.measures.resize(1,numPts);
        for (index_t p = 0; p < numPts; ++p) // for all points
        {
            typename gsAsConstMatrix<T,domDim,tarDim>::Tr jac =
                    gsAsConstMatrix<T,domDim,tarDim>(InOut.values[1].col(p).data(),d, n).transpose();
//            if (tarDim == domDim && tarDim!=-1)
            if ( tarDim!=-1 ? tarDim == domDim : n==d )
                InOut.measures(0,p) = math::abs(jac.determinant());
            else
                InOut.measures(0,p) = math::sqrt( ( jac.transpose()*jac  )
                                                  .determinant() );
        }
    }

    if (InOut.flags & NEED_GRAD_TRANSFORM)
    {
        // domDim<=tarDim makes sense

        InOut.jacInvTr.resize(d*n, numPts);
        for (index_t p=0; p!=numPts; ++p)
        {
            const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);

            if ( tarDim!=-1 ? tarDim == domDim : n==d )
                gsAsMatrix<T,tarDim,domDim>(InOut.jacInvTr.col(p).data(), n, d)
                        = jacT.cramerInverse();
            else
            {
                gsAsMatrix<T,tarDim,domDim>(InOut.jacInvTr.col(p).data(), n, d)
                        = jacT.transpose()*(jacT*jacT.transpose()).cramerInverse();
            }
        }
    }


    // Normal vector of hypersurface
    if (n==d+1 && InOut.flags & NEED_NORMAL)
    {
        GISMO_ASSERT( n == d + 1, "Codimension should be equal to one");

        typename gsMatrix<T,domDim,tarDim>::ColMinorMatrixType   minor;
        InOut.normals.resize(n, numPts);

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

    // Second Fundamantan form of surface
    if (InOut.flags & NEED_2ND_FFORM)
    {
        //domDim=2, tarDim=3
        InOut.fundForms.resize(d*d, numPts);
        const index_t sz = d*(d+1)/2;
        for (index_t p=0; p!=numPts; ++p)
        {
            const gsAsConstMatrix<T,-1,tarDim> ddT(InOut.values[2].col(p).data(), sz, n);
            const T nrm = InOut.normals.col(p).norm();
            if (0!=nrm)
            {
                InOut.fundForms(0,p) = ddT.row(0).dot(InOut.normals.col(p)) / nrm;
                InOut.fundForms(3,p) = ddT.row(1).dot(InOut.normals.col(p)) / nrm;
                InOut.fundForms(1,p) = InOut.fundForms(2,p) =
                    ddT.row(2).dot(InOut.normals.col(p)) / nrm;
            }
        }
    }

    // Outer normal vector
    if ( InOut.flags & NEED_OUTER_NORMAL)
    {
        if (InOut.side==boundary::none)
            gsWarn<< "Computing boundary normal without a valid side.\n";
        const T   sgn = sideOrientation(InOut.side);
        const int dir = InOut.side.direction();
        InOut.outNormals.resize(n,numPts);

        if (tarDim!=-1 ? tarDim == domDim : n==d)
        {
            if ( 1==n ) { InOut.outNormals.setConstant(sgn); return; } // 1D case

            typename gsMatrix<T,domDim,tarDim>::FirstMinorMatrixType minor;
            for (index_t p=0;  p!=numPts; ++p)
            {
                const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
                T alt_sgn = sgn * (T)( //jacT.rows()==jacT.cols() &&
                                    jacT.determinant()<0 ? -1 : 1);
                for (int i = 0; i != n; ++i) //for all components of the normal
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

    /*
    // Curvature of isoparametric curve
    if ( InOut.flags & NEED_CURVATURE)
    {
        //domDim=2, tarDim=3
        const int dir = InOut.side.direction();
        const index_t sz = domDim*(domDim+1)/2;
        InOut.curvature.resize(1,numPts);
        for (index_t p=0; p!=numPts; ++p)
        {
            const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
            const gsAsConstMatrix<T,-1,tarDim> ddT(InOut.values[2].col(p).data(), sz, n);
            const T nrm = ddT.row(dir).norm(); // (d^2/dir^2)G
            if (0!=nrm)
                InOut.curvature.at(p) = jacT.row(dir).cross(ddT.row(dir)).norm() / (nrm*nrm*nrm);
        }
    }
    */
}


// Computes map data out of this map
template <class T>
void gsFunction<T>::computeMap(gsMapData<T> & InOut) const
{
    // Fill function data
    if ( InOut.flags & (NEED_GRAD_TRANSFORM|NEED_MEASURE|NEED_NORMAL|NEED_OUTER_NORMAL) ) //NEED_JACOBIAN
        InOut.flags |= NEED_DERIV;
    if ( InOut.flags & (NEED_2ND_FFORM) ) //NEED_HESSIAN
        InOut.flags |= NEED_DERIV | NEED_DERIV2 | NEED_NORMAL;

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
