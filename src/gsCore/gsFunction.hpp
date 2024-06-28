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

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif
//#include <gsOptimizer/gsGradientDescent.h>
#include <gsOptimizer/gsFunctionAdaptor.h>

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

template<class T>
void gsFunction<T>::invertPoints(const gsMatrix<T> & points,
                                 gsMatrix<T> & result,
                                 const T accuracy, const bool useInitialPoint) const
{
    result.resize(this->domainDim(), points.cols() );
    gsVector<T> arg;
    for ( index_t i = 0; i!= points.cols(); ++i)
    {
        if (useInitialPoint)
            arg = result.col(i);
        else
            arg = this->parameterCenter();
        //arg = _argMinNormOnGrid(16);

        const int iter = this->newtonRaphson(points.col(i), arg, true, accuracy, 250);
        if (iter>100)
            gsWarn<< "Inversion took "<<iter<<" steps for "<< points.col(i).transpose() <<" (result="<< arg.transpose()<< ")\n";

        if (-1==iter)
        {
            gsWarn<< "Inversion failed for: "<< points.col(i).transpose() <<" (result="<< arg.transpose()<< ")\n";
            result.col(i).setConstant( std::numeric_limits<T>::infinity() );
        }
        else
            result.col(i) = arg;
    }
/* // alternative impl using closestPointTo
    result.resize(parDim(), points.cols() );
    gsVector<T> pt, arg;
    for ( index_t i = 0; i!= points.cols(); ++i )
    {
        pt = points.col(i);
        if (useInitialPoint)
            arg = result.col(i);

        this->closestPointTo(pt, arg, accuracy, useInitialPoint);
        if ( (this->eval(arg)-pt).norm()<=accuracy )
            result.col(i) = arg;
        else
        {
            //result.col(i) = arg;
            result.col(i).setConstant( std::numeric_limits<T>::infinity() );
        }
    }
*/
}

template<class T>
void gsFunction<T>::invertPointGrid(gsGridIterator<T,0> & git,
                                    gsMatrix<T> & result, const T accuracy,
                                    const bool useInitialPoint) const
{
    result.resize(this->domainDim(), git.numPoints() );
    gsVector<T> arg;
    auto cw = git.numPointsCwise();
    git.reset();
    index_t i = 0;
    int iter  = -1;
    gsInfo << "Iterations: ";
    for(;git; ++git)
    {
        if (-1==iter && !useInitialPoint)
            arg = this->parameterCenter();
        else if (-1==iter && useInitialPoint)
            arg = result.col(0);
        else
            arg = (i%cw[0]==0 ? result.col(i-cw[0]) : result.col(i-1) );

        iter = this->newtonRaphson(*git, arg, true, accuracy, 500);
        if (100 < iter) gsInfo << " "<< iter;
        if (-1==iter)
        {
            gsWarn<< "--- "<<i<<" Inversion failed for: "<< git->transpose() <<" ("
                  << (i%cw[0]==0 ? result.col(i-cw[0]) : result.col(i-1) ).transpose() << " -> " << arg.transpose()<< " )\n";
            //result.col(i).setConstant( std::numeric_limits<T>::infinity() );
        }
        //else gsWarn<< "+Inversion OK "<< git->transpose() <<" (result="<< arg.transpose()<< ")\n";

        result.col(i++) = arg;
    }
    gsInfo << "\n";
}

template <class T>
template <int mode,int _Dim>
int gsFunction<T>::newtonRaphson_impl(
    const gsVector<T> & value,
    gsVector<T> & arg,
    bool withSupport,
    const T accuracy,
    int max_loop,
    T damping_factor, T scale) const
{
    const index_t n = value.rows();
    const bool squareJac = (n == domainDim());//assumed for _Dim!=-1

    GISMO_ASSERT( arg.size() == domainDim(),
                  "Input argument has wrong dimensions: "<< arg.transpose() );

    gsMatrix<T> supp;
    gsMatrix<T,_Dim,(_Dim==-1?-1:1)> delta , residual;
    gsMatrix<T,_Dim,_Dim> jac;

    if (withSupport && support().cols()==2)
    {
        supp = support();
        GISMO_ASSERT( (arg.array()>=supp.col(0).array()).all() &&
                      (arg.array()<=supp.col(1).array()).all(),
                      "Initial point is outside the domain.\n point = "<<arg<<"\n domain = "<<supp);
    }
    int iter = 1;
    T rnorm[2]; rnorm[0]=1e100;
    //T alpha=.5, beta=.5;
    gsFuncData<T> fd(0==mode?(NEED_VALUE|NEED_DERIV):(NEED_DERIV|NEED_HESSIAN));

    do {
        this->compute(arg,fd);
        residual = (0==mode?fd.values[0]:fd.values[1]);

        residual.noalias() = value - scale*residual;
        rnorm[iter%2] = residual.norm();

        if(rnorm[iter%2] <= accuracy) // residual below threshold
        {
            //gsInfo <<"--- OK: Accuracy "<<rnorm[iter%2]<<" reached after "<<iter <<"iterations.\n";
            return iter;
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

        const T rr = ( 1==iter ? (T)1.51 : rnorm[(iter-1)%2]/rnorm[iter%2] ); //important to start with small step
        damping_factor = rr<1.5 ? math::max((T)0.1 + (rr/99),(rr-(T)0.5)*damping_factor) : math::min((T)1,rr*damping_factor);

        //Line search
        /*
        real_t val = (this->eval(arg + damping_factor*delta)*scale - value).norm();

        for(index_t i = 1; i!=11; ++i)
        {
            val = 0 ;
        }
        */

        //gsInfo << "Newton it " << iter << " arg=" << arg.transpose() << ", f(arg)="
        //       << (0==mode?fd.values[0]:fd.values[1]).transpose() << ", res=" << residual.transpose()
        //       <<" ("<<rr<<" ~ "<<damping_factor<<"), norm=" << rnorm[iter%2] << "\n";

        // update arg
        arg += damping_factor * delta;

        if ( withSupport  && support().cols()==2)
        {
            if ( delta.norm()<accuracy )
            {
                //gsInfo <<"OK: Newton reached boundary of support "<< delta.norm() <<"\n";
                return iter;
            }
            arg = arg.cwiseMax( supp.col(0) ).cwiseMin( supp.col(1) );
        }

    } while (++iter <= max_loop);

    gsWarn <<"--- Newton method did not converge after "<< max_loop
           <<" iterations. Residual norm: "<< rnorm[max_loop%2]<<".\n";

    return -1;
}

template <class T>
gsVector<T> gsFunction<T>::_argMinOnGrid(index_t numpts) const
{
    gsVector<T> result;
    gsMatrix<T> supp = this->support();
    if (0!=supp.size())
    {
        result = this->parameterCenter();
        gsGridIterator<T,CUBE> pt(supp, numpts);//per direction
        T val, mval = this->eval(result).value();
        for(;pt; ++pt)
        {
            val = this->eval(*pt).value();
            if (val < mval )
            {
                mval   = val;
                result = *pt;
            }
        }
    }
    else
    {
        //take random points?
        result.setZero( domainDim() );
    }
    return result;
}

template <class T>
gsVector<T> gsFunction<T>::_argMinNormOnGrid(index_t numpts) const
{
    gsVector<T> result;
    gsMatrix<T> supp = this->support();
    if (0!=supp.size())
    {
        result = this->parameterCenter();
        gsGridIterator<T,CUBE> pt(supp, numpts);//per direction
        T val, mval = this->eval(result).squaredNorm();
        for(;pt; ++pt)
        {
            if ( (val = this->eval(*pt).squaredNorm())<mval )
            {
                mval   = val;
                result = *pt;
            }
        }
    }
    else
    {
        //take random points?
        result.setZero( domainDim() );
    }
    return result;
}

template <class T>
int gsFunction<T>::newtonRaphson(const gsVector<T> & value,
                                 gsVector<T> & arg,
                                 bool withSupport,
                                 const T accuracy,
                                 int max_loop,
                                 T damping_factor) const
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
                                  T damping_factor) const
{
    GISMO_ASSERT(1==targetDim(), "Currently argMin works for scalar functions");
    gsVector<T> result;

    // Initial point
    if ( 0 != init.size() )
        result = give(init);
    else
    {
        result = _argMinOnGrid(20);
    }

    #ifdef gsHLBFGS_ENABLED
    gsFunctionAdaptor<T> fmin(*this);
    // gsIpOpt<T> solver( &fmin );
    //gsGradientDescent<T> solver( &fmin );
    gsHLBFGS<T> solver( &fmin );
    solver.options().setInt("MaxIterations",100);
    solver.options().setInt("Verbose",0);
    //MinStep..1e-12
    solver.solve(result);
    result = solver.currentDesign();
    //gsDebugVar(result);
    return result;
#else
    int dd=domainDim();
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
#endif

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

template<class T>
void gsFunction<T>::recoverPoints(gsMatrix<T> & xyz, gsMatrix<T> & uv, index_t k,
                                  const T accuracy) const
{
    gsVector<index_t> ind(xyz.rows()-1);
    for (index_t i = 0; i!= xyz.rows(); ++i)
        if (i<k) ind[i]=i;
        else if (i>k) ind[i-1]=i;       

    gsMatrix<T> pt = xyz(ind,gsEigen::all);
    gsFuncCoordinate<T> fc(*this, give(ind));

    //find low accuracy closest point
    // uv.resize(this->domainDim(), xyz.cols() );
    // for (index_t i = 0; i!= xyz.cols(); ++i)
    // {
    //     gsSquaredDistance2<T> dist2(fc, pt.col(i));
    //     uv.col(i) = dist2.argMin(accuracy, 10000) ;
    // }

    fc.invertPoints(pt,uv,accuracy,false); //true
    xyz = this->eval(uv);
    //possible check: pt close to xyz
}

template<class T>
void gsFunction<T>::recoverPointGrid(gsGridIterator<T,0> & git,
                                     gsMatrix<T> & xyz, gsMatrix<T> & uv,
                                     index_t k, const T accuracy) const
{
    gsVector<index_t> ind(xyz.rows()-1);
    for (index_t i = 0; i!= xyz.rows(); ++i)
        if (i<k) ind[i]=i;
        else if (i>k) ind[i-1]=i;

    gsMatrix<T> pt = xyz(ind,gsEigen::all);
    gsFuncCoordinate<T> fc(*this, give(ind));

    fc.invertPointGrid(git,uv,accuracy,false); //true
    xyz = this->eval(uv);
}

template <class T>
void gsFunction<T>::eval_component_into(const gsMatrix<T>&,
                                        const index_t,
                                        gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsFunction<T>::parameterCenter( const boxCorner& bc ) const
{
    gsMatrix<T> supp = this->support();
    const index_t dim = supp.rows();
    gsMatrix<T> coordinates(dim,1);
    gsVector<bool> boxPar = bc.parameters(dim);
    for (index_t d=0; d<dim;++d)
    {
        if (boxPar(d))
            coordinates(d,0) = supp(d,1);
        else
            coordinates(d,0) = supp(d,0);
    }
    return coordinates;
}

template<class T>
gsMatrix<T> gsFunction<T>::parameterCenter( const boxSide& bc ) const
{
    gsMatrix<T> supp = this->support();
    const index_t dim = supp.rows();
    gsMatrix<T> coordinates(dim,1);
    const index_t dir = bc.direction();
    for (index_t d=0; d<dim;++d)
    {
        if (d != dir)
            coordinates(d,0) = ( supp(d,1) + supp(d,0) ) / (T)(2);
        else if (bc.parameter())
            coordinates(d,0) = supp(d,1);
        else
            coordinates(d,0) = supp(d,0);
    }
    return coordinates;
}

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
inline void computeAuxiliaryData(const gsFunction<T> &src, gsMapData<T> & InOut, int d, int n)
{
    //GISMO_ASSERT( domDim*tarDim == 1, "Both domDim and tarDim must have the same sign");
    const index_t numPts = InOut.points.cols();
    
    // If the measure on a boundary is requested, calculate the outer normal vector as well
    if ( InOut.side!=boundary::none && (InOut.flags & NEED_MEASURE) ) InOut.flags|=NEED_OUTER_NORMAL;

    // Outer normal vector
    if ( InOut.flags & NEED_OUTER_NORMAL)
    {
        if (InOut.side==boundary::none)
            gsWarn<< "Computing boundary normal without a valid side.\n";
        // gsDebugVar( InOut.side);
        const T   sgn = sideOrientation(InOut.side);
        // gsDebugVar( sgn );
        const int dir = InOut.side.direction();
        InOut.outNormals.resize(n,numPts);

        if (tarDim!=-1 ? tarDim == domDim : n==d)
        {
            if ( 1==n ) { InOut.outNormals.setConstant(sgn); return; } // 1D case

            T det_sgn = 0;
            typename gsMatrix<T,domDim,tarDim>::FirstMinorMatrixType minor;
            // Determine Jacobian determinant's sign, assume constant along boundary
            if (InOut.flags & SAME_ELEMENT )
            {
                for (index_t p=0;  p!=numPts; ++p)
                {
                    const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
                    T detJacTcurr = jacT.determinant();
                    if ( math::abs(detJacTcurr) >= 1e-7 )
                    {
                        det_sgn = detJacTcurr < 0 ? -1 : 1;
                        break;
                    }
                }
                if ( 0 == det_sgn )
                {
                    gsMatrix<T> parameterCenter = src.parameterCenter(InOut.side);
                    T detJacTcurr = src.jacobian(parameterCenter).determinant();
                    det_sgn = detJacTcurr < 0 ? -1 : 1;
                }
            }
            for (index_t p=0;  p!=numPts; ++p)
            {
                const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
                // BUG: When the determinant is really close to zero but negative,
                // the result might be the opposite of what is expected because of alt_sgn
                
                if (! (InOut.flags &SAME_ELEMENT)  )
                {
                    T detJacTcurr = jacT.determinant();
                    det_sgn = math::abs(detJacTcurr) < 1e-7 ? 0 : 
                        ( detJacTcurr < 0 ? -1 : 1 );
                    if ( 0 == det_sgn )
                    {
                        gsMatrix<T> parameterCenter = src.parameterCenter(InOut.side);
                        detJacTcurr = src.jacobian(parameterCenter).determinant();
                        det_sgn = detJacTcurr < 0 ? -1 : 1;
                    }
                }
                GISMO_ENSURE(det_sgn!=0, "Cannot find a non-zero Jacobian determinant.\n" << InOut.points);

                T alt_sgn = sgn *  det_sgn;
                for (int i = 0; i != n; ++i) //for all components of the normal
                {
                    jacT.firstMinor(dir, i, minor);
                    InOut.outNormals(i,p) = alt_sgn * minor.determinant();
                    alt_sgn  *= -1;
                }
            }
        }
        else // lower-dim boundary case, d + 1 == n
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
                InOut.outNormals.col(p)=jacT.transpose()*param/metric.determinant();
                //InOut.outNormals.col(p)=(jacT.transpose()*param).normalized()*jacT.col(!dir).norm();
            }
        }

    }

    // Measure
    if (InOut.flags & NEED_MEASURE)
    {
        InOut.measures.resize(1,numPts);
        if (InOut.side==boundary::none) // If in the domain's interior
        {
            for (index_t p = 0; p < numPts; ++p) // for all points
            {
                typename gsAsConstMatrix<T,domDim,tarDim>::Tr jac =
                        gsAsConstMatrix<T,domDim,tarDim>(InOut.values[1].col(p).data(),d, n).transpose();
                if ( tarDim!=-1 ? tarDim == domDim : n==d )
                    InOut.measures(0,p) = math::abs(jac.determinant());
                else //unequal dimensions
                    InOut.measures(0,p) = math::sqrt( ( jac.transpose()*jac  )
                                                    .determinant() );
            }
        }
        else // If on boundary
        {
            GISMO_ASSERT(d==2, "Only works for boundary curves..");
            const int dir = InOut.side.direction();
            typename gsMatrix<T,domDim,tarDim>::ColMinorMatrixType   minor;
            InOut.measures.resize(1, numPts);
            for (index_t p = 0; p != numPts; ++p) // for all points
            {
                const gsAsConstMatrix<T,domDim,tarDim> jacT(InOut.values[1].col(p).data(), d, n);
                InOut.measures.at(p) = jacT.row(!dir).norm();
            }
            //InOut.measures = InOut.outNormals.colwise().norm(); // problematic on 3d curve boundary
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
    if ( (InOut.flags & NEED_NORMAL) && (tarDim!=-1 ? tarDim == domDim+1 : n==d+1) )
    {
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
        //gsDebugVar(InOut.normals);
    }

    // Second Fundamantan form of surface
    if (InOut.flags & NEED_2ND_FFORM)
    {
        // gsDebugVar(src);
        // gsDebugVar(InOut.normals);
        // gsDebugVar(n);
        // gsDebugVar(d+1);
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
    case 11: computeAuxiliaryData<T,1,1>(*this, InOut, Dim.first, Dim.second); break;
    case 21: computeAuxiliaryData<T,1,2>(*this, InOut, Dim.first, Dim.second); break;
//    case 31: computeAuxiliaryData<T,1,3>(*this, InOut, Dim.first, Dim.second); break;
//    case 41: computeAuxiliaryData<T,1,4>(*this, InOut, Dim.first, Dim.second); break;
    // surfaces
//  case 12: computeAuxiliaryData<T,2,1>(*this, InOut, Dim.first, Dim.second); break;
    case 22: computeAuxiliaryData<T,2,2>(*this, InOut, Dim.first, Dim.second); break;
    case 32: computeAuxiliaryData<T,2,3>(*this, InOut, Dim.first, Dim.second); break;
//    case 42: computeAuxiliaryData<T,2,4>(*this, InOut, Dim.first, Dim.second); break;
// volumes
//  case 13: computeAuxiliaryData<T,3,1>(*this, InOut, Dim.first, Dim.second); break;
//  case 23: computeAuxiliaryData<T,3,2>(*this, InOut, Dim.first, Dim.second); break;
    case 33: computeAuxiliaryData<T,3,3>(*this, InOut, Dim.first, Dim.second); break;
//    case 43: computeAuxiliaryData<T,3,4>(InOut, Dim.first, Dim.second); break;
// 4D bulks
//  case 14: computeAuxiliaryData<T,4,1>(*this, InOut, Dim.first, Dim.second); break;
//  case 24: computeAuxiliaryData<T,4,2>(*this, InOut, Dim.first, Dim.second); break;
//  case 34: computeAuxiliaryData<T,4,3>(*this, InOut, Dim.first, Dim.second); break;
    case 44: computeAuxiliaryData<T,4,4>(*this, InOut, Dim.first, Dim.second); break;
    default: computeAuxiliaryData<T,-1,-1>(*this, InOut, Dim.first, Dim.second); break;
    }

}


} // namespace gismo
