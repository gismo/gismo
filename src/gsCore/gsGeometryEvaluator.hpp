/** @file gsGeometryEvaluator.hpp

    @brief Provides implementation of GeometryEvaluation interface.

    This file is part of the G+Smo library-

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, J. Sogn, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

template<typename T,short_t ParDim>
void secDerToHessian(typename gsMatrix<T,ParDim*(ParDim+1)/2,1>::constRef & secDers,
                     gsMatrix<T,ParDim,ParDim> & hessian)
{
    switch ( ParDim )
    {
    case 1:
        hessian(0,0)=secDers(0,0);
        break;
    case 2:
        hessian(0,0)=secDers(0,0);
        hessian(1,1)=secDers(1,0);
        hessian(0,1)=
        hessian(1,0)=secDers(2,0);
        break;
    case 3:
        hessian(0,0)=secDers(0,0);
        hessian(1,1)=secDers(1,0);
        hessian(2,2)=secDers(2,0);
        hessian(0,1)=
        hessian(1,0)=secDers(3,0);
        hessian(0,2)=
        hessian(2,0)=secDers(4,0);
        hessian(1,2)=
        hessian(2,1)=secDers(5,0);
        break;
    default:
        break;
    }
}

template<typename T,short_t ParDim>
void hessianToSecDer (const gsMatrix<T,ParDim,ParDim> & hessian,
                      typename gsMatrix<T>::Row secDers)
{
    switch ( ParDim )
    {
    case 1:
        secDers(0,0)=hessian(0,0);
        break;
    case 2:
        secDers(0,0)=hessian(0,0);
        secDers(0,1)=hessian(1,1);
        secDers(0,2)=(hessian(1,0)+hessian(0,1)) / 2.0;
        break;
    case 3:
        secDers(0,0)=hessian(0,0);
        secDers(0,1)=hessian(1,1);
        secDers(0,2)=hessian(2,2);
        secDers(0,3)=(hessian(0,1)+hessian(1,0)) / 2.0;
        secDers(0,4)=(hessian(0,2)+hessian(2,0)) / 2.0;
        secDers(0,5)=(hessian(1,2)+hessian(2,1)) / 2.0;
        break;
    default:
        break;
    }
}


/*
// Partial derivatives of the Jacobian matrix
template<typename T,short_t ParDim> void
secDerToJacPartial(const typename gsMatrix<T>::constColumn  & secDer,
                   gsMatrix<T,ParDim,ParDim> * DJac)
{
    switch ( ParDim )
    {
    case 1:
        *DJac = secDer;
        break;
    case 2:
        // Note: stride is 3
        //dx
        DJac[0](0,0) = secDer(0);//G1xx
        DJac[0](0,1) = secDer(2);//G1xy
        DJac[0](1,0) = secDer(3);//G2xx
        DJac[0](1,1) = secDer(5);//G2xy
        //dy
        DJac[1](0,0) = secDer(2);//G1yx
        DJac[1](0,1) = secDer(1);//G1yy
        DJac[1](1,0) = secDer(5);//G2yx
        DJac[1](1,1) = secDer(4);//G2yy
        break;
    case 3:
        // Note: stride is 6
        //du
        DJac[0](0,0) = secDer(0 );//G1uu
        DJac[0](0,1) = secDer(3 );//G1uv
        DJac[0](0,2) = secDer(4 );//G1uw
        DJac[0](1,0) = secDer(6 );//G2uu
        DJac[0](1,1) = secDer(9 );//G2uv
        DJac[0](1,2) = secDer(10);//G2uw
        DJac[0](2,0) = secDer(12);//G3uu
        DJac[0](2,1) = secDer(15);//G3uv
        DJac[0](2,2) = secDer(16);//G3uw
        //dv
        DJac[1](0,0) = secDer(3 );//G1vu
        DJac[1](0,1) = secDer(1 );//G1vv
        DJac[1](0,2) = secDer(5 );//G1vw
        DJac[1](1,0) = secDer(9 );//G2vu
        DJac[1](1,1) = secDer(7 );//G2vv
        DJac[1](1,2) = secDer(11);//G2vw
        DJac[1](2,0) = secDer(15);//G3vu
        DJac[1](2,1) = secDer(13);//G3vv
        DJac[1](2,2) = secDer(17);//G3vw
        //dw
        DJac[2](0,0) = secDer(4 );//G1wu
        DJac[2](0,1) = secDer(5 );//G1wv
        DJac[2](0,2) = secDer(2 );//G1ww
        DJac[2](1,0) = secDer(10);//G2wu
        DJac[2](1,1) = secDer(11);//G2wv
        DJac[2](1,2) = secDer(8 );//G2ww
        DJac[2](2,0) = secDer(16);//G3wu
        DJac[2](2,1) = secDer(17);//G3wv
        DJac[2](2,2) = secDer(14);//G3ww
        break;
    default:
        break;
    }
}
//*/

template<typename T,short_t ParDim,short_t GeoDim>
void secDerToTensor(const typename gsMatrix<T>::constColumn & secDers,
                    gsMatrix<T,ParDim,ParDim> * a)
{
    static const int dim = ParDim*(ParDim+1)/2;
    for(int i=0;i<GeoDim;++i)
        secDerToHessian<T,ParDim>(secDers.template segment<dim>(i*dim),a[i]);
}

// Geometry transformation
template <class T, short_t ParDim, short_t GeoDim>
struct gsGeoTransform
{
    // Integral transformation for volume integrals
    static void getVolumeElements( gsMatrix<T> const & jacobians,
                                   gsVector<T> & result)
    {
        const index_t numPts = jacobians.cols() / ParDim;
        result.resize(numPts);

        for (index_t i = 0; i < numPts; ++i)
        {
            const Eigen::Block<const typename gsMatrix<T>::Base,GeoDim,ParDim> Ji =
                    jacobians.template block<GeoDim,ParDim>(0, i * ParDim);

            result[i] = math::sqrt( (Ji.transpose() * Ji).determinant() );
        }
    }

    // Transformation for parametric gradient
    static void getGradTransform( gsMatrix<T> const & jacobians,
                                  gsMatrix<T> & result)
    {
        const index_t numPts = jacobians.cols() / ParDim;
        result.resize(GeoDim, numPts * ParDim);

        for (index_t i = 0; i < numPts; ++i)
        {
            const Eigen::Block<const typename gsMatrix<T>::Base,GeoDim,ParDim> Ji =
                    jacobians.template block<GeoDim,ParDim>(0, i * ParDim);

            result.template block<GeoDim,ParDim>(0, i*ParDim) =
                Ji *  ( Ji.transpose() * Ji ).inverse().eval(); // temporary here
        }
    }

    // Laplacian transformation for parametric second derivatives
    static void getLaplTransform(index_t k,
                                 gsMatrix<T> const & jacobians,
                                 const typename gsMatrix<T>::Block allDerivs2,
                                 gsMatrix<T> & result)
    {
        //
    }

    // Outer (co-)normal vector at point k
    static void getOuterNormal(index_t k,
                               gsMatrix<T> const & jacobians,
                               boxSide s,
                               gsMatrix<T> & result)
    {
        // Assumes points k on boundary "s"
        // Assumes codim = 0 or 1

        const T   sgn = sideOrientation(s); // (!) * m_jacSign;
        const int dir = s.direction();

        result.resize(GeoDim);

        const gsMatrix<T, GeoDim, ParDim> Jk =
                jacobians.template block<GeoDim,ParDim>(0, k*ParDim);

        T alt_sgn = sgn;

        gsMatrix<T, ParDim-1, ParDim-1> minor;
        for (int i = 0; i != GeoDim; ++i) // for all components of the normal
        {
            Jk.firstMinor(i, dir, minor);
            result[i] = alt_sgn * minor.determinant();
            alt_sgn  *= -1;
        }

        /*
        // tangent vector
        const gsVector<T,GeoDim> tangent  = sgn * jacobians.template block<GeoDim, 1>(0,!dir);
        if ( tangent.squaredNorm() < 1e-10 )
        gsWarn<< "Zero tangent.\n";

        // Manifold outer normal vector
        const gsVector<T,3> outer_normal = jacobians.template block<GeoDim, 1>(0,0).cross(
        jacobians.template block<GeoDim, 1>(0,1) ).normalized();

        // Co-normal vector
        result = tangent.cross( outer_normal ) ;
        */
    }

    static void getNormal(index_t k,
                          gsMatrix<T> const & jacobians,
                          boxSide s,
                          gsMatrix<T> & result)
    {
        //
    }
};

// Geometry transformation : Specialization for co-dimension zero
template <class T, short_t ParDim>
struct gsGeoTransform<T,ParDim,ParDim>
{
    // Integral transformation for volume integrals
    static void getVolumeElements(gsMatrix<T> const & jacobians,
                                  gsVector<T> & result)
    {
        const index_t numPts = jacobians.cols() / ParDim;
        result.resize(numPts);

        for (index_t i = 0; i < numPts; ++i)
        {
            result[i] =
                    math::abs((jacobians.template block<ParDim,ParDim>(0, i*ParDim)).determinant());
        }
    }

    // Transformation for parametric gradient
    static void getGradTransform( gsMatrix<T> const & jacobians,
                                  gsMatrix<T> & result)
    {
        const index_t numPts = jacobians.cols() / ParDim;
        result.resize(ParDim, numPts * ParDim);

        for (index_t i = 0; i < numPts; ++i)
        {
            result.template block<ParDim,ParDim>(0, i*ParDim) =
                    jacobians.template block<ParDim,ParDim>(0, i*ParDim).inverse().transpose();
        }
    }

    // Laplacian transformation at allDerivs2.col(k)
    static void getLaplTransform(index_t k,
                                 gsMatrix<T> const & jacobians,
                                 const typename gsMatrix<T>::Block & allDerivs2,
                                 gsMatrix<T> & result)
    {
        //
    }

    // Outer (co-)normal vector at point k
    static void getOuterNormal(index_t k,
                               gsMatrix<T> const & jacobians,
                               boxSide s,
                               gsMatrix<T> & result)
    {
        // Assumes points k on boundary "s"
        // Assumes codim = 0 or 1

        const T   sgn = sideOrientation(s); // (!) * m_jacSign;
        const int dir = s.direction();

        result.resize(ParDim);

        const gsMatrix<T, ParDim, ParDim> & Jk =
                jacobians.template block<ParDim,ParDim>(0, k*ParDim);

        gsMatrix<T, ParDim-1, ParDim-1> minor;
        T alt_sgn = sgn;
        for (int i = 0; i != ParDim; ++i) // for all components of the normal
        {
            Jk.firstMinor(i, dir, minor);
            result[i] = alt_sgn * minor.determinant();
            alt_sgn  *= -1;
        }
    }

    static void getNormal(index_t k,
                          gsMatrix<T> const & jacobians,
                          boxSide s,
                          gsMatrix<T> & result)
    {
        //
    }

};


template <class T, short_t ParDim, short_t codim> void
gsGenericGeometryEvaluator<T,ParDim,codim>::
outerNormal(index_t k, boxSide s, gsVector<T> & result) const
{
    GISMO_ASSERT(this->m_flags & NEED_JACOBIAN, "Jacobians not computed");
    
    const T   sgn = sideOrientation(s) * m_orientation;
    const int dir = s.direction();
    
    // assumes points u on boundary "s"
    result.resize(GeoDim);
    if (ParDim + 1 == GeoDim) // surface case GeoDim == 3
    {
        const gsMatrix<T,GeoDim, ParDim> Jk = 
            m_jacobians.template block<GeoDim,ParDim>(0, k*ParDim);
        // fixme: generalize to nD
        normal(k,result);
        result = result.normalized().cross( sgn * Jk.template block<GeoDim, 1>(0,!dir) );
        
        /*
          gsDebugVar(result.transpose()); // result 1
          normal(k,result);
          Jk.col(dir) = result.normalized();
          gsMatrix<T, ParDim, ParDim> minor;
          T alt_sgn = sgn;
          for (int i = 0; i != GeoDim; ++i) // for all components of the normal
          {
          Jk.rowMinor(i, minor);
          result[i] = alt_sgn * minor.determinant();
          alt_sgn = -alt_sgn;
          }
          gsDebugVar(result.transpose()); // result 2
        //*/
    }
    else // planar case
    {
        GISMO_ASSERT( ParDim == GeoDim, "Codim different than zero/one");

        if ( 1==GeoDim ) { result[0] = sgn; return; } // 1D case

        const gsMatrix<T, ParDim, ParDim> Jk = 
            m_jacobians.template block<ParDim,ParDim>(0, k*ParDim);

        T alt_sgn = sgn;
        typename gsMatrix<T,ParDim,ParDim>::FirstMinorMatrixType minor;
        for (int i = 0; i != ParDim; ++i) // for all components of the normal
        {
            Jk.firstMinor(i, dir, minor);
            result[i] = alt_sgn * minor.determinant();
            alt_sgn = -alt_sgn;
        }
    }
}

template <class T, short_t ParDim, short_t codim>
void
gsGenericGeometryEvaluator<T,ParDim,codim>::evaluateAt(const gsMatrix<T>& u)
{
    GISMO_ASSERT( m_maxDeriv != -1, "Error in evaluation flags. -1 not supported yet.");

    m_numPts = u.cols();
    m_geo.basis().evalAllDers_into(u, m_maxDeriv, m_basisVals);

    // todo: If we assume all points (u) lie on the same element we may
    // store the actives only once
    m_geo.basis().active_into(u, m_active);

    if (this->m_flags & NEED_VALUE)
        computeValues();
    if (this->m_flags & NEED_JACOBIAN)
        computeJacobians();
    if (this->m_flags & NEED_MEASURE)
        gsGeoTransform<T,ParDim,GeoDim>::getVolumeElements(m_jacobians, m_measures);
    if (this->m_flags & NEED_GRAD_TRANSFORM)
        gsGeoTransform<T,ParDim,GeoDim>::getGradTransform(m_jacobians, m_jacInvs);
    if (this->m_flags & NEED_2ND_DER)
        compute2ndDerivs();
/*
    if (this->m_flags & NEED_DIV)
        divergence(m_div);
    if (this->m_flags & NEED_CURL)
        computeCurl();
    if (this->m_flags & NEED_LAPLACIAN)
        computeLaplacian();
    if (this->m_flags & NEED_NORMAL)
        computeNormal();
*/
}


template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeValues()
{
    const gsMatrix<T> & coefs = m_geo.coefs();
    const gsMatrix<T> & bVals = m_basisVals.front();

    m_values.resize(coefs.cols(), m_numPts);

    for (index_t j=0; j < m_numPts; ++j) // for all evaluation points
    {
        m_values.col(j) =  coefs.row( m_active(0,j) ) * bVals(0,j);

        for ( index_t i=1; i< m_active.rows() ; i++ )   // for all non-zero basis functions
            m_values.col(j)  +=   coefs.row( m_active(i,j) ) * bVals(i,j);
    }
}

template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeJacobians()
{
    const index_t numActive = m_active.rows();
    const gsMatrix<T> & coefs = m_geo.coefs();
    const gsMatrix<T> & bVals = m_basisVals[1];

    m_jacobians.setZero(GeoDim, m_numPts * ParDim);

    for (index_t j = 0; j < m_numPts; ++j)
        for (index_t i = 0; i < numActive; ++i) // for all active basis functions
        {
            m_jacobians.template block<GeoDim,ParDim>(0,j*ParDim).transpose().noalias() +=
                    bVals.template block<ParDim,1>(i*ParDim, j) *
                    coefs.template block<1,GeoDim>(m_active(i,j),0) ;

            // to check: which is faster ? the above or this..
            // m_jacobians.template block<GeoDim,ParDim>(0,j*ParDim).noalias() +=
            //     (bVals.template block<ParDim,1>(i*ParDim, j)
            //      * coefs.template block<1,GeoDim>(m_active(i,j),0) ).transpose();
        }
}


template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::compute2ndDerivs()
{
    const index_t numActive = m_active.rows();
    const gsMatrix<T> & coefs = m_geo.coefs();
    const gsMatrix<T> & bVals = m_basisVals[2];
    // Number of differnt 2ed derivativs combinations
    const index_t numDeriv = ParDim + (ParDim*(ParDim - 1))/2;

    m_2ndDers.setZero(GeoDim*numDeriv, m_numPts);

    gsMatrix<T> reshape; // to be removed later

    for (index_t j = 0; j < m_numPts; ++j)
    {
        for (index_t i = 0; i < numActive; ++i) // for all active basis functions
        {
            reshape.noalias() =
                    bVals.template block<numDeriv,1>(i*numDeriv, j)
                    * coefs.template block<1,GeoDim>(m_active(i,j),0);

            reshape.resize(GeoDim*numDeriv,1);
            m_2ndDers.template block<GeoDim*numDeriv,1>(0,j).noalias() += reshape;
        }
    }
}


template <class T, short_t ParDim, short_t codim> // AM: Not yet tested
void gsGenericGeometryEvaluator<T,ParDim,codim>::
transformValuesHdiv( index_t k,
                     const std::vector<gsMatrix<T> >& allValues,
                     gsMatrix<T>  & result) const
{
    GISMO_ASSERT( static_cast<size_t>(ParDim) == allValues.size(), "The number of components must equal the dimension of the domain.");
    int numA = 0;
    for(int comp = 0; comp < ParDim; ++comp)
        numA += allValues[comp].rows();
    result.resize(GeoDim,numA); // GeoDim x numA

    const T det = this->jacDet(k);
    const typename gsMatrix<T>::constColumns & jac = this->jacobian(k);
    index_t c = 0;
    for(int comp = 0; comp < ParDim; ++comp)  // component
    {
        const typename gsMatrix<T>::constColumn & bVals = allValues[comp].col(k);
        for( index_t j=0; j< allValues[comp].rows() ; ++j) // active of component
        {
            result.col(c++) = ( bVals[j] / det ) * jac.col(comp);
        }
    }
}

//template <class T, short_t ParDim, short_t codim>  // AM: Not yet tested
//void gsGenericGeometryEvaluator<T,ParDim,codim>::
//transformGradsHdiv( index_t k,
//                    const std::vector<gsMatrix<T> >& allValues,
//                    std::vector<gsMatrix<T> > const & allGrads,
//                    gsMatrix<T> & result) const
//{
// /*
//    //Assumptions: GeoDim = ParDim = TargetDim

//    index_t c = 0;
//    for(size_t comp = 0; comp < allValues.size(); ++comp)
//        c += allValues[c].rows();
//    result.resize(GeoDim*ParDim,c);

//    const T det = this->jacDet(k);
//    const typename gsMatrix<T>::constColumn  & secDer = this->deriv2(k);
//    const typename gsMatrix<T>::constColumns & Jac    = this->jacobian(k);
//    const Eigen::Transpose< const typename gsMatrix<T>::constColumns > &
//        invJac = this->gradTransform(k).transpose();

//    gsMatrix<T,GeoDim,ParDim> DJac[ParDim];
//    secDerToJacPartial<ParDim,T>(secDer,DJac);

//    gsVector<T> gradDetJrec(ParDim);
//    for (int i=0; i<ParDim; ++i)
//        gradDetJrec[i] = - ( invJac * DJac[i] ).trace() / det;

//    c = 0;
//    for(int comp = 0; comp < GeoDim; ++comp) // component
//    {
//        const typename gsMatrix<T>::constColumn & bvals = allValues[comp].col(k);
//        const typename gsMatrix<T>::constColumn & bder  = allGrads [comp].col(k);

//        for( index_t j=0; j< bvals.rows() ; ++j) // active of component
//        {
//            gsAsMatrix<T> tGrad(result.col(c).data(), GeoDim, GeoDim);

//            tGrad.noalias() = Jac.col(comp) * bder.segment(j*ParDim,ParDim).transpose() * invJac / det;

//            // tGrad.colwise() += gradDetJrec[i];
//            for (int i=0; i<ParDim; ++i) // result column
//            {
//                tGrad.col(i).array()   += gradDetJrec[i];  //JS2: I think this is incorrect!

//                tGrad.col(i) += ( allValues[comp](j,k) / det ) * DJac[i].col(j);
//                }

//            ++c;// next basis function
//        }
//    }
// //*/
//}

//Verified for 2D
template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::
transformLaplaceHgrad(  index_t k,
                        const gsMatrix<T> & allGrads,
                        const gsMatrix<T> & allHessians,
                        gsMatrix<T> & result) const
{
    gsMatrix<T> hessians;
    transformDeriv2Hgrad(k,allGrads,allHessians,hessians);
    result=hessians.leftCols(ParDim).rowwise().sum(); // trace of Hessian
    result.transposeInPlace();
}

//Verified for 2D
template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::
transformDeriv2Hgrad(  index_t k,
                        const gsMatrix<T> & funcGrad,
                        const gsMatrix<T> & funcSecDir,
                        gsMatrix<T> & result) const
{
    GISMO_ASSERT(
                   (ParDim==1 && (GeoDim==1||GeoDim==2||GeoDim==3))
                || (ParDim==2 && (GeoDim==2||GeoDim==3))
                || (ParDim==3 && GeoDim==3 ), "No implementation for this case");

    // important sizes
    static const int parSecDirSize = ParDim*(ParDim+1)/2;
    static const int fisSecDirSize = GeoDim*(GeoDim+1)/2;

    // allgrads
    const index_t numGrads = funcGrad.rows() / ParDim;

    result.resize(numGrads,fisSecDirSize);

    typename gsMatrix<T,GeoDim,ParDim>::constRef JM1 =
        m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim);
    typename gsMatrix<T,ParDim,GeoDim>::constRef JMT =
        m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim).transpose();

    // First part: J^-T H(u) J^-1
    gsMatrix<T,ParDim,ParDim> parFuncHessian;
    for (index_t i = 0; i < numGrads ; ++i)
    {
        secDerToHessian<T,ParDim>(
            funcSecDir.template block<parSecDirSize,1>(i*parSecDirSize,k), parFuncHessian);
        hessianToSecDer<T,GeoDim>( JM1 * parFuncHessian * JMT, result.row(i) );
    }

    // Second part: sum_i[  J^-T H(G_i) J^-1 ( e_i^T J^-T grad(u) ) ]
    const typename gsMatrix<T>::constColumn  & secDer = this->deriv2(k);
    gsMatrix<T,ParDim,ParDim> DDG[GeoDim]; // Each matrix is the Hessian of a component of the Geometry
    secDerToTensor<T,ParDim,GeoDim>(secDer, DDG);
    gsMatrix<T> HGT(GeoDim,fisSecDirSize);
    for(int i = 0; i < GeoDim; ++i)
        hessianToSecDer<T,GeoDim>(JM1 * DDG[i] * JMT, HGT.row(i));

    // Lastpart: substract part2 from part1
    const gsAsConstMatrix<T,ParDim> grads_k(funcGrad.col(k).data(), ParDim, numGrads);
    result.noalias() -=  grads_k.transpose() * JM1.transpose() * HGT;
    // 1 x d * d x d * d x d * d * s -> 1 x s
}


/*
template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeCurl()
{
    /// NOT TESTED FOR CORRECTNESS
    GISMO_ASSERT( ParDim==3 && ParDim == GeoDim, "Implemented only for 3D geometries in 3D spaces");
    const int numPoints=m_jacobians.cols()/ParDim;

    m_curl.resize(GeoDim,numPoints);
    for (int i=0; i<numPoints;++i)
    {
        const typename gsMatrix<T>::constColumns Jac = this->jacobian(i);

        m_curl(0,i)= Jac(2,3) - Jac(3,2);
        m_curl(1,i)= Jac(3,1) - Jac(1,3);
        m_curl(2,i)= Jac(1,2) - Jac(2,1);
    }
}

template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeLaplacian()
{
    /// NOT TESTED FOR CORRECTNESS
    const int numPoints=m_jacobians.cols()/ParDim;
    m_lap.resize(GeoDim,numPoints);
    for (int i=0; i<numPoints;++i)
    {
        m_lap.col(i)=m_2ndDers.template block<ParDim,GeoDim>(0,i*GeoDim).colwise().sum().transpose();
    }
}

template <class T, short_t ParDim, short_t codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeNormal()
{
    /// NOT TESTED FOR CORRECTNESS
    GISMO_ASSERT( codim==1, "Codimension should be equal to 1");
    const int numPoints=m_jacobians.cols()/ParDim;
    m_normal.resize(GeoDim,numPoints);
    for (int i=0; i<numPoints;++i)
    {
    // TODO remove the extra copy
        gsVector<T> temp;
        normal(i,temp);
        m_normal.col(i)=temp;
    }
}
*/

template<class T>
gsGeometryEvaluator<T> *
evaluator1(unsigned flags, const gsGeometry<T> & g)
{
    switch (g.coDim())
    {
        case 0:
            return new gsGenericGeometryEvaluator<T, 1, 0>(g, flags);
        case 1:
            return new gsGenericGeometryEvaluator<T, 1, 1>(g, flags);
        case 2:
            return new gsGenericGeometryEvaluator<T, 1, 2>(g, flags);
        default:
            GISMO_ERROR("Codimension problem.");
    }
}

template<class T>
gsGeometryEvaluator<T> *
evaluator2(unsigned flags, const gsGeometry<T> & g)
{
    switch (g.coDim())
    {
        case 0:
            return new gsGenericGeometryEvaluator<T, 2, 0>(g, flags);
        case 1:
            return new gsGenericGeometryEvaluator<T, 2, 1>(g, flags);
        case -1:
            return new gsGenericGeometryEvaluator<T, 2, -1>(g, flags);
        default:
            GISMO_ERROR("Codimension problem.");
    }
}

template<class T>
gsGeometryEvaluator<T> *
evaluator3(unsigned flags, const gsGeometry<T> & g)
{
    switch (g.coDim())
    {
        case 0:
            return new gsGenericGeometryEvaluator<T, 3, 0>(g, flags);
        case 1:
            return new gsGenericGeometryEvaluator<T, 3, 1>(g, flags);
        case -2:
            return new gsGenericGeometryEvaluator<T, 3, -2>(g, flags);
        default:
            GISMO_ERROR("Codimension problem.( codim=" << g.coDim());
    }
}

template<class T>
gsGeometryEvaluator<T> *
evaluator4(unsigned flags, const gsGeometry<T> & g)
{
    switch (g.coDim())
    {
        case 0:
            return new gsGenericGeometryEvaluator<T, 4, 0>(g, flags);
        case 1:
            return new gsGenericGeometryEvaluator<T, 4, 1>(g, flags);
        case -3:
            return new gsGenericGeometryEvaluator<T, 4, -3>(g, flags);
        default:
            GISMO_ERROR("Codimension problem, parDim=" << g.parDim()
                                                       << ", coDim=" << g.coDim() << ".");
    }
}

} // namespace gismo
