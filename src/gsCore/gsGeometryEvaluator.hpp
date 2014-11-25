/** @file gsGeometryEvaluator.hpp

    @brief Provides implementation of GeometryEvaluation interface.
    
    This file is part of the G+Smo library-
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, J. Sogn, A. Mantzaflaris
*/

namespace gismo
{

/// Geometry transformation 
template <class T, int ParDim, int GeoDim>
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
                    Ji *  ( Ji.transpose() * Ji ).inverse();
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
                               boundary::side s,
                               gsMatrix<T> & result)
    {
        // Assumes points k on boundary "s"
        // Assumes codim = 0 or 1

        const T   sgn = sideOrientation(s); // (!) * m_jacSign;
        const int dir = direction(s);

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
                          boundary::side s,
                          gsMatrix<T> & result)
    {
        //
    }
};

/// Geometry transformation : Specialization for co-dimension zero
template <class T, int ParDim>
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
                    math::fabs((jacobians.template block<ParDim,ParDim>(0, i*ParDim)).determinant());
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
                               boundary::side s,
                               gsMatrix<T> & result)
    {
        // Assumes points k on boundary "s"
        // Assumes codim = 0 or 1

        const T   sgn = sideOrientation(s); // (!) * m_jacSign;
        const int dir = direction(s);

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
                          boundary::side s,
                          gsMatrix<T> & result)
    {
        //
    }

};


template <class T, int ParDim, int codim>
void
gsGenericGeometryEvaluator<T,ParDim,codim>::evaluateAt(const gsMatrix<T>& u)
{
    GISMO_ASSERT( m_maxDeriv != -1, "Error in evaluation flags. -1 not supported yet.");

    m_numPts = u.cols();
    m_geo.basis().evalAllDers_into(u, m_maxDeriv, m_basisVals);
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
    if (this->m_flags & NEED_DIV)
        divergence(m_div);
    if (this->m_flags & NEED_CURL)
        computeCurl();
    if (this->m_flags & NEED_LAPLACIAN)
        computeLaplacian();
    if (this->m_flags & NEED_NORMAL)
        computeNormal();
}


template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeValues()
{
    const gsMatrix<T> & coefs = m_geo.coefs();

    m_values.resize(coefs.cols(), m_numPts);

    for (index_t j=0; j < m_numPts; ++j)                  // for all evaluation points
    {
        m_values.col(j) =  coefs.row( m_active(0,j) ) * m_basisVals(0,j);
        for ( index_t i=1; i< m_active.rows() ; i++ )   // for all non-zero basis functions
            m_values.col(j)  +=   coefs.row( m_active(i,j) ) * m_basisVals(i,j);
    }
}

template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeJacobians()
{
    const index_t numActive = m_active.rows();
    const gsMatrix<T> & coefs = m_geo.coefs();

    m_jacobians.setZero(GeoDim, m_numPts * ParDim);

    for (index_t j = 0; j < m_numPts; ++j)
        for (index_t i = 0; i < numActive; ++i) // for all active basis functions
        {
            m_jacobians.template block<GeoDim,ParDim>(0,j*ParDim).transpose().noalias() +=
                    m_basisVals.template block<ParDim,1>(numActive+i*ParDim, j) *
                    coefs.template block<1,GeoDim>(m_active(i,j),0) ;

            // to check: which is faster ? the above or this..
            // m_jacobians.template block<GeoDim,ParDim>(0,j*ParDim).noalias() +=
            //     (m_basisVals.template block<ParDim,1>(numActive+i*ParDim, j)
            //      * coefs.template block<1,GeoDim>(m_active(i,j),0) ).transpose();
        }
}




template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::compute2ndDerivs()
{
    const index_t numActive = m_active.rows();
    const gsMatrix<T> & coefs = m_geo.coefs();
    // Number of differnt 2ed derivativs combinations
    const index_t numDeriv = ParDim + (ParDim*(ParDim - 1))/2;
    // Starting index of the second derivatives in m_basisVals
    const index_t der2start = (1 + ParDim) * numActive;

    m_2ndDers.setZero(GeoDim*numDeriv, m_numPts);

    gsMatrix<T> reshape; // to be removed later

    for (index_t j = 0; j < m_numPts; ++j)
    {
        for (index_t i = 0; i < numActive; ++i) // for all active basis functions
        {
            reshape.noalias() =
                    m_basisVals.template block<numDeriv,1>(der2start+i*numDeriv, j)
                    * coefs.template block<1,GeoDim>(m_active(i,j),0);

            reshape.resize(GeoDim*numDeriv,1);
            m_2ndDers.template block<GeoDim*numDeriv,1>(0,j).noalias() += reshape;
        }
    }
}

template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeCurl()
{
    /// NOT TESTED FOR CORRECTNESS
    GISMO_ASSERT( ParDim==3 && ParDim == GeoDim, "Implemented only for 3D geometries in 3D spaces");
    const int numPoints=m_jacobians.cols()/ParDim;
    m_curl.resize(GeoDim,numPoints);
    for (int i=0; i<numPoints;++i)
    {
        m_curl(0,i)=gsGeometryEvaluator<T>::jacobian(i)(2,3)-gsGeometryEvaluator<T>::jacobian(i)(3,2);
        m_curl(1,i)=gsGeometryEvaluator<T>::jacobian(i)(3,1)-gsGeometryEvaluator<T>::jacobian(i)(1,3);
        m_curl(2,i)=gsGeometryEvaluator<T>::jacobian(i)(1,2)-gsGeometryEvaluator<T>::jacobian(i)(2,1);
    }
}

template <class T, int ParDim, int codim>
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

template <class T, int ParDim, int codim>
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

template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::
     transformBasisValueDivergencePreserving(
        index_t k,
        const gsMatrix<T>& allValues, // NumAct X TargetDim
        gsMatrix<T>  & trfValues) const
{
    GISMO_ASSERT(this->m_flags & NEED_MEASURE, "Need determinant and Jacobian");
    GISMO_ASSERT(GeoDim == ParDim, "Implementation assumes geometric dimention and parametric dimention are the same");
    GISMO_ASSERT(GeoDim == allValues.cols(), "Wrong format in allValues (Also assumes that target dimention is equal to geometric dimention");

    index_t TarDim = GeoDim;
    index_t numAct = allValues.rows();

    gsVector<T> componetVector;


    trfValues.setZero(numAct,TarDim);

    for (index_t comp = 0;  comp < TarDim; ++comp)
    {
        componetVector.setZero(TarDim);
        componetVector(comp) = 1.0;
        for ( index_t j=0; j< numAct ; j++ ) // for all active functions
        {
            //JS2: Im not sure if I should use the transposed or not (or inverse)
            trfValues.row(j) += ((m_jacobians.block(0, k*TarDim, TarDim,TarDim) * componetVector) * allValues(j,comp)/m_measures[k]).transpose();
        }
    }
}

template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::
    transformGradientsDivergencePreserving(
        index_t k,
        const gsMatrix<T>& allValues, // NumAct X TargetDim
        std::vector<gsMatrix<T> > const & allGrads_vec,
        std::vector<gsMatrix<T> > & trfGradsK_vec) const
{
    //Assumtions:
    //GeoDim = ParDim = TargetDim
    //numAct is equal for three components

    GISMO_ASSERT(this->m_flags & NEED_GRAD_TRANSFORM, "J^-1 not computed");
    GISMO_ASSERT(this->m_flags & NEED_2ND_DER , "Second derivatives are not computed");
    GISMO_ASSERT(GeoDim == ParDim, "Implementation assumes geometric dimention and parametric dimention are the same");


    trfGradsK_vec.resize(GeoDim);
    gsVector<index_t> numAct_vec(GeoDim);
    gsVector<T> determinantDerivative;
    determinantDerivative.setZero(GeoDim);
    std::vector<gsMatrix<T> > allGradsResize_vec;       //NumAct X ParDim
    std::vector<gsMatrix<T> > allGradsResizeGeoVar_vec; //NumAct X ParDim u1hx u1hy u1hz
    allGradsResizeGeoVar_vec.resize(GeoDim); //NumAct X ParDim u1hx u1hy u1hz
    gsMatrix<T> trfValues; //NumAct X ParDim
    index_t k1 = k*ParDim;                         //Index for Jacobian
    index_t k2 = (ParDim + (ParDim*(ParDim-1))/2); //Index second derivative matrix

    //JS2: Note tested
    if(GeoDim==2)
    {
        //JS2: CODE GOTEN FROM GEOPDES SOFTWARE

        //xu = m_jacobians(0,0+k1)
        //yu = m_jacobians(1,0+k1)
        //xv = m_jacobians(0,1+k1)
        //yv = m_jacobians(1,1+k1)
        //xuu = m_2ndDers(0+0*k2,k)
        //xvv = m_2ndDers(1+0*k2,k)
        //xuv = m_2ndDers(2+0*k2,k)
        //yuu = m_2ndDers(0+1*k2,k)
        //yvv = m_2ndDers(1+1*k2,k)
        //yuv = m_2ndDers(2+1*k2,k)
        //det = m_measures[k]

        //wh = u1h = allValuesT.row(0)
        //zh = u2h = allValuesT.row(1)
        //whu = allGradsResizeGeoVar_vec[0].row(0)
        //whv = allGradsResizeGeoVar_vec[0].row(1)
        //zhu = allGradsResizeGeoVar_vec[1].row(0)
        //zhv = allGradsResizeGeoVar_vec[1].row(1)

        numAct_vec[0] = (allGrads_vec[0].rows() / ParDim);
        numAct_vec[1] = (allGrads_vec[1].rows() / ParDim);

        T detdx = m_jacobians(1,0+k1)*(m_jacobians(1,0+k1)*m_2ndDers(1+0*k2,k) + m_jacobians(0,1+k1)*m_2ndDers(2+1*k2,k) - m_jacobians(0,0+k1)*m_2ndDers(1+1*k2,k));
        detdx  += m_jacobians(1,1+k1)*(m_jacobians(1,1+k1)*m_2ndDers(0+0*k2,k) + m_jacobians(0,0+k1)*m_2ndDers(2+1*k2,k) - m_jacobians(0,1+k1)*m_2ndDers(0+1*k2,k));
        detdx  +=-2*m_jacobians(1,0+k1)*m_jacobians(1,1+k1)*m_2ndDers(2+0*k2,k);
        detdx  /= m_measures[k];
        //detdx = (yu.*(yu.*xvv + xv.*yuv - xu.*yvv) + ...
        //         yv.*(yv.*xuu + xu.*yuv - xv.*yuu) - 2*yu.*yv.*xuv)./det;

        T detdy = m_jacobians(0,0+k1)*(m_jacobians(1,1+k1)*m_2ndDers(2+0*k2,k) + m_jacobians(0,0+k1)*m_2ndDers(1+1*k2,k) - m_jacobians(1,0+k1)*m_2ndDers(1+0*k2,k));
        detdy  += m_jacobians(0,1+k1)*(m_2ndDers(2+0*k2,k)*m_jacobians(1,0+k1) + m_2ndDers(0+1*k2,k)*m_jacobians(0,1+k1) - m_jacobians(1,1+k1)*m_2ndDers(0+0*k2,k));
        detdy  +=-2*m_jacobians(0,0+k1)*m_jacobians(0,1+k1)*m_2ndDers(2+1*k2,k);
        detdy  /= m_measures[k];
        //detdy = (xu.*(yv.*xuv + xu.*yvv - yu.*xvv) + ...
        //         xv.*(xuv.*yu + yuu.*xv - yv.*xuu) - 2*xu.*xv.*yuv)./det;

        trfGradsK_vec[0].setZero(ParDim,numAct_vec[0]);
        trfGradsK_vec[1].setZero(ParDim,numAct_vec[0]);

        trfValues.setZero(numAct_vec[0],ParDim);
        trfValues.col(0).noalias() = (m_jacobians(0,0+k1)*allValues.col(0) + m_jacobians(0,1+k1)*allValues.col(1)); //v1 = xu.*u1h + xv.*u2h = xu.*wh + xv.*zh
        trfValues.col(1).noalias() = (m_jacobians(1,0+k1)*allValues.col(0) + m_jacobians(1,1+k1)*allValues.col(1)); //v2 = yu.*u1h + yv.*u2h = yu.*wh + yv.*zh

        trfValues.transposeInPlace();  //TargetDim X NumAct
        gsMatrix<T> allValuesT =  allValues.transpose();  //TargetDim X NumAct


        gsMatrix<T> tmp1 =allGrads_vec[0].col(k);
        gsMatrix<T> tmp2 =allGrads_vec[1].col(k);
        tmp1.resize(ParDim,numAct_vec[0]);
        tmp2.resize(ParDim,numAct_vec[1]);
        allGradsResize_vec.push_back(tmp1); //NumAct X ParDim
        allGradsResize_vec.push_back(tmp2);
        //allGradsResizeGeoVar_vec[0].transposeInPlace();  //ParDim X NumAct
        //allGradsResizeGeoVar_vec[1].transposeInPlace();  //ParDim X NumAct


        //v1x // 1 X NumAct
        trfGradsK_vec[0].row(0).noalias() = (-detdx*trfValues.row(0) \
                + m_jacobians(1,1+k1)*(m_2ndDers(0+0*k2,k)*allValuesT.row(0) + m_2ndDers(2+0*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(0,0+k1)*allGradsResize_vec[0].row(0) + m_jacobians(0,1+k1)*allGradsResize_vec[1].row(0)) \
                - m_jacobians(1,0+k1)*(m_2ndDers(2+0*k2,k)*allValuesT.row(0) + m_2ndDers(1+0*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(0,0+k1)*allGradsResize_vec[0].row(1) + m_jacobians(0,1+k1)*allGradsResize_vec[1].row(1)))/(m_measures[k]*m_measures[k]);


        //v1y // 1 X NumAct
        trfGradsK_vec[0].row(1).noalias() = (-detdy*trfValues.row(0) \
                + m_jacobians(0,0+k1)*(m_2ndDers(2+0*k2,k)*allValuesT.row(0) + m_2ndDers(1+0*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(0,0+k1)*allGradsResize_vec[0].row(1) + m_jacobians(0,1+k1)*allGradsResize_vec[1].row(1)) \
                - m_jacobians(0,1+k1)*(m_2ndDers(0+0*k2,k)*allValuesT.row(0) + m_2ndDers(2+0*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(0,0+k1)*allGradsResize_vec[0].row(0) + m_jacobians(0,1+k1)*allGradsResize_vec[1].row(0)))/(m_measures[k]*m_measures[k]);

        //v2x // 1 X NumAct
        trfGradsK_vec[1].row(0).noalias() = (-detdx*trfValues.row(1) \
                + m_jacobians(1,1+k1)*(m_2ndDers(0+1*k2,k)*allValuesT.row(0) + m_2ndDers(2+1*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(1,0+k1)*allGradsResize_vec[0].row(0) + m_jacobians(1,1+k1)*allGradsResize_vec[1].row(0)) \
                - m_jacobians(1,0+k1)*(m_2ndDers(2+1*k2,k)*allValuesT.row(0) + m_2ndDers(1+1*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(1,0+k1)*allGradsResize_vec[0].row(1) + m_jacobians(1,1+k1)*allGradsResize_vec[1].row(1)))/(m_measures[k]*m_measures[k]);

        //v2y // 1 X NumAct
        trfGradsK_vec[1].row(1).noalias() = (-detdy*trfValues.row(1) \
                + m_jacobians(0,0+k1)*(m_2ndDers(2+1*k2,k)*allValuesT.row(0) + m_2ndDers(1+1*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(1,0+k1)*allGradsResize_vec[0].row(1) + m_jacobians(1,1+k1)*allGradsResize_vec[1].row(1)) \
                - m_jacobians(0,1+k1)*(m_2ndDers(0+1*k2,k)*allValuesT.row(0) + m_2ndDers(2+1*k2,k)*allValuesT.row(1) + \
                                       m_jacobians(1,0+k1)*allGradsResize_vec[0].row(0) + m_jacobians(1,1+k1)*allGradsResize_vec[1].row(0)))/(m_measures[k]*m_measures[k]);

        /*
        v1x = (-(xu.*wh + xv.*zh).*detdx +...
                    yv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu)...
                   -yu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv))./det2;

        v1y = (-(xu.*wh + xv.*zh).*detdy +...
                    xu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv)...
                   -xv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu))./det2;

        v2x = (-(yu.*wh + yv.*zh).*detdx +...
                    yv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu)...
                   -yu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv))./det2;

        v2y = (-(yu.*wh + yv.*zh).*detdy +...
                    xu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv)...
                   -xv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu))./det2;
        */

    }

    else if(GeoDim==3) //3D case
    {

        gsVector<T> xdu(3), xdv(3), xdw(3), ydu(3), ydv(3), ydw(3), zdu(3), zdv(3), zdw(3);

        for (index_t i = 0; i < GeoDim; ++i)
        {
            numAct_vec[i] = (allGrads_vec[i].rows() / ParDim);

            xdu(i) = m_2ndDers(0+0*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(3+0*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(4+0*k2,k)*m_jacInvs(i,2+k1);
            xdv(i) = m_2ndDers(3+0*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(1+0*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(5+0*k2,k)*m_jacInvs(i,2+k1);
            xdw(i) = m_2ndDers(4+0*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(5+0*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(2+0*k2,k)*m_jacInvs(i,2+k1);
            ydu(i) = m_2ndDers(0+1*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(3+1*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(4+1*k2,k)*m_jacInvs(i,2+k1);
            ydv(i) = m_2ndDers(3+1*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(1+1*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(5+1*k2,k)*m_jacInvs(i,2+k1);
            ydw(i) = m_2ndDers(4+1*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(5+1*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(2+1*k2,k)*m_jacInvs(i,2+k1);
            zdu(i) = m_2ndDers(0+2*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(3+2*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(4+2*k2,k)*m_jacInvs(i,2+k1);
            zdv(i) = m_2ndDers(3+2*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(1+2*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(5+2*k2,k)*m_jacInvs(i,2+k1);
            zdw(i) = m_2ndDers(4+2*k2,k)*m_jacInvs(i,0+k1) + m_2ndDers(5+2*k2,k)*m_jacInvs(i,1+k1) + m_2ndDers(2+2*k2,k)*m_jacInvs(i,2+k1);

            determinantDerivative[i] += xdu(i)*m_jacInvs(0,0+k1)+xdv(i)*m_jacInvs(0,1+k1)+xdw(i)*m_jacInvs(0,2+k1);
            determinantDerivative[i] += ydu(i)*m_jacInvs(1,0+k1)+ydv(i)*m_jacInvs(1,1+k1)+ydw(i)*m_jacInvs(1,2+k1);
            determinantDerivative[i] += zdu(i)*m_jacInvs(2,0+k1)+zdv(i)*m_jacInvs(2,1+k1)+zdw(i)*m_jacInvs(2,2+k1);
            determinantDerivative[i] *= m_measures[k];

            gsMatrix<T> tmp =allGrads_vec[i].col(k);
            tmp.resize(ParDim,numAct_vec[i]);
            allGradsResize_vec.push_back(tmp.transpose());


            allGradsResizeGeoVar_vec[i].setZero(numAct_vec[i], ParDim);
            allGradsResizeGeoVar_vec[i].col(0).noalias() = allGradsResize_vec[i] * m_jacInvs.col(0+k1);
            allGradsResizeGeoVar_vec[i].col(1).noalias() = allGradsResize_vec[i] * m_jacInvs.col(1+k1);
            allGradsResizeGeoVar_vec[i].col(2).noalias() = allGradsResize_vec[i] * m_jacInvs.col(2+k1);
            allGradsResizeGeoVar_vec[i].transposeInPlace();  //ParDim X NumAct
        }

        //trfValues //NumAct X ParDim  #allValues // NumAct X TargetDim
        trfValues.setZero(numAct_vec[0],ParDim);
        trfValues.col(0).noalias() = (m_jacobians(0,0+k1)*allValues.col(0) + m_jacobians(0,1+k1)*allValues.col(1) +m_jacobians(0,2+k1)*allValues.col(2))/m_measures[k];
        trfValues.col(1).noalias() = (m_jacobians(1,0+k1)*allValues.col(0) + m_jacobians(1,1+k1)*allValues.col(1) +m_jacobians(1,2+k1)*allValues.col(2))/m_measures[k];
        trfValues.col(2).noalias() = (m_jacobians(2,0+k1)*allValues.col(0) + m_jacobians(2,1+k1)*allValues.col(1) +m_jacobians(2,2+k1)*allValues.col(2))/m_measures[k];
        //v1 = 1./det.*(xu.*u1h + xv.*u2h + xw.*u3h);
        //v2 = 1./det.*(yu.*u1h + yv.*u2h + yw.*u3h);
        //v3 = 1./det.*(zu.*u1h + zv.*u2h + zw.*u3h);

        trfValues.transposeInPlace();  //TargetDim X NumAct
        gsMatrix<T> allValuesT =  allValues.transpose();  //TargetDim X NumAct

        trfGradsK_vec[0].setZero(ParDim,numAct_vec[0]);
        trfGradsK_vec[1].setZero(ParDim,numAct_vec[0]);
        trfGradsK_vec[2].setZero(ParDim,numAct_vec[0]);

        trfGradsK_vec[0].row(0).noalias() = (-determinantDerivative[0]*trfValues.row(0) \
                + xdu(0)*allValuesT.row(0) + m_jacobians(0,0+k1)*allGradsResizeGeoVar_vec[0].row(0) \
                + xdv(0)*allValuesT.row(1) + m_jacobians(0,1+k1)*allGradsResizeGeoVar_vec[1].row(0) \
                + xdw(0)*allValuesT.row(2) + m_jacobians(0,2+k1)*allGradsResizeGeoVar_vec[2].row(0))/m_measures[k];
        //*/
        trfGradsK_vec[0].row(1).noalias() = (-determinantDerivative[1]*trfValues.row(0) \
                + xdu(1)*allValuesT.row(0) + m_jacobians(0,0+k1)*allGradsResizeGeoVar_vec[0].row(1) \
                + xdv(1)*allValuesT.row(1) + m_jacobians(0,1+k1)*allGradsResizeGeoVar_vec[1].row(1) \
                + xdw(1)*allValuesT.row(2) + m_jacobians(0,2+k1)*allGradsResizeGeoVar_vec[2].row(1))/m_measures[k];

        trfGradsK_vec[0].row(2).noalias() = (-determinantDerivative[2]*trfValues.row(0) \
                + xdu(2)*allValuesT.row(0) + m_jacobians(0,0+k1)*allGradsResizeGeoVar_vec[0].row(2) \
                + xdv(2)*allValuesT.row(1) + m_jacobians(0,1+k1)*allGradsResizeGeoVar_vec[1].row(2) \
                + xdw(2)*allValuesT.row(2) + m_jacobians(0,2+k1)*allGradsResizeGeoVar_vec[2].row(2))/m_measures[k];


        trfGradsK_vec[1].row(0).noalias() = (-determinantDerivative[0]*trfValues.row(1) \
                + ydu(0)*allValuesT.row(0) + m_jacobians(1,0+k1)*allGradsResizeGeoVar_vec[0].row(0) \
                + ydv(0)*allValuesT.row(1) + m_jacobians(1,1+k1)*allGradsResizeGeoVar_vec[1].row(0) \
                + ydw(0)*allValuesT.row(2) + m_jacobians(1,2+k1)*allGradsResizeGeoVar_vec[2].row(0))/m_measures[k];

        trfGradsK_vec[1].row(1).noalias() = (-determinantDerivative[1]*trfValues.row(1) \
                + ydu(1)*allValuesT.row(0) + m_jacobians(1,0+k1)*allGradsResizeGeoVar_vec[0].row(1) \
                + ydv(1)*allValuesT.row(1) + m_jacobians(1,1+k1)*allGradsResizeGeoVar_vec[1].row(1) \
                + ydw(1)*allValuesT.row(2) + m_jacobians(1,2+k1)*allGradsResizeGeoVar_vec[2].row(1))/m_measures[k];

        trfGradsK_vec[1].row(2).noalias() = (-determinantDerivative[2]*trfValues.row(1) \
                + ydu(2)*allValuesT.row(0) + m_jacobians(1,0+k1)*allGradsResizeGeoVar_vec[0].row(2) \
                + ydv(2)*allValuesT.row(1) + m_jacobians(1,1+k1)*allGradsResizeGeoVar_vec[1].row(2) \
                + ydw(2)*allValuesT.row(2) + m_jacobians(1,2+k1)*allGradsResizeGeoVar_vec[2].row(2))/m_measures[k];


        trfGradsK_vec[2].row(0).noalias() = (-determinantDerivative[0]*trfValues.row(2) \
                + zdu(0)*allValuesT.row(0) + m_jacobians(2,0+k1)*allGradsResizeGeoVar_vec[0].row(0) \
                + zdv(0)*allValuesT.row(1) + m_jacobians(2,1+k1)*allGradsResizeGeoVar_vec[1].row(0) \
                + zdw(0)*allValuesT.row(2) + m_jacobians(2,2+k1)*allGradsResizeGeoVar_vec[2].row(0))/m_measures[k];

        trfGradsK_vec[2].row(1).noalias() = (-determinantDerivative[1]*trfValues.row(2) \
                + zdu(1)*allValuesT.row(0) + m_jacobians(2,0+k1)*allGradsResizeGeoVar_vec[0].row(1) \
                + zdv(1)*allValuesT.row(1) + m_jacobians(2,1+k1)*allGradsResizeGeoVar_vec[1].row(1) \
                + zdw(1)*allValuesT.row(2) + m_jacobians(2,2+k1)*allGradsResizeGeoVar_vec[2].row(1))/m_measures[k];

        trfGradsK_vec[2].row(2).noalias() = (-determinantDerivative[2]*trfValues.row(2) \
                + zdu(2)*allValuesT.row(0) + m_jacobians(2,0+k1)*allGradsResizeGeoVar_vec[0].row(2) \
                + zdv(2)*allValuesT.row(1) + m_jacobians(2,1+k1)*allGradsResizeGeoVar_vec[1].row(2) \
                + zdw(2)*allValuesT.row(2) + m_jacobians(2,2+k1)*allGradsResizeGeoVar_vec[2].row(2))/m_measures[k];

    }
    else
    {
        GISMO_ERROR("Geometric dimention not valid");
    }
}





}; // namespace gismo
