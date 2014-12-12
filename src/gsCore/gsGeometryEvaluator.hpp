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

// Geometry transformation 
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


//////////////////////////////////////////////////
//////////////////////////////////////////////////


template <class T, int ParDim, int codim>
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


template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::computeValues()
{
    const gsMatrix<T> & coefs = m_geo.coefs();

    m_values.resize(coefs.cols(), m_numPts);

    for (index_t j=0; j < m_numPts; ++j) // for all evaluation points
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
void gsGenericGeometryEvaluator<T,ParDim,codim>::
     transformValuesHdiv(
        index_t k,
        const gsMatrix<T>& allValues, // NumAct X TargetDim
        gsMatrix<T>  & trfValues) const
{


}

template <class T, int ParDim, int codim>
void gsGenericGeometryEvaluator<T,ParDim,codim>::
    transformGradsHdiv(
        index_t k,
        const gsMatrix<T>& allValues, // NumAct X TargetDim
        std::vector<gsMatrix<T> > const & allGrads_vec,
        std::vector<gsMatrix<T> > & trfGradsK_vec) const
{


}


/*
template <class T, int ParDim, int codim>
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
*/

}; // namespace gismo
