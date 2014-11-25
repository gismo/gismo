/** @file gsGeometryEvaluator.h

    @brief Provides implementation of BasisEvaluator class.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#pragma once

#include <gsCore/gsBasisEvaluator.h>
#include <gsCore/gsGeoTransform.hpp>
#include <gsCore/gsGeometry.h>

namespace gismo {

template <typename T, int ParDim, int TarDim, int AmbDim>
gsBasisEvaluator<T> * makeBasisEvaluatorWithDims (const std::vector<gsBasis<T> *> &basis, unsigned flags, ValueTransformationType geoTrans )
{
    gsBasis<T> * basisPtr[TarDim];
    for (int i=0; i<TarDim;++i)
        basisPtr[i]=basis[i];
    switch (geoTrans)
    {
    case NO_TRANSFORMATION:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoNoTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,flags);
    case INVERSE_COMPOSITION:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoGradPreservingTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,flags);
    case DIV_CONFORMING:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoDivPreservingTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,flags);
    case CURL_CONFORMING:
        return new gsGenericBasisEvaluator<T,ParDim,TarDim,gsGeoCurlPreservingTransform<T,ParDim,TarDim,AmbDim> >(basisPtr,flags);
    default:
        gsWarn<<"I do not know this way to transform to the physical domain!\nUse constants from gsCore/gsBasisEvaluator.h.\n";
        return NULL;
    }
}


template <typename T, int ParDim, int TarDim>
gsBasisEvaluator<T> * makeBasisEvaluatorWithParAndTarDim (const std::vector<gsBasis<T> *> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{

    if (geo == NULL)
    {
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,ParDim>( basis, flags,  geoTrans);
    }
    else if( geo->coefDim() < ParDim)
    {
        gsWarn<<"Impossible to make basis evaluator with a parametrization from R^n to R^m with m<n.\n";
        return NULL;
    }
    else switch (geo->coefDim())
    {
    case 1:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,1>( basis, flags, geoTrans);
    case 2:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,2>( basis, flags, geoTrans);
    case 3:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,3>( basis, flags, geoTrans);
    case 4:
        return makeBasisEvaluatorWithDims<T,ParDim,TarDim,4>( basis, flags, geoTrans);
    default:
        gsWarn<<"Cannot construct basis evaluator for basis on R^n n>=5\nuse the appropriate instantiation of gsGenericBasisEvaluator!\n";
        return NULL;
    }
}


template <typename T, int TarDim>
gsBasisEvaluator<T> * makeBasisEvaluatorWithTarDim (const std::vector<gsBasis<T> *> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    switch (basis[0]->dim())
    {
    case 1:
        return makeBasisEvaluatorWithParAndTarDim<T,1,TarDim>( basis, flags, geo,  geoTrans);
    case 2:
        return makeBasisEvaluatorWithParAndTarDim<T,2,TarDim>( basis, flags, geo,  geoTrans);
    case 3:
        return makeBasisEvaluatorWithParAndTarDim<T,3,TarDim>( basis, flags, geo,  geoTrans);
    case 4:
        return makeBasisEvaluatorWithParAndTarDim<T,4,TarDim>( basis, flags, geo,  geoTrans);
    default:
        gsWarn<<"Cannot construct basis evaluator for basis on R^n n>=5\nuse the appropriate instantiation of gsGenericBasisEvaluator!\n";
        return NULL;
    }
}


template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator ( const gsBasis<T> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    std::vector<gsBasis<T>* >  basis_Array(1);
    basis_Array[0]= (gsBasis<T>*)&basis;
    switch (basis.dim())
    {
    case 1:
        return makeBasisEvaluatorWithParAndTarDim<T,1,1>(basis_Array,flags,geo,geoTrans);
    case 2:
        return makeBasisEvaluatorWithParAndTarDim<T,2,1>(basis_Array,flags,geo,geoTrans);
    case 3:
        return makeBasisEvaluatorWithParAndTarDim<T,3,1>(basis_Array,flags,geo,geoTrans);
    case 4:
        return makeBasisEvaluatorWithParAndTarDim<T,4,1>(basis_Array,flags,geo,geoTrans);
    default:
        gsWarn<<"Cannot construct basis evaluator for basis on R^n n>=5\nuse the appropriate instantiation of gsGenericBasisEvaluator!\n";
        return NULL;
    }
}

template <typename T>
gsBasisEvaluator<T> * makeBasisEvaluator (const std::vector<gsBasis<T> *> &basis, unsigned flags, const gsGeometry<T> *geo, ValueTransformationType geoTrans )
{
    switch (basis.size())
    {
    case 0:
        gsWarn<<"Cannot make evaluator without basis!\n";
        return NULL;
    case 1:
        return makeBasisEvaluatorWithTarDim<T,1>(basis,flags,geo,geoTrans);
    case 2:
        return makeBasisEvaluatorWithTarDim<T,2>(basis,flags,geo,geoTrans);
    case 3:
        return makeBasisEvaluatorWithTarDim<T,3>(basis,flags,geo,geoTrans);
    case 4:
        return makeBasisEvaluatorWithTarDim<T,4>(basis,flags,geo,geoTrans);
    default:
        gsWarn<<"Cannot make evaluator with target dimension >=4.\nUse the appropriate instantiation of gsGenericEvaluator.\n";
        return NULL;
    }
}


template <typename T, int ParDim, int TarDim, typename geometryTransform >
gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::gsGenericBasisEvaluator( gsBasis<T> *(&basis)[TarDim], unsigned flags)
    : gsBasisEvaluator<T>(flags)
{
    for (int i=0;i<TarDim;++i)
        m_basis[i]=basis[i];
    setFlags(flags);
    m_parDim=ParDim;
    m_tarDim=TarDim;
    m_active_shift[0]=0;
    for (int i =1;i<TarDim;++i)
        m_active_shift[i]=m_active_shift[i-1]+m_basis[i-1]->size();
    m_spaceDim=m_active_shift[TarDim-1]+m_basis[TarDim-1]->size();
}

template <typename T, int ParDim, int TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::setFlags   ( unsigned newFlags)
{
    m_flags=geometryTransform::addAuxiliaryFlags(newFlags);
    m_geo_flags=geometryTransform::getGeometryFlags(m_flags);

    if (m_flags & NEED_2ND_DER)
        m_max_deriv=2;
    else if (m_flags & (NEED_GRAD | NEED_JACOBIAN) )
        m_max_deriv=1;
    else if (m_flags & NEED_VALUE)
        m_max_deriv=0;
    else
        m_max_deriv=-1;
}




template <typename T, int ParDim, int TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::evaluateAt ( const gsMatrix<T> &points)
{
    int                 active_num[TarDim];
    gsMatrix<unsigned>  active[TarDim];
    unsigned            tot_active(0);

    for (int i =0;i<TarDim;++i)
    {
        m_basis[i]->active_into(points.col(0),active[i]);
        active_num[i]=active[i].rows();
        tot_active+=active_num[i];
        m_basis[i]->evalAllDers_into(points, m_max_deriv, m_basis_vals[i]);
    }
    m_actives.resize(tot_active,1);

    int size;
    int start = tot_active;
    for (int i =TarDim-1;i>=0;--i)
    {
        size  = active[i].rows();
        start -=size;
        m_actives.block(start,0,size,1)=(active[i].array()+m_active_shift[i]);
    }

    if (this->m_flags & NEED_VALUE)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeValues(this, NULL, m_basis_vals, active_num, m_values);
    if (this->m_flags & NEED_GRAD)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeGrads(this, NULL, m_basis_vals, active_num, m_derivs);
    if (this->m_flags & NEED_JACOBIAN)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeJacobians(this, NULL, m_basis_vals, active_num, m_jacobians);
    if (this->m_flags & NEED_DIV)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeDivs(this,  NULL, m_basis_vals, active_num,m_divs);
    if (this->m_flags & NEED_CURL)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeCurls(this, NULL, m_basis_vals, active_num, m_curls);
    if (this->m_flags & NEED_2ND_DER)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeSecDers(this, NULL, m_basis_vals, active_num, m_2ndDers);
    if (this->m_flags & NEED_LAPLACIAN)
        gsGeoNoTransform<T,ParDim,TarDim,ParDim>::computeLaplacians(this, NULL, m_basis_vals, active_num, m_laps);
}


template <typename T, int ParDim, int TarDim, typename geometryTransform >
void gsGenericBasisEvaluator<T,ParDim,TarDim,geometryTransform>::evaluateAt ( const gsMatrix<T> &points, const gsGeometryEvaluator<T> &geoEval)
{
    int                 active_num[TarDim];
    gsMatrix<unsigned>  active[TarDim];
    unsigned            tot_active(0);

    for (int i =0;i<TarDim;++i)
    {
        m_basis[i]->active_into(points.col(0),active[i]);
        active_num[i]=active[i].rows();
        tot_active+=active_num[i];
        m_basis[i]->evalAllDers_into(points, m_max_deriv, m_basis_vals[i]);
    }
    m_actives.resize(tot_active,1);

    int size;
    int start = tot_active;
    for (int i =TarDim-1;i>=0;--i)
    {
        size  = active[i].rows();
        start -=size;
        m_actives.block(start,0,size,1)=(active[i].array()+m_active_shift[i]);
    }

    if (this->m_flags & NEED_VALUE)
        geometryTransform::computeValues(this, &geoEval, m_basis_vals, active_num, m_values);
    if (this->m_flags & NEED_GRAD)
        geometryTransform::computeGrads(this, &geoEval, m_basis_vals, active_num, m_derivs);
    if (this->m_flags & NEED_JACOBIAN)
        geometryTransform::computeJacobians(this, &geoEval, m_basis_vals, active_num, m_jacobians);
    if (this->m_flags & NEED_DIV)
        geometryTransform::computeDivs(this,  &geoEval, m_basis_vals, active_num, m_divs);
    if (this->m_flags & NEED_CURL)
        geometryTransform::computeCurls(this, &geoEval, m_basis_vals, active_num,m_curls);
    if (this->m_flags & NEED_2ND_DER)
        geometryTransform::computeSecDers(this, &geoEval, m_basis_vals, active_num,m_2ndDers);
    if (this->m_flags & NEED_LAPLACIAN)
        geometryTransform::computeLaplacians(this, &geoEval, m_basis_vals, active_num,m_laps);
}


} // namespace gismo


