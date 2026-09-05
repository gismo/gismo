/** @file gsCHTensorBasis.cpp

    @brief C interface: hierarchical (HB, THB) spline bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// Header-only builds: the B-spline basis header must come first, it
// completes gsTensorBSplineBasis before gsSurfMesh.hpp needs it
#include <gsNurbs/gsBSplineBasis.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsHDomainIterator.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsHSplines/cinterface/gsCHTensorBasis.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

template <int dim>
void gsHTensorBasis_elements_into_impl(gsCBasis * b, bool getKnotBoxes,
                                                     bool getIndexBoxes,
                                                     bool getLevels,
                                                     gsCMatrix*    knotBoxes,
                                                     gsCMatrixInt* indexBoxes,
                                                     gsCVectorInt* levels)
{
    int N   = RICAST_B(b)->numElements();

    auto * el = RICAST_M(knotBoxes);
    auto * bx = RICAST_Mi(indexBoxes);
    auto * lv = RICAST_Vi(levels);

    if (getKnotBoxes) el->resize(dim,2*N);
    if (getIndexBoxes) bx->resize(2*dim+1,N);
    if (getLevels) lv->resize(N);

    auto domain = RICAST_B(b)->domain();
    auto domIt  = domain->beginAll();
    auto domEnd = domain->endAll();

    GISMO_ENSURE((dynamic_cast<gsHDomainIterator<double,dim> *>(domIt.get())),"The domain iterator is not a hierarchical domain iterator");
    gsHTensorBasis<dim,double> * basis = dynamic_cast<gsHTensorBasis<dim,double> *>(RICAST_B(b));

    int id=0;
    gsVector<double,dim> low, upp;
    for (; domIt<domEnd; ++domIt, ++id)
    {
        gsHDomainIterator<double,dim> * domItH = dynamic_cast<gsHDomainIterator<double,dim> *>(domIt.get());
        low = domIt.lowerCorner();
        upp = domIt.upperCorner();
        if (getKnotBoxes)
        {
            el->col(2*id) = low;
            el->col(2*id+1) = upp;
        }
        if (getLevels)
        {
            lv->at(id) = domItH->getLevel();
        }
        if (getIndexBoxes)
        {
            for(int j = 0; j < dim;j++)
            {
                // Convert the parameter coordinates to (unique) knot indices
                const gsKnotVector<real_t> & kv = basis->tensorLevel(domItH->getLevel()).knots(j);
                int k1 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd(),
                                        low[j] ) - 1).uIndex();
                int k2 = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd()+1,
                                        upp[j] ) - 1).uIndex();

                // Trivial cells trigger some refinement
                if ( k1 == k2)
                {
                    if (0!=k1) {--k1;}
                    ++k2;
                }

                // Store the data...
                (*bx)(0,id) = domItH->getLevel();
                (*bx)(1+j,id) = k1;
                (*bx)(1+j+dim,id) = k2;
            }
        }
    }
}

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCBasis* gsTHBSplineBasis1_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<1,double>* >(b);
    return RICAST_CB(new gsTHBSplineBasis<1,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTHBSplineBasis2_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<2,double>* >(b);
    return RICAST_CB(new  gsTHBSplineBasis<2,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTHBSplineBasis3_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<3,double>* >(b);
    return RICAST_CB(new  gsTHBSplineBasis<3,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTHBSplineBasis4_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<4,double>* >(b);
    return RICAST_CB(new  gsTHBSplineBasis<4,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsHBSplineBasis1_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<1,double>* >(b);
    return RICAST_CB(new gsHBSplineBasis<1,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsHBSplineBasis2_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<2,double>* >(b);
    return RICAST_CB(new  gsHBSplineBasis<2,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsHBSplineBasis3_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<3,double>* >(b);
    return RICAST_CB(new  gsHBSplineBasis<3,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsHBSplineBasis4_create(gsCBasis* b)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<4,double>* >(b);
    return RICAST_CB(new  gsHBSplineBasis<4,double>(*basis_ptr,false));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsHTensorBasis_elements_into(gsCBasis * b, bool getKnotBoxes,
                                                             bool getIndexBoxes,
                                                             bool getLevels,
                                                             gsCMatrix*    knotBoxes,
                                                             gsCMatrixInt* indexBoxes,
                                                             gsCVectorInt* levels)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            gsHTensorBasis_elements_into_impl<1>(b,getKnotBoxes,getIndexBoxes,getLevels,knotBoxes,indexBoxes,levels);
            break;
        case 2:
            gsHTensorBasis_elements_into_impl<2>(b,getKnotBoxes,getIndexBoxes,getLevels,knotBoxes,indexBoxes,levels);
            break;
        case 3:
            gsHTensorBasis_elements_into_impl<3>(b,getKnotBoxes,getIndexBoxes,getLevels,knotBoxes,indexBoxes,levels);
            break;
        case 4:
            gsHTensorBasis_elements_into_impl<4>(b,getKnotBoxes,getIndexBoxes,getLevels,knotBoxes,indexBoxes,levels);
            break;
        default:
            GISMO_ERROR("gsHTensorBasis_elements_into: Dimension not supported");
    }
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT int gsHTensorBasis_numLevels(gsCBasis * b)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            return reinterpret_cast< gsHTensorBasis<1,double>* >(b)->numLevels();
        case 2:
            return reinterpret_cast< gsHTensorBasis<2,double>* >(b)->numLevels();
        case 3:
            return reinterpret_cast< gsHTensorBasis<3,double>* >(b)->numLevels();
        case 4:
            return reinterpret_cast< gsHTensorBasis<4,double>* >(b)->numLevels();
        default:
            GISMO_ERROR("gsHTensorBasis_numLevels: domainDim not supported");
    }
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsHTensorBasis_maxLevel(gsCBasis * b)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            return reinterpret_cast< gsHTensorBasis<1,double>* >(b)->maxLevel();
        case 2:
            return reinterpret_cast< gsHTensorBasis<2,double>* >(b)->maxLevel();
        case 3:
            return reinterpret_cast< gsHTensorBasis<3,double>* >(b)->maxLevel();
        case 4:
            return reinterpret_cast< gsHTensorBasis<4,double>* >(b)->maxLevel();
        default:
            GISMO_ERROR("gsHTensorBasis_maxLevel: domainDim not supported");
    }
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsHTensorBasis_levelOf(gsCBasis * b, int i)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            return reinterpret_cast< gsHTensorBasis<1,double>* >(b)->levelOf(i);
        case 2:
            return reinterpret_cast< gsHTensorBasis<2,double>* >(b)->levelOf(i);
        case 3:
            return reinterpret_cast< gsHTensorBasis<3,double>* >(b)->levelOf(i);
        case 4:
            return reinterpret_cast< gsHTensorBasis<4,double>* >(b)->levelOf(i);
        default:
            GISMO_ERROR("gsHTensorBasis_levelOf: domainDim not supported");
    }
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT int gsHTensorBasis_getLevelAtPoint(gsCBasis * b, gsCMatrix * Pt)
{
    GISMO_CAPI_BEGIN
    auto * m = RICAST_M(Pt);
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            return reinterpret_cast< gsHTensorBasis<1,double>* >(b)->getLevelAtPoint(*m);
        case 2:
            return reinterpret_cast< gsHTensorBasis<2,double>* >(b)->getLevelAtPoint(*m);
        case 3:
            return reinterpret_cast< gsHTensorBasis<3,double>* >(b)->getLevelAtPoint(*m);
        case 4:
            return reinterpret_cast< gsHTensorBasis<4,double>* >(b)->getLevelAtPoint(*m);
        default:
            GISMO_ERROR("gsHTensorBasis_getLevelAtPoint: domainDim not supported");
    }
    GISMO_CAPI_END(-1)
}

GISMO_EXPORT gsCBasis * gsHTensorBasis_tensorLevel(gsCBasis * b, int l)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            return RICAST_CB(new  gsBSplineBasis<double>(reinterpret_cast< gsHTensorBasis<1,double>* >(b)->tensorLevel(l)));
        case 2:
            return RICAST_CB(new  gsTensorBSplineBasis<2,double>(reinterpret_cast< gsHTensorBasis<2,double>* >(b)->tensorLevel(l)));
        case 3:
            return RICAST_CB(new  gsTensorBSplineBasis<3,double>(reinterpret_cast< gsHTensorBasis<3,double>* >(b)->tensorLevel(l)));
        case 4:
            return RICAST_CB(new  gsTensorBSplineBasis<4,double>(reinterpret_cast< gsHTensorBasis<4,double>* >(b)->tensorLevel(l)));
        default:
            GISMO_ERROR("gsHTensorBasis_tensorLevel: domainDim not supported");
    }
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsHTensorBasis_treeLeafSize(gsCBasis * b)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            reinterpret_cast< gsHTensorBasis<1,double>* >(b)->tree().leafSize();
        case 2:
            reinterpret_cast< gsHTensorBasis<2,double>* >(b)->tree().leafSize();
        case 3:
            reinterpret_cast< gsHTensorBasis<3,double>* >(b)->tree().leafSize();
        case 4:
            reinterpret_cast< gsHTensorBasis<4,double>* >(b)->tree().leafSize();
        default:
            GISMO_ERROR("gsHTensorBasis_treeLeaveSize: domainDim not supported");
    }
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsHTensorBasis_treePrintLeaves(gsCBasis * b)
{
    GISMO_CAPI_BEGIN
    switch (RICAST_B(b)->domainDim())
    {
        case 1:
            reinterpret_cast< gsHTensorBasis<1,double>* >(b)->tree().printLeaves();
        case 2:
            reinterpret_cast< gsHTensorBasis<2,double>* >(b)->tree().printLeaves();
        case 3:
            reinterpret_cast< gsHTensorBasis<3,double>* >(b)->tree().printLeaves();
        case 4:
            reinterpret_cast< gsHTensorBasis<4,double>* >(b)->tree().printLeaves();
        default:
            GISMO_ERROR("gsHTensorBasis_treePrintLeaves: domainDim not supported");
    }
    GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif
