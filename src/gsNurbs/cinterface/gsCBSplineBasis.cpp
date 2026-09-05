/** @file gsCBSplineBasis.cpp

    @brief C interface: B-spline, tensor-product B-spline and NURBS bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsNurbsBasis.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsNurbs/cinterface/gsCBSplineBasis.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCBasis * gsBSplineBasis_create(gsCKnotVector * KV)
{
    GISMO_CAPI_BEGIN
    auto * KV_ptr = RICAST_KV (KV);
    return RICAST_CB (new gsBSplineBasis<double>(*KV_ptr) );
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTensorBSplineBasis2_create(gsCKnotVector* KV1, gsCKnotVector* KV2)
{
    GISMO_CAPI_BEGIN
    auto * KV1_ptr = RICAST_KV (KV1);
    auto * KV2_ptr = RICAST_KV (KV2);
    return RICAST_CB(new gsTensorBSplineBasis<2,double>(*KV1_ptr,*KV2_ptr) );
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTensorBSplineBasis3_create(gsCKnotVector* KV1, gsCKnotVector* KV2,
                                                    gsCKnotVector* KV3)
{
    GISMO_CAPI_BEGIN
    auto * KV1_ptr = RICAST_KV (KV1);
    auto * KV2_ptr = RICAST_KV (KV2);
    auto * KV3_ptr = RICAST_KV (KV3);
    return RICAST_CB(new gsTensorBSplineBasis<3,double>(*KV1_ptr,*KV2_ptr,*KV3_ptr));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTensorBSplineBasis4_create(gsCKnotVector* KV1, gsCKnotVector* KV2,
                                                    gsCKnotVector* KV3, gsCKnotVector* KV4)
{
    GISMO_CAPI_BEGIN
    auto * KV1_ptr = RICAST_KV (KV1);
    auto * KV2_ptr = RICAST_KV (KV2);
    auto * KV3_ptr = RICAST_KV (KV3);
    auto * KV4_ptr = RICAST_KV (KV4);
    return RICAST_CB(new gsTensorBSplineBasis<4,double>(*KV1_ptr,*KV2_ptr,*KV3_ptr,*KV4_ptr));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis * gsNurbsBasis_create(gsCBasis * b, gsCMatrix * weights)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsBSplineBasis<double>* >(b);
    auto * w = RICAST_M(weights);
    return RICAST_CB(new gsNurbsBasis<double>(basis_ptr,*w));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTensorNurbsBasis2_create(gsCBasis* b, gsCMatrix * weights)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<2,double>* >(b);
    auto * w = RICAST_M(weights);
    return RICAST_CB(new  gsTensorNurbsBasis<2,double>(basis_ptr,*w));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTensorNurbsBasis3_create(gsCBasis* b, gsCMatrix * weights)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<3,double>* >(b);
    auto * w = RICAST_M(weights);
    return RICAST_CB(new  gsTensorNurbsBasis<3,double>(basis_ptr,*w));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCBasis* gsTensorNurbsBasis4_create(gsCBasis* b, gsCMatrix * weights)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<4,double>* >(b);
    auto * w = RICAST_M(weights);
    return RICAST_CB(new  gsTensorNurbsBasis<4,double>(basis_ptr,*w));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCKnotVector * gsBSplineBasis_knots(gsCBasis * b)
{
    GISMO_CAPI_BEGIN
    gsKnotVector<double> * KV= &reinterpret_cast< gsBSplineBasis<double>* >(b)->knots();
    return reinterpret_cast<gsCKnotVector*>(KV);
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCKnotVector * gsTensorBSplineBasis_knots(gsCBasis * b, int dir)
{
    GISMO_CAPI_BEGIN
    gsKnotVector<double> * KV=NULL;
    GISMO_ASSERT(RICAST_B(b)->domainDim()>=dir,"gsTensorBSplineBasis_knots: dir out of range");
    switch (RICAST_B(b)->domainDim())
    {
        case 2:
            KV = &reinterpret_cast< gsTensorBSplineBasis<2,double>* >(b)->knots(dir);
        case 3:
            KV = &reinterpret_cast< gsTensorBSplineBasis<3,double>* >(b)->knots(dir);
        case 4:
            KV = &reinterpret_cast< gsTensorBSplineBasis<4,double>* >(b)->knots(dir);
    }

    return reinterpret_cast<gsCKnotVector*>(KV);
    GISMO_CAPI_END(NULL)
}

#ifdef __cplusplus
}
#endif
