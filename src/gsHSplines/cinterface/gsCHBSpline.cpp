/** @file gsCHBSpline.cpp

    @brief C interface: hierarchical (HB, THB) spline geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

// Header-only builds: the B-spline basis header must come first, it
// completes gsTensorBSplineBasis before gsSurfMesh.hpp needs it
#include <gsNurbs/gsBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsHSplines/cinterface/gsCHBSpline.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCGeometry* gsTHBSpline1_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTHBSplineBasis<1,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTHBSpline<1,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTHBSpline2_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTHBSplineBasis<2,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTHBSpline<2,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTHBSpline3_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTHBSplineBasis<3,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTHBSpline<3,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTHBSpline4_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTHBSplineBasis<4,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTHBSpline<4,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsHBSpline1_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsHBSplineBasis<1,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsHBSpline<1,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsHBSpline2_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsHBSplineBasis<2,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsHBSpline<2,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsHBSpline3_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsHBSplineBasis<3,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsHBSpline<3,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsHBSpline4_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsHBSplineBasis<4,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsHBSpline<4,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

#ifdef __cplusplus
}
#endif
