/** @file gsCBSpline.cpp

    @brief C interface: B-spline, tensor-product B-spline and NURBS geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/cinterface/gsCBSpline.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

using namespace gismo;

#ifdef __cplusplus
extern "C"
{
#endif

GISMO_EXPORT gsCGeometry* gsBSpline_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsBSplineBasis<double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsBSpline<double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorBSpline2_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    gsTensorBSplineBasis<2,double>* basis_ptr = reinterpret_cast< gsTensorBSplineBasis<2,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTensorBSpline<2,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorBSpline3_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    gsTensorBSplineBasis<3,double>* basis_ptr = reinterpret_cast< gsTensorBSplineBasis<3,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTensorBSpline<3,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorBSpline4_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorBSplineBasis<4,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTensorBSpline<4,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsNurbs_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsNurbsBasis<double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsNurbs<double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorNurbs2_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    gsTensorNurbsBasis<2,double>* basis_ptr = reinterpret_cast< gsTensorNurbsBasis<2,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTensorNurbs<2,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorNurbs3_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    gsTensorNurbsBasis<3,double>* basis_ptr = reinterpret_cast< gsTensorNurbsBasis<3,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTensorNurbs<3,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorNurbs4_create(gsCBasis* b, gsCMatrix * coefs)
{
    GISMO_CAPI_BEGIN
    auto * basis_ptr = reinterpret_cast< gsTensorNurbsBasis<4,double>* >(b);
    auto * m = RICAST_M(coefs);
    return RICAST_CG(new gsTensorNurbs<4,double>(*basis_ptr,*m));
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT gsCGeometry* gsTensorBSpline2_slice(gsCGeometry * g, int direction, double parameter)
{
    GISMO_CAPI_BEGIN
    auto * g_ptr = reinterpret_cast< gsTensorBSpline<2,double>* >(g);
    typedef typename gsTensorBSpline<2,double>::BoundaryGeometryType GeometryBdr;
    GeometryBdr * bdr = new GeometryBdr();
    g_ptr->slice(direction, parameter, *bdr);
    return RICAST_CG(bdr);
    GISMO_CAPI_END(NULL)
}

#ifdef __cplusplus
}
#endif
