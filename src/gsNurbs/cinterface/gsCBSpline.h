/** @file gsCBSpline.h

    @brief C interface: B-spline, tensor-product B-spline and NURBS geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCore/cinterface/gsCGeometry.h>

#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCGeometry* gsBSpline_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorBSpline2_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorBSpline3_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorBSpline4_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsNurbs_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorNurbs2_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorNurbs3_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorNurbs4_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTensorBSpline2_slice(gsCGeometry * g, int direction, double parameter);

#ifdef __cplusplus
}
#endif
