/** @file gsCHBSpline.h

    @brief C interface: hierarchical (HB, THB) spline geometries.

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

    GISMO_EXPORT gsCGeometry* gsTHBSpline1_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTHBSpline2_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTHBSpline3_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsTHBSpline4_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsHBSpline1_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsHBSpline2_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsHBSpline3_create(gsCBasis* b, gsCMatrix * coef);
    GISMO_EXPORT gsCGeometry* gsHBSpline4_create(gsCBasis* b, gsCMatrix * coef);

#ifdef __cplusplus
}
#endif
