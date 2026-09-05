/** @file gsCBSplineBasis.h

    @brief C interface: B-spline, tensor-product B-spline and NURBS bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCore/cinterface/gsCBasis.h>
#include <gsDomain/cinterface/gsCKnotVector.h>

#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCBasis* gsBSplineBasis_create(gsCKnotVector* knots);
    GISMO_EXPORT gsCBasis* gsTensorBSplineBasis2_create(gsCKnotVector* KV1, gsCKnotVector* KV2);
    GISMO_EXPORT gsCBasis* gsTensorBSplineBasis3_create(gsCKnotVector* KV1, gsCKnotVector* KV2, gsCKnotVector* KV3);
    GISMO_EXPORT gsCBasis* gsTensorBSplineBasis4_create(gsCKnotVector* KV1, gsCKnotVector* KV2, gsCKnotVector* KV3, gsCKnotVector* KV4);
    GISMO_EXPORT gsCBasis * gsNurbsBasis_create(gsCBasis * b, gsCMatrix * weights);
    GISMO_EXPORT gsCBasis* gsTensorNurbsBasis2_create(gsCBasis* b, gsCMatrix * weights);
    GISMO_EXPORT gsCBasis* gsTensorNurbsBasis3_create(gsCBasis* b, gsCMatrix * weights);
    GISMO_EXPORT gsCBasis* gsTensorNurbsBasis4_create(gsCBasis* b, gsCMatrix * weights);
    GISMO_EXPORT gsCKnotVector * gsBSplineBasis_knots(gsCBasis * b);
    GISMO_EXPORT gsCKnotVector * gsTensorBSplineBasis_knots(gsCBasis * b, int dir);

#ifdef __cplusplus
}
#endif
