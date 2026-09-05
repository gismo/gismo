/** @file gsCHTensorBasis.h

    @brief C interface: hierarchical (HB, THB) spline bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCore/cinterface/gsCBasis.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

    GISMO_EXPORT gsCBasis* gsTHBSplineBasis1_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsTHBSplineBasis2_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsTHBSplineBasis3_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsTHBSplineBasis4_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsHBSplineBasis1_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsHBSplineBasis2_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsHBSplineBasis3_create(gsCBasis* basis);
    GISMO_EXPORT gsCBasis* gsHBSplineBasis4_create(gsCBasis* basis);
    GISMO_EXPORT void gsHTensorBasis_elements_into(gsCBasis * b, bool getKnotBoxes,
                                                                 bool getIndexBoxes,
                                                                 bool getLevels,
                                                                 gsCMatrix*    knotBoxes,
                                                                 gsCMatrixInt* indexBoxes,
                                                                 gsCVectorInt* levels);
    GISMO_EXPORT int gsHTensorBasis_numLevels(gsCBasis * b);
    GISMO_EXPORT int gsHTensorBasis_maxLevel(gsCBasis * b);
    GISMO_EXPORT int gsHTensorBasis_levelOf(gsCBasis * b, int i);
    GISMO_EXPORT int gsHTensorBasis_getLevelAtPoint(gsCBasis * b, gsCMatrix * Pt);
    GISMO_EXPORT gsCBasis * gsHTensorBasis_tensorLevel(gsCBasis * b, int l);
    GISMO_EXPORT void gsHTensorBasis_treeLeafSize(gsCBasis * b);
    GISMO_EXPORT void gsHTensorBasis_treePrintLeaves(gsCBasis * b);

#ifdef __cplusplus
}
#endif
