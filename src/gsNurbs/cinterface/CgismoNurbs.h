/** @file CgismoNurbs.h

    @brief C interface umbrella of the gsNurbs module (NURBS creation).

    Include only what you use: this header exposes the C API of gsNurbs
    without dragging in the whole library. <gsCInterface/Cgismo.h>
    (generated) includes the umbrellas of all modules.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCInterface/gsCTypes.h>
#include <gsCInterface/gsCError.h>

#include <gsNurbs/cinterface/gsCNurbsCreator.h>
#include <gsNurbs/cinterface/gsCBSplineBasis.h>
#include <gsNurbs/cinterface/gsCBSpline.h>
