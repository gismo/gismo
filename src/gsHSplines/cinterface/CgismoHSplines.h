/** @file CgismoHSplines.h

    @brief C interface umbrella of the gsHSplines module (quasi-interpolation).

    Include only what you use: this header exposes the C API of gsHSplines
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

#include <gsHSplines/cinterface/gsCQuasiInterpolate.h>
#include <gsHSplines/cinterface/gsCHTensorBasis.h>
#include <gsHSplines/cinterface/gsCHBSpline.h>
