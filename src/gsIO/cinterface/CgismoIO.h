/** @file CgismoIO.h

    @brief C interface umbrella of the gsIO module (file reading).

    Include only what you use: this header exposes the C API of gsIO
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

#include <gsIO/cinterface/gsCReadFile.h>
