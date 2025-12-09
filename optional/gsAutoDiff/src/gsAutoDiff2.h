/** @file gsAutoDiff2.h

    @brief Automatic differentiation data type for C++

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include "gsAutoDiff2_fwd.h"
#include <cmath>
#include <sstream>
#include <type_traits>
#include <utility>

#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsForwardDeclarations.h>

#include <gsMatrix/gsEigenDeclarations.h>

#include <gsAutoDiff/gsVarAdaptor.h>
#include <gsAutoDiff/gsAutoDiffEigen.h>
#include <gsAutoDiff/gsAutoDiffTraits.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>


