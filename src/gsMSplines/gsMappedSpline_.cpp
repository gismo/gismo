/** @file gsMappedSpline_.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsMSplines/gsMappedSpline.h>
#include <gsMSplines/gsMappedSpline.hpp>

namespace gismo
{

#define INST(D) CLASS_TEMPLATE_INST gsMappedSpline<D,real_t>;
GISMO_DIM_FOREACH(INST)
#undef INST

} // end namespace gismo
