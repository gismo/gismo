/** @file gsRationalTHBSpline.h

    @brief Represents a rational truncated hierarchical B-Spline patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsRationalTHBSpline.h>
#include <gsHSplines/gsRationalTHBSpline.hpp>

namespace gismo
{
#define INST(D) CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSpline<D,real_t> >;
GISMO_DIM_FOREACH(INST)
#undef INST

} // namespace gismo
