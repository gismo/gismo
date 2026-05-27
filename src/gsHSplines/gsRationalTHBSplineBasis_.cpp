/** @file gsRationalTHBSplineBasis.h

    @brief Provides declaration of RationalTHBSplineBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsDimMacro.h>

#include <gsHSplines/gsRationalTHBSplineBasis.h>
#include <gsHSplines/gsRationalTHBSplineBasis.hpp>

namespace gismo
{

#define INST(D) CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSplineBasis<D,real_t> >;
GISMO_DIM_FOREACH(INST)
#undef INST

} // namespace gismo
