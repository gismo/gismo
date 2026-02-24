/** @file gsClosestPoint.cpp

    @brief Provides instantization of gsClosestPoint functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D. Mokris
*/

#include <gsCore/gsTemplateTools.h>

#include <gsParasolid/gsClosestPoint.h>
#include <gsParasolid/gsClosestPoint.hpp>

namespace gismo
{
namespace extensions
{
TEMPLATE_INST PK_VECTOR_t
gsPK_VECTOR(const gsVector<real_t, 3>& gsVector);

TEMPLATE_INST bool
gsClosestParam(const gsTensorBSpline<2, real_t>& gsBSurf,
               const gsVector<real_t, 3>& gsPoint,
               gsVector<real_t, 2>& gsResult);

TEMPLATE_INST bool
gsClosestParam(const gsTensorBSpline<2, real_t>& gsBSurf,
               const gsMatrix<real_t>& gsPoints,
               gsMatrix<real_t>& gsResults);

} // namespace extensions
} // namespace gismo
