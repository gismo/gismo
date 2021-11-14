/** @file gsCompositeBasis_.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>

#include <gsMSplines/gsDPatch.h>
#include <gsMSplines/gsDPatch.hpp>

#include <gsMSplines/gsDPatchGeom.h>
#include <gsMSplines/gsDPatchGeom.hpp>

#include <gsMSplines/gsDPatchBasis.h>
#include <gsMSplines/gsDPatchBasis.hpp>

namespace gismo
{

    CLASS_TEMPLATE_INST gsDPatch<2,real_t> ;
    CLASS_TEMPLATE_INST gsDPatchGeom<2,real_t> ;
    CLASS_TEMPLATE_INST gsDPatchBasis<2,real_t> ;

} // end namespace gismo
