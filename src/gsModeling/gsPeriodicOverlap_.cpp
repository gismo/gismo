/** @file gsPeriodicOverlap_.cpp

    @brief Instatiation of the gsPeriodicOverlap class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsPeriodicOverlap.h>
#include <gsModeling/gsPeriodicOverlap.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsPeriodicOverlap<real_t>;

} // namespace gismo
