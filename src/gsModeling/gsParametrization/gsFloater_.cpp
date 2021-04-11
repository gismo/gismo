/** @file gsFloater_.cpp

    @brief Instatiation of the gsFloater class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, D. Mokris
*/

#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsParametrization/gsFloater.h>
#include <gsModeling/gsParametrization/gsFloater.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsFloater<real_t>;

} // namespace gismo
