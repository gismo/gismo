/** @file gsRegistry_.cpp

    @brief Explicit instantiation of gsRegistry for the types that
    libgismo exports a registry of.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
*/

#define gsRegistry_EXPORT

#include <gsCore/gsTemplateTools.h>

#include <gsSystem/gsFileIo.h>
#include <gsSystem/gsRegistry.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsRegistry<gsFileIo>;

}
