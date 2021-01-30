/** @file gsCoDiPack.h

    @brief Header for CoDiPack package

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Moeller
*/

#pragma once

#include <codi.hpp>
#include <gsCore/gsTemplateTools.h>

namespace codi
{

#if defined(GISMO_BUILD_LIB) && defined(gsCoDiPack_EXPORTS)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif

EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<double, double, LinearIndexHandler<int> >, ChunkVector> > >;

EXTERN_CLASS_TEMPLATE
ActiveReal<ForwardEvaluation<ForwardTapeTypes<double, double> > >;

}
