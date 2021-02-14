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

#ifdef GISMO_BUILD_LIB
#ifdef gsCoDiPack_EXPORTS
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif

namespace codi
{

#ifndef real_t
EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<real_t, real_t, LinearIndexHandler<int> >, ChunkVector> > >;

EXTERN_CLASS_TEMPLATE
ActiveReal<ForwardEvaluation<ForwardTapeTypes<real_t, real_t> > >;
#else
EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<GISMO_COEFF_TYPE, GISMO_COEFF_TYPE, LinearIndexHandler<int> >, ChunkVector> > >;

EXTERN_CLASS_TEMPLATE
ActiveReal<ForwardEvaluation<ForwardTapeTypes<GISMO_COEFF_TYPE, GISMO_COEFF_TYPE> > >;
#endif
}
#endif
