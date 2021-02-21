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

namespace codi {
template<class T>
using codi_real = std::conditional<std::is_base_of<codi::Expression<typename codi::TypeTraits<T>::Real, T>, T>::value,
                                   GISMO_COEFF_TYPE, T>;

// RealReverseGen<real_t>
EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<codi_real<real_t>::type, codi_real<real_t>::type, LinearIndexHandler<int> >, ChunkVector> > >;

// RealForwardGen<real_t>
EXTERN_CLASS_TEMPLATE
ActiveReal<ForwardEvaluation<ForwardTapeTypes<real_t, codi_real<real_t>::type> > >;

// RealReverseIndexGen<real_t>
EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiIndexTape<JacobiIndexTapeTypes<ReverseTapeTypes<real_t, real_t, ReuseIndexHandlerUseCount<int> >, ChunkVector> > >;

// RealReverseGen<real_t, Direction<real_t, 4> >
EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiTape<JacobiTapeTypes<ReverseTapeTypes<real_t, Direction<real_t, 4>, LinearIndexHandler<int> >, ChunkVector > > >;

// RealReverseIndexGen<real_t, Direction<real_t, 4> >
EXTERN_CLASS_TEMPLATE
ActiveReal<JacobiIndexTape<JacobiIndexTapeTypes<ReverseTapeTypes<real_t, Direction<real_t, 4>, ReuseIndexHandlerUseCount<int> >, ChunkVector> > >;

} // namespace codi
#endif
