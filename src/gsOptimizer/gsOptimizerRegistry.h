/** @file gsOptimizerRegistry.h

    @brief The optimizer registry: gsRegistry<gsOptimizer<real_t>>.

    Optional modules (gsOptim, gsHLBFGS, ...) register their optimizers
    by name; users create them via
    gsOptimizerRegistry::get().make("gsOptim-BFGS").

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsRegistry.h>
#include <gsOptimizer/gsOptimizer.h>

namespace gismo
{

typedef gsRegistry<gsOptimizer<real_t> > gsOptimizerRegistry;

#ifdef GISMO_BUILD_LIB
EXTERN_CLASS_TEMPLATE gsRegistry<gsOptimizer<real_t> >;
#endif

} // namespace gismo
