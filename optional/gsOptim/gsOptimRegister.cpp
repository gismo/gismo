/** @file gsOptimRegister.cpp

    @brief Registers the gsOptim optimizers in the optimizer registry.

    Compiled into libgismo when gsOptim is a compiled-in optional module,
    and into gsOptim_module when built as a runtime module
    (GISMO_BUILD_MODULE_LIB=ON); registration is idempotent either way.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsOptim/gsOptim.h>
#include <gsOptimizer/gsOptimizerRegistry.h>

namespace gismo
{

GISMO_EXPORT void gsOptimRegisterAll()
{
    typedef gsOptimizer<real_t> B;
    gsOptimizerRegistry & reg = gsOptimizerRegistry::get();
    reg.add("gsOptim-BFGS"  , gsRegistryFactory<B, gsOptimBFGS<real_t>   >);
    reg.add("gsOptim-LBFGS" , gsRegistryFactory<B, gsOptimLBFGS<real_t>  >);
    reg.add("gsOptim-CG"    , gsRegistryFactory<B, gsOptimCG<real_t>     >);
    reg.add("gsOptim-GD"    , gsRegistryFactory<B, gsOptimGD<real_t>     >);
    reg.add("gsOptim-NM"    , gsRegistryFactory<B, gsOptimNM<real_t>     >);
    reg.add("gsOptim-DE"    , gsRegistryFactory<B, gsOptimDE<real_t>     >);
    reg.add("gsOptim-DEPRMM", gsRegistryFactory<B, gsOptimDEPRMM<real_t> >);
    reg.add("gsOptim-PSO"   , gsRegistryFactory<B, gsOptimPSO<real_t>    >);
    reg.add("gsOptim-PSODV" , gsRegistryFactory<B, gsOptimPSODV<real_t>  >);
    reg.add("gsOptim-SUMT"  , gsRegistryFactory<B, gsOptimSUMT<real_t>   >);
}

namespace {
// Fires on library load (libgismo or gsOptim_module)
struct gsOptimRegistrar { gsOptimRegistrar() { gsOptimRegisterAll(); } };
static const gsOptimRegistrar s_gsOptimRegistrar;
}

} // namespace gismo

// Anchor against linker dead-stripping when gismo is consumed as a
// static archive (referenced via force-undefine, see the modularization
// report S12/S3-A1 anchor mechanism)
extern "C" GISMO_EXPORT void gismo_anchor_gsOptim(void) { }
