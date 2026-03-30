/** @file gsOptimizerRegistry.cpp

    @brief Registry singleton for optional optimizers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsOptimizer/gsOptimizerRegistry.h>

namespace gismo
{

gsOptimizerRegistry::RegistryData& gsOptimizerRegistry::getMutableInstance()
{
    // Heap-allocated to avoid static-destruction order issues.
    // One authoritative definition lives in libgismo.so; MODULE .so files
    // resolve this symbol from the host library (RTLD_LOCAL keeps their own
    // symbols local, but calls to libgismo exports go to the shared copy).
    static RegistryData* instance = new RegistryData();
    return *instance;
}

} // namespace gismo
