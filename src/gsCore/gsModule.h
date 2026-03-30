/** @file gsModule.h

    @brief Module registry utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsModules/gsModuleLoader.h>
#include <gsOptimizer/gsOptimizerRegistry.h>

/**
 * @defgroup gismo_modules Modules
 * @brief Optional modules that register themselves at runtime
 *
 * In MODULE mode (GISMO_BUILD_MODULE_LIB=ON), call gsModuleLoader::loadAll()
 * once at program start to discover and load all installed MODULE .so files.
 * In header-only or compiled-into-libgismo mode, include the module header
 * directly — the static registrar fires automatically.
 *
 * Example (MODULE mode):
 * @code
 * gismo::gsModuleLoader::loadAll();
 * if (gismo::gsOptimizerRegistry::has("gsOptim-BFGS"))
 *     auto opt = gismo::gsOptimizerRegistry::get("gsOptim-BFGS");
 * @endcode
 */
