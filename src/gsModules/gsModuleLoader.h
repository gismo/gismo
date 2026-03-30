/** @file gsModuleLoader.h

    @brief Runtime loader for gismo MODULE (.so/.dll) libraries

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsExport.h>
#include <string>
#include <vector>

namespace gismo
{

/**
 * @brief Loads gismo MODULE shared libraries at runtime.
 *
 * Each MODULE .so exports two C-linkage symbols:
 *  - @c gismo_module_info()   – returns a pointer to a @c gismo_module_info_t
 *                               struct used for the version handshake.
 *  - @c gismo_module_register() – triggers factory registration in the
 *                               relevant registry (e.g. gsOptimizerRegistry).
 *
 * ### Typical usage
 * @code
 * // At program start, before any registry query:
 * gismo::gsModuleLoader::loadAll();
 *
 * // Or load a specific module:
 * gismo::gsModuleLoader::load("/path/to/gsOptim_module.so");
 *
 * // Now use the registry as normal:
 * if (gismo::gsOptimizerRegistry::has("gsOptim-BFGS"))
 *     auto opt = gismo::gsOptimizerRegistry::get("gsOptim-BFGS");
 * @endcode
 *
 * ### Module search order for loadAll()
 *  1. @c GISMO_MODULE_DIR environment variable (if set).
 *  2. Compile-time install prefix: @c GISMO_MODULE_INSTALL_DIR (from gsConfig.h).
 *
 * @ingroup gismo_modules
 */
class GISMO_EXPORT gsModuleLoader
{
public:
    /**
     * @brief Scan the module directory and load all MODULE .so files found.
     *
     * Silently skips files that are not valid gismo modules.
     * Emits a warning for modules that fail the version handshake.
     */
    static void loadAll();

    /**
     * @brief Load a single MODULE .so by full filesystem path.
     *
     * Performs the version handshake (abi_major, coeff_type) and calls
     * gismo_module_register() on success.
     *
     * @param path  Full path to the .so / .dll file.
     * @return true if the module was loaded successfully.
     */
    static bool load(const std::string& path);

    /// Return paths of all successfully loaded modules.
    static std::vector<std::string> loadedModules();

    /// Return paths of modules that failed the version handshake.
    static std::vector<std::string> rejectedModules();

private:
    gsModuleLoader() = delete; // non-constructible utility class
};

} // namespace gismo
