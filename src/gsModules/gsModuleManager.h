/** @file gsModuleManager.h

    @brief Runtime manager for gismo module libraries (.so/.dll).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsModules/gsDyLib.h>
#include <string>
#include <vector>

namespace gismo
{

class gsModuleManager;

/// Singleton function returning the module manager
GISMO_EXPORT gsModuleManager & gsModuleManagerSingleton();

/**
   @brief Loads gismo module libraries at runtime.

   A module library exports two C-linkage symbols:
   - \c gismo_module_info()     : ABI handshake (see gsModuleInfo.h)
   - \c gismo_module_register() : populates the relevant gsRegistry

   Use as:
   \code
   gsModuleManager & mm = gsModuleManager::get();
   mm.loadAll();                       // scan GISMO_MODULE_DIR / install dir
   mm.load("/path/to/gsOptim_module.so");
   \endcode

   A rejected or failed load is reported with gsWarn and recorded in
   rejected(); it never throws. Loaded libraries stay mapped for the
   process lifetime (registered factories must remain callable).

   \ingroup Modules
*/
class GISMO_EXPORT gsModuleManager
{
    friend GISMO_EXPORT gsModuleManager & gsModuleManagerSingleton();
public:

    static gsModuleManager & get() { return gsModuleManagerSingleton(); }

    /// Loads all module libraries found in the module directory
    /// (the \c GISMO_MODULE_DIR environment variable if set, otherwise
    /// the compiled-in install location). Returns the number loaded.
    index_t loadAll();

    /// Loads a single module library by path; true on success
    bool load(const std::string & path);

    /// Paths of successfully loaded modules
    const std::vector<std::string> & loaded()   const { return m_loaded;   }

    /// Names (gismo_module_info_t::module_name) of the loaded modules
    const std::vector<std::string> & loadedNames() const { return m_names; }

    /// True if a runtime module named \a name has been loaded
    bool isLoaded(const std::string & name) const;

    /// Paths of modules that failed the handshake or did not load
    const std::vector<std::string> & rejected() const { return m_rejected; }

private:
    gsModuleManager() { }
    gsModuleManager(const gsModuleManager &);

    std::vector<gsDyLib>     m_libs;
    std::vector<std::string> m_loaded, m_rejected, m_names;
};

} // namespace gismo
