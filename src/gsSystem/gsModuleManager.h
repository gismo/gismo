/** @file gsModuleManager.h

    @brief Run-time manager of gismo modules: built-in ones and
    dynamic ones (.so/.dll).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsSystem/gsDyLib.h>

#include <string>
#include <vector>

namespace gismo
{

class gsRegistryBase;

/**
   @brief Registers the built-in modules and loads dynamic modules at run time.

   A dynamic module is a shared library named <tt>\<name\>.so</tt>
   (<tt>.dll</tt> on Windows) that exports two C-linkage symbols:
   - @c gismo_module_info()     : ABI handshake (see gsModuleInfo.h)
   - @c gismo_module_register() : adds the factories of the module to the
     gsRegistry objects

   The manager is a process-wide singleton, created on first use, which
   also registers the built-in modules. Loading never throws: a failure
   is reported with gsWarn and recorded in rejected().

   @code
   gsModuleManager & mm = gsModuleManager::get();
   mm.load("gsOpennurbs");            // by name, searched in searchPaths()
   mm.load("/path/to/gsOptim.so");    // by path
   mm.loadAll();                      // every module found in searchPaths()
   @endcode

   @ingroup System
*/
class GISMO_EXPORT gsModuleManager
{
public:
    /// The unique manager of the process. It is never destroyed, and
    /// registers the built-in modules when first called.
    static gsModuleManager & get();

    /// Loads the dynamic module @a nameOrPath. An argument that contains a
    /// path separator is a file path used as is; otherwise it is a module
    /// name, whose file is looked for in searchPaths(). Returns true if the
    /// module is present afterwards, including the case where it was
    /// present before; false, with a warning and an entry in rejected(),
    /// if the file is not found, cannot be opened, or fails the handshake
    /// (missing or empty module name, ABI version, scalar or index type).
    bool load(const std::string & nameOrPath);

    /// Loads every module file found in searchPaths() that is not already
    /// present nor already rejected. Scans the directories on every call.
    /// Returns the number of modules newly present.
    index_t loadAll();

    /// Removes the registrations of the dynamically loaded module @a name
    /// from every registry, then closes its library. Returns false, and
    /// does nothing, if @a name is not a dynamically loaded module, in
    /// particular for built-in modules. All objects made by the factories
    /// of the module must have been destroyed before this call.
    bool unload(const std::string & name);

    /// True if @a name is a built-in or a dynamically loaded module
    bool has(const std::string & name) const;

    /// True if @a name is a dynamically loaded module
    bool isLoaded(const std::string & name) const;

    /// Names of the built-in modules followed by those of the loaded modules
    std::vector<std::string> modules()  const;

    /// File paths, or names when no file was found, that failed to load
    std::vector<std::string> rejected() const;

    /// Directories searched for module files: the entries of the environment
    /// variable GISMO_MODULE_PATH (separated by ':', or ';' on Windows),
    /// then the directories given to addSearchPath(), then the build-time
    /// GISMO_MODULE_DIR, then GISMO_MODULE_INSTALL_DIR. Empty entries and
    /// repeated directories are skipped. The environment is read on every call.
    std::vector<std::string> searchPaths() const;

    /// Adds @a dir to the directories returned by searchPaths(); an empty
    /// or an already added @a dir is ignored
    void addSearchPath(const std::string & dir);

    /// Name of the module whose register function is running, otherwise empty
    const std::string & registering() const;

    /// Makes @a r known to the manager; called by the gsRegistryBase constructor
    void enlist(gsRegistryBase * r);

private:
    friend class gsRegistryBase;

    gsModuleManager() { }
    gsModuleManager(const gsModuleManager &);
    gsModuleManager & operator=(const gsModuleManager &);

    /// Forgets @a r; called by the gsRegistryBase destructor
    void delist(gsRegistryBase * r);

    /// A dynamically loaded module
    struct Record
    {
        std::string name;
        std::string path;
        gsDyLib     lib;
    };

    std::vector<std::string>     m_builtin;
    std::vector<Record>          m_records;
    std::vector<std::string>     m_rejected;
    std::vector<std::string>     m_searchPaths;
    std::vector<gsRegistryBase*> m_registries;
    std::string                  m_registering;
};

} // namespace gismo
