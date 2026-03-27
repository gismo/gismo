/** @file gsModuleManager.h

    @brief The module manager for dynamically loaded modules for G+Smo.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsSystem/gsModule.h>
#include <gsSystem/gsDyLib.h>

#define ___YPLUGS_BOOTSTRAP_PROC_NAME __y_plugs_gsModule_entry_point
#define ___YPLUGS_BOOTSTRAP_XSTR(s) ___YPLUGS_BOOTSTRAP_STR(s)
#define ___YPLUGS_BOOTSTRAP_STR(s) #s
#define ___YPLUGS_BOOTSTRAP_PROC_NAME_STR                                      \
	___YPLUGS_BOOTSTRAP_STR(__y_plugs_gsModule_entry_point)

#ifdef _WIN32
#define ___YPLUGS_BOOTSTRAP_EXPORT_SYMBOL _declspec(dllexport)
#else
#define ___YPLUGS_BOOTSTRAP_EXPORT_SYMBOL
#endif
#define YPLUGS_BOOTSRAP_SIGNATURE                                              \
	extern "C" ___YPLUGS_BOOTSTRAP_EXPORT_SYMBOL yba::gsModule*               \
		___YPLUGS_BOOTSTRAP_PROC_NAME

#include <string>

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include "Windows.h"
#else //UNIX
#include <dlfcn.h>
#endif

namespace gismo
{
/// Singleton function returning the gsMpi helper object
GISMO_EXPORT gsMpi & gsMpiSingleton(const int& argc, char** argv);

/** 
    @brief 
    \ingroup System
*/
class gsModuleManager
{
    friend GISMO_EXPORT gsMpi & gsMpiSingleton(const int& argc, char** argv);
public:
    gsModule* load_plugin(const std::string& library_name)
    {
        auto library_iterator = m_libraries.find(library_name);
        if(library_iterator == m_libraries.end())
        {
            m_libraries[library_name] = gsDyLib(library_name);
            library_iterator		= m_libraries.find(library_name);
            assert(library_iterator != m_libraries.end());
        }
        gsDyLib& library = library_iterator->second;

        gsModule* plugin_instance_ptr = library.getPluginInstance();
        m_modules.emplace_back(plugin_instance_ptr);

        assert((uintptr_t)plugin_instance_ptr
               == (uintptr_t)m_modules.back().get());
        plugin_instance_ptr->init();

        return plugin_instance_ptr;
    }

    ~gsModuleManager()
    {
        for(auto& plugin : m_modules)
            plugin->unload();
    }

private:
    gsModuleManager();
    gsModuleManager();
    std::unordered_map<std::string, gsDyLib> m_libraries;
    std::vector<std::unique_ptr<gsModule>> m_modules;
};

} // namespace gismo
