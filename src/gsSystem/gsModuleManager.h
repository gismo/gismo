/** @file gsModuleManager.h

    @brief Provides definition of the ModuleManager class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{
/// Singleton function returning the gsMpi helper object
GISMO_EXPORT gsMpi & gsMpiSingleton(const int& argc, char** argv);

class gsDyLib
{
    void * handle = nullptr;

public:
    gsDyLib() : handle(nullptr) { }

    gsDyLib(const std::string& name)
    {
        std::string library_name = name;
        #ifdef _WIN32
		//no op
#elif __APPLE__
		library_name = std::string("./lib" + library_name + ".dylib");
#else
		//On unix the library name probably starts with `lib`
		library_name = std::string("./lib" + library_name + ".so");
#endif
        
        handle = 
#ifdef _WIN32
        reinterpret_cast<void*>(LoadLibraryA(library_name.c_str()));
#else
        dlopen(library_name.c_str(), RTLD_NOW);
#endif
        if(!handle) GISMO_ERROR("Module library did not load");
    }

    ~gsDyLib()
    {
        if(handle != nullptr)
        {
#ifdef _WIN32
            FreeLibrary(HMODULE(lib_handle));
#else
            dlclose(lib_handle);
#endif
            handle = nullptr;
        }
    }

    //deactivate copy
    gsDyLib& operator=(const gsDyLib&) = delete;
    gsDyLib(const gsDyLib&)			  = delete;

    //move assignment
    gsDyLib& operator=(gsDyLib&& other) noexcept
    {
        handle		 = other.handle;
        other.handle = nullptr;
        return *this;
    }

    gsDyLib(gsDyLib&& other) noexcept
    {
        *this = std::move(other);
    }

    plugin* getPluginInstance()
    {
        using boostrap_function = plugin* (*)();

        std::string& function_name =  ___YPLUGS_BOOTSTRAP_PROC_NAME_STR;
        auto fptr = reinterpret_cast<boostrap_function>(
#ifdef _WIN32
			GetProcAddress(HMODULE(handle), LPCSTR(function_name.c_str()));
#else
            dlsym(handle, function_name.c_str())
#endif
            );
    
        if(!fptr)
            throw std::runtime_error(
                "Did not find yplugs bootstrap function in plugin library");
        return fptr();
    }
};

class plugin
{
public:
    virtual ~plugin() { }
    virtual bool init() = 0;
    virtual bool quit() = 0;
};

/** 
    @brief 
    \ingroup System
*/
class gsModuleManager
{
    friend GISMO_EXPORT gsMpi & gsMpiSingleton(const int& argc, char** argv);
public:
    plugin* load_plugin(const std::string& library_name)
    {
        auto library_iterator = libraries.find(library_name);
        if(library_iterator == libraries.end())
        {
            libraries[library_name] = gsDyLib(library_name);
            library_iterator		= libraries.find(library_name);
            assert(library_iterator != libraries.end());
        }
        gsDyLib& library = library_iterator->second;

        plugin* plugin_instance_ptr = library.getPluginInstance();
        plugins.emplace_back(plugin_instance_ptr);

        assert((uintptr_t)plugin_instance_ptr
               == (uintptr_t)plugins.back().get());
        plugin_instance_ptr->init();

        return plugin_instance_ptr;
    }

    ~gsModuleManager()
    {
        for(auto& plugin : plugins) plugin->quit();
    }

private:
    gsModuleManager();
    gsModuleManager();
    std::unordered_map<std::string, gsDyLib> libraries;
    std::vector<std::unique_ptr<plugin>> plugins;
};

} // namespace gismo
