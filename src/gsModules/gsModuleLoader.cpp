/** @file gsModuleLoader.cpp

    @brief Runtime loader for gismo MODULE (.so/.dll) libraries

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsModules/gsModuleLoader.h>
#include <gsModules/gsModuleInfo.h>
#include <gsCore/gsConfig.h>   // GISMO_MAJOR, GISMO_MINOR, GISMO_COEFF_TYPE_STR
#include <gsCore/gsDebug.h>    // gsWarn, gsInfo

#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>             // std::getenv

// -----------------------------------------------------------------------
// Platform-specific dynamic loading
// -----------------------------------------------------------------------
#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
   using gs_module_handle_t = HMODULE;
#  define gs_dlopen(p)     ::LoadLibraryA(p)
#  define gs_dlsym(h,s)    reinterpret_cast<void*>(::GetProcAddress(h,s))
#  define gs_dlclose(h)    ::FreeLibrary(h)
   static std::string gs_dlerror()
   {
       char buf[256]; buf[0] = '\0';
       ::FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM|FORMAT_MESSAGE_IGNORE_INSERTS,
           NULL, ::GetLastError(), 0, buf, sizeof(buf), NULL);
       return std::string(buf);
   }
#  define GS_MODULE_EXT    ".dll"
#else
#  include <dlfcn.h>
   using gs_module_handle_t = void*;
#  define gs_dlopen(p)     ::dlopen(p, RTLD_NOW | RTLD_LOCAL)
#  define gs_dlsym(h,s)    ::dlsym(h,s)
#  define gs_dlclose(h)    ::dlclose(h)
#  define gs_dlerror()     std::string(::dlerror() ? ::dlerror() : "unknown error")
#  if defined(__APPLE__)
#    define GS_MODULE_EXT  ".dylib"
#  else
#    define GS_MODULE_EXT  ".so"
#  endif
#endif

// -----------------------------------------------------------------------
// Platform-specific directory iteration
// -----------------------------------------------------------------------
#if defined(_WIN32)
#  include <windows.h>
   static std::vector<std::string> gs_list_so_files(const std::string& dir)
   {
       std::vector<std::string> result;
       WIN32_FIND_DATAA fd;
       std::string pattern = dir + "\\*" GS_MODULE_EXT;
       HANDLE h = ::FindFirstFileA(pattern.c_str(), &fd);
       if (h == INVALID_HANDLE_VALUE) return result;
       do {
           if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
               result.push_back(dir + "\\" + fd.cFileName);
       } while (::FindNextFileA(h, &fd));
       ::FindClose(h);
       return result;
   }
#else
#  include <dirent.h>
#  include <sys/stat.h>
    static std::vector<std::string> gs_list_so_files(const std::string& dir)
    {
        std::vector<std::string> result;
        DIR* d = ::opendir(dir.c_str());
        if (!d) return result;
        struct dirent* entry;
        const std::string ext(GS_MODULE_EXT);
        while ((entry = ::readdir(d)) != NULL)
        {
            std::string name(entry->d_name);
            if (name.size() > ext.size() &&
                name.compare(name.size() - ext.size(), ext.size(), ext) == 0 &&
                name.find("libgismo") != 0)
            {
                result.push_back(dir + "/" + name);
            }
        }
        ::closedir(d);
        return result;
    }
#endif

namespace gismo
{

// -----------------------------------------------------------------------
// Internal singleton state
// -----------------------------------------------------------------------
namespace {

struct ModuleLoaderData
{
    std::vector<std::string> loaded;
    std::vector<std::string> rejected;
};

static ModuleLoaderData& loaderData()
{
    static ModuleLoaderData* instance = new ModuleLoaderData();
    return *instance;
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Function pointer typedefs matching the MODULE entry point signatures
// -----------------------------------------------------------------------
using InfoFn     = const gismo_module_info_t*(*)();
using RegisterFn = void(*)();

// -----------------------------------------------------------------------
// gsModuleLoader::load
// -----------------------------------------------------------------------
bool gsModuleLoader::load(const std::string& path)
{
    gs_module_handle_t handle = gs_dlopen(path.c_str());
    if (!handle)
    {
        gsWarn << "gsModuleLoader: cannot open '" << path
               << "': " << gs_dlerror() << "\n";
        return false;
    }

    // --- Step 1: version handshake ---
    void* sym_info = gs_dlsym(handle, "gismo_module_info");
    if (!sym_info)
    {
        gsWarn << "gsModuleLoader: '" << path
               << "' has no 'gismo_module_info' symbol — not a gismo module\n";
        gs_dlclose(handle);
        loaderData().rejected.push_back(path);
        return false;
    }

    InfoFn info_fn = reinterpret_cast<InfoFn>(sym_info);
    const gismo_module_info_t* info = info_fn();

    // ABI major must match exactly
    if (info->abi_major != static_cast<uint32_t>(GISMO_MAJOR))
    {
        gsWarn << "gsModuleLoader: '" << path << "' ABI major mismatch: "
               << "module=" << info->abi_major
               << " host=" << GISMO_MAJOR << " — rejected\n";
        gs_dlclose(handle);
        loaderData().rejected.push_back(path);
        return false;
    }

    // ABI minor: module must not be newer than host
    if (info->abi_minor > static_cast<uint32_t>(GISMO_MINOR))
    {
        gsWarn << "gsModuleLoader: '" << path << "' ABI minor too new: "
               << "module=" << info->abi_minor
               << " host=" << GISMO_MINOR << " — rejected\n";
        gs_dlclose(handle);
        loaderData().rejected.push_back(path);
        return false;
    }

    // Coefficient type must match exactly
    if (std::strcmp(info->coeff_type, GISMO_COEFF_TYPE_STR) != 0)
    {
        gsWarn << "gsModuleLoader: '" << path << "' coeff_type mismatch: "
               << "module='" << info->coeff_type
               << "' host='" << GISMO_COEFF_TYPE_STR << "' — rejected\n";
        gs_dlclose(handle);
        loaderData().rejected.push_back(path);
        return false;
    }

    // --- Step 2: register ---
    void* sym_reg = gs_dlsym(handle, "gismo_module_register");
    if (sym_reg)
    {
        RegisterFn reg_fn = reinterpret_cast<RegisterFn>(sym_reg);
        reg_fn();
    }

    // Note: handle is intentionally NOT closed — the .so must remain mapped
    // so that factory lambdas stored in the registry remain valid.
    loaderData().loaded.push_back(path);
    gsInfo << "gsModuleLoader: loaded module '" << info->module_name
           << "' v" << info->module_version
           << " from '" << path << "'\n";
    return true;
}

// -----------------------------------------------------------------------
// gsModuleLoader::loadAll
// -----------------------------------------------------------------------
void gsModuleLoader::loadAll()
{
    // Search path priority:
    //   1. GISMO_MODULE_DIR environment variable
    //   2. Compile-time install prefix (GISMO_MODULE_INSTALL_DIR from gsConfig.h)
    const char* env_dir = std::getenv("GISMO_MODULE_DIR");
    std::string dir = env_dir ? std::string(env_dir) : std::string(GISMO_MODULE_INSTALL_DIR);

    std::vector<std::string> files = gs_list_so_files(dir);
    for (const std::string& f : files)
        load(f);
}

// -----------------------------------------------------------------------
// Accessors
// -----------------------------------------------------------------------
std::vector<std::string> gsModuleLoader::loadedModules()
{
    return loaderData().loaded;
}

std::vector<std::string> gsModuleLoader::rejectedModules()
{
    return loaderData().rejected;
}

} // namespace gismo
