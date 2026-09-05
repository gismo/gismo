/** @file gsModuleManager.cpp

    @brief Runtime manager for gismo module libraries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsCore/gsConfig.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsModules/gsModuleManager.h>
#include <gsModules/gsModuleInfo.h>

#include <cstring>
#include <algorithm>
#include <cstdlib>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#  define GS_MODULE_EXT ".dll"
#elif defined(__APPLE__)
#  include <dlfcn.h>
#  include <dirent.h>
#  define GS_MODULE_EXT ".dylib"
#else
#  include <dlfcn.h>
#  include <dirent.h>
#  define GS_MODULE_EXT ".so"
#endif

namespace gismo
{

// ---------------------------------------------------------------- gsDyLib

gsDyLib::gsDyLib(const std::string & path)
{
#if defined(_WIN32)
    m_handle = reinterpret_cast<void*>(::LoadLibraryA(path.c_str()));
    if (!m_handle)
    {
        char buf[256]; buf[0] = '\0';
        ::FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM|FORMAT_MESSAGE_IGNORE_INSERTS,
                         NULL, ::GetLastError(), 0, buf, sizeof(buf), NULL);
        m_error = buf;
    }
#else
    m_handle = ::dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!m_handle)
    {
        const char * e = ::dlerror();
        m_error = e ? e : "unknown error";
    }
#endif
}

gsDyLib::~gsDyLib()
{
    if (NULL != m_handle)
#if defined(_WIN32)
        ::FreeLibrary(reinterpret_cast<HMODULE>(m_handle));
#else
        ::dlclose(m_handle);
#endif
}

void * gsDyLib::symbol(const char * name) const
{
    if (NULL == m_handle) return NULL;
#if defined(_WIN32)
    return reinterpret_cast<void*>(
        ::GetProcAddress(reinterpret_cast<HMODULE>(m_handle), name));
#else
    return ::dlsym(m_handle, name);
#endif
}

// ------------------------------------------------------- directory listing

namespace {

std::vector<std::string> listModuleFiles(const std::string & dir)
{
    std::vector<std::string> result;
    const std::string ext(GS_MODULE_EXT);
#if defined(_WIN32)
    WIN32_FIND_DATAA fd;
    HANDLE h = ::FindFirstFileA((dir + "\\*" GS_MODULE_EXT).c_str(), &fd);
    if (INVALID_HANDLE_VALUE == h) return result;
    do
        if ( !(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) )
            result.push_back(dir + "\\" + fd.cFileName);
    while (::FindNextFileA(h, &fd));
    ::FindClose(h);
#else
    DIR * d = ::opendir(dir.c_str());
    if (!d) return result;
    while (struct dirent * entry = ::readdir(d))
    {
        const std::string name(entry->d_name);
        if (name.size() > ext.size() &&
            0 == name.compare(name.size()-ext.size(), ext.size(), ext) &&
            0 != name.find("libgismo"))
            result.push_back(dir + "/" + name);
    }
    ::closedir(d);
#endif
    return result;
}

} // namespace

// -------------------------------------------------------- gsModuleManager

gsModuleManager & gsModuleManagerSingleton()
{
    // Never destroyed: registered factories must stay valid through
    // static destruction of other objects.
    static gsModuleManager * mm = new gsModuleManager;
    return *mm;
}

bool gsModuleManager::load(const std::string & path)
{
    gsDyLib lib(path);
    if ( !lib.isOpen() )
    {
        gsWarn << "gsModuleManager: cannot open '" << path << "': "
               << lib.error() << "\n";
        m_rejected.push_back(path);
        return false;
    }

    typedef const gismo_module_info_t * (*InfoFn)();
    typedef void (*RegisterFn)();

    InfoFn infoFn = reinterpret_cast<InfoFn>(lib.symbol("gismo_module_info"));
    const gismo_module_info_t * info = infoFn ? infoFn() : NULL;
    if ( NULL == info )
    {
        gsWarn << "gsModuleManager: '" << path << "' is not a gismo module "
               << "(no usable gismo_module_info) - rejected\n";
        m_rejected.push_back(path);
        return false;
    }

    // ABI handshake: major exact, minor not newer than host, scalar
    // and index types exact
    if ( info->abi_major != static_cast<uint32_t>(GISMO_MAJOR) ||
         info->abi_minor >  static_cast<uint32_t>(GISMO_MINOR) )
    {
        gsWarn << "gsModuleManager: '" << info->module_name << "' ABI "
               << info->abi_major << "." << info->abi_minor
               << " does not match host " << GISMO_MAJOR << "." << GISMO_MINOR
               << " - rejected\n";
        m_rejected.push_back(path);
        return false;
    }
    if ( 0 != std::strcmp(info->coeff_type, GISMO_COEFF_TYPE_STR) ||
         0 != std::strcmp(info->index_type, GISMO_INDEX_TYPE_STR) )
    {
        gsWarn << "gsModuleManager: '" << info->module_name << "' built for ("
               << info->coeff_type << ", " << info->index_type
               << "), host is (" << GISMO_COEFF_TYPE_STR << ", "
               << GISMO_INDEX_TYPE_STR << ") - rejected\n";
        m_rejected.push_back(path);
        return false;
    }

    if ( RegisterFn regFn =
             reinterpret_cast<RegisterFn>(lib.symbol("gismo_module_register")) )
        regFn();

    // The library must stay mapped: registered factories point into it
    m_libs.push_back(give(lib));
    m_loaded.push_back(path);
    m_names.push_back(info->module_name);
    gsInfo << "gsModuleManager: loaded '" << info->module_name
           << "' v" << info->module_version << " (" << path << ")\n";
    return true;
}

bool gsModuleManager::isLoaded(const std::string & name) const
{
    return std::find(m_names.begin(), m_names.end(), name) != m_names.end();
}

index_t gsModuleManager::loadAll()
{
    const char * env = std::getenv("GISMO_MODULE_DIR");
    const std::string dir = env ? std::string(env)
                                : std::string(GISMO_MODULE_INSTALL_DIR);
    index_t count = 0;
    const std::vector<std::string> files = listModuleFiles(dir);
    for (size_t i = 0; i != files.size(); ++i)
        if ( load(files[i]) ) ++count;
    return count;
}

} // namespace gismo
