/** @file gsModuleManager.cpp

    @brief Run-time manager of gismo modules, and the out-of-line part of
    gsRegistryBase.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsSystem/gsModuleManager.h>
#include <gsSystem/gsRegistry.h>
#include <gsSystem/gsModule.h>
#include <gsSystem/gsModuleInfo.h>
#include <gsSystem/gsModuleConfig.h>
#include <gsIO/gsFileManager.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#else
#  include <dirent.h>
#endif

namespace gismo
{

namespace {

#if defined(_WIN32)
const char * const moduleSuffix = ".dll";
const char          pathListSeparator = ';';
#else
const char * const moduleSuffix = ".so";
const char          pathListSeparator = ':';
#endif

/// Full paths of the module files (suffix @c moduleSuffix, name not
/// starting with "libgismo") directly inside @a dir
std::vector<std::string> listModuleFiles(const std::string & dir)
{
    std::vector<std::string> result;
    const std::string ext(moduleSuffix);
    const char sep = gsFileManager::getNativePathSeparator();
#if defined(_WIN32)
    WIN32_FIND_DATAA fd;
    HANDLE h = ::FindFirstFileA((dir + sep + "*" + ext).c_str(), &fd);
    if (INVALID_HANDLE_VALUE == h) return result;
    do
    {
        const std::string name(fd.cFileName);
        if ( !(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) &&
             0 != name.find("libgismo") )
            result.push_back(dir + sep + name);
    }
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
            result.push_back(dir + sep + name);
    }
    ::closedir(d);
#endif
    return result;
}

/// File name of @a path without directory and without the module suffix
std::string moduleNameOf(const std::string & path)
{
    const size_t pos = path.find_last_of(gsFileManager::getValidPathSeparators());
    std::string name = (std::string::npos == pos) ? path : path.substr(pos+1);
    const std::string ext(moduleSuffix);
    if (name.size() >= ext.size() &&
        0 == name.compare(name.size()-ext.size(), ext.size(), ext))
        name.erase(name.size()-ext.size());
    return name;
}

/// Appends @a s to @a v unless it is empty or already there
void appendUnique(std::vector<std::string> & v, const std::string & s)
{
    if ( !s.empty() && std::find(v.begin(), v.end(), s) == v.end() )
        v.push_back(s);
}

} // namespace

// ------------------------------------------------------------ gsRegistryBase

gsRegistryBase::gsRegistryBase()
{
    gsModuleManager::get().enlist(this);
}

gsRegistryBase::~gsRegistryBase()
{
    gsModuleManager::get().delist(this);
}

// -------------------------------------------------------------- gsModuleManager

gsModuleManager & gsModuleManager::get()
{
    // Never destroyed: registered factories must stay valid through
    // static destruction of other objects. Out of line, so that every
    // module reaches this one instance through libgismo.
    static gsModuleManager * mm = NULL;
    if (NULL == mm)
    {
        // Assigned before the built-in modules register: their registries
        // call get() again, from their base constructor.
        mm = new gsModuleManager;
        for (const char * const * n = gsBuiltinModuleNames(); NULL != *n; ++n)
        {
            mm->m_builtin.push_back(*n);
            mm->m_registering = *n;
            gsRegisterBuiltinModule(*n);
            mm->m_registering.clear();
        }
    }
    return *mm;
}

void gsModuleManager::enlist(gsRegistryBase * r)
{
    if ( std::find(m_registries.begin(), m_registries.end(), r) == m_registries.end() )
        m_registries.push_back(r);
}

void gsModuleManager::delist(gsRegistryBase * r)
{
    m_registries.erase(std::remove(m_registries.begin(), m_registries.end(), r),
                       m_registries.end());
}

const std::string & gsModuleManager::registering() const
{ return m_registering; }

bool gsModuleManager::isLoaded(const std::string & name) const
{
    for (size_t i = 0; i != m_records.size(); ++i)
        if (m_records[i].name == name) return true;
    return false;
}

bool gsModuleManager::has(const std::string & name) const
{
    return std::find(m_builtin.begin(), m_builtin.end(), name) != m_builtin.end()
        || isLoaded(name);
}

std::vector<std::string> gsModuleManager::modules() const
{
    std::vector<std::string> res(m_builtin);
    for (size_t i = 0; i != m_records.size(); ++i)
        res.push_back(m_records[i].name);
    return res;
}

std::vector<std::string> gsModuleManager::rejected() const
{ return m_rejected; }

std::vector<std::string> gsModuleManager::searchPaths() const
{
    std::vector<std::string> res;

    if ( const char * env = std::getenv("GISMO_MODULE_PATH") )
    {
        const std::string list(env);
        size_t start = 0;
        while (start <= list.size())
        {
            size_t end = list.find(pathListSeparator, start);
            if (std::string::npos == end) end = list.size();
            appendUnique(res, list.substr(start, end-start));
            start = end + 1;
        }
    }

    for (size_t i = 0; i != m_searchPaths.size(); ++i)
        appendUnique(res, m_searchPaths[i]);

    appendUnique(res, GISMO_MODULE_DIR);
    appendUnique(res, GISMO_MODULE_INSTALL_DIR);
    return res;
}

void gsModuleManager::addSearchPath(const std::string & dir)
{ appendUnique(m_searchPaths, dir); }

bool gsModuleManager::load(const std::string & nameOrPath)
{
    std::string path;
    if ( std::string::npos != nameOrPath.find_first_of(gsFileManager::getValidPathSeparators()) )
        path = nameOrPath;
    else
    {
        if ( has(nameOrPath) ) return true;

        const std::vector<std::string> dirs = searchPaths();
        for (size_t i = 0; i != dirs.size() && path.empty(); ++i)
        {
            const std::string candidate = dirs[i]
                + gsFileManager::getNativePathSeparator() + nameOrPath + moduleSuffix;
            if ( gsFileManager::fileExists(candidate) )
                path = candidate;
        }
        if ( path.empty() )
        {
            gsWarn << "gsModuleManager: no file for module '" << nameOrPath
                   << "' found in the search paths\n";
            appendUnique(m_rejected, nameOrPath);
            return false;
        }
    }

    for (size_t i = 0; i != m_records.size(); ++i)
        if (m_records[i].path == path) return true;

    gsDyLib lib(path);
    if ( !lib.isOpen() )
    {
        gsWarn << "gsModuleManager: cannot open '" << path << "': "
               << lib.error() << "\n";
        appendUnique(m_rejected, path);
        return false;
    }

    typedef const gismo_module_info_t * (*InfoFn)();
    typedef void (*RegisterFn)();

    InfoFn infoFn = reinterpret_cast<InfoFn>(lib.symbol("gismo_module_info"));
    const gismo_module_info_t * info = infoFn ? infoFn() : NULL;
    if ( NULL == info || NULL == info->module_name || '\0' == info->module_name[0] )
    {
        gsWarn << "gsModuleManager: '" << path << "' is not a gismo module "
               << "(no usable gismo_module_info) - rejected\n";
        appendUnique(m_rejected, path);
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
        appendUnique(m_rejected, path);
        return false;
    }
    if ( NULL == info->coeff_type || NULL == info->index_type ||
         0 != std::strcmp(info->coeff_type, GISMO_COEFF_TYPE_STR) ||
         0 != std::strcmp(info->index_type, GISMO_INDEX_TYPE_STR) )
    {
        gsWarn << "gsModuleManager: '" << info->module_name << "' built for ("
               << (info->coeff_type ? info->coeff_type : "?") << ", "
               << (info->index_type ? info->index_type : "?")
               << "), host is (" << GISMO_COEFF_TYPE_STR << ", "
               << GISMO_INDEX_TYPE_STR << ") - rejected\n";
        appendUnique(m_rejected, path);
        return false;
    }

    // The same module from another file, or a built-in one: keep the first
    const std::string name(info->module_name);
    if ( has(name) ) return true;

    if ( RegisterFn regFn =
             reinterpret_cast<RegisterFn>(lib.symbol("gismo_module_register")) )
    {
        m_registering = name;
        regFn();
        m_registering.clear();
    }

    // The library must stay mapped: registered factories point into it
    Record rec;
    rec.name = name;
    rec.path = path;
    rec.lib  = give(lib);
    m_records.push_back(give(rec));
    return true;
}

index_t gsModuleManager::loadAll()
{
    const size_t before = m_records.size();
    const std::vector<std::string> dirs = searchPaths();
    for (size_t i = 0; i != dirs.size(); ++i)
    {
        const std::vector<std::string> files = listModuleFiles(dirs[i]);
        for (size_t j = 0; j != files.size(); ++j)
        {
            if ( has(moduleNameOf(files[j])) ||
                 std::find(m_rejected.begin(), m_rejected.end(), files[j]) != m_rejected.end() )
                continue;
            load(files[j]);
        }
    }
    return static_cast<index_t>(m_records.size() - before);
}

bool gsModuleManager::unload(const std::string & name)
{
    for (size_t i = 0; i != m_records.size(); ++i)
    {
        if (m_records[i].name != name) continue;

        for (size_t j = 0; j != m_registries.size(); ++j)
            m_registries[j]->removeModule(name);
        m_records.erase(m_records.begin() + i);
        return true;
    }
    return false;
}

} // namespace gismo
