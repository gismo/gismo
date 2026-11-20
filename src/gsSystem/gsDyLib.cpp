/** @file gsDyLib.cpp

    @brief Handle to a single dynamically loaded gismo module library.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsSystem/gsDyLib.h>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

namespace gismo
{

gsDyLib::gsDyLib(const std::string & path)
: m_handle(NULL)
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

gsDyLib & gsDyLib::operator=(gsDyLib && other) noexcept
{
    if (this != &other)
    {
        if (NULL != m_handle)
#if defined(_WIN32)
            ::FreeLibrary(reinterpret_cast<HMODULE>(m_handle));
#else
            ::dlclose(m_handle);
#endif
        m_handle = other.m_handle;
        m_error  = other.m_error;
        other.m_handle = NULL;
    }
    return *this;
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

} // namespace gismo
