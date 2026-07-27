/** @file gsDyLib.h

    @brief RAII wrapper for a dynamically loaded library handle.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsExport.h>
#include <string>

namespace gismo
{

/**
   @brief Move-only RAII owner of a dynamic library handle
   (dlopen / LoadLibrary).

   \ingroup Modules
*/
class GISMO_EXPORT gsDyLib
{
public:
    gsDyLib() : m_handle(NULL) { }

    /// Opens the library at \a path (full file name). On failure the
    /// handle stays NULL and error() returns the reason.
    explicit gsDyLib(const std::string & path);

    ~gsDyLib();

    gsDyLib(const gsDyLib &) = delete;
    gsDyLib & operator=(const gsDyLib &) = delete;

    gsDyLib(gsDyLib && other) : m_handle(other.m_handle), m_error(other.m_error)
    { other.m_handle = NULL; }

    gsDyLib & operator=(gsDyLib && other)
    {
        std::swap(m_handle, other.m_handle);
        std::swap(m_error , other.m_error );
        return *this;
    }

    /// True if a library is open
    bool isOpen() const { return NULL != m_handle; }

    /// Looks up \a symbol; returns NULL if absent
    void * symbol(const char * name) const;

    /// Releases ownership of the handle without closing it (used when
    /// the library must stay mapped for the process lifetime)
    void release() { m_handle = NULL; }

    /// The last error message of open/symbol lookup
    const std::string & error() const { return m_error; }

private:
    void * m_handle;
    mutable std::string m_error;
};

} // namespace gismo
