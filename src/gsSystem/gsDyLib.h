/** @file gsDyLib.h

    @brief Handle to a single dynamically loaded gismo module library.

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
   @brief Owns one loaded module library (.so/.dll).

   The library is opened on construction and closed on destruction.
   The class is move-only, so that gsModuleManager can keep the
   collection of loaded modules alive.

   @ingroup System
*/
class GISMO_EXPORT gsDyLib
{
public:
    /// Creates a handle that holds no library
    gsDyLib() : m_handle(NULL) { }

    /// Opens the library at @a path; check isOpen() afterwards
    explicit gsDyLib(const std::string & path);

    /// Closes the library, if one is held
    ~gsDyLib();

    /// Takes over the library of @a other, which is left empty
    gsDyLib(gsDyLib && other) noexcept
    : m_handle(other.m_handle), m_error(other.m_error)
    { other.m_handle = NULL; }

    /// Closes the library held by this object, then takes over the
    /// library of @a other, which is left empty
    gsDyLib & operator=(gsDyLib && other) noexcept;

    gsDyLib(const gsDyLib &) = delete;
    gsDyLib & operator=(const gsDyLib &) = delete;

    /// True if a library is held
    bool isOpen() const { return NULL != m_handle; }

    /// The loader's error message if loading failed
    const std::string & error() const { return m_error; }

    /// Address of the exported symbol @a name, or NULL if not found
    void * symbol(const char * name) const;

private:
    void *      m_handle;
    std::string m_error;
};

} // namespace gismo
