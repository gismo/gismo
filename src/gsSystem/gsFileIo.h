/** @file gsFileIo.h

    @brief Base class of file-format readers and writers provided by
    modules.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsXml.h>
#include <gsSystem/gsRegistry.h>

#include <string>

namespace gismo
{

/**
   @brief Reader and/or writer of a file format that is provided by a module.

   A module registers a factory of its derived class in
   gsRegistry<gsFileIo>, under the lower-case extension of the format
   without the dot (e.g. "3dm"). gsFileData looks the format up with
   find() when it meets an extension it has no native reader or writer
   for.

   An object made by a module's factory must be destroyed before that
   module is unloaded. gsFileData creates and destroys it within a single
   call.

   @ingroup System
*/
class GISMO_EXPORT gsFileIo
{
public:
    /// Unique pointer to a gsFileIo
    typedef memory::unique_ptr<gsFileIo> uPtr;

    virtual ~gsFileIo();

    /// Returns true if read() is implemented; false by default
    virtual bool canRead()  const;

    /// Returns true if write() is implemented; false by default
    virtual bool canWrite() const;

    /// Reads the file @a fn into the xml tree @a data. Returns true on
    /// success; the default implementation returns false.
    virtual bool read (const std::string & fn, internal::gsXmlTree & data) const;

    /// Writes the xml tree @a data to the file @a fn. Returns true on
    /// success; the default implementation returns false.
    virtual bool write(const std::string & fn, const internal::gsXmlTree & data) const;

    /// Makes the reader/writer registered for the extension @a ext (case
    /// insensitive, with or without leading dot), or returns a null
    /// pointer if there is none. If the extension is not registered,
    /// gsModuleManager::loadAll() runs once and the lookup is repeated.
    static uPtr find(const std::string & ext);
};

} // namespace gismo

// The single instance of gsRegistry<gsFileIo> lives in libgismo
// (gsRegistry_.cpp); every other translation unit links to it.
#ifndef GISMO_BUILD_LIB
#include <gsSystem/gsRegistry.hpp>
#else
#ifndef gsRegistry_EXPORT
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsRegistry<gsFileIo>;
}
#endif
#endif
