/** @file gsOpenNurbsFile.h

    @brief File-format class for the Rhinoceros 3DM format, based on OpenNURBS.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
*/

#pragma once

#include <gsSystem/gsFileIo.h>

namespace gismo
{

/**
   @brief Reads and writes Rhinoceros 3DM files through OpenNURBS.

   The gsOpennurbs module entry registers this class in
   @c gsRegistry<gsFileIo> under the key @c "3dm", so that gsFileData
   reaches it through @c gsFileIo::find("3dm").

   Reading appends to the @c xml root of the tree what
   extensions::gsReadOpenNurbs extracts. Writing handles curves,
   surfaces and meshes, as extensions::gsWriteOpenNurbs does.
*/
class gsOpenNurbsFile : public gsFileIo
{
public:
    /// @brief Reading is supported.
    /// @return true
    bool canRead() const GISMO_OVERRIDE;

    /// @brief Writing is supported.
    /// @return true
    bool canWrite() const GISMO_OVERRIDE;

    /// @brief Reads the 3DM file @a fn into @a data.
    /// @param fn   path of the @c .3dm file
    /// @param data tree to which the extracted objects are appended
    /// @return true on success, false if the file could not be read
    bool read(const std::string & fn, internal::gsXmlTree & data) const GISMO_OVERRIDE;

    /// @brief Writes the curves, surfaces and meshes of @a data to the 3DM file @a fn.
    /// @param fn   path of the @c .3dm file to create
    /// @param data tree holding the objects to write
    /// @return true on success, false if the file could not be written
    bool write(const std::string & fn, const internal::gsXmlTree & data) const GISMO_OVERRIDE;
};

} // namespace gismo
