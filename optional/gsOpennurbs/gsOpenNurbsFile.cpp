/** @file gsOpenNurbsFile.cpp

    @brief Implementation of the OpenNURBS file-format class and the gsOpennurbs module entry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
*/

#include <gsOpennurbs/gsOpenNurbsFile.h>
#include <gsOpennurbs/gsReadOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.h>
#include <gsSystem/gsModule.h>

namespace gismo
{

bool gsOpenNurbsFile::canRead()  const { return true; }
bool gsOpenNurbsFile::canWrite() const { return true; }

bool gsOpenNurbsFile::read(const std::string & fn, internal::gsXmlTree & data) const
{ return extensions::gsReadOpenNurbs(fn.c_str(), data); }

bool gsOpenNurbsFile::write(const std::string & fn, const internal::gsXmlTree & data) const
{ return extensions::gsWriteOpenNurbs(fn.c_str(), data); }

GISMO_MODULE_ENTRY(gsOpennurbs, GISMO_VERSION)
{
    gsRegistry<gsFileIo>::get().add("3dm", &gsRegistryFactory<gsFileIo, gsOpenNurbsFile>);
}

} // namespace gismo
