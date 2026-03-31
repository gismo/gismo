/** @file gsReadOpenNurbsIo.h

    @brief Declaration of function for data input from the Rhinoceros 3DM file format.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsOpennurbs/gsReadOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.h>

#include <gsSystem/gsFileIo.h>

#pragma once

namespace gismo {

class gsOpenNurbsFile : public gsFileIo
{
public:
    gsOpenNurbsFile() : gsFileIo("3dm") { }

    virtual bool read(std::string fname, FileData & data)
    { return extensions::gsReadOpenNurbs(fname,data); }
    
    virtual bool write(std::string fname, FileData & data)
    { return extensions::gsWriteOpenNurbs(fname,data); }
};


} // namespace gismo

