/** @file gsWriteOcct.h

    @brief Reading OpenCascade .brep files

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDebug.h>

#include <gsIO/gsXml.h>

#define gsEndl std::endl


namespace gismo {

namespace extensions {


GISMO_EXPORT bool writeOcctStep(const gsSurface<real_t> & srf, const std::string & name);
GISMO_EXPORT bool writeOcctIges(const gsSurface<real_t> & srf, const std::string & name);

} //extensions

}//namespace gismo
