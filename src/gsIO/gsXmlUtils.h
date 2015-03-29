/** @file gsXmlUtils.h

    @brief Provides declaration of input/output XML utilities struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsXml.h>

#include <string>
#include <cstring>


namespace gismo {

namespace internal {

/* generic utils, will not be exported

/// Helper to read tensor bases
template<class Object>
Object * getTensorBasisFromXml ( gsXmlNode * node);

/// Helper to put tensor bases to XML
template<class Object>
gsXmlNode * putTensorBasisToXml ( Object const & obj, gsXmlTree & data);

/// Helper to read rational bases
template<class Object>
Object * getRationalBasisFromXml ( gsXmlNode * node);

/// Helper to put rational bases to XML
template<class Object>
gsXmlNode * putRationalBasisToXml ( Object const & obj, gsXmlTree & data);

*/

}// end namespace internal

}// end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXmlUtils.hpp)
//#include GISMO_HPP_HEADER(gsXmlGenericUtils.hpp)
#endif
