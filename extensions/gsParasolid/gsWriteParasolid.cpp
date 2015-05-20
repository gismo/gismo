/** @file gsWriteParasolid.cpp

    @brief Provides instantization of gsWriteParasolid functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsTemplateTools.h>

#include <gsParasolid/gsWriteParasolid.h>
#include <gsParasolid/gsWriteParasolid.hpp>

namespace gismo {
namespace extensions {

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsGeometry<real_t> & ggeo, std::string const & filename );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsMultiPatch<real_t> & mp, std::string const & filename );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsMesh<real_t>& mesh, std::string const & filename );


}// extensions

}// gismo
