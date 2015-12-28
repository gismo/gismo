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
createPK_BSURF( const gsTensorBSpline< 2, real_t> & bsp, 
                PK_BSURF_t & bsurf,
                bool closed_u, 
                bool closed_v);

TEMPLATE_INST bool
createPK_BCURVE(const gsBSpline<real_t>& curve, PK_BCURVE_t& bcurve);

TEMPLATE_INST bool
exportTHBsurface( const gsTHBSpline<2,real_t>& surface, PK_ASSEMBLY_t& body );

TEMPLATE_INST bool
exportMesh( const gsMesh<real_t>& mesh, PK_BODY_t& body );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsGeometry<real_t> & ggeo, std::string const & filename );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsMultiPatch<real_t> & mp, std::string const & filename );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsMesh<real_t>& mesh, std::string const & filename );

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsTHBSpline<2, real_t>& thb, std::string const & filename);


}// extensions

}// gismo
