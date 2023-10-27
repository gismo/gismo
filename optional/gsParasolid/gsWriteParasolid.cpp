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

CLASS_TEMPLATE_INST gsTrimData<real_t> ;

TEMPLATE_INST bool
createPK_BSURF( const gsTensorBSpline< 2, real_t> & bsp, 
                PK_BSURF_t & bsurf,
                bool closed_u, 
                bool closed_v);

TEMPLATE_INST bool
createPK_BCURVE(const gsBSpline<real_t>& curve, PK_BCURVE_t& bcurve);

TEMPLATE_INST bool
exportTHBsurface<real_t>( const gsTHBSpline<2,real_t>& surface, PK_ASSEMBLY_t& body );

TEMPLATE_INST bool
exportTHBsurface( const gsTHBSpline<2,real_t>& surface,const std::vector<real_t>&par_boxes, PK_ASSEMBLY_t& body );

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

TEMPLATE_INST bool
gsWriteParasolid<real_t>
( const gsTHBSpline<2, real_t>& thb,const std::vector<real_t>&par_boxes, std::string const & filename);

TEMPLATE_INST bool
gsWritePK_SHEET<real_t>
( const gsTensorBSpline<2, real_t>& tp, std::string const & filename);

TEMPLATE_INST bool
getTrimCurvesAndBoundingBoxes<real_t>
( const gsTHBSpline<2, real_t>& surface,
  const std::vector<real_t>& par_boxes,
  std::vector<gsTrimData<real_t> >& trimdata);

TEMPLATE_INST bool
getParBoxAsIndexBoxInLevel(const gsTHBSplineBasis<2, real_t>& basis,unsigned lvl,
                           const std::vector<real_t>& par_box,std::vector<index_t>& index_box);

TEMPLATE_INST bool
parBoxesIntersect(const std::vector<real_t>& par_boxes);

}// extensions

}// gismo
