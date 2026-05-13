/** @file gsPKUtils.cpp

    @brief Provides instantization of gsPKUtils functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D. Mokris
*/

#include <gsParasolid/gsPKUtils.h>
#include <gsParasolid/gsPKUtils.hpp>

namespace gismo {
namespace extensions {

TEMPLATE_INST bool
createPK_BSURF(const gsTensorBSpline< 2, real_t> & bsp, 
			   PK_BSURF_t & bsurf,
			   bool closed_u, 
			   bool closed_v);

TEMPLATE_INST bool
createPK_BCURVE(const gsBSpline<real_t>& curve, PK_BCURVE_t& bcurve);

TEMPLATE_INST bool
createPK_GEOM(const gsGeometry<real_t> & ggeo, PK_GEOM_t & pgeo);

} // namespace extensions

} // namespace gismo
