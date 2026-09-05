/** @file gsRationalTHBSplineBasis.h

    @brief Provides declaration of RationalTHBSplineBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsRationalTHBSplineBasis.h>
#include <gsHSplines/gsRationalTHBSplineBasis.hpp>
#include <gsIO/gsXmlRegistry.h>

namespace gismo
{
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSplineBasis<1,real_t> >;
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSplineBasis<2,real_t> >;
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSplineBasis<3,real_t> >;
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSplineBasis<4,real_t> >;


// XML dispatch registration (priorities: see gsHSplinesXmlRegistration.h)
GISMO_XML_REGISTER_GET(gsBasis<real_t>, gsRationalTHBSplineBasis<TMPLA2(1,real_t)>)
GISMO_XML_REGISTER_GET(gsBasis<real_t>, gsRationalTHBSplineBasis<TMPLA2(2,real_t)>)
GISMO_XML_REGISTER(gsBasis<real_t>, gsRationalTHBSplineBasis<TMPLA2(3,real_t)>, 230)

} // namespace gismo
