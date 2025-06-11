/** @file gsRationalTHBSpline.h

    @brief Represents a rational truncated hierarchical B-Spline patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsRationalTHBSpline.h>
#include <gsHSplines/gsRationalTHBSpline.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSpline<1,real_t> >;
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSpline<2,real_t> >;
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSpline<3,real_t> >;
    CLASS_TEMPLATE_INST internal::gsXml< gsRationalTHBSpline<4,real_t> >;

} // namespace gismo
