/** @file gsMappedBasis.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>

#include <gsC1Basis/gsC1SplineBase.h>

#include <gsUnstructuredSplines/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/gsApproxC1Spline.hpp>

// #include <gsUnstructuredSplines/gsC1Andrea.h>
// #include <gsUnstructuredSplines/gsC1Andrea.hpp>

namespace gismo
{

// CLASS_TEMPLATE_INST gsC1Basis<1,real_t> ;
CLASS_TEMPLATE_INST gsC1SplineBase<2,real_t> ;
CLASS_TEMPLATE_INST gsApproxC1Spline<2,real_t>;
// CLASS_TEMPLATE_INST gsC1Basis<3,real_t> ;

} // end namespace gismo
