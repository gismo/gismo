/** @file gsApproxC1Edge.hpp

    @brief Creates the (approx) C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#include <gsCore/gsTemplateTools.h>

#include <gsUnstructuredSplines/gsApproxC1Edge.h>
#include <gsUnstructuredSplines/gsApproxC1Edge.hpp>

namespace gismo
{

//CLASS_TEMPLATE_INST gsApproxC1Edge<1,real_t> ;
CLASS_TEMPLATE_INST gsApproxC1Edge<2,real_t> ;
//CLASS_TEMPLATE_INST gsApproxC1Edge<3,real_t> ;

} // end namespace gismo