/** @file gsBlockOp_.cpp

    @brief Simple class create a block preconditioner structure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <gsSolver/gsBlockOp.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsBlockOp<real_t>;
CLASS_TEMPLATE_INST gsBlockOp<std::complex<real_t> >;

} // namespace gismo

