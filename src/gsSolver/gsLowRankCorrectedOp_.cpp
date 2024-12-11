/** @file gsLowRankCorrectedOp_.cpp

    @brief Provides the inverse of \f$Ainv^{-1} + U Q^{-1} V^T\f$ with \f$U Q^{-1} V^T\f$ being a low rank matrix using the Sherman Morrisson Woodburry formula

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#include <gsSolver/gsLowRankCorrectedOp.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsLowRankCorrectedOp<real_t>;

}
