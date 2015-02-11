/** @file gsKronecker.h

    @brief Provides functions and classes for working with Kronecker products of matrices and operators.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither
*/
#pragma once

#include <vector>
#include <gsMatrix/gsMatrix.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/**
  Compute the application of the Kronecker product
    kron(ops[0], ops[1], ..., ops[n-1]) * x
  and store it in \a result without computing the large Kronecker product matrix itself.
*/
void applyKronecker(const std::vector< gsLinearOperator* > & ops, const gsMatrix<>& x, gsMatrix<>& result);

}

