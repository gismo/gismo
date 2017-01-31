/** @file gsSimpleOps_.cpp

    @brief Collection of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, C. Hofreither, A. Mantzaflaris
*/

#include <gsSolver/gsSimpleOps.hpp>

namespace gismo
{
    
TEMPLATE_INST void dampedRichardsonSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1.));

TEMPLATE_INST void jacobiSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

TEMPLATE_INST void dampedJacobiSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(0.5));

TEMPLATE_INST void gaussSeidelSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

TEMPLATE_INST void reverseGaussSeidelSweep(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

TEMPLATE_INST void gaussSeidelSingleBlock(const gsSparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, gsVector<index_t>& DoFs);


} // namespace gismo
