/** @file gsSimplePreconditioners_.cpp

    @brief Collection of some simple preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, C. Hofreither, A. Mantzaflaris, S. Takacs
*/

#include <gsSolver/gsSimplePreconditioners.hpp>

namespace gismo
{

namespace internal
{

TEMPLATE_INST void gaussSeidelSweep(gsSparseMatrix<real_t>::Nested& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
TEMPLATE_INST void reverseGaussSeidelSweep(gsSparseMatrix<real_t>::Nested& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

} // namespace internal

} // namespace gismo
