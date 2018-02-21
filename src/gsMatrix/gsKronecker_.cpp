/** @file gsKronecker_.cpp

    @brief Provides functions for working with Kronecker products of matrices.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/

#include <gsMatrix/gsKronecker.hpp>

namespace gismo
{

TEMPLATE_INST
gsSparseMatrix<real_t> getKroneckerProduct(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& B);

TEMPLATE_INST
gsSparseMatrix<real_t, RowMajor> getKroneckerProduct(const gsSparseMatrix<real_t, RowMajor>& A, const gsSparseMatrix<real_t, RowMajor>& B);

TEMPLATE_INST
gsSparseMatrix<std::complex<real_t> > getKroneckerProduct(const gsSparseMatrix<std::complex<real_t> >& A, const gsSparseMatrix<std::complex<real_t> >& B);

TEMPLATE_INST
gsSparseMatrix<std::complex<real_t>, RowMajor> getKroneckerProduct(const gsSparseMatrix<std::complex<real_t>, RowMajor>& A, const gsSparseMatrix<std::complex<real_t>, RowMajor>& B);


TEMPLATE_INST
gsMatrix<real_t> getKroneckerProduct(const gsMatrix<real_t>& A, const gsMatrix<real_t>& B);

TEMPLATE_INST
gsMatrix<std::complex<real_t> > getKroneckerProduct(const gsMatrix<std::complex<real_t> >& A, const gsMatrix<std::complex<real_t> >& B);


}
