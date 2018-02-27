/** @file gsTensorTools.hpp

    @brief Provides functions for working with Kronecker products of matrices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs

*/
#include <gsTensor/gsTensorTools.h>

namespace gismo
{

template <typename T, int _Options>
gsSparseMatrix<T,_Options> getKroneckerProduct(const gsSparseMatrix<T,_Options>& A, const gsSparseMatrix<T,_Options>& B)
{
    typedef gsSparseMatrix<T,_Options> MatrixType;
    typedef typename MatrixType::InnerIterator InnerIterator;

    const index_t Ar = A.rows(), Ac = A.cols();
    const index_t Br = B.rows(), Bc = B.cols();

    MatrixType result(Ar*Br, Ac*Bc);

    if( Ar*Br == 0 || Ac*Bc == 0 )
        return result;

    index_t nz = A.nonZeros() * B.nonZeros();

    gsSparseEntries<T> se;
    se.reserve(nz);

    for (index_t i=0; i < A.outerSize(); ++i)
    {
        for (index_t j=0; j < B.outerSize(); ++j)
        {
            for (InnerIterator ii(A,i); ii; ++ii)
            {
                for (InnerIterator jj(B,j); jj; ++jj)
                {
                    const T val = ii.value() * jj.value();
                    if (val!=(T)0)
                        se.add(
                            ii.row() * Br + jj.row(),
                            ii.col() * Bc + jj.col(),
                            val
                        );
                }
            }
        }
    }

    result.setFrom(se);
    result.makeCompressed();
    return result;
}

template <typename T>
gsMatrix<T> getKroneckerProduct(const gsMatrix<T>& A, const gsMatrix<T>& B)
{
    const index_t Ar = A.rows(), Ac = A.cols();
    const index_t Br = B.rows(), Bc = B.cols();

    gsMatrix<T> result(Ar*Br, Ac*Bc);

    for (index_t i=0; i < Ar; ++i)
        for (index_t j=0; j < Ac; ++j)
            for (index_t ii=0; ii < Br; ++ii)
                for (index_t jj=0; jj < Bc; ++jj)
                    result(i*Br+ii,j*Bc+jj) = A(i,j) * B(ii,jj);

    return result;
}

} // namespace gismo
