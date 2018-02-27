/** @file gsTensorTools.hpp

    @brief Provides functions for working with Kronecker products of matrices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs

    This file is based on unsupported parts of Eigen:

        // This file is part of Eigen, a lightweight C++ template library
        // for linear algebra.
        //
        // Copyright (C) 2011 Kolja Brix <brix@igpm.rwth-aachen.de>
        // Copyright (C) 2011 Andreas Platen <andiplaten@gmx.de>
        // Copyright (C) 2012 Chen-Pang He <jdh8@ms63.hinet.net>
        //
        // This Source Code Form is subject to the terms of the Mozilla
        // Public License v. 2.0. If a copy of the MPL was not distributed
        // with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    result.resizeNonZeros(0);

    if( Ar*Br == 0 || Ac*Bc == 0 )
        return result;

    // compute number of non-zeros per innervectors of result
    {
        // VectorXi is not necessarily big enough.
        Eigen::VectorXi nnzA = Eigen::VectorXi::Zero(MatrixType::IsRowMajor ? A.rows() : A.cols());
        for (index_t kA=0; kA < A.outerSize(); ++kA)
            for (InnerIterator itA(A,kA); itA; ++itA)
                nnzA(MatrixType::IsRowMajor ? itA.row() : itA.col())++;

        Eigen::VectorXi nnzB = Eigen::VectorXi::Zero(MatrixType::IsRowMajor ? B.rows() : B.cols());
        for (index_t kB=0; kB < B.outerSize(); ++kB)
            for (InnerIterator itB(B,kB); itB; ++itB)
                nnzB(MatrixType::IsRowMajor ? itB.row() : itB.col())++;

        Eigen::Matrix<int,Dynamic,Dynamic,ColMajor> nnzAB = nnzB * nnzA.transpose();
        result.reserve(Eigen::VectorXi::Map(nnzAB.data(), nnzAB.size()));
    }

    for (index_t kA=0; kA < A.outerSize(); ++kA)
    {
        for (index_t kB=0; kB < B.outerSize(); ++kB)
        {
            for (InnerIterator itA(A,kA); itA; ++itA)
            {
                for (InnerIterator itB(B,kB); itB; ++itB)
                {
                    const index_t i = itA.row() * Br + itB.row(),
                                  j = itA.col() * Bc + itB.col();
                    result.insert(i,j) = itA.value() * itB.value();
                }
            }
        }
    }

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
