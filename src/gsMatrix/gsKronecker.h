/** @file gsKronecker.h

    @brief Provides functions for working with Kronecker products of matrices.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, S. Takacs
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsExport.h>

namespace gismo
{

/// @brief Compute the Kronecker product of two sparse matrices as a sparse matrix
///
/// @ingroup Matrix
template <typename T, int _Options>
gsSparseMatrix<T,_Options> getKroneckerProduct(const gsSparseMatrix<T,_Options>& A, const gsSparseMatrix<T,_Options>& B);

/// @brief Compute the Kronecker product of two dense matrices as a dense matrix
///
/// @ingroup Matrix
template <typename T>
gsMatrix<T> getKroneckerProduct(const gsMatrix<T>& A, const gsMatrix<T>& B);

/// @brief Compute the Kronecker product of a vector of matrices
///
/// @ingroup Matrix
template <typename MatrixType>
MatrixType getKroneckerProduct(const std::vector< MatrixType >& matrices)
{
    if ( matrices.size() == 0 )
        return MatrixType();
    else if ( matrices.size() == 1 )
        return matrices[0];
    else
    {
        MatrixType result = getKroneckerProduct(matrices[0], matrices[1]);
        MatrixType tmp;
        for ( unsigned i = 2; i<matrices.size(); ++ i )
        {
            tmp.swap(result);
            result = getKroneckerProduct(tmp, matrices[i]);
        }
        return result;
    }
}

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsKronecker.hpp)
#endif


