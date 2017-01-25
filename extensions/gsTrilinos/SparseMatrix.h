/** @file SparseMatrix.h

    @brief Wrapper for Trilinos/Epetra sparse matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsTrilinos/gsTrilinosForwardDecl.h>

namespace gismo
{

namespace trilinos
{

class SparseMatrixPrivate;

class GISMO_EXPORT SparseMatrix
{
public:

    SparseMatrix();
    
    explicit SparseMatrix(const gsSparseMatrix<real_t,RowMajor> & sp, const int rank = 0);

    ~SparseMatrix();

    std::pair<index_t,index_t> dim() const;

    std::pair<index_t,index_t> mydim() const;

    index_t cols() const;
    
    index_t mycols() const;

    index_t rows() const;
    
    index_t myrows() const;

    index_t nonzeros() const;

    index_t mynonzeros() const;
    
    void copyTo(gsSparseMatrix<real_t,RowMajor> & sp, const int rank = 0) const;

    Epetra_CrsMatrix * get() const;

    Teuchos::RCP<Epetra_CrsMatrix> getRCP() const;
    
    void print() const;

private:
    SparseMatrix(const SparseMatrix& other);
    SparseMatrix& operator=(const SparseMatrix& other);
    
private:

    SparseMatrixPrivate * my;
};


}//namespace trilinos

}// namespace gismo
