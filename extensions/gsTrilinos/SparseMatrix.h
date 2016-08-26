/** @file SparseMatrix.h

    @brief Wrapper for Trilinos/Epetra sparse matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsForwardDeclarations.h>

// FD Trilinos
class Epetra_CrsMatrix;
class Epetra_BlockMap;

namespace gismo
{

namespace trilinos
{


class SparseMatrixPrivate;

class GISMO_EXPORT SparseMatrix
{
public:

    SparseMatrix();
    
    explicit SparseMatrix(const gsSparseMatrix<> & sp);

    ~SparseMatrix();

    //Epetra_BlockMap map() const;
    
    void copyTo(gsSparseMatrix<> & sp) const;

    Epetra_CrsMatrix * get() const;

    void print() const;

private:
    SparseMatrix(const SparseMatrix& other);
    SparseMatrix& operator=(const SparseMatrix& other);
    
private:

    SparseMatrixPrivate * my;
};


}//namespace trilinos

}// namespace gismo
