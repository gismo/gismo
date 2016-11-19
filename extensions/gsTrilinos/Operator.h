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
class Epetra_Operator;
class Epetra_BlockMap;

namespace gismo
{

namespace trilinos
{


class OperatorPrivate;

class GISMO_EXPORT Operator
{
public:

    explicit Operator(const gsSparseMatrix<> & sp);

    ~Operator();
    
    Epetra_Operator * get() const;

    void print() const;

public:

    // int apply(Vector & X, Vector & Y) const;
    // ..

    
private:
    Operator(const Operator& other);
    Operator& operator=(const Operator& other);
    
private:

    OperatorPrivate * my;
};


}//namespace trilinos

}// namespace gismo
