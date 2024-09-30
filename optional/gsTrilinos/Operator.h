/** @file SparseMatrix.h

    @brief Wrapper for Trilinos/Epetra operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsTrilinos/gsTrilinosForwardDecl.h>

namespace gismo
{

namespace trilinos
{


class OperatorPrivate;

class GISMO_EXPORT Operator
{
public:

    explicit Operator(const gsSparseMatrix<real_t> & sp);

    ~Operator();
    
    Epetra_Operator * get() const;

    void print() const;
    
private:
    Operator(const Operator& other);
    Operator& operator=(const Operator& other);
    
private:

    OperatorPrivate * my;
};


}//namespace trilinos

}// namespace gismo
