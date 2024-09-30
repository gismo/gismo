/** @file Vector.h

    @brief Wrapper for Trilinos/Epetra vector

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsTrilinos/SparseMatrix.h>

namespace gismo
{

namespace trilinos
{


class VectorPrivate;

class GISMO_EXPORT Vector
{
public:
    
    Vector();

    explicit Vector(const SparseMatrix & _map);
    
    Vector(const gsVector<real_t> & gsVec, const SparseMatrix & _map, const int rank = 0);
    
    explicit Vector(Epetra_Vector * v_ptr);
        
    ~Vector();

    size_t size() const;

    size_t mySize() const;
    
    void setConstant(const double val);

    void setFrom(const SparseMatrix & _map);
    
    void copyTo(gsVector<real_t> & gsVec, const int rank = 0) const;

    void copyTo(gsMatrix<real_t> & gsVec, const int = 0) const
    {
        gsVector<real_t> tmp;
        copyTo(tmp);
        tmp.swap(gsVec);
    }

    Epetra_MultiVector * get() const;

    Teuchos::RCP<Epetra_MultiVector> getRCP() const;
    
    void print() const;
    
private:
    Vector(const Vector& other);
    Vector& operator=(const Vector& other);
    
private:

    VectorPrivate * my;
};


}//namespace trilinos

}// namespace gismo
