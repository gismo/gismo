/** @file Vector.h

    @brief Wrapper for Trilinos/Epetra vector

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>

#include <gsTrilinos/SparseMatrix.h>


// FD Trilinos
class Epetra_Vector;
class Epetra_BlockMap;

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
    
    Vector(const gsVector<> & gsVec, const SparseMatrix & _map);
    
    explicit Vector(Epetra_Vector * v_ptr);
        
    ~Vector();

    size_t size() const;
    
    void setConstant(const double val);

    void setFrom(const SparseMatrix & _map);
        
    void copyTo(gsVector<real_t> & gsVec) const;

    gsVector<real_t> print() const
    {
        gsVector<real_t> tmp;
        tmp.setConstant(11);
        copyTo(tmp);
        return tmp;
    }

    Epetra_Vector * get() const;

private:
    Vector(const Vector& other);
    Vector& operator=(const Vector& other);
    
private:

    VectorPrivate * my;
};


}//namespace trilinos

}// namespace gismo
