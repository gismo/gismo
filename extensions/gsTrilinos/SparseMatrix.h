

#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

namespace trilinos
{


class SparseMatrixPrivate;

class GISMO_EXPORT SparseMatrix
{
public:

    SparseMatrix();
    
    SparseMatrix(const gsSparseMatrix<> & sp);

    ~SparseMatrix();

private:

    SparseMatrixPrivate * my;
};


}//namespace trilinos

}// namespace gismo
