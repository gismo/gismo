

#include "gsTrilinosHeaders.h"
#include "SparseMatrix.h"

namespace gismo
{

namespace trilinos
{


class SparseMatrixPrivate
{
    friend class SparseMatrix;
    
    SparseMatrixPrivate()
    : column_space_map (new Epetra_Map (0, 0, comm) ),
      matrix (new Epetra_FECrsMatrix(View, *column_space_map, *column_space_map, 0))
    { }

    // for testing
    Epetra_SerialComm comm;

    /// Epetra Trilinos mapping of the matrix columns that assigns
    /// parts of the matrix to the individual processes.
    memory::shared_ptr<Epetra_Map> column_space_map;
    
    /// A sparse matrix object in Trilinos 
    memory::shared_ptr<Epetra_FECrsMatrix> matrix;    
};

SparseMatrix::SparseMatrix() : my(new SparseMatrixPrivate)
{ }
    
SparseMatrix::SparseMatrix(const gsSparseMatrix<> & sp)
: my(new SparseMatrixPrivate)
{
    // fill in
    my->matrix->FillComplete();
}

    
SparseMatrix::~SparseMatrix() { delete my; }


}//namespace trilinos

}// namespace gismo
