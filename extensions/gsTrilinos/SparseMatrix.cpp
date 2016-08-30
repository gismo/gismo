/** @file SparseMatrix.cpp

    @brief Wrapper for Trilinos/Epetra sparse matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsTrilinos/SparseMatrix.h>
#include <gsMpi/gsMpiHelper.h>
#include <gsTrilinos/gsTrilinosHeaders.h>

#include <gsCore/gsLinearAlgebra.h>


namespace gismo
{

namespace trilinos
{


class SparseMatrixPrivate
{
    typedef Epetra_CrsMatrix Epetra_Matrix;
    //Epetra_FECrsMatrix matrix;
    
    friend class SparseMatrix;
    
    /// A sparse matrix object in Trilinos 
    memory::shared_ptr<Epetra_Matrix> matrix;
};

SparseMatrix::SparseMatrix() : my(new SparseMatrixPrivate)
{ }
    
SparseMatrix::SparseMatrix(const gsSparseMatrix<> & sp, const int rank)
: my(new SparseMatrixPrivate)
{
#ifdef HAVE_MPI
        Epetra_MpiComm comm (gsMPIHelper::instance().getCommunicator() );
#else
        Epetra_SerialComm comm;
#endif
        
        // The type of global indices.  You could just set this to int,
        // but we want the example to work for Epetra64 as well.
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
        // Epetra was compiled only with 64-bit global index support, so use
        // 64-bit global indices.
        typedef long long global_ordinal_type;
#else
        // Epetra was compiled with 32-bit global index support.  If
        // EPETRA_NO_64BIT_GLOBAL_INDICES is defined, it does not also
        // support 64-bit indices.
        typedef int global_ordinal_type;
#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

        // In G+Smo sparse matrices are stored as column-major (rows
        // are compressed).  In Epetra we work row-wise (compressed
        // columns), therefore a transposition is needed
        gsSparseMatrix<real_t> AT = sp.transpose();

        // The number of rows and columns in the matrix.
        const global_ordinal_type locRows = AT.rows();
        global_ordinal_type glbRows       = AT.rows();
        comm.Broadcast(&glbRows, 1, 0);
        
        GISMO_ENSURE( comm.MyPID() == 0 || 0 == locRows,
                      "Only Processor 0 can fill in entries");
        
        GISMO_ASSERT( AT.isCompressed(), "Need compressed matrix for now"); // todo

        // Construct a map with all the rows on processor 0
        Epetra_Map map0(glbRows, locRows, 0, comm);

        // Collect the number of nonzero entries per row of AT
        gsVector<global_ordinal_type>  nnzPerRow(locRows);
        for(global_ordinal_type i=0; i!=locRows; ++i)
            nnzPerRow[i] = AT.innerVector(i).nonZeros();

        // This distributed matrix is stored entirely on proccessor 0
        Epetra_CrsMatrix _A0(Copy, map0, nnzPerRow.data(), true);
        
        // Fill in _A0 at processor 0
        int err_code = 0;
        for (global_ordinal_type r = 0; r != locRows; ++r)
        {
            const index_t oind = *(AT.outerIndexPtr()+r);
            err_code = _A0.InsertGlobalValues (r, nnzPerRow[r],
                                               AT.valuePtr()+oind,
                                               AT.innerIndexPtr()+oind);
            GISMO_ASSERT(0 == err_code,
                         "InsertGlobalValues failed with err_code="<<err_code);
        }
        
        err_code = _A0.FillComplete ();
        GISMO_ASSERT(0 == err_code, "FillComplete failed with err_code="<<err_code);

        // Construct a Map that puts approximately the same number of
        // equations on each processor.
        Epetra_Map map(glbRows, 0, comm);

        // We've created a sparse matrix _A0 whose rows live entirely on MPI
        // Process 0.  Now we want to distribute it over all the processes.
        my->matrix.reset( new Epetra_CrsMatrix(Copy, map, true) );

        // Redistribute the data, NOT in place, from matrix _A0 (which lives
        // entirely on Proc 0) to *matrix (which is distributed evenly over
        // the processes).
        Epetra_Export exporter(map0, map);
        err_code = my->matrix->Export(_A0, exporter, Insert);
        err_code = my->matrix->FillComplete();
        err_code = my->matrix->OptimizeStorage();
}

    
SparseMatrix::~SparseMatrix() { delete my; }

/*
Epetra_BlockMap SparseMatrix::map() const
{
    return my->matrix->Map();
}
*/

void SparseMatrix::copyTo(gsSparseMatrix<> & sp, const int rank) const
{
/*
    Epetra_MpiComm comm (gsMPIHelper::instance().getCommunicator() );
    const int myrank = comm.MyPID();
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    const long long sz = my->vec->GlobalLength64();
#else
    const int sz = my->vec->GlobalLength();
#endif
    Epetra_Map map0(sz, rank==myrank ? sz : 0, 0, comm);
    Epetra_Vector tmp(map0);
    Epetra_Export exp(my->vec->Map(), map0);
    (void)tmp.Export(*my->vec, exp, Insert);
    if ( myrank == rank )
    {
        gsVec.resize(sz);
        tmp.ExtractCopy(gsVec.data());
        //my->matrix->ExtractGlobalRowCopy
    }
*/
}

Epetra_CrsMatrix * SparseMatrix::get() const
{
    return my->matrix.get();
}

memory::shared_ptr<Epetra_CrsMatrix> SparseMatrix::getPtr()
{
    return my->matrix;
}

void SparseMatrix::print() const
{
    gsInfo << "Processor No. " << gsMPIHelper::instance().rank() << "\n" << *get();    
}

}//namespace trilinos

}// namespace gismo
