/** @file SparseMatrix.cpp

    @brief Wrapper for Trilinos/Epetra sparse matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsMpi/gsMpiHelper.h>
#include <gsTrilinos/gsTrilinosHeaders.h>
#include <gsTrilinos/SparseMatrix.h>

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

    SparseMatrixPrivate() { }

        
    SparseMatrixPrivate(const gsSparseMatrix<real_t> & A)
    {
        // In gismo sparse matrices are stored as column-major (rows
        // are compressed).  In Epetra we work row-wise (compressed
        // columns), therefore a transposition is needed
        gsSparseMatrix<real_t> AT = A.transpose();

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
        
        // The number of rows and columns in the matrix.
        const global_ordinal_type numGlobalElements =
            static_cast<global_ordinal_type> ( AT.rows() );
        GISMO_ENSURE( comm.MyPID() == 0 || 0 == numGlobalElements,
                      "Only Processor 0 can fill in entries");

        GISMO_ASSERT( AT.isCompressed(), "Need compressed matrix for now"); // todo

        // Construct a Map that puts approximately the same number of
        // equations on each processor.
        Epetra_Map map (numGlobalElements, 0, comm);

        // Create a Epetra sparse matrix whose rows have distribution
        // given by the Map.  The max number of entries per row is
        // given by the numbers in nEntriesPerRow.
        global_ordinal_type  nEntriesPerRow[numGlobalElements];
        for(int i=0; i<numGlobalElements; i++)
            nEntriesPerRow[i] = AT.innerVector(i).nonZeros();
        matrix.reset( new Epetra_CrsMatrix(Copy, map, nEntriesPerRow, true) );

        int err_code = 0;
        
        for (global_ordinal_type r = 0; r != numGlobalElements; ++r )
        {
            const int oind = AT.outerIndexPtr()[r];
            err_code = matrix->InsertGlobalValues (r, nEntriesPerRow[r],
                                                   AT.valuePtr()+oind,
                                                   AT.innerIndexPtr()+oind);
            GISMO_ASSERT(0 == err_code,
                         "matrix->InsertGlobalValues failed with err_code="<<err_code);
        }
        
        err_code = matrix->FillComplete ();

        // If any process failed to insert at least one entry, throw.
        int gl_err = 0;
        (void) comm.MaxAll (&err_code, &gl_err, 1);
        GISMO_ENSURE(0 == gl_err, "Some process failed constructing Trilinos matrix");

        /*        
        // Get the list of global indices that this process owns.  In this
        // case, this is unnecessary, because we know that we created a
        // contiguous Map (see above).  (Thus, we really only need the min
        // and max global index on this process.)  However, in general, we
        // don't know what global indices the Map owns, so if we plan to add
        // entries into the sparse matrix using global indices, we have to
        // get the list of global indices this process owns.
        const int numMyElements = map.NumMyElements ();

        global_ordinal_type * myGlobalElements = NULL;

        #ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
        myGlobalElements = map.MyGlobalElements64 ();
        #else
        myGlobalElements = map.MyGlobalElements ();
        #endif

        // In general, tests like this really should synchronize across all
        // processes.  However, the likely cause for this case is a
        // misconfiguration of Epetra, so we expect it to happen on all
        // processes, if it happens at all.
        if (numMyElements > 0 && myGlobalElements == NULL)
        {
            throw std::logic_error ("Failed to get the list of global indices");
        }
        */
    }
    
    /// A sparse matrix object in Trilinos 
    memory::shared_ptr<Epetra_Matrix> matrix;
};

SparseMatrix::SparseMatrix() : my(new SparseMatrixPrivate)
{ }
    
SparseMatrix::SparseMatrix(const gsSparseMatrix<> & sp)
: my(new SparseMatrixPrivate(sp))
{ }

    
SparseMatrix::~SparseMatrix() { delete my; }

/*
Epetra_BlockMap SparseMatrix::map() const
{
    return my->matrix->Map();
}
*/

void SparseMatrix::copyTo(gsSparseMatrix<> & sp) const
{
    // to do
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
    gsMPIHelper & mpi_helper = gsMPIHelper::instance();
    gsInfo << "Processor No. " << mpi_helper.rank() << "\n" << *get();    
}

}//namespace trilinos

}// namespace gismo
