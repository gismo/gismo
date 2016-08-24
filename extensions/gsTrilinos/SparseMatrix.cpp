/** @file SparseMatrix.cpp

    @brief Wrapper for Trilinos/Epetra sparse matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include "gsTrilinosHeaders.h"
#include "SparseMatrix.h"

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

        
    SparseMatrixPrivate(const gsSparseMatrix<real_t> & originalA)
    {
        // get info from gismo/Eigen matrix
        const int outersize = originalA.outerSize();
        const int innersize = originalA.innerSize();
        
        // In gismo sparse matrices are stored as column-major (rows are compressed).
        // In Epetra we work row-wise (compressed columns)
        GISMO_ASSERT( outersize==originalA.cols() && innersize==originalA.rows(),
                      "Matrix dimension error");
        
        gsSparseMatrix<real_t> AT = originalA.transpose();
        
#ifdef HAVE_MPI
        MPI_Init (&argc, &argv);
        Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
        Epetra_SerialComm comm;
#endif // HAVE_MPI
        
        //const int myRank = comm.MyPID ();
        //const int numProcs = comm.NumProc ();
        
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
            static_cast<global_ordinal_type>(originalA.rows());

        // Construct a Map that puts approximately the same number of
        // equations on each processor.
        const global_ordinal_type indexBase = 0;
        Epetra_Map map (numGlobalElements, indexBase, comm);

        // Get the list of global indices that this process owns.  In this
        // example, this is unnecessary, because we know that we created a
        // contiguous Map (see above).  (Thus, we really only need the min
        // and max global index on this process.)  However, in general, we
        // don't know what global indices the Map owns, so if we plan to add
        // entries into the sparse matrix using global indices, we have to
        // get the list of global indices this process owns.
        const int numMyElements = map.NumMyElements ();

        global_ordinal_type * myGlobalElements = NULL;

#       ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
        myGlobalElements = map.MyGlobalElements64 ();
#       else
        myGlobalElements = map.MyGlobalElements ();
#       endif // EPETRA_NO_32BIT_GLOBAL_INDICES

        // In general, tests like this really should synchronize across all
        // processes.  However, the likely cause for this case is a
        // misconfiguration of Epetra, so we expect it to happen on all
        // processes, if it happens at all.
        if (numMyElements > 0 && myGlobalElements == NULL)
        {
            throw std::logic_error ("Failed to get the list of global indices");
        }

        // Create a Epetra sparse matrix whose rows have distribution
        // given by the Map.  The max number of entries per row is
        // given by the numbers in nEntriesPerRow.
        global_ordinal_type  nEntriesPerRow[numGlobalElements];
        for(int i=0; i<numGlobalElements; i++)
            nEntriesPerRow[i] = AT.innerVector(i).nonZeros();
        
        Epetra_CrsMatrix A (Copy, map, nEntriesPerRow, true);

        // Local error code for use below.
        int lclerr = 0;

        // Fill the sparse matrix, one row at a time.  InsertGlobalValues
        // adds entries to the sparse matrix, using global column indices.
        // It changes both the graph structure and the values.
        std::vector<double> tmpValues;
        std::vector<global_ordinal_type> tmpColumns;
        for (int i = 0; i < numMyElements; ++i)
        {
            int globalrow = myGlobalElements[i];

            tmpValues .resize(nEntriesPerRow[globalrow]);
            tmpColumns.resize(nEntriesPerRow[globalrow]);
            
            int c = 0;
            for (gsSparseMatrix<real_t>::InnerIterator it(AT,globalrow); it; ++it)
            {
                GISMO_ASSERT(c < nEntriesPerRow[globalrow], "Sparse matrix Filling failed");
                tmpValues[c]  = it.value();
                tmpColumns[c] = it.col();
                c++;
            }

            if (lclerr == 0)
            {
                //gsInfo<<"Insert ("<<globalrow <<"," << c <<")\n"; 
                lclerr = A.InsertGlobalValues (globalrow, c, 
                                               tmpValues.data(), tmpColumns.data());
            }
            else
            {
                break;
            }
        }

        // If any process failed to insert at least one entry, throw.
        int gblerr = 0;
        (void) comm.MaxAll (&lclerr, &gblerr, 1);
        if (gblerr != 0)
        {
            throw std::runtime_error ("Some process failed to insert an entry.");
        }

        // Tell the sparse matrix that we are done adding entries to it.
        gblerr = A.FillComplete ();
        if (gblerr != 0)
        {    
            std::ostringstream os;
            os << "A.FillComplete() failed with error code " << gblerr << ".";
            throw std::runtime_error (os.str ());
        }
        
        matrix.reset( new Epetra_CrsMatrix(A) );
    }
            
/*   
     Epetra_SerialComm comm;

     /// Epetra Trilinos mapping of the matrix columns that assigns
     /// parts of the matrix to the individual processes.
     memory::shared_ptr<Epetra_Map> column_space_map;
*/
    
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

}//namespace trilinos

}// namespace gismo
