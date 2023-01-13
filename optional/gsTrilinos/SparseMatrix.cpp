/** @file SparseMatrix.cpp

    @brief Wrapper for Trilinos/Epetra sparse matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsTrilinos/SparseMatrix.h>
#include <gsMpi/gsMpi.h>
#include <gsTrilinos/gsTrilinosHeaders.h>


namespace gismo
{

namespace trilinos
{


struct SparseMatrixPrivate
{
    typedef real_t Scalar;
    typedef util::conditional<util::is_same<Scalar,double>::value, Epetra_CrsMatrix,
                              Tpetra::CrsMatrix<Scalar,int,int> >::type Matrix;
    //Epetra_FECrsMatrix matrix;
    
    typedef Epetra_CrsMatrix Epetra_Matrix;
    
    /// A sparse matrix object in Trilinos 
    Teuchos::RCP<Matrix> matrix;
};

SparseMatrix::SparseMatrix() : my(new SparseMatrixPrivate)
{ }
    
SparseMatrix::SparseMatrix(const gsSparseMatrix<real_t,RowMajor> & sp, const int rank)
: my(new SparseMatrixPrivate)
{
    // In G+Smo sparse matrices are usually stored as column-major
    // (rows are compressed).  In Epetra we work row-wise (compressed
    // columns), therefore any ColMajor sp's have to produce a copy

#ifdef HAVE_MPI
    Epetra_MpiComm comm (gsMpi::init().worldComm());
#else
    Epetra_SerialComm comm;
#endif
        
    // The type of global indices.  You could just set this to int,
    // but we want the example to work for Epetra64 as well.
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    // Epetra was compiled only with 64-bit global index support,
    typedef long long global_ordinal_type;
#else
    // Epetra has 32-bit global indices. If
    // EPETRA_NO_64BIT_GLOBAL_INDICES is defined, then no support
    // 64-bit indices.
    typedef int global_ordinal_type;
#endif

    // The number of rows and columns in the matrix.
    global_ordinal_type glbRows       = sp.rows();
    const global_ordinal_type locRows = glbRows;
    comm.Broadcast(&glbRows, 1, rank);
        
    GISMO_ENSURE( comm.MyPID() == rank || 0 == locRows,
                  "Only Processor "<<rank<<" can fill in entries: "
                  "size of local matrix must be zero on other Processors.");
        
    GISMO_ASSERT( sp.isCompressed(), "Need compressed matrix for now");

    // Construct a map with all the rows on processor "rank"
    Epetra_Map map0(glbRows, locRows, 0, comm);

    // Collect the number of nonzero entries per row of sp
    gsVector<global_ordinal_type>  nnzPerRow(locRows);
    for(global_ordinal_type i=0; i!=locRows; ++i)
        nnzPerRow[i] = sp.innerVector(i).nonZeros();

    // This distributed matrix is located entirely on proccessor "rank"
    // Note: we should have used "View", but there seems to be some bug
    // (known restriction: for each row, we can insert values only once)
    Epetra_CrsMatrix _sp0(Copy, map0, nnzPerRow.data(), true);
        
    // Fill in _sp0 at processor "rank"
    int err_code = 0;
    std::vector<global_ordinal_type> tmp;
    for (global_ordinal_type r = 0; r != locRows; ++r)
    {
        const index_t oind = *(sp.outerIndexPtr()+r);
        tmp.resize(nnzPerRow[r]);
        // convert index_t to global_ordinal_type
        copy_n(sp.innerIndexPtr()+oind, nnzPerRow[r], tmp.data());
        err_code = _sp0.InsertGlobalValues(r, nnzPerRow[r],
                                           sp.valuePtr()+oind,
                                           // in case global_ordinal_type == index_t
                                           //sp.innerIndexPtr()+oind
                                           tmp.data()
                                          );
        GISMO_ASSERT(0 == err_code,
                     "InsertGlobalValues failed with err_code="<<err_code);
    }
        
    err_code = _sp0.FillComplete();
    GISMO_ASSERT(0 == err_code, "FillComplete failed with err_code="<<err_code);

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    Epetra_Map map(glbRows, 0, comm);

    // We've created a sparse matrix _sp0 whose rows live entirely on MPI
    // Process 0.  Now we want to distribute it over all the processes.
    // Redistributin is NOT in place, from matrix _sp0 to *matrix
    Epetra_Export exporter(map0, map);
    //my->matrix.reset( new Epetra_CrsMatrix(_sp0, exporter) ); //bug?
    // /*// Equivalent
    my->matrix.reset( new Epetra_CrsMatrix(Copy, map, true) );
    err_code = my->matrix->Export(_sp0, exporter, Insert);
    err_code = my->matrix->FillComplete();
    err_code = my->matrix->OptimizeStorage();
    GISMO_UNUSED(err_code);
    //*/
}

    
SparseMatrix::~SparseMatrix() { delete my; }

/*
  Epetra_BlockMap SparseMatrix::map() const
  {
  return my->matrix->Map();
  }
*/

std::pair<index_t,index_t> SparseMatrix::dim() const
{
    if (util::is_same<index_t,int>::value)
        return std::make_pair(my->matrix->NumGlobalRows(),
                              my->matrix->NumGlobalCols());
    else
        return std::make_pair(my->matrix->NumGlobalRows64(),
                              my->matrix->NumGlobalCols64());
}

std::pair<index_t,index_t> SparseMatrix::mydim() const
{
    return std::make_pair(my->matrix->NumMyRows(),
                          my->matrix->NumMyCols());
}

index_t SparseMatrix::cols() const
{
    if (util::is_same<index_t,int>::value)
        return my->matrix->NumGlobalCols();
    else
        return my->matrix->NumGlobalCols64();
}

index_t SparseMatrix::mycols() const
{
    return my->matrix->NumMyCols();
}

index_t SparseMatrix::rows() const
{
    if (util::is_same<index_t,int>::value)
        return my->matrix->NumGlobalRows();
    else
        return my->matrix->NumGlobalRows64();
}

index_t SparseMatrix::myrows() const
{
    return my->matrix->NumMyRows();
}

index_t SparseMatrix::nonzeros() const
{
    if (util::is_same<index_t,int>::value)
        return my->matrix->NumGlobalNonzeros();
    else
        return my->matrix->NumGlobalNonzeros64();
}

index_t SparseMatrix::mynonzeros() const
{
    return my->matrix->NumMyNonzeros();
}

void SparseMatrix::copyTo(gsSparseMatrix<real_t,RowMajor> & sp, const int rank) const
{
/*
  Epetra_MpiComm comm (gsMpi::init().worldComm() );

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
  
  tmp.ExtractCopy(gsVec.data());
  //my->matrix->ExtractGlobalRowCopy
  }
*/
}

Epetra_CrsMatrix * SparseMatrix::get() const
{
    return my->matrix.get();
}

Teuchos::RCP<Epetra_CrsMatrix> SparseMatrix::getRCP() const
{
    return my->matrix;
}

void SparseMatrix::print() const
{
    gsInfo << "Processor No. " << gsMpi::init().worldRank() << "\n" << *get();    
}

}//namespace trilinos

}// namespace gismo
