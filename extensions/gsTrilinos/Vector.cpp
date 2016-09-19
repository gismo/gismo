/** @file Vector.cpp

    @brief Wrapper for Trilinos/Epetra vector

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsTrilinos/Vector.h>

#include <gsMpi/gsMpi.h>
#include <gsTrilinos/gsTrilinosHeaders.h>

//#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>


namespace gismo
{

namespace trilinos
{

class VectorPrivate
{
    friend class Vector;
/*   
    Epetra_SerialComm comm;
    /// Epetra Trilinos mapping of the matrix columns that assigns
    /// parts of the matrix to the individual processes.
    memory::shared_ptr<Epetra_Map> column_space_map;
*/
    
    /// A vector object in Trilinos 
    memory::shared_ptr<Epetra_Vector> vec;
};

Vector::Vector(const SparseMatrix & _map)
: my(new VectorPrivate)
{
    my->vec.reset(new Epetra_Vector(_map.get()->OperatorDomainMap()));
}

Vector::Vector(const gsVector<> & gsVec, const SparseMatrix & _map, const int rank)
: my(new VectorPrivate)
{
#   ifdef HAVE_MPI
    Epetra_MpiComm comm (gsMpi::init().worldComm() );
#   else
    Epetra_SerialComm comm;
#   endif

    // The number of rows and columns in the matrix.
    index_t glbRows       = gsVec.rows();
    const index_t locRows = glbRows;
    comm.Broadcast(&glbRows, 1, rank);
    GISMO_ENSURE( comm.MyPID() == rank || 0 == locRows,
                  "Only Processor "<<rank<<" can fill in entries");
    
    // Create a temporary Epetra_Vector on Proc "rank"
    Epetra_Map map0(gsVec.rows(), locRows, 0, comm);
    Epetra_Vector tmp(View, map0, const_cast<real_t *>(gsVec.data()) );

    // Initialize the distributed Epetra_Vector
    const Epetra_Map & map = _map.get()->OperatorRangeMap();
    my->vec.reset( new Epetra_Vector(map) );

    int err_code = 0;
    // Redistribute the vector data
    Epetra_Import importer(map0, map);
    err_code = my->vec->Export(tmp, importer, Insert);
    GISMO_ENSURE(0==err_code, "Something went terribly wrong");
}

Vector::Vector(Epetra_Vector * v_ptr) : my(new VectorPrivate)
{
    my->vec.reset(v_ptr, null_deleter<Epetra_Vector> );
}

Vector::Vector() : my(new VectorPrivate) { }

Vector::~Vector() { delete my; }

void Vector::setConstant(const double val)
{
    my->vec->PutScalar(val);
}

void Vector::setFrom(const SparseMatrix & A)
{
    my->vec.reset( new Epetra_Vector(A.get()->OperatorRangeMap()) );
}

size_t Vector::size() const 
{ 
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    return (my->vec ? my->vec->GlobalLength64() : 0 ); 
#else
    return (my->vec ? my->vec->GlobalLength()   : 0 ); 
#endif
}

size_t Vector::mySize() const 
{ 
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    return (my->vec ? my->vec->MyLength64() : 0 ); 
#else
    return (my->vec ? my->vec->MyLength()   : 0 ); 
#endif
}

void Vector::copyTo(gsVector<real_t> & gsVec, const int rank) const
{
    Epetra_MpiComm comm (gsMpi::init().worldComm());
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
    }
}

Epetra_Vector * Vector::get() const
{
    return my->vec.get();
}

void Vector::print() const
{
    gsInfo << "Processor No. " << gsMpi::init().worldRank() << "\n" << *get();    
}


}//namespace trilinos

}// namespace gismo
