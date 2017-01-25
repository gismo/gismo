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
#include <gsCore/gsLinearAlgebra.h>


namespace gismo
{

namespace trilinos
{

struct VectorPrivate
{
    typedef real_t Scalar;
    typedef util::conditional<util::is_same<Scalar,double>::value, Epetra_MultiVector,
                              Tpetra::MultiVector<Scalar,int,int> >::type MVector;
    
    /// A vector object in Trilinos
    Teuchos::RCP<MVector> vec;
};

Vector::Vector(const SparseMatrix & _map)
: my(new VectorPrivate)
{
    my->vec.reset(new VectorPrivate::MVector(_map.get()->OperatorDomainMap(), 1));
}

Vector::Vector(const gsVector<> & gsVec, const SparseMatrix & _map, const int rank)
: my(new VectorPrivate)
{
#   ifdef HAVE_MPI
    Epetra_MpiComm comm (gsMpi::init().worldComm() );
#   else
    Epetra_SerialComm comm;
#   endif

#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    typedef long long global_ordinal_type;
#else
    typedef int global_ordinal_type;
#endif
    
    // The number of rows and columns in the matrix.
    global_ordinal_type glbRows = gsVec.rows();
    const int locRows = glbRows;
    comm.Broadcast(&glbRows, 1, rank);
    GISMO_ENSURE( comm.MyPID() == rank || 0 == locRows,
                  "Only Processor "<<rank<<" can fill in entries");
    
    // Create a temporary Epetra_Vector on Proc "rank"
    Epetra_Map map0(glbRows, locRows, 0, comm);
    Epetra_Vector tmp(View, map0, const_cast<real_t *>(gsVec.data()) );
    
    // Initialize the distributed Epetra_Vector
    const Epetra_Map & map = _map.get()->OperatorRangeMap();
    my->vec.reset( new VectorPrivate::MVector(map, 1) );

    int err_code = 0;
    // Redistribute the vector data
    Epetra_Import importer(map0, map);
    err_code = my->vec->Export(tmp, importer, Insert);
    GISMO_ENSURE(0==err_code, "Something went terribly wrong");
}

Vector::Vector(Epetra_Vector * v_ptr) : my(new VectorPrivate)
{
    my->vec.reset(v_ptr, false ); // has_ownership==false
}

Vector::Vector() : my(new VectorPrivate) { }

Vector::~Vector() { delete my; }

void Vector::setConstant(const double val)
{
    my->vec->PutScalar(val);
}

void Vector::setFrom(const SparseMatrix & A)
{
    my->vec.reset( new VectorPrivate::MVector(A.get()->OperatorRangeMap(), 1) );
}

size_t Vector::size() const 
{ 
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    return (my->vec.is_null() ? 0 : my->vec->GlobalLength64() ); 
#else
    return (my->vec.is_null() ? 0 : my->vec->GlobalLength() ); 
#endif
}

size_t Vector::mySize() const 
{ 
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    return (my->vec.is_null() ? 0 : my->vec->MyLength64() ); 
#else
    return (my->vec.is_null() ? 0 : my->vec->MyLength() ); 
#endif
}

void Vector::copyTo(gsVector<real_t> & gsVec, const int rank) const
{
#ifdef HAVE_MPI
    Epetra_MpiComm comm (gsMpi::init().worldComm());
#else
    Epetra_SerialComm comm;
#endif
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

Epetra_MultiVector * Vector::get() const
{
    return my->vec.get();
}

Teuchos::RCP<Epetra_MultiVector> Vector::getRCP() const
{
    return my->vec;
}


void Vector::print() const
{
    gsInfo << "Processor No. " << gsMpi::init().worldRank() << "\n" << *get();    
}


}//namespace trilinos

}// namespace gismo
