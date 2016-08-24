

#include "Vector.h"

#include "gsTrilinosHeaders.h"
#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsMemory.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>


namespace gismo
{

namespace trilinos
{

class VectorPrivate
{
    friend class Vector;

    VectorPrivate() { }

    VectorPrivate(Epetra_Vector * v_ptr)
    : vec(v_ptr, null_deleter<Epetra_Vector>)
    { }
            
    VectorPrivate(const Epetra_BlockMap & Map)
    : vec(new Epetra_Vector(Map))
    {
        // set to zero ?
        //vec->PutScalar(0.0);
    }

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
: my(new VectorPrivate(_map.get()->OperatorDomainMap()))
{ }

Vector::Vector(const gsVector<> & gsVec, const SparseMatrix & _map)
: my(new VectorPrivate(_map.get()->OperatorDomainMap()))
{
    // todo: fill in from gsVec to my->vec
}

Vector::Vector(Epetra_Vector * v_ptr) : my(new VectorPrivate(v_ptr))
{ }

Vector::Vector() : my(new VectorPrivate) { }

Vector::~Vector() { delete my; }

void Vector::setConstant(const double val)
{
    my->vec->PutScalar(val);
}

void Vector::setFrom(const SparseMatrix & A)
{
    my->vec.reset( new Epetra_Vector(A.get()->OperatorDomainMap()) );
}

size_t Vector::size() const 
{ 
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
    return (my->vec ? my->vec->GlobalLength64() : 0 ); 
#else
    return (my->vec ? my->vec->GlobalLength()   : 0 ); 
#endif
}

void Vector::copyTo(gsVector<real_t> & gsVec)
{
    gsVec.resize(size());
    my->vec->ExtractCopy(gsVec.data());
}

Epetra_Vector * Vector::get() const
{
    return my->vec.get();
}


}//namespace trilinos

}// namespace gismo
