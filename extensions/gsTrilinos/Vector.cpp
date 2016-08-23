

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
            
    VectorPrivate(const Epetra_Map & Map)
    : vec(new Epetra_Vector(Map))
    {
        // set to one for now
        vec->PutScalar(1.0);
    }

/*   
    Epetra_SerialComm comm;

    /// Epetra Trilinos mapping of the matrix columns that assigns
    /// parts of the matrix to the individual processes.
    memory::shared_ptr<Epetra_Map> column_space_map;
*/
    
    /// A sparse matrix object in Trilinos 
    memory::shared_ptr<Epetra_Vector> vec;
};


Vector::Vector() : my(new VectorPrivate) { }

Vector::~Vector() { delete my; }

void Vector::copyTo(gsVector<real_t> & gsVec)
{
    //gsVec.resize(my->vec->size());
    
}


}//namespace trilinos

}// namespace gismo
