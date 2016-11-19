/** @file Operator.cpp

    @brief Wrapper for Trilinos/Epetra operator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsMpi/gsMpi.h>
#include <gsCore/gsLinearAlgebra.h>
#include "gsTrilinosHeaders.h"

#include <gsTrilinos/Operator.h>
#include <gsTrilinos/SparseMatrix.h>




namespace gismo
{

namespace trilinos
{

/**
   Implementation of Matrix operator
*/
class SparseMatrixEPetraOp : public Epetra_Operator
{
    typedef Teuchos::RCP<Epetra_CrsMatrix> EpetraMatrixPtr;
public:

    SparseMatrixEPetraOp() { }
    
    SparseMatrixEPetraOp(const EpetraMatrixPtr & mptr) : matrix(mptr) { }
    
    /* Methods required to support the Epetra_Operator interface */
    
    //! Returns a character string describing the operator
    const char* Label() const {return("trilinos::Operator");}

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
      does not support transpose use, this method should return a value of -1.

      \param UseTranspose - (In) If true, multiply by the transpose of operator, otherwise just use operator.

      \return Always returns 0.
    */
    int SetUseTranspose(bool UseTranspose_in)
    { return matrix->SetUseTranspose(UseTranspose_in); }

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*!
      \param X - (In) An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
      \param Y -(Out) An Epetra_MultiVector of dimension NumVectors containing result.

      \return Integer error code, set to 0 if successful.
      \pre Filled()==true
      \post Unchanged.
    */
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { return matrix->Apply(X,Y); }

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
      method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(),
      the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
      if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
      be accounted for when doing a triangular solve.

      \param X - (In) An Epetra_MultiVector of dimension NumVectors to solve for.
      \param Y - (Out) An Epetra_MultiVector of dimension NumVectors containing result.

      \return Integer error code, set to 0 if successful.
      \pre Filled()==true
      \post Unchanged.
    */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { return matrix->ApplyInverse(X,Y); }

    //! Returns true because this class can compute an Inf-norm.
    bool HasNormInf() const { return matrix->HasNormInf(); }

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const     { return matrix->UseTranspose(); }

    //! Returns the Epetra_Map object associated with the domain of this matrix operator.
    const Epetra_Map& OperatorDomainMap() const
    { return matrix->OperatorDomainMap(); }

    //! Returns the Epetra_Map object associated with the range of this matrix operator.
    const Epetra_Map& OperatorRangeMap() const
    { return matrix->OperatorRangeMap(); }

    double NormInf() const { return matrix->NormInf(); }

    const Epetra_Comm & Comm() const { return matrix->Comm(); }

protected:
     EpetraMatrixPtr matrix;
};

    
class OperatorPrivate
{
    friend class Operator;

    Teuchos::RCP<Epetra_Operator> op;
};
    
Operator::Operator(const gsSparseMatrix<> & sp)
: my(new OperatorPrivate)
{
    my->op.reset( new SparseMatrixEPetraOp( SparseMatrix(sp).getRCP() ) );
}

Operator::~Operator() { delete my; }

Epetra_Operator * Operator::get() const { return my->op.get(); }

void Operator::print() const
{
    gsMpiComm mpi_helper = gsMpi::init().worldComm();
    gsInfo << "Processor No. " << mpi_helper.rank() << "operator \n";    
}

}//namespace trilinos

}// namespace gismo
