/** @file gsTrilinosEigenProblem.h

    @brief Wrappers for Trilinos parallel eigenvalue solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

// FD Trilinos
class Epetra_MultiVector;
class Epetra_Operator;

namespace gismo
{

namespace trilinos
{

namespace solver
{


enum AnasaziMethod { BlockDavidson = 1, LOBPCG = 2, BlockKrylovSchur = 3 };

class EigenProblemPrivate;

/**
    \brief Computes the eigenvalues of largest magnitude of an
    eigenvalue problem $A x = \lambda x$, using Anasazi's
    implementation of the Block Davidson method.

    \see https://trilinos.org/docs/dev/packages/anasazi/doc/html/classAnasazi_1_1SolverManager.html

*/
class GISMO_EXPORT EigenProblem
{
public:
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;

public:

    explicit EigenProblem(const SparseMatrix & A,
                          const AnasaziMethod & method = BlockKrylovSchur);
    // BlockDavidson

    ~EigenProblem();

    void solve() const;


protected:
        EigenProblemPrivate * my;
};


}// namespace solver
}// namespace trilinos
}// namespace gismo
