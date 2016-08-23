/** @file gsTrilinosSolvers.h

    @brief Wrappers for Trilinos solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

namespace gismo
{

namespace trilinos
{

namespace solver
{


class GISMO_EXPORT AbstractSolver
{
public:
    
    AbstractSolver( Epetra_CrsMatrix *A )
    {
        Problem.SetOperator(A);
    }
    
    Epetra_Vector* solve( Epetra_Vector *b )
    {
        Problem.SetRHS(b);
        Epetra_Vector* x;
        Problem.SetLHS(x);
        solveProblem();
        
        Epetra_MultiVector MV= *Problem.GetLHS();
        return MV(0);
    }


protected:
        virtual void solveProblem() = 0;
    
protected:
    Epetra_LinearProblem Problem;
};








};// namespace solver
};// namespace trilinos
};// namespace gismo


