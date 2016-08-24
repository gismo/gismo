/** @file gsTrilinosSolvers.cpp

    @brief Wrappers for Trilinos solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsTrilinos/gsTrilinosSolvers.h>
#include "gsTrilinosHeaders.h"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#include "AztecOO.h"
#include "AztecOO_Version.h"

#include <Epetra_LinearProblem.h>

namespace gismo
{

namespace trilinos
{

namespace solver
{

struct AbstractSolverPrivate
{
    Epetra_LinearProblem Problem;

    Vector solution;
};


AbstractSolver::AbstractSolver(const SparseMatrix & A)
: my(new AbstractSolverPrivate)
{
    my->Problem.SetOperator(A.get());
    my->solution.setFrom(A); // i.e. A.get()->OperatorDomainMap()
    my->Problem.SetLHS(my->solution.get());
}

AbstractSolver::~AbstractSolver()
{
    delete my;
}

const Vector & AbstractSolver::solve( const Vector & b )
{
    my->Problem.SetRHS(b.get());
    solveProblem(); // virtual call
    
    //Epetra_MultiVector & MV = *my->Problem.GetLHS();
    //return Vector(MV(0));
    
    return my->solution;
}

void AbstractSolver::getSolution( gsVector<> & sol ) const
{
    my->solution.copyTo(sol);
}

void GMRES::solveProblem()
{
    AztecOO Solver;
    Solver.SetProblem(my->Problem);
    Solver.SetAztecOption(AZ_solver, AZ_gmres);
    Solver.SetAztecOption(AZ_output,32);
    //Solver.SetPrecOperator(Prec);
    Solver.Iterate(m_maxIter, m_tolerance);
}

void KLU::solveProblem()
{
    static const char * SolverType = "Amesos_KLU";
    static Amesos Factory;
    const bool IsAvailable = Factory.Query(SolverType);
    GISMO_ENSURE(IsAvailable, "Amesos KLU is not available.\n");
    
    Amesos_BaseSolver * Solver = Factory.Create(SolverType, my->Problem);
    Solver->SymbolicFactorization();
    Solver->NumericFactorization();
    Solver->Solve();
}




};// namespace solver
};// namespace trilinos
};// namespace gismo

