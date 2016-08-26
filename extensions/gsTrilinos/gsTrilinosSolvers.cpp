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

//#include "Ifpack_ConfigDefs.h"
//#include "Ifpack.h"
//#include "Ifpack_AdditiveSchwarz.h"

#include "AztecOO.h"
#include "AztecOO_Version.h"

//#include "Amesos_Superlu.h"

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
    /* 1. Preconditionner */
/*    
    // allocates an IFPACK factory. No data is associated 
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    std::string PrecType = "ILU"; // incomplete LU
    int OverlapLevel = 1; // must be >= 0. ignored for Comm.NumProc() == 1
    Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType, &*A, OverlapLevel) );
    assert(Prec != Teuchos::null);
    
    // specify parameters for ILU
    List.set("fact: drop tolerance", 1e-9);
    List.set("fact: level-of-fill", 1);
    // the combine mode is on the following:
    // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
    // Their meaning is as defined in file Epetra_CombineMode.h   
    List.set("schwarz: combine mode", "Add");
    // sets the parameters
    IFPACK_CHK_ERR(Prec->SetParameters(List));
    
    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(Prec->Initialize());
    
    // Builds the preconditioners, by looking for the values of 
    // the matrix.
    IFPACK_CHK_ERR(Prec->Compute());
    
*/
    /* 2. AztecOO solver / GMRES*/
    
    AztecOO Solver;
    Solver.SetProblem(my->Problem);
    Solver.SetAztecOption(AZ_solver, AZ_gmres);
    Solver.SetAztecOption(AZ_output,32);
    //Solver.SetPrecOperator(Prec);
    Solver.Iterate(m_maxIter, m_tolerance);
}

void KLU::solveProblem()
{
    static const char * SolverType = "Klu";
    static Amesos Factory;
    const bool IsAvailable = Factory.Query(SolverType);
    GISMO_ENSURE(IsAvailable, "Amesos KLU is not available.\n");
    
    Amesos_BaseSolver * Solver = Factory.Create(SolverType, my->Problem);
    Solver->SymbolicFactorization();
    Solver->NumericFactorization();
    Solver->Solve();
}

void SuperLU::solveProblem()
{
//    Amesos_Superlu Solver(my->Problem);
//    Solver.SymbolicFactorization();
//    Solver.NumericFactorization();
//    Solver.Solve();
}



};// namespace solver
};// namespace trilinos
};// namespace gismo

