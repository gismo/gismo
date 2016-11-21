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

//#include "Amesos_Superlu.h"

#include "gsTrilinosHeaders.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"

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

class AbstractSolverBelosPrivate
{
public:

    AbstractSolverBelosPrivate() {}
    AbstractSolverBelosPrivate(const SparseMatrix & A)
    {
      solution.setFrom(A);
    }
    ~AbstractSolverBelosPrivate() {}

    Teuchos::RCP< Belos::LinearProblem<double,Epetra_MultiVector,
                                              Epetra_Operator> > Problem;
    Vector solution;
};


AbstractSolver::AbstractSolver(const SparseMatrix & A, const Vector & b)
: my(new AbstractSolverPrivate)
{
    my->Problem.SetOperator(A.get());
    my->solution.setFrom(A); // i.e. A.get()->OperatorDomainMap()
    my->Problem.SetLHS(my->solution.get());
    my->Problem.SetRHS(b.get());
}

AbstractSolver::AbstractSolver(const SparseMatrix & A, const Vector & b,
                               const std::string solverteuchosUser)
                               : myBelos(new AbstractSolverBelosPrivate(A))
{
    // Note: By default the string variable SolverTeuchosUser = "Belos".
    // This can be adapted to get different values if other solvers with
    // similar implementation structures are considered later on.

    SolverTeuchosUser = solverteuchosUser;

    myBelos->Problem = Teuchos::rcp(new Belos::LinearProblem
                       <double, Epetra_MultiVector, Epetra_Operator> 
                       (A.getRCP(), myBelos->solution.getRCP(), b.getRCP()));

//!!!!!!// If the matrix is symmetric, specify this in the linear problem.
//!!!!!!my->Problem->setHermitian();

      // Tell the program that setting of the linear problem is done.
      // Throw an error if failed.
      bool err_set = myBelos->Problem->setProblem();
   
      GISMO_ASSERT(true == err_set, "Error: Belos Problem couldn't be"
                   " initialized.");

}

AbstractSolver::~AbstractSolver()
{
    if ( my != NULL)
      delete my;

    if (myBelos != NULL)
      delete myBelos;
}

const Vector & AbstractSolver::solve()
{
    solveProblem(); // virtual call

    if( SolverTeuchosUser == "Belos" )
    {
      return myBelos->solution;
    }
    else if( SolverTeuchosUser == "")
    {
      //Epetra_MultiVector & MV = *my->Problem.GetLHS();
      //return Vector(MV(0));
      return my->solution;
    }
    else
    {
      GISMO_ERROR("Fatal error: something went wrong with solver.");
    }
}

void AbstractSolver::getSolution( gsVector<real_t> & sol, const int rank) const
{

    if( SolverTeuchosUser == "Belos" )
    {
      myBelos->solution.copyTo(sol,rank);
    }
    else if( SolverTeuchosUser == "")
    {
      my->solution.copyTo(sol,rank);
    }
    else
    {
      GISMO_ERROR("Fatal error: something went wrong with solver.");
    }
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
    Solver.SetAztecOption(AZ_output,AZ_none);//32
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

void Belos_solver::solveProblem()
{

    // Teuchos::ScalarTraits<double>::magnitudeType tol = 1.0e-5;

    double tol = 1.0e-5;

    Teuchos::ParameterList belosList;
    belosList.set( "Block Size", blocksize );         // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );  // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );    // Relative convergence tolerance requested

  // Create an iterative solver manager.

  Teuchos::RCP< Belos::SolverManager
              <double, Epetra_MultiVector, Epetra_Operator> > Solver = 
              Teuchos::rcp( new Belos::BlockCGSolMgr
              <double,Epetra_MultiVector,Epetra_Operator>
              ((myBelos->Problem), Teuchos::rcp(&belosList,false)) );

  // Perform solve

  Belos::ReturnType ret = Solver->solve();

  // Get the number of iterations for this solve.
  
//  int numIters = Solver->getNumIters();

//  // Compute actual residuals.

//  bool badRes = false;
//  std::vector<double> actual_resids( numrhs );
//  std::vector<double> rhs_norm( numrhs );
//  Epetra_MultiVector resid(*Map, numrhs);

    GISMO_ENSURE(ret == Belos::Converged , "Error: Belos Problem couldn't be"
                 " initialized.");
}


};// namespace solver
};// namespace trilinos
};// namespace gismo

