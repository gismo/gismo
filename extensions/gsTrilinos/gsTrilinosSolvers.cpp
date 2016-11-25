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

//#include "Ifpack_ConfigDefs.h" 
//#include "Ifpack.h" 
//#include "Ifpack_AdditiveSchwarz.h" 

//#include "Amesos_Superlu.h"

#include "gsTrilinosHeaders.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"

namespace gismo
{

namespace trilinos
{

struct DataTypes // General types for use in all solvers
{
    typedef real_t Scalar;
                                 
    typedef conditional<util::is_same<Scalar,double>::value, Epetra_MultiVector,
                        Tpetra::MultiVector<Scalar,int,int> >::type MVector;
    
    typedef conditional<util::is_same<Scalar,double>::value,Epetra_Operator,
                        Tpetra::Operator<Scalar,int,int> >::type    Operator;

    typedef Belos::SolverManager<Scalar, MVector, Operator> SolManager;

    typedef Belos::LinearProblem<Scalar, MVector, Operator> BelosLp;

    typedef Epetra_LinearProblem EpetraLp;

    //typedef Teuchos::ParameterList
};

namespace solver
{

struct AbstractSolverPrivate
{
    // Problem
    DataTypes::EpetraLp Problem;

    // Solution vector
    Vector solution;
};

AbstractSolver::AbstractSolver() : my(NULL)
{

}


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

const Vector & AbstractSolver::solve(const Vector & b)
{
    my->Problem.SetRHS(b.get());
    solveProblem(); // virtual call
    return my->solution; 
}

void AbstractSolver::getSolution( gsVector<real_t> & sol, const int rank) const
{
    my->solution.copyTo(sol, rank);
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



/*   --- Belos--- */

// mode values define different solvers
template<int mode> struct BelosSolManager { };

template<>
struct BelosSolManager<BlockCG> 
{ typedef Belos::BlockCGSolMgr<DataTypes::Scalar , DataTypes::MVector, DataTypes::Operator> type; };

template<>
struct BelosSolManager<BlockGmres> 
{ typedef Belos::BlockGmresSolMgr<DataTypes::Scalar , DataTypes::MVector, DataTypes::Operator> type; };

struct BelosSolverPrivate
{
    Teuchos::RCP<DataTypes::SolManager> Solver;

    Teuchos::ParameterList belosList;
    
    DataTypes::BelosLp Problem;
};

template<int mode>
BelosSolver<mode>::BelosSolver(const SparseMatrix & A
                         //, const std::string solver_teuchosUser
    )
: Base(A), myBelos(new BelosSolverPrivate), blocksize(1), maxiters(500)
{ 
    // Note: By default the string variable SolverTeuchosUser = "Belos". 
    // This can be adapted to get different values if other solvers with 
    // similar implementation structures are considered later on. 
    //SolverTeuchosUser = solver_teuchosUser; 

    // Initialize solver manager
    myBelos->Solver = Teuchos::rcp( new typename BelosSolManager<mode>::type );

    myBelos->Problem.setOperator(A.getRCP());
    myBelos->Problem.setLHS(my->solution.getRCP());
    
	// If the matrix is symmetric, specify this in the linear problem. 
	// my->Problem->setHermitian(); 
}

template<int mode>
BelosSolver<mode>::~BelosSolver()
{
    delete myBelos;
}

template<int mode>
void BelosSolver<mode>::solveProblem()
{
    // Grab right-hand side
    myBelos->Problem.setRHS( Teuchos::rcp(my->Problem.GetRHS(), false) );

    // Tell the program that setting of the linear problem is done. 
    // Throw an error if failed. 
    bool err_set = myBelos->Problem.setProblem(); 
    
    GISMO_ASSERT(true == err_set, "Error: Belos Problem couldn't be" 
                 " initialized."); 

    // Teuchos::ScalarTraits<double>::magnitudeType tol = 1.0e-5;

    double tol = 1.0e-5;
    myBelos->belosList.set( "Block Size", blocksize );         // Blocksize to be used by iterative solver
    myBelos->belosList.set( "Maximum Iterations", maxiters );  // Maximum number of iterations allowed
    myBelos->belosList.set( "Convergence Tolerance", tol );    // Relative convergence tolerance requested

    // Set the solver manager data.
    myBelos->Solver->setProblem   (Teuchos::rcp(&myBelos->Problem  , false));
    myBelos->Solver->setParameters(Teuchos::rcp(&myBelos->belosList, false));

    // Perform solve
    Belos::ReturnType ret = myBelos->Solver->solve();
    
    // Get the number of iterations for this solve.
    
//  int numIters = myBelos->Solver->getNumIters();
    
//  // Compute actual residuals.
    
//  bool badRes = false;
//  std::vector<double> actual_resids( numrhs );
//  std::vector<double> rhs_norm( numrhs );
//  MVector resid(*Map, numrhs);
    
    GISMO_ENSURE(ret == Belos::Converged , "Error: Belos Problem couldn't be"
                 " initialized.");
}

template<int mode>
void BelosSolver<mode>::setBlockSize(int bs)
{
    myBelos->belosList.set( "Block Size", blocksize );
}

//------------------------------------------

CLASS_TEMPLATE_INST BelosSolver<BlockGmres>;
CLASS_TEMPLATE_INST BelosSolver<BlockCG>;

};// namespace solver
};// namespace trilinos
};// namespace gismo

