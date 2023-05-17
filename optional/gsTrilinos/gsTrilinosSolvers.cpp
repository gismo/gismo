/** @file gsTrilinosSolvers.cpp

    @brief Wrappers for Trilinos solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, F. Khatami, M. Moeller
*/

#include <gsTrilinos/gsTrilinosSolvers.h>
#include "gsTrilinosHeaders.h"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#include "AztecOO.h"
#include "AztecOO_Version.h"

#include "Ifpack_ConfigDefs.h" 
#include "Ifpack.h" 
#include "Ifpack_AdditiveSchwarz.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"

#include "BelosBiCGStabSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosFixedPointSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosGmresPolySolMgr.hpp"
#include "BelosLSQRSolMgr.hpp"
#include "BelosMinresSolMgr.hpp"
#include "BelosPCPGSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockStochasticCGSolMgr.hpp"
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosTFQMRSolMgr.hpp"
#ifdef Belos_ENABLE_Experimental
#include "BelosBlockGCRODRSolMgr.hpp"
#endif

namespace gismo
{

namespace trilinos
{

/// General types for use in all solvers
struct DataTypes 
{
    /// Data and index types
    typedef real_t  Scalar;
    typedef index_t Index;

    /// Multi vector type
    typedef util::conditional<util::is_same<Scalar,double>::value,
                              Epetra_MultiVector,
                              Tpetra::MultiVector<Scalar,Index,Index> >::type MVector;

    /// Operator type
    typedef util::conditional<util::is_same<Scalar,double>::value,
                              Epetra_Operator,
                              Tpetra::Operator<Scalar,Index,Index> >::type    Operator;

    /// Belos solver manager type
    typedef Belos::SolverManager<Scalar, MVector, Operator>                   SolManager;

    /// Belos linear problem type
    typedef Belos::LinearProblem<Scalar, MVector, Operator>                   BelosLp;

    /// Epetra linear problem type
    typedef Epetra_LinearProblem                                              EpetraLp;
};

namespace solver
{

/*    --- Abstract solver ---    */

struct AbstractSolverPrivate
{
    // Problem
    DataTypes::EpetraLp Problem;

    // Solution vector
    Vector Solution;
};

/// Constructor (default)
AbstractSolver::AbstractSolver() : my(NULL)
{ }

/// Constructor (sparse matrix)
AbstractSolver::AbstractSolver(const SparseMatrix & A)
: my(new AbstractSolverPrivate)
{
    my->Problem.SetOperator(A.get());
    my->Solution.setFrom(A); // i.e. A.get()->OperatorDomainMap()
    my->Problem.SetLHS(my->Solution.get());
}

/// Destructor
AbstractSolver::~AbstractSolver()
{
    delete my;
}

/// Solves problem for the given a right-hand side vector
const Vector & AbstractSolver::solve(const Vector & b)
{
    my->Problem.SetRHS(b.get());
    solveProblem(); // virtual call
    return my->Solution; 
}

/// Returns solution vector
void AbstractSolver::getSolution(gsVector<real_t> & sol, const int rank) const
{
    my->Solution.copyTo(sol, rank);
}

/// Sets parameters from option list
void AbstractSolver::setOptions(const gsOptionList & opt)
{
    std::vector<gsOptionList::OptionListEntry> options = opt.getAllEntries();
    for(std::vector<gsOptionList::OptionListEntry>::iterator it = options.begin();
        it != options.end(); ++it) {        
        if (it->type == "string") {
            set(it->label, it->val);
        }
        else if (it->type == "int") {
            set(it->label, atoi(it->val.c_str()));
        }
        else if (it->type == "real") {
            set(it->label, atof(it->val.c_str()));
            
        }
        else if (it->type == "bool") {
            set(it->label, (it->val != "0"));
        }
        else {
            GISMO_ERROR("Error : Invalid parameter type.");
        }
    }
}

/*    --- Abstract direct solver ---    */

/// No difference to AbstractSolver

/*    --- Abstract iterative solver ---    */

/// Constructor (default)
AbstractIterativeSolver::AbstractIterativeSolver()
: AbstractSolver::AbstractSolver(), tolerance(1e-6), maxIter(100)
{ }

/// Constructor (sparse matrix)
AbstractIterativeSolver::AbstractIterativeSolver(const SparseMatrix & A)
: AbstractSolver::AbstractSolver(A), tolerance(1e-6), maxIter(100)
{ }

/*    --- Amesos solver ---    */

struct AmesosSolverPrivate
{
    // Amesos solver
    Amesos_BaseSolver * Solver;

    // Parameter list
    Teuchos::ParameterList AmesosList;

    // Solver status
    int Status;
    
    // Constructor
    AmesosSolverPrivate(int solver)
    : solver(solver)
    { }

    // Destructor
    ~AmesosSolverPrivate()
    { }

    // Returns SolverType as std::string
    const std::string SolverType() const
    {
        switch(solver) {
        case AmesosSolvers::Lapack :
            return "Amesos_Lapack";
        case AmesosSolvers::KLU :
            return "Amesos_Klu";
        case AmesosSolvers::Umfpack :
            return "Amesos_Umfpack";
        case AmesosSolvers::Pardiso :
            return "Amesos_Pardiso";
        case AmesosSolvers::Taucs :
            return "Amesos_Taucs";
        case AmesosSolvers::SuperLU :
            return "Amesos_Superlu";
        case AmesosSolvers::SuperLUDist :
            return "Amesos_Superludist";
        case AmesosSolvers::Mumps :
            return "Amesos_Mumps";
        case AmesosSolvers::Dscpack :
            return "Amesos_Superlu";
        default :
            GISMO_ERROR("Error : Invalid Amesos solver");
        }
    }

    // Returns SolverList sublist
    Teuchos::ParameterList & SolverList()
    {
        switch(solver) {
        case AmesosSolvers::Lapack :
            return AmesosList.sublist("Lapack");
        case AmesosSolvers::KLU :
            return AmesosList.sublist("Klu");
        case AmesosSolvers::Umfpack :
            return AmesosList.sublist("Umfpack");
        case AmesosSolvers::Pardiso :
            return AmesosList.sublist("Pardiso");
        case AmesosSolvers::Taucs :
            return AmesosList.sublist("Taucs");
        case AmesosSolvers::SuperLU :
            return AmesosList.sublist("Superlu");
        case AmesosSolvers::SuperLUDist :
            return AmesosList.sublist("Superludist");
        case AmesosSolvers::Mumps :
            return AmesosList.sublist("Mumps");
        case AmesosSolvers::Dscpack :
            return AmesosList.sublist("Superlu");
        default :
            GISMO_ERROR("Error : Invalid Amesos solver");
        }
    }
    
private:

    // Solver type
    const int solver;
};

/// Constructor (sparse matrix)
AmesosSolver::AmesosSolver(const SparseMatrix & A, const int solver)
: Base(A), myAmesos(new AmesosSolverPrivate(solver))
{
    Amesos Factory;
    const std::string SolverType = myAmesos->SolverType();

    // Check availability of solver
    GISMO_ENSURE(Factory.Query(SolverType.c_str()),
                 "Error: Amesos solver " + SolverType + " is not available");

    // Create solver
    myAmesos->Solver = Factory.Create(SolverType.c_str(), my->Problem);
    
    // Initialize parameter list by default values
    myAmesos->Solver->setParameterList(Teuchos::RCP<Teuchos::ParameterList>(
                                           &myAmesos->AmesosList, false));
}

/// Destructor
AmesosSolver::~AmesosSolver()
{
    delete myAmesos;
}

/// Solves problem
void AmesosSolver::solveProblem()
{
    myAmesos->Status = myAmesos->Solver->SymbolicFactorization();
    if (myAmesos->Status != 0)
    {
        gsWarn << "Error: Amesos solver failed in symbolic factorization.\n";
        return;
    }

    myAmesos->Status = myAmesos->Solver->NumericFactorization();
    if (myAmesos->Status != 0)
    {
        gsWarn << "Error: Amesos solver failed in numeric factorization.\n";
        return;
    }

    myAmesos->Status = myAmesos->Solver->Solve();
    if (myAmesos->Status != 0)
    {
        gsWarn << "Error: Amesos solver failed.\n";
    }
}

/// Solves problem (overwrite default behaviour of factorization)
void AmesosSolver::solveProblem(const bool noSymbolicFactorization,
                                const bool noNumericFactorization)
{
    if (!noSymbolicFactorization)
        if (myAmesos->Status != 0)
        {
            gsWarn << "Error: Amesos solver failed in symbolic factorization.\n";
            return;
        }
    
    if (!noNumericFactorization)
        if (myAmesos->Status != 0)
        {
            gsWarn << "Error: Amesos solver failed in numeric factorization.\n";
            return;
        }

    myAmesos->Status = myAmesos->Solver->Solve();
    if (myAmesos->Status != 0)
    {
        gsWarn << "Error: Amesos solver failed.\n";
    }
}

/// Returns valid parameters
std::string AmesosSolver::validParams() const
{
    // Create temporal Amesos solver
    Teuchos::ParameterList AmesosList;
    Amesos Factory;
    Amesos_BaseSolver * Solver = Factory.Create(myAmesos->SolverType().c_str(), my->Problem);
    
    // Set default values
    Solver->setParameterList(Teuchos::RCP<Teuchos::ParameterList>(
                                 &AmesosList, false));

    // Return valid parameters
    std::ostringstream os;
    os << "Valid parameters of the current "
        + myAmesos->SolverType() + " solver: \n"
       << AmesosList;
    return os.str();
}

/// Returns current parameters
std::string AmesosSolver::currentParams() const
{ 
    std::ostringstream os;
    os << "Current parameters of the current "
        + myAmesos->SolverType() + " solver: \n"
       << myAmesos->AmesosList;
    return os.str();
}

/// Sets integer parameters
void AmesosSolver::set(const std::string & name, const int & value)
{
    myAmesos->SolverList().set( name, value );
    myAmesos->Solver->setParameterList(Teuchos::RCP<Teuchos::ParameterList>(
                                           &myAmesos->AmesosList, false));
}

/// Sets bool parameters
void AmesosSolver::set(const std::string & name, const bool & value)
{
    myAmesos->SolverList().set( name, value );
    myAmesos->Solver->setParameterList(Teuchos::RCP<Teuchos::ParameterList>(
                                           &myAmesos->AmesosList, false));
}

/// Sets double parameters
void AmesosSolver::set(const std::string & name, const double & value)
{
    myAmesos->SolverList().set( name, value );
    myAmesos->Solver->setParameterList(Teuchos::RCP<Teuchos::ParameterList>(
                                           &myAmesos->AmesosList, false));
}

/// Sets string parameters
void AmesosSolver::set(const std::string & name, const std::string & value)
{
    myAmesos->SolverList().set( name, value );
    myAmesos->Solver->setParameterList(Teuchos::RCP<Teuchos::ParameterList>(
                                           &myAmesos->AmesosList, false));
}

/// Returns status of the solver
std::string AmesosSolver::status() const
{
    // Redirect the output of Amesos::PrintStatus() from std::cout to
    // an output stream buffer, convert to string and return it
    std::streambuf *old = std::cout.rdbuf();
    std::ostringstream os;

    std::cout.rdbuf(os.rdbuf());
    myAmesos->Solver->PrintStatus();
    std::cout.rdbuf(old);

    std::string s = os.str();
    return s;
}

/// Returns timing of the solver
std::string AmesosSolver::timing() const
{
    // Redirect the output of Amesos::PrintTiming() from std::cout to
    // an output stream buffer, convert to string and return it
    std::streambuf *old = std::cout.rdbuf();
    std::ostringstream os;

    std::cout.rdbuf(os.rdbuf());
    myAmesos->Solver->PrintTiming();
    std::cout.rdbuf(old);

    std::string s = os.str();
    return s;
}

/*    --- Aztec solver ---    */

struct AztecSolverPrivate
{
    // Aztec solver
    AztecOO Solver;
    
    // Returns option as std::string
    std::string AztecOptionName(const int option)
    {
        switch(option) {
        case(AZ_solver)            : return "AZ_solver";
        case(AZ_scaling)           : return "AZ_scaling";
        case(AZ_precond)           : return "AZ_precond";
        case(AZ_conv)              : return "AZ_conv";
        case(AZ_output)            : return "AZ_output";
        case(AZ_pre_calc)          : return "AZ_pre_calc";
        case(AZ_max_iter)          : return "AZ_max_iter";
        case(AZ_poly_ord)          : return "AZ_poly_ord";
        case(AZ_overlap)           : return "AZ_overlap";
        case(AZ_type_overlap)      : return "AZ_type_overlap";
        case(AZ_kspace)            : return "AZ_kspace";
        case(AZ_orthog)            : return "AZ_orthog";
        case(AZ_aux_vec)           : return "AZ_aux_vec";
        case(AZ_reorder)           : return "AZ_reorder";
        case(AZ_keep_info)         : return "AZ_keep_info";
        case(AZ_recursion_level)   : return "AZ_recursion_level";
        case(AZ_print_freq)        : return "AZ_print_freq";
        case(AZ_graph_fill)        : return "AZ_graph_fill";
        case(AZ_subdomain_solve)   : return "AZ_subdomain_solve";
        case(AZ_init_guess)        : return "AZ_init_guess";
        case(AZ_keep_kvecs)        : return "AZ_keep_kvecs";
        case(AZ_apply_kvecs)       : return "AZ_apply_kvecs";
        case(AZ_orth_kvecs)        : return "AZ_orth_kvecs";
        case(AZ_ignore_scaling)    : return "AZ_ignore_scaling";
        case(AZ_check_update_size) : return "AZ_check_update_size";
        case(AZ_extreme)           : return "AZ_extreme";
        case(AZ_diagnostics)       : return "AZ_diagnostics";
            
        default :
            GISMO_ERROR("Error : Invalid Aztec option");
        }
    }

    // Returns parameters as std::string
    std::string AztecParamName(const int param)
    {
        switch(param) {
        case(AZ_tol)              : return "AZ_tol";
        case(AZ_drop)             : return "AZ_drop";
        case(AZ_ilut_fill)        : return "AZ_ilut_fill";
        case(AZ_omega)            : return "AZ_omega";
        case(AZ_rthresh)          : return "AZ_rthresh";
        case(AZ_athresh)          : return "AZ_athresh";
        case(AZ_update_reduction) : return "AZ_update_reduction";
        case(AZ_temp)             : return "AZ_temp";
        case(AZ_ill_cond_thresh)  : return "AZ_ill_cond_thresh";
        case(AZ_weights)          : return "AZ_weights";

        default :
            GISMO_ERROR("Error : Invalid Aztec parameter");
        }
    }
};

/// Constructor (sparse matrix)
AztecSolver::AztecSolver( const SparseMatrix &A,
                          const int solver,
                          const int precond,
                          const int subdomain_solver )
: Base(A), myAztec(new AztecSolverPrivate)
{
    // Set linear problem (yet without RHS)
    myAztec->Solver.SetProblem(my->Problem);
    
    // Set solver
    set(AZ_solver,  solver);

    // Set preconditioner
    set(AZ_precond, precond);

    // Set subdomain solver
    set(AZ_subdomain_solve, subdomain_solver);
    
    // Set default parameters
    set(AZ_max_iter, maxIter);
    set(AZ_tol, tolerance);
}

/// Destructor
AztecSolver::~AztecSolver()
{ 
    delete myAztec;
}

/// Sets Belos preconditioner
void AztecSolver::setPreconditioner(const BelosSolver & Belos)
{
    myAztec->Solver.SetPrecOperator(Belos.getPrecOperator());
}

/// Sets ML preconditioner
void AztecSolver::setPreconditioner(const MLSolver & ML)
{
    myAztec->Solver.SetPrecOperator(ML.getPrecOperator());
}

/// Solves the problem
void AztecSolver::solveProblem()
{
    // Grab right-hand side
    myAztec->Solver.SetRHS(my->Problem.GetRHS());

    // Get configuration from internal arrays
    const int *options = myAztec->Solver.GetAllAztecOptions();
    const double *params = myAztec->Solver.GetAllAztecParams();
        
    // Solve linear problem
    myAztec->Solver.Iterate(options[AZ_max_iter], params[AZ_tol]);
}

/// Returns valid parameters
std::string AztecSolver::validParams() const
{
    // Create temporal Aztec solver with default options/parameters
    AztecOO Solver;
    
    const int *options = Solver.GetAllAztecOptions();
    const double *params = Solver.GetAllAztecParams();

    std::ostringstream os;
    os << "Valid parameters of the current Aztec solver: \n";
                                                                
    for (int AZ_def=0; AZ_def<AZ_FIRST_USER_OPTION; AZ_def++)
        os << myAztec->AztecOptionName(AZ_def) << " = " << options[AZ_def] << "\n";

    for (int AZ_def=0; AZ_def<AZ_FIRST_USER_PARAM; AZ_def++)
        os << myAztec->AztecParamName(AZ_def) << " = " << params[AZ_def] << "\n";
            
    return os.str();
}

/// Returns current parameters
std::string AztecSolver::currentParams() const
{
    const int *options = myAztec->Solver.GetAllAztecOptions();
    const double *params = myAztec->Solver.GetAllAztecParams();
    
    std::ostringstream os;
    os << "Current parameters of the current Aztec solver: \n";

    for (int AZ_def=0; AZ_def<AZ_FIRST_USER_OPTION; AZ_def++)
        os << myAztec->AztecOptionName(AZ_def) << " = " << options[AZ_def] << "\n";

    for (int AZ_def=0; AZ_def<AZ_FIRST_USER_PARAM; AZ_def++)
        os << myAztec->AztecParamName(AZ_def) << " = " << params[AZ_def] << "\n";

    return os.str();
}

/// Returns status of the solver
std::string AztecSolver::status() const
{
    const double *status = myAztec->Solver.GetAztecStatus();
    const int solver = myAztec->Solver.GetAztecOption(AZ_solver);
    
    std::string prefix = std::string("Aztec_")
        + std::string(solver == AZ_cg ? "CG" :
                      solver == AZ_cg_condnum ? "CG_condnum" :
                      solver == AZ_gmres ? "Gmres" :
                      solver == AZ_gmres_condnum ? "Gmres_condnum" :
                      solver == AZ_GMRESR ? "Gmres_recursive" :
                      solver == AZ_cgs ? "CGS" :
                      solver == AZ_tfqmr ? "TFQMR" :
                      solver == AZ_bicgstab ? "BiCGStab" :
                      solver == AZ_lu ? "LU" :
                      solver == AZ_slu ? "SuperLU" :
                      solver == AZ_symmlq ? "SymmLQ" :
                      solver == AZ_fixed_pt ? "Fixed_point" :
                      solver == AZ_analyze ? "Analyze" :
                      "unknown")
        +  std::string(" : ");

    std::ostringstream os;
    os << "----------------------------------------------------------------------------\n"
       << prefix << "Number if iterations AZ_its = " << status[AZ_its] << "\n"
       << prefix << "Iteration termination AZ_why = "
       << (status[AZ_why] == AZ_normal ? "normally" :
           status[AZ_why] == AZ_param ? "invalid option" :
           status[AZ_why] == AZ_breakdown ? "numerical breakdown" :
           status[AZ_why] == AZ_loss ? "precision loss" :
           status[AZ_why] == AZ_ill_cond ? "ill-conditioned Hessenberg matrix in GMRES" :
           status[AZ_why] == AZ_maxits ? "maximum number of iterations reached" :
           "unknown") << "\n"
       << prefix << "Absolute residual AZ_r = " << status[AZ_r] << "\n"
       << prefix << "True scaled residual AZ_scaled_r = " << status[AZ_scaled_r] << "\n"
       << prefix << "Estimated scaled residual AZ_rec_r = " << status[AZ_rec_r] << "\n"
       << "----------------------------------------------------------------------------\n";
    return os.str();
}

/// Returns timing of the solver
std::string AztecSolver::timing() const
{
    const double *status = myAztec->Solver.GetAztecStatus();
    const int solver = myAztec->Solver.GetAztecOption(AZ_solver);

    std::string prefix = std::string("Aztec_")
        + std::string(solver == AZ_cg ? "CG" :
                      solver == AZ_cg_condnum ? "CG_condnum" :
                      solver == AZ_gmres ? "Gmres" :
                      solver == AZ_gmres_condnum ? "Gmres_condnum" :
                      solver == AZ_GMRESR ? "Gmres_recursive" :
                      solver == AZ_cgs ? "CGS" :
                      solver == AZ_tfqmr ? "TFQMR" :
                      solver == AZ_bicgstab ? "BiCGStab" :
                      solver == AZ_lu ? "LU" :
                      solver == AZ_slu ? "SuperLU" :
                      solver == AZ_symmlq ? "SymmLQ" :
                      solver == AZ_fixed_pt ? "Fixed_point" :
                      solver == AZ_analyze ? "Analyze" :
                      "unknown")
        +  std::string(" : ");

    std::ostringstream os;
    os << "----------------------------------------------------------------------------\n"
       << prefix << "Total time spent in Aztec = " << status[AZ_solve_time] << " (s)\n"
       << "----------------------------------------------------------------------------\n";
    return os.str();
}

/// Sets integer paramters
void AztecSolver::set(const std::string & name, const int & value)
{
    Teuchos::ParameterList AztecList;
    AztecList.set( name, value );
    myAztec->Solver.SetParameters(AztecList);
}

/// Sets bool parameters
void AztecSolver::set(const std::string & name, const bool & value)
{
    Teuchos::ParameterList AztecList;
    AztecList.set( name, value );
    myAztec->Solver.SetParameters(AztecList);
}

/// Sets double parameters
void AztecSolver::set(const std::string & name, const double & value)
{
    Teuchos::ParameterList AztecList;
    AztecList.set( name, value );
    myAztec->Solver.SetParameters(AztecList);
}

/// Sets string parameters
void AztecSolver::set(const std::string & name, const std::string & value)
{
    Teuchos::ParameterList AztecList;
    AztecList.set( name, value );
    myAztec->Solver.SetParameters(AztecList);
}

/// Sets Aztec option directly
void AztecSolver::set(const int & option, const int & value)
{
    myAztec->Solver.SetAztecOption( option, value );
}

/// Sets Aztec parameter directly
void AztecSolver::set(const int & param, const double & value)
{
    myAztec->Solver.SetAztecParam( param, value );
}

/// Returns number of iterations
int AztecSolver::numIterations() const
{
    return myAztec->Solver.NumIters();
}

/*   --- Belos solver ---    */

struct BelosSolverPrivate
{
    // Belos solver
    Teuchos::RCP<DataTypes::SolManager> Solver;

    // Parameter list
    Teuchos::ParameterList BelosList;
    
    // Belos problem    
    DataTypes::BelosLp Problem;

    // Belos status
    Belos::ReturnType Status;
    
    // !!!NOTE!!! : Belos solvers are configured by parameter
    // lists. However, the solver manager stores an internal copy of
    // all parameters so that it makes no sense to keep a separate
    // parameter list. The most efficient way to set parameters is to
    // create a temporal parameter list when needed and pass only the
    // parameters needed.
};

/// Constructor (sparse matrix)
BelosSolver::BelosSolver(const SparseMatrix & A,
                               const int solver)
: Base(A), myBelos(new BelosSolverPrivate)
{
    // Initialize solver manager
    switch(solver) {
    case BelosSolvers::BiCGStab :
        myBelos->Solver = Teuchos::rcp( new typename Belos::BiCGStabSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::BlockCG :
        myBelos->Solver = Teuchos::rcp( new typename Belos::BlockCGSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
#ifdef Belos_ENABLE_Experimental
    case BelosSolvers::BlockGCRODR :
        myBelos->Solver = Teuchos::rcp( new typename Belos::BlockGCRODRSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
#endif
    case BelosSolvers::BlockGmres :
        myBelos->Solver = Teuchos::rcp( new typename Belos::BlockGmresSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::FixedPoint :
        myBelos->Solver = Teuchos::rcp( new typename Belos::FixedPointSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::GCRODR :
        myBelos->Solver = Teuchos::rcp( new typename Belos::GCRODRSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::GmresPoly :
        myBelos->Solver = Teuchos::rcp( new typename Belos::GmresPolySolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::LSQR :
        myBelos->Solver = Teuchos::rcp( new typename Belos::LSQRSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::Minres :
        myBelos->Solver = Teuchos::rcp( new typename Belos::MinresSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::PCPG :
        myBelos->Solver = Teuchos::rcp( new typename Belos::PCPGSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::PseudoBlockCG :
        myBelos->Solver = Teuchos::rcp( new typename Belos::PseudoBlockCGSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::PseudoBlockGmres :
        myBelos->Solver = Teuchos::rcp( new typename Belos::PseudoBlockGmresSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::PseudoBlockStochasticCG :
        myBelos->Solver = Teuchos::rcp( new typename Belos::PseudoBlockStochasticCGSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::PseudoBlockTFQMR :
        myBelos->Solver = Teuchos::rcp( new typename Belos::PseudoBlockTFQMRSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::RCG :
        myBelos->Solver = Teuchos::rcp( new typename Belos::RCGSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
    case BelosSolvers::TFQMR :
        myBelos->Solver = Teuchos::rcp( new typename Belos::TFQMRSolMgr<DataTypes::Scalar ,
                                        DataTypes::MVector, DataTypes::Operator> );
        break;
        
    default :
        GISMO_ERROR("Error : Invalid Belos solver");
    }

    // Initialize problem
    myBelos->Solver->setProblem(Teuchos::rcp(&myBelos->Problem  , false));
    myBelos->Problem.setOperator(A.getRCP());
    myBelos->Problem.setLHS(my->Solution.getRCP());
    
    // Set default parameters (via temporal parameter list)
    myBelos->BelosList.set( "Maximum Iterations", maxIter );
    myBelos->BelosList.set( "Convergence Tolerance", tolerance );
    myBelos->Solver->setParameters(Teuchos::RCP<Teuchos::ParameterList>(&myBelos->BelosList, false));
}

/// Destructor
BelosSolver::~BelosSolver()
{
    delete myBelos;
}

int BelosSolver::setPreconditioner(
    const std::string & PrecType, const SparseMatrix &A,
    const bool & leftprec, const int & OverlapLevel ) 
{
    // allocates an IFPACK factory. No data is associated
    // to this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation

    Teuchos::RCP<Epetra_RowMatrix> AA = A.getRCP();

    Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create
                                             (PrecType, &*AA, OverlapLevel) );

    assert(Prec != Teuchos::null);

//    // specify parameters for ICT
//    myBelos->belosList.set("fact: drop tolerance", 1e-9);
//    myBelos->belosList.set("fact: ict level-of-fill", 1.0);
//
//    // the combine mode is on the following:
//    // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
//    // Their meaning is as defined in file Epetra_CombineMode.h
//    myBelos->belosList.set("schwarz: combine mode", "Add");

    // sets the parameters
    //IFPACK_CHK_ERR(Prec->SetParameters(myBelos->belosList));

    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(Prec->Initialize());

    // Builds the preconditioners, by looking for the values of
    // the matrix.
    IFPACK_CHK_ERR(Prec->Compute());

    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = 
                              Teuchos::rcp( new Belos::EpetraPrecOp( Prec ) );

    if (leftprec)
    {
      myBelos->Problem.setLeftPrec( belosPrec );
    }
    else 
    {
      myBelos->Problem.setRightPrec( belosPrec );
    }

    return 0;
}

/// Solves problem
void BelosSolver::solveProblem()
{
    // Grab right-hand side
    myBelos->Problem.setRHS( Teuchos::rcp(my->Problem.GetRHS(), false) );

    // Set problem
    GISMO_ASSERT(myBelos->Problem.setProblem(),
                 "Error: Belos solver failed in initialization."); 

    // Perform solve
    myBelos->Status = myBelos->Solver->solve();
}

/// Returns valid parameters
std::string BelosSolver::validParams() const
{
    std::ostringstream os;
    os << "Valid parameters of the current Belos solver: \n" 
       << *myBelos->Solver->getValidParameters() << "\n";
    return os.str();
}

/// Returns current parameters
std::string BelosSolver::currentParams() const
{
    std::ostringstream os;
    os << "Current parameters of the current Belos solver: \n" 
       << *myBelos->Solver->getCurrentParameters() << "\n";
    return os.str();
}

/// Returns status of the solver
std::string BelosSolver::status() const
{
    std::ostringstream os;
    os << "----------------------------------------------------------------------------\n"
       << "Belos : Number of iterations = " << myBelos->Solver->getNumIters() << "\n"
       << "Belos : Iteration converged = " << (myBelos->Status==Belos::Converged) << "\n"
       << "Belos : Absolute residual = " << myBelos->Solver->achievedTol() << "\n"
       << "Belos : Loss of accuracy detected = " << myBelos->Solver->isLOADetected() << "\n"
       << "----------------------------------------------------------------------------\n";
    return os.str();
}

/// Returns timing of the solver
std::string BelosSolver::timing() const
{
    return "Timing for Belos solver not yet implemented\n";
}

/// Sets integer paramters
void BelosSolver::set(const std::string & name, const int & value)
{
    myBelos->BelosList.set( name, value );
    myBelos->Solver->setParameters(Teuchos::RCP<Teuchos::ParameterList>(&myBelos->BelosList, false));
}

/// Sets bool paramters
void BelosSolver::set(const std::string & name, const bool & value)
{
    myBelos->BelosList.set( name, value );
    myBelos->Solver->setParameters(Teuchos::RCP<Teuchos::ParameterList>(&myBelos->BelosList, false));
}

/// Sets double paramters
void BelosSolver::set(const std::string & name, const double & value)
{
    myBelos->BelosList.set( name, value );
    myBelos->Solver->setParameters(Teuchos::RCP<Teuchos::ParameterList>(&myBelos->BelosList, false));
}

/// Sets string paramters
void BelosSolver::set(const std::string & name, const std::string & value)
{
    myBelos->BelosList.set( name, value );
    myBelos->Solver->setParameters(Teuchos::RCP<Teuchos::ParameterList>(&myBelos->BelosList, false));
}

/// Sets Hermitian problem type
void BelosSolver::setHermitian()
{
    myBelos->Problem.setHermitian();
}

/// Returns number of iterations
int BelosSolver::numIterations() const
{
    return myBelos->Solver->getNumIters();
}

/// Returns pointer to internal preconditioner
Epetra_Operator * BelosSolver::getPrecOperator() const
{
    //TODO   return myBelos->Problem.getOperator().get();
    GISMO_NO_IMPLEMENTATION
}

/*    --- Multi Level (ML) ---    */

struct MLSolverPrivate
{
    // Multi-level preconditioner
    ML_Epetra::MultiLevelPreconditioner* MLPrec;

    // Aztec solver
    AztecOO Solver;

    // Parameter list
    Teuchos::ParameterList MLList;
};

/// Constructor (sparse matrix)
MLSolver::MLSolver( const SparseMatrix &A,
                    const int solver )
: Base(A), myML(new MLSolverPrivate)
{
    // Initialize parameter list
    switch(solver) {
    case MLSolvers::SA :
        ML_Epetra::SetDefaults("SA", myML->MLList);
        break;
    case MLSolvers::NSSA :
        ML_Epetra::SetDefaults("NSSA", myML->MLList);
        break;
    case MLSolvers::DD :
        ML_Epetra::SetDefaults("DD", myML->MLList);
        break;
    case MLSolvers::DDLU :
        ML_Epetra::SetDefaults("DD-LU", myML->MLList);
        break;
    case MLSolvers::DDML :
        ML_Epetra::SetDefaults("DD-ML", myML->MLList);
        break;
    case MLSolvers::DDMLLU :
        ML_Epetra::SetDefaults("DD-ML-LU", myML->MLList);
        break;
    default :
        GISMO_ERROR("Error : Invalid ML solver");
    }

    // Attach preconditioner but do not compute it yet
    Teuchos::RCP<Epetra_RowMatrix> AA = A.getRCP();
    myML->MLPrec = new ML_Epetra::MultiLevelPreconditioner
                                     (*AA.get(), myML->MLList, false);
}

/// Destructor
MLSolver::~MLSolver()
{ 
    delete myML;
}

/// Solves problem
void MLSolver::solveProblem()
{
    myML->Solver.SetProblem(my->Problem);

    if (!myML->MLPrec->IsPreconditionerComputed())
        myML->MLPrec->ComputePreconditioner();
    
    myML->Solver.SetPrecOperator(myML->MLPrec);

    // Get configuration from internal arrays
    const int *options = myML->Solver.GetAllAztecOptions();
    const double *params = myML->Solver.GetAllAztecParams();
        
    // Solve linear problem
    myML->Solver.Iterate(options[AZ_max_iter], params[AZ_tol]);
}

/// Returns valid parameters
std::string MLSolver::validParams() const
{
    // Temporal parameter list
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA", MLList);
    
    std::ostringstream os;
    os << "Valid parameters of the current ML solver: \n";
    MLList.print(os);
    return os.str();
}

/// Returns current parameters
std::string MLSolver::currentParams() const
{
    std::ostringstream os;
    os << "Current parameters of the current ML solver: \n";
    myML->MLList.print(os);
    return os.str();
}

/// Returns status of the solver
std::string MLSolver::status() const
{
    myML->MLPrec->ReportTime();
    return "";
}

/// Returns timing of the solver
std::string MLSolver::timing() const
{
    return "";
}

/// Sets integer paramters
void MLSolver::set(const int & option, const int & value)
{
    myML->Solver.SetAztecOption(option, value);
}

/// Sets integer parameters
void MLSolver::set(const std::string & name, const int & value)
{
    myML->MLList.set( name, value );
    myML->MLPrec->SetParameterList(myML->MLList);
}

/// Sets bool parameters
void MLSolver::set(const std::string & name, const bool & value)
{
    myML->MLList.set( name, value );
    myML->MLPrec->SetParameterList(myML->MLList);
}

/// Sets double parameters
void MLSolver::set(const std::string & name, const double & value)
{
    myML->MLList.set( name, value );
    myML->MLPrec->SetParameterList(myML->MLList);
}

/// Sets string parameters
void MLSolver::set(const std::string & name, const std::string & value)
{
    myML->MLList.set( name, value );
    myML->MLPrec->SetParameterList(myML->MLList);
}

/// Returns number of iterations
int MLSolver::numIterations() const
{
    return 0;
}

/// Returns pointer to internal preconditioner
Epetra_Operator * MLSolver::getPrecOperator() const
{
    if (!myML->MLPrec->IsPreconditionerComputed())
         myML->MLPrec->ComputePreconditioner();
    
    return myML->MLPrec;
}
    
//------------------------------------------

};// namespace solver
};// namespace trilinos
};// namespace gismo

