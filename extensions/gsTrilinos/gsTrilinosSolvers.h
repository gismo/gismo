/** @file gsTrilinosSolvers.h

    @brief Wrappers for Trilinos solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, F. Khatami, M. Moeller
 */

#pragma once

#include <gsTrilinos/Vector.h>
#include <gsIO/gsOptionList.h>

#include <az_aztec_defs.h> // in external

namespace gismo
{

namespace trilinos
{

/** @namespace gismo::trilinos::solver

    @brief This namespace contains wrappers for Trilinos linear
    system solvers and eigenvalue solvers
 */
namespace solver
{

/// Forward declaration
struct AmesosSolver;
struct AztecSolver;
struct BelosSolver;
struct MLSolver;

/// Amesos sparse direct solvers
struct AmesosSolvers
{
    enum {
        Lapack                   = 1,                ///< Lapack (serial)
        KLU                      = 2,                ///< KLU (serial)
        Umfpack                  = 3,                ///< Umfpack (serial)
        Pardiso                  = 4,                ///< Pardiso (serial/OpenMP)
        Taucs                    = 5,                ///< Taucs (serial)
        SuperLU                  = 6,                ///< SuperLU (serial)
        SuperLUDist              = 7,                ///< SuperLU (parallel)
        Mumps                    = 8,                ///< Mumps (parallel)
        Dscpack                  = 9                 ///< Dscpack (parallel)
    };
};

/// Aztec solvers
struct AztecSolvers
{
    enum {
        CG                       = AZ_cg,            ///< CG solver
        CGCondNum                = AZ_cg_condnum,    ///< CG solver with cond. num. estimation
        Gmres                    = AZ_gmres,         ///< GMRES solver
        GmresCondNum             = AZ_gmres_condnum, ///< GMRES solver with cond.num. estimation
        CGS                      = AZ_cgs,           ///< CG squared solver
        TFQMR                    = AZ_tfqmr,         ///< TFQMR solver
        BiCGStab                 = AZ_bicgstab,      ///< BiCGStab solver
        LU                       = AZ_lu             ///< sparse direct LU solver
    };
};

/// Aztec preconditioners
struct AztecPreconds
{
    enum {
        None                     = AZ_none,          ///< No preconditioning
        Jacobi                   = AZ_Jacobi,        ///< Jacobi
        Neumann                  = AZ_Neumann,       ///< Neumann series polynomial
        LS                       = AZ_ls,            ///< Least-squares polynomial
        GS                       = AZ_sym_GS,        ///< Symmetric Gauss-Seidel
        DD                       = AZ_dom_decomp     ///< Non-overlapping domain decomposition
    };
};

/// Aztec subdomain solver
struct AztecSubdomainSolvers
{
    enum {
        LU                       = AZ_lu,            ///< Sparse LU decomposition
        ILUT                     = AZ_ilut,          ///< Sparse ILUT decomposition (Saad)
        ILU                      = AZ_ilu,           ///< Sparse ILU(k) decomposition
        RILU                     = AZ_rilu,          ///< Sparse ILU(k,w) decomposition
        BILU                     = AZ_bilu,          ///< Sparse block ILU(k) decomposition
        ICC                      = AZ_icc            ///< Sparse ICC(k) decomposition
    };
};

/// Belos solvers
struct BelosSolvers
{
    enum {
        BiCGStab                =  1, ///< BiCGStab solver
        BlockCG                 =  2, ///< Block GC solver
        /// BlockGCRODR             =  3, ///< Block Recycling GMRES solver
        BlockGmres              =  4, ///< Block GMRES solver
        FixedPoint              =  5, ///< Fixed-point solver
        GCRODR                  =  6, ///< Recycling GMRES solver
        GmresPoly               =  7, ///< Hybrid block GMRES solver
        LSQR                    =  8, ///< LSQR solver
        Minres                  =  9, ///< Minres solver
        PCPG                    = 10, ///< PCPG solver
        PseudoBlockCG           = 11, ///< Pseudo Block CG solver
        PseudoBlockGmres        = 12, ///< Pseudo Block GMRES solver
        PseudoBlockStochasticCG = 13, ///< Pseudo Block CG solver
        PseudoBlockTFQMR        = 14, ///< Pseudo Block TFQMR solver
        RCG                     = 15, ///< RCG solver
        TFQMR                   = 16  ///< TFQMR solver
    };
};

struct MLSolvers
{
    enum {
        SA                      =  1, ///< classical smoothed aggregation preconditioners
        NSSA                    =  2, ///< Petrov-Galerkin preconditioner for nonsymmetric systems
        DD                      =  3, ///< 2-level domain decomposition preconditioners
        ///< based on aggregation
        DDLU                    =  4, ///< DD with exact LU decompositions on each subdomain
        DDML                    =  5  ///< 3-level domain decomposition preconditioners,
                                      ///< with coarser spaces defined by aggregation;
    };
};

struct AbstractSolverPrivate;

/** @brief Abstract solver base class

    This abstract base class defines the set of attributes and methods
    that must be implemented by all solvers.
 */
class GISMO_EXPORT AbstractSolver
{
 public:

    /// Constructors
    AbstractSolver();
    explicit AbstractSolver(const SparseMatrix & A);

    /// Destructor
    ~AbstractSolver();

    /// Solves problem for the given a right-hand side vector
    const Vector & solve(const Vector & b);

    /// Returns solution vector
    void getSolution(gsVector<real_t> & sol, const int rank = 0) const;

    /// Returns information about valid parameters
    virtual std::string validParams() const = 0;

    /// Returns information about current parameters
    virtual std::string currentParams() const = 0;

    /// Returns status of the solver
    virtual std::string status() const = 0;

    /// Returns timing of the solver
    virtual std::string timing() const = 0;
    
    /// Sets parameters
    virtual void set(const std::string & name, const int & value) = 0;
    virtual void set(const std::string & name, const bool & value) = 0;
    virtual void set(const std::string & name, const double & value) = 0;
    virtual void set(const std::string & name, const std::string & value) = 0;

    /// Set parameters from option list
    void setOptions(const gsOptionList & opt);

    /// Get parameters into option list
    /// void getOptions(const gsOptionList & opt);
    
 protected:
    /// Solves problem
    virtual void solveProblem() = 0;
    
 protected:
    AbstractSolverPrivate * my;
};

/** @brief Abstract direct solver base class
    
    This abstract base class defines the set of additional attributes
    and methods that must be implemented by all direct solvers.
 */
class GISMO_EXPORT AbstractDirectSolver : public AbstractSolver
{
 public:

    // default constructor inherited by default
};

/** @brief Abstract iterative solver base class
    
    This abstract base class defines the set of additional attributes
    and methods that must be implemented by all iterative solvers.
 */
class GISMO_EXPORT AbstractIterativeSolver : public AbstractSolver
{
 public:
    
    /// Constructors
    AbstractIterativeSolver();
    explicit AbstractIterativeSolver(const SparseMatrix & A);
    
 protected:

    /// Returns number of iterations
    virtual int numIterations() const = 0;

protected:
    
    /// Default tolerance for all iterative solvers
    double tolerance;

    /// Default maximum number of iterations for all iterative solvers
    int    maxIter;    
};

/*    --- Amesos solver ---    */

struct AmesosSolverPrivate;

/** @brief Amesos solver class

    This class warps the Trilinos Amesos package
 */
class GISMO_EXPORT AmesosSolver : public AbstractDirectSolver
{
 public:
    typedef AbstractDirectSolver Base;
    
 public:

    /// Constructor (sparse matrix)
    explicit AmesosSolver(const SparseMatrix & A,
                          const int solver = AmesosSolvers::KLU);

    /// Destructor
    ~AmesosSolver();

    /// Returns information about parameters
    std::string validParams() const;
    std::string currentParams() const;

    /// Returns status and timing of solver
    std::string status() const;
    std::string timing() const;
    
    /// Sets parameters
    void set(const std::string & name, const int & value);
    void set(const std::string & name, const bool & value);
    void set(const std::string & name, const double & value);
    void set(const std::string & name, const std::string & value);
    
 private:

    /// Solves problem
    void solveProblem();
    void solveProblem(const bool noSymbolicFactorization,
                      const bool noNumericFactorization);

private:
    
    AmesosSolverPrivate * myAmesos;
};

/*    --- Aztec solver ---    */

struct AztecSolverPrivate;

/** @brief Actez solver class

    This class wraps the Trilinos Actez solver package
 */
class GISMO_EXPORT AztecSolver : public AbstractIterativeSolver
{
 public:
    typedef AbstractIterativeSolver Base;
    
 public:

    /// Constructor (sparse matrix)
    explicit AztecSolver(const SparseMatrix & A,
                         const int solver   = AztecSolvers::Gmres,
                         const int precond  = AztecPreconds::None,
                         const int subdomain_solver = AztecSubdomainSolvers::ILUT);

    /// Destructor
    ~AztecSolver();

    /// Sets Belos solver as preconditioner
    void setPreconditioner(const BelosSolver & Belos);
    
    /// Sets ML solver as preconditioner
    void setPreconditioner(const MLSolver & ML);
    
    /// Returns information about parameters
    std::string validParams() const;
    std::string currentParams() const;

    /// Returns status and timing of solver
    std::string status() const;
    std::string timing() const;
    
    /// Sets parameters
    void set(const std::string & name, const int & value);
    void set(const std::string & name, const bool & value);
    void set(const std::string & name, const double & value);
    void set(const std::string & name, const std::string & value);

    /// Sets parameters/options for for Aztec solver directly
    void set(const int & option, const int & value);
    void set(const int & param,  const double & value);
    
    /// Returns number of iterations
    int numIterations() const;
    
 private:

    /// Solves problem
    void solveProblem();

 private:
    
    AztecSolverPrivate * myAztec;
};

/*   --- Belos solver ---    */

struct BelosSolverPrivate;

/** @brief Belos solver class

    This class warps the Trilinos Belos package
 */
class GISMO_EXPORT BelosSolver : public AbstractIterativeSolver
{
 public:
    typedef AbstractIterativeSolver Base;

    /// Constructor (sparse matrix)
    explicit BelosSolver(const SparseMatrix & A,
                         const int solver   = BelosSolvers::BiCGStab);

    /// Destructor
    ~BelosSolver();

    /// Returns information about parameters
    std::string validParams() const;
    std::string currentParams() const;

    /// Returns status and timing of solver
    std::string status() const;
    std::string timing() const;
    
    /// Sets parameters
    void set(const std::string & name, const int & value);
    void set(const std::string & name, const bool & value);
    void set(const std::string & name, const double & value);
    void set(const std::string & name, const std::string & value);

    /// Hermitian problem type used by iterative solver
    /// : if matrix is symmetric, specifies it in the linear problem
    void setHermitian();

    /// OverlapLevel: must be >= 0. If Comm.NumProc() == 1, it is ignored.
    int setPreconditioner(const std::string & precType, const SparseMatrix &A,
                          const bool & leftprec, const int & OverlapLevel=0);

    //    void setMLPreconditioner(const SparseMatrix &A);

    /// Returns number of iterations
    int numIterations() const;
    
 private:

    /// Solves problem
    void solveProblem();

 private:

    BelosSolverPrivate * myBelos;

public:
    
    /// Returns pointer to preconditioner operator
    Epetra_Operator * getPrecOperator() const;   
};

/*    --- Multi Level (ML) ---    */

struct MLSolverPrivate;

/** @brief ML solver class

    This class warps the Trilinos ML package
 */
class GISMO_EXPORT MLSolver : public AbstractIterativeSolver
{

 public:
    typedef AbstractIterativeSolver Base;

 public:
    
    /// Constructor
    explicit MLSolver(const SparseMatrix &A,
                      const int solver = MLSolvers::SA);

    /// Destructor
    ~MLSolver();

    /// Returns information about parameters
    std::string validParams() const;
    std::string currentParams() const;

    /// Returns status and timing of the solver
    std::string status() const;
    std::string timing() const;
    
    /// Sets parameters
    void set(const std::string & name, const int & value);
    void set(const std::string & name, const bool & value);
    void set(const std::string & name, const double & value);
    void set(const std::string & name, const std::string & value);

    /// Sets options/parameters for AztecOO solver
    void set(const int & option, const int & value);
    
    /// Returns number of iterations
    int numIterations() const;
    
 private:

    /// Solves problem
    void solveProblem();

 private:

    MLSolverPrivate * myML;

 public:

    /// Returns pointer to preconditioner operator
    Epetra_Operator * getPrecOperator() const;
};

};// namespace solver
};// namespace trilinos
};// namespace gismo


