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
class AmesosSolver;
class AztecSolver;
class BelosSolver;
class MLSolver;

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

//constants for solver types - from az_aztec_defs.h
#define AZ_cg               0 /* preconditioned conjugate gradient method     */
#define AZ_gmres            1 /* preconditioned gmres method                  */
#define AZ_cgs              2 /* preconditioned cg squared method             */
#define AZ_tfqmr            3 /* preconditioned transpose-free qmr method     */
#define AZ_bicgstab         4 /* preconditioned stabilized bi-cg method       */
#define AZ_slu              5 /* super LU direct method.                      */
#define AZ_symmlq           6 /* indefinite symmetric like symmlq             */
#define AZ_GMRESR           7 /* recursive GMRES (not supported)              */
#define AZ_fixed_pt         8 /* fixed point iteration                        */
#define AZ_analyze          9 /* fixed point iteration                        */
#define AZ_lu              10 /* sparse LU direct method. Also used for a     */
#define AZ_cg_condnum      11
#define AZ_gmres_condnum   12

//constants for preconditioners - from az_aztec_defs.h
#define AZ_none             0 /* no preconditioning. Note: also used for      */
/* scaling, output, overlap options options     */
#define AZ_Jacobi           1 /* Jacobi preconditioning. Note: also used for  */
/* scaling options                              */
#define AZ_sym_GS           2 /* symmetric Gauss-Siedel preconditioning       */
#define AZ_Neumann          3 /* Neumann series polynomial preconditioning    */
#define AZ_ls               4 /* least-squares polynomial preconditioning     */
#define AZ_ilu              6 /* domain decomp with  ilu in subdomains        */
#define AZ_bilu             7 /* domain decomp with block ilu in subdomains   */
/* #define AZ_lu           10    domain decomp with   lu in subdomains        */
#define AZ_icc              8 /* domain decomp with incomp Choleski in domains*/
#define AZ_ilut             9 /* domain decomp with ilut in subdomains        */
#define AZ_rilu            11 /* domain decomp with rilu in subdomains        */
#define AZ_recursive       12 /* Recursive call to AZ_iterate()               */
#define AZ_smoother        13 /* Recursive call to AZ_iterate()               */
#define AZ_dom_decomp      14 /* Domain decomposition using subdomain solver  */
/* given by options[AZ_subdomain_solve]         */
#define AZ_multilevel      15 /* Do multiplicative domain decomp with coarse  */
/* grid (not supported).                        */
#define AZ_user_precond    16 /*  user's preconditioning */
/* Begin Aztec 2.1 mheroux mod */
#define AZ_bilu_ifp        17 /* dom decomp with bilu using ifpack in subdom  */
/* End Aztec 2.1 mheroux mod */

/// Aztec solvers
struct AztecSolvers
{
    enum {
        CG                       = AZ_cg,            ///< CG solver
        CGCondNum                = AZ_cg_condnum,    ///< CG solver with cond. num. estimation
        Gmres                    = AZ_gmres,         ///< GMRES solver
        GmresCondNum             = AZ_gmres_condnum, ///< GMRES solver with cond.num. estimation
        GmresR                   = AZ_GMRESR,        ///< recursive GMRES
        CGS                      = AZ_cgs,           ///< CG squared solver
        TFQMR                    = AZ_tfqmr,         ///< TFQMR solver
        BiCGStab                 = AZ_bicgstab,      ///< BiCGStab solver
        LU                       = AZ_lu,            ///< Sparse direct LU solver
        SLU                      = AZ_slu,           ///< Super LU direct solver
        SymmLQ                   = AZ_symmlq,        ///< Indefinite symmetric like SymmLQ
        FixedPoint               = AZ_fixed_pt,      ///< Fixed-point iteration
        Analyze                  = AZ_analyze        ///< Fixed-point iteration
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
        DD                       = AZ_dom_decomp,    ///< Non-overlapping domain decomposition
        ILU                      = AZ_ilu,           ///< DD with ILU in subdomains
        BILU                     = AZ_bilu,          ///< DD with block-ILU in subdomains
        ICC                      = AZ_icc,           ///< DD with incomplete Choleski in subdomains
        ILUT                     = AZ_ilut,          ///< DD with ILUT in subdomains
        RILU                     = AZ_rilu,          ///< DD with RILU in subdomains
        Recursive                = AZ_recursive,     ///< Recursive call to iteration
        Smoother                 = AZ_smoother,      ///< Recursive call to iteration
        Multilevel               = AZ_multilevel,    ///< DD with coarse grid
        UserDefined              = AZ_user_precond,  ///< User defined preconditioner
        IfpackBILU               = AZ_bilu_ifp       ///< DD with BILU (Ifpack version)
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

#undef AZ_cg
#undef AZ_gmres
#undef AZ_cgs
#undef AZ_tfqmr
#undef AZ_bicgstab
#undef AZ_slu
#undef AZ_symmlq
#undef AZ_GMRESR
#undef AZ_fixed_pt
#undef AZ_analyze
#undef AZ_lu
#undef AZ_cg_condnum
#undef AZ_gmres_condnum
//
#undef AZ_none
#undef AZ_Jacobi
#undef AZ_sym_GS
#undef AZ_Neumann
#undef AZ_ls
#undef AZ_ilu
#undef AZ_bilu
#undef AZ_icc
#undef AZ_ilut
#undef AZ_rilu
#undef AZ_recursive
#undef AZ_smoother
#undef AZ_dom_decomp
#undef AZ_multilevel
#undef AZ_user_precond
#undef AZ_bilu_ifp


/// Belos solvers
struct BelosSolvers
{
    enum {
        BiCGStab                =  1, ///< BiCGStab solver
        BlockCG                 =  2, ///< Block GC solver
#ifdef Belos_ENABLE_Experimental
        BlockGCRODR             =  3, ///< Block Recycling GMRES solver
#endif
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
        DDML                    =  5, ///< 3-level domain decomposition preconditioners,
                                      ///< with coarser spaces defined by aggregation
        DDMLLU                  =  6  ///< 3-level domain decomposition preconditioners,
                                      ///< with coarser spaces defined by aggregation
                                      ///< and exact LU decomposition on each subdomain
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

    AbstractDirectSolver() { }
    explicit AbstractDirectSolver(const SparseMatrix & A)
        : AbstractSolver(A) { }

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

}// namespace solver
}// namespace trilinos
}// namespace gismo
