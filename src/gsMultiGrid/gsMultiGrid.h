/** @file gsMultiGrid.h

    @brief Multigrid solver for isogeometric discretizations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsSolver/gsPreconditioner.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/** @brief
    Matrix-free multigrid solver isogeometric problems.

    This class implements geometric multigrid for isogeometric patches.
    To use it, pass the operators representing the stiffness matrix on ALL levels
    and the transfer matrices. The smoother HAS TO BE set with setSmoother().

    \par Additional options

    The solver can be configured to use V- or W-cycles using setNumCycles().

    The smoother to use can be set with setSmoother(). Here, a pointer
    to the gsPreconditioner object has to be provided for each grid level.

    The number of pre- and post-smoothing steps can be configured with
    setNumPreSmooth() and setNumPostSmooth(), respectively.

    Full multigrid (FMG) and cascadic multigrid are supported and sometimes
    extremely efficient. It should however be considered
    experimental since the theory is not well understood at this point.

    \ingroup Solver
*/

template<class T>
class gsMultiGridOp : public gsPreconditionerOp<T>
{

public:

    /// Shared pointer for gsMultiGridOp
    typedef memory::shared_ptr<gsMultiGridOp> Ptr;

    /// Unique pointer for gsMultiGridOp
    typedef memory::unique_ptr<gsMultiGridOp> uPtr;

    /// Shared pointer to gsLinearOperator
    typedef typename gsLinearOperator<T>::Ptr OpPtr;

    /// Shared pointer to gsPreconditionerOp
    typedef typename gsPreconditionerOp<T>::Ptr PrecondPtr;

    /// Direct base class
    typedef gsPreconditionerOp<T> Base;

    /// Matrix type
    typedef gsSparseMatrix<T> SpMatrix;

    /// Matrix type
    typedef gsSparseMatrix<T, RowMajor> SpMatrixRowMajor;

    /// Smart pointer to matrix type
    typedef memory::shared_ptr<SpMatrix> SpMatrixPtr;

    /// Smart pointer to matrix type
    typedef memory::shared_ptr<SpMatrixRowMajor> SpMatrixRowMajorPtr;

    /// @brief Constructor
    ///
    /// @param fineMatrix                Stiffness matrix on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    gsMultiGridOp( SpMatrix fineMatrix, std::vector< SpMatrixRowMajor > transferMatrices, OpPtr coarseSolver = OpPtr() );

    /// @brief Constructor
    ///
    /// @param fineMatrix                Stiffness matrix (as smart pointers) on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    gsMultiGridOp( SpMatrixPtr fineMatrix, std::vector< SpMatrixRowMajorPtr > transferMatrices, OpPtr coarseSolver = OpPtr() );

    /// @brief Constructor for a matix-free variant
    ///
    /// @param ops                       Linear operators representing the stiffness matrix on all levels
    /// @param prolong                   Linear operators representing the prolongation operators
    /// @param restrict                  Linear operators representing the restriction operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    gsMultiGridOp( const std::vector<OpPtr>& ops, const std::vector<OpPtr>& prolong, const std::vector<OpPtr>& restrict, OpPtr coarseSolver = OpPtr() );

    /// Make function returning smart pointer
    ///
    /// @param fineMatrix                Stiffness matrix on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    static uPtr make( SpMatrix fineMatrix, std::vector< SpMatrixRowMajor > transferMatrices, OpPtr coarseSolver = OpPtr() )
        { return uPtr( new gsMultiGridOp( give(fineMatrix), give(transferMatrices), give(coarseSolver) ) ); }

    /// Make function returning smart pointer
    ///
    /// @param fineMatrix                Stiffness matrix (as smart pointers) on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    static uPtr make( SpMatrixPtr fineMatrix, std::vector< SpMatrixRowMajorPtr > transferMatrices, OpPtr coarseSolver = OpPtr() )
        { return uPtr( new gsMultiGridOp( give(fineMatrix), give(transferMatrices), give(coarseSolver) ) ); }

    /// Make function returning a shared pointer for a matix-free variant
    ///
    /// @param ops                       Linear operators representing the stiffness matrix on all levels
    /// @param prolong                   Linear operators representing the prolongation operators
    /// @param restrict                  Linear operators representing the restriction operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    static uPtr make( const std::vector<OpPtr>& ops, const std::vector<OpPtr>& prolongation, const std::vector<OpPtr>& restriction, OpPtr coarseSolver = OpPtr() )
    { return uPtr( new gsMultiGridOp( ops, prolongation, restriction, give(coarseSolver) ) ); }

private:
    // Init function that is used by matrix based constructors
    void init( SpMatrixPtr fineMatrix, std::vector< SpMatrixRowMajorPtr > transferMatrices, OpPtr coarseSolver );
    void initCoarseSolver();
public:

    /// Apply smoothing step
    void smoothingStep(index_t level, const gsMatrix<T>& rhs, gsMatrix<T>& x) const;

    /// Apply smoothing step on finest grid
    void smoothingStep(const gsMatrix<T>& rhs, gsMatrix<T>& x) const
    { smoothingStep(finestLevel(), rhs, x); }

    /// Estimates for a smoother I - S^{-1} A the largest eigenvalue of S^{-1} A. Can be used to adjust the damping parameters.
    T estimateLargestEigenvalueOfSmoothedOperator(index_t level, index_t steps = 10)
    { return m_smoother[level]->estimateLargestEigenvalueOfPreconditionedSystem(steps); }

    /// Perform one multigrid cycle on the iterate \a x with right-hand side \a f at the given \a level.
    void multiGridStep(index_t level, const gsMatrix<T>& rhs, gsMatrix<T>& x) const;

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const
    { multiGridStep(finestLevel(), rhs, x); }

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const
    {
        // This should cover most cases:
        std::swap(m_numPreSmooth,m_numPostSmooth);
        step(rhs,x);
        std::swap(m_numPreSmooth,m_numPostSmooth);
    }

    /// Perform one full multigrid cycle and store the resulting solution vector in \a result.
    void fullMultiGrid(const std::vector< gsMatrix<T> >& rhs, const std::vector< gsMatrix<T> >& dirichletIntp, gsMatrix<T>& result) const;

    /// Perform one cascadic multigrid cycle and store the resulting solution vector in \a result.
    void cascadicMultiGrid(const std::vector< gsMatrix<T> >& rhs, const std::vector< gsMatrix<T> >& dirichletIntp, gsMatrix<T>& result) const;

    /// Restrict a vector of coefficients from the level \a lf, given by \a fine, to the next coarser level.
    void restrictVector(index_t lf, const gsMatrix<T>& fine, gsMatrix<T>& coarse) const;

    /// Prolong a vector of coefficients from the level \a lc, given by \a coarse, to the next finer level.
    void prolongVector(index_t lc, const gsMatrix<T>& coarse, gsMatrix<T>& fine) const;

    /// Solve the problem with direct method on coarsest level
    void solveCoarse(const gsMatrix<T>& rhs, gsMatrix<T>& result) const
    { m_coarseSolver->apply( rhs, result ); }

    index_t numLevels() const               { return n_levels;               } ///< Number of levels in the multigrid construction.
    index_t finestLevel() const             { return n_levels - 1;           } ///< The index of the finest level (0-based).

    OpPtr underlyingOp(index_t lvl) const   { return m_ops[lvl];             } ///< Underlying operator (=stiffness matrix) for given level.
    OpPtr underlyingOp() const              { return m_ops[finestLevel()];   } ///< Underlying operator (=stiffness matrix) for finest level.

    /// Set underlying operator (=stiffness matrix) for certain level
    void setUnderlyingOp(index_t lvl, OpPtr op)
    {
        GISMO_ASSERT ( lvl >= 0 && lvl < n_levels, "The given level is not feasible." );
        m_ops[lvl] = op;
    }

    const SpMatrix& matrix(index_t lvl) const;                                  ///< Stiffness matrix for given level.
    const SpMatrix& matrix() const        { return matrix(finestLevel());     } ///< Stiffness matrix for finest level.

    index_t nDofs(index_t lvl) const      { return underlyingOp(lvl)->cols(); } ///< Number of dofs for the given level.
    index_t nDofs()            const      { return nDofs( finestLevel() );    } ///< Number of dofs for the finest level.

    index_t rows() const                  { return underlyingOp()->rows();    }
    index_t cols() const                  { return underlyingOp()->cols();    }

    void setSmoother(index_t lvl, const PrecondPtr& sm);                        ///< Set the smoother
    PrecondPtr smoother(index_t lvl)                { return m_smoother[lvl]; } ///< Get the smoother

    void setCoarseSolver(const OpPtr& sol)          { m_coarseSolver = sol;   } ///< Get the coarse solver
    OpPtr coarseSolver()                            { return m_coarseSolver;  } ///< Get the coarse solver

    void setNumPreSmooth(index_t n)                 { m_numPreSmooth = n;     } ///< Set number of pre-smoothing steps to perform.
    void setNumPostSmooth(index_t n)                { m_numPostSmooth = n;    } ///< Set number of post-smoothing steps to perform.
    void setNumCycles(index_t n)                    { m_numCycles = n;        } ///< Set number of cycles (usually 1 for V-cycle or 2 for W-cycle).
    void setMaxIters(index_t n)                     { m_maxIters = n;         } ///< Set the maximum number of iterations for the member function \a solve
    void setTol(T tol)                              { m_tol = tol;            } ///< Set the error bound for the member function \a solve
    void setCorarseGridCorrectionDamping(T damping) { m_damping = damping;    } ///<Set the damping of for the coarse-grid correction

    static gsOptionList defaultOptions();                                       ///< Returns a list of default options
    virtual void setOptions(const gsOptionList & opt);                          ///< Set the options based on a gsOptionList

protected:

    gsMultiGridOp() {}

private:

    // Number of levels
    index_t n_levels;

    // Underlying operators (=stiffness matrix) on each level
    std::vector< OpPtr > m_ops;

    // Smoothers on each level
    std::vector< PrecondPtr > m_smoother;

    // Prolongation matrix for each grid transition
    std::vector< OpPtr > m_prolong;
    std::vector< OpPtr > m_restrict;

    // solver for the coarsest-grid problem
    OpPtr m_coarseSolver;

    mutable index_t m_numPreSmooth;
    mutable index_t m_numPostSmooth;
    index_t m_numCycles;
    index_t m_maxIters;
    T m_tol;
    T m_damping;

}; // class gsMultiGridOp

}  // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMultiGrid.hpp)
#endif
