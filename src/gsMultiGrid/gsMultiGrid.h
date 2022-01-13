/** @file gsMultiGrid.h

    @brief Multigrid preconditioner

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

/** @brief Multigrid preconditioner
 *
 *  This class implements an generic geometric multigrid framework. This class
 *  implements a preconditioner which can be used, e.g., in an iterative solver,
 *  like \a gsConjugateGradient. If a simple multigrid iteration is required,
 *  use this class as a preconditioner in \a gsGradientMethod with step size 1.
 *
 *  The class expects transfer matrices, which can be constructed with the
 *  class \a gsGridHierarchy.
 *
 *  The intergrid transfers (restriction and prolongation) and the stiffness
 *  matrix are stored as \a gsLinearOperator, which means that also matrix-free
 *  variants are possible. If linear operators are provided, operators
 *  representing the siffness matrix have to be provided for all grid levles.
 *  If the integrid transfers and the stiffness matrix are provided as a sparse
 *  matrix, the matrices for the coarser levels are computed automatically based
 *  on the Galerkin principle.
 *
 *  For all levels, a smoother has to be provided by defining setSmoother().
 *
 *  The solver for the coarsest grid level is created automatically if not
 *  provided by the caller.
 *
 *  The solver can be configured to use V- or W-cycles using setNumCycles().
 *  The number of pre- and post-smoothing steps can be configured with
 *  setNumPreSmooth() and setNumPostSmooth(), respectively.
 *
 *  Also full multigrid and cascadic multigrid algorithms for computing an
 *  initial guess are provided.
 *
 *  @ingroup Solver
**/

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
    typedef gsMatrix<T> Matrix;

    /// Sparse matrix type
    typedef gsSparseMatrix<T> SpMatrix;

    /// Matrix type for transfers
    typedef gsSparseMatrix<T, RowMajor> SpMatrixRowMajor;

    /// Smart pointer to sparse matrix type
    typedef memory::shared_ptr<SpMatrix> SpMatrixPtr;

    /// Smart pointer to matrix type for transfers
    typedef memory::shared_ptr<SpMatrixRowMajor> SpMatrixRowMajorPtr;

    /// @brief Constructor
    ///
    /// @param fineMatrix                Stiffness matrix on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    gsMultiGridOp(
        SpMatrix fineMatrix,
        std::vector<SpMatrixRowMajor> transferMatrices,
        OpPtr coarseSolver = OpPtr()
    );

    /// @brief Constructor
    ///
    /// @param fineMatrix                Stiffness matrix (as smart pointers) on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    gsMultiGridOp(
        SpMatrixPtr fineMatrix,
        std::vector<SpMatrixRowMajorPtr> transferMatrices,
        OpPtr coarseSolver = OpPtr()
    );

    /// @brief Constructor for a matix-free variant
    ///
    /// @param ops                       Linear operators representing the stiffness matrix on all levels
    /// @param prolongation              Linear operators representing the prolongation operators
    /// @param restriction               Linear operators representing the restriction operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    gsMultiGridOp(
        const std::vector<OpPtr>& ops,
        const std::vector<OpPtr>& prolongation,
        const std::vector<OpPtr>& restriction,
        OpPtr coarseSolver = OpPtr()
    );

    /// @brief Make function returning smart pointer
    ///
    /// @param fineMatrix                Stiffness matrix on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    static uPtr make(
        SpMatrix fineMatrix,
        std::vector<SpMatrixRowMajor> transferMatrices,
        OpPtr coarseSolver = OpPtr()
    )
    { return uPtr( new gsMultiGridOp( give(fineMatrix), give(transferMatrices), give(coarseSolver) ) ); }

    /// @brief Make function returning smart pointer
    ///
    /// @param fineMatrix                Stiffness matrix (as smart pointers) on the finest grid
    /// @param transferMatrices          Intergrid transfer matrices representing restriction and prolongation operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    static uPtr make(
        SpMatrixPtr fineMatrix,
        std::vector<SpMatrixRowMajorPtr> transferMatrices,
        OpPtr coarseSolver = OpPtr()
    )
    { return uPtr( new gsMultiGridOp( give(fineMatrix), give(transferMatrices), give(coarseSolver) ) ); }

    /// @brief Make function returning a shared pointer for a matix-free variant
    ///
    /// @param ops                       Linear operators representing the stiffness matrix on all levels
    /// @param prolongation              Linear operators representing the prolongation operators
    /// @param restriction               Linear operators representing the restriction operators
    /// @param coarseSolver              Linear operator representing the exact solver on the coarsest grid level,
    ///                                  defaulted to a direct solver (PartialPivLUSolver)
    static uPtr make(
        const std::vector<OpPtr>& ops,
        const std::vector<OpPtr>& prolongation,
        const std::vector<OpPtr>& restriction,
        OpPtr coarseSolver = OpPtr()
    )
    { return uPtr( new gsMultiGridOp( ops, prolongation, restriction, give(coarseSolver) ) ); }

private:
    /// Init function that is used by matrix based constructors
    void init(SpMatrixPtr fineMatrix, std::vector<SpMatrixRowMajorPtr> transferMatrices);
    /// Init solver on coarsest grid level
    void initCoarseSolver() const;
public:

    /// @brief Apply smoothing steps on the corresponding level
    ///
    /// @param[in]     level   Level
    /// @param[in]     rhs     Right hand side
    /// @param[in,out] x       Current iterate vector
    void smoothingStep(index_t level, const Matrix& rhs, Matrix& x) const;

    /// @brief Apply smoothing steps on finest grid
    ///
    /// @param[in]     rhs     Right hand side
    /// @param[in,out] x       Current iterate vector
    void smoothingStep(const Matrix& rhs, Matrix& x) const
    { smoothingStep(finestLevel(), rhs, x); }

    /// @brief Perform one multigrid cycle, starting from the given level
    ///
    /// @param[in]     level   Level
    /// @param[in]     rhs     Right-hand side
    /// @param[in,out] x       Current iterate vector
    void multiGridStep(index_t level, const Matrix& rhs, Matrix& x) const;

    /// @brief Perform one multigrid cycle
    ///
    /// @param[in]     rhs     Right-hand side
    /// @param[in,out] x       Current iterate vector
    void step(const Matrix& rhs, Matrix& x) const override
    { multiGridStep(finestLevel(), rhs, x); }

    /// @brief Perform one transposed multigrid cycle
    ///
    /// @param[in]     rhs     Right-hand side
    /// @param[in,out] x       Current iterate vector
    ///
    /// This assumes that the stiffness matrix is symmetric, that SymmSmooth is true
    /// and that the smoothers implement stepT correctly.
    void stepT(const Matrix& rhs, Matrix& x) const override
    {
        std::swap(m_numPreSmooth,m_numPostSmooth);
        step(rhs,x);
        std::swap(m_numPreSmooth,m_numPostSmooth);
    }

    /// @brief Perform one full multigrid cycle and store the resulting solution vector in result.
    ///
    /// @param[in]  rhs         The right-hand-sides for all grid levels
    /// @param[in]  fixedValues The values for the fixed dofs (like eliminated Dirichlet values)
    /// @param[out] result      The destination for the result
    void fullMultiGrid(
        const std::vector<Matrix>& rhs,
        const std::vector<Matrix>& dirichletIntp,
        Matrix& result
    ) const;

    /// @brief Perform one cascadic multigrid cycle and store the resulting solution vector in result.
    ///
    /// @param[in]  rhs         The right-hand-sides for all grid levels
    /// @param[in]  fixedValues The values for the fixed dofs (like eliminated Dirichlet values)
    /// @param[out] result      The destination for the result
    void cascadicMultiGrid(
        const std::vector<Matrix>& rhs,
        const std::vector<Matrix>& fixedValues,
        Matrix& result
    ) const;

    /// Restrict a vector of coefficients from the level lf, given by fine, to the next coarser level.
    void restrictVector(index_t lf, const Matrix& fine, Matrix& coarse) const;

    /// Prolong a vector of coefficients from the level lc, given by coarse, to the next finer level.
    void prolongVector(index_t lc, const Matrix& coarse, Matrix& fine) const;

    /// Solve the problem with direct solver on coarsest level
    void solveCoarse(const Matrix& rhs, Matrix& result) const
    {
        if (!m_coarseSolver) initCoarseSolver();
        m_coarseSolver->apply(rhs, result);
    }

    /// Number of levels in the multigrid construction
    index_t numLevels() const                         { return m_nLevels;                 }

    /// The index of the finest level (0-based)
    index_t finestLevel() const                       { return m_nLevels - 1;             }

    /// Set underlying operator (=stiffness matrix) for certain level
    void setUnderlyingOp(index_t lvl, OpPtr op)
    {
        GISMO_ASSERT ( lvl >= 0 && lvl < m_nLevels, "gsMultiGridOp: The given level "<<lvl<<" is not feasible." );
        m_ops[lvl] = op;
    }

    /// Underlying operator (=stiffness matrix) for given level
    OpPtr underlyingOp(index_t lvl) const             { return m_ops[lvl];               }

    /// Underlying operator (=stiffness matrix) for finest level
    OpPtr underlyingOp() const override               { return m_ops[finestLevel()];     }

    /// Stiffness matrix for given level
    const SpMatrix& matrix(index_t lvl) const;

    /// Stiffness matrix for finest level
    const SpMatrix& matrix() const                    { return matrix(finestLevel());     }

    /// Number of dofs for the given level
    index_t nDofs(index_t lvl) const                  { return underlyingOp(lvl)->cols(); }

    /// Number of dofs for the finest level
    index_t nDofs() const                             { return nDofs( finestLevel() );    }

    // For docs, see base class
    index_t rows() const override                     { return underlyingOp()->rows();    }

    // For docs, see base class
    index_t cols() const override                     { return underlyingOp()->cols();    }

    /// Set the smoother
    /// @param lvl  The corresponding level
    /// @param sm   The smoother
    void setSmoother(index_t lvl, const PrecondPtr& sm);

    /// Get the smoother
    /// @param lvl  The corresponding level
    PrecondPtr smoother(index_t lvl) const            { return m_smoother[lvl];           }

    /// Set the solver for the coarsest problem (level 0)
    void setCoarseSolver(const OpPtr& sol)            { m_coarseSolver = sol;             }

    /// Get the solver for the coarsest problem (level 0)
    OpPtr coarseSolver() const                        { return m_coarseSolver;            }

    /// Set number of pre-smoothing steps to perform
    void setNumPreSmooth(index_t n)                   { m_numPreSmooth = n;               }

    /// Get number of pre-smoothing steps to perform
    index_t setNumPreSmooth() const                   { return m_numPreSmooth;            }

    /// Set number of post-smoothing steps to perform
    void setNumPostSmooth(index_t n)                  { m_numPostSmooth = n;              }

    /// Get number of post-smoothing steps to perform
    index_t setNumPostSmooth() const                  { return m_numPostSmooth;           }

    /// Set symmetric post-smoothing option
    /// @param s       Iff true (default), stepT is called for post-smoothing
    void setSymmSmooth(bool s)                        { m_symmSmooth = s;                 }

    /// Get symmetric post-smoothing option
    bool symmSmooth() const                           { return m_symmSmooth;              }

    /// Set number of coarse grid steps to be applied
    /// @param n       Number of levels (1 for V-cycle, 2 for W-cycle)
    void setNumCycles(index_t n)
    {
        m_numCycles.setConstant(m_nLevels-1,n);
        // The direct solver on coarsest level is only invoked once
        m_numCycles[0] = 1;
    }

    /// Set number of coarse grid steps to be applied
    /// @param lvl   Level of coarser grid (multiGridStep on level lvl+1 calls
    ///              n times multiGridStep for level lvl)
    /// @param n     Number of calls (1 for V-cycle, 2 for W-cycle)
    void setNumCycles(index_t lvl, index_t n)
    {
        GISMO_ASSERT(lvl>=0&&lvl<m_nLevels-1, "gsMultiGrid::setNumCycles: Givel level out of bounds.");
        m_numCycles[lvl] = n;
    }

    /// Get number of coarse grid steps to be applied
    index_t numCycles(index_t lvl) const              { return m_numCycles[lvl];           }

    /// Set the damping of for the coarse-grid correction
    void setCorarseGridCorrectionDamping(T damping)   { m_damping = damping;               }

    /// Get the damping of for the coarse-grid correction
    T corarseGridCorrectionDamping() const            { return m_damping;                  }

    /// Returns a list of default options
    static gsOptionList defaultOptions();

    /// Set the options based on a gsOptionList
    void setOptions(const gsOptionList& opt) override;

private:

    /// Number of levels
    index_t m_nLevels;

    /// Underlying operators (=stiffness matrix) on each level
    std::vector<OpPtr> m_ops;

    /// Smoothers on each level
    std::vector<PrecondPtr> m_smoother;

    /// Prolongation matrix for each grid transition
    std::vector<OpPtr> m_prolong;

    /// Restriction matrix for each grid transition
    std::vector<OpPtr> m_restrict;

    /// Solver for the coarsest-grid problem
    mutable OpPtr m_coarseSolver;

    /// Number of pre-smoothing steps
    mutable index_t m_numPreSmooth;

    /// Number of post-smoothing steps
    mutable index_t m_numPostSmooth;

    /// Symmmetric smoothing
    bool m_symmSmooth;

    /// Number of coarse-grid steps (1=V-cycle, 2=W-cycle)
    gsVector<index_t> m_numCycles;

    /// Damping for the coarse-grid correction
    T m_damping;

}; // class gsMultiGridOp

}  // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMultiGrid.hpp)
#endif
