/** @file gsIetiSystem.h

    @brief This class represents a IETI system

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

/** @brief   This class represents a IETI system.
 *
 *  The IETI saddle point system is a system of the form:
 *
 *  \f[
 *     \begin{pmatrix}
 *        \tilde A_1 &            &             &            &  \tilde B_1^\top \\
 *                   & \tilde A_2 &             &            &  \tilde B_2^\top \\
 *                   &            &   \ddots    &            &  \vdots          \\
 *                   &            &             & \tilde A_N &  \tilde B_N^\top \\
 *        \tilde B_1 & \tilde B_2 &   \cdots    & \tilde B_N &     0            \\
 *     \end{pmatrix}
 *  \f]
 *
 *  The corresponding Schur complement is
 *
 *  \f[
 *       \sum_{k=1}^N   \tilde B_k   \tilde A_k^{-1}  \tilde B_k^\top
 *  \f]
 *
 *  For a standard IETI-DP setup, \f$ \tilde A_k \f$ and \f$ \tilde B_k \f$ are
 *  obtained from the original matrices \f$ A_k \f$ and \f$ B_k \f$ by
 *  eliminating the primal dofs or by incorporating a constraint that sets them
 *  to zero.
 *
 *  This class does not have any special treatment for the primal problem of a
 *  IETI-DP solver. Thus, the primal problem is just another subdomain and in
 *  case of IETI-DP, we have N=K+1, where K is the number of patches.
 *
 *  The matrices \f$ \tilde A_k \f$ are stored in a vector that can be accessed
 *  (read and write) via \ref localMatrixOp. To allow certain matrix-free
 *  variants, each of them is stored as \a gsLinearOperator.
 *
 *  The inverses \f$ \tilde A_k^{-1} \f$ are stored as vector accessable via
 *  \ref localSolverOp. If the matrices \f$ \tilde A_k \f$ are stored as
 *  \a gsMatrixOp< gsSparseMatrix<T> >, LU solvers will be automatically created
 *  on the fly if needed. Otherwise or if the caller wants other local solvers
 *  (like inexact ones), the vector can be populated by the caller.
 *
 *  The matrices \f$ \tilde B_k \f$ are stored in a vector accessible via
 *  \ref jumpMatrix.
 *
 *  The right-hand sides are stored in a vector accessible via \ref localRhs.
 *
 *  @ingroup Solver
**/

template< typename T >
class gsIetiSystem
{
    typedef gsLinearOperator<T>               Op;              ///< Linear operator
    typedef memory::shared_ptr<Op>            OpPtr;           ///< Shared pointer to linear operator
    typedef gsSparseMatrix<T>                 SparseMatrix;    ///< Sparse matrix type
    typedef gsMatrixOp<SparseMatrix>          SparseMatrixOp;  ///< Sparse matrix type as operatpr
    typedef gsSparseMatrix<T,RowMajor>        JumpMatrix;      ///< Sparse matrix type for jumps
    typedef memory::shared_ptr<JumpMatrix>    JumpMatrixPtr;   ///< Shared pointer to sparse matrix type for jumps
    typedef gsMatrix<T>                       Matrix;          ///< Matrix type
public:

    /// @brief Reserves the memory required to store the number of subdomains
    /// @param n Number of subdomains
    void reserve(index_t n);

    /// @briefs Adds a new subdomain
    ///
    /// Subdomain might be, e.g., a patch-local problem or the primal problem
    ///
    /// @param jumpMatrix       The associated jump matrix
    /// @param localMatrixOp    The operator that represents the local stiffness matrix
    /// @param localRhs         The contribution to the right-hand side
    /// @param localSolverOp    The operator that represents a solver for the
    ///                         local problem. This parameter is optional; the
    ///                         solver is created automatically if needed.
    void addSubdomain(JumpMatrixPtr jumpMatrix, OpPtr localMatrixOp,
        Matrix localRhs, OpPtr localSolverOp = OpPtr());

    /// Access the jump matrix
    JumpMatrixPtr&       jumpMatrix(index_t k)           { return m_jumpMatrices[k];   }
    const JumpMatrixPtr& jumpMatrix(index_t k) const     { return m_jumpMatrices[k];   }

    /// Access the local stiffness matrix (as \a gsLinearOperator)
    OpPtr&               localMatrixOp(index_t k)        { return m_localMatrixOps[k]; }
    const OpPtr&         localMatrixOp(index_t k) const  { return m_localMatrixOps[k]; }

    /// Access the local right-hand side
    Matrix&              localRhs(index_t k)             { return m_localRhs[k];       }
    const Matrix&        localRhs(index_t k) const       { return m_localRhs[k];       }

    /// Access the local solver operator
    OpPtr&               localSolverOp(index_t k)        { return m_localSolverOps[k]; }
    const OpPtr&         localSolverOp(index_t k) const  { return m_localSolverOps[k]; }

    /// @brief Returns the number of Lagrange multipliers
    ///
    /// This requires that at least one jump matrix has been set.
    index_t nLagrangeMultipliers() const
    {
        GISMO_ASSERT( m_jumpMatrices.size()>0, "gsIetiSystem: Number of Lagrange multipliers "
            "can only be determined if there are jump matrices.");
        return m_jumpMatrices[0]->rows();
    }

    /// @brief Returns \a gsLinearOperator that represents the Schur complement
    ///        for the IETI problem
    ///
    /// If the local solvers have not been provided, this function will generate
    /// them from the \ref localMatrixOp if they are \a gsMatrixOp<gsSparseMatrix<T>>
    OpPtr schurComplement() const;

    /// @brief Returns the right-hand-side that is required for the Schur
    ///        complement formulation of the IETI problem
    ///
    /// If the local solvers have not been provided, this function will generate
    /// them from the \ref localMatrixOp if they are \a gsMatrixOp<gsSparseMatrix<T>>
    Matrix rhsForSchurComplement() const;

    /// @brief Returns the local solutions for the individual subdomains
    ///
    /// @param multipliers  The Lagrange multipliers previously computed
    ///                     (based on the Schur complement form)
    ///
    /// If the local solvers have not been provided, this function will generate
    /// them from the \ref localMatrixOp if they are \a gsMatrixOp<gsSparseMatrix<T>>
    std::vector<Matrix> constructSolutionFromLagrangeMultipliers(const Matrix& multipliers) const;

    /// @brief Returns \a gsLinearOperator that represents the IETI problem as
    ///        saddle point problem
    OpPtr saddlePointProblem() const;

    /// @brief Returns the right-hand-side that is required for the saddle point
    ///        formulation of the IETI problem
    Matrix rhsForSaddlePoint() const;

    /// @brief Returns the local solutions for the individual subdomains
    ///
    /// @param x  The solution computed using the saddle point system
    std::vector<Matrix> constructSolutionFromSaddlePoint(const Matrix& x) const;


private:
    void setupSparseLUSolvers() const;                ///< Setup solvers if not provided by user

    std::vector<JumpMatrixPtr>  m_jumpMatrices;       ///< Stores the jump matrices
    std::vector<OpPtr>          m_localMatrixOps;     ///< Stores the local matrix ops \f$ \tilde A_k \f$
    std::vector<Matrix>         m_localRhs;           ///< Stores the local right-hand sides
    mutable std::vector<OpPtr>  m_localSolverOps;     ///< Stores the local solvers
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIetiSystem.hpp)
#endif
