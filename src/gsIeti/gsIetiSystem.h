/** @file gsIetiSystem.h

    @brief This class represents a IETI problem. Its algorithms allow to set up a IETI solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

/** @brief   This class represents a IETI problem. Its algorithms allow to set up a IETI solver.
 *
 *  The IETI saddle point system is a system of the form:
 *
 *  \f[
 *       \begin{pmatrix}
 *            \tilde A_1 &            &             &            &  \tilde B_1^\top \\
 *                       & \tilde A_2 &             &            &  \tilde B_2^\top \\
 *                       &            &   \ddots    &            &  \vdots   \\
 *                       &            &             & \tilde A_K &  \tilde B_K^\top \\
 *            \tilde B_1 & \tilde B_2 &   \cdots    & \tilde B_K &     0     \\
 *       \end{pmatrix}
 *  \f]
 *
 *  The correspondung Schur complement is
 *
 *  \f[
 *       \sum_{k=1}^K   \tilde B_k   \tilde A_k^{-1}  \tilde B_k^\top
 *  \f]
 *
 *  For a standard IETI-dp setup, the matrices \f$ \tilde A_k \f$ and \f$ \tilde B_k \f$ are obtained
 *  from the original matrices \f$ A_k \f$ and \f$ B_k \f$ by eliminating the primal dofs (or by
 *  incorporatung a constraint that sets them to zero).
 *
 *  The matrices \f$ \tilde A_k \f$ are stored in the vector \a localMatrixOps . To allow certain
 *  matrix-free variats, they are stored in form of a vector of \a gsLinearOperator s.
 *
 *  The inverses \f$ \tilde A_k^{-1} \f$ are stored in the vector \a localSolverOps . As far as the
 *  matrices \f$ \tilde A_k\f$ are stored as \a gsMatrixOp containing \a gsSparseMatrix , LU solvers can
 *  be automatically created by calling \a setupSparseLUSolvers. Otherwise or if the caller wants other
 *  local solvers (like inexact ones), the vector can be populated by the caller.
 *
 *  The matrices \f$ \tilde B_k \f$ are stored in the vector \a jumpMatrices .
 *
 *  The right-hand sides are stored in the vector \a localRhs .
 *
 *  @note This class does not have any special treatment for the primal problem of a IETI-dp solver.
 *
 *  @ingroup Solver
**/

template< typename T >
class gsIetiSystem
{
    typedef gsLinearOperator<T>               Op;
    typedef memory::shared_ptr<Op>            OpPtr;
    typedef gsSparseMatrix<T>                 SparseMatrix;
    typedef gsMatrixOp< gsSparseMatrix<T> >   SparseMatrixOp;
    typedef gsSparseMatrix<T,RowMajor>        JumpMatrix;
    typedef memory::shared_ptr<JumpMatrix>    JumpMatrixPtr;
    typedef gsMatrix<T>                       Matrix;
public:

    /// @brief Reservs the memory required to store the given number of subdomain
    /// @param n Number of subdomains
    void reserve(index_t n);

    /// @briefs Adds a new subdomain
    ///
    /// Subdomain might be, e.g., a patch-local problem or the primal problem
    ///
    /// @param jumpMatrix       The associated jump matrix
    /// @param localMatrixOp    The operator that represents the local problem
    /// @param localRhs         The contribution to the right-hand side
    void addSubdomain(JumpMatrixPtr jumpMatrix, OpPtr localMatrixOp, Matrix localRhs);

    /// Access the vector of jump matrices
    std::vector<JumpMatrixPtr>&       jumpMatrices()          { return m_jumpMatrices;   }
    const std::vector<JumpMatrixPtr>& jumpMatrices() const    { return m_jumpMatrices;   }

    /// Access the vector of local system matrices (as \a gsLinearOperator)
    std::vector<OpPtr>&               localMatrixOps()        { return m_localMatrixOps; }
    const std::vector<OpPtr>&         localMatrixOps() const  { return m_localMatrixOps; }

    /// Access the vector of local right-hand sides
    std::vector<Matrix>&              localRhs()              { return m_localRhs;       }
    const std::vector<Matrix>&        localRhs() const        { return m_localRhs;       }

    /// Access the vector of local solver operators
    std::vector<OpPtr>&               localSolverOps()        { return m_localSolverOps; }
    const std::vector<OpPtr>&         localSolverOps() const  { return m_localSolverOps; }

    /// @brief Populates the member \a m_localSolverOps
    ///
    /// This function assums that \a m_jumpMatrices and \a m_localMatrixOps have been populated
    /// first. Moreover, it requres that all \a localMatrixOp s are actually of type
    /// \a gsMatrixOp<gsSparseMatrix<T>> .
    void setupSparseLUSolvers();

    /// Returns the number of Lagrange multipliers
    ///
    /// This requires that at least one jump matrix has been set.
    index_t numberOfLagrangeMultipliers() const
    {
        GISMO_ASSERT( m_jumpMatrices.size()>0, "gsIetiSystem: Number of Lagrange multipliers "
            "can only be determined if there are jump matrices.");
        return m_jumpMatrices[0]->rows();
    }

    /// @brief Returns \a gsLinearOperator that represents the IETI problem as saddle point problem
    ///
    /// This requires that the jump matrices (\a jumpMatrices) and the local matrices
    /// (\a localMatrixOps) have been provided.
    OpPtr saddlePointProblem() const;

    /// @brief Returns \a gsLinearOperator that represents the Schur complement for the IETI problem
    ///
    /// This requires that the jump matrices (\a jumpMatrices) and the local solvers
    /// (\a localSolverOps) have been provided. The local solvers can be automatically generated
    /// by calling the member function \a setupSparseLUSolvers .
    OpPtr schurComplement() const;

    /// @brief Returns the right-hand-side that is required for the Schur complement formulation of the IETI problem
    ///
    /// This requires that the jump matrices (\a jumpMatrices) and the local solvers
    /// (\a localSolverOps) have been provided. The local solvers can be automatically generated
    /// by calling the member function \a setupSparseLUSolvers .
    Matrix rhsForSchurComplement() const;

    /// @brief Returns the local solutions for the individual subdomains
    ///
    /// @param multipliers  The Lagrange multipliers previously computed (based on the Schur complement form)
    ///
    /// This requires that the jump matrices (\a jumpMatrices), the local right-hand sides (\a localRhs) and
    /// the local solvers (\a localSolverOps) have been provided. The local solvers can be automatically
    /// generated by calling the member function \a setupSparseLUSolvers .
    std::vector<Matrix> constructSolutionFromLagrangeMultipliers(const gsMatrix<T>& multipliers) const;

private:
    std::vector<JumpMatrixPtr>  m_jumpMatrices;
    std::vector<OpPtr>          m_localMatrixOps;
    std::vector<Matrix>         m_localRhs;
    std::vector<OpPtr>          m_localSolverOps;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIetiSystem.hpp)
#endif
