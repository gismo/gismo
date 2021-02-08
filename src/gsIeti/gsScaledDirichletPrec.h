/** @file gsScaledDirichletPrec.h

    @brief This class represents the sclaed Dirichlet preconditioner.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsUtils/gsSortedVector.h>

namespace gismo
{

/** @brief   This class represents the scaled Dirichlet preconditioner for a IETI problem.
 *
 *  The IETI saddle point system is a system of the form:
 *
 *  Its formal representation is
 *
 *  \f[
 *       \sum_{k=1}^K  \hat B_k  D_k  S_k D_k \hat B_k^\top
 *  \f]
 *
 *  It is a preconditioner for the Schur complement of the IETI system
 *
 *  \f[
 *       \begin{pmatrix}
 *            \tilde A_1 &            &             &            &  \tilde B_1^\top \\
 *                       & \tilde A_2 &             &            &  \tilde B_2^\top \\
 *                       &            &   \ddots    &            &  \vdots   \\
 *                       &            &             & \tilde A_N &  \tilde B_N^\top \\
 *            \tilde B_1 & \tilde B_2 &   \cdots    & \tilde B_N &     0     \\
 *       \end{pmatrix}
 *  \f]
 *
 *  For a standard IETI-dp setup, we additionally have a primal problem, thus N=K+1. In this
 *  case, the matrices \f$ \tilde A_k \f$ and \f$ \tilde B_k \f$ are obtained from the
 *  original matrices \f$ A_k \f$ and \f$ B_k \f$ by eliminating the primal dofs (or by
 *  incorporatung a constraint that sets them to zero).
 *
 *  The matrices \f$ S_k \f$ are stored in the vector \a localSchurOps . To allow certain matrix-free
 *  variats, they are stored in form of a vector of \a gsLinearOperator s. These matrices are the
 *  Schur-complements of the matrices \f$ A_k \f$ with respect to the degrees of freedom on the
 *  skeleton.
 *
 *  The matrices \f$\tilde B_k\f$ are stored in the vector \a jumpMatrices . These matrices usually
 *  differ from the matrices \f$ B_k \f$ from the IETI-system since -- for the preconditioner --
 *  the jump matrices have to be restricted to the Skeleton.
 *
 *  If the matrices \f$ A_k \f$ and \f$ B_k \f$ are given, the function \a restrictToSkeleton allows
 *  to compute the corresponding matrices \f$ S_k \f$ and \f$ \hat B_k \f$. The degrees of freedom
 *  belonging to the skeleton can be specified by the caller. The caller can use the function
 *  \a getSkeletonDofs to extract this information from the jump matrices, i.e., skeleton dofs are
 *  those that are affected by a Lagrange multiplier.
 *
 *  The scaling matrcies \f$ D_k \f$ are stored in the vector \a scalingMatrices. They can be provided
 *  by the caller or generated automatically using \a setupMultiplicityScaling .
 *
 *  @ingroup Solver
**/

template< typename T >
class gsScaledDirichletPrec
{
    typedef gsLinearOperator<T>               Op;
    typedef memory::shared_ptr<Op>            OpPtr;
    typedef gsSparseMatrix<T>                 SparseMatrix;
    typedef gsSparseMatrix<T,RowMajor>        JumpMatrix;
    typedef memory::shared_ptr<JumpMatrix>    JumpMatrixPtr;
    typedef gsMatrix<T>                       Matrix;
public:

    /// @brief Reservs the memory required to store the given number of subdomain
    /// @param n Number of subdomains
    void reserve( index_t n )
    {
        m_jumpMatrices.reserve(n);
        m_localSchurOps.reserve(n);
        m_localScaling.reserve(n);
    }

    /// @briefs Adds a new subdomain
    ///
    /// Subdomain might be, e.g., a patch-local problem or the primal problem
    ///
    /// @param jumpMatrix       The associated jump matrix
    /// @param localSchurOp     The operator that represents the local Schur complement
    ///
    /// @note These tow parameters can also be provided as \a std::pair
    void addSubdomain( JumpMatrix jumpMatrix, OpPtr localSchurOp )
    {
        m_jumpMatrices.push_back(jumpMatrix.moveToPtr());
        m_localSchurOps.push_back(give(localSchurOp));
        m_localScaling.push_back( Matrix() );
    }

    void addSubdomain( std::pair<JumpMatrix,OpPtr> data )
    {
        addSubdomain(give(data.first), give(data.second));
    }

    /// Access the jump matrix
    JumpMatrixPtr&       jumpMatrix(index_t i)           { return m_jumpMatrices[i];  }
    const JumpMatrixPtr& jumpMatrix(index_t i) const     { return m_jumpMatrices[i];  }

    /// Access the local Schur complements operator
    OpPtr&               localSchurOps(index_t i)        { return m_localSchurOps[i]; }
    const OpPtr&         localSchurOps(index_t i) const  { return m_localSchurOps[i]; }

    /// Access the local scaling matrix (as row vector)
    Matrix&              localScaling(index_t i)         { return m_localScaling[i];  }
    const Matrix&        localScaling(index_t i) const   { return m_localScaling[i];  }

    /// @brief Extracts the skeleton dofs from the jump matrix
    ///
    /// This means that a dof is considered to be on the skeleton iff at least one Lagrange
    /// multiplier acts on it. This might lead to other results than the function that is
    /// provided by \a gsIetiMapper.
    static gsSortedVector<index_t> getSkeletonDofs( const JumpMatrix& jm );

    /// Restricts the jump matrix to the given dofs (just takes the corresponding cols)
    static JumpMatrix restrictJumpMatrix( const JumpMatrix& jm, const std::vector<index_t> dofs );

    /// Computes the matrix blocks with respect to the given dofs
    ///
    /// If 0 corresponds to the list of dofs and 1 remains to the others, this function returns
    /// the blocks A00, A10, A01, A11 of A
    static std::vector<SparseMatrix> matrixBlocks( const SparseMatrix & mat, const std::vector<index_t> dofs );

    /// Computes the Schur complement with respect to the given dofs
    ///
    /// Uses sparse Cholesky solver
    static OpPtr schurComplement( const SparseMatrix & mat, const std::vector<index_t> dofs );

    /// Combines \a restrictJumpMatrix and \a schurComplement
    static std::pair<JumpMatrix,OpPtr> restrictToSkeleton(
        const JumpMatrix& jm,
        const SparseMatrix& mat,
        const std::vector<index_t>& dofs
    )
    { return std::pair<JumpMatrix,OpPtr>( restrictJumpMatrix( jm, dofs ), schurComplement( mat, dofs ) ); }

    /// @brief Returns the number of Lagrange multipliers.
    /// This requires that at least one subdomain was defined.
    index_t numberOfLagrangeMultipliers() const
    {
        GISMO_ASSERT( !m_jumpMatrices.empty(), "gsScaledDirichletPrec: Number of Lagrange multipliers "
            "can only be determined if there are jump matrices.");
        return m_jumpMatrices[0]->rows();
    }

    /// @brief This sets up the member vector \a localScaling based on multiplicity scaling
    ///
    /// This requires that the subdomains have been defined first.
    void setupMultiplicityScaling();

    /// @brief This returns the preconditioner as \a gsLinearOperator
    ///
    /// This requires that the subdomains have been defined first.
    OpPtr preconditioner() const;

public:
    std::vector<OpPtr>          m_localSchurOps;
    std::vector<Matrix>         m_localScaling;
    std::vector<JumpMatrixPtr>  m_jumpMatrices;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsScaledDirichletPrec.hpp)
#endif
