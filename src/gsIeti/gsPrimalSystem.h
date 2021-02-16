/** @file gsPrimalSystem.h

    @brief This class represents the primal system for a IETI-dp algorithm

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsMatrixOp.h>
#include <gsMatrix/gsVector.h>

namespace gismo
{


/** @brief   This class represents the primal system for a IETI-dp algorithm
 *
 *  The class \a gsIetiSystem does not know anything about the primal problem.
 *  For that class the primal problem is just another subdomain.
 *
 *  The class at hand allows to handle primal degrees of freedom (dofs). It is
 *  purely algebraic and assumes that the caller provides the constraint matrix
 *  C (or, more precisely, its individual rows) and a mapping (as vector of
 *  indices) that tells the class to which primal dof the corresponding
 *  constraint belongs. So, for example, if a primal dof is a vertex value, then
 *  that vertex gets one identifying index identifying it as one single primal
 *  dof, while there is a constraint for each patch the corner belongs to.
 *
 *  It is assumed that the class receives the matrix \f$ A \f$ and the member
 *  \ref incorporateConstraints constructs the local saddle point matrix:
 *
 *  \f[
 *       \tilde A
 *       =
 *       \begin{pmatrix}
 *            A   &    C^\top \\
 *            C   &
 *       \end{pmatrix}
 *  \f]
 *
 *  This matrix is then to be handed over to the \a gsIetiSystem. (Simultaneously,
 *  the local rhs and the local jump matrix are amended).
 *
 *  The member \ref primalBasis allows to construct the basis for the primal
 *  space by the principle of energy minimization.
 *
 *  The member \ref addContribution add the local contributions to the primal
 *  problem.
 *
 *  The member \ref handleConstraints combines the member
 *  \ref incorporateConstraints, the setup of a sparse LU solver for the local
 *  system, and the members \a primalBasis and \ref addContribution into one
 *  single function.
 *
 *  After going through all patches, the class \a gsIetiSystem gets one more
 *  subdomain: the primal one (which is set up based on \ref jumpMatrix,
 *  \ref localMatrix and \ref localRhs).
 *
 *  After solving, the member \ref distributePrimalSolution distributes the
 *  solution obtained for the primal problem back to the individual patches.
 *
 *  @ingroup Solver
**/

template< typename T >
class gsPrimalSystem
{
private:
    typedef gsLinearOperator<T>               Op;              ///< Linear operator
    typedef memory::shared_ptr<Op>            OpPtr;           ///< Shared pointer to linear operator
    typedef gsSparseMatrix<T>                 SparseMatrix;    ///< Sparse matrix type
    typedef gsSparseMatrix<T,RowMajor>        JumpMatrix;      ///< Sparse matrix type for jumps
    typedef gsSparseVector<T>                 SparseVector;    ///< Sparse vector type
    typedef gsMatrix<T>                       Matrix;          ///< Matrix type
public:

    /// @brief Constructor
    ///
    /// @param nPrimalDofs                Number of primal constraints in total
    gsPrimalSystem(index_t nPrimalDofs);

    /// @brief Incorporates the given constraints in the local system
    ///
    /// @param[in]  primalConstraints     Primal constraints; the given vectors
    ///                                   make up the matrix \f$ C_k \f$
    /// @param[in]  jumpMatrix            Jump matrix \f$ B_k \f$
    /// @param[in]  localMatrix           Local matrix \f$ A_k \f$
    /// @param[in]  localRhs              Local right-hand-side \f$ f_k \f$
    /// @param[out] modifiedJumpMatrix    Jump matrix \f$ \tilde B_k \f$
    /// @param[out] modifiedLocalMatrix   Local matrix \f$ \tilde A_k \f$
    /// @param[out] modifiedLocalRhs      Local right-hand-side \f$ \tilde f_k \f$
    /// @param[out] rhsForBasis           Right-hand-side for basis computation
    ///
    /// The output is for the local matrix,
    ///  \f[
    ///      \tilde{A}_k =
    ///      \begin{pmatrix}
    ///           A_k   &  C_k^\top \\  C_k   &  0
    ///      \end{pmatrix}
    ///  \f]
    /// for the jump matrix
    ///  \f[
    ///      \tilde{B}_k =
    ///      \begin{pmatrix}
    ///           B_k   &   0
    ///      \end{pmatrix}
    ///  \f]
    /// and the right-hand side
    ///  \f[
    ///      \tilde f_k =
    ///      \begin{pmatrix}
    ///           f_k \\ 0
    ///      \end{pmatrix}
    ///  \f]
    ///
    static void incorporateConstraints(
        const std::vector<SparseVector>& primalConstraints,
        bool eliminatePointwiseConstraints,
        const SparseMatrix& localMatrix,
        SparseMatrix& modifiedLocalMatrix,
        SparseMatrix& localEmbedding,
        SparseMatrix& embeddingForBasis,
        Matrix& rhsForBasis
    );

    /// @brief Returns the matrix representation of the energy minimizing primal basis
    ///
    /// @param  localSaddlePointSolver  Solver that realizes \f$ \tilde{A}_k^{-1} \f$
    /// @param  rhsForBasis             The right-hand side required for the computation
    ///                                 of the basis as provided by \ref incorporateConstraints
    /// @param  primalDofIndices        Vector, that contains for every primal
    ///                                 constraint the index of the respective
    ///                                 primal dof
    /// @param  nPrimalDofs             The total number of primal dofs
    static gsSparseMatrix<T> primalBasis(
        OpPtr localSaddlePointSolver,
        const SparseMatrix& embeddingForBasis,
        const Matrix& rhsForBasis,
        const std::vector<index_t>& primalDofIndices,
        index_t nPrimalDofs
    );

    /// @brief Adds contributions for a patch to the data hold in the class
    ///
    /// @param  jumpMatrix        Jump matrix \f$ B_k \f$
    /// @param  localMatrix       Local stiffness matrix \f$ A_k \f$
    /// @param  localRhs          Local right-hand side \f$ f_k \f$
    /// @param  primalBasis       Matrix representation of the primal basis
    /// @param  embedding         The map \f$ \tilde u_k \f$ to \f$ u_k \f$; if not
    ///                           provided, defaulted to (possibly rectangular) identity
    void addContribution(
        const JumpMatrix& jumpMatrix,
        const SparseMatrix& localMatrix,
        const Matrix& localRhs,
        SparseMatrix primalBasis,
        OpPtr embedding = OpPtr()
    );

    /// @brief Convenience function for handling the primal constraints
    ///
    /// @param[in]     primalConstraints  Primal constraints; the given vectors
    ///                                   make up the matrix \f$ C_k \f$
    /// @param[in]     primalDofIndices   Vector, that contains for every primal
    ///                                   constraint the index of the respective
    ///                                   primal dof
    /// @param[in,out] jumpMatrix         Jump matrix \f$ B_k \f$
    /// @param[in,out] localMatrix        Local stiffness matrix \f$ A_k \f$
    /// @param[in,out] localRhs           Local right-hand side \f$ f_k \f$
    ///
    /// The implementation is basically:
    /// @code{.cpp}
    ///     //TODO
    /// @endcode
    void handleConstraints(
        const std::vector<SparseVector>& primalConstraints,
        const std::vector<index_t>& primalDofIndices,
        JumpMatrix& jumpMatrix,
        SparseMatrix& localMatrix,
        Matrix& localRhs
    );

    /// @brief  Distributes the given solution for K+1 subdomains to the K patches
    ///
    /// @param    sol   The solution, first for the K patches, followed by the
    ///                 contribution for the primal dofs. The solutions for the K
    ///                 patches is expected to have first the values for all patch-local
    ///                 degrees of freedom, possibly followed by degrees of freedom
    ///                 from Lagrange-mutlipliers used for enforcing primal constraints.
    /// @returns        The solution for the K patches
    std::vector<Matrix> distributePrimalSolution( std::vector<Matrix> sol );

    /// @brief Returns the jump matrix for the primal problem
    JumpMatrix&                           jumpMatrix()        { return m_jumpMatrix;         }
    const JumpMatrix&                     jumpMatrix() const  { return m_jumpMatrix;         }

    /// @brief Returns the local stiffness matrix for the primal problem
    SparseMatrix&                         localMatrix()       { return m_localMatrix;        }
    const SparseMatrix&                   localMatrix() const { return m_localMatrix;        }

    /// @brief Returns the right-hand-side for the primal problem
    Matrix&                               localRhs()          { return m_localRhs;           }
    const Matrix&                         localRhs() const    { return m_localRhs;           }

    /// @brief Returns the size of the primal problem (number of primal dofs)
    index_t nPrimalDofs() const                               { return m_localMatrix.rows();            }

    /// @brief TODO
    index_t eliminatePointwiseConstraints() const             { return m_eliminatePointwiseConstraints; }

    /// @brief TODO
    void setEliminatePointwiseConstraints(bool v)             { m_eliminatePointwiseConstraints = v;    }

private:
    JumpMatrix                  m_jumpMatrix;   ///< The jump matrix for the primal problem
    SparseMatrix                m_localMatrix;  ///< The overall matrix for the primal problem
    Matrix                      m_localRhs;     ///< The right-hand side for the primal problem
    std::vector<SparseMatrix>   m_primalBases;  ///< The bases for the primal dofs on the patches
    std::vector<OpPtr>          m_embeddings;   ///< For each patch, the map \f$ \tilde u_k \f$ to \f$ u_k \f$
    bool                        m_eliminatePointwiseConstraints; //TODO
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPrimalSystem.hpp)
#endif
