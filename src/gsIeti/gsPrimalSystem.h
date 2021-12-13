/** @file gsPrimalSystem.h

    @brief This class represents the primal system for a IETI-DP algorithm

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


/** @brief   This class represents the primal system for a IETI-DP algorithm
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
 *  It is assumed that the class receives the matrix \f$ A_k \f$ and the vectors
 *  \f$ c_k^{(n)} \f$ that make up the constraints matrix \f$ C_k \f$ such
 *  that the local saddle point formulation reads as follows:
 *
 *  \f[
 *       \tilde A_k
 *       =
 *       \begin{pmatrix}
 *            A_k   &    C_k^\top \\
 *            C_k   &
 *       \end{pmatrix}
 *  \f]
 *
 *  Using \ref setEliminatePointwiseConstraints, it is possible to instruct the
 *  class to eliminate pointwise constraints. A constraint is pointwise iff the vector
 *  \f$ c_k^{(n)} \f$ only has one non-zero entry accordning to the sparsity pattern
 *  of the sparse vector data structure. This is typically the case for vertex values.
 *
 *  Usually, it will be the best only to call \ref handleConstraints, which modifies
 *  the local system matrix, the jump matrix and the right-hand side such that they
 *  can then handed over to the \a gsIetiSystem. That function simultainously collects
 *  the contributions for the primal problem, which is handed over to \a gsIetiSystem
 *  as last subspace.
 *
 *  The member \ref handleConstraints is based on the static members
 *  \ref incorporateConstraints, \ref primalBasis and the member \ref addContribution.
 *  \ref incorporateConstraints sets up the matrix \f$ \tilde A_k \f$,
 *  \f$ R_c \f$ (right-hand-side for basis), \f$ P_c \f$ (embedding for basis) and
 *  \f$ P_0 \f$ (local embedding). The member \ref primalBasis constructs the basis
 *  for the primal space by the principle of energy minimization by calculating
 *  \f$ \Psi = P_c^\top \tilde{A}_k^{-1} R_c \f$. The member \ref addContribution adds
 *  the local contributions to the primal problem. Finally, the matrix \f$ P_0 \f$
 *  is used to modify the jump matrix and the local rhs.
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
    /// @param[in]  eliminatePointwiseConstraints  True iff pointwise constraints (typically
    ///                                   vertex values are eliminated (and not treated as local
    ///                                   saddle point problems). A constraint is pointwise iff
    ///                                   \f$ C_k \f$ has one non-zero entry.
    /// @param[in]  localMatrix           Local matrix \f$ A_k \f$
    /// @param[out] modifiedLocalMatrix   Local matrix \f$ \tilde A_k \f$
    /// @param[out] localEmbedding        Embedding \f$ P_0 \f$. This matrix is used by to modify
    ///                                   the right-hand side and the jump matrix
    ///                                   (cf. \ref handleConstraints)
    /// @param[out] embeddingForBasis     Embedding \f$ P_c \f$ for energy-minimizing basis
    ///                                   (cf. \ref primalBasis)
    /// @param[out] rhsForBasis           Right-hand side \f$ R_c \f$ for energy-minimizing basis
    ///                                   (cf. \ref primalBasis)
    ///
    /// Iff the constraint is not handled by elimination, the output is:
    /// \f[
    ///      \tilde{A}_k =
    ///      \begin{pmatrix}
    ///           A_k   &  C_k^\top \\  C_k   &  0
    ///      \end{pmatrix}
    ///      \quad\text{and}\quad
    ///      P_0 = P_c =
    ///      \begin{pmatrix}
    ///           I   &   0
    ///      \end{pmatrix}
    ///      \quad\text{and}\quad
    ///      R_c =
    ///      \begin{pmatrix}
    ///           0   \\   I
    ///      \end{pmatrix}.
    /// \f]
    ///
    /// Iff the constraint is handled by elimination and if we assume that the first block
    /// component of
    /// \f[
    ///      A_k =
    ///      \begin{pmatrix}
    ///           A_{11} & A_{12} \\  A_{21} & A_{22}
    ///      \end{pmatrix}
    /// \f]
    /// is to be eliminated, then we have:
    /// \f[
    ///      \tilde{A}_k =
    ///      \begin{pmatrix}
    ///           I \\ & A_{22}
    ///      \end{pmatrix}
    ///      \quad\text{and}\quad
    ///      P_0 =
    ///      \begin{pmatrix}
    ///           0   \\ &  I
    ///      \end{pmatrix}
    ///      \quad\text{and}\quad
    ///      P_c =
    ///      \begin{pmatrix}
    ///           I   \\ &  I
    ///      \end{pmatrix}
    ///      \quad\text{and}\quad
    ///      R_c =
    ///      \begin{pmatrix}
    ///           I   \\ -A_{21}
    ///      \end{pmatrix}.
    /// \f]
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
    /// @param  embeddingForBasis       Embedding \f$ P_c \f$ as provided by
    ///                                 \ref incorporateConstraints
    /// @param  rhsForBasis             Embedding \f$ R_c \f$ as provided by
    ///                                 \ref incorporateConstraints
    /// @param  primalDofIndices        Vector, that contains for every primal
    ///                                 constraint the index of the respective
    ///                                 primal dof
    /// @param  nPrimalDofs             The total number of primal dofs (\ref nPrimalDofs())
    ///
    /// The local basis is given by \f$ \Psi = P_c^\top \tilde{A}_k^{-1} R_c \f$.
    ///
    /// @returns a version of \f$ \Psi \f$ where the column indices are changed according
    ///          to primalDofIndices.
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
    ///    gsSparseMatrix<T> modifiedLocalMatrix, localEmbedding, embeddingForBasis;
    ///    gsMatrix<T> rhsForBasis;
    ///
    ///    gsPrimalSystem<T>::incorporateConstraints(
    ///        primalConstraints, this->eliminatePointwiseConstraints(), localMatrix,
    ///        modifiedLocalMatrix, localEmbedding, embeddingForBasis, rhsForBasis);
    ///
    ///    this->addContribution(
    ///        jumpMatrix, localMatrix, localRhs,
    ///        gsPrimalSystem<T>::primalBasis(
    ///            makeSparseLUSolver(modifiedLocalMatrix),
    ///            embeddingForBasis, rhsForBasis, primalDofIndices, this->nPrimalDofs()
    ///        )
    ///    );
    ///
    ///    localMatrix  = give(modifiedLocalMatrix);
    ///    localRhs     = localEmbedding * localRhs;
    ///    jumpMatrix   = jumpMatrix * localEmbedding.transpose();
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

    /// Returns the jump matrix for the primal problem
    JumpMatrix&                           jumpMatrix()        { return m_jumpMatrix;                    }
    const JumpMatrix&                     jumpMatrix() const  { return m_jumpMatrix;                    }

    /// Returns the local stiffness matrix for the primal problem
    SparseMatrix&                         localMatrix()       { return m_localMatrix;                   }
    const SparseMatrix&                   localMatrix() const { return m_localMatrix;                   }

    /// Returns the right-hand-side for the primal problem
    Matrix&                               localRhs()          { return m_localRhs;                      }
    const Matrix&                         localRhs() const    { return m_localRhs;                      }

    /// Returns the size of the primal problem (number of primal dofs)
    index_t nPrimalDofs() const                               { return m_localMatrix.rows();            }

    /// Returns true iff \ref handleConstraints will eliminate pointwise constraints (typically vertex values)
    index_t eliminatePointwiseConstraints() const             { return m_eliminatePointwiseConstraints; }

    /// Iff true, \ref handleConstraints will eliminate pointwise constraints (typically vertex values)
    void setEliminatePointwiseConstraints(bool v)             { m_eliminatePointwiseConstraints = v;    }

private:
    JumpMatrix                  m_jumpMatrix;   ///< The jump matrix for the primal problem
    SparseMatrix                m_localMatrix;  ///< The overall matrix for the primal problem
    Matrix                      m_localRhs;     ///< The right-hand side for the primal problem
    std::vector<SparseMatrix>   m_primalBases;  ///< The bases for the primal dofs on the patches
    std::vector<OpPtr>          m_embeddings;   ///< For each patch, the map \f$ \tilde u_k \f$ to \f$ u_k \f$
    bool                        m_eliminatePointwiseConstraints; ///< \ref handleConstraints will eliminate pointwise constraints
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPrimalSystem.hpp)
#endif
