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
    typedef gsLinearOperator<T>               Op;
    typedef memory::shared_ptr<Op>            OpPtr;
    typedef gsSparseMatrix<T>                 SparseMatrix;
    typedef gsSparseMatrix<T,RowMajor>        JumpMatrix;
    typedef gsSparseVector<T>                 SparseVector;
    typedef gsMatrix<T>                       Matrix;
public:

    /// @brief Initialize the object
    ///
    /// @param nPrimalDofs                Number of primal constraints in total
    void init(index_t nPrimalDofs);

    /// @brief Incorporates the given constraints in the local system
    ///
    /// @param[in]     primalConstraints  Primal constraints; the given vectors
    ///                                   make up the matrix C
    /// @param[in,out] jumpMatrix         Jump matrix B
    /// @param[in,out] localMatrix        Local matrix A
    /// @param[in,out] localRhs           Local right-hand-side f
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
    ///      \underline{\tilde f}_k =
    ///      \begin{pmatrix}
    ///           \underline f_k \\ 0
    ///      \end{pmatrix}
    ///  \f]
    ///
    static void incorporateConstraints(
        const std::vector<SparseVector>& primalConstraints,
        JumpMatrix& jumpMatrix,
        SparseMatrix& localMatrix,
        Matrix& localRhs
    );

    /// @brief Returns the matrix representation of the energy minimizing primal basis
    ///
    /// @param  localSaddlePointSolver  Solver that realizes \f$ \tilde{A}^{-1} \f$
    /// @param  primalDofIndices        Vector, that contains for every primal
    ///                                 constraint the index of the respective
    ///                                 primal dof
    /// @param  nPrimalDofs             The total number of primal dofs
    ///
    static gsSparseMatrix<T> primalBasis(
        OpPtr localSaddlePointSolver,
        const std::vector<index_t>& primalDofIndices,
        index_t nPrimalDofs
    );

    /// @brief Adds contributions for a patch to the data hold in the class
    ///
    /// @param[in,out]  jumpMatrix    Jump matrix B_k
    /// @param[in,out]  localMatrix   Local matrix A_k
    /// @param[in,out]  localRhs      Local right-hand-side f_k
    /// @param[in]      primalBasis   Matrix representation of the primal basis
    ///
    void addContribution(
        const JumpMatrix& jumpMatrix,
        const SparseMatrix& localMatrix,
        const Matrix& localRhs,
        SparseMatrix primalBasis
    );

    /// @brief Convenience function for handling the primal constraints

    /// @param[in]     primalConstraints  Primal constraints; the given vectors
    ///                                   make up the matrix C_k
    /// @param[in]     primalDofIndices   Vector, that contains for every primal
    ///                                   constraint the index of the respective
    ///                                   primal dof
    /// @param[in,out] jumpMatrix         Jump matrix B_k
    /// @param[in,out] localMatrix        Local matrix A_k
    /// @param[in,out] localRhs           Local right-hand-side f_k
    ///
    void handleConstraints(
        const std::vector<SparseVector>& primalConstraints,
        const std::vector<index_t>& primalDofIndices,
        JumpMatrix& jumpMatrix,
        SparseMatrix& localMatrix,
        Matrix& localRhs
    )
    {
        incorporateConstraints(primalConstraints,jumpMatrix,localMatrix,localRhs);
        addContribution(jumpMatrix,localMatrix,localRhs,
            primalBasis(makeSparseLUSolver(localMatrix),primalDofIndices,nPrimalDofs())
        );
    }

    /// @brief  Distributes the given solution for K+1 subdomains to the K patches
    ///
    /// @param    sol   The solution, first for the K patches, followed by the
    ///                 contribution for the primal dofs
    /// @returns        The solution for the K patches
    std::vector<Matrix> distributePrimalSolution( std::vector<Matrix> sol );

    /// @brief Returns the jump matrix for the primal problem
    JumpMatrix&                           jumpMatrix()        { return m_jumpMatrix;         }
    const JumpMatrix&                     jumpMatrix() const  { return m_jumpMatrix;         }

    /// @brief Returns the local matrix for the primal problem
    SparseMatrix&                         localMatrix()       { return m_localMatrix;        }
    const SparseMatrix&                   localMatrix() const { return m_localMatrix;        }

    /// @brief Returns the right-hand-side for the primal problem
    Matrix&                               localRhs()          { return m_localRhs;           }
    const Matrix&                         localRhs() const    { return m_localRhs;           }

    /// @brief Returns the size of the primal problem (number of primal dofs)
    index_t nPrimalDofs() const                               { return m_localMatrix.rows(); }

private:
    JumpMatrix                  m_jumpMatrix;   ///< The jump matrix for the primal problem
    SparseMatrix                m_localMatrix;  ///< The overall matrix for the primal problem
    Matrix                      m_localRhs;     ///< The right-hand side for the primal problem
    std::vector<SparseMatrix>   m_primalBases;  ///< The bases for the primal dofs on the patches
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPrimalSystem.hpp)
#endif
