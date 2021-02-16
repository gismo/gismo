/** @file gsPrimalSystem.hpp

    @brief This class represents the primal system for a IETI-dp algorithm

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

namespace gismo
{

template <class T>
gsPrimalSystem<T>::gsPrimalSystem(index_t nPrimalDofs)
{
    this->m_localMatrix.resize(nPrimalDofs,nPrimalDofs);
    this->m_localRhs.setZero(nPrimalDofs,1);
}

template <class T>
void gsPrimalSystem<T>::incorporateConstraints(
        const std::vector<SparseVector>& primalConstraints,
        const JumpMatrix& jumpMatrix,
        const SparseMatrix& localMatrix,
        const Matrix& localRhs,
        JumpMatrix& modifiedJumpMatrix,
        SparseMatrix& modifiedLocalMatrix,
        Matrix& modifiedLocalRhs,
        Matrix& rhsForBasis
    )
{
    const index_t localDofs = localMatrix.rows();
    const index_t nrPrimalConstraints = primalConstraints.size();
    if (nrPrimalConstraints==0) return;

    gsSparseEntries<T> seLocalMatrix;
    seLocalMatrix.reserve( localMatrix.nonZeros() + nrPrimalConstraints * localDofs );

    for (index_t i=0; i<localDofs; ++i)
        for (typename SparseMatrix::InnerIterator it(localMatrix,i); it; ++it)
            seLocalMatrix.add(it.row(), it.col(), it.value());

    for (index_t i=0; i<nrPrimalConstraints; ++i)
        for (typename SparseVector::InnerIterator it(primalConstraints[i]); it; ++it)
        {
            seLocalMatrix.add(it.row(), localDofs+i, it.value());
            seLocalMatrix.add(localDofs+i, it.row(), it.value());
        }

    modifiedLocalMatrix.clear();
    modifiedLocalMatrix.resize(localDofs+nrPrimalConstraints, localDofs+nrPrimalConstraints);
    modifiedLocalMatrix.setFrom(seLocalMatrix);
    modifiedLocalMatrix.makeCompressed();

    GISMO_ASSERT( localRhs.rows() == localDofs,
        "gsPrimalSystem::incorporateConstraint: Right-hand side does not have the expected number of columns;");

    modifiedLocalRhs.setZero(localDofs+nrPrimalConstraints, localRhs.cols());
    modifiedLocalRhs.topRows(localDofs) = localRhs;

    GISMO_ASSERT( jumpMatrix.cols() == localDofs,
        "gsPrimalSystem::incorporateConstraint: Jump matrix does not have the expected number of columns;");

    modifiedJumpMatrix = jumpMatrix;
    modifiedJumpMatrix.conservativeResize(jumpMatrix.rows(), localDofs+nrPrimalConstraints);

    rhsForBasis.setZero(localDofs+nrPrimalConstraints,nrPrimalConstraints);
    for (index_t i=0; i<nrPrimalConstraints; ++i)
        rhsForBasis(localDofs+i,i) = 1;

}

template <class T>
typename gsPrimalSystem<T>::SparseMatrix
gsPrimalSystem<T>::primalBasis(
        OpPtr localSaddlePointSolver,
        const Matrix& rhsForBasis,
        const std::vector<index_t>& primalDofIndices,
        index_t nPrimalDofs
    )
{
    const index_t nrPrimalConstraints = primalDofIndices.size();

    GISMO_ASSERT( nrPrimalConstraints<=nPrimalDofs, "gsPrimalSystem::primalBasis: "
        "There are more local constraints that there are constraints in total. "
        "Forgot to call gsPrimalSystem::init()?" );

    const index_t localDofs = localSaddlePointSolver->rows() - nrPrimalConstraints;

    SparseMatrix result( localDofs, nPrimalDofs );

    if (nPrimalDofs==0) return result;

    Matrix tmp;
    localSaddlePointSolver->apply(rhsForBasis, tmp);

    gsSparseEntries<T> se_result;
    se_result.reserve(localDofs*nrPrimalConstraints);
    for (index_t i=0; i<localDofs; ++i)
        for (index_t j=0; j<nrPrimalConstraints; ++j)
        {
            GISMO_ASSERT( primalDofIndices[j]>=0 && primalDofIndices[j]<nPrimalDofs,
                "gsPrimalSystem::primalBasis: Invalid index.");
            se_result.add(i,primalDofIndices[j],tmp(i,j));
        }

    result.setFrom(se_result);

    return result;
}

template <class T>
void gsPrimalSystem<T>::addContribution(
        const JumpMatrix& jumpMatrix,
        const SparseMatrix& localMatrix,
        const Matrix& localRhs,
        SparseMatrix primalBasis
    )
{
    if (m_jumpMatrix.rows()==0&&m_jumpMatrix.cols()==0)
        m_jumpMatrix.resize(jumpMatrix.rows(),nPrimalDofs());

    GISMO_ASSERT( primalBasis.cols() == m_jumpMatrix.cols(),
        "gsPrimalSystem::incorporate: The given problem size does not match the stored primal problem size. "
        "Forgot to call gsPrimalSystem::init()?" );

    m_localMatrix     += primalBasis.transpose() * localMatrix * primalBasis;
    m_localRhs        += primalBasis.transpose() * localRhs;
    m_jumpMatrix      += JumpMatrix(jumpMatrix * primalBasis);
    m_primalBases.push_back(give(primalBasis));
}

template <class T>
void gsPrimalSystem<T>::handleConstraints(
        const std::vector<SparseVector>& primalConstraints,
        const std::vector<index_t>& primalDofIndices,
        JumpMatrix& jumpMatrix,
        SparseMatrix& localMatrix,
        Matrix& localRhs
    )
{
    JumpMatrix modifiedJumpMatrix;
    SparseMatrix modifiedLocalMatrix;
    Matrix modifiedLocalRhs, rhsForBasis;

    incorporateConstraints(primalConstraints,jumpMatrix,localMatrix,localRhs,
        modifiedJumpMatrix,modifiedLocalMatrix,modifiedLocalRhs,rhsForBasis);

    addContribution(
        jumpMatrix, localMatrix, localRhs,
        primalBasis( makeSparseLUSolver(modifiedLocalMatrix),
            rhsForBasis, primalDofIndices, nPrimalDofs()
        )
    );

    jumpMatrix   = give(modifiedJumpMatrix);
    localMatrix  = give(modifiedLocalMatrix);
    localRhs     = give(modifiedLocalRhs);
}

template <class T>
std::vector<typename gsPrimalSystem<T>::Matrix>
gsPrimalSystem<T>::distributePrimalSolution( std::vector<Matrix> sol )
{
    const index_t sz = this->m_primalBases.size();

    // If the primal problem is empty, there might just not be any primal subdomain
    if (static_cast<index_t>(sol.size())==sz && this->m_jumpMatrix.cols()==0)
        return sol;

    GISMO_ASSERT(static_cast<index_t>(sol.size())==sz+1, "gsPrimalSystem::distributePrimalSolution "
        "expects that there is one more subdomain that patches.");

    for (index_t i=0; i<sz; ++i)
    {
        GISMO_ASSERT( sol[i].rows() >= this->m_primalBases[i].rows()
            && this->m_primalBases[i].cols() == sol.back().rows()
            && sol.back().cols() == sol[i].cols(),
            "gsPrimalSystem::distributePrimalSolution: Dimensions do not agree: "
            << sol[i].rows() << ">=" << this->m_primalBases[i].rows() << "&&"
            << this->m_primalBases[i].cols() << "==" << sol.back().rows() << "&&"
            << sol.back().cols() << "==" << sol[i].cols() << " ( i=" << i << "). "
            << "This method assumes the primal subspace to be the last one." );

        sol[i].conservativeResize( this->m_primalBases[i].rows(), Eigen::NoChange );
        sol[i] += this->m_primalBases[i] * sol.back();
    }

    sol.pop_back();

    return sol;
}


} // namespace gismo
