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
    : m_localMatrix(nPrimalDofs,nPrimalDofs), m_eliminatePointwiseConstraints(false)
{
    this->m_localRhs.setZero(nPrimalDofs,1);
}

template <class T>
void gsPrimalSystem<T>::incorporateConstraints(
        const std::vector<SparseVector>& primalConstraints,
        bool eliminatePointwiseConstraints,
        const SparseMatrix& localMatrix,
        SparseMatrix& modifiedLocalMatrix,
        SparseMatrix& localEmbedding,
        SparseMatrix& embeddingForBasis,
        Matrix& rhsForBasis
    )
{
    const index_t localDofs = localMatrix.rows();
    const index_t nrPrimalConstraints = primalConstraints.size();

    // Which dofs should we eliminate?
    index_t nElimDofs = 0;
    std::vector<bool> eliminatedDof(localDofs, false);
    std::vector<bool> eliminatedConstaint(nrPrimalConstraints, false);

    {
        gsSparseEntries<T> seModifiedLocalMatrix;
        seModifiedLocalMatrix.reserve( localMatrix.nonZeros() + nrPrimalConstraints * localDofs );

        if (eliminatePointwiseConstraints)
        {
            for (index_t i=0; i!=nrPrimalConstraints; ++i)
            {
                if (primalConstraints[i].nonZeros() == 1)
                {
                    typename SparseVector::InnerIterator it(primalConstraints[i]);
                    // mark dof and constraint as to be eliminated
                    eliminatedDof[it.row()] = true;
                    eliminatedConstaint[i] = true;
                    // set diagonal matrix of new matrix
                    seModifiedLocalMatrix.add(it.row(),it.row(),(T)1);
                    ++nElimDofs;
                }
             }
        }

        // Copy entries of local matrix if not eliminated
        for (index_t i=0; i!=localDofs; ++i)
            for (typename SparseMatrix::InnerIterator it(localMatrix,i); it; ++it)
                if ( !eliminatedDof[it.row()] && !eliminatedDof[it.col()]  )
                    seModifiedLocalMatrix.add(it.row(), it.col(), it.value());

        // Make saddle point for non-eliminated primal constraints
        for (index_t i=0, j=0; i!=nrPrimalConstraints; ++i)
        {
            if ( !eliminatedConstaint[i] )
            {
                for (typename SparseVector::InnerIterator it(primalConstraints[i]); it; ++it)
                {
                    if ( !eliminatedDof[it.row()] )
                    {
                        seModifiedLocalMatrix.add(it.row(), localDofs+j, it.value());
                        seModifiedLocalMatrix.add(localDofs+j, it.row(), it.value());
                    }
                }
                ++j;
            }
        }

        modifiedLocalMatrix.clear();
        modifiedLocalMatrix.resize(localDofs+nrPrimalConstraints-nElimDofs, localDofs+nrPrimalConstraints-nElimDofs);
        modifiedLocalMatrix.setFrom(seModifiedLocalMatrix);
        modifiedLocalMatrix.makeCompressed();
    }

    // Compute the embedding matrices
    {
        localEmbedding.clear();
        localEmbedding.resize(localDofs+nrPrimalConstraints-nElimDofs,localDofs);
        embeddingForBasis.clear();
        embeddingForBasis.resize(localDofs+nrPrimalConstraints-nElimDofs,localDofs);
        gsSparseEntries<T> seLocalEmbedding, seEmbeddingForBasis;
        seLocalEmbedding.reserve(localDofs-nElimDofs);
        seEmbeddingForBasis.reserve(localDofs);
        for (index_t i=0; i!=localDofs; ++i)
        {
            if ( !eliminatedDof[i] )
                seLocalEmbedding.add(i,i,(T)1);
            seEmbeddingForBasis.add(i,i,(T)1);
        }
        localEmbedding.setFrom(seLocalEmbedding);
        localEmbedding.makeCompressed();
        embeddingForBasis.setFrom(seEmbeddingForBasis);
        embeddingForBasis.makeCompressed();
    }

    // Compute rhs for basis
    rhsForBasis.setZero(localDofs+nrPrimalConstraints-nElimDofs,nrPrimalConstraints);
    for (index_t i=0, j=0; i!=nrPrimalConstraints; ++i)
    {
        if ( eliminatedConstaint[i] )
        {
            typename SparseVector::InnerIterator it(primalConstraints[i]);
            const index_t idx = it.row();
            for (typename SparseMatrix::InnerIterator it2(localMatrix, idx); it2; ++it2)
            {
                if ( !eliminatedDof[it2.row()] )
                    rhsForBasis(it2.row(),i) = - it2.value();
            }
            rhsForBasis(idx,i) = (T)1;
        }
        else
        {
            rhsForBasis(localDofs+j,i) = (T)1;
            ++j;
        }
    }
}

template <class T>
typename gsPrimalSystem<T>::SparseMatrix
gsPrimalSystem<T>::primalBasis(
        OpPtr localSaddlePointSolver,
        const SparseMatrix& embeddingForBasis,
        const Matrix& rhsForBasis,
        const std::vector<index_t>& primalDofIndices,
        index_t nPrimalDofs
    )
{
    const index_t nrPrimalConstraints = primalDofIndices.size();
    const index_t localDofs = embeddingForBasis.cols();

    GISMO_ASSERT( nrPrimalConstraints<=nPrimalDofs, "gsPrimalSystem::primalBasis: "
        "There are more local constraints that there are constraints in total. "
        "Forgot to call gsPrimalSystem::init()?" );

    SparseMatrix result( localDofs, nPrimalDofs );

    if (nPrimalDofs==0) return result;

    Matrix tmp;
    localSaddlePointSolver->apply(rhsForBasis, tmp);
    Matrix basis = embeddingForBasis.transpose() * tmp;

    gsSparseEntries<T> se_result;
    se_result.reserve(localDofs*nrPrimalConstraints);
    for (index_t i=0; i<localDofs; ++i)
        for (index_t j=0; j<nrPrimalConstraints; ++j)
        {
            GISMO_ASSERT( primalDofIndices[j]>=0 && primalDofIndices[j]<nPrimalDofs,
                "gsPrimalSystem::primalBasis: Invalid index.");
            se_result.add(i,primalDofIndices[j],basis(i,j));
        }

    result.setFrom(se_result);
    return result;
}

template <class T>
void gsPrimalSystem<T>::addContribution(
        const JumpMatrix& jumpMatrix,
        const SparseMatrix& localMatrix,
        const Matrix& localRhs,
        SparseMatrix primalBasis,
        OpPtr embedding
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
    m_embeddings.push_back(give(embedding));
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
    SparseMatrix modifiedLocalMatrix, localEmbedding, embeddingForBasis;
    Matrix rhsForBasis;

    incorporateConstraints(primalConstraints,eliminatePointwiseConstraints(),
        localMatrix,
        modifiedLocalMatrix,localEmbedding,embeddingForBasis,rhsForBasis);

    addContribution(
        jumpMatrix, localMatrix, localRhs,
        primalBasis(
            makeSparseLUSolver(modifiedLocalMatrix),
            embeddingForBasis, rhsForBasis, primalDofIndices, nPrimalDofs()
        )
    );

    localMatrix  = give(modifiedLocalMatrix);
    localRhs     = localEmbedding * localRhs;
    jumpMatrix   = jumpMatrix * localEmbedding.transpose();
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

        if (m_embeddings[i])
        {
            Matrix tmp;
            m_embeddings[i]->apply(sol[i],tmp);
            sol[i].swap(tmp);
        }
        else
            sol[i].conservativeResize( this->m_primalBases[i].rows(), Eigen::NoChange );
        sol[i] += this->m_primalBases[i] * sol.back();
    }

    sol.pop_back();

    return sol;
}


} // namespace gismo
