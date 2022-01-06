/** @file gsScaledDirichletPrec.hpp

    @brief This class represents the sclaed Dirichlet preconditioner.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsAdditiveOp.h>

namespace gismo
{

template <class T>
gsSortedVector<index_t>
gsScaledDirichletPrec<T>::skeletonDofs( const JumpMatrix& jm )
{
    gsSortedVector<index_t> result;
    for (index_t i=0; i<jm.outerSize(); ++i)
        for (typename JumpMatrix::InnerIterator it(jm, i); it; ++it)
            result.push_sorted_unique(it.col());
    return result;
}

template <class T>
typename gsScaledDirichletPrec<T>::JumpMatrix
gsScaledDirichletPrec<T>::restrictJumpMatrix( const JumpMatrix& jm, const std::vector<index_t> dofs )
{
    gsVector<index_t> reverse;
    reverse.setZero( jm.cols() );
    const index_t sz = dofs.size();
    for (index_t i=0; i<sz; ++i)
        reverse[dofs[i]] = i+1;

    gsSparseEntries<T> triplets;
    triplets.reserve(jm.nonZeros());
    for (index_t i=0; i<jm.outerSize(); ++i)
        for (typename JumpMatrix::InnerIterator it(jm, i); it; ++it)
            if (reverse[it.col()] > 0)
            {
                triplets.add(
                    it.row(),
                    reverse[it.col()]-1,
                    it.value()
                );
            }

    JumpMatrix result(jm.rows(), dofs.size());
    result.setFrom(triplets);
    return result;
}

template <class T>
typename gsScaledDirichletPrec<T>::Blocks
gsScaledDirichletPrec<T>::matrixBlocks( const SparseMatrix& mat, const std::vector<index_t> dofs )
{
    gsVector<index_t> reverse;
    reverse.setZero( mat.cols() );
    const index_t sz = dofs.size();
    for (index_t i=0; i<sz; ++i)
        reverse[dofs[i]] = i+1;
    index_t j=0;
    for (index_t i=0; i<mat.cols(); ++i)
        if (reverse[i]==0)
            reverse[i] = --j;

    gsSparseEntries<T> se_A00, se_A01, se_A10, se_A11;
    se_A00.reserve( 2 * mat.nonZeros() * dofs.size() / mat.rows() );
    se_A01.reserve( 2 * mat.nonZeros() * dofs.size() / mat.rows() );
    se_A10.reserve( 2 * mat.nonZeros() * dofs.size() / mat.rows() );
    se_A11.reserve( mat.nonZeros() );
    for (index_t i=0; i<mat.outerSize(); ++i)
        for (typename SparseMatrix::InnerIterator it(mat, i); it; ++it)
        {
            if (reverse[it.row()] > 0 && reverse[it.col()] > 0)
            {
                se_A00.add(
                    reverse[it.row()]-1,
                    reverse[it.col()]-1,
                    it.value()
                );
            }
            else if (reverse[it.row()] > 0 && reverse[it.col()] < 0)
            {
                se_A10.add(
                    reverse[it.row()]-1,
                    -reverse[it.col()]-1,
                    it.value()
                );
            }
            else if (reverse[it.row()] < 0 && reverse[it.col()] > 0)
            {
                se_A01.add(
                    -reverse[it.row()]-1,
                    reverse[it.col()]-1,
                    it.value()
                );
            }
            else //if (reverse[it.col()] < 0 && reverse[it.row()] < 0)
            {
                se_A11.add(
                    -reverse[it.row()]-1,
                    -reverse[it.col()]-1,
                    it.value()
                );
            }
        }

    Blocks result;

    result.A00.resize(           dofs.size(),            dofs.size());
    result.A00.setFrom(se_A00);
    result.A01.resize(mat.rows()-dofs.size(),            dofs.size());
    result.A01.setFrom(se_A01);
    result.A10.resize(           dofs.size(), mat.rows()-dofs.size());
    result.A10.setFrom(se_A10);
    result.A11.resize(mat.rows()-dofs.size(), mat.rows()-dofs.size());
    result.A11.setFrom(se_A11);

    return result;

}

template <class T>
typename gsScaledDirichletPrec<T>::OpPtr
gsScaledDirichletPrec<T>::schurComplement( Blocks matrixBlocks, OpPtr solver )
{
    matrixBlocks.A01 *= -1;
    return gsSumOp<T>::make(
        makeMatrixOp(matrixBlocks.A00.moveToPtr()),
        gsProductOp<T>::make(
            makeMatrixOp(matrixBlocks.A01.moveToPtr()),
            give(solver),
            makeMatrixOp(matrixBlocks.A10.moveToPtr())
        )
    );
}

template <class T>
void gsScaledDirichletPrec<T>::setupMultiplicityScaling()
{
    const index_t pnr = m_jumpMatrices.size();

    for (index_t k=0; k<pnr; ++k)
    {
        const index_t sz = m_jumpMatrices[k]->cols();
        gsSparseMatrix<T> sc;
        sc.resize(sz, sz);

        for (index_t i=0; i<sz; ++i)
          sc(i,i) = 1;

        JumpMatrix & jm = *(m_jumpMatrices[k]);

        for (index_t i=0; i<jm.outerSize(); ++i){
            for (typename JumpMatrix::InnerIterator it(jm, i); it; ++it)
            {
                const index_t c = it.col();
                sc(c,c) += 1;
            }
        }

        for (index_t i=0; i<sz; ++i)
            sc(i,i) = T(1)/sc(i,i);

        m_localScalingOps[k] = makeMatrixOp(sc.moveToPtr());
    }
    m_localScalingTransOps = m_localScalingOps;
}

template <class T>
void gsScaledDirichletPrec<T>::setupDeluxeScaling()
{
    // Assert that the jumpmatrices are registered
    GISMO_ASSERT(m_jumpMatrices.size() != 0, "gsScaledDirichletPrec<T>::setupDeluxeScaling(...), Forgot to call addSubdomain(...)?");
    std::vector<std::vector<SparseMatrix> > restrictions(m_jumpMatrices.size());

    // Construct restriction matrices to the faces, assume that there are no Lagrange multipliers for the corners
    std::vector<typename gsSumOp<>::Ptr > scalingOps(m_jumpMatrices.size()), scalingTransOps(m_jumpMatrices.size());
    for(index_t k = 0; k < (index_t)scalingOps.size(); ++k)
    {
        scalingOps[k] = gsSumOp<>::make();
        scalingTransOps[k] = gsSumOp<>::make();
    }

    // Always two patches are connected by one Lagrange multiplier
    gsMatrix<index_t> lambda2Patches(2, m_jumpMatrices[0]->rows()); // Assume that all jump matrices have the same rows
    for (index_t i=0; i<m_jumpMatrices[0]->outerSize(); ++i){
        for (size_t k = 0; k < m_jumpMatrices.size(); ++k)
        for (typename JumpMatrix::InnerIterator it(*m_jumpMatrices[k], i); it; ++it)
        {
            index_t c = it.col();
            if((*m_jumpMatrices[k])(i, c) == 1)
                lambda2Patches(0, i) = k;
            else if((*m_jumpMatrices[k])(i, c) == -1)
                lambda2Patches(1, i) = k;
        }
    }

    for (index_t col = 0; col < lambda2Patches.cols(); ++col)
    {
        index_t patch1 = lambda2Patches(0, col);
        index_t patch2 = lambda2Patches(1, col);

        SparseMatrix localSchurRestriction1, localSchurRestriction2;
        gsSparseEntries<T> se1_restriction, se2_restriction;
        JumpMatrix & jm1 = *(m_jumpMatrices[patch1]);
        JumpMatrix & jm2 = *(m_jumpMatrices[patch2]);
        se1_restriction.reserve( jm1.cols() );
        se2_restriction.reserve( jm2.cols() );

        index_t r = 0;
        for (index_t i=0; i<jm1.outerSize(); ++i){
            for (typename JumpMatrix::InnerIterator it(jm1, i); it; ++it)
            {
                const index_t c = it.col();
                GISMO_ASSERT(jm1.col(c).norm() == 1,"Deluxe scaling only allowed if there are no Lagrange multipliers for the corners!");
                for(typename JumpMatrix::InnerIterator it2(jm2, i); it2; ++it2)
                {
                    const index_t c2 = it2.col();
                    if(jm1(i, c) != 0 && jm2(i, c2) != 0) // the corresponding dofs are connected
                    {
                        se1_restriction.add(r, c, T(1));
                        se2_restriction.add(r, c2, T(1));
                        r++;
                        break;
                    }
                }
            }
        }

        localSchurRestriction1.resize(r, m_localSchurOps[patch1]->cols());
        localSchurRestriction2.resize(r, m_localSchurOps[patch2]->cols());
        localSchurRestriction1.setFrom(se1_restriction);
        localSchurRestriction2.setFrom(se2_restriction);

        // Restrict the A_{\Gamma\Gamma}^{k} and A_{\Gamma\Gamma}^{l} to the interface to obtain \Sum_i R_{iF} A_{\Gamma\Gamma}^{i} R_{iF}^T,
        // where R_{iF} denotes the restriction to the interface F
        SparseMatrix AFF = localSchurRestriction1 * m_blocks[patch1].A00 * localSchurRestriction1.transpose() + localSchurRestriction2 * m_blocks[patch2].A00 * localSchurRestriction2.transpose();
        SparseMatrix AFI1 = localSchurRestriction1 * m_blocks[patch1].A10;
        SparseMatrix AFI2 =  localSchurRestriction2 * m_blocks[patch2].A10;
        SparseMatrix AIF1 = m_blocks[patch1].A01 * localSchurRestriction1.transpose();
        SparseMatrix AIF2 = m_blocks[patch2].A01 * localSchurRestriction2.transpose();

        gsSparseEntries<T> se_blockMat;
        se_blockMat.reserve(m_blocks[patch1].A00.rows() + m_blocks[patch1].A11.rows() + m_blocks[patch1].A11.cols() + m_blocks[patch2].A11.cols());
        for (index_t i=0; i<AFF.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(AFF, i); it; ++it)
            {
                se_blockMat.add(i, it.col(), AFF(i, it.col()));
            }
        }

        for (index_t i=0; i<AFI1.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(AFI1, i); it; ++it)
            {
                se_blockMat.add(i, it.col() + AFF.cols(), AFI1(i, it.col()));
            }
        }

        for (index_t i=0; i<AFI2.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(AFI2, i); it; ++it)
            {
                se_blockMat.add(i, it.col() + AFF.cols() + AFI1.cols(), AFI2(i, it.col()));
            }
        }

        for (index_t i=0; i<AIF1.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(AIF1, i); it; ++it)
            {
                se_blockMat.add(i + AFF.rows(), it.col(), AIF1(i, it.col()));
            }
        }

        for (index_t i=0; i<AIF2.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(AIF2, i); it; ++it)
            {
                se_blockMat.add(i + AFF.rows() + AIF1.rows(), it.col(), AIF2(i, it.col()));
            }
        }

        for (index_t i=0; i<m_blocks[patch1].A11.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(m_blocks[patch1].A11, i); it; ++it)
            {
                se_blockMat.add(i + AFF.rows(), it.col() + AFF.cols(), m_blocks[patch1].A11(i, it.col()));
            }
        }

        for (index_t i=0; i<m_blocks[patch2].A11.outerSize(); ++i){
            for (typename SparseMatrix::InnerIterator it(m_blocks[patch2].A11, i); it; ++it)
            {
                se_blockMat.add(i + AFF.rows() + AIF1.rows(), it.col() + AFF.cols() + AFI1.cols(), m_blocks[patch2].A11(i, it.col()));
            }
        }

        gsSparseEntries<T> se_I; se_I.reserve(AFF.cols());
        for (index_t i=0; i<AFF.rows(); ++i) {
            se_I.add(i, i, T(1));
        }

        SparseMatrix blockMat(AFF.rows() + AIF1.rows() + AIF2.rows(), AFF.cols() + AFI1.cols() + AFI2.cols());
        blockMat.setFrom(se_blockMat);
        SparseMatrix I(AFF.cols(), AFF.rows() + AIF1.rows() + AIF2.rows());
        I.setFrom(se_I);

        // Restrict the Schur complements on the patches to the interfaces
        typename gsAdditiveOp<>::Ptr SFk = gsAdditiveOp<>::make(std::vector<SparseMatrix>{localSchurRestriction1},
                                                                std::vector<OpPtr>{m_localSchurOps[patch1]});

        typename gsAdditiveOp<>::Ptr SFl = gsAdditiveOp<>::make(std::vector<SparseMatrix>{localSchurRestriction2},
                                                                 std::vector<OpPtr>{m_localSchurOps[patch2]});

        // Build the inverse of the sum of Schur complements
        typename gsAdditiveOp<>::Ptr embedding = gsAdditiveOp<>::make(std::vector<SparseMatrix>{I},
                                                                 std::vector<OpPtr>{makeSparseCholeskySolver(blockMat)});

        // Build the scaling matrices
        typename gsProductOp<>::Ptr prod1 = gsProductOp<>::make(SFk, embedding);
        typename gsProductOp<>::Ptr prod2 = gsProductOp<>::make(SFl, embedding);
        typename gsProductOp<>::Ptr prod1trans = gsProductOp<>::make(embedding, SFk);
        typename gsProductOp<>::Ptr prod2trans = gsProductOp<>::make(embedding, SFl);

        scalingOps[patch1]->addOperator(gsAdditiveOp<>::make(std::vector<SparseMatrix>{SFk->getTransfers()[0]->transpose()}, std::vector<OpPtr>{prod1}));
        scalingOps[patch2]->addOperator(gsAdditiveOp<>::make(std::vector<SparseMatrix>{SFl->getTransfers()[0]->transpose()}, std::vector<OpPtr>{prod2}));

        scalingTransOps[patch1]->addOperator(gsAdditiveOp<>::make(std::vector<SparseMatrix>{SFk->getTransfers()[0]->transpose()}, std::vector<OpPtr>{prod1trans}));
        scalingTransOps[patch2]->addOperator(gsAdditiveOp<>::make(std::vector<SparseMatrix>{SFl->getTransfers()[0]->transpose()}, std::vector<OpPtr>{prod2trans}));
    }

    for(index_t k = 0; k < (index_t)scalingOps.size(); ++k)
    {
        m_localScalingOps[k] = std::move(scalingOps[k]);
        m_localScalingTransOps[k] = std::move(scalingTransOps[k]);
    }
}

template <class T>
typename gsScaledDirichletPrec<T>::OpPtr
gsScaledDirichletPrec<T>::preconditioner() const
{
    const index_t pnr = m_jumpMatrices.size();

    for (index_t i=0; i<pnr; ++i)
    {
        GISMO_ASSERT( m_localScalingOps[0] != nullptr && m_localScalingOps[pnr-1] != nullptr &&
                              m_localScalingTransOps[0] != nullptr && m_localScalingTransOps[pnr-1] != nullptr,
            "gsScaledDirichletPrec::preconditioner needs the local scaling operators given. "
            "Forgot to call setup...Scaling()?" );
    }

    typename gsAdditiveOp<T>::Ptr result = gsAdditiveOp<T>::make();

    for (index_t i=0; i<pnr; ++i)
    {
        typename gsProductOp<T>::Ptr local = gsProductOp<T>::make();
        local->addOperator(m_localScalingTransOps[i]);
        local->addOperator(m_localSchurOps[i]);
        local->addOperator(m_localScalingOps[i]);
        result->addOperator(m_jumpMatrices[i],local);
    }

    return result;
}


} // namespace gismo
