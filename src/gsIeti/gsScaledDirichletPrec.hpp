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
#include <gsCore/gsMultiBasis.h>
#include <gsSolver/gsBlockOp.h>

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
        const index_t sz = m_localSchurOps[k]->rows();
        //gsMatrix<T> & sc = m_localScaling[k];
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
void gsScaledDirichletPrec<T>::setupDeluxeScaling(const std::vector<boundaryInterface>& ifaces) // TODO: Consider entry != 0 for the corners
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

    for (std::vector<boundaryInterface>::const_iterator it = ifaces.begin(); it != ifaces.end(); it++)
    {
        patchSide side1 = it->first();
        patchSide side2 = it->second();

        SparseMatrix localSchurRestriction1, localSchurRestriction2;
        gsSparseEntries<T> se1_restriction, se2_restriction;
        JumpMatrix & jm1 = *(m_jumpMatrices[side1.patch]);
        JumpMatrix & jm2 = *(m_jumpMatrices[side2.patch]);
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

        localSchurRestriction1.resize(r, m_localSchurOps[side1.patch]->cols());
        localSchurRestriction2.resize(r, m_localSchurOps[side2.patch]->cols());
        localSchurRestriction1.setFrom(se1_restriction);
        localSchurRestriction2.setFrom(se2_restriction);

        // Restrict the Schur complements S^k and S^l to the interface to obtain \Sum_i R_{iF} S^i R_{iF}^T,
        // where R_{iF} denotes the restriction to the interface F
        typename gsAdditiveOp<>::Ptr ifaceRestriction = gsAdditiveOp<>::make(std::vector<gsSparseMatrix<T, RowMajor>>{localSchurRestriction1, localSchurRestriction2},
                                                                std::vector<OpPtr>{m_localSchurOps[side1.patch], m_localSchurOps[side2.patch]}
                                                                );

        // Build the Schur complement matrices restricted to the interface, i.e., S^k_F + S^l_F
        gsMatrix<T> result;
        ifaceRestriction->toMatrix(result); // TODO:Might be costly if matrix dofs on faces or edges grow, better solution??

        // Build the deluxe scaling matrices (S^k_F + S^l_F)^{-1}S^k_F and (S^k_F + S^l_F)^{-1}S^l_F, respectively
        // They sum up to the identity!
        typename gsLinearOperator<T>::Ptr SkFSlFInv = makePartialPivLUSolver(result);
        typename gsProductOp<>::Ptr deluxeMatk = gsProductOp<>::make( gsProductOp<>::make(
                makeMatrixOp(ifaceRestriction->getTransfers()[0]->transpose()), ifaceRestriction->getOps()[0], makeMatrixOp(ifaceRestriction->getTransfers()[0])
                        ),
                                                                SkFSlFInv
                                                                );
        typename gsProductOp<>::Ptr deluxeMatl = gsProductOp<>::make( gsProductOp<>::make(
                makeMatrixOp(ifaceRestriction->getTransfers()[1]->transpose()), ifaceRestriction->getOps()[1], makeMatrixOp(ifaceRestriction->getTransfers()[1])
                ),
                                                                      SkFSlFInv
                                                                 );

        // Build also the transposed deluxe scaling matrices
        typename gsProductOp<>::Ptr deluxeMatTransk = gsProductOp<>::make(SkFSlFInv, gsProductOp<>::make(
                makeMatrixOp(ifaceRestriction->getTransfers()[0]->transpose()), ifaceRestriction->getOps()[0], makeMatrixOp(ifaceRestriction->getTransfers()[0])
                                                                      )
        );
        typename gsProductOp<>::Ptr deluxeMatTransl = gsProductOp<>::make(SkFSlFInv, gsProductOp<>::make(
                makeMatrixOp(ifaceRestriction->getTransfers()[1]->transpose()), ifaceRestriction->getOps()[1], makeMatrixOp(ifaceRestriction->getTransfers()[1])
                                                                      )
        );

        // Add the deluxe scaling matrices to a sum operator for each patch
        scalingOps[side1.patch]->addOperator(gsAdditiveOp<>::make(std::vector<gsSparseMatrix<T, RowMajor>>{ifaceRestriction->getTransfers()[0]->transpose()}, std::vector<OpPtr>{deluxeMatk}));
        scalingOps[side2.patch]->addOperator(gsAdditiveOp<>::make(std::vector<gsSparseMatrix<T, RowMajor>>{ifaceRestriction->getTransfers()[1]->transpose()}, std::vector<OpPtr>{deluxeMatl}));

        scalingTransOps[side1.patch]->addOperator(gsAdditiveOp<>::make(std::vector<gsSparseMatrix<T, RowMajor>>{ifaceRestriction->getTransfers()[0]->transpose()}, std::vector<OpPtr>{deluxeMatTransk}));
        scalingTransOps[side2.patch]->addOperator(gsAdditiveOp<>::make(std::vector<gsSparseMatrix<T, RowMajor>>{ifaceRestriction->getTransfers()[1]->transpose()}, std::vector<OpPtr>{deluxeMatTransl}));
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
