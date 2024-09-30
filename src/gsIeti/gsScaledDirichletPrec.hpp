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
        const index_t sz = m_localSchurOps[k]->rows();
        gsMatrix<T> & sc = m_localScaling[k];
        sc.resize(sz, 1);

        for (index_t i=0; i<sz; ++i)
          sc(i,0) = 1;

        JumpMatrix & jm = *(m_jumpMatrices[k]);

        for (index_t i=0; i<jm.outerSize(); ++i){
            for (typename JumpMatrix::InnerIterator it(jm, i); it; ++it)
            {
                const index_t c = it.col();
                sc(c,0) += 1;
            }
        }
    }
}

template <class T>
typename gsScaledDirichletPrec<T>::OpPtr
gsScaledDirichletPrec<T>::preconditioner() const
{
    const index_t pnr = m_jumpMatrices.size();

    for (index_t i=0; i<pnr; ++i)
    {
        GISMO_ASSERT( m_localScaling[i].rows() > 0 && m_localScaling[i].cols() == 1,
            "gsScaledDirichletPrec::preconditioner needs the local scaling matrices given. "
            "Forgot to call setupMultiplicityScaling()?" );
    }

    std::vector<OpPtr> scalingOps;
    scalingOps.reserve(pnr);

    for (index_t k=0; k<pnr; ++k)
    {
        const index_t sz = m_localSchurOps[k]->rows();
        gsSparseMatrix<T> scaling(sz,sz);
        for (index_t i=0; i<sz; ++i)
            scaling(i,i) = (T)1/m_localScaling[k](i,0);
        scalingOps.push_back(makeMatrixOp(scaling.moveToPtr()));
    }

    typename gsAdditiveOp<T>::Ptr result = gsAdditiveOp<T>::make();

    for (index_t i=0; i<pnr; ++i)
    {
        typename gsProductOp<T>::Ptr local = gsProductOp<T>::make();
        local->addOperator(scalingOps[i]);
        local->addOperator(m_localSchurOps[i]);
        local->addOperator(scalingOps[i]);
        result->addOperator(m_jumpMatrices[i],local);
    }

    return result;
}


} // namespace gismo
