/** @file gsPrimalSystem.hpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsIeti/gsPrimalSystem.h>

#include <gsSolver/gsMatrixOp.h>


#define DEBUGVAR(a) gsInfo << "  " << #a << ": " << a << std::endl
#define DEBUGMATRIX(a) gsInfo << "  " << #a << ": " << a.rows() << " x " << a.cols() << std::endl


namespace gismo
{

template <class T>
void gsPrimalSystem<T>::init(index_t primalProblemSize, index_t nrLagrangeMultipliers)
{
    this->primalProblemSize = primalProblemSize;
    jumpMatrix.resize(nrLagrangeMultipliers,primalProblemSize);
    localMatrix.resize(primalProblemSize,primalProblemSize);
    localRhs.setZero(primalProblemSize,1);
}

template <class T>
void gsPrimalSystem<T>::extendLocalSystem(
        const std::vector<gsSparseVector<T>>& primalConstraints,
        const std::vector<index_t>& primalConstraintsMapper,
        gsSparseMatrix<T,RowMajor>& jumpMatrix,
        gsSparseMatrix<T>& localMatrix,
        gsMatrix<T>& localRhs,
        index_t primalProblemSize
    )
{
    const index_t localDofs = localMatrix.rows();
    const index_t nrPrimalConstraints = primalProblemSize; //TODO: only active ones
    if (nrPrimalConstraints==0) return;

    gsVector<index_t> activator;
    activator.setZero(nrPrimalConstraints);

    localMatrix.conservativeResize(localDofs+nrPrimalConstraints, localDofs+nrPrimalConstraints);

    for (index_t k=0; k<primalConstraints.size(); ++k)
    {
        const index_t col = primalConstraintsMapper[k];
        //localMatrix.block(localDofs+col,0,1,localDofs) = primalConstraints[k];
        //localMatrix.block(0,localDofs+col,localDofs,1) = primalConstraints[k].transpose();
        for (index_t j=0; j<primalConstraints[k].outerSize(); ++j)
            for (typename gsSparseVector<T>::InnerIterator it(primalConstraints[k], j); it; ++it)
            {
                localMatrix(it.row(), localDofs+col) = it.value();
                localMatrix(localDofs+col, it.row()) = it.value();
                activator[col] = 1;
            }
    }
    for (index_t k=0; k<nrPrimalConstraints; ++k)
    {
        if (activator[k]==0)
            localMatrix(localDofs+k,localDofs+k) = 1;  //shut off...
        // TODO: instead of this, we should just make the problem smaller
    }

    localMatrix.makeCompressed();

    GISMO_ASSERT( localRhs.rows() == localDofs,
        "gsPrimalSystem::extendLocalSystem: Right-hand side does not have the expected number of columns;");

    localRhs.conservativeResize(localDofs+nrPrimalConstraints, Eigen::NoChange);
    localRhs.bottomRows(nrPrimalConstraints).setZero();

    GISMO_ASSERT( jumpMatrix.cols() == localDofs,
        "gsPrimalSystem::extendLocalSystem: Jump matrix does not have the expected number of columns;");

    jumpMatrix.conservativeResize(jumpMatrix.rows(), localDofs+nrPrimalConstraints);
    jumpMatrix.makeCompressed();

}

template <class T>
gsSparseMatrix<T> gsPrimalSystem<T>::primalBasis(
        typename gsLinearOperator<T>::Ptr localSaddlePointSolver,
        const std::vector<index_t>& primalConstraintsMapper,
        index_t primalProblemSize
    )
{
    //const index_t nrPrimalConstraints = primalConstraintsMapper.size();
    const index_t nrPrimalConstraints = primalProblemSize;
    const index_t localDofs = localSaddlePointSolver->rows() - primalProblemSize;

    gsSparseMatrix<T> result( localDofs, nrPrimalConstraints );
    if (nrPrimalConstraints==0) return result;

    gsMatrix<T> id;
    id.setZero(localDofs+primalProblemSize,nrPrimalConstraints);

    for (index_t i=0; i<nrPrimalConstraints; ++i)
        id(localDofs+i,i) = 1;

    gsMatrix<T> tmp;
    localSaddlePointSolver->apply(id, tmp);

    for (index_t i=0; i<localDofs; ++i)
        for (index_t j=0; j<nrPrimalConstraints; ++j)
            result(i,j) = tmp(i,j);

    return result;
}

template <class T>
void gsPrimalSystem<T>::incorporate(
        const gsSparseMatrix<T,RowMajor>& jumpMatrix,
        const gsSparseMatrix<T>& localMatrix,
        const gsMatrix<T>& localRhs,
        const gsSparseMatrix<T>& primalBasis
    )
{
    GISMO_ASSERT( primalBasis.cols() == primalProblemSize,
        "gsPrimalSystem::incorporate: The given problem size does not match the stored primal problem size" );
    const index_t sz = primalBasis.rows();
    this->localMatrix     += primalBasis.transpose() * localMatrix.block(0,0,sz,sz) * primalBasis;
    this->localRhs        += primalBasis.transpose() * localRhs.topRows(sz);
    gsSparseMatrix<T,RowMajor> tmp(jumpMatrix.leftCols(sz) * primalBasis);
    this->jumpMatrix      += tmp;
}

template <class T>
void gsPrimalSystem<T>::handleConstraints(
        const std::vector< gsSparseVector<T> >& primalConstraints,
        const std::vector<index_t>& primalConstraintsMapper,
        gsSparseMatrix<T,RowMajor>& jumpMatrix,
        gsSparseMatrix<T>& localMatrix,
        gsMatrix<T>& localRhs
    )
{
    extendLocalSystem(
        primalConstraints,
        primalConstraintsMapper,
        jumpMatrix,
        localMatrix,
        localRhs,
        primalProblemSize
    );

    incorporate(
        jumpMatrix,
        localMatrix,
        localRhs,
        primalBasis( makeSparseLUSolver( localMatrix ), primalConstraintsMapper, primalProblemSize )
    );
}

} // namespace gismo
