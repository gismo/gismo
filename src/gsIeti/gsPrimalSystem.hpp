/** @file gsPrimalSystem.hpp

    @brief This class represents the primal system and allows to incorporate the primal constraints

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
    this->m_jumpMatrix.resize(nrLagrangeMultipliers,primalProblemSize);
    this->m_localMatrix.resize(primalProblemSize,primalProblemSize);
    this->m_localRhs.setZero(primalProblemSize,1);
}

template <class T>
void gsPrimalSystem<T>::incorporateConstraints(
        const std::vector<gsSparseVector<T>>& primalConstraints,
        gsSparseMatrix<T,RowMajor>& jumpMatrix,
        gsSparseMatrix<T>& localMatrix,
        gsMatrix<T>& localRhs
    )
{
    const index_t localDofs = localMatrix.rows();
    const index_t nrPrimalConstraints = primalConstraints.size();
    if (nrPrimalConstraints==0) return;

    localMatrix.conservativeResize(localDofs+nrPrimalConstraints, localDofs+nrPrimalConstraints);

    for (index_t i=0; i<nrPrimalConstraints; ++i)
    {
        //localMatrix.block(localDofs+i,0,1,localDofs) = primalConstraints[i];
        //localMatrix.block(0,localDofs+i,localDofs,1) = primalConstraints[i].transpose();
        for (index_t j=0; j<primalConstraints[i].outerSize(); ++j)
            for (typename gsSparseVector<T>::InnerIterator it(primalConstraints[i], j); it; ++it)
            {
                localMatrix(it.row(), localDofs+i) = it.value();
                localMatrix(localDofs+i, it.row()) = it.value();
            }
    }

    localMatrix.makeCompressed();

    GISMO_ASSERT( localRhs.rows() == localDofs,
        "gsPrimalSystem::incorporateConstraint: Right-hand side does not have the expected number of columns;");

    localRhs.conservativeResize(localDofs+nrPrimalConstraints, Eigen::NoChange);
    localRhs.bottomRows(nrPrimalConstraints).setZero();

    GISMO_ASSERT( jumpMatrix.cols() == localDofs,
        "gsPrimalSystem::incorporateConstraint: Jump matrix does not have the expected number of columns;");

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
    const index_t nrPrimalConstraints = primalConstraintsMapper.size();

    GISMO_ASSERT( nrPrimalConstraints<=primalProblemSize, "gsPrimalSystem::primalBasis: "
        "There are more local constrains that there are constraints in total. "
        "Forgot to call gsPrimalSystem::init()?" );

    const index_t localDofs = localSaddlePointSolver->rows() - nrPrimalConstraints;

    gsSparseMatrix<T> result( localDofs, primalProblemSize );

    if (primalProblemSize==0) return result;

    gsMatrix<T> id;
    id.setZero(localDofs+nrPrimalConstraints,nrPrimalConstraints);

    for (index_t i=0; i<nrPrimalConstraints; ++i)
    {
        GISMO_ASSERT( primalConstraintsMapper[i]>=0 && primalConstraintsMapper[i]<primalProblemSize,
            "gsPrimalSystem::primalBasis: Invalid index.");
        id(localDofs+i,i) = 1;
    }

    gsMatrix<T> tmp;
    localSaddlePointSolver->apply(id, tmp);

    for (index_t i=0; i<localDofs; ++i)
        for (index_t j=0; j<nrPrimalConstraints; ++j)
            result(i,primalConstraintsMapper[j]) = tmp(i,j);

    return result;
}

template <class T>
void gsPrimalSystem<T>::addContribution(
        const gsSparseMatrix<T,RowMajor>& jumpMatrix,
        const gsSparseMatrix<T>& localMatrix,
        const gsMatrix<T>& localRhs,
        gsSparseMatrix<T> primalBasis
    )
{
    GISMO_ASSERT( primalBasis.cols() == m_jumpMatrix.cols(),
        "gsPrimalSystem::incorporate: The given problem size does not match the stored primal problem size. "
        "Forgot to call gsPrimalSystem::init()?" );
    const index_t sz = primalBasis.rows();
    m_localMatrix     += primalBasis.transpose() * localMatrix.block(0,0,sz,sz) * primalBasis;
    m_localRhs        += primalBasis.transpose() * localRhs.topRows(sz);
    m_jumpMatrix      += gsSparseMatrix<T,RowMajor>(jumpMatrix.leftCols(sz) * primalBasis);
    m_primalBases.push_back(give(primalBasis));
}

template <class T>
std::vector< gsMatrix<T> > gsPrimalSystem<T>::distributePrimalSolution( std::vector< gsMatrix<T> > sol )
{
    const index_t sz = this->m_primalBases.size();

    // If the primal problem is empty, there might just not be any primal subdomain
    if (sol.size()==sz && this->m_jumpMatrix.cols()==0)
        return sol;

    GISMO_ASSERT(sol.size()==sz+1, "gsPrimalSystem::distributePrimalSolution expects that there "
        "is one more subdomain that patches.");

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
