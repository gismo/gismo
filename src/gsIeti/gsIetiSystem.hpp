/** @file gsIetiSystem.hpp

    @brief This class represents a IETI system

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsAdditiveOp.h>

namespace gismo
{

template<class T>
void gsIetiSystem<T>::reserve(index_t n)
{
    this->m_localMatrixOps.reserve(n);
    this->m_localRhs.reserve(n);
    this->m_localSolverOps.reserve(n);
    this->m_jumpMatrices.reserve(n);
}

template<class T>
void gsIetiSystem<T>::addSubdomain(JumpMatrixPtr jumpMatrix, OpPtr localMatrixOp, Matrix localRhs, OpPtr localSolverOp)
{
    GISMO_ASSERT( jumpMatrix->cols() == localMatrixOp->cols()
        && localMatrixOp->rows() == localMatrixOp->cols()
        && localRhs.rows() == localMatrixOp->cols(),
        "gsIetiSystem::addPatch: Dimensions do not agree." );

    GISMO_ASSERT( this->m_jumpMatrices.empty() || this->m_jumpMatrices[0]->rows() == jumpMatrix->rows(),
        "gsIetiSystem::addPatch: The dimension of the given transfer matrix "
        "is incompatible to the transfer matrices given previously." );

    this->m_jumpMatrices.push_back(give(jumpMatrix));
    this->m_localMatrixOps.push_back(give(localMatrixOp));
    this->m_localRhs.push_back(give(localRhs));
    this->m_localSolverOps.push_back(give(localSolverOp));
}

template<class T>
void gsIetiSystem<T>::setupSparseLUSolvers() const
{
    const size_t sz = this->m_localSolverOps.size();
    for (size_t i=0; i<sz; ++i)
    {
        if (!m_localSolverOps[i]) // If not yet provided...
        {
            SparseMatrixOp* matop = dynamic_cast<SparseMatrixOp*>(this->m_localMatrixOps[i].get());
            GISMO_ENSURE( matop, "gsIetiSystem::setupSparseLUSolvers The local solvers can only "
              "be computed on the fly if the local systems in localMatrixOps are of type "
              "gsMatrixOp<gsSparseMatrix<T>>. Please provide solvers via members .addSubdomain "
              "or .solverOp" );
            this->m_localSolverOps[i] = makeSparseLUSolver(SparseMatrix(matop->matrix()));
        }
    }
}

template<class T>
typename gsIetiSystem<T>::OpPtr gsIetiSystem<T>::saddlePointProblem() const
{
    const size_t sz = this->m_localMatrixOps.size();
    typename gsBlockOp<T>::Ptr result = gsBlockOp<T>::make( sz+1, sz+1 );
    for (size_t i=0; i<sz; ++i)
    {
        result->addOperator( i, i, m_localMatrixOps[i] );
        // We hope that the transposed operator does not outlive the non-transposed one.
        result->addOperator( i, sz, makeMatrixOp( this->m_jumpMatrices[i]->transpose() ) );
        result->addOperator( sz, i, makeMatrixOp( this->m_jumpMatrices[i] ) );
    }
    return result;
}

template<class T>
typename gsIetiSystem<T>::OpPtr gsIetiSystem<T>::schurComplement() const
{
    setupSparseLUSolvers();
    return gsAdditiveOp<T>::make( this->m_jumpMatrices, this->m_localSolverOps );
}


template<class T>
gsMatrix<T> gsIetiSystem<T>::rhsForSchurComplement() const
{
    setupSparseLUSolvers();
    Matrix result;
    result.setZero( this->nLagrangeMultipliers(), this->m_localRhs[0].cols());
    const index_t numPatches = this->m_jumpMatrices.size();
    for (index_t i=0; i<numPatches; ++i)
    {
        Matrix tmp;
        this->m_localSolverOps[i]->apply( this->m_localRhs[i], tmp );
        result += *(this->m_jumpMatrices[i]) * tmp;
    }
    return result;
}

template<class T>
std::vector< gsMatrix<T> > gsIetiSystem<T>::constructSolutionFromLagrangeMultipliers(const Matrix& multipliers) const
{
    setupSparseLUSolvers();

    const index_t numPatches = this->m_jumpMatrices.size();
    std::vector<Matrix> result;
    result.resize(numPatches);
    for (index_t i=0; i<numPatches; ++i)
    {
        this->m_localSolverOps[i]->apply( this->m_localRhs[i]-this->m_jumpMatrices[i]->transpose()*multipliers, result[i] );
    }
    return result;
}

template<class T>
gsMatrix<T> gsIetiSystem<T>::rhsForSaddlePoint() const
{
    const index_t sz = m_localMatrixOps.size();
    index_t rows = nLagrangeMultipliers();
    for (index_t k=0; k<sz; ++k)
        rows += m_localRhs[k].rows();
    Matrix result;
    result.setZero(rows,m_localRhs[0].cols());
    index_t i=0;
    for (index_t k=0; k<sz; ++k)
    {
        const index_t l = m_localRhs[k].rows();
        result.middleRows(i, l) = m_localRhs[k];
        i += l;
    }
    return result;
}

template<class T>
std::vector< gsMatrix<T> > gsIetiSystem<T>::constructSolutionFromSaddlePoint(const Matrix& x) const
{
    const index_t sz = m_localMatrixOps.size();
    std::vector<Matrix> result;
    result.reserve(sz);
    index_t i=0;
    for (index_t k=0; k<sz; ++k)
    {
        const index_t l = m_localRhs[k].rows();
        result.push_back( x.middleRows(i, l) );
        i += l;
    }
    return result;
}


} // namespace gismo
