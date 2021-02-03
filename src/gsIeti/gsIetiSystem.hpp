/** @file gsIetiSystem.hpp

    @brief This class represents a IETI problem. Its algorithms allow to set up a IETI solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsIeti/gsIetiSystem.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsSumOp.h>
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
void gsIetiSystem<T>::addSubdomain(JumpMatrixPtr jumpMatrix, OpPtr localMatrixOp, Matrix localRhs)
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
}

template<class T>
void gsIetiSystem<T>::setupSparseLUSolvers()
{
    const size_t sz = this->m_localMatrixOps.size();
    this->m_localSolverOps.clear();
    this->m_localSolverOps.reserve(sz);
    for (size_t i=0; i<sz; ++i)
    {
        SparseMatrixOp* matop = dynamic_cast<SparseMatrixOp*>(this->m_localMatrixOps[i].get());
        GISMO_ENSURE( matop, "gsIetiSystem::setupSparseLUSolvers requires the "
          "local systems in localMatrixOps to by of type gsMatrixOp<gsSparseMatrix<T>>." );
        this->m_localSolverOps.push_back(makeSparseLUSolver(SparseMatrix(matop->matrix())));
    }
}

template<class T>
typename gsIetiSystem<T>::OpPtr gsIetiSystem<T>::saddlePointProblem() const
{
    GISMO_ASSERT( this->m_jumpMatrices.size() == this->m_localMatrixOps.size(),
        "gsIeti: The number of jump matrices must match the number of local problems." );
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
    GISMO_ASSERT( this->m_localMatrixOps.empty() || this->m_jumpMatrices.size() == this->m_localMatrixOps.size(),
        "gsIeti: The number of restriction operators must match the number of local problems." );
    GISMO_ASSERT( this->m_localSolverOps.size() == this->m_jumpMatrices.size(),
        "gsIetiSystem::schurComplement() requires solvers for the local subproblems. "
        "Forgot to call setupSparseLUSolvers?" );
    return gsAdditiveOp<T>::make( this->m_jumpMatrices, this->m_localSolverOps );
}


template<class T>
gsMatrix<T> gsIetiSystem<T>::rhsForSchurComplement() const
{
    GISMO_ASSERT( this->m_localMatrixOps.empty() || this->m_jumpMatrices.size() == this->m_localMatrixOps.size(),
        "gsIeti: The number of restriction operators must match the number of local problems." );
    GISMO_ASSERT( this->m_localSolverOps.size() == this->m_jumpMatrices.size(),
        "gsIetiSystem::rhsForSchurComplement() requires solvers for the local subproblems. "
        "Forgot to call setupSparseLUSolvers?" );
    GISMO_ASSERT( this->m_localRhs.size() == this->m_jumpMatrices.size(),
        "gsIetiSystem::rhsForSchurComplement() requires the right-hand sides for the local subproblems." );

    gsMatrix<T> result;
    result.setZero( this->numberOfLagrangeMultipliers(), this->m_localRhs[0].cols());
    const index_t numPatches = this->m_jumpMatrices.size();
    for (index_t i=0; i<numPatches; ++i)
    {
        gsMatrix<T> tmp;
        this->m_localSolverOps[i]->apply( this->m_localRhs[i], tmp );
        result += *(this->m_jumpMatrices[i]) * tmp;
    }
    return result;
}

template<class T>
std::vector< gsMatrix<T> > gsIetiSystem<T>::constructSolutionFromLagrangeMultipliers(const gsMatrix<T>& multipliers) const
{
    GISMO_ASSERT( this->m_localMatrixOps.empty() || this->m_jumpMatrices.size() == this->m_localMatrixOps.size(),
        "gsIeti: The number of restriction operators must match the number of local problems." );
    GISMO_ASSERT( this->m_localSolverOps.size() == this->m_jumpMatrices.size(),
        "gsIetiSystem::constructSolutionFromLagrangeMultipliers() requires solvers for the local subproblems. "
        "Forgot to call setupSparseLUSolvers?" );
    GISMO_ASSERT( this->m_localRhs.size() == this->m_jumpMatrices.size(),
        "gsIetiSystem::constructSolutionFromLagrangeMultipliers() requires the right-hand sides for the local subproblems." );

    const index_t numPatches = this->m_jumpMatrices.size();
    std::vector< gsMatrix<T> > result;
    result.resize(numPatches);
    for (index_t i=0; i<numPatches; ++i)
    {
        this->m_localSolverOps[i]->apply( this->m_localRhs[i]-this->m_jumpMatrices[i]->transpose()*multipliers, result[i] );
    }
    return result;
}

} // namespace gismo
