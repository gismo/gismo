/** @file gsIetiSystem.h

    @brief This class represents a IETI problem. Its algorithms allow to set up a IETI solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsAdditiveOp.h>

#define DEBUGVAR(a) gsInfo << "  " << #a << ": " << a << std::endl
#define DEBUGMATRIX(a) gsInfo << "  " << #a << ": " << a.rows() << " x " << a.cols() << std::endl

namespace gismo
{

/** @brief   This class represents a IETI problem. Its algorithms allow to set up a IETI solver.

    \ingroup Solver
*/
template< typename T >
class gsIetiSystem
{
    typedef typename gsLinearOperator<T>::Ptr OpPtr;
    typedef gsSparseMatrix<T,RowMajor> Transfer;
    typedef memory::shared_ptr<Transfer> TransferPtr;
    typedef gsMatrix<T> Matrix;
public:

    void setupSparseLUSolvers()
    {
        const size_t sz = localMatrixOps.size();
        localSolverOps.clear();
        localSolverOps.reserve(sz);
        for (size_t i=0; i<sz; ++i)
        {
            gsMatrixOp< gsSparseMatrix<T> >* matop
              = dynamic_cast< gsMatrixOp< gsSparseMatrix<T> >* >(localMatrixOps[i].get());
            GISMO_ENSURE( matop, "gsIetiSystem::setupSparseLUSolvers requires the "
              "local systems in localMatrixOps to by of type gsMatrixOp< gsSparseMatrix<T> >." );
            localSolverOps.push_back(makeSparseLUSolver(gsSparseMatrix<T>(matop->matrix())));
        }
    }

    OpPtr saddlePointProblem() const
    {
        GISMO_ASSERT( jumpMatrices.size() == localMatrixOps.size(),
            "gsIeti: The number of restriction operators must match the number of local problems." );
        const size_t sz = localMatrixOps.size();
        typename gsBlockOp<T>::Ptr result = gsBlockOp<T>::make( sz+1, sz+1 );
        for (size_t i=0; i<sz; ++i)
        {
            result->addOperator( i, i, localMatrixOps[i] );
            // We hope that the transposed operator does not outlive the non-transposed one.
            result->addOperator( i, sz, makeMatrixOp( jumpMatrices[i]->transpose() ) );
            result->addOperator( sz, i, makeMatrixOp( jumpMatrices[i] ) );
        }
        return result;
    }

    index_t numberOfLagrangeMultipliers() const
    {
        GISMO_ASSERT( jumpMatrices.size()>0, "gsIetiSystem: Number of Lagrange multipliers "
            "can only be determined if there are jump matrices.");
        return jumpMatrices[0]->rows();
    }

    OpPtr schurComplement() const
    {
        GISMO_ASSERT( jumpMatrices.size() == localMatrixOps.size(),
            "gsIeti: The number of restriction operators must match the number of local problems." );
        GISMO_ASSERT( localSolverOps.size() == localMatrixOps.size(),
            "gsIetiSystem::schurComplement() requires solvers for the local subproblems." );
        return gsAdditiveOp<T>::make( jumpMatrices, localSolverOps );
    }

    gsMatrix<T> rhsForSchurComplement() const
    {
        GISMO_ASSERT( localSolverOps.size() == jumpMatrices.size(),
            "gsIetiSystem::rhsForSchurComplement() requires solvers for the local subproblems." );
        GISMO_ASSERT( localRhs.size() == jumpMatrices.size(),
            "gsIetiSystem::rhsForSchurComplement() requires the right-hand sides for the local subproblems." );
        gsMatrix<T> result;
        result.setZero( numberOfLagrangeMultipliers(), localRhs[0].cols());
        const index_t numPatches = jumpMatrices.size();
        for (index_t i=0; i<numPatches; ++i)
        {
            gsMatrix<T> tmp;
            localSolverOps[i]->apply( localRhs[i], tmp );
            result += *(jumpMatrices[i]) * tmp;
        }
        return result;
    }

    std::vector< gsMatrix<T> > constructSolutionFromLagrangeMultipliers(const gsMatrix<T>& multipliers) const
    {
        GISMO_ASSERT( localSolverOps.size() == jumpMatrices.size(),
            "gsIetiSystem::rhsForSchurComplement() requires solvers for the local subproblems." );
        const index_t numPatches = jumpMatrices.size();
        std::vector< gsMatrix<T> > result;
        result.reserve(numPatches);
        for (index_t i=0; i<numPatches; ++i)
        {
            gsMatrix<T> tmp;
            localSolverOps[i]->apply( localRhs[i]-jumpMatrices[i]->transpose()*multipliers, tmp );
            result.push_back(tmp);
        }
        return result;
    }

    void reserve( index_t n )
    {
        localMatrixOps.reserve(n);
        localRhs.reserve(n);
        localSolverOps.reserve(n);
        jumpMatrices.reserve(n);
    }

    struct subdomain {
        TransferPtr jumpMatrix;
        OpPtr       localMatrixOp;
        Matrix      localRhs;
        subdomain(TransferPtr _jumpMatrix, OpPtr _localMatrixOp, Matrix _localRhs)
            : jumpMatrix(_jumpMatrix), localMatrixOp(_localMatrixOp), localRhs(_localRhs) {}
    };

    void addPatch(subdomain data)
    {
        jumpMatrices.push_back(give(data.jumpMatrix));
        localMatrixOps.push_back(give(data.localMatrixOp));
        localRhs.push_back(give(data.localRhs));
    }

public:
    std::vector< TransferPtr >  jumpMatrices;
    std::vector< OpPtr >        localMatrixOps;
    std::vector< Matrix >       localRhs;
    std::vector< OpPtr >        localSolverOps;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIetiSystem.hpp)
#endif
