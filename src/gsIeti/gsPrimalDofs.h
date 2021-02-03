/** @file gsPrimalDofs.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsMatrix/gsVector.h>

namespace gismo
{
#define DEBUGVAR(a) gsInfo << "  " << #a << ": " << a << std::endl
#define DEBUGMATRIX(a) gsInfo << "  " << #a << ": " << a.rows() << " x " << a.cols() << std::endl

template< typename T >
class gsPrimalDofs
{
public:

    void init(index_t primalProblemSize, index_t nrLagrangeMultipliers)
    {
        this->primalProblemSize = primalProblemSize;
        jumpMatrix.resize(nrLagrangeMultipliers,primalProblemSize);
        localMatrix.resize(primalProblemSize,primalProblemSize);
        localRhs.setZero(primalProblemSize,1);
    }

    void incorporate( const std::vector< gsSparseVector<T> >& primalConstraints, const std::vector<index_t>& primalConstraintsMapper,
        gsSparseMatrix<T,RowMajor>& jumpMatrix, gsSparseMatrix<T>& localMatrix, gsMatrix<T>& localRhs)
    {
        const index_t localDofs = localRhs.rows();
        const index_t nrPrimalConstraints = primalProblemSize; //TODO: only active ones

        if (nrPrimalConstraints==0) return;

        gsVector<index_t> activator;
        activator.setZero(nrPrimalConstraints);

        gsSparseMatrix<> localSaddlePointProblem(localDofs+nrPrimalConstraints,localDofs+nrPrimalConstraints);
        for (index_t j=0; j<localMatrix.outerSize(); ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(localMatrix, j); it; ++it)
                localSaddlePointProblem(it.row(), it.col()) = it.value();
        for (index_t k=0; k<primalConstraints.size(); ++k)
            for (index_t j=0; j<primalConstraints[k].outerSize(); ++j)
                for (typename gsSparseVector<T>::InnerIterator it(primalConstraints[k], j); it; ++it)
                {
                    GISMO_ASSERT( it.col() == 0, "Vec??"<<it.col() );
                    const index_t col = primalConstraintsMapper[k];
                    localSaddlePointProblem(it.row(), localDofs+col) = it.value();
                    localSaddlePointProblem(localDofs+col, it.row()) = it.value();
                    activator[col] = 1;
                }
        for (index_t k=0; k<nrPrimalConstraints; ++k)
        {
            if (activator[k]==0)
                localSaddlePointProblem(localDofs+k,localDofs+k) = 1;  //shut off...
        }

        gsMatrix<> rhsForLocalSaddlePointProblem;
        rhsForLocalSaddlePointProblem.setZero(localRhs.rows()+nrPrimalConstraints,1);
        rhsForLocalSaddlePointProblem.topRows(localRhs.rows()) = localRhs;

        gsSparseMatrix<T,RowMajor> extendedJumpMatrix(jumpMatrix.rows(),localDofs+nrPrimalConstraints);
        for (index_t j=0; j<jumpMatrix.outerSize(); ++j)
            for (typename gsSparseMatrix<T,RowMajor>::InnerIterator it(jumpMatrix, j); it; ++it)
                extendedJumpMatrix(it.row(), it.col()) = it.value();


        gsLinearOperator<>::Ptr lu = makeSparseLUSolver(localSaddlePointProblem);

        gsSparseMatrix<> primalBasis(localDofs,nrPrimalConstraints);
        {
            gsMatrix<> id;
            id.setZero(localDofs+nrPrimalConstraints,nrPrimalConstraints);
            for (index_t i=0; i<nrPrimalConstraints; ++i)
                id(localDofs+i,i) = 1;
            gsMatrix<> tmp;
            lu->apply(id, tmp);
            for (index_t i=0; i<localDofs; ++i)
                for (index_t j=0; j<nrPrimalConstraints; ++j)
                    primalBasis(i,j) = tmp(i,j);
        }

        this->localMatrix     += primalBasis.transpose() * localMatrix * primalBasis;
        this->localRhs        += primalBasis.transpose() * localRhs;
        this->jumpMatrix      += gsSparseMatrix<real_t,RowMajor>(jumpMatrix * primalBasis);

        // Write back:

        jumpMatrix   = extendedJumpMatrix;
        localMatrix  = localSaddlePointProblem;
        localRhs     = rhsForLocalSaddlePointProblem;

    }


public:
    index_t primalProblemSize;
    gsSparseMatrix<T, RowMajor> jumpMatrix;
    gsSparseMatrix<T> localMatrix;
    gsMatrix<T> localRhs;



};

} // namespace gismo

//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsPrimalDofs.hpp)
//#endif
