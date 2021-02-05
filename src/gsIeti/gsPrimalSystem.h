/** @file gsPrimalSystem.h

    @brief This class represents the primal system and allows to incorporate the primal constraints

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsMatrix/gsVector.h>

namespace gismo
{

template< typename T >
class gsPrimalSystem
{
public:

    void init(index_t primalProblemSize, index_t nrLagrangeMultipliers);

    static void extendLocalSystem(
        const std::vector<gsSparseVector<T>>& primalConstraints,
        gsSparseMatrix<T,RowMajor>& jumpMatrix,
        gsSparseMatrix<T>& localMatrix,
        gsMatrix<T>& localRhs
    );

    static gsSparseMatrix<T> primalBasis(
        typename gsLinearOperator<T>::Ptr localSaddlePointSolver,
        const std::vector<index_t>& primalConstraintsMapper,
        index_t primalProblemSize
    );

    void incorporate(
        const gsSparseMatrix<T,RowMajor>& jumpMatrix,
        const gsSparseMatrix<T>& localMatrix,
        const gsMatrix<T>& localRhs,
        gsSparseMatrix<T> primalBasis
    );

    void handleConstraints(
        const std::vector< gsSparseVector<T> >& primalConstraints,
        const std::vector<index_t>& primalConstraintsMapper,
        gsSparseMatrix<T,RowMajor>& jumpMatrix,
        gsSparseMatrix<T>& localMatrix,
        gsMatrix<T>& localRhs
    );

    std::vector< gsMatrix<T> > distributePrimalSolution( std::vector< gsMatrix<T> > sol );

public:
    index_t primalProblemSize;
    gsSparseMatrix<T, RowMajor> jumpMatrix;
    gsSparseMatrix<T> localMatrix;
    gsMatrix<T> localRhs;
    std::vector< gsSparseMatrix<T> > primalBases;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPrimalSystem.hpp)
#endif
