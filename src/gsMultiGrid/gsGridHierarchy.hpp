/** @file gsGridHierarchy.hpp

    @brief Coarsening algorithms for knot vectors and bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsIO/gsOptionList.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsCore/gsMultiBasis.h>

namespace gismo
{

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByRefinement(
    gsMultiBasis<T> mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& options,
    index_t levels,
    index_t numberOfKnotsToBeInserted,
    index_t multiplicityOfKnotsToBeInserted,
    index_t unk
    )
{
    gsGridHierarchy<T> result;
    result.m_mBases.resize(levels);
    result.m_transferMatrices.resize(levels-1);
    result.m_mBases[0] = give(mBasis);
    for ( index_t i=1; i<levels; ++i )
    {
        result.m_mBases[i] = result.m_mBases[i-1];
        result.m_mBases[i].uniformRefine_withTransfer(
            result.m_transferMatrices[i-1],
            boundaryConditions,
            options,
            numberOfKnotsToBeInserted,
            multiplicityOfKnotsToBeInserted,
            unk
        );
    }
    return result;
}

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByCoarsening(
    gsMultiBasis<T> mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& options,
    index_t levels,
    index_t degreesOfFreedom,
    index_t unk
    )
{
    gsGridHierarchy<T> result;

    result.m_mBases.push_back(give(mBasis));

    index_t lastSize = result.m_mBases[0].totalSize();

    for (index_t i = 0; i < levels-1 && lastSize > degreesOfFreedom; ++i)
    {
        gsSparseMatrix<T, RowMajor> transferMatrix;
        gsMultiBasis<T> coarseMBasis = result.m_mBases[i];
        coarseMBasis.uniformCoarsen_withTransfer(
            transferMatrix,
            boundaryConditions,
            options,
            1,
            unk
        );

        index_t newSize = coarseMBasis.totalSize();
        // If the number of dofs could not be decreased, then cancel. However, if only the number
        // of levels was specified, then this should be ignored (the caller might need to have a
        // fixed number of levels).
        if (lastSize <= newSize && degreesOfFreedom > 0)
             break;
        lastSize = newSize;

        result.m_mBases.push_back(give(coarseMBasis));
        result.m_transferMatrices.push_back(give(transferMatrix));
    }

    std::reverse( result.m_mBases.begin(), result.m_mBases.end() );
    std::reverse( result.m_transferMatrices.begin(), result.m_transferMatrices.end() );

    return result;
}

} // namespace gismo
