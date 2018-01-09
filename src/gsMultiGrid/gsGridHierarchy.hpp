/** @file gsGridHierarchy.hpp

    @brief Coarsening algorithms for knot vectors and bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsMultiGrid/gsGridHierarchy.h>

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
    index_t multiplicityOfKnotsToBeInserted
    )
{
    gsGridHierarchy<T> result;
    result.m_boundaryConditions = boundaryConditions,
    result.m_options = options,
    result.m_mBases.resize(levels);
    result.m_transferMatrices.resize(levels-1);
    result.m_mBases[0] = give(mBasis);
    for ( index_t i=1; i<levels; ++i )
        uniformRefine_withTransfer(
            result.m_mBases[i-1],
            result.m_boundaryConditions,
            result.m_options,
            numberOfKnotsToBeInserted,
            multiplicityOfKnotsToBeInserted,
            result.m_mBases[i],
            result.m_transferMatrices[i-1]
        );
    return result;
}

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByCoarsening(
    gsMultiBasis<T> mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& options,
    index_t levels,
    index_t degreesOfFreedom
    )
{
    gsGridHierarchy<T> result;
    result.m_boundaryConditions = boundaryConditions,
    result.m_options = options,

    result.m_mBases.push_back(give(mBasis));

    index_t lastSize = result.m_mBases[0].totalSize();

    for (int i = 0; i < levels && lastSize > degreesOfFreedom; ++i)
    {
        gsMultiBasis<T> coarseMBasis;
        gsSparseMatrix<T, RowMajor> transferMatrix;
        coarsenMultiBasis_withTransfer(
            result.m_mBases[i],
            boundaryConditions,
            options,
            coarseMBasis,
            transferMatrix
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
