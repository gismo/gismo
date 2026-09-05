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
#include <gsHSplines/gsHTensorBasis.h>
#include <gsCore/gsOptionList.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsContainers/gsMultiBasis.h>
#include <gsAssembler/gsDofMapperCreator.h>

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
gsGridHierarchy<T> gsGridHierarchy<T>::buildByRefinement(
    gsMultiBasis<T> mBasis,
    index_t numComponents,
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

    bool glueInterfaces = iFace::glue==options.askInt("InterfaceStrategy", 1);
    bool eliminateDirichlet = dirichlet::elimination==options.askInt("DirichletStrategy",11);
    gsDofMapper coarseMapper, fineMapper;
    std::vector< gsSparseMatrix<T, RowMajor> > localTransferMatrices(result.m_mBases[0].nBases());
    for ( index_t i=1; i<levels; ++i )
    {
        result.m_mBases[i] = result.m_mBases[i-1];
        gsMultiBasis<T> & basis = result.m_mBases[i];
        // Store the mapper for the current level
        coarseMapper = eliminateDirichlet
            ? createMapper(basis, boundaryConditions, numComponents, unk, glueInterfaces)
            : createMapper(basis, numComponents, glueInterfaces);
        coarseMapper.finalize();
        // Refine each base and store the local transfer matrices
        for (size_t j = 0; j < basis.nBases(); ++j)
            basis.basis(j).uniformRefine_withTransfer(localTransferMatrices[j],numberOfKnotsToBeInserted,multiplicityOfKnotsToBeInserted);
        // Store the fine mapper
        fineMapper = eliminateDirichlet
            ? createMapper(basis, boundaryConditions, numComponents, unk, glueInterfaces)
            : createMapper(basis, numComponents, glueInterfaces);
        fineMapper.finalize();
        // Compute the global transfer matrix
        basis.combineTransferMatrices(localTransferMatrices,coarseMapper,fineMapper,result.m_transferMatrices[i-1]);
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

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByCoarsening(
    gsMultiBasis<T> mBasis,
    index_t numComponents,
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

    bool glueInterfaces = iFace::glue==options.askInt("InterfaceStrategy", 1);
    bool eliminateDirichlet = dirichlet::elimination==options.askInt("DirichletStrategy",11);
    std::vector< gsSparseMatrix<T, RowMajor> > localTransferMatrices(result.m_mBases[0].nBases());
    gsDofMapper coarseMapper, fineMapper;
    for (index_t i = 0; i < levels-1 && lastSize > degreesOfFreedom; ++i)
    {
        gsSparseMatrix<T, RowMajor> transferMatrix;
        gsMultiBasis<T> coarseMBasis = result.m_mBases[i];
        // Store the mapper for the current level
        fineMapper = eliminateDirichlet
            ? createMapper(coarseMBasis, boundaryConditions, numComponents, unk, glueInterfaces)
            : createMapper(coarseMBasis, numComponents, glueInterfaces);
        fineMapper.finalize();
        // Coarsen each base and store the local transfer matrices
        for (size_t j = 0; j < coarseMBasis.nBases(); ++j)
            coarseMBasis.basis(j).uniformCoarsen_withTransfer(localTransferMatrices[j],1);
        // Store the coarse mapper
        coarseMapper = eliminateDirichlet
            ? createMapper(coarseMBasis, boundaryConditions, numComponents, unk, glueInterfaces)
            : createMapper(coarseMBasis, numComponents, glueInterfaces);
        coarseMapper.finalize();
        // Compute the global transfer matrix
        coarseMBasis.combineTransferMatrices(localTransferMatrices,coarseMapper,fineMapper,transferMatrix);

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

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByHierarchicalCoarsening(
    gsMultiBasis<T> mBasis,
    index_t numComponents,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& options,
    index_t unk
    )
{
    // Check all bases
    index_t maxLevel = 0;
    for (size_t i = 0; i < mBasis.nBases(); ++i)
        if (gsHTensorBasis<1,T>* basis = dynamic_cast<gsHTensorBasis<1,T>*>(&mBasis.basis(i)))
            maxLevel = math::max(maxLevel, (index_t)basis->maxLevel());
        else if (gsHTensorBasis<2,T>* basis = dynamic_cast<gsHTensorBasis<2,T>*>(&mBasis.basis(i)))
            maxLevel = math::max(maxLevel, (index_t)basis->maxLevel());
        else if (gsHTensorBasis<3,T>* basis = dynamic_cast<gsHTensorBasis<3,T>*>(&mBasis.basis(i)))
            maxLevel = math::max(maxLevel, (index_t)basis->maxLevel());
        else if (gsHTensorBasis<4,T>* basis = dynamic_cast<gsHTensorBasis<4,T>*>(&mBasis.basis(i)))
            maxLevel = math::max(maxLevel, (index_t)basis->maxLevel());
        else
            GISMO_ERROR("Basis " << i << " must be hierarchical.");
    maxLevel = math::min(maxLevel, (index_t)options.askInt("Levels", 2));


    gsGridHierarchy<T> result;

    result.m_mBases.push_back(give(mBasis));

    index_t lastSize = result.m_mBases[0].totalSize();

    bool glueInterfaces = iFace::glue==options.askInt("InterfaceStrategy", 1);
    bool eliminateDirichlet = dirichlet::elimination==options.askInt("DirichletStrategy",11);
    std::vector< gsSparseMatrix<T, RowMajor> > localTransferMatrices(result.m_mBases[0].nBases());
    gsDofMapper coarseMapper, fineMapper;
    for (index_t i = 0; i != maxLevel; i++)
    {
        gsSparseMatrix<T, RowMajor> transferMatrix;
        gsMultiBasis<T> coarseMBasis = result.m_mBases[i];
        // Store the mapper for the current level
        fineMapper = eliminateDirichlet
            ? createMapper(coarseMBasis, boundaryConditions, numComponents, unk, glueInterfaces)
            : createMapper(coarseMBasis, numComponents, glueInterfaces);
        fineMapper.finalize();
        gsDebugVar(fineMapper.freeSize());
        // Coarsen each base and store the local transfer matrices
        for (size_t j = 0; j < coarseMBasis.nBases(); ++j)
            if (coarseMBasis.basis(j).domainDim()==1)
                static_cast<gsHTensorBasis<1,T>*>(&coarseMBasis.basis(j))->unrefineFinestLevel_withTransfer(localTransferMatrices[j]);
            else if (coarseMBasis.basis(j).domainDim()==2)
                static_cast<gsHTensorBasis<2,T>*>(&coarseMBasis.basis(j))->unrefineFinestLevel_withTransfer(localTransferMatrices[j]);
            else if (coarseMBasis.basis(j).domainDim()==3)
                static_cast<gsHTensorBasis<3,T>*>(&coarseMBasis.basis(j))->unrefineFinestLevel_withTransfer(localTransferMatrices[j]);
            else if (coarseMBasis.basis(j).domainDim()==4)
                static_cast<gsHTensorBasis<4,T>*>(&coarseMBasis.basis(j))->unrefineFinestLevel_withTransfer(localTransferMatrices[j]);
            else
                // This should never happen due to the check above
                GISMO_ERROR("Basis " << j << " must be hierarchical.");
        // Store the coarse mapper
        coarseMapper = eliminateDirichlet
            ? createMapper(coarseMBasis, boundaryConditions, numComponents, unk, glueInterfaces)
            : createMapper(coarseMBasis, numComponents, glueInterfaces);
        coarseMapper.finalize();
        gsDebugVar(coarseMapper.freeSize());
        // Compute the global transfer matrix
        coarseMBasis.combineTransferMatrices(localTransferMatrices,coarseMapper,fineMapper,transferMatrix);

        index_t newSize = coarseMBasis.totalSize();
        // If the number of dofs could not be decreased, then cancel. However, if only the number
        // of levels was specified, then this should be ignored (the caller might need to have a
        // fixed number of levels).
        if (lastSize <= newSize)
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
