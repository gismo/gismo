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

#include <vector>
#include <iterator>

namespace gismo
{

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByRefinement(
    gsMultiBasis<T> mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    index_t levels,
    index_t numberOfKnotsToBeInserted,
    index_t multiplicityOfKnotsToBeInserted
    )
{
    gsGridHierarchy<T> result;
    result.m_boundaryConditions = boundaryConditions,
    result.m_assemblerOptions = assemblerOptions,
    result.m_mBases.resize(levels);
    result.m_transferMatrices.resize(levels-1);
    result.m_localTransferMatrices.resize(levels-1);
    result.m_mBases[0] = give(mBasis);
    for ( index_t i=1; i<levels; ++i )
        uniformRefine_withTransfer(
            result.m_mBases[i-1],
            result.m_boundaryConditions,
            result.m_assemblerOptions,
            numberOfKnotsToBeInserted,
            multiplicityOfKnotsToBeInserted,
            result.m_mBases[i],
            result.m_transferMatrices[i-1],
            result.m_localTransferMatrices[i-1]
        );
    return result;
}

template <typename T>
index_t sizeOfMultiBasis( const gsMultiBasis<T>& mb )
{
    const index_t nBases = mb.nBases();
    index_t result = 0;
    for (index_t i = 0; i<nBases; ++i)
        result += mb[i].size();
    return result;
}

template <typename T>
gsGridHierarchy<T> gsGridHierarchy<T>::buildByCoarsening(
    gsMultiBasis<T> mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    index_t levels,
    index_t degreesOfFreedom
    )
{
    gsGridHierarchy<T> result;
    result.m_boundaryConditions = boundaryConditions,
    result.m_assemblerOptions = assemblerOptions,

    result.m_mBases.push_back(give(mBasis));

    index_t lastSize = sizeOfMultiBasis(result.m_mBases[0]);

    for (int i = 0; i < levels && lastSize > degreesOfFreedom; ++i)
    {
        gsMultiBasis<T> coarseMBasis;
        gsSparseMatrix<T, RowMajor> transferMatrix;
        std::vector< gsSparseMatrix<T, RowMajor> > localTransferMatrices;
        coarsenMultiBasis_withTransfer(
            result.m_mBases[i],
            boundaryConditions,
            assemblerOptions,
            coarseMBasis,
            transferMatrix,
            localTransferMatrices
        );

        index_t newSize = sizeOfMultiBasis(coarseMBasis);
        // If the number of dofs could not be decreased, then cancel. However, if only the number
        // of levels was specified, then this should be ignored (the caller might need to have a
        // fixed number of levels).
        if (lastSize <= newSize && degreesOfFreedom > 0)
             break;
        lastSize = newSize;

        result.m_mBases.push_back(give(coarseMBasis));
        result.m_transferMatrices.push_back(give(transferMatrix));
        result.m_localTransferMatrices.push_back(give(localTransferMatrices));
    }

    std::reverse( result.m_mBases.begin(), result.m_mBases.end() );
    std::reverse( result.m_transferMatrices.begin(), result.m_transferMatrices.end() );
    std::reverse( result.m_localTransferMatrices.begin(), result.m_localTransferMatrices.end() );

    return result;
}

/*template <typename T>
typename gsMultiGridOp<T>::uPtr gsGridHierarchy<T>::getMultiGrid( gsSparseMatrix<T> systemMatrix, const gsOptionList& opt, const gsMultiPatch<T>* mp )
    { return getMultiGrid( systemMatrix.moveToPtr(), opt, mp ); }

template <typename T>
typename gsMultiGridOp<T>::uPtr gsGridHierarchy<T>::getMultiGrid( typename gsSparseMatrix<T>::Ptr systemMatrix, const gsOptionList& opt, const gsMultiPatch<T>* mp )
{
    GISMO_ASSERT( m_mBases.size() > 0 && m_mBases.size() == m_transferMatrices.size() + 1,
                    "There is a problem with the number of multi bases or the number of transfer matrices.");

    typedef typename gsSparseMatrix<T, RowMajor>::Ptr TransferMatrixPtr;

    const size_t sz = m_transferMatrices.size();

    // TODO: This makes a copy, which is not a big deal, but could be avoided if shared pointers were
    // used internally
    std::vector<TransferMatrixPtr> transferMatrixPtrs(sz);
    for (size_t i = 0; i < sz; ++i)
        transferMatrixPtrs[i] = TransferMatrixPtr(new gsSparseMatrix<T, RowMajor>(m_transferMatrices[i]));

    typename gsMultiGridOp<T>::uPtr mg = gsMultiGridOp<T>::make( give(systemMatrix), give(transferMatrixPtrs) );
    mg->setOptions(opt);
    gsSmootherFactory smf;
    smf.setOptions(opt);
    for (index_t i = 1; i < mg->numLevels(); ++i)
    {
        mg->setSmoother(
            i,
            smf.getSmootherForLevel( i, mg->matrixPtr(i), &(m_mBases[i]), mp, &m_boundaryConditions, &m_assemblerOptions )
        );
    }
    return mg;
}*/

template <typename T>
void uniformRefine_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    index_t refineKnots,
    index_t mult,
    gsMultiBasis<T> & refinedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix,
    std::vector< gsSparseMatrix<T, RowMajor> >& localTransferMatrices
    )
{
    // initialize
    const index_t nBases = mBasis.nBases();
    std::vector< gsBasis<T>* > refinedBases(nBases);
    localTransferMatrices.resize(nBases);

    // setup patchwise without boundary conditions
    for (index_t i=0; i<nBases; ++i )
    {
        refinedBases[i] = (gsBasis<T>*)mBasis[i].clone().release();
        refinedBases[i]->uniformRefine_withTransfer( localTransferMatrices[i], refineKnots, mult );
    }

    // setup the new multibasis object with refined bases and old topology
    refinedMBasis = gsMultiBasis<>(refinedBases, mBasis.topology());

    // get dof dofMappers
    gsDofMapper fineMapper, coarseMapper;

    mBasis.getMapper(
            (dirichlet::strategy)assemblerOptions.askInt("DirichletStrategy",11),
            (iFace    ::strategy)assemblerOptions.askInt("InterfaceStrategy", 1),
            boundaryConditions,
            coarseMapper,
            0
    );

    refinedMBasis.getMapper(
            (dirichlet::strategy)assemblerOptions.askInt("DirichletStrategy",11),
            (iFace    ::strategy)assemblerOptions.askInt("InterfaceStrategy", 1),
            boundaryConditions,
            fineMapper,
            0
    );

    // restrict to free dofs
    combineTransferMatrices( localTransferMatrices, coarseMapper, fineMapper, transferMatrix );

}

template <typename T>
gsKnotVector<T> coarsenKnotVector(const gsKnotVector<T>& kv, std::vector<T>& removedKnots)
{
    std::vector<int> removeIdx;
    std::vector<T> coarseKnots;

    removedKnots.clear();
    removedKnots.reserve( kv.size() / 2 );
    removeIdx.reserve( kv.size() / 2 );

    // determine indices and values of knots to be removed
    const int first = kv.degree() + 1;              // first non-boundary knot
    const int last  = kv.size() - kv.degree() - 2;  // last non-boundary knot
    for (int i = first; i <= last; i += 2)
    {
        removeIdx.push_back( i );
        removedKnots.push_back( kv[i] );
    }

    // copy non-removed knots into coarseKnots
    copyIfNotIndexed( kv.get(), removeIdx, coarseKnots );

    return gsKnotVector<T>( give(coarseKnots), kv.degree() );
}

// Helper function allowing coarsenBasis doing its work
template <unsigned d, typename T>
typename gsTensorBSplineBasis<d, T>::uPtr coarsenTensorBasis(const gsTensorBSplineBasis<d, T>& b, std::vector< std::vector<T> > & removedKnots)
{
    removedKnots.clear();
    removedKnots.resize(d);
    std::vector< gsBasis<T>* > coarseB(d);

    for (unsigned i = 0; i < d; ++i)
        coarseB[i] = new gsBSplineBasis<T>( coarsenKnotVector( b.component(i).knots(), removedKnots[i] ) );

    return gsTensorBSplineBasis<d,T>::make( coarseB );
}

template <unsigned d, typename T>
inline typename gsTensorBSplineBasis<d, T>::uPtr _coarsenTensorBasis(const gsBasis<T>& b, std::vector< std::vector<T> > & removedKnots)
{
    const gsTensorBSplineBasis<d, T> * tb = dynamic_cast<const gsTensorBSplineBasis<d,T>*>(&b);
    if( !tb )
        GISMO_ERROR ("Coarsening only works for gsTensorBSplineBasis.");
    return coarsenTensorBasis<d, T>(*tb,removedKnots);
}

template <typename T>
typename gsBasis<T>::uPtr coarsenBasis(const gsBasis<T>& b, std::vector< std::vector<T> > & removedKnots)
{
    GISMO_ASSERT( (index_t)(b.dim()) == (index_t)(removedKnots.size()), "The dimensions do not agree." );
    switch (b.dim())
    {
        case 1: return _coarsenTensorBasis<1,T>(b, removedKnots);
        case 2: return _coarsenTensorBasis<2,T>(b, removedKnots);
        case 3: return _coarsenTensorBasis<3,T>(b, removedKnots);
        case 4: return _coarsenTensorBasis<4,T>(b, removedKnots);
        default: GISMO_ERROR ("Coarsening is instanciated only for up to 4 dimensions.");
    }
}

// Helper function allowing coarsenBasis_withTransfer doing its work
template <unsigned d, typename T>
typename gsTensorBSplineBasis<d, T>::uPtr coarsenTensorBasis_withTransfer(const gsBasis<T>& b, gsSparseMatrix<T, RowMajor>& transferMatrix)
{
    std::vector< std::vector<T> > removedKnots;
    typename gsTensorBSplineBasis<d, T>::uPtr result = _coarsenTensorBasis<d, T>( b, removedKnots );

    typename gsTensorBSplineBasis<d, T>::uPtr tmp = result->clone();
    tmp->refine_withTransfer( transferMatrix, removedKnots );

    return result;
}

template <typename T>
typename gsBasis<T>::uPtr coarsenBasis_withTransfer(const gsBasis<T>& b, gsSparseMatrix<T, RowMajor>& transferMatrix)
{
    switch (b.dim())
    {
        case 1: return coarsenTensorBasis_withTransfer<1,T>(b, transferMatrix);
        case 2: return coarsenTensorBasis_withTransfer<2,T>(b, transferMatrix);
        case 3: return coarsenTensorBasis_withTransfer<3,T>(b, transferMatrix);
        case 4: return coarsenTensorBasis_withTransfer<4,T>(b, transferMatrix);
        default: GISMO_ERROR ("Coarsening is instanciated only for up to 4 dimensions.");
    }
}

template <typename T>
void coarsenMultiBasis_withTransfer(
    const gsMultiBasis<T>& mBasis,
    const gsBoundaryConditions<T>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    gsMultiBasis<T> & coarsenedMBasis,
    gsSparseMatrix<T, RowMajor>& transferMatrix,
    std::vector< gsSparseMatrix<T, RowMajor> >& localTransferMatrices
    )
{
    // initialize
    const index_t nBases = mBasis.nBases();
    std::vector< gsBasis<T>* > coarsenedBases(nBases);
    localTransferMatrices.resize(nBases);

    // setup patchwise without boundary conditions
    for (index_t i=0; i<nBases; ++i )
        coarsenedBases[i] = coarsenBasis_withTransfer( mBasis[i], localTransferMatrices[i] ).release();

    // setup the new multibasis object with coarsened bases and old topology
    coarsenedMBasis = gsMultiBasis<>(coarsenedBases, mBasis.topology());

    // get dof dofMappers
    gsDofMapper fineMapper, coarseMapper;

    mBasis.getMapper(
            (dirichlet::strategy)assemblerOptions.askInt("DirichletStrategy",11),
            (iFace    ::strategy)assemblerOptions.askInt("InterfaceStrategy", 1),
            boundaryConditions,
            fineMapper,
            0
    );

    coarsenedMBasis.getMapper(
            (dirichlet::strategy)assemblerOptions.askInt("DirichletStrategy",11),
            (iFace    ::strategy)assemblerOptions.askInt("InterfaceStrategy", 1),
            boundaryConditions,
            coarseMapper,
            0
    );

    // restrict to free dofs
    combineTransferMatrices( localTransferMatrices, coarseMapper, fineMapper, transferMatrix );
}

template <typename T>
struct take_first {
    T operator() (const T& a, const T&b) { return a; }
};

template <typename T>
void combineTransferMatrices(
    const std::vector< gsSparseMatrix<T, RowMajor> >& localTransferMatrices,
    const gsDofMapper& coarseMapper,
    const gsDofMapper& fineMapper,
    gsSparseMatrix<T, RowMajor>& transferMatrix
    )
{
    const index_t nBases = localTransferMatrices.size();

    index_t nonzeros = 0;
    index_t trRows = 0;
    index_t trCols = 0;

    for (index_t j=0; j<nBases; ++j)
    {
        nonzeros += localTransferMatrices[j].nonZeros();
        trRows += localTransferMatrices[j].rows();
        trCols += localTransferMatrices[j].cols();
    }

    gsSparseEntries<T> entries;
    entries.reserve( nonzeros );

    for (index_t j=0; j<nBases; ++j)
    {
        for (index_t k=0; k < localTransferMatrices[j].outerSize(); ++k)
        {
            for (typename gsSparseMatrix<T, RowMajor>::iterator it(localTransferMatrices[j],k); it; ++it)
            {
                const index_t coarse_dof_idx = coarseMapper.index(it.col(),j);
                const index_t   fine_dof_idx = fineMapper.index(it.row(),j);

                if (coarseMapper.is_free_index(coarse_dof_idx) && fineMapper.is_free_index(fine_dof_idx))
                    entries.add(fine_dof_idx, coarse_dof_idx, it.value());
            }
        }
    }

    transferMatrix.resize(fineMapper.freeSize(), coarseMapper.freeSize()); 
    transferMatrix.setFromTriplets(entries.begin(), entries.end(), take_first<T>());
    transferMatrix.makeCompressed();
}
    

template <typename Cont>
void copyIfNotIndexed(const Cont& data, const std::vector<index_t>& indices, Cont& result)
{
    result.clear();
    result.reserve( data.size() - indices.size() );

    // Copy blocks between two indices at a time.
    typename Cont::const_iterator itBlockBegin = data.begin();
    for (std::vector<index_t>::const_iterator idxit = indices.begin(); idxit != indices.end(); ++idxit)
    {
        typename Cont::const_iterator itBlockEnd = data.begin() + *idxit;
        GISMO_ASSERT( itBlockBegin <= itBlockEnd, "Indices are not in increasing order." );
        std::copy(itBlockBegin, itBlockEnd, std::back_inserter(result));
        itBlockBegin = itBlockEnd + 1;
    }

    // Copy last block.
    std::copy(itBlockBegin, data.end(), std::back_inserter(result));
}

} // namespace gismo
