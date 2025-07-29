/** @file gsIetiMapper.hpp

    @brief Algorithms that help with assembling the matrices required for IETI-solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsAssembler/gsGenericAssembler.h>

/*    Concerning the status flag m_status:
 *       (m_status&1)!=0    means that the object has been initialized by calling init or the value constructor
 *       (m_status&2)!=0    means that the jump matrices have been computed
 *       (m_status&4)!=0    means that corners have been set up as primal constraints
 *       (m_status&flag)!=0 for flag = 8, 16,... means that edges, faces, ... have been set up as primal constraints
 */

namespace gismo
{

template <class T>
void gsIetiMapper<T>::init(
        const gsMultiBasis<T>& multiBasis,
        gsDofMapper dofMapperGlobal,
        const Matrix& fixedPart
    )
{
    GISMO_ASSERT( dofMapperGlobal.componentsSize() == 1, "gsIetiMapper::init: "
        "Got only 1 multi basis, so a gsDofMapper with only 1 component is expected." );
    GISMO_ASSERT( dofMapperGlobal.numPatches() == multiBasis.nBases(), "gsIetiMapper::init: "
        "Number of patches does not agree." );

    const index_t nPatches = dofMapperGlobal.numPatches();
    m_multiBasis = &multiBasis;
    m_dofMapperGlobal = give(dofMapperGlobal);
    m_dofMapperLocal.clear();
    m_dofMapperLocal.resize(nPatches);
    m_fixedPart.clear();
    m_fixedPart.resize(nPatches);
    m_jumpMatrices.clear();
    m_nPrimalDofs = 0;
    m_primalConstraints.clear();
    m_primalConstraints.resize(nPatches);
    m_primalDofIndices.clear();
    m_primalDofIndices.resize(nPatches);
    m_status = 1;

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t nDofs = m_dofMapperGlobal.patchSize(k);
        GISMO_ASSERT( nDofs==m_multiBasis->piece(k).size(), "gsIetiMapper::init: "
            "The mapper for patch "<<k<<" has not as many dofs as the corresponding basis." );

        m_dofMapperLocal[k].setIdentity(1,nDofs);

        // Eliminate boundary dofs (we do not consider the full floating case).
        for (index_t i=0; i<nDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx))
                m_dofMapperLocal[k].eliminateDof(i,0);
        }
        m_dofMapperLocal[k].finalize();

        const index_t szFixedPart = m_dofMapperLocal[k].boundarySize();
        m_fixedPart[k].setZero(szFixedPart,1);
        for (index_t i=0; i<nDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx))
            {
                const index_t globalBoundaryIdx = m_dofMapperGlobal.bindex(i,k);
                const index_t localBoundaryIdx = m_dofMapperLocal[k].bindex(i,0);
                m_fixedPart[k](localBoundaryIdx,0) = fixedPart(globalBoundaryIdx,0);
            }
        }
    }
}


template <class T>
typename gsIetiMapper<T>::Matrix
gsIetiMapper<T>::constructGlobalSolutionFromLocalSolutions( const std::vector<Matrix>& localContribs )
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    GISMO_ASSERT( nPatches == static_cast<index_t>(localContribs.size()),
        "gsIetiMapper::constructGlobalSolutionFromLocalSolutions; The number of local contributions does "
        "not argee with the number of patches." );

    Matrix result;
    result.setZero( m_dofMapperGlobal.freeSize(), localContribs[0].cols() );

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t sz=m_dofMapperLocal[k].size();
        for (index_t i=0; i<sz; ++i)
        {
            if (m_dofMapperLocal[k].is_free(i,0) && m_dofMapperGlobal.is_free(i,k))
                result.row(m_dofMapperGlobal.index(i,k)) = localContribs[k].row(m_dofMapperLocal[k].index(i,0));
        }
    }
    return result;
}

namespace {

struct dof_helper {
    index_t globalIndex;
    index_t patch;
    index_t localIndex;
    bool operator<(const dof_helper& other) const
    { return globalIndex < other.globalIndex; }
};

template<class T>
gsSparseVector<T> makeUnitVector(index_t len, index_t comp)
{
    GISMO_ENSURE(comp<len, "Wrong dimensions.");
    gsSparseVector<T> vec(len);
    vec[comp] = 1;
    //gsInfo << "makeUnitVector("<<len<<", "<<comp<<")\n";
    return vec;
}

template<class SparseMatrixOrVector>
std::vector<index_t> findUnitVectorRow(const SparseMatrixOrVector& v, index_t index, bool transpose = false)
{
    //gsInfo << "fundUnitVectorRow(" << v.toDense() << ",\nindex=" << index << ", transpose=" << transpose << ") -> ";
    gsVector<index_t> data;
    data.setZero(transpose?v.cols():v.rows());
    for (index_t i=0; i<v.outerSize(); ++i)
    {
        for (typename SparseMatrixOrVector::InnerIterator it(v, i); it; ++it)
            if (it.value() != 0)
                data[transpose?it.col():it.row()] |= ((transpose?it.row():it.col()) == index) ? 1 : 2;
    }
    std::vector<index_t> result;
    for (index_t i=0; i<data.rows(); ++i)
        if (data[i]==1)
            result.push_back(i);
    return result;
}

} // namespace

template <class T>
void gsIetiMapper<T>::cornersAsPrimals()
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( !(m_status&4), "gsIetiMapper::cornersAsPrimals: This function has already been called." );
    m_status |= 4;

    const index_t nPatches = m_dofMapperLocal.size();

    // Construct all corners
    std::vector<dof_helper> corners;
    const index_t dim = m_multiBasis->dim();
    corners.reserve((1<<dim)*nPatches);
    // Add corners on all patches
    for (index_t k=0; k<nPatches; ++k)
    {
        for (boxCorner it = boxCorner::getFirst(dim); it!=boxCorner::getEnd(dim); ++it)
        {
            const index_t idx = (*m_multiBasis)[k].functionAtCorner(it);
            dof_helper dh;
            dh.globalIndex = m_dofMapperGlobal.index( idx, k );
            dh.patch = k;
            dh.localIndex = m_dofMapperLocal[k].index( idx, 0 );
            if (m_dofMapperGlobal.is_free_index(dh.globalIndex))
                corners.push_back(dh);
        }
    }

    // Sort corners to collapse corners with same global index
    std::sort(corners.begin(), corners.end());

    // Create data
    index_t lastIndex = -1;
    const index_t sz = corners.size();
    for (index_t i=0; i<sz; ++i)
    {
        if (lastIndex!=corners[i].globalIndex)
        {
            lastIndex = corners[i].globalIndex;
            if (i+1<sz&&corners[i+1].globalIndex==corners[i].globalIndex)
                ++m_nPrimalDofs;
            else
                continue; // Ignore corners that are not shared
        }
        const index_t cornerIndex = m_nPrimalDofs - 1;
        const index_t patch       = corners[i].patch;
        const index_t localIndex  = corners[i].localIndex;

        m_primalConstraints[patch].push_back(makeUnitVector<T>(m_dofMapperLocal[patch].freeSize(),localIndex));
        m_primalDofIndices[patch].push_back(cornerIndex);
    }

}

template <class T>
void gsIetiMapper<T>::declareDofAsPrimal( index_t patch, index_t index, bool checkUnique )
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    if (checkUnique)
    {
        for (size_t i=0; i<m_primalConstraints[patch].size(); ++i)
            if (findUnitVectorRow(m_primalConstraints[patch][i], index, true).size()>0)
                return;
    }

    //if (!m_dofMapperLocal[patch].is_free(index, 0)) return;
    std::pair<index_t, index_t> tmp = m_dofMapperLocal[patch].anyPreImage(index);
    GISMO_ASSERT( tmp.first == 0, "Unexpected patch index, expected 0, got " << tmp.first);
    const index_t indexInLocalBasis = tmp.second;
    const index_t globalIndex = m_dofMapperGlobal.index(indexInLocalBasis, patch);

    std::vector<std::pair<index_t, index_t>> localIndices;
    m_dofMapperGlobal.preImage(globalIndex,localIndices);
    for (size_t i=0; i<localIndices.size(); ++i)
    {
        const index_t k=localIndices[i].first;
        const index_t idx=m_dofMapperLocal[k].index(localIndices[i].second,0);
        m_primalConstraints[k].push_back(makeUnitVector<T>(m_dofMapperLocal[k].freeSize(),idx));
        m_primalDofIndices[k].push_back(m_nPrimalDofs);
    }
    ++m_nPrimalDofs;
}

template <class T>
gsSparseVector<T> gsIetiMapper<T>::assembleAverage(
    const gsGeometry<T>& geo,
    const gsBasis<T>& basis,
    const gsDofMapper& dm,
    boxComponent bc
)
{
    gsMatrix<index_t> indices;

    gsMatrix<T> moments = gsGenericAssembler<T>(
        *(geo.component(bc)),
        *(basis.componentBasis_withIndices(bc, indices, false))
    ).assembleMoments(
        gsConstantFunction<T>(1,geo.targetDim())
    );

    SparseVector constraint( dm.freeSize() );
    T sum = (T)0;
    const index_t sz = moments.size();
    GISMO_ASSERT( sz == indices.size(), "Internal error." );
    for (index_t i=0; i<sz; ++i)
    {
        const index_t idx = dm.index( indices(i,0), 0 );
        if (dm.is_free_index(idx))
        {
            constraint[idx] = moments(i,0);
            sum += moments(i,0);
        }
    }
    return constraint / sum;

}


template <class T>
void gsIetiMapper<T>::interfaceAveragesAsPrimals( const gsMultiPatch<T>& geo, const short_t d )
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( d>0, "gsIetiMapper::interfaceAveragesAsPrimals cannot handle corners." );
    GISMO_ASSERT( d<=m_multiBasis->dim(), "gsIetiMapper::interfaceAveragesAsPrimals: "
        "Interfaces cannot have larger dimension than considered object." );
    GISMO_ASSERT( (index_t)(geo.nPatches()) == m_multiBasis->nPieces(),
        "gsIetiMapper::interfaceAveragesAsPrimals: The given geometry does not fit.");
    GISMO_ASSERT( geo.parDim() == m_multiBasis->dim(),
        "gsIetiMapper::interfaceAveragesAsPrimals: The given geometry does not fit.");

    const unsigned flag = 1<<(2+d);
    GISMO_ASSERT( !(m_status&flag), "gsIetiMapper::interfaceAveragesAsPrimals: This function has "
        " already been called for d="<<d );
    m_status |= flag;

    std::vector< std::vector<patchComponent> > components = geo.allComponents();
    const index_t nComponents = components.size();
    for (index_t n=0; n<nComponents; ++n)
    {
        const index_t sz = components[n].size();
        if ( components[n][0].dim() == d && ( sz > 1 || m_multiBasis->dim() == d ))
        {
            index_t used = 0;
            for (index_t i=0; i<sz; ++i)
            {
                const index_t k = components[n][i].patch();
                gsSparseVector<T> constr = assembleAverage(
                    geo[k],
                    (*m_multiBasis)[k],
                    m_dofMapperLocal[k],
                    components[n][i]
                );
                if ( constr.nonZeros() > 0 )
                {
                    m_primalConstraints[k].push_back(give(constr));
                    m_primalDofIndices[k].push_back(m_nPrimalDofs);
                    ++used;
                }
            }
            GISMO_ASSERT( used==0 || used == sz, "Internal error." );
            if (used)
                ++m_nPrimalDofs;
        }
    }
}


template <class T>
void gsIetiMapper<T>::customPrimalConstraints( std::vector< std::pair<index_t,SparseVector> > data )
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    const index_t sz = data.size();
    for (index_t i=0; i<sz; ++i)
    {
        const index_t patch = data[i].first;
        m_primalConstraints[patch].push_back(give(data[i].second));
        m_primalDofIndices[patch].push_back(m_nPrimalDofs);
    }
    ++m_nPrimalDofs;
}

template <class T>
std::vector<index_t> gsIetiMapper<T>::skeletonDofs( const index_t patch ) const
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    std::vector<index_t> result;
    const index_t patchSize = m_dofMapperGlobal.patchSize(patch);
    const index_t dim = m_multiBasis->dim();
    result.reserve(2*dim*std::pow(patchSize,(1.0-dim)/dim));
    for (index_t i=0; i<patchSize; ++i)
        if ( m_dofMapperGlobal.is_coupled(i,patch) )
            result.push_back(m_dofMapperLocal[patch].index(i,0));
    return result;
}

template <class T>
void gsIetiMapper<T>::computeJumpMatrices( bool fullyRedundant, bool excludeCorners, const std::vector<index_t>& exclude )
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    GISMO_ASSERT( !(m_status&2), "gsIetiMapper::computeJumpMatrices: This function has already been called." );
    m_status |= 2;

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    const index_t coupledSize = m_dofMapperGlobal.coupledSize();

    // Find the groups of to be coupled indices
    std::vector< std::vector< std::pair<index_t,index_t> > > coupling;
    coupling.resize(coupledSize);

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t patchSize = m_dofMapperGlobal.patchSize(k);
        for (index_t i=0; i<patchSize; ++i)
        {
            const index_t globalIndex = m_dofMapperGlobal.index(i,k);
            if ( m_dofMapperGlobal.is_coupled_index(globalIndex) )
            {
                const index_t coupledIndex = m_dofMapperGlobal.cindex(i,k);
                const index_t localIndex = m_dofMapperLocal[k].index(i,0);
                coupling[coupledIndex].push_back(
                    std::pair<index_t,index_t>(k,localIndex)
                );
            }
        }
    }

    // Erease data for corners if so desired
    if (excludeCorners)
    {
        const index_t dim = m_multiBasis->dim();
        for (index_t k=0; k<nPatches; ++k)
        {
            for (boxCorner it = boxCorner::getFirst(dim); it!=boxCorner::getEnd(dim); ++it)
            {
                const index_t idx = (*m_multiBasis)[k].functionAtCorner(it);
                if ( m_dofMapperGlobal.is_coupled(idx, k) )
                {
                    const index_t coupledIndex = m_dofMapperGlobal.cindex(idx,k);
                    coupling[coupledIndex].clear();
                }
            }
        }
    }

    for (std::size_t i=0; i<exclude.size(); ++i)
    {
        const std::pair<index_t, index_t> tmp = m_dofMapperGlobal.anyPreImage(exclude[i]);
        const index_t patch = tmp.first;
        const index_t idx = tmp.second;
        if ( m_dofMapperGlobal.is_coupled(idx, patch) )
        {
            const index_t coupledIndex = m_dofMapperGlobal.cindex(idx, patch);
            coupling[coupledIndex].clear();
        }
    }

    // Compute the number of Lagrange multipliers
    index_t numLagrangeMult = 0;
    for (index_t i=0; i<coupledSize; ++i)
    {
        const index_t n = coupling[i].size();
        if (n==0)
            continue;
        GISMO_ASSERT( n!=1, "Coupled dof that is only coupled to itself.");
        if (fullyRedundant)
            numLagrangeMult += (n * (n-1))/2;
        else
            numLagrangeMult += n-1;
    }

    // Compute the jump matrices
    std::vector< gsSparseEntries<T> > jumpMatrices_se(nPatches);
    for (index_t i=0; i<nPatches; ++i)
    {
        const index_t dim = m_multiBasis->dim();
        jumpMatrices_se[i].reserve(std::pow(m_dofMapperLocal[i].freeSize(),(1.0-dim)/dim));
    }

    index_t multiplier = 0;
    for (index_t i=0; i<coupledSize; ++i)
    {
        const index_t n = coupling[i].size();
        if (n==0)
            continue;
        const index_t maxIndex = fullyRedundant ? (n-1) : 1;
        for (index_t j1=0; j1<maxIndex; ++j1)
        {
            const index_t patch1 = coupling[i][j1].first;
            const index_t localMappedIndex1 = coupling[i][j1].second;
            for (index_t j2=j1+1; j2<n; ++j2)
            {
                const index_t patch2 = coupling[i][j2].first;
                const index_t localMappedIndex2 = coupling[i][j2].second;
                jumpMatrices_se[patch1].add(multiplier,localMappedIndex1,(T)1);
                jumpMatrices_se[patch2].add(multiplier,localMappedIndex2,(T)-1);
                ++multiplier;
            }
        }
    }
    GISMO_ASSERT( multiplier == numLagrangeMult, "gsIetiMapper::computeJumpMatrices: Internal error: "
        << multiplier << "!=" << numLagrangeMult );

    m_jumpMatrices.clear();
    for (index_t i=0; i<nPatches; ++i)
    {
        m_jumpMatrices.push_back(JumpMatrix(numLagrangeMult, m_dofMapperLocal[i].freeSize()));
        m_jumpMatrices[i].setFrom(jumpMatrices_se[i]);
    }

}

template <class T>
gsMatrix<T> gsIetiMapper<T>::incorporateFixedPart(index_t k, const gsMatrix<T>& localSolution) const
{
    GISMO_ASSERT(localSolution.cols() == m_fixedPart[k].cols(), "gsIetiMapper::incorporateFixedPart: Dimension missmatch.");
    const std::size_t sz = m_dofMapperLocal[k].totalSize();
    Matrix coeffs(sz, localSolution.cols());
    for (std::size_t i = 0; i < sz; ++i)
    {
        if ( m_dofMapperLocal[k].is_free(i, 0) )
            coeffs.row(i) = localSolution.row( m_dofMapperLocal[k].index(i, 0) );
        else
            coeffs.row(i) = m_fixedPart[k].row( m_dofMapperLocal[k].bindex(i, 0) );
    }
    return coeffs;
}


} // namespace gismo
