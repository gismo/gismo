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
    GISMO_ASSERT( dofMapperGlobal.componentsSize() == 1, "gsIetiMapper<T>::init: "
        "Got only 1 multi basis, so a gsDofMapper with only 1 component is expected." );
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
        const index_t patchDofs = m_dofMapperGlobal.patchSize(k);
        m_dofMapperLocal[k].setIdentity(1,patchDofs);

        for (index_t i=0; i<patchDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx))
                m_dofMapperLocal[k].eliminateDof(i,0);
        }
        m_dofMapperLocal[k].finalize();

        const index_t szFixedPart = m_dofMapperLocal[k].boundarySize();
        m_fixedPart[k].setZero(szFixedPart,1);
        for (index_t i=0; i<patchDofs; ++i)
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
void gsIetiMapper<T>::setFixedPart(const Matrix &fixedPart)
{
    index_t nPatches = m_dofMapperGlobal.numPatches();
    for (index_t k = 0; k < nPatches; ++k) {
        const index_t szFixedPart = m_dofMapperLocal[k].boundarySize();
        m_fixedPart[k].setZero(szFixedPart,1);
        for (size_t i=0; i<m_dofMapperLocal[k].patchSize(0); ++i)
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

namespace{
struct dof_helper {
    index_t globalIndex;
    index_t patch;
    index_t localIndex;
    bool operator<(const dof_helper& other) const
    { return globalIndex < other.globalIndex; }
};
}

template <class T>
void gsIetiMapper<T>::cornersAsPrimals()
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( !(m_status&4), "gsIetiMapper::cornersAsPrimals: This function has already been called." );
    m_status |= 4;

    const index_t nPatches = m_dofMapperLocal.size();

    // Construct all corners
    const short_t dim = m_multiBasis->dim();
    std::vector<dof_helper> corners;
    corners.reserve(math::pow(2,dim)*nPatches);

    for (index_t k=0; k<nPatches; ++k){
        for (boxCorner it = boxCorner::getFirst(dim); it!=boxCorner::getEnd(dim); ++it)
        {
            const index_t idx = (*m_multiBasis)[k].functionAtCorner(it);
            dof_helper dh;
            dh.globalIndex = m_dofMapperGlobal.index( idx, k );
            dh.patch = k;
            dh.localIndex = m_dofMapperLocal[k].index( idx, 0 );
            if (m_dofMapperGlobal.is_free_index(dh.globalIndex))
            {
                corners.push_back(dh);
                for (size_t n=0; n<m_artificialEdges[k].size(); ++n)
                {
                    const patchSide s = m_artificialEdges[k][n].second;
                    const patchSide t = m_artificialEdges[k][n].first;
                    index_t idxAsArtificial = dgFindCorrespondingExtraIndex(s, t, idx);
                    dof_helper artificialDh;
                    artificialDh.globalIndex = m_dofMapperGlobal.index(idxAsArtificial, s.patch);
                    artificialDh.patch = s.patch;
                    artificialDh.localIndex = (m_dofMapperLocal[s.patch].index( idxAsArtificial,0));
                    corners.push_back(artificialDh);
                }
            }
        }
    }
    std::sort(corners.begin(), corners.end());

    // Create data
    index_t lastIndex=-1;
    const index_t sz = corners.size();
    //First add the "real" primal constraints ....
    for (index_t i=0; i<sz; ++i) {
        if (lastIndex != corners[i].globalIndex) {
            lastIndex = corners[i].globalIndex;
            ++m_nPrimalDofs;
        }
        const index_t cornerIndex = m_nPrimalDofs - 1;
        const index_t patch       = corners[i].patch;
        const index_t localIndex  = corners[i].localIndex;

        SparseVector constr(m_dofMapperLocal[patch].freeSize());
        constr[localIndex] = 1;

        m_primalConstraints[patch].push_back(give(constr));
        m_primalDofIndices[patch].push_back(cornerIndex);
    }

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
void gsIetiMapper<T>::interfaceAveragesAsPrimals(const gsMultiPatch<T>& geo, const short_t d)
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( d>0, "gsIetiMapper::interfaceAveragesAsPrimals cannot handle corners." );
    GISMO_ASSERT( d<m_multiBasis->dim(), "gsIetiMapper::interfaceAveragesAsPrimals: "
        "Interfaces must have smaller dimension than considered object." );
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
        if ( sz > 1 && components[n][0].dim() == d )
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
void gsIetiMapper<T>::customPrimalConstraints(std::vector< std::pair<index_t,SparseVector> > data)
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
    result.reserve(4*math::sqrt(patchSize));
    for (index_t i=0; i<patchSize; ++i)
    {
        if ( m_dofMapperGlobal.is_coupled(i,patch) )
            result.push_back(m_dofMapperLocal[patch].index(i,0));
    }
    return result;
}

template <class T>
void gsIetiMapper<T>::computeJumpMatrices(bool fullyRedundant, bool excludeCorners)
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    GISMO_ASSERT( !(m_status&2), "gsIetiMapper::computeJumpMatrices: This function has already been called." );
    m_status |= 2;

    if(m_status&32)
        unfinalizeMapper();


    const index_t nPatches = m_dofMapperGlobal.numPatches();
    const index_t coupledSize = m_dofMapperGlobal.coupledSize();

    // Find the groups of to be coupled indices
    std::vector< std::vector< std::pair<index_t,index_t> > > coupling;
    coupling.resize(coupledSize);

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t patchSize = m_multiBasis->basis(k).size();
        for (index_t i=0; i<patchSize; ++i)
        {
            const index_t globalIndex = m_dofMapperGlobal.index(i,k);
            if ( m_dofMapperGlobal.is_coupled_index(globalIndex) )
            {
                const index_t coupledIndex = m_dofMapperGlobal.cindex(i,k);
                coupling[coupledIndex].push_back(
                    std::pair<index_t,index_t>(k,i)
                );
            }
        }
    }

    // Erease data for corners if so desired
    if (excludeCorners)
    {
        const index_t dim = m_multiBasis->dim();
        for (index_t k=0; k<nPatches; ++k)
            for (boxCorner it = boxCorner::getFirst(dim); it!=boxCorner::getEnd(dim); ++it)
            {
                const index_t idx = (*m_multiBasis)[k].functionAtCorner(it);
                const index_t globalIndex = m_dofMapperGlobal.index(idx,k);
                if ( m_dofMapperGlobal.is_coupled_index(globalIndex) )
                {
                    const index_t coupledIndex = m_dofMapperGlobal.cindex(idx,k);
                    coupling[coupledIndex].clear();
                }
            }
    }

    // Compute the number of Lagrange multipliers
    index_t numLagrangeMult = 0;
    for (index_t i=0; i<coupledSize; ++i)
    {
        const index_t n = coupling[i].size();
        GISMO_ASSERT( n>1 || excludeCorners, "gsIetiMapper::computeJumpMatrices:"
            "Found a coupled dof that is not coupled to any other dof." );
        if (fullyRedundant)
            numLagrangeMult += (n * (n-1))/2;
        else
            numLagrangeMult += n-1;
    }

    // Compute the jump matrices
    std::vector< gsSparseEntries<T> > jumpMatrices_se(nPatches);
    for (index_t i=0; i<nPatches; ++i)
        jumpMatrices_se[i].reserve( math::sqrt( m_dofMapperLocal[i].freeSize() ) );

    index_t multiplier = 0;
    for (index_t i=0; i<coupledSize; ++i)
    {
        const index_t n = coupling[i].size();
        index_t maxIndex = fullyRedundant ? (n-1) : 1;
        for (index_t j1=0; j1<maxIndex; ++j1)
        {
            const index_t patch1 = coupling[i][j1].first;
            const index_t localIndex1 = coupling[i][j1].second;
            const index_t localMappedIndex1 = m_dofMapperLocal[patch1].index(localIndex1,0);
            for (index_t j2=j1+1; j2<n; ++j2)
            {
                const index_t patch2 = coupling[i][j2].first;
                const index_t localIndex2 = coupling[i][j2].second;
                const index_t localMappedIndex2 = m_dofMapperLocal[patch2].index(localIndex2,0);
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
        //gsInfo << "jumpmat\n" << m_jumpMatrices[i].toDense() << "\n";
    }

}

template<class T>
void gsIetiMapper<T>::fullDG()
{
    index_t nPatches = m_dofMapperGlobal.numPatches();
    m_artificialEdges.clear();
    m_artificialEdges.resize(nPatches);
    m_artificialDofsPerSide.clear();
    m_artificialDofsPerSide.resize(nPatches);

    gsVector<index_t> artificialDofs(nPatches); artificialDofs.setZero();
    gsDofMapper dGMapper;
    std::vector<gsDofMapper> localdGmappers;

    gsVector<index_t> localSpaceDofs(nPatches); localSpaceDofs.setZero();
    for (index_t k = 0; k < nPatches; ++k) {
        m_artificialDofsPerSide[k] = {0,0,0,0,0};
        localSpaceDofs[k] += m_dofMapperGlobal.patchSize(k);
    }
    const gsBoxTopology& topology = m_multiBasis->topology();
    for(typename gsBoxTopology::const_iiterator it  = topology.iBegin();it<topology.iEnd();it++)
    {
        const patchSide side1 = it->first();
        const patchSide side2 = it->second();
        m_artificialEdges[side1.patch].push_back(std::pair<patchSide,patchSide>(side1, side2));
        m_artificialEdges[side2.patch].push_back(std::pair<patchSide,patchSide>(side2, side1));

        m_artificialDofsPerSide[side1.patch][side1.index()]=m_multiBasis->basis(side2.patch).boundary(side2.index()).rows();
        m_artificialDofsPerSide[side2.patch][side2.index()]=m_multiBasis->basis(side1.patch).boundary(side1.index()).rows();

        artificialDofs[side1.patch] += m_multiBasis->basis(side2.patch).boundary(side2.index()).rows();
        artificialDofs[side2.patch] += m_multiBasis->basis(side1.patch).boundary(side1.index()).rows();

        localSpaceDofs[side1.patch] += m_multiBasis->basis(side2.patch).boundary(side2.index()).rows();
        localSpaceDofs[side2.patch] += m_multiBasis->basis(side1.patch).boundary(side1.index()).rows();
    }

    dGMapper = gsDofMapper(localSpaceDofs);

    // Match the "real" and "artificial" dofs
    for(typename gsBoxTopology::const_iiterator it  = topology.iBegin();it<topology.iEnd();it++)
    {
        const patchSide side1 = it->first();
        const patchSide side2 = it->second();

        gsMatrix<index_t> bnd1 = m_multiBasis->basis(side1.patch).boundary(side1.index());
        gsMatrix<index_t> bnd2 = m_multiBasis->basis(side2.patch).boundary(side2.index());

        for (index_t i = 0; i<bnd1.rows(); ++i)
            dGMapper.matchDof(side1.patch, (bnd1)(i,0), side2.patch,m_dofMapperGlobal.patchSize(side2.patch)+dgOffset(side2.patch, side2.side())+i);

        for (index_t i = 0; i<bnd2.rows(); ++i)
            dGMapper.matchDof(side2.patch, (bnd2)(i,0), side1.patch, m_dofMapperGlobal.patchSize(side1.patch)+dgOffset(side1.patch, side1.side())+i);
    }

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t patchDofs = m_dofMapperGlobal.patchSize(k);
        localdGmappers[k].setIdentity(1,patchDofs+artificialDofs[k]);

        for (index_t i=0; i<patchDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx)){
                localdGmappers[k].eliminateDof(i,0);
                dGMapper.eliminateDof(i, k);
            }
        }

        for (size_t n=0; n<m_artificialEdges[k].size(); ++n)
        {
            const patchSide s = m_artificialEdges[k][n].second;
            const patchSide t = m_artificialEdges[k][n].first;
            const index_t additionalDofs = m_multiBasis->basis(s.patch).boundary(s.index()).rows();
            for (index_t i = 0; i < additionalDofs; ++i) {
                const index_t idx = m_dofMapperGlobal.index(m_multiBasis->basis(s.patch).boundary(s.index())(i,0),s.patch);
                if (m_dofMapperGlobal.is_boundary_index(idx)){
                    localdGmappers[k].eliminateDof(patchDofs + dgOffset(k, t.side()) + i,0);
                    dGMapper.eliminateDof(patchDofs + dgOffset(k, t.side()) + i,k);
                }
            }
        }
        localdGmappers[k].finalize();
    }

    dGMapper.finalize();
    m_dofMapperGlobal = dGMapper;
    m_dofMapperLocal = localdGmappers;
}

template <class T>
void gsIetiMapper<T>::registerArtificialIface(const index_t patch, const boxSide s1, const patchSide ps2)
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    if(!(m_status&32))
        dgInit();

    if(!(m_status&64))
        m_status |= 64;



    gsDofMapper newlocal;
    gsVector<index_t> artificialDofs(1); artificialDofs.setZero();


    m_artificialEdges[patch].push_back(std::pair<patchSide,patchSide>(patchSide(patch, s1), ps2));
    m_artificialDofsPerSide[patch][s1.index()]=m_multiBasis->basis(ps2.patch).boundary(ps2.index()).rows();
    artificialDofs(0) += m_multiBasis->basis(ps2.patch).boundary(ps2.index()).rows();

    const index_t patchDofs = m_dofMapperLocal[patch].patchSize(0);
    newlocal.setIdentity(1,patchDofs+artificialDofs(0));

    for (index_t i=0; i<m_multiBasis->basis(patch).size(); ++i)
    {
        const index_t idx = m_dofMapperGlobal.index(i,patch);
        if (m_dofMapperGlobal.is_boundary_index(idx))
            newlocal.eliminateDof(i,0);
    }

    for (size_t n=0; n<m_artificialEdges[patch].size(); ++n)
    {
        const patchSide s = m_artificialEdges[patch][n].second;
        const patchSide t = m_artificialEdges[patch][n].first;
        const index_t additionalDofs = m_multiBasis->basis(s.patch).boundary(s.index()).rows();
        for (index_t i = 0; i < additionalDofs; ++i) {
            const index_t idx = m_dofMapperGlobal.index(m_multiBasis->basis(s.patch).boundary(s.index())(i,0),s.patch);
            if (m_dofMapperGlobal.is_boundary_index(idx)){
                newlocal.eliminateDof(m_multiBasis->basis(patch).size() + dgOffset(patch, t.side()) + i,0);
            }
        }
    }

    m_dofMapperLocal[patch] = newlocal;
}

template <class T>
void gsIetiMapper<T>::unfinalizeMapper() {
    index_t nPatches = m_dofMapperGlobal.numPatches();
    gsVector<index_t> artificialDofs(nPatches); artificialDofs.setZero();
    gsDofMapper dGMapper;

    gsVector<index_t> localSpaceDofs(nPatches); localSpaceDofs.setZero();
    for (index_t k = 0; k < nPatches; ++k) {
        m_dofMapperLocal[k].finalize();
        localSpaceDofs[k] += m_dofMapperLocal[k].patchSize(0);
    }

    dGMapper = gsDofMapper(localSpaceDofs);

    const gsBoxTopology& topology = m_multiBasis->topology();
    // Match the "real" and "artificial" dofs
    for(typename gsBoxTopology::const_iiterator it  = topology.iBegin();it<topology.iEnd();it++)
    {
        const patchSide side1 = it->first();
        const patchSide side2 = it->second();

        gsMatrix<index_t> bnd1 = m_multiBasis->basis(side1.patch).boundary(side1.index());
        gsMatrix<index_t> bnd2 = m_multiBasis->basis(side2.patch).boundary(side2.index());

        for (index_t i = 0; i<bnd1.rows(); ++i)
            dGMapper.matchDof(side1.patch, (bnd1)(i,0), side2.patch,m_dofMapperGlobal.patchSize(side2.patch)+dgOffset(side2.patch, side2.side())+i);

        for (index_t i = 0; i<bnd2.rows(); ++i)
            dGMapper.matchDof(side2.patch, (bnd2)(i,0), side1.patch, m_dofMapperGlobal.patchSize(side1.patch)+dgOffset(side1.patch, side1.side())+i);
    }

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t patchDofs = m_dofMapperGlobal.patchSize(k);
        for (index_t i=0; i<patchDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx)){
                dGMapper.eliminateDof(i, k);
            }
        }

        for (size_t n=0; n<m_artificialEdges[k].size(); ++n)
        {
            const patchSide s = m_artificialEdges[k][n].second;
            const patchSide t = m_artificialEdges[k][n].first;
            const index_t additionalDofs = m_multiBasis->basis(s.patch).boundary(s.index()).rows();
            for (index_t i = 0; i < additionalDofs; ++i) {
                const index_t idx = m_dofMapperGlobal.index(m_multiBasis->basis(s.patch).boundary(s.index())(i,0),s.patch);
                if (m_dofMapperGlobal.is_boundary_index(idx)){
                    dGMapper.eliminateDof(patchDofs + dgOffset(k, t.side()) + i,k);
                }
            }
        }
    }

    dGMapper.finalize();
    m_dofMapperGlobal = dGMapper;
}

template <class T>
void gsIetiMapper<T>::dgInit()
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( !(m_status&32), "gsIetiMapper: The dG components have already been initialized" );
    m_status |= 32;

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    m_artificialEdges.clear();
    m_artificialEdges.resize(nPatches);
    m_artificialDofsPerSide.clear();
    m_artificialDofsPerSide.resize(nPatches);
    for (index_t k = 0; k < nPatches; ++k)
        m_artificialDofsPerSide[k] = {0,0,0,0,0};

}
} // namespace gismo
