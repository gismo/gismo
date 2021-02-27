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
 *       (m_status&2)!=0    means that artificial interfaces have been registered
 *       (m_status&4)!=0    means that setupMappers has been called
 *       (m_status&8)!=0    means that the jump matrices have been computed
 *       (m_status&16)!=0   means that corners have been set up as primal constraints
 *       (m_status&flag)!=0 for flag = 32, 64,... means that edges, faces, ... have been set up as primal constraints
 */

namespace gismo
{

template <class T>
void gsIetiMapper<T>::init(
        const gsMultiBasis<T>& multiBasis,
        gsDofMapper dofMapperGlobal,
        Matrix fixedPart
    )
{
    GISMO_ASSERT( dofMapperGlobal.componentsSize() == 1, "gsIetiMapper<T>::init: "
        "Got only 1 multi basis, so a gsDofMapper with only 1 component is expected." );
    const index_t nPatches = dofMapperGlobal.numPatches();
    m_multiBasis = &multiBasis;
    m_dofMapperGlobal = give(dofMapperGlobal);
    m_dofMapperLocal.clear();
    m_fixedPartGlobal = give(fixedPart);
    m_fixedPart.clear();
    m_jumpMatrices.clear();
    m_nPrimalDofs = 0;
    m_primalConstraints.clear();
    m_primalConstraints.resize(nPatches);
    m_primalDofIndices.clear();
    m_primalDofIndices.resize(nPatches);
    m_artificialIfaces.clear();
    m_artificialIfaces.resize(nPatches);
    m_status = 1;
}

template <class T>
void gsIetiMapper<T>::setupMappers()
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( !(m_status&4), "gsIetiMapper:: setupMappers has already been called." );
    m_status |= 4;

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    m_dofMapperLocal.reserve(nPatches);
    m_fixedPart.reserve(nPatches);
    for (index_t k=0; k<nPatches; ++k)
    {
        // The first patch is the k th patch. The following "patches" are the
        // artificial interfaces.
        const index_t nDofs = m_dofMapperGlobal.patchSize(k);
        const index_t nArtIf = m_artificialIfaces[k].size();
        gsVector<index_t> sizes(1+nArtIf);
        sizes[0] = nDofs;
        for (index_t l=0; l<nArtIf; ++l)
            sizes[l+1] = m_artificialIfaces[k][l].ifaceIndices.rows();
        m_dofMapperLocal.push_back(gsDofMapper(sizes));

        // Eliminate boundary dofs (we do not consider the full floating case).
        for (index_t i=0; i<nDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx))
                m_dofMapperLocal[k].eliminateDof(i,0);
        }
        for (index_t l=0; l<nArtIf; ++l)
        {
            for (index_t i=0; i<m_artificialIfaces[k][l].ifaceIndices.rows(); ++i)
            {
                const index_t kk = m_artificialIfaces[k][l].artificialIface.patch;
                const index_t ii = m_artificialIfaces[k][l].ifaceIndices[i];
                const index_t idx = m_dofMapperGlobal.index(ii,kk);
                if (m_dofMapperGlobal.is_boundary_index(idx))
                    m_dofMapperLocal[k].eliminateDof(i,l+1);
            }
        }

        m_dofMapperLocal[k].finalize();

        const index_t szFixedPart = m_dofMapperLocal[k].boundarySize();
        m_fixedPart.push_back(Matrix(szFixedPart,1));
        m_fixedPart[k].setZero();
        for (index_t i=0; i<nDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx))
            {
                const index_t globalBoundaryIdx = m_dofMapperGlobal.bindex(i,k);
                const index_t localBoundaryIdx = m_dofMapperLocal[k].bindex(i,0);
                m_fixedPart[k](localBoundaryIdx,0) = m_fixedPartGlobal(globalBoundaryIdx,0);
            }
        }
        for (index_t l=0; l<nArtIf; ++l)
        {
            for (index_t i=0; i<m_artificialIfaces[k][l].ifaceIndices.rows(); ++i)
            {
                const index_t kk = m_artificialIfaces[k][l].artificialIface.patch;
                const index_t ii = m_artificialIfaces[k][l].ifaceIndices[i];
                const index_t idx = m_dofMapperGlobal.index(ii,kk);
                if (m_dofMapperGlobal.is_boundary_index(idx))
                {
                    const index_t globalBoundaryIdx = m_dofMapperGlobal.bindex(ii,kk);
                    const index_t localBoundaryIdx = m_dofMapperLocal[k].bindex(i,l+1);
                    m_fixedPart[k](localBoundaryIdx,0) = m_fixedPartGlobal(globalBoundaryIdx,0);
                }
            }
        }
    }
    // Cleanup
    m_fixedPartGlobal.resize(0,0);
}


template <class T>
typename gsIetiMapper<T>::Matrix
gsIetiMapper<T>::constructGlobalSolutionFromLocalSolutions( const std::vector<Matrix>& localContribs )
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );

    if (!(m_status&4))
        setupMappers();

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    GISMO_ASSERT( nPatches == static_cast<index_t>(localContribs.size()),
        "gsIetiMapper::constructGlobalSolutionFromLocalSolutions; The number of local contributions does "
        "not argee with the number of patches." );

    Matrix result;
    result.setZero( m_dofMapperGlobal.freeSize(), localContribs[0].cols() );

    // We are never extracting the solution from artificial interfaces
    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t sz=m_dofMapperLocal[k].size();
        for (index_t i=0; i<sz; ++i)
        {
            // There is an asignment. This means that if there are several values, we just take the last one.
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

template<typename Container, typename Element>
index_t indexOf( const Container& c, Element e )
{
    const Element* it = std::lower_bound(c.begin(),c.end(),e);
    if (it==c.end() || *it>e) return -1;
    return it-c.begin();
}
}

template <class T>
void gsIetiMapper<T>::cornersAsPrimals()
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( !(m_status&16), "gsIetiMapper::cornersAsPrimals: This function has already been called." );
    m_status |= 16;

    if (!(m_status&4))
        setupMappers();

    const index_t nPatches = m_dofMapperLocal.size();

    // Construct all corners
    std::vector<dof_helper> corners;
    const index_t dim = m_multiBasis->dim();
    corners.reserve((1<<dim)*nPatches);
    // Add corners on all patches
    for (index_t k=0; k<nPatches; ++k)
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
    // Are there corners on the artificial interfaces as well?
    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t nArtIfaces = m_artificialIfaces[k].size();
        for (index_t l=0; l<nArtIfaces; ++l)
        {
            const index_t kk = m_artificialIfaces[k][l].artificialIface.patch;
            std::vector<boxCorner> cornersOnSide;
            m_artificialIfaces[k][l].artificialIface.getContainedCorners(dim,cornersOnSide);
            for (index_t i=0; i<cornersOnSide.size(); ++i)
            {
                const boxCorner& it = cornersOnSide[i];
                const index_t idx = (*m_multiBasis)[kk].functionAtCorner(it);
                const index_t idxOnIface = indexOf(m_artificialIfaces[k][l].ifaceIndices, idx);
                if (idxOnIface > -1)
                {
                    dof_helper dh;
                    dh.globalIndex = m_dofMapperGlobal.index( idx, kk );
                    dh.patch = k;
                    dh.localIndex = m_dofMapperLocal[k].index(idxOnIface, l+1);
                    if (m_dofMapperGlobal.is_free_index(dh.globalIndex))
                    //if (m_dofMapperLocal[kk].is_free(idx,0))
                    //if (m_dofMapperLocal[k].is_free(idxOnIface,l+1))

                        corners.push_back(dh);
                }
            }
        }
    }
    // Sort corners to collapse corners with same global index
    std::sort(corners.begin(), corners.end());

    // Create data
    index_t lastIndex=-1;
    const index_t sz = corners.size();
    for (index_t i=0; i<sz; ++i)
    {
        if (lastIndex!=corners[i].globalIndex)
        {
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

    GISMO_ASSERT( !(m_status&2), "gsIetiMapper::interfaceAveragesAsPrimals "
        "is not implemented for artificial ifaces." ); //TODO

    const unsigned flag = 1<<(4+d);
    GISMO_ASSERT( !(m_status&flag), "gsIetiMapper::interfaceAveragesAsPrimals: This function has "
        " already been called for d="<<d );
    m_status |= flag;

    if (!(m_status&4))
        setupMappers();

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

    if (!(m_status&4))
        setupMappers();

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
std::vector<index_t> gsIetiMapper<T>::skeletonDofs(const index_t patch) const
{
    GISMO_ASSERT( m_status&4, "gsIetiMapper::skeletonDofs can only be called after "
        "jump matrices have been computed." );

    std::vector<index_t> result;
    const index_t patchSize = m_dofMapperGlobal.patchSize(patch);
    const index_t dim = m_multiBasis->dim();
    result.reserve(2*dim*std::pow(patchSize,(1.0-dim)/dim));
    // All dofs on patch that are coupled
    for (index_t i=0; i<patchSize; ++i)
        if ( m_dofMapperGlobal.is_coupled(i,patch) )
            result.push_back( m_dofMapperLocal[patch].index(i,0) );
    // All dofs on artificial iface
    for (index_t l=0; l<m_artificialIfaces[patch].size(); ++l)
    {
        const index_t nDofs = m_artificialIfaces[patch][l].ifaceIndices.size();
        for (index_t i=0; i<nDofs; ++i)
            if ( m_dofMapperLocal[patch].is_free(i,l+1) )
                result.push_back( m_dofMapperLocal[patch].index(i,l+1) );
    }
    // All dofs on patch that are part of artificial iface of other patch
    for (index_t kk=0; kk<m_artificialIfaces.size(); ++kk)
    {
        for (index_t l=0; l<m_artificialIfaces[kk].size(); ++l)
        {
            // Does that artificial interface refer back?
            if ( patch == m_artificialIfaces[kk][l].artificialIface.patch )
            {
                const index_t nDofs = m_artificialIfaces[kk][l].ifaceIndices.size();
                for (index_t i=0; i<nDofs; ++i)
                {
                    const index_t idx = m_artificialIfaces[kk][l].ifaceIndices[i];
                    if ( m_dofMapperLocal[patch].is_free(idx,0) )
                        result.push_back(m_dofMapperLocal[patch].index(idx,0));
                }
            }
        }
    }
    return result;
}

template <class T>
void gsIetiMapper<T>::computeJumpMatrices(bool fullyRedundant, bool excludeCorners)
{
    GISMO_ASSERT( m_status&1, "gsIetiMapper: The class has not been initialized." );
    GISMO_ASSERT( !(m_status&8), "gsIetiMapper::computeJumpMatrices: This function has already been called." );
    m_status |= 8;

    if (!(m_status&4))
        setupMappers();

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    const index_t coupledSize = m_dofMapperGlobal.coupledSize();

    std::vector< std::vector< std::pair<index_t,index_t> > > coupling;
    coupling.resize(coupledSize);

    // Find the groups of to be coupled indices in the global mapper
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
    // Additionally, the artificial ifaces are coupled with the iface from
    // the patch where they were taken from
    {
        index_t reserve = 0;
        for (index_t k=0; k<nPatches; ++k)
            for (index_t l=0; l<m_artificialIfaces[k].size(); ++l)
                reserve += m_artificialIfaces[k][l].ifaceIndices.size();
        coupling.resize(coupledSize+reserve);
    }
    // TODO: Think about strategy: The following joins iface with artificial iface,
    // which is basically minimally connection; never fully redundant.
    for (index_t k=0; k<nPatches; ++k)
    {
        for (index_t l=0; l<m_artificialIfaces[k].size(); ++l)
        {
            const index_t nDofs = m_artificialIfaces[k][l].ifaceIndices.size();
            for (index_t i=0; i<nDofs; ++i)
            {
                std::vector< std::pair<index_t,index_t> > coupling_pair(2);
                // Binding artificial iface
                coupling_pair[0].first = k;
                coupling_pair[0].second = m_dofMapperLocal[k].index(i,l+1);
                // and corresponding real iface
                const index_t idx = m_artificialIfaces[k][l].ifaceIndices[i];
                const index_t kk = m_artificialIfaces[k][l].artificialIface.patch;
                coupling_pair[1].first = kk;
                coupling_pair[1].second = m_dofMapperLocal[kk].index(idx,0);
                coupling.push_back( give(coupling_pair) );
            }
        }
    }

    // Erease data for corners if so desired
    if (excludeCorners)
    {
        GISMO_ASSERT( !(m_status&2), "gsIetiMapper::computeJumpMatrices: "
            "The option excludeCorners is not implemented for artificial ifaces." ); //TODO
        const index_t dim = m_multiBasis->dim();
        for (index_t k=0; k<nPatches; ++k)
        {
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
    {
        const index_t dim = m_multiBasis->dim();
        jumpMatrices_se[i].reserve(std::pow(m_dofMapperLocal[i].freeSize(),(1.0-dim)/dim));
    }

    index_t multiplier = 0;
    for (index_t i=0; i<coupledSize; ++i)
    {
        const index_t n = coupling[i].size();
        const index_t maxIndex = fullyRedundant ? (n-1) : 1;
        for (index_t j1=0; j1<maxIndex; ++j1)
        {
            for (index_t j2=j1+1; j2<n; ++j2)
            {
                const index_t patch1 = coupling[i][j1].first;
                const index_t localIndex1 = coupling[i][j1].second;
                const index_t patch2 = coupling[i][j2].first;
                const index_t localIndex2 = coupling[i][j2].second;
                jumpMatrices_se[patch1].add(multiplier,localIndex1,(T)1);
                jumpMatrices_se[patch2].add(multiplier,localIndex2,(T)-1);
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
void gsIetiMapper<T>::registerArtificialIface(patchSide realIface, patchSide artificialIface)
{
    m_status |= 2;
    ArtificialIface ai;
    ai.realIface = realIface;
    ai.artificialIface = artificialIface;
    ai.ifaceIndices = (*m_multiBasis)[realIface.patch].boundary(artificialIface);
    m_artificialIfaces[realIface.patch].push_back(give(ai));
}

template <class T>
void gsIetiMapper<T>::registerAllArtificialIfaces()
{
    const gsBoxTopology& top = m_multiBasis->topology();
    const short_t dim        = m_multiBasis->domainDim();
    const index_t nPatches   = m_dofMapperGlobal.numPatches();
    for (index_t k=0; k<nPatches; ++k)
    {
        for (boxSide it = boxSide::getFirst(dim); it!=boxSide::getEnd(dim); ++it) {
            patchSide realIface(k,it);
            patchSide artificialIface;
            if(top.getNeighbour(realIface, artificialIface))
                registerArtificialIface(realIface, artificialIface);
        }
    }
}

} // namespace gismo
