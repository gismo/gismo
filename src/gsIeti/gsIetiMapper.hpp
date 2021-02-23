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
#include <gsAssembler/gsVisitorDg.h>

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
        const Matrix& fixedPart,
        const bool dG
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

    m_dG = dG;
    m_artificialEdges.clear();
    m_artificialEdges.resize(nPatches);
    m_artificialDofsPerSide.clear();
    m_artificialDofsPerSide.resize(nPatches);
    gsVector<index_t> artificialDofs(nPatches); artificialDofs.setZero();

    if(m_dG){
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

        m_dofMapperGlobaldG = gsDofMapper(localSpaceDofs);

        // Match the "real" and "artificial" dofs
        for(typename gsBoxTopology::const_iiterator it  = topology.iBegin();it<topology.iEnd();it++)
        {
            const patchSide side1 = it->first();
            const patchSide side2 = it->second();

            gsMatrix<index_t> bnd1 = m_multiBasis->basis(side1.patch).boundary(side1.index());
            gsMatrix<index_t> bnd2 = m_multiBasis->basis(side2.patch).boundary(side2.index());

            for (index_t i = 0; i<bnd1.rows(); ++i)
                m_dofMapperGlobaldG.matchDof(side1.patch, (bnd1)(i,0), side2.patch,m_dofMapperGlobal.patchSize(side2.patch)+dgOffset(side2.patch, side2.side())+i);

            for (index_t i = 0; i<bnd2.rows(); ++i)
                m_dofMapperGlobaldG.matchDof(side2.patch, (bnd2)(i,0), side1.patch, m_dofMapperGlobal.patchSize(side1.patch)+dgOffset(side1.patch, side1.side())+i);
        }
    }

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t patchDofs = m_dofMapperGlobal.patchSize(k);
        m_dofMapperLocal[k].setIdentity(1,patchDofs+artificialDofs[k]);

        for (index_t i=0; i<patchDofs; ++i)
        {
            const index_t idx = m_dofMapperGlobal.index(i,k);
            if (m_dofMapperGlobal.is_boundary_index(idx)){
                m_dofMapperLocal[k].eliminateDof(i,0);
                m_dofMapperGlobaldG.eliminateDof(i, k);
            }
        }

        gsInfo << "\n patch: " << k<<"\n";
        for (size_t n=0; n<m_artificialEdges[k].size(); ++n)
        {
            const patchSide s = m_artificialEdges[k][n].second;
            const patchSide t = m_artificialEdges[k][n].first;
            const index_t additionalDofs = m_multiBasis->basis(s.patch).boundary(s.index()).rows();
            for (index_t i = 0; i < additionalDofs; ++i) {
                const index_t idx = m_dofMapperGlobal.index(m_multiBasis->basis(s.patch).boundary(s.index())(i,0),s.patch);
                if (m_dofMapperGlobal.is_boundary_index(idx)){
                    m_dofMapperLocal[k].eliminateDof(patchDofs + dgOffset(k, t.side()) + i,0);
                    m_dofMapperGlobaldG.eliminateDof(patchDofs + dgOffset(k, t.side()) + i,k);
                }
            }
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

    if(m_dG)
        m_dofMapperGlobaldG.finalize();


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
    std::vector<index_t> localIndex;
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
    std::vector<dof_helper> corners;
    corners.reserve(4*nPatches);
    const short_t dim = m_multiBasis->dim();
    for (index_t k=0; k<nPatches; ++k){
        for (boxCorner it = boxCorner::getFirst(dim); it!=boxCorner::getEnd(dim); ++it)
        {
            const index_t idx = (*m_multiBasis)[k].functionAtCorner(it);
            dof_helper dh;
            dh.globalIndex = m_dofMapperGlobaldG.index( idx, k );
            dh.patch = k;
            dh.localIndex.push_back(m_dofMapperLocal[k].index( idx, 0 ));
            if (m_dofMapperGlobaldG.is_free_index(dh.globalIndex))
            {
                std::vector<boxSide> side;
                it.getContainingSides(dim, side);
                for (typename std::vector<boxSide>::iterator  bit = side.begin();  bit != side.end(); ++bit) {
                    std::vector<std::pair<patchSide,patchSide> >::iterator it = std::find_if( m_artificialEdges[k].begin(), m_artificialEdges[k].end(),
                                                                                            [bit](const std::pair<patchSide, patchSide>& element){ return element.first.side() == *bit;} );

                    if(it != m_artificialEdges[k].end()) {
                        patchSide artificialEdge = it->second;
                        //index_t artificial_corner = m_dofMapperGlobal.patchSize(k) + dgOffset(artificialEdge.patch, artificialEdge.side());
                        index_t artificial_corner = m_dofMapperGlobal.patchSize(k) + dgOffset(k, it->first.side());

                        gsMatrix<index_t> bndBasis = m_multiBasis->basis(artificialEdge.patch).boundary(artificialEdge);
                        const index_t candidate = bndBasis(bndBasis.rows()-1,0);
                        index_t c = m_dofMapperGlobaldG.index( candidate, artificialEdge.patch );
                        std::vector<std::pair<index_t, index_t> > patchDofs;
                        m_dofMapperGlobaldG.preImage(c, patchDofs);

                        for (size_t n = 0; n < patchDofs.size(); ++n)
                            if (patchDofs[n].first == k)
                                if(m_dofMapperLocal[k].index( patchDofs[n].second,0) < m_dofMapperLocal[k].freeSize()){
                                    artificial_corner += bndBasis.rows() - 1;
                                    break;
                                }

                        // The corner can only be the first or the last entry of the boundary basis
                        dh.localIndex.push_back((m_dofMapperLocal[k].index( artificial_corner,0)));
                    }

                }
                corners.push_back(dh);
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
        const index_t localIndex  = corners[i].localIndex[0];

        SparseVector constr(m_dofMapperLocal[patch].freeSize());
        constr[localIndex] = 1;
        m_primalConstraints[patch].push_back(give(constr));
        m_primalDofIndices[patch].push_back(cornerIndex);
    }

    // .... and then the "artificial" ones
    for (index_t i=0;i<sz; ++i){
        const index_t patch = corners[i].patch;
        for (size_t  j = 1; j < corners[i].localIndex.size(); ++j) {
            SparseVector constr(m_dofMapperLocal[patch].freeSize());
            const index_t localIndex = corners[i].localIndex[j];
            constr[localIndex] = 1;
            m_primalConstraints[patch].push_back(give(constr));
        }
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

    const index_t nPatches = m_dofMapperGlobal.numPatches();
    const index_t coupledSize = m_dofMapperGlobaldG.coupledSize();

    // Find the groups of to be coupled indices
    std::vector< std::vector< std::pair<index_t,index_t> > > coupling;
    coupling.resize(coupledSize);

    for (index_t k=0; k<nPatches; ++k)
    {
        const index_t patchSize = m_dofMapperGlobaldG.patchSize(k);
        for (index_t i=0; i<patchSize; ++i)
        {
            const index_t globalIndex = m_dofMapperGlobaldG.index(i,k);
            if ( m_dofMapperGlobaldG.is_coupled_index(globalIndex) )
            {
                const index_t coupledIndex = m_dofMapperGlobaldG.cindex(i,k);
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
                const index_t globalIndex = m_dofMapperGlobaldG.index(idx,k);
                if ( m_dofMapperGlobaldG.is_coupled_index(globalIndex) )
                {
                    const index_t coupledIndex = m_dofMapperGlobaldG.cindex(idx,k);
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
void gsIetiMapper<T>::adddGInterfaceContributions(gsAssembler<T>* localassembler, gsSparseMatrix<T>& localmatrix, gsMatrix<T>& localrhs,
                                                  const index_t patch, const gsMultiPatch<T>& domain, const gsOptionList options) const
{
    GISMO_ASSERT(m_dG, "gsIETIMapper::adddGInterfaceContributions: this is only allowed in dG mode");
    gsMatrix<index_t> *actives1, *actives2;
    gsMatrix<index_t> activesExtra1, activesExtra2;

    gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights

    gsQuadRule<T> QuRule;

    for (size_t e = 0; e < m_artificialEdges[patch].size(); ++e)
    {
        patchSide side1, side2;
        if(m_multiBasis->basis(patch).numElements(m_artificialEdges[patch][e].second.side()) > m_multiBasis->basis(patch).numElements(m_artificialEdges[patch][e].first.side())) {
            side2 = m_artificialEdges[patch][e].first;
            side1 = m_artificialEdges[patch][e].second;
        }
        else {
            side1 = m_artificialEdges[patch][e].first;
            side2 = m_artificialEdges[patch][e].second;
        }

        gsVisitorDg<T> dg(localassembler->pde());

        boundaryInterface bi(side1, side2, domain.domainDim());
        gsRemapInterface<T> interfaceMap(domain, *m_multiBasis, bi);

        const index_t patch1      = side1.patch;
        const index_t patch2      = side2.patch;
        const gsBasis<T> & B1 = m_multiBasis->basis(patch1);
        const gsBasis<T> & B2 = m_multiBasis->basis(patch2);

        // Initialize
        dg.initialize(B1, B2, bi, options, QuRule);

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt1 = interfaceMap.makeDomainIterator();
        typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator(bi.second().side()); // for the correct penalty term in the Visitor
        typename gsBasis<T>::domainIter domExtra = B1.makeDomainIterator(bi.first().side());

        // iterate over all boundary grid cells on the "left"
        for (; domIt1->good(); domIt1->next() )
        {
            QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            dg.evaluate(B1, domain[patch1], B2, domain[patch2], quNodes1, quNodes2);

            // Assemble on element
            dg.assemble(*domExtra,*domIt2, quWeights); // the iterator is only used for the calculation of the cell size

            //extract the actives and prepare them for the IETI_locToGlob map
            dg.getActives(actives1, actives2);
            prepareActives(bi,*actives1,*actives2,activesExtra1 );


            //do the map
            dg.localToGlobalIETI(m_dofMapperLocal[patch1],m_fixedPart[patch1],
                                 activesExtra1, localmatrix, localrhs);

        }
    }

}

template<class T>
void gsIetiMapper<T>::prepareActives(const boundaryInterface & bi, gsMatrix<index_t>& actives1, gsMatrix<index_t>& actives2,
                                          gsMatrix<index_t>& activesExtra1) const {
    index_t patch1 = bi.first().patch;
    patchSide side1 = bi.first();
    patchSide side2 = bi.second();
    const index_t n1 = actives1.rows();
    const index_t n2 = actives2.rows();

    actives1.conservativeResize(n1 + n2, actives1.cols());
    actives2.conservativeResize(n2 + n1, actives2.cols());
    actives1.bottomRows(n2).setZero();
    actives2.bottomRows(n1).setZero();

    activesExtra1.resize(m_artificialDofsPerSide[patch1][side1.index()],1);
    index_t iter1 = 0;
    for (index_t i = 0; i < n2; i++) {
        index_t k = dgFindCorrespondingExtraIndex(side1, side2, actives2(i,0));
        if (k >= 0) {
            (actives1)(n1+i,0)= k;
            activesExtra1(iter1, 0) = i;
            iter1++;
        }
    }
    //shrink to the appropriate size, since not all extra dofs are active on an element!
    activesExtra1.conservativeResize(iter1, activesExtra1.cols());
}



} // namespace gismo
