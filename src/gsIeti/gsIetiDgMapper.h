/** @file gsIetiDgMapper.h

    @brief TODO

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, R. Schneckenleitner
*/

#include <gsIeti/gsIetiMapper.h>

namespace gismo {

template<class T=real_t>
class gsIetiDgMapper : protected gsIetiMapper<T> {
    typedef gsIetiMapper<T> Base;

public:
    gsIetiDgMapper( const gsMultiPatch<T>& mp, const gsMultiBasis<T>& mb, gsDofMapper dm, const gsMatrix<T>& fixedPart )
        : m_status(0), m_artificialIfaces(mb.nBases()), m_multiPatch(&mp), m_multiBasis(&mb), m_dofMapperOrig(give(dm)), m_fixedPart(fixedPart) {}

    /// Data structure for artificial interfaces, as used in dG settings
    struct ArtificialIface {
       /// The real interface (local patch to which the artificial interface is assigned+side)
       patchSide           realIface;
       /// The artificial interface (foreign patch from which the artificial interface taken+side)
       patchSide           artificialIface;
       /// The indicies of the basis of the foreign patch that belong to artificial interface
       /// Vector is assumed to be sorted
       gsVector<index_t>   ifaceIndices;
    };

    /// @brief Registers an artificial interface, as used in dG settings
    /// @param realIface        The real interface (local patch to which the artificial
    ///                         interface is assigned+side)
    /// @param artificialIface  The artificial interface (foreign patch from which the
    ///                         artificial interface taken+side)
    void registerArtificialIface(patchSide realIface, patchSide artificialIface)
    {
        GISMO_ASSERT( m_status==0, "gsIetiDgMapper: Cannot register artificial interfaces after requesting data." );
        ArtificialIface ai;
        ai.realIface = realIface;
        ai.artificialIface = artificialIface;
        ai.ifaceIndices = (*m_multiBasis)[artificialIface.patch].boundary(artificialIface);
        m_artificialIfaces[realIface.patch].push_back(give(ai));
    }

    /// Calls \ref registerArtificialIface for all interface known to the underlying box topology
    void registerAllArtificialIfaces()
    {
        GISMO_ASSERT( m_status==0, "gsIetiDgMapper: Cannot register artificial interfaces after requesting data." );
        const gsBoxTopology& top = m_multiBasis->topology();
        const short_t dim        = m_multiBasis->domainDim();
        const index_t nPatches   = m_dofMapperOrig.numPatches();
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

    template<typename Container, typename Element>
    static index_t indexOf( const Container& c, Element e ) //TODO: unused
    {
        const Element* it = std::lower_bound(c.begin(),c.end(),e);
        if (it==c.end() || *it>e) return -1;
        return it-c.begin();
    }

    /// Setup of the dof mappers
    void finalize()
    {
        GISMO_ASSERT( m_status==0, "gsIetiDgMapper: Already finalized." );
        m_status = 1;

        // TODO: we can get rid of this assumption if we incorporate dofs that
        // have been coupled into m_dofMapperOrig
        GISMO_ASSERT( m_dofMapperOrig.coupledSize() == 0, "Not implemented." );

        const index_t nPatches = m_dofMapperOrig.numPatches();

        m_dofMapperLocal2.clear();
        m_dofMapperLocal2.reserve(nPatches);

        // Subdomain = patch + artificial ifaces
        gsVector<index_t> subdomainSizes(nPatches);
        std::vector< gsVector<index_t> > sizes(nPatches);
        for (index_t k=0; k<nPatches; ++k)
        {
            const index_t nArtIf = m_artificialIfaces[k].size();
            sizes[k].resize(1+nArtIf);
            sizes[k][0]       = m_dofMapperOrig.patchSize(k);
            subdomainSizes[k] = m_dofMapperOrig.patchSize(k);
            for (index_t l=0; l<nArtIf; ++l)
            {
                sizes[k][l+1]      = m_dofMapperOrig.patchSize( m_artificialIfaces[k][l].artificialIface.patch );
                subdomainSizes[k] += m_artificialIfaces[k][l].ifaceIndices.rows();
            }
            m_dofMapperLocal2.push_back( gsDofMapper(sizes[k]) );
        }

        m_dofMapperMod = gsDofMapper(subdomainSizes);

        for (index_t k=0; k<nPatches; ++k)
        {
            const index_t patchDofs = m_dofMapperOrig.patchSize(k);
            // Eliminate boundary dofs
            for (index_t i=0; i<patchDofs; ++i)
            {
                const index_t idx = m_dofMapperOrig.index(i,k);
                if (m_dofMapperOrig.is_boundary_index(idx))
                {
                    m_dofMapperMod.eliminateDof(i,k);
                    m_dofMapperLocal2[k].eliminateDof(i,0);
                }
            }
            // Match or eliminate dofs on artificial interfaces
            const index_t nArtIf = m_artificialIfaces[k].size();
            {
                index_t idx = patchDofs;
                for (index_t l=0; l<nArtIf; ++l)
                {
                    for (index_t i=0; i<m_artificialIfaces[k][l].ifaceIndices.rows(); ++i)
                    {
                        const index_t kk = m_artificialIfaces[k][l].artificialIface.patch;
                        const index_t ii = m_artificialIfaces[k][l].ifaceIndices[i];
                        m_dofMapperMod.matchDof( k,idx, kk,ii );

                        const index_t global_idx = m_dofMapperOrig.index(ii,kk);
                        if (m_dofMapperOrig.is_boundary_index(global_idx))
                        {
                            //m_dofMapperMod.eliminateDof(idx,k); // TODO: necessary?
                            m_dofMapperLocal2[k].eliminateDof(ii,l+1);
                        }
                        ++idx;
                    }
                }
            }
            // Eliminate all dofs on m_dofMapperLocal2 that are not on the artificial interface
            for (index_t l=0; l<nArtIf; ++l)
            {
                index_t last_idx = -1;
                for (index_t i=0; i<m_artificialIfaces[k][l].ifaceIndices.rows(); ++i)
                {
                    const index_t idx = m_artificialIfaces[k][l].ifaceIndices[i];
                    GISMO_ASSERT (idx>last_idx, "Not sorted; k="<<k<<"; l="<<l);
                    for (index_t j=last_idx+1; j<idx; ++j)
                        m_dofMapperLocal2[k].eliminateDof(j,l+1);
                    last_idx = idx;
                }
                GISMO_ASSERT (sizes[k][l+1]>last_idx, "Not sorted; k="<<k<<"; l="<<l);
                for (index_t j=last_idx+1; j<sizes[k][l+1]; ++j)
                    m_dofMapperLocal2[k].eliminateDof(j,l+1);
            }
            m_dofMapperLocal2[k].finalize();
        }
        m_dofMapperMod.finalize();

        Base::init( *m_multiBasis, dofMapperMod(), m_fixedPart );

        // Populate Base::m_artificialDofInfo
        /*const index_t nDofs = Base::m_dofMapperGlobal.freeSize();
        gsMatrix<index_t> dofs(nDofs,2);   // Has information (patch, localIndex) for each global dof
        dofs.setZero();

        for (index_t k=0; k<nPatches; ++k)
        {
            // Here we only consider the real values
            const index_t sz = Base::m_multiBasis->piece(k).size();
            for (index_t i=0; i<sz; ++i)
            {
                const index_t globalIndex = Base::m_dofMapperGlobal.index(i,k);
                if (Base::m_dofMapperGlobal.is_free_index(globalIndex))
                {
                    GISMO_ASSERT( dofs(globalIndex,0) == 0, "Internal error.");
                    dofs(globalIndex,0) = k;
                    dofs(globalIndex,1) = Base::m_dofMapperLocal2[k].index(i,0) + 1;
                }
            }
        }

        Base::m_artificialDofInfo.resize( nPatches );
        for (index_t k=0; k<nPatches; ++k)
        {
            const index_t sz  = Base::m_multiBasis->piece(k).size();
            const index_t sz2 = Base::m_dofMapperGlobal.patchSize(k);
            for (index_t i=sz; i<sz2; ++i)
            {
                const index_t globalIndex       = Base::m_dofMapperGlobal.index(i,k);
                if (Base::m_dofMapperGlobal.is_free_index(globalIndex))
                {
                    const index_t otherPatch        = dofs(globalIndex,0);
                    const index_t indexOnOtherPatch = dofs(globalIndex,1) - 1;
                    GISMO_ASSERT( indexOnOtherPatch>=0, "Internal error." );
                    gsVector<index_t> & which = Base::m_artificialDofInfo[otherPatch][k];
                    if (which.rows() == 0)
                        which.setZero( Base::m_dofMapperLocal[otherPatch].freeSize(), 1 );
                    which[indexOnOtherPatch] = i + 1;
                }
            }
        }*/

    }

    /// @brief Returns artificial interfaces for given patch
    const std::vector<ArtificialIface>& artificialIfaces(index_t k)
    {
        if (!m_status) finalize();
        return m_artificialIfaces[k];
    }

    /// @brief The local dof mappers
    const std::vector<gsDofMapper>& dofMappersLocal()
    {
        if (!m_status) finalize();
        return m_dofMapperLocal2;
    }

    const gsDofMapper& dofMapperLocal(index_t k)
    {
        if (!m_status) finalize();
        return m_dofMapperLocal2[k];
    }


    gsMultiPatch<T> multiPatchLocal(index_t patch)
    {
        if (!m_status) finalize();

        const std::vector<ArtificialIface>& artIfaces = m_artificialIfaces[patch];
        std::vector<gsGeometry<T>*> mp_local_data;
        mp_local_data.reserve(1+2*m_multiPatch->geoDim());
        mp_local_data.push_back((*m_multiPatch)[patch].clone().release());
        for (size_t i=0; i<artIfaces.size(); ++i)
             mp_local_data.push_back((*m_multiPatch)[artIfaces[i].artificialIface.patch].clone().release());
        gsMultiPatch<T> mp_local(mp_local_data);
        mp_local.computeTopology();
        return mp_local;
    }

    gsMultiBasis<T> multiBasisLocal(index_t patch)
    {
        if (!m_status) finalize();

        const std::vector<ArtificialIface>& artIfaces = m_artificialIfaces[patch];
        std::vector<gsBasis<T>*> mb_local_data;
        mb_local_data.reserve(1+2*m_multiPatch->geoDim());
        mb_local_data.push_back((*m_multiBasis)[patch].clone().release());
        for (size_t i=0; i<artIfaces.size(); ++i)
             mb_local_data.push_back((*m_multiBasis)[artIfaces[i].artificialIface.patch].clone().release());
        gsMultiBasis<T> mb_local(mb_local_data,multiPatchLocal(patch));
        return mb_local;
    }

    const gsDofMapper& dofMapperMod()
    {
        if (!m_status) finalize();
        return m_dofMapperMod;
    }

    gsMatrix<T> fixedPart(index_t k)
    {
        gsMatrix<T> result(dofMapperLocal(k).boundarySize(),1);
        result.setZero();
        result.topRows(Base::fixedPart(k).rows()) = Base::fixedPart(k); // TODO
        return result;
    }

    //using Base::cornersAsPrimals;
    //using Base::interfaceAveragesAsPrimals;
    using Base::customPrimalConstraints;
    using Base::computeJumpMatrices;
    using Base::skeletonDofs;
    using Base::initFeSpace;
    using Base::constructGlobalSolutionFromLocalSolutions;
    using Base::nLagrangeMultipliers;
    using Base::nPrimalDofs;
    using Base::primalConstraints;
    using Base::primalDofIndices;
    using Base::jumpMatrix;
    //using Base::dofMapperGlobal;
    //using Base::dofMapperLocal;
    //using Base::fixedPart;
    //using Base::multiBasis;

    const gsMultiBasis<T>& multiBasis() const { return *m_multiBasis; }
    const gsMultiBasis<T>& multiPatch() const { return *m_multiPatch; }

    void cornersAsPrimals( bool includeIsolated=false )
    {
        const index_t primals_old = nPrimalDofs();
        Base::cornersAsPrimals(true);
        const index_t primals_new = nPrimalDofs();
        transferToAi(primals_old,primals_new);
        if (!includeIsolated) removeIsolated(primals_old,primals_new);
    }

    void interfaceAveragesAsPrimals( const gsMultiPatch<T>& geo, short_t d, bool includeIsolated=false )
    {
        const index_t primals_old = nPrimalDofs();
        Base::interfaceAveragesAsPrimals(geo,d,true);
        const index_t primals_new = nPrimalDofs();
        transferToAi(primals_old,primals_new);
        if (!includeIsolated) removeIsolated(primals_old,primals_new);
    }

    void removeIsolated(index_t primals_old, index_t primals_new)
    {
        std::vector<index_t> count(primals_new);
        std::vector< std::vector<index_t> >& primalDofIndices = Base::m_primalDofIndices;
        for (index_t k=0; k<primalDofIndices.size(); ++k)
            for (index_t l=0; l<primalDofIndices[k].size(); ++l)
                if (primalDofIndices[k][l]<primals_new)
                    count[primalDofIndices[k][l]]++;

        std::vector<index_t> newIndices(primals_new);
        for (index_t i=0, j=0; i<primals_new; ++i)
        {
            if (i>=primals_old && count[i]==1)
                newIndices[i] = -1;
            else
            {
                newIndices[i] = j;
                ++j;
            }
        }

        for (index_t k=0; k<primalDofIndices.size(); ++k)
            for (index_t l=0; l<primalDofIndices[k].size(); ++l)
            {
                const index_t newIndex = newIndices[primalDofIndices[k][l]];
                if (newIndex==-1)
                {
                    gsInfo << "Removed a constraint from patch "<<k<<"\n";
                    primalDofIndices[k].erase(primalDofIndices[k].begin()+l);
                    Base::m_primalConstraints[k].erase(Base::m_primalConstraints[k].begin()+l);
                    --l;
                }
                else
                {
                    primalDofIndices[k][l] = newIndex;
                }
            }
    }


protected:
    void transferToAi(index_t primals_old, index_t primals_new )
    {
        for (size_t k=0; k<m_artificialIfaces.size(); ++k)
            for (size_t i=0; i<m_artificialIfaces[k].size(); ++i)
            {
                ArtificialIface& ai = m_artificialIfaces[k][i];
                /// The real interface (local patch to which the artificial interface is assigned+side)
                //ai.realIface;
                /// The artificial interface (foreign patch from which the artificial interface taken+side)
                //ai.artificialIface;
                /// The indicies of the basis of the foreign patch that belong to artificial interface
                /// Vector is assumed to be sorted
                //ai.ifaceIndices; //gsVector<index_t>
                const index_t patch = ai.artificialIface.patch;
                for (size_t j=0; j<Base::m_primalDofIndices[patch].size(); ++j)
                {
                    if (primals_old<=Base::m_primalDofIndices[patch][j]
                        && Base::m_primalDofIndices[patch][j]<primals_new)
                    {
                        gsInfo <<"Constraint #"<<j<<" of patch "<<patch<<" (global primal dof # "
                                <<Base::m_primalDofIndices[patch][j]<<") might needed to be added to patch # "
                                <<k<<"(=="<<ai.realIface.patch<<")\n";

                        std::vector<index_t> ifaceIndicesMapped(ai.ifaceIndices.size());
                        for (index_t ii=0; ii<ai.ifaceIndices.size(); ++ii)
                        {
                            ifaceIndicesMapped[ii] = m_dofMapperLocal2[patch].index(ai.ifaceIndices[ii],0/**/);
                        }

                        bool b = true;
                        for (typename gsSparseVector<T>::InnerIterator it(Base::m_primalConstraints[patch][j],0);
                            it; ++it)
                        {
                            bool contains = (std::find(ifaceIndicesMapped.begin(),ifaceIndicesMapped.end(),it.row()) != ifaceIndicesMapped.end());
                            b &= contains;
                        }
                        if (b)
                        {
                            gsInfo << "Yes it, should be added!\n";
                            /*index_t OFFSET = Base::m_dofMapperLocal2[k].freeSize();
                            for (size_t ii=i; ii<m_artificialIfaces[k].size(); ++ii)
                            {
                                OFFSET -= m_artificialIfaces[k][ii].ifaceIndices.rows();
                            }*/
                            gsSparseVector<T> newConstraint(m_dofMapperLocal2[k].freeSize());
                            for (typename gsSparseVector<T>::InnerIterator it(Base::m_primalConstraints[patch][j],0);
                                it; ++it)
                            {
                                index_t idx = std::find(ifaceIndicesMapped.begin(),ifaceIndicesMapped.end(),it.row()) - ifaceIndicesMapped.begin();

                                {
                                std::vector<std::pair<index_t,index_t> >  result;
                                m_dofMapperLocal2[patch].preImage(it.row(),result);
                                GISMO_ASSERT(result.size()==1, "result.size()="<<result.size());
                                GISMO_ASSERT(result[0].first == 0, "result[0].first="<<result[0].first);
                                idx = result[0].second; ////////
                                }

                                index_t idx_transfered = /*Self::*/m_dofMapperLocal2[k].index(idx,i+1); //TODO: is this correct?
/*   gsInfo<<"m_multiBasis->piece(k).size()="<<m_multiBasis->piece(k).size()<<"\n";
   gsInfo<<"m_dofMapperOrig.patchSize(k)="<<m_dofMapperOrig.patchSize(k)<<"\n";
   gsInfo<<"m_dofMapperMod.patchSize(k)="<<m_dofMapperMod.patchSize(k)<<"\n";
   gsInfo<<"m_dofMapperLocal2[k].patchSize(0)="<<m_dofMapperLocal2[k].patchSize(0)<<"\n";
   gsInfo<<"m_dofMapperLocal2[k].patchSize(1)="<<m_dofMapperLocal2[k].patchSize(1)<<"\n";
   gsInfo<<"m_dofMapperLocal2[k].patchSize(2)="<<m_dofMapperLocal2[k].patchSize(2)<<"\n";
   gsInfo<<"idx="<<idx<<"\n";
   gsInfo<<"idx_transfered="<<idx_transfered<<"\n";
   gsInfo<<"i="<<i<<"\n";
   gsInfo<<"it.row()="<<it.row()<<"\n";
   gsInfo<<"it.value()="<<it.value()<<"\n";
   gsInfo<<"newConstraint:"<<newConstraint.rows()<<"x"<<newConstraint.cols()<<"\n";*/
                                newConstraint[idx_transfered]=it.value();
                            }
                            Base::m_primalConstraints[k].push_back( give(newConstraint) );
                            Base::m_primalDofIndices[k].push_back( Base::m_primalDofIndices[patch][j] );
                            gsInfo << "Yup, it now is added!\n";
                        }
                        else
                            gsInfo << "No problem.\n";
                    }
                }
            }
        /*
                    // If there artificial dofs, we have to find all pre-images, which are then mapped back
                    std::vector< std::pair<index_t,index_t> > preImages;
                    m_dofMapperGlobal.preImage(dh.globalIndex, preImages);
                    gsInfo << "Found " << preImages.size() << " pre-images.\n";
                    for (size_t i=0; i<preImages.size(); ++i)
                    {
                        dof_helper dh2;
                        dh2.globalIndex = dh.globalIndex;
                        dh2.patch = preImages[i].first;
                        dh2.localIndex = m_dofMapperLocal2[preImages[i].first].index( preImages[i].second, 0 );
                        corners.push_back(give(dh2));
                    }
        */
    }

private:
    index_t                                       m_status;              ///< TODO: docs
    std::vector< std::vector<ArtificialIface> >   m_artificialIfaces;    ///< TODO: docs
    const gsMultiPatch<T>*                        m_multiPatch;          ///< Pointer to the respective multipatch
    const gsMultiBasis<T>*                        m_multiBasis;          ///< Pointer to the respective multibasis
    gsDofMapper                                   m_dofMapperOrig;       ///< The original dof mapper
    gsDofMapper                                   m_dofMapperMod;        ///< The modified dof mapper
    std::vector<gsDofMapper>                      m_dofMapperLocal2;     ///< The modified dof mapper, including the dofs for the neighboring patches
    gsMatrix<T>                                   m_fixedPart;           ///< TODO: docs
};

} // namespace gismo
