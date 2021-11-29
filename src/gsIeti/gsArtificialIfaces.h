/** @file gsArttficialIfaces.h

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
class gsArtificialIfaces {

public:
    gsArtificialIfaces( const gsMultiPatch<T>& mp, const gsMultiBasis<T>& mb, gsDofMapper dm )
        : m_status(0), m_artificialIfaces(mb.nBases()), m_multiPatch(&mp), m_multiBasis(&mb), m_dofMapperOrig(give(dm)) {}

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
        GISMO_ASSERT( m_status==0, "gsArtificialIfaces: Cannot register artificial interfaces after requesting data." );
        ArtificialIface ai;
        ai.realIface = realIface;
        ai.artificialIface = artificialIface;
        ai.ifaceIndices = (*m_multiBasis)[artificialIface.patch].boundary(artificialIface);
        m_artificialIfaces[realIface.patch].push_back(give(ai));
    }

    /// Calls \ref registerArtificialIface for all interface known to the underlying box topology
    void registerAllArtificialIfaces()
    {
        GISMO_ASSERT( m_status==0, "gsArtificialIfaces: Cannot register artificial interfaces after requesting data." );
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
        GISMO_ASSERT( m_status==0, "gsArtificialIfaces: Already finalized." );
        m_status = 1;

        // TODO: we can get rid of this assumption if we incorporate dofs that
        // have been coupled into m_dofMapperOrig
        GISMO_ASSERT( m_dofMapperOrig.coupledSize() == 0, "Not implemented." );

        const index_t nPatches = m_dofMapperOrig.numPatches();

        m_dofMapperLocal.clear();
        m_dofMapperLocal.reserve(nPatches);

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
            m_dofMapperLocal.push_back( gsDofMapper(sizes[k]) );
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
                    m_dofMapperLocal[k].eliminateDof(i,0);
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
                            m_dofMapperLocal[k].eliminateDof(ii,l+1);
                        }
                        ++idx;
                    }
                }
            }
            // Eliminate all dofs on m_dofMapperLocal that are not on the artificial interface
            for (index_t l=0; l<nArtIf; ++l)
            {
                index_t last_idx = -1;
                for (index_t i=0; i<m_artificialIfaces[k][l].ifaceIndices.rows(); ++i)
                {
                    const index_t idx = m_artificialIfaces[k][l].ifaceIndices[i];
                    GISMO_ASSERT (idx>last_idx, "Not sorted; k="<<k<<"; l="<<l);
                    for (index_t j=last_idx+1; j<idx; ++j)
                        m_dofMapperLocal[k].eliminateDof(j,l+1);
                    last_idx = idx;
                }
                GISMO_ASSERT (sizes[k][l+1]>last_idx, "Not sorted; k="<<k<<"; l="<<l);
                for (index_t j=last_idx+1; j<sizes[k][l+1]; ++j)
                    m_dofMapperLocal[k].eliminateDof(j,l+1);
            }
            m_dofMapperLocal[k].finalize();
        }
        m_dofMapperMod.finalize();
    }

    gsSparseMatrix<index_t> subdomainToPatch()
    {
        if (!m_status) finalize();
        const index_t nPatches = m_artificialIfaces.size();
        gsSparseMatrix<index_t> result(nPatches,nPatches);
        result.setIdentity();
        for (index_t k=0; k<nPatches; ++k)
        {
            index_t idx = 1;
            const index_t nArtIf = m_artificialIfaces[k].size();
            for (index_t i=0; i<nArtIf; ++i)
            {
                index_t rr = m_artificialIfaces[k][i].realIface.patch;
                index_t aa = m_artificialIfaces[k][i].artificialIface.patch;
                result(rr,aa) = idx;
                ++idx;
            }
        }

        return result;
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
        return m_dofMapperLocal;
    }

    const gsDofMapper& dofMapperLocal(index_t k)
    {
        if (!m_status) finalize();
        return m_dofMapperLocal[k];
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

private:
    index_t                                       m_status;
    std::vector< std::vector<ArtificialIface> >   m_artificialIfaces;    ///< TODO: docs
    const gsMultiPatch<T>*                        m_multiPatch;          ///< Pointer to the respective multipatch
    const gsMultiBasis<T>*                        m_multiBasis;          ///< Pointer to the respective multibasis
    gsDofMapper                                   m_dofMapperOrig;       ///< The original dof mapper
    gsDofMapper                                   m_dofMapperMod;        ///< The modified dof mapper
    std::vector<gsDofMapper>                      m_dofMapperLocal;      ///< The modified dof mapper, including the dofs for the neighboring patches

};

} // namespace gismo
