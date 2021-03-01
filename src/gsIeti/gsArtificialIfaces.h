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
    gsArtificialIfaces( const gsMultiBasis<T>& mb, gsDofMapper dm, gsMatrix<T> fd )
        : m_artificialIfaces(mb.nBases()), m_multiBasis(&mb), m_dofMapperOrig(give(dm)), m_fixedDofsOrig(give(fd)) {}

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
        ArtificialIface ai;
        ai.realIface = realIface;
        ai.artificialIface = artificialIface;
        ai.ifaceIndices = (*m_multiBasis)[realIface.patch].boundary(artificialIface);
        m_artificialIfaces[realIface.patch].push_back(give(ai));
    }

    /// Calls \ref registerArtificialIface for all interface known to the underlying box topology
    void registerAllArtificialIfaces()
    {
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
    static index_t indexOf( const Container& c, Element e )
    {
        const Element* it = std::lower_bound(c.begin(),c.end(),e);
        if (it==c.end() || *it>e) return -1;
        return it-c.begin();
    }

    /// Setup of the dof mappers
    void finalize()
    {
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
                            m_dofMapperMod.eliminateDof(idx,k); // TODO: necessary?
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
                    for (index_t j=last_idx+1; j<idx; ++j)
                        m_dofMapperLocal[k].eliminateDof(j,l+1);
                    last_idx = idx;
                }
                for (index_t j=last_idx+1; j<sizes[k][l+1]; ++j)
                    m_dofMapperLocal[k].eliminateDof(j,l+1);
            }
            m_dofMapperLocal[k].finalize();
        }

        m_dofMapperMod.finalize();

        m_fixedDofsModLocal.resize(nPatches);
        for (index_t k=0; k<nPatches; ++k)
            m_fixedDofsModLocal[k].setZero(m_dofMapperLocal[k].boundarySize(),1);
        // TODO: compute m_fixedDofsModLocal (currently this is only correct for hom. Dir. bc).

    }


public:

    /// @brief Returns artificial interfaces for given patch
    const std::vector<ArtificialIface>& artificialIfaces(index_t k) const  { return m_artificialIfaces[k];        }

    /// @brief Reference to the multi basis object being passed to constructur or \ref init
    const gsMultiBasis<T>& multiBasis() const                              { return *m_multiBasis;                }

    /// @brief The original dof mapper
    const gsDofMapper& dofMapperOrig() const                               { return m_dofMapperOrig;              }

    /// @brief The modified dof mapper
    const gsDofMapper& dofMapperMod() const                                { return m_dofMapperMod;               }

    /// @brief The local dof mapper
    const gsDofMapper& dofMapperLocal(index_t patch) const                 { return m_dofMapperLocal[patch];      }

    /// @brief The original dof mapper
    const gsMatrix<T>& fixedDofsOrig() const                               { return m_fixedDofsOrig;              }

    /// @brief The modified dof mapper
    const gsMatrix<T>& fixedDofs(index_t patch) const                      { return m_fixedDofsModLocal[patch];   }

private:
    std::vector< std::vector<ArtificialIface> >   m_artificialIfaces;    ///< TODO: docs
    const gsMultiBasis<T>*                        m_multiBasis;          ///< Pointer to the respective multibasis
    gsDofMapper                                   m_dofMapperOrig;       ///< The original dof mapper
    gsDofMapper                                   m_dofMapperMod;        ///< The modified dof mapper, including the dofs for the art ifaces
    std::vector<gsDofMapper>                      m_dofMapperLocal;      ///< The modified dof mapper, including the dofs for the neighboring patches
    gsMatrix<T>                                   m_fixedDofsOrig;       ///< TODO: docs
    std::vector< gsMatrix<T> >                    m_fixedDofsModLocal;   ///< TODO: docs


};

} // namespace gismo
