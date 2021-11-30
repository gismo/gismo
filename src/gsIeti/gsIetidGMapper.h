/** @file gsIetidGMapper.h

    @brief Algorithms that help with assembling the matrices required for IETI-solvers
    in the context of discontinuous Galerkin discretizations

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, R. Schneckenleitner
*/

#include <gsIeti/gsIetiMapper.h>

namespace gismo {

// unused:
template <typename Container, typename Element>
static index_t indexOf( const Container& c, Element e ) //TODO: unused
{
    const Element* it = std::lower_bound(c.begin(),c.end(),e);
    if (it==c.end() || *it>e) return -1;
    return it-c.begin();
}

/** @brief Ieti Mapper for discontinuous Galerkin discretizations
 *
 *  Algorithms that help with assembling the matrices required for IETI-Solvers
 *  in the context of discontinuous Galerkin discretizations
 *
 *  See also: \a gsIetiMapper
 *
 *  @ingroup Solver
*/
template <class T>
class gsIetidGMapper : protected gsIetiMapper<T> {
    // Use protected inherence in order not to need virtual functions
    typedef gsIetiMapper<T> Base;

public:
    /// Data structure for artificial interfaces
    struct ArtificialIface {
       /// The artificial interface (foreign patch from which the artificial interface taken+side)
       patchSide           takenFrom;
       /// The real interface (local patch to which the artificial interface is assigned+side)
       patchSide           assignedTo;
    };

    /// @brief Default constructor
    gsIetidGMapper() {}

    /// @brief Create the ieti mapper
    ///
    /// @param multiBasis       The corresponding \a gsMultiBasis
    /// @param dofMapperGlobal  The dof mapper for the global problem
    /// @param fixedPart        The function values for the eliminated dofs
    /// @param artIfaces        A vector of all artificial interfaces to be incorporated
    ///
    /// The multibasis object has to outlive the ieti mapper
    gsIetidGMapper(
        const gsMultiBasis<T>& multiBasis,
        gsDofMapper dofMapperGlobal,
        const gsMatrix<T>& fixedPart,
        const std::vector<ArtificialIface>& artIfaces
    )
    { init(multiBasis, dofMapperGlobal, fixedPart, artIfaces); }

    /// @brief Init the ieti mapper after default construction
    ///
    /// @param multiBasis       The corresponding \a gsMultiBasis
    /// @param dofMapperGlobal  The dof mapper for the global problem
    /// @param fixedPart        The function values for the eliminated dofs
    /// @param artIfaces        A vector of all artificial interfaces to be incorporated
    ///
    /// The multibasis object has to outlive the ieti mapper
    void init(
        const gsMultiBasis<T>& multiBasis,
        gsDofMapper dofMapperGlobal,
        const gsMatrix<T>& fixedPart,
        const std::vector<ArtificialIface>& artIfaces
    );

    /// @brief Gives a vector that contains all possible artificial interfaces for the given discretization
    ///
    /// @param mb     The corresponding multi basis
    static std::vector<ArtificialIface> allArtificialIfaces( const gsMultiBasis<T>& mb );

    /// @brief Returns artificial interfaces for given patch
    ///
    /// @param k      The coresponding patch number
    const std::vector<ArtificialIface>& artificialIfaces(index_t k) const
    {
        GISMO_ASSERT( m_status, "gsIetidGMapper: After calling the default constructur, data must be given via member init." );
        return m_artificialIfaces[k];
    }

    /// @brief Returns local dof mapper for given patch
    ///
    /// @param k      The coresponding patch number
    const gsDofMapper& dofMapperLocal(index_t k) const
    {
        GISMO_ASSERT( m_status, "gsIetidGMapper: After calling the default constructur, data must be given via member init." );
        return m_dofMapperLocal2[k];
    }

    void localSpaces(const gsMultiPatch<T>& mp, index_t patch, gsMultiPatch<T>& mp_local, gsMultiBasis<T>& mb_local) const
    {
        GISMO_ASSERT( m_status, "gsIetidGMapper: After calling the default constructur, data must be given via member init." );

        const std::vector<ArtificialIface>& artIfaces = m_artificialIfaces[patch];
        std::vector<gsGeometry<T>*> mp_local_data;
        std::vector<gsBasis<T>*> mb_local_data;
        mp_local_data.reserve(1+2*mp.geoDim());
        mb_local_data.reserve(1+2*mp.geoDim());
        mp_local_data.push_back(mp[patch].clone().release());
        mb_local_data.push_back((*m_multiBasis)[patch].clone().release());
        for (size_t i=0; i<artIfaces.size(); ++i)
    {
             mp_local_data.push_back(mp[artIfaces[i].takenFrom.patch].clone().release());
             mb_local_data.push_back((*m_multiBasis)[artIfaces[i].takenFrom.patch].clone().release());
        }
        mp_local = gsMultiPatch<T>(mp_local_data);
        mp_local.computeTopology();
        mb_local = gsMultiBasis<T>(mb_local_data,mp_local);
    }

    gsMatrix<T> fixedPart(index_t k) const
    {
        GISMO_ASSERT( Base::m_status, "gsIetidGMapper: After calling the default constructur, data must be given via member init." );
        gsMatrix<T> result(dofMapperLocal(k).boundarySize(),1);
        result.setZero();
        result.topRows(Base::fixedPart(k).rows()) = Base::fixedPart(k); // TODO: this is not correct for the artificial interfaces
        return result;
    }

    void cornersAsPrimals( bool includeIsolated=false )
    {
        const index_t primals_old = nPrimalDofs();
        Base::cornersAsPrimals(true);
        const index_t primals_new = nPrimalDofs();
        transferConstraintsToArtificialIfaces(primals_old,primals_new);
        if (!includeIsolated) removeIsolatedPrimalConstraints(primals_old,primals_new);
    }

    void interfaceAveragesAsPrimals( const gsMultiPatch<T>& geo, short_t d, bool includeIsolated=false )
    {
        const index_t primals_old = nPrimalDofs();
        Base::interfaceAveragesAsPrimals(geo,d,true);
        const index_t primals_new = nPrimalDofs();
        transferConstraintsToArtificialIfaces(primals_old,primals_new);
        if (!includeIsolated) removeIsolatedPrimalConstraints(primals_old,primals_new);
    }

public:
    // Make members of base clase accessible
    // Functions that are overwritten here have been commented out.

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
    using Base::multiBasis;

private:
    // Internal functions
    void transferConstraintsToArtificialIfaces(index_t primals_old, index_t primals_new );

    void removeIsolatedPrimalConstraints(index_t primals_old, index_t primals_new);

private:
    // Data members and accessors to them
    using Base::m_status;
    using Base::m_multiBasis;

    std::vector< std::vector<ArtificialIface> >   m_artificialIfaces;    ///< TODO: docs
    /// The indicies of the basis of the foreign patch that belong to artificial interface
    /// Vector is assumed to be sorted. Same order as in m_artificialIfaces
    std::vector< std::vector< gsVector<index_t> > > m_ifaceIndices;      ///< TODO: docs
    gsDofMapper                                   m_dofMapperOrig;       ///< The original dof mapper
    gsDofMapper                                   m_dofMapperMod;        ///< The modified dof mapper
    std::vector<gsDofMapper>                      m_dofMapperLocal2;     ///< The modified dof mapper, including the dofs for the neighboring patches
};

} // namespace gismo

// To be moved to hpp file later //////////////

namespace gismo {

template <class T>
std::vector<typename gsIetidGMapper<T>::ArtificialIface> gsIetidGMapper<T>::allArtificialIfaces( const gsMultiBasis<T>& mb )
{
    std::vector<ArtificialIface> result;
    const gsBoxTopology& top = mb.topology();
    const short_t dim        = mb.domainDim();
    const index_t nPatches   = mb.nPieces();
    result.reserve(nPatches*dim*4);
    for (index_t k=0; k<nPatches; ++k)
    {
        for (boxSide it = boxSide::getFirst(dim); it!=boxSide::getEnd(dim); ++it)
        {
            ArtificialIface artIface;
            artIface.assignedTo = patchSide(k,it);
            if(top.getNeighbour(artIface.assignedTo, artIface.takenFrom))
                result.push_back(artIface);
        }
    }
    return result;
}

template <class T>
void gsIetidGMapper<T>::init( const gsMultiBasis<T>& mb, gsDofMapper dm, const gsMatrix<T>& fixedPart, const std::vector<ArtificialIface>& artIfaces )
{
    GISMO_ASSERT( !m_status,
        "gsIetidGMapper::init cannot be called after object was initialized by previous call or value constructor." );

    m_artificialIfaces.resize(mb.nBases());
    m_ifaceIndices.resize(mb.nBases());
    m_dofMapperOrig = dm;

    // We store them ordered
    index_t sz = artIfaces.size();
    for (index_t i=0; i<sz; ++i)
    {
        m_artificialIfaces[ artIfaces[i].assignedTo.patch ].push_back( artIfaces[i] );
        m_ifaceIndices    [ artIfaces[i].assignedTo.patch ].push_back(
            mb[artIfaces[i].takenFrom.patch].boundary(artIfaces[i].takenFrom)
        );
    }

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
            sizes[k][l+1]      = m_dofMapperOrig.patchSize( m_artificialIfaces[k][l].takenFrom.patch );
            subdomainSizes[k] += m_ifaceIndices[k][l].rows();
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
                for (index_t i=0; i<m_ifaceIndices[k][l].rows(); ++i)
                {
                    const index_t kk = m_artificialIfaces[k][l].takenFrom.patch;
                    const index_t ii = m_ifaceIndices[k][l][i];
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
            for (index_t i=0; i<m_ifaceIndices[k][l].rows(); ++i)
            {
                const index_t idx = m_ifaceIndices[k][l][i];
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

    Base::init( mb, m_dofMapperMod, fixedPart );

}


template <class T>
void gsIetidGMapper<T>::transferConstraintsToArtificialIfaces(index_t primals_old, index_t primals_new )
{
    for (size_t k=0; k<m_artificialIfaces.size(); ++k)
    {
        for (size_t i=0; i<m_artificialIfaces[k].size(); ++i)
        {
            ArtificialIface& ai = m_artificialIfaces[k][i];
            const index_t patch = ai.takenFrom.patch;
            for (size_t j=0; j<Base::m_primalDofIndices[patch].size(); ++j)
            {
                if (primals_old<=Base::m_primalDofIndices[patch][j]
                    && Base::m_primalDofIndices[patch][j]<primals_new)
                {
                    //gsInfo <<"Constraint #"<<j<<" of patch "<<patch<<" (global primal dof # "
                    //        <<Base::m_primalDofIndices[patch][j]<<") might needed to be added to patch # "
                    //        <<k<<"(=="<<ai.assignedTo.patch<<")\n";

                    std::vector<index_t> ifaceIndicesMapped(m_ifaceIndices[k][i].size());
                    for (index_t ii=0; ii<m_ifaceIndices[k][i].size(); ++ii)
                    {
                        ifaceIndicesMapped[ii] = m_dofMapperLocal2[patch].index(m_ifaceIndices[k][i][ii],0/**/);
                    }

                    bool b = true;
                    for (typename gsSparseVector<T>::InnerIterator it(Base::m_primalConstraints[patch][j],0);
                        it; ++it)
                    {
                        b &= (std::find(ifaceIndicesMapped.begin(),ifaceIndicesMapped.end(),it.row()) != ifaceIndicesMapped.end());
                    }
                    if (b)
                    {
                        //gsInfo << "Yes it, should be added!\n";
                        gsSparseVector<T> newConstraint(m_dofMapperLocal2[k].freeSize());
                        for (typename gsSparseVector<T>::InnerIterator it(Base::m_primalConstraints[patch][j],0);
                            it; ++it)
                        {
                            index_t idx;
                            {
                                std::vector<std::pair<index_t,index_t> > result;
                                m_dofMapperLocal2[patch].preImage(it.row(),result);
                                GISMO_ASSERT(result.size()==1, "result.size()="<<result.size());
                                GISMO_ASSERT(result[0].first == 0, "result[0].first="<<result[0].first);
                                idx = result[0].second;
                            }

                            index_t idx_transfered = m_dofMapperLocal2[k].index(idx,i+1);
                            newConstraint[idx_transfered]=it.value();
                        }
                        Base::m_primalConstraints[k].push_back( give(newConstraint) );
                        Base::m_primalDofIndices[k].push_back( Base::m_primalDofIndices[patch][j] );
                        //gsInfo << "Yup, it now is added!\n";
                    }
                    //else
                        //gsInfo << "No problem.\n";
                }
            }
        }
    }
}

template <class T>
void gsIetidGMapper<T>::removeIsolatedPrimalConstraints(index_t primals_old, index_t primals_new)
{
    std::vector<index_t> count(primals_new);
    std::vector< std::vector<index_t> >& primalDofIndices = Base::m_primalDofIndices;
    for (size_t k=0; k<primalDofIndices.size(); ++k)
        for (size_t l=0; l<primalDofIndices[k].size(); ++l)
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

    for (size_t k=0; k<primalDofIndices.size(); ++k)
        for (size_t l=0; l<primalDofIndices[k].size(); ++l)
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


} // namespace gismo
