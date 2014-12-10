/** @file gsDofMapper.h

    @brief Provides the gsDofMapper class for re-indexing DoFs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#include <gsCore/gsDofMapper.h>


namespace gismo 
{

gsDofMapper::gsDofMapper() : 
m_shift(0), m_numFreeDofs(0), m_curElimId(-1), m_curCouplingId(1)
{ 
    m_offset.resize(1,0);
}


void gsDofMapper::localToGlobal(const gsMatrix<unsigned>& locals,
                                index_t patchIndex,
                                gsMatrix<unsigned>& globals) const
{
    GISMO_ASSERT( locals.cols() == 1, "localToGlobal: Expecting one column of locals");
    const index_t numActive = locals.rows();
    
    globals.resize(numActive,1);
    
    for (index_t i = 0; i < numActive; ++i)
        globals(i,0) = MAPPER_PATCH_DOF(locals(i,0), patchIndex)+m_shift;
}

// This function can not have enough information to do its job if dim>2
// GISMO_DEPRECATED
void gsDofMapper::matchInterface( index_t k1, index_t k2,
                                     const gsMatrix<unsigned> & b1,
                                     const gsMatrix<unsigned> & b2,
                                     const gsVector<bool> &orient )
{
    // Boundaries must be conforming (matching)
    const index_t sz = b1.size();
    if ( sz != b2.size() )
    {
        gsWarn<<"gsDofMapper: Problem: non-conforming interface "<<
                "("    <<k1<<","<< //i.first().side<<
                ")<->("<<k2<<","<< //i.second().side<<
                ") ~ ("<<sz<<","<<b2.size() <<").\n";

        gsWarn<< b1.transpose() << "\n";
        gsWarn<< b2.transpose() << "\n";

        GISMO_ERROR( "Non-conforming boundaries.\n" );
    }

    if( orient.size() )
    {
        if ( orient[0] == true ) // assumes 1D side
        {
            for ( index_t k=0; k<sz; ++k)
                this->matchDof( k1, b1(k,0), k2, b2(k,0) );
        }
        else
        {
            for ( index_t k=0; k<sz; ++k)
                this->matchDof( k1, b1(k,0), k2, b2(sz-k-1,0) );
        }
    }
    else
        for ( index_t k=0; k<sz; ++k)
            this->matchDof( k1, b1(k,0), k2, b2(k,0) );
}

void gsDofMapper::colapseDofs(index_t k, const gsMatrix<unsigned> & b )
{
    const index_t last = b.size()-1;
    for ( index_t k=0; k<last; ++k)
    {
        this->matchDof( k, b(k,0), k, b(k+1,0) );
    }
}





void gsDofMapper::matchDof( index_t u, index_t i, index_t v, index_t j )
{
    index_t d1 = MAPPER_PATCH_DOF(i,u);
    index_t d2 = MAPPER_PATCH_DOF(j,v);

    // make sure that d1 <= d2, simplifies implementation
    if (d1 > d2)
    {
        std::swap(d1, d2);
        std::swap(u, v);
        std::swap(i, j);
    }

    if (d1 < 0)         // first dof is eliminated
    {
        if (d2 < 0)         mergeDofsGlobally( d1, d2 );        // both are eliminated, merge their indices
        else if (d2 == 0)   MAPPER_PATCH_DOF(j,v) = d1;             // second is free, eliminate it along with first
        else /* d2 > 0*/    replaceDofGlobally( d2, d1 );       // second is coupling, eliminate all instances of it
    }
    else if (d1 == 0)   // first dof is a free dof
    {
        if (d2 == 0)        MAPPER_PATCH_DOF(i,u) = MAPPER_PATCH_DOF(j,v) = m_curCouplingId++;  // both are free, assign them a new coupling id
        else if (d2 > 0)    MAPPER_PATCH_DOF(i,u) = d2;             // second is coupling, add first to the same coupling group
        else                GISMO_ERROR("Something went terribly wrong");
    }
    else /* d1 > 0 */   // first dof is a coupling dof
    {
        GISMO_ASSERT(d2 > 0, "Something went terribly wrong");
        mergeDofsGlobally( d1, d2 );                            // both are coupling dofs, merge them
    }

    // if we merged two different non-eliminated dofs, we lost one free dof
    if ( (d1 != d2 && (d1 >= 0 || d2 >= 0) ) || (d1 == 0 && d2 == 0) )
        --m_numFreeDofs;
}


void gsDofMapper::markBoundary( index_t k, const gsMatrix<unsigned> & boundaryDofs )
{
    for (index_t i = 0; i < boundaryDofs.rows(); ++i)
    {
        eliminateDof( boundaryDofs(i,0), k );
    }

    // TO DO: save inverse, eg boundaryDofs(i,0)
}

void gsDofMapper::eliminateDof( index_t i, index_t k )
{
    const index_t old = MAPPER_PATCH_DOF(i,k);
    if (old == 0)       // regular free dof
    {
        --m_numFreeDofs;
        MAPPER_PATCH_DOF(i,k) = m_curElimId--;
    }
    else if (old > 0)   // coupling dof
    {
        --m_numFreeDofs;
        replaceDofGlobally( old, m_curElimId-- );
    }
    // else: old < 0: already an eliminated dof, nothing to do
}

void gsDofMapper::finalize()
{
    GISMO_ASSERT(m_curElimId!=0, "Error in gsDofMapper::finalize() called twince.");
    
    index_t curFreeDof = 0;                 // free dofs start at 0
    index_t curElimDof = m_numFreeDofs;     // eliminated dofs start after free dofs
    std::vector<index_t> couplingDofs(m_curCouplingId - 1, -1);
    std::vector<index_t> elimDofs(-m_curElimId - 1, -1);

    // Coupled dofs start after standard dofs
    index_t curCplDof = m_numFreeDofs - couplingDofs.size();

    for (std::size_t k = 0; k < m_dofs.size(); ++k)
    {
        const index_t dofType = m_dofs[k];

        if (dofType == 0)       // standard dof
            m_dofs[k] = curFreeDof++;
        else if (dofType < 0)   // eliminated dof
        {
            const index_t id = -dofType - 1;
            if (elimDofs[id] < 0)
                elimDofs[id] = curElimDof++;
            m_dofs[k] = elimDofs[id];
        }
        else // dofType > 0     // coupling dof
        {
            const index_t id = dofType - 1;
            if (couplingDofs[id] < 0)
                //couplingDofs[id] = curFreeDof++;
                couplingDofs[id] =   curCplDof++;
            m_dofs[k] = couplingDofs[id];
        }
    }
    m_numElimDofs = curElimDof - m_numFreeDofs;

    GISMO_ASSERT(curCplDof == m_numFreeDofs,
                 "gsDofMapper::finalize() - computed number of coupling dofs does not match allocated number");
    GISMO_ASSERT(curFreeDof + static_cast<index_t>(couplingDofs.size()) == m_numFreeDofs,
                 "gsDofMapper::finalize() - computed number of free dofs does not match allocated number");

    m_curElimId = 0;// Only equal to zero after finalize is called.

}

void gsDofMapper::print() const
{
    gsInfo<<"Dofs: "<< this->size() <<"\n";
    gsInfo<<" free: "<< this->freeSize() <<"\n";
    gsInfo<<" elim: "<< this->boundarySize() <<"\n";
}



void gsDofMapper::setIdentity(index_t nPatches, size_t nDofs)
{
    m_curElimId     = -1;
    m_curCouplingId =  1;
    m_numFreeDofs = nDofs;

    // Initialize all offsets to zero
    m_offset.resize(nPatches, 0);

    m_dofs.resize( m_numFreeDofs, 0);

    finalize();
}

void gsDofMapper::initPatchDofs(const gsVector<index_t> & patchDofSizes)
{
    m_curElimId     = -1;
    m_curCouplingId =  1;

    const size_t nPatches = patchDofSizes.size();

    // Initialize offsets and dof holder
    m_offset.reserve( nPatches );
    m_offset.push_back(0);
    for (size_t k = 1; k < nPatches; ++k)
    {
        m_offset.push_back( m_offset.back() + patchDofSizes[k-1] );
    }

    m_numFreeDofs = m_offset.back() + patchDofSizes[nPatches-1];

    m_dofs.resize( m_numFreeDofs, 0);
}



void gsDofMapper::replaceDofGlobally(index_t oldIdx, index_t newIdx)
{
    std::replace( m_dofs.begin(), m_dofs.end(), oldIdx, newIdx );
}

void gsDofMapper::mergeDofsGlobally(index_t dof1, index_t dof2)
{
    if (dof1 != dof2)
    {
        // replace the larger by the smaller for more consistent numbering.
        if (dof1 < dof2)
            std::swap(dof1, dof2);

        replaceDofGlobally(dof1, dof2);
    }
}


index_t gsDofMapper::coupledSize() const
{ 
    // Property: coupled (eliminated or not )DoFs appear more than once in the mapping.
    GISMO_ENSURE(m_curElimId==0, "finalize() was not called on gsDofMapper");

    index_t count = 0;
    std::map<index_t,index_t> CountMap;

    for (std::vector<index_t>::const_iterator it = m_dofs.begin(); it != m_dofs.end(); ++it)
        CountMap[*it]++;

    for (std::map<index_t,index_t>::const_iterator it = CountMap.begin(); it != CountMap.end(); ++it)
        if ( it->first > m_numFreeDofs && it->second > 1 )
            count++;

    return count; 
}


void gsDofMapper::setShift (index_t shift)
{
    m_shift=shift;
}

} // namespace gismo


