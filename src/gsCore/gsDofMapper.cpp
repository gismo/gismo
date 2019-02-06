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
m_shift(0), m_numFreeDofs(0), m_numCpldDofs(1), m_curElimId(-1)
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

    /* //Testing second overload of localToGlobal
    index_t nf;
    gsMatrix<unsigned> tmp;
    localToGlobal(locals, patchIndex, tmp, nf);
    for (index_t i = 0; i != numActive; ++i)
        globals.at(tmp(i,0)) = tmp(i,1);
    return;
    */

    for (index_t i = 0; i < numActive; ++i)
        globals(i,0) = index(locals(i,0), patchIndex);
}

void gsDofMapper::localToGlobal(const gsMatrix<unsigned>& locals,
                                index_t patchIndex,
                                gsMatrix<unsigned>& globals,
                                index_t & numFree) const
{
    GISMO_ASSERT( locals.cols() == 1, "localToGlobal: Expecting one column of locals");
    GISMO_ASSERT( &locals != &globals, "localToGlobal: Inplace not supported");
    const index_t numActive = locals.rows();
    globals.resize(numActive, 2);

    numFree = 0;
    index_t bot = numActive;
    for (index_t i = 0; i != numActive; ++i)
    {
        const index_t ii = index(locals(i,0), patchIndex);
        if ( is_free_index(ii) )
        {
            globals(numFree  , 0) = i ;
            globals(numFree++, 1) = ii;
        }
        else // is_boundary_index(ii)
        {
            globals(--bot, 0) = i ;
            globals(  bot, 1) = ii;
        }
    }

    //GISMO_ASSERT(numFree == bot, "Something went wrong in localToGlobal");
}

gsVector<index_t> gsDofMapper::asVector() const
{
    gsVector<index_t> v(m_dofs.size());
    for(index_t i = 0; i!= v.size(); ++i)
        v[i] = m_dofs[i]+m_shift;
    return v;
}

gsVector<index_t> gsDofMapper::inverseAsVector() const
{
    GISMO_ASSERT(isPermutation(), "This dofMapper is not 1-1");
    gsVector<index_t> v(size());
    for(index_t i = 0; i!= v.size(); ++i)
        v[m_dofs[i]] = i;
    return v;
}

void gsDofMapper::colapseDofs(index_t k, const gsMatrix<unsigned> & b )
{
    const index_t last = b.size()-1;

    for ( index_t l=0; l!=last; ++l)
    {
        this->matchDof( k, b(l,0), k, b(l+1,0) );
    }
}

void gsDofMapper::matchDof( index_t u, index_t i, index_t v, index_t j )
{
    //GISMO_ASSERT(u < m_offset.size() && i< m_offset[u+1] - m_offset[u] ,"Invalid patch-local dof "<<i<<" in matchDof.");
    //GISMO_ASSERT(v < m_offset.size() && j< m_offset[v+1] - m_offset[v] ,"Invalid patch-local dof "<<i<<" in matchDof.");
    //TODO: add check for correct range [by adding a last offset = number of total dofs]

    GISMO_ASSERT(static_cast<size_t>(u)<numPatches(), "Invalid patch index "<< u <<" >= "<< numPatches() );
    GISMO_ASSERT(static_cast<size_t>(v)<numPatches(), "Invalid patch index "<< v <<" >= "<< numPatches() );

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
        if (d2 < 0)
            mergeDofsGlobally( d1, d2 );  // both are eliminated, merge their indices
        else if (d2 == 0)
            MAPPER_PATCH_DOF(j,v) = d1;   // second is free, eliminate it along with first
        else /* d2 > 0*/
            replaceDofGlobally( d2, d1 ); // second is coupling, eliminate all instances of it
    }
    else if (d1 == 0)   // first dof is a free dof
    {
        if (d2 == 0)
        {
            MAPPER_PATCH_DOF(i,u) = MAPPER_PATCH_DOF(j,v) = m_numCpldDofs++;  // both are free, assign them a new coupling id
            if (u==v && i==j) return;
        }
        else if (d2 > 0)
            MAPPER_PATCH_DOF(i,u) = d2;   // second is coupling, add first to the same coupling group
        else
            GISMO_ERROR("Something went terribly wrong");
    }
    else /* d1 > 0 */   // first dof is a coupling dof
    {
        GISMO_ASSERT(d2 > 0, "Something went terribly wrong");
        mergeDofsGlobally( d1, d2 );      // both are coupling dofs, merge them
    }

    // if we merged two different non-eliminated dofs, we lost one free dof
    if ( (d1 != d2 && (d1 >= 0 || d2 >= 0) ) || (d1 == 0 && d2 == 0) )
        --m_numFreeDofs;
}

void gsDofMapper::matchDofs(index_t u, const gsMatrix<unsigned> & b1,
                            index_t v,const gsMatrix<unsigned> & b2)
{
    const index_t sz = b1.size();
    GISMO_ASSERT( sz == b2.size(), "Waiting for same number of DoFs");
    for ( index_t k=0; k<sz; ++k)
        this->matchDof( u, b1(k,0), v, b2(k,0) );
}

void gsDofMapper::markCoupled( index_t i, index_t k )
{
    matchDof(k,i,k,i);
}

void gsDofMapper::markTagged( index_t i, index_t k )
{
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );

    //see gsSortedVector::push_sorted_unique
    index_t t = index(i,k);
    std::vector<index_t>::iterator pos = std::lower_bound(m_tagged.begin(), m_tagged.end(), t );

    if ( pos == m_tagged.end() || *pos != t )// If not found
        m_tagged.insert(pos, t);
}

void gsDofMapper::markBoundary( index_t k, const gsMatrix<unsigned> & boundaryDofs )
{
    for (index_t i = 0; i < boundaryDofs.rows(); ++i)
    {
        eliminateDof( boundaryDofs(i,0), k );
    }

    // TO DO: save inverse, eg boundaryDofs(i,0)
}

void gsDofMapper::markCoupledAsTagged()
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    m_tagged.reserve(m_tagged.size()+m_numCpldDofs);
    for(int i=0; i< m_numCpldDofs;++i)
        m_tagged.push_back(m_numFreeDofs-m_numCpldDofs+i);

    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique(m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );
}

void gsDofMapper::eliminateDof( index_t i, index_t k )
{
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
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
    GISMO_ASSERT(m_curElimId!=0, "Error in gsDofMapper::finalize() called twice.");

    // For assigning coupling and eliminated dofs to continuous
    // indices (-1 = unassigned)
    std::vector<index_t> couplingDofs(m_numCpldDofs -1, -1);
    std::vector<index_t> elimDofs    (-m_curElimId - 1, -1);

    // Free dofs start at 0
    index_t curFreeDof = 0;
    // Coupling dofs start after standard dofs (=num of zeros in m_dofs)
    index_t curCplDof = std::count(m_dofs.begin(), m_dofs.end(), 0);

    // Eliminated dofs start after free dofs
    index_t curElimDof = m_numFreeDofs;
    // Devise number of coupled dofs (m_numCpldDofs was used as
    // coupling id up to here)
    m_numCpldDofs = m_numFreeDofs - curCplDof;

    /*// For debugging: counting the number of coupled and boundary dofs
    std::vector<index_t> alldofs = m_dofs;
    std::sort( alldofs.begin(), alldofs.end() );
    alldofs.erase( std::unique( alldofs.begin(), alldofs.end() ), alldofs.end() );
    const index_t numCoupled =
    std::count_if( alldofs.begin(), alldofs.end(),
                          std::bind2nd(std::greater<index_t>(), 0) );
    const index_t numBoundary =
    std::count_if( alldofs.begin(), alldofs.end(),
                          std::bind2nd(std::less<index_t>(), 0) );
    */

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
                couplingDofs[id] = curCplDof++;
            m_dofs[k] = couplingDofs[id];
        }
    }

    GISMO_ASSERT(curCplDof == m_numFreeDofs,
                 "gsDofMapper::finalize() - computed number of coupling "
                 "dofs does not match allocated number, "<<curCplDof<<"!="<<m_numFreeDofs);
    GISMO_ASSERT(curFreeDof + m_numCpldDofs == m_numFreeDofs,
                 "gsDofMapper::finalize() - computed number of free dofs "
                 "does not match allocated number");

    // Devise number of eliminated dofs
    m_numElimDofs = curElimDof - m_numFreeDofs;

    m_curElimId = 0;// Only equal to zero after finalize is called.
}

std::ostream& gsDofMapper::print( std::ostream& os ) const
{
    os<<" Dofs: "<< this->size() <<"\n";
    os<<" free: "<< this->freeSize() <<"\n";
    os<<" coupled: "<< this->coupledSize() <<"\n";
    os<<" tagged: "<< this->taggedSize() <<"\n";
    os<<" elim: "<< this->boundarySize() <<"\n";
    return os;
}

void gsDofMapper::setIdentity(index_t nPatches, size_t nDofs)
{
    m_curElimId   = -1;
    m_numFreeDofs = nDofs;
    m_numElimDofs = 0;
    m_numCpldDofs =  1;

    m_shift = m_bshift = 0;

    // Initialize all offsets to zero
    m_offset.resize(nPatches, 0);

    m_dofs.resize( m_numFreeDofs, 0);
}

void gsDofMapper::permuteFreeDofs(const gsVector<index_t>& permutation)
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(m_numFreeDofs == permutation.size(), "permutation size does not match number of free dofs");
    //GISMO_ASSERT(m_tagged.empty(), "you cannot permute the dofVector twice, combine the permutation");
    //We could check here, that permutation is indeed a permutation

    //make a copy of the old ordering, easiest way to implement the permutation. Inplace reordering is quite hard.
    std::vector<index_t> dofs = m_dofs;
    //first use a tempory coupled vector, other wise is_coupled_index() would not work correctly
    std::vector<index_t> tagged_permuted;
    tagged_permuted.reserve(taggedSize()); //reserve enough memory

    for(index_t i=0; i<(index_t)dofs.size();++i)
    {
        index_t idx = dofs[i];
        if(is_free_index(idx))
        {
            m_dofs[i] = permutation[idx];
            //fill  bookkeeping for tagged dofs
            if(is_tagged_index(idx))
                tagged_permuted.push_back(m_dofs[i]);
        }
        else if(is_tagged_index(idx)) //Take care about eliminated tagged dofs
            tagged_permuted.push_back(idx);
    }
    m_tagged.swap(tagged_permuted);

    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique (m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );

    //coupled dofs cannot be tracked anymore
    m_numCpldDofs = 0;
}


void gsDofMapper::initPatchDofs(const gsVector<index_t> & patchDofSizes)
{
    m_curElimId   = -1;
    m_numElimDofs = 0;
    m_numCpldDofs =  1;
    m_shift = m_bshift = 0;

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

void gsDofMapper::preImage(const index_t gl,
                           std::vector<std::pair<index_t,index_t> > & result) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    typedef std::vector<index_t>::const_iterator citer;

    result.clear();

    std::size_t cur = 0;//local offsetted index
    for (citer it = m_dofs.begin(); it != m_dofs.end(); ++it, ++cur)
    {
        if ( *it == gl )
        {
            // Get the patch index of "cur" by "un-offsetting"
            const index_t patch = std::upper_bound(m_offset.begin(), m_offset.end(), cur)
                                - m_offset.begin() - 1;

            // Found a patch-dof pair
            result.push_back( std::make_pair(patch, cur - m_offset[patch] - m_shift) );
        }
    }
}

bool gsDofMapper::indexOnPatch(const index_t gl, const index_t k) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    const citer istart = m_dofs.begin()+m_offset[k];
    const citer iend   = istart + patchSize(k);
    return (std::find(istart, iend, gl)!=iend);
}

index_t gsDofMapper::boundarySizeWithDuplicates() const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");

    const index_t s = m_numFreeDofs + m_shift - 1;
    return std::count_if(m_dofs.begin(), m_dofs.end(),
                         std::bind2nd(std::greater<index_t>(), s) );
}


index_t gsDofMapper::coupledSize() const
{
    GISMO_ENSURE(m_curElimId==0, "finalize() was not called on gsDofMapper");
    return m_numCpldDofs;
/*// Implementation without saving this number:
    // Property: coupled (eliminated or not) DoFs appear more than once in the mapping.
    GISMO_ENSURE(m_curElimId==0, "finalize() was not called on gsDofMapper");

    std::vector<index_t> CountMap(m_numFreeDofs,0);

    // Count number of appearances of each free DoF
    for (std::vector<index_t>::const_iterator it = m_dofs.begin(); it != m_dofs.end(); ++it)
        if ( *it < m_numFreeDofs )
            CountMap[*it]++;

    // Count the number of freeDoFs that appear more than once
    return std::count_if( CountMap.begin(), CountMap.end(),
                          std::bind2nd(std::greater<index_t>(), 1) );
*/
}

index_t gsDofMapper::taggedSize() const
{
    return m_tagged.size();
}

template<class Predicate, class Iterator>
gsVector<index_t> gsDofMapper::find_impl(Iterator istart, Iterator iend, Predicate pred)
{
    gsVector<index_t> rvo(std::count_if(istart, iend, pred));
    index_t * c = rvo.data();
    Iterator cur = std::find_if(istart, iend, pred);
    while( cur!=iend )
    {
        *(c++) = std::distance(istart,cur);
        cur = std::find_if(cur+1, iend, pred);
    }
    return rvo;
}

gsVector<index_t> gsDofMapper::findBoundary(const index_t k) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    const index_t s = m_numFreeDofs + m_shift - 1;
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs.begin() + m_offset[k];
    citer iend   = istart + patchSize(k);
    return find_impl(istart, iend, std::bind2nd(std::greater<index_t>(),s));
}

gsVector<index_t> gsDofMapper::findFree(const index_t k) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    const index_t s = m_numFreeDofs + m_shift;
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs.begin() + m_offset[k];
    citer iend   = istart + patchSize(k);
    return find_impl(istart, iend, std::bind2nd(std::less<index_t>(),s));
}

struct _isBetween
{
    _isBetween(const index_t l, const index_t u) : _l(l), _u(u) { }
    index_t _l, _u;
    bool operator()(const index_t i) { return  (i < _u) && (i > _l); }
};

gsVector<index_t> gsDofMapper::findCoupled(const index_t k, const index_t j) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );

    if (k==j) return gsVector<index_t>();

    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs.begin() + m_offset[k];
    citer iend   = istart + patchSize(k);

    const index_t l = m_numFreeDofs+m_shift-m_numCpldDofs-1;
    const index_t u = m_numFreeDofs+m_shift;
    if (-1==j)
        return find_impl(istart, iend, _isBetween(l,u) );
    else
    {
        citer istartj = m_dofs.begin() + m_offset[j];
        citer iendj   = istartj + patchSize(j);
        std::list<index_t> v;
        citer cur = std::find_if(istart, iend, _isBetween(l,u));
        while( cur!=iend )
        {
            if ( std::find(istartj,iendj,*cur)!=iendj )
                v.push_back( std::distance(istart,cur) );
            cur = std::find_if(cur+1, iend, _isBetween(l,u));
        }

        gsVector<index_t> res;
        res.resize(v.size());
        index_t * a = res.data();
        for( std::list<index_t>::const_iterator it = v.begin();
             it!=v.end(); ++it) *(a++) = *it;
        return res;
    }
}

gsVector<index_t> gsDofMapper::findFreeUncoupled(const index_t k) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    const citer istart = m_dofs.begin() + m_offset[k];
    const citer iend   = istart + patchSize(k);
    return find_impl(istart, iend,
                     _isBetween(m_shift-1, m_numFreeDofs+m_shift-m_numCpldDofs) );
}

gsVector<index_t> gsDofMapper::findTagged(const index_t k) const
{
    GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs.begin() + m_offset[k];
    citer iend   = istart + patchSize(k);

    std::vector<index_t> si;
//    std::set_intersection<index_t>(istart, iend,m_tagged.begin(), m_tagged.end(),
//                                   std::back_inserter(si));
    gsVector<index_t> rvo;
    //rvo.swap(si);
    return rvo;
}

void gsDofMapper::setShift (index_t shift)
{
    m_shift=shift;
}

void gsDofMapper::setBoundaryShift (index_t shift)
{
    m_bshift=shift;
}

} // namespace gismo
