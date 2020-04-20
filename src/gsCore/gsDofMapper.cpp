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
  m_offset(1,0), m_shift(0), m_numFreeDofs(1,0), m_numElimDofs(1,0), 
  m_numCpldDofs(1,0), m_curElimId(-1)
{ }

void gsDofMapper::localToGlobal(const gsMatrix<index_t>& locals,
                                index_t patchIndex,
                                gsMatrix<index_t>& globals,
                                index_t comp) const
{
    GISMO_ASSERT( locals.cols() == 1, "localToGlobal: Expecting one column of locals, got " << locals.cols() << ".");
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
        globals(i,0) = index(locals(i,0), patchIndex, comp);
}

void gsDofMapper::localToGlobal2(const gsMatrix<index_t>& locals,
                                 index_t patchIndex,
                                 gsMatrix<index_t>& globals,
                                 index_t & numFree,
                                 index_t comp) const
{
    GISMO_ASSERT( locals.cols() == 1, "localToGlobal: Expecting one column of locals");
    GISMO_ASSERT( &locals != &globals, "localToGlobal: Inplace not supported");
    const index_t numActive = locals.rows();
    globals.resize(numActive, 2);

    numFree = 0;
    index_t bot = numActive;
    for (index_t i = 0; i != numActive; ++i)
    {
      const index_t ii = index(locals(i,0), patchIndex, comp);
      if ( is_free_index(ii))
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

gsVector<index_t> gsDofMapper::asVector(index_t comp) const
{
  gsVector<index_t> v(m_dofs[comp].size());
  for(size_t j = 0; j!= m_dofs[comp].size(); ++j)
    v[j] = m_dofs[comp][j] + m_shift;
  return v;
}

void gsDofMapper::colapseDofs(index_t k, const gsMatrix<unsigned> & b,
			      index_t comp)
{
    const index_t last = b.size()-1;
    for ( index_t l=0; l!=last; ++l)
    {
        this->matchDof( k, b(l,0), k, b(l+1,0), comp);
    }
}

void gsDofMapper::matchDof(index_t u, index_t i,
			   index_t v, index_t j, index_t comp)
{
    //GISMO_ASSERT(u < m_offset.size() && i< m_offset[u+1] - m_offset[u] ,"Invalid patch-local dof "<<i<<" in matchDof.");
    //GISMO_ASSERT(v < m_offset.size() && j< m_offset[v+1] - m_offset[v] ,"Invalid patch-local dof "<<i<<" in matchDof.");
    //TODO: add check for correct range [by adding a last offset = number of total dofs]

    GISMO_ASSERT(static_cast<size_t>(u)<numPatches(), "Invalid patch index "<< u <<" >= "<< numPatches() );
    GISMO_ASSERT(static_cast<size_t>(v)<numPatches(), "Invalid patch index "<< v <<" >= "<< numPatches() );

    index_t d1 = MAPPER_PATCH_DOF(i,u,comp);
    index_t d2 = MAPPER_PATCH_DOF(j,v,comp);

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
	  mergeDofsGlobally(d1, d2, comp);  // both are eliminated, merge their indices
        else if (d2 == 0)
            MAPPER_PATCH_DOF(j,v, comp) = d1;   // second is free, eliminate it along with first
        else /* d2 > 0*/
            replaceDofGlobally(d2, d1, comp); // second is coupling, eliminate all instances of it
    }
    else if (d1 == 0)   // first dof is a free dof
    {
        if (d2 == 0)
        {
            MAPPER_PATCH_DOF(i,u,comp) = MAPPER_PATCH_DOF(j,v,comp) = m_numCpldDofs[1+comp]++;  // both are free, assign them a new coupling id
            if (u==v && i==j) return;
        }
        else if (d2 > 0)
            MAPPER_PATCH_DOF(i,u,comp) = d2;   // second is coupling, add first to the same coupling group
        else
            GISMO_ERROR("Something went terribly wrong");
    }
    else /* d1 > 0 */   // first dof is a coupling dof
    {
        GISMO_ASSERT(d2 > 0, "Something went terribly wrong");
        mergeDofsGlobally( d1, d2, comp);      // both are coupling dofs, merge them
    }

    // if we merged two different non-eliminated dofs, we lost one free dof
    if ( (d1 != d2 && (d1 >= 0 || d2 >= 0) ) || (d1 == 0 && d2 == 0) )
        --m_numFreeDofs[1+comp];
}

void gsDofMapper::matchDofs(index_t u, const gsMatrix<index_t> & b1,
                            index_t v,const gsMatrix<index_t> & b2,
			                index_t comp)
{
    const index_t sz = b1.size();
    GISMO_ASSERT( sz == b2.size(), "Waiting for same number of DoFs");
    for ( index_t k=0; k<sz; ++k)
      this->matchDof( u, b1(k,0), v, b2(k,0), comp );
}

void gsDofMapper::markCoupled(index_t i, index_t k, index_t comp)
{
    matchDof(k,i,k,i,comp);
}

void gsDofMapper::markTagged( index_t i, index_t k, index_t comp)
{
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );

    //see gsSortedVector::push_sorted_unique
    index_t t = index(i,k,comp);
    std::vector<index_t>::iterator pos = std::lower_bound(m_tagged.begin(), m_tagged.end(), t );

    if ( pos == m_tagged.end() || *pos != t )// If not found
        m_tagged.insert(pos, t);
}


void gsDofMapper::markBoundary(index_t k, const gsMatrix<index_t> & boundaryDofs, index_t comp)
{
    for (index_t i = 0; i < boundaryDofs.rows(); ++i)
      eliminateDof( boundaryDofs.at(i), k, comp );
}

void gsDofMapper::markCoupledAsTagged()
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    m_tagged.reserve(m_tagged.size()+m_numCpldDofs.back());
    std::vector<index_t>::const_iterator fr = m_numFreeDofs.begin()+1;
    std::vector<index_t>::const_iterator el = m_numElimDofs.begin();
    for(std::vector<index_t>::const_iterator cp = 
	  m_numCpldDofs.begin()+1;
	cp!=m_numCpldDofs.end(); ++cp, ++el, ++fr)
      {
	const index_t nd = *fr + *el;
	for(int i=0; i< *cp;++i)
	  m_tagged.push_back(nd+i);
      }
    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique(m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );
}

void gsDofMapper::eliminateDof( index_t i, index_t k, index_t comp)
{
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    const index_t old = MAPPER_PATCH_DOF(i,k,comp);
    if (old == 0)       // regular free dof
    {
        --m_numFreeDofs[comp+1];
        MAPPER_PATCH_DOF(i,k,comp) = m_curElimId--;
    }
    else if (old > 0)   // coupling dof
    {
        --m_numFreeDofs[comp+1];
        replaceDofGlobally( old, m_curElimId--, comp);//superfluous ElimId
    }
    // else: old < 0: already an eliminated dof, nothing to do
}

void gsDofMapper::finalize()
{
    GISMO_ASSERT(m_curElimId!=0, "Error in gsDofMapper::finalize() called twice.");

    for (size_t c = 0; c!=m_dofs.size(); ++c)
      {
	finalizeComp(c);
	
	//off-set
	m_curElimId -= m_numElimDofs[c+1];
	m_numFreeDofs[c+1] += m_numFreeDofs[c];
	m_numElimDofs[c+1] += m_numElimDofs[c];
	m_numCpldDofs[c+1] += m_numCpldDofs[c];
      }

    if ( 1!=m_dofs.size() )
      for (size_t c = 0; c!=m_dofs.size(); ++c)
	{
	  std::vector<index_t> & dofs = m_dofs[c];
	  for(std::vector<index_t>::iterator j = 
		dofs.begin(); j!= dofs.end(); ++j)
	    *j =  (*j<m_numFreeDofs[c+1]+m_numElimDofs[c] ?
		   *j - m_numElimDofs[c]                  :
		   *j - m_numFreeDofs[c+1] + m_numFreeDofs.back()
		   );
	}

    // Only bigger or equal to zero after finalize is called.
    m_curElimId = m_numFreeDofs.back();
    //if () m_curElimId += m_numElimDofs.back();
}

void gsDofMapper::finalizeComp(const index_t comp)
{
    std::vector<index_t> & dofs = m_dofs[comp];

    // For assigning coupling and eliminated dofs to continuous
    // indices (-1 = unassigned)
    std::vector<index_t> couplingDofs(m_numCpldDofs[comp+1] -1, -1);
    //std::vector<index_t> elimDofs(-m_curElimId - 1 - m_numElimDofs[comp], -1);
    std::map<index_t,index_t> elimDofs;
    // Free dofs start at offset
    index_t curFreeDof = m_numFreeDofs[comp]+m_numElimDofs[comp];
    // Eliminated dofs start after free dofs plus previous components
    index_t curElimDof = m_numFreeDofs[comp+1] + curFreeDof;

    // Coupling dofs start after standard dofs (=num of zeros in dofs)
    index_t curCplDof = std::count(dofs.begin(), dofs.end(), 0);
    // Devise number of coupled dofs (m_numCpldDofs was used as
    // coupling id up to here)
    m_numCpldDofs[comp+1] = m_numFreeDofs[comp+1] - curCplDof;
    curCplDof += curFreeDof; //off-set

    /*// For debugging: counting the number of coupled and boundary dofs
    std::vector<index_t> alldofs = dofs;
    std::sort( alldofs.begin(), alldofs.end() );
    alldofs.erase( std::unique( alldofs.begin(), alldofs.end() ), alldofs.end() );
    const index_t numCoupled =
    std::count_if( alldofs.begin(), alldofs.end(),
                          GS_BIND2ND(std::greater<index_t>(), 0) );
    const index_t numBoundary =
    std::count_if( alldofs.begin(), alldofs.end(),
                          GS_BIND2ND(std::less<index_t>(), 0) );
    */

    for (size_t k = 0; k < dofs.size(); ++k)
    {
        const index_t dofType = dofs[k];

        if (dofType == 0)       // standard dof
            dofs[k] = curFreeDof++;
        else if (dofType < 0)   // eliminated dof
        {
            const index_t id = -dofType-1;
	    if (elimDofs.find(id)==elimDofs.end())
	      elimDofs[id] = curElimDof++;
            dofs[k] = elimDofs[id];
        }
        else // dofType > 0     // coupling dof
        {
            const index_t id = dofType - 1;
            if (couplingDofs[id] < 0)
                couplingDofs[id] = curCplDof++;
            dofs[k] = couplingDofs[id];
        }
    }

    // Devise number of eliminated dofs
    m_numElimDofs[comp+1] = curElimDof - curCplDof;

    /*
    gsDebugVar(m_numFreeDofs[comp+1]);
    gsDebugVar(m_numCpldDofs[comp+1]);
    gsDebugVar(m_numElimDofs[comp+1]);
    */

    curCplDof -= m_numFreeDofs[comp]+m_numElimDofs[comp]; //de-off-set
    GISMO_ASSERT(curCplDof == m_numFreeDofs[1+comp],
                 "gsDofMapper::finalize() - computed number of coupling "
                 "dofs does not match allocated number, "<<curCplDof<<"!="<<m_numFreeDofs[comp+1]);

    curFreeDof -= m_numFreeDofs[comp]+m_numElimDofs[comp];//de-off-set
    GISMO_ASSERT(curFreeDof + m_numCpldDofs[comp+1] == m_numFreeDofs[comp+1],
                 "gsDofMapper::finalize() - computed number of free dofs "
                 "does not match allocated number");
}

std::ostream& gsDofMapper::print( std::ostream& os ) const
{
  os<<" Dofs: "<< this->size() 
    <<"\n components: "<< m_dofs.size()<<"\n";
    os<<" free: "<< this->freeSize() <<"\n";
    os<<" coupled: "<< this->coupledSize() <<"\n";
    os<<" tagged: "<< this->taggedSize() <<"\n";
    os<<" elim: "<< this->boundarySize() <<"\n";
    if ( 1!=m_dofs.size() )
      {
	os<<" Free per comp: "<< gsAsConstVector<index_t>(m_numFreeDofs).transpose() <<"\n";
	os<<" Elim per comp: "<< gsAsConstVector<index_t>(m_numElimDofs).transpose() <<"\n";
	os<<" Cpld per comp: "<< gsAsConstVector<index_t>(m_numCpldDofs).transpose() <<"\n";
      }

    return os;
}

  void gsDofMapper::setIdentity(index_t nPatches, size_t nDofs,
				size_t nComp)
{
    m_curElimId   = -1;
    m_shift = m_bshift = 0;
    m_numFreeDofs.assign(nComp+1,nDofs); m_numFreeDofs.front()=0;
    m_numElimDofs.assign(nComp+1,0);
    m_numCpldDofs.assign(nComp+1,1); m_numCpldDofs.front()=0;
    
    //todo: check nDofs%nPatches==0 and initialize correctly
    m_offset.resize(nPatches, 0);

    m_dofs.resize(nComp, std::vector<index_t>(nDofs, 0));
}

  void gsDofMapper::permuteFreeDofs(const gsVector<index_t>& permutation, index_t comp)
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(m_numFreeDofs[comp+1] == permutation.size(), "permutation size does not match number of free dofs");
    //GISMO_ASSERT(m_tagged.empty(), "you cannot permute the dofVector twice, combine the permutation");
    //We could check here, that permutation is indeed a permutation

    //make a copy of the old ordering, easiest way to implement the permutation. Inplace reordering is quite hard.
    std::vector<index_t> dofs = m_dofs[comp];
    //first use a tempory coupled vector, other wise is_coupled_index() would not work correctly
    std::vector<index_t> tagged_permuted;
    tagged_permuted.reserve(taggedSize()); //reserve enough memory

    for(index_t i=0; i<(index_t)dofs.size();++i)
    {
        const index_t idx = dofs[i];
        if(is_free_index(idx))
        {
            m_dofs[comp][i] = permutation[idx];
            //fill  bookkeeping for tagged dofs
            if(is_tagged_index(idx))
                tagged_permuted.push_back(m_dofs[comp][i]);
        }
        else if(is_tagged_index(idx)) //Take care about eliminated tagged dofs
            tagged_permuted.push_back(idx);
    }
    m_tagged.swap(tagged_permuted);

    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique (m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );

    //coupled dofs cannot be tracked anymore on this component
    m_numCpldDofs[comp+1]=m_numCpldDofs[comp];
    for(std::vector<index_t>::iterator s=m_numCpldDofs.begin()+comp+2;
	s<m_numCpldDofs.end(); ++s)
      *s += (*s-1);
}


  void gsDofMapper::initPatchDofs(const gsVector<index_t> & patchDofSizes, index_t nComp)
{
    m_curElimId   = -1;
    m_shift = m_bshift = 0;
    m_numElimDofs.assign(nComp+1,0);
    m_numCpldDofs.assign(nComp+1,1); m_numCpldDofs.front()=0;

    const size_t nPatches = patchDofSizes.size();

    // Initialize offsets and dof holder
    m_offset.reserve( nPatches );
    m_offset.push_back(0);
    for (size_t k = 1; k < nPatches; ++k)
    {
        m_offset.push_back( m_offset.back() + patchDofSizes[k-1] );
    }

    m_numFreeDofs.assign(nComp+1,
    m_offset.back() + patchDofSizes[nPatches-1]);
    m_numFreeDofs.front()=0;

    m_dofs.resize(nComp, std::vector<index_t>(m_numFreeDofs.back(), 0));
}

void gsDofMapper::replaceDofGlobally(index_t oldIdx, index_t newIdx)
{
  for(size_t i = 0; i!= m_dofs.size(); ++i)
    std::replace(m_dofs[i].begin(), m_dofs[i].end(), oldIdx, newIdx );
}

  void gsDofMapper::replaceDofGlobally(index_t oldIdx, index_t newIdx, index_t comp)
{
    std::vector<index_t> & dofs = m_dofs[comp];
    std::replace(dofs.begin(), dofs.end(), oldIdx, newIdx );
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

  void gsDofMapper::mergeDofsGlobally(index_t dof1, index_t dof2, index_t comp)
{
    if (dof1 != dof2)
    {
        if (dof1 < dof2)
            std::swap(dof1, dof2);
        replaceDofGlobally(dof1, dof2, comp);
    }
}

void gsDofMapper::preImage(const index_t gl,
                           std::vector<std::pair<index_t,index_t> > & result) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    typedef std::vector<index_t>::const_iterator citer;
    const std::vector<index_t> & dofs = m_dofs[componentOf(gl)];
    result.clear();
    size_t cur = 0;//local offsetted index

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
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

gsVector<index_t> gsDofMapper::inverseAsVector(index_t comp) const
{
    GISMO_ASSERT(isPermutation(), "This dofMapper is not 1-1");
    gsVector<index_t> v(size());
      for(size_t j = 0; j!= m_dofs[comp].size(); ++j)
	v[ m_dofs[comp][j] ] = j;
    return v;
}

std::map<index_t,index_t>
gsDofMapper::inverseOnPatch(const index_t k) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );

    std::map<index_t,index_t> inv;
    //inv.reserve(patchSize(k));
    typedef std::vector<index_t>::const_iterator citer;

    for(size_t i = 0; i!= m_dofs.size(); ++i)
      {
	citer it = m_dofs[i].begin()+m_offset[k];
	for(size_t j = 0; j!= m_dofs[i].size(); ++j,++it)
	  inv[*it]=j;
      }
    return inv;
}

bool gsDofMapper::indexOnPatch(const index_t gl, const index_t k) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    const std::vector<index_t> & dofs = m_dofs[componentOf(gl)];
    const citer istart = dofs.begin()+m_offset[k];
    const citer iend   = istart + patchSize(k);
    return (std::find(istart, iend, gl)!=iend);
}

index_t gsDofMapper::boundarySizeWithDuplicates() const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");

    index_t res = 0;
    for (size_t i = 0; i!= m_dofs.size(); ++i)
      res += std::count_if(m_dofs[i].begin(), m_dofs[i].end(),
			   GS_BIND2ND(std::greater<index_t>(),
					freeSize(i) - 1) );
    return res;
}

index_t gsDofMapper::coupledSize() const
{
    GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    return m_numCpldDofs.back();
/*// Implementation without saving this number:
    // Property: coupled (eliminated or not) DoFs appear more than once in the mapping.
    GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");

    std::vector<index_t> CountMap(m_numFreeDofs,0);

    // Count number of appearances of each free DoF
    for (std::vector<index_t>::const_iterator it = m_dofs.begin(); it != m_dofs.end(); ++it)
        if ( *it < m_numFreeDofs )
            CountMap[*it]++;

    // Count the number of freeDoFs that appear more than once
    return std::count_if( CountMap.begin(), CountMap.end(),
                          GS_BIND2ND(std::greater<index_t>(), 1) );
*/
}

index_t gsDofMapper::taggedSize() const
{
    return m_tagged.size();
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
