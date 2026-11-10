/** @file gsDofMapper.h

    @brief Provides the gsDofMapper class for re-indexing DoFs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#include <gsCore/gsDofMapper.h>
#include <limits>


namespace gismo
{

gsDofMapper::gsDofMapper() :
  m_offset(2,0), m_nPatches(1), m_shift(0), m_numFreeDofs(1,0), m_numElimDofs(1,0),
  m_numCpldDofs(1,0), m_curElimId(-1)
{ }

gsDofMapper::gsDofMapper(const gsVector<index_t> & patchDofSizes, index_t nComp)
{
    const gsMatrix<index_t> sz = patchDofSizes.replicate(1, nComp);
    initPatchDofs(sz);
}

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

    GISMO_ASSERT(comp>-1,"Component is invalid");
    for (index_t i = 0; i < numActive; ++i)
        globals(i,0) = index(locals(i,0), patchIndex, comp);
}

void gsDofMapper::localToGlobal2(const gsMatrix<index_t>& locals,
                                 index_t patchIndex,
                                 gsMatrix<index_t>& globals,
                                 index_t & numFree,
                                 index_t comp) const
{
    GISMO_ASSERT(comp>-1,"Component is invalid");
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
    GISMO_ASSERT(comp>-1,"Component is invalid");
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
    GISMO_ASSERT(static_cast<size_t>(u)<numPatches(), "Invalid patch index "<< u <<" >= "<< numPatches() );
    GISMO_ASSERT(static_cast<size_t>(v)<numPatches(), "Invalid patch index "<< v <<" >= "<< numPatches() );

    if (-1==comp)
    {
        for (index_t c = 0; static_cast<size_t>(c) != componentsSize(); ++c)
            GISMO_ENSURE(patchSize(u,c)==patchSize(u,0) && patchSize(v,c)==patchSize(v,0),
                         "matchDof: comp=-1 needs patches "<<u<<" and "<<v<<" to have the same size in every component; call matchDof once per component");
        for (index_t c = 0; static_cast<size_t>(c) != componentsSize(); ++c)
            matchDof(u,i,v,j,c);
        return;
    }

    GISMO_ASSERT(off(u,comp)+i < off(m_nPatches,comp), "Invalid patch-local dof "<<i<<" of patch "<<u<<" in matchDof.");
    GISMO_ASSERT(off(v,comp)+j < off(m_nPatches,comp), "Invalid patch-local dof "<<j<<" of patch "<<v<<" in matchDof.");

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
    GISMO_ASSERT(comp>-1,"Component is invalid");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );

    //see gsSortedVector::push_sorted_unique
    index_t t = MAPPER_PATCH_DOF(i,k,comp);
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
    for(size_t c = 0; c != m_dofs.size(); ++c)
      {
	const index_t cp = m_numCpldDofs[c+1] - m_numCpldDofs[c];
	for(index_t i = m_numFreeDofs[c+1]-cp; i < m_numFreeDofs[c+1]; ++i)
	  m_tagged.push_back(i);
      }
    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique(m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );
}

void gsDofMapper::eliminateDof( index_t i, index_t k, index_t comp)
{
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    GISMO_ASSERT(static_cast<size_t>(comp)<componentsSize() || comp == -1, "Invalid component index "<< comp <<" >= "<< componentsSize() );

    if (-1==comp)
    {
        for (index_t c = 0; static_cast<size_t>(c) != componentsSize(); ++c)
            GISMO_ENSURE(patchSize(k,c)==patchSize(k,0),
                         "eliminateDof: comp=-1 needs patch "<<k<<" to have the same size in every component; call eliminateDof once per component");
        for (index_t c = 0; static_cast<size_t>(c) != componentsSize(); ++c)
            eliminateDof(i,k,c);
        return;
    }

    GISMO_ASSERT(off(k,comp)+i < off(m_nPatches,comp), "Invalid patch-local dof "<<i<<" of patch "<<k<<" in eliminateDof.");

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
    GISMO_ASSERT(m_curElimId<0, "Error in gsDofMapper::finalize() called twice.");

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
    m_shift = 0;
    m_numFreeDofs.assign(nComp+1,nDofs); m_numFreeDofs.front()=0;
    m_numElimDofs.assign(nComp+1,0);
    m_numCpldDofs.assign(nComp+1,1); m_numCpldDofs.front()=0;

    //todo: check nDofs%nPatches==0 and initialize correctly
    m_nPatches = nPatches;
    m_offset.assign(nComp*(nPatches+1), 0);
    for (size_t c = 0; c != nComp; ++c)
        m_offset[c*(nPatches+1)+nPatches] = nDofs;

    m_dofs.assign(nComp, std::vector<index_t>(nDofs, 0));
    m_tagged.clear();
}

  void gsDofMapper::initPatchDofs(const gsMatrix<index_t> & sz)
{
    GISMO_ENSURE(sz.rows()>0 && sz.cols()>0 && sz.minCoeff()>=0,
                 "initPatchDofs: Expected a non-empty nPatches x nComp size matrix with non-negative entries, got "
                 << sz.rows() << " x " << sz.cols() << ".");

    const index_t nP = sz.rows();
    const index_t nC = sz.cols();

    // Global indices span all components; coupling ids count from 1
    const index_t maxDofs = std::numeric_limits<index_t>::max() - 1;
    index_t total = 0;
    for (index_t c = 0; c != nC; ++c)
        for (index_t k = 0; k != nP; ++k)
        {
            GISMO_ENSURE(sz(k,c) <= maxDofs - total,
                         "initPatchDofs: The total number of dofs exceeds the range of index_t ("<<maxDofs<<").");
            total += sz(k,c);
        }

    m_curElimId   = -1;
    m_shift = 0;
    m_numElimDofs.assign(nC+1,0);
    m_numCpldDofs.assign(nC+1,1); m_numCpldDofs.front()=0;

    m_nPatches = nP;
    m_offset.assign(nC*(nP+1), 0);
    m_numFreeDofs.assign(nC+1, 0);
    m_dofs.assign(nC, std::vector<index_t>());
    for (index_t c = 0; c != nC; ++c)
    {
        for (index_t k = 0; k != nP; ++k)
            m_offset[c*(nP+1)+k+1] = m_offset[c*(nP+1)+k] + sz(k,c);
        m_numFreeDofs[c+1] = static_cast<index_t>(off(nP,c));
        m_dofs[c].assign(off(nP,c), 0);
    }
}

void gsDofMapper::replaceDofGlobally(index_t oldIdx, index_t newIdx)
{
  for(size_t i = 0; i!= m_dofs.size(); ++i)
    std::replace(m_dofs[i].begin(), m_dofs[i].end(), oldIdx, newIdx );
}

  void gsDofMapper::replaceDofGlobally(index_t oldIdx, index_t newIdx, index_t comp)
{
    GISMO_ASSERT(comp>-1,"Component is invalid");
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
    const index_t g = gl - m_shift;
    const index_t c = componentOf(gl);
    const std::vector<index_t> & dofs = m_dofs[c];
    const std::vector<size_t>::const_iterator ob = m_offset.begin() + c*(m_nPatches+1);
    result.clear();
    size_t cur = 0;//local offsetted index

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
    {
        if ( *it == g )
        {
            // Get the patch index of "cur" by "un-offsetting"
            const index_t patch = std::upper_bound(ob, ob + m_nPatches, cur) - ob - 1;

            // Found a patch-dof pair
            result.push_back( std::make_pair(patch, cur - off(patch,c)) );
        }
    }
}

std::pair<index_t,index_t> gsDofMapper::anyPreImage(const index_t gl) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    typedef std::vector<index_t>::const_iterator citer;
    const index_t g = gl - m_shift;
    const index_t c = componentOf(gl);
    const std::vector<index_t> & dofs = m_dofs[c];
    const std::vector<size_t>::const_iterator ob = m_offset.begin() + c*(m_nPatches+1);
    size_t cur = 0;//local offsetted index

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
    {
        if ( *it == g )
        {
            // Get the patch index of "cur" by "un-offsetting"
            const index_t patch = std::upper_bound(ob, ob + m_nPatches, cur) - ob - 1;

            // Found a patch-dof pair
            return std::make_pair(patch, cur - off(patch,c));
        }
    }
    GISMO_ERROR("The global index "<< gl <<" is not valid");
}

std::vector<std::pair<index_t,index_t> > gsDofMapper::anyPreImages(index_t comp) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    typedef std::vector<index_t>::const_iterator citer;
    const std::vector<index_t> & dofs = m_dofs[comp];
    const std::vector<size_t>::const_iterator ob = m_offset.begin() + comp*(m_nPatches+1);
    size_t cur = 0;//local offsetted index

    std::vector<std::pair<index_t,index_t> > result(size(), {-1,0});

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
    {
        if ( -1 == result[*it].first )
        {
            // Get the patch index of "cur" by "un-offsetting"
            const index_t patch = std::upper_bound(ob, ob + m_nPatches, cur) - ob - 1;

            // Found a patch-dof pair
            result[*it] = std::make_pair(patch, cur - off(patch,comp));
        }
    }
    return result;
}

gsVector<index_t> gsDofMapper::inverseAsVector(index_t comp) const
{
    GISMO_ASSERT(comp>-1,"Component is invalid");
    GISMO_ASSERT(isPermutation(), "This dofMapper is not 1-1");
    gsVector<index_t> v;
    v.setConstant(size(), -1);
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
        citer it = m_dofs[i].begin()+off(k, static_cast<index_t>(i));
        for(size_t j = 0; j!= patchSize(k, static_cast<index_t>(i)); ++j,++it)
            inv[*it + m_shift]=j;
    }
    return inv;
}

bool gsDofMapper::indexOnPatch(const index_t gl, const index_t k, index_t & local) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    const index_t c = componentOf(gl);
    const std::vector<index_t> & dofs = m_dofs[c];
    const citer istart = dofs.begin()+off(k,c);
    const citer iend   = istart + patchSize(k,c);
    auto it = std::find(istart, iend, gl - m_shift);
    if (iend==it) return false;
    local = std::distance(istart,it);
    return true;
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
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    const index_t s = m_numFreeDofs.back() - 1;
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[0].begin() + off(k,0);
    citer iend   = istart + patchSize(k);
    return find_impl(istart, iend, GS_BIND2ND(std::greater<index_t>(),s));
}

gsVector<index_t> gsDofMapper::findFree(const index_t k) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    const index_t s = m_numFreeDofs.back();
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[0].begin() + off(k,0);
    citer iend   = istart + patchSize(k);
    return find_impl(istart, iend, GS_BIND2ND(std::less<index_t>(),s));
}

namespace {
struct _isBetween
{
    _isBetween(const index_t l, const index_t u) : _l(l), _u(u) { }
    index_t _l, _u;
    bool operator()(const index_t i) { return  (i < _u) && (i > _l); }
};
} // end anonymous namespace

gsVector<index_t> gsDofMapper::findCoupled(const index_t k, const index_t j) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );

    if (k==j) return gsVector<index_t>();

    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[0].begin() + off(k,0);
    citer iend   = istart + patchSize(k);

    const index_t l = m_numFreeDofs[1]-(m_numCpldDofs[1]-m_numCpldDofs[0])-1;
    const index_t u = m_numFreeDofs[1];
    if (-1==j)
        return find_impl(istart, iend, _isBetween(l,u) );
    else
    {
        citer istartj = m_dofs[0].begin() + off(j,0);
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
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    const citer istart = m_dofs[0].begin() + off(k,0);
    const citer iend   = istart + patchSize(k);
    return find_impl(istart, iend,
                     _isBetween(-1, m_numFreeDofs[1]-(m_numCpldDofs[1]-m_numCpldDofs[0])) );
}

gsVector<index_t> gsDofMapper::findTagged(const index_t k) const
{
    GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
    GISMO_ASSERT(static_cast<size_t>(k)<numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[0].begin() + off(k,0);
    citer iend   = istart + patchSize(k);
    return find_impl(istart, iend,
                     [this](const index_t v)
                     { return std::binary_search(m_tagged.begin(), m_tagged.end(), v); });
}
void gsDofMapper::setShift (index_t shift)
{
    m_shift=shift;
}

void gsDofMapper::addShift (index_t shift)
{
    m_shift+=shift;
}

#ifdef GISMO_WITH_PYBIND11
namespace py = pybind11;
void pybind11_init_gsDofMapper(py::module &m)
{
    using Class = gsDofMapper;
    py::class_<Class>(m, "gsDofMapper")

    // Member functions
    .def("asVector",                    &Class::asVector,                   "Returns a vector taking flat local indices to global")
    .def("inverseAsVector",             &Class::inverseAsVector,            "Returns a vector taking global indices to flat local")
    .def("colapseDofs",                 &Class::colapseDofs,                "Calls matchDof() for all dofs on the given patch side i ps. Thus, the whole set of dofs collapses to a single global dof")
    .def("matchDof", &Class::matchDof, "Couples dof \a i of patch \a u with dof \a j of patch \a v such that they refer to the same global dof at component \a comp.")
    .def("matchDofs", &Class::matchDofs, "Couples dofs \a b1 of patch \a u with dofs \a b2 of patch \a v one by one such that they refer to the same global dof.")
    .def("markCoupled", &Class::markCoupled, "Mark the local dof \a i of patch \a k as coupled.")
    .def("markTagged", &Class::markTagged, "Mark a local dof \a i of patch \a k as tagged")
    .def("markCoupledAsTagged", &Class::markCoupledAsTagged, "Mark all coupled dofs as tagged")
    .def("markBoundary", &Class::markBoundary, "Mark the local dofs \a boundaryDofs of patch \a k as eliminated. ")
    .def("eliminateDof", &Class::eliminateDof, "Mark the local dof \a i of patch \a k as eliminated.")
    .def("finalize", &Class::finalize, "Must be called after all boundaries and interfaces have been marked to set up the dof numbering.")
    .def("isFinalized", &Class::isFinalized, "Checks whether finalize() has been called.")
    .def("isPermutation", &Class::isPermutation, "Returns true iff the mapper is a permuatation")

    .def("setIdentity", &Class::setIdentity, "Set this mapping to be the identity")
    .def("setShift", &Class::setShift, "Set the shift amount for the global numbering")
    .def("addShift", &Class::addShift, "Add a shift amount to the global numbering")

    .def("index", &Class::index, "Returns the global dof index associated to local dof \a i of patch \a k.")
    .def("bindex", &Class::bindex, "Returns the boundary index of local dof \a i of patch \a k.")
    .def("cindex", &Class::cindex, "Returns the coupled dof index")
    .def("tindex", &Class::tindex, "Returns the tagged dof index")
    .def("global_to_bindex", &Class::global_to_bindex, "Returns the boundary index of global dof \a gl.")
    .def("is_free_index", &Class::is_free_index, "Returns true if global dof \a gl is not eliminated.")
    .def("is_free", &Class::is_free, "Returns true if local dof \a i of patch \a k is not eliminated.")
    .def("is_boundary_index", &Class::is_boundary_index, "Returns true if global dof \a gl is eliminated")
    .def("is_boundary", &Class::is_boundary, "Returns true if local dof \a i of patch \a k is eliminated.")
    .def("is_coupled", &Class::is_coupled, "Returns true if local dof \a i of patch \a k is coupled.")
    .def("is_coupled_index", &Class::is_coupled_index, "Returns true if \a gl is a coupled dof.")
    .def("is_tagged", &Class::is_tagged, "Returns true if local dof \a i of patch \a k is tagged.")
    .def("is_tagged_index", &Class::is_tagged_index, "Returns true if \a gl is a tagged dof.")
    .def("numComponents", &Class::numComponents, "Returns the number of components present in the mapper")
    .def("size", static_cast<index_t (Class::*)()        const > (&Class::size), "Returns the total number of dofs (free and eliminated).")
    .def("size", static_cast<index_t (Class::*)(index_t) const > (&Class::size), "Returns the total number of dofs (free and eliminated).")
    .def("freeSize", static_cast<index_t (Class::*)()        const > (&Class::freeSize), "Returns the number of free (not eliminated) dofs.")
    .def("freeSize", static_cast<index_t (Class::*)(index_t) const > (&Class::freeSize), "Returns the number of free (not eliminated) dofs.")
    .def("coupledSize", &Class::coupledSize, "Returns the number of coupled (not eliminated) dofs.")
    .def("taggedSize", &Class::taggedSize, "Returns the number of tagged dofs.")
    .def("boundarySize", &Class::boundarySize, "Returns the number of eliminated dofs.")

    .def("offset", &Class::offset, "Returns the offset of patch \a k in component \a c", py::arg("k"), py::arg("c") = 0)
    .def("numPatches", &Class::numPatches, "Returns the number of patches present underneath the mapper")
    .def("mapSize", &Class::mapSize, "Returns the total number of patch-local degrees of freedom that are being mapped")
    .def("componentsSize", &Class::componentsSize, "Returns the components size")
    .def("patchSize", &Class::patchSize, "Returns the total number of patch-local DoFs that live on patch \a k for component \a c")
    .def("totalSize", &Class::totalSize, "Returns the total size of the mapper")
    .def("indexOnPatch", static_cast<bool (Class::*)(index_t,index_t) const > (&Class::indexOnPatch), "For \a gl being a global index, this function returns true whenever \a gl corresponds to patch \a k")
    .def("__str__",
         [] (Class & self)
        {
            std::ostringstream os;
            self.print(os);
            return os.str();
        },
        "Returns a string with information about the object.");
}
#endif
} // namespace gismo
