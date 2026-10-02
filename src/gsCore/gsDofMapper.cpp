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

namespace {

/// The largest dof count a gsDofMapper can represent.
///
/// Callers hand in sizes as \c size_t, but every count the mapper keeps
/// ends up in an \c index_t: \c m_numFreeDofs / \c m_numElimDofs hold the
/// per-component totals, \c finalize() accumulates those into running sums
/// in the same type, and the global indices it hands out are \c index_t as
/// well.  \c index_t is a build-time configurable type that is plain \c int
/// by default and may be narrower still, so a size that fits a \c size_t is
/// no evidence that it fits the mapper.  Sizes must therefore be rejected
/// against this bound -- before anything is narrowed or allocated.
inline size_t maxDofCount()
{
    // sizeof(index_t) <= sizeof(size_t): index_t's signed maximum is exactly
    // representable as a size_t.  Otherwise every size_t fits an index_t and
    // size_t's own maximum is the only real bound.
    return sizeof(index_t) <= sizeof(size_t)
        ? static_cast<size_t>(std::numeric_limits<index_t>::max())
        : std::numeric_limits<size_t>::max();
}

/// The largest component count a gsDofMapper can represent.
///
/// numComponents() and every component argument are \c index_t, and
/// m_numFreeDofs & co. are indexed by component+1 in that type, so a count
/// handed in as \c size_t must fit it.  Zero-sized components pass every
/// dof-count check, so this is not implied by maxDofCount() being checked.
inline size_t maxComponentCount()
{
    return maxDofCount();
}

} // anonymous namespace

gsDofMapper::gsDofMapper() :
  m_nPatches(1), m_layout(PatchConcatenated), m_hasDistinctComponentSpaces(false),
  m_uniformComponents(true),
  m_shift(0), m_bshift(0), m_numFreeDofs(1,0), m_numElimDofs(1,0),
  m_numCpldDofs(1,0), m_curElimId(-1)
{
    checkInvariants();
}

gsDofMapper::gsDofMapper(const std::vector<gsVector<index_t> > & patchDofSizes,
                          bool hasDistinctComponentSpaces)
{
    initRaggedPatchDofs(patchDofSizes, hasDistinctComponentSpaces);
}

void gsDofMapper::checkInvariants() const
{
#ifndef NDEBUG
    GISMO_ASSERT(m_offset.size() == m_dofs.size()*(m_nPatches+1),
                 "gsDofMapper: offset table size "<<m_offset.size()<<" does not match "
                 <<m_dofs.size()<<" components x "<<(m_nPatches+1)<<" (nPatches+1).");
    GISMO_ASSERT(m_numFreeDofs.size() == m_dofs.size()+1,
                 "gsDofMapper: m_numFreeDofs has the wrong size.");
    GISMO_ASSERT(m_numElimDofs.size() == m_dofs.size()+1,
                 "gsDofMapper: m_numElimDofs has the wrong size.");
    GISMO_ASSERT(m_numCpldDofs.size() == m_dofs.size()+1,
                 "gsDofMapper: m_numCpldDofs has the wrong size.");

    for (size_t c = 0; c != m_dofs.size(); ++c)
    {
        const index_t cc = static_cast<index_t>(c);
        GISMO_ASSERT(offAt(cc,0) == 0,
                     "gsDofMapper: offset table does not start at 0 for component "<<c<<".");
        GISMO_ASSERT(offAt(cc, static_cast<index_t>(m_nPatches)) == m_dofs[c].size(),
                     "gsDofMapper: offset sentinel does not match storage size for component "
                     <<c<<": "<<offAt(cc, static_cast<index_t>(m_nPatches))<<" != "<<m_dofs[c].size()<<".");

        if (PatchConcatenated == m_layout)
        {
            for (size_t k = 0; k != m_nPatches; ++k)
                GISMO_ASSERT(offAt(cc,static_cast<index_t>(k)) <= offAt(cc,static_cast<index_t>(k+1)),
                             "gsDofMapper: offsets are not monotone for component "<<c<<" at patch "<<k<<".");
        }
        else // GlobalIdentity
        {
            for (size_t k = 0; k != m_nPatches; ++k)
                GISMO_ASSERT(offAt(cc,static_cast<index_t>(k)) == 0,
                             "gsDofMapper: aliased-identity component "<<c
                             <<" has a nonzero offset on real patch "<<k<<".");
        }
    }

    GISMO_ASSERT(m_uniformComponents == computeUniformComponents(),
                 "gsDofMapper: the cached hasUniformComponents() is stale.");

    if (m_curElimId >= 0) // finalized: the count vectors are now cumulative prefix sums
    {
        for (size_t c = 0; c != m_dofs.size(); ++c)
        {
            GISMO_ASSERT(m_numFreeDofs[c] <= m_numFreeDofs[c+1],
                         "gsDofMapper: m_numFreeDofs is not monotone after finalize().");
            GISMO_ASSERT(m_numElimDofs[c] <= m_numElimDofs[c+1],
                         "gsDofMapper: m_numElimDofs is not monotone after finalize().");
            GISMO_ASSERT(m_numCpldDofs[c] <= m_numCpldDofs[c+1],
                         "gsDofMapper: m_numCpldDofs is not monotone after finalize().");
        }
    }
#endif // NDEBUG
}

bool gsDofMapper::computeUniformComponents() const
{
    if (m_hasDistinctComponentSpaces)
        return false;

    const index_t nComp = numComponents();
    if (nComp <= 1)
        return true;

    if (GlobalIdentity == m_layout)
    {
        for (index_t c = 1; c != nComp; ++c)
            if (totalSize(c) != totalSize(0))
                return false;
    }
    else // PatchConcatenated
    {
        for (index_t p = 0; p != static_cast<index_t>(m_nPatches); ++p)
            for (index_t c = 1; c != nComp; ++c)
                if (patchSize(p,c) != patchSize(p,0))
                    return false;
    }
    return true;
}

void gsDofMapper::localToGlobal(const gsMatrix<index_t>& locals,
                                index_t patchIndex,
                                gsMatrix<index_t>& globals,
                                index_t comp) const
{
    // Debug-only, like index() itself (see dofAt()): this runs once per
    // element and component in every assembler, and even O(1) checks per
    // call are measurable against a loop over a handful of local dofs.
    // Empty means no rows: the loop below runs over rows(), so an N x 0
    // matrix would be read at (i,0) in storage it does not have.
    GISMO_ASSERT( locals.cols() == 1 || locals.rows() == 0,
                  "localToGlobal: Expecting one column of locals, got " << locals.cols() << ".");
    GISMO_ASSERT( validPatch(patchIndex), "localToGlobal: invalid patch "<<patchIndex
                  <<", the mapper has "<<m_nPatches<<" patches.");
    GISMO_ASSERT( validComponent(comp), "localToGlobal: invalid component "<<comp
                  <<", the mapper has "<<m_dofs.size()<<" components.");
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
    // Debug-only, as in localToGlobal().
    GISMO_ASSERT( locals.cols() == 1 || locals.rows() == 0,
                  "localToGlobal2: Expecting one column of locals, got " << locals.cols() << ".");
    // globals is resized to two columns before locals is read, so an aliased
    // call would destroy its own input.
    GISMO_ASSERT( &locals != &globals, "localToGlobal2: Inplace not supported");
    GISMO_ASSERT( validPatch(patchIndex), "localToGlobal2: invalid patch "<<patchIndex
                  <<", the mapper has "<<m_nPatches<<" patches.");
    GISMO_ASSERT( validComponent(comp), "localToGlobal2: invalid component "<<comp
                  <<", the mapper has "<<m_dofs.size()<<" components.");
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
    ensureComponent(comp, "asVector");
  gsVector<index_t> v(m_dofs[comp].size());
  for(size_t j = 0; j!= m_dofs[comp].size(); ++j)
    v[j] = m_dofs[comp][j] + m_shift;
  return v;
}

void gsDofMapper::colapseDofs(index_t k, const gsMatrix<unsigned> & b,
			      index_t comp)
{
    ensureNotFinalized("colapseDofs");
    ensureComponentOrAll(comp, "colapseDofs");
    ensurePatch(k, "colapseDofs");
    // b(l,0) addresses the first column only, so a matrix with more than one
    // column would be read past its first column's end, and one with rows but
    // no column at all would be read in storage it does not have.
    GISMO_ENSURE(b.cols() == 1 || b.rows() == 0,
                 "gsDofMapper::colapseDofs: Expecting one column of dofs, got "<<b.cols()<<".");
    // Fewer than two dofs collapse to nothing.  (b.size()-1 as a loop bound
    // is -1 for an empty b, which the loop below would never reach.)
    const index_t n = b.rows();
    if (n < 2) return;

    // Every entry is validated before the first dof is matched, so that a bad
    // entry cannot leave the collapse half-applied.  The entries are unsigned:
    // compared as such, a value beyond index_t's range is rejected instead of
    // being narrowed onto a valid-looking index.
    for (index_t l = 0; l != n; ++l)
        for (index_t c = (-1 == comp ? 0 : comp);
             c != (-1 == comp ? numComponents() : comp+1); ++c)
            GISMO_ENSURE(static_cast<size_t>(b(l,0)) < localCount(k,c),
                         "gsDofMapper::colapseDofs: invalid local dof "<<b(l,0)<<" of patch "<<k
                         <<", component "<<c<<localRangeInfo(k,c));

    for (index_t l = 0; l+1 != n; ++l)
    {
        const index_t i = static_cast<index_t>(b(l,0)), j = static_cast<index_t>(b(l+1,0));
        if (-1 == comp)
            for (index_t c = 0; c != numComponents(); ++c)
                matchDofImpl(k, i, k, j, c);
        else
            matchDofImpl(k, i, k, j, comp);
    }
}

void gsDofMapper::matchDof(index_t u, index_t i,
			   index_t v, index_t j, index_t comp)
{
    ensureNotFinalized("matchDof");
    ensureComponentOrAll(comp, "matchDof");
    ensurePatch(u, "matchDof");
    ensurePatch(v, "matchDof");
    // Both dofs, in every targeted component, before anything is matched.
    ensureLocal(i, u, comp, "matchDof");
    ensureLocal(j, v, comp, "matchDof");

    if (-1==comp)
    {
        for (index_t c = 0; c != numComponents(); ++c)
            matchDofImpl(u,i,v,j,c);
        return;
    }
    matchDofImpl(u,i,v,j,comp);
}

void gsDofMapper::matchDofImpl(index_t u, index_t i,
                               index_t v, index_t j, index_t comp)
{
    index_t d1 = dofAt(i,u,comp);
    index_t d2 = dofAt(j,v,comp);

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
            dofAt(j,v, comp) = d1;   // second is free, eliminate it along with first
        else /* d2 > 0*/
            replaceDofGlobally(d2, d1, comp); // second is coupling, eliminate all instances of it
    }
    else if (d1 == 0)   // first dof is a free dof
    {
        if (d2 == 0)
        {
            dofAt(i,u,comp) = dofAt(j,v,comp) = ++m_numCpldDofs[1+comp];  // both are free, assign them a new coupling id
            if (u==v && i==j) return;
        }
        else if (d2 > 0)
            dofAt(i,u,comp) = d2;   // second is coupling, add first to the same coupling group
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
    ensureNotFinalized("matchDofs");
    ensureComponentOrAll(comp, "matchDofs");
    ensurePatch(u, "matchDofs");
    ensurePatch(v, "matchDofs");
    const index_t sz = b1.size();
    GISMO_ENSURE( sz == b2.size(), "gsDofMapper::matchDofs: Waiting for same number of DoFs, got "
                  <<sz<<" and "<<b2.size()<<".");
    // b(k,0) addresses the first column only; see colapseDofs().  Unlike
    // there, the loops run over size() entries, so an empty matrix of any
    // shape is safe.
    GISMO_ENSURE( (b1.cols() == 1 || sz == 0) && (b2.cols() == 1 || sz == 0),
                  "gsDofMapper::matchDofs: Expecting one column of dofs, got "
                  <<b1.cols()<<" and "<<b2.cols()<<".");
    // Every pair, in every targeted component, before the first is matched.
    for ( index_t k=0; k<sz; ++k)
    {
        ensureLocal(b1(k,0), u, comp, "matchDofs");
        ensureLocal(b2(k,0), v, comp, "matchDofs");
    }

    for ( index_t k=0; k<sz; ++k)
    {
        if (-1 == comp)
            for (index_t c = 0; c != numComponents(); ++c)
                matchDofImpl(u, b1(k,0), v, b2(k,0), c);
        else
            matchDofImpl(u, b1(k,0), v, b2(k,0), comp);
    }
}

void gsDofMapper::markCoupled(index_t i, index_t k, index_t comp)
{
    matchDof(k,i,k,i,comp);
}

void gsDofMapper::markTagged( index_t i, index_t k, index_t comp)
{
    ensureFinalized("markTagged");
    ensureComponent(comp, "markTagged");
    ensurePatch(k, "markTagged");
    ensureLocal(i, k, comp, "markTagged");

    // Tags are stored without the shift, like every other stored dof value,
    // so that markCoupledAsTagged(), permuteFreeDofs(), tindex() and
    // findTagged() -- which all work on stored values -- agree with it, and
    // so that a later setShift() does not detach the tags from their dofs.
    //see gsSortedVector::push_sorted_unique
    index_t t = dofAt(i,k,comp);
    std::vector<index_t>::iterator pos = std::lower_bound(m_tagged.begin(), m_tagged.end(), t );

    if ( pos == m_tagged.end() || *pos != t )// If not found
        m_tagged.insert(pos, t);
}


void gsDofMapper::markBoundary(index_t k, const gsMatrix<index_t> & boundaryDofs, index_t comp)
{
    ensureNotFinalized("markBoundary");
    ensureComponentOrAll(comp, "markBoundary");
    ensurePatch(k, "markBoundary");
    // The dofs are read as boundaryDofs.at(i) for i < rows(): a row vector
    // would lose every entry but its first, and a matrix with rows but no
    // column would be read in storage it does not have.
    GISMO_ENSURE(boundaryDofs.cols() == 1 || boundaryDofs.rows() == 0,
                 "gsDofMapper::markBoundary: Expecting one column of dofs, got "
                 <<boundaryDofs.cols()<<".");
    // Every dof, in every targeted component, before the first is eliminated.
    for (index_t i = 0; i < boundaryDofs.rows(); ++i)
        ensureLocal(boundaryDofs.at(i), k, comp, "markBoundary");

    // Dof by dof, and every component of a dof before the next one: the
    // elimination ids are drawn from one counter shared by all components.
    for (index_t i = 0; i < boundaryDofs.rows(); ++i)
    {
        if (-1 == comp)
            for (index_t c = 0; c != numComponents(); ++c)
                eliminateDofImpl( boundaryDofs.at(i), k, c );
        else
            eliminateDofImpl( boundaryDofs.at(i), k, comp );
    }
}

void gsDofMapper::markCoupledAsTagged()
{
    ensureFinalized("markCoupledAsTagged");
    m_tagged.reserve(m_tagged.size()+m_numCpldDofs.back());
    // The coupled dofs of a component occupy the top of that component's own
    // free block, [m_numFreeDofs[c+1]-nc, m_numFreeDofs[c+1]), where nc is
    // the component's own coupled count -- a prefix difference, since
    // m_numCpldDofs is cumulative after finalize().  This is exactly the band
    // is_coupled_index() tests; the eliminated blocks all sit above it.
    for (size_t c = 0; c+1 != m_numCpldDofs.size(); ++c)
    {
        const index_t nc = m_numCpldDofs[c+1] - m_numCpldDofs[c];
        for (index_t i = 0; i != nc; ++i)
            m_tagged.push_back(m_numFreeDofs[c+1] - nc + i);
    }
    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique(m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );
}

void gsDofMapper::eliminateDof( index_t i, index_t k, index_t comp)
{
    ensureNotFinalized("eliminateDof");
    ensureComponentOrAll(comp, "eliminateDof");
    ensurePatch(k, "eliminateDof");
    ensureLocal(i, k, comp, "eliminateDof");

    if (-1==comp)
    {
        for (index_t c = 0; c != numComponents(); ++c)
            eliminateDofImpl(i,k,c);
        return;
    }
    eliminateDofImpl(i,k,comp);
}

void gsDofMapper::eliminateDofImpl( index_t i, index_t k, index_t comp)
{
    const index_t old = dofAt(i,k,comp);
    if (old == 0)       // regular free dof
    {
        --m_numFreeDofs[comp+1];
        dofAt(i,k,comp) = m_curElimId--;
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
    GISMO_ENSURE(m_curElimId<0, "Error in gsDofMapper::finalize() called twice.");
    checkInvariants();

    for (size_t c = 0; c!=m_dofs.size(); ++c)
      {
	finalizeComp(c);

	//off-set
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

    checkInvariants();
}

void gsDofMapper::finalizeComp(const index_t comp)
{
    std::vector<index_t> & dofs = m_dofs[comp];

    // For assigning coupling and eliminated dofs to continuous
    // indices (-1 = unassigned)
    std::vector<index_t> couplingDofs(m_numCpldDofs[comp+1], -1);
    //std::vector<index_t> elimDofs(-m_curElimId - 1 - m_numElimDofs[comp], -1);
    std::map<index_t,index_t> elimDofs;
    // Free dofs start at offset
    index_t curFreeDof = m_numFreeDofs[comp]+m_numElimDofs[comp];
    // Eliminated dofs start after free dofs plus previous components
    index_t curElimDof = m_numFreeDofs[comp+1] + curFreeDof;

    // Coupling dofs start after standard dofs (=num of zeros in dofs)
    index_t curCplDof = std::count(dofs.begin(), dofs.end(), 0);
    // Devise number of coupled dofs (m_numCpldDofs counted the coupling
    // ids handed out up to here)
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
            const index_t id = -(dofType+1); // dofType may be the lowest index_t
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
    os<<" patches: "<< m_nPatches <<"\n";
    os<<" layout: "<< (GlobalIdentity==m_layout ? "global-identity" : "patch-concatenated") <<"\n";
    os<<" distinct component spaces: "<< (m_hasDistinctComponentSpaces ? "yes" : "no") <<"\n";
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

void gsDofMapper::setIdentity(index_t nPatches, size_t nDofs, size_t nComp)
{
    // Checked here as well as in the overload below, because the vector
    // built for it is allocated first.
    GISMO_ENSURE(nComp <= maxComponentCount(),
                 "setIdentity: "<<nComp<<" components exceed the largest representable "
                 "component count ("<<maxComponentCount()<<").");
    setIdentity(nPatches, std::vector<size_t>(nComp, nDofs));
}

void gsDofMapper::setIdentity(index_t nPatches, const std::vector<size_t> & dofsPerComponent)
{
    GISMO_ENSURE(nPatches > 0, "setIdentity: Expected at least one patch, got " << nPatches << ".");
    GISMO_ENSURE(!dofsPerComponent.empty(), "setIdentity: Expected at least one component.");

    const size_t nComp = dofsPerComponent.size();
    const size_t np = static_cast<size_t>(nPatches);
    GISMO_ENSURE(nComp <= maxComponentCount(),
                 "setIdentity: "<<nComp<<" components exceed the largest representable "
                 "component count ("<<maxComponentCount()<<").");

    // Validate before allocating or narrowing anything: each component's
    // total must be representable, and so must their sum, which finalize()
    // accumulates into m_numFreeDofs (see maxDofCount()).
    size_t total = 0;
    for (size_t c = 0; c != nComp; ++c)
    {
        GISMO_ENSURE(dofsPerComponent[c] <= maxDofCount(),
                     "setIdentity: dof count "<<dofsPerComponent[c]<<" of component "<<c
                     <<" exceeds the largest representable dof count ("<<maxDofCount()<<").");
        GISMO_ENSURE(total <= maxDofCount() - dofsPerComponent[c],
                     "setIdentity: cumulative dof count over components exceeds the largest "
                     "representable dof count ("<<maxDofCount()<<") at component "<<c<<".");
        total += dofsPerComponent[c];
    }

    m_curElimId   = -1;
    m_shift = m_bshift = 0;
    // setIdentity() is a full reset of a possibly already-populated mapper,
    // so every derived member has to go -- including the tag list, whose
    // entries are indices of the numbering being discarded here.
    m_tagged.clear();
    m_numFreeDofs.assign(nComp+1,0);
    m_numElimDofs.assign(nComp+1,0);
    m_numCpldDofs.assign(nComp+1,0);
    for (size_t c = 0; c != nComp; ++c)
        m_numFreeDofs[c+1] = static_cast<index_t>(dofsPerComponent[c]);

    m_nPatches = np;
    m_layout = GlobalIdentity;
    m_hasDistinctComponentSpaces = false;

    // Aliased layout: every real patch offset is zero, only the sentinel
    // (index m_nPatches) carries that component's identity total.
    m_offset.assign(nComp * (np+1), 0);
    m_dofs.assign(nComp, std::vector<index_t>());
    for (size_t c = 0; c != nComp; ++c)
    {
        m_offset[c*(np+1) + np] = dofsPerComponent[c];
        m_dofs[c].assign(dofsPerComponent[c], 0);
    }

    m_uniformComponents = computeUniformComponents();
    checkInvariants();
}

void gsDofMapper::permuteFreeDofs(const gsVector<index_t>& permutation, index_t comp)
{
    ensureFinalized("permuteFreeDofs");
    ensureComponent(comp, "permuteFreeDofs");
    // The permutation is component-local: it permutes component comp's own
    // free block among itself.  m_numFreeDofs is a cumulative prefix sum
    // after finalize(), so the block's length is the prefix difference and
    // its first index is m_numFreeDofs[comp], by which the stored values are
    // rebased before they index the permutation.
    const index_t base  = m_numFreeDofs[comp];
    const index_t nFree = m_numFreeDofs[comp+1] - base;
    GISMO_ENSURE(nFree == permutation.size(), "gsDofMapper::permuteFreeDofs: permutation size "
                 <<permutation.size()<<" does not match the number of free dofs "<<nFree
                 <<" of component "<<comp);
    //GISMO_ASSERT(m_tagged.empty(), "you cannot permute the dofVector twice, combine the permutation");

    // Validated in full before anything is rewritten.  A value out of range
    // would store an index outside the component's free block -- or outside
    // the mapper altogether -- and a repeated value would map two free dofs
    // onto one index and leave another index without a preimage.
    {
        std::vector<bool> seen(nFree, false);
        for (index_t p = 0; p != nFree; ++p)
        {
            const index_t val = permutation[p];
            GISMO_ENSURE(val >= 0 && val < nFree, "gsDofMapper::permuteFreeDofs: permutation value "
                         <<val<<" at position "<<p<<" is out of range [0,"<<nFree<<").");
            GISMO_ENSURE(!seen[val], "gsDofMapper::permuteFreeDofs: permutation value "
                         <<val<<" occurs more than once; not a permutation.");
            seen[val] = true;
        }
    }

    //make a copy of the old ordering, easiest way to implement the permutation. Inplace reordering is quite hard.
    std::vector<index_t> dofs = m_dofs[comp];

    for(index_t i=0; i<(index_t)dofs.size();++i)
    {
        // idx is a stored, unshifted value, so it is classified against the
        // unshifted free range; is_free_index() expects a shifted index.
        const index_t idx = dofs[i];
        if(idx < m_curElimId)
        {
            GISMO_ASSERT(idx-base >= 0 && idx-base < nFree,
                         "gsDofMapper::permuteFreeDofs: free index "<<idx
                         <<" is outside component "<<comp<<"'s free block ["
                         <<base<<","<<base+nFree<<").");
            m_dofs[comp][i] = base + permutation[idx-base];
        }
    }

    // Only the tags inside this component's free block are moved by a
    // component-local permutation; every other tag -- another component's,
    // or an eliminated dof's -- is kept as it is, so the list is remapped in
    // place rather than rebuilt from this component's storage.
    for(std::vector<index_t>::iterator t = m_tagged.begin(); t != m_tagged.end(); ++t)
        if (*t >= base && *t < base + nFree)
            *t = base + permutation[*t - base];

    //sort and delete the duplicated ones
    std::sort(m_tagged.begin(),m_tagged.end());
    std::vector<index_t>::iterator it = std::unique (m_tagged.begin(),m_tagged.end());
    m_tagged.resize( std::distance(m_tagged.begin(),it) );

    // The coupled dofs of this component cannot be tracked anymore, so its
    // own coupled count drops to zero.  m_numCpldDofs is cumulative, so
    // every later prefix loses exactly that count.
    const index_t nCpld = m_numCpldDofs[comp+1] - m_numCpldDofs[comp];
    m_numCpldDofs[comp+1] = m_numCpldDofs[comp];
    for(std::vector<index_t>::iterator s=m_numCpldDofs.begin()+comp+2;
	s<m_numCpldDofs.end(); ++s)
      *s -= nCpld;
}


void gsDofMapper::initPatchDofs(const gsVector<index_t> & patchDofSizes, index_t nComp)
{
    GISMO_ENSURE( nComp > 0, "initPatchDofs: Expected at least one component, got " << nComp << ".");

    const size_t nPatches = patchDofSizes.size();
    GISMO_ENSURE( nPatches > 0, "initPatchDofs: Expected at least one patch, got " << nPatches << ".");

    for (size_t k = 0; k != nPatches; ++k)
        GISMO_ENSURE( patchDofSizes[k] >= 0,
                      "initPatchDofs: Negative patch dof size "<<patchDofSizes[k]
                      <<" at patch "<<k<<".");

    // Build the (single, shared-across-components) offset row first, with
    // a cumulative-overflow check, before touching any member state.
    std::vector<size_t> row(nPatches+1, 0);
    for (size_t k = 0; k < nPatches; ++k)
    {
        const size_t sz = static_cast<size_t>(patchDofSizes[k]);
        GISMO_ENSURE( row[k] <= maxDofCount() - sz,
                      "initPatchDofs: cumulative patch dof size exceeds the largest "
                      "representable dof count ("<<maxDofCount()<<") at patch "<<k<<"." );
        row[k+1] = row[k] + sz;
    }

    // The same size is broadcast to every component, and finalize() sums the
    // per-component totals into m_numFreeDofs, so that sum must fit too.
    GISMO_ENSURE( 0 == row.back() ||
                  static_cast<size_t>(nComp) <= maxDofCount() / row.back(),
                  "initPatchDofs: total dof count over "<<nComp<<" components exceeds the "
                  "largest representable dof count ("<<maxDofCount()<<")." );

    m_curElimId   = -1;
    m_shift = m_bshift = 0;
    m_numElimDofs.assign(nComp+1,0);
    m_numCpldDofs.assign(nComp+1,0);

    m_nPatches = nPatches;
    m_layout = PatchConcatenated;
    m_hasDistinctComponentSpaces = false;

    m_offset.assign(static_cast<size_t>(nComp) * (nPatches+1), 0);
    for (index_t c = 0; c != nComp; ++c)
        std::copy(row.begin(), row.end(), m_offset.begin() + static_cast<size_t>(c)*(nPatches+1));

    m_numFreeDofs.assign(nComp+1, static_cast<index_t>(row.back()));
    m_numFreeDofs.front()=0;

    m_dofs.assign(nComp, std::vector<index_t>(row.back(), 0));

    m_uniformComponents = computeUniformComponents();
    checkInvariants();
}

void gsDofMapper::initRaggedPatchDofs(const std::vector<gsVector<index_t> > & patchDofSizes,
                                       bool hasDistinctComponentSpaces)
{
    const size_t nComp = patchDofSizes.size();
    GISMO_ENSURE( nComp > 0, "gsDofMapper: Expected at least one component, got 0.");
    GISMO_ENSURE( nComp <= maxComponentCount(),
                  "gsDofMapper: "<<nComp<<" components exceed the largest representable "
                  "component count ("<<maxComponentCount()<<").");

    const size_t nPatches = static_cast<size_t>(patchDofSizes.front().size());
    GISMO_ENSURE( nPatches > 0, "gsDofMapper: Expected at least one patch, got 0.");

    for (size_t c = 0; c != nComp; ++c)
        GISMO_ENSURE( static_cast<size_t>(patchDofSizes[c].size()) == nPatches,
                      "gsDofMapper: component "<<c<<" reports "<<patchDofSizes[c].size()
                      <<" patches, expected "<<nPatches<<" (every component must share the "
                      "same patch count).");

    for (size_t c = 0; c != nComp; ++c)
        for (size_t k = 0; k != nPatches; ++k)
            GISMO_ENSURE( patchDofSizes[c][k] >= 0,
                          "gsDofMapper: negative patch dof size "<<patchDofSizes[c][k]
                          <<" at component "<<c<<", patch "<<k<<".");

    // Build every component's offset row first, with a cumulative-overflow
    // check, before touching any member state / allocating m_dofs.
    std::vector<std::vector<size_t> > rows(nComp, std::vector<size_t>(nPatches+1, 0));
    size_t total = 0;
    for (size_t c = 0; c != nComp; ++c)
    {
        for (size_t k = 0; k != nPatches; ++k)
        {
            const size_t sz = static_cast<size_t>(patchDofSizes[c][k]);
            GISMO_ENSURE( rows[c][k] <= maxDofCount() - sz,
                          "gsDofMapper: cumulative patch dof size in component "<<c
                          <<" exceeds the largest representable dof count ("<<maxDofCount()
                          <<") at patch "<<k<<"." );
            rows[c][k+1] = rows[c][k] + sz;
        }

        // finalize() accumulates the per-component totals into m_numFreeDofs,
        // so their sum must be representable as well.
        GISMO_ENSURE( total <= maxDofCount() - rows[c].back(),
                      "gsDofMapper: total dof count over components exceeds the largest "
                      "representable dof count ("<<maxDofCount()<<") at component "<<c<<"." );
        total += rows[c].back();
    }

    m_curElimId   = -1;
    m_shift = m_bshift = 0;
    m_numElimDofs.assign(nComp+1,0);
    m_numCpldDofs.assign(nComp+1,0);

    m_nPatches = nPatches;
    m_layout = PatchConcatenated;
    m_hasDistinctComponentSpaces = hasDistinctComponentSpaces;

    m_offset.assign(nComp * (nPatches+1), 0);
    m_dofs.assign(nComp, std::vector<index_t>());
    m_numFreeDofs.assign(nComp+1, 0);
    for (size_t c = 0; c != nComp; ++c)
    {
        std::copy(rows[c].begin(), rows[c].end(), m_offset.begin() + c*(nPatches+1));
        m_dofs[c].assign(rows[c].back(), 0);
        m_numFreeDofs[c+1] = static_cast<index_t>(rows[c].back());
    }

    m_uniformComponents = computeUniformComponents();
    checkInvariants();
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
    ensureFinalized("preImage");
    typedef std::vector<index_t>::const_iterator citer;
    // gl is shifted and the stored values are not, so the search is for the
    // unshifted value.  The returned dofs are patch-local indices, which the
    // global shift has no part in.
    const index_t comp = componentOf(gl);
    const index_t g    = gl - m_shift;
    const std::vector<index_t> & dofs = m_dofs[comp];
    result.clear();
    size_t cur = 0;//local offsetted index

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
    {
        if ( *it == g )
        {
            if (GlobalIdentity == m_layout)
            {
                // Aliased storage: every patch shares the same range;
                // patch 0 is the canonical stored preimage.
                result.push_back( std::make_pair(index_t(0), static_cast<index_t>(cur)) );
            }
            else
            {
                // Get the patch index of "cur" by "un-offsetting"
                const index_t patch = static_cast<index_t>(
                    std::upper_bound(offBegin(comp), offEnd(comp), cur) - offBegin(comp) - 1);

                // Found a patch-dof pair
                result.push_back( std::make_pair(patch, static_cast<index_t>(cur - offAt(comp,patch))) );
            }
        }
    }
}

std::pair<index_t,index_t> gsDofMapper::anyPreImage(const index_t gl) const
{
    ensureFinalized("anyPreImage");
    typedef std::vector<index_t>::const_iterator citer;
    // Same index conventions as preImage().
    const index_t comp = componentOf(gl);
    const index_t g    = gl - m_shift;
    const std::vector<index_t> & dofs = m_dofs[comp];
    size_t cur = 0;//local offsetted index

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
    {
        if ( *it == g )
        {
            if (GlobalIdentity == m_layout)
                return std::make_pair(index_t(0), static_cast<index_t>(cur));

            // Get the patch index of "cur" by "un-offsetting"
            const index_t patch = static_cast<index_t>(
                std::upper_bound(offBegin(comp), offEnd(comp), cur) - offBegin(comp) - 1);

            // Found a patch-dof pair
            return std::make_pair(patch, static_cast<index_t>(cur - offAt(comp,patch)));
        }
    }
    GISMO_ERROR("The global index "<< gl <<" is not valid");
}

std::vector<std::pair<index_t,index_t> > gsDofMapper::anyPreImages(index_t comp) const
{
    // result[*it] is indexed by stored values, which are indices in
    // [0,size()) only once finalized.
    ensureFinalized("anyPreImages");
    ensureComponent(comp, "anyPreImages");
    typedef std::vector<index_t>::const_iterator citer;
    const std::vector<index_t> & dofs = m_dofs[comp];
    size_t cur = 0;//local offsetted index

    // One entry per global index, at the unshifted position: the stored
    // values are exactly the unshifted indices, all of them below size().
    // A dof that this component does not own gets a sentinel in both slots,
    // since 0 is a valid patch-local index.
    std::vector<std::pair<index_t,index_t> > result(size(), std::make_pair(index_t(-1),index_t(-1)));

    for (citer it = dofs.begin(); it != dofs.end(); ++it, ++cur)
    {
        if ( -1 == result[*it].first )
        {
            if (GlobalIdentity == m_layout)
            {
                result[*it] = std::make_pair(index_t(0), static_cast<index_t>(cur));
            }
            else
            {
                // Get the patch index of "cur" by "un-offsetting"
                const index_t patch = static_cast<index_t>(
                    std::upper_bound(offBegin(comp), offEnd(comp), cur) - offBegin(comp) - 1);

                // Found a patch-dof pair
                result[*it] = std::make_pair(patch, static_cast<index_t>(cur - offAt(comp,patch)));
            }
        }
    }
    return result;
}

gsVector<index_t> gsDofMapper::inverseAsVector(index_t comp) const
{
    ensureFinalized("inverseAsVector");
    ensureComponent(comp, "inverseAsVector");
    GISMO_ASSERT(isPermutation(), "This dofMapper is not 1-1");
    // Every position that is not the image of a local dof of this component
    // -- every other component's block, in particular -- holds -1.
    gsVector<index_t> v = gsVector<index_t>::Constant(size(), -1);
      for(size_t j = 0; j!= m_dofs[comp].size(); ++j)
	v[ m_dofs[comp][j] ] = j;
    return v;
}

std::map<index_t,index_t>
gsDofMapper::inverseOnPatch(const index_t k) const
{
    ensureFinalized("inverseOnPatch");
    ensurePatch(k, "inverseOnPatch");

    std::map<index_t,index_t> inv;
    //inv.reserve(patchSize(k));
    typedef std::vector<index_t>::const_iterator citer;

    for(size_t i = 0; i!= m_dofs.size(); ++i)
    {
        const index_t c = static_cast<index_t>(i);
        // Only the dofs that live on patch k, as many as this component has
        // there.  Under the aliased layout patchSize() is the component's
        // global total on every patch, so the same expression yields the
        // complete inverse.
        //
        // The keys are global indices, so they carry the shift like index().
        citer it = m_dofs[i].begin() + offAt(c, k);
        const size_t n = patchSize(k, c);
        for(size_t j = 0; j!= n; ++j,++it)
            inv[*it + m_shift]=j;
    }
    return inv;
}

bool gsDofMapper::indexOnPatch(const index_t gl, const index_t k, index_t & local) const
{
    ensureFinalized("indexOnPatch");
    ensurePatch(k, "indexOnPatch");
    typedef std::vector<index_t>::const_iterator citer;
    // gl is shifted, the stored values are not.  An index outside the
    // mapper's range lives on no patch at all.
    if (!shiftedInRange(gl, size())) return false;
    const index_t g = gl - m_shift;
    const index_t comp = componentOfUnshifted(g);
    const std::vector<index_t> & dofs = m_dofs[comp];
    const citer istart = dofs.begin() + offAt(comp, k);
    const citer iend   = istart + patchSize(k, comp);
    auto it = std::find(istart, iend, g);
    if (iend==it) return false;
    local = std::distance(istart,it);
    return true;
}

index_t gsDofMapper::boundarySizeWithDuplicates() const
{
    ensureFinalized("boundarySizeWithDuplicates");

    // Eliminated dofs of every component are numbered above the free dofs of
    // ALL components, so the threshold is the global free count.
    const index_t s = m_numFreeDofs.back() - 1;
    index_t res = 0;
    for (size_t i = 0; i!= m_dofs.size(); ++i)
      res += std::count_if(m_dofs[i].begin(), m_dofs[i].end(),
			   GS_BIND2ND(std::greater<index_t>(), s) );
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

gsVector<index_t> gsDofMapper::findBoundary(const index_t k, const index_t comp) const
{
    ensureFinalized("findBoundary");
    ensurePatch(k, "findBoundary");
    ensureComponent(comp, "findBoundary");
    // Eliminated dofs of every component are numbered above the free dofs of
    // all components, so this threshold is component-independent.  Like
    // every threshold in the find* queries it is compared with the stored,
    // unshifted values, so the shift has no part in it.
    const index_t s = m_numFreeDofs.back() - 1;
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[comp].begin() + offAt(comp,k);
    citer iend   = istart + patchSize(k,comp);
    return find_impl(istart, iend, GS_BIND2ND(std::greater<index_t>(),s));
}

gsVector<index_t> gsDofMapper::findFree(const index_t k, const index_t comp) const
{
    ensureFinalized("findFree");
    ensurePatch(k, "findFree");
    ensureComponent(comp, "findFree");
    const index_t s = m_numFreeDofs.back();
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[comp].begin() + offAt(comp,k);
    citer iend   = istart + patchSize(k,comp);
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

gsVector<index_t> gsDofMapper::findCoupled(const index_t k, const index_t j,
                                           const index_t comp) const
{
    ensureFinalized("findCoupled");
    ensurePatch(k, "findCoupled");
    ensureComponent(comp, "findCoupled");

    GISMO_ENSURE(-1 == j || validPatch(j), "gsDofMapper::findCoupled: invalid patch "<<j
                 <<", expected -1 (any) or a value in [0,"<<m_nPatches<<").");
    if (k==j) return gsVector<index_t>();

    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[comp].begin() + offAt(comp,k);
    citer iend   = istart + patchSize(k,comp);

    // The coupled dofs of a component sit at the top of that component's own
    // free block; both bounds are prefix values of that component, not the
    // last component's totals.
    const index_t nCpld = m_numCpldDofs[comp+1] - m_numCpldDofs[comp];
    const index_t l = m_numFreeDofs[comp+1]-nCpld-1;
    const index_t u = m_numFreeDofs[comp+1];
    if (-1==j)
        return find_impl(istart, iend, _isBetween(l,u) );
    else
    {
        citer istartj = m_dofs[comp].begin() + offAt(comp,j);
        citer iendj   = istartj + patchSize(j,comp);
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

gsVector<index_t> gsDofMapper::findFreeUncoupled(const index_t k, const index_t comp) const
{
    ensureFinalized("findFreeUncoupled");
    ensurePatch(k, "findFreeUncoupled");
    ensureComponent(comp, "findFreeUncoupled");
    typedef std::vector<index_t>::const_iterator citer;
    const citer istart = m_dofs[comp].begin() + offAt(comp,k);
    const citer iend   = istart + patchSize(k,comp);
    // Below this component's own coupled band and at or above the start of
    // its own free block.
    const index_t nCpld = m_numCpldDofs[comp+1] - m_numCpldDofs[comp];
    return find_impl(istart, iend,
                     _isBetween(m_numFreeDofs[comp]-1,
                                m_numFreeDofs[comp+1]-nCpld) );
}

gsVector<index_t> gsDofMapper::findTagged(const index_t k, const index_t comp) const
{
    ensureFinalized("findTagged");
    ensurePatch(k, "findTagged");
    ensureComponent(comp, "findTagged");
    typedef std::vector<index_t>::const_iterator citer;
    citer istart = m_dofs[comp].begin() + offAt(comp,k);
    citer iend   = istart + patchSize(k,comp);
    std::list<index_t> si;
    std::set_intersection(istart, iend, m_tagged.begin(),
                          m_tagged.end(), std::back_inserter(si));
    gsVector<index_t> rvo;
    si.assign(si.begin(),si.end());
    return rvo;
}

void gsDofMapper::componentFailed(index_t c, const char * where) const
{
    GISMO_ENSURE(validComponent(c), "gsDofMapper::"<<where<<": invalid component "<<c
                 <<", the mapper has "<<m_dofs.size()<<" components.");
}

void gsDofMapper::componentOrAllFailed(index_t c, const char * where) const
{
    GISMO_ENSURE(-1 == c || validComponent(c), "gsDofMapper::"<<where<<": invalid component "<<c
                 <<", expected -1 (all) or a value in [0,"<<m_dofs.size()<<").");
}

void gsDofMapper::patchFailed(index_t k, const char * where) const
{
    GISMO_ENSURE(validPatch(k), "gsDofMapper::"<<where<<": invalid patch "<<k
                 <<", the mapper has "<<m_nPatches<<" patches.");
}

void gsDofMapper::localFailed(index_t i, index_t k, index_t c, const char * where) const
{
    GISMO_ENSURE(validLocal(i,k,c), "gsDofMapper::"<<where<<": invalid local dof "<<i
                 <<" of patch "<<k<<", component "<<c<<localRangeInfo(k,c));
}

void gsDofMapper::finalizedFailed(const char * where) const
{
    GISMO_ENSURE(m_curElimId>=0, "gsDofMapper::"<<where<<": finalize() was not called.");
}

void gsDofMapper::notFinalizedFailed(const char * where) const
{
    GISMO_ENSURE(m_curElimId<0, "gsDofMapper::"<<where<<": the mapper is already finalized.");
}

void gsDofMapper::ensureShift(index_t shift, bool boundary, const char * where) const
{
    // Stored indices lie in [0,n) with n = size() (global) or boundarySize()
    // (boundary numbering) once finalized, and n <= mapSize() before that,
    // since every index has at least one preimage.  mapSize() itself was
    // checked against index_t at construction.  The shifted range, including
    // its one-past-the-end value, must fit: shift + n <= max.  A negative
    // shift only moves the range down from shift >= min, which always fits.
    const index_t n = m_curElimId >= 0 ? (boundary ? boundarySize() : size())
                                       : static_cast<index_t>(mapSize());
    GISMO_ENSURE(shift <= std::numeric_limits<index_t>::max() - n,
                 "gsDofMapper::"<<where<<": shift "<<shift<<" would take the "<<n<<" indices of this "
                 "mapper past the largest representable index ("<<std::numeric_limits<index_t>::max()<<").");
}

void gsDofMapper::setShift (index_t shift)
{
    ensureShift(shift, false, "setShift");
    m_shift=shift;
}

void gsDofMapper::addShift (index_t shift)
{
    GISMO_ENSURE(shift >= 0 ? m_shift <= std::numeric_limits<index_t>::max() - shift
                            : m_shift >= std::numeric_limits<index_t>::min() - shift,
                 "gsDofMapper::addShift: adding "<<shift<<" to the shift "<<m_shift<<" overflows.");
    ensureShift(m_shift + shift, false, "addShift");
    m_shift+=shift;
}

void gsDofMapper::setBoundaryShift (index_t shift)
{
    ensureShift(shift, true, "setBoundaryShift");
    m_bshift=shift;
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
    // Following requires binding of gsBoxTopology
    // .def("setMatchingInterfaces",       &Class::setMatchingInterfaces,      "Called to initialize the gsDofMapper with matching interfaces after m_bases have already been set")
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

    // gsDofMapper::setIdentity is overloaded (scalar total, and
    // per-component totals): explicit casts are required here so that
    // the correct overload is bound -- taking &Class::setIdentity
    // directly is ambiguous once there is more than one overload.
    .def("setIdentity", static_cast<void (Class::*)(index_t,size_t,size_t)>(&Class::setIdentity), "Set this mapping to be the identity")
    .def("setIdentity", static_cast<void (Class::*)(index_t,const std::vector<size_t>&)>(&Class::setIdentity), "Set this mapping to be the identity, with a dof total per component")
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

    // pybind11 does not inherit C++ default arguments: offset() gained a
    // defaulted component parameter, so the default has to be restated here
    // or every existing one-argument Python call breaks.
    .def("offset", &Class::offset, "Returns the offset corresponding to patch \a k for component \a c",
         py::arg("k"), py::arg("c") = 0)
    .def("numPatches", &Class::numPatches, "Returns the number of patches present underneath the mapper")
    .def("mapSize", &Class::mapSize, "Returns the total number of patch-local degrees of freedom that are being mapped")
    .def("componentsSize", &Class::componentsSize, "Returns the components size")
    .def("hasDistinctComponentSpaces", &Class::hasDistinctComponentSpaces, "Returns whether this mapper was declared to be built from distinct per-component bases")
    .def("hasUniformComponents", &Class::hasUniformComponents, "Returns whether every component can be treated as one shared space")
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
