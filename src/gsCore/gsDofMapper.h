/** @file gsDofMapper.h

    @brief Provides the gsDofMapper class for re-indexing DoFs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsExport.h>

namespace gismo
{

/** @brief Maintains a mapping from patch-local dofs to global dof indices
    and allows the elimination of individual dofs.

    A \em dof (degree of freedom) is, roughly speaking, an unknown in a
    discretization of a PDE. However, some dofs may be eliminated before
    solving the system and won't actually translate into unknowns.
    An example are dofs on Dirichlet boundaries.

    This class creates a mapping between an arbitrary number of
    per-patch local dofs to an enumeration of global dofs.
    Furthermore, dofs can also be marked as eliminated.

    This is a a many-to-one mapping: many patch-local dofs are mapped
    to a single global dfo (index).
    Every global dof gets a unique number, forming a continuous range
    starting from 0. This range has length gsDofMapper::size().
    The dofs are numbered in the following order:

    - first the standard \em free (non-eliminated) dofs, ie. dofs that
      are not coupled on the boundary. For the standard dofs there is
      unique pre-image pair (patch,localdof).

    - then the \em free dofs which are coupled with other dofs For the
      coupled dofs there is a list of pre-image pairs of the form
      (patch,localdof).

      Upto here we get all the dofs which are \em free (number:
      gsDofMapper::freeSize() ). Then a final group follows:

    - then the dofs that are on Dirichlet boundaries (number:
      gsDofMapper::boundarySize()). These dofs might have a unique
      pre-image or not.

    The boundary (eg. eliminated) dofs have their own 0-based
    numbering. The index of an boundary global dof in this numbering
    can be queried with gsDofMapper::bindex().

    The object must be finalized before it is used,
    i.e. gsDofMapper::finalize() has to be called once before use.

    Each component of the mapper carries its own patch-local storage
    (see initPatchDofs()/the ragged constructor), so the number of
    local dofs on a given patch may differ from one component to
    another ("ragged" storage).  There are two kinds of layout
    (see layout()):

    - patch-concatenated: ordinary storage, one contiguous range of
      local dofs per (patch,component);
    - global identity / aliased: produced by setIdentity(), where the
      local index is already a component-global index and does not
      depend on the patch argument.  See the per-method documentation
      below for the exact aliased semantics.

    Whether a mapper was built from more than one distinct
    per-component basis object is a separately \em declared property
    (hasDistinctComponentSpaces()), never inferred from sizes or
    offsets: two components can have identical per-patch sizes while
    still being built from different bases (this is precisely the
    Raviart-Thomas situation on an isotropic mesh).

    \ingroup Core

*/
class GISMO_EXPORT gsDofMapper
{
public:

    /// The two ways a mapper's per-(patch,component) storage may be
    /// laid out.  A mapper never infers which one it has from a
    /// coincidental offset pattern -- the layout is fixed at
    /// construction and carried through reset/swap.
    enum Layout
    {
        PatchConcatenated = 0, ///< ordinary patch-local storage (the default)
        GlobalIdentity    = 1  ///< setIdentity()-built aliased/global storage
    };

    /// Default empty constructor
    gsDofMapper();

    /**
     * @brief construct a dof mapper with a given number of dofs per
     * patch
     *
     * @param patchDofSizes
     */
    gsDofMapper(const gsVector<index_t> &patchDofSizes, index_t nComp = 1)
    {
        initPatchDofs(patchDofSizes, nComp);
    }

    /**
     * @brief Construct a mapper from one patch-dof-size vector per
     * component (patch-concatenated layout).  Every component must
     * report the same number of patches; local sizes may otherwise
     * differ freely between components and between patches ("ragged"
     * storage).
     *
     * @param patchDofSizes           one entry per component, each a
     *                                vector of per-patch local dof counts
     * @param hasDistinctComponentSpaces
     *                                must be declared explicitly by the
     *                                caller, never inferred: true when
     *                                the components were built from
     *                                genuinely different basis objects
     *                                (e.g. a Raviart-Thomas component
     *                                pair), false when they happen to
     *                                share one basis.  See the class
     *                                documentation: a square-mesh
     *                                Raviart-Thomas mapper has identical
     *                                per-component sizes even though its
     *                                bases differ, so this cannot be
     *                                derived from \a patchDofSizes.
     */
    gsDofMapper(const std::vector<gsVector<index_t> > & patchDofSizes,
                bool hasDistinctComponentSpaces);

    void swap(gsDofMapper & other)
    {
        m_dofs  .swap(other.m_dofs);
        m_offset.swap(other.m_offset);
        std::swap(m_nPatches, other.m_nPatches);
        std::swap(m_layout,   other.m_layout);
        std::swap(m_hasDistinctComponentSpaces, other.m_hasDistinctComponentSpaces);
        std::swap(m_uniformComponents, other.m_uniformComponents);

        std::swap(m_shift      , other.m_shift);
        std::swap(m_bshift     , other.m_bshift);
        std::swap(m_numFreeDofs, other.m_numFreeDofs);
        std::swap(m_numElimDofs, other.m_numElimDofs);
        std::swap(m_numCpldDofs, other.m_numCpldDofs);
        std::swap(m_curElimId  , other.m_curElimId);
        std::swap(m_tagged     , other.m_tagged);
    }

private:

    /// Initialize by vector of DoF indices and dimension
    void initPatchDofs(const gsVector<index_t> & patchDofSizes,
		       index_t nComp = 1);

    /// Initialize by one patch-dof-size vector per component (ragged,
    /// patch-concatenated layout).
    void initRaggedPatchDofs(const std::vector<gsVector<index_t> > & patchDofSizes,
                              bool hasDistinctComponentSpaces);

    // Flat-offset-table accessors (the only way this class touches
    // m_offset): m_offset is one row of m_nPatches+1 entries per
    // component, laid out consecutively.
    size_t offAt(index_t c, index_t k) const
    { return m_offset[static_cast<size_t>(c)*(m_nPatches+1)+static_cast<size_t>(k)]; }

    std::vector<size_t>::const_iterator offBegin(index_t c) const
    { return m_offset.begin() + static_cast<size_t>(c)*(m_nPatches+1); }

    std::vector<size_t>::const_iterator offEnd(index_t c) const
    { return offBegin(c) + (m_nPatches+1); }

    /// Read-write access to the stored value of local dof \a i of patch
    /// \a k in component \a c (setup-time encoding: 0 = free, negative =
    /// eliminated, positive = coupling id -- see the m_dofs comment below).
    /// The only way this class' own code touches m_dofs/m_offset together;
    /// overloaded on constness instead of macro-expanded so both mutating
    /// setup code and const query methods can use one accessor.
    ///
    /// The bounds check is debug-only.  dofAt() runs once per local dof per
    /// element in every assembler, so the per-dof accessors built on it
    /// (index(), bindex(), cindex(), tindex(), freeIndex(), and the
    /// is_free/is_boundary/is_coupled/is_tagged queries taking a local dof)
    /// and the per-element localToGlobal()/localToGlobal2() are not checked
    /// in Release builds; every other public entry point validates its
    /// arguments with the ensure*() helpers below before reaching it.
    inline index_t & dofAt(index_t i, index_t k, index_t c)
    {
        GISMO_ASSERT(validLocal(i,k,c), "gsDofMapper: invalid local dof "<<i<<" of patch "<<k
                     <<", component "<<c<<localRangeInfo(k,c));
        return m_dofs[c][offAt(c,k)+i];
    }

    inline index_t dofAt(index_t i, index_t k, index_t c) const
    {
        GISMO_ASSERT(validLocal(i,k,c), "gsDofMapper: invalid local dof "<<i<<" of patch "<<k
                     <<", component "<<c<<localRangeInfo(k,c));
        return m_dofs[c][offAt(c,k)+i];
    }

    // --- argument validation -------------------------------------------
    //
    // The valid*() predicates never touch storage for an invalid argument.
    // The ensure*() helpers throw std::runtime_error (GISMO_ENSURE) in
    // Release builds too: they guard memory safety, so they must not
    // compile away.  \a where names the public method for the message.

    bool validComponent(index_t c) const
    { return c >= 0 && static_cast<size_t>(c) < m_dofs.size(); }

    bool validPatch(index_t k) const
    { return k >= 0 && static_cast<size_t>(k) < m_nPatches; }

    /// Number of valid local indices of (patch \a k, component \a c), for
    /// arguments already known to be valid.  Patch-specific under the
    /// patch-concatenated layout, so that an oversized local index cannot
    /// spill into the next patch's storage; the component's global identity
    /// total on every patch under the aliased layout.
    size_t localCount(index_t k, index_t c) const
    {
        return GlobalIdentity == m_layout ? m_dofs[c].size()
                                          : offAt(c,k+1) - offAt(c,k);
    }

    bool validLocal(index_t i, index_t k, index_t c) const
    {
        return validComponent(c) && validPatch(k) &&
            i >= 0 && static_cast<size_t>(i) < localCount(k,c);
    }

    /// Message tail describing the valid local range of (\a k, \a c), or
    /// why there is none.
    std::string localRangeInfo(index_t k, index_t c) const
    {
        std::ostringstream os;
        if (!validComponent(c))
            os << ": the mapper has " << m_dofs.size() << " components.";
        else if (!validPatch(k))
            os << ": the mapper has " << m_nPatches << " patches.";
        else
            os << ": the valid range is [0," << localCount(k,c) << ").";
        return os.str();
    }

    // Each ensure*() below is an inline O(1) test; the message and the throw
    // live in the matching *Failed() function in gsDofMapper.cpp, called
    // only for an invalid argument.  Formatting the message inline makes
    // the helper too large to inline, so every checked query, including
    // the inline size()/patchSize()/offset() family, would pay an
    // out-of-line call for a check that is one comparison.

    void ensureComponent(index_t c, const char * where) const
    { if (!validComponent(c)) componentFailed(c, where); }

    /// As ensureComponent(), but also accepts the broadcast value -1.
    void ensureComponentOrAll(index_t c, const char * where) const
    { if (-1 != c && !validComponent(c)) componentOrAllFailed(c, where); }

    void ensurePatch(index_t k, const char * where) const
    { if (!validPatch(k)) patchFailed(k, where); }

    /// Local dof \a i of patch \a k must exist in component \a c, or, for
    /// \a c == -1, in every component: a broadcast is validated for all its
    /// targets before any of them is touched.
    void ensureLocal(index_t i, index_t k, index_t c, const char * where) const
    {
        if (-1 == c)
        {
            for (index_t cc = 0; static_cast<size_t>(cc) != m_dofs.size(); ++cc)
                ensureLocal(i, k, cc, where);
            return;
        }
        if (!validLocal(i,k,c)) localFailed(i, k, c, where);
    }

    void ensureFinalized(const char * where) const
    { if (m_curElimId<0) finalizedFailed(where); }

    /// Setup mutators (matchDof(s), markCoupled, eliminateDof, markBoundary,
    /// colapseDofs) rewrite the setup-time encoding of m_dofs, which
    /// finalize() replaces by the final numbering: applied afterwards they
    /// would misread global indices as coupling or elimination ids.
    void ensureNotFinalized(const char * where) const
    { if (m_curElimId>=0) notFinalizedFailed(where); }

    // Throw for the argument the matching ensure*() rejected.
    void componentFailed(index_t c, const char * where) const;
    void componentOrAllFailed(index_t c, const char * where) const;
    void patchFailed(index_t k, const char * where) const;
    void localFailed(index_t i, index_t k, index_t c, const char * where) const;
    void finalizedFailed(const char * where) const;
    void notFinalizedFailed(const char * where) const;

    /// Validates a new global (\a boundary false) or boundary (\a boundary
    /// true) shift: every index this mapper hands out, and one past the last
    /// of them (lastIndex(), firstIndex()+freeSize()), must be representable.
    /// Before finalize() the final size is unknown; mapSize() bounds it.
    void ensureShift(index_t shift, bool boundary, const char * where) const;

    /// Single-component setup mutators, for arguments already validated.
    void matchDofImpl(index_t u, index_t i, index_t v, index_t j, index_t comp);
    void eliminateDofImpl(index_t i, index_t k, index_t comp);

    /// True iff the shifted global index \a gl lies in [m_shift, m_shift+n).
    /// Written so that nothing overflows for any \a gl and any shift, which
    /// the obvious "g = gl - m_shift; 0 <= g && g < n" does not: that
    /// subtraction is signed overflow for \a gl near the bottom of index_t
    /// with a positive shift, or near the top with a negative one.  Once
    /// gl >= m_shift the true difference lies in [0, 2^digits), which the
    /// unsigned subtraction represents exactly.  After it returns true,
    /// gl - m_shift is in [0,n) and safe to compute.  n < 0 (an unfinalized
    /// m_curElimId) is an empty range.
    bool shiftedInRange(index_t gl, index_t n) const
    {
        typedef std::make_unsigned<index_t>::type uindex_t;
        return n > 0 && gl >= m_shift &&
            static_cast<uindex_t>(gl) - static_cast<uindex_t>(m_shift) < static_cast<uindex_t>(n);
    }

    /// Component of the UNSHIFTED index \a g, i.e. of a stored dof value,
    /// which must lie in [0,size()).  After finalize() the free blocks of
    /// all components come first and their eliminated blocks follow, so
    /// \a g is looked up in whichever of the two prefix vectors covers it.
    index_t componentOfUnshifted(index_t g) const
    {
        return (g<m_numFreeDofs.back() ?
                std::distance(m_numFreeDofs.begin(), std::upper_bound(m_numFreeDofs.begin(), m_numFreeDofs.end(), g))
              : std::distance(m_numElimDofs.begin(),std::upper_bound(m_numElimDofs.begin(), m_numElimDofs.end(), g-m_numFreeDofs.back())) ) - 1;
    }

    /// The value hasUniformComponents() caches, computed from the sizes and
    /// the declared distinctness in O(nPatches*nComp).
    bool computeUniformComponents() const;

    /// Debug-only structural invariant check (GISMO_ASSERT-based, so it
    /// compiles away entirely under NDEBUG -- this is not a release-mode
    /// guard).  Invoked after construction/reset and before/after
    /// finalize().
    void checkInvariants() const;

public:

    /// Returns a vector taking flat local indices to global
    gsVector<index_t> asVector(index_t comp = 0) const;

    /** \brief Returns a vector taking global indices to flat local

        The inverse of asVector(\a comp) in the full (unshifted) global
        index space: size() entries, with -1 in every position that is not
        the image of a local dof of component \a comp.  The entry of global
        index \a gl is at position \a gl minus the shift, as in
        anyPreImages().

        Assumes that the mapper is a permutation
    */
    gsVector<index_t> inverseAsVector(index_t comp = 0) const;

    /// Called to initialize the gsDofMapper with matching interfaces
    /// after m_bases have already been set
    void setMatchingInterfaces(const gsBoxTopology & mp);

    // The setup mutators below (colapseDofs, matchDof, matchDofs,
    // markCoupled, markBoundary, eliminateDof) may only be called before
    // finalize().  Their component argument is -1 (every component) or a
    // valid component, and each local dof must lie in [0, patchSize(k,c))
    // for every component it targets.  All arguments of a call -- every
    // entry of a batch, every component of a broadcast -- are validated
    // before the first change, so a call that throws leaves the mapper as
    // it was.  The checks throw in Release builds too.

    /** \brief Calls matchDof() for all dofs on the given patch side
     * \a i ps. Thus, the whole set of dofs collapses to a single
     * global dof
     *
     */
    void colapseDofs(index_t k, const gsMatrix<unsigned> & b, index_t comp = 0);

    /// \brief Couples dof \a i of patch \a u with dof \a j of patch
    /// \a v such that they refer to the same global dof at component
    /// \a comp.
    void matchDof( index_t u, index_t i, index_t v, index_t j, index_t comp = 0);

    /// \brief Couples dofs \a b1 of patch \a u with dofs \a b2 of patch
    /// \a v one by one such that they refer to the same global dof.
    void matchDofs(index_t u, const gsMatrix<index_t> & b1,
                   index_t v, const gsMatrix<index_t> & b2,
		           index_t comp = 0);

    /// Mark the local dof \a i of patch \a k as coupled.
    void markCoupled(index_t i, index_t k, index_t comp = 0);

    /// Mark a local dof \a i of patch \a k as tagged
    void markTagged(index_t i, index_t k, index_t comp = 0);

    /// Mark all coupled dofs as tagged
    void markCoupledAsTagged();

    /// Mark the local dofs \a boundaryDofs of patch \a k as eliminated.
    // to do: put k at the end
    void markBoundary(index_t k, const gsMatrix<index_t> & boundaryDofs, index_t comp = 0);

    /// Mark the local dof \a i of patch \a k as eliminated.
    void eliminateDof(index_t i, index_t k, index_t comp = 0);

    /// \brief Must be called after all boundaries and interfaces have
    /// been marked to set up the dof numbering.  Must be called exactly
    /// once; a second call throws.
    void finalize();

    /// \brief Checks whether finalize() has been called.
    bool isFinalized() const { return m_curElimId>=0; }

    /// \brief Returns true iff the mapper is a permuatation
    bool isPermutation() const { return static_cast<size_t>(size())==mapSize(); }

    /// \brief Print summary
    std::ostream& print( std::ostream& os = gsInfo ) const;

    ///\brief Set this mapping to be the identity
    void setIdentity(index_t nPatches, size_t nDofs, size_t nComp = 1);

    ///\brief Set this mapping to be the identity, with a (possibly
    /// unequal) total dof count per component.  The scalar-size overload
    /// delegates here by broadcasting one total to every component.
    void setIdentity(index_t nPatches, const std::vector<size_t> & dofsPerComponent);

    ///\brief Set the shift amount for the global numbering
    ///
    /// Throws if the shifted indices, up to and including one past the
    /// last, would not be representable as index_t.  Before finalize() the
    /// bound is taken from mapSize(), which no final size can exceed.
    void setShift(index_t shift);

    ///\brief Add a shift amount to the global numbering (same bound as
    /// setShift()).
    void addShift(index_t shift);

    /// \brief Permutes the mapped free indices according to permutation, i.e.,  dofs_perm[idx] = dofs_old[permutation[idx]]
    ///
    /// \a permutation must be a permutation of [0, n) with n the number of
    /// free dofs of component \a comp; anything else throws before the
    /// mapper is changed.
    ///
    /// \warning Applying a permutation makes the functions regarding coupled dofs (cindex, is_coupled_index,.. ) invalid.
    /// The dofs are still coupled, but you have no way of extracting them. If you need this functions, first call
    /// markCoupledAsTagged() and then use the corresponding functions for tagged dofs.
    void permuteFreeDofs(const gsVector<index_t>& permutation, index_t comp = 0);

    ///\brief Returns the smallest global index present in component \a comp
    ///
    /// After finalize() the free blocks of all components come first,
    /// component by component, and the eliminated blocks of all components
    /// follow them.  There are three cases:
    ///
    /// - the component has free dofs: the start of its own free block;
    /// - it has none but has eliminated ones: the start of its own
    ///   eliminated block, which lies above every component's free block;
    /// - it owns no dof at all (an empty component, or a mapper that has
    ///   not been finalized): no index is present, and the start of its
    ///   free block is reported.
    ///
    /// \a comp may also equal numComponents(), which reports the end of the
    /// last component's free block; in particular firstIndex() is valid on
    /// a default-constructed mapper.
    index_t firstIndex(index_t comp = 0) const
    {
        GISMO_ENSURE(comp >= 0 && static_cast<size_t>(comp) < m_numFreeDofs.size(),
                     "gsDofMapper::firstIndex: invalid component "<<comp<<", expected a value in [0,"
                     <<m_numFreeDofs.size()-1<<"].");
        if (static_cast<size_t>(comp)+1 < m_numFreeDofs.size() &&
            m_numFreeDofs[comp+1] == m_numFreeDofs[comp] &&   // no free dof
            m_numElimDofs[comp+1] != m_numElimDofs[comp])     // but eliminated ones
            return m_numFreeDofs.back() + m_numElimDofs[comp] + m_shift;
        return m_numFreeDofs[comp] + m_shift;
    }

    ///\brief Returns one past the biggest value of the free indices
    index_t lastIndex() const { return m_shift + freeSize(); }

    ///\brief Set the shift amount for the boundary numbering (same
    /// bound as setShift(), with boundarySize() in place of size()).
    void setBoundaryShift(index_t shift);

    /** \brief Computes the global indices of the input local indices
     *
     * \param[in] locals a column matrix with the local indices
     * \param[in] patchIndex the index of the patch where the local indices belong to
     * \param[out] globals the global indices of the patch
     *
     * \note On the assembly hot path: like index(), its arguments are
     * checked in debug builds only.
     */
    void localToGlobal(const gsMatrix<index_t>& locals,
                       index_t patchIndex,
                       gsMatrix<index_t>& globals,
		               index_t comp = 0) const;

    /** \brief Computes the global indices of the input local indices
     *
     * \param[in] locals a column matrix with the local indices
     * \param[in] patchIndex the index of the patch where the local indices belong to
     * \param[out] globals the local-global correspondance
     * \param[out] numFree the number of free indices in \a local
     *
     * \note On the assembly hot path: like index(), its arguments are
     * checked in debug builds only.  \a locals and \a globals must be
     * distinct objects.
     */
    void localToGlobal2(const gsMatrix<index_t>& locals,
                        index_t patchIndex,
                        gsMatrix<index_t>& globals,
                        index_t & numFree,
		                index_t comp = 0) const;

    /** \brief Returns the index associated to local dof \a i of patch \a k without shifts.
     *
     * \note This method only works after all interfaces and boundaries have been
     * marked and finalize() has been called.
     */

    inline index_t freeIndex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return dofAt(i,k,c);
    }

    /// \brief Returns the component that the global dof \a gl belongs to.
    ///
    /// \a gl is a global index as returned by index(), i.e. including the
    /// shift \a s set by setShift(), and must lie in [s, s+size());
    /// anything else throws.
    index_t componentOf(index_t gl) const
    {
         GISMO_ASSERT(m_curElimId>=0,"finalize() was not called on gsDofMapper");
         GISMO_ENSURE(shiftedInRange(gl, size()), "gsDofMapper::componentOf(): global index "<<gl
                      <<" is outside the mapper's range (shift "<<m_shift<<", size "<<size()<<")");
         return componentOfUnshifted(gl - m_shift);
    }

    /** \brief Returns the global dof index associated to local dof \a i of patch \a k.
     *
     * \note This method only works after all interfaces and boundaries have been
     * marked and finalize() has been called.
     *
     * \note Like the other per-dof accessors (bindex(), cindex(), tindex(),
     * freeIndex()) this checks its arguments in debug builds only: it is on
     * the assembly hot path.  \a i must lie in [0, patchSize(k,c)).
     */

    inline index_t index(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return dofAt(i,k,c)+m_shift;
    }

    /// @brief Returns the boundary index of local dof \a i of patch \a k.
    ///
    /// Produces undefined results if local dof (i,k) does not lie on the boundary.
    inline index_t bindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return dofAt(i,k,c) - m_numFreeDofs.back()
            //- m_numElimDofs[c]
            + m_bshift;
    }

    /// Returns true iff all DoFs are considered as free
    bool allFree() const
    { return m_numFreeDofs.back()+m_numElimDofs.back()==m_curElimId; }

    /// @brief Returns the coupled dof index
    ///
    /// The coupled dofs of component \a c occupy the top of that component's
    /// own free block, and are numbered consecutively across components, so
    /// the shift applied here is the CUMULATIVE coupled prefix -- unlike
    /// is_coupled_index(), which needs the component's own count to locate
    /// the band.
    inline index_t cindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return dofAt(i,k,c) - m_numFreeDofs[c+1]
	  + m_numCpldDofs[c+1];
    }

    /// @brief Returns the tagged dof index
    ///
    /// Tags are stored without the shift (see getTagged()), so the stored
    /// value dofAt() is what is searched for.
    inline index_t tindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return std::distance(m_tagged.begin(),std::lower_bound(m_tagged.begin(),m_tagged.end(),dofAt(i,k,c)));
    }

    /// @brief Returns the boundary index of global dof \a gl.
    ///
    /// Produces undefined results if dof \a gl does not lie on the boundary.
    inline index_t global_to_bindex(index_t gl) const
    {
        GISMO_ASSERT( is_boundary_index( gl ),
                      "global_to_bindex(): dof "<<gl<<" is not on the boundary");
        return gl - m_numFreeDofs.back() - m_shift + m_bshift;
    //gl -= m_numFreeDofs.back() + m_shift;
	//const index_t c = std::distance(m_numElimDofs.begin(),
    //    std::upper_bound(m_numElimDofs.begin(), m_numElimDofs.end(),gl)) -1;
 	//return gl + m_bshift; // - m_numElimDofs[c]
    }

    // The is_*_index predicates take a shifted global index, as returned by
    // index(), and answer false for any index outside the shifted range
    // [m_shift, m_shift+size()) instead of classifying it: an index below
    // the shift belongs to no dof of this mapper, free or otherwise.

    /// Returns true if global dof \a gl is not eliminated.
    inline bool is_free_index(index_t gl) const
    {
      return shiftedInRange(gl, m_curElimId);
    }

    /// Returns true if local dof \a i of patch \a k is not eliminated.
    inline bool is_free( index_t i, index_t k = 0, index_t c = 0) const
    { return is_free_index( index(i, k, c) ); }

    /// Returns true if global dof \a gl is eliminated
    inline bool is_boundary_index( index_t gl ) const
    {
      return shiftedInRange(gl, m_numFreeDofs.back() + m_numElimDofs.back()) &&
          gl - m_shift >= m_numFreeDofs.back();
    }

    /// Returns true if local dof \a i of patch \a k is eliminated.
    inline bool is_boundary(index_t i, index_t k = 0, index_t c = 0) const
    {return is_boundary_index( index(i, k, c) );}

    /// Returns true if local dof \a i of patch \a k is coupled.
    inline bool is_coupled( index_t i, index_t k = 0, index_t c = 0) const
    { return  is_coupled_index( index(i, k, c) ); }

    /// Returns true if \a gl is a coupled dof.
    inline bool is_coupled_index(index_t gl) const
    {
      // Coupled dofs are free, so anything outside the free range is
      // answered here; it also keeps componentOfUnshifted() in range.
      if (!shiftedInRange(gl, m_numFreeDofs.back())) return false;
      const index_t g = gl - m_shift;
      const index_t gc = componentOfUnshifted(g);
      const index_t vv = m_numFreeDofs[gc+1];
      // The coupled dofs of a component sit at the top of that component's
      // own free block, so the band is that component's own coupled count
      // wide: the difference of the cumulative prefix m_numCpldDofs.
      const index_t nc = m_numCpldDofs[gc+1] - m_numCpldDofs[gc];
      return  (g < vv &&      // is a free dof of component gc, and
               g >= vv - nc); // lies in its coupled band
    }

    /// Returns true if local dof \a i of patch \a k is tagged.
    inline bool is_tagged(index_t i, index_t k = 0, index_t c = 0) const
    { return  is_tagged_index( index(i, k, c) ); }

    /// Returns true if \a gl is a tagged dof.
    inline bool is_tagged_index(index_t gl) const
    {
        // Tags are stored without the shift; see getTagged().  Every tag is
        // a valid index, so an out-of-range gl is not tagged.
        return shiftedInRange(gl, m_numFreeDofs.back() + m_numElimDofs.back()) &&
            std::binary_search(m_tagged.begin(),m_tagged.end(),gl - m_shift);
    }

    /// Returns the number of components present in the mapper
    inline index_t numComponents() const
    { return static_cast<index_t>(m_dofs.size()); }

    /// Returns the total number of dofs (free and eliminated).
    inline index_t size() const
    {
        GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return freeSize() + boundarySize();
    }

    /// Returns the total number of dofs (free and eliminated).
    inline index_t size(index_t comp) const
    {
        GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        ensureComponent(comp, "size");
        return m_numFreeDofs[comp+1]-m_numFreeDofs[comp]
	  + m_numElimDofs[comp+1]-m_numElimDofs[comp];
    }

    /// Returns the number of free (not eliminated) dofs.
    inline index_t freeSize() const
    {
      return m_curElimId;
    }

    inline index_t freeSize(index_t comp) const
    {
      ensureComponent(comp, "freeSize");
      return m_numFreeDofs[comp+1]-m_numFreeDofs[comp] +
	(allFree() ? m_numElimDofs[comp+1]-m_numElimDofs[comp] : 0 );
    }

    /// \brief Returns the sorted vector of tagged dofs.
    ///
    /// The entries are stored WITHOUT the shift, i.e. they are freeIndex()
    /// values and not index() values, so that tagging is unaffected by a
    /// later setShift() or addShift().  Add the shift to compare them with
    /// the global indices returned by index().
    const std::vector<index_t> & getTagged() const { return m_tagged; }

    /// Returns the number of coupled (not eliminated) dofs.
    index_t coupledSize() const;

    /// Returns the number of tagged dofs.
    index_t taggedSize() const;

    /// Returns the number of eliminated dofs.
    inline index_t boundarySize() const
    {
      GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return m_numElimDofs.back();
    }

    index_t boundarySizeWithDuplicates() const;

    /// Returns the offset corresponding to patch \a k for component \a c.
    /// Zero for every real patch under the global-identity/aliased layout
    /// (see Layout).
    size_t offset(index_t k, index_t c = 0) const
    {
        ensureComponent(c, "offset");
        ensurePatch(k, "offset");
        return offAt(c,k);
    }

    /// Returns the number of patches present underneath the mapper
    size_t numPatches() const {return m_nPatches;}

    /// \brief Returns the total number of patch-local degrees of
    /// freedom that are being mapped
    size_t mapSize() const
    {
        size_t s = 0;
        for (size_t c = 0; c != m_dofs.size(); ++c)
            s += m_dofs[c].size();
        return s;
    }

    size_t componentsSize() const {return m_dofs.size();}

    /// Returns the storage layout of this mapper (see Layout).
    /// Declared at construction and never inferred from an observed
    /// offset pattern.
    Layout layout() const { return m_layout; }

    /// Returns true if this mapper was declared, at construction, to
    /// have been built from more than one distinct per-component basis
    /// object. This is a DECLARED property, never inferred from sizes
    /// or offsets: a square-mesh Raviart-Thomas mapper is
    /// size-indistinguishable from a uniform one, yet its components are
    /// built from different bases.
    bool hasDistinctComponentSpaces() const { return m_hasDistinctComponentSpaces; }

    /// \brief Returns true if every component can be treated as one shared
    /// space, as a single-basis-per-space evaluator such as the expression
    /// assembler requires: the component spaces were not declared distinct,
    /// and every component has the same cardinality as component 0
    /// (patchSize(p,c)==patchSize(p,0) on every patch under the
    /// patch-concatenated layout, totalSize(c)==totalSize(0) under the
    /// global-identity layout).
    ///
    /// O(1): component sizes are fixed at construction, so the answer is
    /// computed there and cached, which lets callers check it on every
    /// assembly call.
    bool hasUniformComponents() const { return m_uniformComponents; }

    /// \brief Returns the total number of patch-local DoFs
    /// that live on patch \a k for component \a c.  Under the
    /// global-identity/aliased layout this is the component-global
    /// identity total for every patch (see Layout and
    /// setIdentity()), not just the last one.
    size_t patchSize(const index_t k, const index_t c = 0) const
    {
        ensureComponent(c, "patchSize");
        ensurePatch(k, "patchSize");
        return localCount(k,c);
    }

    size_t totalSize(const index_t c = 0) const
    {
        ensureComponent(c, "totalSize");
        return m_dofs[c].size();
    }

    /// \brief For \a gl being a global index, this function returns a
    /// vector of pairs (patch,dof) that contains all the pairs which
    /// map to \a gl
    ///
    /// \a gl includes the shift, as returned by index(), and must be a
    /// valid index (see componentOf()); the returned dofs are patch-local
    /// and never shifted.
    void preImage(index_t gl, std::vector<std::pair<index_t,index_t> > & result) const;

    /// \brief For \a gl being a global index, this function returns a
    /// pair (patch,dof) that maps to \a gl
    ///
    /// Same index conventions as preImage().
    std::pair<index_t,index_t> anyPreImage(index_t gl) const;

    /// \brief For all global index, this function assigns
    /// a pair (patch,dof) that maps to that global index
    ///
    /// The result has exactly size() entries, one per global index, at
    /// position \a gl minus the shift.  Entries of global indices that
    /// belong to another component than \a comp are (-1,-1).
    std::vector<std::pair<index_t,index_t> > anyPreImages(index_t comp = 0) const;

    /// \brief Produces the inverse of the mapping on patch \a k
    /// assuming that the map is invertible on that patch
    ///
    /// The result covers every component, keyed by global index (including
    /// the shift, as returned by index()) and valued by the patch-local
    /// index within the component that global index belongs to.  Only the
    /// dofs that live on patch \a k are considered; under the
    /// global-identity/aliased layout every patch carries the complete
    /// inverse, so the result is the same for every patch.
    std::map<index_t,index_t> inverseOnPatch(const index_t k) const;

    /// \brief For \a gl being a global index, this function returns
    /// true whenever \a gl corresponds to patch \a k
    ///
    /// \a gl includes the shift, as returned by index(); an index outside
    /// the mapper's range corresponds to no patch.
    bool indexOnPatch(const index_t gl, const index_t k, index_t & local) const;

    inline bool indexOnPatch(const index_t gl, const index_t k) const
    {
        index_t local;
        return indexOnPatch(gl, k, local);
    }

    /// \brief For \a n being an index which is already offset, it
    /// returns the global index where it is mapped to by the dof
    /// mapper.  Walks the (small) component list linearly rather than
    /// assuming every component has the same storage size, which does
    /// not hold for ragged storage.
    ///
    /// \a n must lie in [0, mapSize()).  The bound is checked by the walk
    /// itself -- running out of components means \a n was too large -- so
    /// the check costs nothing beyond the walk.
    inline index_t mapIndex(index_t n) const
    {
        const index_t n0 = n;
        GISMO_ENSURE(n >= 0, "gsDofMapper::mapIndex: negative index "<<n<<".");
        size_t c = 0;
        while (c != m_dofs.size() && static_cast<size_t>(n) >= m_dofs[c].size())
        {
            n -= static_cast<index_t>(m_dofs[c].size());
            ++c;
        }
        GISMO_ENSURE(c != m_dofs.size(), "gsDofMapper::mapIndex: index "<<n0
                     <<" is outside [0,"<<mapSize()<<").");
        return m_dofs[c][n] + m_shift;
    }

    /// \brief Returns all boundary dofs on patch k of component \a comp
    /// (local dof indices)
    gsVector<index_t> findBoundary(const index_t k, const index_t comp = 0) const;

    /// \brief Returns all free dofs on patch k of component \a comp
    /// (local dof indices)
    gsVector<index_t> findFree(const index_t k, const index_t comp = 0) const;

    /// \brief Returns all coupled dofs on patch k of component \a comp
    /// (local dof indices).  If \a j is a patch index, only the dofs that
    /// patch k shares with patch j are returned.
    gsVector<index_t> findCoupled(const index_t k, const index_t j = -1,
                                  const index_t comp = 0) const;

    /// \brief Returns all free, not coupled dofs on patch k of component
    /// \a comp (local dof indices)
    gsVector<index_t> findFreeUncoupled(const index_t k, const index_t comp = 0) const;

    /// \brief Returns all tagged dofs on patch k of component \a comp
    /// (local dof indices)
    gsVector<index_t> findTagged(const index_t k, const index_t comp = 0) const;

private:

    template<class Predicate, class Iterator>
      static gsVector<index_t> find_impl(Iterator istart, Iterator iend, Predicate pred);

    void finalizeComp(const index_t comp);

    // replace all references to oldIdx by newIdx
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx);
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx, index_t comp);

    void mergeDofsGlobally(index_t dof1, index_t dof2);
    void mergeDofsGlobally(index_t dof1, index_t dof2, index_t comp);

// Data members
private:

    // m_dofs/m_patchDofs stores for each patch the mapping from local to global dofs.
    //
    // During setup, the entries have a different meaning:
    //   0        -- regular free dof
    //   negative -- an eliminated dof
    //   positive -- a coupling dof
    // For nonzero entries, the value is an id which identifies the eliminated/coupling
    // group of the dof. Dofs with the same id will get the same dof index in the
    // final numbering stage in finalize().

    // Representation of each component as a single vector plus
    // offsets for patch-local indices
    std::vector<std::vector<index_t> >  m_dofs;

    /// Number of patches, shared by all components: ragged storage means
    /// different sizes per (component,patch), never different patch sets.
    size_t m_nPatches;

    /// Flat per-component offset table, row stride m_nPatches+1:
    /// entry [c*(m_nPatches+1)+k] is the start offset of patch k's local
    /// dofs within m_dofs[c] (k==m_nPatches is the sentinel, equal to
    /// m_dofs[c].size()).  Touch only through offAt()/offBegin()/offEnd().
    std::vector<size_t> m_offset;

    /// Storage layout (patch-concatenated or global-identity/aliased).
    /// Declared at construction, never inferred; see Layout.
    Layout m_layout;

    /// Declared (never inferred) at construction: true if this mapper
    /// was built from more than one distinct per-component basis object.
    bool m_hasDistinctComponentSpaces;

    /// hasUniformComponents(), cached at construction: nothing changes the
    /// component sizes or the declared distinctness afterwards.
    bool m_uniformComponents;

    /// Shifting of the global index (zero by default)
    ///
    /// Everything stored in this class -- m_dofs, m_tagged, the count
    /// vectors -- lives in the unshifted index space.  The shift is added
    /// only to the global indices a public method returns and subtracted
    /// from the global indices a public method takes, so changing it never
    /// requires touching the stored state.
    index_t m_shift;

    /// Shifting of the boundary index (zero by default)
    index_t m_bshift;

    /// Offsets of free dofs, nComp+1
    std::vector<index_t> m_numFreeDofs;
    /// Offsets of eliminated dofs, nComp+1
    std::vector<index_t> m_numElimDofs;
    /// Offsets of coupled dofs, nComp+1.  During setup, entry c+1 is instead
    /// the number of coupling ids handed out in component c; the ids are
    /// 1..count, so the count never exceeds the component's dof count.
    std::vector<index_t> m_numCpldDofs;

    // used during setup: running id for current eliminated dof
    // After finalize() is called m_curElimId takes positive value.
    index_t m_curElimId;

    /// Stores the tagged indices, sorted and unshifted
    std::vector<index_t> m_tagged;

}; // class gsDofMapper

/// Print (as string) a dofmapper structure
inline std::ostream& operator<<( std::ostream& os, const gsDofMapper& b )
{
    return b.print( os );
}

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsDofMapper
   */
  void pybind11_init_gsDofMapper(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo
