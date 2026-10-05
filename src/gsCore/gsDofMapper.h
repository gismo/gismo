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

#include <algorithm>
#include <unordered_map>
#include <limits>

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

    Storage modes.  By default (storage::dense) the local-to-global table
    holds one entry per local dof of every component, i.e. O(N_c) memory
    for a component with N_c local dofs.  A mapper may instead be created
    in storage::sparse mode, in which only the \em marked positions
    (interface, coupled, eliminated and collapsed dofs) are stored and the
    regular positions are numbered by counting; the memory is then
    O(M_c), M_c being the number of marked positions of component c.  The
    mode is chosen at construction or by setIdentity() (a full reset), and
    permuteFreeDofs() converts a sparse mapper to dense storage.  A sparse
    mapper returns, for every query, exactly what the dense mapper built by
    the same calls returns.  A single lookup in sparse mode costs O(log M_c)
    after finalize(); before finalize() it is one hash-map find (expected
    O(1)).  After localize() a sparse mapper holds a run table instead:
    O(#runs) memory, with #runs at most the number of maximal runs of
    local regular dofs plus M_c, and O(log #runs) lookup.

    \ingroup Core

*/
class GISMO_EXPORT gsDofMapper
{
private:

    /// One maximal run of a localized sparse component: positions
    /// [start, start+len) carry the unshifted values [val, val+len).
    /// Interleaved so that a lookup hit reads a single 3*sizeof(index_t)
    /// record. A component with R runs answers one position in
    /// O(log R) (binary search on start) and an ascending batch of n
    /// positions in O(log R + n) (see index_into).
    struct Run
    {
        index_t start, len, val;
    };

    struct DofUnionFind
    {
        /// Return the canonical setup-time label for the given value.
        index_t representative(index_t value);

        /// Merge the equivalence classes containing two labels.
        void unite(index_t first, index_t second);

        /// Return the allocated storage used by this disjoint-set forest.
        size_t nBytes() const;

    private:
        /// Find or create the node corresponding to a setup-time label.
        index_t nodeFor(index_t value);

        /// Find a forest root and apply path compression.
        index_t find(index_t node);

        std::unordered_map<index_t,index_t> m_nodes;
        std::vector<index_t>               m_parent;
        std::vector<unsigned char>         m_rank;
        std::vector<index_t>               m_label;
    };

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

    /// Storage of the local-to-global table (see the class documentation).
    enum class storage { dense, /**< one stored entry per position, O(N_c) per component */ sparse /**< only marked positions are stored (see the class documentation) */ };

    /// Returns the storage mode fixed at construction (or at the last full
    /// reset: setIdentity(), or permuteFreeDofs() which densifies).
    storage storageMode() const { return m_storage; }

    /// Default empty constructor
    gsDofMapper();

    /**
     * @brief construct a dof mapper with a given number of dofs per
     * patch
     *
     * @param patchDofSizes
     * @param nComp number of components
     * @param st    storage mode of the local-to-global table
     */
    gsDofMapper(const gsVector<index_t> &patchDofSizes, index_t nComp = 1,
                storage st = storage::dense)
    {
        initPatchDofs(patchDofSizes, nComp, st);
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
     * @param st                      storage mode of the local-to-global
     *                                table
     */
    gsDofMapper(const std::vector<gsVector<index_t> > & patchDofSizes,
                bool hasDistinctComponentSpaces, storage st = storage::dense);

    void swap(gsDofMapper & other)
    {
        m_dofs  .swap(other.m_dofs);
        std::swap(m_storage, other.m_storage);
        m_marked.swap(other.m_marked);
        m_keys  .swap(other.m_keys);
        m_vals  .swap(other.m_vals);
        m_regBase.swap(other.m_regBase);
        m_runs.swap(other.m_runs);
        std::swap(m_localizedSparse, other.m_localizedSparse);
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
        m_unionFind.swap(other.m_unionFind);
    }

private:

    /// Initialize by vector of DoF indices and dimension
    void initPatchDofs(const gsVector<index_t> & patchDofSizes,
                       index_t nComp = 1, storage st = storage::dense);

    /// Initialize by one patch-dof-size vector per component (ragged,
    /// patch-concatenated layout).
    void initRaggedPatchDofs(const std::vector<gsVector<index_t> > & patchDofSizes,
                             bool hasDistinctComponentSpaces,
                             storage st = storage::dense);

    /// Sets the storage mode of a full reset to \a st and empties every
    /// container of the other mode; for sparse storage the setup maps are
    /// created empty, one per component (\a nComp of them).  m_dofs is
    /// left sized to \a nComp with empty inner vectors; dense callers
    /// size them afterwards.
    void resetStorage(storage st, size_t nComp);

    /// Converts a sparse mapper to dense storage (O(sum_c N_c) memory).
    /// Requires a finalized mapper.
    void densify();

    // Flat-offset-table accessors (the only way this class touches
    // m_offset): m_offset is one row of m_nPatches+1 entries per
    // component, laid out consecutively.
    size_t offAt(index_t c, index_t k) const
    { return m_offset[static_cast<size_t>(c)*(m_nPatches+1)+static_cast<size_t>(k)]; }

    std::vector<size_t>::const_iterator offBegin(index_t c) const
    { return m_offset.begin() + static_cast<size_t>(c)*(m_nPatches+1); }

    std::vector<size_t>::const_iterator offEnd(index_t c) const
    { return offBegin(c) + (m_nPatches+1); }

    /// Logical number of positions N_c of component \a c, i.e. the end
    /// sentinel of its offset row.  Valid in both layouts and both storage
    /// modes, unlike m_dofs[c].size(), which is empty in sparse mode.
    size_t compSize(index_t c) const
    { return offAt(c, static_cast<index_t>(m_nPatches)); }

    /// Setup-time value (see the comment on the data members) of local dof
    /// \a i of patch \a k in component \a c.  In sparse mode this is a
    /// pure lookup: an absent position reads as 0 and no entry is created.
    inline index_t setupValue(index_t i, index_t k, index_t c) const
    {
        GISMO_ASSERT(validLocal(i,k,c), "gsDofMapper: invalid local dof "<<i<<" of patch "<<k
                     <<", component "<<c<<localRangeInfo(k,c));
        const size_t p = offAt(c,k)+i;
        if (storage::dense == m_storage) return m_dofs[c][p];
        const std::unordered_map<index_t,index_t>::const_iterator it =
            m_marked[c].find(static_cast<index_t>(p));
        return m_marked[c].end() == it ? 0 : it->second;
    }

    /// Stores the nonzero setup-time value \a v at local dof \a i of patch
    /// \a k in component \a c.
    inline void setSetupValue(index_t i, index_t k, index_t c, index_t v)
    {
        GISMO_ASSERT(validLocal(i,k,c), "gsDofMapper: invalid local dof "<<i<<" of patch "<<k
                     <<", component "<<c<<localRangeInfo(k,c));
        const size_t p = offAt(c,k)+i;
        if (storage::dense == m_storage) m_dofs[c][p] = v;
        else m_marked[c][static_cast<index_t>(p)] = v;
    }

    /// Value stored for position \a p of component \a c in sparse mode.
    /// Before finalize() this is the setup-time value (0 if absent);
    /// afterwards a binary search over the marked positions, a hit giving the
    /// stored final id and a miss the id of the (p - r)-th regular position, r
    /// being the number of marked positions below \a p.  O(log M_c).  After
    /// localize() the value comes from the run table through the inline
    /// runValue(): O(log #runs).  This out-of-line function itself serves the
    /// setup path and the finalized non-localized path; callers on a hot path
    /// test m_localizedSparse and call runValue() directly.
    index_t sparseValue(index_t c, size_t p) const;

    /// Value at position \a p of component \a c of a localized sparse mapper:
    /// the run with the largest start <= p, remoteDof() if no run covers p.
    /// O(log #runs).  Requires m_localizedSparse.
    inline index_t runValue(index_t c, size_t p) const
    {
        const index_t q = static_cast<index_t>(p);
        const std::vector<Run> & runs = m_runs[c];
        std::vector<Run>::const_iterator it = std::upper_bound(runs.begin(), runs.end(), q,
            [](index_t x, const Run & r) { return x < r.start; });
        if (it == runs.begin()) return remoteDof();
        --it;
        const index_t d = q - it->start;
        return d < it->len ? it->val + d : remoteDof();
    }

    /// Value stored for position \a p of component \a c (either mode).
    inline index_t valueAtPos(index_t c, size_t p) const
    {
        if (storage::dense == m_storage) return m_dofs[c][p];
        if (m_localizedSparse) return runValue(c, p);
        return sparseValue(c, p);
    }

    /// Calls f(p, value) for the positions p in [\a pBegin, \a pEnd) of
    /// component \a c in increasing order, until f returns true.  Sparse
    /// finalized storage is walked by merging with the marked positions (or,
    /// after localize(), with the runs): O(pEnd - pBegin + log M_c)
    /// (O(pEnd - pBegin + log #runs)), no temporary storage.  Before
    /// finalize() the walk does one hash lookup per position: O(pEnd - pBegin)
    /// expected.  Dense storage is a plain loop.
    template<class F>
    void forEachValue(index_t c, size_t pBegin, size_t pEnd, F f) const;

    /// Post-finalize read of the stored (unshifted) value of local dof \a i
    /// of patch \a k in component \a c.
    ///
    /// The bounds check is debug-only.  dofAt() runs once per local dof per
    /// element in every assembler, so the per-dof accessors built on it
    /// (index(), bindex(), cindex(), tindex(), freeIndex(), and the
    /// is_free/is_boundary/is_coupled/is_tagged queries taking a local dof)
    /// and the per-element localToGlobal()/localToGlobal2() are not checked
    /// in Release builds; every other public entry point validates its
    /// arguments with the ensure*() helpers below before reaching it.
    ///
    /// Cost: dense storage O(1); localized sparse storage O(log #runs),
    /// inline; finalized non-localized sparse storage O(log M_c), out of line.
    inline index_t dofAt(index_t i, index_t k, index_t c) const
    {
        GISMO_ASSERT(validLocal(i,k,c), "gsDofMapper: invalid local dof "<<i<<" of patch "<<k
                     <<", component "<<c<<localRangeInfo(k,c));
        const size_t p = offAt(c,k)+i;
        if (storage::dense == m_storage) return m_dofs[c][p];
        if (m_localizedSparse) return runValue(c, p);
        return sparseValue(c, p);
    }

    // --- argument validation -------------------------------------------
    //
    // The valid*() predicates never touch storage for an invalid argument.
    // The ensure*() helpers throw std::runtime_error (GISMO_ENSURE) in
    // Release builds too: they guard memory safety, so they must not
    // compile away.  \a where names the public method for the message.

    bool validComponent(index_t c) const
    { return c >= 0 && c < numComponents(); }

    bool validPatch(index_t k) const
    { return k >= 0 && static_cast<size_t>(k) < m_nPatches; }

    /// Number of valid local indices of (patch \a k, component \a c), for
    /// arguments already known to be valid.  Patch-specific under the
    /// patch-concatenated layout, so that an oversized local index cannot
    /// spill into the next patch's storage; the component's global identity
    /// total on every patch under the aliased layout.
    size_t localCount(index_t k, index_t c) const
    {
        return GlobalIdentity == m_layout ? compSize(c)
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
            os << ": the mapper has " << numComponents() << " components.";
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
            for (index_t cc = 0; cc != numComponents(); ++cc)
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

        Requires isPermutation(); otherwise a std::runtime_error is thrown,
        e.g. for coupled interfaces, or for a localized mapper with remote
        positions (see localize()).

        Complexity: O(size() + compSize(\a comp)) for dense storage; sparse
        storage adds O(log M_c) or O(log #runs) for the walk start, as in
        forEachValue().
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
    ///
    /// This is a full reset: the storage mode becomes \a st, so calling it
    /// with the default argument on a sparse mapper returns the mapper to
    /// dense storage.
    void setIdentity(index_t nPatches, size_t nDofs, size_t nComp = 1,
                     storage st = storage::dense);

    ///\brief Set this mapping to be the identity, with a (possibly
    /// unequal) total dof count per component.  The scalar-size overload
    /// delegates here by broadcasting one total to every component.
    ///
    /// This is a full reset: the storage mode becomes \a st, so calling it
    /// with the default argument on a sparse mapper returns the mapper to
    /// dense storage.
    void setIdentity(index_t nPatches, const std::vector<size_t> & dofsPerComponent,
                     storage st = storage::dense);

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
    /// A mapper in sparse storage is converted to dense storage (O(sum_c N_c)
    /// memory) once the permutation has been validated; storageMode() reports
    /// storage::dense afterwards.
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

    /** \brief Global indices of the local dofs in column \a col of \a act
     *  (patch \a patch, component \a comp), written to \a out.
     *
     *  out[i] == index(act(i,col), patch, comp) for i in [0, act.rows()),
     *  the shift included and remote dofs reported as remoteDof() plus the
     *  shift.  finalize() must have been called.  \a out must hold
     *  act.rows() entries and must not alias \a act.
     *
     *  Complexity: O(n) for dense storage, n = act.rows().  For sparse
     *  storage one binary search locates the first entry, after which a
     *  cursor into the (sorted) key or run table advances by galloping:
     *  O(log R + n + sum log g), g being the number of table entries
     *  skipped between consecutive inputs, which is O(log R + n) for the
     *  lexicographic actives of a tensor B-spline basis (R = number of
     *  marked positions before localize(), number of runs after).  Any
     *  input order is correct; an entry smaller than its predecessor
     *  restarts the cursor with a binary search, O(log R) extra per descent.
     *
     *  \note On the assembly hot path: arguments are checked in debug
     *  builds only.
     */
    void index_into(const gsMatrix<index_t> & act, index_t col,
                    index_t patch, index_t comp, index_t * out) const;

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
        GISMO_ASSERT( is_boundary_index( gl ) && !is_remote_index( gl ),
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

    /**
       @brief Turns the mapper into a rank-local one.

       The free dofs listed in \a localDofs (indices as returned by
       index(), sorted ascending, unique) are renumbered to
       0,...,localDofs.size()-1 (plus the shift), in the given order.
       All other free dofs become \em remote: index() returns a value
       for which is_remote_index() is true, and they must not be used
       for assembly or evaluation. Eliminated dofs keep their boundary
       numbering, so the fixed (Dirichlet) values remain valid.

       After the call, freeSize() is the number of local dofs, so the
       matrices, right-hand sides and solution vectors of an assembler
       using this mapper are rank-local. \a localDofs is the local to
       global map of the free dofs.

       Typical use: \a localDofs are all free dofs active on the
       elements of the rank (owned and ghost dofs).

       The numbering is the same in both storage modes, and localizing an
       already localized mapper is allowed (the free dofs are then the
       current local ones).  Remote dofs are treated as follows.
       - index() returns remoteDof() plus the shift, and is_remote_index()
         is true for it;
       - is_free_index(), is_boundary_index(), is_coupled_index() and
         is_tagged_index() are false for it;
       - eliminated dofs keep their boundary numbering, moved down to start
         at the new free size, so global_to_bindex() and the fixed values
         remain valid;
       - the whole-array queries that compare stored values with a
         threshold see remote positions as values above every free id:
         findBoundary() lists them and boundarySizeWithDuplicates() counts
         them, whereas findFree(), findFreeUncoupled() and findCoupled()
         exclude them.

       In sparse storage the result is held as a run table: per component
       the maximal runs of consecutive positions carrying consecutive local
       ids, every other position being remote.  Memory is O(#runs), with
       #runs at most the number of maximal runs of local regular dofs plus
       M_c, the number of marked positions of component c.  A lookup costs
       O(log #runs), and no loop runs over all positions.  Time is O(n_c + M_c + M_c log n)
       for the first localization and O(#runs log n + n) for a repeated one,
       n being localDofs.size() and n_c the part of it in the regular id
       range of component c.
    */
    void localize(const std::vector<index_t> & localDofs);

    /// Returns true if global dof \a gl was dropped by localize()
    inline bool is_remote_index(index_t gl) const
    { return gl - m_shift >= remoteDof(); }

    /// Returns true if local dof \a i of patch \a k was dropped by localize()
    inline bool is_remote(index_t i, index_t k = 0, index_t c = 0) const
    { return is_remote_index( index(i, k, c) ); }

    /// The (shift-less) index stored for dofs dropped by localize()
    static index_t remoteDof()
    { return std::numeric_limits<index_t>::max() / 2; }

    /// Returns the number of components present in the mapper
    inline index_t numComponents() const
    { return static_cast<index_t>(m_dofs.size()); }  // m_dofs is sized to the component count in both storage modes

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
        for (index_t c = 0; c != numComponents(); ++c)
            s += compSize(c);
        return s;
    }

    size_t componentsSize() const { return static_cast<size_t>(numComponents()); }

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
        return compSize(c);
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
    ///
    /// On a localized mapper (see localize()) the result still has size()
    /// entries, that is, the local free count plus the eliminated count.
    /// Slot \a gl minus the shift of every local free dof and every
    /// eliminated dof of \a comp holds a (patch, patch-local index) pair
    /// that maps to it, the one at the lowest position.  Remote positions are
    /// skipped and contribute to no slot; a remote dof has no slot, as it has
    /// no global index on this rank.
    ///
    /// Complexity: O(size() + compSize(\a comp) log numPatches()), plus the
    /// storage-dependent walk start of forEachValue().
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
        index_t c = 0;
        while (c != numComponents() && static_cast<size_t>(n) >= compSize(c))
        {
            n -= static_cast<index_t>(compSize(c));
            ++c;
        }
        GISMO_ENSURE(c != numComponents(), "gsDofMapper::mapIndex: index "<<n0
                     <<" is outside [0,"<<mapSize()<<").");
        return valueAtPos(c, static_cast<size_t>(n)) + m_shift;
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

    /// \brief Returns the capacity in bytes of this mapper's containers
    /// (not their logical size).  The hash maps of sparse setup are
    /// estimated from their bucket count and stored pairs, without node
    /// links or allocator overhead, so the result is an estimate and not a
    /// heap measurement.
    inline size_t nBytes() const
    {
        size_t bytes = sizeof(*this);

        bytes += m_dofs.capacity() * sizeof(std::vector<index_t>);
        for (size_t c = 0; c != m_dofs.size(); ++c)
            bytes += m_dofs[c].capacity() * sizeof(index_t);

        bytes += m_marked.capacity() * sizeof(std::unordered_map<index_t,index_t>);
        for (const std::unordered_map<index_t,index_t> & mk : m_marked)
        {
            bytes += mk.bucket_count() * sizeof(void *);
            bytes += mk.size() * sizeof(std::pair<const index_t,index_t>);
        }
        bytes += m_keys.capacity() * sizeof(std::vector<index_t>);
        for (const std::vector<index_t> & v : m_keys)
            bytes += v.capacity() * sizeof(index_t);
        bytes += m_vals.capacity() * sizeof(std::vector<index_t>);
        for (const std::vector<index_t> & v : m_vals)
            bytes += v.capacity() * sizeof(index_t);
        bytes += m_regBase.capacity() * sizeof(index_t);
        bytes += m_runs.capacity() * sizeof(std::vector<Run>);
        for (const std::vector<Run> & v : m_runs)
            bytes += v.capacity() * sizeof(Run);

        bytes += m_offset.capacity() * sizeof(size_t);

        bytes += m_numFreeDofs.capacity() * sizeof(index_t);
        bytes += m_numElimDofs.capacity() * sizeof(index_t);
        bytes += m_numCpldDofs.capacity() * sizeof(index_t);
        bytes += m_tagged.capacity()      * sizeof(index_t);

        bytes += m_unionFind.capacity() * sizeof(DofUnionFind);
        for (const DofUnionFind & uf : m_unionFind)
            bytes += uf.nBytes();

        return bytes;
    }

private:

    template<class Predicate, class Iterator>
      static gsVector<index_t> find_impl(Iterator istart, Iterator iend, Predicate pred);

    /// Local indices on patch \a k of component \a comp whose stored value
    /// satisfies \a pred; the sparse-storage counterpart of find_impl().
    template<class Predicate>
    gsVector<index_t> findSparse(const index_t k, const index_t comp, Predicate pred) const;

    void finalizeComp(const index_t comp);

    /// Sparse-storage counterpart of finalizeComp(): replays the dense
    /// numbering of component \a comp from its marked positions.
    /// O(M_c log M_c) time and O(M_c) memory.
    void finalizeCompSparse(const index_t comp);

    // Merge the equivalence classes represented by oldIdx and newIdx.
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx);
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx, index_t comp);

    void mergeDofsGlobally(index_t dof1, index_t dof2);
    void mergeDofsGlobally(index_t dof1, index_t dof2, index_t comp);

    /// Return the current representative of a setup-time equivalence label.
    index_t canonicalDof(index_t dof, index_t comp);

    /// Reset the setup-time union-find state for the requested components.
    void resetUnionFind(size_t nComp)
    {
        m_unionFind.clear();
        m_unionFind.resize(nComp);
    }

// Data members
private:

    // m_dofs/m_patchDofs stores for each patch the mapping from local to global dofs
    // (dense storage only; in sparse storage every inner vector is empty, while the
    // outer vector stays sized to the component count, which is what numComponents()
    // reports in both modes).
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

    /// Storage mode, fixed by construction or by a full reset.
    storage m_storage;

    // Sparse storage.  Positions are per component, p = offAt(c,k)+i in [0,N_c).
    //
    // During setup, m_marked[c] maps a marked position to its setup-time value
    // (never 0: an absent position is a regular free dof).  finalize() replaces
    // it by the sorted marked positions m_keys[c] with their final unshifted
    // ids m_vals[c], and by m_regBase[c], the id of the first regular position;
    // the regular positions are numbered consecutively in position order.
    // m_marked is released by finalize().  localize() releases m_keys, m_vals
    // and m_regBase in turn and replaces them by a run table.
    std::vector<std::unordered_map<index_t,index_t> > m_marked;
    std::vector<std::vector<index_t> >  m_keys, m_vals;
    std::vector<index_t>                m_regBase;

    // Sparse storage after localize(): per component, maximal runs of
    // consecutive positions with consecutive values, sorted by start position.
    // A run never mixes free values (< m_curElimId) with eliminated ones.  A
    // position covered by no run is remote (remoteDof()).  m_keys, m_vals and
    // m_regBase are empty in this state.
    std::vector<std::vector<Run> >      m_runs;
    bool                                m_localizedSparse;

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

    /// Setup-time equivalence classes, one disjoint-set forest per component.
    /// This is released by finalize() after the flat mapper has been relabeled.
    std::vector<DofUnionFind> m_unionFind;

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
