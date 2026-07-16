/** @file gsIndexSubDomain.h

    @brief A subdomain defined by a list of global element indices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsSubDomain.h>
#include <algorithm>
#include <map>
#include <tuple>

namespace gismo
{

// Forward declarations
template<class T> class gsIndexSubDomain;
template<class T> class gsIndexSubDomainFilterIterator;
template<class T> class gsIndexSubDomainPatchView;


/**
   @brief Iterator for gsIndexSubDomain.

   Walks only the subset of elements specified by the index list, by
   maintaining a cached parent-domain iterator that advances incrementally.

   Invariant: after every operation, m_parentIt sits at parent element
   m_indices[m_id]. The constructor establishes this by seeking forward
   from position 0 to m_indices[0].

   \ingroup Domain
*/
template<class T>
class gsIndexSubDomainIterator : public gsDomainIterator<T>
{
    typedef gsDomainIterator<T> Base;
    typedef typename gsDomainIterator<T>::uPtr uPtr;

    const gsIndexSubDomain<T>&  m_subdomain;
    const std::vector<index_t>& m_indices;  // reference into subdomain's sorted list
    gsDomainIteratorWrapper<T>  m_parentIt; // cached; always at m_indices[m_id]

public:
    explicit gsIndexSubDomainIterator(const gsIndexSubDomain<T>& subdomain);

    // Default copy ctor: m_parentIt clones its internal iterator (already at the right position)
    gsIndexSubDomainIterator(const gsIndexSubDomainIterator&) = default;

    uPtr clone() const override
    { return uPtr(new gsIndexSubDomainIterator(*this)); }

    gsVector<T> lowerCorner() const override { return m_parentIt.lowerCorner(); }
    gsVector<T> upperCorner() const override { return m_parentIt.upperCorner(); }

    // Delegate to the parent iterator so cross-patch composite domains report the right patch
    index_t patchIndex() const override { return m_parentIt.patchIndex(); }

    // Return per-patch local id (matches gsCompositeDomainIterator::localId() semantics)
    size_t localId() const override { return m_parentIt.localId(); }

private:
    void next() override { next(1); }

    // Jump the parent by the delta between the two sorted indices so that one
    // bulk advance replaces k individual steps. Called before nextId(k).
    void next(index_t increment) override
    {
        const index_t target = static_cast<index_t>(this->m_id) + increment;
        if (target < static_cast<index_t>(m_indices.size()))
            m_parentIt += static_cast<index_t>(m_indices[target])
                        - static_cast<index_t>(m_indices[this->m_id]);
        // if target >= size we are stepping to end; parent position is irrelevant
    }

    void reset() override;

}; // class gsIndexSubDomainIterator


/**
   @brief A sub-domain defined by a sorted list of global element indices.

   Constructed from a shared_ptr to a parent gsDomain and a (possibly
   unsorted, possibly duplicate) vector of global element indices.  The
   constructor sorts and deduplicates the index vector.

   The subdomain co-owns the parent via shared_ptr, so the parent lifetime
   is automatically managed — no dangling references.

   The primary use case is restricting gsExprAssembler integration to a
   METIS-computed partition:
   \code
   typename gsDomain<T>::Ptr dom = mb.domain();
   auto sub = memory::make_shared(new gsIndexSubDomain<T>(dom, myIndices));
   A.setIntegrationDomain(sub);
   A.assemble(igrad(u), igrad(v) * meas(G));
   \endcode

   \ingroup Domain
*/
template<class T>
class gsIndexSubDomain : public gsSubDomain<T>
{
    typedef gsSubDomain<T> Base;

    typename gsDomain<T>::Ptr m_parentPtr; // co-owns the parent domain
    std::vector<index_t>      m_indices;   // sorted, unique, global element IDs

    // Cache of subdomain(k) results, keyed by patch. subdomain() is called
    // twice per boundary/interface loop (once for beginBdr, once for
    // endBdr; see gsExprAssembler.h call sites) and repeatedly across
    // assemble() calls. Pre-populated for every patch eagerly at
    // construction (Fix 1: read-only during assembly -- no OpenMP race);
    // subdomain() below is a pure lookup, never inserts.
    std::map<index_t, typename gsDomain<T>::Ptr> m_subdomainCache;

    // Whether the parent domain provides an analytic elementIndex() fast
    // path (probed once at construction against the parent's very first
    // element -- a structural property of the domain *type*, not of the
    // data, so one probe determines it for the whole domain). When true,
    // m_volumeBoxes below is never built; volumeElementAt() always resolves
    // through the fast path.
    bool m_hasFastElementIndex = false;

    // Built eagerly at construction (Fix 1: read-only during assembly -- no
    // OpenMP race), but only when the parent lacks an elementIndex() fast
    // path (e.g. THB domains). (patch, lowerCorner, upperCorner) per global
    // element id, covering the *whole* parent domain (not just m_indices).
    // Used for point-location of boundary/interface elements as a fallback;
    // see volumeElementAt().
    std::vector<std::tuple<index_t, gsVector<T>, gsVector<T>>> m_volumeBoxes;

    // Called once from the constructor (serial by definition). Probes
    // elementIndex() capability and, in debug builds, validates it against
    // the beginAll() reference for every element -- a elementIndex()
    // override returning ids out of beginAll() order would silently
    // corrupt boundary/interface ownership filtering, so this is a
    // correctness invariant, not just a sanity check.
    void initVolumeLookup()
    {
        auto it  = m_parentPtr->beginAll();
        auto end = m_parentPtr->endAll();
        if (it == end) { m_hasFastElementIndex = true; return; } // empty parent: nothing to probe or scan
        m_hasFastElementIndex =
            m_parentPtr->elementIndex(static_cast<index_t>(it.patchIndex()), it.centerPoint()) >= 0;
        if (!m_hasFastElementIndex)
            m_volumeBoxes.resize(m_parentPtr->numElements());
#ifndef NDEBUG
        for (; it != end; ++it)
        {
            if (!m_hasFastElementIndex)
                m_volumeBoxes[it.id()] = std::make_tuple(
                    static_cast<index_t>(it.patchIndex()), it.lowerCorner(), it.upperCorner());
            const index_t fast = m_parentPtr->elementIndex(
                static_cast<index_t>(it.patchIndex()), it.centerPoint());
            GISMO_ASSERT(fast < 0 || fast == static_cast<index_t>(it.id()),
                "gsIndexSubDomain: elementIndex() fast lookup disagrees with "
                "the beginAll() reference at global element "<<it.id()<<
                " (patch "<<it.patchIndex()<<"): got "<<fast<<", expected "<<
                it.id()<<". A gsDomain::elementIndex() override is returning "
                "ids out of beginAll() order.");
            // Round-trip check for the id-based path (Fix 3): a local id
            // (as returned by a raw per-patch iterator's id()) converted via
            // globalElementId() must land back on this same global id -- the
            // invariant gsIndexSubDomainFilterIterator::isOwned() (volume
            // case) and gsIndexSubDomainPatchView::numElements() rely on.
            const index_t roundTrip = m_parentPtr->globalElementId(
                static_cast<index_t>(it.patchIndex()), static_cast<index_t>(it.localId()));
            GISMO_ASSERT(roundTrip == static_cast<index_t>(it.id()),
                "gsIndexSubDomain: globalElementId(patch, localId) disagrees "
                "with the beginAll() reference at global element "<<it.id()<<
                " (patch "<<it.patchIndex()<<", local id "<<it.localId()<<
                "): got "<<roundTrip<<", expected "<<it.id()<<".");
        }
#else
        if (!m_hasFastElementIndex)
            for (; it != end; ++it)
                m_volumeBoxes[it.id()] = std::make_tuple(
                    static_cast<index_t>(it.patchIndex()), it.lowerCorner(), it.upperCorner());
#endif
    }

public:
    typedef typename gsDomain<T>::Ptr      Ptr;
    typedef typename gsDomain<T>::iterator iterator;

    /**
       @brief Construct from a parent domain (shared ownership) and element indices.

       @param parentPtr  Shared pointer to the full domain.  The subdomain
                         increments the reference count so the parent stays
                         alive as long as any subdomain references it.
       @param indices    Global element indices (unsorted/duplicate allowed;
                         constructor normalises them).
    */
    gsIndexSubDomain(typename gsDomain<T>::Ptr parentPtr, std::vector<index_t> indices)
    : Base(), m_parentPtr(give(parentPtr)), m_indices(give(indices))
    {
        GISMO_ASSERT(m_parentPtr, "Parent domain pointer must not be null.");
        std::sort(m_indices.begin(), m_indices.end());
        m_indices.erase(std::unique(m_indices.begin(), m_indices.end()),
                        m_indices.end());
        initVolumeLookup();
        // Eagerly pre-populate subdomain() for every patch (Fix 1: makes
        // subdomain() a pure, race-free lookup during assembly).
        for (index_t k = 0; k < static_cast<index_t>(m_parentPtr->nPieces()); ++k)
            m_subdomainCache[k] =
                memory::make_shared(new gsIndexSubDomainPatchView<T>(*this, k));
    }

    const gsDomain<T>& parentDomain() const override { return *m_parentPtr; }

    const std::vector<index_t>& elementIndices() const override { return m_indices; }

    bool contains(index_t elementIndex) const override
    {
        return std::binary_search(m_indices.begin(), m_indices.end(), elementIndex);
    }

    size_t numElements() const override { return m_indices.size(); }

    // Quadrature rule selection in gsExprEvaluator queries subdomain(patchIndex)
    // and nPieces(). A plain delegation to the parent (as before) would hand
    // back the *unrestricted* per-patch domain, so boundary/interface
    // iteration on the result would walk every element of the patch instead
    // of only the ones owned by this subdomain (they get inserted nranks times
    // under ADD_VALUES). Wrap the parent's per-patch domain in a filtering
    // view instead; volume queries (degree/dim/boundingBox/...) are still
    // forwarded straight through.
    typename gsDomain<T>::Ptr subdomain(index_t k) const override
    {
        auto cacheIt = m_subdomainCache.find(k);
        GISMO_ASSERT(cacheIt != m_subdomainCache.end(),
            "gsIndexSubDomain::subdomain: patch "<<k<<" out of range -- "
            "all patches are pre-populated at construction.");
        return cacheIt->second;
    }

    size_t nPieces() const override { return m_parentPtr->nPieces(); }

    /// Locates the global (parent beginAll()-order) element id of the volume
    /// element on \a patch containing the point \a u. Used to test whether a
    /// boundary/interface element (which lives in a different index space)
    /// is owned by this subdomain. Uses the parent domain's own
    /// elementIndex() analytic fast path when available (see
    /// gsDomain::elementIndex, capability probed once at construction);
    /// falls back to an O(parent numElements()) linear scan over the
    /// eagerly-built m_volumeBoxes for domain types that don't provide one.
    index_t volumeElementAt(index_t patch, const gsVector<T>& u) const
    {
        if (m_hasFastElementIndex)
        {
            const index_t fast = m_parentPtr->elementIndex(patch, u);
            GISMO_ASSERT(fast >= 0, "gsIndexSubDomain: elementIndex() fast "
                "path unexpectedly failed for a domain type that reported "
                "having one at construction.");
            return fast;
        }

        for (size_t id = 0; id < m_volumeBoxes.size(); ++id)
        {
            const auto & box = m_volumeBoxes[id];
            if (std::get<0>(box) != patch) continue;
            const gsVector<T> & lo = std::get<1>(box);
            const gsVector<T> & hi = std::get<2>(box);
            bool inside = true;
            for (index_t d = 0; d < u.size(); ++d)
                if (u[d] < lo[d] || u[d] > hi[d]) { inside = false; break; }
            if (inside) return static_cast<index_t>(id);
        }
        GISMO_ERROR("gsIndexSubDomain: could not locate the volume element "
                    "containing the given point on patch "<<patch<<".");
    }

    /// Maps a boundary-element iterator to the global id of the volume
    /// element it borders. Tries the purely-integer
    /// gsDomainIterator::adjacentVolumeLocalId() path first (no
    /// floating-point nudge, no risk of the nudge landing on the wrong side
    /// of a knot for thin elements); falls back to nudging the boundary
    /// element's center point half a cell into the interior along the
    /// side's normal direction and locating it via volumeElementAt().
    index_t volumeElementIdFromBoundary(const gsDomainIteratorWrapper<T> & bdrIt) const
    {
        const index_t localId = bdrIt.adjacentVolumeLocalId();
        if (localId >= 0)
            return m_parentPtr->globalElementId(bdrIt.patchIndex(), localId);

        gsVector<T> u = bdrIt.centerPoint();
        const boxSide bs = bdrIt.side();
        const short_t dir = bs.direction();
        const T eps = static_cast<T>(0.5) * bdrIt.getPerpendicularCellSize();
        u[dir] += (bs.parameter() ? -eps : eps);
        return volumeElementAt(bdrIt.patchIndex(), u);
    }

    iterator beginAll() const override
    {
        if (m_indices.empty())
            return this->endAll(); // degenerate: empty subdomain
        return iterator(new gsIndexSubDomainIterator<T>(*this));
    }

    short_t degree(short_t i = 0) const override { return m_parentPtr->degree(i); }
    short_t dim()                  const override { return m_parentPtr->dim(); }

    gsMatrix<T> boundingBox() const override { return m_parentPtr->boundingBox(); }
    gsMesh<T>   mesh()        const override { return m_parentPtr->mesh(); }

    std::ostream& print(std::ostream& os) const override
    {
        os << "gsIndexSubDomain with " << m_indices.size() << " / "
           << m_parentPtr->numElements() << " elements.";
        return os;
    }

}; // class gsIndexSubDomain


/**
   @brief Filtering iterator used by gsIndexSubDomainPatchView.

   Wraps a raw iterator range from the parent's per-patch domain (either
   volume or boundary elements) and skips any element not owned by the
   gsIndexSubDomain, per gsIndexSubDomain::contains(). For boundary elements
   ownership is decided via gsIndexSubDomain::volumeElementIdFromBoundary()
   (boundary elements live in a different index space than volume elements).

   \ingroup Domain
*/
template<class T>
class gsIndexSubDomainFilterIterator : public gsDomainIterator<T>
{
    typedef gsDomainIterator<T> Base;
    typedef typename Base::uPtr uPtr;

    const gsIndexSubDomain<T>& m_owner;
    gsDomainIteratorWrapper<T> m_it;
    gsDomainIteratorWrapper<T> m_end;
    bool m_boundary; // true: m_it walks boundary elements; false: volume elements

public:
    gsIndexSubDomainFilterIterator(const gsIndexSubDomain<T>& owner,
                                    gsDomainIteratorWrapper<T> it,
                                    gsDomainIteratorWrapper<T> end,
                                    bool isBoundary)
    : Base(0, (isBoundary && it != end) ? it.side() : boxSide(boundary::none)),
      m_owner(owner), m_it(give(it)), m_end(give(end)), m_boundary(isBoundary)
    {
        if (m_it != m_end)
            this->setPatchIndex(m_it.patchIndex());
        skipToOwned();
    }

    gsIndexSubDomainFilterIterator(const gsIndexSubDomainFilterIterator&) = default;

    uPtr clone() const override
    { return uPtr(new gsIndexSubDomainFilterIterator(*this)); }

    gsVector<T> lowerCorner() const override { return m_it.lowerCorner(); }
    gsVector<T> upperCorner() const override { return m_it.upperCorner(); }
    index_t patchIndex() const override { return m_it.patchIndex(); }
    size_t localId() const override { return m_it.localId(); }
    const T getPerpendicularCellSize() const override { return m_it.getPerpendicularCellSize(); }
    bool isBoundaryElement() const override { return m_it.isBoundaryElement(); }

private:
    void next() override { next(1); }

    void next(index_t increment) override
    {
        for (index_t i = 0; i < increment; ++i)
        {
            if (m_it == m_end) return;
            ++m_it;
            skipToOwned();
        }
    }

    // Not supported: the underlying per-patch boundary iterators (e.g.
    // gsTensorDomainBoundaryIterator) don't implement clone(), so there is
    // no way to stash a rewindable copy of the start position -- and, like
    // the base gsDomainIterator::reset() default, this path is never
    // exercised by the assembler's own boundary/interface iteration, which
    // always constructs a fresh begin/end pair instead of resetting one.
    void reset() override { GISMO_NO_IMPLEMENTATION }

    void skipToOwned()
    {
        while (m_it != m_end && !isOwned())
            ++m_it;
    }

    bool isOwned() const
    {
        // Note: m_it.id() (when m_it walks a raw per-patch domain, as
        // returned by e.g. gsCompositeDomain::subdomain(k)) is the *local*
        // element id within that patch, not the global id m_owner.contains()
        // expects (gsCompositeDomainIterator::localId() vs id()).
        //
        // Boundary case: volumeElementIdFromBoundary() locates the volume
        // id (via adjacentVolumeLocalId() when available, else point
        // location -- see gsIndexSubDomain.h).
        //
        // Volume case: m_it.id() IS the local id, and
        // gsDomain::globalElementId() converts local -> global directly, no
        // point location needed at all (Fix 3 -- round-trip-validated
        // against beginAll() in gsIndexSubDomain::initVolumeLookup()).
        const index_t volId = m_boundary
            ? m_owner.volumeElementIdFromBoundary(m_it)
            : m_owner.parentDomain().globalElementId(
                  static_cast<index_t>(m_it.patchIndex()), static_cast<index_t>(m_it.id()));
        return m_owner.contains(volId);
    }

}; // class gsIndexSubDomainFilterIterator


/**
   @brief Per-patch, ownership-filtered view of a gsIndexSubDomain.

   Returned by gsIndexSubDomain::subdomain(k): wraps the parent's raw
   per-patch domain and filters volume/boundary iteration down to the
   elements owned by the gsIndexSubDomain (see gsIndexSubDomainFilterIterator).
   Geometric queries that don't depend on ownership (degree, dim,
   boundingBox, mesh) are forwarded straight to the raw per-patch domain.

   \ingroup Domain
*/
template<class T>
class gsIndexSubDomainPatchView : public gsDomain<T>
{
    const gsIndexSubDomain<T>& m_owner;
    typename gsDomain<T>::Ptr  m_patchDomain; // raw, unfiltered per-patch domain

    // Lazily filled, unlike gsIndexSubDomain's own caches (Fix 1a): the
    // side/patch pair isn't known until this view is constructed, so
    // eager-for-every-side isn't practical here. INVARIANT (load-bearing,
    // not just documentation): these must only ever be first-touched from
    // serial code -- gsExprAssembler's per-(patch,side) warm-up loop before
    // #pragma omp parallel (see _warmupSubdomainViews(), Fix 1b) -- never
    // from inside the parallel region itself. Threads inside the parallel
    // region must only read through an already-warm cache.
    mutable std::map<short_t, size_t> m_bdrCountCache;
    mutable bool   m_numElementsCached = false;
    mutable size_t m_numElementsCache  = 0;

public:
    typedef typename gsDomain<T>::Ptr      Ptr;
    typedef typename gsDomain<T>::iterator iterator;

    gsIndexSubDomainPatchView(const gsIndexSubDomain<T>& owner, index_t patch)
    : m_owner(owner), m_patchDomain(owner.parentDomain().subdomain(patch))
    {
        this->setPatchIndex(patch);
    }

    iterator beginAll() const override
    {
        return iterator(new gsIndexSubDomainFilterIterator<T>(
            m_owner, m_patchDomain->beginAll(), m_patchDomain->endAll(), false));
    }

    iterator beginBdr(const boxSide bs = boundary::all) const override
    {
        return iterator(new gsIndexSubDomainFilterIterator<T>(
            m_owner, m_patchDomain->beginBdr(bs), m_patchDomain->endBdr(bs), true));
    }

    iterator endBdr(const boxSide bs = boundary::all) const override
    {
        return iterator(new gsDomainIteratorEnd<T>(this->numElementsBdr(bs)));
    }

    size_t numElements() const override
    {
        if (m_numElementsCached) return m_numElementsCache;
        size_t n = 0;
        auto it  = m_patchDomain->beginAll();
        auto end = m_patchDomain->endAll();
        // it.id() is local to the raw per-patch domain (see isOwned() in
        // gsIndexSubDomainFilterIterator); globalElementId() converts it to
        // the global id directly, no point location needed (Fix 3).
        for (; it != end; ++it)
            if (m_owner.contains(m_owner.parentDomain().globalElementId(
                    static_cast<index_t>(it.patchIndex()), static_cast<index_t>(it.id())))) ++n;
        m_numElementsCache = n;
        m_numElementsCached = true;
        return n;
    }

    size_t numElementsBdr(boxSide const & bs = boundary::all) const override
    {
        auto cacheIt = m_bdrCountCache.find(bs.index());
        if (cacheIt != m_bdrCountCache.end()) return cacheIt->second;
        size_t n = 0;
        auto it  = m_patchDomain->beginBdr(bs);
        auto end = m_patchDomain->endBdr(bs);
        for (; it != end; ++it)
            if (m_owner.contains(m_owner.volumeElementIdFromBoundary(it))) ++n;
        m_bdrCountCache[bs.index()] = n;
        return n;
    }

    short_t degree(short_t i = 0) const override { return m_patchDomain->degree(i); }
    short_t dim()                  const override { return m_patchDomain->dim(); }

    gsMatrix<T> boundingBox() const override { return m_patchDomain->boundingBox(); }
    gsMesh<T>   mesh()        const override { return m_patchDomain->mesh(); }

    std::ostream& print(std::ostream& os) const override
    {
        os << "gsIndexSubDomainPatchView on patch " << this->patchIndex();
        return os;
    }

}; // class gsIndexSubDomainPatchView


// -----------------------------------------------------------------------
// gsIndexSubDomainIterator out-of-line definitions
// (placed after gsIndexSubDomain so that subdomain.parentDomain() /
//  subdomain.elementIndices() are fully visible)
// -----------------------------------------------------------------------

template<class T>
gsIndexSubDomainIterator<T>::gsIndexSubDomainIterator(
        const gsIndexSubDomain<T>& subdomain)
: Base(0),
  m_subdomain(subdomain),
  m_indices(subdomain.elementIndices()),
  m_parentIt(subdomain.parentDomain().beginAll())
{
    // Establish invariant: m_parentIt must sit at m_indices[0], not at parent element 0.
    // Without this seek, every partition except the one starting at index 0 iterates
    // the wrong elements silently.
    if (!m_indices.empty())
        m_parentIt += static_cast<index_t>(m_indices[0]);
}

template<class T>
void gsIndexSubDomainIterator<T>::reset()
{
    this->resetId();
    m_parentIt = m_subdomain.parentDomain().beginAll();
    if (!m_indices.empty())
        m_parentIt += static_cast<index_t>(m_indices[0]);
}

} // namespace gismo
