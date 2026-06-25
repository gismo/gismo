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

namespace gismo
{

// Forward declaration needed by gsIndexSubDomainIterator
template<class T> class gsIndexSubDomain;


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
   auto sub = memory::make_shared<gsIndexSubDomain<T>>(dom, myIndices);
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
    }

    const gsDomain<T>& parentDomain() const override { return *m_parentPtr; }

    const std::vector<index_t>& elementIndices() const override { return m_indices; }

    bool contains(index_t elementIndex) const override
    {
        return std::binary_search(m_indices.begin(), m_indices.end(), elementIndex);
    }

    size_t numElements() const override { return m_indices.size(); }

    // Quadrature rule selection in gsExprEvaluator queries subdomain(patchIndex)
    // and nPieces().  Delegate both to the parent so multi-patch evaluations
    // (patchIndex > 0) don't hit the base-class assert.
    typename gsDomain<T>::Ptr subdomain(index_t k) const override
    { return m_parentPtr->subdomain(k); }

    size_t nPieces() const override { return m_parentPtr->nPieces(); }

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
