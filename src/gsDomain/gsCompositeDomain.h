/** @file gsDomainContainer.h

    @brief Container for multiple domains.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomain.h>
#include <gsDomain/gsDomainIterator.h>

#include <gsUtils/gsCombinatorics.h>

#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

template<class T> class gsCompositeDomain;

/**
   @brief A domain typically coming from a multipatch basis/geometry


   patchIndex() : The index of the patch where we are at
   subdomainIndex() : The index of the subdomain where we are at
   id(): The global element index of the current element
   localId(): The index of the current element in its containing subdomain
   
   
   side(): The index of side where this elemement lies (if applicable)
 */
template <class T>
class gsCompositeDomainIterator : public gsDomainIterator<T>
{
    typedef gsDomainIterator<T> Base;
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef std::vector<domainPtr> domainContainer;
    typedef typename gsDomainIterator<T>::uPtr domainIter;

    domainContainer m_domains;
    std::vector<size_t> m_numElOffset; //offsets
    gsDomainIteratorWrapper<T> m_cur;

    index_t m_sid; //< Composite (sub-domain) id

public:
    explicit gsCompositeDomainIterator(index_t _id = 0) : Base(_id) { }

    gsCompositeDomainIterator(domainContainer _dom)
    : Base(), m_domains(give(_dom))
    {
        GISMO_ASSERT(!m_domains.empty(), "Empty..");
        m_numElOffset.reserve(m_domains.size()+1);
        m_numElOffset.push_back(0);
        m_sid = 0;
        for( auto & sd : m_domains )
        {
            m_numElOffset.push_back(m_numElOffset.back()+sd->numElements());
        }
        m_cur = m_domains.front()->beginAll();
    }

    gsCompositeDomainIterator(const gsCompositeDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsCompositeDomainIterator(*this)); }

    virtual ~gsCompositeDomainIterator() { }

    
    index_t subdomainIndex() const override { return m_sid; }

    virtual size_t localId() const { return m_cur.id(); }

    index_t patchIndex() const override { return m_cur.patchIndex(); }
    
private:

    void next() override
    {
        // previous
        //note: we cannot rely on this->id()
        //if (m_cur.id() + 1 == m_numEl[m_sid+1]-m_numEl[m_sid])

        if ( this->id() + 1 == m_numElOffset[m_sid+1])
        {
            ++m_sid;
            if ((size_t)m_sid<m_domains.size())
            {
                m_cur = m_domains[m_sid]->beginAll();
            }
            else
                return;
        }
        else
            ++m_cur;
        return;
    }

    /// \brief Proceeds to the next element (skipping \p increment elements).
    void next(index_t increment) override
    {
        const size_t pos = this->id() + increment;
        if ( pos < m_numElOffset[m_sid+1])
        {
            m_cur += increment;
            return;
        }

        //note: could end at min( --m_numElOffset.end(), m_numElOffset.begin() + increment)
        auto it = --std::upper_bound(m_numElOffset.begin()+m_sid, --m_numElOffset.end(), pos);
        m_sid = it - m_numElOffset.begin();
        m_cur  = m_domains[m_sid]->beginAll();
        m_cur +=  pos - m_numElOffset[m_sid];
        return;
    }

    virtual gsVector<T> lowerCorner() const override
    { return m_cur.lowerCorner(); }

    virtual gsVector<T> upperCorner() const override
    { return m_cur.upperCorner(); }

}; // class gsCompositeDomainIterator


template<class T>
class gsCompositeDomain : public gsDomain<T>
{
private:
    typedef gsDomain<T> Base;
    typedef typename Base::Ptr Ptr;
    typedef typename Base::iterator iterator;
    typedef std::vector<Ptr> domainContainer;

    domainContainer m_domains;

    const gsBoxTopology * m_topology;

    // Cumulative element count before patch i, i.e. the same offsets
    // gsCompositeDomainIterator recomputes on every beginAll() call
    // (m_numElOffset there) -- cached here once, eagerly, at construction
    // time so elementIndex()/globalElementId() below can convert a per-patch
    // local id into the global beginAll()-order id without re-walking every
    // patch's elements. Built eagerly (not lazily) so that subdomain()
    // callers can query it from inside an OpenMP parallel region without a
    // data race (see gsIndexSubDomain.h / gsExprAssembler.h Fix 1).
    std::vector<size_t> m_numElOffset;

    void buildNumElOffset()
    {
        m_numElOffset.assign(1, 0);
        for (const Ptr & sd : m_domains)
            m_numElOffset.push_back(m_numElOffset.back() + sd->numElements());
    }

public:

    /** @brief Constructor from a \ref gsMultiBasis.
     *
     * @param multiBasis The multiBasis from which the domains are extracted.
     */
    //gsCompositeDomain(const gsFunctionSet<T> & multiBasis)
    gsCompositeDomain(const gsMultiBasis<T> & multiBasis)
    : Base(), m_domains(multiBasis.nPieces()), m_topology(&multiBasis.topology())
    {
        for (index_t i = 0; i != multiBasis.nPieces(); ++i)
        {
            m_domains[i] = multiBasis.basis(i).domain();
            m_domains[i]->setPatchIndex(i);
        }
        buildNumElOffset();
    }

    /** @brief Constructor from a \ref gsMultipatch.
     *
     * @param mp The multipatch from which the domains are extracted.
     */
    //gsCompositeDomain(const gsFunctionSet<T> & multiBasis)
    gsCompositeDomain(const gsMultiPatch<T> & mp)
    : Base(), m_domains(mp.nPieces()), m_topology(&mp.topology())
    {
        for (index_t i = 0; i != mp.nPieces(); ++i)
        {
            m_domains[i] = mp.patch(i).basis().domain();
            m_domains[i]->setPatchIndex(i);
        }
        buildNumElOffset();
    }

    // Caller is responsible for tagging patch indices on the supplied domains
    // via setPatchIndex(); untagged domains default to patchIndex()==0.
    gsCompositeDomain(domainContainer domains)
    : Base(), m_domains(give(domains)) { buildNumElOffset(); }

    // void insert(Ptr other);

    Ptr subdomain(index_t k) const override { return m_domains[k]; }

    // returns the index of the first element on subdomain \a k in the global numbering
    //size_t offset(size_t k) const;

    size_t nPieces() const override { return m_domains.size(); }

    const domainContainer & subdomains() const { return m_domains;}

    iterator beginAll() const  override { return iterator(new gsCompositeDomainIterator<T>(m_domains)); }

    /// See \ref gsDomain.h for documentation.
    size_t numElements() const override
    {
        size_t sz = 0;
        for (size_t i = 0; i < m_domains.size(); ++i)
            sz += m_domains[i]->numElements();
        return sz;
    }

    /** @brief Degree (maximum) of the domain
    */
    short_t degree(short_t i = 0) const override
    {
        GISMO_ASSERT(m_domains.size(), "Empty composite domain.");
        short_t result = m_domains[0]->degree(i);
        for (size_t p = 0; p < m_domains.size(); ++p)
            if (m_domains[p]->degree(i) > result )
                result = m_domains[p]->degree(i);
        return result;
    }

    /// See \ref gsDomain.h for documentation.
    short_t dim() const override { return m_domains.front()->dim(); }

    /// Fast element lookup; see \ref gsDomain::elementIndex. Delegates to
    /// the per-patch domain's own elementIndex() for the local id, then
    /// converts it to the global beginAll()-order id via globalElementId()
    /// below. Returns -1 (falls back to the caller's own point-location) if
    /// \a patch's own domain has no fast lookup.
    index_t elementIndex(index_t patch, const gsVector<T>& u) const override
    {
        GISMO_ASSERT(patch>=0 && static_cast<size_t>(patch)<m_domains.size(),
                     "Patch index "<<patch<<" out of range.");
        const index_t local = m_domains[patch]->elementIndex(0, u);
        if (local < 0) return -1;
        return globalElementId(patch, local);
    }

    /// See \ref gsDomain::globalElementId. Converts a per-patch local
    /// (beginAll()-order) element id into the global beginAll()-order id of
    /// this composite domain, by adding the patch's cumulative-element
    /// offset (same convention gsCompositeDomainIterator uses to go the
    /// other way). \a m_numElOffset is built eagerly at construction (Fix
    /// 1), so this is a pure lookup -- safe to call from multiple threads.
    index_t globalElementId(index_t patch, index_t localId) const override
    {
        GISMO_ASSERT(patch>=0 && static_cast<size_t>(patch)<m_domains.size(),
                     "Patch index "<<patch<<" out of range.");
        return static_cast<index_t>(m_numElOffset[patch]) + localId;
    }

    /// See \ref gsDomain.h for documentation.
    gsMatrix<T> boundingBox() const override
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \ref gsDomain.h for documentation.
    virtual gsMesh<T> mesh() const override
    {
        // gsMesh<T> mesh;
        // mesh.setDimension(d);
        // mesh.setBasis(m_basis);
        // return mesh;
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \ref gsDomain.h for documentation.
    std::ostream &print(std::ostream &os) const override
    {
        os << "Domain container with " << m_domains.size() << " domains.";
        return os;
    }

    const gsBoxTopology & topology() const { return *m_topology; }

}; // class gsCompositeDomain


} // namespace gismo
