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
#include <gsDomain/gsDecompositionStrategy.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiPatch.h>

#include <gsUtils/gsCombinatorics.h>

#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

template<class T> class gsCompositeDomain;

/**
   @brief A domain typically coming from a multipatch basis/geometry
 */

template <class T>
class gsCompositeDomainIterator : public gsDomainIterator<T>
{
    typedef gsDomainIterator<T> Base;
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef std::vector<domainPtr> domainContainer;
    typedef typename gsDomainIterator<T>::uPtr domainIter;

    domainContainer m_domains;
    typename domainContainer::const_iterator m_curDomain;
    std::vector<size_t> m_offset; //offsets
    typename std::vector<size_t>::const_iterator m_curOffset;
    gsDomainIteratorWrapper<T> m_cur;

public:
    explicit gsCompositeDomainIterator(index_t _id = 0) : Base(_id) { }

    gsCompositeDomainIterator(const domainContainer& _dom)
    : Base(), m_domains(_dom)
    {
        GISMO_ASSERT(!m_domains.empty(), "Empty..");
        m_curDomain = m_domains.begin();

        m_offset.reserve(m_domains.size()+1);
        m_offset.push_back(0);
        for( auto & sd : m_domains )
            m_offset.push_back(m_offset.back()+sd->numElements());
        m_curOffset = m_offset.begin();
        m_cur = (*m_curDomain)->beginAll();
        this->m_pside.patch = m_cur.patch();
    }

    gsCompositeDomainIterator(const gsCompositeDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsCompositeDomainIterator(*this)); }

    virtual ~gsCompositeDomainIterator() { }

    index_t patch() const override
    {
        if (m_curDomain == m_domains.end())
            return -1;
        return m_cur.patch();
    }

    size_t id() const override
    {
        if (m_curDomain == m_domains.end())
            return m_offset.back();
        return *m_curOffset + m_cur.id();
    }

    index_t subdomainIndex() const
    {
        if (m_curDomain == m_domains.end())
            return -1;
        return std::distance(m_domains.begin(), m_curDomain);
    }

    virtual size_t localId() const override { return m_cur.localId(); }

private:

    void next() override
    {
        if (m_curDomain == m_domains.end()) return;
        
        ++m_cur; // Always increment the child iterator first
        
        // Check if the incremented iterator is still within the current child domain
        if (m_cur != (*m_curDomain)->endAll()) // Using comparison with end()
        {
            this->m_pside.patch = m_cur.patch(); // Update patch ID from the current element
        }
        else // m_cur is no longer good, meaning it reached the end of the child domain
        {
            // Advance to the next child domain
            ++m_curOffset;
            ++m_curDomain;
            if (m_curDomain != m_domains.end())
            {
                m_cur = (*m_curDomain)->beginAll(); // Reinitialize m_cur for the new child domain
                this->m_pside.patch = m_cur.patch(); // Update patch ID from the first element of the new child domain
            } else {
                this->m_pside.patch = -1; // End of iteration
            }
        }
    }

    /// \brief Proceeds to the next element (skipping \p increment elements).
    void next(index_t increment) override
    {
        GISMO_ASSERT(increment >= 0, "Negative increment is not supported.");
        if (increment == 0)
            return;
            
        if (m_curDomain == m_domains.end()) return;

        const size_t target_global_pos = (*m_curOffset + m_cur.id()) + increment;

        const size_t total_elements = m_offset.back();
        if (target_global_pos >= total_elements)
        {
            // Move iterator to the end
            m_curDomain = m_domains.end();
            m_curOffset = m_offset.end() - 1;
            this->m_pside.patch = -1; // End of iteration
            return;
        }

        // Find which subdomain the target position falls into
        // std::upper_bound takes a range [first, last)
        // Here, m_offset.end() should be correct as it includes the sentinel total_elements
        auto it_offset = std::upper_bound(m_offset.begin(), m_offset.end(), target_global_pos);
        --it_offset; // Decrement to get the start offset of the correct subdomain

        // Update current domain and offset
        m_curOffset = it_offset;
        m_curDomain = m_domains.begin() + std::distance(m_offset.begin(), it_offset);
        
        // Adjust local position within the current subdomain
        index_t local_pos = target_global_pos - *m_curOffset;
        m_cur = (*m_curDomain)->beginAll(); // Reinitialize m_cur for the new child domain
        m_cur += local_pos; // Advance m_cur to the exact local position

        this->m_pside.patch = m_cur.patch(); // Update patch ID after m_cur is fully advanced
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
    typedef std::vector<Ptr> domainContainer; //domain

    domainContainer m_domains;
    const gsBoxTopology * m_topology;
    
public:

    /** @brief Empty constructor
     */
    gsCompositeDomain(index_t patchId = -1)
    : Base(patchId), m_domains(), m_topology(nullptr)
    {
    }

    /** @brief Constructor from a \ref gsMultiBasis.
     *
     * @param multiBasis The multiBasis from which the domains are extracted.
     */
    gsCompositeDomain(const gsMultiBasis<T> & multiBasis)
    : Base(-1), m_domains(multiBasis.nPieces()), m_topology(&multiBasis.topology())
    {
        for (index_t i = 0; i != multiBasis.nPieces(); ++i)
        {
            m_domains[i] = multiBasis.basis(i).domain();
            m_domains[i]->setPatchId(i);
        }
    }

    /** @brief Constructor from a \ref gsMultipatch.
     *
     * @param mp The multipatch from which the domains are extracted.
     */
    gsCompositeDomain(const gsMultiPatch<T> & mp)
    : Base(-1), m_domains(mp.nPieces()), m_topology(&mp.topology())
    {
        for (index_t i = 0; i != mp.nPieces(); ++i)
        {
            m_domains[i] = mp.patch(i).basis().domain();
            m_domains[i]->setPatchId(i);
        }
    }

    /** @brief Add a domain to the composite domain
     *
     * @param domain The domain to add
     */
    void addDomain(Ptr domain)
    {
        m_domains.push_back(domain);
    }

    void addDomain(Ptr domain, index_t patchID)
    {
        domain->setPatchId(patchID);
        m_domains.push_back(domain);
    }

    Ptr subdomain(index_t k) const override { return m_domains[k]; }
    index_t getGlobalPatchId(index_t k) const { return m_domains[k]->patchId(); }

    gsCompositeDomain<T> patchDomain(index_t patchId) const
    {
        gsCompositeDomain<T> result;
        for (const auto& domain : m_domains)
        {
            if (domain->patchId() == patchId)
                result.addDomain(domain);
        }
        return result;
    }

    size_t nPieces() const override { return m_domains.size(); }

    const domainContainer & subdomains() const { return m_domains;}

    iterator beginAll() const  override { return iterator(new gsCompositeDomainIterator<T>(m_domains)); }

    /// See \ref gsDomain.h for documentation.
    size_t numElements() const override
    {
        size_t sz = 0;
        for (const auto & d : m_domains)
            sz += d->numElements();
        return sz;
    }

    /** @brief Degree (maximum) of the domain
    */
    short_t degree(short_t i = 0) const override
    {
        GISMO_ASSERT(!m_domains.empty(), "Empty composite domain.");
        short_t result = m_domains[0]->degree(i);
        for (size_t p = 1; p < m_domains.size(); ++p)
            if (m_domains[p]->degree(i) > result )
                result = m_domains[p]->degree(i);
        return result;
    }

    /// See \ref gsDomain.h for documentation.
    short_t dim() const override { return m_domains.front()->dim(); }

    /// See \ref gsDomain.h for documentation.
    gsMatrix<T> boundingBox() const override
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \ref gsDomain.h for documentation.
    virtual gsMesh<T> mesh() const override
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \ref gsDomain.h for documentation.
    std::ostream &print(std::ostream &os) const override
    {
        os << "Domain container with " << m_domains.size() << " domains.";
        return os;
    }

    const gsBoxTopology & topology() const { return *m_topology; }

    virtual typename gsDomain<T>::Ptr decompose(index_t npieces) const override
    {
        return this->decompose(npieces, decompositionStrategy::tensor);
    }

    /// @brief Decompose a multi-patch domain into a specified number of subdomains.
    ///
    /// This method aims to decompose the current composite domain, which consists of multiple
    /// patches (subdomains), into a target number of `npieces`. The behavior depends on
    /// the relationship between the number of existing patches (`npatches`) and `npieces`,
    /// as well as the chosen `decompositionStrategy`.
    ///
    /// @param npieces The desired number of subdomains after decomposition.
    /// @param strategy The decomposition strategy to use.
    ///                 - `decompositionStrategy::tensor`:
    ///                   - If `npatches == 0`: Returns an empty composite domain.
    ///                   - If `npatches == npieces`: Each original patch is returned as a subdomain.
    ///                   - If `npatches < npieces`: The `npieces` are distributed as evenly as possible
    ///                     among the original patches. Each original patch is then recursively
    ///                     decomposed into its assigned number of pieces using its own `decompose`
    ///                     method (e.g., `gsTensorDomain::decompose` for tensor patches).
    ///                   - If `npatches > npieces`: The original patches are merged into `npieces`
    ///                     new composite subdomains. Patches are assigned to the new subdomains
    ///                     to balance the total number of elements as evenly as possible.
    ///                 - `decompositionStrategy::localOptimalBalancing`:
    ///                   - If `npatches == 0` or `numElements() == 0`: Returns an empty composite domain.
    ///                   - Distributes `npieces` among original patches based on their proportional
    ///                     element count, aiming for local balance. Each original patch is then
    ///                     split into index subdomains (using `gsIndexSubDomain`) according to its
    ///                     assigned number of pieces. This strategy ensures each original patch
    ///                     contributes at least one piece, if possible, and aims to balance the workload
    ///                     across patches based on their size. If `npieces < npatches`, the effective
    ///                     number of resulting pieces will be `npatches` (one for each original patch),
    ///                     as pieces cannot be reduced below one per original patch in this strategy.
    ///                 - `decompositionStrategy::optimalBalancing`:
    ///                   - If `npatches == 0` or `numElements() == 0`: Returns an empty composite domain.
    ///                   - Caps `npieces` at the total number of elements if `npieces` is larger.
    ///                   - All elements from all original patches are collected and then distributed
    ///                     globally as evenly as possible into `npieces` resulting subdomains.
    ///                     Each resulting subdomain is a composite domain containing `gsIndexSubDomain`
    ///                     objects that reference specific elements from the original patches.
    ///                     This provides the most balanced element distribution globally.
    /// @return A shared pointer to a new `gsCompositeDomain` containing the decomposed subdomains.
    virtual typename gsDomain<T>::Ptr decompose(index_t npieces, decompositionStrategy strategy) const override;

}; // class gsCompositeDomain

} // namespace gismo

#include <gsDomain/gsCompositeDomain.hpp>
