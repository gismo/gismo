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
    explicit gsCompositeDomainIterator(index_t _id = 0);

    gsCompositeDomainIterator(const domainContainer& _dom);

    gsCompositeDomainIterator(const gsCompositeDomainIterator & other);
    domainIter clone() const override;

    virtual ~gsCompositeDomainIterator();

    index_t patch() const override;

    size_t id() const override;

    index_t subdomainIndex() const;

    virtual size_t localId() const override;

private:

    void next() override;

    /// \brief Proceeds to the next element (skipping \p increment elements).
    void next(index_t increment) override;

    virtual gsVector<T> lowerCorner() const override;

    virtual gsVector<T> upperCorner() const override;

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
    std::vector<Ptr> m_dependencies;
    const gsBoxTopology * m_topology;
    
public:

    /** @brief Empty constructor
     */
    gsCompositeDomain(index_t patchId = -1);

    /** @brief Constructor from a \ref gsMultiBasis.
     *
     * @param multiBasis The multiBasis from which the domains are extracted.
     */
    gsCompositeDomain(const gsMultiBasis<T> & multiBasis);

    /** @brief Constructor from a \ref gsMultipatch.
     *
     * @param mp The multipatch from which the domains are extracted.
     */
    gsCompositeDomain(const gsMultiPatch<T> & mp);

    /** @brief Add a domain to the composite domain
     *
     * @param domain The domain to add
     */
    void addDomain(Ptr domain);

    void addDomain(Ptr domain, index_t patchID);

    /// Keep a shared ownership of a domain alive for the lifetime of this composite
    void keepAlive(Ptr domain);

    Ptr subdomain(index_t k) const override;
    index_t getGlobalPatchId(index_t k) const;

    gsCompositeDomain<T> patchDomain(index_t patchId) const;

    size_t nPieces() const override;

    const domainContainer & subdomains() const;

    iterator beginAll() const  override;

    /// See \ref gsDomain.h for documentation.
    size_t numElements() const override;

    /** @brief Degree (maximum) of the domain
    */
    short_t degree(short_t i = 0) const override;

    /// See \ref gsDomain.h for documentation.
    short_t dim() const override;

    /// See \ref gsDomain.h for documentation.
    gsMatrix<T> boundingBox() const override;

    /// See \ref gsDomain.h for documentation.
    virtual gsMesh<T> mesh() const override;

    /// See \ref gsDomain.h for documentation.
    std::ostream &print(std::ostream &os) const override;

    const gsBoxTopology & topology() const;

    /// @brief Decompose a multi-patch domain into a specified number of subdomains.
    ///
    /// This method aims to decompose the current composite domain, which consists of multiple
    /// patches (subdomains), into a target number of `npieces`. The behavior depends on
    /// the relationship between the number of existing patches (`npatches`) and `npieces`.
    ///
    /// @param npieces The desired number of subdomains after decomposition.
    /// @return A shared pointer to the decomposed domain.
    virtual typename gsDomain<T>::Ptr decompose(size_t npieces) const override;

}; // class gsCompositeDomain

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCompositeDomain.hpp)
#endif