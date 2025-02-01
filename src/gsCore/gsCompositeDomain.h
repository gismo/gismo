/** @file gsDomainContainer.h

    @brief Container for multiple domains.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomain.h>
#include <gsCore/gsDomainIterator.h>

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
    std::vector<size_t> m_numEl;
    size_t m_index;
    gsDomainIteratorWrapper<T> m_cur;
    typedef gsDomainIterator<T> Base;
    using Base::m_domain;

public:
    explicit gsCompositeDomainIterator(index_t _id = 0) : Base(_id), m_index(0) { }

    gsCompositeDomainIterator( const gsCompositeDomain<T>& _dom)
    : Base(_dom), m_index(1)
    {
        m_numEl.reserve(_dom.numPieces()+1);
        m_numEl.push_back(0);
        for( auto & sd : _dom.subdomains() )
            m_numEl.push_back(m_numEl.back()+sd->numElements());
        m_cur = _dom.subdomain(0)->beginAll();
    }

    virtual ~gsCompositeDomainIterator() { }

private:

    bool next()
    {
        ++m_cur;
        if (this->id() + 1 == m_numEl[m_index])
        {
            ++m_index;
            if (m_index<m_numEl.size())
            {
                m_cur = m_domain->subdomain(m_index-1)->beginAll();
                m_cur.get()->setPatch(m_index-1);
            }
            else
                return false;
        }
        return true;
    }

    /// \brief Proceeds to the next element (skipping \p increment elements).
    bool next(index_t increment)
    {
        const size_t pos = increment + this->id();
        if ( pos < m_numEl[m_index])
        {
            m_cur += increment;
            return true;
        }

        auto it = std::lower_bound(m_numEl.begin()+m_index, m_numEl.end(), pos);
        m_index = it - m_numEl.begin();
        m_cur = m_domain->subdomain(m_index-1)->beginAll();
        m_cur +=  pos - m_numEl[m_index];
        return true;
    }

    virtual gsVector<T> lowerCorner() const
    { return m_cur.lowerCorner(); }

    virtual gsVector<T> upperCorner() const
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

public:

    /** @brief Constructor from a \ref gsMultiBasis.
     *
     * @param multiBasis The multiBasis from which the domains are extracted.
     */
    //gsCompositeDomain(const gsFunctionSet<T> & multiBasis)
    gsCompositeDomain(const gsMultiBasis<T> & multiBasis)
    : Base(), m_domains(multiBasis.size()), m_topology(&multiBasis.topology())
    {
        for (index_t i = 0; i != multiBasis.nPieces(); ++i)
            m_domains[i] = multiBasis.basis(i).domain();
    }

    Ptr subdomain(index_t k) const { return m_domains[k]; }

    // returns the index of the first element on subdomain \a k in the global numbering
    //size_t offset(size_t k) const;
                
    size_t numPieces() const { return m_domains.size(); }

    const domainContainer & subdomains() const { return m_domains;}

    iterator beginAll() const { return new gsCompositeDomainIterator<T>(*this); }

    /// See \ref gsDomain.h for documentation.
    size_t numElements() const override
    {
        size_t size = 0;
        for (size_t i = 0; i < m_domains.size(); ++i)
            size += m_domains[i]->numElements();
        return size;
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

}; // class gsCompositeDomain


} // namespace gismo
