/** @file gsIndexSubDomain.h

    @brief A subdomain defined by a set of element indices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsSubDomain.h>
#include <gsDomain/gsSubDomainIterator.h>

namespace gismo
{

// Forward declaration
template<class T> class gsIndexSubDomainIterator;

/**
    @brief A subdomain defined by an explicit set of element indices.
    
    This class creates a subdomain from a parent domain by specifying
    which elements (by their global indices in the parent) belong to
    this subdomain.
    
    The iterator for this subdomain filters the parent domain's iterator
    to only yield elements that are in the index set.
    
    \ingroup Domain
*/
template<class T>
class gsIndexSubDomain : public gsSubDomain<T>
{
public:
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef gsDomainIteratorWrapper<T> domainIter;

    /// Constructor
    /// @param parent The parent domain
    /// @param indices Vector of element indices in the parent domain
    gsIndexSubDomain(const gsDomain<T>& parent, 
                     const std::vector<index_t>& indices)
        : m_parent(parent), m_indices(indices)
    {
        // Sort and unique the indices
        std::sort(m_indices.begin(), m_indices.end());
        m_indices.erase(std::unique(m_indices.begin(), m_indices.end()), 
                        m_indices.end());
    }

    /// Constructor with move semantics for indices
    gsIndexSubDomain(const gsDomain<T>& parent, 
                     std::vector<index_t>&& indices)
        : m_parent(parent), m_indices(give(indices))
    {
        // Sort and unique the indices
        std::sort(m_indices.begin(), m_indices.end());
        m_indices.erase(std::unique(m_indices.begin(), m_indices.end()), 
                        m_indices.end());
    }

    virtual domainIter beginAll() const override;
    
    domainIter endAll() const override
    {
        return domainIter(new gsDomainIteratorEnd<T>(m_indices.size()));
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    size_t numElements() const override
    {
        return m_indices.size();
    }

    size_t numElementsBdr(const boxSide& s = boundary::none) const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    short_t dim() const override
    {
        return m_parent.dim();
    }

    short_t degree(short_t i) const override
    {
        return m_parent.degree(i);
    }

    gsMatrix<T> boundingBox() const override
    {
        return m_parent.boundingBox();
    }

    gsMesh<T> mesh() const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Return the parent domain
    const gsDomain<T>& parentDomain() const override { return m_parent; }

    unsigned getGlobalPatchId() const { return 0; } // Default for gsIndexSubDomain

    /// Get the element indices
    const std::vector<index_t>& elementIndices() const override { return m_indices; }

    /// Check if element is in this subdomain
    bool contains(index_t elementIndex) const override
    {
        return std::binary_search(m_indices.begin(), m_indices.end(), elementIndex);
    }

private:
    const gsDomain<T>& m_parent;
    std::vector<index_t> m_indices;
};


/**
    @brief Iterator for gsIndexSubDomain.
    
    Iterates through the parent domain but only yields elements that
    are in the index set of the subdomain.
    
    \ingroup Domain
*/
template<class T>
class gsIndexSubDomainIterator : public gsSubDomainIterator<T>
{
    typedef gsSubDomainIterator<T> Base;
    typedef typename gsDomainIterator<T>::uPtr domainIterUPtr;

public:
    explicit gsIndexSubDomainIterator(const gsIndexSubDomain<T>& subdomain)
        : Base(subdomain, 0), m_subdomain(subdomain),
          m_indices(subdomain.elementIndices()),
          m_current(0)
    {
    }

    gsIndexSubDomainIterator(const gsIndexSubDomainIterator& other) = default;

    domainIterUPtr clone() const override
    {
        return domainIterUPtr(new gsIndexSubDomainIterator(*this));
    }

    virtual ~gsIndexSubDomainIterator() { }

    bool good() const
    {
        return this->m_id < (index_t)m_indices.size();
    }

    index_t numElements() const
    {
        return m_indices.size();
    }

private:

    size_t localId() const override 
    { 
        return (this->m_id < (index_t)m_indices.size()) ? m_indices[this->m_id] : 0; 
    }

    void next() override
    {
        // Don't do anything here - let nextId() handle the increment
    }

    virtual gsVector<T> lowerCorner() const override
    {
        auto parent_it = m_subdomain.parentDomain().beginAll();
        parent_it += this->localId();
        return parent_it.lowerCorner();
    }

    virtual gsVector<T> upperCorner() const override
    {
        auto parent_it = m_subdomain.parentDomain().beginAll();
        parent_it += this->localId();
        return parent_it.upperCorner();
    }

private:
    const gsIndexSubDomain<T>& m_subdomain;
    const std::vector<index_t>& m_indices;
    index_t m_current;
};

// Implementation
template<class T>
typename gsIndexSubDomain<T>::domainIter
gsIndexSubDomain<T>::beginAll() const
{
    return domainIter(new gsIndexSubDomainIterator<T>(*this));
}

} // namespace gismo
