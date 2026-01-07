/** @file gsHSubDomain.h

    @brief A subdomain of a hierarchical domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsSubDomain.h>
#include <gsDomain/gsSubDomainIterator.h>
#include <gsDomain/gsTensorSubDomain.h>
#include <gsDomain/gsCompositeDomain.h>
#include <gsDomain/gsHDomain.h>

namespace gismo
{

// Forward declaration
template<short_t d, class T, class Z> class gsHSubDomainIterator;

/**
    @brief A subdomain of a hierarchical domain at a specific level.
    
    A hierarchical subdomain is defined by selecting a specific level
    from a gsHDomain and a range of tensor indices at that level.
    The subdomain is implemented using a gsTensorSubDomain constructed
    on the level's tensor domain.
    
    \ingroup HSplines
*/
template<short_t d, class T, class Z = index_t>
class gsHSubDomain : public gsSubDomain<T>
{
public:
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef gsHDomain<d, T, Z> hDomain;
    typedef gsHTree<d, Z> hTree;
    typedef typename hTree::const_literator leafIterator;

    /// Constructor selecting a specific level and tensor range
    /// @param parent The parent gsHDomain
    /// @param basis The hierarchical tensor basis
    /// @param level The hierarchical level (0 = coarsest)
    /// @param ranges Vector of pairs (start, end) for tensor indices in each dimension
    gsHSubDomain(const gsHDomain<d, T, Z>& parent, 
                 const gsHTensorBasis<d,T>& basis,
                 short_t level,
                 const std::vector<std::pair<index_t, index_t>>& ranges)
        : m_parent(parent), m_basis(basis), m_level(level), m_ranges(ranges),
          m_leafIndices()
    {
        GISMO_ASSERT(level >= 0 && level < (short_t)m_basis.nLevels(),
                     "Level " << level << " out of range [0," << m_basis.nLevels()-1 << "]");
        GISMO_ASSERT(ranges.size() == (size_t)d,
                     "Number of ranges must equal dimension " << d);
        
        constructTensorSubDomains();
    }

    /// Constructor selecting a specific level and specific leaf indices
    /// @param parent The parent gsHDomain
    /// @param basis The hierarchical tensor basis
    /// @param level The hierarchical level
    /// @param leafIndices Indices of leaves at this level to include
    gsHSubDomain(const gsHDomain<d, T, Z>& parent, 
                 const gsHTensorBasis<d,T>& basis,
                 short_t level,
                 const std::vector<index_t>& leafIndices)
        : m_parent(parent), m_basis(basis), m_level(level), m_ranges(),
          m_leafIndices(leafIndices)
    {
        GISMO_ASSERT(level >= 0 && level < (short_t)m_basis.nLevels(),
                     "Level " << level << " out of range [0," << m_basis.nLevels()-1 << "]");
        
        constructTensorSubDomainsFromLeaves();
    }

    virtual domainIter beginAll() const override;

    domainIter endAll() const override
    {
        if (!m_tensorSubDomains.empty())
            return m_tensorSubDomains[0]->endAll();
        return domainIter(new gsDomainIteratorEnd<T>(0));
    }

    size_t numElements() const override
    {
        size_t total = 0;
        for (const auto& sd : m_tensorSubDomains)
            total += sd->numElements();
        return total;
    }

    size_t numElementsBdr(const boxSide& s = boundary::none) const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    short_t dim() const override { return d; }

    short_t degree(short_t i) const override
    {
        return m_basis.degree(i);
    }

    gsMatrix<T> boundingBox() const override
    {
        return m_basis.support();
    }

    gsMesh<T> mesh() const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Return the parent domain
    const gsDomain<T>& parentDomain() const override { return m_parent; }

    /// Get the hierarchical level
    short_t level() const { return m_level; }

    /// Get the tensor subdomains
    const std::vector<domainPtr>& tensorSubDomains() const { return m_tensorSubDomains; }

private:

    void constructTensorSubDomains()
    {
        // Collect all leaves at this level
        leafIterator it = m_parent.tree().beginLeafIterator();
        std::vector<std::pair<leafIterator, std::vector<std::pair<index_t, index_t>>>> selectedLeaves;
        
        while (it.good())
        {
            if (it.level() == m_level)
            {
                // Check if this leaf intersects with our ranges
                bool intersects = true;
                std::vector<std::pair<index_t, index_t>> leafRanges(d);
                
                for (short_t i = 0; i < d; ++i)
                {
                    index_t ll = it.lowerCorner()[i];
                    index_t uu = it.upperCorner()[i];
                    
                    if (m_basis.manualLevels())
                    {
                        m_basis._diadicIndexToKnotIndex(m_level, i, ll);
                        m_basis._diadicIndexToKnotIndex(m_level, i, uu);
                    }
                    
                    // Check intersection
                    if (uu <= m_ranges[i].first || ll >= m_ranges[i].second)
                    {
                        intersects = false;
                        break;
                    }
                    
                    // Clamp to range
                    leafRanges[i] = {
                        math::max(ll, m_ranges[i].first),
                        math::min(uu, m_ranges[i].second)
                    };
                }
                
                if (intersects)
                {
                    selectedLeaves.push_back({it, leafRanges});
                }
            }
            it.next();
        }
        
        // Create tensor subdomains for selected leaves
        // Get tensor domain at this level
        const auto& tensorBasis = m_basis.tensorLevel(m_level);
        auto tensorDomain = tensorBasis.domain();
        
        for (const auto& leafEntry : selectedLeaves)
        {
            const auto& leaf = leafEntry.first;
            const auto& leafRanges = leafEntry.second;
            auto subdomain = gismo::memory::make_shared<gsTensorSubDomain<d,T>>(
                *tensorDomain, leafRanges);
            m_tensorSubDomains.push_back(subdomain);
        }
    }

    void constructTensorSubDomainsFromLeaves()
    {
        // Collect all leaves at this level
        leafIterator it = m_parent.tree().beginLeafIterator();
        index_t leafCounter = 0;
        const auto& tensorBasis = m_basis.tensorLevel(m_level);
        auto tensorDomain = tensorBasis.domain();
        
        while (it.good())
        {
            if (it.level() == m_level)
            {
                // Check if this leaf is in our selection
                if (std::find(m_leafIndices.begin(), m_leafIndices.end(), leafCounter) 
                    != m_leafIndices.end())
                {
                    std::vector<std::pair<index_t, index_t>> leafRanges(d);
                    
                    for (short_t i = 0; i < d; ++i)
                    {
                        index_t ll = it.lowerCorner()[i];
                        index_t uu = it.upperCorner()[i];
                        
                        if (m_basis.manualLevels())
                        {
                            m_basis._diadicIndexToKnotIndex(m_level, i, ll);
                            m_basis._diadicIndexToKnotIndex(m_level, i, uu);
                        }
                        
                        leafRanges[i] = {ll, uu};
                    }
                    
                    auto subdomain = gismo::memory::make_shared<gsTensorSubDomain<d,T>>(
                        *tensorDomain, leafRanges);
                    m_tensorSubDomains.push_back(subdomain);
                }
                leafCounter++;
            }
            it.next();
        }
    }

private:
    const gsHDomain<d, T, Z>& m_parent;
    const gsHTensorBasis<d,T>& m_basis;
    short_t m_level;
    std::vector<std::pair<index_t, index_t>> m_ranges;
    std::vector<index_t> m_leafIndices;
    std::vector<domainPtr> m_tensorSubDomains;
};


/**
    @brief Iterator for hierarchical subdomains.
    
    Iterates through elements at a specific hierarchical level
    within the given tensor ranges.
    
    \ingroup HSplines
*/
template<short_t d, class T, class Z = long long>
class gsHSubDomainIterator : public gsSubDomainIterator<T>
{
    typedef gsSubDomainIterator<T> Base;
    typedef typename gsDomainIterator<T>::uPtr domainIterUPtr;

public:
    explicit gsHSubDomainIterator(const gsHSubDomain<d, T, Z>& subdomain)
        : Base(0), m_subdomains(subdomain.tensorSubDomains()), 
          m_currentIdx(0)
    {
        if (!m_subdomains.empty())
            m_currentIter = m_subdomains[0]->beginAll();
    }

    gsHSubDomainIterator(const gsHSubDomainIterator& other) = default;

    domainIterUPtr clone() const override
    {
        return domainIterUPtr(new gsHSubDomainIterator(*this));
    }

    virtual ~gsHSubDomainIterator() { }

private:

    size_t localId() const override
    {
        return m_currentIter.id();
    }

    void next() override
    {
        ++m_currentIter;
        if (!m_currentIter.good() && m_currentIdx + 1 < (index_t)m_subdomains.size())
        {
            ++m_currentIdx;
            m_currentIter = m_subdomains[m_currentIdx]->beginAll();
        }
    }

    void next(index_t increment) override
    {
        for (index_t i = 0; i < increment; ++i)
            next();
    }

private:
    std::vector<typename gsDomain<T>::Ptr> m_subdomains;
    index_t m_currentIdx;
    gsDomainIteratorWrapper<T> m_currentIter;
};

// Implementation
template<short_t d, class T, class Z>
typename gsHSubDomain<d, T, Z>::domainIter
gsHSubDomain<d, T, Z>::beginAll() const
{
    return domainIter(new gsHSubDomainIterator<d, T, Z>(*this));
}

} // namespace gismo

