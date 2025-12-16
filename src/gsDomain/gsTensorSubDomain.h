/** @file gsTensorSubDomain.h

    @brief A subdomain defined by tensor-product rectangular ranges.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsSubDomain.h>
#include <gsDomain/gsSubDomainIterator.h>
#include <gsDomain/gsTensorDomain.h>

namespace gismo
{

// Forward declaration
template<class T, int D> class gsTensorSubDomainIterator;

/**
    @brief A subdomain of a tensor product domain defined by axis-aligned ranges.
    
    This class defines a subdomain as a tensor product of ranges in each
    parametric direction. For example, in 2D this would be a rectangular region
    defined by [u_min, u_max] x [v_min, v_max] in parameter space.
    
    The element ranges are specified as knot line indices.
    
    \ingroup Domain
*/
template<class T, int D>
class gsTensorSubDomain : public gsSubDomain<T>
{
public:
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename gsKnotVector<T>::const_uiterator knotIter;
    
    // Forward declaration
    template<class T2, int D2>
    friend class gsTensorSubDomainIterator;
    
    /// Range for a single parametric direction [min_knot_index, max_knot_index)
    struct Range
    {
        index_t start;  ///< Starting knot line index
        index_t end;    ///< Ending knot line index (exclusive)
        
        Range(index_t s = 0, index_t e = 0) : start(s), end(e) { }
        
        index_t size() const { return end - start; }
    };

    /// Constructor
    /// @param parent The parent gsTensorDomain
    /// @param ranges Vector of D ranges, one for each parametric direction
    /// @param patchId The global patch ID this subdomain represents
    gsTensorSubDomain(const gsTensorDomain<T, D>& parent,
                      const std::vector<Range>& ranges,
                      index_t patchId)
        : gsSubDomain<T>(patchId), m_parent(parent), m_ranges(ranges)
    {
        GISMO_ASSERT(ranges.size() == D, 
                     "Number of ranges must match domain dimension");
        gsInfo<<"Made subdomain with ranges: ";
        for (const auto& r : m_ranges)
            gsInfo<<"["<<r.start<<","<<r.end<<") ";
        gsInfo<<"\n";
        // Compute number of elements
        m_numElements = 1;
        for (const auto& r : m_ranges)
            m_numElements *= r.size();
    }

    virtual domainIter beginAll() const override;
    
    domainIter endAll() const override
    {
        return domainIter(new gsDomainIteratorEnd<T>(m_numElements));
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    size_t numElements() const override { return m_numElements; }

    size_t numElementsBdr(const boxSide& s = boundary::none) const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    short_t dim() const override { return D; }

    short_t degree(short_t i) const override
    {
        gsInfo<<m_parent;
        return m_parent.degree(i);
    }
    
    gsMatrix<T> boundingBox() const override
    {
        gsMatrix<T> result(D, 2);
        // This is approximate - would need access to knot vectors to be exact
        return m_parent.boundingBox();
    }

    gsMesh<T> mesh() const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Return the parent domain
    const gsDomain<T>& parentDomain() const override { return m_parent; }

    /// Get the ranges
    const std::vector<Range>& ranges() const { return m_ranges; }

    /// Convert tensor index (indices in each direction) to global element index
    index_t tensorIndexToGlobal(const std::vector<index_t>& tensorIdx) const
    {
        GISMO_ASSERT(tensorIdx.size() == D, "Tensor index dimension mismatch");
        index_t result = 0;
        index_t stride = 1;
        
        for (int d = D - 1; d >= 0; --d)
        {
            result += (m_ranges[d].start + tensorIdx[d]) * stride;
            
            // Get size of this direction in parent
            auto parentSize = m_parent.component(d)->numElements();
            stride *= parentSize;
        }
        return result;
    }

    /// Check if element (by global index) belongs to this subdomain
    bool contains(index_t elementIndex) const override
    {
        // Convert global index to tensor indices
        std::vector<index_t> tensorIdx(D);
        index_t remaining = elementIndex;
        
        for (int d = D - 1; d >= 0; --d)
        {
            auto parentSize = m_parent.component(d)->numElements();
            tensorIdx[d] = remaining % parentSize;
            remaining /= parentSize;
            
            // Check if in range
            if (tensorIdx[d] < m_ranges[d].start || tensorIdx[d] >= m_ranges[d].end)
                return false;
        }
        return true;
    }

    /// Get all element indices in this subdomain
    const std::vector<index_t>& elementIndices() const override
    {
        if (m_cachedIndices.empty())
            m_cachedIndices = computeElementIndices();
        return m_cachedIndices;
    }

private:
    std::vector<index_t> computeElementIndices() const
    {
        std::vector<index_t> result;
        result.reserve(m_numElements);
        
        // Create a multi-index iterator over the ranges
        // Iterate through all combinations of the tensor indices
        std::vector<index_t> indices(D);
        for (int d = 0; d < D; ++d)
            indices[d] = m_ranges[d].start;
        
        while (true) {
            // Compute global index from tensor indices
            index_t globalIdx = tensorIndexToGlobal(indices);
            result.push_back(globalIdx);
            
            // Increment the multi-index
            int d = D - 1;
            while (d >= 0) {
                indices[d]++;
                if (indices[d] < m_ranges[d].end)
                    break;
                indices[d] = m_ranges[d].start;
                d--;
            }
            if (d < 0) break; // We've wrapped around all dimensions
        }
        
        return result;
    }

    const gsTensorDomain<T, D>& m_parent;
    std::vector<Range> m_ranges;
    size_t m_numElements;
    mutable std::vector<index_t> m_cachedIndices;
};


/**
    @brief Fast iterator for gsTensorSubDomain exploiting tensor structure.
    
    This iterator efficiently traverses a tensor subdomain by iterating
    through the ranges in each direction without checking all elements.
    
    \ingroup Domain
*/
template<class T, int D>
class gsTensorSubDomainIterator : public gsSubDomainIterator<T>
{
    typedef gsSubDomainIterator<T> Base;
    typedef typename gsDomainIterator<T>::uPtr domainIterUPtr;
    typedef typename gsTensorSubDomain<T, D>::Range Range;

public:
    explicit gsTensorSubDomainIterator(const gsTensorSubDomain<T, D>& subdomain)
        : Base(subdomain), m_subdomain(subdomain), m_ranges(subdomain.ranges()),
          m_tensorIdx(D, 0), m_currentElement(0)
    {
        // Initialize tensor indices to the start of each range
        for (int d = 0; d < D; ++d)
            m_tensorIdx[d] = m_ranges[d].start;
    }

    gsTensorSubDomainIterator(const gsTensorSubDomainIterator& other) = default;

    domainIterUPtr clone() const override
    {
        return domainIterUPtr(new gsTensorSubDomainIterator(*this));
    }

    virtual ~gsTensorSubDomainIterator() { }

    bool good() const
    {
        return m_currentElement < (index_t)m_subdomain.numElements();
    }

    index_t numElements() const
    {
        return m_subdomain.numElements();
    }

    index_t patch() const override
    {
        return Base::patch(); // Delegate to base class implementation
    }

private:

    size_t localId() const override
    {
        return m_subdomain.tensorIndexToGlobal(m_tensorIdx);
    }

    void next() override
    {
        ++m_currentElement;
        
        if (!good())
            return;
            
        // Increment the rightmost (fastest varying) index
        for (int d = D - 1; d >= 0; --d)
        {
            ++m_tensorIdx[d];
            if (m_tensorIdx[d] < m_ranges[d].end)
                break;
            // Reset and carry to next dimension
            m_tensorIdx[d] = m_ranges[d].start;
        }
    }

    void next(index_t increment) override
    {
        m_currentElement += increment;
        
        if (!good())
            return;
            
        // Convert current element counter to tensor indices
        index_t remaining = m_currentElement;
        std::vector<index_t> sizes(D);
        for (int d = 0; d < D; ++d)
            sizes[d] = m_ranges[d].size();
        
        for (int d = D - 1; d >= 0; --d)
        {
            m_tensorIdx[d] = m_ranges[d].start + (remaining % sizes[d]);
            remaining /= sizes[d];
        }
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> result(D);
        for (short_t dim = 0; dim < D; ++dim) {
            const gsTensorDomain<T, D>* parentTensorDomain = 
                dynamic_cast<const gsTensorDomain<T, D>*>(&m_subdomain.parentDomain());
            GISMO_ASSERT(parentTensorDomain != nullptr, "Parent domain is not a gsTensorDomain!");

            const auto& knotVectorPtr = parentTensorDomain->component(dim);
            const gsKnotVector<T>* kv = dynamic_cast<const gsKnotVector<T>*>(knotVectorPtr.get());
            GISMO_ASSERT(kv != nullptr, "Component is not a gsKnotVector!");
            
            // Get the lower bound of the element from the knot vector
            // m_tensorIdx[dim] is the 0-based element index within the relevant range
            result[dim] = *(kv->domainUBegin() + m_tensorIdx[dim]);
        }
        return result;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> result(D);
        for (short_t dim = 0; dim < D; ++dim) {
            const gsTensorDomain<T, D>* parentTensorDomain = 
                dynamic_cast<const gsTensorDomain<T, D>*>(&m_subdomain.parentDomain());
            GISMO_ASSERT(parentTensorDomain != nullptr, "Parent domain is not a gsTensorDomain!");

            const auto& knotVectorPtr = parentTensorDomain->component(dim);
            const gsKnotVector<T>* kv = dynamic_cast<const gsKnotVector<T>*>(knotVectorPtr.get());
            GISMO_ASSERT(kv != nullptr, "Component is not a gsKnotVector!");
            
            // Get the upper bound of the element from the knot vector
            // m_tensorIdx[dim] is the 0-based element index within the relevant range
            result[dim] = *(kv->domainUBegin() + m_tensorIdx[dim] + 1);
        }
        return result;
    }

private:
    const gsTensorSubDomain<T, D>& m_subdomain;
    const std::vector<Range>& m_ranges;
    std::vector<index_t> m_tensorIdx;
    index_t m_currentElement;
};

// Implementation
template<class T, int D>
typename gsTensorSubDomain<T, D>::domainIter
gsTensorSubDomain<T, D>::beginAll() const
{
    return domainIter(new gsTensorSubDomainIterator<T, D>(*this));
}

} // namespace gismo