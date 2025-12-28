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
template<class T,short_t d> class gsTensorSubDomainIterator;

/**
    @brief A subdomain of a tensor product domain defined by axis-aligned ranges.
    
    This class defines a subdomain as a tensor product of ranges in each
    parametric direction. For example, in 2D this would be a rectangular region
    defined by [u_min, u_max] x [v_min, v_max] in parameter space.
    
    The element ranges are specified as knot line indices.
    
    \ingroup Domain
*/
template<short_t d,class T>
class gsTensorSubDomain : public gsSubDomain<T>
{
public:
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename gsKnotVector<T>::const_uiterator knotIter;
    
    // Forward declaration
    template<short_t d_, class T_>
    friend class gsTensorDomain;

    template<class T_, short_t d_>
    friend class gsTensorSubDomainIterator;
    
    /// Range for a single parametric direction [min_knot_index, max_knot_index)
    struct Range
    {
        index_t start;  ///< Starting knot line index
        index_t end;    ///< Ending knot line index (exclusive)
        
        Range(index_t s = 0, index_t e = 0);
        
        index_t size() const;
    };

    /// Constructor
    /// @param parent The parent gsTensorDomain
    /// @param ranges Vector of D ranges, one for each parametric direction
    /// @param patchId The global patch ID this subdomain represents
    gsTensorSubDomain(const gsTensorDomain<d,T>& parent,
                      const std::vector<Range>& ranges,
                      index_t patchId,
                      memory::shared_ptr<gsTensorDomain<d,T>> parentPtr = memory::shared_ptr<gsTensorDomain<d,T>>());

    virtual domainIter beginAll() const override;
    
    domainIter endAll() const override;

    domainIter beginBdr(const boxSide bs) const override;

    size_t numElements() const override;

    size_t numElementsBdr(const boxSide& s = boundary::none) const override;

    short_t dim() const override;

    short_t degree(short_t i) const override;
    
    gsMatrix<T> boundingBox() const override;

    gsMesh<T> mesh() const override;

    /// Return the parent domain
    const gsDomain<T>& parentDomain() const override;

    /// Return the tensor parent domain (must exist for tensor-based operations)
    const gsTensorDomain<d,T>* tensorParent() const;

    /// Get the ranges
    const std::vector<Range>& ranges() const;

    /// Convert tensor index (indices in each direction) to global element index
    index_t tensorIndexToGlobal(const std::vector<index_t>& tensorIdx) const;

    /// Check if element (by global index) belongs to this subdomain
    bool contains(index_t elementIndex) const override;

    /// Get all element indices in this subdomain
    const std::vector<index_t>& elementIndices() const override;

    static typename gismo::gsDomain<T>::Ptr decompose(const gismo::gsTensorDomain<d,T>& domain,
                                      const std::vector<Range>& ranges,
                                      size_t npieces,
                                      index_t patchId,
                                      memory::shared_ptr<gsTensorDomain<d,T>> parentPtr = memory::shared_ptr<gsTensorDomain<d,T>>());

    typename gismo::gsDomain<T>::Ptr decompose(size_t npieces) const override;

private:
    std::vector<index_t> computeElementIndices() const;

    const gismo::gsTensorDomain<d,T>& m_parent;
    const gsTensorDomain<d,T>* m_tensorParent;
    memory::shared_ptr<gsTensorDomain<d,T>> m_parentPtr;
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
template<class T,short_t d>
class gsTensorSubDomainIterator : public gsSubDomainIterator<T>
{
    typedef gsSubDomainIterator<T> Base;
    typedef typename gsDomainIterator<T>::uPtr domainIterUPtr;
    typedef typename gsTensorSubDomain<d,T>::Range Range;

public:
    explicit gsTensorSubDomainIterator(const gsTensorSubDomain<d,T>& subdomain);

    gsTensorSubDomainIterator(const gsTensorSubDomainIterator& other);

    domainIterUPtr clone() const override;

    virtual ~gsTensorSubDomainIterator();

    bool good() const;

    index_t numElements() const;

    index_t patch() const override;

private:

    size_t localId() const override;

    void next() override;

    void next(index_t increment) override;

    gsVector<T> lowerCorner() const override;

    gsVector<T> upperCorner() const override;

private:
    const gsTensorSubDomain<d,T>& m_subdomain;
    const std::vector<Range>& m_ranges;
    std::vector<index_t> m_tensorIdx;
    index_t m_currentElement;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorSubDomain.hpp)
#endif