/** @file gsTensorSubDomain.hpp

    @brief Implementation of gsTensorSubDomain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsTensorSubDomain.h>
#include <gsDomain/gsKdNode.h>    // Needed for gsKdNode
#include <gsDomain/gsCompositeDomain.h>
#include <queue>

namespace gismo
{

template<short_t d,class T>
gsTensorSubDomain<d,T>::Range::Range(index_t s, index_t e) : start(s), end(e) { }

template<short_t d,class T>
index_t gsTensorSubDomain<d,T>::Range::size() const { return end - start; }

template<short_t d,class T>
gsTensorSubDomain<d,T>::gsTensorSubDomain(const gsTensorDomain<d,T>& parent,
                      const std::vector<Range>& ranges,
                      index_t patchId,
                      memory::shared_ptr<gsTensorDomain<d,T>> parentPtr)
    : gsSubDomain<T>(patchId), m_parent(parent), m_tensorParent(&parent), m_parentPtr(parentPtr), m_ranges(ranges)
{
    GISMO_ASSERT(ranges.size() == d, 
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

template<short_t d,class T>
typename gsTensorSubDomain<d,T>::domainIter gsTensorSubDomain<d,T>::beginAll() const
{
    return domainIter(new gsTensorSubDomainIterator<T,d>(*this));
}

template<short_t d,class T>
typename gsTensorSubDomain<d,T>::domainIter gsTensorSubDomain<d,T>::endAll() const
{
    return domainIter(new gsDomainIteratorEnd<T>(m_numElements));
}

template<short_t d,class T>
typename gsTensorSubDomain<d,T>::domainIter gsTensorSubDomain<d,T>::beginBdr(const boxSide bs) const
{
    GISMO_NO_IMPLEMENTATION
}

template<short_t d,class T>
size_t gsTensorSubDomain<d,T>::numElements() const { return m_numElements; }

template<short_t d,class T>
size_t gsTensorSubDomain<d,T>::numElementsBdr(const boxSide& s) const
{
    GISMO_NO_IMPLEMENTATION
}

template<short_t d,class T>
short_t gsTensorSubDomain<d,T>::dim() const { return d; }

template<short_t d,class T>
short_t gsTensorSubDomain<d,T>::degree(short_t i) const
{
    gsInfo<<m_parent;
    return m_parent.degree(i);
}

template<short_t d,class T>
gsMatrix<T> gsTensorSubDomain<d,T>::boundingBox() const
{
    gsMatrix<T> result(d, 2);
    // This is approximate - would need access to knot vectors to be exact
    return m_parent.boundingBox();
}

template<short_t d,class T>
gsMesh<T> gsTensorSubDomain<d,T>::mesh() const
{
    GISMO_NO_IMPLEMENTATION
}

template<short_t d,class T>
const gsDomain<T>& gsTensorSubDomain<d,T>::parentDomain() const { return m_parentPtr ? *m_parentPtr : m_parent; }

template<short_t d,class T>
const gsTensorDomain<d,T>* gsTensorSubDomain<d,T>::tensorParent() const { return m_tensorParent; }

template<short_t d,class T>
const std::vector<typename gsTensorSubDomain<d,T>::Range>& gsTensorSubDomain<d,T>::ranges() const { return m_ranges; }

template<short_t d,class T>
index_t gsTensorSubDomain<d,T>::tensorIndexToGlobal(const std::vector<index_t>& tensorIdx) const
{
    GISMO_ASSERT(tensorIdx.size() == d, "Tensor index dimension mismatch");
    index_t result = 0;
    index_t stride = 1;
    
    for (short_t dim = d - 1; dim >= 0; --dim)
    {
        result += (m_ranges[dim].start + tensorIdx[dim]) * stride;
        
        auto parentSize = m_parent.component(dim)->numElements();
        stride *= parentSize;
    }
    return result;
}

template<short_t d,class T>
bool gsTensorSubDomain<d,T>::contains(index_t elementIndex) const
{
    std::vector<index_t> tensorIdx(d);
    index_t remaining = elementIndex;
    
    for (short_t dim = d - 1; dim >= 0; --dim)
    {
        auto parentSize = m_parent.component(dim)->numElements();
        tensorIdx[dim] = remaining % parentSize;
        remaining /= parentSize;
        
        if (tensorIdx[dim] < m_ranges[dim].start || tensorIdx[dim] >= m_ranges[dim].end)
            return false;
    }
    return true;
}

template<short_t d,class T>
const std::vector<index_t>& gsTensorSubDomain<d,T>::elementIndices() const
{
    if (m_cachedIndices.empty())
        m_cachedIndices = computeElementIndices();
    return m_cachedIndices;
}

template<short_t d,class T>
typename gsDomain<T>::Ptr gsTensorSubDomain<d,T>::decompose(const gsTensorDomain<d,T>& domain,
                                   const std::vector<Range>& ranges,
                                   size_t npieces,
                                   index_t patchId,
                                   memory::shared_ptr<gsTensorDomain<d,T>> parentPtr)
{
    struct CompareNode
    {
        bool operator()(const gsKdNode<d,index_t>* a,
                        const gsKdNode<d,index_t>* b) const
        {
            typename gsKdNode<d,index_t>::point sizeA = a->uppCorner() - a->lowCorner();
            typename gsKdNode<d,index_t>::point sizeB = b->uppCorner() - b->lowCorner();
            return sizeA.prod() < sizeB.prod();
        }
    };

    auto composite = memory::make_shared(new gsCompositeDomain<T>());

    gsVector<index_t,d> low, upp;
    size_t totalElements = 1;
    for (short_t dim = 0; dim < d; ++dim) 
    {
        low[dim] = ranges[dim].start;
        upp[dim] = ranges[dim].end;
        totalElements *= (ranges[dim].end - ranges[dim].start);
    }

    if (npieces <= 0 || totalElements == 0) 
    {
        return composite;
    }

    if (npieces == 1) 
    {
        composite->addDomain(
            memory::make_shared(new gsTensorSubDomain<d,T>(domain, ranges, patchId, parentPtr)));
        return composite;
    }

    gsKdNode<d,index_t>* root = new gsKdNode<d,index_t>(low, upp);

    GISMO_ASSERT(npieces<=totalElements, "Requested more pieces than elements.");

    std::priority_queue<gsKdNode<d,index_t>*, std::vector<gsKdNode<d,index_t>*>, CompareNode> nodeQueue;
    nodeQueue.push(root);

    while (nodeQueue.size() < npieces) 
    {
        gsKdNode<d,index_t>* currentNode = nodeQueue.top();
        nodeQueue.pop();

        auto size = currentNode->uppCorner() - currentNode->lowCorner();
        int splitAxis = 0;
        if (d > 1) {
            splitAxis = (size[0] >= size[1]) ? 0 : 1;
            for(int i = 2; i < d; ++i) {
                if (size[i] > size[splitAxis]) {
                    splitAxis = i;
                }
            }
        }

        index_t midIndex = (currentNode->lowCorner()[splitAxis] + currentNode->uppCorner()[splitAxis]) / 2;

        if (midIndex == currentNode->lowCorner()[splitAxis] || midIndex == currentNode->uppCorner()[splitAxis]) {
            nodeQueue.push(currentNode);
            if (nodeQueue.size() == npieces) break;
            continue;
        }
        
        currentNode->split(splitAxis, midIndex);
        
        nodeQueue.push(currentNode->left);
        nodeQueue.push(currentNode->right);
    }

    std::vector<gsKdNode<d,index_t>*> leafNodes;
    while (!nodeQueue.empty()) 
    {
        leafNodes.push_back(nodeQueue.top());
        nodeQueue.pop();
    }

    for (const auto& leafNode : leafNodes) 
    {
        std::vector<typename gsTensorSubDomain<d,T>::Range> leafRanges;
        for (short_t dim = 0; dim < d; ++dim) {
            leafRanges.push_back(typename gsTensorSubDomain<d,T>::Range(
                leafNode->lowCorner()[dim], leafNode->uppCorner()[dim]));
        }
        auto subdomain = memory::make_shared(new gsTensorSubDomain<d,T>(domain, leafRanges, patchId, parentPtr));
        composite->addDomain(subdomain);
    }

    delete root;
    
    return composite;
}

template<short_t d,class T>
typename gsDomain<T>::Ptr gsTensorSubDomain<d,T>::decompose(size_t npieces) const
{
    return decompose(m_parent, m_ranges, npieces, this->patchId(), m_parentPtr);
}

template<short_t d,class T>
std::vector<index_t> gsTensorSubDomain<d,T>::computeElementIndices() const
{
    std::vector<index_t> result;
    result.reserve(m_numElements);
    
    std::vector<index_t> indices(d);
    for (short_t dim = 0; dim < d; ++dim)
        indices[dim] = m_ranges[dim].start;
    
    while (true) 
    {
        index_t globalIdx = tensorIndexToGlobal(indices);
        result.push_back(globalIdx);
        
        short_t dim = d - 1;
        while (dim >= 0) {
            indices[dim]++;
            if (indices[dim] < m_ranges[dim].end)
                break;
            indices[dim] = m_ranges[dim].start;
            dim--;
        }
        if (dim < 0) break;
    }
    
    return result;
}


template<class T,short_t d>
gsTensorSubDomainIterator<T,d>::gsTensorSubDomainIterator(const gsTensorSubDomain<d,T>& subdomain)
    : Base(subdomain), m_subdomain(subdomain), m_ranges(subdomain.ranges()),
      m_tensorIdx(d, 0), m_currentElement(0)
{
    for (short_t dim = 0; dim < d; ++dim)
        m_tensorIdx[dim] = m_ranges[dim].start;
}

template<class T,short_t d>
gsTensorSubDomainIterator<T,d>::gsTensorSubDomainIterator(const gsTensorSubDomainIterator& other) = default;

template<class T,short_t d>
typename gsTensorSubDomainIterator<T,d>::domainIterUPtr gsTensorSubDomainIterator<T,d>::clone() const
{
    return domainIterUPtr(new gsTensorSubDomainIterator(*this));
}

template<class T,short_t d>
gsTensorSubDomainIterator<T,d>::~gsTensorSubDomainIterator() { }

template<class T,short_t d>
bool gsTensorSubDomainIterator<T,d>::good() const
{
    return m_currentElement < (index_t)m_subdomain.numElements();
}

template<class T,short_t d>
index_t gsTensorSubDomainIterator<T,d>::numElements() const
{
    return m_subdomain.numElements();
}

template<class T,short_t d>
index_t gsTensorSubDomainIterator<T,d>::patch() const
{
    return Base::patch();
}

template<class T,short_t d>
size_t gsTensorSubDomainIterator<T,d>::localId() const
{
    return m_subdomain.tensorIndexToGlobal(m_tensorIdx);
}

template<class T,short_t d>
void gsTensorSubDomainIterator<T,d>::next()
{
    ++m_currentElement;
    
    if (!good())
        return;
        
    for (short_t dim = d - 1; dim >= 0; --dim)
    {
        ++m_tensorIdx[dim];
        if (m_tensorIdx[dim] < m_ranges[dim].end)
            break;
        m_tensorIdx[dim] = m_ranges[dim].start;
    }
}

template<class T,short_t d>
void gsTensorSubDomainIterator<T,d>::next(index_t increment)
{
    m_currentElement += increment;
    
    if (!good())
        return;
        
    index_t remaining = m_currentElement;
    std::vector<index_t> sizes(d);
    for (short_t dim = 0; dim < d; ++dim)
        sizes[dim] = m_ranges[dim].size();
    
    for (short_t dim = d - 1; dim >= 0; --dim)
    {
        m_tensorIdx[dim] = m_ranges[dim].start + (remaining % sizes[dim]);
        remaining /= sizes[dim];
    }
}

template<class T,short_t d>
gsVector<T> gsTensorSubDomainIterator<T,d>::lowerCorner() const
{
    gsVector<T> result(d);
    for (short_t dim = 0; dim < d; ++dim) {
        const gsTensorDomain<d,T>* parentTensorDomain = m_subdomain.tensorParent();
        GISMO_ASSERT(parentTensorDomain != nullptr, "Parent domain is not a gsTensorDomain!");

        const auto& knotVectorPtr = parentTensorDomain->component(dim);
        const gsKnotVector<T>* kv = dynamic_cast<const gsKnotVector<T>*>(knotVectorPtr.get());
        GISMO_ASSERT(kv != nullptr, "Component is not a gsKnotVector!");
        
        result[dim] = *(kv->domainUBegin() + m_tensorIdx[dim]);
    }
    return result;
}

template<class T,short_t d>
gsVector<T> gsTensorSubDomainIterator<T,d>::upperCorner() const
{
    gsVector<T> result(d);
    for (short_t dim = 0; dim < d; ++dim) {
        const gsTensorDomain<d,T>* parentTensorDomain = m_subdomain.tensorParent();
        GISMO_ASSERT(parentTensorDomain != nullptr, "Parent domain is not a gsTensorDomain!");

        const auto& knotVectorPtr = parentTensorDomain->component(dim);
        const gsKnotVector<T>* kv = dynamic_cast<const gsKnotVector<T>*>(knotVectorPtr.get());
        GISMO_ASSERT(kv != nullptr, "Component is not a gsKnotVector!");
        
        result[dim] = *(kv->domainUBegin() + m_tensorIdx[dim] + 1);
    }
    return result;
}

} // namespace gismo

