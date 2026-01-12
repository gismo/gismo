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
gsTensorSubDomain<d,T>::gsTensorSubDomain(const gsTensorDomain<d,T>& parent,
                        const gsVector<index_t> & start,
                        const gsVector<index_t> & end,
                        index_t patchId,
                        memory::shared_ptr<gsTensorDomain<d,T>> parentPtr)
    : gsSubDomain<T>(patchId < 0 ? 0 : patchId), m_parent(parent), m_tensorParent(&parent), m_parentPtr(parentPtr), m_start(start), m_end(end)
{
    GISMO_ASSERT(start.rows() == d, "Start vector dimension mismatch");
    GISMO_ASSERT(end.rows() == d, "End vector dimension mismatch");
}

template<short_t d,class T>
typename gsTensorSubDomain<d,T>::domainIter gsTensorSubDomain<d,T>::beginAll() const
{
    return domainIter(new gsTensorSubDomainIterator<T,d>(*this));
}

template<short_t d,class T>
typename gsTensorSubDomain<d,T>::domainIter gsTensorSubDomain<d,T>::endAll() const
{
    return domainIter(new gsDomainIteratorEnd<T>(this->numElements()));
}

template<short_t d,class T>
typename gsTensorSubDomain<d,T>::domainIter gsTensorSubDomain<d,T>::beginBdr(const boxSide bs) const
{
    GISMO_NO_IMPLEMENTATION
}

template<short_t d,class T>
size_t gsTensorSubDomain<d,T>::numElements() const
{
    // Compute number of elements
    size_t numElements = 1;
    for (short_t i = 0; i < d; ++i)
    {
        GISMO_ASSERT(m_start[i] >= 0 && static_cast<size_t>(m_end[i]) <= m_parent.component(i)->numElements(),
                     "Subdomain range out of bounds");
        GISMO_ASSERT(m_start[i] < m_end[i], "Invalid subdomain range");
        numElements *= (static_cast<size_t>(m_end[i]) - static_cast<size_t>(m_start[i]));
    }
    return numElements;
}

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
const gsVector<index_t> & gsTensorSubDomain<d,T>::start() const { return m_start; }

template<short_t d,class T>
const gsVector<index_t> & gsTensorSubDomain<d,T>::end() const { return m_end; }

template<short_t d,class T>
index_t gsTensorSubDomain<d,T>::tensorIndexToGlobal(const gsVector<index_t>& tensorIdx) const
{
    GISMO_ASSERT(tensorIdx.size() == d, "Tensor index dimension mismatch");
    // Unpack lexicographical index
    index_t globalIndex = 0;
    index_t stride = 1;
    for (short_t dim = 0; dim < d; ++dim)
    {
        index_t idxInDim = tensorIdx[dim];
        GISMO_ASSERT(idxInDim >= m_start[dim] && idxInDim <= m_end[dim],
                     "Tensor index out of subdomain range");
        globalIndex += idxInDim * stride;
        stride *= m_parent.component(dim)->numElements();
    }
    return globalIndex;
}

template<short_t d,class T>
bool gsTensorSubDomain<d,T>::contains(index_t elementIndex) const
{
    gsVector<index_t> tensorIdx(d);
    // Convert global index to tensor index (unpack lexicographical index)
    index_t remainder = elementIndex;
    for (short_t dim = 0; dim < d; ++dim)
    {
        index_t dimSize = m_parent.component(dim)->numElements();
        tensorIdx[dim] = remainder % dimSize;
        remainder /= dimSize;

        if (tensorIdx[dim] < m_start[dim] || tensorIdx[dim] >= m_end[dim])
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
                                    const gsVector<index_t> & start,
                                    const gsVector<index_t> & end,
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
        low[dim] = start[dim];
        upp[dim] = end[dim];
        totalElements *= (end[dim] - start[dim]);
    }

    if (npieces <= 0 || totalElements == 0) 
    {
        return composite;
    }

    if (npieces == 1) 
    {
        composite->addDomain(
            memory::make_shared(new gsTensorSubDomain<d,T>(domain, start,end, patchId, parentPtr)));
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
        if (d > 1) 
        {
            splitAxis = (size[0] >= size[1]) ? 0 : 1;
            for(int i = 2; i < d; ++i) 
            {
                if (size[i] > size[splitAxis]) 
                {
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
        gsVector<index_t> leaf_low(d), leaf_upp(d);
        for (short_t dim = 0; dim < d; ++dim) 
        {
            leaf_low[dim] = leafNode->lowCorner()[dim];
            leaf_upp[dim] = leafNode->uppCorner()[dim];
        }
        auto subdomain = memory::make_shared(new gsTensorSubDomain<d,T>(domain, leaf_low, leaf_upp, patchId, parentPtr));
        composite->addDomain(subdomain);
    }

    delete root;
    
    return composite;
}

template<short_t d,class T>
typename gsDomain<T>::Ptr gsTensorSubDomain<d,T>::decompose(size_t npieces) const
{
    return decompose(m_parent, m_start,m_end, npieces, this->patchId(), m_parentPtr);
}

template<short_t d,class T>
std::vector<index_t> gsTensorSubDomain<d,T>::computeElementIndices() const
{
    std::vector<index_t> result;
    result.reserve(this->numElements());
    
    gsVector<index_t> cur = m_start;
    
    // Compute lexicographical order
    bool good = true;
    do
    {
        index_t globalIdx = tensorIndexToGlobal(cur);
        result.push_back(globalIdx);
        good = nextLexicographic(cur, m_start, m_end);
    } while (good);

    return result;
}


template<class T,short_t d>
gsTensorSubDomainIterator<T,d>::gsTensorSubDomainIterator(const gsTensorSubDomain<d,T>& subdomain)
: 
Base(subdomain), 
m_subdomain(subdomain), 
m_start(subdomain.start()), 
m_end(subdomain.end()), 
m_currentElement(0)
{
    // Initialize ranges
    m_tensorIdx.resize(d);
    for (short_t dim = 0; dim < d; ++dim)
    {
        m_tensorIdx[dim] = m_start[dim];
    }
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
    if (!good()) // needed??
        return;
    
    nextLexicographic(m_tensorIdx, m_start, m_end);
}

template<class T,short_t d>
void gsTensorSubDomainIterator<T,d>::next(index_t increment)
{
    m_currentElement += increment;
    
    if (!good())
        return;

    for (index_t i = 0; i < increment; ++i)
        nextLexicographic(m_tensorIdx, m_start, m_end);
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

