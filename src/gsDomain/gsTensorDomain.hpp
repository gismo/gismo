/** @file gsTensorDomain.hpp

    @brief Provides implementation of gsTensorDomain methods.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsTensorSubDomain.h> // Keep this for gsTensorSubDomain references
#include <gsDomain/gsCompositeDomain.h> // Needed for gsCompositeDomain

namespace gismo
{

template<class T, int D>
typename gsDomain<T>::Ptr
gsTensorDomain<T, D>::decompose(index_t npieces, decompositionStrategy strategy) const
{
    if (strategy != decompositionStrategy::tensor)
    {
        GISMO_ERROR("This decomposition strategy is not implemented for gsTensorDomain!");
    }

    auto composite = memory::make_shared(new gsCompositeDomain<T>());
    
    index_t totalElements = this->numElements();

    // Edge cases for npieces
    if (npieces <= 0 || totalElements == 0) {
        return composite; // Return empty composite domain
    }

    if (npieces == 1) {
        std::vector<typename gsTensorSubDomain<T, D>::Range> ranges;
        for (short_t dim = 0; dim < D; ++dim) {
            ranges.push_back(typename gsTensorSubDomain<T, D>::Range(0, m_knotVectors[dim]->numElements()));
        }
        auto subdomain = memory::make_shared(new gsTensorSubDomain<T, D>(*this, ranges, this->patchId()));
        composite->addDomain(subdomain);
        return composite;
    }
    
    // Determine the effective number of pieces to create, capped by total elements
    index_t effective_npieces = std::min(npieces, totalElements);
    
    // Get element counts per dimension
    std::vector<index_t> elementsPerDim(D);
    for (short_t dim = 0; dim < D; ++dim)
        elementsPerDim[dim] = m_knotVectors[dim]->numElements();
    
    std::vector<index_t> numSplits(D, 1);
    index_t current_product_of_splits = 1;

    // Greedily distribute splits to reach effective_npieces
    // Repeatedly split the dimension that currently has the largest "chunk" size (elements / splits)
    // until we have effective_npieces subdomains or cannot split further.
    while (current_product_of_splits < effective_npieces) {
        short_t bestDim = -1;
        real_t maxChunkSize = -1.0; 

        for (short_t dim = 0; dim < D; ++dim) {
            // Check if this dimension can still be split further
            if (numSplits[dim] < elementsPerDim[dim]) { 
                // Calculate chunk size if we were to split this dimension one more time
                real_t potentialChunkSize = static_cast<real_t>(elementsPerDim[dim]) / (numSplits[dim] + 1);
                if (potentialChunkSize > maxChunkSize) {
                    maxChunkSize = potentialChunkSize;
                    bestDim = dim;
                }
            }
        }

        if (bestDim != -1) {
            numSplits[bestDim]++;
            current_product_of_splits = 1; // Recalculate product after split
            for (short_t dim = 0; dim < D; ++dim) {
                current_product_of_splits *= numSplits[dim];
            }
        } else {
            // Cannot split any more dimensions, we've reached the maximum possible tensor product subdomains
            break;
        }
    }
    
    // Create ranges for each dimension based on numSplits
    typedef typename gsTensorSubDomain<T, D>::Range Range;
    std::vector<std::vector<Range>> dimRanges(D);
    for (short_t dim = 0; dim < D; ++dim) {
        index_t elemPerDim = elementsPerDim[dim];
        index_t splits = numSplits[dim]; // Use the determined number of splits
        
        for (index_t s = 0; s < splits; ++s) {
            index_t start = (s * elemPerDim) / splits;
            index_t end = ((s + 1) * elemPerDim) / splits;
            dimRanges[dim].push_back(Range(start, end));
        }
    }
    
    // Create all tensor subdomain combinations
    std::function<void(short_t, std::vector<Range>&)> 
        createSubdomains = [&](short_t dim, std::vector<Range>& current) {
        if (dim == D) {
            auto subdomain = memory::make_shared(
                new gsTensorSubDomain<T, D>(*this, current, this->patchId())); 
            composite->addDomain(subdomain);
            return;
        }
        
        current.resize(D); // Ensure current has enough space
        for (const auto& range : dimRanges[dim]) {
            current[dim] = range;
            createSubdomains(dim + 1, current);
        }
    };
    
    std::vector<Range> current(D);
    createSubdomains(0, current);
    
    return composite;
}

template<class T, int D>
typename gsDomain<T>::Ptr
gsTensorDomain<T, D>::decompose(index_t npieces) const
{
    return this->decompose(npieces, decompositionStrategy::tensor);
}


template<class T, int D>
size_t gsTensorDomain<T, D>::numElements() const
{
    size_t nElem = 1;
    for (short_t dim = 0; dim < D; ++dim)
        nElem *= m_knotVectors[dim]->numElements();
    return nElem;
}

} // namespace gismo
