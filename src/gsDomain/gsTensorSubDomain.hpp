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

namespace gismo
{

template<class T, int D>
typename gsTensorSubDomain<T, D>::domainIter
gsTensorSubDomain<T, D>::beginAll() const
{
    return domainIter(new gsTensorSubDomainIterator<T, D>(*this));
}

template<class T, int D>
std::vector<index_t>
gsTensorSubDomain<T, D>::computeElementIndices() const
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

} // namespace gismo
