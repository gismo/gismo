/** @file gsHSubDomain.hpp

    @brief Implementation of gsHSubDomain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsHSubDomain.h>

namespace gismo
{

template<short_t d, class T, class Z>
typename gsHSubDomain<d, T, Z>::domainIter
gsHSubDomain<d, T, Z>::beginAll() const
{
    return domainIter(new gsHSubDomainIterator<d, T, Z>(*this));
}

template<short_t d, class T, class Z>
void gsHSubDomain<d, T, Z>::computeElementIndices()
{
    m_elementIndices.clear();
    
    // Get all elements from parent domain
    size_t numElements = m_parent.numElements();
    
    // For now, collect all elements - this is a simplified version
    // In a more sophisticated implementation, we would access the tree structure
    // to determine which elements belong to our level(s)
    for (size_t i = 0; i < numElements; ++i)
        m_elementIndices.push_back(i);
    
    // Sort indices for binary search in contains()
    std::sort(m_elementIndices.begin(), m_elementIndices.end());
}

} // namespace gismo
