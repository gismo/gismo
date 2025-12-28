/** @file gsDomain.hpp

    @brief Provides implementation of gsDomain methods.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsDomain.h>
#include <gsDomain/gsCompositeDomain.h> // Keep this for gsTensorSubDomain references
#include <gsDomain/gsIndexSubDomain.h> // Needed for gsCompositeDomain

namespace gismo
{

template<class T>
typename gsDomain<T>::Ptr
gsDomain<T>::decompose(size_t npieces) const
{
    // Base implementation: divide the total number of elements equally (assuming 1 piece)
    memory::shared_ptr<gsCompositeDomain<T>> result = memory::make_shared(new gsCompositeDomain<T>());
    index_t totalElements = this->numElements();
    index_t elementsPerPiece = totalElements / npieces;
    index_t remainingElements = totalElements % npieces;
    index_t pieces = static_cast<index_t>(npieces);
    for (index_t i = 0; i < pieces; ++i)
    {
        std::vector<index_t> elementIndices;
        elementIndices.reserve(elementsPerPiece + (i < remainingElements ? 1 : 0));
        index_t startIdx = i * elementsPerPiece + std::min(i, remainingElements);
        index_t endIdx = startIdx + elementsPerPiece + (i < remainingElements ? 1 : 0);
        for (index_t j = startIdx; j < endIdx; ++j)
        {
            elementIndices.push_back(j);
        }
        memory::shared_ptr<gsIndexSubDomain<T>> subdomain = memory::make_shared(new gsIndexSubDomain<T>(*this, elementIndices));
        result->addDomain(subdomain);
    }
    return result;
}


} // namespace gismo
