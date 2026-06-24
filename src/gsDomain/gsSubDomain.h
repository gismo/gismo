/** @file gsSubDomain.h

    @brief Abstract base for a sub-domain of a parent domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsDomain.h>

namespace gismo
{

/**
   @brief Abstract base class for a sub-domain of a parent domain.

   A gsSubDomain represents a subset of elements from a parent gsDomain,
   identified by a set of global element indices. Concrete subclasses
   provide the iteration mechanism over the selected elements.

   The primary concrete implementation is gsIndexSubDomain, which stores
   a sorted list of global element indices and iterates them by walking the
   parent domain's iterator.

   \ingroup Domain
*/
template<class T>
class gsSubDomain : public gsDomain<T>
{
public:

    /// Returns the parent domain from which this subdomain is drawn
    virtual const gsDomain<T>& parentDomain() const = 0;

    /// Returns the sorted list of global element indices belonging to this subdomain
    virtual const std::vector<index_t>& elementIndices() const = 0;

    /// Returns true if the given global element index belongs to this subdomain
    virtual bool contains(index_t elementIndex) const = 0;

}; // class gsSubDomain

} // namespace gismo
