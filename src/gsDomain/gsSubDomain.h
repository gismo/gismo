/** @file gsSubDomain.h

    @brief Abstract base class for subdomains of a domain.

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
    @brief Abstract base class representing a subdomain of a gsDomain.
    
    A subdomain is a subset of elements from a parent domain. It provides
    iteration and querying functionality for its elements.
    
    \ingroup Domain
*/
template<class T>
class gsSubDomain : public gsDomain<T>
{
public:
    typedef typename gsDomain<T>::Ptr domainPtr;
    typedef typename gsDomain<T>::uPtr domainUPtr;
    typedef memory::shared_ptr<gsSubDomain<T>> Ptr;
    typedef memory::unique_ptr<gsSubDomain<T>> uPtr;

    explicit gsSubDomain(index_t patchId = -1) 
    : 
    gsDomain<T>(patchId) 
    {}

    explicit gsSubDomain(const gsDomain<T>& parent, index_t patchId = -1)
    : 
    gsDomain<T>(patchId) 
    {}

    virtual ~gsSubDomain() { }

    /// Return the parent domain
    virtual const gsDomain<T>& parentDomain() const = 0;

    /// Return number of elements in this subdomain
    size_t numElements() const override = 0;

    /// Return begin iterator
    virtual typename gsDomain<T>::iterator beginAll() const = 0;
    
    /// Return end iterator
    virtual typename gsDomain<T>::iterator endAll() const = 0;

    /// Get the element indices of this subdomain in the parent domain
    virtual const std::vector<index_t>& elementIndices() const = 0;

    /// Check if an element (by global index) belongs to this subdomain
    virtual bool contains(index_t elementIndex) const = 0;
};

} // namespace gismo
