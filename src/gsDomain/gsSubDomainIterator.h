/** @file gsSubDomainIterator.h

    @brief Abstract base class for subdomain iterators.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>

namespace gismo
{

/**
    @brief Abstract base class for iterators over subdomains.
    
    Inherits from gsDomainIterator and provides the base interface
    for specialized subdomain iterators.
    
    \ingroup Domain
*/
template<class T>
class gsSubDomainIterator : public gsDomainIterator<T>
{
public:
    typedef gsDomainIterator<T> Base;
    
    virtual ~gsSubDomainIterator() { }

    virtual index_t patch() const override { return m_subdomain.patchId(); }

protected:
    /// Protected constructor
    explicit gsSubDomainIterator(const gsSubDomain<T>& subdomain, index_t _id = 0) 
        : Base(_id, subdomain.patchId()), m_subdomain(subdomain) { }

private:
    const gsSubDomain<T>& m_subdomain;
};

} // namespace gismo
