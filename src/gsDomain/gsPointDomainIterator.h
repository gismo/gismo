/** @file gsPointDomainIterator.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>

namespace gismo
{

template<class T>
class gsPointDomainIterator : public gsDomainIterator<T>
{
private:
    typedef typename std::vector<T>::const_iterator  uiter;

    typedef typename gsDomainIterator<T>::uPtr domainIter;
public:

    gsPointDomainIterator(const gsPointDomain<T> & domain)
    :
    gsDomainIterator<T>(), m_domain(domain)
    { }

    domainIter clone() const override { return domainIter(new gsPointDomainIterator(*this)); }

    /// Documentation in gsDomainIterator.h
    void next() override
    {
        /* Transparent, since only the ID needs to be updated (delegated to operator++ in gsDomainIterator) */
    }

    /// Documentation in gsDomainIterator.h
    void next(index_t increment) override
    {
        /* Transparent, since only the ID needs to be updated (delegated to operator+= in gsDomainIterator) */
    }

    /// Documentation in gsDomainIterator.h
    void prev() override
    {
        /* Transparent, since only the ID needs to be updated (delegated to operator-- in gsDomainIterator) */
    }

    /// Documentation in gsDomainIterator.h
    void prev(index_t decrement) override
    {
        /* Transparent, since only the ID needs to be updated (delegated to operator-= in gsDomainIterator) */
    }

    gsVector<T> lowerCorner() const override
    { return m_domain.points().col(m_id); }

    gsVector<T> upperCorner() const override
    { return m_domain.points().col(m_id); }

// Data members
protected:
    using gsDomainIterator<T>::m_id;
    const gsPointDomain<T> & m_domain;

}; // class gsPointDomainIterator


} // namespace gismo
