/** @file gsTensorDomainIterator.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTensorDomain.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

template<class T>
class gsKnotDomainIterator : public gsDomainIterator<T>
{
private:
    typedef typename gsKnotVector<T>::const_uiterator knotIterator;
    typedef typename gsDomainIterator<T>::uPtr domainIter;

    // Data members
    knotIterator m_it, m_itEnd;

public:

    gsKnotDomainIterator(const gsKnotVector<T> & _knots, bool start = true)
    :
    gsDomainIterator<T>(start ? 0 : _knots.numElements()),
    m_it(start ? _knots.domainUBegin() : _knots.domainUEnd()),
    m_itEnd(_knots.domainUEnd())
    {

    }

    gsKnotDomainIterator(const gsKnotDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsKnotDomainIterator(*this)); }

    // Documentation in gsDomainIterator.h
    void next() override
    {
        ++m_it;
    }

    // Documentation in gsDomainIterator.h
    void next(index_t increment) override
    {
        m_it += increment;
    }

    // Documentation in gsDomainIterator.h
    void prev() override
    {
        --m_it;
    }

    // Documentation in gsDomainIterator.h
    void prev(index_t decrement) override
    {
        m_it -= decrement;
    }

    // Documentation in gsDomainIterator.h
    void reset() override
    {
        m_it.reset();
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower;
        lower.resize(1);
        lower[0] = m_it.value();
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper;
        upper.resize(1);
        upper[0] = (m_it+1).value();
        return upper;
    }

    bool isBoundaryElement() const override
    {
        return ( 0==m_it.uIndex() || m_it+1==m_itEnd);
    }

    index_t domainDim() const {return 1;}

}; // class gsKnotDomainIterator


} // namespace gismo
