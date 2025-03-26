/** @file gsBreaksIterator.h

    @brief Iterator over the elements of 1D elements defined by breaks in a std container

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTensorDomain.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

template<class T>
class gsBreaksIterator : public gsDomainIterator<T>
{
private:
    // todo: generalize for any kind of standard iterator
    typedef typename std::vector<T>::const_iterator domainIter;

public:

    explicit gsBreaksIterator(const std::vector<T> & _knots, bool start = true)
    :
    gsDomainIterator<T>(start ? 0 : _knots.size()-1),
    m_it(start ? _knots.begin() : _knots.end()-1),
    m_itBegin(_knots.begin()),
    m_itEnd(_knots.end()-1)
    {
    }

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
        m_it = m_itBegin;
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower;
        lower.resize(1);
        lower[0] = *m_it;
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper;
        upper.resize(1);
        upper[0] = *(m_it+1);
        return upper;
    }

    bool isBoundaryElement() const override
    {
        return ( m_itBegin==m_it || m_it+1==m_itEnd);
    }

    index_t domainDim() const {return 1;}

    static gsDomainIteratorWrapper<T> make(const std::vector<T> & knots, bool start = true)
    {
        return gsDomainIteratorWrapper<T>(new gsBreaksIterator<T>(knots, start));
    }

// Data members
private:
    domainIter m_it, m_itBegin, m_itEnd;

}; // class gsBreaksIterator


} // namespace gismo
