/** @file gsTensorDomainIterator.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomainIterator.h>
#include <gsTensor/gsTensorDomain.h>
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

    gsBreaksIterator(const std::vector<T> & _knots, bool start = true)
    :
    gsDomainIterator<T>(start ? 0 : _knots.size()),
    m_it(start ? _knots.begin() : _knots.end()-1),
    m_itBegin(_knots.begin()),
    m_itEnd(_knots.end()-1)
    {
    }

    // Documentation in gsDomainIterator.h
    bool next() override
    {
        ++m_it;
        return (m_it != m_itEnd);
    }

    // Documentation in gsDomainIterator.h
    bool next(index_t increment) override
    {
        m_it += increment;
        return (m_it < m_itEnd);
    }

    // Documentation in gsDomainIterator.h
    bool prev() override
    {
        --m_it;
        // WARNING: This does not check whether the iterator passes the beginning
        return (m_it >= m_itBegin);
    }

    // Documentation in gsDomainIterator.h
    bool prev(index_t decrement) override
    {
        m_it -= decrement;
        // WARNING: This does not check whether the iterator passes the beginning
        return (m_it >= m_itBegin);
    }

    // Documentation in gsDomainIterator.h
    void reset()
    {
        m_it = m_itBegin;
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T,1> lower;
        lower[0] = *m_it;
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T,1> upper;
        upper[0] = *(m_it+1);
        return upper;
    }

    bool isBoundaryElement() const
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
    gsVector<T,1> lower, upper;

    domainIter m_it, m_itBegin, m_itEnd;

}; // class gsBreaksIterator


} // namespace gismo
