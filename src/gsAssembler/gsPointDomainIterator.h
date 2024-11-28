/** @file gsPointDomainIterator.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomainIterator.h>

namespace gismo
{

template<class T>
class gsPointDomainIterator : public gsDomainIterator<T>
{
private:
    typedef typename std::vector<T>::const_iterator  uiter;

public:

    gsPointDomainIterator(const gsPointDomain<T,D> & domain)
    :
    d( domain.dim() ),
    N( domain.points().cols() ),
    m_id(0),
    lower  ( domain.points().col(0) ),
    upper  ( domain.points().col(0) )
    {
        center  = domain.points().col(0);
        m_isGood = (domain.points().col(0)!=domain.points().col(N-1));
    }

    // Documentation in gsDomainIterator.h
    bool next()
    {
        m_isGood = m_isGood && (m_id < N-1);
        if (m_isGood)
        {
            ++m_id;
            lower = upper = center = domain.points().col(m_id);
        }
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    bool next(index_t increment)
    {
        m_isGood = m_isGood && (m_id + increment < N);
        if (m_isGood)
        {
            m_id += increment;
            lower = upper = center = domain.points().col(m_id);
        }
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    void reset()
    {
        m_id = 0;
        lower = upper = center = domain.points().col(m_id);
    }

    const gsVector<T> & lowerCorner() const
    { return lower; }

    const gsVector<T> & upperCorner() const
    { return upper; }

    index_t domainDim() const {return d;}

// Data members
public:
    using gsDomainIterator<T>::center;

protected:
    using gsDomainIterator<T>::m_id;
    using gsDomainIterator<T>::m_isGood;

private:
    // the dimension of the parameter space
    int d;

    // parameter coordinates of current grid cell
    gsVector<T> lower, upper;

}; // class gsPointDomainIterator


} // namespace gismo
