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
class gsKnotDomainIterator : public gsDomainIterator<T>
{
private:
    typedef typename gsKnotVector<T>::const_uiterator domainIter;

public:

    gsKnotDomainIterator(const gsKnotVector<T> & _knots, bool start = true)
    :
    m_it(start ? _knots.domainUBegin() : _knots.domainUEnd()-1),
    m_itEnd(_knots.domainUEnd()-1)
    {

    }

    // Documentation in gsDomainIterator.h
    bool next() override
    {
        ++m_it;
        m_isGood = (m_it != m_itEnd);
        gsDebugVar(m_it.value());
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    bool next(index_t increment) override
    {
        m_it += increment;
        m_isGood = (m_it < m_itEnd);
        gsDebugVar(m_it.value());
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    bool prev() override
    {
        --m_it;
        // WARNING: This does not check whether the iterator passes the beginning
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    bool prev(index_t decrement) override
    {
        m_it -= decrement;
        // WARNING: This does not check whether the iterator passes the beginning
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    void reset()
    {
        m_it.reset();
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T,1> lower;
        lower[0] = m_it.value();
        gsDebugVar(lower[0]);
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T,1> upper;
        upper[0] = (m_it+1).value();
        return upper;
    }

    bool isBoundaryElement() const
    {
        return ( 0==m_it.uIndex() || m_it+1==m_itEnd);
    }

    index_t domainDim() const {return 1;}

// Data members
protected:
    using gsDomainIterator<T>::m_isGood;

private:
    gsVector<T,1> lower, upper;

    domainIter m_it, m_itEnd;

}; // class gsKnotDomainIterator


} // namespace gismo
