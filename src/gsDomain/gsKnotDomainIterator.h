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
    knotIterator m_it, m_itBegin, m_itEnd;

public:

    gsKnotDomainIterator(const gsKnotVector<T> & _knots, bool start = true)
    :
    gsDomainIterator<T>(start ? 0 : _knots.numElements()),
    m_it(start ? _knots.domainUBegin() : _knots.domainUEnd()),
    m_itBegin(_knots.domainUBegin()),
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
        m_it = m_itBegin;
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


/** @brief Iterator over the boundary elements of a one-dimensional
 *  (knot vector) domain.
 *
 *  The boundary of a 1D domain is a set of points, so every element visited
 *  here is degenerate: its lower and upper corner coincide with a domain
 *  endpoint. This mirrors gsTensorDomainBoundaryIterator, which likewise
 *  collapses the element extent along the fixed direction.
 *
 *  getPerpendicularCellSize() returns the length of the adjacent knot span --
 *  the h entering Nitsche-type boundary penalties -- again matching the
 *  higher-dimensional iterator.
 */
template<class T>
class gsKnotDomainBoundaryIterator : public gsDomainIterator<T>
{
private:
    typedef typename gsKnotVector<T>::const_uiterator knotIterator;
    typedef typename gsDomainIterator<T>::uPtr domainIter;

    // Data members: the endpoint(s) of the domain and the length of the knot
    // span adjacent to each. Two entries for boundary::all, one otherwise.
    std::vector<T> m_values, m_sizes;
    index_t m_pos;

public:

    gsKnotDomainBoundaryIterator(const gsKnotVector<T> & _knots,
                                 const boxSide & s, bool start = true)
    :
    gsDomainIterator<T>(0, s),
    m_pos(0)
    {
        const bool west = (boundary::all == s || !s.parameter());
        const bool east = (boundary::all == s ||  s.parameter());

        const knotIterator ub = _knots.domainUBegin();
        const knotIterator ue = _knots.domainUEnd();

        if (west)
        {
            m_values.push_back( ub.value() );
            m_sizes .push_back( (ub+1).value() - ub.value() );
        }
        if (east)
        {
            m_values.push_back( ue.value() );
            m_sizes .push_back( ue.value() - (ue-1).value() );
        }

        if (!start)
        {
            m_pos = static_cast<index_t>(m_values.size());
            this->m_id = m_pos;
        }
    }

    gsKnotDomainBoundaryIterator(const gsKnotDomainBoundaryIterator & other) = default;
    domainIter clone() const override
    { return domainIter(new gsKnotDomainBoundaryIterator(*this)); }

    // Documentation in gsDomainIterator.h
    void next() override { ++m_pos; }

    // Documentation in gsDomainIterator.h
    void next(index_t increment) override { m_pos += increment; }

    // Documentation in gsDomainIterator.h
    void prev() override { --m_pos; }

    // Documentation in gsDomainIterator.h
    void prev(index_t decrement) override { m_pos -= decrement; }

    // Documentation in gsDomainIterator.h
    void reset() override { m_pos = 0; }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower(1);
        lower[0] = m_values[m_pos];
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper(1);
        upper[0] = m_values[m_pos];
        return upper;
    }

    const T getPerpendicularCellSize() const override
    { return m_sizes[m_pos]; }

    bool isBoundaryElement() const override { return true; }

    index_t domainDim() const {return 1;}

}; // class gsKnotDomainBoundaryIterator


} // namespace gismo
