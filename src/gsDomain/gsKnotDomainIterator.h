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
        // Threaded loops step past the end when there are fewer elements
        // than threads; stop at the end, as the tensor iterators do.
        m_it += std::min<index_t>(increment, m_itEnd - m_it);
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

/**
 * @brief Iterator over the boundary of a knot-vector domain on side
 * west or east: one element, the end point of the domain on that side.
 *
 * The element has zero length (lowerCorner() == upperCorner()), like the
 * fixed direction of gsTensorDomainBoundaryIterator in higher dimensions;
 * gsQuadRule::mapTo maps a rule with one node in that direction onto the
 * point with weight 1.
 */
template<class T>
class gsKnotDomainBoundaryIterator : public gsDomainIterator<T>
{
private:
    typedef typename gsDomainIterator<T>::uPtr domainIter;

    // Data members
    T m_point; // end knot of the domain on the side
    T m_width; // length of the knot span next to it

public:

    gsKnotDomainBoundaryIterator(const gsKnotVector<T> & _knots, const boxSide & s)
    : gsDomainIterator<T>(0, s)
    {
        GISMO_ENSURE(s == boundary::west || s == boundary::east,
                     "A knot vector domain has the boundary sides west and east, got "<< s);
        if (s.parameter())
        {
            m_point = *_knots.domainUEnd();
            m_width = m_point - *(_knots.domainUEnd() - 1);
        }
        else
        {
            m_point = *_knots.domainUBegin();
            m_width = *(_knots.domainUBegin() + 1) - m_point;
        }
    }

    gsKnotDomainBoundaryIterator(const gsKnotDomainBoundaryIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsKnotDomainBoundaryIterator(*this)); }

    // A single element: the loops compare element ids with endBdr(), so
    // stepping changes the id only
    void next() override { }
    void next(index_t) override { }
    void prev() override { }
    void prev(index_t) override { }
    void reset() override { }

    gsVector<T> lowerCorner() const override
    { return gsVector<T>::Constant(1, m_point); }

    gsVector<T> upperCorner() const override
    { return gsVector<T>::Constant(1, m_point); }

    const T getPerpendicularCellSize() const override { return m_width; }

    bool isBoundaryElement() const override { return true; }

}; // class gsKnotDomainBoundaryIterator


} // namespace gismo
