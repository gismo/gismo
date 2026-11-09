/** @file gsKnotDomainIterator.h

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
    bool m_bdr;

public:

    gsKnotDomainIterator(const gsKnotVector<T> & _knots, bool start = true)
    :
    gsDomainIterator<T>(start ? 0 : _knots.numElements()),
    m_it(start ? _knots.domainUBegin() : _knots.domainUEnd()),
    m_itBegin(_knots.domainUBegin()),
    m_itEnd(_knots.domainUEnd()),
    m_bdr(false)
    {

    }

    /** @brief Boundary mode: iterates over the boundary side \a s of the
        1-D domain.

        The boundary of an interval is its end point, so the range holds
        exactly one degenerate element [x_s, x_s], with x_s the first
        (boundary::west) or last (boundary::east) point of the domain. The
        one-point Gauss rule (node 0, weight 2) of the fixed direction,
        mapped onto it with gsQuadRule::mapTo (which replaces a zero element
        size by 0.5), yields the point x_s with weight 1.

        The range is [begin,end) with ids 0 and 1 (\a start true or false).
        Any other side is not a side of a 1-D domain: the range is then
        empty (begin and end both have id 0).

        The position of the iterator is fixed; only the id advances.
    */
    gsKnotDomainIterator(const gsKnotVector<T> & _knots, const boxSide & s,
                         bool start = true)
    :
    gsDomainIterator<T>((start || !isSide(s)) ? 0 : 1, s),
    m_it(s.index()==boundary::east ? _knots.domainUEnd() : _knots.domainUBegin()),
    m_itBegin(_knots.domainUBegin()),
    m_itEnd(_knots.domainUEnd()),
    m_bdr(true)
    {
        this->setPatchIndex(_knots.patchIndex());
    }

    gsKnotDomainIterator(const gsKnotDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsKnotDomainIterator(*this)); }

    // Documentation in gsDomainIterator.h
    void next() override
    {
        if (m_bdr) return;
        ++m_it;
    }

    // Documentation in gsDomainIterator.h
    void next(index_t increment) override
    {
        if (m_bdr) return;
        m_it += increment;
    }

    // Documentation in gsDomainIterator.h
    void prev() override
    {
        if (m_bdr) return;
        --m_it;
    }

    // Documentation in gsDomainIterator.h
    void prev(index_t decrement) override
    {
        if (m_bdr) return;
        m_it -= decrement;
    }

    // Documentation in gsDomainIterator.h
    void reset() override
    {
        if (m_bdr) return;
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
        upper[0] = m_bdr ? m_it.value() : (m_it+1).value();
        return upper;
    }

    bool isBoundaryElement() const override
    {
        return m_bdr || ( 0==m_it.uIndex() || m_it+1==m_itEnd);
    }

    /// In boundary mode: the length of the knot span adjacent to the end point.
    const T getPerpendicularCellSize() const override
    {
        if (!m_bdr)
            return gsDomainIterator<T>::getPerpendicularCellSize();
        return this->side().parameter()
            ? m_itEnd.value() - (m_itEnd-1).value()
            : (m_itBegin+1).value() - m_itBegin.value();
    }

    index_t domainDim() const {return 1;}

private:
    static bool isSide(const boxSide & s)
    { return s.index()==boundary::west || s.index()==boundary::east; }

}; // class gsKnotDomainIterator


} // namespace gismo
