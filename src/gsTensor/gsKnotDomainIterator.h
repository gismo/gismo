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
class gsKnotDomainIterator : public gsTensorDomainIterator<T,1>
{
private:
    typedef typename gsKnotVector<T>::const_uiterator domainIter;

public:

    gsKnotDomainIterator(const gsKnotVector<T> & _knots)
    : m_knotvector(&_knots), m_it(_knots.domainUBegin())
    {

    }

    // Documentation in gsDomainIterator.h
    bool next()
    {
        ++m_it;
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    bool next(index_t increment)
    {
        m_it += increment;
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    void reset()
    {
        m_it.reset();
    }

    const gsVector<T> & lowerCorner() const
    { return lower; }

    const gsVector<T> & upperCorner() const
    { return upper; }

    bool isBoundaryElement() const
    {
        return ( 0==m_it.uIndex() || curElement+1==meshEnd);
    }

    index_t domainDim() const {return 1;}

private:

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void update()
    {
        lower[0]  = curElement.value();
        upper[0]  = (curElement+1).value();
        center  = (lower + upper) / (T)(2);
    }


//    size_t numElements() const
//    {
//
//    }

// Data members
public:
    using gsDomainIterator<T>::center;

protected:
    using gsDomainIterator<T>::m_isGood;

private:
    gsKnotVector<T> * m_knotvector;

    domainIter m_it;

}; // class gsKnotDomainIterator


} // namespace gismo
