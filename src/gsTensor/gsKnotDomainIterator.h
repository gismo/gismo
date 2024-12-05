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

    gsKnotDomainIterator(const gsKnotVector<T> & domain)
    : lower(1),
      upper(1)
    {
        // compute breaks and mesh size
        meshEnd    = domain->domainEnd();
        curElement = meshStart = domain->domainBegin();

        if (meshEnd == meshStart)
            m_isGood = false;

        if (m_isGood)
            update();
    }

    // Documentation in gsDomainIterator.h
    bool next()
    {
        ++curElement;
        m_isGood = (curElement != meshEnd);

        if (m_isGood)
        {
            update();
            ++m_id; //increment id
        }
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    bool next(index_t increment)
    {
        curElement += increment;
        m_isGood = (curElement < meshEnd);

        if (m_isGood)
        {
            update();
            m_id += increment; //increment id
        }
        return m_isGood;
    }

    // Documentation in gsDomainIterator.h
    void reset()
    {
        m_id = 0;
        curElement = meshStart;

        m_isGood = (meshEnd != meshStart);
        if (m_isGood)
            update();
    }

    /// return the tensor index of the current element
    gsVector<unsigned, 1> index() const
    {
        gsVector<unsigned, 1> curr_index(1);
        curr_index[0]  = curElement.uIndex();
        return curr_index;
    }

    void getVertices(gsMatrix<T>& result)
    {
        result.resize( 1, 2);
        result(0,0) = curElement.value();
        result(0,1) = (curElement+1).value();
    }

    const gsVector<T> & lowerCorner() const
    { return lower; }

    const gsVector<T> & upperCorner() const
    { return upper; }

    bool isBoundaryElement() const
    {
        return (curElement==meshStart || curElement+1==meshEnd);
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
    using gsDomainIterator<T>::m_id;
    using gsDomainIterator<T>::m_isGood;

private:
    // Extent of the tensor grid
    gsVector<domainIter, 1> meshStart, meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<domainIter, 1> curElement;

    // parameter coordinates of current grid cell
    gsVector<T> lower, upper;

public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
}; // class gsTensorDomainIterator


} // namespace gismo
