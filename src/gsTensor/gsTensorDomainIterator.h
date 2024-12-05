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
// Documentation in gsDomainIterator.h
// Class which enables iteration over all elements of a tensor product parameter domain

/**
 * @brief Re-implements gsDomainIterator for iteration over all elements of a <b>tensor product</b> parameter domain.\n
 * <em>See gsDomainIterator for more detailed documentation and an example of the typical use!!!</em>
 *
 * \ingroup Tensor
 */

template<class T, int D>
class gsTensorDomainIterator : public gsDomainIterator<T>
{
private:
    typedef typename gsDomainIterator<T>::uPtr domainIter;

public:

    gsTensorDomainIterator(const gsTensorDomain<T,D> & domain)
    : lower  ( gsVector<T, D>::Zero(D) ),
      upper  ( gsVector<T, D>::Zero(D) )
    {
        center  = gsVector<T, D>::Zero(D);

        // compute breaks and mesh size
        meshStart.resize(D);
        meshEnd.resize(D);
        curElement.resize(D);

        for (int i=0; i < D; ++i)
        {
            meshEnd[i]    = domain.component(i)->end();
            curElement[i] = meshStart[i] = domain.component(i)->begin();

            if (meshEnd[i] == meshStart[i])
                m_isGood = false;
        }

        if (m_isGood)
            update();
    }

    // gsTensorDomainIterator(const gsTensorDomain<T,D> & domain, boxSide s=boundary::none)
    // : lower  ( gsVector<T, D>::Zero(D) ),
    //   upper  ( gsVector<T, D>::Zero(D) )
    // {
    //     center  = gsVector<T, D>::Zero(D);

    //     // compute breaks and mesh size
    //     meshStart.resize(D);
    //     meshEnd.resize(D);
    //     curElement.resize(D);

    //     for (int i=0; i < D; ++i)
    //     {
    //         meshEnd[i]    = domain.component(i)->end();
    //         curElement[i] = meshStart[i] = domain.component(i)->begin();

    //         if (meshEnd[i] == meshStart[i])
    //             m_isGood = false;
    //     }

    //     if (s!=boundary::none)
    //     {
    //         auto par = s.parameter();
    //         auto dir = s.direction();
    //         // Fixed direction
    //         meshEnd[dir]    = ( par ? domain.component(dir)->end() - 1 : domain.component(dir)->begin() + 1 );
    //         curElement[dir] =
    //         meshStart[dir]  = ( par ? domain.component(dir)->end() - 2 : domain.component(dir)->begin()     );
    //         // tindex = curElement[dir] - domain.component(dir)->begin();
    //     }


    //     if (m_isGood)
    //         update();
    // }

    // Documentation in gsDomainIterator.h
    bool next()
    {
        this->advance();
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
        for (index_t i = 0; i < increment; i++)
            this->advance();
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
        for (index_t i = 0; i < D; ++i)
            curElement[i]->reset();

        m_isGood = ( meshEnd.array() != curElement.array() ).all() ;
        if (m_isGood)
            update();
    }

    /// return the tensor index of the current element
    gsVector<unsigned, D> index() const
    {
        gsVector<unsigned, D> curr_index(D);
        for (int i = 0; i < D; ++i)
            curr_index[i]  = curElement[i]->index();
        return curr_index;
    }

    void getVertices(gsMatrix<T>& result)
    {
        result.resize( D, 1 << D);

        gsVector<T,D> v, l, u;
        l.setZero();
        u.setOnes();
        v.setZero();
        int r = 0;
        do {
            for ( int i = 0; i< D; ++i)
                result(i,r) = ( v[i] ? upper[i] : lower[i] );
        }
        while ( nextCubeVertex(v, l, u) );
    }

    const gsVector<T> & lowerCorner() const
    { return lower; }

    const gsVector<T> & upperCorner() const
    { return upper; }

    bool isBoundaryElement() const
    {
        for (int i = 0; i< D; ++i)
            if ( curElement[i]->isBoundaryElement() )
                return true;
        return false;
    }

    index_t domainDim() const {return D;}

private:

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void update()
    {
        for (int i = 0; i < D; ++i)
        {
            lower[i]  = curElement[i]->lowerCorner();
            upper[i]  = curElement[i]->upperCorner();
            center[i] = curElement[i]->center;
        }
    }

    void advance()
    {
        for (index_t i = 0; i < D; ++i)
        {
            // increase current dimension
            if (++(*curElement[i]) == meshEnd[i])     // current dimension exhausted ?
            {
                if (i == D - 1)         // was it the last one?
                    m_isGood = false;       // then all elements exhausted
                else
                    curElement[i]->reset();  // otherwise, reset this and increase the next dimension
            }
            else
                m_isGood = true;            // current dimension not yet exhausted, return current vector
        }
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
    gsVector<domainIter, D> meshStart, meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<domainIter, D> curElement;

    // parameter coordinates of current grid cell
    gsVector<T> lower, upper;

public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
}; // class gsTensorDomainIterator

} // namespace gismo
