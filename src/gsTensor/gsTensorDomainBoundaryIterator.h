/** @file gsTensorDomainBoundaryIterator.h

    @brief Iterator over the boundary elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomainIterator.h>
#include <gsTensor/gsTensorDomain.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

/**
 * @brief Re-implements gsDomainIterator for iteration over all elements of the boundary of a <b>tensor product</b> parameter domain.\n
 * <em>See gsDomainIterator for more detailed documentation and an example of the typical use!!!</em>
 *
 * \ingroup Tensor
 */

// Class which enables iteration over all elements of a tensor product parameter domain
// Documentation in gsDomainIterator.h

template<class T, int D, typename uiter>
//template<class T, int D=Dynamic, typename uiter = typename gsKnotVector<T>::const_uiterator>
class gsTensorDomainBoundaryIterator : public gsDomainIterator<T>
{
public:

    gsTensorDomainBoundaryIterator(const gsTensorDomain<T,D> & domain, const boxSide & s)
    :
    gsDomainIterator<T>(domain, s),
    d( domain.dim() ),
    lower  ( gsVector<T, D>::Zero(d) ),
    upper  ( gsVector<T, D>::Zero(d) )
    {
        center  = gsVector<T, D>::Zero(d);
        par = s.parameter();
        dir = s.direction();
        meshStart.resize(d);
        meshEnd.resize(d);
        curElement.resize(d);

        for (int i=0; i < dir; ++i)
        {
            meshEnd[i]   = domain.component(i)->end();
            meshStart[i] = curElement[i] = domain.component(i)->begin();

            if (meshEnd[i] == curElement[i])
                m_isGood = false;
        }

        // Fixed direction
        meshEnd[dir]    = ( par ? domain.component(dir)->end() - 1 : domain.component(dir)->begin() + 1 );
        curElement[dir] =
        meshStart[dir]  = ( par ? domain.component(dir)->end() - 2 : domain.component(dir)->begin()     );
        tindex = curElement[dir] - domain.component(dir)->begin();

        for (int i=dir+1; i < d; ++i)
        {
            meshEnd[i]   = domain.component(i)->end() - 1;
            meshStart[i] = curElement[i] = domain.component(i)->begin();

            if (meshEnd[i] == curElement[i])
                m_isGood = false;
        }

        if (m_isGood)
            update();
    }


    gsTensorDomainBoundaryIterator( const gsBasis<T>& b, const boxSide & s )
    :
    gsTensorDomainBoundaryIterator(static_cast<const gsTensorDomain<T,D>&>(b.domain()), s)
    {
    }

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element; returns true if end not reached yet
    bool next()
    {
        m_isGood = m_isGood && nextLexicographic(curElement, meshStart, meshEnd);
        if (m_isGood)
            update();
        return m_isGood;
    }

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element (skipping #increment elements);
    // returns true if end not reached yet
    bool next(index_t increment)
    {
        for (index_t i = 0; i < increment; i++)
            m_isGood = m_isGood && nextLexicographic(curElement, meshStart, meshEnd);
        if (m_isGood)
            update();
        return m_isGood;
    }

    /// Resets the iterator implementation copied from the constructor
    /// note that it fails for sides containing 1 element in any direction
    /// do not know the rationale for it
    void reset()
    {
        curElement=meshStart;
        m_isGood = true;
        for(int i=0; i < d; ++i)
        {
            if (i!=dir && curElement[i]==meshEnd[i])
                m_isGood=false;
        }
        if (m_isGood)
            update();
    }

    /// Return the tensor index of the current element
    gsVector<unsigned, D> index() const
    {
        gsVector<unsigned, D> curr_index(d);
        for (int i = 0; i < dir; ++i)
            curr_index[i]  = curElement[i] - meshStart[i];
        for (int i = dir+1; i < d; ++i)
            curr_index[i]  = curElement[i] - meshStart[i];
        curr_index[dir]  = tindex;
        return curr_index;
    }

    const gsVector<T> & lowerCorner() const
    { return lower; }

    const gsVector<T> & upperCorner() const
    { return upper; }

    const T getPerpendicularCellSize() const
    {
        return *(curElement[dir]+1) - *curElement[dir];
    }

    GISMO_DEPRECATED
    void adjacent( const gsVector<bool> & orient,
                   gsDomainIterator<T>  & other )
    {
        GISMO_NO_IMPLEMENTATION
        // // 2D only for now

        // gsTensorDomainBoundaryIterator & other_ =
        //     static_cast< gsTensorDomainBoundaryIterator &>(other);

        // int a1 = !dir;
        // int a2 = !other_.dir;

        // other_.curElement[a2] = std::lower_bound(
        //     other_.breaks[a2].begin(), other_.breaks[a2].end(),
        //     orient[0] ? *curElement[a1] : *(curElement[a1]+1) );
        // other_.update();
    }

    /// Function to set the breakpoints in direction \a i manually
    void setBreaks(std::vector<T> newBreaks, index_t i) // i: direction
    {
        GISMO_NO_IMPLEMENTATION
        // breaks[i].swap(newBreaks);
        // meshEnd[i]   = breaks[i].end() - 1;
        // meshStart[i] = curElement[i] = breaks[i].begin();
        // reset();
    }

private:

    /// Computes lower, upper and center point of the current element, and computes the
    /// active functions. Plus some additional, boundary-iterator-specific things.
    void update()
    {
        for (int i = 0; i < dir ; ++i)
        {
            lower[i]  = curElement[i]->lowerCorner();
            upper[i]  = curElement[i]->upperCorner();
            center[i] = curElement[i]->center;
        }
        lower[dir]  =
        upper[dir]  =
        center[dir] = (par ? curElement[dir]->upperCorner() : curElement[dir]->lowerCorner() );
        for (int i = dir+1; i < d; ++i)
        {
            lower[i]  = curElement[i]->lowerCorner();
            upper[i]  = curElement[i]->upperCorner();
            center[i] = curElement[i]->center;
        }

        //gsDebug<<"lower: "<< lower.transpose() <<", upper="<<upper.transpose() <<"\n";
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

// Data members
protected:
    using gsDomainIterator<T>::m_isGood;
    using gsDomainIterator<T>::center;
    using gsDomainIterator<T>::m_side;

private:

    // the dimension of the parameter space
    short_t d;

    // Boundary parameters
    short_t  dir;
    bool par;
    unsigned tindex;


    // First mesh-line on the tensor grid
    gsVector<uiter, D> meshStart;

    // Last mesh-line on the tensor grid
    gsVector<uiter, D> meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<uiter, D> curElement;

    // parameter coordinates of current grid cell
    gsVector<T> lower, upper;

public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
}; // class gsTensorDomainBoundaryIterator


} // namespace gismo
