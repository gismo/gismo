/** @file gsTensorDomainBoundaryIterator.h

    @brief Iterator over the boundary elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTensorDomain.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsDomain/gsBreaksIterator.h>

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
class gsTensorDomainBoundaryIterator : public gsDomainIterator<T>
{
    typedef gsDomainIteratorWrapper<T> domainIterWrapper;
    typedef typename gsDomain<T>::uPtr domainPtr;
    using typename  gsDomainIterator<T>::uPtr;
    using typename  gsDomainIterator<T>::Ptr;
public:

    explicit gsTensorDomainBoundaryIterator(const gsTensorDomain<T,D> & domain,
                                            const boxSide & s)
    : gsDomainIterator<T>(0, s),
      d( domain.dim() )
    {
        par = s.parameter();
        dir = s.direction();
        meshStart.resize(d);
        meshEnd.resize(d);
        curElement.resize(d);

        for (int i=0; i < dir; ++i)
        {
            meshEnd[i]    = give(domain.component(i)->endAll());
            meshStart[i]  = give(domain.component(i)->beginAll());
            curElement[i] = give(domain.component(i)->beginAll());
        }

        // Fixed direction
        if (par)
        {
            meshEnd[dir]    = give(domain.component(dir)->endAll());
            curElement[dir] = give(domain.component(dir)->endAll());
            curElement[dir]-=1;
            meshStart[dir]  = give(domain.component(dir)->endAll());
            meshStart[dir] -=1; //note: ending value
        }
        else
        {
            meshEnd[dir]    = give(domain.component(dir)->beginAll());
            meshEnd[dir]   +=1;
            curElement[dir] = give(domain.component(dir)->beginAll());
            meshStart[dir]  = give(domain.component(dir)->beginAll());
        }

        tindex = curElement[dir] - domain.component(dir)->beginAll();

        for (int i=dir+1; i < d; ++i)
        {
            meshEnd[i]    = give(domain.component(i)->endAll());
            meshStart[i]  = give(domain.component(i)->beginAll());
            curElement[i] = give(domain.component(i)->beginAll());
        }
    }


    gsTensorDomainBoundaryIterator( const gsBasis<T>& b, const boxSide & s )
    :
    gsTensorDomainBoundaryIterator(static_cast<const gsTensorDomain<T,D>&>(*b.domain()), s)
    { }

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element; returns true if end not reached yet
    void next() override
    {
        nextLexicographicIter(curElement, meshEnd, dir);
    }

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element (skipping #increment elements);
    // returns true if end not reached yet
    void next(index_t increment) override
    {
        bool isGood(true);
        for (index_t i = 0; i < increment; i++)
            isGood = isGood && nextLexicographicIter(curElement, meshEnd, dir);
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

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower;
        lower.resize(d);
        for (short_t i = 0; i < dir ; ++i)
            lower[i]  = curElement[i].lowerCorner().value();
        lower[dir]  = (par ? curElement[dir].upperCorner().value() : curElement[dir].lowerCorner().value() );
        for (short_t i = dir+1; i < d; ++i)
            lower[i]  = curElement[i].lowerCorner().value();
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper;
        upper.resize(d);
        for (short_t i = 0; i < dir ; ++i)
            upper[i]  = curElement[i].upperCorner().value();
        upper[dir]  = (par ? curElement[dir].upperCorner().value() : curElement[dir].lowerCorner().value() );
        for (short_t i = dir+1; i < d; ++i)
            upper[i]  = curElement[i].upperCorner().value();
        return upper;
    }

    const T getPerpendicularCellSize() const override
    {
        return curElement[dir].upperCorner().value() - curElement[dir].lowerCorner().value();
    }

    GISMO_DEPRECATED
    void adjacent( const gsVector<bool> & /* orient */,
                   gsDomainIterator<T>  & /* other */ ) override
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
    void setBreaks(const std::vector<T> & newBreaks, index_t i) // i: direction
    {
        GISMO_ASSERT(i!=dir, "Changing non-boundary breakpoints is not supported.");
        meshStart[i] = give(gsBreaksIterator<T>::make(newBreaks, true));
        curElement[i] = give(gsBreaksIterator<T>::make(newBreaks, true));
        meshEnd[i]   = give(gsBreaksIterator<T>::make(newBreaks, false));

        // Note: reset() has a bug, therefore we do not call it at all
        //reset();
    }


// Data members
private:

    // the dimension of the parameter space
    short_t d;

    // Boundary parameters
    short_t  dir;
    bool par;
    unsigned tindex;


    // First mesh-line on the tensor grid
    gsVector<domainIterWrapper, D> meshStart;

    // Last mesh-line on the tensor grid
    gsVector<domainIterWrapper, D> meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<domainIterWrapper, D> curElement;

public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
}; // class gsTensorDomainBoundaryIterator


} // namespace gismo
