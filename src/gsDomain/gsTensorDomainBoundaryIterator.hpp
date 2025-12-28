/** @file gsTensorDomainBoundaryIterator.hpp

    @brief Iterator over the boundary elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsTensorDomainBoundaryIterator.h>

namespace gismo
{

template<class T, int D, typename uiter>
gsTensorDomainBoundaryIterator<T,D,uiter>::gsTensorDomainBoundaryIterator(const gsTensorDomain<D,T> & domain,
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

template<class T, int D, typename uiter>
gsTensorDomainBoundaryIterator<T,D,uiter>::gsTensorDomainBoundaryIterator( const gsBasis<T>& b, const boxSide & s )
    : gsTensorDomainBoundaryIterator(static_cast<const gsTensorDomain<D,T>&>(*b.domain()), s)
{ }

template<class T, int D, typename uiter>
void gsTensorDomainBoundaryIterator<T,D,uiter>::next()
{
    nextLexicographicIter(curElement, meshEnd, dir);
}

template<class T, int D, typename uiter>
void gsTensorDomainBoundaryIterator<T,D,uiter>::next(index_t increment)
{
    bool isGood(true);
    for (index_t i = 0; i < increment; i++)
        isGood = isGood && nextLexicographicIter(curElement, meshEnd, dir);
}

template<class T, int D, typename uiter>
gsVector<unsigned, D> gsTensorDomainBoundaryIterator<T,D,uiter>::index() const
{
    gsVector<unsigned, D> curr_index(d);
    for (int i = 0; i < dir; ++i)
        curr_index[i]  = curElement[i] - meshStart[i];
    for (int i = dir+1; i < d; ++i)
        curr_index[i]  = curElement[i] - meshStart[i];
    curr_index[dir]  = tindex;
    return curr_index;
}

template<class T, int D, typename uiter>
gsVector<T> gsTensorDomainBoundaryIterator<T,D,uiter>::lowerCorner() const
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

template<class T, int D, typename uiter>
gsVector<T> gsTensorDomainBoundaryIterator<T,D,uiter>::upperCorner() const
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

template<class T, int D, typename uiter>
const T gsTensorDomainBoundaryIterator<T,D,uiter>::getPerpendicularCellSize() const
{
    return curElement[dir].upperCorner().value() - curElement[dir].lowerCorner().value();
}

template<class T, int D, typename uiter>
void gsTensorDomainBoundaryIterator<T,D,uiter>::adjacent( const gsVector<bool> & /* orient */,
               gsDomainIterator<T>  & /* other */ )
{
    GISMO_NO_IMPLEMENTATION
}

template<class T, int D, typename uiter>
void gsTensorDomainBoundaryIterator<T,D,uiter>::setBreaks(const std::vector<T> & newBreaks, index_t i) // i: direction
{
    GISMO_ASSERT(i!=dir, "Changing non-boundary breakpoints is not supported.");
    meshStart[i] = give(gsBreaksIterator<T>::make(newBreaks, true));
    curElement[i] = give(gsBreaksIterator<T>::make(newBreaks, true));
    meshEnd[i]   = give(gsBreaksIterator<T>::make(newBreaks, false));
}

} // namespace gismo
