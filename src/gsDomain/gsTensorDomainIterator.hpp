/** @file gsTensorDomainIterator.hpp

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsTensorDomainIterator.h>

namespace gismo
{

template<class T, int D>
gsTensorDomainIterator<T,D>::gsTensorDomainIterator(const gsTensorDomain<D,T> & domain, index_t patchId)
    : gsDomainIterator<T>(0, patchId, boundary::none) // Call base constructor with patch ID
{
    for (int i=0; i < D; ++i)
    {
        meshEnd[i]    = give(domain.component(i)->endAll()  );
        meshStart[i]  = give(domain.component(i)->beginAll());
        curElement[i] = give(domain.component(i)->beginAll());
    }
}

template<class T, int D>
gsTensorDomainIterator<T,D>::gsTensorDomainIterator(const gsTensorDomainIterator<T,D> & other) = default;

template<class T, int D>
typename gsTensorDomainIterator<T,D>::domainIter gsTensorDomainIterator<T,D>::clone() const { return domainIter(new gsTensorDomainIterator(*this)); }

template<class T, int D>
void gsTensorDomainIterator<T,D>::next()
{
    nextLexicographicIter(curElement, meshEnd);
}

template<class T, int D>
void gsTensorDomainIterator<T,D>::next(index_t increment)
{
    bool isGood(true);
    for (index_t i = 0; i < increment; i++)
        isGood = isGood && nextLexicographicIter(curElement, meshEnd);
}

template<class T, int D>
void gsTensorDomainIterator<T,D>::reset()
{
    for (index_t i = 0; i < D; ++i)
        curElement[i].reset();
}

template<class T, int D>
gsVector<unsigned, D> gsTensorDomainIterator<T,D>::index() const
{
    gsVector<unsigned, D> curr_index(D);
    for (int i = 0; i < D; ++i)
        curr_index[i]  = curElement[i].id();
    return curr_index;
}

template<class T, int D>
void gsTensorDomainIterator<T,D>::getVertices(gsMatrix<T>& result)
{
    result.resize( D, 1 << D);

    const gsVector<T> lower = lowerCorner();
    const gsVector<T> upper = upperCorner();
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

template<class T, int D>
gsVector<T> gsTensorDomainIterator<T,D>::lowerCorner() const
{
    gsVector<T> lower(D);
    for (short_t i = 0; i < D ; ++i)
        lower[i]  = curElement[i].lowerCorner().value();
    return lower;
}

template<class T, int D>
gsVector<T> gsTensorDomainIterator<T,D>::upperCorner() const
{
    gsVector<T> upper(D);
    for (short_t i = 0; i < D ; ++i)
        upper[i]  = curElement[i].upperCorner().value();
    return upper;
}

template<class T, int D>
bool gsTensorDomainIterator<T,D>::isBoundaryElement() const
{
    for (int i = 0; i< D; ++i)
        if ( curElement[i].isBoundaryElement() )
            return true;
    return false;
}

template<class T, int D>
index_t gsTensorDomainIterator<T,D>::domainDim() const {return D;}

} // namespace gismo
