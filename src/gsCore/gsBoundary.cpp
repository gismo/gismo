/** @file gsBoundary.cpp

    @brief Provides implementation of functions related to interfaces and boundaries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/


#include <gsCore/gsBoundary.h>
#include <gsTensor/gsTensorTools.h>

namespace gismo {


void boxSide::getContainedCorners(short_t dim, std::vector<boxCorner> &corners) const
{
    GISMO_ASSERT(dim>=0, "Dimension must be non negative");
    corners.clear();
    corners.reserve( 1ULL<<(dim-1) );
    const short_t dir = direction();
    const bool    par = parameter();
    for (boxCorner c=boxCorner::getFirst(dim); c<boxCorner::getEnd(dim);++c)
    {
        if (c.parameters(dim)(dir) == par)
            corners.push_back(c);
    }
}

void patchSide::getContainedCorners(short_t dim, std::vector<patchCorner> &corners) const
{
    std::vector<boxCorner> tmp;
    boxSide::getContainedCorners(dim, tmp);
    corners.clear();
    corners.reserve(tmp.size());
    for (std::vector<boxCorner>::iterator it=tmp.begin(); it<tmp.end(); ++it)
        corners.push_back( patchCorner(patch, *it) );
}

boxComponent::boxComponent( boxSide b, short_t total_dim ) : m_total_dim(total_dim)
{
    const index_t d = (b.m_index-1)/2;
    const index_t o = (b.m_index-1)%2;
    //m_index = (o+1)*3^d
    m_index = o+1;
    for (index_t c = 0; c<d; ++c)
        m_index *= 3;
}

boxComponent::boxComponent( boxCorner b, short_t total_dim ) : m_total_dim(total_dim)
{
    b.m_index -= 1;
    m_index = 0;
    index_t c = 1;
    for (short_t d=0; d<total_dim; ++d)
    {
        m_index += (1+b.m_index%2) * c;
        c *= 3;
        b.m_index /= 2;
    }
}

short_t boxComponent::dim() const
{
        index_t result = 0, tmp = m_index;
        for (short_t i=0; i<m_total_dim; ++i)
        {
            if (tmp%3 == 0)
                ++result;
            tmp /= 3;
        }
        return result;
}

std::vector<boxCorner> boxComponent::containedCorners() const
{
    const index_t result_sz = 1u<<dim();

    std::vector<boxCorner> result;
    result.reserve(result_sz);

    for (index_t i=0; i<result_sz; ++i)
    {
            index_t idx = 0;
            index_t ii=i, jj=m_index, c=1;
            for (index_t d=0; d<m_total_dim; ++d)
            {
                if (jj%3==0)
                {
                    idx += c*(ii%2);
                    ii /= 2; jj /= 3;
                }
                else
                {
                    idx += c*(jj%3-1);
                    jj /= 3;
                }
                c *= 2;
            }
            result.push_back(boxCorner(1+idx));
    }
    return result;
}

std::vector<boxSide> boxComponent::containingSides() const
{
    const index_t d = dim();

    std::vector<boxSide> result;
    result.reserve(d);
    index_t data = m_index;
    for (index_t i=0; i<d; ++i)
    {
        const index_t loc = data%3;
        if (loc)
            result.push_back(boxSide(loc+2*d));
        data /= 3;
    }
    return result;
}

boxSide boxComponent::asSide() const
{
    GISMO_ENSURE( dim() == m_total_dim-1, "This is not a side." );
    index_t d = 0, idx = m_index;
    while (idx>0)
    {
        if (idx%3)
            return boxSide( idx%3 + 2*d );
        idx /= 3;
        ++d;
    }
    return boxSide(0);
}

boxCorner boxComponent::asCorner() const
{
    GISMO_ENSURE( dim() == 0, "This is not a corner." );
    index_t factor = 1, idx = m_index, result = 0;
    while (idx>0)
    {
        result += ( idx%3 - 1 ) * factor;
        idx /= 3;
        factor *= 3;
    }
    return boxCorner(1+result);
}

boxComponent::location boxComponent::locationForDirection(index_t direction) const
{
    GISMO_ASSERT( direction >= 0 && direction < m_total_dim, "Out of bounds" );
    index_t result = m_index;
    for (index_t i=0; i<direction; ++i)
        result /= 3;
    return location(result%3);
}

void boxComponent::setLocationForDirection(index_t direction, boxComponent::location par)
{
    const index_t diff = par - locationForDirection(direction);
    if (diff)
    {
        index_t factor = 1;
        for (index_t i=0; i<direction; ++i)
            factor *= 3;
        m_index += diff * factor;
    }
}

boxComponent boxComponent::opposite() const
{
    boxComponent result(*this);
    for (index_t i=0; i<m_total_dim; ++i)
    {
        location loc = locationForDirection(i);
        if (loc == begin)
            result.setLocationForDirection(i,end);
        else if (loc == end)
            result.setLocationForDirection(i,begin);
        //if loc == interor, do nothing.
    }
    return result;
}

std::vector<patchCorner> patchComponent::containedCorners() const
{
    std::vector<boxCorner> tmp = boxComponent::containedCorners();
    std::vector<patchCorner> result;
    const index_t sz = tmp.size();
    result.reserve(sz);
    for (index_t i=0; i<sz; ++i)
        result.push_back(patchCorner(m_patch,tmp[i]));
    return result;
}

std::vector<patchSide> patchComponent::containingSides() const
{
    std::vector<boxSide> tmp = boxComponent::containingSides();
    std::vector<patchSide> result;
    const index_t sz = tmp.size();
    result.reserve(sz);
    for (index_t i=0; i<sz; ++i)
        result.push_back(patchSide(m_patch,tmp[i]));
    return result;
}

void boundaryInterface::faceData(gsVector<bool> & flip, gsVector<index_t> & perm) const
{
    const index_t d = directionMap.size();
    flip.resize(d-1);
    perm.resize(d-1);
    const index_t s1 = ps1.direction();
    const index_t s2 = ps2.direction();
    index_t c = 0;

    for (index_t k = 0; k!=s1; ++k)
    {
        flip[c] = directionOrientation[k];
        perm[c] = (directionMap[k] < s2 ? directionMap[k] : directionMap[k]-1);
        c++;
    }
    for (index_t k = s1+1; k!=d; ++k)
    {
        flip[c] = directionOrientation[k];
        perm[c] = (directionMap[k] < s2 ? directionMap[k] : directionMap[k]-1);
        c++;
    }
}

void boundaryInterface::cornerMap(gsVector<index_t> & cmap) const
{
    gsVector<bool>    flip;
    gsVector<index_t> perm;
    faceData(flip, perm);
    cubeIsometry(flip, perm, cmap);
}


void boundaryInterface::reorderCorners(gsMatrix<index_t> & boundary) const
{
    gsVector<index_t> cmap;
    cornerMap(cmap);
    boundary = cmap.asPermutation() * boundary;
}


} //namespace gismo
