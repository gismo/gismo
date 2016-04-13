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


void boxSide::getContainedCorners (int dim, std::vector<boxCorner> &corners) const
{
    GISMO_ASSERT(dim>=0, "Dimension must be non negative");
    corners.resize(0);
    corners.reserve( 1U<<(dim-1) );
    const index_t dir = direction();
    const bool    par = parameter();
    for (boxCorner c=boxCorner::getFirst(dim); c<boxCorner::getEnd(dim);++c)
    {
        if (c.parameters(dim)(dir) == par)
            corners.push_back(c);
    }
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


void boundaryInterface::reorderCorners(gsMatrix<unsigned> & boundary) const
{
    gsVector<index_t> cmap;
    boundary = cmap.asPermutation() * boundary;
}


} //namespace gismo


