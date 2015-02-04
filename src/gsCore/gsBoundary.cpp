/** @file gsBoundary.cpp

    @brief Provides implementation of functions related to interfaces and boundaries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/


#include <gsCore/gsBoundary.h>
#include <gsTensor/gsTensorTools.h>

namespace gismo {


void boxSide::getContainedCorners (int dim, std::vector<boxCorner> &corners) const
{
    GISMO_ASSERT(dim>=0, "Dimension must be non negative");
    corners.resize(0);
    corners.reserve(1<<(dim-1));
    const index_t dir = direction();
    const bool    par = parameter();
    for (boxCorner c=boxCorner::getFirst(dim); c<boxCorner::getEnd(dim);++c)
    {
        if (c.parameters(dim)(dir) == par)
            corners.push_back(c);
    }
}


void boundaryInterface::matchDofs(gsVector<int>    bSize,
                                  gsMatrix<unsigned> & b1,
                                  const gsMatrix<unsigned> & b2) const
{
    if (b1.size() == 1)
        return;

    // Get structure of the interface dofs
    const index_t s1 = first() .direction();
    const index_t s2 = second().direction();
    const gsVector<bool>    & dirOr = dirOrientation();
    const gsVector<index_t> & bMap  = dirMap();
    const index_t  d = bMap.size();

    gsVector<int>  bPerm(d-1);
    index_t c = 0;
    for (index_t k = 0; k<d; ++k )
    {
        if ( k == s1 ) // skip ?
            continue;

        if ( ! dirOr[k] ) // flip ?
            flipTensorVector(c, bSize, b1);
        
        bPerm[c] = ( bMap[k] < s2 ? bMap[k] : bMap[k]-1 );
        c++;
    }
    
    // Permute
    permuteTensorVector<unsigned,-1>(bPerm, bSize, b1);
}


} //namespace gismo


