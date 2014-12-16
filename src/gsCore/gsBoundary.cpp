/** @file gsBoundary.cpp

    @brief Provides implementation of functions related to interfaces and boundaries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/


#include <gsCore/gsBoundary.h>

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

}


