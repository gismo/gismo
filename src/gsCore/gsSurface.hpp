/** @file gsSurface.hpp

    @brief Provides implementation of Surface common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsTensor/gsGridIterator.h>
#include <gsUtils/gsMesh/gsMesh.h>


namespace gismo
{

template<class T>
void gsSurface<T>::toMesh(gsMesh<T> & msh, int npoints) const
{
    const gsMatrix<T> param     = this->parameterRange();
    gsMatrix<T> cp;
    gsGridIterator<T,CUBE> pIter(param, npoints);
    for(; pIter; ++pIter)
    {
        this->eval_into( *pIter, cp);
        msh.addVertex( cp );
    }

    const typename gsGridIterator<T,CUBE>::point_index np = pIter.numPointsCwise();
    for(index_t j = 0; j + 1 != np[1]; j++)
        for(index_t i= 0; i + 1 != np[0]; i++)
        {
            const index_t ind1 = j * np[0] + i;
            const index_t ind2 = ind1 + np[0];
            //msh.addFace(ind1, ind1+1, ind2+1, ind2);
            msh.addFace(ind1  , ind1+1, ind2+1);
            msh.addFace(ind2+1, ind2  , ind1  );
        }
}

/*
   gsVector<unsigned> v;
    v.setZero(2);
    do
    {
        msh.addFace(v[0], v[0]+1, v[1]+1, v[1]);
    }
    while( nextLexicographic(v, np) );
*/

} // namespace gismo
