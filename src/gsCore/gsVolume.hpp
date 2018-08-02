/** @file gsVolume.hpp

    @brief Provides implementation of Volume common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsPointGrid.h>
#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{


template<class T>
void gsVolume<T>::toMesh(gsMesh<T> & msh, int npoints) const
{
    // OLD CODE NEVER USED gsTensorGridIter has been superseded

    // const gsMatrix<T> param        = this->parameterRange();
    // const gsVector<unsigned> np    = uniformSampleCount(param.col(0), param.col(1), npoints );

    // gsMatrix<T> cp; // Curve point

    // gsTensorGridIter<3,T> grid(param, np);

    // for(; grid.isGood(); grid.next() )
    // {
    //     this->eval_into(grid.currPoint, cp);
    //     msh.addVertex( cp );
    // }

    // gsVector<unsigned> v;
    // v.setZero(3);
    // do
    // {
    //     msh.addFace(v[0], v[0]+1, v[1]+1, v[1]);
    // }
    // while( nextLexicographic(v, np) )
}

template<class T>
gsGeometryEvaluator<T> *
gsVolume<T>::evaluator(unsigned flags) const
{
    switch ( this->coDim() )
    {
    case 0:
        return new gsGenericGeometryEvaluator<T,3,0 >(*this, flags);
    case 1:
        return new gsGenericGeometryEvaluator<T,3,1 >(*this, flags);
    case -2:
        return new gsGenericGeometryEvaluator<T,3,-2>(*this, flags);
    default:
        GISMO_ERROR("Codimension problem.( codim="<<this->coDim() );
    }
}

} // namespace gismo
