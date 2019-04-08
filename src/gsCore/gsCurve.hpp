/** @file gsCurve.hpp

    @brief Provides implementation of Curve common operations.

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
short_t gsCurve<T>::degree() const
{
    return this->basis().degree(0);
}

template<class T>
void gsCurve<T>::toMesh(gsMesh<T> & msh, int npoints) const
{
    const gsMatrix<T> param        = this->parameterRange();
    const gsMatrix<T> paramSamples = gsPointGrid<T>(param(0,0), param(0,1), npoints);
    //const bool        zzero        = ( this->geoDim()==2 );

    gsMatrix<T> cp; // Curve point

    this->eval_into(paramSamples.col(0), cp);

    typename gsMesh<T>::VertexHandle v1,
    //v0 = msh.addVertex( cp(0,0), cp(1,0), zzero ? T(0) : cp(2,0) );
    v0 = msh.addVertex( cp );

    for ( int i = 1; i<npoints; ++i)
    {
        this->eval_into(paramSamples.col(i), cp);
        //v1 = msh.addVertex( cp(0,0), cp(1,0), zzero ? T(0) : cp(2,0) );
        v1 = msh.addVertex( cp );
        msh.addEdge(v0, v1);
        v0 = v1;
    }
}

} // namespace gismo
