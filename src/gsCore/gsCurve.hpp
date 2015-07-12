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
#include <gsCore/gsBasis.h>

namespace gismo
{

template<class T> 
int 
gsCurve<T>::degree() const 
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

template<class T>
unsigned gsCurve<T>::functionAtCorner(boxCorner const & c) const
{
    GISMO_ASSERT(c<boxCorner(3),"Invalid endpoint of curve.");
    return ( c == 1 ? 0 : this->coefsSize()-1);
}

template<class T> 
gsGeometryEvaluator<T> *
gsCurve<T>::evaluator(unsigned flags) const
{
    switch ( this->coDim() )
    {
    case 0:
        return new gsGenericGeometryEvaluator<T,1,0 >(*this, flags);
    case 1:
        return new gsGenericGeometryEvaluator<T,1,1>(*this, flags);
    case 2:
        return new gsGenericGeometryEvaluator<T,1,2>(*this, flags);
    default:
        GISMO_ERROR("Codimension problem.");
    }
}


}; // namespace gismo

