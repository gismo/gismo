/** @file gsHeVertex.h

    @brief Provides gsHeVertex class for a vertex of a gsHeMesh

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsMesh/gsHeMeshElement.h>
#include <gsUtils/gsMesh/gsHalfEdge.h>


namespace gismo {

template <class T>  class gsHalfEdge;

template <class T> 
class gsHeVertex  : public gsVertex<T>
{
public:
    typedef gsVertex<T> Base;
    typedef T scalar_t;
    //typedef typename gsVertex<T>::MeshElement MeshElement;
    typedef typename gsHeMeshElement<T>::gsHalfEdgeHandle gsHalfEdgeHandle;
public:
    gsHeVertex(T x, T y, T z = 0) : Base(x,y,z) { hed=0; };

    virtual ~gsHeVertex(){ };


public:

    // halfedge out of this vertex
    gsHalfEdgeHandle hed;
};

};// namespace gismo


