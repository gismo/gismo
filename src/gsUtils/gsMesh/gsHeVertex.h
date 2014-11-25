
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


