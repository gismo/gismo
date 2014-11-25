#ifndef _CELL_H_
#define _CELL_H_

#include <gsUtils/gsMesh/gsMesh.h>
#include <gsUtils/gsMesh/gsMeshElement.h>

namespace gismo {



template <class Vertex>
class gsCell : public gsMeshElement<Vertex>
{
public:
    typedef gsMeshElement<Vertex > MeshElement;
    typedef typename gsMeshElement<Vertex>::scalar_t scalar_t;
    typedef typename gsMeshElement<Vertex>::gsVertexHandle gsVertexHandle;
    typedef typename gsMeshElement<Vertex>::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename gsMeshElement<Vertex>::gsFaceHandle gsFaceHandle;
    typedef typename gsMeshElement<Vertex>::gsCellHandle gsCellHandle;

public:
    gsCell() : MeshElement() { };

    explicit gsCell(int i) : MeshElement(i) { };

    ~gsCell() { };    

public:
    gsFaceHandle boundary;
    gsCellHandle next;
};

};// namespace gismo

#endif
