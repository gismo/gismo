/** @file gsCell.h

    @brief Provides gsCell class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

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
