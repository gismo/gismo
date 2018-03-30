/** @file gsMeshElement.h

    @brief Provides gsMeshElement class - a vertex, edge, face or cell of a gsMesh

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>


namespace gismo {


template <class T>
class gsMeshElement
{
public:
    typedef T scalar_t;
    typedef gsVertex<T>        * gsVertexHandle;
    typedef gsEdge<T>          * gsEdgeHandle;
    typedef gsFace<T>          * gsFaceHandle;
    typedef gsCell<T>          * gsCellHandle;

    typedef gsHalfEdge<T>      * gsHalfEdgeHandle;
    typedef gsHeVertex<T>      * gsHeVertexHandle;
    typedef gsHalfFace<T>      * gsHalfFaceHandle;

public:
    explicit gsMeshElement(int i = 0) : id(i)
    { }

    virtual ~gsMeshElement() { }
    
    int getId() const   { return id; }
    void setId(int i)   { id=i; }
 
public:

    static gsVertexHandle makeVertex( scalar_t x, scalar_t y, scalar_t z = 0)
    { return new gsVertex<T>(x,y,z); }

    static gsVertexHandle makeVertex( gsVector<T> const & u )
    { return new gsVertex<T>(u); }
    
    static gsFaceHandle makeFace( std::vector<gsVertexHandle> const & vert)
    { return new gsFace<T>(vert); }

    static gsFaceHandle makeFace(gsVertexHandle v0, gsVertexHandle v1, 
                                 gsVertexHandle v2, gsVertexHandle v3)
    { return new gsFace<T>(v0,v1,v2,v3); }

    static gsFaceHandle makeFace(gsVertexHandle v0, gsVertexHandle v1, gsVertexHandle v2)
    { return new gsFace<T>(v0,v1,v2); }
        
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const 
    { os<<"gsMeshElement\n"; return os; }

    friend std::ostream& operator<<(std::ostream& os, const gsMeshElement& e) 
    { return e.print(os); }
    
private:
    int id;
};


} // namespace gismo

