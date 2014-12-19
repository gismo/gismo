/** @file gsHeMeshElement.h

    @brief Provides gsHeMeshElement class - a vertex, edge or face of a gsHeMesh.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <vector>

#include <gsCore/gsForwardDeclarations.h>



namespace gismo {

template <class T > class gsCell;

template <class T > class gsHeVertex;
template <class T > class gsHalfFace;
template <class T > class gsHalfEdge;


template <class T >
class gsHeMeshElement
{
public:
    typedef T scalar_t;

    typedef gsHalfEdge<T>      * gsHalfEdgeHandle;
    typedef gsHeVertex<T>      * gsHeVertexHandle;
    typedef gsHalfFace<T>      * gsHalfFaceHandle;
    typedef gsCell<T>          * gsCellHandle;

public:
    explicit gsHeMeshElement(int i = 0) : id(i)
    { }

    int getId() const   { return id; }
    void setId(int i)   { id=i; }
 
public:
        
    static gsHeVertexHandle makeHeVertex( scalar_t x, scalar_t y, scalar_t z = 0)
    { return new gsHeVertex<T>(x,y,z); }

    static gsHalfEdgeHandle makeHalfEdge( gsHeVertexHandle source, gsHalfFaceHandle f)
    { return new gsHalfEdge<T>(source,f, 0, true); }

    static gsHalfFaceHandle makeHalfFace( const std::vector<gsHalfEdgeHandle> & hedges)
    { return new gsHalfFace<T>(hedges); }

    static gsCellHandle makeCell( const std::vector<gsHalfFaceHandle> & hfaces)
    { return new gsCell<T>(hfaces); }
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const 
    { os<<"gsHeMeshElement\n"; return os; }

    friend std::ostream& operator<<(std::ostream& os, const gsHeMeshElement& e)
    { return e.print(os); }
    
private:
    int id;
};


};// namespace gismo

