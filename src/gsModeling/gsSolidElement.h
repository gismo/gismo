/** @file gsSolidElement.h

    @brief Provides gsSolidElement class - interface for an element (vertex,
    edge or face) of a solid

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

#pragma once

#include <gsCore/gsExport.h>

namespace gismo {

template <class T > class gsSolidHeVertex;
template <class T > class gsSolidHalfEdge;
template <class T > class gsSolidHalfFace;
template <class T > class gsVolumeBlock;

template <class T >
class GISMO_DEFAULT_VIS gsSolidElement
{
public:
    typedef T scalar_t;
    typedef gsSolidHeVertex<T>      * gsSolidHeVertexHandle;
    typedef gsSolidHalfEdge<T>      * gsSolidHalfEdgeHandle;
    typedef gsSolidHalfFace<T>      * gsSolidHalfFaceHandle;
    typedef gsVolumeBlock<T>        * gsVolumeHandle;

private:
    int id;
   
public:
    explicit gsSolidElement(int i = 0) : id(i)
    { }

    virtual ~gsSolidElement() { }

    int getId() const   { return id; }
    void setId(int i)   { id=i; }
 
public:

    static gsSolidHeVertexHandle makeHeVertex( scalar_t x, scalar_t y, scalar_t z = 0)
        { return new gsSolidHeVertex<T>(x,y,z); }
//     gsSolidHalfFaceHandle makeFace( std::vector<gsVertexHandle> const & vert)
//         { return new gsFace<T>(vert); }
//     gsSolidHalfFaceHandle makeFace(gsVertexHandle v0, gsVertexHandle v1, gsVertexHandle v2,  gsVertexHandle v3)
//         { return new gsFace<T>(v0,v1,v2,v3); }
//     gsSolidHalfFaceHandle makeFace(gsVertexHandle v0, gsVertexHandle v1, gsVertexHandle v2)
//         { return new gsFace<T>(v0,v1,v2); }
    
    static gsSolidHalfEdgeHandle makeSolidHalfEdge( gsSolidHeVertexHandle source, gsSolidHalfFaceHandle f)
        { return new gsSolidHalfEdge<T>(source,f, 0, true); }
//     gsSolidHalfFaceHandle makeHalfFace( std::vector<gsSolidHalfEdgeHandle> const & hedges)
//         { return new gsHalfFace<T>(hedges); }
    static gsVolumeHandle makeVolume( const std::vector<gsSolidHalfFaceHandle>& hfaces)
        { return new gsVolumeBlock<T>(hfaces); }

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const { os<<"gsSolidElement\n"; return os; }
    
    friend std::ostream& operator<<(std::ostream& os, const gsSolidElement& e) { return e.print(os); }
};


} // namespace gismo
