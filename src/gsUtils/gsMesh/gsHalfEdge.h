/** @file gsHalfEdge.h

    @brief Provides gsHalfEdge class for an edge of a gsHeMesh

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

//http://www.cs.purdue.edu/homes/cmh/distribution/books/geo.html
// http://www.cs.utah.edu/~xchen/euler-doc/

#ifndef _HALFEDGE_H_
#define _HALFEDGE_H_

#include <gsUtils/gsMesh/gsHeMeshElement.h>
#include <gsUtils/gsMesh/gsHalfEdge.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo {


enum gsHEdSign{ PLUS=true, MINUS=false };

template <class T>
class gsHalfEdge : public gsHeMeshElement<T>
{
public:
    typedef gsHeMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsHeVertexHandle gsVertexHandle;
    typedef typename MeshElement::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename MeshElement::gsHalfFaceHandle gsFaceHandle;
    typedef gsVector3d<T> gsVector;
    typedef gsVector * gsVectorHandle;     

public:
    gsHalfEdge(gsVertexHandle const & v, gsFaceHandle const & f, int const& i, bool const& conv) : MeshElement(i), source(v), face(f), is_convex(conv) 
        { 
            mate = 0;
            next = 0;
            prev = 0;	    
        };

    gsHalfEdge() { };
    virtual ~gsHalfEdge(){ };

    std::ostream &print(std::ostream &os) const
        {
            //os<<"gsHalfEdge from "<< *source ;
            return os;
        };

public:

    //gsEdgeHandle edge;#
    gsHalfEdgeHandle mate;

    gsVertexHandle source;
    gsHalfEdgeHandle next, prev;

    gsFaceHandle face;// parent face
    
    // Additional properties
    bool is_convex;
    gsKnotVector<T>* kv; // Knot vector and degree of the trimming line (in parameter domain, not the trimming curve in the trimming surface)
    gsMatrix<T>* cp; // Control points of  the trimming line     
    
};

// template <class T>
// addHe(typename gsHeMeshElement<Vertex>::gsVertexHandle *v, gsHalfEdge<T> *h, gsHEdSign sign)
//  {
//     gsHalfEdge<T> *he;
//     if (h->edge == 0)
//         he = h;
//     else {
//         he = new gsHalfEdge<T>();
//         h->prev->next = he;
//         he->prev = h->prev;
//         h->prev = he;
//         he->next = h;
//     }

//     he->edge = e;
//     he->origin = v;
//     he->loop = h->loop;
//     if(sign == PLUS)
//         e->hed1 = he;
//     else 
//         e->hed2 = he;

//     return he;
// }

// template <class T>
// gsHalfEdge<T>* delHe(gsHalfEdge<T> *he)
//  {
//     if(he->edge == 0) {
//         delete he;
//     } else if (he->next == he) {
//         he->edge = 0;
//         return he;
//     }
//     else {
//         he->prev->next = he->next;
//         he->next->prev = he->prev;
//         delete he;
//         return he->prev;
//     }
// }

};//namespace gismo

#endif
