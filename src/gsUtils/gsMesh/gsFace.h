/** @file gsFace.h

    @brief Provides gsFace class for a face of a gsMesh.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D. Mayer
*/

#pragma once

#include <gsUtils/gsMesh/gsMeshElement.h>
#include <gsCore/gsLinearAlgebra.h>


namespace gismo 
{

template <class T> 
class gsFace  : public gsMeshElement<T>
{
public:
    /// Shared pointer for gsFace
    typedef memory::shared_ptr< gsFace > Ptr;

    /// Unique pointer for gsFace
    typedef memory::unique_ptr< gsFace > uPtr;

    typedef gsMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsVertexHandle gsVertexHandle;
    typedef typename MeshElement::gsEdgeHandle gsEdgeHandle;
    typedef typename MeshElement::gsFaceHandle gsFaceHandle;

public:
    gsFace() : MeshElement() { }

    virtual ~gsFace() { }

    gsFace(std::vector<gsVertexHandle> const & vert ) : MeshElement() 
        { 
            vertices = vert;
            for ( typename std::vector<gsVertexHandle>::iterator 
                  it = vertices.begin(); it!= vertices.end(); ++it)
	      (*it)->addFace(this);
        }

    gsFaceHandle handle() { return static_cast<gsFaceHandle>(this); }

    gsFace(gsVertexHandle const & v0, gsVertexHandle const & v1, gsVertexHandle const & v2) : MeshElement() 
        { 
            vertices.push_back(v0);
            vertices.push_back(v1);
            vertices.push_back(v2);
            v0->addFace(this);
            v1->addFace(this);
            v2->addFace(this);
        }

    gsFace(gsVertexHandle const & v0, gsVertexHandle const & v1, gsVertexHandle const & v2, gsVertexHandle const & v3) : MeshElement() 
        { 
            vertices.push_back(v0);
            vertices.push_back(v1);
            vertices.push_back(v2);
            vertices.push_back(v3);
            v0->addFace( handle() );
            v1->addFace(this);
            v2->addFace(this);
            v3->addFace(this);
        }

    // clone function
    //GISMO_CLONE_FUNCTION(gsFace)
    uPtr clone() const { return uPtr(new gsFace(*this)); }

    bool operator< (gsFace const & rhs) const
    {
        return ( Xless<T>(this->vertices[0],rhs.vertices[0])||
           ( (this->vertices[0]->x() == rhs.vertices[0]->x() &&
              this->vertices[0]->y() == rhs.vertices[0]->y() &&
              this->vertices[0]->z() == rhs.vertices[0]->z() ) &&
                Xless<T>(this->vertices[1],rhs.vertices[1]) )||
           ( (this->vertices[0]->x() == rhs.vertices[0]->x() &&
              this->vertices[0]->y() == rhs.vertices[0]->y() &&
              this->vertices[0]->z() == rhs.vertices[0]->z() &&
              this->vertices[1]->x() == rhs.vertices[1]->x() &&
              this->vertices[1]->y() == rhs.vertices[1]->y() &&
              this->vertices[1]->z() == rhs.vertices[1]->z() ) &&
                Xless<T>(this->vertices[2],rhs.vertices[2]) ) );
    }
    bool operator == (gsFace const & rhs) const
    {
        return(
           this->vertices[0]->x() == rhs.vertices[0]->x() &&
           this->vertices[0]->y() == rhs.vertices[0]->y() &&
           this->vertices[0]->z() == rhs.vertices[0]->z() &&
           this->vertices[1]->x() == rhs.vertices[1]->x() &&
           this->vertices[1]->y() == rhs.vertices[1]->y() &&
           this->vertices[1]->z() == rhs.vertices[1]->z() &&
           this->vertices[2]->x() == rhs.vertices[2]->x() &&
           this->vertices[2]->y() == rhs.vertices[2]->y() &&
           this->vertices[2]->z() == rhs.vertices[2]->z() );
    }
    bool operator != (gsFace const & rhs)const
    {
        return !(*this==rhs);
    }
    inline void addVertex(gsVertexHandle const & v)
        {
            vertices.push_back( v );
            v->addFace(this);
            return v;
        }

    void move(scalar_t const&  dx, scalar_t const&  dy, scalar_t const&  dz) 
        {
            //for ( typename 
            //it->move(dx,dy,dz)
        }

    std::ostream &print(std::ostream &os) const
        {
            os<<"gsFace: ";
            for ( typename std::vector<gsVertexHandle>::const_iterator
                  it = vertices.begin(); it!= vertices.end(); ++it)
                os<< (*it)->getId()<<" ";
            os<<"\n";
            return os;
        }
    

     gsVector3d<T> orthogonalVector()
        {
         gsVector3d<T> result;
         gsVector3d<T> lhs(((*vertices[1]).x()-(*vertices[0]).x()),
                 ((*vertices[1]).y()-(*vertices[0]).y()),
                 ((*vertices[1]).z()-(*vertices[0]).z()));
         gsVector3d<T> rhs(((*vertices[2]).x()-(*vertices[0]).x()),
                 ((*vertices[2]).y()-(*vertices[0]).y()),
                 ((*vertices[2]).z()-(*vertices[0]).z()));
         result = rhs.cross(lhs);
         return result;
        }

public:

    //List of vetrices on this face
    std::vector<gsVertexHandle> vertices;
    std::vector<gsFaceHandle> nFaces;
    std::vector<gsEdgeHandle> nEdges;
    int faceIdentity;
};

} // namespace gismo


