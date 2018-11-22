/** @file gsSolidHeVertex.h

    @brief Provides gsSolidHeVertex - a vertex of a gsSolid.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

#pragma once

#include <gsModeling/gsSolidElement.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo {

template <class T>  class gsSolidHalfEdge;

template <class T> 
class gsSolidHeVertex  : public gsSolidElement<T>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef T scalar_t;
    typedef gsSolidElement<T> SolidElement;
    typedef typename SolidElement::gsSolidHeVertexHandle gsSolidHeVertexHandle;
    typedef typename SolidElement::gsSolidHalfEdgeHandle gsSolidHalfEdgeHandle;
    typedef typename SolidElement::gsSolidHalfFaceHandle gsSolidHalfFaceHandle;    
    
// Data members    
public:
    gsVector3d<T> coords; // TODO: should be inherited instead, makes it more naturally
    // halfedge out of this vertex
    gsSolidHalfEdgeHandle hed;    
    
public:
    gsSolidHeVertex(scalar_t x, scalar_t y, scalar_t z = 0) : SolidElement(), coords(x,y,z), hed(0) { }
    gsSolidHeVertex(scalar_t x, scalar_t y, scalar_t z, int i) : SolidElement(i), coords(x,y,z), hed(0) { }

    virtual ~gsSolidHeVertex(){ }
    
    /// translate the vertex towards a displacement verctor (dx,dy,dz)
    void move(scalar_t dx, scalar_t dy, scalar_t dz) {coords.x() += dx;coords.y() += dy;coords.z() += dz;}
    
//--------------------------------------------------------------------------
// const members
public:
    gsVector3d<T> getCoordinate() const {return coords;}
    T   x () const { return coords(0); }
    T   y () const { return coords(1); }
    T   z () const { return coords(2); }
    
    std::ostream &print(std::ostream &os) const; 
    
    /// Getting pointers to all faces incident to the current vertex
    std::vector<gsSolidHalfFaceHandle> getHalfFaces() const;    
    
    /// Getting pointer to the HE connecting two given vertex IDs, if the output is a NULL pointer then there is no such a HE
    gsSolidHalfEdgeHandle getHalfEdge(gsSolidHeVertexHandle anotherVertex) const;
    
    /// determine if there is a half-edge from this vertex to \a anotherVertex
    bool hasHalfEdge(gsSolidHeVertexHandle anotherVertex) const;
    
    /// Getting pointers to the faces that contains two vertices given by their IDs
    std::vector< gsSolidHalfFaceHandle > getFacesContaining2Vertices(gsSolidHeVertexHandle anotherVertex, bool const & IsBoundaryHalfEdge=false) const;          
    
    /// Return pointer to the half-edge on a specified face that has this vertex as the source (or destination, if test==true)
    gsSolidHalfEdgeHandle getHalfEdgeOnFace(gsSolidHalfFaceHandle f, bool dest);

    /// Collect all pointers that emanate from this vertex (has this vertex as their sources)
    std::vector<gsSolidHalfEdgeHandle> halfEdges() const;

    /// Check if the current HeVertex has similar coordinates to other HeVertex \param other within some tolerance
    /// Please use this member with caution
    bool isEquiv(gsSolidHeVertexHandle other, T tol) const
    { using std::abs;if ( abs(x()-other->x())<=tol && abs(y()-other->y())<=tol && abs(z()-other->z())<=tol  ) return true; return false; }
};

//=============================================================================
// SOURCE
//=============================================================================
template <class T>
std::ostream &gsSolidHeVertex<T>::print(std::ostream &os) const 
{
    os<<"Vertex " << this->getId() << "( "<<coords.x() << " " << coords.y() << " " << coords.z() << " )\n";
    return os;      
}

template <class T>
std::vector<typename gsSolidHeVertex<T>::gsSolidHalfFaceHandle> gsSolidHeVertex<T>::getHalfFaces() const
{
  // See gsSolid<T>::getHEwrtTwoVertexIDs for algorithm explaination
  std::vector< gsSolidHalfFaceHandle > faceList;
  gsSolidHalfEdgeHandle he1 = hed;
  gsSolidHalfEdgeHandle instant_he = he1;
  do
  {
    faceList.push_back(instant_he->face);
    instant_he = instant_he->prev->mate;
  } while (instant_he!= he1);
  return faceList;  
}

template <class T>
typename gsSolidHeVertex<T>::gsSolidHalfEdgeHandle gsSolidHeVertex<T>::getHalfEdge(gsSolidHeVertexHandle anotherVertex) const
{
  // Theoretical background: ID1: the current vertex, ID2: = anotherVertex. Let he1 is the HE emanating from ID1, then he1->prev and he1 are the two HEs of 1 face containing ID1. 
  // Furthermore: he1->prev->mate and he1->prev->mate->pre are the two HEs of the adjacent face containing ID1 in a CCW fashion.
  // Continue to do that way, we can access all HEs emanating or shooting at ID1
  const gsSolidHalfEdgeHandle he1 = hed;
  gsSolidHalfEdgeHandle instant_he = he1;
  do
  {
    // check if target of instant_he is ID2
    if (instant_he->target()==anotherVertex)
      return instant_he;
    instant_he = instant_he->prev->mate;
  } while (instant_he!=he1);

  GISMO_ERROR("ERROR:gsSolidHeVertex.h: No HE is found while it is supposed to be existent");
}

template <class T>
bool gsSolidHeVertex<T>::hasHalfEdge(gsSolidHeVertexHandle anotherVertex) const
{
  gsSolidHalfEdgeHandle he1 = hed;
  gsSolidHalfEdgeHandle instant_he = he1;
  do
  {
    if (instant_he->target()==anotherVertex) {return true;};
    instant_he = instant_he->prev->mate;
  } while (instant_he!=he1);
  return false;
}

template <class T>
std::vector< typename gsSolidHeVertex<T>::gsSolidHalfFaceHandle > 
// *isBoundaryHalfEdge* = is an existing half-edge
gsSolidHeVertex<T>::getFacesContaining2Vertices(gsSolidHeVertexHandle anotherVertex, bool const & IsBoundaryHalfEdge) const
{
  std::vector<gsSolidHalfFaceHandle> faceList;
  if (IsBoundaryHalfEdge)
  {
    gsSolidHalfEdgeHandle he = getHalfEdge(anotherVertex);
    // as this is a boundary edge, there are exactly two faces incident to the edge
    faceList.push_back(he->face);
    faceList.push_back(he->mate->face);
  }
  else
  {
    // as this is NOT a boundary edge, there are only ONE face containing the edge
    std::vector<gsSolidHalfFaceHandle> faceList1;
    std::vector<gsSolidHalfFaceHandle> faceList2;
    faceList1  = getHalfFaces();
    faceList2 = anotherVertex->getHalfFaces();
    // find intersection of faceList and faceList2
    for (unsigned i=0;i!=faceList1.size();i++)
    {
      for (unsigned j=0;j!=faceList2.size();j++)
      {
	if (faceList1[i]==faceList2[j]) {faceList.push_back(faceList1[i]); break;};
      };
    };
    
  };
  return faceList;
}

template<class T>
typename gsSolidHeVertex<T>::gsSolidHalfEdgeHandle gsSolidHeVertex<T>::getHalfEdgeOnFace(typename gsSolidHeVertex<T>::gsSolidHalfFaceHandle f, bool dest)
{
  gsSolidHalfEdgeHandle currentEdge = hed;
  while(1)
  {
    if(!dest && currentEdge->face == f) return currentEdge;
    currentEdge = currentEdge->mate;
    if(dest && currentEdge->face == f) return currentEdge;
    currentEdge = currentEdge->next;
    assert(currentEdge != hed); // didn't find the face
  }
}

template<class T>
std::vector<typename gsSolidHeVertex<T>::gsSolidHalfEdgeHandle> gsSolidHeVertex<T>::halfEdges() const
{
  gsSolidHalfEdgeHandle he1 = hed;
  gsSolidHalfEdgeHandle instant_he = he1;
  std::vector<gsSolidHalfEdgeHandle> heSet;
  do
  {
    heSet.push_back(instant_he);
    instant_he = instant_he->prev->mate;
  } while (instant_he!=he1);
  return heSet;
}

} // namespace gismo


