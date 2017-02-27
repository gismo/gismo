/** @file gsSolidHalfFace.h

    @brief Provides gsSolidHalfFace - a (half-)face of a gsSolid.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley
*/

#pragma once


#include <gsModeling/gsSolidElement.h>

namespace gismo {



template <class T>
class gsSolidHalfFace : public gsSolidElement<T>
{
public:
    typedef gsSolidElement<T> SolidElement;
    typedef T Scalar_t;
    typedef typename gsSolidElement<T>::gsSolidHeVertexHandle gsSolidHeVertexHandle;
    typedef typename gsSolidElement<T>::gsSolidHalfEdgeHandle gsSolidHalfEdgeHandle;
    typedef typename gsSolidElement<T>::gsSolidHalfFaceHandle gsSolidHalfFaceHandle;
    typedef typename gsSolidElement<T>::gsVolumeHandle gsVolumeHandle;    

public:
    gsVolumeHandle vol; // Point to the volume that the face belongs to
    gsSolidHalfFaceHandle mate;
    // Every half-edge represents an edge loop
    // loop[0] is the outer loop oriented CCW
    std::vector<gsSolidHalfEdgeHandle >  loop;    
    gsTrimSurface<T> * surf;
    
public:
    gsSolidHalfFace() : SolidElement() { vol = 0; }

    explicit gsSolidHalfFace(int i) : SolidElement(i) { vol = 0; }

    virtual ~gsSolidHalfFace() { delete surf; }
    
    /// set values of the members: next, prev, face of the half-edges in the boundary of a face
    void setBoundary(std::vector<gsSolidHalfEdgeHandle> const &  hedges);

//----------------------------------------------------------------------------
// const members
public:  
    /// loop numbers
    int loopN() const {return loop.size();}

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const; 
    
    /// Getting pointers to Vertices of the OUTER loop of the current face    
    std::vector<gsSolidHeVertexHandle> getVertices() const;

    /// Getting pointers to HEs of the OUTER loop a given face
    std::vector<gsSolidHalfEdgeHandle> getHalfEdges() const;  
    
    gsSolidHalfEdgeHandle getHalfEdgeFromLoopOrder(unsigned loopNumber,int curveNumber) const
    {
        gsSolidHalfEdgeHandle he=*(loop.begin()+loopNumber);
        for ( int i=0; i != surf->domain().loop(loopNumber).size(); ++i)
        {
            if (i==curveNumber) break;
            he=he->next;
        }
        return he;
    }
    
    gsSolidHalfEdgeHandle getHalfEdgeFromBoundaryOrder(unsigned curveNumber) const {return getHalfEdgeFromLoopOrder(0, curveNumber);}
    
    unsigned nCurvesOfOneLoop(unsigned loopNumber) { return surf->domain().loop(loopNumber).curves().size(); }
    
    /// find the index, within this face, of edge \a e.
    int indexOfEdge(gsSolidHalfEdgeHandle e) const;
    
    /// find the index, within this face, of vertex \a v.
    int indexOfVertex(gsSolidHeVertexHandle v) const;
    
    /// find out if this face contains a given vertex \a v.
    bool containsVertex(const gsSolidHeVertexHandle v) const;

}; // end class

//=============================================================================
// SOURCE
//=============================================================================
template <class T>
std::ostream &operator<<(std::ostream &os, const gsSolidHalfFace<T>& obj)
{return obj.print(os); }

template <class T>
void gsSolidHalfFace<T>::setBoundary(std::vector<gsSolidHalfEdgeHandle> const &  hedges)
{
  typename std::vector<gsSolidHalfEdgeHandle>::const_iterator p = hedges.begin();
  loop.push_back(*p);
  (*p)->prev = hedges.back();
  (*p)->face = this;
  for ( typename std::vector<gsSolidHalfEdgeHandle>::const_iterator 
          it = hedges.begin()+1; it!= hedges.end(); ++it)
  {   
      (*it)->face = this;
      (*it)->prev = *p;
      (*p)->next  = *it;
      p++;
  }
  (*p)->next = hedges.front();
}

template <class T>
std::ostream &gsSolidHalfFace<T>::print(std::ostream &os) const 
{          
  os<<"\nHalf face ID = "<< this->getId() << " with "; 
  os<<"Coordinates of vertices of face 0: \n";
  std::vector<gsSolidHeVertexHandle> fvert = getVertices();
  for (typename std::vector<gsSolidHeVertexHandle>::const_iterator it=fvert.begin();it!=fvert.end();++it) os<<**it;
  os<<"\n";
  return os;       
}

template <class T>
std::vector<typename gsSolidHalfFace<T>::gsSolidHeVertexHandle> gsSolidHalfFace<T>::getVertices() const
{
  std::vector<gsSolidHeVertexHandle> vert;      
    gsSolidHalfEdgeHandle edge1 = loop[0];
    vert.push_back(edge1->source);
    gsSolidHalfEdgeHandle current_edge=edge1->next;
    while ( current_edge != edge1 )
    {
      vert.push_back(current_edge->source);
      current_edge=current_edge->next;
    }
  return vert;
}

template <class T>
std::vector< typename gsSolidHalfFace<T>::gsSolidHalfEdgeHandle > gsSolidHalfFace<T>::getHalfEdges() const
{
  std::vector<gsSolidHalfEdgeHandle> HEs;      
    gsSolidHalfEdgeHandle edge1 = loop[0];
    HEs.push_back(edge1);
    gsSolidHalfEdgeHandle current_edge=edge1->next;
    while ( current_edge != edge1 )
    {
      HEs.push_back(current_edge);
      current_edge=current_edge->next;
    }
  return HEs;
}

template <class T>
int gsSolidHalfFace<T>::indexOfEdge(gsSolidHalfEdgeHandle e) const
{
  int idx = 0;
  gsSolidHalfEdgeHandle tempEdge = this->loop[0];
  while(tempEdge != e)
  {
    tempEdge = tempEdge->next;
    idx++;
    assert(tempEdge != this->loop[0]); // went round the loop without finding the edge
  }
  return idx;
}

template <class T>
int gsSolidHalfFace<T>::indexOfVertex(gsSolidHeVertexHandle v) const
{
  int idx = 0;
  gsSolidHalfEdgeHandle tempEdge = this->loop[0];
  while(tempEdge->source != v)
  {
    tempEdge = tempEdge->next;
    idx++;
    assert(tempEdge != this->loop[0]); // went round the loop without finding the vertex
  }
  return idx;
}

template <class T>
bool gsSolidHalfFace<T>::containsVertex(const gsSolidHeVertexHandle v) const
{
  gsSolidHalfEdgeHandle tempEdge = this->loop[0];
  do
  {
    tempEdge = tempEdge->next;
    if (tempEdge->source == v) { return true; };
  } while (tempEdge != this->loop[0]);
  return false;
}


} // namespace gismo
