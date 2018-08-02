/** @file gsEdge.h

    @brief Provides gsEdge class for an edge of a gsMesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mayer
*/

#pragma once

#include <gsUtils/gsMesh/gsMeshElement.h>

namespace gismo {

template <class T>
class gsEdge : public gsMeshElement<T>
{
public:
    typedef gsMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsVertexHandle gsVertexHandle;
    typedef typename MeshElement::gsEdgeHandle gsEdgeHandle;
    typedef typename MeshElement::gsFaceHandle gsFaceHandle;
    typedef gsVector3d<T> gsVector;
    typedef gsVector * gsVectorHandle;

public:
  gsEdge() { } // : MeshElement()

  gsEdge(gsVertexHandle const & v0, gsVertexHandle const & v1 ):
    source(v0), target(v1)
    {
//      faceId1=0;
//      faceId2=0;
      // next = 0;
      // prev = 0;
    }


  virtual ~gsEdge(){ }

  std::ostream &print(std::ostream &os) const
    {
      os<<"gsEdge from "<< *source<<" to "<<*target ;
      return os;
    };


  void orderVertices()
    {
      if ( Xless<T>(source,target) )
	   std::swap(source,target);
    }

// Lex compare
bool operator< (gsEdge const & rhs) const
{
  return ( Xless<T>(this->source,rhs.source)

          ||
       ( (this->source->x() == rhs.source->x() &&
          this->source->y() == rhs.source->y() &&
          this->source->z() == rhs.source->z() ) &&
         Xless<T>(this->target,rhs.target) ));
}

bool operator == (gsEdge const & rhs) const
{
    return (((this->source->x())==rhs.source->x())&&
            ((this->source->y())==rhs.source->y())&&
            ((this->source->z())==rhs.source->z())&&
            ((this->target->x())==rhs.target->x())&&
            ((this->target->y())==rhs.target->y())&&
            ((this->target->z())==rhs.target->z()));
}
bool operator> (gsEdge const & rhs) const
{
    return !(*this<rhs)&&!(*this==rhs);
}
bool operator != (gsEdge const & rhs) const
{
    return !(*this==rhs);
}

public:

  gsVertexHandle source, target;
  std::vector<gsFaceHandle> nFaces;
  bool sharp;
  std::vector<int> numPatches;
};

} //namespace gismo
