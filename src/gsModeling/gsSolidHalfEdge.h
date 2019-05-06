/** @file gsSolidHalfEdge.h

    @brief Provides gsSolidHalfEdge - a half-edge of a gsSolid

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley, J. Speh
*/

#pragma once

#include <gsModeling/gsSolidElement.h>

namespace gismo {

template <class T>
class gsSolidHalfEdge : public gsSolidElement<T>
{
public:
    typedef gsSolidElement<T> SolidElement;
    typedef typename SolidElement::scalar_t scalar_t;
    typedef gsSolidHeVertex<T> Vertex;
    typedef gsSolidHalfEdge<T> HalfEdge;
    typedef gsSolidHalfFace<T> HalfFace;
/*    typedef gsVector3d<T> gsVector;
    typedef gsVector * gsVectorHandle; */ 
    
public:

    HalfEdge* mate;

    Vertex* source;
    HalfEdge* next;
    HalfEdge* prev;
    HalfFace* face;// parent face
        
    //gsGeometry<T> * trim_curve;    
    bool is_convex;
private:
  int loopNum; // loop number
  
/// Accessors
public:
  int loopN() const {return loopNum;}
  Vertex* target() const {return mate->source;}

public:
    gsSolidHalfEdge() : SolidElement() { eps  = 2.220446049250313e-16; }


    gsSolidHalfEdge(Vertex* v, HalfFace* f, int i, bool conv)
        : SolidElement(i), source(v), face(f), is_convex(conv) 
    { 
        mate = 0;
        next = 0;
        prev = 0;
        loopNum = 0;
        //trim_curve = 0;
        eps = 2.220446049250313e-16;
    }

    gsSolidHalfEdge(Vertex* v, HalfFace* f, int i, bool conv, int loopN) :
        SolidElement(i), source(v), face(f), is_convex(conv), loopNum(loopN)
    {
        mate = 0;
        next = 0;
        prev = 0;
        eps = 2.220446049250313e-16;
    }


    explicit gsSolidHalfEdge(int i) : SolidElement(i) { }

    /// check if two HEs are "equivalent", ie., if their sources have the same coordinates, and their targets have the same coordinates
    bool isEquiv(HalfEdge* other, T tolFactor=1e-8) const
    {
        using std::abs;
        T tol = tolFactor*eps;
        return source->isEquiv(other->source, tol) && target()->isEquiv(other->target(), tol);
    }

    virtual ~gsSolidHalfEdge() { }
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const {      
        os<<"\ngsSolidHalfEdge number: " << this->getId() << " source:" << *source << " target: " << *target() << std::endl;
      return os; }
      
    T eps;
    /// Get the index of the corresponding trimming curve in the gsTrimSurface of the face
    // TODO: get rid of the following, using the indexOfEdge instead
    int trimLoopInd(T tolerance){ GISMO_UNUSED(tolerance); return this->face->indexOfEdge(this); }
    int trimLoopInd();

    /// Moves along edge "n" times
    HalfEdge* moveAlongEdge(int n = 1)
    {
        HalfEdge* edge = this;
        for (int times = 0; times < n; times++)
        {
            edge = edge->next;
        }

        return edge;
    }
};

//================================================================
// Source
//================================================================

/// Print (as string) operator to be used by all mesh elements
template<class T>
std::ostream &operator<<(std::ostream &os, const gsSolidHalfEdge<T>& me)
{return me.print(os); }

template<class T>
int gsSolidHalfEdge<T>::trimLoopInd()
{
    return trimLoopInd(eps);
}


} // namespace gismo
  
