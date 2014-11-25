#ifndef _HALFFACE_H_
#define _HALFFACE_H_

#include <list>
#include <gsUtils/gsMesh/gsHeMesh.h>
#include <gsUtils/gsMesh/gsHeMeshElement.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo {
    
template <class T>
class gsHalfFace : public gsHeMeshElement<T>
{
public:
    typedef gsHeMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsHeVertexHandle gsVertexHandle;
    typedef typename MeshElement::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename MeshElement::gsHalfFaceHandle gsHalfFaceHandle;
    typedef typename MeshElement::gsCellHandle gsCellHandle;
    typedef gsVector3d<T> gsVector;
    typedef gsVector * gsVectorHandle;
    
public:

    gsHalfFaceHandle mate;
    gsHalfFaceHandle next; //prev;
    gsHalfEdgeHandle boundary; // todo: make a *loop* member as done with gsSolid
    gsCellHandle parent;    

public:
    gsHalfFace(): MeshElement() { };

    gsHalfFace(gsHalfEdgeHandle const & he, int const& i): MeshElement(i), boundary(he){ };
    
    gsHalfFace(std::vector<gsVertexHandle> const & verts, int const& i): MeshElement(i){}
    
    gsHalfFace(std::vector<gsHalfEdgeHandle> const &  hedges, int const& i)
        : MeshElement(i), boundary(hedges[0])
        {
            this->setBoundary( hedges );
        };

    virtual ~gsHalfFace();

    // set values of the members: next, prev, face of the half-edges in the boundary of a face
    void setBoundary(std::vector<gsHalfEdgeHandle> const &  hedges)
        {
            typename std::vector<gsHalfEdgeHandle>::const_iterator p = hedges.begin();
	    this->boundary = *p;
            (*p)->prev = hedges.back();
            (*p)->face = this;	    
            for ( typename std::vector<gsHalfEdgeHandle>::const_iterator 
                      it = hedges.begin()+1; it!= hedges.end(); ++it)
            {   
                (*it)->face = this;
                (*it)->prev = *p;
                (*p)->next  = *it;
                p++;
            }
            (*p)->next = hedges.front();           	    	 
        };

    std::ostream &print(std::ostream &os) const
        {
            os<<"gsHalfFace \n" ;
            return os;
        };
 
};

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


template <typename T>
gsHalfFace<T>::~gsHalfFace()
{
    // delete boundary half-edges
}
    

};// namespace gismo
#endif
