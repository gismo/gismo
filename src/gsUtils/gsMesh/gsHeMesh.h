/** @file gsHeMesh.h

    @brief Provides gsHeMesh class - half-edge data structure

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen
*/

// Filetype
// http://www.riken.jp/brict/Yoshizawa/Research/PLYformat/PLYformat.html

#ifndef _MESH_H_
#define _MESH_H_

#include <set>

#include <gsUtils/gsMesh/gsBoundingBox.h>
#include <gsUtils/gsMesh/gsHeVertex.h>
#include <gsUtils/gsMesh/gsHalfEdge.h>
#include <gsUtils/gsMesh/gsHalfFace.h>
#include <gsUtils/gsMesh/gsCell.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo {

template <class T >
class gsHeMesh : public gsHeMeshElement<T>
{
public:
    typedef gsHeMeshElement<T> MeshElement;
    typedef typename gsHeMeshElement<T>::scalar_t scalar_t;
    typedef typename gsHeMeshElement<T>::gsHeVertexHandle gsVertexHandle;
    typedef typename gsHeMeshElement<T>::gsHalfEdgeHandle gsHalfEdgeHandle;
    typedef typename gsHeMeshElement<T>::gsHalfFaceHandle gsHalfFaceHandle;
    typedef typename gsHeMeshElement<T>::gsCellHandle gsCellHandle;
    typedef gsVector3d<T>* gsVectorHandle;
    typedef gsMatrix<T> gsMatrixT;
    
public:

    std::set< gsVector<int>  > halfFace_indices;

    int numVertices;
    int numHalfEdges;
    int numHalfFaces;
    int numCells;

    //std::map<int, gsVertexHandle> vertices;
    std::vector<gsVertexHandle > vertex;
    std::vector<gsHalfEdgeHandle >  edge;
    std::vector<gsHalfFaceHandle >  face;
    std::vector<gsCellHandle >      cell;

    gsBoundingBox<scalar_t> bb;

    bool manifold;
    bool initialized;    

public:
    gsHeMesh(int const & i = 0) : MeshElement(i), initialized(false) 
        {
            numHalfFaces = 0;
            numVertices = 0;
            numCells = 0;
            numHalfEdges = 0;
             
            manifold = true;
        };
        
    ~gsHeMesh();
    
    void addVertex(gsVertexHandle v);   

    // add coords to gsHeVertex, not yet hed
    gsVertexHandle addVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z=0)
        {
            gsVertexHandle v = this->makeHeVertex(x,y,z);
            vertex.push_back(v );
            v->setId(numVertices++);
            return v;
        };

    void addHalfEdge(gsHalfEdgeHandle he)
        {
            edge.push_back(he);
            he->setId(numHalfEdges++);
        };

    void addEdge(gsVertexHandle const& v1, gsVertexHandle const& v2)
        {
            gsHalfEdgeHandle e1 = makeHalfEdge(v1);
            gsHalfEdgeHandle e2 = makeHalfEdge(v2);
            e1->mate = e2;
            e2->mate = e1;
            edge.push_back( e1 );
            e1->setId(numHalfEdges++);
            edge.push_back( e2 );
            e2->setId(numHalfEdges++);
        };

    gsHalfFaceHandle addHalfFace(std::vector<gsVertexHandle> V);
    
    //------------------------------------------------------------------------------------------------------
    void addCell(gsCellHandle *c) 
        {
            cell.push_back(c);
            c->setId(numCells++);
        };
    
    //------------------------------------------------------------------------------------------------------
    gsHalfEdgeHandle findHalfEdge( gsVertexHandle const & s)
        {
            gsHalfEdgeHandle h = 0;
            for ( typename std::set<gsHalfEdgeHandle>::iterator 
                      it = edge.begin(); it!= edge.end(); ++it)
                if ( it->source == s ) 
                {
                    h = *it;
                    break;
                }
            return h;
        };
    
	   std::ostream &print(std::ostream &os) const
        {
            os<<"gsHeMesh with "<<numVertices<<" vertices and "<<numHalfFaces<<" half-faces.\n";
            for ( typename std::vector<gsVertexHandle>::const_iterator
                      it = vertex.begin(); it!= vertex.end(); ++it)
                os<< **it ;
            return os;
        };

    //------------------------------------------------------------------------------------------------------
    // Set up default values for the spline information of each trimming curve attached to each HE	
    void setDefaultTrimmingLine();	
    
    //------------------------------------------------------------------------------------------------------
    // Assigning mates for each HE	
    void setHeMate();	    
    
    //------------------------------------------------------------------------------------------------------
    //generating full control points from four corners of a patch 
    gsMatrixT* simpleFaceCP(std::vector< gsVector3d<T>* > corner, gsKnotVector<T>* kv1, gsKnotVector<T>* kv2);    

};
 
//////////////////////////////////////////////////////
// Source
//////////////////////////////////////////////////////

template<class T>
void gsHeMesh<T>::addVertex(gsVertexHandle v){
    vertex.push_back(v);
    if(!initialized) {
	bb.pMax.x() = bb.pMin.x() = v->x();
	bb.pMax.y() = bb.pMin.y() = v->y();
	bb.pMax.z() = bb.pMin.z() = v->z();
	initialized = true;
    }
    else {
	if (bb.pMax.x() < v->x())
	    bb.pMax.x() = v->x();
	if (bb.pMin.x() > v->x())
	    bb.pMin.x() = v->x();
	
	if (bb.pMax.y < v->y())
	    bb.pMax.y = v->y();
	if (bb.pMin.y > v->y())
	    bb.pMin.y = v->y();
	
	if (bb.pMax.z() < v->z())
	    bb.pMax.z() = v->z();
	if (bb.pMin.z() > v->z())
	    bb.pMin.z() = v->z();
    }
    numVertices++;
};

template<class T>
gsHeMesh<T>::~gsHeMesh()
{
    //apaga os vertices

    typename std::vector<gsVertexHandle>::iterator vIter;
    for(vIter = vertex.begin(); vIter != vertex.end(); vIter++) {
        delete *vIter;
    }
    vertex.clear();

    //apaga as arestas
    typename std::vector<gsHalfEdgeHandle>::iterator eIter;
    for(eIter = edge.begin(); eIter != edge.end(); eIter++) {
        delete *eIter;
    }
    edge.clear();
    
    //apaga as faces
    typename std::vector<gsHalfFaceHandle>::iterator fIter;
    for(fIter = face.begin(); fIter != face.end(); fIter++) {
        delete *fIter;
    }
    face.clear();
};

    
template<class T>
void gsHeMesh<T>::setHeMate()
    {      
      // Method in the meantime: consider each pair of halfedges and see if they are mates
      // Todo: consider each pair of faces 
      unsigned int noMate(0); // number of mates 
      for (typename std::vector<gsHalfEdgeHandle>::iterator it1=edge.begin(); it1!=edge.end()-1; ++it1)
      {
	for (typename std::vector<gsHalfEdgeHandle>::iterator it2=it1+1; it2!=edge.end(); ++it2)
	{
	  // a pair of half edges are mates iff the source of one of them is the target of the orther
      gsVector3d<T> source1, source2, target1, target2;
	  source1 = (*(*(*it1)).source).coords;
	  source2 = (*(*(*it2)).source).coords;
	  target1 = (*it1)->next->source->coords;
	  target2 = (*it2)->next->source->coords;
	  if (source1[0]==target2[0] && source1[1]==target2[1] && source1[2]==target2[2]
	    && source2[0]==target1[0] && source2[1]==target1[1] && source2[2]==target1[2])
	  {
	    noMate++;
	    (*it1)->mate = *it2;
	    (*it2)->mate = *it1;
	  }
	}
      }	  
      // check if the number of mates is the same as the number of assignments
      if (2*noMate!=edge.size())
      {
	std::cout <<"\n"<<"The number of assignments of HE mates (="<< noMate <<") is NOT equal to half of the HE number (="<< edge.size()/2 <<"), this is most likely because of the wrong order of the vertices of a face "<<"\n";
	exit(1);
      }
    }; 

template<class T>
typename gsHeMesh<T>::gsHalfFaceHandle gsHeMesh<T>::addHalfFace(std::vector<gsVertexHandle> V)
{
    gsHalfFaceHandle f = new gsHalfFace<T>();
    //gsHalfEdgeHandle he;
    std::vector<gsHalfEdgeHandle> H;
    for ( typename std::vector<gsVertexHandle>::iterator
          it = V.begin(); it!= V.end(); ++it)
    {
    // initialize the half-edge associated with the vertex *it and the face f
    gsHalfEdgeHandle tempHe = this->makeHalfEdge(*it,f);
    tempHe->setId(numHalfEdges++);
    H.push_back( tempHe );
    edge.push_back( tempHe );
    // assign the incident halfedge associated with the vertex *it to *it
    // fist check if hed is already assigned to avoid doubling assignments
    if (!(*it)->hed)
    (*it)->hed = tempHe;
    }
    f->setBoundary( H );
    face.push_back(f);
    f->setId(numHalfFaces++);
    return f;
}

};// namespace gismo
#endif
