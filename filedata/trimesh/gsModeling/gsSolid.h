/** @file gsSolid.h

    @brief Provides declaration of gsSolid class, a boundary-represented solid

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley
*/

#pragma once
 
#include <gsCore/gsLinearAlgebra.h>

#include <gsModeling/gsSolidElement.h>

#include <gsModeling/gsSolidHeVertex.h>
#include <gsModeling/gsSolidHalfEdge.h>
#include <gsModeling/gsSolidHalfFace.h>
#include <gsModeling/gsVolumeBlock.h>
#include <gsUtils/gsMesh/gsBoundingBox.h>

namespace gismo {

/// @brief Class for representing a solid made up of vertices, edges, faces, and volumes.
///
/// \ingroup Modeling
template <class T>
class GISMO_EXPORT gsSolid : public gsSolidElement<T>
{
public: 
    typedef memory::shared_ptr<gsSolid> Ptr;
    typedef memory::unique_ptr<gsSolid> uPtr;
 
    typedef gsSolidElement<T> SolidElement;
    typedef typename gsSolidElement<T>::scalar_t scalar_t;    
    typedef typename gsSolidElement<T>::gsSolidHeVertexHandle gsSolidHeVertexHandle;
    typedef typename gsSolidElement<T>::gsSolidHalfEdgeHandle gsSolidHalfEdgeHandle;
    typedef typename gsSolidElement<T>::gsSolidHalfFaceHandle gsSolidHalfFaceHandle;
    typedef typename gsSolidElement<T>::gsVolumeHandle        gsVolumeHandle;
    typedef gsMatrix<T> gsMatrixT;
    
    /// Iterators
    typedef typename std::vector<gsSolidHalfFaceHandle*>::iterator face_iterator;
    typedef typename std::vector<gsGeometry<T> *>::const_iterator  const_face_iterator;

    typedef typename std::vector<gsSolidHalfEdgeHandle*>::iterator       edge_iterator;
    typedef typename std::vector<gsSolidHalfEdgeHandle*>::const_iterator const_edge_iterator;

    typedef typename std::vector<gsSolidHeVertexHandle*>::iterator       vertex_iterator;
    typedef typename std::vector<gsSolidHeVertexHandle*>::const_iterator const_vertex_iterator;

    typedef typename std::vector<gsVolumeHandle*>::iterator       volume_iterator;
    typedef typename std::vector<gsVolumeHandle*>::const_iterator const_volume_iterator;

public:

    unsigned numVertices;
    unsigned numHalfEdges;
    unsigned numHalfFaces;
    unsigned numVolumes;

    std::vector<gsSolidHeVertexHandle > vertex;
    std::vector<gsSolidHalfEdgeHandle >   edge;
    std::vector<gsSolidHalfFaceHandle >   face;
    std::vector<gsVolumeHandle >        volume;            
    
    bool manifold;
    bool initialized;       
    
protected:
    gsBoundingBox<T> bb;
    bool FaceCCW; // true if vertices of the outer loop of a face are ordered in a CCW fashion when viewing from infinity

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Constructors
public:

    gsSolid(int const & i = 0) : SolidElement(i), initialized(false) 
        {
            numHalfFaces = 0;
            numVertices = 0;
            numVolumes = 0;
            numHalfEdges = 0;             
            manifold = true;
            FaceCCW = true;
        }

    virtual ~gsSolid();

private:    // disable copying (it can be done but is not implemented)
    gsSolid( const gsSolid& );
    gsSolid& operator=( const gsSolid& );
    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Accessors= iterators on vertices, faces, edges, boundingbox
public:
    inline int nVolumes()  	const { return volume.size();}

    inline int nHalfFaces() const { return face.size();}

    inline int nVertices()	const { return vertex.size();}

    inline size_t nHalfEdges()    const { return edge.size();}

    inline bool isFaceCCW()    const { return FaceCCW;}

    inline gsSolidHeVertexHandle getVertexFromID(int const & _id) const 
    {return *(vertex.begin()+_id);}

    inline gsSolidHalfEdgeHandle getHalfEdgeFromID(int const & _id) const 
    {return *(edge.begin()+_id);}

    inline gsSolidHalfFaceHandle getHalfFaceFromID(int const & _id) const 
    {return *(face.begin()+_id);}

    gsBoundingBox<T> getBoundingBox() const {return bb;}
    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Non-constant members
public:
    /// add coords to gsHeVertex, not yet pointers to HEs
    void addHeVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z=0);

    /// add one face as a trimmed surface, the order of the vertices
    /// in *V* must be either CCW or CW when viewing from infinity for
    /// all faces
    gsSolidHalfFaceHandle addFace(std::vector<gsSolidHeVertexHandle> V);

    gsSolidHalfFaceHandle addFace(std::vector<gsSolidHeVertexHandle> V, 
                                  gsTrimSurface<T> * tsurf);

    /// Add a half face to the solid
    //
    /// \param loopV  a vector of vectors of vertices of one loop, the first vector is the outer loop
    /// \param tsurf  a (pointer to) trimmed surface
    gsSolidHalfFaceHandle addFace(std::vector< std::vector<gsSolidHeVertexHandle> > loopV, 
                                  gsTrimSurface<T> * tsurf);

    gsSolidHalfFaceHandle addFace_PlanarPolygon(std::vector<gsSolidHeVertexHandle> V);

    /// add one face as a trimmed surface, the order of the vertices in *V* must be either CCW or CW when viewing from infinity for all faces 
    gsSolidHalfFaceHandle addFace_4Vertices(gsSolidHeVertexHandle v0,
                                            gsSolidHeVertexHandle v1,
                                            gsSolidHeVertexHandle v2,
                                            gsSolidHeVertexHandle v3)
    {
        std::vector<gsSolidHeVertexHandle> V(4);
        V[0] = v0; V[1] = v1; V[2] = v2; V[3] = v3;
        return addFace(V);
    }
    
    /// Split a face \a f along a given spline. Returns the new face created
    /// as a result. The new face will have the edges from \a startVertex to
    /// \a endVertex and the reverse of \a domainSpline. The original face
    /// \a f will retain the edges from \a endVertex to \a startVertex and
    /// gain \a domainSpline (forwards).
    gsSolidHalfFaceHandle splitFace(gsSolidHalfFaceHandle f, 
                                    gsSolidHeVertexHandle startVertex, 
                                    gsSolidHeVertexHandle endVertex, 
                                    gsBSpline<T> *domainSpline);
    
    /// add a volume using its handle 
    void addVolume(gsVolumeHandle vol) 
        {
            volume.push_back(vol);
            vol->setId(numVolumes++);
        }
	
    /// add a volume using handles of its half faces
    void addVolume(std::vector<gsSolidHalfFaceHandle> hfaces) 
        {
	    gsVolumeHandle vol = this->makeVolume(hfaces);
            addVolume(vol);
        }
    
    /// Assigning mates for each HE	
    void setHeMate();    

    /// Define (TODO: detect automatically) nonconvex edges
    std::vector<gsSolidHalfEdgeHandle> detectNonConvexEdges(std::vector<int> const & ncEdgeV1, std::vector<int> const & ncEdgeV2);

    /// @brief Starting with a specified face, chase round all the faces
    /// connected to it and move them to a new volume.
    gsVolumeBlock<T> *newVolume(gsSolidHalfFaceHandle startingFace);

    /// sanity checks on the graph structure of the solid
    void checkStructure(bool checkVerts = false) const;

    /// @brief Make a face that can be used as the mate of the given face, add
    /// both faces and split off a new volume.
    gsSolidHalfFace<T> *addFaceWithMate(const std::vector<gsSolidHeVertexHandle> &Verts, gsTrimSurface<T> *surf);

    /// \brief Insert a new vertex to an edge of the volume.
    /// \param he       the half-edge into which to insert the vertex
    /// \note This insertion only affects the volume containing the HE.
    void insertNewVertex(gsSolidHalfEdgeHandle he);
    	
    /// @brief If there are impeding edges that make the solid inseparatable along a given HE \a he, this routine
    /// will create a new vertex on each impeding edges.
    void handleImpedingEdges(gsSolidHalfEdgeHandle he);

    //void addPatch(gsPatch<T> *f) { this->addFace(f) ;};
    //void addCurve(gsCurve<T> *f) { this->addEdge(f) ;};      

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Const members    
public:
    /// Check if the solid is separatable along a given edge \a he by collecting
    /// all "impeding edges" which connect the two faces incident to \a he.
    std::vector<gsSolidHalfEdgeHandle> impedingEdges(gsSolidHalfEdgeHandle he) const;

    /// plot edge graph	
    gsMultiPatch<T> plotEdgeGraph();
    
    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const;               
}; 
    
} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSolid.hpp)
#endif
