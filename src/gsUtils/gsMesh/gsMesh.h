/** @file gsMesh.h

    @brief Provides declaration of the Mesh class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mayer
*/

#pragma once

#include <gsUtils/gsMesh/gsMeshElement.h>
#include <gsUtils/gsMesh/gsVertex.h>
#include <gsUtils/gsMesh/gsFace.h>
#include <gsUtils/gsMesh/gsEdge.h>
#include <gsUtils/gsSortedVector.h>


namespace gismo {

/**
   \brief Class Representing a triangle mesh with 3D vertices.

   \ingroup Utils
*/
template <class T>
class GISMO_EXPORT gsMesh : public gsMeshElement<T>
{
public:
    typedef memory::shared_ptr<gsMesh> Ptr;
    typedef memory::unique_ptr<gsMesh> uPtr;
    typedef gsMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsVertexHandle VertexHandle;
    typedef typename MeshElement::gsFaceHandle FaceHandle;
    typedef typename MeshElement::gsEdgeHandle EdgeHandle;
    //typedef typename gsMeshElement<Vertex>::gsCellHandle gsCellHandle;
    typedef gsEdge<T>                         Edge;
    typedef gsVertex<T>                       Vertex;

public:

    gsMesh() : MeshElement()
    { }

    gsMesh(const gsMesh<T> & mesh) : MeshElement()
    {
        this->operator=(mesh);
    }

    gsMesh& operator=(const gsMesh& other)
    {
        if (this!=&other)
        {
            // Assert a already cleared gsMesh (getId)
            cloneAll(other.m_vertex, m_vertex);  // fine, new pointers for all vertices

            // copy all pointers to it's original counterpart in face
            cloneAll(other.m_face, m_face);
            for (size_t i = 0; i < other.m_face.size(); ++i)
            {
                for (size_t j = 0; j != 3; ++j)
                {
                    GISMO_ASSERT(m_vertex[other.m_face[i]->vertices[j]->getId()]->getId() == other.m_face[i]->vertices[j]->getId(), "gsMesh(const gsMesh<T> & mesh): getId() of vertex and face don't match");
                    m_face[i]->vertices[j] = m_vertex[other.m_face[i]->vertices[j]->getId()];
                }
            }

            // iterate over all edges and make them new
            m_edge = other.m_edge;
            for (size_t i = 0; i != other.m_edge.size(); ++i)
            {
                GISMO_ASSERT(other.m_edge[i].source->getId() == m_vertex[i]->getId(), "gsMesh(const gsMesh<T> & mesh): getId() of vertex and edge.source don't match");
                m_edge[i].source = m_vertex[other.m_edge[i].source->getId()];

                GISMO_ASSERT(other.m_edge[i].source->getId() == m_vertex[i]->getId(), "gsMesh(const gsMesh<T> & mesh): getId() of vertex and edge.target don't match");
                m_edge[i].target = m_vertex[other.m_edge[i].target->getId()];
            }
        }
        return *this;
    }

    /// Create a mesh corresponding to the parametric domain of the
    /// basis \a basis.
    /// @param basis Basis whose parametric domain should be meshed.
    /// @param midPts Number of intermediate vertices to be inserted
    /// to each element side. Cf. gsMesh::addLine(VertexHandle,
    /// VertexHandle, int).
    explicit gsMesh(const gsBasis<T> & basis, int midPts = 0);

    virtual ~gsMesh();

    // void addVertex(gsVertexHandle v){
    //     vertices.push_back(v);
    //     if(!initialized) {
    //         bb.pMax.x() = bb.pMin.x() = v->x();
    //         bb.pMax.y() = bb.pMin.y() = v->y();
    //         bb.pMax.z() = bb.pMin.z() = v->z();
    //         initialized = true;
    //     }
    //     else {
    //         if (bb.pMax.x() < v->x())
    //             bb.pMax.x() = v->x();
    //         if (bb.pMin.x() > v->x())
    //             bb.pMin.x() = v->x();

    //         if (bb.pMax.y < v->y())
    //             bb.pMax.y = v->y();
    //         if (bb.pMin.y > v->y())
    //             bb.pMin.y = v->y();

    //         if (bb.pMax.z() < v->z())
    //             bb.pMax.z() = v->z();
    //         if (bb.pMin.z() > v->z())
    //             bb.pMin.z() = v->z();
    //     }
    //     numVertices++;
    // };

    VertexHandle addVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z=0);

    VertexHandle addVertex(gsVector<T> const & u);

    void addEdge(VertexHandle v0, VertexHandle v1);

    void addEdge(int const vind0, int const vind1);

    void addEdge(gsVector<T> const & u0,
                 gsVector<T> const & u1 );

    FaceHandle addFace(std::vector<VertexHandle> const & vert);

    FaceHandle addFace(VertexHandle const & v0, VertexHandle const & v1,
                       VertexHandle const & v2);

    FaceHandle addFace(VertexHandle const & v0, VertexHandle const & v1,
                       VertexHandle const & v2,  VertexHandle const & v3);

    FaceHandle addFace(std::vector<int> const & vert);

    FaceHandle addFace(const int v0, const int v1, const int v2);

    FaceHandle addFace(const int v0, const int v1, const int v2, const int v3);

    /// Add to the mesh a list of vertices connected with edges
    /// One edges are added between two successive vertices in the input.
    /// Vertices of the edges are assumed to be inexistent vertices.
    /// \param points matrix containing in the columns a list of points
    // works for 2D vertices for now
    void addLine(gsMatrix<T> const & points);

    /// Inserts a straight line in the mesh, between \a v0 and \a v1,
    /// with \a midPts intermediate vertices distributed linearly
    /// between \a v0 and \a v1 (used for plotting).
    void addLine(VertexHandle v0, VertexHandle v1, int midPts = 0);

    std::ostream &print(std::ostream &os) const;

//    bool getTurnDirection(gsVector3d<T> A, gsVector3d<T> B, gsVector3d<T> C)
//    {
//        gsVector3d<T> vec1 = B - A ;
//        gsVector3d<T> vec2 = C - B ;
//        gsVector3d<T> normal = vec2.cross( vec1 ) ;
//        if (conditionedAngle( vec1,  vec2,  normal) >= EIGEN_PI)
//        {
//            return 1;
//        }
//        else
//        {
//            return 0;
//        }
//    }

    /** \brief reorders the vertices of all faces of an .stl mesh,
     * such that only 1 vertex is used instead of #(adjacent
     * triangles) vertices
     */
    gsMesh& cleanMesh();

    gsMesh& reserve(size_t vertex, size_t face, size_t edge);

    size_t numVertices() const { return m_vertex.size(); }
    size_t numEdges()    const { return m_edge.size(); }
    size_t numFaces()    const { return m_face.size(); }

    const std::vector<VertexHandle > & vertices() const
    { return m_vertex; }

    const std::vector<Edge > & edges() const
    { return m_edge; }

    const std::vector<FaceHandle > & faces() const
    { return m_face; }

    const Vertex & vertex(size_t i) const { return *m_vertex[i]; }
    Vertex & vertex(size_t i) { return *m_vertex[i]; }

public: //protected: -- todo

    std::vector<VertexHandle > m_vertex;
    std::vector<FaceHandle >  m_face;
    gsSortedVector<Edge> m_edge;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMesh.hpp)
#endif
