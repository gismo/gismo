/** @file gsSurfMesh.h

    @brief Half edge mesh structure

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsXml.h>

#include <gsMesh2/gsSurfMeshTopology.h>

namespace gismo
{

/** \brief A halfedge data structure for polygonal meshes with vertex positions.

    gsSurfMesh adds geometry to gsSurfMeshTopology: it stores a position of type
    \c Point per vertex and implements every operation that needs those
    positions (normals, edge lengths, subdivision, spline extraction, I/O).

    Everything purely combinatorial -- the Vertex / Halfedge / Edge / Face
    handles, the connectivity, the property containers, the iterators and
    circulators, and operations such as add_face(), collapse(), flip() or
    garbage_collection() -- is inherited unchanged from gsSurfMeshTopology and
    compiled only once into the library rather than once per \c Scalar.

    The nested types of the base are inherited, so \c gsSurfMesh<T>::Vertex,
    \c gsSurfMesh<T>::Halfedge, \c gsSurfMesh<T>::Vertex_property<X> and so on
    all keep working exactly as before.

    \tparam Scalar coefficient type of the vertex positions

    \sa gsSurfMeshTopology
*/
template <class Scalar>
class GISMO_EXPORT gsSurfMesh : public gsSurfMeshTopology
{
public:

    /// Non-templated topology base
    typedef gsSurfMeshTopology Base;

    using Self = gsSurfMesh<Scalar>;

    /// Point type
    typedef gsVector3d<Scalar> Point;

    /// Normal type
    typedef Point Normal;

    // the base declares overloads of these taking a Vertex; bring them in so
    // that the Point-taking overloads added below do not hide them
    using Base::split;
    using Base::insert_vertex;
    using Base::quad_split;

public:

    /// \name Construct, destruct, assignment
    //@{

    /// default constructor
    gsSurfMesh();

    /// Constructor from a set of points
    gsSurfMesh(const gsMatrix<Scalar> & pts);

    // destructor (is virtual, since we inherit from Geometry_representation)
    virtual ~gsSurfMesh();

    /// copy constructor: copies \c rhs to \c *this. performs a deep copy of all properties.
    gsSurfMesh(const gsSurfMesh& rhs) : Base() { operator=(rhs); }

    /// assign \c rhs to \c *this. performs a deep copy of all properties.
    gsSurfMesh& operator=(const gsSurfMesh& rhs);

    /// move \c rhs to \c *this
    gsSurfMesh& operator=(gsSurfMesh&& rhs) noexcept;

    /// assign \c rhs to \c *this. does not copy custom properties.
    gsSurfMesh& assign(const gsSurfMesh& rhs);

    //@}

public:

    /// read mesh from file \c filename. file extension determines file type.
    /// \sa write(const std::string& filename)
    bool read(const std::string& filename);

    /// write mesh to file \c filename. file extensions determines file type.
    /// \sa read(const std::string& filename)
    bool write(const std::string& filename) const;

public:

    /// add a new vertex with position \c p
    Vertex add_vertex(const Point& p);

    /// Add a mesh to the current one. The two meshes should be distinct from each other
    ///
    /// It retuns the vector of vertices that will be the
    /// map of local vertices'ids in the \c subMesh with the global vertices' id
    /// of out current mesh.
    ///
    /// \param subMesh Distinct mesh to be added in the current one
    gsVector<Vertex> add_mesh(gsSurfMesh& subMesh);

public:

    /** Split the face \c f by first adding point \c p to the mesh and then
     inserting edges between \c p and the vertices of \c f. For a triangle
     this is a standard one-to-three split.
     \sa split(Face, Vertex)
     */
    Vertex split(Face f, const Point& p) { Vertex v=add_vertex(p); split(f,v); return v; }

    /** Split the edge \c e by first adding point \c p to the mesh and then
     connecting it to the two vertices of the adjacent triangles that are
     opposite to edge \c e. Returns the halfedge pointing to \c p that is
     created by splitting the existing edge \c e.

     \attention This function is only valid for triangle meshes.
     \sa split(Edge, Vertex)
     */
    Halfedge split(Edge e, const Point& p) { return split(e, add_vertex(p)); }

    /** Subdivide the edge \c e = (v0,v1) by splitting it into the two edge
     (v0,p) and (p,v1). Note that this function does not introduce any
     other edge or faces. It simply splits the edge. Returns halfedge that
     points to \c p.
     \sa insert_vertex(Edge, Vertex)
     \sa insert_vertex(Halfedge, Vertex)
     */
    Halfedge insert_vertex(Edge e, const Point& p)
    {
        return insert_vertex(halfedge(e,0), add_vertex(p));
    }

    /// Quad-split uniformly at the half of each edge respectively.
    void quad_split();

    /**  Quad-split at uniform \c w positions on each edge respectively.

     Depending on the \c w, for each face, each each is splitted into \c w + 1 edges.
        * For w = 1: 1 Face ---> 4 Faces.
        * For w = 2: 1 Face ---> 9 Faces.

     \param w: Option regarding the edge split in each edge. For now works for \c w <= 2.
    */
    void quad_split(index_t w);

    /// Angle between halfedges
    ///
    ///
    /// \param h1: Halfedge 1
    /// \param h2: Halfedge 2
    real_t angle(Halfedge h1, Halfedge h2);

    /// Add property "v:halfedge" of normalized halfedges to all faces of the mesh.
    ///
    /// If the mesh is plotted, this property shows in a every face a halfedge.
    ///
    /// The property will be added to the vertices.
    /// close to vertices not on vertices, to be shown on face.
    void display_halfedge();

    /// Flips halfedge orientation i.e., if it is CW, becomes CCW.
    Self flip_orientation();

public:  // mesh operations related to subdivision schemes

    /// Augment mesh boundaries for boundary control on dual subdivision schemes.
    ///
    /// Necessary all faces that touch boundary be quads.
    ///
    /// In practice is the implementation of A. Nashri 1987 - "Polyhedral Subdivision Methods for Free-Form Surfaces"
    ///
    void polyhedral_modification_boundary();

    /// creates in-place the dual mesh (Dual-graph) for 2-manifolds without boundaries using barycentric method.
    void dual_mesh_inplace();

    /// returns dual mesh (Dual-graph) creation for 2-manifolds without boundaries using barycentric method.
    Self dual_mesh();

    /// calculate barycenter of a face
    Point face_barycenter(Face f);

public:

    /// position of a vertex (read only)
    const Point& position(Vertex v) const { return vpoint_[v]; }

    /// position of a vertex
    Point& position(Vertex v) { return vpoint_[v]; }

    /// compute the length of edge \c e.
    Scalar edge_length(Edge e) const;

    /// Vertex point coordinates
    const Vertex_property<Point> & points() const
    { return const_cast<Vertex_property<Point>&>(vpoint_); }

    /// Vertex point coordinates
    Vertex_property<Point> & points() { return vpoint_; }

    /// vector of vertex positions
    std::vector<Point>& pointsVec() { return vpoint_.vector(); }

    /// vector of vertex positions
    const std::vector<Point>& pointsVec() const { return vpoint_.vector(); }

    /// compute face normals by calling compute_face_normal(Face) for each face.
    void update_face_normals();

    /// compute normal vector of face \c f.
    Normal compute_face_normal(Face f) const;

    /// compute vertex normals by calling compute_vertex_normal(Vertex) for each vertex.
    void update_vertex_normals();

    /// compute normal vector of vertex \c v.
    Normal compute_vertex_normal(Vertex v) const;

public: // Operations related to b-splines and mesh

    /// Generate linear tensor-product patches (possibly merging faces)
    gsMultiPatch<Scalar> linear_patches() const;

    void mergeDoubleVertices();

public: // Extract spline functions

    /// Creates spline representations of mesh in a given degree
    gsMultiPatch<Scalar> asSpline(int deg) const;

    /// Creates a patch of spefic degree from a given halfedge in a mesh
    memory::unique_ptr<gsGeometry<Scalar> > asPatch(Halfedge h, int deg) const;

private:

    template<class S>
    friend std::ostream& operator<<(std::ostream& os, const gsSurfMesh<S> & sm)
    {
        os<<"gsSurfMesh with "<<sm.n_vertices()<<" vertices, "<<sm.n_edges()<<
            " edges and "<<sm.n_faces()<<" faces.\n";
        return os;
    }

private:

    Vertex_property<Point>   vpoint_;
    Vertex_property<Normal>  vnormal_;
    Face_property<Normal>    fnormal_;
};


namespace internal
{

template <class Scalar>
class GISMO_EXPORT gsXml< gsSurfMesh<Scalar> >
{
private:
    gsXml();
public:
    static std::string tag () { return "Mesh"; }
    static std::string type() { return ""; } //no inheritance

    GSXML_GET_POINTER(gsSurfMesh<Scalar>);
    static void get_into(gsXmlNode * node, gsSurfMesh<Scalar> & result);
    static gsXmlNode * put (const gsSurfMesh<Scalar> & obj, gsXmlTree & data);
};

}//namespace internal

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSurfMesh.hpp)
#endif
