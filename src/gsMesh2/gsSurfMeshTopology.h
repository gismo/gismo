/** @file gsSurfMeshTopology.h

    @brief Half edge mesh topology: the scalar-independent (non-templated) base
    of gsSurfMesh, holding connectivity, handles, iterators and properties.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsMeshTopology.h>

namespace gismo
{

/** \brief A halfedge data structure for polygonal meshes: topology only.

    This class carries everything about a polygonal \em surface mesh that does
    not depend on the coordinate type.  The handles, the connectivity, the
    property containers, the iterators and circulators and the incremental
    construction come from gsMeshTopology; what is added here is the reading of
    those entities as the cells of a surface, together with the Euler operators
    that only make sense under that reading: triangulate(), collapse(), flip(),
    split() and their helpers.

    Vertex positions and every operation needing them live in the derived
    gsSurfMesh<Scalar>.

    Being non-templated, all of this compiles exactly once into the library
    instead of once per scalar type.

    \sa gsMeshTopology, gsSurfMesh, gsVolMeshTopology
*/
class GISMO_EXPORT gsSurfMeshTopology : public gsMeshTopology
{
public:

    /// Dimension-agnostic topology base
    typedef gsMeshTopology Base;

    using Self = gsSurfMeshTopology;

public:

    /// \name Construct, destruct, assignment
    ///
    /// The mesh holds no state beyond the base, so the compiler-generated
    /// members are correct: they forward to gsMeshTopology, which performs the
    /// deep copy of the property containers and re-binds the property handles.
    //@{

    gsSurfMeshTopology() : Base() {}

    virtual ~gsSurfMeshTopology() {}

    //@}

public:

    /// returns whether the mesh a triangle mesh. this function simply tests
    /// each face, and therefore is not very efficient.
    /// \sa trianglate(), triangulate(Face)
    bool is_triangle_mesh() const;

    /// returns whether the mesh a quad mesh. this function simply tests
    /// each face, and therefore is not very efficient.
    bool is_quad_mesh() const;

    /// triangulate the entire mesh, by calling triangulate(Face) for each face.
    /// \sa trianglate(Face)
    void triangulate();

    /// triangulate the face \c f
    /// \sa trianglate()
    void triangulate(Face f);

    /// returns whether collapsing the halfedge \c h is topologically legal.
    /// \attention This function is only valid for triangle meshes.
    bool is_collapse_ok(Halfedge h);

    /** Collapse the halfedge \c h by moving its start vertex into its target
     vertex. For non-boundary halfedges this function removes one vertex, three
     edges, and two faces. For boundary halfedges it removes one vertex, two
     edges and one face.
     \attention This function is only valid for triangle meshes.
     \attention Halfedge collapses might lead to invalid faces. Call
     is_collapse_ok(Halfedge) to be sure the collapse is legal.
     \attention The removed items are only marked as deleted. You have
     to call garbage_collection() to finally remove them.
     */
    void collapse(Halfedge h);


    /** Split the face \c f by inserting edges between \c p and the vertices
     of \c f. For a triangle this is a standard one-to-three split.
     \sa split(Face, const Point&)
     */
    void split(Face f, Vertex v);



    /** \brief Triangle - split a face connecting vertex \a v with the vertices in vector \c edgeverts

        The vector should contain a vertex for each edge of a face.
        The halfedge orientation in each trinagle is the same as in the original face.

        \param edgeverts: Vector with vertices. We need one vertex per edge.
        \param f: Face to be splitted
        \param v: Inserted vertex for split.
    */
    void split_to_triangles(std::vector<Vertex>& edgeverts, Face f, Vertex v);


    /// Quad-split face connecting vertex \a v, starting from corner
    /// \a s of the face
    /// \a f is assumed to have 8 vertices, and contains halfedge \a s
    void quad_split(Face f, Vertex v, Halfedge s);

    /** Split the edge \c e by connecting vertex \c v it to the two
     vertices of the adjacent triangles that are opposite to edge \c
     e. Returns the halfedge pointing to \c p that is created by splitting
     the existing edge \c e.

     \attention This function is only valid for triangle meshes.
     \sa split(Edge, Point)
     */
    Halfedge split(Edge e, Vertex v);


    /** Subdivide the edge \c e = (v0,v1) by splitting it into the two edge
     (v0,v) and (v,v1). Note that this function does not introduce any other
     edge or faces. It simply splits the edge. Returns halfedge that points to \c p.
     \sa insert_vertex(Edge, Point)
     \sa insert_vertex(Halfedge, Vertex)
     */
    Halfedge insert_vertex(Edge e, Vertex v)
    {
        return insert_vertex(halfedge(e,0), v);
    }

    /** Subdivide the edge \c e = (v0,v1) by splitting it into the two edge
     (v0,v) and (v,v1). Note that this function does not introduce any other
     edge or faces. It simply splits the edge. Returns halfedge that points to \c p.
     \sa insert_vertex(Edge, Point)
     \sa insert_vertex(Edge, Vertex)
     */
    Halfedge insert_vertex(Halfedge h, Vertex v);


    /// insert edge between the to-vertices v0 of h0 and v1 of h1.
    /// returns the new halfedge from v0 to v1.
    /// \attention h0 and h1 have to belong to the same face
    Halfedge insert_edge(Halfedge h0, Halfedge h1);


    /** Check whether flipping edge \c e is topologically
     \attention This function is only valid for triangle meshes.
     \sa flip(Edge)
     */
    bool is_flip_ok(Edge e) const;

    /** Flip edge \c e: Remove edge \c e and add an edge between the two vertices
     opposite to edge \c e of the two incident triangles.
     \attention This function is only valid for triangle meshes.
     \sa is_flip_ok(Edge)
     */
    void flip(Edge e);

    /// Mesh statistics
    ///
    /// Print info on mesh as:
    ///    * Number of vertices, faces and edges.
    ///    * Number of extraordinary vertices, with minimum and maximum valence for interior and boundary.
    ///    * Number of extraordinary faces, with minimum and maximum valence for interior and boundary.
    ///
    /// \param eoc_verbose: If true, returns the number of each extraordinary case (EV,EF)
    void mesh_statistics(bool eoc_verbose = false);

public:

    unsigned int hcount(Vertex v, const Halfedge_property<bool>  & prop) const;

    /// returns the sum of all valences of faces in the mesh
    unsigned int face_valence_sum() const;

    inline unsigned int num_faces(Self::Vertex v) const
    {
        unsigned int cnt = 0;
        for (auto f: faces(v)) { ++cnt; GISMO_UNUSED(f); }
        return cnt;
    }

    /// Returns the degree of a vertex.
    /*
     Let k(v) = #incident faces
     Let n(v) = #incident edges
     Degree d(v) := if (v interior) d(v)= k=n else d(v)= 2*k=2*(n-1)
     Therefore:
     - Corner vertex:  d==2
     - Regular vertex (boundary or interior): d==4
     - Irregular vertex (boundary or interior): d!=4 and d!=2
     Examples:
     - Corner vertex            : k=1, n=2,   d=2
     - Regular interior vertex  : k=4, n=4,   d=4
     - Regular boundary vertex  : k=2, n=3,   d=4
     - Irregular interior vertex: k=k, n=k,   d=k   (i.e. 3,5,6,7,..)
     - Irregular boundary vertex: k=k, n=k+1, d=2*k (i.e. 6,8,10,12,..)
    */
    inline unsigned int vertex_degree(Self::Vertex v) const
    {
        const unsigned int k = num_faces(v);
        return is_boundary(v) ? 2*k : k;
    }

    /// returns whether \c v is an extraordinary vertex, i.e. neither a corner
    /// (degree 2) nor a regular vertex (degree 4)
    inline bool is_regular(Vertex v) const
    {
        int d = vertex_degree(v);
        return d!=4 && d!=2;
    }

    // Returns true if there is a halfedge with hflag set to true emenating from vertex \a v
    bool has_flag(Vertex v, const Halfedge_property<bool> & hflag) const;

protected: //------------------------------------------------ helper functions

    /// Helper for halfedge collapse
    void remove_edge(Halfedge h);

    /// Helper for halfedge collapse
    void remove_loop(Halfedge h);

    friend std::ostream& operator<<(std::ostream& os, const gsSurfMeshTopology & sm)
    {
        os<<"gsSurfMeshTopology with "<<sm.n_vertices()<<" vertices, "<<sm.n_edges()<<
            " edges and "<<sm.n_faces()<<" faces.\n";
        return os;
    }
};

} // namespace gismo
