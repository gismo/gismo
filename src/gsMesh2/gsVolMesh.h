/** @file gsVolMesh.h

    @brief Half-face mesh structure for polyhedral volume meshes, with vertex
    positions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#include <gsMesh2/gsVolMeshTopology.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

/** \brief A half-face data structure for polyhedral volume meshes with vertex
    positions.

    gsVolMesh adds geometry to gsVolMeshTopology in exactly the way gsSurfMesh
    adds it to gsSurfMeshTopology: it stores a position of type \c Point per
    geometric Vertex and implements the operations that need those positions
    (edge lengths, barycenters, half-face normals, cell volumes, extraction of
    the boundary surface, I/O).

    Everything purely combinatorial -- the two tiers of handles, the
    \f$\beta_3\f$ relation, the iterators and circulators, add_cell() and
    garbage_collection() -- is inherited from gsVolMeshTopology and compiled
    only once into the library rather than once per \c Scalar.

    \tparam Scalar coefficient type of the vertex positions

    \sa gsVolMeshTopology, gsSurfMesh
*/
template <class Scalar>
class GISMO_EXPORT gsVolMesh : public gsVolMeshTopology
{
public:

    /// Non-templated topology base
    typedef gsVolMeshTopology Base;

    using Self = gsVolMesh<Scalar>;

    /// Point type
    typedef gsVector3d<Scalar> Point;

    /// Normal type
    typedef Point Normal;

    /// The surface mesh type that boundary_mesh() and cell_mesh() return
    typedef gsSurfMesh<Scalar> Surface;

    // the base declares add_vertex() without a position; bring it in so that
    // the Point-taking overload added below does not hide it
    using Base::add_vertex;

public:

    /// \name Construct, destruct, assignment
    //@{

    /// default constructor
    gsVolMesh();

    /// constructor from a set of points, one per column
    explicit gsVolMesh(const gsMatrix<Scalar> & pts);

    virtual ~gsVolMesh();

    /// copy constructor: copies \c rhs to \c *this. performs a deep copy of all properties.
    gsVolMesh(const gsVolMesh& rhs) : Base() { operator=(rhs); }

    /// assign \c rhs to \c *this. performs a deep copy of all properties.
    gsVolMesh& operator=(const gsVolMesh& rhs);

    /// move \c rhs to \c *this
    gsVolMesh& operator=(gsVolMesh&& rhs) noexcept;

    /// assign \c rhs to \c *this. does not copy custom properties.
    gsVolMesh& assign(const gsVolMesh& rhs);

    //@}

public:

    /// add a new vertex with position \c p
    Vertex add_vertex(const Point& p);

public:

    /// \name Geometry
    //@{

    /// position of a vertex (read only)
    const Point& position(Vertex v) const { return vpoint_[v]; }

    /// position of a vertex
    Point& position(Vertex v) { return vpoint_[v]; }

    /// vertex point coordinates
    const Vertex_property<Point> & points() const
    { return const_cast<Vertex_property<Point>&>(vpoint_); }

    /// vertex point coordinates
    Vertex_property<Point> & points() { return vpoint_; }

    /// vector of vertex positions
    std::vector<Point>& pointsVec() { return vpoint_.vector(); }

    /// vector of vertex positions
    const std::vector<Point>& pointsVec() const { return vpoint_.vector(); }

    /// the length of edge \c e
    Scalar edge_length(Edge e) const;

    /// the barycenter of the vertices of half-face \c hf
    Point barycenter(Halfface hf) const;

    /// the barycenter of the vertices of face \c f
    Point barycenter(Face f) const { return barycenter(halfface(f,0)); }

    /// the barycenter of the vertices of cell \c c
    Point barycenter(Cell c) const;

    /** the area-weighted normal of half-face \c hf, computed with Newell's
        method, so that it is well defined for non-planar polygons.

        Half-faces are oriented outwards with respect to their cell, so the
        normal points out of cell(hf).
    */
    Normal normal(Halfface hf) const;

    /** the signed volume of cell \c c.

        Computed by the divergence theorem over the boundary of the cell, with
        each half-face triangulated as a fan.  The value is positive for a cell
        whose half-faces are consistently oriented outwards, which is what
        add_cell() produces, so a negative volume signals an inverted cell.
    */
    Scalar volume(Cell c) const;

    /// the sum of volume(Cell) over all cells
    Scalar volume() const;

    //@}

public:

    /// \name Surface views
    //@{

    /** the boundary of the mesh, as a standalone gsSurfMesh.

        The faces are oriented outwards.  The returned mesh carries a vertex
        property \c "v:volvertex" of type \c index_t holding, for each of its
        vertices, the index of the gsVolMesh vertex it came from.
    */
    Surface boundary_mesh() const;

    /// the boundary of the single cell \c c, as a standalone gsSurfMesh,
    /// carrying the same \c "v:volvertex" back-map as boundary_mesh()
    Surface cell_mesh(Cell c) const;

    //@}

public:

    /// \name Input / output
    //@{

    /** read the mesh from \c filename; the extension determines the file type.

        Supported: \c .msh (Gmsh, ASCII versions 2.2 and 4.1), \c .vtk (legacy
        VTK, ASCII unstructured grid) and \c .vtu (VTK XML, inline ASCII).
        On failure the mesh is left untouched.
        \sa write(const std::string&)
    */
    bool read(const std::string& filename);

    /** write the mesh to \c filename; the extension determines the file type.

        Supported: \c .msh (Gmsh ASCII 2.2), \c .vtk (legacy VTK) and \c .vtu
        (VTK XML) for the volume mesh itself, and any extension gsSurfMesh can
        write (\c .off, \c .obj, \c .stl) for its boundary surface.

        Only \c .vtu can express a general polyhedron; the other two formats
        report the cells they had to leave out.
        \sa read(const std::string&)
    */
    bool write(const std::string& filename) const;

    //@}

private:

    template<class S>
    friend std::ostream& operator<<(std::ostream& os, const gsVolMesh<S> & vm)
    {
        os<<"gsVolMesh with "<<vm.n_vertices()<<" vertices, "<<vm.n_edges()<<
            " edges, "<<vm.n_faces()<<" faces and "<<vm.n_cells()<<" cells.\n";
        return os;
    }

private:

    Vertex_property<Point> vpoint_;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsVolMesh.hpp)
#endif
