/** @file gsVolCatmullClark.h

    @brief Catmull-Clark subdivision on a volumetric mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsVolSubdivisionScheme.h>

namespace gismo
{

/** \brief Catmull-Clark subdivision for polyhedral volume meshes.

    The solid analogue of gsCatmullClark, after MacCracken and Joy,
    <em>Free-form deformations with lattices of arbitrary topology</em>,
    SIGGRAPH 1996.

    <b>Topology.</b> Every (cell, corner) pair becomes one new cell. Walking the
    darts around a corner inside its cell gives the cyclic sequence of edges
    \f$e_i\f$ and faces \f$f_i\f$ at that corner, with \f$f_i\f$ bounded by
    \f$e_i\f$ and \f$e_{i+1}\f$. The new cell is bounded by

    - \f$k\f$ quads around the vertex, \f$(V, E_i, F_i, E_{i+1})\f$, and
    - \f$k\f$ quads around the cell point, \f$(E_i, F_{i-1}, C, F_i)\f$,

    where \f$V, E, F, C\f$ are the vertex, edge, face and cell points. That is
    \f$2k\f$ faces over \f$2k+2\f$ vertices, which for \f$k=3\f$ is exactly a
    hexahedron.

    The refined mesh is therefore all-hexahedral precisely when every corner of
    every cell has three incident faces -- true for tetrahedra, hexahedra and
    prisms, and false at the apex of a pyramid, which has four and yields an
    8-faced cell.

    A corner of valence \f$k\f$ produces a cell in which both the vertex point
    and the cell point again have valence \f$k\f$, so a non-trivalent corner is
    \em not repaired by further refinement: it persists as two poles of the same
    valence. A mesh of tetrahedra, hexahedra and prisms becomes hexahedral after
    one step; one containing pyramids never does.

    <b>Geometry.</b> The interior masks are MacCracken and Joy's:

    - cell point: the centroid of the cell's vertices;
    - face point: the average of the face's vertices and of the cell points of
      the cells sharing it;
    - edge point: the average of the edge's two endpoints, the face points of
      its faces and the cell points of its cells;
    - vertex point: \f$(A + 3B + 3C + V)/8\f$, with \f$A\f$ the average of the
      incident cell points, \f$B\f$ of the incident face points and \f$C\f$ of
      the incident edge midpoints.

    \attention Boundary entities currently use these same interior masks. The
    published scheme instead applies the ordinary surface Catmull-Clark masks
    there, so that the boundary converges to the Catmull-Clark limit surface.
    Until that is implemented <b>the boundary of the mesh shrinks</b> with every
    step and does not agree with gsCatmullClark applied to boundary_mesh().
    The interior is unaffected.

    Unlike the surface schemes this one does not edit the mesh in place: a
    3-map cannot be refined by local surgery, so the refined mesh is built from
    scratch and moved over the target.

    \tparam Scalar coefficient type of the vertex positions

    \sa gsCatmullClark, gsVolSubdivisionScheme, gsVolMesh
*/
template <class Scalar=real_t>
class GISMO_EXPORT gsVolCatmullClark : public gsVolSubdivisionScheme<Scalar>
{
public:
    typedef gsVolSubdivisionScheme<Scalar> Base;

    typedef typename gsVolMesh<Scalar>::Point    Point;
    typedef typename gsVolMesh<Scalar>::Vertex   Vertex;
    typedef typename gsVolMesh<Scalar>::Edge     Edge;
    typedef typename gsVolMesh<Scalar>::Face     Face;
    typedef typename gsVolMesh<Scalar>::Cell     Cell;
    typedef typename gsVolMesh<Scalar>::Corner   Corner;
    typedef typename gsVolMesh<Scalar>::Halfedge Halfedge;
    typedef typename gsVolMesh<Scalar>::Halfface Halfface;

public: // Constructors

    /// \brief Constructor with a mesh to target.
    explicit gsVolCatmullClark(gsVolMesh<Scalar>* mesh = nullptr)
    : gsVolSubdivisionScheme<Scalar>()
    {
        this->assign(mesh);
    }

    /// Apply one subdivision step to \a mesh
    static void apply(gsVolMesh<Scalar>& mesh, size_t repetitions = 1)
    {
        gsVolCatmullClark<Scalar> cc(&mesh);
        cc.subdivide(repetitions);
    }

    // declaring check_mesh() below would otherwise hide the base overload
    // that takes a mesh and targets it
    using Base::check_mesh;

    /// \brief The scheme needs a mesh with at least one cell; there is no
    /// further restriction, cells of any polyhedral shape are refined.
    typename Base::gsSubdivisionMeshValidity check_mesh() GISMO_OVERRIDE;

protected:

    void subdivide_impl() GISMO_OVERRIDE;

}; // class gsVolCatmullClark

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsVolCatmullClark.hpp)
#endif
