/** @file gsFreeformSubdivision.h

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#pragma once

#include <gsMesh2/gsSurfMesh.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsMesh2/gsSubdivisionScheme.h>

namespace gismo
{

/// Data class that contains the additional information per face needed to perform freeform subdivision.
/// Represents a 5x5 grid of Bezier control points on the face.
/// The grid is supposed to be laid out as follows:
///
///    V0        E1         V1
///      +-----------------+
///      |  0  1  2  3  4  |
///      |                 |
///      |  5  6  7  8  9  |
///  E0  |                 |  E2
///      | 10 11 12 13 14  |
///      |                 |
///      | 15 16 17 18 19  |
///      |                 |
///      | 20 21 22 23 24  |
///      +-----------------+
///    V3         E3        V2
///
/// where C0-C3 are the vertices in the order mesh.vertices(face) traverses them and E0-E3 are the (half)edges in the order mesh.halfedges(face) traverses them.
/// In particular, the position of control points 0, 4, 20, 24 should correspond to the position of the vertices (if present).
class GISMO_EXPORT gsFreeformFaceData {

friend class gsFreeformSubdivision;

typedef gsSurfMesh::Point Point;
typedef gsSurfMesh::Vertex Vertex;
typedef gsSurfMesh::Face Face;
typedef gsSurfMesh::Halfedge Halfedge;
typedef gsSurfMesh::Edge Edge;

private: // consts that tell us about the orientation of the data
    // The control points for each (half)edge.
    const static std::vector<std::array<size_t, 5>> gs_FREEFORMDATA_EDGES;

    // The control points for each (half)edge, one layer farther in.
    const static std::vector<std::array<size_t, 3>> gs_FREEFORMDATA_INNEREDGES;

private: // members
    // The 25 bezier control points.
    std::array<gismo::gsVector3d<real_t>, 25> control_points;
    // A back reference to the face this data belongs to.
    Face face;

public: // Contructors
    /// Default contstructor with zeroed inner control points and empty face reference.
    gsFreeformFaceData()
    : control_points(),
      face(0) {}

    /// Default contstructor with zeroed inner control points.
    gsFreeformFaceData(Face face)
    : control_points(),
      face(face) {}

    /// Constructor that takes a mesh and chooses the control points of this patch to achieve C1 smoothness. 
    gsFreeformFaceData(const gsSurfMesh& mesh, gsSurfMesh::Face face);

public: // Control point accessors

    /// Returns an array containing the 5 points on the given edge of this face, ordered in the direction of the half edge.
    /// If the given half-edge does not belong to this face, this will fail in debug mode and return an array of zero vectors.
    std::array<gismo::gsVector3d<>, 5> edge_control_points(
      const gsSurfMesh& mesh,
      Halfedge hedge
    );

    /// Returns an array containing the 3 points on the given edge one layer in on this face, ordered in the direction of the half edge.
    /// If the given half-edge does not belong to this face, this will fail in debug mode and return an array of zero vectors.
    std::array<gismo::gsVector3d<>, 3> edge_inner_control_points(
      const gsSurfMesh& mesh,
      Halfedge hedge
    );

    /// Returns the corner point on the given vertex of this face.
    /// If the given vertex does not belong to this face, this will fail in debug mode and return an array of zero vectors.
    gismo::gsVector3d<> vertex_control_point(
      const gsSurfMesh& mesh,
      Vertex v
    );

    /// Returns the corner point on the given vertex of this face, one layer in.
    /// If the given vertex does not belong to this face, this will fail in debug mode and return an array of zero vectors.
    gismo::gsVector3d<> vertex_inner_control_point(
      const gsSurfMesh& mesh,
      Vertex v
    );

public: // Conversions
    /// Returns a Bezier patch corresponding to these control points.
    gismo::gsTensorBSpline<2, real_t> patch();
  
};//namespace internal

/// class for subdivision schemes in polygonal meshes.
class GISMO_EXPORT gsFreeformSubdivision : public gsSubdivisionScheme
{

public: // Constructors
  /// Default constructor.
  /// Catmull-Clark has no special options.
  gsFreeformSubdivision() : gsSubdivisionScheme() {}

public:
  void subdivide(gsSurfMesh &mesh) override;
  // gsSubdivisionMeshValidity valid_mesh(const gsSurfMesh &mesh) override;

  /// Takes a given mesh with free form data and makes it C1 by adjusting the outer 16 control points of each bezier patch.
  void make_c1(gsSurfMesh &mesh);
  
};//namespace internal

} // namespace gismo

