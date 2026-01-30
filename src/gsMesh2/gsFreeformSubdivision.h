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
#include <vector>

namespace gismo
{

/// Data class that contains the additional information per face needed to perform freeform subdivision.
/// Represents a 5x5 grid of Bezier control points on the face.
/// The grid is supposed to be laid out as follows:
///
///    V0        E1         V1
///      +-----------------+
///      | 00 01 02 03 04  |
///      |                 |
///      | 10 11 12 13 14  |
///  E0  |                 |  E2
///      | 20 21 22 23 24  |
///      |                 |
///      | 30 31 32 33 34  |
///      |                 |
///      | 40 41 42 43 44  |
///      +-----------------+
///    V3         E3        V2
///
/// where C0-C3 are the vertices in the order mesh.vertices(face) traverses them and E0-E3 are the (half)edges in the order mesh.halfedges(face) traverses them.
/// In particular, the position of control points 0, 4, 20, 24 should correspond to the position of the vertices (if present).
class GISMO_EXPORT gsFreeformFaceData {

friend class gsFreeformSubdivision;

using Point = gsSurfMesh::Point;
using Vertex = gsSurfMesh::Vertex;
using Face = gsSurfMesh::Face;
using Halfedge = gsSurfMesh::Halfedge;
using Edge = gsSurfMesh::Edge;

private: // members
    // The 25 bezier control points.
    gsMatrix<gismo::gsVector3d<real_t>, 5, 5> control_points;
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

    /// Returns a vector containing the control points along to a halfedge, `inset` rows of control points offset to the middle, and always in the same direction as the halfedge.
    /// I.e. with inset 0, returns the entire row or column on this halfedge.
    /// I.e. with inset 1, one row/column farther in and without the first and last element.
    /// I.e. with inset 2, two rows/columns farther in and without the first, second, second-to-last and last element.
    /// If the given half edge does not belong to this face, this will fail in debug mode and return an empty vector.
    /// If the inset is to large (greater than half the size of the control net), this will return an empty vector.
    /// Pointers returned by this function are valid as long as the underlying mesh remains unchanged and should be regenerated after.
    std::vector<gsVector3d<>*> edge_control_points(
      gsSurfMesh& mesh,
      Halfedge hedge,
      size_t inset
    );

    /// Returns a the control point closest to the given vertex, in the `inset`th outer most layer of the control net.
    /// I.e. with inset 0, return the correct corner of the control net associated to this vertex.
    /// I.e. with inset 1, one row & column farther in.
    /// If the given vertex does not belong to this face, this will fail in debug mode or return a nullpointer.
    /// If the inset is to large (greater than half the size of the control net), this will return an empty vector.
    /// Pointers returned by this function are valid as long as the underlying mesh remains unchanged and should be regenerated after.
    gsVector3d<>* vertex_control_point(
      gsSurfMesh& mesh,
      Vertex v,
      size_t inset
    ) ;

public: // Conversions
    /// Returns a Bezier patch corresponding to these control points.
    const gismo::gsTensorBSpline<2, real_t> patch() const;
  
};//namespace internal

/// class for subdivision schemes in polygonal meshes.
class GISMO_EXPORT gsFreeformSubdivision : public gsSubdivisionScheme
{

public: // Constructors
  /// Default constructor.
  /// Catmull-Clark has no special options.
  gsFreeformSubdivision() : gsSubdivisionScheme() {}

private:
  /// Helper function
  /// Splits the given matrix, interpreted as a control net of a Bezier patch, into two control nets of Bezier patches of equal degree that in their union again form the original net.
  /// This is done using the algorithm of deCasteljau.
  /// The net is assumed to be quadratic in size.
  /// The net is split horizontally, in the 'row' direction of the matrix.
  static std::array<gsMatrix<gsVector3d<>, Dynamic, Dynamic>, 2> deCasteljau(const gsMatrix<gsVector3d<>, Dynamic, Dynamic>& control_net);

  /// Rotates a matrix around its center, returning the rotated matrix.
  static gsMatrix<gsVector3d<>, Dynamic, Dynamic> rotate(const gsMatrix<gsVector3d<>, Dynamic, Dynamic>& mat);
public:
  void subdivide(gsSurfMesh &mesh) override;
  // gsSubdivisionMeshValidity valid_mesh(const gsSurfMesh &mesh) override;

  /// Takes a given mesh with free form data and makes it C1 by adjusting the outer 16 control points of each bezier patch.
  void make_c1(gsSurfMesh &mesh);
  
};//namespace internal


} // namespace gismo

