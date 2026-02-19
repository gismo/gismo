/** @file gsFreeformSubdivision.h

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#pragma once

#include <gsMesh2/gsSubdivisionScheme.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

/// \brief Subdivision algorithm based on freeform splines.
///
/// Class for subdivision schemes based on freeform spline control nets on
/// quadrangular meshes. Also provides other support functions for working with
/// such meshes.
template <size_t N, size_t D>
class GISMO_EXPORT gsFreeformSubdivision : public gsSubdivisionScheme
{

public: // Constructors
    /// \brief Default constructor.
    ///
    /// Default constructor. Sets no options and leaves the targeted mesh as a
    /// nullpointer.
    gsFreeformSubdivision() : gsSubdivisionScheme() {}

    /// \brief Constructor with a mesh to target.
    ///
    /// Constructor that accepts a mesh to be targeted by this constructor.
    /// Sets no options.
    gsFreeformSubdivision(gsSurfMesh* mesh) : gsSubdivisionScheme(mesh) {}

private: // Helper functions
    /// \brief Performs deCasteljau algorithm to split a matrix into two.
    ///
    /// Splits the given matrix, interpreted as a control net of a Bezier patch,
    /// into two control nets of Bezier patches of equal degree that in their
    /// union again form the original net. This is done using the algorithm of
    /// deCasteljau. The net is assumed to be quadratic in size. The net is
    /// split horizontally, in the 'row' direction of the matrix.
    ///
    /// \param control_net The original matrix to be split.
    static std::array<gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>, 2>
    deCasteljau(
        const gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>& control_net);

    /// \brief Loads a model patch of the given valence and type.
    ///
    /// Looks in the `filedata/freeformsubdivision` folder for model patches for
    /// the EV subdivision and loads them to a gismo Bezier patch. The patches
    /// will combine in such a way that the patches `fine_1` to `fine_4` form a
    /// perfect partition of patch `coarse` as follows:
    /// ```
    /// fine_4   fine_3
    /// fine_1   fine_2
    /// ```
    /// With the u/v-Direction of `coarse` pointing left/up from the bottom left
    /// and the u direction of each `fine_v` pointing outward from the central
    /// shared point and having the same orientation (e.g. `fine_1` points u
    /// leftwards and v downwards).
    ///
    /// \param valence The valence of the desired patch. Valid values are 3, 5,
    /// 6, 7, 8, 9, 10.
    /// \param v The subtype of the patch. Valid values are "coarse", "fine_1",
    /// "fine_2", "fine_3", "fine_4".
    static gismo::gsTensorBSpline<2, real_t>
    load_model_patch(int valence, std::string subtype);

    /// \brief Re-orients the faces of the given mesh.
    ///
    /// This method changes the assigned halfedge of each face in the given mesh
    /// in such a way that, if the face has an extraordinary vertex, the
    /// halfedge of that face is the one pointing to that vertex. This will
    /// cause the first vertex of that face to be the EV. If any faces in the
    /// mesh have multiple EVs, there is no guarantee on which one will be
    /// chosen and it is suggested you change your mesh so that EVs have a
    /// larger distance between them.
    ///
    /// \param mesh The mesh to be re-oriented.
    void orient_faces();

    /// \brief Orders a vector of 4 faces such that the face with a given vertex
    /// is first.
    ///
    /// This method accepts a vector of four elements and a vertex. It will
    /// search each of the faces for the given vertex and then rotate the 4
    /// elements until a face containing that vertex is first.
    ///
    /// \param first_vertex The target vertex that determines the face order.
    /// Should be contained in at least one of the faces.
    /// \param faces The faces to be ordered. Should have exactly 4 elements or
    /// unexpected results may occur.
    std::array<Face, 4> order_faces(Vertex first_vertex,
                                    std::array<Face, 4> faces);

public:
    gsSubdivisionScheme::gsSubdivisionMeshValidity check_mesh() override;

    void subdivide() override;

    /// \brief Initializes C0 control nets on all faces.
    ///
    /// Takes a given mesh without freeform data and initializes a freeform data
    /// control net with C0 continuity on each face.
    void initialize_data();

    /// \brief Turns a $C^0$ set of control nets into a $C^s$ set.
    ///
    /// Takes a given mesh with freeform data and makes it $C^s$ by adjusting
    /// the outer layer of control points of each bezier patch. This only causes
    /// $C^s$ at edges and ordinary vertices. No guarantee is made for
    /// extraordinary vertices.
    ///
    /// \param degree The degree of smoothness desired. As of now, only $C^1$ is
    /// supported.
    void smooth(size_t degree);

    /// \brief Converts to a Gismo multipatch object.
    ///
    /// Converts a given mesh with freeform data into a multipatch that can be
    /// easily displayed by e.g. Paraview. Each face and its control net are
    /// converted to one appropriately sized patch.
    gsMultiPatch<> multipatch();

}; // namespace internal

/// \brief Manages a control net on a face of a surf mesh.
///
/// Data class that contains the additional information per face needed to
/// perform freeform subdivision. Represents a 5x5 grid of Bezier control points
/// on the face. The grid is supposed to be laid out as follows:
///
/// ```
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
/// ```
///
/// where `C0`-`C3` are the vertices in the order `mesh.vertices(face)`
/// traverses them and `E0`-`E3` are the (half)edges in the order
/// `mesh.halfedges(face)` traverses them. In particular, the position of
/// control points `00`, `04`, `40`, `44` should correspond to the position of
/// the vertices (if present).
template <size_t N, size_t D> class GISMO_EXPORT gsFreeformFaceData
{
    template <size_t M, size_t E> friend class gsFreeformSubdivision;

    using Point = gsSurfMesh::Point;
    using Vertex = gsSurfMesh::Vertex;
    using Face = gsSurfMesh::Face;
    using Halfedge = gsSurfMesh::Halfedge;
    using Edge = gsSurfMesh::Edge;

private: // members
    // The 25 bezier control points.
    gsMatrix<gismo::gsVector<real_t, D>, N, N> control_points;
    // A back reference to the face this data belongs to.
    Face face;

private: // helpers
public:  // Contructors
    /// \brief Default constructor with everything empty.
    ///
    /// Default constructor. All control points are the zero vector, and the
    /// back reference to the face is empty.
    gsFreeformFaceData() : control_points(), face(0)
    {
        for (size_t i = 0; i < N * N; ++i)
        {
            control_points(i) = gsVector<real_t, D>::Zero(D);
        }
    }

    /// \brief Zero constructor with a back reference.
    ///
    /// Basic constructor. All control points are the zero vector, and the
    /// back reference points to the given face.
    gsFreeformFaceData(Face face) : control_points(), face(face)
    {
        for (size_t i = 0; i < N * N; ++i)
        {
            control_points(i) = gsVector<real_t, D>::Zero(D);
        }
    }

    /// \brief Full constructor with back reference and control net.
    ///
    /// Basic constructor. All control points are the zero vector, and the
    /// back reference points to the given face.
    gsFreeformFaceData(
        gsMatrix<gismo::gsVector<real_t, D>, N, N> control_points, Face face)
        : control_points(control_points), face(face)
    {
    }

    /// \brief Basic constructor that creates a C0 mesh.
    ///
    /// Constructor that takes a mesh and chooses the control points of this
    /// evenly distributed over the the given face, creating a C0 mesh.
    gsFreeformFaceData(const gsSurfMesh& mesh, gsSurfMesh::Face face);

public: // Control point accessors
    /// \brief Returns control points along the edge of the face.
    ///
    /// Returns a matrix with pointers to all vectors in the control net.
    /// The matrix is oriented in such a way that the control point with indices
    /// `(0,0)` is right at the from-vertex of the given halfedge and the first
    /// row `(0,0)` to `(0, N-1)` follows in the direction of that halfedge,
    /// ending at its to-vertex with indices `(0, N-1)`.
    ///
    /// \param mesh The mesh this control net lives in.
    /// \param hedge The halfedge we orient this net on.
    gsMatrix<gsVector<real_t, D>*> control_points_oriented(gsSurfMesh& mesh,
                                                           Halfedge hedge);

public: // Conversions
    /// Returns a Bezier patch corresponding to these control points.
    const gismo::gsTensorBSpline<2, real_t> patch() const;

}; // namespace internal

} // namespace gismo
