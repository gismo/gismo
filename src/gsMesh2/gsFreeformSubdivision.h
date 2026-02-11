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
#include <vector>

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
    /// Default constructor. Only constructor, as this scheme has no special
    /// options.
    gsFreeformSubdivision() : gsSubdivisionScheme() {}

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

    /// \brief Checks if a vertex is ordinary.
    ///
    /// Checks if the given vertex has is ordinary, i.e.
    ///  - has exactly 4 neighbors OR
    ///  - lies on the boundary of the mesh.
    ///
    /// \param mesh The surrounding mesh this vertex belongs to.
    /// \param v The vertex to be analyzed.
    static bool is_ordinary(const gsSurfMesh& mesh, const Vertex& v);

    static gismo::gsTensorBSpline<2, real_t> load_patch(int valence,
                                             std::string subtype);



public:
    gsSubdivisionScheme::gsSubdivisionMeshValidity
    valid_mesh(const gsSurfMesh&) override;

    void subdivide(gsSurfMesh& mesh) override;

    /// \brief Initializes C0 control nets on all faces.
    ///
    /// Takes a given mesh without freeform data and initializes a freeform data
    /// control net with C0 continuity on each face.
    void initialize_data(gsSurfMesh& mesh);

    /// \brief Turns a $C^0$ set of control nets into a $C^s$ set.
    ///
    /// Takes a given mesh with freeform data and makes it $C^s$ by adjusting the
    /// outer layer of control points of each bezier patch. This only causes $C^s$
    /// at edges and ordinary vertices. No guarantee is made for extraordinary
    /// vertices.
    ///
    /// \param degree The degree of smoothness desired. As of now, only $C^1$ is supported.
    void smooth(gsSurfMesh& mesh, size_t degree);

    /// \brief Converts to a Gismo multipatch object.
    ///
    /// Converts a given mesh with freeform data into a multipatch that can be
    /// easily displayed by e.g. Paraview. Each face and its control net are
    /// converted to one appropriately sized patch.
    gsMultiPatch<> multipatch(const gsSurfMesh& mesh);

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
    /// row `(0,0)` to `(N-1,0)` follows in the direction of that halfedge,
    /// ending at its to-vertex with indices `(N-1,0)`.
    ///
    /// \param mesh The mesh this control net lives in.
    /// \param hedge The halfedge we orient this net on.
    gsMatrix<gsVector<real_t, D>*> control_points_oriented(gsSurfMesh& mesh,
                                                           Halfedge hedge);

public: // Conversions
    /// Returns a Bezier patch corresponding to these control points.
    const gismo::gsTensorBSpline<2, real_t> patch() const;

}; // namespace internal

/// \brief Matrix rotation.
///
/// Rotates a leftwise matrix around its center, returning the rotated matrix
/// without changing the original.
///
/// \param mat The matrix to be rotated.
template <class T>
gsMatrix<T, Dynamic, Dynamic> rotate_l(const gsMatrix<T, Dynamic, Dynamic>& mat)
{
    gsMatrix<T, Dynamic, Dynamic> res;
    res.resize(mat.cols(), mat.rows());
    for (int i = 0; i < res.rows(); ++i)
    {
        for (int j = 0; j < res.cols(); ++j)
        {
            res(i, j) = mat(j, mat.cols() - 1 - i);
        }
    }
    return res;
}

/// \brief Matrix rotation.
///
/// Rotates a matrix rightwise around its center, returning the rotated matrix
/// without changing the original.
///
/// \param mat The matrix to be rotated.
template <class T>
gsMatrix<T, Dynamic, Dynamic> rotate_r(const gsMatrix<T, Dynamic, Dynamic>& mat)
{
    gsMatrix<T, Dynamic, Dynamic> res;
    res.resize(mat.cols(), mat.rows());
    for (int i = 0; i < res.rows(); ++i)
    {
        for (int j = 0; j < res.cols(); ++j)
        {
            res(i, j) = mat(mat.rows() - 1 - j, i);
        }
    }
    return res;
}

} // namespace gismo
