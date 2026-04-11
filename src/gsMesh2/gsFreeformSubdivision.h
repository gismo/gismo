/** @file gsFreeformSubdivision.h

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#pragma once

#include "gsCore/gsFunctionExpr.h"
#include "gsIO/gsParaviewCollection.h"
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
template <size_t N>
class GISMO_EXPORT gsFreeformSubdivision : public gsSubdivisionScheme
{
private: // Space Dimension
    size_t D;

public: // Constructors
    /// \brief Default constructor.
    ///
    /// Default constructor. Sets no options and leaves the targeted mesh as a
    /// nullpointer.
    gsFreeformSubdivision(size_t D) : gsFreeformSubdivision(nullptr, D) {}

    /// \brief Constructor with a mesh to target.
    ///
    /// Constructs the object targeting the given mesh and registers the
    /// following options:
    /// - \c optimize_fit (switch, default \c false): when active, fits around
    ///   extraordinary vertices by optimising a functional instead of using
    ///   linear constraints loaded from file.
    /// - \c weighted_fit (switch, default \c false): when active, weights the
    ///   least-squares fit around extraordinary vertices using a per-sample
    ///   weight vector loaded from
    ///   \c filedata/freeform/val\<v\>_weights.xml.
    /// - \c model_patch_path (string, default \c "freeform/bubble/"): path
    ///   to the directory containing the model patch \c .xml files.
    ///
    /// \param mesh Pointer to the \c gsSurfMesh to be targeted by this object.
    gsFreeformSubdivision(gsSurfMesh* mesh, size_t D)
        : gsSubdivisionScheme(mesh), D(D)
    {
        m_options.addSwitch("optimize_fit",
                            "When active, fits around EVs by optimizing with "
                            "respect to a functional instead of using linear "
                            "constraints loaded from a file.",
                            false);
        m_options.addSwitch("weighted_fit",
                            "When active, weights the least-squares EV fit "
                            "using a per-sample weight vector loaded from "
                            "`filedata/freeform/val<v>_weights.xml`.",
                            false);
        m_options.addString("model_patch_path", "Path to the model patches.",
                            "freeform/bubble/");
    }

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
    /// \return An array of two control nets representing the left and right
    /// halves of the split patch.
    static std::array<gsMatrix<gsVector<real_t>, Dynamic, Dynamic>, 2>
    deCasteljau(
        const gsMatrix<gsVector<real_t>, Dynamic, Dynamic>& control_net);

    /// \brief Loads a model patch of the given valence and type.
    ///
    /// Reads the file \c Val<valence>Fct1.xml from the directory given by the
    /// \c model_patch_path option, selects the patch corresponding to
    /// \c subtype, drops the third coordinate (returns a 2D patch), and
    /// applies the rotation or scaling needed so that the \f$(0,0)\f$
    /// parameter point of every \c fine_* patch coincides with the central
    /// shared point of the four fine patches. The patches combine as follows:
    /// ```
    /// fine_4   fine_3
    /// fine_1   fine_2
    /// ```
    /// With the u/v-direction of \c coarse pointing right/up from the bottom
    /// left and the u-direction of each \c fine_* pointing outward from the
    /// central shared point with the same orientation (e.g. \c fine_1 points
    /// u leftwards and v downwards). The \c coarse patch is additionally
    /// scaled by a factor of 2.
    ///
    /// \param valence The valence of the desired patch. Valid values are 3, 5,
    /// 6, 7, 8, 9, 10.
    /// \param subtype The subtype of the patch. Valid values are \c "coarse",
    /// \c "fine_1", \c "fine_2", \c "fine_3", \c "fine_4".
    /// \return The loaded 2D Bézier patch as a \c gsTensorBSpline<2>.
    gismo::gsTensorBSpline<2, real_t> load_model_patch(int valence,
                                                       std::string subtype);

    /// \brief Re-orients the faces of the given mesh.
    ///
    /// This method changes the assigned halfedge of each face in the given mesh
    /// in such a way that, if the face has an extraordinary vertex, the
    /// halfedge of that face is the one pointing to that vertex. This will
    /// cause the first vertex of that face to be the EV. If any faces in the
    /// mesh have multiple EVs, there is no guarantee on which one will be
    /// chosen and it is suggested you change your mesh so that EVs have a
    /// larger distance between them. The control points of each affected face
    /// are rotated in lockstep with its halfedge so that their orientation
    /// remains consistent.
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
    /// \return The reordered array of faces, with the face containing
    /// \c first_vertex rotated to index 0.
    std::array<Face, 4> order_faces(Vertex first_vertex,
                                    std::array<Face, 4> faces);

    /// \brief Returns whether \c v is an ordinary vertex.
    ///
    /// A vertex is ordinary if it has valence 4 or lies on the boundary.
    ///
    /// \param mesh The surface mesh the vertex belongs to.
    /// \param v    The vertex to test.
    /// \return \c true if \c v is ordinary, \c false otherwise.
    inline bool is_ordinary(const gsSurfMesh& mesh, const Vertex& v) const
    {
        return mesh.valence(v) == 4 || mesh.is_boundary(v);
    }

    /// \brief Fits control points around an EV using linear constraints.
    ///
    /// Solves the regularised least-squares system \f$A x \approx
    /// \text{target}\f$ subject to linear constraints loaded from file, using a
    /// Tikhonov-augmented block system.
    ///
    /// \param A       Coefficient matrix of the fitting system.
    /// \param target  Right-hand side of the fitting system (D columns, one per
    ///                coordinate).
    /// \param valence Valence of the extraordinary vertex being fitted.
    /// \return Solution matrix \f$x\f$ satisfying the constraints.
    gsMatrix<real_t> fit_ev(gsMatrix<real_t> A, gsMatrix<real_t> target,
                            size_t valence);

    /// \brief Fits control points around an EV by kernel-space optimisation.
    ///
    /// Solves \f$A x = \text{target}\f$ via least squares and then optimises
    /// within the null space of \f$A\f$ to equalise the solution coefficients,
    /// minimising variation around the extraordinary vertex.
    ///
    /// \param A       Coefficient matrix of the fitting system.
    /// \param target  Right-hand side of the fitting system (D columns, one per
    ///                coordinate).
    /// \param valence Valence of the extraordinary vertex being fitted.
    /// \return Optimised solution matrix \f$x\f$.
    gsMatrix<real_t> fit_ev_opt(gsMatrix<real_t> A, gsMatrix<real_t> target,
                                size_t valence);

public:
    /// \brief Validates the mesh and returns its validity status.
    ///
    /// Checks whether every face of the targeted mesh has exactly 4 vertices
    /// (quad-only requirement). Returns \c UNDETERMINED (not \c VALID) when
    /// the mesh passes, since further geometric checks are not performed.
    /// Returns \c INVALID and emits a warning on the first non-quad face found.
    ///
    /// \return \c gsSubdivisionMeshValidity::UNDETERMINED if all faces are
    /// quads, \c gsSubdivisionMeshValidity::INVALID otherwise.
    gsSubdivisionScheme::gsSubdivisionMeshValidity check_mesh() override;

    /// \brief Performs one step of freeform subdivision on the mesh.
    ///
    /// Calls \c orient_faces(), then applies \c quad_split() to the underlying
    /// \c gsSurfMesh. For each original face the four child faces receive new
    /// \f$N \times N\f$ control nets according to one of two rules:
    /// - **Ordinary vertices** (all four corner valences equal 4, or on the
    ///   boundary): the de Casteljau algorithm is applied twice — once in each
    ///   parameter direction — to split the parent control net into four child
    ///   nets exactly.
    /// - **Extraordinary vertex** (one corner has valence \f$\neq 4\f$ and is
    ///   not on the boundary): model patches for the given valence are loaded
    ///   via \c load_model_patch(), the child control nets are obtained by
    ///   sampling the parent patch at model-patch parameter locations found
    ///   via Newton–Raphson, and fitting an \f$N \times N\f$ Bézier patch to
    ///   those samples.
    void subdivide() override;

    /// Use only garbage-collected meshes with contiguous indices!
    void basis_data(gsMultiPatch<>& multi_patch, gsMultiBasis<>& multi_basis,
                    gsMappedBasis<2>& mapped_basis);

    /// \brief Solves the Laplace-Beltrami problem on the mesh.
    ///
    /// Replaces the last coordinate of this patch with an approximate solution
    /// to the Laplace-Beltrami problem.  The discrete space is built using a
    /// \c gsMappedBasis with an identity mapping matrix, which is equivalent
    /// to the standard multi-basis and can later be replaced by a non-trivial
    /// mapping to impose patch-coupling or smoothness constraints.
    ///
    /// \param rhs The right hand side function.
    void laplace_beltrami(gsFunctionExpr<> rhs);

    /// \brief Adds another dimension to the control point vecs using with a
    /// least-squares fit to a function of the existing coordinates.
    ///
    /// For each face, samples the existing coordinates of the Bézier
    /// patch at \f$N^2\f$ parameter points, evaluates \c function at those
    /// positions, fits a new \f$N \times N\f$ Bézier patch to the function
    /// values by least squares, and writes the resulting control
    /// coefficients back as the new last coordinate of each control point.
    ///
    /// \param function A real-valued function in \f$D-1\f$ real variables.
    void fit_last_coordinate_to_function(gsFunctionExpr<> function);

    /// \brief Computes the approximation error of the freeform patches against
    /// a reference function.
    ///
    /// For each face, samples the Bézier patch at \c samples_per_face squared
    /// parameter points. At each sample the last coordinate of the patch is
    /// compared to the value of \c function evaluated at the first \f$D-1\f$
    /// coordinates. Returns a vector whose first entry is the \f$L^\infty\f$
    /// error and whose second entry is the root-mean-square \f$L^2\f$ error
    /// over all samples.
    ///
    /// \param function       A real-valued reference function in \f$D-1\f$
    /// variables.
    /// \param samples_per_face Number of sample points per parameter direction
    /// per face.
    /// \return A vector containing the \f$L^\infty\f$ and rms-\f$L^2\f$ errors.
    gsVector<real_t, 2> error(gsFunctionExpr<> function,
                              size_t samples_per_face);

    /// \brief Writes a per-face approximation-error field to Paraview.
    ///
    /// For each face, evaluates the pointwise absolute error between the last
    /// coordinate of the Bézier patch and \c function evaluated at the first
    /// \f$D-1\f$ coordinates. The error is rescaled to [0, 1] using \c
    /// max_error (which can be obtained via \c error()) fitted by a spline of
    /// degree \f$2(N-1)\f$, and written as a Paraview \c .vts file. If
    /// \c collection is non-null the generated file is registered in it at the
    /// given \c timestep.
    ///
    /// \note Only available when \c D >= 3. The spatial geometry shown in
    /// Paraview uses the first \f$\min(3, D-1)\f$ coordinates; for \c D > 4
    /// coordinates 4 and beyond are not visualised.
    ///
    /// \param function    A real-valued reference function in \f$D-1\f$
    /// variables.
    /// \param max_error   Value used to normalise the error to [0, 1].
    /// \param name        Base path/name for the output \c .vts files.
    /// \param collection  Optional Paraview collection to register the files
    /// in.
    /// \param timestep    Timestep index used when registering in \c
    /// collection.
    void write_paraview_error(gsFunctionExpr<real_t> function, real_t max_error,
                              std::string name,
                              gsParaviewCollection* collection = nullptr,
                              size_t timestep = 0);

    /// \brief Initializes the targeted mesh from an xml file.
    ///
    /// Clears the targeted mesh, then loads all \c gsTensorBSpline<2> objects
    /// from the \c .xml file at \c filepath. Each patch must have an
    /// \f$N \times N\f$ control net of \f$D\f$-dimensional vectors. For each
    /// patch a quad face is added to the mesh, with shared corner vertices
    /// detected by coordinate comparison using a tolerance of \f$10^{-10}\f$.
    /// The control net is stored in the \c bezier_points face property with
    /// the BSpline index transposition \f$(i,j) \to (j,i)\f$ applied to match
    /// the internal face layout convention.
    ///
    /// \param filepath Path to the \c .xml file to load.
    void initialize_data_xml(std::string filepath);

    /// \brief Initializes the targeted mesh from an off file.
    ///
    /// Clears the targeted mesh, loads the \c .off file at \c filepath
    /// (which contains 3-dimensional point data), adds the
    /// \c bezier_points face property, and initializes each face's control
    /// net via bilinear interpolation of its four corner positions. The first
    /// \f$\min(D,3)\f$ coordinates are taken from the mesh; any coordinates
    /// beyond index 2 are initialised to zero.
    ///
    /// \param filepath Path to the \c .off file to load.
    void initialize_data_off(std::string filepath);

    /// \brief Initializes the targeted mesh from a file, dispatching on
    /// extension.
    ///
    /// Detects the file format from the extension of \c filepath and calls the
    /// appropriate loader:
    /// - \c .xml → \ref initialize_data_xml
    /// - \c .off → \ref initialize_data_off
    ///
    /// If the extension is not recognised, a warning is emitted and the mesh
    /// is left unchanged.
    ///
    /// \param filepath Path to the file to load (\c .xml or \c .off).
    void initialize_data(std::string filepath);

    /// \brief Initializes the targeted mesh from a file, dispatching on
    /// extension, and sets the dimension of the mesh.
    ///
    /// First sets the dimension of this freeform subdivision to D.
    /// Then detects the file format from the extension of \c filepath and calls
    /// the appropriate loader:
    /// - \c .xml → \ref initialize_data_xml
    /// - \c .off → \ref initialize_data_off
    ///
    /// If the extension is not recognised, a warning is emitted and the mesh
    /// is left unchanged.
    ///
    /// \param filepath Path to the file to load (\c .xml or \c .off).
    /// \param D The dimension of the mesh.
    void initialize_data(std::string filepath, size_t D);

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

    /// \brief Turns a $C^0$ set of control nets into a $C^s$ set.
    ///
    /// Takes a given mesh with freeform data and makes it $C^s$ by adjusting
    /// the outer layer of control points of each bezier patch. This only causes
    /// $C^s$ at edges and ordinary vertices. No guarantee is made for
    /// extraordinary vertices.
    /// Features return arguments that return the coefficient matrix of the
    /// fitting for each extraordinary vertex, as well as the matrix containing
    /// the values of the outer control points of patches around that vertex.
    ///
    /// \param degree The degree of smoothness desired. As of now, only $C^1$ is
    /// supported.
    /// \param ev_coefficients      Output: coefficient matrices of the EV fit,
    ///                             one entry per extraordinary vertex.
    /// \param ev_coefficients_outer Output: matrices of the outer control point
    ///                              values around each extraordinary vertex.
    void smooth(size_t degree, std::vector<gsMatrix<real_t>>& ev_coefficients,
                std::vector<gsMatrix<real_t>>& ev_coefficients_outer);

    /// \brief Converts to a Gismo multipatch object.
    ///
    /// Converts the targeted mesh with freeform data into a multipatch that can
    /// be easily displayed by e.g. Paraview. Each face and its control net are
    /// converted to one appropriately sized patch.
    ///
    /// \return A \c gsMultiPatch<> containing one patch per face.
    gsMultiPatch<> multipatch();

    /// \brief Writes the targeted mesh with freeform data to a Paraview file.
    ///
    /// Writes the targeted mesh with freeform data to a paraview file for easy
    /// viewing. Can also register the files to given paraview collections at
    /// the given timestep.
    ///
    /// \param name             Base path/name for the output \c .vts files.
    /// \param collection       Optional Paraview collection to register the
    ///                         surface patch files in.
    /// \param cnet_collection  Optional Paraview collection to register the
    ///                         control net files in. Only used when
    ///                         \c control_net is true.
    /// \param timestep         Timestep index used when registering in
    ///                         \c collection and \c cnet_collection.
    /// \param control_net      When set to true, the control net of the mesh is
    ///                         also written to paraview and, if
    ///                         \c cnet_collection is provided, registered
    ///                         there.
    void write_paraview(std::string name,
                        gsParaviewCollection* collection = nullptr,
                        gsParaviewCollection* cnet_collection = nullptr,
                        size_t timestep = 0, bool control_net = false);

}; // namespace internal

/// \brief Manages a control net on a face of a surf mesh.
///
/// Data class that contains the additional information per face needed to
/// perform freeform subdivision. Represents an \f$N \times N\f$ grid of
/// Bézier control points on the face. The grid is supposed to be laid out
/// as follows:
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
/// where `V0`–`V3` are the vertices in the order `mesh.vertices(face)`
/// traverses them and `E0`–`E3` are the (half)edges in the order
/// `mesh.halfedges(face)` traverses them. In particular, the position of
/// control points `00`, `04`, `40`, `44` should correspond to the position of
/// the vertices (if present).
template <size_t N> class GISMO_EXPORT gsFreeformFaceData
{
    template <size_t M> friend class gsFreeformSubdivision;

    using Point = gsSurfMesh::Point;
    using Vertex = gsSurfMesh::Vertex;
    using Face = gsSurfMesh::Face;
    using Halfedge = gsSurfMesh::Halfedge;
    using Edge = gsSurfMesh::Edge;

private: // members
    /// The \f$N \times N\f$ matrix of \f$D\f$-dimensional Bézier control
    /// points of this face patch.
    gsMatrix<gismo::gsVector<real_t>, N, N> control_points;
    /// A back-reference to the face of the mesh this data belongs to.
    Face face;
    /// Dimension of the space
    size_t D;

public: // Contructors
    /// \brief Default constructor with everything empty.
    ///
    /// Default constructor. All control points are the zero vector, and the
    /// back reference to the face is empty.
    gsFreeformFaceData(size_t D) : control_points(), face(0), D(D)
    {
        for (size_t i = 0; i < N * N; ++i)
        {
            control_points(i) = gsVector<real_t>::Zero(D);
        }
    }

    /// \brief Zero constructor with a back reference.
    ///
    /// All control points are set to the zero vector, and the back reference
    /// points to the given face.
    ///
    /// \param face The face this data is associated with.
    gsFreeformFaceData(Face face, size_t D) : control_points(), face(face), D(D)
    {
        for (size_t i = 0; i < N * N; ++i)
        {
            control_points(i) = gsVector<real_t>::Zero(D);
        }
    }

    /// \brief Full constructor with back reference and control net.
    ///
    /// Constructs the face data with the given control net and face reference.
    ///
    /// \param control_points The \f$N \times N\f$ matrix of \f$D\f$-dimensional
    ///                       Bézier control points for this face.
    /// \param face           The face this data is associated with.
    gsFreeformFaceData(gsMatrix<gismo::gsVector<real_t>, N, N> control_points,
                       Face face)
        : control_points(control_points), face(face)
    {
        this->D = control_points(0, 0).size();
    }

    /// \brief Basic constructor that creates a $C^0$ control net.
    ///
    /// Constructs the face data with control points distributed by bilinear
    /// interpolation over the four corners of the given mesh face. The first
    /// \f$\min(D,3)\f$ coordinates of each control point are taken from the
    /// 3-dimensional mesh vertex positions; any coordinates beyond index 2
    /// are initialised to zero.
    ///
    /// \param mesh The surface mesh containing the face.
    /// \param face The face over which to distribute the control points.
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
    /// \return An \f$N \times N\f$ matrix of pointers into the control net,
    /// oriented so that
    /// row 0 runs from the from-vertex to the to-vertex of \c hedge.
    gsMatrix<gsVector<real_t>*> control_points_oriented(gsSurfMesh& mesh,
                                                        Halfedge hedge);

public: // Conversions
    /// \brief Returns a Bézier patch corresponding to these control points.
    ///
    /// \return A \c gsTensorBSpline<2> of degree N-1 in each direction,
    /// built from the stored N×N control net.
    const gismo::gsTensorBSpline<2, real_t> patch() const;

}; // namespace internal

} // namespace gismo
