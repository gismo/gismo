/** @file gsFrogSplines.h

    @brief Classes for FROG Splines on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#pragma once

#include "gsCore/gsFunctionExpr.h"
#include "gsCore/gsMultiPatch.h"
#include "gsIO/gsParaviewCollection.h"
#include "gsMatrix/gsVector.h"
#include <gsMesh2/gsSubdivisionScheme.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

/// \brief Manager class for FROG splines.
///
/// This class manages a quadrangular mesh and basises of FROG splines, allowing
/// for subdivision or fitting operations.
template <size_t N>
class GISMO_EXPORT gsFrogSplines : public gsSubdivisionScheme
{
private: // Space Dimension
    size_t D;

    using gsPatch = gsTensorBSpline<2>;

public: // Constructors
    /// \brief Default constructor.
    ///
    /// Default constructor. Sets no options and leaves the targeted mesh as a
    /// nullpointer.
    gsFrogSplines(size_t D) : gsFrogSplines(nullptr, D) {}

    /// \brief Constructor with a mesh to target.
    ///
    /// Constructs the object targeting the given mesh and registers the
    /// following options:
    /// - \c optimize_fit (switch, default \c false): when active, chooses
    ///   extraordinary-vertex coefficient representatives by minimising the
    ///   smoothing functional used by \c fit_ev_opt(); otherwise it drives the
    ///   linear functionals loaded from \c Val<v>Constraints.xml to zero.
    /// - \c weighted_fit (switch, default \c false): when active, weights the
    ///   least-squares fit around extraordinary vertices using a per-sample
    ///   weight vector loaded from
    ///   \c filedata/frog/val\<v\>_weights.xml.
    /// - \c frog_dir (string, default \c "frog/bubble/"): path
    ///   to the directory containing a set of frog spline generating
    ///   functions for each required valence.
    ///
    /// \param mesh Pointer to the \c gsSurfMesh to be targeted by this object.
    gsFrogSplines(gsSurfMesh* mesh, size_t D) : gsSubdivisionScheme(mesh), D(D)
    {
        m_options.addSwitch("optimize_fit",
                            "When active, chooses EV coefficient "
                            "representatives by minimizing the functional "
                            "used by fit_ev_opt(); otherwise it drives the "
                            "linear functionals loaded from file to zero.",
                            false);
        m_options.addSwitch("weighted_fit",
                            "When active, weights the least-squares EV fit "
                            "using a per-sample weight vector loaded from "
                            "`filedata/frog/val<v>_weights.xml`.",
                            false);
        m_options.addString("frog_dir", "Path to the directory containing a set of frog spline generating functions for each required valence.",
                            "frog/bubble/");
    }

public: // Data initialization
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
    /// First sets the dimension of this frog spline subdivision to D.
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

    /// \brief Initializes the targeted mesh from an xml file.
    ///
    /// Reads a \c gsMultiPatch from the \c .xml file at \c filepath and
    /// delegates all mesh and face-data construction to
    /// \ref initialize_data_multipatch.
    ///
    /// \param filepath Path to the \c .xml file to load.
    void initialize_data_xml(std::string filepath);

    /// \brief Initializes the targeted mesh from a gsMultiPatch.
    ///
    /// Builds the combinatorial mesh and Bézier face data from \c mpatch.
    /// The method proceeds in two steps:
    ///
    /// **Topology construction.**
    /// For each patch, only the four corner coefficients (indices \c 0,
    /// \c N-1, \c N*(N-1), \c N*N-1) are extracted to form a bilinear proxy
    /// patch. All proxies are collected into a temporary \c gsMultiPatch,
    /// whose topology is computed and then converted to the internal
    /// \c gsSurfMesh via \c toMesh(). This yields one quad face per input
    /// patch with correct vertex-sharing between adjacent patches.
    ///
    /// **Face data initialisation.**
    /// A \c bezier_points face property is added to the mesh, typed as
    /// \c gsPatch with an \f$N \times N\f$ Bézier basis and a zero
    /// \f$D\f$-column coefficient matrix as default. Each face is then
    /// assigned the corresponding patch from \c mpatch: the first
    /// \f$\min(D, \mathrm{geoDim})\f$ coefficient columns are copied and
    /// any remaining columns are left as zero.
    ///
    /// \param mpatch A \c gsMultiPatch whose patch ordering matches the
    ///               desired face ordering of the mesh. Each patch must have
    ///               at least \f$N \times N\f$ control points.
    void initialize_data_multipatch(gsMultiPatch<> mpatch);

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

    /// \brief Scales all vertex positions and Bézier control points.
    ///
    /// Applies a per-axis scaling to every vertex position in the mesh and to
    /// every Bézier control point stored in the face data. Each coordinate
    /// \f$x_i\f$ is multiplied by the corresponding entry \c factors[i].
    ///
    /// \param factors A 3-vector whose components are the scale factors for
    ///                the x, y, and z axes respectively.
    void scale(gsVector3d<> factors);

private: // Helper functions
    /// \brief Rotates a control net clockwise (equivalent to a
    /// counter-clockwise
    ///        geometric rotation of the surface patch).
    ///
    /// The control net is stored as a flat \f$(n^2 \times D)\f$ coefficient
    /// matrix in row-major order (flat index \f$i \cdot n + j\f$ for grid
    /// position \f$(i,j)\f$). The grid position \f$(i,j)\f$ is mapped to
    /// \f$(n-1-j,\, i)\f$.
    ///
    /// \note Follows the same naming convention as \c gsMatrix::rotate_cw(),
    ///       where "cw" refers to the rotation of the index grid, not the
    ///       resulting geometric orientation.
    ///
    /// \param coefs The flat coefficient matrix to rotate.
    /// \param n     Side length of the square grid.
    /// \return A new coefficient matrix with the rotated layout.
    static gsMatrix<> rotate_coefs_cw(const gsMatrix<>& coefs, size_t n);

    /// \brief Rotates a control net counter-clockwise (equivalent to a
    /// clockwise
    ///        geometric rotation of the surface patch).
    ///
    /// The control net is stored as a flat \f$(n^2 \times D)\f$ coefficient
    /// matrix in row-major order. The grid position \f$(i,j)\f$ is mapped to
    /// \f$(j,\, n-1-i)\f$.
    ///
    /// \note Follows the same naming convention as \c gsMatrix::rotate_ccw(),
    ///       where "ccw" refers to the rotation of the index grid, not the
    ///       resulting geometric orientation.
    ///
    /// \param coefs The flat coefficient matrix to rotate.
    /// \param n     Side length of the square grid.
    /// \return A new coefficient matrix with the rotated layout.
    static gsMatrix<> rotate_coefs_ccw(const gsMatrix<>& coefs, size_t n);

    /// \brief Loads a model patch of the given valence and type.
    ///
    /// Reads the file \c Val<valence>Fct1.xml from the directory given by the
    /// \c frog_dir option, selects the patch corresponding to
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

    /// \brief Converts to a Gismo multipatch object.
    ///
    /// Converts the targeted mesh with patch data into a multipatch that can
    /// be easily displayed by e.g. Paraview. Each face and its control net are
    /// converted to one appropriately sized patch.
    ///
    /// \return A \c gsMultiPatch<> containing one patch per face.
    gsMultiPatch<> multipatch();

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

public: // Subdivision
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

    /// \brief Performs one step of frog-based subdivision on the mesh.
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

public: // Fitters
    /// \brief Smooths the geometry by L2-projecting it onto the mapped basis.
    ///
    /// Builds the mapped basis via \c c1_basis(), projects the current
    /// multipatch geometry onto that space in the \f$L^2\f$ sense, and then
    /// post-processes each extraordinary-vertex coefficient block in the kernel
    /// of the corresponding blending functions. When \c optimize_fit is
    /// enabled, the post-processing minimises the equalisation functional used
    /// by \c fit_ev_opt(); otherwise it enforces the linear functionals loaded
    /// from \c Val<v>Constraints.xml. In both cases the kernel basis is loaded
    /// from \c Val<v>Kernel.xml. The final projected control points are written
    /// back to the per-face control nets.
    ///
    /// \param degree Requested smoothness degree. Currently only implemented
    /// for C1.
    /// \return Matrix of mapped coefficients with one row per global mapped
    /// degree of freedom and one column per geometric coordinate.
    gsMatrix<real_t> smooth(size_t degree);

    /// \brief Appends a scalar field as a new control-point coordinate.
    ///
    /// Evaluates \c function on the current geometry and fits the resulting
    /// scalar values with the frog patch representation. The fitted values
    /// are written into a new last coordinate of every control point, and the
    /// ambient dimension \c D is increased by one.
    ///
    /// \param function A real-valued function of the current geometric
    ///                 coordinates.
    void fit_function(gsFunctionExpr<> function);

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

public: // Basises
    /// \brief Builds the mapped \f$C^1\f$ basis used by smoothing and PDE
    /// solves.
    ///
    /// Converts the current frog control-net mesh into a \c gsMultiPatch
    /// and its underlying \c gsMultiBasis, then assembles a sparse
    /// local-to-global mapper whose image defines the \f$C^1\f$ spline space.
    /// The construction proceeds in four stages:
    /// - interior extraordinary vertices reserve degrees of freedom for the EV
    ///   blending-function and map all local control points inside the EV
    ///   support through the corresponding model-patch coefficient rows,
    /// - ordinary interior control points receive free global degrees of
    ///   freedom,
    /// - interior edge points are expressed by the \f$C^1\f$ averaging relation
    ///   between the two adjacent inner rows, while boundary edge points remain
    ///   free,
    /// - corner points are handled as follows:
    ///   - non-boundary vertices are averaged from the surrounding
    ///     \f$(1,1)\f$ inner points,
    ///   - a boundary corner with one adjacent patch remains free,
    ///   - a boundary corner with two adjacent patches is averaged from the
    ///     two neighboring control points on the boundary edges,
    ///   - boundary corners with more than two adjacent patches are averaged
    ///     from one neighboring edge control point per outgoing edge
    ///     (equivalently \f$\#\mathrm{patches}+1\f$ points).
    ///
    /// The resulting global degrees of freedom are ordered as follows:
    /// - for each extraordinary vertex in mesh-vertex iteration order:
    ///   - the central frog blending function,
    ///   - the remaining frog blending functions around that vertex,
    ///   - the additional free control-point degrees of freedom around that EV,
    /// - ordinary inner points in face order,
    /// - edge points in face/halfedge order,
    /// - corner points.
    ///
    /// \warning Use only garbage-collected meshes with contiguous face and
    /// vertex indices.
    ///
    /// \param multi_patch Output multipatch containing one Bézier patch per
    /// mesh face.
    /// \param multi_basis Output tensor-product basis underlying \p
    /// multi_patch.
    /// \param mapped_basis Output mapped basis that encodes the \f$C^1\f$
    /// coupling relations.
    /// \param constraint_matrix Output sparse matrix encoding linear
    /// constraints on global DOFs near extraordinary vertices.
    /// \param use_generating_system If true (default), use all fitting functions
    /// including linearly dependent ones near EVs (generating system). If false,
    /// remove linearly dependent functions to obtain a proper basis.
    void c1_basis(gsMultiPatch<>& multi_patch, gsMultiBasis<>& multi_basis,
                  gsMappedBasis<2>& mapped_basis, bool use_generating_system = true);

public: // errors & output
    /// \brief Computes the approximation error of the patch patches against
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

    /// \brief Writes the targeted mesh with patch data to a Paraview file.
    ///
    /// Writes the targeted mesh with patch data to a paraview file for easy
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

}; // class gsFrogSplines

} // namespace gismo
