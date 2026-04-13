/** @file gsFreeformSubdivision.cpp

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#include "gsAssembler/gsAssemblerOptions.h"
#include "gsAssembler/gsExprAssembler.h"
#include "gsCore/gsField.h"
#include "gsCore/gsFunctionExpr.h"
#include "gsCore/gsMultiBasis.h"
#include "gsIO/gsParaviewCollection.h"
#include "gsIO/gsWriteParaview.h"
#include "gsMSplines/gsMappedBasis.h"
#include "gsMSplines/gsMappedSpline.h"
#include "gsMatrix/gsMatrix.h"
#include "gsMatrix/gsSparseSolver.h"
#include "gsMatrix/gsVector.h"
#include "gsPde/gsBoundaryConditions.h"
#include "gsUtils/gsL2Projection.h"
#include <cmath>
#include <cstdlib>
#include <gismo.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsMultiPatch.h>
#include <gsIO/gsFileData.h>
#include <gsIO/gsFileData.hpp>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSubdivisionScheme.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsModeling/gsFitting.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <memory>

namespace gismo
{

template <size_t N>
gsFreeformFaceData<N>::gsFreeformFaceData(const gsSurfMesh& mesh,
                                          gsSurfMesh::Face face)
    : control_points(), face(face)
{
    // Create a vector of the 4 corners, copying the 3D mesh vertex positions
    // into D-dimensional vectors. The first min(D,3) coordinates are taken
    // from the mesh; any remaining coordinates beyond index 2 are zeroed.
    std::vector<gsVector<real_t>> corners;
    corners.reserve(4);
    for (auto const& v : mesh.vertices(face))
    {
        gsVector<real_t> corner = gsVector<real_t>::Zero(D);
        const gsSurfMesh::Point& p = mesh.position(v);
        for (size_t k = 0; k < std::min<size_t>(D, 3); ++k)
            corner(k) = p[k];
        corners.emplace_back(corner);
    }
    assert(corners.size() == 4);

    // Choose the control points (N*N total) as appropriate linear combinations
    // of the corners.
    real_t denom = real_t((N - 1) * (N - 1));
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            real_t n_r = real_t(N);
            real_t i_r = real_t(i);
            real_t j_r = real_t(j);
            this->control_points(i, j) =
                corners[0] * ((n_r - 1) - j_r) * ((n_r - 1) - i_r) / denom +
                corners[1] * j_r * ((n_r - 1) - i_r) / denom +
                corners[3] * ((n_r - 1) - j_r) * i_r / denom +
                corners[2] * j_r * i_r / denom;
        }
    }
}

template <size_t N>
gsMatrix<gsVector<real_t>*>
gsFreeformFaceData<N>::control_points_oriented(gsSurfMesh& mesh, Halfedge hedge)
{

    gsMatrix<gsVector<real_t>*> result;
    result.resize(control_points.rows(), control_points.cols());
    for (int i = 0; i < control_points.rows(); ++i)
    {
        for (int j = 0; j < control_points.cols(); ++j)
        {
            result(i, j) = &control_points(i, j);
        }
    }
    result = result.rotate_ccw();
    // find the edge on the face
    for (auto const& he : mesh.halfedges(face))
    {
        if (he == hedge)
            break;

        result = result.rotate_cw();
    }
    return result;
}

template <size_t N>
const gismo::gsTensorBSpline<2, real_t> gsFreeformFaceData<N>::patch() const
{
    // Create a spline basis for a normal bezier patch.
    gsKnotVector<> kv(0, 1, 0, N);
    gsTensorBSplineBasis<2, real_t> basis(kv, kv);
    // Create a coefficient matrix out of the control points.
    // Technically, you could use just [i] and one loop here since the elements
    // of a matrix are laid out row-wise, but this might be clearer to read.
    gsMatrix<> coeffs(N * N, D);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            int total_index = i * N + j;
            coeffs.row(total_index) = control_points(i, j);
        }
    }

    return gsTensorBSpline<2>(basis, coeffs);
}

template <size_t N>
std::array<gsMatrix<gsVector<real_t>, Dynamic, Dynamic>, 2>
gsFreeformSubdivision<N>::deCasteljau(
    const gsMatrix<gsVector<real_t>, Dynamic, Dynamic>& control_net)
{

    // Create the 3d data vector and ensure it has the right size.
    std::vector<gsMatrix<gsVector<real_t>, Dynamic, Dynamic>> points;
    points.resize(N);
    // The first layer is just the starting points
    points[0] = control_net;
    // each further layer is one shorter than the previous, but just as wide.
    for (size_t k = 1; k < N; ++k)
    {
        points[k].resize(N - k, N);
    }

    // now construct each layer from the previous one by linear combination of
    // adjacent points into the next layer.
    for (size_t k = 1; k < N; ++k)
    {
        for (size_t i = 0; i < N - k; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                points[k](i, j) =
                    (points[k - 1](i, j) + points[k - 1](i + 1, j)) * 0.5;
            }
        }
    }

    // finally collect the first vertical layer and last-in-each-row diagonal
    // layer into two result matrices.
    gsMatrix<gsVector<real_t>, Dynamic, Dynamic> result1;
    gsMatrix<gsVector<real_t>, Dynamic, Dynamic> result2;
    result1.resize(N, N);
    result2.resize(N, N);

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            result1(i, j) = points[i](0, j);
            result2(i, j) = points[(N - 1) - i](i, j);
        }
    }

    return {result1, result2};
}

template <size_t N>
gismo::gsTensorBSpline<2, real_t>
gsFreeformSubdivision<N>::load_model_patch(int valence, std::string subtype)
{
    // Load all patches from Val<valence>Fct1.xml
    auto path = m_options.getString("model_patch_path") + "Val" +
                std::to_string(valence) + "Fct1.xml";
    auto patches =
        gsFileData<real_t>(path).getAll<gsTensorBSpline<2, real_t>>();

    // Select the correct patch index based on subtype
    size_t patch_index(0);
    if (subtype == "fine_1" || subtype == "coarse")
        patch_index = 0;
    else if (subtype == "fine_2")
        patch_index = valence + 2;
    else if (subtype == "fine_3")
        patch_index = valence + 1;
    else if (subtype == "fine_4")
        patch_index = valence + 0;

    // Drop the z-coordinate: keep only the first 2 columns of the coef matrix
    // Then split into two 1-column matrices
    gsMatrix<real_t> coefsx = patches[patch_index]->coefs().col(0);
    gsMatrix<real_t> coefsy = patches[patch_index]->coefs().col(1);
    // reshape each to a quadratic matrix
    coefsx = coefsx.reshape(N, N);
    coefsy = coefsy.reshape(N, N);

    // use rotation/scaling based on subtype to rotate correctly
    // in a way that the (0,0) point for all fine patches is the central point.
    if (subtype == "coarse")
    {
        coefsx *= 2.0;
        coefsy *= 2.0;
    }
    else if (subtype == "fine_1")
    {
        coefsx = coefsx.rotate_ccw().rotate_ccw();
        coefsy = coefsy.rotate_ccw().rotate_ccw();
    }
    else if (subtype == "fine_2")
    {
        coefsx = coefsx.rotate_cw();
        coefsy = coefsy.rotate_cw();
    }
    else if (subtype == "fine_3")
    {
        coefsx = coefsx;
        coefsy = coefsy;
    }
    else if (subtype == "fine_4")
    {
        coefsx = coefsx.rotate_ccw();
        coefsy = coefsy.rotate_ccw();
    }

    // re-reshape back to a 1-column matrix
    coefsx = coefsx.reshape(N * N, 1);
    coefsy = coefsy.reshape(N * N, 1);

    // recombine into a 2-column matrix
    gsMatrix<real_t> coefs(N * N, 2);
    coefs.col(0) = coefsx;
    coefs.col(1) = coefsy;

    // build the tensor spline
    return gsTensorBSpline<2>(patches[patch_index]->basis(), give(coefs));
}

template <size_t N>
std::array<gsSurfMesh::Face, 4>
gsFreeformSubdivision<N>::order_faces(Vertex first_vertex,
                                      std::array<gsSurfMesh::Face, 4> faces)
{
    auto& mesh = *m_mesh;
    size_t first_face(0);
    for (size_t i = 0; i < 4; ++i)
    {
        // get the face
        Face f(faces[i]);
        // go through the vertices of this face and check if it has the
        // searched-for vertex
        for (auto const& v : mesh.vertices(f))
        {
            if (v == first_vertex)
            {
                first_face = i;
                break;
            }
        }
    }

    // Collate the faces into a correctly ordered array as well.
    std::array<Face, 4> children_faces_ordered;
    for (int i = 0; i < 4; ++i)
    {
        children_faces_ordered[i] = faces[(i + first_face) % 4];
    }

    return children_faces_ordered;
}

template <size_t N> void gsFreeformSubdivision<N>::orient_faces()
{
    // Get data
    auto& mesh = *m_mesh;
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points"));

    // As a pre-step, we make sure that for each face, if it has an
    // extraordinary vertex, that vertex is the top left one
    for (Face f : mesh.faces())
    {
        // first, check if this face has an EV at all
        bool has_ev(false);
        for (Vertex v : mesh.vertices(f))
        {
            has_ev |= !is_ordinary(mesh, v);
        }

        // if it does, rotate its first halfedge around until it points to the
        // EV
        if (has_ev)
        {
            while (is_ordinary(mesh, mesh.to_vertex(mesh.halfedge(f))))
            {
                // rotate the edge
                mesh.set_halfedge(f, mesh.next_halfedge(mesh.halfedge(f)));
                // also rotate the control points
                gsMatrix<gsVector<real_t>, Dynamic, Dynamic> old_points(
                    face_data_vec.vector()[f.idx()].control_points);
                face_data_vec.vector()[f.idx()].control_points =
                    old_points.rotate_cw();
            }
        }
    }
}

template <size_t N> void gsFreeformSubdivision<N>::subdivide()
{
    auto& mesh = *m_mesh;
    // First, make sure all faces are correctly oriented with the EV as their
    // first vertex.
    // Remember the number of faces before the subdivision.
    size_t n(mesh.faces_size());

    orient_faces();
    // Remember the first vertex of each face (this is where the control nets of
    // each face data are oriented on).
    std::vector<Vertex> first_vertices;
    for (Face f : mesh.faces())
    {
        // remember the first vertex
        first_vertices.emplace_back(mesh.to_vertex(mesh.halfedge(f)));
    }

    // Do the quad split of the (abstract) base faces.
    mesh.quad_split();

    // Get face data
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points"));

    // Now fix the data on each face.
    for (size_t initial_face = 0; initial_face < n; ++initial_face)
    {
        Vertex fv = first_vertices[initial_face];

        // The indices of the child faces can be calculated from the guarantees
        // made by `quad_split`.
        std::array<Face, 4> child_faces = {
            Face(initial_face), Face(n + initial_face * 3),
            Face(n + initial_face * 3 + 1), Face(n + initial_face * 3 + 2)};

        // Check if we have an EV. If there is, it must be the first vertex
        // thanks to our preprocessing.
        if (!is_ordinary(mesh, fv))
        {
            // === EXTRAORDINARY VERTICES ===
            // via lots of fitting

            // load coarse matrix
            auto valence = mesh.valence(fv);
            auto coarse_model = load_model_patch(valence, "coarse");
            // and remember the associated face data
            auto coarse_patch = face_data_vec.vector()[initial_face].patch();

            // Collate the faces into a correctly ordered array.
            auto children_ordered =
                order_faces(first_vertices[initial_face], child_faces);

            // Now fit each face with a new control net.
            for (size_t f_idx = 0; f_idx < 4; ++f_idx)
            {
                // first, load the appropriate fine model control points
                auto fine_model = load_model_patch(
                    valence, "fine_" + std::to_string(f_idx + 1));

                // get sample points
                gsMatrix<> samples(D, N * N);
                gsMatrix<> params(2, N * N);

                for (size_t i = 0; i < N * N; ++i)
                {

                    // The parameter of the sample point.
                    params.col(i) = gsVector<real_t, 2>::vec(
                        real_t(std::floor(i % N)) / real_t(N - 1),
                        real_t(std::floor(i / N)) / real_t(N - 1));
                    // The point in the geometry on the fine model patch.
                    gsVector<real_t, 2> point = fine_model.eval(params.col(i));

                    gsVector<real_t> closest_point =
                        gsVector<real_t>::vec(0.5, 0.5);
                    // Get the parameters of the same point in the coarse
                    // geometry model via Newton-Raphson. Note that internally,
                    // the tolerance is squared, so this is a tolerance of
                    // 1e-12.
                    coarse_model.closestPointTo(point, closest_point, 1e-6,
                                                true);

                    // Sample the old control net.
                    samples.col(i) = coarse_patch.eval(closest_point);
                    // if(f_idx == 0)
                    //     gsInfo << i << ": " << closest_point.transpose() <<
                    //     "\n";
                }

                // Fit a NxN Bezier patch to these samples
                gsKnotVector<> kv(0, 1, 0, N);
                gsTensorBSplineBasis<2> basis(kv, kv);
                gsFitting<> fitter(params, samples, basis);
                fitter.compute(0.0);
                gsGeometry<>* result = fitter.result();

                // Extract control points
                const gsMatrix<>& coeffs = result->coefs();

                // Reshape the control points
                gsMatrix<gsVector<real_t>, Dynamic, Dynamic> control_points;
                control_points.resize(N, N);
                for (size_t i = 0; i < N * N; ++i)
                {
                    control_points(i / N, i % N) = coeffs.row(i);
                }

                // The fine patches have their u direction pointing outwards
                // (first point is the center) while our system has the first
                // halfedge pointing outwards, so the first point is on the
                // edge. So we do a rotation here.
                control_points = control_points.rotate_cw();

                // Now that we have the new control net, update the face data
                // with the correct face and that net.
                auto& data =
                    face_data_vec.vector()[children_ordered[f_idx].idx()];
                data.face = children_ordered[f_idx];
                data.control_points = control_points;
            }
        }
        else
        {
            // === ORDINARY VERTICES ===
            // via deCasteljau

            // Get the face data and store it in a temporary dynamic 2d array.
            gsMatrix<gsVector<real_t>, Dynamic, Dynamic> control_net(
                face_data_vec.vector()[initial_face].control_points);

            // now control_net is a (n+1)*(n+1) matrix of control points (degree
            // n) Perform deCasteljau once to divide into two (n+1)*(n+1)
            // matrices of control points.
            auto const first_split = deCasteljau(control_net);

            // Perform deCasteljau again on both of them, to get 4 (n+1)*(n+1)
            // matrices of control points In between, we need to transpose so we
            // now divide in the other direction.
            auto top_split = deCasteljau(first_split[0].transpose());
            auto bot_split = deCasteljau(first_split[1].transpose());

            // Re-transpose
            top_split[0] = top_split[0].transpose().eval();
            top_split[1] = top_split[1].transpose().eval();
            bot_split[0] = bot_split[0].transpose().eval();
            bot_split[1] = bot_split[1].transpose().eval();

            // Rotate based on position 0, 1, 2, and 3 times.
            bot_split[1] = bot_split[1].rotate_cw();
            bot_split[0] = bot_split[0].rotate_cw().rotate_cw();
            top_split[0] = top_split[0].rotate_ccw();

            // Collate all these matrices in the correct order into an array.
            std::array<gsMatrix<gsVector<real_t>, Dynamic, Dynamic>, 4> arr = {
                top_split[0], top_split[1], bot_split[1], bot_split[0]};

            // Collate the faces into a correctly ordered array.
            auto children_ordered =
                order_faces(first_vertices[initial_face], child_faces);

            // Correct back references of face data and give them the correct
            // control points.
            for (size_t f = 0; f < 4; ++f)
            {
                auto& data = face_data_vec.vector()[children_ordered[f].idx()];
                data.face = children_ordered[f];
                data.control_points = arr[f];
            }
        }
    }
};
template <size_t N>
gsMatrix<real_t> gsFreeformSubdivision<N>::fit_ev_opt(gsMatrix<real_t> A,
                                                      gsMatrix<real_t> target,
                                                      size_t valence)
{
    // Just a direct least-squares solve with no regularisation.
    gsMatrix<real_t> solution = A.colPivHouseholderQr().solve(target);

    auto Apv = A.fullPivLu();
    Apv.setThreshold(1e-8);
    gsMatrix<> K = Apv.kernel();
    // gsWrite(K, "../filedata/" + m_options.getString("model_patch_path") +
    //                "Val" + std::to_string(valence) + "Kernel.xml");

    gsMatrix<> diff(2 * valence, 2 * valence + 1);
    diff.setZero();
    for (size_t i = 0; i < 2 * valence; i++)
    {
        diff(i, 0) = 1.0;
        diff(i, i + 1) = -1.0;
    }

    gsMatrix<> w = (diff * K).colPivHouseholderQr().solve(-diff * solution);

    gsInfo << "Initial fitting error: " << (A * solution - target).norm()
           << "\n";
    gsInfo << "Final fitting error: "
           << (A * (solution + K * w) - target).norm() << "\n";

    return solution + K * w;
}

template <size_t N>
gsMatrix<real_t> gsFreeformSubdivision<N>::fit_ev(gsMatrix<real_t> A,
                                                  gsMatrix<real_t> target,
                                                  size_t valence)
{
    gsMatrix<real_t> weights;
    if (m_options.getSwitch("weighted_fit"))
    {
        auto _readFile = gsReadFile<>(
            "freeform/val" + std::to_string(valence) + "_weights.xml", weights);
    }
    else
    {
        weights = gsVector<real_t>::Ones(A.rows());
    }
    auto W = weights.col(0).asDiagonal();
    // Solve the least-squares system A * solution = target
    // with Tikhonov regularization and additional constraints by
    // building the augmented system:
    //```
    //    [ A^T*W*A + lambda*I     C^T ] [ x ] = [ A^T*W*target ]
    //    [     C                  0   ] [ y ]   [     0        ]
    //```
    //  and solving
    size_t function_count(2 * valence + 1);
    // Regularization parameter
    // Set to 0 to allow a second smoothing to have 0 error.
    real_t lambda = 0.0;
    // Load the constraints into a matrix
    // Number of constraints should be the number of fitting functions
    // (`function_count = 2 * valence + 1`) minus the dimension of the
    // space (`v + 3`), i.e. `v-2`.
    gsMatrix<real_t> constraints;
    auto _readFile =
        gsReadFile<>(m_options.getString("model_patch_path") + "Val" +
                         std::to_string(valence) + "Constraints.xml",
                     constraints);
    size_t constraint_count = constraints.rows();

    // Build the matrix & target
    gsMatrix<real_t> augmented_A(function_count + constraint_count,
                                 function_count + constraint_count);
    gsMatrix<real_t> augmented_target(function_count + constraint_count, D);
    // Zero it first
    augmented_A.setZero();
    // Top left: Tikhonov system
    augmented_A.topLeftCorner(function_count, function_count) =
        A.transpose() * W * A +
        lambda * gsMatrix<real_t>::Identity(function_count, function_count);

    // Top right & Bottom right: Constraints
    augmented_A.topRightCorner(function_count, constraint_count) =
        constraints.transpose();
    augmented_A.bottomLeftCorner(constraint_count, function_count) =
        constraints;

    // Bottom right: Zero (Lagrange multiplier block)
    augmented_A.bottomRightCorner(constraint_count, constraint_count).setZero();

    // Top of augmented target: target via Tikhonov
    augmented_target.topRows(function_count) = A.transpose() * W * target;
    // Bottom of augmented target: zero, to ensure constraints are
    // fulfilled
    augmented_target.bottomRows(constraint_count).setZero();

    // Actually solve the system
    gsMatrix<real_t> augmented_solution =
        augmented_A.fullPivHouseholderQr().solve(augmented_target);

    // Extract the solution without the zeroes from the augmented system
    gsMatrix<real_t> solution = augmented_solution.topRows(function_count);

    gsInfo << "Total fitting error: " << (A * solution - target).norm() << "\n";
    gsInfo << "Constraint error: " << (constraints * solution).norm() << "\n";

    return solution;
}

template <size_t N>
void gsFreeformSubdivision<N>::smooth(
    size_t degree, std::vector<gsMatrix<real_t>>& ev_coefficients,
    std::vector<gsMatrix<real_t>>& ev_coefficients_outer)
{
    // TODO: Retrieve coefficients
}

template <size_t N> void gsFreeformSubdivision<N>::smooth(size_t degree)
{
    GISMO_ASSERT(degree == 1, "Only C1 smoothing supported.");
    // TODO: Maybe work with the same Assember/Solver as in fit_function and
    // laplace_beltrami here?

    gsMultiPatch<> multi_patch;
    gsMultiBasis<> multi_basis;
    gsMappedBasis<2> mapped_basis;
    this->c1_basis(multi_patch, multi_basis, mapped_basis);

    gsMatrix<real_t> coefficients;
    gsL2Projection<real_t>::project(multi_basis, mapped_basis, multi_patch,
                                    coefficients);

    // The L2 projection returns a flattened vector with one block per
    // coordinate direction; gsMappedSpline expects one row per mapped DoF.
    coefficients = coefficients.reshape(mapped_basis.size(), D);
    gsMappedSpline<2, real_t> solSpline(mapped_basis, coefficients);
    gsMultiPatch<> solField = solSpline.exportToPatches();

    // Write the solution back.
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        m_mesh->get_face_property<gsFreeformFaceData<N>>("bezier_points"));
    for (size_t k = 0; k < face_data_vec.vector().size(); ++k)
    {
        for (size_t i = 0; i < N; ++i)
        {

            for (size_t j = 0; j < N; ++j)
            {
                face_data_vec[k].control_points(i, j) =
                    solField.patch(k).coefs().row(i * N + j);
            }
        }
    }
}

template <size_t N>
void gsFreeformSubdivision<N>::c1_basis(gsMultiPatch<>& multi_patch,
                                        gsMultiBasis<>& multi_basis,
                                        gsMappedBasis<2>& mapped_basis)
{
    auto& mesh = *m_mesh;

    // normal multi patch and multi basis
    multi_patch = this->multipatch();
    multi_patch.computeTopology();
    multi_basis = gsMultiBasis<>(multi_patch);

    const size_t nPatches = multi_patch.nPatches();
    const index_t nLocalDofs = (index_t)(nPatches * N * N);
    GISMO_ASSERT(nLocalDofs == (index_t)multi_basis.totalSize(),
                 "Local DOF count mismatch");

    // ---- Helpers ------------------------------------------------

    // Position of halfedge h in face f's CCW halfedge list (0-based).
    auto get_k = [&](Face f, Halfedge h) -> index_t
    {
        index_t k = 0;
        for (Halfedge he : mesh.halfedges(f))
        {
            if (he == h)
                return k;
            ++k;
        }
        GISMO_ERROR("Halfedge not found in face");
    };

    // Flat local DOF index for patch p at oriented position (i,j).
    auto local_dof_rotated = [&](index_t p, index_t k, index_t i,
                                 index_t j) -> index_t
    {
        index_t rot_i = i;
        index_t rot_j = j;
        switch (k % 4)
        {
        case 0:
            rot_i = N - 1 - j;
            rot_j = i;
            break;
        case 1:
            break;
        case 2:
            rot_i = j;
            rot_j = N - 1 - i;
            break;
        default:
            rot_i = N - 1 - i;
            rot_j = N - 1 - j;
            break;
        }

        return p * (index_t)(N * N) + rot_i * (index_t)N + rot_j;
    };

    // ---- Working state ------------------------------------------

    // pre_mapper[ldof] is the sparse row of the mapper matrix
    // represented as a map from global DOF index to coefficient.
    std::vector<std::map<index_t, real_t>> pre_mapper(nLocalDofs);
    index_t global_dof_count = 0;

    // handled[ldof] becomes true once a mapping has been assigned.
    std::vector<bool> handled(nLocalDofs, false);

    // ================================================================
    // Phase 1 — EV vertices
    //
    // For each extraordinary vertex v (valence != 4, interior):
    //   • Assign 2*v+1 EV global DOFs.
    //   • For every control point in the 4v-patch ring that lies
    //     inside the EV support (all fitting functions non-zero),
    //     set its mapper row to the corresponding row of A_control
    //     (i.e. the fitting-function coefficient vector).
    //   • For inner control points outside the EV support ("outer
    //     inner" points) assign a free global DOF.
    //   Boundary points that are outside the EV support are left for
    //   the C1 constraint pass below.
    // ================================================================
    for (Vertex v : mesh.vertices())
    {
        if (is_ordinary(mesh, v))
            continue;

        const size_t valence = mesh.valence(v);
        const size_t patches_count = 4 * valence;
        const size_t function_count = 2 * valence + 1;

        // Collect the 4*valence (face, orienting-halfedge) pairs in
        // the same order as smooth(): valence inner patches first,
        // then 3*valence outer patches (three per inner halfedge).
        std::vector<Face> ev_faces;
        std::vector<Halfedge> ev_halfedges;
        ev_faces.reserve(patches_count);
        ev_halfedges.reserve(patches_count);

        for (Halfedge h : mesh.halfedges(v))
        {
            ev_faces.push_back(mesh.face(h));
            ev_halfedges.push_back(h);
        }

        for (Halfedge h_orig : mesh.halfedges(v))
        {
            Halfedge h = mesh.next_halfedge(h_orig);
            h = mesh.opposite_halfedge(h);
            h = mesh.next_halfedge(h);
            ev_faces.push_back(mesh.face(h));
            ev_halfedges.push_back(h);

            h = mesh.next_halfedge(h);
            h = mesh.next_halfedge(h);
            h = mesh.opposite_halfedge(h);
            ev_faces.push_back(mesh.face(h));
            ev_halfedges.push_back(h);

            h = mesh.prev_halfedge(h);
            h = mesh.opposite_halfedge(h);
            h = mesh.prev_halfedge(h);
            ev_faces.push_back(mesh.face(h));
            ev_halfedges.push_back(h);
        }

        // Load the 2v+1 model-patch fitting functions.
        std::vector<std::vector<std::unique_ptr<gsTensorBSpline<2, real_t>>>>
            fitting_functions;
        fitting_functions.reserve(function_count);
        for (size_t i = 0; i < function_count; ++i)
        {
            fitting_functions.emplace_back(
                gsFileData<real_t>(m_options.getString("model_patch_path") +
                                   "Val" + std::to_string(valence) + "Fct" +
                                   std::to_string(i) + ".xml")
                    .getAll<gsTensorBSpline<2, real_t>>());
        }

        // Reserve 2v+1 contiguous global DOFs for this EV.
        const index_t ev_dof_start = global_dof_count;
        global_dof_count += (index_t)function_count;

        for (size_t p = 0; p < patches_count; ++p)
        {
            const index_t patch_idx = ev_faces[p].idx();

            const index_t k = get_k(ev_faces[p], ev_halfedges[p]);

            for (size_t ux = 0; ux < N; ++ux)
            {
                for (size_t vx = 0; vx < N; ++vx)
                {
                    // Inside EV support = all fitting functions non-zero here.
                    const bool in_support = std::all_of(
                        fitting_functions.begin(), fitting_functions.end(),
                        [&](const auto& ff)
                        {
                            return ff[p]->coef((index_t)(ux * N + vx), 2) !=
                                   real_t(0);
                        });

                    // The values vx and ux are rotated to fit the EV, so we
                    // rotate them back to get the correct local DOF
                    const index_t ldof =
                        local_dof_rotated(patch_idx, k, vx, ux);

                    if (in_support && !handled[ldof])
                    {
                        // Map this point via the fitting-function coefficients
                        // (the A_control row) to the EV global DOFs.
                        for (size_t j = 0; j < function_count; ++j)
                        {
                            const real_t coeff = fitting_functions[j][p]->coef(
                                (index_t)(ux * N + vx), 2);
                            if (coeff != real_t(0))
                                pre_mapper[ldof][ev_dof_start + j] = coeff;
                        }
                        handled[ldof] = true;
                    }
                    else if (
                        // not in support of ev basis functions
                        !in_support &&
                        // but also not on the edge of a patch
                        ux > 0 && vx > 0 && ux < N - 1 && vx < N - 1 &&
                        // and not already handled
                        !handled[ldof])
                    {
                        // Outer inner point: free global DOF.
                        pre_mapper[ldof][global_dof_count] = real_t(1);
                        ++global_dof_count;
                        handled[ldof] = true;
                    }
                    // All other points are boundary non-support points and are
                    // handled in Phase 3/4.
                }
            }
        }
    } // end EV loop

    // ================================================================
    // Phase 2 — Ordinary inner points
    //
    // Every interior control point (1 <= i,j <= N-2) not yet handled
    // (i.e. not in any EV support or outer-inner position) gets its
    // own free global DOF.  These are the 9 independent DOFs per
    // ordinary patch.
    // ================================================================
    for (size_t p = 0; p < nPatches; ++p)
    {
        for (size_t i = 1; i < N - 1; ++i)
        {
            for (size_t j = 1; j < N - 1; ++j)
            {
                const index_t ldof = local_dof_rotated(p, 1, i, j);
                if (!handled[ldof])
                {
                    pre_mapper[ldof][global_dof_count] = real_t(1);
                    ++global_dof_count;
                    handled[ldof] = true;
                }
            }
        }
    }

    // ================================================================
    // Phase 3 — Interior edge points (boundary rows/columns, not
    //           corners, i.e. oriented (0,j) for j=1..N-2).
    //
    // For each interior mesh edge the C1 condition is
    //   edge_pt = 0.5 * inner_this + 0.5 * inner_adj
    // where inner_this is oriented (1,j) in the current patch and
    // inner_adj is oriented (1,N-1-j) in the adjacent patch.
    //
    // The boundary points of the EV block that are inside the EV
    // support are already handled (Phase 1); the !handled guard skips
    // them.  Edges on the mesh boundary get a free global DOF because
    // there is no adjacent patch to couple to.
    // ================================================================
    for (index_t p = 0; p < (index_t)nPatches; ++p)
    {
        Face f(p);
        index_t k = 0;
        for (Halfedge h : mesh.halfedges(f))
        {
            for (index_t j = 1; j < (index_t)N - 1; ++j)
            {
                const index_t ldof = local_dof_rotated(p, k, 0, j);

                if (handled[ldof])
                    continue;

                if (mesh.is_boundary(mesh.opposite_halfedge(h)))
                {
                    // Mesh boundary: unconstrained free DOF.
                    pre_mapper[ldof][global_dof_count] = real_t(1);
                    ++global_dof_count;
                }
                else
                {
                    // C1 constraint: average of the two adjacent inner pts.
                    const index_t inner_ldof = local_dof_rotated(p, k, 1, j);

                    const Halfedge h_opp = mesh.opposite_halfedge(h);
                    const Face f_adj = mesh.face(h_opp);
                    const index_t p_adj = f_adj.idx();
                    const index_t k_adj = get_k(f_adj, h_opp);

                    const index_t inner_adj_ldof =
                        local_dof_rotated(p_adj, k_adj, 1, N - 1 - j);

                    for (auto& kv1 : pre_mapper[inner_ldof])
                        pre_mapper[ldof][kv1.first] += real_t(0.5) * kv1.second;
                    for (auto& kv2 : pre_mapper[inner_adj_ldof])
                        pre_mapper[ldof][kv2.first] += real_t(0.5) * kv2.second;
                }
                handled[ldof] = true;
            }
            ++k;
        }
    }

    // ================================================================
    // Phase 4 — Corner points (oriented (0,0) of each halfedge).
    //
    // Two rules depending on whether v_corner is interior or on the
    // mesh boundary:
    //
    // Interior vertex (not mesh-boundary):
    //   Set to the average of (1,1) points of adjacent patches.
    //
    // Mesh-boundary vertex (including interior corners, i.e. vertices
    // where N >= 2 patches meet but a gap is present on the boundary):
    //   All patches sharing this vertex receive the SAME single free
    //   global DOF for their corner, enforcing C0 without involving
    //   any inner-point DOFs.  The DOF remains free so that
    //   gsDofMapper can eliminate it via markBoundary.
    //   Using a shared DOF (rather than the (1,1)-average) is
    //   essential: the average would only constrain the SUM of the
    //   inner DOFs, leaving their difference unconstrained and causing
    //   the second row of control points to escape when Dirichlet BCs
    //   fix the combined value.
    //
    // Points already handled by the EV support (Phase 1) are skipped.
    // ================================================================

    // One global DOF per distinct boundary vertex, shared across all
    // patches whose corner lands on that vertex.
    std::map<Vertex, index_t> boundary_corner_dof;

    for (Face f : mesh.faces())
    {
        index_t k = 0;
        for (Halfedge h : mesh.halfedges(f))
        {
            const index_t ldof = local_dof_rotated(f.idx(), k, 0, 0);

            if (!handled[ldof])
            {
                const Vertex v_corner = mesh.from_vertex(h);

                if (mesh.is_boundary(v_corner))
                {
                    // Boundary vertex: all patches at this vertex share
                    // one free global DOF so their corners are C0-coupled
                    // without constraining any inner-point DOFs.
                    auto res =
                        boundary_corner_dof.emplace(v_corner, global_dof_count);
                    if (res.second) // newly inserted → allocate the DOF
                        ++global_dof_count;
                    pre_mapper[ldof][res.first->second] = real_t(1);
                }
                else
                {

                    // Interior vertex: average of (1,1) inner points of all
                    // surrounding faces.
                    std::map<index_t, real_t> row;
                    index_t face_count = 0;

                    for (Halfedge h_out : mesh.halfedges(v_corner))
                    {
                        if (mesh.is_boundary(h_out))
                            continue;

                        const index_t p_adj = mesh.face(h_out).idx();
                        const index_t k_adj = get_k(mesh.face(h_out), h_out);

                        // The (1,1) inner point in the frame oriented at
                        // v_corner.
                        const index_t inner_ldof =
                            local_dof_rotated(p_adj, k_adj, 1, 1);

                        for (auto& kv : pre_mapper[inner_ldof])
                            row[kv.first] += kv.second;
                        ++face_count;
                    }

                    for (auto& kv : row)
                        pre_mapper[ldof][kv.first] =
                            kv.second / (real_t)face_count;
                }

                // this has now been handled
                handled[ldof] = true;
            }
            ++k;
        }
    }

    // Sanity check.
    for (index_t i = 0; i < nLocalDofs; ++i)
    {
        if (!handled[i])
        {
            gsWarn << "basis_data: local DOF " << i
                   << " was not handled — assigning free DOF.\n";
            pre_mapper[i][global_dof_count] = real_t(1);
            ++global_dof_count;
        }
    }

    gsInfo << "Total degrees of freedom: " << global_dof_count << "\n";

    // ================================================================
    // Assemble the sparse mapper matrix  (nLocalDofs x nGlobalDofs).
    // ================================================================
    gsSparseMatrix<> mapper_mat(nLocalDofs, global_dof_count);
    {
        gsSparseEntries<real_t> triplets;
        triplets.reserve(nLocalDofs * 4);
        for (index_t row = 0; row < nLocalDofs; ++row)
            for (auto& kv : pre_mapper[row])
                if (std::abs(kv.second) > real_t(1e-15))
                    triplets.add(row, kv.first, kv.second);
        mapper_mat.setFrom(triplets);
    }
    mapper_mat.makeCompressed();
    mapped_basis.init(multi_basis, mapper_mat);
}

template <size_t N>
void gsFreeformSubdivision<N>::laplace_beltrami(gsFunctionExpr<real_t> rhs)
{
    gsMultiPatch<> multi_patch;
    gsMultiBasis<> multi_basis;
    gsMappedBasis<2> mapped_basis;
    this->c1_basis(multi_patch, multi_basis, mapped_basis);

    // Set up the expression assembler.
    gsExprAssembler<> A(1, 1);
    // Must be called before any computation; sets the integration domain.
    A.setIntegrationElements(multi_basis);
    // The geometry map that parametrizes the surface over the patches $\Omega$.
    auto G = A.getMap(multi_patch);
    auto u = A.getSpace(mapped_basis);
    // Pull back the right hand side function $f$ to $\Omega$.
    // rhs is defined in the 2-D parametric domain, so evaluate it there
    // (not at physical points via G).
    auto ff = A.getCoeff(rhs);

    // Set homogeneous Dirichlet BCs on every boundary side.  Without
    // computeTopology() above, bBegin()==bEnd() and the system would have an
    // unconstrained constant null space on every patch, causing CG divergence.
    gsBoundaryConditions<> bc;
    bc.setGeoMap(multi_patch);
    for (auto it = multi_patch.bBegin(); it != multi_patch.bEnd(); ++it)
        bc.addCondition(*it, condition_type::dirichlet, nullptr);
    u.setup(bc, dirichlet::l2Projection, -1);
    A.initSystem();

    // Weak form stiffness matrix and RHS
    A.assemble(                                   // expressions:
        igrad(u, G) * igrad(u, G).tr() * meas(G), // stiffness
        u * ff * meas(G)                          // $\int f \cdot v d\Gamma$
    );

    // Solver and stuff
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute(A.matrix());
    gsMatrix<> solVector = solver.solve(A.rhs());
    auto solution = A.getSolution(u, solVector);

    // Extract the solution into a gsMappedSpline, then export to patches.
    gsMappedSpline<2, real_t> solSpline;
    solution.extract(solSpline);
    gsMultiPatch<> solField = solSpline.exportToPatches();

    // Write the solution back into the last entry.
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        m_mesh->get_face_property<gsFreeformFaceData<N>>("bezier_points"));
    for (size_t k = 0; k < face_data_vec.vector().size(); ++k)
    {
        for (size_t i = 0; i < N; ++i)
        {

            for (size_t j = 0; j < N; ++j)
            {
                face_data_vec[k].control_points(i, j).conservativeResize(D + 1);
                face_data_vec[k].control_points(i, j)(D) =
                    solField.patch(k).coefs()(i * N + j, 0);
            }
        }
        face_data_vec[k].D++;
    }

    D++;
}

template <size_t N>
void gsFreeformSubdivision<N>::fit_function(gsFunctionExpr<real_t> function)
{
    // TODO: Output coefficients again

    gsMultiPatch<> multi_patch;
    gsMultiBasis<> multi_basis;
    gsMappedBasis<2> mapped_basis;
    this->c1_basis(multi_patch, multi_basis, mapped_basis);

    // gsMatrix<real_t> coefficients;
    // real_t err = gsL2Projection<real_t>::project(
    //     mapped_basis,
    //     multi_basis,
    //     multi_patch,
    //     function,
    //     coefficients
    // );
    // gsMappedSpline<2, real_t> solSpline(mapped_basis, coefficients);
    // gsMultiPatch<> solField = solSpline.exportToPatches();

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(multi_basis);

    auto G = A.getMap(multi_patch);
    auto u = A.getSpace(mapped_basis, 1);
    auto ff = A.getCoeff(function, G);

    u.setup();
    A.initSystem();
    A.assemble(u * u.tr(), u * ff);

    gsSparseSolver<real_t>::LU solver;
    solver.compute(A.matrix());
    gsMatrix<> coefficients = solver.solve(A.rhs());
    auto solution = A.getSolution(u, coefficients);

    // Extract the solution into a gsMappedSpline, then export to patches.
    gsMappedSpline<2, real_t> solSpline;
    solution.extract(solSpline);
    gsMultiPatch<> solField = solSpline.exportToPatches();

    // Write the solution back into the last entry.
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        m_mesh->get_face_property<gsFreeformFaceData<N>>("bezier_points"));
    for (size_t k = 0; k < face_data_vec.vector().size(); ++k)
    {
        for (size_t i = 0; i < N; ++i)
        {

            for (size_t j = 0; j < N; ++j)
            {
                face_data_vec[k].control_points(i, j).conservativeResize(D + 1);
                face_data_vec[k].control_points(i, j)(D) =
                    solField.patch(k).coefs()(i * N + j, 0);
            }
        }
        face_data_vec[k].D++;
    }

    D++;
}

template <size_t N>
void gsFreeformSubdivision<N>::write_paraview_error(
    gsFunctionExpr<real_t> function, real_t max_error, std::string name,
    gsParaviewCollection* collection, size_t timestep)
{
    // Prepare data
    auto& mesh = *m_mesh;
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points"));

    // Number of faces needs to be counted for registration in the collection.
    size_t face_counter(0);

    for (auto f : mesh.faces())
    {
        // read out the control points and create a patch
        auto patch = face_data_vec[f.idx()].patch();

        // Sample at S*S points (S > N) for a fit that can capture errors.
        const size_t S = 2 * N;
        gsMatrix<> error_samples(1, S * S);
        gsMatrix<> params(2, S * S);

        for (size_t i = 0; i < S * S; ++i)
        {
            // Get the parameters at the sample point.
            params.col(i) = gsVector<real_t, 2>::vec(
                real_t(std::floor(i % S)) / real_t(S - 1),
                real_t(std::floor(i / S)) / real_t(S - 1));

            // Get the value of the patch
            gsVector<real_t> point = patch.eval(params.col(i));

            // Use the function to find the error
            error_samples(0, i) =
                std::abs(function.eval(point.topRows(D - 1))(0) - point(D - 1));
        }

        // rescale error to be in [0,1].
        error_samples /= max_error;

        // Fit a 1D B-spline patch to the sampled error field.
        // Use degree 2*(N-1) to better capture the error shape.
        gsKnotVector<> kv_err(0, 1, 0, S - 1);
        gsTensorBSplineBasis<2, real_t> fit_basis(kv_err, kv_err);
        gsFitting<> fitter(params, error_samples, fit_basis);
        fitter.compute(0.0);
        gsGeometry<>* result = fitter.result();

        // Now create a B-spline patch using the first min(3, D-1) coordinates
        // of the original patch. For D == 3 this gives a 2D geometry; for D > 3
        // this gives a full 3D geometry. gsWriteParaview supports both.
        size_t geo_dim = (D - 1 < 3) ? D - 1 : 3;
        gsTensorBSpline<2, real_t> patchGeo(
            patch.basis(), patch.coefs().leftCols(geo_dim).eval());

        // Combine the two objects into an error field.
        gsField<> field(patchGeo, *result);

        // Write this to a file.
        gsWriteParaview(field, name + "_" + std::to_string(face_counter), 1000);

        // If requested, also register it in the given collection at the given
        // time step.
        if (collection != nullptr)
        {
            std::string basename(name.substr(name.rfind('/') + 1));
            collection->addPart(basename + "_" + std::to_string(face_counter) +
                                    "0.vts",
                                timestep);
        }

        ++face_counter;
    }
}

template <size_t N>
gsVector<real_t, 2>
gsFreeformSubdivision<N>::error(gsFunctionExpr<real_t> function,
                                size_t samples_per_face)
{

    auto& mesh = *m_mesh;
    size_t spf = samples_per_face;
    gsProperty<gsFreeformFaceData<N>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points"));

    real_t error_linf(0.);
    real_t error_l2(0.);
    real_t error_count(0);

    for (auto f : mesh.faces())
    {
        // Look at the patch for this face.
        auto patch = face_data_vec[f.idx()].patch();

        // Now sample this patch
        for (size_t i = 0; i < spf * spf; ++i)
        {
            // Get the value of the patch
            gsVector<real_t> point = patch.eval(gsVector<real_t, 2>::vec(
                real_t(std::floor(i % spf)) / real_t(spf - 1),
                real_t(std::floor(i / spf)) / real_t(spf - 1)));

            // Compare its last coordinate to the value of the function applied
            // to the first coordinates and collate update the error
            // receptables.
            real_t err =
                abs(point(D - 1) - function.eval(point.topRows(D - 1))(0));
            error_linf = std::max(error_linf, err);
            error_l2 += err * err;
            ++error_count;
        }
    }

    gsVector<real_t> error =
        gsVector<real_t>::vec(error_linf, sqrt(error_l2 / real_t(error_count)));

    gsInfo << "Error L0: " << error(0) << ".\n";
    gsInfo << "Error L2: " << error(1) << ".\n";

    return error;
}

template <size_t N>
void gsFreeformSubdivision<N>::initialize_data_xml(std::string filepath)
{
    auto& mesh = *m_mesh;
    // Clear the mesh
    mesh = gsSurfMesh();

    // Load the TensorBSplinePatches
    std::vector<std::unique_ptr<gsTensorBSpline<2, real_t>>> patches =
        gsFileData<real_t>(filepath).getAll<gsTensorBSpline<2, real_t>>();

    // Initialize the property.
    auto bezier_points = mesh.add_face_property(std::string("bezier_points"),
                                                gsFreeformFaceData<N>(D));

    // Map from corner positions to vertex indices for detecting shared
    // vertices. We use a tolerance-based comparison for floating-point
    // coordinates.
    struct VecLess
    {
        bool operator()(const gsVector<>& a, const gsVector<>& b) const
        {
            GISMO_ASSERT(a.size() == b.size(), "Size mismatch in VecLess");
            return std::lexicographical_compare(a.data(), a.data() + a.size(),
                                                b.data(), b.data() + b.size());
        }
    };
    std::map<gsVector<>, gsSurfMesh::Vertex, VecLess> cornerMap;
    const real_t tolerance = 1e-10;

    auto findOrCreateVertex =
        [&](const gsMatrix<real_t>& point) -> gsSurfMesh::Vertex
    {
        // Round coordinates for map lookup, ensure there are at least 3
        // dimensions
        gsVector<> key(point.size());
        for (index_t i = 0; i < point.size(); ++i)
        {
            key[i] = std::round(point(i) / tolerance) * tolerance;
        }

        auto it = cornerMap.find(key);
        if (it != cornerMap.end())
        {
            return it->second;
        }
        else
        {
            // Mesh vertices are always 3D; pad with zeros if the file has
            // fewer than 3 dimensions.
            gsSurfMesh::Point pt(0, 0, 0);
            for (index_t k = 0; k < std::min(point.size(), (index_t)3); ++k)
                pt[k] = point(k);
            gsSurfMesh::Vertex v = mesh.add_vertex(pt);
            cornerMap[key] = v;
            return v;
        }
    };

    // Process each patch
    for (size_t patchIdx = 0; patchIdx < patches.size(); ++patchIdx)
    {
        const auto& patch = patches[patchIdx];

        // Get control points and dimensions
        const gsMatrix<real_t>& coefs = patch->coefs();

        // Extract corner control points (lexicographic indexing: i + j*n_u).
        // BSpline corners: (0,0), (N-1,0), (N-1,N-1), (0,N-1) in (u,v)
        // coordinates. Map to mesh vertices based on their physical positions.
        std::vector<gsSurfMesh::Vertex> corners(4);
        corners[0] =
            findOrCreateVertex(coefs.row(0 + 0 * N)); // BSpline (0,0) → v0
        corners[1] = findOrCreateVertex(
            coefs.row((N - 1) + 0 * N)); // BSpline (N-1,0) → v1
        corners[2] = findOrCreateVertex(
            coefs.row((N - 1) + (N - 1) * N)); // BSpline (N-1,N-1) → v2
        corners[3] = findOrCreateVertex(
            coefs.row(0 + (N - 1) * N)); // BSpline (0,N-1) → v3

        // Add face
        gsSurfMesh::Face f =
            mesh.add_quad(corners[0], corners[1], corners[2], corners[3]);

        // Build control points matrix for gsFreeformFaceData
        // The first halfedge goes from v3→v0 (corners[3]→corners[0])
        // Control point layout should be:
        //   faceControlPoints(0,0) near v0 = BSpline (0,0)
        //   faceControlPoints(0,N-1) near v1 = BSpline (N-1,0)
        //   faceControlPoints(N-1,0) near v3 = BSpline (0,N-1)
        //   faceControlPoints(N-1,N-1) near v2 = BSpline (N-1,N-1)
        // So the mapping is: faceControlPoints(i,j) = BSpline(j, i)
        gsMatrix<gsVector<real_t>, N, N> faceControlPoints;

        const size_t fileDim = (size_t)coefs.cols();
        const size_t copyDim = std::min(fileDim, D);
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                // Map face matrix (i,j) to B-spline (u,v) = (j,i).
                index_t linearIdx = j + i * N;
                // Copy the first copyDim coordinates; zero-fill any remainder
                // up to D if the file dimension is smaller than D.
                gsVector<real_t> point = gsVector<real_t>::Zero(D);
                for (size_t k = 0; k < copyDim; ++k)
                    point(k) = coefs(linearIdx, k);
                faceControlPoints(i, j) = point;
            }
        }

        // Create gsFreeformFaceData with control points and face back reference
        bezier_points[f] = gsFreeformFaceData<N>(faceControlPoints, f);
    }
}

template <size_t N>
void gsFreeformSubdivision<N>::initialize_data_off(std::string filepath)
{
    auto& mesh = *m_mesh;
    // Clear the mesh
    mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(filepath, mesh);
    // Initialize the property.
    mesh.add_face_property(std::string("bezier_points"),
                           gsFreeformFaceData<N>(D));

    // Get the data. It will be empty and non-valid at this point.
    gsProperty<gsFreeformFaceData<N>> patch_data =
        mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points");

    // Each patch is now initalized with basic face data.
    for (auto f : mesh.faces())
    {
        patch_data.vector()[f.idx()] = gsFreeformFaceData<N>(mesh, f);
    }
}

template <size_t N>
void gsFreeformSubdivision<N>::initialize_data(std::string filepath)
{
    std::string xml(".xml");
    std::string off(".off");
    // Check the filetype to be loaded.
    if (std::equal(filepath.begin() + filepath.size() - xml.size(),
                   filepath.end(), xml.begin()))
    {
        gsInfo << "Loading xml\n";
        initialize_data_xml(filepath);
    }
    else if (std::equal(filepath.begin() + filepath.size() - off.size(),
                        filepath.end(), off.begin()))
    {
        gsInfo << "Loading off\n";
        initialize_data_off(filepath);
    }
    else
    {
        gsWarn << "Unsupported Filetype! Mesh not initialized.\n";
    }
    m_mesh->garbage_collection();
}

template <size_t N>
void gsFreeformSubdivision<N>::initialize_data(std::string filepath, size_t D)
{
    this->D = D;
    initialize_data(filepath);
}

template <size_t N>
void gsFreeformSubdivision<N>::write_paraview(
    std::string name, gsParaviewCollection* collection,
    gsParaviewCollection* cnet_collection, size_t timestep, bool control_net)
{
    auto mp(multipatch());
    std::string basename(name.substr(name.rfind('/') + 1));

    gsWriteParaview(mp, name, 1000, false, control_net);

    if (collection == nullptr)
        return;
    // Register each patch's .vts file in the time series collection
    for (size_t j = 0; j < mp.nPatches(); ++j)
    {
        collection->addPart(basename + "_" + std::to_string(j) + ".vts",
                            timestep, "", j);
        if (control_net && cnet_collection != nullptr)
        {
            // Also register the control net in the collection
            cnet_collection->addPart(basename + "_" + std::to_string(j) +
                                         "_cnet.vtp",
                                     timestep, "", j);
        }
    }
}

template <size_t N> gsMultiPatch<> gsFreeformSubdivision<N>::multipatch()
{
    auto& mesh = *m_mesh;
    gsMultiPatch<> patch;

    // Get the vector containing all the face data.
    gsProperty<gsFreeformFaceData<N>> face_data_vec =
        mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points");

    // For each face, convert its control net to a patch and add it to the
    // multipatch. Order doesn't matter.
    for (auto face : mesh.faces())
    {
        patch.addPatch(face_data_vec.vector()[face.idx()].patch());
    }

    return patch;
}

template <size_t N>
gsSubdivisionScheme::gsSubdivisionMeshValidity
gsFreeformSubdivision<N>::check_mesh()
{
    auto& mesh = *m_mesh;
    for (Face f : mesh.faces())
    {
        size_t count(0);
        for ([[maybe_unused]] Vertex v : mesh.vertices(f))
        {
            ++count;
        }
        if (count != 4)
        {
            gsWarn << "This mesh has at least one non-quadrangular face "
                      "(Vertex count "
                   << count << ").";
            return gsSubdivisionScheme::gsSubdivisionMeshValidity::INVALID;
        }
    }
    return gsSubdivisionScheme::gsSubdivisionMeshValidity::UNDETERMINED;
}

template class gsFreeformSubdivision<5>;
template class gsFreeformFaceData<5>;
template class gsFreeformSubdivision<6>;
template class gsFreeformFaceData<6>;
template class gsFreeformSubdivision<9>;
template class gsFreeformFaceData<9>;

} // namespace gismo
