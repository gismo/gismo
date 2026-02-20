/** @file gsFreeformSubdivision.cpp

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#include <cmath>
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

template <size_t N, size_t D>
gsFreeformFaceData<N, D>::gsFreeformFaceData(const gsSurfMesh& mesh,
                                             gsSurfMesh::Face face)
    : control_points(), face(face)
{
    // Create a vector of the 4 corners.
    std::vector<gismo::gsVector3d<>> corners;
    corners.reserve(4);
    for (auto const& v : mesh.vertices(face))
    {
        corners.emplace_back(mesh.position(v));
    }
    assert(points.size() == 4);

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

template <size_t N, size_t D>
gsMatrix<gsVector<real_t, D>*>
gsFreeformFaceData<N, D>::control_points_oriented(gsSurfMesh& mesh,
                                                  Halfedge hedge)
{

    gsMatrix<gsVector<real_t, D>*> result;
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

template <size_t N, size_t D>
const gismo::gsTensorBSpline<2, real_t> gsFreeformFaceData<N, D>::patch() const
{
    // Create a spline basis for a normal bezier patch.
    gsKnotVector<> kv1(0, 1, 0, N);
    gsKnotVector<> kv2(0, 1, 0, N);
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
    // Create a coefficient matrix out of the control points.
    // Technically, you could use just [i] and one loop here since the elements
    // of a matrix are layed out row-wise, but this might be clearer to read.
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

template <size_t N, size_t D>
std::array<gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>, 2>
gsFreeformSubdivision<N, D>::deCasteljau(
    const gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>& control_net)
{

    // Create the 3d data vector and ensure it has the right size.
    std::vector<gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>> points;
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
    gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic> result1;
    gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic> result2;
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

template <size_t N, size_t D>
gismo::gsTensorBSpline<2, real_t>
gsFreeformSubdivision<N, D>::load_model_patch(int valence, std::string subtype)
{
    // build the filepath
    auto path = "freeformSubdivision/control_net_d" + std::to_string(N) + "_v" +
                std::to_string(valence) + "_" + subtype + ".xml";

    // load the file containing the correct matrix of control points
    gsFileData<real_t> filedata(path);
    auto mat = filedata.getFirst<gsMatrix<real_t>>();

    // create a basis vector etc.
    gsKnotVector<> kv1(0, 1, 0, N);
    gsKnotVector<> kv2(0, 1, 0, N);
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    return gsTensorBSpline<2>(basis, *mat);
}

template <size_t N, size_t D>
std::array<gsSurfMesh::Face, 4>
gsFreeformSubdivision<N, D>::order_faces(Vertex first_vertex,
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

template <size_t N, size_t D> void gsFreeformSubdivision<N, D>::orient_faces()
{
    // Get data
    auto mesh = *m_mesh;
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points"));

    // As a pre-step, we make sure that for each face, if it has an
    // extraordinary vertex, that vertex is the top left one
    for (Face f : mesh.faces())
    {
        // first, check if this face has an EV at all
        bool has_ev(false);
        for (Vertex v : mesh.vertices(f))
        {
            has_ev |= !mesh.is_ordinary(v);
        }

        // if it does, rotate its first halfedge around until it points to the
        // EV
        if (has_ev)
        {
            while (mesh.is_ordinary(mesh.to_vertex(mesh.halfedge(f))))
            {
                // rotate the edge
                mesh.set_halfedge(f, mesh.next_halfedge(mesh.halfedge(f)));
                // also rotate the control points
                gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic> old_points(
                    face_data_vec.vector()[f.idx()].control_points);
                face_data_vec.vector()[f.idx()].control_points =
                    old_points.rotate_cw();
            }
        }
    }
}

template <size_t N, size_t D> void gsFreeformSubdivision<N, D>::subdivide()
{
    auto& mesh = *m_mesh;
    // First, make sure all faces are correctly oriented with the EV as their
    // first vertex.
    orient_faces();

    // Remember the first vertex of each face (this is where the control nets of
    // each face data are oriented on).
    std::vector<Vertex> first_vertices;
    for (Face f : mesh.faces())
    {
        // remember the first vertex
        first_vertices.emplace_back(mesh.to_vertex(mesh.halfedge(f)));
    }
    // Remember the number of faces before the subdivision.
    size_t n(first_vertices.size());

    // Do the quad split of the (abstract) base faces.
    mesh.quad_split();

    // Get face data
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points"));

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
        if (!mesh.is_ordinary(fv))
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
                    // the tolerance is squared, so this is a tolerance of 1e-4.
                    coarse_model.closestPointTo(point, closest_point, 1e-2,
                                                true);

                    // Sample the old control net.
                    samples.col(i) = coarse_patch.eval(closest_point);
                }

                // Fit a NxN Bezier patch to these samples
                gsKnotVector<> kv1(0, 1, 0, N);
                gsKnotVector<> kv2(0, 1, 0, N);
                gsTensorBSplineBasis<2> basis(kv1, kv2);
                gsFitting<> fitter(params, samples, basis);
                fitter.compute(0.0);
                gsGeometry<>* result = fitter.result();

                // Extract control points
                const gsMatrix<>& coeffs = result->coefs();

                // Reshape the control points
                gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic> control_points;
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
            gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic> control_net(
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
            std::array<gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>, 4> arr =
                {top_split[0], top_split[1], bot_split[1], bot_split[0]};

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

template <size_t N, size_t D>
void gsFreeformSubdivision<N, D>::smooth(size_t degree)
{
    auto& mesh = *m_mesh;
    // Ensure we have a high enough degree
    if (degree + 1 > N / 2)
    {
        gsWarn
            << "Degree of Bezier control net to small for this smoothness.\n";
        return;
    }

    // Currently, only C^1 is supported.
    if (degree != 1)
    {
        gsWarn << "Currently, only C^1 smoothing is supported.\n";
        return;
    }

    // Cache face data
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points"));

    // First, correct each vertex
    for (const Vertex& v : mesh.vertices())
    {
        // Check for EVs on any face adjacent to this vertex.
        bool has_ev(false);
        for (Face f : mesh.faces(v))
        {
            for (Vertex v_face : mesh.vertices(f))
            {
                has_ev |= !mesh.is_ordinary(v_face);
            }
        }

        if (!mesh.is_ordinary(v))
        {
            size_t valence = mesh.valence(v);
            size_t patches = 4 * valence;

            // Collect all the control points.
            // The `4 * valence` faces will be ordered as follows, with the `o`
            // indicating its (0,0) control point and the arrow its primary
            // u-direction towards the (0,N-1) control point:
            // ```
            //    +-----+  +-----+  +-----+  +-----+
            //    |    ^|  |    ^|  |     |  |     |
            //    |  5 ||  |  4 ||  |  15 |  |  14 |
            //    |    o|  |    o|  |o--> |  |o--> |
            //    +-----+  +-----+  +-----+  +-----+
            //    +-----+  +-----+  +-----+  +-----+
            //    |    ^|  |    ^|  |     |  |     |
            //    |  6 ||  |  0 ||  |  3  |  |  13 |
            //    |    o|  |    o|  |o--> |  |o--> |
            //    +-----+  +-----+  +-----+  +-----+
            //    +-----+  +-----+  +-----+  +-----+
            //    | <--o|  | <--o|  |o    |  |o    |
            //    |  7  |  |  1  |  || 2  |  || 12 |
            //    |     |  |     |  |v    |  |v    |
            //    +-----+  +-----+  +-----+  +-----+
            //    +-----+  +-----+  +-----+  +-----+
            //    | <--o|  | <--o|  |o    |  |o    |
            //    |  8  |  |  9  |  || 10 |  || 11 |
            //    |     |  |     |  |v    |  |v    |
            //    +-----+  +-----+  +-----+  +-----+
            //
            // ```
            // TODO: And then we transpose
            std::vector<gsMatrix<gsVector<real_t, D>*>> control_nets;

            // First, collect the inner nets.
            for (Halfedge h : mesh.halfedges(v))
            {
                gsMatrix<gsVector<real_t, D>*> cp =
                    face_data_vec[mesh.face(h).idx()]
                        .control_points_oriented(mesh, h)
                        .transpose();
                control_nets.emplace_back(cp);
            }

            // Then, go through the halfedges again to get the outer nets.
            for (Halfedge h : mesh.halfedges(v))
            {
                // Follow halfedges to the next patch
                h = mesh.next_halfedge(h);
                h = mesh.opposite_halfedge(h);
                h = mesh.next_halfedge(h);
                // Get its oriented control net
                gsMatrix<gsVector<real_t, D>*> cp1 =
                    face_data_vec[mesh.face(h).idx()]
                        .control_points_oriented(mesh, h)
                        .transpose();
                // Place it in the matrix.
                control_nets.emplace_back(cp1);

                // Repeat two more times.
                h = mesh.next_halfedge(h);
                h = mesh.next_halfedge(h);
                h = mesh.opposite_halfedge(h);
                gsMatrix<gsVector<real_t, D>*> cp2 =
                    face_data_vec[mesh.face(h).idx()]
                        .control_points_oriented(mesh, h)
                        .transpose();
                control_nets.emplace_back(cp2);

                h = mesh.prev_halfedge(h);
                h = mesh.opposite_halfedge(h);
                h = mesh.prev_halfedge(h);
                gsMatrix<gsVector<real_t, D>*> cp3 =
                    face_data_vec[mesh.face(h).idx()]
                        .control_points_oriented(mesh, h)
                        .transpose();
                control_nets.emplace_back(cp3);
            }

            // Now, use this data to construct patches that we can sample.
            std::vector<gsTensorBSpline<2, real_t>> target_patches;
            for (size_t p = 0; p < patches; ++p)
            {
                // Knot vectors for a tensor product bezier patch
                gsKnotVector<> kv1(0, 1, 0, N);
                gsKnotVector<> kv2(0, 1, 0, N);
                gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
                // Create a coefficient matrix out of the control points.
                gsMatrix<> coeffs(N * N, D);
                for (size_t i = 0; i < N; ++i)
                {
                    for (size_t j = 0; j < N; ++j)
                    {
                        size_t total_index = i * N + j;
                        coeffs.row(total_index) = *control_nets[p](i, j);
                    }
                }

                // Create the patch and place it in the vector.
                target_patches.emplace_back(basis, coeffs);
            }

            // This is the number of total points that need to be moved.
            size_t point_count(valence * (N * N + 2 * N + 2 * 2 + N * 2));
            // This is the number of total points we are fitting to, excluding
            // doubles.
            // size_t sample_count(1 + valence * (N-1) * (N-1) + valence * 2 *
            // (N-1) + valence * 2 + valence * (N-1));
            size_t sample_count(1 + 3 * valence);

            // The number of fitting functions.
            // Note these are linearly dependent and actually the dimension of
            // this space is only valence + 3.
            size_t function_count(2 * valence + 1);

            // We now load the patches of these fitting functions.
            std::vector<
                std::vector<std::unique_ptr<gsTensorBSpline<2, real_t>>>>
                fitting_functions;
            fitting_functions.reserve(2 * valence + 1);

            for (size_t i = 1; i <= function_count; ++i)
            {
                // Construct the filepath for the ith basis function and
                // load it with gismo utilities
                fitting_functions.emplace_back(
                    gsFileData<real_t>(
                        "freeformSubdivision/fitting_functions/Val" +
                        std::to_string(valence) + "Fct" + std::to_string(i) +
                        ".xml")
                        .getAll<gsTensorBSpline<2, real_t>>());
            }

            // For the least-squares fit we now build:
            // - The matrix `A` of dimension `sample_count x fit_functions`, in
            // which the entry A(i,j) is the value of the j-th fitting function
            // at the ith sample point. . The matrix `A_cp` of dimension
            // `sample_count x fit_functions`, in which the entry A(i,j) is the
            // ith control point of meshes from the jth fitting function.
            // - The matrix `target` of dimension `sample_count x D`, in which
            // the row b(i,-) is the value of the target patch at fitting point
            // i
            //
            // The control points and sampling points are arranged in the same
            // way.
            gsMatrix<real_t> A(sample_count, function_count);
            gsMatrix<real_t> A_cp(point_count, function_count);
            gsMatrix<real_t> target(sample_count, D);

            {
                // the sampling index - this is incremented whenever we sample a
                // point
                size_t i_s(0);
                // The point index
                size_t i_p(0);
                // Iterating over patches
                for (size_t p = 0; p < patches; ++p)
                {
                    // Iterating over the N * N sample points on that patch by
                    // index, u dimension.
                    for (size_t ux = 0; ux < N; ++ux)
                    {
                        // Iterating over the N * N sample points on that patch
                        // by index, v dimension.
                        for (size_t vx = 0; vx < N; ++vx)
                        {
                            // control point of fitting function at this point
                            // If this is zero, we want to skip it since we only
                            // fit within the support of the fitting functions.
                            // TODO: Check all, not just 1.
                            real_t cp =
                                fitting_functions[0][p]->coef(vx * N + ux, 2);

                            if (cp == 0.0)
                                continue;

                            // First, log the z-coordinate of the control points
                            // into A_cp.
                            for (size_t j = 0; j < function_count; ++j)
                            {
                                A_cp(i_p, j) = fitting_functions[j][p]->coef(
                                    vx * N + ux, 2);
                            }
                            i_p++;

                            // Then, check if this is a sample point
                            if (!((p == 0 && vx == 0 &&
                                   ux == 0) || //
                                               // (ux > 0 && vx < N-1)
                                  (p < valence && ux == N - 1 &&
                                   (vx == 0 || vx == N - 1)) || //
                                  (p < valence && ux == 2 && vx == 2)))
                            {
                                continue;
                            }

                            // Transform the indices into real coordinates for
                            // the sampling points.
                            real_t u = real_t(ux) / real_t(N - 1);
                            real_t v = real_t(vx) / real_t(N - 1);

                            // For each function, save its z-value into A as
                            // described above
                            for (size_t j = 0; j < function_count; ++j)
                            {
                                // sample
                                gsVector<real_t, D> val =
                                    fitting_functions[j][p]->eval(
                                        gsVector<real_t>::vec(u, v));
                                A(i_s, j) = val.z();
                            }

                            // Lastly, save the sample value of the target patch
                            gsVector<real_t, D> val = target_patches[p].eval(
                                gsVector<real_t>::vec(u, v));
                            target.row(i_s) = val;

                            i_s++;
                        }
                    }
                }
                // gsInfo << "Total samples: " << i_s << "\n";
                // gsInfo << "Total points: " << i_p << "\n";
            }

            // Solve the least-squares system A * solution = target
            // with Tikhonov regularization and additional constraints by
            // building the augmented system:
            //```
            //`   /                         \  /  \    /           \
            //`  | A^T*A + lambda*I     C^T | | x | = | A^T*target |
            //`  |     C                0   | | y |   |     0      |
            //`  \                         /  \  /    \            /
            //```
            //  and solving

            // Regularization parameter
            real_t lambda = 1e-8;
            // Load the constraints into a matrix
            // Number of constraints should be the number of fitting functions
            // (`function_count = 2 * valence + 1`) minus the dimension of the
            // space (`v + 3`), i.e. `v-2`.
            gsMatrix<real_t> constraints;
            auto _readFile =
                gsReadFile<>("freeformSubdivision/fitting_functions/Val" +
                                 std::to_string(valence) + "Constraints.xml",
                             constraints);
            size_t constraint_count = constraints.rows();

            // Build the matrix & target
            gsMatrix<real_t> augmented_A(function_count + constraint_count,
                                         function_count + constraint_count);
            gsMatrix<real_t> augmented_target(function_count + constraint_count,
                                              target.cols());
            // Zero it first
            augmented_A.setZero();
            // Top left: Tikhonov system
            augmented_A.topLeftCorner(function_count, function_count) =
                A.transpose() * A +
                lambda *
                    gsMatrix<real_t>::Identity(function_count, function_count);

            // Top right & Bottom right: Constraints
            augmented_A.topRightCorner(function_count, constraint_count) =
                constraints.transpose();
            augmented_A.bottomLeftCorner(constraint_count, function_count) =
                constraints;

            // Bottom left: Zero
            augmented_A.bottomRightCorner(constraint_count, constraint_count).setZero();

            // Top of augmented target: target via Tikhonov
            augmented_target.topRows(function_count) = A.transpose() * target;
            // Bottom of augmented target: zero, to ensure constraints are fulfilled
            augmented_target.bottomRows(constraint_count).setZero();

            // Acutally solve the system
            gsMatrix<real_t> augmented_solution =
                augmented_A.colPivHouseholderQr().solve(augmented_target);

            // Extract the solution without the zeroes from the augmented system
            gsMatrix<real_t> solution = augmented_solution.topRows(function_count);

            gsInfo << "Total fiting error: " << (A * solution - target).norm()
                   << "\n";

            // Print the coefficients
            gsInfo << "Fitting coefficients: \n";
            for (int i = 0; i < solution.rows(); ++i)
            {
                gsInfo << solution.row(i) << "\n";
            }

            // Now, the coefficients in `solution` give a linear combination of
            // the fitting functions that approximates the original target
            // patches. By multiplying with A_cp, we get the same linear
            // combination of the control points of the fitting functions,
            // resulting in a collection of control points that we can write
            // back into the control nets.
            auto new_values = A_cp * solution;

            // We need to make sure to use the same ordering as we did when
            // sampling.
            {
                size_t i(0);
                for (size_t p = 0; p < patches; ++p)
                {
                    for (size_t ux = 0; ux < N; ++ux)
                    {
                        for (size_t vx = 0; vx < N; ++vx)
                        {
                            real_t cp =
                                fitting_functions[0][p]->coef(vx * N + ux, 2);

                            // Makre sure to skip all points that we skipped
                            // when collecting A_cp, those that lie outside of
                            // the support of the fitting functions.
                            if (cp == 0.0)
                                continue;

                            *control_nets[p](vx, ux) = new_values.row(i);

                            i++;
                        }
                    }
                }
            }
        }
        else
        {

            if (has_ev)
                continue;

            // first, collect all the control points around this vertex. They
            // will be arrayed like this:
            // ```
            //      ...   ...     ...   ...
            // ... (1,1) (1,0)   (0,1) (1,1) ...
            // ... (0,1) (0,0)   (0,0) (1,0) ...
            // ...             V
            // ... (1,0) (0,0)   (0,0) (0,1) ...
            // ... (1,1) (0,1)   (1,0) (1,1) ...
            //      ...   ...     ...   ...
            // ```
            // The points directly neighboring another across different patches
            // should be equal (e.g. top left (1,0) and bottom left (0,1)) as
            // the C0 condition.
            std::vector<gsMatrix<gsVector<real_t, D>*>> control_points_faces;
            size_t insert_index(0);
            for (Halfedge h : mesh.halfedges(v))
            {
                // we do skip non-existent faces
                if (mesh.is_boundary(h))
                {
                    // When a boundary is found, reset the insert-index.
                    // This ensures that, if we are on the boundary we still
                    // have a continuous series of faces.
                    insert_index = 0;
                    continue;
                }
                // Insert the face into our list at the appropriate position.
                control_points_faces.insert(
                    control_points_faces.begin() + insert_index++,
                    face_data_vec[mesh.face(h).idx()].control_points_oriented(
                        mesh, h));
            }
            // Now, each different point needs a row in the matrix. We have 4
            // base points with one patch in the top left. Each additional patch
            // adds another 2 points. Technically, the last patch (if connected)
            // adds only 1, but we can simply overdetermine the system a bit
            // since we use least squares anyways.
            size_t rows(2 + 2 * control_points_faces.size());

            // Now create a matrix that represents these C1 equations, i.e. each
            // point on a boundary is colinear with the ones on either side.
            // This results one row for each point, writing it as a linear
            // combination of the 4 free points in the top left patch. The
            // indices of the points above in the matrix are as follows:
            // ```
            // 2 3 3 4
            // 1 0 0 5
            // 9 0 0 5
            // 8 7 7 6
            // ```
            // where 9 has to be equal to 1, but this will happen automatically.
            auto matrix = gsMatrix<real_t>(rows, 4);
            // The first four equations just say that the first four points are
            // equal to themselves.
            matrix.row(0) << 1., 0., 0., 0.;
            matrix.row(1) << 0., 1., 0., 0.;
            matrix.row(2) << 0., 0., 1., 0.;
            matrix.row(3) << 0., 0., 0., 1.;

            // For each new face, its first vertex (with an even index e.g. 4 or
            // 8 above) depends on the two previous vertices and its second
            // vertex (with an odd index, e.g. 5 or 7 above) depends on the
            // center 0 and the second vertex of the preprevious face.
            for (size_t i = 4; i < rows; ++i)
            {
                if (i % 2 == 0)
                    matrix.row(i) = 2. * matrix.row(i - 1) - matrix.row(i - 2);
                else
                    matrix.row(i) = 2. * matrix.row(0) - matrix.row(i - 4);
            }

            // Now, for each of these 9 points, we want to find its desired
            // value by looking at the respective value of the old (non-smooth)
            // control net.
            // TODO: Does this overvalue the 0/9 point?
            gsMatrix<real_t> target_matrix(rows, D);
            target_matrix.setZero();
            target_matrix.row(0) = control_points_faces[0](0, 0)->transpose();
            target_matrix.row(1) = control_points_faces[0](0, 1)->transpose();

            for (size_t i = 0; i < control_points_faces.size(); ++i)
            {
                target_matrix.row(2 * i + 2) =
                    control_points_faces[i](1, 1)->transpose();
                target_matrix.row(2 * i + 3) =
                    control_points_faces[i](1, 0)->transpose();
            }

            // Now do a least squares fit.
            // I.e. we are searching for values for the 4 free points
            // (transformed into 9 points via `matrix` that are thus C1 smooth)
            // such that the squared distance of all 9 points to their previous
            // values (given in `target_matrix`) is minimal.
            gsMatrix<real_t> solution =
                matrix.colPivHouseholderQr().solve(target_matrix);

            auto new_points = matrix * solution;

            // Now re-assign the correct solution rows back to the points.
            // The linear combinations are the same as above, so the result will
            // be C1 smooth.
            for (size_t i = 0; i < control_points_faces.size(); ++i)
            {
                *(control_points_faces[i](0, 0)) = new_points.row(0);
                *(control_points_faces[i](0, 1)) = new_points.row(2 * i + 1);
                *(control_points_faces[i](1, 1)) = new_points.row(2 * i + 2);
                *(control_points_faces[i](1, 0)) = new_points.row(2 * i + 3);
            }
        }
    }
    // End of edge correction

    // Now, correct the remaining part of each edge that isn't surrounding a
    // vertex.
    for (const Edge& e : mesh.edges())
    {
        Halfedge halfedge0 = mesh.halfedge(e, 0);
        Halfedge halfedge1 = mesh.halfedge(e, 1);

        // If we are on the boundary of the mesh, nothing needs to be done.
        if (mesh.is_boundary(halfedge0) || mesh.is_boundary(halfedge1))
            continue;

        // Get the faces. These must be valid by the checks above.
        Face face0 = mesh.face(halfedge0);
        Face face1 = mesh.face(halfedge1);

        // Skip this edge if either of these faces has an EV.
        bool has_ev(false);
        for (Vertex v : mesh.vertices(face0))
        {
            has_ev |= !mesh.is_ordinary(v);
        }
        for (Vertex v : mesh.vertices(face1))
        {
            has_ev |= !mesh.is_ordinary(v);
        }

        if (has_ev)
            continue;

        // Get the control points. They will be arrayed like this:
        // ```
        // (1,0) (1,1) ... (1,i)... (1,N-2) (1,N-1)
        // (0,0) (0,1) ... (0,i)... (0,N-2) (0,N-1)
        // -------------- halfedge0 -------------->
        // <------------- halfedge1 ---------------
        // (0,N-1) (0,N-2) ... (0,i) ... (0,1) (0,0)
        // (1,N-1) (1,N-2) ... (1,i) ... (1,1) (1,0)
        //
        // ```
        // Note that in each column, the middle points should be equal by C0
        // conditions.
        auto cp0 =
            face_data_vec[face0.idx()].control_points_oriented(mesh, halfedge0);
        auto cp1 =
            face_data_vec[face1.idx()].control_points_oriented(mesh, halfedge1);

        // Now, for each column above, we need to ensure the three points are
        // colinear.
        for (size_t i = 2; i < N - 2; ++i)
        {
            // Extract the column.
            gsVector<real_t, D>* i0 = cp0(1, i);
            gsVector<real_t, D>* m0 = cp0(0, i);
            gsVector<real_t, D>* m1 = cp1(0, N - 1 - i);
            gsVector<real_t, D>* i1 = cp1(1, N - 1 - i);

            // Ensure the mesh was at least C0 and warn the user if it is not.
            if ((*m0 - *m1).squaredNorm() > 1e-5)
            {
                gsWarn << face0 << " and " << face1 << " along edges "
                       << halfedge0 << " and " << halfedge1 << " have values ("
                       << m0->transpose() << ") and (" << m1->transpose()
                       << ").\n";
                continue;
            }

            // Now, of these three points, two will be 'free' and the last will
            // be determined. The following matrix generates all three points
            // from the free ones.
            gsMatrix<real_t, D, 2> matrix;
            matrix.row(0) << 1., 0.;
            matrix.row(1) << 0.5, 0.5;
            matrix.row(2) << 0., 1.;

            // As a target matrix, use the old values of the points.
            gsMatrix<real_t> target_matrix(3, D);
            target_matrix.setZero();
            target_matrix.row(0) = i0->transpose();
            target_matrix.row(1) = m0->transpose();
            target_matrix.row(2) = i1->transpose();

            // We now do a least squares fit, i.e. search for two free points
            // such that the three generated ( and thus colinear) points wil lbe
            // as close as possible to the original three.
            gsMatrix<real_t> solution =
                matrix.colPivHouseholderQr().solve(target_matrix);

            // And finally re-assign all four values.
            *i0 = solution.row(0);
            *m0 = (solution.row(0) + solution.row(1)) * 0.5;
            *m1 = (solution.row(0) + solution.row(1)) * 0.5;
            *i1 = solution.row(1);
        }
    }
    // End of vertex correction
}

template <size_t N, size_t D>
void gsFreeformSubdivision<N, D>::initialize_data_off(std::string filepath)
{
    auto& mesh = *m_mesh;
    // Clear the mesh
    mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(filepath, mesh);
    // Initialize the property.
    mesh.add_face_property(std::string("bezier_points"),
                           gsFreeformFaceData<N, D>());

    // Get the data. It will be empty and non-valid at this point.
    gsProperty<gsFreeformFaceData<N, D>> patch_data =
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points");

    // Each patch is now initalized with basic face data.
    for (auto f : mesh.faces())
    {
        patch_data.vector()[f.idx()] = gsFreeformFaceData<N, D>(mesh, f);
    }
}

template <size_t N, size_t D>
void gsFreeformSubdivision<N, D>::initialize_data_xml(std::string filepath)
{
    auto& mesh = *m_mesh;
    // Clear the mesh
    mesh = gsSurfMesh();

    // Load the TensorBSplinePatches
    std::vector<std::unique_ptr<gsTensorBSpline<2, real_t>>> patches =
        gsFileData<real_t>(filepath).getAll<gsTensorBSpline<2, real_t>>();

    // Initialize the property.
    auto bezier_points = mesh.add_face_property(std::string("bezier_points"),
                                                gsFreeformFaceData<N, D>());

    // Map from corner positions to vertex indices for detecting shared vertices
    // We use a tolerance-based comparison for floating point coordinates
    std::map<std::array<real_t, D>, gsSurfMesh::Vertex> cornerMap;
    const real_t tolerance = 1e-10;

    auto findOrCreateVertex =
        [&](const gsMatrix<real_t>& point) -> gsSurfMesh::Vertex
    {
        // Round coordinates for map lookup
        std::array<real_t, D> key;
        for (size_t i = 0; i < D; ++i)
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
            gsSurfMesh::Vertex v = mesh.add_vertex(gsSurfMesh::Point(point));
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

        // Extract corner control points (lexicographic indexing: i + j*n_u)
        // BSpline corners: (0,0), (N-1,0), (N-1,N-1), (0,N-1) in (u,v)
        // coordinates Map to mesh vertices based on their physical positions
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
        gsMatrix<gsVector<real_t, D>, N, N> faceControlPoints;

        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                // Map face matrix (i,j) to BSpline (u,v) = (j,i)
                index_t linearIdx = j + i * N;
                // Extract the row as a 3D point and assign it
                auto point = coefs.row(linearIdx);
                faceControlPoints(i, j) = point;
            }
        }

        // Create gsFreeformFaceData with control points and face back reference
        bezier_points[f] = gsFreeformFaceData<N, D>(faceControlPoints, f);
    }
}

template <size_t N, size_t D>
gsMultiPatch<> gsFreeformSubdivision<N, D>::multipatch()
{
    auto& mesh = *m_mesh;
    gsMultiPatch<> patch;

    // Get the vector containing all the face data.
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec =
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points");

    // For each face, convert its control net to a patch and add it to the
    // multipatch. Order doesn't matter.
    for (auto face : mesh.faces())
    {
        patch.addPatch(face_data_vec.vector()[face.idx()].patch());
    }

    return patch;
}

template <size_t N, size_t D>
gsSubdivisionScheme::gsSubdivisionMeshValidity
gsFreeformSubdivision<N, D>::check_mesh()
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

template class gsFreeformSubdivision<5, 3>;
template class gsFreeformFaceData<5, 3>;
template class gsFreeformSubdivision<5, 4>;
template class gsFreeformFaceData<5, 4>;
template class gsFreeformSubdivision<6, 3>;
template class gsFreeformFaceData<6, 3>;
template class gsFreeformSubdivision<9, 3>;
template class gsFreeformFaceData<9, 3>;

} // namespace gismo
