/** @file gsFreeformSubdivision.cpp

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsDebug.h"
#include <gsCore/gsMultiPatch.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSubdivisionScheme.h>
#include <gsNurbs/gsTensorBSpline.h>

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
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            real_t denom = real_t((N - 1) * (N - 1));
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
    result = rotate_r(result);
    // find the edge on the face
    for (auto const& he : mesh.halfedges(face))
    {
        if (he == hedge)
            break;

        result = rotate_l(result);
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
    gsMatrix<> coeffs(N * N, 3);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            int total_index = i * N + j;
            coeffs(total_index, 0) = control_points(i, j).x();
            coeffs(total_index, 1) = D >= 2 ? control_points(i, j).y() : 0.0;
            coeffs(total_index, 2) = D >= 3 ? control_points(i, j).z() : 0.0;
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
void gsFreeformSubdivision<N, D>::subdivide(gsSurfMesh& mesh)
{
    // As a pre-step, we make sure that for each face, if it has an
    // extraordinary vertex, that vertex is the top left one
    for (Face f : mesh.faces())
    {
        // first, check if this face has an EV at all
        bool has_ev(false);
        for (Vertex v : mesh.vertices(f))
        {
            has_ev = has_ev || !is_ordinary(mesh, v);
        }

        // if it does, rotate its first halfedge around until it points to the
        // EV
        if (has_ev)
        {
            while (is_ordinary(mesh, mesh.to_vertex(mesh.halfedge(f))))
            {
                mesh.set_halfedge(f, mesh.next_halfedge(mesh.halfedge(f)));
            }
        }
    }

    // Remember the first vertex of each face (this is where the control nets of
    // each face data are oriented on).
    std::vector<Vertex> first_vertices;
    // List to keep track of extraordinary vertices on each face.
    // An invalid vertex with idx -1 is used to signify no EV on a face.
    // Multiple EVs on a face should not happen, use a better mesh.
    std::vector<Vertex> extraordinary_vertices;
    for (Face f : mesh.faces())
    {
        // remember the first vertex
        first_vertices.emplace_back(mesh.to_vertex(mesh.halfedge(f)));

        // reserve a spot for possible EVs by placing an invalid vertex that
        // will also serve as a marker for faces without EVs.
        extraordinary_vertices.emplace_back(-1);

        // check if there is an EV
        for (const Vertex& v : mesh.vertices(f))
        {
            // if an EV is found, replace the invalid vertex with it
            if (!is_ordinary(mesh, v))
            {
                extraordinary_vertices.back() = v;
            }
        }
    }

    // Split each face into 4 and get info about the way they were split.
    std::map<gsSurfMesh::Face, std::vector<gsSurfMesh::Face>> face_map =
        mesh.quad_split();

    // Get face data
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points"));

    // Now fix the data on each face.
    for (auto const& parent_to_children_faces : face_map)
    {
        if (extraordinary_vertices[parent_to_children_faces.first.idx()]
                .is_valid())
        {
            Vertex ev =
                extraordinary_vertices[parent_to_children_faces.first.idx()];
            Vertex fv = first_vertices[parent_to_children_faces.first.idx()];
            gsInfo << parent_to_children_faces.first << " has an EV: " << ev
                   << ".\n";
            gsInfo << parent_to_children_faces.first
                   << " has first vertex: " << fv << ".\n";
        }
        // else //un-comment this once ev subdivision is done
        {
            // === ORDINARY VERTICES ===
            // via deCasteljau

            // Get the face data and store it in a temporary dynamic 2d array.
            gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic> control_net(
                face_data_vec.vector()[parent_to_children_faces.first.idx()]
                    .control_points);

            // now control_net is a (n+1)*(n+1) matrix of control points (degree
            // n) Perform deCasteljau once to divide into two (n+1)*(n+1)
            // matrices of control points.
            auto const first_split = this->deCasteljau(control_net);

            // Perform deCasteljau again on both of them, to get 4 (n+1)*(n+1)
            // matrices of control points In between, we need to transpose so we
            // now divide in the other direction.
            auto top_split = this->deCasteljau(first_split[0].transpose());
            auto bot_split = this->deCasteljau(first_split[1].transpose());

            // re-transpose
            top_split[0] = top_split[0].transpose().eval();
            top_split[1] = top_split[1].transpose().eval();
            bot_split[0] = bot_split[0].transpose().eval();
            bot_split[1] = bot_split[1].transpose().eval();

            // rotate
            bot_split[1] = rotate_l(bot_split[1]);
            bot_split[0] = rotate_l(rotate_l(bot_split[0]));
            top_split[0] = rotate_l(rotate_l(rotate_l(top_split[0])));

            // Collate all these matrices in the correct order into an array.
            std::array<gsMatrix<gsVector<real_t, D>, Dynamic, Dynamic>, 4> arr =
                {top_split[0], top_split[1], bot_split[1], bot_split[0]};

            // find the new face that contains the top_left vertex
            Vertex first_vertex =
                first_vertices[parent_to_children_faces.first.idx()];
            size_t first_face(0);
            for (size_t i = 0; i < 4; ++i)
            {
                // get the face
                Face f(parent_to_children_faces.second[i]);
                // go through the vertices of this face and check if it has the
                // searched-for vertex
                for (auto const& v : mesh.vertices(f))
                {
                    if (v == first_vertex)
                        first_face = i;
                }
            }

            // Collate the faces into a correctly ordered array as well.
            std::array<Face, 4> children_faces_ordered;
            for (int i = 0; i < 4; ++i)
            {
                children_faces_ordered[i] =
                    parent_to_children_faces.second[(i + first_face) % 4];
            }

            // Correct back references of face data and give them the correct
            // control points.
            for (size_t f = 0; f < 4; ++f)
            {
                auto data =
                    &face_data_vec.vector()[children_faces_ordered[f].idx()];
                data->face = children_faces_ordered[f];
                data->control_points = arr[f];
            }
        }
    }
};

template <size_t N, size_t D>
bool gsFreeformSubdivision<N, D>::is_ordinary(const gsSurfMesh& mesh,
                                              const Vertex& v)
{
    auto count(0);
    for ([[maybe_unused]] const auto he : mesh.halfedges(v))
    {
        ++count;
    }
    return count == 4 || mesh.is_boundary(v);
}

template <size_t N, size_t D>
void gsFreeformSubdivision<N, D>::smooth(gsSurfMesh& mesh, size_t degree)
{
    // Ensure we have enough degree
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

    // Get face data
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec(
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points"));

    // First, correct each vertex
    for (const Vertex& v : mesh.vertices())
    {
        // ignore EVs
        if (!is_ordinary(mesh, v))
            continue;

        // first, collect all the control points in a square of side length 1
        // around this vertex. Note that each control point on a boundary is
        // represented up to twice, the one in the center even up to four times.
        // In the end we get up to 9 distinct points
        // Points are arranged like this:
        // 0 1 4
        // 2 3 5
        // 8 7 6
        std::vector<gsMatrix<gsVector<real_t, D>*>> control_points_faces;
        for (Halfedge h : mesh.halfedges(v))
        {
            // we do skip non-existent faces
            if (mesh.is_boundary(h))
                continue;
            control_points_faces.emplace_back(
                face_data_vec[mesh.face(h).idx()].control_points_oriented(mesh,
                                                                          h));
        }

        // we have 4 base equations at 1 face, then two more per face, except
        // only 1 more for the last face.
        size_t eqs(2 + 2 * control_points_faces.size() +
                   (control_points_faces.size() == 4 ? -1 : 0));

        // create a matrix that represents the C1 equations, i.e. each point on
        // a boundary is colinear with the ones on either side. This results in
        // multiple equations for the center point.
        // In the end, we will have four free points and all others should be
        // dependent.
        auto matrix = gsMatrix<real_t>(eqs, 4);
        // The first four equations just say that the first four points are
        // equal to themselves.
        matrix.row(0) << 1., 0., 0., 0.;
        matrix.row(1) << 0., 1., 0., 0.;
        matrix.row(2) << 0., 0., 1., 0.;
        matrix.row(3) << 0., 0., 0., 1.;

        // if a second face is present, its points must be colinear with the
        // first face
        if (control_points_faces.size() >= 2)
        {
            matrix.row(4) << -1., 2., 0., 0.;
            matrix.row(5) << 0., 0., -1., 2.;
        }

        // if a third face is present, its points must be colinear with the
        // second face, which are in turn expressed via the first face.
        if (control_points_faces.size() >= 3)
        {
            matrix.row(6) << 1., -2., -2., 4.;
            matrix.row(7) << 0., -1., 0., 2.;
        }

        // if a fourth face is present, its points must be colinear with the
        // first and third face.
        if (control_points_faces.size() >= 4)
        {
            matrix.row(8) << -1., 0., 2., 0.;
        }

        // Now do a least squares fit.
        // I.e. we are searching for the set of 9 points (or rather 4 points
        // determining 9) that have the minimum distance to the original 9
        // points while fulfilling all equations (this latter part is achieved
        // by only searching for 4 points).

        // Set the target points as row vectors of a matrix.
        // The number of rows depends on the faces present.
        gsMatrix<real_t> target_matrix(eqs, D);
        target_matrix.setZero();
        target_matrix.row(0) = control_points_faces[0](1, 1)->transpose();
        target_matrix.row(1) = control_points_faces[0](1, 0)->transpose();
        target_matrix.row(2) = control_points_faces[0](0, 1)->transpose();
        target_matrix.row(3) = control_points_faces[0](0, 0)->transpose();

        if (control_points_faces.size() >= 2)
        {
            target_matrix.row(4) = control_points_faces[1](1, 1)->transpose();
            target_matrix.row(5) = control_points_faces[1](1, 0)->transpose();
        }

        if (control_points_faces.size() >= 3)
        {
            target_matrix.row(6) = control_points_faces[2](1, 1)->transpose();
            target_matrix.row(7) = control_points_faces[2](1, 0)->transpose();
        }

        if (control_points_faces.size() >= 4)
        {
            target_matrix.row(8) = control_points_faces[3](1, 1)->transpose();
        }

        // actually perform the least squares
        gsMatrix<real_t> solution =
            matrix.colPivHouseholderQr().solve(target_matrix);

        // assign the solutions, again depending on faces present
        *(control_points_faces[0](1, 1)) = solution.row(0);
        *(control_points_faces[0](1, 0)) = solution.row(1);
        *(control_points_faces[0](0, 1)) = solution.row(2);
        *(control_points_faces[0](0, 0)) = solution.row(3);

        if (control_points_faces.size() >= 2)
        {
            *(control_points_faces[1](0, 0)) = solution.row(3);
            *(control_points_faces[1](0, 1)) = solution.row(1);
            *(control_points_faces[1](1, 0)) =
                2. * solution.row(3) - solution.row(2);
            *(control_points_faces[1](1, 1)) =
                2. * solution.row(1) - solution.row(0);
        }

        if (control_points_faces.size() >= 3)
        {
            *(control_points_faces[2](0, 0)) = solution.row(3);
            *(control_points_faces[2](0, 1)) =
                2. * solution.row(3) - solution.row(2);
            *(control_points_faces[2](1, 0)) =
                2. * solution.row(3) - solution.row(1);
            *(control_points_faces[2](1, 1)) =
                4. * solution.row(3) - 2. * solution.row(2) -
                2. * solution.row(1) + solution.row(0);
        }

        if (control_points_faces.size() >= 4)
        {
            *(control_points_faces[3](0, 0)) = solution.row(3);
            *(control_points_faces[3](0, 1)) =
                2. * solution.row(3) - solution.row(1);
            *(control_points_faces[3](1, 0)) = solution.row(2);
            *(control_points_faces[3](1, 1)) =
                2. * solution.row(2) - solution.row(0);
        }
    }
    // End of edge correction

    // Now, correct the remaining part of each edge that isn't surrounding a
    // vertex.
    for (const Edge& e : mesh.edges())
    {
        auto halfedge0 = mesh.halfedge(e, 0);
        auto halfedge1 = mesh.halfedge(e, 1);

        // if we are on the boundary of the mesh, nothing needs to be done
        if (mesh.is_boundary(halfedge0) || mesh.is_boundary(halfedge1))
            continue;

        // Get the faces. These must be valid or else we would have continued
        // above.
        auto face0 = mesh.face(halfedge0);
        auto face1 = mesh.face(halfedge1);

        // Get the control points correctly oriented with respect to these
        // halfedges.
        auto cp0 =
            face_data_vec[face0.idx()].control_points_oriented(mesh, halfedge0);
        auto cp1 =
            face_data_vec[face1.idx()].control_points_oriented(mesh, halfedge1);

        // correct all points but the first and last two
        for (size_t i = 2; i < N - 2; ++i)
        {
            // we get a triple of points across boundary. The cetner point is of
            // course represented in both meshes.
            gsVector<real_t, D>* i0 = cp0(1, i);
            gsVector<real_t, D>* m0 = cp0(0, i);
            gsVector<real_t, D>* m1 = cp1(0, N - 1 - i);
            gsVector<real_t, D>* i1 = cp1(1, N - 1 - i);

            // make sure the user already put in a C0 mesh
            if (*m0 != *m1)
                gsWarn << "Points along an edge are not equal, mesh might be "
                          "corrupted.\n";

            // generate the matrix that describes the C1 equation for these
            // points, i.e. the center point is colinear with the outer two
            auto matrix = gsMatrix<real_t, 3, 2>({1., 0.5, 0., 0., 0.5, 1.});

            // We now do a least squares fit, i.e. search for three points as
            // close as possible to the original three that are colinear. We do
            // this by putting the target vectors as row of a matrix, getting
            // our solutions back as rows as well.

            gsMatrix<real_t> target_matrix(3, D);
            target_matrix.setZero();
            target_matrix.row(0) = i0->transpose();
            target_matrix.row(1) = m0->transpose();
            target_matrix.row(2) = i1->transpose();

            gsMatrix<real_t> solution =
                matrix.colPivHouseholderQr().solve(target_matrix);

            *i0 = solution.row(0);
            *m0 = (solution.row(0) + solution.row(1)) * 0.5;
            *m1 = (solution.row(0) + solution.row(1)) * 0.5;
            *i1 = solution.row(1);
        }
    }
    // End of vertex correction
}

template <size_t N, size_t D>
void gsFreeformSubdivision<N, D>::initialize_data(gsSurfMesh& mesh)
{
    mesh.add_face_property(std::string("bezier_points"),
                           gsFreeformFaceData<N, D>());
    gsProperty<gsFreeformFaceData<N, D>> patch_data =
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points");
    for (auto f : mesh.faces())
    {
        patch_data.vector()[f.idx()] = gsFreeformFaceData<N, D>(mesh, f);
    }
}

template <size_t N, size_t D>
gsMultiPatch<> gsFreeformSubdivision<N, D>::multipatch(const gsSurfMesh& mesh)
{
    gsMultiPatch<> patch;

    // get the vector containing all the face data
    gsProperty<gsFreeformFaceData<N, D>> face_data_vec =
        mesh.get_face_property<gsFreeformFaceData<N, D>>("bezier_points");

    // for each face, convert its control net to a patch and add it to the
    // multipatch. Order doesn't matter.
    for (auto face : mesh.faces())
    {
        patch.addPatch(face_data_vec.vector()[face.idx()].patch());
    }

    return patch;
}

template <size_t N, size_t D>
gsSubdivisionScheme::gsSubdivisionMeshValidity
gsFreeformSubdivision<N, D>::valid_mesh(const gsSurfMesh& mesh)
{
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
