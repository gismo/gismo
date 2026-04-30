/** @file gsFrogSplines.cpp

    @brief Classes for FROG Splines on a quadrangular mesh.

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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <gismo.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsMultiPatch.h>
#include <gsIO/gsFileData.h>
#include <gsIO/gsFileData.hpp>
#include <gsMesh2/gsFrogSplines.h>
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
void gsFrogSplines<N>::initialize_data(std::string filepath, size_t D)
{
    this->D = D;
    initialize_data(filepath);
}

template <size_t N>
void gsFrogSplines<N>::initialize_data(std::string filepath)
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
void gsFrogSplines<N>::initialize_data_xml(std::string filepath)
{
    auto& mesh = *m_mesh;

    gsMultiPatch<> patch;

    // Load the patches
    auto _file = gsReadFile<>(filepath, patch);

    // Convert each patch into a degree 1 (2x2) patch for mesh calculation
    gsMultiPatch<> corner_mp;
    for (size_t p = 0; p < patch.nPatches(); ++p)
    {
        auto& geo = patch.patch(p);
        auto& coefs = geo.coefs();

        // Extract 4 corners (coef index = i + j*n1)
        gsMatrix<real_t> cc(4, D);
        cc.row(0) = coefs.row(0);             // (0,    0   )
        cc.row(1) = coefs.row(N - 1);        // (N-1, 0   )
        cc.row(2) = coefs.row(N * (N - 1)); // (0,    N-1)
        cc.row(3) = coefs.row(N * N - 1);   // (N-1, N-1)

        // Clamped linear knot vector: {0,0,1,1}
        gsKnotVector<real_t> kv(0.0, 1.0, 0, 2);
        corner_mp.addPatch(gsTensorBSpline<2, real_t>(kv, kv, cc));
    }
    // Calculate the topology and the mesh
    corner_mp.computeTopology();
    // This mesh turn each patch into (n-1)^2 many patches, so by reducing n to 2 first we get one face per patch with correct connectivity.
    // This overrides the old mesh if present.
    mesh = corner_mp.toMesh();


    // Prepare the basis 
    gsKnotVector<real_t> kv_default(0.0, 1.0, 0, N);
    gsTensorBSplineBasis<2, real_t> basis_default(kv_default, kv_default);

    // Initialize the property on the mesh.
    auto face_data_vec = mesh.add_face_property(
        std::string("bezier_points"),
        gsPatch(basis_default, gsMatrix<>::Zero(N * N, D)));

    // Now associate each face with its patch.
    for (Face f : mesh.faces())
    {
        auto coefs = patch.patch(f.idx()).coefs();

        gsMatrix<> full_coefs = gsMatrix<>::Zero(N * N, D);
        const index_t n_cols = std::min((index_t)D, coefs.cols());
        full_coefs.leftCols(n_cols) = coefs.leftCols(n_cols);
        face_data_vec[f] = gsPatch(basis_default, full_coefs);
    }
}

template <size_t N>
void gsFrogSplines<N>::initialize_data_off(std::string filepath)
{
    auto& mesh = *m_mesh;
    // Clear the mesh
    mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(filepath, mesh);
    // Initialize the property.
    gsKnotVector<real_t> kv_default(0.0, 1.0, 0, N);
    gsTensorBSplineBasis<2, real_t> basis_default(kv_default, kv_default);
    mesh.add_face_property(std::string("bezier_points"),
                           gsPatch(basis_default, gsMatrix<>::Zero(N * N, D)));

    // Get the data.
    gsProperty<gsPatch> patch_data =
        mesh.get_face_property<gsPatch>("bezier_points");

    // Each patch is now initialized with a bilinear interpolation of its
    // corner vertices.
    for (auto f : mesh.faces())
    {
        // Collect the 4 corner positions in CCW order.
        std::vector<gsVector<real_t>> corners;
        for (Vertex v : mesh.vertices(f))
            corners.push_back(mesh.position(v));

        // Bilinear interpolation over an N×N grid.
        gsMatrix<real_t> coefs(N * N, D);
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                real_t u = real_t(j) / real_t(N - 1);
                real_t v = real_t(i) / real_t(N - 1);
                gsVector<real_t> pt = (1 - u) * (1 - v) * corners[0] +
                                      u * (1 - v) * corners[1] +
                                      (1 - u) * v * corners[3] +
                                      u * v * corners[2];
                coefs.row(i * N + j) = pt.transpose().leftCols(D);
            }
        }

        patch_data.vector()[f.idx()] = gsPatch(basis_default, coefs);
    }
}

template <size_t N> void gsFrogSplines<N>::scale(gsVector3d<> factors)
{
    auto& mesh = *m_mesh;
    gsProperty<gsPatch> face_data_vec(
        mesh.get_face_property<gsPatch>("bezier_points"));

    size_t n = std::min((size_t)factors.size(), (size_t)3);
    for (auto v : mesh.vertices())
    {
        mesh.position(v).head(n).array() *= factors.array().head(n);
    }

    for (auto& patch : face_data_vec.vector())
    {
        for (index_t c = 0;
             c < std::min<index_t>((index_t)factors.size(), patch.coefs().cols());
             ++c)
            patch.coefs().col(c) *= factors[c];
    }
}

template <size_t N>
gsMatrix<> gsFrogSplines<N>::rotate_coefs_cw(const gsMatrix<>& coefs,
                                                      size_t n)
{
    // CCW rotation: grid position (i,j) maps to (n-1-j, i).
    // Flat index i*n+j maps to (n-1-j)*n+i.
    gsMatrix<> result(coefs.rows(), coefs.cols());
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            result.row((n - 1 - j) * n + i) = coefs.row(i * n + j);
    return result;
}

template <size_t N>
gsMatrix<> gsFrogSplines<N>::rotate_coefs_ccw(const gsMatrix<>& coefs,
                                                       size_t n)
{
    // CW rotation: grid position (i,j) maps to (j, n-1-i).
    // Flat index i*n+j maps to j*n+(n-1-i).
    gsMatrix<> result(coefs.rows(), coefs.cols());
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            result.row(j * n + (n - 1 - i)) = coefs.row(i * n + j);
    return result;
}

template <size_t N>
gismo::gsTensorBSpline<2, real_t>
gsFrogSplines<N>::load_model_patch(int valence, std::string subtype)
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

template <size_t N> void gsFrogSplines<N>::orient_faces()
{
    // Get data
    auto& mesh = *m_mesh;
    gsProperty<gsPatch> face_data_vec(
        mesh.get_face_property<gsPatch>("bezier_points"));

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
                // also rotate the control points (CCW to match next_halfedge)
                gsPatch& p = face_data_vec.vector()[f.idx()];
                p.coefs() = rotate_coefs_cw(p.coefs(), N);
            }
        }
    }
}

template <size_t N> gsMultiPatch<> gsFrogSplines<N>::multipatch()
{
    auto& mesh = *m_mesh;
    gsMultiPatch<> patch;

    // Get the vector containing all the face data.
    gsProperty<gsPatch> face_data_vec =
        mesh.get_face_property<gsPatch>("bezier_points");

    // For each face, convert its control net to a patch and add it to the
    // multipatch. Order doesn't matter.
    for (auto face : mesh.faces())
    {
        patch.addPatch(face_data_vec.vector()[face.idx()]);
    }

    return patch;
}

template <size_t N>
std::array<gsSurfMesh::Face, 4>
gsFrogSplines<N>::order_faces(Vertex first_vertex,
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

template <size_t N>
gsSubdivisionScheme::gsSubdivisionMeshValidity
gsFrogSplines<N>::check_mesh()
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

template <size_t N> void gsFrogSplines<N>::subdivide()
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
    gsProperty<gsPatch> face_data_vec(
        mesh.get_face_property<gsPatch>("bezier_points"));

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
            // and remember the associated face data as a copy (not a
            // reference) to guard against aliasing when a child face
            // overwrites the parent's property slot.
            gsPatch coarse_patch = face_data_vec.vector()[initial_face];

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
                }

                // Fit a NxN Bezier patch to these samples
                gsKnotVector<> kv(0, 1, 0, N);
                gsTensorBSplineBasis<2> basis(kv, kv);
                gsFitting<> fitter(params, samples, basis);
                fitter.compute(0.0);
                gsGeometry<>* result = fitter.result();

                // The fine patches have their u direction pointing outwards
                // (first point is the center) while our system has the first
                // halfedge pointing outwards, so the first point is on the
                // edge. So we do a CCW rotation here.
                auto& data =
                    face_data_vec.vector()[children_ordered[f_idx].idx()];
                data.coefs() = rotate_coefs_cw(result->coefs(), N);
            }
        }
        else
        {
            // === ORDINARY VERTICES ===
            // Subdivide using gsPatch::splitAt to perform Bezier knot insertion

            // Copy the parent patch to avoid aliasing.
            const gsPatch src = face_data_vec.vector()[initial_face];

            // Split in direction 1 (i/row) then direction 0 (j/col) at 0.5
            // to obtain the four quadrant sub-patches.
            gsPatch top, bot, tl, tr, bl, br;
            src.splitAt(1, real_t(0.5), top, bot);
            top.splitAt(0, real_t(0.5), tl, tr);
            bot.splitAt(0, real_t(0.5), bl, br);

            // Collate the faces into a correctly ordered array.
            auto children_ordered =
                order_faces(first_vertices[initial_face], child_faces);

            // Write quadrant coefficients with the appropriate rotations.
            face_data_vec.vector()[children_ordered[0].idx()].coefs() =
                rotate_coefs_ccw(tl.coefs(), N);
            face_data_vec.vector()[children_ordered[1].idx()].coefs() =
                tr.coefs();
            face_data_vec.vector()[children_ordered[2].idx()].coefs() =
                rotate_coefs_cw(br.coefs(), N);
            face_data_vec.vector()[children_ordered[3].idx()].coefs() =
                rotate_coefs_ccw(rotate_coefs_ccw(bl.coefs(), N), N);
        }
    }
};

template <size_t N>
gsMatrix<real_t> gsFrogSplines<N>::smooth(size_t degree)
{
    GISMO_ASSERT(degree == 1, "Only C1 smoothing supported.");

    // ===============================
    // Phase 1: Project to smooth
    // ===============================

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

    // ===============================
    // Phase 2: Optimize Coefficients
    // ===============================

    auto& mesh = *m_mesh;
    const bool optimize_fit = m_options.getSwitch("optimize_fit");
    const std::string model_patch_path =
        m_options.getString("model_patch_path");

    std::vector<index_t> ev_coef_starts;
    std::vector<index_t> ev_coef_ends;
    ev_coef_starts.push_back(0);

    for (Vertex v : mesh.vertices())
    {
        if (is_ordinary(mesh, v))
            continue;

        ev_coef_ends.push_back(ev_coef_starts[ev_coef_starts.size() - 1] +
                               mesh.valence(v) * 15 + 1);
        ev_coef_starts.push_back(ev_coef_ends[ev_coef_ends.size() - 1] +
                                 20 * mesh.valence(v));
    }
    ev_coef_starts.pop_back();

    size_t ev_index = 0;
    for (Vertex v : mesh.vertices())
    {
        if (is_ordinary(mesh, v))
            continue;

        const size_t valence = mesh.valence(v);
        const index_t function_count = 15 * valence + 1;
        const index_t ev_coef_start = ev_coef_starts[ev_index];

        gsMatrix<real_t> kernel;
        auto _file = gsReadFile<>(model_patch_path + "Val" +
                                      std::to_string(valence) + "Kernel.xml",
                                  kernel);

        gsMatrix<real_t> functional;
        if (optimize_fit)
        {
            functional.resize(15 * valence, function_count);
            functional.setZero();
            for (size_t i = 0; i < function_count-1; ++i)
            {
                functional(i, 0) = real_t(1);
                functional(i, static_cast<index_t>(i) + 1) = real_t(-1);
            }
        }
        else
        {
            gsReadFile<>(model_patch_path + "Val" + std::to_string(valence) +
                             "Constraints.xml",
                         functional);
        }

        const gsMatrix<real_t> ev_coefficients =
            coefficients.block(ev_coef_start, 0, function_count, D);
        const gsMatrix<real_t> kernel_weights =
            (functional * kernel)
                .colPivHouseholderQr()
                .solve(-functional * ev_coefficients);

        coefficients.block(ev_coef_start, 0, function_count, D) =
            ev_coefficients + kernel * kernel_weights;

        ++ev_index;
    }

    gsMappedSpline<2, real_t> solSpline(mapped_basis, coefficients);
    gsMultiPatch<> solField = solSpline.exportToPatches();

    // Write the solution back.
    gsProperty<gsPatch> face_data_vec(
        m_mesh->get_face_property<gsPatch>("bezier_points"));
    for (size_t k = 0; k < face_data_vec.vector().size(); ++k)
        face_data_vec[k].coefs() = solField.patch(k).coefs();

    return coefficients;
}

template <size_t N>
void gsFrogSplines<N>::fit_function(gsFunctionExpr<real_t> function)
{
    gsMultiPatch<> multi_patch;
    gsMultiBasis<> multi_basis;
    gsMappedBasis<2> mapped_basis;
    this->c1_basis(multi_patch, multi_basis, mapped_basis);

    // Use the expression assembler to build a system
    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(multi_basis);

    auto G = A.getMap(multi_patch);
    auto u = A.getSpace(mapped_basis, 1);
    auto ff = A.getCoeff(function, G);

    u.setup();
    A.initSystem();
    // Equation: Int(v * u) = Int(v * f) e.g. u = f
    A.assemble(u * u.tr(), u * ff);

    // Solve this system
    gsSparseSolver<real_t>::LU solver;
    solver.compute(A.matrix());
    gsMatrix<> coefficients = solver.solve(A.rhs());
    auto solution = A.getSolution(u, coefficients);

    // Extract the solution into a gsMappedSpline, then export to patches.
    gsMappedSpline<2, real_t> solSpline;
    solution.extract(solSpline);
    gsMultiPatch<> solField = solSpline.exportToPatches();

    // Write the solution back into the last entry.
    gsProperty<gsPatch> face_data_vec(
        m_mesh->get_face_property<gsPatch>("bezier_points"));
    for (size_t k = 0; k < face_data_vec.vector().size(); ++k)
    {
        gsPatch& p = face_data_vec[k];
        p.coefs().conservativeResize(p.coefs().rows(), D + 1);
        p.coefs().col(D) = solField.patch(k).coefs().col(0);
    }

    D++;
}

template <size_t N>
void gsFrogSplines<N>::laplace_beltrami(gsFunctionExpr<real_t> rhs)
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
    // Evaluate the right hand side function $f$ at physical points on the
    // surface (mapped through G).  Using getCoeff(rhs) without G would
    // evaluate rhs at parametric (u,v) coordinates, which are identical for
    // every patch and would force a rotationally symmetric solution.
    auto ff = A.getCoeff(rhs, G);

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
    gsProperty<gsPatch> face_data_vec(
        m_mesh->get_face_property<gsPatch>("bezier_points"));
    for (size_t k = 0; k < face_data_vec.vector().size(); ++k)
    {
        gsPatch& p = face_data_vec[k];
        p.coefs().conservativeResize(p.coefs().rows(), D + 1);
        p.coefs().col(D) = solField.patch(k).coefs().col(0);
    }

    D++;
}

template <size_t N>
void gsFrogSplines<N>::c1_basis(gsMultiPatch<>& multi_patch,
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
        const size_t function_count = 15 * valence + 1;

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
                    // Inside EV support = any fitting functions non-zero here.
                    const bool in_support = std::any_of(
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
    // Corner rules:
    //
    // Non-boundary vertex:
    //   Set to the average of the adjacent (1,1) inner points.
    //
    // Boundary vertex with 1 adjacent patch:
    //   The corner remains a free DOF.
    //
    // Boundary vertex with 2 adjacent patches:
    //   Set to the average of the two neighboring control points on the
    //   boundary edges.
    //
    // Boundary vertex with more than 2 adjacent patches:
    //   Set to the average of one neighboring edge point per outgoing edge
    //   (thus adjacent_patch_count + 1 contributors).
    //
    // Points already handled by the EV support (Phase 1) are skipped.
    // ================================================================

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
                const size_t adjacent_patch_count = mesh.num_faces(v_corner);

                if (!mesh.is_boundary(v_corner))
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
                else if (adjacent_patch_count == 1)
                {
                    // This is a vertex on an outer corner
                    // It is thus completely free.
                    pre_mapper[ldof][global_dof_count] = real_t(1);
                    ++global_dof_count;
                }
                else if (adjacent_patch_count == 2)
                {
                    // This is a vertex on a corner of a patch that lies along
                    // the edge of the boundary.
                    // It is the average of the just-adjacent control points
                    // along the two edges that lie on the boundary.
                    for (Halfedge h_out : mesh.halfedges(v_corner))
                    {
                        auto h_in = mesh.opposite_halfedge(h_out);
                        if (!mesh.is_boundary(h_out) && mesh.is_boundary(h_in))
                        {
                            auto edge_ldof = local_dof_rotated(
                                mesh.face(h_out).idx(),
                                get_k(mesh.face(h_out), h_out), 0, 1);
                            for (auto& kv2 : pre_mapper[edge_ldof])
                                pre_mapper[ldof][kv2.first] +=
                                    real_t(0.5) * kv2.second;
                        }

                        if (mesh.is_boundary(h_out) && !mesh.is_boundary(h_in))
                        {
                            auto edge_ldof = local_dof_rotated(
                                mesh.face(h_in).idx(),
                                get_k(mesh.face(h_in), h_in), 0, N - 2);
                            for (auto& kv2 : pre_mapper[edge_ldof])
                                pre_mapper[ldof][kv2.first] +=
                                    real_t(0.5) * kv2.second;
                        }
                    }
                }
                else
                {
                    // All others are non-free vertices on patch corners that
                    // touch the boundary but are only interior corners of the
                    // boundary.

                    real_t factor = 1. / ((real_t)adjacent_patch_count + 1);

                    for (Halfedge h_out : mesh.halfedges(v_corner))
                    {
                        auto h_in = mesh.opposite_halfedge(h_out);
                        if (!mesh.is_boundary(h_out))
                        {
                            auto edge_ldof = local_dof_rotated(
                                mesh.face(h_out).idx(),
                                get_k(mesh.face(h_out), h_out), 0, 1);
                            for (auto& kv2 : pre_mapper[edge_ldof])
                                pre_mapper[ldof][kv2.first] +=
                                    factor * kv2.second;
                        }
                        else if (!mesh.is_boundary(h_in))
                        {
                            auto edge_ldof = local_dof_rotated(
                                mesh.face(h_in).idx(),
                                get_k(mesh.face(h_in), h_in), 0, N - 2);
                            for (auto& kv2 : pre_mapper[edge_ldof])
                                pre_mapper[ldof][kv2.first] +=
                                    factor * kv2.second;
                        }
                    }
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
gsVector<real_t, 2>
gsFrogSplines<N>::error(gsFunctionExpr<real_t> function,
                                size_t samples_per_face)
{

    auto& mesh = *m_mesh;
    size_t spf = samples_per_face;
    gsProperty<gsPatch> face_data_vec(
        mesh.get_face_property<gsPatch>("bezier_points"));

    real_t error_linf(0.);
    real_t error_l2(0.);
    real_t error_count(0);

    for (auto f : mesh.faces())
    {
        // Look at the patch for this face.
        const gsPatch& patch = face_data_vec[f.idx()];

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

    gsInfo << "Error LI: " << error(0) << ".\n";
    gsInfo << "Error L2: " << error(1) << ".\n";

    return error;
}

template <size_t N>
void gsFrogSplines<N>::write_paraview_error(
    gsFunctionExpr<real_t> function, real_t max_error, std::string name,
    gsParaviewCollection* collection, size_t timestep)
{
    // Prepare data
    auto& mesh = *m_mesh;
    gsProperty<gsPatch> face_data_vec(
        mesh.get_face_property<gsPatch>("bezier_points"));

    // Number of faces needs to be counted for registration in the collection.
    size_t face_counter(0);

    for (auto f : mesh.faces())
    {
        // read out the control points and create a patch
        const gsPatch& patch = face_data_vec[f.idx()];

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
void gsFrogSplines<N>::write_paraview(
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



template class gsFrogSplines<5>;
template class gsFrogSplines<6>;
template class gsFrogSplines<9>;

} // namespace gismo
