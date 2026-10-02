/** @file poisson2_nitsche_immersed_geometry_example.cpp

    @brief Poisson equation on an immersed 3D domain (a triangle mesh loaded
           from a Wavefront .obj file) inside [0,1]^3, solved with the Finite
           Cell Method.

    Pipeline (Poisson analogue of calculate_volume_immersed_example.cpp):

      1.-4. Identical setup to calculate_volume_immersed_example.cpp: load and
            rescale the .obj mesh into [0,1]^3, build the mesh's
            signed-distance level set, classify the background grid with
            gsImplicitTrimmedDomain<3>.
      5. Assemble the Poisson system explicitly (no gsExprAssembler), reusing
         the same cell loop / Algoim quadrature (+ uniform-grid fallback on
         re-entrant cut cells where Algoim returns zero nodes) as the volume
         example:
           stiffness_ij = int_Omega grad(phi_i) . grad(phi_j)
           rhs_i        = int_Omega f * phi_i
         plus a tiny full-background regularization term so DOFs whose
         support never touches Omega do not leave the matrix singular.
      6. Weak Dirichlet / Nitsche terms on the immersed boundary {phi==0}
         (on by default, --noNitsche to disable), using the Algoim surface
         quadrature rule and the level-set-gradient normal
         n = grad(phi)/|grad(phi)| (finite differences via
         gsMeshSignedDist::deriv_into).
      7. Per refinement level: solve, compute L2/H1 errors against the
         manufactured solution over {phi<0}, and (--plot) export the
         quadrature points actually used (same 4 x N [xyz; weight] layout as
         the volume example), tagged with the refinement level.

    Usage:
    ./bin/poisson2_nitsche_immersed_geometry_example -f obj/spot.obj -r 2 -g 1e3 --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsDomain/gsMeshLevelSet.h>

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
//  main
// =============================================================================
int main(int argc, char * argv[])
{
    std::string filename =
        "obj/spot.obj";
    std::string out = gsFileManager::getCanonicRepresentation(
        gsFileManager::getPath(std::string(__FILE__), true) + "../output_poisson_immersed3d");
    index_t numRefine  = 3;
    index_t numElevate = 0;
    real_t  fill       = 0.9;  // fraction of [0,1]^3 spanned by the largest extent
    real_t  gamma      = 1e3;  // Nitsche penalty parameter (scaled by 1/h)
    bool    plot       = false;
    bool    noNitsche  = false;

    gsCmdLine cmd("Poisson equation on an immersed 3D .obj domain via the Finite Cell Method.");
    cmd.addString("f", "file",            "Wavefront .obj mesh file",            filename);
    cmd.addString("o", "output",          "Output folder for ParaView files",    out);
    cmd.addInt   ("r", "uniformRefine",   "Number of uniform refinement steps",  numRefine);
    cmd.addInt   ("e", "degreeElevation", "Number of degree elevation steps",    numElevate);
    cmd.addReal  ("",  "fill",            "Fill fraction of [0,1]^3 (0..1)",      fill);
    cmd.addReal  ("g", "gamma",           "Nitsche penalty parameter",           gamma);
    cmd.addSwitch("plot",      "Create ParaView output",                        plot);
    cmd.addSwitch("noNitsche", "Disable the weak/Nitsche BC (interior forcing only)", noNitsche);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    out = gsFileManager::getCanonicRepresentation(out);

    // -------------------------------------------------------------------------
    // 1. Load the triangle mesh
    // -------------------------------------------------------------------------
    gsSurfMesh mesh;
    if (!gsReadSurfMesh(filename, mesh))
    {
        gsWarn << "Failed to read a triangle mesh from: " << filename << "\n";
        return EXIT_FAILURE;
    }
    const std::string outputStem = gsFileManager::getBasename(filename);
    const std::size_t nVert = mesh.n_vertices();
    const std::size_t nTri  = mesh.n_faces();
    gsInfo << "Loaded mesh '" << filename << "': "
           << nVert << " vertices, " << nTri << " triangles.\n";

    // -------------------------------------------------------------------------
    // 2. Rescale / center the mesh into the unit parametric box [0,1]^3
    //    mapped = (p - center) * scale + 0.5,  scale = fill / maxExtent
    // -------------------------------------------------------------------------
    const real_t scale = gsNormalizeToUnitBox(mesh, fill);
    GISMO_UNUSED(scale);

    // -------------------------------------------------------------------------
    // 3. Analytic level set on the rescaled mesh
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(3, 2);
    bbox.col(0).setZero();
    bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(mesh, bbox);

    // -------------------------------------------------------------------------
    // 4. Background box [0,1]^3 (identity map: parameter space == physical space)
    // -------------------------------------------------------------------------
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    const index_t deg = dbasis.maxCwiseDegree();

    // Manufactured solution: f = -Delta(u_exact), imposed weakly on {phi==0}.
    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);
    gsFunctionExpr<> f_rhs("3*pi^2*sin(pi*x)*sin(pi*y)*sin(pi*z)", 3);

    gsFileManager::mkdir(out);

    gsInfo << "\nDegree of basis: " << deg << "\n";
    gsInfo << "Solving the Poisson equation on the immersed domain in [0,1]^3.\n";
    if (noNitsche)
        gsInfo << "Weak/Nitsche BC disabled: only the interior forcing is imposed.\n";
    gsInfo << "\nDoFs (dot1=assembled, dot2=solved, dot3=error): ";

    // Number of uniform sample points per direction used as fallback when
    // Algoim returns 0 nodes (re-entrant / multi-crossing cut cells).
    const index_t nFallback = 4;
    const int D = 3;

    gsVector<real_t> l2err(numRefine + 1), h1err(numRefine + 1);

    // -------------------------------------------------------------------------
    // 5. Refinement loop: classify cells, assemble, solve, compute errors
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        gsTensorBSplineBasis<3, real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<3, real_t> *>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        real_t hmax = 0;
        for (std::size_t p = 0; p != dbasis.nBases(); ++p)
            hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());

        // Cut-cell classification of the background grid.
        gsImplicitTrimmedDomain<3, real_t> tr_domain(impl_fun, *tbsPtr);

        // Quadrature rules -- exactly the same objects as the volume example.
        gsGaussRule<real_t>         gauss(gsVector<index_t, 3>::Constant(deg + 1));
        gsAlgoimGenericRule<real_t> algoimRule(impl_fun, *tbsPtr);            // {phi<0} volume

        gsOptionList surfOpts = gsAlgoimGenericRule<real_t>::defaultOptions();
        surfOpts.setInt("dim", D);
        gsAlgoimGenericRule<real_t> algoimSurf(impl_fun, *tbsPtr, surfOpts);  // {phi==0} surface

        // Cut-cell volume quadrature = Algoim rule, with the same
        // uniform-grid fallback as the volume example when Algoim returns
        // zero nodes (re-entrant / multi-crossing cells).
        auto cutCellPoints = [&](const gsVector<real_t> & lo_c, const gsVector<real_t> & hi_c,
                                  gsMatrix<real_t> & pts, gsVector<real_t> & wts)
        {
            algoimRule.mapTo(lo_c, hi_c, pts, wts);
            if (wts.size() > 0) return;

            const index_t nTotal = nFallback * nFallback * nFallback;
            gsMatrix<real_t> spts(3, nTotal);
            index_t col = 0;
            for (index_t ii = 0; ii < nFallback; ++ii)
            for (index_t jj = 0; jj < nFallback; ++jj)
            for (index_t kk = 0; kk < nFallback; ++kk, ++col)
            {
                spts(0,col) = lo_c[0] + (ii+0.5)/nFallback*(hi_c[0]-lo_c[0]);
                spts(1,col) = lo_c[1] + (jj+0.5)/nFallback*(hi_c[1]-lo_c[1]);
                spts(2,col) = lo_c[2] + (kk+0.5)/nFallback*(hi_c[2]-lo_c[2]);
            }
            gsMatrix<real_t> fbVals;
            impl_fun.eval_into(spts, fbVals);
            std::vector<index_t> insideCols;
            for (index_t ci = 0; ci < nTotal; ++ci)
                if (fbVals(0,ci) < 0) insideCols.push_back(ci);
            if (insideCols.empty()) { pts.resize(3,0); wts.resize(0); return; }

            const real_t cellVol = (hi_c - lo_c).prod();
            const real_t w = cellVol / static_cast<real_t>(nTotal);
            pts.resize(3, static_cast<index_t>(insideCols.size()));
            wts.resize(static_cast<index_t>(insideCols.size()));
            for (std::size_t ci = 0; ci < insideCols.size(); ++ci)
            {
                pts.col(ci) = spts.col(insideCols[ci]);
                wts[ci] = w;
            }
        };

        // ---------------------------------------------------------------
        // Explicit assembly. The geometry map is the identity on [0,1]^3,
        // so parametric gradients equal physical gradients (Jacobian = I).
        // ---------------------------------------------------------------
        const index_t nDofs = tbsPtr->size();
        gsVector<real_t> F_vec(nDofs);
        F_vec.setZero();
        gsSparseEntries<real_t> triplets;

        // Accumulate stiffness (+ load) from one batch of volume points.
        auto addVolume = [&](const gsMatrix<real_t> & pts, const gsVector<real_t> & wts,
                              real_t coeff = (real_t)1, bool withRhs = true)
        {
            if (pts.cols() == 0) return;
            gsMatrix<index_t> act;
            gsMatrix<real_t>  bv, bd, fv;
            tbsPtr->active_into(pts, act);
            tbsPtr->eval_into  (pts, bv);
            tbsPtr->deriv_into (pts, bd);
            if (withRhs) f_rhs.eval_into(pts, fv);
            const index_t na = act.rows();
            for (index_t q = 0; q < pts.cols(); ++q)
            {
                const real_t w = wts(q) * coeff;
                for (index_t i = 0; i < na; ++i)
                {
                    const index_t di = act(i,q);
                    if (withRhs)
                        F_vec(di) += w * bv(i,q) * fv(0,q);
                    for (index_t j = 0; j < na; ++j)
                    {
                        real_t kij = 0;
                        for (int d = 0; d < D; ++d)
                            kij += bd(i*D+d,q) * bd(j*D+d,q);
                        triplets.add(di, act(j,q), w * kij);
                    }
                }
            }
        };

        gsMatrix<real_t> pts_tmp; gsVector<real_t> wts_tmp;

        // Tiny full-background regularization: keeps K nonsingular for DOFs
        // whose support never touches {phi<0}.
        {
            const real_t alpha = 1e-10;
            for (auto it = tbsPtr->domain()->beginAll(); it != tbsPtr->domain()->endAll(); ++it)
            {
                gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
                addVolume(pts_tmp, wts_tmp, alpha, /*withRhs=*/false);
            }
        }

        // int_Omega grad(u).grad(v)  and  int_Omega f.v,  Omega = {phi<0}.
        for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addVolume(pts_tmp, wts_tmp);
        }
        for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
        {
            cutCellPoints(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addVolume(pts_tmp, wts_tmp);
        }

        // ---------------------------------------------------------------
        // 6. Weak (Nitsche) immersed boundary condition on {phi==0}.
        //    Normal: n = grad(phi)/|grad(phi)| (finite differences via
        //    gsMeshSignedDist). Disable with --noNitsche if this ever proves
        //    fragile -- without it, the interior forcing is still solved
        //    correctly, but u=0 is not imposed on the immersed boundary, so
        //    the solution will not converge to u_exact.
        // ---------------------------------------------------------------
        if (!noNitsche)
        {
            for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
            {
                gsMatrix<real_t> spts; gsVector<real_t> swts;
                algoimSurf.mapTo(it.lowerCorner(), it.upperCorner(), spts, swts);
                if (spts.cols() == 0) continue;

                gsMatrix<index_t> act;
                gsMatrix<real_t>  bv, bd, phi_d, gDv;
                tbsPtr->active_into(spts, act);
                tbsPtr->eval_into  (spts, bv);
                tbsPtr->deriv_into (spts, bd);
                impl_fun.deriv_into(spts, phi_d);   // 3 x numSurfPts (finite differences)
                u_exact.eval_into  (spts, gDv);      // 1 x numSurfPts
                const index_t na = act.rows();

                for (index_t q = 0; q < spts.cols(); ++q)
                {
                    const real_t w = swts(q);
                    real_t nNorm2 = 0;
                    for (int d = 0; d < D; ++d)
                        nNorm2 += phi_d(d,q) * phi_d(d,q);
                    const real_t nNorm = std::sqrt(nNorm2);
                    if (nNorm < 1e-14) continue;

                    for (index_t i = 0; i < na; ++i)
                    {
                        const index_t di = act(i,q);
                        const real_t  vi = bv(i,q);
                        real_t gi_n = 0;
                        for (int d = 0; d < D; ++d)
                            gi_n += bd(i*D+d,q) * phi_d(d,q) / nNorm;

                        // RHS: -(grad(phi_i).n)g_D + (gamma/h)phi_i g_D
                        F_vec(di) += w * (-gi_n * gDv(0,q) + (gamma/hmax) * vi * gDv(0,q));

                        for (index_t j = 0; j < na; ++j)
                        {
                            const index_t dj = act(j,q);
                            const real_t  vj = bv(j,q);
                            real_t gj_n = 0;
                            for (int d = 0; d < D; ++d)
                                gj_n += bd(j*D+d,q) * phi_d(d,q) / nNorm;

                            // -(gi.n)vj - vi(gj.n) + (gamma/h)vi.vj
                            const real_t kij = -gi_n*vj - vi*gj_n + (gamma/hmax)*vi*vj;
                            triplets.add(di, dj, w * kij);
                        }
                    }
                }
            }
        }

        // Assemble sparse matrix and solve.
        gsSparseMatrix<real_t> K(nDofs, nDofs);
        K.setFrom(triplets);
        triplets.clear();

        gsInfo << nDofs << "." << std::flush;

        gsVector<real_t> solVector(nDofs);
        {
#           ifdef GISMO_WITH_PARDISO
                gsSparseSolver<real_t>::PardisoLDLT slvr;
#           else
                gsSparseSolver<real_t>::CGDiagonal slvr;
#           endif
            slvr.compute(K);
            solVector = slvr.solve(F_vec);
        }
        gsInfo << "." << std::flush;

        // Error norms, integrated over {phi<0} with the same quadrature.
        real_t l2sq = 0, h1sq = 0;
        auto addError = [&](const gsMatrix<real_t> & pts, const gsVector<real_t> & wts)
        {
            if (pts.cols() == 0) return;
            gsMatrix<index_t> act;
            gsMatrix<real_t>  bv, bd, uExV, uExD;
            tbsPtr->active_into(pts, act);
            tbsPtr->eval_into  (pts, bv);
            tbsPtr->deriv_into (pts, bd);
            u_exact.eval_into  (pts, uExV);        // 1 x numPts
            u_exact.deriv_into (pts, uExD);        // 3 x numPts
            const index_t na = act.rows();
            for (index_t q = 0; q < pts.cols(); ++q)
            {
                const real_t w = wts(q);
                real_t uh = 0;
                for (index_t i = 0; i < na; ++i)
                    uh += bv(i,q) * solVector(act(i,q));
                const real_t dl2 = uh - uExV(0,q);
                l2sq += w * dl2 * dl2;
                for (int d = 0; d < D; ++d)
                {
                    real_t guh_d = 0;
                    for (index_t i = 0; i < na; ++i)
                        guh_d += bd(i*D+d,q) * solVector(act(i,q));
                    const real_t dh1 = guh_d - uExD(d,q);
                    h1sq += w * dh1 * dh1;
                }
            }
        };

        for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addError(pts_tmp, wts_tmp);
        }
        for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
        {
            cutCellPoints(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addError(pts_tmp, wts_tmp);
        }
        l2err[r] = math::sqrt(l2sq);
        h1err[r] = l2err[r] + math::sqrt(h1sq);
        gsInfo << ". " << std::flush;

        // ---------------------------------------------------------------
        // 7. ParaView output: the quadrature points actually used for the
        //    volume integrals (same 4 x N [xyz; weight] layout as the
        //    volume example), tagged with the refinement level.
        // ---------------------------------------------------------------
        if (plot)
        {
            gsMatrix<real_t> pts, phys, quad(4, 0);
            gsVector<real_t> wts;

            for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
            {
                cutCellPoints(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0) continue;
                mp.patch(0).eval_into(pts, phys);
                const index_t c = quad.cols();
                quad.conservativeResize(4, c + pts.cols());
                quad.block(0, c, 3, pts.cols()) = phys;
                quad.row(3).segment(c, pts.cols()) = wts.transpose();
            }
            for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
            {
                gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0) continue;
                mp.patch(0).eval_into(pts, phys);
                const index_t c = quad.cols();
                quad.conservativeResize(4, c + pts.cols());
                quad.block(0, c, 3, pts.cols()) = phys;
                quad.row(3).segment(c, pts.cols()) = wts.transpose();
            }
            if (quad.cols() > 0)
                gsWriteParaviewPoints(quad, out + "/quadrature_points_r" + std::to_string(r));

            if (r == numRefine)
            {
            gsWriteParaview(mesh, out + "/" + outputStem + "_geometry");

                gsMesh<real_t> bgMesh(dbasis.basis(0));
                gsWriteParaview(bgMesh, out + "/" + outputStem + "_background_mesh");
            }
        }

        gsInfo << " (h=" << hmax << ", gamma/h=" << gamma/hmax << ")\n";
    }

    gsInfo << "\nL2 error: " << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
    gsInfo << "H1 error: " << std::scientific << h1err.transpose() << "\n";

    if (numRefine > 0)
    {
        gsInfo << "\nEoC (L2): " << std::fixed << std::setprecision(2)
               << (l2err.head(numRefine).array() / l2err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
        gsInfo << "EoC (H1): " << std::fixed << std::setprecision(2)
               << (h1err.head(numRefine).array() / h1err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
    }

    if (plot)
        gsInfo << "\nParaView files written to: " << out << "\n";

    return EXIT_SUCCESS;
}
