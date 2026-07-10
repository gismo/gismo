/** @file poisson2_nitsche_immersed_geometry_octree_example.cpp

    @brief Poisson equation on an immersed 3D domain (a triangle mesh loaded
           from a Wavefront .obj file) inside [0,1]^3, solved with the Finite
           Cell Method, using Algoim combined with recursive octree
           subdivision of the cut cells for the immersed quadrature.

    This is a copy of poisson2_nitsche_immersed_geometry_example.cpp with a
    single addition: the cut-cell volume quadrature. Instead of calling Algoim
    once on the whole cut cell, the cut cell (as classified by
    gsImplicitTrimmedDomain<3>) is recursively subdivided into an octree of
    depth --quadDepth, and Algoim (+ the uniform-grid fallback) is applied on
    each leaf box. This gives more reliable / cheaper immersed quadrature on
    geometrically difficult cut cells.

      - --quadDepth 0 : exactly the original behaviour (Algoim on the full cell).
      - --quadDepth 1 : each cut cell is split into 8 leaf boxes.
      - --quadDepth 2 : each cut cell is split into 8^2 leaf boxes. etc.

    This is adaptive QUADRATURE only. The B-spline basis, the trimmed-domain
    classification, the solver, and the Poisson problem are unchanged. Only
    cut cells (beginBdr(boundary::none)) are subdivided.

    Usage:
    ./bin/poisson2_nitsche_immersed_geometry_octree_example -f spot.obj -r 2 -g 1e3 --quadDepth 1 --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include "gsMeshLevelSet.h"

#include <cmath>
#include <functional>
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
        "/Users/lucasventavinuela/gismo_gmsh/optional/gsGmsh/filedata/spot.obj";
    std::string out = gsFileManager::getCanonicRepresentation(
        gsFileManager::getPath(std::string(__FILE__), true) + "../output_poisson_immersed3d");
    index_t numRefine  = 3;
    index_t numElevate = 0;
    real_t  fill       = 0.9;  // fraction of [0,1]^3 spanned by the largest extent
    real_t  gamma      = 1e3;  // Nitsche penalty parameter (scaled by 1/h)
    index_t quadDepth  = 0;    // octree subdivision depth for cut-cell quadrature
    bool    plot       = false;
    bool    noNitsche  = false;

    gsCmdLine cmd("Poisson equation on an immersed 3D .obj domain via the Finite Cell Method (octree cut-cell quadrature).");
    cmd.addString("f", "file",            "Wavefront .obj mesh file",            filename);
    cmd.addString("o", "output",          "Output folder for ParaView files",    out);
    cmd.addInt   ("r", "uniformRefine",   "Number of uniform refinement steps",  numRefine);
    cmd.addInt   ("e", "degreeElevation", "Number of degree elevation steps",    numElevate);
    cmd.addReal  ("",  "fill",            "Fill fraction of [0,1]^3 (0..1)",      fill);
    cmd.addReal  ("g", "gamma",           "Nitsche penalty parameter",           gamma);
    cmd.addInt   ("",  "quadDepth",       "Octree subdivision depth for cut-cell quadrature", quadDepth);
    cmd.addSwitch("plot",      "Create ParaView output",                        plot);
    cmd.addSwitch("noNitsche", "Disable the weak/Nitsche BC (interior forcing only)", noNitsche);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    out = gsFileManager::getCanonicRepresentation(out);

    // -------------------------------------------------------------------------
    // 1. Load the triangle mesh
    // -------------------------------------------------------------------------
    std::vector<double> verts;
    std::vector<int>    tris;
    if (!loadObjMesh(filename, verts, tris))
    {
        gsWarn << "Failed to read a triangle mesh from: " << filename << "\n";
        return EXIT_FAILURE;
    }
    const std::string outputStem = fileStem(filename);
    const std::size_t nVert = verts.size() / 3;
    const std::size_t nTri  = tris.size()  / 3;
    gsInfo << "Loaded mesh '" << filename << "': "
           << nVert << " vertices, " << nTri << " triangles.\n";

    // -------------------------------------------------------------------------
    // 2. Rescale / center the mesh into the unit parametric box [0,1]^3
    //    mapped = (p - center) * scale + 0.5,  scale = fill / maxExtent
    // -------------------------------------------------------------------------
    Vec3 lo = {{ verts[0], verts[1], verts[2] }};
    Vec3 hi = lo;
    for (std::size_t i = 0; i < nVert; ++i)
        for (int d = 0; d < 3; ++d)
        {
            lo[d] = std::min(lo[d], verts[3 * i + d]);
            hi[d] = std::max(hi[d], verts[3 * i + d]);
        }
    const Vec3 center = {{ 0.5 * (lo[0] + hi[0]),
                           0.5 * (lo[1] + hi[1]),
                           0.5 * (lo[2] + hi[2]) }};
    const double extent = std::max(hi[0] - lo[0],
                          std::max(hi[1] - lo[1], hi[2] - lo[2]));
    const double scale  = fill / extent;
    for (std::size_t i = 0; i < nVert; ++i)
        for (int d = 0; d < 3; ++d)
            verts[3 * i + d] = (verts[3 * i + d] - center[d]) * scale + 0.5;

    // -------------------------------------------------------------------------
    // 3. Analytic level set on the rescaled mesh
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(3, 2);
    bbox.col(0).setZero();
    bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(verts, tris, bbox);

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
    gsInfo << "Cut-cell quadrature octree depth: " << quadDepth << "\n";
    gsInfo << "Solving the Poisson equation on the immersed domain in [0,1]^3.\n";
    if (noNitsche)
        gsInfo << "Weak/Nitsche BC disabled: only the interior forcing is imposed.\n";
    gsInfo << "\nDoFs (dot1=assembled, dot2=solved, dot3=error): ";

    // Number of uniform sample points per direction used as fallback when
    // Algoim returns 0 nodes (re-entrant / multi-crossing cut cells).
    const index_t nFallback = 4;
    const int D = 3;

    gsVector<real_t> l2err(numRefine + 1), h1err(numRefine + 1);

    // ParaView collections spanning all refinement levels: each category is
    // written into its own subfolder and gathered by a single .pvd, where the
    // refinement level r plays the role of a time step.
    std::unique_ptr<gsParaviewCollection> colInterior, colCut, colAll,
                                          colBg, colSol, colExact;
    if (plot)
    {
        colInterior.reset(new gsParaviewCollection(out + "/points_interior/interior"));
        colCut     .reset(new gsParaviewCollection(out + "/points_cutcells/cutcells"));
        colAll     .reset(new gsParaviewCollection(out + "/points_all/all"));
        colBg      .reset(new gsParaviewCollection(out + "/background/background"));
        colSol     .reset(new gsParaviewCollection(out + "/solution/solution"));
        colExact   .reset(new gsParaviewCollection(out + "/exact/exact"));
    }

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

        // ---------------------------------------------------------------
        // Per-leaf cut-cell volume quadrature = Algoim rule on a box, with
        // the same uniform-grid fallback as the volume example when Algoim
        // returns zero nodes (re-entrant / multi-crossing cells).
        // ---------------------------------------------------------------
        auto leafCutCellPoints = [&](const gsVector<real_t> & lo_c, const gsVector<real_t> & hi_c,
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

        // Append a (pts,wts) block to the accumulators.
        auto appendPoints = [&](gsMatrix<real_t> & ptsAll, gsVector<real_t> & wtsAll,
                                const gsMatrix<real_t> & ptsNew, const gsVector<real_t> & wtsNew)
        {
            if (ptsNew.cols() == 0) return;
            const index_t oldCols = ptsAll.cols();
            ptsAll.conservativeResize(3, oldCols + ptsNew.cols());
            ptsAll.block(0, oldCols, 3, ptsNew.cols()) = ptsNew;

            const index_t oldSize = wtsAll.size();
            wtsAll.conservativeResize(oldSize + wtsNew.size());
            wtsAll.segment(oldSize, wtsNew.size()) = wtsNew;
        };

        // Recursive octree subdivision of a cut cell. At depth==quadDepth the
        // per-leaf Algoim(+fallback) rule is applied; otherwise the box is
        // split at its midpoint into 8 children.
        std::function<void(const gsVector<real_t>&, const gsVector<real_t>&, index_t,
                           gsMatrix<real_t>&, gsVector<real_t>&)> collectCutCellPoints;
        collectCutCellPoints = [&](const gsVector<real_t> & lo_c, const gsVector<real_t> & hi_c,
                                   index_t depth, gsMatrix<real_t> & ptsAll, gsVector<real_t> & wtsAll)
        {
            if (depth >= quadDepth)
            {
                gsMatrix<real_t> ptsLeaf; gsVector<real_t> wtsLeaf;
                leafCutCellPoints(lo_c, hi_c, ptsLeaf, wtsLeaf);
                appendPoints(ptsAll, wtsAll, ptsLeaf, wtsLeaf);
                return;
            }

            const gsVector<real_t> mid = (real_t)0.5 * (lo_c + hi_c);
            gsVector<real_t> childLo(3), childHi(3);
            for (int ix = 0; ix < 2; ++ix)
            for (int iy = 0; iy < 2; ++iy)
            for (int iz = 0; iz < 2; ++iz)
            {
                childLo[0] = ix == 0 ? lo_c[0] : mid[0];
                childHi[0] = ix == 0 ? mid[0] : hi_c[0];
                childLo[1] = iy == 0 ? lo_c[1] : mid[1];
                childHi[1] = iy == 0 ? mid[1] : hi_c[1];
                childLo[2] = iz == 0 ? lo_c[2] : mid[2];
                childHi[2] = iz == 0 ? mid[2] : hi_c[2];
                collectCutCellPoints(childLo, childHi, depth + 1, ptsAll, wtsAll);
            }
        };

        // Cut-cell quadrature entry point: octree over the cut cell, Algoim on
        // the leaves. For quadDepth==0 this is one Algoim call on the full cell
        // (identical to the non-octree example).
        auto cutCellPoints = [&](const gsVector<real_t> & lo_c, const gsVector<real_t> & hi_c,
                                 gsMatrix<real_t> & pts, gsVector<real_t> & wts)
        {
            pts.resize(3, 0);
            wts.resize(0);
            collectCutCellPoints(lo_c, hi_c, 0, pts, wts);
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

        index_t nQuad = 0;

        // int_Omega grad(u).grad(v)  and  int_Omega f.v,  Omega = {phi<0}.
        for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addVolume(pts_tmp, wts_tmp);
            nQuad += pts_tmp.cols();

        }
        for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
        {
            cutCellPoints(it.lowerCorner(), it.upperCorner(), pts_tmp, wts_tmp);
            addVolume(pts_tmp, wts_tmp);
            nQuad += pts_tmp.cols();
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
        // 7. ParaView output, one subfolder + .pvd per category:
        //     - interior / cut-cell / all quadrature points (4 x N clouds),
        //     - the background mesh,
        //     - the numerical and exact solution as point clouds sampled at
        //       the quadrature points inside {phi<0} (4 x N: xyz + value).
        //   The refinement level r is used as the .pvd time step.
        // ---------------------------------------------------------------
        if (plot)
        {
            const std::string rs = std::to_string(r);
            gsMatrix<real_t> pts, phys, quadCut(4, 0), quadIn(4, 0), quadAll(4, 0),
                             solNum(4, 0), solEx(4, 0);
            gsVector<real_t> wts;

            // Append physical points + a scalar (weight or value) as 4 x N.
            auto append = [&](gsMatrix<real_t> & acc, const gsMatrix<real_t> & phys,
                              const gsVector<real_t> & scalar)
            {
                const index_t c = acc.cols();
                acc.conservativeResize(4, c + phys.cols());
                acc.block(0, c, 3, phys.cols()) = phys;
                acc.row(3).segment(c, phys.cols()) = scalar.transpose();
            };

            // Evaluate u_h and u_exact at parametric pts and append them as
            // point clouds (value in row 3).
            auto appendSolution = [&](const gsMatrix<real_t> & pts, const gsMatrix<real_t> & phys)
            {
                gsMatrix<index_t> act;
                gsMatrix<real_t>  bv, uEx;
                tbsPtr->active_into(pts, act);
                tbsPtr->eval_into  (pts, bv);
                u_exact.eval_into  (pts, uEx);
                gsVector<real_t> uh(pts.cols());
                for (index_t q = 0; q < pts.cols(); ++q)
                {
                    real_t v = 0;
                    for (index_t i = 0; i < act.rows(); ++i)
                        v += bv(i,q) * solVector(act(i,q));
                    uh[q] = v;
                }
                append(solNum, phys, uh);
                append(solEx,  phys, uEx.row(0).transpose());
            };

            // Cut-cell quadrature points (also collected into "all").
            for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
            {
                cutCellPoints(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0) continue;
                mp.patch(0).eval_into(pts, phys);
                append(quadCut, phys, wts);
                append(quadAll, phys, wts);
                appendSolution(pts, phys);
            }
            // Interior (fully inside {phi<0}) Gauss quadrature points.
            for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
            {
                gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0) continue;
                mp.patch(0).eval_into(pts, phys);
                append(quadIn, phys, wts);
                append(quadAll, phys, wts);
                appendSolution(pts, phys);
            }
            if (quadCut.cols() > 0)
            {
                gsWriteParaviewPoints(quadCut, out + "/points_cutcells/cutcells_r" + rs);
                colCut->addPart("cutcells_r" + rs + ".vtp", r);
            }
            if (quadIn.cols() > 0)
            {
                gsWriteParaviewPoints(quadIn, out + "/points_interior/interior_r" + rs);
                colInterior->addPart("interior_r" + rs + ".vtp", r);
            }
            if (quadAll.cols() > 0)
            {
                gsWriteParaviewPoints(quadAll, out + "/points_all/all_r" + rs);
                colAll->addPart("all_r" + rs + ".vtp", r);
            }

            // Background mesh at every refinement level.
            gsMesh<real_t> bgMesh(dbasis.basis(0));
            gsWriteParaview(bgMesh, out + "/background/background_r" + rs);
            colBg->addPart("background_r" + rs + ".vtp", r);

            // Numerical and exact solution as point clouds.
            if (solNum.cols() > 0)
            {
                gsWriteParaviewPoints(solNum, out + "/solution/solution_r" + rs);
                colSol->addPart("solution_r" + rs + ".vtp", r);
                gsWriteParaviewPoints(solEx, out + "/exact/exact_r" + rs);
                colExact->addPart("exact_r" + rs + ".vtp", r);
            }
        }

        gsInfo << " (h=" << hmax << ", gamma/h=" << gamma/hmax
               << ", quadPts=" << nQuad << ")\n";
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
    {
        // Immersed geometry (the input triangle mesh), written once.
        gsMesh<real_t> geomMesh;
        for (std::size_t i = 0; i < nVert; ++i)
            geomMesh.addVertex(verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]);
        for (std::size_t t = 0; t < nTri; ++t)
            geomMesh.addFace(tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]);
        gsWriteParaview(geomMesh, out + "/" + outputStem + "_geometry");

        colInterior->save(); colCut->save(); colAll->save();
        colBg->save(); colSol->save(); colExact->save();

        gsInfo << "\nParaView files written to: " << out << "\n";
    }

    return EXIT_SUCCESS;
}