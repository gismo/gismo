/** @file volume_cow_fcm_octree_adaptive_example.cpp

    @brief Volume of a 3D triangle mesh via the Finite Cell Method with
           adaptive octree cut-cell quadrature (same strategy as
           poisson2_nitsche_immersed_geometry_octree_adaptive_example.cpp).

    Each background cut cell is subdivided adaptively using the ~1-Lipschitz
    signed-distance test (radius = half the box diagonal):
      phi(center) < -radius  => child fully INSIDE  -> exact box volume
      phi(center) >  radius  => child fully OUTSIDE -> skipped
      otherwise              => child still CUT     -> recurse / Algoim at max depth

    --quadDepth 0 reproduces the original single-Algoim-per-cell behaviour.

    Performance notes:
      * Fully-interior sub-cells contribute their exact box volume (no Gauss
        points needed to integrate the constant 1) -> far fewer evaluations.
      * The cut cells are collected once and integrated in parallel (OpenMP);
        Algoim and the signed-distance evaluations are const / read-only, so
        the loop is thread-safe.

    Outputs (mirrors the Poisson adaptive example):
      out/points_interior/  interior quadrature/box points (+ .pvd)
      out/points_cutcells/  Algoim cut-cell points        (+ .pvd)
      out/points_all/       union of the above            (+ .pvd)
      out/background/        background mesh per level     (+ .pvd)
      out/<stem>_geometry.*  the immersed triangle mesh (once)
      out/results/volume.txt the computed volume per refinement level

    Usage:
    ./bin/volume_cow_fcm_octree_adaptive_example -f spot.obj -r 3 --quadDepth 2 --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include "gsMeshLevelSet.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

// Free helpers (definitions at the end of this file).
static int  boxClass(const gsFunction<real_t> & phi,
                     const gsVector<real_t> & lo, const gsVector<real_t> & hi);
static void leafCutCellPoints(const gsAlgoimGenericRule<real_t> & algoim,
                              const gsFunction<real_t> & phi, index_t nFallback,
                              const gsVector<real_t> & lo, const gsVector<real_t> & hi,
                              gsMatrix<real_t> & pts, gsVector<real_t> & wts);
static void collectAdaptive(const gsAlgoimGenericRule<real_t> & algoim,
                            const gsFunction<real_t> & phi,
                            index_t quadDepth, index_t nFallback,
                            const gsVector<real_t> & lo, const gsVector<real_t> & hi,
                            index_t depth,
                            gsMatrix<real_t> & insCtr, gsVector<real_t> & insVol,
                            gsMatrix<real_t> & cutP,   gsVector<real_t> & cutW);
static void appendCloud(gsMatrix<real_t> & cloud, const gsMatrix<real_t> & phys,
                        const gsVector<real_t> & scalar);

int main(int argc, char * argv[])
{
    std::string filename =
        "/Users/lucasventavinuela/gismo_gmsh/optional/gsGmsh/filedata/spot.obj";
    std::string out = "output_volume_adaptive";
    index_t numRefine  = 3;
    real_t  fill       = 0.9;
    index_t quadDepth  = 0;
    index_t classifySamples = 5;
    bool    plot       = false;

    gsCmdLine cmd("Volume of a 3D .obj mesh via FCM + adaptive octree cut-cell quadrature.");
    cmd.addString("f", "file",          "Wavefront .obj mesh file",           filename);
    cmd.addString("o", "output",        "Output folder for ParaView files",   out);
    cmd.addInt   ("r", "uniformRefine", "Number of uniform refinement steps", numRefine);
    cmd.addReal  ("",  "fill",          "Fill fraction of [0,1]^3 (0..1)",    fill);
    cmd.addInt   ("",  "quadDepth",     "Adaptive octree depth for cut cells", quadDepth);
    cmd.addInt   ("c",  "classifySamples", "Cell-classification samples per direction", classifySamples);
    cmd.addSwitch("plot", "Create ParaView output",                           plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // -------------------------------------------------------------------------
    // 1. Load and rescale the mesh into [0,1]^3
    // -------------------------------------------------------------------------
    std::vector<double> verts;
    std::vector<int>    tris;
    if (!loadObjMesh(filename, verts, tris))
    { gsWarn << "Failed to read: " << filename << "\n"; return EXIT_FAILURE; }

    const std::string stem = fileStem(filename);
    const std::size_t nVert = verts.size() / 3;
    const std::size_t nTri  = tris.size()  / 3;
    gsInfo << "Loaded '" << filename << "': " << nVert << " verts, " << nTri << " tris.\n";

    Vec3 lo = {{verts[0],verts[1],verts[2]}}, hi = lo;
    for (std::size_t i = 0; i < nVert; ++i)
        for (int d = 0; d < 3; ++d)
        { lo[d] = std::min(lo[d],verts[3*i+d]); hi[d] = std::max(hi[d],verts[3*i+d]); }
    const Vec3 center = {{0.5*(lo[0]+hi[0]), 0.5*(lo[1]+hi[1]), 0.5*(lo[2]+hi[2])}};
    const double extent = std::max(hi[0]-lo[0], std::max(hi[1]-lo[1], hi[2]-lo[2]));
    const double scale  = fill / extent;
    for (std::size_t i = 0; i < nVert; ++i)
        for (int d = 0; d < 3; ++d)
            verts[3*i+d] = (verts[3*i+d] - center[d]) * scale + 0.5;
    const double bboxVolPhys = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);

    // -------------------------------------------------------------------------
    // 2. Level set and background B-spline box [0,1]^3
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(3,2); bbox.col(0).setZero(); bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(verts, tris, bbox);

    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsMultiBasis<> dbasis(mp, true);
    const index_t deg = dbasis.maxCwiseDegree();

    gsFileManager::mkdir(out);

    gsInfo << "\nAdaptive octree depth: " << quadDepth << "\n";
    gsInfo << "Cell-classification samples: " << classifySamples << " per direction\n";
    gsInfo << std::setw(8)  << "refine"
           << std::setw(10) << "cells/dim"
           << std::setw(18) << "vol (param)"
           << std::setw(18) << "vol (phys)"
           << std::setw(12) << "quadPts" << "\n";

    const double invScale3 = 1.0 / (scale*scale*scale);
    const index_t nFallback = 4;

    // ParaView collections spanning all refinement levels (one .pvd each).
    std::unique_ptr<gsParaviewCollection> colInterior, colCut, colAll, colBg;
    if (plot)
    {
        colInterior.reset(new gsParaviewCollection(out + "/points_interior/interior"));
        colCut     .reset(new gsParaviewCollection(out + "/points_cutcells/cutcells"));
        colAll     .reset(new gsParaviewCollection(out + "/points_all/all"));
        colBg      .reset(new gsParaviewCollection(out + "/background/background"));
    }

    // Per-level history for the results file.
    std::vector<index_t> cellsHist, quadHist;
    std::vector<real_t>  volParHist, volPhyHist;

    // -------------------------------------------------------------------------
    // 3. Refinement loop
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();
        gsTensorBSplineBasis<3,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<3,real_t>*>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        gsImplicitTrimmedDomain<3,real_t> tr_domain(impl_fun, *tbsPtr, classifySamples);
        gsAlgoimGenericRule<real_t> algoimRule(impl_fun, *tbsPtr);

        // Fully-interior background cells: add their exact box volume.
        real_t volInterior = 0;
        gsMatrix<real_t> intCtr(3,0);
        gsVector<real_t> intVol;
        for (auto it = tr_domain.beginInterior(); it != tr_domain.end<InteriorSign>(); ++it)
        {
            const gsVector<real_t> l = it.lowerCorner(), u = it.upperCorner();
            const real_t v = (u - l).prod();
            volInterior += v;
            if (plot)
            {
                const index_t c = intCtr.cols();
                intCtr.conservativeResize(3, c + 1);
                intCtr.col(c) = (real_t)0.5 * (l + u);
                intVol.conservativeResize(c + 1);
                intVol[c] = v;
            }
        }

        // Collect cut cells once so the expensive octree/Algoim work can run
        // over a random-access container (and in parallel).
        std::vector<gsVector<real_t>> cutLo, cutHi;
        for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
        { cutLo.push_back(it.lowerCorner()); cutHi.push_back(it.upperCorner()); }
        const int nCut = static_cast<int>(cutLo.size());

        // Global point clouds (only touched under a critical section when plotting).
        gsMatrix<real_t> quadIn(4,0), quadCut(4,0), quadAll(4,0);

        real_t  volCut = 0;
        index_t nQuad  = 0;

#       pragma omp parallel for schedule(dynamic) reduction(+:volCut,nQuad)
        for (int i = 0; i < nCut; ++i)
        {
            gsMatrix<real_t> insCtr(3,0), cutP(3,0);
            gsVector<real_t> insVol, cutW;
            collectAdaptive(algoimRule, impl_fun, quadDepth, nFallback,
                            cutLo[i], cutHi[i], 0, insCtr, insVol, cutP, cutW);

            volCut += insVol.sum() + cutW.sum();
            nQuad  += insCtr.cols() + cutP.cols();

            if (plot)
            {
                gsMatrix<real_t> physIn, physCut;
                mp.patch(0).eval_into(insCtr, physIn);
                mp.patch(0).eval_into(cutP,   physCut);
#               pragma omp critical
                {
                    appendCloud(quadIn,  physIn,  insVol);
                    appendCloud(quadCut, physCut, cutW);
                    appendCloud(quadAll, physIn,  insVol);
                    appendCloud(quadAll, physCut, cutW);
                }
            }
        }

        const real_t  volume = volInterior + volCut;
        const index_t cellsPerDim = tbsPtr->component(0).numElements();

        gsInfo << std::setw(8)  << r
               << std::setw(10) << cellsPerDim
               << std::setw(18) << std::fixed << std::setprecision(6) << volume
               << std::setw(18) << std::fixed << std::setprecision(6) << volume*invScale3
               << std::setw(12) << nQuad << "\n";

        cellsHist.push_back(cellsPerDim);
        quadHist .push_back(nQuad);
        volParHist.push_back(volume);
        volPhyHist.push_back(volume*invScale3);

        // --- ParaView output ---
        if (plot)
        {
            const std::string rs = std::to_string(r);

            // Interior background cells (one point per cell center).
            gsMatrix<real_t> physInt;
            mp.patch(0).eval_into(intCtr, physInt);
            appendCloud(quadIn,  physInt, intVol);
            appendCloud(quadAll, physInt, intVol);

            if (quadIn.cols() > 0)
            {
                gsWriteParaviewPoints(quadIn, out + "/points_interior/interior_r" + rs);
                colInterior->addPart("interior_r" + rs + ".vtp", r);
            }
            if (quadCut.cols() > 0)
            {
                gsWriteParaviewPoints(quadCut, out + "/points_cutcells/cutcells_r" + rs);
                colCut->addPart("cutcells_r" + rs + ".vtp", r);
            }
            if (quadAll.cols() > 0)
            {
                gsWriteParaviewPoints(quadAll, out + "/points_all/all_r" + rs);
                colAll->addPart("all_r" + rs + ".vtp", r);
            }

            gsMesh<real_t> bgMesh(dbasis.basis(0));
            gsWriteParaview(bgMesh, out + "/background/background_r" + rs);
            colBg->addPart("background_r" + rs + ".vtp", r);
        }
    }

    // -------------------------------------------------------------------------
    // 4. Results file
    // -------------------------------------------------------------------------
    gsFileManager::mkdir(out + "/results");
    {
        std::ofstream txt((out + "/results/volume.txt").c_str());
        txt << "# Volume of '" << filename << "' via FCM + adaptive octree\n";
        txt << "# mesh: " << nVert << " vertices, " << nTri << " triangles\n";
        txt << "# adaptive octree depth (quadDepth): " << quadDepth << "\n";
        txt << "# cell-classification samples per direction: " << classifySamples << "\n";
        txt << "# rescaling factor (param/phys): " << scale << "\n\n";
        txt << std::setw(8)  << "refine" << std::setw(10) << "cells/dim"
            << std::setw(20) << "vol(param)" << std::setw(20) << "vol(phys)"
            << std::setw(12) << "quadPts" << "\n";
        for (std::size_t k = 0; k < volParHist.size(); ++k)
            txt << std::setw(8)  << k << std::setw(10) << cellsHist[k]
                << std::setw(20) << std::fixed << std::setprecision(8) << volParHist[k]
                << std::setw(20) << std::fixed << std::setprecision(8) << volPhyHist[k]
                << std::setw(12) << quadHist[k] << "\n";
        const std::size_t last = volParHist.size() - 1;
        txt << "\n=== Result (finest level r=" << last << ") ===\n";
        txt << "Volume (parametric [0,1]^3): " << std::setprecision(8) << volParHist[last] << "\n";
        txt << "Volume (physical units):     " << std::setprecision(8) << volPhyHist[last] << "\n";
        txt << "Percentage of [0,1]^3 box:   " << std::setprecision(2) << 100.0*volParHist[last] << " %\n";
        txt << "Percentage of object bbox:   " << std::setprecision(2)
            << 100.0*volPhyHist[last]/bboxVolPhys << " %\n";
    }
    gsInfo << "\nResults written to: " << out << "/results/volume.txt\n";

    // -------------------------------------------------------------------------
    // 5. Finalize ParaView output
    // -------------------------------------------------------------------------
    if (plot)
    {
        gsMesh<real_t> geomMesh;
        for (std::size_t i = 0; i < nVert; ++i)
            geomMesh.addVertex(verts[3*i], verts[3*i+1], verts[3*i+2]);
        for (std::size_t t = 0; t < nTri; ++t)
            geomMesh.addFace(tris[3*t], tris[3*t+1], tris[3*t+2]);
        gsWriteParaview(geomMesh, out + "/" + stem + "_geometry");

        colInterior->save(); colCut->save(); colAll->save(); colBg->save();
        gsInfo << "ParaView files written to: " << out << "\n";
    }

    return EXIT_SUCCESS;
}

// =============================================================================
//  Helper definitions
// =============================================================================

// Signed-distance box classification (phi is ~1-Lipschitz):
//   +1 : box fully outside {phi<0}   (phi(center) >  radius)
//   -1 : box fully inside  {phi<0}   (phi(center) < -radius)
//    0 : box still cut / uncertain
// radius = half the box diagonal. Conservative: uncertain boxes are refined.
static int boxClass(const gsFunction<real_t> & phi,
                    const gsVector<real_t> & lo, const gsVector<real_t> & hi)
{
    gsMatrix<real_t> c(3,1), v;
    c.col(0) = (real_t)0.5 * (lo + hi);
    phi.eval_into(c, v);
    const real_t radius = (real_t)0.5 * (hi - lo).norm();
    if (v(0,0) >  radius) return  1;
    if (v(0,0) < -radius) return -1;
    return 0;
}

// Algoim {phi<0} volumetric rule on a leaf box, with a uniform-grid midpoint
// fallback when Algoim returns zero nodes (re-entrant / multi-crossing cells).
static void leafCutCellPoints(const gsAlgoimGenericRule<real_t> & algoim,
                              const gsFunction<real_t> & phi, index_t nFallback,
                              const gsVector<real_t> & lo, const gsVector<real_t> & hi,
                              gsMatrix<real_t> & pts, gsVector<real_t> & wts)
{
    algoim.mapTo(lo, hi, pts, wts);
    if (wts.size() > 0) return;

    const index_t nTotal = nFallback*nFallback*nFallback;
    gsMatrix<real_t> spts(3, nTotal);
    index_t col = 0;
    for (index_t ii = 0; ii < nFallback; ++ii)
    for (index_t jj = 0; jj < nFallback; ++jj)
    for (index_t kk = 0; kk < nFallback; ++kk, ++col)
    {
        spts(0,col) = lo[0] + (ii+0.5)/nFallback*(hi[0]-lo[0]);
        spts(1,col) = lo[1] + (jj+0.5)/nFallback*(hi[1]-lo[1]);
        spts(2,col) = lo[2] + (kk+0.5)/nFallback*(hi[2]-lo[2]);
    }
    gsMatrix<real_t> fbVals;
    phi.eval_into(spts, fbVals);
    std::vector<index_t> ins;
    for (index_t ci = 0; ci < nTotal; ++ci)
        if (fbVals(0,ci) < 0) ins.push_back(ci);
    if (ins.empty()) { pts.resize(3,0); wts.resize(0); return; }

    const real_t w = (hi - lo).prod() / static_cast<real_t>(nTotal);
    pts.resize(3, static_cast<index_t>(ins.size()));
    wts.resize(    static_cast<index_t>(ins.size()));
    for (std::size_t ci = 0; ci < ins.size(); ++ci)
    { pts.col(ci) = spts.col(ins[ci]); wts[ci] = w; }
}

// Adaptive octree over one background cut cell. Fully-inside children add their
// exact box volume (center + volume, no Gauss points); fully-outside children
// are skipped; still-cut children are subdivided until quadDepth, where Algoim
// (+ fallback) is applied. quadDepth==0 => Algoim on the whole cell.
static void collectAdaptive(const gsAlgoimGenericRule<real_t> & algoim,
                            const gsFunction<real_t> & phi,
                            index_t quadDepth, index_t nFallback,
                            const gsVector<real_t> & lo, const gsVector<real_t> & hi,
                            index_t depth,
                            gsMatrix<real_t> & insCtr, gsVector<real_t> & insVol,
                            gsMatrix<real_t> & cutP,   gsVector<real_t> & cutW)
{
    if (depth >= quadDepth)
    {
        gsMatrix<real_t> p; gsVector<real_t> w;
        leafCutCellPoints(algoim, phi, nFallback, lo, hi, p, w);
        if (p.cols() > 0)
        {
            const index_t oc = cutP.cols(), os = cutW.size();
            cutP.conservativeResize(3, oc + p.cols());
            cutP.block(0, oc, 3, p.cols()) = p;
            cutW.conservativeResize(os + w.size());
            cutW.segment(os, w.size()) = w;
        }
        return;
    }

    const gsVector<real_t> mid = (real_t)0.5 * (lo + hi);
    gsVector<real_t> clo(3), chi(3);
    for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
    for (int iz = 0; iz < 2; ++iz)
    {
        clo[0] = ix==0 ? lo[0] : mid[0]; chi[0] = ix==0 ? mid[0] : hi[0];
        clo[1] = iy==0 ? lo[1] : mid[1]; chi[1] = iy==0 ? mid[1] : hi[1];
        clo[2] = iz==0 ? lo[2] : mid[2]; chi[2] = iz==0 ? mid[2] : hi[2];

        const int cls = boxClass(phi, clo, chi);
        if (cls > 0) continue;                         // outside: skip
        if (cls < 0)                                   // inside: exact box volume
        {
            const index_t c = insCtr.cols();
            insCtr.conservativeResize(3, c + 1);
            insCtr.col(c) = (real_t)0.5 * (clo + chi);
            insVol.conservativeResize(c + 1);
            insVol[c] = (chi - clo).prod();
            continue;
        }
        collectAdaptive(algoim, phi, quadDepth, nFallback,
                        clo, chi, depth+1, insCtr, insVol, cutP, cutW);
    }
}

// Append physical points + a scalar (row 3) to a 4 x N point cloud.
static void appendCloud(gsMatrix<real_t> & cloud, const gsMatrix<real_t> & phys,
                        const gsVector<real_t> & scalar)
{
    if (phys.cols() == 0) return;
    const index_t c = cloud.cols();
    cloud.conservativeResize(4, c + phys.cols());
    cloud.block(0, c, 3, phys.cols()) = phys;
    cloud.row(3).segment(c, phys.cols()) = scalar.transpose();
}
