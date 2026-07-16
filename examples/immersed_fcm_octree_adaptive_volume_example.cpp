/** @file immersed_fcm_octree_adaptive_volume_example.cpp

    @brief Volume of a 3D triangle mesh via the Finite Cell Method with
           adaptive octree cut-cell quadrature.

    Cell classification uses the Separating Axis Theorem (SAT) against the
    actual triangle mesh instead of signed-distance point sampling.  This
    prevents thin features (e.g. cow ears) from being silently dropped.

    Classification logic (background cell / octree sub-cell):
      SAT hits a triangle         => CUT    -> recurse / Algoim at max depth
      no hit, phi(center) < 0    => INSIDE  -> add exact box volume
      no hit, phi(center) >= 0   => OUTSIDE -> skip

    Performance:
      * Per-cell triangle candidate lists are pre-filtered by AABB overlap,
        so SAT only runs against the small subset that could touch each cell.
      * Fully-inside sub-cells contribute their exact volume (no quadrature).
      * Cut cells are processed in parallel (OpenMP).

    Outputs (same folder structure as the Poisson adaptive example):
      out/points_interior/  interior box centers (+ .pvd)
      out/points_cutcells/  Algoim cut points    (+ .pvd)
      out/points_all/       union of the above   (+ .pvd)
      out/background/       background mesh per level (+ .pvd)
      out/<stem>_geometry.* the immersed triangle mesh (once)
      out/results/volume.txt computed volume per refinement level

    Usage:
    ./bin/immersed_fcm_octree_adaptive_volume_example -f spot.obj -r 3 --quadDepth 2 --plot

    This file is part of the G+Smo library.
    Mozilla Public License, v. 2.0.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include "gsGmsh/gsMeshLevelSet.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace gismo;

// Per-triangle AABB (pre-computed once; used to filter SAT candidates per cell).
struct TriBBox { double lo[3], hi[3]; };

// Free helpers (definitions at the end of this file).
static bool triBoxSAT  (const Vec3& a, const Vec3& b, const Vec3& c,
                        const gsVector<real_t>& lo, const gsVector<real_t>& hi);
static bool boxHitsMesh(const std::vector<double>& verts, const std::vector<int>& tris,
                        const std::vector<int>& cIdx,
                        const gsVector<real_t>& lo, const gsVector<real_t>& hi);
static int  boxClassSAT(const std::vector<double>& verts, const std::vector<int>& tris,
                        const std::vector<int>& cIdx, const gsFunction<real_t>& phi,
                        const gsVector<real_t>& lo, const gsVector<real_t>& hi);
static void leafCutCellPoints(const gsAlgoimGenericRule<real_t>& algoim,
                              const gsFunction<real_t>& phi, index_t nFallback,
                              const gsVector<real_t>& lo, const gsVector<real_t>& hi,
                              gsMatrix<real_t>& pts, gsVector<real_t>& wts);
static void collectAdaptive(const gsAlgoimGenericRule<real_t>& algoim,
                            const gsFunction<real_t>& phi,
                            const std::vector<double>& verts, const std::vector<int>& tris,
                            const std::vector<int>& cIdx,
                            index_t quadDepth, index_t nFallback,
                            const gsVector<real_t>& lo, const gsVector<real_t>& hi,
                            index_t depth,
                            gsMatrix<real_t>& insCtr, gsVector<real_t>& insVol,
                            gsMatrix<real_t>& cutP,   gsVector<real_t>& cutW);
static void appendCloud(gsMatrix<real_t>& cloud, const gsMatrix<real_t>& phys,
                        const gsVector<real_t>& scalar);

int main(int argc, char * argv[])
{
    std::string filename =
        "/Users/lucasventavinuela/gismo_gmsh/optional/gsGmsh/filedata/spot.obj";
    std::string out = "output_volume_adaptive";
    index_t numRefine  = 3;
    real_t  fill       = 0.9;
    index_t quadDepth  = 0;
    std::string indicator = "integralChange";
    real_t  indicatorTol = 1e-2;
    bool    plot       = false;

    gsCmdLine cmd("Volume of a 3D .obj mesh via FCM + adaptive octree (SAT classification).");
    cmd.addString("f", "file",          "Wavefront .obj mesh file",           filename);
    cmd.addString("o", "output",        "Output folder for ParaView files",   out);
    cmd.addInt   ("r", "uniformRefine", "Number of uniform refinement steps", numRefine);
    cmd.addReal  ("",  "fill",          "Fill fraction of [0,1]^3 (0..1)",    fill);
    cmd.addInt   ("",  "quadDepth",     "Adaptive octree depth for cut cells", quadDepth);
    cmd.addString("",  "indicator",     "Adaptive indicator: uniform, fallback or "
                                        "integralChange",                      indicator);
    cmd.addReal  ("",  "indicatorTol",  "integralChange acceptance tolerance", indicatorTol);
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
    // 2. Pre-compute per-triangle AABBs (done once; filter candidates per cell)
    // -------------------------------------------------------------------------
    std::vector<TriBBox> triBBoxes(nTri);
    for (std::size_t t = 0; t < nTri; ++t)
    {
        const int ia = tris[3*t], ib = tris[3*t+1], ic = tris[3*t+2];
        for (int d = 0; d < 3; ++d)
        {
            triBBoxes[t].lo[d] = std::min({verts[3*ia+d], verts[3*ib+d], verts[3*ic+d]});
            triBBoxes[t].hi[d] = std::max({verts[3*ia+d], verts[3*ib+d], verts[3*ic+d]});
        }
    }

    // -------------------------------------------------------------------------
    // 3. Level set and background B-spline box [0,1]^3
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(3,2); bbox.col(0).setZero(); bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(verts, tris, bbox);

    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsMultiBasis<> dbasis(mp, true);

    gsFileManager::mkdir(out);

    gsInfo << "\nAdaptive octree depth : " << quadDepth << "\n";
    gsInfo << std::setw(8)  << "refine"
           << std::setw(10) << "cells/dim"
           << std::setw(18) << "vol (param)"
           << std::setw(18) << "vol (phys)"
            << std::setw(12) << "quadPts"
            << std::setw(12) << "bdryPts"
            << std::setw(12) << "subBoxes" << "\n";

    const double invScale3 = 1.0 / (scale*scale*scale);
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
    std::vector<index_t> cellsHist, quadHist, boundaryQuadHist, interiorSubQuadHist;
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

        gsOptionList adaptiveOptions = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
        adaptiveOptions.setInt("maxDepth", quadDepth);
        adaptiveOptions.setInt("nFallback", 4);
        adaptiveOptions.setString("indicator", indicator);
        adaptiveOptions.setReal("indicatorTol", indicatorTol);
        adaptiveOptions.setSwitch("analyticInterior", true);
        gsAlgoimAdaptiveRule<real_t> adaptiveRule(impl_fun, *tbsPtr, adaptiveOptions);

        // Element break-points (unique knot values per direction).
        const std::vector<real_t> bx = tbsPtr->knots(0).breaks();
        const std::vector<real_t> by = tbsPtr->knots(1).breaks();
        const std::vector<real_t> bz = tbsPtr->knots(2).breaks();

        // Classify every background element with SAT + phi.
        real_t volInterior = 0;
        gsMatrix<real_t> intCtr(3,0);
        gsVector<real_t> intVol;
        std::vector<gsVector<real_t>> cutLo, cutHi;
        std::vector<std::vector<int>> cutTris;   // per-cut-cell candidate triangles

        for (std::size_t i = 0; i+1 < bx.size(); ++i)
        for (std::size_t j = 0; j+1 < by.size(); ++j)
        for (std::size_t k = 0; k+1 < bz.size(); ++k)
        {
            gsVector<real_t> lo(3), hi(3);
            lo << bx[i], by[j], bz[k];
            hi << bx[i+1], by[j+1], bz[k+1];

            // AABB pre-filter: only triangles whose bbox overlaps this cell.
            std::vector<int> cIdx;
            cIdx.reserve(32);
            for (std::size_t t = 0; t < nTri; ++t)
                if (triBBoxes[t].hi[0] >= lo[0] && triBBoxes[t].lo[0] <= hi[0] &&
                    triBBoxes[t].hi[1] >= lo[1] && triBBoxes[t].lo[1] <= hi[1] &&
                    triBBoxes[t].hi[2] >= lo[2] && triBBoxes[t].lo[2] <= hi[2])
                    cIdx.push_back(static_cast<int>(t));

            const int cls = boxClassSAT(verts, tris, cIdx, impl_fun, lo, hi);
            if (cls > 0) continue;   // outside: skip
            if (cls < 0)             // inside: exact volume
            {
                const real_t v = (hi - lo).prod();
                volInterior += v;
                if (plot)
                {
                    const index_t c = intCtr.cols();
                    intCtr.conservativeResize(3, c + 1);
                    intCtr.col(c) = 0.5*(lo + hi);
                    intVol.conservativeResize(c + 1);
                    intVol[c] = v;
                }
                continue;
            }
            cutLo.push_back(lo);
            cutHi.push_back(hi);
            cutTris.push_back(std::move(cIdx));
        }
        const int nCut = static_cast<int>(cutLo.size());

        // Global point clouds (only touched under a critical section when plotting).
        gsMatrix<real_t> quadIn(4,0), quadCut(4,0), quadAll(4,0);

        real_t  volCut = 0;
        index_t nQuad  = 0;
        index_t nBoundaryQuad = 0;
        index_t nInteriorSubQuad = 0;


#       pragma omp parallel for schedule(dynamic) reduction(+:volCut,nQuad,nBoundaryQuad,nInteriorSubQuad)
        for (int i = 0; i < nCut; ++i)
        {
            gsMatrix<real_t> insCtr(3,0), cutP(3,0);
            gsVector<real_t> insVol, cutW;
            adaptiveRule.mapToSeparated(
                cutLo[i], cutHi[i], insCtr, insVol, cutP, cutW,
                [&](const gsVector<real_t> & childLo,
                    const gsVector<real_t> & childHi) -> int
                {
                    return boxClassSAT(verts, tris, cutTris[i], impl_fun,
                                       childLo, childHi);
                });

            volCut += insVol.sum() + cutW.sum();
            nQuad  += insCtr.cols() + cutP.cols();

            // Inside your OpenMP loop, accumulate them separately:
            nBoundaryQuad    += cutP.cols();
            nInteriorSubQuad += insCtr.cols();

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
             << std::setw(12) << nQuad
             << std::setw(12) << nBoundaryQuad
             << std::setw(12) << nInteriorSubQuad << "\n";

        cellsHist.push_back(cellsPerDim);
        quadHist .push_back(nQuad);
         boundaryQuadHist.push_back(nBoundaryQuad);
         interiorSubQuadHist.push_back(nInteriorSubQuad);
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
        txt << "# rescaling factor (param/phys): " << scale << "\n\n";
        txt << std::setw(8)  << "refine" << std::setw(10) << "cells/dim"
            << std::setw(20) << "vol(param)" << std::setw(20) << "vol(phys)"
            << std::setw(12) << "quadPts"
            << std::setw(12) << "bdryPts"
            << std::setw(12) << "subBoxes" << "\n";
        for (std::size_t k = 0; k < volParHist.size(); ++k)
            txt << std::setw(8)  << k << std::setw(10) << cellsHist[k]
                << std::setw(20) << std::fixed << std::setprecision(8) << volParHist[k]
                << std::setw(20) << std::fixed << std::setprecision(8) << volPhyHist[k]
                << std::setw(12) << quadHist[k]
                << std::setw(12) << boundaryQuadHist[k]
                << std::setw(12) << interiorSubQuadHist[k] << "\n";
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
//  Helper implementations
// =============================================================================

// Single SAT axis test: project triangle vertices + box radius onto axis ax.
// Returns true if this axis separates (=> no intersection along this axis).
static inline bool satAxis(const Vec3& ax,
                           const Vec3& v0, const Vec3& v1, const Vec3& v2,
                           const Vec3& h)
{
    const double p0 = dot(ax,v0), p1 = dot(ax,v1), p2 = dot(ax,v2);
    const double r  = h[0]*std::abs(ax[0]) + h[1]*std::abs(ax[1]) + h[2]*std::abs(ax[2]);
    return std::min({p0,p1,p2}) > r || std::max({p0,p1,p2}) < -r;
}

// Möller (1997) triangle–AABB intersection test using 13 SAT axes.
static bool triBoxSAT(const Vec3& a, const Vec3& b, const Vec3& c,
                      const gsVector<real_t>& lo, const gsVector<real_t>& hi)
{
    const Vec3 bc = {{0.5*(lo[0]+hi[0]), 0.5*(lo[1]+hi[1]), 0.5*(lo[2]+hi[2])}};
    const Vec3 h  = {{0.5*(hi[0]-lo[0]), 0.5*(hi[1]-lo[1]), 0.5*(hi[2]-lo[2])}};
    const Vec3 va = sub(a,bc), vb = sub(b,bc), vc = sub(c,bc);
    const Vec3 e0 = sub(vb,va), e1 = sub(vc,vb), e2 = sub(va,vc);

    // 3 AABB face normals.
    for (int d = 0; d < 3; ++d)
    { Vec3 ax = {{0,0,0}}; ax[d] = 1.0; if (satAxis(ax,va,vb,vc,h)) return false; }

    // Triangle face normal.
    if (satAxis(cross(e0,e1), va,vb,vc,h)) return false;

    // 9 cross products: each box axis × each triangle edge.
    const Vec3 baxes[3] = {{{1,0,0}}, {{0,1,0}}, {{0,0,1}}};
    const Vec3 edges[3] = {e0, e1, e2};
    for (int ai = 0; ai < 3; ++ai)
        for (int ej = 0; ej < 3; ++ej)
            if (satAxis(cross(baxes[ai],edges[ej]), va,vb,vc,h)) return false;

    return true;  // no separating axis found -> intersection
}

// Returns true if box [lo,hi] intersects any triangle in the candidate list cIdx.
static bool boxHitsMesh(const std::vector<double>& verts, const std::vector<int>& tris,
                        const std::vector<int>& cIdx,
                        const gsVector<real_t>& lo, const gsVector<real_t>& hi)
{
    for (int ti : cIdx)
    {
        const int ia = tris[3*ti], ib = tris[3*ti+1], ic = tris[3*ti+2];
        const Vec3 a = {{verts[3*ia],verts[3*ia+1],verts[3*ia+2]}};
        const Vec3 b = {{verts[3*ib],verts[3*ib+1],verts[3*ib+2]}};
        const Vec3 c = {{verts[3*ic],verts[3*ic+1],verts[3*ic+2]}};
        if (triBoxSAT(a,b,c,lo,hi)) return true;
    }
    return false;
}

// Hybrid cell classification (SAT + Lipschitz radius test):
//   SAT hit                         =>  0  (cut: contains a surface triangle)
//   no SAT hit, phi(ctr) >  radius  => +1  (confidently outside)
//   no SAT hit, phi(ctr) < -radius  => -1  (confidently inside -> exact volume)
//   no SAT hit, |phi(ctr)| <= radius =>  0  (near boundary, no direct hit -> Algoim)
//
// Why both tests?
//   SAT alone: cells just inside the boundary (no triangle) get phi(ctr)<0 -> full
//              box volume added, but part of the cell may be outside -> overestimates.
//   Radius alone: misses thin-feature cells whose surface triangles clip a corner
//              of the cell without the center being close to any triangle.
//   Combined: SAT catches thin features; radius keeps near-boundary cells honest.
static int boxClassSAT(const std::vector<double>& verts, const std::vector<int>& tris,
                       const std::vector<int>& cIdx, const gsFunction<real_t>& phi,
                       const gsVector<real_t>& lo, const gsVector<real_t>& hi)
{
    if (boxHitsMesh(verts, tris, cIdx, lo, hi)) return 0;   // SAT: definitely cut
    gsMatrix<real_t> c(3,1), v;
    c.col(0) = 0.5*(lo + hi);
    phi.eval_into(c, v);
    const real_t radius = (real_t)0.5 * (hi - lo).norm();  // half-diagonal
    if (v(0,0) >  radius) return  1;   // confidently outside
    if (v(0,0) < -radius) return -1;   // confidently inside
    return 0;   // near boundary but no direct triangle hit -> send to Algoim
}

// Algoim {phi<0} volumetric rule on a leaf box, with uniform midpoint fallback
// when Algoim returns zero nodes (re-entrant / multi-crossing cells).
static void leafCutCellPoints(const gsAlgoimGenericRule<real_t>& algoim,
                              const gsFunction<real_t>& phi, index_t nFallback,
                              const gsVector<real_t>& lo, const gsVector<real_t>& hi,
                              gsMatrix<real_t>& pts, gsVector<real_t>& wts)
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
    gsMatrix<real_t> vals; phi.eval_into(spts, vals);
    std::vector<index_t> ins;
    for (index_t ci = 0; ci < nTotal; ++ci)
        if (vals(0,ci) < 0) ins.push_back(ci);
    if (ins.empty()) { pts.resize(3,0); wts.resize(0); return; }
    const real_t w = (hi - lo).prod() / static_cast<real_t>(nTotal);
    pts.resize(3, static_cast<index_t>(ins.size()));
    wts.resize(   static_cast<index_t>(ins.size()));
    for (std::size_t ci = 0; ci < ins.size(); ++ci)
    { pts.col(ci) = spts.col(ins[ci]); wts[ci] = w; }
}

// Adaptive octree over one background cut cell.
// Uses SAT-based classification so thin features are never silently pruned.
// cIdx: candidate triangles pre-filtered by the parent cell's AABB.
static void collectAdaptive(const gsAlgoimGenericRule<real_t>& algoim,
                            const gsFunction<real_t>& phi,
                            const std::vector<double>& verts, const std::vector<int>& tris,
                            const std::vector<int>& cIdx,
                            index_t quadDepth, index_t nFallback,
                            const gsVector<real_t>& lo, const gsVector<real_t>& hi,
                            index_t depth,
                            gsMatrix<real_t>& insCtr, gsVector<real_t>& insVol,
                            gsMatrix<real_t>& cutP,   gsVector<real_t>& cutW)
{
    if (depth >= quadDepth)
    {
        gsMatrix<real_t> p; gsVector<real_t> w;
        leafCutCellPoints(algoim, phi, nFallback, lo, hi, p, w);
        if (p.cols() > 0)
        {
            const index_t oc = cutP.cols();
            cutP.conservativeResize(3, oc + p.cols());
            cutP.block(0, oc, 3, p.cols()) = p;
            cutW.conservativeResize(oc + w.size());
            cutW.segment(oc, w.size()) = w;
        }
        return;
    }

    const gsVector<real_t> mid = 0.5*(lo + hi);
    gsVector<real_t> clo(3), chi(3);
    for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
    for (int iz = 0; iz < 2; ++iz)
    {
        clo[0] = ix ? mid[0] : lo[0];  chi[0] = ix ? hi[0] : mid[0];
        clo[1] = iy ? mid[1] : lo[1];  chi[1] = iy ? hi[1] : mid[1];
        clo[2] = iz ? mid[2] : lo[2];  chi[2] = iz ? hi[2] : mid[2];

        // SAT-based: cannot falsely prune cells containing thin triangles.
        const int cls = boxClassSAT(verts, tris, cIdx, phi, clo, chi);
        if (cls > 0) continue;   // outside: skip
        if (cls < 0)             // inside: exact volume
        {
            const index_t c = insCtr.cols();
            insCtr.conservativeResize(3, c + 1);
            insCtr.col(c) = 0.5*(clo + chi);
            insVol.conservativeResize(c + 1);
            insVol[c] = (chi - clo).prod();
            continue;
        }
        collectAdaptive(algoim, phi, verts, tris, cIdx, quadDepth, nFallback,
                        clo, chi, depth+1, insCtr, insVol, cutP, cutW);
    }
}

// Append physical points + scalar values as a 4×N point cloud.
static void appendCloud(gsMatrix<real_t>& cloud, const gsMatrix<real_t>& phys,
                        const gsVector<real_t>& scalar)
{
    if (phys.cols() == 0) return;
    const index_t c = cloud.cols();
    cloud.conservativeResize(4, c + phys.cols());
    cloud.block(0, c, 3, phys.cols()) = phys;
    cloud.row(3).segment(c, phys.cols()) = scalar.transpose();
}
