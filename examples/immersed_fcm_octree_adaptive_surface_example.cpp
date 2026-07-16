/** @file immersed_fcm_octree_adaptive_surface_example.cpp

    @brief Area and perimeter of the Spain coastline via the Finite Cell Method
           with adaptive quad-tree cut-cell quadrature (Algoim).

    The physical domain (Spain) is given as an *immersed* boundary: a closed
    coastline polygon read either from a Gmsh 2D mesh (.msh, MSH 4.1) or from a
    plain text file of ordered "x y" vertices.  We do NOT use the conforming
    Gmsh triangulation for FEM.  Instead we build a signed-distance level set

        phi(x,y) < 0  ->  inside  Spain
        phi(x,y) = 0  ->  on the coastline
        phi(x,y) > 0  ->  outside Spain

    and embed it in a simple Cartesian tensor-product B-spline background mesh.
    The Finite Cell Method machinery then classifies cells and integrates:

            - Area      : interior background cells (exact) + adaptive quad-tree over
                                        cut cells (interior sub-cells exact + Algoim {phi<0} leaves)
            - Perimeter : Algoim surface quadrature ({phi=0}) over cut background cells

        Output counters follow the 3D adaptive volume example:
            - quadPts : quadrature points created inside adaptive cut-cell processing
            - bdryPts : Algoim/fallback points on cut leaves
            - subBoxes: fully-inside adaptive subboxes

    This mirrors area_fcm_example.cpp but replaces the analytic level set
    (gsFunctionExpr) with the geometry-based SpainLevelSet.

    Usage:
      ./immersed_fcm_octree_adaptive_surface_example [-f spain_mesh.msh] [-r 5] [--quadDepth 2]
                          [-o output_spain_fcm] [--plot]

    Examples:
      ./build/bin/immersed_fcm_octree_adaptive_surface_example -r 5 --quadDepth 1 --plot
      ./build/bin/immersed_fcm_octree_adaptive_surface_example -f coastline_points.txt -r 6 --quadDepth 2

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Venta-Viñuela
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace gismo;

// =============================================================================
//  Closed coastline polygon: ordered vertices + edges + robust in-out test.
// =============================================================================

// A straight coastline segment between two consecutive polygon vertices.
struct Seg { real_t x0, y0, x1, y1; };

// Build the closed list of edges from ordered vertices (last joins first).
static std::vector<Seg> polygonSegments(const std::vector<std::pair<real_t,real_t> >& pts)
{
    std::vector<Seg> segs;
    const size_t n = pts.size();
    segs.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        const std::pair<real_t,real_t>& a = pts[i];
        const std::pair<real_t,real_t>& b = pts[(i + 1) % n];
        segs.push_back(Seg{a.first, a.second, b.first, b.second});
    }
    return segs;
}

// Closest point on segment s to (px,py), written to (cx,cy).
static inline void closestOnSeg(real_t px, real_t py, const Seg& s,
                                real_t& cx, real_t& cy)
{
    const real_t dx = s.x1 - s.x0, dy = s.y1 - s.y0;
    const real_t len2 = dx * dx + dy * dy;
    real_t t = (len2 > 0) ? ((px - s.x0) * dx + (py - s.y0) * dy) / len2 : (real_t)0;
    t = std::max<real_t>(0, std::min<real_t>(1, t));
    cx = s.x0 + t * dx;
    cy = s.y0 + t * dy;
}

// Robust ray-casting (crossing-number / even-odd) point-in-polygon test.
static inline bool pointInside(real_t px, real_t py, const std::vector<Seg>& segs)
{
    int cn = 0;
    for (size_t i = 0; i < segs.size(); ++i)
    {
        const Seg& s = segs[i];
        const bool crosses = (s.y0 <= py) != (s.y1 <= py);
        if (crosses)
        {
            const real_t xint = s.x0 + (py - s.y0) / (s.y1 - s.y0) * (s.x1 - s.x0);
            if (px < xint) ++cn;
        }
    }
    return (cn & 1) != 0;
}

// =============================================================================
//  SpainLevelSet: signed-distance level set of the closed coastline polygon.
//  Plugs straight into gsImplicitTrimmedDomain / gsAlgoimGenericRule
//  (negative inside, zero on the coastline, positive outside).
// =============================================================================
template<class T>
class SpainLevelSet : public gsFunction<T>
{
public:
    GISMO_CLONE_FUNCTION(SpainLevelSet)

    SpainLevelSet(std::vector<Seg> segs, const gsMatrix<T>& bbox)
        : m_segs(give(segs)), m_bbox(bbox) {}

    short_t     domainDim() const override { return 2; }
    short_t     targetDim() const override { return 1; }
    gsMatrix<T> support()   const override { return m_bbox; }

    // Signed distance to the coastline: -d inside, +d outside.
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k < u.cols(); ++k)
        {
            const T x = u(0, k), y = u(1, k);
            T best = std::numeric_limits<T>::max();
            for (size_t i = 0; i < m_segs.size(); ++i)
            {
                T cx, cy;
                closestOnSeg(x, y, m_segs[i], cx, cy);
                const T dx = x - cx, dy = y - cy;
                best = std::min(best, dx * dx + dy * dy);
            }
            const T dist = math::sqrt(best);
            result(0, k) = pointInside(x, y, m_segs) ? -dist : dist;
        }
    }

private:
    std::vector<Seg> m_segs;
    gsMatrix<T>      m_bbox;   // 2 x 2: col 0 = lower, col 1 = upper
};

// =============================================================================
//  Coastline sources
// =============================================================================

// (1) Plain text file with one ordered "x y" vertex per line.
static std::vector<std::pair<real_t,real_t> > loadPolygonTxt(const std::string& file)
{
    std::vector<std::pair<real_t,real_t> > pts;
    std::ifstream in(file.c_str());
    GISMO_ENSURE(in.good(), "Cannot open polygon text file: " << file);
    real_t x, y;
    while (in >> x >> y) pts.emplace_back(x, y);
    return pts;
}

// (2) Gmsh 2D mesh (MSH 4.1): extract the ordered coastline loop by chaining
//     the 1-D boundary line elements (entityDim == 1, elemType == 1).
//
//     We store every node coordinate, collect all boundary edges, then walk
//     the (degree-2) adjacency graph to recover a single ordered closed loop.
static std::vector<std::pair<real_t,real_t> > loadPolygonGmsh(const std::string& file)
{
    std::ifstream in(file.c_str());
    GISMO_ENSURE(in.good(), "Cannot open Gmsh mesh file: " << file);

    std::string line;
    // Version check (this reader targets the MSH 4.x block layout).
    real_t version = 0;
    while (std::getline(in, line))
        if (line.rfind("$MeshFormat", 0) == 0)
        {
            in >> version;
            std::getline(in, line); // finish the format line
            break;
        }
    GISMO_ENSURE(version >= 4.0,
                 "loadPolygonGmsh expects MSH 4.x; got version " << version
                 << ". Export as a plain 'x y' text file and use -f *.txt.");

    // ---- $Nodes : tag -> (x,y) ---------------------------------------------
    std::unordered_map<index_t, std::pair<real_t,real_t> > coord;
    while (std::getline(in, line))
        if (line.rfind("$Nodes", 0) == 0) break;

    index_t numBlocks = 0, numNodes = 0, minTag = 0, maxTag = 0;
    {
        std::getline(in, line);
        std::istringstream hs(line);
        hs >> numBlocks >> numNodes >> minTag >> maxTag;
    }
    coord.reserve(static_cast<size_t>(numNodes) * 2);
    for (index_t b = 0; b < numBlocks; ++b)
    {
        index_t edim, etag, parametric, nInBlock;
        std::getline(in, line);
        std::istringstream bh(line);
        bh >> edim >> etag >> parametric >> nInBlock;

        std::vector<index_t> tags(nInBlock);
        for (index_t i = 0; i < nInBlock; ++i)
        {
            std::getline(in, line);
            std::istringstream ts(line);
            ts >> tags[i];
        }
        for (index_t i = 0; i < nInBlock; ++i)
        {
            std::getline(in, line);
            std::istringstream cs(line);
            real_t x, y, z;
            cs >> x >> y >> z;
            coord[tags[i]] = std::make_pair(x, y);
        }
    }

    // ---- $Elements : collect boundary edges (dim==1, type==1) --------------
    while (std::getline(in, line))
        if (line.rfind("$Elements", 0) == 0) break;

    index_t eNumBlocks = 0, eNumElems = 0, eMin = 0, eMax = 0;
    {
        std::getline(in, line);
        std::istringstream hs(line);
        hs >> eNumBlocks >> eNumElems >> eMin >> eMax;
    }

    std::vector<std::pair<index_t,index_t> > edges;
    for (index_t b = 0; b < eNumBlocks; ++b)
    {
        index_t edim, etag, etype, nInBlock;
        std::getline(in, line);
        std::istringstream bh(line);
        bh >> edim >> etag >> etype >> nInBlock;

        const bool boundaryLines = (edim == 1 && etype == 1); // 2-node line
        for (index_t i = 0; i < nInBlock; ++i)
        {
            std::getline(in, line);
            if (!boundaryLines) continue;
            std::istringstream es(line);
            index_t etagId, n1, n2;
            es >> etagId >> n1 >> n2;
            edges.emplace_back(n1, n2);
        }
    }
    GISMO_ENSURE(!edges.empty(),
                 "No boundary line elements (dim=1) found in " << file);

    // ---- Chain edges into a single ordered closed loop ---------------------
    std::unordered_map<index_t, std::vector<index_t> > adj;
    adj.reserve(edges.size() * 2);
    for (size_t i = 0; i < edges.size(); ++i)
    {
        adj[edges[i].first ].push_back(edges[i].second);
        adj[edges[i].second].push_back(edges[i].first );
    }

    const index_t start = edges.front().first;
    std::vector<std::pair<real_t,real_t> > pts;
    pts.reserve(edges.size() + 1);

    index_t prev = -1, cur = start;
    for (size_t guard = 0; guard <= edges.size() + 1; ++guard)
    {
        std::unordered_map<index_t, std::pair<real_t,real_t> >::const_iterator
            it = coord.find(cur);
        GISMO_ENSURE(it != coord.end(),
                     "Boundary node " << cur << " missing from $Nodes.");
        pts.push_back(it->second);

        const std::vector<index_t>& nb = adj[cur];
        index_t next = -1;
        for (size_t j = 0; j < nb.size(); ++j)
            if (nb[j] != prev) { next = nb[j]; break; }

        if (next == -1 || next == start) break; // closed the loop
        prev = cur;
        cur  = next;
    }
    return pts;
}

// =============================================================================
//  FCM helpers (identical in spirit to area_fcm_example.cpp)
// =============================================================================

// Append 2D points (2xN) + scalar weights as a 4xN cloud (x, y, z=0, scalar).
static void appendCloud2D(gsMatrix<real_t>& cloud,
                          const gsMatrix<real_t>& pts2d,
                          const gsVector<real_t>& scalar)
{
    if (pts2d.cols() == 0) return;
    const index_t c = cloud.cols();
    cloud.conservativeResize(4, c + pts2d.cols());
    cloud.block(0, c, 2, pts2d.cols()) = pts2d;
    cloud.row(2).segment(c, pts2d.cols()).setZero();
    cloud.row(3).segment(c, pts2d.cols()) = scalar.transpose();
}

// Tensor Gauss rule on a fully-inside box.  These points are used only where
// no Algoim cut-cell rule is needed; the integral is still exactly the box area.
static void appendInteriorGauss2D(index_t nGauss,
                                  const gsVector<real_t>& lo,
                                  const gsVector<real_t>& hi,
                                  gsMatrix<real_t>& pts,
                                  gsVector<real_t>& wts)
{
    gsGaussRule<real_t> rule(gsVector<index_t,2>::Constant(nGauss));
    gsMatrix<real_t> qPts;
    gsVector<real_t> qWts;
    rule.mapTo(lo, hi, qPts, qWts);

    const index_t c = pts.cols();
    pts.conservativeResize(2, c + qPts.cols());
    pts.block(0, c, 2, qPts.cols()) = qPts;
    wts.conservativeResize(c + qWts.size());
    wts.segment(c, qWts.size()) = qWts;
}

// Lipschitz radius cell classification (signed distance is 1-Lipschitz):
//   phi(center) >  radius  => +1 (outside)
//   phi(center) < -radius  => -1 (inside -> exact area)
//   otherwise              =>  0 (cut -> recurse / Algoim)
static int boxClassPhi2D(const gsFunction<real_t>& phi,
                         const gsVector<real_t>& lo,
                         const gsVector<real_t>& hi)
{
    gsMatrix<real_t> c(2, 1), v;
    c.col(0) = 0.5*(lo + hi);
    phi.eval_into(c, v);
    const real_t radius = (real_t)0.5 * (hi - lo).norm();
    if (v(0,0) >  radius) return  1;
    if (v(0,0) < -radius) return -1;
    return 0;
}

// Algoim {phi<0} area rule on a leaf box, with a uniform midpoint fallback
// when Algoim returns zero nodes (e.g. re-entrant / multi-crossing cells).
static void leafCutCellPoints2D(const gsAlgoimGenericRule<real_t>& volRule,
                                const gsFunction<real_t>& phi,
                                index_t nFallback,
                                const gsVector<real_t>& lo, const gsVector<real_t>& hi,
                                gsMatrix<real_t>& pts, gsVector<real_t>& wts)
{
    volRule.mapTo(lo, hi, pts, wts);
    if (wts.size() > 0) return;

    const index_t nTotal = nFallback * nFallback;
    gsMatrix<real_t> spts(2, nTotal);
    index_t col = 0;
    for (index_t ii = 0; ii < nFallback; ++ii)
    for (index_t jj = 0; jj < nFallback; ++jj, ++col)
    {
        spts(0,col) = lo[0] + (ii + (real_t)0.5) / nFallback * (hi[0] - lo[0]);
        spts(1,col) = lo[1] + (jj + (real_t)0.5) / nFallback * (hi[1] - lo[1]);
    }

    gsMatrix<real_t> vals;
    phi.eval_into(spts, vals);
    std::vector<index_t> inside;
    for (index_t ci = 0; ci < nTotal; ++ci)
        if (vals(0,ci) < 0) inside.push_back(ci);

    if (inside.empty()) { pts.resize(2,0); wts.resize(0); return; }

    const real_t w = (hi - lo).prod() / static_cast<real_t>(nTotal);
    pts.resize(2, static_cast<index_t>(inside.size()));
    wts.resize(   static_cast<index_t>(inside.size()));
    for (std::size_t ci = 0; ci < inside.size(); ++ci)
    {
        pts.col(static_cast<index_t>(ci)) = spts.col(inside[ci]);
        wts[static_cast<index_t>(ci)] = w;
    }
}

// Adaptive quad-tree over one background cut cell.
static void collectAdaptive2D(const gsAlgoimGenericRule<real_t>& volRule,
                              const gsFunction<real_t>& phi,
                              index_t quadDepth, index_t nFallback,
                              const gsVector<real_t>& lo, const gsVector<real_t>& hi,
                              index_t depth,
                              index_t& nSubBoxes,
                              gsMatrix<real_t>& subIntCtr, gsVector<real_t>& subIntArea,
                              gsMatrix<real_t>& leafPts,   gsVector<real_t>& leafWts)
{
    if (depth >= quadDepth)
    {
        gsMatrix<real_t> pts; gsVector<real_t> wts;
        leafCutCellPoints2D(volRule, phi, nFallback, lo, hi, pts, wts);
        if (pts.cols() > 0)
        {
            const index_t c = leafPts.cols();
            leafPts.conservativeResize(2, c + pts.cols());
            leafPts.block(0, c, 2, pts.cols()) = pts;
            leafWts.conservativeResize(c + wts.size());
            leafWts.segment(c, wts.size()) = wts;
        }
        return;
    }

    const gsVector<real_t> mid = 0.5*(lo + hi);
    gsVector<real_t> clo(2), chi(2);
    for (int ix = 0; ix < 2; ++ix)
    for (int iy = 0; iy < 2; ++iy)
    {
        clo[0] = ix ? mid[0] : lo[0];  chi[0] = ix ? hi[0] : mid[0];
        clo[1] = iy ? mid[1] : lo[1];  chi[1] = iy ? hi[1] : mid[1];

        const int cls = boxClassPhi2D(phi, clo, chi);
        if (cls > 0) continue;   // outside: skip
        if (cls < 0)             // inside: exact area
        {
            ++nSubBoxes;
            appendInteriorGauss2D(nFallback, clo, chi, subIntCtr, subIntArea);
            continue;
        }
        collectAdaptive2D(volRule, phi, quadDepth, nFallback, clo, chi, depth + 1,
                          nSubBoxes,
                          subIntCtr, subIntArea, leafPts, leafWts);
    }
}

// =============================================================================
//  main
// =============================================================================
int main(int argc, char* argv[])
{
    std::string file      =
        "/Users/lucasventavinuela/gismo_gmsh/optional/gsGmsh/filedata/spain_mesh.msh";
    std::string outDir    = "output_spain_fcm";
    index_t     numRef    = 5;
    index_t     degree    = 2;
    index_t     quadDepth = 1;
    std::string indicator = "integralChange";
    real_t      indicatorTol = 1e-2;
    real_t      fill      = (real_t)0.9;
    bool        plot      = false;

    gsCmdLine cmd("2D FCM: area and perimeter of the Spain coastline "
                  "(immersed boundary) via adaptive quad-tree + Algoim.");
    cmd.addString("f", "file",      "Coastline source: Gmsh .msh (MSH 4.x) or "
                                    "text file of ordered 'x y' vertices",     file);
    cmd.addString("o", "output",    "Output folder for results/plots",          outDir);
    cmd.addInt   ("r", "refine",    "Number of uniform refinements",            numRef);
    cmd.addInt   ("e", "degree",    "B-spline degree of the background mesh",    degree);
    cmd.addInt   ("",  "quadDepth", "Adaptive quad-tree depth for cut cells",   quadDepth);
    cmd.addString("",  "indicator", "Adaptive indicator: uniform, fallback or "
                                    "integralChange",                          indicator);
    cmd.addReal  ("",  "indicatorTol", "integralChange acceptance tolerance",  indicatorTol);
    cmd.addReal  ("",  "fill",      "Fill fraction of [0,1]^2 (0..1)",          fill);
    cmd.addSwitch("plot", "Write ParaView output",                              plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // -------------------------------------------------------------------------
    // 1. Read the coastline as an ordered closed polygon.
    // -------------------------------------------------------------------------
    const bool isMsh = file.size() >= 4 &&
                       file.compare(file.size() - 4, 4, ".msh") == 0;

    std::vector<std::pair<real_t,real_t> > pts =
        isMsh ? loadPolygonGmsh(file) : loadPolygonTxt(file);

    GISMO_ENSURE(pts.size() >= 3, "Coastline needs at least 3 vertices (got "
                                  << pts.size() << ").");

    // -------------------------------------------------------------------------
    // 2. Compute physical bounding box, then rescale into parametric [0,1]^2.
    //    Mirrors the normalisation used in volume_fcm_octree_adaptive_example.
    // -------------------------------------------------------------------------
    real_t xmin =  std::numeric_limits<real_t>::max(), ymin = xmin;
    real_t xmax = -std::numeric_limits<real_t>::max(), ymax = xmax;
    for (size_t i = 0; i < pts.size(); ++i)
    {
        xmin = std::min(xmin, pts[i].first);  xmax = std::max(xmax, pts[i].first);
        ymin = std::min(ymin, pts[i].second); ymax = std::max(ymax, pts[i].second);
    }
    const real_t cx     = (real_t)0.5 * (xmin + xmax);
    const real_t cy     = (real_t)0.5 * (ymin + ymax);
    const real_t extent = std::max(xmax - xmin, ymax - ymin);
    const real_t scale  = fill / extent;

    // Rescale every vertex into [0,1]^2 (longest axis maps to [0.5-fill/2, 0.5+fill/2]).
    for (size_t i = 0; i < pts.size(); ++i)
    {
        pts[i].first  = (pts[i].first  - cx) * scale + (real_t)0.5;
        pts[i].second = (pts[i].second - cy) * scale + (real_t)0.5;
    }

    // Build polygon segments and reference values from the rescaled coordinates.
    std::vector<Seg> segs = polygonSegments(pts);

    real_t exactArea = 0, exactPerim = 0;
    for (size_t i = 0; i < segs.size(); ++i)
    {
        const Seg& s = segs[i];
        exactArea  += s.x0 * s.y1 - s.x1 * s.y0;
        exactPerim += math::sqrt((s.x1 - s.x0)*(s.x1 - s.x0)
                                + (s.y1 - s.y0)*(s.y1 - s.y0));
    }
    exactArea = (real_t)0.5 * std::abs(exactArea);

    // Back-conversion factors: param -> physical units.
    const real_t invScale2 = (real_t)1.0 / (scale * scale);
    const real_t invScale  = (real_t)1.0 / scale;

    gsInfo << "Coastline source   : " << file << "\n"
           << "Vertices           : " << pts.size() << "\n"
           << "Edges              : " << segs.size() << "\n"
           << "Physical bbox      : [" << xmin << ", " << xmax << "] x ["
                                      << ymin << ", " << ymax << "]\n"
           << "fill / scale       : " << fill << " / " << scale << "\n"
           << "Ref area  (param)  : " << exactArea           << "\n"
           << "Ref area  (phys)   : " << exactArea * invScale2 << " (shoelace)\n"
           << "Ref perim (param)  : " << exactPerim           << "\n"
           << "Ref perim (phys)   : " << exactPerim * invScale << " (polygon edge length)\n"
           << "degree             : " << degree << "\n"
           << "quadDepth          : " << quadDepth << "\n\n";

    // -------------------------------------------------------------------------
    // 3. Background mesh on [0,1]^2; level set works in parametric coordinates.
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(2, 2);
    bbox.col(0) << (real_t)0, (real_t)0;
    bbox.col(1) << (real_t)1, (real_t)1;

    SpainLevelSet<real_t> phi(segs, bbox);

    gsInfo << std::setw(6)  << "ref"
           << std::setw(12) << "cells/dim"
           << std::setw(16) << "area(par)"
           << std::setw(16) << "area(phys)"
           << std::setw(14) << "perim(par)"
           << std::setw(14) << "perim(phys)"
           << std::setw(12) << "intCells"
           << std::setw(12) << "quadPts"
           << std::setw(12) << "bdryPts"
           << std::setw(12) << "subBoxes"
            << std::setw(10) << "fallback"
           << std::setw(13) << "err(area)"
           << std::setw(13) << "err(perim)\n";

    gsFileManager::mkdir(outDir);
    gsFileManager::mkdir(outDir + "/results");

    std::ofstream fout((outDir + "/results/area.txt").c_str());
    fout << "# 2D area and perimeter of the Spain coastline via FCM "
         << "(immersed boundary) + adaptive quad-tree + Algoim\n"
         << "# source: " << file << "\n"
         << "# quadDepth = " << quadDepth << "\n"
         << "# NOTE: (par) = parametric [0,1]^2 coords; "
            "(phys) = original physical units.\n\n"
         << std::setw(6)  << "#ref" << std::setw(12) << "cells/dim"
         << std::setw(16) << "area(par)" << std::setw(16) << "area(phys)"
         << std::setw(14) << "perim(par)" << std::setw(14) << "perim(phys)"
         << std::setw(12) << "intCells"
         << std::setw(12) << "quadPts"
         << std::setw(12) << "bdryPts"
         << std::setw(12) << "subBoxes"
         << std::setw(10) << "fallback"
         << std::setw(16) << "err_area" << std::setw(16) << "err_perim\n";

    // ParaView collections spanning all refinement levels.
    std::unique_ptr<gsParaviewCollection> colInterior, colCut, colAll, colBg;
    if (plot)
    {
        gsFileManager::mkdir(outDir + "/points_interior");
        gsFileManager::mkdir(outDir + "/points_cutcells");
        gsFileManager::mkdir(outDir + "/points_all");
        gsFileManager::mkdir(outDir + "/background");
        colInterior.reset(new gsParaviewCollection(outDir + "/points_interior/interior"));
        colCut     .reset(new gsParaviewCollection(outDir + "/points_cutcells/cutcells"));
        colAll     .reset(new gsParaviewCollection(outDir + "/points_all/all"));
        colBg      .reset(new gsParaviewCollection(outDir + "/background/background"));

        // Level set field (phi=0 isoline is the coastline).
        gsFileManager::mkdir(outDir + "/levelset");
        gsWriteParaview(phi, bbox, outDir + "/levelset/spain_levelset", 40000);
    }

    // -------------------------------------------------------------------------
    // 4. Refinement loop (Finite Cell Method).
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRef; ++r)
    {
        const int nElem = 1 << r;  // 2^r elements per direction over [0,1]^2
        gsKnotVector<> kvx((real_t)0, (real_t)1, nElem - 1,
                   static_cast<short_t>(degree + 1), 1,
                   static_cast<short_t>(degree));
        gsKnotVector<> kvy((real_t)0, (real_t)1, nElem - 1,
                   static_cast<short_t>(degree + 1), 1,
                   static_cast<short_t>(degree));
        gsTensorBSplineBasis<2,real_t> bkgBasis(kvx, kvy);
        const index_t nFallback = bkgBasis.maxDegree() + 1;

        gsImplicitTrimmedDomain<2,real_t> tr_domain(phi, bkgBasis);

        gsOptionList adaptiveOptions = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
        adaptiveOptions.setInt("maxDepth", quadDepth);
        adaptiveOptions.setInt("nFallback", nFallback);
        adaptiveOptions.setString("indicator", indicator);
        adaptiveOptions.setReal("indicatorTol", indicatorTol);
        adaptiveOptions.setSwitch("analyticInterior", true);
        gsAlgoimAdaptiveRule<real_t> volumeRule(phi, bkgBasis, adaptiveOptions);

        gsOptionList surfaceOptions = adaptiveOptions;
        surfaceOptions.setInt("dim", 2);
        gsAlgoimAdaptiveRule<real_t> surfaceRule(phi, bkgBasis, surfaceOptions);

        // -- Background interior cells: exact area (no quadrature needed)
        real_t           areaInterior = 0;
        index_t          nInterior    = 0;
        index_t          nInteriorQPs = 0;
        gsMatrix<real_t> intCtr(2, 0);
        gsVector<real_t> intArea;

        for (auto it = tr_domain.beginInterior();
             it != tr_domain.end<InteriorSign>(); ++it)
        {
            const real_t cellArea = (it.upperCorner() - it.lowerCorner()).prod();
            areaInterior += cellArea;
            ++nInterior;
            nInteriorQPs += nFallback * nFallback;
            if (plot)
                appendInteriorGauss2D(nFallback, it.lowerCorner(), it.upperCorner(),
                                      intCtr, intArea);
        }

        // -- Background cut cells: adaptive quad-tree + Algoim at leaves
        real_t          areaCut   = 0;
        real_t          perim     = 0;
        index_t         nSubInteriorQPs = 0;
        index_t         nBoundaryQPs    = 0;
        index_t         nSubBoxes       = 0;

        gsMatrix<real_t> quadIn (4, 0);
        gsMatrix<real_t> quadCut(4, 0);
        gsMatrix<real_t> quadAll(4, 0);

        if (plot) appendCloud2D(quadIn, intCtr, intArea);

        for (auto it = tr_domain.beginBdr(boundary::none);
             it != tr_domain.endBdr(boundary::none); ++it)
        {
            gsMatrix<real_t> subIntCtr(2, 0), leafPts(2, 0);
            gsVector<real_t> subIntArea, leafWts;

            const auto classifier = [&](const gsVector<real_t> & childLo,
                                        const gsVector<real_t> & childHi) -> int
            {
                return boxClassPhi2D(phi, childLo, childHi);
            };

            volumeRule.mapToSeparated(it.lowerCorner(), it.upperCorner(),
                                      subIntCtr, subIntArea,
                                      leafPts, leafWts, classifier);

            areaCut   += subIntArea.sum() + leafWts.sum();
            nSubInteriorQPs += subIntCtr.cols();
            nBoundaryQPs    += leafPts.cols();

            if (plot)
            {
                appendCloud2D(quadIn,  subIntCtr, subIntArea);
                appendCloud2D(quadCut, leafPts,   leafWts);
                appendCloud2D(quadAll, subIntCtr, subIntArea);
                appendCloud2D(quadAll, leafPts,   leafWts);
            }

            gsMatrix<real_t> surfInterior(2, 0), surfCut(2, 0);
            gsVector<real_t> surfInteriorW, surfCutW;
            surfaceRule.mapToSeparated(it.lowerCorner(), it.upperCorner(),
                                       surfInterior, surfInteriorW,
                                       surfCut, surfCutW, classifier);
            perim += surfCutW.sum();
        }

        if (plot) appendCloud2D(quadAll, intCtr, intArea);

        nSubBoxes = volumeRule.stats().nSubBoxes;   // accumulated over the loop
        const real_t area = areaInterior + areaCut;
    const index_t nQuadPts = nSubInteriorQPs + nBoundaryQPs;

        // -- Console output
        gsInfo << std::setw(6)  << r
               << std::setw(12) << nElem
               << std::setw(16) << std::fixed << std::setprecision(6) << area
               << std::setw(16) << std::fixed << std::setprecision(6) << area  * invScale2
               << std::setw(14) << std::fixed << std::setprecision(6) << perim
               << std::setw(14) << std::fixed << std::setprecision(6) << perim * invScale
               << std::setw(12) << nInterior
               << std::setw(12) << nQuadPts
               << std::setw(12) << nBoundaryQPs
               << std::setw(12) << nSubBoxes
               << std::setw(10) << nFallback
               << std::setw(13) << std::scientific << std::setprecision(2)
               << std::abs(area  - exactArea)
               << std::setw(13) << std::abs(perim - exactPerim) << "\n";

        // -- Results file
        fout   << std::setw(6)  << r
               << std::setw(12) << nElem
               << std::setw(16) << std::fixed << std::setprecision(8) << area
               << std::setw(16) << std::fixed << std::setprecision(8) << area  * invScale2
               << std::setw(14) << std::fixed << std::setprecision(8) << perim
               << std::setw(14) << std::fixed << std::setprecision(8) << perim * invScale
               << std::setw(12) << nInterior
               << std::setw(12) << nQuadPts
               << std::setw(12) << nBoundaryQPs
               << std::setw(12) << nSubBoxes
               << std::setw(10) << nFallback
               << std::setw(16) << std::scientific << std::setprecision(3)
               << std::abs(area  - exactArea)
               << std::setw(16) << std::abs(perim - exactPerim) << "\n";

        // -- ParaView output
        if (plot)
        {
            const std::string rs = std::to_string(r);

            if (quadIn.cols() > 0)
            {
                gsWriteParaviewPoints(quadIn,
                    outDir + "/points_interior/interior_r" + rs);
                colInterior->addPart("interior_r" + rs + ".vtp", r);
            }
            if (quadCut.cols() > 0)
            {
                gsWriteParaviewPoints(quadCut,
                    outDir + "/points_cutcells/cutcells_r" + rs);
                colCut->addPart("cutcells_r" + rs + ".vtp", r);
            }
            if (quadAll.cols() > 0)
            {
                gsWriteParaviewPoints(quadAll,
                    outDir + "/points_all/all_r" + rs);
                colAll->addPart("all_r" + rs + ".vtp", r);
            }

            gsMesh<real_t> bgMesh(bkgBasis);
            gsWriteParaview(bgMesh, outDir + "/background/background_r" + rs);
            colBg->addPart("background_r" + rs + ".vtp", r);
        }
    }

    fout.close();
    gsInfo << "\nResults written to: " << outDir << "/results/area.txt\n";

    if (plot)
    {
        colInterior->save(); colCut->save(); colAll->save(); colBg->save();
        gsInfo << "ParaView files written to: " << outDir << "\n";
    }

    return 0;
}
