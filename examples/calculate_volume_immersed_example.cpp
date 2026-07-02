/** @file calculate_volume_immersed_example.cpp

    @brief Volume of a 3D triangle mesh (cow.obj) via the Finite Cell Method.

    Pipeline (3D analogue of poisson2_nitsche_immersed_example.cpp):

      1. Read a triangle mesh from a Wavefront .obj file (the cow).
      2. Rescale/center the mesh into the unit parametric box [0,1]^3.
      3. Build an analytic level set phi(x,y,z) from the triangle soup:
           - distance : minimal point-to-triangle distance to the surface.
           - sign     : generalized (solid-angle) winding number; phi < 0
                        inside, phi > 0 outside (the convention used by
                        gsImplicitTrimmedDomain).
      4. Classify the background cells of [0,1]^3 with
         gsImplicitTrimmedDomain<3> (interior / cut / exterior).
      5. Integrate the constant 1 over the immersed domain using the Algoim
         quadrature rule on cut cells (phi < 0 volumetric rule) and a full
         Gauss rule on interior cells.  The result is the volume of the cow
         in the parametric box; dividing by the rescaling factor recovers the
         physical volume.
      6. (--plot) Export the rescaled cow surface together with the Algoim
         volumetric quadrature points to ParaView.

    Usage:
    ./bin/calculate_volume_immersed_example [-f cow.obj] [-o output_cow3d] [-r 3] [-e 0] [--fill 0.9] [--plot]
    ./build/bin/calculate_volume_immersed_example -f optional/gsGmsh/filedata/spot.obj -o output_spot -r 3 --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
//  Wavefront .obj loader (vertices + triangulated faces)
// =============================================================================
//  Reads "v x y z" vertex lines and "f ..." face lines.  Faces may use the
//  v, v/vt, v//vn or v/vt/vn notations; only the vertex index is used.  Faces
//  with more than three corners are triangulated with a simple fan.
static bool loadObjMesh(const std::string & filename,
                        std::vector<double> & verts,   // 3 * Nv, row-major (x,y,z)
                        std::vector<int>    & tris)     // 3 * Nt, zero-based
{
    std::ifstream in(filename.c_str());
    if (!in.is_open())
        return false;

    verts.clear();
    tris.clear();

    std::string line;
    while (std::getline(in, line))
    {
        if (line.size() < 2) continue;
        if (line[0] == 'v' && line[1] == ' ')
        {
            std::istringstream ss(line.substr(2));
            double x, y, z;
            ss >> x >> y >> z;
            verts.push_back(x);
            verts.push_back(y);
            verts.push_back(z);
        }
        else if (line[0] == 'f' && line[1] == ' ')
        {
            std::istringstream ss(line.substr(2));
            std::vector<int> face;
            std::string tok;
            while (ss >> tok)
            {
                // Keep only the part before the first '/'.
                const std::size_t slash = tok.find('/');
                if (slash != std::string::npos)
                    tok = tok.substr(0, slash);
                if (tok.empty()) continue;
                int idx = std::atoi(tok.c_str());
                if (idx < 0) // negative = relative index
                    idx = static_cast<int>(verts.size() / 3) + idx;
                else
                    idx -= 1; // .obj is 1-based
                face.push_back(idx);
            }
            // Fan triangulation of the (convex) polygon.
            for (std::size_t k = 2; k < face.size(); ++k)
            {
                tris.push_back(face[0]);
                tris.push_back(face[k - 1]);
                tris.push_back(face[k]);
            }
        }
    }
    return !verts.empty() && !tris.empty();
}

static std::string fileStem(const std::string & filename)
{
    const std::size_t slash = filename.find_last_of("/\\");
    const std::size_t begin = (slash == std::string::npos) ? 0 : slash + 1;
    const std::size_t dot = filename.find_last_of('.');
    const std::size_t end = (dot == std::string::npos || dot < begin) ? filename.size() : dot;
    return filename.substr(begin, end - begin);
}

// =============================================================================
//  Geometry helpers (plain doubles)
// =============================================================================
typedef std::array<double, 3> Vec3;

static inline Vec3 sub(const Vec3 & a, const Vec3 & b)
{ return {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}}; }
static inline double dot(const Vec3 & a, const Vec3 & b)
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static inline Vec3 cross(const Vec3 & a, const Vec3 & b)
{ return {{a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}}; }
static inline double norm(const Vec3 & a) { return std::sqrt(dot(a, a)); }

// Squared distance from point p to triangle (a,b,c) (Ericson, RTCD).
static double sqDistPointTriangle(const Vec3 & p, const Vec3 & a,
                                  const Vec3 & b, const Vec3 & c)
{
    const Vec3 ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) { const Vec3 d = sub(p, a); return dot(d, d); }

    const Vec3 bp = sub(p, b);
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) { const Vec3 d = sub(p, b); return dot(d, d); }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        const Vec3 q = {{a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]}};
        const Vec3 d = sub(p, q); return dot(d, d);
    }

    const Vec3 cp = sub(p, c);
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) { const Vec3 d = sub(p, c); return dot(d, d); }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        const Vec3 q = {{a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]}};
        const Vec3 d = sub(p, q); return dot(d, d);
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const Vec3 bc = sub(c, b);
        const Vec3 q = {{b[0] + w * bc[0], b[1] + w * bc[1], b[2] + w * bc[2]}};
        const Vec3 d = sub(p, q); return dot(d, d);
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    const Vec3 q = {{a[0] + ab[0] * v + ac[0] * w,
                     a[1] + ab[1] * v + ac[1] * w,
                     a[2] + ab[2] * v + ac[2] * w}};
    const Vec3 d = sub(p, q); return dot(d, d);
}

// Signed solid angle subtended by triangle (a,b,c) at point p
// (Van Oosterom & Strackee).  The sum over a closed mesh equals
// 4*pi * (generalized winding number).
static double signedSolidAngle(const Vec3 & p, const Vec3 & a,
                               const Vec3 & b, const Vec3 & c)
{
    const Vec3 va = sub(a, p), vb = sub(b, p), vc = sub(c, p);
    const double la = norm(va), lb = norm(vb), lc = norm(vc);
    const double numer = dot(va, cross(vb, vc));
    const double denom = la * lb * lc
                       + dot(va, vb) * lc
                       + dot(vb, vc) * la
                       + dot(vc, va) * lb;
    return 2.0 * std::atan2(numer, denom);
}

// =============================================================================
//  gsMeshSignedDist<T>
//  Signed-distance level set of a closed triangle mesh.
//
//  Convention (matches gsImplicitTrimmedDomain):
//    phi < 0  ->  inside    (active / interior)
//    phi > 0  ->  outside
//    phi = 0  ->  on the surface (cut cell)
// =============================================================================
template<class T>
class gsMeshSignedDist : public gsFunction<T>
{
public:
    GISMO_CLONE_FUNCTION(gsMeshSignedDist)

    gsMeshSignedDist(std::vector<double> verts, std::vector<int> tris,
                     const gsMatrix<T> & bbox)
        : m_verts(give(verts)), m_tris(give(tris)), m_bbox(bbox)
    {}

    short_t     domainDim() const override { return 3; }
    short_t     targetDim() const override { return 1; }
    gsMatrix<T> support()   const override { return m_bbox; }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        const std::size_t nt = m_tris.size() / 3;
        result.resize(1, u.cols());
        for (index_t k = 0; k < u.cols(); ++k)
        {
            const Vec3 p = {{ static_cast<double>(u(0, k)),
                              static_cast<double>(u(1, k)),
                              static_cast<double>(u(2, k)) }};

            double best  = std::numeric_limits<double>::max();
            double omega = 0.0;
            for (std::size_t t = 0; t < nt; ++t)
            {
                const int ia = m_tris[3 * t + 0];
                const int ib = m_tris[3 * t + 1];
                const int ic = m_tris[3 * t + 2];
                const Vec3 a = {{ m_verts[3 * ia], m_verts[3 * ia + 1], m_verts[3 * ia + 2] }};
                const Vec3 b = {{ m_verts[3 * ib], m_verts[3 * ib + 1], m_verts[3 * ib + 2] }};
                const Vec3 c = {{ m_verts[3 * ic], m_verts[3 * ic + 1], m_verts[3 * ic + 2] }};

                const double d2 = sqDistPointTriangle(p, a, b, c);
                if (d2 < best) best = d2;
                omega += signedSolidAngle(p, a, b, c);
            }

            const double pi      = 3.14159265358979323846;
            const double dist    = std::sqrt(best);
            const double winding = omega / (4.0 * pi);
            // |winding| ~ 1 inside, ~ 0 outside (sign-robust to mesh orientation)
            const bool inside = std::abs(winding) > 0.5;
            result(0, k) = inside ? -static_cast<T>(dist) : static_cast<T>(dist);
        }
    }

private:
    std::vector<double> m_verts; // 3 * Nv
    std::vector<int>    m_tris;  // 3 * Nt (zero-based)
    gsMatrix<T>         m_bbox;  // 3 x 2 : col 0 = lower, col 1 = upper
};

// =============================================================================
//  main
// =============================================================================
int main(int argc, char * argv[])
{
    std::string filename =
        "/Users/lucasventavinuela/gismo_gmsh/optional/gsGmsh/filedata/cow.obj";
    std::string out = gsFileManager::getCanonicRepresentation(
        gsFileManager::getPath(std::string(__FILE__), true) + "../output_cow3d");
    index_t numRefine  = 3;
    index_t numElevate = 0;
    real_t  fill       = 0.9;  // fraction of [0,1]^3 spanned by the largest extent
    bool    plot       = false;

    gsCmdLine cmd("Volume of a 3D triangle mesh from an .obj file via the Finite Cell Method.");
    cmd.addString("f", "file",            "Wavefront .obj mesh file",            filename);
    cmd.addString("o", "output",          "Output folder for ParaView files",    out);
    cmd.addInt   ("r", "uniformRefine",   "Number of uniform refinement steps",  numRefine);
    cmd.addInt   ("e", "degreeElevation", "Number of degree elevation steps",    numElevate);
    cmd.addReal  ("",  "fill",            "Fill fraction of [0,1]^3 (0..1)",      fill);
    cmd.addSwitch("plot", "Create ParaView output",                              plot);
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

    gsInfo << "Physical bounding box: ["
           << lo[0] << ", " << lo[1] << ", " << lo[2] << "] - ["
           << hi[0] << ", " << hi[1] << ", " << hi[2] << "]\n";
    gsInfo << "Rescaling factor (param/phys): " << scale << "\n";

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

    // -------------------------------------------------------------------------
    // 5. Output directory
    // -------------------------------------------------------------------------
    gsFileManager::mkdir(out);

    gsInfo << "\nDegree of basis: " << dbasis.maxCwiseDegree() << "\n";
    gsInfo << "Computing the volume of the immersed domain (cow) in [0,1]^3.\n\n";
    gsInfo << std::setw(8) << "refine" << std::setw(10) << "cells/dim"
           << std::setw(18) << "vol (param)" << std::setw(18) << "vol (phys)" << "\n";

    const double invScale3 = 1.0 / (scale * scale * scale);

    // Number of uniform sample points per direction used as fallback when Algoim
    // returns 0 nodes (re-entrant / multi-crossing cut cells).  n^3 = 64 pts.
    const index_t nFallback = 4;

    // -------------------------------------------------------------------------
    // 6. Refinement loop: classify cells and integrate 1 over the cow
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        gsTensorBSplineBasis<3, real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<3, real_t> *>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        // Cut-cell classification of the background grid.
        gsImplicitTrimmedDomain<3, real_t> tr_domain(impl_fun, *tbsPtr);

        // Algoim volumetric rule for {phi < 0} on cut cells (default dim = -1).
        gsAlgoimGenericRule<real_t> algoimRule(impl_fun, *tbsPtr);

        real_t volume = 0;
        gsMatrix<real_t> pts;
        gsVector<real_t> wts;

        // Interior cells (sign < 0) contribute their full (unit-Jacobian) volume.
        for (auto it = tr_domain.beginInterior();
             it != tr_domain.end<InteriorSign>(); ++it)
            volume += (it.upperCorner() - it.lowerCorner()).prod();

        // Cut cells (sign == 0): Algoim volumetric quadrature of {phi < 0}.
        // When Algoim returns no nodes (re-entrant / multi-crossing cell),
        // fall back to a uniform grid fraction estimate.
        for (auto it = tr_domain.beginBdr(boundary::none);
             it != tr_domain.endBdr(boundary::none); ++it)
        {
            algoimRule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            if (wts.size() > 0)
            {
                volume += wts.sum();
            }
            else
            {
                const gsVector<real_t> & lo_c = it.lowerCorner();
                const gsVector<real_t> & hi_c = it.upperCorner();
                const real_t cellVol = (hi_c - lo_c).prod();
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
                const index_t nInside = static_cast<index_t>((fbVals.array() < 0).count());
                volume += cellVol * static_cast<real_t>(nInside) / nTotal;
            }
        }

        const index_t cellsPerDim = tbsPtr->component(0).numElements();
        gsInfo << std::setw(8) << r << std::setw(10) << cellsPerDim
               << std::setw(18) << std::fixed << std::setprecision(6) << volume
               << std::setw(18) << std::fixed << std::setprecision(6)
               << volume * invScale3 << "\n";

        if (r == numRefine)
        {
            const real_t bboxVol = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
            gsInfo << "\n=== Result (finest level) ===\n";
            gsInfo << "Volume of object (parametric [0,1]^3): "
                   << std::fixed << std::setprecision(6) << volume << "\n";
            gsInfo << "Volume of object (physical units):    "
                   << std::fixed << std::setprecision(6) << volume * invScale3 << "\n";
            gsInfo << "Percentage of [0,1]^3 box occupied:   "
                   << std::fixed << std::setprecision(2) << 100.0 * volume << " %\n";
            gsInfo << "Percentage of object bounding box:    "
                   << std::fixed << std::setprecision(2)
                   << 100.0 * (volume * invScale3) / bboxVol << " %\n\n";
        }

        // ---------------------------------------------------------------------
        // 7. ParaView output (final refinement level only)
        // ---------------------------------------------------------------------
        if (plot && r == numRefine)
        {
            // 7a. Rescaled input geometry as a triangle mesh.
            gsMesh<real_t> geomMesh;
            for (std::size_t i = 0; i < nVert; ++i)
                geomMesh.addVertex(verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]);
            for (std::size_t t = 0; t < nTri; ++t)
                geomMesh.addFace(tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]);
            gsWriteParaview(geomMesh, out + "/" + outputStem + "_geometry");

            // 7b. Background grid.
            gsMesh<real_t> bgMesh(dbasis.basis(0));
            gsWriteParaview(bgMesh, out + "/" + outputStem + "_background_mesh");

            // 7c. Algoim volumetric quadrature points covering the WHOLE cow.
            //     Cut cells (sign == 0) use the Algoim {phi < 0} rule; interior
            //     cells (sign < 0) use a full Gauss rule so the solid body is
            //     filled with points and no region is left empty.
            //     Row 3 stores the quadrature weight.
            gsMatrix<real_t> quad(4, 0), phys;

            for (auto it = tr_domain.beginBdr(boundary::none);
                 it != tr_domain.endBdr(boundary::none); ++it)
            {
                algoimRule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0)
                {
                    // Algoim failure (re-entrant cell): fall back to uniform grid.
                    const gsVector<real_t> & lo_c = it.lowerCorner();
                    const gsVector<real_t> & hi_c = it.upperCorner();
                    const real_t cellVol = (hi_c - lo_c).prod();
                    const index_t nTotal = nFallback * nFallback * nFallback;
                    gsMatrix<real_t> spts(3, nTotal);
                    index_t scol = 0;
                    for (index_t ii = 0; ii < nFallback; ++ii)
                    for (index_t jj = 0; jj < nFallback; ++jj)
                    for (index_t kk = 0; kk < nFallback; ++kk, ++scol)
                    {
                        spts(0,scol) = lo_c[0] + (ii+0.5)/nFallback*(hi_c[0]-lo_c[0]);
                        spts(1,scol) = lo_c[1] + (jj+0.5)/nFallback*(hi_c[1]-lo_c[1]);
                        spts(2,scol) = lo_c[2] + (kk+0.5)/nFallback*(hi_c[2]-lo_c[2]);
                    }
                    gsMatrix<real_t> fbVals;
                    impl_fun.eval_into(spts, fbVals);
                    std::vector<index_t> insideCols;
                    for (index_t ci = 0; ci < nTotal; ++ci)
                        if (fbVals(0,ci) < 0) insideCols.push_back(ci);
                    if (insideCols.empty()) continue;
                    pts.resize(3, static_cast<index_t>(insideCols.size()));
                    wts.resize(static_cast<index_t>(insideCols.size()));
                    // const real_t w = cellVol / static_cast<real_t>(insideCols.size());
                    const real_t w = cellVol / static_cast<real_t>(nTotal);
                    for (index_t ci = 0; ci < static_cast<index_t>(insideCols.size()); ++ci)
                    {
                        pts.col(ci) = spts.col(insideCols[ci]);
                        wts[ci] = w;
                    }
                }
                mp.patch(0).eval_into(pts, phys);
                const index_t c = quad.cols();
                quad.conservativeResize(4, c + pts.cols());
                quad.block(0, c, 3, pts.cols()) = phys;
                quad.row(3).segment(c, pts.cols()) = wts.transpose();
            }

            gsGaussRule<real_t> gauss(gsVector<index_t, 3>::Constant(
                dbasis.maxCwiseDegree() + 1));
            for (auto it = tr_domain.beginInterior();
                 it != tr_domain.end<InteriorSign>(); ++it)
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
            {
                const std::string base = out + "/" + outputStem + "_quadrature_points";
                gsInfo << "\nExporting " << quad.cols()
                       << " quadrature points to " << base << ".vtp\n";
                gsWriteParaviewPoints(quad, base);
            }

            gsInfo << "ParaView files written to: " << out << "\n";
        }
    }

    return EXIT_SUCCESS;
}
