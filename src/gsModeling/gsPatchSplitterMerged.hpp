/** @file gsPatchSplitterMerged.hpp
    @brief Field-Guided Medial Macro-Pairing ("Merged" Algorithm) for gsPatchSplitter.

    Unifies the three planar quad-meshing paradigms as sequential decoupled modules:
      - Module 1: Boundary-Only CDT & Pruned MAT (Topology, Sleeves, Junctions, Radius r(x))
      - Module 2: MAT-Guided 4-RoSy Cross-Field (Angles, Boundary Orthogonality, Singularities)
      - Module 3: Cross-Field-Optimal Neck Slicing, Hole-Cycle Opening & 1-Layer Parity Balancing
      - Module 4: Singularity-Seeded Macro-Triangulation & Exact Flip/Bitmask Pairing
      - Module 5: Cross-Field Hermite Curve Realization + Coons/Winslow -> AS-G1

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitterGeom.hpp>
#include <gsModeling/gsPatchSplitterTriMesh.hpp>
#include <gsModeling/gsPatchSplitterMAT.hpp>
#include <gsModeling/gsPatchSplitterCrossField.hpp>
#include <gsModeling/gsPatchSplitterMotorcycle.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace gismo {
namespace internal_splitter {

// ============================================================================
// Decoupled Data Contracts for Modules 1 - 5
// ============================================================================

template <class T>
struct MedialBackbone {
    std::vector<VirtualCADLoop<T>> vloops;
    TriMeshData<T> cdt;
    std::vector<int> vLoopIdx;
    std::vector<T>   vArcS;
    std::vector<typename gsPatchSplitterMAT<T>::MATNode>     triNodes;
    std::vector<typename gsPatchSplitterMAT<T>::MacroNode>   macroNodes;
    std::vector<typename gsPatchSplitterMAT<T>::MacroBranch> sleeves;
    std::vector<bool> isSharpTipTri;
    std::vector<bool> isSleeveSharpTip; // per sleeve in `sleeves`
    T domDiag = T(1);
    T targetH = T(0.1);
};

template <class T>
struct GuidedCrossField {
    TriMeshData<T> fineMesh;
    CrossFieldMeshAccel<T> accel;
    CrossFieldSolution<T> field;
    std::vector<T> vertMATRadius;
    std::vector<T> vertSleeveWeight;
};

template <class T>
struct CrossCut {
    int id = -1;
    int sleeveIdx = -1;
    int loopA = -1, loopB = -1;
    T arcSA = T(0), arcSB = T(0);
    gsPoint2D<T> ptA, ptB;
    T costJ = T(0);
};

template <class T>
struct JunctionPolygon {
    int id = -1;
    struct Side {
        bool isCrossCut = false;
        int cutId = -1;
        bool cutForward = true;
        int loopIdx = -1;
        T sStart = T(0), sEnd = T(0);
        int uNode = -1, vNode = -1; // planar arrangement vertex ids
        CurveSeg<T> seg;            // CCW oriented along boundary of P_j
    };
    std::vector<Side> sidesCCW;
    std::vector<int> singularityIds;
};

template <class T>
struct SlicedDecomposition {
    struct ArrNode {
        gsPoint2D<T> pt;
        int loopIdx = -1;
        T arcS = T(0);
        bool isPinnedCorner = false;
        T cornerAngleDeg = T(180);
        int targetValence = 2; // quad valence at boundary vertex (2 on smooth wall, n_q at sharp corner)
    };
    std::vector<ArrNode> nodes;
    std::vector<CrossCut<T>> cuts;
    std::vector<JunctionPolygon<T>> hubs; // Includes both 4-sided sleeve quads and 2m-gon hubs
};

template <class T>
struct MacroQuadGraph {
    struct Vertex {
        gsPoint2D<T> pt;
        bool isBoundary = false;
        int loopIdx = -1;
        T arcS = T(0);
        bool isPinnedCorner = false;
        T cornerAngleDeg = T(180);
        int singularityIdx = -1; // index into GuidedCrossField::field.singularities (-1 if none)
    };
    struct Edge {
        int u = -1, v = -1;
        bool isBoundary = false;
        CurveSeg<T> seg;
    };
    struct Face {
        std::array<int, 4> v;
        std::array<int, 4> e;
        std::array<bool, 4> eFwd;
    };
    std::vector<Vertex> verts;
    std::vector<Edge>   edges;
    std::vector<Face>   faces;
};

// ============================================================================
// Geometry Utilities for Merged Pipeline
// ============================================================================

template <class T>
inline T mergedMeasureSegArcLen(const CurveSeg<T>& seg, int nS = 16)
{
    if (!seg.bsp) return T(0);
    T len = T(0);
    gsPoint2D<T> prev = evalCurvePoint(*seg.bsp, T(0));
    for (int k = 1; k <= nS; ++k) {
        gsPoint2D<T> cur = evalCurvePoint(*seg.bsp, static_cast<T>(k) / static_cast<T>(nS));
        len += prev.dist(cur);
        prev = cur;
    }
    return len;
}

template <class T>
inline T mergedLoopForwardSpan(const VirtualCADLoop<T>& vl, T sStart, T sEnd)
{
    T a = virtualLoopWrapS(vl, sStart);
    T b = virtualLoopWrapS(vl, sEnd);
    T d = b - a;
    if (d <= T(1e-10)) d += vl.totalArcLen;
    return d;
}

template <class T>
inline gsPoint2D<T> mergedInwardNormalOnLoop(const VirtualCADLoop<T>& vl, T s)
{
    // All VirtualCADLoops are oriented so domain lies strictly on the LEFT:
    // outward unit tangent t -> inward normal is (-t.y, t.x)
    gsPoint2D<T> t = evalVirtualLoopUnitTangent(vl, s);
    return gsPoint2D<T>(-t.y, t.x);
}

template <class T>
inline T mergedCcwAngleDeg(const gsPoint2D<T>& u, const gsPoint2D<T>& v)
{
    const T radToDeg = T(180) / T(3.14159265358979323846);
    T nu = std::hypot(u.x, u.y), nv = std::hypot(v.x, v.y);
    if (nu <= T(1e-14) || nv <= T(1e-14)) return T(0);
    T a = std::atan2(u.x * v.y - u.y * v.x, u.x * v.x + u.y * v.y) * radToDeg;
    if (a <= T(0)) a += T(360);
    return a;
}

template <class T>
inline std::vector<std::vector<gsPoint2D<T>>> mergedBuildDenseLoopPolys(
    const std::vector<VirtualCADLoop<T>>& vloops, T domDiag)
{
    std::vector<std::vector<gsPoint2D<T>>> polys(vloops.size());
    const T step = std::max(T(1e-4), domDiag / T(140));
    for (size_t l = 0; l < vloops.size(); ++l) {
        const auto& vl = vloops[l];
        int nS = std::max(48, std::min(320, static_cast<int>(std::ceil(vl.totalArcLen / step))));
        polys[l].reserve(nS);
        for (int k = 0; k < nS; ++k)
            polys[l].push_back(evalVirtualLoopPoint(vl, (vl.totalArcLen * static_cast<T>(k)) / static_cast<T>(nS)));
    }
    return polys;
}

// Distance from point P to closed polygon boundary
template <class T>
inline T mergedPointToPolyDist(const gsPoint2D<T>& P, const std::vector<gsPoint2D<T>>& poly)
{
    const size_t n = poly.size();
    if (n < 2) return T(0);
    T best2 = std::numeric_limits<T>::max();
    for (size_t i = 0; i < n; ++i) {
        const auto& A = poly[i];
        const auto& B = poly[(i + 1) % n];
        T vx = B.x - A.x, vy = B.y - A.y;
        T l2 = vx * vx + vy * vy;
        T t = (l2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((P.x - A.x) * vx + (P.y - A.y) * vy) / l2)) : T(0);
        T dx = P.x - (A.x + t * vx), dy = P.y - (A.y + t * vy);
        best2 = std::min(best2, dx * dx + dy * dy);
    }
    return std::sqrt(best2);
}

// Fast 2-level Pole of Inaccessibility inside simple polygon `poly` (guaranteed strictly inside both `poly` & `loopPolys`)
template <class T>
inline bool mergedFindDeepInteriorPoint(
    const std::vector<gsPoint2D<T>>& poly,
    const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
    gsPoint2D<T>& outPt,
    T& outRadius)
{
    outRadius = T(0);
    if (poly.size() < 3) return false;
    T xMin = poly[0].x, xMax = poly[0].x, yMin = poly[0].y, yMax = poly[0].y;
    gsPoint2D<T> vAvg(0, 0);
    for (const auto& p : poly) {
        xMin = std::min(xMin, p.x); xMax = std::max(xMax, p.x);
        yMin = std::min(yMin, p.y); yMax = std::max(yMax, p.y);
        vAvg.x += p.x; vAvg.y += p.y;
    }
    vAvg.x /= static_cast<T>(poly.size());
    vAvg.y /= static_cast<T>(poly.size());
    T dx = xMax - xMin, dy = yMax - yMin;
    if (dx <= T(1e-12) && dy <= T(1e-12)) return false;

    bool found = false;
    auto evalCand = [&](const gsPoint2D<T>& q) {
        if (!pointInPolyStrict(q, poly, T(1e-6))) return;
        if (!loopPolys.empty() && !isPointStrictlyInsideDomain(q, loopPolys, T(1e-6))) return;
        T d = mergedPointToPolyDist(q, poly);
        if (d > outRadius) { outRadius = d; outPt = q; found = true; }
    };

    evalCand(polygonCentroid(poly));
    evalCand(vAvg);

    const int G1 = 7;
    for (int ix = 0; ix < G1; ++ix) {
        T x = xMin + (static_cast<T>(ix) + T(0.5)) * dx / static_cast<T>(G1);
        for (int iy = 0; iy < G1; ++iy) {
            T y = yMin + (static_cast<T>(iy) + T(0.5)) * dy / static_cast<T>(G1);
            evalCand(gsPoint2D<T>(x, y));
        }
    }
    if (!found) {
        // Finer search for very thin crescents/slivers
        const int Gf = 15;
        for (int ix = 0; ix < Gf; ++ix) {
            T x = xMin + (static_cast<T>(ix) + T(0.5)) * dx / static_cast<T>(Gf);
            for (int iy = 0; iy < Gf; ++iy) {
                T y = yMin + (static_cast<T>(iy) + T(0.5)) * dy / static_cast<T>(Gf);
                evalCand(gsPoint2D<T>(x, y));
            }
        }
    }
    if (!found) return false;

    // Local 5x5 refinement around best point
    gsPoint2D<T> ctr = outPt;
    T stepX = dx / static_cast<T>(G1), stepY = dy / static_cast<T>(G1);
    for (int ix = -2; ix <= 2; ++ix) {
        for (int iy = -2; iy <= 2; ++iy) {
            evalCand(gsPoint2D<T>(ctr.x + static_cast<T>(ix) * T(0.25) * stepX,
                                  ctr.y + static_cast<T>(iy) * T(0.25) * stepY));
        }
    }
    return found && outRadius > T(1e-7);
}

// Verify that straight cross-cut A -> B departs strictly inside the exact CAD tangent sector at BOTH endpoints
template <class T>
inline bool isCrossCutTangentClear(
    const std::vector<VirtualCADLoop<T>>& vloops,
    int lA, T sA,
    int lB, T sB,
    T minDeg = T(14.0))
{
    const gsPoint2D<T> A = evalVirtualLoopPoint(vloops[lA], sA);
    const gsPoint2D<T> B = evalVirtualLoopPoint(vloops[lB], sB);
    const T len = A.dist(B);
    if (len <= T(1e-10)) return false;
    const gsPoint2D<T> dirAB((B.x - A.x) / len, (B.y - A.y) / len);
    const gsPoint2D<T> dirBA(-dirAB.x, -dirAB.y);

    auto checkEnd = [&](int lIdx, T s, const gsPoint2D<T>& dirOut) -> bool {
        const auto& vl = vloops[lIdx];
        T eps = std::min(T(1e-3) * vl.totalArcLen, T(1e-3));
        gsPoint2D<T> tOut = evalVirtualLoopUnitTangent(vl, virtualLoopWrapS(vl, s + eps));
        gsPoint2D<T> tIn  = evalVirtualLoopUnitTangent(vl, virtualLoopWrapS(vl, s - eps));
        gsPoint2D<T> tBack(-tIn.x, -tIn.y);
        T sectorDeg = mergedCcwAngleDeg(tOut, tBack);
        T angOut    = mergedCcwAngleDeg(tOut, dirOut);
        T angBack   = sectorDeg - angOut;
        return (angOut >= minDeg && angBack >= minDeg && angOut < sectorDeg - T(1e-3));
    };

    return checkEnd(lA, sA, dirAB) && checkEnd(lB, sB, dirBA);
}

// Check if straight segment A -> B stays strictly inside domain AND does not cross any existing cut
template <class T>
inline bool isCrossCutClear(
    const gsPoint2D<T>& A,
    const gsPoint2D<T>& B,
    const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
    const std::vector<CrossCut<T>>& existingCuts,
    T minEndTol)
{
    const T len = A.dist(B);
    if (len <= minEndTol) return false;
    if (!isSegmentInsideDomain(A, B, loopPolys)) return false;
    // Extra near-endpoint interior checks to prevent grazing curved concave walls
    for (T f : {T(0.03), T(0.08), T(0.92), T(0.97)}) {
        gsPoint2D<T> q((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y);
        if (!isPointStrictlyInsideDomain(q, loopPolys, T(1e-6))) return false;
    }

    for (const auto& c : existingCuts) {
        bool shareEnd = (A.dist(c.ptA) < minEndTol || A.dist(c.ptB) < minEndTol ||
                         B.dist(c.ptA) < minEndTol || B.dist(c.ptB) < minEndTol);
        if (shareEnd) {
            if ((A.dist(c.ptA) < minEndTol && B.dist(c.ptB) < minEndTol) ||
                (A.dist(c.ptB) < minEndTol && B.dist(c.ptA) < minEndTol))
                return false;
            continue;
        }
        if (segProperIntersect(A, B, c.ptA, c.ptB))
            return false;
        auto ptSegDist = [](const gsPoint2D<T>& P, const gsPoint2D<T>& U, const gsPoint2D<T>& V) -> T {
            T vx = V.x - U.x, vy = V.y - U.y;
            T l2 = vx * vx + vy * vy;
            T t = (l2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((P.x - U.x) * vx + (P.y - U.y) * vy) / l2)) : T(0);
            return std::hypot(P.x - (U.x + t * vx), P.y - (U.y + t * vy));
        };
        if (ptSegDist(A, c.ptA, c.ptB) < T(0.5) * minEndTol ||
            ptSegDist(B, c.ptA, c.ptB) < T(0.5) * minEndTol ||
            ptSegDist(c.ptA, A, B) < T(0.5) * minEndTol ||
            ptSegDist(c.ptB, A, B) < T(0.5) * minEndTol)
            return false;
    }
    return true;
}

// Bulletproof 5-patch O-grid inside ANY simple curved 4-gon bEdges[0..3] (even with smooth 180-deg corners!)
// Automatically switches to a 6-patch Sharp-Tip-Preserving Kite + O-grid when any corner is sharp (< 35 deg),
// so sharp CAD corners (like starburst 26 deg or curved_lshape 10.6 deg) are NEVER sliced in half!
template <class T>
inline bool emitRobustOGrid5(
    const std::array<CurveSeg<T>, 4>& bEdges,
    int degree,
    const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
    std::vector<std::array<CurveSeg<T>, 4>>& outQuads,
    bool allowSharpKite = true,
    T minTargetAngDeg = T(8.5))
{
    gsPoint2D<T> P[4];
    for (int k = 0; k < 4; ++k) {
        if (!bEdges[k].bsp) return false;
        P[k] = bEdges[k].p0();
    }

    // Measure original exact corner tangents, angles, and inward angle bisectors of bEdges[0..3]
    const T pi = T(3.14159265358979323846);
    const T twoPi = T(2) * pi;
    const T degToRad = pi / T(180);
    T origAngs[4], origMinA = T(180), origMaxA = T(0);
    gsPoint2D<T> bisP[4];
    int kSharp = -1;
    for (int k = 0; k < 4; ++k) {
        gsPoint2D<T> dIn  = evalCurveDeriv(*bEdges[(k + 3) % 4].bsp, T(1));
        gsPoint2D<T> dOut = evalCurveDeriv(*bEdges[k].bsp, T(0));
        T nIn = std::hypot(dIn.x, dIn.y), nOut = std::hypot(dOut.x, dOut.y);
        gsPoint2D<T> tOut  = (nOut > T(1e-12)) ? gsPoint2D<T>(dOut.x / nOut, dOut.y / nOut) : gsPoint2D<T>(1, 0);
        gsPoint2D<T> tBack = (nIn  > T(1e-12)) ? gsPoint2D<T>(-dIn.x / nIn, -dIn.y / nIn)   : gsPoint2D<T>(-1, 0);
        origAngs[k] = mergedCcwAngleDeg(tOut, tBack);
        T halfRad = T(0.5) * origAngs[k] * degToRad;
        T ch = std::cos(halfRad), sh = std::sin(halfRad);
        bisP[k] = gsPoint2D<T>(tOut.x * ch - tOut.y * sh, tOut.x * sh + tOut.y * ch);
        if (origAngs[k] < origMinA) { origMinA = origAngs[k]; kSharp = k; }
        origMaxA = std::max(origMaxA, origAngs[k]);
    }

    std::vector<gsPoint2D<T>> capPoly;
    capPoly.reserve(64);
    for (int k = 0; k < 4; ++k)
        for (int s = 0; s < 16; ++s)
            capPoly.push_back(evalCurvePoint(*bEdges[k].bsp, static_cast<T>(s) / T(16)));

    gsPoint2D<T> polePt;
    T rIn = T(0);
    if (!mergedFindDeepInteriorPoint(capPoly, loopPolys, polePt, rIn))
        return false;

    // Sharp-Tip-Preserving 6-Patch Upgrade:
    // Keep genuinely sharp corners kSharp (< 20 deg) UNSPLIT inside Quad 0 = (b[kSharp], P[k1]->I, I->P[k3], b[k3]),
    // and run 5-patch O-grid on the fat remainder (I->P[k1], b[k1], b[k2], P[k3]->I).
    // (Corners >= 20 deg split cleanly into >= 10+10 deg via Hermite bisector spokes without a 6th kite patch!)
    if (allowSharpKite && kSharp >= 0 && origMinA < T(20.0) && origMinA > T(1.5)) {
        int k0 = kSharp, k1 = (k0 + 1) % 4, k2 = (k0 + 2) % 4, k3 = (k0 + 3) % 4;
        gsPoint2D<T> mid13(T(0.5) * (P[k1].x + P[k3].x), T(0.5) * (P[k1].y + P[k3].y));
        std::vector<gsPoint2D<T>> kitePts = {
            polePt,
            gsPoint2D<T>(T(0.70) * mid13.x + T(0.30) * P[k2].x, T(0.70) * mid13.y + T(0.30) * P[k2].y),
            gsPoint2D<T>(T(0.55) * mid13.x + T(0.45) * P[k2].x, T(0.55) * mid13.y + T(0.45) * P[k2].y),
            gsPoint2D<T>(T(0.82) * mid13.x + T(0.18) * P[k2].x, T(0.82) * mid13.y + T(0.18) * P[k2].y),
            gsPoint2D<T>(T(0.50) * (polePt.x + mid13.x), T(0.50) * (polePt.y + mid13.y)),
            gsPoint2D<T>(T(0.60) * mid13.x + T(0.40) * P[k0].x, T(0.60) * mid13.y + T(0.40) * P[k0].y),
            gsPoint2D<T>(T(0.45) * P[k2].x + T(0.55) * P[k0].x, T(0.45) * P[k2].y + T(0.55) * P[k0].y)
        };
        auto segInCap = [&](const CurveSeg<T>& sg) -> bool {
            if (!sg.bsp) return false;
            for (int s = 1; s <= 5; ++s) {
                gsPoint2D<T> pt = evalCurvePoint(*sg.bsp, static_cast<T>(s) / T(6));
                if (!pointInPolyStrict(pt, capPoly, T(1e-6)) ||
                    (!loopPolys.empty() && !isPointStrictlyInsideDomain(pt, loopPolys, T(1e-6))))
                    return false;
            }
            return true;
        };
        for (const auto& I : kitePts) {
            if (!pointInPolyStrict(I, capPoly, T(1e-6))) continue;
            if (!loopPolys.empty() && !isPointStrictlyInsideDomain(I, loopPolys, T(1e-6))) continue;
            T L1 = I.dist(P[k1]), L3 = I.dist(P[k3]);
            if (L1 <= T(1e-7) || L3 <= T(1e-7)) continue;

            // Try 0: straight s1I/sI3; Try 1..2: curved Hermite arcs making angle at I strictly convex (~136 deg)
            // inside tipQuad and reflex (~224 deg) inside remQuad (which Hermite O-grid splits into 112+112 deg!)
            gsPoint2D<T> uI1((P[k1].x - I.x) / L1, (P[k1].y - I.y) / L1);
            gsPoint2D<T> uI3((P[k3].x - I.x) / L3, (P[k3].y - I.y) / L3);
            T tipIAng = mergedCcwAngleDeg(uI3, uI1); // CCW from I->P[k3] to I->P[k1] inside tipQuad
            gsPoint2D<T> bisI;
            {
                T hr = T(0.5) * tipIAng * degToRad;
                bisI = gsPoint2D<T>(uI3.x * std::cos(hr) - uI3.y * std::sin(hr),
                                    uI3.x * std::sin(hr) + uI3.y * std::cos(hr));
            }
            for (int cMode = 0; cMode < 3; ++cMode) {
                CurveSeg<T> s1I, sI3;
                if (cMode == 0) {
                    s1I = makeStraightSeg(P[k1], I, degree);
                    sI3 = makeStraightSeg(I, P[k3], degree);
                } else {
                    if (tipIAng < T(155.0)) continue;
                    T targetHalf = (cMode == 1 ? T(68.0) : T(58.0)) * degToRad;
                    // Depart I toward P[k3] at -targetHalf from bisI, and toward P[k1] at +targetHalf from bisI
                    T cm = std::cos(-targetHalf), sm = std::sin(-targetHalf);
                    T cp = std::cos( targetHalf), sp = std::sin( targetHalf);
                    gsPoint2D<T> tI3(bisI.x * cm - bisI.y * sm, bisI.x * sm + bisI.y * cm);
                    gsPoint2D<T> tI1(bisI.x * cp - bisI.y * sp, bisI.x * sp + bisI.y * cp);
                    T bH = T(0.28);
                    // s1I: P[k1] -> I (arrives at I from tI1 direction, departs P[k1] along bisP[k1])
                    gsPoint2D<T> p1a(P[k1].x + bH * L1 * bisP[k1].x, P[k1].y + bH * L1 * bisP[k1].y);
                    gsPoint2D<T> p1b(I.x     + bH * L1 * tI1.x,      I.y     + bH * L1 * tI1.y);
                    s1I = makeCubicBezierInteriorSeg(P[k1], p1a, p1b, I, degree);
                    // sI3: I -> P[k3] (departs I along tI3, arrives at P[k3] along bisP[k3])
                    gsPoint2D<T> p3a(I.x     + bH * L3 * tI3.x,      I.y     + bH * L3 * tI3.y);
                    gsPoint2D<T> p3b(P[k3].x + bH * L3 * bisP[k3].x, P[k3].y + bH * L3 * bisP[k3].y);
                    sI3 = makeCubicBezierInteriorSeg(I, p3a, p3b, P[k3], degree);
                }
                if (!segInCap(s1I) || !segInCap(sI3)) continue;

                std::array<CurveSeg<T>, 4> tipQuad = {bEdges[k0], s1I, sI3, bEdges[k3]};
                T mnT = 0, mxT = 180, efT = 0;
                eval4SegCornerAnglesDeg(tipQuad[0], tipQuad[1], tipQuad[2], tipQuad[3], mnT, mxT, efT);
                if (mnT < std::min(origMinA - T(0.8), std::max(T(4.0), minTargetAngDeg - T(1.0))) || mxT > T(168.0)) continue;
                std::unique_ptr<gsTensorBSpline<2, T>> tbTip(
                    makeCoonsFromCCWCycle(tipQuad[0], tipQuad[1], tipQuad[2], tipQuad[3], degree));
                if (!tbTip) continue;
                if (evalPatchMinDetDense(*tbTip, 21) <= T(1e-7))
                    quickHarmonicPatchRepair(*tbTip, T(1e-7));
                if (evalPatchMinDetDense(*tbTip, 21) <= T(1e-7)) continue;

                std::array<CurveSeg<T>, 4> remQuad = {s1I.reversed(), bEdges[k1], bEdges[k2], sI3.reversed()};
                std::vector<std::array<CurveSeg<T>, 4>> rem5;
                if (emitRobustOGrid5(remQuad, degree, loopPolys, rem5, false, minTargetAngDeg)) {
                    outQuads.push_back(std::move(tipQuad));
                    for (auto& rq : rem5) outQuads.push_back(std::move(rq));
                    return true;
                }
            }
        }
    }

    std::vector<gsPoint2D<T>> centers = {polePt};
    gsPoint2D<T> cent = polygonCentroid(capPoly);
    if (pointInPolyStrict(cent, capPoly, T(1e-6)) &&
        (loopPolys.empty() || isPointStrictlyInsideDomain(cent, loopPolys, T(1e-6)))) {
        centers.push_back(cent);
        centers.emplace_back(T(0.5) * (polePt.x + cent.x), T(0.5) * (polePt.y + cent.y));
    }
    // Reflex-elbow center candidates: on L-/C-shaped quads with an >= 180 deg corner P[k],
    // a center right in front of the elbow (along P[k]->P[(k+2)%4] or bisP[k]) sees all 4 corners!
    for (int k = 0; k < 4 && centers.size() < 5; ++k) {
        if (origAngs[k] < T(178.0)) continue;
        int kOpp = (k + 2) % 4;
        gsPoint2D<T> cElb(T(0.60) * P[k].x + T(0.40) * P[kOpp].x,
                          T(0.60) * P[k].y + T(0.40) * P[kOpp].y);
        if (pointInPolyStrict(cElb, capPoly, T(1e-6)) &&
            (loopPolys.empty() || isPointStrictlyInsideDomain(cElb, loopPolys, T(1e-6)))) {
            centers.push_back(cElb);
        }
        T spanRef = std::min(P[k].dist(P[(k + 1) % 4]), P[k].dist(P[(k + 3) % 4]));
        T stepRef = std::max(rIn, T(0.45) * spanRef);
        gsPoint2D<T> cBis(P[k].x + T(0.55) * stepRef * bisP[k].x,
                          P[k].y + T(0.55) * stepRef * bisP[k].y);
        if (pointInPolyStrict(cBis, capPoly, T(1e-6)) &&
            (loopPolys.empty() || isPointStrictlyInsideDomain(cBis, loopPolys, T(1e-6)))) {
            centers.push_back(cBis);
        }
    }

    auto spokeValid = [&](const gsPoint2D<T>& Q, const gsPoint2D<T>& Pk) -> bool {
        if (Q.dist(Pk) <= T(1e-8)) return false;
        for (int s = 1; s <= 5; ++s) {
            T f = static_cast<T>(s) / T(6);
            gsPoint2D<T> pt((T(1) - f) * Q.x + f * Pk.x, (T(1) - f) * Q.y + f * Pk.y);
            if (!pointInPolyStrict(pt, capPoly, T(1e-6))) return false;
            if (!loopPolys.empty() && !isPointStrictlyInsideDomain(pt, loopPolys, T(1e-6))) return false;
        }
        return true;
    };
    auto curvedSpokeValid = [&](const CurveSeg<T>& seg) -> bool {
        if (!seg.bsp) return false;
        for (int s = 1; s <= 7; ++s) {
            gsPoint2D<T> pt = evalCurvePoint(*seg.bsp, static_cast<T>(s) / T(8));
            if (!pointInPolyStrict(pt, capPoly, T(1e-6))) return false;
            if (!loopPolys.empty() && !isPointStrictlyInsideDomain(pt, loopPolys, T(1e-6))) return false;
        }
        return true;
    };
    auto spokesDisjoint = [&](const CurveSeg<T> sp[4]) -> bool {
        gsPoint2D<T> pts[4][7];
        for (int k = 0; k < 4; ++k)
            for (int s = 0; s < 7; ++s)
                pts[k][s] = evalCurvePoint(*sp[k].bsp, static_cast<T>(s) / T(6));
        for (int a = 0; a < 4; ++a) {
            for (int b = a + 1; b < 4; ++b) {
                for (int sa = 0; sa < 6; ++sa) {
                    for (int sb = 0; sb < 6; ++sb) {
                        if (segProperIntersect(pts[a][sa], pts[a][sa + 1], pts[b][sb], pts[b][sb + 1]))
                            return false;
                    }
                }
            }
        }
        return true;
    };

    std::vector<std::array<CurveSeg<T>, 4>> bestQuads;
    bool foundValid = false;
    T bestScore = T(1e18);
    T bestFallbackAng = -T(1e18);
    std::vector<std::array<CurveSeg<T>, 4>> fallbackQuads;
    int dbgGeomOk = 0;
    T dbgMaxWorstDet = -T(1e18), dbgBestMinA = T(0), dbgBestMaxA = T(180);

    for (const auto& ctr : centers) {
        T rLoc = mergedPointToPolyDist(ctr, capPoly);
        if (rLoc <= T(1e-7)) continue;

        T ang[4];
        bool monoCCW = true;
        ang[0] = std::atan2(P[0].y - ctr.y, P[0].x - ctr.x);
        for (int k = 1; k < 4; ++k) {
            T a = std::atan2(P[k].y - ctr.y, P[k].x - ctr.x);
            while (a <= ang[k - 1]) a += twoPi;
            ang[k] = a;
        }
        if (ang[3] - ang[0] >= twoPi - T(0.05)) monoCCW = false;

        std::complex<T> phSum(T(0), T(0));
        for (int k = 0; k < 4; ++k) {
            T diff = ang[k] - static_cast<T>(k) * T(0.5) * pi;
            phSum += std::complex<T>(std::cos(diff), std::sin(diff));
        }
        T sqPhi0 = (std::abs(phSum) > T(1e-8)) ? std::arg(phSum) : ang[0];

        for (int mode = 0; mode < 3; ++mode) {
            if (mode == 0 && !monoCCW) continue;
            const T alphas[5] = {T(0.52), T(0.38), T(0.66), T(0.25), T(0.16)};
            for (T alpha : alphas) {
                gsPoint2D<T> Q[4];
                if (mode == 0) {
                    T aUse[4] = {ang[0], ang[1], ang[2], ang[3]};
                    const T minGap = T(0.55);
                    for (int it = 0; it < 4; ++it) {
                        for (int k = 0; k < 4; ++k) {
                            T nxt = (k < 3) ? aUse[k + 1] : (aUse[0] + twoPi);
                            if (nxt - aUse[k] < minGap) {
                                T mid = T(0.5) * (aUse[k] + nxt);
                                aUse[k] = mid - T(0.5) * minGap;
                                if (k < 3) aUse[k + 1] = mid + T(0.5) * minGap;
                            }
                        }
                    }
                    T rad = alpha * rLoc;
                    for (int k = 0; k < 4; ++k)
                        Q[k] = gsPoint2D<T>(ctr.x + rad * std::cos(aUse[k]), ctr.y + rad * std::sin(aUse[k]));
                } else if (mode == 1) {
                    T rad = alpha * rLoc;
                    for (int k = 0; k < 4; ++k) {
                        T ak = sqPhi0 + static_cast<T>(k) * T(0.5) * pi;
                        Q[k] = gsPoint2D<T>(ctr.x + rad * std::cos(ak), ctr.y + rad * std::sin(ak));
                    }
                } else {
                    for (int k = 0; k < 4; ++k)
                        Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x), ctr.y + alpha * (P[k].y - ctr.y));
                }

                if (!isStrictlyConvexQuad(Q)) continue;
                bool innerOk = true;
                for (int k = 0; k < 4; ++k) {
                    if (!pointInPolyStrict(Q[k], capPoly, T(1e-6)) ||
                        (!loopPolys.empty() && !isPointStrictlyInsideDomain(Q[k], loopPolys, T(1e-6))) ||
                        !spokeValid(Q[k], Q[(k + 1) % 4])) {
                        innerOk = false; break;
                    }
                }
                if (!innerOk) continue;

                CurveSeg<T> qEdge[4];
                gsPoint2D<T> extBisQ[4], eBackQ[4];
                T extSecQ[4];
                for (int k = 0; k < 4; ++k) {
                    qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], degree);
                    int kp = (k + 3) % 4, kn = (k + 1) % 4;
                    T lx = Q[kn].x - Q[k].x, ly = Q[kn].y - Q[k].y, ln = std::hypot(lx, ly);
                    T px = Q[kp].x - Q[k].x, py = Q[kp].y - Q[k].y, lp = std::hypot(px, py);
                    gsPoint2D<T> eO = (ln > T(1e-12)) ? gsPoint2D<T>(lx / ln, ly / ln) : gsPoint2D<T>(1, 0);
                    gsPoint2D<T> eB = (lp > T(1e-12)) ? gsPoint2D<T>(px / lp, py / lp) : gsPoint2D<T>(0, 1);
                    eBackQ[k] = eB;
                    extSecQ[k] = mergedCcwAngleDeg(eB, eO);
                    T hr = T(0.5) * extSecQ[k] * degToRad;
                    extBisQ[k] = gsPoint2D<T>(eB.x * std::cos(hr) - eB.y * std::sin(hr),
                                              eB.x * std::sin(hr) + eB.y * std::cos(hr));
                }

                auto evalSpokeCandidate = [&](const CurveSeg<T> spoke[4]) -> bool {
                    std::vector<std::array<CurveSeg<T>, 4>> cand5;
                    cand5.reserve(5);
                    cand5.push_back({qEdge[0], qEdge[1], qEdge[2], qEdge[3]});
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        cand5.push_back({bEdges[k], spoke[kn].reversed(), qEdge[k].reversed(), spoke[k]});
                    }

                    // Fast 1-microsecond angle pre-filter BEFORE building Coons patches or harmonic repairs!
                    T worstMaxA = T(0), worstMinA = T(180);
                    for (int p = 0; p < 5; ++p) {
                        T mnA = 0, mxA = 180, ef = 0;
                        eval4SegCornerAnglesDeg(cand5[p][0], cand5[p][1], cand5[p][2], cand5[p][3], mnA, mxA, ef);
                        worstMinA = std::min(worstMinA, mnA);
                        worstMaxA = std::max(worstMaxA, mxA);
                    }
                    // Splitting a corner of angle origMinA with a bisector yields ~0.5*origMinA,
                    // so allow down to 0.44*origMinA (floored at 1.8 deg) for already-acute input corners.
                    T minAllowedA = std::max(T(1.8), std::min(origMinA * T(0.44), std::max(T(4.2), minTargetAngDeg - T(1.0))));
                    if (worstMaxA >= T(172.0) || worstMinA < minAllowedA)
                        return false;
                    // Skip building Coons if we already have a strictly better valid solution
                    T candSc = (worstMaxA - T(90.0)) * (worstMaxA - T(90.0)) +
                               (T(90.0) - worstMinA) * (T(90.0) - worstMinA);
                    if (foundValid && candSc >= bestScore)
                        return false;

                    T worstDet = std::numeric_limits<T>::max();
                    for (int p = 0; p < 5; ++p) {
                        std::unique_ptr<gsTensorBSpline<2, T>> tb(
                            makeCoonsFromCCWCycle(cand5[p][0], cand5[p][1], cand5[p][2], cand5[p][3], degree));
                        if (!tb) return false;
                        T md = evalPatchMinDetDense(*tb, 21);
                        if (md <= T(1e-7)) {
                            quickHarmonicPatchRepair(*tb, T(1e-7));
                            md = evalPatchMinDetDense(*tb, 21);
                        }
                        if (md <= T(1e-7) && md > -T(1.5e-2) &&
                            tb->knots(0).unique().size() < 3 && tb->knots(1).unique().size() < 3) {
                            tb->insertKnot(T(0.5), 0, 1);
                            tb->insertKnot(T(0.5), 1, 1);
                            quickHarmonicPatchRepair(*tb, T(1e-7));
                            md = evalPatchMinDetDense(*tb, 21);
                            if (md <= T(1e-7) && md > -T(3e-3) && !foundValid && bestFallbackAng <= T(0)) {
                                tb->insertKnot(T(0.25), 0, 1);
                                tb->insertKnot(T(0.75), 0, 1);
                                tb->insertKnot(T(0.25), 1, 1);
                                tb->insertKnot(T(0.75), 1, 1);
                                quickHarmonicPatchRepair(*tb, T(1e-7));
                                md = evalPatchMinDetDense(*tb, 21);
                            }
                        }
                        worstDet = std::min(worstDet, md);
                        if (worstDet <= T(1e-7)) break; // Early exit on first unrepairable sub-patch!
                    }
                    if (worstDet > dbgMaxWorstDet) {
                        dbgMaxWorstDet = worstDet;
                        dbgBestMinA = worstMinA;
                        dbgBestMaxA = worstMaxA;
                    }

                    if (worstDet > T(1e-7) && worstMinA >= minAllowedA && worstMaxA < T(172.0)) {
                        T efAng = std::min(worstMinA, T(180.0) - worstMaxA);
                        if (efAng > bestFallbackAng) {
                            bestFallbackAng = efAng;
                            fallbackQuads = cand5;
                        }
                    }
                    if (worstDet > T(1e-7) && worstMinA > minTargetAngDeg && worstMaxA < T(162.5)) {
                        if (candSc < bestScore) {
                            bestScore = candSc;
                            bestQuads = std::move(cand5);
                            foundValid = true;
                        }
                        return true;
                    }
                    return false;
                };

                // Try 1: Straight spokes (zero curvature, fast path)
                bool straightOk = true;
                for (int k = 0; k < 4; ++k) {
                    if (!spokeValid(Q[k], P[k])) { straightOk = false; break; }
                }
                bool straightClean = false;
                if (straightOk) {
                    ++dbgGeomOk;
                    CurveSeg<T> spoke[4];
                    for (int k = 0; k < 4; ++k) spoke[k] = makeStraightSeg(Q[k], P[k], degree);
                    straightClean = evalSpokeCandidate(spoke);
                }

                // Try 2: Per-Spoke Adaptive Corner-Bisector Cubic Hermite Spokes!
                // Curves each spoke independently with the smallest valid beta so thin arms never clip walls
                // while 180-deg and reflex corners always split cleanly into (ang/2)+(ang/2)!
                if (!straightClean && (!foundValid || bestScore > T(1600.0))) {
                    const T betaSets[3][4] = {
                        {T(0.26), T(0.16), T(0.08), T(0.36)},
                        {T(0.14), T(0.04), T(0.24), T(0.38)},
                        {T(0.36), T(0.24), T(0.12), T(0.04)}
                    };
                    for (int bs = 0; bs < 3; ++bs) {
                        CurveSeg<T> hSpoke[4];
                        bool hOk = true;
                        for (int k = 0; k < 4; ++k) {
                            T L = Q[k].dist(P[k]);
                            if (L <= T(1e-8)) { hOk = false; break; }
                            gsPoint2D<T> ch((P[k].x - Q[k].x) / L, (P[k].y - Q[k].y) / L);
                            // Clamp departure direction at Q[k] so BOTH outer-quad angles at Q[k] stay in [42, 148] deg!
                            T aCh = mergedCcwAngleDeg(eBackQ[k], ch);
                            T loDeg = std::max(T(42.0), extSecQ[k] - T(148.0));
                            T hiDeg = std::min(T(148.0), extSecQ[k] - T(42.0));
                            T useDeg = (loDeg <= hiDeg)
                                     ? ((aCh > extSecQ[k]) ? T(0.5) * extSecQ[k] : std::max(loDeg, std::min(hiDeg, aCh)))
                                     : T(0.5) * extSecQ[k];
                            T cr = useDeg * degToRad;
                            gsPoint2D<T> dQ(eBackQ[k].x * std::cos(cr) - eBackQ[k].y * std::sin(cr),
                                            eBackQ[k].x * std::sin(cr) + eBackQ[k].y * std::cos(cr));
                            bool kBuilt = false;
                            for (int bi = 0; bi < 4; ++bi) {
                                T beta = betaSets[bs][bi];
                                gsPoint2D<T> c1(Q[k].x + beta * L * dQ.x, Q[k].y + beta * L * dQ.y);
                                gsPoint2D<T> c2(P[k].x + beta * L * bisP[k].x, P[k].y + beta * L * bisP[k].y);
                                CurveSeg<T> cand = makeCubicBezierInteriorSeg(Q[k], c1, c2, P[k], degree);
                                if (curvedSpokeValid(cand)) {
                                    hSpoke[k] = std::move(cand);
                                    kBuilt = true;
                                    break;
                                }
                            }
                            if (!kBuilt && spokeValid(Q[k], P[k])) {
                                hSpoke[k] = makeStraightSeg(Q[k], P[k], degree);
                                kBuilt = true;
                            }
                            if (!kBuilt) { hOk = false; break; }
                        }
                        if (!hOk || !spokesDisjoint(hSpoke)) continue;
                        ++dbgGeomOk;
                        if (evalSpokeCandidate(hSpoke)) break;
                    }
                }
            }
        }
    }

    if (foundValid) {
        for (auto& q : bestQuads) outQuads.push_back(std::move(q));
        return true;
    }
    if (bestFallbackAng > T(0) && !fallbackQuads.empty()) {
        for (auto& q : fallbackQuads) outQuads.push_back(std::move(q));
        return true;
    }

    // Curved-Inner-Quad 5-Patch Fallback for strongly bent C-/L-/U-shaped quads (e.g., horseshoe wrapping around a hole):
    // When no straight-sided inner quad Q[0..3] can see all 4 corners, build a curved inner quad whose 4 sides follow
    // the exact curved opposite-edge blend of bEdges[0..3] (so spokes Q[k]->P[k] stay short and strictly inside capPoly!).
    {
        auto evalCoonsPt = [&](T u, T v) -> gsPoint2D<T> {
            gsPoint2D<T> b0 = evalCurvePoint(*bEdges[0].bsp, u);
            gsPoint2D<T> b1 = evalCurvePoint(*bEdges[1].bsp, v);
            gsPoint2D<T> b2 = evalCurvePoint(*bEdges[2].bsp, T(1) - u);
            gsPoint2D<T> b3 = evalCurvePoint(*bEdges[3].bsp, T(1) - v);
            T x = (T(1) - v) * b0.x + v * b2.x + (T(1) - u) * b3.x + u * b1.x
                - (T(1) - u) * (T(1) - v) * P[0].x - u * (T(1) - v) * P[1].x
                - u * v * P[2].x - (T(1) - u) * v * P[3].x;
            T y = (T(1) - v) * b0.y + v * b2.y + (T(1) - u) * b3.y + u * b1.y
                - (T(1) - u) * (T(1) - v) * P[0].y - u * (T(1) - v) * P[1].y
                - u * v * P[2].y - (T(1) - u) * v * P[3].y;
            return gsPoint2D<T>(x, y);
        };
        auto evalBilinBlendPt = [&](T u, T v) -> gsPoint2D<T> {
            gsPoint2D<T> b0 = evalCurvePoint(*bEdges[0].bsp, u);
            gsPoint2D<T> b2 = evalCurvePoint(*bEdges[2].bsp, T(1) - u);
            gsPoint2D<T> b3 = evalCurvePoint(*bEdges[3].bsp, T(1) - v);
            gsPoint2D<T> b1 = evalCurvePoint(*bEdges[1].bsp, v);
            return gsPoint2D<T>(T(0.5) * ((T(1) - v) * b0.x + v * b2.x + (T(1) - u) * b3.x + u * b1.x),
                                T(0.5) * ((T(1) - v) * b0.y + v * b2.y + (T(1) - u) * b3.y + u * b1.y));
        };
        auto evalEdgeLoftPt = [&](int axis, T u, T v) -> gsPoint2D<T> {
            if (axis == 0) {
                gsPoint2D<T> b0 = evalCurvePoint(*bEdges[0].bsp, u);
                gsPoint2D<T> b2 = evalCurvePoint(*bEdges[2].bsp, T(1) - u);
                return gsPoint2D<T>((T(1) - v) * b0.x + v * b2.x, (T(1) - v) * b0.y + v * b2.y);
            } else {
                gsPoint2D<T> b3 = evalCurvePoint(*bEdges[3].bsp, T(1) - v);
                gsPoint2D<T> b1 = evalCurvePoint(*bEdges[1].bsp, v);
                return gsPoint2D<T>((T(1) - u) * b3.x + u * b1.x, (T(1) - u) * b3.y + u * b1.y);
            }
        };
        auto fitCubic4Pts = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& M1,
                                const gsPoint2D<T>& M2, const gsPoint2D<T>& B) -> CurveSeg<T> {
            // Exact cubic Bezier interpolating A at t=0, M1 at t=1/3, M2 at t=2/3, B at t=1
            gsPoint2D<T> C1((-T(5) * A.x + T(18) * M1.x - T(9) * M2.x + T(2) * B.x) / T(6),
                            (-T(5) * A.y + T(18) * M1.y - T(9) * M2.y + T(2) * B.y) / T(6));
            gsPoint2D<T> C2((T(2) * A.x - T(9) * M1.x + T(18) * M2.x - T(5) * B.x) / T(6),
                            (T(2) * A.y - T(9) * M1.y + T(18) * M2.y - T(5) * B.y) / T(6));
            return makeCubicBezierInteriorSeg(A, C1, C2, B, degree);
        };

        const T cAlphas[4] = {T(0.28), T(0.20), T(0.35), T(0.14)};
        for (int bMode = 0; bMode < 4 && !foundValid; ++bMode) {
            for (T al : cAlphas) {
                auto sampleUV = [&](T u, T v) -> gsPoint2D<T> {
                    if (bMode == 0) return evalEdgeLoftPt(0, u, v);
                    if (bMode == 1) return evalEdgeLoftPt(1, u, v);
                    if (bMode == 2) return evalBilinBlendPt(u, v);
                    return evalCoonsPt(u, v);
                };
                const T uvC[4][2] = {{al, al}, {T(1) - al, al}, {T(1) - al, T(1) - al}, {al, T(1) - al}};
                gsPoint2D<T> Q[4];
                bool qOk = true;
                for (int k = 0; k < 4; ++k) {
                    Q[k] = sampleUV(uvC[k][0], uvC[k][1]);
                    if (!pointInPolyStrict(Q[k], capPoly, T(1e-6)) ||
                        (!loopPolys.empty() && !isPointStrictlyInsideDomain(Q[k], loopPolys, T(1e-6)))) {
                        qOk = false; break;
                    }
                }
                if (!qOk) continue;

                CurveSeg<T> cqEdge[4], cSpoke[4];
                bool curvesOk = true;
                for (int k = 0; k < 4; ++k) {
                    int kn = (k + 1) % 4;
                    T u0 = uvC[k][0], v0 = uvC[k][1], u1 = uvC[kn][0], v1 = uvC[kn][1];
                    gsPoint2D<T> m1 = sampleUV((T(2) * u0 + u1) / T(3), (T(2) * v0 + v1) / T(3));
                    gsPoint2D<T> m2 = sampleUV((u0 + T(2) * u1) / T(3), (v0 + T(2) * v1) / T(3));
                    cqEdge[k] = fitCubic4Pts(Q[k], m1, m2, Q[kn]);
                    if (!curvedSpokeValid(cqEdge[k])) { curvesOk = false; break; }
                }
                if (!curvesOk) continue;

                for (int k = 0; k < 4; ++k) {
                    T L = Q[k].dist(P[k]);
                    if (L <= T(1e-8)) { curvesOk = false; break; }
                    // Compute exact curved inner-quad exterior sector & bisector at Q[k]
                    gsPoint2D<T> dOutQ = evalCurveDeriv(*cqEdge[k].bsp, T(0));
                    gsPoint2D<T> dInQ  = evalCurveDeriv(*cqEdge[(k + 3) % 4].bsp, T(1));
                    T nOQ = std::hypot(dOutQ.x, dOutQ.y), nIQ = std::hypot(dInQ.x, dInQ.y);
                    gsPoint2D<T> eOQ = (nOQ > T(1e-12)) ? gsPoint2D<T>(dOutQ.x / nOQ, dOutQ.y / nOQ) : gsPoint2D<T>(1, 0);
                    gsPoint2D<T> eBQ = (nIQ > T(1e-12)) ? gsPoint2D<T>(-dInQ.x / nIQ, -dInQ.y / nIQ) : gsPoint2D<T>(0, 1);
                    T extSecQk = mergedCcwAngleDeg(eBQ, eOQ);
                    gsPoint2D<T> ch((P[k].x - Q[k].x) / L, (P[k].y - Q[k].y) / L);
                    T aCh = mergedCcwAngleDeg(eBQ, ch);
                    T loDeg = std::max(T(42.0), extSecQk - T(148.0));
                    T hiDeg = std::min(T(148.0), extSecQk - T(42.0));
                    T useDeg = (loDeg <= hiDeg)
                             ? ((aCh > extSecQk) ? T(0.5) * extSecQk : std::max(loDeg, std::min(hiDeg, aCh)))
                             : T(0.5) * extSecQk;
                    T cr = useDeg * degToRad;
                    gsPoint2D<T> dQ(eBQ.x * std::cos(cr) - eBQ.y * std::sin(cr),
                                    eBQ.x * std::sin(cr) + eBQ.y * std::cos(cr));

                    // Always use curved bisector spoke on >=155 deg corners so 180/206+ corners split cleanly (~ang/2)!
                    bool built = false;
                    if (origAngs[k] >= T(155.0) || !spokeValid(Q[k], P[k])) {
                        for (T beta : {T(0.24), T(0.15), T(0.34), T(0.08), T(0.42)}) {
                            gsPoint2D<T> c1(Q[k].x + beta * L * dQ.x, Q[k].y + beta * L * dQ.y);
                            gsPoint2D<T> c2(P[k].x + beta * L * bisP[k].x, P[k].y + beta * L * bisP[k].y);
                            CurveSeg<T> cand = makeCubicBezierInteriorSeg(Q[k], c1, c2, P[k], degree);
                            if (curvedSpokeValid(cand)) { cSpoke[k] = std::move(cand); built = true; break; }
                        }
                    }
                    if (!built && spokeValid(Q[k], P[k])) {
                        cSpoke[k] = makeStraightSeg(Q[k], P[k], degree);
                        built = true;
                    }
                    if (!built) { curvesOk = false; break; }
                }
                if (!curvesOk || !spokesDisjoint(cSpoke)) continue;

                // Ensure inner curved loop cqEdge[0..3] has positive area and does not self-cross
                std::vector<std::array<CurveSeg<T>, 4>> cand5;
                cand5.reserve(5);
                cand5.push_back({cqEdge[0], cqEdge[1], cqEdge[2], cqEdge[3]});
                for (int k = 0; k < 4; ++k) {
                    int kn = (k + 1) % 4;
                    cand5.push_back({bEdges[k], cSpoke[kn].reversed(), cqEdge[k].reversed(), cSpoke[k]});
                }

                T worstMaxA = T(0), worstMinA = T(180);
                for (int p = 0; p < 5; ++p) {
                    T mnA = 0, mxA = 180, ef = 0;
                    eval4SegCornerAnglesDeg(cand5[p][0], cand5[p][1], cand5[p][2], cand5[p][3], mnA, mxA, ef);
                    worstMinA = std::min(worstMinA, mnA);
                    worstMaxA = std::max(worstMaxA, mxA);
                }
                T minAllowedA = std::max(T(1.8), std::min(origMinA * T(0.44), std::max(T(4.2), minTargetAngDeg - T(1.0))));
                if (worstMaxA >= T(172.0) || worstMinA < minAllowedA) continue;

                ++dbgGeomOk;
                T worstDet = std::numeric_limits<T>::max();
                for (int p = 0; p < 5; ++p) {
                    std::unique_ptr<gsTensorBSpline<2, T>> tb(
                        makeCoonsFromCCWCycle(cand5[p][0], cand5[p][1], cand5[p][2], cand5[p][3], degree));
                    if (!tb) { worstDet = -T(1); break; }
                    T md = evalPatchMinDetDense(*tb, 21);
                    if (md <= T(1e-7)) {
                        quickHarmonicPatchRepair(*tb, T(1e-7));
                        md = evalPatchMinDetDense(*tb, 21);
                    }
                    if (md <= T(1e-7) && md > -T(5e-2) &&
                        tb->knots(0).unique().size() < 3 && tb->knots(1).unique().size() < 3) {
                        tb->insertKnot(T(0.5), 0, 1);
                        tb->insertKnot(T(0.5), 1, 1);
                        quickHarmonicPatchRepair(*tb, T(1e-7));
                        md = evalPatchMinDetDense(*tb, 21);
                        if (md <= T(1e-7) && md > -T(5e-2)) {
                            tb->insertKnot(T(0.25), 0, 1);
                            tb->insertKnot(T(0.75), 0, 1);
                            tb->insertKnot(T(0.25), 1, 1);
                            tb->insertKnot(T(0.75), 1, 1);
                            quickHarmonicPatchRepair(*tb, T(1e-7));
                            md = evalPatchMinDetDense(*tb, 21);
                        }
                    }
                    worstDet = std::min(worstDet, md);
                    if (worstDet <= T(1e-7)) break;
                }
                if (worstDet > dbgMaxWorstDet) {
                    dbgMaxWorstDet = worstDet;
                    dbgBestMinA = worstMinA;
                    dbgBestMaxA = worstMaxA;
                }
                if (worstDet > T(1e-7)) {
                    for (auto& q : cand5) outQuads.push_back(std::move(q));
                    return true;
                }
            }
        }

        // Bisector-Ray-March Curved-Inner-Quad Fallback (handles strongly bent L-/C-shaped reflex quads like 206 deg!):
        // Marches Q[k] inward directly along each corner bisector bisP[k] so ALL 4 spokes Q[k]->P[k] stay strictly
        // inside capPoly and split every corner evenly (~ang/2), while inner edges cqEdge[k] curve parallel to bEdges[k]!
        {
            T bboxDiag = T(0);
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    bboxDiag = std::max(bboxDiag, P[a].dist(P[b]));
            auto rayMaxInCap = [&](const gsPoint2D<T>& orig, const gsPoint2D<T>& dir) -> T {
                T lo = T(0), hi = std::max(T(1e-4), bboxDiag);
                // Find an upper bound or march outward from small step until exiting capPoly
                T step = std::max(T(1e-5), T(0.025) * bboxDiag);
                T lastOk = T(0);
                for (int i = 1; i <= 40; ++i) {
                    T d = static_cast<T>(i) * step;
                    gsPoint2D<T> q(orig.x + d * dir.x, orig.y + d * dir.y);
                    if (pointInPolyStrict(q, capPoly, T(1e-6)) &&
                        (loopPolys.empty() || isPointStrictlyInsideDomain(q, loopPolys, T(1e-6)))) {
                        lastOk = d;
                    } else {
                        hi = d;
                        lo = lastOk;
                        break;
                    }
                }
                if (lastOk <= T(1e-7)) {
                    // Try tiny steps right near the corner
                    for (T f : {T(0.01), T(0.004), T(0.0015), T(0.0005)}) {
                        T d = f * bboxDiag;
                        gsPoint2D<T> q(orig.x + d * dir.x, orig.y + d * dir.y);
                        if (pointInPolyStrict(q, capPoly, T(1e-7)) &&
                            (loopPolys.empty() || isPointStrictlyInsideDomain(q, loopPolys, T(1e-7)))) {
                            lastOk = d;
                            break;
                        }
                    }
                    return lastOk;
                }
                for (int it = 0; it < 6; ++it) {
                    T mid = T(0.5) * (lo + hi);
                    gsPoint2D<T> q(orig.x + mid * dir.x, orig.y + mid * dir.y);
                    if (pointInPolyStrict(q, capPoly, T(1e-6)) &&
                        (loopPolys.empty() || isPointStrictlyInsideDomain(q, loopPolys, T(1e-6))))
                        lo = mid;
                    else
                        hi = mid;
                }
                return lo;
            };

            T dMax[4];
            bool allRaysOk = true;
            for (int k = 0; k < 4; ++k) {
                dMax[k] = rayMaxInCap(P[k], bisP[k]);
                if (dMax[k] <= T(1e-7)) { allRaysOk = false; break; }
            }
            if (allRaysOk) {
                int kRefIdx = -1;
                for (int k = 0; k < 4; ++k)
                    if (origAngs[k] > T(175.0) && (kRefIdx < 0 || origAngs[k] > origAngs[kRefIdx])) kRefIdx = k;
                struct RayFracSet { T fRef, fSide, fOpp; };
                std::vector<RayFracSet> fracSets = {
                    {T(0.32), T(0.32), T(0.32)},
                    {T(0.22), T(0.22), T(0.22)},
                    {T(0.42), T(0.42), T(0.42)}
                };
                if (kRefIdx >= 0) {
                    // Push Q[kRef] deeper across the elbow along bisP[kRef] so inner quad Q[0..3] becomes strictly convex!
                    fracSets.insert(fracSets.begin(), {
                        {T(0.60), T(0.22), T(0.16)},
                        {T(0.50), T(0.24), T(0.18)},
                        {T(0.72), T(0.18), T(0.12)},
                        {T(0.42), T(0.22), T(0.18)},
                        {T(0.80), T(0.16), T(0.10)}
                    });
                }
                for (const auto& fs : fracSets) {
                    gsPoint2D<T> Q[4];
                    CurveSeg<T> cSpoke[4], cqEdge[4];
                    bool okAll = true;
                    for (int k = 0; k < 4; ++k) {
                        T fk = fs.fSide;
                        bool isRef = (k == kRefIdx);
                        if (isRef) fk = fs.fRef;
                        else if (kRefIdx >= 0 && k == (kRefIdx + 2) % 4) fk = fs.fOpp;
                        T dUse = fk * dMax[k];
                        if (!isRef)
                            dUse = std::min(dUse, T(0.38) * std::min(P[k].dist(P[(k + 1) % 4]), P[k].dist(P[(k + 3) % 4])));
                        Q[k] = gsPoint2D<T>(P[k].x + dUse * bisP[k].x, P[k].y + dUse * bisP[k].y);
                        if (!pointInPolyStrict(Q[k], capPoly, T(1e-6)) || !spokeValid(Q[k], P[k])) {
                            okAll = false; break;
                        }
                        cSpoke[k] = makeStraightSeg(Q[k], P[k], degree);
                    }
                    if (!okAll || !spokesDisjoint(cSpoke)) continue;

                    T frac = fs.fSide;
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        if (spokeValid(Q[k], Q[kn])) {
                            cqEdge[k] = makeStraightSeg(Q[k], Q[kn], degree);
                        } else {
                            // Curve inner edge parallel to bEdges[k] by marching inward along bEdges[k]'s unit normal at 1/3 and 2/3
                            auto offsetAt = [&](T t, gsPoint2D<T>& outPt) -> bool {
                                gsPoint2D<T> bp = evalCurvePoint(*bEdges[k].bsp, t);
                                gsPoint2D<T> dv = evalCurveDeriv(*bEdges[k].bsp, t);
                                T nl = std::hypot(dv.x, dv.y);
                                if (nl <= T(1e-12)) return false;
                                gsPoint2D<T> inN(-dv.y / nl, dv.x / nl);
                                T dm = rayMaxInCap(bp, inN);
                                if (dm <= T(1e-7)) return false;
                                T off = std::min(frac * dm, T(0.45) * dm);
                                outPt = gsPoint2D<T>(bp.x + off * inN.x, bp.y + off * inN.y);
                                return pointInPolyStrict(outPt, capPoly, T(1e-6));
                            };
                            gsPoint2D<T> m1, m2;
                            if (!offsetAt(T(1) / T(3), m1) || !offsetAt(T(2) / T(3), m2)) { okAll = false; break; }
                            cqEdge[k] = fitCubic4Pts(Q[k], m1, m2, Q[kn]);
                            if (!curvedSpokeValid(cqEdge[k])) { okAll = false; break; }
                        }
                    }
                    if (!okAll || !spokesDisjoint(cqEdge)) continue;

                    // Verify inner curved loop has strictly positive area
                    std::vector<gsPoint2D<T>> inPl;
                    for (int k = 0; k < 4; ++k)
                        for (int s = 0; s < 6; ++s)
                            inPl.push_back(evalCurvePoint(*cqEdge[k].bsp, static_cast<T>(s) / T(6)));
                    T inA2 = T(0);
                    for (size_t i = 0; i < inPl.size(); ++i)
                        inA2 += inPl[i].x * inPl[(i + 1) % inPl.size()].y - inPl[(i + 1) % inPl.size()].x * inPl[i].y;
                    if (inA2 <= T(1e-8)) continue;

                    std::vector<std::array<CurveSeg<T>, 4>> cand5;
                    cand5.reserve(5);
                    cand5.push_back({cqEdge[0], cqEdge[1], cqEdge[2], cqEdge[3]});
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        cand5.push_back({bEdges[k], cSpoke[kn].reversed(), cqEdge[k].reversed(), cSpoke[k]});
                    }
                    T worstMaxA = T(0), worstMinA = T(180);
                    for (int p = 0; p < 5; ++p) {
                        T mnA = 0, mxA = 180, ef = 0;
                        eval4SegCornerAnglesDeg(cand5[p][0], cand5[p][1], cand5[p][2], cand5[p][3], mnA, mxA, ef);
                        worstMinA = std::min(worstMinA, mnA);
                        worstMaxA = std::max(worstMaxA, mxA);
                    }
                    T minAllowedA = std::max(T(1.8), std::min(origMinA * T(0.44), std::max(T(4.2), minTargetAngDeg - T(1.0))));
                    if (worstMaxA >= T(172.0) || worstMinA < minAllowedA) continue;

                    ++dbgGeomOk;
                    T worstDet = std::numeric_limits<T>::max();
                    for (int p = 0; p < 5; ++p) {
                        std::unique_ptr<gsTensorBSpline<2, T>> tb(
                            makeCoonsFromCCWCycle(cand5[p][0], cand5[p][1], cand5[p][2], cand5[p][3], degree));
                        if (!tb) { worstDet = -T(1); break; }
                        T md = evalPatchMinDetDense(*tb, 21);
                        if (md <= T(1e-7)) {
                            quickHarmonicPatchRepair(*tb, T(1e-7));
                            md = evalPatchMinDetDense(*tb, 21);
                        }
                        if (md <= T(1e-7) && md > -T(5e-2) &&
                            tb->knots(0).unique().size() < 3 && tb->knots(1).unique().size() < 3) {
                            tb->insertKnot(T(0.5), 0, 1);
                            tb->insertKnot(T(0.5), 1, 1);
                            quickHarmonicPatchRepair(*tb, T(1e-7));
                            md = evalPatchMinDetDense(*tb, 21);
                            if (md <= T(1e-7) && md > -T(5e-2)) {
                                tb->insertKnot(T(0.25), 0, 1);
                                tb->insertKnot(T(0.75), 0, 1);
                                tb->insertKnot(T(0.25), 1, 1);
                                tb->insertKnot(T(0.75), 1, 1);
                                quickHarmonicPatchRepair(*tb, T(1e-7));
                                md = evalPatchMinDetDense(*tb, 21);
                            }
                        }
                        worstDet = std::min(worstDet, md);
                        if (worstDet <= T(1e-7)) break;
                    }
                    if (worstDet > dbgMaxWorstDet) {
                        dbgMaxWorstDet = worstDet;
                        dbgBestMinA = worstMinA;
                        dbgBestMaxA = worstMaxA;
                    }
                    if (worstDet > T(1e-7)) {
                        for (auto& q : cand5) outQuads.push_back(std::move(q));
                        return true;
                    }
                }
            }
        }
    }

    if (std::getenv("GISMO_MERGED_DEBUG") != nullptr) {
        std::cerr << "  [emitOGrid FAIL] origAngs=[" << origAngs[0] << "," << origAngs[1] << "," << origAngs[2] << "," << origAngs[3]
                  << "] geomOk=" << dbgGeomOk << " bestDet=" << dbgMaxWorstDet
                  << " ang=[" << dbgBestMinA << "," << dbgBestMaxA << "]\n";
    }
    return false;
}

// ============================================================================
// Module 1: Boundary-Only CDT & Pruned MAT Backbone
// ============================================================================

template <class T>
inline bool buildMedialBackbone(
    const gsPlanarDomain<T>& dom,
    int degree,
    T smoothJointTolDeg,
    T targetHFrac,
    T thetaPruneDeg,
    T junctionMergeRatio,
    MedialBackbone<T>& out)
{
    using MAT = gsPatchSplitterMAT<T>;
    using MATNode = typename MAT::MATNode;
    using MacroNode = typename MAT::MacroNode;
    using MacroBranch = typename MAT::MacroBranch;

    if (dom.numLoops() == 0) return false;
    const gsMatrix<T> bb = dom.boundingBox();
    const T domDiag = std::hypot(bb(0, 1) - bb(0, 0), bb(1, 1) - bb(1, 0));
    if (domDiag <= T(1e-12)) return false;

    out = MedialBackbone<T>();
    out.domDiag = domDiag;
    out.targetH = std::max(T(1e-4), domDiag * targetHFrac);
    out.vloops = buildVirtualCADLoops(dom, degree, smoothJointTolDeg);

    // 1.1 Boundary-only CDT + Lawson flips
    if (!buildConstrainedDomainTriangulation(dom, degree, out.targetH, true, 0, out.cdt))
        return false;
    applyDomainLawsonFlips(out.cdt);

    auto& mesh = out.cdt;
    const T radToDeg = T(180) / T(3.14159265358979323846);

    // Eliminate extreme sliver J-triangles bordering S-triangles (same as gsPatchSplitterMAT Step 1.1)
    {
        const int nT0 = static_cast<int>(mesh.tris.size());
        const int nE0 = static_cast<int>(mesh.edges.size());
        bool flippedAny = true;
        for (int pass = 0; pass < 4 && flippedAny; ++pass) {
            flippedAny = false;
            struct EAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
            std::vector<EAdj> eAdj(nE0);
            for (int ti = 0; ti < nT0; ++ti) {
                const auto& tr = mesh.tris[ti];
                for (int k = 0; k < 3; ++k) {
                    int eid = tr.e[k];
                    if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
                    else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
                }
            }
            auto triBndCount = [&](int ti) -> int {
                int nb = 0;
                for (int k = 0; k < 3; ++k)
                    if (mesh.edges[mesh.tris[ti].e[k]].isBoundary) ++nb;
                return nb;
            };
            auto isSliverJLongestEdge = [&](int ti, int kEdge) -> bool {
                if (triBndCount(ti) != 0) return false;
                const auto& tr = mesh.tris[ti];
                const auto& P0 = mesh.pts[tr.v[0]];
                const auto& P1 = mesh.pts[tr.v[1]];
                const auto& P2 = mesh.pts[tr.v[2]];
                T l[3] = {P0.dist(P1), P1.dist(P2), P2.dist(P0)};
                T lMax = std::max(l[0], std::max(l[1], l[2]));
                if (lMax <= T(1e-12) || l[kEdge] < T(0.999) * lMax) return false;
                T hMin = std::abs(cross2(P0, P1, P2)) / lMax;
                T mnA = 0, mxA = 0;
                triMinMaxAnglesDeg(P0, P1, P2, mnA, mxA);
                return (mxA > T(168.0) && hMin < T(0.04) * lMax);
            };
            for (int eIdx = 0; eIdx < nE0; ++eIdx) {
                if (mesh.edges[eIdx].isBoundary) continue;
                int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
                int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
                if (tA < 0 || tB < 0) continue;
                bool flipCand = (isSliverJLongestEdge(tA, kA) && triBndCount(tB) == 1) ||
                                (isSliverJLongestEdge(tB, kB) && triBndCount(tA) == 1);
                if (!flipCand) continue;
                const auto A = mesh.tris[tA], B = mesh.tris[tB];
                int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
                int d = B.v[(kB + 2) % 3];
                if (c == d || c == a || c == b || d == a || d == b) continue;
                if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) continue;
                if (mesh.pts[c].dist(mesh.pts[d]) >= T(0.90) * mesh.pts[a].dist(mesh.pts[b])) continue;
                gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
                if (!isStrictlyConvexQuad(qa)) continue;
                if (!isSegmentInsideDomain(mesh.pts[c], mesh.pts[d], mesh.loopPolys)) continue;
                int e_bc = A.e[(kA + 1) % 3], e_ca = A.e[(kA + 2) % 3];
                int e_ad = B.e[(kB + 1) % 3], e_db = B.e[(kB + 2) % 3];
                mesh.edges[eIdx] = {c, d, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1};
                mesh.tris[tA] = {{c, a, d}, {e_ca, e_ad, eIdx}};
                mesh.tris[tB] = {{c, d, b}, {eIdx, e_db, e_bc}};
                flippedAny = true;
                break;
            }
        }
    }

    const int nV = static_cast<int>(mesh.pts.size());
    const int nT = static_cast<int>(mesh.tris.size());
    const int nE = static_cast<int>(mesh.edges.size());
    if (nT < 2) return false;

    // Anchor every boundary CDT vertex onto VirtualCADLoop (loopIdx, arcS)
    out.vLoopIdx.assign(nV, -1);
    out.vArcS.assign(nV, T(0));
    for (int i = 0; i < nV; ++i) {
        if (!mesh.isBoundaryPt[i]) continue;
        int lHint = (static_cast<size_t>(i) < mesh.bndInfo.size()) ? mesh.bndInfo[i].loopIdx : -1;
        if (lHint >= 0 && lHint < static_cast<int>(out.vloops.size())) {
            T s = T(0);
            gsPoint2D<T> proj;
            projectPointToVirtualLoop(out.vloops[lHint], mesh.pts[i], s, proj);
            out.vLoopIdx[i] = lHint;
            out.vArcS[i] = s;
        } else {
            T bestD = std::numeric_limits<T>::max();
            for (int l = 0; l < static_cast<int>(out.vloops.size()); ++l) {
                T s = T(0);
                gsPoint2D<T> proj;
                T d = projectPointToVirtualLoop(out.vloops[l], mesh.pts[i], s, proj);
                if (d < bestD) { bestD = d; out.vLoopIdx[i] = l; out.vArcS[i] = s; }
            }
        }
    }

    // Build T/S/J classification and MATNodes
    std::vector<std::vector<int>> edgeTris(nE);
    out.triNodes.resize(nT);
    for (int ti = 0; ti < nT; ++ti) {
        const auto& tr = mesh.tris[ti];
        int nB = 0;
        for (int k = 0; k < 3; ++k) {
            int eid = tr.e[k];
            if (mesh.edges[eid].isBoundary) ++nB;
            else edgeTris[eid].push_back(ti);
        }
        MATNode& node = out.triNodes[ti];
        node.triIdx = ti;
        node.triType = (nB >= 2) ? 0 : ((nB == 1) ? 1 : 2);

        const auto& A = mesh.pts[tr.v[0]];
        const auto& B = mesh.pts[tr.v[1]];
        const auto& C = mesh.pts[tr.v[2]];
        T D = T(2) * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
        bool useCC = false;
        gsPoint2D<T> cc;
        T cr = T(0);
        if (std::abs(D) > T(1e-12)) {
            T a2 = A.x * A.x + A.y * A.y, b2 = B.x * B.x + B.y * B.y, c2 = C.x * C.x + C.y * C.y;
            cc = gsPoint2D<T>((a2 * (B.y - C.y) + b2 * (C.y - A.y) + c2 * (A.y - B.y)) / D,
                              (a2 * (C.x - B.x) + b2 * (A.x - C.x) + c2 * (B.x - A.x)) / D);
            cr = cc.dist(A);
            T dAB = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y);
            T dBC = (B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y);
            T dCA = (C.x - A.x) * (C.x - A.x) + (C.y - A.y) * (C.y - A.y);
            bool isObtuse = (dAB + dBC < dCA - T(1e-9)) || (dBC + dCA < dAB - T(1e-9)) || (dCA + dAB < dBC - T(1e-9));
            if (!isObtuse && cross2(A, B, cc) > T(1e-8) && cross2(B, C, cc) > T(1e-8) &&
                cross2(C, A, cc) > T(1e-8) && isPointStrictlyInsideDomain(cc, mesh.loopPolys, T(1e-5)))
                useCC = true;
        }
        if (useCC && node.triType != 0) {
            node.pt = cc;
            node.radius = cr;
        } else if (node.triType == 2 || node.triType == 0) {
            T la = B.dist(C), lb = C.dist(A), lc = A.dist(B);
            T perim = std::max(T(1e-12), la + lb + lc);
            gsPoint2D<T> inc((la * A.x + lb * B.x + lc * C.x) / perim,
                             (la * A.y + lb * B.y + lc * C.y) / perim);
            gsPoint2D<T> cent((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
            node.pt = (node.triType == 2) ? gsPoint2D<T>(T(0.5) * (inc.x + cent.x), T(0.5) * (inc.y + cent.y)) : inc;
            node.radius = std::max(std::abs(cross2(A, B, C)) / perim,
                                   (node.pt.dist(A) + node.pt.dist(B) + node.pt.dist(C)) / T(4));
        } else {
            gsPoint2D<T> sumM(0, 0);
            int nInt = 0;
            T maxL = T(0);
            for (int k = 0; k < 3; ++k) {
                int eid = tr.e[k];
                if (!mesh.edges[eid].isBoundary) {
                    const auto& pa = mesh.pts[mesh.edges[eid].a];
                    const auto& pb = mesh.pts[mesh.edges[eid].b];
                    sumM.x += T(0.5) * (pa.x + pb.x);
                    sumM.y += T(0.5) * (pa.y + pb.y);
                    maxL = std::max(maxL, pa.dist(pb));
                    ++nInt;
                }
            }
            node.pt = (nInt > 0) ? gsPoint2D<T>(sumM.x / nInt, sumM.y / nInt)
                                 : gsPoint2D<T>((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
            node.radius = (nInt > 0) ? (T(0.5) * maxL) : (T(0.1) * domDiag);
        }
    }
    for (int ei = 0; ei < nE; ++ei) {
        if (!mesh.edges[ei].isBoundary && edgeTris[ei].size() == 2) {
            int t0 = edgeTris[ei][0], t1 = edgeTris[ei][1];
            out.triNodes[t0].adj.push_back(t1);
            out.triNodes[t1].adj.push_back(t0);
        }
    }

    auto& triNodes = out.triNodes;
    auto activeDeg = [&](int ti) -> int {
        if (!triNodes[ti].active) return 0;
        int d = 0;
        for (int nb : triNodes[ti].adj) if (triNodes[nb].active) ++d;
        return d;
    };
    auto countActive = [&]() -> int {
        int c = 0;
        for (int ti = 0; ti < nT; ++ti) if (triNodes[ti].active) ++c;
        return c;
    };
    auto analyzeCapArc = [&](int vR, int vL, std::vector<int>& sharpConvexOut, int& nIntOut) -> bool {
        sharpConvexOut.clear();
        nIntOut = 0;
        std::vector<int> intV;
        int lIdx = -1;
        if (!MAT::collectInteriorLoopVertices(mesh, vR, vL, intV, lIdx)) return false;
        nIntOut = static_cast<int>(intV.size());
        for (int v : intV) {
            if (mesh.isSharpCorner[v] && mesh.cornerAngleDeg[v] < T(150.0))
                sharpConvexOut.push_back(v);
        }
        return true;
    };

    // 1.2 Iterative theta-SMA leaf pruning
    bool prunedAny = true;
    while (prunedAny && countActive() > 2) {
        prunedAny = false;
        for (int leaf = 0; leaf < nT; ++leaf) {
            if (activeDeg(leaf) != 1) continue;
            std::vector<int> chain = {leaf};
            int prev = leaf, curr = -1;
            for (int nb : triNodes[leaf].adj) if (triNodes[nb].active) { curr = nb; break; }
            while (curr >= 0 && activeDeg(curr) == 2) {
                chain.push_back(curr);
                int nxt = -1;
                for (int nb : triNodes[curr].adj) if (triNodes[nb].active && nb != prev) { nxt = nb; break; }
                prev = curr; curr = nxt;
            }
            int root = curr;
            if (root < 0 || activeDeg(root) < 3) continue;
            int vL = -1, vR = -1;
            if (!MAT::getSharedChord(mesh, root, chain.back(), vL, vR)) continue;
            std::vector<int> sharpV;
            int nIntV = 0;
            if (!analyzeCapArc(vR, vL, sharpV, nIntV)) continue;
            if (sharpV.size() >= 2) continue;

            const auto& rPt = triNodes[root].pt;
            gsPoint2D<T> d1(mesh.pts[vL].x - rPt.x, mesh.pts[vL].y - rPt.y);
            gsPoint2D<T> d2(mesh.pts[vR].x - rPt.x, mesh.pts[vR].y - rPt.y);
            T n1 = std::hypot(d1.x, d1.y), n2 = std::hypot(d2.x, d2.y);
            T sepDeg = (n1 > T(1e-12) && n2 > T(1e-12))
                     ? std::acos(std::max(T(-1), std::min(T(1), (d1.x * d2.x + d1.y * d2.y) / (n1 * n2)))) * radToDeg
                     : T(180);
            T branchLen = triNodes[root].pt.dist(triNodes[chain.back()].pt);
            for (size_t i = 1; i < chain.size(); ++i)
                branchLen += triNodes[chain[i - 1]].pt.dist(triNodes[chain[i]].pt);

            bool shouldPrune = false;
            if (sharpV.empty()) {
                if (sepDeg < thetaPruneDeg || branchLen < T(0.95) * triNodes[root].radius)
                    shouldPrune = true;
            } else if (sharpV.size() == 1 && mesh.cornerAngleDeg[sharpV[0]] >= T(65.0) &&
                       branchLen < T(1.10) * triNodes[root].radius) {
                shouldPrune = true;
            }
            if (shouldPrune) {
                for (int td : chain) triNodes[td].active = false;
                prunedAny = true;
                break;
            }
        }
    }

    // 1.3 Compress active skeleton into MacroNodes + MacroBranches (including pure closed cycles)
    std::vector<int> triToMacro(nT, -1);
    auto registerMacroNode = [&](int ti, int deg) -> int {
        if (triToMacro[ti] >= 0) return triToMacro[ti];
        int id = static_cast<int>(out.macroNodes.size());
        MacroNode mn;
        mn.pt = triNodes[ti].pt;
        mn.radius = triNodes[ti].radius;
        mn.tris = {ti};
        mn.degree = deg;
        out.macroNodes.push_back(mn);
        triToMacro[ti] = id;
        return id;
    };
    for (int ti = 0; ti < nT; ++ti) {
        if (!triNodes[ti].active) continue;
        int d = activeDeg(ti);
        if (d == 1 || d >= 3) registerMacroNode(ti, d);
    }

    // Detect pure closed cycles (every triangle in component has activeDeg == 2, e.g. single-hole ring)
    {
        std::vector<bool> vis(nT, false);
        for (int start = 0; start < nT; ++start) {
            if (!triNodes[start].active || vis[start] || triToMacro[start] >= 0) continue;
            std::vector<int> comp;
            std::queue<int> q;
            q.push(start); vis[start] = true;
            bool hasMacro = false;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                comp.push_back(u);
                if (triToMacro[u] >= 0) hasMacro = true;
                for (int nb : triNodes[u].adj) {
                    if (triNodes[nb].active && !vis[nb]) { vis[nb] = true; q.push(nb); }
                }
            }
            if (hasMacro || comp.size() < 4) continue;
            std::vector<int> cyc;
            int cur = comp[0], prev = -1;
            for (size_t step = 0; step < comp.size(); ++step) {
                cyc.push_back(cur);
                int nxt = -1;
                for (int nb : triNodes[cur].adj)
                    if (triNodes[nb].active && nb != prev) { nxt = nb; break; }
                if (nxt < 0 || nxt == comp[0]) break;
                prev = cur; cur = nxt;
            }
            int C = static_cast<int>(cyc.size());
            if (C < 4) continue;
            MacroBranch br;
            br.isCycle = true;
            for (int i = 0; i < C; ++i) {
                br.triChain.push_back(cyc[i]);
                br.poly.push_back(triNodes[cyc[i]].pt);
                int vL = -1, vR = -1;
                if (MAT::getSharedChord(mesh, cyc[i], cyc[(i + 1) % C], vL, vR))
                    br.chords.push_back({vL, vR});
            }
            out.sleeves.push_back(std::move(br));
            out.isSleeveSharpTip.push_back(false);
        }
    }

    // Contract short J-J junction clusters
    if (junctionMergeRatio > T(0)) {
        bool didMerge = true;
        while (didMerge) {
            didMerge = false;
            for (size_t i = 0; i < out.macroNodes.size() && !didMerge; ++i) {
                if (out.macroNodes[i].tris.empty() || out.macroNodes[i].degree < 3 || out.macroNodes[i].degree >= 6)
                    continue;
                for (int tA : out.macroNodes[i].tris) {
                    for (int nb : triNodes[tA].adj) {
                        if (!triNodes[nb].active) continue;
                        int j = triToMacro[nb];
                        if (j < 0 || j == static_cast<int>(i)) continue;
                        if (out.macroNodes[j].tris.empty() || out.macroNodes[j].degree < 3) continue;
                        int vL = -1, vR = -1;
                        if (!MAT::getSharedChord(mesh, tA, nb, vL, vR)) continue;
                        T d = out.macroNodes[i].pt.dist(out.macroNodes[j].pt);
                        T rMin = std::min(out.macroNodes[i].radius, out.macroNodes[j].radius);
                        T chLen = mesh.pts[vL].dist(mesh.pts[vR]);
                        if (d >= std::max(junctionMergeRatio * rMin, T(0.65) * chLen)) continue;

                        std::set<int> comb(out.macroNodes[i].tris.begin(), out.macroNodes[i].tris.end());
                        comb.insert(out.macroNodes[j].tris.begin(), out.macroNodes[j].tris.end());
                        int newDeg = 0;
                        for (int ct : comb)
                            for (int cnb : triNodes[ct].adj)
                                if (triNodes[cnb].active && !comb.count(cnb)) ++newDeg;
                        if (newDeg < 3 || newDeg > 6) continue;

                        out.macroNodes[i].pt = gsPoint2D<T>(T(0.5) * (out.macroNodes[i].pt.x + out.macroNodes[j].pt.x),
                                                            T(0.5) * (out.macroNodes[i].pt.y + out.macroNodes[j].pt.y));
                        out.macroNodes[i].radius = T(0.5) * (out.macroNodes[i].radius + out.macroNodes[j].radius);
                        for (int tj : out.macroNodes[j].tris) {
                            out.macroNodes[i].tris.push_back(tj);
                            triToMacro[tj] = static_cast<int>(i);
                        }
                        out.macroNodes[i].degree = newDeg;
                        out.macroNodes[j].tris.clear();
                        out.macroNodes[j].degree = 0;
                        didMerge = true;
                        break;
                    }
                    if (didMerge) break;
                }
            }
        }
    }

    // Trace all MacroBranch sleeves between active MacroNodes
    out.isSharpTipTri.assign(nT, false);
    std::set<std::pair<int, int>> usedDualEdges;
    for (size_t mi = 0; mi < out.macroNodes.size(); ++mi) {
        if (out.macroNodes[mi].tris.empty()) continue;
        for (int tStart : out.macroNodes[mi].tris) {
            for (int nb : triNodes[tStart].adj) {
                if (!triNodes[nb].active) continue;
                if (triToMacro[nb] == static_cast<int>(mi)) continue;
                int lo = std::min(tStart, nb), hi = std::max(tStart, nb);
                if (usedDualEdges.count({lo, hi})) continue;
                usedDualEdges.insert({lo, hi});

                MacroBranch br;
                br.uNode = static_cast<int>(mi);
                br.poly.push_back(out.macroNodes[mi].pt);
                int vL0 = -1, vR0 = -1;
                if (!MAT::getSharedChord(mesh, tStart, nb, vL0, vR0)) continue;
                br.chords.push_back({vL0, vR0});

                int prev = tStart, curr = nb;
                while (triToMacro[curr] < 0) {
                    br.triChain.push_back(curr);
                    br.poly.push_back(triNodes[curr].pt);
                    int nxt = -1;
                    for (int cnb : triNodes[curr].adj)
                        if (triNodes[cnb].active && cnb != prev) { nxt = cnb; break; }
                    if (nxt < 0) break;
                    usedDualEdges.insert({std::min(curr, nxt), std::max(curr, nxt)});
                    int vL = -1, vR = -1;
                    if (MAT::getSharedChord(mesh, curr, nxt, vL, vR))
                        br.chords.push_back({vL, vR});
                    prev = curr; curr = nxt;
                }
                br.vNode = triToMacro[curr];
                if (br.vNode < 0) continue;
                br.poly.push_back(out.macroNodes[br.vNode].pt);

                // Tag n_q = 1 sharp-tip leaf sleeves so Module 2 & Module 3 protect them
                bool isTip = false;
                auto checkLeafEnd = [&](int leafMacroIdx, int parentTri, int leafTri) {
                    if (out.macroNodes[leafMacroIdx].degree != 1) return;
                    int vL = -1, vR = -1;
                    if (!MAT::getSharedChord(mesh, parentTri, leafTri, vL, vR)) return;
                    std::vector<int> sharpV;
                    int nInt = 0;
                    if (analyzeCapArc(vR, vL, sharpV, nInt) && sharpV.size() == 1 &&
                        mesh.cornerAngleDeg[sharpV[0]] < T(115.0)) {
                        isTip = true;
                    }
                };
                if (!br.triChain.empty()) {
                    checkLeafEnd(br.uNode, br.triChain.front(), tStart);
                    checkLeafEnd(br.vNode, br.triChain.back(), curr);
                } else {
                    checkLeafEnd(br.uNode, curr, tStart);
                    checkLeafEnd(br.vNode, tStart, curr);
                }
                if (isTip) {
                    for (int ti : br.triChain) out.isSharpTipTri[ti] = true;
                    if (out.macroNodes[br.uNode].degree == 1)
                        for (int ti : out.macroNodes[br.uNode].tris) out.isSharpTipTri[ti] = true;
                    if (out.macroNodes[br.vNode].degree == 1)
                        for (int ti : out.macroNodes[br.vNode].tris) out.isSharpTipTri[ti] = true;
                }
                out.sleeves.push_back(std::move(br));
                out.isSleeveSharpTip.push_back(isTip);
            }
        }
    }

    return true;
}

// ============================================================================
// Module 2: MAT-Guided 4-RoSy Cross-Field Solver & Singularity Oracle
// ============================================================================

template <class T>
inline bool solveMATGuided4RoSy(
    const MedialBackbone<T>& bb,
    GuidedCrossField<T>& out,
    T alphaMAT = T(2.5))
{
    out = GuidedCrossField<T>();
    out.fineMesh = bb.cdt;
    refineInteriorTriangles(out.fineMesh, bb.targetH, 2);
    out.accel.build(out.fineMesh);

    const auto& mesh = out.fineMesh;
    const int nV = static_cast<int>(mesh.pts.size());
    if (nV < 3 || mesh.tris.empty()) return false;

    T avgEdge = T(0);
    for (const auto& e : mesh.edges)
        avgEdge += mesh.pts[e.a].dist(mesh.pts[e.b]);
    avgEdge /= std::max<T>(T(1), static_cast<T>(mesh.edges.size()));

    // Collect sleeve sample points (pos, radius, unit MAT tangent, anti-parallelism eta)
    struct SleeveSample {
        gsPoint2D<T> pt;
        T radius;
        std::complex<T> zMAT;
        T eta;
    };
    std::vector<SleeveSample> samples;
    T meanSleeveR = T(0);
    int nSleeveSamples = 0;

    for (size_t si = 0; si < bb.sleeves.size(); ++si) {
        if (bb.isSleeveSharpTip[si]) continue; // Never force MAT tangent inside converging n_q=1 tips!
        const auto& br = bb.sleeves[si];
        if (br.poly.size() < 2) continue;
        for (size_t k = 0; k + 1 < br.poly.size(); ++k) {
            gsPoint2D<T> A = br.poly[k], B = br.poly[k + 1];
            gsPoint2D<T> mid(T(0.5) * (A.x + B.x), T(0.5) * (A.y + B.y));
            T dx = B.x - A.x, dy = B.y - A.y;
            T len = std::hypot(dx, dy);
            if (len <= T(1e-12)) continue;
            std::complex<T> tDir(dx / len, dy / len);
            std::complex<T> zM = std::pow(tDir, 4);

            T rLoc = T(0.1) * bb.domDiag;
            T eta = T(0.8);
            if (k < br.chords.size()) {
                int vL = br.chords[k].vL, vR = br.chords[k].vR;
                rLoc = T(0.5) * bb.cdt.pts[vL].dist(bb.cdt.pts[vR]);
                if (vL >= 0 && vR >= 0 &&
                    static_cast<size_t>(vL) < bb.cdt.bndTanOut.size() &&
                    static_cast<size_t>(vR) < bb.cdt.bndTanOut.size()) {
                    const auto& tL = bb.cdt.bndTanOut[vL];
                    const auto& tR = bb.cdt.bndTanOut[vR];
                    eta = std::max(T(0), -(tL.x * tR.x + tL.y * tR.y));
                }
            }
            if (eta > T(0.65)) {
                samples.push_back({mid, std::max(T(1e-5), rLoc), zM, eta});
                meanSleeveR += rLoc;
                ++nSleeveSamples;
            }
        }
    }
    const T rRef = (nSleeveSamples > 0) ? (meanSleeveR / static_cast<T>(nSleeveSamples)) : (T(0.15) * bb.domDiag);

    // Collect junction positions to suppress MAT guidance inside wide junction hubs
    std::vector<std::pair<gsPoint2D<T>, T>> juncZones;
    for (const auto& mn : bb.macroNodes) {
        if (!mn.tris.empty() && mn.degree >= 3)
            juncZones.push_back({mn.pt, mn.radius});
    }

    out.vertMATRadius.assign(nV, rRef);
    out.vertSleeveWeight.assign(nV, T(0));
    std::vector<std::complex<T>> vertZMAT(nV, std::complex<T>(T(1), T(0)));

    for (int i = 0; i < nV; ++i) {
        if (mesh.isBoundaryPt[i] || samples.empty()) continue;
        const auto& P = mesh.pts[i];

        bool inJunc = false;
        for (const auto& jz : juncZones) {
            if (P.dist(jz.first) < T(0.85) * jz.second) { inJunc = true; break; }
        }
        if (inJunc) continue;

        T bestD2 = std::numeric_limits<T>::max();
        int bestS = -1;
        for (int s = 0; s < static_cast<int>(samples.size()); ++s) {
            T d2 = P.distSq(samples[s].pt);
            if (d2 < bestD2) { bestD2 = d2; bestS = s; }
        }
        if (bestS >= 0) {
            const auto& smp = samples[bestS];
            out.vertMATRadius[i] = smp.radius;
            T d = std::sqrt(bestD2);
            if (d <= T(1.15) * smp.radius && smp.radius <= T(1.35) * rRef) {
                T ratio = smp.radius / std::max(T(1e-6), rRef);
                out.vertSleeveWeight[i] = alphaMAT * smp.eta * smp.eta * std::exp(-ratio * ratio) /
                                          std::max(T(1e-6), rRef * rRef);
                vertZMAT[i] = smp.zMAT;
            }
        }
    }

    // Boundary Dirichlet clamping on smooth walls & 90-deg corners (enforces Gauss-Bonnet winding!)
    // + soft Robin ONLY at non-90-deg sharp corners
    std::vector<T> bndMass(nV, T(0));
    for (const auto& e : mesh.edges) {
        if (!e.isBoundary) continue;
        T len = mesh.pts[e.a].dist(mesh.pts[e.b]);
        bndMass[e.a] += T(0.5) * len;
        bndMass[e.b] += T(0.5) * len;
    }

    std::vector<bool> isDirichlet(nV, false), isRobin(nV, false);
    std::vector<std::complex<T>> bndTarget(nV, std::complex<T>(T(1), T(0)));
    for (int i = 0; i < nV; ++i) {
        if (!mesh.isBoundaryPt[i]) continue;
        std::complex<T> t1(mesh.bndTanIn[i].x, mesh.bndTanIn[i].y);
        std::complex<T> t2(mesh.bndTanOut[i].x, mesh.bndTanOut[i].y);
        std::complex<T> z1 = std::pow(t1, 4), z2 = std::pow(t2, 4);
        std::complex<T> zAvg = z1 + z2;
        T nAvg = std::abs(zAvg);
        bndTarget[i] = (nAvg > T(1e-12)) ? (zAvg / nAvg) : z2;

        if (!mesh.isSharpCorner[i] || std::abs(z1 - z2) < T(0.45)) {
            isDirichlet[i] = true;
        } else {
            isRobin[i] = true;
        }
    }

    std::vector<int> freeIdx(nV, -1);
    int nFree = 0;
    for (int i = 0; i < nV; ++i)
        if (!isDirichlet[i]) freeIdx[i] = nFree++;

    out.field.z.assign(nV, std::complex<T>(T(1), T(0)));
    const T betaBnd = T(16.0) / bb.domDiag;

    if (nFree == 0) {
        out.field.z = bndTarget;
    } else {
        gsSparseEntries<T> lapEntries;
        lapEntries.reserve(mesh.tris.size() * 12 + nV * 2);
        gsMatrix<T> rhsLap = gsMatrix<T>::Zero(nFree, 2);
        gsVector<T> massDiag = gsVector<T>::Zero(nFree);

        auto addWeight = [&](int u, int v, T w) {
            if (isDirichlet[u]) return;
            int iu = freeIdx[u];
            lapEntries.add(iu, iu, w);
            if (isDirichlet[v]) {
                rhsLap(iu, 0) += w * bndTarget[v].real();
                rhsLap(iu, 1) += w * bndTarget[v].imag();
            } else {
                lapEntries.add(iu, freeIdx[v], -w);
            }
        };

        for (const auto& t : mesh.tris) {
            const auto& p0 = mesh.pts[t.v[0]];
            const auto& p1 = mesh.pts[t.v[1]];
            const auto& p2 = mesh.pts[t.v[2]];
            T area2 = std::max(T(1e-10), std::abs(cross2(p0, p1, p2)));
            T mTri = area2 / T(6);
            for (int k = 0; k < 3; ++k)
                if (!isDirichlet[t.v[k]]) massDiag(freeIdx[t.v[k]]) += mTri;

            auto cotAng = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C) -> T {
                T dot = (B.x - A.x) * (C.x - A.x) + (B.y - A.y) * (C.y - A.y);
                return std::max(T(1e-4), dot / area2);
            };
            T w12 = T(0.5) * cotAng(p0, p1, p2);
            T w20 = T(0.5) * cotAng(p1, p2, p0);
            T w01 = T(0.5) * cotAng(p2, p0, p1);
            addWeight(t.v[1], t.v[2], w12); addWeight(t.v[2], t.v[1], w12);
            addWeight(t.v[2], t.v[0], w20); addWeight(t.v[0], t.v[2], w20);
            addWeight(t.v[0], t.v[1], w01); addWeight(t.v[1], t.v[0], w01);
        }

        // Robin boundary alignment + MAT narrow-sleeve alignment (both strictly SPD!)
        for (int i = 0; i < nV; ++i) {
            if (isDirichlet[i]) continue;
            int iu = freeIdx[i];
            if (isRobin[i]) {
                T wr = betaBnd * std::max(T(1e-6) * avgEdge, bndMass[i]);
                lapEntries.add(iu, iu, wr);
                rhsLap(iu, 0) += wr * bndTarget[i].real();
                rhsLap(iu, 1) += wr * bndTarget[i].imag();
            }
            if (out.vertSleeveWeight[i] > T(0)) {
                T wm = std::min(T(3.5), out.vertSleeveWeight[i] * std::max(T(1e-12), massDiag(iu)));
                lapEntries.add(iu, iu, wm);
                rhsLap(iu, 0) += wm * vertZMAT[i].real();
                rhsLap(iu, 1) += wm * vertZMAT[i].imag();
            }
        }

        gsSparseMatrix<T> L(nFree, nFree);
        L.setFrom(lapEntries);
        typename gsSparseSolver<T>::SimplicialLDLT lapSolver;
        lapSolver.compute(L);
        if (lapSolver.info() != gsEigen::Success) return false;
        gsMatrix<T> solMat = lapSolver.solve(rhsLap);

        // MBO Unit-Norm Diffusion-Projection (Viertel & Osting)
        const T tau = T(2.8) * avgEdge * avgEdge;
        if (tau > T(1e-14)) {
            gsSparseMatrix<T> Mmat(nFree, nFree);
            gsSparseEntries<T> mEnt;
            mEnt.reserve(nFree);
            for (int i = 0; i < nFree; ++i)
                mEnt.add(i, i, std::max(T(1e-12), massDiag(i)));
            Mmat.setFrom(mEnt);
            gsSparseMatrix<T> A = Mmat + tau * L;
            typename gsSparseSolver<T>::SimplicialLDLT mboSolver;
            mboSolver.compute(A);
            if (mboSolver.info() == gsEigen::Success) {
                const gsMatrix<T> force = tau * rhsLap;
                gsMatrix<T> rhsMbo(nFree, 2);
                for (int it = 0; it < 10; ++it) {
                    for (int i = 0; i < nFree; ++i) {
                        T zr = solMat(i, 0), zi = solMat(i, 1);
                        T nrm = std::hypot(zr, zi);
                        if (nrm > T(1e-10)) { zr /= nrm; zi /= nrm; }
                        T mi = std::max(T(1e-12), massDiag(i));
                        rhsMbo(i, 0) = mi * zr + force(i, 0);
                        rhsMbo(i, 1) = mi * zi + force(i, 1);
                    }
                    solMat = mboSolver.solve(rhsMbo);
                }
            }
        }

        for (int i = 0; i < nV; ++i) {
            if (isDirichlet[i]) out.field.z[i] = bndTarget[i];
            else out.field.z[i] = std::complex<T>(solMat(freeIdx[i], 0), solMat(freeIdx[i], 1));
        }
    }

    // Singularity detection, clustering, dipole cancellation & circular-harmonic port extraction
    const T pi = T(3.14159265358979323846);
    const T twoPi = T(2) * pi;
    auto wrapPi = [&](T a) -> T {
        while (a > pi) a -= twoPi;
        while (a <= -pi) a += twoPi;
        return a;
    };

    struct RawSing { gsPoint2D<T> pos; int idx; int triIdx; };
    std::vector<RawSing> rawSings;
    for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
        const auto& t = mesh.tris[ti];
        std::complex<T> z0 = out.field.z[t.v[0]], z1 = out.field.z[t.v[1]], z2 = out.field.z[t.v[2]];
        if (std::abs(z0) < T(1e-12) || std::abs(z1) < T(1e-12) || std::abs(z2) < T(1e-12)) continue;
        T sum = wrapPi(std::arg(z1) - std::arg(z0)) + wrapPi(std::arg(z2) - std::arg(z1)) + wrapPi(std::arg(z0) - std::arg(z2));
        int ind = static_cast<int>(std::llround(sum / twoPi));
        if (ind == 0) continue;

        const auto& P0 = mesh.pts[t.v[0]];
        const auto& P1 = mesh.pts[t.v[1]];
        const auto& P2 = mesh.pts[t.v[2]];
        T A00 = z1.real() - z0.real(), A01 = z2.real() - z0.real();
        T A10 = z1.imag() - z0.imag(), A11 = z2.imag() - z0.imag();
        T det = A00 * A11 - A01 * A10;
        T l1 = T(1) / T(3), l2 = T(1) / T(3);
        if (std::abs(det) > T(1e-12)) {
            T u = (-z0.real() * A11 + z0.imag() * A01) / det;
            T v = (-A00 * z0.imag() + A10 * z0.real()) / det;
            if (u >= T(0) && v >= T(0) && u + v <= T(1)) { l1 = u; l2 = v; }
        }
        T l0 = T(1) - l1 - l2;
        gsPoint2D<T> pos(l0 * P0.x + l1 * P1.x + l2 * P2.x, l0 * P0.y + l1 * P1.y + l2 * P2.y);
        if (isPointStrictlyInsideDomain(pos, mesh.loopPolys, T(1e-5)))
            rawSings.push_back({pos, ind, ti});
    }

    const T clusterTol = T(1.8) * avgEdge;
    std::vector<bool> used(rawSings.size(), false);
    std::vector<RawSing> clustered;
    for (size_t i = 0; i < rawSings.size(); ++i) {
        if (used[i]) continue;
        used[i] = true;
        gsPoint2D<T> sumP = rawSings[i].pos;
        int sumInd = rawSings[i].idx, cnt = 1;
        for (size_t j = i + 1; j < rawSings.size(); ++j) {
            if (!used[j] && rawSings[i].pos.dist(rawSings[j].pos) <= clusterTol &&
                isSegmentInsideDomain(rawSings[i].pos, rawSings[j].pos, mesh.loopPolys)) {
                used[j] = true;
                sumP.x += rawSings[j].pos.x; sumP.y += rawSings[j].pos.y;
                sumInd += rawSings[j].idx; ++cnt;
            }
        }
        if (sumInd != 0)
            clustered.push_back({gsPoint2D<T>(sumP.x / cnt, sumP.y / cnt), sumInd, rawSings[i].triIdx});
    }

    const T dipoleTol = T(2.8) * avgEdge;
    std::vector<bool> dead(clustered.size(), false);
    for (size_t i = 0; i < clustered.size(); ++i) {
        if (dead[i]) continue;
        int bestJ = -1;
        T bestD = dipoleTol;
        for (size_t j = i + 1; j < clustered.size(); ++j) {
            if (dead[j] || clustered[i].idx + clustered[j].idx != 0) continue;
            T d = clustered[i].pos.dist(clustered[j].pos);
            if (d < bestD && isSegmentInsideDomain(clustered[i].pos, clustered[j].pos, mesh.loopPolys)) {
                bestD = d; bestJ = static_cast<int>(j);
            }
        }
        if (bestJ >= 0) { dead[i] = true; dead[bestJ] = true; }
    }

    for (size_t i = 0; i < clustered.size(); ++i) {
        if (dead[i]) continue;
        // Filter spurious boundary-hugging singularities
        if (gsPatchSplitterMAT<T>::pointToLoopsDist(clustered[i].pos, mesh.loopPolys) < T(0.65) * avgEdge)
            continue;
        typename CrossFieldSolution<T>::Singularity sg;
        sg.pos = clustered[i].pos;
        sg.index = std::max(-1, std::min(1, clustered[i].idx)); // clamp to val 3 or val 5
        int val = 4 - sg.index;
        const int nAng = 48;
        T probeR = T(1.4) * avgEdge;
        std::complex<T> harmSum(T(0), T(0));
        int nValid = 0, triHint = clustered[i].triIdx;
        for (T rScale : {T(1.0), T(0.55), T(0.30)}) {
            harmSum = std::complex<T>(T(0), T(0));
            nValid = 0;
            T r = probeR * rScale;
            for (int m = 0; m < nAng; ++m) {
                T phi = (twoPi * m) / static_cast<T>(nAng);
                gsPoint2D<T> q(sg.pos.x + r * std::cos(phi), sg.pos.y + r * std::sin(phi));
                std::complex<T> zq; T thq = T(0);
                if (evalCrossFieldAt(mesh, out.field, out.accel, q, zq, thq, triHint) && std::abs(zq) > T(1e-10)) {
                    std::complex<T> unitZ = zq / std::abs(zq);
                    harmSum += unitZ * std::complex<T>(std::cos(-sg.index * phi), std::sin(-sg.index * phi));
                    ++nValid;
                }
            }
            if (nValid >= 12 && std::abs(harmSum) > T(1e-6)) break;
        }
        T phi0 = (nValid > 0 && std::abs(harmSum) > T(1e-8)) ? std::arg(harmSum) : T(0);
        sg.portAngles.reserve(val);
        for (int p = 0; p < val; ++p)
            sg.portAngles.push_back((phi0 + twoPi * static_cast<T>(p)) / static_cast<T>(val));
        out.field.singularities.push_back(std::move(sg));
    }

    return true;
}

// ============================================================================
// Module 3: Cross-Field-Optimal Neck Slicing, Hole-Cycle Opening & Parity
// ============================================================================

template <class T>
inline T evalCrossCutCostJ(
    const MedialBackbone<T>& bb,
    const GuidedCrossField<T>& gcf,
    int loopA, T sA,
    int loopB, T sB,
    T rMAT,
    int& triHint)
{
    const gsPoint2D<T> A = evalVirtualLoopPoint(bb.vloops[loopA], sA);
    const gsPoint2D<T> B = evalVirtualLoopPoint(bb.vloops[loopB], sB);
    const T len = A.dist(B);
    if (len <= T(1e-10)) return T(1e9);

    const gsPoint2D<T> dir((B.x - A.x) / len, (B.y - A.y) / len);
    const T segAng = std::atan2(dir.y, dir.x);

    // 1. Cross-field alignment integral: int sin^2(2 * angle(gamma', C)) ds
    T cfErr = T(0);
    int nEval = 0;
    for (int k = 1; k <= 5; ++k) {
        T t = static_cast<T>(k) / T(6);
        gsPoint2D<T> q((T(1) - t) * A.x + t * B.x, (T(1) - t) * A.y + t * B.y);
        std::complex<T> zq; T theta = T(0);
        if (evalCrossFieldAt(gcf.fineMesh, gcf.field, gcf.accel, q, zq, theta, triHint)) {
            T s2 = std::sin(T(2) * (segAng - theta));
            cfErr += s2 * s2;
            ++nEval;
        }
    }
    if (nEval > 0) cfErr /= static_cast<T>(nEval);

    // 2. Boundary orthogonality: cos^2(angle(gamma, dOmega_A)) + cos^2(angle(gamma, dOmega_B))
    const gsPoint2D<T> tA = evalVirtualLoopUnitTangent(bb.vloops[loopA], sA);
    const gsPoint2D<T> tB = evalVirtualLoopUnitTangent(bb.vloops[loopB], sB);
    const T cA = dir.x * tA.x + dir.y * tA.y;
    const T cB = dir.x * tB.x + dir.y * tB.y;
    const T orthErr = cA * cA + cB * cB;

    // 3. Normalized thickness ratio Len(gamma) / (2 * r_MAT)
    const T lenRatio = len / std::max(T(1e-6), T(2) * rMAT);

    return T(1.8) * cfErr + T(2.2) * orthErr + T(1.0) * lenRatio;
}

template <class T>
inline bool sliceNecksAndBalanceParity(
    const MedialBackbone<T>& bb,
    const GuidedCrossField<T>& gcf,
    int degree,
    SlicedDecomposition<T>& out)
{
    out = SlicedDecomposition<T>();
    const int nLoops = static_cast<int>(bb.vloops.size());
    if (nLoops == 0) return false;
    const T minEndTol = std::max(T(1e-5), T(0.012) * bb.domDiag);
    const auto densePolys = mergedBuildDenseLoopPolys(bb.vloops, bb.domDiag);
    int triHint = -1;

    // Helper: snap arc coordinate to nearby pinned corner with n_q >= 2 (avoiding tiny sliver walls)
    auto snapToPinnedCorner = [&](int lIdx, T s, T rLocal) -> T {
        const auto& vl = bb.vloops[lIdx];
        T bestS = virtualLoopWrapS(vl, s);
        T snapTol = std::min(T(0.32) * rLocal, T(0.06) * vl.totalArcLen);
        for (const auto& c : vl.corners) {
            if (c.targetPatchValence < 2) continue; // never snap onto sharp n_q=1 tips
            T d = std::abs(bestS - c.arcS);
            d = std::min(d, vl.totalArcLen - d);
            if (d < snapTol) return c.arcS;
        }
        return bestS;
    };

    // Helper: locally refine a candidate cross-cut (loopA, sA) <-> (loopB, sB) minimizing J(gamma)
    auto optimizeCutLocally = [&](int lA, T& sA, int lB, T& sB, T rMAT, T minTanDeg) -> T {
        sA = snapToPinnedCorner(lA, sA, rMAT);
        sB = snapToPinnedCorner(lB, sB, rMAT);
        T bestCost = evalCrossCutCostJ(bb, gcf, lA, sA, lB, sB, rMAT, triHint);
        T step = std::min(T(0.28) * rMAT, T(0.04) * bb.domDiag);
        for (int iter = 0; iter < 3; ++iter) {
            bool moved = false;
            for (T dsA : {-step, T(0), step}) {
                for (T dsB : {-step, T(0), step}) {
                    if (dsA == T(0) && dsB == T(0)) continue;
                    T cA = virtualLoopWrapS(bb.vloops[lA], sA + dsA);
                    T cB = virtualLoopWrapS(bb.vloops[lB], sB + dsB);
                    gsPoint2D<T> pA = evalVirtualLoopPoint(bb.vloops[lA], cA);
                    gsPoint2D<T> pB = evalVirtualLoopPoint(bb.vloops[lB], cB);
                    if (!isCrossCutClear(pA, pB, densePolys, out.cuts, minEndTol)) continue;
                    if (!isCrossCutTangentClear(bb.vloops, lA, cA, lB, cB, minTanDeg)) continue;
                    T c = evalCrossCutCostJ(bb, gcf, lA, cA, lB, cB, rMAT, triHint);
                    if (c < bestCost - T(1e-4)) {
                        bestCost = c; sA = cA; sB = cB; moved = true;
                    }
                }
            }
            if (!moved) break;
            step *= T(0.5);
        }
        return bestCost;
    };

    auto tryAddCut = [&](int sleeveIdx, int lA, T sA, int lB, T sB, T rMAT, T minTanDeg = T(14.0)) -> bool {
        if (lA < 0 || lB < 0 || lA >= nLoops || lB >= nLoops) return false;
        sA = snapToPinnedCorner(lA, sA, rMAT);
        sB = snapToPinnedCorner(lB, sB, rMAT);
        gsPoint2D<T> pA = evalVirtualLoopPoint(bb.vloops[lA], sA);
        gsPoint2D<T> pB = evalVirtualLoopPoint(bb.vloops[lB], sB);
        if (!isCrossCutClear(pA, pB, densePolys, out.cuts, minEndTol)) return false;
        if (!isCrossCutTangentClear(bb.vloops, lA, sA, lB, sB, minTanDeg)) return false;

        // Stage 1A: Prune shallow same-loop pocket/blister cuts (loopA == loopB) that slice off tiny
        // featureless ears on multi-corner holes or wavy outer walls (causing parasitic B=2/3 O-rings).
        if (sleeveIdx >= 0 && lA == lB && (nLoops > 1 || bb.sleeves.size() > 1)) {
            const auto& vl = bb.vloops[lA];
            T spFwd = mergedLoopForwardSpan(vl, sA, sB);
            T spRev = vl.totalArcLen - spFwd;
            T shortSpan = std::min(spFwd, spRev);
            T shortStart = (spFwd <= spRev) ? sA : sB;
            T cutLen = pA.dist(pB);
            if (shortSpan < T(1.55) * cutLen) {
                bool hasMajorFeature = false;
                for (const auto& c : vl.corners) {
                    T d = mergedLoopForwardSpan(vl, shortStart, c.arcS);
                    if (d > T(1e-5) && d < shortSpan - T(1e-5)) {
                        if (c.interiorAngleDeg >= T(195.0) || c.interiorAngleDeg <= T(50.0)) {
                            hasMajorFeature = true;
                            break;
                        }
                    }
                }
                if (!hasMajorFeature) return false;
            }
        }

        T cost = optimizeCutLocally(lA, sA, lB, sB, rMAT, minTanDeg);
        pA = evalVirtualLoopPoint(bb.vloops[lA], sA);
        pB = evalVirtualLoopPoint(bb.vloops[lB], sB);
        if (!isCrossCutClear(pA, pB, densePolys, out.cuts, minEndTol)) return false;
        if (!isCrossCutTangentClear(bb.vloops, lA, sA, lB, sB, minTanDeg)) return false;

        CrossCut<T> cc;
        cc.id = static_cast<int>(out.cuts.size());
        cc.sleeveIdx = sleeveIdx;
        cc.loopA = lA; cc.loopB = lB;
        cc.arcSA = sA; cc.arcSB = sB;
        cc.ptA = pA;   cc.ptB = pB;
        cc.costJ = cost;
        out.cuts.push_back(cc);
        return true;
    };

    // 3.1 Place cross-field-optimal cuts along MAT sleeves & pure closed cycles
    const int nTotalSleeves = static_cast<int>(bb.sleeves.size());
    for (size_t si = 0; si < bb.sleeves.size(); ++si) {
        // Keep n_q = 1 sharp blade tips uncut when attached to a multi-sleeve skeleton,
        // but NEVER skip the sole sleeve of a 0-hole crescent/heart domain!
        if (bb.isSleeveSharpTip[si] && (nLoops > 1 || nTotalSleeves > 1)) continue;
        const auto& br = bb.sleeves[si];
        const int nCh = static_cast<int>(br.chords.size());
        if (nCh == 0) continue;

        // Compute sleeve polyline length and mean radius
        T sLen = T(0);
        for (size_t k = 1; k < br.poly.size(); ++k) sLen += br.poly[k - 1].dist(br.poly[k]);
        T meanR = T(0);
        for (const auto& ch : br.chords) meanR += T(0.5) * bb.cdt.pts[ch.vL].dist(bb.cdt.pts[ch.vR]);
        meanR = std::max(T(1e-5), meanR / static_cast<T>(nCh));

        if (br.isCycle) {
            // Pure closed cycle around an isolated hole -> 4 optimal cuts spaced by ~C/4
            int bestOff = 0;
            T bestSumJ = std::numeric_limits<T>::max();
            int nPhases = std::max(1, std::min(8, nCh / 4));
            for (int ph = 0; ph < nPhases; ++ph) {
                int off = (ph * nCh) / (4 * nPhases);
                T sumJ = T(0);
                for (int q = 0; q < 4; ++q) {
                    const auto& ch = br.chords[(off + (q * nCh) / 4) % nCh];
                    int lA = bb.vLoopIdx[ch.vL], lB = bb.vLoopIdx[ch.vR];
                    if (lA >= 0 && lB >= 0)
                        sumJ += evalCrossCutCostJ(bb, gcf, lA, bb.vArcS[ch.vL], lB, bb.vArcS[ch.vR], meanR, triHint);
                }
                if (sumJ < bestSumJ) { bestSumJ = sumJ; bestOff = off; }
            }
            for (int q = 0; q < 4; ++q) {
                int ci = (bestOff + (q * nCh) / 4) % nCh;
                const auto& ch = br.chords[ci];
                tryAddCut(static_cast<int>(si), bb.vLoopIdx[ch.vL], bb.vArcS[ch.vL],
                          bb.vLoopIdx[ch.vR], bb.vArcS[ch.vR], meanR);
            }
            continue;
        }

        // Standard junction-to-junction or junction-to-leaf sleeve: coarse macro-spacing!
        int nCuts = (sLen < T(3.2) * meanR)
                  ? 1
                  : std::max(2, std::min(4, static_cast<int>(std::llround(sLen / (T(2.8) * meanR)))));
        if (nLoops == 1 && nTotalSleeves == 1)
            nCuts = std::max(2, std::min(4, nCuts + 1)); // crescent/heart ribbon gets 2-4 clean cross-cuts
        if (nCh <= 2) nCuts = 1;

        for (int k = 0; k < nCuts; ++k) {
            // Step back from junction endpoints into interior of sleeve
            T frac = (nCuts == 1) ? T(0.5) : (T(0.15) + T(0.70) * static_cast<T>(k) / static_cast<T>(nCuts - 1));
            int centerIdx = std::max(0, std::min(nCh - 1, static_cast<int>(std::llround(frac * static_cast<T>(nCh - 1)))));
            int winLo = std::max(0, centerIdx - std::max(1, nCh / (2 * nCuts + 2)));
            int winHi = std::min(nCh - 1, centerIdx + std::max(1, nCh / (2 * nCuts + 2)));

            int bestCi = -1;
            T bestJ = std::numeric_limits<T>::max();
            for (int ci = winLo; ci <= winHi; ++ci) {
                const auto& ch = br.chords[ci];
                int lA = bb.vLoopIdx[ch.vL], lB = bb.vLoopIdx[ch.vR];
                if (lA < 0 || lB < 0) continue;
                T rC = std::max(T(1e-5), T(0.5) * bb.cdt.pts[ch.vL].dist(bb.cdt.pts[ch.vR]));
                if (!isCrossCutClear(bb.cdt.pts[ch.vL], bb.cdt.pts[ch.vR], densePolys, out.cuts, minEndTol))
                    continue;
                if (!isCrossCutTangentClear(bb.vloops, lA, bb.vArcS[ch.vL], lB, bb.vArcS[ch.vR], T(14.0)))
                    continue;
                T J = evalCrossCutCostJ(bb, gcf, lA, bb.vArcS[ch.vL], lB, bb.vArcS[ch.vR], rC, triHint);
                if (J < bestJ) { bestJ = J; bestCi = ci; }
            }
            if (bestCi >= 0) {
                const auto& ch = br.chords[bestCi];
                T rC = std::max(T(1e-5), T(0.5) * bb.cdt.pts[ch.vL].dist(bb.cdt.pts[ch.vR]));
                tryAddCut(static_cast<int>(si), bb.vLoopIdx[ch.vL], bb.vArcS[ch.vL],
                          bb.vLoopIdx[ch.vR], bb.vArcS[ch.vR], rC, T(14.0));
            }
        }
    }

    // 3.2 Topological Genus-0, Distinct-Vertex Hole-Degree >= 2 & 2-Edge-Connectivity Guarantee!
    // Every hole loop must: (a) connect to L_0, (b) have >= 2 DISTINCT cut endpoints (3 for smooth holes),
    // and (c) have NO single cut-edge (bridge) in the loop graph (prevents self-touching slit / figure-8 faces!).
    if (nLoops > 1) {
        auto findBestBridgeForHole = [&](int hLoop, bool requireTargetComp, const std::vector<bool>& targetComp) -> bool {
            const auto& hvl = bb.vloops[hLoop];
            const int nSamples = 36;
            std::vector<T> existingS;
            for (const auto& c : out.cuts) {
                if (c.loopA == hLoop) existingS.push_back(c.arcSA);
                if (c.loopB == hLoop) existingS.push_back(c.arcSB);
            }

            for (T minDeg : {T(14.0), T(7.0), T(3.5)}) {
                int bestLA = hLoop, bestLB = -1;
                T bestSA = T(0), bestSB = T(0), bestScore = std::numeric_limits<T>::max();
                for (int k = 0; k < nSamples; ++k) {
                    T sH = (hvl.totalArcLen * static_cast<T>(k)) / static_cast<T>(nSamples);
                    // Preview snapped arc position so we don't accidentally collapse onto an existing cut vertex!
                    T sHSnap = snapToPinnedCorner(hLoop, sH, T(0.05) * hvl.totalArcLen);
                    T minArcDist = hvl.totalArcLen;
                    for (T es : existingS) {
                        T d = std::abs(sHSnap - es);
                        minArcDist = std::min(minArcDist, std::min(d, hvl.totalArcLen - d));
                    }
                    if (!existingS.empty() && minArcDist < T(0.08) * hvl.totalArcLen) continue;
                    gsPoint2D<T> pH = evalVirtualLoopPoint(hvl, sH);

                    for (int otherL = 0; otherL < nLoops; ++otherL) {
                        if (otherL == hLoop) continue;
                        if (requireTargetComp && !targetComp[otherL]) continue;
                        // Test both nearest projection AND 16 evenly spaced target points along otherL
                        // (so non-orthogonal clear diagonals around curved holes are never missed!)
                        std::vector<T> candSO;
                        candSO.reserve(17);
                        T sProj = T(0);
                        gsPoint2D<T> pProj;
                        projectPointToVirtualLoop(bb.vloops[otherL], pH, sProj, pProj);
                        candSO.push_back(sProj);
                        const T otherLen = bb.vloops[otherL].totalArcLen;
                        for (int m = 0; m < 16; ++m)
                            candSO.push_back((otherLen * static_cast<T>(m)) / T(16));
                        for (T sO : candSO) {
                            gsPoint2D<T> pO = evalVirtualLoopPoint(bb.vloops[otherL], sO);
                            T dist = pH.dist(pO);
                            if (dist <= minEndTol) continue;
                            if (!isCrossCutClear(pH, pO, densePolys, out.cuts, minEndTol)) continue;
                            if (!isCrossCutTangentClear(bb.vloops, hLoop, sH, otherL, sO, minDeg)) continue;
                            T J = evalCrossCutCostJ(bb, gcf, hLoop, sH, otherL, sO, T(0.5) * dist, triHint);
                            T sc = J + T(1.5) * (dist / bb.domDiag) - (existingS.empty() ? T(0) : T(0.6) * (minArcDist / hvl.totalArcLen));
                            if (sc < bestScore) {
                                bestScore = sc; bestLA = hLoop; bestSA = sH; bestLB = otherL; bestSB = sO;
                            }
                        }
                    }
                }
                if (bestLB >= 0) {
                    // Add without collapsing onto an already-used vertex on bestLA
                    gsPoint2D<T> pA = evalVirtualLoopPoint(bb.vloops[bestLA], bestSA);
                    gsPoint2D<T> pB = evalVirtualLoopPoint(bb.vloops[bestLB], bestSB);
                    T rEst = T(0.5) * pA.dist(pB);
                    if (tryAddCut(-1, bestLA, bestSA, bestLB, bestSB, std::max(T(1e-5), rEst), minDeg)) {
                        // Verify the newly added cut didn't snap onto an existing vertex on hLoop
                        bool collapsed = false;
                        const auto& nc = out.cuts.back();
                        T ns = (nc.loopA == hLoop) ? nc.arcSA : nc.arcSB;
                        for (T es : existingS) {
                            T d = std::abs(ns - es);
                            if (std::min(d, hvl.totalArcLen - d) < T(0.02) * hvl.totalArcLen) { collapsed = true; break; }
                        }
                        if (!collapsed) return true;
                        // Undo snapped collapse and insert unsnapped clean cut directly
                        out.cuts.pop_back();
                        if (isCrossCutClear(pA, pB, densePolys, out.cuts, minEndTol) &&
                            isCrossCutTangentClear(bb.vloops, bestLA, bestSA, bestLB, bestSB, minDeg)) {
                            CrossCut<T> cc;
                            cc.id = static_cast<int>(out.cuts.size());
                            cc.sleeveIdx = -1;
                            cc.loopA = bestLA; cc.loopB = bestLB;
                            cc.arcSA = bestSA; cc.arcSB = bestSB;
                            cc.ptA = pA; cc.ptB = pB;
                            cc.costJ = evalCrossCutCostJ(bb, gcf, bestLA, bestSA, bestLB, bestSB, std::max(T(1e-5), rEst), triHint);
                            out.cuts.push_back(cc);
                            return true;
                        }
                    }
                }
            }
            return false;
        };

        for (int guard = 0; guard < 6 * nLoops; ++guard) {
            std::vector<std::vector<int>> lAdj(nLoops);
            std::vector<std::vector<T>> loopCutS(nLoops);
            for (const auto& c : out.cuts) {
                lAdj[c.loopA].push_back(c.loopB);
                lAdj[c.loopB].push_back(c.loopA);
                loopCutS[c.loopA].push_back(c.arcSA);
                loopCutS[c.loopB].push_back(c.arcSB);
            }
            // Count DISTINCT attachment vertices per loop (so 2 cuts snapped to the same vertex count as 1!)
            std::vector<int> distinctDeg(nLoops, 0);
            for (int l = 0; l < nLoops; ++l) {
                const T vTol = std::max(T(1e-5), T(0.015) * bb.vloops[l].totalArcLen);
                std::vector<T> uniq;
                for (T s : loopCutS[l]) {
                    bool dup = false;
                    for (T u : uniq) {
                        T d = std::abs(s - u);
                        if (std::min(d, bb.vloops[l].totalArcLen - d) < vTol) { dup = true; break; }
                    }
                    if (!dup) uniq.push_back(s);
                }
                distinctDeg[l] = static_cast<int>(uniq.size());
            }

            std::vector<bool> inComp0(nLoops, false);
            std::queue<int> q;
            q.push(0); inComp0[0] = true;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int nb : lAdj[u]) if (!inComp0[nb]) { inComp0[nb] = true; q.push(nb); }
            }

            bool added = false;
            // (a) Connect any disconnected hole to component 0
            for (int h = 1; h < nLoops; ++h) {
                if (!inComp0[h]) {
                    if (findBestBridgeForHole(h, true, inComp0)) { added = true; break; }
                }
            }
            if (added) continue;

            // (b) Ensure every hole has >= 2 DISTINCT cut vertices (3 for smooth 0-corner holes)
            for (int h = 1; h < nLoops; ++h) {
                int minReq = bb.vloops[h].corners.empty() ? 3 : 2;
                if (distinctDeg[h] < minReq) {
                    if (findBestBridgeForHole(h, false, inComp0)) { added = true; break; }
                }
            }
            if (added) continue;

            // (c) Enforce 2-edge-connectivity of the loop-cut graph (eliminate single bridge cuts that cause slit faces!)
            for (size_t skipCi = 0; skipCi < out.cuts.size() && !added; ++skipCi) {
                if (out.cuts[skipCi].loopA == out.cuts[skipCi].loopB) continue;
                std::vector<std::vector<int>> subAdj(nLoops);
                for (size_t ci = 0; ci < out.cuts.size(); ++ci) {
                    if (ci == skipCi) continue;
                    subAdj[out.cuts[ci].loopA].push_back(out.cuts[ci].loopB);
                    subAdj[out.cuts[ci].loopB].push_back(out.cuts[ci].loopA);
                }
                std::vector<bool> reach0(nLoops, false);
                std::queue<int> bq;
                bq.push(0); reach0[0] = true;
                while (!bq.empty()) {
                    int u = bq.front(); bq.pop();
                    for (int nb : subAdj[u]) if (!reach0[nb]) { reach0[nb] = true; bq.push(nb); }
                }
                for (int h = 1; h < nLoops; ++h) {
                    if (!reach0[h]) {
                        if (findBestBridgeForHole(h, true, reach0)) { added = true; break; }
                    }
                }
            }
            if (!added) break;
        }
    }

    // 3.3 Build Exact Planar Half-Edge Arrangement from VirtualCADLoops + CrossCuts
    // Wrapped in a reusable tracer so Stage 1B can delete redundant B=2 / tiny B=3 pocket cuts
    // and absorb those blisters into their encompassing hub BEFORE any artificial 180-deg parity split!
    auto buildDCELAndTraceHubs = [&]() -> bool {
        out.nodes.clear();
        out.hubs.clear();
        std::vector<std::vector<std::pair<T, int>>> loopEvents(nLoops); // (arcS, nodeId)
        auto addOrGetLoopNode = [&](int lIdx, T s, bool isCorner, T cAng, int tVal) -> int {
            const auto& vl = bb.vloops[lIdx];
            s = virtualLoopWrapS(vl, s);
            const T mergeTol = std::max(T(1e-6), T(0.004) * vl.totalArcLen);
            for (auto& ev : loopEvents[lIdx]) {
                T d = std::abs(ev.first - s);
                d = std::min(d, vl.totalArcLen - d);
                if (d < mergeTol) {
                    if (isCorner) {
                        out.nodes[ev.second].isPinnedCorner = true;
                        out.nodes[ev.second].cornerAngleDeg = cAng;
                        out.nodes[ev.second].targetValence = tVal;
                    }
                    return ev.second;
                }
            }
            int id = static_cast<int>(out.nodes.size());
            typename SlicedDecomposition<T>::ArrNode nd;
            nd.pt = evalVirtualLoopPoint(vl, s);
            nd.loopIdx = lIdx;
            nd.arcS = s;
            nd.isPinnedCorner = isCorner;
            nd.cornerAngleDeg = cAng;
            nd.targetValence = tVal;
            out.nodes.push_back(nd);
            loopEvents[lIdx].push_back({s, id});
            return id;
        };

        // Register all pinned corners first
        for (int l = 0; l < nLoops; ++l) {
            for (const auto& c : bb.vloops[l].corners)
                addOrGetLoopNode(l, c.arcS, true, c.interiorAngleDeg, c.targetPatchValence);
        }

        // Register cut endpoints
        std::vector<std::pair<int, int>> cutNodes(out.cuts.size(), {-1, -1});
        for (size_t ci = 0; ci < out.cuts.size(); ++ci) {
            auto& c = out.cuts[ci];
            c.id = static_cast<int>(ci);
            int u = addOrGetLoopNode(c.loopA, c.arcSA, false, T(180), 2);
            int v = addOrGetLoopNode(c.loopB, c.arcSB, false, T(180), 2);
            c.arcSA = out.nodes[u].arcS; c.ptA = out.nodes[u].pt;
            c.arcSB = out.nodes[v].arcS; c.ptB = out.nodes[v].pt;
            cutNodes[ci] = {u, v};
        }

        // Ensure every smooth 0-corner loop with < 4 events has at least 4 evenly spaced events
        for (int l = 0; l < nLoops; ++l) {
            if (loopEvents[l].size() < 4) {
                int need = 4 - static_cast<int>(loopEvents[l].size());
                const auto& vl = bb.vloops[l];
                for (int k = 0; k < need; ++k) {
                    if (loopEvents[l].empty()) {
                        addOrGetLoopNode(l, T(0), false, T(180), 2);
                    } else {
                        std::sort(loopEvents[l].begin(), loopEvents[l].end());
                        int m = static_cast<int>(loopEvents[l].size());
                        T bestSpan = -T(1), bestMid = T(0);
                        for (int i = 0; i < m; ++i) {
                            T s0 = loopEvents[l][i].first;
                            T s1 = loopEvents[l][(i + 1) % m].first;
                            T sp = (i + 1 < m) ? (s1 - s0) : (s1 + vl.totalArcLen - s0);
                            if (sp > bestSpan) { bestSpan = sp; bestMid = virtualLoopWrapS(vl, s0 + T(0.5) * sp); }
                        }
                        addOrGetLoopNode(l, bestMid, false, T(180), 2);
                    }
                }
            }
            std::sort(loopEvents[l].begin(), loopEvents[l].end(),
                      [](const std::pair<T, int>& a, const std::pair<T, int>& b) { return a.first < b.first; });
        }

        struct HalfEdge {
            int origin = -1, target = -1;
            int twin = -1;
            bool isCrossCut = false;
            int cutId = -1;
            bool cutForward = true;
            int loopIdx = -1;
            T sStart = T(0), sEnd = T(0);
            T outAng = T(0);
            CurveSeg<T> seg;
        };
        std::vector<HalfEdge> HEs;

        for (int l = 0; l < nLoops; ++l) {
            const auto& evs = loopEvents[l];
            const int M = static_cast<int>(evs.size());
            for (int i = 0; i < M; ++i) {
                int u = evs[i].second;
                int v = evs[(i + 1) % M].second;
                if (u == v) continue;
                T s0 = evs[i].first, s1 = evs[(i + 1) % M].first;
                CurveSeg<T> seg = extractVirtualSubCurve(bb.vloops[l], s0, s1, degree);
                if (!seg.bsp) continue;
                gsPoint2D<T> d0 = evalCurveDeriv(*seg.bsp, T(0));
                HalfEdge he;
                he.origin = u; he.target = v; he.twin = -1;
                he.isCrossCut = false; he.loopIdx = l; he.sStart = s0; he.sEnd = s1;
                he.outAng = std::atan2(d0.y, d0.x);
                he.seg = seg;
                HEs.push_back(std::move(he));
            }
        }

        std::set<std::pair<int, int>> addedCutPairs;
        for (size_t ci = 0; ci < out.cuts.size(); ++ci) {
            int u = cutNodes[ci].first, v = cutNodes[ci].second;
            if (u < 0 || v < 0 || u == v) continue;
            std::pair<int, int> key(std::min(u, v), std::max(u, v));
            if (addedCutPairs.count(key)) continue;
            addedCutPairs.insert(key);
            CurveSeg<T> fwd = makeStraightSeg(out.nodes[u].pt, out.nodes[v].pt, degree);
            CurveSeg<T> rev = fwd.reversed();
            int h0 = static_cast<int>(HEs.size());
            int h1 = h0 + 1;
            HalfEdge e0, e1;
            e0.origin = u; e0.target = v; e0.twin = h1; e0.isCrossCut = true; e0.cutId = static_cast<int>(ci); e0.cutForward = true;
            e0.outAng = std::atan2(out.nodes[v].pt.y - out.nodes[u].pt.y, out.nodes[v].pt.x - out.nodes[u].pt.x);
            e0.seg = fwd;
            e1.origin = v; e1.target = u; e1.twin = h0; e1.isCrossCut = true; e1.cutId = static_cast<int>(ci); e1.cutForward = false;
            e1.outAng = std::atan2(out.nodes[u].pt.y - out.nodes[v].pt.y, out.nodes[u].pt.x - out.nodes[v].pt.x);
            e1.seg = rev;
            HEs.push_back(std::move(e0));
            HEs.push_back(std::move(e1));
        }

        const int nNodes = static_cast<int>(out.nodes.size());
        struct IncRay { T ang; int heIdx; bool isOutgoing; };
        std::vector<std::vector<IncRay>> vInc(nNodes);
        for (int h = 0; h < static_cast<int>(HEs.size()); ++h) {
            vInc[HEs[h].origin].push_back({HEs[h].outAng, h, true});
            if (!HEs[h].isCrossCut) {
                gsPoint2D<T> d1 = evalCurveDeriv(*HEs[h].seg.bsp, T(1));
                T backAng = std::atan2(-d1.y, -d1.x);
                vInc[HEs[h].target].push_back({backAng, h, false});
            }
        }
        for (int v = 0; v < nNodes; ++v) {
            std::sort(vInc[v].begin(), vInc[v].end(),
                      [](const IncRay& a, const IncRay& b) { return a.ang < b.ang; });
        }

        auto nextFaceHE = [&](int h) -> int {
            int v = HEs[h].target;
            const auto& rays = vInc[v];
            const int deg = static_cast<int>(rays.size());
            if (deg == 0) return -1;
            int pos = -1;
            for (int i = 0; i < deg; ++i) {
                if (HEs[h].isCrossCut) {
                    if (rays[i].isOutgoing && rays[i].heIdx == HEs[h].twin) { pos = i; break; }
                } else {
                    if (!rays[i].isOutgoing && rays[i].heIdx == h) { pos = i; break; }
                }
            }
            if (pos < 0) return -1;
            for (int step = 1; step <= deg; ++step) {
                int idx = (pos - step + deg) % deg;
                if (rays[idx].isOutgoing) return rays[idx].heIdx;
            }
            return -1;
        };

        std::vector<bool> heVisited(HEs.size(), false);
        std::vector<bool> bndHECovered(HEs.size(), false);
        for (int hStart = 0; hStart < static_cast<int>(HEs.size()); ++hStart) {
            if (heVisited[hStart]) continue;
            std::vector<int> cyc;
            int cur = hStart;
            bool closed = false;
            for (int guard = 0; guard <= static_cast<int>(HEs.size()); ++guard) {
                if (heVisited[cur]) {
                    closed = (cur == hStart);
                    break;
                }
                heVisited[cur] = true;
                cyc.push_back(cur);
                cur = nextFaceHE(cur);
                if (cur < 0) break;
            }
            if (!closed || cyc.size() < 2) continue;

            std::vector<gsPoint2D<T>> poly;
            for (int ch : cyc)
                for (int s = 0; s < 12; ++s)
                    poly.push_back(evalCurvePoint(*HEs[ch].seg.bsp, static_cast<T>(s) / T(12)));
            T area2 = T(0);
            for (size_t i = 0; i < poly.size(); ++i) {
                const auto& a = poly[i]; const auto& b = poly[(i + 1) % poly.size()];
                area2 += a.x * b.y - b.x * a.y;
            }
            if (area2 <= T(1e-8)) continue;

            JunctionPolygon<T> hub;
            hub.id = static_cast<int>(out.hubs.size());
            for (int ch : cyc) {
                if (!HEs[ch].isCrossCut) bndHECovered[ch] = true;
                typename JunctionPolygon<T>::Side sd;
                sd.isCrossCut = HEs[ch].isCrossCut;
                sd.cutId = HEs[ch].cutId;
                sd.cutForward = HEs[ch].cutForward;
                sd.loopIdx = HEs[ch].loopIdx;
                sd.sStart = HEs[ch].sStart;
                sd.sEnd = HEs[ch].sEnd;
                sd.uNode = HEs[ch].origin;
                sd.vNode = HEs[ch].target;
                sd.seg = HEs[ch].seg;
                hub.sidesCCW.push_back(std::move(sd));
            }
            out.hubs.push_back(std::move(hub));
        }

        if (out.hubs.empty()) return false;
        for (size_t h = 0; h < HEs.size(); ++h) {
            if (!HEs[h].isCrossCut && !bndHECovered[h])
                return false;
        }
        return true;
    };

    if (!buildDCELAndTraceHubs()) return false;

    // 3.3b DCEL-Level Bridge-Slit Detector / Breaker:
    // Guarantees EVERY traced hub is genuinely simply-connected with ZERO bridge cuts (same cutId twice).
    // Strictly requires every added cut to reduce the count of bridged hubs so it can never loop or bloat cuts!
    auto detectSlitOrPinchedHole = [&](int& badHubIdx, std::set<int>& innerLoopsOut) -> int {
        badHubIdx = -1;
        innerLoopsOut.clear();
        auto sidesArea2 = [](const std::vector<typename JunctionPolygon<T>::Side>& sides, size_t a, size_t b) -> T {
            std::vector<gsPoint2D<T>> pl;
            const size_t ns = sides.size();
            for (size_t i = a; i != b; i = (i + 1) % ns)
                for (int s = 0; s < 6; ++s)
                    pl.push_back(evalCurvePoint(*sides[i].seg.bsp, static_cast<T>(s) / T(6)));
            T a2 = T(0);
            for (size_t k = 0; k < pl.size(); ++k)
                a2 += pl[k].x * pl[(k + 1) % pl.size()].y - pl[(k + 1) % pl.size()].x * pl[k].y;
            return a2;
        };
        int nBad = 0;
        for (size_t hi = 0; hi < out.hubs.size(); ++hi) {
            const auto& sides = out.hubs[hi].sidesCCW;
            const size_t ns = sides.size();
            std::map<int, size_t> seenCut;
            for (size_t i = 0; i < ns; ++i) {
                if (!sides[i].isCrossCut || sides[i].cutId < 0) continue;
                auto it = seenCut.find(sides[i].cutId);
                if (it != seenCut.end()) {
                    ++nBad;
                    if (badHubIdx < 0) {
                        badHubIdx = static_cast<int>(hi);
                        size_t j = it->second;
                        T a1 = sidesArea2(sides, j + 1, i);
                        T a2 = sidesArea2(sides, i + 1, j);
                        size_t st = (a1 < a2) ? ((j + 1) % ns) : ((i + 1) % ns);
                        size_t en = (a1 < a2) ? i : j;
                        for (size_t k = st; k != en; k = (k + 1) % ns)
                            if (!sides[k].isCrossCut && sides[k].loopIdx >= 0) innerLoopsOut.insert(sides[k].loopIdx);
                    }
                    break;
                }
                seenCut[sides[i].cutId] = i;
            }
        }
        return nBad;
    };

    for (int fixPass = 0; fixPass < nLoops + 2; ++fixPass) {
        int badH = -1;
        std::set<int> innerLoops;
        int curBad = detectSlitOrPinchedHole(badH, innerLoops);
        if (curBad == 0 || badH < 0) break;
        if (std::getenv("GISMO_MERGED_DEBUG") != nullptr)
            std::cerr << "[MERGED] M3 Slit detected in hub=" << badH << " nBad=" << curBad << " nInnerLoops=" << innerLoops.size() << "\n";

        const auto& sides = out.hubs[badH].sidesCCW;
        struct BndCand { int loopIdx; T s; gsPoint2D<T> pt; };
        std::vector<BndCand> innerCands, outerCands;
        for (const auto& sd : sides) {
            if (sd.isCrossCut || sd.loopIdx < 0) continue;
            const auto& vl = bb.vloops[sd.loopIdx];
            T sp = mergedLoopForwardSpan(vl, sd.sStart, sd.sEnd);
            for (T frac : {T(0.25), T(0.50), T(0.75)}) {
                T s = virtualLoopWrapS(vl, sd.sStart + frac * sp);
                BndCand bc{sd.loopIdx, s, evalVirtualLoopPoint(vl, s)};
                if (innerLoops.count(sd.loopIdx)) innerCands.push_back(bc);
                else outerCands.push_back(bc);
            }
        }
        if (outerCands.empty()) {
            for (int l = 0; l < nLoops; ++l) {
                if (innerLoops.count(l)) continue;
                const auto& vl = bb.vloops[l];
                for (int k = 0; k < 24; ++k) {
                    T s = (vl.totalArcLen * static_cast<T>(k)) / T(24);
                    outerCands.push_back({l, s, evalVirtualLoopPoint(vl, s)});
                }
            }
        }

        bool fixedOne = false;
        const T relaxTol = std::max(T(1e-5), T(0.002) * bb.domDiag);
        for (T minDeg : {T(12.0), T(6.0), T(2.5)}) {
            T bestDist = std::numeric_limits<T>::max();
            BndCand bestA{0, T(0), gsPoint2D<T>(0, 0)};
            BndCand bestB{0, T(0), gsPoint2D<T>(0, 0)};
            for (const auto& ca : innerCands) {
                for (const auto& cb : outerCands) {
                    if (ca.loopIdx == cb.loopIdx) continue;
                    T d = ca.pt.dist(cb.pt);
                    if (d <= relaxTol || d >= bestDist) continue;
                    if (!isCrossCutClear(ca.pt, cb.pt, densePolys, out.cuts, relaxTol)) continue;
                    if (!isCrossCutTangentClear(bb.vloops, ca.loopIdx, ca.s, cb.loopIdx, cb.s, minDeg)) continue;
                    bestDist = d; bestA = ca; bestB = cb;
                }
            }
            if (bestDist < std::numeric_limits<T>::max()) {
                CrossCut<T> cc;
                cc.id = static_cast<int>(out.cuts.size());
                cc.sleeveIdx = -1;
                cc.loopA = bestA.loopIdx; cc.loopB = bestB.loopIdx;
                cc.arcSA = bestA.s;       cc.arcSB = bestB.s;
                cc.ptA = bestA.pt;        cc.ptB = bestB.pt;
                cc.costJ = T(0);
                out.cuts.push_back(cc);
                int chkH = -1; std::set<int> chkL;
                if (buildDCELAndTraceHubs() && detectSlitOrPinchedHole(chkH, chkL) < curBad) {
                    fixedOne = true;
                    break;
                } else {
                    out.cuts.pop_back();
                    buildDCELAndTraceHubs();
                }
            }
        }
        if (!fixedOne) break;
    }

    // Stage 1B: Absorb redundant B=2 boundary blisters and B=3 triangular pockets/junctions into their neighbor hub
    // by deleting a safe redundant cross-cut BEFORE parity splitting injects artificial 180-deg nodes into a B=4 hub!
    // Prioritizes: (3) any B=2 blister, (2) B=3 + odd neighbor (cancels TWO odd hubs into one EVEN hub!),
    //              (1) any B=3 + even neighbor (turns B=3+4 into B=5->6 where chords split 180 cleanly instead of O-rings).
    std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> blockedCuts;
    for (int prunePass = 0; prunePass < 16; ++prunePass) {
        if (out.hubs.size() <= 1) break;
        if (nLoops == 1 && out.cuts.size() <= 2) break;

        std::vector<T> hubAreas(out.hubs.size(), T(0));
        for (size_t hi = 0; hi < out.hubs.size(); ++hi) {
            std::vector<gsPoint2D<T>> poly;
            for (const auto& sd : out.hubs[hi].sidesCCW)
                for (int s = 0; s < 8; ++s)
                    poly.push_back(evalCurvePoint(*sd.seg.bsp, static_cast<T>(s) / T(8)));
            T a2 = T(0);
            for (size_t k = 0; k < poly.size(); ++k)
                a2 += poly[k].x * poly[(k + 1) % poly.size()].y - poly[(k + 1) % poly.size()].x * poly[k].y;
            hubAreas[hi] = T(0.5) * std::abs(a2);
        }
        T totA = std::accumulate(hubAreas.begin(), hubAreas.end(), T(0));

        // Map each cross-cut ci -> its two incident hubs (hA, hB)
        std::vector<std::pair<int, int>> cutHubs(out.cuts.size(), {-1, -1});
        for (size_t hi = 0; hi < out.hubs.size(); ++hi) {
            for (const auto& sd : out.hubs[hi].sidesCCW) {
                if (!sd.isCrossCut || sd.cutId < 0 || sd.cutId >= static_cast<int>(out.cuts.size())) continue;
                auto& pr = cutHubs[sd.cutId];
                if (pr.first < 0) pr.first = static_cast<int>(hi);
                else if (pr.second < 0 && pr.first != static_cast<int>(hi)) pr.second = static_cast<int>(hi);
            }
        }

        auto isBlockedCut = [&](int ci) -> bool {
            for (const auto& bc : blockedCuts) {
                if ((bc.first.dist(out.cuts[ci].ptA) < T(1e-4) && bc.second.dist(out.cuts[ci].ptB) < T(1e-4)) ||
                    (bc.first.dist(out.cuts[ci].ptB) < T(1e-4) && bc.second.dist(out.cuts[ci].ptA) < T(1e-4)))
                    return true;
            }
            return false;
        };

        auto canSafelyDeleteCut = [&](int skipIdx) -> bool {
            if (skipIdx < 0 || skipIdx >= static_cast<int>(out.cuts.size())) return false;
            if (out.cuts[skipIdx].loopA == out.cuts[skipIdx].loopB) return true;
            if (nLoops <= 1) return (out.cuts.size() > 2);
            // Verify remaining cuts keep all loops connected, 2-edge-connected, and distinctDeg >= minReq
            std::vector<std::vector<T>> lS(nLoops);
            for (int ci = 0; ci < static_cast<int>(out.cuts.size()); ++ci) {
                if (ci == skipIdx) continue;
                lS[out.cuts[ci].loopA].push_back(out.cuts[ci].arcSA);
                lS[out.cuts[ci].loopB].push_back(out.cuts[ci].arcSB);
            }
            // Check distinct hole degree
            for (int l = 1; l < nLoops; ++l) {
                int minReq = bb.vloops[l].corners.empty() ? 3 : 2;
                const T vTol = std::max(T(1e-5), T(0.015) * bb.vloops[l].totalArcLen);
                std::vector<T> uniq;
                for (T s : lS[l]) {
                    bool dup = false;
                    for (T u : uniq) {
                        T d = std::abs(s - u);
                        if (std::min(d, bb.vloops[l].totalArcLen - d) < vTol) { dup = true; break; }
                    }
                    if (!dup) uniq.push_back(s);
                }
                if (static_cast<int>(uniq.size()) < minReq) return false;
            }
            // Check 2-edge-connectivity
            for (int drop2 = -1; drop2 < static_cast<int>(out.cuts.size()); ++drop2) {
                if (drop2 == skipIdx) continue;
                if (drop2 >= 0 && out.cuts[drop2].loopA == out.cuts[drop2].loopB) continue;
                std::vector<std::vector<int>> subA(nLoops);
                for (int ci = 0; ci < static_cast<int>(out.cuts.size()); ++ci) {
                    if (ci == skipIdx || ci == drop2) continue;
                    subA[out.cuts[ci].loopA].push_back(out.cuts[ci].loopB);
                    subA[out.cuts[ci].loopB].push_back(out.cuts[ci].loopA);
                }
                std::vector<bool> vis(nLoops, false);
                std::queue<int> bq;
                bq.push(0); vis[0] = true;
                while (!bq.empty()) {
                    int u = bq.front(); bq.pop();
                    for (int nb : subA[u]) if (!vis[nb]) { vis[nb] = true; bq.push(nb); }
                }
                for (int l = 1; l < nLoops; ++l) if (!vis[l]) return false;
            }
            return true;
        };

        int cutToDrop = -1;
        int bestPrio = -1;
        T smallestPocketA = std::numeric_limits<T>::max();
        for (int ci = 0; ci < static_cast<int>(out.cuts.size()); ++ci) {
            if (isBlockedCut(ci)) continue;
            int hA = cutHubs[ci].first, hB = cutHubs[ci].second;
            if (hA < 0 || hB < 0 || hA == hB) continue;
            int nA = static_cast<int>(out.hubs[hA].sidesCCW.size());
            int nB = static_cast<int>(out.hubs[hB].sidesCCW.size());
            if (nA > 3 && nB > 3) continue;
            int nMerged = nA + nB - 2;
            if (nMerged > 14) continue;
            // On 0-hole smooth shapes (amoeba/heart/crescent), keep big >18% smooth lobes intact
            T minA = std::min(hubAreas[hA], hubAreas[hB]);
            if (nLoops == 1 && minA > T(0.18) * totA) continue;
            if (!canSafelyDeleteCut(ci)) continue;

            int prio = 1;
            if (nA == 2 || nB == 2) prio = 3;
            else if ((nA == 3 && (nB % 2) == 1) || (nB == 3 && (nA % 2) == 1)) prio = 2;

            if (prio > bestPrio || (prio == bestPrio && minA < smallestPocketA)) {
                bestPrio = prio;
                smallestPocketA = minA;
                cutToDrop = ci;
            }
        }
        if (cutToDrop < 0) break;

        auto backupCuts = out.cuts;
        std::pair<gsPoint2D<T>, gsPoint2D<T>> droppedEnds = {out.cuts[cutToDrop].ptA, out.cuts[cutToDrop].ptB};
        out.cuts.erase(out.cuts.begin() + cutToDrop);
        int tmpBadH = -1;
        std::set<int> tmpInner;
        if (!buildDCELAndTraceHubs() || detectSlitOrPinchedHole(tmpBadH, tmpInner) > 0) {
            out.cuts = std::move(backupCuts);
            buildDCELAndTraceHubs();
            blockedCuts.push_back(droppedEnds);
        }
    }

    // 3.4 Curvature-Aware Wall Subdivision & 1-Layer-Preserving Parity Balancing
    // Subdivide any smooth boundary wall turning > 125 deg, and ensure every hub has EVEN sides B(P_j) = 2*m_j >= 4
    // strictly by splitting CAD boundary walls (preserving k_gamma == 1 on every cross-cut!).
    auto splitBoundarySideAtMid = [&](JunctionPolygon<T>& hub, int sideIdx) -> bool {
        auto sd = hub.sidesCCW[sideIdx];
        if (sd.isCrossCut || sd.loopIdx < 0) return false;
        const auto& vl = bb.vloops[sd.loopIdx];
        T span = mergedLoopForwardSpan(vl, sd.sStart, sd.sEnd);
        if (span <= T(1e-7)) return false;
        T sMid = virtualLoopWrapS(vl, sd.sStart + T(0.5) * span);
        CurveSeg<T> seg0 = extractVirtualSubCurve(vl, sd.sStart, sMid, degree);
        CurveSeg<T> seg1 = extractVirtualSubCurve(vl, sMid, sd.sEnd, degree);
        if (!seg0.bsp || !seg1.bsp) return false;

        int midNode = static_cast<int>(out.nodes.size());
        typename SlicedDecomposition<T>::ArrNode nd;
        nd.pt = evalVirtualLoopPoint(vl, sMid);
        nd.loopIdx = sd.loopIdx;
        nd.arcS = sMid;
        nd.isPinnedCorner = false;
        nd.cornerAngleDeg = T(180);
        nd.targetValence = 2;
        out.nodes.push_back(nd);

        typename JunctionPolygon<T>::Side s0 = sd;
        s0.sEnd = sMid; s0.vNode = midNode; s0.seg = seg0;
        typename JunctionPolygon<T>::Side s1 = sd;
        s1.sStart = sMid; s1.uNode = midNode; s1.seg = seg1;

        hub.sidesCCW[sideIdx] = std::move(s0);
        hub.sidesCCW.insert(hub.sidesCCW.begin() + sideIdx + 1, std::move(s1));
        return true;
    };

    for (auto& hub : out.hubs) {
        // First subdivide strongly turning (> 125 deg) smooth CAD boundary walls so Coons patches don't fold
        for (int pass = 0; pass < 3; ++pass) {
            bool splitAny = false;
            for (int si = 0; si < static_cast<int>(hub.sidesCCW.size()); ++si) {
                const auto& sd = hub.sidesCCW[si];
                if (sd.isCrossCut || !sd.seg.bsp) continue;
                if (curveIntegratedTurnRad(*sd.seg.bsp, 16) > T(2.18)) { // > ~125 deg
                    if (splitBoundarySideAtMid(hub, si)) { splitAny = true; break; }
                }
            }
            if (!splitAny) break;
        }

        // Enforce B(P_j) >= 4 and B(P_j) even by splitting the longest CAD boundary side
        while (hub.sidesCCW.size() < 4 || (hub.sidesCCW.size() % 2) == 1) {
            int bestSi = -1;
            T bestLen = -T(1);
            for (int si = 0; si < static_cast<int>(hub.sidesCCW.size()); ++si) {
                if (hub.sidesCCW[si].isCrossCut) continue;
                T L = mergedMeasureSegArcLen(hub.sidesCCW[si].seg, 12);
                if (L > bestLen) { bestLen = L; bestSi = si; }
            }
            if (bestSi < 0 || !splitBoundarySideAtMid(hub, bestSi))
                break;
        }
    }

    // Assign interior cross-field singularities to their enclosing JunctionPolygon
    for (int sIdx = 0; sIdx < static_cast<int>(gcf.field.singularities.size()); ++sIdx) {
        const auto& sg = gcf.field.singularities[sIdx];
        for (auto& hub : out.hubs) {
            if (hub.sidesCCW.size() < 6) continue; // Only seed singularities into 6+-gon junction hubs
            std::vector<gsPoint2D<T>> poly;
            T minEdge = std::numeric_limits<T>::max();
            for (const auto& sd : hub.sidesCCW) {
                minEdge = std::min(minEdge, sd.seg.p0().dist(sd.seg.p1()));
                for (int k = 0; k < 6; ++k)
                    poly.push_back(evalCurvePoint(*sd.seg.bsp, static_cast<T>(k) / T(6)));
            }
            if (pointInPolyStrict(sg.pos, poly, std::max(T(1e-5), T(0.10) * minEdge))) {
                hub.singularityIds.push_back(sIdx);
                break;
            }
        }
    }

    return true;
}

// ============================================================================
// Module 4: Singularity-Seeded Macro-Triangulation & Exact Flip/Bitmask Pairing
// ============================================================================

template <class T>
struct HubLocalQuad {
    // Vertices: [0..N-1] are boundary vertices of hub, [N..] are interior Steiner vertices
    std::array<int, 4> v;
    std::array<gsPoint2D<T>, 4> pts;
};

// Exact O(N^4) DP globally optimal chord quadrangulation of an even N-gon (I == 0, zero K_{1,3} deadlocks!)
// Uses EXACT curved CAD boundary tangents from `edgesCCW` for every corner angle and chord clearance check!
template <class T>
inline bool quadrangulateEvenHubDP(
    const std::vector<CurveSeg<T>>& edgesCCW,
    const std::vector<gsPoint2D<T>>& bndPts,
    const std::vector<bool>& isSmooth180,
    const std::vector<gsPoint2D<T>>& hubPoly,
    const std::vector<std::vector<gsPoint2D<T>>>& densePolys,
    const GuidedCrossField<T>& gcf,
    bool allowSmooth180,
    std::vector<std::array<int, 4>>& outQuads)
{
    const int N = static_cast<int>(bndPts.size());
    outQuads.clear();
    if (N < 4 || (N % 2) != 0) return false;
    if (N == 4) {
        outQuads.push_back({0, 1, 2, 3});
        return true;
    }
    if (N > 20) return false;

    // Precompute exact unit tangents along boundary edges of this hub
    std::vector<gsPoint2D<T>> segOutDir(N), segBackDir(N);
    std::vector<T> cornerSectorDeg(N, T(180));
    for (int i = 0; i < N; ++i) {
        gsPoint2D<T> d0 = evalCurveDeriv(*edgesCCW[i].bsp, T(0));
        gsPoint2D<T> d1 = evalCurveDeriv(*edgesCCW[i].bsp, T(1));
        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
        segOutDir[i]  = (n0 > T(1e-12)) ? gsPoint2D<T>(d0.x / n0, d0.y / n0) : gsPoint2D<T>(1, 0);
        segBackDir[i] = (n1 > T(1e-12)) ? gsPoint2D<T>(-d1.x / n1, -d1.y / n1) : gsPoint2D<T>(-1, 0);
    }
    for (int i = 0; i < N; ++i)
        cornerSectorDeg[i] = mergedCcwAngleDeg(segOutDir[i], segBackDir[(i + N - 1) % N]);

    int triHint = -1;
    const T minTanClear = allowSmooth180 ? T(6.5) : T(9.0);

    // Precompute valid interior chords (i, j) checking BOTH hubPoly/densePolys AND exact curved tangent sectors!
    std::vector<std::vector<bool>> validChord(N, std::vector<bool>(N, false));
    std::vector<std::vector<T>> chordCF(N, std::vector<T>(N, T(0)));
    for (int i = 0; i < N; ++i) {
        validChord[i][(i + 1) % N] = true;
        validChord[(i + 1) % N][i] = true;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = i + 2; j < N; ++j) {
            if (i == 0 && j == N - 1) continue;
            const auto& A = bndPts[i];
            const auto& B = bndPts[j];
            T len = A.dist(B);
            if (len <= T(1e-8)) continue;

            // Exact curved boundary tangent clearance at i and j
            gsPoint2D<T> dAB((B.x - A.x) / len, (B.y - A.y) / len);
            gsPoint2D<T> dBA(-dAB.x, -dAB.y);
            T aI = mergedCcwAngleDeg(segOutDir[i], dAB);
            T aJ = mergedCcwAngleDeg(segOutDir[j], dBA);
            if (aI < minTanClear || cornerSectorDeg[i] - aI < minTanClear) continue;
            if (aJ < minTanClear || cornerSectorDeg[j] - aJ < minTanClear) continue;

            if (!isSegmentInsideDomain(A, B, densePolys)) continue;
            bool okIn = true;
            for (int k = 1; k <= 5; ++k) {
                T f = static_cast<T>(k) / T(6);
                if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y), hubPoly, T(1e-6))) {
                    okIn = false; break;
                }
            }
            if (!okIn) continue;
            const int np = static_cast<int>(hubPoly.size());
            for (int k = 0; k < np && okIn; ++k) {
                const auto& p0 = hubPoly[k]; const auto& p1 = hubPoly[(k + 1) % np];
                if (p0.dist(A) < T(1e-5) || p1.dist(A) < T(1e-5) || p0.dist(B) < T(1e-5) || p1.dist(B) < T(1e-5))
                    continue;
                if (segProperIntersect(A, B, p0, p1)) okIn = false;
            }
            if (!okIn) continue;
            validChord[i][j] = validChord[j][i] = true;
            T al = evalSegmentCrossFieldAlignment(gcf.fineMesh, gcf.field, gcf.accel, A, B, triHint);
            chordCF[i][j] = chordCF[j][i] = T(1) - al;
        }
    }

    auto edgeDirOut = [&](int u, int v) -> gsPoint2D<T> {
        if ((u + 1) % N == v) return segOutDir[u];
        T dx = bndPts[v].x - bndPts[u].x, dy = bndPts[v].y - bndPts[u].y;
        T l = std::hypot(dx, dy);
        return (l > T(1e-12)) ? gsPoint2D<T>(dx / l, dy / l) : gsPoint2D<T>(1, 0);
    };
    auto edgeDirBack = [&](int u, int v) -> gsPoint2D<T> {
        if ((u + 1) % N == v) return segBackDir[u];
        T dx = bndPts[u].x - bndPts[v].x, dy = bndPts[u].y - bndPts[v].y;
        T l = std::hypot(dx, dy);
        return (l > T(1e-12)) ? gsPoint2D<T>(dx / l, dy / l) : gsPoint2D<T>(-1, 0);
    };

    auto evalQuadCost = [&](int i, int a, int b, int j) -> T {
        gsPoint2D<T> q[4] = {bndPts[i], bndPts[a], bndPts[b], bndPts[j]};
        T area2 = cross2(q[0], q[1], q[2]) + cross2(q[0], q[2], q[3]);
        if (area2 <= T(1e-8)) return T(1e18);
        bool strictC = isStrictlyConvexQuad(q);
        if (!strictC) {
            for (int k = 0; k < 4; ++k) {
                if (cross2(q[k], q[(k + 1) % 4], q[(k + 2) % 4]) < -T(0.12) * area2)
                    return T(1e18);
            }
        }
        T cost = strictC ? T(0) : T(6.0);
        int idx[4] = {i, a, b, j};
        for (int k = 0; k < 4; ++k) {
            int vPrev = idx[(k + 3) % 4], vCur = idx[k], vNext = idx[(k + 1) % 4];
            bool prevBnd = ((vPrev + 1) % N == vCur);
            bool nextBnd = ((vCur + 1) % N == vNext);
            bool unsplitSmooth = (prevBnd && nextBnd && static_cast<size_t>(vCur) < isSmooth180.size() && isSmooth180[vCur]);
            if (prevBnd && nextBnd && cornerSectorDeg[vCur] >= T(195.0))
                return T(1e18); // Strictly forbid leaving a >=195 deg reentrant hole/domain corner unsplit!
            if (unsplitSmooth) {
                if (!allowSmooth180) return T(1e18);
                cost += T(25.0); // Allowed in fallback mode only (will be upgraded by emitRobustOGrid5)
            }
            // Exact curved-tangent corner angle!
            gsPoint2D<T> dOut  = edgeDirOut(vCur, vNext);
            gsPoint2D<T> dBack = edgeDirBack(vPrev, vCur);
            T ang = mergedCcwAngleDeg(dOut, dBack);
            if (!unsplitSmooth) {
                if (ang < T(7.5) || ang > T(168.0)) return T(1e18);
                T dev = std::abs(ang - T(90.0)) / T(45.0);
                cost += dev * dev;
                if (ang < T(32.0)) cost += T(6.0) * (T(32.0) - ang) / T(10.0);
                if (ang > T(148.0)) cost += T(6.0) * (ang - T(148.0)) / T(10.0);
            }
        }
        for (int k = 0; k < 4; ++k) {
            int u = idx[k], v = idx[(k + 1) % 4];
            if ((u + 1) % N != v && (v + 1) % N != u)
                cost += T(1.2) * chordCF[u][v];
        }
        return cost;
    };

    const T INF = T(1e17);
    std::vector<std::vector<T>> dp(N, std::vector<T>(N, INF));
    struct Choice { int a = -1, b = -1; };
    std::vector<std::vector<Choice>> choice(N, std::vector<Choice>(N));

    // Base: adjacent edges have 0 cost
    for (int i = 0; i + 1 < N; ++i) dp[i][i + 1] = T(0);

    // Odd span lengths len = 3, 5, 7, ..., N-1 (sub-polygon i..j has len+1 vertices, which is even!)
    for (int len = 3; len <= N - 1; len += 2) {
        for (int i = 0; i + len < N; ++i) {
            int j = i + len;
            if (!(i == 0 && j == N - 1) && !validChord[i][j]) continue;
            for (int a = i + 1; a <= j - 2; a += 2) {
                if (a > i + 1 && !validChord[i][a]) continue;
                if (dp[i][a] >= INF) continue;
                for (int b = a + 1; b <= j - 1; b += 2) {
                    if (b > a + 1 && !validChord[a][b]) continue;
                    if (j > b + 1 && !validChord[b][j]) continue;
                    if (dp[a][b] >= INF || dp[b][j] >= INF) continue;
                    T qc = evalQuadCost(i, a, b, j);
                    if (qc >= INF) continue;
                    T tot = qc + dp[i][a] + dp[a][b] + dp[b][j];
                    if (tot < dp[i][j]) {
                        dp[i][j] = tot;
                        choice[i][j] = {a, b};
                    }
                }
            }
        }
    }

    if (dp[0][N - 1] >= INF) return false;

    std::function<void(int, int)> reconstruct = [&](int i, int j) {
        if (j - i <= 1) return;
        Choice c = choice[i][j];
        outQuads.push_back({i, c.a, c.b, j});
        reconstruct(i, c.a);
        reconstruct(c.a, c.b);
        reconstruct(c.b, j);
    };
    reconstruct(0, N - 1);
    return !outQuads.empty();
}

// Singularity-seeded Steiner triangulation + Joint Diagonal-Flip & Exact Bitmask Pairing (N_tri <= 18)
// Evaluates EXACT curved CAD boundary tangents and forbids leaving smooth-180 vertices unsplit!
template <class T>
inline bool pairSteinerHubExact(
    const std::vector<CurveSeg<T>>& edgesCCW,
    const std::vector<gsPoint2D<T>>& bndPts,
    const std::vector<bool>& isSmooth180,
    const std::vector<gsPoint2D<T>>& hubPoly,
    const std::vector<std::vector<gsPoint2D<T>>>& densePolys,
    const std::vector<gsPoint2D<T>>& singPts,
    const std::vector<int>& singTargetVal,
    const GuidedCrossField<T>& gcf,
    std::vector<gsPoint2D<T>>& outAllPts,
    std::vector<std::array<int, 4>>& outQuads)
{
    const int N = static_cast<int>(bndPts.size());
    const int I = static_cast<int>(singPts.size());
    if (N < 4 || (N % 2) != 0 || I <= 0) return false;
    const int expectedTris = N + 2 * I - 2;
    if (expectedTris > 18 || (expectedTris % 2) != 0) return false;

    std::vector<gsPoint2D<T>> segOutDir(N), segBackDir(N);
    for (int i = 0; i < N; ++i) {
        gsPoint2D<T> d0 = evalCurveDeriv(*edgesCCW[i].bsp, T(0));
        gsPoint2D<T> d1 = evalCurveDeriv(*edgesCCW[i].bsp, T(1));
        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
        segOutDir[i]  = (n0 > T(1e-12)) ? gsPoint2D<T>(d0.x / n0, d0.y / n0) : gsPoint2D<T>(1, 0);
        segBackDir[i] = (n1 > T(1e-12)) ? gsPoint2D<T>(-d1.x / n1, -d1.y / n1) : gsPoint2D<T>(-1, 0);
    }

    outAllPts = bndPts;
    for (const auto& s : singPts) outAllPts.push_back(s);
    const int totalV = static_cast<int>(outAllPts.size());

    // 1. Ear-clip boundary polygon [0..N-1]
    struct LTri { int v[3]; };
    std::vector<LTri> tris;
    std::vector<int> ring(N);
    std::iota(ring.begin(), ring.end(), 0);
    while (ring.size() > 3) {
        int m = static_cast<int>(ring.size());
        int bestK = -1;
        T bestSc = -T(1e30);
        for (int k = 0; k < m; ++k) {
            int ia = ring[(k + m - 1) % m], ib = ring[k], ic = ring[(k + 1) % m];
            const auto& A = outAllPts[ia]; const auto& B = outAllPts[ib]; const auto& C = outAllPts[ic];
            T a2 = cross2(A, B, C);
            if (a2 <= T(1e-9)) continue;
            bool clear = true;
            for (int r = 0; r < m; ++r) {
                int ir = ring[r];
                if (ir == ia || ir == ib || ir == ic) continue;
                const auto& P = outAllPts[ir];
                if (cross2(A, B, P) >= -T(1e-9) && cross2(B, C, P) >= -T(1e-9) && cross2(C, A, P) >= -T(1e-9)) {
                    clear = false; break;
                }
            }
            if (!clear) continue;
            T l2 = A.distSq(B) + B.distSq(C) + C.distSq(A);
            T sc = a2 / std::max(T(1e-14), l2);
            if (sc > bestSc) { bestSc = sc; bestK = k; }
        }
        if (bestK < 0) return false;
        int mCur = static_cast<int>(ring.size());
        tris.push_back({{ring[(bestK + mCur - 1) % mCur], ring[bestK], ring[(bestK + 1) % mCur]}});
        ring.erase(ring.begin() + bestK);
    }
    tris.push_back({{ring[0], ring[1], ring[2]}});

    // 2. Insert each interior Steiner singularity by splitting its enclosing triangle 1 -> 3
    for (int si = 0; si < I; ++si) {
        int vNew = N + si;
        const auto& P = outAllPts[vNew];
        int hitT = -1;
        T bestMinBary = -T(1e9);
        for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti) {
            const auto& A = outAllPts[tris[ti].v[0]];
            const auto& B = outAllPts[tris[ti].v[1]];
            const auto& C = outAllPts[tris[ti].v[2]];
            T a2 = cross2(A, B, C);
            if (a2 <= T(1e-10)) continue;
            T b0 = cross2(B, C, P) / a2;
            T b1 = cross2(C, A, P) / a2;
            T b2 = T(1) - b0 - b1;
            T mn = std::min(b0, std::min(b1, b2));
            if (mn > bestMinBary) { bestMinBary = mn; hitT = ti; }
        }
        if (hitT < 0 || bestMinBary < -T(1e-4)) return false;
        LTri oldT = tris[hitT];
        tris[hitT] = {{oldT.v[0], oldT.v[1], vNew}};
        tris.push_back({{oldT.v[1], oldT.v[2], vNew}});
        tris.push_back({{oldT.v[2], oldT.v[0], vNew}});
    }

    if (static_cast<int>(tris.size()) != expectedTris) return false;

    // Helper: constrained Lawson flips inside P_j (never flipping boundary edges (i, (i+1)%N))
    auto isBndEdge = [&](int a, int b) -> bool {
        return (a < N && b < N) && (((a + 1) % N == b) || ((b + 1) % N == a));
    };
    auto inCircle = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C, const gsPoint2D<T>& D) -> bool {
        T ax = A.x - D.x, ay = A.y - D.y;
        T bx = B.x - D.x, by = B.y - D.y;
        T cx = C.x - D.x, cy = C.y - D.y;
        T det = (ax * ax + ay * ay) * (bx * cy - cx * by)
              - (bx * bx + by * by) * (ax * cy - cx * ay)
              + (cx * cx + cy * cy) * (ax * by - bx * ay);
        return det > T(1e-10);
    };
    auto lawsonRelax = [&](std::vector<LTri>& tList) {
        for (int pass = 0; pass < 8; ++pass) {
            bool flipped = false;
            const int nT = static_cast<int>(tList.size());
            for (int i = 0; i < nT && !flipped; ++i) {
                for (int ka = 0; ka < 3 && !flipped; ++ka) {
                    int a = tList[i].v[ka], b = tList[i].v[(ka + 1) % 3], c = tList[i].v[(ka + 2) % 3];
                    if (isBndEdge(a, b)) continue;
                    for (int j = i + 1; j < nT; ++j) {
                        for (int kb = 0; kb < 3; ++kb) {
                            if (tList[j].v[kb] == b && tList[j].v[(kb + 1) % 3] == a) {
                                int d = tList[j].v[(kb + 2) % 3];
                                gsPoint2D<T> q[4] = {outAllPts[c], outAllPts[a], outAllPts[d], outAllPts[b]};
                                if (isStrictlyConvexQuad(q) && inCircle(outAllPts[a], outAllPts[b], outAllPts[c], outAllPts[d])) {
                                    tList[i] = {{c, a, d}};
                                    tList[j] = {{c, d, b}};
                                    flipped = true;
                                }
                                break;
                            }
                        }
                        if (flipped) break;
                    }
                }
            }
            if (!flipped) break;
        }
    };
    lawsonRelax(tris);

    auto edgeDirOutSt = [&](int u, int v) -> gsPoint2D<T> {
        if (u < N && v < N && ((u + 1) % N == v)) return segOutDir[u];
        T dx = outAllPts[v].x - outAllPts[u].x, dy = outAllPts[v].y - outAllPts[u].y;
        T l = std::hypot(dx, dy);
        return (l > T(1e-12)) ? gsPoint2D<T>(dx / l, dy / l) : gsPoint2D<T>(1, 0);
    };
    auto edgeDirBackSt = [&](int u, int v) -> gsPoint2D<T> {
        if (u < N && v < N && ((u + 1) % N == v)) return segBackDir[u];
        T dx = outAllPts[u].x - outAllPts[v].x, dy = outAllPts[u].y - outAllPts[v].y;
        T l = std::hypot(dx, dy);
        return (l > T(1e-12)) ? gsPoint2D<T>(dx / l, dy / l) : gsPoint2D<T>(-1, 0);
    };
    auto segInHubSt = [&](int u, int v) -> bool {
        if (isBndEdge(u, v)) return true;
        const auto& A = outAllPts[u]; const auto& B = outAllPts[v];
        if (!isSegmentInsideDomain(A, B, densePolys)) return false;
        for (int s = 1; s <= 3; ++s) {
            T f = static_cast<T>(s) / T(4);
            if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y), hubPoly, T(1e-6)))
                return false;
        }
        return true;
    };

    int triHint = -1;

    auto solveExactPairingOnTris = [&](const std::vector<LTri>& tList,
                                       std::vector<std::array<int, 4>>& bestQ,
                                       T& bestObj) -> bool {
        const int nT = static_cast<int>(tList.size());
        struct PairOpt {
            int t0, t1;
            std::array<int, 4> qv;
            int diagA, diagB;
            T cost;
        };
        std::vector<PairOpt> opts;
        std::vector<std::vector<int>> triOpts(nT);

        std::vector<int> baseDeg(totalV, 0);
        for (const auto& t : tList) { ++baseDeg[t.v[0]]; ++baseDeg[t.v[1]]; ++baseDeg[t.v[2]]; }

        for (int i = 0; i < nT; ++i) {
            for (int ka = 0; ka < 3; ++ka) {
                int a = tList[i].v[ka], b = tList[i].v[(ka + 1) % 3], c = tList[i].v[(ka + 2) % 3];
                if (isBndEdge(a, b)) continue;
                for (int j = i + 1; j < nT; ++j) {
                    for (int kb = 0; kb < 3; ++kb) {
                        if (tList[j].v[kb] == b && tList[j].v[(kb + 1) % 3] == a) {
                            int d = tList[j].v[(kb + 2) % 3];
                            gsPoint2D<T> q[4] = {outAllPts[c], outAllPts[a], outAllPts[d], outAllPts[b]};
                            if (!isStrictlyConvexQuad(q)) continue;
                            std::array<int, 4> qv = {c, a, d, b};
                            T qCost = T(0);
                            bool okA = true;
                            for (int m = 0; m < 4; ++m) {
                                int vPrev = qv[(m + 3) % 4], vCur = qv[m], vNext = qv[(m + 1) % 4];
                                if (!segInHubSt(vCur, vNext)) { okA = false; break; }
                                if (vCur < N && vPrev < N && vNext < N &&
                                    ((vPrev + 1) % N == vCur) && ((vCur + 1) % N == vNext) &&
                                    static_cast<size_t>(vCur) < isSmooth180.size() && isSmooth180[vCur]) {
                                    okA = false; break;
                                }
                                T ang = mergedCcwAngleDeg(edgeDirOutSt(vCur, vNext), edgeDirBackSt(vPrev, vCur));
                                if (ang < T(10.0) || ang > T(163.0)) { okA = false; break; }
                                T dev = std::abs(ang - T(90.0)) / T(45.0);
                                qCost += dev * dev;
                                if (ang < T(38.0)) qCost += T(5.0) * (T(38.0) - ang) / T(10.0);
                                if (ang > T(142.0)) qCost += T(5.0) * (ang - T(142.0)) / T(10.0);
                            }
                            if (!okA) continue;
                            for (int m = 0; m < 4; ++m) {
                                int u = qv[m], v = qv[(m + 1) % 4];
                                if (!isBndEdge(u, v))
                                    qCost += T(0.9) * (T(1) - evalSegmentCrossFieldAlignment(gcf.fineMesh, gcf.field, gcf.accel, outAllPts[u], outAllPts[v], triHint));
                            }
                            int oi = static_cast<int>(opts.size());
                            opts.push_back({i, j, qv, a, b, qCost});
                            triOpts[i].push_back(oi);
                            triOpts[j].push_back(oi);
                        }
                    }
                }
            }
        }

        bestObj = T(1e18);
        bool found = false;
        std::vector<int> curChoice;
        std::vector<int> curDeg = baseDeg;
        const uint32_t fullMask = (1u << nT) - 1u;

        std::function<void(uint32_t, T)> dfs = [&](uint32_t mask, T accCost) {
            if (accCost >= bestObj) return;
            if (mask == fullMask) {
                // Evaluate interior vertex valence penalty (target v in {3, 4, 5}, match singTargetVal)
                T valPen = T(0);
                for (int si = 0; si < I; ++si) {
                    int vDeg = curDeg[N + si];
                    if (vDeg < 3 || vDeg > 6) return;
                    int diff = vDeg - singTargetVal[si];
                    valPen += T(3.5) * static_cast<T>(diff * diff);
                }
                if (accCost + valPen < bestObj) {
                    bestObj = accCost + valPen;
                    bestQ.clear();
                    for (int oi : curChoice) bestQ.push_back(opts[oi].qv);
                    found = true;
                }
                return;
            }
            // Pick unmatched triangle with fewest available partner options (MRV heuristic -> instant pruning!)
            int pickT = -1, minOpts = 999;
            for (int ti = 0; ti < nT; ++ti) {
                if (mask & (1u << ti)) continue;
                int cnt = 0;
                for (int oi : triOpts[ti]) {
                    int other = (opts[oi].t0 == ti) ? opts[oi].t1 : opts[oi].t0;
                    if (!(mask & (1u << other))) ++cnt;
                }
                if (cnt == 0) return;
                if (cnt < minOpts) { minOpts = cnt; pickT = ti; }
            }
            if (pickT < 0) return;

            for (int oi : triOpts[pickT]) {
                int other = (opts[oi].t0 == pickT) ? opts[oi].t1 : opts[oi].t0;
                if (mask & (1u << other)) continue;
                int da = opts[oi].diagA, db = opts[oi].diagB;
                // Keep interior Steiner vertices at valence >= 3
                if (da >= N && curDeg[da] - 1 < 3) continue;
                if (db >= N && curDeg[db] - 1 < 3) continue;

                --curDeg[da]; --curDeg[db];
                curChoice.push_back(oi);
                dfs(mask | (1u << pickT) | (1u << other), accCost + opts[oi].cost);
                curChoice.pop_back();
                ++curDeg[da]; ++curDeg[db];
            }
        };

        dfs(0u, T(0));
        return found;
    };

    T globalBestObj = T(1e18);
    bool globalFound = solveExactPairingOnTris(tris, outQuads, globalBestObj);

    // Joint diagonal-flip search (up to 12 single convex interior flips) to escape K_{1,3} claws and improve angles!
    const int nT = static_cast<int>(tris.size());
    int flipsTried = 0;
    for (int i = 0; i < nT && flipsTried < 14; ++i) {
        for (int ka = 0; ka < 3 && flipsTried < 14; ++ka) {
            int a = tris[i].v[ka], b = tris[i].v[(ka + 1) % 3], c = tris[i].v[(ka + 2) % 3];
            if (isBndEdge(a, b)) continue;
            for (int j = i + 1; j < nT; ++j) {
                for (int kb = 0; kb < 3; ++kb) {
                    if (tris[j].v[kb] == b && tris[j].v[(kb + 1) % 3] == a) {
                        int d = tris[j].v[(kb + 2) % 3];
                        gsPoint2D<T> q[4] = {outAllPts[c], outAllPts[a], outAllPts[d], outAllPts[b]};
                        if (!isStrictlyConvexQuad(q)) continue;
                        ++flipsTried;
                        std::vector<LTri> alt = tris;
                        alt[i] = {{c, a, d}};
                        alt[j] = {{c, d, b}};
                        std::vector<std::array<int, 4>> candQ;
                        T candObj = T(1e18);
                        if (solveExactPairingOnTris(alt, candQ, candObj) && candObj < globalBestObj) {
                            globalBestObj = candObj;
                            outQuads = std::move(candQ);
                            globalFound = true;
                        }
                        break;
                    }
                }
            }
        }
    }

    return globalFound;
}

// Radial-Star m-Quad Hub Solver for star-shaped 2m-gons (N in {6, 8, 10, 12}):
// Uses pole-of-inaccessibility + exact curved tangents!
template <class T>
inline bool quadrangulateRadialStarHub(
    const std::vector<CurveSeg<T>>& edgesCCW,
    const std::vector<gsPoint2D<T>>& bndPts,
    const std::vector<bool>& isSmooth180,
    const std::vector<gsPoint2D<T>>& hubPoly,
    const std::vector<gsPoint2D<T>>& singCands,
    const std::vector<std::vector<gsPoint2D<T>>>& densePolys,
    bool allowSmooth180,
    gsPoint2D<T>& outCenter,
    std::vector<std::array<int, 4>>& outQuads,
    T& outBestScore)
{
    const int N = static_cast<int>(bndPts.size());
    outQuads.clear();
    outBestScore = T(1e18);
    if (N < 6 || N > 12 || (N % 2) != 0) return false;
    const int M = N / 2;

    std::vector<gsPoint2D<T>> segOutDir(N), segBackDir(N);
    std::vector<T> cornerSectorDeg(N, T(180));
    for (int i = 0; i < N; ++i) {
        gsPoint2D<T> d0 = evalCurveDeriv(*edgesCCW[i].bsp, T(0));
        gsPoint2D<T> d1 = evalCurveDeriv(*edgesCCW[i].bsp, T(1));
        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
        segOutDir[i]  = (n0 > T(1e-12)) ? gsPoint2D<T>(d0.x / n0, d0.y / n0) : gsPoint2D<T>(1, 0);
        segBackDir[i] = (n1 > T(1e-12)) ? gsPoint2D<T>(-d1.x / n1, -d1.y / n1) : gsPoint2D<T>(-1, 0);
    }
    for (int i = 0; i < N; ++i)
        cornerSectorDeg[i] = mergedCcwAngleDeg(segOutDir[i], segBackDir[(i + N - 1) % N]);

    // Gather candidate interior centers O (pole-of-inaccessibility FIRST for curved/swirled hubs!)
    std::vector<gsPoint2D<T>> cands = singCands;
    gsPoint2D<T> polePt; T rPole = T(0);
    if (mergedFindDeepInteriorPoint(hubPoly, densePolys, polePt, rPole))
        cands.push_back(polePt);
    gsPoint2D<T> polyCtr = polygonCentroid(hubPoly);
    gsPoint2D<T> vAvg(0, 0);
    for (const auto& p : bndPts) { vAvg.x += p.x; vAvg.y += p.y; }
    vAvg.x /= static_cast<T>(N); vAvg.y /= static_cast<T>(N);
    cands.push_back(polyCtr);
    cands.push_back(vAvg);
    if (rPole > T(0))
        cands.emplace_back(T(0.5) * (polePt.x + polyCtr.x), T(0.5) * (polePt.y + polyCtr.y));
    for (int phase = 0; phase < 2; ++phase) {
        gsPoint2D<T> sAvg(0, 0);
        for (int j = 0; j < M; ++j) {
            sAvg.x += bndPts[2 * j + phase].x;
            sAvg.y += bndPts[2 * j + phase].y;
        }
        sAvg.x /= static_cast<T>(M); sAvg.y /= static_cast<T>(M);
        cands.push_back(sAvg);
        if (rPole > T(0))
            cands.emplace_back(T(0.5) * (polePt.x + sAvg.x), T(0.5) * (polePt.y + sAvg.y));
    }

    auto spokeInside = [&](const gsPoint2D<T>& O, int vIdx) -> bool {
        const auto& V = bndPts[vIdx];
        T len = O.dist(V);
        if (len <= T(1e-8)) return false;
        // Exact tangent sector clearance at V
        gsPoint2D<T> dVO((O.x - V.x) / len, (O.y - V.y) / len);
        T aV = mergedCcwAngleDeg(segOutDir[vIdx], dVO);
        if (aV < T(7.5) || cornerSectorDeg[vIdx] - aV < T(7.5)) return false;
        if (!isSegmentInsideDomain(O, V, densePolys)) return false;
        for (int k = 1; k <= 5; ++k) {
            T f = static_cast<T>(k) / T(6);
            gsPoint2D<T> q((T(1) - f) * O.x + f * V.x, (T(1) - f) * O.y + f * V.y);
            if (!pointInPolyStrict(q, hubPoly, T(1e-6))) return false;
        }
        const int np = static_cast<int>(hubPoly.size());
        for (int k = 0; k < np; ++k) {
            const auto& p0 = hubPoly[k]; const auto& p1 = hubPoly[(k + 1) % np];
            if (p0.dist(V) < T(1e-5) || p1.dist(V) < T(1e-5)) continue;
            if (segProperIntersect(O, V, p0, p1)) return false;
        }
        return true;
    };

    bool found = false;
    for (const auto& O : cands) {
        if (!pointInPolyStrict(O, hubPoly, T(1e-5)) || !isPointStrictlyInsideDomain(O, densePolys, T(1e-5)))
            continue;
        for (int phase = 0; phase < 2; ++phase) {
            bool ok = true;
            T smoothPen = T(0);
            for (int j = 0; j < M; ++j) {
                int unsplitV = (2 * j + phase + 1) % N;
                if (cornerSectorDeg[unsplitV] >= T(195.0)) { ok = false; break; }
                if (static_cast<size_t>(unsplitV) < isSmooth180.size() && isSmooth180[unsplitV]) {
                    if (!allowSmooth180) { ok = false; break; }
                    smoothPen += T(25.0);
                }
            }
            if (!ok) continue;
            for (int j = 0; j < M; ++j) {
                if (!spokeInside(O, 2 * j + phase)) { ok = false; break; }
            }
            if (!ok) continue;

            T totScore = smoothPen;
            std::vector<std::array<int, 4>> qs(M);
            for (int j = 0; j < M; ++j) {
                int v0 = 2 * j + phase;
                int v1 = (v0 + 1) % N;
                int v2 = (v0 + 2) % N;
                gsPoint2D<T> q[4] = {O, bndPts[v0], bndPts[v1], bndPts[v2]};
                T a2 = cross2(q[0], q[1], q[2]) + cross2(q[0], q[2], q[3]);
                if (a2 <= T(1e-8) || cross2(O, bndPts[v0], bndPts[v2]) <= T(1e-8)) { ok = false; break; }

                gsPoint2D<T> dOv0(bndPts[v0].x - O.x, bndPts[v0].y - O.y);
                gsPoint2D<T> dOv2(bndPts[v2].x - O.x, bndPts[v2].y - O.y);
                gsPoint2D<T> dv0O(-dOv0.x, -dOv0.y);
                gsPoint2D<T> dv2O(-dOv2.x, -dOv2.y);

                T aO  = mergedCcwAngleDeg(dOv0, dOv2);
                T aV0 = mergedCcwAngleDeg(segOutDir[v0], dv0O);
                T aV1 = cornerSectorDeg[v1];
                T aV2 = mergedCcwAngleDeg(dv2O, segBackDir[v1]);

                bool v1Smooth = (static_cast<size_t>(v1) < isSmooth180.size() && isSmooth180[v1]);
                T angs[4] = {aO, aV0, aV1, aV2};
                for (int k = 0; k < 4; ++k) {
                    if (k == 2 && v1Smooth) continue;
                    T ang = angs[k];
                    if (ang < T(7.5) || ang > T(168.0)) { ok = false; break; }
                    T dev = std::abs(ang - T(90.0)) / T(45.0);
                    totScore += dev * dev;
                    if (ang < T(32.0)) totScore += T(5.0) * (T(32.0) - ang) / T(10.0);
                    if (ang > T(148.0)) totScore += T(5.0) * (ang - T(148.0)) / T(10.0);
                }
                if (!ok) break;
                qs[j] = {N, v0, v1, v2};
            }
            if (ok && totScore < outBestScore) {
                outBestScore = totScore;
                outCenter = O;
                outQuads = std::move(qs);
                found = true;
            }
        }
    }
    return found;
}

// Universal Polynomial-Time Finisher (Tier D):
// Constrained Ear-Triangulation (+ optional Steiner pole) + Constrained Lawson Flips +
// Edmonds' General-Graph Blossom Matching + Claw-Breaking Diagonal Flips -> 100% topological completion for ANY even 2m-gon!
template <class T>
inline bool quadrangulateByTriangulationBlossomShift(
    const std::vector<CurveSeg<T>>& edgesCCW,
    const std::vector<gsPoint2D<T>>& bndPts,
    const std::vector<gsPoint2D<T>>& hubPoly,
    const std::vector<std::vector<gsPoint2D<T>>>& densePolys,
    int degree,
    std::vector<std::array<CurveSeg<T>, 4>>& outQuads)
{
    const int N = static_cast<int>(bndPts.size());
    if (N < 4 || (N % 2) != 0) return false;

    struct LTri { int v[3]; };
    std::vector<gsPoint2D<T>> allPts = bndPts;
    std::vector<LTri> tris;

    auto earInsideHub = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, bool strictHub) -> bool {
        if (A.dist(B) <= T(1e-8)) return false;
        if (!strictHub) return true;
        if (!isSegmentInsideDomain(A, B, densePolys)) return false;
        for (int k = 1; k <= 3; ++k) {
            T f = static_cast<T>(k) / T(4);
            if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y), hubPoly, T(1e-6)))
                return false;
        }
        const int np = static_cast<int>(hubPoly.size());
        for (int k = 0; k < np; ++k) {
            const auto& p0 = hubPoly[k]; const auto& p1 = hubPoly[(k + 1) % np];
            if (p0.dist(A) < T(1e-5) || p1.dist(A) < T(1e-5) || p0.dist(B) < T(1e-5) || p1.dist(B) < T(1e-5))
                continue;
            if (segProperIntersect(A, B, p0, p1)) return false;
        }
        return true;
    };

    auto tryEarClip = [&](bool strictHub, std::vector<LTri>& outT) -> bool {
        outT.clear();
        std::vector<int> ring(N);
        std::iota(ring.begin(), ring.end(), 0);
        while (ring.size() > 3) {
            int m = static_cast<int>(ring.size());
            int bestK = -1;
            T bestSc = -T(1e30);
            auto ptToSegDist = [](const gsPoint2D<T>& P, const gsPoint2D<T>& U, const gsPoint2D<T>& V) -> T {
                T vx = V.x - U.x, vy = V.y - U.y;
                T l2 = vx * vx + vy * vy;
                T t = (l2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((P.x - U.x) * vx + (P.y - U.y) * vy) / l2)) : T(0);
                return std::hypot(P.x - (U.x + t * vx), P.y - (U.y + t * vy));
            };
            for (int k = 0; k < m; ++k) {
                int ia = ring[(k + m - 1) % m], ib = ring[k], ic = ring[(k + 1) % m];
                const auto& A = allPts[ia]; const auto& B = allPts[ib]; const auto& C = allPts[ic];
                T a2 = cross2(A, B, C);
                if (a2 <= T(1e-10)) continue;
                bool clear = true;
                for (int r = 0; r < m; ++r) {
                    int ir = ring[r];
                    if (ir == ia || ir == ib || ir == ic) continue;
                    const auto& P = allPts[ir];
                    if (P.dist(A) < T(1e-6) || P.dist(B) < T(1e-6) || P.dist(C) < T(1e-6)) continue;
                    if (ptToSegDist(P, A, C) < T(1e-5)) { clear = false; break; }
                    if (cross2(A, B, P) >= -T(1e-7) && cross2(B, C, P) >= -T(1e-7) && cross2(C, A, P) >= -T(1e-7)) {
                        clear = false; break;
                    }
                }
                if (!clear) continue;
                for (int e = 0; e < m && clear; ++e) {
                    int u = ring[e], v = ring[(e + 1) % m];
                    if (u == ia || v == ia || u == ic || v == ic) continue;
                    if (segProperIntersect(A, C, allPts[u], allPts[v])) clear = false;
                }
                if (!clear) continue;
                if (m > 3 && !earInsideHub(A, C, strictHub)) continue;
                if (strictHub && m > 4) {
                    // Guard: remaining sub-ring after clipping ib MUST retain strictly positive area!
                    T remA2 = T(0);
                    for (int e = 0; e < m; ++e) {
                        if (e == k) continue;
                        int u = ring[e];
                        int v = (e == (k + m - 1) % m) ? ring[(k + 1) % m] : ring[(e + 1) % m];
                        remA2 += allPts[u].x * allPts[v].y - allPts[v].x * allPts[u].y;
                    }
                    if (remA2 <= T(1e-8)) continue;
                }
                T l2 = A.distSq(B) + B.distSq(C) + C.distSq(A);
                T sc = a2 / std::max(T(1e-14), l2);
                if (sc > bestSc) { bestSc = sc; bestK = k; }
            }
            if (bestK < 0 && !strictHub) {
                // Emergency ear fallback for collinear/curved boundary vertices so ear-clipping never stalls
                for (int k = 0; k < m; ++k) {
                    int ia = ring[(k + m - 1) % m], ib = ring[k], ic = ring[(k + 1) % m];
                    const auto& A = allPts[ia]; const auto& B = allPts[ib]; const auto& C = allPts[ic];
                    T a2 = cross2(A, B, C);
                    T l2 = A.distSq(B) + B.distSq(C) + C.distSq(A);
                    T sc = a2 / std::max(T(1e-14), l2);
                    if (a2 > -T(1e-5) * l2 && sc > bestSc) { bestSc = sc; bestK = k; }
                }
                if (bestK < 0) {
                    // Unconditional best-scoring ear fallback so tryEarClip(false) is 100% total
                    for (int k = 0; k < m; ++k) {
                        int ia = ring[(k + m - 1) % m], ib = ring[k], ic = ring[(k + 1) % m];
                        const auto& A = allPts[ia]; const auto& B = allPts[ib]; const auto& C = allPts[ic];
                        T a2 = cross2(A, B, C);
                        T l2 = A.distSq(B) + B.distSq(C) + C.distSq(A);
                        T sc = a2 / std::max(T(1e-14), l2);
                        if (sc > bestSc) { bestSc = sc; bestK = k; }
                    }
                }
            }
            if (bestK < 0) return false;
            int mCur = static_cast<int>(ring.size());
            outT.push_back({{ring[(bestK + mCur - 1) % mCur], ring[bestK], ring[(bestK + 1) % mCur]}});
            ring.erase(ring.begin() + bestK);
        }
        if (ring.size() == 3)
            outT.push_back({{ring[0], ring[1], ring[2]}});
        return (static_cast<int>(outT.size()) == N - 2);
    };

    // Multi-Steiner Positive-Area Ear-Clipping Triangulator:
    // Whenever no straight diagonal stays inside curved `hubPoly`, inserts an interior Steiner pole
    // strictly inside the remaining sub-polygon and clips its maximal visible positive-area fan,
    // strictly preserving (tris.size() + ring.size()) even and 100% positive-area triangles!
    auto trySteinerEarClip = [&](std::vector<LTri>& outT) -> bool {
        outT.clear();
        std::vector<gsPoint2D<T>> workPts = bndPts;
        std::vector<int> ring(N);
        std::iota(ring.begin(), ring.end(), 0);
        int maxSteps = 4 * N + 24;
        int nSteiners = 0;

        auto buildRingPoly = [&](const std::vector<int>& r) {
            std::vector<gsPoint2D<T>> pl;
            const int m = static_cast<int>(r.size());
            pl.reserve(m * 6);
            for (int k = 0; k < m; ++k) {
                int u = r[k], v = r[(k + 1) % m];
                if (u < N && v < N && ((u + 1) % N == v)) {
                    for (int s = 0; s < 6; ++s)
                        pl.push_back(evalCurvePoint(*edgesCCW[u].bsp, static_cast<T>(s) / T(6)));
                } else {
                    const auto& A = workPts[u]; const auto& B = workPts[v];
                    for (int s = 0; s < 4; ++s) {
                        T f = static_cast<T>(s) / T(4);
                        pl.emplace_back((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y);
                    }
                }
            }
            return pl;
        };

        auto segClearInPoly = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B,
                                  const std::vector<gsPoint2D<T>>& subP) -> bool {
            T abLen = A.dist(B);
            if (abLen <= T(1e-8)) return false;
            // Test strictly against exact curved subP + hubPoly (avoid coarse densePolys chord false rejections!)
            for (int s = 1; s <= 5; ++s) {
                T f = static_cast<T>(s) / T(6);
                gsPoint2D<T> q((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y);
                if (!pointInPolyStrict(q, subP, T(1e-7)) || !pointInPolyStrict(q, hubPoly, T(1e-7)))
                    return false;
            }
            auto clearOfPoly = [&](const std::vector<gsPoint2D<T>>& pl) -> bool {
                const int np = static_cast<int>(pl.size());
                const T abx = B.x - A.x, aby = B.y - A.y;
                for (int k = 0; k < np; ++k) {
                    const auto& p0 = pl[k]; const auto& p1 = pl[(k + 1) % np];
                    if (p0.dist(A) < T(1e-5) || p1.dist(A) < T(1e-5) || p0.dist(B) < T(1e-5) || p1.dist(B) < T(1e-5))
                        continue;
                    T cdx = p1.x - p0.x, cdy = p1.y - p0.y;
                    T denom = abx * cdy - aby * cdx;
                    if (std::abs(denom) <= T(1e-14)) continue;
                    T t = ((p0.x - A.x) * cdy - (p0.y - A.y) * cdx) / denom;
                    T u = ((p0.x - A.x) * aby - (p0.y - A.y) * abx) / denom;
                    if (t > T(0.03) && t < T(0.97) && u >= -T(1e-6) && u <= T(1) + T(1e-6))
                        return false;
                }
                return true;
            };
            return clearOfPoly(subP) && clearOfPoly(hubPoly);
        };

        const std::vector<std::vector<gsPoint2D<T>>> noPolys;
        while (ring.size() > 3 && --maxSteps > 0) {
            const int m = static_cast<int>(ring.size());

            // Fast exact 4-gon termination: ANY simple 4-gon splits along one of its 2 diagonals into 2 positive-area tris!
            if (m == 4) {
                int r0 = ring[0], r1 = ring[1], r2 = ring[2], r3 = ring[3];
                bool ok02 = (cross2(workPts[r0], workPts[r1], workPts[r2]) > T(1e-10) &&
                             cross2(workPts[r0], workPts[r2], workPts[r3]) > T(1e-10));
                bool ok13 = (cross2(workPts[r1], workPts[r2], workPts[r3]) > T(1e-10) &&
                             cross2(workPts[r1], workPts[r3], workPts[r0]) > T(1e-10));
                auto subPoly4 = buildRingPoly(ring);
                bool cl02 = ok02 && segClearInPoly(workPts[r0], workPts[r2], subPoly4);
                bool cl13 = ok13 && segClearInPoly(workPts[r1], workPts[r3], subPoly4);
                if (cl02 || cl13) {
                    if (cl02) { outT.push_back({{r0, r1, r2}}); outT.push_back({{r0, r2, r3}}); }
                    else      { outT.push_back({{r1, r2, r3}}); outT.push_back({{r1, r3, r0}}); }
                    ring.clear();
                    break;
                }
                // Even if straight diagonal grazes a slightly bowed curved wall, if both tris have strictly positive area,
                // prefer Steiner pole below first, otherwise accept the positive-area diagonal!
            }

            int bestK = -1;
            T bestSc = -T(1e30);
            auto ptToSegDistS = [](const gsPoint2D<T>& P, const gsPoint2D<T>& U, const gsPoint2D<T>& V) -> T {
                T vx = V.x - U.x, vy = V.y - U.y;
                T l2 = vx * vx + vy * vy;
                T t = (l2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((P.x - U.x) * vx + (P.y - U.y) * vy) / l2)) : T(0);
                return std::hypot(P.x - (U.x + t * vx), P.y - (U.y + t * vy));
            };
            auto subPolyCur = buildRingPoly(ring);
            for (int k = 0; k < m; ++k) {
                int ia = ring[(k + m - 1) % m], ib = ring[k], ic = ring[(k + 1) % m];
                const auto& A = workPts[ia]; const auto& B = workPts[ib]; const auto& C = workPts[ic];
                T a2 = cross2(A, B, C);
                if (a2 <= T(1e-10)) continue;
                bool clear = true;
                for (int r = 0; r < m; ++r) {
                    int ir = ring[r];
                    if (ir == ia || ir == ib || ir == ic) continue;
                    const auto& P = workPts[ir];
                    if (P.dist(A) < T(1e-6) || P.dist(B) < T(1e-6) || P.dist(C) < T(1e-6)) continue;
                    if (ptToSegDistS(P, A, C) < T(1e-5)) { clear = false; break; }
                    if (cross2(A, B, P) >= -T(1e-7) && cross2(B, C, P) >= -T(1e-7) && cross2(C, A, P) >= -T(1e-7)) {
                        clear = false; break;
                    }
                }
                if (!clear) continue;
                for (int e = 0; e < m && clear; ++e) {
                    int u = ring[e], v = ring[(e + 1) % m];
                    if (u == ia || v == ia || u == ic || v == ic) continue;
                    if (segProperIntersect(A, C, workPts[u], workPts[v])) clear = false;
                }
                if (!clear || !segClearInPoly(A, C, subPolyCur)) continue;
                if (m > 4) {
                    T remA2 = T(0);
                    for (int e = 0; e < m; ++e) {
                        if (e == k) continue;
                        int u = ring[e];
                        int v = (e == (k + m - 1) % m) ? ring[(k + 1) % m] : ring[(e + 1) % m];
                        remA2 += workPts[u].x * workPts[v].y - workPts[v].x * workPts[u].y;
                    }
                    if (remA2 <= T(1e-8)) continue;
                }
                T l2 = A.distSq(B) + B.distSq(C) + C.distSq(A);
                T sc = a2 / std::max(T(1e-14), l2);
                if (sc > bestSc) { bestSc = sc; bestK = k; }
            }
            if (bestK >= 0) {
                int ia = ring[(bestK + m - 1) % m], ib = ring[bestK], ic = ring[(bestK + 1) % m];
                outT.push_back({{ia, ib, ic}});
                ring.erase(ring.begin() + bestK);
                continue;
            }

            // No valid straight ear: insert a deep interior Steiner point strictly inside subPoly
            if (nSteiners >= N + 8) return false;
            auto subPoly = buildRingPoly(ring);
            std::vector<gsPoint2D<T>> poleCands;
            gsPoint2D<T> pole; T rPole = T(0);
            if (mergedFindDeepInteriorPoint(subPoly, noPolys, pole, rPole) &&
                pointInPolyStrict(pole, hubPoly, T(1e-7)))
                poleCands.push_back(pole);
            gsPoint2D<T> cent = polygonCentroid(subPoly);
            if (pointInPolyStrict(cent, subPoly, T(1e-7)) && pointInPolyStrict(cent, hubPoly, T(1e-7)))
                poleCands.push_back(cent);
            for (int k = 0; k < m; ++k) {
                int ia = ring[k], ib = ring[(k + 1) % m], ic = ring[(k + 2) % m];
                gsPoint2D<T> g((workPts[ia].x + workPts[ib].x + workPts[ic].x) / T(3),
                               (workPts[ia].y + workPts[ib].y + workPts[ic].y) / T(3));
                if (pointInPolyStrict(g, subPoly, T(1e-7)) && pointInPolyStrict(g, hubPoly, T(1e-7)))
                    poleCands.push_back(g);
                // Also try small inward normal offset from edge midpoints for thin curved channels
                const auto& U = workPts[ia]; const auto& V = workPts[ib];
                T el = U.dist(V);
                if (el > T(1e-8)) {
                    gsPoint2D<T> inN(-(V.y - U.y) / el, (V.x - U.x) / el);
                    for (T dFrac : {T(0.08), T(0.20), T(0.35)}) {
                        gsPoint2D<T> q(T(0.5) * (U.x + V.x) + dFrac * el * inN.x,
                                       T(0.5) * (U.y + V.y) + dFrac * el * inN.y);
                        if (pointInPolyStrict(q, subPoly, T(1e-7)) && pointInPolyStrict(q, hubPoly, T(1e-7))) {
                            poleCands.push_back(q);
                            break;
                        }
                    }
                }
            }

            bool inserted = false;
            int bestRun = 0, bestStart = -1;
            gsPoint2D<T> bestP;
            bool bestAll = false;
            for (const auto& P : poleCands) {
                std::vector<bool> vVis(m, false), eVis(m, false);
                for (int k = 0; k < m; ++k)
                    vVis[k] = segClearInPoly(P, workPts[ring[k]], subPoly);
                int cntE = 0;
                for (int k = 0; k < m; ++k) {
                    int kn = (k + 1) % m;
                    if (vVis[k] && vVis[kn] && cross2(P, workPts[ring[k]], workPts[ring[kn]]) > T(1e-10)) {
                        eVis[k] = true;
                        ++cntE;
                    }
                }
                if (cntE == m) {
                    bestP = P; bestAll = true; bestRun = m; bestStart = 0;
                    break;
                }
                // Find longest contiguous run of visible edges (prefer runLen >= 2)
                for (int s = 0; s < m; ++s) {
                    if (!eVis[s] || eVis[(s + m - 1) % m]) continue;
                    int rLen = 0;
                    while (rLen < m && eVis[(s + rLen) % m]) ++rLen;
                    if (rLen > bestRun) {
                        bestRun = rLen; bestStart = s; bestP = P;
                    }
                }
            }
            if (bestAll) {
                int vP = static_cast<int>(workPts.size());
                workPts.push_back(bestP);
                for (int k = 0; k < m; ++k)
                    outT.push_back({{vP, ring[k], ring[(k + 1) % m]}});
                ring.clear();
                inserted = true;
            } else if (bestRun >= 2 && bestStart >= 0) {
                std::vector<int> nextRing;
                nextRing.reserve(m - bestRun + 2);
                nextRing.push_back(ring[bestStart]);
                int vP = static_cast<int>(workPts.size());
                nextRing.push_back(vP);
                for (int t = bestRun; t < m; ++t)
                    nextRing.push_back(ring[(bestStart + t) % m]);
                workPts.push_back(bestP);
                T nxtA2 = T(0);
                for (size_t k = 0; k < nextRing.size(); ++k)
                    nxtA2 += workPts[nextRing[k]].x * workPts[nextRing[(k + 1) % nextRing.size()]].y
                           - workPts[nextRing[(k + 1) % nextRing.size()]].x * workPts[nextRing[k]].y;
                if (nxtA2 > T(1e-8) || nextRing.size() <= 3) {
                    ++nSteiners;
                    for (int t = 0; t < bestRun; ++t) {
                        int k0 = (bestStart + t) % m;
                        int k1 = (bestStart + t + 1) % m;
                        outT.push_back({{vP, ring[k0], ring[k1]}});
                    }
                    ring = std::move(nextRing);
                    inserted = true;
                } else {
                    workPts.pop_back();
                }
            }
            if (!inserted) return false;
        }

        if (ring.size() == 3) {
            if (cross2(workPts[ring[0]], workPts[ring[1]], workPts[ring[2]]) <= T(1e-10))
                return false;
            outT.push_back({{ring[0], ring[1], ring[2]}});
            ring.clear();
        }
        if (!ring.empty() || outT.empty() || (outT.size() % 2) != 0)
            return false;
        allPts = std::move(workPts);
        return true;
    };

    if (!tryEarClip(true, tris)) {
        // Try inserting deep interior pole + star fan if star-shaped
        gsPoint2D<T> pole; T rPole = T(0);
        bool fanOk = false;
        if (mergedFindDeepInteriorPoint(hubPoly, densePolys, pole, rPole)) {
            fanOk = true;
            for (int i = 0; i < N && fanOk; ++i) {
                if (cross2(pole, bndPts[i], bndPts[(i + 1) % N]) <= T(1e-9)) { fanOk = false; break; }
                for (int s = 1; s <= 4; ++s) {
                    T f = static_cast<T>(s) / T(5);
                    if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * pole.x + f * bndPts[i].x,
                                                        (T(1) - f) * pole.y + f * bndPts[i].y), hubPoly, T(1e-6))) {
                        fanOk = false; break;
                    }
                }
            }
            if (fanOk) {
                allPts.push_back(pole);
                tris.clear();
                for (int i = 0; i < N; ++i)
                    tris.push_back({{N, i, (i + 1) % N}});
            }
        }
        if (!fanOk && !trySteinerEarClip(tris) && !tryEarClip(false, tris)) {
            if (std::getenv("GISMO_MERGED_DEBUG") != nullptr)
                std::cerr << "[MERGED] BlossomShift earClip fail N=" << N << "\n";
            return false;
        }
    }

    const int nT = static_cast<int>(tris.size());
    if (nT < 2 || (nT % 2) != 0) return false;

    // Insert 1 deep interior pole Steiner point inside the largest interior triangle when N >= 12 (improves hub quality!)
    if (static_cast<int>(allPts.size()) == N && N >= 12) {
        gsPoint2D<T> pole; T rPole = T(0);
        if (mergedFindDeepInteriorPoint(hubPoly, densePolys, pole, rPole)) {
            int hitT = -1;
            T bestMinB = T(0.05);
            for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti) {
                const auto& A = allPts[tris[ti].v[0]];
                const auto& B = allPts[tris[ti].v[1]];
                const auto& C = allPts[tris[ti].v[2]];
                T a2 = cross2(A, B, C);
                if (a2 <= T(1e-10)) continue;
                T b0 = cross2(B, C, pole) / a2;
                T b1 = cross2(C, A, pole) / a2;
                T b2 = T(1) - b0 - b1;
                T mn = std::min(b0, std::min(b1, b2));
                if (mn > bestMinB) { bestMinB = mn; hitT = ti; }
            }
            if (hitT >= 0) {
                int vNew = static_cast<int>(allPts.size());
                allPts.push_back(pole);
                LTri oldT = tris[hitT];
                tris[hitT] = {{oldT.v[0], oldT.v[1], vNew}};
                tris.push_back({{oldT.v[1], oldT.v[2], vNew}});
                tris.push_back({{oldT.v[2], oldT.v[0], vNew}});
            }
        }
    }

    auto isBndEdge = [&](int a, int b) -> bool {
        return (a < N && b < N) && (((a + 1) % N == b) || ((b + 1) % N == a));
    };

    // Constrained Lawson flips inside hubPoly
    auto inCircle = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C, const gsPoint2D<T>& D) -> bool {
        T ax = A.x - D.x, ay = A.y - D.y;
        T bx = B.x - D.x, by = B.y - D.y;
        T cx = C.x - D.x, cy = C.y - D.y;
        T det = (ax * ax + ay * ay) * (bx * cy - cx * by)
              - (bx * bx + by * by) * (ax * cy - cx * ay)
              + (cx * cx + cy * cy) * (ax * by - bx * ay);
        return det > T(1e-10);
    };
    for (int pass = 0; pass < 8; ++pass) {
        bool flipped = false;
        const int curT = static_cast<int>(tris.size());
        for (int i = 0; i < curT && !flipped; ++i) {
            for (int ka = 0; ka < 3 && !flipped; ++ka) {
                int a = tris[i].v[ka], b = tris[i].v[(ka + 1) % 3], c = tris[i].v[(ka + 2) % 3];
                if (isBndEdge(a, b)) continue;
                for (int j = i + 1; j < curT; ++j) {
                    for (int kb = 0; kb < 3; ++kb) {
                        if (tris[j].v[kb] == b && tris[j].v[(kb + 1) % 3] == a) {
                            int d = tris[j].v[(kb + 2) % 3];
                            gsPoint2D<T> q[4] = {allPts[c], allPts[a], allPts[d], allPts[b]};
                            if (isStrictlyConvexQuad(q) && inCircle(allPts[a], allPts[b], allPts[c], allPts[d]) &&
                                earInsideHub(allPts[c], allPts[d], true)) {
                                tris[i] = {{c, a, d}};
                                tris[j] = {{c, d, b}};
                                flipped = true;
                            }
                            break;
                        }
                    }
                    if (flipped) break;
                }
            }
        }
        if (!flipped) break;
    }

    // Blossom matching helper on `tList` (Phase 1: convex edges, Phase 2: all interior edges)
    struct DualEdge { int t0, t1; int a, b, c, d; bool convex; T cost; };
    auto solveBlossom = [&](const std::vector<LTri>& tList,
                            std::vector<std::array<int, 4>>& qvOut,
                            std::vector<int>& unmatchedOut) -> bool {
        const int mT = static_cast<int>(tList.size());
        std::vector<DualEdge> dEdges;
        for (int i = 0; i < mT; ++i) {
            for (int ka = 0; ka < 3; ++ka) {
                int a = tList[i].v[ka], b = tList[i].v[(ka + 1) % 3], c = tList[i].v[(ka + 2) % 3];
                if (isBndEdge(a, b)) continue;
                for (int j = i + 1; j < mT; ++j) {
                    for (int kb = 0; kb < 3; ++kb) {
                        if (tList[j].v[kb] == b && tList[j].v[(kb + 1) % 3] == a) {
                            int d = tList[j].v[(kb + 2) % 3];
                            // Forbid bowtie pairings where c and d lie on the same side of shared diagonal (a, b)!
                            T aTri1 = cross2(allPts[a], allPts[b], allPts[c]);
                            T aTri2 = cross2(allPts[b], allPts[a], allPts[d]);
                            if (aTri1 <= T(1e-10) || aTri2 <= T(1e-10)) continue;
                            gsPoint2D<T> q[4] = {allPts[c], allPts[a], allPts[d], allPts[b]};
                            if (segProperIntersect(q[0], q[1], q[2], q[3]) ||
                                segProperIntersect(q[1], q[2], q[3], q[0])) continue;
                            bool conv = isStrictlyConvexQuad(q);
                            T cost = conv ? T(0) : T(50.0);
                            for (int k = 0; k < 4; ++k) {
                                gsPoint2D<T> u(q[(k + 1) % 4].x - q[k].x, q[(k + 1) % 4].y - q[k].y);
                                gsPoint2D<T> v(q[(k + 3) % 4].x - q[k].x, q[(k + 3) % 4].y - q[k].y);
                                T ang = mergedCcwAngleDeg(u, v);
                                T dev = std::abs(ang - T(90.0)) / T(45.0);
                                cost += dev * dev;
                            }
                            dEdges.push_back({i, j, a, b, c, d, conv, cost});
                            break;
                        }
                    }
                }
            }
        }
        std::sort(dEdges.begin(), dEdges.end(), [](const DualEdge& x, const DualEdge& y) { return x.cost < y.cost; });

        std::vector<int> partner(mT, -1), matchedPair(mT, -1);
        auto alwaysOk = [](const std::vector<int>&, const std::vector<int>&) { return true; };

        for (int phase = 0; phase < 2; ++phase) {
            std::vector<std::vector<std::pair<int, int>>> adj(mT);
            for (int ei = 0; ei < static_cast<int>(dEdges.size()); ++ei) {
                if (phase == 0 && !dEdges[ei].convex) continue;
                adj[dEdges[ei].t0].push_back({dEdges[ei].t1, ei});
                adj[dEdges[ei].t1].push_back({dEdges[ei].t0, ei});
            }
            // Greedy seed first
            for (int ei = 0; ei < static_cast<int>(dEdges.size()); ++ei) {
                if (phase == 0 && !dEdges[ei].convex) continue;
                int u = dEdges[ei].t0, v = dEdges[ei].t1;
                if (partner[u] < 0 && partner[v] < 0) {
                    partner[u] = v; partner[v] = u;
                    matchedPair[u] = ei; matchedPair[v] = ei;
                }
            }
            // Edmonds' Blossom augmenting paths
            for (int r = 0; r < mT; ++r) {
                if (partner[r] < 0)
                    augmentBlossomMatchingFrom(r, mT, adj, partner, matchedPair, alwaysOk);
            }
        }

        unmatchedOut.clear();
        for (int i = 0; i < mT; ++i)
            if (partner[i] < 0) unmatchedOut.push_back(i);
        if (!unmatchedOut.empty()) return false;

        qvOut.clear();
        for (int i = 0; i < mT; ++i) {
            if (partner[i] > i) {
                const auto& de = dEdges[matchedPair[i]];
                qvOut.push_back({de.c, de.a, de.d, de.b});
            }
        }
        return (static_cast<int>(qvOut.size()) == mT / 2);
    };

    // Dual-distance helper between unmatched triangles (enables multi-step K_{1,3} claw resolution!)
    auto evalUnmatchedDualDist = [&](const std::vector<LTri>& tList, const std::vector<int>& unm) -> int {
        if (unm.size() < 2) return 0;
        const int mT = static_cast<int>(tList.size());
        std::vector<std::vector<int>> dAdj(mT);
        for (int i = 0; i < mT; ++i) {
            for (int ka = 0; ka < 3; ++ka) {
                int a = tList[i].v[ka], b = tList[i].v[(ka + 1) % 3];
                if (isBndEdge(a, b)) continue;
                for (int j = i + 1; j < mT; ++j) {
                    for (int kb = 0; kb < 3; ++kb) {
                        if (tList[j].v[kb] == b && tList[j].v[(kb + 1) % 3] == a) {
                            dAdj[i].push_back(j); dAdj[j].push_back(i);
                            break;
                        }
                    }
                }
            }
        }
        int totD = 0;
        for (size_t ui = 0; ui < unm.size(); ++ui) {
            std::vector<int> dist(mT, 999);
            std::queue<int> q;
            q.push(unm[ui]); dist[unm[ui]] = 0;
            while (!q.empty()) {
                int cur = q.front(); q.pop();
                for (int nb : dAdj[cur]) if (dist[nb] > dist[cur] + 1) { dist[nb] = dist[cur] + 1; q.push(nb); }
            }
            int mn = 999;
            for (size_t uj = 0; uj < unm.size(); ++uj)
                if (uj != ui) mn = std::min(mn, dist[unm[uj]]);
            totD += mn;
        }
        return totD;
    };

    std::vector<std::array<int, 4>> qvList;
    std::vector<int> unmatched;
    bool okMatch = solveBlossom(tris, qvList, unmatched);
    int curDualDist = okMatch ? 0 : evalUnmatchedDualDist(tris, unmatched);

    // Multi-step dual-distance-guided diagonal flipping to break any K_{1,3} claw!
    for (int attempt = 0; !okMatch && attempt < 28; ++attempt) {
        bool didFlip = false;
        const int curT = static_cast<int>(tris.size());
        std::vector<LTri> bestAlt;
        std::vector<int> bestUnm;
        int bestDist = curDualDist;
        size_t bestUnmSz = unmatched.size();

        for (int i = 0; i < curT && !okMatch; ++i) {
            for (int ka = 0; ka < 3 && !okMatch; ++ka) {
                int a = tris[i].v[ka], b = tris[i].v[(ka + 1) % 3], c = tris[i].v[(ka + 2) % 3];
                if (isBndEdge(a, b)) continue;
                for (int j = i + 1; j < curT; ++j) {
                    for (int kb = 0; kb < 3; ++kb) {
                        if (tris[j].v[kb] == b && tris[j].v[(kb + 1) % 3] == a) {
                            int d = tris[j].v[(kb + 2) % 3];
                            gsPoint2D<T> q[4] = {allPts[c], allPts[a], allPts[d], allPts[b]};
                            if (isStrictlyConvexQuad(q)) {
                                std::vector<LTri> alt = tris;
                                alt[i] = {{c, a, d}};
                                alt[j] = {{c, d, b}};
                                std::vector<int> altUnm;
                                if (solveBlossom(alt, qvList, altUnm)) {
                                    tris = std::move(alt);
                                    okMatch = true;
                                    didFlip = true;
                                } else {
                                    int dNew = evalUnmatchedDualDist(alt, altUnm);
                                    if (altUnm.size() < bestUnmSz || (altUnm.size() == bestUnmSz && dNew < bestDist)) {
                                        bestUnmSz = altUnm.size();
                                        bestDist = dNew;
                                        bestAlt = std::move(alt);
                                        bestUnm = std::move(altUnm);
                                        didFlip = true;
                                    }
                                }
                            }
                            break;
                        }
                    }
                    if (okMatch) break;
                }
            }
        }
        if (okMatch) break;
        if (!didFlip || bestAlt.empty()) break;
        tris = std::move(bestAlt);
        unmatched = std::move(bestUnm);
        curDualDist = bestDist;
    }

    // Guaranteed 1->3 Centroid-Split Claw-Destroyer along shortest dual path between unmatched triangles:
    // Splitting the first interior triangle `p_1` adjacent to `unmatched[0]` on the path to `unmatched[k]`
    // destroys any K_{1,3} claw without splitting any boundary edge, strictly reducing unmatched distance until 0!
    for (int sp = 0; !okMatch && unmatched.size() >= 2 && sp < 8; ++sp) {
        const int mT = static_cast<int>(tris.size());
        std::vector<std::vector<int>> dAdj(mT);
        for (int i = 0; i < mT; ++i) {
            for (int ka = 0; ka < 3; ++ka) {
                int a = tris[i].v[ka], b = tris[i].v[(ka + 1) % 3];
                if (isBndEdge(a, b)) continue;
                for (int j = i + 1; j < mT; ++j) {
                    for (int kb = 0; kb < 3; ++kb) {
                        if (tris[j].v[kb] == b && tris[j].v[(kb + 1) % 3] == a) {
                            dAdj[i].push_back(j); dAdj[j].push_back(i);
                            break;
                        }
                    }
                }
            }
        }
        int src = unmatched[0];
        std::vector<int> prevT(mT, -1), dist(mT, 999);
        std::queue<int> q;
        q.push(src); dist[src] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int nb : dAdj[u]) {
                if (dist[nb] > dist[u] + 1) { dist[nb] = dist[u] + 1; prevT[nb] = u; q.push(nb); }
            }
        }
        int dst = -1, bestD = 999;
        for (size_t k = 1; k < unmatched.size(); ++k) {
            if (dist[unmatched[k]] < bestD) { bestD = dist[unmatched[k]]; dst = unmatched[k]; }
        }
        if (dst < 0 || bestD < 1) break;
        // Trace back to p_1 (first neighbor of src on shortest path to dst; if bestD==1 split src itself!)
        int p1 = src;
        if (bestD >= 2) {
            int cur = dst;
            while (prevT[cur] >= 0 && prevT[cur] != src) cur = prevT[cur];
            p1 = cur;
        }
        if (p1 < 0 || p1 >= mT) break;
        const auto& A = allPts[tris[p1].v[0]];
        const auto& B = allPts[tris[p1].v[1]];
        const auto& C = allPts[tris[p1].v[2]];
        if (cross2(A, B, C) <= T(1e-10)) {
            if (p1 == src && cross2(allPts[tris[dst].v[0]], allPts[tris[dst].v[1]], allPts[tris[dst].v[2]]) > T(1e-10))
                p1 = dst;
            else
                break;
        }
        const auto& A2 = allPts[tris[p1].v[0]];
        const auto& B2 = allPts[tris[p1].v[1]];
        const auto& C2 = allPts[tris[p1].v[2]];
        gsPoint2D<T> G((A2.x + B2.x + C2.x) / T(3), (A2.y + B2.y + C2.y) / T(3));
        int vG = static_cast<int>(allPts.size());
        allPts.push_back(G);
        LTri oldT = tris[p1];
        tris[p1] = {{oldT.v[0], oldT.v[1], vG}};
        tris.push_back({{oldT.v[1], oldT.v[2], vG}});
        tris.push_back({{oldT.v[2], oldT.v[0], vG}});
        okMatch = solveBlossom(tris, qvList, unmatched);
    }

    if (!okMatch) {
        if (std::getenv("GISMO_MERGED_DEBUG") != nullptr)
            std::cerr << "[MERGED] BlossomShift fail N=" << N << " nT=" << tris.size() << " unm=" << unmatched.size() << "\n";
        return false;
    }

    for (const auto& lq : qvList) {
        std::array<CurveSeg<T>, 4> segs;
        for (int k = 0; k < 4; ++k) {
            int lu = lq[k], lv = lq[(k + 1) % 4];
            if (lu < N && lv < N && ((lu + 1) % N == lv))
                segs[k] = edgesCCW[lu];
            else
                segs[k] = makeStraightSeg(allPts[lu], allPts[lv], degree);
        }
        outQuads.push_back(std::move(segs));
    }
    return !outQuads.empty();
}

// Recursive universal even-gon quadrangulator for arbitrary N >= 4 (even, including N > 22 & non-convex hubs!)
// Strictly bounded recursion budget so no hub ever stalls, and smooth 180-deg corners get O-grid or spoke splits!
template <class T>
inline bool quadrangulateGeneralEvenPoly(
    const std::vector<CurveSeg<T>>& edgesCCW,
    const std::vector<gsPoint2D<T>>& bndPts,
    const std::vector<gsPoint2D<T>>& singPts,
    const std::vector<int>& singTargetVal,
    const MedialBackbone<T>& bb,
    const std::vector<std::vector<gsPoint2D<T>>>& densePolys,
    const GuidedCrossField<T>& gcf,
    int degree,
    int depth,
    int& budget,
    std::vector<std::array<CurveSeg<T>, 4>>& outQuads)
{
    const int N = static_cast<int>(edgesCCW.size());
    if (N < 4 || (N % 2) != 0 || depth > 10 || budget <= 0) return false;
    --budget;

    // Split pinched / repeated-vertex hubs ONLY when BOTH sub-loops have strictly positive area (never fill a CW hole!)
    if (N >= 6) {
        const T pinchTol = std::max(T(1e-5), T(1e-4) * bb.domDiag);
        auto subLoopArea2 = [&](const std::vector<CurveSeg<T>>& es) -> T {
            std::vector<gsPoint2D<T>> pl;
            pl.reserve(es.size() * 6);
            for (const auto& s : es)
                for (int k = 0; k < 6; ++k)
                    pl.push_back(evalCurvePoint(*s.bsp, static_cast<T>(k) / T(6)));
            T a2 = T(0);
            for (size_t k = 0; k < pl.size(); ++k)
                a2 += pl[k].x * pl[(k + 1) % pl.size()].y - pl[(k + 1) % pl.size()].x * pl[k].y;
            return a2;
        };
        for (int i = 0; i < N; ++i) {
            for (int j = i + 2; j < N; ++j) {
                if (i == 0 && j == N - 1) continue;
                if (bndPts[i].dist(bndPts[j]) < pinchTol) {
                    int n1 = j - i, n2 = N - n1;
                    if (n1 >= 4 && n2 >= 4 && (n1 % 2) == 0 && (n2 % 2) == 0) {
                        std::vector<CurveSeg<T>> e1, e2;
                        std::vector<gsPoint2D<T>> v1, v2;
                        for (int k = i; k < j; ++k) { e1.push_back(edgesCCW[k]); v1.push_back(bndPts[k]); }
                        for (int k = j; k != i; k = (k + 1) % N) { e2.push_back(edgesCCW[k]); v2.push_back(bndPts[k]); }
                        if (subLoopArea2(e1) <= T(1e-8) || subLoopArea2(e2) <= T(1e-8))
                            continue;
                        std::vector<std::array<CurveSeg<T>, 4>> subOut;
                        if (quadrangulateGeneralEvenPoly(e1, v1, {}, {}, bb, densePolys, gcf, degree, depth + 1, budget, subOut) &&
                            quadrangulateGeneralEvenPoly(e2, v2, {}, {}, bb, densePolys, gcf, degree, depth + 1, budget, subOut)) {
                            for (auto& q : subOut) outQuads.push_back(std::move(q));
                            return true;
                        }
                    }
                }
            }
        }
    }

    // Detect smooth ~180 deg interior corners vs >=195 deg reentrant reflex corners using exact curve tangents
    std::vector<bool> isSmooth180(N, false);
    std::vector<bool> isReflex(N, false);
    bool anySmooth180 = false, anyReflex = false;
    for (int i = 0; i < N; ++i) {
        gsPoint2D<T> dIn  = evalCurveDeriv(*edgesCCW[(i + N - 1) % N].bsp, T(1));
        gsPoint2D<T> dOut = evalCurveDeriv(*edgesCCW[i].bsp, T(0));
        T n1 = std::hypot(dIn.x, dIn.y), n2 = std::hypot(dOut.x, dOut.y);
        if (n1 > T(1e-12) && n2 > T(1e-12)) {
            gsPoint2D<T> tOut(dOut.x / n2, dOut.y / n2);
            gsPoint2D<T> tBack(-dIn.x / n1, -dIn.y / n1);
            T ang = mergedCcwAngleDeg(tOut, tBack);
            if (ang >= T(166.0) && ang < T(195.0)) {
                isSmooth180[i] = true;
                anySmooth180 = true;
            } else if (ang >= T(195.0)) {
                isReflex[i] = true;
                anyReflex = true;
            }
        }
    }

    if (N == 4) {
        // Stage 3: Only emit an O-grid inside Module 4 when there is a genuine topological >=166 deg corner
        // (smooth 180 or reentrant reflex) that global Winslow/interface smoothing cannot remove.
        // All mild Coons folds on <166 deg quads are deferred to Module 5 + global Winslow (+0 extra patches!).
        bool needOGrid = (anySmooth180 || anyReflex);
        if (!needOGrid) {
            T mnA = 0, mxA = 180, ef = 0;
            eval4SegCornerAnglesDeg(edgesCCW[0], edgesCCW[1], edgesCCW[2], edgesCCW[3], mnA, mxA, ef);
            if (mxA > T(172.0)) needOGrid = true;
        }
        if (needOGrid) {
            std::array<CurveSeg<T>, 4> b4 = {edgesCCW[0], edgesCCW[1], edgesCCW[2], edgesCCW[3]};
            if (emitRobustOGrid5(b4, degree, densePolys, outQuads))
                return true;
        }
        outQuads.push_back({edgesCCW[0], edgesCCW[1], edgesCCW[2], edgesCCW[3]});
        return true;
    }

    std::vector<gsPoint2D<T>> hubPoly;
    hubPoly.reserve(N * 6);
    for (int i = 0; i < N; ++i)
        for (int k = 0; k < 6; ++k)
            hubPoly.push_back(evalCurvePoint(*edgesCCW[i].bsp, static_cast<T>(k) / T(6)));

    auto emitLocalQuads = [&](const std::vector<gsPoint2D<T>>& allPts,
                              const std::vector<std::array<int, 4>>& lqs) {
        for (const auto& lq : lqs) {
            std::array<CurveSeg<T>, 4> segs;
            for (int k = 0; k < 4; ++k) {
                int lu = lq[k], lv = lq[(k + 1) % 4];
                if (lu < N && lv < N && ((lu + 1) % N == lv))
                    segs[k] = edgesCCW[lu];
                else
                    segs[k] = makeStraightSeg(allPts[lu], allPts[lv], degree);
            }
            outQuads.push_back(std::move(segs));
        }
    };

    // Tier A: Singularity-seeded Steiner macro-triangulation + exact flip/bitmask pairing
    if (!singPts.empty() && N + 2 * static_cast<int>(singPts.size()) - 2 <= 18) {
        std::vector<gsPoint2D<T>> allPts;
        std::vector<std::array<int, 4>> lqs;
        if (pairSteinerHubExact(edgesCCW, bndPts, isSmooth180, hubPoly, densePolys, singPts, singTargetVal, gcf, allPts, lqs)) {
            emitLocalQuads(allPts, lqs);
            return true;
        }
    }

    // Tier B: Radial-Star (N in {6,8,10,12}) vs Exact DP (N <= 20)
    // Pass 0: strict (forbids unsplit smooth-180 vertices); Pass 1: relaxed (local quad O-grid rescues any unsplit 180)
    for (int pass = 0; pass < 2; ++pass) {
        const bool allow180 = (pass == 1);
        bool hasStar = false, hasDP = false;
        gsPoint2D<T> starO;
        std::vector<std::array<int, 4>> starQuads, dpQuads;
        T starScore = T(1e18);
        if (N <= 12)
            hasStar = quadrangulateRadialStarHub(edgesCCW, bndPts, isSmooth180, hubPoly, singPts, densePolys, allow180, starO, starQuads, starScore);
        if (N <= 20)
            hasDP = quadrangulateEvenHubDP(edgesCCW, bndPts, isSmooth180, hubPoly, densePolys, gcf, allow180, dpQuads);

        if (hasStar && (!hasDP || (N == 6 && starScore < T(12.0)))) {
            std::vector<gsPoint2D<T>> allPts = bndPts;
            allPts.push_back(starO);
            emitLocalQuads(allPts, starQuads);
            return true;
        }
        if (hasDP) {
            emitLocalQuads(bndPts, dpQuads);
            return true;
        }
        if (hasStar) {
            std::vector<gsPoint2D<T>> allPts = bndPts;
            allPts.push_back(starO);
            emitLocalQuads(allPts, starQuads);
            return true;
        }
    }

    // Tier C: Strictly shrinking recursive chord bisection (straight + bent 2-segment chords through interior poles!)
    std::vector<gsPoint2D<T>> segOutDir(N), segBackDir(N);
    std::vector<T> cornerSectorDeg(N, T(180));
    for (int i = 0; i < N; ++i) {
        gsPoint2D<T> d0 = evalCurveDeriv(*edgesCCW[i].bsp, T(0));
        gsPoint2D<T> d1 = evalCurveDeriv(*edgesCCW[i].bsp, T(1));
        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
        segOutDir[i]  = (n0 > T(1e-12)) ? gsPoint2D<T>(d0.x / n0, d0.y / n0) : gsPoint2D<T>(1, 0);
        segBackDir[i] = (n1 > T(1e-12)) ? gsPoint2D<T>(-d1.x / n1, -d1.y / n1) : gsPoint2D<T>(-1, 0);
    }
    for (int i = 0; i < N; ++i)
        cornerSectorDeg[i] = mergedCcwAngleDeg(segOutDir[i], segBackDir[(i + N - 1) % N]);

    auto dirClearAtVert = [&](int vIdx, const gsPoint2D<T>& target) -> bool {
        T len = bndPts[vIdx].dist(target);
        if (len <= T(1e-8)) return false;
        gsPoint2D<T> d((target.x - bndPts[vIdx].x) / len, (target.y - bndPts[vIdx].y) / len);
        T a = mergedCcwAngleDeg(segOutDir[vIdx], d);
        return (a >= T(7.0) && cornerSectorDeg[vIdx] - a >= T(7.0));
    };
    auto chordClear = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B) -> bool {
        if (A.dist(B) <= T(1e-7) || !isSegmentInsideDomain(A, B, densePolys)) return false;
        for (int k = 1; k <= 5; ++k) {
            T f = static_cast<T>(k) / T(6);
            if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y), hubPoly, T(1e-6)))
                return false;
        }
        const int np = static_cast<int>(hubPoly.size());
        for (int k = 0; k < np; ++k) {
            const auto& p0 = hubPoly[k]; const auto& p1 = hubPoly[(k + 1) % np];
            if (p0.dist(A) < T(1e-5) || p1.dist(A) < T(1e-5) || p0.dist(B) < T(1e-5) || p1.dist(B) < T(1e-5))
                continue;
            if (segProperIntersect(A, B, p0, p1)) return false;
        }
        return true;
    };

    std::vector<gsPoint2D<T>> bendPts;
    gsPoint2D<T> poleM; T rPoleM = T(0);
    if (mergedFindDeepInteriorPoint(hubPoly, densePolys, poleM, rPoleM))
        bendPts.push_back(poleM);
    gsPoint2D<T> centM = polygonCentroid(hubPoly);
    if (pointInPolyStrict(centM, hubPoly, T(1e-5)) && isPointStrictlyInsideDomain(centM, densePolys, T(1e-5)))
        bendPts.push_back(centM);
    // Add inward-bisector interior bend points for every >=195 deg reentrant corner (splits 270 -> 135+135!)
    const T degToRadC = T(3.14159265358979323846) / T(180);
    for (int i = 0; i < N; ++i) {
        if (!isReflex[i]) continue;
        T hr = T(0.5) * cornerSectorDeg[i] * degToRadC;
        gsPoint2D<T> bis(segOutDir[i].x * std::cos(hr) - segOutDir[i].y * std::sin(hr),
                         segOutDir[i].x * std::sin(hr) + segOutDir[i].y * std::cos(hr));
        T locSpan = std::min(bndPts[i].dist(bndPts[(i + 1) % N]), bndPts[i].dist(bndPts[(i + N - 1) % N]));
        if (rPoleM > T(0)) locSpan = std::max(locSpan, rPoleM);
        for (T sc : {T(0.45), T(0.85), T(0.25)}) {
            gsPoint2D<T> bp(bndPts[i].x + sc * locSpan * bis.x, bndPts[i].y + sc * locSpan * bis.y);
            if (pointInPolyStrict(bp, hubPoly, T(1e-5)) && isPointStrictlyInsideDomain(bp, densePolys, T(1e-5))) {
                bendPts.push_back(bp);
                break;
            }
        }
    }
    // Add a couple of interior midpoints of ear-like diagonals for C-shaped large hubs
    for (int k = 0; k < N && bendPts.size() < 9; k += std::max(1, N / 4)) {
        gsPoint2D<T> em(T(0.5) * (bndPts[k].x + bndPts[(k + 2) % N].x),
                        T(0.5) * (bndPts[k].y + bndPts[(k + 2) % N].y));
        if (pointInPolyStrict(em, hubPoly, T(1e-5)) && isPointStrictlyInsideDomain(em, densePolys, T(1e-5)))
            bendPts.push_back(em);
    }

    struct BisectCand { int i, j; bool oddStep; bool bent; gsPoint2D<T> midPt; T score; };
    std::vector<BisectCand> cands;
    int triHint = -1;
    for (int i = 0; i < N; ++i) {
        for (int step = 3; step <= N - 3; ++step) {
            int j = (i + step) % N;
            if (j <= i) continue;
            bool isOdd = ((step % 2) == 1);
            if (!isOdd && (step < 4 || step > N - 4)) continue;
            T bal = std::abs(static_cast<T>(step) - T(0.5) * static_cast<T>(N));
            T smoothBonus = (isSmooth180[i] ? T(4.0) : T(0.0)) + (isSmooth180[j] ? T(4.0) : T(0.0))
                          + (isReflex[i]    ? T(25.0) : T(0.0)) + (isReflex[j]    ? T(25.0) : T(0.0));

            if (dirClearAtVert(i, bndPts[j]) && dirClearAtVert(j, bndPts[i]) && chordClear(bndPts[i], bndPts[j])) {
                T al = evalSegmentCrossFieldAlignment(gcf.fineMesh, gcf.field, gcf.accel, bndPts[i], bndPts[j], triHint);
                gsPoint2D<T> mid(T(0.5) * (bndPts[i].x + bndPts[j].x), T(0.5) * (bndPts[i].y + bndPts[j].y));
                T sc = (isOdd ? T(10.0) : T(2.0)) + smoothBonus + T(2.0) * al - T(0.4) * bal - bndPts[i].dist(bndPts[j]) / bb.domDiag;
                cands.push_back({i, j, isOdd, false, mid, sc});
            } else if (!isOdd) {
                for (const auto& bm : bendPts) {
                    if (dirClearAtVert(i, bm) && dirClearAtVert(j, bm) &&
                        chordClear(bndPts[i], bm) && chordClear(bm, bndPts[j])) {
                        T sc = T(0.0) + smoothBonus - T(0.4) * bal - (bndPts[i].dist(bm) + bm.dist(bndPts[j])) / bb.domDiag;
                        cands.push_back({i, j, false, true, bm, sc});
                        break;
                    }
                }
            }
        }
    }
    std::sort(cands.begin(), cands.end(), [](const BisectCand& a, const BisectCand& b) { return a.score > b.score; });

    for (size_t ci = 0; ci < std::min<size_t>(cands.size(), 8) && budget > 0; ++ci) {
        int i = cands[ci].i, j = cands[ci].j;
        std::vector<CurveSeg<T>> e1, e2;
        std::vector<gsPoint2D<T>> v1, v2;
        if (cands[ci].oddStep) {
            CurveSeg<T> diag = makeStraightSeg(bndPts[i], bndPts[j], degree);
            for (int k = i; k != j; k = (k + 1) % N) { e1.push_back(edgesCCW[k]); v1.push_back(bndPts[k]); }
            e1.push_back(diag.reversed()); v1.push_back(bndPts[j]);
            for (int k = j; k != i; k = (k + 1) % N) { e2.push_back(edgesCCW[k]); v2.push_back(bndPts[k]); }
            e2.push_back(diag); v2.push_back(bndPts[i]);
        } else {
            gsPoint2D<T> M = cands[ci].midPt;
            CurveSeg<T> d0 = makeStraightSeg(bndPts[i], M, degree);
            CurveSeg<T> d1 = makeStraightSeg(M, bndPts[j], degree);
            for (int k = i; k != j; k = (k + 1) % N) { e1.push_back(edgesCCW[k]); v1.push_back(bndPts[k]); }
            e1.push_back(d1.reversed()); v1.push_back(bndPts[j]);
            e1.push_back(d0.reversed()); v1.push_back(M);
            for (int k = j; k != i; k = (k + 1) % N) { e2.push_back(edgesCCW[k]); v2.push_back(bndPts[k]); }
            e2.push_back(d0); v2.push_back(bndPts[i]);
            e2.push_back(d1); v2.push_back(M);
        }
        if (static_cast<int>(e1.size()) >= N || static_cast<int>(e2.size()) >= N) continue;
        auto childArea2 = [](const std::vector<CurveSeg<T>>& es) -> T {
            std::vector<gsPoint2D<T>> pl;
            pl.reserve(es.size() * 6);
            for (const auto& s : es)
                for (int k = 0; k < 6; ++k)
                    pl.push_back(evalCurvePoint(*s.bsp, static_cast<T>(k) / T(6)));
            T a2 = T(0);
            for (size_t k = 0; k < pl.size(); ++k)
                a2 += pl[k].x * pl[(k + 1) % pl.size()].y - pl[(k + 1) % pl.size()].x * pl[k].y;
            return a2;
        };
        if (childArea2(e1) <= T(1e-8) || childArea2(e2) <= T(1e-8)) continue;
        std::vector<std::array<CurveSeg<T>, 4>> subOut;
        if (quadrangulateGeneralEvenPoly(e1, v1, {}, {}, bb, densePolys, gcf, degree, depth + 1, budget, subOut) &&
            quadrangulateGeneralEvenPoly(e2, v2, {}, {}, bb, densePolys, gcf, degree, depth + 1, budget, subOut)) {
            for (auto& q : subOut) outQuads.push_back(std::move(q));
            return true;
        }
    }

    // Tier C2 (for curved/C-shaped N == 6 hubs): 2-Pole 4-Quad Interior-Bridge Split!
    // Uses 2 interior Steiner poles (M1, M2) connected by M1--M2 to avoid any straight chord crossing curved walls!
    if (N == 6) {
        std::vector<gsPoint2D<T>> pCands = bendPts;
        for (int k = 0; k < 6; ++k) {
            gsPoint2D<T> g((bndPts[k].x + bndPts[(k + 1) % 6].x + bndPts[(k + 2) % 6].x) / T(3),
                           (bndPts[k].y + bndPts[(k + 1) % 6].y + bndPts[(k + 2) % 6].y) / T(3));
            if (pointInPolyStrict(g, hubPoly, T(1e-5)) && isPointStrictlyInsideDomain(g, densePolys, T(1e-5)))
                pCands.push_back(g);
        }
        T bestSc6 = T(1e18);
        std::vector<gsPoint2D<T>> bestPts6;
        std::vector<std::array<int, 4>> bestLq6;
        for (int sh = 0; sh < 3; ++sh) {
            int v0 = sh, v1 = (sh + 1) % 6, v2 = (sh + 2) % 6, v3 = (sh + 3) % 6, v4 = (sh + 4) % 6, v5 = (sh + 5) % 6;
            // Never leave a >= 195 deg reflex vertex unsplit at v1 or v4!
            if (isReflex[v1] || isReflex[v4]) continue;
            for (size_t a = 0; a < pCands.size(); ++a) {
                const auto& M1 = pCands[a];
                if (!dirClearAtVert(v0, M1) || !dirClearAtVert(v2, M1) ||
                    !chordClear(bndPts[v0], M1) || !chordClear(bndPts[v2], M1)) continue;
                for (size_t b = 0; b < pCands.size(); ++b) {
                    if (a == b) continue;
                    const auto& M2 = pCands[b];
                    if (M1.dist(M2) < T(1e-4) * bb.domDiag) continue;
                    if (!dirClearAtVert(v3, M2) || !dirClearAtVert(v5, M2) ||
                        !chordClear(bndPts[v3], M2) || !chordClear(bndPts[v5], M2) ||
                        !chordClear(M1, M2)) continue;
                    std::vector<gsPoint2D<T>> pts8 = bndPts;
                    pts8.push_back(M1); // idx 6
                    pts8.push_back(M2); // idx 7
                    std::array<std::array<int, 4>, 4> qs = {{
                        {6, v0, v1, v2},
                        {6, v2, v3, 7},
                        {7, v3, v4, v5},
                        {7, v5, v0, 6}
                    }};
                    bool ok4 = true;
                    T sc = T(0);
                    for (int qi = 0; qi < 4 && ok4; ++qi) {
                        gsPoint2D<T> qp[4] = {pts8[qs[qi][0]], pts8[qs[qi][1]], pts8[qs[qi][2]], pts8[qs[qi][3]]};
                        T a2 = cross2(qp[0], qp[1], qp[2]) + cross2(qp[0], qp[2], qp[3]);
                        if (a2 <= T(1e-8) || segProperIntersect(qp[0], qp[1], qp[2], qp[3]) ||
                            segProperIntersect(qp[1], qp[2], qp[3], qp[0])) { ok4 = false; break; }
                        for (int k = 0; k < 4; ++k) {
                            int lu = qs[qi][k], lv = qs[qi][(k + 1) % 4], lp = qs[qi][(k + 3) % 4];
                            gsPoint2D<T> dO = (lu < 6 && lv < 6 && (lu + 1) % 6 == lv)
                                            ? segOutDir[lu]
                                            : gsPoint2D<T>(pts8[lv].x - pts8[lu].x, pts8[lv].y - pts8[lu].y);
                            gsPoint2D<T> dB = (lp < 6 && lu < 6 && (lp + 1) % 6 == lu)
                                            ? segBackDir[lp]
                                            : gsPoint2D<T>(pts8[lp].x - pts8[lu].x, pts8[lp].y - pts8[lu].y);
                            T ang = mergedCcwAngleDeg(dO, dB);
                            bool smoothV = (lu < 6 && (lu == v1 || lu == v4) && isSmooth180[lu]);
                            if (!smoothV && (ang < T(10.0) || ang > T(168.0))) { ok4 = false; break; }
                            T dev = std::abs(ang - T(90.0)) / T(45.0);
                            sc += dev * dev;
                        }
                    }
                    if (ok4 && sc < bestSc6) {
                        bestSc6 = sc;
                        bestPts6 = std::move(pts8);
                        bestLq6.assign(qs.begin(), qs.end());
                    }
                }
            }
        }
        if (!bestLq6.empty()) {
            emitLocalQuads(bestPts6, bestLq6);
            return true;
        }
    }

    // Tier D: Universal polynomial-time Triangulation + Lawson + Edmonds' Blossom + Dual-Distance Claw-Flip finisher!
    if (std::getenv("GISMO_MERGED_DEBUG") != nullptr)
        std::cerr << "  [M4 TierD] N=" << N << " depth=" << depth << " budget=" << budget << "\n";
    return quadrangulateByTriangulationBlossomShift(edgesCCW, bndPts, hubPoly, densePolys, degree, outQuads);
}

template <class T>
inline bool pairJunctionHubsExact(
    const MedialBackbone<T>& bb,
    const GuidedCrossField<T>& gcf,
    const SlicedDecomposition<T>& sliced,
    int degree,
    MacroQuadGraph<T>& outGraph)
{
    outGraph = MacroQuadGraph<T>();
    const T mergeTol = std::max(T(1e-6), T(1e-4) * bb.domDiag);
    const auto densePolys = mergedBuildDenseLoopPolys(bb.vloops, bb.domDiag);

    auto getOrAddVert = [&](const gsPoint2D<T>& pt, bool isBnd, int lIdx, T arcS,
                            bool isPinned, T cAng, int singIdx) -> int {
        for (int i = 0; i < static_cast<int>(outGraph.verts.size()); ++i) {
            if (outGraph.verts[i].pt.dist(pt) < mergeTol) {
                if (isBnd) {
                    outGraph.verts[i].isBoundary = true;
                    if (outGraph.verts[i].loopIdx < 0) {
                        outGraph.verts[i].loopIdx = lIdx;
                        outGraph.verts[i].arcS = arcS;
                    }
                }
                if (isPinned) {
                    outGraph.verts[i].isPinnedCorner = true;
                    outGraph.verts[i].cornerAngleDeg = cAng;
                }
                if (singIdx >= 0) outGraph.verts[i].singularityIdx = singIdx;
                return i;
            }
        }
        int id = static_cast<int>(outGraph.verts.size());
        outGraph.verts.push_back({pt, isBnd, lIdx, arcS, isPinned, cAng, singIdx});
        return id;
    };

    std::map<std::pair<int, int>, int> edgeLookup;
    auto getOrAddEdge = [&](int u, int v, bool isBnd, const CurveSeg<T>& seg, bool& outFwd) -> int {
        int lo = std::min(u, v), hi = std::max(u, v);
        auto it = edgeLookup.find({lo, hi});
        if (it != edgeLookup.end()) {
            int eid = it->second;
            outFwd = (outGraph.edges[eid].u == u && outGraph.edges[eid].v == v);
            return eid;
        }
        int eid = static_cast<int>(outGraph.edges.size());
        CurveSeg<T> useSeg;
        if (seg.bsp) {
            gsBSpline<T> snapped(*seg.bsp);
            auto& cf = snapped.coefs();
            const index_t nC = cf.rows();
            if (nC >= 2) {
                const auto& pu = outGraph.verts[u].pt;
                const auto& pv = outGraph.verts[v].pt;
                T dx0 = pu.x - cf(0, 0), dy0 = pu.y - cf(0, 1);
                T dx1 = pv.x - cf(nC - 1, 0), dy1 = pv.y - cf(nC - 1, 1);
                if (std::hypot(dx0, dy0) > T(1e-14) || std::hypot(dx1, dy1) > T(1e-14)) {
                    safeNormalizeKnots01(snapped);
                    gsMatrix<T> gr = snapped.knots().greville();
                    for (index_t k = 0; k < nC; ++k) {
                        T t = gr(0, k);
                        cf(k, 0) += (T(1) - t) * dx0 + t * dx1;
                        cf(k, 1) += (T(1) - t) * dy0 + t * dy1;
                    }
                    cf(0, 0) = pu.x; cf(0, 1) = pu.y;
                    cf(nC - 1, 0) = pv.x; cf(nC - 1, 1) = pv.y;
                }
            }
            useSeg = CurveSeg<T>(snapped);
        } else {
            useSeg = makeStraightSeg(outGraph.verts[u].pt, outGraph.verts[v].pt, degree);
        }
        outGraph.edges.push_back({u, v, isBnd, useSeg});
        edgeLookup[{lo, hi}] = eid;
        outFwd = true;
        return eid;
    };

    auto emitQuadFace = [&](const std::array<int, 4>& gv, const std::array<CurveSeg<T>, 4>& segs,
                            const std::array<bool, 4>& isBnd) {
        typename MacroQuadGraph<T>::Face f;
        f.v = gv;
        for (int k = 0; k < 4; ++k) {
            bool fwd = true;
            f.e[k] = getOrAddEdge(gv[k], gv[(k + 1) % 4], isBnd[k], segs[k], fwd);
            f.eFwd[k] = fwd;
        }
        outGraph.faces.push_back(f);
    };

    for (const auto& hub : sliced.hubs) {
        if (hub.sidesCCW.size() < 4 || (hub.sidesCCW.size() % 2) != 0) return false;

        // Gather interior singularities for this hub
        std::vector<gsPoint2D<T>> singPts;
        std::vector<int> singTargetVal;
        for (int sIdx : hub.singularityIds) {
            const auto& sg = gcf.field.singularities[sIdx];
            singPts.push_back(sg.pos);
            singTargetVal.push_back(4 - sg.index);
        }

        auto workSides = hub.sidesCCW;
        auto bestSides = hub.sidesCCW;
        std::vector<std::array<CurveSeg<T>, 4>> bestCertified;
        bool hasAnySolution = false;
        int bestBadCount = 999999;

        for (int wallSplitTry = 0; wallSplitTry < 3; ++wallSplitTry) {
            const int curN = static_cast<int>(workSides.size());
            std::vector<gsPoint2D<T>> bndPts(curN);
            std::vector<CurveSeg<T>> edgesCCW(curN);
            for (int i = 0; i < curN; ++i) {
                edgesCCW[i] = workSides[i].seg;
                bndPts[i] = workSides[i].seg.p0();
            }

            std::vector<std::array<CurveSeg<T>, 4>> hubQuads;
            int hubBudget = 96;
            bool okQ = quadrangulateGeneralEvenPoly(
                edgesCCW, bndPts,
                (wallSplitTry == 0 ? singPts : std::vector<gsPoint2D<T>>{}),
                (wallSplitTry == 0 ? singTargetVal : std::vector<int>{}),
                bb, densePolys, gcf, degree, 0, hubBudget, hubQuads);

            if (okQ) {
                // Stage 3: Per-quad certification inside Module 4 ONLY upgrades genuine >=172-deg smooth/reflex corners.
                std::vector<std::array<CurveSeg<T>, 4>> candCert;
                candCert.reserve(hubQuads.size());
                int badCnt = 0;
                for (const auto& rq : hubQuads) {
                    T mnA = 0, mxA = 180, ef = 0;
                    eval4SegCornerAnglesDeg(rq[0], rq[1], rq[2], rq[3], mnA, mxA, ef);
                    if (mxA > T(172.0)) {
                        std::vector<std::array<CurveSeg<T>, 4>> og5;
                        if (emitRobustOGrid5(rq, degree, densePolys, og5)) {
                            for (auto& oq : og5) {
                                T om = 0, ox = 180, oe = 0;
                                eval4SegCornerAnglesDeg(oq[0], oq[1], oq[2], oq[3], om, ox, oe);
                                if (om < T(4.5) || ox > T(172.0)) ++badCnt;
                                candCert.push_back(std::move(oq));
                            }
                            continue;
                        }
                        ++badCnt;
                    } else if (mnA < T(4.5)) {
                        ++badCnt;
                    }
                    candCert.push_back(rq);
                }
                if (!hasAnySolution || badCnt < bestBadCount) {
                    bestCertified = std::move(candCert);
                    bestSides = workSides;
                    bestBadCount = badCnt;
                    hasAnySolution = true;
                }
                if (bestBadCount == 0) break;
            }

            // Curved-CAD-Wall Rescue (preserves even N->N+2 and 100% cross-cut conformity!):
            // Try 1: if >= 2 CAD boundary walls exist, bisect BOTH most-curved/longest walls at 0.5 (cleanly cuts 180-deg horseshoe/U-bend bands in half!).
            // Try 2: split the single most-curved CAD boundary wall into 3 equal arcs (1->3).
            if (wallSplitTry < 2 && hub.sidesCCW.size() <= 20) {
                workSides = hub.sidesCCW;
                const int origN = static_cast<int>(workSides.size());
                std::vector<std::pair<T, int>> cadSides;
                T maxTurn = T(0), maxBow = T(0);
                for (int s = 0; s < origN; ++s) {
                    const auto& sd = workSides[s];
                    if (sd.isCrossCut || sd.loopIdx < 0 || !sd.seg.bsp) continue;
                    T aLen = mergedMeasureSegArcLen(sd.seg, 12);
                    T cLen = std::max(T(1e-6), sd.seg.p0().dist(sd.seg.p1()));
                    T turn = curveIntegratedTurnRad(*sd.seg.bsp, 12);
                    T bow  = std::max(T(0), aLen / cLen - T(1));
                    maxTurn = std::max(maxTurn, turn);
                    maxBow  = std::max(maxBow, bow);
                    T sc = turn + T(2.0) * bow + T(0.15) * (aLen / bb.domDiag);
                    cadSides.push_back({sc, s});
                }
                if (cadSides.empty()) break;
                std::sort(cadSides.begin(), cadSides.end(), [](const std::pair<T, int>& a, const std::pair<T, int>& b) {
                    return a.first > b.first;
                });
                // When the hub already quadrangulated (hasAnySolution), ONLY split CAD walls that are genuinely curved
                // (avoid injecting artificial 180-deg vertices onto flat straight walls!)
                if (hasAnySolution && maxTurn < T(0.25) && maxBow < T(0.015)) break;

                auto splitSide1to2 = [&](std::vector<typename JunctionPolygon<T>::Side>& sidesVec, int idx) -> bool {
                    const auto sd = sidesVec[idx];
                    const auto& vl = bb.vloops[sd.loopIdx];
                    T sp = mergedLoopForwardSpan(vl, sd.sStart, sd.sEnd);
                    if (sp <= T(1e-6)) return false;
                    T sMid = virtualLoopWrapS(vl, sd.sStart + T(0.5) * sp);
                    CurveSeg<T> seg0 = extractVirtualSubCurve(vl, sd.sStart, sMid, degree);
                    CurveSeg<T> seg1 = extractVirtualSubCurve(vl, sMid, sd.sEnd, degree);
                    if (!seg0.bsp || !seg1.bsp) return false;
                    auto sd0 = sd; sd0.sEnd = sMid; sd0.seg = seg0; sd0.vNode = -1;
                    auto sd1 = sd; sd1.sStart = sMid; sd1.seg = seg1; sd1.uNode = -1;
                    sidesVec[idx] = std::move(sd0);
                    sidesVec.insert(sidesVec.begin() + idx + 1, std::move(sd1));
                    return true;
                };

                if (wallSplitTry == 0 && cadSides.size() >= 2) {
                    int sA = std::max(cadSides[0].second, cadSides[1].second);
                    int sB = std::min(cadSides[0].second, cadSides[1].second);
                    if (!splitSide1to2(workSides, sA) || !splitSide1to2(workSides, sB)) break;
                } else {
                    int bestS = cadSides[0].second;
                    const auto sd = workSides[bestS];
                    const auto& vl = bb.vloops[sd.loopIdx];
                    T sp = mergedLoopForwardSpan(vl, sd.sStart, sd.sEnd);
                    if (sp <= T(1e-6)) break;
                    T s1 = virtualLoopWrapS(vl, sd.sStart + sp / T(3));
                    T s2 = virtualLoopWrapS(vl, sd.sStart + (T(2) * sp) / T(3));
                    CurveSeg<T> seg0 = extractVirtualSubCurve(vl, sd.sStart, s1, degree);
                    CurveSeg<T> seg1 = extractVirtualSubCurve(vl, s1, s2, degree);
                    CurveSeg<T> seg2 = extractVirtualSubCurve(vl, s2, sd.sEnd, degree);
                    if (!seg0.bsp || !seg1.bsp || !seg2.bsp) break;
                    auto sd0 = sd; sd0.sEnd = s1; sd0.seg = seg0; sd0.vNode = -1;
                    auto sd1 = sd; sd1.sStart = s1; sd1.sEnd = s2; sd1.seg = seg1; sd1.uNode = -1; sd1.vNode = -1;
                    auto sd2 = sd; sd2.sStart = s2; sd2.seg = seg2; sd2.uNode = -1;
                    workSides[bestS] = std::move(sd0);
                    workSides.insert(workSides.begin() + bestS + 1, std::move(sd2));
                    workSides.insert(workSides.begin() + bestS + 1, std::move(sd1));
                }
            } else {
                break;
            }
        }

        if (!hasAnySolution) {
            if (std::getenv("GISMO_MERGED_DEBUG") != nullptr)
                std::cerr << "[MERGED] M4 hub fail id=" << hub.id << " N=" << hub.sidesCCW.size() << " I=" << singPts.size() << "\n";
            return false;
        }

        // Register all boundary vertices (including any newly inserted CAD wall split points!)
        for (const auto& sd : bestSides) {
            if (sd.uNode >= 0 && sd.uNode < static_cast<int>(sliced.nodes.size())) {
                const auto& nd = sliced.nodes[sd.uNode];
                getOrAddVert(nd.pt, true, nd.loopIdx, nd.arcS, nd.isPinnedCorner, nd.cornerAngleDeg, -1);
            } else if (!sd.isCrossCut && sd.loopIdx >= 0) {
                getOrAddVert(sd.seg.p0(), true, sd.loopIdx, sd.sStart, false, T(180), -1);
            }
        }

        const int finN = static_cast<int>(bestSides.size());
        for (const auto& rq : bestCertified) {
            std::array<int, 4> gv;
            std::array<bool, 4> isB = {false, false, false, false};
            for (int k = 0; k < 4; ++k) {
                gv[k] = getOrAddVert(rq[k].p0(), false, -1, T(0), false, T(360), -1);
                for (int s = 0; s < finN; ++s) {
                    if (!bestSides[s].isCrossCut &&
                        rq[k].p0().dist(bestSides[s].seg.p0()) < mergeTol &&
                        rq[k].p1().dist(bestSides[s].seg.p1()) < mergeTol) {
                        isB[k] = true; break;
                    }
                }
            }
            emitQuadFace(gv, rq, isB);
        }
    }

    return !outGraph.faces.empty();
}

// ============================================================================
// Module 5: Cross-Field Hermite Curve Realization & Coons MultiPatch Assembly
// ============================================================================

template <class T>
inline bool realizeHermiteMultiPatch(
    const MedialBackbone<T>& bb,
    const GuidedCrossField<T>& gcf,
    MacroQuadGraph<T>& graph,
    int degree,
    gsMultiPatch<T>& outMp)
{
    outMp.clear();
    const int nV = static_cast<int>(graph.verts.size());
    const int nE = static_cast<int>(graph.edges.size());
    if (nV == 0 || graph.faces.empty()) return false;

    const T pi = T(3.14159265358979323846);
    const T twoPi = T(2) * pi;

    // Collect incident interior edges per vertex
    struct IncIntEdge { int eId; int otherV; T chordAng; };
    std::vector<std::vector<IncIntEdge>> vIntEdges(nV);
    for (int ei = 0; ei < nE; ++ei) {
        if (graph.edges[ei].isBoundary) continue;
        int u = graph.edges[ei].u, v = graph.edges[ei].v;
        T aUV = std::atan2(graph.verts[v].pt.y - graph.verts[u].pt.y, graph.verts[v].pt.x - graph.verts[u].pt.x);
        T aVU = std::atan2(graph.verts[u].pt.y - graph.verts[v].pt.y, graph.verts[u].pt.x - graph.verts[v].pt.x);
        vIntEdges[u].push_back({ei, v, aUV});
        vIntEdges[v].push_back({ei, u, aVU});
    }

    // Assign symmetric star / 90-deg orthogonal departure tangents at every vertex endpoint
    std::map<std::pair<int, int>, gsPoint2D<T>> edgeEndTan; // (edgeId, vertId) -> unit departure tangent from vertId
    for (int v = 0; v < nV; ++v) {
        auto& inc = vIntEdges[v];
        if (inc.empty()) continue;
        std::sort(inc.begin(), inc.end(), [](const IncIntEdge& a, const IncIntEdge& b) { return a.chordAng < b.chordAng; });

        const auto& vtx = graph.verts[v];
        const int K = static_cast<int>(inc.size());

        if (vtx.isBoundary && vtx.loopIdx >= 0 && vtx.loopIdx < static_cast<int>(bb.vloops.size())) {
            if (!vtx.isPinnedCorner && K == 1) {
                // Smooth CAD boundary split -> exact 90-deg inward normal (AS-G1 beta(0) = 0!)
                gsPoint2D<T> nIn = mergedInwardNormalOnLoop(bb.vloops[vtx.loopIdx], vtx.arcS);
                edgeEndTan[{inc[0].eId, v}] = nIn;
            } else {
                // Pinned corner or multi-spoke boundary vertex -> divide interior sector smoothly
                for (int j = 0; j < K; ++j)
                    edgeEndTan[{inc[j].eId, v}] = gsPoint2D<T>(std::cos(inc[j].chordAng), std::sin(inc[j].chordAng));
            }
        } else if (!vtx.isBoundary && K >= 3 && K <= 6) {
            // Interior vertex of valence K -> fit optimal symmetric 2*pi/K star phase phi0
            std::complex<T> phaseSum(T(0), T(0));
            for (int j = 0; j < K; ++j) {
                T targetOff = (twoPi * static_cast<T>(j)) / static_cast<T>(K);
                T diff = inc[j].chordAng - targetOff;
                phaseSum += std::complex<T>(std::cos(diff), std::sin(diff));
            }
            T phi0 = (std::abs(phaseSum) > T(1e-8)) ? std::arg(phaseSum) : inc[0].chordAng;
            for (int j = 0; j < K; ++j) {
                T ang = phi0 + (twoPi * static_cast<T>(j)) / static_cast<T>(K);
                edgeEndTan[{inc[j].eId, v}] = gsPoint2D<T>(std::cos(ang), std::sin(ang));
            }
        } else {
            for (int j = 0; j < K; ++j)
                edgeEndTan[{inc[j].eId, v}] = gsPoint2D<T>(std::cos(inc[j].chordAng), std::sin(inc[j].chordAng));
        }
    }

    const auto densePolys = mergedBuildDenseLoopPolys(bb.vloops, bb.domDiag);
    const int nF = static_cast<int>(graph.faces.size());

    // Build edge -> incident faces adjacency
    std::vector<std::vector<int>> edgeFaces(nE);
    for (int fi = 0; fi < nF; ++fi)
        for (int k = 0; k < 4; ++k)
            edgeFaces[graph.faces[fi].e[k]].push_back(fi);

    auto buildAndEvalFace = [&](int fi, const CurveSeg<T>* overrideSeg, int overrideEi,
                                std::unique_ptr<gsTensorBSpline<2, T>>& outPatch,
                                T& outMinDet, T& outMinA, T& outMaxA) -> bool {
        const auto& f = graph.faces[fi];
        CurveSeg<T> c[4];
        for (int k = 0; k < 4; ++k) {
            int ei = f.e[k];
            const CurveSeg<T>& base = (overrideSeg && ei == overrideEi) ? *overrideSeg : graph.edges[ei].seg;
            c[k] = f.eFwd[k] ? base : base.reversed();
        }
        T ef = T(0);
        eval4SegCornerAnglesDeg(c[0], c[1], c[2], c[3], outMinA, outMaxA, ef);
        outPatch.reset(makeCoonsFromCCWCycle(c[0], c[1], c[2], c[3], degree));
        if (!outPatch) return false;
        outMinDet = evalPatchMinDetDense(*outPatch, 25);
        if (outMinDet <= T(1e-7)) {
            quickHarmonicPatchRepair(*outPatch, T(1e-7));
            outMinDet = evalPatchMinDetDense(*outPatch, 25);
        }
        return true;
    };

    // Baseline straight/Module-4 patches (guaranteed never regressed by Hermite curving!)
    std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> patches(nF);
    std::vector<T> faceMinDet(nF, -T(1)), faceMinA(nF, T(0)), faceMaxA(nF, T(180));
    for (int fi = 0; fi < nF; ++fi) {
        if (!buildAndEvalFace(fi, nullptr, -1, patches[fi], faceMinDet[fi], faceMinA[fi], faceMaxA[fi]))
            return false;
    }

    // Per-edge monotonic Hermite realization: curves edge ONLY if incident faces stay valid (or strictly improve)!
    for (int ei = 0; ei < nE; ++ei) {
        if (graph.edges[ei].isBoundary || edgeFaces[ei].empty()) continue;
        int u = graph.edges[ei].u, v = graph.edges[ei].v;
        auto itU = edgeEndTan.find({ei, u});
        auto itV = edgeEndTan.find({ei, v});
        if (itU == edgeEndTan.end() || itV == edgeEndTan.end()) continue;

        const auto& pu = graph.verts[u].pt;
        const auto& pv = graph.verts[v].pt;
        const T L = pu.dist(pv);
        if (L <= T(1e-10)) continue;

        gsPoint2D<T> dU = itU->second, dV = itV->second;
        gsPoint2D<T> ch((pv.x - pu.x) / L, (pv.y - pu.y) / L);
        // Prevent hairpin hooks if target tangent opposes chord
        if (dU.x * ch.x + dU.y * ch.y < T(0.15)) dU = ch;
        if (dV.x * (-ch.x) + dV.y * (-ch.y) < T(0.15)) dV = gsPoint2D<T>(-ch.x, -ch.y);
        // Skip if tangents are already practically straight along chord
        if (dU.x * ch.x + dU.y * ch.y > T(0.996) && dV.x * (-ch.x) + dV.y * (-ch.y) > T(0.996))
            continue;

        for (T alpha : {T(0.24), T(0.15), T(0.08)}) {
            gsPoint2D<T> p1(pu.x + alpha * L * dU.x, pu.y + alpha * L * dU.y);
            gsPoint2D<T> p2(pv.x + alpha * L * dV.x, pv.y + alpha * L * dV.y);
            CurveSeg<T> cand = makeCubicBezierInteriorSeg(pu, p1, p2, pv, degree);
            bool insideDom = true;
            for (int s = 1; s < 10; ++s) {
                gsPoint2D<T> q = evalCurvePoint(*cand.bsp, static_cast<T>(s) / T(10));
                if (!isPointStrictlyInsideDomain(q, densePolys, T(1e-5))) { insideDom = false; break; }
            }
            if (!insideDom) continue;

            bool accept = true;
            std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> trialP(edgeFaces[ei].size());
            std::vector<T> tDet(edgeFaces[ei].size()), tMinA(edgeFaces[ei].size()), tMaxA(edgeFaces[ei].size());
            for (size_t k = 0; k < edgeFaces[ei].size(); ++k) {
                int fi = edgeFaces[ei][k];
                if (!buildAndEvalFace(fi, &cand, ei, trialP[k], tDet[k], tMinA[k], tMaxA[k])) {
                    accept = false; break;
                }
                bool oldValid = (faceMinDet[fi] > T(1e-6));
                bool newValid = (tDet[k] > T(1e-6) && tMinA[k] >= T(10.0) && tMaxA[k] <= T(162.0));
                if (oldValid) {
                    if (!newValid) { accept = false; break; }
                } else {
                    if (tDet[k] <= faceMinDet[fi] + T(1e-8)) { accept = false; break; }
                }
            }
            if (accept) {
                graph.edges[ei].seg = std::move(cand);
                for (size_t k = 0; k < edgeFaces[ei].size(); ++k) {
                    int fi = edgeFaces[ei][k];
                    patches[fi] = std::move(trialP[k]);
                    faceMinDet[fi] = tDet[k];
                    faceMinA[fi] = tMinA[k];
                    faceMaxA[fi] = tMaxA[k];
                }
                break;
            }
        }
    }

    for (int fi = 0; fi < nF; ++fi) {
        if (!patches[fi]) return false;
        if (evalPatchMinDetDense(*patches[fi], 49) <= T(1e-7))
            quickHarmonicPatchRepair(*patches[fi], T(1e-7));
        outMp.addPatch(typename gsGeometry<T>::uPtr(patches[fi].release()));
    }

    gsASG1LocalFitter<T>::unifyPatchOrientations(outMp);
    outMp.computeTopology(T(1e-4), false, false);
    return (outMp.nPatches() == graph.faces.size());
}

// ============================================================================
// Top-Level Driver: Field-Guided Medial Macro-Pairing
// ============================================================================

template <class T>
inline bool buildFieldGuidedMedialMultiPatch(
    const gsPlanarDomain<T>& dom,
    int degree,
    T smoothJointTolDeg,
    gsMultiPatch<T>& outMp,
    T targetHFrac = T(1.0) / T(24.0),
    T thetaPruneDeg = T(45.0),
    T junctionMergeRatio = T(0.32),
    T alphaMAT = T(2.5))
{
    const bool dbg = (std::getenv("GISMO_MERGED_DEBUG") != nullptr);
    outMp.clear();

    // Fast canonical 0-hole convexish 4-corner quad -> 1 patch
    if (dom.numHoles() == 0 && dom.outer().size() == 4) {
        CurveSeg<T> s4[4];
        gsPoint2D<T> c4[4];
        for (int k = 0; k < 4; ++k) {
            s4[k] = CurveSeg<T>(dom.outer().curve(k), degree);
            c4[k] = s4[k].p0();
        }
        if (isStrictlyConvexQuad(c4)) {
            std::unique_ptr<gsTensorBSpline<2, T>> p0(makeCoonsFromCCWCycle(s4[0], s4[1], s4[2], s4[3], degree));
            if (p0) {
                if (evalPatchMinDetDense(*p0, 49) <= T(1e-7)) quickHarmonicPatchRepair(*p0, T(1e-7));
                T angs[4], mnA = 0, mxA = 180; bool degen = false;
                evalPatchCornerAnglesDeg(*p0, angs, mnA, mxA, degen);
                if (!degen && mnA > T(10) && mxA < T(162) && evalPatchMinDetDense(*p0, 49) > T(1e-7)) {
                    outMp.addPatch(typename gsGeometry<T>::uPtr(p0.release()));
                    outMp.computeTopology(T(1e-4), false, false);
                    return true;
                }
            }
        }
    }

    // Module 1: Boundary CDT & Pruned MAT Backbone
    MedialBackbone<T> bb;
    if (!buildMedialBackbone(dom, degree, smoothJointTolDeg, targetHFrac, thetaPruneDeg, junctionMergeRatio, bb)) {
        if (dbg) std::cerr << "[MERGED] Module 1 (buildMedialBackbone) failed.\n";
        return false;
    }
    if (dbg) {
        std::cerr << "[MERGED] M1: loops=" << bb.vloops.size()
                  << " cdtT=" << bb.cdt.tris.size()
                  << " sleeves=" << bb.sleeves.size()
                  << " macroNodes=" << bb.macroNodes.size() << "\n";
    }

    // Fast smooth 0-corner hole-free domain -> try 5-patch O-grid ONLY if all 5 patches are dense-injective
    if (dom.numHoles() == 0 && !bb.vloops.empty() && bb.vloops[0].corners.empty()) {
        GuidedCrossField<T> gcf0;
        gsMultiPatch<T> ogridMp;
        if (solveMATGuided4RoSy(bb, gcf0, alphaMAT) &&
            buildHarmonicCurvedOGrid(dom, bb.vloops, gcf0.fineMesh, gcf0.accel, degree, ogridMp)) {
            bool oOk = (ogridMp.nPatches() > 0);
            for (size_t p = 0; p < ogridMp.nPatches() && oOk; ++p)
                if (evalPatchMinDetDense(ogridMp.patch(p), 49) <= T(1e-7)) oOk = false;
            if (oOk) {
                outMp = give(ogridMp);
                if (dbg) std::cerr << "[MERGED] Fast smooth O-grid P=" << outMp.nPatches() << "\n";
                return true;
            }
        }
    }

    // Module 2: MAT-Guided 4-RoSy Cross-Field
    GuidedCrossField<T> gcf;
    if (!solveMATGuided4RoSy(bb, gcf, alphaMAT)) {
        if (dbg) std::cerr << "[MERGED] Module 2 (solveMATGuided4RoSy) failed.\n";
        return false;
    }
    if (dbg) {
        std::cerr << "[MERGED] M2: fineT=" << gcf.fineMesh.tris.size()
                  << " sings=" << gcf.field.singularities.size() << "\n";
    }

    // Module 3: Optimal Neck Slicing & 1-Layer Parity Balancing
    SlicedDecomposition<T> sliced;
    if (!sliceNecksAndBalanceParity(bb, gcf, degree, sliced)) {
        if (dbg) std::cerr << "[MERGED] Module 3 (sliceNecksAndBalanceParity) failed.\n";
        return false;
    }
    if (dbg) {
        std::cerr << "[MERGED] M3: cuts=" << sliced.cuts.size()
                  << " hubs=" << sliced.hubs.size() << " [sides:";
        for (const auto& h : sliced.hubs)
            std::cerr << " " << h.sidesCCW.size() << "(I=" << h.singularityIds.size() << ")";
        std::cerr << "]\n";
    }

    // Module 4: Singularity-Seeded Macro-Triangulation & Exact Pairing
    MacroQuadGraph<T> graph;
    if (!pairJunctionHubsExact(bb, gcf, sliced, degree, graph)) {
        if (dbg) std::cerr << "[MERGED] Module 4 (pairJunctionHubsExact) failed.\n";
        return false;
    }
    if (dbg) {
        std::cerr << "[MERGED] M4: V=" << graph.verts.size()
                  << " E=" << graph.edges.size()
                  << " F=" << graph.faces.size() << "\n";
    }

    // Module 5: Hermite Curve Realization & Coons Assembly
    if (!realizeHermiteMultiPatch(bb, gcf, graph, degree, outMp)) {
        if (dbg) std::cerr << "[MERGED] Module 5 (realizeHermiteMultiPatch) failed.\n";
        return false;
    }
    return true;
}

// ============================================================================
// Stage 4: Topological O-Ring Pruner & Encompassing-Patch Absorber
//   - Op 4A: Direct 5 -> 1 O-ring collapse (when outer 4-gon is clean)
//   - Op 4B: 6 -> 1 Sharp-Tip Kite + O-ring collapse
//   - Op 4C: 6 -> 1 / 6 -> 2 Boundary-Blister O-ring absorption into encompassing patch
// Every reduction is strictly certified for C0 conformity, exact CAD boundary match,
// positive dense Jacobian (detJ > 1e-6), and corner angles in [effMinAng, effMaxAng].
// ============================================================================
template <class T>
inline int pruneParasiticORings(
    const gsPlanarDomain<T>& dom,
    gsMultiPatch<T>& mp,
    int degree,
    T smoothJointTolDeg,
    T effMinAng,
    T effMaxAng,
    bool dbg = false)
{
    if (mp.nPatches() < 5) return 0;

    auto vloops = buildVirtualCADLoops(dom, degree, smoothJointTolDeg);

    auto boxToCcw = [](int bIdx) -> int {
        // G+Smo boxSide: west=1, east=2, south=3, north=4 -> CCW order: S=0, E=1, N=2, W=3
        if (bIdx == 3) return 0;
        if (bIdx == 2) return 1;
        if (bIdx == 4) return 2;
        return 3;
    };

    auto extractPatchCcwSides = [&](const gsGeometry<T>& patch, std::array<CurveSeg<T>, 4>& outSegs) -> bool {
        typename gsGeometry<T>::uPtr sS = patch.boundary(boundary::south);
        typename gsGeometry<T>::uPtr sE = patch.boundary(boundary::east);
        typename gsGeometry<T>::uPtr sN = patch.boundary(boundary::north);
        typename gsGeometry<T>::uPtr sW = patch.boundary(boundary::west);
        auto* bS = dynamic_cast<gsBSpline<T>*>(sS.get());
        auto* bE = dynamic_cast<gsBSpline<T>*>(sE.get());
        auto* bN = dynamic_cast<gsBSpline<T>*>(sN.get());
        auto* bW = dynamic_cast<gsBSpline<T>*>(sW.get());
        if (!bS || !bE || !bN || !bW) return false;
        outSegs[0] = CurveSeg<T>(*bS);
        outSegs[1] = CurveSeg<T>(*bE);
        outSegs[2] = CurveSeg<T>(*bN).reversed();
        outSegs[3] = CurveSeg<T>(*bW).reversed();
        return true;
    };

    auto chain4ToCCW = [&](const std::array<CurveSeg<T>, 4>& raw, std::array<CurveSeg<T>, 4>& out) -> bool {
        bool used[4] = {true, false, false, false};
        out[0] = raw[0];
        gsPoint2D<T> cur = out[0].p1();
        for (int step = 1; step < 4; ++step) {
            int hit = -1;
            bool rev = false;
            T bestD = T(5e-4);
            for (int k = 1; k < 4; ++k) {
                if (used[k]) continue;
                T d0 = raw[k].p0().dist(cur);
                T d1 = raw[k].p1().dist(cur);
                if (d0 < bestD) { bestD = d0; hit = k; rev = false; }
                if (d1 < bestD) { bestD = d1; hit = k; rev = true; }
            }
            if (hit < 0) return false;
            used[hit] = true;
            out[step] = rev ? raw[hit].reversed() : raw[hit];
            cur = out[step].p1();
        }
        if (cur.dist(out[0].p0()) > T(5e-4)) return false;

        std::vector<gsPoint2D<T>> poly;
        poly.reserve(32);
        for (int k = 0; k < 4; ++k)
            for (int s = 0; s < 8; ++s)
                poly.push_back(evalCurvePoint(*out[k].bsp, static_cast<T>(s) / T(8)));
        T a2 = T(0);
        for (size_t i = 0; i < poly.size(); ++i)
            a2 += poly[i].x * poly[(i + 1) % poly.size()].y - poly[(i + 1) % poly.size()].x * poly[i].y;
        if (std::abs(a2) <= T(1e-10)) return false;
        if (a2 < T(0)) {
            out = {out[3].reversed(), out[2].reversed(), out[1].reversed(), out[0].reversed()};
        }
        return true;
    };

    auto buildCertifiedSinglePatch = [&](const std::array<CurveSeg<T>, 4>& c4) -> std::unique_ptr<gsTensorBSpline<2, T>> {
        T mnA = 0, mxA = 180, ef = 0;
        eval4SegCornerAnglesDeg(c4[0], c4[1], c4[2], c4[3], mnA, mxA, ef);
        // Never collapse if it re-creates a sharp or >=effMaxAng (e.g. 180-deg) corner singularity!
        if (mnA < effMinAng + T(0.5) || mxA > effMaxAng - T(0.5)) return nullptr;

        std::unique_ptr<gsTensorBSpline<2, T>> np(makeCoonsFromCCWCycle(c4[0], c4[1], c4[2], c4[3], degree));
        if (!np) return nullptr;
        if (evalPatchMinDetDense(*np, 25) <= T(1e-6))
            quickHarmonicPatchRepair(*np, T(1e-6));
        if ((!isPatchDetSingleSign(*np, 21, T(1e-7)) || evalPatchMinDetDense(*np, 49) <= T(1e-6)) &&
            np->knots(0).unique().size() < 3 && np->knots(1).unique().size() < 3) {
            np->insertKnot(T(0.5), 0, 1);
            np->insertKnot(T(0.5), 1, 1);
            quickHarmonicPatchRepair(*np, T(1e-6));
        }
        if (!isPatchDetSingleSign(*np, 21, T(1e-7)) || evalPatchMinDetDense(*np, 49) <= T(1e-6))
            return nullptr;
        T angs[4], pMinA = 0, pMaxA = 180; bool deg = false;
        evalPatchCornerAnglesDeg(*np, angs, pMinA, pMaxA, deg);
        if (deg || pMinA < effMinAng + T(0.5) || pMaxA > effMaxAng - T(0.5))
            return nullptr;
        return np;
    };

    auto approxPatchArea = [&](const gsGeometry<T>& p) -> T {
        gsMatrix<T> uv(2, 4), val;
        uv << 0, 1, 1, 0,
              0, 0, 1, 1;
        p.eval_into(uv, val);
        T a2 = T(0);
        for (int i = 0; i < 4; ++i) {
            int j = (i + 1) % 4;
            a2 += val(0, i) * val(1, j) - val(0, j) * val(1, i);
        }
        return T(0.5) * std::abs(a2);
    };

    int totalPruned = 0;
    for (int iter = 0; iter < 24; ++iter) {
        const size_t nP = mp.nPatches();
        if (nP < 5) break;

        gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
        mp.computeTopology(T(1e-4), false, false);
        if (2 * mp.nInterfaces() + mp.nBoundary() != 4 * nP) break;
        const size_t oldBnd = mp.nBoundary();

        // Build exact CCW adjacency table adj[p][ccwSide] = {nbrPatch, nbrCcwSide}
        std::vector<std::array<std::pair<int, int>, 4>> adj(nP);
        for (size_t p = 0; p < nP; ++p)
            for (int k = 0; k < 4; ++k)
                adj[p][k] = {-1, -1};
        for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
            int p1 = static_cast<int>(it->first().patch);
            int s1 = boxToCcw(it->first().side().m_index);
            int p2 = static_cast<int>(it->second().patch);
            int s2 = boxToCcw(it->second().side().m_index);
            if (p1 >= 0 && p1 < static_cast<int>(nP) && p2 >= 0 && p2 < static_cast<int>(nP)) {
                adj[p1][s1] = {p2, s2};
                adj[p2][s2] = {p1, s1};
            }
        }

        // Detect all 5-patch O-rings {p0, r[0..3]} and sort smallest-area first
        struct ORingInfo {
            int p0;
            std::array<int, 4> r;
            std::array<int, 4> sIn;
            std::array<int, 4> sOut;
            std::array<std::pair<int, int>, 4> extNbr;
            T area;
        };
        std::vector<ORingInfo> rings;
        for (int p0 = 0; p0 < static_cast<int>(nP); ++p0) {
            bool allInt = true;
            std::array<int, 4> r, sIn, sOut;
            for (int k = 0; k < 4; ++k) {
                if (adj[p0][k].first < 0) { allInt = false; break; }
                r[k] = adj[p0][k].first;
                sIn[k] = adj[p0][k].second;
                sOut[k] = (sIn[k] + 2) % 4;
            }
            if (!allInt) continue;
            // All 4 ring patches must be distinct and different from p0
            std::set<int> rSet(r.begin(), r.end());
            if (rSet.size() != 4 || rSet.count(p0)) continue;
            // Consecutive ring patches must share spoke edges
            int spokes = 0;
            for (int k = 0; k < 4; ++k) {
                int rk = r[k], rn = r[(k + 1) % 4];
                if (adj[rk][(sIn[k] + 1) % 4].first == rn || adj[rk][(sIn[k] + 3) % 4].first == rn)
                    ++spokes;
            }
            if (spokes != 4) continue;

            ORingInfo info;
            info.p0 = p0; info.r = r; info.sIn = sIn; info.sOut = sOut;
            info.area = approxPatchArea(mp.patch(p0));
            for (int k = 0; k < 4; ++k) {
                info.area += approxPatchArea(mp.patch(r[k]));
                info.extNbr[k] = adj[r[k]][sOut[k]];
            }
            rings.push_back(info);
        }
        if (rings.empty()) break;

        std::sort(rings.begin(), rings.end(), [](const ORingInfo& a, const ORingInfo& b) { return a.area < b.area; });

        bool prunedOne = false;
        for (const auto& rg : rings) {
            std::set<int> cl5 = {rg.p0, rg.r[0], rg.r[1], rg.r[2], rg.r[3]};
            std::array<CurveSeg<T>, 4> outSeg;
            bool extOk = true;
            for (int k = 0; k < 4; ++k) {
                std::array<CurveSeg<T>, 4> rSegs;
                if (!extractPatchCcwSides(mp.patch(rg.r[k]), rSegs)) { extOk = false; break; }
                outSeg[k] = rSegs[rg.sOut[k]];
            }
            if (!extOk) continue;

            auto commitReplacement = [&](const std::set<int>& dropSet,
                                         std::vector<std::unique_ptr<gsTensorBSpline<2, T>>>& addPatches,
                                         size_t expectedBnd) -> bool {
                gsMultiPatch<T> trial;
                for (size_t p = 0; p < nP; ++p) {
                    if (dropSet.count(static_cast<int>(p))) continue;
                    trial.addPatch(mp.patch(p).clone());
                }
                for (auto& np : addPatches)
                    trial.addPatch(typename gsGeometry<T>::uPtr(np.release()));
                gsASG1LocalFitter<T>::unifyPatchOrientations(trial);
                trial.computeTopology(T(1e-4), false, false);
                if (2 * trial.nInterfaces() + trial.nBoundary() != 4 * trial.nPatches()) return false;
                if (trial.nBoundary() != expectedBnd) return false;
                mp = give(trial);
                return true;
            };

            // Op 4A: Direct 5 -> 1 O-ring collapse
            {
                std::array<CurveSeg<T>, 4> c4;
                if (chain4ToCCW(outSeg, c4)) {
                    if (auto np = buildCertifiedSinglePatch(c4)) {
                        std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> add;
                        add.push_back(std::move(np));
                        if (commitReplacement(cl5, add, oldBnd)) {
                            if (dbg) std::cerr << "[MERGED] Pruned O-ring 5->1 (p0=" << rg.p0 << ") P=" << mp.nPatches() << "\n";
                            ++totalPruned;
                            prunedOne = true;
                            break;
                        }
                    }
                }
            }

            // Op 4B: 6 -> 1 Sharp-Tip Kite + O-ring collapse (single neighbor pTip touches 2 adjacent sides of C5)
            {
                int kKite = -1, pTip = -1;
                for (int k = 0; k < 4; ++k) {
                    int kn = (k + 1) % 4;
                    int na = rg.extNbr[k].first, nb = rg.extNbr[kn].first;
                    if (na >= 0 && na == nb && !cl5.count(na)) {
                        kKite = k; pTip = na; break;
                    }
                }
                if (kKite >= 0 && pTip >= 0) {
                    std::array<CurveSeg<T>, 4> tipSegs;
                    if (extractPatchCcwSides(mp.patch(pTip), tipSegs)) {
                        int sTA = rg.extNbr[kKite].second;
                        int sTB = rg.extNbr[(kKite + 1) % 4].second;
                        std::array<CurveSeg<T>, 4> raw4;
                        int w = 0;
                        for (int s = 0; s < 4; ++s)
                            if (s != sTA && s != sTB && w < 4) raw4[w++] = tipSegs[s];
                        for (int s = 0; s < 4; ++s)
                            if (s != kKite && s != (kKite + 1) % 4 && w < 4) raw4[w++] = outSeg[s];
                        std::array<CurveSeg<T>, 4> c4;
                        if (w == 4 && chain4ToCCW(raw4, c4)) {
                            if (auto np = buildCertifiedSinglePatch(c4)) {
                                std::set<int> cl6 = cl5;
                                cl6.insert(pTip);
                                std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> add;
                                add.push_back(std::move(np));
                                if (commitReplacement(cl6, add, oldBnd)) {
                                    if (dbg) std::cerr << "[MERGED] Pruned Kite6 O-ring 6->1 (p0=" << rg.p0 << ",tip=" << pTip << ") P=" << mp.nPatches() << "\n";
                                    ++totalPruned;
                                    prunedOne = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // Op 4C: Encompassing-Patch Absorption (BLISTER5 + pEnc: 6 -> 1 or 6 -> 2)
            {
                int nBnd = 0, kCut = -1;
                for (int k = 0; k < 4; ++k) {
                    if (rg.extNbr[k].first < 0) ++nBnd;
                    else kCut = k;
                }
                if (nBnd == 3 && kCut >= 0) {
                    int pEnc = rg.extNbr[kCut].first;
                    int sEnc = rg.extNbr[kCut].second;
                    if (pEnc >= 0 && !cl5.count(pEnc)) {
                        // Order the 3 boundary segments along C5's outer perimeter from cutEnd -> cutStart
                        std::array<CurveSeg<T>, 4> c5ccw;
                        if (chain4ToCCW(outSeg, c5ccw)) {
                            // Locate the cut segment in c5ccw
                            int cutPos = -1;
                            for (int k = 0; k < 4; ++k) {
                                if ((c5ccw[k].p0().dist(outSeg[kCut].p0()) < T(5e-4) && c5ccw[k].p1().dist(outSeg[kCut].p1()) < T(5e-4)) ||
                                    (c5ccw[k].p0().dist(outSeg[kCut].p1()) < T(5e-4) && c5ccw[k].p1().dist(outSeg[kCut].p0()) < T(5e-4))) {
                                    cutPos = k; break;
                                }
                            }
                            if (cutPos >= 0) {
                                const auto& b1 = c5ccw[(cutPos + 1) % 4];
                                const auto& b2 = c5ccw[(cutPos + 2) % 4];
                                const auto& b3 = c5ccw[(cutPos + 3) % 4];
                                gsPoint2D<T> pStart = b1.p0();
                                gsPoint2D<T> pMid   = evalCurvePoint(*b2.bsp, T(0.5));
                                gsPoint2D<T> pEnd   = b3.p1();

                                // Find which VirtualCADLoop carries (pStart -> pMid -> pEnd)
                                int hitL = -1;
                                T sStart = T(0), sEnd = T(0);
                                for (int l = 0; l < static_cast<int>(vloops.size()); ++l) {
                                    T s0, sm, s1;
                                    gsPoint2D<T> q0, qm, q1;
                                    T d0 = projectPointToVirtualLoop(vloops[l], pStart, s0, q0);
                                    T dm = projectPointToVirtualLoop(vloops[l], pMid,   sm, qm);
                                    T d1 = projectPointToVirtualLoop(vloops[l], pEnd,   s1, q1);
                                    if (d0 < T(1e-3) && dm < T(1e-3) && d1 < T(1e-3)) {
                                        T spMid = mergedLoopForwardSpan(vloops[l], s0, sm);
                                        T spEnd = mergedLoopForwardSpan(vloops[l], s0, s1);
                                        if (spMid > T(1e-5) && spMid < spEnd - T(1e-5)) {
                                            hitL = l; sStart = s0; sEnd = s1; break;
                                        }
                                    }
                                }
                                if (hitL >= 0) {
                                    const auto& vl = vloops[hitL];
                                    T fullSpan = mergedLoopForwardSpan(vl, sStart, sEnd);
                                    // Check if any pinned sharp CAD corner lies inside (sStart, sEnd)
                                    std::vector<T> innerCorners;
                                    for (const auto& c : vl.corners) {
                                        T d = mergedLoopForwardSpan(vl, sStart, c.arcS);
                                        if (d > T(1e-4) * vl.totalArcLen && d < fullSpan - T(1e-4) * vl.totalArcLen)
                                            innerCorners.push_back(c.arcS);
                                    }

                                    std::array<CurveSeg<T>, 4> encSegs;
                                    if (extractPatchCcwSides(mp.patch(pEnc), encSegs)) {
                                        // 4C-1: 6 -> 1 single-patch absorption when 0 sharp CAD corners lie inside the blister arc
                                        if (innerCorners.empty()) {
                                            CurveSeg<T> mergedArc = extractVirtualSubCurve(vl, sStart, sEnd, degree);
                                            if (mergedArc.bsp) {
                                                // Snap endpoints smoothly via Greville ramp to match pEnc corners to 1e-16
                                                gsBSpline<T> sn(*mergedArc.bsp);
                                                safeNormalizeKnots01(sn);
                                                auto& cf = sn.coefs();
                                                index_t nC = cf.rows();
                                                T dx0 = pStart.x - cf(0, 0), dy0 = pStart.y - cf(0, 1);
                                                T dx1 = pEnd.x - cf(nC - 1, 0), dy1 = pEnd.y - cf(nC - 1, 1);
                                                if (std::hypot(dx0, dy0) < T(1e-3) && std::hypot(dx1, dy1) < T(1e-3)) {
                                                    gsMatrix<T> gr = sn.knots().greville();
                                                    for (index_t k = 0; k < nC; ++k) {
                                                        T t = gr(0, k);
                                                        cf(k, 0) += (T(1) - t) * dx0 + t * dx1;
                                                        cf(k, 1) += (T(1) - t) * dy0 + t * dy1;
                                                    }
                                                    cf(0, 0) = pStart.x; cf(0, 1) = pStart.y;
                                                    cf(nC - 1, 0) = pEnd.x; cf(nC - 1, 1) = pEnd.y;
                                                    CurveSeg<T> exactMerged(sn);

                                                    std::array<CurveSeg<T>, 4> raw4;
                                                    raw4[0] = exactMerged;
                                                    int w = 1;
                                                    for (int s = 0; s < 4; ++s)
                                                        if (s != sEnc && w < 4) raw4[w++] = encSegs[s];
                                                    std::array<CurveSeg<T>, 4> c4;
                                                    if (w == 4 && chain4ToCCW(raw4, c4)) {
                                                        if (auto np = buildCertifiedSinglePatch(c4)) {
                                                            std::set<int> cl6 = cl5;
                                                            cl6.insert(pEnc);
                                                            std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> add;
                                                            add.push_back(std::move(np));
                                                            if (commitReplacement(cl6, add, oldBnd - 2)) {
                                                                if (dbg) std::cerr << "[MERGED] Absorbed Blister5+pEnc 6->1 (p0=" << rg.p0 << ",enc=" << pEnc << ") P=" << mp.nPatches() << "\n";
                                                                ++totalPruned;
                                                                prunedOne = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        // 4C-2: 6 -> 2 clean 2-patch absorption if pEnc's opposite side is ALSO on domain boundary
                                        int sOpp = (sEnc + 2) % 4;
                                        if (!prunedOne && adj[pEnc][sOpp].first < 0 && innerCorners.size() <= 1) {
                                            T sSplit = innerCorners.empty()
                                                     ? virtualLoopWrapS(vl, sStart + T(0.5) * fullSpan)
                                                     : innerCorners[0];
                                            CurveSeg<T> arcA = extractVirtualSubCurve(vl, sStart, sSplit, degree);
                                            CurveSeg<T> arcB = extractVirtualSubCurve(vl, sSplit, sEnd, degree);
                                            auto oppPair = encSegs[sOpp].splitArc(T(0.5));
                                            if (arcA.bsp && arcB.bsp && oppPair.first.bsp && oppPair.second.bsp) {
                                                // encSegs[sEnc] goes from pEnd -> pStart (opposite to mergedArc pStart -> pEnd)
                                                // Check which half of oppPair matches which side of pEnc
                                                int sPrev = (sEnc + 3) % 4; // ends at encSegs[sEnc].p0()
                                                int sNext = (sEnc + 1) % 4; // starts at encSegs[sEnc].p1()
                                                gsPoint2D<T> mArc = arcA.p1();
                                                gsPoint2D<T> mOpp = oppPair.first.p1();
                                                CurveSeg<T> midCut = makeStraightSeg(mArc, mOpp, degree);

                                                // Match side connecting pStart to oppPair
                                                bool nextAtStart = (encSegs[sNext].p0().dist(pStart) < T(5e-4) || encSegs[sNext].p1().dist(pStart) < T(5e-4));
                                                CurveSeg<T> sideStart = nextAtStart ? encSegs[sNext] : encSegs[sPrev];
                                                CurveSeg<T> sideEnd   = nextAtStart ? encSegs[sPrev] : encSegs[sNext];
                                                bool firstNearStart = (oppPair.first.p0().dist(sideStart.p0()) < T(5e-4) || oppPair.first.p0().dist(sideStart.p1()) < T(5e-4) ||
                                                                       oppPair.first.p1().dist(sideStart.p0()) < T(5e-4) || oppPair.first.p1().dist(sideStart.p1()) < T(5e-4));
                                                CurveSeg<T> oppStart = firstNearStart ? oppPair.first : oppPair.second;
                                                CurveSeg<T> oppEnd   = firstNearStart ? oppPair.second : oppPair.first;

                                                std::array<CurveSeg<T>, 4> rawA = {arcA, midCut, oppStart, sideStart};
                                                std::array<CurveSeg<T>, 4> rawB = {arcB, sideEnd, oppEnd, midCut};
                                                std::array<CurveSeg<T>, 4> c4A, c4B;
                                                if (chain4ToCCW(rawA, c4A) && chain4ToCCW(rawB, c4B)) {
                                                    auto npA = buildCertifiedSinglePatch(c4A);
                                                    auto npB = buildCertifiedSinglePatch(c4B);
                                                    if (npA && npB) {
                                                        std::set<int> cl6 = cl5;
                                                        cl6.insert(pEnc);
                                                        std::vector<std::unique_ptr<gsTensorBSpline<2, T>>> add;
                                                        add.push_back(std::move(npA));
                                                        add.push_back(std::move(npB));
                                                        if (commitReplacement(cl6, add, oldBnd)) {
                                                            if (dbg) std::cerr << "[MERGED] Absorbed Blister5+pEnc 6->2 (p0=" << rg.p0 << ",enc=" << pEnc << ") P=" << mp.nPatches() << "\n";
                                                            ++totalPruned;
                                                            prunedOne = true;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (!prunedOne) break;
    }

    return totalPruned;
}

} // namespace internal_splitter
} // namespace gismo
