/** @file gsPatchSplitterTemplates.hpp
    @brief Structured topological decomposition templates for gsPatchSplitter:
           Sharp-corner block rings, cross-8 quad rings, bilinear rescue rings,
           proximity rings, and 4-sided multi-hole sliced grids.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitterGeom.hpp>
#include <gsModeling/gsPatchSplitterTriMesh.hpp>

#include <vector>
#include <array>
#include <functional>
#include <queue>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <chrono>

namespace gismo {
namespace internal_splitter {

// Shared helper returning the exterior CCW tangent cone [aStart, aStart + dCCW] at the junction prevSeg -> nextSeg.
template <class T>
inline bool computeExteriorCone(
    const CurveSeg<T>& prevSeg,
    const CurveSeg<T>& nextSeg,
    const gsPoint2D<T>& dir,
    T& aStart, T& dCCW, T& aDir, T& minEdgeLen)
{
    auto dIn  = evalCurveDeriv(*prevSeg.bsp, T(1));
    auto dOut = evalCurveDeriv(*nextSeg.bsp, T(0));
    T nIn = std::hypot(dIn.x, dIn.y), nOut = std::hypot(dOut.x, dOut.y);
    T lenD = std::hypot(dir.x, dir.y);
    if (nIn <= T(1e-14) || nOut <= T(1e-14) || lenD <= T(1e-14))
        return false;
    const T twoPi = T(2) * T(3.14159265358979323846);
    gsPoint2D<T> tBack(-dIn.x / nIn, -dIn.y / nIn);
    gsPoint2D<T> tOut(dOut.x / nOut, dOut.y / nOut);
    aStart = std::atan2(tBack.y, tBack.x);
    dCCW = std::atan2(tOut.y, tOut.x) - aStart;
    while (dCCW <= T(0)) dCCW += twoPi;
    while (dCCW > twoPi) dCCW -= twoPi;
    aDir = std::atan2(dir.y / lenD, dir.x / lenD) - aStart;
    while (aDir < T(0)) aDir += twoPi;
    while (aDir >= twoPi) aDir -= twoPi;
    T lPrev = std::hypot(prevSeg.p1().x - prevSeg.p0().x, prevSeg.p1().y - prevSeg.p0().y);
    T lNext = std::hypot(nextSeg.p1().x - nextSeg.p0().x, nextSeg.p1().y - nextSeg.p0().y);
    minEdgeLen = std::min(lPrev, lNext);
    return true;
}


// Build a spoke CurveSeg H -> P; curves gently into the exterior safe cone only when a
// straight spoke grazes a hole tangent (< 14 deg) or creates a > 172 deg angle at a sharp corner.
template <class T>
inline CurveSeg<T> makeExteriorConeSafeSpoke(
    const CurveSeg<T>& prevSeg,
    const CurveSeg<T>& nextSeg,
    const gsPoint2D<T>& H,
    const gsPoint2D<T>& P,
    const std::vector<gsPoint2D<T>>& holePoly,
    int degree,
    bool& wasCurved)
{
    wasCurved = false;
    auto mkLin = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        std::unique_ptr<gsBSpline<T>> ptr(makeStraightBSpline(a, b, degree));
        return CurveSeg<T>(*ptr);
    };
    gsPoint2D<T> diff(P.x - H.x, P.y - H.y);
    T lenD = std::hypot(diff.x, diff.y);
    T aStart = 0, dCCW = 0, aDir = 0, minEdgeLen = 0;
    if (lenD > T(1e-14) && computeExteriorCone(prevSeg, nextSeg, diff, aStart, dCCW, aDir, minEdgeLen)) {
        const T pi = T(3.14159265358979323846);
        const T trigMargin = T(0.32);
        if (aDir < trigMargin || aDir > dCCW - trigMargin || aDir > pi - trigMargin || (dCCW - aDir) > pi - trigMargin) {
            const T targetMargin = T(0.36);
            T loSafe = std::max(targetMargin, dCCW - pi + targetMargin);
            T hiSafe = std::min(pi - targetMargin, dCCW - targetMargin);
            T targetRel = (loSafe < hiSafe) ? std::max(loSafe, std::min(hiSafe, aDir)) : (T(0.5) * dCCW);
            T absAng = aStart + targetRel;
            gsPoint2D<T> tanDir(std::cos(absAng), std::sin(absAng));
            // Smooth convex (zero-inflection) cubic Bezier: compact c1 along safe tangent, c2 = midpoint(c1, P).
            T tLen = std::min(T(0.18) * lenD, std::max(T(0.06) * lenD, T(0.85) * minEdgeLen));
            gsPoint2D<T> c1(H.x + tLen * tanDir.x, H.y + tLen * tanDir.y);
            gsPoint2D<T> c2(T(0.5) * (c1.x + P.x), T(0.5) * (c1.y + P.y));
            gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
            gsMatrix<T> cf(4, 2);
            cf << H.x, H.y, c1.x, c1.y, c2.x, c2.y, P.x, P.y;
            gsBSpline<T> bsp(kv, cf);
            if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
            bool curveClear = true;
            for (int k = 1; k < 16; ++k) {
                gsPoint2D<T> q = evalCurvePoint(bsp, static_cast<T>(k) / T(16));
                if (pointInPolyStrict(q, holePoly, T(1e-7))) { curveClear = false; break; }
            }
            if (curveClear) {
                wasCurved = true;
                return CurveSeg<T>(bsp);
            }
        }
    }
    return mkLin(H, P);
}

// Adjust an interior cap point relative to hole vertex vm so the segment vm -> pt
// stays inside the safe exterior convex cone (both quad corner angles at vm < 166 deg).
template <class T>
inline gsPoint2D<T> clampPointToExteriorCone(
    const CurveSeg<T>& prevSeg,
    const CurveSeg<T>& nextSeg,
    const gsPoint2D<T>& vm,
    const gsPoint2D<T>& cand)
{
    gsPoint2D<T> diff(cand.x - vm.x, cand.y - vm.y);
    T lenD = std::hypot(diff.x, diff.y);
    T aStart = 0, dCCW = 0, aDir = 0, minEdgeLen = 0;
    if (lenD <= T(1e-14) || !computeExteriorCone(prevSeg, nextSeg, diff, aStart, dCCW, aDir, minEdgeLen))
        return cand;
    const T pi = T(3.14159265358979323846);
    T loSafe = std::max(T(0.32), dCCW - pi + T(0.32));
    T hiSafe = std::min(pi - T(0.32), dCCW - T(0.32));
    T targetRel = (loSafe < hiSafe) ? std::max(loSafe, std::min(hiSafe, aDir)) : (T(0.5) * dCCW);
    if (std::abs(targetRel - aDir) < T(1e-6))
        return cand;
    T absAng = aStart + targetRel;
    return gsPoint2D<T>(vm.x + lenD * std::cos(absAng), vm.y + lenD * std::sin(absAng));
}

// Build a sharp-corner-preserving conforming quad ring inside a 4-sided block [cS, cE, cN, cW]
// for an arbitrary even-segment CCW hole (or K=3 split to 4).
// All hole corners become patch vertices; outer block edges cS, cE, cN, cW remain unsplit.
template <class T>
inline bool buildSharpCornerBlockRing(
    const CurveSeg<T>& cS, const CurveSeg<T>& cE,
    const CurveSeg<T>& cN, const CurveSeg<T>& cW,
    const std::vector<CurveSeg<T>>& hCCW,
    const int anchorIdx[4],
    const std::vector<gsPoint2D<T>>& holePoly,
    std::vector<std::unique_ptr<gsGeometry<T>>>& outPatches,
    int degree,
    T outWeight = T(0.62),
    std::string* failReason = nullptr,
    T spokeSplitT = T(0.5))
{
    const int K = static_cast<int>(hCCW.size());
    // Compute sectors: sector s runs from anchorIdx[s] to anchorIdx[(s+1)%4]
    std::vector<std::vector<CurveSeg<T>>> sectors(4);
    for (int s = 0; s < 4; ++s) {
        int i0 = anchorIdx[s];
        int i1 = anchorIdx[(s + 1) % 4];
        int cur = i0;
        for (int step = 0; step < K; ++step) {
            sectors[s].push_back(hCCW[cur]);
            cur = (cur + 1) % K;
            if (cur == i1) break;
        }
        if (sectors[s].empty()) return false;
    }

    const gsPoint2D<T> P[4] = { cS.p0(), cS.p1(), cE.p1(), cN.p0() }; // SW, SE, NE, NW
    const gsPoint2D<T> H[4] = {
        hCCW[anchorIdx[0]].p0(), hCCW[anchorIdx[1]].p0(),
        hCCW[anchorIdx[2]].p0(), hCCW[anchorIdx[3]].p0()
    };

    // Quick visibility check on anchor spokes
    for (int s = 0; s < 4; ++s) {
        if (!isHoleSpokeClearOfHole(H[s], P[s], holePoly)) {
            if (failReason) *failReason = "spokeBlocked_s" + std::to_string(s);
            return false;
        }
    }

    int totalRem = 0;
    for (int s = 0; s < 4; ++s) totalRem += static_cast<int>(sectors[s].size()) - 1;
    if (totalRem % 2 != 0) return false;

    // Helper: create a straight spoke CurveSeg A -> B
    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    struct QuadGraph {
        std::vector<gsPoint2D<T>> pts;
        std::vector<bool> pinned;
        struct QItem {
            int vid[4];
            CurveSeg<T> e[4]; // cycle: v0-(e0)->v1-(e1)->v2-(e2)->v3-(e3)->v0
            bool straight[4];
        };
        std::vector<QItem> items;

        int addVertex(const gsPoint2D<T>& p, bool isPinned) {
            for (size_t i = 0; i < pts.size(); ++i) {
                T dx = pts[i].x - p.x, dy = pts[i].y - p.y;
                if (dx * dx + dy * dy < T(1e-20)) return static_cast<int>(i);
            }
            pts.push_back(p);
            pinned.push_back(isPinned);
            return static_cast<int>(pts.size()) - 1;
        }
        int addFreeVertex(const gsPoint2D<T>& p) {
            pts.push_back(p);
            pinned.push_back(false);
            return static_cast<int>(pts.size()) - 1;
        }
        void addQuad(const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                     const CurveSeg<T>& e23, const CurveSeg<T>& e30,
                     int v0 = -1, int v1 = -1, int v2 = -1, int v3 = -1) {
            QItem q;
            q.vid[0] = (v0 >= 0) ? v0 : addVertex(e01.p0(), true);
            q.vid[1] = (v1 >= 0) ? v1 : addVertex(e12.p0(), true);
            q.vid[2] = (v2 >= 0) ? v2 : addVertex(e23.p0(), true);
            q.vid[3] = (v3 >= 0) ? v3 : addVertex(e30.p0(), true);
            q.e[0] = e01; q.e[1] = e12; q.e[2] = e23; q.e[3] = e30;
            for (int m = 0; m < 4; ++m)
                q.straight[m] = (!pinned[q.vid[m]] || !pinned[q.vid[(m + 1) % 4]]);
            items.push_back(q);
        }
    } qg;

    // Quad helper storing graph record
    auto addCCWQuad = [&](const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                          const CurveSeg<T>& e23, const CurveSeg<T>& e30,
                          int v0 = -1, int v1 = -1, int v2 = -1, int v3 = -1) {
        qg.addQuad(e01, e12, e23, e30, v0, v1, v2, v3);
    };

    auto projU = [](const gsPoint2D<T>& q, const gsPoint2D<T>& Tl, const gsPoint2D<T>& Tr, T uMin, T uMax) -> T {
        T dx = Tr.x - Tl.x, dy = Tr.y - Tl.y;
        T len2 = dx * dx + dy * dy;
        if (len2 <= T(1e-20)) return T(0.5) * (uMin + uMax);
        T u = ((q.x - Tl.x) * dx + (q.y - Tl.y) * dy) / len2;
        return std::max(uMin, std::min(uMax, u));
    };
    auto lerpPt = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b, T u) -> gsPoint2D<T> {
        return gsPoint2D<T>((T(1) - u) * a.x + u * b.x, (T(1) - u) * a.y + u * b.y);
    };
    CurveSeg<T> outerTop[4] = { cS, cE, cN.reversed(), cW.reversed() };
    std::vector<gsPoint2D<T>> blockPoly;
    blockPoly.reserve(32);
    for (int s = 0; s < 4; ++s) {
        for (int k = 0; k < 8; ++k)
            blockPoly.push_back(evalCurvePoint(*outerTop[s].bsp, static_cast<T>(k) / T(8)));
    }
    CurveSeg<T> sp[4];
    bool spCurved[4] = {false, false, false, false};
    int anchorVid[4] = {-1, -1, -1, -1};
    for (int s = 0; s < 4; ++s) {
        const int idx = anchorIdx[s];
        const auto& prevSeg = hCCW[(idx + K - 1) % K];
        const auto& nextSeg = hCCW[idx];
        sp[s] = makeExteriorConeSafeSpoke(prevSeg, nextSeg, H[s], P[s], holePoly, degree, spCurved[s]);
        anchorVid[s] = qg.addVertex(H[s], true);
    }

    struct SectorFrame {
        std::vector<CurveSeg<T>> holeEdges; // size 1 or 2 after Step A
        CurveSeg<T> topEdge;                // T_left -> T_right
        CurveSeg<T> leftSpoke;              // B_left -> T_left (outward)
        CurveSeg<T> rightSpoke;             // B_right -> T_right (outward)
        int vidTl = -1, vidTr = -1;
    };

    auto placeCapMidpoint = [&](const SectorFrame& f, gsPoint2D<T>& Qout) -> bool {
        const gsPoint2D<T> Tl = qg.pts[f.vidTl], Tr = qg.pts[f.vidTr];
        const CurveSeg<T>& segIn  = f.holeEdges[0];
        const CurveSeg<T>& segOut = f.holeEdges[1];
        gsPoint2D<T> vm = segIn.p1();
        T u = projU(vm, Tl, Tr, T(0.20), T(0.80));
        gsPoint2D<T> T_out = evalCurvePoint(*f.topEdge.bsp, u);
        Qout = gsPoint2D<T>(outWeight * T_out.x + (T(1) - outWeight) * vm.x,
                            outWeight * T_out.y + (T(1) - outWeight) * vm.y);
        Qout = clampPointToExteriorCone(segIn, segOut, vm, Qout);
        if (pointInPolyStrict(Qout, holePoly, T(1e-7)))
            return false;
        if (!isHoleSpokeClearOfHole(vm, Qout, holePoly))
            return false;
        return true;
    };

    SectorFrame fr[4];
    for (int s = 0; s < 4; ++s) {
        fr[s].holeEdges = sectors[s];
        fr[s].topEdge = outerTop[s];
        fr[s].leftSpoke = sp[s];
        fr[s].rightSpoke = sp[(s + 1) % 4];
        fr[s].vidTl = qg.addVertex(fr[s].topEdge.p0(), true);
        fr[s].vidTr = qg.addVertex(fr[s].topEdge.p1(), true);
    }

    // Track sectors whose size==2 needs half-sector closure at a pre-split spoke midpoint M.
    struct PendingTwo {
        bool active = false;
        bool mOnRight = false;
        int vM = -1;
        CurveSeg<T> spLo, spHi;
    } pend2[4];

    auto peelSize4RightM = [&](int s, int vM, const gsPoint2D<T>& M,
                               const CurveSeg<T>& spLo, const CurveSeg<T>& spHi) -> bool {
        gsPoint2D<T> Tl = qg.pts[fr[s].vidTl], Tr = qg.pts[fr[s].vidTr];
        gsPoint2D<T> v0 = fr[s].holeEdges[0].p0(), v1 = fr[s].holeEdges[0].p1();
        T u1 = projU(v1, Tl, Tr, T(0.18), T(0.55));
        gsPoint2D<T> T1 = evalCurvePoint(*fr[s].topEdge.bsp, u1);
        gsPoint2D<T> p1(outWeight * T1.x + (T(1) - outWeight) * v1.x,
                        outWeight * T1.y + (T(1) - outWeight) * v1.y);
        if (fr[s].holeEdges.size() >= 2) {
            gsPoint2D<T> p1c = clampPointToExteriorCone(fr[s].holeEdges[0], fr[s].holeEdges[1], v1, p1);
            if (!pointInPolyStrict(p1c, holePoly, T(1e-7)) &&
                pointInPolyStrict(p1c, blockPoly, T(1e-7)) &&
                isHoleSpokeClearOfHole(v1, p1c, holePoly))
                p1 = p1c;
        }
        if (pointInPolyStrict(p1, holePoly, T(1e-7)) ||
            !pointInPolyStrict(p1, blockPoly, T(1e-7)) ||
            !isHoleSpokeClearOfHole(v1, p1, holePoly) ||
            !isHoleSpokeClearOfHole(p1, M, holePoly))
            return false;
        int vp1 = qg.addFreeVertex(p1);
        int vv0 = qg.addVertex(v0, true), vv1 = qg.addVertex(v1, true);
        CurveSeg<T> e_p1M = mkLine(p1, M), e_v1p1 = mkLine(v1, p1), e_p1Tl = mkLine(p1, Tl);
        addCCWQuad(e_p1M, spHi, fr[s].topEdge.reversed(), e_p1Tl.reversed(), vp1, vM, fr[s].vidTr, fr[s].vidTl);
        addCCWQuad(fr[s].holeEdges[0], e_v1p1, e_p1Tl, fr[s].leftSpoke.reversed(), vv0, vv1, vp1, fr[s].vidTl);
        fr[s].holeEdges.erase(fr[s].holeEdges.begin());
        fr[s].topEdge = e_p1M;
        fr[s].leftSpoke = e_v1p1;
        fr[s].rightSpoke = spLo;
        fr[s].vidTl = vp1;
        fr[s].vidTr = vM;
        return true;
    };

    auto peelSize4LeftM = [&](int s, int vM, const gsPoint2D<T>& M,
                              const CurveSeg<T>& spLo, const CurveSeg<T>& spHi) -> bool {
        gsPoint2D<T> Tl = qg.pts[fr[s].vidTl], Tr = qg.pts[fr[s].vidTr];
        gsPoint2D<T> v3 = fr[s].holeEdges.back().p0(), v4 = fr[s].holeEdges.back().p1();
        T u2 = projU(v3, Tl, Tr, T(0.45), T(0.82));
        gsPoint2D<T> T2 = evalCurvePoint(*fr[s].topEdge.bsp, u2);
        gsPoint2D<T> p2(outWeight * T2.x + (T(1) - outWeight) * v3.x,
                        outWeight * T2.y + (T(1) - outWeight) * v3.y);
        if (fr[s].holeEdges.size() >= 2) {
            const size_t nh = fr[s].holeEdges.size();
            gsPoint2D<T> p2c = clampPointToExteriorCone(fr[s].holeEdges[nh - 2], fr[s].holeEdges.back(), v3, p2);
            if (!pointInPolyStrict(p2c, holePoly, T(1e-7)) &&
                pointInPolyStrict(p2c, blockPoly, T(1e-7)) &&
                isHoleSpokeClearOfHole(v3, p2c, holePoly))
                p2 = p2c;
        }
        if (pointInPolyStrict(p2, holePoly, T(1e-7)) ||
            !pointInPolyStrict(p2, blockPoly, T(1e-7)) ||
            !isHoleSpokeClearOfHole(v3, p2, holePoly) ||
            !isHoleSpokeClearOfHole(M, p2, holePoly))
            return false;
        int vp2 = qg.addFreeVertex(p2);
        int vv3 = qg.addVertex(v3, true), vv4 = qg.addVertex(v4, true);
        CurveSeg<T> e_Mp2 = mkLine(M, p2), e_v3p2 = mkLine(v3, p2), e_p2Tr = mkLine(p2, Tr);
        addCCWQuad(e_Mp2, e_p2Tr, fr[s].topEdge.reversed(), spHi.reversed(), vM, vp2, fr[s].vidTr, fr[s].vidTl);
        addCCWQuad(fr[s].holeEdges.back(), fr[s].rightSpoke, e_p2Tr.reversed(), e_v3p2.reversed(), vv3, vv4, fr[s].vidTr, vp2);
        fr[s].holeEdges.pop_back();
        fr[s].topEdge = e_Mp2;
        fr[s].leftSpoke = spLo;
        fr[s].rightSpoke = e_v3p2;
        fr[s].vidTl = vM;
        fr[s].vidTr = vp2;
        return true;
    };

    // Pre-pair all even sectors (size 2 or 4) BEFORE Step A:
    // - Adjacent pair (a, b=(a+1)%4): split shared spoke sp[b] at M
    // - Opposite pair (a, b=(a+2)%4): split bridge oppMid=(a+1)%4 spokes at M1, M2 and emit bridge top strip
    // - Any size==4 sector in a pair peels 1 edge into its split M, becoming size==3 for Step A!
    {
        const T tSplit = std::max(T(0.15), std::min(T(0.85), spokeSplitT));
        std::vector<int> evenSec;
        for (int s = 0; s < 4; ++s)
            if (fr[s].holeEdges.size() % 2 == 0) evenSec.push_back(s);

        auto bindRightLeft = [&](int u, int v, int vM, const gsPoint2D<T>& M,
                                 const CurveSeg<T>& spLo, const CurveSeg<T>& spHi) -> bool {
            if (fr[u].holeEdges.size() == 4) {
                if (!peelSize4RightM(u, vM, M, spLo, spHi)) return false;
            } else {
                pend2[u] = {true, true, vM, spLo, spHi};
            }
            if (fr[v].holeEdges.size() == 4) {
                if (!peelSize4LeftM(v, vM, M, spLo, spHi)) return false;
            } else {
                pend2[v] = {true, false, vM, spLo, spHi};
            }
            return true;
        };

        auto pairAdjacent = [&](int a, int b) -> bool {
            gsBSpline<T> p0, p1;
            splitBSplineExact(*fr[a].rightSpoke.bsp, tSplit, p0, p1);
            CurveSeg<T> spLo(p0), spHi(p1);
            gsPoint2D<T> M = spLo.p1();
            int spIdx = (a + 1) % 4;
            int vM = spCurved[spIdx] ? qg.addVertex(M, true) : qg.addFreeVertex(M);
            return bindRightLeft(a, b, vM, M, spLo, spHi);
        };

        if (evenSec.size() == 2) {
            int a = evenSec[0], b = evenSec[1];
            if ((a + 1) % 4 == b) {
                if (!pairAdjacent(a, b)) return false;
            } else if ((b + 1) % 4 == a) {
                if (!pairAdjacent(b, a)) return false;
            } else {
                int mid = (a + 1) % 4;
                gsBSpline<T> am0, am1, mb0, mb1;
                splitBSplineExact(*fr[mid].leftSpoke.bsp, tSplit, am0, am1);
                splitBSplineExact(*fr[mid].rightSpoke.bsp, tSplit, mb0, mb1);
                CurveSeg<T> spAM_lo(am0), spAM_hi(am1), spMB_lo(mb0), spMB_hi(mb1);
                gsPoint2D<T> M1 = spAM_lo.p1(), M2 = spMB_lo.p1();
                int vM1 = spCurved[mid] ? qg.addVertex(M1, true) : qg.addFreeVertex(M1);
                int vM2 = spCurved[(mid + 1) % 4] ? qg.addVertex(M2, true) : qg.addFreeVertex(M2);
                CurveSeg<T> e_M1M2 = mkLine(M1, M2);
                addCCWQuad(e_M1M2, spMB_hi, fr[mid].topEdge.reversed(), spAM_hi.reversed(),
                           vM1, vM2, fr[mid].vidTr, fr[mid].vidTl);
                fr[mid].topEdge = e_M1M2;
                fr[mid].leftSpoke = spAM_lo;
                fr[mid].rightSpoke = spMB_lo;
                fr[mid].vidTl = vM1;
                fr[mid].vidTr = vM2;
                if (fr[a].holeEdges.size() == 4) {
                    if (!peelSize4RightM(a, vM1, M1, spAM_lo, spAM_hi)) return false;
                } else {
                    pend2[a] = {true, true, vM1, spAM_lo, spAM_hi};
                }
                if (fr[b].holeEdges.size() == 4) {
                    if (!peelSize4LeftM(b, vM2, M2, spMB_lo, spMB_hi)) return false;
                } else {
                    pend2[b] = {true, false, vM2, spMB_lo, spMB_hi};
                }
            }
        } else if (evenSec.size() == 4) {
            if (!pairAdjacent(0, 1) || !pairAdjacent(2, 3))
                return false;
        }
    }

    // Step A: Recursively apply 3->1 pillow cap in any sector with >= 3 hole edges until every sector has 1 or 2 edges.
    for (int s = 0; s < 4; ++s) {
        while (fr[s].holeEdges.size() >= 3) {
            gsPoint2D<T> Tl = qg.pts[fr[s].vidTl];
            gsPoint2D<T> Tr = qg.pts[fr[s].vidTr];
            gsPoint2D<T> v1 = fr[s].holeEdges.front().p1();
            gsPoint2D<T> v2 = fr[s].holeEdges.back().p0();

            // Place p1, p2 on a chord parallel to Tl--Tr strictly in front of both v1 and v2;
            // fall back to curve-point interpolation if topEdge curves past the straight chord (e.g., bent starburst rays).
            T dxT = Tr.x - Tl.x, dyT = Tr.y - Tl.y;
            T lenT = std::hypot(dxT, dyT);
            if (lenT <= T(1e-14)) return false;
            gsPoint2D<T> nIn(-dyT / lenT, dxT / lenT); // CCW interior normal (domain on left of Tl->Tr)
            T h1 = (v1.x - Tl.x) * nIn.x + (v1.y - Tl.y) * nIn.y;
            T h2 = (v2.x - Tl.x) * nIn.x + (v2.y - Tl.y) * nIn.y;
            T hMin = std::min(h1, h2);
            T u1 = projU(v1, Tl, Tr, T(0.15), T(0.68));
            T u2 = projU(v2, Tl, Tr, std::min(T(0.85), u1 + T(0.12)), T(0.85));
            gsPoint2D<T> p1, p2;
            if (hMin > T(1e-6)) {
                T hCap = (T(1) - outWeight) * hMin;
                gsPoint2D<T> T_left13 = lerpPt(Tl, Tr, u1);
                gsPoint2D<T> T_rght23 = lerpPt(Tl, Tr, u2);
                p1 = gsPoint2D<T>(T_left13.x + hCap * nIn.x, T_left13.y + hCap * nIn.y);
                p2 = gsPoint2D<T>(T_rght23.x + hCap * nIn.x, T_rght23.y + hCap * nIn.y);
            } else {
                gsPoint2D<T> T1 = evalCurvePoint(*fr[s].topEdge.bsp, u1);
                gsPoint2D<T> T2 = evalCurvePoint(*fr[s].topEdge.bsp, u2);
                T wC = T(0.65) * outWeight;
                p1 = gsPoint2D<T>((T(1) - wC) * v1.x + wC * T1.x, (T(1) - wC) * v1.y + wC * T1.y);
                p2 = gsPoint2D<T>((T(1) - wC) * v2.x + wC * T2.x, (T(1) - wC) * v2.y + wC * T2.y);
            }
            {
                const size_t nh = fr[s].holeEdges.size();
                gsPoint2D<T> p1c = clampPointToExteriorCone(fr[s].holeEdges[0], fr[s].holeEdges[1], v1, p1);
                gsPoint2D<T> p2c = clampPointToExteriorCone(fr[s].holeEdges[nh - 2], fr[s].holeEdges.back(), v2, p2);
                if (!pointInPolyStrict(p1c, holePoly, T(1e-7)) &&
                    pointInPolyStrict(p1c, blockPoly, T(1e-7)) &&
                    isHoleSpokeClearOfHole(v1, p1c, holePoly))
                    p1 = p1c;
                if (!pointInPolyStrict(p2c, holePoly, T(1e-7)) &&
                    pointInPolyStrict(p2c, blockPoly, T(1e-7)) &&
                    isHoleSpokeClearOfHole(v2, p2c, holePoly))
                    p2 = p2c;
            }

            if (pointInPolyStrict(p1, holePoly, T(1e-7)) || pointInPolyStrict(p2, holePoly, T(1e-7)) ||
                !pointInPolyStrict(p1, blockPoly, T(1e-7)) || !pointInPolyStrict(p2, blockPoly, T(1e-7)) ||
                !isHoleSpokeClearOfHole(v1, p1, holePoly) || !isHoleSpokeClearOfHole(v2, p2, holePoly) ||
                !isHoleSpokeClearOfHole(p1, p2, holePoly)) {
                if (failReason && hMin <= T(1e-6)) *failReason = "hMin<=0";
                return false;
            }

            int vp1 = qg.addFreeVertex(p1);
            int vp2 = qg.addFreeVertex(p2);
            int vBl = qg.addVertex(fr[s].holeEdges.front().p0(), true);
            int vv1 = qg.addVertex(v1, true);
            int vv2 = qg.addVertex(v2, true);
            int vBr = qg.addVertex(fr[s].holeEdges.back().p1(), true);

            CurveSeg<T> top_p1p2 = mkLine(p1, p2);
            CurveSeg<T> left_v1p1 = mkLine(v1, p1);
            CurveSeg<T> right_v2p2 = mkLine(v2, p2);
            CurveSeg<T> capL = mkLine(p1, Tl);
            CurveSeg<T> capR = mkLine(p2, Tr);

            addCCWQuad(top_p1p2, capR, fr[s].topEdge.reversed(), capL.reversed(), vp1, vp2, fr[s].vidTr, fr[s].vidTl);
            addCCWQuad(fr[s].holeEdges.front(), left_v1p1, capL, fr[s].leftSpoke.reversed(), vBl, vv1, vp1, fr[s].vidTl);
            addCCWQuad(fr[s].holeEdges.back(), fr[s].rightSpoke, capR.reversed(), right_v2p2.reversed(), vv2, vBr, fr[s].vidTr, vp2);

            std::vector<CurveSeg<T>> mid(fr[s].holeEdges.begin() + 1, fr[s].holeEdges.end() - 1);
            fr[s].holeEdges = std::move(mid);
            fr[s].topEdge = top_p1p2;
            fr[s].leftSpoke = left_v1p1;
            fr[s].rightSpoke = right_v2p2;
            fr[s].vidTl = vp1;
            fr[s].vidTr = vp2;
        }
    }

    // Emit remaining sectors cleanly (pend2 sectors close into pre-split M, unhandled sectors have size 1).
    for (int s = 0; s < 4; ++s) {
        if (pend2[s].active) {
            if (fr[s].holeEdges.size() != 2) return false;
            gsPoint2D<T> Q;
            if (!placeCapMidpoint(fr[s], Q)) return false;
            int vQ = qg.addFreeVertex(Q);
            int vBl = qg.addVertex(fr[s].holeEdges[0].p0(), true);
            gsPoint2D<T> vm = fr[s].holeEdges[0].p1();
            int vvm = qg.addVertex(vm, true);
            int vBr = qg.addVertex(fr[s].holeEdges[1].p1(), true);
            int vM = pend2[s].vM;
            gsPoint2D<T> M = qg.pts[vM];
            CurveSeg<T> e_vmQ = mkLine(vm, Q);
            if (pend2[s].mOnRight) {
                CurveSeg<T> e_QTl = mkLine(Q, qg.pts[fr[s].vidTl]), e_QM = mkLine(Q, M);
                addCCWQuad(fr[s].holeEdges[0], e_vmQ, e_QTl, fr[s].leftSpoke.reversed(), vBl, vvm, vQ, fr[s].vidTl);
                addCCWQuad(fr[s].holeEdges[1], pend2[s].spLo, e_QM.reversed(), e_vmQ.reversed(), vvm, vBr, vM, vQ);
                addCCWQuad(e_QM, pend2[s].spHi, fr[s].topEdge.reversed(), e_QTl.reversed(), vQ, vM, fr[s].vidTr, fr[s].vidTl);
            } else {
                CurveSeg<T> e_QTr = mkLine(Q, qg.pts[fr[s].vidTr]), e_MQ = mkLine(M, Q);
                addCCWQuad(fr[s].holeEdges[0], e_vmQ, e_MQ.reversed(), pend2[s].spLo.reversed(), vBl, vvm, vQ, vM);
                addCCWQuad(fr[s].holeEdges[1], fr[s].rightSpoke, e_QTr.reversed(), e_vmQ.reversed(), vvm, vBr, fr[s].vidTr, vQ);
                addCCWQuad(e_MQ, e_QTr, fr[s].topEdge.reversed(), pend2[s].spHi.reversed(), vM, vQ, fr[s].vidTr, fr[s].vidTl);
            }
        } else {
            if (fr[s].holeEdges.size() != 1) return false;
            int vBl = qg.addVertex(fr[s].holeEdges[0].p0(), true);
            int vBr = qg.addVertex(fr[s].holeEdges[0].p1(), true);
            addCCWQuad(fr[s].holeEdges[0], fr[s].rightSpoke, fr[s].topEdge.reversed(), fr[s].leftSpoke.reversed(),
                       vBl, vBr, fr[s].vidTr, fr[s].vidTl);
        }
    }

    // Guarded Laplacian + convexity hill-climb on free interior vertices before Coons generation
    {
        const int NV = static_cast<int>(qg.pts.size());
        std::vector<bool> hasCurvedAnchor(NV, false);
        for (int s = 0; s < 4; ++s) {
            if (spCurved[s] && anchorVid[s] >= 0 && anchorVid[s] < NV)
                hasCurvedAnchor[anchorVid[s]] = true;
        }
        std::vector<std::vector<int>> adj(NV);
        std::vector<std::vector<int>> incQ(NV);
        for (size_t qi_idx = 0; qi_idx < qg.items.size(); ++qi_idx) {
            const auto& qi = qg.items[qi_idx];
            for (int m = 0; m < 4; ++m) {
                adj[qi.vid[m]].push_back(qi.vid[(m + 1) % 4]);
                adj[qi.vid[(m + 1) % 4]].push_back(qi.vid[m]);
                incQ[qi.vid[m]].push_back(static_cast<int>(qi_idx));
            }
        }
        auto minIncSin = [&](int v) -> T {
            T minS = T(1e18);
            for (int qi_idx : incQ[v]) {
                const auto& qi = qg.items[qi_idx];
                gsPoint2D<T> cp[4] = { qg.pts[qi.vid[3]], qg.pts[qi.vid[2]], qg.pts[qi.vid[1]], qg.pts[qi.vid[0]] };
                int vIdx[4] = { qi.vid[3], qi.vid[2], qi.vid[1], qi.vid[0] };
                for (int k = 0; k < 4; ++k) {
                    if (hasCurvedAnchor[vIdx[(k + 1) % 4]]) continue;
                    const auto& A = cp[k];
                    const auto& B = cp[(k + 1) % 4];
                    const auto& C = cp[(k + 2) % 4];
                    T cr = cross2(A, B, C);
                    T l1 = std::hypot(B.x - A.x, B.y - A.y);
                    T l2 = std::hypot(C.x - B.x, C.y - B.y);
                    T sA = (l1 > T(1e-12) && l2 > T(1e-12)) ? (cr / (l1 * l2)) : T(-1);
                    minS = std::min(minS, sA);
                }
            }
            return minS;
        };
        auto allQuadsConvex = [&](std::string* qErr = nullptr) -> bool {
            for (size_t idx = 0; idx < qg.items.size(); ++idx) {
                const auto& qi = qg.items[idx];
                gsPoint2D<T> cp[4] = { qg.pts[qi.vid[3]], qg.pts[qi.vid[2]], qg.pts[qi.vid[1]], qg.pts[qi.vid[0]] };
                int vIdx[4] = { qi.vid[3], qi.vid[2], qi.vid[1], qi.vid[0] };
                for (int k = 0; k < 4; ++k) {
                    const auto& A = cp[k];
                    const auto& B = cp[(k + 1) % 4];
                    const auto& C = cp[(k + 2) % 4];
                    T cr = cross2(A, B, C);
                    T l1 = std::hypot(B.x - A.x, B.y - A.y);
                    T l2 = std::hypot(C.x - B.x, C.y - B.y);
                    T sinA = (l1 > T(1e-12) && l2 > T(1e-12)) ? (cr / (l1 * l2)) : T(0);
                    // Only relax straight-chord angle threshold at hole anchors equipped with a curved bisector spoke
                    // so evaluateRingPatches inspects the actual curved Coons Jacobian.
                    T minReq = hasCurvedAnchor[vIdx[(k + 1) % 4]] ? T(-0.60) : T(1e-3);
                    if (sinA <= minReq) {
                        if (qErr) {
                            std::ostringstream oss;
                            oss << "q=" << idx << "/k=" << k << " sin=" << sinA
                                << " [(" << cp[0].x << "," << cp[0].y << "),("
                                << cp[1].x << "," << cp[1].y << "),("
                                << cp[2].x << "," << cp[2].y << "),("
                                << cp[3].x << "," << cp[3].y << ")]";
                            *qErr = oss.str();
                        }
                        return false;
                    }
                }
            }
            return true;
        };
        auto isFreeCandClear = [&](int v, const gsPoint2D<T>& cand) -> bool {
            if (pointInPolyStrict(cand, holePoly, T(1e-7))) return false;
            if (!pointInPolyStrict(cand, blockPoly, T(1e-7))) return false;
            for (int nb : adj[v]) {
                if (!isHoleSpokeClearOfHole(qg.pts[nb], cand, holePoly)) return false;
            }
            return true;
        };
        for (int iter = 0; iter < 8; ++iter) {
            bool anyMoved = false;
            for (int v = 0; v < NV; ++v) {
                if (qg.pinned[v] || adj[v].empty()) continue;
                gsPoint2D<T> avg(0, 0);
                for (int nb : adj[v]) { avg.x += qg.pts[nb].x; avg.y += qg.pts[nb].y; }
                avg.x /= adj[v].size(); avg.y /= adj[v].size();
                gsPoint2D<T> oldP = qg.pts[v];
                if (oldP.dist(avg) < T(1e-7)) continue;
                T oldMinS = minIncSin(v);
                bool moved = false;
                for (T alpha : {T(1.0), T(0.5)}) {
                    gsPoint2D<T> cand((T(1) - alpha) * oldP.x + alpha * avg.x,
                                      (T(1) - alpha) * oldP.y + alpha * avg.y);
                    qg.pts[v] = cand;
                    T newMinS = minIncSin(v);
                    if (((oldMinS > T(0.05) && newMinS > T(0.05)) || (newMinS > oldMinS + T(1e-4))) &&
                        isFreeCandClear(v, cand))
                    {
                        moved = true;
                        anyMoved = true;
                        break;
                    }
                }
                if (!moved) qg.pts[v] = oldP;
            }
            if (!anyMoved) break;
        }
        // Targeted directional + coordinate convexity/angle rescue if any free vertex touches a shallow/concave corner
        const T dirX[8] = { T(1), T(-1), T(0), T(0), T(0.70710678), T(0.70710678), T(-0.70710678), T(-0.70710678) };
        const T dirY[8] = { T(0), T(0), T(1), T(-1), T(0.70710678), T(-0.70710678), T(0.70710678), T(-0.70710678) };
        for (int pass = 0; pass < 2; ++pass) {
            bool anyImproved = false;
            for (int v = 0; v < NV; ++v) {
                if (qg.pinned[v] || adj[v].empty()) continue;
                T curMin = minIncSin(v);
                if (curMin >= T(0.26)) continue;
                gsPoint2D<T> oldP = qg.pts[v];
                gsPoint2D<T> bestP = oldP;
                T bestS = curMin;
                T sumLen = T(0);
                for (int nb : adj[v]) {
                    sumLen += std::hypot(qg.pts[nb].x - oldP.x, qg.pts[nb].y - oldP.y);
                    for (T frac : {T(0.20), T(0.40)}) {
                        gsPoint2D<T> cand((T(1) - frac) * oldP.x + frac * qg.pts[nb].x,
                                          (T(1) - frac) * oldP.y + frac * qg.pts[nb].y);
                        qg.pts[v] = cand;
                        T s = minIncSin(v);
                        if (s > bestS + T(1e-4) && isFreeCandClear(v, cand)) { bestS = s; bestP = cand; }
                    }
                }
                T avgLen = sumLen / static_cast<T>(adj[v].size());
                if (avgLen > T(1e-12)) {
                    for (T sf : {T(0.14), T(0.05)}) {
                        T st = sf * avgLen;
                        for (int d = 0; d < 8; ++d) {
                            gsPoint2D<T> cand(oldP.x + st * dirX[d], oldP.y + st * dirY[d]);
                            qg.pts[v] = cand;
                            T s = minIncSin(v);
                            if (s > bestS + T(1e-4) && isFreeCandClear(v, cand)) { bestS = s; bestP = cand; }
                        }
                    }
                }
                if (bestS > curMin + T(1e-4)) anyImproved = true;
                qg.pts[v] = bestP;
            }
            if (!anyImproved) break;
        }
        if (!allQuadsConvex(failReason)) return false;
    }

    // Regenerate straight edges connected to relaxed vertices and emit Coons patches
    for (auto& qi : qg.items) {
        for (int m = 0; m < 4; ++m) {
            if (qi.straight[m])
                qi.e[m] = mkLine(qg.pts[qi.vid[m]], qg.pts[qi.vid[(m + 1) % 4]]);
        }
        auto r0 = qi.e[3].reversed(), r1 = qi.e[2].reversed();
        outPatches.emplace_back(makeCoonsPatchFrom4Curves(*r0.bsp, *r1.bsp, *qi.e[1].bsp, *qi.e[0].bsp, degree));
    }
    return true;
}


// Compute minimum positive distance t > 0 from origin along unit direction `dir` to polygon boundary.
template <class T>
inline T rayDistToPoly(
    const gsPoint2D<T>& origin,
    const gsPoint2D<T>& dir,
    const std::vector<gsPoint2D<T>>& poly)
{
    T bestT = T(1e18);
    const int N = static_cast<int>(poly.size());
    for (int i = 0; i < N; ++i) {
        const gsPoint2D<T>& A = poly[i];
        const gsPoint2D<T>& B = poly[(i + 1) % N];
        T ax = A.x - origin.x, ay = A.y - origin.y;
        T sx = B.x - A.x, sy = B.y - A.y;
        T denom = dir.x * sy - dir.y * sx;
        if (std::abs(denom) <= T(1e-14)) continue;
        T t = (ax * sy - ay * sx) / denom;
        T u = (ax * dir.y - ay * dir.x) / denom;
        if (t > T(1e-9) && u >= -T(1e-6) && u <= T(1) + T(1e-6))
            bestT = std::min(bestT, t);
    }
    return bestT;
}

// Blend interior CPs toward Greville-Harmonic (first) and corner-bilinear map when a convex Coons patch overshoots detJ.
template <class T>
inline bool bilinearBlendRescueRing(
    std::vector<std::unique_ptr<gsGeometry<T>>>& cand,
    const std::vector<gsPoint2D<T>>& holePoly,
    T& sc, T& md, T& ma,
    std::string* failMsg = nullptr)
{
    bool fixedAny = false;
    for (auto& gp : cand) {
        T mn = 0, mx = 0, ab = 0, en = 0;
        evalPatchJacobianAndWinslow(*gp, 21, mn, mx, ab, en);
        if (mn <= T(2e-7) && mx >= -T(2e-7)) {
            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(gp.get())) {
                if (quickHarmonicPatchRepair(*tb, T(2e-7))) {
                    fixedAny = true;
                    continue;
                }
                const int nu = tb->basis().component(0).size();
                const int nv = tb->basis().component(1).size();
                if (nu >= 3 && nv >= 3 && mn > T(-0.02)) {
                    gsMatrix<T> origCf = tb->coefs();
                    gsMatrix<T> gu = tb->basis().component(0).knots().greville();
                    gsMatrix<T> gv = tb->basis().component(1).knots().greville();
                    T P00[2] = {origCf(0, 0), origCf(0, 1)};
                    T P10[2] = {origCf(nu - 1, 0), origCf(nu - 1, 1)};
                    T P01[2] = {origCf((nv - 1) * nu, 0), origCf((nv - 1) * nu, 1)};
                    T P11[2] = {origCf(nu - 1 + (nv - 1) * nu, 0), origCf(nu - 1 + (nv - 1) * nu, 1)};
                    bool fixedThis = false;
                    for (T al : {T(0.25), T(0.50), T(0.75), T(1.0)}) {
                        gsMatrix<T> cf = origCf;
                        for (int j = 1; j < nv - 1; ++j) {
                            T eta = gv(0, j);
                            for (int i = 1; i < nu - 1; ++i) {
                                T xi = gu(0, i);
                                for (int d = 0; d < 2; ++d) {
                                    T bil = (1-xi)*(1-eta)*P00[d] + xi*(1-eta)*P10[d]
                                          + (1-xi)*eta*P01[d] + xi*eta*P11[d];
                                    cf(i + j * nu, d) = (1 - al) * origCf(i + j * nu, d) + al * bil;
                                }
                            }
                        }
                        tb->setCoefs(cf);
                        T tmn = 0, tmx = 0, tab = 0, ten = 0;
                        evalPatchJacobianAndWinslow(*tb, 15, tmn, tmx, tab, ten);
                        if (tmn > T(1e-9)) { fixedThis = true; fixedAny = true; break; }
                    }
                    if (!fixedThis) tb->setCoefs(origCf);
                }
            }
        }
    }
    if (!fixedAny) return false;
    return evaluateRingPatches(cand, holePoly, sc, md, ma, failMsg);
}

// Forward declarations for helpers defined below
template <class T>
inline bool verifyPatchSetLocal(
    std::vector<std::unique_ptr<gsGeometry<T>>>& patches,
    const std::vector<std::vector<gsPoint2D<T>>>& holePolys,
    const std::vector<gsPoint2D<T>>* cellPolyPtr,
    T minAngTol,
    T maxAngTol,
    T& totalEnergyOut);

template <class T>
inline bool buildHoleSmoothingCollar(
    const gsCurveLoop<T>& hole,
    const std::vector<gsPoint2D<T>>* enclosingPolyPtr,
    int degree,
    T minAngTol,
    T maxAngTol,
    std::vector<std::unique_ptr<gsGeometry<T>>>& collarPatches,
    std::vector<CurveSeg<T>>& outerCollarCCW,
    std::vector<int>& outerEdgeToPatch);

template <class T>
inline bool splitOneCollarFlankForEvenness(
    std::vector<std::unique_ptr<gsGeometry<T>>>& collarPatches,
    std::vector<CurveSeg<T>>& outerCollarCCW,
    std::vector<int>& outerEdgeToPatch,
    int degree);

template <class T>
inline bool runSharpBlockRingSearch(
    const CurveSeg<T>& bcS, const CurveSeg<T>& bcE,
    const CurveSeg<T>& bcN, const CurveSeg<T>& bcW,
    const gsPoint2D<T> bP[4],
    const std::vector<CurveSeg<T>>& workH,
    const std::vector<gsPoint2D<T>>& holePoly,
    int degree,
    const std::vector<T>& splitTs,
    const std::vector<T>& oWeights,
    bool earlyStop,
    std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
    T& bestScoreOut, T& bestMaxAngOut,
    int* dbgVisOut = nullptr, int* dbgCapOut = nullptr, int* dbgEvalOut = nullptr,
    std::string* firstCapOut = nullptr, std::string* firstFailOut = nullptr)
{
    const int M = static_cast<int>(workH.size());
    std::vector<std::vector<bool>> vis(4, std::vector<bool>(M, false));
    std::vector<std::vector<T>> distMat(4, std::vector<T>(M, T(1e18)));
    for (int s = 0; s < 4; ++s) {
        for (int i = 0; i < M; ++i) {
            vis[s][i] = isHoleSpokeClearOfHole(workH[i].p0(), bP[s], holePoly);
            if (vis[s][i])
                distMat[s][i] = workH[i].p0().dist(bP[s]);
        }
    }

    // Keep only top-4 closest visible anchors per block corner when M > 6 to avoid combinatorial blowup
    if (M > 6) {
        for (int s = 0; s < 4; ++s) {
            std::vector<int> idxs;
            for (int i = 0; i < M; ++i) if (vis[s][i]) idxs.push_back(i);
            if (idxs.size() > 4) {
                std::sort(idxs.begin(), idxs.end(), [&](int i, int j) { return distMat[s][i] < distMat[s][j]; });
                std::vector<bool> keep(M, false);
                for (size_t k = 0; k < 4; ++k) keep[idxs[k]] = true;
                vis[s] = std::move(keep);
            }
        }
    }

    struct AnchorTuple { int a[4]; int maxSec; T alignCost; };
    std::vector<AnchorTuple> tuples;
    const int maxL = (M >= 8) ? 4 : 3;
    for (int a0 = 0; a0 < M; ++a0) {
        if (!vis[0][a0]) continue;
        for (int l0 = 1; l0 <= maxL; ++l0) {
            int a1 = (a0 + l0) % M;
            if (!vis[1][a1]) continue;
            for (int l1 = 1; l1 <= maxL; ++l1) {
                int a2 = (a1 + l1) % M;
                if (!vis[2][a2]) continue;
                for (int l2 = 1; l2 <= maxL; ++l2) {
                    int l3 = M - (l0 + l1 + l2);
                    if (l3 < 1 || l3 > maxL) continue;
                    int a3 = (a2 + l2) % M;
                    if (!vis[3][a3]) continue;
                    int maxSec = std::max(std::max(l0, l1), std::max(l2, l3));
                    T ac = distMat[0][a0] + distMat[1][a1] + distMat[2][a2] + distMat[3][a3];
                    tuples.push_back({{a0, a1, a2, a3}, maxSec, ac});
                }
            }
        }
    }
    std::sort(tuples.begin(), tuples.end(), [](const AnchorTuple& x, const AnchorTuple& y) {
        if (x.maxSec != y.maxSec) return x.maxSec < y.maxSec;
        return x.alignCost < y.alignCost;
    });
    if (tuples.size() > 16) tuples.resize(16);

    bestScoreOut = T(1e18);
    bestMaxAngOut = T(180);
    bool found = false;
    int ringEvals = 0;
    const int maxRingEvals = earlyStop ? 32 : 48;
    for (const auto& tp : tuples) {
        if (ringEvals >= maxRingEvals) break;
        const int* a = tp.a;
        const int maxSec = tp.maxSec;
        if (dbgVisOut) ++(*dbgVisOut);

        for (T st : splitTs) {
            for (T ow : oWeights) {
                if (ringEvals >= maxRingEvals) break;
                ++ringEvals;
                std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                std::string cmsg;
                if (!buildSharpCornerBlockRing(bcS, bcE, bcN, bcW, workH, a, holePoly, cand, degree, ow, &cmsg, st)) {
                    if (firstCapOut && (firstCapOut->empty() || (maxSec == 2 && firstCapOut->find("max2:") == std::string::npos)))
                        *firstCapOut = (maxSec == 2 ? "max2:" : "") + cmsg;
                    if (dbgCapOut) ++(*dbgCapOut);
                    continue;
                }
                T sc = 0, md = 0, ma = 0;
                std::string fmsg;
                bool pass = evaluateRingPatches(cand, holePoly, sc, md, ma, (firstFailOut && firstFailOut->empty()) ? &fmsg : nullptr);
                if (!pass)
                    pass = bilinearBlendRescueRing(cand, holePoly, sc, md, ma, (firstFailOut && firstFailOut->empty()) ? &fmsg : nullptr);
                if (pass) {
                    T totalMetric = sc + T(5) * static_cast<T>(maxSec) + (ma > T(163.5) ? T(1e5) : T(0));
                    if (!found || totalMetric < bestScoreOut) {
                        found = true; bestScoreOut = totalMetric; bestMaxAngOut = ma; bestOut = std::move(cand);
                    }
                    if (earlyStop && bestMaxAngOut <= T(158))
                        return true;
                } else {
                    if (firstFailOut && firstFailOut->empty()) *firstFailOut = fmsg;
                    if (dbgEvalOut) ++(*dbgEvalOut);
                }
            }
        }
    }
    return found;
}

// Build an 8-quad conforming ring inside a 4-sided block [cS, cE, cN, cW] for a K=12 cross hole
// (4 concave ~90 deg notches spaced by 3 edges each: 4 notch quads + 4 arm-tip quads).
template <class T>
inline bool tryBuildCross8QuadRing(
    const CurveSeg<T>& cS, const CurveSeg<T>& cE,
    const CurveSeg<T>& cN, const CurveSeg<T>& cW,
    const std::vector<CurveSeg<T>>& hCCW,
    const std::vector<gsPoint2D<T>>& holePoly,
    std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
    T& bestScoreOut, T& bestMaxAngOut,
    int degree)
{
    if (hCCW.size() != 12) return false;
    const T pi = T(3.14159265358979323846);
    std::vector<int> notches;
    for (int i = 0; i < 12; ++i) {
        const auto& prev = hCCW[(i + 11) % 12];
        const auto& next = hCCW[i];
        gsPoint2D<T> d(next.p1().x - prev.p0().x, next.p1().y - prev.p0().y);
        if (std::hypot(d.x, d.y) <= T(1e-12)) d = gsPoint2D<T>(T(1), T(0));
        T aStart = 0, dCCW = 0, aDir = 0, minEdge = 0;
        if (computeExteriorCone(prev, next, d, aStart, dCCW, aDir, minEdge) &&
            dCCW < (T(150) * pi / T(180)))
        {
            notches.push_back(i);
        }
    }
    if (notches.size() != 4) return false;
    for (int s = 0; s < 4; ++s) {
        if ((notches[(s + 1) % 4] - notches[s] + 12) % 12 != 3)
            return false;
    }

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    const gsPoint2D<T> P[4] = { cS.p0(), cS.p1(), cE.p1(), cN.p0() };
    const CurveSeg<T> outerTop[4] = { cS, cE, cN.reversed(), cW.reversed() };

    bool found = false;
    bestScoreOut = T(1e18);
    bestMaxAngOut = T(180);

    for (int rot = 0; rot < 4; ++rot) {
        CurveSeg<T> eIn[4], eOut[4], eTip[4], spL[4], spR[4];
        bool clearOk = true;
        for (int s = 0; s < 4; ++s) {
            int i0 = notches[(s + rot) % 4];
            eIn[s]  = hCCW[(i0 + 11) % 12]; // T2_{s-1} -> V_s
            eOut[s] = hCCW[i0];             // V_s -> T1_s
            eTip[s] = hCCW[(i0 + 1) % 12];  // T1_s -> T2_s
            if (!isHoleSpokeClearOfHole(eIn[s].p0(),  P[s], holePoly) ||
                !isHoleSpokeClearOfHole(eOut[s].p1(), P[s], holePoly)) {
                clearOk = false;
                break;
            }
            spL[s] = mkLine(P[s], eIn[s].p0());
            spR[s] = mkLine(P[s], eOut[s].p1());
        }
        if (!clearOk) continue;

        std::vector<std::unique_ptr<gsGeometry<T>>> cand;
        cand.reserve(8);
        for (int s = 0; s < 4; ++s) {
            // Notch quad N_s: (P00=T1_s, P10=V_s, P11=T2_{s-1}, P01=P_s)
            cand.emplace_back(makeCoonsPatchFrom4Curves(
                *eOut[s].reversed().bsp, *eIn[s].reversed().bsp,
                *spL[s].bsp, *spR[s].reversed().bsp, degree));
            // Arm-tip quad A_s: (P00=P_s, P10=P_{s+1}, P11=T2_s, P01=T1_s)
            cand.emplace_back(makeCoonsPatchFrom4Curves(
                *outerTop[s].bsp, *spL[(s + 1) % 4].bsp,
                *eTip[s].bsp, *spR[s].bsp, degree));
        }

        T sc = 0, md = 0, ma = 0;
        if (evaluateRingPatches(cand, holePoly, sc, md, ma)) {
            if (!found || sc < bestScoreOut) {
                found = true;
                bestScoreOut = sc;
                bestMaxAngOut = ma;
                bestOut = std::move(cand);
            }
        }
    }
    return found;
}

// Build conforming 4-patch quad ring inside a 4-sided block [oS, oE, oN, oW] around an arbitrary hole.
template <class T>
inline void buildBlockQuadRing(
    const gsCurve<T>& oS, const gsCurve<T>& oE,
    const gsCurve<T>& oN, const gsCurve<T>& oW,
    const gsCurveLoop<T>& hole,
    gsMultiPatch<T>& mp,
    int degree)
{
    auto dist2 = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        T dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    };

    // Auto-orient the 4 outer block curves consistently:
    // cS: P0(SW)->P1(SE), cE: P1->P2(NE), cW: P0->P3(NW), cN: P3->P2
    CurveSeg<T> cS(oS, degree), cE(oE, degree), cN(oN, degree), cW(oW, degree);
    gsPoint2D<T> P0 = cS.p0(), P1 = cS.p1();
    if (dist2(cE.p0(), P1) > dist2(cE.p1(), P1)) cE = cE.reversed();
    gsPoint2D<T> P2 = cE.p1();
    if (dist2(cW.p0(), P0) > dist2(cW.p1(), P0)) cW = cW.reversed();
    gsPoint2D<T> P3 = cW.p1();
    if (dist2(cN.p0(), P3) > dist2(cN.p1(), P3)) cN = cN.reversed();
    const gsPoint2D<T> P[4] = { P0, P1, P2, P3 };

    // Orient hole curves strictly CCW and sample polygon for inside-hole checks
    std::vector<CurveSeg<T>> hCCW = orientLoopSegments(hole, degree, true);
    std::vector<gsPoint2D<T>> holePoly;
    holePoly.reserve(hCCW.size() * 12);
    for (const auto& seg : hCCW) {
        gsPoint2D<T> a = seg.p0(), b = seg.p1();
        gsPoint2D<T> pm = evalCurvePoint(*seg.bsp, T(0.5));
        gsPoint2D<T> pq = evalCurvePoint(*seg.bsp, T(0.25));
        T dMid = std::hypot(pm.x - T(0.5) * (a.x + b.x), pm.y - T(0.5) * (a.y + b.y));
        T dQtr = std::hypot(pq.x - (T(0.75) * a.x + T(0.25) * b.x), pq.y - (T(0.75) * a.y + T(0.25) * b.y));
        if (dMid < T(1e-6) && dQtr < T(1e-6)) {
            holePoly.push_back(a);
        } else {
            for (int k = 0; k < 12; ++k)
                holePoly.push_back(evalCurvePoint(*seg.bsp, static_cast<T>(k) / T(12)));
        }
    }
    const std::vector<T> hTurns = computeLoopVertexTurnDegs(hCCW);

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    // Helper: assemble and evaluate a 4-patch ring candidate from 4 sector chains inside an arbitrary 4-sided block
    auto tryBuildFrom4SectorsOnBlock = [&](const CurveSeg<T>& bcS, const CurveSeg<T>& bcE,
                                           const CurveSeg<T>& bcN, const CurveSeg<T>& bcW,
                                           const gsPoint2D<T> bP[4],
                                           const std::vector<CurveSeg<T>> sec[4],
                                           std::vector<std::unique_ptr<gsGeometry<T>>>& cand,
                                           bool& hasInteriorSharp, T& score, T& minDet, T& maxAng) -> bool {
        cand.clear();
        hasInteriorSharp = false;
        const T radToDeg = T(180) / T(3.14159265358979323846);
        for (int s = 0; s < 4; ++s) {
            if (sec[s].empty()) return false;
            for (size_t j = 1; j < sec[s].size(); ++j) {
                auto d0 = evalCurveDeriv(*sec[s][j - 1].bsp, T(1));
                auto d1 = evalCurveDeriv(*sec[s][j].bsp, T(0));
                T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
                if (n0 > T(1e-14) && n1 > T(1e-14)) {
                    T cr = (d0.x * d1.y - d0.y * d1.x) / (n0 * n1);
                    T dt = (d0.x * d1.x + d0.y * d1.y) / (n0 * n1);
                    if (std::abs(std::atan2(cr, dt)) * radToDeg > T(28)) {
                        hasInteriorSharp = true;
                        return false;
                    }
                }
            }
        }
        gsPoint2D<T> H[4];
        for (int s = 0; s < 4; ++s) H[s] = sec[s].front().p0();
        for (int s = 0; s < 4; ++s) {
            if (!isHoleSpokeClearOfHole(H[s], bP[s], holePoly))
                return false;
        }
        std::shared_ptr<gsBSpline<T>> hJoin[4];
        for (int s = 0; s < 4; ++s) hJoin[s] = joinCurveSegsExact(sec[s], degree);

        // hJoin[0]: H0->H1 (S north), hJoin[1]: H1->H2 (E west),
        // hJoin[2]: H2->H3 (needs reversed -> H3->H2 for N south),
        // hJoin[3]: H3->H0 (needs reversed -> H0->H3 for W east)
        auto hN = CurveSeg<T>(*hJoin[2]).reversed().bsp;
        auto hW = CurveSeg<T>(*hJoin[3]).reversed().bsp;

        CurveSeg<T> spOut[4], spIn[4];
        for (int s = 0; s < 4; ++s) {
            bool wc = false;
            spOut[s] = makeExteriorConeSafeSpoke(sec[(s + 3) % 4].back(), sec[s].front(), H[s], bP[s], holePoly, degree, wc);
            spIn[s]  = spOut[s].reversed();
        }

        cand.emplace_back(makeCoonsPatchFrom4Curves(*bcS.bsp, *spIn[1].bsp, *hJoin[0], *spIn[0].bsp, degree));
        cand.emplace_back(makeCoonsPatchFrom4Curves(*spOut[1].bsp, *bcE.bsp, *spOut[2].bsp, *hJoin[1], degree));
        cand.emplace_back(makeCoonsPatchFrom4Curves(*hN, *spOut[2].bsp, *bcN.bsp, *spOut[3].bsp, degree));
        cand.emplace_back(makeCoonsPatchFrom4Curves(*spIn[0].bsp, *hW, *spIn[3].bsp, *bcW.bsp, degree));
        std::string fMsg;
        bool ok = evaluateRingPatches(cand, holePoly, score, minDet, maxAng, &fMsg);
        if (!ok)
            ok = bilinearBlendRescueRing(cand, holePoly, score, minDet, maxAng, &fMsg);
        if (!ok && std::getenv("GISMO_SPLITTER_DEBUG") && hCCW.size() <= 4)
            std::cerr << "    [4SEC fail K=" << hCCW.size() << "] " << fMsg << "\n";
        return ok;
    };
    auto tryBuildFrom4Sectors = [&](const std::vector<CurveSeg<T>> sec[4],
                                    std::vector<std::unique_ptr<gsGeometry<T>>>& cand,
                                    bool& hasInteriorSharp, T& score, T& minDet, T& maxAng) -> bool {
        return tryBuildFrom4SectorsOnBlock(cS, cE, cN, cW, P, sec, cand, hasInteriorSharp, score, minDet, maxAng);
    };

    std::vector<std::unique_ptr<gsGeometry<T>>> bestSmooth4;
    T bestSmoothScore = T(1e18), bestSmoothMaxAng = T(180);
    bool foundSmooth4 = false;

    const int K = static_cast<int>(hCCW.size());
    auto considerSplits = [&](const std::vector<CurveSeg<T>>& ring) {
        const int M = static_cast<int>(ring.size());
        if (M < 4) return;
        if (M == 4) {
            for (int rot = 0; rot < 4; ++rot) {
                std::vector<CurveSeg<T>> sec[4];
                for (int s = 0; s < 4; ++s) sec[s] = { ring[(rot + s) % 4] };
                std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                bool hasSharp = false; T sc = 0, md = 0, ma = 0;
                if (tryBuildFrom4Sectors(sec, cand, hasSharp, sc, md, ma)) {
                    if (!hasSharp && (!foundSmooth4 || sc < bestSmoothScore)) {
                        foundSmooth4 = true; bestSmoothScore = sc; bestSmoothMaxAng = ma; bestSmooth4 = std::move(cand);
                    }
                }
            }
            return;
        }
        // Choose 4 anchor indices aligned with P[0..3] rays from hole centroid + neighbor shifts
        gsPoint2D<T> holeCtr(0, 0);
        for (const auto& hp : holePoly) { holeCtr.x += hp.x; holeCtr.y += hp.y; }
        holeCtr.x /= holePoly.size(); holeCtr.y /= holePoly.size();
        int baseAnchor[4];
        for (int s = 0; s < 4; ++s) {
            T angP = std::atan2(P[s].y - holeCtr.y, P[s].x - holeCtr.x);
            T bestDiff = T(1e9); int bestIdx = 0;
            for (int i = 0; i < M; ++i) {
                gsPoint2D<T> v = ring[i].p0();
                T a = std::atan2(v.y - holeCtr.y, v.x - holeCtr.x);
                T d = std::abs(a - angP);
                if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
                if (d < bestDiff) { bestDiff = d; bestIdx = i; }
            }
            baseAnchor[s] = bestIdx;
        }
        // Try small offsets around baseAnchor while preserving strictly increasing cyclic order
        for (int d0 = -1; d0 <= 1; ++d0)
        for (int d1 = -1; d1 <= 1; ++d1)
        for (int d2 = -1; d2 <= 1; ++d2)
        for (int d3 = -1; d3 <= 1; ++d3) {
            int a[4] = {
                (baseAnchor[0] + d0 + M) % M, (baseAnchor[1] + d1 + M) % M,
                (baseAnchor[2] + d2 + M) % M, (baseAnchor[3] + d3 + M) % M
            };
            // Check distinct and CCW cyclic ordering
            bool distinct = true;
            for (int x = 0; x < 4; ++x)
                for (int y = x + 1; y < 4; ++y)
                    if (a[x] == a[y]) distinct = false;
            if (!distinct) continue;
            int stepSum = 0;
            std::vector<CurveSeg<T>> sec[4];
            for (int s = 0; s < 4; ++s) {
                int len = (a[(s + 1) % 4] - a[s] + M) % M;
                if (len <= 0) { stepSum = -1; break; }
                stepSum += len;
                for (int j = 0; j < len; ++j) sec[s].push_back(ring[(a[s] + j) % M]);
            }
            if (stepSum != M) continue;
            std::vector<std::unique_ptr<gsGeometry<T>>> cand;
            bool hasSharp = false; T sc = 0, md = 0, ma = 0;
            if (tryBuildFrom4Sectors(sec, cand, hasSharp, sc, md, ma)) {
                if (!hasSharp && (!foundSmooth4 || sc < bestSmoothScore)) {
                    foundSmooth4 = true; bestSmoothScore = sc; bestSmoothMaxAng = ma; bestSmooth4 = std::move(cand);
                }
            }
        }
    };

    if (K == 1) {
        // Split single closed curve into 4 equal-arc segments starting near angle toward P0
        gsPoint2D<T> holeCtr(0, 0);
        for (const auto& hp : holePoly) { holeCtr.x += hp.x; holeCtr.y += hp.y; }
        holeCtr.x /= holePoly.size(); holeCtr.y /= holePoly.size();
        T targetAng = std::atan2(P0.y - holeCtr.y, P0.x - holeCtr.x);
        T bestU = T(0), bestDiff = T(1e9);
        for (int k = 0; k < 64; ++k) {
            T u = static_cast<T>(k) / T(64);
            gsPoint2D<T> q = evalCurvePoint(*hCCW[0].bsp, u);
            T d = std::abs(std::atan2(q.y - holeCtr.y, q.x - holeCtr.x) - targetAng);
            if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
            if (d < bestDiff) { bestDiff = d; bestU = u; }
        }
        std::vector<CurveSeg<T>> r4;
        CurveSeg<T> shifted = hCCW[0];
        if (bestU > T(1e-4) && bestU < T(1) - T(1e-4)) {
            auto sp = hCCW[0].splitAtParam(bestU);
            auto joined = joinCurveSegsExact(std::vector<CurveSeg<T>>{sp.second, sp.first}, degree);
            shifted = CurveSeg<T>(*joined);
        }
        CurveSeg<T> rem = shifted;
        for (int k = 0; k < 3; ++k) {
            auto sp = rem.splitArc(T(1) / static_cast<T>(4 - k));
            r4.push_back(sp.first);
            rem = sp.second;
        }
        r4.push_back(rem);
        considerSplits(r4);
    } else if (K == 2) {
        std::vector<CurveSeg<T>> r4;
        for (int i = 0; i < 2; ++i) {
            auto sp = hCCW[i].splitArc(T(0.5));
            r4.push_back(sp.first);
            r4.push_back(sp.second);
        }
        considerSplits(r4);
    } else if (K == 3) {
        // Split each single edge choice in turn (at t=0.5) -> 4-segment ring keeping all 3 sharp corners as vertices
        for (int eSplit = 0; eSplit < 3; ++eSplit) {
            std::vector<CurveSeg<T>> r4;
            for (int i = 0; i < 3; ++i) {
                if (i == eSplit) {
                    auto sp = hCCW[i].splitArc(T(0.5));
                    r4.push_back(sp.first);
                    r4.push_back(sp.second);
                } else {
                    r4.push_back(hCCW[i]);
                }
            }
            considerSplits(r4);
        }
    } else if (K == 4) {
        considerSplits(hCCW);
        // Also test midpoint-shifted 4-ring if hole corners are smooth (<= 25 deg turns)
        bool allSmooth = true;
        for (T t : hTurns) if (t > T(25)) { allSmooth = false; break; }
        if (allSmooth) {
            std::vector<CurveSeg<T>> halves;
            for (int i = 0; i < 4; ++i) {
                auto sp = hCCW[i].splitArc(T(0.5));
                halves.push_back(sp.first);
                halves.push_back(sp.second);
            }
            std::vector<CurveSeg<T>> rot4;
            for (int i = 0; i < 4; ++i) {
                auto j = joinCurveSegsExact(std::vector<CurveSeg<T>>{halves[(2 * i + 1) % 8], halves[(2 * i + 2) % 8]}, degree);
                rot4.emplace_back(*j);
            }
            considerSplits(rot4);
            // Also test radial angle-aligned 4-way split toward block corners P[0..3]
            auto joinedFull = joinCurveSegsExact(hCCW, degree);
            gsPoint2D<T> hc(0, 0);
            for (const auto& hp : holePoly) { hc.x += hp.x; hc.y += hp.y; }
            hc.x /= holePoly.size(); hc.y /= holePoly.size();
            const int NS = 128;
            std::vector<T> sampleU(NS), sampleAng(NS);
            for (int k = 0; k < NS; ++k) {
                sampleU[k] = static_cast<T>(k) / static_cast<T>(NS);
                gsPoint2D<T> q = evalCurvePoint(*joinedFull, sampleU[k]);
                sampleAng[k] = std::atan2(q.y - hc.y, q.x - hc.x);
            }
            T cutU[4];
            for (int s = 0; s < 4; ++s) {
                T tAng = std::atan2(P[s].y - hc.y, P[s].x - hc.x);
                T bestD = T(1e9); int bi = 0;
                for (int k = 0; k < NS; ++k) {
                    T d = std::abs(sampleAng[k] - tAng);
                    if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
                    if (d < bestD) { bestD = d; bi = k; }
                }
                cutU[s] = sampleU[bi];
            }
            // Shift parameter origin to cutU[0] and sort remaining cuts
            CurveSeg<T> fullSeg(*joinedFull);
            if (cutU[0] > T(1e-4) && cutU[0] < T(1) - T(1e-4)) {
                auto sp0 = fullSeg.splitAtParam(cutU[0]);
                fullSeg = CurveSeg<T>(*joinCurveSegsExact(std::vector<CurveSeg<T>>{sp0.second, sp0.first}, degree));
            }
            // Locate P[1..3] along shifted fullSeg and split at exact radial params
            T uShift[3];
            bool uOk = true;
            for (int s = 1; s < 4; ++s) {
                T tAng = std::atan2(P[s].y - hc.y, P[s].x - hc.x);
                T bestD = T(1e9), bu = T(0);
                for (int k = 1; k < NS; ++k) {
                    T u = static_cast<T>(k) / static_cast<T>(NS);
                    gsPoint2D<T> q = evalCurvePoint(*fullSeg.bsp, u);
                    T d = std::abs(std::atan2(q.y - hc.y, q.x - hc.x) - tAng);
                    if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
                    if (d < bestD) { bestD = d; bu = u; }
                }
                uShift[s - 1] = bu;
            }
            if (uShift[0] > T(0.05) && uShift[0] < uShift[1] - T(0.05) &&
                uShift[1] < uShift[2] - T(0.05) && uShift[2] < T(0.95))
            {
                auto sp0 = fullSeg.splitAtParam(uShift[0]);
                T u1l = (uShift[1] - uShift[0]) / (T(1) - uShift[0]);
                auto sp1 = sp0.second.splitAtParam(u1l);
                T u2l = (uShift[2] - uShift[1]) / (T(1) - uShift[1]);
                auto sp2 = sp1.second.splitAtParam(u2l);
                std::vector<CurveSeg<T>> rRadial = {sp0.first, sp1.first, sp2.first, sp2.second};
                {
                    std::vector<CurveSeg<T>> sec[4];
                    for (int s = 0; s < 4; ++s) sec[s] = { rRadial[s] };
                    std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                    bool hasSharp = false; T sc = 0, md = 0, ma = 0;
                    if (tryBuildFrom4Sectors(sec, cand, hasSharp, sc, md, ma)) {
                        if (!hasSharp && (!foundSmooth4 || sc < bestSmoothScore)) {
                            foundSmooth4 = true; bestSmoothScore = sc; bestSmoothMaxAng = ma; bestSmooth4 = std::move(cand);
                        }
                    }
                }
            }
            (void)uOk;
        }
    } else {
        // K >= 5: first test 4-sector ring (succeeds as smooth 4-patch when turns are small/smooth)
        considerSplits(hCCW);
    }

    if (foundSmooth4 && bestSmoothMaxAng <= T(163.5)) {
        if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSmooth4\n";
        for (auto& p : bestSmooth4)
            mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
        return;
    }

    gsPoint2D<T> holeCtr(0, 0);
    for (const auto& hp : holePoly) { holeCtr.x += hp.x; holeCtr.y += hp.y; }
    if (!holePoly.empty()) { holeCtr.x /= holePoly.size(); holeCtr.y /= holePoly.size(); }
    T holeR = T(0);
    for (const auto& hp : holePoly)
        holeR = std::max(holeR, std::hypot(hp.x - holeCtr.x, hp.y - holeCtr.y));
    T maxBlockDist = T(0);
    for (int s = 0; s < 4; ++s)
        maxBlockDist = std::max(maxBlockDist, std::hypot(P[s].x - holeCtr.x, P[s].y - holeCtr.y));

    auto runSharpSearchOnBlock = [&](const CurveSeg<T>& bcS, const CurveSeg<T>& bcE,
                                     const CurveSeg<T>& bcN, const CurveSeg<T>& bcW,
                                     const gsPoint2D<T> bP[4],
                                     const std::vector<CurveSeg<T>>& workH,
                                     const std::vector<T>& splitTs,
                                     const std::vector<T>& oWeights,
                                     bool earlyStop,
                                     std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
                                     T& bestScoreOut, T& bestMaxAngOut,
                                     int* dbgVisOut = nullptr, int* dbgCapOut = nullptr, int* dbgEvalOut = nullptr,
                                     std::string* firstCapOut = nullptr, std::string* firstFailOut = nullptr) -> bool
    {
        return runSharpBlockRingSearch(
            bcS, bcE, bcN, bcW, bP, workH, holePoly, degree, splitTs, oWeights, earlyStop,
            bestOut, bestScoreOut, bestMaxAngOut, dbgVisOut, dbgCapOut, dbgEvalOut, firstCapOut, firstFailOut);
    };

    auto tryBuildBufferedRing = [&](std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
                                    T& bestMaxAngOut,
                                    bool allowSharpFallback = true) -> bool
    {
        if (holeR <= T(1e-12)) return false;
        const T pi = T(3.14159265358979323846);
        std::vector<gsPoint2D<T>> blockPoly;
        blockPoly.reserve(32);
        const CurveSeg<T> outerSegs[4] = { cS, cE, cN.reversed(), cW.reversed() };
        for (int s = 0; s < 4; ++s) {
            for (int k = 0; k < 8; ++k)
                blockPoly.push_back(evalCurvePoint(*outerSegs[s].bsp, static_cast<T>(k) / T(8)));
        }
        if (!pointInPolyStrict(holeCtr, blockPoly, T(1e-7))) return false;

        T angP[4];
        for (int s = 0; s < 4; ++s)
            angP[s] = std::atan2(P[s].y - holeCtr.y, P[s].x - holeCtr.x);

        T angSq[4];
        bool haveNotchSq = false;
        if (K == 12) {
            std::vector<int> notches;
            for (int i = 0; i < 12; ++i) {
                const auto& prev = hCCW[(i + 11) % 12];
                const auto& next = hCCW[i];
                gsPoint2D<T> d(next.p1().x - prev.p0().x, next.p1().y - prev.p0().y);
                T aStart = 0, dCCW = 0, aDir = 0, minEdge = 0;
                if (computeExteriorCone(prev, next, d, aStart, dCCW, aDir, minEdge) &&
                    dCCW < (T(150) * pi / T(180)))
                    notches.push_back(i);
            }
            if (notches.size() == 4) {
                T nAng[4];
                for (int s = 0; s < 4; ++s) {
                    gsPoint2D<T> v = hCCW[notches[s]].p0();
                    nAng[s] = std::atan2(v.y - holeCtr.y, v.x - holeCtr.x);
                }
                int bestRot = 0;
                T bestAlign = -T(1e9);
                for (int r = 0; r < 4; ++r) {
                    T al = T(0);
                    for (int s = 0; s < 4; ++s) al += std::cos(nAng[(s + r) % 4] - angP[s]);
                    if (al > bestAlign) { bestAlign = al; bestRot = r; }
                }
                for (int s = 0; s < 4; ++s) angSq[s] = nAng[(s + bestRot) % 4];
                haveNotchSq = true;
            }
        }
        if (!haveNotchSq) {
            T sumC = T(0), sumS = T(0);
            for (int s = 0; s < 4; ++s) {
                T a = angP[s] - static_cast<T>(s) * (pi / T(2));
                sumC += std::cos(a);
                sumS += std::sin(a);
            }
            T baseA = std::atan2(sumS, sumC);
            for (int s = 0; s < 4; ++s)
                angSq[s] = baseA + static_cast<T>(s) * (pi / T(2));
        }

        auto wrapPi = [&](T a) -> T {
            while (a >  pi) a -= T(2) * pi;
            while (a < -pi) a += T(2) * pi;
            return a;
        };

        bool foundBuf = false;
        size_t bestCount = std::numeric_limits<size_t>::max();
        T bestScore = T(1e18);
        bestMaxAngOut = T(180);

        auto evalBox = [&](const gsPoint2D<T> B[4]) -> bool {
            for (int s = 0; s < 4; ++s) {
                if (pointInPolyStrict(B[s], holePoly, T(1e-6)) ||
                    !pointInPolyStrict(B[s], blockPoly, T(1e-6)) ||
                    !isHoleSpokeClearOfHole(B[s], P[s], holePoly))
                    return false;
            }
            for (int s = 0; s < 4; ++s) {
                if (cross2(B[s], B[(s + 1) % 4], B[(s + 2) % 4]) <= T(1e-10))
                    return false;
                gsPoint2D<T> mid(T(0.5) * (B[s].x + B[(s + 1) % 4].x),
                                 T(0.5) * (B[s].y + B[(s + 1) % 4].y));
                if (!pointInPolyStrict(mid, blockPoly, T(1e-6)) ||
                    !isHoleSpokeClearOfHole(B[s], B[(s + 1) % 4], holePoly))
                    return false;
            }

            CurveSeg<T> bS = mkLine(B[0], B[1]);
            CurveSeg<T> bE = mkLine(B[1], B[2]);
            CurveSeg<T> bN = mkLine(B[3], B[2]);
            CurveSeg<T> bW = mkLine(B[0], B[3]);

            std::vector<std::unique_ptr<gsGeometry<T>>> outer4;
            outer4.reserve(4);
            CurveSeg<T> spOut[4], spIn[4];
            for (int s = 0; s < 4; ++s) {
                spOut[s] = mkLine(B[s], P[s]);
                spIn[s]  = spOut[s].reversed();
            }
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*cS.bsp, *spIn[1].bsp, *bS.bsp, *spIn[0].bsp, degree));
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*spOut[1].bsp, *cE.bsp, *spOut[2].bsp, *bE.bsp, degree));
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*bN.bsp, *spOut[2].bsp, *cN.bsp, *spOut[3].bsp, degree));
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*spIn[0].bsp, *bW.bsp, *spIn[3].bsp, *cW.bsp, degree));

            T outSc = 0, outMd = 0, outMa = 0;
            if (!evaluateRingPatches(outer4, holePoly, outSc, outMd, outMa))
                return false;

            std::vector<std::unique_ptr<gsGeometry<T>>> innerCand;
            T inSc = 0, inMa = T(180);
            bool inOk = false;

            if (K == 12) {
                inOk = tryBuildCross8QuadRing(bS, bE, bN, bW, hCCW, holePoly, innerCand, inSc, inMa, degree);
            } else if (K == 3 || K == 4) {
                auto tryRing4 = [&](const std::vector<CurveSeg<T>>& r4) {
                    for (int rot = 0; rot < 4; ++rot) {
                        std::vector<CurveSeg<T>> sec[4];
                        for (int s = 0; s < 4; ++s) sec[s] = { r4[(rot + s) % 4] };
                        std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                        bool hs = false; T sc = 0, md = 0, ma = 0;
                        if (tryBuildFrom4SectorsOnBlock(bS, bE, bN, bW, B, sec, cand, hs, sc, md, ma) && !hs) {
                            if (!inOk || sc < inSc) {
                                inOk = true; inSc = sc; inMa = ma; innerCand = std::move(cand);
                            }
                        }
                    }
                };
                if (K == 4) {
                    tryRing4(hCCW);
                } else {
                    for (int eSplit = 0; eSplit < 3; ++eSplit) {
                        std::vector<CurveSeg<T>> r4;
                        for (int i = 0; i < 3; ++i) {
                            if (i == eSplit) {
                                auto sp = hCCW[i].splitArc(T(0.5));
                                r4.push_back(sp.first);
                                r4.push_back(sp.second);
                            } else {
                                r4.push_back(hCCW[i]);
                            }
                        }
                        tryRing4(r4);
                    }
                }
            }
            if (!inOk && allowSharpFallback && K >= 3) {
                std::vector<CurveSeg<T>> wH;
                if (K == 3 || K == 4) {
                    for (const auto& s : hCCW) {
                        auto sp = s.splitArc(T(0.5));
                        wH.push_back(sp.first);
                        wH.push_back(sp.second);
                    }
                } else {
                    wH = hCCW;
                    if (wH.size() % 2 != 0) {
                        int longestIdx = 0;
                        T bestL = -T(1);
                        for (size_t i = 0; i < wH.size(); ++i) {
                            T d = std::hypot(wH[i].p1().x - wH[i].p0().x, wH[i].p1().y - wH[i].p0().y);
                            if (d > bestL) { bestL = d; longestIdx = static_cast<int>(i); }
                        }
                        auto sp = wH[longestIdx].splitArc(T(0.5));
                        wH[longestIdx] = sp.first;
                        wH.insert(wH.begin() + longestIdx + 1, sp.second);
                    }
                }
                inOk = runSharpSearchOnBlock(
                    bS, bE, bN, bW, B, wH,
                    {T(0.5)}, {T(0.62), T(0.45)}, true,
                    innerCand, inSc, inMa);
            }

            if (inOk) {
                size_t totN = outer4.size() + innerCand.size();
                T totMa = std::max(outMa, inMa);
                T totSc = outSc + inSc;
                if (!foundBuf || totN < bestCount || (totN == bestCount && totSc < bestScore)) {
                    foundBuf = true;
                    bestCount = totN;
                    bestScore = totSc;
                    bestMaxAngOut = totMa;
                    bestOut.clear();
                    bestOut.reserve(totN);
                    for (auto& p : outer4) bestOut.push_back(std::move(p));
                    for (auto& p : innerCand) bestOut.push_back(std::move(p));
                }
                if (totMa <= T(158) && ((K == 12 && totN <= 12) || (K <= 4 && totN <= 8) || (K >= 5 && totN <= 20)))
                    return true;
            }
            return false;
        };

        const std::vector<T> scaleList = allowSharpFallback
            ? std::vector<T>{T(1.90), T(2.25)}
            : std::vector<T>{T(1.90), T(2.25), T(1.60), T(2.65)};
        const std::vector<T> blendList = allowSharpFallback
            ? std::vector<T>{T(0.0), T(0.50)}
            : std::vector<T>{T(0.0), T(0.35), T(0.70), T(1.0)};

        for (T scale : scaleList) {
            for (T blend : blendList) {
                gsPoint2D<T> B[4];
                bool raysOk = true;
                for (int s = 0; s < 4; ++s) {
                    T th = angSq[s] + blend * wrapPi(angP[s] - angSq[s]);
                    gsPoint2D<T> dir(std::cos(th), std::sin(th));
                    T dWall = rayDistToPoly(holeCtr, dir, blockPoly);
                    if (dWall > T(1e9)) { raysOk = false; break; }
                    T dist = std::min(scale * holeR, T(0.76) * dWall);
                    B[s] = gsPoint2D<T>(holeCtr.x + dist * dir.x, holeCtr.y + dist * dir.y);
                }
                if (raysOk && evalBox(B))
                    return true;
            }
        }
        // Also test block-proportioned lerp boxes for tall/skewed cells where 1:1 boxes fail outer4 aspect checks
        const std::vector<T> lerpList = allowSharpFallback
            ? std::vector<T>{T(0.48), T(0.58)}
            : std::vector<T>{T(0.48), T(0.58), T(0.38)};
        for (T al : lerpList) {
            gsPoint2D<T> B[4];
            for (int s = 0; s < 4; ++s)
                B[s] = gsPoint2D<T>((T(1) - al) * holeCtr.x + al * P[s].x,
                                    (T(1) - al) * holeCtr.y + al * P[s].y);
            if (evalBox(B))
                return true;
        }
        return foundBuf;
    };

    // 1. Direct 8-quad cross ring for K=12 crosses when enclosing block has clean visibility & angles
    if (K == 12) {
        std::vector<std::unique_ptr<gsGeometry<T>>> cross8;
        T cSc = 0, cMa = T(180);
        if (tryBuildCross8QuadRing(cS, cE, cN, cW, hCCW, holePoly, cross8, cSc, cMa, degree) &&
            cMa <= T(163.5))
        {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked direct Cross8Quad P=8\n";
            for (auto& p : cross8)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
    }

    // 2. For K=3,4 when direct 4-sector ring failed or exceeded 163.5 deg, try compact 8-patch buffered ring before edge doubling
    if (K == 3 || K == 4) {
        std::vector<std::unique_ptr<gsGeometry<T>>> buf8;
        T bMa = T(180);
        if (tryBuildBufferedRing(buf8, bMa, false) && buf8.size() <= 8 && bMa <= T(163.5)) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked buffered K=" << K << " P=" << buf8.size() << "\n";
            for (auto& p : buf8)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
    }

    // 3. For K=12 crosses on tall/skewed blocks, try 12-patch buffered cross ring before general sharp search
    if (K == 12) {
        std::vector<std::unique_ptr<gsGeometry<T>>> buf12;
        T bMa = T(180);
        if (tryBuildBufferedRing(buf12, bMa, false) && buf12.size() <= 12 && bMa <= T(163.5)) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked buffered Cross P=" << buf12.size() << "\n";
            for (auto& p : buf12)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
    }

    // 4. For K >= 5 with sharp corners (or K=3,4 fallback), build sharp/midpoint-aligned structured ring
    if (K >= 3) {
        std::vector<CurveSeg<T>> workH;
        if (K == 3 || K == 4) {
            for (const auto& s : hCCW) {
                auto sp = s.splitArc(T(0.5));
                workH.push_back(sp.first);
                workH.push_back(sp.second);
            }
        } else {
            workH = hCCW;
            if (workH.size() % 2 != 0) {
                int longestIdx = 0;
                T bestL = -T(1);
                for (size_t i = 0; i < workH.size(); ++i) {
                    T d = std::hypot(workH[i].p1().x - workH[i].p0().x, workH[i].p1().y - workH[i].p0().y);
                    if (d > bestL) { bestL = d; longestIdx = static_cast<int>(i); }
                }
                auto sp = workH[longestIdx].splitArc(T(0.5));
                workH[longestIdx] = sp.first;
                workH.insert(workH.begin() + longestIdx + 1, sp.second);
            }
        }

        std::vector<std::unique_ptr<gsGeometry<T>>> bestSharp;
        T bestSharpScore = T(1e18), bestSharpMaxAng = T(180);
        int dbgVisCount = 0, dbgCapFail = 0, dbgEvalFail = 0;
        std::string firstCapFail, firstFail;
        bool foundSharp = runSharpSearchOnBlock(
            cS, cE, cN, cW, P, workH, {T(0.5)},
            {T(0.62), T(0.45), T(0.78)}, true,
            bestSharp, bestSharpScore, bestSharpMaxAng,
            &dbgVisCount, &dbgCapFail, &dbgEvalFail, &firstCapFail, &firstFail);
        if (!foundSharp && holeR > T(1e-12) && maxBlockDist > T(2.2) * holeR) {
            foundSharp = runSharpSearchOnBlock(
                cS, cE, cN, cW, P, workH, {T(0.25)},
                {T(0.62), T(0.45), T(0.78)}, true,
                bestSharp, bestSharpScore, bestSharpMaxAng,
                &dbgVisCount, &dbgCapFail, &dbgEvalFail, &firstCapFail, &firstFail);
        }

        if (foundSharp && (!foundSmooth4 || bestSharpMaxAng <= T(163.5) || bestSharpMaxAng <= bestSmoothMaxAng)) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSharp P=" << bestSharp.size() << "\n";
            for (auto& p : bestSharp)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
        {
            std::vector<std::unique_ptr<gsGeometry<T>>> bufBest;
            T bufMaxAng = T(180);
            if (tryBuildBufferedRing(bufBest, bufMaxAng, true) &&
                (!foundSmooth4 || bufMaxAng <= T(163.5) || bufMaxAng <= bestSmoothMaxAng))
            {
                if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked buffered fallback P=" << bufBest.size() << "\n";
                for (auto& p : bufBest)
                    mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
                return;
            }
        }
        if (foundSmooth4) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSmooth4 (backup)\n";
            for (auto& p : bestSmooth4)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
        if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] foundSharp=false for K=" << K << " (visQuad=" << dbgVisCount << ", capFail=" << dbgCapFail << "[" << firstCapFail << "], evalFail=" << dbgEvalFail << "[" << firstFail << "])\n";
    }

    if (foundSmooth4) {
        if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSmooth4 (backup)\n";
        for (auto& p : bestSmooth4)
            mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
        return;
    }

    // 5. Collar rescue: place a conforming 2-spoke-per-spike hole smoothing collar around sharp/twisted holes
    {
        std::vector<gsPoint2D<T>> blockPolyCollar;
        blockPolyCollar.reserve(32);
        const CurveSeg<T> outerSegsCol[4] = { cS, cE, cN.reversed(), cW.reversed() };
        for (int s = 0; s < 4; ++s)
            for (int k = 0; k < 8; ++k)
                blockPolyCollar.push_back(evalCurvePoint(*outerSegsCol[s].bsp, static_cast<T>(k) / T(8)));

        std::vector<std::unique_ptr<gsGeometry<T>>> collarP;
        std::vector<CurveSeg<T>> colOuter;
        std::vector<int> edgeToPatch;
        if (buildHoleSmoothingCollar(hole, &blockPolyCollar, degree, T(4), T(163.5), collarP, colOuter, edgeToPatch)) {
            if ((colOuter.size() % 2) == 1)
                splitOneCollarFlankForEvenness(collarP, colOuter, edgeToPatch, degree);
            if ((colOuter.size() % 2) == 0 && colOuter.size() >= 4) {
                std::vector<gsPoint2D<T>> colPoly;
                colPoly.reserve(colOuter.size() * 8);
                for (const auto& cs : colOuter)
                    for (int k = 0; k < 8; ++k)
                        colPoly.push_back(evalCurvePoint(*cs.bsp, static_cast<T>(k) / T(8)));

                std::vector<std::unique_ptr<gsGeometry<T>>> outerRing;
                T rSc = 0, rMa = T(180);
                bool ringOk = false;
                if (colOuter.size() == 4) {
                    for (int rot = 0; rot < 4; ++rot) {
                        CurveSeg<T> spOut[4], spIn[4];
                        std::shared_ptr<gsBSpline<T>> hJ[4];
                        bool visOk = true;
                        for (int s = 0; s < 4; ++s) {
                            const CurveSeg<T>& seg = colOuter[(rot + s) % 4];
                            const CurveSeg<T>& prevSeg = colOuter[(rot + s + 3) % 4];
                            hJ[s] = joinCurveSegsExact(std::vector<CurveSeg<T>>{seg}, degree);
                            if (!isHoleSpokeClearOfHole(seg.p0(), P[s], colPoly)) { visOk = false; break; }
                            bool wc = false;
                            spOut[s] = makeExteriorConeSafeSpoke(prevSeg, seg, seg.p0(), P[s], colPoly, degree, wc);
                            spIn[s]  = spOut[s].reversed();
                        }
                        if (!visOk) continue;
                        auto hN = CurveSeg<T>(*hJ[2]).reversed().bsp;
                        auto hW = CurveSeg<T>(*hJ[3]).reversed().bsp;
                        std::vector<std::unique_ptr<gsGeometry<T>>> cand4;
                        cand4.emplace_back(makeCoonsPatchFrom4Curves(*cS.bsp, *spIn[1].bsp, *hJ[0], *spIn[0].bsp, degree));
                        cand4.emplace_back(makeCoonsPatchFrom4Curves(*spOut[1].bsp, *cE.bsp, *spOut[2].bsp, *hJ[1], degree));
                        cand4.emplace_back(makeCoonsPatchFrom4Curves(*hN, *spOut[2].bsp, *cN.bsp, *spOut[3].bsp, degree));
                        cand4.emplace_back(makeCoonsPatchFrom4Curves(*spIn[0].bsp, *hW, *spIn[3].bsp, *cW.bsp, degree));
                        T sc = 0, md = 0, ma = 0;
                        if (evaluateRingPatches(cand4, colPoly, sc, md, ma) ||
                            bilinearBlendRescueRing(cand4, colPoly, sc, md, ma)) {
                            if (!ringOk || sc < rSc) {
                                ringOk = true; rSc = sc; rMa = ma; outerRing = std::move(cand4);
                            }
                        }
                    }
                }
                if (!ringOk) {
                    ringOk = runSharpBlockRingSearch(
                        cS, cE, cN, cW, P, colOuter, colPoly, degree,
                        {T(0.5)}, {T(0.62), T(0.45), T(0.78)}, true,
                        outerRing, rSc, rMa);
                }
                if (ringOk && rMa <= T(163.5)) {
                    for (auto& p : collarP) mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
                    for (auto& p : outerRing) mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
                    return;
                }
            }
        }
    }

    // All valid ring constructions failed: leave mp unchanged so callers fail fast without untangling folded patches.
    if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] FAILED for K=" << K << "\n";
    (void)holeCtr;
}

// Conforming quad ring between a 4-curve outer loop and an arbitrary inner hole
template <class T>
inline void buildQuadRingByProximity(
    const gsCurveLoop<T>& outer,
    const gsCurveLoop<T>& hole,
    gsMultiPatch<T>& mp,
    int degree)
{
    if (outer.size() == 4 && hole.size() == 4) {
        // Check if 1-to-1 curve proximity gives valid single-sign Coons ring (ideal for curved 4-arc annulus/plate)
        gsMultiPatch<T> fastMp;
        for (int i = 0; i < 4; ++i) {
            const auto& oCrv = outer.curve(i);
            auto pA = startPoint(oCrv);
            auto pB = endPoint(oCrv);
            int best_h = -1;
            bool best_rev = false;
            T best_dist = T(1e18);
            for (int h = 0; h < 4; ++h) {
                const auto& hCrv = hole.curve(h);
                auto hA = startPoint(hCrv);
                auto hB = endPoint(hCrv);
                T d_fwd = (pA.x - hA.x)*(pA.x - hA.x) + (pA.y - hA.y)*(pA.y - hA.y) +
                          (pB.x - hB.x)*(pB.x - hB.x) + (pB.y - hB.y)*(pB.y - hB.y);
                if (d_fwd < best_dist) { best_dist = d_fwd; best_h = h; best_rev = false; }
                T d_rev = (pA.x - hB.x)*(pA.x - hB.x) + (pA.y - hB.y)*(pA.y - hB.y) +
                          (pB.x - hA.x)*(pB.x - hA.x) + (pB.y - hA.y)*(pB.y - hA.y);
                if (d_rev < best_dist) { best_dist = d_rev; best_h = h; best_rev = true; }
            }
            std::unique_ptr<gsCurve<T>> hCrv(hole.curve(best_h).clone());
            if (best_rev) safeReverseCurve(*hCrv);
            auto hA = startPoint(*hCrv);
            auto hB = endPoint(*hCrv);
            std::unique_ptr<gsBSpline<T>> sp_S(makeStraightBSpline(hA, pA, degree));
            std::unique_ptr<gsBSpline<T>> sp_N(makeStraightBSpline(hB, pB, degree));
            auto patch = makeCoonsPatchFrom4Curves(*sp_S, oCrv, *sp_N, *hCrv, degree);
            fastMp.addPatch(typename gsGeometry<T>::uPtr(patch));
        }
        bool ok = (fastMp.nPatches() == 4);
        if (ok) {
            fastMp.computeTopology(T(1e-5), false, false);
            if (fastMp.topology().nInterfaces() != 4) ok = false;
        }
        for (size_t p = 0; ok && p < fastMp.nPatches(); ++p) {
            if (!isPatchDetSingleSign(fastMp.patch(p), 15, T(1e-8))) ok = false;
        }
        if (ok) {
            for (size_t p = 0; p < fastMp.nPatches(); ++p)
                mp.addPatch(typename gsGeometry<T>::uPtr(fastMp.patch(p).clone().release()));
            return;
        }
    }

    if (outer.size() == 4) {
        std::vector<CurveSeg<T>> oCCW = orientLoopSegments(outer, degree, true);
        buildBlockQuadRing(*oCCW[0].bsp, *oCCW[1].bsp, *oCCW[2].reversed().bsp, *oCCW[3].reversed().bsp, hole, mp, degree);
    }
}

// Create an exact or cubic-interpolated CurveSeg connecting points at iso-parameter on patch.
template <class T>
inline CurveSeg<T> extractBlockCutCurve(
    const gsGeometry<T>& refPatch,
    T u0, T v0, T u1, T v1, int degree)
{
    gsMatrix<T> uv(2, 4), pts;
    for (int k = 0; k < 4; ++k) {
        T t = static_cast<T>(k) / T(3);
        uv(0, k) = (T(1) - t) * u0 + t * u1;
        uv(1, k) = (T(1) - t) * v0 + t * v1;
    }
    refPatch.eval_into(uv, pts);
    gsPoint2D<T> q0(pts(0, 0), pts(1, 0)), q1(pts(0, 1), pts(1, 1));
    gsPoint2D<T> q2(pts(0, 2), pts(1, 2)), q3(pts(0, 3), pts(1, 3));

    // Check if essentially straight
    auto distToChord = [&](const gsPoint2D<T>& q) -> T {
        T dx = q3.x - q0.x, dy = q3.y - q0.y;
        T L = std::hypot(dx, dy);
        if (L <= T(1e-14)) return T(0);
        return std::abs(dx * (q.y - q0.y) - dy * (q.x - q0.x)) / L;
    };
    if (distToChord(q1) < T(1e-5) && distToChord(q2) < T(1e-5)) {
        std::unique_ptr<gsBSpline<T>> bsp(makeStraightBSpline(q0, q3, std::max(3, degree)));
        return CurveSeg<T>(*bsp);
    }

    // Closed-form cubic Bezier interpolating q0, q1(1/3), q2(2/3), q3(1)
    gsPoint2D<T> p1((-T(5) * q0.x + T(18) * q1.x - T(9) * q2.x + T(2) * q3.x) / T(6),
                    (-T(5) * q0.y + T(18) * q1.y - T(9) * q2.y + T(2) * q3.y) / T(6));
    gsPoint2D<T> p2((T(2) * q0.x - T(9) * q1.x + T(18) * q2.x - T(5) * q3.x) / T(6),
                    (T(2) * q0.y - T(9) * q1.y + T(18) * q2.y - T(5) * q3.y) / T(6));
    gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
    gsMatrix<T> coefs(4, 2);
    coefs << q0.x, q0.y, p1.x, p1.y, p2.x, p2.y, q3.x, q3.y;
    gsBSpline<T> bsp(kv, coefs);
    if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
    return CurveSeg<T>(bsp);
}

// Decompose an arbitrary 4-sided macro-block [cS(P0->P1), cE(P1->P2), cN(P3->P2), cW(P0->P3)]
// enclosing 0, 1, or N holes into a conforming multipatch without T-junctions.
template <class T>
inline bool build4SidedMultiHoleGrid(
    const CurveSeg<T>& cS, const CurveSeg<T>& cE,
    const CurveSeg<T>& cN, const CurveSeg<T>& cW,
    const std::vector<const gsCurveLoop<T>*>& holes,
    gsMultiPatch<T>& mp,
    int degree)
{
    auto dist2 = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        T dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    };
    CurveSeg<T> oS = cS, oE = cE, oN = cN, oW = cW;
    gsPoint2D<T> P0 = oS.p0(), P1 = oS.p1();
    if (dist2(oE.p0(), P1) > dist2(oE.p1(), P1)) oE = oE.reversed();
    if (dist2(oW.p0(), P0) > dist2(oW.p1(), P0)) oW = oW.reversed();
    gsPoint2D<T> P3 = oW.p1();
    if (dist2(oN.p0(), P3) > dist2(oN.p1(), P3)) oN = oN.reversed();

    if (holes.empty()) {
        mp.addPatch(typename gsGeometry<T>::uPtr(
            makeCoonsPatchFrom4Curves(*oS.bsp, *oE.bsp, *oN.bsp, *oW.bsp, degree)));
        return true;
    }
    if (holes.size() == 1) {
        buildBlockQuadRing(*oS.bsp, *oE.bsp, *oN.bsp, *oW.bsp, *holes[0], mp, degree);
        return true;
    }

    // Build smooth 4-sided reference Coons map over [0,1]^2
    std::unique_ptr<gsTensorBSpline<2, T>> refPatch(
        makeCoonsPatchFrom4Curves(*oS.bsp, *oE.bsp, *oN.bsp, *oW.bsp, degree));

    // Sample refPatch on 41x41 grid to invert physical points to (u, v) in [0,1]^2
    const int G = 41;
    std::vector<gsPoint2D<T>> gPts(G * G);
    {
        gsMatrix<T> uv(2, G * G), val;
        for (int j = 0; j < G; ++j)
            for (int i = 0; i < G; ++i) {
                uv(0, j * G + i) = static_cast<T>(i) / (G - 1);
                uv(1, j * G + i) = static_cast<T>(j) / (G - 1);
            }
        refPatch->eval_into(uv, val);
        for (int idx = 0; idx < G * G; ++idx)
            gPts[idx] = gsPoint2D<T>(val(0, idx), val(1, idx));
    }

    auto invertPointUV = [&](const gsPoint2D<T>& q) -> gsPoint2D<T> {
        T bestD = T(1e18); int bi = 0, bj = 0;
        for (int j = 0; j < G; ++j)
            for (int i = 0; i < G; ++i) {
                T d = dist2(q, gPts[j * G + i]);
                if (d < bestD) { bestD = d; bi = i; bj = j; }
            }
        T u = static_cast<T>(bi) / (G - 1), v = static_cast<T>(bj) / (G - 1);
        gsMatrix<T> uv(2, 1), pt, jac;
        for (int it = 0; it < 6; ++it) {
            uv(0, 0) = u; uv(1, 0) = v;
            refPatch->eval_into(uv, pt);
            refPatch->jacobian_into(uv, jac);
            T rx = q.x - pt(0, 0), ry = q.y - pt(1, 0);
            T det = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
            if (std::abs(det) < T(1e-14)) break;
            T du = ( jac(1, 1) * rx - jac(0, 1) * ry) / det;
            T dv = (-jac(1, 0) * rx + jac(0, 0) * ry) / det;
            u = std::max(T(0), std::min(T(1), u + du));
            v = std::max(T(0), std::min(T(1), v + dv));
            if (du * du + dv * dv < T(1e-16)) break;
        }
        return gsPoint2D<T>(u, v);
    };

    const int nH = static_cast<int>(holes.size());
    struct HoleUVBox {
        T u0 = 1, u1 = 0, v0 = 1, v1 = 0;
        T uc = 0.5, vc = 0.5;
        std::vector<gsPoint2D<T>> poly;
    };
    std::vector<HoleUVBox> hBoxes(nH);
    for (int h = 0; h < nH; ++h) {
        const auto& hl = *holes[h];
        for (int c = 0; c < hl.size(); ++c) {
            gsMatrix<T> pr = hl.curve(c).parameterRange();
            for (int k = 0; k < 8; ++k) {
                T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(k) / T(8));
                gsPoint2D<T> p = evalCurvePoint(hl.curve(c), t);
                hBoxes[h].poly.push_back(p);
                gsPoint2D<T> uv = invertPointUV(p);
                hBoxes[h].u0 = std::min(hBoxes[h].u0, uv.x);
                hBoxes[h].u1 = std::max(hBoxes[h].u1, uv.x);
                hBoxes[h].v0 = std::min(hBoxes[h].v0, uv.y);
                hBoxes[h].v1 = std::max(hBoxes[h].v1, uv.y);
            }
        }
        hBoxes[h].uc = T(0.5) * (hBoxes[h].u0 + hBoxes[h].u1);
        hBoxes[h].vc = T(0.5) * (hBoxes[h].v0 + hBoxes[h].v1);
    }

    // 1. Find all gap-midpoint candidate separator cuts in u and v that do not pierce any hole
    std::vector<std::pair<T, T>> uIvals(nH), vIvals(nH);
    std::vector<gsPoint2D<T>> uvCenters(nH);
    for (int h = 0; h < nH; ++h) {
        uIvals[h] = {hBoxes[h].u0, hBoxes[h].u1};
        vIvals[h] = {hBoxes[h].v0, hBoxes[h].v1};
        uvCenters[h] = {hBoxes[h].uc, hBoxes[h].vc};
    }
    std::vector<T> uCuts = SeparatorGrid1D<T>::findClearCuts(uIvals, T(0.01));
    std::vector<T> vCuts = SeparatorGrid1D<T>::findClearCuts(vIvals, T(0.01));

    SeparatorGrid1D<T>::pruneRedundantCuts(uCuts, vCuts, uIvals, vIvals, uvCenters);

    if (SeparatorGrid1D<T>::gridMaxItems(uCuts, vCuts, uvCenters) <= 1 && (!uCuts.empty() || !vCuts.empty())) {
        const int nx = static_cast<int>(uCuts.size()) + 1;
        const int ny = static_cast<int>(vCuts.size()) + 1;
        std::vector<T> us = {T(0)};
        us.insert(us.end(), uCuts.begin(), uCuts.end());
        us.push_back(T(1));
        std::vector<T> vs = {T(0)};
        vs.insert(vs.end(), vCuts.begin(), vCuts.end());
        vs.push_back(T(1));

        std::vector<CurveSeg<T>> sSegs = splitSegMulti(oS, uCuts);
        std::vector<CurveSeg<T>> nSegs = splitSegMulti(oN, uCuts);
        std::vector<CurveSeg<T>> wSegs = splitSegMulti(oW, vCuts);
        std::vector<CurveSeg<T>> eSegs = splitSegMulti(oE, vCuts);

        // Precompute exact interior vertical cuts (at u = us[i] for i=1..nx-1) and horizontal cuts
        std::vector<std::vector<CurveSeg<T>>> vEdges(nx + 1, std::vector<CurveSeg<T>>(ny));
        std::vector<std::vector<CurveSeg<T>>> hEdges(ny + 1, std::vector<CurveSeg<T>>(nx));
        for (int j = 0; j < ny; ++j) {
            vEdges[0][j]  = wSegs[j];
            vEdges[nx][j] = eSegs[j];
        }
        for (int i = 0; i < nx; ++i) {
            hEdges[0][i]  = sSegs[i];
            hEdges[ny][i] = nSegs[i];
        }
        for (int i = 1; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                vEdges[i][j] = extractBlockCutCurve(*refPatch, us[i], vs[j], us[i], vs[j + 1], degree);
        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                hEdges[j][i] = extractBlockCutCurve(*refPatch, us[i], vs[j], us[i + 1], vs[j], degree);

        // Assign holes to grid cells
        std::vector<const gsCurveLoop<T>*> cellHole(nx * ny, nullptr);
        for (int h = 0; h < nH; ++h) {
            int ix = static_cast<int>(std::upper_bound(uCuts.begin(), uCuts.end(), hBoxes[h].uc) - uCuts.begin());
            int iy = static_cast<int>(std::upper_bound(vCuts.begin(), vCuts.end(), hBoxes[h].vc) - vCuts.begin());
            cellHole[iy * nx + ix] = holes[h];
        }

        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const auto& cBot = hEdges[iy][ix];
                const auto& cRgt = vEdges[ix + 1][iy];
                const auto& cTop = hEdges[iy + 1][ix];
                const auto& cLft = vEdges[ix][iy];
                if (cellHole[iy * nx + ix]) {
                    buildBlockQuadRing(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, *cellHole[iy * nx + ix], mp, degree);
                } else {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsPatchFrom4Curves(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, degree)));
                }
            }
        }
        return true;
    }

    // 2. For 2 holes with overlapping axis bboxes, try tilted straight cut in (u, v) space
    if (nH == 2) {
        const T cands[5] = {T(0.25), T(0.38), T(0.50), T(0.62), T(0.75)};
        for (int axis = 0; axis < 2; ++axis) {
            for (T a0 : cands) {
                for (T a1 : cands) {
                    auto ptAt = [&](T u, T v) -> gsPoint2D<T> {
                        gsMatrix<T> uv(2, 1), out; uv(0, 0) = u; uv(1, 0) = v;
                        refPatch->eval_into(uv, out);
                        return gsPoint2D<T>(out(0, 0), out(1, 0));
                    };
                    gsPoint2D<T> A = (axis == 0) ? ptAt(a0, T(0)) : ptAt(T(0), a0);
                    gsPoint2D<T> B = (axis == 0) ? ptAt(a1, T(1)) : ptAt(T(1), a1);
                    auto sideSign = [&](const std::vector<gsPoint2D<T>>& poly) -> int {
                        int sgn = 0;
                        for (const auto& p : poly) {
                            T cr = (B.x - A.x) * (p.y - A.y) - (B.y - A.y) * (p.x - A.x);
                            int s = (cr > T(1e-3)) ? 1 : ((cr < -T(1e-3)) ? -1 : 0);
                            if (s == 0) return 0;
                            if (sgn == 0) sgn = s;
                            else if (s != sgn) return 0;
                        }
                        return sgn;
                    };
                    int s0 = sideSign(hBoxes[0].poly), s1 = sideSign(hBoxes[1].poly);
                    if (s0 != 0 && s1 != 0 && s0 != s1) {
                        if (axis == 0) {
                            auto sSp = oS.splitAtParam(a0);
                            auto nSp = oN.splitAtParam(a1);
                            CurveSeg<T> cut = extractBlockCutCurve(*refPatch, a0, T(0), a1, T(1), degree);
                            const gsCurveLoop<T>* hLeft  = (s0 > 0) ? holes[0] : holes[1];
                            const gsCurveLoop<T>* hRight = (s0 > 0) ? holes[1] : holes[0];
                            buildBlockQuadRing(*sSp.first.bsp, *cut.bsp, *nSp.first.bsp, *oW.bsp, *hLeft, mp, degree);
                            buildBlockQuadRing(*sSp.second.bsp, *oE.bsp, *nSp.second.bsp, *cut.bsp, *hRight, mp, degree);
                        } else {
                            auto wSp = oW.splitAtParam(a0);
                            auto eSp = oE.splitAtParam(a1);
                            CurveSeg<T> cut = extractBlockCutCurve(*refPatch, T(0), a0, T(1), a1, degree);
                            const gsCurveLoop<T>* hBot = (s0 < 0) ? holes[0] : holes[1];
                            const gsCurveLoop<T>* hTop = (s0 < 0) ? holes[1] : holes[0];
                            buildBlockQuadRing(*oS.bsp, *eSp.first.bsp, *cut.bsp, *wSp.first.bsp, *hBot, mp, degree);
                            buildBlockQuadRing(*cut.bsp, *eSp.second.bsp, *oN.bsp, *wSp.second.bsp, *hTop, mp, degree);
                        }
                        return true;
                    }
                }
            }
        }
    }

    // 3. For 3 holes in a non-grid triangle (e.g. plate_4corner_3holes), try 3-way Steiner partition into 3 quads
    if (nH == 3) {
        // Rotate 4 block sides so P[k] is central-split corner
        CurveSeg<T> edgesCCW[4] = { oS, oE, oN.reversed(), oW.reversed() };
        gsPoint2D<T> ctr[3];
        for (int h = 0; h < 3; ++h) {
            ctr[h] = polygonCentroid(hBoxes[h].poly);
        }
        gsPoint2D<T> O_base((ctr[0].x + ctr[1].x + ctr[2].x) / T(3),
                            (ctr[0].y + ctr[1].y + ctr[2].y) / T(3));
        for (int k = 0; k < 4; ++k) {
            const CurveSeg<T>& eK    = edgesCCW[k];       // P_k -> P_{k+1}
            const CurveSeg<T>& eNext = edgesCCW[(k+1)%4]; // P_{k+1} -> P_{k+2}
            const CurveSeg<T>& ePrev = edgesCCW[(k+2)%4]; // P_{k+2} -> P_{k+3}
            const CurveSeg<T>& eLast = edgesCCW[(k+3)%4]; // P_{k+3} -> P_k
            gsPoint2D<T> P_opp = eNext.p1(); // P_{k+2}

            for (T t1 : {T(0.4), T(0.5), T(0.6)}) {
                for (T t2 : {T(0.4), T(0.5), T(0.6)}) {
                    auto spA = eK.splitAtParam(t1);     // P_k -> M1 -> P_{k+1}
                    auto spB = eLast.splitAtParam(t2);  // P_{k+3} -> M2 -> P_k
                    gsPoint2D<T> M1 = spA.first.p1(), M2 = spB.first.p1();
                    for (T w : {T(0.0), T(0.15)}) {
                        gsPoint2D<T> O((T(1) - w) * O_base.x + w * P_opp.x,
                                       (T(1) - w) * O_base.y + w * P_opp.y);
                        bool oInHole = false;
                        for (int h = 0; h < 3; ++h)
                            if (pointInPolyStrict(O, hBoxes[h].poly, T(1e-4))) { oInHole = true; break; }
                        if (oInHole) continue;

                        std::vector<gsPoint2D<T>> qPoly[3] = {
                            {spA.first.p0(), M1, O, M2},
                            {M1, eK.p1(), P_opp, O},
                            {O, P_opp, eLast.p0(), M2}
                        };
                        gsPoint2D<T> qq[3][4];
                        bool convOk = true;
                        for (int q = 0; q < 3; ++q) {
                            for (int v = 0; v < 4; ++v) qq[q][v] = qPoly[q][v];
                            if (!isStrictlyConvexQuad(qq[q])) { convOk = false; break; }
                        }
                        if (!convOk) continue;

                        int match[3] = {-1, -1, -1};
                        bool allClear = true;
                        for (int q = 0; q < 3; ++q) {
                            for (int h = 0; h < 3; ++h) {
                                bool allIn = true;
                                for (const auto& hp : hBoxes[h].poly)
                                    if (!pointInPolyStrict(hp, qPoly[q], T(1e-4))) { allIn = false; break; }
                                if (allIn) {
                                    if (match[q] != -1) allClear = false;
                                    match[q] = h;
                                }
                            }
                            if (match[q] == -1) { allClear = false; break; }
                        }
                        if (!allClear || match[0] == match[1] || match[0] == match[2] || match[1] == match[2])
                            continue;

                        std::unique_ptr<gsBSpline<T>> s1(makeStraightBSpline(M1, O, degree));
                        std::unique_ptr<gsBSpline<T>> s2(makeStraightBSpline(O, P_opp, degree));
                        std::unique_ptr<gsBSpline<T>> s3(makeStraightBSpline(M2, O, degree));
                        std::unique_ptr<gsBSpline<T>> s3_rev(s3->clone()); safeReverseCurve(*s3_rev);
                        std::unique_ptr<gsBSpline<T>> s2_rev(s2->clone()); safeReverseCurve(*s2_rev);

                        buildBlockQuadRing(*spA.first.bsp, *s1, *s3_rev, *spB.second.reversed().bsp, *holes[match[0]], mp, degree);
                        buildBlockQuadRing(*spA.second.bsp, *eNext.bsp, *s2_rev, *s1, *holes[match[1]], mp, degree);
                        buildBlockQuadRing(*s2, *ePrev.bsp, *spB.first.bsp, *s3, *holes[match[2]], mp, degree);
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

// Decompose a 4-sided outer domain [oCCW] enclosing >= 2 holes that separate cleanly along
// physical XY axes inside an axis-aligned inner box Q, joined to oCCW by a conforming O-grid frame.
template <class T>
inline bool build4SidedInnerAxisGrid(
    const std::vector<CurveSeg<T>>& oCCW,
    const std::vector<const gsCurveLoop<T>*>& holes,
    gsMultiPatch<T>& mp,
    int degree)
{
    if (oCCW.size() != 4 || holes.size() < 2) return false;
    const int nH = static_cast<int>(holes.size());

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    std::vector<gsPoint2D<T>> outerPoly;
    outerPoly.reserve(32);
    T oxMin = T(1e18), oxMax = -T(1e18), oyMin = T(1e18), oyMax = -T(1e18);
    for (int s = 0; s < 4; ++s) {
        for (int k = 0; k < 8; ++k) {
            gsPoint2D<T> q = evalCurvePoint(*oCCW[s].bsp, static_cast<T>(k) / T(8));
            outerPoly.push_back(q);
            oxMin = std::min(oxMin, q.x); oxMax = std::max(oxMax, q.x);
            oyMin = std::min(oyMin, q.y); oyMax = std::max(oyMax, q.y);
        }
    }
    const T outerDiag = std::hypot(oxMax - oxMin, oyMax - oyMin);
    if (outerDiag <= T(1e-12)) return false;

    struct HoleXYBox { T x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18, xc = 0, yc = 0; };
    std::vector<HoleXYBox> hBoxes(nH);
    T hxMin = T(1e18), hxMax = -T(1e18), hyMin = T(1e18), hyMax = -T(1e18);
    for (int h = 0; h < nH; ++h) {
        const auto& hl = *holes[h];
        for (int c = 0; c < hl.size(); ++c) {
            gsMatrix<T> pr = hl.curve(c).parameterRange();
            for (int k = 0; k < 8; ++k) {
                T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(k) / T(8));
                gsPoint2D<T> p = evalCurvePoint(hl.curve(c), t);
                hBoxes[h].x0 = std::min(hBoxes[h].x0, p.x); hBoxes[h].x1 = std::max(hBoxes[h].x1, p.x);
                hBoxes[h].y0 = std::min(hBoxes[h].y0, p.y); hBoxes[h].y1 = std::max(hBoxes[h].y1, p.y);
            }
        }
        hBoxes[h].xc = T(0.5) * (hBoxes[h].x0 + hBoxes[h].x1);
        hBoxes[h].yc = T(0.5) * (hBoxes[h].y0 + hBoxes[h].y1);
        hxMin = std::min(hxMin, hBoxes[h].x0); hxMax = std::max(hxMax, hBoxes[h].x1);
        hyMin = std::min(hyMin, hBoxes[h].y0); hyMax = std::max(hyMax, hBoxes[h].y1);
    }
    const T hSpan = std::max(hxMax - hxMin, hyMax - hyMin);
    if (hSpan <= T(1e-12)) return false;

    std::vector<std::pair<T, T>> xIvals(nH), yIvals(nH);
    std::vector<gsPoint2D<T>> xyCenters(nH);
    T loAllX = T(1e18), hiAllX = -T(1e18);
    T loAllY = T(1e18), hiAllY = -T(1e18);
    for (int h = 0; h < nH; ++h) {
        xIvals[h] = {hBoxes[h].x0, hBoxes[h].x1};
        yIvals[h] = {hBoxes[h].y0, hBoxes[h].y1};
        xyCenters[h] = {hBoxes[h].xc, hBoxes[h].yc};
        loAllX = std::min(loAllX, hBoxes[h].x0); hiAllX = std::max(hiAllX, hBoxes[h].x1);
        loAllY = std::min(loAllY, hBoxes[h].y0); hiAllY = std::max(hiAllY, hBoxes[h].y1);
    }
    T minGapX = T(0.02) * std::max(hiAllX - loAllX, T(1e-6));
    T minGapY = T(0.02) * std::max(hiAllY - loAllY, T(1e-6));
    std::vector<T> xCuts = SeparatorGrid1D<T>::findClearCuts(xIvals, minGapX);
    std::vector<T> yCuts = SeparatorGrid1D<T>::findClearCuts(yIvals, minGapY);

    SeparatorGrid1D<T>::pruneRedundantCuts(xCuts, yCuts, xIvals, yIvals, xyCenters);

    if (SeparatorGrid1D<T>::gridMaxItems(xCuts, yCuts, xyCenters) > 1 || (xCuts.empty() && yCuts.empty()))
        return false;

    const int nx = static_cast<int>(xCuts.size()) + 1;
    const int ny = static_cast<int>(yCuts.size()) + 1;
    const int P_in = 2 * (nx + ny);

    auto splitOptsForCount = [](int n) -> std::vector<std::vector<T>> {
        if (n <= 1) return {{}};
        if (n == 2) return {{T(0.50)}, {T(0.35)}, {T(0.65)}};
        if (n == 3) return {
            {T(1) / T(3), T(2) / T(3)},
            {T(0.18), T(0.50)},
            {T(0.50), T(0.82)},
            {T(0.25), T(0.75)}
        };
        std::vector<T> uni(n - 1);
        for (int k = 1; k < n; ++k) uni[k - 1] = static_cast<T>(k) / static_cast<T>(n);
        return {uni};
    };

    auto ptSegDist = [](const gsPoint2D<T>& p, const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> T {
        T dx = b.x - a.x, dy = b.y - a.y;
        T L2 = dx * dx + dy * dy;
        T t = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2)) : T(0);
        return std::hypot(a.x + t * dx - p.x, a.y + t * dy - p.y);
    };
    auto minDistToOuter = [&](const gsPoint2D<T>& p) -> T {
        T mn = T(1e18);
        const int no = static_cast<int>(outerPoly.size());
        for (int i = 0; i < no; ++i)
            mn = std::min(mn, ptSegDist(p, outerPoly[i], outerPoly[(i + 1) % no]));
        return mn;
    };

    bool foundFrame = false;
    T bestFrameScore = T(1e18);
    std::vector<std::unique_ptr<gsGeometry<T>>> bestFrame;
    std::vector<T> bestXs, bestYs;

    for (T padFrac : {T(0.075), T(0.095), T(0.055), T(0.12)}) {
        T pad = padFrac * hSpan;
        T qx0 = hxMin - pad, qx1 = hxMax + pad;
        T qy0 = hyMin - pad, qy1 = hyMax + pad;
        gsPoint2D<T> qc[4] = {{qx0, qy0}, {qx1, qy0}, {qx1, qy1}, {qx0, qy1}};
        bool boxClear = true;
        for (int c = 0; c < 4; ++c) {
            if (!pointInPolyStrict(qc[c], outerPoly, T(1e-3)) ||
                minDistToOuter(qc[c]) < T(0.025) * outerDiag)
            {
                boxClear = false; break;
            }
            gsPoint2D<T> mid(T(0.5) * (qc[c].x + qc[(c + 1) % 4].x),
                             T(0.5) * (qc[c].y + qc[(c + 1) % 4].y));
            if (!pointInPolyStrict(mid, outerPoly, T(1e-3)) ||
                minDistToOuter(mid) < T(0.025) * outerDiag)
            {
                boxClear = false; break;
            }
        }
        if (!boxClear) continue;

        std::vector<T> xs = {qx0};
        xs.insert(xs.end(), xCuts.begin(), xCuts.end());
        xs.push_back(qx1);
        std::vector<T> ys = {qy0};
        ys.insert(ys.end(), yCuts.begin(), yCuts.end());
        ys.push_back(qy1);

        std::vector<CurveSeg<T>> iSeg;
        iSeg.reserve(P_in);
        for (int i = 0; i < nx; ++i)      iSeg.push_back(mkLine({xs[i], qy0},     {xs[i + 1], qy0}));
        for (int j = 0; j < ny; ++j)      iSeg.push_back(mkLine({qx1,   ys[j]},   {qx1, ys[j + 1]}));
        for (int i = nx; i > 0; --i)      iSeg.push_back(mkLine({xs[i], qy1},     {xs[i - 1], qy1}));
        for (int j = ny; j > 0; --j)      iSeg.push_back(mkLine({qx0,   ys[j]},   {qx0, ys[j - 1]}));

        std::vector<gsPoint2D<T>> qPoly = {qc[0], qc[1], qc[2], qc[3]};

        std::vector<std::array<int, 4>> sideModes = {{nx, ny, nx, ny}};
        if (nx != ny) sideModes.push_back({ny, nx, ny, nx});

        for (const auto& nSide : sideModes) {
            auto opts0 = splitOptsForCount(nSide[0]);
            auto opts1 = splitOptsForCount(nSide[1]);
            for (const auto& ts0 : opts0) {
                auto s0 = splitSegMulti(oCCW[0], ts0);
                auto s2 = splitSegMulti(oCCW[2], ts0);
                for (const auto& ts1 : opts1) {
                    auto s1 = splitSegMulti(oCCW[1], ts1);
                    auto s3 = splitSegMulti(oCCW[3], ts1);
                    std::vector<CurveSeg<T>> oSeg;
                    oSeg.reserve(P_in);
                    oSeg.insert(oSeg.end(), s0.begin(), s0.end());
                    oSeg.insert(oSeg.end(), s1.begin(), s1.end());
                    oSeg.insert(oSeg.end(), s2.begin(), s2.end());
                    oSeg.insert(oSeg.end(), s3.begin(), s3.end());
                    if (static_cast<int>(oSeg.size()) != P_in) continue;

                    for (int shift = 0; shift < P_in; ++shift) {
                        std::vector<std::unique_ptr<gsGeometry<T>>> frameCand;
                        frameCand.reserve(P_in);
                        for (int k = 0; k < P_in; ++k) {
                            const auto& oE = oSeg[k];
                            const auto& iE = iSeg[(k + shift) % P_in];
                            CurveSeg<T> spL = mkLine(oE.p0(), iE.p0());
                            CurveSeg<T> spR = mkLine(oE.p1(), iE.p1());
                            frameCand.emplace_back(makeCoonsPatchFrom4Curves(*oE.bsp, *spR.bsp, *iE.bsp, *spL.bsp, degree));
                        }
                        T sc = 0, md = 0, ma = 0;
                        if (evaluateRingPatches(frameCand, qPoly, sc, md, ma)) {
                            if (!foundFrame || sc < bestFrameScore) {
                                foundFrame = true;
                                bestFrameScore = sc;
                                bestFrame = std::move(frameCand);
                                bestXs = xs;
                                bestYs = ys;
                            }
                        }
                    }
                }
            }
        }
        if (foundFrame) break;
    }

    if (!foundFrame) return false;

    for (auto& p : bestFrame)
        mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));

    std::vector<const gsCurveLoop<T>*> cellHole(nx * ny, nullptr);
    for (int h = 0; h < nH; ++h) {
        int ix = static_cast<int>(std::upper_bound(xCuts.begin(), xCuts.end(), hBoxes[h].xc) - xCuts.begin());
        int iy = static_cast<int>(std::upper_bound(yCuts.begin(), yCuts.end(), hBoxes[h].yc) - yCuts.begin());
        cellHole[iy * nx + ix] = holes[h];
    }

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            CurveSeg<T> cBot = mkLine({bestXs[ix],     bestYs[iy]},     {bestXs[ix + 1], bestYs[iy]});
            CurveSeg<T> cRgt = mkLine({bestXs[ix + 1], bestYs[iy]},     {bestXs[ix + 1], bestYs[iy + 1]});
            CurveSeg<T> cTop = mkLine({bestXs[ix],     bestYs[iy + 1]}, {bestXs[ix + 1], bestYs[iy + 1]});
            CurveSeg<T> cLft = mkLine({bestXs[ix],     bestYs[iy]},     {bestXs[ix],     bestYs[iy + 1]});
            if (cellHole[iy * nx + ix]) {
                buildBlockQuadRing(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, *cellHole[iy * nx + ix], mp, degree);
            } else {
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsPatchFrom4Curves(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, degree)));
            }
        }
    }
    return true;
}


// ============================================================================
// Generalized Smooth Offset O-Grid, K-Sector Hole Ring, and Recursive Cell Splitter
// ============================================================================

template <class T>
inline bool isPolySimpleCCW(const std::vector<gsPoint2D<T>>& poly, T minEdgeLen = T(1e-9))
{
    const int n = static_cast<int>(poly.size());
    if (n < 3) return false;
    T area2 = T(0);
    for (int i = 0; i < n; ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % n];
        if (a.dist(b) <= minEdgeLen) return false;
        area2 += a.x * b.y - b.x * a.y;
    }
    if (area2 <= T(1e-10)) return false;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;
            if (segProperIntersect(poly[i], poly[(i + 1) % n], poly[j], poly[(j + 1) % n]))
                return false;
        }
    }
    return true;
}


// Verify candidate patches without full domain boundary Hausdorff search (checks single-sign detJ,
// repairing folded Coons interiors via Greville-Harmonic/Winslow, corner angles, and hole/cell clearance).
template <class T>
inline bool verifyPatchSetLocal(
    std::vector<std::unique_ptr<gsGeometry<T>>>& patches,
    const std::vector<std::vector<gsPoint2D<T>>>& holePolys,
    const std::vector<gsPoint2D<T>>* cellPolyPtr,
    T minAngTol,
    T maxAngTol,
    T& totalEnergyOut)
{
    if (patches.empty()) return false;
    totalEnergyOut = T(0);
    const T minEff = std::min(minAngTol, T(15));
    const T maxEff = std::max(maxAngTol, T(163.5));

    gsMatrix<T> uInner(2, 9), valInner;
    const T pts3[3] = {T(0.22), T(0.5), T(0.78)};
    int col = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            uInner(0, col) = pts3[a];
            uInner(1, col) = pts3[b];
            ++col;
        }

    for (auto& gp : patches) {
        if (!gp) return false;
        if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(gp.get())) {
            T mn = 0, mx = 0, ma = 0, en = 0;
            evalPatchJacobianAndWinslow(*tb, 21, mn, mx, ma, en);
            if (mx < -T(1e-7)) {
                safeReverseTensorBSpline(*tb, 0);
                evalPatchJacobianAndWinslow(*tb, 21, mn, mx, ma, en);
            }
            if (mn <= T(1e-7)) {
                if (!quickHarmonicPatchRepair(*tb, T(1e-7)))
                    return false;
                evalPatchJacobianAndWinslow(*tb, 21, mn, mx, ma, en);
                if (mn <= T(1e-7)) return false;
            }
            totalEnergyOut += en;
        } else {
            T mn = 0, mx = 0, ma = 0, en = 0;
            evalPatchJacobianAndWinslow(*gp, 21, mn, mx, ma, en);
            if (mn <= T(1e-7) && mx >= -T(1e-7)) return false;
            totalEnergyOut += en;
        }

        T angs[4], mnA = 0, mxA = 180;
        bool degen = false;
        evalPatchCornerAnglesDeg(*gp, angs, mnA, mxA, degen);
        if (degen || mnA < minEff || mxA > maxEff) return false;

        if (!holePolys.empty()) {
            gp->eval_into(uInner, valInner);
            for (int k = 0; k < 9; ++k) {
                gsPoint2D<T> q(valInner(0, k), valInner(1, k));
                for (const auto& hp : holePolys) {
                    if (pointInPolyStrict(q, hp, T(1e-7)))
                        return false;
                }
            }
        }
    }
    return true;
}

template <class T>
inline bool buildHoleSmoothingCollar(
    const gsCurveLoop<T>& hole,
    const std::vector<gsPoint2D<T>>* enclosingPolyPtr,
    int degree,
    T minAngTol,
    T maxAngTol,
    std::vector<std::unique_ptr<gsGeometry<T>>>& collarPatches,
    std::vector<CurveSeg<T>>& outerCollarCCW,
    std::vector<int>& outerEdgeToPatch)
{
    collarPatches.clear();
    outerCollarCCW.clear();
    outerEdgeToPatch.clear();

    std::vector<CurveSeg<T>> hRaw = orientLoopSegments(hole, degree, true);
    if (hRaw.empty()) return false;

    std::vector<CurveSeg<T>> hWork;
    for (const auto& s : hRaw) {
        T turn = curveIntegratedTurnRad(*s.bsp, 24);
        if (turn > T(1.95)) {
            int nSub = std::min(12, std::max(2, static_cast<int>(std::ceil(turn / T(1.30)))));
            CurveSeg<T> rem = s;
            for (int k = 0; k < nSub - 1; ++k) {
                auto sp = rem.splitArc(T(1) / static_cast<T>(nSub - k));
                hWork.push_back(sp.first);
                rem = sp.second;
            }
            hWork.push_back(rem);
        } else {
            hWork.push_back(s);
        }
    }

    const int Mw = static_cast<int>(hWork.size());
    if (Mw < 2) return false;

    std::vector<gsPoint2D<T>> holePoly;
    holePoly.reserve(Mw * 12);
    for (const auto& s : hWork)
        for (int k = 0; k < 12; ++k)
            holePoly.push_back(evalCurvePoint(*s.bsp, static_cast<T>(k) / T(12)));
    if (holePoly.size() < 4) return false;
    {
        T a2 = 0;
        for (size_t i = 0; i < holePoly.size(); ++i) {
            const auto& p = holePoly[i];
            const auto& q = holePoly[(i + 1) % holePoly.size()];
            a2 += p.x * q.y - q.x * p.y;
        }
        if (a2 < 0) std::reverse(holePoly.begin(), holePoly.end());
    }

    const T pi = T(3.14159265358979323846);
    const T radToDeg = T(180) / pi;
    std::vector<gsPoint2D<T>> H(Mw);
    std::vector<T> aStart(Mw, T(0)), dCCW(Mw, pi), hTheta(Mw, T(180)), adjLen(Mw, T(1));
    for (int i = 0; i < Mw; ++i) {
        H[i] = hWork[i].p0();
        const auto& prevS = hWork[(i + Mw - 1) % Mw];
        const auto& nextS = hWork[i];
        gsPoint2D<T> dummy(1, 0);
        T aDir = 0, minEL = 0;
        if (computeExteriorCone(prevS, nextS, dummy, aStart[i], dCCW[i], aDir, minEL))
            hTheta[i] = dCCW[i] * radToDeg;
        T l0 = prevS.p0().dist(evalCurvePoint(*prevS.bsp, T(0.5))) + evalCurvePoint(*prevS.bsp, T(0.5)).dist(prevS.p1());
        T l1 = nextS.p0().dist(evalCurvePoint(*nextS.bsp, T(0.5))) + evalCurvePoint(*nextS.bsp, T(0.5)).dist(nextS.p1());
        adjLen[i] = std::max(T(1e-6), std::min(l0, l1));
    }

    auto rayClr = [&](const gsPoint2D<T>& orig, T ang, T minEdge) -> T {
        gsPoint2D<T> dir(std::cos(ang), std::sin(ang));
        T eps = std::min(T(1e-4) * minEdge, T(1e-5));
        gsPoint2D<T> oEps(orig.x + eps * dir.x, orig.y + eps * dir.y);
        T dH = rayDistToPoly(oEps, dir, holePoly);
        T dE = enclosingPolyPtr ? rayDistToPoly(oEps, dir, *enclosingPolyPtr) : T(1e18);
        return std::min(dH, dE);
    };

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        return makeStraightSeg(a, b, degree);
    };
    const std::vector<std::vector<gsPoint2D<T>>> hVec = {holePoly};

    for (T alpha : {T(0.24), T(0.16), T(0.10), T(0.05), T(0.32)}) {
        struct CornerNode {
            bool splitTip = false;
            gsPoint2D<T> q0, mPt, q1;
        };
        std::vector<CornerNode> cn(Mw);
        std::vector<gsPoint2D<T>> ringPts;
        std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> spokeSegs;
        bool ptsOk = true;

        for (int i = 0; i < Mw && ptsOk; ++i) {
            if (hTheta[i] > T(280)) {
                cn[i].splitTip = true;
                T th0 = aStart[i] + dCCW[i] / T(3);
                T thM = aStart[i] + T(0.5) * dCCW[i];
                T th1 = aStart[i] + (T(2) * dCCW[i]) / T(3);
                T c0 = rayClr(H[i], th0, adjLen[i]);
                T cM = rayClr(H[i], thM, adjLen[i]);
                T c1 = rayClr(H[i], th1, adjLen[i]);
                T d = std::min(alpha * std::min(c0, c1), T(0.36) * adjLen[i]);
                if (d <= T(1e-5)) { ptsOk = false; break; }
                T dM = std::min(T(1.16) * d, T(0.70) * cM);
                if (dM <= T(0.65) * d) { ptsOk = false; break; }
                cn[i].q0  = gsPoint2D<T>(H[i].x + d * std::cos(th0), H[i].y + d * std::sin(th0));
                cn[i].mPt = gsPoint2D<T>(H[i].x + dM * std::cos(thM), H[i].y + dM * std::sin(thM));
                cn[i].q1  = gsPoint2D<T>(H[i].x + d * std::cos(th1), H[i].y + d * std::sin(th1));
                ringPts.push_back(cn[i].q0);
                ringPts.push_back(cn[i].mPt);
                ringPts.push_back(cn[i].q1);
                spokeSegs.push_back({H[i], cn[i].q0});
                spokeSegs.push_back({H[i], cn[i].q1});
            } else {
                cn[i].splitTip = false;
                T th0 = aStart[i] + T(0.5) * dCCW[i];
                T c0 = rayClr(H[i], th0, adjLen[i]);
                T d = std::min(alpha * c0, T(0.38) * adjLen[i]);
                if (d <= T(1e-5)) { ptsOk = false; break; }
                cn[i].q0 = gsPoint2D<T>(H[i].x + d * std::cos(th0), H[i].y + d * std::sin(th0));
                cn[i].q1 = cn[i].q0;
                ringPts.push_back(cn[i].q0);
                spokeSegs.push_back({H[i], cn[i].q0});
            }
        }
        if (!ptsOk) continue;

        for (const auto& rp : ringPts) {
            if (pointInPolyStrict(rp, holePoly, T(1e-7))) { ptsOk = false; break; }
            if (enclosingPolyPtr && !pointInPolyStrict(rp, *enclosingPolyPtr, T(1e-7))) { ptsOk = false; break; }
        }
        if (!ptsOk || !isPolySimpleCCW(ringPts, T(1e-10))) continue;

        for (size_t i = 0; i < spokeSegs.size() && ptsOk; ++i) {
            if (!isHoleSpokeClearOfHole(spokeSegs[i].first, spokeSegs[i].second, holePoly)) { ptsOk = false; break; }
            if (enclosingPolyPtr && !isSegmentInsidePoly(spokeSegs[i].first, spokeSegs[i].second, *enclosingPolyPtr)) { ptsOk = false; break; }
            for (size_t j = i + 1; j < spokeSegs.size() && ptsOk; ++j) {
                if (spokeSegs[i].first.dist(spokeSegs[j].first) < T(1e-11)) continue;
                if (segProperIntersect(spokeSegs[i].first, spokeSegs[i].second, spokeSegs[j].first, spokeSegs[j].second))
                    ptsOk = false;
            }
        }
        const int nR = static_cast<int>(ringPts.size());
        for (int e = 0; e < nR && ptsOk; ++e) {
            const auto& A = ringPts[e];
            const auto& B = ringPts[(e + 1) % nR];
            if (!isHoleSpokeClearOfHole(A, B, holePoly)) { ptsOk = false; break; }
            if (enclosingPolyPtr && !isSegmentInsidePoly(A, B, *enclosingPolyPtr)) { ptsOk = false; break; }
        }
        if (!ptsOk) continue;

        std::vector<std::unique_ptr<gsGeometry<T>>> patches;
        std::vector<CurveSeg<T>> colOut;
        std::vector<int> e2p;
        for (int i = 0; i < Mw; ++i) {
            int inxt = (i + 1) % Mw;
            CurveSeg<T> oFlank = mkLine(cn[i].q1, cn[inxt].q0);
            CurveSeg<T> spR = mkLine(cn[inxt].q0, H[inxt]);
            CurveSeg<T> spL = mkLine(H[i], cn[i].q1);
            int pIdx = static_cast<int>(patches.size());
            patches.emplace_back(makeCoonsFromCCWCycle(oFlank, spR, hWork[i].reversed(), spL, degree));
            colOut.push_back(oFlank);
            e2p.push_back(pIdx);

            if (cn[inxt].splitTip) {
                CurveSeg<T> oA = mkLine(cn[inxt].q0, cn[inxt].mPt);
                CurveSeg<T> oB = mkLine(cn[inxt].mPt, cn[inxt].q1);
                patches.emplace_back(makeCoonsFromCCWCycle(
                    mkLine(H[inxt], cn[inxt].q0), oA, oB, mkLine(cn[inxt].q1, H[inxt]), degree));
                colOut.push_back(oA);
                e2p.push_back(-1);
                colOut.push_back(oB);
                e2p.push_back(-1);
            }
        }

        T totEn = 0;
        if (verifyPatchSetLocal(patches, hVec, enclosingPolyPtr, std::min(minAngTol, T(6)), std::max(maxAngTol, T(163.5)), totEn)) {
            collarPatches = std::move(patches);
            outerCollarCCW = std::move(colOut);
            outerEdgeToPatch = std::move(e2p);
            return true;
        }
    }
    return false;
}

template <class T>
inline bool splitOneCollarFlankForEvenness(
    std::vector<std::unique_ptr<gsGeometry<T>>>& collarPatches,
    std::vector<CurveSeg<T>>& outerCollarCCW,
    std::vector<int>& outerEdgeToPatch,
    int degree)
{
    int bestEi = -1;
    T bestLen = -T(1);
    for (int i = 0; i < static_cast<int>(outerCollarCCW.size()); ++i) {
        int pi = outerEdgeToPatch[i];
        if (pi < 0 || pi >= static_cast<int>(collarPatches.size())) continue;
        T len = outerCollarCCW[i].p0().dist(outerCollarCCW[i].p1());
        if (len > bestLen) { bestLen = len; bestEi = i; }
    }
    if (bestEi < 0) return false;
    int pIdx = outerEdgeToPatch[bestEi];
    auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(collarPatches[pIdx].get());
    if (!tb) return false;

    const CurveSeg<T>& seg = outerCollarCCW[bestEi];
    T uSplit = arcLengthSplitParam(*seg.bsp, T(0.5));
    auto spSeg = seg.splitAtParam(uSplit);

    safeNormalizeKnotVector01(tb->knots(0));
    safeNormalizeKnotVector01(tb->knots(1));
    gsMatrix<T> par(2, 2), pts;
    par << T(0), T(1), T(0), T(0);
    tb->eval_into(par, pts);
    gsPoint2D<T> s0(pts(0, 0), pts(1, 0)), s1(pts(0, 1), pts(1, 1));
    T uPatchSplit = uSplit;
    bool revU = false;
    if (s0.dist(seg.p1()) < T(1e-6) && s1.dist(seg.p0()) < T(1e-6)) {
        uPatchSplit = T(1) - uSplit;
        revU = true;
    } else if (!(s0.dist(seg.p0()) < T(1e-6) && s1.dist(seg.p1()) < T(1e-6))) {
        return false;
    }

    gsTensorBSpline<2, T> left, right;
    tb->splitAt(0, uPatchSplit, left, right);
    safeNormalizeKnotVector01(left.knots(0));
    safeNormalizeKnotVector01(left.knots(1));
    safeNormalizeKnotVector01(right.knots(0));
    safeNormalizeKnotVector01(right.knots(1));
    (void)degree;

    collarPatches[pIdx].reset(new gsTensorBSpline<2, T>(revU ? right : left));
    int newPIdx = static_cast<int>(collarPatches.size());
    collarPatches.emplace_back(new gsTensorBSpline<2, T>(revU ? left : right));

    outerCollarCCW[bestEi] = spSeg.first;
    outerCollarCCW.insert(outerCollarCCW.begin() + bestEi + 1, spSeg.second);
    outerEdgeToPatch.insert(outerEdgeToPatch.begin() + bestEi + 1, newPIdx);
    return true;
}

template <class T>
inline bool meshHoleFreePolygon(
    const std::vector<CurveSeg<T>>& bndSegs,
    const std::vector<bool>& bndIsCut,
    int degree,
    T minAngTol,
    T maxAngTol,
    T domDiag,
    int depth,
    int& solverNodes,
    int maxSolverNodes,
    std::vector<std::unique_ptr<gsGeometry<T>>>& outP)
{
    if (depth > 16 || ++solverNodes > maxSolverNodes || bndSegs.size() != bndIsCut.size())
        return false;
    outP.clear();

    struct HalfSegLocal {
        CurveSeg<T> seg;
        bool isCutEdge = false;
    };
    const T radToDeg = T(180) / T(3.14159265358979323846);
    auto evalVertAng = [&](const CurveSeg<T>& prevS, const CurveSeg<T>& nextS) -> T {
        auto dIn = evalCurveDeriv(*prevS.bsp, T(1));
        auto dOut = evalCurveDeriv(*nextS.bsp, T(0));
        T n1 = std::hypot(dIn.x, dIn.y), n2 = std::hypot(dOut.x, dOut.y);
        if (n1 <= T(1e-14) || n2 <= T(1e-14)) return T(180);
        T cr = (dIn.x * dOut.y - dIn.y * dOut.x) / (n1 * n2);
        T dt = (dIn.x * dOut.x + dIn.y * dOut.y) / (n1 * n2);
        T ang = T(180) - std::atan2(cr, dt) * radToDeg;
        if (ang <= T(0)) ang += T(360);
        return ang;
    };
    auto unitTan = [](const CurveSeg<T>& seg, T u) -> gsPoint2D<T> {
        gsPoint2D<T> d = evalCurveDeriv(*seg.bsp, u);
        T n = std::hypot(d.x, d.y);
        return (n > T(1e-14)) ? gsPoint2D<T>(d.x / n, d.y / n) : gsPoint2D<T>(1, 0);
    };
    auto evalCornerCuts = [&](const gsPoint2D<T>& tInFwd, const gsPoint2D<T>& tOutFwd) -> T {
        T cr = tInFwd.x * tOutFwd.y - tInFwd.y * tOutFwd.x;
        T dt = tInFwd.x * tOutFwd.x + tInFwd.y * tOutFwd.y;
        T ang = T(180) - std::atan2(cr, dt) * radToDeg;
        if (ang <= T(0)) ang += T(360);
        return ang;
    };
    auto sampleLocalPoly = [](const std::vector<HalfSegLocal>& bnd, int nPerSeg = 12) {
        std::vector<gsPoint2D<T>> poly;
        poly.reserve(bnd.size() * nPerSeg);
        for (const auto& hs : bnd)
            for (int k = 0; k < nPerSeg; ++k)
                poly.push_back(evalCurvePoint(*hs.seg.bsp, static_cast<T>(k) / static_cast<T>(nPerSeg)));
        return poly;
    };
    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        return makeStraightSeg(a, b, degree);
    };

    std::vector<HalfSegLocal> work(bndSegs.size());
    for (size_t i = 0; i < bndSegs.size(); ++i)
        work[i] = {bndSegs[i], bndIsCut[i]};

    while (work.size() > 4) {
        const int K = static_cast<int>(work.size());
        int mergeIdx = -1;
        for (int i = 0; i < K; ++i) {
            int ip = (i + K - 1) % K;
            if (!work[ip].isCutEdge && !work[i].isCutEdge) {
                T ang = evalVertAng(work[ip].seg, work[i].seg);
                if (std::abs(ang - T(180)) <= T(14) &&
                    curveIntegratedTurnRad(*work[ip].seg.bsp) + curveIntegratedTurnRad(*work[i].seg.bsp) <= T(1.85))
                {
                    mergeIdx = i;
                    break;
                }
            }
        }
        if (mergeIdx < 0) break;
        int ip = (mergeIdx + K - 1) % K;
        auto jPtr = joinCurveSegsExact(std::vector<CurveSeg<T>>{work[ip].seg, work[mergeIdx].seg}, degree);
        if (!jPtr) break;
        work[ip].seg = CurveSeg<T>(*jPtr);
        work.erase(work.begin() + mergeIdx);
    }

    // Special 4-patch acute-tooth cap for K==3 triangles with 1 cut base edge and <38 deg opposite tip
    if (work.size() == 3) {
        int cutCnt = 0, cutIdx = -1;
        for (int i = 0; i < 3; ++i)
            if (work[i].isCutEdge) { ++cutCnt; cutIdx = i; }
        if (cutCnt == 1) {
            int i0 = cutIdx, i1 = (cutIdx + 1) % 3, i2 = (cutIdx + 2) % 3;
            T tipAng = evalVertAng(work[i1].seg, work[i2].seg);
            if (tipAng < T(38)) {
                const std::vector<gsPoint2D<T>> cPoly3 = sampleLocalPoly(work, 12);
                gsPoint2D<T> r2 = work[i0].seg.p0(), r1 = work[i0].seg.p1(), tip = work[i2].seg.p0();
                auto spA0 = work[i1].seg.splitArc(T(1) / T(3));
                auto spA1 = spA0.second.splitArc(T(0.5));
                CurveSeg<T> sA = spA0.first, sB = spA1.first, sC = spA1.second;
                auto spN = work[i2].seg.splitArc(T(0.5));
                CurveSeg<T> sN0 = spN.first, sN1 = spN.second;
                gsPoint2D<T> M1 = sA.p1(), M2 = sB.p1(), N1 = sN0.p1();

                std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> qCands = {
                    { gsPoint2D<T>(T(0.42) * M1.x + T(0.42) * N1.x + T(0.16) * r2.x,
                                   T(0.42) * M1.y + T(0.42) * N1.y + T(0.16) * r2.y),
                      gsPoint2D<T>(T(0.46) * M2.x + T(0.44) * N1.x + T(0.10) * tip.x,
                                   T(0.46) * M2.y + T(0.44) * N1.y + T(0.10) * tip.y) },
                    { gsPoint2D<T>(T(0.28) * r1.x + T(0.34) * r2.x + T(0.38) * tip.x,
                                   T(0.28) * r1.y + T(0.34) * r2.y + T(0.38) * tip.y),
                      gsPoint2D<T>(T(0.16) * r1.x + T(0.26) * r2.x + T(0.58) * tip.x,
                                   T(0.16) * r1.y + T(0.26) * r2.y + T(0.58) * tip.y) }
                };
                for (const auto& qq : qCands) {
                    const gsPoint2D<T>& Q1 = qq.first;
                    const gsPoint2D<T>& Q2 = qq.second;
                    if (!pointInPolyStrict(Q1, cPoly3, T(1e-6)) || !pointInPolyStrict(Q2, cPoly3, T(1e-6)))
                        continue;
                    std::vector<std::unique_ptr<gsGeometry<T>>> triCap;
                    triCap.emplace_back(makeCoonsFromCCWCycle(sC, sN0, mkLine(N1, Q2), mkLine(Q2, M2), degree));
                    triCap.emplace_back(makeCoonsFromCCWCycle(sB, mkLine(M2, Q2), mkLine(Q2, Q1), mkLine(Q1, M1), degree));
                    triCap.emplace_back(makeCoonsFromCCWCycle(mkLine(Q2, N1), sN1, mkLine(r2, Q1), mkLine(Q1, Q2), degree));
                    triCap.emplace_back(makeCoonsFromCCWCycle(sA, mkLine(M1, Q1), mkLine(Q1, r2), work[i0].seg, degree));
                    T en = 0;
                    if (verifyPatchSetLocal(triCap, {}, &cPoly3, std::min(minAngTol, tipAng * T(0.85)), maxAngTol, en)) {
                        outP = std::move(triCap);
                        return true;
                    }
                }
            }
        }
    }

    if (work.size() % 2 == 1) {
        int splitIdx = -1;
        T maxLen = -T(1);
        for (size_t i = 0; i < work.size(); ++i) {
            if (work[i].isCutEdge) continue;
            T d = work[i].seg.p0().dist(work[i].seg.p1());
            if (d > maxLen) { maxLen = d; splitIdx = static_cast<int>(i); }
        }
        if (splitIdx < 0) return false;
        auto sp = work[splitIdx].seg.splitArc(T(0.5));
        work[splitIdx] = {sp.first, false};
        work.insert(work.begin() + splitIdx + 1, {sp.second, false});
    }

    if (work.size() == 2) {
        if (work[0].isCutEdge || work[1].isCutEdge) return false;
        std::vector<HalfSegLocal> expanded;
        for (int i = 0; i < 2; ++i) {
            auto sp = work[i].seg.splitArc(T(0.5));
            expanded.push_back({sp.first, false});
            expanded.push_back({sp.second, false});
        }
        work = std::move(expanded);
    }

    const int K = static_cast<int>(work.size());
    if (K < 4 || K % 2 != 0) return false;
    const std::vector<gsPoint2D<T>> cPoly = sampleLocalPoly(work, 12);

    if (K == 4) {
        T angs[4];
        bool needSplit180 = false;
        const T maxLeafA = std::min(maxAngTol, T(163.5));
        for (int i = 0; i < 4; ++i) {
            angs[i] = evalVertAng(work[(i + 3) % 4].seg, work[i].seg);
            if (angs[i] > maxLeafA || angs[i] < std::min(minAngTol, T(6)))
                needSplit180 = true;
        }
        if (!needSplit180) {
            std::vector<std::unique_ptr<gsGeometry<T>>> one;
            one.emplace_back(makeCoonsPatchFrom4Curves(
                *work[0].seg.bsp, *work[1].seg.bsp,
                *work[2].seg.reversed().bsp, *work[3].seg.reversed().bsp, degree));
            T en = 0;
            if (verifyPatchSetLocal(one, {}, &cPoly, minAngTol, maxAngTol, en)) {
                outP = std::move(one);
                return true;
            }
        }

        for (int base = 0; base < 2; ++base) {
            int s1 = (base + 1) % 4, s3 = (base + 3) % 4;
            if (work[s1].isCutEdge || work[s3].isCutEdge) continue;
            for (int M : {2, 3, 4}) {
                CurveSeg<T> remR = work[s1].seg;
                CurveSeg<T> remL = work[s3].seg.reversed();
                std::vector<CurveSeg<T>> rSegs(M), lSegs(M);
                for (int j = 0; j < M - 1; ++j) {
                    T f = T(1) / static_cast<T>(M - j);
                    auto sr = remR.splitArc(f); rSegs[j] = sr.first; remR = sr.second;
                    auto sl = remL.splitArc(f); lSegs[j] = sl.first; remL = sl.second;
                }
                rSegs[M - 1] = remR; lSegs[M - 1] = remL;
                std::vector<CurveSeg<T>> hSegs(M + 1);
                hSegs[0] = work[base].seg;
                hSegs[M] = work[(base + 2) % 4].seg.reversed();
                bool okH = true;
                for (int j = 1; j < M; ++j) {
                    if (!isSegmentInsidePoly(lSegs[j].p0(), rSegs[j].p0(), cPoly)) { okH = false; break; }
                    hSegs[j] = mkLine(lSegs[j].p0(), rSegs[j].p0());
                }
                if (!okH) continue;
                std::vector<std::unique_ptr<gsGeometry<T>>> strip;
                for (int j = 0; j < M; ++j) {
                    strip.emplace_back(makeCoonsPatchFrom4Curves(
                        *hSegs[j].bsp, *rSegs[j].bsp, *hSegs[j + 1].bsp, *lSegs[j].bsp, degree));
                }
                T en = 0;
                if (verifyPatchSetLocal(strip, {}, &cPoly, minAngTol, maxAngTol, en)) {
                    outP = std::move(strip);
                    return true;
                }
            }
        }

        gsPoint2D<T> P[4] = {work[0].seg.p0(), work[1].seg.p0(), work[2].seg.p0(), work[3].seg.p0()};
        gsPoint2D<T> ctr = polygonCentroid(cPoly);
        if (!pointInPolyStrict(ctr, cPoly, T(1e-5)))
            ctr = gsPoint2D<T>(T(0.25)*(P[0].x+P[1].x+P[2].x+P[3].x), T(0.25)*(P[0].y+P[1].y+P[2].y+P[3].y));
        if (pointInPolyStrict(ctr, cPoly, T(1e-5))) {
            for (T alpha : {T(0.45), T(0.32), T(0.56)}) {
                gsPoint2D<T> Q[4];
                bool qOk = true;
                for (int k = 0; k < 4; ++k) {
                    Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x), ctr.y + alpha * (P[k].y - ctr.y));
                    if (!pointInPolyStrict(Q[k], cPoly, T(1e-5)) || !isSegmentInsidePoly(Q[k], P[k], cPoly)) {
                        qOk = false; break;
                    }
                }
                if (!qOk || !isStrictlyConvexQuad(Q)) continue;
                std::vector<std::unique_ptr<gsGeometry<T>>> og;
                og.emplace_back(makeCoonsFromCCWCycle(mkLine(Q[0], Q[1]), mkLine(Q[1], Q[2]), mkLine(Q[2], Q[3]), mkLine(Q[3], Q[0]), degree));
                for (int k = 0; k < 4; ++k) {
                    int kn = (k + 1) % 4;
                    og.emplace_back(makeCoonsPatchFrom4Curves(
                        *mkLine(Q[k], Q[kn]).bsp, *mkLine(Q[kn], P[kn]).bsp,
                        *work[k].seg.bsp, *mkLine(Q[k], P[k]).bsp, degree));
                }
                T en = 0;
                if (verifyPatchSetLocal(og, {}, &cPoly, minAngTol, maxAngTol, en)) {
                    outP = std::move(og);
                    return true;
                }
            }
        }
        return false;
    }

    std::vector<T> angs(K);
    std::vector<bool> mustBreak(K, false);
    for (int i = 0; i < K; ++i) {
        angs[i] = evalVertAng(work[(i + K - 1) % K].seg, work[i].seg);
        if (angs[i] > T(154)) mustBreak[i] = true;
    }

    struct Chord { int a, b; CurveSeg<T> seg; T score; };
    std::vector<Chord> chords;
    for (int a = 0; a < K; ++a) {
        gsPoint2D<T> A = work[a].seg.p0();
        gsPoint2D<T> tAin = unitTan(work[(a + K - 1) % K].seg, T(1));
        gsPoint2D<T> tAout = unitTan(work[a].seg, T(0));
        for (int step = 3; step <= K - 3; step += 2) {
            int b = (a + step) % K;
            if (b < a && ((a - b) % 2 == 1) && (K - step == step)) continue;
            gsPoint2D<T> B = work[b].seg.p0();
            if (!isSegmentInsidePoly(A, B, cPoly)) continue;
            CurveSeg<T> cSeg = mkLine(A, B);
            gsPoint2D<T> tCutA = unitTan(cSeg, T(0));
            gsPoint2D<T> tCutB = unitTan(cSeg, T(1));
            gsPoint2D<T> tBin = unitTan(work[(b + K - 1) % K].seg, T(1));
            gsPoint2D<T> tBout = unitTan(work[b].seg, T(0));

            gsPoint2D<T> tCutArev(-tCutA.x, -tCutA.y);
            gsPoint2D<T> tCutBrev(-tCutB.x, -tCutB.y);
            T c1_B = evalCornerCuts(tBin, tCutBrev);
            T c1_A = evalCornerCuts(tCutArev, tAout);
            T c2_A = evalCornerCuts(tAin, tCutA);
            T c2_B = evalCornerCuts(tCutB, tBout);
            T mnC = std::min(std::min(c1_A, c1_B), std::min(c2_A, c2_B));
            T mxC = std::max(std::max(c1_A, c1_B), std::max(c2_A, c2_B));
            if (mnC < std::min(minAngTol, T(14)) || mxC > std::min(maxAngTol, T(162))) continue;

            T sc = std::min(mnC, T(180) - mxC) - T(8) * (A.dist(B) / std::max(T(1e-6), domDiag));
            if (mustBreak[a]) sc += T(45);
            if (mustBreak[b]) sc += T(45);
            chords.push_back({a, b, cSeg, sc});
        }
    }
    std::sort(chords.begin(), chords.end(), [](const Chord& x, const Chord& y) { return x.score > y.score; });

    for (size_t ci = 0; ci < std::min<size_t>(chords.size(), 4); ++ci) {
        int a = chords[ci].a, b = chords[ci].b;
        std::vector<CurveSeg<T>> s1, s2;
        std::vector<bool> c1, c2;
        for (int k = a; k != b; k = (k + 1) % K) { s1.push_back(work[k].seg); c1.push_back(work[k].isCutEdge); }
        s1.push_back(chords[ci].seg.reversed()); c1.push_back(true);
        for (int k = b; k != a; k = (k + 1) % K) { s2.push_back(work[k].seg); c2.push_back(work[k].isCutEdge); }
        s2.push_back(chords[ci].seg); c2.push_back(true);

        std::vector<std::unique_ptr<gsGeometry<T>>> p1, p2;
        if (meshHoleFreePolygon(s1, c1, degree, minAngTol, maxAngTol, domDiag, depth + 1, solverNodes, maxSolverNodes, p1) &&
            meshHoleFreePolygon(s2, c2, degree, minAngTol, maxAngTol, domDiag, depth + 1, solverNodes, maxSolverNodes, p2)) {
            outP = std::move(p1);
            for (auto& g : p2) outP.push_back(std::move(g));
            return true;
        }
    }

    gsPoint2D<T> ctr = polygonCentroid(cPoly);
    if (pointInPolyStrict(ctr, cPoly, T(1e-5))) {
        for (int par = 0; par < 2; ++par) {
            bool visOk = true;
            for (int s = 0; s < K / 2 && visOk; ++s) {
                int idx = (par + 2 * s) % K;
                if (!isSegmentInsidePoly(ctr, work[idx].seg.p0(), cPoly)) visOk = false;
            }
            if (!visOk) continue;
            std::vector<std::unique_ptr<gsGeometry<T>>> kp;
            for (int s = 0; s < K / 2; ++s) {
                int a = (par + 2 * s) % K, b = (a + 1) % K, c = (a + 2) % K;
                kp.emplace_back(makeCoonsPatchFrom4Curves(
                    *mkLine(ctr, work[a].seg.p0()).bsp,
                    *work[a].seg.bsp,
                    *work[b].seg.reversed().bsp,
                    *mkLine(ctr, work[c].seg.p0()).bsp,
                    degree));
            }
            T en = 0;
            if (verifyPatchSetLocal(kp, {}, &cPoly, minAngTol, maxAngTol, en)) {
                outP = std::move(kp);
                return true;
            }
        }
    }
    return false;
}

// Triangulate a simple CCW straight polygon pts[0..N-1] and pair adjacent triangles into quads
// (preferring 0 leftover triangles and maximizing minimum quad corner angle).
template <class T>
inline bool triangulateAndPairSimplePolygon(
    const std::vector<gsPoint2D<T>>& pts,
    T minAngTol,
    T maxAngTol,
    std::vector<std::array<int, 4>>& outQuads,
    std::vector<std::array<int, 3>>& outRemTris,
    T& bestMinAngOut)
{
    const int N = static_cast<int>(pts.size());
    outQuads.clear();
    outRemTris.clear();
    bestMinAngOut = T(0);
    if (N < 3) return false;
    if (N == 3) {
        outRemTris.push_back({0, 1, 2});
        return true;
    }

    auto quadCornerMinMax = [&](int a, int b, int c, int d, T& mn, T& mx) -> bool {
        gsPoint2D<T> q[4] = {pts[a], pts[b], pts[c], pts[d]};
        if (!isStrictlyConvexQuad(q)) return false;
        const T r2d = T(180) / T(3.14159265358979323846);
        mn = T(360); mx = T(0);
        for (int k = 0; k < 4; ++k) {
            const auto& p0 = q[(k + 3) % 4];
            const auto& pk = q[k];
            const auto& p1 = q[(k + 1) % 4];
            T ux = p1.x - pk.x, uy = p1.y - pk.y;
            T vx = p0.x - pk.x, vy = p0.y - pk.y;
            T nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
            if (nu <= T(1e-12) || nv <= T(1e-12)) return false;
            T ang = std::atan2((ux * vy - uy * vx) / (nu * nv), (ux * vx + uy * vy) / (nu * nv)) * r2d;
            if (ang <= T(0)) ang += T(360);
            mn = std::min(mn, ang);
            mx = std::max(mx, ang);
        }
        return (mn >= minAngTol && mx <= maxAngTol);
    };

    bool foundAny = false;
    size_t bestRem = N + 10;
    T bestScore = -T(1e9);

    for (int mode = 0; mode < 3; ++mode) {
        std::vector<int> ring(N);
        for (int i = 0; i < N; ++i) ring[i] = i;
        std::vector<std::array<int, 3>> tris;
        auto ptInTri = [&](int ia, int ib, int ic, int ip) -> bool {
            return cross2(pts[ia], pts[ib], pts[ip]) >= -T(1e-9) &&
                   cross2(pts[ib], pts[ic], pts[ip]) >= -T(1e-9) &&
                   cross2(pts[ic], pts[ia], pts[ip]) >= -T(1e-9);
        };
        bool earOk = true;
        while (ring.size() > 3) {
            int m = static_cast<int>(ring.size());
            int bestI = -1;
            T bestS = -T(1e30);
            for (int i = 0; i < m; ++i) {
                int ia = ring[(i + m - 1) % m], ib = ring[i], ic = ring[(i + 1) % m];
                T a2 = cross2(pts[ia], pts[ib], pts[ic]);
                if (a2 <= T(1e-8)) continue;
                bool clr = true;
                for (int k = 0; k < m; ++k) {
                    int ik = ring[k];
                    if (ik == ia || ik == ib || ik == ic) continue;
                    if (ptInTri(ia, ib, ic, ik)) { clr = false; break; }
                }
                if (!clr) continue;
                T dAC2 = (pts[ia].x - pts[ic].x) * (pts[ia].x - pts[ic].x) +
                         (pts[ia].y - pts[ic].y) * (pts[ia].y - pts[ic].y);
                T sc = T(0);
                if (mode == 0) {
                    sc = T(1) / (dAC2 + T(1e-14));
                } else if (mode == 1) {
                    T l2 = (pts[ia].x - pts[ib].x) * (pts[ia].x - pts[ib].x) + (pts[ia].y - pts[ib].y) * (pts[ia].y - pts[ib].y) +
                           (pts[ib].x - pts[ic].x) * (pts[ib].x - pts[ic].x) + (pts[ib].y - pts[ic].y) * (pts[ib].y - pts[ic].y) + dAC2;
                    sc = a2 / (l2 + T(1e-14));
                } else {
                    T mnA = 0, mxA = 180;
                    triMinMaxAnglesDeg(pts[ia], pts[ib], pts[ic], mnA, mxA);
                    sc = mnA;
                }
                if (sc > bestS) { bestS = sc; bestI = i; }
            }
            if (bestI < 0) { earOk = false; break; }
            int ia = ring[(bestI + m - 1) % m], ib = ring[bestI], ic = ring[(bestI + 1) % m];
            tris.push_back({ia, ib, ic});
            ring.erase(ring.begin() + bestI);
        }
        if (!earOk || ring.size() != 3) continue;
        if (cross2(pts[ring[0]], pts[ring[1]], pts[ring[2]]) <= T(1e-8)) continue;
        tris.push_back({ring[0], ring[1], ring[2]});

        const int nT = static_cast<int>(tris.size());
        // Apply Lawson flips inside simple polygon if mode == 2 (creates alternate dual graph when snake has forks)
        if (mode == 2 && nT >= 2) {
            for (int pass = 0; pass < nT * 2; ++pass) {
                bool flipped = false;
                for (int t0 = 0; t0 < nT && !flipped; ++t0) {
                    for (int k0 = 0; k0 < 3 && !flipped; ++k0) {
                        int a = tris[t0][k0], b = tris[t0][(k0 + 1) % 3], c = tris[t0][(k0 + 2) % 3];
                        for (int t1 = t0 + 1; t1 < nT && !flipped; ++t1) {
                            for (int k1 = 0; k1 < 3; ++k1) {
                                if (tris[t1][k1] == b && tris[t1][(k1 + 1) % 3] == a) {
                                    int d = tris[t1][(k1 + 2) % 3];
                                    gsPoint2D<T> q4[4] = {pts[c], pts[a], pts[d], pts[b]};
                                    if (!isStrictlyConvexQuad(q4)) continue;
                                    T mn0, mx0, mn1, mx1, nm0, nx0, nm1, nx1;
                                    triMinMaxAnglesDeg(pts[a], pts[b], pts[c], mn0, mx0);
                                    triMinMaxAnglesDeg(pts[b], pts[a], pts[d], mn1, mx1);
                                    triMinMaxAnglesDeg(pts[c], pts[a], pts[d], nm0, nx0);
                                    triMinMaxAnglesDeg(pts[c], pts[d], pts[b], nm1, nx1);
                                    if (std::min(nm0, nm1) > std::min(mn0, mn1) + T(0.5)) {
                                        tris[t0] = {c, a, d};
                                        tris[t1] = {c, d, b};
                                        flipped = true;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                if (!flipped) break;
            }
        }

        // Build valid pair edges between adjacent triangles
        struct PE { int t0, t1; std::array<int, 4> q; T sc; };
        std::vector<PE> pedges;
        for (int t0 = 0; t0 < nT; ++t0) {
            for (int k0 = 0; k0 < 3; ++k0) {
                int a = tris[t0][k0], b = tris[t0][(k0 + 1) % 3], c = tris[t0][(k0 + 2) % 3];
                for (int t1 = t0 + 1; t1 < nT; ++t1) {
                    for (int k1 = 0; k1 < 3; ++k1) {
                        if (tris[t1][k1] == b && tris[t1][(k1 + 1) % 3] == a) {
                            int d = tris[t1][(k1 + 2) % 3];
                            T mn = 0, mx = 180;
                            if (quadCornerMinMax(c, a, d, b, mn, mx)) {
                                T eff = std::min(mn, T(180) - mx);
                                pedges.push_back({t0, t1, {c, a, d, b}, eff});
                            }
                        }
                    }
                }
            }
        }

        // Exact DFS maximum matching (nT <= 14 is instantaneous)
        std::vector<std::vector<int>> tAdj(nT);
        for (int pi = 0; pi < static_cast<int>(pedges.size()); ++pi) {
            tAdj[pedges[pi].t0].push_back(pi);
            tAdj[pedges[pi].t1].push_back(pi);
        }
        std::vector<bool> used(nT, false);
        std::vector<int> cur, bestMatch;
        int bestCard = -1;
        T bestMatchSc = -T(1e9);
        int nodes = 0;
        std::function<void(int, T)> dfs = [&](int ti, T minEdgeSc) {
            if (++nodes > 6000) return;
            while (ti < nT && used[ti]) ++ti;
            if (ti == nT) {
                int card = static_cast<int>(cur.size());
                if (card > bestCard || (card == bestCard && minEdgeSc > bestMatchSc)) {
                    bestCard = card;
                    bestMatchSc = minEdgeSc;
                    bestMatch = cur;
                }
                return;
            }
            int rem = 0;
            for (int j = ti; j < nT; ++j) if (!used[j]) ++rem;
            if (static_cast<int>(cur.size()) + rem / 2 < bestCard) return;

            for (int pi : tAdj[ti]) {
                const auto& pe = pedges[pi];
                int other = (pe.t0 == ti) ? pe.t1 : pe.t0;
                if (used[other]) continue;
                used[ti] = used[other] = true;
                cur.push_back(pi);
                dfs(ti + 1, std::min(minEdgeSc, pe.sc));
                cur.pop_back();
                used[ti] = used[other] = false;
                if (nodes > 6000) return;
            }
            used[ti] = true;
            dfs(ti + 1, minEdgeSc);
            used[ti] = false;
        };
        dfs(0, T(180));

        std::vector<bool> matched(nT, false);
        std::vector<std::array<int, 4>> mQuads;
        for (int pi : bestMatch) {
            matched[pedges[pi].t0] = matched[pedges[pi].t1] = true;
            mQuads.push_back(pedges[pi].q);
        }
        std::vector<std::array<int, 3>> rTris;
        for (int ti = 0; ti < nT; ++ti)
            if (!matched[ti]) rTris.push_back(tris[ti]);

        if (!foundAny || rTris.size() < bestRem || (rTris.size() == bestRem && bestMatchSc > bestScore)) {
            foundAny = true;
            bestRem = rTris.size();
            bestScore = bestMatchSc;
            outQuads = std::move(mQuads);
            outRemTris = std::move(rTris);
            bestMinAngOut = (bestCard > 0) ? bestMatchSc : minAngTol;
        }
        if (bestRem == 0 && bestScore > T(28)) break;
    }
    return foundAny;
}

// Smooth / Irregular 0-Hole Offset O-Grid (solves smooth non-star shapes like amoeba_smooth_irregular).
template <class T>
inline bool buildSmoothOffsetOGrid(
    const gsPlanarDomain<T>& dom,
    int degree,
    gsMultiPatch<T>& outMp,
    T minAngTol = T(12),
    T maxAngTol = T(162))
{
    if (dom.numLoops() != 1 || dom.outer().size() == 0) return false;
    std::vector<CurveSeg<T>> oCCW = orientLoopSegments(dom.outer(), degree, true);
    auto joined = joinCurveSegsExact(oCCW, degree);
    if (!joined) return false;
    CurveSeg<T> fullLoop(*joined);

    const std::vector<gsPoint2D<T>> outerPoly = sampleCurveLoopPoly(dom.outer(), 32);
    if (outerPoly.size() < 6) return false;

    const int NS = 256;
    std::vector<T> uGrid(NS + 1), sArc(NS + 1, T(0)), sBend(NS + 1, T(0));
    std::vector<gsPoint2D<T>> sPts(NS + 1), sTans(NS + 1);
    for (int i = 0; i <= NS; ++i) {
        T u = static_cast<T>(i) / static_cast<T>(NS);
        uGrid[i] = u;
        sPts[i] = evalCurvePoint(*fullLoop.bsp, u);
        gsPoint2D<T> d = evalCurveDeriv(*fullLoop.bsp, u);
        T nm = std::hypot(d.x, d.y);
        sTans[i] = (nm > T(1e-14)) ? gsPoint2D<T>(d.x / nm, d.y / nm) : gsPoint2D<T>(1, 0);
    }
    T totalLen = T(0);
    for (int i = 1; i <= NS; ++i) {
        T ds = sPts[i].dist(sPts[i - 1]);
        sArc[i] = sArc[i - 1] + ds;
        totalLen += ds;
    }
    if (totalLen <= T(1e-12)) return false;
    const T approxRad = totalLen / (T(2) * T(3.14159265358979323846));
    for (int i = 1; i <= NS; ++i) {
        T ds = sPts[i].dist(sPts[i - 1]);
        T cr = sTans[i - 1].x * sTans[i].y - sTans[i - 1].y * sTans[i].x;
        T dt = sTans[i - 1].x * sTans[i].x + sTans[i - 1].y * sTans[i].y;
        T dAng = std::abs(std::atan2(cr, dt));
        sBend[i] = sBend[i - 1] + ds + T(0.35) * approxRad * dAng;
    }

    auto paramAtCumFrac = [&](const std::vector<T>& cum, T frac) -> T {
        T f = frac - std::floor(frac);
        if (f < T(0)) f += T(1);
        T target = f * cum.back();
        for (int i = 1; i <= NS; ++i) {
            if (cum[i] >= target) {
                T span = cum[i] - cum[i - 1];
                T al = (span > T(1e-14)) ? (target - cum[i - 1]) / span : T(0);
                return (static_cast<T>(i - 1) + al) / static_cast<T>(NS);
            }
        }
        return T(0);
    };

    bool foundBest = false;
    T bestCost = T(1e30);
    gsMultiPatch<T> bestMp;

    const T minEff = std::min(minAngTol, T(12));
    const T maxEff = std::min(T(163), std::max(maxAngTol, T(158)));

    for (int N : {6, 8, 10, 12, 14, 16, 20, 24}) {
        if (foundBest && bestMp.nPatches() <= static_cast<size_t>(N + N / 2))
            break;

        struct Scheme { const std::vector<T>* cum; T phase; };
        const Scheme schemes[3] = {
            {&sArc,  T(0)},
            {&sArc,  T(0.5) / static_cast<T>(N)},
            {&sBend, T(0)}
        };

        for (const auto& sch : schemes) {
            std::vector<T> uStations(N);
            for (int k = 0; k < N; ++k)
                uStations[k] = paramAtCumFrac(*sch.cum, sch.phase + static_cast<T>(k) / static_cast<T>(N));

            // Shift start of fullLoop to uStations[0] while keeping exact curve geometry
            CurveSeg<T> shifted = fullLoop;
            std::vector<T> relSplits(N - 1);
            T u0 = uStations[0];
            if (u0 > T(1e-4) && u0 < T(1) - T(1e-4)) {
                auto sp0 = fullLoop.splitAtParam(u0);
                auto sj = joinCurveSegsExact(std::vector<CurveSeg<T>>{sp0.second, sp0.first}, degree);
                if (!sj) continue;
                shifted = CurveSeg<T>(*sj);
                for (int k = 1; k < N; ++k) {
                    T uRel = uStations[k] - u0;
                    if (uRel < T(0)) uRel += T(1);
                    relSplits[k - 1] = uRel;
                }
            } else {
                for (int k = 1; k < N; ++k) relSplits[k - 1] = uStations[k];
            }
            std::sort(relSplits.begin(), relSplits.end());
            bool splitsOk = (relSplits.front() > T(0.02) && relSplits.back() < T(0.98));
            for (int k = 1; k < N - 1 && splitsOk; ++k)
                if (relSplits[k] - relSplits[k - 1] < T(0.02)) splitsOk = false;
            if (!splitsOk) continue;

            std::vector<CurveSeg<T>> bEdges = splitSegMulti(shifted, relSplits);
            if (static_cast<int>(bEdges.size()) != N) continue;

            std::vector<gsPoint2D<T>> P(N), inNorm(N);
            std::vector<T> clr(N, T(0)), lfs(N, totalLen), adjLen(N, T(0));
            for (int k = 0; k < N; ++k) P[k] = bEdges[k].p0();
            const int nOP = static_cast<int>(outerPoly.size());
            for (int k = 0; k < N; ++k) {
                const auto& prev = bEdges[(k + N - 1) % N];
                const auto& curr = bEdges[k];
                auto dIn = evalCurveDeriv(*prev.bsp, T(1));
                auto dOut = evalCurveDeriv(*curr.bsp, T(0));
                T n1 = std::hypot(dIn.x, dIn.y), n2 = std::hypot(dOut.x, dOut.y);
                gsPoint2D<T> t1 = (n1 > T(1e-14)) ? gsPoint2D<T>(dIn.x / n1, dIn.y / n1) : gsPoint2D<T>(1, 0);
                gsPoint2D<T> t2 = (n2 > T(1e-14)) ? gsPoint2D<T>(dOut.x / n2, dOut.y / n2) : gsPoint2D<T>(1, 0);
                gsPoint2D<T> ni(-t1.y - t2.y, t1.x + t2.x);
                T nl = std::hypot(ni.x, ni.y);
                inNorm[k] = (nl > T(1e-14)) ? gsPoint2D<T>(ni.x / nl, ni.y / nl) : gsPoint2D<T>(-t2.y, t2.x);
                clr[k] = std::min(rayDistToPoly(P[k], inNorm[k], outerPoly), totalLen);
                adjLen[k] = std::min(P[k].dist(P[(k + 1) % N]), P[k].dist(P[(k + N - 1) % N]));
                // Local feature size: distance to non-adjacent parts of outerPoly so narrow necks/fingers stay simple
                int nearIdx = 0;
                T nearD = T(1e30);
                for (int oi = 0; oi < nOP; ++oi) {
                    T d = P[k].dist(outerPoly[oi]);
                    if (d < nearD) { nearD = d; nearIdx = oi; }
                }
                const int skipWin = std::max(2, nOP / std::max(N, 6));
                T minNonAdj = totalLen;
                for (int oi = 0; oi < nOP; ++oi) {
                    int diff = std::abs(oi - nearIdx);
                    diff = std::min(diff, nOP - diff);
                    if (diff <= skipWin) continue;
                    minNonAdj = std::min(minNonAdj, P[k].dist(outerPoly[oi]));
                }
                lfs[k] = minNonAdj;
            }

            for (T alpha : {T(0.22), T(0.14), T(0.08), T(0.04), T(0.32)}) {
                std::vector<gsPoint2D<T>> Q(N);
                bool qInitOk = true;
                for (int k = 0; k < N; ++k) {
                    T dOff = std::min(alpha * std::min(clr[k], T(0.90) * lfs[k]), T(0.42) * adjLen[k]);
                    Q[k] = gsPoint2D<T>(P[k].x + dOff * inNorm[k].x, P[k].y + dOff * inNorm[k].y);
                    if (!pointInPolyStrict(Q[k], outerPoly, T(1e-5))) { qInitOk = false; break; }
                }
                if (!qInitOk) continue;

                // Laplacian smooth Q while staying inside outerPoly
                for (int sw = 0; sw < 4; ++sw) {
                    std::vector<gsPoint2D<T>> Qnew = Q;
                    for (int k = 0; k < N; ++k) {
                        const auto& Qp = Q[(k + N - 1) % N];
                        const auto& Qn = Q[(k + 1) % N];
                        gsPoint2D<T> avg(T(0.5) * Q[k].x + T(0.25) * (Qp.x + Qn.x),
                                         T(0.5) * Q[k].y + T(0.25) * (Qp.y + Qn.y));
                        if (pointInPolyStrict(avg, outerPoly, T(1e-5)) &&
                            isSegmentInsidePoly(P[k], avg, outerPoly) &&
                            P[k].dist(avg) > T(0.06) * clr[k])
                            Qnew[k] = avg;
                    }
                    Q = std::move(Qnew);
                }

                if (!isPolySimpleCCW(Q)) continue;
                bool spokesOk = true;
                for (int k = 0; k < N && spokesOk; ++k) {
                    if (!isSegmentInsidePoly(P[k], Q[k], outerPoly) ||
                        !isSegmentInsidePoly(Q[k], Q[(k + 1) % N], outerPoly))
                        spokesOk = false;
                }
                for (int i = 0; i < N && spokesOk; ++i) {
                    for (int j = i + 1; j < N && spokesOk; ++j) {
                        if (segProperIntersect(P[i], Q[i], P[j], Q[j])) spokesOk = false;
                    }
                    for (int j = 0; j < N && spokesOk; ++j) {
                        if (j == i || (j + 1) % N == i) continue;
                        if (segProperIntersect(P[i], Q[i], Q[j], Q[(j + 1) % N])) spokesOk = false;
                    }
                }
                if (!spokesOk) continue;

                auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
                    return makeStraightSeg(a, b, degree);
                };

                // Option 1: Pairwise ear-clipping quadification of inner even polygon Q (total = 3N/2 - 1 patches)
                std::vector<std::array<int, 4>> iQuads;
                std::vector<std::array<int, 3>> iRemTris;
                T pairMinA = 0;
                bool pairedClean = triangulateAndPairSimplePolygon(Q, minEff, maxEff, iQuads, iRemTris, pairMinA) && iRemTris.empty();
                if (pairedClean) {
                    std::vector<std::unique_ptr<gsGeometry<T>>> patches;
                    patches.reserve(N + iQuads.size());
                    for (int k = 0; k < N; ++k) {
                        int kn = (k + 1) % N;
                        CurveSeg<T> spL = mkLine(Q[k], P[k]);
                        CurveSeg<T> spR = mkLine(Q[kn], P[kn]);
                        CurveSeg<T> iEd = mkLine(Q[k], Q[kn]);
                        patches.emplace_back(makeCoonsPatchFrom4Curves(
                            *iEd.bsp, *spR.bsp, *bEdges[k].bsp, *spL.bsp, degree));
                    }
                    for (const auto& q : iQuads) {
                        patches.emplace_back(makeCoonsPatchFrom4Curves(
                            *mkLine(Q[q[0]], Q[q[1]]).bsp,
                            *mkLine(Q[q[1]], Q[q[2]]).bsp,
                            *mkLine(Q[q[3]], Q[q[2]]).bsp,
                            *mkLine(Q[q[0]], Q[q[3]]).bsp,
                            degree));
                    }
                    T totEn = 0;
                    if (verifyPatchSetLocal(patches, {}, &outerPoly, minEff, maxEff, totEn)) {
                        T cost = T(2) * static_cast<T>(patches.size()) + T(0.01) * totEn;
                        if (!foundBest || cost < bestCost) {
                            foundBest = true; bestCost = cost;
                            bestMp.clear();
                            for (auto& gp : patches) bestMp.addPatch(typename gsGeometry<T>::uPtr(gp.release()));
                        }
                    }
                }

                // Option 2: Star-kites inside Q from centroid C_Q (N/2 interior quads, total = 3N/2 patches)
                bool starOk = false;
                gsPoint2D<T> ctr = polygonCentroid(Q);
                if (pointInPolyStrict(ctr, outerPoly, T(1e-5)) && pointInPolyStrict(ctr, Q, T(1e-5))) {
                    bool allVis = true;
                    for (int k = 0; k < N && allVis; ++k)
                        if (!isSegmentInsidePoly(ctr, Q[k], Q) && !isSegmentInsidePoly(ctr, Q[k], outerPoly))
                            allVis = false;
                    if (allVis) {
                        for (int par = 0; par < 2; ++par) {
                            std::vector<std::unique_ptr<gsGeometry<T>>> patches;
                            patches.reserve(N + N / 2);
                            for (int k = 0; k < N; ++k) {
                                int kn = (k + 1) % N;
                                CurveSeg<T> spL = mkLine(Q[k], P[k]);
                                CurveSeg<T> spR = mkLine(Q[kn], P[kn]);
                                CurveSeg<T> iEd = mkLine(Q[k], Q[kn]);
                                patches.emplace_back(makeCoonsPatchFrom4Curves(
                                    *iEd.bsp, *spR.bsp, *bEdges[k].bsp, *spL.bsp, degree));
                            }
                            for (int s = 0; s < N / 2; ++s) {
                                int a = (par + 2 * s) % N;
                                int b = (a + 1) % N;
                                int c = (a + 2) % N;
                                patches.emplace_back(makeCoonsPatchFrom4Curves(
                                    *mkLine(ctr, Q[a]).bsp,
                                    *mkLine(Q[a], Q[b]).bsp,
                                    *mkLine(Q[c], Q[b]).bsp,
                                    *mkLine(ctr, Q[c]).bsp,
                                    degree));
                            }
                            T totEn = 0;
                            if (verifyPatchSetLocal(patches, {}, &outerPoly, minEff, maxEff, totEn)) {
                                starOk = true;
                                T cost = T(2) * static_cast<T>(patches.size()) + T(0.01) * totEn;
                                if (!foundBest || cost < bestCost) {
                                    foundBest = true; bestCost = cost;
                                    bestMp.clear();
                                    for (auto& gp : patches) bestMp.addPatch(typename gsGeometry<T>::uPtr(gp.release()));
                                }
                            }
                        }
                    }
                }

                // Option 3: Recursive chord meshing on inner polygon Q (keeping Q edges unsplit) when ear-pairing/star failed
                if (!pairedClean && !starOk) {
                    std::vector<CurveSeg<T>> qSegs;
                    std::vector<bool> qIsCut(N, true);
                    qSegs.reserve(N);
                    for (int k = 0; k < N; ++k) qSegs.push_back(mkLine(Q[k], Q[(k + 1) % N]));
                    int qNodes = 0;
                    std::vector<std::unique_ptr<gsGeometry<T>>> qInner;
                    if (meshHoleFreePolygon(qSegs, qIsCut, degree, minEff, maxEff, totalLen, 0, qNodes, 220, qInner)) {
                        std::vector<std::unique_ptr<gsGeometry<T>>> patches;
                        patches.reserve(N + qInner.size());
                        for (int k = 0; k < N; ++k) {
                            int kn = (k + 1) % N;
                            CurveSeg<T> spL = mkLine(Q[k], P[k]);
                            CurveSeg<T> spR = mkLine(Q[kn], P[kn]);
                            patches.emplace_back(makeCoonsPatchFrom4Curves(
                                *qSegs[k].bsp, *spR.bsp, *bEdges[k].bsp, *spL.bsp, degree));
                        }
                        for (auto& qp : qInner) patches.push_back(std::move(qp));
                        T totEn = 0;
                        if (verifyPatchSetLocal(patches, {}, &outerPoly, minEff, maxEff, totEn)) {
                            T cost = T(2) * static_cast<T>(patches.size()) + T(0.01) * totEn;
                            if (!foundBest || cost < bestCost) {
                                foundBest = true; bestCost = cost;
                                bestMp.clear();
                                for (auto& gp : patches) bestMp.addPatch(typename gsGeometry<T>::uPtr(gp.release()));
                            }
                        }
                    }
                }
            }
        }
    }

    if (foundBest) {
        bestMp.computeTopology(T(1e-4), false, false);
        outMp = give(bestMp);
        return true;
    }
    return false;
}

// General K-sector hole ring inside a K >= 3 CCW cell, supporting multi-spoke tip trisection (>285 deg)
// and matching arbitrary (K, M) hole/cell edge counts by splitting only !isCutEdge outer segments or hole segments.
template <class T>
inline bool buildGeneralKSectorHoleRing(
    const std::vector<CurveSeg<T>>& cellCCW,
    const gsCurveLoop<T>& hole,
    const std::vector<gsPoint2D<T>>& cellPoly,
    std::vector<std::unique_ptr<gsGeometry<T>>>& outPatches,
    int degree,
    T minAngTol,
    T maxAngTol,
    const std::vector<bool>& cellIsCut = {})
{
    const int K = static_cast<int>(cellCCW.size());
    if (K < 3) return false;
    std::vector<bool> isCut(K, false);
    if (cellIsCut.size() == static_cast<size_t>(K)) isCut = cellIsCut;

    std::vector<CurveSeg<T>> hRaw = orientLoopSegments(hole, degree, true);
    if (hRaw.empty()) return false;
    std::vector<CurveSeg<T>> hCCW;
    for (const auto& s : hRaw) {
        T turn = curveIntegratedTurnRad(*s.bsp, 24);
        if (turn > T(1.95)) {
            int nSub = std::min(12, std::max(2, static_cast<int>(std::ceil(turn / T(1.30)))));
            CurveSeg<T> rem = s;
            for (int k = 0; k < nSub - 1; ++k) {
                auto sp = rem.splitArc(T(1) / static_cast<T>(nSub - k));
                hCCW.push_back(sp.first);
                rem = sp.second;
            }
            hCCW.push_back(rem);
        } else {
            hCCW.push_back(s);
        }
    }
    std::vector<gsPoint2D<T>> holePoly = sampleCurveLoopPoly(hole, 16);
    if (holePoly.size() < 4 || hCCW.empty()) return false;
    {
        T a2 = 0;
        for (size_t i = 0; i < holePoly.size(); ++i) {
            const auto& p = holePoly[i];
            const auto& q = holePoly[(i + 1) % holePoly.size()];
            a2 += p.x * q.y - q.x * p.y;
        }
        if (a2 < 0) std::reverse(holePoly.begin(), holePoly.end());
    }
    const std::vector<std::vector<gsPoint2D<T>>> hVec = {holePoly};

    bool triedK4Collar = false;
    if (K == 4) {
        triedK4Collar = true;
        gsMultiPatch<T> ringMp;
        buildBlockQuadRing(*cellCCW[0].bsp, *cellCCW[1].bsp, *cellCCW[2].reversed().bsp, *cellCCW[3].reversed().bsp,
                           hole, ringMp, degree);
        if (ringMp.nPatches() >= 4) {
            bool cutPreserved = true;
            for (int s = 0; s < 4 && cutPreserved; ++s) {
                if (!isCut[s]) continue;
                gsPoint2D<T> eA = cellCCW[s].p0(), eB = cellCCW[s].p1();
                bool foundSide = false;
                for (size_t p = 0; p < ringMp.nPatches() && !foundSide; ++p) {
                    const gsMatrix<T> par0 = (gsMatrix<T>(2, 4) << 0, 1, 1, 0, 0, 0, 1, 1).finished();
                    gsMatrix<T> corners;
                    ringMp.patch(p).eval_into(par0, corners);
                    for (int sd = 0; sd < 4; ++sd) {
                        int sdN = (sd + 1) % 4;
                        gsPoint2D<T> c0(corners(0, sd), corners(1, sd));
                        gsPoint2D<T> c1(corners(0, sdN), corners(1, sdN));
                        if ((c0.dist(eA) < T(1e-5) && c1.dist(eB) < T(1e-5)) ||
                            (c0.dist(eB) < T(1e-5) && c1.dist(eA) < T(1e-5))) {
                            foundSide = true;
                            break;
                        }
                    }
                }
                if (!foundSide) cutPreserved = false;
            }
            if (cutPreserved) {
                std::vector<std::unique_ptr<gsGeometry<T>>> cands;
                cands.reserve(ringMp.nPatches());
                for (size_t p = 0; p < ringMp.nPatches(); ++p)
                    cands.emplace_back(ringMp.patch(p).clone().release());
                T en = 0;
                if (verifyPatchSetLocal(cands, hVec, &cellPoly, minAngTol, maxAngTol, en)) {
                    outPatches = std::move(cands);
                    return true;
                }
            }
        }
    }

    const T radToDeg = T(180) / T(3.14159265358979323846);
    std::vector<gsPoint2D<T>> P(K);
    for (int s = 0; s < K; ++s) P[s] = cellCCW[s].p0();

    auto solveRingForInnerLoop = [&](const std::vector<CurveSeg<T>>& inLoopCCW,
                                     const std::vector<gsPoint2D<T>>& inLoopPoly,
                                     std::vector<std::unique_ptr<gsGeometry<T>>>& ringOut) -> bool {
        const int M = static_cast<int>(inLoopCCW.size());
        const std::vector<std::vector<gsPoint2D<T>>> inVec = {inLoopPoly};
        std::vector<T> hTheta(M, T(180));
        int nSharpH = 0;
        for (int i = 0; i < M; ++i) {
            auto d0 = evalCurveDeriv(*inLoopCCW[(i + M - 1) % M].bsp, T(1));
            auto d1 = evalCurveDeriv(*inLoopCCW[i].bsp, T(0));
            T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
            if (n0 > T(1e-14) && n1 > T(1e-14)) {
                T cr = (d0.x * d1.y - d0.y * d1.x) / (n0 * n1);
                T dt = (d0.x * d1.x + d0.y * d1.y) / (n0 * n1);
                T turnDeg = std::atan2(cr, dt) * radToDeg;
                if (std::abs(turnDeg) > T(38)) ++nSharpH;
                T extAng = T(180) + turnDeg;
                if (extAng <= T(0)) extAng += T(360);
                if (extAng > T(360)) extAng -= T(360);
                hTheta[i] = extAng;
            }
        }

        auto trySpokeLayout = [&](const std::vector<int>& mBase) -> bool {
            std::vector<CurveSeg<T>> cellExp = cellCCW;
            std::vector<bool> cutExp = isCut;
            std::vector<CurveSeg<T>> hExp = inLoopCCW;
            std::vector<int> mExp = mBase;

            auto neededK = [&]() -> int {
                int sumS = 0;
                for (int mv : mExp) sumS += mv;
                return 2 * sumS - static_cast<int>(hExp.size());
            };

            int kNeed = neededK();
            while (static_cast<int>(cellExp.size()) < kNeed) {
                int bestIdx = -1;
                T bestLen = -T(1);
                for (int i = 0; i < static_cast<int>(cellExp.size()); ++i) {
                    if (cutExp[i]) continue;
                    T len = cellExp[i].p0().dist(cellExp[i].p1());
                    if (len <= T(1e-14))
                        len = T(2) * cellExp[i].p0().dist(evalCurvePoint(*cellExp[i].bsp, T(0.5)));
                    if (len > bestLen) { bestLen = len; bestIdx = i; }
                }
                if (bestIdx < 0) return false;
                T tMid = arcLengthSplitParam(*cellExp[bestIdx].bsp, T(0.5));
                auto sp = cellExp[bestIdx].splitAtParam(tMid);
                cellExp[bestIdx] = sp.first;
                cutExp[bestIdx] = false;
                cellExp.insert(cellExp.begin() + bestIdx + 1, sp.second);
                cutExp.insert(cutExp.begin() + bestIdx + 1, false);
            }

            while (kNeed < static_cast<int>(cellExp.size())) {
                int bestIdx = -1;
                T bestLen = -T(1);
                for (int i = 0; i < static_cast<int>(hExp.size()); ++i) {
                    T len = hExp[i].p0().dist(hExp[i].p1());
                    if (len <= T(1e-14))
                        len = T(2) * hExp[i].p0().dist(evalCurvePoint(*hExp[i].bsp, T(0.5)));
                    if (len > bestLen) { bestLen = len; bestIdx = i; }
                }
                if (bestIdx < 0) return false;
                T tMid = arcLengthSplitParam(*hExp[bestIdx].bsp, T(0.5));
                auto sp = hExp[bestIdx].splitAtParam(tMid);
                hExp[bestIdx] = sp.first;
                hExp.insert(hExp.begin() + bestIdx + 1, sp.second);
                mExp.insert(mExp.begin() + bestIdx + 1, 1);
                kNeed = neededK();
            }

            const int Ke = static_cast<int>(cellExp.size());
            const int Me = static_cast<int>(hExp.size());
            if (Ke != kNeed || Ke < 3 || Me < 2) return false;

            std::vector<gsPoint2D<T>> Pe(Ke);
            std::vector<T> cAng(Ke, T(90));
            for (int s = 0; s < Ke; ++s) {
                Pe[s] = cellExp[s].p0();
                auto din = evalCurveDeriv(*cellExp[(s + Ke - 1) % Ke].bsp, T(1));
                auto dout = evalCurveDeriv(*cellExp[s].bsp, T(0));
                T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
                if (n1 > T(1e-14) && n2 > T(1e-14)) {
                    T cr = (din.x * dout.y - din.y * dout.x) / (n1 * n2);
                    T dt = (din.x * dout.x + din.y * dout.y) / (n1 * n2);
                    T ang = T(180) - std::atan2(cr, dt) * radToDeg;
                    if (ang <= T(0)) ang += T(360);
                    cAng[s] = ang;
                }
            }

            std::vector<int> rotOrder(Ke);
            std::iota(rotOrder.begin(), rotOrder.end(), 0);
            if (Ke > 10) {
                gsPoint2D<T> h0 = hExp[0].p0();
                std::sort(rotOrder.begin(), rotOrder.end(), [&](int ra, int rb) {
                    return Pe[ra].dist(h0) < Pe[rb].dist(h0);
                });
                rotOrder.resize(8);
            }

            T bestEn = T(1e30);
            std::vector<std::unique_ptr<gsGeometry<T>>> bestP;
            bool okAny = false;

            for (int rot : rotOrder) {
                std::vector<int> vStart(Me), vEnd(Me);
                int curA = rot;
                bool wedgeAngOk = true;
                for (int j = 0; j < Me; ++j) {
                    vStart[j] = curA;
                    if (mExp[j] == 2) {
                        int midV = (curA + 1) % Ke;
                        if (cAng[midV] > std::min(maxAngTol, T(162))) {
                            wedgeAngOk = false;
                            break;
                        }
                        curA = (curA + 2) % Ke;
                        vEnd[j] = curA;
                        curA = (curA + 1) % Ke;
                    } else {
                        vEnd[j] = curA;
                        curA = (curA + 1) % Ke;
                    }
                }
                if (!wedgeAngOk) continue;

                // Curved-spoke builder: falls back to bowed cubic Bezier spoke when straight line-of-sight
                // clips a curved hole or concave cell wall.
                auto buildClearSpoke = [&](const CurveSeg<T>& hp, const CurveSeg<T>& hn,
                                           const gsPoint2D<T>& H, const gsPoint2D<T>& Pdst,
                                           CurveSeg<T>& outSeg, std::vector<gsPoint2D<T>>& outPoly) -> bool {
                    auto sample16 = [](const CurveSeg<T>& seg) -> std::vector<gsPoint2D<T>> {
                        std::vector<gsPoint2D<T>> pts(17);
                        for (int k = 0; k <= 16; ++k)
                            pts[k] = evalCurvePoint(*seg.bsp, static_cast<T>(k) / T(16));
                        return pts;
                    };
                    if (isHoleSpokeClearOfHole(H, Pdst, inLoopPoly) && isSegmentInsidePoly(H, Pdst, cellPoly)) {
                        bool wc = false;
                        outSeg = makeExteriorConeSafeSpoke(hp, hn, H, Pdst, inLoopPoly, degree, wc);
                        outPoly = sample16(outSeg);
                        bool ok = true;
                        for (int k = 1; k < 16; ++k)
                            if (pointInPolyStrict(outPoly[k], inLoopPoly, T(1e-6)) || !pointInPolyStrict(outPoly[k], cellPoly, T(1e-6)))
                                { ok = false; break; }
                        if (ok) return true;
                    }
                    // Try curved cubic Bezier spokes bowing around concave hole/outer walls
                    T L = H.dist(Pdst);
                    if (L <= T(1e-10)) return false;
                    gsPoint2D<T> d((Pdst.x - H.x) / L, (Pdst.y - H.y) / L);
                    gsPoint2D<T> perp(-d.y, d.x);
                    for (T bow : {T(0.20), T(-0.20), T(0.35), T(-0.35), T(0.50), T(-0.50)}) {
                        gsPoint2D<T> c1(H.x + (L / T(3)) * d.x + bow * L * perp.x,
                                        H.y + (L / T(3)) * d.y + bow * L * perp.y);
                        gsPoint2D<T> c2(H.x + (T(2) * L / T(3)) * d.x + bow * L * perp.x,
                                        H.y + (T(2) * L / T(3)) * d.y + bow * L * perp.y);
                        gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
                        gsMatrix<T> cf(4, 2);
                        cf << H.x, H.y, c1.x, c1.y, c2.x, c2.y, Pdst.x, Pdst.y;
                        gsBSpline<T> bsp(kv, cf);
                        if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
                        CurveSeg<T> cand(bsp);
                        std::vector<gsPoint2D<T>> poly = sample16(cand);
                        bool ok = true;
                        for (int k = 1; k < 16; ++k) {
                            if (pointInPolyStrict(poly[k], inLoopPoly, T(1e-6)) ||
                                !pointInPolyStrict(poly[k], cellPoly, T(1e-6))) {
                                ok = false; break;
                            }
                        }
                        if (ok) {
                            // Verify sub-segments don't cross hole or cell boundary
                            for (int k = 1; k + 2 < 16 && ok; ++k) {
                                for (size_t e = 0; e < inLoopPoly.size() && ok; ++e)
                                    if (segProperIntersect(poly[k], poly[k + 1], inLoopPoly[e], inLoopPoly[(e + 1) % inLoopPoly.size()]))
                                        ok = false;
                                for (size_t e = 0; e < cellPoly.size() && ok; ++e)
                                    if (segProperIntersect(poly[k], poly[k + 1], cellPoly[e], cellPoly[(e + 1) % cellPoly.size()]))
                                        ok = false;
                            }
                        }
                        if (ok) {
                            outSeg = cand;
                            outPoly = std::move(poly);
                            return true;
                        }
                    }
                    return false;
                };

                auto spokesCross = [](const std::vector<gsPoint2D<T>>& pa, const std::vector<gsPoint2D<T>>& pb) -> bool {
                    for (size_t i = 0; i + 1 < pa.size(); ++i) {
                        for (size_t j = 0; j + 1 < pb.size(); ++j) {
                            if ((i == 0 && j == 0 && pa[0].dist(pb[0]) < T(1e-6)) ||
                                (i + 2 == pa.size() && j + 2 == pb.size() && pa.back().dist(pb.back()) < T(1e-6)))
                                continue;
                            if (segProperIntersect(pa[i], pa[i + 1], pb[j], pb[j + 1]))
                                return true;
                        }
                    }
                    return false;
                };

                std::vector<CurveSeg<T>> spA(Me), spB(Me);
                bool visOk = true;
                std::vector<std::vector<gsPoint2D<T>>> spPolys;
                spPolys.reserve(2 * Me);
                for (int j = 0; j < Me && visOk; ++j) {
                    gsPoint2D<T> Hj = hExp[j].p0();
                    const auto& hp = hExp[(j + Me - 1) % Me];
                    const auto& hn = hExp[j];
                    std::vector<gsPoint2D<T>> polyA;
                    if (!buildClearSpoke(hp, hn, Hj, Pe[vStart[j]], spA[j], polyA)) {
                        visOk = false; break;
                    }
                    spPolys.push_back(std::move(polyA));
                    if (mExp[j] == 2) {
                        std::vector<gsPoint2D<T>> polyB;
                        if (!buildClearSpoke(hp, hn, Hj, Pe[vEnd[j]], spB[j], polyB)) {
                            visOk = false; break;
                        }
                        spPolys.push_back(std::move(polyB));
                    } else {
                        spB[j] = spA[j];
                    }
                }
                for (size_t i = 0; i < spPolys.size() && visOk; ++i) {
                    for (size_t j = i + 1; j < spPolys.size() && visOk; ++j) {
                        if (spokesCross(spPolys[i], spPolys[j]))
                            visOk = false;
                    }
                }
                if (!visOk) continue;

                std::vector<std::unique_ptr<gsGeometry<T>>> cands;
                cands.reserve(Ke);
                for (int j = 0; j < Me; ++j) {
                    int jn = (j + 1) % Me;
                    if (mExp[j] == 2) {
                        int midE = (vStart[j] + 1) % Ke;
                        cands.emplace_back(makeCoonsFromCCWCycle(
                            cellExp[vStart[j]], cellExp[midE], spB[j].reversed(), spA[j], degree));
                    }
                    cands.emplace_back(makeCoonsPatchFrom4Curves(
                        *cellExp[vEnd[j]].bsp, *spA[jn].reversed().bsp, *hExp[j].bsp, *spB[j].reversed().bsp, degree));
                }
                T en = 0;
                if (verifyPatchSetLocal(cands, inVec, &cellPoly, minAngTol, maxAngTol, en) && (!okAny || en < bestEn)) {
                    okAny = true; bestEn = en; bestP = std::move(cands);
                }
            }
            if (okAny) {
                ringOut = std::move(bestP);
                return true;
            }
            return false;
        };

        std::vector<int> mMulti(M, 1), mSingle(M, 1);
        bool hasTip285 = false;
        for (int i = 0; i < M; ++i) {
            if (hTheta[i] > T(285)) { mMulti[i] = 2; hasTip285 = true; }
        }
        if (hasTip285 && trySpokeLayout(mMulti))
            return true;
        if (trySpokeLayout(mSingle))
            return true;

        if (nSharpH == 0) {
            auto hFullPtr = joinCurveSegsExact(inLoopCCW, degree);
            if (!hFullPtr) return false;
            CurveSeg<T> hFull(*hFullPtr);
            gsPoint2D<T> hCtr = polygonCentroid(inLoopPoly);

            const int NS = 128;
            std::vector<gsPoint2D<T>> sPts(NS);
            std::vector<T> sAng(NS);
            for (int i = 0; i < NS; ++i) {
                sPts[i] = evalCurvePoint(*hFull.bsp, static_cast<T>(i) / static_cast<T>(NS));
                sAng[i] = std::atan2(sPts[i].y - hCtr.y, sPts[i].x - hCtr.x);
            }

            const T twoPi = T(2) * T(3.14159265358979323846);
            auto angDiff = [&](T a, T b) -> T {
                T d = std::abs(a - b);
                return (d > T(3.14159265358979323846)) ? (twoPi - d) : d;
            };

            std::vector<T> uBase(K);
            for (int pickMode = 0; pickMode < 2; ++pickMode) {
                for (int s = 0; s < K; ++s) {
                    if (pickMode == 0) {
                        T bestD = T(1e30); int bi = 0;
                        for (int i = 0; i < NS; ++i) {
                            T d = sPts[i].dist(P[s]);
                            if (d < bestD && isHoleSpokeClearOfHole(sPts[i], P[s], inLoopPoly)) {
                                bestD = d; bi = i;
                            }
                        }
                        uBase[s] = static_cast<T>(bi) / static_cast<T>(NS);
                    } else {
                        T targetA = std::atan2(P[s].y - hCtr.y, P[s].x - hCtr.x);
                        T bestD = T(1e30); int bi = 0;
                        for (int i = 0; i < NS; ++i) {
                            T d = angDiff(sAng[i], targetA);
                            if (d < bestD) { bestD = d; bi = i; }
                        }
                        uBase[s] = static_cast<T>(bi) / static_cast<T>(NS);
                    }
                }

                T u0 = uBase[0];
                std::vector<T> relU(K - 1);
                bool cyclicOk = true;
                for (int s = 1; s < K; ++s) {
                    T ur = uBase[s] - u0;
                    if (ur < T(0)) ur += T(1);
                    relU[s - 1] = ur;
                }
                for (int s = 0; s < K - 1; ++s) {
                    if (relU[s] < T(0.03) || relU[s] > T(0.97)) cyclicOk = false;
                    if (s > 0 && relU[s] - relU[s - 1] < T(0.03)) cyclicOk = false;
                }
                if (!cyclicOk) continue;

                CurveSeg<T> hShifted = hFull;
                if (u0 > T(1e-4) && u0 < T(1) - T(1e-4)) {
                    auto sp0 = hFull.splitAtParam(u0);
                    auto jPtr = joinCurveSegsExact(std::vector<CurveSeg<T>>{sp0.second, sp0.first}, degree);
                    if (!jPtr) continue;
                    hShifted = CurveSeg<T>(*jPtr);
                }
                std::vector<CurveSeg<T>> hSec = splitSegMulti(hShifted, relU);
                if (static_cast<int>(hSec.size()) != K) continue;

                std::vector<CurveSeg<T>> sp(K);
                bool visOk = true;
                for (int s = 0; s < K; ++s) {
                    gsPoint2D<T> H = hSec[s].p0();
                    if (!isHoleSpokeClearOfHole(H, P[s], inLoopPoly) || !isSegmentInsidePoly(H, P[s], cellPoly)) {
                        visOk = false; break;
                    }
                    bool wc = false;
                    sp[s] = makeExteriorConeSafeSpoke(hSec[(s + K - 1) % K], hSec[s], H, P[s], inLoopPoly, degree, wc);
                }
                for (int i = 0; i < K && visOk; ++i)
                    for (int j = i + 1; j < K && visOk; ++j)
                        if (segProperIntersect(sp[i].p0(), sp[i].p1(), sp[j].p0(), sp[j].p1())) visOk = false;
                if (!visOk) continue;

                std::vector<std::unique_ptr<gsGeometry<T>>> cands;
                cands.reserve(K);
                for (int s = 0; s < K; ++s) {
                    int sn = (s + 1) % K;
                    cands.emplace_back(makeCoonsPatchFrom4Curves(
                        *cellCCW[s].bsp, *sp[sn].reversed().bsp, *hSec[s].bsp, *sp[s].reversed().bsp, degree));
                }
                T en = 0;
                if (verifyPatchSetLocal(cands, inVec, &cellPoly, minAngTol, maxAngTol, en)) {
                    ringOut = std::move(cands);
                    return true;
                }
            }
        }
        return false;
    };

    if (solveRingForInnerLoop(hCCW, holePoly, outPatches))
        return true;

    // Collar fallback for sharp spikes (>280 deg) or strongly twisted holes inside any K-cell
    std::vector<std::unique_ptr<gsGeometry<T>>> collarP;
    std::vector<CurveSeg<T>> colOuter;
    std::vector<int> edgeToPatch;
    if (buildHoleSmoothingCollar(hole, &cellPoly, degree, minAngTol, maxAngTol, collarP, colOuter, edgeToPatch)) {
        auto tryCombineCollar = [&](std::vector<std::unique_ptr<gsGeometry<T>>>& cp,
                                    const std::vector<CurveSeg<T>>& co) -> bool {
            std::vector<gsPoint2D<T>> colPoly;
            colPoly.reserve(co.size() * 8);
            for (const auto& cs : co)
                for (int k = 0; k < 8; ++k)
                    colPoly.push_back(evalCurvePoint(*cs.bsp, static_cast<T>(k) / T(8)));

            std::vector<std::unique_ptr<gsGeometry<T>>> outerRing;
            if (solveRingForInnerLoop(co, colPoly, outerRing)) {
                outPatches.clear();
                for (auto& g : cp) outPatches.push_back(std::move(g));
                for (auto& g : outerRing) outPatches.push_back(std::move(g));
                T en = 0;
                return verifyPatchSetLocal(outPatches, hVec, &cellPoly, minAngTol, maxAngTol, en);
            }
            if (K == 4 && !triedK4Collar && (co.size() % 2) == 0 && co.size() >= 4) {
                const gsPoint2D<T> P4[4] = {cellCCW[0].p0(), cellCCW[1].p0(), cellCCW[2].p0(), cellCCW[3].p0()};
                T rSc = 0, rMa = T(180);
                if (runSharpBlockRingSearch(
                        cellCCW[0], cellCCW[1], cellCCW[2].reversed(), cellCCW[3].reversed(),
                        P4, co, colPoly, degree,
                        {T(0.5)}, {T(0.62), T(0.45), T(0.78)}, true,
                        outerRing, rSc, rMa)) {
                    outPatches.clear();
                    for (auto& g : cp) outPatches.push_back(std::move(g));
                    for (auto& g : outerRing) outPatches.push_back(std::move(g));
                    T en = 0;
                    return verifyPatchSetLocal(outPatches, hVec, &cellPoly, minAngTol, maxAngTol, en);
                }
            }
            return false;
        };

        if (tryCombineCollar(collarP, colOuter))
            return true;
        if ((colOuter.size() % 2) == 1 && splitOneCollarFlankForEvenness(collarP, colOuter, edgeToPatch, degree)) {
            if (tryCombineCollar(collarP, colOuter))
                return true;
        }
    }

    return false;
}

// Generalized rotation-invariant recursive cell splitter for arbitrary 0-hole and multi-hole domains.
// Guarantees 1-to-1 conforming interfaces by never splitting interior cut edges (`isCutEdge == true`).
template <class T>
inline bool computeRecursiveCellSplitting(
    const gsPlanarDomain<T>& dom,
    int degree,
    T weightPatches,
    T weightDistortion,
    T minAngTol,
    T maxAngTol,
    gsMultiPatch<T>& outMp,
    std::string& outLabel,
    std::function<bool(const std::vector<CurveSeg<T>>&,
                       const std::vector<bool>&,
                       const gsCurveLoop<T>*,
                       std::vector<std::unique_ptr<gsGeometry<T>>>&)> leafFallback = nullptr)
{
    if (dom.numLoops() == 0 || dom.outer().size() == 0) return false;
    const gsMatrix<T> bbox = dom.boundingBox();
    const T domDiag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    if (domDiag <= T(1e-12)) return false;

    const int nHoles = dom.numHoles();
    std::vector<std::vector<gsPoint2D<T>>> holePolys(nHoles);
    std::vector<gsPoint2D<T>> holeCtrs(nHoles);
    std::vector<T> holeRadii(nHoles, T(0));
    for (int h = 0; h < nHoles; ++h) {
        const auto& hl = dom.loop(h + 1);
        for (int c = 0; c < hl.size(); ++c) {
            gsMatrix<T> pr = hl.curve(c).parameterRange();
            for (int s = 0; s < 12; ++s) {
                T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(s) / T(12));
                holePolys[h].push_back(evalCurvePoint(hl.curve(c), t));
            }
        }
        holeCtrs[h] = polygonCentroid(holePolys[h]);
        for (const auto& hp : holePolys[h])
            holeRadii[h] = std::max(holeRadii[h], hp.dist(holeCtrs[h]));
    }

    struct HalfSeg {
        CurveSeg<T> seg;
        bool isCutEdge = false;
    };
    struct Cell {
        std::vector<HalfSeg> boundary;
        std::vector<int> holeIds;
    };

    const T radToDeg = T(180) / T(3.14159265358979323846);
    auto sampleCellPoly = [](const std::vector<HalfSeg>& bnd, int nPerSeg = 12) -> std::vector<gsPoint2D<T>> {
        std::vector<gsPoint2D<T>> poly;
        poly.reserve(bnd.size() * nPerSeg);
        for (const auto& hs : bnd) {
            for (int k = 0; k < nPerSeg; ++k)
                poly.push_back(evalCurvePoint(*hs.seg.bsp, static_cast<T>(k) / static_cast<T>(nPerSeg)));
        }
        return poly;
    };

    auto evalCellVertexAngleDeg = [&](const CurveSeg<T>& prevSeg, const CurveSeg<T>& nextSeg) -> T {
        auto dIn = evalCurveDeriv(*prevSeg.bsp, T(1));
        auto dOut = evalCurveDeriv(*nextSeg.bsp, T(0));
        T n1 = std::hypot(dIn.x, dIn.y), n2 = std::hypot(dOut.x, dOut.y);
        if (n1 <= T(1e-14) || n2 <= T(1e-14)) return T(180);
        T cr = (dIn.x * dOut.y - dIn.y * dOut.x) / (n1 * n2);
        T dt = (dIn.x * dOut.x + dIn.y * dOut.y) / (n1 * n2);
        T ang = T(180) - std::atan2(cr, dt) * radToDeg;
        if (ang <= T(0)) ang += T(360);
        return ang;
    };

    auto evalCornerBetweenCuts = [&](const gsPoint2D<T>& tInFwd, const gsPoint2D<T>& tOutFwd) -> T {
        T cr = tInFwd.x * tOutFwd.y - tInFwd.y * tOutFwd.x;
        T dt = tInFwd.x * tOutFwd.x + tInFwd.y * tOutFwd.y;
        T ang = T(180) - std::atan2(cr, dt) * radToDeg;
        if (ang <= T(0)) ang += T(360);
        return ang;
    };

    auto unitTanAt = [](const CurveSeg<T>& seg, T u) -> gsPoint2D<T> {
        gsPoint2D<T> d = evalCurveDeriv(*seg.bsp, u);
        T n = std::hypot(d.x, d.y);
        return (n > T(1e-14)) ? gsPoint2D<T>(d.x / n, d.y / n) : gsPoint2D<T>(1, 0);
    };

    auto segMinDistPoly = [](const std::vector<gsPoint2D<T>>& segPts, const std::vector<gsPoint2D<T>>& poly, bool checkCross) -> T {
        const int np = static_cast<int>(poly.size());
        const int ns = static_cast<int>(segPts.size());
        if (checkCross) {
            for (int a = 0; a + 1 < ns; ++a) {
                for (int i = 0; i < np; ++i) {
                    if (segProperIntersect(segPts[a], segPts[a + 1], poly[i], poly[(i + 1) % np]))
                        return T(-1);
                }
            }
        }
        T minD2 = T(1e30);
        for (int a = 0; a < ns; ++a) {
            const auto& p = segPts[a];
            if (pointInPolyStrict(p, poly, T(1e-7))) return T(-1);
            for (int i = 0; i < np; ++i) {
                const auto& A = poly[i];
                const auto& B = poly[(i + 1) % np];
                T dx = B.x - A.x, dy = B.y - A.y;
                T L2 = dx * dx + dy * dy;
                T t = (L2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((p.x - A.x) * dx + (p.y - A.y) * dy) / L2)) : T(0);
                T qx = A.x + t * dx - p.x, qy = A.y + t * dy - p.y;
                minD2 = std::min(minD2, qx * qx + qy * qy);
            }
        }
        return std::sqrt(minD2);
    };

    auto sampleCurvePts = [](const CurveSeg<T>& seg, int M = 16) -> std::vector<gsPoint2D<T>> {
        std::vector<gsPoint2D<T>> res(M + 1);
        gsMatrix<T> u(1, M + 1), val;
        for (int i = 0; i <= M; ++i) u(0, i) = static_cast<T>(i) / static_cast<T>(M);
        seg.bsp->eval_into(u, val);
        for (int i = 0; i <= M; ++i) res[i] = gsPoint2D<T>(val(0, i), val(1, i));
        return res;
    };

    auto makeCurvedBezierSeg = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& Mctrl, const gsPoint2D<T>& B) -> CurveSeg<T> {
        gsKnotVector<T> kv(T(0), T(1), 0, 3, 1, 2);
        gsMatrix<T> cf(3, 2);
        cf << A.x, A.y, Mctrl.x, Mctrl.y, B.x, B.y;
        gsBSpline<T> bsp(kv, cf);
        if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
        return CurveSeg<T>(bsp);
    };

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        return makeStraightSeg(a, b, degree);
    };

    // Build initial root Cell (keeping corners intact, subdividing outer curves whose turning > 2.15 rad)
    Cell rootCell;
    {
        std::vector<CurveSeg<T>> oCCW = orientLoopSegments(dom.outer(), degree, true);
        if (oCCW.empty()) return false;
        for (const auto& cs : oCCW) {
            T turn = curveIntegratedTurnRad(*cs.bsp, 24);
            if (turn > T(2.15)) {
                int nSub = std::min(8, std::max(2, static_cast<int>(std::ceil(turn / T(1.35)))));
                CurveSeg<T> rem = cs;
                for (int k = 0; k < nSub - 1; ++k) {
                    auto sp = rem.splitArc(T(1) / static_cast<T>(nSub - k));
                    rootCell.boundary.push_back({sp.first, false});
                    rem = sp.second;
                }
                rootCell.boundary.push_back({rem, false});
            } else {
                rootCell.boundary.push_back({cs, false});
            }
        }
        for (int h = 0; h < nHoles; ++h)
            rootCell.holeIds.push_back(h);
    }

    int solverNodes = 0;
    int maxSolverNodes = 450;
    int leafCalls = 0;
    const int maxLeafCalls = 24;
    const auto tStart = std::chrono::steady_clock::now();
    auto timeExpired = [&]() -> bool {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count() > 3.5;
    };

    // Recursive hole-free quad mesher delegating to standalone meshHoleFreePolygon
    auto meshHoleFree = [&](const Cell& cell, int depth, std::vector<std::unique_ptr<gsGeometry<T>>>& outP) -> bool {
        if (timeExpired()) return false;
        std::vector<CurveSeg<T>> segs;
        std::vector<bool> cuts;
        segs.reserve(cell.boundary.size());
        cuts.reserve(cell.boundary.size());
        for (const auto& hs : cell.boundary) {
            segs.push_back(hs.seg);
            cuts.push_back(hs.isCutEdge);
        }
        return meshHoleFreePolygon(segs, cuts, degree, minAngTol, maxAngTol, domDiag, depth, solverNodes, maxSolverNodes, outP);
    };

    // Main recursive cell solver for arbitrary hole count >= 0
    std::function<bool(const Cell&, int, std::vector<std::unique_ptr<gsGeometry<T>>>&)> solveCell;
    solveCell = [&](const Cell& cell, int depth, std::vector<std::unique_ptr<gsGeometry<T>>>& outP) -> bool {
        if (depth > 14 || ++solverNodes > maxSolverNodes || timeExpired()) return false;
        outP.clear();

        if (cell.holeIds.empty()) {
            if (meshHoleFree(cell, depth, outP))
                return true;
            if (depth > 0 && leafFallback) {
                std::vector<CurveSeg<T>> leafCCW;
                std::vector<bool> leafIsCut;
                for (const auto& hs : cell.boundary) {
                    leafCCW.push_back(hs.seg);
                    leafIsCut.push_back(hs.isCutEdge);
                }
                std::vector<std::unique_ptr<gsGeometry<T>>> fbP;
                if (leafFallback(leafCCW, leafIsCut, nullptr, fbP) && !fbP.empty()) {
                    outP = std::move(fbP);
                    return true;
                }
            }
            return false;
        }

        const int K = static_cast<int>(cell.boundary.size());
        const std::vector<gsPoint2D<T>> cPoly = sampleCellPoly(cell.boundary, 12);

        // Compute corner angles & count reflex corners
        std::vector<T> angDeg(K);
        std::vector<int> reflexIdx;
        for (int i = 0; i < K; ++i) {
            angDeg[i] = evalCellVertexAngleDeg(cell.boundary[(i + K - 1) % K].seg, cell.boundary[i].seg);
            if (angDeg[i] > T(194)) reflexIdx.push_back(i);
        }

        // Direct 1-hole leaf solve when cell has 1 hole (and is mostly convex, <=1 reflex, or K <= 5)
        if (cell.holeIds.size() == 1 && (reflexIdx.size() <= 1 || K <= 5) && leafCalls < maxLeafCalls) {
            int hid = cell.holeIds[0];
            std::vector<HalfSeg> compBnd = cell.boundary;
            while (compBnd.size() > 4) {
                const int KC = static_cast<int>(compBnd.size());
                int mIdx = -1;
                for (int i = 0; i < KC; ++i) {
                    int ip = (i + KC - 1) % KC;
                    if (!compBnd[ip].isCutEdge && !compBnd[i].isCutEdge) {
                        T ang = evalCellVertexAngleDeg(compBnd[ip].seg, compBnd[i].seg);
                        if (std::abs(ang - T(180)) <= T(12) &&
                            curveIntegratedTurnRad(*compBnd[ip].seg.bsp) + curveIntegratedTurnRad(*compBnd[i].seg.bsp) <= T(1.85)) {
                            mIdx = i; break;
                        }
                    }
                }
                if (mIdx < 0) break;
                int ip = (mIdx + KC - 1) % KC;
                auto jPtr = joinCurveSegsExact(std::vector<CurveSeg<T>>{compBnd[ip].seg, compBnd[mIdx].seg}, degree);
                if (!jPtr) break;
                compBnd[ip].seg = CurveSeg<T>(*jPtr);
                compBnd.erase(compBnd.begin() + mIdx);
            }
            const auto& useBnd = (static_cast<int>(compBnd.size()) < K) ? compBnd : cell.boundary;
            std::vector<CurveSeg<T>> leafCCW;
            std::vector<bool> leafIsCut;
            leafCCW.reserve(useBnd.size());
            leafIsCut.reserve(useBnd.size());
            for (const auto& hs : useBnd) {
                leafCCW.push_back(hs.seg);
                leafIsCut.push_back(hs.isCutEdge);
            }
            const std::vector<gsPoint2D<T>> leafPoly =
                (static_cast<int>(compBnd.size()) < K) ? sampleCellPoly(compBnd, 12) : cPoly;
            ++leafCalls;
            solverNodes += 12;
            std::vector<std::unique_ptr<gsGeometry<T>>> leafP;
            if (buildGeneralKSectorHoleRing(leafCCW, dom.loop(hid + 1), leafPoly, leafP, degree, minAngTol, maxAngTol, leafIsCut)) {
                outP = std::move(leafP);
                return true;
            }
            if (depth > 0 && leafFallback && !timeExpired()) {
                std::vector<std::unique_ptr<gsGeometry<T>>> fbP;
                if (leafFallback(leafCCW, leafIsCut, &dom.loop(hid + 1), fbP) && !fbP.empty()) {
                    outP = std::move(fbP);
                    return true;
                }
            }
            if (timeExpired() || solverNodes > maxSolverNodes) return false;
        }

        // Pre-sample all K boundary segments once (16 intervals each) for fast ray-casting & polygon slicing
        std::vector<std::vector<gsPoint2D<T>>> segPts(K);
        for (int i = 0; i < K; ++i)
            segPts[i] = sampleCurvePts(cell.boundary[i].seg, 16);

        // Generate candidate boundary anchors on Cell (vertices + interior points on !isCutEdge sides)
        struct Anchor {
            int edgeIdx;
            T u;
            bool isVertex;
            gsPoint2D<T> pt;
            gsPoint2D<T> tIn;  // arriving at pt along boundary
            gsPoint2D<T> tOut; // leaving pt along boundary
        };
        std::vector<Anchor> anchors;
        for (int i = 0; i < K; ++i) {
            anchors.push_back({
                i, T(0), true,
                cell.boundary[i].seg.p0(),
                unitTanAt(cell.boundary[(i + K - 1) % K].seg, T(1)),
                unitTanAt(cell.boundary[i].seg, T(0))
            });
        }

        const size_t maxAnchors = 44;
        auto addEdgeAnchor = [&](int ei, T u) {
            if (anchors.size() >= maxAnchors) return;
            if (cell.boundary[ei].isCutEdge) return;
            if (u < T(0.16) || u > T(0.84)) return;
            gsPoint2D<T> p = evalCurvePoint(*cell.boundary[ei].seg.bsp, u);
            for (const auto& ex : anchors) {
                if (!ex.isVertex && ex.edgeIdx == ei && std::abs(ex.u - u) < T(0.08)) return;
                if (ex.pt.dist(p) < T(0.04) * domDiag) return;
            }
            gsPoint2D<T> t = unitTanAt(cell.boundary[ei].seg, u);
            anchors.push_back({ei, u, false, p, t, t});
        };

        // 1) Ray intersections along reflex bisectors and hole-pair bisectors onto !isCutEdge boundaries (using cached segPts)
        auto projectRayToBoundary = [&](const gsPoint2D<T>& orig, const gsPoint2D<T>& dir) {
            T bestT = T(1e18);
            int bestE = -1;
            T bestU = T(0.5);
            for (int i = 0; i < K; ++i) {
                if (cell.boundary[i].isCutEdge) continue;
                const auto& spts = segPts[i];
                for (int s = 0; s < 16; ++s) {
                    const auto& A = spts[s];
                    const auto& B = spts[s + 1];
                    T sx = B.x - A.x, sy = B.y - A.y;
                    T ax = A.x - orig.x, ay = A.y - orig.y;
                    T denom = dir.x * sy - dir.y * sx;
                    if (std::abs(denom) <= T(1e-12)) continue;
                    T t = (ax * sy - ay * sx) / denom;
                    T w = (ax * dir.y - ay * dir.x) / denom;
                    if (t > T(0.03) * domDiag && t < bestT && w >= T(0) && w <= T(1)) {
                        bestT = t;
                        bestE = i;
                        bestU = (static_cast<T>(s) + w) / T(16);
                    }
                }
            }
            if (bestE >= 0) addEdgeAnchor(bestE, bestU);
        };

        for (int r : reflexIdx) {
            gsPoint2D<T> tPrev = unitTanAt(cell.boundary[(r + K - 1) % K].seg, T(1));
            gsPoint2D<T> tNext = unitTanAt(cell.boundary[r].seg, T(0));
            gsPoint2D<T> bis(tPrev.x - tNext.x, tPrev.y - tNext.y);
            T bn = std::hypot(bis.x, bis.y);
            if (bn > T(1e-12)) {
                bis.x /= bn; bis.y /= bn;
                projectRayToBoundary(cell.boundary[r].seg.p0(), bis);
            }
        }
        if (cell.holeIds.size() >= 2) {
            const int nCH = static_cast<int>(cell.holeIds.size());
            struct HPair { int a, b; T d; };
            std::vector<HPair> pairs;
            for (int a = 0; a < nCH; ++a)
                for (int b = a + 1; b < nCH; ++b)
                    pairs.push_back({cell.holeIds[a], cell.holeIds[b], holeCtrs[cell.holeIds[a]].dist(holeCtrs[cell.holeIds[b]])});
            std::sort(pairs.begin(), pairs.end(), [](const HPair& x, const HPair& y) { return x.d < y.d; });
            for (size_t pi = 0; pi < std::min<size_t>(pairs.size(), 6); ++pi) {
                int ha = pairs[pi].a, hb = pairs[pi].b;
                gsPoint2D<T> M(T(0.5) * (holeCtrs[ha].x + holeCtrs[hb].x),
                               T(0.5) * (holeCtrs[ha].y + holeCtrs[hb].y));
                gsPoint2D<T> diff(holeCtrs[hb].x - holeCtrs[ha].x, holeCtrs[hb].y - holeCtrs[ha].y);
                T dl = std::hypot(diff.x, diff.y);
                if (dl > T(1e-12)) {
                    gsPoint2D<T> perp(-diff.y / dl, diff.x / dl);
                    projectRayToBoundary(M, perp);
                    projectRayToBoundary(M, gsPoint2D<T>(-perp.x, -perp.y));
                }
            }
        }

        // 2) Arc-length midpoints & projections on !isCutEdge segments (capped by maxAnchors)
        for (int i = 0; i < K && anchors.size() < maxAnchors; ++i) {
            if (cell.boundary[i].isCutEdge) continue;
            addEdgeAnchor(i, arcLengthSplitParam(*cell.boundary[i].seg.bsp, T(0.5)));
            if (K <= 6) {
                addEdgeAnchor(i, arcLengthSplitParam(*cell.boundary[i].seg.bsp, T(0.30)));
                addEdgeAnchor(i, arcLengthSplitParam(*cell.boundary[i].seg.bsp, T(0.70)));
            }
        }
        if (anchors.size() < maxAnchors) {
            for (int i = 0; i < K && anchors.size() < maxAnchors; ++i) {
                if (cell.boundary[i].isCutEdge) continue;
                for (int r : reflexIdx) {
                    if (r == i || r == (i + 1) % K) continue;
                    addEdgeAnchor(i, cell.boundary[i].seg.closestParam(cell.boundary[r].seg.p0()));
                }
                for (int hid : cell.holeIds) {
                    addEdgeAnchor(i, cell.boundary[i].seg.closestParam(holeCtrs[hid]));
                }
            }
        }

        // Evaluate all valid anchor pairs lazily (zero gsBSpline allocations until top maxTry selection)
        struct CandMeta {
            int ia, ib;
            T offS = T(0);
            gsPoint2D<T> mCtrl;
            std::vector<int> holesL, holesR;
            int tier = 0;
            T score = T(0);
        };
        std::vector<CandMeta> cands;
        const int nA = static_cast<int>(anchors.size());

        auto sameSegEdge = [&](const Anchor& A, const Anchor& B) -> bool {
            if (A.edgeIdx == B.edgeIdx) return true;
            if (A.isVertex && !B.isVertex && ((B.edgeIdx + 1) % K == A.edgeIdx)) return true;
            if (B.isVertex && !A.isVertex && ((A.edgeIdx + 1) % K == B.edgeIdx)) return true;
            if (A.isVertex && B.isVertex) {
                int df = (B.edgeIdx - A.edgeIdx + K) % K;
                if (df <= 1 || df >= K - 1) return true;
            }
            return false;
        };

        auto isHoleInPolyMajority = [&](int hid, const std::vector<gsPoint2D<T>>& poly) -> bool {
            const auto& hp = holePolys[hid];
            int nh = static_cast<int>(hp.size());
            int inCnt = 0;
            for (int k = 0; k < 4; ++k) {
                if (pointInPolyStrict(hp[(k * nh) / 4], poly, T(1e-8))) ++inCnt;
            }
            return inCnt >= 2;
        };

        auto buildPolyLR = [&](const Anchor& A, const Anchor& B,
                               const std::vector<gsPoint2D<T>>& cutPts,
                               std::vector<gsPoint2D<T>>& polyL,
                               std::vector<gsPoint2D<T>>& polyR) {
            polyL.clear(); polyR.clear();
            auto appendArc = [&](int eIdx, int sStart, int sEnd, std::vector<gsPoint2D<T>>& dst) {
                for (int s = sStart; s < sEnd; ++s) dst.push_back(segPts[eIdx][s]);
            };
            int sA = A.isVertex ? 0 : std::max(1, std::min(15, static_cast<int>(std::round(A.u * T(16)))));
            int sB = B.isVertex ? 0 : std::max(1, std::min(15, static_cast<int>(std::round(B.u * T(16)))));
            // polyL: A -> ... -> B along cell boundary, then B -> A along cutPts reversed
            if (!A.isVertex) {
                polyL.push_back(A.pt);
                appendArc(A.edgeIdx, sA + 1, 16, polyL);
            } else {
                appendArc(A.edgeIdx, 0, 16, polyL);
            }
            for (int e = (A.edgeIdx + 1) % K; e != B.edgeIdx; e = (e + 1) % K)
                appendArc(e, 0, 16, polyL);
            if (!B.isVertex)
                appendArc(B.edgeIdx, 0, sB, polyL);
            for (int s = static_cast<int>(cutPts.size()) - 1; s >= 1; --s)
                polyL.push_back(cutPts[s]);

            // polyR: B -> ... -> A along cell boundary, then A -> B along cutPts forward
            if (!B.isVertex) {
                polyR.push_back(B.pt);
                appendArc(B.edgeIdx, sB + 1, 16, polyR);
            } else {
                appendArc(B.edgeIdx, 0, 16, polyR);
            }
            for (int e = (B.edgeIdx + 1) % K; e != A.edgeIdx; e = (e + 1) % K)
                appendArc(e, 0, 16, polyR);
            if (!A.isVertex)
                appendArc(A.edgeIdx, 0, sA, polyR);
            for (int s = 0; s + 1 < static_cast<int>(cutPts.size()); ++s)
                polyR.push_back(cutPts[s]);
        };

        std::vector<gsPoint2D<T>> cutPts(15), bestCutPts(15), polyL, polyR;
        for (int ia = 0; ia < nA; ++ia) {
            const auto& A = anchors[ia];
            for (int ib = ia + 1; ib < nA; ++ib) {
                const auto& B = anchors[ib];
                if (sameSegEdge(A, B)) continue;

                T lenAB = A.pt.dist(B.pt);
                if (lenAB < T(0.02) * domDiag) continue;

                gsPoint2D<T> dirAB((B.pt.x - A.pt.x) / lenAB, (B.pt.y - A.pt.y) / lenAB);
                gsPoint2D<T> perpAB(-dirAB.y, dirAB.x);
                gsPoint2D<T> midAB(T(0.5) * (A.pt.x + B.pt.x), T(0.5) * (A.pt.y + B.pt.y));

                // Fast angle screen for straight cut before polygon/hole tests
                auto checkCutAngles = [&](const gsPoint2D<T>& tCutA, const gsPoint2D<T>& tCutB,
                                          T& mnOut, T& mxOut) -> bool {
                    gsPoint2D<T> tCutArev(-tCutA.x, -tCutA.y);
                    gsPoint2D<T> tCutBrev(-tCutB.x, -tCutB.y);
                    T angL_B = evalCornerBetweenCuts(B.tIn, tCutBrev);
                    T angL_A = evalCornerBetweenCuts(tCutArev, A.tOut);
                    T angR_A = evalCornerBetweenCuts(A.tIn, tCutA);
                    T angR_B = evalCornerBetweenCuts(tCutB, B.tOut);
                    mnOut = std::min(std::min(angL_A, angL_B), std::min(angR_A, angR_B));
                    mxOut = std::max(std::max(angL_A, angL_B), std::max(angR_A, angR_B));
                    return (mnOut >= std::min(minAngTol, T(18)) && mxOut <= T(158));
                };

                T chosenOffS = T(0);
                gsPoint2D<T> chosenMctrl = midAB;
                T bestClear = -T(1);
                T bestMnAng = T(0), bestMxAng = T(180);
                bool foundCutGeom = false;
                bool straightHitOuter = false;

                for (T offS : {T(0), T(0.15), T(-0.15)}) {
                    if (offS != T(0) && (straightHitOuter || (foundCutGeom && bestClear > T(0.025) * domDiag)))
                        break;
                    gsPoint2D<T> Mctrl = (offS == T(0))
                        ? midAB
                        : gsPoint2D<T>(midAB.x + offS * lenAB * perpAB.x, midAB.y + offS * lenAB * perpAB.y);
                    gsPoint2D<T> tCutA = dirAB, tCutB = dirAB;
                    if (offS != T(0)) {
                        T dAx = Mctrl.x - A.pt.x, dAy = Mctrl.y - A.pt.y;
                        T nA0 = std::hypot(dAx, dAy);
                        T dBx = B.pt.x - Mctrl.x, dBy = B.pt.y - Mctrl.y;
                        T nB0 = std::hypot(dBx, dBy);
                        if (nA0 <= T(1e-12) || nB0 <= T(1e-12)) continue;
                        tCutA = gsPoint2D<T>(dAx / nA0, dAy / nA0);
                        tCutB = gsPoint2D<T>(dBx / nB0, dBy / nB0);
                    }
                    T mnA = 0, mxA = 180;
                    if (!checkCutAngles(tCutA, tCutB, mnA, mxA)) continue;

                    for (int s = 0; s <= 14; ++s) {
                        T t = static_cast<T>(s) / T(14);
                        if (offS == T(0)) {
                            cutPts[s] = gsPoint2D<T>((T(1) - t) * A.pt.x + t * B.pt.x,
                                                     (T(1) - t) * A.pt.y + t * B.pt.y);
                        } else {
                            T om = T(1) - t;
                            cutPts[s] = gsPoint2D<T>(om * om * A.pt.x + T(2) * om * t * Mctrl.x + t * t * B.pt.x,
                                                     om * om * A.pt.y + T(2) * om * t * Mctrl.y + t * t * B.pt.y);
                        }
                    }

                    bool insideC = true;
                    for (int s = 1; s < 14; ++s) {
                        if (!pointInPolyStrict(cutPts[s], cPoly, T(1e-6))) { insideC = false; break; }
                    }
                    if (!insideC) {
                        if (offS == T(0)) straightHitOuter = true;
                        continue;
                    }
                    for (int s = 1; s + 2 < 14 && insideC; ++s) {
                        for (size_t e = 0; e < cPoly.size(); ++e) {
                            if (segProperIntersect(cutPts[s], cutPts[s + 1], cPoly[e], cPoly[(e + 1) % cPoly.size()])) {
                                insideC = false; break;
                            }
                        }
                    }
                    if (!insideC) {
                        if (offS == T(0)) straightHitOuter = true;
                        continue;
                    }

                    T minHClr = domDiag;
                    bool hitsHole = false;
                    for (int hid : cell.holeIds) {
                        T hd = segMinDistPoly(cutPts, holePolys[hid], true);
                        if (hd <= T(1e-3) * domDiag) { hitsHole = true; break; }
                        minHClr = std::min(minHClr, hd);
                    }
                    if (!hitsHole && minHClr > bestClear) {
                        bestClear = minHClr;
                        chosenOffS = offS;
                        chosenMctrl = Mctrl;
                        bestMnAng = mnA;
                        bestMxAng = mxA;
                        bestCutPts = cutPts;
                        foundCutGeom = true;
                    }
                }
                if (!foundCutGeom) continue;

                // Compute sub-cell edge counts & orig-edge counts analytically without splitting B-splines
                int cntL = 1, cntR = 1; // +1 for the cut edge itself
                int origL = 0, origR = 0;
                if (!A.isVertex) { ++cntL; ++origL; ++cntR; ++origR; }
                else             { ++cntL; if (!cell.boundary[A.edgeIdx].isCutEdge) ++origL; }
                for (int e = (A.edgeIdx + 1) % K; e != B.edgeIdx; e = (e + 1) % K) {
                    ++cntL;
                    if (!cell.boundary[e].isCutEdge) ++origL;
                }
                if (!B.isVertex) { ++cntL; ++origL; ++cntR; ++origR; }
                else             { ++cntR; if (!cell.boundary[B.edgeIdx].isCutEdge) ++origR; }
                for (int e = (B.edgeIdx + 1) % K; e != A.edgeIdx; e = (e + 1) % K) {
                    ++cntR;
                    if (!cell.boundary[e].isCutEdge) ++origR;
                }
                if (cntL < 3 || cntR < 3) continue;

                // Partition holes using lightweight polyL/polyR built from cached segPts + bestCutPts
                buildPolyLR(A, B, bestCutPts, polyL, polyR);
                std::vector<int> holesL, holesR;
                for (int hid : cell.holeIds) {
                    bool inL = isHoleInPolyMajority(hid, polyL);
                    bool inR = isHoleInPolyMajority(hid, polyR);
                    if (inL && !inR)      holesL.push_back(hid);
                    else if (inR && !inL) holesR.push_back(hid);
                    else                  { holesL.clear(); holesR.clear(); break; }
                }
                if (holesL.size() + holesR.size() != cell.holeIds.size()) continue;

                // Parity guard: never create a 0-hole cell whose edges are ALL isCutEdge == true with ODD count
                if (holesL.empty() && origL == 0 && (cntL % 2 != 0 || cntL < 4)) continue;
                if (holesR.empty() && origR == 0 && (cntR % 2 != 0 || cntR < 4)) continue;

                // Require strict progress
                bool splitsHoles = (!holesL.empty() && !holesR.empty());
                bool fixesReflex = (A.isVertex && angDeg[A.edgeIdx] > T(194)) || (B.isVertex && angDeg[B.edgeIdx] > T(194));
                bool simplifiesK = (cntL < K && cntR < K) ||
                                   (cell.holeIds.size() == 1 && ((holesL.size() == 1 && cntL == 4) ||
                                                                 (holesR.size() == 1 && cntR == 4)));

                if (cell.holeIds.size() >= 2 && !splitsHoles && !fixesReflex && !simplifiesK) continue;
                if (cell.holeIds.size() == 1 && !fixesReflex && !simplifiesK) continue;

                bool chopsZeroHoleQuad = (K >= 10 && reflexIdx.size() >= 2) &&
                                         ((holesL.empty() && cntL == 4) || (holesR.empty() && cntR == 4));
                int tier = (splitsHoles || chopsZeroHoleQuad) ? 0 : 1;
                int hDiff = std::abs(static_cast<int>(holesL.size()) - static_cast<int>(holesR.size()));
                T sc = (splitsHoles ? T(100) - T(12) * hDiff : T(0))
                     + (chopsZeroHoleQuad ? T(95) : T(0))
                     + (fixesReflex ? T(65) : T(0))
                     + T(40) * (bestClear / std::max(T(1e-6), domDiag))
                     + T(0.35) * std::min(bestMnAng, T(180) - bestMxAng)
                     - T(22) * (lenAB / domDiag)
                     + (A.isVertex ? T(6) : T(0)) + (B.isVertex ? T(6) : T(0));
                if (holesL.size() == 1 && cntL == 4) sc += T(15);
                if (holesR.size() == 1 && cntR == 4) sc += T(15);
                if ((holesL.size() >= 2 && origL < 2) || (holesR.size() >= 2 && origR < 2))
                    sc -= T(45);

                cands.push_back({ia, ib, chosenOffS, chosenMctrl, std::move(holesL), std::move(holesR), tier, sc});
            }
        }

        if (cands.empty()) return false;
        std::sort(cands.begin(), cands.end(), [](const CandMeta& x, const CandMeta& y) {
            if (x.tier != y.tier) return x.tier < y.tier;
            return x.score > y.score;
        });

        auto materializeSubCells = [&](const CandMeta& cm, Cell& cL, Cell& cR) -> bool {
            const auto& A = anchors[cm.ia];
            const auto& B = anchors[cm.ib];
            CurveSeg<T> cutSeg = (cm.offS == T(0))
                ? mkLine(A.pt, B.pt)
                : makeCurvedBezierSeg(A.pt, cm.mCtrl, B.pt);

            const auto& seq = cell.boundary;
            int idxA = A.edgeIdx, idxB = B.edgeIdx;
            std::vector<HalfSeg> exp;
            int posA = -1, posB = -1;
            for (int i = 0; i < K; ++i) {
                if (i == idxA && !A.isVertex) {
                    auto sp = seq[i].seg.splitAtParam(A.u);
                    exp.push_back({sp.first, false});
                    posA = static_cast<int>(exp.size());
                    exp.push_back({sp.second, false});
                } else if (i == idxB && !B.isVertex) {
                    auto sp = seq[i].seg.splitAtParam(B.u);
                    exp.push_back({sp.first, false});
                    posB = static_cast<int>(exp.size());
                    exp.push_back({sp.second, false});
                } else {
                    if (i == idxA && A.isVertex) posA = static_cast<int>(exp.size());
                    if (i == idxB && B.isVertex) posB = static_cast<int>(exp.size());
                    exp.push_back(seq[i]);
                }
            }
            if (posA < 0 || posB < 0 || posA == posB) return false;
            const int KE = static_cast<int>(exp.size());
            for (int k = posA; k != posB; k = (k + 1) % KE) cL.boundary.push_back(exp[k]);
            cL.boundary.push_back({cutSeg.reversed(), true});
            for (int k = posB; k != posA; k = (k + 1) % KE) cR.boundary.push_back(exp[k]);
            cR.boundary.push_back({cutSeg, true});
            cL.holeIds = cm.holesL;
            cR.holeIds = cm.holesR;
            return (cL.boundary.size() >= 3 && cR.boundary.size() >= 3);
        };

        const size_t maxTry = (depth == 0) ? 5 : ((depth == 1) ? 3 : 2);
        bool foundSol = false;
        T bestCost = T(1e30);
        std::vector<std::unique_ptr<gsGeometry<T>>> bestSol;

        for (size_t ci = 0; ci < std::min(cands.size(), maxTry); ++ci) {
            if (timeExpired()) break;
            if (depth == 0)
                maxSolverNodes = std::min(450, solverNodes + 140);
            Cell cL, cR;
            if (!materializeSubCells(cands[ci], cL, cR)) continue;
            std::vector<std::unique_ptr<gsGeometry<T>>> pL, pR;
            if (!solveCell(cL, depth + 1, pL)) continue;
            if (!solveCell(cR, depth + 1, pR)) continue;

            std::vector<std::unique_ptr<gsGeometry<T>>> comb;
            comb.reserve(pL.size() + pR.size());
            for (auto& g : pL) comb.push_back(std::move(g));
            for (auto& g : pR) comb.push_back(std::move(g));
            T cost = weightPatches * static_cast<T>(comb.size());
            if (!foundSol || cost < bestCost) {
                foundSol = true;
                bestCost = cost;
                bestSol = std::move(comb);
                if (bestCost <= weightPatches * T(6 * cell.holeIds.size() + 6) || (depth == 0 && ci >= 1))
                    break;
            }
        }

        if (foundSol) {
            outP = std::move(bestSol);
            return true;
        }
        return false;
    };

    std::vector<std::unique_ptr<gsGeometry<T>>> allPatches;
    if (!solveCell(rootCell, 0, allPatches) || allPatches.empty())
        return false;

    outMp.clear();
    for (auto& gp : allPatches)
        outMp.addPatch(typename gsGeometry<T>::uPtr(gp.release()));
    outLabel = "RecursiveCellDecomposition";
    (void)weightDistortion;
    return true;
}



} // namespace internal_splitter
} // namespace gismo
