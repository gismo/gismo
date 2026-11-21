/** @file gsPatchSplitterTriMesh.hpp
    @brief Constrained Delaunay Triangulation, edge-based Lawson flips, interior smoothing,
           curvature-based refinement, obtuse bisection, and quad-union pairing for gsPatchSplitter.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitterGeom.hpp>

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <iostream>

namespace gismo {
namespace internal_splitter {

template <class T>
struct TriMeshData {
    struct Tri { int v[3]; int e[3]; };
    struct Quad { int v[4]; int e[4]; };
    struct EdgeRec {
        int a = -1, b = -1;
        CurveSeg<T> seg;
        CurveSeg<T> h0, h1;
        gsPoint2D<T> mid;
        bool splitDone = false;
        bool isBoundary = false;
        int loopIdx = -1;
        T splitFrac = T(0.5);
    };

    struct BoundaryVertexInfo {
        int loopIdx = -1;
        int curveIdx = -1;
        T paramU = T(0);
    };

    struct BoundaryLoopRing {
        std::vector<int> vertices; // ordered vertex indices along the loop
        std::vector<int> edges;    // edges[i] is edgeId connecting vertices[i] -> vertices[(i+1)%M]
    };

    std::vector<gsPoint2D<T>> pts;
    std::vector<bool> isBoundaryPt;
    std::vector<bool> isSmoothBnd;
    std::vector<bool> isSharpCorner;
    std::vector<T> cornerAngleDeg;
    std::vector<gsPoint2D<T>> bndTanIn;
    std::vector<gsPoint2D<T>> bndTanOut;
    std::vector<BoundaryVertexInfo> bndInfo;
    std::vector<BoundaryLoopRing> loopRings;
    std::vector<EdgeRec> edges;
    std::vector<Tri> tris;
    std::vector<Quad> quads;
    std::vector<std::vector<gsPoint2D<T>>> loopPolys;
};

template <class T>
inline bool isPointStrictlyInsideDomain(const gsPoint2D<T>& p,
                                        const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
                                        T tol = T(1e-6))
{
    if (loopPolys.empty()) return false;
    if (!pointInPolyStrict(p, loopPolys[0], tol)) return false;
    for (size_t h = 1; h < loopPolys.size(); ++h) {
        // Point must not be inside or on hole boundary
        const auto& hp = loopPolys[h];
        int n = static_cast<int>(hp.size());
        T tol2 = tol * tol;
        for (int i = 0; i < n; ++i) {
            const auto& a = hp[i];
            const auto& b = hp[(i + 1) % n];
            T dx = b.x - a.x, dy = b.y - a.y;
            T L2 = dx * dx + dy * dy;
            if (L2 > T(1e-20)) {
                T t = std::max(T(0), std::min(T(1), ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2));
                T qx = a.x + t * dx - p.x;
                T qy = a.y + t * dy - p.y;
                if (qx * qx + qy * qy <= tol2) return false;
            }
        }
        bool inHole = false;
        for (int i = 0, j = n - 1; i < n; j = i++) {
            if (((hp[i].y > p.y) != (hp[j].y > p.y)) &&
                (p.x < (hp[j].x - hp[i].x) * (p.y - hp[i].y) / (hp[j].y - hp[i].y + T(1e-30)) + hp[i].x))
                inHole = !inHole;
        }
        if (inHole) return false;
    }
    return true;
}

template <class T>
inline bool isSegmentInsideDomain(const gsPoint2D<T>& a, const gsPoint2D<T>& b,
                                  const std::vector<std::vector<gsPoint2D<T>>>& loopPolys)
{
    const T abx = b.x - a.x, aby = b.y - a.y;
    for (const auto& poly : loopPolys) {
        int m = static_cast<int>(poly.size());
        for (int i = 0; i < m; ++i) {
            const auto& c = poly[i];
            const auto& d = poly[(i + 1) % m];
            const T cdx = d.x - c.x, cdy = d.y - c.y;
            const T denom = abx * cdy - aby * cdx;
            if (std::abs(denom) <= T(1e-14)) continue;
            const T t = ((c.x - a.x) * cdy - (c.y - a.y) * cdx) / denom;
            const T u = ((c.x - a.x) * aby - (c.y - a.y) * abx) / denom;
            if (t > T(0.02) && t < T(0.98) && u >= -T(1e-6) && u <= T(1) + T(1e-6))
                return false;
        }
    }
    for (int k = 1; k <= 5; ++k) {
        T t = static_cast<T>(k) / T(6);
        gsPoint2D<T> mid((T(1) - t) * a.x + t * b.x, (T(1) - t) * a.y + t * b.y);
        if (!isPointStrictlyInsideDomain(mid, loopPolys, T(1e-6)))
            return false;
    }
    return true;
}

// Check if vector d lies strictly inside CCW sector from ray u to ray v around origin
template <class T>
inline bool vecStrictlyInSector(const gsPoint2D<T>& u, const gsPoint2D<T>& v, const gsPoint2D<T>& d)
{
    T cu = u.x * d.y - u.y * d.x;
    T cv = d.x * v.y - d.y * v.x;
    T uv = u.x * v.y - u.y * v.x;
    if (uv >= T(0))
        return (cu > T(1e-9)) && (cv > T(1e-9));
    else
        return (cu > T(1e-9)) || (cv > T(1e-9));
}

// Compute minimum and maximum interior angles (deg) of triangle (A, B, C).
template <class T>
inline void triMinMaxAnglesDeg(const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
                               T& minAngOut, T& maxAngOut)
{
    const T radToDeg = T(180) / T(3.14159265358979323846);
    auto angAt = [&](const gsPoint2D<T>& P, const gsPoint2D<T>& Q, const gsPoint2D<T>& R) -> T {
        T ux = Q.x - P.x, uy = Q.y - P.y;
        T vx = R.x - P.x, vy = R.y - P.y;
        T lu = std::hypot(ux, uy), lv = std::hypot(vx, vy);
        if (lu <= T(1e-14) || lv <= T(1e-14)) return T(0);
        T c = std::max(T(-1), std::min(T(1), (ux * vx + uy * vy) / (lu * lv)));
        return std::acos(c) * radToDeg;
    };
    T a0 = angAt(A, B, C), a1 = angAt(B, C, A), a2 = angAt(C, A, B);
    minAngOut = std::min(a0, std::min(a1, a2));
    maxAngOut = std::max(a0, std::max(a1, a2));
}

template <class T>
inline bool buildConstrainedDomainTriangulation(
    const gsPlanarDomain<T>& dom,
    int degree,
    T targetH,
    bool refineEdges,
    int earMode,
    TriMeshData<T>& outMesh,
    T maxBendRad = T(0.44),
    bool refineHolePerim = true,
    bool addFeatureProj = false,
    bool ensureEvenParity = false)
{
    if (dom.numLoops() == 0) return false;
    const int nLoops = dom.numLoops();

    outMesh = TriMeshData<T>();
    outMesh.loopPolys.resize(nLoops);
    for (int l = 0; l < nLoops; ++l)
        outMesh.loopPolys[l] = sampleCurveLoopPoly(dom.loop(l), 28);

    // Prepare boundary segments for each loop (ensuring >= 4 segments per loop)
    struct LoopSeg {
        CurveSeg<T> seg;
        int loopIdx = 0;
        int curveIdx = 0;
        T paramU0 = T(0);
        T paramU1 = T(1);
    };
    std::vector<std::vector<LoopSeg>> loopSegs(nLoops);
    std::vector<std::vector<LoopSeg>> rawLoops(nLoops);
    const gsMatrix<T> domBbox = dom.boundingBox();
    const T domDiag = std::hypot(domBbox(0, 1) - domBbox(0, 0), domBbox(1, 1) - domBbox(1, 0));

    for (int l = 0; l < nLoops; ++l) {
        const auto& lp = dom.loop(l);
        const int K = lp.size();
        if (K == 0) return false;

        std::vector<LoopSeg> raw;
        if (K == 1) {
            auto sH = CurveSeg<T>(lp.curve(0), degree).splitArc(T(0.5));
            auto s1 = sH.first.splitArc(T(0.5));
            auto s2 = sH.second.splitArc(T(0.5));
            raw = { {s1.first, l, 0, T(0.0), T(0.25)},
                    {s1.second, l, 0, T(0.25), T(0.5)},
                    {s2.first, l, 0, T(0.5), T(0.75)},
                    {s2.second, l, 0, T(0.75), T(1.0)} };
        } else if (K == 2) {
            auto s0 = CurveSeg<T>(lp.curve(0), degree).splitArc(T(0.5));
            auto s1 = CurveSeg<T>(lp.curve(1), degree).splitArc(T(0.5));
            raw = { {s0.first, l, 0, T(0.0), T(0.5)},
                    {s0.second, l, 0, T(0.5), T(1.0)},
                    {s1.first, l, 1, T(0.0), T(0.5)},
                    {s1.second, l, 1, T(0.5), T(1.0)} };
        } else {
            for (int i = 0; i < K; ++i)
                raw.push_back({CurveSeg<T>(lp.curve(i), degree), l, i, T(0.0), T(1.0)});
        }
        rawLoops[l] = std::move(raw);
    }

    // Collect feature points (all hole vertices and outer vertices) for orthogonal projection splits
    std::vector<gsPoint2D<T>> featurePts;
    if (addFeatureProj && nLoops > 1) {
        for (int l = 0; l < nLoops; ++l)
            for (const auto& rseg : rawLoops[l])
                featurePts.push_back(rseg.seg.p0());
    }

    for (int l = 0; l < nLoops; ++l) {
        const auto& raw = rawLoops[l];

        // Gather split parameters per raw segment from uniform/integrated-bend refinement + feature projections
        std::vector<std::vector<T>> segSplits(raw.size());
        if (domDiag > T(1e-12)) {
            T loopPerim = T(0);
            std::vector<T> segLens(raw.size(), T(0));
            std::vector<T> segTurns(raw.size(), T(0));
            for (size_t si = 0; si < raw.size(); ++si) {
                const auto& bsp = *raw[si].seg.bsp;
                const int NS = 24;
                gsMatrix<T> uMat(1, NS + 1), ptsMat;
                for (int ks = 0; ks <= NS; ++ks)
                    uMat(0, ks) = static_cast<T>(ks) / static_cast<T>(NS);
                bsp.eval_into(uMat, ptsMat);
                T len = T(0), totTurn = T(0);
                gsPoint2D<T> prevTan(0, 0);
                bool hasPrevTan = false;
                for (int ks = 1; ks <= NS; ++ks) {
                    T dx = ptsMat(0, ks) - ptsMat(0, ks - 1);
                    T dy = ptsMat(1, ks) - ptsMat(1, ks - 1);
                    T ds = std::hypot(dx, dy);
                    len += ds;
                    if (ds > T(1e-14)) {
                        gsPoint2D<T> curTan(dx / ds, dy / ds);
                        if (hasPrevTan) {
                            T cr = prevTan.x * curTan.y - prevTan.y * curTan.x;
                            T dt = prevTan.x * curTan.x + prevTan.y * curTan.y;
                            totTurn += std::abs(std::atan2(cr, dt));
                        }
                        prevTan = curTan;
                        hasPrevTan = true;
                    }
                }
                segLens[si] = len;
                segTurns[si] = totTurn;
                loopPerim += len;
            }
            if (refineEdges && targetH > T(1e-8)) {
                T localH = targetH;
                if (l > 0 && loopPerim > T(1e-8) && refineHolePerim)
                    localH = std::min(targetH, loopPerim / T(8));

                for (size_t si = 0; si < raw.size(); ++si) {
                    int nLen  = static_cast<int>(std::ceil(segLens[si] / localH));
                    int nBend = static_cast<int>(std::ceil(segTurns[si] / maxBendRad));
                    int nSub = std::max(1, std::min(24, std::max(nLen, nBend)));
                    for (int k = 1; k < nSub; ++k) {
                        T uArc = arcLengthSplitParam(*raw[si].seg.bsp, static_cast<T>(k) / static_cast<T>(nSub));
                        segSplits[si].push_back(uArc);
                    }
                }
            } else {
                // Even in coarse pass, subdivide S-curves/hairpins whose integrated turning exceeds ~106 deg
                for (size_t si = 0; si < raw.size(); ++si) {
                    if (segTurns[si] > T(1.85)) {
                        int nSub = std::max(2, std::min(16, static_cast<int>(std::ceil(segTurns[si] / T(1.15)))));
                        for (int k = 1; k < nSub; ++k) {
                            T uArc = arcLengthSplitParam(*raw[si].seg.bsp, static_cast<T>(k) / static_cast<T>(nSub));
                            segSplits[si].push_back(uArc);
                        }
                    }
                }
            }
        }

        if (addFeatureProj && nLoops > 1 && domDiag > T(1e-12)) {
            for (size_t si = 0; si < raw.size(); ++si) {
                const auto& s = raw[si].seg;
                gsPoint2D<T> pm = evalCurvePoint(*s.bsp, T(0.5));
                T sLen = s.p0().dist(pm) + pm.dist(s.p1());
                if (sLen < (l == 0 ? T(0.16) : T(0.12)) * domDiag) continue;

                // Pick up to 2 best nearly-perpendicular interior feature projections on this segment
                std::vector<std::pair<T, T>> projCands; // (dist, u)
                for (const auto& q : featurePts) {
                    if (q.dist(s.p0()) < T(0.06) * domDiag || q.dist(s.p1()) < T(0.06) * domDiag)
                        continue;
                    T u = s.closestParam(q);
                    if (u < T(0.16) || u > T(0.84)) continue;
                    gsPoint2D<T> spPt = evalCurvePoint(*s.bsp, u);
                    T dist = q.dist(spPt);
                    if (dist < T(1e-4) || dist > T(0.35) * domDiag) continue;
                    auto td = evalCurveDeriv(*s.bsp, u);
                    T tn = std::hypot(td.x, td.y);
                    if (tn <= T(1e-14)) continue;
                    T cosPar = std::abs((q.x - spPt.x) * td.x + (q.y - spPt.y) * td.y) / (dist * tn);
                    if (cosPar > T(0.32)) continue;
                    if (!isSegmentInsideDomain(q, spPt, outMesh.loopPolys)) continue;
                    projCands.push_back({dist, u});
                }
                std::sort(projCands.begin(), projCands.end(),
                          [](const std::pair<T, T>& a, const std::pair<T, T>& b) { return a.first < b.first; });
                int added = 0;
                for (const auto& pc : projCands) {
                    if (added >= 2) break;
                    bool tooClose = false;
                    for (T ex : segSplits[si]) {
                        if (std::abs(ex - pc.second) < T(0.14)) { tooClose = true; break; }
                    }
                    if (!tooClose) {
                        segSplits[si].push_back(pc.second);
                        ++added;
                    }
                }
            }
        }

        for (size_t si = 0; si < raw.size(); ++si) {
            auto& spParams = segSplits[si];
            if (spParams.empty()) {
                loopSegs[l].push_back(raw[si]);
                continue;
            }
            std::sort(spParams.begin(), spParams.end());
            std::vector<T> filt;
            for (T u : spParams) {
                if (u <= T(0.08) || u >= T(0.92)) continue;
                if (filt.empty() || u - filt.back() >= T(0.08))
                    filt.push_back(u);
            }
            if (filt.empty()) {
                loopSegs[l].push_back(raw[si]);
                continue;
            }
            CurveSeg<T> rem = raw[si].seg;
            T prevU = T(0);
            T u0 = raw[si].paramU0, u1 = raw[si].paramU1;
            for (T uAbs : filt) {
                T locU = (uAbs - prevU) / std::max(T(1e-9), T(1) - prevU);
                locU = std::max(T(0.05), std::min(T(0.95), locU));
                auto sp = rem.splitAtParam(locU);
                T curU = u0 + uAbs * (u1 - u0);
                loopSegs[l].push_back({sp.first, l, raw[si].curveIdx, u0 + prevU * (u1 - u0), curU});
                rem = sp.second;
                prevU = uAbs;
            }
            loopSegs[l].push_back({rem, l, raw[si].curveIdx, u0 + prevU * (u1 - u0), u1});
        }

        // Rebuild loopPolys directly from loopSegs so every ring vertex coincides with loopPolys
        outMesh.loopPolys[l].clear();
        for (const auto& ls : loopSegs[l]) {
            for (int k = 0; k < 16; ++k) {
                T u = static_cast<T>(k) / T(16);
                outMesh.loopPolys[l].push_back(evalCurvePoint(*ls.seg.bsp, u));
            }
        }
    }

    if (ensureEvenParity) {
        size_t totalBnd = 0;
        int bestL = -1, bestSi = -1;
        T bestScore = T(-1);
        for (int l = 0; l < nLoops; ++l) {
            totalBnd += loopSegs[l].size();
            for (int si = 0; si < static_cast<int>(loopSegs[l].size()); ++si) {
                T arcLen = T(0);
                gsPoint2D<T> prevPt = evalCurvePoint(*loopSegs[l][si].seg.bsp, T(0));
                for (int ks = 1; ks <= 8; ++ks) {
                    gsPoint2D<T> curPt = evalCurvePoint(*loopSegs[l][si].seg.bsp, static_cast<T>(ks) / T(8));
                    arcLen += prevPt.dist(curPt);
                    prevPt = curPt;
                }
                T sc = arcLen * (l == 0 ? T(1.05) : T(1.0));
                if (sc > bestScore) { bestScore = sc; bestL = l; bestSi = si; }
            }
        }
        if ((totalBnd % 2) == 1 && bestL >= 0 && bestSi >= 0) {
            const auto orig = loopSegs[bestL][bestSi];
            T tMid = arcLengthSplitParam(*orig.seg.bsp, T(0.5));
            auto sp = orig.seg.splitAtParam(tMid);
            T uMid = orig.paramU0 + tMid * (orig.paramU1 - orig.paramU0);
            LoopSeg firstHalf = {sp.first, orig.loopIdx, orig.curveIdx, orig.paramU0, uMid};
            LoopSeg secondHalf = {sp.second, orig.loopIdx, orig.curveIdx, uMid, orig.paramU1};
            loopSegs[bestL][bestSi] = firstHalf;
            loopSegs[bestL].insert(loopSegs[bestL].begin() + bestSi + 1, secondHalf);

            outMesh.loopPolys[bestL].clear();
            for (const auto& ls : loopSegs[bestL]) {
                for (int k = 0; k < 16; ++k) {
                    T u = static_cast<T>(k) / T(16);
                    outMesh.loopPolys[bestL].push_back(evalCurvePoint(*ls.seg.bsp, u));
                }
            }
        }
    }

    // Populate vertices and boundary edges for every loop
    struct RingNode {
        int ptIdx;
        CurveSeg<T> outSeg;
        int edgeId;
        gsPoint2D<T> q25, q50, q75;
    };
    auto makeRingNode = [](int ptIdx, const CurveSeg<T>& seg, int eid) -> RingNode {
        RingNode rn{ptIdx, seg, eid, seg.p0(), seg.p0(), seg.p1()};
        if (seg.bsp) {
            gsMatrix<T> ev, par(1, 3);
            par(0, 0) = T(0.25); par(0, 1) = T(0.5); par(0, 2) = T(0.75);
            seg.bsp->eval_into(par, ev);
            rn.q25 = gsPoint2D<T>(ev(0, 0), ev(1, 0));
            rn.q50 = gsPoint2D<T>(ev(0, 1), ev(1, 1));
            rn.q75 = gsPoint2D<T>(ev(0, 2), ev(1, 2));
        }
        return rn;
    };
    std::vector<std::vector<RingNode>> rings(nLoops);
    std::map<std::pair<int, int>, int> edgeMap;

    auto addUndirEdge = [&](int a, int b, const CurveSeg<T>& seg, bool isBdr, int loopIdx) -> int {
        int lo = std::min(a, b), hi = std::max(a, b);
        auto it = edgeMap.find({lo, hi});
        if (it != edgeMap.end()) return it->second;
        int id = static_cast<int>(outMesh.edges.size());
        // EdgeRec stores seg oriented a -> b
        outMesh.edges.push_back({a, b, seg, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, isBdr, loopIdx});
        edgeMap[{lo, hi}] = id;
        return id;
    };

    const T radToDeg = T(180) / T(3.14159265358979323846);
    outMesh.loopRings.clear();
    outMesh.loopRings.resize(nLoops);
    for (int l = 0; l < nLoops; ++l) {
        int M = static_cast<int>(loopSegs[l].size());
        outMesh.loopRings[l].vertices.reserve(M);
        outMesh.loopRings[l].edges.reserve(M);
        int baseIdx = static_cast<int>(outMesh.pts.size());
        T signedArea2 = T(0);
        {
            const auto& poly = outMesh.loopPolys[l];
            for (size_t k = 0; k < poly.size(); ++k) {
                size_t j = (k + 1) % poly.size();
                signedArea2 += poly[k].x * poly[j].y - poly[j].x * poly[k].y;
            }
        }
        // Outer loop CCW (area > 0) -> domain on left; Hole loop CW (area < 0) -> domain on left
        const bool domainOnLeft = (l == 0) ? (signedArea2 >= T(0)) : (signedArea2 <= T(0));
        for (int i = 0; i < M; ++i) {
            outMesh.pts.push_back(loopSegs[l][i].seg.p0());
            outMesh.isBoundaryPt.push_back(true);
            outMesh.bndInfo.push_back({l, loopSegs[l][i].curveIdx, loopSegs[l][i].paramU0});

            const auto& prevBsp = *loopSegs[l][(i + M - 1) % M].seg.bsp;
            const auto& nextBsp = *loopSegs[l][i].seg.bsp;
            auto din  = evalCurveDeriv(prevBsp, T(1));
            auto dout = evalCurveDeriv(nextBsp, T(0));
            T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
            gsPoint2D<T> tIn  = (n1 > T(1e-14)) ? gsPoint2D<T>(din.x / n1, din.y / n1) : gsPoint2D<T>(1, 0);
            gsPoint2D<T> tOut = (n2 > T(1e-14)) ? gsPoint2D<T>(dout.x / n2, dout.y / n2) : gsPoint2D<T>(1, 0);
            T cr = tIn.x * tOut.y - tIn.y * tOut.x;
            T dt = tIn.x * tOut.x + tIn.y * tOut.y;
            T turnRad = std::atan2(cr, dt);
            bool smooth = (std::abs(turnRad) <= T(0.10));
            // Interior domain angle at vertex i when domain lies on the left of the oriented boundary
            T intAngDeg = domainOnLeft ? (T(180) - turnRad * radToDeg) : (T(180) + turnRad * radToDeg);
            outMesh.isSmoothBnd.push_back(smooth);
            outMesh.isSharpCorner.push_back(std::abs(turnRad) > T(0.45)); // > ~25 deg turn
            outMesh.cornerAngleDeg.push_back(intAngDeg);
            outMesh.bndTanIn.push_back(tIn);
            outMesh.bndTanOut.push_back(tOut);
        }
        for (int i = 0; i < M; ++i) {
            int u = baseIdx + i;
            int v = baseIdx + (i + 1) % M;
            int eid = addUndirEdge(u, v, loopSegs[l][i].seg, true, l);
            rings[l].push_back(makeRingNode(u, loopSegs[l][i].seg, eid));
            outMesh.loopRings[l].vertices.push_back(u);
            outMesh.loopRings[l].edges.push_back(eid);
        }
    }

    // Bridge holes into outer ring (ordered by descending max_x)
    std::vector<RingNode> activeRing = rings[0];
    std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> bridges;

    if (nLoops > 1) {
        struct HoleOrder { int loopIdx; T maxX; int rightPos; };
        std::vector<HoleOrder> holes;
        for (int l = 1; l < nLoops; ++l) {
            T mx = -std::numeric_limits<T>::max();
            int rp = 0;
            for (int k = 0; k < static_cast<int>(rings[l].size()); ++k) {
                T x = outMesh.pts[rings[l][k].ptIdx].x;
                if (x > mx) { mx = x; rp = k; }
            }
            holes.push_back({l, mx, rp});
        }
        std::sort(holes.begin(), holes.end(),
                  [](const HoleOrder& a, const HoleOrder& b) { return a.maxX > b.maxX; });

        for (const auto& ho : holes) {
            const auto& hRing = rings[ho.loopIdx];
            int hSize = static_cast<int>(hRing.size());
            int bestHPos = -1, bestOuterPos = -1;

            // Score candidate bridges in all directions by distance and departure angle (avoid grazing < 14 deg slits)
            std::vector<int> ringCount(outMesh.pts.size(), 0);
            for (const auto& rn : activeRing) ++ringCount[rn.ptIdx];

            T bestCost = std::numeric_limits<T>::max();
            for (int pass = 0; pass < 3 && bestHPos < 0; ++pass) {
                for (int hp = 0; hp < hSize; ++hp) {
                    int hmIdx = hRing[hp].ptIdx;
                    const auto& M = outMesh.pts[hmIdx];
                    const auto& mIn  = outMesh.bndTanIn[hmIdx];
                    const auto& mOut = outMesh.bndTanOut[hmIdx];

                    const int actSz = static_cast<int>(activeRing.size());
                    for (int op = 0; op < actSz; ++op) {
                        int ovIdx = activeRing[op].ptIdx;
                        const auto& V = outMesh.pts[ovIdx];
                        if (pass == 0 && V.x < M.x - T(1e-6)) continue;

                        T d = std::hypot(V.x - M.x, V.y - M.y);
                        if (d <= T(1e-12) || d >= bestCost) continue;

                        gsPoint2D<T> dirMV(V.x - M.x, V.y - M.y);
                        gsPoint2D<T> dirVM(-dirMV.x, -dirMV.y);
                        if (ringCount[ovIdx] > 1) {
                            const auto& vNext = outMesh.pts[activeRing[(op + 1) % actSz].ptIdx];
                            const auto& vPrev = outMesh.pts[activeRing[(op + actSz - 1) % actSz].ptIdx];
                            if (!vecStrictlyInSector(gsPoint2D<T>(vNext.x - V.x, vNext.y - V.y),
                                                     gsPoint2D<T>(vPrev.x - V.x, vPrev.y - V.y), dirVM))
                                continue;
                        }

                        T ux = dirMV.x / d, uy = dirMV.y / d;
                        T s1 = std::abs(ux * mIn.y - uy * mIn.x);
                        T s2 = std::abs(ux * mOut.y - uy * mOut.x);
                        T s3 = T(1), s4 = T(1);
                        if (outMesh.isBoundaryPt[ovIdx]) {
                            const auto& vIn  = outMesh.bndTanIn[ovIdx];
                            const auto& vOut = outMesh.bndTanOut[ovIdx];
                            s3 = std::abs(ux * vIn.y - uy * vIn.x);
                            s4 = std::abs(ux * vOut.y - uy * vOut.x);
                        }
                        T minSin = std::min(std::min(s1, s2), std::min(s3, s4));
                        if (pass <= 1 && minSin < T(0.24)) continue;

                        T angFactor = (minSin < T(0.24))
                            ? T(8.0)
                            : (T(1.0) + T(2.0) * (T(1.0) - minSin) * (T(1.0) - minSin));
                        T cost = d * angFactor;
                        if (cost >= bestCost) continue;

                        if (!isSegmentInsideDomain(M, V, outMesh.loopPolys)) continue;

                        bool crossBridge = false;
                        for (const auto& br : bridges) {
                            if (segProperIntersect(M, V, br.first, br.second)) {
                                crossBridge = true;
                                break;
                            }
                        }
                        if (crossBridge) continue;

                        bestCost = cost;
                        bestHPos = hp;
                        bestOuterPos = op;
                    }
                }
            }

            if (bestHPos < 0 || bestOuterPos < 0) {
                int splitHPos = -1, splitActOp = -1;
                T splitParam = T(0.5), bestSplitCost = std::numeric_limits<T>::max();
                const int actSz = static_cast<int>(activeRing.size());
                for (int hp = 0; hp < hSize; ++hp) {
                    int hmIdx = hRing[hp].ptIdx;
                    const auto& M = outMesh.pts[hmIdx];
                    for (int op = 0; op < actSz; ++op) {
                        const int eOld = activeRing[op].edgeId;
                        if (eOld < 0 || eOld >= static_cast<int>(outMesh.edges.size()) || !outMesh.edges[eOld].isBoundary)
                            continue;
                        if (!activeRing[op].outSeg.bsp) continue;
                        T uCand = std::max(T(0.20), std::min(T(0.80), activeRing[op].outSeg.closestParam(M)));
                        gsPoint2D<T> V = evalCurvePoint(*activeRing[op].outSeg.bsp, uCand);
                        T d = std::hypot(V.x - M.x, V.y - M.y);
                        if (d <= T(1e-12) || d >= bestSplitCost) continue;
                        if (!isSegmentInsideDomain(M, V, outMesh.loopPolys)) continue;
                        bool crossBridge = false;
                        for (const auto& br : bridges) {
                            if (segProperIntersect(M, V, br.first, br.second)) { crossBridge = true; break; }
                        }
                        if (crossBridge) continue;
                        bestSplitCost = d;
                        splitHPos = hp;
                        splitActOp = op;
                        splitParam = uCand;
                    }
                }
                if (splitHPos < 0 || splitActOp < 0) return false;

                int eOld = activeRing[splitActOp].edgeId;
                int uIdx = activeRing[splitActOp].ptIdx;
                int vIdxOld = activeRing[(splitActOp + 1) % actSz].ptIdx;
                auto sp = activeRing[splitActOp].outSeg.splitAtParam(splitParam);
                gsPoint2D<T> midPt = sp.first.p1();
                int newIdx = static_cast<int>(outMesh.pts.size());
                outMesh.pts.push_back(midPt);
                outMesh.isBoundaryPt.push_back(true);
                outMesh.isSmoothBnd.push_back(true);
                outMesh.isSharpCorner.push_back(false);
                outMesh.cornerAngleDeg.push_back(T(180));
                auto din = evalCurveDeriv(*sp.first.bsp, T(1));
                T dn = std::hypot(din.x, din.y);
                gsPoint2D<T> tMid = (dn > T(1e-14)) ? gsPoint2D<T>(din.x / dn, din.y / dn) : gsPoint2D<T>(1, 0);
                outMesh.bndTanIn.push_back(tMid);
                outMesh.bndTanOut.push_back(tMid);

                int loopIdxOld = outMesh.edges[eOld].loopIdx;
                bool fwdOld = (outMesh.edges[eOld].a == uIdx);
                outMesh.edges[eOld] = fwdOld
                    ? typename TriMeshData<T>::EdgeRec{uIdx, newIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, loopIdxOld}
                    : typename TriMeshData<T>::EdgeRec{newIdx, uIdx, sp.first.reversed(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, loopIdxOld};
                edgeMap.erase({std::min(uIdx, vIdxOld), std::max(uIdx, vIdxOld)});
                edgeMap[{std::min(uIdx, newIdx), std::max(uIdx, newIdx)}] = eOld;
                int eSecond = addUndirEdge(newIdx, vIdxOld, sp.second, true, loopIdxOld);
                outMesh.bndInfo.push_back({loopIdxOld, -1, T(0)});
                if (loopIdxOld >= 0 && loopIdxOld < static_cast<int>(outMesh.loopRings.size())) {
                    auto& lr = outMesh.loopRings[loopIdxOld];
                    for (size_t ri = 0; ri < lr.edges.size(); ++ri) {
                        if (lr.edges[ri] == eOld) {
                            lr.vertices.insert(lr.vertices.begin() + static_cast<long>(ri + 1), newIdx);
                            lr.edges.insert(lr.edges.begin() + static_cast<long>(fwdOld ? (ri + 1) : ri), eSecond);
                            break;
                        }
                    }
                }
                activeRing[splitActOp] = makeRingNode(uIdx, sp.first, eOld);
                activeRing.insert(activeRing.begin() + splitActOp + 1, makeRingNode(newIdx, sp.second, eSecond));
                bestHPos = splitHPos;
                bestOuterPos = splitActOp + 1;
            }

            int vIdx = activeRing[bestOuterPos].ptIdx;
            int mIdx = hRing[bestHPos].ptIdx;
            CurveSeg<T> segVM = makeStraightSeg(outMesh.pts[vIdx], outMesh.pts[mIdx], degree);
            CurveSeg<T> segMV = segVM.reversed();
            int bridgeEid = addUndirEdge(vIdx, mIdx, segVM, false, -1);
            bridges.push_back({outMesh.pts[vIdx], outMesh.pts[mIdx]});

            // Splice hole ring into activeRing at bestOuterPos:
            // ... -> [V, segVM] -> [M, hRing...] -> [M, segMV] -> [V, oldOutSeg] -> ...
            std::vector<RingNode> newRing;
            newRing.reserve(activeRing.size() + hRing.size() + 2);
            for (int k = 0; k < bestOuterPos; ++k)
                newRing.push_back(activeRing[k]);

            newRing.push_back(makeRingNode(vIdx, segVM, bridgeEid));
            for (int k = 0; k < hSize; ++k)
                newRing.push_back(hRing[(bestHPos + k) % hSize]);
            newRing.push_back(makeRingNode(mIdx, segMV, bridgeEid));

            for (int k = bestOuterPos; k < static_cast<int>(activeRing.size()); ++k)
                newRing.push_back(activeRing[k]);

            activeRing = std::move(newRing);
        }
    }

    auto getOrAddDiag = [&](int a, int b) -> int {
        int lo = std::min(a, b), hi = std::max(a, b);
        auto it = edgeMap.find({lo, hi});
        if (it != edgeMap.end()) return it->second;
        return addUndirEdge(lo, hi, makeStraightSeg(outMesh.pts[lo], outMesh.pts[hi], degree), false, -1);
    };

    auto pointInTri = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B,
                          const gsPoint2D<T>& C, const gsPoint2D<T>& P) -> bool {
        return cross2(A, B, P) >= -1e-9 && cross2(B, C, P) >= -1e-9 && cross2(C, A, P) >= -1e-9;
    };

    std::vector<int> triCount(outMesh.pts.size(), 0);
    int subDivs = 0;

    while (activeRing.size() > 3) {
        int m = static_cast<int>(activeRing.size());
        int bestIdx = -1;
        T bestScore = -T(1e30);

        for (int pass = 0; pass < 2 && bestIdx < 0; ++pass) {
            for (int idx = 0; idx < m; ++idx) {
                int iPrev = (idx + m - 1) % m;
                int iNext = (idx + 1) % m;
                int ia = activeRing[iPrev].ptIdx;
                int ib = activeRing[idx].ptIdx;
                int ic = activeRing[iNext].ptIdx;

                if (ia == ic || ia == ib || ib == ic) continue;
                if (pass == 0 && outMesh.isSmoothBnd[ib] && triCount[ib] == 0) continue;

                const auto& A = outMesh.pts[ia];
                const auto& B = outMesh.pts[ib];
                const auto& C = outMesh.pts[ic];
                T area2 = cross2(A, B, C);
                if (area2 <= T(1e-8)) continue;

                // Ensure diagonal A->C does not cross any boundary or active ring chord
                if (!isSegmentInsideDomain(A, C, outMesh.loopPolys)) continue;

                bool crossesRing = false;
                for (int k = 0; k < m; ++k) {
                    int ka = activeRing[k].ptIdx;
                    int kb = activeRing[(k + 1) % m].ptIdx;
                    if (ka == ia || ka == ic || kb == ia || kb == ic) continue;
                    if (segProperIntersect(A, C, outMesh.pts[ka], outMesh.pts[kb])) {
                        crossesRing = true;
                        break;
                    }
                }
                if (crossesRing) continue;

                // If vertex ib appears multiple times in activeRing (bridge endpoint), check wedge clearance
                bool wedgeOk = true;
                gsPoint2D<T> rayC(C.x - B.x, C.y - B.y);
                gsPoint2D<T> rayA(A.x - B.x, A.y - B.y);
                for (int k = 0; k < m && wedgeOk; ++k) {
                    if (k == idx || activeRing[k].ptIdx != ib) continue;
                    int kp = activeRing[(k + m - 1) % m].ptIdx;
                    int kn = activeRing[(k + 1) % m].ptIdx;
                    gsPoint2D<T> r1(outMesh.pts[kp].x - B.x, outMesh.pts[kp].y - B.y);
                    gsPoint2D<T> r2(outMesh.pts[kn].x - B.x, outMesh.pts[kn].y - B.y);
                    if (vecStrictlyInSector(rayC, rayA, r1) || vecStrictlyInSector(rayC, rayA, r2))
                        wedgeOk = false;
                }
                if (!wedgeOk) continue;

                bool clear = true;
                for (int k = 0; k < m && clear; ++k) {
                    int ik = activeRing[k].ptIdx;
                    if (ik == ia || ik == ib || ik == ic) continue;
                    if (pointInTri(A, B, C, outMesh.pts[ik])) clear = false;
                }
                for (int k = 0; k < m && clear; ++k) {
                    if (k == iPrev || k == idx) continue;
                    int ka = activeRing[k].ptIdx;
                    int kb = activeRing[(k + 1) % m].ptIdx;
                    // Skip the reverse copy of a bridge edge coincident with AB or BC
                    if ((ka == ia && kb == ib) || (ka == ib && kb == ia) ||
                        (ka == ib && kb == ic) || (ka == ic && kb == ib))
                        continue;
                    if (pointInTri(A, B, C, activeRing[k].q25) ||
                        pointInTri(A, B, C, activeRing[k].q50) ||
                        pointInTri(A, B, C, activeRing[k].q75)) {
                        clear = false;
                        break;
                    }
                }
                if (!clear) continue;

                T sc = T(0);
                T dAC2 = (A.x - C.x) * (A.x - C.x) + (A.y - C.y) * (A.y - C.y);
                if (earMode == 0) {
                    T l2 = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y)
                         + (B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y) + dAC2;
                    sc = area2 / (l2 + T(1e-14));
                } else if (earMode == 1) {
                    sc = T(1) / (dAC2 + T(1e-14));
                } else {
                    T lAB = std::hypot(A.x - B.x, A.y - B.y);
                    T lBC = std::hypot(B.x - C.x, B.y - C.y);
                    sc = area2 / ((lAB * lBC) * (dAC2 + T(1e-14)) + T(1e-14));
                }
                if (outMesh.isSmoothBnd[ib]) sc *= T(0.2);
                if (sc > bestScore) {
                    bestScore = sc;
                    bestIdx = idx;
                }
            }
        }

        if (bestIdx < 0) {
            if (subDivs >= std::max(48, 2 * static_cast<int>(activeRing.size()))) return false;
            int splitK = -1;
            T maxSag = -T(1);
            int longestK = -1;
            T maxLen = T(0);
            for (int k = 0; k < m; ++k) {
                int u = activeRing[k].ptIdx;
                int v = activeRing[(k + 1) % m].ptIdx;
                int eOld = activeRing[k].edgeId;
                if (!outMesh.edges[eOld].isBoundary) continue;
                const gsPoint2D<T>& mid = activeRing[k].q50;
                T chord = outMesh.pts[u].dist(outMesh.pts[v]);
                if (chord > maxLen) { maxLen = chord; longestK = k; }
                T sag = (chord > T(1e-12)) ? std::abs(cross2(outMesh.pts[u], outMesh.pts[v], mid)) / (chord * chord + T(1e-12)) : T(0);
                if (sag > maxSag) { maxSag = sag; splitK = k; }
            }
            if (splitK < 0 || maxSag < T(1e-5)) splitK = longestK;
            if (splitK < 0) return false;

            int u = activeRing[splitK].ptIdx;
            int v = activeRing[(splitK + 1) % m].ptIdx;
            int eOld = activeRing[splitK].edgeId;
            auto sp = activeRing[splitK].outSeg.splitArc(T(0.5));
            gsPoint2D<T> midPt = sp.first.p1();
            int newIdx = static_cast<int>(outMesh.pts.size());
            outMesh.pts.push_back(midPt);
            outMesh.isBoundaryPt.push_back(true);
            outMesh.isSmoothBnd.push_back(true);
            outMesh.isSharpCorner.push_back(false);
            outMesh.cornerAngleDeg.push_back(T(180));
            auto din = evalCurveDeriv(*sp.first.bsp, T(1));
            T dn = std::hypot(din.x, din.y);
            gsPoint2D<T> tMid = (dn > T(1e-14)) ? gsPoint2D<T>(din.x / dn, din.y / dn) : gsPoint2D<T>(1, 0);
            outMesh.bndTanIn.push_back(tMid);
            outMesh.bndTanOut.push_back(tMid);
            triCount.push_back(0);

            int loopIdxOld = outMesh.edges[eOld].loopIdx;
            bool fwdOld = (outMesh.edges[eOld].a == u);
            outMesh.edges[eOld] = fwdOld
                ? typename TriMeshData<T>::EdgeRec{u, newIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, loopIdxOld}
                : typename TriMeshData<T>::EdgeRec{newIdx, u, sp.first.reversed(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, loopIdxOld};
            edgeMap.erase({std::min(u, v), std::max(u, v)});
            edgeMap[{std::min(u, newIdx), std::max(u, newIdx)}] = eOld;
            int eSecond = addUndirEdge(newIdx, v, sp.second, true, loopIdxOld);

            outMesh.bndInfo.push_back({loopIdxOld, -1, T(0)});
            if (loopIdxOld >= 0 && loopIdxOld < static_cast<int>(outMesh.loopRings.size())) {
                auto& lr = outMesh.loopRings[loopIdxOld];
                for (size_t ri = 0; ri < lr.edges.size(); ++ri) {
                    if (lr.edges[ri] == eOld) {
                        lr.vertices.insert(lr.vertices.begin() + static_cast<long>(ri + 1), newIdx);
                        lr.edges.insert(lr.edges.begin() + static_cast<long>(fwdOld ? (ri + 1) : ri), eSecond);
                        break;
                    }
                }
            }

            activeRing[splitK] = makeRingNode(u, sp.first, eOld);
            activeRing.insert(activeRing.begin() + splitK + 1, makeRingNode(newIdx, sp.second, eSecond));
            ++subDivs;
            continue;
        }

        int iPrev = (bestIdx + m - 1) % m;
        int iNext = (bestIdx + 1) % m;
        int ia = activeRing[iPrev].ptIdx;
        int ib = activeRing[bestIdx].ptIdx;
        int ic = activeRing[iNext].ptIdx;

        int e0 = activeRing[iPrev].edgeId;
        int e1 = activeRing[bestIdx].edgeId;
        int e2 = getOrAddDiag(ia, ic);
        outMesh.tris.push_back({{ia, ib, ic}, {e0, e1, e2}});
        ++triCount[ia];
        ++triCount[ib];
        ++triCount[ic];

        const auto& er2 = outMesh.edges[e2];
        CurveSeg<T> diagSeg = (er2.a == ia && er2.b == ic) ? er2.seg : er2.seg.reversed();
        activeRing[iPrev] = makeRingNode(ia, diagSeg, e2);
        activeRing.erase(activeRing.begin() + bestIdx);
    }

    if (activeRing.size() != 3) return false;
    int ia = activeRing[0].ptIdx, ib = activeRing[1].ptIdx, ic = activeRing[2].ptIdx;
    if (cross2(outMesh.pts[ia], outMesh.pts[ib], outMesh.pts[ic]) <= T(1e-8)) return false;
    outMesh.tris.push_back({{ia, ib, ic}, {activeRing[0].edgeId, activeRing[1].edgeId, activeRing[2].edgeId}});
    return true;
}

// Extract exact CAD boundary curve between two boundary vertices vStart -> vEnd on loop l
template <class T>
inline CurveSeg<T> extractExactCADBoundaryBetween(
    const TriMeshData<T>& mesh,
    int vStart, int vEnd,
    int degree,
    bool forward = true)
{
    if (vStart < 0 || vEnd < 0 ||
        static_cast<size_t>(vStart) >= mesh.pts.size() ||
        static_cast<size_t>(vEnd) >= mesh.pts.size()) {
        return CurveSeg<T>();
    }
    if (vStart == vEnd) {
        return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
    }
    if (mesh.loopRings.empty()) {
        return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
    }

    int lStart = (static_cast<size_t>(vStart) < mesh.bndInfo.size()) ? mesh.bndInfo[vStart].loopIdx : -1;
    int lEnd   = (static_cast<size_t>(vEnd) < mesh.bndInfo.size()) ? mesh.bndInfo[vEnd].loopIdx : -1;
    if (lStart < 0 || lStart != lEnd || lStart >= static_cast<int>(mesh.loopRings.size())) {
        return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
    }

    const auto& ring = mesh.loopRings[lStart];
    int M = static_cast<int>(ring.vertices.size());
    if (M < 2) {
        return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
    }

    int iStart = -1, iEnd = -1;
    for (int i = 0; i < M; ++i) {
        if (ring.vertices[i] == vStart) iStart = i;
        if (ring.vertices[i] == vEnd)   iEnd = i;
    }
    if (iStart < 0 || iEnd < 0) {
        return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
    }

    std::vector<CurveSeg<T>> segs;
    if (forward) {
        int i = iStart;
        while (i != iEnd) {
            int eid = ring.edges[i];
            const auto& er = mesh.edges[eid];
            int u = ring.vertices[i];
            int v = ring.vertices[(i + 1) % M];
            CurveSeg<T> s = (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
            segs.push_back(s);
            i = (i + 1) % M;
        }
    } else {
        int i = iStart;
        while (i != iEnd) {
            int prevI = (i - 1 + M) % M;
            int eid = ring.edges[prevI];
            const auto& er = mesh.edges[eid];
            int u = ring.vertices[prevI];
            int v = ring.vertices[i];
            CurveSeg<T> s = (er.a == u && er.b == v) ? er.seg.reversed() : er.seg;
            segs.push_back(s);
            i = prevI;
        }
    }

    if (segs.empty()) {
        return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
    }
    if (segs.size() == 1) {
        return segs[0];
    }
    auto joined = joinCurveSegsExact(segs, degree);
    if (joined) {
        return CurveSeg<T>(*joined);
    }
    return makeStraightSeg(mesh.pts[vStart], mesh.pts[vEnd], degree);
}

template <class T>
inline int minRequiredBoundaryTriValence(const TriMeshData<T>& mesh, int v)
{
    if (v < 0 || v >= static_cast<int>(mesh.pts.size()) || !mesh.isBoundaryPt[v])
        return 1;
    if (mesh.cornerAngleDeg[v] > T(325))
        return 4;
    if (mesh.cornerAngleDeg[v] > T(275))
        return 3;
    if (mesh.isSmoothBnd[v] || mesh.cornerAngleDeg[v] > T(155))
        return 2;
    return 1;
}

template <class T>
inline void applyDomainLawsonFlips(TriMeshData<T>& mesh)
{
    const int nT = static_cast<int>(mesh.tris.size());
    const int nE = static_cast<int>(mesh.edges.size());
    if (nT < 2 || nE == 0) return;

    struct EdgeAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
    std::vector<EdgeAdj> eAdj(nE);
    std::vector<int> vertDeg(mesh.pts.size(), 0);
    for (int ti = 0; ti < nT; ++ti) {
        const auto& tr = mesh.tris[ti];
        for (int k = 0; k < 3; ++k) {
            ++vertDeg[tr.v[k]];
            int eid = tr.e[k];
            if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
            else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
        }
    }

    auto updateAdjSlot = [&](int eid, int oldTri, int newTri, int newK) {
        if (eAdj[eid].t0 == oldTri) { eAdj[eid].t0 = newTri; eAdj[eid].k0 = newK; }
        else if (eAdj[eid].t1 == oldTri) { eAdj[eid].t1 = newTri; eAdj[eid].k1 = newK; }
    };

    auto boundaryEdgeClearsDiagSide = [&](int eid, const gsPoint2D<T>& C, const gsPoint2D<T>& D, T reqSign) -> bool {
        const auto& er = mesh.edges[eid];
        if (!er.isBoundary || !er.seg.bsp) return true;
        for (T u : {T(0.25), T(0.5), T(0.75)}) {
            gsPoint2D<T> q = evalCurvePoint(*er.seg.bsp, u);
            if (reqSign * cross2(C, D, q) <= T(1e-8))
                return false;
        }
        return true;
    };

    std::queue<int> qEdges;
    std::vector<bool> inQ(nE, false);
    for (int eIdx = 0; eIdx < nE; ++eIdx) {
        if (!mesh.edges[eIdx].isBoundary && eAdj[eIdx].t0 >= 0 && eAdj[eIdx].t1 >= 0) {
            qEdges.push(eIdx);
            inQ[eIdx] = true;
        }
    }

    const int maxFlips = std::max(32, nT * 12);
    int flips = 0;
    while (!qEdges.empty() && flips < maxFlips) {
        int eIdx = qEdges.front();
        qEdges.pop();
        inQ[eIdx] = false;
        if (mesh.edges[eIdx].isBoundary) continue;

        int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
        int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
        if (tA < 0 || tB < 0) continue;

        const auto A = mesh.tris[tA];
        const auto B = mesh.tris[tB];
        int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
        int d = B.v[(kB + 2) % 3];
        if (c == d || c == a || c == b || d == a || d == b) continue;
        if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) continue;

        gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
        if (!isStrictlyConvexQuad(qa)) continue;
        if (!isSegmentInsideDomain(mesh.pts[c], mesh.pts[d], mesh.loopPolys)) continue;

        // Delaunay max-min-angle criterion
        T curMin0, curMax0, curMin1, curMax1;
        triMinMaxAnglesDeg(mesh.pts[a], mesh.pts[b], mesh.pts[c], curMin0, curMax0);
        triMinMaxAnglesDeg(mesh.pts[b], mesh.pts[a], mesh.pts[d], curMin1, curMax1);
        T newMin0, newMax0, newMin1, newMax1;
        triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[a], mesh.pts[d], newMin0, newMax0);
        triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[d], mesh.pts[b], newMin1, newMax1);

        T curMin = std::min(curMin0, curMin1);
        T newMin = std::min(newMin0, newMin1);
        T curMax = std::max(curMax0, curMax1);
        T newMax = std::max(newMax0, newMax1);
        if (newMin <= curMin + T(1e-3) && !(std::abs(newMin - curMin) <= T(1e-3) && newMax < curMax - T(1e-2)))
            continue;

        if (vertDeg[a] <= minRequiredBoundaryTriValence(mesh, a) ||
            vertDeg[b] <= minRequiredBoundaryTriValence(mesh, b))
            continue;

        int e_bc = A.e[(kA + 1) % 3];
        int e_ca = A.e[(kA + 2) % 3];
        int e_ad = B.e[(kB + 1) % 3];
        int e_db = B.e[(kB + 2) % 3];
        if (!boundaryEdgeClearsDiagSide(e_ca, mesh.pts[c], mesh.pts[d], T(-1)) ||
            !boundaryEdgeClearsDiagSide(e_ad, mesh.pts[c], mesh.pts[d], T(-1)) ||
            !boundaryEdgeClearsDiagSide(e_db, mesh.pts[c], mesh.pts[d], T(1)) ||
            !boundaryEdgeClearsDiagSide(e_bc, mesh.pts[c], mesh.pts[d], T(1)))
            continue;

        mesh.edges[eIdx] = {c, d, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1};
        mesh.tris[tA] = {{c, a, d}, {e_ca, e_ad, eIdx}};
        mesh.tris[tB] = {{c, d, b}, {eIdx, e_db, e_bc}};

        eAdj[eIdx] = {tA, 2, tB, 0};
        updateAdjSlot(e_ca, tA, tA, 0);
        updateAdjSlot(e_ad, tB, tA, 1);
        updateAdjSlot(e_db, tB, tB, 1);
        updateAdjSlot(e_bc, tA, tB, 2);

        --vertDeg[a]; --vertDeg[b];
        ++vertDeg[c]; ++vertDeg[d];
        ++flips;

        for (int ne : {e_ca, e_ad, e_db, e_bc}) {
            if (!mesh.edges[ne].isBoundary && !inQ[ne] && eAdj[ne].t0 >= 0 && eAdj[ne].t1 >= 0) {
                qEdges.push(ne);
                inQ[ne] = true;
            }
        }
    }
}

template <class T>
inline bool curvedBoundaryEdgeClearsVertex(const TriMeshData<T>& mesh, int eid, int u, int v, const gsPoint2D<T>& opp)
{
    const auto& er = mesh.edges[eid];
    if (!er.isBoundary || !er.seg.bsp) return true;
    const gsPoint2D<T>& pu = mesh.pts[u];
    const gsPoint2D<T>& pv = mesh.pts[v];
    if (cross2(pu, pv, opp) <= T(1e-7)) return false;
    for (T s : {T(0.25), T(0.5), T(0.75)}) {
        gsPoint2D<T> q = evalCurvePoint(*er.seg.bsp, s);
        if (cross2(pv, opp, q) <= T(1e-8) || cross2(opp, pu, q) <= T(1e-8))
            return false;
    }
    return true;
}

template <class T>
inline void smoothInteriorVertices(TriMeshData<T>& mesh, int sweeps = 4)
{
    const int nV = static_cast<int>(mesh.pts.size());
    std::vector<std::vector<int>> adjTris(nV);
    std::vector<std::vector<int>> adjVerts(nV);
    for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
        const auto& t = mesh.tris[ti];
        for (int k = 0; k < 3; ++k) {
            adjTris[t.v[k]].push_back(ti);
            adjVerts[t.v[k]].push_back(t.v[(k + 1) % 3]);
            adjVerts[t.v[k]].push_back(t.v[(k + 2) % 3]);
        }
    }
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (int v = 0; v < nV; ++v) {
            if (mesh.isBoundaryPt[v] || adjVerts[v].empty()) continue;
            gsPoint2D<T> avg(0, 0);
            for (int nb : adjVerts[v]) { avg.x += mesh.pts[nb].x; avg.y += mesh.pts[nb].y; }
            avg.x /= static_cast<T>(adjVerts[v].size());
            avg.y /= static_cast<T>(adjVerts[v].size());
            if (!isPointStrictlyInsideDomain(avg, mesh.loopPolys, T(1e-5))) continue;

            T oldMinAng = T(180);
            for (int ti : adjTris[v]) {
                const auto& t = mesh.tris[ti];
                T mn = T(0), mx = T(180);
                triMinMaxAnglesDeg(mesh.pts[t.v[0]], mesh.pts[t.v[1]], mesh.pts[t.v[2]], mn, mx);
                oldMinAng = std::min(oldMinAng, mn);
            }
            const T reqMinAng = std::min(T(8.0), std::max(T(1.5), oldMinAng - T(0.2)));

            gsPoint2D<T> old = mesh.pts[v];
            mesh.pts[v] = avg;
            bool ok = true;
            for (int ti : adjTris[v]) {
                const auto& t = mesh.tris[ti];
                if (cross2(mesh.pts[t.v[0]], mesh.pts[t.v[1]], mesh.pts[t.v[2]]) <= T(1e-7)) {
                    ok = false; break;
                }
                T mn = T(0), mx = T(180);
                triMinMaxAnglesDeg(mesh.pts[t.v[0]], mesh.pts[t.v[1]], mesh.pts[t.v[2]], mn, mx);
                if (mn < reqMinAng) { ok = false; break; }
                for (int k = 0; k < 3; ++k) {
                    if (t.v[(k + 2) % 3] == v) {
                        int eid = t.e[k];
                        if (mesh.edges[eid].isBoundary) {
                            if (!curvedBoundaryEdgeClearsVertex(mesh, eid, t.v[k], t.v[(k + 1) % 3], avg) ||
                                !isSegmentInsideDomain(mesh.pts[t.v[k]], avg, mesh.loopPolys) ||
                                !isSegmentInsideDomain(mesh.pts[t.v[(k + 1) % 3]], avg, mesh.loopPolys)) {
                                ok = false; break;
                            }
                        }
                    }
                }
                if (!ok) break;
            }
            if (!ok) mesh.pts[v] = old;
        }
    }
    applyDomainLawsonFlips(mesh);
}

template <class T>
inline void refineInteriorTriangles(TriMeshData<T>& mesh, T targetH, int maxRounds = 3)
{
    const T areaThresh = T(0.65) * targetH * targetH;
    for (int round = 0; round < maxRounds; ++round) {
        int nOrig = static_cast<int>(mesh.tris.size());
        bool added = false;
        for (int ti = 0; ti < nOrig; ++ti) {
            const auto t = mesh.tris[ti];
            int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
            T a2 = cross2(mesh.pts[v0], mesh.pts[v1], mesh.pts[v2]);
            if (a2 * T(0.5) <= areaThresh) continue;

            gsPoint2D<T> c((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x) / T(3),
                           (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y) / T(3));
            if (!isPointStrictlyInsideDomain(c, mesh.loopPolys, T(1e-5))) continue;
            if (cross2(mesh.pts[v0], mesh.pts[v1], c) <= T(1e-7) ||
                cross2(mesh.pts[v1], mesh.pts[v2], c) <= T(1e-7) ||
                cross2(mesh.pts[v2], mesh.pts[v0], c) <= T(1e-7))
                continue;
            T mn0, mx0, mn1, mx1, mn2, mx2;
            triMinMaxAnglesDeg(mesh.pts[v0], mesh.pts[v1], c, mn0, mx0);
            triMinMaxAnglesDeg(mesh.pts[v1], mesh.pts[v2], c, mn1, mx1);
            triMinMaxAnglesDeg(mesh.pts[v2], mesh.pts[v0], c, mn2, mx2);
            if (std::min(mn0, std::min(mn1, mn2)) < T(4.0)) continue;
            if (!isSegmentInsideDomain(mesh.pts[v0], c, mesh.loopPolys) ||
                !isSegmentInsideDomain(mesh.pts[v1], c, mesh.loopPolys) ||
                !isSegmentInsideDomain(mesh.pts[v2], c, mesh.loopPolys))
                continue;
            if (!curvedBoundaryEdgeClearsVertex(mesh, t.e[0], v0, v1, c) ||
                !curvedBoundaryEdgeClearsVertex(mesh, t.e[1], v1, v2, c) ||
                !curvedBoundaryEdgeClearsVertex(mesh, t.e[2], v2, v0, c))
                continue;

            int vc = static_cast<int>(mesh.pts.size());
            mesh.pts.push_back(c);
            mesh.isBoundaryPt.push_back(false);
            mesh.isSmoothBnd.push_back(false);
            mesh.isSharpCorner.push_back(false);
            mesh.cornerAngleDeg.push_back(T(360));
            mesh.bndTanIn.push_back(gsPoint2D<T>(1, 0));
            mesh.bndTanOut.push_back(gsPoint2D<T>(1, 0));

            int e_0c = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({v0, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
            int e_1c = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({v1, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
            int e_2c = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({v2, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});

            mesh.tris[ti] = {{v0, v1, vc}, {t.e[0], e_1c, e_0c}};
            mesh.tris.push_back({{v1, v2, vc}, {t.e[1], e_2c, e_1c}});
            mesh.tris.push_back({{v2, v0, vc}, {t.e[2], e_0c, e_2c}});
            added = true;
        }
        if (!added) break;
        applyDomainLawsonFlips(mesh);
    }

    // Safe Laplacian smoothing of interior Steiner vertices
    smoothInteriorVertices(mesh, 4);
}

// Rivara longest-edge bisection for triangles whose 1-to-3 quadification still has angles > maxQuadAngleDeg
template <class T>
inline void refineNeedleCircumcenters(TriMeshData<T>& mesh, T targetH, int maxRounds = 2, T minTriAngThresh = T(16), bool allowCentroidFallback = false, int maxNewPts = 36)
{
    const T radToDeg = T(180) / T(3.14159265358979323846);
    int totalAdded = 0;
    for (int round = 0; round < maxRounds && totalAdded < maxNewPts; ++round) {
        int nOrig = static_cast<int>(mesh.tris.size());
        bool added = false;
        for (int ti = 0; ti < nOrig && totalAdded < maxNewPts; ++ti) {
            if (mesh.tris.size() > 420) break;
            const auto t = mesh.tris[ti];
            int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
            const auto& P0 = mesh.pts[v0];
            const auto& P1 = mesh.pts[v1];
            const auto& P2 = mesh.pts[v2];
            T a2 = cross2(P0, P1, P2);
            if (a2 <= T(1e-10)) continue;

            T l0 = std::hypot(P1.x - P0.x, P1.y - P0.y);
            T l1 = std::hypot(P2.x - P1.x, P2.y - P1.y);
            T l2 = std::hypot(P0.x - P2.x, P0.y - P2.y);
            T lMin = std::min(l0, std::min(l1, l2));
            T lMax = std::max(l0, std::max(l1, l2));
            if (lMax <= T(2.3) * lMin || lMax <= T(0.28) * targetH) continue;

            auto angAt = [&](const gsPoint2D<T>& Q, const gsPoint2D<T>& U, const gsPoint2D<T>& V, T lu, T lv) -> T {
                if (lu < T(1e-14) || lv < T(1e-14)) return T(0);
                T ct = ((U.x - Q.x) * (V.x - Q.x) + (U.y - Q.y) * (V.y - Q.y)) / (lu * lv);
                return std::acos(std::max(T(-1), std::min(T(1), ct))) * radToDeg;
            };
            T angs[3] = { angAt(P0, P1, P2, l0, l2), angAt(P1, P2, P0, l1, l0), angAt(P2, P0, P1, l2, l1) };
            int kMin = 0;
            if (angs[1] < angs[kMin]) kMin = 1;
            if (angs[2] < angs[kMin]) kMin = 2;
            int vMinPt = t.v[kMin];
            bool atSharpCorner = mesh.isSharpCorner[vMinPt] &&
                                 mesh.edges[t.e[kMin]].isBoundary &&
                                 mesh.edges[t.e[(kMin + 2) % 3]].isBoundary;
            if (atSharpCorner || angs[kMin] >= minTriAngThresh) continue;

            gsPoint2D<T> G((P0.x + P1.x + P2.x) / T(3), (P0.y + P1.y + P2.y) / T(3));
            T lsum = l0 + l1 + l2; // opposite: l1->P0, l2->P1, l0->P2
            gsPoint2D<T> O_inc = (lsum > T(1e-14))
                ? gsPoint2D<T>((l1 * P0.x + l2 * P1.x + l0 * P2.x) / lsum,
                               (l1 * P0.y + l2 * P1.y + l0 * P2.y) / lsum)
                : G;

            std::vector<gsPoint2D<T>> steinerCands;
            T ax = P1.x - P0.x, ay = P1.y - P0.y;
            T bx = P2.x - P0.x, by = P2.y - P0.y;
            T D = T(2) * (ax * by - ay * bx);
            bool ccValid = false;
            gsPoint2D<T> ccRaw = G;
            if (std::abs(D) > T(1e-14)) {
                T aSq = ax * ax + ay * ay, bSq = bx * bx + by * by;
                ccRaw = gsPoint2D<T>(P0.x + (by * aSq - ay * bSq) / D,
                                     P0.y + (ax * bSq - bx * aSq) / D);
                T h0 = cross2(P0, P1, ccRaw) / (l0 + T(1e-14));
                T h1 = cross2(P1, P2, ccRaw) / (l1 + T(1e-14));
                T h2 = cross2(P2, P0, ccRaw) / (l2 + T(1e-14));
                if (std::min(h0, std::min(h1, h2)) >= T(0.25) * lMin) {
                    steinerCands.push_back(ccRaw);
                    ccValid = true;
                }
            }
            if (!ccValid && !allowCentroidFallback) continue;
            if (!ccValid) {
                // Ungor off-center blends along G->ccRaw and O_inc->G staying strictly inside triangle
                steinerCands.push_back(O_inc);
                steinerCands.push_back(G);
                for (T w : {T(0.25), T(0.5), T(0.75)}) {
                    steinerCands.emplace_back((T(1) - w) * O_inc.x + w * ccRaw.x,
                                              (T(1) - w) * O_inc.y + w * ccRaw.y);
                    steinerCands.emplace_back((T(1) - w) * G.x + w * ccRaw.x,
                                              (T(1) - w) * G.y + w * ccRaw.y);
                }
            }

            gsPoint2D<T> bestPt = G;
            T bestMinChildAng = -T(1);
            bool foundPt = false;
            for (const auto& cand : steinerCands) {
                T h0 = cross2(P0, P1, cand) / (l0 + T(1e-14));
                T h1 = cross2(P1, P2, cand) / (l1 + T(1e-14));
                T h2 = cross2(P2, P0, cand) / (l2 + T(1e-14));
                if (cross2(P0, P1, cand) <= T(1e-7) || cross2(P1, P2, cand) <= T(1e-7) || cross2(P2, P0, cand) <= T(1e-7)) continue;
                if (std::min(h0, std::min(h1, h2)) < T(0.16) * lMin) continue;
                if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-5))) continue;
                if (!isSegmentInsideDomain(P0, cand, mesh.loopPolys) ||
                    !isSegmentInsideDomain(P1, cand, mesh.loopPolys) ||
                    !isSegmentInsideDomain(P2, cand, mesh.loopPolys))
                    continue;
                if (!curvedBoundaryEdgeClearsVertex(mesh, t.e[0], v0, v1, cand) ||
                    !curvedBoundaryEdgeClearsVertex(mesh, t.e[1], v1, v2, cand) ||
                    !curvedBoundaryEdgeClearsVertex(mesh, t.e[2], v2, v0, cand))
                    continue;
                T mn0, mx0, mn1, mx1, mn2, mx2;
                triMinMaxAnglesDeg(P0, P1, cand, mn0, mx0);
                triMinMaxAnglesDeg(P1, P2, cand, mn1, mx1);
                triMinMaxAnglesDeg(P2, P0, cand, mn2, mx2);
                T mnAll = std::min(mn0, std::min(mn1, mn2));
                if (!ccValid && mnAll <= angs[kMin] + T(2.0)) continue;
                if (!foundPt || mnAll > bestMinChildAng) {
                    bestMinChildAng = mnAll;
                    bestPt = cand;
                    foundPt = true;
                    if (ccValid) break;
                }
            }
            if (!foundPt) continue;

            int vc = static_cast<int>(mesh.pts.size());
            mesh.pts.push_back(bestPt);
            mesh.isBoundaryPt.push_back(false);
            mesh.isSmoothBnd.push_back(false);
            mesh.isSharpCorner.push_back(false);
            mesh.cornerAngleDeg.push_back(T(360));
            mesh.bndTanIn.push_back(gsPoint2D<T>(1, 0));
            mesh.bndTanOut.push_back(gsPoint2D<T>(1, 0));

            int e_0c = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({v0, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
            int e_1c = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({v1, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
            int e_2c = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({v2, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});

            mesh.tris[ti] = {{v0, v1, vc}, {t.e[0], e_1c, e_0c}};
            mesh.tris.push_back({{v1, v2, vc}, {t.e[1], e_2c, e_1c}});
            mesh.tris.push_back({{v2, v0, vc}, {t.e[2], e_0c, e_2c}});
            added = true;
            ++totalAdded;
        }
        if (!added) break;
        applyDomainLawsonFlips(mesh);
    }
    smoothInteriorVertices(mesh, 4);
}

// Shared helper that tries splitting one chosen edge of triangle `badTi` (and its mate if interior).
template <class T>
inline bool trySplitTriangleEdges(
    TriMeshData<T>& mesh,
    int badTi,
    const int order[3],
    int degree,
    T minAllowedChildAng)
{
    const auto bt = mesh.tris[badTi];
    auto segClears = [&](const CurveSeg<T>& seg, const gsPoint2D<T>& pStart,
                         const gsPoint2D<T>& pEnd, const gsPoint2D<T>& qOpp) -> bool {
        if (!seg.bsp) return true;
        for (T s : {T(0.25), T(0.5), T(0.75)}) {
            gsPoint2D<T> q = evalCurvePoint(*seg.bsp, s);
            if (cross2(pEnd, qOpp, q) <= T(1e-8) || cross2(qOpp, pStart, q) <= T(1e-8))
                return false;
        }
        return true;
    };

    for (int oi = 0; oi < 3; ++oi) {
        int kMax = order[oi];
        int u = bt.v[kMax], v = bt.v[(kMax + 1) % 3], opp = bt.v[(kMax + 2) % 3];
        int eIdx = bt.e[kMax];
        int e_vo = bt.e[(kMax + 1) % 3];
        int e_ou = bt.e[(kMax + 2) % 3];

        auto oldEr = mesh.edges[eIdx];
        if (!oldEr.seg.bsp)
            oldEr.seg = makeStraightSeg(mesh.pts[oldEr.a], mesh.pts[oldEr.b], degree);

        int mateTi = -1, mateK = -1;
        for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
            if (ti == badTi) continue;
            for (int k = 0; k < 3; ++k) {
                if (mesh.tris[ti].e[k] == eIdx) { mateTi = ti; mateK = k; break; }
            }
            if (mateTi >= 0) break;
        }
        int mOpp = (mateTi >= 0) ? mesh.tris[mateTi].v[(mateK + 2) % 3] : -1;

        T tMid = arcLengthSplitParam(*oldEr.seg.bsp, T(0.5));
        std::vector<T> tCands = {tMid};
        T tOpp = std::max(T(0.20), std::min(T(0.80), oldEr.seg.closestParam(mesh.pts[opp])));
        if (std::abs(tOpp - tMid) > T(0.03)) tCands.push_back(tOpp);
        if (mOpp >= 0) {
            T tMate = std::max(T(0.20), std::min(T(0.80), oldEr.seg.closestParam(mesh.pts[mOpp])));
            T tAvg = T(0.5) * (tOpp + tMate);
            if (std::abs(tMate - tMid) > T(0.03)) tCands.push_back(tMate);
            if (std::abs(tAvg - tMid) > T(0.03))  tCands.push_back(tAvg);
        }

        std::pair<CurveSeg<T>, CurveSeg<T>> bestSp;
        T bestSplitScore = -T(1e9);
        bool foundSplit = false;
        for (size_t ci = 0; ci < tCands.size(); ++ci) {
            auto spTry = oldEr.seg.splitAtParam(tCands[ci]);
            gsPoint2D<T> sPt = spTry.first.p1();
            if (cross2(mesh.pts[u], sPt, mesh.pts[opp]) <= T(1e-7) ||
                cross2(sPt, mesh.pts[v], mesh.pts[opp]) <= T(1e-7))
                continue;
            if (!curvedBoundaryEdgeClearsVertex(mesh, e_ou, opp, u, sPt) ||
                !curvedBoundaryEdgeClearsVertex(mesh, e_vo, v, opp, sPt))
                continue;
            if (oldEr.isBoundary) {
                if (!isSegmentInsideDomain(mesh.pts[opp], sPt, mesh.loopPolys)) continue;
                const bool fwdUV = (oldEr.a == u);
                const auto& segUS = fwdUV ? spTry.first : spTry.second.reversed();
                const auto& segSV = fwdUV ? spTry.second : spTry.first.reversed();
                if (!segClears(segUS, mesh.pts[u], sPt, mesh.pts[opp]) ||
                    !segClears(segSV, sPt, mesh.pts[v], mesh.pts[opp]))
                    continue;
            }
            if (mateTi >= 0) {
                int mu = mesh.tris[mateTi].v[mateK], mv = mesh.tris[mateTi].v[(mateK + 1) % 3];
                int e_mvo = mesh.tris[mateTi].e[(mateK + 1) % 3];
                int e_mou = mesh.tris[mateTi].e[(mateK + 2) % 3];
                if (cross2(mesh.pts[mu], sPt, mesh.pts[mOpp]) <= T(1e-7) ||
                    cross2(sPt, mesh.pts[mv], mesh.pts[mOpp]) <= T(1e-7))
                    continue;
                if (!curvedBoundaryEdgeClearsVertex(mesh, e_mou, mOpp, mu, sPt) ||
                    !curvedBoundaryEdgeClearsVertex(mesh, e_mvo, mv, mOpp, sPt))
                    continue;
                if (oldEr.isBoundary) {
                    if (!isSegmentInsideDomain(mesh.pts[mOpp], sPt, mesh.loopPolys)) continue;
                }
            }
            T mnA, mxA, mnB, mxB;
            triMinMaxAnglesDeg(mesh.pts[u], sPt, mesh.pts[opp], mnA, mxA);
            triMinMaxAnglesDeg(sPt, mesh.pts[v], mesh.pts[opp], mnB, mxB);
            T mnAll = std::min(mnA, mnB);
            if (mateTi >= 0) {
                int mu = mesh.tris[mateTi].v[mateK], mv = mesh.tris[mateTi].v[(mateK + 1) % 3];
                T mnC, mxC, mnD, mxD;
                triMinMaxAnglesDeg(mesh.pts[mu], sPt, mesh.pts[mOpp], mnC, mxC);
                triMinMaxAnglesDeg(sPt, mesh.pts[mv], mesh.pts[mOpp], mnD, mxD);
                mnAll = std::min(mnAll, std::min(mnC, mnD));
            }
            if (mnAll < minAllowedChildAng) continue;
            T sc = mnAll + (ci == 0 ? T(1.5) : T(0));
            if (!foundSplit || sc > bestSplitScore) {
                bestSplitScore = sc;
                bestSp = spTry;
                foundSplit = true;
            }
        }
        if (!foundSplit) continue;
        auto& sp = bestSp;
        gsPoint2D<T> midPt = sp.first.p1();

        int mIdx = static_cast<int>(mesh.pts.size());
        mesh.pts.push_back(midPt);
        mesh.isBoundaryPt.push_back(oldEr.isBoundary);
        mesh.isSmoothBnd.push_back(oldEr.isBoundary);
        mesh.isSharpCorner.push_back(false);
        mesh.cornerAngleDeg.push_back(oldEr.isBoundary ? T(180) : T(360));
        auto din = evalCurveDeriv(*sp.first.bsp, T(1));
        T dn = std::hypot(din.x, din.y);
        gsPoint2D<T> tMidDir = (dn > T(1e-14)) ? gsPoint2D<T>(din.x / dn, din.y / dn) : gsPoint2D<T>(1, 0);
        mesh.bndTanIn.push_back(tMidDir);
        mesh.bndTanOut.push_back(tMidDir);

        int oldA = oldEr.a, oldB = oldEr.b;
        mesh.edges[eIdx] = {oldA, mIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, oldEr.isBoundary, oldEr.loopIdx};
        int e_mb = static_cast<int>(mesh.edges.size());
        mesh.edges.push_back({mIdx, oldB, sp.second, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, oldEr.isBoundary, oldEr.loopIdx});

        auto halfFor = [&](int endPt) { return (endPt == oldA) ? eIdx : e_mb; };
        int e_om = static_cast<int>(mesh.edges.size());
        mesh.edges.push_back({opp, mIdx, makeStraightSeg(mesh.pts[opp], midPt, degree),
                              CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
        mesh.tris[badTi] = {{u, mIdx, opp}, {halfFor(u), e_om, e_ou}};
        mesh.tris.push_back({{mIdx, v, opp}, {halfFor(v), e_vo, e_om}});

        if (mateTi >= 0) {
            const auto mt = mesh.tris[mateTi];
            int mu = mt.v[mateK], mv = mt.v[(mateK + 1) % 3];
            int e_mvo = mt.e[(mateK + 1) % 3];
            int e_mou = mt.e[(mateK + 2) % 3];
            int e_mom = static_cast<int>(mesh.edges.size());
            mesh.edges.push_back({mOpp, mIdx, makeStraightSeg(mesh.pts[mOpp], midPt, degree),
                                  CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
            mesh.tris[mateTi] = {{mu, mIdx, mOpp}, {halfFor(mu), e_mom, e_mou}};
            mesh.tris.push_back({{mIdx, mv, mOpp}, {halfFor(mv), e_mvo, e_mom}});
        }

        applyDomainLawsonFlips(mesh);
        return true;
    }
    return false;
}

template <class T>
inline void fixUndervalentAndReflexCorners(TriMeshData<T>& mesh, int degree, int maxSplits = 12)
{
    auto edgeDirAt = [&](int ei, int fromV, int toV, bool atEnd) -> gsPoint2D<T> {
        const auto& er = mesh.edges[ei];
        if (er.isBoundary && er.seg.bsp) {
            bool fwd = (er.a == fromV);
            T u = (atEnd == fwd) ? T(1) : T(0);
            gsPoint2D<T> dv = evalCurveDeriv(*er.seg.bsp, u);
            T sgn = ((fwd && !atEnd) || (!fwd && atEnd)) ? T(1) : T(-1);
            T nm = std::hypot(dv.x, dv.y);
            if (nm > T(1e-14)) return gsPoint2D<T>(sgn * dv.x / nm, sgn * dv.y / nm);
        }
        T chSign = atEnd ? T(-1) : T(1);
        T dx = chSign * (mesh.pts[toV].x - mesh.pts[fromV].x);
        T dy = chSign * (mesh.pts[toV].y - mesh.pts[fromV].y);
        T nm = std::hypot(dx, dy);
        return (nm > T(1e-14)) ? gsPoint2D<T>(dx / nm, dy / nm) : gsPoint2D<T>(1, 0);
    };
    const T radToDeg = T(180) / T(3.14159265358979323846);

    for (int iter = 0; iter < maxSplits; ++iter) {
        if (mesh.tris.size() > 500) break;
        std::vector<int> vertDeg(mesh.pts.size(), 0);
        for (const auto& tr : mesh.tris) {
            ++vertDeg[tr.v[0]]; ++vertDeg[tr.v[1]]; ++vertDeg[tr.v[2]];
        }
        struct Cand { int ti, cornerIdx; T score; };
        std::vector<Cand> cands;
        for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
            const auto& t = mesh.tris[ti];
            for (int idx = 0; idx < 3; ++idx) {
                int vk = t.v[idx];
                if (vertDeg[vk] < minRequiredBoundaryTriValence(mesh, vk)) {
                    cands.push_back({ti, idx, T(200) + mesh.cornerAngleDeg[vk]});
                    continue;
                }
                int eOut = t.e[idx], eIn = t.e[(idx + 2) % 3];
                if (mesh.isBoundaryPt[vk] && (mesh.edges[eOut].isBoundary != mesh.edges[eIn].isBoundary)) {
                    int vPrev = t.v[(idx + 2) % 3], vNext = t.v[(idx + 1) % 3];
                    gsPoint2D<T> tOut = edgeDirAt(eOut, vk, vNext, false);
                    gsPoint2D<T> tBack = edgeDirAt(eIn, vPrev, vk, true);
                    T cAng = std::atan2(tOut.x * tBack.y - tOut.y * tBack.x,
                                        tOut.x * tBack.x + tOut.y * tBack.y) * radToDeg;
                    if (cAng <= T(0)) cAng += T(360);
                    if (cAng > T(160) && cAng < T(215))
                        cands.push_back({ti, idx, cAng});
                }
            }
        }
        if (cands.empty()) break;
        std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.score > b.score; });

        bool didSplit = false;
        for (const auto& c : cands) {
            int order[3] = { (c.cornerIdx + 1) % 3, c.cornerIdx, (c.cornerIdx + 2) % 3 };
            if (trySplitTriangleEdges(mesh, c.ti, order, degree, T(1.5))) {
                didSplit = true;
                break;
            }
        }
        if (!didSplit) break;
    }
}

template <class T>
inline void bisectObtuseTriangles(TriMeshData<T>& mesh, int degree, T maxQuadAngleDeg = T(146), int maxSplits = 28, T minQuadAngleDeg = T(14.5))
{
    auto edgeDirAt = [&](int ei, int fromV, int toV, bool atEnd) -> gsPoint2D<T> {
        const auto& er = mesh.edges[ei];
        if (er.isBoundary && er.seg.bsp) {
            bool fwd = (er.a == fromV);
            T u = (atEnd == fwd) ? T(1) : T(0);
            gsPoint2D<T> dv = evalCurveDeriv(*er.seg.bsp, u);
            T sgn = ((fwd && !atEnd) || (!fwd && atEnd)) ? T(1) : T(-1);
            T nm = std::hypot(dv.x, dv.y);
            if (nm > T(1e-14)) return gsPoint2D<T>(sgn * dv.x / nm, sgn * dv.y / nm);
        }
        T chSign = atEnd ? T(-1) : T(1);
        T dx = chSign * (mesh.pts[toV].x - mesh.pts[fromV].x);
        T dy = chSign * (mesh.pts[toV].y - mesh.pts[fromV].y);
        T nm = std::hypot(dx, dy);
        return (nm > T(1e-14)) ? gsPoint2D<T>(dx / nm, dy / nm) : gsPoint2D<T>(1, 0);
    };

    auto triQuadMaxAngle = [&](const typename TriMeshData<T>::Tri& t, const std::vector<int>& vertDeg, int* outLowValIdx = nullptr) -> T {
        if (outLowValIdx) *outLowValIdx = -1;
        T forcedScore = T(-1);
        for (int idx = 0; idx < 3; ++idx) {
            int vk = t.v[idx];
            if (vertDeg[vk] < minRequiredBoundaryTriValence(mesh, vk)) {
                T sc = T(200) + mesh.cornerAngleDeg[vk];
                if (sc > forcedScore) {
                    forcedScore = sc;
                    if (outLowValIdx) *outLowValIdx = idx;
                }
            }
        }
        if (forcedScore > T(0)) return forcedScore;

        const auto& A = mesh.pts[t.v[0]];
        const auto& B = mesh.pts[t.v[1]];
        const auto& C = mesh.pts[t.v[2]];
        auto midOf = [&](int eIdx, int oppV) -> gsPoint2D<T> {
            const auto& er = mesh.edges[eIdx];
            if (er.seg.bsp) {
                if (er.isBoundary && oppV >= 0) {
                    T uPrj = std::max(T(0.14), std::min(T(0.86), er.seg.closestParam(mesh.pts[oppV])));
                    T uMid = er.seg.splitArc(T(0.5)).first.p1().dist(mesh.pts[oppV]);
                    gsPoint2D<T> pPrj = evalCurvePoint(*er.seg.bsp, uPrj);
                    if (pPrj.dist(mesh.pts[oppV]) < uMid * T(0.95))
                        return pPrj;
                }
                return er.seg.splitArc(T(0.5)).first.p1();
            }
            return gsPoint2D<T>(T(0.5) * (mesh.pts[er.a].x + mesh.pts[er.b].x),
                                T(0.5) * (mesh.pts[er.a].y + mesh.pts[er.b].y));
        };
        gsPoint2D<T> m01 = midOf(t.e[0], t.v[2]);
        gsPoint2D<T> m12 = midOf(t.e[1], t.v[0]);
        gsPoint2D<T> m20 = midOf(t.e[2], t.v[1]);
        gsPoint2D<T> O;
        if (!pickOptimalTri1to3Center(A, B, C, m01, m12, m20, O)) return T(180);

        auto quadScore = [&](const gsPoint2D<T> q[4], bool cornerSharp) -> T {
            const T radToDeg = T(180) / T(3.14159265358979323846);
            T mx = T(0), mn = T(180);
            for (int k = 0; k < 4; ++k) {
                const auto& p0 = q[(k + 3) % 4];
                const auto& pk = q[k];
                const auto& p1 = q[(k + 1) % 4];
                T ux = p1.x - pk.x, uy = p1.y - pk.y;
                T vx = p0.x - pk.x, vy = p0.y - pk.y;
                T nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
                if (nu <= T(1e-14) || nv <= T(1e-14)) return T(180);
                T ang = std::atan2((ux * vy - uy * vx) / (nu * nv),
                                   (ux * vx + uy * vy) / (nu * nv)) * radToDeg;
                if (ang <= T(0)) ang += T(360);
                if (ang > mx) mx = ang;
                if (!(k == 0 && cornerSharp) && ang < mn) mn = ang;
            }
            if (mn < minQuadAngleDeg)
                return std::max(mx, maxQuadAngleDeg + (minQuadAngleDeg - mn) * T(3));
            return mx;
        };
        auto isCornerSharp = [&](int idxInTri) -> bool {
            int vk = t.v[idxInTri];
            return mesh.isSharpCorner[vk] &&
                   mesh.edges[t.e[idxInTri]].isBoundary &&
                   mesh.edges[t.e[(idxInTri + 2) % 3]].isBoundary;
        };
        gsPoint2D<T> q0[4] = {A, m01, O, m20};
        gsPoint2D<T> q1[4] = {B, m12, O, m01};
        gsPoint2D<T> q2[4] = {C, m20, O, m12};
        T bestAng = std::max(quadScore(q0, isCornerSharp(0)),
                             std::max(quadScore(q1, isCornerSharp(1)), quadScore(q2, isCornerSharp(2))));

        const T radToDeg = T(180) / T(3.14159265358979323846);
        for (int idx = 0; idx < 3; ++idx) {
            int eOut = t.e[idx], eIn = t.e[(idx + 2) % 3];
            if (!mesh.edges[eOut].isBoundary && !mesh.edges[eIn].isBoundary) continue;
            if (mesh.edges[eOut].isBoundary && mesh.edges[eIn].isBoundary) continue;
            int vPrev = t.v[(idx + 2) % 3], vk = t.v[idx], vNext = t.v[(idx + 1) % 3];
            gsPoint2D<T> tOut = edgeDirAt(eOut, vk, vNext, false);
            gsPoint2D<T> tBack = edgeDirAt(eIn, vPrev, vk, true);
            T cr = tOut.x * tBack.y - tOut.y * tBack.x;
            T dt = tOut.x * tBack.x + tOut.y * tBack.y;
            T cAng = std::atan2(cr, dt) * radToDeg;
            if (cAng <= T(0)) cAng += T(360);
            if (cAng > T(158) && cAng < T(215) && cAng > bestAng) bestAng = cAng;
        }
        return bestAng;
    };

    for (int iter = 0; iter < maxSplits; ++iter) {
        if (mesh.tris.size() > 500) break;
        std::vector<int> vertDeg(mesh.pts.size(), 0);
        for (const auto& tr : mesh.tris) {
            ++vertDeg[tr.v[0]]; ++vertDeg[tr.v[1]]; ++vertDeg[tr.v[2]];
        }
        struct BadCand { int ti, lvi; T ang; };
        std::vector<BadCand> cands;
        for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
            int lvi = -1;
            T a = triQuadMaxAngle(mesh.tris[ti], vertDeg, &lvi);
            if (a > maxQuadAngleDeg)
                cands.push_back({ti, lvi, a});
        }
        if (cands.empty()) break;
        std::sort(cands.begin(), cands.end(), [](const BadCand& x, const BadCand& y) { return x.ang > y.ang; });

        bool didSplit = false;
        for (const auto& cand : cands) {
            const auto bt = mesh.tris[cand.ti];
            int order[3] = {0, 1, 2};
            T lens[3];
            for (int k = 0; k < 3; ++k)
                lens[k] = mesh.pts[bt.v[k]].dist(mesh.pts[bt.v[(k + 1) % 3]]);
            std::sort(order, order + 3, [&](int a, int b) { return lens[a] > lens[b]; });

            T oldMn0, oldMx0;
            triMinMaxAnglesDeg(mesh.pts[bt.v[0]], mesh.pts[bt.v[1]], mesh.pts[bt.v[2]], oldMn0, oldMx0);
            T minAllowedChildAng = std::max(T(3.5), std::min(T(10.0), T(0.65) * oldMn0));

            int forceIdx = cand.lvi;
            if (forceIdx < 0 && oldMx0 > T(150)) {
                // Put edge opposite the widest vertex angle first
                T bestCornerA = -T(1);
                for (int idx = 0; idx < 3; ++idx) {
                    const auto& Q = mesh.pts[bt.v[idx]];
                    const auto& U = mesh.pts[bt.v[(idx + 1) % 3]];
                    const auto& V = mesh.pts[bt.v[(idx + 2) % 3]];
                    T lu = Q.dist(U), lv = Q.dist(V);
                    if (lu > T(1e-14) && lv > T(1e-14)) {
                        T ca = std::acos(std::max(T(-1), std::min(T(1), ((U.x-Q.x)*(V.x-Q.x)+(U.y-Q.y)*(V.y-Q.y))/(lu*lv)))) * T(180) / T(3.14159265358979323846);
                        if (ca > bestCornerA) { bestCornerA = ca; forceIdx = idx; }
                    }
                }
            }
            if (forceIdx >= 0) {
                int oppEdge = (forceIdx + 1) % 3;
                for (int pos = 0; pos < 3; ++pos) {
                    if (order[pos] == oppEdge && pos > 0) {
                        std::swap(order[0], order[pos]);
                        break;
                    }
                }
                minAllowedChildAng = T(1.5);
            }
            if (trySplitTriangleEdges(mesh, cand.ti, order, degree, minAllowedChildAng)) {
                didSplit = true;
                break;
            }
        }
        if (!didSplit) break;
    }
}


template <class T>
inline bool evalMeshQuadUnion(
    const TriMeshData<T>& mesh,
    const int qv[4], const int qe[4],
    T minPairAng, T maxPairAng,
    T& effMinOut, T& scoreOut,
    bool allowCurvedInteriorChord = false)
{
    gsPoint2D<T> qc[4] = {mesh.pts[qv[0]], mesh.pts[qv[1]], mesh.pts[qv[2]], mesh.pts[qv[3]]};
    const bool strictConv = isStrictlyConvexQuad(qc);
    if (!strictConv) {
        if (!allowCurvedInteriorChord) return false;
        // For curved-chord pairing, allow mild non-convexity ONLY if signed area is positive
        // and at least one edge is interior so a curved Bezier chord can restore strict convexity.
        T area2 = cross2(qc[0], qc[1], qc[2]) + cross2(qc[0], qc[2], qc[3]);
        if (area2 <= T(1e-7)) return false;
        bool hasInt = false;
        for (int k = 0; k < 4; ++k) if (!mesh.edges[qe[k]].isBoundary) hasInt = true;
        if (!hasInt) return false;
    }

    for (int k = 0; k < 4; ++k) {
        const auto& er = mesh.edges[qe[k]];
        if (!er.isBoundary || !er.seg.bsp) continue;
        const gsPoint2D<T>& p1 = qc[(k + 1) % 4];
        const gsPoint2D<T>& p2 = qc[(k + 2) % 4];
        const gsPoint2D<T>& p3 = qc[(k + 3) % 4];
        const gsPoint2D<T>& p0 = qc[k];
        bool e1Int = !mesh.edges[qe[(k + 1) % 4]].isBoundary;
        bool e2Int = !mesh.edges[qe[(k + 2) % 4]].isBoundary;
        bool e3Int = !mesh.edges[qe[(k + 3) % 4]].isBoundary;
        for (T u : {T(0.25), T(0.5), T(0.75)}) {
            gsPoint2D<T> sPt = evalCurvePoint(*er.seg.bsp, u);
            if ((!allowCurvedInteriorChord || !e1Int) && cross2(p1, p2, sPt) <= T(1e-8)) return false;
            if ((!allowCurvedInteriorChord || !e2Int) && cross2(p2, p3, sPt) <= T(1e-8)) return false;
            if ((!allowCurvedInteriorChord || !e3Int) && cross2(p3, p0, sPt) <= T(1e-8)) return false;
        }
    }

    auto edgeDir = [&](int ei, int fromV, int toV, bool atEnd) -> gsPoint2D<T> {
        const auto& er = mesh.edges[ei];
        if (er.seg.bsp && (er.isBoundary || allowCurvedInteriorChord)) {
            bool fwd = (er.a == fromV);
            T u = (atEnd == fwd) ? T(1) : T(0);
            gsPoint2D<T> dv = evalCurveDeriv(*er.seg.bsp, u);
            T sgn = ((fwd && !atEnd) || (!fwd && atEnd)) ? T(1) : T(-1);
            T nm = std::hypot(dv.x, dv.y);
            if (nm > T(1e-14)) return gsPoint2D<T>(sgn * dv.x / nm, sgn * dv.y / nm);
        }
        T chSign = atEnd ? T(-1) : T(1);
        T dx = chSign * (mesh.pts[toV].x - mesh.pts[fromV].x);
        T dy = chSign * (mesh.pts[toV].y - mesh.pts[fromV].y);
        T nm = std::hypot(dx, dy);
        return (nm > T(1e-14)) ? gsPoint2D<T>(dx / nm, dy / nm) : gsPoint2D<T>(1, 0);
    };

    const T radToDeg = T(180) / T(3.14159265358979323846);
    T minAng = T(360), maxAng = T(0), maxDev90 = T(0);
    T lenMin = T(1e30), lenMax = T(0);
    for (int k = 0; k < 4; ++k) {
        int kn = (k + 1) % 4, kp = (k + 3) % 4;
        gsPoint2D<T> tOut  = edgeDir(qe[k],  qv[k],  qv[kn], false);
        gsPoint2D<T> tBack = edgeDir(qe[kp], qv[kp], qv[k],  true);
        T cr = tOut.x * tBack.y - tOut.y * tBack.x;
        T dt = tOut.x * tBack.x + tOut.y * tBack.y;
        T ang = std::atan2(cr, dt) * radToDeg;
        if (ang <= T(0)) ang += T(360);
        bool bothBnd = mesh.edges[qe[k]].isBoundary && mesh.edges[qe[kp]].isBoundary;
        T locMin = (allowCurvedInteriorChord && !bothBnd) ? std::max(T(3), minPairAng - T(6)) : minPairAng;
        T locMax = (allowCurvedInteriorChord && !bothBnd) ? std::min(T(178), maxPairAng + T(14)) : maxPairAng;
        if (ang < locMin || ang > locMax) return false;
        minAng = std::min(minAng, ang);
        maxAng = std::max(maxAng, ang);
        maxDev90 = std::max(maxDev90, std::abs(ang - T(90)));
        T eLen = mesh.pts[qv[k]].dist(mesh.pts[qv[kn]]);
        lenMin = std::min(lenMin, eLen);
        lenMax = std::max(lenMax, eLen);
    }

    effMinOut = std::min(minAng, T(180) - maxAng);
    T arRatio = (lenMax > T(1e-14)) ? (lenMin / lenMax) : T(1);
    scoreOut = effMinOut + T(2.25) * arRatio - T(0.15) * maxDev90 - (strictConv ? T(0) : T(12));
    return true;
}

// Co-optimize interior diagonals where a convex flip unlocks a substantially better adjacent quad union.
template <class T>
inline void flipForBetterQuadUnions(TriMeshData<T>& mesh, T minPairAng, T maxPairAng)
{
    const int nT = static_cast<int>(mesh.tris.size());
    const int nE = static_cast<int>(mesh.edges.size());
    if (nT < 2 || nE == 0) return;

    for (int pass = 0; pass < 2; ++pass) {
        struct EdgeAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
        std::vector<EdgeAdj> eAdj(nE);
        std::vector<int> vertDeg(mesh.pts.size(), 0);
        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            for (int k = 0; k < 3; ++k) {
                ++vertDeg[tr.v[k]];
                int eid = tr.e[k];
                if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
                else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
            }
        }

        auto edgePairEff = [&](int eIdx, int tA, int kA, int tB, int kB) -> T {
            const auto& A = mesh.tris[tA];
            const auto& B = mesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) return T(0);
            int qv[4] = {c, a, d, b};
            int qe[4] = {A.e[(kA + 2) % 3], B.e[(kB + 1) % 3], B.e[(kB + 2) % 3], A.e[(kA + 1) % 3]};
            T eff = 0, sc = 0;
            return evalMeshQuadUnion(mesh, qv, qe, minPairAng, maxPairAng, eff, sc) ? eff : T(0);
        };

        std::vector<T> triBestEff(nT, T(0));
        for (int eIdx = 0; eIdx < nE; ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = eAdj[eIdx].t0, tB = eAdj[eIdx].t1;
            if (tA < 0 || tB < 0) continue;
            T eff = edgePairEff(eIdx, tA, eAdj[eIdx].k0, tB, eAdj[eIdx].k1);
            triBestEff[tA] = std::max(triBestEff[tA], eff);
            triBestEff[tB] = std::max(triBestEff[tB], eff);
        }

        bool anyFlip = false;
        for (int eIdx = 0; eIdx < nE; ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
            int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
            if (tA < 0 || tB < 0) continue;
            T curBest = std::max(triBestEff[tA], triBestEff[tB]);
            if (curBest >= T(32)) continue;

            const auto A = mesh.tris[tA], B = mesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) continue;
            if (vertDeg[a] <= minRequiredBoundaryTriValence(mesh, a) ||
                vertDeg[b] <= minRequiredBoundaryTriValence(mesh, b))
                continue;

            gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
            if (!isStrictlyConvexQuad(qa)) continue;
            if (!isSegmentInsideDomain(mesh.pts[c], mesh.pts[d], mesh.loopPolys)) continue;

            T mn0, mx0, mn1, mx1;
            triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[a], mesh.pts[d], mn0, mx0);
            triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[d], mesh.pts[b], mn1, mx1);
            if (std::min(mn0, mn1) < T(10) || std::max(mx0, mx1) > T(158)) continue;

            int e_bc = A.e[(kA + 1) % 3], e_ca = A.e[(kA + 2) % 3];
            int e_ad = B.e[(kB + 1) % 3], e_db = B.e[(kB + 2) % 3];

            bool bndSafe = true;
            for (int chk : {e_ca, e_ad, e_db, e_bc}) {
                const auto& er = mesh.edges[chk];
                if (!er.isBoundary || !er.seg.bsp) continue;
                T reqSign = (chk == e_ca || chk == e_ad) ? T(-1) : T(1);
                for (T u : {T(0.25), T(0.5), T(0.75)}) {
                    if (reqSign * cross2(mesh.pts[c], mesh.pts[d], evalCurvePoint(*er.seg.bsp, u)) <= T(1e-8)) {
                        bndSafe = false; break;
                    }
                }
                if (!bndSafe) break;
            }
            if (!bndSafe) continue;

            // Temporarily apply flip and inspect union options around the 5 cavity edges
            auto oldEdge = mesh.edges[eIdx];
            auto oldTA = mesh.tris[tA], oldTB = mesh.tris[tB];
            mesh.edges[eIdx] = {c, d, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1};
            mesh.tris[tA] = {{c, a, d}, {e_ca, e_ad, eIdx}};
            mesh.tris[tB] = {{c, d, b}, {eIdx, e_db, e_bc}};

            T newBest = T(0);
            auto testEdgeSlot = [&](int eid, int triIdx, int slotIdx) {
                if (mesh.edges[eid].isBoundary) return;
                int ot = (eAdj[eid].t0 == triIdx) ? eAdj[eid].t1 : eAdj[eid].t0;
                int ok = (eAdj[eid].t0 == triIdx) ? eAdj[eid].k1 : eAdj[eid].k0;
                if (ot < 0) {
                    if (eid == eIdx) newBest = std::max(newBest, edgePairEff(eIdx, tA, 2, tB, 0));
                    return;
                }
                newBest = std::max(newBest, edgePairEff(eid, triIdx, slotIdx, ot, ok));
            };
            testEdgeSlot(eIdx, tA, 2);
            testEdgeSlot(e_ca, tA, 0);
            testEdgeSlot(e_ad, tA, 1);
            testEdgeSlot(e_db, tB, 1);
            testEdgeSlot(e_bc, tB, 2);

            bool accept = (curBest <= T(1e-6) && newBest >= T(16)) || (newBest > curBest + T(6));
            if (accept) {
                eAdj[eIdx] = {tA, 2, tB, 0};
                if (eAdj[e_ca].t0 == tA) eAdj[e_ca].k0 = 0; else eAdj[e_ca].k1 = 0;
                if (eAdj[e_ad].t0 == tB) { eAdj[e_ad].t0 = tA; eAdj[e_ad].k0 = 1; } else { eAdj[e_ad].t1 = tA; eAdj[e_ad].k1 = 1; }
                if (eAdj[e_db].t0 == tB) eAdj[e_db].k0 = 1; else eAdj[e_db].k1 = 1;
                if (eAdj[e_bc].t0 == tA) { eAdj[e_bc].t0 = tB; eAdj[e_bc].k0 = 2; } else { eAdj[e_bc].t1 = tB; eAdj[e_bc].k1 = 2; }
                --vertDeg[a]; --vertDeg[b]; ++vertDeg[c]; ++vertDeg[d];
                triBestEff[tA] = std::max(triBestEff[tA], newBest);
                triBestEff[tB] = std::max(triBestEff[tB], newBest);
                anyFlip = true;
            } else {
                mesh.edges[eIdx] = oldEdge;
                mesh.tris[tA] = oldTA;
                mesh.tris[tB] = oldTB;
            }
        }
        if (!anyFlip) break;
    }
}

// Construct a clamped cubic Bezier interior CurveSeg elevated to `degree`.
template <class T>
inline CurveSeg<T> makeCubicBezierInteriorSeg(
    const gsPoint2D<T>& q0,
    const gsPoint2D<T>& p1,
    const gsPoint2D<T>& p2,
    const gsPoint2D<T>& q3,
    int degree)
{
    gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
    gsMatrix<T> coefs(4, 2);
    coefs << q0.x, q0.y, p1.x, p1.y, p2.x, p2.y, q3.x, q3.y;
    gsBSpline<T> bsp(kv, coefs);
    if (bsp.degree() < degree)
        bsp.degreeElevate(degree - bsp.degree());
    return CurveSeg<T>(bsp);
}

// Build a smooth cubic Hermite interior spoke O -> M departing O along unit direction `tO`
// and arriving at M along incoming normal `tM_inward` (pointing from M back toward interior),
// guaranteed to stay strictly inside `loopPolys`. Falls back to straight if clearance fails.
template <class T>
inline CurveSeg<T> makeSymmetricHermiteSpoke(
    const gsPoint2D<T>& O,
    const gsPoint2D<T>& M,
    const gsPoint2D<T>& tO,
    const gsPoint2D<T>& tM_inward,
    const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
    int degree)
{
    const T L = O.dist(M);
    if (L <= T(1e-10)) return makeStraightSeg(O, M, degree);
    T nO = std::hypot(tO.x, tO.y), nM = std::hypot(tM_inward.x, tM_inward.y);
    if (nO <= T(1e-12) || nM <= T(1e-12)) return makeStraightSeg(O, M, degree);
    gsPoint2D<T> dO(tO.x / nO, tO.y / nO);
    gsPoint2D<T> dM(tM_inward.x / nM, tM_inward.y / nM);
    // Avoid creating hairpin loops if tangent points strongly backward (> 115 deg from chord)
    gsPoint2D<T> ch((M.x - O.x) / L, (M.y - O.y) / L);
    if (dO.x * ch.x + dO.y * ch.y < T(-0.35)) dO = ch;
    if (dM.x * (-ch.x) + dM.y * (-ch.y) < T(-0.35)) dM = gsPoint2D<T>(-ch.x, -ch.y);

    for (T alpha : {T(0.32), T(0.22), T(0.13)}) {
        gsPoint2D<T> p1(O.x + alpha * L * dO.x, O.y + alpha * L * dO.y);
        gsPoint2D<T> p2(M.x + alpha * L * dM.x, M.y + alpha * L * dM.y);
        CurveSeg<T> cs = makeCubicBezierInteriorSeg(O, p1, p2, M, degree);
        bool ok = true;
        for (int s = 1; s < 10; ++s) {
            gsPoint2D<T> pt = evalCurvePoint(*cs.bsp, static_cast<T>(s) / T(10));
            if (!isPointStrictlyInsideDomain(pt, loopPolys, T(1e-5))) { ok = false; break; }
        }
        if (ok) return cs;
    }
    return makeStraightSeg(O, M, degree);
}

// Edmonds' general-graph Blossom augmenting-path search from an unmatched vertex `root`.
// Returns true and commits one augmenting path (subject to `canCommit` valence guard) if found.
template <class CanCommitFn>
inline bool augmentBlossomMatchingFrom(
    int root,
    int nV,
    const std::vector<std::vector<std::pair<int, int>>>& adj, // (neighborV, pairEdgeIdx)
    std::vector<int>& partner,
    std::vector<int>& matchedPair,
    CanCommitFn&& canCommit)
{
    if (partner[root] >= 0) return false;

    std::vector<int> p(nV, -1);
    std::vector<int> pEdge(nV, -1);
    std::vector<int> base(nV);
    for (int i = 0; i < nV; ++i) base[i] = i;
    std::vector<bool> used(nV, false);
    std::vector<bool> blossom(nV, false);

    auto lca = [&](int a, int b) -> int {
        std::vector<bool> vis(nV, false);
        for (;;) {
            a = base[a];
            vis[a] = true;
            if (partner[a] < 0) break;
            a = p[partner[a]];
        }
        for (;;) {
            b = base[b];
            if (vis[b]) return b;
            b = p[partner[b]];
        }
    };

    auto markPath = [&](int v, int b, int children, int edgeToChild) {
        while (base[v] != b) {
            blossom[base[v]] = blossom[base[partner[v]]] = true;
            p[v] = children;
            pEdge[v] = edgeToChild;
            children = partner[v];
            edgeToChild = matchedPair[v];
            v = p[partner[v]];
        }
    };

    std::queue<int> q;
    used[root] = true;
    q.push(root);
    int foundEnd = -1;

    while (!q.empty() && foundEnd < 0) {
        int v = q.front();
        q.pop();
        for (const auto& nb : adj[v]) {
            int to = nb.first;
            int pi = nb.second;
            if (base[v] == base[to] || partner[v] == to) continue;
            if (to == root || (partner[to] >= 0 && p[partner[to]] >= 0)) {
                int curbase = lca(v, to);
                std::fill(blossom.begin(), blossom.end(), false);
                markPath(v, curbase, to, pi);
                markPath(to, curbase, v, pi);
                for (int i = 0; i < nV; ++i) {
                    if (blossom[base[i]]) {
                        base[i] = curbase;
                        if (!used[i]) {
                            used[i] = true;
                            q.push(i);
                        }
                    }
                }
            } else if (p[to] < 0) {
                p[to] = v;
                pEdge[to] = pi;
                if (partner[to] < 0) {
                    foundEnd = to;
                    break;
                }
                int m = partner[to];
                used[m] = true;
                q.push(m);
            }
        }
    }

    if (foundEnd < 0) return false;

    // Reconstruct added/removed pair edges along the augmenting path
    std::vector<int> addPairs, remPairs;
    std::vector<std::tuple<int, int, int>> flips; // (u, v, addPi)
    int cur = foundEnd;
    int guard = 0;
    while (cur >= 0 && guard++ <= nV) {
        int pv = p[cur];
        if (pv < 0) return false;
        int addPi = pEdge[cur];
        int ppv = partner[pv];
        if (ppv >= 0) remPairs.push_back(matchedPair[pv]);
        addPairs.push_back(addPi);
        flips.emplace_back(pv, cur, addPi);
        cur = ppv;
    }
    if (guard > nV || !canCommit(addPairs, remPairs))
        return false;

    for (const auto& f : flips) {
        int u = std::get<0>(f), v = std::get<1>(f), addPi = std::get<2>(f);
        partner[u] = v;
        partner[v] = u;
        matchedPair[u] = addPi;
        matchedPair[v] = addPi;
    }
    return true;
}


} // namespace internal_splitter
} // namespace gismo
