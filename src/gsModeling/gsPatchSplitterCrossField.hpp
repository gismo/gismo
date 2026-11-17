/** @file gsPatchSplitterCrossField.hpp
    @brief Constrained Multi-Hole Triangulation, 4-RoSy Harmonic Cross-Field Solver,
           Separatrix Tracing, Strip Collapse, and Conforming Arrangement Quadification.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

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
    };

    std::vector<gsPoint2D<T>> pts;
    std::vector<bool> isBoundaryPt;
    std::vector<bool> isSmoothBnd;
    std::vector<bool> isSharpCorner;
    std::vector<T> cornerAngleDeg;
    std::vector<gsPoint2D<T>> bndTanIn;
    std::vector<gsPoint2D<T>> bndTanOut;
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

template <class T>
inline bool buildConstrainedDomainTriangulation(
    const gsPlanarDomain<T>& dom,
    int degree,
    T targetH,
    bool refineEdges,
    int earMode,
    TriMeshData<T>& outMesh,
    T maxBendRad = T(0.44))
{
    if (dom.numLoops() == 0) return false;
    const int nLoops = dom.numLoops();

    outMesh = TriMeshData<T>();
    outMesh.loopPolys.resize(nLoops);
    for (int l = 0; l < nLoops; ++l)
        outMesh.loopPolys[l] = sampleCurveLoopPoly(dom.loop(l), 28);

    // Prepare boundary segments for each loop (ensuring >= 4 segments per loop)
    struct LoopSeg { CurveSeg<T> seg; int loopIdx; };
    std::vector<std::vector<CurveSeg<T>>> loopSegs(nLoops);

    for (int l = 0; l < nLoops; ++l) {
        const auto& lp = dom.loop(l);
        const int K = lp.size();
        if (K == 0) return false;

        std::vector<CurveSeg<T>> raw;
        if (K == 1) {
            auto sH = CurveSeg<T>(lp.curve(0), degree).splitArc(T(0.5));
            auto s1 = sH.first.splitArc(T(0.5));
            auto s2 = sH.second.splitArc(T(0.5));
            raw = {s1.first, s1.second, s2.first, s2.second};
        } else if (K == 2) {
            auto s0 = CurveSeg<T>(lp.curve(0), degree).splitArc(T(0.5));
            auto s1 = CurveSeg<T>(lp.curve(1), degree).splitArc(T(0.5));
            raw = {s0.first, s0.second, s1.first, s1.second};
        } else {
            for (int i = 0; i < K; ++i)
                raw.emplace_back(lp.curve(i), degree);
        }

        if (refineEdges && targetH > T(1e-8)) {
            T loopPerim = T(0);
            std::vector<T> segLens(raw.size(), T(0));
            for (size_t si = 0; si < raw.size(); ++si) {
                gsPoint2D<T> pm = evalCurvePoint(*raw[si].bsp, T(0.5));
                segLens[si] = raw[si].p0().dist(pm) + pm.dist(raw[si].p1());
                loopPerim += segLens[si];
            }
            T localH = targetH;
            if (l > 0 && loopPerim > T(1e-8))
                localH = std::min(targetH, loopPerim / T(8));

            for (size_t si = 0; si < raw.size(); ++si) {
                const auto& s = raw[si];
                T len = segLens[si];
                auto d0 = evalCurveDeriv(*s.bsp, T(0));
                auto d1 = evalCurveDeriv(*s.bsp, T(1));
                T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
                T bend = T(0);
                if (n0 > T(1e-14) && n1 > T(1e-14)) {
                    T dot = std::max(T(-1), std::min(T(1), (d0.x * d1.x + d0.y * d1.y) / (n0 * n1)));
                    bend = std::acos(dot);
                }
                int nLen  = static_cast<int>(std::ceil(len / localH));
                int nBend = static_cast<int>(std::ceil(bend / maxBendRad));
                int nSub = std::max(1, std::min(12, std::max(nLen, nBend)));
                if (nSub == 1) {
                    loopSegs[l].push_back(s);
                } else {
                    CurveSeg<T> rem = s;
                    for (int k = 0; k < nSub - 1; ++k) {
                        T frac = T(1) / static_cast<T>(nSub - k);
                        auto sp = rem.splitArc(frac);
                        loopSegs[l].push_back(sp.first);
                        rem = sp.second;
                    }
                    loopSegs[l].push_back(rem);
                }
            }
        } else {
            loopSegs[l] = std::move(raw);
        }

        // Rebuild loopPolys directly from loopSegs so every ring vertex coincides with loopPolys
        outMesh.loopPolys[l].clear();
        for (const auto& seg : loopSegs[l]) {
            for (int k = 0; k < 16; ++k) {
                T u = static_cast<T>(k) / T(16);
                outMesh.loopPolys[l].push_back(evalCurvePoint(*seg.bsp, u));
            }
        }
    }

    // Populate vertices and boundary edges for every loop
    struct RingNode { int ptIdx; CurveSeg<T> outSeg; int edgeId; };
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
    for (int l = 0; l < nLoops; ++l) {
        int M = static_cast<int>(loopSegs[l].size());
        int baseIdx = static_cast<int>(outMesh.pts.size());
        for (int i = 0; i < M; ++i) {
            outMesh.pts.push_back(loopSegs[l][i].p0());
            outMesh.isBoundaryPt.push_back(true);

            const auto& prevBsp = *loopSegs[l][(i + M - 1) % M].bsp;
            const auto& nextBsp = *loopSegs[l][i].bsp;
            auto din  = evalCurveDeriv(prevBsp, T(1));
            auto dout = evalCurveDeriv(nextBsp, T(0));
            T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
            gsPoint2D<T> tIn  = (n1 > T(1e-14)) ? gsPoint2D<T>(din.x / n1, din.y / n1) : gsPoint2D<T>(1, 0);
            gsPoint2D<T> tOut = (n2 > T(1e-14)) ? gsPoint2D<T>(dout.x / n2, dout.y / n2) : gsPoint2D<T>(1, 0);
            T cr = tIn.x * tOut.y - tIn.y * tOut.x;
            T dt = tIn.x * tOut.x + tIn.y * tOut.y;
            T turnRad = std::atan2(cr, dt);
            bool smooth = (std::abs(turnRad) <= T(0.10));
            // Interior domain angle at vertex i:
            // Outer loop is CCW -> interior angle = pi - turnRad; Hole loop is CW -> interior angle = pi + turnRad
            T intAngDeg = (l == 0) ? (T(180) - turnRad * radToDeg) : (T(180) + turnRad * radToDeg);
            outMesh.isSmoothBnd.push_back(smooth);
            outMesh.isSharpCorner.push_back(std::abs(turnRad) > T(0.45)); // > ~25 deg turn
            outMesh.cornerAngleDeg.push_back(intAngDeg);
            outMesh.bndTanIn.push_back(tIn);
            outMesh.bndTanOut.push_back(tOut);
        }
        for (int i = 0; i < M; ++i) {
            int u = baseIdx + i;
            int v = baseIdx + (i + 1) % M;
            int eid = addUndirEdge(u, v, loopSegs[l][i], true, l);
            rings[l].push_back({u, loopSegs[l][i], eid});
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

            // Score candidate bridges by distance and departure angle (avoid grazing < 14 deg slits)
            T bestCost = std::numeric_limits<T>::max();
            for (int pass = 0; pass < 2 && bestHPos < 0; ++pass) {
                for (int hp = 0; hp < hSize; ++hp) {
                    int hmIdx = hRing[hp].ptIdx;
                    const auto& M = outMesh.pts[hmIdx];
                    const auto& mIn  = outMesh.bndTanIn[hmIdx];
                    const auto& mOut = outMesh.bndTanOut[hmIdx];

                    for (int op = 0; op < static_cast<int>(activeRing.size()); ++op) {
                        int ovIdx = activeRing[op].ptIdx;
                        const auto& V = outMesh.pts[ovIdx];
                        if (pass == 0 && V.x < M.x - T(1e-6)) continue;

                        T d = std::hypot(V.x - M.x, V.y - M.y);
                        if (d <= T(1e-12) || d >= bestCost) continue;

                        T ux = (V.x - M.x) / d, uy = (V.y - M.y) / d;
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
                        if (pass == 0 && minSin < T(0.24)) continue;

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

            if (bestHPos < 0 || bestOuterPos < 0) return false;

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

            newRing.push_back({vIdx, segVM, bridgeEid});
            for (int k = 0; k < hSize; ++k)
                newRing.push_back(hRing[(bestHPos + k) % hSize]);
            newRing.push_back({mIdx, segMV, bridgeEid});

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
                    const auto& bsp = *activeRing[k].outSeg.bsp;
                    gsMatrix<T> ev, par(1, 3);
                    par(0, 0) = T(0.25); par(0, 1) = T(0.5); par(0, 2) = T(0.75);
                    bsp.eval_into(par, ev);
                    for (int s = 0; s < 3; ++s) {
                        if (pointInTri(A, B, C, gsPoint2D<T>(ev(0, s), ev(1, s)))) {
                            clear = false;
                            break;
                        }
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
            if (subDivs >= 24) return false;
            int splitK = -1;
            T maxSag = -T(1);
            int longestK = -1;
            T maxLen = T(0);
            for (int k = 0; k < m; ++k) {
                int u = activeRing[k].ptIdx;
                int v = activeRing[(k + 1) % m].ptIdx;
                int eOld = activeRing[k].edgeId;
                if (!outMesh.edges[eOld].isBoundary) continue;
                auto sp = activeRing[k].outSeg.splitArc(T(0.5));
                gsPoint2D<T> mid = sp.first.p1();
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
            outMesh.edges[eOld] = {u, newIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, loopIdxOld};
            edgeMap.erase({std::min(u, v), std::max(u, v)});
            edgeMap[{std::min(u, newIdx), std::max(u, newIdx)}] = eOld;
            int eSecond = addUndirEdge(newIdx, v, sp.second, true, loopIdxOld);

            activeRing[splitK].outSeg = sp.first;
            activeRing[splitK].edgeId = eOld;
            activeRing.insert(activeRing.begin() + splitK + 1, {newIdx, sp.second, eSecond});
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
        activeRing[iPrev].outSeg = (er2.a == ia && er2.b == ic) ? er2.seg : er2.seg.reversed();
        activeRing[iPrev].edgeId = e2;
        activeRing.erase(activeRing.begin() + bestIdx);
    }

    if (activeRing.size() != 3) return false;
    int ia = activeRing[0].ptIdx, ib = activeRing[1].ptIdx, ic = activeRing[2].ptIdx;
    if (cross2(outMesh.pts[ia], outMesh.pts[ib], outMesh.pts[ic]) <= T(1e-8)) return false;
    outMesh.tris.push_back({{ia, ib, ic}, {activeRing[0].edgeId, activeRing[1].edgeId, activeRing[2].edgeId}});
    return true;
}

template <class T>
inline void applyDomainLawsonFlips(TriMeshData<T>& mesh)
{
    const int nT = static_cast<int>(mesh.tris.size());
    if (nT < 2) return;
    for (int pass = 0; pass < nT * 2; ++pass) {
        bool flipped = false;
        for (int eIdx = 0; eIdx < static_cast<int>(mesh.edges.size()); ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = -1, kA = -1, tB = -1, kB = -1;
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].e[k] == eIdx) {
                        if (tA < 0) { tA = ti; kA = k; }
                        else        { tB = ti; kB = k; }
                    }
                }
            }
            if (tA < 0 || tB < 0) continue;
            const auto A = mesh.tris[tA];
            const auto B = mesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
            if (!isStrictlyConvexQuad(qa)) continue;
            if (!isSegmentInsideDomain(mesh.pts[c], mesh.pts[d], mesh.loopPolys)) continue;

            auto minTriSin = [&](int i0, int i1, int i2) -> T {
                const auto& P0 = mesh.pts[i0];
                const auto& P1 = mesh.pts[i1];
                const auto& P2 = mesh.pts[i2];
                T cr2 = cross2(P0, P1, P2);
                if (cr2 <= T(1e-8)) return -T(1);
                T l2 = (P0.x - P1.x) * (P0.x - P1.x) + (P0.y - P1.y) * (P0.y - P1.y)
                     + (P1.x - P2.x) * (P1.x - P2.x) + (P1.y - P2.y) * (P1.y - P2.y)
                     + (P2.x - P0.x) * (P2.x - P0.x) + (P2.y - P0.y) * (P2.y - P0.y);
                return cr2 / (l2 + T(1e-14));
            };
            if (std::min(minTriSin(c, a, d), minTriSin(c, d, b)) <= std::min(minTriSin(a, b, c), minTriSin(b, a, d)) + T(1e-4))
                continue;

            int degA = 0, degB = 0;
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].v[k] == a) ++degA;
                    if (mesh.tris[ti].v[k] == b) ++degB;
                }
            }
            if ((mesh.isSmoothBnd[a] && degA <= 2) || (mesh.isSmoothBnd[b] && degB <= 2)) continue;

            int e_bc = A.e[(kA + 1) % 3];
            int e_ca = A.e[(kA + 2) % 3];
            int e_ad = -1, e_db = -1;
            for (int k = 0; k < 3; ++k) {
                int u = B.v[k], v = B.v[(k + 1) % 3];
                if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
            }
            if (e_ad < 0 || e_db < 0) continue;

            mesh.edges[eIdx] = {c, d, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1};
            mesh.tris[tA] = {{c, a, d}, {e_ca, e_ad, eIdx}};
            mesh.tris[tB] = {{c, d, b}, {eIdx, e_db, e_bc}};
            flipped = true;
            break;
        }
        if (!flipped) break;
    }
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
            gsPoint2D<T> old = mesh.pts[v];
            mesh.pts[v] = avg;
            bool ok = isPointStrictlyInsideDomain(avg, mesh.loopPolys, T(1e-6));
            for (int ti : adjTris[v]) {
                if (!ok) break;
                const auto& t = mesh.tris[ti];
                if (cross2(mesh.pts[t.v[0]], mesh.pts[t.v[1]], mesh.pts[t.v[2]]) <= T(1e-8))
                    ok = false;
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
            if (!isPointStrictlyInsideDomain(c, mesh.loopPolys, T(1e-6))) continue;

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
inline void refineNeedleCircumcenters(TriMeshData<T>& mesh, T targetH, int maxRounds = 2, T minTriAngThresh = T(16), bool allowCentroidFallback = false)
{
    const T radToDeg = T(180) / T(3.14159265358979323846);
    for (int round = 0; round < maxRounds; ++round) {
        int nOrig = static_cast<int>(mesh.tris.size());
        bool added = false;
        for (int ti = 0; ti < nOrig; ++ti) {
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
            if (lMax <= T(2.3) * lMin || lMax <= T(0.35) * targetH) continue;

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

            gsPoint2D<T> cc((P0.x + P1.x + P2.x) / T(3),
                            (P0.y + P1.y + P2.y) / T(3));
            T ax = P1.x - P0.x, ay = P1.y - P0.y;
            T bx = P2.x - P0.x, by = P2.y - P0.y;
            T D = T(2) * (ax * by - ay * bx);
            bool ccValid = false;
            if (std::abs(D) > T(1e-14)) {
                T aSq = ax * ax + ay * ay, bSq = bx * bx + by * by;
                gsPoint2D<T> cand(P0.x + (by * aSq - ay * bSq) / D,
                                  P0.y + (ax * bSq - bx * aSq) / D);
                T h0 = cross2(P0, P1, cand) / (l0 + T(1e-14));
                T h1 = cross2(P1, P2, cand) / (l1 + T(1e-14));
                T h2 = cross2(P2, P0, cand) / (l2 + T(1e-14));
                if (std::min(h0, std::min(h1, h2)) >= T(0.25) * lMin) {
                    cc = cand;
                    ccValid = true;
                }
            }
            if (!ccValid && !allowCentroidFallback) continue;
            if (!isPointStrictlyInsideDomain(cc, mesh.loopPolys, T(1e-6))) continue;

            int vc = static_cast<int>(mesh.pts.size());
            mesh.pts.push_back(cc);
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
    smoothInteriorVertices(mesh, 4);
}

template <class T>
inline void bisectObtuseTriangles(TriMeshData<T>& mesh, int degree, T maxQuadAngleDeg = T(146), int maxSplits = 28, T minQuadAngleDeg = T(14.5))
{
    auto triQuadMaxAngle = [&](const typename TriMeshData<T>::Tri& t) -> T {
        const auto& A = mesh.pts[t.v[0]];
        const auto& B = mesh.pts[t.v[1]];
        const auto& C = mesh.pts[t.v[2]];
        auto midOf = [&](int eIdx) -> gsPoint2D<T> {
            const auto& er = mesh.edges[eIdx];
            if (er.seg.bsp) return er.seg.splitArc(T(0.5)).first.p1();
            return gsPoint2D<T>(T(0.5) * (mesh.pts[er.a].x + mesh.pts[er.b].x),
                                T(0.5) * (mesh.pts[er.a].y + mesh.pts[er.b].y));
        };
        gsPoint2D<T> m01 = midOf(t.e[0]);
        gsPoint2D<T> m12 = midOf(t.e[1]);
        gsPoint2D<T> m20 = midOf(t.e[2]);
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
        return std::max(quadScore(q0, isCornerSharp(0)),
                        std::max(quadScore(q1, isCornerSharp(1)), quadScore(q2, isCornerSharp(2))));
    };

    for (int iter = 0; iter < maxSplits; ++iter) {
        if (mesh.tris.size() > 500) break;
        int badTi = -1;
        T worstAng = maxQuadAngleDeg;
        for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
            T a = triQuadMaxAngle(mesh.tris[ti]);
            if (a > worstAng) { worstAng = a; badTi = ti; }
        }
        if (badTi < 0) break;

        const auto bt = mesh.tris[badTi];
        int order[3] = {0, 1, 2};
        T lens[3];
        for (int k = 0; k < 3; ++k)
            lens[k] = mesh.pts[bt.v[k]].dist(mesh.pts[bt.v[(k + 1) % 3]]);
        std::sort(order, order + 3, [&](int a, int b) { return lens[a] > lens[b]; });

        bool didSplit = false;
        for (int oi = 0; oi < 3 && !didSplit; ++oi) {
            int kMax = order[oi];
            int u = bt.v[kMax], v = bt.v[(kMax + 1) % 3], opp = bt.v[(kMax + 2) % 3];
            int eIdx = bt.e[kMax];
            int e_vo = bt.e[(kMax + 1) % 3];
            int e_ou = bt.e[(kMax + 2) % 3];

            auto oldEr = mesh.edges[eIdx];
            if (!oldEr.seg.bsp)
                oldEr.seg = makeStraightSeg(mesh.pts[oldEr.a], mesh.pts[oldEr.b], degree);
            auto sp = oldEr.seg.splitArc(T(0.5));
            gsPoint2D<T> midPt = sp.first.p1();

            if (oldEr.isBoundary && !isSegmentInsideDomain(mesh.pts[opp], midPt, mesh.loopPolys))
                continue;

            int mateTi = -1, mateK = -1;
            for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
                if (ti == badTi) continue;
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].e[k] == eIdx) { mateTi = ti; mateK = k; break; }
                }
                if (mateTi >= 0) break;
            }
            if (mateTi >= 0) {
                int mOpp = mesh.tris[mateTi].v[(mateK + 2) % 3];
                if (oldEr.isBoundary && !isSegmentInsideDomain(mesh.pts[mOpp], midPt, mesh.loopPolys))
                    continue;
            }

            int mIdx = static_cast<int>(mesh.pts.size());
            mesh.pts.push_back(midPt);
            mesh.isBoundaryPt.push_back(oldEr.isBoundary);
            mesh.isSmoothBnd.push_back(oldEr.isBoundary);
            mesh.isSharpCorner.push_back(false);
            mesh.cornerAngleDeg.push_back(oldEr.isBoundary ? T(180) : T(360));
            auto din = evalCurveDeriv(*sp.first.bsp, T(1));
            T dn = std::hypot(din.x, din.y);
            gsPoint2D<T> tMid = (dn > T(1e-14)) ? gsPoint2D<T>(din.x / dn, din.y / dn) : gsPoint2D<T>(1, 0);
            mesh.bndTanIn.push_back(tMid);
            mesh.bndTanOut.push_back(tMid);

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
                int mu = mt.v[mateK], mv = mt.v[(mateK + 1) % 3], mOpp = mt.v[(mateK + 2) % 3];
                int e_mvo = mt.e[(mateK + 1) % 3];
                int e_mou = mt.e[(mateK + 2) % 3];
                int e_mom = static_cast<int>(mesh.edges.size());
                mesh.edges.push_back({mOpp, mIdx, makeStraightSeg(mesh.pts[mOpp], midPt, degree),
                                      CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
                mesh.tris[mateTi] = {{mu, mIdx, mOpp}, {halfFor(mu), e_mom, e_mou}};
                mesh.tris.push_back({{mIdx, mv, mOpp}, {halfFor(mv), e_mvo, e_mom}});
            }

            applyDomainLawsonFlips(mesh);
            didSplit = true;
        }
        if (!didSplit) break;
    }
}

// ============================================================================
// 4-RoSy Harmonic Cross-Field Solver (L z = 0) & Singularity Extraction
// ============================================================================

template <class T>
struct CrossFieldSolution {
    std::vector<std::complex<T>> z;
    struct Singularity {
        gsPoint2D<T> pos;
        int index = 0; // +1 (val 3), -1 (val 5), +4 (O-grid disk)
        std::vector<T> portAngles;
    };
    std::vector<Singularity> singularities;
};

template <class T>
inline CrossFieldSolution<T> solveHarmonic4RoSy(const TriMeshData<T>& mesh)
{
    CrossFieldSolution<T> sol;
    const int nV = static_cast<int>(mesh.pts.size());
    sol.z.assign(nV, std::complex<T>(T(1), T(0)));
    if (nV < 3 || mesh.tris.empty()) return sol;

    std::vector<bool> isDirichlet(nV, false);
    std::vector<std::complex<T>> dirVal(nV, std::complex<T>(T(1), T(0)));

    for (int i = 0; i < nV; ++i) {
        if (!mesh.isBoundaryPt[i]) continue;
        std::complex<T> t1(mesh.bndTanIn[i].x, mesh.bndTanIn[i].y);
        std::complex<T> t2(mesh.bndTanOut[i].x, mesh.bndTanOut[i].y);
        std::complex<T> z1 = std::pow(t1, 4);
        std::complex<T> z2 = std::pow(t2, 4);

        if (!mesh.isSharpCorner[i]) {
            isDirichlet[i] = true;
            std::complex<T> avg = z1 + z2;
            T nrm = std::abs(avg);
            dirVal[i] = (nrm > T(1e-12)) ? (avg / nrm) : z2;
        } else if (std::abs(z1 - z2) < T(0.35)) {
            // 90 deg or 270 deg corner: incident tangents agree in 4-RoSy space
            isDirichlet[i] = true;
            std::complex<T> avg = z1 + z2;
            T nrm = std::abs(avg);
            dirVal[i] = (nrm > T(1e-12)) ? (avg / nrm) : z2;
        }
    }

    std::vector<int> freeIdx(nV, -1);
    int nFree = 0;
    for (int i = 0; i < nV; ++i)
        if (!isDirichlet[i]) freeIdx[i] = nFree++;

    if (nFree == 0) {
        sol.z = dirVal;
    } else {
        gsSparseEntries<T> entries;
        entries.reserve(mesh.tris.size() * 12);
        gsMatrix<T> rhs = gsMatrix<T>::Zero(nFree, 2);

        auto addWeight = [&](int u, int v, T w) {
            if (isDirichlet[u]) return;
            int iu = freeIdx[u];
            entries.add(iu, iu, w);
            if (isDirichlet[v]) {
                rhs(iu, 0) += w * dirVal[v].real();
                rhs(iu, 1) += w * dirVal[v].imag();
            } else {
                entries.add(iu, freeIdx[v], -w);
            }
        };

        for (const auto& t : mesh.tris) {
            const auto& p0 = mesh.pts[t.v[0]];
            const auto& p1 = mesh.pts[t.v[1]];
            const auto& p2 = mesh.pts[t.v[2]];
            T area2 = std::max(T(1e-10), std::abs(cross2(p0, p1, p2)));

            auto cotAngle = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C) -> T {
                T dot = (B.x - A.x) * (C.x - A.x) + (B.y - A.y) * (C.y - A.y);
                return std::max(T(1e-4), dot / area2);
            };
            T w12 = T(0.5) * cotAngle(p0, p1, p2);
            T w20 = T(0.5) * cotAngle(p1, p2, p0);
            T w01 = T(0.5) * cotAngle(p2, p0, p1);

            addWeight(t.v[1], t.v[2], w12);
            addWeight(t.v[2], t.v[1], w12);
            addWeight(t.v[2], t.v[0], w20);
            addWeight(t.v[0], t.v[2], w20);
            addWeight(t.v[0], t.v[1], w01);
            addWeight(t.v[1], t.v[0], w01);
        }

        gsSparseMatrix<T> L(nFree, nFree);
        L.setFrom(entries);
        typename gsSparseSolver<T>::SimplicialLDLT solver;
        solver.compute(L);
        gsMatrix<T> solMat = solver.solve(rhs);

        for (int i = 0; i < nV; ++i) {
            if (isDirichlet[i]) {
                sol.z[i] = dirVal[i];
            } else {
                int ii = freeIdx[i];
                sol.z[i] = std::complex<T>(solMat(ii, 0), solMat(ii, 1));
            }
        }
    }

    // Extract triangle singularities via angle winding
    auto wrapPi = [](T a) -> T {
        const T pi = T(3.14159265358979323846);
        while (a > pi) a -= T(2) * pi;
        while (a <= -pi) a += T(2) * pi;
        return a;
    };

    struct RawSing { gsPoint2D<T> pos; int idx; };
    std::vector<RawSing> rawSings;
    const T twoPi = T(2) * T(3.14159265358979323846);

    for (const auto& t : mesh.tris) {
        std::complex<T> z0 = sol.z[t.v[0]], z1 = sol.z[t.v[1]], z2 = sol.z[t.v[2]];
        if (std::abs(z0) < T(1e-12) || std::abs(z1) < T(1e-12) || std::abs(z2) < T(1e-12))
            continue;
        T a0 = std::arg(z0), a1 = std::arg(z1), a2 = std::arg(z2);
        T sum = wrapPi(a1 - a0) + wrapPi(a2 - a1) + wrapPi(a0 - a2);
        int ind = static_cast<int>(std::llround(sum / twoPi));
        if (ind == 0) continue;

        // Solve linear zero z(u,v) = 0 in triangle barycentrics
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
        gsPoint2D<T> pos(l0 * P0.x + l1 * P1.x + l2 * P2.x,
                         l0 * P0.y + l1 * P1.y + l2 * P2.y);
        if (isPointStrictlyInsideDomain(pos, mesh.loopPolys, T(1e-5)))
            rawSings.push_back({pos, ind});
    }

    // Cluster nearby singularities
    T avgEdge = T(0);
    for (const auto& e : mesh.edges)
        avgEdge += mesh.pts[e.a].dist(mesh.pts[e.b]);
    avgEdge /= std::max<T>(T(1), static_cast<T>(mesh.edges.size()));
    const T clusterTol = T(1.8) * avgEdge;

    std::vector<bool> used(rawSings.size(), false);
    for (size_t i = 0; i < rawSings.size(); ++i) {
        if (used[i]) continue;
        used[i] = true;
        gsPoint2D<T> sumP = rawSings[i].pos;
        int sumInd = rawSings[i].idx;
        int cnt = 1;
        for (size_t j = i + 1; j < rawSings.size(); ++j) {
            if (!used[j] && rawSings[i].pos.dist(rawSings[j].pos) <= clusterTol) {
                used[j] = true;
                sumP.x += rawSings[j].pos.x;
                sumP.y += rawSings[j].pos.y;
                sumInd += rawSings[j].idx;
                ++cnt;
            }
        }
        if (sumInd != 0) {
            typename CrossFieldSolution<T>::Singularity sg;
            sg.pos = gsPoint2D<T>(sumP.x / cnt, sumP.y / cnt);
            sg.index = sumInd;
            sol.singularities.push_back(sg);
        }
    }

    return sol;
}

// Evaluate cross-field angle theta in [0, pi/2) at arbitrary point inside mesh
template <class T>
inline bool evalCrossFieldAngleAt(const TriMeshData<T>& mesh,
                                  const CrossFieldSolution<T>& cf,
                                  const gsPoint2D<T>& p,
                                  T& thetaOut,
                                  int& triOut)
{
    for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
        const auto& t = mesh.tris[ti];
        const auto& A = mesh.pts[t.v[0]];
        const auto& B = mesh.pts[t.v[1]];
        const auto& C = mesh.pts[t.v[2]];
        T area2 = cross2(A, B, C);
        if (std::abs(area2) < T(1e-14)) continue;
        T l0 = cross2(B, C, p) / area2;
        T l1 = cross2(C, A, p) / area2;
        T l2 = T(1) - l0 - l1;
        if (l0 >= -T(1e-4) && l1 >= -T(1e-4) && l2 >= -T(1e-4)) {
            std::complex<T> z = l0 * cf.z[t.v[0]] + l1 * cf.z[t.v[1]] + l2 * cf.z[t.v[2]];
            if (std::abs(z) < T(1e-12)) return false;
            thetaOut = T(0.25) * std::arg(z);
            triOut = ti;
            return true;
        }
    }
    return false;
}

// ============================================================================
// Planar Arrangement Construction & Bulletproof Conforming Quadification
// ============================================================================

template <class T>
inline bool buildCrossFieldArrangementCandidate(
    const gsPlanarDomain<T>& dom,
    int degree,
    T stripCollapseRatio,
    T angleThresholdDeg,
    gsMultiPatch<T>& outMp,
    std::string& methodLabel)
{
    const gsMatrix<T> bbox = dom.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    if (diag <= T(1e-12)) return false;
    const T targetH = diag / T(22);

    TriMeshData<T> fineMesh;
    if (!buildConstrainedDomainTriangulation(dom, degree, targetH, true, 0, fineMesh))
        return false;
    applyDomainLawsonFlips(fineMesh);
    refineInteriorTriangles(fineMesh, targetH, 2);

    CrossFieldSolution<T> cf = solveHarmonic4RoSy(fineMesh);

    // Check if an O-grid +4 singularity (or 4x +1 / 2x +2 valence-3 cluster on a smooth disk/oval) occurred
    bool hasDiskCluster = false;
    gsPoint2D<T> diskCenter(0, 0);
    int totalSingIdx = 0;
    gsPoint2D<T> posSum(0, 0);
    int posCnt = 0;
    for (const auto& sg : cf.singularities) {
        totalSingIdx += sg.index;
        if (sg.index > 0) {
            posSum.x += sg.pos.x * sg.index;
            posSum.y += sg.pos.y * sg.index;
            posCnt += sg.index;
        }
        if (sg.index >= 2)
            hasDiskCluster = true;
    }
    int nSharpCorners = 0;
    for (size_t i = 0; i < fineMesh.pts.size(); ++i)
        if (fineMesh.isBoundaryPt[i] && fineMesh.isSharpCorner[i]) ++nSharpCorners;

    gsPoint2D<T> polyCtr = polygonCentroid(fineMesh.loopPolys[0]);
    if (posCnt > 0)
        diskCenter = gsPoint2D<T>(posSum.x / posCnt, posSum.y / posCnt);
    else
        diskCenter = polyCtr;

    if (!hasDiskCluster && dom.numHoles() == 0 && nSharpCorners == 0)
        hasDiskCluster = true;

    // For smooth 0-hole domains where harmonic cross-field detects a +4 center zero (or valence-3 cluster),
    // construct the canonical 5-patch butterfly O-grid directly aligned with principal cross directions.
    if (dom.numHoles() == 0 && dom.outer().size() <= 4 && nSharpCorners == 0 && hasDiskCluster) {
        std::vector<CurveSeg<T>> bEdges;
        const auto& outer = dom.outer();
        if (outer.size() == 1) {
            auto sH = CurveSeg<T>(outer.curve(0), degree).splitArc(T(0.5));
            auto s1 = sH.first.splitArc(T(0.5));
            auto s2 = sH.second.splitArc(T(0.5));
            bEdges = {s1.first, s1.second, s2.first, s2.second};
        } else if (outer.size() == 2) {
            auto s0 = CurveSeg<T>(outer.curve(0), degree).splitArc(T(0.5));
            auto s1 = CurveSeg<T>(outer.curve(1), degree).splitArc(T(0.5));
            bEdges = {s0.first, s0.second, s1.first, s1.second};
        } else if (outer.size() == 4) {
            for (int i = 0; i < 4; ++i) bEdges.emplace_back(outer.curve(i), degree);
        }

        if (bEdges.size() == 4) {
            gsPoint2D<T> P[4];
            for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
            gsPoint2D<T> pAvg((P[0].x + P[1].x + P[2].x + P[3].x) * T(0.25),
                              (P[0].y + P[1].y + P[2].y + P[3].y) * T(0.25));
            std::vector<gsPoint2D<T>> cands = {diskCenter, polyCtr, pAvg};
            const T alphas[3] = {T(0.45), T(0.35), T(0.55)};

            T bestMaxAng = T(360);
            gsMultiPatch<T> bestOGrid;
            for (const auto& ctr : cands) {
                if (!isPointStrictlyInsideDomain(ctr, fineMesh.loopPolys, T(1e-5))) continue;
                for (T alpha : alphas) {
                    gsPoint2D<T> Q[4];
                    bool ok = true;
                    for (int k = 0; k < 4; ++k) {
                        Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                            ctr.y + alpha * (P[k].y - ctr.y));
                        if (!isPointStrictlyInsideDomain(Q[k], fineMesh.loopPolys, T(1e-5))) { ok = false; break; }
                    }
                    if (!ok || !isStrictlyConvexQuad(Q)) continue;

                    CurveSeg<T> qEdge[4], spoke[4];
                    for (int k = 0; k < 4; ++k) {
                        qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], degree);
                        spoke[k] = makeStraightSeg(Q[k], P[k], degree);
                    }
                    gsMultiPatch<T> candMp;
                    candMp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], degree)));
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        candMp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(qEdge[k], spoke[kn], bEdges[k].reversed(), spoke[k].reversed(), degree)));
                    }

                    bool allOk = true;
                    T curMaxA = T(0);
                    for (size_t p = 0; p < candMp.nPatches() && allOk; ++p) {
                        if (!isPatchDetSingleSign(candMp.patch(p), 17, T(1e-7))) { allOk = false; break; }
                        T angs[4], mnA = 0, mxA = 180;
                        bool degen = false;
                        evalPatchCornerAnglesDeg(candMp.patch(p), angs, mnA, mxA, degen);
                        if (degen || mxA > T(152)) allOk = false;
                        curMaxA = std::max(curMaxA, mxA);
                    }
                    if (allOk && curMaxA < bestMaxAng) {
                        bestMaxAng = curMaxA;
                        bestOGrid = give(candMp);
                    }
                }
            }

            if (bestOGrid.nPatches() == 5) {
                outMp = give(bestOGrid);
                methodLabel = "HarmonicCrossField";
                return true;
            }
        }
    }

    // General Cross-Field Separatrix Tracing from singularities and concave corners
    struct RaySeed {
        gsPoint2D<T> start;
        gsPoint2D<T> dir;
        int sourceId;
    };
    std::vector<RaySeed> seeds;

    // 1. Seed from concave boundary corners (> angleThresholdDeg)
    for (int i = 0; i < static_cast<int>(fineMesh.pts.size()); ++i) {
        if (!fineMesh.isBoundaryPt[i] || !fineMesh.isSharpCorner[i]) continue;
        if (fineMesh.cornerAngleDeg[i] <= angleThresholdDeg) continue;

        gsPoint2D<T> tBack(-fineMesh.bndTanIn[i].x, -fineMesh.bndTanIn[i].y);
        gsPoint2D<T> tFwd(fineMesh.bndTanOut[i].x, fineMesh.bndTanOut[i].y);
        // On CCW outer loop, interior sweeps CCW from tFwd to tBack
        T angStart = std::atan2(tFwd.y, tFwd.x);
        T sweepRad = fineMesh.cornerAngleDeg[i] * T(3.14159265358979323846) / T(180);
        int nPorts = std::max(1, static_cast<int>(std::floor(fineMesh.cornerAngleDeg[i] / T(90) + T(0.5))) - 1);
        for (int p = 1; p <= nPorts; ++p) {
            T ang = angStart + sweepRad * (static_cast<T>(p) / static_cast<T>(nPorts + 1));
            seeds.push_back({fineMesh.pts[i], gsPoint2D<T>(std::cos(ang), std::sin(ang)), -(i + 1)});
        }
    }

    // 2. Seed from interior singularities (valence-3 -> 3 ports, valence-5 -> 5 ports)
    const T pi = T(3.14159265358979323846);
    for (size_t s = 0; s < cf.singularities.size(); ++s) {
        const auto& sg = cf.singularities[s];
        int val = (sg.index > 0) ? 3 : 5;
        T baseAng = T(0);
        T thetaLoc = T(0);
        int triDummy = -1;
        gsPoint2D<T> probe(sg.pos.x + T(1.5) * targetH, sg.pos.y);
        if (evalCrossFieldAngleAt(fineMesh, cf, probe, thetaLoc, triDummy))
            baseAng = thetaLoc;
        for (int p = 0; p < val; ++p) {
            T ang = baseAng + (T(2) * pi * p) / static_cast<T>(val);
            seeds.push_back({sg.pos, gsPoint2D<T>(std::cos(ang), std::sin(ang)), static_cast<int>(s + 1)});
        }
    }

    // Trace each separatrix through the cross field until hitting boundary or snapping to singularity/corner
    struct TracedPath {
        std::vector<gsPoint2D<T>> pts;
        int hitBndLoop = -1;
        int hitBndEdge = -1;
        gsPoint2D<T> hitPt;
        int sourceId = 0;
    };
    std::vector<TracedPath> paths;
    const T stepLen = T(0.35) * targetH;
    const int maxSteps = 200;

    for (const auto& sd : seeds) {
        TracedPath tp;
        tp.sourceId = sd.sourceId;
        tp.pts.push_back(sd.start);
        gsPoint2D<T> cur = sd.start;
        gsPoint2D<T> dir = sd.dir;
        bool terminated = false;

        for (int st = 0; st < maxSteps && !terminated; ++st) {
            // Update direction from cross-field (branch closest to current dir)
            T th = T(0);
            int triDummy = -1;
            if (st > 2 && evalCrossFieldAngleAt(fineMesh, cf, cur, th, triDummy)) {
                T bestDot = -T(2);
                gsPoint2D<T> bestDir = dir;
                for (int k = 0; k < 4; ++k) {
                    T a = th + k * (pi * T(0.5));
                    gsPoint2D<T> cand(std::cos(a), std::sin(a));
                    T d = cand.x * dir.x + cand.y * dir.y;
                    if (d > bestDot) { bestDot = d; bestDir = cand; }
                }
                if (bestDot > T(0.5)) dir = bestDir;
            }

            gsPoint2D<T> nxt(cur.x + stepLen * dir.x, cur.y + stepLen * dir.y);

            // Check intersection with boundary loop polygons
            T bestT = T(2);
            gsPoint2D<T> bestHit;
            int bestLoop = -1;
            for (int l = 0; l < static_cast<int>(fineMesh.loopPolys.size()); ++l) {
                const auto& poly = fineMesh.loopPolys[l];
                int np = static_cast<int>(poly.size());
                for (int i = 0; i < np; ++i) {
                    if (st == 0 && sd.sourceId < 0 && cur.dist(poly[i]) < T(1e-6)) continue;
                    const auto& A = poly[i];
                    const auto& B = poly[(i + 1) % np];
                    T rx = nxt.x - cur.x, ry = nxt.y - cur.y;
                    T sx = B.x - A.x, sy = B.y - A.y;
                    T denom = rx * sy - ry * sx;
                    if (std::abs(denom) < T(1e-14)) continue;
                    T t = ((A.x - cur.x) * sy - (A.y - cur.y) * sx) / denom;
                    T u = ((A.x - cur.x) * ry - (A.y - cur.y) * rx) / denom;
                    if (t > T(1e-4) && t <= T(1) + T(1e-4) && u >= -T(1e-4) && u <= T(1) + T(1e-4)) {
                        if (t < bestT) {
                            bestT = t;
                            bestHit = gsPoint2D<T>(cur.x + t * rx, cur.y + t * ry);
                            bestLoop = l;
                        }
                    }
                }
            }

            if (bestLoop >= 0) {
                tp.pts.push_back(bestHit);
                tp.hitBndLoop = bestLoop;
                tp.hitPt = bestHit;
                terminated = true;
                break;
            }

            // Snap to existing singularity if close after leaving source
            if (st > 3) {
                for (size_t s = 0; s < cf.singularities.size(); ++s) {
                    if (sd.sourceId == static_cast<int>(s + 1)) continue;
                    if (nxt.dist(cf.singularities[s].pos) < T(1.1) * targetH) {
                        tp.pts.push_back(cf.singularities[s].pos);
                        tp.hitPt = cf.singularities[s].pos;
                        terminated = true;
                        break;
                    }
                }
            }
            if (terminated) break;

            tp.pts.push_back(nxt);
            cur = nxt;
        }
        if (terminated && tp.pts.size() >= 2)
            paths.push_back(std::move(tp));
    }

    // Topological Strip Collapse: merge boundary hits that lie within stripCollapseRatio * diag / 4
    const T snapDist = stripCollapseRatio * targetH * T(3.0);
    for (size_t i = 0; i < paths.size(); ++i) {
        if (paths[i].hitBndLoop < 0) continue;
        for (size_t j = i + 1; j < paths.size(); ++j) {
            if (paths[j].hitBndLoop == paths[i].hitBndLoop &&
                paths[i].hitPt.dist(paths[j].hitPt) < snapDist) {
                gsPoint2D<T> mid(T(0.5) * (paths[i].hitPt.x + paths[j].hitPt.x),
                                 T(0.5) * (paths[i].hitPt.y + paths[j].hitPt.y));
                paths[i].hitPt = mid;
                paths[j].hitPt = mid;
                paths[i].pts.back() = mid;
                paths[j].pts.back() = mid;
            }
        }
    }

    // If no singularities or separatrices were extracted (e.g., convex quad or clean polygon),
    // return false so structural/triangulation candidate evaluation takes over seamlessly.
    if (paths.empty()) return false;

    // ------------------------------------------------------------------------
    // Build Planar Arrangement from Exact Boundary Loops + Traced Separatrices
    // ------------------------------------------------------------------------
    struct ArrVertex {
        gsPoint2D<T> pt;
        bool isBoundary = false;
    };
    struct ArrEdge {
        int u = -1, v = -1;
        CurveSeg<T> seg; // oriented u -> v
        bool isBoundary = false;
        int loopIdx = -1;
    };

    std::vector<ArrVertex> arrVerts;
    std::vector<ArrEdge> arrEdges;

    auto addArrVertex = [&](const gsPoint2D<T>& p, bool isBnd, T mergeTol) -> int {
        for (size_t i = 0; i < arrVerts.size(); ++i) {
            if (arrVerts[i].pt.dist(p) <= mergeTol) {
                if (isBnd) arrVerts[i].isBoundary = true;
                return static_cast<int>(i);
            }
        }
        int id = static_cast<int>(arrVerts.size());
        arrVerts.push_back({p, isBnd});
        return id;
    };

    const T vMergeTol = T(0.15) * targetH;

    // Collect boundary loops as CurveSegs and track split parameters from path hits
    const int nLoops = dom.numLoops();
    std::vector<std::vector<CurveSeg<T>>> loopBaseSegs(nLoops);
    std::vector<std::vector<std::vector<T>>> loopSplitUs(nLoops);
    for (int l = 0; l < nLoops; ++l) {
        const auto& lp = dom.loop(l);
        int K = lp.size();
        if (K == 1) {
            auto sH = CurveSeg<T>(lp.curve(0), degree).splitArc(T(0.5));
            auto s1 = sH.first.splitArc(T(0.5));
            auto s2 = sH.second.splitArc(T(0.5));
            loopBaseSegs[l] = {s1.first, s1.second, s2.first, s2.second};
        } else if (K == 2) {
            auto s0 = CurveSeg<T>(lp.curve(0), degree).splitArc(T(0.5));
            auto s1 = CurveSeg<T>(lp.curve(1), degree).splitArc(T(0.5));
            loopBaseSegs[l] = {s0.first, s0.second, s1.first, s1.second};
        } else {
            for (int c = 0; c < K; ++c)
                loopBaseSegs[l].emplace_back(lp.curve(c), degree);
        }
        loopSplitUs[l].resize(loopBaseSegs[l].size());
    }

    // Snap boundary hits onto exact boundary CurveSegs
    struct ChordCut { int vStart = -1, vEnd = -1; gsPoint2D<T> pStart, pEnd; };
    std::vector<ChordCut> chordCuts;

    for (auto& path : paths) {
        if (path.hitBndLoop >= 0) {
            int l = path.hitBndLoop;
            int bestC = -1;
            T bestU = T(0.5);
            T bestD = std::numeric_limits<T>::max();
            for (size_t c = 0; c < loopBaseSegs[l].size(); ++c) {
                T u = loopBaseSegs[l][c].closestParam(path.hitPt);
                gsPoint2D<T> q = evalCurvePoint(*loopBaseSegs[l][c].bsp, u);
                T d = q.dist(path.hitPt);
                if (d < bestD) { bestD = d; bestC = static_cast<int>(c); bestU = u; }
            }
            if (bestC >= 0) {
                if (bestU < T(0.08)) {
                    path.hitPt = loopBaseSegs[l][bestC].p0();
                } else if (bestU > T(0.92)) {
                    path.hitPt = loopBaseSegs[l][bestC].p1();
                } else {
                    path.hitPt = evalCurvePoint(*loopBaseSegs[l][bestC].bsp, bestU);
                    loopSplitUs[l][bestC].push_back(bestU);
                }
            }
        }
        if (path.pts.front().dist(path.hitPt) > vMergeTol) {
            if (isSegmentInsideDomain(path.pts.front(), path.hitPt, fineMesh.loopPolys)) {
                chordCuts.push_back({-1, -1, path.pts.front(), path.hitPt});
            }
        }
    }

    if (chordCuts.empty()) return false;

    // Split boundary segments at hit parameters and add boundary arrangement edges
    for (int l = 0; l < nLoops; ++l) {
        for (size_t c = 0; c < loopBaseSegs[l].size(); ++c) {
            auto& us = loopSplitUs[l][c];
            std::sort(us.begin(), us.end());
            std::vector<T> cleanUs;
            for (T u : us) {
                if (u > T(0.05) && u < T(0.95) &&
                    (cleanUs.empty() || u - cleanUs.back() > T(0.05)))
                    cleanUs.push_back(u);
            }
            std::vector<CurveSeg<T>> pieces;
            CurveSeg<T> rem = loopBaseSegs[l][c];
            T prevU = T(0);
            for (T u : cleanUs) {
                T locU = (u - prevU) / std::max(T(1e-9), T(1) - prevU);
                auto sp = rem.splitAtParam(locU);
                pieces.push_back(sp.first);
                rem = sp.second;
                prevU = u;
            }
            pieces.push_back(rem);
            for (const auto& seg : pieces) {
                int va = addArrVertex(seg.p0(), true, vMergeTol);
                int vb = addArrVertex(seg.p1(), true, vMergeTol);
                if (va != vb)
                    arrEdges.push_back({va, vb, seg, true, l});
            }
        }
    }

    // Intersect interior chord cuts cleanly so crossings become degree-4 vertices (0 T-junctions)
    std::vector<std::vector<T>> cutSplits(chordCuts.size());
    for (size_t i = 0; i < chordCuts.size(); ++i) {
        for (size_t j = i + 1; j < chordCuts.size(); ++j) {
            const auto& A = chordCuts[i].pStart;
            const auto& B = chordCuts[i].pEnd;
            const auto& C = chordCuts[j].pStart;
            const auto& D = chordCuts[j].pEnd;
            T rx = B.x - A.x, ry = B.y - A.y;
            T sx = D.x - C.x, sy = D.y - C.y;
            T denom = rx * sy - ry * sx;
            if (std::abs(denom) < T(1e-12)) continue;
            T t = ((C.x - A.x) * sy - (C.y - A.y) * sx) / denom;
            T u = ((C.x - A.x) * ry - (C.y - A.y) * rx) / denom;
            if (t > T(0.02) && t < T(0.98) && u > T(0.02) && u < T(0.98)) {
                cutSplits[i].push_back(t);
                cutSplits[j].push_back(u);
            }
        }
    }

    std::set<std::pair<int, int>> existingPairs;
    for (const auto& ae : arrEdges)
        existingPairs.insert({std::min(ae.u, ae.v), std::max(ae.u, ae.v)});

    for (size_t i = 0; i < chordCuts.size(); ++i) {
        auto& ts = cutSplits[i];
        ts.push_back(T(0));
        ts.push_back(T(1));
        std::sort(ts.begin(), ts.end());
        std::vector<T> uniq;
        for (T t : ts) {
            if (uniq.empty() || t - uniq.back() > T(0.03))
                uniq.push_back(t);
        }
        if (uniq.size() < 2) uniq = {T(0), T(1)};
        uniq.front() = T(0);
        uniq.back() = T(1);

        const auto& A = chordCuts[i].pStart;
        const auto& B = chordCuts[i].pEnd;
        for (size_t k = 0; k + 1 < uniq.size(); ++k) {
            gsPoint2D<T> pa((T(1) - uniq[k]) * A.x + uniq[k] * B.x,
                            (T(1) - uniq[k]) * A.y + uniq[k] * B.y);
            gsPoint2D<T> pb((T(1) - uniq[k + 1]) * A.x + uniq[k + 1] * B.x,
                            (T(1) - uniq[k + 1]) * A.y + uniq[k + 1] * B.y);
            int va = addArrVertex(pa, false, vMergeTol);
            int vb = addArrVertex(pb, false, vMergeTol);
            if (va == vb) continue;
            auto key = std::make_pair(std::min(va, vb), std::max(va, vb));
            if (existingPairs.count(key)) continue;
            existingPairs.insert(key);
            arrEdges.push_back({va, vb, makeStraightSeg(arrVerts[va].pt, arrVerts[vb].pt, degree), false, -1});
        }
    }

    // Prune degree-1 interior vertices (dangling rays that did not terminate)
    bool pruned = true;
    while (pruned) {
        pruned = false;
        std::vector<int> deg(arrVerts.size(), 0);
        for (const auto& e : arrEdges) { ++deg[e.u]; ++deg[e.v]; }
        for (size_t i = 0; i < arrEdges.size(); ++i) {
            if (arrEdges[i].isBoundary) continue;
            if (deg[arrEdges[i].u] <= 1 || deg[arrEdges[i].v] <= 1) {
                arrEdges.erase(arrEdges.begin() + static_cast<long>(i));
                pruned = true;
                break;
            }
        }
    }

    if (arrEdges.size() < 4) return false;

    // Build half-edge structure to extract minimal bounded faces
    struct HalfEdge {
        int origin = -1, target = -1;
        int arrEdgeIdx = -1;
        bool forward = true;
        int twin = -1;
        int next = -1;
        T angle = T(0);
    };

    const int nE = static_cast<int>(arrEdges.size());
    std::vector<HalfEdge> halfEdges(2 * nE);
    std::vector<std::vector<int>> outHalf(arrVerts.size());

    for (int i = 0; i < nE; ++i) {
        const auto& e = arrEdges[i];
        int h0 = 2 * i, h1 = 2 * i + 1;
        const auto& pu = arrVerts[e.u].pt;
        const auto& pv = arrVerts[e.v].pt;
        halfEdges[h0] = {e.u, e.v, i, true,  h1, -1, std::atan2(pv.y - pu.y, pv.x - pu.x)};
        halfEdges[h1] = {e.v, e.u, i, false, h0, -1, std::atan2(pu.y - pv.y, pu.x - pv.x)};
        outHalf[e.u].push_back(h0);
        outHalf[e.v].push_back(h1);
    }

    for (size_t v = 0; v < arrVerts.size(); ++v) {
        auto& oh = outHalf[v];
        std::sort(oh.begin(), oh.end(), [&](int a, int b) {
            return halfEdges[a].angle < halfEdges[b].angle;
        });
    }

    for (int h = 0; h < 2 * nE; ++h) {
        int v = halfEdges[h].target;
        int tw = halfEdges[h].twin;
        const auto& oh = outHalf[v];
        auto it = std::find(oh.begin(), oh.end(), tw);
        if (it == oh.end() || oh.empty()) return false;
        size_t pos = static_cast<size_t>(it - oh.begin());
        size_t prevPos = (pos + oh.size() - 1) % oh.size();
        halfEdges[h].next = oh[prevPos];
    }

    std::vector<bool> heVisited(2 * nE, false);
    std::vector<std::vector<int>> faces;
    for (int h = 0; h < 2 * nE; ++h) {
        if (heVisited[h]) continue;
        std::vector<int> cycle;
        int cur = h;
        bool okCycle = true;
        for (int guard = 0; guard < 2 * nE + 2; ++guard) {
            if (heVisited[cur]) { okCycle = (cur == h); break; }
            heVisited[cur] = true;
            cycle.push_back(cur);
            cur = halfEdges[cur].next;
            if (cur < 0) { okCycle = false; break; }
        }
        if (!okCycle || cycle.size() < 3) continue;

        // Compute signed area & centroid to keep only bounded interior faces inside domain
        std::vector<gsPoint2D<T>> poly;
        poly.reserve(cycle.size());
        for (int ch : cycle) poly.push_back(arrVerts[halfEdges[ch].origin].pt);
        T area2 = T(0);
        for (size_t k = 0; k < poly.size(); ++k) {
            const auto& p = poly[k];
            const auto& q = poly[(k + 1) % poly.size()];
            area2 += p.x * q.y - q.x * p.y;
        }
        if (area2 <= T(1e-7)) continue; // Skip exterior unbounded face & hole faces

        gsPoint2D<T> ctr = polygonCentroid(poly);
        if (!isPointStrictlyInsideDomain(ctr, fineMesh.loopPolys, T(1e-6))) continue;
        faces.push_back(std::move(cycle));
    }

    if (faces.empty()) return false;

    auto getArrSeg = [&](int heIdx) -> CurveSeg<T> {
        const auto& he = halfEdges[heIdx];
        const auto& ae = arrEdges[he.arrEdgeIdx];
        return he.forward ? ae.seg : ae.seg.reversed();
    };

    bool allQuads = true;
    for (const auto& f : faces)
        if (f.size() != 4) { allQuads = false; break; }

    outMp.clear();
    if (allQuads) {
        for (const auto& f : faces) {
            outMp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(getArrSeg(f[0]), getArrSeg(f[1]),
                                      getArrSeg(f[2]), getArrSeg(f[3]), degree)));
        }
        methodLabel = "HarmonicCrossField";
        return true;
    }

    // Bulletproof Catmull-Clark midpoint + centroid quadification of every arrangement N-gon face
    struct SplitSegs { CurveSeg<T> h0, h1; gsPoint2D<T> mid; };
    std::vector<SplitSegs> edgeHalves(nE);
    for (int i = 0; i < nE; ++i) {
        auto sp = arrEdges[i].seg.splitArc(T(0.5));
        edgeHalves[i].h0 = sp.first;
        edgeHalves[i].h1 = sp.second;
        edgeHalves[i].mid = sp.first.p1();
    }

    auto getArrHalves = [&](int heIdx, CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
        const auto& he = halfEdges[heIdx];
        const auto& eh = edgeHalves[he.arrEdgeIdx];
        outMid = eh.mid;
        if (he.forward) {
            out0 = eh.h0;
            out1 = eh.h1;
        } else {
            out0 = eh.h1.reversed();
            out1 = eh.h0.reversed();
        }
    };

    for (const auto& f : faces) {
        int nf = static_cast<int>(f.size());
        std::vector<CurveSeg<T>> h0(nf), h1(nf);
        std::vector<gsPoint2D<T>> mids(nf), vertsF(nf);
        for (int k = 0; k < nf; ++k) {
            vertsF[k] = arrVerts[halfEdges[f[k]].origin].pt;
            getArrHalves(f[k], h0[k], h1[k], mids[k]);
        }
        gsPoint2D<T> O = polygonCentroid(mids);
        if (!isPointStrictlyInsideDomain(O, fineMesh.loopPolys, T(1e-6)))
            O = polygonCentroid(vertsF);

        std::vector<CurveSeg<T>> spokes(nf);
        for (int k = 0; k < nf; ++k)
            spokes[k] = makeStraightSeg(O, mids[k], degree);

        for (int k = 0; k < nf; ++k) {
            int kp = (k + nf - 1) % nf;
            outMp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spokes[kp], h1[kp], h0[k], spokes[k].reversed(), degree)));
        }
    }

    methodLabel = "HarmonicCrossField";
    return outMp.nPatches() > 0;
}

} // namespace internal_splitter
} // namespace gismo
