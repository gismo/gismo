/** @file gsPatchSplitterMotorcycle.hpp
    @brief Minimal-Patch Cross-Field Layout Engine:
           - Stage 0: Virtual G1-Consolidated CAD Boundary Loops
           - Stage 1: Curvature-Regularized Robin 4-RoSy + Singularity Dipole Cancellation
           - Stage 2A: Harmonic-Gradient Curved 5-Patch O-Grid (smooth non-star H=0 domains)
           - Stage 2B: Simultaneous Motorcycle T-Mesh Partitioning + Hole-Bridge Guarantee
           - Stage 3: Native Constrained T-Mesh Integer Quantization + Macro-Block Emission

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitterCrossField.hpp>
#include <gsModeling/gsSpringPatch.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace gismo {
namespace internal_splitter {

// ============================================================================
// Stage 0: Virtual G1-Consolidated CAD Boundary Loops
// ============================================================================

template <class T>
struct VirtualCADLoop {
    int loopIdx = -1;
    bool isHole = false;
    T totalArcLen = T(0);

    static constexpr int kArcSamples = 32;

    struct Piece {
        CurveSeg<T> seg;
        T s0 = T(0);
        T s1 = T(0);
        std::array<T, kArcSamples + 1> sTable{}; // cumulative local arc length on [0, 1]
    };
    std::vector<Piece> pieces;

    struct PinnedCorner {
        T arcS = T(0);
        int pieceIdx = -1;
        gsPoint2D<T> pt;
        T interiorAngleDeg = T(180);
        int targetPatchValence = 1;
    };
    std::vector<PinnedCorner> corners;
};

template <class T>
inline T virtualLoopWrapS(const VirtualCADLoop<T>& vl, T s)
{
    if (vl.totalArcLen <= T(1e-14)) return T(0);
    T r = std::fmod(s, vl.totalArcLen);
    if (r < T(0)) r += vl.totalArcLen;
    if (r >= vl.totalArcLen) r -= vl.totalArcLen;
    return r;
}

template <class T>
inline T virtualPieceArcToParam(const typename VirtualCADLoop<T>::Piece& pc, T localS)
{
    const int N = VirtualCADLoop<T>::kArcSamples;
    const T len = pc.sTable[N];
    if (len <= T(1e-14)) return T(0);
    if (localS <= T(0)) return T(0);
    if (localS >= len) return T(1);
    auto it = std::upper_bound(pc.sTable.begin(), pc.sTable.end(), localS);
    int idx = static_cast<int>(it - pc.sTable.begin()) - 1;
    idx = std::max(0, std::min(N - 1, idx));
    T sA = pc.sTable[idx], sB = pc.sTable[idx + 1];
    T frac = (sB > sA + T(1e-14)) ? (localS - sA) / (sB - sA) : T(0);
    return (static_cast<T>(idx) + frac) / static_cast<T>(N);
}

template <class T>
inline T virtualPieceParamToArc(const typename VirtualCADLoop<T>::Piece& pc, T u)
{
    const int N = VirtualCADLoop<T>::kArcSamples;
    u = std::max(T(0), std::min(T(1), u));
    T fu = u * static_cast<T>(N);
    int idx = std::min(N - 1, std::max(0, static_cast<int>(std::floor(fu))));
    T frac = fu - static_cast<T>(idx);
    return pc.sTable[idx] + frac * (pc.sTable[idx + 1] - pc.sTable[idx]);
}

template <class T>
inline int findVirtualPieceAtS(const VirtualCADLoop<T>& vl, T sWrapped, T& outLocalU)
{
    outLocalU = T(0);
    if (vl.pieces.empty()) return -1;
    for (size_t i = 0; i < vl.pieces.size(); ++i) {
        if (sWrapped < vl.pieces[i].s1 - T(1e-11) || i + 1 == vl.pieces.size()) {
            outLocalU = virtualPieceArcToParam<T>(vl.pieces[i], sWrapped - vl.pieces[i].s0);
            return static_cast<int>(i);
        }
    }
    outLocalU = T(1);
    return static_cast<int>(vl.pieces.size()) - 1;
}

template <class T>
inline gsPoint2D<T> evalVirtualLoopPoint(const VirtualCADLoop<T>& vl, T s)
{
    T u = T(0);
    int pi = findVirtualPieceAtS(vl, virtualLoopWrapS(vl, s), u);
    if (pi < 0 || !vl.pieces[pi].seg.bsp) return gsPoint2D<T>(0, 0);
    return evalCurvePoint(*vl.pieces[pi].seg.bsp, u);
}

template <class T>
inline gsPoint2D<T> evalVirtualLoopUnitTangent(const VirtualCADLoop<T>& vl, T s)
{
    T u = T(0);
    int pi = findVirtualPieceAtS(vl, virtualLoopWrapS(vl, s), u);
    if (pi < 0 || !vl.pieces[pi].seg.bsp) return gsPoint2D<T>(1, 0);
    gsPoint2D<T> d = evalCurveDeriv(*vl.pieces[pi].seg.bsp, u);
    T n = std::hypot(d.x, d.y);
    return (n > T(1e-14)) ? gsPoint2D<T>(d.x / n, d.y / n) : gsPoint2D<T>(1, 0);
}

template <class T>
inline T projectPointToVirtualLoop(const VirtualCADLoop<T>& vl,
                                   const gsPoint2D<T>& q,
                                   T& outArcS,
                                   gsPoint2D<T>& outPt)
{
    outArcS = T(0);
    outPt = q;
    T bestD = std::numeric_limits<T>::max();
    for (size_t i = 0; i < vl.pieces.size(); ++i) {
        const auto& pc = vl.pieces[i];
        if (!pc.seg.bsp) continue;
        T u = pc.seg.closestParam(q);
        gsPoint2D<T> p = evalCurvePoint(*pc.seg.bsp, u);
        T d = p.dist(q);
        if (d < bestD) {
            bestD = d;
            outPt = p;
            outArcS = virtualLoopWrapS(vl, pc.s0 + virtualPieceParamToArc<T>(pc, u));
        }
    }
    return bestD;
}

template <class T>
inline std::vector<VirtualCADLoop<T>> buildVirtualCADLoops(
    const gsPlanarDomain<T>& dom,
    int degree,
    T sharpThreshDeg = T(22.0))
{
    const int nLoops = dom.numLoops();
    std::vector<VirtualCADLoop<T>> vloops(nLoops);
    const T radToDeg = T(180) / T(3.14159265358979323846);

    for (int l = 0; l < nLoops; ++l) {
        auto& vl = vloops[l];
        vl.loopIdx = l;
        vl.isHole = (l > 0);
        const auto& lp = dom.loop(l);
        const int nC = lp.size();
        if (nC == 0) continue;

        std::vector<CurveSeg<T>> rawSegs;
        rawSegs.reserve(nC);
        for (int c = 0; c < nC; ++c)
            rawSegs.emplace_back(lp.curve(c), degree);

        // Enforce consistent loop orientation: outer (l==0) CCW (area > 0), holes (l>0) CW (area < 0)
        // so the domain interior is ALWAYS strictly to the LEFT of every loop!
        std::vector<gsPoint2D<T>> poly;
        poly.reserve(nC * 16);
        for (int c = 0; c < nC; ++c)
            for (int k = 0; k < 16; ++k)
                poly.push_back(evalCurvePoint(*rawSegs[c].bsp, static_cast<T>(k) / T(16)));
        T area2 = T(0);
        for (size_t i = 0; i < poly.size(); ++i) {
            const auto& a = poly[i];
            const auto& b = poly[(i + 1) % poly.size()];
            area2 += a.x * b.y - b.x * a.y;
        }
        bool needFlip = (l == 0) ? (area2 < T(0)) : (area2 > T(0));
        if (needFlip) {
            std::reverse(rawSegs.begin(), rawSegs.end());
            for (auto& s : rawSegs) s = s.reversed();
        }

        // Build arc-length tables
        vl.pieces.resize(nC);
        T cumS = T(0);
        const int N = VirtualCADLoop<T>::kArcSamples;
        for (int c = 0; c < nC; ++c) {
            auto& pc = vl.pieces[c];
            pc.seg = rawSegs[c];
            pc.s0 = cumS;
            pc.sTable[0] = T(0);
            gsPoint2D<T> prev = evalCurvePoint(*pc.seg.bsp, T(0));
            for (int k = 1; k <= N; ++k) {
                gsPoint2D<T> cur = evalCurvePoint(*pc.seg.bsp, static_cast<T>(k) / static_cast<T>(N));
                pc.sTable[k] = pc.sTable[k - 1] + prev.dist(cur);
                prev = cur;
            }
            cumS += pc.sTable[N];
            pc.s1 = cumS;
        }
        vl.totalArcLen = cumS;

        // Classify true geometric corners (dissolving smooth G1 joints where |ang - 180| < sharpThreshDeg)
        for (int c = 0; c < nC; ++c) {
            int cPrev = (c + nC - 1) % nC;
            gsPoint2D<T> dIn = evalCurveDeriv(*vl.pieces[cPrev].seg.bsp, T(1));
            gsPoint2D<T> dOut = evalCurveDeriv(*vl.pieces[c].seg.bsp, T(0));
            T n1 = std::hypot(dIn.x, dIn.y), n2 = std::hypot(dOut.x, dOut.y);
            if (n1 <= T(1e-14) || n2 <= T(1e-14)) continue;
            dIn.x /= n1; dIn.y /= n1;
            dOut.x /= n2; dOut.y /= n2;
            T cr = dIn.x * dOut.y - dIn.y * dOut.x;
            T dt = dIn.x * dOut.x + dIn.y * dOut.y;
            T turnDeg = std::atan2(cr, dt) * radToDeg; // >0 left turn (convex), <0 right turn (reflex)
            T intAngDeg = T(180) - turnDeg;
            if (intAngDeg <= T(0)) intAngDeg += T(360);
            if (intAngDeg > T(360)) intAngDeg -= T(360);

            if (std::abs(intAngDeg - T(180)) >= sharpThreshDeg) {
                typename VirtualCADLoop<T>::PinnedCorner pc;
                pc.arcS = vl.pieces[c].s0;
                pc.pieceIdx = c;
                pc.pt = vl.pieces[c].seg.p0();
                pc.interiorAngleDeg = intAngDeg;
                pc.targetPatchValence = std::max(1, static_cast<int>(std::floor((intAngDeg + T(45)) / T(90))));
                vl.corners.push_back(pc);
            }
        }
    }
    return vloops;
}

// Extract exact CAD sub-curve moving FORWARD along VirtualCADLoop from arc-length sStart to sEnd
template <class T>
inline CurveSeg<T> extractVirtualSubCurve(const VirtualCADLoop<T>& vl,
                                          T sStart,
                                          T sEnd,
                                          int degree)
{
    if (vl.pieces.empty() || vl.totalArcLen <= T(1e-14))
        return CurveSeg<T>();

    const T L = vl.totalArcLen;
    T a = virtualLoopWrapS(vl, sStart);
    T b = virtualLoopWrapS(vl, sEnd);
    T span = (b > a + T(1e-10)) ? (b - a) : (b + L - a);
    if (span <= T(1e-10)) return CurveSeg<T>();

    auto slicePiece = [&](const CurveSeg<T>& seg, T u0, T u1) -> CurveSeg<T> {
        u0 = std::max(T(0), std::min(T(1), u0));
        u1 = std::max(T(0), std::min(T(1), u1));
        if (u1 <= u0 + T(1e-6)) return CurveSeg<T>();
        if (u0 <= T(1e-5) && u1 >= T(1) - T(1e-5)) return seg;
        CurveSeg<T> cur = seg;
        if (u0 > T(1e-5)) {
            cur = cur.splitAtParam(u0).second;
            u1 = (u1 - u0) / std::max(T(1e-12), T(1) - u0);
        }
        if (u1 < T(1) - T(1e-5))
            cur = cur.splitAtParam(u1).first;
        return cur;
    };

    std::vector<CurveSeg<T>> outSegs;
    T curS = a;
    T rem = span;
    const int maxSteps = static_cast<int>(vl.pieces.size()) * 2 + 4;
    for (int guard = 0; guard < maxSteps && rem > T(1e-8); ++guard) {
        T uStart = T(0);
        int pi = findVirtualPieceAtS(vl, virtualLoopWrapS(vl, curS), uStart);
        if (pi < 0) break;
        const auto& pc = vl.pieces[pi];
        T pieceLen = pc.s1 - pc.s0;
        T avail = std::max(T(0), pc.s1 - virtualLoopWrapS(vl, curS));
        if (avail <= T(1e-10)) {
            curS = virtualLoopWrapS(vl, pc.s1 + T(1e-9));
            continue;
        }
        T take = std::min(rem, avail);
        T uEnd = (take >= avail - T(1e-9) || pieceLen <= T(1e-12))
               ? T(1)
               : virtualPieceArcToParam<T>(pc, (virtualLoopWrapS(vl, curS) - pc.s0) + take);
        CurveSeg<T> sl = slicePiece(pc.seg, uStart, uEnd);
        if (sl.bsp) outSegs.push_back(sl);
        curS = virtualLoopWrapS(vl, curS + take);
        rem -= take;
    }

    if (outSegs.empty()) return CurveSeg<T>();
    if (outSegs.size() == 1) return outSegs[0];
    auto joined = joinCurveSegsExact(outSegs, degree);
    return joined ? CurveSeg<T>(*joined) : outSegs[0];
}

// Concatenate K curve segments onto exact uniform macro-parameter breakpoints i/K in [0, 1]
// so exact knot-splitting of a macro-patch at u = i/K recovers every sub-segment bit-for-bit!
template <class T>
inline CurveSeg<T> joinCurveSegsUniformParam(const std::vector<CurveSeg<T>>& segs, int degree)
{
    if (segs.empty()) return CurveSeg<T>();
    if (segs.size() == 1) {
        auto copy = std::make_shared<gsBSpline<T>>(*segs[0].bsp);
        safeNormalizeKnots01(*copy);
        if (copy->degree() < degree)
            copy->degreeElevate(degree - copy->degree());
        return CurveSeg<T>(*copy);
    }

    int d = degree;
    for (const auto& s : segs)
        if (s.bsp) d = std::max(d, s.bsp->degree());

    const int K = static_cast<int>(segs.size());
    std::vector<std::unique_ptr<gsBSpline<T>>> pieces(K);
    for (int i = 0; i < K; ++i) {
        pieces[i].reset(new gsBSpline<T>(*segs[i].bsp));
        safeNormalizeKnots01(*pieces[i]);
        if (pieces[i]->degree() < d)
            pieces[i]->degreeElevate(d - pieces[i]->degree());
    }

    std::vector<T> allKnots;
    for (int m = 0; m < d + 1; ++m) allKnots.push_back(T(0));
    int totalCoefs = 1;
    for (int i = 0; i < K; ++i) {
        const auto& kv = pieces[i]->knots();
        const T s0 = static_cast<T>(i) / static_cast<T>(K);
        const T s1 = static_cast<T>(i + 1) / static_cast<T>(K);
        const index_t kvSz = static_cast<index_t>(kv.size());
        for (index_t j = d + 1; j < kvSz - (d + 1); ++j) {
            T u = std::max(s0 + T(1e-12), std::min(s1 - T(1e-12), s0 + (s1 - s0) * kv[j]));
            allKnots.push_back(std::max(allKnots.back(), u));
        }
        if (i + 1 < K) {
            for (int m = 0; m < d; ++m) allKnots.push_back(s1);
        }
        totalCoefs += static_cast<int>(pieces[i]->coefs().rows()) - 1;
    }
    for (int m = 0; m < d + 1; ++m) allKnots.push_back(T(1));

    gsMatrix<T> coefs(totalCoefs, 2);
    coefs.row(0) = pieces[0]->coefs().row(0);
    int rowIdx = 1;
    for (int i = 0; i < K; ++i) {
        const auto& pc = pieces[i]->coefs();
        const index_t nr = pc.rows();
        for (index_t r = 1; r < nr - 1; ++r)
            coefs.row(rowIdx++) = pc.row(r);
        if (i + 1 < K) {
            const auto& nc = pieces[i + 1]->coefs();
            coefs.row(rowIdx++) = T(0.5) * (pc.row(nr - 1) + nc.row(0));
        } else {
            coefs.row(rowIdx++) = pc.row(nr - 1);
        }
    }
    gsKnotVector<T> newKv(static_cast<short_t>(d), allKnots.begin(), allKnots.end());
    return CurveSeg<T>(gsBSpline<T>(newKv, coefs));
}

// ============================================================================
// Stage 1: Curvature-Regularized Soft-Boundary (Robin) 4-RoSy + Dipole Pruning
// ============================================================================

template <class T>
inline CrossFieldSolution<T> solveRegularized4RoSy(
    const TriMeshData<T>& mesh,
    const std::vector<VirtualCADLoop<T>>& vloops,
    const CrossFieldMeshAccel<T>* accelPtr = nullptr)
{
    CrossFieldSolution<T> sol;
    const int nV = static_cast<int>(mesh.pts.size());
    sol.z.assign(nV, std::complex<T>(T(1), T(0)));
    if (nV < 3 || mesh.tris.empty()) return sol;

    CrossFieldMeshAccel<T> localAccel;
    if (!accelPtr) {
        localAccel.build(mesh);
        accelPtr = &localAccel;
    }

    T minX = mesh.pts[0].x, maxX = mesh.pts[0].x, minY = mesh.pts[0].y, maxY = mesh.pts[0].y;
    for (const auto& p : mesh.pts) {
        minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
    }
    const T diag = std::max(T(1e-9), std::hypot(maxX - minX, maxY - minY));

    T avgEdge = T(0);
    for (const auto& e : mesh.edges)
        avgEdge += mesh.pts[e.a].dist(mesh.pts[e.b]);
    avgEdge /= std::max<T>(T(1), static_cast<T>(mesh.edges.size()));

    // Compute 1D boundary lumped mass per boundary vertex
    std::vector<T> bndMass(nV, T(0));
    for (const auto& e : mesh.edges) {
        if (!e.isBoundary) continue;
        T len = mesh.pts[e.a].dist(mesh.pts[e.b]);
        bndMass[e.a] += T(0.5) * len;
        bndMass[e.b] += T(0.5) * len;
    }

    // Gaussian arc-length smoothed boundary tangent on smooth spans (eliminates wavy-boundary dipoles!)
    auto evalSmoothedBndZ = [&](int vIdx, bool& isPinnedCornerVert) -> std::complex<T> {
        isPinnedCornerVert = false;
        std::complex<T> t1(mesh.bndTanIn[vIdx].x, mesh.bndTanIn[vIdx].y);
        std::complex<T> t2(mesh.bndTanOut[vIdx].x, mesh.bndTanOut[vIdx].y);
        std::complex<T> zRaw = std::pow(t1, 4) + std::pow(t2, 4);
        T nRaw = std::abs(zRaw);
        zRaw = (nRaw > T(1e-12)) ? (zRaw / nRaw) : std::pow(t2, 4);

        if (vIdx >= static_cast<int>(mesh.bndInfo.size())) return zRaw;
        int lIdx = mesh.bndInfo[vIdx].loopIdx;
        if (lIdx < 0 || lIdx >= static_cast<int>(vloops.size()) || vloops[lIdx].pieces.empty())
            return zRaw;

        const auto& vl = vloops[lIdx];
        T s0 = T(0);
        gsPoint2D<T> proj;
        projectPointToVirtualLoop(vl, mesh.pts[vIdx], s0, proj);

        // Distance to nearest pinned sharp corner along the loop
        T minCornerDist = vl.totalArcLen;
        for (const auto& c : vl.corners) {
            T ds = std::abs(s0 - c.arcS);
            ds = std::min(ds, vl.totalArcLen - ds);
            minCornerDist = std::min(minCornerDist, ds);
        }
        if (minCornerDist < T(1.5) * avgEdge) {
            isPinnedCornerVert = true;
            return zRaw;
        }

        // Smooth tangent in complex 4-RoSy space along the smooth span (never crossing a pinned corner!)
        T sigma = std::min(T(0.075) * diag, T(0.45) * minCornerDist);
        if (sigma <= T(0.8) * avgEdge) return zRaw;

        std::complex<T> zSum(T(0), T(0));
        const int nK = 9;
        for (int k = -nK; k <= nK; ++k) {
            T ds = (static_cast<T>(k) / static_cast<T>(nK)) * T(2.0) * sigma;
            if (std::abs(ds) >= minCornerDist * T(0.95)) continue;
            T w = std::exp(-T(0.5) * (ds / sigma) * (ds / sigma));
            gsPoint2D<T> tk = evalVirtualLoopUnitTangent(vl, s0 + ds);
            std::complex<T> zk = std::pow(std::complex<T>(tk.x, tk.y), 4);
            zSum += w * zk;
        }
        T nrm = std::abs(zSum);
        return (nrm > T(1e-10)) ? (zSum / nrm) : zRaw;
    };

    std::vector<bool> isDirichlet(nV, false);
    std::vector<bool> isRobin(nV, false);
    std::vector<std::complex<T>> bndTarget(nV, std::complex<T>(T(1), T(0)));

    for (int i = 0; i < nV; ++i) {
        if (!mesh.isBoundaryPt[i]) continue;
        bool nearPinned = false;
        std::complex<T> zTgt = evalSmoothedBndZ(i, nearPinned);
        bndTarget[i] = zTgt;

        std::complex<T> t1(mesh.bndTanIn[i].x, mesh.bndTanIn[i].y);
        std::complex<T> t2(mesh.bndTanOut[i].x, mesh.bndTanOut[i].y);
        std::complex<T> z1 = std::pow(t1, 4), z2 = std::pow(t2, 4);

        if (nearPinned) {
            // Hard Dirichlet only at near-90/180/270 corners; free non-orthogonal sharp corners
            if (std::abs(z1 - z2) < T(0.45))
                isDirichlet[i] = true;
            else
                isRobin[i] = true;
        } else {
            // Smooth boundary interior: soft Robin alignment allows cubic splines to absorb boundary waves!
            isRobin[i] = true;
        }
    }

    std::vector<int> freeIdx(nV, -1);
    int nFree = 0;
    for (int i = 0; i < nV; ++i)
        if (!isDirichlet[i]) freeIdx[i] = nFree++;

    const T betaBnd = T(14.0) / diag;

    if (nFree == 0) {
        sol.z = bndTarget;
    } else {
        gsSparseEntries<T> lapEntries;
        lapEntries.reserve(mesh.tris.size() * 12 + nV);
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
                if (!isDirichlet[t.v[k]])
                    massDiag(freeIdx[t.v[k]]) += mTri;

            auto cotAngle = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C) -> T {
                T dot = (B.x - A.x) * (C.x - A.x) + (B.y - A.y) * (C.y - A.y);
                return std::max(T(1e-4), dot / area2);
            };
            T w12 = T(0.5) * cotAngle(p0, p1, p2);
            T w20 = T(0.5) * cotAngle(p1, p2, p0);
            T w01 = T(0.5) * cotAngle(p2, p0, p1);

            addWeight(t.v[1], t.v[2], w12); addWeight(t.v[2], t.v[1], w12);
            addWeight(t.v[2], t.v[0], w20); addWeight(t.v[0], t.v[2], w20);
            addWeight(t.v[0], t.v[1], w01); addWeight(t.v[1], t.v[0], w01);
        }

        // Add Robin boundary alignment penalty betaBnd * M_bnd * (z - z_smooth)
        for (int i = 0; i < nV; ++i) {
            if (!isRobin[i] || isDirichlet[i]) continue;
            int iu = freeIdx[i];
            T wr = betaBnd * std::max(T(1e-6) * avgEdge, bndMass[i]);
            lapEntries.add(iu, iu, wr);
            rhsLap(iu, 0) += wr * bndTarget[i].real();
            rhsLap(iu, 1) += wr * bndTarget[i].imag();
        }

        gsSparseMatrix<T> L(nFree, nFree);
        L.setFrom(lapEntries);
        typename gsSparseSolver<T>::SimplicialLDLT lapSolver;
        lapSolver.compute(L);
        gsMatrix<T> solMat = lapSolver.solve(rhsLap);

        // MBO Unit-Norm Diffusion-Projection (Viertel & Osting)
        const T tau = T(2.8) * avgEdge * avgEdge;
        if (tau > T(1e-14)) {
            gsSparseMatrix<T> Mmat(nFree, nFree);
            {
                gsSparseEntries<T> mEnt;
                mEnt.reserve(nFree);
                for (int i = 0; i < nFree; ++i)
                    mEnt.add(i, i, std::max(T(1e-12), massDiag(i)));
                Mmat.setFrom(mEnt);
            }
            gsSparseMatrix<T> A = Mmat + tau * L;
            typename gsSparseSolver<T>::SimplicialLDLT mboSolver;
            mboSolver.compute(A);
            if (mboSolver.info() == gsEigen::Success) {
                const gsMatrix<T> bndForce = tau * rhsLap;
                gsMatrix<T> rhsMbo(nFree, 2);
                for (int it = 0; it < 10; ++it) {
                    for (int i = 0; i < nFree; ++i) {
                        T zr = solMat(i, 0), zi = solMat(i, 1);
                        T nrm = std::hypot(zr, zi);
                        if (nrm > T(1e-10)) { zr /= nrm; zi /= nrm; }
                        T mi = std::max(T(1e-12), massDiag(i));
                        rhsMbo(i, 0) = mi * zr + bndForce(i, 0);
                        rhsMbo(i, 1) = mi * zi + bndForce(i, 1);
                    }
                    solMat = mboSolver.solve(rhsMbo);
                }
            }
        }

        for (int i = 0; i < nV; ++i) {
            if (isDirichlet[i]) {
                sol.z[i] = bndTarget[i];
            } else {
                int ii = freeIdx[i];
                sol.z[i] = std::complex<T>(solMat(ii, 0), solMat(ii, 1));
            }
        }
    }

    // Singularity extraction + extended dipole cancellation + boundary-skin medial repulsion
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
        std::complex<T> z0 = sol.z[t.v[0]], z1 = sol.z[t.v[1]], z2 = sol.z[t.v[2]];
        if (std::abs(z0) < T(1e-12) || std::abs(z1) < T(1e-12) || std::abs(z2) < T(1e-12))
            continue;
        T a0 = std::arg(z0), a1 = std::arg(z1), a2 = std::arg(z2);
        T sum = wrapPi(a1 - a0) + wrapPi(a2 - a1) + wrapPi(a0 - a2);
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
        gsPoint2D<T> pos(l0 * P0.x + l1 * P1.x + l2 * P2.x,
                         l0 * P0.y + l1 * P1.y + l2 * P2.y);
        if (isPointStrictlyInsideDomain(pos, mesh.loopPolys, T(1e-5)))
            rawSings.push_back({pos, ind, ti});
    }

    // Cluster nearby singularities
    const T clusterTol = T(2.0) * avgEdge;
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
                sumInd += rawSings[j].idx;
                ++cnt;
            }
        }
        if (sumInd != 0)
            clustered.push_back({gsPoint2D<T>(sumP.x / cnt, sumP.y / cnt), sumInd, rawSings[i].triIdx});
    }

    // Extended opposite-sign dipole cancellation (closest visible +1/-1 pairs first!)
    const T dipoleTol = std::max(T(3.5) * avgEdge, T(0.20) * diag);
    struct DipPair { T d; size_t i; size_t j; };
    std::vector<DipPair> dipCands;
    for (size_t i = 0; i < clustered.size(); ++i) {
        for (size_t j = i + 1; j < clustered.size(); ++j) {
            if (clustered[i].idx + clustered[j].idx != 0) continue;
            T d = clustered[i].pos.dist(clustered[j].pos);
            if (d < dipoleTol && isSegmentInsideDomain(clustered[i].pos, clustered[j].pos, mesh.loopPolys))
                dipCands.push_back({d, i, j});
        }
    }
    std::sort(dipCands.begin(), dipCands.end(), [](const DipPair& a, const DipPair& b) { return a.d < b.d; });
    std::vector<bool> dead(clustered.size(), false);
    for (const auto& dp : dipCands) {
        if (!dead[dp.i] && !dead[dp.j]) { dead[dp.i] = true; dead[dp.j] = true; }
    }

    // Helper: distance to nearest boundary point & outward vector
    auto closestBndDistAndPt = [&](const gsPoint2D<T>& q, gsPoint2D<T>& bPt) -> T {
        T bestD = std::numeric_limits<T>::max();
        bPt = q;
        for (const auto& poly : mesh.loopPolys) {
            int np = static_cast<int>(poly.size());
            for (int k = 0; k < np; ++k) {
                const auto& A = poly[k];
                const auto& B = poly[(k + 1) % np];
                T dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy;
                T t = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((q.x - A.x) * dx + (q.y - A.y) * dy) / L2)) : T(0);
                gsPoint2D<T> c(A.x + t * dx, A.y + t * dy);
                T d = q.dist(c);
                if (d < bestD) { bestD = d; bPt = c; }
            }
        }
        return bestD;
    };

    // Boundary-skin medial repulsion + circular-harmonic port extraction
    const T minSkinClear = T(0.05) * diag;
    for (size_t i = 0; i < clustered.size(); ++i) {
        if (dead[i]) continue;
        typename CrossFieldSolution<T>::Singularity sg;
        sg.pos = clustered[i].pos;
        sg.index = clustered[i].idx;

        gsPoint2D<T> bPt;
        T dBnd = closestBndDistAndPt(sg.pos, bPt);
        if (dBnd < minSkinClear && dBnd > T(1e-9)) {
            gsPoint2D<T> inDir((sg.pos.x - bPt.x) / dBnd, (sg.pos.y - bPt.y) / dBnd);
            T bestClear = dBnd;
            gsPoint2D<T> bestP = sg.pos;
            for (T step : {T(0.02), T(0.035), T(0.05), T(0.065)}) {
                gsPoint2D<T> cand(sg.pos.x + step * diag * inDir.x, sg.pos.y + step * diag * inDir.y);
                if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-5)) ||
                    !isSegmentInsideDomain(sg.pos, cand, mesh.loopPolys)) break;
                gsPoint2D<T> tmpB;
                T dc = closestBndDistAndPt(cand, tmpB);
                if (dc > bestClear) { bestClear = dc; bestP = cand; }
            }
            sg.pos = bestP;
        }

        int val = std::max(2, std::min(8, 4 - sg.index));
        if (sg.index < 4) {
            const int nAng = 48;
            T probeR = T(1.5) * avgEdge;
            std::complex<T> harmSum(T(0), T(0));
            int nValid = 0;
            int triHint = clustered[i].triIdx;
            for (T rScale : {T(1.0), T(0.55), T(0.30)}) {
                harmSum = std::complex<T>(T(0), T(0));
                nValid = 0;
                T r = probeR * rScale;
                for (int m = 0; m < nAng; ++m) {
                    T phi = (twoPi * m) / static_cast<T>(nAng);
                    gsPoint2D<T> q(sg.pos.x + r * std::cos(phi), sg.pos.y + r * std::sin(phi));
                    std::complex<T> zq;
                    T thq = T(0);
                    if (evalCrossFieldAt(mesh, sol, *accelPtr, q, zq, thq, triHint)) {
                        T nrm = std::abs(zq);
                        if (nrm > T(1e-10)) {
                            std::complex<T> unitZ = zq / nrm;
                            std::complex<T> rot(std::cos(-sg.index * phi), std::sin(-sg.index * phi));
                            harmSum += unitZ * rot;
                            ++nValid;
                        }
                    }
                }
                if (nValid >= 12 && std::abs(harmSum) > T(1e-6)) break;
            }
            T phi0 = (nValid > 0 && std::abs(harmSum) > T(1e-8)) ? std::arg(harmSum) : T(0);
            sg.portAngles.reserve(val);
            for (int p = 0; p < val; ++p)
                sg.portAngles.push_back((phi0 + twoPi * static_cast<T>(p)) / static_cast<T>(val));
        }
        sol.singularities.push_back(std::move(sg));
    }

    return sol;
}

// ============================================================================
// Stage 2A: Poisson-Torsion Curved 5-Patch O-Grid (Smooth H=0 Non-Star Domains)
// ============================================================================

template <class T>
inline bool buildHarmonicCurvedOGrid(
    const gsPlanarDomain<T>& dom,
    const std::vector<VirtualCADLoop<T>>& vloops,
    const TriMeshData<T>& fineMesh,
    const CrossFieldMeshAccel<T>& accel,
    int degree,
    gsMultiPatch<T>& outMp)
{
    if (dom.numHoles() != 0 || vloops.empty() || !vloops[0].corners.empty())
        return false;

    const auto& vl = vloops[0];
    if (vl.pieces.empty() || vl.totalArcLen <= T(1e-12)) return false;

    const int nV = static_cast<int>(fineMesh.pts.size());
    const int nT = static_cast<int>(fineMesh.tris.size());
    if (nV < 10 || nT < 8 || fineMesh.loopPolys.empty()) return false;

    T minX = fineMesh.pts[0].x, maxX = fineMesh.pts[0].x, minY = fineMesh.pts[0].y, maxY = fineMesh.pts[0].y;
    for (const auto& p : fineMesh.pts) {
        minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
    }
    const T diag = std::max(T(1e-9), std::hypot(maxX - minX, maxY - minY));

    T avgEdge = T(0);
    for (const auto& e : fineMesh.edges)
        avgEdge += fineMesh.pts[e.a].dist(fineMesh.pts[e.b]);
    avgEdge /= std::max<T>(T(1), static_cast<T>(fineMesh.edges.size()));

    // 1. Solve smooth Poisson torsion potential: -Delta(phi) = 1 inside Omega, phi = 0 on boundary.
    //    Unlike a singular point-source Laplace pole, Poisson torsion has NO singularity, stays smooth
    //    along elongated/non-star lobes (e.g. amoeba), and on a unit disk equals phi_norm(r) = 1 - r^2!
    std::vector<int> freeIdx(nV, -1);
    int nFree = 0;
    for (int i = 0; i < nV; ++i)
        if (!fineMesh.isBoundaryPt[i]) freeIdx[i] = nFree++;
    if (nFree < 4) return false;

    gsSparseEntries<T> lapEnt;
    lapEnt.reserve(nT * 12);
    gsVector<T> rhs = gsVector<T>::Zero(nFree);
    auto addW = [&](int u, int v, T w) {
        if (fineMesh.isBoundaryPt[u]) return;
        int iu = freeIdx[u];
        lapEnt.add(iu, iu, w);
        if (!fineMesh.isBoundaryPt[v]) lapEnt.add(iu, freeIdx[v], -w);
    };
    for (const auto& t : fineMesh.tris) {
        const auto& p0 = fineMesh.pts[t.v[0]];
        const auto& p1 = fineMesh.pts[t.v[1]];
        const auto& p2 = fineMesh.pts[t.v[2]];
        T area2 = std::max(T(1e-10), std::abs(cross2(p0, p1, p2)));
        T mTri = area2 / T(6);
        for (int k = 0; k < 3; ++k)
            if (!fineMesh.isBoundaryPt[t.v[k]])
                rhs(freeIdx[t.v[k]]) += mTri;

        auto cotA = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C) -> T {
            return std::max(T(1e-4), ((B.x - A.x) * (C.x - A.x) + (B.y - A.y) * (C.y - A.y)) / area2);
        };
        T w12 = T(0.5) * cotA(p0, p1, p2), w20 = T(0.5) * cotA(p1, p2, p0), w01 = T(0.5) * cotA(p2, p0, p1);
        addW(t.v[1], t.v[2], w12); addW(t.v[2], t.v[1], w12);
        addW(t.v[2], t.v[0], w20); addW(t.v[0], t.v[2], w20);
        addW(t.v[0], t.v[1], w01); addW(t.v[1], t.v[0], w01);
    }

    gsSparseMatrix<T> L(nFree, nFree);
    L.setFrom(lapEnt);
    typename gsSparseSolver<T>::SimplicialLDLT solver;
    solver.compute(L);
    if (solver.info() != gsEigen::Success) return false;
    gsVector<T> phiFree = solver.solve(rhs);

    T phiMax = T(0);
    for (int i = 0; i < nFree; ++i)
        phiMax = std::max(phiMax, phiFree(i));
    if (phiMax <= T(1e-14)) return false;

    std::vector<T> phi(nV, T(0));
    for (int i = 0; i < nV; ++i)
        if (!fineMesh.isBoundaryPt[i])
            phi[i] = std::max(T(0), std::min(T(1), phiFree(freeIdx[i]) / phiMax));

    // Area-weighted vertex gradient gV = grad(phi) for smooth C0 RK2 streamline tracing
    std::vector<gsPoint2D<T>> gradV(nV, gsPoint2D<T>(0, 0));
    std::vector<T> areaV(nV, T(0));
    for (int ti = 0; ti < nT; ++ti) {
        const auto& t = fineMesh.tris[ti];
        const auto& p0 = fineMesh.pts[t.v[0]];
        const auto& p1 = fineMesh.pts[t.v[1]];
        const auto& p2 = fineMesh.pts[t.v[2]];
        T a2 = cross2(p0, p1, p2);
        if (std::abs(a2) <= T(1e-14)) continue;
        T gx = (phi[t.v[0]] * (p1.y - p2.y) + phi[t.v[1]] * (p2.y - p0.y) + phi[t.v[2]] * (p0.y - p1.y)) / a2;
        T gy = (phi[t.v[0]] * (p2.x - p1.x) + phi[t.v[1]] * (p0.x - p2.x) + phi[t.v[2]] * (p1.x - p0.x)) / a2;
        T wA = std::abs(a2);
        for (int k = 0; k < 3; ++k) {
            gradV[t.v[k]].x += wA * gx;
            gradV[t.v[k]].y += wA * gy;
            areaV[t.v[k]] += wA;
        }
    }
    for (int i = 0; i < nV; ++i) {
        if (areaV[i] > T(1e-14)) { gradV[i].x /= areaV[i]; gradV[i].y /= areaV[i]; }
    }

    // Fast triangle-walk evaluator for phi and grad(phi)
    auto evalPhiAndGrad = [&](const gsPoint2D<T>& p, T& outPhi, gsPoint2D<T>& outG, int& triHint) -> bool {
        auto baryEval = [&](int ti, T& l0, T& l1, T& l2) -> bool {
            if (ti < 0 || ti >= nT) return false;
            const auto& t = fineMesh.tris[ti];
            T a2 = cross2(fineMesh.pts[t.v[0]], fineMesh.pts[t.v[1]], fineMesh.pts[t.v[2]]);
            if (std::abs(a2) < T(1e-14)) return false;
            l0 = cross2(fineMesh.pts[t.v[1]], fineMesh.pts[t.v[2]], p) / a2;
            l1 = cross2(fineMesh.pts[t.v[2]], fineMesh.pts[t.v[0]], p) / a2;
            l2 = T(1) - l0 - l1;
            return true;
        };
        auto finish = [&](int ti, T l0, T l1, T l2) -> bool {
            const auto& t = fineMesh.tris[ti];
            outPhi = l0 * phi[t.v[0]] + l1 * phi[t.v[1]] + l2 * phi[t.v[2]];
            outG = gsPoint2D<T>(l0 * gradV[t.v[0]].x + l1 * gradV[t.v[1]].x + l2 * gradV[t.v[2]].x,
                                l0 * gradV[t.v[0]].y + l1 * gradV[t.v[1]].y + l2 * gradV[t.v[2]].y);
            triHint = ti;
            return true;
        };
        if (triHint >= 0 && triHint < nT) {
            int cur = triHint, prev = -1;
            for (int st = 0; st < 64; ++st) {
                T l0 = 0, l1 = 0, l2 = 0;
                if (!baryEval(cur, l0, l1, l2)) break;
                if (l0 >= -T(1e-4) && l1 >= -T(1e-4) && l2 >= -T(1e-4)) return finish(cur, l0, l1, l2);
                int eo[3] = {1, 2, 0}; T vs[3] = {l0, l1, l2};
                if (vs[1] < vs[0]) { std::swap(vs[0], vs[1]); std::swap(eo[0], eo[1]); }
                if (vs[2] < vs[1]) { std::swap(vs[1], vs[2]); std::swap(eo[1], eo[2]); }
                if (vs[1] < vs[0]) { std::swap(vs[0], vs[1]); std::swap(eo[0], eo[1]); }
                int nxt = -1;
                for (int m = 0; m < 3; ++m) {
                    if (vs[m] >= -T(1e-4)) break;
                    int nb = accel.triAdj[cur][eo[m]];
                    if (nb >= 0 && nb != prev) { nxt = nb; break; }
                }
                if (nxt < 0) break;
                prev = cur; cur = nxt;
            }
        }
        for (int ti = 0; ti < nT; ++ti) {
            const auto& bb = accel.triBBox[ti];
            if (p.x < bb[0] || p.x > bb[1] || p.y < bb[2] || p.y > bb[3]) continue;
            T l0 = 0, l1 = 0, l2 = 0;
            if (baryEval(ti, l0, l1, l2) && l0 >= -T(1e-4) && l1 >= -T(1e-4) && l2 >= -T(1e-4))
                return finish(ti, l0, l1, l2);
        }
        return false;
    };

    // 2. Marching-triangles level-set contour extractor for phi == phiIso (guaranteed simple CCW loop!)
    auto extractClosedIsoLoop = [&](T phiIso, std::vector<gsPoint2D<T>>& outLoop) -> bool {
        outLoop.clear();
        struct Seg2 { gsPoint2D<T> a, b; bool used = false; };
        std::vector<Seg2> segs;
        for (int ti = 0; ti < nT; ++ti) {
            const auto& t = fineMesh.tris[ti];
            T v[3] = {phi[t.v[0]], phi[t.v[1]], phi[t.v[2]]};
            int above[3] = {v[0] > phiIso ? 1 : 0, v[1] > phiIso ? 1 : 0, v[2] > phiIso ? 1 : 0};
            int sum = above[0] + above[1] + above[2];
            if (sum == 0 || sum == 3) continue;

            int odd = -1;
            int target = (sum == 1) ? 1 : 0;
            for (int k = 0; k < 3; ++k) if (above[k] == target) { odd = k; break; }
            int k1 = (odd + 1) % 3, k2 = (odd + 2) % 3;
            auto interp = [&](int i, int j) -> gsPoint2D<T> {
                T dv = v[j] - v[i];
                T al = (std::abs(dv) > T(1e-14)) ? (phiIso - v[i]) / dv : T(0.5);
                al = std::max(T(0), std::min(T(1), al));
                return gsPoint2D<T>((T(1) - al) * fineMesh.pts[t.v[i]].x + al * fineMesh.pts[t.v[j]].x,
                                    (T(1) - al) * fineMesh.pts[t.v[i]].y + al * fineMesh.pts[t.v[j]].y);
            };
            gsPoint2D<T> A = interp(odd, k1);
            gsPoint2D<T> B = interp(odd, k2);
            if (A.dist(B) <= T(1e-10)) continue;
            if (sum == 1) segs.push_back({A, B, false});
            else          segs.push_back({B, A, false});
        }
        if (segs.size() < 4) return false;

        const T snapTol = T(1e-4) * diag;
        std::vector<gsPoint2D<T>> bestLoop;
        for (size_t s0 = 0; s0 < segs.size(); ++s0) {
            if (segs[s0].used) continue;
            segs[s0].used = true;
            std::vector<gsPoint2D<T>> chain = {segs[s0].a, segs[s0].b};
            bool closed = false;
            for (size_t step = 0; step < segs.size(); ++step) {
                if (chain.size() >= 4 && chain.back().dist(chain.front()) <= snapTol) {
                    chain.pop_back();
                    closed = true;
                    break;
                }
                int bestNext = -1;
                T bestD = snapTol;
                for (size_t j = 0; j < segs.size(); ++j) {
                    if (segs[j].used) continue;
                    T d = chain.back().dist(segs[j].a);
                    if (d < bestD) { bestD = d; bestNext = static_cast<int>(j); }
                }
                if (bestNext < 0) break;
                segs[bestNext].used = true;
                chain.push_back(segs[bestNext].b);
            }
            if (closed && chain.size() > bestLoop.size())
                bestLoop = std::move(chain);
        }
        if (bestLoop.size() < 8) return false;

        T area2 = T(0);
        for (size_t i = 0; i < bestLoop.size(); ++i) {
            const auto& a = bestLoop[i];
            const auto& b = bestLoop[(i + 1) % bestLoop.size()];
            area2 += a.x * b.y - b.x * a.y;
        }
        if (area2 < T(0)) std::reverse(bestLoop.begin(), bestLoop.end());
        outLoop = std::move(bestLoop);
        return true;
    };

    // 3. Evaluate 5-patch curved O-grid candidates across Poisson level sets
    //    using BOTH outward -grad(phi) streamlines from the core AND inward +grad(phi) streamlines from the outer loop!
    const bool oDbg = (std::getenv("GISMO_CF_DEBUG") != nullptr);
    int dbgTraced = 0, dbgOrderOk = 0, dbgAngOk = 0, dbgDetOk = 0;
    T bestScore = -T(1e9);
    gsMultiPatch<T> bestMp;
    T bestFallbackScore = -T(1e9);
    gsMultiPatch<T> bestFallbackMp;

    const T isoLevels[5] = {T(0.68), T(0.55), T(0.42), T(0.78), T(0.28)};
    const T stepLen = T(0.25) * avgEdge;
    const auto& outerPoly = fineMesh.loopPolys[0];
    const int nOuterPoly = static_cast<int>(outerPoly.size());

    auto evalCandidate5 = [&](const std::array<T, 4>& sBnd,
                              const std::array<T, 4>& isoS,
                              const std::array<gsPoint2D<T>, 4>& P,
                              const std::array<gsPoint2D<T>, 4>& Q,
                              std::array<std::vector<gsPoint2D<T>>, 4>& spokePolys,
                              const std::function<std::vector<gsPoint2D<T>>(T, T)>& extractIsoSubPoly,
                              T isoTot) {
        ++dbgTraced;
        // Verify cyclic CCW order on both outer loop and core loop
        T sumB = T(0), sumI = T(0);
        for (int k = 0; k < 4; ++k) {
            T db = sBnd[(k + 1) % 4] - sBnd[k];
            if (db <= T(0)) db += vl.totalArcLen;
            if (db < T(0.05) * vl.totalArcLen || db > T(0.72) * vl.totalArcLen) return;
            sumB += db;

            T di = isoS[(k + 1) % 4] - isoS[k];
            if (di <= T(0)) di += isoTot;
            if (di < T(0.05) * isoTot || di > T(0.72) * isoTot) return;
            sumI += di;
        }
        if (std::abs(sumB - vl.totalArcLen) > T(0.05) * vl.totalArcLen ||
            std::abs(sumI - isoTot) > T(0.05) * isoTot) return;
        ++dbgOrderOk;

        CurveSeg<T> bEdge[4], cCore[4], spoke[4];
        for (int k = 0; k < 4; ++k) {
            bEdge[k] = extractVirtualSubCurve(vl, sBnd[k], sBnd[(k + 1) % 4], degree);
            if (!bEdge[k].bsp) return;
            std::vector<gsPoint2D<T>> corePoly = extractIsoSubPoly(isoS[k], isoS[(k + 1) % 4]);
            corePoly.front() = Q[k];
            corePoly.back() = Q[(k + 1) % 4];
            cCore[k] = fitStreamlinePolylineToCurveSeg(corePoly, degree, fineMesh.loopPolys);
            spokePolys[k].front() = Q[k];
            spokePolys[k].back() = P[k];
            spoke[k] = fitStreamlinePolylineToCurveSeg(spokePolys[k], degree, fineMesh.loopPolys);
        }

        gsMultiPatch<T> candMp;
        auto addSafePatch = [&](const CurveSeg<T>& e0, const CurveSeg<T>& e1,
                                const CurveSeg<T>& e2, const CurveSeg<T>& e3,
                                bool refineU, bool refineV) -> bool {
            gsTensorBSpline<2, T>* tb = makeCoonsFromCCWCycle(e0, e1, e2, e3, degree);
            if (!tb) return false;
            // Insert conforming interior knots along coarse directions so harmonic/Spring maps have real interior DOFs
            // to untangle strongly curved amoeba lobes (e.g., 35x4 -> 35x7)!
            if (refineU && tb->knots(0).unique().size() < 5) {
                for (T u : {T(0.25), T(0.5), T(0.75)}) tb->insertKnot(u, 0, 1);
            }
            if (refineV && tb->knots(1).unique().size() < 5) {
                for (T v : {T(0.25), T(0.5), T(0.75)}) tb->insertKnot(v, 1, 1);
            }
            if (evalPatchMinDetDense(*tb, 33) <= T(1e-7)) {
                quickHarmonicPatchRepair(*tb, T(1e-7));
                if (evalPatchMinDetDense(*tb, 33) <= T(1e-7)) {
                    try {
                        gsMultiPatch<T> bndC;
                        bndC.addPatch(tb->boundary(boundary::south));
                        bndC.addPatch(tb->boundary(boundary::east));
                        bndC.addPatch(tb->boundary(boundary::north));
                        bndC.addPatch(tb->boundary(boundary::west));
                        gsSpringPatch<T> sp(bndC);
                        const gsGeometry<T>& res = sp.compute();
                        if (const auto* sTb = dynamic_cast<const gsTensorBSpline<2, T>*>(&res)) {
                            if (evalPatchMinDetDense(*sTb, 33) > evalPatchMinDetDense(*tb, 33))
                                *tb = *sTb;
                        }
                    } catch (...) {}
                }
            }
            candMp.addPatch(typename gsGeometry<T>::uPtr(tb));
            return true;
        };

        bool ok5 = addSafePatch(cCore[0], cCore[1], cCore[2], cCore[3], true, true);
        for (int k = 0; k < 4 && ok5; ++k)
            ok5 = addSafePatch(bEdge[k], spoke[(k + 1) % 4].reversed(), cCore[k].reversed(), spoke[k], true, true);
        if (!ok5 || candMp.nPatches() != 5) return;

        T minDetAll = std::numeric_limits<T>::max();
        T minEfAll = T(180);
        int nNeg = 0;
        for (size_t p = 0; p < 5; ++p) {
            T angs[4], mnA = 0, mxA = 180; bool deg = false;
            evalPatchCornerAnglesDeg(candMp.patch(p), angs, mnA, mxA, deg);
            if (deg || mnA < T(8.0) || mxA > T(168.0)) return;
            minEfAll = std::min(minEfAll, std::min(mnA, T(180) - mxA));
            T d49 = evalPatchMinDetDense(candMp.patch(p), 49);
            minDetAll = std::min(minDetAll, d49);
            if (d49 <= T(1e-8)) ++nNeg;
        }
        ++dbgAngOk;
        T maxER = evalMultiPatchMaxEdgeRatio(candMp);
        if (nNeg == 0) {
            ++dbgDetOk;
            T sc = minEfAll + T(4.0) * std::log(std::max(T(1e-9), minDetAll)) - T(4.0) * std::max(T(0), maxER - T(3.2));
            if (sc > bestScore) { bestScore = sc; bestMp = give(candMp); }
        } else {
            // Keep best near-valid 5-patch candidate so evaluateAndRegister's Winslow + coupled relaxation can untangle it!
            T fsc = minEfAll - T(10.0) * static_cast<T>(nNeg) + T(20.0) * minDetAll - T(2.0) * std::max(T(0), maxER - T(4.0));
            if (fsc > bestFallbackScore) { bestFallbackScore = fsc; bestFallbackMp = give(candMp); }
        }
    };

    for (T phiIso : isoLevels) {
        std::vector<gsPoint2D<T>> isoLoop;
        if (!extractClosedIsoLoop(phiIso, isoLoop)) continue;

        const size_t M = isoLoop.size();
        std::vector<T> isoCum(M + 1, T(0));
        for (size_t i = 0; i < M; ++i)
            isoCum[i + 1] = isoCum[i] + isoLoop[i].dist(isoLoop[(i + 1) % M]);
        const T isoTot = isoCum[M];
        if (isoTot <= T(1e-8)) continue;

        auto sampleIsoAt = [&](T s) -> gsPoint2D<T> {
            T sw = std::fmod(s, isoTot);
            if (sw < T(0)) sw += isoTot;
            for (size_t i = 0; i < M; ++i) {
                if (isoCum[i + 1] >= sw || i + 1 == M) {
                    T ds = isoCum[i + 1] - isoCum[i];
                    T al = (ds > T(1e-14)) ? (sw - isoCum[i]) / ds : T(0);
                    const auto& A = isoLoop[i];
                    const auto& B = isoLoop[(i + 1) % M];
                    return gsPoint2D<T>((T(1) - al) * A.x + al * B.x, (T(1) - al) * A.y + al * B.y);
                }
            }
            return isoLoop[0];
        };

        auto projectToIsoLoop = [&](const gsPoint2D<T>& q, T& outIsoS, gsPoint2D<T>& outPt) {
            T bestD = std::numeric_limits<T>::max();
            outIsoS = T(0); outPt = isoLoop[0];
            for (size_t i = 0; i < M; ++i) {
                const auto& A = isoLoop[i];
                const auto& B = isoLoop[(i + 1) % M];
                T dx = B.x - A.x, dy = B.y - A.y, L2 = dx * dx + dy * dy;
                T t = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((q.x - A.x) * dx + (q.y - A.y) * dy) / L2)) : T(0);
                gsPoint2D<T> c(A.x + t * dx, A.y + t * dy);
                T d = q.dist(c);
                if (d < bestD) { bestD = d; outPt = c; outIsoS = isoCum[i] + t * (isoCum[i + 1] - isoCum[i]); }
            }
        };

        auto extractIsoSubPoly = [&](T sA, T sB) -> std::vector<gsPoint2D<T>> {
            T len = (sB > sA + T(1e-10)) ? (sB - sA) : (sB + isoTot - sA);
            std::vector<gsPoint2D<T>> sub;
            const int nSub = 16;
            sub.reserve(nSub + 1);
            for (int k = 0; k <= nSub; ++k)
                sub.push_back(sampleIsoAt(sA + len * (static_cast<T>(k) / static_cast<T>(nSub))));
            return sub;
        };

        // Strategy A: Trace OUTWARD along -grad(phi) from 4 equal-arc-length points on the core loop!
        // (Harmonic/Poisson outward streamlines naturally hit concave bays and separate convex lobes!)
        for (int ph = 0; ph < 6; ++ph) {
            std::array<T, 4> isoS{}, sBnd{};
            std::array<gsPoint2D<T>, 4> P{}, Q{};
            std::array<std::vector<gsPoint2D<T>>, 4> spokePolys;
            bool traceOk = true;
            for (int k = 0; k < 4 && traceOk; ++k) {
                isoS[k] = std::fmod((static_cast<T>(ph) / T(24) + static_cast<T>(k) * T(0.25)) * isoTot, isoTot);
                Q[k] = sampleIsoAt(isoS[k]);
                gsPoint2D<T> cur = Q[k];
                std::vector<gsPoint2D<T>> outPts = {Q[k]};
                int triHint = -1;
                bool hitBnd = false;
                for (int st = 0; st < 400; ++st) {
                    T phVal = T(0);
                    gsPoint2D<T> g1(0, 0);
                    if (!evalPhiAndGrad(cur, phVal, g1, triHint)) break;
                    T gn = std::hypot(g1.x, g1.y);
                    if (gn <= T(1e-12)) break;
                    // Outward direction = -grad(phi)
                    gsPoint2D<T> d1(-g1.x / gn, -g1.y / gn);
                    gsPoint2D<T> mid(cur.x + T(0.5) * stepLen * d1.x, cur.y + T(0.5) * stepLen * d1.y);
                    T phMid = T(0); gsPoint2D<T> g2;
                    gsPoint2D<T> d2 = d1;
                    if (evalPhiAndGrad(mid, phMid, g2, triHint)) {
                        T gn2 = std::hypot(g2.x, g2.y);
                        if (gn2 > T(1e-12)) d2 = gsPoint2D<T>(-g2.x / gn2, -g2.y / gn2);
                    }
                    gsPoint2D<T> nxt(cur.x + stepLen * d2.x, cur.y + stepLen * d2.y);

                    // Check boundary intersection [cur, nxt]
                    T bestT = T(2);
                    gsPoint2D<T> hitPt = nxt;
                    for (int bi = 0; bi < nOuterPoly; ++bi) {
                        const auto& A = outerPoly[bi];
                        const auto& B = outerPoly[(bi + 1) % nOuterPoly];
                        T rx = nxt.x - cur.x, ry = nxt.y - cur.y;
                        T sx = B.x - A.x, sy = B.y - A.y;
                        T den = rx * sy - ry * sx;
                        if (std::abs(den) < T(1e-14)) continue;
                        T t = ((A.x - cur.x) * sy - (A.y - cur.y) * sx) / den;
                        T u = ((A.x - cur.x) * ry - (A.y - cur.y) * rx) / den;
                        if (t > T(1e-4) && t <= T(1) + T(1e-4) && u >= -T(1e-4) && u <= T(1) + T(1e-4)) {
                            if (t < bestT) { bestT = t; hitPt = gsPoint2D<T>(cur.x + t * rx, cur.y + t * ry); }
                        }
                    }
                    if (bestT <= T(1) + T(1e-4) || !isPointStrictlyInsideDomain(nxt, fineMesh.loopPolys, T(1e-6))) {
                        projectPointToVirtualLoop(vl, hitPt, sBnd[k], P[k]);
                        outPts.push_back(P[k]);
                        hitBnd = true;
                        break;
                    }
                    outPts.push_back(nxt);
                    cur = nxt;
                }
                if (!hitBnd) { traceOk = false; break; }
                spokePolys[k] = std::move(outPts);
            }
            if (traceOk)
                evalCandidate5(sBnd, isoS, P, Q, spokePolys, extractIsoSubPoly, isoTot);
        }

        // Strategy B: Trace INWARD along +grad(phi) from 4 uniform outer boundary anchors
        for (int ph = 0; ph < 4; ++ph) {
            std::array<T, 4> isoS{}, sBnd{};
            std::array<gsPoint2D<T>, 4> P{}, Q{};
            std::array<std::vector<gsPoint2D<T>>, 4> spokePolys;
            bool traceOk = true;
            for (int k = 0; k < 4 && traceOk; ++k) {
                sBnd[k] = virtualLoopWrapS(vl, (static_cast<T>(ph) * T(0.0625) + static_cast<T>(k) * T(0.25)) * vl.totalArcLen);
                P[k] = evalVirtualLoopPoint(vl, sBnd[k]);
                gsPoint2D<T> tanB = evalVirtualLoopUnitTangent(vl, sBnd[k]);
                gsPoint2D<T> inN(-tanB.y, tanB.x);
                gsPoint2D<T> cur(P[k].x + T(0.15) * avgEdge * inN.x, P[k].y + T(0.15) * avgEdge * inN.y);
                if (!isPointStrictlyInsideDomain(cur, fineMesh.loopPolys, T(1e-6)))
                    cur = gsPoint2D<T>(P[k].x + T(0.05) * avgEdge * inN.x, P[k].y + T(0.05) * avgEdge * inN.y);

                std::vector<gsPoint2D<T>> inwardPts = {P[k]};
                int triHint = -1;
                bool hitIso = false;
                for (int st = 0; st < 400; ++st) {
                    T phVal = T(0);
                    gsPoint2D<T> g1(0, 0);
                    if (!evalPhiAndGrad(cur, phVal, g1, triHint)) break;
                    if (phVal >= phiIso) {
                        projectToIsoLoop(cur, isoS[k], Q[k]);
                        inwardPts.push_back(Q[k]);
                        hitIso = true;
                        break;
                    }
                    T gn = std::hypot(g1.x, g1.y);
                    if (gn <= T(1e-12)) break;
                    g1.x /= gn; g1.y /= gn;
                    gsPoint2D<T> mid(cur.x + T(0.5) * stepLen * g1.x, cur.y + T(0.5) * stepLen * g1.y);
                    T phMid = T(0); gsPoint2D<T> g2 = g1;
                    if (evalPhiAndGrad(mid, phMid, g2, triHint)) {
                        T gn2 = std::hypot(g2.x, g2.y);
                        if (gn2 > T(1e-12)) { g2.x /= gn2; g2.y /= gn2; } else g2 = g1;
                    }
                    gsPoint2D<T> nxt(cur.x + stepLen * g2.x, cur.y + stepLen * g2.y);
                    T phNxt = T(0); gsPoint2D<T> gTmp;
                    if (evalPhiAndGrad(nxt, phNxt, gTmp, triHint) && phNxt >= phiIso) {
                        T denom = std::max(T(1e-12), phNxt - phVal);
                        T frac = std::max(T(0), std::min(T(1), (phiIso - phVal) / denom));
                        gsPoint2D<T> hitPt(cur.x + frac * (nxt.x - cur.x), cur.y + frac * (nxt.y - cur.y));
                        projectToIsoLoop(hitPt, isoS[k], Q[k]);
                        inwardPts.push_back(Q[k]);
                        hitIso = true;
                        break;
                    }
                    inwardPts.push_back(cur);
                    cur = nxt;
                }
                if (!hitIso) { traceOk = false; break; }
                std::reverse(inwardPts.begin(), inwardPts.end());
                spokePolys[k] = std::move(inwardPts);
            }
            if (traceOk)
                evalCandidate5(sBnd, isoS, P, Q, spokePolys, extractIsoSubPoly, isoTot);
        }
    }

    if (oDbg) {
        std::cerr << "[CF_DBG] Stage2A stats: traced=" << dbgTraced << " orderOk=" << dbgOrderOk
                  << " angOk=" << dbgAngOk << " detOk=" << dbgDetOk << "\n";
    }
    if (bestMp.nPatches() == 5) {
        outMp = give(bestMp);
        return true;
    }
    if (bestFallbackMp.nPatches() == 5) {
        outMp = give(bestFallbackMp);
        return true;
    }
    return false;
}

// ============================================================================
// Stage 2B: Simultaneous Motorcycle Graph -> Rectangular T-Mesh Partitioning
// ============================================================================

template <class T>
struct TMeshGraph {
    struct Node {
        gsPoint2D<T> pt;
        bool isBoundary = false;
        int loopIdx = -1;
        T loopArcS = T(0);
        bool isPinnedCorner = false;
        bool isTJunction = false;
    };
    struct Segment {
        int u = -1, v = -1;
        CurveSeg<T> seg;
        bool isBoundary = false;
        int loopIdx = -1;
        T s0 = T(0), s1 = T(0);
        T arcLen = T(0);
        int quantLen = 1;
    };
    struct HalfEdge {
        int origin = -1, target = -1, segIdx = -1;
        bool forward = true;
        int twin = -1, next = -1;
        T angle = T(0);
    };
    struct RectFace {
        // 4 CCW sides: 0=South, 1=East, 2=North, 3=West
        std::array<std::vector<int>, 4> sideHEs;
        int corners[4] = {-1, -1, -1, -1};
        T Lu = T(1), Lv = T(1);
    };
    std::vector<Node> nodes;
    std::vector<Segment> segs;
    std::vector<HalfEdge> halfEdges;
    std::vector<RectFace> faces;
};

template <class T>
inline T measureCurveSegArcLen(const CurveSeg<T>& seg, int nS = 16)
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
inline bool buildSimultaneousMotorcycleTMesh(
    const gsPlanarDomain<T>& dom,
    const std::vector<VirtualCADLoop<T>>& vloops,
    const TriMeshData<T>& fineMesh,
    const CrossFieldSolution<T>& cf,
    const CrossFieldMeshAccel<T>& accel,
    int degree,
    int reflexMode,
    TMeshGraph<T>& outTMesh)
{
    outTMesh = TMeshGraph<T>();
    const int nLoops = static_cast<int>(vloops.size());
    if (nLoops == 0 || fineMesh.pts.empty() || fineMesh.loopPolys.empty()) return false;

    const gsMatrix<T> bbox = dom.boundingBox();
    const T diag = std::max(T(1e-9), std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0)));

    T avgEdge = T(0);
    for (const auto& e : fineMesh.edges)
        avgEdge += fineMesh.pts[e.a].dist(fineMesh.pts[e.b]);
    avgEdge = std::max(diag / T(40), avgEdge / std::max<T>(T(1), static_cast<T>(fineMesh.edges.size())));

    const T pi = T(3.14159265358979323846);
    const T radToDeg = T(180) / pi;

    struct Ray {
        int id = -1;
        int sourceId = 0;         // >0 singularity, <0 pinned corner/bridge
        int startLoop = -1;
        T startArcS = T(0);
        int endLoop = -1;
        T endArcS = T(0);
        int targetSourceId = 0;
        std::vector<gsPoint2D<T>> pts;
        gsPoint2D<T> dir;
        int triHint = -1;
        T arcLen = T(0);
        bool active = true;
        bool alive = true;
        bool terminated = false;
    };

    struct SourceInfo {
        int sourceId = 0;
        gsPoint2D<T> pt;
        int loopIdx = -1;
        T arcS = T(0);
    };

    std::vector<Ray> rays;
    std::vector<SourceInfo> allSources;
    int nextCornerSrcId = -1;

    // 1. Seed mandatory rays from VirtualCADLoop PinnedCorners
    for (int l = 0; l < nLoops; ++l) {
        const auto& vl = vloops[l];
        for (const auto& pc : vl.corners) {
            int sid = nextCornerSrcId--;
            allSources.push_back({sid, pc.pt, l, pc.arcS});
            if (pc.interiorAngleDeg <= T(152.0)) continue;

            // Domain interior is strictly to the LEFT of every VirtualCADLoop!
            gsPoint2D<T> tOut = evalVirtualLoopUnitTangent(vl, pc.arcS + T(1e-4) * vl.totalArcLen);
            T angStart = std::atan2(tOut.y, tOut.x);
            T sweepRad = pc.interiorAngleDeg * pi / T(180);

            int nRays = 1;
            if (reflexMode == 1 && pc.interiorAngleDeg <= T(285.0))
                nRays = 1;
            else
                nRays = std::max(1, pc.targetPatchValence - 1);

            for (int p = 1; p <= nRays; ++p) {
                T ang = angStart + sweepRad * (static_cast<T>(p) / static_cast<T>(nRays + 1));
                Ray r;
                r.id = static_cast<int>(rays.size());
                r.sourceId = sid;
                r.startLoop = l;
                r.startArcS = pc.arcS;
                r.pts.push_back(pc.pt);
                r.dir = gsPoint2D<T>(std::cos(ang), std::sin(ang));
                rays.push_back(std::move(r));
            }
        }
    }

    // 2. Seed rays from surviving interior singularities
    for (size_t s = 0; s < cf.singularities.size(); ++s) {
        const auto& sg = cf.singularities[s];
        int sid = static_cast<int>(s + 1);
        allSources.push_back({sid, sg.pos, -1, T(0)});
        for (T ang : sg.portAngles) {
            Ray r;
            r.id = static_cast<int>(rays.size());
            r.sourceId = sid;
            r.pts.push_back(sg.pos);
            r.dir = gsPoint2D<T>(std::cos(ang), std::sin(ang));
            rays.push_back(std::move(r));
        }
    }

    auto pickCrossBranch = [&](const gsPoint2D<T>& q, const gsPoint2D<T>& refDir, int& triHint, gsPoint2D<T>& outDir) -> bool {
        std::complex<T> zq;
        T th = T(0);
        if (!evalCrossFieldAt(fineMesh, cf, accel, q, zq, th, triHint)) return false;
        T bestDot = -T(2);
        gsPoint2D<T> bestD = refDir;
        for (int k = 0; k < 4; ++k) {
            T a = th + k * (pi * T(0.5));
            gsPoint2D<T> cand(std::cos(a), std::sin(a));
            T dt = cand.x * refDir.x + cand.y * refDir.y;
            if (dt > bestDot) { bestDot = dt; bestD = cand; }
        }
        if (bestDot > T(0.35)) { outDir = bestD; return true; }
        return false;
    };

    const T stepLen = T(0.28) * avgEdge;
    const int maxSteps = 440;
    const T maxArcLen = T(2.6) * diag;
    const T tipSnapTol = T(1.55) * avgEdge;
    const T srcSnapTol = T(1.15) * avgEdge;

    // Reusable synchronous motorcycle propagator (works on any subset of active rays!)
    auto propagateActiveRays = [&]() {
        for (int st = 0; st < maxSteps; ++st) {
            int nActive = 0;
            for (const auto& r : rays) if (r.alive && r.active) ++nActive;
            if (nActive == 0) break;

            for (size_t i = 0; i < rays.size(); ++i) {
                auto& r = rays[i];
                if (!r.alive || !r.active) continue;
                if (r.arcLen >= maxArcLen) { r.active = false; continue; }

                const gsPoint2D<T> cur = r.pts.back();
                if (r.arcLen >= T(0.85) * avgEdge) {
                    gsPoint2D<T> k1 = r.dir;
                    pickCrossBranch(cur, r.dir, r.triHint, k1);
                    gsPoint2D<T> mid(cur.x + T(0.5) * stepLen * k1.x, cur.y + T(0.5) * stepLen * k1.y);
                    gsPoint2D<T> k2 = k1;
                    if (pickCrossBranch(mid, k1, r.triHint, k2)) r.dir = k2;
                    else r.dir = k1;
                }
                gsPoint2D<T> nxt(cur.x + stepLen * r.dir.x, cur.y + stepLen * r.dir.y);

                // Rule A1: Anti-parallel active tip fusion (head-on channel capture -> single smooth separatrix!)
                if (r.arcLen >= T(0.8) * avgEdge) {
                    int fuseJ = -1;
                    T bestFD = tipSnapTol;
                    for (size_t j = 0; j < rays.size(); ++j) {
                        if (j == i || !rays[j].alive || !rays[j].active) continue;
                        if (rays[j].sourceId == r.sourceId && r.arcLen < T(2.5) * avgEdge) continue;
                        T dt = r.dir.x * rays[j].dir.x + r.dir.y * rays[j].dir.y;
                        if (dt > -T(0.68)) continue;
                        T d = cur.dist(rays[j].pts.back());
                        if (d < bestFD && isSegmentInsideDomain(cur, rays[j].pts.back(), fineMesh.loopPolys)) {
                            bestFD = d;
                            fuseJ = static_cast<int>(j);
                        }
                    }
                    if (fuseJ >= 0) {
                        auto& rj = rays[fuseJ];
                        for (int k = static_cast<int>(rj.pts.size()) - 1; k >= 0; --k)
                            r.pts.push_back(rj.pts[k]);
                        r.endLoop = rj.startLoop;
                        r.endArcS = rj.startArcS;
                        r.targetSourceId = rj.sourceId;
                        r.active = false;
                        r.terminated = true;
                        rj.active = false;
                        rj.alive = false;
                        continue;
                    }
                }

                // Rule A2: Capture nearby source node (singularity or pinned corner)
                if (r.arcLen >= T(1.1) * avgEdge) {
                    int hitSrc = -1;
                    T bestSD = srcSnapTol;
                    for (size_t si = 0; si < allSources.size(); ++si) {
                        if (allSources[si].sourceId == r.sourceId) continue;
                        T d = nxt.dist(allSources[si].pt);
                        if (d < bestSD && isSegmentInsideDomain(cur, allSources[si].pt, fineMesh.loopPolys)) {
                            bestSD = d;
                            hitSrc = static_cast<int>(si);
                        }
                    }
                    if (hitSrc >= 0) {
                        r.pts.push_back(allSources[hitSrc].pt);
                        r.endLoop = allSources[hitSrc].loopIdx;
                        r.endArcS = allSources[hitSrc].arcS;
                        r.targetSourceId = allSources[hitSrc].sourceId;
                        r.active = false;
                        r.terminated = true;
                        continue;
                    }
                }

                // Rule C: Check boundary loop crossing
                T bestBndT = T(2);
                gsPoint2D<T> bestBndHit = nxt;
                int bestBndLoop = -1;
                for (int l = 0; l < static_cast<int>(fineMesh.loopPolys.size()); ++l) {
                    const auto& poly = fineMesh.loopPolys[l];
                    const int np = static_cast<int>(poly.size());
                    for (int k = 0; k < np; ++k) {
                        const auto& A = poly[k];
                        const auto& B = poly[(k + 1) % np];
                        if (r.arcLen < T(0.6) * avgEdge && l == r.startLoop &&
                            (r.pts.front().dist(A) < T(1.2) * avgEdge || r.pts.front().dist(B) < T(1.2) * avgEdge))
                            continue;
                        T rx = nxt.x - cur.x, ry = nxt.y - cur.y;
                        T sx = B.x - A.x, sy = B.y - A.y;
                        T den = rx * sy - ry * sx;
                        if (std::abs(den) < T(1e-14)) continue;
                        T t = ((A.x - cur.x) * sy - (A.y - cur.y) * sx) / den;
                        T u = ((A.x - cur.x) * ry - (A.y - cur.y) * rx) / den;
                        if (t > T(1e-4) && t <= T(1) + T(1e-4) && u >= -T(1e-4) && u <= T(1) + T(1e-4)) {
                            if (t < bestBndT) {
                                bestBndT = t;
                                bestBndHit = gsPoint2D<T>(cur.x + t * rx, cur.y + t * ry);
                                bestBndLoop = l;
                            }
                        }
                    }
                }

                // Rule B: Stop-on-Trail collision against ALL live trails (creates T-junction, NEVER crosses!)
                T bestTrT = bestBndT;
                gsPoint2D<T> bestTrHit;
                bool hitTrail = false;
                if (r.arcLen >= T(0.5) * avgEdge) {
                    T rx = nxt.x - cur.x, ry = nxt.y - cur.y;
                    T rLen = std::hypot(rx, ry);
                    for (size_t j = 0; j < rays.size(); ++j) {
                        if (!rays[j].alive || rays[j].pts.size() < 2) continue;
                        const auto& pj = rays[j].pts;
                        size_t maxSeg = (j == i) ? (pj.size() > 9 ? pj.size() - 9 : 0) : (pj.size() - 1);
                        for (size_t k = 0; k < maxSeg; ++k) {
                            const auto& C = pj[k];
                            const auto& D = pj[k + 1];
                            if (j != i && rays[j].sourceId == r.sourceId && r.arcLen < T(2.0) * avgEdge &&
                                (C.dist(r.pts.front()) < T(1.5) * avgEdge || D.dist(r.pts.front()) < T(1.5) * avgEdge))
                                continue;
                            T sx = D.x - C.x, sy = D.y - C.y;
                            T sLen = std::hypot(sx, sy);
                            if (sLen <= T(1e-12) || rLen <= T(1e-12)) continue;
                            T den = rx * sy - ry * sx;
                            // Reject grazing co-linear overlaps; stop on transverse hits (>= ~20 deg)
                            if (std::abs(den) < T(0.32) * rLen * sLen) continue;
                            T t = ((C.x - cur.x) * sy - (C.y - cur.y) * sx) / den;
                            T u = ((C.x - cur.x) * ry - (C.y - cur.y) * rx) / den;
                            if (t > T(1e-4) && t <= bestTrT && u >= -T(1e-3) && u <= T(1) + T(1e-3)) {
                                bestTrT = t;
                                bestTrHit = gsPoint2D<T>(cur.x + t * rx, cur.y + t * ry);
                                hitTrail = true;
                            }
                        }
                    }
                }

                if (hitTrail) {
                    r.pts.push_back(bestTrHit);
                    r.active = false;
                    r.terminated = true;
                } else if (bestBndLoop >= 0) {
                    T hitS = T(0);
                    gsPoint2D<T> projPt = bestBndHit;
                    if (bestBndLoop < nLoops)
                        projectPointToVirtualLoop(vloops[bestBndLoop], bestBndHit, hitS, projPt);
                    r.pts.push_back(projPt);
                    r.endLoop = bestBndLoop;
                    r.endArcS = hitS;
                    r.active = false;
                    r.terminated = true;
                } else {
                    if (!isPointStrictlyInsideDomain(nxt, fineMesh.loopPolys, T(1e-6))) {
                        // Project onto closest loop
                        int clL = 0; T clS = T(0); gsPoint2D<T> clP = nxt; T clD = T(1e30);
                        for (int l = 0; l < nLoops; ++l) {
                            T sTmp = T(0); gsPoint2D<T> pTmp;
                            T d = projectPointToVirtualLoop(vloops[l], cur, sTmp, pTmp);
                            if (d < clD) { clD = d; clL = l; clS = sTmp; clP = pTmp; }
                        }
                        r.pts.push_back(clP);
                        r.endLoop = clL;
                        r.endArcS = clS;
                        r.active = false;
                        r.terminated = true;
                    } else {
                        r.pts.push_back(nxt);
                        r.arcLen += stepLen;
                    }
                }
            }
        }
    };

    propagateActiveRays();

    // Remove unterminated rays (didn't hit boundary, trail, or source)
    for (auto& r : rays)
        if (!r.terminated || r.pts.size() < 2) r.alive = false;

    // 3. Topological Hole-Bridge & Connectedness Guarantee
    auto getLoopHitArcs = [&](int l) -> std::vector<T> {
        std::vector<T> hs;
        for (const auto& r : rays) {
            if (!r.alive) continue;
            if (r.startLoop == l) hs.push_back(r.startArcS);
            if (r.endLoop == l)   hs.push_back(r.endArcS);
        }
        return hs;
    };

    auto arcDistLoop = [&](int l, T a, T b) -> T {
        T d = std::abs(a - b);
        return std::min(d, vloops[l].totalArcLen - d);
    };

    // Launch a bridge ray from vloops[l] at arcS along cross-field/inward normal (stopping on trails/boundaries)
    auto launchHoleBridgeFromS = [&](int l, T sLaunch) -> bool {
        const auto& vl = vloops[l];
        gsPoint2D<T> p0 = evalVirtualLoopPoint(vl, sLaunch);
        gsPoint2D<T> t0 = evalVirtualLoopUnitTangent(vl, sLaunch);
        gsPoint2D<T> inN(-t0.y, t0.x); // left normal points strictly INTO domain!

        Ray r;
        r.id = static_cast<int>(rays.size());
        r.sourceId = nextCornerSrcId--;
        r.startLoop = l;
        r.startArcS = sLaunch;
        r.pts.push_back(p0);
        r.dir = inN;
        int th = -1;
        gsPoint2D<T> cfDir;
        gsPoint2D<T> probe(p0.x + T(0.4) * avgEdge * inN.x, p0.y + T(0.4) * avgEdge * inN.y);
        if (pickCrossBranch(probe, inN, th, cfDir) && (cfDir.x * inN.x + cfDir.y * inN.y > T(0.55)))
            r.dir = cfDir;
        rays.push_back(std::move(r));
        propagateActiveRays();
        auto& lr = rays.back();
        // Reject self-loop back to same hole if very short
        if (!lr.terminated || lr.pts.size() < 2 ||
            (lr.endLoop == l && arcDistLoop(l, lr.startArcS, lr.endArcS) < T(0.12) * vl.totalArcLen)) {
            lr.alive = false;
            return false;
        }
        return true;
    };

    // Rule H1: Every hole loop l >= 1 needs >= 2 well-spaced incident rays
    for (int l = 1; l < nLoops; ++l) {
        const auto& vl = vloops[l];
        if (vl.totalArcLen <= T(1e-10)) continue;
        for (int attempt = 0; attempt < 4; ++attempt) {
            std::vector<T> hits = getLoopHitArcs(l);
            // Count distinct well-spaced hits (>= 0.14 * L apart)
            std::vector<T> uniq;
            for (T s : hits) {
                bool dup = false;
                for (T u : uniq) if (arcDistLoop(l, s, u) < T(0.14) * vl.totalArcLen) { dup = true; break; }
                if (!dup) uniq.push_back(s);
            }
            if (uniq.size() >= 2) break;

            // Choose best launch point on vloops[l] maximizing spacing from existing hits + cross-field alignment
            T bestS = T(0);
            T bestScore = -T(1e9);
            const int nCand = 24;
            for (int c = 0; c < nCand; ++c) {
                T sc = (static_cast<T>(c) / static_cast<T>(nCand)) * vl.totalArcLen;
                T minSp = vl.totalArcLen;
                for (T u : uniq) minSp = std::min(minSp, arcDistLoop(l, sc, u));
                if (!uniq.empty() && minSp < T(0.22) * vl.totalArcLen) continue;
                gsPoint2D<T> p = evalVirtualLoopPoint(vl, sc);
                gsPoint2D<T> t = evalVirtualLoopUnitTangent(vl, sc);
                gsPoint2D<T> n(-t.y, t.x);
                gsPoint2D<T> pr(p.x + T(0.5) * avgEdge * n.x, p.y + T(0.5) * avgEdge * n.y);
                if (!isPointStrictlyInsideDomain(pr, fineMesh.loopPolys, T(1e-6))) continue;
                std::complex<T> zq; T th = T(0); int triH = -1;
                T align = T(0);
                if (evalCrossFieldAt(fineMesh, cf, accel, pr, zq, th, triH) && std::abs(zq) > T(1e-10)) {
                    std::complex<T> n4 = std::pow(std::complex<T>(n.x, n.y), 4);
                    align = (zq.real() * n4.real() + zq.imag() * n4.imag()) / std::abs(zq);
                }
                T val = T(2.0) * align + (uniq.empty() ? T(0) : (minSp / vl.totalArcLen));
                if (val > bestScore) { bestScore = val; bestS = sc; }
            }
            if (bestScore > -T(1e8)) {
                if (!launchHoleBridgeFromS(l, bestS)) {
                    // Fallback: shortest visible straight bridge to any other loop
                    gsPoint2D<T> pA = evalVirtualLoopPoint(vl, bestS);
                    T bestD = T(1e30);
                    gsPoint2D<T> bestB;
                    int bestOL = -1; T bestOS = T(0);
                    for (int ol = 0; ol < nLoops; ++ol) {
                        if (ol == l) continue;
                        for (int k = 0; k < 24; ++k) {
                            T so = (static_cast<T>(k) / T(24)) * vloops[ol].totalArcLen;
                            gsPoint2D<T> pB = evalVirtualLoopPoint(vloops[ol], so);
                            T d = pA.dist(pB);
                            if (d < bestD && isSegmentInsideDomain(pA, pB, fineMesh.loopPolys)) {
                                bestD = d; bestB = pB; bestOL = ol; bestOS = so;
                            }
                        }
                    }
                    if (bestOL >= 0) {
                        Ray rb;
                        rb.id = static_cast<int>(rays.size());
                        rb.sourceId = nextCornerSrcId--;
                        rb.startLoop = l; rb.startArcS = bestS;
                        rb.endLoop = bestOL; rb.endArcS = bestOS;
                        rb.pts = {pA, bestB};
                        rb.active = false; rb.alive = true; rb.terminated = true;
                        rays.push_back(std::move(rb));
                    }
                }
            }
        }
    }

    // Rule H2: Guarantee ALL loops are in ONE single connected component with the outer loop (0)!
    for (int connIter = 0; connIter < nLoops; ++connIter) {
        // Build Union-Find over loops [0..nLoops-1] + live rays [nLoops .. nLoops+nRays-1]
        const int nR = static_cast<int>(rays.size());
        std::vector<int> uf(nLoops + nR);
        std::iota(uf.begin(), uf.end(), 0);
        auto ufFind = [&](int x) {
            while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
            return x;
        };
        auto ufUnite = [&](int a, int b) {
            a = ufFind(a); b = ufFind(b);
            if (a != b) uf[a] = b;
        };
        for (int i = 0; i < nR; ++i) {
            if (!rays[i].alive) continue;
            if (rays[i].startLoop >= 0) ufUnite(rays[i].startLoop, nLoops + i);
            if (rays[i].endLoop >= 0)   ufUnite(rays[i].endLoop, nLoops + i);
        }
        // Connect rays that share a source or intersect/touch each other
        for (int i = 0; i < nR; ++i) {
            if (!rays[i].alive) continue;
            for (int j = i + 1; j < nR; ++j) {
                if (!rays[j].alive) continue;
                if (rays[i].sourceId == rays[j].sourceId ||
                    rays[i].pts.front().dist(rays[j].pts.front()) < T(0.8) * avgEdge ||
                    rays[i].pts.back().dist(rays[j].pts.back()) < T(0.8) * avgEdge ||
                    rays[i].pts.back().dist(rays[j].pts.front()) < T(0.8) * avgEdge ||
                    rays[i].pts.front().dist(rays[j].pts.back()) < T(0.8) * avgEdge) {
                    ufUnite(nLoops + i, nLoops + j);
                    continue;
                }
                // Check if ray endpoint touches other ray polyline (T-junction)
                auto touches = [&](const gsPoint2D<T>& q, const std::vector<gsPoint2D<T>>& pl) -> bool {
                    for (size_t k = 0; k + 1 < pl.size(); ++k) {
                        T dx = pl[k + 1].x - pl[k].x, dy = pl[k + 1].y - pl[k].y, L2 = dx * dx + dy * dy;
                        T t = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((q.x - pl[k].x) * dx + (q.y - pl[k].y) * dy) / L2)) : T(0);
                        if (std::hypot(q.x - (pl[k].x + t * dx), q.y - (pl[k].y + t * dy)) < T(0.35) * avgEdge)
                            return true;
                    }
                    return false;
                };
                if (touches(rays[i].pts.back(), rays[j].pts) || touches(rays[j].pts.back(), rays[i].pts))
                    ufUnite(nLoops + i, nLoops + j);
            }
        }

        int root0 = ufFind(0);
        int disLoop = -1;
        for (int l = 1; l < nLoops; ++l)
            if (ufFind(l) != root0) { disLoop = l; break; }
        if (disLoop < 0) break;

        // Shortest straight interior bridge between disLoop and any loop currently in root0
        T bestD = T(1e30);
        gsPoint2D<T> bestA, bestB;
        int bestOL = -1; T bestSA = T(0), bestSB = T(0);
        for (int ka = 0; ka < 24; ++ka) {
            T sa = (static_cast<T>(ka) / T(24)) * vloops[disLoop].totalArcLen;
            gsPoint2D<T> pA = evalVirtualLoopPoint(vloops[disLoop], sa);
            for (int ol = 0; ol < nLoops; ++ol) {
                if (ufFind(ol) != root0) continue;
                for (int kb = 0; kb < 28; ++kb) {
                    T sb = (static_cast<T>(kb) / T(28)) * vloops[ol].totalArcLen;
                    gsPoint2D<T> pB = evalVirtualLoopPoint(vloops[ol], sb);
                    T d = pA.dist(pB);
                    if (d < bestD && isSegmentInsideDomain(pA, pB, fineMesh.loopPolys)) {
                        bestD = d; bestA = pA; bestB = pB; bestOL = ol; bestSA = sa; bestSB = sb;
                    }
                }
            }
        }
        if (bestOL < 0) break;
        Ray rb;
        rb.id = static_cast<int>(rays.size());
        rb.sourceId = nextCornerSrcId--;
        rb.startLoop = disLoop; rb.startArcS = bestSA;
        rb.endLoop = bestOL; rb.endArcS = bestSB;
        rb.pts = {bestA, bestB};
        rb.active = false; rb.alive = true; rb.terminated = true;
        rays.push_back(std::move(rb));
    }

    // Collect live rays
    std::vector<Ray> liveRays;
    for (auto& r : rays)
        if (r.alive && r.pts.size() >= 2 && r.pts.front().dist(r.pts.back()) > T(1e-5) * diag)
            liveRays.push_back(std::move(r));
    if (liveRays.empty()) return false;

    // 4. Exact Planar Arrangement: Boundary VirtualCADLoops + Live Rays (with T-Junction & Crossing Splits)
    TMeshGraph<T> G;
    const T nodeMergeTol = T(0.22) * avgEdge;

    auto addNode = [&](const gsPoint2D<T>& pt, bool isBnd, int lIdx, T arcS, bool isPinned, T mergeTol) -> int {
        for (size_t i = 0; i < G.nodes.size(); ++i) {
            if (isBnd != G.nodes[i].isBoundary) continue;
            if (isBnd && G.nodes[i].loopIdx != lIdx) continue;
            if (G.nodes[i].pt.dist(pt) <= mergeTol) {
                if (isPinned) G.nodes[i].isPinnedCorner = true;
                return static_cast<int>(i);
            }
        }
        int id = static_cast<int>(G.nodes.size());
        typename TMeshGraph<T>::Node nd;
        nd.pt = pt; nd.isBoundary = isBnd; nd.loopIdx = lIdx; nd.loopArcS = arcS; nd.isPinnedCorner = isPinned;
        G.nodes.push_back(nd);
        return id;
    };

    // Collect all boundary split arc-lengths on each VirtualCADLoop
    std::vector<std::vector<std::pair<T, bool>>> loopSplitPts(nLoops); // (arcS, isPinned)
    for (int l = 0; l < nLoops; ++l) {
        for (const auto& pc : vloops[l].corners)
            loopSplitPts[l].push_back({pc.arcS, true});
    }
    for (const auto& r : liveRays) {
        if (r.startLoop >= 0 && r.startLoop < nLoops)
            loopSplitPts[r.startLoop].push_back({virtualLoopWrapS(vloops[r.startLoop], r.startArcS), false});
        if (r.endLoop >= 0 && r.endLoop < nLoops)
            loopSplitPts[r.endLoop].push_back({virtualLoopWrapS(vloops[r.endLoop], r.endArcS), false});
    }

    // Deduplicate boundary split points (keeping pinned corners fixed!) and ensure >= 3 vertices per loop
    std::vector<std::vector<std::pair<T, bool>>> cleanLoopPts(nLoops);
    std::vector<std::vector<int>> cleanLoopNodeIds(nLoops);
    for (int l = 0; l < nLoops; ++l) {
        const auto& vl = vloops[l];
        const T L = vl.totalArcLen;
        if (L <= T(1e-12)) continue;
        const T minGap = std::min(T(0.03) * diag, T(0.08) * L);

        // First keep all pinned corners
        std::vector<std::pair<T, bool>> kept;
        for (const auto& sp : loopSplitPts[l]) {
            if (!sp.second) continue;
            bool dup = false;
            for (const auto& k : kept)
                if (arcDistLoop(l, sp.first, k.first) < T(1e-4) * L) { dup = true; break; }
            if (!dup) kept.push_back(sp);
        }
        // Then add non-pinned ray hits not too close to existing points
        for (const auto& sp : loopSplitPts[l]) {
            if (sp.second) continue;
            bool dup = false;
            for (const auto& k : kept)
                if (arcDistLoop(l, sp.first, k.first) < minGap) { dup = true; break; }
            if (!dup) kept.push_back(sp);
        }
        // Ensure smooth loop has at least 4 nodes (or bisect > 135 deg arcs on 2-node smooth loops)
        if (kept.empty()) {
            for (int k = 0; k < 4; ++k) kept.push_back({static_cast<T>(k) * T(0.25) * L, false});
        } else if (kept.size() == 1) {
            T s0 = kept[0].first;
            for (int k = 1; k < 4; ++k) kept.push_back({virtualLoopWrapS(vl, s0 + static_cast<T>(k) * T(0.25) * L), false});
        } else if (kept.size() == 2) {
            std::sort(kept.begin(), kept.end());
            T s0 = kept[0].first, s1 = kept[1].first;
            kept.push_back({virtualLoopWrapS(vl, T(0.5) * (s0 + s1)), false});
            kept.push_back({virtualLoopWrapS(vl, s1 + T(0.5) * (L - (s1 - s0))), false});
        }
        std::sort(kept.begin(), kept.end(), [](const std::pair<T, bool>& a, const std::pair<T, bool>& b) { return a.first < b.first; });
        cleanLoopPts[l] = kept;

        for (const auto& kp : kept) {
            gsPoint2D<T> pt = evalVirtualLoopPoint(vl, kp.first);
            int nid = addNode(pt, true, l, kp.first, kp.second, T(1e-5) * diag);
            cleanLoopNodeIds[l].push_back(nid);
        }

        const size_t nK = kept.size();
        for (size_t k = 0; k < nK; ++k) {
            size_t kn = (k + 1) % nK;
            T sA = kept[k].first, sB = kept[kn].first;
            CurveSeg<T> cSeg = extractVirtualSubCurve(vl, sA, sB, degree);
            if (!cSeg.bsp) continue;
            typename TMeshGraph<T>::Segment sg;
            sg.u = cleanLoopNodeIds[l][k];
            sg.v = cleanLoopNodeIds[l][kn];
            sg.seg = cSeg;
            sg.isBoundary = true;
            sg.loopIdx = l;
            sg.s0 = sA;
            sg.s1 = sB;
            sg.arcLen = (sB > sA) ? (sB - sA) : (sB + L - sA);
            if (sg.u != sg.v) G.segs.push_back(std::move(sg));
        }
    }

    // Snap every liveRay boundary endpoint to its exact kept boundary node
    auto snapRayEndToLoop = [&](int l, T s, gsPoint2D<T>& pt) -> int {
        if (l < 0 || l >= nLoops || cleanLoopPts[l].empty()) return -1;
        int bestK = 0;
        T bestD = T(1e30);
        for (size_t k = 0; k < cleanLoopPts[l].size(); ++k) {
            T d = arcDistLoop(l, s, cleanLoopPts[l][k].first);
            if (d < bestD) { bestD = d; bestK = static_cast<int>(k); }
        }
        int nid = cleanLoopNodeIds[l][bestK];
        pt = G.nodes[nid].pt;
        return nid;
    };

    const size_t nLR = liveRays.size();
    std::vector<int> rayStartNode(nLR, -1), rayEndNode(nLR, -1);
    for (size_t i = 0; i < nLR; ++i) {
        if (liveRays[i].startLoop >= 0)
            rayStartNode[i] = snapRayEndToLoop(liveRays[i].startLoop, liveRays[i].startArcS, liveRays[i].pts.front());
        if (liveRays[i].endLoop >= 0)
            rayEndNode[i] = snapRayEndToLoop(liveRays[i].endLoop, liveRays[i].endArcS, liveRays[i].pts.back());
    }

    // Compute cumulative arc length & ALL pairwise trail T-junctions / crossings for exact planarity
    std::vector<std::vector<T>> polyCumS(nLR);
    for (size_t i = 0; i < nLR; ++i) {
        polyCumS[i].resize(liveRays[i].pts.size(), T(0));
        for (size_t k = 1; k < liveRays[i].pts.size(); ++k)
            polyCumS[i][k] = polyCumS[i][k - 1] + liveRays[i].pts[k - 1].dist(liveRays[i].pts[k]);
    }

    std::vector<std::vector<T>> raySplits(nLR);
    for (size_t i = 0; i < nLR; ++i) {
        const auto& Pi = liveRays[i].pts;
        const T Li = polyCumS[i].back();
        if (Li <= T(1e-12)) continue;

        // Check if endpoint of ray i lands on the interior of ray j (T-junction!)
        for (int endFlag = 0; endFlag < 2; ++endFlag) {
            if (endFlag == 0 && rayStartNode[i] >= 0) continue;
            if (endFlag == 1 && rayEndNode[i] >= 0) continue;
            const gsPoint2D<T>& q = (endFlag == 0) ? Pi.front() : Pi.back();
            for (size_t j = 0; j < nLR; ++j) {
                if (j == i) continue;
                const auto& Pj = liveRays[j].pts;
                const T Lj = polyCumS[j].back();
                if (Lj <= T(1e-12)) continue;
                T bestD = T(0.45) * avgEdge;
                T bestSj = -T(1);
                gsPoint2D<T> bestPt = q;
                for (size_t kj = 0; kj + 1 < Pj.size(); ++kj) {
                    T dx = Pj[kj + 1].x - Pj[kj].x, dy = Pj[kj + 1].y - Pj[kj].y, L2 = dx * dx + dy * dy;
                    if (L2 <= T(1e-20)) continue;
                    T u = std::max(T(0), std::min(T(1), ((q.x - Pj[kj].x) * dx + (q.y - Pj[kj].y) * dy) / L2));
                    gsPoint2D<T> proj(Pj[kj].x + u * dx, Pj[kj].y + u * dy);
                    T d = q.dist(proj);
                    if (d < bestD) {
                        bestD = d;
                        bestPt = proj;
                        bestSj = (polyCumS[j][kj] + u * std::sqrt(L2)) / Lj;
                    }
                }
                if (bestSj >= T(0)) {
                    if (bestSj > T(0.03) && bestSj < T(0.97)) {
                        raySplits[j].push_back(bestSj);
                        if (endFlag == 0) liveRays[i].pts.front() = bestPt;
                        else              liveRays[i].pts.back()  = bestPt;
                    } else if (bestSj <= T(0.03)) {
                        if (endFlag == 0) { liveRays[i].pts.front() = Pj.front(); if (rayStartNode[j] >= 0) rayStartNode[i] = rayStartNode[j]; }
                        else              { liveRays[i].pts.back()  = Pj.front(); if (rayStartNode[j] >= 0) rayEndNode[i]   = rayStartNode[j]; }
                    } else {
                        if (endFlag == 0) { liveRays[i].pts.front() = Pj.back(); if (rayEndNode[j] >= 0) rayStartNode[i] = rayEndNode[j]; }
                        else              { liveRays[i].pts.back()  = Pj.back(); if (rayEndNode[j] >= 0) rayEndNode[i]   = rayEndNode[j]; }
                    }
                }
            }
        }

        // Safety net: also split any remaining interior crossings (e.g. from hole-bridge segments)
        for (size_t j = i + 1; j < nLR; ++j) {
            const auto& Pj = liveRays[j].pts;
            const T Lj = polyCumS[j].back();
            if (Lj <= T(1e-12)) continue;
            for (size_t ki = 0; ki + 1 < Pi.size(); ++ki) {
                const auto& A = Pi[ki]; const auto& B = Pi[ki + 1];
                T rx = B.x - A.x, ry = B.y - A.y, lI = std::hypot(rx, ry);
                if (lI <= T(1e-12)) continue;
                for (size_t kj = 0; kj + 1 < Pj.size(); ++kj) {
                    const auto& C = Pj[kj]; const auto& D = Pj[kj + 1];
                    T sx = D.x - C.x, sy = D.y - C.y, lJ = std::hypot(sx, sy);
                    if (lJ <= T(1e-12)) continue;
                    T den = rx * sy - ry * sx;
                    if (std::abs(den) < T(0.18) * lI * lJ) continue;
                    T t = ((C.x - A.x) * sy - (C.y - A.y) * sx) / den;
                    T u = ((C.x - A.x) * ry - (C.y - A.y) * rx) / den;
                    if (t >= T(0) && t <= T(1) && u >= T(0) && u <= T(1)) {
                        T si = (polyCumS[i][ki] + t * lI) / Li;
                        T sj = (polyCumS[j][kj] + u * lJ) / Lj;
                        if (si > T(0.03) && si < T(0.97) && sj > T(0.03) && sj < T(0.97)) {
                            raySplits[i].push_back(si);
                            raySplits[j].push_back(sj);
                        }
                    }
                }
            }
        }
    }

    // Recompute polyCumS after any endpoint snaps (T-junction or boundary) before extracting sub-polylines!
    for (size_t i = 0; i < nLR; ++i) {
        polyCumS[i].assign(liveRays[i].pts.size(), T(0));
        for (size_t k = 1; k < liveRays[i].pts.size(); ++k)
            polyCumS[i][k] = polyCumS[i][k - 1] + liveRays[i].pts[k - 1].dist(liveRays[i].pts[k]);
    }

    auto extractSubPoly = [&](const std::vector<gsPoint2D<T>>& pts, const std::vector<T>& cumS,
                              T f0, T f1) -> std::vector<gsPoint2D<T>> {
        const T tot = cumS.back();
        const T s0 = f0 * tot, s1 = f1 * tot;
        auto sampleS = [&](T st) -> gsPoint2D<T> {
            for (size_t k = 1; k < pts.size(); ++k) {
                if (cumS[k] >= st) {
                    T ds = cumS[k] - cumS[k - 1];
                    T al = (ds > T(1e-14)) ? (st - cumS[k - 1]) / ds : T(0);
                    return gsPoint2D<T>((T(1) - al) * pts[k - 1].x + al * pts[k].x,
                                        (T(1) - al) * pts[k - 1].y + al * pts[k].y);
                }
            }
            return pts.back();
        };
        std::vector<gsPoint2D<T>> sub;
        sub.push_back(sampleS(s0));
        for (size_t k = 1; k + 1 < pts.size(); ++k)
            if (cumS[k] > s0 + T(1e-9) && cumS[k] < s1 - T(1e-9))
                sub.push_back(pts[k]);
        sub.push_back(sampleS(s1));
        return sub;
    };

    std::set<std::pair<int, int>> existingPairs;
    for (const auto& s : G.segs)
        existingPairs.insert({std::min(s.u, s.v), std::max(s.u, s.v)});

    for (size_t i = 0; i < nLR; ++i) {
        auto& ts = raySplits[i];
        ts.push_back(T(0)); ts.push_back(T(1));
        std::sort(ts.begin(), ts.end());
        std::vector<T> uniq;
        for (T t : ts)
            if (uniq.empty() || t - uniq.back() > T(0.02)) uniq.push_back(t);
        if (uniq.size() < 2) uniq = {T(0), T(1)};
        uniq.front() = T(0); uniq.back() = T(1);

        for (size_t k = 0; k + 1 < uniq.size(); ++k) {
            std::vector<gsPoint2D<T>> sub = extractSubPoly(liveRays[i].pts, polyCumS[i], uniq[k], uniq[k + 1]);
            int va = (k == 0 && rayStartNode[i] >= 0) ? rayStartNode[i]
                   : addNode(sub.front(), false, -1, T(0), false, nodeMergeTol);
            int vb = (k + 2 == uniq.size() && rayEndNode[i] >= 0) ? rayEndNode[i]
                   : addNode(sub.back(), false, -1, T(0), false, nodeMergeTol);
            if (va == vb) continue;
            auto key = std::make_pair(std::min(va, vb), std::max(va, vb));
            if (existingPairs.count(key)) continue;
            existingPairs.insert(key);
            sub.front() = G.nodes[va].pt;
            sub.back()  = G.nodes[vb].pt;
            CurveSeg<T> fitted = fitStreamlinePolylineToCurveSeg(sub, degree, fineMesh.loopPolys);
            typename TMeshGraph<T>::Segment sg;
            sg.u = va; sg.v = vb; sg.seg = fitted; sg.isBoundary = false; sg.loopIdx = -1;
            sg.arcLen = measureCurveSegArcLen(fitted);
            G.segs.push_back(std::move(sg));
        }
    }

    auto pruneInteriorSpurs = [&]() {
        bool pruned = true;
        while (pruned) {
            pruned = false;
            std::vector<int> deg(G.nodes.size(), 0);
            for (const auto& s : G.segs) { ++deg[s.u]; ++deg[s.v]; }
            for (size_t i = 0; i < G.segs.size(); ++i) {
                if (G.segs[i].isBoundary) continue;
                if (deg[G.segs[i].u] <= 1 || deg[G.segs[i].v] <= 1) {
                    G.segs.erase(G.segs.begin() + static_cast<long>(i));
                    pruned = true;
                    break;
                }
            }
        }
    };
    auto segExists = [&](int u, int v) -> bool {
        if (u == v) return true;
        for (const auto& s : G.segs)
            if ((s.u == u && s.v == v) || (s.u == v && s.v == u)) return true;
        return false;
    };
    auto splitGraphSegAtFrac = [&](int segIdx, T frac) -> int {
        frac = std::max(T(0.05), std::min(T(0.95), frac));
        const auto orig = G.segs[segIdx];
        CurveSeg<T> leftSeg, rightSeg;
        T midS = T(0);
        T leftArc = T(0), rightArc = T(0);
        if (orig.isBoundary && orig.loopIdx >= 0 && orig.loopIdx < nLoops) {
            const auto& vl = vloops[orig.loopIdx];
            leftArc  = frac * orig.arcLen;
            rightArc = orig.arcLen - leftArc;
            midS = virtualLoopWrapS(vl, orig.s0 + leftArc);
            leftSeg  = extractVirtualSubCurve(vl, orig.s0, midS, degree);
            rightSeg = extractVirtualSubCurve(vl, midS, orig.s1, degree);
        } else {
            auto sp = orig.seg.splitArc(frac);
            leftSeg = sp.first;
            rightSeg = sp.second;
            leftArc = measureCurveSegArcLen(leftSeg);
            rightArc = measureCurveSegArcLen(rightSeg);
        }
        if (!leftSeg.bsp || !rightSeg.bsp) return -1;
        gsPoint2D<T> mPt = leftSeg.p1();
        typename TMeshGraph<T>::Node nd;
        nd.pt = mPt;
        nd.isBoundary = orig.isBoundary;
        nd.loopIdx = orig.loopIdx;
        nd.loopArcS = midS;
        nd.isTJunction = true;
        int mId = static_cast<int>(G.nodes.size());
        G.nodes.push_back(nd);

        G.segs[segIdx].v = mId;
        G.segs[segIdx].seg = leftSeg;
        G.segs[segIdx].s1 = midS;
        G.segs[segIdx].arcLen = leftArc;

        typename TMeshGraph<T>::Segment rsg;
        rsg.u = mId; rsg.v = orig.v;
        rsg.seg = rightSeg;
        rsg.isBoundary = orig.isBoundary;
        rsg.loopIdx = orig.loopIdx;
        rsg.s0 = midS; rsg.s1 = orig.s1;
        rsg.arcLen = rightArc;
        G.segs.push_back(std::move(rsg));
        return mId;
    };

    pruneInteriorSpurs();

    // 4B. Graph-Level Planar Hole-Bridge & Connectedness Guarantee DIRECTLY on G!
    //     Guarantees EVERY hole loop l >= 1 has >= 2 incident interior segments AND is connected to outer loop 0,
    //     splitting any crossed existing segment at the exact earliest hit so planarity is 100% preserved!
    auto addPlanarGraphBridge = [&](int uNode, int vNode) -> bool {
        if (uNode < 0 || vNode < 0 || uNode == vNode || segExists(uNode, vNode)) return false;
        const gsPoint2D<T> A = G.nodes[uNode].pt;
        const gsPoint2D<T> B = G.nodes[vNode].pt;
        if (A.dist(B) <= T(1e-5) * diag) return false;
        if (!isSegmentInsideDomain(A, B, fineMesh.loopPolys)) return false;

        // Check earliest crossing against ALL existing segments of G (12 polyline samples per segment)
        T bestT = T(1) + T(1e-4);
        int hitSeg = -1;
        T hitFrac = T(0.5);
        T rx = B.x - A.x, ry = B.y - A.y;
        const int nS = 12;
        for (size_t ei = 0; ei < G.segs.size(); ++ei) {
            const auto& sg = G.segs[ei];
            if (sg.u == uNode || sg.v == uNode || sg.u == vNode || sg.v == vNode) continue;
            gsPoint2D<T> prev = evalCurvePoint(*sg.seg.bsp, T(0));
            for (int k = 1; k <= nS; ++k) {
                gsPoint2D<T> cur = evalCurvePoint(*sg.seg.bsp, static_cast<T>(k) / static_cast<T>(nS));
                T sx = cur.x - prev.x, sy = cur.y - prev.y;
                T den = rx * sy - ry * sx;
                if (std::abs(den) > T(1e-14)) {
                    T t = ((prev.x - A.x) * sy - (prev.y - A.y) * sx) / den;
                    T u = ((prev.x - A.x) * ry - (prev.y - A.y) * rx) / den;
                    if (t > T(1e-4) && t < bestT && u >= -T(1e-3) && u <= T(1) + T(1e-3)) {
                        bestT = t;
                        hitSeg = static_cast<int>(ei);
                        hitFrac = (static_cast<T>(k - 1) + std::max(T(0), std::min(T(1), u))) / static_cast<T>(nS);
                    }
                }
                prev = cur;
            }
        }

        int finalV = vNode;
        if (hitSeg >= 0) {
            if (hitFrac < T(0.10))      finalV = G.segs[hitSeg].u;
            else if (hitFrac > T(0.90)) finalV = G.segs[hitSeg].v;
            else                        finalV = splitGraphSegAtFrac(hitSeg, hitFrac);
        }
        if (finalV < 0 || finalV == uNode || segExists(uNode, finalV)) return false;
        if (!isSegmentInsideDomain(A, G.nodes[finalV].pt, fineMesh.loopPolys)) return false;

        CurveSeg<T> cs = makeStraightSeg(A, G.nodes[finalV].pt, degree);
        typename TMeshGraph<T>::Segment ns;
        ns.u = uNode; ns.v = finalV; ns.seg = cs; ns.isBoundary = false; ns.loopIdx = -1;
        ns.arcLen = A.dist(G.nodes[finalV].pt);
        G.segs.push_back(std::move(ns));
        return true;
    };

    for (int gBridgeIter = 0; gBridgeIter < 4 * nLoops + 4; ++gBridgeIter) {
        const int nN = static_cast<int>(G.nodes.size());
        std::vector<int> uf(nN);
        std::iota(uf.begin(), uf.end(), 0);
        auto ufF = [&](int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
        auto ufU = [&](int a, int b) { a = ufF(a); b = ufF(b); if (a != b) uf[a] = b; };

        std::vector<int> loopIntDeg(nLoops, 0);
        std::vector<int> nodeIntDeg(nN, 0);
        std::vector<int> loopRepNode(nLoops, -1);
        for (int ni = 0; ni < nN; ++ni)
            if (G.nodes[ni].isBoundary && G.nodes[ni].loopIdx >= 0 && G.nodes[ni].loopIdx < nLoops && loopRepNode[G.nodes[ni].loopIdx] < 0)
                loopRepNode[G.nodes[ni].loopIdx] = ni;
        for (const auto& s : G.segs) {
            ufU(s.u, s.v);
            if (!s.isBoundary) {
                ++nodeIntDeg[s.u]; ++nodeIntDeg[s.v];
                if (G.nodes[s.u].isBoundary && G.nodes[s.u].loopIdx >= 0) ++loopIntDeg[G.nodes[s.u].loopIdx];
                if (G.nodes[s.v].isBoundary && G.nodes[s.v].loopIdx >= 0) ++loopIntDeg[G.nodes[s.v].loopIdx];
            }
        }
        int root0 = (loopRepNode[0] >= 0) ? ufF(loopRepNode[0]) : -1;
        int needLoop = -1;
        bool needConnectOuter = false;
        for (int l = 1; l < nLoops; ++l) {
            if (loopRepNode[l] < 0) continue;
            if (root0 >= 0 && ufF(loopRepNode[l]) != root0) { needLoop = l; needConnectOuter = true; break; }
        }
        if (needLoop < 0) {
            for (int l = 1; l < nLoops; ++l) {
                if (loopRepNode[l] >= 0 && loopIntDeg[l] < 2) { needLoop = l; break; }
            }
        }
        if (needLoop < 0) break;

        // Collect candidate source nodes on needLoop (prefer unused boundary nodes spaced from existing bridges!)
        std::vector<int> srcNodes;
        for (int ni = 0; ni < nN; ++ni)
            if (G.nodes[ni].isBoundary && G.nodes[ni].loopIdx == needLoop && nodeIntDeg[ni] == 0)
                srcNodes.push_back(ni);
        if (srcNodes.empty()) {
            for (int ni = 0; ni < nN; ++ni)
                if (G.nodes[ni].isBoundary && G.nodes[ni].loopIdx == needLoop)
                    srcNodes.push_back(ni);
        }

        int bestU = -1, bestV = -1;
        T bestD = T(1e30);
        for (int u : srcNodes) {
            for (int v = 0; v < nN; ++v) {
                if (v == u) continue;
                if (G.nodes[v].isBoundary && G.nodes[v].loopIdx == needLoop) continue;
                if (needConnectOuter && root0 >= 0 && ufF(v) != root0) continue;
                if (segExists(u, v)) continue;
                T d = G.nodes[u].pt.dist(G.nodes[v].pt);
                if (d < bestD && isSegmentInsideDomain(G.nodes[u].pt, G.nodes[v].pt, fineMesh.loopPolys)) {
                    bestD = d; bestU = u; bestV = v;
                }
            }
        }
        if (bestU < 0 && needConnectOuter) {
            // Relax: bridge to ANY visible node not in needLoop's component (multi-hop bridges toward outer!)
            int rootL = ufF(loopRepNode[needLoop]);
            for (int u : srcNodes) {
                for (int v = 0; v < nN; ++v) {
                    if (v == u || ufF(v) == rootL || segExists(u, v)) continue;
                    T d = G.nodes[u].pt.dist(G.nodes[v].pt);
                    if (d < bestD && isSegmentInsideDomain(G.nodes[u].pt, G.nodes[v].pt, fineMesh.loopPolys)) {
                        bestD = d; bestU = u; bestV = v;
                    }
                }
            }
        }
        if (bestU < 0 || !addPlanarGraphBridge(bestU, bestV)) break;
    }

    if (G.segs.size() < 4) return false;

    // 5. Planar Half-Edge Face Extraction & Automatic RectFace Classification/Rectification
    auto buildHalfEdgesAndExtractRawFaces = [&](std::vector<std::vector<int>>& rawFaces) -> bool {
        pruneInteriorSpurs();
        const int nE = static_cast<int>(G.segs.size());
        G.halfEdges.assign(2 * nE, typename TMeshGraph<T>::HalfEdge());
        std::vector<std::vector<int>> outHalf(G.nodes.size());

        for (int i = 0; i < nE; ++i) {
            const auto& e = G.segs[i];
            int h0 = 2 * i, h1 = 2 * i + 1;
            gsPoint2D<T> d0 = evalCurveDeriv(*e.seg.bsp, T(0));
            gsPoint2D<T> d1 = evalCurveDeriv(*e.seg.bsp, T(1));
            T a0 = (std::hypot(d0.x, d0.y) > T(1e-12))
                 ? std::atan2(d0.y, d0.x)
                 : std::atan2(G.nodes[e.v].pt.y - G.nodes[e.u].pt.y, G.nodes[e.v].pt.x - G.nodes[e.u].pt.x);
            T a1 = (std::hypot(d1.x, d1.y) > T(1e-12))
                 ? std::atan2(-d1.y, -d1.x)
                 : std::atan2(G.nodes[e.u].pt.y - G.nodes[e.v].pt.y, G.nodes[e.u].pt.x - G.nodes[e.v].pt.x);
            G.halfEdges[h0] = {e.u, e.v, i, true,  h1, -1, a0};
            G.halfEdges[h1] = {e.v, e.u, i, false, h0, -1, a1};
            outHalf[e.u].push_back(h0);
            outHalf[e.v].push_back(h1);
        }

        for (size_t v = 0; v < G.nodes.size(); ++v) {
            auto& oh = outHalf[v];
            std::sort(oh.begin(), oh.end(), [&](int a, int b) { return G.halfEdges[a].angle < G.halfEdges[b].angle; });
        }

        for (int h = 0; h < 2 * nE; ++h) {
            int v = G.halfEdges[h].target;
            int tw = G.halfEdges[h].twin;
            const auto& oh = outHalf[v];
            auto it = std::find(oh.begin(), oh.end(), tw);
            if (it == oh.end() || oh.empty()) return false;
            size_t pos = static_cast<size_t>(it - oh.begin());
            G.halfEdges[h].next = oh[(pos + oh.size() - 1) % oh.size()];
        }

        std::vector<bool> vis(2 * nE, false);
        rawFaces.clear();
        for (int h = 0; h < 2 * nE; ++h) {
            if (vis[h]) continue;
            std::vector<int> cyc;
            int cur = h;
            bool okCyc = true;
            for (int guard = 0; guard < 2 * nE + 2; ++guard) {
                if (vis[cur]) { okCyc = (cur == h); break; }
                vis[cur] = true;
                cyc.push_back(cur);
                cur = G.halfEdges[cur].next;
                if (cur < 0) { okCyc = false; break; }
            }
            if (!okCyc || cyc.size() < 2) continue;

            // Exact topological interior-face criterion:
            // Every VirtualCADLoop is oriented with domain interior strictly to the LEFT (outer CCW, holes CW).
            // Therefore, exterior and hole-interior faces ALWAYS contain backward boundary half-edges (!he.forward),
            // whereas EVERY genuine interior domain face contains ZERO backward boundary half-edges!
            bool hasBackBnd = false;
            for (int ch : cyc) {
                const auto& he = G.halfEdges[ch];
                const auto& sg = G.segs[he.segIdx];
                if (sg.isBoundary && !he.forward) { hasBackBnd = true; break; }
            }
            if (hasBackBnd) continue;

            std::vector<gsPoint2D<T>> poly;
            poly.reserve(cyc.size() * 6);
            for (int ch : cyc) {
                const auto& he = G.halfEdges[ch];
                const auto& sg = G.segs[he.segIdx];
                for (int s = 0; s < 6; ++s) {
                    T u = he.forward ? (static_cast<T>(s) / T(6)) : (T(1) - static_cast<T>(s) / T(6));
                    poly.push_back(evalCurvePoint(*sg.seg.bsp, u));
                }
            }
            T area2 = T(0);
            for (size_t k = 0; k < poly.size(); ++k) {
                const auto& a = poly[k];
                const auto& b = poly[(k + 1) % poly.size()];
                area2 += a.x * b.y - b.x * a.y;
            }
            if (area2 <= T(1e-7)) continue;
            rawFaces.push_back(std::move(cyc));
        }
        return !rawFaces.empty();
    };

    auto getHECurve = [&](int heIdx) -> CurveSeg<T> {
        const auto& he = G.halfEdges[heIdx];
        const auto& sg = G.segs[he.segIdx];
        return he.forward ? sg.seg : sg.seg.reversed();
    };

    // Helper: exact straight ray-cast from vertex `startK` along `dir` against face boundary `cyc`
    // (by Jordan Curve Theorem, guaranteed to hit a non-adjacent edge for any strictly inward ray!)
    auto castStraightFaceChord = [&](const std::vector<int>& cyc, int startK, const gsPoint2D<T>& dir) -> bool {
        const int M = static_cast<int>(cyc.size());
        const int vStart = G.halfEdges[cyc[startK]].origin;
        const gsPoint2D<T> P0 = G.nodes[vStart].pt;
        const int nS = 20;

        T bestT = T(1e30);
        int hitK = -1;
        T hitFracHE = T(0.5);
        for (int k = 0; k < M; ++k) {
            if (k == startK || (k + 1) % M == startK) continue;
            CurveSeg<T> cs = getHECurve(cyc[k]);
            gsPoint2D<T> prev = evalCurvePoint(*cs.bsp, T(0));
            for (int s = 1; s <= nS; ++s) {
                gsPoint2D<T> cur = evalCurvePoint(*cs.bsp, static_cast<T>(s) / static_cast<T>(nS));
                T sx = cur.x - prev.x, sy = cur.y - prev.y;
                T den = dir.x * sy - dir.y * sx;
                if (std::abs(den) > T(1e-14)) {
                    T t = ((prev.x - P0.x) * sy - (prev.y - P0.y) * sx) / den;
                    T u = ((prev.x - P0.x) * dir.y - (prev.y - P0.y) * dir.x) / den;
                    if (t > T(1e-4) * diag && t < bestT && u >= -T(1e-3) && u <= T(1) + T(1e-3)) {
                        bestT = t;
                        hitK = k;
                        hitFracHE = (static_cast<T>(s - 1) + std::max(T(0), std::min(T(1), u))) / static_cast<T>(nS);
                    }
                }
                prev = cur;
            }
        }
        if (hitK < 0) return false;
        const auto& heHit = G.halfEdges[cyc[hitK]];
        int vEnd = -1;
        if (hitFracHE < T(0.14) && heHit.origin != vStart && !segExists(vStart, heHit.origin)) {
            vEnd = heHit.origin;
        } else if (hitFracHE > T(0.86) && heHit.target != vStart && !segExists(vStart, heHit.target)) {
            vEnd = heHit.target;
        } else {
            T segFrac = heHit.forward ? hitFracHE : (T(1) - hitFracHE);
            vEnd = splitGraphSegAtFrac(heHit.segIdx, segFrac);
        }
        if (vEnd < 0 || vEnd == vStart || segExists(vStart, vEnd)) return false;
        if (!isSegmentInsideDomain(P0, G.nodes[vEnd].pt, fineMesh.loopPolys)) return false;
        CurveSeg<T> cs = makeStraightSeg(P0, G.nodes[vEnd].pt, degree);
        typename TMeshGraph<T>::Segment ns;
        ns.u = vStart; ns.v = vEnd; ns.seg = cs; ns.isBoundary = false;
        ns.arcLen = P0.dist(G.nodes[vEnd].pt);
        G.segs.push_back(std::move(ns));
        return true;
    };

    // Helper: trace a curved cross-field ray strictly inside face `cyc` from vertex `startK` until it hits
    // another segment of `cyc`, splitting the hit segment (or snapping to its endpoint) and inserting a curved chord!
    auto traceFaceChordFromVertex = [&](const std::vector<int>& cyc, int startK, const gsPoint2D<T>& initDir) -> bool {
        const int M = static_cast<int>(cyc.size());
        const int vStart = G.halfEdges[cyc[startK]].origin;
        const gsPoint2D<T> P0 = G.nodes[vStart].pt;

        const int nS = 20;
        std::vector<std::vector<gsPoint2D<T>>> ePolys(M);
        for (int k = 0; k < M; ++k) {
            CurveSeg<T> cs = getHECurve(cyc[k]);
            ePolys[k].resize(nS + 1);
            for (int s = 0; s <= nS; ++s)
                ePolys[k][s] = evalCurvePoint(*cs.bsp, static_cast<T>(s) / static_cast<T>(nS));
        }

        const T ds = T(0.22) * avgEdge;
        gsPoint2D<T> cur = P0;
        gsPoint2D<T> dir = initDir;
        std::vector<gsPoint2D<T>> poly = {P0};
        int triH = -1;
        T arc = T(0);

        auto commitHit = [&](int hitK, T hitFracHE) -> bool {
            if (hitK < 0 || hitK == startK || (hitK + 1) % M == startK) return false;
            const auto& heHit = G.halfEdges[cyc[hitK]];
            int vEnd = -1;
            if (hitFracHE < T(0.14) && heHit.origin != vStart && !segExists(vStart, heHit.origin)) {
                vEnd = heHit.origin;
            } else if (hitFracHE > T(0.86) && heHit.target != vStart && !segExists(vStart, heHit.target)) {
                vEnd = heHit.target;
            } else {
                T segFrac = heHit.forward ? hitFracHE : (T(1) - hitFracHE);
                vEnd = splitGraphSegAtFrac(heHit.segIdx, segFrac);
            }
            if (vEnd < 0 || vEnd == vStart || segExists(vStart, vEnd)) return false;
            poly.push_back(G.nodes[vEnd].pt);
            CurveSeg<T> fitted = fitStreamlinePolylineToCurveSeg(poly, degree, fineMesh.loopPolys);
            typename TMeshGraph<T>::Segment ns;
            ns.u = vStart; ns.v = vEnd; ns.seg = fitted; ns.isBoundary = false;
            ns.arcLen = measureCurveSegArcLen(fitted);
            G.segs.push_back(std::move(ns));
            return true;
        };

        for (int st = 0; st < 280 && arc < T(1.8) * diag; ++st) {
            if (arc >= T(0.70) * avgEdge) {
                gsPoint2D<T> k1 = dir;
                pickCrossBranch(cur, dir, triH, k1);
                gsPoint2D<T> mid(cur.x + T(0.5) * ds * k1.x, cur.y + T(0.5) * ds * k1.y);
                gsPoint2D<T> k2 = k1;
                if (pickCrossBranch(mid, k1, triH, k2)) dir = k2;
                else dir = k1;
            }
            gsPoint2D<T> nxt(cur.x + ds * dir.x, cur.y + ds * dir.y);

            T bestT = T(2);
            int hitK = -1;
            T hitFracHE = T(0.5);
            T rx = nxt.x - cur.x, ry = nxt.y - cur.y;
            for (int k = 0; k < M; ++k) {
                if (st < 3 && (k == startK || (k + 1) % M == startK)) continue;
                for (int s = 0; s < nS; ++s) {
                    const auto& A = ePolys[k][s];
                    const auto& B = ePolys[k][s + 1];
                    if (st == 0 && (A.dist(P0) < T(1e-4) * diag || B.dist(P0) < T(1e-4) * diag)) continue;
                    T sx = B.x - A.x, sy = B.y - A.y;
                    T den = rx * sy - ry * sx;
                    if (std::abs(den) < T(1e-14)) continue;
                    T t = ((A.x - cur.x) * sy - (A.y - cur.y) * sx) / den;
                    T u = ((A.x - cur.x) * ry - (A.y - cur.y) * rx) / den;
                    if (t > T(1e-4) && t <= T(1) + T(1e-4) && u >= -T(1e-3) && u <= T(1) + T(1e-3)) {
                        if (t < bestT) {
                            bestT = t;
                            hitK = k;
                            hitFracHE = (static_cast<T>(s) + std::max(T(0), std::min(T(1), u))) / static_cast<T>(nS);
                        }
                    }
                }
            }
            if (hitK >= 0) return commitHit(hitK, hitFracHE);

            if (!isPointStrictlyInsideDomain(nxt, fineMesh.loopPolys, T(1e-6))) {
                // Project cur onto nearest non-adjacent edge of cyc instead of dropping!
                int clK = -1; T clFrac = T(0.5); T clD = T(1.5) * ds;
                for (int k = 0; k < M; ++k) {
                    if (k == startK || (k + 1) % M == startK) continue;
                    for (int s = 0; s < nS; ++s) {
                        const auto& A = ePolys[k][s]; const auto& B = ePolys[k][s + 1];
                        T sx = B.x - A.x, sy = B.y - A.y, L2 = sx * sx + sy * sy;
                        T u = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((nxt.x - A.x) * sx + (nxt.y - A.y) * sy) / L2)) : T(0);
                        T d = std::hypot(nxt.x - (A.x + u * sx), nxt.y - (A.y + u * sy));
                        if (d < clD) { clD = d; clK = k; clFrac = (static_cast<T>(s) + u) / static_cast<T>(nS); }
                    }
                }
                return (clK >= 0) ? commitHit(clK, clFrac) : false;
            }
            poly.push_back(nxt);
            arc += ds;
            cur = nxt;
        }
        return false;
    };

    // Classify corners of a raw face cycle (returns corner indices into `cyc` and interior angles at every vertex)
    auto classifyFaceCorners = [&](const std::vector<int>& cyc,
                                   std::vector<int>& cornerPos,
                                   std::vector<T>& vertAngles) {
        const int M = static_cast<int>(cyc.size());
        vertAngles.assign(M, T(180));
        cornerPos.clear();
        for (int k = 0; k < M; ++k) {
            CurveSeg<T> cIn  = getHECurve(cyc[(k + M - 1) % M]);
            CurveSeg<T> cOut = getHECurve(cyc[k]);
            gsPoint2D<T> dIn  = evalCurveDeriv(*cIn.bsp,  T(1));
            gsPoint2D<T> dOut = evalCurveDeriv(*cOut.bsp, T(0));
            T n1 = std::hypot(dIn.x, dIn.y), n2 = std::hypot(dOut.x, dOut.y);
            if (n1 <= T(1e-14) || n2 <= T(1e-14)) { vertAngles[k] = T(90); cornerPos.push_back(k); continue; }
            gsPoint2D<T> tOut(dOut.x / n2, dOut.y / n2);
            gsPoint2D<T> tBack(-dIn.x / n1, -dIn.y / n1);
            T cr = tOut.x * tBack.y - tOut.y * tBack.x;
            T dt = tOut.x * tBack.x + tOut.y * tBack.y;
            T ang = std::atan2(cr, dt) * radToDeg;
            if (ang <= T(0)) ang += T(360);
            vertAngles[k] = ang;
            // True corner if angle deviates from 180 deg by >= 36 deg (i.e., < 144 or > 216)
            if (ang < T(144.0) || ang > T(216.0))
                cornerPos.push_back(k);
        }
    };

    const bool mDbg = (std::getenv("GISMO_CF_DEBUG") != nullptr);

    // Rectify loop: convert ALL raw faces into exact 4-corner non-reflex RectFaces!
    std::vector<std::vector<int>> rawFaces;
    for (int rectIter = 0; rectIter < 120; ++rectIter) {
        if (!buildHalfEdgesAndExtractRawFaces(rawFaces)) return false;

        bool modified = false;
        for (size_t fi = 0; fi < rawFaces.size() && !modified; ++fi) {
            const auto& cyc = rawFaces[fi];
            const int M = static_cast<int>(cyc.size());
            std::vector<int> cPos;
            std::vector<T> vAng;
            classifyFaceCorners(cyc, cPos, vAng);

            // 1. Resolve any reflex corner (> 214 deg) FIRST (straight chord, else curved RK2 cross-field chord!)
            int reflexK = -1;
            T worstReflex = T(214.0);
            for (int k : cPos) {
                if (vAng[k] > worstReflex) { worstReflex = vAng[k]; reflexK = k; }
            }

            if (reflexK >= 0) {
                int vStart = G.halfEdges[cyc[reflexK]].origin;
                gsPoint2D<T> P0 = G.nodes[vStart].pt;
                CurveSeg<T> cOut = getHECurve(cyc[reflexK]);
                gsPoint2D<T> dOut = evalCurveDeriv(*cOut.bsp, T(0));
                T no = std::hypot(dOut.x, dOut.y);
                T angStart = (no > T(1e-12)) ? std::atan2(dOut.y, dOut.x) : T(0);
                T bisAng = angStart + T(0.5) * vAng[reflexK] * pi / T(180);
                gsPoint2D<T> bDir(std::cos(bisAng), std::sin(bisAng));

                std::vector<gsPoint2D<T>> fPoly;
                for (int ch : cyc) {
                    CurveSeg<T> cs = getHECurve(ch);
                    for (int s = 0; s < 6; ++s) fPoly.push_back(evalCurvePoint(*cs.bsp, static_cast<T>(s) / T(6)));
                }
                auto segInsideFace = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B) -> bool {
                    if (A.dist(B) <= T(1e-5) * diag || !isSegmentInsideDomain(A, B, fineMesh.loopPolys)) return false;
                    for (int s = 1; s <= 5; ++s) {
                        T f = static_cast<T>(s) / T(6);
                        if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y), fPoly, T(1e-6)))
                            return false;
                    }
                    const int np = static_cast<int>(fPoly.size());
                    for (int s = 0; s < np; ++s) {
                        const auto& pA = fPoly[s]; const auto& pB = fPoly[(s + 1) % np];
                        if (pA.dist(A) < T(1e-4) * diag || pB.dist(A) < T(1e-4) * diag ||
                            pA.dist(B) < T(1e-4) * diag || pB.dist(B) < T(1e-4) * diag) continue;
                        if (segProperIntersect(A, B, pA, pB)) return false;
                    }
                    return true;
                };

                int bestTargetNode = -1, bestSplitSeg = -1;
                T bestSplitFrac = T(0.5), bestScore = -T(1e9);
                for (int k = 0; k < M; ++k) {
                    if (k == reflexK || (k + 1) % M == reflexK || (reflexK + 1) % M == k) continue;
                    int vT = G.halfEdges[cyc[k]].origin;
                    if (vT == vStart || !segInsideFace(P0, G.nodes[vT].pt)) continue;
                    T d = P0.dist(G.nodes[vT].pt);
                    gsPoint2D<T> uDir((G.nodes[vT].pt.x - P0.x) / d, (G.nodes[vT].pt.y - P0.y) / d);
                    T al = uDir.x * bDir.x + uDir.y * bDir.y;
                    if (al < T(0.10)) continue;
                    T sc = T(2.5) * al - d / diag;
                    if (sc > bestScore) { bestScore = sc; bestTargetNode = vT; bestSplitSeg = -1; }
                }
                for (int k = 0; k < M; ++k) {
                    if (k == reflexK || (k + 1) % M == reflexK) continue;
                    int segId = G.halfEdges[cyc[k]].segIdx;
                    for (T fr : {T(0.5), T(0.35), T(0.65)}) {
                        gsPoint2D<T> PM = evalCurvePoint(*G.segs[segId].seg.bsp, fr);
                        if (!segInsideFace(P0, PM)) continue;
                        T d = P0.dist(PM);
                        gsPoint2D<T> uDir((PM.x - P0.x) / d, (PM.y - P0.y) / d);
                        T al = uDir.x * bDir.x + uDir.y * bDir.y;
                        if (al < T(0.25)) continue;
                        T sc = T(2.0) * al - d / diag - T(0.15);
                        if (sc > bestScore) { bestScore = sc; bestTargetNode = -1; bestSplitSeg = segId; bestSplitFrac = fr; }
                    }
                }

                if (bestTargetNode >= 0 && !segExists(vStart, bestTargetNode)) {
                    CurveSeg<T> cs = makeStraightSeg(P0, G.nodes[bestTargetNode].pt, degree);
                    typename TMeshGraph<T>::Segment ns;
                    ns.u = vStart; ns.v = bestTargetNode; ns.seg = cs; ns.isBoundary = false; ns.arcLen = P0.dist(G.nodes[bestTargetNode].pt);
                    G.segs.push_back(std::move(ns));
                    modified = true;
                } else if (bestSplitSeg >= 0) {
                    int mId = splitGraphSegAtFrac(bestSplitSeg, bestSplitFrac);
                    if (mId >= 0 && !segExists(vStart, mId)) {
                        CurveSeg<T> cs = makeStraightSeg(P0, G.nodes[mId].pt, degree);
                        typename TMeshGraph<T>::Segment ns;
                        ns.u = vStart; ns.v = mId; ns.seg = cs; ns.isBoundary = false; ns.arcLen = P0.dist(G.nodes[mId].pt);
                        G.segs.push_back(std::move(ns));
                        modified = true;
                    }
                }
                if (!modified) {
                    // Curved RK2 cross-field chord fallback, then exact analytical straight ray-cast!
                    for (T frac : {T(0.5), T(0.33), T(0.67), T(0.25), T(0.75)}) {
                        T a = angStart + frac * vAng[reflexK] * pi / T(180);
                        gsPoint2D<T> rDir(std::cos(a), std::sin(a));
                        if (traceFaceChordFromVertex(cyc, reflexK, rDir) || castStraightFaceChord(cyc, reflexK, rDir)) {
                            modified = true;
                            break;
                        }
                    }
                }
                continue;
            }

            const int K = static_cast<int>(cPos.size());
            if (K == 4) continue;

            // 2. K <= 2 (lens/digon) with M < 4: split longest segment(s) so face has >= 4 vertices
            if (K <= 2 && M < 4) {
                int bestSeg = G.halfEdges[cyc[0]].segIdx;
                T bestL = -T(1);
                for (int ch : cyc) {
                    int sid = G.halfEdges[ch].segIdx;
                    if (G.segs[sid].arcLen > bestL) { bestL = G.segs[sid].arcLen; bestSeg = sid; }
                }
                if (splitGraphSegAtFrac(bestSeg, T(0.5)) >= 0) {
                    modified = true;
                    continue;
                }
            }

            // 3. K >= 5 (convex 5-gon, 6-gon, ...): split via straight or curved cross-field chord!
            if (K >= 5) {
                std::vector<gsPoint2D<T>> fPoly;
                for (int ch : cyc) {
                    CurveSeg<T> cs = getHECurve(ch);
                    for (int s = 0; s < 6; ++s) fPoly.push_back(evalCurvePoint(*cs.bsp, static_cast<T>(s) / T(6)));
                }
                auto segInFace = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B) -> bool {
                    if (A.dist(B) <= T(1e-5) * diag || !isSegmentInsideDomain(A, B, fineMesh.loopPolys)) return false;
                    for (int s = 1; s <= 5; ++s) {
                        T f = static_cast<T>(s) / T(6);
                        if (!pointInPolyStrict(gsPoint2D<T>((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y), fPoly, T(1e-6)))
                            return false;
                    }
                    const int np = static_cast<int>(fPoly.size());
                    for (int s = 0; s < np; ++s) {
                        const auto& pA = fPoly[s]; const auto& pB = fPoly[(s + 1) % np];
                        if (pA.dist(A) < T(1e-4) * diag || pB.dist(A) < T(1e-4) * diag ||
                            pA.dist(B) < T(1e-4) * diag || pB.dist(B) < T(1e-4) * diag) continue;
                        if (segProperIntersect(A, B, pA, pB)) return false;
                    }
                    return true;
                };

                if (K >= 6) {
                    int bestU = -1, bestV = -1;
                    T bestSc = -T(1e9);
                    int triH = -1;
                    for (int a = 0; a < K; ++a) {
                        for (int step = 2; step <= K - 2; ++step) {
                            int b = (a + step) % K;
                            if (b <= a) continue;
                            int vA = G.halfEdges[cyc[cPos[a]]].origin;
                            int vB = G.halfEdges[cyc[cPos[b]]].origin;
                            if (vA == vB || !segInFace(G.nodes[vA].pt, G.nodes[vB].pt)) continue;
                            T al = evalSegmentCrossFieldAlignment(fineMesh, cf, accel, G.nodes[vA].pt, G.nodes[vB].pt, triH);
                            T evenBonus = ((step % 2) == 1) ? T(0.6) : T(0);
                            T sc = al + evenBonus - G.nodes[vA].pt.dist(G.nodes[vB].pt) / diag;
                            if (sc > bestSc) { bestSc = sc; bestU = vA; bestV = vB; }
                        }
                    }
                    if (bestU >= 0 && !segExists(bestU, bestV)) {
                        CurveSeg<T> cs = makeStraightSeg(G.nodes[bestU].pt, G.nodes[bestV].pt, degree);
                        typename TMeshGraph<T>::Segment ns;
                        ns.u = bestU; ns.v = bestV; ns.seg = cs; ns.isBoundary = false; ns.arcLen = G.nodes[bestU].pt.dist(G.nodes[bestV].pt);
                        G.segs.push_back(std::move(ns));
                        modified = true;
                        continue;
                    }
                }

                // K == 5 (or fallback for K >= 6): straight chord from corner to opposite side midpoint (1 T-node -> two 4-gons!)
                int bestCornerNode = -1, bestOppSeg = -1;
                T bestSc = -T(1e9);
                int triH = -1;
                for (int a = 0; a < K; ++a) {
                    int vA = G.halfEdges[cyc[cPos[a]]].origin;
                    int oppC0 = cPos[(a + 2) % K];
                    int oppC1 = cPos[(a + 3) % K];
                    int k = oppC0;
                    for (int step = 0; step < M && k != oppC1; ++step, k = (k + 1) % M) {
                        int segId = G.halfEdges[cyc[k]].segIdx;
                        gsPoint2D<T> pM = evalCurvePoint(*G.segs[segId].seg.bsp, T(0.5));
                        if (!segInFace(G.nodes[vA].pt, pM)) continue;
                        T al = evalSegmentCrossFieldAlignment(fineMesh, cf, accel, G.nodes[vA].pt, pM, triH);
                        T sc = al - T(0.5) * G.nodes[vA].pt.dist(pM) / diag;
                        if (sc > bestSc) { bestSc = sc; bestCornerNode = vA; bestOppSeg = segId; }
                    }
                }
                if (bestCornerNode >= 0 && bestOppSeg >= 0) {
                    int mId = splitGraphSegAtFrac(bestOppSeg, T(0.5));
                    if (mId >= 0 && !segExists(bestCornerNode, mId)) {
                        CurveSeg<T> cs = makeStraightSeg(G.nodes[bestCornerNode].pt, G.nodes[mId].pt, degree);
                        typename TMeshGraph<T>::Segment ns;
                        ns.u = bestCornerNode; ns.v = mId; ns.seg = cs; ns.isBoundary = false;
                        ns.arcLen = G.nodes[bestCornerNode].pt.dist(G.nodes[mId].pt);
                        G.segs.push_back(std::move(ns));
                        modified = true;
                        continue;
                    }
                }
                // Curved RK2 + exact analytical ray-cast fallback for strongly bent K >= 5 faces
                for (int a = 0; a < K && !modified; ++a) {
                    int ck = cPos[a];
                    CurveSeg<T> cOut = getHECurve(cyc[ck]);
                    gsPoint2D<T> dOut = evalCurveDeriv(*cOut.bsp, T(0));
                    T no = std::hypot(dOut.x, dOut.y);
                    if (no <= T(1e-12)) continue;
                    T a0 = std::atan2(dOut.y, dOut.x);
                    for (T fr : {T(0.5), T(0.35), T(0.65)}) {
                        T ang = a0 + fr * vAng[ck] * pi / T(180);
                        gsPoint2D<T> rDir(std::cos(ang), std::sin(ang));
                        if (traceFaceChordFromVertex(cyc, ck, rDir) || castStraightFaceChord(cyc, ck, rDir)) {
                            modified = true;
                            break;
                        }
                    }
                }
                if (modified) continue;
            }

            // 4. K == 3 (triangle, WITH OR WITHOUT T-junctions on its 3 sides!)
            if (K == 3) {
                std::vector<int> heToFace(G.halfEdges.size(), -1);
                for (int f2 = 0; f2 < static_cast<int>(rawFaces.size()); ++f2)
                    for (int ch : rawFaces[f2]) heToFace[ch] = f2;

                int killSeg = -1;
                for (int ch : cyc) {
                    int segId = G.halfEdges[ch].segIdx;
                    if (G.segs[segId].isBoundary) continue;
                    int of = heToFace[G.halfEdges[ch].twin];
                    if (of >= 0 && of != static_cast<int>(fi)) {
                        std::vector<int> oC; std::vector<T> oA;
                        classifyFaceCorners(rawFaces[of], oC, oA);
                        if (oC.size() == 3) { killSeg = segId; break; }
                    }
                }
                if (killSeg >= 0) {
                    G.segs.erase(G.segs.begin() + killSeg);
                    modified = true;
                    continue;
                }

                // Local 1-to-3 T-star: reuse existing T-node on any side with cnt >= 2, split single-seg sides at 0.5!
                std::vector<gsPoint2D<T>> fPoly;
                for (int ch : cyc) {
                    CurveSeg<T> cs = getHECurve(ch);
                    for (int s = 0; s < 6; ++s) fPoly.push_back(evalCurvePoint(*cs.bsp, static_cast<T>(s) / T(6)));
                }
                gsPoint2D<T> V3[3], M3[3];
                int existingMidNode[3] = {-1, -1, -1};
                int splitSegId[3] = {-1, -1, -1};
                for (int s = 0; s < 3; ++s) {
                    int c0 = cPos[s], c1 = cPos[(s + 1) % 3];
                    int cnt = (c1 - c0 + M) % M;
                    V3[s] = G.nodes[G.halfEdges[cyc[c0]].origin].pt;
                    if (cnt == 1) {
                        splitSegId[s] = G.halfEdges[cyc[c0]].segIdx;
                        M3[s] = evalCurvePoint(*G.segs[splitSegId[s]].seg.bsp, T(0.5));
                    } else {
                        int midK = (c0 + cnt / 2) % M;
                        existingMidNode[s] = G.halfEdges[cyc[midK]].origin;
                        M3[s] = G.nodes[existingMidNode[s]].pt;
                    }
                }

                // Candidate interior centers O strictly inside fPoly (including inward-normal probes for curved/thin triangles!)
                gsPoint2D<T> O_opt((M3[0].x + M3[1].x + M3[2].x) / T(3), (M3[0].y + M3[1].y + M3[2].y) / T(3));
                pickOptimalTri1to3Center(V3[0], V3[1], V3[2], M3[0], M3[1], M3[2], O_opt);
                std::vector<gsPoint2D<T>> oCands = {
                    O_opt,
                    gsPoint2D<T>((M3[0].x + M3[1].x + M3[2].x) / T(3), (M3[0].y + M3[1].y + M3[2].y) / T(3)),
                    gsPoint2D<T>((V3[0].x + V3[1].x + V3[2].x) / T(3), (V3[0].y + V3[1].y + V3[2].y) / T(3)),
                    polygonCentroid(fPoly)
                };
                for (int s = 0; s < 3; ++s) {
                    gsPoint2D<T> oppV = V3[(s + 2) % 3];
                    for (T w : {T(0.20), T(0.35), T(0.50)})
                        oCands.emplace_back((T(1) - w) * M3[s].x + w * oppV.x, (T(1) - w) * M3[s].y + w * oppV.y);
                }
                gsPoint2D<T> chosenO = O_opt;
                bool foundO = false;
                for (const auto& Oc : oCands) {
                    if (!pointInPolyStrict(Oc, fPoly, T(1e-7)) || !isPointStrictlyInsideDomain(Oc, fineMesh.loopPolys, T(1e-7)))
                        continue;
                    bool spokesOk = true;
                    for (int s = 0; s < 3 && spokesOk; ++s) {
                        for (int k = 1; k <= 3; ++k) {
                            T f = static_cast<T>(k) / T(4);
                            gsPoint2D<T> q((T(1) - f) * Oc.x + f * M3[s].x, (T(1) - f) * Oc.y + f * M3[s].y);
                            if (!pointInPolyStrict(q, fPoly, T(1e-7))) { spokesOk = false; break; }
                        }
                    }
                    if (spokesOk) { chosenO = Oc; foundO = true; break; }
                }
                if (foundO) {
                    int midIds[3] = {existingMidNode[0], existingMidNode[1], existingMidNode[2]};
                    for (int s = 0; s < 3; ++s) {
                        if (midIds[s] < 0 && splitSegId[s] >= 0)
                            midIds[s] = splitGraphSegAtFrac(splitSegId[s], T(0.5));
                    }
                    if (midIds[0] >= 0 && midIds[1] >= 0 && midIds[2] >= 0) {
                        typename TMeshGraph<T>::Node oNd;
                        oNd.pt = chosenO; oNd.isBoundary = false;
                        int oId = static_cast<int>(G.nodes.size());
                        G.nodes.push_back(oNd);
                        for (int s = 0; s < 3; ++s) {
                            CurveSeg<T> spk = makeStraightSeg(chosenO, G.nodes[midIds[s]].pt, degree);
                            typename TMeshGraph<T>::Segment ns;
                            ns.u = oId; ns.v = midIds[s]; ns.seg = spk; ns.isBoundary = false;
                            ns.arcLen = chosenO.dist(G.nodes[midIds[s]].pt);
                            G.segs.push_back(std::move(ns));
                        }
                        modified = true;
                        continue;
                    }
                }
                // Fallback for razor-thin curved sliver triangle: dissolve longest interior edge (3+4 -> 5, then 5 -> 2 quads!)
                int dissolveSeg = -1;
                T bestLenD = -T(1);
                for (int ch : cyc) {
                    int segId = G.halfEdges[ch].segIdx;
                    if (G.segs[segId].isBoundary) continue;
                    int of = heToFace[G.halfEdges[ch].twin];
                    if (of >= 0 && of != static_cast<int>(fi) && G.segs[segId].arcLen > bestLenD) {
                        bestLenD = G.segs[segId].arcLen;
                        dissolveSeg = segId;
                    }
                }
                if (dissolveSeg >= 0) {
                    G.segs.erase(G.segs.begin() + dissolveSeg);
                    modified = true;
                    continue;
                }
            }
        }
        if (!modified) break;
    }

    // Build final RectFaces from rectified rawFaces
    if (!buildHalfEdgesAndExtractRawFaces(rawFaces)) {
        if (mDbg) std::cerr << "[CF_DBG] MotoFail: rawFaces empty\n";
        return false;
    }
    // Verify 100% forward boundary half-edge coverage (catches any unbridged hole or dropped face early!)
    {
        std::vector<bool> heUsed(G.halfEdges.size(), false);
        for (const auto& cyc : rawFaces)
            for (int ch : cyc) heUsed[ch] = true;
        int missBnd = 0;
        for (size_t h = 0; h < G.halfEdges.size(); ++h)
            if (G.halfEdges[h].forward && G.segs[G.halfEdges[h].segIdx].isBoundary && !heUsed[h])
                ++missBnd;
        if (missBnd > 0) {
            if (mDbg) std::cerr << "[CF_DBG] MotoFail: missing " << missBnd << " forward boundary HEs\n";
            return false;
        }
    }

    G.faces.clear();
    for (size_t fi = 0; fi < rawFaces.size(); ++fi) {
        const auto& cyc = rawFaces[fi];
        const int M = static_cast<int>(cyc.size());
        std::vector<int> cPos;
        std::vector<T> vAng;
        classifyFaceCorners(cyc, cPos, vAng);

        if (cPos.size() < 4 && M >= 4) {
            std::vector<std::pair<T, int>> cand;
            std::set<int> usedC(cPos.begin(), cPos.end());
            for (int k = 0; k < M; ++k)
                if (!usedC.count(k)) cand.push_back({std::abs(vAng[k] - T(180)), k});
            std::sort(cand.begin(), cand.end(), [](const std::pair<T, int>& a, const std::pair<T, int>& b) { return a.first > b.first; });
            for (const auto& cp : cand) {
                if (cPos.size() == 4) break;
                cPos.push_back(cp.second);
            }
            std::sort(cPos.begin(), cPos.end());
        }
        if (cPos.size() > 4) {
            std::vector<std::pair<T, int>> scored;
            for (int k : cPos) scored.push_back({std::abs(vAng[k] - T(90)), k});
            std::sort(scored.begin(), scored.end());
            bool mildDrop = true;
            for (size_t i = 4; i < scored.size(); ++i)
                if (vAng[scored[i].second] < T(115.0) || vAng[scored[i].second] > T(240.0)) mildDrop = false;
            if (!mildDrop) {
                if (mDbg) {
                    std::cerr << "[CF_DBG] MotoFail: face " << fi << "/" << rawFaces.size()
                              << " M=" << M << " K=" << cPos.size() << " angs=[";
                    for (int k : cPos) std::cerr << vAng[k] << " ";
                    std::cerr << "]\n";
                }
                return false;
            }
            cPos.clear();
            for (int i = 0; i < 4; ++i) cPos.push_back(scored[i].second);
            std::sort(cPos.begin(), cPos.end());
        }
        if (cPos.size() != 4) {
            if (mDbg) {
                std::cerr << "[CF_DBG] MotoFail: face " << fi << "/" << rawFaces.size()
                          << " M=" << M << " K=" << cPos.size() << "\n";
            }
            return false;
        }

        typename TMeshGraph<T>::RectFace rf;
        T sideLen[4] = {T(0), T(0), T(0), T(0)};
        for (int s = 0; s < 4; ++s) {
            int k0 = cPos[s], k1 = cPos[(s + 1) % 4];
            rf.corners[s] = G.halfEdges[cyc[k0]].origin;
            int k = k0;
            for (int guard = 0; guard < M && k != k1; ++guard, k = (k + 1) % M) {
                rf.sideHEs[s].push_back(cyc[k]);
                sideLen[s] += G.segs[G.halfEdges[cyc[k]].segIdx].arcLen;
            }
            if (rf.sideHEs[s].empty()) return false;
        }
        rf.Lu = std::max(T(1e-6), T(0.5) * (sideLen[0] + sideLen[2]));
        rf.Lv = std::max(T(1e-6), T(0.5) * (sideLen[1] + sideLen[3]));
        G.faces.push_back(std::move(rf));
    }

    if (G.faces.empty()) return false;
    outTMesh = std::move(G);
    return true;
}

// ============================================================================
// Stage 3: Native Constrained T-Mesh Integer Quantization + Macro-Block Emission
// ============================================================================

inline int64_t tmeshGcd64(int64_t a, int64_t b)
{
    a = (a < 0) ? -a : a;
    b = (b < 0) ? -b : b;
    while (b != 0) { int64_t t = a % b; a = b; b = t; }
    return (a == 0) ? 1 : a;
}

struct TMeshRat {
    int64_t n = 0, d = 1;
    TMeshRat() = default;
    explicit TMeshRat(int64_t num, int64_t den = 1) : n(num), d(den) { norm(); }
    void norm() {
        if (d < 0) { n = -n; d = -d; }
        if (n == 0) { d = 1; return; }
        int64_t g = tmeshGcd64(n, d);
        n /= g; d /= g;
    }
    bool isZero() const { return n == 0; }
    double toDouble() const { return static_cast<double>(n) / static_cast<double>(d); }
    TMeshRat operator+(const TMeshRat& o) const { return TMeshRat(n * o.d + o.n * d, d * o.d); }
    TMeshRat operator-(const TMeshRat& o) const { return TMeshRat(n * o.d - o.n * d, d * o.d); }
    TMeshRat operator*(const TMeshRat& o) const { return TMeshRat(n * o.n, d * o.d); }
    TMeshRat operator/(const TMeshRat& o) const { return TMeshRat(n * o.d, d * o.n); }
    TMeshRat operator-() const { return TMeshRat(-n, d); }
};

template <class T>
inline bool quantizeTMeshMinimalPatches(
    TMeshGraph<T>& tmesh,
    T lambdaAR,
    int maxPatchesCap)
{
    const int F = static_cast<int>(tmesh.faces.size());
    const int E = static_cast<int>(tmesh.segs.size());
    if (F == 0 || E == 0) return false;

    // Build exact opposite-side incidence matrix C (2F x E)
    // Row 2*f:     sum_{e in South} n_e - sum_{e in North} n_e = 0
    // Row 2*f + 1: sum_{e in East}  n_e - sum_{e in West}  n_e = 0
    const int R = 2 * F;
    std::vector<std::vector<double>> M(R, std::vector<double>(E, 0.0));
    std::vector<std::array<std::vector<int>, 4>> faceSegs(F);
    for (int f = 0; f < F; ++f) {
        for (int s = 0; s < 4; ++s) {
            for (int he : tmesh.faces[f].sideHEs[s]) {
                int eid = tmesh.halfEdges[he].segIdx;
                faceSegs[f][s].push_back(eid);
            }
        }
        for (int eid : faceSegs[f][0]) M[2 * f][eid] += 1.0;
        for (int eid : faceSegs[f][2]) M[2 * f][eid] -= 1.0;
        for (int eid : faceSegs[f][1]) M[2 * f + 1][eid] += 1.0;
        for (int eid : faceSegs[f][3]) M[2 * f + 1][eid] -= 1.0;
    }

    // Exact rational Gauss-Jordan RREF (with safe double fallback if any denominator explodes)
    bool useExact = true;
    std::vector<std::vector<TMeshRat>> MR(R, std::vector<TMeshRat>(E));
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < E; ++c)
            MR[r][c] = TMeshRat(static_cast<int64_t>(std::llround(M[r][c])));

    std::vector<int> pivCols, freeCols;
    int pivRow = 0;
    for (int col = 0; col < E && pivRow < R && useExact; ++col) {
        int sel = -1;
        for (int r = pivRow; r < R; ++r)
            if (!MR[r][col].isZero()) { sel = r; break; }
        if (sel < 0) continue;
        std::swap(MR[pivRow], MR[sel]);
        TMeshRat pv = MR[pivRow][col];
        for (int c = col; c < E; ++c) {
            MR[pivRow][c] = MR[pivRow][c] / pv;
            if (std::abs(MR[pivRow][c].n) > 1000000000LL || MR[pivRow][c].d > 1000000000LL) { useExact = false; break; }
        }
        for (int r = 0; r < R && useExact; ++r) {
            if (r == pivRow || MR[r][col].isZero()) continue;
            TMeshRat fac = MR[r][col];
            for (int c = col; c < E; ++c) {
                MR[r][c] = MR[r][c] - fac * MR[pivRow][c];
                if (std::abs(MR[r][c].n) > 1000000000LL || MR[r][c].d > 1000000000LL) { useExact = false; break; }
            }
        }
        pivCols.push_back(col);
        ++pivRow;
    }

    std::vector<std::vector<double>> B; // n_piv[i] = sum_j B[i][j] * n_free[j]
    int nPiv = 0, d = 0;

    if (useExact) {
        std::vector<bool> isPiv(E, false);
        for (int c : pivCols) isPiv[c] = true;
        for (int c = 0; c < E; ++c) if (!isPiv[c]) freeCols.push_back(c);
        nPiv = static_cast<int>(pivCols.size());
        d = static_cast<int>(freeCols.size());
        B.assign(nPiv, std::vector<double>(d, 0.0));
        for (int i = 0; i < nPiv; ++i)
            for (int j = 0; j < d; ++j)
                B[i][j] = (-MR[i][freeCols[j]]).toDouble();
    } else {
        // Double-precision RREF fallback with partial pivoting
        pivCols.clear(); freeCols.clear();
        std::vector<std::vector<double>> MD = M;
        pivRow = 0;
        for (int col = 0; col < E && pivRow < R; ++col) {
            int sel = -1;
            double bestAbs = 1e-9;
            for (int r = pivRow; r < R; ++r)
                if (std::abs(MD[r][col]) > bestAbs) { bestAbs = std::abs(MD[r][col]); sel = r; }
            if (sel < 0) continue;
            std::swap(MD[pivRow], MD[sel]);
            double pv = MD[pivRow][col];
            for (int c = col; c < E; ++c) MD[pivRow][c] /= pv;
            for (int r = 0; r < R; ++r) {
                if (r == pivRow || std::abs(MD[r][col]) < 1e-12) continue;
                double fac = MD[r][col];
                for (int c = col; c < E; ++c) MD[r][c] -= fac * MD[pivRow][c];
            }
            pivCols.push_back(col);
            ++pivRow;
        }
        std::vector<bool> isPiv(E, false);
        for (int c : pivCols) isPiv[c] = true;
        for (int c = 0; c < E; ++c) if (!isPiv[c]) freeCols.push_back(c);
        nPiv = static_cast<int>(pivCols.size());
        d = static_cast<int>(freeCols.size());
        B.assign(nPiv, std::vector<double>(d, 0.0));
        for (int i = 0; i < nPiv; ++i)
            for (int j = 0; j < d; ++j)
                B[i][j] = -MD[i][freeCols[j]];
    }

    // Immediate infeasibility check (spiral T-cycle where a pivot has all B[i][j] <= 0)
    for (int i = 0; i < nPiv; ++i) {
        bool hasPos = false;
        for (int j = 0; j < d; ++j) if (B[i][j] > 1e-7) { hasPos = true; break; }
        if (!hasPos) return false;
    }

    // Evaluate full integer assignment from free variables, verifying exact C * n == 0 and n >= 1
    auto evalFreeVec = [&](const std::vector<int>& nFree,
                           std::vector<int>& nFull,
                           int& totalPatchesOut,
                           double& costOut) -> bool {
        nFull.assign(E, 0);
        for (int j = 0; j < d; ++j) {
            if (nFree[j] < 1) return false;
            nFull[freeCols[j]] = nFree[j];
        }
        for (int i = 0; i < nPiv; ++i) {
            double val = 0.0;
            for (int j = 0; j < d; ++j) val += B[i][j] * static_cast<double>(nFree[j]);
            int iv = static_cast<int>(std::llround(val));
            if (iv < 1 || std::abs(val - static_cast<double>(iv)) > 1e-6) return false;
            nFull[pivCols[i]] = iv;
        }
        // Double-check exact opposite-side equality on every face
        int totP = 0;
        double arCost = 0.0;
        for (int f = 0; f < F; ++f) {
            int ns = 0, nn = 0, ne = 0, nw = 0;
            for (int e : faceSegs[f][0]) ns += nFull[e];
            for (int e : faceSegs[f][2]) nn += nFull[e];
            for (int e : faceSegs[f][1]) ne += nFull[e];
            for (int e : faceSegs[f][3]) nw += nFull[e];
            if (ns != nn || ne != nw || ns < 1 || ne < 1) return false;
            totP += ns * ne;
            double r = (static_cast<double>(ns) * static_cast<double>(tmesh.faces[f].Lv)) /
                       std::max(1e-12, static_cast<double>(ne) * static_cast<double>(tmesh.faces[f].Lu));
            arCost += (r + 1.0 / std::max(1e-9, r) - 2.0);
        }
        totalPatchesOut = totP;
        costOut = static_cast<double>(totP) + static_cast<double>(lambdaAR) * arCost;
        return true;
    };

    // Determine exact LCM multiplier for rational entries in B (usually 1 or 2)
    int64_t lcmDen = 1;
    if (useExact) {
        for (int i = 0; i < nPiv; ++i) {
            for (int j = 0; j < d; ++j) {
                int64_t den = MR[i][freeCols[j]].d;
                if (den > 1 && lcmDen <= 12)
                    lcmDen = (lcmDen / tmeshGcd64(lcmDen, den)) * den;
            }
        }
    }

    // Step 3b: Find feasible initial integer seed via greedy feasibility repair + continuous relaxation
    std::vector<int> bestFree(d, 1), bestFull;
    double bestCost = 1e30;
    int bestPatches = 1000000;
    bool haveBest = false;

    auto registerSeed = [&](const std::vector<int>& candFree) {
        std::vector<int> full; int totP = 0; double c = 0.0;
        if (evalFreeVec(candFree, full, totP, c) && c < bestCost) {
            bestCost = c; bestPatches = totP; bestFree = candFree; bestFull = full; haveBest = true;
        }
    };

    auto greedyRepairFree = [&](std::vector<int> curFree) {
        for (int it = 0; it < 200; ++it) {
            int worstPiv = -1;
            double worstVal = 1.0 - 1e-6;
            for (int i = 0; i < nPiv; ++i) {
                double val = 0.0;
                for (int j = 0; j < d; ++j) val += B[i][j] * static_cast<double>(curFree[j]);
                if (val < worstVal) { worstVal = val; worstPiv = i; }
            }
            if (worstPiv < 0) break;
            // Increment free variable with largest positive contribution to worstPiv and fewest negative side-effects
            int bestJ = -1;
            double bestGain = -1e9;
            for (int j = 0; j < d; ++j) {
                if (B[worstPiv][j] <= 1e-7) continue;
                double gain = 4.0 * B[worstPiv][j];
                for (int i = 0; i < nPiv; ++i)
                    if (B[i][j] < -1e-7) gain += B[i][j];
                if (gain > bestGain) { bestGain = gain; bestJ = j; }
            }
            if (bestJ < 0 || curFree[bestJ] >= 12) break;
            ++curFree[bestJ];
        }
        registerSeed(curFree);
        if (lcmDen > 1 && lcmDen <= 6) {
            std::vector<int> scaled = curFree;
            for (int& v : scaled) v *= static_cast<int>(lcmDen);
            registerSeed(scaled);
        }
    };

    // Seed 1: uniform 1s repaired
    greedyRepairFree(std::vector<int>(d, 1));

    // Seed 2: physical-length-proportional seed (handles elongated channels nicely!)
    {
        std::vector<double> lens(E);
        for (int e = 0; e < E; ++e) lens[e] = std::max(1e-9, static_cast<double>(tmesh.segs[e].arcLen));
        std::vector<double> sortedL = lens;
        std::sort(sortedL.begin(), sortedL.end());
        double refL = sortedL[sortedL.size() / 4]; // 25th percentile length
        for (double scale : {1.0, 1.5, 2.2}) {
            std::vector<int> sFree(d, 1);
            for (int j = 0; j < d; ++j)
                sFree[j] = std::max(1, std::min(6, static_cast<int>(std::llround(lens[freeCols[j]] / std::max(1e-9, scale * refL)))));
            greedyRepairFree(sFree);
        }
    }

    // Step 3c: Bounded Branch-and-Bound for small nullspace dimension (d <= 14)
    if (d > 0 && d <= 14) {
        std::vector<int> uBound(d, 4);
        for (int j = 0; j < d; ++j) {
            int ub = haveBest ? std::min(6, std::max(3, bestFree[j] + 2)) : 5;
            uBound[j] = ub;
        }
        // Order free variables by coupling count
        std::vector<int> order(d);
        std::iota(order.begin(), order.end(), 0);
        std::vector<int> coup(d, 0);
        for (int i = 0; i < nPiv; ++i)
            for (int j = 0; j < d; ++j)
                if (std::abs(B[i][j]) > 1e-7) ++coup[j];
        std::sort(order.begin(), order.end(), [&](int a, int b) { return coup[a] > coup[b]; });

        // Precompute remaining max positive contribution per pivot at each recursion depth
        std::vector<std::vector<double>> remMax(d + 1, std::vector<double>(nPiv, 0.0));
        for (int dep = d - 1; dep >= 0; --dep) {
            int j = order[dep];
            for (int i = 0; i < nPiv; ++i) {
                double add = (B[i][j] > 0.0) ? (B[i][j] * uBound[j]) : (B[i][j] * 1.0);
                remMax[dep][i] = remMax[dep + 1][i] + add;
            }
        }

        std::vector<int> curFree(d, 1);
        std::vector<double> curPiv(nPiv, 0.0);
        int nodeCount = 0;
        const int maxNodes = 6000;

        std::function<void(int)> dfs = [&](int dep) {
            if (++nodeCount > maxNodes) return;
            if (dep == d) {
                registerSeed(curFree);
                return;
            }
            int j = order[dep];
            for (int val = 1; val <= uBound[j]; ++val) {
                curFree[j] = val;
                bool feasible = true;
                for (int i = 0; i < nPiv; ++i) {
                    curPiv[i] += B[i][j] * val;
                    if (curPiv[i] + remMax[dep + 1][i] < 1.0 - 1e-5) feasible = false;
                }
                if (feasible) dfs(dep + 1);
                for (int i = 0; i < nPiv; ++i) curPiv[i] -= B[i][j] * val;
            }
        };
        dfs(0);
    }

    if (!haveBest) return false;

    // Step 3d: Discrete Coordinate + Pairwise Local Search (refines bestFree for any d!)
    for (int iter = 0; iter < 48; ++iter) {
        bool improved = false;
        // Single-variable +-1 moves
        for (int j = 0; j < d; ++j) {
            for (int delta : {-1, 1}) {
                if (bestFree[j] + delta < 1 || bestFree[j] + delta > 10) continue;
                std::vector<int> cand = bestFree;
                cand[j] += delta;
                std::vector<int> full; int totP = 0; double c = 0.0;
                if (evalFreeVec(cand, full, totP, c) && c < bestCost - 1e-9) {
                    bestFree = cand; bestFull = full; bestPatches = totP; bestCost = c;
                    improved = true;
                }
            }
        }
        // Pairwise (+1, -1) or (-1, -1) moves (essential for coupled T-junctions!)
        if (!improved && d <= 48) {
            for (int j = 0; j < d && !improved; ++j) {
                for (int k = j + 1; k < d && !improved; ++k) {
                    const int deltas[4][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}};
                    for (int m = 0; m < 4; ++m) {
                        int nj = bestFree[j] + deltas[m][0];
                        int nk = bestFree[k] + deltas[m][1];
                        if (nj < 1 || nj > 10 || nk < 1 || nk > 10) continue;
                        std::vector<int> cand = bestFree;
                        cand[j] = nj; cand[k] = nk;
                        std::vector<int> full; int totP = 0; double c = 0.0;
                        if (evalFreeVec(cand, full, totP, c) && c < bestCost - 1e-9) {
                            bestFree = cand; bestFull = full; bestPatches = totP; bestCost = c;
                            improved = true;
                            break;
                        }
                    }
                }
            }
        }
        if (!improved) break;
    }

    if (bestPatches <= 0 || bestPatches > maxPatchesCap + 40)
        return false;

    for (int e = 0; e < E; ++e)
        tmesh.segs[e].quantLen = bestFull[e];
    return true;
}

template <class T>
inline void splitMacroPatchUniformGrid(
    const gsTensorBSpline<2, T>& macroTb,
    int Nu,
    int Nv,
    gsMultiPatch<T>& outMp)
{
    auto addNorm = [&](const gsTensorBSpline<2, T>& bsp) {
        auto* cp = new gsTensorBSpline<2, T>(bsp);
        safeNormalizeKnotVector01(cp->knots(0));
        safeNormalizeKnotVector01(cp->knots(1));
        outMp.addPatch(typename gsGeometry<T>::uPtr(cp));
    };

    gsTensorBSpline<2, T> base(macroTb);
    safeNormalizeKnotVector01(base.knots(0));
    safeNormalizeKnotVector01(base.knots(1));

    if (Nu <= 1 && Nv <= 1) {
        addNorm(base);
        return;
    }

    std::vector<gsTensorBSpline<2, T>> uStrips;
    gsTensorBSpline<2, T> remU = base;
    for (int i = 1; i < Nu; ++i) {
        T xi = static_cast<T>(i) / static_cast<T>(Nu);
        gsTensorBSpline<2, T> left, right;
        remU.splitAt(0, xi, left, right);
        uStrips.push_back(std::move(left));
        remU = std::move(right);
    }
    uStrips.push_back(std::move(remU));

    for (auto& strip : uStrips) {
        gsTensorBSpline<2, T> remV = strip;
        for (int j = 1; j < Nv; ++j) {
            T eta = static_cast<T>(j) / static_cast<T>(Nv);
            gsTensorBSpline<2, T> bot, top;
            remV.splitAt(1, eta, bot, top);
            addNorm(bot);
            remV = std::move(top);
        }
        addNorm(remV);
    }
}

template <class T>
inline bool emitQuantizedTMeshPatches(
    const TMeshGraph<T>& tmesh,
    const std::vector<VirtualCADLoop<T>>& vloops,
    int degree,
    gsMultiPatch<T>& outMp)
{
    outMp.clear();
    const int E = static_cast<int>(tmesh.segs.size());
    const int F = static_cast<int>(tmesh.faces.size());
    if (E == 0 || F == 0) return false;

    // 1. Pre-split every undirected segment e into exact quantLen canonical sub-curves (u -> v)
    std::vector<std::vector<CurveSeg<T>>> segSubs(E);
    for (int e = 0; e < E; ++e) {
        const auto& sg = tmesh.segs[e];
        const int ne = std::max(1, sg.quantLen);
        if (ne == 1) {
            segSubs[e].push_back(sg.seg);
        } else if (sg.isBoundary && sg.loopIdx >= 0 && sg.loopIdx < static_cast<int>(vloops.size())) {
            // Exact CAD slicing along VirtualCADLoop (exact endpoints sg.s0 and sg.s1!)
            const auto& vl = vloops[sg.loopIdx];
            for (int k = 0; k < ne; ++k) {
                T sa = (k == 0) ? sg.s0 : virtualLoopWrapS(vl, sg.s0 + sg.arcLen * (static_cast<T>(k) / static_cast<T>(ne)));
                T sb = (k + 1 == ne) ? sg.s1 : virtualLoopWrapS(vl, sg.s0 + sg.arcLen * (static_cast<T>(k + 1) / static_cast<T>(ne)));
                CurveSeg<T> piece = extractVirtualSubCurve(vl, sa, sb, degree);
                if (!piece.bsp) return false;
                segSubs[e].push_back(piece);
            }
        } else {
            // Equal arc-length exact B-spline splitting of interior curve
            std::vector<T> ts;
            ts.reserve(ne - 1);
            for (int k = 1; k < ne; ++k)
                ts.push_back(arcLengthSplitParam(*sg.seg.bsp, static_cast<T>(k) / static_cast<T>(ne)));
            segSubs[e] = splitSegMulti(sg.seg, ts);
            if (static_cast<int>(segSubs[e].size()) != ne) return false;
        }
    }

    auto collectSideSubCurves = [&](const std::vector<int>& sideHEs, bool reverseEntireSide) -> std::vector<CurveSeg<T>> {
        std::vector<CurveSeg<T>> out;
        for (int heIdx : sideHEs) {
            const auto& he = tmesh.halfEdges[heIdx];
            const auto& subs = segSubs[he.segIdx];
            if (he.forward) {
                for (const auto& s : subs) out.push_back(s);
            } else {
                for (int k = static_cast<int>(subs.size()) - 1; k >= 0; --k)
                    out.push_back(subs[k].reversed());
            }
        }
        if (reverseEntireSide) {
            std::reverse(out.begin(), out.end());
            for (auto& s : out) s = s.reversed();
        }
        return out;
    };

    // 2. Emit each quantized RectFace (1x1 directly, Nu x Nv via macro-block smoothing + exact uniform knot split!)
    for (int f = 0; f < F; ++f) {
        const auto& rf = tmesh.faces[f];
        std::vector<CurveSeg<T>> cS_vec = collectSideSubCurves(rf.sideHEs[0], false); // C0 -> C1
        std::vector<CurveSeg<T>> cE_vec = collectSideSubCurves(rf.sideHEs[1], false); // C1 -> C2
        std::vector<CurveSeg<T>> cN_vec = collectSideSubCurves(rf.sideHEs[2], true);  // C3 -> C2 (increasing u)
        std::vector<CurveSeg<T>> cW_vec = collectSideSubCurves(rf.sideHEs[3], true);  // C0 -> C3 (increasing v)

        const int Nu = static_cast<int>(cS_vec.size());
        const int Nv = static_cast<int>(cE_vec.size());
        if (Nu == 0 || Nv == 0 || static_cast<int>(cN_vec.size()) != Nu || static_cast<int>(cW_vec.size()) != Nv)
            return false;

        CurveSeg<T> mS = joinCurveSegsUniformParam(cS_vec, degree);
        CurveSeg<T> mE = joinCurveSegsUniformParam(cE_vec, degree);
        CurveSeg<T> mN = joinCurveSegsUniformParam(cN_vec, degree);
        CurveSeg<T> mW = joinCurveSegsUniformParam(cW_vec, degree);
        if (!mS.bsp || !mE.bsp || !mN.bsp || !mW.bsp) return false;

        std::unique_ptr<gsTensorBSpline<2, T>> macroTb(
            makeCoonsPatchFrom4Curves(*mS.bsp, *mE.bsp, *mN.bsp, *mW.bsp, degree));
        if (!macroTb) return false;

        if (evalPatchMinDetDense(*macroTb, 33) <= T(1e-7)) {
            quickHarmonicPatchRepair(*macroTb, T(1e-7));
            if (evalPatchMinDetDense(*macroTb, 33) <= T(1e-7)) {
                try {
                    gsMultiPatch<T> bndC;
                    bndC.addPatch(macroTb->boundary(boundary::south));
                    bndC.addPatch(macroTb->boundary(boundary::east));
                    bndC.addPatch(macroTb->boundary(boundary::north));
                    bndC.addPatch(macroTb->boundary(boundary::west));
                    gsSpringPatch<T> sp(bndC);
                    const gsGeometry<T>& res = sp.compute();
                    if (const auto* sTb = dynamic_cast<const gsTensorBSpline<2, T>*>(&res)) {
                        if (evalPatchMinDetDense(*sTb, 33) > evalPatchMinDetDense(*macroTb, 33))
                            *macroTb = *sTb;
                    }
                } catch (...) {}
            }
        }

        splitMacroPatchUniformGrid(*macroTb, Nu, Nv, outMp);
    }

    return outMp.nPatches() > 0;
}

} // namespace internal_splitter
} // namespace gismo
