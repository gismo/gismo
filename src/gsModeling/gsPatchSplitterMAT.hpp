/** @file gsPatchSplitterMAT.hpp
    @brief 2D Delaunay-Dual Medial Axis Transform (MAT) Quad Decomposition for gsPatchSplitter.

    Constructs a 2D Medial Axis Transform from a boundary-only constrained Delaunay
    triangulation (CDT), prunes spurious branches via theta-SMA object-angle filtering,
    contracts clustered junctions, smooths medial branches with bend-station insertion,
    and emits an all-quad conforming gsMultiPatch<T> preserving exact CAD boundaries.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova
*/

#pragma once

#include <gsModeling/gsPatchSplitterGeom.hpp>
#include <gsModeling/gsPatchSplitterTriMesh.hpp>
#include <gsModeling/gsPatchSplitterTemplates.hpp>
#include <gsCore/gsMultiPatch.h>
#include <gsModeling/gsPlanarDomain.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <memory>
#include <iostream>
#include <fstream>
#include <queue>
#include <utility>

namespace gismo {
namespace internal_splitter {

// ============================================================================
// Legacy VTK Debug Exporters (Guarded by GISMO_MAT_DEBUG)
// ============================================================================

template <class T>
inline void exportCDTToVTK(const TriMeshData<T>& mesh, const std::vector<int>& triType, const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) return;
    int nP = static_cast<int>(mesh.pts.size());
    int nT = static_cast<int>(mesh.tris.size());
    out << "# vtk DataFile Version 3.0\nCDT Triangles\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << nP << " float\n";
    for (const auto& p : mesh.pts)
        out << static_cast<float>(p.x) << " " << static_cast<float>(p.y) << " 0\n";
    out << "CELLS " << nT << " " << (4 * nT) << "\n";
    for (const auto& tr : mesh.tris)
        out << "3 " << tr.v[0] << " " << tr.v[1] << " " << tr.v[2] << "\n";
    out << "CELL_TYPES " << nT << "\n";
    for (int i = 0; i < nT; ++i) out << "5\n";
    out << "CELL_DATA " << nT << "\nSCALARS triType int 1\nLOOKUP_TABLE default\n";
    for (int t : triType) out << t << "\n";
}

template <class T>
inline void exportSkeletonToVTK(
    const std::vector<std::vector<gsPoint2D<T>>>& polylines,
    const std::vector<T>& polylineRadii,
    const std::vector<gsPoint2D<T>>& junctions,
    const std::vector<T>& junctionRadii,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) return;
    std::vector<gsPoint2D<T>> allPts;
    std::vector<std::vector<int>> cellConn;
    std::vector<int> branchIds;
    std::vector<float> cellRadii;

    int curBranchId = 0;
    for (size_t i = 0; i < polylines.size(); ++i) {
        const auto& poly = polylines[i];
        if (poly.size() < 2) continue;
        std::vector<int> conn;
        for (const auto& pt : poly) {
            conn.push_back(static_cast<int>(allPts.size()));
            allPts.push_back(pt);
        }
        cellConn.push_back(conn);
        branchIds.push_back(curBranchId++);
        cellRadii.push_back(i < polylineRadii.size() ? static_cast<float>(polylineRadii[i]) : 0.0f);
    }
    for (size_t i = 0; i < junctions.size(); ++i) {
        std::vector<int> conn = {static_cast<int>(allPts.size())};
        allPts.push_back(junctions[i]);
        cellConn.push_back(conn);
        branchIds.push_back(curBranchId++);
        cellRadii.push_back(i < junctionRadii.size() ? static_cast<float>(junctionRadii[i]) : 0.0f);
    }
    int totalCellEntries = 0;
    for (const auto& c : cellConn) totalCellEntries += 1 + static_cast<int>(c.size());

    out << "# vtk DataFile Version 3.0\nMedial Skeleton\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << allPts.size() << " float\n";
    for (const auto& p : allPts)
        out << static_cast<float>(p.x) << " " << static_cast<float>(p.y) << " 0\n";
    out << "CELLS " << cellConn.size() << " " << totalCellEntries << "\n";
    for (const auto& c : cellConn) {
        out << c.size();
        for (int idx : c) out << " " << idx;
        out << "\n";
    }
    out << "CELL_TYPES " << cellConn.size() << "\n";
    for (const auto& c : cellConn) out << (c.size() == 1 ? "1\n" : "4\n");
    out << "CELL_DATA " << cellConn.size() << "\n";
    out << "SCALARS branchId int 1\nLOOKUP_TABLE default\n";
    for (int bid : branchIds) out << bid << "\n";
    out << "SCALARS radius float 1\nLOOKUP_TABLE default\n";
    for (float rad : cellRadii) out << rad << "\n";
}

template <class T>
inline void exportSpokesToVTK(
    const std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>>& spokes,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) return;
    int nS = static_cast<int>(spokes.size());
    out << "# vtk DataFile Version 3.0\nMedial Spokes\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << (2 * nS) << " float\n";
    for (const auto& sp : spokes) {
        out << static_cast<float>(sp.first.x) << " " << static_cast<float>(sp.first.y) << " 0\n";
        out << static_cast<float>(sp.second.x) << " " << static_cast<float>(sp.second.y) << " 0\n";
    }
    out << "CELLS " << nS << " " << (3 * nS) << "\n";
    for (int i = 0; i < nS; ++i)
        out << "2 " << (2 * i) << " " << (2 * i + 1) << "\n";
    out << "CELL_TYPES " << nS << "\n";
    for (int i = 0; i < nS; ++i) out << "3\n";
}

// ============================================================================
// 2D Medial Axis Transform (MAT) Quad Decomposition Engine
// ============================================================================

template <class T>
class gsPatchSplitterMAT {
public:
    struct MATNode {
        gsPoint2D<T> pt;
        T radius = 0;
        int triIdx = -1;
        int triType = 1; // 0: T (leaf), 1: S (sleeve), 2: J (junction)
        std::vector<int> adj; // adjacent CDT triangle indices across interior edges
        bool active = true;
    };

    struct Chord {
        int vL = -1; // boundary vertex on LEFT when moving forward along branch
        int vR = -1; // boundary vertex on RIGHT when moving forward along branch
    };

    struct MacroBranch {
        int uNode = -1; // index into macroNodes
        int vNode = -1; // index into macroNodes
        std::vector<int> triChain;       // ordered active triangles strictly between uTri and vTri
        std::vector<Chord> chords;       // chords[0] exits uTri; chords.back() enters vTri
        std::vector<gsPoint2D<T>> poly;  // medial polyline from uNode.pt to vNode.pt
        bool isCycle = false;
    };

    struct MacroNode {
        gsPoint2D<T> pt;
        T radius = 0;
        std::vector<int> tris;  // 1 CDT triangle (or 2 if contracted J-J pair)
        int degree = 0;         // active macro degree (1=leaf, 2=station, >=3=junction)
    };

    static CurveSeg<T> getCADBoundarySegment(
        const gsPlanarDomain<T>& /*dom*/,
        int degree,
        const TriMeshData<T>& mesh,
        int vStart, int vEnd)
    {
        return extractExactCADBoundaryBetween(mesh, vStart, vEnd, degree, true);
    }

    // ------------------------------------------------------------------------
    // Internal Geometry Helpers
    // ------------------------------------------------------------------------

    // Return (vLeft, vRight) when crossing from CCW triangle tFrom into tTo across their shared edge.
    static bool getSharedChord(const TriMeshData<T>& mesh, int tFrom, int tTo, int& vLeft, int& vRight)
    {
        const auto& A = mesh.tris[tFrom];
        const auto& B = mesh.tris[tTo];
        for (int ka = 0; ka < 3; ++ka) {
            for (int kb = 0; kb < 3; ++kb) {
                if (A.e[ka] == B.e[kb]) {
                    vRight = A.v[ka];
                    vLeft  = A.v[(ka + 1) % 3];
                    return true;
                }
            }
        }
        return false;
    }

    // Collect boundary vertices strictly between vStart and vEnd walking forward on their loop ring.
    // Returns false if vStart and vEnd are not on the same boundary loop.
    static bool collectInteriorLoopVertices(
        const TriMeshData<T>& mesh,
        int vStart, int vEnd,
        std::vector<int>& outInteriorVerts,
        int& outLoopIdx)
    {
        outInteriorVerts.clear();
        outLoopIdx = -1;
        if (vStart < 0 || vEnd < 0 || vStart == vEnd) return false;
        if (static_cast<size_t>(vStart) >= mesh.bndInfo.size() ||
            static_cast<size_t>(vEnd) >= mesh.bndInfo.size())
            return false;
        int l0 = mesh.bndInfo[vStart].loopIdx;
        int l1 = mesh.bndInfo[vEnd].loopIdx;
        if (l0 < 0 || l0 != l1 || l0 >= static_cast<int>(mesh.loopRings.size()))
            return false;

        const auto& verts = mesh.loopRings[l0].vertices;
        const int M = static_cast<int>(verts.size());
        int i0 = -1, i1 = -1;
        for (int i = 0; i < M; ++i) {
            if (verts[i] == vStart) i0 = i;
            if (verts[i] == vEnd)   i1 = i;
        }
        if (i0 < 0 || i1 < 0) return false;

        outLoopIdx = l0;
        int cur = (i0 + 1) % M;
        int guard = 0;
        while (cur != i1 && guard++ < M) {
            outInteriorVerts.push_back(verts[cur]);
            cur = (cur + 1) % M;
        }
        return (cur == i1);
    }

    // Fit a variation-diminishing, spoke-cone-safe cubic Bezier to an interior medial polyline Pa -> Pb.
    // Replaces Runge-prone Lagrange interpolation with convex-hull attenuation + spoke-bisector candidates,
    // picking the curve that maximizes the worst corner-angle margin across both flanking quads!
    static CurveSeg<T> fitMedialCurveConeSafe(
        const std::vector<gsPoint2D<T>>& poly,
        int degree,
        const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
        const CurveSeg<T>& bndL,
        const CurveSeg<T>& bndR,
        const CurveSeg<T>& spLa_Pa,
        const CurveSeg<T>& spPa_Ra,
        const CurveSeg<T>& spPb_Lb,
        const CurveSeg<T>& spRb_Pb)
    {
        const gsPoint2D<T>& A = poly.front();
        const gsPoint2D<T>& B = poly.back();
        const T chordLen = A.dist(B);
        CurveSeg<T> straight = makeStraightSeg(A, B, degree);
        if (chordLen <= T(1e-12))
            return straight;

        auto makeCubic = [&](const gsPoint2D<T>& p1, const gsPoint2D<T>& p2) -> CurveSeg<T> {
            gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
            gsMatrix<T> coefs(4, 2);
            coefs << A.x, A.y, p1.x, p1.y, p2.x, p2.y, B.x, B.y;
            gsBSpline<T> bsp(kv, coefs);
            if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
            return CurveSeg<T>(bsp);
        };

        auto curveInside = [&](const CurveSeg<T>& seg) -> bool {
            for (int k = 1; k < 12; ++k) {
                gsPoint2D<T> q = evalCurvePoint(*seg.bsp, static_cast<T>(k) / T(12));
                if (!isPointStrictlyInsideDomain(q, loopPolys, T(1e-6))) return false;
            }
            return true;
        };

        auto evalCandScore = [&](const CurveSeg<T>& cand) -> T {
            T mnL = 0, mxL = 180, efL = -180;
            T mnR = 0, mxR = 180, efR = -180;
            eval4SegCornerAnglesDeg(cand, spPb_Lb, bndL, spLa_Pa, mnL, mxL, efL);
            eval4SegCornerAnglesDeg(spPa_Ra, bndR, spRb_Pb, cand.reversed(), mnR, mxR, efR);
            return std::min(efL, efR);
        };

        CurveSeg<T> bestSeg = straight;
        T bestScore = curveInside(straight) ? evalCandScore(straight) : -T(1e8);

        // 1. Variation-diminishing convex-hull cubic candidates along arc-length polyline samples (zero Runge overshoot)
        if (poly.size() >= 3) {
            std::vector<T> s(poly.size(), T(0));
            for (size_t i = 1; i < poly.size(); ++i)
                s[i] = s[i - 1] + poly[i - 1].dist(poly[i]);
            const T totS = s.back();
            if (totS > T(1e-12)) {
                auto sampleAt = [&](T frac) -> gsPoint2D<T> {
                    T target = frac * totS;
                    for (size_t i = 1; i < poly.size(); ++i) {
                        if (s[i] >= target) {
                            T ds = s[i] - s[i - 1];
                            T al = (ds > T(1e-14)) ? (target - s[i - 1]) / ds : T(0);
                            return gsPoint2D<T>((T(1) - al) * poly[i - 1].x + al * poly[i].x,
                                                (T(1) - al) * poly[i - 1].y + al * poly[i].y);
                        }
                    }
                    return B;
                };
                gsPoint2D<T> q1 = sampleAt(T(1) / T(3));
                gsPoint2D<T> q2 = sampleAt(T(2) / T(3));
                gsPoint2D<T> c1Lin((T(2) * A.x + B.x) / T(3), (T(2) * A.y + B.y) / T(3));
                gsPoint2D<T> c2Lin((A.x + T(2) * B.x) / T(3), (A.y + T(2) * B.y) / T(3));
                for (T lam : {T(1.0), T(0.70), T(0.45), T(0.25)}) {
                    gsPoint2D<T> p1((T(1) - lam) * c1Lin.x + lam * q1.x, (T(1) - lam) * c1Lin.y + lam * q1.y);
                    gsPoint2D<T> p2((T(1) - lam) * c2Lin.x + lam * q2.x, (T(1) - lam) * c2Lin.y + lam * q2.y);
                    CurveSeg<T> cand = makeCubic(p1, p2);
                    if (!curveInside(cand)) continue;
                    T sc = evalCandScore(cand);
                    if (sc > bestScore) { bestScore = sc; bestSeg = cand; }
                }
            }
        }

        // 2. Spoke-bisector-guided cubic candidates whenever bestScore is still tight (< 26 deg margin, i.e. > 154 deg or < 26 deg)
        if (bestScore < T(26.0)) {
            const T twoPi = T(2) * T(3.14159265358979323846);
            auto normVec = [](const gsPoint2D<T>& v) -> gsPoint2D<T> {
                T n = std::hypot(v.x, v.y);
                return (n > T(1e-14)) ? gsPoint2D<T>(v.x / n, v.y / n) : gsPoint2D<T>(1, 0);
            };
            // At Pa: forward medial cone sweeps CCW from dRa (Pa->vRa) to dLa (Pa->vLa)
            gsPoint2D<T> dRa = normVec(evalCurveDeriv(*spPa_Ra.bsp, T(0)));
            gsPoint2D<T> dLaIn = evalCurveDeriv(*spLa_Pa.bsp, T(1));
            gsPoint2D<T> dLa = normVec(gsPoint2D<T>(-dLaIn.x, -dLaIn.y));
            T aRa = std::atan2(dRa.y, dRa.x);
            T dConeA = std::atan2(dLa.y, dLa.x) - aRa;
            while (dConeA <= T(1e-9)) dConeA += twoPi;
            while (dConeA > twoPi)    dConeA -= twoPi;
            T aBisA = aRa + T(0.5) * dConeA;
            gsPoint2D<T> tBisA(std::cos(aBisA), std::sin(aBisA));

            // At Pb: backward medial cone (Pb->Pa) sweeps CCW from dLb (Pb->vLb) to dRb (Pb->vRb)
            gsPoint2D<T> dLb = normVec(evalCurveDeriv(*spPb_Lb.bsp, T(0)));
            gsPoint2D<T> dRbIn = evalCurveDeriv(*spRb_Pb.bsp, T(1));
            gsPoint2D<T> dRb = normVec(gsPoint2D<T>(-dRbIn.x, -dRbIn.y));
            T aLb = std::atan2(dLb.y, dLb.x);
            T dConeB = std::atan2(dRb.y, dRb.x) - aLb;
            while (dConeB <= T(1e-9)) dConeB += twoPi;
            while (dConeB > twoPi)    dConeB -= twoPi;
            T aBisB = aLb + T(0.5) * dConeB;
            gsPoint2D<T> tBisB(std::cos(aBisB), std::sin(aBisB));

            for (T f : {T(0.28), T(0.20), T(0.13)}) {
                gsPoint2D<T> p1(A.x + f * chordLen * tBisA.x, A.y + f * chordLen * tBisA.y);
                gsPoint2D<T> p2(B.x + f * chordLen * tBisB.x, B.y + f * chordLen * tBisB.y);
                CurveSeg<T> cand = makeCubic(p1, p2);
                if (!curveInside(cand)) continue;
                T sc = evalCandScore(cand);
                if (sc > bestScore) { bestScore = sc; bestSeg = cand; }
            }
        }

        return bestSeg;
    }

    // Build a conforming 5-patch O-grid inside 4 CCW edges e[0..3] with corners P[0..3].
    static bool emitOGrid5Inside(
        const CurveSeg<T> bEdges[4],
        int degree,
        const std::vector<std::vector<gsPoint2D<T>>>& loopPolys,
        std::vector<gsTensorBSpline<2, T>*>& outPatches)
    {
        gsPoint2D<T> P[4];
        for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
        gsPoint2D<T> pAvg(T(0.25) * (P[0].x + P[1].x + P[2].x + P[3].x),
                          T(0.25) * (P[0].y + P[1].y + P[2].y + P[3].y));

        // Sample boundary polygon of this 4-sided cap to compute accurate centroid
        std::vector<gsPoint2D<T>> capPoly;
        capPoly.reserve(32);
        for (int k = 0; k < 4; ++k)
            for (int s = 0; s < 8; ++s)
                capPoly.push_back(evalCurvePoint(*bEdges[k].bsp, static_cast<T>(s) / T(8)));
        gsPoint2D<T> ctrPoly = polygonCentroid(capPoly);

        std::vector<gsPoint2D<T>> cands = {ctrPoly, pAvg};
        const T alphas[4] = {T(0.42), T(0.32), T(0.50), T(0.24)};

        std::vector<gsTensorBSpline<2, T>*> bestFive;
        T bestMaxAng = T(360);

        // Axis vectors u = P[2] - P[0], v = P[3] - P[1] for isotropic rhombus inner-quad candidates
        gsPoint2D<T> d02(P[2].x - P[0].x, P[2].y - P[0].y);
        gsPoint2D<T> d13(P[3].x - P[1].x, P[3].y - P[1].y);
        T len02 = std::hypot(d02.x, d02.y), len13 = std::hypot(d13.x, d13.y);

        for (const auto& ctr : cands) {
            if (!isPointStrictlyInsideDomain(ctr, loopPolys, T(1e-6))) continue;
            for (int mode = 0; mode < 2; ++mode) {
                for (T alpha : alphas) {
                    gsPoint2D<T> Q[4];
                    if (mode == 0) {
                        for (int k = 0; k < 4; ++k) {
                            Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                                ctr.y + alpha * (P[k].y - ctr.y));
                        }
                    } else {
                        if (len02 <= T(1e-10) || len13 <= T(1e-10)) continue;
                        gsPoint2D<T> uDir(d02.x / len02, d02.y / len02);
                        gsPoint2D<T> vDir(-uDir.y, uDir.x);
                        if (vDir.x * d13.x + vDir.y * d13.y < T(0)) { vDir.x = -vDir.x; vDir.y = -vDir.y; }
                        T ru = alpha * T(0.5) * std::min( len02, len13 );
                        T rv = alpha * T(0.5) * std::min( len13, T(1.5) * len02 );
                        Q[0] = gsPoint2D<T>(ctr.x - ru * uDir.x, ctr.y - ru * uDir.y);
                        Q[1] = gsPoint2D<T>(ctr.x - rv * vDir.x, ctr.y - rv * vDir.y);
                        Q[2] = gsPoint2D<T>(ctr.x + ru * uDir.x, ctr.y + ru * uDir.y);
                        Q[3] = gsPoint2D<T>(ctr.x + rv * vDir.x, ctr.y + rv * vDir.y);
                    }
                    bool ok = true;
                    for (int k = 0; k < 4; ++k) {
                        gsPoint2D<T> midQP(T(0.5) * (Q[k].x + P[k].x), T(0.5) * (Q[k].y + P[k].y));
                        gsPoint2D<T> qtrQP(T(0.25) * Q[k].x + T(0.75) * P[k].x, T(0.25) * Q[k].y + T(0.75) * P[k].y);
                        gsPoint2D<T> midQQ(T(0.5) * (Q[k].x + Q[(k + 1) % 4].x),
                                           T(0.5) * (Q[k].y + Q[(k + 1) % 4].y));
                        // Keep Q[k] and spokes strictly inside BOTH the whole domain AND this 4-sided cap itself,
                        // preventing any inner O-grid quad from crossing the cap's interior spoke edges into neighbors.
                        if (!isPointStrictlyInsideDomain(Q[k], loopPolys, T(1e-6)) ||
                            !isPointStrictlyInsideDomain(midQP, loopPolys, T(1e-6)) ||
                            !isPointStrictlyInsideDomain(qtrQP, loopPolys, T(1e-6)) ||
                            !isSegmentInsideDomain(Q[k], Q[(k + 1) % 4], loopPolys) ||
                            !pointInPolyStrict(Q[k], capPoly, T(1e-5)) ||
                            !pointInPolyStrict(midQP, capPoly, T(1e-5)) ||
                            !pointInPolyStrict(qtrQP, capPoly, T(1e-5)) ||
                            !pointInPolyStrict(midQQ, capPoly, T(1e-5))) {
                            ok = false; break;
                        }
                        // Ensure Q stays strictly on the interior side of the cap's two straight spoke edges P[0]->P[1] and P[3]->P[0]
                        if (cross2(P[0], P[1], Q[k]) <= T(1e-7) ||
                            cross2(P[3], P[0], Q[k]) <= T(1e-7)) {
                            ok = false; break;
                        }
                    }
                    if (!ok || !isStrictlyConvexQuad(Q)) continue;

                    CurveSeg<T> qEdge[4], spoke[4];
                    for (int k = 0; k < 4; ++k) {
                        qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], degree);
                        spoke[k] = makeStraightSeg(Q[k], P[k], degree);
                    }

                    std::vector<gsTensorBSpline<2, T>*> trial;
                    trial.reserve(5);
                    trial.push_back(makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], degree));
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        trial.push_back(makeCoonsFromCCWCycle(
                            bEdges[k], spoke[kn].reversed(), qEdge[k].reversed(), spoke[k], degree));
                    }

                    bool trialOk = true;
                    T worstQuality = T(180);
                    for (auto* patch : trial) {
                        if (!patch) { trialOk = false; break; }
                        if (!isPatchDetSingleSign(*patch, 15, T(1e-7)))
                            quickHarmonicPatchRepair(*patch, T(1e-7));
                        if (!isPatchDetSingleSign(*patch, 15, T(1e-7))) { trialOk = false; break; }
                        T angs[4], mnA = 0, mxA = 180;
                        bool degen = false;
                        evalPatchCornerAnglesDeg(*patch, angs, mnA, mxA, degen);
                        if (degen || mnA < T(14) || mxA > T(160)) { trialOk = false; break; }
                        worstQuality = std::min(worstQuality, std::min(mnA, T(180) - mxA));
                    }
                    // Maximize worstQuality (stored inversely in bestMaxAng = 180 - worstQuality)
                    T score = T(180) - worstQuality;
                    if (trialOk && score < bestMaxAng) {
                        for (auto* old : bestFive) delete old;
                        bestFive = std::move(trial);
                        bestMaxAng = score;
                    } else {
                        for (auto* p : trial) delete p;
                    }
                }
            }
        }

        if (bestFive.size() == 5) {
            for (auto* p : bestFive) outPatches.push_back(p);
            return true;
        }
        return false;
    }

    // ------------------------------------------------------------------------
    // Main Entry Point: True 2D Delaunay-Dual MAT Quad Decomposition
    // ------------------------------------------------------------------------
    static bool buildMATMultiPatch(
        const gsPlanarDomain<T>& dom,
        int degree,
        gsMultiPatch<T>& outMp,
        T targetHFrac = T(1.0) / T(22.0),
        T thetaPruneDeg = T(45.0),
        T junctionMergeRatio = T(0.32))
    {
        if (dom.numLoops() == 0) return false;
        outMp.clear();

        const bool dbg = (std::getenv("GISMO_MAT_DEBUG") != nullptr);
        const auto bb = dom.boundingBox();
        const T domDiag = std::hypot(bb(0, 1) - bb(0, 0), bb(1, 1) - bb(1, 0));
        if (domDiag <= T(1e-12)) return false;
        const T targetH = std::max(T(1e-4), domDiag * targetHFrac);
        const T radToDeg = T(180) / T(3.14159265358979323846);

        // --------------------------------------------------------------------
        // Step 1.1: Boundary-Only CDT & T/S/J Classification
        // --------------------------------------------------------------------
        TriMeshData<T> mesh;
        if (!buildConstrainedDomainTriangulation(dom, degree, targetH, true, 0, mesh))
            return false;
        applyDomainLawsonFlips(mesh);

        const int nT = static_cast<int>(mesh.tris.size());
        const int nE = static_cast<int>(mesh.edges.size());
        if (nT < 2) return false;

        std::vector<int> triType(nT, 1);
        std::vector<std::vector<int>> edgeTris(nE);
        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            int nB = 0;
            for (int k = 0; k < 3; ++k) {
                int eid = tr.e[k];
                if (mesh.edges[eid].isBoundary) ++nB;
                else edgeTris[eid].push_back(ti);
            }
            triType[ti] = (nB >= 2) ? 0 : ((nB == 1) ? 1 : 2);
        }

        std::vector<MATNode> triNodes(nT);
        int countT = 0, countS = 0, countJ = 0;
        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            const auto& A = mesh.pts[tr.v[0]];
            const auto& B = mesh.pts[tr.v[1]];
            const auto& C = mesh.pts[tr.v[2]];
            MATNode& node = triNodes[ti];
            node.triIdx = ti;
            node.triType = triType[ti];
            if (node.triType == 0) ++countT;
            else if (node.triType == 1) ++countS;
            else ++countJ;

            // Compute circumcenter; clamp obtuse/out-of-domain triangles safely inside domain
            T D = T(2) * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
            bool useCC = false;
            gsPoint2D<T> cc;
            T cr = T(0);
            if (std::abs(D) > T(1e-12)) {
                T a2 = A.x * A.x + A.y * A.y;
                T b2 = B.x * B.x + B.y * B.y;
                T c2 = C.x * C.x + C.y * C.y;
                cc = gsPoint2D<T>((a2 * (B.y - C.y) + b2 * (C.y - A.y) + c2 * (A.y - B.y)) / D,
                                  (a2 * (C.x - B.x) + b2 * (A.x - C.x) + c2 * (B.x - A.x)) / D);
                cr = cc.dist(A);
                T dAB = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y);
                T dBC = (B.x - C.x) * (B.x - C.x) + (B.y - C.y) * (B.y - C.y);
                T dCA = (C.x - A.x) * (C.x - A.x) + (C.y - A.y) * (C.y - A.y);
                bool isObtuse = (dAB + dBC < dCA - T(1e-9)) || (dBC + dCA < dAB - T(1e-9)) || (dCA + dAB < dBC - T(1e-9));
                if (!isObtuse &&
                    cross2(A, B, cc) > T(1e-8) && cross2(B, C, cc) > T(1e-8) && cross2(C, A, cc) > T(1e-8) &&
                    isPointStrictlyInsideDomain(cc, mesh.loopPolys, T(1e-5))) {
                    useCC = true;
                }
            }

            if (useCC && node.triType != 0) {
                node.pt = cc;
                node.radius = cr;
            } else if (node.triType == 2 || node.triType == 0) {
                // J-triangle or T-triangle: incenter is guaranteed strictly inside the triangle
                // (never on the interior chord midpoint, avoiding 180-deg collinear leaf corners!)
                T la = B.dist(C), lb = C.dist(A), lc = A.dist(B);
                T perim = std::max(T(1e-12), la + lb + lc);
                node.pt = gsPoint2D<T>((la * A.x + lb * B.x + lc * C.x) / perim,
                                       (la * A.y + lb * B.y + lc * C.y) / perim);
                node.radius = std::abs(cross2(A, B, C)) / perim;
            } else {
                // S/T-triangle: chordal midpoint (average of interior edge midpoints)
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
                if (nInt > 0) {
                    node.pt = gsPoint2D<T>(sumM.x / nInt, sumM.y / nInt);
                    node.radius = T(0.5) * maxL;
                } else {
                    node.pt = gsPoint2D<T>((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
                    node.radius = T(0.2) * domDiag;
                }
            }
        }

        for (int ei = 0; ei < nE; ++ei) {
            if (!mesh.edges[ei].isBoundary && edgeTris[ei].size() == 2) {
                int t0 = edgeTris[ei][0], t1 = edgeTris[ei][1];
                triNodes[t0].adj.push_back(t1);
                triNodes[t1].adj.push_back(t0);
            }
        }

        if (dbg) {
            std::cout << "[MAT] CDT: V=" << mesh.pts.size() << " T=" << nT
                      << " (T=" << countT << ", S=" << countS << ", J=" << countJ << ")\n";
            exportCDTToVTK(mesh, triType, "mat_debug_cdt.vtk");
        }

        // --------------------------------------------------------------------
        // Step 1.2: Fast Canonical 0-Hole Special Cases (1-Patch Quad / Smooth O-Grid)
        // --------------------------------------------------------------------
        if (dom.numLoops() == 1) {
            int nOuter = static_cast<int>(dom.outer().size());
            if (nOuter == 4) {
                CurveSeg<T> bsn[4];
                for (int k = 0; k < 4; ++k) bsn[k] = CurveSeg<T>(dom.outer().curve(k), degree);
                gsPoint2D<T> c4[4] = {bsn[0].p0(), bsn[1].p0(), bsn[2].p0(), bsn[3].p0()};
                if (isStrictlyConvexQuad(c4)) {
                    auto patch = makeCoonsFromCCWCycle(bsn[0], bsn[1], bsn[2], bsn[3], degree);
                    if (patch) {
                        T minDet = 0, maxDet = 0, minAbs = 0, en = 0;
                        evalPatchJacobianAndWinslow(*patch, 15, minDet, maxDet, minAbs, en);
                        T angs[4], mnA = 0, mxA = 180;
                        bool degen = false;
                        evalPatchCornerAnglesDeg(*patch, angs, mnA, mxA, degen);
                        if (minDet > T(1e-5) && !degen && mnA > T(10) && mxA < T(155)) {
                            outMp.addPatch(typename gsGeometry<T>::uPtr(patch));
                            outMp.computeTopology(T(1e-4), false, false);
                            if (dbg) std::cout << "[MAT] 1-patch convex quad.\n";
                            return true;
                        }
                        delete patch;
                    }
                }
            }

            int sharpCount = 0;
            for (bool sc : mesh.isSharpCorner) if (sc) ++sharpCount;
            if (sharpCount == 0 && !mesh.loopRings.empty() && mesh.loopRings[0].vertices.size() >= 4) {
                const auto& ringV = mesh.loopRings[0].vertices;
                const int M = static_cast<int>(ringV.size());
                int st[4] = {0, M / 4, M / 2, (3 * M) / 4};
                CurveSeg<T> bArc[4];
                for (int k = 0; k < 4; ++k)
                    bArc[k] = extractExactCADBoundaryBetween(mesh, ringV[st[k]], ringV[st[(k + 1) % 4]], degree, true);
                std::vector<gsTensorBSpline<2, T>*> ogridPatches;
                if (emitOGrid5Inside(bArc, degree, mesh.loopPolys, ogridPatches)) {
                    for (auto* p : ogridPatches)
                        outMp.addPatch(typename gsGeometry<T>::uPtr(p));
                    outMp.computeTopology(T(1e-4), false, false);
                    if (dbg) std::cout << "[MAT] Smooth 5-patch O-grid.\n";
                    return true;
                }
                if (buildSmoothOffsetOGrid(dom, degree, outMp)) {
                    if (dbg) std::cout << "[MAT] Smooth offset O-grid (P=" << outMp.nPatches() << ").\n";
                    return true;
                }
            }
        }

        // --------------------------------------------------------------------
        // Step 1.3: Iterative Theta-SMA Leaf Pruning on the CDT Dual Graph
        // --------------------------------------------------------------------
        auto activeDeg = [&](int ti) -> int {
            if (!triNodes[ti].active) return 0;
            int d = 0;
            for (int nb : triNodes[ti].adj)
                if (triNodes[nb].active) ++d;
            return d;
        };

        auto countActive = [&]() -> int {
            int cnt = 0;
            for (int ti = 0; ti < nT; ++ti) if (triNodes[ti].active) ++cnt;
            return cnt;
        };

        // Helper: inspect the boundary arc covered by a leaf triangle/branch entered via chord (vL, vR)
        auto analyzeCapArc = [&](int vR, int vL, std::vector<int>& sharpConvexOut, int& nInteriorVertsOut) -> bool {
            sharpConvexOut.clear();
            nInteriorVertsOut = 0;
            std::vector<int> intVerts;
            int loopIdx = -1;
            if (!collectInteriorLoopVertices(mesh, vR, vL, intVerts, loopIdx))
                return false;
            nInteriorVertsOut = static_cast<int>(intVerts.size());
            for (int v : intVerts) {
                if (mesh.isSharpCorner[v] && mesh.cornerAngleDeg[v] < T(150.0))
                    sharpConvexOut.push_back(v);
            }
            return true;
        };

        // Pass A: Prune spurious hair branches (root at activeDeg==3 junction, leaf at activeDeg==1)
        bool prunedAny = true;
        while (prunedAny && countActive() > 2) {
            prunedAny = false;
            for (int leaf = 0; leaf < nT; ++leaf) {
                if (activeDeg(leaf) != 1) continue;

                // Walk from leaf back until hitting a node with activeDeg != 2
                std::vector<int> chain = {leaf};
                int prev = leaf;
                int curr = -1;
                for (int nb : triNodes[leaf].adj) if (triNodes[nb].active) { curr = nb; break; }
                while (curr >= 0 && activeDeg(curr) == 2) {
                    chain.push_back(curr);
                    int nxt = -1;
                    for (int nb : triNodes[curr].adj)
                        if (triNodes[nb].active && nb != prev) { nxt = nb; break; }
                    prev = curr;
                    curr = nxt;
                }
                int root = curr;
                if (root < 0 || activeDeg(root) < 3) continue; // only prune hairs attached to junctions

                // Entry chord from root into chain.back()
                int vL = -1, vR = -1;
                if (!getSharedChord(mesh, root, chain.back(), vL, vR)) continue;

                std::vector<int> sharpVerts;
                int nIntV = 0;
                if (!analyzeCapArc(vR, vL, sharpVerts, nIntV)) continue;

                // Never prune a branch whose cap arc would swallow >= 2 sharp corners (keep real arms/teeth!)
                if (sharpVerts.size() >= 2) continue;

                // Compute separation angle theta at root between (root->vL) and (root->vR)
                const auto& rPt = triNodes[root].pt;
                gsPoint2D<T> d1(mesh.pts[vL].x - rPt.x, mesh.pts[vL].y - rPt.y);
                gsPoint2D<T> d2(mesh.pts[vR].x - rPt.x, mesh.pts[vR].y - rPt.y);
                T n1 = std::hypot(d1.x, d1.y), n2 = std::hypot(d2.x, d2.y);
                T sepAngleDeg = T(180);
                if (n1 > T(1e-12) && n2 > T(1e-12)) {
                    T cosA = std::max(T(-1), std::min(T(1), (d1.x * d2.x + d1.y * d2.y) / (n1 * n2)));
                    sepAngleDeg = std::acos(cosA) * radToDeg;
                }

                // Compute polyline length of this hair branch
                T branchLen = triNodes[root].pt.dist(triNodes[chain.back()].pt);
                for (size_t i = 1; i < chain.size(); ++i)
                    branchLen += triNodes[chain[i - 1]].pt.dist(triNodes[chain[i]].pt);

                bool shouldPrune = false;
                if (sharpVerts.empty()) {
                    // Pure smooth-boundary discretization hair
                    if (sepAngleDeg < thetaPruneDeg || branchLen < T(0.85) * triNodes[root].radius)
                        shouldPrune = true;
                } else if (sharpVerts.size() == 1) {
                    // Single sharp convex corner hair: ONLY prune if `root` also has another short
                    // leaf branch (so pruning merges two corner hairs of a blunt rectangular arm end
                    // into a 2-corner cap Subcase B, NEVER swallowing a lone side corner into a branch wall!).
                    T cornerAng = mesh.cornerAngleDeg[sharpVerts[0]];
                    if (cornerAng >= T(65.0) && branchLen < T(1.15) * triNodes[root].radius) {
                        bool hasSiblingLeafHair = false;
                        for (int rnb : triNodes[root].adj) {
                            if (!triNodes[rnb].active || rnb == chain.back()) continue;
                            int p = root, c = rnb;
                            T sibLen = triNodes[root].pt.dist(triNodes[c].pt);
                            int steps = 0;
                            while (c >= 0 && activeDeg(c) == 2 && steps++ < 16) {
                                int nx = -1;
                                for (int x : triNodes[c].adj) if (triNodes[x].active && x != p) { nx = x; break; }
                                if (nx < 0) break;
                                sibLen += triNodes[c].pt.dist(triNodes[nx].pt);
                                p = c; c = nx;
                            }
                            if (c >= 0 && activeDeg(c) == 1 && sibLen < T(1.35) * triNodes[root].radius) {
                                hasSiblingLeafHair = true;
                                break;
                            }
                        }
                        if (hasSiblingLeafHair)
                            shouldPrune = true;
                    }
                }

                if (shouldPrune) {
                    // Verify that pruning keeps active graph non-empty and valid
                    for (int tDel : chain) triNodes[tDel].active = false;
                    prunedAny = true;
                    break;
                }
            }
        }

        // Pass B: Retreat any active leaf whose chord (vL, vR) sits right ON a sharp corner or has 0 interior vertices
        bool retreated = true;
        while (retreated && countActive() > 2) {
            retreated = false;
            for (int leaf = 0; leaf < nT; ++leaf) {
                if (activeDeg(leaf) != 1) continue;
                int parent = -1;
                for (int nb : triNodes[leaf].adj) if (triNodes[nb].active) { parent = nb; break; }
                if (parent < 0) continue;

                int vL = -1, vR = -1;
                if (!getSharedChord(mesh, parent, leaf, vL, vR)) continue;
                std::vector<int> sharpV;
                int nIntV = 0;
                if (!analyzeCapArc(vR, vL, sharpV, nIntV)) continue;

                // If the leaf arc has 0 interior vertices (degenerate 1-edge cap),
                // both endpoints are sharp corners on a straight end wall, OR
                // it is a shallow 1-2 vertex smooth cap (avoiding 10:1 flat smooth sliver caps), retreat leaf into parent
                bool bothSharpEnd = (mesh.isSharpCorner[vR] && mesh.isSharpCorner[vL] && sharpV.empty());
                bool shallowSmoothCap = (sharpV.empty() && nIntV <= 2);
                if (nIntV == 0 || bothSharpEnd || shallowSmoothCap) {
                    // Only retreat if parent still has another neighbor and won't accumulate > 2 sharp corners
                    if (activeDeg(parent) == 2) {
                        int grandParent = -1;
                        for (int nb : triNodes[parent].adj)
                            if (triNodes[nb].active && nb != leaf) { grandParent = nb; break; }
                        if (grandParent >= 0) {
                            int gvL = -1, gvR = -1;
                            if (getSharedChord(mesh, grandParent, parent, gvL, gvR)) {
                                std::vector<int> gSharp;
                                int gInt = 0;
                                int maxAllowedSharp = shallowSmoothCap ? 0 : 2;
                                if (analyzeCapArc(gvR, gvL, gSharp, gInt) &&
                                    static_cast<int>(gSharp.size()) <= maxAllowedSharp) {
                                    triNodes[leaf].active = false;
                                    retreated = true;
                                    break;
                                }
                            }
                        }
                    } else if (activeDeg(parent) == 3 && nIntV == 0) {
                        triNodes[leaf].active = false;
                        retreated = true;
                        break;
                    }
                }
            }
        }

        // Pass C (Retreat-to-Stub): For short sharp-corner hair branches (1 corner or short 2-corner stub)
        // attached to a junction `root`, deactivate the redundant sleeve triangles chain[0..size-2] and keep
        // ONLY the entry triangle chain.back() as a 0-sleeve degree-1 stub leaf. Its end-cap then attaches
        // directly to `root` (1 clean kite quad per corner instead of 3-5 quads!) while preserving exact sharp corners!
        {
            for (int leaf = 0; leaf < nT; ++leaf) {
                if (activeDeg(leaf) != 1) continue;
                std::vector<int> chain = {leaf};
                int prev = leaf, curr = -1;
                for (int nb : triNodes[leaf].adj) if (triNodes[nb].active) { curr = nb; break; }
                while (curr >= 0 && activeDeg(curr) == 2) {
                    chain.push_back(curr);
                    int nxt = -1;
                    for (int nb : triNodes[curr].adj)
                        if (triNodes[nb].active && nb != prev) { nxt = nb; break; }
                    prev = curr;
                    curr = nxt;
                }
                int root = curr;
                if (root < 0 || activeDeg(root) < 3 || chain.size() < 2) continue;

                int vL = -1, vR = -1;
                if (!getSharedChord(mesh, root, chain.back(), vL, vR)) continue;
                std::vector<int> sharpVerts;
                int nIntV = 0;
                if (!analyzeCapArc(vR, vL, sharpVerts, nIntV)) continue;
                if (sharpVerts.size() != 1 && sharpVerts.size() != 2) continue;

                T branchLen = triNodes[root].pt.dist(triNodes[chain.back()].pt);
                for (size_t i = 1; i < chain.size(); ++i)
                    branchLen += triNodes[chain[i - 1]].pt.dist(triNodes[chain[i]].pt);
                T maxLenMult = (sharpVerts.size() == 1) ? T(1.95) : T(1.10);
                if (branchLen >= maxLenMult * triNodes[root].radius) continue;

                const auto& rPt = triNodes[root].pt;
                const auto& pL = mesh.pts[vL];
                const auto& pR = mesh.pts[vR];
                if (cross2(pR, pL, rPt) <= T(1e-7)) continue;
                gsPoint2D<T> d1(pL.x - rPt.x, pL.y - rPt.y);
                gsPoint2D<T> d2(pR.x - rPt.x, pR.y - rPt.y);
                T n1 = std::hypot(d1.x, d1.y), n2 = std::hypot(d2.x, d2.y);
                if (n1 <= T(1e-12) || n2 <= T(1e-12)) continue;
                // Avoid retreating a very long thin needle tip or creating > 152 deg angles at vL/vR
                if (sharpVerts.size() == 1) {
                    const auto& pC = mesh.pts[sharpVerts[0]];
                    T tipDist = rPt.dist(pC);
                    T chordW  = pL.dist(pR);
                    if (chordW > T(1e-12) && tipDist > T(2.8) * chordW) continue;
                    // Check straight kite corner angles at pR and pL (rPt -> pR -> pC and pC -> pL -> rPt)
                    auto cornerDeg = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C) -> T {
                        T ux = A.x - B.x, uy = A.y - B.y, vx = C.x - B.x, vy = C.y - B.y;
                        T lu = std::hypot(ux, uy), lv = std::hypot(vx, vy);
                        if (lu <= T(1e-12) || lv <= T(1e-12)) return T(180);
                        T cs = std::max(T(-1), std::min(T(1), (ux * vx + uy * vy) / (lu * lv)));
                        return std::acos(cs) * radToDeg;
                    };
                    T aR = cornerDeg(rPt, pR, pC);
                    T aL = cornerDeg(pC, pL, rPt);
                    if (aR < T(22.0) || aR > T(152.0) || aL < T(22.0) || aL > T(152.0)) continue;
                }
                T cosA = std::max(T(-1), std::min(T(1), (d1.x * d2.x + d1.y * d2.y) / (n1 * n2)));
                T sepDeg = std::acos(cosA) * radToDeg;
                if (sepDeg < T(24.0) || sepDeg > T(148.0)) continue;
                gsPoint2D<T> midL(T(0.5) * (rPt.x + pL.x), T(0.5) * (rPt.y + pL.y));
                gsPoint2D<T> midR(T(0.5) * (rPt.x + pR.x), T(0.5) * (rPt.y + pR.y));
                if (!isPointStrictlyInsideDomain(midL, mesh.loopPolys, T(1e-6)) ||
                    !isPointStrictlyInsideDomain(midR, mesh.loopPolys, T(1e-6)))
                    continue;

                // Deactivate sleeve triangles, leaving chain.back() as a 0-sleeve degree-1 stub leaf on root
                for (size_t i = 0; i + 1 < chain.size(); ++i)
                    triNodes[chain[i]].active = false;
            }
        }

        // --------------------------------------------------------------------
        // Step 1.4: Compress Active Skeleton into MacroNodes & MacroBranches
        // --------------------------------------------------------------------
        std::vector<int> triToMacro(nT, -1);
        std::vector<MacroNode> macroNodes;

        auto registerMacroNode = [&](int ti, int deg) -> int {
            if (triToMacro[ti] >= 0) return triToMacro[ti];
            int id = static_cast<int>(macroNodes.size());
            MacroNode mn;
            mn.pt = triNodes[ti].pt;
            mn.radius = triNodes[ti].radius;
            mn.tris = {ti};
            mn.degree = deg;
            macroNodes.push_back(mn);
            triToMacro[ti] = id;
            return id;
        };

        for (int ti = 0; ti < nT; ++ti) {
            if (!triNodes[ti].active) continue;
            int d = activeDeg(ti);
            if (d == 1 || d >= 3)
                registerMacroNode(ti, d);
        }

        // Handle pure closed cycles (active components where every active triangle has activeDeg == 2,
        // e.g. single-hole plate / wavy ring annulus) by placing 4 well-aligned stations per cycle.
        {
            std::vector<bool> vis(nT, false);
            for (int start = 0; start < nT; ++start) {
                if (!triNodes[start].active || vis[start] || triToMacro[start] >= 0) continue;
                // Check if this connected component has any macro node
                std::vector<int> comp;
                std::queue<int> q;
                q.push(start);
                vis[start] = true;
                bool hasMacro = false;
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    comp.push_back(u);
                    if (triToMacro[u] >= 0) hasMacro = true;
                    for (int nb : triNodes[u].adj) {
                        if (triNodes[nb].active && !vis[nb]) {
                            vis[nb] = true;
                            q.push(nb);
                        }
                    }
                }
                if (hasMacro || comp.size() < 4) continue;

                // Order the pure cycle sequentially
                std::vector<int> cyc;
                cyc.reserve(comp.size());
                int cur = comp[0], prev = -1;
                for (size_t step = 0; step < comp.size(); ++step) {
                    cyc.push_back(cur);
                    int nxt = -1;
                    for (int nb : triNodes[cur].adj) {
                        if (triNodes[nb].active && nb != prev) { nxt = nb; break; }
                    }
                    if (nxt < 0 || nxt == comp[0]) break;
                    prev = cur;
                    cur = nxt;
                }
                int C = static_cast<int>(cyc.size());
                if (C < 4) continue;

                // Pick 4 stations roughly spaced by C/4, preferring chords touching sharp corners
                int bestOff = 0;
                int bestSharpScore = -1;
                for (int off = 0; off < std::max(1, C / 4); ++off) {
                    int sc = 0;
                    for (int s = 0; s < 4; ++s) {
                        int idx = (off + (s * C) / 4) % C;
                        int nxt = (idx + 1) % C;
                        int vL = -1, vR = -1;
                        if (getSharedChord(mesh, cyc[idx], cyc[nxt], vL, vR)) {
                            if (mesh.isSharpCorner[vL]) ++sc;
                            if (mesh.isSharpCorner[vR]) ++sc;
                        }
                    }
                    if (sc > bestSharpScore) { bestSharpScore = sc; bestOff = off; }
                }
                for (int s = 0; s < 4; ++s) {
                    int ti = cyc[(bestOff + (s * C) / 4) % C];
                    registerMacroNode(ti, 2);
                }
            }
        }

        // If skeleton collapsed completely (0 or 1 active macro nodes), bail cleanly
        if (macroNodes.size() < 2) {
            if (dbg) std::cout << "[MAT] Collapsed to < 2 macro nodes; skipping.\n";
            return false;
        }

        // Contract short J-J edges (0 sleeve triangles between adjacent junction MacroNodes, including
        // multi-triangle clusters up to valence 6) whenever a star-kernel point sees all outer half-edges.
        if (junctionMergeRatio > T(0)) {
            bool didMerge = true;
            while (didMerge) {
                didMerge = false;
                for (size_t i = 0; i < macroNodes.size() && !didMerge; ++i) {
                    if (macroNodes[i].tris.empty() || macroNodes[i].degree < 3 || macroNodes[i].degree >= 6)
                        continue;
                    for (int tA : macroNodes[i].tris) {
                        for (int nb : triNodes[tA].adj) {
                            if (!triNodes[nb].active) continue;
                            int j = triToMacro[nb];
                            if (j < 0 || j == static_cast<int>(i)) continue;
                            if (macroNodes[j].tris.empty() || macroNodes[j].degree < 3) continue;
                            int vL = -1, vR = -1;
                            if (!getSharedChord(mesh, tA, nb, vL, vR)) continue;

                            T d = macroNodes[i].pt.dist(macroNodes[j].pt);
                            T rMin = std::min(macroNodes[i].radius, macroNodes[j].radius);
                            T chordLen = mesh.pts[vL].dist(mesh.pts[vR]);
                            if (d >= std::max(junctionMergeRatio * rMin, T(0.75) * chordLen)) continue;

                            // Build combined triangle set and verify new active degree <= 6
                            std::set<int> combTris(macroNodes[i].tris.begin(), macroNodes[i].tris.end());
                            combTris.insert(macroNodes[j].tris.begin(), macroNodes[j].tris.end());
                            int newDeg = 0;
                            for (int ct : combTris) {
                                for (int cnb : triNodes[ct].adj) {
                                    if (triNodes[cnb].active && !combTris.count(cnb))
                                        ++newDeg;
                                }
                            }
                            if (newDeg < 3 || newDeg > 6) continue;

                            // Collect all outer CCW directed half-edges (a -> b) of the merged cluster
                            std::vector<std::pair<int, int>> extHalfEdges;
                            gsPoint2D<T> vSum(0, 0);
                            int vCnt = 0;
                            for (int ct : combTris) {
                                const auto& tr = mesh.tris[ct];
                                for (int k = 0; k < 3; ++k) {
                                    int eid = tr.e[k];
                                    bool isInternal = false;
                                    if (!mesh.edges[eid].isBoundary) {
                                        for (int ot : edgeTris[eid]) {
                                            if (ot != ct && combTris.count(ot)) { isInternal = true; break; }
                                        }
                                    }
                                    if (!isInternal) {
                                        int va = tr.v[k], vb = tr.v[(k + 1) % 3];
                                        extHalfEdges.push_back({va, vb});
                                        vSum.x += mesh.pts[va].x;
                                        vSum.y += mesh.pts[va].y;
                                        ++vCnt;
                                    }
                                }
                            }
                            if (extHalfEdges.size() < 4 || vCnt == 0) continue;

                            gsPoint2D<T> avgN(T(0.5) * (macroNodes[i].pt.x + macroNodes[j].pt.x),
                                              T(0.5) * (macroNodes[i].pt.y + macroNodes[j].pt.y));
                            gsPoint2D<T> midC(T(0.5) * (mesh.pts[vL].x + mesh.pts[vR].x),
                                              T(0.5) * (mesh.pts[vL].y + mesh.pts[vR].y));
                            gsPoint2D<T> ctrQ(vSum.x / T(vCnt), vSum.y / T(vCnt));
                            gsPoint2D<T> midPt = avgN;
                            bool foundKernel = false;
                            for (const auto& cand : {avgN, ctrQ, midC,
                                                     gsPoint2D<T>(T(0.5) * (avgN.x + ctrQ.x), T(0.5) * (avgN.y + ctrQ.y)),
                                                     gsPoint2D<T>(T(0.5) * (avgN.x + midC.x), T(0.5) * (avgN.y + midC.y))}) {
                                if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-5))) continue;
                                bool allPos = true;
                                for (const auto& he : extHalfEdges) {
                                    if (cross2(mesh.pts[he.first], mesh.pts[he.second], cand) <= T(1e-7)) {
                                        allPos = false; break;
                                    }
                                }
                                if (allPos) {
                                    midPt = cand;
                                    foundKernel = true;
                                    break;
                                }
                            }
                            if (!foundKernel) continue;

                            // Contract j into i
                            macroNodes[i].pt = midPt;
                            macroNodes[i].radius = T(0.5) * (macroNodes[i].radius + macroNodes[j].radius);
                            for (int tj : macroNodes[j].tris) {
                                macroNodes[i].tris.push_back(tj);
                                triToMacro[tj] = static_cast<int>(i);
                            }
                            macroNodes[i].degree = newDeg;
                            macroNodes[j].tris.clear();
                            macroNodes[j].degree = 0;
                            didMerge = true;
                            break;
                        }
                        if (didMerge) break;
                    }
                }
            }
        }

        // Trace all MacroBranches between active MacroNodes
        std::vector<MacroBranch> branches;
        std::set<std::pair<int, int>> usedDualEdges;

        for (size_t mi = 0; mi < macroNodes.size(); ++mi) {
            if (macroNodes[mi].tris.empty()) continue;
            for (int tStart : macroNodes[mi].tris) {
                for (int nb : triNodes[tStart].adj) {
                    if (!triNodes[nb].active) continue;
                    if (triToMacro[nb] == static_cast<int>(mi)) continue; // internal edge of contracted J-J
                    int lo = std::min(tStart, nb), hi = std::max(tStart, nb);
                    if (usedDualEdges.count({lo, hi})) continue;
                    usedDualEdges.insert({lo, hi});

                    MacroBranch br;
                    br.uNode = static_cast<int>(mi);
                    br.poly.push_back(macroNodes[mi].pt);

                    int vL0 = -1, vR0 = -1;
                    if (!getSharedChord(mesh, tStart, nb, vL0, vR0)) continue;
                    br.chords.push_back({vL0, vR0});

                    int prev = tStart;
                    int curr = nb;
                    while (triToMacro[curr] < 0) {
                        br.triChain.push_back(curr);
                        br.poly.push_back(triNodes[curr].pt);
                        int nxt = -1;
                        for (int cnb : triNodes[curr].adj) {
                            if (triNodes[cnb].active && cnb != prev) { nxt = cnb; break; }
                        }
                        if (nxt < 0) break;
                        usedDualEdges.insert({std::min(curr, nxt), std::max(curr, nxt)});
                        int vL = -1, vR = -1;
                        if (getSharedChord(mesh, curr, nxt, vL, vR))
                            br.chords.push_back({vL, vR});
                        prev = curr;
                        curr = nxt;
                    }
                    br.vNode = triToMacro[curr];
                    if (br.vNode < 0) continue;
                    br.poly.push_back(macroNodes[br.vNode].pt);
                    branches.push_back(std::move(br));
                }
            }
        }

        if (branches.empty()) return false;

        // Smooth interior medial polylines of every branch (4 constrained Laplacian sweeps)
        for (auto& br : branches) {
            if (br.poly.size() < 3) continue;
            int np = static_cast<int>(br.poly.size());
            for (int iter = 0; iter < 4; ++iter) {
                for (int i = 1; i < np - 1; ++i) {
                    gsPoint2D<T> sm(T(0.5) * br.poly[i].x + T(0.25) * (br.poly[i - 1].x + br.poly[i + 1].x),
                                    T(0.5) * br.poly[i].y + T(0.25) * (br.poly[i - 1].y + br.poly[i + 1].y));
                    if (isPointStrictlyInsideDomain(sm, mesh.loopPolys, T(1e-5)))
                        br.poly[i] = sm;
                }
            }
        }

        // Position every active degree-1 leaf MacroNode along the medial direction behind its entry chord
        // midpoint M so PL never sits on M (no 180-deg cap corner), never enters leafTri (no >180-deg dart),
        // and reaches ~110..150 deg even on wide, short sleeves where |pAdj - M| << Lc.
        for (size_t mi = 0; mi < macroNodes.size(); ++mi) {
            auto& mn = macroNodes[mi];
            if (mn.tris.empty() || mn.degree != 1) continue;
            int leafTri = mn.tris[0];
            int parentTri = -1;
            for (int nb : triNodes[leafTri].adj)
                if (triNodes[nb].active) { parentTri = nb; break; }
            if (parentTri < 0) continue;
            if (triToMacro[parentTri] >= 0) {
                mn.pt = macroNodes[triToMacro[parentTri]].pt;
                for (auto& br : branches) {
                    if (br.uNode == static_cast<int>(mi) && !br.poly.empty()) br.poly.front() = mn.pt;
                    if (br.vNode == static_cast<int>(mi) && !br.poly.empty()) br.poly.back()  = mn.pt;
                }
                continue;
            }
            int vL = -1, vR = -1;
            if (!getSharedChord(mesh, parentTri, leafTri, vL, vR)) continue;
            const gsPoint2D<T>& pL = mesh.pts[vL];
            const gsPoint2D<T>& pR = mesh.pts[vR];
            gsPoint2D<T> M(T(0.5) * (pL.x + pR.x), T(0.5) * (pL.y + pR.y));
            T Lc = pL.dist(pR);
            if (Lc <= T(1e-12)) continue;

            // Locate the incident branch and adjacent smoothed medial points pAdj, pNext
            int incBrIdx = -1;
            bool leafAtFront = true;
            gsPoint2D<T> pAdj = triNodes[parentTri].pt;
            gsPoint2D<T> pNext = pAdj;
            bool hasNext = false;
            for (size_t bi = 0; bi < branches.size(); ++bi) {
                const auto& br = branches[bi];
                if (br.uNode == static_cast<int>(mi) && br.poly.size() >= 2) {
                    incBrIdx = static_cast<int>(bi);
                    leafAtFront = true;
                    pAdj = br.poly[1];
                    if (br.poly.size() >= 3) { pNext = br.poly[2]; hasNext = true; }
                    break;
                }
                if (br.vNode == static_cast<int>(mi) && br.poly.size() >= 2) {
                    incBrIdx = static_cast<int>(bi);
                    leafAtFront = false;
                    pAdj = br.poly[br.poly.size() - 2];
                    if (br.poly.size() >= 3) { pNext = br.poly[br.poly.size() - 3]; hasNext = true; }
                    break;
                }
            }

            gsPoint2D<T> bestPL(T(0.60) * M.x + T(0.40) * pAdj.x, T(0.60) * M.y + T(0.40) * pAdj.y);
            T bestScore = std::numeric_limits<T>::max();
            T bestAngDeg = T(180);
            bool usedDeepNext = false;

            auto evalLeafCand = [&](const gsPoint2D<T>& cand, bool isDeep) {
                if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-6))) return;
                if (cross2(pR, pL, cand) <= T(1e-7)) return;
                gsPoint2D<T> midL(T(0.5) * (cand.x + pL.x), T(0.5) * (cand.y + pL.y));
                gsPoint2D<T> midR(T(0.5) * (cand.x + pR.x), T(0.5) * (cand.y + pR.y));
                if (!isPointStrictlyInsideDomain(midL, mesh.loopPolys, T(1e-6)) ||
                    !isPointStrictlyInsideDomain(midR, mesh.loopPolys, T(1e-6)))
                    return;
                gsPoint2D<T> dL(pL.x - cand.x, pL.y - cand.y);
                gsPoint2D<T> dR(pR.x - cand.x, pR.y - cand.y);
                T nL = std::hypot(dL.x, dL.y), nR = std::hypot(dR.x, dR.y);
                if (nL <= T(1e-12) || nR <= T(1e-12)) return;
                T cosA = std::max(T(-1), std::min(T(1), (dL.x * dR.x + dL.y * dR.y) / (nL * nR)));
                T angDeg = std::acos(cosA) * radToDeg;
                // Smooth objective: target ~120 deg, strongly penalize > 152 deg, never drop all candidates
                T sc = std::abs(angDeg - T(120)) + (angDeg > T(152) ? T(12) * (angDeg - T(152)) : T(0))
                     + (angDeg < T(75) ? T(4) * (T(75) - angDeg) : T(0))
                     + (isDeep ? T(3) : T(0));
                if (sc < bestScore) {
                    bestScore = sc;
                    bestAngDeg = angDeg;
                    bestPL = cand;
                    usedDeepNext = isDeep;
                }
            };

            // Pass 1: candidates along M -> pAdj (including up to 0.80 when pAdj is close to M)
            for (T al : {T(0.32), T(0.45), T(0.58), T(0.70), T(0.80), T(0.24), T(0.16)}) {
                gsPoint2D<T> cand((T(1) - al) * M.x + al * pAdj.x, (T(1) - al) * M.y + al * pAdj.y);
                evalLeafCand(cand, false);
            }

            // Pass 2: if sleeve is wide & short (bestAng > 150 deg), search deeper toward pNext or inward chord normal
            if (bestAngDeg > T(150) && hasNext) {
                for (T al : {T(0.25), T(0.36), T(0.48), T(0.60), T(0.72)}) {
                    gsPoint2D<T> cand((T(1) - al) * M.x + al * pNext.x, (T(1) - al) * M.y + al * pNext.y);
                    gsPoint2D<T> midN(T(0.5) * (cand.x + pNext.x), T(0.5) * (cand.y + pNext.y));
                    if (isPointStrictlyInsideDomain(midN, mesh.loopPolys, T(1e-6)))
                        evalLeafCand(cand, true);
                }
                // Inward unit normal to chord pR -> pL (points into sleeve since cross2(pR, pL, cand) > 0)
                gsPoint2D<T> nIn(-(pL.y - pR.y) / Lc, (pL.x - pR.x) / Lc);
                for (T hf : {T(0.18), T(0.24), T(0.30)}) {
                    gsPoint2D<T> cand(M.x + hf * Lc * nIn.x, M.y + hf * Lc * nIn.y);
                    // Only accept if cand stays strictly before pNext along M -> pNext
                    T dotRem = (pNext.x - cand.x) * (pNext.x - M.x) + (pNext.y - cand.y) * (pNext.y - M.y);
                    gsPoint2D<T> midN(T(0.5) * (cand.x + pNext.x), T(0.5) * (cand.y + pNext.y));
                    if (dotRem > T(1e-8) && isPointStrictlyInsideDomain(midN, mesh.loopPolys, T(1e-6)))
                        evalLeafCand(cand, true);
                }
            }

            mn.pt = bestPL;
            for (auto& br : branches) {
                if (br.uNode == static_cast<int>(mi) && !br.poly.empty()) br.poly.front() = mn.pt;
                if (br.vNode == static_cast<int>(mi) && !br.poly.empty()) br.poly.back()  = mn.pt;
            }
            // If PL moved deeper toward pNext, re-interpolate pAdj between PL and pNext so polyline stays monotone
            if (usedDeepNext && incBrIdx >= 0 && branches[incBrIdx].poly.size() >= 3) {
                auto& bp = branches[incBrIdx].poly;
                size_t adjIdx = leafAtFront ? 1 : (bp.size() - 2);
                gsPoint2D<T> newAdj(T(0.5) * (bestPL.x + pNext.x), T(0.5) * (bestPL.y + pNext.y));
                if (isPointStrictlyInsideDomain(newAdj, mesh.loopPolys, T(1e-6)))
                    bp[adjIdx] = newAdj;
            }
        }

        // Export skeleton debug VTK if requested
        if (dbg) {
            std::vector<std::vector<gsPoint2D<T>>> dbgPolys;
            std::vector<T> dbgPRad;
            std::vector<gsPoint2D<T>> dbgJunc;
            std::vector<T> dbgJRad;
            for (const auto& br : branches) {
                dbgPolys.push_back(br.poly);
                dbgPRad.push_back(T(0.5) * (macroNodes[br.uNode].radius + macroNodes[br.vNode].radius));
            }
            int nJunc = 0, nLeaf = 0;
            for (const auto& mn : macroNodes) {
                if (mn.tris.empty()) continue;
                if (mn.degree >= 3) { ++nJunc; dbgJunc.push_back(mn.pt); dbgJRad.push_back(mn.radius); }
                else if (mn.degree == 1) ++nLeaf;
            }
            std::cout << "[MAT] Pruned skeleton: branches=" << branches.size()
                      << " junctions=" << nJunc << " leaves=" << nLeaf << "\n";
            exportSkeletonToVTK(dbgPolys, dbgPRad, dbgJunc, dbgJRad, "mat_debug_skeleton.vtk");
        }

        // --------------------------------------------------------------------
        // Step 1.5: Pre-compute Branch Stations (with Reflex Fan Splitting)
        //           and Boundary-Cone-Safe Curved MAT Spokes
        // --------------------------------------------------------------------
        struct Station {
            gsPoint2D<T> pt;
            int vL, vR;
            int polyIdx;
        };

        auto isSpokeInteriorClear = [&](const gsPoint2D<T>& medPt, const gsPoint2D<T>& bndPt) -> bool {
            for (T t : {T(0.25), T(0.50), T(0.75)}) {
                gsPoint2D<T> q((T(1) - t) * medPt.x + t * bndPt.x, (T(1) - t) * medPt.y + t * bndPt.y);
                if (!isPointStrictlyInsideDomain(q, mesh.loopPolys, T(1e-6))) return false;
            }
            return true;
        };

        std::vector<std::vector<Station>> allBranchStations(branches.size());
        for (size_t bi = 0; bi < branches.size(); ++bi) {
            const auto& br = branches[bi];
            const int nChords = static_cast<int>(br.chords.size());
            if (nChords == 0) continue;
            // If a leaf MacroNode is directly adjacent to a junction MacroNode (0 sleeve triangles),
            // its end-cap is emitted directly from the junction MacroNode (avoiding a valence-2 Case 1 diamond!).
            if (br.triChain.empty() && (macroNodes[br.uNode].degree == 1 || macroNodes[br.vNode].degree == 1))
                continue;

            auto& stations = allBranchStations[bi];
            stations.push_back({macroNodes[br.uNode].pt, br.chords.front().vL, br.chords.front().vR, 0});

            if (nChords >= 3) {
                const T avgRad = std::max(T(1e-6), T(0.5) * (macroNodes[br.uNode].radius + macroNodes[br.vNode].radius));
                T accumLen = T(0);
                T accumTurn = T(0);
                gsPoint2D<T> prevDir(0, 0);
                bool hasDir = false;

                for (int ci = 1; ci < nChords - 1; ++ci) {
                    int pi = std::min(ci, static_cast<int>(br.poly.size()) - 2);
                    gsPoint2D<T> curPt = br.poly[pi];
                    gsPoint2D<T> step(curPt.x - br.poly[pi - 1].x, curPt.y - br.poly[pi - 1].y);
                    T stepL = std::hypot(step.x, step.y);
                    accumLen += stepL;
                    if (stepL > T(1e-12)) {
                        gsPoint2D<T> curDir(step.x / stepL, step.y / stepL);
                        if (hasDir) {
                            T cosT = std::max(T(-1), std::min(T(1), prevDir.x * curDir.x + prevDir.y * curDir.y));
                            accumTurn += std::acos(cosT) * radToDeg;
                        }
                        prevDir = curDir;
                        hasDir = true;
                    }

                    int vLi = br.chords[ci].vL;
                    int vRi = br.chords[ci].vR;
                    const auto& prevSt = stations.back();
                    const auto& endCh  = br.chords.back();

                    // 1. Force a station when chord hits a sharp CAD corner (including FIRST contact with a true reflex corner > 205 deg)
                    bool hitSharpL = mesh.isSharpCorner[vLi] && vLi != prevSt.vL &&
                                     (vLi != endCh.vL || (mesh.cornerAngleDeg[vLi] > T(205.0) && vRi != prevSt.vR && vRi != endCh.vR));
                    bool hitSharpR = mesh.isSharpCorner[vRi] && vRi != prevSt.vR &&
                                     (vRi != endCh.vR || (mesh.cornerAngleDeg[vRi] > T(205.0) && vLi != prevSt.vL && vLi != endCh.vL));
                    bool hitSharpCorner = hitSharpL || hitSharpR;

                    // 2. Reflex-corner fan splitting: ONLY when one wall is pinned at a true reflex sharp corner
                    //    (cornerAngleDeg > 198 deg, NOT smooth 180 deg walls!) and the spoke around R has swept > 68 deg.
                    bool reflexFanSplit = false;
                    auto checkReflexFan = [&](int vPin, int vOpp, int prevOpp, int endOpp) -> bool {
                        if (vOpp == prevOpp || vOpp == endOpp) return false;
                        if (!mesh.isSharpCorner[vPin] || mesh.cornerAngleDeg[vPin] <= T(198.0)) return false;
                        const auto& R = mesh.pts[vPin];
                        gsPoint2D<T> u0(prevSt.pt.x - R.x, prevSt.pt.y - R.y);
                        gsPoint2D<T> u1(curPt.x - R.x, curPt.y - R.y);
                        T l0 = std::hypot(u0.x, u0.y), l1 = std::hypot(u1.x, u1.y);
                        if (l0 <= T(1e-12) || l1 <= T(1e-12)) return false;
                        T c = std::max(T(-1), std::min(T(1), (u0.x * u1.x + u0.y * u1.y) / (l0 * l1)));
                        T sweepDeg = std::acos(c) * radToDeg;
                        T threshDeg = (mesh.cornerAngleDeg[vPin] > T(240.0)) ? T(65.0) : T(78.0);
                        return sweepDeg > threshDeg;
                    };
                    if (vLi == prevSt.vL)
                        reflexFanSplit = checkReflexFan(vLi, vRi, prevSt.vR, endCh.vR);
                    else if (vRi == prevSt.vR)
                        reflexFanSplit = checkReflexFan(vRi, vLi, prevSt.vL, endCh.vL);

                    bool wantSplit = hitSharpCorner || reflexFanSplit ||
                                     (accumTurn > T(65.0)) || (accumLen > T(3.6) * avgRad);
                    int minRem = (hitSharpCorner || reflexFanSplit) ? 1 : 2;
                    if (wantSplit && (nChords - 1 - ci) >= minRem) {
                        bool bothAdvance = (vLi != prevSt.vL && vRi != prevSt.vR) &&
                                           (vLi != endCh.vL  && vRi != endCh.vR);
                        bool oneAdvances = (vLi != prevSt.vL || vRi != prevSt.vR) &&
                                           (vLi != endCh.vL  || vRi != endCh.vR);
                        if (bothAdvance || ((hitSharpCorner || reflexFanSplit) && oneAdvances)) {
                            gsPoint2D<T> candPt = curPt;
                            if (!isPointStrictlyInsideDomain(candPt, mesh.loopPolys, T(1e-5)) ||
                                !isSpokeInteriorClear(candPt, mesh.pts[vLi]) ||
                                !isSpokeInteriorClear(candPt, mesh.pts[vRi]))
                            {
                                candPt = gsPoint2D<T>(T(0.5) * (mesh.pts[vLi].x + mesh.pts[vRi].x),
                                                      T(0.5) * (mesh.pts[vLi].y + mesh.pts[vRi].y));
                            }
                            if (isPointStrictlyInsideDomain(candPt, mesh.loopPolys, T(1e-5)) &&
                                isSpokeInteriorClear(candPt, mesh.pts[vLi]) &&
                                isSpokeInteriorClear(candPt, mesh.pts[vRi]))
                            {
                                stations.push_back({candPt, vLi, vRi, pi});
                                accumLen = T(0);
                                accumTurn = T(0);
                            }
                        }
                    }
                }
            }

            stations.push_back({macroNodes[br.vNode].pt, br.chords.back().vL, br.chords.back().vR,
                                static_cast<int>(br.poly.size()) - 1});
        }

        // Pre-register every mesh-vertex spoke (vBnd, P) so all spokes incident to each boundary vertex
        // are coordinated into a single non-crossing, boundary-cone-safe curved Bezier registry (exact C^0!).
        std::map<int, std::vector<gsPoint2D<T>>> bndSpokeTargets;
        auto registerBndSpoke = [&](int vBnd, const gsPoint2D<T>& medPt) {
            if (vBnd < 0 || static_cast<size_t>(vBnd) >= mesh.pts.size()) return;
            auto& vec = bndSpokeTargets[vBnd];
            for (const auto& q : vec)
                if (q.dist(medPt) < T(1e-8)) return;
            vec.push_back(medPt);
        };

        for (const auto& stations : allBranchStations) {
            for (const auto& st : stations) {
                registerBndSpoke(st.vL, st.pt);
                registerBndSpoke(st.vR, st.pt);
            }
        }
        for (size_t mi = 0; mi < macroNodes.size(); ++mi) {
            const auto& mn = macroNodes[mi];
            if (mn.tris.empty() || mn.degree != 1) continue;
            int leafTri = mn.tris[0];
            int parentTri = -1;
            for (int nb : triNodes[leafTri].adj)
                if (triNodes[nb].active) { parentTri = nb; break; }
            if (parentTri < 0) continue;
            int vL = -1, vR = -1;
            if (!getSharedChord(mesh, parentTri, leafTri, vL, vR)) continue;
            registerBndSpoke(vR, mn.pt);
            registerBndSpoke(vL, mn.pt);
            std::vector<int> sharpV;
            int nIntV = 0;
            if (analyzeCapArc(vR, vL, sharpV, nIntV) && sharpV.size() >= 3)
                registerBndSpoke(sharpV[sharpV.size() / 2], mn.pt);
        }

        // Determine loop orientation defensively (outer CCW => domainOnLeft, hole CW => domainOnLeft)
        std::vector<bool> loopDomainOnLeft(mesh.loopPolys.size(), true);
        for (size_t l = 0; l < mesh.loopPolys.size(); ++l) {
            const auto& poly = mesh.loopPolys[l];
            T area2 = T(0);
            for (size_t k = 0; k < poly.size(); ++k) {
                size_t kn = (k + 1) % poly.size();
                area2 += poly[k].x * poly[kn].y - poly[kn].x * poly[k].y;
            }
            loopDomainOnLeft[l] = (l == 0) ? (area2 >= T(0)) : (area2 <= T(0));
        }

        // Build curved cubic Bezier spoke B -> P with exact departure tangent tB at B and exact arrival tangent uBP at P.
        // Optional neighbor rays (prevP, nextP) prevent any curved spoke from bulging across an adjacent spoke at B!
        auto buildCurvedSpokeBndToMed = [&](const gsPoint2D<T>& B, const gsPoint2D<T>& P,
                                            const gsPoint2D<T>& tB, T maxT1Abs,
                                            const gsPoint2D<T>* prevP, const gsPoint2D<T>* nextP) -> CurveSeg<T> {
            T L = B.dist(P);
            if (L <= T(1e-12)) return makeStraightSeg(B, P, degree);
            gsPoint2D<T> uBP((P.x - B.x) / L, (P.y - B.y) / L);
            const T mults[5] = {T(1.0), T(0.72), T(0.48), T(0.30), T(0.16)};
            for (T m : mults) {
                T t1 = std::min(T(0.24) * L, m * maxT1Abs);
                T t2 = std::min(T(0.28) * L, T(1.30) * t1);
                gsPoint2D<T> c1(B.x + t1 * tB.x, B.y + t1 * tB.y);
                gsPoint2D<T> c2(P.x - t2 * uBP.x, P.y - t2 * uBP.y);
                gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
                gsMatrix<T> cf(4, 2);
                cf << B.x, B.y, c1.x, c1.y, c2.x, c2.y, P.x, P.y;
                gsBSpline<T> bsp(kv, cf);
                if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
                bool insideOk = true;
                for (int ks = 1; ks < 12; ++ks) {
                    gsPoint2D<T> q = evalCurvePoint(bsp, static_cast<T>(ks) / T(12));
                    if (!isPointStrictlyInsideDomain(q, mesh.loopPolys, T(1e-6))) { insideOk = false; break; }
                    // Stay strictly between adjacent spokes (prevP is CW, nextP is CCW from P around B)
                    if (prevP && cross2(B, *prevP, P) > T(1e-8) && cross2(B, *prevP, q) <= T(1e-9)) {
                        insideOk = false; break;
                    }
                    if (nextP && cross2(B, P, *nextP) > T(1e-8) && cross2(B, q, *nextP) <= T(1e-9)) {
                        insideOk = false; break;
                    }
                }
                if (insideOk) return CurveSeg<T>(bsp);
            }
            return makeStraightSeg(B, P, degree);
        };

        // Compute coordinated arrival tangents for every boundary vertex and cache CurveSeg(vBnd -> P)
        const T piVal = T(3.14159265358979323846);
        const T twoPi = T(2) * piVal;
        struct CachedSpoke { gsPoint2D<T> P; CurveSeg<T> seg; };
        std::map<int, std::vector<CachedSpoke>> bndSpokeCache;

        for (const auto& kvPair : bndSpokeTargets) {
            int vBnd = kvPair.first;
            const auto& ptsVec = kvPair.second;
            const int K = static_cast<int>(ptsVec.size());
            if (K == 0) continue;
            const gsPoint2D<T>& B = mesh.pts[vBnd];

            // Interior CCW cone [aStart, aStart + dCCW] at vBnd
            const int lIdx = (static_cast<size_t>(vBnd) < mesh.bndInfo.size()) ? mesh.bndInfo[vBnd].loopIdx : 0;
            const bool domOnLeft = (lIdx >= 0 && static_cast<size_t>(lIdx) < loopDomainOnLeft.size())
                                 ? loopDomainOnLeft[lIdx] : true;
            const gsPoint2D<T>& tOut = mesh.bndTanOut[vBnd];
            const gsPoint2D<T>& tIn  = mesh.bndTanIn[vBnd];
            gsPoint2D<T> tBack(-tIn.x, -tIn.y);
            const gsPoint2D<T>& ray0 = domOnLeft ? tOut : tBack;
            const gsPoint2D<T>& ray1 = domOnLeft ? tBack : tOut;
            T aStart = std::atan2(ray0.y, ray0.x);
            T dCCW   = std::atan2(ray1.y, ray1.x) - aStart;
            while (dCCW <= T(1e-9)) dCCW += twoPi;
            while (dCCW > twoPi)    dCCW -= twoPi;

            struct SpokeItem { gsPoint2D<T> P; T rawAng; T len; };
            std::vector<SpokeItem> items;
            items.reserve(K);
            const T midExt = T(0.5) * (dCCW + twoPi);
            T minSpokeLen = std::numeric_limits<T>::max();
            for (const auto& P : ptsVec) {
                T len = B.dist(P);
                minSpokeLen = std::min(minSpokeLen, len);
                T ang = std::atan2(P.y - B.y, P.x - B.x) - aStart;
                while (ang < T(0)) ang += twoPi;
                while (ang >= twoPi) ang -= twoPi;
                // Map exterior overshoot smoothly to negative (near ray0) vs > dCCW (near ray1)
                if (ang > midExt) ang -= twoPi;
                items.push_back({P, ang, len});
            }
            std::sort(items.begin(), items.end(),
                      [](const SpokeItem& a, const SpokeItem& b) { return a.rawAng < b.rawAng; });

            const T avgGap   = dCCW / T(K + 1);
            const T wallMin  = std::min(T(0.46), T(0.42) * avgGap);  // >= ~26.4 deg from CAD walls
            const T spokeMin = std::min(T(0.18), T(0.35) * avgGap);  // >= ~10.3 deg between adjacent spokes
            const T gMax     = std::max(T(2.54), avgGap + T(0.03));  // <= ~145.5 deg whenever feasible

            bool alreadySafe = true;
            for (int j = 0; j < K; ++j) {
                T reqMin = (j == 0) ? wallMin : spokeMin;
                T gapPrev = (j == 0) ? items[0].rawAng : (items[j].rawAng - items[j - 1].rawAng);
                if (gapPrev < reqMin || gapPrev > gMax) { alreadySafe = false; break; }
            }
            if (alreadySafe && ((dCCW - items.back().rawAng) < wallMin || (dCCW - items.back().rawAng) > gMax))
                alreadySafe = false;

            std::vector<T> targetAng(K);
            for (int j = 0; j < K; ++j)
                targetAng[j] = std::max(wallMin, std::min(dCCW - wallMin, items[j].rawAng));

            if (!alreadySafe) {
                for (int sweep = 0; sweep < 40; ++sweep) {
                    for (int j = 0; j < K; ++j) {
                        T lo = (j == 0) ? wallMin : (targetAng[j - 1] + spokeMin);
                        T hi = (j == 0) ? gMax    : (targetAng[j - 1] + gMax);
                        targetAng[j] = std::max(lo, std::min(hi, targetAng[j]));
                    }
                    for (int j = K - 1; j >= 0; --j) {
                        T hi = (j == K - 1) ? (dCCW - wallMin) : (targetAng[j + 1] - spokeMin);
                        T lo = (j == K - 1) ? (dCCW - gMax)    : (targetAng[j + 1] - gMax);
                        targetAng[j] = std::max(lo, std::min(hi, targetAng[j]));
                    }
                }
                bool projOk = true;
                for (int j = 0; j < K; ++j) {
                    T reqMin = (j == 0) ? wallMin : spokeMin;
                    T gp = (j == 0) ? targetAng[0] : (targetAng[j] - targetAng[j - 1]);
                    if (gp < reqMin - T(1e-3) || gp > gMax + T(1e-3)) projOk = false;
                }
                if ((dCCW - targetAng.back()) < wallMin - T(1e-3) || (dCCW - targetAng.back()) > gMax + T(1e-3))
                    projOk = false;
                if (!projOk) {
                    for (int j = 0; j < K; ++j)
                        targetAng[j] = avgGap * T(j + 1);
                }
            }

            auto& cacheVec = bndSpokeCache[vBnd];
            cacheVec.reserve(K);
            for (int j = 0; j < K; ++j) {
                T diffAng = std::abs(targetAng[j] - items[j].rawAng);
                CurveSeg<T> seg;
                if (diffAng <= T(0.02)) {
                    seg = makeStraightSeg(B, items[j].P, degree);
                } else {
                    T absA = aStart + targetAng[j];
                    gsPoint2D<T> tB(std::cos(absA), std::sin(absA));
                    // Cap handle length by 0.35 * minSpokeLen so a long curved spoke never sweeps past a short neighbor spoke
                    T maxT1Abs = (K == 1) ? (T(0.22) * items[j].len)
                                          : std::min(T(0.16) * items[j].len, T(0.35) * minSpokeLen);
                    const gsPoint2D<T>* prevP = (j > 0) ? &items[j - 1].P : nullptr;
                    const gsPoint2D<T>* nextP = (j + 1 < K) ? &items[j + 1].P : nullptr;
                    seg = buildCurvedSpokeBndToMed(B, items[j].P, tB, maxT1Abs, prevP, nextP);
                }
                cacheVec.push_back({items[j].P, seg});
            }
        }

        auto getSpokeBndToMed = [&](int vBnd, const gsPoint2D<T>& medPt) -> CurveSeg<T> {
            auto it = bndSpokeCache.find(vBnd);
            if (it != bndSpokeCache.end()) {
                for (const auto& cs : it->second)
                    if (cs.P.dist(medPt) < T(1e-7)) return cs.seg;
            }
            return makeStraightSeg(mesh.pts[vBnd], medPt, degree);
        };
        auto getSpokeMedToBnd = [&](const gsPoint2D<T>& medPt, int vBnd) -> CurveSeg<T> {
            return getSpokeBndToMed(vBnd, medPt).reversed();
        };

        // Helper for spokes landing on a smooth interior parameter uAtM of an oriented CAD boundary curve bndSeg
        auto makeMidCurveSpokeMedToBnd = [&](const gsPoint2D<T>& medPt, const CurveSeg<T>& bndSeg, T uAtM) -> CurveSeg<T> {
            gsPoint2D<T> M = evalCurvePoint(*bndSeg.bsp, uAtM);
            gsPoint2D<T> der = evalCurveDeriv(*bndSeg.bsp, uAtM);
            T nd = std::hypot(der.x, der.y);
            if (nd <= T(1e-12)) return makeStraightSeg(medPt, M, degree);
            T aStart = std::atan2(der.y / nd, der.x / nd);
            T ang = std::atan2(medPt.y - M.y, medPt.x - M.x) - aStart;
            while (ang < T(0)) ang += twoPi;
            while (ang >= twoPi) ang -= twoPi;
            if (ang > T(1.5) * piVal) ang -= twoPi;
            const T margin = T(0.48); // ~27.5 deg from boundary tangent on both sides
            if (ang >= margin && ang <= piVal - margin)
                return makeStraightSeg(medPt, M, degree);
            T targetA = std::max(margin, std::min(piVal - margin, ang));
            gsPoint2D<T> tB(std::cos(aStart + targetA), std::sin(aStart + targetA));
            return buildCurvedSpokeBndToMed(M, medPt, tB, T(0.20) * M.dist(medPt), nullptr, nullptr).reversed();
        };

        // --------------------------------------------------------------------
        // Step 1.6: True MAT Quad Emission (Branches + Bend Stations + Caps)
        // --------------------------------------------------------------------
        std::vector<gsTensorBSpline<2, T>*> quadPatches;
        std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> debugSpokes;
        bool emitFailed = false;

        auto pushQuad = [&](const CurveSeg<T>& e0, const CurveSeg<T>& e1,
                            const CurveSeg<T>& e2, const CurveSeg<T>& e3) -> bool {
            auto* patch = makeCoonsFromCCWCycle(e0, e1, e2, e3, degree);
            if (!patch) return false;
            if (evalPatchMinDetDense(*patch, 49) <= T(1e-7)) {
                quickHarmonicPatchRepair(*patch, T(1e-7));
                // If a stiff 4-CP transverse direction opposes a highly curved multi-span wall (>=12 CPs),
                // insert one midpoint knot along the stiff direction and re-repair interior CPs.
                if (evalPatchMinDetDense(*patch, 49) <= T(1e-7)) {
                    int nu = patch->basis().component(0).size();
                    int nv = patch->basis().component(1).size();
                    if (nu == degree + 1 && nv >= 10) {
                        patch->insertKnot(T(0.5), 0, 1);
                        quickHarmonicPatchRepair(*patch, T(1e-7));
                    } else if (nv == degree + 1 && nu >= 10) {
                        patch->insertKnot(T(0.5), 1, 1);
                        quickHarmonicPatchRepair(*patch, T(1e-7));
                    }
                }
            }
            quadPatches.push_back(patch);
            return true;
        };

        // Emit conforming quads for ONE sub-segment of a medial branch:
        // from medial point P_a (chord vL_a, vR_a) to P_b (chord vL_b, vR_b).
        auto emitBranchSegment = [&](
            const gsPoint2D<T>& Pa, int vLa, int vRa,
            const gsPoint2D<T>& Pb, int vLb, int vRb,
            const std::vector<gsPoint2D<T>>& segPoly) -> bool
        {
            const gsPoint2D<T>& pLa = mesh.pts[vLa];
            const gsPoint2D<T>& pRa = mesh.pts[vRa];
            const gsPoint2D<T>& pLb = mesh.pts[vLb];
            const gsPoint2D<T>& pRb = mesh.pts[vRb];

            debugSpokes.push_back({Pa, pLa});
            debugSpokes.push_back({Pa, pRa});
            debugSpokes.push_back({Pb, pLb});
            debugSpokes.push_back({Pb, pRb});

            // Case 1: Both chords identical (0-sleeve edge between two unmerged J-triangles)
            // -> 1 conforming diamond quad: Pa -> vR -> Pb -> vL
            if (vLa == vLb && vRa == vRb) {
                CurveSeg<T> e0 = getSpokeMedToBnd(Pa, vRa);
                CurveSeg<T> e1 = getSpokeBndToMed(vRa, Pb);
                CurveSeg<T> e2 = getSpokeMedToBnd(Pb, vLa);
                CurveSeg<T> e3 = getSpokeBndToMed(vLa, Pa);
                return pushQuad(e0, e1, e2, e3);
            }

            // Case 2: Left wall pinned at reflex/pivot vertex R = pLa (vLa == vLb), Right wall advances (vRa -> vRb)
            // -> 3-Quad Medial Y-Split: keeps medial spine Pa -> Pm -> Pb strictly separating left and right spokes
            //    at Pa and Pb (never putting pLa -> Pa -> pRa in the same quad!), with Pm optimized so
            //    all 3 quads (Left: Pa->Pm->Pb->R, Right1: Pa->vRa->MR->Pm, Right2: Pm->MR->vRb->Pb) have valid angles!
            if (vLa == vLb && vRa != vRb) {
                CurveSeg<T> bndR = extractExactCADBoundaryBetween(mesh, vRa, vRb, degree, true);
                gsPoint2D<T> Mab(T(0.5) * (Pa.x + Pb.x), T(0.5) * (Pa.y + Pb.y));
                T uSplit = std::max(T(0.25), std::min(T(0.75), bndR.closestParam(Mab)));
                auto spR = bndR.splitAtParam(uSplit);
                gsPoint2D<T> MR = spR.first.p1();

                CurveSeg<T> spPa_Ra = getSpokeMedToBnd(Pa, vRa);
                CurveSeg<T> spRb_Pb = getSpokeBndToMed(vRb, Pb);
                CurveSeg<T> spPb_R  = getSpokeMedToBnd(Pb, vLa);
                CurveSeg<T> spR_Pa  = getSpokeBndToMed(vLa, Pa);

                gsPoint2D<T> Pm(Mab.x + T(0.28) * (Mab.x - pLa.x), Mab.y + T(0.28) * (Mab.y - pLa.y));
                T bestEff = -T(1e9);
                std::vector<gsPoint2D<T>> pmCands;
                // Median-ray extension from R=pLa through Mab: mathematically guaranteed inside cone (R->Pa, R->Pb)!
                for (T gam : {T(0.25), T(0.38), T(0.52), T(0.68), T(0.85), T(0.16), T(0.10)})
                    pmCands.emplace_back(Mab.x + gam * (Mab.x - pLa.x), Mab.y + gam * (Mab.y - pLa.y));
                for (T beta : {T(0.28), T(0.18), T(0.38), T(0.52), T(0.64)}) {
                    for (T sShift : {T(0.0), T(-0.16), T(0.16)}) {
                        gsPoint2D<T> base((T(1) - std::abs(sShift)) * Mab.x + std::abs(sShift) * (sShift > 0 ? Pb.x : Pa.x),
                                          (T(1) - std::abs(sShift)) * Mab.y + std::abs(sShift) * (sShift > 0 ? Pb.y : Pa.y));
                        pmCands.emplace_back((T(1) - beta) * base.x + beta * MR.x, (T(1) - beta) * base.y + beta * MR.y);
                    }
                }
                for (const auto& cand : pmCands) {
                    gsPoint2D<T> mMR(T(0.5) * (cand.x + MR.x), T(0.5) * (cand.y + MR.y));
                    if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-6)) ||
                        !isPointStrictlyInsideDomain(mMR, mesh.loopPolys, T(1e-6)))
                        continue;
                    if (cross2(Pa, cand, Pb) <= T(1e-8) ||
                        cross2(cand, Pb, pLa) <= T(1e-8) ||
                        cross2(pLa, Pa, cand) <= T(1e-8))
                        continue;
                    CurveSeg<T> m1 = makeStraightSeg(Pa, cand, degree);
                    CurveSeg<T> m2 = makeStraightSeg(cand, Pb, degree);
                    CurveSeg<T> sM = makeMidCurveSpokeMedToBnd(cand, bndR, uSplit);
                    T mn0, mx0, ef0, mn1, mx1, ef1, mn2, mx2, ef2;
                    eval4SegCornerAnglesDeg(m1, m2, spPb_R, spR_Pa, mn0, mx0, ef0);
                    eval4SegCornerAnglesDeg(spPa_Ra, spR.first, sM.reversed(), m1.reversed(), mn1, mx1, ef1);
                    eval4SegCornerAnglesDeg(sM, spR.second, spRb_Pb, m2.reversed(), mn2, mx2, ef2);
                    T eff = std::min(ef0, std::min(ef1, ef2));
                    if (eff > bestEff) { bestEff = eff; Pm = cand; }
                }
                debugSpokes.push_back({Pm, MR});
                CurveSeg<T> spPm_MR = makeMidCurveSpokeMedToBnd(Pm, bndR, uSplit);
                std::vector<gsPoint2D<T>> poly1 = {Pa, Pm}, poly2 = {Pm, Pb};
                CurveSeg<T> med1 = fitMedialCurveConeSafe(
                    poly1, degree, mesh.loopPolys, spPb_R, spR.first, spR_Pa, spPa_Ra, makeStraightSeg(Pm, Pb, degree), spPm_MR.reversed());
                CurveSeg<T> med2 = fitMedialCurveConeSafe(
                    poly2, degree, mesh.loopPolys, spR_Pa, spR.second, med1, spPm_MR, spPb_R, spRb_Pb);

                // Left Quad (CCW):   Pa -> Pm -> Pb -> R
                // Right Quad 1 (CCW): Pa -> vRa -> MR -> Pm
                // Right Quad 2 (CCW): Pm -> MR -> vRb -> Pb
                return pushQuad(med1, med2, spPb_R, spR_Pa) &&
                       pushQuad(spPa_Ra, spR.first, spPm_MR.reversed(), med1.reversed()) &&
                       pushQuad(spPm_MR, spR.second, spRb_Pb, med2.reversed());
            }

            // Case 3: Right wall pinned at reflex/pivot vertex R = pRa (vRa == vRb), Left wall advances (vLb -> vLa)
            // -> Symmetric 3-Quad Medial Y-Split with Pm along median ray from R=pRa through Mab / toward ML!
            if (vRa == vRb && vLa != vLb) {
                CurveSeg<T> bndL = extractExactCADBoundaryBetween(mesh, vLb, vLa, degree, true);
                gsPoint2D<T> Mab(T(0.5) * (Pa.x + Pb.x), T(0.5) * (Pa.y + Pb.y));
                T uSplit = std::max(T(0.25), std::min(T(0.75), bndL.closestParam(Mab)));
                auto spL = bndL.splitAtParam(uSplit);
                gsPoint2D<T> ML = spL.first.p1();

                CurveSeg<T> spPa_R  = getSpokeMedToBnd(Pa, vRa);
                CurveSeg<T> spR_Pb  = getSpokeBndToMed(vRa, Pb);
                CurveSeg<T> spPb_Lb = getSpokeMedToBnd(Pb, vLb);
                CurveSeg<T> spLa_Pa = getSpokeBndToMed(vLa, Pa);

                gsPoint2D<T> Pm(Mab.x + T(0.28) * (Mab.x - pRa.x), Mab.y + T(0.28) * (Mab.y - pRa.y));
                T bestEff = -T(1e9);
                std::vector<gsPoint2D<T>> pmCands;
                for (T gam : {T(0.25), T(0.38), T(0.52), T(0.68), T(0.85), T(0.16), T(0.10)})
                    pmCands.emplace_back(Mab.x + gam * (Mab.x - pRa.x), Mab.y + gam * (Mab.y - pRa.y));
                for (T beta : {T(0.28), T(0.18), T(0.38), T(0.52), T(0.64)}) {
                    for (T sShift : {T(0.0), T(-0.16), T(0.16)}) {
                        gsPoint2D<T> base((T(1) - std::abs(sShift)) * Mab.x + std::abs(sShift) * (sShift > 0 ? Pb.x : Pa.x),
                                          (T(1) - std::abs(sShift)) * Mab.y + std::abs(sShift) * (sShift > 0 ? Pb.y : Pa.y));
                        pmCands.emplace_back((T(1) - beta) * base.x + beta * ML.x, (T(1) - beta) * base.y + beta * ML.y);
                    }
                }
                for (const auto& cand : pmCands) {
                    gsPoint2D<T> mML(T(0.5) * (cand.x + ML.x), T(0.5) * (cand.y + ML.y));
                    if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-6)) ||
                        !isPointStrictlyInsideDomain(mML, mesh.loopPolys, T(1e-6)))
                        continue;
                    if (cross2(Pb, cand, Pa) <= T(1e-8) ||
                        cross2(cand, Pa, pRa) <= T(1e-8) ||
                        cross2(pRa, Pb, cand) <= T(1e-8))
                        continue;
                    CurveSeg<T> m1 = makeStraightSeg(Pa, cand, degree);
                    CurveSeg<T> m2 = makeStraightSeg(cand, Pb, degree);
                    CurveSeg<T> sM = makeMidCurveSpokeMedToBnd(cand, bndL, uSplit);
                    T mn0, mx0, ef0, mn1, mx1, ef1, mn2, mx2, ef2;
                    eval4SegCornerAnglesDeg(spPa_R, spR_Pb, m2.reversed(), m1.reversed(), mn0, mx0, ef0);
                    eval4SegCornerAnglesDeg(m1, sM, spL.second, spLa_Pa, mn1, mx1, ef1);
                    eval4SegCornerAnglesDeg(m2, spPb_Lb, spL.first, sM.reversed(), mn2, mx2, ef2);
                    T eff = std::min(ef0, std::min(ef1, ef2));
                    if (eff > bestEff) { bestEff = eff; Pm = cand; }
                }
                debugSpokes.push_back({Pm, ML});
                CurveSeg<T> spPm_ML = makeMidCurveSpokeMedToBnd(Pm, bndL, uSplit);
                // Note: for Case 3, med1 runs Pa -> Pm with Left Quad 1 on left and Right Quad on right
                std::vector<gsPoint2D<T>> poly1 = {Pa, Pm}, poly2 = {Pm, Pb};
                CurveSeg<T> med1 = fitMedialCurveConeSafe(
                    poly1, degree, mesh.loopPolys, spL.second, spR_Pb, spLa_Pa, spPa_R, spPm_ML, makeStraightSeg(Pb, Pm, degree));
                CurveSeg<T> med2 = fitMedialCurveConeSafe(
                    poly2, degree, mesh.loopPolys, spL.first, spPa_R, spPm_ML.reversed(), med1.reversed(), spPb_Lb, spR_Pb);

                // Right Quad (CCW):  Pa -> R -> Pb -> Pm
                // Left Quad 1 (CCW): Pa -> Pm -> ML -> vLa
                // Left Quad 2 (CCW): Pm -> Pb -> vLb -> ML
                return pushQuad(spPa_R, spR_Pb, med2.reversed(), med1.reversed()) &&
                       pushQuad(med1, spPm_ML, spL.second, spLa_Pa) &&
                       pushQuad(med2, spPb_Lb, spL.first, spPm_ML.reversed());
            }

            // Case 4: Standard branch segment (both walls advance: vLa != vLb and vRa != vRb)
            // -> 2 conforming quads flanking the variation-diminishing, spoke-cone-safe medial curve Pa -> Pb
            CurveSeg<T> bndL = extractExactCADBoundaryBetween(mesh, vLb, vLa, degree, true);
            CurveSeg<T> bndR = extractExactCADBoundaryBetween(mesh, vRa, vRb, degree, true);

            CurveSeg<T> spPb_Lb = getSpokeMedToBnd(Pb, vLb);
            CurveSeg<T> spLa_Pa = getSpokeBndToMed(vLa, Pa);
            CurveSeg<T> spPa_Ra = getSpokeMedToBnd(Pa, vRa);
            CurveSeg<T> spRb_Pb = getSpokeBndToMed(vRb, Pb);

            CurveSeg<T> med = fitMedialCurveConeSafe(
                segPoly, degree, mesh.loopPolys, bndL, bndR, spLa_Pa, spPa_Ra, spPb_Lb, spRb_Pb);

            // Left quad (CCW):  Pa -> Pb -> vLb -> vLa
            // Right quad (CCW): Pa -> vRa -> vRb -> Pb
            return pushQuad(med, spPb_Lb, bndL, spLa_Pa) &&
                   pushQuad(spPa_Ra, bndR, spRb_Pb, med.reversed());
        };

        // Emit all MacroBranches using pre-computed stations & boundary-cone-safe curved spokes
        for (size_t bi = 0; bi < branches.size() && !emitFailed; ++bi) {
            const auto& br = branches[bi];
            const auto& stations = allBranchStations[bi];
            for (size_t si = 0; si + 1 < stations.size(); ++si) {
                const auto& stA = stations[si];
                const auto& stB = stations[si + 1];
                std::vector<gsPoint2D<T>> subPoly;
                subPoly.push_back(stA.pt);
                for (int pi = stA.polyIdx + 1; pi < stB.polyIdx && pi < static_cast<int>(br.poly.size()); ++pi)
                    subPoly.push_back(br.poly[pi]);
                subPoly.push_back(stB.pt);

                if (!emitBranchSegment(stA.pt, stA.vL, stA.vR, stB.pt, stB.vL, stB.vR, subPoly)) {
                    emitFailed = true;
                    break;
                }
            }
        }

        // Emit end-caps for every active degree-1 leaf MacroNode
        for (size_t mi = 0; mi < macroNodes.size() && !emitFailed; ++mi) {
            const auto& mn = macroNodes[mi];
            if (mn.tris.empty() || mn.degree != 1) continue;

            int leafTri = mn.tris[0];
            int parentTri = -1;
            for (int nb : triNodes[leafTri].adj)
                if (triNodes[nb].active) { parentTri = nb; break; }
            if (parentTri < 0) continue;

            // Crossing from parentTri INTO leafTri gives (vL, vR) facing into the leaf
            int vL = -1, vR = -1;
            if (!getSharedChord(mesh, parentTri, leafTri, vL, vR)) { emitFailed = true; break; }

            std::vector<int> sharpV;
            int nIntV = 0;
            if (!analyzeCapArc(vR, vL, sharpV, nIntV)) { emitFailed = true; break; }

            const gsPoint2D<T>& PL = mn.pt;

            if (sharpV.size() == 1) {
                // Subcase A: Single sharp tip C -> 1 kite quad (PL -> vR -> C -> vL)
                int vC = sharpV[0];
                CurveSeg<T> sp0  = getSpokeMedToBnd(PL, vR);
                CurveSeg<T> bnd1 = extractExactCADBoundaryBetween(mesh, vR, vC, degree, true);
                CurveSeg<T> bnd2 = extractExactCADBoundaryBetween(mesh, vC, vL, degree, true);
                CurveSeg<T> sp3  = getSpokeBndToMed(vL, PL);
                if (!pushQuad(sp0, bnd1, bnd2, sp3)) emitFailed = true;
            } else if (sharpV.size() == 2) {
                // Subcase B: Blunt rectangular end (2 corners C1, C2) -> 2 corner quads
                // (reducePolychordStrips collapses the centerline across M_end into a 1-wide arm!)
                int vC1 = sharpV[0];
                int vC2 = sharpV[1];
                CurveSeg<T> endWall = extractExactCADBoundaryBetween(mesh, vC1, vC2, degree, true);
                T uMid = arcLengthSplitParam(*endWall.bsp, T(0.5));
                auto spW = endWall.splitAtParam(uMid);
                gsPoint2D<T> Mend = spW.first.p1();
                debugSpokes.push_back({PL, Mend});

                CurveSeg<T> spPL_R   = getSpokeMedToBnd(PL, vR);
                CurveSeg<T> bndR_C1  = extractExactCADBoundaryBetween(mesh, vR, vC1, degree, true);
                CurveSeg<T> spM_PL   = makeMidCurveSpokeMedToBnd(PL, endWall, uMid).reversed();
                CurveSeg<T> bndC2_L  = extractExactCADBoundaryBetween(mesh, vC2, vL, degree, true);
                CurveSeg<T> spL_PL   = getSpokeBndToMed(vL, PL);

                // Quad 1 (CCW): PL -> vR -> C1 -> Mend
                // Quad 2 (CCW): PL -> Mend -> C2 -> vL
                if (!pushQuad(spPL_R, bndR_C1, spW.first, spM_PL) ||
                    !pushQuad(spM_PL.reversed(), spW.second, bndC2_L, spL_PL))
                    emitFailed = true;
            } else if (sharpV.empty()) {
                // Subcase C: Smooth rounded end-cap (0 sharp corners) -> 5-patch O-grid cap (0 180-deg corners!)
                CurveSeg<T> fullArc = extractExactCADBoundaryBetween(mesh, vR, vL, degree, true);
                auto spArc = fullArc.splitArc(T(0.5));
                CurveSeg<T> capEdges[4] = {
                    getSpokeMedToBnd(PL, vR),
                    spArc.first,
                    spArc.second,
                    getSpokeBndToMed(vL, PL)
                };
                if (!emitOGrid5Inside(capEdges, degree, mesh.loopPolys, quadPatches)) {
                    // Fallback: single quad if O-grid didn't fit (e.g. very shallow gentle arc)
                    if (!pushQuad(capEdges[0], capEdges[1], capEdges[2], capEdges[3]))
                        emitFailed = true;
                }
            } else {
                // Subcase D (rare: >= 3 sharp corners): emit 2 quads via middle corner
                int midIdx = static_cast<int>(sharpV.size()) / 2;
                int vMid = sharpV[midIdx];
                CurveSeg<T> spPL_R  = getSpokeMedToBnd(PL, vR);
                CurveSeg<T> spM_PL  = getSpokeBndToMed(vMid, PL);
                CurveSeg<T> spL_PL  = getSpokeBndToMed(vL, PL);
                int vC1 = sharpV[0];
                int vC2 = sharpV.back();
                if (!pushQuad(spPL_R,
                              extractExactCADBoundaryBetween(mesh, vR, vC1, degree, true),
                              extractExactCADBoundaryBetween(mesh, vC1, vMid, degree, true),
                              spM_PL) ||
                    !pushQuad(spM_PL.reversed(),
                              extractExactCADBoundaryBetween(mesh, vMid, vC2, degree, true),
                              extractExactCADBoundaryBetween(mesh, vC2, vL, degree, true),
                              spL_PL))
                    emitFailed = true;
            }
        }

        if (dbg) {
            exportSpokesToVTK(debugSpokes, "mat_debug_spokes.vtk");
            std::cout << "[MAT] Emitted " << quadPatches.size()
                      << " MAT quads (emitFailed=" << emitFailed << ").\n";
        }

        if (emitFailed || quadPatches.empty()) {
            for (auto* p : quadPatches) delete p;
            return false;
        }

        for (auto* p : quadPatches)
            outMp.addPatch(typename gsGeometry<T>::uPtr(p));

        outMp.computeTopology(T(1e-4), false, false);
        const bool topoOk = (2 * outMp.nInterfaces() + outMp.nBoundary() == 4 * outMp.nPatches());
        if (!topoOk) {
            if (dbg) std::cerr << "[MAT] Non-conforming topology (P=" << outMp.nPatches()
                               << ", IF=" << outMp.nInterfaces()
                               << ", BND=" << outMp.nBoundary() << ").\n";
            outMp.clear();
            return false;
        }
        return true;
    }
};

} // namespace internal_splitter
} // namespace gismo
