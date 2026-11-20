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

    // Fit a straight segment or smooth cubic Bezier to an interior medial polyline, staying inside domain.
    static CurveSeg<T> fitMedialCurve(
        const std::vector<gsPoint2D<T>>& poly,
        int degree,
        const std::vector<std::vector<gsPoint2D<T>>>& loopPolys)
    {
        const gsPoint2D<T>& A = poly.front();
        const gsPoint2D<T>& B = poly.back();
        const T chordLen = A.dist(B);
        if (poly.size() <= 2 || chordLen <= T(1e-12))
            return makeStraightSeg(A, B, degree);

        // Check deviation from straight chord
        T maxDev = T(0);
        for (size_t i = 1; i + 1 < poly.size(); ++i) {
            T dev = std::abs(cross2(A, B, poly[i])) / chordLen;
            maxDev = std::max(maxDev, dev);
        }
        if (maxDev < T(0.03) * chordLen && isSegmentInsideDomain(A, B, loopPolys))
            return makeStraightSeg(A, B, degree);

        // Compute arc-length parameterization of polyline and sample at 0, 1/3, 2/3, 1
        std::vector<T> s(poly.size(), T(0));
        for (size_t i = 1; i < poly.size(); ++i)
            s[i] = s[i - 1] + poly[i - 1].dist(poly[i]);
        const T totS = s.back();
        if (totS <= T(1e-12))
            return makeStraightSeg(A, B, degree);

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

        gsPoint2D<T> q0 = A, q1 = sampleAt(T(1) / T(3)), q2 = sampleAt(T(2) / T(3)), q3 = B;
        gsPoint2D<T> p1((-T(5) * q0.x + T(18) * q1.x - T(9) * q2.x + T(2) * q3.x) / T(6),
                        (-T(5) * q0.y + T(18) * q1.y - T(9) * q2.y + T(2) * q3.y) / T(6));
        gsPoint2D<T> p2((T(2) * q0.x - T(9) * q1.x + T(18) * q2.x - T(5) * q3.x) / T(6),
                        (T(2) * q0.y - T(9) * q1.y + T(18) * q2.y - T(5) * q3.y) / T(6));

        gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
        gsMatrix<T> coefs(4, 2);
        coefs << q0.x, q0.y, p1.x, p1.y, p2.x, p2.y, q3.x, q3.y;
        gsBSpline<T> bsp(kv, coefs);
        if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());

        bool insideOk = true;
        for (int k = 1; k < 12; ++k) {
            gsPoint2D<T> q = evalCurvePoint(bsp, static_cast<T>(k) / T(12));
            if (!isPointStrictlyInsideDomain(q, loopPolys, T(1e-6))) { insideOk = false; break; }
        }
        if (insideOk) return CurveSeg<T>(bsp);
        return makeStraightSeg(A, B, degree);
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
        const T alphas[4] = {T(0.42), T(0.35), T(0.50), T(0.28)};

        std::vector<gsTensorBSpline<2, T>*> bestFive;
        T bestMaxAng = T(360);

        for (const auto& ctr : cands) {
            if (!isPointStrictlyInsideDomain(ctr, loopPolys, T(1e-6))) continue;
            for (T alpha : alphas) {
                gsPoint2D<T> Q[4];
                bool ok = true;
                for (int k = 0; k < 4; ++k) {
                    Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                        ctr.y + alpha * (P[k].y - ctr.y));
                    if (!isPointStrictlyInsideDomain(Q[k], loopPolys, T(1e-6))) { ok = false; break; }
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
                T maxA = T(0);
                for (auto* patch : trial) {
                    if (!patch || !isPatchDetSingleSign(*patch, 15, T(1e-7))) { trialOk = false; break; }
                    T angs[4], mnA = 0, mxA = 180;
                    bool degen = false;
                    evalPatchCornerAnglesDeg(*patch, angs, mnA, mxA, degen);
                    if (degen || mxA > T(162)) { trialOk = false; break; }
                    maxA = std::max(maxA, mxA);
                }
                if (trialOk && maxA < bestMaxAng) {
                    for (auto* old : bestFive) delete old;
                    bestFive = std::move(trial);
                    bestMaxAng = maxA;
                } else {
                    for (auto* p : trial) delete p;
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

            if (useCC) {
                node.pt = cc;
                node.radius = cr;
            } else if (node.triType == 2) {
                // J-triangle: incenter is guaranteed strictly inside the triangle
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
                    if (sepAngleDeg < thetaPruneDeg || branchLen < T(0.65) * triNodes[root].radius)
                        shouldPrune = true;
                } else if (sharpVerts.size() == 1) {
                    // Single sharp convex corner hair (e.g. corner of a rectangular arm end)
                    // Prune short corner hairs so blunt ends stay clean 2-corner caps instead of Y-splits
                    T cornerAng = mesh.cornerAngleDeg[sharpVerts[0]];
                    if (cornerAng >= T(65.0) && branchLen < T(1.15) * triNodes[root].radius)
                        shouldPrune = true;
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

                // If the leaf arc has 0 interior vertices (degenerate 1-edge cap) OR
                // both endpoints are sharp corners on a straight end wall, retreat leaf into parent
                bool bothSharpEnd = (mesh.isSharpCorner[vR] && mesh.isSharpCorner[vL] && sharpV.empty());
                if (nIntV == 0 || bothSharpEnd) {
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
                                if (analyzeCapArc(gvR, gvL, gSharp, gInt) && gSharp.size() <= 2) {
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

        // Contract very short J-J edges (0 sleeve triangles between two degree-3 J-triangles)
        // into a single degree-4 junction node when the merged quad is strictly convex.
        if (junctionMergeRatio > T(0)) {
            std::vector<bool> merged(macroNodes.size(), false);
            for (size_t i = 0; i < macroNodes.size(); ++i) {
                if (merged[i] || macroNodes[i].degree != 3 || macroNodes[i].tris.size() != 1) continue;
                int tA = macroNodes[i].tris[0];
                for (int nb : triNodes[tA].adj) {
                    int j = triToMacro[nb];
                    if (j <= static_cast<int>(i) || j < 0) continue;
                    if (merged[j] || macroNodes[j].degree != 3 || macroNodes[j].tris.size() != 1) continue;
                    T d = macroNodes[i].pt.dist(macroNodes[j].pt);
                    T rMin = std::min(macroNodes[i].radius, macroNodes[j].radius);
                    if (d >= junctionMergeRatio * rMin) continue;

                    // Verify merged 4-gon is strictly convex and midpoint has clear spokes
                    int vL = -1, vR = -1;
                    if (!getSharedChord(mesh, tA, nb, vL, vR)) continue;
                    int oppA = -1, oppB = -1;
                    for (int k = 0; k < 3; ++k) {
                        if (mesh.tris[tA].v[k] != vL && mesh.tris[tA].v[k] != vR) oppA = mesh.tris[tA].v[k];
                        if (mesh.tris[nb].v[k] != vL && mesh.tris[nb].v[k] != vR) oppB = mesh.tris[nb].v[k];
                    }
                    if (oppA < 0 || oppB < 0) continue;
                    gsPoint2D<T> q4[4] = {mesh.pts[oppA], mesh.pts[vR], mesh.pts[oppB], mesh.pts[vL]};
                    if (!isStrictlyConvexQuad(q4)) continue;

                    gsPoint2D<T> midPt(T(0.5) * (macroNodes[i].pt.x + macroNodes[j].pt.x),
                                       T(0.5) * (macroNodes[i].pt.y + macroNodes[j].pt.y));
                    if (!isPointStrictlyInsideDomain(midPt, mesh.loopPolys, T(1e-5))) continue;
                    if (cross2(q4[0], q4[1], midPt) <= T(1e-8) || cross2(q4[1], q4[2], midPt) <= T(1e-8) ||
                        cross2(q4[2], q4[3], midPt) <= T(1e-8) || cross2(q4[3], q4[0], midPt) <= T(1e-8))
                        continue;

                    // Contract j into i
                    macroNodes[i].pt = midPt;
                    macroNodes[i].radius = T(0.5) * (macroNodes[i].radius + macroNodes[j].radius);
                    macroNodes[i].tris.push_back(nb);
                    macroNodes[i].degree = 4;
                    triToMacro[nb] = static_cast<int>(i);
                    merged[j] = true;
                    macroNodes[j].tris.clear();
                    macroNodes[j].degree = 0;
                    break;
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
        // Step 1.5: True MAT Quad Emission (Branches + Bend Stations + Caps)
        // --------------------------------------------------------------------
        std::vector<gsTensorBSpline<2, T>*> quadPatches;
        std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> debugSpokes;
        bool emitFailed = false;

        auto pushQuad = [&](const CurveSeg<T>& e0, const CurveSeg<T>& e1,
                            const CurveSeg<T>& e2, const CurveSeg<T>& e3) -> bool {
            auto* patch = makeCoonsFromCCWCycle(e0, e1, e2, e3, degree);
            if (!patch) return false;
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
                CurveSeg<T> e0 = makeStraightSeg(Pa, pRa, degree);
                CurveSeg<T> e1 = makeStraightSeg(pRa, Pb, degree);
                CurveSeg<T> e2 = makeStraightSeg(Pb, pLa, degree);
                CurveSeg<T> e3 = makeStraightSeg(pLa, Pa, degree);
                return pushQuad(e0, e1, e2, e3);
            }

            // Case 2: Left wall collapsed to a single reflex corner R (vLa == vLb), Right wall advances
            // -> Split right CAD boundary at midpoint M_R and connect R -> M_R => 2 conforming quads!
            if (vLa == vLb && vRa != vRb) {
                CurveSeg<T> bndR = extractExactCADBoundaryBetween(mesh, vRa, vRb, degree, true);
                auto spR = bndR.splitArc(T(0.5));
                gsPoint2D<T> MR = spR.first.p1();
                if (!isSegmentInsideDomain(pLa, MR, mesh.loopPolys))
                    return false;
                debugSpokes.push_back({pLa, MR});

                CurveSeg<T> spPa_Ra = makeStraightSeg(Pa, pRa, degree);
                CurveSeg<T> spMR_R  = makeStraightSeg(MR, pLa, degree);
                CurveSeg<T> spR_Pa  = makeStraightSeg(pLa, Pa, degree);
                CurveSeg<T> spRb_Pb = makeStraightSeg(pRb, Pb, degree);
                CurveSeg<T> spPb_R  = makeStraightSeg(Pb, pLa, degree);

                // Quad 1 (CCW): Pa -> vRa -> MR -> R
                // Quad 2 (CCW): R -> MR -> vRb -> Pb
                return pushQuad(spPa_Ra, spR.first, spMR_R, spR_Pa) &&
                       pushQuad(spMR_R.reversed(), spR.second, spRb_Pb, spPb_R);
            }

            // Case 3: Right wall collapsed to a single reflex corner R (vRa == vRb), Left wall advances
            // -> Split left CAD boundary (vLb -> vLa forward) at midpoint M_L and connect R -> M_L => 2 quads!
            if (vRa == vRb && vLa != vLb) {
                CurveSeg<T> bndL = extractExactCADBoundaryBetween(mesh, vLb, vLa, degree, true);
                auto spL = bndL.splitArc(T(0.5));
                gsPoint2D<T> ML = spL.first.p1();
                if (!isSegmentInsideDomain(pRa, ML, mesh.loopPolys))
                    return false;
                debugSpokes.push_back({pRa, ML});

                CurveSeg<T> spPa_R  = makeStraightSeg(Pa, pRa, degree);
                CurveSeg<T> spR_ML  = makeStraightSeg(pRa, ML, degree);
                CurveSeg<T> spLa_Pa = makeStraightSeg(pLa, Pa, degree);
                CurveSeg<T> spR_Pb  = makeStraightSeg(pRa, Pb, degree);
                CurveSeg<T> spPb_Lb = makeStraightSeg(Pb, pLb, degree);

                // Quad 1 (CCW): Pa -> R -> ML -> vLa
                // Quad 2 (CCW): R -> Pb -> vLb -> ML
                return pushQuad(spPa_R, spR_ML, spL.second, spLa_Pa) &&
                       pushQuad(spR_Pb, spPb_Lb, spL.first, spR_ML.reversed());
            }

            // Case 4: Standard branch segment (both walls advance: vLa != vLb and vRa != vRb)
            // -> 2 conforming quads flanking the smooth medial curve Pa -> Pb
            CurveSeg<T> med  = fitMedialCurve(segPoly, degree, mesh.loopPolys);
            CurveSeg<T> bndL = extractExactCADBoundaryBetween(mesh, vLb, vLa, degree, true);
            CurveSeg<T> bndR = extractExactCADBoundaryBetween(mesh, vRa, vRb, degree, true);

            CurveSeg<T> spPb_Lb = makeStraightSeg(Pb, pLb, degree);
            CurveSeg<T> spLa_Pa = makeStraightSeg(pLa, Pa, degree);
            CurveSeg<T> spPa_Ra = makeStraightSeg(Pa, pRa, degree);
            CurveSeg<T> spRb_Pb = makeStraightSeg(pRb, Pb, degree);

            // Left quad (CCW):  Pa -> Pb -> vLb -> vLa
            // Right quad (CCW): Pa -> vRa -> vRb -> Pb
            return pushQuad(med, spPb_Lb, bndL, spLa_Pa) &&
                   pushQuad(spPa_Ra, bndR, spRb_Pb, med.reversed());
        };

        // Process all MacroBranches (inserting intermediate bend/aspect stations if curved or long)
        for (const auto& br : branches) {
            if (emitFailed) break;
            const int nChords = static_cast<int>(br.chords.size());
            if (nChords == 0) continue;

            // Choose station indices in [0 .. nChords-1] along br.chords
            // Station 0 is always macroNodes[br.uNode].pt with chord[0].
            // Last station is always macroNodes[br.vNode].pt with chord[nChords - 1].
            struct Station {
                gsPoint2D<T> pt;
                int vL, vR;
                int polyIdx;
            };
            std::vector<Station> stations;
            stations.push_back({macroNodes[br.uNode].pt, br.chords.front().vL, br.chords.front().vR, 0});

            if (nChords >= 3) {
                const T avgRad = std::max(T(1e-6), T(0.5) * (macroNodes[br.uNode].radius + macroNodes[br.vNode].radius));
                T accumLen = T(0);
                T accumTurn = T(0);
                gsPoint2D<T> prevDir(0, 0);
                bool hasDir = false;

                for (int ci = 1; ci < nChords - 1; ++ci) {
                    // poly index corresponding to chord ci is ci (poly[0]=u, poly[1..m]=triChain, poly[m+1]=v)
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

                    // Split if branch turned > 55 deg OR exceeded 2.8 * radius, and remaining tail is not tiny
                    bool wantSplit = (accumTurn > T(55.0)) || (accumLen > T(2.8) * avgRad);
                    if (wantSplit && (nChords - 1 - ci) >= 2) {
                        int vLi = br.chords[ci].vL;
                        int vRi = br.chords[ci].vR;
                        // Avoid splitting on an identical chord to previous station or end chord
                        if ((vLi != stations.back().vL || vRi != stations.back().vR) &&
                            (vLi != br.chords.back().vL || vRi != br.chords.back().vR))
                        {
                            // Choose station point strictly inside domain with clear spokes to vLi & vRi
                            gsPoint2D<T> candPt = curPt;
                            if (!isPointStrictlyInsideDomain(candPt, mesh.loopPolys, T(1e-5)) ||
                                !isSegmentInsideDomain(candPt, mesh.pts[vLi], mesh.loopPolys) ||
                                !isSegmentInsideDomain(candPt, mesh.pts[vRi], mesh.loopPolys))
                            {
                                candPt = gsPoint2D<T>(T(0.5) * (mesh.pts[vLi].x + mesh.pts[vRi].x),
                                                      T(0.5) * (mesh.pts[vLi].y + mesh.pts[vRi].y));
                            }
                            if (isPointStrictlyInsideDomain(candPt, mesh.loopPolys, T(1e-5)) &&
                                isSegmentInsideDomain(candPt, mesh.pts[vLi], mesh.loopPolys) &&
                                isSegmentInsideDomain(candPt, mesh.pts[vRi], mesh.loopPolys))
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
            const gsPoint2D<T>& pR = mesh.pts[vR];
            const gsPoint2D<T>& pL = mesh.pts[vL];

            if (sharpV.size() == 1) {
                // Subcase A: Single sharp tip C -> 1 kite quad (PL -> vR -> C -> vL)
                int vC = sharpV[0];
                CurveSeg<T> sp0  = makeStraightSeg(PL, pR, degree);
                CurveSeg<T> bnd1 = extractExactCADBoundaryBetween(mesh, vR, vC, degree, true);
                CurveSeg<T> bnd2 = extractExactCADBoundaryBetween(mesh, vC, vL, degree, true);
                CurveSeg<T> sp3  = makeStraightSeg(pL, PL, degree);
                if (!pushQuad(sp0, bnd1, bnd2, sp3)) emitFailed = true;
            } else if (sharpV.size() == 2) {
                // Subcase B: Blunt rectangular end (2 corners C1, C2) -> 2 corner quads
                // (reducePolychordStrips collapses the centerline across M_end into a 1-wide arm!)
                int vC1 = sharpV[0];
                int vC2 = sharpV[1];
                CurveSeg<T> endWall = extractExactCADBoundaryBetween(mesh, vC1, vC2, degree, true);
                auto spW = endWall.splitArc(T(0.5));
                gsPoint2D<T> Mend = spW.first.p1();
                debugSpokes.push_back({PL, Mend});

                CurveSeg<T> spPL_R   = makeStraightSeg(PL, pR, degree);
                CurveSeg<T> bndR_C1  = extractExactCADBoundaryBetween(mesh, vR, vC1, degree, true);
                CurveSeg<T> spM_PL   = makeStraightSeg(Mend, PL, degree);
                CurveSeg<T> bndC2_L  = extractExactCADBoundaryBetween(mesh, vC2, vL, degree, true);
                CurveSeg<T> spL_PL   = makeStraightSeg(pL, PL, degree);

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
                    makeStraightSeg(PL, pR, degree),
                    spArc.first,
                    spArc.second,
                    makeStraightSeg(pL, PL, degree)
                };
                if (!emitOGrid5Inside(capEdges, degree, mesh.loopPolys, quadPatches)) {
                    // Fallback: single quad if O-grid didn't fit (e.g. very shallow gentle arc)
                    if (!pushQuad(capEdges[0], capEdges[1], capEdges[2], capEdges[3]))
                        emitFailed = true;
                }
            } else {
                // Subcase D (rare: >= 3 sharp corners): fan radial quads from PL to consecutive corners
                // Split middle edges so total segments are even if needed, or emit 2 quads via middle corner
                int midIdx = static_cast<int>(sharpV.size()) / 2;
                int vMid = sharpV[midIdx];
                CurveSeg<T> spPL_R  = makeStraightSeg(PL, pR, degree);
                CurveSeg<T> spM_PL  = makeStraightSeg(mesh.pts[vMid], PL, degree);
                CurveSeg<T> spL_PL  = makeStraightSeg(pL, PL, degree);
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
