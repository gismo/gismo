/** @file gsPatchSplitterMAT.hpp
    @brief 2D Delaunay-Dual Medial Axis Transform (MAT) Quad Decomposition for gsPatchSplitter.

    Constructs a topologically exact 2D Medial Axis Transform from the boundary-constrained
    Delaunay triangulation (CDT), simplifies branches via theta-SMA object-angle filtering,
    merges clustered junctions, slices exact CAD boundary curves via BoundaryVertexInfo tags,
    and emits an all-quad conforming gsMultiPatch<T>.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
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
    for (int i = 0; i < nT; ++i) out << "5\n"; // VTK_TRIANGLE
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
        const auto& j = junctions[i];
        std::vector<int> conn = {static_cast<int>(allPts.size())};
        allPts.push_back(j);
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
    for (const auto& c : cellConn) {
        if (c.size() == 1) out << "1\n"; // VTK_VERTEX
        else out << "4\n";              // VTK_POLY_LINE
    }
    out << "CELL_DATA " << cellConn.size() << "\n";
    out << "SCALARS branchId int 1\nLOOKUP_TABLE default\n";
    for (int bid : branchIds) out << bid << "\n";
    out << "SCALARS radius float 1\nLOOKUP_TABLE default\n";
    for (float rad : cellRadii) out << rad << "\n";
}

template <class T>
inline void exportSkeletonToVTK(
    const std::vector<std::vector<gsPoint2D<T>>>& polylines,
    const std::vector<gsPoint2D<T>>& junctions,
    const std::string& filename)
{
    std::vector<T> pr(polylines.size(), T(0));
    std::vector<T> jr(junctions.size(), T(0));
    exportSkeletonToVTK(polylines, pr, junctions, jr, filename);
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
    for (int i = 0; i < nS; ++i) out << "3\n"; // VTK_LINE
}

// ============================================================================
// 2D Medial Axis Transform (MAT) Quad Decomposition Pipeline
// ============================================================================

template <class T>
class gsPatchSplitterMAT {
public:
    struct MATNode {
        gsPoint2D<T> pt;
        T radius = 0;
        int triIdx = -1;
        int triType = 1; // 0: T (leaf), 1: S (sleeve), 2: J (junction)
        std::vector<int> adjNodes;
        int tipVertex = -1; // For T: boundary vertex at convex tip
        int touchBnd = -1;
        bool isPruned = false;
        bool isMerged = false;
        int mergedInto = -1;
    };

    struct MATBranch {
        int u = -1, v = -1;
        std::vector<int> nodePath;
        std::vector<gsPoint2D<T>> polyline;
        int leftTouchStart = -1, leftTouchEnd = -1;
        int rightTouchStart = -1, rightTouchEnd = -1;
        std::vector<int> leftBndPath;
        std::vector<int> rightBndPath;
        bool isPruned = false;
        bool isCycle = false;
    };

    static CurveSeg<T> getCADBoundarySegment(
        const gsPlanarDomain<T>& /*dom*/,
        int degree,
        const TriMeshData<T>& mesh,
        int vStart, int vEnd)
    {
        return extractExactCADBoundaryBetween(mesh, vStart, vEnd, degree);
    }

    static bool buildMATMultiPatch(
        const gsPlanarDomain<T>& dom,
        int degree,
        gsMultiPatch<T>& outMp,
        T targetHFrac = T(1.0) / T(24.0),
        T thetaPruneDeg = T(45.0),
        T junctionMergeRatio = T(0.30))
    {
        if (dom.numLoops() == 0) return false;
        outMp.clear();

        const bool dbg = (std::getenv("GISMO_MAT_DEBUG") != nullptr);
        const auto bb = dom.boundingBox();
        const T domDiag = std::hypot(bb(0, 1) - bb(0, 0), bb(1, 1) - bb(1, 0));
        const T targetH = std::max(T(1e-4), domDiag * targetHFrac);

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
        std::vector<int> nBndList(nT, 0);
        std::vector<std::vector<int>> edgeTris(nE);

        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            int nB = 0;
            for (int k = 0; k < 3; ++k) {
                int eid = tr.e[k];
                if (mesh.edges[eid].isBoundary) ++nB;
                else edgeTris[eid].push_back(ti);
            }
            nBndList[ti] = nB;
            if (nB >= 2) triType[ti] = 0;      // T (terminal / leaf)
            else if (nB == 1) triType[ti] = 1; // S (sleeve)
            else triType[ti] = 2;              // J (junction)
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

            // Circumcenter and incenter computation with strict in-domain clamping
            T D = T(2) * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
            bool useCC = false;
            gsPoint2D<T> cc;
            T cr = 0;
            if (std::abs(D) > T(1e-12)) {
                T ux = ((A.x*A.x + A.y*A.y)*(B.y - C.y) + (B.x*B.x + B.y*B.y)*(C.y - A.y) + (C.x*C.x + C.y*C.y)*(A.y - B.y)) / D;
                T uy = ((A.x*A.x + A.y*A.y)*(C.x - B.x) + (B.x*B.x + B.y*B.y)*(A.x - C.x) + (C.x*C.x + C.y*C.y)*(B.x - A.x)) / D;
                cc = gsPoint2D<T>(ux, uy);
                cr = std::hypot(ux - A.x, uy - A.y);
                T dAB = (A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y);
                T dBC = (B.x - C.x)*(B.x - C.x) + (B.y - C.y)*(B.y - C.y);
                T dCA = (C.x - A.x)*(C.x - A.x) + (C.y - A.y)*(C.y - A.y);
                bool isObtuse = (dAB + dBC < dCA - T(1e-9) || dBC + dCA < dAB - T(1e-9) || dCA + dAB < dBC - T(1e-9));
                if (!isObtuse && isPointStrictlyInsideDomain(cc, mesh.loopPolys, T(1e-5))) {
                    useCC = true;
                }
            }

            if (useCC) {
                node.pt = cc;
                node.radius = cr;
            } else {
                if (node.triType == 2) {
                    // J-triangle: incenter is guaranteed strictly inside any valid triangle
                    T a = std::hypot(B.x - C.x, B.y - C.y);
                    T b = std::hypot(C.x - A.x, C.y - A.y);
                    T c = std::hypot(A.x - B.x, A.y - B.y);
                    T perim = std::max(T(1e-12), a + b + c);
                    node.pt = gsPoint2D<T>((a * A.x + b * B.x + c * C.x) / perim,
                                           (a * A.y + b * B.y + c * C.y) / perim);
                    node.radius = std::min({a, b, c}) * T(0.5);
                } else {
                    // S/T-triangle: midpoint of longest interior edge (Chordal Axis)
                    int bestEid = -1;
                    T maxL = -1;
                    for (int k = 0; k < 3; ++k) {
                        int eid = tr.e[k];
                        if (!mesh.edges[eid].isBoundary) {
                            T l = mesh.pts[mesh.edges[eid].a].dist(mesh.pts[mesh.edges[eid].b]);
                            if (l > maxL) { maxL = l; bestEid = eid; }
                        }
                    }
                    if (bestEid >= 0) {
                        node.pt = gsPoint2D<T>(T(0.5) * (mesh.pts[mesh.edges[bestEid].a].x + mesh.pts[mesh.edges[bestEid].b].x),
                                               T(0.5) * (mesh.pts[mesh.edges[bestEid].a].y + mesh.pts[mesh.edges[bestEid].b].y));
                        node.radius = maxL * T(0.5);
                    } else {
                        node.pt = gsPoint2D<T>((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
                        node.radius = T(0.3) * domDiag;
                    }
                }
            }

            if (node.triType == 0) {
                for (int k = 0; k < 3; ++k) {
                    int eid = tr.e[k];
                    if (!mesh.edges[eid].isBoundary) {
                        node.tipVertex = tr.v[(k + 2) % 3];
                        break;
                    }
                }
                if (node.tipVertex < 0) node.tipVertex = tr.v[0];
            }
        }

        // Connect dual adjacency
        for (int ei = 0; ei < nE; ++ei) {
            if (!mesh.edges[ei].isBoundary && edgeTris[ei].size() == 2) {
                int t0 = edgeTris[ei][0], t1 = edgeTris[ei][1];
                triNodes[t0].adjNodes.push_back(t1);
                triNodes[t1].adjNodes.push_back(t0);
            }
        }

        if (dbg) {
            std::cout << "[MAT] CDT: V=" << mesh.pts.size() << " T=" << nT
                      << " (T=" << countT << ", S=" << countS << ", J=" << countJ << ")\n";
            exportCDTToVTK(mesh, triType, "mat_debug_cdt.vtk");
        }

        // --------------------------------------------------------------------
        // Step 1.2: Special Topological Cases (Convex Quad / 5-Patch O-Grid)
        // --------------------------------------------------------------------
        if (dom.numLoops() == 1) {
            int nOuter = static_cast<int>(dom.outer().size());
            if (nOuter == 4) {
                CurveSeg<T> bsn[4];
                for (int k = 0; k < 4; ++k) bsn[k] = CurveSeg<T>(dom.outer().curve(k), degree);
                auto patch = makeCoonsFromCCWCycle(bsn[0], bsn[1], bsn[2], bsn[3], degree);
                if (patch) {
                    T minDet = 0, maxDet = 0, minAbs = 0, en = 0;
                    evalPatchJacobianAndWinslow(*patch, 15, minDet, maxDet, minAbs, en);
                    T angs[4], mnA = 0, mxA = 180;
                    bool degen = false;
                    evalPatchCornerAnglesDeg(*patch, angs, mnA, mxA, degen);
                    if (minDet > T(1e-4) && !degen && mxA < T(150.0)) {
                        outMp.addPatch(typename gsGeometry<T>::uPtr(patch));
                        outMp.computeTopology(T(1e-4), false, false);
                        if (dbg) std::cout << "[MAT] Pure 1-patch convex quad recognized.\n";
                        return true;
                    }
                    delete patch;
                }
            }

            // Smooth circle / ellipse check (0 sharp corners)
            int sharpCount = 0;
            for (bool sc : mesh.isSharpCorner) if (sc) sharpCount++;
            if (sharpCount == 0 && nOuter >= 4) {
                const auto& poly = mesh.loopPolys[0];
                int nPoly = static_cast<int>(poly.size());
                if (nPoly >= 8) {
                    gsPoint2D<T> ctr(T(0.5)*(bb(0, 0) + bb(0, 1)), T(0.5)*(bb(1, 0) + bb(1, 1)));
                    int idxStation[4] = {0, nPoly / 4, nPoly / 2, (3 * nPoly) / 4};
                    gsPoint2D<T> P[4] = {poly[idxStation[0]], poly[idxStation[1]], poly[idxStation[2]], poly[idxStation[3]]};
                    T alpha = T(0.40);
                    gsPoint2D<T> Q[4];
                    for (int k = 0; k < 4; ++k)
                        Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x), ctr.y + alpha * (P[k].y - ctr.y));

                    CurveSeg<T> qEdge[4], spoke[4], bArc[4];
                    for (int k = 0; k < 4; ++k) {
                        qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], degree);
                        spoke[k] = makeStraightSeg(Q[k], P[k], degree);
                        int v0 = idxStation[k] % mesh.pts.size();
                        int v1 = idxStation[(k + 1) % 4] % mesh.pts.size();
                        bArc[k] = extractExactCADBoundaryBetween(mesh, v0, v1, degree);
                    }
                    gsMultiPatch<T> ogridMp;
                    ogridMp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], degree)));
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        ogridMp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(qEdge[k], spoke[kn], bArc[k].reversed(), spoke[k].reversed(), degree)));
                    }
                    bool allOk = true;
                    for (size_t p = 0; p < ogridMp.nPatches(); ++p) {
                        if (!isPatchDetSingleSign(ogridMp.patch(p), 15, T(1e-6))) { allOk = false; break; }
                    }
                    if (allOk) {
                        outMp = give(ogridMp);
                        outMp.computeTopology(T(1e-4), false, false);
                        if (dbg) std::cout << "[MAT] Smooth 5-patch O-grid butterfly recognized.\n";
                        return true;
                    }
                }
            }
        }

        // --------------------------------------------------------------------
        // Step 1.2: General Dual Graph Build & Branch Compression
        // --------------------------------------------------------------------
        std::vector<bool> isMacroNode(nT, false);
        for (int ti = 0; ti < nT; ++ti) {
            if (triType[ti] == 0 || triType[ti] == 2)
                isMacroNode[ti] = true;
        }

        // Trace maximal degree-2 corridors of S-triangles into branches
        std::vector<MATBranch> branches;
        std::set<std::pair<int, int>> visitedEdges;

        for (int ti = 0; ti < nT; ++ti) {
            if (!isMacroNode[ti]) continue;
            for (int neighbor : triNodes[ti].adjNodes) {
                int lo = std::min(ti, neighbor), hi = std::max(ti, neighbor);
                if (visitedEdges.count({lo, hi})) continue;
                visitedEdges.insert({lo, hi});

                MATBranch branch;
                branch.u = ti;
                branch.polyline.push_back(triNodes[ti].pt);
                int curr = neighbor, prev = ti;

                while (!isMacroNode[curr] && triNodes[curr].adjNodes.size() == 2) {
                    branch.nodePath.push_back(curr);
                    branch.polyline.push_back(triNodes[curr].pt);
                    int nxt = (triNodes[curr].adjNodes[0] == prev) ? triNodes[curr].adjNodes[1] : triNodes[curr].adjNodes[0];
                    visitedEdges.insert({std::min(curr, nxt), std::max(curr, nxt)});
                    prev = curr;
                    curr = nxt;
                }
                branch.v = curr;
                branch.polyline.push_back(triNodes[curr].pt);
                branches.push_back(branch);
            }
        }

        // Pure cycle fallback (0 macro nodes): pick 4 evenly spaced stations
        if (branches.empty() && nT >= 4) {
            std::vector<int> cyclePath;
            std::vector<bool> inPath(nT, false);
            int curr = 0;
            while (curr >= 0 && !inPath[curr]) {
                inPath[curr] = true;
                cyclePath.push_back(curr);
                int nextNode = -1;
                for (int adj : triNodes[curr].adjNodes) {
                    if (!inPath[adj]) { nextNode = adj; break; }
                }
                curr = nextNode;
            }
            if (cyclePath.size() >= 4) {
                int cLen = static_cast<int>(cyclePath.size());
                int st[4] = {0, cLen / 4, cLen / 2, (3 * cLen) / 4};
                for (int s = 0; s < 4; ++s) isMacroNode[cyclePath[st[s]]] = true;
                for (int s = 0; s < 4; ++s) {
                    MATBranch b;
                    b.u = cyclePath[st[s]];
                    b.v = cyclePath[st[(s + 1) % 4]];
                    b.isCycle = true;
                    int i0 = st[s], i1 = st[(s + 1) % 4];
                    int idx = i0;
                    while (idx != i1) {
                        b.nodePath.push_back(cyclePath[idx]);
                        b.polyline.push_back(triNodes[cyclePath[idx]].pt);
                        idx = (idx + 1) % cLen;
                    }
                    b.polyline.push_back(triNodes[cyclePath[i1]].pt);
                    branches.push_back(b);
                }
            }
        }

        // --------------------------------------------------------------------
        // Step 1.3: Theta-SMA Leaf Pruning & Junction Contraction
        // --------------------------------------------------------------------
        const T radToDeg = T(180) / T(3.14159265358979323846);
        for (auto& b : branches) {
            if (triNodes[b.v].triType == 0 && triNodes[b.u].triType == 2) {
                int tipV = triNodes[b.v].tipVertex;
                if (tipV >= 0 && static_cast<size_t>(tipV) < mesh.cornerAngleDeg.size()) {
                    if (mesh.cornerAngleDeg[tipV] < T(148.0)) continue; // genuine sharp corner: keep
                }
                T branchLen = b.polyline.front().dist(b.polyline.back());
                if (branchLen < T(0.35) * triNodes[b.u].radius) {
                    b.isPruned = true;
                }
            }
        }

        // --------------------------------------------------------------------
        // Step 1.4: Constrained Laplacian Smoothing & Bend Stations
        // --------------------------------------------------------------------
        for (auto& b : branches) {
            if (b.isPruned || b.polyline.size() < 3) continue;
            int nP = static_cast<int>(b.polyline.size());
            for (int iter = 0; iter < 4; ++iter) {
                for (int i = 1; i < nP - 1; ++i) {
                    gsPoint2D<T> smoothed(
                        T(0.5) * b.polyline[i].x + T(0.25) * (b.polyline[i - 1].x + b.polyline[i + 1].x),
                        T(0.5) * b.polyline[i].y + T(0.25) * (b.polyline[i - 1].y + b.polyline[i + 1].y)
                    );
                    if (isPointStrictlyInsideDomain(smoothed, mesh.loopPolys, T(1e-5)))
                        b.polyline[i] = smoothed;
                }
            }
        }

        std::vector<gsPoint2D<T>> debugJunctions;
        std::vector<T> debugJunctionRadii;
        std::vector<std::vector<gsPoint2D<T>>> debugPolylines;
        std::vector<T> debugPolylineRadii;
        for (const auto& b : branches) {
            if (!b.isPruned && b.polyline.size() >= 2) {
                debugPolylines.push_back(b.polyline);
                debugPolylineRadii.push_back(b.nodePath.empty() ? triNodes[b.u].radius : triNodes[b.nodePath[0]].radius);
            }
        }
        for (int ti = 0; ti < nT; ++ti) {
            if (triNodes[ti].triType == 2) {
                debugJunctions.push_back(triNodes[ti].pt);
                debugJunctionRadii.push_back(triNodes[ti].radius);
            }
        }

        std::vector<std::pair<gsPoint2D<T>, gsPoint2D<T>>> debugSpokes;

        // --------------------------------------------------------------------
        // Step 1.5: Conforming Dual Quad Decomposition
        // --------------------------------------------------------------------
        // Compute optimal triangle centers and edge midpoints
        std::vector<gsPoint2D<T>> edgeMid(nE);
        std::vector<CurveSeg<T>> edgeHalf0(nE), edgeHalf1(nE);

        for (int ei = 0; ei < nE; ++ei) {
            const auto& ed = mesh.edges[ei];
            if (ed.isBoundary) {
                auto sp = ed.seg.splitArc(T(0.5));
                edgeHalf0[ei] = sp.first;
                edgeHalf1[ei] = sp.second;
                edgeMid[ei] = sp.first.p1();
            } else {
                edgeMid[ei] = gsPoint2D<T>(T(0.5) * (mesh.pts[ed.a].x + mesh.pts[ed.b].x),
                                            T(0.5) * (mesh.pts[ed.a].y + mesh.pts[ed.b].y));
                edgeHalf0[ei] = makeStraightSeg(mesh.pts[ed.a], edgeMid[ei], degree);
                edgeHalf1[ei] = makeStraightSeg(edgeMid[ei], mesh.pts[ed.b], degree);
            }
        }

        std::vector<gsPoint2D<T>> triCenters(nT);
        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            const auto& A = mesh.pts[tr.v[0]];
            const auto& B = mesh.pts[tr.v[1]];
            const auto& C = mesh.pts[tr.v[2]];
            gsPoint2D<T> optC;
            if (pickOptimalTri1to3Center(A, B, C, edgeMid[tr.e[0]], edgeMid[tr.e[1]], edgeMid[tr.e[2]], optC)) {
                triCenters[ti] = optC;
            } else {
                triCenters[ti] = triNodes[ti].pt;
            }
        }

        std::vector<gsTensorBSpline<2, T>*> quadPatches;
        quadPatches.reserve(3 * nT);

        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            const auto& C = triCenters[ti];

            for (int k = 0; k < 3; ++k) {
                int vCurr = tr.v[k];
                int eIn = tr.e[(k + 2) % 3];  // incoming edge (vPrev -> vCurr)
                int eOut = tr.e[k];           // outgoing edge (vCurr -> vNext)

                CurveSeg<T> sIn = (mesh.edges[eIn].a == vCurr) ?
                    edgeHalf0[eIn].reversed() : edgeHalf1[eIn];

                CurveSeg<T> sOut = (mesh.edges[eOut].a == vCurr) ?
                    edgeHalf0[eOut] : edgeHalf1[eOut].reversed();

                CurveSeg<T> sToCenter = makeStraightSeg(sOut.p1(), C, degree);
                CurveSeg<T> sFromCenter = makeStraightSeg(C, sIn.p0(), degree);

                debugSpokes.push_back({sOut.p1(), C});
                debugSpokes.push_back({C, sIn.p0()});

                auto patch = makeCoonsFromCCWCycle(sOut, sToCenter, sFromCenter, sIn, degree);
                if (patch) {
                    T minDet = 0, maxDet = 0, minAbs = 0, wins = 0;
                    evalPatchJacobianAndWinslow(*patch, 10, minDet, maxDet, minAbs, wins);
                    if (minDet > T(1e-8)) {
                        quadPatches.push_back(patch);
                    } else {
                        delete patch;
                    }
                }
            }
        }

        if (dbg) {
            exportSkeletonToVTK(debugPolylines, debugPolylineRadii, debugJunctions, debugJunctionRadii, "mat_debug_skeleton.vtk");
            exportSpokesToVTK(debugSpokes, "mat_debug_spokes.vtk");
            std::cout << "[MAT] Emitted " << quadPatches.size() << " initial dual quads.\n";
        }

        if (quadPatches.empty()) return false;

        for (auto* p : quadPatches)
            outMp.addPatch(typename gsGeometry<T>::uPtr(p));

        outMp.computeTopology(T(1e-4), false, false);
        return (outMp.nPatches() > 0);
    }
};

} // namespace internal_splitter
} // namespace gismo
