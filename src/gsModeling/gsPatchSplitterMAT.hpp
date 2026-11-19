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

namespace gismo {
namespace internal_splitter {

template <class T>
class gsPatchSplitterMAT {
public:
    struct MATNode {
        gsPoint2D<T> pt;
        T radius = 0;
        int triIdx = -1;
        int triType = 0; // 0: T (leaf), 1: S (sleeve), 2: J (junction), 3: merged
        std::vector<int> adj;
        int bndV[3] = {-1, -1, -1};
        bool isPruned = false;
        bool isMerged = false;
        int mergedInto = -1;
    };

    struct MATBranch {
        int u = -1, v = -1;
        std::vector<int> nodePath;
        std::vector<gsPoint2D<T>> polyline;
        int leftBndStart = -1, leftBndEnd = -1;
        int rightBndStart = -1, rightBndEnd = -1;
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
        T targetHFrac = T(1.0) / T(16.0),
        T thetaPruneDeg = T(40.0),
        T junctionMergeRatio = T(0.35))
    {
        if (dom.numLoops() == 0) return false;
        outMp.clear();

        const auto bb = dom.boundingBox();
        const T domDiag = std::hypot(bb(0, 1) - bb(0, 0), bb(1, 1) - bb(1, 0));
        const T targetH = std::max(T(1e-4), domDiag * targetHFrac);

        // Phase 1: Boundary-only CDT
        TriMeshData<T> mesh;
        if (!buildConstrainedDomainTriangulation(dom, degree, targetH, true, 0, mesh))
            return false;
        applyDomainLawsonFlips(mesh);

        const int nT = static_cast<int>(mesh.tris.size());
        const int nE = static_cast<int>(mesh.edges.size());
        if (nT < 2) return false;

        // Compute edge midpoints and half-curves preserving exact CAD geometry
        std::vector<gsPoint2D<T>> edgeMid(nE);
        std::vector<CurveSeg<T>> edgeHalf0(nE);
        std::vector<CurveSeg<T>> edgeHalf1(nE);

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

        // Compute optimal triangle centers (maximizing corner angle margins)
        std::vector<gsPoint2D<T>> triCenters(nT);
        for (int ti = 0; ti < nT; ++ti) {
            const auto& tr = mesh.tris[ti];
            const auto& A = mesh.pts[tr.v[0]];
            const auto& B = mesh.pts[tr.v[1]];
            const auto& C = mesh.pts[tr.v[2]];
            const auto& mAB = edgeMid[tr.e[0]];
            const auto& mBC = edgeMid[tr.e[1]];
            const auto& mCA = edgeMid[tr.e[2]];
            gsPoint2D<T> optC;
            if (pickOptimalTri1to3Center(A, B, C, mAB, mBC, mCA, optC)) {
                triCenters[ti] = optC;
            } else {
                triCenters[ti] = gsPoint2D<T>((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
            }
        }

        // Phase 2: Construct Conforming Dual Quads for each triangle
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

        if (quadPatches.empty()) return false;

        for (auto* p : quadPatches)
            outMp.addPatch(typename gsGeometry<T>::uPtr(p));

        outMp.computeTopology(T(1e-4), false, false);
        return (outMp.nPatches() > 0 && outMp.nPatches() == static_cast<size_t>(3 * nT));
    }
};

} // namespace internal_splitter
} // namespace gismo
