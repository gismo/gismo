/** @file gsPatchSplitterCrossField.hpp
    @brief 4-RoSy Harmonic Cross-Field Solver, Separatrix Tracing, Strip Collapse,
           and Conforming Arrangement Quadification for gsPatchSplitter<T>.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitterTriMesh.hpp>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseSolver.h>

#include <complex>

namespace gismo {
namespace internal_splitter {

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
    auto testTri = [&](int ti) -> bool {
        if (ti < 0 || ti >= static_cast<int>(mesh.tris.size())) return false;
        const auto& t = mesh.tris[ti];
        const auto& A = mesh.pts[t.v[0]];
        const auto& B = mesh.pts[t.v[1]];
        const auto& C = mesh.pts[t.v[2]];
        T area2 = cross2(A, B, C);
        if (std::abs(area2) < T(1e-14)) return false;
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
        return false;
    };

    if (triOut >= 0 && testTri(triOut)) return true;

    for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
        if (ti == triOut) continue;
        if (testTri(ti)) return true;
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

        gsPoint2D<T> tFwd(fineMesh.bndTanOut[i].x, fineMesh.bndTanOut[i].y);
        const gsPoint2D<T>& tIn = fineMesh.bndTanIn[i];
        T cr = tIn.x * tFwd.y - tIn.y * tFwd.x;
        T dt = tIn.x * tFwd.x + tIn.y * tFwd.y;
        T turnRad = std::atan2(cr, dt);
        T angIfLeft = T(180) - turnRad * (T(180) / T(3.14159265358979323846));
        T sweepSign = (std::abs(fineMesh.cornerAngleDeg[i] - angIfLeft) < T(1)) ? T(1) : T(-1);
        T angStart = std::atan2(tFwd.y, tFwd.x);
        T sweepRad = fineMesh.cornerAngleDeg[i] * T(3.14159265358979323846) / T(180);
        int nPorts = std::max(1, static_cast<int>(std::floor(fineMesh.cornerAngleDeg[i] / T(90) + T(0.5))) - 1);
        for (int p = 1; p <= nPorts; ++p) {
            T ang = angStart + sweepSign * sweepRad * (static_cast<T>(p) / static_cast<T>(nPorts + 1));
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
        int curTri = -1;

        for (int st = 0; st < maxSteps && !terminated; ++st) {
            // Update direction from cross-field (branch closest to current dir)
            T th = T(0);
            if (st > 2 && evalCrossFieldAngleAt(fineMesh, cf, cur, th, curTri)) {
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
