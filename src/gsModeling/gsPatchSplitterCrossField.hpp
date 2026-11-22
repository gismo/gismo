/** @file gsPatchSplitterCrossField.hpp
    @brief MBO Unit-Norm 4-RoSy Harmonic Cross-Field Solver, O(1) Triangle-Walk Evaluator,
           Exact Singularity Port Extraction, RK2 Curved Separatrix Arrangement with
           Corridor Merging & Local Even-Gon Resolution, and Cross-Field-Aligned
           Coarse Blossom-Quad Decomposition with Dual-Path Zipper & Curved Hermite Rescue.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitterTriMesh.hpp>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseSolver.h>

#include <array>
#include <complex>
#include <functional>
#include <queue>
#include <tuple>

namespace gismo {
namespace internal_splitter {

// ============================================================================
// Fast Triangle-Adjacency Accelerator for O(1) Cross-Field Point Location
// ============================================================================

template <class T>
struct CrossFieldMeshAccel {
    std::vector<std::array<int, 3>> triAdj;   // triAdj[ti][k] = neighbor tri across edge k (-1 if boundary)
    std::vector<std::array<T, 4>>   triBBox;  // {xmin, xmax, ymin, ymax}

    void build(const TriMeshData<T>& mesh)
    {
        const int nT = static_cast<int>(mesh.tris.size());
        const int nE = static_cast<int>(mesh.edges.size());
        triAdj.assign(nT, {-1, -1, -1});
        triBBox.resize(nT);

        std::vector<std::pair<int, int>> eFirst(nE, {-1, -1});
        for (int ti = 0; ti < nT; ++ti) {
            const auto& t = mesh.tris[ti];
            const auto& A = mesh.pts[t.v[0]];
            const auto& B = mesh.pts[t.v[1]];
            const auto& C = mesh.pts[t.v[2]];
            triBBox[ti] = {
                std::min(A.x, std::min(B.x, C.x)) - T(1e-7),
                std::max(A.x, std::max(B.x, C.x)) + T(1e-7),
                std::min(A.y, std::min(B.y, C.y)) - T(1e-7),
                std::max(A.y, std::max(B.y, C.y)) + T(1e-7)
            };
            for (int k = 0; k < 3; ++k) {
                int eid = t.e[k];
                if (eid < 0 || eid >= nE) continue;
                if (eFirst[eid].first < 0) {
                    eFirst[eid] = {ti, k};
                } else {
                    int ot = eFirst[eid].first;
                    int ok = eFirst[eid].second;
                    triAdj[ti][k] = ot;
                    triAdj[ot][ok] = ti;
                }
            }
        }
    }
};

// ============================================================================
// 4-RoSy MBO Unit-Norm Harmonic Cross-Field Solver & Singularity Extraction
// ============================================================================

template <class T>
struct CrossFieldSolution {
    std::vector<std::complex<T>> z;
    struct Singularity {
        gsPoint2D<T> pos;
        int index = 0; // +1 (val 3), -1 (val 5), +4 (O-grid disk), valence = 4 - index
        std::vector<T> portAngles;
    };
    std::vector<Singularity> singularities;
};

// O(1) amortized triangle-walk evaluation of complex 4-RoSy field z and cross angle theta in [-pi/4, pi/4]
template <class T>
inline bool evalCrossFieldAt(const TriMeshData<T>& mesh,
                             const CrossFieldSolution<T>& cf,
                             const CrossFieldMeshAccel<T>& accel,
                             const gsPoint2D<T>& p,
                             std::complex<T>& zOut,
                             T& thetaOut,
                             int& triHint)
{
    const int nT = static_cast<int>(mesh.tris.size());
    if (nT == 0 || cf.z.size() != mesh.pts.size()) return false;

    auto baryEval = [&](int ti, T& l0, T& l1, T& l2) -> bool {
        if (ti < 0 || ti >= nT) return false;
        const auto& t = mesh.tris[ti];
        const auto& A = mesh.pts[t.v[0]];
        const auto& B = mesh.pts[t.v[1]];
        const auto& C = mesh.pts[t.v[2]];
        T area2 = cross2(A, B, C);
        if (std::abs(area2) < T(1e-14)) return false;
        l0 = cross2(B, C, p) / area2;
        l1 = cross2(C, A, p) / area2;
        l2 = T(1) - l0 - l1;
        return true;
    };

    auto finalizeHit = [&](int ti, T l0, T l1, T l2) -> bool {
        const auto& t = mesh.tris[ti];
        zOut = l0 * cf.z[t.v[0]] + l1 * cf.z[t.v[1]] + l2 * cf.z[t.v[2]];
        if (std::abs(zOut) < T(1e-12)) return false;
        thetaOut = T(0.25) * std::arg(zOut);
        triHint = ti;
        return true;
    };

    if (triHint >= 0 && triHint < nT && static_cast<int>(accel.triAdj.size()) == nT) {
        int cur = triHint;
        int prev = -1;
        for (int step = 0; step < 64; ++step) {
            T l0 = 0, l1 = 0, l2 = 0;
            if (!baryEval(cur, l0, l1, l2)) break;
            if (l0 >= -T(1e-4) && l1 >= -T(1e-4) && l2 >= -T(1e-4))
                return finalizeHit(cur, l0, l1, l2);
            int edgeOrder[3] = {1, 2, 0};
            T vals[3] = {l0, l1, l2};
            if (vals[1] < vals[0]) { std::swap(vals[0], vals[1]); std::swap(edgeOrder[0], edgeOrder[1]); }
            if (vals[2] < vals[1]) { std::swap(vals[1], vals[2]); std::swap(edgeOrder[1], edgeOrder[2]); }
            if (vals[1] < vals[0]) { std::swap(vals[0], vals[1]); std::swap(edgeOrder[0], edgeOrder[1]); }

            int nxt = -1;
            for (int m = 0; m < 3; ++m) {
                if (vals[m] >= -T(1e-4)) break;
                int nb = accel.triAdj[cur][edgeOrder[m]];
                if (nb >= 0 && nb != prev) { nxt = nb; break; }
            }
            if (nxt < 0) break;
            prev = cur;
            cur = nxt;
        }
    }

    const bool useBB = (static_cast<int>(accel.triBBox.size()) == nT);
    for (int ti = 0; ti < nT; ++ti) {
        if (useBB) {
            const auto& bb = accel.triBBox[ti];
            if (p.x < bb[0] || p.x > bb[1] || p.y < bb[2] || p.y > bb[3]) continue;
        }
        T l0 = 0, l1 = 0, l2 = 0;
        if (baryEval(ti, l0, l1, l2) && l0 >= -T(1e-4) && l1 >= -T(1e-4) && l2 >= -T(1e-4))
            return finalizeHit(ti, l0, l1, l2);
    }
    return false;
}

template <class T>
inline bool evalCrossFieldAngleAt(const TriMeshData<T>& mesh,
                                  const CrossFieldSolution<T>& cf,
                                  const gsPoint2D<T>& p,
                                  T& thetaOut,
                                  int& triOut)
{
    CrossFieldMeshAccel<T> emptyAccel;
    std::complex<T> zDummy;
    return evalCrossFieldAt(mesh, cf, emptyAccel, p, zDummy, thetaOut, triOut);
}

// Cross-field alignment score of directed segment A -> B in [-1, +1]
template <class T>
inline T evalSegmentCrossFieldAlignment(const TriMeshData<T>& fineMesh,
                                        const CrossFieldSolution<T>& cf,
                                        const CrossFieldMeshAccel<T>& accel,
                                        const gsPoint2D<T>& A,
                                        const gsPoint2D<T>& B,
                                        int& triHint)
{
    T dx = B.x - A.x, dy = B.y - A.y;
    if (dx * dx + dy * dy <= T(1e-24)) return T(0);
    gsPoint2D<T> mid(T(0.5) * (A.x + B.x), T(0.5) * (A.y + B.y));
    std::complex<T> zVal;
    T theta = T(0);
    if (!evalCrossFieldAt(fineMesh, cf, accel, mid, zVal, theta, triHint))
        return T(0);
    T segAng = std::atan2(dy, dx);
    return std::cos(T(4) * (segAng - theta));
}

template <class T>
inline CrossFieldSolution<T> solveHarmonic4RoSy(const TriMeshData<T>& mesh,
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
        } else if (std::abs(z1 - z2) < T(0.40)) {
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

    T avgEdge = T(0);
    for (const auto& e : mesh.edges)
        avgEdge += mesh.pts[e.a].dist(mesh.pts[e.b]);
    avgEdge /= std::max<T>(T(1), static_cast<T>(mesh.edges.size()));

    if (nFree == 0) {
        sol.z = dirVal;
    } else {
        gsSparseEntries<T> lapEntries;
        lapEntries.reserve(mesh.tris.size() * 12);
        gsMatrix<T> rhsLap = gsMatrix<T>::Zero(nFree, 2);
        gsVector<T> massDiag = gsVector<T>::Zero(nFree);

        auto addWeight = [&](int u, int v, T w) {
            if (isDirichlet[u]) return;
            int iu = freeIdx[u];
            lapEntries.add(iu, iu, w);
            if (isDirichlet[v]) {
                rhsLap(iu, 0) += w * dirVal[v].real();
                rhsLap(iu, 1) += w * dirVal[v].imag();
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
            for (int k = 0; k < 3; ++k) {
                if (!isDirichlet[t.v[k]])
                    massDiag(freeIdx[t.v[k]]) += mTri;
            }

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

        // Step 1: Initial harmonic Dirichlet solve L_ff z_free = -L_fb z_bnd
        gsSparseMatrix<T> L(nFree, nFree);
        L.setFrom(lapEntries);
        typename gsSparseSolver<T>::SimplicialLDLT lapSolver;
        lapSolver.compute(L);
        gsMatrix<T> solMat = lapSolver.solve(rhsLap);

        // Step 2: MBO Unit-Norm Diffusion-Projection (Viertel & Osting)
        const T tau = T(2.5) * avgEdge * avgEdge;
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
                sol.z[i] = dirVal[i];
            } else {
                int ii = freeIdx[i];
                sol.z[i] = std::complex<T>(solMat(ii, 0), solMat(ii, 1));
            }
        }
    }

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

    // Cluster nearby singularities & cancel opposite-sign dipoles
    const T clusterTol = T(1.8) * avgEdge;
    std::vector<bool> used(rawSings.size(), false);
    std::vector<RawSing> clustered;
    for (size_t i = 0; i < rawSings.size(); ++i) {
        if (used[i]) continue;
        used[i] = true;
        gsPoint2D<T> sumP = rawSings[i].pos;
        int sumInd = rawSings[i].idx;
        int cnt = 1;
        for (size_t j = i + 1; j < rawSings.size(); ++j) {
            if (!used[j] && rawSings[i].pos.dist(rawSings[j].pos) <= clusterTol &&
                isSegmentInsideDomain(rawSings[i].pos, rawSings[j].pos, mesh.loopPolys)) {
                used[j] = true;
                sumP.x += rawSings[j].pos.x;
                sumP.y += rawSings[j].pos.y;
                sumInd += rawSings[j].idx;
                ++cnt;
            }
        }
        if (sumInd != 0)
            clustered.push_back({gsPoint2D<T>(sumP.x / cnt, sumP.y / cnt), sumInd, rawSings[i].triIdx});
    }

    const T dipoleTol = T(2.6) * avgEdge;
    std::vector<bool> dead(clustered.size(), false);
    for (size_t i = 0; i < clustered.size(); ++i) {
        if (dead[i]) continue;
        int bestJ = -1;
        T bestD = dipoleTol;
        for (size_t j = i + 1; j < clustered.size(); ++j) {
            if (dead[j] || clustered[i].idx + clustered[j].idx != 0) continue;
            T d = clustered[i].pos.dist(clustered[j].pos);
            if (d < bestD && isSegmentInsideDomain(clustered[i].pos, clustered[j].pos, mesh.loopPolys)) {
                bestD = d;
                bestJ = static_cast<int>(j);
            }
        }
        if (bestJ >= 0) { dead[i] = true; dead[bestJ] = true; }
    }

    // Exact circular-harmonic port extraction: phi_p = (phi0 + 2*pi*p) / (4 - k)
    for (size_t i = 0; i < clustered.size(); ++i) {
        if (dead[i]) continue;
        typename CrossFieldSolution<T>::Singularity sg;
        sg.pos = clustered[i].pos;
        sg.index = clustered[i].idx;
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
// Polyline-to-Curved-BSpline Fitter & Local Even-Gon Chord Splitting Helpers
// ============================================================================

template <class T>
inline CurveSeg<T> fitStreamlinePolylineToCurveSeg(
    const std::vector<gsPoint2D<T>>& poly,
    int degree,
    const std::vector<std::vector<gsPoint2D<T>>>& loopPolys)
{
    if (poly.size() < 2)
        return CurveSeg<T>();
    const gsPoint2D<T>& A = poly.front();
    const gsPoint2D<T>& B = poly.back();
    const T chord = A.dist(B);
    CurveSeg<T> straight = makeStraightSeg(A, B, degree);
    if (poly.size() <= 2 || chord <= T(1e-12))
        return straight;

    T maxDev = T(0);
    for (size_t i = 1; i + 1 < poly.size(); ++i)
        maxDev = std::max(maxDev, std::abs(cross2(A, B, poly[i])) / chord);
    if (maxDev < T(0.015) * chord && isSegmentInsideDomain(A, B, loopPolys))
        return straight;

    std::vector<T> s(poly.size(), T(0));
    for (size_t i = 1; i < poly.size(); ++i)
        s[i] = s[i - 1] + poly[i - 1].dist(poly[i]);
    const T totS = s.back();
    if (totS <= T(1e-12)) return straight;

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
    gsPoint2D<T> p1Int((-T(5) * A.x + T(18) * q1.x - T(9) * q2.x + T(2) * B.x) / T(6),
                       (-T(5) * A.y + T(18) * q1.y - T(9) * q2.y + T(2) * B.y) / T(6));
    gsPoint2D<T> p2Int((T(2) * A.x - T(9) * q1.x + T(18) * q2.x - T(5) * B.x) / T(6),
                       (T(2) * A.y - T(9) * q1.y + T(18) * q2.y - T(5) * B.y) / T(6));
    gsPoint2D<T> p1Lin((T(2) * A.x + B.x) / T(3), (T(2) * A.y + B.y) / T(3));
    gsPoint2D<T> p2Lin((A.x + T(2) * B.x) / T(3), (A.y + T(2) * B.y) / T(3));

    auto curveInside = [&](const CurveSeg<T>& seg) -> bool {
        for (int k = 1; k < 12; ++k) {
            gsPoint2D<T> pt = evalCurvePoint(*seg.bsp, static_cast<T>(k) / T(12));
            if (!isPointStrictlyInsideDomain(pt, loopPolys, T(1e-6))) return false;
        }
        return true;
    };

    const std::pair<gsPoint2D<T>, gsPoint2D<T>> basePairs[2] = {{p1Int, p2Int}, {q1, q2}};
    for (int bi = 0; bi < 2; ++bi) {
        for (T lam : {T(1.0), T(0.75), T(0.50), T(0.25)}) {
            gsPoint2D<T> c1((T(1) - lam) * p1Lin.x + lam * basePairs[bi].first.x,
                            (T(1) - lam) * p1Lin.y + lam * basePairs[bi].first.y);
            gsPoint2D<T> c2((T(1) - lam) * p2Lin.x + lam * basePairs[bi].second.x,
                            (T(1) - lam) * p2Lin.y + lam * basePairs[bi].second.y);
            CurveSeg<T> cand = makeCubicBezierInteriorSeg(A, c1, c2, B, degree);
            if (curveInside(cand))
                return cand;
        }
    }
    return straight;
}

template <class T>
inline bool splitEvenFaceToQuadsRecursive(
    const std::vector<CurveSeg<T>>& edgesCCW,
    const std::vector<gsPoint2D<T>>& vertsCCW,
    const TriMeshData<T>& fineMesh,
    const CrossFieldSolution<T>& cf,
    const CrossFieldMeshAccel<T>& accel,
    int degree,
    std::vector<std::array<CurveSeg<T>, 4>>& outQuads)
{
    const int N = static_cast<int>(edgesCCW.size());
    if (N < 4 || (N % 2) != 0) return false;
    if (N == 4) {
        T mnA = 0, mxA = 180, ef = 0;
        eval4SegCornerAnglesDeg(edgesCCW[0], edgesCCW[1], edgesCCW[2], edgesCCW[3], mnA, mxA, ef);
        if (mnA < T(8.0) || mxA > T(168.0)) return false;
        outQuads.push_back({edgesCCW[0], edgesCCW[1], edgesCCW[2], edgesCCW[3]});
        return true;
    }

    std::vector<gsPoint2D<T>> facePoly;
    facePoly.reserve(N * 8);
    for (int i = 0; i < N; ++i)
        for (int s = 0; s < 8; ++s)
            facePoly.push_back(evalCurvePoint(*edgesCCW[i].bsp, static_cast<T>(s) / T(8)));

    const T radToDeg = T(180) / T(3.14159265358979323846);
    auto vecAng = [&](const gsPoint2D<T>& u, const gsPoint2D<T>& v) -> T {
        T nu = std::hypot(u.x, u.y), nv = std::hypot(v.x, v.y);
        if (nu <= T(1e-14) || nv <= T(1e-14)) return T(0);
        T a = std::atan2(u.x * v.y - u.y * v.x, u.x * v.x + u.y * v.y) * radToDeg;
        if (a <= T(0)) a += T(360);
        return a;
    };
    auto segDirOut = [](const CurveSeg<T>& s) -> gsPoint2D<T> {
        auto d = evalCurveDeriv(*s.bsp, T(0));
        T n = std::hypot(d.x, d.y);
        return (n > T(1e-14)) ? gsPoint2D<T>(d.x / n, d.y / n) : gsPoint2D<T>(1, 0);
    };
    auto segDirBack = [](const CurveSeg<T>& s) -> gsPoint2D<T> {
        auto d = evalCurveDeriv(*s.bsp, T(1));
        T n = std::hypot(d.x, d.y);
        return (n > T(1e-14)) ? gsPoint2D<T>(-d.x / n, -d.y / n) : gsPoint2D<T>(-1, 0);
    };

    struct DiagCand { int i, j; T score; };
    std::vector<DiagCand> cands;
    int triHint = -1;

    for (int i = 0; i < N; ++i) {
        for (int step = 3; step <= N - 3; step += 2) {
            int j = (i + step) % N;
            if (j < i) continue;
            const gsPoint2D<T>& A = vertsCCW[i];
            const gsPoint2D<T>& B = vertsCCW[j];
            T L = A.dist(B);
            if (L <= T(1e-8)) continue;
            if (!isSegmentInsideDomain(A, B, fineMesh.loopPolys)) continue;

            bool insideFace = true;
            for (int k = 1; k <= 5; ++k) {
                T f = static_cast<T>(k) / T(6);
                gsPoint2D<T> q((T(1) - f) * A.x + f * B.x, (T(1) - f) * A.y + f * B.y);
                if (!pointInPolyStrict(q, facePoly, T(1e-6))) { insideFace = false; break; }
            }
            if (!insideFace) continue;
            int np = static_cast<int>(facePoly.size());
            for (int k = 0; k < np && insideFace; ++k) {
                const auto& p0 = facePoly[k];
                const auto& p1 = facePoly[(k + 1) % np];
                if (p0.dist(A) < T(1e-5) || p1.dist(A) < T(1e-5) ||
                    p0.dist(B) < T(1e-5) || p1.dist(B) < T(1e-5)) continue;
                if (segProperIntersect(A, B, p0, p1)) insideFace = false;
            }
            if (!insideFace) continue;

            gsPoint2D<T> dAB((B.x - A.x) / L, (B.y - A.y) / L);
            gsPoint2D<T> dBA(-dAB.x, -dAB.y);
            T a1 = vecAng(segDirOut(edgesCCW[i]), dAB);
            T a2 = vecAng(dAB, segDirBack(edgesCCW[(i + N - 1) % N]));
            T a3 = vecAng(segDirOut(edgesCCW[j]), dBA);
            T a4 = vecAng(dBA, segDirBack(edgesCCW[(j + N - 1) % N]));
            T mnA = std::min(std::min(a1, a2), std::min(a3, a4));
            T mxA = std::max(std::max(a1, a2), std::max(a3, a4));
            if (mnA < T(12.0) || mxA > T(162.0)) continue;

            T effMin = std::min(mnA, T(180) - mxA);
            T cfAl = evalSegmentCrossFieldAlignment(fineMesh, cf, accel, A, B, triHint);
            cands.push_back({i, j, effMin + T(8.0) * cfAl});
        }
    }

    std::sort(cands.begin(), cands.end(), [](const DiagCand& a, const DiagCand& b) { return a.score > b.score; });

    for (size_t ci = 0; ci < std::min<size_t>(cands.size(), 6); ++ci) {
        int i = cands[ci].i, j = cands[ci].j;
        CurveSeg<T> diag = makeStraightSeg(vertsCCW[i], vertsCCW[j], degree);

        std::vector<CurveSeg<T>> e1, e2;
        std::vector<gsPoint2D<T>> v1, v2;
        for (int k = i; k != j; k = (k + 1) % N) { e1.push_back(edgesCCW[k]); v1.push_back(vertsCCW[k]); }
        e1.push_back(diag.reversed()); v1.push_back(vertsCCW[j]);

        for (int k = j; k != i; k = (k + 1) % N) { e2.push_back(edgesCCW[k]); v2.push_back(vertsCCW[k]); }
        e2.push_back(diag); v2.push_back(vertsCCW[i]);

        std::vector<std::array<CurveSeg<T>, 4>> subOut;
        if (splitEvenFaceToQuadsRecursive(e1, v1, fineMesh, cf, accel, degree, subOut) &&
            splitEvenFaceToQuadsRecursive(e2, v2, fineMesh, cf, accel, degree, subOut)) {
            for (auto& q : subOut) outQuads.push_back(std::move(q));
            return true;
        }
    }
    return false;
}

// ============================================================================
// Engine 1: Curved-Streamline MBO Cross-Field Separatrix Arrangement
// ============================================================================

template <class T>
inline bool buildCrossFieldSeparatrixFromSolved(
    const gsPlanarDomain<T>& dom,
    const TriMeshData<T>& fineMesh,
    const CrossFieldSolution<T>& cf,
    const CrossFieldMeshAccel<T>& accel,
    int degree,
    T stripCollapseRatio,
    T angleThresholdDeg,
    int maxPatchesCap,
    gsMultiPatch<T>& outMp,
    std::string& methodLabel)
{
    const gsMatrix<T> bbox = dom.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    if (diag <= T(1e-12)) return false;

    T avgEdge = T(0);
    for (const auto& e : fineMesh.edges)
        avgEdge += fineMesh.pts[e.a].dist(fineMesh.pts[e.b]);
    const T targetH = std::max(diag / T(32), avgEdge / std::max<T>(T(1), static_cast<T>(fineMesh.edges.size())));

    int nSharpCorners = 0;
    int nReflexCorners = 0;
    for (size_t i = 0; i < fineMesh.pts.size(); ++i) {
        if (fineMesh.isBoundaryPt[i] && fineMesh.isSharpCorner[i]) {
            ++nSharpCorners;
            if (fineMesh.cornerAngleDeg[i] > angleThresholdDeg) ++nReflexCorners;
        }
    }

    // Fast Path 1A: 0-Hole Convexish 4-Curve Quad -> 1 Patch
    if (dom.numHoles() == 0 && dom.outer().size() == 4 && nReflexCorners == 0) {
        CurveSeg<T> s4[4];
        gsPoint2D<T> c4[4];
        for (int k = 0; k < 4; ++k) {
            s4[k] = CurveSeg<T>(dom.outer().curve(k), degree);
            c4[k] = s4[k].p0();
        }
        if (isStrictlyConvexQuad(c4)) {
            std::unique_ptr<gsTensorBSpline<2, T>> p0(
                makeCoonsFromCCWCycle(s4[0], s4[1], s4[2], s4[3], degree));
            if (p0) {
                if (evalPatchMinDetDense(*p0, 49) <= T(1e-7))
                    quickHarmonicPatchRepair(*p0, T(1e-7));
                T angs[4], mnA = 0, mxA = 180;
                bool degen = false;
                evalPatchCornerAnglesDeg(*p0, angs, mnA, mxA, degen);
                if (!degen && mnA > T(10) && mxA < T(162) && evalPatchMinDetDense(*p0, 49) > T(1e-7)) {
                    outMp.clear();
                    outMp.addPatch(typename gsGeometry<T>::uPtr(p0.release()));
                    methodLabel = "HarmonicCrossField";
                    return true;
                }
            }
        }
    }

    // Fast Path 1B: Smooth 0-Corner Hole-Free Domain -> 5-Patch Butterfly O-Grid
    gsPoint2D<T> polyCtr = polygonCentroid(fineMesh.loopPolys[0]);
    gsPoint2D<T> diskCenter = polyCtr;
    {
        gsPoint2D<T> posSum(0, 0);
        int posCnt = 0;
        for (const auto& sg : cf.singularities) {
            if (sg.index > 0) {
                posSum.x += sg.pos.x * sg.index;
                posSum.y += sg.pos.y * sg.index;
                posCnt += sg.index;
            }
        }
        if (posCnt > 0)
            diskCenter = gsPoint2D<T>(posSum.x / posCnt, posSum.y / posCnt);
    }

    if (dom.numHoles() == 0 && nSharpCorners == 0 && !fineMesh.loopRings.empty()) {
        const auto& ringV = fineMesh.loopRings[0].vertices;
        const int M = static_cast<int>(ringV.size());
        if (M >= 4) {
            T bestMaxAng = T(360);
            gsMultiPatch<T> bestOGrid;
            const int nPhases = (dom.outer().size() <= 4) ? 2 : 6;
            for (int ph = 0; ph < nPhases; ++ph) {
                std::vector<CurveSeg<T>> bEdges;
                if (ph == 0 && dom.outer().size() <= 4 && (dom.outer().size() == 1 || dom.outer().size() == 2 || dom.outer().size() == 4)) {
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
                    } else {
                        for (int i = 0; i < 4; ++i) bEdges.emplace_back(outer.curve(i), degree);
                    }
                } else {
                    int off = (ph * M) / (4 * std::max(1, nPhases));
                    int st[4] = {off % M, (off + M / 4) % M, (off + M / 2) % M, (off + (3 * M) / 4) % M};
                    bEdges.resize(4);
                    for (int k = 0; k < 4; ++k)
                        bEdges[k] = extractExactCADBoundaryBetween(fineMesh, ringV[st[k]], ringV[st[(k + 1) % 4]], degree, true);
                }
                if (bEdges.size() != 4) continue;

                gsPoint2D<T> P[4];
                for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
                gsPoint2D<T> pAvg((P[0].x + P[1].x + P[2].x + P[3].x) * T(0.25),
                                  (P[0].y + P[1].y + P[2].y + P[3].y) * T(0.25));
                std::vector<gsPoint2D<T>> cands = {diskCenter, polyCtr, pAvg};
                const T alphas[4] = {T(0.45), T(0.35), T(0.55), T(0.26)};

                for (const auto& ctr : cands) {
                    if (!isPointStrictlyInsideDomain(ctr, fineMesh.loopPolys, T(1e-5))) continue;
                    for (T alpha : alphas) {
                        gsPoint2D<T> Q[4];
                        bool ok = true;
                        for (int k = 0; k < 4; ++k) {
                            Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                                ctr.y + alpha * (P[k].y - ctr.y));
                            if (!isPointStrictlyInsideDomain(Q[k], fineMesh.loopPolys, T(1e-5)) ||
                                !isSegmentInsideDomain(Q[k], P[k], fineMesh.loopPolys)) {
                                ok = false; break;
                            }
                        }
                        if (!ok || !isStrictlyConvexQuad(Q)) continue;
                        for (int k = 0; k < 4; ++k)
                            if (!isSegmentInsideDomain(Q[k], Q[(k + 1) % 4], fineMesh.loopPolys)) { ok = false; break; }
                        if (!ok) continue;

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
                            if (evalPatchMinDetDense(candMp.patch(p), 49) <= T(1e-7)) {
                                if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&candMp.patch(p)))
                                    quickHarmonicPatchRepair(*tb, T(1e-7));
                            }
                            if (evalPatchMinDetDense(candMp.patch(p), 49) <= T(1e-7)) { allOk = false; break; }
                            T angs[4], mnA = 0, mxA = 180;
                            bool degen = false;
                            evalPatchCornerAnglesDeg(candMp.patch(p), angs, mnA, mxA, degen);
                            if (degen || mnA < T(12) || mxA > T(154)) allOk = false;
                            curMaxA = std::max(curMaxA, mxA);
                        }
                        if (allOk && curMaxA < bestMaxAng) {
                            bestMaxAng = curMaxA;
                            bestOGrid = give(candMp);
                        }
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

    // Skip heavy separatrix arrangement on ultra-complex >= 4-hole domains where
    // O(S^2) separatrix crossings create dense arrangements (Engine 2 handles them cleanly).
    if (dom.numHoles() >= 4)
        return false;

    // Step 1B: Seed Separatrices from Concave Corners & Exact Singularity Ports
    struct RaySeed {
        gsPoint2D<T> start;
        gsPoint2D<T> dir;
        int sourceId;
    };
    std::vector<RaySeed> seeds;
    std::vector<std::pair<int, gsPoint2D<T>>> allSources;

    const T pi = T(3.14159265358979323846);
    const T radToDeg = T(180) / pi;

    for (int i = 0; i < static_cast<int>(fineMesh.pts.size()); ++i) {
        if (!fineMesh.isBoundaryPt[i] || !fineMesh.isSharpCorner[i]) continue;
        int sid = -(i + 1);
        allSources.push_back({sid, fineMesh.pts[i]});
        if (fineMesh.cornerAngleDeg[i] <= angleThresholdDeg) continue;

        gsPoint2D<T> tFwd(fineMesh.bndTanOut[i].x, fineMesh.bndTanOut[i].y);
        const gsPoint2D<T>& tIn = fineMesh.bndTanIn[i];
        T cr = tIn.x * tFwd.y - tIn.y * tFwd.x;
        T dt = tIn.x * tFwd.x + tIn.y * tFwd.y;
        T turnRad = std::atan2(cr, dt);
        T angIfLeft = T(180) - turnRad * radToDeg;
        T sweepSign = (std::abs(fineMesh.cornerAngleDeg[i] - angIfLeft) < T(1)) ? T(1) : T(-1);
        T angStart = std::atan2(tFwd.y, tFwd.x);
        T sweepRad = fineMesh.cornerAngleDeg[i] * pi / T(180);
        int nPorts = std::max(1, static_cast<int>(std::floor(fineMesh.cornerAngleDeg[i] / T(90) + T(0.5))) - 1);
        for (int p = 1; p <= nPorts; ++p) {
            T ang = angStart + sweepSign * sweepRad * (static_cast<T>(p) / static_cast<T>(nPorts + 1));
            seeds.push_back({fineMesh.pts[i], gsPoint2D<T>(std::cos(ang), std::sin(ang)), sid});
        }
    }

    for (size_t s = 0; s < cf.singularities.size(); ++s) {
        const auto& sg = cf.singularities[s];
        int sid = static_cast<int>(s + 1);
        allSources.push_back({sid, sg.pos});
        for (T ang : sg.portAngles)
            seeds.push_back({sg.pos, gsPoint2D<T>(std::cos(ang), std::sin(ang)), sid});
    }
    if (seeds.empty()) return false;

    // Step 1C: RK2 Curved Streamline Tracing
    struct TracedPath {
        std::vector<gsPoint2D<T>> pts;
        int hitBndLoop = -1;
        gsPoint2D<T> hitPt;
        int sourceId = 0;
        int targetSourceId = 0;
        bool alive = true;
    };

    auto pickNearestBranch = [&](const gsPoint2D<T>& q, const gsPoint2D<T>& refDir, int& triHint, gsPoint2D<T>& outDir) -> bool {
        std::complex<T> zq;
        T th = T(0);
        if (!evalCrossFieldAt(fineMesh, cf, accel, q, zq, th, triHint)) return false;
        T bestDot = -T(2);
        gsPoint2D<T> bestD = refDir;
        for (int k = 0; k < 4; ++k) {
            T a = th + k * (pi * T(0.5));
            gsPoint2D<T> cand(std::cos(a), std::sin(a));
            T d = cand.x * refDir.x + cand.y * refDir.y;
            if (d > bestDot) { bestDot = d; bestD = cand; }
        }
        if (bestDot > T(0.35)) { outDir = bestD; return true; }
        return false;
    };

    std::vector<TracedPath> paths;
    const T stepLen = T(0.30) * targetH;
    const int maxSteps = 320;
    const T maxArcLen = T(2.4) * diag;
    const T snapSingTol = T(1.25) * targetH;
    const T snapCornerTol = T(1.05) * targetH;

    for (const auto& sd : seeds) {
        TracedPath tp;
        tp.sourceId = sd.sourceId;
        tp.pts.push_back(sd.start);
        gsPoint2D<T> cur = sd.start;
        gsPoint2D<T> dir = sd.dir;
        bool terminated = false;
        int curTri = -1;
        T arcLen = T(0);

        for (int st = 0; st < maxSteps && !terminated && arcLen < maxArcLen; ++st) {
            if (arcLen >= T(0.9) * targetH) {
                gsPoint2D<T> k1 = dir;
                pickNearestBranch(cur, dir, curTri, k1);
                gsPoint2D<T> midPt(cur.x + T(0.5) * stepLen * k1.x, cur.y + T(0.5) * stepLen * k1.y);
                gsPoint2D<T> k2 = k1;
                if (pickNearestBranch(midPt, k1, curTri, k2)) dir = k2;
                else dir = k1;
            }

            gsPoint2D<T> nxt(cur.x + stepLen * dir.x, cur.y + stepLen * dir.y);

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

            if (arcLen >= T(1.1) * targetH) {
                int bestTarget = 0;
                T bestD = T(1e30);
                gsPoint2D<T> bestPos;
                for (const auto& src : allSources) {
                    if (src.first == sd.sourceId) continue;
                    T tol = (src.first > 0) ? snapSingTol : snapCornerTol;
                    T d = nxt.dist(src.second);
                    if (d < tol && d < bestD && isSegmentInsideDomain(cur, src.second, fineMesh.loopPolys)) {
                        bestD = d; bestTarget = src.first; bestPos = src.second;
                    }
                }
                if (bestTarget != 0) {
                    tp.pts.push_back(bestPos);
                    tp.hitPt = bestPos;
                    tp.targetSourceId = bestTarget;
                    if (bestTarget < 0) {
                        int bvi = -bestTarget - 1;
                        if (bvi >= 0 && bvi < static_cast<int>(fineMesh.bndInfo.size()))
                            tp.hitBndLoop = fineMesh.bndInfo[bvi].loopIdx;
                    }
                    terminated = true;
                    break;
                }
            }

            tp.pts.push_back(nxt);
            arcLen += stepLen;
            cur = nxt;
        }
        if (terminated && tp.pts.size() >= 2)
            paths.push_back(std::move(tp));
    }
    if (paths.empty()) return false;

    // Step 1D: True Anti-Parallel Corridor Merging & Reciprocal Deduplication
    const T corridorTol = std::max(T(1.4) * targetH, stripCollapseRatio * std::min(T(0.09) * diag, T(4.0) * targetH));
    for (size_t i = 0; i < paths.size(); ++i) {
        if (!paths[i].alive || paths[i].targetSourceId == 0) continue;
        for (size_t j = i + 1; j < paths.size(); ++j) {
            if (!paths[j].alive) continue;
            if ((paths[j].sourceId == paths[i].targetSourceId && paths[j].targetSourceId == paths[i].sourceId) ||
                (paths[j].sourceId == paths[i].sourceId && paths[j].targetSourceId == paths[i].targetSourceId)) {
                paths[j].alive = false;
            }
        }
    }

    auto pathUnitDirAt = [](const std::vector<gsPoint2D<T>>& pts, size_t idx) -> gsPoint2D<T> {
        size_t i0 = (idx > 0) ? (idx - 1) : 0;
        size_t i1 = (idx + 1 < pts.size()) ? (idx + 1) : (pts.size() - 1);
        T dx = pts[i1].x - pts[i0].x, dy = pts[i1].y - pts[i0].y;
        T n = std::hypot(dx, dy);
        return (n > T(1e-14)) ? gsPoint2D<T>(dx / n, dy / n) : gsPoint2D<T>(1, 0);
    };

    for (size_t i = 0; i < paths.size(); ++i) {
        if (!paths[i].alive || paths[i].targetSourceId != 0) continue;
        for (size_t j = i + 1; j < paths.size(); ++j) {
            if (!paths[j].alive || paths[j].targetSourceId != 0) continue;
            if (paths[i].sourceId == paths[j].sourceId) continue;

            const gsPoint2D<T>& Si = paths[i].pts.front();
            const gsPoint2D<T>& Sj = paths[j].pts.front();
            size_t bestKi = 0;
            T minDi = T(1e30), minDj = T(1e30);
            for (size_t k = 1; k < paths[i].pts.size(); ++k) {
                T d = paths[i].pts[k].dist(Sj);
                if (d < minDi) { minDi = d; bestKi = k; }
            }
            for (size_t k = 1; k < paths[j].pts.size(); ++k) {
                T d = paths[j].pts[k].dist(Si);
                if (d < minDj) minDj = d;
            }
            if (minDi > corridorTol || minDj > corridorTol) continue;

            gsPoint2D<T> dirI = pathUnitDirAt(paths[i].pts, bestKi);
            gsPoint2D<T> dirJ = pathUnitDirAt(paths[j].pts, 0);
            if (dirI.x * dirJ.x + dirI.y * dirJ.y > -T(0.60)) continue;

            std::vector<gsPoint2D<T>> mergedPts(paths[i].pts.begin(), paths[i].pts.begin() + static_cast<long>(bestKi));
            mergedPts.push_back(Sj);
            if (mergedPts.size() >= 2 && isSegmentInsideDomain(mergedPts[mergedPts.size() - 2], Sj, fineMesh.loopPolys)) {
                paths[i].pts = std::move(mergedPts);
                paths[i].hitPt = Sj;
                paths[i].targetSourceId = paths[j].sourceId;
                paths[i].hitBndLoop = -1;
                if (paths[j].sourceId < 0) {
                    int bvi = -paths[j].sourceId - 1;
                    if (bvi >= 0 && bvi < static_cast<int>(fineMesh.bndInfo.size()))
                        paths[i].hitBndLoop = fineMesh.bndInfo[bvi].loopIdx;
                }
                paths[j].alive = false;
                break;
            }
        }
    }

    // Step 1E: Consistent Exact Boundary Splitting & Curved Polyline Arrangement
    struct ArrVertex { gsPoint2D<T> pt; bool isBoundary = false; int loopIdx = -1; };
    struct ArrEdge { int u = -1, v = -1; CurveSeg<T> seg; bool isBoundary = false; int loopIdx = -1; };
    std::vector<ArrVertex> arrVerts;
    std::vector<ArrEdge> arrEdges;

    const T vMergeTol = T(0.18) * targetH;
    auto addArrVertex = [&](const gsPoint2D<T>& p, bool isBnd, int lIdx, T mergeTol) -> int {
        for (size_t i = 0; i < arrVerts.size(); ++i) {
            if (arrVerts[i].pt.dist(p) <= mergeTol) {
                if (isBnd) {
                    arrVerts[i].isBoundary = true;
                    if (arrVerts[i].loopIdx < 0) arrVerts[i].loopIdx = lIdx;
                }
                return static_cast<int>(i);
            }
        }
        int id = static_cast<int>(arrVerts.size());
        arrVerts.push_back({p, isBnd, lIdx});
        return id;
    };

    const int nLoops = dom.numLoops();
    std::vector<std::vector<CurveSeg<T>>> loopBaseSegs(nLoops);
    std::vector<std::vector<std::vector<T>>> loopSplitUs(nLoops);
    std::vector<int> loopHitCount(nLoops, 0);
    for (int l = 0; l < nLoops; ++l) {
        const auto& lp = dom.loop(l);
        for (int c = 0; c < lp.size(); ++c)
            loopBaseSegs[l].emplace_back(lp.curve(c), degree);
        loopSplitUs[l].resize(lp.size());
    }

    struct LivePolyRec {
        std::vector<gsPoint2D<T>> pts;
        int startLoop = -1, startCrv = -1; T startU = T(0);
        int endLoop = -1, endCrv = -1;     T endU = T(0);
    };
    std::vector<LivePolyRec> livePolys;

    auto projectToLoop = [&](int l, gsPoint2D<T>& pt, int& outCrv, T& outU) {
        outCrv = -1; outU = T(0);
        if (l < 0 || l >= nLoops) return;
        int bestC = -1;
        T bestU = T(0.5), bestD = std::numeric_limits<T>::max();
        for (size_t c = 0; c < loopBaseSegs[l].size(); ++c) {
            T u = loopBaseSegs[l][c].closestParam(pt);
            gsPoint2D<T> q = evalCurvePoint(*loopBaseSegs[l][c].bsp, u);
            T d = q.dist(pt);
            if (d < bestD) { bestD = d; bestC = static_cast<int>(c); bestU = u; }
        }
        if (bestC >= 0) {
            ++loopHitCount[l];
            outCrv = bestC;
            if (bestU < T(0.06)) {
                outU = T(0);
                pt = loopBaseSegs[l][bestC].p0();
            } else if (bestU > T(0.94)) {
                outU = T(1);
                pt = loopBaseSegs[l][bestC].p1();
            } else {
                outU = bestU;
                pt = evalCurvePoint(*loopBaseSegs[l][bestC].bsp, bestU);
                loopSplitUs[l][bestC].push_back(bestU);
            }
        }
    };

    for (auto& path : paths) {
        if (!path.alive || path.pts.size() < 2) continue;
        LivePolyRec rec;
        rec.pts = path.pts;
        if (path.sourceId < 0) {
            int bvi = -path.sourceId - 1;
            if (bvi >= 0 && bvi < static_cast<int>(fineMesh.bndInfo.size())) {
                rec.startLoop = fineMesh.bndInfo[bvi].loopIdx;
                projectToLoop(rec.startLoop, rec.pts.front(), rec.startCrv, rec.startU);
            }
        }
        if (path.hitBndLoop >= 0) {
            rec.endLoop = path.hitBndLoop;
            projectToLoop(rec.endLoop, rec.pts.back(), rec.endCrv, rec.endU);
        }
        if (rec.pts.front().dist(rec.pts.back()) > vMergeTol)
            livePolys.push_back(std::move(rec));
    }

    // Bridge any unhit hole loop cleanly
    for (int l = 1; l < nLoops; ++l) {
        if (loopHitCount[l] >= 2) continue;
        const auto& hPoly = fineMesh.loopPolys[l];
        if (hPoly.size() < 4) continue;
        gsPoint2D<T> hCtr = polygonCentroid(hPoly);
        int needed = 2 - loopHitCount[l];
        std::vector<gsPoint2D<T>> usedDirs;
        for (int b = 0; b < needed; ++b) {
            T bestDist = std::numeric_limits<T>::max();
            gsPoint2D<T> bestA, bestB;
            int bestOtherLoop = -1;
            const int nH = static_cast<int>(hPoly.size());
            const int strideH = std::max(1, nH / 16);
            for (int hi = 0; hi < nH; hi += strideH) {
                const auto& pH = hPoly[hi];
                gsPoint2D<T> rH(pH.x - hCtr.x, pH.y - hCtr.y);
                T nR = std::hypot(rH.x, rH.y);
                if (nR > T(1e-12)) { rH.x /= nR; rH.y /= nR; }
                bool tooAligned = false;
                for (const auto& ud : usedDirs)
                    if (ud.x * rH.x + ud.y * rH.y > T(0.35)) { tooAligned = true; break; }
                if (tooAligned) continue;

                for (int ol = 0; ol < nLoops; ++ol) {
                    if (ol == l) continue;
                    const auto& oPoly = fineMesh.loopPolys[ol];
                    const int nO = static_cast<int>(oPoly.size());
                    const int strideO = std::max(1, nO / 24);
                    for (int oi = 0; oi < nO; oi += strideO) {
                        T d = pH.dist(oPoly[oi]);
                        if (d < bestDist && d > vMergeTol && isSegmentInsideDomain(pH, oPoly[oi], fineMesh.loopPolys)) {
                            bestDist = d; bestA = pH; bestB = oPoly[oi]; bestOtherLoop = ol;
                        }
                    }
                }
            }
            if (bestOtherLoop >= 0) {
                LivePolyRec rec;
                rec.startLoop = l; rec.endLoop = bestOtherLoop;
                projectToLoop(l, bestA, rec.startCrv, rec.startU);
                projectToLoop(bestOtherLoop, bestB, rec.endCrv, rec.endU);
                if (isSegmentInsideDomain(bestA, bestB, fineMesh.loopPolys)) {
                    if (rec.startU > T(0.03) && rec.startU < T(0.97))
                        loopSplitUs[l][rec.startCrv].push_back(rec.startU);
                    if (rec.endU > T(0.03) && rec.endU < T(0.97))
                        loopSplitUs[bestOtherLoop][rec.endCrv].push_back(rec.endU);
                    rec.pts = {bestA, bestB};
                    livePolys.push_back(std::move(rec));
                    gsPoint2D<T> ud(bestA.x - hCtr.x, bestA.y - hCtr.y);
                    T nU = std::hypot(ud.x, ud.y);
                    if (nU > T(1e-12)) usedDirs.emplace_back(ud.x / nU, ud.y / nU);
                }
            }
        }
    }

    if (livePolys.empty()) return false;

    for (int l = 0; l < nLoops; ++l) {
        size_t totalVerts = loopBaseSegs[l].size();
        for (const auto& us : loopSplitUs[l]) totalVerts += us.size();
        if (totalVerts < 3 && !loopBaseSegs[l].empty()) {
            if (loopBaseSegs[l].size() == 1) {
                loopSplitUs[l][0].push_back(T(0.25));
                loopSplitUs[l][0].push_back(T(0.50));
                loopSplitUs[l][0].push_back(T(0.75));
            } else if (loopBaseSegs[l].size() == 2) {
                loopSplitUs[l][0].push_back(T(0.50));
                loopSplitUs[l][1].push_back(T(0.50));
            }
        }
    }

    // Deduplicate split parameters and snap EVERY livePoly boundary endpoint to the EXACT kept split vertex!
    std::vector<std::vector<std::vector<T>>> cleanSplitUs(nLoops);
    for (int l = 0; l < nLoops; ++l) {
        cleanSplitUs[l].resize(loopBaseSegs[l].size());
        for (size_t c = 0; c < loopBaseSegs[l].size(); ++c) {
            auto& us = loopSplitUs[l][c];
            std::sort(us.begin(), us.end());
            auto& cl = cleanSplitUs[l][c];
            for (T u : us) {
                if (u > T(0.06) && u < T(0.94) && (cl.empty() || u - cl.back() > T(0.06)))
                    cl.push_back(u);
            }
            std::vector<CurveSeg<T>> pieces = splitSegMulti(loopBaseSegs[l][c], cl);
            for (const auto& seg : pieces) {
                int va = addArrVertex(seg.p0(), true, l, T(1e-5) * diag);
                int vb = addArrVertex(seg.p1(), true, l, T(1e-5) * diag);
                if (va != vb)
                    arrEdges.push_back({va, vb, seg, true, l});
            }
        }
    }

    auto snapEndpointToCleanSplit = [&](int l, int c, T u, gsPoint2D<T>& pt) {
        if (l < 0 || l >= nLoops || c < 0 || c >= static_cast<int>(loopBaseSegs[l].size())) return;
        T bestU = T(0), bestDiff = std::abs(u - T(0));
        if (std::abs(u - T(1)) < bestDiff) { bestU = T(1); bestDiff = std::abs(u - T(1)); }
        for (T cu : cleanSplitUs[l][c]) {
            if (std::abs(u - cu) < bestDiff) { bestU = cu; bestDiff = std::abs(u - cu); }
        }
        pt = (bestU <= T(1e-12)) ? loopBaseSegs[l][c].p0()
           : ((bestU >= T(1) - T(1e-12)) ? loopBaseSegs[l][c].p1()
              : evalCurvePoint(*loopBaseSegs[l][c].bsp, bestU));
    };
    for (auto& rec : livePolys) {
        if (rec.startLoop >= 0) snapEndpointToCleanSplit(rec.startLoop, rec.startCrv, rec.startU, rec.pts.front());
        if (rec.endLoop >= 0)   snapEndpointToCleanSplit(rec.endLoop,   rec.endCrv,   rec.endU,   rec.pts.back());
    }

    // Intersect actual curved polylines
    const size_t nP = livePolys.size();
    std::vector<std::vector<T>> polyCumS(nP);
    for (size_t i = 0; i < nP; ++i) {
        polyCumS[i].resize(livePolys[i].pts.size(), T(0));
        for (size_t k = 1; k < livePolys[i].pts.size(); ++k)
            polyCumS[i][k] = polyCumS[i][k - 1] + livePolys[i].pts[k - 1].dist(livePolys[i].pts[k]);
    }

    std::vector<std::vector<T>> polySplits(nP);
    for (size_t i = 0; i < nP; ++i) {
        const auto& Pi = livePolys[i].pts;
        const T Li = polyCumS[i].back();
        if (Li <= T(1e-12)) continue;
        for (size_t j = i + 1; j < nP; ++j) {
            const auto& Pj = livePolys[j].pts;
            const T Lj = polyCumS[j].back();
            if (Lj <= T(1e-12)) continue;
            for (size_t ki = 0; ki + 1 < Pi.size(); ++ki) {
                const auto& A = Pi[ki];
                const auto& B = Pi[ki + 1];
                T rx = B.x - A.x, ry = B.y - A.y;
                T lenI = std::hypot(rx, ry);
                if (lenI <= T(1e-12)) continue;
                for (size_t kj = 0; kj + 1 < Pj.size(); ++kj) {
                    const auto& C = Pj[kj];
                    const auto& D = Pj[kj + 1];
                    T sx = D.x - C.x, sy = D.y - C.y;
                    T lenJ = std::hypot(sx, sy);
                    if (lenJ <= T(1e-12)) continue;
                    T denom = rx * sy - ry * sx;
                    if (std::abs(denom) < T(0.22) * lenI * lenJ) continue;
                    T t = ((C.x - A.x) * sy - (C.y - A.y) * sx) / denom;
                    T u = ((C.x - A.x) * ry - (C.y - A.y) * rx) / denom;
                    if (t >= T(0) && t <= T(1) && u >= T(0) && u <= T(1)) {
                        T si = (polyCumS[i][ki] + t * lenI) / Li;
                        T sj = (polyCumS[j][kj] + u * lenJ) / Lj;
                        if (si > T(0.05) && si < T(0.95) && sj > T(0.05) && sj < T(0.95)) {
                            polySplits[i].push_back(si);
                            polySplits[j].push_back(sj);
                        }
                    }
                }
            }
        }
    }

    auto extractSubPolyline = [&](const std::vector<gsPoint2D<T>>& pts, const std::vector<T>& cumS,
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
    for (const auto& ae : arrEdges)
        existingPairs.insert({std::min(ae.u, ae.v), std::max(ae.u, ae.v)});

    for (size_t i = 0; i < nP; ++i) {
        auto& ts = polySplits[i];
        ts.push_back(T(0));
        ts.push_back(T(1));
        std::sort(ts.begin(), ts.end());
        std::vector<T> uniq;
        for (T t : ts)
            if (uniq.empty() || t - uniq.back() > T(0.025)) uniq.push_back(t);
        if (uniq.size() < 2) uniq = {T(0), T(1)};
        uniq.front() = T(0);
        uniq.back() = T(1);

        for (size_t k = 0; k + 1 < uniq.size(); ++k) {
            std::vector<gsPoint2D<T>> sub = extractSubPolyline(livePolys[i].pts, polyCumS[i], uniq[k], uniq[k + 1]);
            int va = addArrVertex(sub.front(), false, -1, vMergeTol);
            int vb = addArrVertex(sub.back(), false, -1, vMergeTol);
            if (va == vb) continue;
            auto key = std::make_pair(std::min(va, vb), std::max(va, vb));
            if (existingPairs.count(key)) continue;
            existingPairs.insert(key);
            sub.front() = arrVerts[va].pt;
            sub.back() = arrVerts[vb].pt;
            CurveSeg<T> fitted = fitStreamlinePolylineToCurveSeg(sub, degree, fineMesh.loopPolys);
            arrEdges.push_back({va, vb, fitted, false, -1});
        }
    }

    // Prune degree-1 interior rays
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

    // Dissolve smooth valence-2 boundary vertices (preserving >= 3 vertices per loop)
    bool dissolved = true;
    while (dissolved) {
        dissolved = false;
        std::vector<std::vector<int>> vInc(arrVerts.size());
        std::vector<int> loopEdgeCnt(nLoops, 0);
        for (int ei = 0; ei < static_cast<int>(arrEdges.size()); ++ei) {
            vInc[arrEdges[ei].u].push_back(ei);
            vInc[arrEdges[ei].v].push_back(ei);
            if (arrEdges[ei].isBoundary && arrEdges[ei].loopIdx >= 0 && arrEdges[ei].loopIdx < nLoops)
                ++loopEdgeCnt[arrEdges[ei].loopIdx];
        }
        for (int vi = 0; vi < static_cast<int>(arrVerts.size()); ++vi) {
            if (!arrVerts[vi].isBoundary || vInc[vi].size() != 2) continue;
            int e0 = vInc[vi][0], e1 = vInc[vi][1];
            if (!arrEdges[e0].isBoundary || !arrEdges[e1].isBoundary) continue;
            int lIdx = arrEdges[e0].loopIdx;
            if (lIdx < 0 || lIdx != arrEdges[e1].loopIdx || loopEdgeCnt[lIdx] <= 3) continue;

            // Ensure consistent orientation vA -> vi -> vB along the boundary loop
            if (arrEdges[e0].v != vi && arrEdges[e1].v == vi) std::swap(e0, e1);
            if (arrEdges[e0].v != vi || arrEdges[e1].u != vi) continue;
            int vA = arrEdges[e0].u, vB = arrEdges[e1].v;
            if (vA == vB) continue;

            auto dIn = evalCurveDeriv(*arrEdges[e0].seg.bsp, T(1));
            auto dOut = evalCurveDeriv(*arrEdges[e1].seg.bsp, T(0));
            T nIn = std::hypot(dIn.x, dIn.y), nOut = std::hypot(dOut.x, dOut.y);
            if (nIn <= T(1e-12) || nOut <= T(1e-12)) continue;
            if ((dIn.x * dOut.x + dIn.y * dOut.y) / (nIn * nOut) < T(0.99)) continue;

            auto joinedPtr = joinCurveSegsExact(std::vector<CurveSeg<T>>{arrEdges[e0].seg, arrEdges[e1].seg}, degree);
            if (!joinedPtr) continue;
            arrEdges[e0] = {vA, vB, CurveSeg<T>(*joinedPtr), true, lIdx};
            arrEdges.erase(arrEdges.begin() + e1);
            dissolved = true;
            break;
        }
    }

    if (arrEdges.size() < 4) return false;

    struct HalfEdge {
        int origin = -1, target = -1, arrEdgeIdx = -1;
        bool forward = true;
        int twin = -1, next = -1;
        T angle = T(0);
    };

    auto extractFaces = [&](std::vector<HalfEdge>& halfEdges, std::vector<std::vector<int>>& faces) -> bool {
        const int nE = static_cast<int>(arrEdges.size());
        halfEdges.assign(2 * nE, HalfEdge());
        std::vector<std::vector<int>> outHalf(arrVerts.size());

        for (int i = 0; i < nE; ++i) {
            const auto& e = arrEdges[i];
            int h0 = 2 * i, h1 = 2 * i + 1;
            gsPoint2D<T> d0 = evalCurveDeriv(*e.seg.bsp, T(0));
            gsPoint2D<T> d1 = evalCurveDeriv(*e.seg.bsp, T(1));
            T a0 = (std::hypot(d0.x, d0.y) > T(1e-12))
                 ? std::atan2(d0.y, d0.x)
                 : std::atan2(arrVerts[e.v].pt.y - arrVerts[e.u].pt.y, arrVerts[e.v].pt.x - arrVerts[e.u].pt.x);
            T a1 = (std::hypot(d1.x, d1.y) > T(1e-12))
                 ? std::atan2(-d1.y, -d1.x)
                 : std::atan2(arrVerts[e.u].pt.y - arrVerts[e.v].pt.y, arrVerts[e.u].pt.x - arrVerts[e.v].pt.x);
            halfEdges[h0] = {e.u, e.v, i, true,  h1, -1, a0};
            halfEdges[h1] = {e.v, e.u, i, false, h0, -1, a1};
            outHalf[e.u].push_back(h0);
            outHalf[e.v].push_back(h1);
        }

        for (size_t v = 0; v < arrVerts.size(); ++v) {
            auto& oh = outHalf[v];
            std::sort(oh.begin(), oh.end(), [&](int a, int b) { return halfEdges[a].angle < halfEdges[b].angle; });
        }

        for (int h = 0; h < 2 * nE; ++h) {
            int v = halfEdges[h].target;
            int tw = halfEdges[h].twin;
            const auto& oh = outHalf[v];
            auto it = std::find(oh.begin(), oh.end(), tw);
            if (it == oh.end() || oh.empty()) return false;
            size_t pos = static_cast<size_t>(it - oh.begin());
            halfEdges[h].next = oh[(pos + oh.size() - 1) % oh.size()];
        }

        std::vector<bool> heVisited(2 * nE, false);
        faces.clear();
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

            std::vector<gsPoint2D<T>> poly;
            poly.reserve(cycle.size() * 4);
            for (int ch : cycle) {
                const auto& he = halfEdges[ch];
                const auto& ae = arrEdges[he.arrEdgeIdx];
                for (int s = 0; s < 4; ++s) {
                    T u = he.forward ? (static_cast<T>(s) / T(4)) : (T(1) - static_cast<T>(s) / T(4));
                    poly.push_back(evalCurvePoint(*ae.seg.bsp, u));
                }
            }
            T area2 = T(0);
            for (size_t k = 0; k < poly.size(); ++k) {
                const auto& p = poly[k];
                const auto& q = poly[(k + 1) % poly.size()];
                area2 += p.x * q.y - q.x * p.y;
            }
            if (area2 <= T(1e-7)) continue;
            gsPoint2D<T> ctr = polygonCentroid(poly);
            if (!isPointStrictlyInsideDomain(ctr, fineMesh.loopPolys, T(1e-6))) continue;
            faces.push_back(std::move(cycle));
        }
        return !faces.empty();
    };

    std::vector<HalfEdge> halfEdges;
    std::vector<std::vector<int>> faces;
    if (!extractFaces(halfEdges, faces)) return false;

    // Merge adjacent odd-sided faces (3+3 -> 4, 3+5 -> 6, 5+5 -> 8)
    for (int mPass = 0; mPass < 16; ++mPass) {
        std::vector<int> heToFace(halfEdges.size(), -1);
        for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi)
            for (int ch : faces[fi]) heToFace[ch] = fi;

        int killEdge = -1;
        for (int fi = 0; fi < static_cast<int>(faces.size()) && killEdge < 0; ++fi) {
            if ((faces[fi].size() % 2) == 0) continue;
            for (int ch : faces[fi]) {
                int ei = halfEdges[ch].arrEdgeIdx;
                if (arrEdges[ei].isBoundary) continue;
                int of = heToFace[halfEdges[ch].twin];
                if (of >= 0 && of != fi && (faces[of].size() % 2) == 1) {
                    killEdge = ei; break;
                }
            }
        }
        if (killEdge < 0) break;
        arrEdges.erase(arrEdges.begin() + killEdge);
        if (!extractFaces(halfEdges, faces)) return false;
    }

    auto getArrSeg = [&](int heIdx) -> CurveSeg<T> {
        const auto& he = halfEdges[heIdx];
        const auto& ae = arrEdges[he.arrEdgeIdx];
        return he.forward ? ae.seg : ae.seg.reversed();
    };

    // Zero-Catmull-Clark emission if all faces are even (4-gons direct, 2k-gons locally chord-split!)
    bool allEven = true;
    for (const auto& f : faces)
        if ((f.size() % 2) != 0) { allEven = false; break; }

    if (allEven) {
        std::vector<std::array<CurveSeg<T>, 4>> quadList;
        bool splitOk = true;
        for (const auto& f : faces) {
            if (f.size() == 4) {
                quadList.push_back({getArrSeg(f[0]), getArrSeg(f[1]), getArrSeg(f[2]), getArrSeg(f[3])});
            } else {
                std::vector<CurveSeg<T>> fEdges;
                std::vector<gsPoint2D<T>> fVerts;
                for (int ch : f) {
                    fEdges.push_back(getArrSeg(ch));
                    fVerts.push_back(arrVerts[halfEdges[ch].origin].pt);
                }
                if (!splitEvenFaceToQuadsRecursive(fEdges, fVerts, fineMesh, cf, accel, degree, quadList)) {
                    splitOk = false; break;
                }
            }
        }
        if (splitOk && !quadList.empty() && static_cast<int>(quadList.size()) <= maxPatchesCap) {
            outMp.clear();
            for (const auto& q : quadList)
                outMp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(q[0], q[1], q[2], q[3], degree)));
            methodLabel = "HarmonicCrossField";
            return true;
        }
    }

    // Bounded midpoint-centroid fallback ONLY if compact (<= min(maxPatchesCap, 48))
    int ccCount = 0;
    for (const auto& f : faces) ccCount += static_cast<int>(f.size());
    if (ccCount == 0 || ccCount > std::min(maxPatchesCap, 48)) return false;

    const int nE = static_cast<int>(arrEdges.size());
    struct SplitSegs { CurveSeg<T> h0, h1; gsPoint2D<T> mid; };
    std::vector<SplitSegs> edgeHalves(nE);
    for (int i = 0; i < nE; ++i) {
        auto sp = arrEdges[i].seg.splitArc(T(0.5));
        edgeHalves[i] = {sp.first, sp.second, sp.first.p1()};
    }

    auto getArrHalves = [&](int heIdx, CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
        const auto& he = halfEdges[heIdx];
        const auto& eh = edgeHalves[he.arrEdgeIdx];
        outMid = eh.mid;
        if (he.forward) { out0 = eh.h0; out1 = eh.h1; }
        else            { out0 = eh.h1.reversed(); out1 = eh.h0.reversed(); }
    };

    outMp.clear();
    for (const auto& f : faces) {
        int nf = static_cast<int>(f.size());
        std::vector<CurveSeg<T>> h0(nf), h1(nf);
        std::vector<gsPoint2D<T>> mids(nf), vertsF(nf);
        for (int k = 0; k < nf; ++k) {
            vertsF[k] = arrVerts[halfEdges[f[k]].origin].pt;
            getArrHalves(f[k], h0[k], h1[k], mids[k]);
        }
        gsPoint2D<T> O = polygonCentroid(mids);
        if (nf == 3) {
            gsPoint2D<T> optO;
            if (pickOptimalTri1to3CenterCurved(
                    vertsF[0], vertsF[1], vertsF[2],
                    h0[0], h1[0], h0[1], h1[1], h0[2], h1[2],
                    arrEdges[halfEdges[f[0]].arrEdgeIdx].isBoundary,
                    arrEdges[halfEdges[f[1]].arrEdgeIdx].isBoundary,
                    arrEdges[halfEdges[f[2]].arrEdgeIdx].isBoundary, optO))
                O = optO;
        } else if (!isPointStrictlyInsideDomain(O, fineMesh.loopPolys, T(1e-6))) {
            O = polygonCentroid(vertsF);
        }
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

// ============================================================================
// Engine 2: MBO Cross-Field-Aligned Coarse Blossom-Quad Decomposition
// ============================================================================

template <class T>
inline void flipEdgesForCrossField(
    TriMeshData<T>& mesh,
    const TriMeshData<T>& fineMesh,
    const CrossFieldSolution<T>& cf,
    const CrossFieldMeshAccel<T>& accel,
    T minPairAng,
    T maxPairAng)
{
    const int nT = static_cast<int>(mesh.tris.size());
    const int nE = static_cast<int>(mesh.edges.size());
    if (nT < 2 || nE == 0) return;

    int triHint = -1;
    auto quadCFAlign = [&](const int qv[4], int diagA, int diagB) -> T {
        T sum = T(0);
        for (int k = 0; k < 4; ++k)
            sum += evalSegmentCrossFieldAlignment(fineMesh, cf, accel, mesh.pts[qv[k]], mesh.pts[qv[(k + 1) % 4]], triHint);
        T dAl = evalSegmentCrossFieldAlignment(fineMesh, cf, accel, mesh.pts[diagA], mesh.pts[diagB], triHint);
        return T(0.25) * sum - T(0.35) * dAl;
    };

    for (int pass = 0; pass < 3; ++pass) {
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

        auto edgePairScore = [&](int tA, int kA, int tB, int kB) -> T {
            const auto& A = mesh.tris[tA];
            const auto& B = mesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) return T(-1e6);
            int qv[4] = {c, a, d, b};
            int qe[4] = {A.e[(kA + 2) % 3], B.e[(kB + 1) % 3], B.e[(kB + 2) % 3], A.e[(kA + 1) % 3]};
            T eff = 0, sc = 0;
            if (!evalMeshQuadUnion(mesh, qv, qe, minPairAng, maxPairAng, eff, sc, false))
                return T(0);
            return eff + T(7.5) * quadCFAlign(qv, a, b);
        };

        std::vector<T> triBestSc(nT, T(0));
        for (int eIdx = 0; eIdx < nE; ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = eAdj[eIdx].t0, tB = eAdj[eIdx].t1;
            if (tA < 0 || tB < 0) continue;
            T sc = edgePairScore(tA, eAdj[eIdx].k0, tB, eAdj[eIdx].k1);
            triBestSc[tA] = std::max(triBestSc[tA], sc);
            triBestSc[tB] = std::max(triBestSc[tB], sc);
        }

        bool anyFlip = false;
        for (int eIdx = 0; eIdx < nE; ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
            int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
            if (tA < 0 || tB < 0) continue;

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

            T curBest = std::max(triBestSc[tA], triBestSc[tB]);
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
                if (ot < 0 || eid == eIdx) {
                    if (eid == eIdx) newBest = std::max(newBest, edgePairScore(tA, 2, tB, 0));
                    return;
                }
                newBest = std::max(newBest, edgePairScore(triIdx, slotIdx, ot, ok));
            };
            testEdgeSlot(eIdx, tA, 2);
            testEdgeSlot(e_ca, tA, 0);
            testEdgeSlot(e_ad, tA, 1);
            testEdgeSlot(e_db, tB, 1);
            testEdgeSlot(e_bc, tB, 2);

            if ((curBest <= T(1e-6) && newBest >= T(12)) || (newBest > curBest + T(5.0))) {
                eAdj[eIdx] = {tA, 2, tB, 0};
                if (eAdj[e_ca].t0 == tA) eAdj[e_ca].k0 = 0; else eAdj[e_ca].k1 = 0;
                if (eAdj[e_ad].t0 == tB) { eAdj[e_ad].t0 = tA; eAdj[e_ad].k0 = 1; } else { eAdj[e_ad].t1 = tA; eAdj[e_ad].k1 = 1; }
                if (eAdj[e_db].t0 == tB) eAdj[e_db].k0 = 1; else eAdj[e_db].k1 = 1;
                if (eAdj[e_bc].t0 == tA) { eAdj[e_bc].t0 = tB; eAdj[e_bc].k0 = 2; } else { eAdj[e_bc].t1 = tB; eAdj[e_bc].k1 = 2; }
                --vertDeg[a]; --vertDeg[b]; ++vertDeg[c]; ++vertDeg[d];
                triBestSc[tA] = std::max(triBestSc[tA], newBest);
                triBestSc[tB] = std::max(triBestSc[tB], newBest);
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

template <class T>
inline void injectCrossFieldSingularities(
    TriMeshData<T>& mesh,
    const CrossFieldSolution<T>& cf,
    T minClearDist,
    int maxInject = 10)
{
    int injected = 0;
    for (const auto& sg : cf.singularities) {
        if (injected >= maxInject) break;
        if (!isPointStrictlyInsideDomain(sg.pos, mesh.loopPolys, T(1e-4))) continue;
        bool tooClose = false;
        for (const auto& pt : mesh.pts)
            if (pt.dist(sg.pos) < minClearDist) { tooClose = true; break; }
        if (tooClose) continue;

        int hitTi = -1;
        for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
            const auto& t = mesh.tris[ti];
            const auto& A = mesh.pts[t.v[0]];
            const auto& B = mesh.pts[t.v[1]];
            const auto& C = mesh.pts[t.v[2]];
            if (cross2(A, B, C) <= T(1e-10)) continue;
            if (cross2(A, B, sg.pos) > T(1e-7) &&
                cross2(B, C, sg.pos) > T(1e-7) &&
                cross2(C, A, sg.pos) > T(1e-7)) {
                T mn0, mx0, mn1, mx1, mn2, mx2;
                triMinMaxAnglesDeg(A, B, sg.pos, mn0, mx0);
                triMinMaxAnglesDeg(B, C, sg.pos, mn1, mx1);
                triMinMaxAnglesDeg(C, A, sg.pos, mn2, mx2);
                if (std::min(mn0, std::min(mn1, mn2)) >= T(10.0) &&
                    isSegmentInsideDomain(A, sg.pos, mesh.loopPolys) &&
                    isSegmentInsideDomain(B, sg.pos, mesh.loopPolys) &&
                    isSegmentInsideDomain(C, sg.pos, mesh.loopPolys) &&
                    curvedBoundaryEdgeClearsVertex(mesh, t.e[0], t.v[0], t.v[1], sg.pos) &&
                    curvedBoundaryEdgeClearsVertex(mesh, t.e[1], t.v[1], t.v[2], sg.pos) &&
                    curvedBoundaryEdgeClearsVertex(mesh, t.e[2], t.v[2], t.v[0], sg.pos)) {
                    hitTi = ti; break;
                }
            }
        }
        if (hitTi < 0) continue;

        const auto t = mesh.tris[hitTi];
        int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
        int vc = static_cast<int>(mesh.pts.size());
        mesh.pts.push_back(sg.pos);
        mesh.isBoundaryPt.push_back(false);
        mesh.isSmoothBnd.push_back(false);
        mesh.isSharpCorner.push_back(false);
        mesh.cornerAngleDeg.push_back(T(360));
        mesh.bndTanIn.push_back(gsPoint2D<T>(1, 0));
        mesh.bndTanOut.push_back(gsPoint2D<T>(1, 0));

        int e0c = static_cast<int>(mesh.edges.size());
        mesh.edges.push_back({v0, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
        int e1c = static_cast<int>(mesh.edges.size());
        mesh.edges.push_back({v1, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
        int e2c = static_cast<int>(mesh.edges.size());
        mesh.edges.push_back({v2, vc, CurveSeg<T>(), CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});

        mesh.tris[hitTi] = {{v0, v1, vc}, {t.e[0], e1c, e0c}};
        mesh.tris.push_back({{v1, v2, vc}, {t.e[1], e2c, e1c}});
        mesh.tris.push_back({{v2, v0, vc}, {t.e[2], e0c, e2c}});
        ++injected;
    }
    if (injected > 0)
        applyDomainLawsonFlips(mesh);
}

// ============================================================================
// Cross-Field Multi-Patch Geometry Conditioning Helpers:
//   1. evalMultiPatchMaxEdgeRatio  : worst corner-chord edge ratio across patches
//   2. relaxInteriorVertices       : conformal Greville-weighted 1-ring relaxation
//      of strictly interior shared vertices (e.g., ellipse O-grid & pinched hubs)
//   3. relaxBadCornerAngles        : conformal interior-interface tangent-handle
//      rotation to rescue near-degenerate corners (< 24 deg or > 150 deg)
// ============================================================================

template <class T>
inline T evalMultiPatchMaxEdgeRatio(const gsMultiPatch<T>& mp)
{
    T maxR = T(1);
    gsMatrix<T> u(2, 4), ev;
    u << T(0), T(1), T(1), T(0),
         T(0), T(0), T(1), T(1);
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        mp.patch(p).eval_into(u, ev);
        T mn = std::numeric_limits<T>::max(), mx = T(0);
        for (int k = 0; k < 4; ++k) {
            int kn = (k + 1) % 4;
            T L = std::hypot(ev(0, kn) - ev(0, k), ev(1, kn) - ev(1, k));
            mn = std::min(mn, L);
            mx = std::max(mx, L);
        }
        if (mn > T(1e-12))
            maxR = std::max(maxR, mx / mn);
        else
            maxR = std::max(maxR, T(1e3));
    }
    return maxR;
}

template <class T>
inline std::vector<T> computeNormalizedGreville(const gsKnotVector<T>& kv)
{
    gsMatrix<T> g = kv.greville();
    const int n = static_cast<int>(g.cols());
    std::vector<T> out(n, T(0));
    if (n <= 1) return out;
    T a = g(0, 0), b = g(0, n - 1), span = b - a;
    if (std::abs(span) <= T(1e-14)) {
        for (int i = 0; i < n; ++i) out[i] = static_cast<T>(i) / static_cast<T>(n - 1);
        return out;
    }
    for (int i = 0; i < n; ++i)
        out[i] = std::max(T(0), std::min(T(1), (g(0, i) - a) / span));
    return out;
}

template <class T>
inline bool relaxInteriorVertices(gsMultiPatch<T>& mp, int numSweeps = 5)
{
    const int nP = static_cast<int>(mp.nPatches());
    if (nP < 3) return false;

    mp.computeTopology(T(1e-4), false, false);
    if (2 * mp.nInterfaces() + mp.nBoundary() != 4 * mp.nPatches()) return false;

    std::vector<gsTensorBSpline<2, T>*> tbs(nP, nullptr);
    std::vector<int> nu(nP, 0), nv(nP, 0);
    std::vector<std::vector<T>> gu(nP), gv(nP);
    for (int p = 0; p < nP; ++p) {
        tbs[p] = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tbs[p]) return false;
        nu[p] = tbs[p]->basis().size(0);
        nv[p] = tbs[p]->basis().size(1);
        gu[p] = computeNormalizedGreville(tbs[p]->knots(0));
        gv[p] = computeNormalizedGreville(tbs[p]->knots(1));
    }

    auto cornerPos = [&](int p, int ci, int cj) -> gsPoint2D<T> {
        int i = ci ? (nu[p] - 1) : 0;
        int j = cj ? (nv[p] - 1) : 0;
        const auto& cf = tbs[p]->coefs();
        return gsPoint2D<T>(cf(i + j * nu[p], 0), cf(i + j * nu[p], 1));
    };

    // Cluster corners into unique mesh vertices
    struct IncCorner { int p; int ci; int cj; };
    std::vector<gsPoint2D<T>> verts;
    std::vector<std::vector<IncCorner>> vInc;
    std::vector<bool> vFixed;
    std::vector<std::array<int, 4>> pVertId(nP); // (0,0)=0, (1,0)=1, (1,1)=2, (0,1)=3

    auto findOrAddV = [&](const gsPoint2D<T>& q) -> int {
        for (size_t i = 0; i < verts.size(); ++i)
            if (verts[i].dist(q) < T(1e-4)) return static_cast<int>(i);
        int id = static_cast<int>(verts.size());
        verts.push_back(q);
        vInc.emplace_back();
        vFixed.push_back(false);
        return id;
    };

    const int cIdx[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (int p = 0; p < nP; ++p) {
        for (int c = 0; c < 4; ++c) {
            int vid = findOrAddV(cornerPos(p, cIdx[c][0], cIdx[c][1]));
            pVertId[p][c] = vid;
            vInc[vid].push_back({p, cIdx[c][0], cIdx[c][1]});
        }
        // Mark vertices touching any domain boundary side as fixed
        for (int s = 1; s <= 4; ++s) {
            if (!mp.isInterface(patchSide(p, static_cast<boxSide>(s)))) {
                if (s == 1) { vFixed[pVertId[p][0]] = vFixed[pVertId[p][3]] = true; }
                if (s == 2) { vFixed[pVertId[p][1]] = vFixed[pVertId[p][2]] = true; }
                if (s == 3) { vFixed[pVertId[p][0]] = vFixed[pVertId[p][1]] = true; }
                if (s == 4) { vFixed[pVertId[p][3]] = vFixed[pVertId[p][2]] = true; }
            }
        }
    }

    // Also fix endpoints of any interface whose Greville abscissae do not match in the ACTUAL
    // endpoint orientation (preserves exact C0 conformity on non-symmetric knot vectors!)
    auto sideCPPos = [&](int p, int s, int k) -> gsPoint2D<T> {
        int loc = (s == 1) ? (0 + k * nu[p])
                : ((s == 2) ? ((nu[p] - 1) + k * nu[p])
                   : ((s == 3) ? (k + 0 * nu[p]) : (k + (nv[p] - 1) * nu[p])));
        return gsPoint2D<T>(tbs[p]->coefs()(loc, 0), tbs[p]->coefs()(loc, 1));
    };
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
        int p1 = static_cast<int>(it->first().patch), s1 = it->first().side().m_index;
        int p2 = static_cast<int>(it->second().patch), s2 = it->second().side().m_index;
        const auto& g1 = (s1 <= 2) ? gv[p1] : gu[p1];
        const auto& g2 = (s2 <= 2) ? gv[p2] : gu[p2];
        bool okMatch = (g1.size() == g2.size() && g1.size() >= 2);
        if (okMatch) {
            const int n = static_cast<int>(g1.size());
            gsPoint2D<T> a0 = sideCPPos(p1, s1, 0), a1 = sideCPPos(p1, s1, n - 1);
            gsPoint2D<T> b0 = sideCPPos(p2, s2, 0), b1 = sideCPPos(p2, s2, n - 1);
            bool flipped = (a0.dist(b1) + a1.dist(b0)) < (a0.dist(b0) + a1.dist(b1));
            if ((flipped ? (a0.dist(b1) + a1.dist(b0)) : (a0.dist(b0) + a1.dist(b1))) > T(1e-3))
                okMatch = false;
            for (int k = 0; k < n && okMatch; ++k) {
                T g2Match = flipped ? (T(1) - g2[n - 1 - k]) : g2[k];
                if (std::abs(g1[k] - g2Match) > T(1e-6))
                    okMatch = false;
            }
        }
        if (!okMatch) {
            auto markSide = [&](int p, int s) {
                if (s == 1) { vFixed[pVertId[p][0]] = vFixed[pVertId[p][3]] = true; }
                if (s == 2) { vFixed[pVertId[p][1]] = vFixed[pVertId[p][2]] = true; }
                if (s == 3) { vFixed[pVertId[p][0]] = vFixed[pVertId[p][1]] = true; }
                if (s == 4) { vFixed[pVertId[p][3]] = vFixed[pVertId[p][2]] = true; }
            };
            markSide(p1, s1);
            markSide(p2, s2);
        }
    }

    auto evalIncScore = [&](int vid, int gridN, T& minDetOut) -> T {
        minDetOut = std::numeric_limits<T>::max();
        T minEf = T(180);
        T maxAR = T(1);
        for (const auto& ic : vInc[vid]) {
            T angs[4], mnA = 0, mxA = 180; bool deg = false;
            evalPatchCornerAnglesDeg(*tbs[ic.p], angs, mnA, mxA, deg);
            if (deg) return -T(1e5);
            minEf = std::min(minEf, std::min(mnA, T(180) - mxA));
            gsPoint2D<T> c0 = cornerPos(ic.p, 0, 0), c1 = cornerPos(ic.p, 1, 0);
            gsPoint2D<T> c2 = cornerPos(ic.p, 1, 1), c3 = cornerPos(ic.p, 0, 1);
            T es[4] = {c0.dist(c1), c1.dist(c2), c2.dist(c3), c3.dist(c0)};
            T mnE = std::min(std::min(es[0], es[1]), std::min(es[2], es[3]));
            T mxE = std::max(std::max(es[0], es[1]), std::max(es[2], es[3]));
            if (mnE > T(1e-12)) maxAR = std::max(maxAR, mxE / mnE);
        }
        for (const auto& ic : vInc[vid]) {
            T dN = evalPatchMinDetDense(*tbs[ic.p], gridN);
            minDetOut = std::min(minDetOut, dN);
            if (dN <= T(1e-8)) return -T(1e6) + dN;
        }
        return minEf - T(2.5) * std::max(T(0), maxAR - T(3.0)) + T(4.0) * std::log(std::max(T(1e-9), minDetOut));
    };

    bool anyMoved = false;
    for (int sw = 0; sw < numSweeps; ++sw) {
        bool sweepMoved = false;
        for (size_t vid = 0; vid < verts.size(); ++vid) {
            if (vFixed[vid] || vInc[vid].size() < 3) continue;

            // Compute 1-ring neighbor vertices (deduplicated)
            std::vector<int> nbs;
            for (const auto& ic : vInc[vid]) {
                int cSelf = (ic.ci == 0 && ic.cj == 0) ? 0 : ((ic.ci == 1 && ic.cj == 0) ? 1 : ((ic.ci == 1 && ic.cj == 1) ? 2 : 3));
                int nb1 = pVertId[ic.p][(cSelf + 1) % 4];
                int nb2 = pVertId[ic.p][(cSelf + 3) % 4];
                for (int nb : {nb1, nb2})
                    if (nb != static_cast<int>(vid) && std::find(nbs.begin(), nbs.end(), nb) == nbs.end())
                        nbs.push_back(nb);
            }
            if (nbs.empty()) continue;

            gsPoint2D<T> target(0, 0);
            for (int nb : nbs) { target.x += verts[nb].x; target.y += verts[nb].y; }
            target.x /= static_cast<T>(nbs.size());
            target.y /= static_cast<T>(nbs.size());
            gsPoint2D<T> diff(target.x - verts[vid].x, target.y - verts[vid].y);
            if (std::hypot(diff.x, diff.y) <= T(1e-6)) continue;

            T curDet = 0;
            T curSc = evalIncScore(static_cast<int>(vid), 17, curDet);
            if (curDet <= T(1e-8)) continue;

            // Backup incident patch coefs
            std::vector<gsMatrix<T>> saved;
            saved.reserve(vInc[vid].size());
            for (const auto& ic : vInc[vid]) saved.push_back(tbs[ic.p]->coefs());

            T bestAl = T(0), bestSc = curSc;
            for (T al : {T(0.65), T(0.45), T(0.28), T(0.15), T(0.07)}) {
                T dx = al * diff.x, dy = al * diff.y;
                for (size_t k = 0; k < vInc[vid].size(); ++k) {
                    const auto& ic = vInc[vid][k];
                    gsMatrix<T>& cf = tbs[ic.p]->coefs();
                    cf = saved[k];
                    for (int j = 0; j < nv[ic.p]; ++j) {
                        T wv = ic.cj ? gv[ic.p][j] : (T(1) - gv[ic.p][j]);
                        for (int i = 0; i < nu[ic.p]; ++i) {
                            T wu = ic.ci ? gu[ic.p][i] : (T(1) - gu[ic.p][i]);
                            T w = wu * wv;
                            cf(i + j * nu[ic.p], 0) += w * dx;
                            cf(i + j * nu[ic.p], 1) += w * dy;
                        }
                    }
                }
                T trDet = 0;
                T trSc = evalIncScore(static_cast<int>(vid), 17, trDet);
                if (trDet > T(1e-8) && trSc > bestSc + T(0.05)) {
                    // Verify on denser 25x25 grid before accepting
                    bool ok25 = true;
                    for (const auto& ic : vInc[vid])
                        if (evalPatchMinDetDense(*tbs[ic.p], 25) <= T(1e-8)) { ok25 = false; break; }
                    if (ok25) {
                        bestSc = trSc;
                        bestAl = al;
                        break;
                    }
                }
            }

            if (bestAl > T(0)) {
                verts[vid].x += bestAl * diff.x;
                verts[vid].y += bestAl * diff.y;
                sweepMoved = true;
                anyMoved = true;
            } else {
                for (size_t k = 0; k < vInc[vid].size(); ++k)
                    tbs[vInc[vid][k].p]->coefs() = saved[k];
            }
        }
        if (!sweepMoved) break;
    }
    return anyMoved;
}

template <class T>
inline bool relaxBadCornerAngles(gsMultiPatch<T>& mp, T effMinAng, T effMaxAng, int numSweeps = 3)
{
    const int nP = static_cast<int>(mp.nPatches());
    if (nP < 2) return false;

    mp.computeTopology(T(1e-4), false, false);
    if (2 * mp.nInterfaces() + mp.nBoundary() != 4 * mp.nPatches()) return false;

    std::vector<gsTensorBSpline<2, T>*> tbs(nP, nullptr);
    std::vector<int> nu(nP, 0), nv(nP, 0);
    std::vector<std::vector<T>> gu(nP), gv(nP);
    for (int p = 0; p < nP; ++p) {
        tbs[p] = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tbs[p]) return false;
        nu[p] = tbs[p]->basis().size(0);
        nv[p] = tbs[p]->basis().size(1);
        gu[p] = computeNormalizedGreville(tbs[p]->knots(0));
        gv[p] = computeNormalizedGreville(tbs[p]->knots(1));
    }

    auto sideLen = [&](int p, int s) -> int { return (s <= 2) ? nv[p] : nu[p]; };
    auto sideCP = [&](int p, int s, int k) -> int {
        if (s == 1) return 0 + k * nu[p];
        if (s == 2) return (nu[p] - 1) + k * nu[p];
        if (s == 3) return k + 0 * nu[p];
        return k + (nv[p] - 1) * nu[p];
    };
    auto cpPt = [&](int p, int loc) -> gsPoint2D<T> {
        return gsPoint2D<T>(tbs[p]->coefs()(loc, 0), tbs[p]->coefs()(loc, 1));
    };

    auto applyHandleDelta = [&](int p, int s, int kEdge, T dx, T dy) {
        gsMatrix<T>& cf = tbs[p]->coefs();
        if (s == 1 || s == 2) {
            // West (i=0) or East (i=nu-1), fixed row j = kEdge
            for (int i = 0; i < nu[p]; ++i) {
                T w = (s == 1) ? (T(1) - gu[p][i]) : gu[p][i];
                cf(i + kEdge * nu[p], 0) += w * dx;
                cf(i + kEdge * nu[p], 1) += w * dy;
            }
        } else {
            // South (j=0) or North (j=nv-1), fixed col i = kEdge
            for (int j = 0; j < nv[p]; ++j) {
                T w = (s == 3) ? (T(1) - gv[p][j]) : gv[p][j];
                cf(kEdge + j * nu[p], 0) += w * dx;
                cf(kEdge + j * nu[p], 1) += w * dy;
            }
        }
    };

    auto patchAngleMargin = [&](int p) -> T {
        T angs[4], mnA = 0, mxA = 180; bool deg = false;
        evalPatchCornerAnglesDeg(*tbs[p], angs, mnA, mxA, deg);
        if (deg) return -T(1e4);
        T mMin = mnA - effMinAng;
        T mMax = effMaxAng - mxA;
        T margin = std::min(mMin, mMax);
        T ef = std::min(mnA, T(180) - mxA);
        return margin + T(0.05) * ef;
    };

    const T targetMin = std::max(effMinAng + T(8.0), T(24.0));
    const T targetMax = std::min(effMaxAng - T(8.0), T(150.0));
    const T degToRad = T(3.14159265358979323846) / T(180);
    const T rotDegs[10] = {T(-32), T(-22), T(-14), T(-8), T(-4), T(4), T(8), T(14), T(22), T(32)};

    bool anyFixed = false;
    for (int sw = 0; sw < numSweeps; ++sw) {
        bool sweepChanged = false;
        for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
            int p1 = static_cast<int>(it->first().patch), s1 = it->first().side().m_index;
            int p2 = static_cast<int>(it->second().patch), s2 = it->second().side().m_index;
            if (p1 == p2) continue;
            int n1 = sideLen(p1, s1), n2 = sideLen(p2, s2);
            if (n1 != n2 || n1 < 4) continue;

            // Verify matching endpoint alignment AND orientation-aware Greville equality!
            gsPoint2D<T> a0 = cpPt(p1, sideCP(p1, s1, 0)), a1 = cpPt(p1, sideCP(p1, s1, n1 - 1));
            gsPoint2D<T> b0 = cpPt(p2, sideCP(p2, s2, 0)), b1 = cpPt(p2, sideCP(p2, s2, n2 - 1));
            bool flipped = (a0.dist(b1) + a1.dist(b0)) < (a0.dist(b0) + a1.dist(b1));
            if ((flipped ? (a0.dist(b1) + a1.dist(b0)) : (a0.dist(b0) + a1.dist(b1))) > T(1e-3)) continue;
            const auto& g1 = (s1 <= 2) ? gv[p1] : gu[p1];
            const auto& g2 = (s2 <= 2) ? gv[p2] : gu[p2];
            bool grevOk = true;
            for (int k = 0; k < n1; ++k) {
                T g2m = flipped ? (T(1) - g2[n1 - 1 - k]) : g2[k];
                if (std::abs(g1[k] - g2m) > T(1e-6)) { grevOk = false; break; }
            }
            if (!grevOk) continue;

            // Quick trigger check: skip if both patches already have comfortable angles
            T a1s[4], mn1 = 0, mx1 = 180, a2s[4], mn2 = 0, mx2 = 180;
            bool d1 = false, d2 = false;
            evalPatchCornerAnglesDeg(*tbs[p1], a1s, mn1, mx1, d1);
            evalPatchCornerAnglesDeg(*tbs[p2], a2s, mn2, mx2, d2);
            if (!d1 && !d2 && mn1 >= targetMin && mx1 <= targetMax && mn2 >= targetMin && mx2 <= targetMax)
                continue;

            // Smooth multi-CP band width (1 CP for single-span degree-3; ~25% of edge for multi-knot curves!)
            const int kBand = std::min(n1 - 2, std::max(1, n1 / 4));

            // Try rotating the tangent handle smoothly at end 0 and end 1
            for (int endIdx = 0; endIdx < 2; ++endIdx) {
                int k0 = (endIdx == 0) ? 0 : (n1 - 1);
                int kStep = (endIdx == 0) ? 1 : -1;

                T baseSc = std::min(patchAngleMargin(p1), patchAngleMargin(p2));
                T baseMinDet = std::min(evalPatchMinDetDense(*tbs[p1], 33),
                                        evalPatchMinDetDense(*tbs[p2], 33));

                gsPoint2D<T> P0 = cpPt(p1, sideCP(p1, s1, k0));
                gsMatrix<T> save1 = tbs[p1]->coefs();
                gsMatrix<T> save2 = tbs[p2]->coefs();

                auto applySmoothBand = [&](T rdeg, T scale) {
                    tbs[p1]->coefs() = save1;
                    tbs[p2]->coefs() = save2;
                    for (int b = 1; b <= kBand; ++b) {
                        int k1 = k0 + b * kStep;
                        int k2 = flipped ? (n1 - 1 - k1) : k1;
                        T frac = static_cast<T>(b - 1) / static_cast<T>(std::max(1, kBand));
                        T decay = (T(1) - frac) * (T(1) - frac);
                        T rad = (rdeg * decay) * degToRad;
                        T ca = std::cos(rad), sa = std::sin(rad);
                        T sc = T(1) + decay * (scale - T(1));
                        int loc1 = sideCP(p1, s1, k1);
                        T hx = save1(loc1, 0) - P0.x;
                        T hy = save1(loc1, 1) - P0.y;
                        T dx = sc * (ca * hx - sa * hy) - hx;
                        T dy = sc * (sa * hx + ca * hy) - hy;
                        applyHandleDelta(p1, s1, k1, dx, dy);
                        applyHandleDelta(p2, s2, k2, dx, dy);
                    }
                };

                T bestSc = baseSc;
                T bestDeg = T(0), bestScale = T(1);
                bool found = false;

                for (T rdeg : rotDegs) {
                    for (T sc : {T(1.0), T(0.85)}) {
                        applySmoothBand(rdeg, sc);
                        // Fast angle filter FIRST (microseconds) before dense 33x33 Jacobian check!
                        T trSc = std::min(patchAngleMargin(p1), patchAngleMargin(p2));
                        if (trSc <= bestSc + T(0.25)) continue;

                        T td1 = evalPatchMinDetDense(*tbs[p1], 33);
                        T td2 = evalPatchMinDetDense(*tbs[p2], 33);
                        T trMinDet = std::min(td1, td2);
                        bool detOk = (baseMinDet > T(1e-8)) ? (trMinDet > T(1e-8)) : (trMinDet >= baseMinDet - T(1e-9));
                        if (!detOk) continue;

                        bestSc = trSc;
                        bestDeg = rdeg; bestScale = sc;
                        found = true;
                    }
                }

                if (found) {
                    applySmoothBand(bestDeg, bestScale);
                    sweepChanged = true;
                    anyFixed = true;
                } else {
                    tbs[p1]->coefs() = save1;
                    tbs[p2]->coefs() = save2;
                }
            }
        }
        if (!sweepChanged) break;
    }
    return anyFixed;
}

template <class T, class AcceptCandidateFn>
inline void runCrossFieldBlossomQuadEngine(
    const gsPlanarDomain<T>& dom,
    const TriMeshData<T>& fineMesh,
    const CrossFieldSolution<T>& cf,
    const CrossFieldMeshAccel<T>& accel,
    int degree,
    T minAngTol,
    T maxAngTol,
    int maxPatchesCap,
    AcceptCandidateFn&& acceptCandidateFn)
{
    const gsMatrix<T> bbox = dom.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    if (diag <= T(1e-12)) return;

    const T effMinAng = std::min(minAngTol, std::max(T(4.0), minAngTol - T(3.5)));
    const T effMaxAng = std::max(maxAngTol, T(164.2));
    const T relMinPair = std::max(T(3.0), effMinAng - T(3.5));
    const T relMaxPair = std::min(T(175.0), effMaxAng + T(10.0));

    struct BQConfig {
        bool refineEdges;
        T hDiv;
        T maxBendRad;
        bool refineHolePerim;
        bool addFeatureProj;
        bool injectSings;
        bool needleRefine;
        int bisectSplits;
        int earMode;
    };
    // Progressive coarseness ladder: NEVER run obtuse bisection on the first 4 configs so
    // N_tri stays at the theoretical minimum B + 2*V_int + 2*(H - 1)!
    const BQConfig configs[8] = {
        {false, T(3.0), T(1.20), false, false, false, false, 0,  0}, // 0: Theoretical minimum coarse!
        {false, T(3.0), T(1.05), false, false, false, true,  0,  0}, // 1: Ultra-coarse + needle Steiner rescue
        {true,  T(2.4), T(0.95), false, false, false, false, 0,  0}, // 2: Ultra-coarse curved boundary (minimal swirl/pinch!)
        {true,  T(2.8), T(0.82), false, true,  false, true,  0,  0}, // 3: Coarse + feature proj
        {true,  T(3.2), T(0.76), false, true,  true,  true,  0,  0}, // 4: Coarse + MBO singularity Steiner points
        {true,  T(4.4), T(0.60), false, true,  false, true,  6,  0}, // 5: Moderate curvature (swirls/amoebas)
        {true,  T(5.6), T(0.48), true,  true,  true,  true,  10, 0}, // 6: Medium-fine + MBO singularities
        {true,  T(7.2), T(0.38), true,  false, false, true,  14, 1}  // 7: Fine rescue for extreme star/swirl domains
    };

    for (int ci = 0; ci < 8; ++ci) {
        const auto& cfg = configs[ci];
        TriMeshData<T> mesh;
        if (!buildConstrainedDomainTriangulation(
                dom, degree, diag / cfg.hDiv, cfg.refineEdges, cfg.earMode, mesh,
                cfg.maxBendRad, cfg.refineHolePerim, cfg.addFeatureProj, /*ensureEvenParity=*/true))
            continue;
        if (mesh.tris.empty()) continue;

        applyDomainLawsonFlips(mesh);
        if (cfg.injectSings && !cf.singularities.empty())
            injectCrossFieldSingularities(mesh, cf, diag / T(18), 8);
        if (cfg.needleRefine)
            refineNeedleCircumcenters(mesh, diag / cfg.hDiv, 2, T(13.0), true, (ci <= 3) ? 8 : 16);
        fixUndervalentAndReflexCorners(mesh, degree, (ci <= 2) ? 4 : 8);
        if (cfg.bisectSplits > 0)
            bisectObtuseTriangles(mesh, degree, T(152.0), cfg.bisectSplits, std::max(T(10.0), effMinAng));

        // Restore even triangle parity if a boundary split flipped parity
        if ((mesh.tris.size() % 2) == 1) {
            int bestTi = -1, bestK = -1;
            T bestLen = T(-1);
            for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
                for (int k = 0; k < 3; ++k) {
                    int eid = mesh.tris[ti].e[k];
                    if (!mesh.edges[eid].isBoundary) continue;
                    T L = mesh.pts[mesh.edges[eid].a].dist(mesh.pts[mesh.edges[eid].b]);
                    if (L > bestLen) { bestLen = L; bestTi = ti; bestK = k; }
                }
            }
            if (bestTi >= 0) {
                int ord[3] = {bestK, (bestK + 1) % 3, (bestK + 2) % 3};
                trySplitTriangleEdges(mesh, bestTi, ord, degree, T(3.0));
            }
        }

        // Rotate interior diagonals to align quad edges with the MBO 4-RoSy cross-field
        flipEdgesForCrossField(mesh, fineMesh, cf, accel, effMinAng, effMaxAng);

        const int nT = static_cast<int>(mesh.tris.size());
        if (nT == 0 || nT / 2 > maxPatchesCap) continue;

        for (auto& er : mesh.edges) {
            if (!er.seg.bsp)
                er.seg = makeStraightSeg(mesh.pts[er.a], mesh.pts[er.b], degree);
        }

        // Build Cross-Field-Weighted Pairing Graph
        struct PairEdge { int t0, t1; typename TriMeshData<T>::Quad q; T effMin; T score; bool relaxed; };
        std::vector<PairEdge> pairEdges;
        int triHint = -1;

        for (int eIdx = 0; eIdx < static_cast<int>(mesh.edges.size()); ++eIdx) {
            if (mesh.edges[eIdx].isBoundary) continue;
            int tA = -1, kA = -1, tB = -1, kB = -1;
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].e[k] == eIdx) {
                        if (tA < 0) { tA = ti; kA = k; } else { tB = ti; kB = k; }
                    }
                }
            }
            if (tA < 0 || tB < 0) continue;
            const auto& A = mesh.tris[tA];
            const auto& B = mesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            if (B.v[kB] != b || B.v[(kB + 1) % 3] != a) continue;
            int qv[4] = {c, a, d, b};
            int qe[4] = {A.e[(kA + 2) % 3], B.e[(kB + 1) % 3], B.e[(kB + 2) % 3], A.e[(kA + 1) % 3]};

            T cfSum = T(0);
            for (int k = 0; k < 4; ++k)
                cfSum += evalSegmentCrossFieldAlignment(fineMesh, cf, accel, mesh.pts[qv[k]], mesh.pts[qv[(k + 1) % 4]], triHint);
            T cfDiag = evalSegmentCrossFieldAlignment(fineMesh, cf, accel, mesh.pts[a], mesh.pts[b], triHint);
            T cfBonus = T(8.0) * (T(0.25) * cfSum - T(0.35) * cfDiag);

            T effMin = 0, sc = 0;
            if (evalMeshQuadUnion(mesh, qv, qe, effMinAng, effMaxAng, effMin, sc, false)) {
                pairEdges.push_back({tA, tB, {{c, a, d, b}, {qe[0], qe[1], qe[2], qe[3]}}, effMin, sc + cfBonus, false});
            } else if (evalMeshQuadUnion(mesh, qv, qe, relMinPair, relMaxPair, effMin, sc, true)) {
                pairEdges.push_back({tA, tB, {{c, a, d, b}, {qe[0], qe[1], qe[2], qe[3]}}, effMin, sc + cfBonus - T(28), true});
            }
        }
        std::sort(pairEdges.begin(), pairEdges.end(),
                  [](const PairEdge& p, const PairEdge& q) { return p.score > q.score; });

        std::vector<int> initBaseDeg(mesh.pts.size(), 0);
        for (const auto& t : mesh.tris) { ++initBaseDeg[t.v[0]]; ++initBaseDeg[t.v[1]]; ++initBaseDeg[t.v[2]]; }

        // Upgraded cross-field + orthogonal-arrival curved interior edge optimizer (supports asymmetric weights!)
        auto curveInteriorEdgesForQuads = [&](std::vector<typename TriMeshData<T>::EdgeRec>& edges,
                                              const std::vector<typename TriMeshData<T>::Quad>& qList,
                                              const std::vector<bool>& needSp) {
            std::vector<std::vector<int>> eQuads(edges.size());
            for (int qi = 0; qi < static_cast<int>(qList.size()); ++qi)
                for (int k = 0; k < 4; ++k) eQuads[qList[qi].e[k]].push_back(qi);

            auto getSeg = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                return (edges[eIdx].a == u && edges[eIdx].b == v) ? edges[eIdx].seg : edges[eIdx].seg.reversed();
            };
            auto evalIncMin = [&](int eIdx) -> T {
                T mnAll = T(180);
                for (int qi : eQuads[eIdx]) {
                    const auto& q = qList[qi];
                    if (needSp[q.e[0]] || needSp[q.e[1]]) continue;
                    T mnA = 0, mxA = 180, ef = 0;
                    eval4SegCornerAnglesDeg(
                        getSeg(q.e[0], q.v[0], q.v[1]), getSeg(q.e[1], q.v[1], q.v[2]),
                        getSeg(q.e[2], q.v[2], q.v[3]), getSeg(q.e[3], q.v[3], q.v[0]),
                        mnA, mxA, ef);
                    std::unique_ptr<gsTensorBSpline<2, T>> qp(
                        makeCoonsFromCCWCycle(getSeg(q.e[0], q.v[0], q.v[1]), getSeg(q.e[1], q.v[1], q.v[2]),
                                              getSeg(q.e[2], q.v[2], q.v[3]), getSeg(q.e[3], q.v[3], q.v[0]), degree));
                    if (!qp || !isPatchDetSingleSign(*qp, 13, T(1e-8))) ef -= T(200);
                    mnAll = std::min(mnAll, ef);
                }
                return mnAll;
            };

            int locHint = -1;
            for (int sweep = 0; sweep < 2; ++sweep) {
                bool anyMoved = false;
                for (int ei = 0; ei < static_cast<int>(edges.size()); ++ei) {
                    if (edges[ei].isBoundary || needSp[ei] || eQuads[ei].empty()) continue;
                    T baseSc = evalIncMin(ei);
                    if (baseSc >= T(32.0)) continue;

                    int va = edges[ei].a, vb = edges[ei].b;
                    const auto& A = mesh.pts[va];
                    const auto& B = mesh.pts[vb];
                    T L = A.dist(B);
                    if (L <= T(1e-9)) continue;
                    gsPoint2D<T> dAB((B.x - A.x) / L, (B.y - A.y) / L);
                    gsPoint2D<T> perp(-dAB.y, dAB.x);

                    auto cfDirAt = [&](int vIdx, const gsPoint2D<T>& fallback) -> gsPoint2D<T> {
                        if (mesh.isBoundaryPt[vIdx]) {
                            gsPoint2D<T> tAvg(mesh.bndTanIn[vIdx].x + mesh.bndTanOut[vIdx].x,
                                              mesh.bndTanIn[vIdx].y + mesh.bndTanOut[vIdx].y);
                            T tn = std::hypot(tAvg.x, tAvg.y);
                            if (tn > T(1e-12)) {
                                gsPoint2D<T> n1(-tAvg.y / tn, tAvg.x / tn), n2(tAvg.y / tn, -tAvg.x / tn);
                                gsPoint2D<T> nUse = (n1.x * fallback.x + n1.y * fallback.y >= n2.x * fallback.x + n2.y * fallback.y) ? n1 : n2;
                                gsPoint2D<T> probe(mesh.pts[vIdx].x + T(0.05) * L * nUse.x,
                                                   mesh.pts[vIdx].y + T(0.05) * L * nUse.y);
                                if (isPointStrictlyInsideDomain(probe, mesh.loopPolys, T(1e-6))) return nUse;
                            }
                        }
                        std::complex<T> zq; T th = T(0);
                        if (evalCrossFieldAt(fineMesh, cf, accel, mesh.pts[vIdx], zq, th, locHint)) {
                            const T pi2 = T(1.57079632679489662);
                            T bestD = -T(2); gsPoint2D<T> bestV = fallback;
                            for (int k = 0; k < 4; ++k) {
                                gsPoint2D<T> cd(std::cos(th + k * pi2), std::sin(th + k * pi2));
                                T d = cd.x * fallback.x + cd.y * fallback.y;
                                if (d > bestD) { bestD = d; bestV = cd; }
                            }
                            if (bestD > T(0.25)) return bestV;
                        }
                        return fallback;
                    };

                    gsPoint2D<T> nA = cfDirAt(va, dAB);
                    gsPoint2D<T> nB = cfDirAt(vb, gsPoint2D<T>(-dAB.x, -dAB.y));
                    CurveSeg<T> oldSeg = edges[ei].seg, bestSeg = oldSeg;
                    T bestSc = baseSc;

                    std::vector<CurveSeg<T>> segCands;
                    const T wPairs[7][2] = {
                        {T(0.45), T(0.45)}, {T(0.75), T(0.75)}, {T(0.95), T(0.95)},
                        {T(0.75), T(0.0)},  {T(0.0),  T(0.75)}, {T(0.95), T(0.25)}, {T(0.25), T(0.95)}
                    };
                    for (int wi = 0; wi < 7; ++wi) {
                        gsPoint2D<T> tA((T(1) - wPairs[wi][0]) * dAB.x + wPairs[wi][0] * nA.x,
                                        (T(1) - wPairs[wi][0]) * dAB.y + wPairs[wi][0] * nA.y);
                        gsPoint2D<T> tB((T(1) - wPairs[wi][1]) * (-dAB.x) + wPairs[wi][1] * nB.x,
                                        (T(1) - wPairs[wi][1]) * (-dAB.y) + wPairs[wi][1] * nB.y);
                        T na = std::hypot(tA.x, tA.y), nb = std::hypot(tB.x, tB.y);
                        if (na > T(1e-12) && nb > T(1e-12)) {
                            tA.x /= na; tA.y /= na; tB.x /= nb; tB.y /= nb;
                            for (T hFrac : {T(0.33), T(0.22)}) {
                                segCands.push_back(makeCubicBezierInteriorSeg(
                                    A, gsPoint2D<T>(A.x + hFrac * L * tA.x, A.y + hFrac * L * tA.y),
                                    gsPoint2D<T>(B.x + hFrac * L * tB.x, B.y + hFrac * L * tB.y), B, degree));
                            }
                        }
                    }
                    for (T bow : {T(-0.14), T(0.14), T(-0.24), T(0.24), T(-0.34), T(0.34)}) {
                        segCands.push_back(makeCubicBezierInteriorSeg(
                            A, gsPoint2D<T>(A.x + (L / T(3)) * dAB.x + bow * L * perp.x, A.y + (L / T(3)) * dAB.y + bow * L * perp.y),
                            gsPoint2D<T>(A.x + (T(2) * L / T(3)) * dAB.x + bow * L * perp.x, A.y + (T(2) * L / T(3)) * dAB.y + bow * L * perp.y),
                            B, degree));
                    }

                    for (const auto& cs : segCands) {
                        bool okIn = true;
                        for (int s = 1; s < 10; ++s) {
                            if (!isPointStrictlyInsideDomain(evalCurvePoint(*cs.bsp, static_cast<T>(s) / T(10)), mesh.loopPolys, T(1e-5))) {
                                okIn = false; break;
                            }
                        }
                        if (!okIn) continue;
                        edges[ei].seg = cs;
                        T sc = evalIncMin(ei);
                        if (sc > bestSc + T(1.0)) { bestSc = sc; bestSeg = cs; anyMoved = true; }
                    }
                    edges[ei].seg = bestSeg;
                }
                if (!anyMoved) break;
            }
        };

        // Helper: run Blossom matching with optional banned pair set and relaxed flag
        auto runMatching = [&](bool allowRelaxed, const std::vector<bool>& bannedPair,
                               std::vector<typename TriMeshData<T>::Quad>& joinedQuads,
                               std::vector<int>& quadPairIdx,
                               std::vector<typename TriMeshData<T>::Tri>& remTris) {
            std::vector<int> baseDeg = initBaseDeg;
            std::vector<int> matchedPair(nT, -1), partnerTri(nT, -1);
            std::vector<std::vector<std::pair<int, int>>> triAdjStrict(nT), triAdjAll(nT);
            for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                if (bannedPair[pi]) continue;
                const auto& pe = pairEdges[pi];
                if (pe.relaxed && !allowRelaxed) continue;
                triAdjAll[pe.t0].push_back({pe.t1, pi});
                triAdjAll[pe.t1].push_back({pe.t0, pi});
                if (!pe.relaxed) {
                    triAdjStrict[pe.t0].push_back({pe.t1, pi});
                    triAdjStrict[pe.t1].push_back({pe.t0, pi});
                }
            }

            for (int tier = 0; tier < 2; ++tier) {
                if (tier == 1 && !allowRelaxed) break;
                for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                    if (bannedPair[pi]) continue;
                    const auto& pe = pairEdges[pi];
                    if (tier == 0 && pe.relaxed) continue;
                    if (matchedPair[pe.t0] >= 0 || matchedPair[pe.t1] >= 0) continue;
                    int a = pe.q.v[1], b = pe.q.v[3];
                    if (baseDeg[a] <= minRequiredBoundaryTriValence(mesh, a) ||
                        baseDeg[b] <= minRequiredBoundaryTriValence(mesh, b))
                        continue;
                    matchedPair[pe.t0] = pi; matchedPair[pe.t1] = pi;
                    partnerTri[pe.t0] = pe.t1; partnerTri[pe.t1] = pe.t0;
                    --baseDeg[a]; --baseDeg[b];
                }
            }

            auto runBlossom = [&](const std::vector<std::vector<std::pair<int, int>>>& useAdj) {
                for (int ti = 0; ti < nT; ++ti) {
                    if (matchedPair[ti] >= 0) continue;
                    augmentBlossomMatchingFrom(
                        ti, nT, useAdj, partnerTri, matchedPair,
                        [&](const std::vector<int>& addPairs, const std::vector<int>& remPairs) -> bool {
                            std::map<int, int> delta;
                            for (int pi : remPairs) { ++delta[pairEdges[pi].q.v[1]]; ++delta[pairEdges[pi].q.v[3]]; }
                            for (int pi : addPairs) { --delta[pairEdges[pi].q.v[1]]; --delta[pairEdges[pi].q.v[3]]; }
                            for (const auto& kv : delta)
                                if (baseDeg[kv.first] + kv.second < minRequiredBoundaryTriValence(mesh, kv.first))
                                    return false;
                            for (const auto& kv : delta) baseDeg[kv.first] += kv.second;
                            return true;
                        });
                }
            };
            runBlossom(triAdjStrict);
            if (allowRelaxed) runBlossom(triAdjAll);

            joinedQuads.clear();
            quadPairIdx.clear();
            remTris.clear();
            for (int ti = 0; ti < nT; ++ti) {
                if (matchedPair[ti] >= 0) {
                    if (partnerTri[ti] > ti) {
                        joinedQuads.push_back(pairEdges[matchedPair[ti]].q);
                        quadPairIdx.push_back(matchedPair[ti]);
                    }
                } else {
                    remTris.push_back(mesh.tris[ti]);
                }
            }
        };

        // Helper: try Dual-Path Zipper (Global Shortest-Valid-Path-First) + Localized Polychord Split
        auto emitMixedWithRescues = [&](const std::vector<typename TriMeshData<T>::Quad>& joinedQuads,
                                        const std::vector<typename TriMeshData<T>::Tri>& remTris) -> bool {
            // 1. Try Global Shortest-Valid-Path Dual-Path Unmatched-Triangle Zipper (minimal patch increase!)
            if (!remTris.empty() && (remTris.size() % 2) == 0 && !joinedQuads.empty()) {
                const int nR = static_cast<int>(remTris.size());
                const int nQ = static_cast<int>(joinedQuads.size());
                const int nF = nQ + nR;
                std::vector<typename TriMeshData<T>::EdgeRec> zEdges = mesh.edges;
                std::vector<bool> noSp(zEdges.size(), false);
                curveInteriorEdgesForQuads(zEdges, joinedQuads, noSp);

                struct FaceSide { int fId; int side; };
                std::vector<std::vector<FaceSide>> edgeInc(zEdges.size());
                for (int qi = 0; qi < nQ; ++qi)
                    for (int s = 0; s < 4; ++s) edgeInc[joinedQuads[qi].e[s]].push_back({qi, s});
                for (int ri = 0; ri < nR; ++ri)
                    for (int s = 0; s < 3; ++s) edgeInc[remTris[ri].e[s]].push_back({nQ + ri, s});

                struct DualAdj { int toF; int eId; int mySide; int toSide; };
                std::vector<std::vector<DualAdj>> dualAdj(nF);
                for (int eid = 0; eid < static_cast<int>(zEdges.size()); ++eid) {
                    if (zEdges[eid].isBoundary || edgeInc[eid].size() != 2) continue;
                    auto a = edgeInc[eid][0], b = edgeInc[eid][1];
                    dualAdj[a.fId].push_back({b.fId, eid, a.side, b.side});
                    dualAdj[b.fId].push_back({a.fId, eid, b.side, a.side});
                }

                auto getZSeg = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                    return (zEdges[eIdx].a == u && zEdges[eIdx].b == v) ? zEdges[eIdx].seg : zEdges[eIdx].seg.reversed();
                };
                auto quadCtr = [&](int qIdx) -> gsPoint2D<T> {
                    const auto& q = joinedQuads[qIdx];
                    return gsPoint2D<T>(
                        T(0.25) * (mesh.pts[q.v[0]].x + mesh.pts[q.v[1]].x + mesh.pts[q.v[2]].x + mesh.pts[q.v[3]].x),
                        T(0.25) * (mesh.pts[q.v[0]].y + mesh.pts[q.v[1]].y + mesh.pts[q.v[2]].y + mesh.pts[q.v[3]].y));
                };

                // Two passes: strict angle tolerance first, then relaxed (for relaxBadCornerAngles rescue!)
                for (int zMode = 0; zMode < 2; ++zMode) {
                    const T minZipA = (zMode == 0) ? std::max(T(3.5), effMinAng - T(3.0)) : T(2.0);
                    const T maxZipA = (zMode == 0) ? std::min(T(172.0), effMaxAng + T(6.0)) : T(176.5);

                    auto pickEndM = [&](int triF, int triSide, int qF, gsPoint2D<T>& outM) -> bool {
                        const auto& t = remTris[triF - nQ];
                        int u = t.v[triSide], v = t.v[(triSide + 1) % 3], opp = t.v[(triSide + 2) % 3];
                        gsPoint2D<T> ctr = quadCtr(qF);
                        const auto& q = joinedQuads[qF];
                        gsPoint2D<T> qp[4] = {mesh.pts[q.v[0]], mesh.pts[q.v[1]], mesh.pts[q.v[2]], mesh.pts[q.v[3]]};
                        T bestSc = -T(1e9); bool ok = false;
                        const T uFracs[3] = {T(0.50), T(0.36), T(0.64)};
                        const T alFracs[6] = {T(0.22), T(0.32), T(0.14), T(0.42), T(0.52), T(0.08)};
                        for (T uf : uFracs) {
                            gsPoint2D<T> base((T(1) - uf) * mesh.pts[u].x + uf * mesh.pts[v].x,
                                              (T(1) - uf) * mesh.pts[u].y + uf * mesh.pts[v].y);
                            for (T al : alFracs) {
                                gsPoint2D<T> cand(base.x + al * (ctr.x - base.x), base.y + al * (ctr.y - base.y));
                                if (!isPointStrictlyInsideDomain(cand, mesh.loopPolys, T(1e-5))) continue;
                                bool inQ = true;
                                for (int k = 0; k < 4; ++k) if (cross2(qp[k], qp[(k + 1) % 4], cand) <= T(1e-7)) { inQ = false; break; }
                                if (!inQ) continue;
                                gsPoint2D<T> tq[4] = {mesh.pts[opp], mesh.pts[u], cand, mesh.pts[v]};
                                if (!isStrictlyConvexQuad(tq)) continue;
                                T mnA, mxA, efA;
                                eval4SegCornerAnglesDeg(getZSeg(t.e[(triSide + 2) % 3], opp, u),
                                                        makeStraightSeg(mesh.pts[u], cand, degree),
                                                        makeStraightSeg(cand, mesh.pts[v], degree),
                                                        getZSeg(t.e[(triSide + 1) % 3], v, opp), mnA, mxA, efA);
                                if (mnA < minZipA || mxA > maxZipA) continue;
                                if (!ok || efA > bestSc) { bestSc = efA; outM = cand; ok = true; }
                            }
                        }
                        return ok;
                    };

                    auto buildPathQuads = [&](const std::vector<int>& pathF, const std::vector<int>& pathE,
                                              const std::vector<int>& sideOut, const std::vector<int>& sideIn,
                                              std::vector<std::array<CurveSeg<T>, 4>>& outPQuads) -> bool {
                        outPQuads.clear();
                        const int L = static_cast<int>(pathE.size());
                        if (L == 1) {
                            const auto& A = remTris[pathF[0] - nQ];
                            const auto& B = remTris[pathF[1] - nQ];
                            int kA = sideOut[0], kB = sideIn[0];
                            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3], d = B.v[(kB + 2) % 3];
                            outPQuads.push_back({
                                getZSeg(A.e[(kA + 2) % 3], c, a), getZSeg(B.e[(kB + 1) % 3], a, d),
                                getZSeg(B.e[(kB + 2) % 3], d, b), getZSeg(A.e[(kA + 1) % 3], b, c)
                            });
                            return true;
                        }
                        std::vector<gsPoint2D<T>> M(L);
                        if (!pickEndM(pathF[0], sideOut[0], pathF[1], M[0]) ||
                            !pickEndM(pathF[L], sideIn[L - 1], pathF[L - 1], M[L - 1]))
                            return false;
                        for (int i = 1; i < L - 1; ++i)
                            M[i] = evalCurvePoint(*zEdges[pathE[i]].seg.bsp, T(0.5));

                        auto pushEndTriQuad = [&](int triF, int s, const gsPoint2D<T>& mPt) {
                            const auto& t = remTris[triF - nQ];
                            int u = t.v[s], v = t.v[(s + 1) % 3], opp = t.v[(s + 2) % 3];
                            outPQuads.push_back({
                                getZSeg(t.e[(s + 2) % 3], opp, u),
                                makeStraightSeg(mesh.pts[u], mPt, degree),
                                makeStraightSeg(mPt, mesh.pts[v], degree),
                                getZSeg(t.e[(s + 1) % 3], v, opp)
                            });
                        };
                        pushEndTriQuad(pathF[0], sideOut[0], M[0]);
                        pushEndTriQuad(pathF[L], sideIn[L - 1], M[L - 1]);

                        for (int i = 1; i < L; ++i) {
                            const auto& q = joinedQuads[pathF[i]];
                            int sA = sideIn[i - 1], sB = sideOut[i];
                            const auto& mA = M[i - 1];
                            const auto& mB = M[i];
                            if (sA == sB) return false;
                            int diff = (sB - sA + 4) % 4;
                            if (diff == 2) {
                                int r = sA;
                                int v0 = q.v[r], v1 = q.v[(r + 1) % 4], v2 = q.v[(r + 2) % 4], v3 = q.v[(r + 3) % 4];
                                int e1 = q.e[(r + 1) % 4], e3 = q.e[(r + 3) % 4];
                                gsPoint2D<T> ca[4] = {mesh.pts[v0], mA, mB, mesh.pts[v3]};
                                gsPoint2D<T> cb[4] = {mA, mesh.pts[v1], mesh.pts[v2], mB};
                                if (!isStrictlyConvexQuad(ca) || !isStrictlyConvexQuad(cb) ||
                                    !isSegmentInsideDomain(mA, mB, mesh.loopPolys)) return false;
                                CurveSeg<T> sMid = makeStraightSeg(mA, mB, degree);
                                outPQuads.push_back({makeStraightSeg(mesh.pts[v0], mA, degree), sMid,
                                                     makeStraightSeg(mB, mesh.pts[v3], degree), getZSeg(e3, v3, v0)});
                                outPQuads.push_back({makeStraightSeg(mA, mesh.pts[v1], degree), getZSeg(e1, v1, v2),
                                                     makeStraightSeg(mesh.pts[v2], mB, degree), sMid.reversed()});
                            } else {
                                int r = (diff == 1) ? sA : sB;
                                const auto& m0 = (diff == 1) ? mA : mB;
                                const auto& m1 = (diff == 1) ? mB : mA;
                                int v0 = q.v[r], v1 = q.v[(r + 1) % 4], v2 = q.v[(r + 2) % 4], v3 = q.v[(r + 3) % 4];
                                int e2 = q.e[(r + 2) % 4], e3 = q.e[(r + 3) % 4];
                                gsPoint2D<T> ctr = quadCtr(pathF[i]);
                                gsPoint2D<T> oCands[4] = {
                                    gsPoint2D<T>((m0.x + m1.x + mesh.pts[v3].x) / T(3), (m0.y + m1.y + mesh.pts[v3].y) / T(3)),
                                    ctr,
                                    gsPoint2D<T>(T(0.5) * (ctr.x + mesh.pts[v3].x), T(0.5) * (ctr.y + mesh.pts[v3].y)),
                                    gsPoint2D<T>(T(0.25) * (m0.x + m1.x + mesh.pts[v0].x + mesh.pts[v2].x),
                                                 T(0.25) * (m0.y + m1.y + mesh.pts[v0].y + mesh.pts[v2].y))
                                };
                                gsPoint2D<T> bestO = ctr; T bestEf = -T(1); bool foundO = false;
                                for (int oi = 0; oi < 4; ++oi) {
                                    const auto& O = oCands[oi];
                                    if (!isPointStrictlyInsideDomain(O, mesh.loopPolys, T(1e-5))) continue;
                                    gsPoint2D<T> c0[4] = {mesh.pts[v0], m0, O, mesh.pts[v3]};
                                    gsPoint2D<T> c1[4] = {m0, mesh.pts[v1], m1, O};
                                    gsPoint2D<T> c2[4] = {O, m1, mesh.pts[v2], mesh.pts[v3]};
                                    if (!isStrictlyConvexQuad(c0) || !isStrictlyConvexQuad(c1) || !isStrictlyConvexQuad(c2)) continue;
                                    CurveSeg<T> sm0 = makeStraightSeg(m0, O, degree), sm1 = makeStraightSeg(O, m1, degree), sv3 = makeStraightSeg(O, mesh.pts[v3], degree);
                                    T mn0, mx0, ef0, mn1, mx1, ef1, mn2, mx2, ef2;
                                    eval4SegCornerAnglesDeg(makeStraightSeg(mesh.pts[v0], m0, degree), sm0, sv3, getZSeg(e3, v3, v0), mn0, mx0, ef0);
                                    eval4SegCornerAnglesDeg(makeStraightSeg(m0, mesh.pts[v1], degree), makeStraightSeg(mesh.pts[v1], m1, degree), sm1.reversed(), sm0.reversed(), mn1, mx1, ef1);
                                    eval4SegCornerAnglesDeg(sm1, makeStraightSeg(m1, mesh.pts[v2], degree), getZSeg(e2, v2, v3), sv3.reversed(), mn2, mx2, ef2);
                                    T efM = std::min(ef0, std::min(ef1, ef2));
                                    if (!foundO || efM > bestEf) { bestEf = efM; bestO = O; foundO = true; }
                                }
                                if (!foundO) return false;
                                CurveSeg<T> sm0 = makeStraightSeg(m0, bestO, degree), sm1 = makeStraightSeg(bestO, m1, degree), sv3 = makeStraightSeg(bestO, mesh.pts[v3], degree);
                                outPQuads.push_back({makeStraightSeg(mesh.pts[v0], m0, degree), sm0, sv3, getZSeg(e3, v3, v0)});
                                outPQuads.push_back({makeStraightSeg(m0, mesh.pts[v1], degree), makeStraightSeg(mesh.pts[v1], m1, degree), sm1.reversed(), sm0.reversed()});
                                outPQuads.push_back({sm1, makeStraightSeg(m1, mesh.pts[v2], degree), getZSeg(e2, v2, v3), sv3.reversed()});
                            }
                        }
                        return true;
                    };

                    std::vector<bool> faceUsed(nF, false);
                    std::vector<std::array<CurveSeg<T>, 4>> zipQuads;
                    int zippedTris = 0;

                    while (zippedTris < nR) {
                        // Find the globally shortest VALID dual path between ANY two unused unmatched triangles!
                        struct CandPath {
                            int len;
                            std::vector<int> pathF, pathE, sideOut, sideIn;
                        };
                        std::vector<CandPath> candPaths;
                        for (int ri = 0; ri < nR; ++ri) {
                            int sF = nQ + ri;
                            if (faceUsed[sF]) continue;
                            std::vector<int> prevF(nF, -1), prevEdge(nF, -1), prevSideFrom(nF, -1), prevSideTo(nF, -1), dist(nF, -1);
                            std::queue<int> bfs;
                            prevF[sF] = sF; dist[sF] = 0;
                            bfs.push(sF);
                            int foundTargets = 0;
                            while (!bfs.empty() && foundTargets < 3) {
                                int u = bfs.front(); bfs.pop();
                                if (dist[u] >= 12) continue;
                                for (const auto& nb : dualAdj[u]) {
                                    int v = nb.toF;
                                    if (faceUsed[v] || prevF[v] >= 0) continue;
                                    if (u == sF && v >= nQ) {
                                        const auto& A = remTris[u - nQ];
                                        const auto& B = remTris[v - nQ];
                                        int kA = nb.mySide, kB = nb.toSide;
                                        gsPoint2D<T> qc[4] = {
                                            mesh.pts[A.v[(kA + 2) % 3]], mesh.pts[A.v[kA]],
                                            mesh.pts[B.v[(kB + 2) % 3]], mesh.pts[A.v[(kA + 1) % 3]]
                                        };
                                        if (!isStrictlyConvexQuad(qc)) continue;
                                    }
                                    prevF[v] = u; prevEdge[v] = nb.eId; prevSideFrom[v] = nb.mySide; prevSideTo[v] = nb.toSide;
                                    dist[v] = dist[u] + 1;
                                    if (v >= nQ) {
                                        if (v > sF) {
                                            CandPath cp;
                                            cp.len = dist[v];
                                            for (int cur = v; cur != sF; cur = prevF[cur]) {
                                                cp.pathF.push_back(cur); cp.pathE.push_back(prevEdge[cur]);
                                                cp.sideOut.push_back(prevSideFrom[cur]); cp.sideIn.push_back(prevSideTo[cur]);
                                            }
                                            cp.pathF.push_back(sF);
                                            std::reverse(cp.pathF.begin(), cp.pathF.end());
                                            std::reverse(cp.pathE.begin(), cp.pathE.end());
                                            std::reverse(cp.sideOut.begin(), cp.sideOut.end());
                                            std::reverse(cp.sideIn.begin(), cp.sideIn.end());
                                            candPaths.push_back(std::move(cp));
                                            ++foundTargets;
                                        }
                                        continue; // don't route through another unmatched triangle
                                    }
                                    bfs.push(v);
                                }
                            }
                        }

                        std::sort(candPaths.begin(), candPaths.end(),
                                  [](const CandPath& a, const CandPath& b) { return a.len < b.len; });

                        bool committed = false;
                        std::vector<std::array<CurveSeg<T>, 4>> pQuads;
                        for (const auto& cp : candPaths) {
                            if (buildPathQuads(cp.pathF, cp.pathE, cp.sideOut, cp.sideIn, pQuads)) {
                                for (const auto& q : pQuads) zipQuads.push_back(q);
                                for (int f : cp.pathF) faceUsed[f] = true;
                                zippedTris += 2;
                                committed = true;
                                break;
                            }
                        }
                        if (!committed) break;
                    }

                    if (zippedTris == nR) {
                        int estP = static_cast<int>(zipQuads.size());
                        for (int qi = 0; qi < nQ; ++qi) if (!faceUsed[qi]) ++estP;
                        if (estP > 0 && estP <= maxPatchesCap + 40) {
                            gsMultiPatch<T> mp;
                            for (int qi = 0; qi < nQ; ++qi) {
                                if (faceUsed[qi]) continue;
                                const auto& q = joinedQuads[qi];
                                mp.addPatch(typename gsGeometry<T>::uPtr(
                                    makeCoonsFromCCWCycle(getZSeg(q.e[0], q.v[0], q.v[1]), getZSeg(q.e[1], q.v[1], q.v[2]),
                                                          getZSeg(q.e[2], q.v[2], q.v[3]), getZSeg(q.e[3], q.v[3], q.v[0]), degree)));
                            }
                            for (const auto& zq : zipQuads)
                                mp.addPatch(typename gsGeometry<T>::uPtr(
                                    makeCoonsFromCCWCycle(zq[0], zq[1], zq[2], zq[3], degree)));
                            if (acceptCandidateFn(mp, "HarmonicCrossField"))
                                return true;
                        }
                    }
                }
            }

            // 2. Selective Polychord-Propagated Split with Adaptive Split-Param & Hermite Spokes
            if (remTris.empty()) return false;
            std::vector<typename TriMeshData<T>::EdgeRec> edges = mesh.edges;
            std::vector<bool> needSplit(edges.size(), false);
            for (const auto& t : remTris)
                needSplit[t.e[0]] = needSplit[t.e[1]] = needSplit[t.e[2]] = true;
            bool changed = true;
            while (changed) {
                changed = false;
                for (const auto& q : joinedQuads) {
                    if (needSplit[q.e[0]] != needSplit[q.e[2]]) { needSplit[q.e[0]] = needSplit[q.e[2]] = true; changed = true; }
                    if (needSplit[q.e[1]] != needSplit[q.e[3]]) { needSplit[q.e[1]] = needSplit[q.e[3]] = true; changed = true; }
                }
            }

            int estP = 3 * static_cast<int>(remTris.size());
            for (const auto& q : joinedQuads) {
                bool s0 = needSplit[q.e[0]], s1 = needSplit[q.e[1]];
                estP += (s0 && s1) ? 4 : ((s0 || s1) ? 2 : 1);
            }
            // Allow headroom for downstream reducePolychordStrips to collapse strips back under maxPatchesCap!
            if (estP > maxPatchesCap + 80) return false;

            curveInteriorEdgesForQuads(edges, joinedQuads, needSplit);
            for (size_t ei = 0; ei < edges.size(); ++ei) {
                if (needSplit[ei]) {
                    auto sp = edges[ei].seg.splitArc(T(0.5));
                    edges[ei].h0 = sp.first; edges[ei].h1 = sp.second;
                    edges[ei].mid = sp.first.p1(); edges[ei].splitDone = true;
                }
            }

            auto getSeg = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                return (edges[eIdx].a == u && edges[eIdx].b == v) ? edges[eIdx].seg : edges[eIdx].seg.reversed();
            };
            auto getHalves = [&](int eIdx, int u, int v, CurveSeg<T>& o0, CurveSeg<T>& o1, gsPoint2D<T>& om) {
                om = edges[eIdx].mid;
                if (edges[eIdx].a == u && edges[eIdx].b == v) { o0 = edges[eIdx].h0; o1 = edges[eIdx].h1; }
                else                                          { o0 = edges[eIdx].h1.reversed(); o1 = edges[eIdx].h0.reversed(); }
            };
            auto chooseQuad4Ctr = [&](const typename TriMeshData<T>::Quad& q,
                                      const CurveSeg<T>& e0a, const CurveSeg<T>& e0b, const gsPoint2D<T>& m0,
                                      const CurveSeg<T>& e1a, const CurveSeg<T>& e1b, const gsPoint2D<T>& m1,
                                      const CurveSeg<T>& e2a, const CurveSeg<T>& e2b, const gsPoint2D<T>& m2,
                                      const CurveSeg<T>& e3a, const CurveSeg<T>& e3b, const gsPoint2D<T>& m3,
                                      gsPoint2D<T>& outO) -> bool {
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                gsPoint2D<T> Om((m0.x + m1.x + m2.x + m3.x) * T(0.25), (m0.y + m1.y + m2.y + m3.y) * T(0.25));
                gsPoint2D<T> Ov((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x + mesh.pts[v3].x) * T(0.25),
                                (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y + mesh.pts[v3].y) * T(0.25));
                std::vector<gsPoint2D<T>> oCands = {Om, Ov,
                    gsPoint2D<T>(T(0.5) * (Om.x + Ov.x), T(0.5) * (Om.y + Ov.y))};
                for (int v : {v0, v1, v2, v3})
                    oCands.emplace_back(T(0.75) * Om.x + T(0.25) * mesh.pts[v].x,
                                        T(0.75) * Om.y + T(0.25) * mesh.pts[v].y);
                bool found = false; T bestEf = -T(1e9);
                for (const auto& Oc : oCands) {
                    gsPoint2D<T> c0[4] = {mesh.pts[v0], m0, Oc, m3};
                    gsPoint2D<T> c1[4] = {mesh.pts[v1], m1, Oc, m0};
                    gsPoint2D<T> c2[4] = {mesh.pts[v2], m2, Oc, m1};
                    gsPoint2D<T> c3[4] = {mesh.pts[v3], m3, Oc, m2};
                    if (!isStrictlyConvexQuad(c0) || !isStrictlyConvexQuad(c1) ||
                        !isStrictlyConvexQuad(c2) || !isStrictlyConvexQuad(c3)) continue;
                    auto s0 = makeStraightSeg(Oc, m0, degree), s1 = makeStraightSeg(Oc, m1, degree);
                    auto s2 = makeStraightSeg(Oc, m2, degree), s3 = makeStraightSeg(Oc, m3, degree);
                    T mn, mx, f0, f1, f2, f3;
                    eval4SegCornerAnglesDeg(e0a, s0.reversed(), s3, e3b, mn, mx, f0);
                    eval4SegCornerAnglesDeg(e1a, s1.reversed(), s0, e0b, mn, mx, f1);
                    eval4SegCornerAnglesDeg(e2a, s2.reversed(), s1, e1b, mn, mx, f2);
                    eval4SegCornerAnglesDeg(e3a, s3.reversed(), s2, e2b, mn, mx, f3);
                    T ef = std::min(std::min(f0, f1), std::min(f2, f3));
                    if (!found || ef > bestEf) { bestEf = ef; outO = Oc; found = true; }
                }
                return found;
            };

            auto evalFaceMinEf = [&](bool isTri, int fIdx) -> T {
                if (isTri) {
                    const auto& t = remTris[fIdx];
                    int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                    CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b; gsPoint2D<T> m0, m1, m2, O;
                    getHalves(t.e[0], v0, v1, e0a, e0b, m0);
                    getHalves(t.e[1], v1, v2, e1a, e1b, m1);
                    getHalves(t.e[2], v2, v0, e2a, e2b, m2);
                    if (!pickOptimalTri1to3CenterCurved(
                            mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], e0a, e0b, e1a, e1b, e2a, e2b,
                            edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O))
                        return -T(100);
                    auto s0 = makeStraightSeg(O, m0, degree), s1 = makeStraightSeg(O, m1, degree), s2 = makeStraightSeg(O, m2, degree);
                    T mn, mx, f0, f1, f2;
                    eval4SegCornerAnglesDeg(e0a, s0.reversed(), s2, e2b, mn, mx, f0);
                    eval4SegCornerAnglesDeg(e1a, s1.reversed(), s0, e0b, mn, mx, f1);
                    eval4SegCornerAnglesDeg(e2a, s2.reversed(), s1, e1b, mn, mx, f2);
                    return std::min(f0, std::min(f1, f2));
                }
                const auto& q = joinedQuads[fIdx];
                bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                if (!sp0 && !sp1) {
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(getSeg(q.e[0], v0, v1), getSeg(q.e[1], v1, v2),
                                            getSeg(q.e[2], v2, v3), getSeg(q.e[3], v3, v0), mn, mx, ef);
                    return ef;
                } else if (sp0 && !sp1) {
                    CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m0, m2;
                    getHalves(q.e[0], v0, v1, e0a, e0b, m0);
                    getHalves(q.e[2], v2, v3, e2a, e2b, m2);
                    auto sMid = makeStraightSeg(m0, m2, degree);
                    T mn, mx, f0, f1;
                    eval4SegCornerAnglesDeg(e0a, sMid, e2b, getSeg(q.e[3], v3, v0), mn, mx, f0);
                    eval4SegCornerAnglesDeg(e0b, getSeg(q.e[1], v1, v2), e2a, sMid.reversed(), mn, mx, f1);
                    return std::min(f0, f1);
                } else if (!sp0 && sp1) {
                    CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m1, m3;
                    getHalves(q.e[1], v1, v2, e1a, e1b, m1);
                    getHalves(q.e[3], v3, v0, e3a, e3b, m3);
                    auto sMid = makeStraightSeg(m1, m3, degree);
                    T mn, mx, f0, f1;
                    eval4SegCornerAnglesDeg(getSeg(q.e[0], v0, v1), e1a, sMid, e3b, mn, mx, f0);
                    eval4SegCornerAnglesDeg(sMid.reversed(), e1b, getSeg(q.e[2], v2, v3), e3a, mn, mx, f1);
                    return std::min(f0, f1);
                } else {
                    CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m0, m1, m2, m3, O;
                    getHalves(q.e[0], v0, v1, e0a, e0b, m0);
                    getHalves(q.e[1], v1, v2, e1a, e1b, m1);
                    getHalves(q.e[2], v2, v3, e2a, e2b, m2);
                    getHalves(q.e[3], v3, v0, e3a, e3b, m3);
                    if (!chooseQuad4Ctr(q, e0a, e0b, m0, e1a, e1b, m1, e2a, e2b, m2, e3a, e3b, m3, O))
                        return -T(100);
                    auto s0 = makeStraightSeg(O, m0, degree), s1 = makeStraightSeg(O, m1, degree);
                    auto s2 = makeStraightSeg(O, m2, degree), s3 = makeStraightSeg(O, m3, degree);
                    T mn, mx, f0, f1, f2, f3;
                    eval4SegCornerAnglesDeg(e0a, s0.reversed(), s3, e3b, mn, mx, f0);
                    eval4SegCornerAnglesDeg(e1a, s1.reversed(), s0, e0b, mn, mx, f1);
                    eval4SegCornerAnglesDeg(e2a, s2.reversed(), s1, e1b, mn, mx, f2);
                    eval4SegCornerAnglesDeg(e3a, s3.reversed(), s2, e2b, mn, mx, f3);
                    return std::min(std::min(f0, f1), std::min(f2, f3));
                }
            };

            // Optimize split parameter & interior midpoint position on tight split edges
            struct IncF { bool isTri; int idx; };
            std::vector<std::vector<IncF>> eInc(edges.size());
            for (int ti = 0; ti < static_cast<int>(remTris.size()); ++ti)
                for (int k = 0; k < 3; ++k) eInc[remTris[ti].e[k]].push_back({true, ti});
            for (int qi = 0; qi < static_cast<int>(joinedQuads.size()); ++qi)
                for (int k = 0; k < 4; ++k) eInc[joinedQuads[qi].e[k]].push_back({false, qi});

            for (size_t ei = 0; ei < edges.size(); ++ei) {
                if (!needSplit[ei] || eInc[ei].empty()) continue;
                auto evalInc = [&]() -> T {
                    T mn = T(180);
                    for (const auto& f : eInc[ei]) mn = std::min(mn, evalFaceMinEf(f.isTri, f.idx));
                    return mn;
                };
                T bestSc = evalInc();
                if (bestSc >= T(26.0)) continue;
                CurveSeg<T> bestH0 = edges[ei].h0, bestH1 = edges[ei].h1;
                gsPoint2D<T> bestMid = edges[ei].mid;

                std::vector<T> uCands = {T(0.18), T(0.28), T(0.38), T(0.62), T(0.72), T(0.82)};
                for (const auto& f : eInc[ei]) {
                    if (f.isTri) {
                        const auto& t = remTris[f.idx];
                        for (int k = 0; k < 3; ++k)
                            if (t.e[k] == static_cast<int>(ei))
                                uCands.push_back(std::max(T(0.14), std::min(T(0.86), edges[ei].seg.closestParam(mesh.pts[t.v[(k + 2) % 3]]))));
                    }
                }
                for (T u : uCands) {
                    auto sp = edges[ei].seg.splitAtParam(u);
                    edges[ei].h0 = sp.first; edges[ei].h1 = sp.second; edges[ei].mid = sp.first.p1();
                    T sc = evalInc();
                    if (sc > bestSc + T(1.0)) { bestSc = sc; bestH0 = edges[ei].h0; bestH1 = edges[ei].h1; bestMid = edges[ei].mid; }
                }
                if (!edges[ei].isBoundary) {
                    const auto& A = mesh.pts[edges[ei].a];
                    const auto& B = mesh.pts[edges[ei].b];
                    T L = A.dist(B);
                    if (L > T(1e-9)) {
                        gsPoint2D<T> perp(-(B.y - A.y) / L, (B.x - A.x) / L);
                        gsPoint2D<T> baseM = bestMid;
                        for (T off : {T(-0.18), T(-0.09), T(0.09), T(0.18)}) {
                            gsPoint2D<T> cM(baseM.x + off * L * perp.x, baseM.y + off * L * perp.y);
                            if (!isPointStrictlyInsideDomain(cM, mesh.loopPolys, T(1e-5)) ||
                                !isSegmentInsideDomain(A, cM, mesh.loopPolys) ||
                                !isSegmentInsideDomain(cM, B, mesh.loopPolys)) continue;
                            edges[ei].mid = cM;
                            edges[ei].h0 = makeStraightSeg(A, cM, degree);
                            edges[ei].h1 = makeStraightSeg(cM, B, degree);
                            T sc = evalInc();
                            if (sc > bestSc + T(1.0)) { bestSc = sc; bestH0 = edges[ei].h0; bestH1 = edges[ei].h1; bestMid = cM; }
                        }
                    }
                }
                edges[ei].h0 = bestH0; edges[ei].h1 = bestH1; edges[ei].mid = bestMid;
            }

            // Build with straight spokes first, then symmetric curved Hermite spokes if needed!
            for (int spokeMode = 0; spokeMode < 2; ++spokeMode) {
                const bool useHermite = (spokeMode == 1);
                auto inwardNormMid = [&](int eIdx, const gsPoint2D<T>& mPt, const gsPoint2D<T>& refIn) -> gsPoint2D<T> {
                    gsPoint2D<T> toIn(refIn.x - mPt.x, refIn.y - mPt.y);
                    T ln = std::hypot(toIn.x, toIn.y);
                    gsPoint2D<T> fb = (ln > T(1e-12)) ? gsPoint2D<T>(toIn.x / ln, toIn.y / ln) : gsPoint2D<T>(1, 0);
                    if (!edges[eIdx].isBoundary || !edges[eIdx].h0.bsp) return fb;
                    gsPoint2D<T> td = evalCurveDeriv(*edges[eIdx].h0.bsp, T(1));
                    T tn = std::hypot(td.x, td.y);
                    if (tn <= T(1e-12)) return fb;
                    gsPoint2D<T> n1(-td.y / tn, td.x / tn), n2(td.y / tn, -td.x / tn);
                    return (n1.x * fb.x + n1.y * fb.y >= n2.x * fb.x + n2.y * fb.y) ? n1 : n2;
                };
                auto mkSpoke = [&](const gsPoint2D<T>& O, const gsPoint2D<T>& mPt, int eIdx, const gsPoint2D<T>& tO) -> CurveSeg<T> {
                    if (!useHermite) return makeStraightSeg(O, mPt, degree);
                    return makeSymmetricHermiteSpoke(O, mPt, tO, inwardNormMid(eIdx, mPt, O), mesh.loopPolys, degree);
                };

                gsMultiPatch<T> mp;
                bool okBuild = true;
                for (const auto& t : remTris) {
                    int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                    CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b; gsPoint2D<T> m0, m1, m2, O;
                    getHalves(t.e[0], v0, v1, e0a, e0b, m0);
                    getHalves(t.e[1], v1, v2, e1a, e1b, m1);
                    getHalves(t.e[2], v2, v0, e2a, e2b, m2);
                    if (!pickOptimalTri1to3CenterCurved(
                            mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], e0a, e0b, e1a, e1b, e2a, e2b,
                            edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O)) {
                        okBuild = false; break;
                    }
                    const T twoPi3 = T(2.09439510239319549);
                    T a0 = std::atan2(m0.y - O.y, m0.x - O.x), a1 = std::atan2(m1.y - O.y, m1.x - O.x), a2 = std::atan2(m2.y - O.y, m2.x - O.x);
                    T sx = std::cos(a0) + std::cos(a1 - twoPi3) + std::cos(a2 - T(2) * twoPi3);
                    T sy = std::sin(a0) + std::sin(a1 - twoPi3) + std::sin(a2 - T(2) * twoPi3);
                    T bA = (std::hypot(sx, sy) > T(1e-9)) ? std::atan2(sy, sx) : a0;
                    auto s0 = mkSpoke(O, m0, t.e[0], gsPoint2D<T>(std::cos(bA), std::sin(bA)));
                    auto s1 = mkSpoke(O, m1, t.e[1], gsPoint2D<T>(std::cos(bA + twoPi3), std::sin(bA + twoPi3)));
                    auto s2 = mkSpoke(O, m2, t.e[2], gsPoint2D<T>(std::cos(bA + T(2) * twoPi3), std::sin(bA + T(2) * twoPi3)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, s0.reversed(), s2, e2b, degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e1a, s1.reversed(), s0, e0b, degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e2a, s2.reversed(), s1, e1b, degree)));
                }
                for (const auto& q : joinedQuads) {
                    if (!okBuild) break;
                    bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                    int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                    if (!sp0 && !sp1) {
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(getSeg(q.e[0], v0, v1), getSeg(q.e[1], v1, v2),
                                                  getSeg(q.e[2], v2, v3), getSeg(q.e[3], v3, v0), degree)));
                    } else if (sp0 && !sp1) {
                        CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m0, m2;
                        getHalves(q.e[0], v0, v1, e0a, e0b, m0);
                        getHalves(q.e[2], v2, v3, e2a, e2b, m2);
                        auto sMid = useHermite
                            ? makeSymmetricHermiteSpoke(m0, m2, inwardNormMid(q.e[0], m0, m2), inwardNormMid(q.e[2], m2, m0), mesh.loopPolys, degree)
                            : makeStraightSeg(m0, m2, degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, sMid, e2b, getSeg(q.e[3], v3, v0), degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0b, getSeg(q.e[1], v1, v2), e2a, sMid.reversed(), degree)));
                    } else if (!sp0 && sp1) {
                        CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m1, m3;
                        getHalves(q.e[1], v1, v2, e1a, e1b, m1);
                        getHalves(q.e[3], v3, v0, e3a, e3b, m3);
                        auto sMid = useHermite
                            ? makeSymmetricHermiteSpoke(m1, m3, inwardNormMid(q.e[1], m1, m3), inwardNormMid(q.e[3], m3, m1), mesh.loopPolys, degree)
                            : makeStraightSeg(m1, m3, degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(getSeg(q.e[0], v0, v1), e1a, sMid, e3b, degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(sMid.reversed(), e1b, getSeg(q.e[2], v2, v3), e3a, degree)));
                    } else {
                        CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m0, m1, m2, m3, O;
                        getHalves(q.e[0], v0, v1, e0a, e0b, m0);
                        getHalves(q.e[1], v1, v2, e1a, e1b, m1);
                        getHalves(q.e[2], v2, v3, e2a, e2b, m2);
                        getHalves(q.e[3], v3, v0, e3a, e3b, m3);
                        if (!chooseQuad4Ctr(q, e0a, e0b, m0, e1a, e1b, m1, e2a, e2b, m2, e3a, e3b, m3, O))
                            O = gsPoint2D<T>((m0.x + m1.x + m2.x + m3.x) * T(0.25), (m0.y + m1.y + m2.y + m3.y) * T(0.25));
                        const T pi2 = T(1.57079632679489662);
                        T b0 = std::atan2(m0.y - O.y, m0.x - O.x), b1 = std::atan2(m1.y - O.y, m1.x - O.x);
                        T b2 = std::atan2(m2.y - O.y, m2.x - O.x), b3 = std::atan2(m3.y - O.y, m3.x - O.x);
                        T qx = std::cos(b0) + std::cos(b1 - pi2) + std::cos(b2 - T(2) * pi2) + std::cos(b3 - T(3) * pi2);
                        T qy = std::sin(b0) + std::sin(b1 - pi2) + std::sin(b2 - T(2) * pi2) + std::sin(b3 - T(3) * pi2);
                        T qA = (std::hypot(qx, qy) > T(1e-9)) ? std::atan2(qy, qx) : b0;
                        auto s0 = mkSpoke(O, m0, q.e[0], gsPoint2D<T>(std::cos(qA), std::sin(qA)));
                        auto s1 = mkSpoke(O, m1, q.e[1], gsPoint2D<T>(std::cos(qA + pi2), std::sin(qA + pi2)));
                        auto s2 = mkSpoke(O, m2, q.e[2], gsPoint2D<T>(std::cos(qA + T(2) * pi2), std::sin(qA + T(2) * pi2)));
                        auto s3 = mkSpoke(O, m3, q.e[3], gsPoint2D<T>(std::cos(qA + T(3) * pi2), std::sin(qA + T(3) * pi2)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, s0.reversed(), s3, e3b, degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e1a, s1.reversed(), s0, e0b, degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e2a, s2.reversed(), s1, e1b, degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e3a, s3.reversed(), s2, e2b, degree)));
                    }
                }
                if (okBuild && acceptCandidateFn(mp, "HarmonicCrossField"))
                    return true;
            }
            return false;
        };

        // Try Matching A (strict + relaxed, with up to 3 bad-quad bans if 1-2 quads fail!)
        std::vector<bool> banned(pairEdges.size(), false);
        for (int attempt = 0; attempt < 3; ++attempt) {
            std::vector<typename TriMeshData<T>::Quad> joinedQuads;
            std::vector<int> quadPairIdx;
            std::vector<typename TriMeshData<T>::Tri> remTris;
            runMatching(true, banned, joinedQuads, quadPairIdx, remTris);

            if (remTris.empty() && !joinedQuads.empty() && static_cast<int>(joinedQuads.size()) <= maxPatchesCap) {
                std::vector<typename TriMeshData<T>::EdgeRec> edges = mesh.edges;
                std::vector<bool> noSplit(edges.size(), false);
                curveInteriorEdgesForQuads(edges, joinedQuads, noSplit);

                auto getSeg = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                    return (edges[eIdx].a == u && edges[eIdx].b == v) ? edges[eIdx].seg : edges[eIdx].seg.reversed();
                };
                gsMultiPatch<T> mp;
                std::vector<int> badPairIds;
                for (size_t qi = 0; qi < joinedQuads.size(); ++qi) {
                    const auto& q = joinedQuads[qi];
                    auto* qp = makeCoonsFromCCWCycle(
                        getSeg(q.e[0], q.v[0], q.v[1]), getSeg(q.e[1], q.v[1], q.v[2]),
                        getSeg(q.e[2], q.v[2], q.v[3]), getSeg(q.e[3], q.v[3], q.v[0]), degree);
                    T angs[4], mnA = 0, mxA = 180; bool deg = false;
                    evalPatchCornerAnglesDeg(*qp, angs, mnA, mxA, deg);
                    if (deg || mnA < effMinAng || mxA > effMaxAng || !isPatchDetSingleSign(*qp, 15, T(1e-8)))
                        badPairIds.push_back(quadPairIdx[qi]);
                    mp.addPatch(typename gsGeometry<T>::uPtr(qp));
                }
                if (acceptCandidateFn(mp, "HarmonicCrossField"))
                    return;
                if (badPairIds.empty() || badPairIds.size() > 3)
                    break;
                for (int pid : badPairIds) banned[pid] = true;
            } else {
                if (emitMixedWithRescues(joinedQuads, remTris))
                    return;
                break;
            }
        }

        // Try Matching B (strict-only pairs -> cleaner angles + zipper/localized split!)
        {
            std::vector<bool> noBan(pairEdges.size(), false);
            std::vector<typename TriMeshData<T>::Quad> joinedQuads;
            std::vector<int> quadPairIdx;
            std::vector<typename TriMeshData<T>::Tri> remTris;
            runMatching(false, noBan, joinedQuads, quadPairIdx, remTris);
            if (!remTris.empty()) {
                if (emitMixedWithRescues(joinedQuads, remTris))
                    return;
            }
        }
    }
}

template <class T>
inline bool buildCrossFieldArrangementCandidate(
    const gsPlanarDomain<T>& dom,
    int degree,
    T stripCollapseRatio,
    T angleThresholdDeg,
    gsMultiPatch<T>& outMp,
    std::string& methodLabel,
    int maxPatchesCap = 140)
{
    const gsMatrix<T> bbox = dom.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    if (diag <= T(1e-12)) return false;
    const T targetH = diag / ((dom.numHoles() >= 4) ? T(26) : T(22));

    TriMeshData<T> fineMesh;
    if (!buildConstrainedDomainTriangulation(dom, degree, targetH, true, 0, fineMesh))
        return false;
    applyDomainLawsonFlips(fineMesh);
    refineInteriorTriangles(fineMesh, targetH, 2);

    CrossFieldMeshAccel<T> accel;
    accel.build(fineMesh);
    CrossFieldSolution<T> cf = solveHarmonic4RoSy(fineMesh, &accel);

    return buildCrossFieldSeparatrixFromSolved(
        dom, fineMesh, cf, accel, degree, stripCollapseRatio, angleThresholdDeg, maxPatchesCap, outMp, methodLabel);
}

} // namespace internal_splitter
} // namespace gismo
