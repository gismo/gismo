/** @file gsPatchSplitterPolychord.hpp
    @brief Global Polychord Optimization, Tangent-Angle Diagnostics, Relaxation,
           and ArgyrisSplines Formatting for gsPatchSplitter<T>.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsSpringPatch.h>
#ifdef gsHLBFGS_ENABLED
#include <gsModeling/gsBarrierPatch.h>
#endif

#include <queue>
#include <numeric>

namespace gismo {
namespace internal_splitter {

// ============================================================================
// Corner Tangent Angle & Aspect Ratio Diagnostics
// ============================================================================

template <class T>
inline void evalPatchCornerAnglesDeg(const gsGeometry<T>& patch,
                                     T anglesDeg[4],
                                     T& minAngDeg,
                                     T& maxAngDeg,
                                     bool& hasDegenerateCorner)
{
    const gsMatrix<T> supp = patch.support();
    const T u0 = supp(0, 0), u1 = supp(0, 1);
    const T v0 = supp(1, 0), v1 = supp(1, 1);
    const T du = u1 - u0, dv = v1 - v0;

    // Corner order: 0=(u0,v0) SW, 1=(u1,v0) SE, 2=(u1,v1) NE, 3=(u0,v1) NW
    // Inward tangent signs (su, sv) so t1 x t2 = +detJ when detJ > 0
    const T cu[4] = {u0, u1, u1, u0};
    const T cv[4] = {v0, v0, v1, v1};
    const bool t1IsU[4] = {true, false, true, false};
    const T s1[4] = {T(1), T(1), T(-1), T(-1)};
    const T s2[4] = {T(1), T(-1), T(-1), T(1)};

    // Batch evaluate center (col 0) + 4 corners (cols 1..4)
    gsMatrix<T> pts(2, 5), dBatch;
    pts(0, 0) = T(0.5) * (u0 + u1); pts(1, 0) = T(0.5) * (v0 + v1);
    for (int c = 0; c < 4; ++c) {
        pts(0, c + 1) = cu[c];
        pts(1, c + 1) = cv[c];
    }
    patch.deriv_into(pts, dBatch);

    const T detMid = dBatch(0, 0) * dBatch(3, 0) - dBatch(2, 0) * dBatch(1, 0);
    const T orientSign = (detMid >= T(0)) ? T(1) : T(-1);

    minAngDeg = T(360);
    maxAngDeg = T(0);
    hasDegenerateCorner = false;
    const T radToDeg = T(180) / T(3.14159265358979323846);

    for (int c = 0; c < 4; ++c) {
        int col = c + 1;
        T ux = dBatch(0, col), vx = dBatch(1, col);
        T uy = dBatch(2, col), vy = dBatch(3, col);
        T t1x = t1IsU[c] ? s1[c] * ux : s1[c] * vx;
        T t1y = t1IsU[c] ? s1[c] * uy : s1[c] * vy;
        T t2x = t1IsU[c] ? s2[c] * vx : s2[c] * ux;
        T t2y = t1IsU[c] ? s2[c] * vy : s2[c] * uy;

        T n1 = std::hypot(t1x, t1y);
        T n2 = std::hypot(t2x, t2y);

        // If a corner control point collapses (zero tangent), probe slightly inward
        if (n1 < T(1e-10) || n2 < T(1e-10)) {
            const T eps = T(1e-3);
            gsMatrix<T> pt(2, 1), der;
            pt(0, 0) = (cu[c] == u0) ? (u0 + eps * du) : (u1 - eps * du);
            pt(1, 0) = (cv[c] == v0) ? (v0 + eps * dv) : (v1 - eps * dv);
            patch.deriv_into(pt, der);
            ux = der(0, 0); vx = der(1, 0);
            uy = der(2, 0); vy = der(3, 0);
            t1x = t1IsU[c] ? s1[c] * ux : s1[c] * vx;
            t1y = t1IsU[c] ? s1[c] * uy : s1[c] * vy;
            t2x = t1IsU[c] ? s2[c] * vx : s2[c] * ux;
            t2y = t1IsU[c] ? s2[c] * vy : s2[c] * uy;
            n1 = std::hypot(t1x, t1y);
            n2 = std::hypot(t2x, t2y);
        }

        if (n1 < T(1e-12) || n2 < T(1e-12)) {
            anglesDeg[c] = T(180);
            hasDegenerateCorner = true;
        } else {
            const T dot = (t1x * t2x + t1y * t2y) / (n1 * n2);
            const T cr  = orientSign * (t1x * t2y - t1y * t2x) / (n1 * n2);
            T ang = std::atan2(cr, dot) * radToDeg;
            if (ang <= T(-1e-6))
                ang += T(360);
            else if (ang < T(0))
                ang = T(0);
            anglesDeg[c] = ang;
            if (ang > T(175) || ang < T(5))
                hasDegenerateCorner = true;
        }
        minAngDeg = std::min(minAngDeg, anglesDeg[c]);
        maxAngDeg = std::max(maxAngDeg, anglesDeg[c]);
    }
}

template <class T>
inline void evalPatchLengthsAndAR(const gsGeometry<T>& patch,
                                  T& Lu,
                                  T& Lv,
                                  T& aspectRatio)
{
    const gsMatrix<T> supp = patch.support();
    const T u0 = supp(0, 0), u1 = supp(0, 1);
    const T v0 = supp(1, 0), v1 = supp(1, 1);
    const int nS = 16;

    auto lineArcLen = [&](bool alongU, T fixedFrac) -> T {
        gsMatrix<T> pts(2, nS + 1), vals;
        for (int i = 0; i <= nS; ++i) {
            T t = static_cast<T>(i) / static_cast<T>(nS);
            if (alongU) {
                pts(0, i) = u0 + t * (u1 - u0);
                pts(1, i) = v0 + fixedFrac * (v1 - v0);
            } else {
                pts(0, i) = u0 + fixedFrac * (u1 - u0);
                pts(1, i) = v0 + t * (v1 - v0);
            }
        }
        patch.eval_into(pts, vals);
        T len = T(0);
        for (int i = 0; i < nS; ++i)
            len += std::hypot(vals(0, i + 1) - vals(0, i), vals(1, i + 1) - vals(1, i));
        return len;
    };

    const T fracs[3] = {T(0), T(0.5), T(1)};
    Lu = T(0);
    Lv = T(0);
    for (int k = 0; k < 3; ++k) {
        Lu += lineArcLen(true, fracs[k]);
        Lv += lineArcLen(false, fracs[k]);
    }
    Lu /= T(3);
    Lv /= T(3);

    const T eps = T(1e-12);
    aspectRatio = std::max(Lu / std::max(Lv, eps), Lv / std::max(Lu, eps));
}

// ============================================================================
// Polychord Graph Decomposition (Union-Find + BFS Orientation Propagation)
// ============================================================================

struct PolychordData {
    int nPatches = 0;
    int nChords = 0;
    // Per half-side id = 4*p + (side.m_index - 1) where west=0, east=1, south=2, north=3
    std::vector<int> sideToChord;
    std::vector<bool> sideFlipped;
    std::vector<bool> chordOrientable;
    // Per patch: chord id for u-split (south/north edges) and v-split (west/east edges)
    std::vector<int> patchChordU;
    std::vector<int> patchChordV;
    std::vector<bool> patchFlipU;
    std::vector<bool> patchFlipV;
    std::vector<std::vector<int>> chordHalfSides;
};

template <class T>
inline PolychordData buildPolychords(const gsMultiPatch<T>& mp)
{
    PolychordData pd;
    pd.nPatches = static_cast<int>(mp.nPatches());
    const int nHalf = 4 * pd.nPatches;
    if (nHalf == 0) return pd;

    struct Link { int to; bool flip; };
    std::vector<std::vector<Link>> adj(nHalf);

    // Opposite sides within each patch: west(0)<->east(1), south(2)<->north(3)
    for (int p = 0; p < pd.nPatches; ++p) {
        int b = 4 * p;
        adj[b + 0].push_back({b + 1, false});
        adj[b + 1].push_back({b + 0, false});
        adj[b + 2].push_back({b + 3, false});
        adj[b + 3].push_back({b + 2, false});
    }

    // Shared interfaces across patches
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
        const boundaryInterface& ifc = *it;
        const patchSide ps1 = ifc.first();
        const patchSide ps2 = ifc.second();
        int id1 = 4 * static_cast<int>(ps1.patch) + (ps1.side().m_index - 1);
        int id2 = 4 * static_cast<int>(ps2.patch) + (ps2.side().m_index - 1);
        const short_t tanDir1 = 1 - ps1.direction();
        const bool flipped = !ifc.dirOrientation(ps1, tanDir1);
        adj[id1].push_back({id2, flipped});
        adj[id2].push_back({id1, flipped});
    }

    pd.sideToChord.assign(nHalf, -1);
    pd.sideFlipped.assign(nHalf, false);

    for (int s = 0; s < nHalf; ++s) {
        if (pd.sideToChord[s] >= 0) continue;
        int cid = pd.nChords++;
        bool orientOk = true;
        std::vector<int> members;

        std::queue<int> q;
        pd.sideToChord[s] = cid;
        pd.sideFlipped[s] = false;
        q.push(s);

        while (!q.empty()) {
            int cur = q.front();
            q.pop();
            members.push_back(cur);
            for (const auto& lk : adj[cur]) {
                bool expectedFlip = pd.sideFlipped[cur] ^ lk.flip;
                if (pd.sideToChord[lk.to] < 0) {
                    pd.sideToChord[lk.to] = cid;
                    pd.sideFlipped[lk.to] = expectedFlip;
                    q.push(lk.to);
                } else if (pd.sideFlipped[lk.to] != expectedFlip) {
                    orientOk = false;
                }
            }
        }
        pd.chordOrientable.push_back(orientOk);
        pd.chordHalfSides.push_back(std::move(members));
    }

    pd.patchChordU.resize(pd.nPatches);
    pd.patchChordV.resize(pd.nPatches);
    pd.patchFlipU.resize(pd.nPatches);
    pd.patchFlipV.resize(pd.nPatches);

    for (int p = 0; p < pd.nPatches; ++p) {
        // u-splitting cuts south(2) & north(3); v-splitting cuts west(0) & east(1)
        pd.patchChordU[p] = pd.sideToChord[4 * p + 2];
        pd.patchFlipU[p]  = pd.sideFlipped[4 * p + 2];
        pd.patchChordV[p] = pd.sideToChord[4 * p + 0];
        pd.patchFlipV[p]  = pd.sideFlipped[4 * p + 0];
    }
    return pd;
}

template <class T>
inline T sideArcLengthFracParam(const gsGeometry<T>& patch, boxSide side, T frac, int nSamples = 64)
{
    const gsMatrix<T> supp = patch.support();
    const short_t nDir = side.direction();
    const short_t tDir = 1 - nDir;
    const T fixedVal = supp(nDir, side.parameter() ? 1 : 0);
    const T t0 = supp(tDir, 0), t1 = supp(tDir, 1);

    gsMatrix<T> pts(2, nSamples + 1), vals;
    for (int i = 0; i <= nSamples; ++i) {
        T t = static_cast<T>(i) / static_cast<T>(nSamples);
        pts(nDir, i) = fixedVal;
        pts(tDir, i) = t0 + t * (t1 - t0);
    }
    patch.eval_into(pts, vals);

    std::vector<T> s(nSamples + 1, T(0));
    for (int i = 1; i <= nSamples; ++i)
        s[i] = s[i - 1] + std::hypot(vals(0, i) - vals(0, i - 1), vals(1, i) - vals(1, i - 1));

    if (s.back() <= T(1e-14)) return frac;
    const T target = frac * s.back();
    for (int i = 1; i <= nSamples; ++i) {
        if (s[i] >= target) {
            T ds = s[i] - s[i - 1];
            T alpha = (ds > T(1e-14)) ? (target - s[i - 1]) / ds : T(0);
            return (static_cast<T>(i - 1) + alpha) / static_cast<T>(nSamples);
        }
    }
    return frac;
}

template <class T>
inline void splitTensorPatchGrid(
    const gsGeometry<T>& geom,
    const std::vector<T>& uSplitsNorm,
    const std::vector<T>& vSplitsNorm,
    gsMultiPatch<T>& outMp)
{
    auto addNormalized = [&](const gsTensorBSpline<2, T>& bsp) {
        auto copy = new gsTensorBSpline<2, T>(bsp);
        copy->knots(0).affineTransformTo(T(0), T(1));
        copy->knots(1).affineTransformTo(T(0), T(1));
        outMp.addPatch(typename gsGeometry<T>::uPtr(copy));
    };
    auto addNormalizedNurbs = [&](const gsTensorNurbs<2, T>& nrb) {
        auto copy = new gsTensorNurbs<2, T>(nrb);
        copy->knots(0).affineTransformTo(T(0), T(1));
        copy->knots(1).affineTransformTo(T(0), T(1));
        outMp.addPatch(typename gsGeometry<T>::uPtr(copy));
    };

    if (const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&geom)) {
        gsTensorBSpline<2, T> base(*tb);
        base.knots(0).affineTransformTo(T(0), T(1));
        base.knots(1).affineTransformTo(T(0), T(1));

        std::vector<gsTensorBSpline<2, T>> uStrips;
        gsTensorBSpline<2, T> remU = base;
        for (T xi : uSplitsNorm) {
            T lo = remU.knots(0).sbegin().value();
            T hi = (--remU.knots(0).send()).value();
            if (xi > lo + T(1e-7) && xi < hi - T(1e-7)) {
                gsTensorBSpline<2, T> left, right;
                remU.splitAt(0, xi, left, right);
                uStrips.push_back(std::move(left));
                remU = std::move(right);
            }
        }
        uStrips.push_back(std::move(remU));

        for (auto& strip : uStrips) {
            gsTensorBSpline<2, T> remV = strip;
            for (T eta : vSplitsNorm) {
                T lo = remV.knots(1).sbegin().value();
                T hi = (--remV.knots(1).send()).value();
                if (eta > lo + T(1e-7) && eta < hi - T(1e-7)) {
                    gsTensorBSpline<2, T> bot, top;
                    remV.splitAt(1, eta, bot, top);
                    addNormalized(bot);
                    remV = std::move(top);
                }
            }
            addNormalized(remV);
        }
        return;
    }

    if (const auto* tn = dynamic_cast<const gsTensorNurbs<2, T>*>(&geom)) {
        gsTensorNurbs<2, T> base(*tn);
        base.knots(0).affineTransformTo(T(0), T(1));
        base.knots(1).affineTransformTo(T(0), T(1));

        std::vector<gsTensorNurbs<2, T>> uStrips;
        gsTensorNurbs<2, T> remU = base;
        for (T xi : uSplitsNorm) {
            T lo = remU.knots(0).sbegin().value();
            T hi = (--remU.knots(0).send()).value();
            if (xi > lo + T(1e-7) && xi < hi - T(1e-7)) {
                gsTensorNurbs<2, T> left, right;
                remU.splitAt(0, xi, left, right);
                uStrips.push_back(std::move(left));
                remU = std::move(right);
            }
        }
        uStrips.push_back(std::move(remU));

        for (auto& strip : uStrips) {
            gsTensorNurbs<2, T> remV = strip;
            for (T eta : vSplitsNorm) {
                T lo = remV.knots(1).sbegin().value();
                T hi = (--remV.knots(1).send()).value();
                if (eta > lo + T(1e-7) && eta < hi - T(1e-7)) {
                    gsTensorNurbs<2, T> bot, top;
                    remV.splitAt(1, eta, bot, top);
                    addNormalizedNurbs(bot);
                    remV = std::move(top);
                }
            }
            addNormalizedNurbs(remV);
        }
        return;
    }

    outMp.addPatch(geom.clone());
}

} // namespace internal_splitter

// ============================================================================
// Stage 2: Global Polychord Aspect-Ratio Optimization
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::optimizePolychordSplits(gsMultiPatch<T>& mp) const
{
    using namespace internal_splitter;
    if (mp.nPatches() == 0) return false;

    mp.computeTopology(T(1e-4), false, false);
    PolychordData pd = buildPolychords(mp);
    const int nP = pd.nPatches;
    const int nC = pd.nChords;
    if (nC == 0) return true;

    std::vector<T> Lu(nP, T(1)), Lv(nP, T(1)), AR(nP, T(1));
    for (int p = 0; p < nP; ++p)
        evalPatchLengthsAndAR(mp.patch(p), Lu[p], Lv[p], AR[p]);

    const T optLambda = static_cast<T>(m_options.getReal("PatchCountWeight"));
    const T lambdaPatch = (m_hasDomain && !m_hasInputMp) ? std::max(optLambda, T(0.55)) : optLambda;
    const T targetAR = std::max(T(0.5), static_cast<T>(m_options.getReal("TargetAspectRatio")));
    const int maxDiv = 16;

    auto evalObjective = [&](const std::vector<int>& n) -> T {
        T energy = T(0);
        for (int p = 0; p < nP; ++p) {
            T hu = Lu[p] / static_cast<T>(n[pd.patchChordU[p]]);
            T hv = Lv[p] / static_cast<T>(n[pd.patchChordV[p]]);
            T r = (hu / std::max(hv, T(1e-12))) / targetAR;
            T invR = T(1) / std::max(r, T(1e-12));
            energy += (r + invR - T(2));
            energy += lambdaPatch * static_cast<T>(n[pd.patchChordU[p]] * n[pd.patchChordV[p]]);
        }
        return energy;
    };

    std::vector<int> nSplit(nC, 1);
    T curE = evalObjective(nSplit);

    // Coordinate descent over integer polychord counts
    for (int iter = 0; iter < 64; ++iter) {
        bool improved = false;
        for (int c = 0; c < nC; ++c) {
            // Try increasing nSplit[c]
            if (nSplit[c] < maxDiv) {
                ++nSplit[c];
                T tryE = evalObjective(nSplit);
                if (tryE < curE - T(1e-10)) {
                    curE = tryE;
                    improved = true;
                    continue;
                }
                --nSplit[c];
            }
            // Try decreasing nSplit[c]
            if (nSplit[c] > 1) {
                --nSplit[c];
                T tryE = evalObjective(nSplit);
                if (tryE < curE - T(1e-10)) {
                    curE = tryE;
                    improved = true;
                    continue;
                }
                ++nSplit[c];
            }
        }
        if (!improved) break;
    }

    bool anySplit = false;
    for (int c = 0; c < nC; ++c)
        if (nSplit[c] > 1) { anySplit = true; break; }
    if (!anySplit) return true;

    // Compute synchronized split parameters in [0,1] for each polychord
    std::vector<std::vector<T>> chordParams(nC);
    for (int c = 0; c < nC; ++c) {
        if (nSplit[c] <= 1) continue;
        chordParams[c].reserve(nSplit[c] - 1);
        if (!pd.chordOrientable[c] || pd.chordHalfSides[c].empty()) {
            for (int k = 1; k < nSplit[c]; ++k)
                chordParams[c].push_back(static_cast<T>(k) / static_cast<T>(nSplit[c]));
        } else {
            // Use mean arc-length parameter across representative side(s) of the polychord
            int repSideId = pd.chordHalfSides[c][0];
            int repP = repSideId / 4;
            boxSide repS = static_cast<boxSide>((repSideId % 4) + 1);
            bool repFlip = pd.sideFlipped[repSideId];
            for (int k = 1; k < nSplit[c]; ++k) {
                T f = static_cast<T>(k) / static_cast<T>(nSplit[c]);
                T fUse = repFlip ? (T(1) - f) : f;
                T uRep = sideArcLengthFracParam(mp.patch(repP), repS, fUse);
                T uCanon = repFlip ? (T(1) - uRep) : uRep;
                chordParams[c].push_back(uCanon);
            }
            std::sort(chordParams[c].begin(), chordParams[c].end());
        }
    }

    auto orientSplits = [](const std::vector<T>& canon, bool flipped) -> std::vector<T> {
        if (!flipped) return canon;
        std::vector<T> rev;
        rev.reserve(canon.size());
        for (auto it = canon.rbegin(); it != canon.rend(); ++it)
            rev.push_back(T(1) - *it);
        return rev;
    };

    gsMultiPatch<T> refinedMp;
    for (int p = 0; p < nP; ++p) {
        std::vector<T> uSplits = orientSplits(chordParams[pd.patchChordU[p]], pd.patchFlipU[p]);
        std::vector<T> vSplits = orientSplits(chordParams[pd.patchChordV[p]], pd.patchFlipV[p]);
        splitTensorPatchGrid(mp.patch(p), uSplits, vSplits, refinedMp);
    }

    refinedMp.computeTopology(T(1e-4), false, false);
    mp = give(refinedMp);
    return true;
}

// ============================================================================
// Stage 2a: Conforming Dual-Sheet / Polychord Strip Reduction (Zero Hanging Nodes)
// ============================================================================

namespace internal_splitter {
struct PolychordSheetPair { int p, sP, q, sQ; };
}

template <class T>
bool gsPatchSplitter<T>::reducePolychordStrips(gsMultiPatch<T>& mp)
{
    using namespace internal_splitter;
    if (mp.nPatches() <= 2) return false;

    const T domMinAng = m_hasDomain ? computeDomainMinCornerAngleDeg(m_domain) : T(90);
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    T curMinA = T(360), curMaxA = T(0);
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T angs[4], mnA = 0, mxA = 180;
        bool degen = false;
        evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
        curMinA = std::min(curMinA, mnA);
        curMaxA = std::max(curMaxA, mxA);
    }
    const T effMinAng = std::min(minAngTol, std::max(T(8), curMinA - T(0.5)));
    const T effMaxAng = std::max(maxAngTol, std::min(T(164.2), curMaxA + T(0.2)));

    bool anyCollapsed = false;

    auto sidePoint = [&](size_t p, int sIdx, T t) -> gsPoint2D<T> {
        const gsMatrix<T> supp = mp.patch(p).support();
        T u0 = supp(0, 0), u1 = supp(0, 1), v0 = supp(1, 0), v1 = supp(1, 1);
        gsMatrix<T> par(2, 1), val;
        if (sIdx == 0)      { par(0, 0) = u0; par(1, 0) = v0 + t * (v1 - v0); }
        else if (sIdx == 1) { par(0, 0) = u1; par(1, 0) = v0 + t * (v1 - v0); }
        else if (sIdx == 2) { par(0, 0) = u0 + t * (u1 - u0); par(1, 0) = v0; }
        else                { par(0, 0) = u0 + t * (u1 - u0); par(1, 0) = v1; }
        mp.patch(p).eval_into(par, val);
        return gsPoint2D<T>(val(0, 0), val(1, 0));
    };

    auto transDirOut = [&](size_t p, int trSide, const gsPoint2D<T>& atPt) -> gsPoint2D<T> {
        gsPoint2D<T> p0 = sidePoint(p, trSide, T(0));
        gsPoint2D<T> p1 = sidePoint(p, trSide, T(1));
        T tA = (atPt.dist(p0) < atPt.dist(p1)) ? T(0.05) : T(0.95);
        gsPoint2D<T> pStep = sidePoint(p, trSide, tA);
        T dx = pStep.x - atPt.x, dy = pStep.y - atPt.y;
        T len = std::hypot(dx, dy);
        return (len > T(1e-14)) ? gsPoint2D<T>(dx / len, dy / len) : gsPoint2D<T>(0, 0);
    };

    auto orientLeft = [](gsTensorBSpline<2, T>& tb, int sideIdx) {
        tb.knots(0).affineTransformTo(T(0), T(1));
        tb.knots(1).affineTransformTo(T(0), T(1));
        if (sideIdx == 1)      tb.reverse(0);
        else if (sideIdx == 2) tb.swapDirections(0, 1);
        else if (sideIdx == 3) { tb.swapDirections(0, 1); tb.reverse(0); }
    };
    auto orientRight = [](gsTensorBSpline<2, T>& tb, int sideIdx) {
        tb.knots(0).affineTransformTo(T(0), T(1));
        tb.knots(1).affineTransformTo(T(0), T(1));
        if (sideIdx == 0)      tb.reverse(0);
        else if (sideIdx == 2) { tb.swapDirections(0, 1); tb.reverse(0); }
        else if (sideIdx == 3) tb.swapDirections(0, 1);
    };

    auto extractSideSeg = [&](const gsTensorBSpline<2, T>& tb, boxSide bs) -> CurveSeg<T> {
        typename gsGeometry<T>::uPtr b = tb.boundary(bs);
        auto* bsp = dynamic_cast<gsBSpline<T>*>(b.get());
        if (bsp)
            return CurveSeg<T>(*bsp);
        gsMatrix<T> p0 = tb.eval(tb.support().col(0)), p1 = tb.eval(tb.support().col(1));
        return makeStraightSeg(gsPoint2D<T>(p0(0, 0), p0(1, 0)), gsPoint2D<T>(p1(0, 0), p1(1, 0)), m_degree);
    };

    auto joinCurveSegsExact = [&](const CurveSeg<T>& c0, const CurveSeg<T>& c1) -> CurveSeg<T> {
        if (!c0.bsp || !c1.bsp) return c0.bsp ? c0 : c1;
        gsBSpline<T> b0(*c0.bsp), b1(*c1.bsp);
        short_t d = std::max(b0.degree(), b1.degree());
        if (b0.degree() < d) b0.degreeElevate(d - b0.degree());
        if (b1.degree() < d) b1.degreeElevate(d - b1.degree());
        b0.knots().affineTransformTo(T(0), T(0.5));
        b1.knots().affineTransformTo(T(0.5), T(1));

        const int n0 = b0.coefs().rows();
        const int n1 = b1.coefs().rows();
        std::vector<T> kts;
        kts.reserve(static_cast<size_t>(n0 + n1 + d));
        const auto& k0 = b0.knots();
        const auto& k1 = b1.knots();
        for (int i = 0; i < static_cast<int>(k0.size()) - (d + 1); ++i)
            kts.push_back(k0[i]);
        for (int r = 0; r < d; ++r)
            kts.push_back(T(0.5));
        for (int i = d + 1; i < static_cast<int>(k1.size()); ++i)
            kts.push_back(k1[i]);

        gsMatrix<T> cfs(n0 + n1 - 1, 2);
        cfs.topRows(n0) = b0.coefs();
        cfs.bottomRows(n1 - 1) = b1.coefs().bottomRows(n1 - 1);
        cfs.row(n0 - 1) = T(0.5) * (b0.coefs().bottomRows(1) + b1.coefs().topRows(1));

        gsKnotVector<T> kv(kts, static_cast<short_t>(d));
        return CurveSeg<T>(gsBSpline<T>(kv, cfs));
    };

    for (int outerPass = 0; outerPass < 32 && mp.nPatches() > 2; ++outerPass) {
        mp.computeTopology(T(1e-4), false, false);
        const int nP = static_cast<int>(mp.nPatches());
        const auto& ifcs = mp.topology().interfaces();
        if (ifcs.empty()) break;

        // ifcAtSide[4*p + s] = interface index, or -1 if boundary
        std::vector<int> ifcAtSide(4 * nP, -1);
        for (int ii = 0; ii < static_cast<int>(ifcs.size()); ++ii) {
            const auto& ps1 = ifcs[ii].first();
            const auto& ps2 = ifcs[ii].second();
            ifcAtSide[4 * static_cast<int>(ps1.patch) + (ps1.side().m_index - 1)] = ii;
            ifcAtSide[4 * static_cast<int>(ps2.patch) + (ps2.side().m_index - 1)] = ii;
        }

        std::vector<bool> triedIfc(ifcs.size(), false);
        bool collapsedOne = false;

        for (int startIi = 0; startIi < static_cast<int>(ifcs.size()); ++startIi) {
            if (triedIfc[startIi]) continue;

            std::vector<PolychordSheetPair> pairs;
            std::map<int, int> patchRole; // patch -> pairIdx
            std::queue<int> qPairs;

            const auto& sIfc = ifcs[startIi];
            int p0 = static_cast<int>(sIfc.first().patch), sP0 = sIfc.first().side().m_index - 1;
            int q0 = static_cast<int>(sIfc.second().patch), sQ0 = sIfc.second().side().m_index - 1;
            if (p0 == q0) continue;

            pairs.push_back({p0, sP0, q0, sQ0});
            patchRole[p0] = 0;
            patchRole[q0] = 0;
            qPairs.push(0);
            triedIfc[startIi] = true;

            bool sheetOk = true;
            while (!qPairs.empty() && sheetOk) {
                const PolychordSheetPair pr = pairs[qPairs.front()];
                qPairs.pop();

                int trP[2] = { (pr.sP < 2) ? 2 : 0, (pr.sP < 2) ? 3 : 1 };
                int trQ[2] = { (pr.sQ < 2) ? 2 : 0, (pr.sQ < 2) ? 3 : 1 };

                for (int endIdx = 0; endIdx < 2 && sheetOk; ++endIdx) {
                    int sTp = trP[endIdx];
                    gsPoint2D<T> cornerP0 = sidePoint(pr.p, pr.sP, T(endIdx));
                    // Match transverse side on q sharing cornerP0
                    int sTq = -1;
                    for (int candIdx = 0; candIdx < 2; ++candIdx) {
                        int candQ = trQ[candIdx];
                        if (cornerP0.dist(sidePoint(pr.q, candQ, T(0))) < T(1e-3) ||
                            cornerP0.dist(sidePoint(pr.q, candQ, T(1))) < T(1e-3))
                        {
                            sTq = candQ;
                            break;
                        }
                    }
                    if (sTq < 0) { sheetOk = false; break; }

                    int iiP = ifcAtSide[4 * pr.p + sTp];
                    int iiQ = ifcAtSide[4 * pr.q + sTq];
                    if ((iiP < 0) != (iiQ < 0)) { sheetOk = false; break; }

                    if (iiP < 0 && iiQ < 0) {
                        // Both boundary: ensure smooth (collinear) join >= 162 deg, never erase a corner
                        gsPoint2D<T> dP = transDirOut(pr.p, sTp, cornerP0);
                        gsPoint2D<T> dQ = transDirOut(pr.q, sTq, cornerP0);
                        if (dP.x * dQ.x + dP.y * dQ.y > -T(0.951)) {
                            sheetOk = false;
                            break;
                        }
                        continue;
                    }

                    // Both interior: neighbors must share a matching interface through cornerP0
                    const auto& ifP = ifcs[iiP];
                    int pNext = (static_cast<int>(ifP.first().patch) == pr.p && (ifP.first().side().m_index - 1) == sTp)
                              ? static_cast<int>(ifP.second().patch) : static_cast<int>(ifP.first().patch);
                    const auto& ifQ = ifcs[iiQ];
                    int qNext = (static_cast<int>(ifQ.first().patch) == pr.q && (ifQ.first().side().m_index - 1) == sTq)
                              ? static_cast<int>(ifQ.second().patch) : static_cast<int>(ifQ.first().patch);
                    if (pNext == qNext || pNext == pr.p || qNext == pr.q) { sheetOk = false; break; }

                    int nextIi = -1, sNextP = -1, sNextQ = -1;
                    for (int spCand = 0; spCand < 4; ++spCand) {
                        int candIi = ifcAtSide[4 * pNext + spCand];
                        if (candIi < 0) continue;
                        const auto& candIf = ifcs[candIi];
                        int otherP = (static_cast<int>(candIf.first().patch) == pNext && (candIf.first().side().m_index - 1) == spCand)
                                   ? static_cast<int>(candIf.second().patch) : static_cast<int>(candIf.first().patch);
                        int otherS = (static_cast<int>(candIf.first().patch) == pNext && (candIf.first().side().m_index - 1) == spCand)
                                   ? (candIf.second().side().m_index - 1) : (candIf.first().side().m_index - 1);
                        if (otherP != qNext) continue;
                        if (cornerP0.dist(sidePoint(pNext, spCand, T(0))) < T(1e-3) ||
                            cornerP0.dist(sidePoint(pNext, spCand, T(1))) < T(1e-3))
                        {
                            nextIi = candIi; sNextP = spCand; sNextQ = otherS;
                            break;
                        }
                    }
                    if (nextIi < 0) { sheetOk = false; break; }
                    triedIfc[nextIi] = true;

                    auto itP = patchRole.find(pNext);
                    auto itQ = patchRole.find(qNext);
                    if (itP != patchRole.end() || itQ != patchRole.end()) {
                        if (itP == patchRole.end() || itQ == patchRole.end() || itP->second != itQ->second) {
                            sheetOk = false; break;
                        }
                        const auto& ex = pairs[itP->second];
                        if (!((ex.p == pNext && ex.q == qNext && ex.sP == sNextP && ex.sQ == sNextQ) ||
                              (ex.p == qNext && ex.q == pNext && ex.sP == sNextQ && ex.sQ == sNextP)))
                        {
                            sheetOk = false; break;
                        }
                    } else {
                        int newIdx = static_cast<int>(pairs.size());
                        pairs.push_back({pNext, sNextP, qNext, sNextQ});
                        patchRole[pNext] = newIdx;
                        patchRole[qNext] = newIdx;
                        qPairs.push(newIdx);
                    }
                }
            }

            if (!sheetOk || pairs.empty()) continue;

            // Construct merged candidate multi-patch
            gsMultiPatch<T> candMp;
            for (int p = 0; p < nP; ++p) {
                if (patchRole.find(p) == patchRole.end())
                    candMp.addPatch(mp.patch(p).clone());
            }

            bool buildOk = true;
            for (const auto& pr : pairs) {
                const auto* tbP = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(pr.p));
                const auto* tbQ = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(pr.q));
                if (!tbP || !tbQ) { buildOk = false; break; }

                gsTensorBSpline<2, T> A(*tbP), B(*tbQ);
                orientRight(A, pr.sP);
                orientLeft(B, pr.sQ);

                gsMatrix<T> pA0(2, 1), pB0(2, 1), pB1(2, 1);
                pA0(0, 0) = T(1); pA0(1, 0) = T(0);
                pB0(0, 0) = T(0); pB0(1, 0) = T(0);
                pB1(0, 0) = T(0); pB1(1, 0) = T(1);
                gsMatrix<T> vA0, vB0, vB1;
                A.eval_into(pA0, vA0);
                B.eval_into(pB0, vB0);
                B.eval_into(pB1, vB1);
                if ((vA0 - vB1).norm() < (vA0 - vB0).norm())
                    B.reverse(1);

                CurveSeg<T> cS = joinCurveSegsExact(extractSideSeg(A, boundary::south), extractSideSeg(B, boundary::south));
                CurveSeg<T> cE = extractSideSeg(B, boundary::east);
                CurveSeg<T> cN = joinCurveSegsExact(extractSideSeg(A, boundary::north), extractSideSeg(B, boundary::north)).reversed();
                CurveSeg<T> cW = extractSideSeg(A, boundary::west).reversed();

                gsTensorBSpline<2, T>* merged = makeCoonsFromCCWCycle(cS, cE, cN, cW, m_degree);
                if (!merged) { buildOk = false; break; }
                candMp.addPatch(typename gsGeometry<T>::uPtr(merged));
            }
            if (!buildOk || candMp.nPatches() >= mp.nPatches()) continue;

            // Quick local angle + Jacobian filter on merged patches before full validation
            gsASG1LocalFitter<T>::unifyPatchOrientations(candMp);
            optimizeWinslowBarrier(candMp);
            candMp.computeTopology(T(1e-4), false, false);
            if (2 * candMp.nInterfaces() + candMp.nBoundary() != 4 * candMp.nPatches()) continue;

            bool localOk = true;
            const size_t nOrigKept = static_cast<size_t>(nP) - 2 * pairs.size();
            for (size_t p = nOrigKept; p < candMp.nPatches(); ++p) {
                T angs[4], mnA = 0, mxA = 180; bool degen = false;
                evalPatchCornerAnglesDeg(candMp.patch(p), angs, mnA, mxA, degen);
                if (degen || mnA < minAngTol || mxA > maxAngTol) { localOk = false; break; }
                T md = 0, mxd = 0, ma = 0, en = 0;
                evalPatchJacobianAndWinslow(candMp.patch(p), 13, md, mxd, ma, en);
                if (md <= T(1e-8) && mxd >= -T(1e-8)) { localOk = false; break; }
            }
            if (!localOk) continue;

            if (m_hasDomain && !validateMultiPatch(m_domain, candMp, T(1e-4), nullptr, effMinAng, effMaxAng)) {
                if (candMp.nPatches() <= 300) {
                    optimizeWinslowBarrier(candMp);
                    candMp.computeTopology(T(1e-4), false, false);
                    if (!validateMultiPatch(m_domain, candMp, T(1e-4), nullptr, effMinAng, effMaxAng))
                        continue;
                } else {
                    continue;
                }
            }

            mp = give(candMp);
            anyCollapsed = true;
            collapsedOne = true;
            break;
        }
        if (!collapsedOne) break;
    }
    return anyCollapsed;
}

// ============================================================================
// Knot Vector Harmonization Across Polychords
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::harmonizeMultiPatchKnots(gsMultiPatch<T>& mp) const
{
    using namespace internal_splitter;
    if (mp.nPatches() == 0) return false;
    mp.computeTopology(T(1e-4), false, false);

    // Ensure [0,1]^2 parameter domain on all B-spline patches
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
            tb->knots(0).affineTransformTo(T(0), T(1));
            tb->knots(1).affineTransformTo(T(0), T(1));
        }
    }

    PolychordData pd = buildPolychords(mp);
    for (int c = 0; c < pd.nChords; ++c) {
        struct Entry { int patch; int dir; bool flip; };
        std::vector<Entry> entries;
        for (int p = 0; p < pd.nPatches; ++p) {
            if (pd.patchChordU[p] == c) entries.push_back({p, 0, pd.patchFlipU[p]});
            if (pd.patchChordV[p] == c) entries.push_back({p, 1, pd.patchFlipV[p]});
        }
        if (entries.size() <= 1) continue;

        short_t maxDeg = 1;
        for (const auto& e : entries) {
            if (const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(e.patch)))
                maxDeg = std::max(maxDeg, tb->degree(e.dir));
        }

        // Elevate degrees along this polychord
        for (const auto& e : entries) {
            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(e.patch))) {
                if (tb->degree(e.dir) < maxDeg)
                    tb->degreeElevate(maxDeg - tb->degree(e.dir), e.dir);
            }
        }

        // Compute canonical union of internal knots along the polychord with 1e-10 tolerance deduplication
        const T knotTol = T(1e-10);
        struct KnotCluster { T u; int mult; };
        std::vector<KnotCluster> canonKnots;
        bool initUnion = false;
        for (const auto& e : entries) {
            if (const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(e.patch))) {
                initUnion = true;
                gsKnotVector<T> kv = tb->knots(e.dir);
                if (pd.chordOrientable[c] && e.flip) kv.reverse();
                kv.affineTransformTo(T(0), T(1));
                std::vector<KnotCluster> localClusters;
                for (index_t ki = 0; ki < kv.size(); ++ki) {
                    T u = kv[ki];
                    if (u <= knotTol || u >= T(1) - knotTol) continue;
                    bool found = false;
                    for (auto& lc : localClusters) {
                        if (std::abs(lc.u - u) < knotTol) { ++lc.mult; found = true; break; }
                    }
                    if (!found) localClusters.push_back({u, 1});
                }
                for (const auto& lc : localClusters) {
                    int clampedMult = std::min<int>(lc.mult, maxDeg);
                    bool found = false;
                    for (auto& cc : canonKnots) {
                        if (std::abs(cc.u - lc.u) < knotTol) {
                            cc.mult = std::max(cc.mult, clampedMult);
                            found = true;
                            break;
                        }
                    }
                    if (!found) canonKnots.push_back({lc.u, clampedMult});
                }
            }
        }
        if (!initUnion || canonKnots.empty()) continue;

        // Insert missing knots into each patch along this polychord
        for (const auto& e : entries) {
            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(e.patch))) {
                const auto& kv = tb->knots(e.dir);
                for (const auto& cc : canonKnots) {
                    T targetU = (pd.chordOrientable[c] && e.flip) ? (T(1) - cc.u) : cc.u;
                    if (targetU <= knotTol || targetU >= T(1) - knotTol) continue;
                    int curMult = 0;
                    for (index_t ki = 0; ki < kv.size(); ++ki)
                        if (std::abs(kv[ki] - targetU) < knotTol) ++curMult;
                    int addMult = cc.mult - curMult;
                    if (addMult > 0)
                        tb->insertKnot(targetU, e.dir, addMult);
                }
            }
        }
    }

    mp.computeTopology(T(1e-4), false, false);
    return true;
}

// ============================================================================
// Relaxation & Native G+Smo Parameterization
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::relaxAndParameterize(gsMultiPatch<T>& mp, bool runBarrier)
{
    if (mp.nPatches() == 0) return false;
    harmonizeMultiPatchKnots(mp);

    // Snapshot pre-barrier valid state so relaxation never introduces a negative/zero Jacobian fold
    bool hadValidDet = (mp.nPatches() > 0);
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T mn = 0, mx = 0, ma = 0, en = 0;
        internal_splitter::evalPatchJacobianAndWinslow(mp.patch(p), 17, mn, mx, ma, en);
        if (mn <= T(1e-8) && mx >= -T(1e-8)) { hadValidDet = false; break; }
    }
    gsMultiPatch<T> preRelaxMp;
    if (hadValidDet) preRelaxMp = mp;

    // Check if any patch needs SpringPatch recovery (e.g. folded Coons interior)
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T minDet = 0, maxDet = 0, minAbs = 0, en = 0;
        internal_splitter::evalPatchJacobianAndWinslow(mp.patch(p), 17, minDet, maxDet, minAbs, en);
        if (minDet <= T(1e-7) && maxDet >= -T(1e-7)) {
            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
                // Try native gsSpringPatch from the 4 boundary curves
                gsMultiPatch<T> bndCurves;
                bndCurves.addPatch(tb->boundary(boundary::south));
                bndCurves.addPatch(tb->boundary(boundary::east));
                bndCurves.addPatch(tb->boundary(boundary::north));
                bndCurves.addPatch(tb->boundary(boundary::west));
                try {
                    gsSpringPatch<T> spring(bndCurves);
                    const gsGeometry<T>& res = spring.compute();
                    T sMinDet = 0, sMaxDet = 0, sMinAbs = 0, sEn = 0;
                    internal_splitter::evalPatchJacobianAndWinslow(res, 17, sMinDet, sMaxDet, sMinAbs, sEn);
                    if ((sMinDet > T(1e-7)) || (sMaxDet < -T(1e-7)) || (sMinDet > minDet)) {
                        if (const auto* sTb = dynamic_cast<const gsTensorBSpline<2, T>*>(&res)) {
                            *tb = *sTb;
                        }
                    }
                } catch (...) {
                    // Keep existing patch if spring construction fails
                }
            }
        }
    }

    if (runBarrier) {
#ifdef gsHLBFGS_ENABLED
        if (m_options.getInt("ParameterizationMethod") == 1) {
            try {
                gsBarrierPatch<2, T> bp(mp, true);
                bp.compute();
                mp = bp.result();
            } catch (...) {
                optimizeWinslowBarrier(mp);
            }
        } else {
            optimizeWinslowBarrier(mp);
        }
#else
        optimizeWinslowBarrier(mp);
#endif
    }

    if (hadValidDet) {
        bool brokeDet = false;
        for (size_t p = 0; p < mp.nPatches(); ++p) {
            T mn = 0, mx = 0, ma = 0, en = 0;
            internal_splitter::evalPatchJacobianAndWinslow(mp.patch(p), 21, mn, mx, ma, en);
            if (mn <= T(1e-8) && mx >= -T(1e-8)) { brokeDet = true; break; }
        }
        if (brokeDet) mp = give(preRelaxMp);
    }

    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(T(1e-4), false, false);
    return true;
}

// ============================================================================
// ArgyrisSplines Post-Processing & Formatting
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::formatForArgyris(gsMultiPatch<T>& mp) const
{
    if (mp.nPatches() == 0) return false;

    const short_t targetDeg = static_cast<short_t>(std::max(1, m_options.getInt("Degree")));
    const short_t targetReg = static_cast<short_t>(std::max(1, std::min<int>(targetDeg - 2, m_options.getInt("Regularity"))));
    const int minSpans = std::max(1, m_options.getInt("MinKnotSpans"));
    const int mult = std::max<int>(1, targetDeg - targetReg);

    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(T(1e-4), false, false);

    // Elevate degree uniformly if below targetDeg
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
            tb->knots(0).affineTransformTo(T(0), T(1));
            tb->knots(1).affineTransformTo(T(0), T(1));
            for (short_t d = 0; d < 2; ++d) {
                if (tb->degree(d) < targetDeg)
                    tb->degreeElevate(targetDeg - tb->degree(d), d);
            }
        }
    }

    harmonizeMultiPatchKnots(mp);

    // Ensure minimum number of knot spans along every patch of every polychord
    if (minSpans > 1) {
        internal_splitter::PolychordData pd = internal_splitter::buildPolychords(mp);
        std::vector<int> chordMinSpans(pd.nChords, 1000000);
        for (int p = 0; p < pd.nPatches; ++p) {
            if (const auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p))) {
                int su = static_cast<int>(tb->knots(0).unique().size()) - 1;
                int sv = static_cast<int>(tb->knots(1).unique().size()) - 1;
                chordMinSpans[pd.patchChordU[p]] = std::min(chordMinSpans[pd.patchChordU[p]], su);
                chordMinSpans[pd.patchChordV[p]] = std::min(chordMinSpans[pd.patchChordV[p]], sv);
            }
        }

        bool insertedAny = false;
        for (int c = 0; c < pd.nChords; ++c) {
            int curSpans = chordMinSpans[c];
            if (curSpans >= minSpans) continue;
            insertedAny = true;

            // Compute uniform knots to insert along chord c so every patch gets >= minSpans
            std::vector<T> targetKnots;
            for (int s = 1; s < minSpans; ++s) {
                T u = static_cast<T>(s) / static_cast<T>(minSpans);
                targetKnots.push_back(u);
            }

            for (int p = 0; p < pd.nPatches; ++p) {
                if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p))) {
                    for (int dir = 0; dir < 2; ++dir) {
                        int pc = (dir == 0) ? pd.patchChordU[p] : pd.patchChordV[p];
                        if (pc != c) continue;
                        for (T u : targetKnots) {
                            const auto& kv = tb->knots(dir);
                            int curMult = 0;
                            for (index_t ki = 0; ki < kv.size(); ++ki)
                                if (std::abs(kv[ki] - u) < T(1e-10)) ++curMult;
                            int addMult = mult - curMult;
                            if (addMult > 0)
                                tb->insertKnot(u, dir, addMult);
                        }
                    }
                }
            }
        }
        if (insertedAny)
            harmonizeMultiPatchKnots(mp);
    }

    if (m_options.getSwitch("ApplyASG1Fitting") && targetDeg >= 3) {
        const index_t innerKnots = std::max<index_t>(2, minSpans - 1);
        auto tryFitMode = [&](typename gsASG1LocalFitter<T>::FittingMode mode, index_t nInner) -> bool {
            gsASG1LocalFitter<T> fitter(mp, targetDeg, targetReg, nInner);
            fitter.setMode(mode);
            fitter.fit();
            if (!fitter.checkASG1(T(1e-3)))
                return false;
            gsMultiPatch<T> fitted = fitter.result();
            gsASG1LocalFitter<T>::unifyPatchOrientations(fitted);
            fitted.computeTopology(T(1e-4), false, false);
            for (size_t p = 0; p < fitted.nPatches(); ++p) {
                if (auto* tb = dynamic_cast<const gsTensorBSpline<2, T>*>(&fitted.patch(p))) {
                    T mnDet = 0, mxDet = 0, mnAbs = 0, en = 0;
                    internal_splitter::evalPatchJacobianAndWinslow(*tb, 21, mnDet, mxDet, mnAbs, en);
                    if ((mnDet <= T(1e-9) && mxDet >= -T(1e-9)) ||
                        !internal_splitter::isPatchDetSingleSign(*tb, 21, T(1e-7)) ||
                        !internal_splitter::isPatchDetSingleSign(*tb, 20, T(1e-7)))
                        return false;
                    T angs[4], minA = 0, maxA = 180;
                    bool degen = false;
                    internal_splitter::evalPatchCornerAnglesDeg(*tb, angs, minA, maxA, degen);
                    if (degen || minA < T(4.0) || maxA > T(164.9))
                        return false;
                }
            }
            // Reject fit if boundary patch corners drifted significantly
            if (fitted.nPatches() == mp.nPatches()) {
                for (const auto& bs : mp.topology().boundaries()) {
                    const size_t p = static_cast<size_t>(bs.patch);
                    gsMatrix<T> par(2, 2);
                    if (bs.side() == boundary::west)      { par << T(0), T(0), T(0), T(1); }
                    else if (bs.side() == boundary::east) { par << T(1), T(1), T(0), T(1); }
                    else if (bs.side() == boundary::south){ par << T(0), T(1), T(0), T(0); }
                    else                                  { par << T(0), T(1), T(1), T(1); }
                    gsMatrix<T> e0, e1;
                    mp.patch(p).eval_into(par, e0);
                    fitted.patch(p).eval_into(par, e1);
                    T bDiff = (e0 - e1).colwise().norm().maxCoeff();
                    if (bDiff > T(1.8e-3))
                        return false;
                }
            }
            if (m_hasDomain && !verifyBoundaryMatch(m_domain, fitted, T(1.9e-3)))
                return false;
            mp = give(fitted);
            return true;
        };
        auto runWithMode = [&](typename gsASG1LocalFitter<T>::FittingMode mode,
                               std::initializer_list<index_t> extras) -> bool {
            for (index_t extra : extras) {
                const index_t nInner = innerKnots + extra;
                const index_t n1D = 4 + 2 * nInner;
                if (mode == gsASG1LocalFitter<T>::FittingMode::Global &&
                    static_cast<index_t>(mp.nPatches()) * n1D * n1D > 3200)
                    continue;
                if (tryFitMode(mode, nInner))
                    return true;
            }
            return false;
        };
        const std::initializer_list<index_t> allExtras = {index_t(0), index_t(1), index_t(2), index_t(3), index_t(4), index_t(5), index_t(6), index_t(8), index_t(14), index_t(22)};
        if (mp.nPatches() <= 20) {
            if (!runWithMode(gsASG1LocalFitter<T>::FittingMode::Global, allExtras))
                runWithMode(gsASG1LocalFitter<T>::FittingMode::Local, allExtras);
        } else if (mp.nPatches() <= 50) {
            if (!runWithMode(gsASG1LocalFitter<T>::FittingMode::Local, {index_t(0), index_t(1), index_t(2), index_t(3), index_t(4)}) &&
                !runWithMode(gsASG1LocalFitter<T>::FittingMode::Global, {index_t(0), index_t(1), index_t(2)}))
                runWithMode(gsASG1LocalFitter<T>::FittingMode::Local, {index_t(5), index_t(6), index_t(8), index_t(14), index_t(22)});
        } else {
            runWithMode(gsASG1LocalFitter<T>::FittingMode::Local, allExtras);
        }
    }

    return true;
}

// ============================================================================
// Quality Report Diagnostics
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::QualityReport
gsPatchSplitter<T>::analyzeQuality(const gsMultiPatch<T>& result) const
{
    QualityReport rep;
    rep.nPatches = static_cast<int>(result.nPatches());
    if (rep.nPatches == 0) return rep;

    gsMultiPatch<T> mpCopy(result);
    mpCopy.computeTopology(T(1e-4), false, false);
    rep.nInterfaces = static_cast<int>(mpCopy.nInterfaces());
    rep.nBoundarySides = static_cast<int>(mpCopy.nBoundary());

    // Check for unmatched interior sides (hanging nodes / T-junctions)
    int accountedSides = 2 * rep.nInterfaces + rep.nBoundarySides;
    rep.hasHangingNodes = (accountedSides != 4 * rep.nPatches);

    rep.isConforming = !rep.hasHangingNodes;
    gsMatrix<T> uMat(1, 1), pt1, pt2;
    for (const auto& ifc : mpCopy.interfaces()) {
        auto b1 = mpCopy.patch(ifc.first().patch).boundary(ifc.first().side());
        auto b2 = mpCopy.patch(ifc.second().patch).boundary(ifc.second().side());
        gsMatrix<T> pr1 = b1->parameterRange(), pr2 = b2->parameterRange();
        uMat(0, 0) = T(0.5) * (pr1(0, 0) + pr1(0, 1));
        b1->eval_into(uMat, pt1);
        uMat(0, 0) = T(0.5) * (pr2(0, 0) + pr2(0, 1));
        b2->eval_into(uMat, pt2);
        if (std::hypot(pt1(0, 0) - pt2(0, 0), pt1(1, 0) - pt2(1, 0)) > T(1e-4))
            rep.isConforming = false;
    }

    // Interior vertex valences via gsBoxTopology corner grouping
    std::vector<std::vector<patchCorner>> cornerGroups;
    for (size_t p = 0; p < mpCopy.nPatches(); ++p) {
        for (boxCorner c = boxCorner::getFirst(2); c < boxCorner::getEnd(2); ++c) {
            patchCorner pc(p, c);
            std::vector<patchCorner> grp;
            mpCopy.getCornerList(pc, grp);
            std::sort(grp.begin(), grp.end());
            if (std::find(cornerGroups.begin(), cornerGroups.end(), grp) == cornerGroups.end())
                cornerGroups.push_back(grp);
        }
    }
    for (const auto& grp : cornerGroups) {
        bool isBnd = false;
        for (const auto& pc : grp) {
            std::vector<patchSide> sides;
            pc.getContainingSides(2, sides);
            for (const auto& ps : sides)
                if (mpCopy.isBoundary(ps)) { isBnd = true; break; }
            if (isBnd) break;
        }
        if (!isBnd) {
            int val = static_cast<int>(grp.size());
            rep.intValences.push_back(val);
            if (val != 4) ++rep.nExtraordinaryVertices;
        }
    }

    rep.minDetJ = std::numeric_limits<T>::max();
    rep.maxDetJ = -std::numeric_limits<T>::max();
    rep.minAspectRatio = std::numeric_limits<T>::max();
    rep.maxAspectRatio = T(0);
    rep.minAngleDeg = T(360);
    rep.maxAngleDeg = T(0);

    for (size_t p = 0; p < mpCopy.nPatches(); ++p) {
        T pMinDet = 0, pMaxDet = 0, pMinAbs = 0, pEn = 0;
        internal_splitter::evalPatchJacobianAndWinslow(mpCopy.patch(p), 21, pMinDet, pMaxDet, pMinAbs, pEn);
        rep.minDetJ = std::min(rep.minDetJ, pMinDet);
        rep.maxDetJ = std::max(rep.maxDetJ, pMaxDet);

        T lu = 0, lv = 0, ar = 1;
        internal_splitter::evalPatchLengthsAndAR(mpCopy.patch(p), lu, lv, ar);
        rep.minAspectRatio = std::min(rep.minAspectRatio, ar);
        rep.maxAspectRatio = std::max(rep.maxAspectRatio, ar);

        T angs[4], minA = 0, maxA = 0;
        bool degen = false;
        internal_splitter::evalPatchCornerAnglesDeg(mpCopy.patch(p), angs, minA, maxA, degen);
        rep.minAngleDeg = std::min(rep.minAngleDeg, minA);
        rep.maxAngleDeg = std::max(rep.maxAngleDeg, maxA);
    }

    rep.isValid = rep.isConforming && !rep.hasHangingNodes && (rep.minDetJ > T(1e-9)) && (rep.maxAngleDeg < T(165));
    return rep;
}

template <class T>
void gsPatchSplitter<T>::printQualityReport(const gsMultiPatch<T>& result) const
{
    QualityReport r = analyzeQuality(result);
    gsInfo << "=== gsPatchSplitter Quality Report ===\n"
           << "  Patches:               " << r.nPatches << "\n"
           << "  Interfaces / Boundary: " << r.nInterfaces << " / " << r.nBoundarySides << "\n"
           << "  Hanging Nodes:         " << (r.hasHangingNodes ? "YES (FAILED)" : "NO (0 T-junctions)") << "\n"
           << "  Conforming Topology:   " << (r.isConforming ? "YES" : "NO") << "\n"
           << "  Extraordinary Vertices:" << r.nExtraordinaryVertices << " (interior valences:";
    if (r.intValences.empty()) gsInfo << " none";
    else for (int v : r.intValences) gsInfo << " " << v;
    gsInfo << ")\n"
           << "  Min / Max det(J):      " << r.minDetJ << " / " << r.maxDetJ << "\n"
           << "  Aspect Ratio Range:    [" << r.minAspectRatio << ", " << r.maxAspectRatio << "]\n"
           << "  Corner Angle Range:    [" << r.minAngleDeg << " deg, " << r.maxAngleDeg << " deg]\n"
           << "  Overall Valid:         " << (r.isValid ? "PASS" : "FAIL") << "\n"
           << "======================================\n";
}

} // namespace gismo
