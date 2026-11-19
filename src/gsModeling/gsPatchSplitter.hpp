/** @file gsPatchSplitter.hpp
    @brief Implementation of gsPatchSplitter<T>: Planar Quadrilateral Domain Decomposition.

    Decomposes 2D planar domains (gsPlanarDomain) into conforming multi-patch
    tensor B-spline geometries (gsMultiPatch) suitable for isogeometric analysis
    and downstream AS-G1 reparametrization (gsASG1LocalFitter / gsAsG1GeometryGenerator).

    Key Properties & Implemented Pipeline:
    1. Exact CAD boundary preservation: outer boundary curves (and hole curves)
       are retained verbatim or subdivided via exact B-spline splitting (splitAt)
       and knot harmonization without perturbing boundary control points.
    2. Hole-free decomposition: evaluates structural quad layouts (convex 1-patch
       preservation, O-grid, ribbon, reflex L-cut, star fans/kites, comb spine-teeth)
       and triangulation-based 1-to-3 Catmull-Clark quadification fallback, selecting
       the lowest-cost valid candidate by patch count and distortion energy.
    3. Domains with holes: quad-ring and block-ring decomposition templates for
       single-hole and supported multi-hole topologies.
    4. Parameterization & validation: Coons blending with interior-only, barrier-guarded
       smoothing; checks conformity, exact boundary match, per-patch nonzero
       single-sign Jacobian determinant, strictly convex patch corners, and absence
       of self-intersections or patch overlaps.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsModeling/gsPatchSplitter.h>
#include <gsModeling/gsPlanarDomainGenerator.h>
#include <gsModeling/gsPatchSplitterGeom.hpp>
#include <gsModeling/gsPatchSplitterTriMesh.hpp>
#include <gsModeling/gsASG1LocalFitter.h>
#include <gsModeling/gsCoonsPatch.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsReadFile.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseSolver.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>

namespace gismo {

// ============================================================================
// Constructor & Ingestion
// ============================================================================

template <class T>
gsOptionList gsPatchSplitter<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("Strategy", "0: Auto/Hybrid, 1: CrossField, 2: ConcaveCut/MAT, 3: PolychordOnly", 0);
    opt.addReal("TargetAspectRatio", "Target patch aspect ratio Lu/Lv", 1.0);
    opt.addReal("MaxAspectRatio", "Preferred maximum patch aspect ratio", 2.0);
    opt.addReal("PatchCountWeight", "Penalty weight per patch in polychord optimization", 0.18);
    opt.addReal("AngleThresholdDeg", "Boundary corner angle above which corner is split", 135.0);
    opt.addReal("MinCornerAngleDeg", "Minimum admissible patch corner tangent angle (deg)", 12.0);
    opt.addReal("MaxCornerAngleDeg", "Maximum admissible patch corner tangent angle (deg)", 155.0);
    opt.addReal("StripCollapseRatio", "Relative distance to snap parallel separatrices", 0.25);
    opt.addInt("Degree", "Target B-spline degree (>= 3 for ArgyrisSplines)", 3);
    opt.addInt("Regularity", "Target spline regularity r (1 <= r <= p - 2)", 1);
    opt.addInt("MinKnotSpans", "Minimum knot spans per patch direction (>= 3 for ArgyrisSplines)", 3);
    opt.addSwitch("EnablePolychordSplitting", "Run global polychord aspect-ratio optimization in compute()", true);
    opt.addSwitch("ApplyASG1Fitting", "Run gsASG1LocalFitter in formatForArgyris()", false);
    opt.addInt("ParameterizationMethod", "0: Coons + Spring + Winslow, 1: Coons + gsBarrierPatch (if enabled)", 0);
    return opt;
}

template <class T>
gsPatchSplitter<T>::gsPatchSplitter()
    : m_options(defaultOptions())
{
}

template <class T>
gsPatchSplitter<T>::gsPatchSplitter(const gsPlanarDomain<T>& domain)
    : m_domain(domain),
      m_hasDomain(true),
      m_options(defaultOptions())
{
}

template <class T>
gsPatchSplitter<T>::gsPatchSplitter(const gsMultiPatch<T>& inputMp)
    : m_inputMp(inputMp),
      m_hasInputMp(true),
      m_options(defaultOptions())
{
}

template <class T>
void gsPatchSplitter<T>::setDomain(const gsPlanarDomain<T>& domain)
{
    m_domain = domain;
    m_hasDomain = true;
    m_hasInputMp = false;
    m_cachedDomainMinCornerAngle = -1;
    m_bestCandidate = SplittingCandidate();
}

template <class T>
void gsPatchSplitter<T>::setMultiPatch(const gsMultiPatch<T>& inputMp)
{
    m_inputMp = inputMp;
    m_hasInputMp = true;
    m_bestCandidate = SplittingCandidate();
}

template <class T>
void gsPatchSplitter<T>::loadDomainFromXml(const std::string& xml_file)
{
    gsReadFile<T> rf(xml_file);
    memory::unique_ptr<gsPlanarDomain<T>> pDom = rf;
    if (pDom) {
        m_domain = *pDom;
        m_hasDomain = true;
        m_hasInputMp = false;
        m_cachedDomainMinCornerAngle = -1;
        m_bestCandidate = SplittingCandidate();
    } else {
        GISMO_ERROR("XML file " << xml_file << " does not contain a gsPlanarDomain.");
    }
}

template <class T>
void gsPatchSplitter<T>::setStrategy(Strategy strategy)
{
    m_strategy = strategy;
    m_options.setInt("Strategy", static_cast<int>(strategy));
}

template <class T>
void gsPatchSplitter<T>::setDegree(int degree)
{
    m_degree = std::max(1, degree);
    m_options.setInt("Degree", m_degree);
}

template <class T>
void gsPatchSplitter<T>::setWeights(T weightPatches, T weightDistortion)
{
    m_weightPatches = weightPatches;
    m_weightDistortion = weightDistortion;
}

// ============================================================================
// Helper: Curve Handling & Exact CAD B-Spline Splitting
// ============================================================================

namespace internal_splitter {
// Shared helper returning the exterior CCW tangent cone [aStart, aStart + dCCW] at the junction prevSeg -> nextSeg.
template <class T>
inline bool computeExteriorCone(
    const CurveSeg<T>& prevSeg,
    const CurveSeg<T>& nextSeg,
    const gsPoint2D<T>& dir,
    T& aStart, T& dCCW, T& aDir, T& minEdgeLen)
{
    auto dIn  = evalCurveDeriv(*prevSeg.bsp, T(1));
    auto dOut = evalCurveDeriv(*nextSeg.bsp, T(0));
    T nIn = std::hypot(dIn.x, dIn.y), nOut = std::hypot(dOut.x, dOut.y);
    T lenD = std::hypot(dir.x, dir.y);
    if (nIn <= T(1e-14) || nOut <= T(1e-14) || lenD <= T(1e-14))
        return false;
    const T twoPi = T(2) * T(3.14159265358979323846);
    gsPoint2D<T> tBack(-dIn.x / nIn, -dIn.y / nIn);
    gsPoint2D<T> tOut(dOut.x / nOut, dOut.y / nOut);
    aStart = std::atan2(tBack.y, tBack.x);
    dCCW = std::atan2(tOut.y, tOut.x) - aStart;
    while (dCCW <= T(0)) dCCW += twoPi;
    while (dCCW > twoPi) dCCW -= twoPi;
    aDir = std::atan2(dir.y / lenD, dir.x / lenD) - aStart;
    while (aDir < T(0)) aDir += twoPi;
    while (aDir >= twoPi) aDir -= twoPi;
    T lPrev = std::hypot(prevSeg.p1().x - prevSeg.p0().x, prevSeg.p1().y - prevSeg.p0().y);
    T lNext = std::hypot(nextSeg.p1().x - nextSeg.p0().x, nextSeg.p1().y - nextSeg.p0().y);
    minEdgeLen = std::min(lPrev, lNext);
    return true;
}

// Build a spoke CurveSeg H -> P; curves gently into the exterior safe cone only when a
// straight spoke grazes a hole tangent (< 14 deg) or creates a > 172 deg angle at a sharp corner.
template <class T>
inline CurveSeg<T> makeExteriorConeSafeSpoke(
    const CurveSeg<T>& prevSeg,
    const CurveSeg<T>& nextSeg,
    const gsPoint2D<T>& H,
    const gsPoint2D<T>& P,
    const std::vector<gsPoint2D<T>>& holePoly,
    int degree,
    bool& wasCurved)
{
    wasCurved = false;
    auto mkLin = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        std::unique_ptr<gsBSpline<T>> ptr(makeStraightBSpline(a, b, degree));
        return CurveSeg<T>(*ptr);
    };
    gsPoint2D<T> diff(P.x - H.x, P.y - H.y);
    T lenD = std::hypot(diff.x, diff.y);
    T aStart = 0, dCCW = 0, aDir = 0, minEdgeLen = 0;
    if (lenD > T(1e-14) && computeExteriorCone(prevSeg, nextSeg, diff, aStart, dCCW, aDir, minEdgeLen)) {
        const T pi = T(3.14159265358979323846);
        const T trigMargin = T(0.32);
        if (aDir < trigMargin || aDir > dCCW - trigMargin || aDir > pi - trigMargin || (dCCW - aDir) > pi - trigMargin) {
            const T targetMargin = T(0.36);
            T loSafe = std::max(targetMargin, dCCW - pi + targetMargin);
            T hiSafe = std::min(pi - targetMargin, dCCW - targetMargin);
            T targetRel = (loSafe < hiSafe) ? std::max(loSafe, std::min(hiSafe, aDir)) : (T(0.5) * dCCW);
            T absAng = aStart + targetRel;
            gsPoint2D<T> tanDir(std::cos(absAng), std::sin(absAng));
            // Smooth convex (zero-inflection) cubic Bezier: compact c1 along safe tangent, c2 = midpoint(c1, P).
            T tLen = std::min(T(0.18) * lenD, std::max(T(0.06) * lenD, T(0.85) * minEdgeLen));
            gsPoint2D<T> c1(H.x + tLen * tanDir.x, H.y + tLen * tanDir.y);
            gsPoint2D<T> c2(T(0.5) * (c1.x + P.x), T(0.5) * (c1.y + P.y));
            gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
            gsMatrix<T> cf(4, 2);
            cf << H.x, H.y, c1.x, c1.y, c2.x, c2.y, P.x, P.y;
            gsBSpline<T> bsp(kv, cf);
            if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
            bool curveClear = true;
            for (int k = 1; k < 16; ++k) {
                gsPoint2D<T> q = evalCurvePoint(bsp, static_cast<T>(k) / T(16));
                if (pointInPolyStrict(q, holePoly, T(1e-7))) { curveClear = false; break; }
            }
            if (curveClear) {
                wasCurved = true;
                return CurveSeg<T>(bsp);
            }
        }
    }
    return mkLin(H, P);
}

// Adjust an interior cap point relative to hole vertex vm so the segment vm -> pt
// stays inside the safe exterior convex cone (both quad corner angles at vm < 166 deg).
template <class T>
inline gsPoint2D<T> clampPointToExteriorCone(
    const CurveSeg<T>& prevSeg,
    const CurveSeg<T>& nextSeg,
    const gsPoint2D<T>& vm,
    const gsPoint2D<T>& cand)
{
    gsPoint2D<T> diff(cand.x - vm.x, cand.y - vm.y);
    T lenD = std::hypot(diff.x, diff.y);
    T aStart = 0, dCCW = 0, aDir = 0, minEdgeLen = 0;
    if (lenD <= T(1e-14) || !computeExteriorCone(prevSeg, nextSeg, diff, aStart, dCCW, aDir, minEdgeLen))
        return cand;
    const T pi = T(3.14159265358979323846);
    T loSafe = std::max(T(0.32), dCCW - pi + T(0.32));
    T hiSafe = std::min(pi - T(0.32), dCCW - T(0.32));
    T targetRel = (loSafe < hiSafe) ? std::max(loSafe, std::min(hiSafe, aDir)) : (T(0.5) * dCCW);
    if (std::abs(targetRel - aDir) < T(1e-6))
        return cand;
    T absAng = aStart + targetRel;
    return gsPoint2D<T>(vm.x + lenD * std::cos(absAng), vm.y + lenD * std::sin(absAng));
}

// Build a sharp-corner-preserving conforming quad ring inside a 4-sided block [cS, cE, cN, cW]
// for an arbitrary even-segment CCW hole (or K=3 split to 4).
// All hole corners become patch vertices; outer block edges cS, cE, cN, cW remain unsplit.
template <class T>
inline bool buildSharpCornerBlockRing(
    const CurveSeg<T>& cS, const CurveSeg<T>& cE,
    const CurveSeg<T>& cN, const CurveSeg<T>& cW,
    const std::vector<CurveSeg<T>>& hCCW,
    const int anchorIdx[4],
    const std::vector<gsPoint2D<T>>& holePoly,
    std::vector<std::unique_ptr<gsGeometry<T>>>& outPatches,
    int degree,
    T outWeight = T(0.62),
    std::string* failReason = nullptr,
    T spokeSplitT = T(0.5))
{
    const int K = static_cast<int>(hCCW.size());
    // Compute sectors: sector s runs from anchorIdx[s] to anchorIdx[(s+1)%4]
    std::vector<std::vector<CurveSeg<T>>> sectors(4);
    for (int s = 0; s < 4; ++s) {
        int i0 = anchorIdx[s];
        int i1 = anchorIdx[(s + 1) % 4];
        int cur = i0;
        for (int step = 0; step < K; ++step) {
            sectors[s].push_back(hCCW[cur]);
            cur = (cur + 1) % K;
            if (cur == i1) break;
        }
        if (sectors[s].empty()) return false;
    }

    const gsPoint2D<T> P[4] = { cS.p0(), cS.p1(), cE.p1(), cN.p0() }; // SW, SE, NE, NW
    const gsPoint2D<T> H[4] = {
        hCCW[anchorIdx[0]].p0(), hCCW[anchorIdx[1]].p0(),
        hCCW[anchorIdx[2]].p0(), hCCW[anchorIdx[3]].p0()
    };

    // Quick visibility check on anchor spokes
    for (int s = 0; s < 4; ++s) {
        if (!isHoleSpokeClearOfHole(H[s], P[s], holePoly)) {
            if (failReason) *failReason = "spokeBlocked_s" + std::to_string(s);
            return false;
        }
    }

    int totalRem = 0;
    for (int s = 0; s < 4; ++s) totalRem += static_cast<int>(sectors[s].size()) - 1;
    if (totalRem % 2 != 0) return false;

    // Helper: create a straight spoke CurveSeg A -> B
    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    struct QuadGraph {
        std::vector<gsPoint2D<T>> pts;
        std::vector<bool> pinned;
        struct QItem {
            int vid[4];
            CurveSeg<T> e[4]; // cycle: v0-(e0)->v1-(e1)->v2-(e2)->v3-(e3)->v0
            bool straight[4];
        };
        std::vector<QItem> items;

        int addVertex(const gsPoint2D<T>& p, bool isPinned) {
            for (size_t i = 0; i < pts.size(); ++i) {
                T dx = pts[i].x - p.x, dy = pts[i].y - p.y;
                if (dx * dx + dy * dy < T(1e-20)) return static_cast<int>(i);
            }
            pts.push_back(p);
            pinned.push_back(isPinned);
            return static_cast<int>(pts.size()) - 1;
        }
        int addFreeVertex(const gsPoint2D<T>& p) {
            pts.push_back(p);
            pinned.push_back(false);
            return static_cast<int>(pts.size()) - 1;
        }
        void addQuad(const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                     const CurveSeg<T>& e23, const CurveSeg<T>& e30,
                     int v0 = -1, int v1 = -1, int v2 = -1, int v3 = -1) {
            QItem q;
            q.vid[0] = (v0 >= 0) ? v0 : addVertex(e01.p0(), true);
            q.vid[1] = (v1 >= 0) ? v1 : addVertex(e12.p0(), true);
            q.vid[2] = (v2 >= 0) ? v2 : addVertex(e23.p0(), true);
            q.vid[3] = (v3 >= 0) ? v3 : addVertex(e30.p0(), true);
            q.e[0] = e01; q.e[1] = e12; q.e[2] = e23; q.e[3] = e30;
            for (int m = 0; m < 4; ++m)
                q.straight[m] = (!pinned[q.vid[m]] || !pinned[q.vid[(m + 1) % 4]]);
            items.push_back(q);
        }
    } qg;

    // Quad helper storing graph record
    auto addCCWQuad = [&](const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                          const CurveSeg<T>& e23, const CurveSeg<T>& e30,
                          int v0 = -1, int v1 = -1, int v2 = -1, int v3 = -1) {
        qg.addQuad(e01, e12, e23, e30, v0, v1, v2, v3);
    };

    auto projU = [](const gsPoint2D<T>& q, const gsPoint2D<T>& Tl, const gsPoint2D<T>& Tr, T uMin, T uMax) -> T {
        T dx = Tr.x - Tl.x, dy = Tr.y - Tl.y;
        T len2 = dx * dx + dy * dy;
        if (len2 <= T(1e-20)) return T(0.5) * (uMin + uMax);
        T u = ((q.x - Tl.x) * dx + (q.y - Tl.y) * dy) / len2;
        return std::max(uMin, std::min(uMax, u));
    };
    auto lerpPt = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b, T u) -> gsPoint2D<T> {
        return gsPoint2D<T>((T(1) - u) * a.x + u * b.x, (T(1) - u) * a.y + u * b.y);
    };
    CurveSeg<T> outerTop[4] = { cS, cE, cN.reversed(), cW.reversed() };
    std::vector<gsPoint2D<T>> blockPoly;
    blockPoly.reserve(32);
    for (int s = 0; s < 4; ++s) {
        for (int k = 0; k < 8; ++k)
            blockPoly.push_back(evalCurvePoint(*outerTop[s].bsp, static_cast<T>(k) / T(8)));
    }
    CurveSeg<T> sp[4];
    bool spCurved[4] = {false, false, false, false};
    int anchorVid[4] = {-1, -1, -1, -1};
    for (int s = 0; s < 4; ++s) {
        const int idx = anchorIdx[s];
        const auto& prevSeg = hCCW[(idx + K - 1) % K];
        const auto& nextSeg = hCCW[idx];
        sp[s] = makeExteriorConeSafeSpoke(prevSeg, nextSeg, H[s], P[s], holePoly, degree, spCurved[s]);
        anchorVid[s] = qg.addVertex(H[s], true);
    }

    struct SectorFrame {
        std::vector<CurveSeg<T>> holeEdges; // size 1 or 2 after Step A
        CurveSeg<T> topEdge;                // T_left -> T_right
        CurveSeg<T> leftSpoke;              // B_left -> T_left (outward)
        CurveSeg<T> rightSpoke;             // B_right -> T_right (outward)
        int vidTl = -1, vidTr = -1;
    };

    auto placeCapMidpoint = [&](const SectorFrame& f, gsPoint2D<T>& Qout) -> bool {
        const gsPoint2D<T> Tl = qg.pts[f.vidTl], Tr = qg.pts[f.vidTr];
        const CurveSeg<T>& segIn  = f.holeEdges[0];
        const CurveSeg<T>& segOut = f.holeEdges[1];
        gsPoint2D<T> vm = segIn.p1();
        T u = projU(vm, Tl, Tr, T(0.20), T(0.80));
        gsPoint2D<T> T_out = evalCurvePoint(*f.topEdge.bsp, u);
        Qout = gsPoint2D<T>(outWeight * T_out.x + (T(1) - outWeight) * vm.x,
                            outWeight * T_out.y + (T(1) - outWeight) * vm.y);
        Qout = clampPointToExteriorCone(segIn, segOut, vm, Qout);
        if (pointInPolyStrict(Qout, holePoly, T(1e-7)))
            return false;
        if (!isHoleSpokeClearOfHole(vm, Qout, holePoly))
            return false;
        return true;
    };

    SectorFrame fr[4];
    for (int s = 0; s < 4; ++s) {
        fr[s].holeEdges = sectors[s];
        fr[s].topEdge = outerTop[s];
        fr[s].leftSpoke = sp[s];
        fr[s].rightSpoke = sp[(s + 1) % 4];
        fr[s].vidTl = qg.addVertex(fr[s].topEdge.p0(), true);
        fr[s].vidTr = qg.addVertex(fr[s].topEdge.p1(), true);
    }

    // Track sectors whose size==2 needs half-sector closure at a pre-split spoke midpoint M.
    struct PendingTwo {
        bool active = false;
        bool mOnRight = false;
        int vM = -1;
        CurveSeg<T> spLo, spHi;
    } pend2[4];

    auto peelSize4RightM = [&](int s, int vM, const gsPoint2D<T>& M,
                               const CurveSeg<T>& spLo, const CurveSeg<T>& spHi) -> bool {
        gsPoint2D<T> Tl = qg.pts[fr[s].vidTl], Tr = qg.pts[fr[s].vidTr];
        gsPoint2D<T> v0 = fr[s].holeEdges[0].p0(), v1 = fr[s].holeEdges[0].p1();
        T u1 = projU(v1, Tl, Tr, T(0.18), T(0.55));
        gsPoint2D<T> T1 = evalCurvePoint(*fr[s].topEdge.bsp, u1);
        gsPoint2D<T> p1(outWeight * T1.x + (T(1) - outWeight) * v1.x,
                        outWeight * T1.y + (T(1) - outWeight) * v1.y);
        if (fr[s].holeEdges.size() >= 2) {
            gsPoint2D<T> p1c = clampPointToExteriorCone(fr[s].holeEdges[0], fr[s].holeEdges[1], v1, p1);
            if (!pointInPolyStrict(p1c, holePoly, T(1e-7)) &&
                pointInPolyStrict(p1c, blockPoly, T(1e-7)) &&
                isHoleSpokeClearOfHole(v1, p1c, holePoly))
                p1 = p1c;
        }
        if (pointInPolyStrict(p1, holePoly, T(1e-7)) ||
            !pointInPolyStrict(p1, blockPoly, T(1e-7)) ||
            !isHoleSpokeClearOfHole(v1, p1, holePoly) ||
            !isHoleSpokeClearOfHole(p1, M, holePoly))
            return false;
        int vp1 = qg.addFreeVertex(p1);
        int vv0 = qg.addVertex(v0, true), vv1 = qg.addVertex(v1, true);
        CurveSeg<T> e_p1M = mkLine(p1, M), e_v1p1 = mkLine(v1, p1), e_p1Tl = mkLine(p1, Tl);
        addCCWQuad(e_p1M, spHi, fr[s].topEdge.reversed(), e_p1Tl.reversed(), vp1, vM, fr[s].vidTr, fr[s].vidTl);
        addCCWQuad(fr[s].holeEdges[0], e_v1p1, e_p1Tl, fr[s].leftSpoke.reversed(), vv0, vv1, vp1, fr[s].vidTl);
        fr[s].holeEdges.erase(fr[s].holeEdges.begin());
        fr[s].topEdge = e_p1M;
        fr[s].leftSpoke = e_v1p1;
        fr[s].rightSpoke = spLo;
        fr[s].vidTl = vp1;
        fr[s].vidTr = vM;
        return true;
    };

    auto peelSize4LeftM = [&](int s, int vM, const gsPoint2D<T>& M,
                              const CurveSeg<T>& spLo, const CurveSeg<T>& spHi) -> bool {
        gsPoint2D<T> Tl = qg.pts[fr[s].vidTl], Tr = qg.pts[fr[s].vidTr];
        gsPoint2D<T> v3 = fr[s].holeEdges.back().p0(), v4 = fr[s].holeEdges.back().p1();
        T u2 = projU(v3, Tl, Tr, T(0.45), T(0.82));
        gsPoint2D<T> T2 = evalCurvePoint(*fr[s].topEdge.bsp, u2);
        gsPoint2D<T> p2(outWeight * T2.x + (T(1) - outWeight) * v3.x,
                        outWeight * T2.y + (T(1) - outWeight) * v3.y);
        if (fr[s].holeEdges.size() >= 2) {
            const size_t nh = fr[s].holeEdges.size();
            gsPoint2D<T> p2c = clampPointToExteriorCone(fr[s].holeEdges[nh - 2], fr[s].holeEdges.back(), v3, p2);
            if (!pointInPolyStrict(p2c, holePoly, T(1e-7)) &&
                pointInPolyStrict(p2c, blockPoly, T(1e-7)) &&
                isHoleSpokeClearOfHole(v3, p2c, holePoly))
                p2 = p2c;
        }
        if (pointInPolyStrict(p2, holePoly, T(1e-7)) ||
            !pointInPolyStrict(p2, blockPoly, T(1e-7)) ||
            !isHoleSpokeClearOfHole(v3, p2, holePoly) ||
            !isHoleSpokeClearOfHole(M, p2, holePoly))
            return false;
        int vp2 = qg.addFreeVertex(p2);
        int vv3 = qg.addVertex(v3, true), vv4 = qg.addVertex(v4, true);
        CurveSeg<T> e_Mp2 = mkLine(M, p2), e_v3p2 = mkLine(v3, p2), e_p2Tr = mkLine(p2, Tr);
        addCCWQuad(e_Mp2, e_p2Tr, fr[s].topEdge.reversed(), spHi.reversed(), vM, vp2, fr[s].vidTr, fr[s].vidTl);
        addCCWQuad(fr[s].holeEdges.back(), fr[s].rightSpoke, e_p2Tr.reversed(), e_v3p2.reversed(), vv3, vv4, fr[s].vidTr, vp2);
        fr[s].holeEdges.pop_back();
        fr[s].topEdge = e_Mp2;
        fr[s].leftSpoke = spLo;
        fr[s].rightSpoke = e_v3p2;
        fr[s].vidTl = vM;
        fr[s].vidTr = vp2;
        return true;
    };

    // Pre-pair all even sectors (size 2 or 4) BEFORE Step A:
    // - Adjacent pair (a, b=(a+1)%4): split shared spoke sp[b] at M
    // - Opposite pair (a, b=(a+2)%4): split bridge oppMid=(a+1)%4 spokes at M1, M2 and emit bridge top strip
    // - Any size==4 sector in a pair peels 1 edge into its split M, becoming size==3 for Step A!
    {
        const T tSplit = std::max(T(0.15), std::min(T(0.85), spokeSplitT));
        std::vector<int> evenSec;
        for (int s = 0; s < 4; ++s)
            if (fr[s].holeEdges.size() % 2 == 0) evenSec.push_back(s);

        auto bindRightLeft = [&](int u, int v, int vM, const gsPoint2D<T>& M,
                                 const CurveSeg<T>& spLo, const CurveSeg<T>& spHi) -> bool {
            if (fr[u].holeEdges.size() == 4) {
                if (!peelSize4RightM(u, vM, M, spLo, spHi)) return false;
            } else {
                pend2[u] = {true, true, vM, spLo, spHi};
            }
            if (fr[v].holeEdges.size() == 4) {
                if (!peelSize4LeftM(v, vM, M, spLo, spHi)) return false;
            } else {
                pend2[v] = {true, false, vM, spLo, spHi};
            }
            return true;
        };

        auto pairAdjacent = [&](int a, int b) -> bool {
            gsBSpline<T> p0, p1;
            splitBSplineExact(*fr[a].rightSpoke.bsp, tSplit, p0, p1);
            CurveSeg<T> spLo(p0), spHi(p1);
            gsPoint2D<T> M = spLo.p1();
            int spIdx = (a + 1) % 4;
            int vM = spCurved[spIdx] ? qg.addVertex(M, true) : qg.addFreeVertex(M);
            return bindRightLeft(a, b, vM, M, spLo, spHi);
        };

        if (evenSec.size() == 2) {
            int a = evenSec[0], b = evenSec[1];
            if ((a + 1) % 4 == b) {
                if (!pairAdjacent(a, b)) return false;
            } else if ((b + 1) % 4 == a) {
                if (!pairAdjacent(b, a)) return false;
            } else {
                int mid = (a + 1) % 4;
                gsBSpline<T> am0, am1, mb0, mb1;
                splitBSplineExact(*fr[mid].leftSpoke.bsp, tSplit, am0, am1);
                splitBSplineExact(*fr[mid].rightSpoke.bsp, tSplit, mb0, mb1);
                CurveSeg<T> spAM_lo(am0), spAM_hi(am1), spMB_lo(mb0), spMB_hi(mb1);
                gsPoint2D<T> M1 = spAM_lo.p1(), M2 = spMB_lo.p1();
                int vM1 = spCurved[mid] ? qg.addVertex(M1, true) : qg.addFreeVertex(M1);
                int vM2 = spCurved[(mid + 1) % 4] ? qg.addVertex(M2, true) : qg.addFreeVertex(M2);
                CurveSeg<T> e_M1M2 = mkLine(M1, M2);
                addCCWQuad(e_M1M2, spMB_hi, fr[mid].topEdge.reversed(), spAM_hi.reversed(),
                           vM1, vM2, fr[mid].vidTr, fr[mid].vidTl);
                fr[mid].topEdge = e_M1M2;
                fr[mid].leftSpoke = spAM_lo;
                fr[mid].rightSpoke = spMB_lo;
                fr[mid].vidTl = vM1;
                fr[mid].vidTr = vM2;
                if (fr[a].holeEdges.size() == 4) {
                    if (!peelSize4RightM(a, vM1, M1, spAM_lo, spAM_hi)) return false;
                } else {
                    pend2[a] = {true, true, vM1, spAM_lo, spAM_hi};
                }
                if (fr[b].holeEdges.size() == 4) {
                    if (!peelSize4LeftM(b, vM2, M2, spMB_lo, spMB_hi)) return false;
                } else {
                    pend2[b] = {true, false, vM2, spMB_lo, spMB_hi};
                }
            }
        } else if (evenSec.size() == 4) {
            if (!pairAdjacent(0, 1) || !pairAdjacent(2, 3))
                return false;
        }
    }

    // Step A: Recursively apply 3->1 pillow cap in any sector with >= 3 hole edges until every sector has 1 or 2 edges.
    for (int s = 0; s < 4; ++s) {
        while (fr[s].holeEdges.size() >= 3) {
            gsPoint2D<T> Tl = qg.pts[fr[s].vidTl];
            gsPoint2D<T> Tr = qg.pts[fr[s].vidTr];
            gsPoint2D<T> v1 = fr[s].holeEdges.front().p1();
            gsPoint2D<T> v2 = fr[s].holeEdges.back().p0();

            // Place p1, p2 on a chord parallel to Tl--Tr strictly in front of both v1 and v2;
            // fall back to curve-point interpolation if topEdge curves past the straight chord (e.g., bent starburst rays).
            T dxT = Tr.x - Tl.x, dyT = Tr.y - Tl.y;
            T lenT = std::hypot(dxT, dyT);
            if (lenT <= T(1e-14)) return false;
            gsPoint2D<T> nIn(-dyT / lenT, dxT / lenT); // CCW interior normal (domain on left of Tl->Tr)
            T h1 = (v1.x - Tl.x) * nIn.x + (v1.y - Tl.y) * nIn.y;
            T h2 = (v2.x - Tl.x) * nIn.x + (v2.y - Tl.y) * nIn.y;
            T hMin = std::min(h1, h2);
            T u1 = projU(v1, Tl, Tr, T(0.15), T(0.68));
            T u2 = projU(v2, Tl, Tr, std::min(T(0.85), u1 + T(0.12)), T(0.85));
            gsPoint2D<T> p1, p2;
            if (hMin > T(1e-6)) {
                T hCap = (T(1) - outWeight) * hMin;
                gsPoint2D<T> T_left13 = lerpPt(Tl, Tr, u1);
                gsPoint2D<T> T_rght23 = lerpPt(Tl, Tr, u2);
                p1 = gsPoint2D<T>(T_left13.x + hCap * nIn.x, T_left13.y + hCap * nIn.y);
                p2 = gsPoint2D<T>(T_rght23.x + hCap * nIn.x, T_rght23.y + hCap * nIn.y);
            } else {
                gsPoint2D<T> T1 = evalCurvePoint(*fr[s].topEdge.bsp, u1);
                gsPoint2D<T> T2 = evalCurvePoint(*fr[s].topEdge.bsp, u2);
                T wC = T(0.65) * outWeight;
                p1 = gsPoint2D<T>((T(1) - wC) * v1.x + wC * T1.x, (T(1) - wC) * v1.y + wC * T1.y);
                p2 = gsPoint2D<T>((T(1) - wC) * v2.x + wC * T2.x, (T(1) - wC) * v2.y + wC * T2.y);
            }
            {
                const size_t nh = fr[s].holeEdges.size();
                gsPoint2D<T> p1c = clampPointToExteriorCone(fr[s].holeEdges[0], fr[s].holeEdges[1], v1, p1);
                gsPoint2D<T> p2c = clampPointToExteriorCone(fr[s].holeEdges[nh - 2], fr[s].holeEdges.back(), v2, p2);
                if (!pointInPolyStrict(p1c, holePoly, T(1e-7)) &&
                    pointInPolyStrict(p1c, blockPoly, T(1e-7)) &&
                    isHoleSpokeClearOfHole(v1, p1c, holePoly))
                    p1 = p1c;
                if (!pointInPolyStrict(p2c, holePoly, T(1e-7)) &&
                    pointInPolyStrict(p2c, blockPoly, T(1e-7)) &&
                    isHoleSpokeClearOfHole(v2, p2c, holePoly))
                    p2 = p2c;
            }

            if (pointInPolyStrict(p1, holePoly, T(1e-7)) || pointInPolyStrict(p2, holePoly, T(1e-7)) ||
                !pointInPolyStrict(p1, blockPoly, T(1e-7)) || !pointInPolyStrict(p2, blockPoly, T(1e-7)) ||
                !isHoleSpokeClearOfHole(v1, p1, holePoly) || !isHoleSpokeClearOfHole(v2, p2, holePoly) ||
                !isHoleSpokeClearOfHole(p1, p2, holePoly)) {
                if (failReason && hMin <= T(1e-6)) *failReason = "hMin<=0";
                return false;
            }

            int vp1 = qg.addFreeVertex(p1);
            int vp2 = qg.addFreeVertex(p2);
            int vBl = qg.addVertex(fr[s].holeEdges.front().p0(), true);
            int vv1 = qg.addVertex(v1, true);
            int vv2 = qg.addVertex(v2, true);
            int vBr = qg.addVertex(fr[s].holeEdges.back().p1(), true);

            CurveSeg<T> top_p1p2 = mkLine(p1, p2);
            CurveSeg<T> left_v1p1 = mkLine(v1, p1);
            CurveSeg<T> right_v2p2 = mkLine(v2, p2);
            CurveSeg<T> capL = mkLine(p1, Tl);
            CurveSeg<T> capR = mkLine(p2, Tr);

            addCCWQuad(top_p1p2, capR, fr[s].topEdge.reversed(), capL.reversed(), vp1, vp2, fr[s].vidTr, fr[s].vidTl);
            addCCWQuad(fr[s].holeEdges.front(), left_v1p1, capL, fr[s].leftSpoke.reversed(), vBl, vv1, vp1, fr[s].vidTl);
            addCCWQuad(fr[s].holeEdges.back(), fr[s].rightSpoke, capR.reversed(), right_v2p2.reversed(), vv2, vBr, fr[s].vidTr, vp2);

            std::vector<CurveSeg<T>> mid(fr[s].holeEdges.begin() + 1, fr[s].holeEdges.end() - 1);
            fr[s].holeEdges = std::move(mid);
            fr[s].topEdge = top_p1p2;
            fr[s].leftSpoke = left_v1p1;
            fr[s].rightSpoke = right_v2p2;
            fr[s].vidTl = vp1;
            fr[s].vidTr = vp2;
        }
    }

    // Emit remaining sectors cleanly (pend2 sectors close into pre-split M, unhandled sectors have size 1).
    for (int s = 0; s < 4; ++s) {
        if (pend2[s].active) {
            if (fr[s].holeEdges.size() != 2) return false;
            gsPoint2D<T> Q;
            if (!placeCapMidpoint(fr[s], Q)) return false;
            int vQ = qg.addFreeVertex(Q);
            int vBl = qg.addVertex(fr[s].holeEdges[0].p0(), true);
            gsPoint2D<T> vm = fr[s].holeEdges[0].p1();
            int vvm = qg.addVertex(vm, true);
            int vBr = qg.addVertex(fr[s].holeEdges[1].p1(), true);
            int vM = pend2[s].vM;
            gsPoint2D<T> M = qg.pts[vM];
            CurveSeg<T> e_vmQ = mkLine(vm, Q);
            if (pend2[s].mOnRight) {
                CurveSeg<T> e_QTl = mkLine(Q, qg.pts[fr[s].vidTl]), e_QM = mkLine(Q, M);
                addCCWQuad(fr[s].holeEdges[0], e_vmQ, e_QTl, fr[s].leftSpoke.reversed(), vBl, vvm, vQ, fr[s].vidTl);
                addCCWQuad(fr[s].holeEdges[1], pend2[s].spLo, e_QM.reversed(), e_vmQ.reversed(), vvm, vBr, vM, vQ);
                addCCWQuad(e_QM, pend2[s].spHi, fr[s].topEdge.reversed(), e_QTl.reversed(), vQ, vM, fr[s].vidTr, fr[s].vidTl);
            } else {
                CurveSeg<T> e_QTr = mkLine(Q, qg.pts[fr[s].vidTr]), e_MQ = mkLine(M, Q);
                addCCWQuad(fr[s].holeEdges[0], e_vmQ, e_MQ.reversed(), pend2[s].spLo.reversed(), vBl, vvm, vQ, vM);
                addCCWQuad(fr[s].holeEdges[1], fr[s].rightSpoke, e_QTr.reversed(), e_vmQ.reversed(), vvm, vBr, fr[s].vidTr, vQ);
                addCCWQuad(e_MQ, e_QTr, fr[s].topEdge.reversed(), pend2[s].spHi.reversed(), vM, vQ, fr[s].vidTr, fr[s].vidTl);
            }
        } else {
            if (fr[s].holeEdges.size() != 1) return false;
            int vBl = qg.addVertex(fr[s].holeEdges[0].p0(), true);
            int vBr = qg.addVertex(fr[s].holeEdges[0].p1(), true);
            addCCWQuad(fr[s].holeEdges[0], fr[s].rightSpoke, fr[s].topEdge.reversed(), fr[s].leftSpoke.reversed(),
                       vBl, vBr, fr[s].vidTr, fr[s].vidTl);
        }
    }

    // Guarded Laplacian + convexity hill-climb on free interior vertices before Coons generation
    {
        const int NV = static_cast<int>(qg.pts.size());
        std::vector<bool> hasCurvedAnchor(NV, false);
        for (int s = 0; s < 4; ++s) {
            if (spCurved[s] && anchorVid[s] >= 0 && anchorVid[s] < NV)
                hasCurvedAnchor[anchorVid[s]] = true;
        }
        std::vector<std::vector<int>> adj(NV);
        std::vector<std::vector<int>> incQ(NV);
        for (size_t qi_idx = 0; qi_idx < qg.items.size(); ++qi_idx) {
            const auto& qi = qg.items[qi_idx];
            for (int m = 0; m < 4; ++m) {
                adj[qi.vid[m]].push_back(qi.vid[(m + 1) % 4]);
                adj[qi.vid[(m + 1) % 4]].push_back(qi.vid[m]);
                incQ[qi.vid[m]].push_back(static_cast<int>(qi_idx));
            }
        }
        auto minIncSin = [&](int v) -> T {
            T minS = T(1e18);
            for (int qi_idx : incQ[v]) {
                const auto& qi = qg.items[qi_idx];
                gsPoint2D<T> cp[4] = { qg.pts[qi.vid[3]], qg.pts[qi.vid[2]], qg.pts[qi.vid[1]], qg.pts[qi.vid[0]] };
                int vIdx[4] = { qi.vid[3], qi.vid[2], qi.vid[1], qi.vid[0] };
                for (int k = 0; k < 4; ++k) {
                    if (hasCurvedAnchor[vIdx[(k + 1) % 4]]) continue;
                    const auto& A = cp[k];
                    const auto& B = cp[(k + 1) % 4];
                    const auto& C = cp[(k + 2) % 4];
                    T cr = cross2(A, B, C);
                    T l1 = std::hypot(B.x - A.x, B.y - A.y);
                    T l2 = std::hypot(C.x - B.x, C.y - B.y);
                    T sA = (l1 > T(1e-12) && l2 > T(1e-12)) ? (cr / (l1 * l2)) : T(-1);
                    minS = std::min(minS, sA);
                }
            }
            return minS;
        };
        auto allQuadsConvex = [&](std::string* qErr = nullptr) -> bool {
            for (size_t idx = 0; idx < qg.items.size(); ++idx) {
                const auto& qi = qg.items[idx];
                gsPoint2D<T> cp[4] = { qg.pts[qi.vid[3]], qg.pts[qi.vid[2]], qg.pts[qi.vid[1]], qg.pts[qi.vid[0]] };
                int vIdx[4] = { qi.vid[3], qi.vid[2], qi.vid[1], qi.vid[0] };
                for (int k = 0; k < 4; ++k) {
                    const auto& A = cp[k];
                    const auto& B = cp[(k + 1) % 4];
                    const auto& C = cp[(k + 2) % 4];
                    T cr = cross2(A, B, C);
                    T l1 = std::hypot(B.x - A.x, B.y - A.y);
                    T l2 = std::hypot(C.x - B.x, C.y - B.y);
                    T sinA = (l1 > T(1e-12) && l2 > T(1e-12)) ? (cr / (l1 * l2)) : T(0);
                    // Only relax straight-chord angle threshold at hole anchors equipped with a curved bisector spoke
                    // so evaluateRingPatches inspects the actual curved Coons Jacobian.
                    T minReq = hasCurvedAnchor[vIdx[(k + 1) % 4]] ? T(-0.60) : T(1e-3);
                    if (sinA <= minReq) {
                        if (qErr) {
                            std::ostringstream oss;
                            oss << "q=" << idx << "/k=" << k << " sin=" << sinA
                                << " [(" << cp[0].x << "," << cp[0].y << "),("
                                << cp[1].x << "," << cp[1].y << "),("
                                << cp[2].x << "," << cp[2].y << "),("
                                << cp[3].x << "," << cp[3].y << ")]";
                            *qErr = oss.str();
                        }
                        return false;
                    }
                }
            }
            return true;
        };
        auto isFreeCandClear = [&](int v, const gsPoint2D<T>& cand) -> bool {
            if (pointInPolyStrict(cand, holePoly, T(1e-7))) return false;
            if (!pointInPolyStrict(cand, blockPoly, T(1e-7))) return false;
            for (int nb : adj[v]) {
                if (!isHoleSpokeClearOfHole(qg.pts[nb], cand, holePoly)) return false;
            }
            return true;
        };
        for (int iter = 0; iter < 15; ++iter) {
            for (int v = 0; v < NV; ++v) {
                if (qg.pinned[v] || adj[v].empty()) continue;
                gsPoint2D<T> avg(0, 0);
                for (int nb : adj[v]) { avg.x += qg.pts[nb].x; avg.y += qg.pts[nb].y; }
                avg.x /= adj[v].size(); avg.y /= adj[v].size();
                gsPoint2D<T> oldP = qg.pts[v];
                T oldMinS = minIncSin(v);
                bool moved = false;
                for (T alpha : {T(1.0), T(0.6), T(0.3)}) {
                    gsPoint2D<T> cand((T(1) - alpha) * oldP.x + alpha * avg.x,
                                      (T(1) - alpha) * oldP.y + alpha * avg.y);
                    qg.pts[v] = cand;
                    T newMinS = minIncSin(v);
                    if (((oldMinS > T(0.05) && newMinS > T(0.05)) || (newMinS > oldMinS + T(1e-4))) &&
                        isFreeCandClear(v, cand))
                    {
                        moved = true;
                        break;
                    }
                }
                if (!moved) qg.pts[v] = oldP;
            }
        }
        // Targeted directional + coordinate convexity/angle rescue if any free vertex touches a shallow/concave corner
        const T dirX[8] = { T(1), T(-1), T(0), T(0), T(0.70710678), T(0.70710678), T(-0.70710678), T(-0.70710678) };
        const T dirY[8] = { T(0), T(0), T(1), T(-1), T(0.70710678), T(-0.70710678), T(0.70710678), T(-0.70710678) };
        for (int pass = 0; pass < 4; ++pass) {
            for (int v = 0; v < NV; ++v) {
                if (qg.pinned[v] || adj[v].empty()) continue;
                T curMin = minIncSin(v);
                if (curMin >= T(0.28)) continue;
                gsPoint2D<T> oldP = qg.pts[v];
                gsPoint2D<T> bestP = oldP;
                T bestS = curMin;
                T sumLen = T(0);
                for (int nb : adj[v]) {
                    sumLen += std::hypot(qg.pts[nb].x - oldP.x, qg.pts[nb].y - oldP.y);
                    for (T frac : {T(0.15), T(0.30), T(0.45)}) {
                        gsPoint2D<T> cand((T(1) - frac) * oldP.x + frac * qg.pts[nb].x,
                                          (T(1) - frac) * oldP.y + frac * qg.pts[nb].y);
                        qg.pts[v] = cand;
                        T s = minIncSin(v);
                        if (s > bestS + T(1e-4) && isFreeCandClear(v, cand)) { bestS = s; bestP = cand; }
                    }
                }
                T avgLen = sumLen / static_cast<T>(adj[v].size());
                if (avgLen > T(1e-12)) {
                    for (T sf : {T(0.14), T(0.05)}) {
                        T st = sf * avgLen;
                        for (int d = 0; d < 8; ++d) {
                            gsPoint2D<T> cand(oldP.x + st * dirX[d], oldP.y + st * dirY[d]);
                            qg.pts[v] = cand;
                            T s = minIncSin(v);
                            if (s > bestS + T(1e-4) && isFreeCandClear(v, cand)) { bestS = s; bestP = cand; }
                        }
                    }
                }
                qg.pts[v] = bestP;
            }
        }
        if (!allQuadsConvex(failReason)) return false;
    }

    // Regenerate straight edges connected to relaxed vertices and emit Coons patches
    for (auto& qi : qg.items) {
        for (int m = 0; m < 4; ++m) {
            if (qi.straight[m])
                qi.e[m] = mkLine(qg.pts[qi.vid[m]], qg.pts[qi.vid[(m + 1) % 4]]);
        }
        auto r0 = qi.e[3].reversed(), r1 = qi.e[2].reversed();
        outPatches.emplace_back(makeCoonsPatchFrom4Curves(*r0.bsp, *r1.bsp, *qi.e[1].bsp, *qi.e[0].bsp, degree));
    }
    return true;
}

// Compute minimum positive distance t > 0 from origin along unit direction `dir` to polygon boundary.
template <class T>
inline T rayDistToPoly(
    const gsPoint2D<T>& origin,
    const gsPoint2D<T>& dir,
    const std::vector<gsPoint2D<T>>& poly)
{
    T bestT = T(1e18);
    const int N = static_cast<int>(poly.size());
    for (int i = 0; i < N; ++i) {
        const gsPoint2D<T>& A = poly[i];
        const gsPoint2D<T>& B = poly[(i + 1) % N];
        T ax = A.x - origin.x, ay = A.y - origin.y;
        T sx = B.x - A.x, sy = B.y - A.y;
        T denom = dir.x * sy - dir.y * sx;
        if (std::abs(denom) <= T(1e-14)) continue;
        T t = (ax * sy - ay * sx) / denom;
        T u = (ax * dir.y - ay * dir.x) / denom;
        if (t > T(1e-9) && u >= -T(1e-6) && u <= T(1) + T(1e-6))
            bestT = std::min(bestT, t);
    }
    return bestT;
}

// Blend interior CPs toward corner-bilinear map when a convex Coons patch slightly overshoots detJ.
template <class T>
inline bool bilinearBlendRescueRing(
    std::vector<std::unique_ptr<gsGeometry<T>>>& cand,
    const std::vector<gsPoint2D<T>>& holePoly,
    T& sc, T& md, T& ma,
    std::string* failMsg = nullptr)
{
    bool fixedAny = false;
    for (auto& gp : cand) {
        T mn = 0, mx = 0, ab = 0, en = 0;
        evalPatchJacobianAndWinslow(*gp, 15, mn, mx, ab, en);
        if (mn > T(-0.006) && mn <= T(1e-9)) {
            if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(gp.get())) {
                const int nu = tb->basis().component(0).size();
                const int nv = tb->basis().component(1).size();
                if (nu >= 3 && nv >= 3) {
                    gsMatrix<T> origCf = tb->coefs();
                    gsMatrix<T> gu = tb->basis().component(0).knots().greville();
                    gsMatrix<T> gv = tb->basis().component(1).knots().greville();
                    T P00[2] = {origCf(0, 0), origCf(0, 1)};
                    T P10[2] = {origCf(nu - 1, 0), origCf(nu - 1, 1)};
                    T P01[2] = {origCf((nv - 1) * nu, 0), origCf((nv - 1) * nu, 1)};
                    T P11[2] = {origCf(nu - 1 + (nv - 1) * nu, 0), origCf(nu - 1 + (nv - 1) * nu, 1)};
                    bool fixedThis = false;
                    for (T al : {T(0.25), T(0.50), T(0.75), T(1.0)}) {
                        gsMatrix<T> cf = origCf;
                        for (int j = 1; j < nv - 1; ++j) {
                            T eta = gv(0, j);
                            for (int i = 1; i < nu - 1; ++i) {
                                T xi = gu(0, i);
                                for (int d = 0; d < 2; ++d) {
                                    T bil = (1-xi)*(1-eta)*P00[d] + xi*(1-eta)*P10[d]
                                          + (1-xi)*eta*P01[d] + xi*eta*P11[d];
                                    cf(i + j * nu, d) = (1 - al) * origCf(i + j * nu, d) + al * bil;
                                }
                            }
                        }
                        tb->setCoefs(cf);
                        T tmn = 0, tmx = 0, tab = 0, ten = 0;
                        evalPatchJacobianAndWinslow(*tb, 15, tmn, tmx, tab, ten);
                        if (tmn > T(1e-9)) { fixedThis = true; fixedAny = true; break; }
                    }
                    if (!fixedThis) tb->setCoefs(origCf);
                }
            }
        }
    }
    if (!fixedAny) return false;
    return evaluateRingPatches(cand, holePoly, sc, md, ma, failMsg);
}

// Build an 8-quad conforming ring inside a 4-sided block [cS, cE, cN, cW] for a K=12 cross hole
// (4 concave ~90 deg notches spaced by 3 edges each: 4 notch quads + 4 arm-tip quads).
template <class T>
inline bool tryBuildCross8QuadRing(
    const CurveSeg<T>& cS, const CurveSeg<T>& cE,
    const CurveSeg<T>& cN, const CurveSeg<T>& cW,
    const std::vector<CurveSeg<T>>& hCCW,
    const std::vector<gsPoint2D<T>>& holePoly,
    std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
    T& bestScoreOut, T& bestMaxAngOut,
    int degree)
{
    if (hCCW.size() != 12) return false;
    const T pi = T(3.14159265358979323846);
    std::vector<int> notches;
    for (int i = 0; i < 12; ++i) {
        const auto& prev = hCCW[(i + 11) % 12];
        const auto& next = hCCW[i];
        gsPoint2D<T> d(next.p1().x - prev.p0().x, next.p1().y - prev.p0().y);
        if (std::hypot(d.x, d.y) <= T(1e-12)) d = gsPoint2D<T>(T(1), T(0));
        T aStart = 0, dCCW = 0, aDir = 0, minEdge = 0;
        if (computeExteriorCone(prev, next, d, aStart, dCCW, aDir, minEdge) &&
            dCCW < (T(150) * pi / T(180)))
        {
            notches.push_back(i);
        }
    }
    if (notches.size() != 4) return false;
    for (int s = 0; s < 4; ++s) {
        if ((notches[(s + 1) % 4] - notches[s] + 12) % 12 != 3)
            return false;
    }

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    const gsPoint2D<T> P[4] = { cS.p0(), cS.p1(), cE.p1(), cN.p0() };
    const CurveSeg<T> outerTop[4] = { cS, cE, cN.reversed(), cW.reversed() };

    bool found = false;
    bestScoreOut = T(1e18);
    bestMaxAngOut = T(180);

    for (int rot = 0; rot < 4; ++rot) {
        CurveSeg<T> eIn[4], eOut[4], eTip[4], spL[4], spR[4];
        bool clearOk = true;
        for (int s = 0; s < 4; ++s) {
            int i0 = notches[(s + rot) % 4];
            eIn[s]  = hCCW[(i0 + 11) % 12]; // T2_{s-1} -> V_s
            eOut[s] = hCCW[i0];             // V_s -> T1_s
            eTip[s] = hCCW[(i0 + 1) % 12];  // T1_s -> T2_s
            if (!isHoleSpokeClearOfHole(eIn[s].p0(),  P[s], holePoly) ||
                !isHoleSpokeClearOfHole(eOut[s].p1(), P[s], holePoly)) {
                clearOk = false;
                break;
            }
            spL[s] = mkLine(P[s], eIn[s].p0());
            spR[s] = mkLine(P[s], eOut[s].p1());
        }
        if (!clearOk) continue;

        std::vector<std::unique_ptr<gsGeometry<T>>> cand;
        cand.reserve(8);
        for (int s = 0; s < 4; ++s) {
            // Notch quad N_s: (P00=T1_s, P10=V_s, P11=T2_{s-1}, P01=P_s)
            cand.emplace_back(makeCoonsPatchFrom4Curves(
                *eOut[s].reversed().bsp, *eIn[s].reversed().bsp,
                *spL[s].bsp, *spR[s].reversed().bsp, degree));
            // Arm-tip quad A_s: (P00=P_s, P10=P_{s+1}, P11=T2_s, P01=T1_s)
            cand.emplace_back(makeCoonsPatchFrom4Curves(
                *outerTop[s].bsp, *spL[(s + 1) % 4].bsp,
                *eTip[s].bsp, *spR[s].bsp, degree));
        }

        T sc = 0, md = 0, ma = 0;
        if (evaluateRingPatches(cand, holePoly, sc, md, ma)) {
            if (!found || sc < bestScoreOut) {
                found = true;
                bestScoreOut = sc;
                bestMaxAngOut = ma;
                bestOut = std::move(cand);
            }
        }
    }
    return found;
}

// Build conforming 4-patch quad ring inside a 4-sided block [oS, oE, oN, oW] around an arbitrary hole.
template <class T>
inline void buildBlockQuadRing(
    const gsCurve<T>& oS, const gsCurve<T>& oE,
    const gsCurve<T>& oN, const gsCurve<T>& oW,
    const gsCurveLoop<T>& hole,
    gsMultiPatch<T>& mp,
    int degree)
{
    auto dist2 = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        T dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    };

    // Auto-orient the 4 outer block curves consistently:
    // cS: P0(SW)->P1(SE), cE: P1->P2(NE), cW: P0->P3(NW), cN: P3->P2
    CurveSeg<T> cS(oS, degree), cE(oE, degree), cN(oN, degree), cW(oW, degree);
    gsPoint2D<T> P0 = cS.p0(), P1 = cS.p1();
    if (dist2(cE.p0(), P1) > dist2(cE.p1(), P1)) cE = cE.reversed();
    gsPoint2D<T> P2 = cE.p1();
    if (dist2(cW.p0(), P0) > dist2(cW.p1(), P0)) cW = cW.reversed();
    gsPoint2D<T> P3 = cW.p1();
    if (dist2(cN.p0(), P3) > dist2(cN.p1(), P3)) cN = cN.reversed();
    const gsPoint2D<T> P[4] = { P0, P1, P2, P3 };

    // Orient hole curves strictly CCW and sample polygon for inside-hole checks
    std::vector<CurveSeg<T>> hCCW = orientLoopSegments(hole, degree, true);
    std::vector<gsPoint2D<T>> holePoly;
    holePoly.reserve(hCCW.size() * 12);
    for (const auto& seg : hCCW) {
        gsPoint2D<T> a = seg.p0(), b = seg.p1();
        gsPoint2D<T> pm = evalCurvePoint(*seg.bsp, T(0.5));
        gsPoint2D<T> pq = evalCurvePoint(*seg.bsp, T(0.25));
        T dMid = std::hypot(pm.x - T(0.5) * (a.x + b.x), pm.y - T(0.5) * (a.y + b.y));
        T dQtr = std::hypot(pq.x - (T(0.75) * a.x + T(0.25) * b.x), pq.y - (T(0.75) * a.y + T(0.25) * b.y));
        if (dMid < T(1e-6) && dQtr < T(1e-6)) {
            holePoly.push_back(a);
        } else {
            for (int k = 0; k < 12; ++k)
                holePoly.push_back(evalCurvePoint(*seg.bsp, static_cast<T>(k) / T(12)));
        }
    }
    const std::vector<T> hTurns = computeLoopVertexTurnDegs(hCCW);

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    // Helper: assemble and evaluate a 4-patch ring candidate from 4 sector chains inside an arbitrary 4-sided block
    auto tryBuildFrom4SectorsOnBlock = [&](const CurveSeg<T>& bcS, const CurveSeg<T>& bcE,
                                           const CurveSeg<T>& bcN, const CurveSeg<T>& bcW,
                                           const gsPoint2D<T> bP[4],
                                           const std::vector<CurveSeg<T>> sec[4],
                                           std::vector<std::unique_ptr<gsGeometry<T>>>& cand,
                                           bool& hasInteriorSharp, T& score, T& minDet, T& maxAng) -> bool {
        cand.clear();
        hasInteriorSharp = false;
        const T radToDeg = T(180) / T(3.14159265358979323846);
        for (int s = 0; s < 4; ++s) {
            if (sec[s].empty()) return false;
            for (size_t j = 1; j < sec[s].size(); ++j) {
                auto d0 = evalCurveDeriv(*sec[s][j - 1].bsp, T(1));
                auto d1 = evalCurveDeriv(*sec[s][j].bsp, T(0));
                T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
                if (n0 > T(1e-14) && n1 > T(1e-14)) {
                    T cr = (d0.x * d1.y - d0.y * d1.x) / (n0 * n1);
                    T dt = (d0.x * d1.x + d0.y * d1.y) / (n0 * n1);
                    if (std::abs(std::atan2(cr, dt)) * radToDeg > T(28)) {
                        hasInteriorSharp = true;
                        return false;
                    }
                }
            }
        }
        gsPoint2D<T> H[4];
        for (int s = 0; s < 4; ++s) H[s] = sec[s].front().p0();
        for (int s = 0; s < 4; ++s) {
            if (!isHoleSpokeClearOfHole(H[s], bP[s], holePoly))
                return false;
        }
        std::shared_ptr<gsBSpline<T>> hJoin[4];
        for (int s = 0; s < 4; ++s) hJoin[s] = joinCurveSegsExact(sec[s], degree);

        // hJoin[0]: H0->H1 (S north), hJoin[1]: H1->H2 (E west),
        // hJoin[2]: H2->H3 (needs reversed -> H3->H2 for N south),
        // hJoin[3]: H3->H0 (needs reversed -> H0->H3 for W east)
        std::unique_ptr<gsBSpline<T>> hN(new gsBSpline<T>(*hJoin[2])); hN->reverse(); hN->knots().affineTransformTo(T(0), T(1));
        std::unique_ptr<gsBSpline<T>> hW(new gsBSpline<T>(*hJoin[3])); hW->reverse(); hW->knots().affineTransformTo(T(0), T(1));

        CurveSeg<T> spOut[4], spIn[4];
        for (int s = 0; s < 4; ++s) {
            bool wc = false;
            spOut[s] = makeExteriorConeSafeSpoke(sec[(s + 3) % 4].back(), sec[s].front(), H[s], bP[s], holePoly, degree, wc);
            spIn[s]  = spOut[s].reversed();
        }

        cand.emplace_back(makeCoonsPatchFrom4Curves(*bcS.bsp, *spIn[1].bsp, *hJoin[0], *spIn[0].bsp, degree));
        cand.emplace_back(makeCoonsPatchFrom4Curves(*spOut[1].bsp, *bcE.bsp, *spOut[2].bsp, *hJoin[1], degree));
        cand.emplace_back(makeCoonsPatchFrom4Curves(*hN, *spOut[2].bsp, *bcN.bsp, *spOut[3].bsp, degree));
        cand.emplace_back(makeCoonsPatchFrom4Curves(*spIn[0].bsp, *hW, *spIn[3].bsp, *bcW.bsp, degree));
        std::string fMsg;
        bool ok = evaluateRingPatches(cand, holePoly, score, minDet, maxAng, &fMsg);
        if (!ok)
            ok = bilinearBlendRescueRing(cand, holePoly, score, minDet, maxAng, &fMsg);
        if (!ok && std::getenv("GISMO_SPLITTER_DEBUG") && hCCW.size() <= 4)
            std::cerr << "    [4SEC fail K=" << hCCW.size() << "] " << fMsg << "\n";
        return ok;
    };
    auto tryBuildFrom4Sectors = [&](const std::vector<CurveSeg<T>> sec[4],
                                    std::vector<std::unique_ptr<gsGeometry<T>>>& cand,
                                    bool& hasInteriorSharp, T& score, T& minDet, T& maxAng) -> bool {
        return tryBuildFrom4SectorsOnBlock(cS, cE, cN, cW, P, sec, cand, hasInteriorSharp, score, minDet, maxAng);
    };

    std::vector<std::unique_ptr<gsGeometry<T>>> bestSmooth4;
    T bestSmoothScore = T(1e18), bestSmoothMaxAng = T(180);
    bool foundSmooth4 = false;

    const int K = static_cast<int>(hCCW.size());
    auto considerSplits = [&](const std::vector<CurveSeg<T>>& ring) {
        const int M = static_cast<int>(ring.size());
        if (M < 4) return;
        if (M == 4) {
            for (int rot = 0; rot < 4; ++rot) {
                std::vector<CurveSeg<T>> sec[4];
                for (int s = 0; s < 4; ++s) sec[s] = { ring[(rot + s) % 4] };
                std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                bool hasSharp = false; T sc = 0, md = 0, ma = 0;
                if (tryBuildFrom4Sectors(sec, cand, hasSharp, sc, md, ma)) {
                    if (!hasSharp && (!foundSmooth4 || sc < bestSmoothScore)) {
                        foundSmooth4 = true; bestSmoothScore = sc; bestSmoothMaxAng = ma; bestSmooth4 = std::move(cand);
                    }
                }
            }
            return;
        }
        // Choose 4 anchor indices aligned with P[0..3] rays from hole centroid + neighbor shifts
        gsPoint2D<T> holeCtr(0, 0);
        for (const auto& hp : holePoly) { holeCtr.x += hp.x; holeCtr.y += hp.y; }
        holeCtr.x /= holePoly.size(); holeCtr.y /= holePoly.size();
        int baseAnchor[4];
        for (int s = 0; s < 4; ++s) {
            T angP = std::atan2(P[s].y - holeCtr.y, P[s].x - holeCtr.x);
            T bestDiff = T(1e9); int bestIdx = 0;
            for (int i = 0; i < M; ++i) {
                gsPoint2D<T> v = ring[i].p0();
                T a = std::atan2(v.y - holeCtr.y, v.x - holeCtr.x);
                T d = std::abs(a - angP);
                if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
                if (d < bestDiff) { bestDiff = d; bestIdx = i; }
            }
            baseAnchor[s] = bestIdx;
        }
        // Try small offsets around baseAnchor while preserving strictly increasing cyclic order
        for (int d0 = -1; d0 <= 1; ++d0)
        for (int d1 = -1; d1 <= 1; ++d1)
        for (int d2 = -1; d2 <= 1; ++d2)
        for (int d3 = -1; d3 <= 1; ++d3) {
            int a[4] = {
                (baseAnchor[0] + d0 + M) % M, (baseAnchor[1] + d1 + M) % M,
                (baseAnchor[2] + d2 + M) % M, (baseAnchor[3] + d3 + M) % M
            };
            // Check distinct and CCW cyclic ordering
            bool distinct = true;
            for (int x = 0; x < 4; ++x)
                for (int y = x + 1; y < 4; ++y)
                    if (a[x] == a[y]) distinct = false;
            if (!distinct) continue;
            int stepSum = 0;
            std::vector<CurveSeg<T>> sec[4];
            for (int s = 0; s < 4; ++s) {
                int len = (a[(s + 1) % 4] - a[s] + M) % M;
                if (len <= 0) { stepSum = -1; break; }
                stepSum += len;
                for (int j = 0; j < len; ++j) sec[s].push_back(ring[(a[s] + j) % M]);
            }
            if (stepSum != M) continue;
            std::vector<std::unique_ptr<gsGeometry<T>>> cand;
            bool hasSharp = false; T sc = 0, md = 0, ma = 0;
            if (tryBuildFrom4Sectors(sec, cand, hasSharp, sc, md, ma)) {
                if (!hasSharp && (!foundSmooth4 || sc < bestSmoothScore)) {
                    foundSmooth4 = true; bestSmoothScore = sc; bestSmoothMaxAng = ma; bestSmooth4 = std::move(cand);
                }
            }
        }
    };

    if (K == 1) {
        // Split single closed curve into 4 equal-arc segments starting near angle toward P0
        gsPoint2D<T> holeCtr(0, 0);
        for (const auto& hp : holePoly) { holeCtr.x += hp.x; holeCtr.y += hp.y; }
        holeCtr.x /= holePoly.size(); holeCtr.y /= holePoly.size();
        T targetAng = std::atan2(P0.y - holeCtr.y, P0.x - holeCtr.x);
        T bestU = T(0), bestDiff = T(1e9);
        for (int k = 0; k < 64; ++k) {
            T u = static_cast<T>(k) / T(64);
            gsPoint2D<T> q = evalCurvePoint(*hCCW[0].bsp, u);
            T d = std::abs(std::atan2(q.y - holeCtr.y, q.x - holeCtr.x) - targetAng);
            if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
            if (d < bestDiff) { bestDiff = d; bestU = u; }
        }
        std::vector<CurveSeg<T>> r4;
        CurveSeg<T> shifted = hCCW[0];
        if (bestU > T(1e-4) && bestU < T(1) - T(1e-4)) {
            auto sp = hCCW[0].splitAtParam(bestU);
            auto joined = joinCurveSegsExact(std::vector<CurveSeg<T>>{sp.second, sp.first}, degree);
            shifted = CurveSeg<T>(*joined);
        }
        CurveSeg<T> rem = shifted;
        for (int k = 0; k < 3; ++k) {
            auto sp = rem.splitArc(T(1) / static_cast<T>(4 - k));
            r4.push_back(sp.first);
            rem = sp.second;
        }
        r4.push_back(rem);
        considerSplits(r4);
    } else if (K == 2) {
        std::vector<CurveSeg<T>> r4;
        for (int i = 0; i < 2; ++i) {
            auto sp = hCCW[i].splitArc(T(0.5));
            r4.push_back(sp.first);
            r4.push_back(sp.second);
        }
        considerSplits(r4);
    } else if (K == 3) {
        // Split each single edge choice in turn (at t=0.5) -> 4-segment ring keeping all 3 sharp corners as vertices
        for (int eSplit = 0; eSplit < 3; ++eSplit) {
            std::vector<CurveSeg<T>> r4;
            for (int i = 0; i < 3; ++i) {
                if (i == eSplit) {
                    auto sp = hCCW[i].splitArc(T(0.5));
                    r4.push_back(sp.first);
                    r4.push_back(sp.second);
                } else {
                    r4.push_back(hCCW[i]);
                }
            }
            considerSplits(r4);
        }
    } else if (K == 4) {
        considerSplits(hCCW);
        // Also test midpoint-shifted 4-ring if hole corners are smooth (<= 25 deg turns)
        bool allSmooth = true;
        for (T t : hTurns) if (t > T(25)) { allSmooth = false; break; }
        if (allSmooth) {
            std::vector<CurveSeg<T>> halves;
            for (int i = 0; i < 4; ++i) {
                auto sp = hCCW[i].splitArc(T(0.5));
                halves.push_back(sp.first);
                halves.push_back(sp.second);
            }
            std::vector<CurveSeg<T>> rot4;
            for (int i = 0; i < 4; ++i) {
                auto j = joinCurveSegsExact(std::vector<CurveSeg<T>>{halves[(2 * i + 1) % 8], halves[(2 * i + 2) % 8]}, degree);
                rot4.emplace_back(*j);
            }
            considerSplits(rot4);
            // Also test radial angle-aligned 4-way split toward block corners P[0..3]
            auto joinedFull = joinCurveSegsExact(hCCW, degree);
            gsPoint2D<T> hc(0, 0);
            for (const auto& hp : holePoly) { hc.x += hp.x; hc.y += hp.y; }
            hc.x /= holePoly.size(); hc.y /= holePoly.size();
            const int NS = 128;
            std::vector<T> sampleU(NS), sampleAng(NS);
            for (int k = 0; k < NS; ++k) {
                sampleU[k] = static_cast<T>(k) / static_cast<T>(NS);
                gsPoint2D<T> q = evalCurvePoint(*joinedFull, sampleU[k]);
                sampleAng[k] = std::atan2(q.y - hc.y, q.x - hc.x);
            }
            T cutU[4];
            for (int s = 0; s < 4; ++s) {
                T tAng = std::atan2(P[s].y - hc.y, P[s].x - hc.x);
                T bestD = T(1e9); int bi = 0;
                for (int k = 0; k < NS; ++k) {
                    T d = std::abs(sampleAng[k] - tAng);
                    if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
                    if (d < bestD) { bestD = d; bi = k; }
                }
                cutU[s] = sampleU[bi];
            }
            // Shift parameter origin to cutU[0] and sort remaining cuts
            CurveSeg<T> fullSeg(*joinedFull);
            if (cutU[0] > T(1e-4) && cutU[0] < T(1) - T(1e-4)) {
                auto sp0 = fullSeg.splitAtParam(cutU[0]);
                fullSeg = CurveSeg<T>(*joinCurveSegsExact(std::vector<CurveSeg<T>>{sp0.second, sp0.first}, degree));
            }
            // Locate P[1..3] along shifted fullSeg and split at exact radial params
            T uShift[3];
            bool uOk = true;
            for (int s = 1; s < 4; ++s) {
                T tAng = std::atan2(P[s].y - hc.y, P[s].x - hc.x);
                T bestD = T(1e9), bu = T(0);
                for (int k = 1; k < NS; ++k) {
                    T u = static_cast<T>(k) / static_cast<T>(NS);
                    gsPoint2D<T> q = evalCurvePoint(*fullSeg.bsp, u);
                    T d = std::abs(std::atan2(q.y - hc.y, q.x - hc.x) - tAng);
                    if (d > T(3.14159265358979323846)) d = T(2) * T(3.14159265358979323846) - d;
                    if (d < bestD) { bestD = d; bu = u; }
                }
                uShift[s - 1] = bu;
            }
            if (uShift[0] > T(0.05) && uShift[0] < uShift[1] - T(0.05) &&
                uShift[1] < uShift[2] - T(0.05) && uShift[2] < T(0.95))
            {
                auto sp0 = fullSeg.splitAtParam(uShift[0]);
                T u1l = (uShift[1] - uShift[0]) / (T(1) - uShift[0]);
                auto sp1 = sp0.second.splitAtParam(u1l);
                T u2l = (uShift[2] - uShift[1]) / (T(1) - uShift[1]);
                auto sp2 = sp1.second.splitAtParam(u2l);
                std::vector<CurveSeg<T>> rRadial = {sp0.first, sp1.first, sp2.first, sp2.second};
                {
                    std::vector<CurveSeg<T>> sec[4];
                    for (int s = 0; s < 4; ++s) sec[s] = { rRadial[s] };
                    std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                    bool hasSharp = false; T sc = 0, md = 0, ma = 0;
                    if (tryBuildFrom4Sectors(sec, cand, hasSharp, sc, md, ma)) {
                        if (!hasSharp && (!foundSmooth4 || sc < bestSmoothScore)) {
                            foundSmooth4 = true; bestSmoothScore = sc; bestSmoothMaxAng = ma; bestSmooth4 = std::move(cand);
                        }
                    }
                }
            }
            (void)uOk;
        }
    } else {
        // K >= 5: first test 4-sector ring (succeeds as smooth 4-patch when turns are small/smooth)
        considerSplits(hCCW);
    }

    if (foundSmooth4 && bestSmoothMaxAng <= T(163.5)) {
        if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSmooth4\n";
        for (auto& p : bestSmooth4)
            mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
        return;
    }

    gsPoint2D<T> holeCtr(0, 0);
    for (const auto& hp : holePoly) { holeCtr.x += hp.x; holeCtr.y += hp.y; }
    if (!holePoly.empty()) { holeCtr.x /= holePoly.size(); holeCtr.y /= holePoly.size(); }
    T holeR = T(0);
    for (const auto& hp : holePoly)
        holeR = std::max(holeR, std::hypot(hp.x - holeCtr.x, hp.y - holeCtr.y));
    T maxBlockDist = T(0);
    for (int s = 0; s < 4; ++s)
        maxBlockDist = std::max(maxBlockDist, std::hypot(P[s].x - holeCtr.x, P[s].y - holeCtr.y));

    auto runSharpSearchOnBlock = [&](const CurveSeg<T>& bcS, const CurveSeg<T>& bcE,
                                     const CurveSeg<T>& bcN, const CurveSeg<T>& bcW,
                                     const gsPoint2D<T> bP[4],
                                     const std::vector<CurveSeg<T>>& workH,
                                     const std::vector<T>& splitTs,
                                     const std::vector<T>& oWeights,
                                     bool earlyStop,
                                     std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
                                     T& bestScoreOut, T& bestMaxAngOut,
                                     int* dbgVisOut = nullptr, int* dbgCapOut = nullptr, int* dbgEvalOut = nullptr,
                                     std::string* firstCapOut = nullptr, std::string* firstFailOut = nullptr) -> bool
    {
        const int M = static_cast<int>(workH.size());
        std::vector<std::vector<bool>> vis(4, std::vector<bool>(M, false));
        for (int s = 0; s < 4; ++s)
            for (int i = 0; i < M; ++i)
                vis[s][i] = isHoleSpokeClearOfHole(workH[i].p0(), bP[s], holePoly);

        bestScoreOut = T(1e18);
        bestMaxAngOut = T(180);
        bool found = false;
        const int maxL = (M >= 8) ? 4 : 3;
        for (int a0 = 0; a0 < M; ++a0) {
            if (!vis[0][a0]) continue;
            for (int l0 = 1; l0 <= maxL; ++l0) {
                int a1 = (a0 + l0) % M;
                if (!vis[1][a1]) continue;
                for (int l1 = 1; l1 <= maxL; ++l1) {
                    int a2 = (a1 + l1) % M;
                    if (!vis[2][a2]) continue;
                    for (int l2 = 1; l2 <= maxL; ++l2) {
                        int l3 = M - (l0 + l1 + l2);
                        if (l3 < 1 || l3 > maxL) continue;
                        int a3 = (a2 + l2) % M;
                        if (!vis[3][a3]) continue;
                        int a[4] = {a0, a1, a2, a3};
                        int maxSec = std::max(std::max(l0, l1), std::max(l2, l3));
                        if (dbgVisOut) ++(*dbgVisOut);

                        for (T st : splitTs) {
                            for (T ow : oWeights) {
                                std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                                std::string cmsg;
                                if (!buildSharpCornerBlockRing(bcS, bcE, bcN, bcW, workH, a, holePoly, cand, degree, ow, &cmsg, st)) {
                                    if (firstCapOut && (firstCapOut->empty() || (maxSec == 2 && firstCapOut->find("max2:") == std::string::npos)))
                                        *firstCapOut = (maxSec == 2 ? "max2:" : "") + cmsg;
                                    if (dbgCapOut) ++(*dbgCapOut);
                                    continue;
                                }
                                T sc = 0, md = 0, ma = 0;
                                std::string fmsg;
                                bool pass = evaluateRingPatches(cand, holePoly, sc, md, ma, (firstFailOut && firstFailOut->empty()) ? &fmsg : nullptr);
                                if (!pass)
                                    pass = bilinearBlendRescueRing(cand, holePoly, sc, md, ma, (firstFailOut && firstFailOut->empty()) ? &fmsg : nullptr);
                                if (pass) {
                                    T totalMetric = sc + T(5) * static_cast<T>(maxSec) + (ma > T(163.5) ? T(1e5) : T(0));
                                    if (!found || totalMetric < bestScoreOut) {
                                        found = true; bestScoreOut = totalMetric; bestMaxAngOut = ma; bestOut = std::move(cand);
                                    }
                                    if (earlyStop && bestMaxAngOut <= T(154))
                                        return true;
                                } else {
                                    if (firstFailOut && firstFailOut->empty()) *firstFailOut = fmsg;
                                    if (dbgEvalOut) ++(*dbgEvalOut);
                                }
                            }
                        }
                    }
                }
            }
        }
        return found;
    };

    auto tryBuildBufferedRing = [&](std::vector<std::unique_ptr<gsGeometry<T>>>& bestOut,
                                    T& bestMaxAngOut,
                                    bool allowSharpFallback = true) -> bool
    {
        if (holeR <= T(1e-12)) return false;
        const T pi = T(3.14159265358979323846);
        std::vector<gsPoint2D<T>> blockPoly;
        blockPoly.reserve(32);
        const CurveSeg<T> outerSegs[4] = { cS, cE, cN.reversed(), cW.reversed() };
        for (int s = 0; s < 4; ++s) {
            for (int k = 0; k < 8; ++k)
                blockPoly.push_back(evalCurvePoint(*outerSegs[s].bsp, static_cast<T>(k) / T(8)));
        }
        if (!pointInPolyStrict(holeCtr, blockPoly, T(1e-7))) return false;

        T angP[4];
        for (int s = 0; s < 4; ++s)
            angP[s] = std::atan2(P[s].y - holeCtr.y, P[s].x - holeCtr.x);

        T angSq[4];
        bool haveNotchSq = false;
        if (K == 12) {
            std::vector<int> notches;
            for (int i = 0; i < 12; ++i) {
                const auto& prev = hCCW[(i + 11) % 12];
                const auto& next = hCCW[i];
                gsPoint2D<T> d(next.p1().x - prev.p0().x, next.p1().y - prev.p0().y);
                T aStart = 0, dCCW = 0, aDir = 0, minEdge = 0;
                if (computeExteriorCone(prev, next, d, aStart, dCCW, aDir, minEdge) &&
                    dCCW < (T(150) * pi / T(180)))
                    notches.push_back(i);
            }
            if (notches.size() == 4) {
                T nAng[4];
                for (int s = 0; s < 4; ++s) {
                    gsPoint2D<T> v = hCCW[notches[s]].p0();
                    nAng[s] = std::atan2(v.y - holeCtr.y, v.x - holeCtr.x);
                }
                int bestRot = 0;
                T bestAlign = -T(1e9);
                for (int r = 0; r < 4; ++r) {
                    T al = T(0);
                    for (int s = 0; s < 4; ++s) al += std::cos(nAng[(s + r) % 4] - angP[s]);
                    if (al > bestAlign) { bestAlign = al; bestRot = r; }
                }
                for (int s = 0; s < 4; ++s) angSq[s] = nAng[(s + bestRot) % 4];
                haveNotchSq = true;
            }
        }
        if (!haveNotchSq) {
            T sumC = T(0), sumS = T(0);
            for (int s = 0; s < 4; ++s) {
                T a = angP[s] - static_cast<T>(s) * (pi / T(2));
                sumC += std::cos(a);
                sumS += std::sin(a);
            }
            T baseA = std::atan2(sumS, sumC);
            for (int s = 0; s < 4; ++s)
                angSq[s] = baseA + static_cast<T>(s) * (pi / T(2));
        }

        auto wrapPi = [&](T a) -> T {
            while (a >  pi) a -= T(2) * pi;
            while (a < -pi) a += T(2) * pi;
            return a;
        };

        bool foundBuf = false;
        size_t bestCount = std::numeric_limits<size_t>::max();
        T bestScore = T(1e18);
        bestMaxAngOut = T(180);

        auto evalBox = [&](const gsPoint2D<T> B[4]) -> bool {
            for (int s = 0; s < 4; ++s) {
                if (pointInPolyStrict(B[s], holePoly, T(1e-6)) ||
                    !pointInPolyStrict(B[s], blockPoly, T(1e-6)) ||
                    !isHoleSpokeClearOfHole(B[s], P[s], holePoly))
                    return false;
            }
            for (int s = 0; s < 4; ++s) {
                if (cross2(B[s], B[(s + 1) % 4], B[(s + 2) % 4]) <= T(1e-10))
                    return false;
                gsPoint2D<T> mid(T(0.5) * (B[s].x + B[(s + 1) % 4].x),
                                 T(0.5) * (B[s].y + B[(s + 1) % 4].y));
                if (!pointInPolyStrict(mid, blockPoly, T(1e-6)) ||
                    !isHoleSpokeClearOfHole(B[s], B[(s + 1) % 4], holePoly))
                    return false;
            }

            CurveSeg<T> bS = mkLine(B[0], B[1]);
            CurveSeg<T> bE = mkLine(B[1], B[2]);
            CurveSeg<T> bN = mkLine(B[3], B[2]);
            CurveSeg<T> bW = mkLine(B[0], B[3]);

            std::vector<std::unique_ptr<gsGeometry<T>>> outer4;
            outer4.reserve(4);
            CurveSeg<T> spOut[4], spIn[4];
            for (int s = 0; s < 4; ++s) {
                spOut[s] = mkLine(B[s], P[s]);
                spIn[s]  = spOut[s].reversed();
            }
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*cS.bsp, *spIn[1].bsp, *bS.bsp, *spIn[0].bsp, degree));
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*spOut[1].bsp, *cE.bsp, *spOut[2].bsp, *bE.bsp, degree));
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*bN.bsp, *spOut[2].bsp, *cN.bsp, *spOut[3].bsp, degree));
            outer4.emplace_back(makeCoonsPatchFrom4Curves(*spIn[0].bsp, *bW.bsp, *spIn[3].bsp, *cW.bsp, degree));

            T outSc = 0, outMd = 0, outMa = 0;
            if (!evaluateRingPatches(outer4, holePoly, outSc, outMd, outMa))
                return false;

            std::vector<std::unique_ptr<gsGeometry<T>>> innerCand;
            T inSc = 0, inMa = T(180);
            bool inOk = false;

            if (K == 12) {
                inOk = tryBuildCross8QuadRing(bS, bE, bN, bW, hCCW, holePoly, innerCand, inSc, inMa, degree);
            } else if (K == 3 || K == 4) {
                auto tryRing4 = [&](const std::vector<CurveSeg<T>>& r4) {
                    for (int rot = 0; rot < 4; ++rot) {
                        std::vector<CurveSeg<T>> sec[4];
                        for (int s = 0; s < 4; ++s) sec[s] = { r4[(rot + s) % 4] };
                        std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                        bool hs = false; T sc = 0, md = 0, ma = 0;
                        if (tryBuildFrom4SectorsOnBlock(bS, bE, bN, bW, B, sec, cand, hs, sc, md, ma) && !hs) {
                            if (!inOk || sc < inSc) {
                                inOk = true; inSc = sc; inMa = ma; innerCand = std::move(cand);
                            }
                        }
                    }
                };
                if (K == 4) {
                    tryRing4(hCCW);
                } else {
                    for (int eSplit = 0; eSplit < 3; ++eSplit) {
                        std::vector<CurveSeg<T>> r4;
                        for (int i = 0; i < 3; ++i) {
                            if (i == eSplit) {
                                auto sp = hCCW[i].splitArc(T(0.5));
                                r4.push_back(sp.first);
                                r4.push_back(sp.second);
                            } else {
                                r4.push_back(hCCW[i]);
                            }
                        }
                        tryRing4(r4);
                    }
                }
            }
            if (!inOk && allowSharpFallback && K >= 3) {
                std::vector<CurveSeg<T>> wH;
                if (K == 3 || K == 4) {
                    for (const auto& s : hCCW) {
                        auto sp = s.splitArc(T(0.5));
                        wH.push_back(sp.first);
                        wH.push_back(sp.second);
                    }
                } else {
                    wH = hCCW;
                    if (wH.size() % 2 != 0) {
                        int longestIdx = 0;
                        T bestL = -T(1);
                        for (size_t i = 0; i < wH.size(); ++i) {
                            T d = std::hypot(wH[i].p1().x - wH[i].p0().x, wH[i].p1().y - wH[i].p0().y);
                            if (d > bestL) { bestL = d; longestIdx = static_cast<int>(i); }
                        }
                        auto sp = wH[longestIdx].splitArc(T(0.5));
                        wH[longestIdx] = sp.first;
                        wH.insert(wH.begin() + longestIdx + 1, sp.second);
                    }
                }
                inOk = runSharpSearchOnBlock(
                    bS, bE, bN, bW, B, wH,
                    {T(0.5)}, {T(0.62), T(0.45), T(0.78)}, true,
                    innerCand, inSc, inMa);
            }

            if (inOk) {
                size_t totN = outer4.size() + innerCand.size();
                T totMa = std::max(outMa, inMa);
                T totSc = outSc + inSc;
                if (!foundBuf || totN < bestCount || (totN == bestCount && totSc < bestScore)) {
                    foundBuf = true;
                    bestCount = totN;
                    bestScore = totSc;
                    bestMaxAngOut = totMa;
                    bestOut.clear();
                    bestOut.reserve(totN);
                    for (auto& p : outer4) bestOut.push_back(std::move(p));
                    for (auto& p : innerCand) bestOut.push_back(std::move(p));
                }
                if (totMa <= T(156) && ((K == 12 && totN <= 12) || (K <= 4 && totN <= 8) || (K >= 5 && totN <= 20)))
                    return true;
            }
            return false;
        };

        const std::vector<T> scaleList = allowSharpFallback
            ? std::vector<T>{T(1.90), T(2.25), T(1.60)}
            : std::vector<T>{T(1.90), T(2.25), T(1.60), T(2.65)};
        const std::vector<T> blendList = allowSharpFallback
            ? std::vector<T>{T(0.0), T(0.50)}
            : std::vector<T>{T(0.0), T(0.35), T(0.70), T(1.0)};

        for (T scale : scaleList) {
            for (T blend : blendList) {
                gsPoint2D<T> B[4];
                bool raysOk = true;
                for (int s = 0; s < 4; ++s) {
                    T th = angSq[s] + blend * wrapPi(angP[s] - angSq[s]);
                    gsPoint2D<T> dir(std::cos(th), std::sin(th));
                    T dWall = rayDistToPoly(holeCtr, dir, blockPoly);
                    if (dWall > T(1e9)) { raysOk = false; break; }
                    T dist = std::min(scale * holeR, T(0.76) * dWall);
                    B[s] = gsPoint2D<T>(holeCtr.x + dist * dir.x, holeCtr.y + dist * dir.y);
                }
                if (raysOk && evalBox(B))
                    return true;
            }
        }
        // Also test block-proportioned lerp boxes for tall/skewed cells where 1:1 boxes fail outer4 aspect checks
        for (T al : {T(0.48), T(0.58), T(0.38)}) {
            gsPoint2D<T> B[4];
            for (int s = 0; s < 4; ++s)
                B[s] = gsPoint2D<T>((T(1) - al) * holeCtr.x + al * P[s].x,
                                    (T(1) - al) * holeCtr.y + al * P[s].y);
            if (evalBox(B))
                return true;
        }
        return foundBuf;
    };

    // 1. Direct 8-quad cross ring for K=12 crosses when enclosing block has clean visibility & angles
    if (K == 12) {
        std::vector<std::unique_ptr<gsGeometry<T>>> cross8;
        T cSc = 0, cMa = T(180);
        if (tryBuildCross8QuadRing(cS, cE, cN, cW, hCCW, holePoly, cross8, cSc, cMa, degree) &&
            cMa <= T(163.5))
        {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked direct Cross8Quad P=8\n";
            for (auto& p : cross8)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
    }

    // 2. For K=3,4 when direct 4-sector ring failed or exceeded 163.5 deg, try compact 8-patch buffered ring before edge doubling
    if (K == 3 || K == 4) {
        std::vector<std::unique_ptr<gsGeometry<T>>> buf8;
        T bMa = T(180);
        if (tryBuildBufferedRing(buf8, bMa, false) && buf8.size() <= 8 && bMa <= T(163.5)) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked buffered K=" << K << " P=" << buf8.size() << "\n";
            for (auto& p : buf8)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
    }

    // 3. For K=12 crosses on tall/skewed blocks, try 12-patch buffered cross ring before general sharp search
    if (K == 12) {
        std::vector<std::unique_ptr<gsGeometry<T>>> buf12;
        T bMa = T(180);
        if (tryBuildBufferedRing(buf12, bMa, false) && buf12.size() <= 12 && bMa <= T(163.5)) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked buffered Cross P=" << buf12.size() << "\n";
            for (auto& p : buf12)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
    }

    // 4. For K >= 5 with sharp corners (or K=3,4 fallback), build sharp/midpoint-aligned structured ring
    if (K >= 3) {
        std::vector<CurveSeg<T>> workH;
        if (K == 3 || K == 4) {
            for (const auto& s : hCCW) {
                auto sp = s.splitArc(T(0.5));
                workH.push_back(sp.first);
                workH.push_back(sp.second);
            }
        } else {
            workH = hCCW;
            if (workH.size() % 2 != 0) {
                int longestIdx = 0;
                T bestL = -T(1);
                for (size_t i = 0; i < workH.size(); ++i) {
                    T d = std::hypot(workH[i].p1().x - workH[i].p0().x, workH[i].p1().y - workH[i].p0().y);
                    if (d > bestL) { bestL = d; longestIdx = static_cast<int>(i); }
                }
                auto sp = workH[longestIdx].splitArc(T(0.5));
                workH[longestIdx] = sp.first;
                workH.insert(workH.begin() + longestIdx + 1, sp.second);
            }
        }

        std::vector<std::unique_ptr<gsGeometry<T>>> bestSharp;
        T bestSharpScore = T(1e18), bestSharpMaxAng = T(180);
        int dbgVisCount = 0, dbgCapFail = 0, dbgEvalFail = 0;
        std::string firstCapFail, firstFail;
        bool foundSharp = runSharpSearchOnBlock(
            cS, cE, cN, cW, P, workH, {T(0.5)},
            {T(0.45), T(0.62), T(0.78), T(0.90)}, false,
            bestSharp, bestSharpScore, bestSharpMaxAng,
            &dbgVisCount, &dbgCapFail, &dbgEvalFail, &firstCapFail, &firstFail);
        if (!foundSharp && holeR > T(1e-12) && maxBlockDist > T(2.2) * holeR) {
            foundSharp = runSharpSearchOnBlock(
                cS, cE, cN, cW, P, workH, {T(0.25)},
                {T(0.45), T(0.62), T(0.78), T(0.90)}, false,
                bestSharp, bestSharpScore, bestSharpMaxAng,
                &dbgVisCount, &dbgCapFail, &dbgEvalFail, &firstCapFail, &firstFail);
        }

        if (foundSharp && (!foundSmooth4 || bestSharpMaxAng <= T(163.5) || bestSharpMaxAng <= bestSmoothMaxAng)) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSharp P=" << bestSharp.size() << "\n";
            for (auto& p : bestSharp)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
        {
            std::vector<std::unique_ptr<gsGeometry<T>>> bufBest;
            T bufMaxAng = T(180);
            if (tryBuildBufferedRing(bufBest, bufMaxAng, true) &&
                (!foundSmooth4 || bufMaxAng <= T(163.5) || bufMaxAng <= bestSmoothMaxAng))
            {
                if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked buffered fallback P=" << bufBest.size() << "\n";
                for (auto& p : bufBest)
                    mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
                return;
            }
        }
        if (foundSmooth4) {
            if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSmooth4 (backup)\n";
            for (auto& p : bestSmooth4)
                mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
            return;
        }
        if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] foundSharp=false for K=" << K << " (visQuad=" << dbgVisCount << ", capFail=" << dbgCapFail << "[" << firstCapFail << "], evalFail=" << dbgEvalFail << "[" << firstFail << "])\n";
    }

    if (foundSmooth4) {
        if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] picked foundSmooth4 (backup)\n";
        for (auto& p : bestSmooth4)
            mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));
        return;
    }

    // Ultimate fallback (only reached if all geometric spoke checks failed):
    // Still use exact joined hole curves instead of straight lines!
    if (std::getenv("GISMO_SPLITTER_DEBUG")) std::cerr << "  [RING] HIT ULTIMATE FALLBACK for K=" << K << "\n";
    (void)holeCtr;
    std::vector<CurveSeg<T>> workH = hCCW;
    while (workH.size() < 4) {
        std::vector<CurveSeg<T>> nextH;
        for (auto& s : workH) {
            auto sp = s.splitArc(T(0.5));
            nextH.push_back(sp.first);
            nextH.push_back(sp.second);
        }
        workH = std::move(nextH);
    }
    const int M = static_cast<int>(workH.size());
    int anchors[4] = {0, M / 4, M / 2, (3 * M) / 4};
    std::vector<CurveSeg<T>> sec[4];
    for (int s = 0; s < 4; ++s) {
        int len = (anchors[(s + 1) % 4] - anchors[s] + M) % M;
        for (int j = 0; j < len; ++j) sec[s].push_back(workH[(anchors[s] + j) % M]);
    }
    std::shared_ptr<gsBSpline<T>> hJ[4];
    gsPoint2D<T> H[4];
    for (int s = 0; s < 4; ++s) {
        hJ[s] = joinCurveSegsExact(sec[s], degree);
        H[s] = sec[s].front().p0();
    }
    std::unique_ptr<gsBSpline<T>> hN(new gsBSpline<T>(*hJ[2])); hN->reverse(); hN->knots().affineTransformTo(T(0), T(1));
    std::unique_ptr<gsBSpline<T>> hW(new gsBSpline<T>(*hJ[3])); hW->reverse(); hW->knots().affineTransformTo(T(0), T(1));
    std::unique_ptr<gsBSpline<T>> sp0(makeStraightBSpline(P0, H[0], degree));
    std::unique_ptr<gsBSpline<T>> sp1i(makeStraightBSpline(P1, H[1], degree));
    std::unique_ptr<gsBSpline<T>> sp1o(makeStraightBSpline(H[1], P1, degree));
    std::unique_ptr<gsBSpline<T>> sp2o(makeStraightBSpline(H[2], P2, degree));
    std::unique_ptr<gsBSpline<T>> sp3o(makeStraightBSpline(H[3], P3, degree));
    std::unique_ptr<gsBSpline<T>> sp3i(makeStraightBSpline(P3, H[3], degree));
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*cS.bsp, *sp1i, *hJ[0], *sp0, degree)));
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*sp1o, *cE.bsp, *sp2o, *hJ[1], degree)));
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*hN, *sp2o, *cN.bsp, *sp3o, degree)));
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*sp0, *hW, *sp3i, *cW.bsp, degree)));
}

// Conforming quad ring between a 4-curve outer loop and an arbitrary inner hole
template <class T>
inline void buildQuadRingByProximity(
    const gsCurveLoop<T>& outer,
    const gsCurveLoop<T>& hole,
    gsMultiPatch<T>& mp,
    int degree)
{
    if (outer.size() == 4 && hole.size() == 4) {
        // Check if 1-to-1 curve proximity gives valid single-sign Coons ring (ideal for curved 4-arc annulus/plate)
        gsMultiPatch<T> fastMp;
        for (int i = 0; i < 4; ++i) {
            const auto& oCrv = outer.curve(i);
            auto pA = startPoint(oCrv);
            auto pB = endPoint(oCrv);
            int best_h = -1;
            bool best_rev = false;
            T best_dist = T(1e18);
            for (int h = 0; h < 4; ++h) {
                const auto& hCrv = hole.curve(h);
                auto hA = startPoint(hCrv);
                auto hB = endPoint(hCrv);
                T d_fwd = (pA.x - hA.x)*(pA.x - hA.x) + (pA.y - hA.y)*(pA.y - hA.y) +
                          (pB.x - hB.x)*(pB.x - hB.x) + (pB.y - hB.y)*(pB.y - hB.y);
                if (d_fwd < best_dist) { best_dist = d_fwd; best_h = h; best_rev = false; }
                T d_rev = (pA.x - hB.x)*(pA.x - hB.x) + (pA.y - hB.y)*(pA.y - hB.y) +
                          (pB.x - hA.x)*(pB.x - hA.x) + (pB.y - hA.y)*(pB.y - hA.y);
                if (d_rev < best_dist) { best_dist = d_rev; best_h = h; best_rev = true; }
            }
            std::unique_ptr<gsCurve<T>> hCrv(hole.curve(best_h).clone());
            if (best_rev) hCrv->reverse();
            auto hA = startPoint(*hCrv);
            auto hB = endPoint(*hCrv);
            std::unique_ptr<gsBSpline<T>> sp_S(makeStraightBSpline(hA, pA, degree));
            std::unique_ptr<gsBSpline<T>> sp_N(makeStraightBSpline(hB, pB, degree));
            auto patch = makeCoonsPatchFrom4Curves(*sp_S, oCrv, *sp_N, *hCrv, degree);
            fastMp.addPatch(typename gsGeometry<T>::uPtr(patch));
        }
        bool ok = (fastMp.nPatches() == 4);
        if (ok) {
            fastMp.computeTopology(T(1e-5), false, false);
            if (fastMp.topology().nInterfaces() != 4) ok = false;
        }
        for (size_t p = 0; ok && p < fastMp.nPatches(); ++p) {
            if (!isPatchDetSingleSign(fastMp.patch(p), 15, T(1e-8))) ok = false;
        }
        if (ok) {
            for (size_t p = 0; p < fastMp.nPatches(); ++p)
                mp.addPatch(typename gsGeometry<T>::uPtr(fastMp.patch(p).clone().release()));
            return;
        }
    }

    if (outer.size() == 4) {
        std::vector<CurveSeg<T>> oCCW = orientLoopSegments(outer, degree, true);
        buildBlockQuadRing(*oCCW[0].bsp, *oCCW[1].bsp, *oCCW[2].reversed().bsp, *oCCW[3].reversed().bsp, hole, mp, degree);
    }
}

// Create an exact or cubic-interpolated CurveSeg connecting points at iso-parameter on patch.
template <class T>
inline CurveSeg<T> extractBlockCutCurve(
    const gsGeometry<T>& refPatch,
    T u0, T v0, T u1, T v1, int degree)
{
    gsMatrix<T> uv(2, 4), pts;
    for (int k = 0; k < 4; ++k) {
        T t = static_cast<T>(k) / T(3);
        uv(0, k) = (T(1) - t) * u0 + t * u1;
        uv(1, k) = (T(1) - t) * v0 + t * v1;
    }
    refPatch.eval_into(uv, pts);
    gsPoint2D<T> q0(pts(0, 0), pts(1, 0)), q1(pts(0, 1), pts(1, 1));
    gsPoint2D<T> q2(pts(0, 2), pts(1, 2)), q3(pts(0, 3), pts(1, 3));

    // Check if essentially straight
    auto distToChord = [&](const gsPoint2D<T>& q) -> T {
        T dx = q3.x - q0.x, dy = q3.y - q0.y;
        T L = std::hypot(dx, dy);
        if (L <= T(1e-14)) return T(0);
        return std::abs(dx * (q.y - q0.y) - dy * (q.x - q0.x)) / L;
    };
    if (distToChord(q1) < T(1e-5) && distToChord(q2) < T(1e-5)) {
        std::unique_ptr<gsBSpline<T>> bsp(makeStraightBSpline(q0, q3, std::max(3, degree)));
        return CurveSeg<T>(*bsp);
    }

    // Closed-form cubic Bezier interpolating q0, q1(1/3), q2(2/3), q3(1)
    gsPoint2D<T> p1((-T(5) * q0.x + T(18) * q1.x - T(9) * q2.x + T(2) * q3.x) / T(6),
                    (-T(5) * q0.y + T(18) * q1.y - T(9) * q2.y + T(2) * q3.y) / T(6));
    gsPoint2D<T> p2((T(2) * q0.x - T(9) * q1.x + T(18) * q2.x - T(5) * q3.x) / T(6),
                    (T(2) * q0.y - T(9) * q1.y + T(18) * q2.y - T(5) * q3.y) / T(6));
    gsKnotVector<T> kv(T(0), T(1), 0, 4, 1, 3);
    gsMatrix<T> coefs(4, 2);
    coefs << q0.x, q0.y, p1.x, p1.y, p2.x, p2.y, q3.x, q3.y;
    gsBSpline<T> bsp(kv, coefs);
    if (bsp.degree() < degree) bsp.degreeElevate(degree - bsp.degree());
    return CurveSeg<T>(bsp);
}

// Decompose an arbitrary 4-sided macro-block [cS(P0->P1), cE(P1->P2), cN(P3->P2), cW(P0->P3)]
// enclosing 0, 1, or N holes into a conforming multipatch without T-junctions.
template <class T>
inline bool build4SidedMultiHoleGrid(
    const CurveSeg<T>& cS, const CurveSeg<T>& cE,
    const CurveSeg<T>& cN, const CurveSeg<T>& cW,
    const std::vector<const gsCurveLoop<T>*>& holes,
    gsMultiPatch<T>& mp,
    int degree)
{
    auto dist2 = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
        T dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    };
    CurveSeg<T> oS = cS, oE = cE, oN = cN, oW = cW;
    gsPoint2D<T> P0 = oS.p0(), P1 = oS.p1();
    if (dist2(oE.p0(), P1) > dist2(oE.p1(), P1)) oE = oE.reversed();
    if (dist2(oW.p0(), P0) > dist2(oW.p1(), P0)) oW = oW.reversed();
    gsPoint2D<T> P3 = oW.p1();
    if (dist2(oN.p0(), P3) > dist2(oN.p1(), P3)) oN = oN.reversed();

    if (holes.empty()) {
        mp.addPatch(typename gsGeometry<T>::uPtr(
            makeCoonsPatchFrom4Curves(*oS.bsp, *oE.bsp, *oN.bsp, *oW.bsp, degree)));
        return true;
    }
    if (holes.size() == 1) {
        buildBlockQuadRing(*oS.bsp, *oE.bsp, *oN.bsp, *oW.bsp, *holes[0], mp, degree);
        return true;
    }

    // Build smooth 4-sided reference Coons map over [0,1]^2
    std::unique_ptr<gsTensorBSpline<2, T>> refPatch(
        makeCoonsPatchFrom4Curves(*oS.bsp, *oE.bsp, *oN.bsp, *oW.bsp, degree));

    // Sample refPatch on 41x41 grid to invert physical points to (u, v) in [0,1]^2
    const int G = 41;
    std::vector<gsPoint2D<T>> gPts(G * G);
    {
        gsMatrix<T> uv(2, G * G), val;
        for (int j = 0; j < G; ++j)
            for (int i = 0; i < G; ++i) {
                uv(0, j * G + i) = static_cast<T>(i) / (G - 1);
                uv(1, j * G + i) = static_cast<T>(j) / (G - 1);
            }
        refPatch->eval_into(uv, val);
        for (int idx = 0; idx < G * G; ++idx)
            gPts[idx] = gsPoint2D<T>(val(0, idx), val(1, idx));
    }

    auto invertPointUV = [&](const gsPoint2D<T>& q) -> gsPoint2D<T> {
        T bestD = T(1e18); int bi = 0, bj = 0;
        for (int j = 0; j < G; ++j)
            for (int i = 0; i < G; ++i) {
                T d = dist2(q, gPts[j * G + i]);
                if (d < bestD) { bestD = d; bi = i; bj = j; }
            }
        T u = static_cast<T>(bi) / (G - 1), v = static_cast<T>(bj) / (G - 1);
        gsMatrix<T> uv(2, 1), pt, jac;
        for (int it = 0; it < 6; ++it) {
            uv(0, 0) = u; uv(1, 0) = v;
            refPatch->eval_into(uv, pt);
            refPatch->jacobian_into(uv, jac);
            T rx = q.x - pt(0, 0), ry = q.y - pt(1, 0);
            T det = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
            if (std::abs(det) < T(1e-14)) break;
            T du = ( jac(1, 1) * rx - jac(0, 1) * ry) / det;
            T dv = (-jac(1, 0) * rx + jac(0, 0) * ry) / det;
            u = std::max(T(0), std::min(T(1), u + du));
            v = std::max(T(0), std::min(T(1), v + dv));
            if (du * du + dv * dv < T(1e-16)) break;
        }
        return gsPoint2D<T>(u, v);
    };

    const int nH = static_cast<int>(holes.size());
    struct HoleUVBox {
        T u0 = 1, u1 = 0, v0 = 1, v1 = 0;
        T uc = 0.5, vc = 0.5;
        std::vector<gsPoint2D<T>> poly;
    };
    std::vector<HoleUVBox> hBoxes(nH);
    for (int h = 0; h < nH; ++h) {
        const auto& hl = *holes[h];
        for (int c = 0; c < hl.size(); ++c) {
            gsMatrix<T> pr = hl.curve(c).parameterRange();
            for (int k = 0; k < 8; ++k) {
                T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(k) / T(8));
                gsPoint2D<T> p = evalCurvePoint(hl.curve(c), t);
                hBoxes[h].poly.push_back(p);
                gsPoint2D<T> uv = invertPointUV(p);
                hBoxes[h].u0 = std::min(hBoxes[h].u0, uv.x);
                hBoxes[h].u1 = std::max(hBoxes[h].u1, uv.x);
                hBoxes[h].v0 = std::min(hBoxes[h].v0, uv.y);
                hBoxes[h].v1 = std::max(hBoxes[h].v1, uv.y);
            }
        }
        hBoxes[h].uc = T(0.5) * (hBoxes[h].u0 + hBoxes[h].u1);
        hBoxes[h].vc = T(0.5) * (hBoxes[h].v0 + hBoxes[h].v1);
    }

    // 1. Find all gap-midpoint candidate separator cuts in u and v that do not pierce any hole
    std::vector<std::pair<T, T>> uIvals(nH), vIvals(nH);
    std::vector<gsPoint2D<T>> uvCenters(nH);
    for (int h = 0; h < nH; ++h) {
        uIvals[h] = {hBoxes[h].u0, hBoxes[h].u1};
        vIvals[h] = {hBoxes[h].v0, hBoxes[h].v1};
        uvCenters[h] = {hBoxes[h].uc, hBoxes[h].vc};
    }
    std::vector<T> uCuts = SeparatorGrid1D<T>::findClearCuts(uIvals, T(0.01));
    std::vector<T> vCuts = SeparatorGrid1D<T>::findClearCuts(vIvals, T(0.01));

    SeparatorGrid1D<T>::pruneRedundantCuts(uCuts, vCuts, uIvals, vIvals, uvCenters);

    if (SeparatorGrid1D<T>::gridMaxItems(uCuts, vCuts, uvCenters) <= 1 && (!uCuts.empty() || !vCuts.empty())) {
        const int nx = static_cast<int>(uCuts.size()) + 1;
        const int ny = static_cast<int>(vCuts.size()) + 1;
        std::vector<T> us = {T(0)};
        us.insert(us.end(), uCuts.begin(), uCuts.end());
        us.push_back(T(1));
        std::vector<T> vs = {T(0)};
        vs.insert(vs.end(), vCuts.begin(), vCuts.end());
        vs.push_back(T(1));

        std::vector<CurveSeg<T>> sSegs = splitSegMulti(oS, uCuts);
        std::vector<CurveSeg<T>> nSegs = splitSegMulti(oN, uCuts);
        std::vector<CurveSeg<T>> wSegs = splitSegMulti(oW, vCuts);
        std::vector<CurveSeg<T>> eSegs = splitSegMulti(oE, vCuts);

        // Precompute exact interior vertical cuts (at u = us[i] for i=1..nx-1) and horizontal cuts
        std::vector<std::vector<CurveSeg<T>>> vEdges(nx + 1, std::vector<CurveSeg<T>>(ny));
        std::vector<std::vector<CurveSeg<T>>> hEdges(ny + 1, std::vector<CurveSeg<T>>(nx));
        for (int j = 0; j < ny; ++j) {
            vEdges[0][j]  = wSegs[j];
            vEdges[nx][j] = eSegs[j];
        }
        for (int i = 0; i < nx; ++i) {
            hEdges[0][i]  = sSegs[i];
            hEdges[ny][i] = nSegs[i];
        }
        for (int i = 1; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                vEdges[i][j] = extractBlockCutCurve(*refPatch, us[i], vs[j], us[i], vs[j + 1], degree);
        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                hEdges[j][i] = extractBlockCutCurve(*refPatch, us[i], vs[j], us[i + 1], vs[j], degree);

        // Assign holes to grid cells
        std::vector<const gsCurveLoop<T>*> cellHole(nx * ny, nullptr);
        for (int h = 0; h < nH; ++h) {
            int ix = static_cast<int>(std::upper_bound(uCuts.begin(), uCuts.end(), hBoxes[h].uc) - uCuts.begin());
            int iy = static_cast<int>(std::upper_bound(vCuts.begin(), vCuts.end(), hBoxes[h].vc) - vCuts.begin());
            cellHole[iy * nx + ix] = holes[h];
        }

        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const auto& cBot = hEdges[iy][ix];
                const auto& cRgt = vEdges[ix + 1][iy];
                const auto& cTop = hEdges[iy + 1][ix];
                const auto& cLft = vEdges[ix][iy];
                if (cellHole[iy * nx + ix]) {
                    buildBlockQuadRing(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, *cellHole[iy * nx + ix], mp, degree);
                } else {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsPatchFrom4Curves(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, degree)));
                }
            }
        }
        return true;
    }

    // 2. For 2 holes with overlapping axis bboxes, try tilted straight cut in (u, v) space
    if (nH == 2) {
        const T cands[5] = {T(0.25), T(0.38), T(0.50), T(0.62), T(0.75)};
        for (int axis = 0; axis < 2; ++axis) {
            for (T a0 : cands) {
                for (T a1 : cands) {
                    auto ptAt = [&](T u, T v) -> gsPoint2D<T> {
                        gsMatrix<T> uv(2, 1), out; uv(0, 0) = u; uv(1, 0) = v;
                        refPatch->eval_into(uv, out);
                        return gsPoint2D<T>(out(0, 0), out(1, 0));
                    };
                    gsPoint2D<T> A = (axis == 0) ? ptAt(a0, T(0)) : ptAt(T(0), a0);
                    gsPoint2D<T> B = (axis == 0) ? ptAt(a1, T(1)) : ptAt(T(1), a1);
                    auto sideSign = [&](const std::vector<gsPoint2D<T>>& poly) -> int {
                        int sgn = 0;
                        for (const auto& p : poly) {
                            T cr = (B.x - A.x) * (p.y - A.y) - (B.y - A.y) * (p.x - A.x);
                            int s = (cr > T(1e-3)) ? 1 : ((cr < -T(1e-3)) ? -1 : 0);
                            if (s == 0) return 0;
                            if (sgn == 0) sgn = s;
                            else if (s != sgn) return 0;
                        }
                        return sgn;
                    };
                    int s0 = sideSign(hBoxes[0].poly), s1 = sideSign(hBoxes[1].poly);
                    if (s0 != 0 && s1 != 0 && s0 != s1) {
                        if (axis == 0) {
                            auto sSp = oS.splitAtParam(a0);
                            auto nSp = oN.splitAtParam(a1);
                            CurveSeg<T> cut = extractBlockCutCurve(*refPatch, a0, T(0), a1, T(1), degree);
                            const gsCurveLoop<T>* hLeft  = (s0 > 0) ? holes[0] : holes[1];
                            const gsCurveLoop<T>* hRight = (s0 > 0) ? holes[1] : holes[0];
                            buildBlockQuadRing(*sSp.first.bsp, *cut.bsp, *nSp.first.bsp, *oW.bsp, *hLeft, mp, degree);
                            buildBlockQuadRing(*sSp.second.bsp, *oE.bsp, *nSp.second.bsp, *cut.bsp, *hRight, mp, degree);
                        } else {
                            auto wSp = oW.splitAtParam(a0);
                            auto eSp = oE.splitAtParam(a1);
                            CurveSeg<T> cut = extractBlockCutCurve(*refPatch, T(0), a0, T(1), a1, degree);
                            const gsCurveLoop<T>* hBot = (s0 < 0) ? holes[0] : holes[1];
                            const gsCurveLoop<T>* hTop = (s0 < 0) ? holes[1] : holes[0];
                            buildBlockQuadRing(*oS.bsp, *eSp.first.bsp, *cut.bsp, *wSp.first.bsp, *hBot, mp, degree);
                            buildBlockQuadRing(*cut.bsp, *eSp.second.bsp, *oN.bsp, *wSp.second.bsp, *hTop, mp, degree);
                        }
                        return true;
                    }
                }
            }
        }
    }

    // 3. For 3 holes in a non-grid triangle (e.g. plate_4corner_3holes), try 3-way Steiner partition into 3 quads
    if (nH == 3) {
        // Rotate 4 block sides so P[k] is central-split corner
        CurveSeg<T> edgesCCW[4] = { oS, oE, oN.reversed(), oW.reversed() };
        gsPoint2D<T> ctr[3];
        for (int h = 0; h < 3; ++h) {
            ctr[h] = polygonCentroid(hBoxes[h].poly);
        }
        gsPoint2D<T> O_base((ctr[0].x + ctr[1].x + ctr[2].x) / T(3),
                            (ctr[0].y + ctr[1].y + ctr[2].y) / T(3));
        for (int k = 0; k < 4; ++k) {
            const CurveSeg<T>& eK    = edgesCCW[k];       // P_k -> P_{k+1}
            const CurveSeg<T>& eNext = edgesCCW[(k+1)%4]; // P_{k+1} -> P_{k+2}
            const CurveSeg<T>& ePrev = edgesCCW[(k+2)%4]; // P_{k+2} -> P_{k+3}
            const CurveSeg<T>& eLast = edgesCCW[(k+3)%4]; // P_{k+3} -> P_k
            gsPoint2D<T> P_opp = eNext.p1(); // P_{k+2}

            for (T t1 : {T(0.4), T(0.5), T(0.6)}) {
                for (T t2 : {T(0.4), T(0.5), T(0.6)}) {
                    auto spA = eK.splitAtParam(t1);     // P_k -> M1 -> P_{k+1}
                    auto spB = eLast.splitAtParam(t2);  // P_{k+3} -> M2 -> P_k
                    gsPoint2D<T> M1 = spA.first.p1(), M2 = spB.first.p1();
                    for (T w : {T(0.0), T(0.15)}) {
                        gsPoint2D<T> O((T(1) - w) * O_base.x + w * P_opp.x,
                                       (T(1) - w) * O_base.y + w * P_opp.y);
                        bool oInHole = false;
                        for (int h = 0; h < 3; ++h)
                            if (pointInPolyStrict(O, hBoxes[h].poly, T(1e-4))) { oInHole = true; break; }
                        if (oInHole) continue;

                        std::vector<gsPoint2D<T>> qPoly[3] = {
                            {spA.first.p0(), M1, O, M2},
                            {M1, eK.p1(), P_opp, O},
                            {O, P_opp, eLast.p0(), M2}
                        };
                        gsPoint2D<T> qq[3][4];
                        bool convOk = true;
                        for (int q = 0; q < 3; ++q) {
                            for (int v = 0; v < 4; ++v) qq[q][v] = qPoly[q][v];
                            if (!isStrictlyConvexQuad(qq[q])) { convOk = false; break; }
                        }
                        if (!convOk) continue;

                        int match[3] = {-1, -1, -1};
                        bool allClear = true;
                        for (int q = 0; q < 3; ++q) {
                            for (int h = 0; h < 3; ++h) {
                                bool allIn = true;
                                for (const auto& hp : hBoxes[h].poly)
                                    if (!pointInPolyStrict(hp, qPoly[q], T(1e-4))) { allIn = false; break; }
                                if (allIn) {
                                    if (match[q] != -1) allClear = false;
                                    match[q] = h;
                                }
                            }
                            if (match[q] == -1) { allClear = false; break; }
                        }
                        if (!allClear || match[0] == match[1] || match[0] == match[2] || match[1] == match[2])
                            continue;

                        std::unique_ptr<gsBSpline<T>> s1(makeStraightBSpline(M1, O, degree));
                        std::unique_ptr<gsBSpline<T>> s2(makeStraightBSpline(O, P_opp, degree));
                        std::unique_ptr<gsBSpline<T>> s3(makeStraightBSpline(M2, O, degree));
                        std::unique_ptr<gsBSpline<T>> s3_rev(s3->clone()); s3_rev->reverse();
                        std::unique_ptr<gsBSpline<T>> s2_rev(s2->clone()); s2_rev->reverse();

                        buildBlockQuadRing(*spA.first.bsp, *s1, *s3_rev, *spB.second.reversed().bsp, *holes[match[0]], mp, degree);
                        buildBlockQuadRing(*spA.second.bsp, *eNext.bsp, *s2_rev, *s1, *holes[match[1]], mp, degree);
                        buildBlockQuadRing(*s2, *ePrev.bsp, *spB.first.bsp, *s3, *holes[match[2]], mp, degree);
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

// Decompose a 4-sided outer domain [oCCW] enclosing >= 2 holes that separate cleanly along
// physical XY axes inside an axis-aligned inner box Q, joined to oCCW by a conforming O-grid frame.
template <class T>
inline bool build4SidedInnerAxisGrid(
    const std::vector<CurveSeg<T>>& oCCW,
    const std::vector<const gsCurveLoop<T>*>& holes,
    gsMultiPatch<T>& mp,
    int degree)
{
    if (oCCW.size() != 4 || holes.size() < 2) return false;
    const int nH = static_cast<int>(holes.size());

    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        return makeStraightSeg(a, b, degree);
    };

    std::vector<gsPoint2D<T>> outerPoly;
    outerPoly.reserve(32);
    T oxMin = T(1e18), oxMax = -T(1e18), oyMin = T(1e18), oyMax = -T(1e18);
    for (int s = 0; s < 4; ++s) {
        for (int k = 0; k < 8; ++k) {
            gsPoint2D<T> q = evalCurvePoint(*oCCW[s].bsp, static_cast<T>(k) / T(8));
            outerPoly.push_back(q);
            oxMin = std::min(oxMin, q.x); oxMax = std::max(oxMax, q.x);
            oyMin = std::min(oyMin, q.y); oyMax = std::max(oyMax, q.y);
        }
    }
    const T outerDiag = std::hypot(oxMax - oxMin, oyMax - oyMin);
    if (outerDiag <= T(1e-12)) return false;

    struct HoleXYBox { T x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18, xc = 0, yc = 0; };
    std::vector<HoleXYBox> hBoxes(nH);
    T hxMin = T(1e18), hxMax = -T(1e18), hyMin = T(1e18), hyMax = -T(1e18);
    for (int h = 0; h < nH; ++h) {
        const auto& hl = *holes[h];
        for (int c = 0; c < hl.size(); ++c) {
            gsMatrix<T> pr = hl.curve(c).parameterRange();
            for (int k = 0; k < 8; ++k) {
                T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(k) / T(8));
                gsPoint2D<T> p = evalCurvePoint(hl.curve(c), t);
                hBoxes[h].x0 = std::min(hBoxes[h].x0, p.x); hBoxes[h].x1 = std::max(hBoxes[h].x1, p.x);
                hBoxes[h].y0 = std::min(hBoxes[h].y0, p.y); hBoxes[h].y1 = std::max(hBoxes[h].y1, p.y);
            }
        }
        hBoxes[h].xc = T(0.5) * (hBoxes[h].x0 + hBoxes[h].x1);
        hBoxes[h].yc = T(0.5) * (hBoxes[h].y0 + hBoxes[h].y1);
        hxMin = std::min(hxMin, hBoxes[h].x0); hxMax = std::max(hxMax, hBoxes[h].x1);
        hyMin = std::min(hyMin, hBoxes[h].y0); hyMax = std::max(hyMax, hBoxes[h].y1);
    }
    const T hSpan = std::max(hxMax - hxMin, hyMax - hyMin);
    if (hSpan <= T(1e-12)) return false;

    std::vector<std::pair<T, T>> xIvals(nH), yIvals(nH);
    std::vector<gsPoint2D<T>> xyCenters(nH);
    T loAllX = T(1e18), hiAllX = -T(1e18);
    T loAllY = T(1e18), hiAllY = -T(1e18);
    for (int h = 0; h < nH; ++h) {
        xIvals[h] = {hBoxes[h].x0, hBoxes[h].x1};
        yIvals[h] = {hBoxes[h].y0, hBoxes[h].y1};
        xyCenters[h] = {hBoxes[h].xc, hBoxes[h].yc};
        loAllX = std::min(loAllX, hBoxes[h].x0); hiAllX = std::max(hiAllX, hBoxes[h].x1);
        loAllY = std::min(loAllY, hBoxes[h].y0); hiAllY = std::max(hiAllY, hBoxes[h].y1);
    }
    T minGapX = T(0.02) * std::max(hiAllX - loAllX, T(1e-6));
    T minGapY = T(0.02) * std::max(hiAllY - loAllY, T(1e-6));
    std::vector<T> xCuts = SeparatorGrid1D<T>::findClearCuts(xIvals, minGapX);
    std::vector<T> yCuts = SeparatorGrid1D<T>::findClearCuts(yIvals, minGapY);

    SeparatorGrid1D<T>::pruneRedundantCuts(xCuts, yCuts, xIvals, yIvals, xyCenters);

    if (SeparatorGrid1D<T>::gridMaxItems(xCuts, yCuts, xyCenters) > 1 || (xCuts.empty() && yCuts.empty()))
        return false;

    const int nx = static_cast<int>(xCuts.size()) + 1;
    const int ny = static_cast<int>(yCuts.size()) + 1;
    const int P_in = 2 * (nx + ny);

    auto splitOptsForCount = [](int n) -> std::vector<std::vector<T>> {
        if (n <= 1) return {{}};
        if (n == 2) return {{T(0.50)}, {T(0.35)}, {T(0.65)}};
        if (n == 3) return {
            {T(1) / T(3), T(2) / T(3)},
            {T(0.18), T(0.50)},
            {T(0.50), T(0.82)},
            {T(0.25), T(0.75)}
        };
        std::vector<T> uni(n - 1);
        for (int k = 1; k < n; ++k) uni[k - 1] = static_cast<T>(k) / static_cast<T>(n);
        return {uni};
    };

    auto ptSegDist = [](const gsPoint2D<T>& p, const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> T {
        T dx = b.x - a.x, dy = b.y - a.y;
        T L2 = dx * dx + dy * dy;
        T t = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2)) : T(0);
        return std::hypot(a.x + t * dx - p.x, a.y + t * dy - p.y);
    };
    auto minDistToOuter = [&](const gsPoint2D<T>& p) -> T {
        T mn = T(1e18);
        const int no = static_cast<int>(outerPoly.size());
        for (int i = 0; i < no; ++i)
            mn = std::min(mn, ptSegDist(p, outerPoly[i], outerPoly[(i + 1) % no]));
        return mn;
    };

    bool foundFrame = false;
    T bestFrameScore = T(1e18);
    std::vector<std::unique_ptr<gsGeometry<T>>> bestFrame;
    std::vector<T> bestXs, bestYs;

    for (T padFrac : {T(0.075), T(0.095), T(0.055), T(0.12)}) {
        T pad = padFrac * hSpan;
        T qx0 = hxMin - pad, qx1 = hxMax + pad;
        T qy0 = hyMin - pad, qy1 = hyMax + pad;
        gsPoint2D<T> qc[4] = {{qx0, qy0}, {qx1, qy0}, {qx1, qy1}, {qx0, qy1}};
        bool boxClear = true;
        for (int c = 0; c < 4; ++c) {
            if (!pointInPolyStrict(qc[c], outerPoly, T(1e-3)) ||
                minDistToOuter(qc[c]) < T(0.025) * outerDiag)
            {
                boxClear = false; break;
            }
            gsPoint2D<T> mid(T(0.5) * (qc[c].x + qc[(c + 1) % 4].x),
                             T(0.5) * (qc[c].y + qc[(c + 1) % 4].y));
            if (!pointInPolyStrict(mid, outerPoly, T(1e-3)) ||
                minDistToOuter(mid) < T(0.025) * outerDiag)
            {
                boxClear = false; break;
            }
        }
        if (!boxClear) continue;

        std::vector<T> xs = {qx0};
        xs.insert(xs.end(), xCuts.begin(), xCuts.end());
        xs.push_back(qx1);
        std::vector<T> ys = {qy0};
        ys.insert(ys.end(), yCuts.begin(), yCuts.end());
        ys.push_back(qy1);

        std::vector<CurveSeg<T>> iSeg;
        iSeg.reserve(P_in);
        for (int i = 0; i < nx; ++i)      iSeg.push_back(mkLine({xs[i], qy0},     {xs[i + 1], qy0}));
        for (int j = 0; j < ny; ++j)      iSeg.push_back(mkLine({qx1,   ys[j]},   {qx1, ys[j + 1]}));
        for (int i = nx; i > 0; --i)      iSeg.push_back(mkLine({xs[i], qy1},     {xs[i - 1], qy1}));
        for (int j = ny; j > 0; --j)      iSeg.push_back(mkLine({qx0,   ys[j]},   {qx0, ys[j - 1]}));

        std::vector<gsPoint2D<T>> qPoly = {qc[0], qc[1], qc[2], qc[3]};

        std::vector<std::array<int, 4>> sideModes = {{nx, ny, nx, ny}};
        if (nx != ny) sideModes.push_back({ny, nx, ny, nx});

        for (const auto& nSide : sideModes) {
            auto opts0 = splitOptsForCount(nSide[0]);
            auto opts1 = splitOptsForCount(nSide[1]);
            for (const auto& ts0 : opts0) {
                auto s0 = splitSegMulti(oCCW[0], ts0);
                auto s2 = splitSegMulti(oCCW[2], ts0);
                for (const auto& ts1 : opts1) {
                    auto s1 = splitSegMulti(oCCW[1], ts1);
                    auto s3 = splitSegMulti(oCCW[3], ts1);
                    std::vector<CurveSeg<T>> oSeg;
                    oSeg.reserve(P_in);
                    oSeg.insert(oSeg.end(), s0.begin(), s0.end());
                    oSeg.insert(oSeg.end(), s1.begin(), s1.end());
                    oSeg.insert(oSeg.end(), s2.begin(), s2.end());
                    oSeg.insert(oSeg.end(), s3.begin(), s3.end());
                    if (static_cast<int>(oSeg.size()) != P_in) continue;

                    for (int shift = 0; shift < P_in; ++shift) {
                        std::vector<std::unique_ptr<gsGeometry<T>>> frameCand;
                        frameCand.reserve(P_in);
                        for (int k = 0; k < P_in; ++k) {
                            const auto& oE = oSeg[k];
                            const auto& iE = iSeg[(k + shift) % P_in];
                            CurveSeg<T> spL = mkLine(oE.p0(), iE.p0());
                            CurveSeg<T> spR = mkLine(oE.p1(), iE.p1());
                            frameCand.emplace_back(makeCoonsPatchFrom4Curves(*oE.bsp, *spR.bsp, *iE.bsp, *spL.bsp, degree));
                        }
                        T sc = 0, md = 0, ma = 0;
                        if (evaluateRingPatches(frameCand, qPoly, sc, md, ma)) {
                            if (!foundFrame || sc < bestFrameScore) {
                                foundFrame = true;
                                bestFrameScore = sc;
                                bestFrame = std::move(frameCand);
                                bestXs = xs;
                                bestYs = ys;
                            }
                        }
                    }
                }
            }
        }
        if (foundFrame) break;
    }

    if (!foundFrame) return false;

    for (auto& p : bestFrame)
        mp.addPatch(typename gsGeometry<T>::uPtr(p.release()));

    std::vector<const gsCurveLoop<T>*> cellHole(nx * ny, nullptr);
    for (int h = 0; h < nH; ++h) {
        int ix = static_cast<int>(std::upper_bound(xCuts.begin(), xCuts.end(), hBoxes[h].xc) - xCuts.begin());
        int iy = static_cast<int>(std::upper_bound(yCuts.begin(), yCuts.end(), hBoxes[h].yc) - yCuts.begin());
        cellHole[iy * nx + ix] = holes[h];
    }

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            CurveSeg<T> cBot = mkLine({bestXs[ix],     bestYs[iy]},     {bestXs[ix + 1], bestYs[iy]});
            CurveSeg<T> cRgt = mkLine({bestXs[ix + 1], bestYs[iy]},     {bestXs[ix + 1], bestYs[iy + 1]});
            CurveSeg<T> cTop = mkLine({bestXs[ix],     bestYs[iy + 1]}, {bestXs[ix + 1], bestYs[iy + 1]});
            CurveSeg<T> cLft = mkLine({bestXs[ix],     bestYs[iy]},     {bestXs[ix],     bestYs[iy + 1]});
            if (cellHole[iy * nx + ix]) {
                buildBlockQuadRing(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, *cellHole[iy * nx + ix], mp, degree);
            } else {
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsPatchFrom4Curves(*cBot.bsp, *cRgt.bsp, *cTop.bsp, *cLft.bsp, degree)));
            }
        }
    }
    return true;
}

} // namespace internal_splitter
} // namespace gismo

#include <gsModeling/gsPatchSplitterPolychord.hpp>
#include <gsModeling/gsPatchSplitterCrossField.hpp>

namespace gismo {

template <class T>
T gsPatchSplitter<T>::getDomainMinCornerAngleDeg() const
{
    if (m_cachedDomainMinCornerAngle < 0 && m_hasDomain) {
        m_cachedDomainMinCornerAngle = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
    }
    return (m_cachedDomainMinCornerAngle >= 0) ? m_cachedDomainMinCornerAngle : T(90);
}

// ============================================================================
// Approach 1: Boundary-Aligned 4-RoSy Harmonic Cross-Field Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeCrossFieldSplitting()
{
    using namespace internal_splitter;
    SplittingCandidate bestCand;
    bestCand.methodUsed = "HarmonicCrossField";
    bestCand.isValid = false;

    if (m_domain.numLoops() == 0) return bestCand;

    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    // 1. Solve 4-RoSy harmonic cross-field + separatrix tracing + strip collapse
    {
        gsMultiPatch<T> cfMp;
        std::string label;
        const T stripRatio = static_cast<T>(m_options.getReal("StripCollapseRatio"));
        const T angThresh  = static_cast<T>(m_options.getReal("AngleThresholdDeg"));
        if (buildCrossFieldArrangementCandidate(m_domain, m_degree, stripRatio, angThresh, cfMp, label)) {
            optimizeWinslowBarrier(cfMp);
            cfMp.computeTopology(T(1e-5), false, false);
            if (validateMultiPatch(m_domain, cfMp, T(1e-4), nullptr, minAngTol, maxAngTol)) {
                T minDet = 0, minA = 0, maxA = 180;
                T energy = evaluateDistortionEnergy(cfMp, &minDet, &minA, &maxA);
                bestCand.multiPatch = give(cfMp);
                bestCand.numPatches = static_cast<int>(bestCand.multiPatch.nPatches());
                bestCand.minJacobianDet = minDet;
                bestCand.distortionEnergy = energy;
                bestCand.minCornerAngleDeg = minA;
                bestCand.maxCornerAngleDeg = maxA;
                bestCand.methodUsed = label;
                bestCand.isValid = true;
            }
        }
    }

    // 2. Compare against conforming triangulation-joining decomposition (hole-free or multi-hole)
    SplittingCandidate triCand = (m_domain.numHoles() == 0)
        ? computeHoleFreeTriangulationSplitting()
        : computeMultiHoleTriangulationSplitting();

    if (triCand.isValid) {
        if (!bestCand.isValid) {
            bestCand = std::move(triCand);
        } else {
            T qBest = m_weightPatches * bestCand.numPatches + m_weightDistortion * bestCand.distortionEnergy;
            T qTri  = m_weightPatches * triCand.numPatches + m_weightDistortion * triCand.distortionEnergy;
            if (qTri < qBest) bestCand = std::move(triCand);
        }
    }

    if (bestCand.isValid)
        bestCand.methodUsed = "HarmonicCrossField";

    return bestCand;
}

// ============================================================================
// Hole-Free Structural & Triangulation-Quad Decomposition Helpers
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeHoleFreeStructuralSplitting()
{
    using namespace internal_splitter;

    SplittingCandidate bestCand;
    bestCand.methodUsed = "MedialAxisSkeleton";
    bestCand.isValid = false;

    const auto& outer = m_domain.outer();
    const int K = outer.size();
    if (K == 0) return bestCand;

    std::vector<CurveSeg<T>> origSegs;
    origSegs.reserve(K);
    for (int i = 0; i < K; ++i)
        origSegs.emplace_back(outer.curve(i), m_degree);

    std::vector<gsPoint2D<T>> outerPoly = sampleCurveLoopPoly(outer, 24);
    gsPoint2D<T> centroid = polygonCentroid(outerPoly);

    // Analyze outer vertices & sharp/reflex angles
    std::vector<gsPoint2D<T>> verts(K);
    std::vector<bool> isSharp(K, false), isReflex(K, false);
    int nSharp = 0, nReflex = 0;
    for (int i = 0; i < K; ++i) {
        verts[i] = origSegs[i].p0();
        const auto& prevBsp = *origSegs[(i + K - 1) % K].bsp;
        const auto& nextBsp = *origSegs[i].bsp;
        auto din  = evalCurveDeriv(prevBsp, T(1));
        auto dout = evalCurveDeriv(nextBsp, T(0));
        T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
        if (n1 > 1e-14 && n2 > 1e-14) {
            T cr = (din.x * dout.y - din.y * dout.x) / (n1 * n2);
            T dt = (din.x * dout.x + din.y * dout.y) / (n1 * n2);
            T ang = std::atan2(cr, dt); // CCW: >0 convex turn, <0 reflex turn
            if (std::abs(ang) > T(0.10)) { // > ~5.7 deg
                isSharp[i] = true;
                ++nSharp;
                if (ang < -T(0.10)) {
                    isReflex[i] = true;
                    ++nReflex;
                }
            }
        }
    }

    auto tryFinalize = [&](gsMultiPatch<T>& mp, const std::string& label, SplittingCandidate& out) -> bool {
        if (mp.nPatches() == 0) return false;
        optimizeWinslowBarrier(mp);
        mp.computeTopology(T(1e-5), false, false);
        const T domMinAng = getDomainMinCornerAngleDeg();
        const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
        const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
        if (!validateMultiPatch(m_domain, mp, T(1e-5), nullptr, minAngTol, maxAngTol)) return false;
        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
        out.multiPatch = mp;
        out.numPatches = static_cast<int>(mp.nPatches());
        out.minJacobianDet = minDet;
        out.distortionEnergy = energy;
        out.minCornerAngleDeg = minA;
        out.maxCornerAngleDeg = maxA;
        out.methodUsed = label;
        out.isValid = true;
        return true;
    };

    auto registerIfBetter = [&](const SplittingCandidate& cand) {
        if (!cand.isValid) return;
        if (!bestCand.isValid) { bestCand = cand; return; }
        T qNew = m_weightPatches * cand.numPatches + m_weightDistortion * cand.distortionEnergy;
        T qOld = m_weightPatches * bestCand.numPatches + m_weightDistortion * bestCand.distortionEnergy;
        if (qNew < qOld) bestCand = cand;
    };

    // ------------------------------------------------------------------------
    // 1. Single-patch check: "convexish" 4-curve hole-free quad -> 1 patch
    //    (A quad is convexish iff no corner lies in the convex hull of the other 3,
    //     which is equivalent to isStrictlyConvexQuad(c4); also allow convex domain
    //     corner angles up to ~174 deg with Winslow/refined-CP untangling for S-ducts.)
    // ------------------------------------------------------------------------
    if (K == 4) {
        gsPoint2D<T> c4[4] = {verts[0], verts[1], verts[2], verts[3]};
        if (isStrictlyConvexQuad(c4)) {
            const T radToDeg = T(180) / T(3.14159265358979323846);
            T domMinC = T(180), domMaxC = T(0);
            for (int i = 0; i < 4; ++i) {
                auto dPrev = internal_splitter::evalCurveDeriv(*origSegs[(i + 3) % 4].bsp, T(1));
                auto dCurr = internal_splitter::evalCurveDeriv(*origSegs[i].bsp, T(0));
                T nP = std::hypot(dPrev.x, dPrev.y), nC = std::hypot(dCurr.x, dCurr.y);
                if (nP > T(1e-12) && nC > T(1e-12)) {
                    T ux = -dPrev.x / nP, uy = -dPrev.y / nP;
                    T vx = dCurr.x / nC, vy = dCurr.y / nC;
                    T ang = std::atan2(vx * uy - vy * ux, ux * vx + uy * vy) * radToDeg;
                    if (ang <= T(0)) ang += T(360);
                    domMinC = std::min(domMinC, ang);
                    domMaxC = std::max(domMaxC, ang);
                }
            }
            if (domMinC > T(3) && domMaxC < T(174)) {
                const T localMinAng = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinC);
                const T localMaxAng = std::max(static_cast<T>(m_options.getReal("MaxCornerAngleDeg")), std::min(T(175), domMaxC + T(4)));
                for (int rLevel = 0; rLevel <= 2; ++rLevel) {
                    gsMultiPatch<T> mp;
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(origSegs[0], origSegs[1], origSegs[2], origSegs[3], m_degree)));
                    if (rLevel > 0) {
                        if (auto* tb = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(0)))
                            tb->uniformRefine(rLevel);
                    }
                    optimizeWinslowBarrier(mp);
                    mp.computeTopology(T(1e-5), false, false);
                    if (validateMultiPatch(m_domain, mp, T(1e-5), nullptr, localMinAng, localMaxAng)) {
                        T minDet = 0, minA = 0, maxA = 180;
                        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
                        SplittingCandidate cand;
                        cand.multiPatch = give(mp);
                        cand.numPatches = 1;
                        cand.minJacobianDet = minDet;
                        cand.distortionEnergy = energy;
                        cand.minCornerAngleDeg = minA;
                        cand.maxCornerAngleDeg = maxA;
                        cand.methodUsed = "SingleConvexQuad";
                        cand.isValid = true;
                        return cand;
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 2. O-Grid (5 patches): ideal for 1-corner / 2-corner shapes & 4-station shapes
    // ------------------------------------------------------------------------
    {
        std::vector<CurveSeg<T>> bEdges;
        if (K == 1) {
            auto sHalf = origSegs[0].splitArc(T(0.5));
            auto sQ1   = sHalf.first.splitArc(T(0.5));
            auto sQ2   = sHalf.second.splitArc(T(0.5));
            bEdges = {sQ1.first, sQ1.second, sQ2.first, sQ2.second};
        } else if (K == 2) {
            auto s0 = origSegs[0].splitArc(T(0.5));
            auto s1 = origSegs[1].splitArc(T(0.5));
            bEdges = {s0.first, s0.second, s1.first, s1.second};
        } else if (K == 4) {
            bEdges = origSegs;
        }

        if (bEdges.size() == 4) {
            gsPoint2D<T> P[4];
            for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
            gsPoint2D<T> pAvg((P[0].x + P[1].x + P[2].x + P[3].x) * T(0.25),
                              (P[0].y + P[1].y + P[2].y + P[3].y) * T(0.25));

            std::vector<gsPoint2D<T>> centers = {centroid, pAvg};
            // Also test shifting center away from any reflex station (e.g. heart cleft)
            for (int k = 0; k < 4; ++k) {
                gsPoint2D<T> opp = P[(k + 2) % 4];
                centers.emplace_back(T(0.75) * pAvg.x + T(0.25) * opp.x,
                                     T(0.75) * pAvg.y + T(0.25) * opp.y);
            }

            T alphas[3] = {T(0.45), T(0.35), T(0.55)};
            for (const auto& ctr : centers) {
                if (!pointInPolyStrict(ctr, outerPoly, T(1e-5))) continue;
                for (T alpha : alphas) {
                    gsPoint2D<T> Q[4];
                    bool qInside = true;
                    for (int k = 0; k < 4; ++k) {
                        Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                            ctr.y + alpha * (P[k].y - ctr.y));
                        if (!pointInPolyStrict(Q[k], outerPoly, T(1e-5))) { qInside = false; break; }
                    }
                    if (!qInside || !isStrictlyConvexQuad(Q)) continue;

                    CurveSeg<T> qEdge[4], spoke[4];
                    for (int k = 0; k < 4; ++k) {
                        qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], m_degree);
                        spoke[k] = makeStraightSeg(Q[k], P[k], m_degree);
                    }

                    gsMultiPatch<T> mp;
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], m_degree)));
                    for (int k = 0; k < 4; ++k) {
                        int kn = (k + 1) % 4;
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(qEdge[k], spoke[kn], bEdges[k].reversed(), spoke[k].reversed(), m_degree)));
                    }

                    SplittingCandidate cand;
                    if (tryFinalize(mp, "MedialOGrid", cand)) {
                        registerIfBetter(cand);
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3. Generalized Ribbon / Spine Decomposition (for 2-curve crescents/lenses)
    // ------------------------------------------------------------------------
    if (K == 2) {
        int mSegsList[2] = {3, 4};
        T tipShifts[2] = {T(0.18), T(0.28)};
        for (int mSegs : mSegsList) {
            std::vector<CurveSeg<T>> botSegs(mSegs), topSegs(mSegs);
            CurveSeg<T> rem0 = origSegs[0];
            CurveSeg<T> rem1 = origSegs[1].reversed(); // Both run V0 -> V1
            for (int s = 0; s < mSegs - 1; ++s) {
                T frac = T(1) / T(mSegs - s);
                auto sp0 = rem0.splitArc(frac);
                botSegs[s] = sp0.first; rem0 = sp0.second;
                auto sp1 = rem1.splitArc(frac);
                topSegs[s] = sp1.first; rem1 = sp1.second;
            }
            botSegs[mSegs - 1] = rem0;
            topSegs[mSegs - 1] = rem1;

            int nMid = mSegs - 1;
            for (T delta : tipShifts) {
                std::vector<gsPoint2D<T>> M(nMid);
                for (int i = 0; i < nMid; ++i) {
                    auto pb = botSegs[i].p1();
                    auto pt = topSegs[i].p1();
                    M[i] = gsPoint2D<T>(T(0.5) * (pb.x + pt.x), T(0.5) * (pb.y + pt.y));
                }
                if (nMid >= 2) {
                    gsPoint2D<T> M0_orig = M[0], Mlast_orig = M[nMid - 1];
                    M[0]        = gsPoint2D<T>((1 - delta) * M0_orig.x + delta * M[1].x,
                                               (1 - delta) * M0_orig.y + delta * M[1].y);
                    M[nMid - 1] = gsPoint2D<T>((1 - delta) * Mlast_orig.x + delta * M[nMid - 2].x,
                                               (1 - delta) * Mlast_orig.y + delta * M[nMid - 2].y);
                }

                std::vector<CurveSeg<T>> bSpoke(nMid), tSpoke(nMid), spine(nMid - 1);
                for (int i = 0; i < nMid; ++i) {
                    bSpoke[i] = makeStraightSeg(botSegs[i].p1(), M[i], m_degree);
                    tSpoke[i] = makeStraightSeg(M[i], topSegs[i].p1(), m_degree);
                }
                for (int i = 0; i < nMid - 1; ++i)
                    spine[i] = makeStraightSeg(M[i], M[i + 1], m_degree);

                gsMultiPatch<T> mp;
                // Left tip kite: (V0, B1, M0, T1)
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(botSegs[0], bSpoke[0], tSpoke[0], topSegs[0].reversed(), m_degree)));
                // Middle pairs
                for (int i = 0; i < nMid - 1; ++i) {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(botSegs[i + 1], bSpoke[i + 1], spine[i].reversed(), bSpoke[i].reversed(), m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(spine[i], tSpoke[i + 1], topSegs[i + 1].reversed(), tSpoke[i].reversed(), m_degree)));
                }
                // Right tip kite: (M_last, B_last, V1, T_last)
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(bSpoke[nMid - 1].reversed(), botSegs[mSegs - 1],
                                          topSegs[mSegs - 1].reversed(), tSpoke[nMid - 1].reversed(), m_degree)));

                SplittingCandidate cand;
                if (tryFinalize(mp, "MedialRibbon", cand)) {
                    registerIfBetter(cand);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3b. Curved Channel / Serpentine Ribbon Decomposition (for K == 4 ducts)
    //     Splits along opposite boundary curves into M stations, bisecting any
    //     obtuse (> 140 deg) corners cleanly with Catmull-Clark edge halving.
    // ------------------------------------------------------------------------
    if (K == 4) {
        std::vector<gsPoint2D<T>> densePoly = sampleCurveLoopPoly(outer, 96);
        auto segStrictlyInsideDense = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> bool {
            const int np = static_cast<int>(densePoly.size());
            const T abx = b.x - a.x, aby = b.y - a.y;
            for (int i = 0; i < np; ++i) {
                const auto& c = densePoly[i];
                const auto& d = densePoly[(i + 1) % np];
                const T cdx = d.x - c.x, cdy = d.y - c.y;
                const T denom = abx * cdy - aby * cdx;
                if (std::abs(denom) <= T(1e-14)) continue;
                const T t = ((c.x - a.x) * cdy - (c.y - a.y) * cdx) / denom;
                const T u = ((c.x - a.x) * aby - (c.y - a.y) * abx) / denom;
                if (t > T(0.015) && t < T(0.985) && u >= -T(1e-6) && u <= T(1) + T(1e-6))
                    return false;
            }
            const T sampleT[5] = {T(0.1), T(0.25), T(0.5), T(0.75), T(0.9)};
            for (T t : sampleT) {
                gsPoint2D<T> mid((T(1) - t) * a.x + t * b.x, (T(1) - t) * a.y + t * b.y);
                if (!pointInPolyStrict(mid, densePoly, T(1e-6))) return false;
            }
            return true;
        };

        auto emitTri1to3 = [&](gsMultiPatch<T>& mp,
                               const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
                               const CurveSeg<T>& eAB0, const CurveSeg<T>& eAB1, const gsPoint2D<T>& mAB,
                               const CurveSeg<T>& eBC0, const CurveSeg<T>& eBC1, const gsPoint2D<T>& mBC,
                               const CurveSeg<T>& eCA0, const CurveSeg<T>& eCA1, const gsPoint2D<T>& mCA) -> bool
        {
            gsPoint2D<T> O_mid((mAB.x + mBC.x + mCA.x) / T(3), (mAB.y + mBC.y + mCA.y) / T(3));
            gsPoint2D<T> O_tri((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
            auto check3 = [&](const gsPoint2D<T>& Oc) -> bool {
                gsPoint2D<T> q0[4] = {A, mAB, Oc, mCA};
                gsPoint2D<T> q1[4] = {B, mBC, Oc, mAB};
                gsPoint2D<T> q2[4] = {C, mCA, Oc, mBC};
                return isStrictlyConvexQuad(q0) && isStrictlyConvexQuad(q1) && isStrictlyConvexQuad(q2);
            };
            gsPoint2D<T> O = O_mid;
            if (!check3(O)) {
                if (check3(O_tri)) O = O_tri;
                else return false;
            }
            auto s0 = makeStraightSeg(O, mAB, m_degree);
            auto s1 = makeStraightSeg(O, mBC, m_degree);
            auto s2 = makeStraightSeg(O, mCA, m_degree);
            mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eAB0, s0.reversed(), s2, eCA1, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eBC0, s1.reversed(), s0, eAB1, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eCA0, s2.reversed(), s1, eBC1, m_degree)));
            return true;
        };

        for (int base = 0; base < 2; ++base) {
            CurveSeg<T> botSeg   = origSegs[base];
            CurveSeg<T> rightSeg = origSegs[(base + 1) % 4];
            CurveSeg<T> topSeg   = origSegs[(base + 2) % 4].reversed(); // L[M] -> R[M]
            CurveSeg<T> leftSeg  = origSegs[(base + 3) % 4].reversed(); // L[0] -> L[M]

            const int mVals[4] = {4, 6, 8, 10};
            for (int mIdx = 0; mIdx < 4; ++mIdx) {
                const int M = mVals[mIdx];
                std::vector<CurveSeg<T>> lSegs(M), rSegs(M);
                {
                    CurveSeg<T> remL = leftSeg, remR = rightSeg;
                    for (int j = 0; j < M - 1; ++j) {
                        T frac = T(1) / static_cast<T>(M - j);
                        auto spL = remL.splitArc(frac);
                        auto spR = remR.splitArc(frac);
                        lSegs[j] = spL.first;  remL = spL.second;
                        rSegs[j] = spR.first;  remR = spR.second;
                    }
                    lSegs[M - 1] = remL;
                    rSegs[M - 1] = remR;
                }

                std::vector<gsPoint2D<T>> Lpts(M + 1), Rpts(M + 1);
                std::vector<CurveSeg<T>> hSegs(M + 1);
                Lpts[0] = botSeg.p0(); Rpts[0] = botSeg.p1(); hSegs[0] = botSeg;
                Lpts[M] = topSeg.p0(); Rpts[M] = topSeg.p1(); hSegs[M] = topSeg;

                bool cutsValid = true;
                for (int j = 1; j < M; ++j) {
                    Lpts[j] = lSegs[j].p0();
                    Rpts[j] = rSegs[j].p0();
                    if (!segStrictlyInsideDense(Lpts[j], Rpts[j])) {
                        cutsValid = false;
                        break;
                    }
                    hSegs[j] = makeStraightSeg(Lpts[j], Rpts[j], m_degree);
                }
                if (!cutsValid) continue;

                // Option A: Direct M-quad channel strip
                {
                    gsMultiPatch<T> mpA;
                    for (int j = 0; j < M; ++j) {
                        mpA.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(hSegs[j], rSegs[j], hSegs[j + 1].reversed(), lSegs[j].reversed(), m_degree)));
                    }
                    SplittingCandidate candA;
                    if (tryFinalize(mpA, "CurvedChannelStrip", candA))
                        registerIfBetter(candA);
                }

                // Option B: Conforming diagonal split of obtuse (> 140 deg) end/strip quads + Catmull-Clark midpoint quadification
                std::vector<CurveSeg<T>> h0(M + 1), h1(M + 1), l0(M), l1(M), r0(M), r1(M);
                std::vector<gsPoint2D<T>> hm(M + 1), lm(M), rm(M);
                for (int j = 0; j <= M; ++j) {
                    auto sp = hSegs[j].splitArc(T(0.5));
                    h0[j] = sp.first; h1[j] = sp.second; hm[j] = h0[j].p1();
                }
                for (int j = 0; j < M; ++j) {
                    auto spL = lSegs[j].splitArc(T(0.5));
                    auto spR = rSegs[j].splitArc(T(0.5));
                    l0[j] = spL.first; l1[j] = spL.second; lm[j] = l0[j].p1();
                    r0[j] = spR.first; r1[j] = spR.second; rm[j] = r0[j].p1();
                }

                gsMultiPatch<T> mpB;
                bool okB = true;
                for (int j = 0; j < M && okB; ++j) {
                    std::unique_ptr<gsTensorBSpline<2, T>> testP(
                        makeCoonsFromCCWCycle(hSegs[j], rSegs[j], hSegs[j + 1].reversed(), lSegs[j].reversed(), m_degree));
                    T angs[4], minA = 0, maxA = 180;
                    bool degen = false;
                    evalPatchCornerAnglesDeg(*testP, angs, minA, maxA, degen);

                    gsPoint2D<T> P[4] = {Lpts[j], Rpts[j], Rpts[j + 1], Lpts[j + 1]};
                    int bestObtuse = -1;
                    T worstA = T(140);
                    for (int c = 0; c < 4; ++c) {
                        if (angs[c] > worstA) { worstA = angs[c]; bestObtuse = c; }
                    }

                    if (bestObtuse >= 0 && !degen) {
                        if (bestObtuse == 1 || bestObtuse == 3) {
                            // Diagonal P1 -> P3
                            if (!segStrictlyInsideDense(P[1], P[3])) { okB = false; break; }
                            auto dSeg = makeStraightSeg(P[1], P[3], m_degree);
                            auto dSp  = dSeg.splitArc(T(0.5));
                            gsPoint2D<T> dm = dSp.first.p1();
                            if (!emitTri1to3(mpB, P[0], P[1], P[3],
                                             h0[j], h1[j], hm[j],
                                             dSp.first, dSp.second, dm,
                                             l1[j].reversed(), l0[j].reversed(), lm[j]) ||
                                !emitTri1to3(mpB, P[1], P[2], P[3],
                                             r0[j], r1[j], rm[j],
                                             h1[j + 1].reversed(), h0[j + 1].reversed(), hm[j + 1],
                                             dSp.second.reversed(), dSp.first.reversed(), dm)) {
                                okB = false;
                            }
                        } else {
                            // Diagonal P0 -> P2
                            if (!segStrictlyInsideDense(P[0], P[2])) { okB = false; break; }
                            auto dSeg = makeStraightSeg(P[0], P[2], m_degree);
                            auto dSp  = dSeg.splitArc(T(0.5));
                            gsPoint2D<T> dm = dSp.first.p1();
                            if (!emitTri1to3(mpB, P[0], P[1], P[2],
                                             h0[j], h1[j], hm[j],
                                             r0[j], r1[j], rm[j],
                                             dSp.second.reversed(), dSp.first.reversed(), dm) ||
                                !emitTri1to3(mpB, P[0], P[2], P[3],
                                             dSp.first, dSp.second, dm,
                                             h1[j + 1].reversed(), h0[j + 1].reversed(), hm[j + 1],
                                             l1[j].reversed(), l0[j].reversed(), lm[j])) {
                                okB = false;
                            }
                        }
                    } else {
                        gsPoint2D<T> qCorners[4] = {Lpts[j], Rpts[j], Rpts[j + 1], Lpts[j + 1]};
                        if (!isStrictlyConvexQuad(qCorners)) { okB = false; break; }
                        gsPoint2D<T> O((hm[j].x + rm[j].x + hm[j + 1].x + lm[j].x) * T(0.25),
                                       (hm[j].y + rm[j].y + hm[j + 1].y + lm[j].y) * T(0.25));
                        auto sBot   = makeStraightSeg(O, hm[j], m_degree);
                        auto sRight = makeStraightSeg(O, rm[j], m_degree);
                        auto sTop   = makeStraightSeg(O, hm[j + 1], m_degree);
                        auto sLeft  = makeStraightSeg(O, lm[j], m_degree);
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(h0[j], sBot.reversed(), sLeft, l0[j].reversed(), m_degree)));
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(h1[j], r0[j], sRight.reversed(), sBot, m_degree)));
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(sRight, r1[j], h1[j + 1].reversed(), sTop.reversed(), m_degree)));
                        mpB.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(sLeft.reversed(), sTop, h0[j + 1].reversed(), l1[j].reversed(), m_degree)));
                    }
                }
                if (okB && mpB.nPatches() > 0) {
                    SplittingCandidate candB;
                    if (tryFinalize(mpB, "CurvedChannelSplit", candB))
                        registerIfBetter(candB);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 3c. Single-Reflex Dart / Boomerang Decomposition (K == 4, nReflex == 1)
    //     Splits edges incident to opposite corner opp=(r+2)%4 at frac t,
    //     creating Quad(opp, mA, Vr, mB) + 2 Triangles, then 1-step Catmull-Clark -> 10 quads.
    // ------------------------------------------------------------------------
    if (K == 4 && nReflex == 1) {
        int r = -1;
        for (int i = 0; i < 4; ++i) if (isReflex[i]) { r = i; break; }
        if (r >= 0) {
            int opp = (r + 2) % 4;
            int v1  = (opp + 1) % 4;
            int v3  = (r + 1) % 4;
            const CurveSeg<T>& segOA = origSegs[opp];     // opp -> v1
            const CurveSeg<T>& seg1r = origSegs[v1];      // v1  -> Vr
            const CurveSeg<T>& segr3 = origSegs[r];       // Vr  -> v3
            const CurveSeg<T>& seg3O = origSegs[v3];      // v3  -> opp

            std::vector<gsPoint2D<T>> densePoly = sampleCurveLoopPoly(outer, 48);
            const T fracs[4] = {T(0.50), T(0.35), T(0.45), T(0.25)};
            for (T frac : fracs) {
                auto spA = segOA.splitArc(frac);          // opp -> mA -> v1
                auto spB = seg3O.splitArc(T(1) - frac);   // v3  -> mB -> opp
                gsPoint2D<T> mA = spA.first.p1();
                gsPoint2D<T> mB = spB.first.p1();
                gsPoint2D<T> Vr = verts[r];
                gsPoint2D<T> Vopp = verts[opp];
                gsPoint2D<T> V1 = verts[v1];
                gsPoint2D<T> V3 = verts[v3];

                if (!isSegmentInsidePoly(mA, Vr, densePoly) || !isSegmentInsidePoly(mB, Vr, densePoly)) continue;

                auto cA = makeStraightSeg(mA, Vr, m_degree);
                auto cB = makeStraightSeg(mB, Vr, m_degree);
                auto spEA0 = spA.first.splitArc(T(0.5));
                auto spEA1 = spA.second.splitArc(T(0.5));
                auto sp1r  = seg1r.splitArc(T(0.5));
                auto spr3  = segr3.splitArc(T(0.5));
                auto spEB0 = spB.first.splitArc(T(0.5));
                auto spEB1 = spB.second.splitArc(T(0.5));
                auto spCA  = cA.splitArc(T(0.5));
                auto spCB  = cB.splitArc(T(0.5));

                auto m_eA0 = spEA0.first.p1(), m_eA1 = spEA1.first.p1();
                auto m_1r  = sp1r.first.p1(),  m_r3  = spr3.first.p1();
                auto m_eB0 = spEB0.first.p1(), m_eB1 = spEB1.first.p1();
                auto m_cA  = spCA.first.p1(),  m_cB  = spCB.first.p1();

                auto addTri1to3 = [&](gsMultiPatch<T>& mpOut,
                                      const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
                                      const CurveSeg<T>& eAB0, const CurveSeg<T>& eAB1, const gsPoint2D<T>& mAB,
                                      const CurveSeg<T>& eBC0, const CurveSeg<T>& eBC1, const gsPoint2D<T>& mBC,
                                      const CurveSeg<T>& eCA0, const CurveSeg<T>& eCA1, const gsPoint2D<T>& mCA) -> bool
                {
                    gsPoint2D<T> O;
                    if (!internal_splitter::pickOptimalTri1to3Center(A, B, C, mAB, mBC, mCA, O))
                        return false;
                    auto s0 = makeStraightSeg(O, mAB, m_degree);
                    auto s1 = makeStraightSeg(O, mBC, m_degree);
                    auto s2 = makeStraightSeg(O, mCA, m_degree);
                    mpOut.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eAB0, s0.reversed(), s2, eCA1, m_degree)));
                    mpOut.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eBC0, s1.reversed(), s0, eAB1, m_degree)));
                    mpOut.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(eCA0, s2.reversed(), s1, eBC1, m_degree)));
                    return true;
                };

                // First try 4-triangle split (diagonal Vopp -> Vr) into 12 quads
                if (isSegmentInsidePoly(Vopp, Vr, densePoly)) {
                    auto cDiag = makeStraightSeg(Vopp, Vr, m_degree);
                    auto spDiag = cDiag.splitArc(T(0.5));
                    auto m_diag = spDiag.first.p1();
                    gsMultiPatch<T> mp12;
                    bool ok4 =
                        addTri1to3(mp12, Vopp, mA, Vr,
                                   spEA0.first, spEA0.second, m_eA0,
                                   spCA.first, spCA.second, m_cA,
                                   spDiag.second.reversed(), spDiag.first.reversed(), m_diag) &&
                        addTri1to3(mp12, Vopp, Vr, mB,
                                   spDiag.first, spDiag.second, m_diag,
                                   spCB.second.reversed(), spCB.first.reversed(), m_cB,
                                   spEB1.first, spEB1.second, m_eB1) &&
                        addTri1to3(mp12, mA, V1, Vr,
                                   spEA1.first, spEA1.second, m_eA1,
                                   sp1r.first, sp1r.second, m_1r,
                                   spCA.second.reversed(), spCA.first.reversed(), m_cA) &&
                        addTri1to3(mp12, Vr, V3, mB,
                                   spr3.first, spr3.second, m_r3,
                                   spEB0.first, spEB0.second, m_eB0,
                                   spCB.first, spCB.second, m_cB);
                    if (ok4) {
                        SplittingCandidate cand12;
                        if (tryFinalize(mp12, "DartSplit12", cand12))
                            registerIfBetter(cand12);
                    }
                }

                gsPoint2D<T> qCorners[4] = {Vopp, mA, Vr, mB};
                if (!isStrictlyConvexQuad(qCorners)) continue;

                gsMultiPatch<T> mpD;
                gsPoint2D<T> Oq((m_eA0.x + m_cA.x + m_cB.x + m_eB1.x) * T(0.25),
                                (m_eA0.y + m_cA.y + m_cB.y + m_eB1.y) * T(0.25));
                auto qS0 = makeStraightSeg(Oq, m_eA0, m_degree);
                auto qS1 = makeStraightSeg(Oq, m_cA, m_degree);
                auto qS2 = makeStraightSeg(Oq, m_cB, m_degree);
                auto qS3 = makeStraightSeg(Oq, m_eB1, m_degree);
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(spEA0.first, qS0.reversed(), qS3, spEB1.second, m_degree)));
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(spEA0.second, spCA.first, qS1.reversed(), qS0, m_degree)));
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(qS1, spCA.second, spCB.second.reversed(), qS2.reversed(), m_degree)));
                mpD.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(qS3.reversed(), qS2, spCB.first.reversed(), spEB1.first, m_degree)));

                bool okTriA = addTri1to3(mpD, mA, V1, Vr,
                                         spEA1.first, spEA1.second, m_eA1,
                                         sp1r.first, sp1r.second, m_1r,
                                         spCA.second.reversed(), spCA.first.reversed(), m_cA);
                bool okTriB = addTri1to3(mpD, Vr, V3, mB,
                                         spr3.first, spr3.second, m_r3,
                                         spEB0.first, spEB0.second, m_eB0,
                                         spCB.first, spCB.second, m_cB);
                if (okTriA && okTriB) {
                    SplittingCandidate candD;
                    if (tryFinalize(mpD, "DartSplit10", candD))
                        registerIfBetter(candD);
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // 4. L-Shape Single-Reflex Cut (3 patches, generalized to any reflex index)
    // ------------------------------------------------------------------------
    if (K == 6 && nReflex == 1) {
        int r = -1;
        for (int i = 0; i < K; ++i) if (isReflex[i]) { r = i; break; }
        if (r >= 0) {
            int eA = (r + 2) % 6;
            int eB = (r + 3) % 6;
            auto spA = origSegs[eA].splitArc(T(0.5));
            auto spB = origSegs[eB].splitArc(T(0.5));
            gsPoint2D<T> Vr = verts[r];
            auto cutA = makeStraightSeg(spA.first.p1(), Vr, m_degree);
            auto cutB = makeStraightSeg(spB.first.p1(), Vr, m_degree);

            gsMultiPatch<T> mp;
            // Corner block: (V_{eB+0.5}, V_{r+4}, V_r, V_{eA+0.5}) -> wait trace edges cleanly:
            // spB.first ends at mB; spB.second starts at mB and ends at V_{eB+1}=V_{r+4}
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spB.second, origSegs[(r + 4) % 6], origSegs[(r + 5) % 6], cutB.reversed(), m_degree)));
            // Second arm block: (mA, cutA, cutB_rev, spB.first) -> wait mB -> mA:
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spB.first, cutB, cutA.reversed(), spA.second, m_degree)));
            // Wait check cycle order: spA.second (mA->V_{eA+1}=mB? No, eB=eA+1! So spA.second ends at V_{eB}=start of spB.first!):
            // Corner at mB: spB.first (V_{eB} -> mB), cutB (mB -> Vr), cutA.reversed() (Vr -> mA), spA.second (mA -> V_{eB})! That IS CCW!
            // Third block: (V_r, origSegs[r], origSegs[(r+1)%6], spA.first, cutA):
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(origSegs[r], origSegs[(r + 1) % 6], spA.first, cutA, m_degree)));

            SplittingCandidate cand;
            if (tryFinalize(mp, "ReflexCut3Patch", cand)) {
                registerIfBetter(cand);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 5a. Alternating Star Kites (K even >= 6, K/2 patches, spokes to reflex vertices)
    // ------------------------------------------------------------------------
    if (K >= 6 && (K % 2 == 0) && nReflex == K / 2 && pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        int startIdx = isReflex[0] ? 0 : 1;
        bool altOk = true;
        for (int i = 0; i < K; ++i) {
            bool expectRef = ((i % 2) == startIdx);
            if (isReflex[i] != expectRef) { altOk = false; break; }
        }
        if (altOk) {
            std::vector<CurveSeg<T>> spokes(K / 2);
            for (int j = 0; j < K / 2; ++j) {
                int vIdx = (startIdx + 2 * j) % K;
                spokes[j] = makeStraightSeg(centroid, verts[vIdx], m_degree);
            }
            gsMultiPatch<T> mp;
            for (int j = 0; j < K / 2; ++j) {
                int jn = (j + 1) % (K / 2);
                int e0 = (startIdx + 2 * j) % K;
                int e1 = (e0 + 1) % K;
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(spokes[j], origSegs[e0], origSegs[e1], spokes[jn].reversed(), m_degree)));
            }
            SplittingCandidate cand;
            if (tryFinalize(mp, "AlternatingStarKites", cand)) {
                registerIfBetter(cand);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 5b. Centroid-to-Midpoint Star Fan (K patches, for 0-reflex star domains)
    // ------------------------------------------------------------------------
    if (K >= 3 && nReflex == 0 && pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        std::vector<CurveSeg<T>> h0(K), h1(K), spokes(K);
        for (int i = 0; i < K; ++i) {
            auto sp = origSegs[i].splitArc(T(0.5));
            h0[i] = sp.first;
            h1[i] = sp.second;
            spokes[i] = makeStraightSeg(centroid, h0[i].p1(), m_degree);
        }
        gsMultiPatch<T> mp;
        for (int i = 0; i < K; ++i) {
            int iprev = (i + K - 1) % K;
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(spokes[iprev], h1[iprev], h0[i], spokes[i].reversed(), m_degree)));
        }
        SplittingCandidate cand;
        if (tryFinalize(mp, "CentroidStarFan", cand)) {
            registerIfBetter(cand);
        }
    }

    // ------------------------------------------------------------------------
    // 5c. Generalized Reflex-Sector Star Kite Decomposition
    //     Connects an interior center C to all reflex vertices; sectors of
    //     length s=1,2 become single kite quads, and longer sectors subdivide
    //     via interior spokes to boundary edge midpoints (e.g. 7-star -> 4 quads).
    // ------------------------------------------------------------------------
    if (nReflex >= 2) {
        std::vector<int> rIdx;
        for (int i = 0; i < K; ++i) {
            if (isReflex[i]) rIdx.push_back(i);
        }
        gsPoint2D<T> rAvg(0, 0);
        for (int ri : rIdx) {
            rAvg.x += verts[ri].x;
            rAvg.y += verts[ri].y;
        }
        rAvg.x /= static_cast<T>(rIdx.size());
        rAvg.y /= static_cast<T>(rIdx.size());

        std::vector<gsPoint2D<T>> starCenters = {rAvg, centroid};
        for (const auto& ctr : starCenters) {
            if (!pointInPolyStrict(ctr, outerPoly, T(1e-5))) continue;

            gsMultiPatch<T> mp;
            bool okStar = true;
            std::vector<CurveSeg<T>> rSpokes(rIdx.size());
            for (size_t a = 0; a < rIdx.size(); ++a) {
                if (!isSegmentInsidePoly(ctr, verts[rIdx[a]], outerPoly)) {
                    okStar = false;
                    break;
                }
                rSpokes[a] = makeStraightSeg(ctr, verts[rIdx[a]], m_degree);
            }
            if (!okStar) continue;

            for (size_t a = 0; a < rIdx.size() && okStar; ++a) {
                size_t b = (a + 1) % rIdx.size();
                int rA = rIdx[a];
                int rB = rIdx[b];
                int span = (rB - rA + K) % K;
                if (span == 0) span = K;

                if (span == 1) {
                    auto sp = origSegs[rA].splitArc(T(0.5));
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(rSpokes[a], sp.first, sp.second, rSpokes[b].reversed(), m_degree)));
                } else if (span == 2) {
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(rSpokes[a], origSegs[rA], origSegs[(rA + 1) % K], rSpokes[b].reversed(), m_degree)));
                } else {
                    // Walk sector in 2-edge steps; if 3 edges remain at end, split middle edge.
                    int pos = 0;
                    CurveSeg<T> curSpoke = rSpokes[a];
                    while (span - pos > 3) {
                        int e0 = (rA + pos) % K;
                        int e1 = (rA + pos + 1) % K;
                        int vNext = (rA + pos + 2) % K;
                        if (!isSegmentInsidePoly(ctr, verts[vNext], outerPoly)) {
                            okStar = false;
                            break;
                        }
                        CurveSeg<T> nextSpoke = makeStraightSeg(ctr, verts[vNext], m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(curSpoke, origSegs[e0], origSegs[e1], nextSpoke.reversed(), m_degree)));
                        curSpoke = nextSpoke;
                        pos += 2;
                    }
                    if (!okStar) break;
                    int rem = span - pos;
                    if (rem == 2) {
                        int e0 = (rA + pos) % K;
                        int e1 = (rA + pos + 1) % K;
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(curSpoke, origSegs[e0], origSegs[e1], rSpokes[b].reversed(), m_degree)));
                    } else if (rem == 3) {
                        int e0 = (rA + pos) % K;
                        int e1 = (rA + pos + 1) % K;
                        int e2 = (rA + pos + 2) % K;
                        auto spMid = origSegs[e1].splitArc(T(0.5));
                        gsPoint2D<T> mPt = spMid.first.p1();
                        if (!isSegmentInsidePoly(ctr, mPt, outerPoly)) {
                            okStar = false;
                            break;
                        }
                        CurveSeg<T> mSpoke = makeStraightSeg(ctr, mPt, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(curSpoke, origSegs[e0], spMid.first, mSpoke.reversed(), m_degree)));
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(mSpoke, spMid.second, origSegs[e2], rSpokes[b].reversed(), m_degree)));
                    } else {
                        okStar = false;
                    }
                }
            }
            if (!okStar) continue;

            SplittingCandidate cand;
            if (tryFinalize(mp, "ReflexStarKites", cand)) {
                registerIfBetter(cand);
            }
        }
    }

    // ------------------------------------------------------------------------
    // 6. Comb / Fingered-Channel Spine Decomposition
    //    Detects T-tooth comb domains (K = 4*T, T >= 2) and tries both:
    //    (a) 3T - 3 conforming patches (diagonal end teeth -> 9 patches for T=4)
    //    (b) 3T - 1 conforming patches (orthogonal grid -> 11 patches for T=4)
    // ------------------------------------------------------------------------
    if (K >= 8 && (K % 4 == 0) && nReflex == 2 * (K / 4 - 1)) {
        const int Tcnt = K / 4;
        for (int kBot = 0; kBot < K; ++kBot) {
            auto idxOf = [&](int off) -> int { return (kBot + off) % K; };
            bool patOk = !isReflex[idxOf(0)] && !isReflex[idxOf(1)] &&
                         !isReflex[idxOf(2)] && !isReflex[idxOf(3)];
            for (int m = 0; m < Tcnt - 1 && patOk; ++m) {
                int baseOff = 4 + 4 * m;
                if (!isReflex[idxOf(baseOff)] || !isReflex[idxOf(baseOff + 1)] ||
                    isReflex[idxOf(baseOff + 2)] || isReflex[idxOf(baseOff + 3)]) {
                    patOk = false;
                }
            }
            if (!patOk) continue;

            // Collect valley floor points ordered right-to-left (R[t], L[t] for t = 0..Tcnt-2)
            // and split parameters along bottom curve eBot = origSegs[kBot] (ordered ascending left-to-right).
            const CurveSeg<T>& eBot = origSegs[idxOf(0)];
            const CurveSeg<T>& eR   = origSegs[idxOf(1)];
            const CurveSeg<T>& eL   = origSegs[idxOf(K - 1)];

            std::vector<gsPoint2D<T>> rPts(Tcnt - 1), lPts(Tcnt - 1);
            std::vector<T> rParams(Tcnt - 1), lParams(Tcnt - 1);
            for (int t = 0; t < Tcnt - 1; ++t) {
                rPts[t] = verts[idxOf(4 + 4 * t)];
                lPts[t] = verts[idxOf(5 + 4 * t)];
                rParams[t] = eBot.closestParam(rPts[t]);
                lParams[t] = eBot.closestParam(lPts[t]);
            }

            // Build ordered valley vertices left-to-right: L_{T-2}, R_{T-2}, ..., L_0, R_0
            std::vector<gsPoint2D<T>> topNodes;
            std::vector<T> botSplitU;
            for (int t = Tcnt - 2; t >= 0; --t) {
                topNodes.push_back(lPts[t]);
                botSplitU.push_back(lParams[t]);
                topNodes.push_back(rPts[t]);
                botSplitU.push_back(rParams[t]);
            }

            bool monoU = !botSplitU.empty() && botSplitU.front() > T(0.02) && botSplitU.back() < T(0.98);
            for (size_t i = 1; i < botSplitU.size() && monoU; ++i) {
                if (botSplitU[i] <= botSplitU[i - 1] + T(0.02)) monoU = false;
            }
            if (!monoU) continue;

            auto splitMulti = [&](const CurveSeg<T>& crv, const std::vector<T>& us) -> std::vector<CurveSeg<T>> {
                std::vector<CurveSeg<T>> segs;
                CurveSeg<T> rem = crv;
                T prevU = T(0);
                for (T u : us) {
                    T localU = (u - prevU) / std::max(T(1e-9), T(1) - prevU);
                    auto sp = rem.splitAtParam(localU);
                    segs.push_back(sp.first);
                    rem = sp.second;
                    prevU = u;
                }
                segs.push_back(rem);
                return segs;
            };

            // (a) Diagonal-end comb decomposition: 3*Tcnt - 3 patches (no side splits)
            {
                std::vector<T> innerU;
                for (size_t i = 1; i + 1 < botSplitU.size(); ++i)
                    innerU.push_back(botSplitU[i]);
                auto botSegs9 = splitMulti(eBot, innerU); // size = 2*Tcnt - 3

                bool diagOk = true;
                gsPoint2D<T> vBotL = verts[idxOf(0)];
                gsPoint2D<T> vBotR = verts[idxOf(1)];
                if (!isSegmentInsidePoly(vBotL, topNodes.front(), outerPoly) ||
                    !isSegmentInsidePoly(topNodes.back(), vBotR, outerPoly)) {
                    diagOk = false;
                }
                for (size_t i = 1; i + 1 < topNodes.size() && diagOk; ++i) {
                    if (!isSegmentInsidePoly(botSegs9[i - 1].p1(), topNodes[i], outerPoly))
                        diagOk = false;
                }

                if (diagOk) {
                    gsMultiPatch<T> mp;
                    CurveSeg<T> diagL = makeStraightSeg(vBotL, topNodes.front(), m_degree);
                    CurveSeg<T> diagR = makeStraightSeg(topNodes.back(), vBotR, m_degree);

                    // Leftmost tooth (t = Tcnt - 1) and rightmost tooth (t = 0)
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(diagL, origSegs[idxOf(K - 3)], origSegs[idxOf(K - 2)], eL, m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(diagR, eR, origSegs[idxOf(2)], origSegs[idxOf(3)], m_degree)));

                    // Tooth bases left-to-right between topNodes[2*j-1] (= R_{...}) and topNodes[2*j] (= L_{...})
                    int nMidTop = static_cast<int>(topNodes.size()) - 1; // 2*Tcnt - 3
                    std::vector<CurveSeg<T>> spineTop(nMidTop);
                    for (int m = 0; m < nMidTop; ++m) {
                        if (m % 2 == 0) {
                            // Valley floor right-to-left index: gapIdx = (Tcnt - 2) - m/2
                            int gapIdx = (Tcnt - 2) - m / 2;
                            spineTop[m] = origSegs[idxOf(4 + 4 * gapIdx)].reversed();
                        } else {
                            // Middle tooth base across inner tooth t = (Tcnt - 1) - (m + 1)/2
                            int toothIdx = (Tcnt - 1) - (m + 1) / 2;
                            spineTop[m] = makeStraightSeg(topNodes[m], topNodes[m + 1], m_degree);
                            int eDropR = idxOf(3 + 4 * toothIdx);
                            int eTop   = idxOf(2 + 4 * toothIdx);
                            int eRiseL = idxOf(5 + 4 * (toothIdx - 1));
                            mp.addPatch(typename gsGeometry<T>::uPtr(
                                makeCoonsFromCCWCycle(spineTop[m], origSegs[eRiseL], origSegs[eTop], origSegs[eDropR], m_degree)));
                        }
                    }

                    // Spine quads left-to-right (m = 0 .. nMidTop - 1)
                    std::vector<CurveSeg<T>> pillars(nMidTop + 1);
                    pillars[0] = diagL.reversed();
                    for (int m = 1; m < nMidTop; ++m)
                        pillars[m] = makeStraightSeg(botSegs9[m - 1].p1(), topNodes[m], m_degree);
                    pillars[nMidTop] = diagR;

                    for (int m = 0; m < nMidTop; ++m) {
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(botSegs9[m], pillars[m + 1], spineTop[m].reversed(), pillars[m].reversed(), m_degree)));
                    }

                    SplittingCandidate cand;
                    if (tryFinalize(mp, "CombDiagonal", cand)) {
                        registerIfBetter(cand);
                    }
                }
            }

            // (b) Orthogonal-spine comb decomposition: 3*Tcnt - 1 patches
            {
                T uL = eL.closestParam(topNodes.front());
                T uR = eR.closestParam(topNodes.back());
                if (uL > T(0.05) && uL < T(0.95) && uR > T(0.05) && uR < T(0.95)) {
                    auto spL = eL.splitAtParam(uL); // spL.first: v(K-1)->midL, spL.second: midL->v(0)
                    auto spR = eR.splitAtParam(uR); // spR.first: v(1)->midR, spR.second: midR->v(2)
                    auto botSegs11 = splitMulti(eBot, botSplitU); // size = 2*Tcnt - 1

                    gsPoint2D<T> midL = spL.first.p1();
                    gsPoint2D<T> midR = spR.first.p1();

                    std::vector<gsPoint2D<T>> allTopNodes;
                    allTopNodes.push_back(midL);
                    allTopNodes.insert(allTopNodes.end(), topNodes.begin(), topNodes.end());
                    allTopNodes.push_back(midR);

                    int nCols = 2 * Tcnt - 1;
                    bool gridOk = isSegmentInsidePoly(midL, topNodes.front(), outerPoly) &&
                                  isSegmentInsidePoly(topNodes.back(), midR, outerPoly);
                    for (int m = 1; m < nCols && gridOk; ++m) {
                        if (!isSegmentInsidePoly(botSegs11[m - 1].p1(), allTopNodes[m], outerPoly))
                            gridOk = false;
                    }

                    if (gridOk) {
                        gsMultiPatch<T> mp;
                        std::vector<CurveSeg<T>> spineTop(nCols);
                        for (int m = 0; m < nCols; ++m) {
                            if (m % 2 == 1) {
                                int gapIdx = (Tcnt - 2) - (m - 1) / 2;
                                spineTop[m] = origSegs[idxOf(4 + 4 * gapIdx)].reversed();
                            } else {
                                spineTop[m] = makeStraightSeg(allTopNodes[m], allTopNodes[m + 1], m_degree);
                                int toothIdx = (Tcnt - 1) - m / 2;
                                if (toothIdx == Tcnt - 1) {
                                    mp.addPatch(typename gsGeometry<T>::uPtr(
                                        makeCoonsFromCCWCycle(spineTop[m], origSegs[idxOf(K - 3)], origSegs[idxOf(K - 2)], spL.first, m_degree)));
                                } else if (toothIdx == 0) {
                                    mp.addPatch(typename gsGeometry<T>::uPtr(
                                        makeCoonsFromCCWCycle(spineTop[m], spR.second, origSegs[idxOf(2)], origSegs[idxOf(3)], m_degree)));
                                } else {
                                    int eDropR = idxOf(3 + 4 * toothIdx);
                                    int eTop   = idxOf(2 + 4 * toothIdx);
                                    int eRiseL = idxOf(5 + 4 * (toothIdx - 1));
                                    mp.addPatch(typename gsGeometry<T>::uPtr(
                                        makeCoonsFromCCWCycle(spineTop[m], origSegs[eRiseL], origSegs[eTop], origSegs[eDropR], m_degree)));
                                }
                            }
                        }

                        std::vector<CurveSeg<T>> pillars(nCols + 1);
                        pillars[0] = spL.second.reversed();
                        for (int m = 1; m < nCols; ++m)
                            pillars[m] = makeStraightSeg(botSegs11[m - 1].p1(), allTopNodes[m], m_degree);
                        pillars[nCols] = spR.first;

                        for (int m = 0; m < nCols; ++m) {
                            mp.addPatch(typename gsGeometry<T>::uPtr(
                                makeCoonsFromCCWCycle(botSegs11[m], pillars[m + 1], spineTop[m].reversed(), pillars[m].reversed(), m_degree)));
                        }

                        SplittingCandidate cand;
                        if (tryFinalize(mp, "CombOrthogonal", cand)) {
                            registerIfBetter(cand);
                        }
                    }
                }
            }
            break;
        }
    }

    return bestCand;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeHoleFreeTriangulationSplitting()
{
    using namespace internal_splitter;

    SplittingCandidate bestCand;
    bestCand.methodUsed = "HarmonicCrossField";
    bestCand.isValid = false;

    const auto& outer = m_domain.outer();
    const int K = outer.size();
    if (K == 0) return bestCand;

    std::vector<gsPoint2D<T>> outerPoly = sampleCurveLoopPoly(outer, 24);
    gsPoint2D<T> centroid = polygonCentroid(outerPoly);

    auto candCost = [&](const SplittingCandidate& c) -> T {
        return m_weightPatches * T(c.numPatches) + m_weightDistortion * c.distortionEnergy;
    };
    auto updateBest = [&](const SplittingCandidate& c) {
        if (!c.isValid) return;
        if (!bestCand.isValid || candCost(c) < candCost(bestCand))
            bestCand = c;
    };

    auto tryFinalize = [&](gsMultiPatch<T>& mp, SplittingCandidate& out) -> bool {
        if (mp.nPatches() == 0) return false;
        optimizeWinslowBarrier(mp);
        mp.computeTopology(T(1e-5), false, false);
        const T domMinAng = getDomainMinCornerAngleDeg();
        const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
        const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
        if (!validateMultiPatch(m_domain, mp, T(1e-5), nullptr, minAngTol, maxAngTol))
            return false;
        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
        out.multiPatch = mp;
        out.numPatches = static_cast<int>(mp.nPatches());
        out.minJacobianDet = minDet;
        out.distortionEnergy = energy;
        out.minCornerAngleDeg = minA;
        out.maxCornerAngleDeg = maxA;
        out.methodUsed = "HarmonicCrossField";
        out.isValid = true;
        return true;
    };

    struct Tri { int v[3]; int e[3]; };  // edges: v0->v1, v1->v2, v2->v0
    struct Quad { int v[4]; int e[4]; }; // edges: v0->v1, v1->v2, v2->v3, v3->v0
    struct EdgeRec {
        int a, b;
        CurveSeg<T> seg; // oriented a -> b
        CurveSeg<T> h0, h1;
        gsPoint2D<T> mid;
        bool splitDone = false;
        bool isBoundary = false;
    };
    struct MixedMesh {
        std::vector<gsPoint2D<T>> pts;
        std::vector<bool> isSmooth;
        std::vector<EdgeRec> edges;
        std::vector<Tri> tris;
        std::vector<Quad> quads;
    };

    auto ensureSplit = [&](MixedMesh& mesh, int eIdx) {
        auto& er = mesh.edges[eIdx];
        if (!er.splitDone) {
            auto sp = er.seg.splitArc(T(0.5));
            er.h0 = sp.first;
            er.h1 = sp.second;
            er.mid = er.h0.p1();
            er.splitDone = true;
        }
    };

    auto getSegDir = [&](const MixedMesh& mesh, int eIdx, int u, int v) -> CurveSeg<T> {
        const auto& er = mesh.edges[eIdx];
        return (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
    };

    auto getHalvesDir = [&](const MixedMesh& mesh, int eIdx, int u, int v,
                            CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
        const auto& er = mesh.edges[eIdx];
        outMid = er.mid;
        if (er.a == u && er.b == v) {
            out0 = er.h0;
            out1 = er.h1;
        } else {
            out0 = er.h1.reversed();
            out1 = er.h0.reversed();
        }
    };

    // Shared quadifier for any conforming Tri+Quad mixed mesh given edge-split flags.
    // Every triangle must have all 3 edges split; every quad must have opposite pairs matched.
    auto emitQuadPatches = [&](MixedMesh& mesh,
                               const std::vector<bool>& needSplit,
                               SplittingCandidate& candOut) -> bool
    {
        for (size_t i = 0; i < mesh.edges.size(); ++i) {
            if (needSplit[i]) ensureSplit(mesh, static_cast<int>(i));
        }

        gsMultiPatch<T> mp;

        // 1. Emit unmatched triangles (1-to-3 Catmull-Clark split)
        for (const auto& t : mesh.tris) {
            if (!needSplit[t.e[0]] || !needSplit[t.e[1]] || !needSplit[t.e[2]]) return false;
            int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
            CurveSeg<T> e01_a, e01_b, e12_a, e12_b, e20_a, e20_b;
            gsPoint2D<T> m01, m12, m20;
            getHalvesDir(mesh, t.e[0], v0, v1, e01_a, e01_b, m01);
            getHalvesDir(mesh, t.e[1], v1, v2, e12_a, e12_b, m12);
            getHalvesDir(mesh, t.e[2], v2, v0, e20_a, e20_b, m20);

            gsPoint2D<T> O;
            if (!internal_splitter::pickOptimalTri1to3Center(mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                                                             m01, m12, m20, O))
                return false;

            auto s01 = makeStraightSeg(O, m01, m_degree);
            auto s12 = makeStraightSeg(O, m12, m_degree);
            auto s20 = makeStraightSeg(O, m20, m_degree);
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(e01_a, s01.reversed(), s20, e20_b, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(e12_a, s12.reversed(), s01, e01_b, m_degree)));
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(e20_a, s20.reversed(), s12, e12_b, m_degree)));
        }

        // 2. Emit quads (0-split, 2-split along one opposite pair, or 4-split)
        for (const auto& q : mesh.quads) {
            bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
            bool sp2 = needSplit[q.e[2]], sp3 = needSplit[q.e[3]];
            if (sp0 != sp2 || sp1 != sp3) return false;
            int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];

            if (!sp0 && !sp1) {
                gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                if (!isStrictlyConvexQuad(c)) return false;
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(getSegDir(mesh, q.e[0], v0, v1),
                                          getSegDir(mesh, q.e[1], v1, v2),
                                          getSegDir(mesh, q.e[2], v2, v3),
                                          getSegDir(mesh, q.e[3], v3, v0),
                                          m_degree)));
            } else if (sp0 && !sp1) {
                CurveSeg<T> e0a, e0b, e2a, e2b;
                gsPoint2D<T> m01, m23;
                getHalvesDir(mesh, q.e[0], v0, v1, e0a, e0b, m01);
                getHalvesDir(mesh, q.e[2], v2, v3, e2a, e2b, m23);
                gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return false;
                auto sMid = makeStraightSeg(m01, m23, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e0a, sMid, e2b, getSegDir(mesh, q.e[3], v3, v0), m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e0b, getSegDir(mesh, q.e[1], v1, v2), e2a, sMid.reversed(), m_degree)));
            } else if (!sp0 && sp1) {
                CurveSeg<T> e1a, e1b, e3a, e3b;
                gsPoint2D<T> m12, m30;
                getHalvesDir(mesh, q.e[1], v1, v2, e1a, e1b, m12);
                getHalvesDir(mesh, q.e[3], v3, v0, e3a, e3b, m30);
                gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return false;
                auto sMid = makeStraightSeg(m12, m30, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(getSegDir(mesh, q.e[0], v0, v1), e1a, sMid, e3b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(sMid.reversed(), e1b, getSegDir(mesh, q.e[2], v2, v3), e3a, m_degree)));
            } else {
                CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b;
                gsPoint2D<T> m01, m12, m23, m30;
                getHalvesDir(mesh, q.e[0], v0, v1, e0a, e0b, m01);
                getHalvesDir(mesh, q.e[1], v1, v2, e1a, e1b, m12);
                getHalvesDir(mesh, q.e[2], v2, v3, e2a, e2b, m23);
                getHalvesDir(mesh, q.e[3], v3, v0, e3a, e3b, m30);

                gsPoint2D<T> O_mid((m01.x + m12.x + m23.x + m30.x) / T(4),
                                   (m01.y + m12.y + m23.y + m30.y) / T(4));
                gsPoint2D<T> O_v((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x + mesh.pts[v3].x) / T(4),
                                 (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y + mesh.pts[v3].y) / T(4));
                auto check4 = [&](const gsPoint2D<T>& Oc) -> bool {
                    gsPoint2D<T> q0[4] = {mesh.pts[v0], m01, Oc, m30};
                    gsPoint2D<T> q1[4] = {mesh.pts[v1], m12, Oc, m01};
                    gsPoint2D<T> q2[4] = {mesh.pts[v2], m23, Oc, m12};
                    gsPoint2D<T> q3[4] = {mesh.pts[v3], m30, Oc, m23};
                    return isStrictlyConvexQuad(q0) && isStrictlyConvexQuad(q1)
                        && isStrictlyConvexQuad(q2) && isStrictlyConvexQuad(q3);
                };
                gsPoint2D<T> O = O_mid;
                if (!check4(O)) {
                    if (check4(O_v)) {
                        O = O_v;
                    } else {
                        bool found = false;
                        for (T w = T(0.1); w <= T(0.91); w += T(0.1)) {
                            gsPoint2D<T> Otry((1 - w) * O_mid.x + w * O_v.x,
                                              (1 - w) * O_mid.y + w * O_v.y);
                            if (check4(Otry)) { O = Otry; found = true; break; }
                        }
                        if (!found) return false;
                    }
                }
                auto s01 = makeStraightSeg(O, m01, m_degree);
                auto s12 = makeStraightSeg(O, m12, m_degree);
                auto s23 = makeStraightSeg(O, m23, m_degree);
                auto s30 = makeStraightSeg(O, m30, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e0a, s01.reversed(), s30, e3b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e1a, s12.reversed(), s01, e0b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e2a, s23.reversed(), s12, e1b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(
                    makeCoonsFromCCWCycle(e3a, s30.reversed(), s23, e2b, m_degree)));
            }
        }
        return tryFinalize(mp, candOut);
    };

    // Opposite-edge parity propagation across quads starting from unmatched triangles.
    auto propagateOppositeSplits = [&](const MixedMesh& mesh, std::vector<bool>& needSplit) {
        needSplit.assign(mesh.edges.size(), false);
        for (const auto& t : mesh.tris) {
            needSplit[t.e[0]] = true;
            needSplit[t.e[1]] = true;
            needSplit[t.e[2]] = true;
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& q : mesh.quads) {
                if (needSplit[q.e[0]] != needSplit[q.e[2]]) {
                    needSplit[q.e[0]] = needSplit[q.e[2]] = true;
                    changed = true;
                }
                if (needSplit[q.e[1]] != needSplit[q.e[3]]) {
                    needSplit[q.e[1]] = needSplit[q.e[3]] = true;
                    changed = true;
                }
            }
        }
    };

    // Ensure no smooth (G1 ~180 deg) boundary vertex has triangle degree == 1 by bisecting
    // the opposite interior edge (Rivara 1-edge bisection preserves conforming mesh).
    auto fixSmoothDegreeOneVertices = [&](MixedMesh& mesh) {
        for (int iter = 0; iter < 8; ++iter) {
            std::vector<int> vDeg(mesh.pts.size(), 0);
            for (const auto& t : mesh.tris) {
                ++vDeg[t.v[0]]; ++vDeg[t.v[1]]; ++vDeg[t.v[2]];
            }
            int badTri = -1, badCorner = -1;
            for (int ti = 0; ti < static_cast<int>(mesh.tris.size()) && badTri < 0; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    int vi = mesh.tris[ti].v[k];
                    if (vDeg[vi] != 1) continue;
                    int vn = mesh.tris[ti].v[(k + 1) % 3];
                    int vp = mesh.tris[ti].v[(k + 2) % 3];
                    T dx1 = mesh.pts[vn].x - mesh.pts[vi].x, dy1 = mesh.pts[vn].y - mesh.pts[vi].y;
                    T dx2 = mesh.pts[vp].x - mesh.pts[vi].x, dy2 = mesh.pts[vp].y - mesh.pts[vi].y;
                    T l1 = std::hypot(dx1, dy1), l2 = std::hypot(dx2, dy2);
                    T cosA = (l1 > T(1e-14) && l2 > T(1e-14)) ? (dx1 * dx2 + dy1 * dy2) / (l1 * l2) : T(1);
                    bool obtuseCorner = (cosA < T(-0.766)); // > 140 deg
                    if (mesh.isSmooth[vi] || obtuseCorner) {
                        int oppEdge = mesh.tris[ti].e[(k + 1) % 3];
                        if (!mesh.edges[oppEdge].isBoundary) {
                            badTri = ti;
                            badCorner = k;
                            break;
                        }
                    }
                }
            }
            if (badTri < 0) break;

            Tri t0 = mesh.tris[badTri];
            int v = t0.v[badCorner];
            int n = t0.v[(badCorner + 1) % 3];
            int p = t0.v[(badCorner + 2) % 3];
            int e_vn = t0.e[badCorner];
            int e_np = t0.e[(badCorner + 1) % 3];
            int e_pv = t0.e[(badCorner + 2) % 3];

            // Find neighbor triangle across e_np
            int mateTri = -1, mateEdgePos = -1;
            for (int ti = 0; ti < static_cast<int>(mesh.tris.size()); ++ti) {
                if (ti == badTri) continue;
                for (int k = 0; k < 3; ++k) {
                    if (mesh.tris[ti].e[k] == e_np) {
                        mateTri = ti;
                        mateEdgePos = k;
                        break;
                    }
                }
                if (mateTri >= 0) break;
            }
            if (mateTri < 0) break;

            Tri t1 = mesh.tris[mateTri];
            int w = t1.v[(mateEdgePos + 2) % 3];
            int e_pw = t1.e[(mateEdgePos + 1) % 3]; // edge from n->w or p->w depending on winding
            int e_wn = t1.e[(mateEdgePos + 2) % 3];
            // Wait: in t1, e[mateEdgePos] goes p -> n (opposite to t0's n -> p),
            // so t1.v[mateEdgePos] == p, t1.v[(mateEdgePos+1)%3] == n, w == t1.v[(mateEdgePos+2)%3].
            if (t1.v[mateEdgePos] != p || t1.v[(mateEdgePos + 1) % 3] != n) {
                // Adjust edge names to match (p, n, w) ordering
                e_pw = t1.e[(mateEdgePos + 1) % 3];
                e_wn = t1.e[(mateEdgePos + 2) % 3];
            } else {
                e_pw = t1.e[(mateEdgePos + 2) % 3];
                e_wn = t1.e[(mateEdgePos + 1) % 3];
            }

            ensureSplit(mesh, e_np);
            int mIdx = static_cast<int>(mesh.pts.size());
            mesh.pts.push_back(mesh.edges[e_np].mid);
            mesh.isSmooth.push_back(false);

            CurveSeg<T> s_nm, s_mp;
            gsPoint2D<T> midTmp;
            getHalvesDir(mesh, e_np, n, p, s_nm, s_mp, midTmp);

            auto addEdge = [&](int a, int b, const CurveSeg<T>& s) -> int {
                int id = static_cast<int>(mesh.edges.size());
                mesh.edges.push_back({a, b, s, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false});
                return id;
            };
            int e_nm = addEdge(n, mIdx, s_nm);
            int e_mp = addEdge(mIdx, p, s_mp);
            int e_vm = addEdge(v, mIdx, makeStraightSeg(mesh.pts[v], mesh.pts[mIdx], m_degree));
            int e_wm = addEdge(w, mIdx, makeStraightSeg(mesh.pts[w], mesh.pts[mIdx], m_degree));

            // Replace badTri and mateTri cleanly
            mesh.tris[badTri]  = Tri{{v, n, mIdx}, {e_vn, e_nm, e_vm}};
            mesh.tris[mateTri] = Tri{{p, mIdx, w}, {e_mp, e_wm, e_pw}};
            mesh.tris.push_back(Tri{{p, v, mIdx}, {e_pv, e_vm, e_mp}});
            mesh.tris.push_back(Tri{{mIdx, n, w}, {e_nm, e_wn, e_wm}});
        }
    };

    // Lawson Delaunay edge flips on interior straight edges inside outerPoly.
    auto applyLawsonDelaunayFlips = [&](MixedMesh& mesh) {
        int nT = static_cast<int>(mesh.tris.size());
        int nE = static_cast<int>(mesh.edges.size());
        if (nT < 2 || nE == 0) return;
        for (int pass = 0; pass < nT * 2; ++pass) {
            struct EdgeAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
            std::vector<EdgeAdj> eAdj(nE);
            std::vector<int> vertDeg(mesh.pts.size(), 0);
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    ++vertDeg[mesh.tris[ti].v[k]];
                    int eid = mesh.tris[ti].e[k];
                    if (eid >= 0 && eid < nE) {
                        if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
                        else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
                    }
                }
            }
            bool flipped = false;
            for (int eIdx = 0; eIdx < nE; ++eIdx) {
                if (mesh.edges[eIdx].isBoundary) continue;
                int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
                int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
                if (tA < 0 || tB < 0) continue;
                Tri A = mesh.tris[tA], B = mesh.tris[tB];
                int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
                int d = B.v[(kB + 2) % 3];
                gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
                if (!isStrictlyConvexQuad(qa)) continue;
                if (!isSegmentInsidePoly(mesh.pts[c], mesh.pts[d], outerPoly)) continue;

                T curMin0, curMax0, curMin1, curMax1;
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[a], mesh.pts[b], mesh.pts[c], curMin0, curMax0);
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[b], mesh.pts[a], mesh.pts[d], curMin1, curMax1);
                T newMin0, newMax0, newMin1, newMax1;
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[a], mesh.pts[d], newMin0, newMax0);
                internal_splitter::triMinMaxAnglesDeg(mesh.pts[c], mesh.pts[d], mesh.pts[b], newMin1, newMax1);
                T curMin = std::min(curMin0, curMin1), newMin = std::min(newMin0, newMin1);
                T curMax = std::max(curMax0, curMax1), newMax = std::max(newMax0, newMax1);
                if (newMin <= curMin + T(1e-3) && !(std::abs(newMin - curMin) <= T(1e-3) && newMax < curMax - T(1e-2)))
                    continue;

                // Protect smooth boundary vertices from dropping to deg == 1 unless needed
                if ((mesh.isSmooth[a] && vertDeg[a] <= 2) || (mesh.isSmooth[b] && vertDeg[b] <= 2)) continue;

                // Find edge indices around A and B matching c->a, a->d, d->b, b->c
                int e_bc = A.e[(kA + 1) % 3];
                int e_ca = A.e[(kA + 2) % 3];
                int e_ad = -1, e_db = -1;
                for (int k = 0; k < 3; ++k) {
                    int u = B.v[k], v = B.v[(k + 1) % 3];
                    if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                    if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
                }
                if (e_ad < 0 || e_db < 0) continue;

                mesh.edges[eIdx] = {c, d, makeStraightSeg(mesh.pts[c], mesh.pts[d], m_degree),
                                    CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false};
                mesh.tris[tA] = Tri{{c, a, d}, {e_ca, e_ad, eIdx}};
                mesh.tris[tB] = Tri{{c, d, b}, {eIdx, e_db, e_bc}};
                flipped = true;
                break;
            }
            if (!flipped) break;
        }
    };

    // Given a valid conforming Tri mesh (quads empty), test:
    //   1) direct 1-to-3 quadification
    //   2) triangle joining (bounded backtracking + greedy) -> OddStripFinish AND CatmullClarkMixedFinish
    auto evaluateTriangleJoining = [&](const MixedMesh& triMesh) {
        // Direct 1-to-3 tri split baseline
        {
            MixedMesh m1 = triMesh;
            std::vector<bool> allTrue(m1.edges.size(), true);
            SplittingCandidate c1;
            if (emitQuadPatches(m1, allTrue, c1)) updateBest(c1);
        }

        const int nT = static_cast<int>(triMesh.tris.size());
        if (nT < 2) return;

        const T radToDeg = T(180) / T(3.14159265358979323846);
        const T domMinA = getDomainMinCornerAngleDeg();
        const T qMinTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinA);
        const T qMaxTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

        struct PairEdge {
            int t0, t1;
            Quad q;
            T score;
        };

        struct EdgeAdj { int t0 = -1, k0 = -1, t1 = -1, k1 = -1; };
        std::vector<EdgeAdj> eAdj(triMesh.edges.size());
        for (int ti = 0; ti < nT; ++ti) {
            for (int k = 0; k < 3; ++k) {
                int eid = triMesh.tris[ti].e[k];
                if (eid >= 0 && eid < static_cast<int>(triMesh.edges.size())) {
                    if (eAdj[eid].t0 < 0) { eAdj[eid].t0 = ti; eAdj[eid].k0 = k; }
                    else if (eAdj[eid].t1 < 0) { eAdj[eid].t1 = ti; eAdj[eid].k1 = k; }
                }
            }
        }

        std::vector<PairEdge> pairEdges;
        for (int eIdx = 0; eIdx < static_cast<int>(triMesh.edges.size()); ++eIdx) {
            if (triMesh.edges[eIdx].isBoundary) continue;
            int tA = eAdj[eIdx].t0, kA = eAdj[eIdx].k0;
            int tB = eAdj[eIdx].t1, kB = eAdj[eIdx].k1;
            if (tA < 0 || tB < 0) continue;
            Tri A = triMesh.tris[tA], B = triMesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            gsPoint2D<T> qc[4] = {triMesh.pts[c], triMesh.pts[a], triMesh.pts[d], triMesh.pts[b]};
            T minArea2 = 0;
            if (!isStrictlyConvexQuad(qc, &minArea2)) continue;

            int e_bc = A.e[(kA + 1) % 3];
            int e_ca = A.e[(kA + 2) % 3];
            int e_ad = -1, e_db = -1;
            for (int k = 0; k < 3; ++k) {
                int u = B.v[k], v = B.v[(k + 1) % 3];
                if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
            }
            if (e_ad < 0 || e_db < 0) continue;

            int qv[4] = {c, a, d, b};
            int qe[4] = {e_ca, e_ad, e_db, e_bc};
            auto edgeDir = [&](int ei, int fromV, int toV, bool atEnd) -> gsPoint2D<T> {
                const auto& er = triMesh.edges[ei];
                if (er.isBoundary && er.seg.bsp) {
                    bool fwd = (er.a == fromV);
                    T u = (atEnd == fwd) ? T(1) : T(0);
                    gsPoint2D<T> dv = evalCurveDeriv(*er.seg.bsp, u);
                    T sgn = ((fwd && !atEnd) || (!fwd && atEnd)) ? T(1) : T(-1);
                    T nm = std::hypot(dv.x, dv.y);
                    if (nm > T(1e-14)) return gsPoint2D<T>(sgn * dv.x / nm, sgn * dv.y / nm);
                }
                T chSign = atEnd ? T(-1) : T(1);
                T dx = chSign * (triMesh.pts[toV].x - triMesh.pts[fromV].x);
                T dy = chSign * (triMesh.pts[toV].y - triMesh.pts[fromV].y);
                T nm = std::hypot(dx, dy);
                return (nm > T(1e-14)) ? gsPoint2D<T>(dx / nm, dy / nm) : gsPoint2D<T>(1, 0);
            };
            bool angOk = true;
            T minAng = T(360), maxAng = T(0), maxDev90 = T(0);
            T lenMin = T(1e30), lenMax = T(0);
            for (int k = 0; k < 4; ++k) {
                int kn = (k + 1) % 4, kp = (k + 3) % 4;
                gsPoint2D<T> tOut = edgeDir(qe[k],  qv[k], qv[kn], false);
                gsPoint2D<T> tBack = edgeDir(qe[kp], qv[kp], qv[k], true);
                T cr = tOut.x * tBack.y - tOut.y * tBack.x;
                T dt = tOut.x * tBack.x + tOut.y * tBack.y;
                T ang = std::atan2(cr, dt) * radToDeg;
                if (ang <= T(0)) ang += T(360);
                if (ang < qMinTol - T(1.5) || ang > qMaxTol) { angOk = false; break; }
                minAng = std::min(minAng, ang);
                maxAng = std::max(maxAng, ang);
                maxDev90 = std::max(maxDev90, std::abs(ang - T(90)));
                T eLen = triMesh.pts[qv[k]].dist(triMesh.pts[qv[kn]]);
                lenMin = std::min(lenMin, eLen);
                lenMax = std::max(lenMax, eLen);
            }
            if (!angOk) continue;

            T effMinQ = std::min(minAng, T(180) - maxAng);
            T arRatio = (lenMax > T(1e-14)) ? (lenMin / lenMax) : T(1);
            T sc = effMinQ + T(2.25) * arRatio - T(0.15) * maxDev90;
            pairEdges.push_back({tA, tB, Quad{{c, a, d, b}, {e_ca, e_ad, e_db, e_bc}}, sc});
        }
        if (pairEdges.empty()) return;

        std::sort(pairEdges.begin(), pairEdges.end(),
                  [](const PairEdge& p, const PairEdge& q) { return p.score > q.score; });

        std::vector<int> baseDeg(triMesh.pts.size(), 0);
        for (const auto& t : triMesh.tris) {
            ++baseDeg[t.v[0]]; ++baseDeg[t.v[1]]; ++baseDeg[t.v[2]];
        }

        // Collect distinct valid matchings (greedy + bounded backtracking)
        std::vector<std::vector<int>> candidateMatchings;
        auto pushUnique = [&](const std::vector<int>& m) {
            std::vector<int> ms = m;
            std::sort(ms.begin(), ms.end());
            for (const auto& prev : candidateMatchings)
                if (prev == ms) return;
            candidateMatchings.push_back(ms);
        };

        // 1. Greedy matching
        {
            std::vector<bool> used(nT, false);
            std::vector<int> deg = baseDeg;
            std::vector<int> gMatch;
            for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                const auto& pe = pairEdges[pi];
                if (used[pe.t0] || used[pe.t1]) continue;
                int a = pe.q.v[1], b = pe.q.v[3];
                if ((triMesh.isSmooth[a] && deg[a] <= 2) || (triMesh.isSmooth[b] && deg[b] <= 2))
                    continue;
                used[pe.t0] = used[pe.t1] = true;
                --deg[a]; --deg[b];
                gMatch.push_back(pi);
            }
            if (!gMatch.empty()) pushUnique(gMatch);
        }

        // 2. Bounded backtracking matching up to nT <= 26 to maximize pairs and quality
        if (nT <= 26) {
            std::vector<std::vector<int>> triPairs(nT);
            for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                triPairs[pairEdges[pi].t0].push_back(pi);
                triPairs[pairEdges[pi].t1].push_back(pi);
            }
            std::vector<bool> used(nT, false);
            std::vector<int> deg = baseDeg;
            std::vector<int> cur, bestMaxCard, bestMaxScore;
            int maxCard = 0;
            T maxScoreVal = -T(1);
            int nodes = 0;

            std::function<void(int, T)> dfs = [&](int ti, T curScore) {
                if (++nodes > 4000) return;
                while (ti < nT && used[ti]) ++ti;
                if (ti == nT) {
                    if (static_cast<int>(cur.size()) > maxCard) {
                        maxCard = static_cast<int>(cur.size());
                        bestMaxCard = cur;
                    }
                    T totalVal = T(cur.size()) * T(10) + curScore;
                    if (totalVal > maxScoreVal) {
                        maxScoreVal = totalVal;
                        bestMaxScore = cur;
                    }
                    return;
                }
                int remTris = 0;
                for (int j = ti; j < nT; ++j) if (!used[j]) ++remTris;
                if (static_cast<int>(cur.size()) + remTris / 2 < maxCard) return;

                // Try matching ti with each partner
                for (int pi : triPairs[ti]) {
                    const auto& pe = pairEdges[pi];
                    int other = (pe.t0 == ti) ? pe.t1 : pe.t0;
                    if (used[other]) continue;
                    int a = pe.q.v[1], b = pe.q.v[3];
                    if ((triMesh.isSmooth[a] && deg[a] <= 2) || (triMesh.isSmooth[b] && deg[b] <= 2))
                        continue;
                    used[ti] = used[other] = true;
                    --deg[a]; --deg[b];
                    cur.push_back(pi);
                    dfs(ti + 1, curScore + pe.score);
                    cur.pop_back();
                    ++deg[a]; ++deg[b];
                    used[ti] = used[other] = false;
                    if (nodes > 4000) return;
                }
                // Option: leave ti unmatched
                used[ti] = true;
                dfs(ti + 1, curScore);
                used[ti] = false;
            };
            dfs(0, T(0));
            if (!bestMaxCard.empty()) pushUnique(bestMaxCard);
            if (!bestMaxScore.empty()) pushUnique(bestMaxScore);
        }

        // Evaluate OddStripFinish and CatmullClarkMixedFinish on each matching
        for (const auto& m : candidateMatchings) {
            MixedMesh joined;
            joined.pts = triMesh.pts;
            joined.isSmooth = triMesh.isSmooth;
            joined.edges = triMesh.edges;
            std::vector<bool> usedTri(nT, false);
            for (int pi : m) {
                usedTri[pairEdges[pi].t0] = true;
                usedTri[pairEdges[pi].t1] = true;
                joined.quads.push_back(pairEdges[pi].q);
            }
            for (int ti = 0; ti < nT; ++ti) {
                if (!usedTri[ti]) joined.tris.push_back(triMesh.tris[ti]);
            }

            // Finish A: Opposite-edge parity propagation (minimal quad split)
            {
                MixedMesh jA = joined;
                std::vector<bool> needSplit;
                propagateOppositeSplits(jA, needSplit);
                SplittingCandidate cA;
                if (emitQuadPatches(jA, needSplit, cA)) updateBest(cA);
            }
            // Finish B: Uniform 1-edge Catmull-Clark mixed split fallback
            {
                MixedMesh jB = joined;
                std::vector<bool> allTrue(jB.edges.size(), true);
                SplittingCandidate cB;
                if (emitQuadPatches(jB, allTrue, cB)) updateBest(cB);
            }
        }
    };

    // Build initial curved boundary segs (subdividing K<=2 loops to >= 4 segments)
    std::vector<CurveSeg<T>> baseSegs;
    if (K == 1) {
        auto sH = CurveSeg<T>(outer.curve(0), m_degree).splitArc(T(0.5));
        auto s1 = sH.first.splitArc(T(0.5));
        auto s2 = sH.second.splitArc(T(0.5));
        baseSegs = {s1.first, s1.second, s2.first, s2.second};
    } else if (K == 2) {
        auto s0 = CurveSeg<T>(outer.curve(0), m_degree).splitArc(T(0.5));
        auto s1 = CurveSeg<T>(outer.curve(1), m_degree).splitArc(T(0.5));
        baseSegs = {s0.first, s0.second, s1.first, s1.second};
    } else {
        baseSegs.reserve(K);
        for (int i = 0; i < K; ++i) baseSegs.emplace_back(outer.curve(i), m_degree);
    }

    auto computeSmoothFlags = [&](const std::vector<CurveSeg<T>>& segs) -> std::vector<bool> {
        int M = static_cast<int>(segs.size());
        std::vector<bool> flags(M, false);
        const T radToDeg = T(180) / T(3.14159265358979323846);
        for (int i = 0; i < M; ++i) {
            const auto& prevBsp = *segs[(i + M - 1) % M].bsp;
            const auto& nextBsp = *segs[i].bsp;
            auto din = evalCurveDeriv(prevBsp, T(1));
            auto dout = evalCurveDeriv(nextBsp, T(0));
            T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
            T cr = (n1 > 1e-14 && n2 > 1e-14) ? (din.x * dout.y - din.y * dout.x) / (n1 * n2) : T(0);
            T dt = (n1 > 1e-14 && n2 > 1e-14) ? (din.x * dout.x + din.y * dout.y) / (n1 * n2) : T(1);
            T turnRad = std::atan2(cr, dt);
            T intAngDeg = T(180) - turnRad * radToDeg;
            // Mark smooth (~180 deg) OR obtuse (> 140 deg) vertices so Rivara bisection splits them
            flags[i] = (std::abs(turnRad) <= T(0.10)) || (intAngDeg > T(140));
        }
        return flags;
    };

    // Curved-safe ear clipping supporting sharp-first protection and multiple scoring modes
    auto buildEarTriangulation = [&](int mode, MixedMesh& outMesh) -> bool {
        struct RingNode {
            int ptIdx;
            CurveSeg<T> outSeg; // curve along active boundary leaving this node
        };
        std::vector<RingNode> ring;
        std::vector<gsPoint2D<T>> pts;
        std::vector<bool> isSmooth = computeSmoothFlags(baseSegs);
        std::vector<EdgeRec> edges;
        std::map<std::pair<int, int>, int> edgeMap;

        for (size_t i = 0; i < baseSegs.size(); ++i) {
            int idx = static_cast<int>(i);
            pts.push_back(baseSegs[i].p0());
            ring.push_back({idx, baseSegs[i]});
        }
        for (size_t i = 0; i < baseSegs.size(); ++i) {
            int a = static_cast<int>(i);
            int b = static_cast<int>((i + 1) % baseSegs.size());
            int id = static_cast<int>(edges.size());
            edges.push_back({a, b, baseSegs[i], CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true});
            edgeMap[{std::min(a, b), std::max(a, b)}] = id;
        }

        auto getOrAddDiag = [&](int a, int b) -> int {
            int lo = std::min(a, b), hi = std::max(a, b);
            auto it = edgeMap.find({lo, hi});
            if (it != edgeMap.end()) return it->second;
            int id = static_cast<int>(edges.size());
            edges.push_back({lo, hi, makeStraightSeg(pts[lo], pts[hi], m_degree),
                             CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false});
            edgeMap[{lo, hi}] = id;
            return id;
        };

        std::vector<int> triCount(pts.size(), 0);
        std::vector<Tri> tris;
        int subDivs = 0;

        auto pointInTri = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B,
                              const gsPoint2D<T>& C, const gsPoint2D<T>& P) -> bool {
            return cross2(A, B, P) >= -1e-9 && cross2(B, C, P) >= -1e-9 && cross2(C, A, P) >= -1e-9;
        };

        while (ring.size() > 3) {
            int m = static_cast<int>(ring.size());
            int bestIdx = -1;
            T bestScore = -T(1e30);

            // Two passes: pass 0 avoids leaving any smooth vertex with final deg == 1
            for (int pass = 0; pass < 2 && bestIdx < 0; ++pass) {
                for (int idx = 0; idx < m; ++idx) {
                    int iPrev = (idx + m - 1) % m;
                    int iNext = (idx + 1) % m;
                    int ia = ring[iPrev].ptIdx;
                    int ib = ring[idx].ptIdx;
                    int ic = ring[iNext].ptIdx;

                    if (pass == 0 && isSmooth[ib] && triCount[ib] == 0) continue;

                    const auto& A = pts[ia];
                    const auto& B = pts[ib];
                    const auto& C = pts[ic];
                    T area2 = cross2(A, B, C);
                    if (area2 <= T(1e-8)) continue;

                    if (!isSegmentInsidePoly(A, C, outerPoly)) continue;

                    bool clear = true;
                    for (int k = 0; k < m && clear; ++k) {
                        int ik = ring[k].ptIdx;
                        if (ik == ia || ik == ib || ik == ic) continue;
                        if (pointInTri(A, B, C, pts[ik])) clear = false;
                    }
                    // Also check interior samples of other active ring curves against curved bulge
                    for (int k = 0; k < m && clear; ++k) {
                        if (k == iPrev || k == idx) continue;
                        const auto& bsp = *ring[k].outSeg.bsp;
                        gsMatrix<T> ev;
                        gsMatrix<T> par(1, 3);
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
                    T dAC2 = std::pow(A.x - C.x, 2) + std::pow(A.y - C.y, 2);
                    if (mode == 0) {
                        T l2 = std::pow(A.x - B.x, 2) + std::pow(A.y - B.y, 2)
                             + std::pow(B.x - C.x, 2) + std::pow(B.y - C.y, 2) + dAC2;
                        sc = area2 / (l2 + T(1e-14));
                    } else if (mode == 1) {
                        sc = T(1) / (dAC2 + T(1e-14));
                    } else {
                        T lAB = std::hypot(A.x - B.x, A.y - B.y);
                        T lBC = std::hypot(B.x - C.x, B.y - C.y);
                        sc = area2 / ((lAB * lBC) * (dAC2 + T(1e-14)) + T(1e-14));
                    }
                    if (isSmooth[ib]) sc *= T(0.2);
                    if (sc > bestScore) {
                        bestScore = sc;
                        bestIdx = idx;
                    }
                }
            }

            if (bestIdx < 0) {
                // Subdivide active ring curve with largest chord sagitta (or longest curved edge)
                if (subDivs >= 24) return false;
                int splitK = -1;
                T maxSag = T(0);
                int longestK = -1;
                T maxLen = T(0);
                for (int k = 0; k < m; ++k) {
                    int kNext = (k + 1) % m;
                    auto sp = ring[k].outSeg.splitArc(T(0.5));
                    gsPoint2D<T> mid = sp.first.p1();
                    const auto& P = pts[ring[k].ptIdx];
                    const auto& Q = pts[ring[kNext].ptIdx];
                    T chord = std::hypot(Q.x - P.x, Q.y - P.y);
                    if (chord > maxLen) { maxLen = chord; longestK = k; }
                    T sag = (chord > 1e-12) ? std::abs(cross2(P, Q, mid)) / (chord * chord + T(1e-12)) : T(0);
                    if (sag > maxSag) {
                        maxSag = sag;
                        splitK = k;
                    }
                }
                if (splitK < 0 || maxSag < T(1e-4)) splitK = longestK;
                if (splitK < 0) return false;

                int kNext = (splitK + 1) % m;
                int u = ring[splitK].ptIdx;
                int v = ring[kNext].ptIdx;
                int eOld = edgeMap[{std::min(u, v), std::max(u, v)}];
                auto sp = ring[splitK].outSeg.splitArc(T(0.5));
                gsPoint2D<T> midPt = sp.first.p1();

                int newIdx = static_cast<int>(pts.size());
                pts.push_back(midPt);
                isSmooth.push_back(true);
                triCount.push_back(0);

                // Keep eOld as u -> midPt, add midPt -> v
                bool oldWasBdr = edges[eOld].isBoundary;
                edges[eOld] = {u, newIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, oldWasBdr};
                edgeMap.erase({std::min(u, v), std::max(u, v)});
                edgeMap[{std::min(u, newIdx), std::max(u, newIdx)}] = eOld;
                int eSecond = static_cast<int>(edges.size());
                edges.push_back({newIdx, v, sp.second, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, oldWasBdr});
                edgeMap[{std::min(newIdx, v), std::max(newIdx, v)}] = eSecond;

                ring[splitK].outSeg = sp.first;
                ring.insert(ring.begin() + splitK + 1, {newIdx, sp.second});
                ++subDivs;
                continue;
            }

            int iPrev = (bestIdx + m - 1) % m;
            int iNext = (bestIdx + 1) % m;
            int ia = ring[iPrev].ptIdx;
            int ib = ring[bestIdx].ptIdx;
            int ic = ring[iNext].ptIdx;

            int e0 = edgeMap[{std::min(ia, ib), std::max(ia, ib)}];
            int e1 = edgeMap[{std::min(ib, ic), std::max(ib, ic)}];
            int e2 = getOrAddDiag(ia, ic);
            tris.push_back({{ia, ib, ic}, {e0, e1, e2}});
            ++triCount[ia];
            ++triCount[ib];
            ++triCount[ic];

            const auto& er2 = edges[e2];
            ring[iPrev].outSeg = (er2.a == ia && er2.b == ic) ? er2.seg : er2.seg.reversed();
            ring.erase(ring.begin() + bestIdx);
        }

        if (ring.size() != 3) return false;
        int ia = ring[0].ptIdx, ib = ring[1].ptIdx, ic = ring[2].ptIdx;
        if (cross2(pts[ia], pts[ib], pts[ic]) <= T(1e-8)) return false;
        int e0 = getOrAddDiag(ia, ib);
        int e1 = getOrAddDiag(ib, ic);
        int e2 = getOrAddDiag(ic, ia);
        tris.push_back({{ia, ib, ic}, {e0, e1, e2}});

        outMesh.pts = std::move(pts);
        outMesh.isSmooth = std::move(isSmooth);
        outMesh.edges = std::move(edges);
        outMesh.tris = std::move(tris);
        outMesh.quads.clear();
        return true;
    };

    // Special K==3 reflex triangle -> 2 triangles by cutting reflex corner to opposite midpoint
    if (K == 3) {
        std::vector<gsPoint2D<T>> polyPts = {baseSegs[0].p0(), baseSegs[1].p0(), baseSegs[2].p0()};
        for (int r = 0; r < 3; ++r) {
            int p = (r + 2) % 3, opp = (r + 1) % 3;
            if (cross2(polyPts[p], polyPts[r], polyPts[opp]) < -1e-8) {
                auto spOpp = baseSegs[opp].splitArc(T(0.5));
                MixedMesh mesh;
                mesh.pts = {polyPts[0], polyPts[1], polyPts[2], spOpp.first.p1()};
                mesh.isSmooth = computeSmoothFlags(baseSegs);
                mesh.isSmooth.push_back(false);
                auto addEdge = [&](int a, int b, const CurveSeg<T>& s, bool isBdr) -> int {
                    int id = static_cast<int>(mesh.edges.size());
                    mesh.edges.push_back({a, b, s, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, isBdr});
                    return id;
                };
                int e_r_opp = addEdge(r, opp, baseSegs[r], true);
                int e_opp_m = addEdge(opp, 3, spOpp.first, true);
                int e_m_r   = addEdge(3, r, makeStraightSeg(mesh.pts[3], mesh.pts[r], m_degree), false);
                int e_m_p   = addEdge(3, p, spOpp.second, true);
                int e_p_r   = addEdge(p, r, baseSegs[p], true);
                mesh.tris = {
                    {{r, opp, 3}, {e_r_opp, e_opp_m, e_m_r}},
                    {{r, 3, p},   {e_m_r,   e_m_p,   e_p_r}}
                };
                evaluateTriangleJoining(mesh);
                break;
            }
        }
    }

    // Strategy T1: General curved ear-clipping (3 priority orders) + Rivara smooth fix + Lawson flips + triangle joining
    for (int mode = 0; mode < 3; ++mode) {
        MixedMesh mesh;
        if (!buildEarTriangulation(mode, mesh)) continue;
        fixSmoothDegreeOneVertices(mesh);
        evaluateTriangleJoining(mesh);

        MixedMesh flipped = mesh;
        applyLawsonDelaunayFlips(flipped);
        evaluateTriangleJoining(flipped);
    }

    // Strategy T2: Centroid-Star Triangulation + triangle joining
    if (pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        int M = static_cast<int>(baseSegs.size());
        std::vector<gsPoint2D<T>> pts(M + 1);
        bool allVisible = true;
        for (int i = 0; i < M; ++i) {
            pts[i] = baseSegs[i].p0();
            if (!isSegmentInsidePoly(centroid, pts[i], outerPoly)) allVisible = false;
        }
        pts[M] = centroid;

        if (allVisible && M >= 3) {
            MixedMesh starMesh;
            starMesh.pts = pts;
            starMesh.isSmooth = computeSmoothFlags(baseSegs);
            starMesh.isSmooth.push_back(false);
            for (int i = 0; i < M; ++i) {
                starMesh.edges.push_back({i, (i + 1) % M, baseSegs[i],
                                          CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true});
            }
            for (int i = 0; i < M; ++i) {
                starMesh.edges.push_back({M, i, makeStraightSeg(centroid, pts[i], m_degree),
                                          CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false});
            }
            for (int i = 0; i < M; ++i) {
                int j = (i + 1) % M;
                starMesh.tris.push_back({{M, i, j}, {M + i, i, M + j}});
            }
            evaluateTriangleJoining(starMesh);
        }
    }

    if (!bestCand.isValid) {
        bestCand = computeMultiHoleTriangulationSplitting();
    }

    return bestCand;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeMultiHoleTriangulationSplitting()
{
    using namespace internal_splitter;

    SplittingCandidate bestCand;
    bestCand.methodUsed = "MultiHoleTriangulation";
    bestCand.isValid = false;

    if (m_domain.numLoops() == 0) return bestCand;
    const gsMatrix<T> bbox = m_domain.boundingBox();
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    auto updateBest = [&](const SplittingCandidate& c) {
        if (!c.isValid) return;
        if (!bestCand.isValid) { bestCand = c; return; }
        T qNew = m_weightPatches * T(c.numPatches) + m_weightDistortion * c.distortionEnergy;
        T qOld = m_weightPatches * T(bestCand.numPatches) + m_weightDistortion * bestCand.distortionEnergy;
        if (qNew < qOld) bestCand = c;
    };

    std::function<void(TriMeshData<T>&, bool, bool)> quadifyTriMesh;
    quadifyTriMesh = [&](TriMeshData<T>& mesh, bool coarseRelaxed, bool allowParityRepair) {
        const int nT = static_cast<int>(mesh.tris.size());
        if (nT == 0) return;

        for (auto& er : mesh.edges) {
            if (!er.seg.bsp)
                er.seg = makeStraightSeg(mesh.pts[er.a], mesh.pts[er.b], m_degree);
            if (!er.splitDone) {
                auto sp = er.seg.splitArc(T(0.5));
                er.h0 = sp.first;
                er.h1 = sp.second;
                er.mid = sp.first.p1();
                er.splitDone = true;
                er.splitFrac = T(0.5);
            }
        }

        const T effMinValAng = coarseRelaxed ? std::max(T(8), minAngTol - T(3.5)) : minAngTol;
        const T effMaxValAng = coarseRelaxed ? std::max(maxAngTol, T(164.2))     : maxAngTol;
        const T effMinPairAng = std::max(T(8), minAngTol - T(4.0));
        const T effMaxPairAng = std::max(maxAngTol, T(164.2));

        auto emitFromMixed = [&](const std::vector<typename TriMeshData<T>::Tri>& tris,
                                 const std::vector<typename TriMeshData<T>::Quad>& quads,
                                 const std::vector<bool>& needSplit) {
            std::vector<typename TriMeshData<T>::EdgeRec> edges = mesh.edges;

            auto resplitEdge = [&](int eIdx, T paramU) {
                auto& er = edges[eIdx];
                T u = std::max(T(0.16), std::min(T(0.84), paramU));
                auto sp = er.seg.splitAtParam(u);
                er.h0 = sp.first;
                er.h1 = sp.second;
                er.mid = sp.first.p1();
                er.splitDone = true;
                er.splitFrac = u;
            };
            auto getSegDir = [&](int eIdx, int u, int v) -> CurveSeg<T> {
                const auto& er = edges[eIdx];
                return (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
            };
            auto getHalvesDir = [&](int eIdx, int u, int v, CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
                const auto& er = edges[eIdx];
                outMid = er.mid;
                if (er.a == u && er.b == v) { out0 = er.h0; out1 = er.h1; }
                else                         { out0 = er.h1.reversed(); out1 = er.h0.reversed(); }
            };

            for (const auto& t : tris)
                if (!needSplit[t.e[0]] || !needSplit[t.e[1]] || !needSplit[t.e[2]]) return;
            for (const auto& q : quads)
                if (needSplit[q.e[0]] != needSplit[q.e[2]] || needSplit[q.e[1]] != needSplit[q.e[3]]) return;

            auto chooseQuad4Center = [&](const typename TriMeshData<T>::Quad& q,
                                         const CurveSeg<T>& e0a, const CurveSeg<T>& e0b, const gsPoint2D<T>& m01,
                                         const CurveSeg<T>& e1a, const CurveSeg<T>& e1b, const gsPoint2D<T>& m12,
                                         const CurveSeg<T>& e2a, const CurveSeg<T>& e2b, const gsPoint2D<T>& m23,
                                         const CurveSeg<T>& e3a, const CurveSeg<T>& e3b, const gsPoint2D<T>& m30,
                                         gsPoint2D<T>& outO, T& outEffMin) -> bool {
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                gsPoint2D<T> O_m((m01.x + m12.x + m23.x + m30.x) / T(4), (m01.y + m12.y + m23.y + m30.y) / T(4));
                gsPoint2D<T> O_v((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x + mesh.pts[v3].x) / T(4),
                                 (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y + mesh.pts[v3].y) / T(4));
                std::vector<gsPoint2D<T>> oCands = {O_m, O_v};
                for (T w : {T(0.25), T(0.5), T(0.75)})
                    oCands.emplace_back((T(1) - w) * O_m.x + w * O_v.x, (T(1) - w) * O_m.y + w * O_v.y);
                for (const auto& vp : {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]})
                    oCands.emplace_back(T(0.75) * O_m.x + T(0.25) * vp.x, T(0.75) * O_m.y + T(0.25) * vp.y);

                bool found = false;
                T bestEff = -T(1e9);
                for (const auto& Oc : oCands) {
                    gsPoint2D<T> q0[4] = {mesh.pts[v0], m01, Oc, m30};
                    gsPoint2D<T> q1[4] = {mesh.pts[v1], m12, Oc, m01};
                    gsPoint2D<T> q2[4] = {mesh.pts[v2], m23, Oc, m12};
                    gsPoint2D<T> q3[4] = {mesh.pts[v3], m30, Oc, m23};
                    if (!isStrictlyConvexQuad(q0) || !isStrictlyConvexQuad(q1) ||
                        !isStrictlyConvexQuad(q2) || !isStrictlyConvexQuad(q3))
                        continue;
                    auto s01 = makeStraightSeg(Oc, m01, m_degree);
                    auto s12 = makeStraightSeg(Oc, m12, m_degree);
                    auto s23 = makeStraightSeg(Oc, m23, m_degree);
                    auto s30 = makeStraightSeg(Oc, m30, m_degree);
                    T mnA, mxA, eff0, eff1, eff2, eff3;
                    eval4SegCornerAnglesDeg(e0a, s01.reversed(), s30, e3b, mnA, mxA, eff0);
                    eval4SegCornerAnglesDeg(e1a, s12.reversed(), s01, e0b, mnA, mxA, eff1);
                    eval4SegCornerAnglesDeg(e2a, s23.reversed(), s12, e1b, mnA, mxA, eff2);
                    eval4SegCornerAnglesDeg(e3a, s30.reversed(), s23, e2b, mnA, mxA, eff3);
                    T effMin = std::min(std::min(eff0, eff1), std::min(eff2, eff3));
                    if (!found || effMin > bestEff) {
                        bestEff = effMin;
                        outO = Oc;
                        outEffMin = effMin;
                        found = true;
                    }
                }
                return found;
            };

            auto evalFaceSplitQuality = [&](bool isTri, int fIdx, T& minAngOut, T& maxAngOut, T& effMinOut) -> bool {
                minAngOut = T(360); maxAngOut = T(0); effMinOut = -T(1);
                if (isTri) {
                    const auto& t = tris[fIdx];
                    int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                    CurveSeg<T> e01a, e01b, e12a, e12b, e20a, e20b;
                    gsPoint2D<T> m01, m12, m20, O;
                    getHalvesDir(t.e[0], v0, v1, e01a, e01b, m01);
                    getHalvesDir(t.e[1], v1, v2, e12a, e12b, m12);
                    getHalvesDir(t.e[2], v2, v0, e20a, e20b, m20);
                    if (!internal_splitter::pickOptimalTri1to3CenterCurved(
                            mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                            e01a, e01b, e12a, e12b, e20a, e20b,
                            edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O))
                        return false;
                    auto s01 = makeStraightSeg(O, m01, m_degree);
                    auto s12 = makeStraightSeg(O, m12, m_degree);
                    auto s20 = makeStraightSeg(O, m20, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(e01a, s01.reversed(), s20, e20b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e12a, s12.reversed(), s01, e01b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e20a, s20.reversed(), s12, e12b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                }
                const auto& q = quads[fIdx];
                bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                if (!sp0 && !sp1) {
                    gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                    if (!isStrictlyConvexQuad(c)) return false;
                    eval4SegCornerAnglesDeg(getSegDir(q.e[0], v0, v1), getSegDir(q.e[1], v1, v2),
                                            getSegDir(q.e[2], v2, v3), getSegDir(q.e[3], v3, v0),
                                            minAngOut, maxAngOut, effMinOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                } else if (sp0 && !sp1) {
                    CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m01, m23;
                    getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                    getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                    gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                    gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                    if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return false;
                    auto sMid = makeStraightSeg(m01, m23, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(e0a, sMid, e2b, getSegDir(q.e[3], v3, v0), mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e0b, getSegDir(q.e[1], v1, v2), e2a, sMid.reversed(), mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                } else if (!sp0 && sp1) {
                    CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m12, m30;
                    getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                    getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                    gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                    gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                    if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return false;
                    auto sMid = makeStraightSeg(m12, m30, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(getSegDir(q.e[0], v0, v1), e1a, sMid, e3b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(sMid.reversed(), e1b, getSegDir(q.e[2], v2, v3), e3a, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                } else {
                    CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m01, m12, m23, m30, O;
                    getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                    getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                    getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                    getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                    T dummyEff = 0;
                    if (!chooseQuad4Center(q, e0a, e0b, m01, e1a, e1b, m12, e2a, e2b, m23, e3a, e3b, m30, O, dummyEff))
                        return false;
                    auto s01 = makeStraightSeg(O, m01, m_degree);
                    auto s12 = makeStraightSeg(O, m12, m_degree);
                    auto s23 = makeStraightSeg(O, m23, m_degree);
                    auto s30 = makeStraightSeg(O, m30, m_degree);
                    T mn, mx, ef;
                    eval4SegCornerAnglesDeg(e0a, s01.reversed(), s30, e3b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e1a, s12.reversed(), s01, e0b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e2a, s23.reversed(), s12, e1b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    eval4SegCornerAnglesDeg(e3a, s30.reversed(), s23, e2b, mn, mx, ef);
                    minAngOut = std::min(minAngOut, mn); maxAngOut = std::max(maxAngOut, mx);
                    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
                    return (minAngOut >= effMinValAng && maxAngOut <= effMaxValAng);
                }
            };

            // Build Coons multi-patch from current `edges`, recording per-face/per-edge patch ownership.
            // If `outBadEdges` is passed on failure, populates unique `needSplit` edges bordering bad patches.
            auto tryBuildFromCurrentEdges = [&](std::vector<int>* outBadEdges = nullptr) -> bool {
                gsMultiPatch<T> mp;
                std::vector<std::vector<int>> patchEdges;
                for (const auto& t : tris) {
                    int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                    CurveSeg<T> e01_a, e01_b, e12_a, e12_b, e20_a, e20_b;
                    gsPoint2D<T> m01, m12, m20;
                    getHalvesDir(t.e[0], v0, v1, e01_a, e01_b, m01);
                    getHalvesDir(t.e[1], v1, v2, e12_a, e12_b, m12);
                    getHalvesDir(t.e[2], v2, v0, e20_a, e20_b, m20);

                    gsPoint2D<T> O;
                    if (!internal_splitter::pickOptimalTri1to3CenterCurved(
                            mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                            e01_a, e01_b, e12_a, e12_b, e20_a, e20_b,
                            edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O)) {
                        if (outBadEdges) {
                            for (int k = 0; k < 3; ++k)
                                if (needSplit[t.e[k]]) outBadEdges->push_back(t.e[k]);
                        }
                        return false;
                    }
                    auto s01 = makeStraightSeg(O, m01, m_degree);
                    auto s12 = makeStraightSeg(O, m12, m_degree);
                    auto s20 = makeStraightSeg(O, m20, m_degree);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e01_a, s01.reversed(), s20, e20_b, m_degree)));
                    patchEdges.push_back({t.e[0], t.e[2]});
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e12_a, s12.reversed(), s01, e01_b, m_degree)));
                    patchEdges.push_back({t.e[1], t.e[0]});
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e20_a, s20.reversed(), s12, e12_b, m_degree)));
                    patchEdges.push_back({t.e[2], t.e[1]});
                }
                for (const auto& q : quads) {
                    bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                    int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                    auto pushQBad = [&]() {
                        if (outBadEdges)
                            for (int k = 0; k < 4; ++k)
                                if (needSplit[q.e[k]]) outBadEdges->push_back(q.e[k]);
                    };
                    if (!sp0 && !sp1) {
                        gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                        if (!isStrictlyConvexQuad(c)) { pushQBad(); return false; }
                        mp.addPatch(typename gsGeometry<T>::uPtr(
                            makeCoonsFromCCWCycle(getSegDir(q.e[0], v0, v1), getSegDir(q.e[1], v1, v2),
                                                  getSegDir(q.e[2], v2, v3), getSegDir(q.e[3], v3, v0), m_degree)));
                        patchEdges.push_back({q.e[0], q.e[1], q.e[2], q.e[3]});
                    } else if (sp0 && !sp1) {
                        CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m01, m23;
                        getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                        getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                        gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                        gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                        if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) { pushQBad(); return false; }
                        auto sMid = makeStraightSeg(m01, m23, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, sMid, e2b, getSegDir(q.e[3], v3, v0), m_degree)));
                        patchEdges.push_back({q.e[0], q.e[2]});
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0b, getSegDir(q.e[1], v1, v2), e2a, sMid.reversed(), m_degree)));
                        patchEdges.push_back({q.e[0], q.e[2]});
                    } else if (!sp0 && sp1) {
                        CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m12, m30;
                        getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                        getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                        gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                        gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                        if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) { pushQBad(); return false; }
                        auto sMid = makeStraightSeg(m12, m30, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(getSegDir(q.e[0], v0, v1), e1a, sMid, e3b, m_degree)));
                        patchEdges.push_back({q.e[1], q.e[3]});
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(sMid.reversed(), e1b, getSegDir(q.e[2], v2, v3), e3a, m_degree)));
                        patchEdges.push_back({q.e[1], q.e[3]});
                    } else {
                        CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m01, m12, m23, m30;
                        getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                        getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                        getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                        getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                        gsPoint2D<T> O;
                        T dummyEff = 0;
                        if (!chooseQuad4Center(q, e0a, e0b, m01, e1a, e1b, m12, e2a, e2b, m23, e3a, e3b, m30, O, dummyEff)) {
                            pushQBad();
                            return false;
                        }
                        auto s01 = makeStraightSeg(O, m01, m_degree);
                        auto s12 = makeStraightSeg(O, m12, m_degree);
                        auto s23 = makeStraightSeg(O, m23, m_degree);
                        auto s30 = makeStraightSeg(O, m30, m_degree);
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, s01.reversed(), s30, e3b, m_degree)));
                        patchEdges.push_back({q.e[0], q.e[3]});
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e1a, s12.reversed(), s01, e0b, m_degree)));
                        patchEdges.push_back({q.e[1], q.e[0]});
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e2a, s23.reversed(), s12, e1b, m_degree)));
                        patchEdges.push_back({q.e[2], q.e[1]});
                        mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e3a, s30.reversed(), s23, e2b, m_degree)));
                        patchEdges.push_back({q.e[3], q.e[2]});
                    }
                }
                if (mp.nPatches() == 0) return false;
                bool angPreOk = true;
                for (size_t p = 0; p < mp.nPatches(); ++p) {
                    T angs[4], mnA = 0, mxA = 180;
                    bool degen = false;
                    evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
                    if (degen || mnA < effMinValAng || mxA > effMaxValAng) {
                        angPreOk = false;
                        if (outBadEdges)
                            for (int eid : patchEdges[p])
                                if (needSplit[eid]) outBadEdges->push_back(eid);
                    }
                }
                if (!angPreOk) return false;
                gsMultiPatch<T> preWinslow = mp;
                optimizeWinslowBarrier(mp);
                mp.computeTopology(T(1e-5), false, false);
                std::string valErr;
                bool valOk = validateMultiPatch(m_domain, mp, T(1e-4), &valErr, effMinValAng, effMaxValAng);
                if (!valOk) {
                    preWinslow.computeTopology(T(1e-5), false, false);
                    if (validateMultiPatch(m_domain, preWinslow, T(1e-4), &valErr, effMinValAng, effMaxValAng)) {
                        mp = give(preWinslow);
                        valOk = true;
                    }
                }
                if (!valOk) {
                    if (outBadEdges) {
                        for (size_t p = 0; p < mp.nPatches(); ++p) {
                            T angs[4], mnA = 0, mxA = 180; bool degen = false;
                            evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
                            bool singleSign = isPatchDetSingleSign(mp.patch(p), 17, T(1e-7));
                            if (degen || mnA < effMinValAng || mxA > effMaxValAng || !singleSign) {
                                for (int eid : patchEdges[p])
                                    if (needSplit[eid]) outBadEdges->push_back(eid);
                            }
                        }
                    }
                    return false;
                }
                T minDet = 0, minA = 0, maxA = 180;
                T energy = evaluateDistortionEnergy(mp, &minDet, &minA, &maxA);
                SplittingCandidate cand;
                cand.multiPatch = give(mp);
                cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
                cand.minJacobianDet = minDet;
                cand.distortionEnergy = energy;
                cand.minCornerAngleDeg = minA;
                cand.maxCornerAngleDeg = maxA;
                cand.methodUsed = "MultiHoleTriangulation";
                cand.isValid = true;
                updateBest(cand);
                return true;
            };

            // Step 1: Fast baseline with untouched 0.5 midpoints
            std::vector<int> badEdges;
            if (tryBuildFromCurrentEdges(&badEdges)) return;
            if (badEdges.empty()) return;

            std::sort(badEdges.begin(), badEdges.end());
            badEdges.erase(std::unique(badEdges.begin(), badEdges.end()), badEdges.end());
            if (badEdges.size() > 24) badEdges.resize(24);

            // Step 2: Adaptive split-point rescue ONLY along edges bordering failing faces/patches
            struct IncFace { bool isTri; int idx; };
            std::vector<std::vector<IncFace>> edgeFaces(edges.size());
            for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti)
                for (int k = 0; k < 3; ++k)
                    edgeFaces[tris[ti].e[k]].push_back({true, ti});
            for (int qi = 0; qi < static_cast<int>(quads.size()); ++qi)
                for (int k = 0; k < 4; ++k)
                    edgeFaces[quads[qi].e[k]].push_back({false, qi});

            auto optimizeIncidentSplitEdge = [&](int eIdx) -> bool {
                const auto& inc = edgeFaces[eIdx];
                if (inc.empty()) return false;
                T curMinAll = T(1e9);
                bool curOkAll = true;
                for (const auto& f : inc) {
                    T mnA, mxA, efA;
                    bool okA = evalFaceSplitQuality(f.isTri, f.idx, mnA, mxA, efA);
                    curOkAll = curOkAll && okA;
                    curMinAll = std::min(curMinAll, okA ? efA : (efA - T(200)));
                }

                T oldFrac = edges[eIdx].splitFrac;
                auto oldH0 = edges[eIdx].h0;
                auto oldH1 = edges[eIdx].h1;
                auto oldMid = edges[eIdx].mid;

                std::vector<T> tCands = {arcLengthSplitParam(*edges[eIdx].seg.bsp, T(0.5))};
                for (const auto& f : inc) {
                    if (f.isTri) {
                        const auto& t = tris[f.idx];
                        for (int k = 0; k < 3; ++k)
                            if (t.e[k] == eIdx)
                                tCands.push_back(edges[eIdx].seg.closestParam(mesh.pts[t.v[(k + 2) % 3]]));
                    } else {
                        const auto& q = quads[f.idx];
                        for (int k = 0; k < 4; ++k)
                            if (q.e[k] == eIdx)
                                tCands.push_back(edges[eIdx].seg.closestParam(edges[q.e[(k + 2) % 4]].mid));
                    }
                }
                for (T tc : {T(0.24), T(0.34), T(0.42), T(0.58), T(0.66), T(0.76)})
                    tCands.push_back(tc);

                T bestScore = curMinAll;
                bool bestOk = curOkAll;
                T bestFrac = oldFrac;
                bool improved = false;
                for (T tTry : tCands) {
                    T uClamped = std::max(T(0.18), std::min(T(0.82), tTry));
                    resplitEdge(eIdx, uClamped);
                    T minAll = T(1e9);
                    bool okAll = true;
                    for (const auto& f : inc) {
                        T mnA, mxA, efA;
                        bool okA = evalFaceSplitQuality(f.isTri, f.idx, mnA, mxA, efA);
                        okAll = okAll && okA;
                        minAll = std::min(minAll, okA ? efA : (efA - T(200)));
                    }
                    // On curved boundary edges, also verify the actual Coons Jacobian of incident faces doesn't flip
                    if (okAll && edges[eIdx].isBoundary) {
                        for (const auto& f : inc) {
                            if (f.isTri) {
                                const auto& t = tris[f.idx];
                                int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                                CurveSeg<T> e01a, e01b, e12a, e12b, e20a, e20b; gsPoint2D<T> m01, m12, m20, O;
                                getHalvesDir(t.e[0], v0, v1, e01a, e01b, m01);
                                getHalvesDir(t.e[1], v1, v2, e12a, e12b, m12);
                                getHalvesDir(t.e[2], v2, v0, e20a, e20b, m20);
                                if (!internal_splitter::pickOptimalTri1to3CenterCurved(
                                        mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                                        e01a, e01b, e12a, e12b, e20a, e20b,
                                        edges[t.e[0]].isBoundary, edges[t.e[1]].isBoundary, edges[t.e[2]].isBoundary, O)) {
                                    okAll = false; break;
                                }
                                auto s01 = makeStraightSeg(O, m01, m_degree), s12 = makeStraightSeg(O, m12, m_degree), s20 = makeStraightSeg(O, m20, m_degree);
                                std::unique_ptr<gsTensorBSpline<2, T>> p0(makeCoonsFromCCWCycle(e01a, s01.reversed(), s20, e20b, m_degree));
                                std::unique_ptr<gsTensorBSpline<2, T>> p1(makeCoonsFromCCWCycle(e12a, s12.reversed(), s01, e01b, m_degree));
                                std::unique_ptr<gsTensorBSpline<2, T>> p2(makeCoonsFromCCWCycle(e20a, s20.reversed(), s12, e12b, m_degree));
                                if (!isPatchDetSingleSign(*p0, 11, T(1e-8)) || !isPatchDetSingleSign(*p1, 11, T(1e-8)) || !isPatchDetSingleSign(*p2, 11, T(1e-8))) {
                                    okAll = false; break;
                                }
                            }
                        }
                    }
                    bool accept = (!bestOk && okAll) ||
                                  (okAll == bestOk && minAll > bestScore + T(1.0));
                    if (accept) {
                        bestOk = okAll;
                        bestScore = minAll;
                        bestFrac = uClamped;
                        improved = true;
                    }
                }
                if (improved) {
                    resplitEdge(eIdx, bestFrac);
                    return true;
                }
                edges[eIdx].splitFrac = oldFrac;
                edges[eIdx].h0 = oldH0;
                edges[eIdx].h1 = oldH1;
                edges[eIdx].mid = oldMid;
                return false;
            };

            bool anyMoved = false;
            for (int eIdx : badEdges) {
                if (optimizeIncidentSplitEdge(eIdx)) anyMoved = true;
            }
            if (anyMoved) tryBuildFromCurrentEdges(nullptr);
        };

        // 1. Direct 1-to-3 split baseline
        {
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(mesh.tris, {}, allTrue);
        }

        // 2. Curved-aware max-min-angle triangle pairing + augmenting paths + local swap optimization
        struct PairEdge { int t0, t1; typename TriMeshData<T>::Quad q; T effMin; T score; };
        std::vector<PairEdge> pairEdges;
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
            int e_bc = A.e[(kA + 1) % 3], e_ca = A.e[(kA + 2) % 3];
            int e_ad = B.e[(kB + 1) % 3], e_db = B.e[(kB + 2) % 3];
            int qe[4] = {e_ca, e_ad, e_db, e_bc};
            T effMin = 0, sc = 0;
            if (!evalMeshQuadUnion(mesh, qv, qe, effMinPairAng, effMaxPairAng, effMin, sc))
                continue;
            pairEdges.push_back({tA, tB, {{c, a, d, b}, {e_ca, e_ad, e_db, e_bc}}, effMin, sc});
        }
        std::sort(pairEdges.begin(), pairEdges.end(),
                  [](const PairEdge& p, const PairEdge& q) { return p.score > q.score; });

        std::vector<int> initBaseDeg(mesh.pts.size(), 0);
        for (const auto& t : mesh.tris) { ++initBaseDeg[t.v[0]]; ++initBaseDeg[t.v[1]]; ++initBaseDeg[t.v[2]]; }
        std::vector<std::vector<int>> triPairsAdj(nT);
        for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
            triPairsAdj[pairEdges[pi].t0].push_back(pi);
            triPairsAdj[pairEdges[pi].t1].push_back(pi);
        }
        std::vector<T> triMinAng(nT, T(60));
        for (int ti = 0; ti < nT; ++ti) {
            T mxA = T(0);
            triMinMaxAnglesDeg(mesh.pts[mesh.tris[ti].v[0]],
                               mesh.pts[mesh.tris[ti].v[1]],
                               mesh.pts[mesh.tris[ti].v[2]], triMinAng[ti], mxA);
        }

        auto runMatchingStrategy = [&](T minEffThresh, std::vector<typename TriMeshData<T>::Quad>& joinedQuadsOut,
                                       std::vector<typename TriMeshData<T>::Tri>& remTrisOut) {
            std::vector<int> baseDeg = initBaseDeg;
            std::vector<int> matchedPair(nT, -1), partnerTri(nT, -1);

            for (int tier = 0; tier < 2; ++tier) {
                T thresh = (tier == 0) ? minEffThresh : T(0);
                if (tier == 1 && minEffThresh <= T(0)) break;
                for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                    const auto& pe = pairEdges[pi];
                    if (pe.effMin < thresh) continue;
                    if (matchedPair[pe.t0] >= 0 || matchedPair[pe.t1] >= 0) continue;
                    int a = pe.q.v[1], b = pe.q.v[3];
                    if ((mesh.isSmoothBnd[a] && baseDeg[a] <= 2) || (mesh.isSmoothBnd[b] && baseDeg[b] <= 2))
                        continue;
                    matchedPair[pe.t0] = pi; matchedPair[pe.t1] = pi;
                    partnerTri[pe.t0] = pe.t1; partnerTri[pe.t1] = pe.t0;
                    --baseDeg[a]; --baseDeg[b];
                }
            }

            // Augmenting paths prioritizing sharpest sliver unmatched triangles first
            std::vector<int> orderUnmatched;
            for (int ti = 0; ti < nT; ++ti)
                if (matchedPair[ti] < 0) orderUnmatched.push_back(ti);
            std::sort(orderUnmatched.begin(), orderUnmatched.end(),
                      [&](int a, int b) { return triMinAng[a] < triMinAng[b]; });

            for (int startTi : orderUnmatched) {
                if (matchedPair[startTi] >= 0) continue;
                std::vector<int> prevTri(nT, -1), prevPair(nT, -1);
                std::queue<int> qBfs;
                qBfs.push(startTi);
                prevTri[startTi] = startTi;
                int foundEnd = -1;
                while (!qBfs.empty() && foundEnd < 0) {
                    int u = qBfs.front(); qBfs.pop();
                    for (int pi : triPairsAdj[u]) {
                        const auto& pe = pairEdges[pi];
                        int v = (pe.t0 == u) ? pe.t1 : pe.t0;
                        if (v == partnerTri[u]) continue;
                        int ca = pe.q.v[1], cb = pe.q.v[3];
                        if ((mesh.isSmoothBnd[ca] && baseDeg[ca] <= 2) || (mesh.isSmoothBnd[cb] && baseDeg[cb] <= 2))
                            continue;
                        if (matchedPair[v] < 0) {
                            prevTri[v] = u; prevPair[v] = pi; foundEnd = v; break;
                        }
                        int w = partnerTri[v];
                        if (w >= 0 && prevTri[w] < 0) {
                            prevTri[v] = u; prevPair[v] = pi; prevTri[w] = v; qBfs.push(w);
                        }
                    }
                }
                if (foundEnd >= 0) {
                    int cur = foundEnd;
                    while (cur != startTi) {
                        int u = prevTri[cur], addPi = prevPair[cur];
                        int nextFree = (u == startTi) ? startTi : prevTri[u];
                        if (u != startTi) {
                            const auto& remPe = pairEdges[matchedPair[u]];
                            ++baseDeg[remPe.q.v[1]]; ++baseDeg[remPe.q.v[3]];
                        }
                        const auto& addPe = pairEdges[addPi];
                        matchedPair[u] = addPi; matchedPair[cur] = addPi;
                        partnerTri[u] = cur; partnerTri[cur] = u;
                        --baseDeg[addPe.q.v[1]]; --baseDeg[addPe.q.v[3]];
                        cur = nextFree;
                    }
                }
            }

            // Local 1-edge shift improvement: replace (t0, t1) with higher-angle (t1, freeT2)
            for (int sweep = 0; sweep < 2; ++sweep) {
                bool changed = false;
                for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
                    const auto& newPe = pairEdges[pi];
                    int tA = newPe.t0, tB = newPe.t1;
                    bool aFree = (matchedPair[tA] < 0), bFree = (matchedPair[tB] < 0);
                    if (aFree == bFree) continue;
                    int freeT = aFree ? tA : tB;
                    int occT  = aFree ? tB : tA;
                    int oldPi = matchedPair[occT];
                    const auto& oldPe = pairEdges[oldPi];
                    if (newPe.effMin <= oldPe.effMin + T(3.0)) continue;
                    int occMate = partnerTri[occT];
                    int ca = newPe.q.v[1], cb = newPe.q.v[3];
                    bool sameVerts = (ca == oldPe.q.v[1] && cb == oldPe.q.v[3]) || (ca == oldPe.q.v[3] && cb == oldPe.q.v[1]);
                    if (!sameVerts && ((mesh.isSmoothBnd[ca] && baseDeg[ca] <= 2) || (mesh.isSmoothBnd[cb] && baseDeg[cb] <= 2)))
                        continue;
                    ++baseDeg[oldPe.q.v[1]]; ++baseDeg[oldPe.q.v[3]];
                    matchedPair[occMate] = -1; partnerTri[occMate] = -1;
                    matchedPair[freeT] = pi; matchedPair[occT] = pi;
                    partnerTri[freeT] = occT; partnerTri[occT] = freeT;
                    --baseDeg[ca]; --baseDeg[cb];
                    changed = true;
                }
                if (!changed) break;
            }

            joinedQuadsOut.clear();
            remTrisOut.clear();
            for (int ti = 0; ti < nT; ++ti) {
                if (matchedPair[ti] >= 0) {
                    if (partnerTri[ti] > ti) joinedQuadsOut.push_back(pairEdges[matchedPair[ti]].q);
                } else {
                    remTrisOut.push_back(mesh.tris[ti]);
                }
            }
        };

        auto propagateAndEmit = [&](const std::vector<typename TriMeshData<T>::Quad>& joinedQuads,
                                    const std::vector<typename TriMeshData<T>::Tri>& remTris) -> bool {
            if (joinedQuads.empty()) return false;
            std::vector<bool> needSplit(mesh.edges.size(), false);
            for (const auto& t : remTris)
                needSplit[t.e[0]] = needSplit[t.e[1]] = needSplit[t.e[2]] = true;
            bool changed = true;
            while (changed) {
                changed = false;
                for (const auto& q : joinedQuads) {
                    if (needSplit[q.e[0]] != needSplit[q.e[2]]) {
                        needSplit[q.e[0]] = needSplit[q.e[2]] = true;
                        changed = true;
                    }
                    if (needSplit[q.e[1]] != needSplit[q.e[3]]) {
                        needSplit[q.e[1]] = needSplit[q.e[3]] = true;
                        changed = true;
                    }
                }
            }
            int prevP = bestCand.isValid ? bestCand.numPatches : 999999;
            emitFromMixed(remTris, joinedQuads, needSplit);
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(remTris, joinedQuads, allTrue);
            return bestCand.isValid && bestCand.numPatches < prevP;
        };

        std::vector<typename TriMeshData<T>::Quad> jQuads0, jQuads1;
        std::vector<typename TriMeshData<T>::Tri> rem0, rem1;
        if (!pairEdges.empty()) {
            runMatchingStrategy(T(0), jQuads0, rem0);
            propagateAndEmit(jQuads0, rem0);

            T highThresh = std::max(effMinPairAng + T(8), T(22));
            runMatchingStrategy(highThresh, jQuads1, rem1);
            if (jQuads1.size() != jQuads0.size() || !rem1.empty())
                propagateAndEmit(jQuads1, rem1);

            // Parity repair: if 1-2 unmatched triangles trigger global cascades on a coarse mesh,
            // bisect one boundary edge of the worst unmatched triangle and re-CDT + re-pair.
            if (allowParityRepair && coarseRelaxed && nT <= 90 && !rem0.empty() && rem0.size() <= 3) {
                int candTri = -1, candEdge = -1, candK = -1;
                T bestLen = -T(1);
                for (const auto& rt : rem0) {
                    for (int k = 0; k < 3; ++k) {
                        int eid = rt.e[k];
                        if (!mesh.edges[eid].isBoundary) continue;
                        T el = mesh.pts[rt.v[k]].dist(mesh.pts[rt.v[(k + 1) % 3]]);
                        if (el > bestLen) {
                            bestLen = el; candEdge = eid; candK = k;
                            for (int ti = 0; ti < nT; ++ti)
                                if (mesh.tris[ti].e[0] == rt.e[0] && mesh.tris[ti].e[1] == rt.e[1] && mesh.tris[ti].e[2] == rt.e[2])
                                    candTri = ti;
                        }
                    }
                }
                if (candTri >= 0 && candEdge >= 0) {
                    const auto bt = mesh.tris[candTri];
                    int u = bt.v[candK], v = bt.v[(candK + 1) % 3], opp = bt.v[(candK + 2) % 3];
                    int e_vo = bt.e[(candK + 1) % 3], e_ou = bt.e[(candK + 2) % 3];
                    auto oldEr = mesh.edges[candEdge];
                    T tMid = arcLengthSplitParam(*oldEr.seg.bsp, T(0.5));
                    T tPrj = std::max(T(0.22), std::min(T(0.78), oldEr.seg.closestParam(mesh.pts[opp])));
                    for (T tTry : {tPrj, tMid}) {
                        TriMeshData<T> rMesh = mesh;
                        auto sp = oldEr.seg.splitAtParam(tTry);
                        gsPoint2D<T> sPt = sp.first.p1();
                        if (cross2(rMesh.pts[u], sPt, rMesh.pts[opp]) <= T(1e-9) ||
                            cross2(sPt, rMesh.pts[v], rMesh.pts[opp]) <= T(1e-9))
                            continue;
                        if (!isSegmentInsideDomain(rMesh.pts[opp], sPt, rMesh.loopPolys)) continue;
                        int mIdx = static_cast<int>(rMesh.pts.size());
                        rMesh.pts.push_back(sPt);
                        rMesh.isBoundaryPt.push_back(true);
                        rMesh.isSmoothBnd.push_back(true);
                        rMesh.isSharpCorner.push_back(false);
                        rMesh.cornerAngleDeg.push_back(T(180));
                        auto din = evalCurveDeriv(*sp.first.bsp, T(1));
                        T dn = std::hypot(din.x, din.y);
                        gsPoint2D<T> tDir = (dn > T(1e-14)) ? gsPoint2D<T>(din.x / dn, din.y / dn) : gsPoint2D<T>(1, 0);
                        rMesh.bndTanIn.push_back(tDir); rMesh.bndTanOut.push_back(tDir);
                        int oldA = oldEr.a, oldB = oldEr.b;
                        rMesh.edges[candEdge] = {oldA, mIdx, sp.first, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, oldEr.loopIdx};
                        int e_mb = static_cast<int>(rMesh.edges.size());
                        rMesh.edges.push_back({mIdx, oldB, sp.second, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, true, oldEr.loopIdx});
                        int e_om = static_cast<int>(rMesh.edges.size());
                        rMesh.edges.push_back({opp, mIdx, makeStraightSeg(rMesh.pts[opp], sPt, m_degree),
                                               CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false, false, -1});
                        auto halfFor = [&](int ep) { return (ep == oldA) ? candEdge : e_mb; };
                        rMesh.tris[candTri] = {{u, mIdx, opp}, {halfFor(u), e_om, e_ou}};
                        rMesh.tris.push_back({{mIdx, v, opp}, {halfFor(v), e_vo, e_om}});
                        applyDomainLawsonFlips(rMesh);
                        flipForBetterQuadUnions(rMesh, effMinPairAng, effMaxPairAng);
                        quadifyTriMesh(rMesh, coarseRelaxed, false);
                        break;
                    }
                }
            }
        }
        if (!bestCand.isValid ||
            (coarseRelaxed && nT <= 35 && bestCand.numPatches > 3 * nT && (jQuads0.empty() || nT <= 35)))
        {
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(mesh.tris, {}, allTrue);
        }
    };

    // Coarse-first passes: true CDT + feature projections + quad-union co-optimization + targeted sliver repair
    for (int coarsePass = 0; coarsePass < 5 && (!bestCand.isValid || (coarsePass == 0 && bestCand.numPatches > 29)); ++coarsePass) {
        const bool cRef = (coarsePass > 0);
        const bool cRefineHoles = (coarsePass >= 3);
        const T cH = (coarsePass <= 1) ? (diag * T(0.65))
                   : ((coarsePass == 2) ? (diag * T(0.45))
                   : ((coarsePass == 3) ? (diag * T(0.32)) : (diag * T(0.25))));
        const T cBend = (coarsePass <= 1) ? T(1.10) : T(0.65);
        const int cNeedleRounds = (coarsePass >= 3) ? 2 : 1;

        for (int featPass = 0; featPass < 2; ++featPass) {
            const bool useFeat = (featPass == 1);
            if (useFeat && bestCand.isValid && (coarsePass > 0 || bestCand.numPatches <= 29)) break;

            for (int mode = 0; mode < 3 && (!bestCand.isValid || (coarsePass == 0 && bestCand.numPatches > 29)); ++mode) {
                TriMeshData<T> coarse;
                if (!buildConstrainedDomainTriangulation(m_domain, m_degree, cH, cRef, mode, coarse, cBend, cRefineHoles, useFeat))
                    continue;
                if (coarsePass == 0 && !useFeat)
                    quadifyTriMesh(coarse, true, true);

                TriMeshData<T> cdt = coarse;
                applyDomainLawsonFlips(cdt);
                if (!bestCand.isValid || (coarsePass <= 1 && cdt.tris.size() <= 45 && bestCand.numPatches > 29))
                    quadifyTriMesh(cdt, true, true);

                if (!bestCand.isValid || (coarsePass <= 1 && cdt.tris.size() <= 60 && bestCand.numPatches > 29)) {
                    TriMeshData<T> qOpt = cdt;
                    flipForBetterQuadUnions(qOpt, std::max(T(8), minAngTol - T(4.0)), std::max(maxAngTol, T(164.2)));
                    quadifyTriMesh(qOpt, true, true);
                }

                if (!bestCand.isValid) {
                    TriMeshData<T> mild = cdt;
                    refineNeedleCircumcenters(mild, cH * T(0.55), cNeedleRounds, T(16.0), true, 16);
                    bisectObtuseTriangles(mild, m_degree, T(146), 28);
                    applyDomainLawsonFlips(mild);
                    quadifyTriMesh(mild, true, false);
                }
            }
        }
    }

    for (int refinePass = 0; refinePass < 2 && !bestCand.isValid; ++refinePass) {
        const bool doRef = (refinePass == 1);
        const T hStep = doRef ? (diag / T(6)) : (diag / T(8));
        for (int mode = 0; mode < 3 && (!bestCand.isValid || !doRef); ++mode) {
            if (mode == 2 && bestCand.isValid) continue;
            TriMeshData<T> coarse;
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, hStep, doRef, mode, coarse, T(0.44), true, true))
                continue;
            TriMeshData<T> cdt = coarse;
            applyDomainLawsonFlips(cdt);
            quadifyTriMesh(cdt, false, false);

            if (!bestCand.isValid) {
                TriMeshData<T> qOpt = cdt;
                flipForBetterQuadUnions(qOpt, std::max(T(8), minAngTol - T(4.0)), std::max(maxAngTol, T(164.2)));
                quadifyTriMesh(qOpt, false, false);
            }

            if (!bestCand.isValid) {
                TriMeshData<T> steiner = cdt;
                refineInteriorTriangles(steiner, hStep, 2);
                bisectObtuseTriangles(steiner, m_degree, T(146), 28);
                applyDomainLawsonFlips(steiner);
                quadifyTriMesh(steiner, false, false);
            }
        }
    }

    if (!bestCand.isValid) {
        TriMeshData<T> fine;
        if (buildConstrainedDomainTriangulation(m_domain, m_degree, diag / T(10), true, 0, fine, T(0.44), true, true)) {
            applyDomainLawsonFlips(fine);
            refineInteriorTriangles(fine, diag / T(10), 2);
            bisectObtuseTriangles(fine, m_degree, T(146), 28);
            applyDomainLawsonFlips(fine);
            quadifyTriMesh(fine, false, false);
        }
    }

    // Fallback: Ruppert-style strict-interior circumcenter refinement for acute needles + coarse diag/4 bisect
    if (!bestCand.isValid) {
        for (int fPass = 0; fPass < 7 && !bestCand.isValid; ++fPass) {
            T fH = (fPass == 0 || fPass == 4 || fPass == 6) ? (diag / T(10))
                 : ((fPass == 1) ? (diag / T(5))
                 : ((fPass == 2 || fPass == 5) ? (diag / T(4)) : (diag / T(7))));
            bool fRef = (fPass != 2 && fPass != 5);
            int fMode = (fPass == 1) ? 1 : 0;
            int nRounds = (fPass == 2) ? 2 : ((fPass >= 3) ? 3 : 2);
            bool runNeedle = (fPass != 2);
            T needleAng = (fPass >= 3) ? T(18) : T(16);
            T maxQ = (fPass == 2) ? T(146) : ((fPass >= 3) ? T(142) : T(144));
            int maxSp = (fPass == 2) ? 28 : ((fPass >= 3) ? 56 : 48);
            T minQ = (fPass >= 3) ? T(16.5) : T(14.5);
            bool allowCentroid = (fPass >= 5);
            TriMeshData<T> needle;
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, fH, fRef, fMode, needle, T(0.44), true, false) &&
                !buildConstrainedDomainTriangulation(m_domain, m_degree, fH, fRef, fMode, needle, T(0.44), true, true))
                continue;
            applyDomainLawsonFlips(needle);
            refineInteriorTriangles(needle, fH, nRounds);
            if (runNeedle)
                refineNeedleCircumcenters(needle, fH, nRounds, needleAng, allowCentroid);
            bisectObtuseTriangles(needle, m_degree, maxQ, maxSp, minQ);
            applyDomainLawsonFlips(needle);
            quadifyTriMesh(needle, false, false);
        }
    }

    return bestCand;
}

// ============================================================================
// Approach 2: Medial Axis Transform (MAT) Skeleton Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeMedialAxisSplitting()
{
    SplittingCandidate candidate;
    candidate.methodUsed = "MedialAxisSkeleton";

    if (m_domain.numLoops() == 0) {
        candidate.isValid = false;
        return candidate;
    }

    const auto& outer = m_domain.outer();
    int nHoles = m_domain.numHoles();
    int nOuterCurves = outer.size();

    if (nHoles == 0) {
        return computeHoleFreeStructuralSplitting();
    }

    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T minHoleAngTol = std::min(minAngTol, T(4.0));
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    auto tryFinalizeHole = [&](SplittingCandidate& cand) -> bool {
        const bool dbg = (std::getenv("GISMO_SPLITTER_DEBUG") != nullptr);
        if (cand.multiPatch.nPatches() == 0) { if (dbg) std::cerr << "[DBG] empty multiPatch\n"; return false; }
        cand.multiPatch.computeTopology(T(1e-4), false, false);
        optimizeWinslowBarrier(cand.multiPatch);
        gsASG1LocalFitter<T>::unifyPatchOrientations(cand.multiPatch);
        std::string vErr;
        if (!validateMultiPatch(m_domain, cand.multiPatch, T(1e-4), &vErr, minHoleAngTol, std::max(maxAngTol, T(164.9)))) {
            if (dbg) {
                T md=0, ma=0, xa=180;
                evaluateDistortionEnergy(cand.multiPatch, &md, &ma, &xa);
                std::cerr << "[DBG] validateMultiPatch failed (P=" << cand.multiPatch.nPatches()
                          << ", IF=" << cand.multiPatch.topology().nInterfaces()
                          << ", BND=" << cand.multiPatch.topology().nBoundary()
                          << ", domMin=" << domMinAng << ", tol=[" << minHoleAngTol << "," << std::max(maxAngTol, T(164.9)) << "]"
                          << ", minDet=" << md << ", ang=[" << ma << "," << xa << "]): " << vErr << "\n";
            }
            cand.multiPatch.clear();
            cand.isValid = false;
            return false;
        }
        T minPos20 = std::numeric_limits<T>::max();
        for (size_t p = 0; p < cand.multiPatch.nPatches(); ++p) {
            T pMin = 0, pMax = 0, pAbs = 0, pEn = 0;
            internal_splitter::evalPatchJacobianAndWinslow(cand.multiPatch.patch(p), 20, pMin, pMax, pAbs, pEn);
            minPos20 = std::min(minPos20, pMin);
        }
        cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
        cand.distortionEnergy = evaluateDistortionEnergy(cand.multiPatch, &cand.minJacobianDet,
                                                         &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
        cand.isValid = (minPos20 > T(1e-8) && cand.minJacobianDet > T(1e-8) && cand.maxCornerAngleDeg < T(164.9));
        if (!cand.isValid) {
            if (dbg) std::cerr << "[DBG] post-check failed (P=" << cand.numPatches
                               << " minPos20=" << minPos20 << " minDet=" << cand.minJacobianDet
                               << " ang=[" << cand.minCornerAngleDeg << "," << cand.maxCornerAngleDeg << "])\n";
            cand.multiPatch.clear();
            return false;
        }
        return true;
    };

    using internal_splitter::CurveSeg;
    std::vector<CurveSeg<T>> oCCW = internal_splitter::orientLoopSegments(outer, m_degree, true);
    auto mkLine = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> CurveSeg<T> {
        std::unique_ptr<gsBSpline<T>> ptr(internal_splitter::makeStraightBSpline(a, b, m_degree));
        return CurveSeg<T>(*ptr);
    };
    auto findReflexIdx = [&](const std::vector<CurveSeg<T>>& loop, T minDeg = T(205)) -> std::vector<int> {
        const int M = static_cast<int>(loop.size());
        const T radToDeg = T(180) / T(3.14159265358979323846);
        std::vector<int> res;
        for (int i = 0; i < M; ++i) {
            auto d0 = internal_splitter::evalCurveDeriv(*loop[(i + M - 1) % M].bsp, T(1));
            auto d1 = internal_splitter::evalCurveDeriv(*loop[i].bsp, T(0));
            T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
            if (n0 <= T(1e-14) || n1 <= T(1e-14)) continue;
            T cr = (d0.x * d1.y - d0.y * d1.x) / (n0 * n1);
            T dt = (d0.x * d1.x + d0.y * d1.y) / (n0 * n1);
            T intAng = T(180) - std::atan2(cr, dt) * radToDeg;
            if (intAng > minDeg) res.push_back(i);
        }
        return res;
    };

    std::vector<std::vector<gsPoint2D<T>>> holePolys(nHoles);
    std::vector<gsPoint2D<T>> holeCtrs(nHoles);
    for (int h = 0; h < nHoles; ++h) {
        const auto& hl = m_domain.loop(h + 1);
        for (int c = 0; c < hl.size(); ++c) {
            gsMatrix<T> pr = hl.curve(c).parameterRange();
            for (int s = 0; s < 8; ++s) {
                T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(s) / T(8));
                holePolys[h].push_back(internal_splitter::evalCurvePoint(hl.curve(c), t));
            }
        }
        holeCtrs[h] = internal_splitter::polygonCentroid(holePolys[h]);
    }

    auto assignHolesToPolys = [&](const std::vector<std::vector<gsPoint2D<T>>>& cellPolys)
        -> std::vector<std::vector<const gsCurveLoop<T>*>>
    {
        std::vector<std::vector<const gsCurveLoop<T>*>> res(cellPolys.size());
        for (int h = 0; h < nHoles; ++h) {
            int bestC = -1;
            for (size_t c = 0; c < cellPolys.size(); ++c) {
                if (internal_splitter::pointInPolyStrict(holeCtrs[h], cellPolys[c], T(1e-9))) {
                    bestC = static_cast<int>(c);
                    break;
                }
            }
            if (bestC < 0) {
                T bestD = T(1e18);
                for (size_t c = 0; c < cellPolys.size(); ++c) {
                    gsPoint2D<T> cc = internal_splitter::polygonCentroid(cellPolys[c]);
                    T d = (cc.x - holeCtrs[h].x)*(cc.x - holeCtrs[h].x) + (cc.y - holeCtrs[h].y)*(cc.y - holeCtrs[h].y);
                    if (d < bestD) { bestD = d; bestC = static_cast<int>(c); }
                }
            }
            if (bestC >= 0) res[bestC].push_back(&m_domain.loop(h + 1));
        }
        return res;
    };

    auto addBlockCCW = [&](const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                           const CurveSeg<T>& e23, const CurveSeg<T>& e30,
                           const std::vector<const gsCurveLoop<T>*>& blockHoles) -> bool {
        return internal_splitter::build4SidedMultiHoleGrid(
            e01, e12, e23.reversed(), e30.reversed(), blockHoles, candidate.multiPatch, m_degree);
    };

    auto blockBoundaryPoly = [](const CurveSeg<T>& e01, const CurveSeg<T>& e12,
                                const CurveSeg<T>& e23, const CurveSeg<T>& e30) -> std::vector<gsPoint2D<T>> {
        std::vector<gsPoint2D<T>> poly;
        poly.reserve(32);
        const CurveSeg<T>* arr[4] = { &e01, &e12, &e23, &e30 };
        for (int s = 0; s < 4; ++s) {
            for (int k = 0; k < 8; ++k)
                poly.push_back(internal_splitter::evalCurvePoint(*arr[s]->bsp, static_cast<T>(k) / T(8)));
        }
        return poly;
    };

    // ------------------------------------------------------------------------
    // Case 2: 4 Outer Curves, Any Number of Holes (1, 2, 3, ..., N)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 4) {
        if (nHoles == 1 && m_domain.loop(1).size() == 4) {
            internal_splitter::buildQuadRingByProximity(outer, m_domain.loop(1), candidate.multiPatch, m_degree);
            if (tryFinalizeHole(candidate)) return candidate;
        }
        std::vector<const gsCurveLoop<T>*> allHoles;
        for (int h = 0; h < nHoles; ++h) allHoles.push_back(&m_domain.loop(h + 1));
        for (int rot = 0; rot < 2; ++rot) {
            candidate.multiPatch.clear();
            if (addBlockCCW(oCCW[rot], oCCW[(rot+1)%4], oCCW[(rot+2)%4], oCCW[(rot+3)%4], allHoles)) {
                if (tryFinalizeHole(candidate)) return candidate;
            }
        }
        if (nHoles >= 2) {
            candidate.multiPatch.clear();
            if (internal_splitter::build4SidedInnerAxisGrid(oCCW, allHoles, candidate.multiPatch, m_degree)) {
                if (tryFinalizeHole(candidate)) return candidate;
            }
        }
        candidate.multiPatch.clear();
    }

    // ------------------------------------------------------------------------
    // Case 5 & 6: 6-Corner L-Shape with 1, 3, or Arbitrary Holes
    // ------------------------------------------------------------------------
    if (nOuterCurves == 6) {
        std::vector<int> reflex = findReflexIdx(oCCW);
        if (reflex.size() == 1) {
            int r = reflex[0];
            CurveSeg<T> c[6];
            for (int k = 0; k < 6; ++k) c[k] = oCCW[(r + k) % 6];

            auto sp2 = c[2].splitArc(T(0.5));
            auto sp3 = c[3].splitArc(T(0.5));
            gsPoint2D<T> Vr = c[0].p0(), M2 = sp2.first.p1(), M3 = sp3.first.p1();
            CurveSeg<T> s2 = mkLine(Vr, M2), s3 = mkLine(Vr, M3);

            std::vector<std::vector<gsPoint2D<T>>> polys = {
                blockBoundaryPoly(c[0], c[1], sp2.first, s2.reversed()),
                blockBoundaryPoly(sp2.second, sp3.first, s3, s2),
                blockBoundaryPoly(sp3.second, c[4], c[5], s3.reversed())
            };
            auto assigned = assignHolesToPolys(polys);
            bool ok = true;
            ok = ok && addBlockCCW(c[0], c[1], sp2.first, s2.reversed(), assigned[0]);
            ok = ok && addBlockCCW(sp2.second, sp3.first, s3, s2, assigned[1]);
            ok = ok && addBlockCCW(sp3.second, c[4], c[5], s3.reversed(), assigned[2]);
            if (ok && tryFinalizeHole(candidate)) return candidate;
            candidate.multiPatch.clear();
        }
    }

    // ------------------------------------------------------------------------
    // Case 7: 2-Corner Crescent with Holes (crescent_2corner_3holes, etc.)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 2 && nHoles >= 1) {
        // Orient cBot left->right, cTop right->left so tips are P_L and P_R
        CurveSeg<T> c0 = oCCW[0], c1 = oCCW[1];
        gsPoint2D<T> P_L = c0.p0(), P_R = c0.p1();
        if (P_L.x > P_R.x) { c0 = oCCW[1]; c1 = oCCW[0]; P_L = c0.p0(); P_R = c0.p1(); }
        // Sort holes by centroid.x ascending
        std::vector<int> hOrder(nHoles);
        std::iota(hOrder.begin(), hOrder.end(), 0);
        std::sort(hOrder.begin(), hOrder.end(), [&](int a, int b) { return holeCtrs[a].x < holeCtrs[b].x; });

        if (nHoles == 3) {
            // Find parameter split cuts separating left / mid / right holes cleanly
            auto findCutParam = [&](const CurveSeg<T>& seg, T targetX) -> T {
                T bestU = T(0.5), bestD = T(1e9);
                for (int k = 4; k <= 60; ++k) {
                    T u = static_cast<T>(k) / T(64);
                    auto p = internal_splitter::evalCurvePoint(*seg.bsp, u);
                    T d = std::abs(p.x - targetX);
                    if (d < bestD) { bestD = d; bestU = u; }
                }
                return bestU;
            };
            T xCut1 = T(0.5) * (holeCtrs[hOrder[0]].x + holeCtrs[hOrder[1]].x);
            T xCut2 = T(0.5) * (holeCtrs[hOrder[1]].x + holeCtrs[hOrder[2]].x);
            T xTipL = T(0.5) * (P_L.x + holeCtrs[hOrder[0]].x) - T(0.18) * (holeCtrs[hOrder[0]].x - P_L.x);
            T xTipR = T(0.5) * (holeCtrs[hOrder[2]].x + P_R.x) + T(0.18) * (P_R.x - holeCtrs[hOrder[2]].x);

            T ubL = findCutParam(c0, xTipL), ub1 = findCutParam(c0, xCut1), ub2 = findCutParam(c0, xCut2), ubR = findCutParam(c0, xTipR);
            T utR = findCutParam(c1, xTipR), ut2 = findCutParam(c1, xCut2), ut1 = findCutParam(c1, xCut1), utL = findCutParam(c1, xTipL);
            if (ubL < ub1 && ub1 < ub2 && ub2 < ubR && utR < ut2 && ut2 < ut1 && ut1 < utL) {
                auto bSegs = internal_splitter::splitSegMulti(c0, {ubL, ub1, ub2, ubR});
                auto tSegs = internal_splitter::splitSegMulti(c1, {utR, ut2, ut1, utL});
                // bSegs: [0:PL->B0, 1:B0->B1, 2:B1->B2, 3:B2->B3, 4:B3->PR]
                // tSegs: [0:PR->T3, 1:T3->T2, 2:T2->T1, 3:T1->T0, 4:T0->PL]
                CurveSeg<T> cut1 = mkLine(bSegs[1].p1(), tSegs[2].p1()); // B1 -> T1
                CurveSeg<T> cut2 = mkLine(bSegs[2].p1(), tSegs[1].p1()); // B2 -> T2

                // Middle block [B1->B2, B2->T2, T2->T1, T1->B1] around middle hole hOrder[1]
                bool midOk = addBlockCCW(bSegs[2], cut2, tSegs[2], cut1.reversed(), {&m_domain.loop(hOrder[1] + 1)});

                // Helper for pinched 3-boundary tip-wing around an outer hole (left or right)
                auto buildTipWing = [&](const CurveSeg<T>& segTipB, const CurveSeg<T>& segB,
                                        const CurveSeg<T>& segCut,  const CurveSeg<T>& segT,
                                        const CurveSeg<T>& segTipT, const gsCurveLoop<T>& hl, bool isLeft) -> bool {
                    // Boundary CCW: P_tip ->(segTipB)-> P_b ->(segB)-> P_cutB ->(segCut)-> P_cutT ->(segT)-> P_t ->(segTipT)-> P_tip
                    std::vector<CurveSeg<T>> hCCW = internal_splitter::orientLoopSegments(hl, m_degree, true);
                    // Wait: 3-sector ring around hole touching H_tip, H_bcut, H_tcut leaves H_tip sharing tip quad
                    // Pick vertex H_tip closest to P_tip (must have clear spokes to P_b and P_t)
                    // If K <= 2 (smooth oval/circle), subdivide into 4 equal-arc segments first so 3 sectors cover cleanly.
                    if (hCCW.size() <= 2) {
                        std::vector<CurveSeg<T>> sub;
                        for (const auto& s : hCCW) {
                            auto sp = s.splitArc(T(0.5));
                            sub.push_back(sp.first); sub.push_back(sp.second);
                        }
                        if (sub.size() == 2) {
                            std::vector<CurveSeg<T>> sub4;
                            for (const auto& s : sub) {
                                auto sp = s.splitArc(T(0.5));
                                sub4.push_back(sp.first); sub4.push_back(sp.second);
                            }
                            hCCW = std::move(sub4);
                        } else {
                            hCCW = std::move(sub);
                        }
                    }
                    std::vector<gsPoint2D<T>> hPoly;
                    for (const auto& s : hCCW)
                        for (int k = 0; k < 8; ++k) hPoly.push_back(internal_splitter::evalCurvePoint(*s.bsp, static_cast<T>(k) / T(8)));
                    gsPoint2D<T> P_tip = segTipB.p0(), P_b = segB.p0(), P_t = segTipT.p0();
                    gsPoint2D<T> P_cb = segCut.p0(), P_ct = segCut.p1();
                    const int M = static_cast<int>(hCCW.size());
                    auto distPt = [&](int i, const gsPoint2D<T>& q) {
                        gsPoint2D<T> p = hCCW[i].p0();
                        return (p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y);
                    };
                    // Pick the valid cyclic triple (idxTip, idxB, idxT) whose 4 Coons wing patches pass evaluateRingPatches with lowest distortion
                    std::vector<std::unique_ptr<gsGeometry<T>>> bestWing;
                    T bestScore = T(1e18);
                    bool foundWing = false;
                    for (int iTip = 0; iTip < M; ++iTip) {
                        if (!internal_splitter::isHoleSpokeClearOfHole(hCCW[iTip].p0(), P_b, hPoly) ||
                            !internal_splitter::isHoleSpokeClearOfHole(hCCW[iTip].p0(), P_t, hPoly))
                            continue;
                        for (int dTB = 1; dTB < M - 1; ++dTB) {
                            int iB = (iTip + dTB) % M;
                            if (!internal_splitter::isHoleSpokeClearOfHole(hCCW[iB].p0(), P_cb, hPoly)) continue;
                            for (int dBT = 1; dTB + dBT < M; ++dBT) {
                                int dTT = M - dTB - dBT;
                                int iT = (iB + dBT) % M;
                                if (!internal_splitter::isHoleSpokeClearOfHole(hCCW[iT].p0(), P_ct, hPoly)) continue;
                                std::vector<CurveSeg<T>> s1, s2, s3;
                                for (int j = 0; j < dTB; ++j) s1.push_back(hCCW[(iTip + j) % M]);
                                for (int j = 0; j < dBT; ++j) s2.push_back(hCCW[(iB   + j) % M]);
                                for (int j = 0; j < dTT; ++j) s3.push_back(hCCW[(iT   + j) % M]);
                                auto hasSharpChain = [](const std::vector<CurveSeg<T>>& ch) -> bool {
                                    const T r2d = T(180) / T(3.14159265358979323846);
                                    for (size_t j = 1; j < ch.size(); ++j) {
                                        auto d0 = internal_splitter::evalCurveDeriv(*ch[j - 1].bsp, T(1));
                                        auto d1 = internal_splitter::evalCurveDeriv(*ch[j].bsp, T(0));
                                        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
                                        if (n0 > T(1e-14) && n1 > T(1e-14)) {
                                            T cr = (d0.x * d1.y - d0.y * d1.x) / (n0 * n1);
                                            T dt = (d0.x * d1.x + d0.y * d1.y) / (n0 * n1);
                                            if (std::abs(std::atan2(cr, dt)) * r2d > T(28))
                                                return true;
                                        }
                                    }
                                    return false;
                                };
                                if (hasSharpChain(s1) || hasSharpChain(s2) || hasSharpChain(s3))
                                    continue;
                                auto j1 = internal_splitter::joinCurveSegsExact(s1, m_degree);
                                auto j2 = internal_splitter::joinCurveSegsExact(s2, m_degree);
                                auto j3 = internal_splitter::joinCurveSegsExact(s3, m_degree);
                                CurveSeg<T> h1(*j1), h2(*j2), h3(*j3);
                                bool dmy = false;
                                CurveSeg<T> spTipB = internal_splitter::makeExteriorConeSafeSpoke(hCCW[(iTip + M - 1) % M], hCCW[iTip], hCCW[iTip].p0(), P_b,  hPoly, m_degree, dmy).reversed();
                                CurveSeg<T> spTipT = internal_splitter::makeExteriorConeSafeSpoke(hCCW[(iTip + M - 1) % M], hCCW[iTip], hCCW[iTip].p0(), P_t,  hPoly, m_degree, dmy);
                                CurveSeg<T> spBCut = internal_splitter::makeExteriorConeSafeSpoke(hCCW[(iB   + M - 1) % M], hCCW[iB],   hCCW[iB].p0(),   P_cb, hPoly, m_degree, dmy).reversed();
                                CurveSeg<T> spTCut = internal_splitter::makeExteriorConeSafeSpoke(hCCW[(iT   + M - 1) % M], hCCW[iT],   hCCW[iT].p0(),   P_ct, hPoly, m_degree, dmy);
                                std::vector<std::unique_ptr<gsGeometry<T>>> cand;
                                cand.emplace_back(internal_splitter::makeCoonsPatchFrom4Curves(*segTipB.bsp, *spTipB.bsp, *spTipT.reversed().bsp, *segTipT.reversed().bsp, m_degree));
                                cand.emplace_back(internal_splitter::makeCoonsPatchFrom4Curves(*segB.bsp, *spBCut.bsp, *h1.bsp, *spTipB.bsp, m_degree));
                                cand.emplace_back(internal_splitter::makeCoonsPatchFrom4Curves(*segCut.bsp, *spTCut.reversed().bsp, *h2.bsp, *spBCut.bsp, m_degree));
                                cand.emplace_back(internal_splitter::makeCoonsPatchFrom4Curves(*segT.bsp, *spTipT.reversed().bsp, *h3.bsp, *spTCut.reversed().bsp, m_degree));
                                T sc = 0, md = 0, ma = 0;
                                if (internal_splitter::evaluateRingPatches(cand, hPoly, sc, md, ma)) {
                                    T alignPen = T(0.1) * (distPt(iTip, P_tip) + distPt(iB, P_cb) + distPt(iT, P_ct))
                                               + (ma > T(163.5) ? T(1e5) : T(0));
                                    if (!foundWing || sc + alignPen < bestScore) {
                                        foundWing = true; bestScore = sc + alignPen; bestWing = std::move(cand);
                                    }
                                }
                            }
                        }
                    }
                    if (!foundWing) return false;
                    for (auto& p : bestWing)
                        candidate.multiPatch.addPatch(typename gsGeometry<T>::uPtr(p.release()));
                    (void)isLeft;
                    return true;
                };

                bool lOk = buildTipWing(bSegs[0], bSegs[1], cut1, tSegs[3], tSegs[4], m_domain.loop(hOrder[0] + 1), true);
                bool rOk = buildTipWing(tSegs[0], tSegs[1], cut2.reversed(), bSegs[3], bSegs[4], m_domain.loop(hOrder[2] + 1), false);
                if (midOk && lOk && rOk && tryFinalizeHole(candidate)) return candidate;
            }
            candidate.multiPatch.clear();
        }
    }

    // ------------------------------------------------------------------------
    // Case 8: Star5 with 1 Hole (star_5corner_1hole)
    // ------------------------------------------------------------------------
    if (nHoles == 1 && nOuterCurves == 5) {
        const auto& hole = m_domain.loop(1);
        if (hole.size() == 10) {
            std::vector<gsBSpline<T>> oSegs(10);
            for (int i = 0; i < 5; ++i)
                internal_splitter::splitBSplineExact(outer.curve(i), 0.5, oSegs[2 * i], oSegs[2 * i + 1]);
            for (int i = 0; i < 10; ++i) {
                int h_idx = (9 - i) % 10;
                std::unique_ptr<gsCurve<T>> hCrv(hole.curve(h_idx).clone());
                hCrv->reverse();
                auto pO_s = internal_splitter::startPoint(oSegs[i]);
                auto pO_e = internal_splitter::endPoint(oSegs[i]);
                auto pH_s = internal_splitter::startPoint(*hCrv);
                auto pH_e = internal_splitter::endPoint(*hCrv);
                std::unique_ptr<gsCurve<T>> wsp(internal_splitter::makeStraightBSpline(pO_s, pH_s, m_degree));
                std::unique_ptr<gsCurve<T>> esp(internal_splitter::makeStraightBSpline(pO_e, pH_e, m_degree));
                candidate.multiPatch.addPatch(typename gsGeometry<T>::uPtr(
                    internal_splitter::makeCoonsPatchFrom4Curves(oSegs[i], *esp, *hCrv, *wsp, m_degree)));
            }
            if (tryFinalizeHole(candidate)) return candidate;
            candidate.multiPatch.clear();
        }
    }

    // ------------------------------------------------------------------------
    // Case 9: 12-Corner Domains (Swiss Cross & 6-Point Starburst)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 12) {
        std::vector<int> reflex = findReflexIdx(oCCW);
        // 9a. Swiss Cross (4 reflex inner corners spaced by 3 edges around center)
        if (reflex.size() == 4) {
            int r0 = reflex[0];
            CurveSeg<T> c[12];
            for (int k = 0; k < 12; ++k) c[k] = oCCW[(r0 + k) % 12];
            CurveSeg<T> L[4];
            for (int a = 0; a < 4; ++a)
                L[a] = mkLine(c[3 * a].p0(), c[(3 * (a + 1)) % 12].p0());

            std::vector<std::vector<gsPoint2D<T>>> polys;
            polys.push_back(blockBoundaryPoly(L[0], L[1], L[2], L[3]));
            for (int a = 0; a < 4; ++a)
                polys.push_back(blockBoundaryPoly(c[3*a], c[3*a+1], c[3*a+2], L[a].reversed()));
            auto assigned = assignHolesToPolys(polys);

            bool ok = addBlockCCW(L[0], L[1], L[2], L[3], assigned[0]);
            for (int a = 0; a < 4; ++a)
                ok = ok && addBlockCCW(c[3*a], c[3*a+1], c[3*a+2], L[a].reversed(), assigned[a + 1]);
            if (ok && tryFinalizeHole(candidate)) return candidate;
            candidate.multiPatch.clear();
        }

        // 9b. 6-Point Starburst (6 reflex inner valleys spaced by 2 edges -> 6 kites meeting at O)
        if (reflex.size() == 6) {
            int r0 = reflex[0];
            CurveSeg<T> c[12];
            for (int k = 0; k < 12; ++k) c[k] = oCCW[(r0 + k) % 12];

            // Try origin O at (0,0) and small perturbation candidates; bend spoke cubic if needed to avoid holes
            const gsPoint2D<T> oCands[5] = { {0, 0}, {T(0.03), T(0.02)}, {T(0.02), T(0.04)}, {T(-0.02), T(0.03)}, {T(0.04), T(-0.01)} };
            for (const auto& O : oCands) {
                bool oInHole = false;
                for (int h = 0; h < nHoles; ++h)
                    if (internal_splitter::pointInPolyStrict(O, holePolys[h], T(1e-3))) { oInHole = true; break; }
                if (oInHole) continue;

                CurveSeg<T> rays[6];
                bool raysOk = true;
                for (int k = 0; k < 6; ++k) {
                    gsPoint2D<T> V = c[2 * k].p0();
                    // Check straight ray against all holes; if obstructed, curve quadratic midpoint perpendicular
                    auto rayClear = [&](const gsPoint2D<T>& midCtrl) -> CurveSeg<T> {
                        gsKnotVector<T> kv(T(0), T(1), 0, 3, 1, 2);
                        gsMatrix<T> cf(3, 2);
                        cf << O.x, O.y, midCtrl.x, midCtrl.y, V.x, V.y;
                        gsBSpline<T> bsp(kv, cf);
                        bsp.degreeElevate(std::max(0, m_degree - 2));
                        return CurveSeg<T>(bsp);
                    };
                    auto segClearOfAllHoles = [&](const CurveSeg<T>& seg) -> bool {
                        for (int s = 1; s < 16; ++s) {
                            auto q = internal_splitter::evalCurvePoint(*seg.bsp, static_cast<T>(s) / T(16));
                            for (int h = 0; h < nHoles; ++h)
                                if (internal_splitter::pointInPolyStrict(q, holePolys[h], T(3e-2))) return false;
                        }
                        return true;
                    };
                    CurveSeg<T> testSeg = mkLine(O, V);
                    if (!segClearOfAllHoles(testSeg)) {
                        gsPoint2D<T> mid(T(0.5)*(O.x+V.x), T(0.5)*(O.y+V.y));
                        gsPoint2D<T> dir(V.x - O.x, V.y - O.y);
                        T L = std::hypot(dir.x, dir.y);
                        gsPoint2D<T> perp(-dir.y / std::max(L, T(1e-9)), dir.x / std::max(L, T(1e-9)));
                        bool fixed = false;
                        for (T off : {T(0.10), T(-0.10), T(0.16), T(-0.16), T(0.22), T(-0.22)}) {
                            gsPoint2D<T> mc(mid.x + off * L * perp.x, mid.y + off * L * perp.y);
                            CurveSeg<T> bent = rayClear(mc);
                            if (segClearOfAllHoles(bent)) { testSeg = bent; fixed = true; break; }
                        }
                        if (!fixed) { raysOk = false; break; }
                    }
                    rays[k] = testSeg;
                }
                if (!raysOk) continue;

                std::vector<std::vector<gsPoint2D<T>>> polys;
                for (int k = 0; k < 6; ++k)
                    polys.push_back(blockBoundaryPoly(rays[k], c[2*k], c[2*k+1], rays[(k+1)%6].reversed()));
                auto assigned = assignHolesToPolys(polys);

                bool ok = true;
                candidate.multiPatch.clear();
                for (int k = 0; k < 6; ++k)
                    ok = ok && addBlockCCW(rays[k], c[2*k], c[2*k+1], rays[(k+1)%6].reversed(), assigned[k]);
                if (ok && tryFinalizeHole(candidate)) return candidate;
                candidate.multiPatch.clear();
            }
        }
    }

    // ------------------------------------------------------------------------
    // Case 10: U-Channel 8-Corner (uchannel_8corner_2holes, uchannel_8corner_4holes)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 8 && nHoles >= 1) {
        std::vector<int> reflex = findReflexIdx(oCCW);
        // Find adjacent pair of reflex corners spanning the inner U floor
        int rFloor = -1;
        for (int r : reflex) {
            if (std::find(reflex.begin(), reflex.end(), (r + 1) % 8) != reflex.end()) { rFloor = r; break; }
        }
        if (rFloor >= 0) {
            CurveSeg<T> c[8];
            for (int k = 0; k < 8; ++k) c[k] = oCCW[(rFloor + k) % 8];
            // c[0]: J_R -> J_L (inner floor), c[1]: J_L -> L_topIn, c[2]: L_topIn -> L_topOut,
            // c[3]: L_topOut -> L_botOut, c[4]: L_botOut -> R_botOut (bottom),
            // c[5]: R_botOut -> R_topOut, c[6]: R_topOut -> R_topIn, c[7]: R_topIn -> J_R
            gsPoint2D<T> J_R = c[0].p0(), J_L = c[0].p1();
            auto findParamByX = [&](const CurveSeg<T>& seg, T targetX) -> T {
                T bestU = T(0.5), bestD = T(1e9);
                for (int k = 1; k < 256; ++k) {
                    T u = static_cast<T>(k) / T(256);
                    T d = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, u).x - targetX);
                    if (d < bestD) { bestD = d; bestU = u; }
                }
                T lo = std::max(T(1e-6), bestU - T(1) / T(256));
                T hi = std::min(T(1) - T(1e-6), bestU + T(1) / T(256));
                for (int it = 0; it < 20; ++it) {
                    T m1 = lo + (hi - lo) / T(3), m2 = hi - (hi - lo) / T(3);
                    T d1 = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, m1).x - targetX);
                    T d2 = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, m2).x - targetX);
                    if (d1 < d2) hi = m2; else lo = m1;
                }
                return T(0.5) * (lo + hi);
            };
            auto findParamByY = [&](const CurveSeg<T>& seg, T targetY) -> T {
                T bestU = T(0.5), bestD = T(1e9);
                for (int k = 1; k < 256; ++k) {
                    T u = static_cast<T>(k) / T(256);
                    T d = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, u).y - targetY);
                    if (d < bestD) { bestD = d; bestU = u; }
                }
                T lo = std::max(T(1e-6), bestU - T(1) / T(256));
                T hi = std::min(T(1) - T(1e-6), bestU + T(1) / T(256));
                for (int it = 0; it < 20; ++it) {
                    T m1 = lo + (hi - lo) / T(3), m2 = hi - (hi - lo) / T(3);
                    T d1 = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, m1).y - targetY);
                    T d2 = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, m2).y - targetY);
                    if (d1 < d2) hi = m2; else lo = m1;
                }
                return T(0.5) * (lo + hi);
            };

            T tBL = findParamByX(c[4], J_L.x), tBR = findParamByX(c[4], J_R.x);
            if (tBL > tBR) std::swap(tBL, tBR);
            auto botSegs = internal_splitter::splitSegMulti(c[4], {tBL, tBR}); // [B_left, B_mid, B_right]
            gsPoint2D<T> B_L = botSegs[0].p1(), B_R = botSegs[1].p1();
            CurveSeg<T> cutL = mkLine(B_L, J_L), cutR = mkLine(B_R, J_R);

            // Check min hole y in Left Arm (x < J_L.x) and Right Arm (x > J_R.x)
            T leftMinY = T(1e9), rightMinY = T(1e9);
            for (int h = 0; h < nHoles; ++h) {
                T hMinY = T(1e9);
                for (const auto& p : holePolys[h]) hMinY = std::min(hMinY, p.y);
                if (holeCtrs[h].x < J_L.x) leftMinY = std::min(leftMinY, hMinY);
                if (holeCtrs[h].x > J_R.x) rightMinY = std::min(rightMinY, hMinY);
            }
            T yCutL = (leftMinY  < J_L.y + T(0.03)) ? std::max(B_L.y + T(0.18)*(J_L.y - B_L.y), leftMinY  - T(0.34)) : J_L.y;
            T yCutR = (rightMinY < J_R.y + T(0.03)) ? std::max(B_R.y + T(0.18)*(J_R.y - B_R.y), rightMinY - T(0.34)) : J_R.y;

            T tWL = findParamByY(c[3], yCutL), tWR = findParamByY(c[5], yCutR);
            auto lWall = c[3].splitAtParam(tWL); // first: L_topOut -> W_L, second: W_L -> L_botOut
            auto rWall = c[5].splitAtParam(tWR); // first: R_botOut -> W_R, second: W_R -> R_topOut
            CurveSeg<T> cutArmL = mkLine(lWall.first.p1(), J_L);
            CurveSeg<T> cutArmR = mkLine(J_R, rWall.first.p1());

            std::vector<std::vector<gsPoint2D<T>>> polys = {
                blockBoundaryPoly(botSegs[1], cutR, c[0], cutL.reversed()),          // Base
                blockBoundaryPoly(botSegs[0], cutL, cutArmL.reversed(), lWall.second), // Left-Bot (hole-free)
                blockBoundaryPoly(cutArmL, c[1], c[2], lWall.first),                   // Left-Upper
                blockBoundaryPoly(botSegs[2], rWall.first, cutArmR.reversed(), cutR.reversed()), // Right-Bot (hole-free)
                blockBoundaryPoly(cutArmR, rWall.second, c[6], c[7])                   // Right-Upper
            };
            auto assigned = assignHolesToPolys(polys);

            // Helper: when an upper arm has >= 2 holes and horizontal base cut (both side walls free),
            // slice vertically into ~1:1 square sub-blocks per hole plus optional empty end caps.
            auto addArmCCW = [&](const CurveSeg<T>& bot, const CurveSeg<T>& rgt,
                                 const CurveSeg<T>& top, const CurveSeg<T>& lft,
                                 const std::vector<const gsCurveLoop<T>*>& armHoles) -> bool
            {
                if (armHoles.size() >= 2 && std::abs(bot.p0().y - bot.p1().y) < T(0.02) * std::abs(bot.p1().x - bot.p0().x)) {
                    T w = std::abs(bot.p1().x - bot.p0().x);
                    T yBotA = T(0.5) * (bot.p0().y + bot.p1().y);
                    T yTopA = T(0.5) * (top.p0().y + top.p1().y);
                    T h = yTopA - yBotA;
                    if (w > T(1e-9) && h > T(1.8) * w) {
                        struct ArmHoleInfo { T yMin, yMax, yMid; const gsCurveLoop<T>* ptr; };
                        std::vector<ArmHoleInfo> hInfo;
                        hInfo.reserve(armHoles.size());
                        for (const auto* hp : armHoles) {
                            T hyMin = T(1e18), hyMax = -T(1e18);
                            for (int ci = 0; ci < hp->size(); ++ci) {
                                gsMatrix<T> pr = hp->curve(ci).parameterRange();
                                for (int s = 0; s < 8; ++s) {
                                    T u = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(s) / T(8));
                                    T py = internal_splitter::evalCurvePoint(hp->curve(ci), u).y;
                                    hyMin = std::min(hyMin, py);
                                    hyMax = std::max(hyMax, py);
                                }
                            }
                            hInfo.push_back({hyMin, hyMax, T(0.5) * (hyMin + hyMax), hp});
                        }
                        std::sort(hInfo.begin(), hInfo.end(), [](const ArmHoleInfo& a, const ArmHoleInfo& b) {
                            return a.yMid < b.yMid;
                        });
                        const size_t M = hInfo.size();
                        bool nonOverlap = true;
                        for (size_t i = 0; i + 1 < M; ++i) {
                            if (hInfo[i].yMax + T(0.05) * w >= hInfo[i + 1].yMin) { nonOverlap = false; break; }
                        }
                        if (nonOverlap) {
                            std::vector<T> midCuts(M - 1);
                            for (size_t i = 0; i + 1 < M; ++i)
                                midCuts[i] = T(0.5) * (hInfo[i].yMax + hInfo[i + 1].yMin);

                            std::vector<T> yCuts;
                            T half0 = std::max(T(0.55) * w, midCuts[0] - hInfo[0].yMid);
                            T yLo = hInfo[0].yMid - half0;
                            if (yLo > yBotA + T(0.25) * w && yLo < hInfo[0].yMin - T(0.10) * w)
                                yCuts.push_back(yLo);
                            for (size_t i = 0; i + 1 < M; ++i)
                                yCuts.push_back(midCuts[i]);
                            T halfLast = std::max(T(0.55) * w, hInfo[M - 1].yMid - midCuts[M - 2]);
                            T yHi = hInfo[M - 1].yMid + halfLast;
                            if (yHi < yTopA - T(0.25) * w && yHi > hInfo[M - 1].yMax + T(0.10) * w)
                                yCuts.push_back(yHi);

                            if (!yCuts.empty()) {
                                CurveSeg<T> lftUp = lft.reversed();
                                std::vector<T> tr(yCuts.size()), tl(yCuts.size());
                                for (size_t i = 0; i < yCuts.size(); ++i) {
                                    tr[i] = findParamByY(rgt, yCuts[i]);
                                    tl[i] = findParamByY(lftUp, yCuts[i]);
                                }
                                bool okY = true;
                                for (size_t i = 0; i < yCuts.size(); ++i) {
                                    if (tr[i] <= T(1e-3) || tr[i] >= T(1) - T(1e-3) ||
                                        tl[i] <= T(1e-3) || tl[i] >= T(1) - T(1e-3)) okY = false;
                                    if (i > 0 && (tr[i] <= tr[i - 1] + T(1e-3) || tl[i] <= tl[i - 1] + T(1e-3)))
                                        okY = false;
                                }
                                if (okY) {
                                    auto rSegs = internal_splitter::splitSegMulti(rgt, tr);
                                    auto lSegs = internal_splitter::splitSegMulti(lftUp, tl);
                                    const int N = static_cast<int>(yCuts.size());
                                    std::vector<std::vector<const gsCurveLoop<T>*>> stripHoles(N + 1);
                                    for (const auto& hi : hInfo) {
                                        int hs = 0;
                                        while (hs < N && hi.yMid > yCuts[hs]) ++hs;
                                        stripHoles[hs].push_back(hi.ptr);
                                    }
                                    bool stripOk = true;
                                    for (int i = 0; i <= N; ++i) {
                                        CurveSeg<T> curBot = (i == 0) ? bot : mkLine(lSegs[i - 1].p1(), rSegs[i - 1].p1());
                                        CurveSeg<T> curTopRev = (i == N) ? top : mkLine(rSegs[i].p1(), lSegs[i].p1());
                                        stripOk = stripOk && addBlockCCW(curBot, rSegs[i], curTopRev, lSegs[i].reversed(), stripHoles[i]);
                                    }
                                    return stripOk;
                                }
                            }
                        }
                    }
                }
                return addBlockCCW(bot, rgt, top, lft, armHoles);
            };

            bool ok = true;
            ok = ok && addBlockCCW(botSegs[1], cutR, c[0], cutL.reversed(), assigned[0]);
            ok = ok && addBlockCCW(botSegs[0], cutL, cutArmL.reversed(), lWall.second, assigned[1]);
            ok = ok && addArmCCW(cutArmL, c[1], c[2], lWall.first, assigned[2]);
            ok = ok && addBlockCCW(botSegs[2], rWall.first, cutArmR.reversed(), cutR.reversed(), assigned[3]);
            ok = ok && addArmCCW(cutArmR, rWall.second, c[6], c[7], assigned[4]);
            if (ok && tryFinalizeHole(candidate)) return candidate;
            candidate.multiPatch.clear();
        }
    }

    // ------------------------------------------------------------------------
    // Case 11: Comb 16-Corner with 4 Holes (comb_16corner_4holes, comb_16corner_4holes_curved)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 16 && nHoles >= 1) {
        // Find the back wall curve opposite the teeth (curve whose start and end corners are convex and adjacent walls extend along Y)
        std::vector<int> reflex = findReflexIdx(oCCW);
        if (reflex.size() == 6) {
            // Back wall sits 4 edges before the first reflex valley corner when sorted from right-to-left teeth
            // Identify the two outer corners whose adjacent edges have the two longest straight-Y spans and do NOT touch reflex corners
            int backIdx = -1;
            T minCenterY = T(1e18);
            for (int i = 0; i < 16; ++i) {
                bool isConv0 = (std::find(reflex.begin(), reflex.end(), i) == reflex.end());
                bool isConv1 = (std::find(reflex.begin(), reflex.end(), (i + 1) % 16) == reflex.end());
                if (isConv0 && isConv1) {
                    gsPoint2D<T> m = internal_splitter::evalCurvePoint(*oCCW[i].bsp, T(0.5));
                    if (m.y < minCenterY && oCCW[i].p1().x > oCCW[i].p0().x) {
                        minCenterY = m.y; backIdx = i;
                    }
                }
            }
            if (backIdx >= 0) {
                CurveSeg<T> c[16];
                for (int k = 0; k < 16; ++k) c[k] = oCCW[(backIdx + k) % 16];
                // Split outer side walls c[1] (right) and c[15] (left) at the tooth-floor y level
                T yFloor = c[4].p0().y;
                auto findParamByY = [&](const CurveSeg<T>& seg, T targetY) -> T {
                    T bestU = T(0.5), bestD = T(1e9);
                    for (int k = 1; k < 256; ++k) {
                        T u = static_cast<T>(k) / T(256);
                        T d = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, u).y - targetY);
                        if (d < bestD) { bestD = d; bestU = u; }
                    }
                    return bestU;
                };
                auto findParamByX = [&](const CurveSeg<T>& seg, T targetX) -> T {
                    T bestU = T(0.5), bestD = T(1e9);
                    for (int k = 1; k < 256; ++k) {
                        T u = static_cast<T>(k) / T(256);
                        T d = std::abs(internal_splitter::evalCurvePoint(*seg.bsp, u).x - targetX);
                        if (d < bestD) { bestD = d; bestU = u; }
                    }
                    return bestU;
                };

                auto rWall = c[1].splitAtParam(findParamByY(c[1], yFloor));   // first: bot->yFloor, second: yFloor->top
                auto lWall = c[15].splitAtParam(findParamByY(c[15], yFloor)); // first: top->yFloor, second: yFloor->bot

                // 6 gap corner points ordered left-to-right: G[0..5]
                gsPoint2D<T> G[6] = {
                    c[12].p1(), c[12].p0(),
                    c[8].p1(),  c[8].p0(),
                    c[4].p1(),  c[4].p0()
                };
                // Choose bottom x stations on c[0] so spine vertical/slanted cuts G[k]->B[k] maximize hole clearance
                std::vector<T> bx(6);
                for (int k = 0; k < 6; ++k) bx[k] = G[k].x;
                auto ptSegDist = [](const gsPoint2D<T>& p, const gsPoint2D<T>& a, const gsPoint2D<T>& b) -> T {
                    T dx = b.x - a.x, dy = b.y - a.y;
                    T L2 = dx * dx + dy * dy;
                    T t = (L2 > T(1e-20)) ? std::max(T(0), std::min(T(1), ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2)) : T(0);
                    return std::hypot(a.x + t * dx - p.x, a.y + t * dy - p.y);
                };
                auto segSegDist = [&](const gsPoint2D<T>& a, const gsPoint2D<T>& b,
                                      const gsPoint2D<T>& p, const gsPoint2D<T>& q) -> T {
                    if (internal_splitter::segProperIntersect(a, b, p, q)) return T(0);
                    return std::min(std::min(ptSegDist(a, p, q), ptSegDist(b, p, q)),
                                    std::min(ptSegDist(p, a, b), ptSegDist(q, a, b)));
                };
                auto segMinDistToAllHoles = [&](const gsPoint2D<T>& A, const gsPoint2D<T>& B) -> T {
                    for (int s = 1; s < 16; ++s) {
                        T t = static_cast<T>(s) / T(16);
                        gsPoint2D<T> q((T(1) - t) * A.x + t * B.x, (T(1) - t) * A.y + t * B.y);
                        for (int h = 0; h < nHoles; ++h)
                            if (internal_splitter::pointInPolyStrict(q, holePolys[h], T(1e-7))) return T(-1);
                    }
                    T minD = T(1e18);
                    for (int h = 0; h < nHoles; ++h) {
                        const auto& hp = holePolys[h];
                        const int nh = static_cast<int>(hp.size());
                        for (int i = 0; i < nh; ++i) {
                            T d = segSegDist(A, B, hp[i], hp[(i + 1) % nh]);
                            if (d <= T(1e-7)) return T(-1);
                            minD = std::min(minD, d);
                        }
                    }
                    return minD;
                };
                const T xMinC = std::min(c[0].p0().x, c[0].p1().x);
                const T xMaxC = std::max(c[0].p0().x, c[0].p1().x);
                const T spineH = std::abs(yFloor - c[0].p0().y);
                const T targetClear = std::min(T(0.08), T(0.15) * spineH);
                for (int k = 0; k < 6; ++k) {
                    gsPoint2D<T> Bdef = internal_splitter::evalCurvePoint(*c[0].bsp, findParamByX(c[0], G[k].x));
                    T dDef = segMinDistToAllHoles(G[k], Bdef);
                    if (dDef < targetClear) {
                        T bestX = G[k].x;
                        T bestScore = (dDef > T(0)) ? std::min(dDef, T(0.14)) : T(-1e9);
                        for (int stepIdx = -14; stepIdx <= 14; ++stepIdx) {
                            if (stepIdx == 0) continue;
                            T dx = static_cast<T>(stepIdx) * T(0.025);
                            T candX = G[k].x + dx;
                            if (candX < xMinC + T(0.06) || candX > xMaxC - T(0.06) - T(0.05) * static_cast<T>(5 - k))
                                continue;
                            if (k > 0 && candX < bx[k - 1] + T(0.05)) continue;
                            if (k < 5 && candX > G[k + 1].x + T(0.25)) continue;
                            gsPoint2D<T> Bc = internal_splitter::evalCurvePoint(*c[0].bsp, findParamByX(c[0], candX));
                            T d = segMinDistToAllHoles(G[k], Bc);
                            if (d <= T(0)) continue;
                            T orderPen = (k < 5 && candX > G[k + 1].x - T(0.04)) ? T(0.05) : T(0);
                            T sc = std::min(d, T(0.14)) - T(0.12) * std::abs(dx) - orderPen;
                            if (sc > bestScore) { bestScore = sc; bestX = candX; }
                        }
                        bx[k] = bestX;
                    }
                }
                std::vector<T> bt(6);
                for (int k = 0; k < 6; ++k) bt[k] = findParamByX(c[0], bx[k]);
                bool mono = true;
                for (int k = 1; k < 6; ++k) if (bt[k] <= bt[k-1] + T(1e-3)) mono = false;

                if (mono) {
                    auto bSegs = internal_splitter::splitSegMulti(c[0], bt); // 7 bottom segments
                    CurveSeg<T> spineCut[6];
                    for (int k = 0; k < 6; ++k) spineCut[k] = mkLine(bSegs[k].p1(), G[k]);

                    // 4 tooth base cuts (at yFloor, spanning G[2k-1]->G[2k] and outer walls)
                    CurveSeg<T> tBase0 = mkLine(lWall.first.p1(), G[0]);
                    CurveSeg<T> tBase1 = mkLine(G[1], G[2]);
                    CurveSeg<T> tBase2 = mkLine(G[3], G[4]);
                    CurveSeg<T> tBase3 = mkLine(G[5], rWall.first.p1());

                    // 7 spine strip top curves (left to right: L->G0, G0->G1, G1->G2, G2->G3, G3->G4, G4->G5, G5->R)
                    CurveSeg<T> sTop[7] = {
                        tBase0,
                        c[12].reversed(),
                        tBase1,
                        c[8].reversed(),
                        tBase2,
                        c[4].reversed(),
                        tBase3
                    };

                    // Helper: isolate a ~1:1 square sub-block around a single hole inside a tall free-walled tooth
                    auto addToothCCW = [&](const CurveSeg<T>& bot, const CurveSeg<T>& rgt,
                                           const CurveSeg<T>& top, const CurveSeg<T>& lft,
                                           const std::vector<const gsCurveLoop<T>*>& toothHoles) -> bool
                    {
                        if (toothHoles.size() == 1) {
                            T w = std::abs(bot.p1().x - bot.p0().x);
                            T yBotT = T(0.5) * (bot.p0().y + bot.p1().y);
                            T yTopT = T(0.5) * (top.p0().y + top.p1().y);
                            T h = yTopT - yBotT;
                            if (w > T(1e-9) && h > T(1.65) * w) {
                                const auto& hl = *toothHoles[0];
                                T hyMin = T(1e18), hyMax = -T(1e18);
                                for (int ci = 0; ci < hl.size(); ++ci) {
                                    gsMatrix<T> pr = hl.curve(ci).parameterRange();
                                    for (int s = 0; s < 8; ++s) {
                                        T u = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (static_cast<T>(s) / T(8));
                                        T py = internal_splitter::evalCurvePoint(hl.curve(ci), u).y;
                                        hyMin = std::min(hyMin, py);
                                        hyMax = std::max(hyMax, py);
                                    }
                                }
                                T hyMid = T(0.5) * (hyMin + hyMax);
                                T padY = std::max(T(0.28) * w, T(0.5) * (w - (hyMax - hyMin)));
                                T yLo = hyMin - padY;
                                T yHi = hyMax + padY;
                                std::vector<T> yCuts;
                                if (yLo > yBotT + T(0.18) * w && yLo < hyMin - T(0.05) * w) yCuts.push_back(yLo);
                                if (yHi < yTopT - T(0.18) * w && yHi > hyMax + T(0.05) * w) yCuts.push_back(yHi);
                                if (!yCuts.empty()) {
                                    CurveSeg<T> lftUp = lft.reversed();
                                    std::vector<T> tr(yCuts.size()), tl(yCuts.size());
                                    for (size_t i = 0; i < yCuts.size(); ++i) {
                                        tr[i] = findParamByY(rgt, yCuts[i]);
                                        tl[i] = findParamByY(lftUp, yCuts[i]);
                                    }
                                    bool okY = true;
                                    for (size_t i = 0; i < yCuts.size(); ++i) {
                                        if (tr[i] <= T(1e-3) || tr[i] >= T(1) - T(1e-3) ||
                                            tl[i] <= T(1e-3) || tl[i] >= T(1) - T(1e-3)) okY = false;
                                        if (i > 0 && (tr[i] <= tr[i - 1] + T(1e-3) || tl[i] <= tl[i - 1] + T(1e-3)))
                                            okY = false;
                                    }
                                    if (okY) {
                                        auto rSegs = internal_splitter::splitSegMulti(rgt, tr);
                                        auto lSegs = internal_splitter::splitSegMulti(lftUp, tl);
                                        const int N = static_cast<int>(yCuts.size());
                                        int holeStrip = 0;
                                        while (holeStrip < N && hyMid > yCuts[holeStrip]) ++holeStrip;
                                        bool stripOk = true;
                                        for (int i = 0; i <= N; ++i) {
                                            CurveSeg<T> curBot = (i == 0) ? bot : mkLine(lSegs[i - 1].p1(), rSegs[i - 1].p1());
                                            CurveSeg<T> curTopRev = (i == N) ? top : mkLine(rSegs[i].p1(), lSegs[i].p1());
                                            std::vector<const gsCurveLoop<T>*> sh;
                                            if (i == holeStrip) sh = toothHoles;
                                            stripOk = stripOk && addBlockCCW(curBot, rSegs[i], curTopRev, lSegs[i].reversed(), sh);
                                        }
                                        return stripOk;
                                    }
                                }
                            }
                        }
                        return addBlockCCW(bot, rgt, top, lft, toothHoles);
                    };

                    // 7 spine strips + 4 teeth = 11 macro-blocks
                    std::vector<std::vector<gsPoint2D<T>>> polys;
                    for (int i = 0; i < 7; ++i) {
                        CurveSeg<T> wL = (i == 0) ? lWall.second.reversed() : spineCut[i - 1];
                        CurveSeg<T> wR = (i == 6) ? rWall.first : spineCut[i];
                        polys.push_back(blockBoundaryPoly(bSegs[i], wR, sTop[i].reversed(), wL.reversed()));
                    }
                    polys.push_back(blockBoundaryPoly(tBase0, c[13], c[14], lWall.first));
                    polys.push_back(blockBoundaryPoly(tBase1, c[9],  c[10], c[11]));
                    polys.push_back(blockBoundaryPoly(tBase2, c[5],  c[6],  c[7]));
                    polys.push_back(blockBoundaryPoly(tBase3, rWall.second, c[2], c[3]));

                    auto assigned = assignHolesToPolys(polys);
                    bool ok = true;
                    for (int i = 0; i < 7; ++i) {
                        CurveSeg<T> wL = (i == 0) ? lWall.second.reversed() : spineCut[i - 1];
                        CurveSeg<T> wR = (i == 6) ? rWall.first : spineCut[i];
                        ok = ok && addBlockCCW(bSegs[i], wR, sTop[i].reversed(), wL.reversed(), assigned[i]);
                    }
                    ok = ok && addToothCCW(tBase0, c[13], c[14], lWall.first, assigned[7]);
                    ok = ok && addToothCCW(tBase1, c[9],  c[10], c[11],       assigned[8]);
                    ok = ok && addToothCCW(tBase2, c[5],  c[6],  c[7],        assigned[9]);
                    ok = ok && addToothCCW(tBase3, rWall.second, c[2], c[3],  assigned[10]);

                    if (ok && tryFinalizeHole(candidate)) return candidate;
                    candidate.multiPatch.clear();

                    // Second chance for spine-dominated combs: add a hole-free transition band at yMid
                    // so vertical cuts mx[k] below yMid frame every spine hole in an aligned rectangular cell.
                    T maxHoleY = -T(1e18);
                    for (int h = 0; h < nHoles; ++h)
                        for (const auto& hp : holePolys[h]) maxHoleY = std::max(maxHoleY, hp.y);
                    if (c[0].p0().y < yFloor && maxHoleY < yFloor - T(0.18) * spineH) {
                        T yMid = std::min(yFloor - T(0.16) * spineH, maxHoleY + T(0.12));
                        if (yMid > c[0].p0().y + T(0.20) * spineH && yMid > maxHoleY + T(0.05)) {
                            std::vector<T> mx(6);
                            for (int k = 0; k < 6; ++k) mx[k] = G[k].x;
                            const T vClear = std::min(T(0.12), T(0.20) * spineH);
                            for (int k = 0; k < 6; ++k) {
                                gsPoint2D<T> Bdef = internal_splitter::evalCurvePoint(*c[0].bsp, findParamByX(c[0], G[k].x));
                                gsPoint2D<T> Mdef(Bdef.x, yMid);
                                T dDef = std::min(segMinDistToAllHoles(Bdef, Mdef), segMinDistToAllHoles(Mdef, G[k]));
                                if (dDef < vClear) {
                                    T bestX = G[k].x;
                                    T bestSc = (dDef > T(0)) ? std::min(dDef, vClear) : T(-1e9);
                                    for (int stepIdx = -14; stepIdx <= 14; ++stepIdx) {
                                        if (stepIdx == 0) continue;
                                        T dx = static_cast<T>(stepIdx) * T(0.025);
                                        T candX = G[k].x + dx;
                                        if (candX < xMinC + T(0.08) || candX > xMaxC - T(0.08) - T(0.06) * static_cast<T>(5 - k))
                                            continue;
                                        T minLeft = (k > 0) ? std::min(T(0.08), T(0.30) * (G[k].x - G[k - 1].x)) : T(0.08);
                                        if (k > 0 && candX < mx[k - 1] + minLeft) continue;
                                        T minRight = (k < 5) ? std::min(T(0.14), T(0.45) * (G[k + 1].x - G[k].x)) : T(0.08);
                                        if (k < 5 && candX > G[k + 1].x - minRight) continue;
                                        gsPoint2D<T> Bc = internal_splitter::evalCurvePoint(*c[0].bsp, findParamByX(c[0], candX));
                                        gsPoint2D<T> Mc(Bc.x, yMid);
                                        T d = std::min(segMinDistToAllHoles(Bc, Mc), segMinDistToAllHoles(Mc, G[k]));
                                        if (d <= T(0)) continue;
                                        T sc = std::min(d, vClear) - T(0.12) * std::abs(dx);
                                        if (sc > bestSc) { bestSc = sc; bestX = candX; }
                                    }
                                    mx[k] = bestX;
                                }
                            }
                            std::vector<T> mt(6);
                            for (int k = 0; k < 6; ++k) mt[k] = findParamByX(c[0], mx[k]);
                            bool mono2 = true;
                            for (int k = 1; k < 6; ++k) if (mt[k] <= mt[k - 1] + T(1e-3)) mono2 = false;
                            if (mono2) {
                                auto bSegs2 = internal_splitter::splitSegMulti(c[0], mt);
                                auto lMidSplit = lWall.second.splitAtParam(findParamByY(lWall.second, yMid)); // yFloor->yMid, yMid->yBot
                                auto rMidSplit = rWall.first.splitAtParam(findParamByY(rWall.first, yMid));   // yBot->yMid, yMid->yFloor
                                gsPoint2D<T> M[6];
                                for (int k = 0; k < 6; ++k) M[k] = gsPoint2D<T>(bSegs2[k].p1().x, yMid);
                                CurveSeg<T> mCut[7], loCut[6], upCut[6];
                                for (int k = 0; k < 6; ++k) {
                                    loCut[k] = mkLine(bSegs2[k].p1(), M[k]);
                                    upCut[k] = mkLine(M[k], G[k]);
                                }
                                for (int i = 0; i < 7; ++i) {
                                    gsPoint2D<T> P0 = (i == 0) ? lMidSplit.first.p1() : M[i - 1];
                                    gsPoint2D<T> P1 = (i == 6) ? rMidSplit.first.p1() : M[i];
                                    mCut[i] = mkLine(P0, P1);
                                }
                                std::vector<std::vector<gsPoint2D<T>>> polys2;
                                for (int i = 0; i < 7; ++i) {
                                    CurveSeg<T> loL = (i == 0) ? lMidSplit.second.reversed() : loCut[i - 1];
                                    CurveSeg<T> loR = (i == 6) ? rMidSplit.first : loCut[i];
                                    polys2.push_back(blockBoundaryPoly(bSegs2[i], loR, mCut[i].reversed(), loL.reversed()));
                                }
                                auto assigned2 = assignHolesToPolys(polys2);
                                bool ok2 = true;
                                for (int i = 0; i < 7; ++i) {
                                    CurveSeg<T> loL = (i == 0) ? lMidSplit.second.reversed() : loCut[i - 1];
                                    CurveSeg<T> loR = (i == 6) ? rMidSplit.first : loCut[i];
                                    ok2 = ok2 && addBlockCCW(bSegs2[i], loR, mCut[i].reversed(), loL.reversed(), assigned2[i]);

                                    CurveSeg<T> upL = (i == 0) ? lMidSplit.first.reversed() : upCut[i - 1];
                                    CurveSeg<T> upR = (i == 6) ? rMidSplit.second : upCut[i];
                                    ok2 = ok2 && addBlockCCW(mCut[i], upR, sTop[i].reversed(), upL.reversed(), {});
                                }
                                std::vector<const gsCurveLoop<T>*> noHoles;
                                ok2 = ok2 && addBlockCCW(tBase0, c[13], c[14], lWall.first, noHoles);
                                ok2 = ok2 && addBlockCCW(tBase1, c[9],  c[10], c[11],       noHoles);
                                ok2 = ok2 && addBlockCCW(tBase2, c[5],  c[6],  c[7],        noHoles);
                                ok2 = ok2 && addBlockCCW(tBase3, rWall.second, c[2], c[3],  noHoles);

                                if (ok2 && tryFinalizeHole(candidate)) return candidate;
                                candidate.multiPatch.clear();
                            }
                        }
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // Case 12: General Conforming Multi-Hole Fallback (only if requested directly)
    // ------------------------------------------------------------------------
    if (m_strategy == Strategy::MedialAxisSkeleton)
        return computeMultiHoleTriangulationSplitting();
    return candidate;
}

// ============================================================================
// Approach 3: Topological Preset Template Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeTemplateSplitting()
{
    SplittingCandidate candidate;
    candidate.methodUsed = "TemplateDecomposition";

    if (m_domain.numLoops() == 0) {
        candidate.isValid = false;
        return candidate;
    }

    const auto& outer = m_domain.outer();
    int nHoles = m_domain.numHoles();
    int nOuterCurves = outer.size();
    const T domMinAng = getDomainMinCornerAngleDeg();
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    if (nHoles == 0) {
        candidate.numPatches = 0;
        candidate.isValid = false;
        return candidate;
    }

    // 1-hole with 4 curves
    if (nHoles == 1 && nOuterCurves == 4 && m_domain.loop(1).size() == 4) {
        internal_splitter::buildQuadRingByProximity(outer, m_domain.loop(1), candidate.multiPatch, m_degree);
        candidate.multiPatch.computeTopology(T(1e-5), false, false);
        optimizeWinslowBarrier(candidate.multiPatch);
        if (validateMultiPatch(m_domain, candidate.multiPatch, T(1e-4), nullptr, minAngTol, maxAngTol)) {
            candidate.numPatches = 4;
            candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet,
                                                                &candidate.minCornerAngleDeg, &candidate.maxCornerAngleDeg);
            candidate.isValid = (candidate.minJacobianDet > 0);
            return candidate;
        }
        candidate.multiPatch.clear();
    }

    if (m_strategy == Strategy::TemplateDecomposition) {
        // General conforming multi-hole fallback (guaranteed 0 hanging nodes & exact boundary match)
        candidate = computeMultiHoleTriangulationSplitting();
        if (candidate.isValid)
            candidate.methodUsed = "TemplateDecomposition";
    }
    return candidate;
}

// ============================================================================
// Mode A Helpers (Existing gsMultiPatch<T> -> Coarse Conforming Layout)
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::extractPlanarDomainFromMultiPatch(const gsMultiPatch<T>& mp, gsPlanarDomain<T>& outDom)
{
    using namespace internal_splitter;
    if (mp.nPatches() == 0) return false;

    gsMultiPatch<T> mpCopy(mp);
    mpCopy.computeTopology(T(1e-4), false, false);

    struct BndEdge {
        std::unique_ptr<gsCurve<T>> crv;
        gsPoint2D<T> p0, p1;
        bool used = false;
    };
    std::vector<BndEdge> edges;
    for (auto it = mpCopy.bBegin(); it != mpCopy.bEnd(); ++it) {
        auto g = mpCopy.patch(it->patch).boundary(it->side());
        std::unique_ptr<gsCurve<T>> crv(dynamic_cast<gsCurve<T>*>(g.release()));
        if (!crv) return false;
        gsPoint2D<T> a = startPoint(*crv);
        gsPoint2D<T> b = endPoint(*crv);
        edges.push_back({std::move(crv), a, b, false});
    }
    if (edges.empty()) return false;

    gsMatrix<T> bbox;
    mpCopy.boundingBox(bbox);
    const T diag = std::hypot(bbox(0, 1) - bbox(0, 0), bbox(1, 1) - bbox(1, 0));
    const T tol = std::max(T(1e-6), T(1e-4) * diag);

    std::vector<std::unique_ptr<gsCurveLoop<T>>> loops;
    std::vector<T> loopAreas;

    for (size_t s = 0; s < edges.size(); ++s) {
        if (edges[s].used) continue;
        edges[s].used = true;
        std::vector<gsCurve<T>*> crvs;
        crvs.push_back(edges[s].crv->clone().release());
        gsPoint2D<T> head = edges[s].p0;
        gsPoint2D<T> tail = edges[s].p1;

        while (tail.dist(head) > tol) {
            int best = -1;
            bool rev = false;
            T bestD = tol * T(10);
            for (size_t k = 0; k < edges.size(); ++k) {
                if (edges[k].used) continue;
                T d0 = tail.dist(edges[k].p0);
                T d1 = tail.dist(edges[k].p1);
                if (d0 < bestD) { bestD = d0; best = static_cast<int>(k); rev = false; }
                if (d1 < bestD) { bestD = d1; best = static_cast<int>(k); rev = true; }
            }
            if (best < 0) break;
            edges[best].used = true;
            std::unique_ptr<gsCurve<T>> c(edges[best].crv->clone().release());
            if (rev) {
                c->reverse();
                tail = edges[best].p0;
            } else {
                tail = edges[best].p1;
            }
            crvs.push_back(c.release());
        }
        auto lp = std::unique_ptr<gsCurveLoop<T>>(new gsCurveLoop<T>(crvs));
        auto poly = sampleCurveLoopPoly(*lp, 20);
        T area2 = T(0);
        for (size_t i = 0; i < poly.size(); ++i) {
            const auto& p = poly[i];
            const auto& q = poly[(i + 1) % poly.size()];
            area2 += p.x * q.y - q.x * p.y;
        }
        loopAreas.push_back(std::abs(area2) * T(0.5));
        loops.push_back(std::move(lp));
    }

    if (loops.empty()) return false;
    size_t outerIdx = 0;
    for (size_t i = 1; i < loops.size(); ++i)
        if (loopAreas[i] > loopAreas[outerIdx]) outerIdx = i;

    outDom = gsPlanarDomain<T>(loops[outerIdx].release());
    for (size_t i = 0; i < loops.size(); ++i) {
        if (i == outerIdx) continue;
        outDom.insertHole(loops[i].release());
    }
    return true;
}

template <class T>
bool gsPatchSplitter<T>::splitDegenerateSinglePatchOGrid(const gsMultiPatch<T>& inMp, gsMultiPatch<T>& outMp) const
{
    using namespace internal_splitter;
    if (inMp.nPatches() != 1) return false;

    const auto& p0 = inMp.patch(0);
    std::unique_ptr<gsCurve<T>> cS(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::south).release()));
    std::unique_ptr<gsCurve<T>> cE(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::east).release()));
    std::unique_ptr<gsCurve<T>> cN(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::north).release()));
    std::unique_ptr<gsCurve<T>> cW(dynamic_cast<gsCurve<T>*>(p0.boundary(boundary::west).release()));
    if (!cS || !cE || !cN || !cW) return false;
    cN->reverse();
    cW->reverse();

    int deg = std::max(m_degree, m_options.getInt("Degree"));
    std::vector<CurveSeg<T>> bEdges = {
        CurveSeg<T>(*cS, deg),
        CurveSeg<T>(*cE, deg),
        CurveSeg<T>(*cN, deg),
        CurveSeg<T>(*cW, deg)
    };

    // Ensure CCW ordering of bEdges
    gsPoint2D<T> P[4];
    for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
    std::vector<gsPoint2D<T>> poly;
    for (int k = 0; k < 4; ++k) {
        for (int s = 0; s < 16; ++s) {
            T u = static_cast<T>(s) / T(16);
            poly.push_back(evalCurvePoint(*bEdges[k].bsp, u));
        }
    }
    T area2 = T(0);
    for (size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        area2 += a.x * b.y - b.x * a.y;
    }
    if (area2 < T(0)) {
        std::reverse(bEdges.begin(), bEdges.end());
        for (auto& e : bEdges) e = e.reversed();
        for (int k = 0; k < 4; ++k) P[k] = bEdges[k].p0();
        std::reverse(poly.begin(), poly.end());
    }

    gsPoint2D<T> ctr = polygonCentroid(poly);
    if (!pointInPolyStrict(ctr, poly, T(1e-5))) {
        ctr = gsPoint2D<T>((P[0].x + P[1].x + P[2].x + P[3].x) * T(0.25),
                           (P[0].y + P[1].y + P[2].y + P[3].y) * T(0.25));
    }

    T bestEnergy = std::numeric_limits<T>::max();
    gsMultiPatch<T> bestMp;
    const T alphas[3] = {T(0.45), T(0.35), T(0.55)};
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    for (T alpha : alphas) {
        gsPoint2D<T> Q[4];
        bool ok = true;
        for (int k = 0; k < 4; ++k) {
            Q[k] = gsPoint2D<T>(ctr.x + alpha * (P[k].x - ctr.x),
                                ctr.y + alpha * (P[k].y - ctr.y));
            if (!pointInPolyStrict(Q[k], poly, T(1e-5))) { ok = false; break; }
        }
        if (!ok || !isStrictlyConvexQuad(Q)) continue;

        CurveSeg<T> qEdge[4], spoke[4];
        for (int k = 0; k < 4; ++k) {
            qEdge[k] = makeStraightSeg(Q[k], Q[(k + 1) % 4], deg);
            spoke[k] = makeStraightSeg(Q[k], P[k], deg);
        }

        gsMultiPatch<T> candMp;
        candMp.addPatch(typename gsGeometry<T>::uPtr(
            makeCoonsFromCCWCycle(qEdge[0], qEdge[1], qEdge[2], qEdge[3], deg)));
        for (int k = 0; k < 4; ++k) {
            int kn = (k + 1) % 4;
            candMp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(qEdge[k], spoke[kn], bEdges[k].reversed(), spoke[k].reversed(), deg)));
        }

        gsASG1LocalFitter<T>::unifyPatchOrientations(candMp);
        candMp.computeTopology(T(1e-5), false, false);
        optimizeWinslowBarrier(candMp);
        candMp.computeTopology(T(1e-5), false, false);
        T minDet = 0, minA = 0, maxA = 180;
        T energy = evaluateDistortionEnergy(candMp, &minDet, &minA, &maxA);
        if (minDet > T(1e-7) && maxA <= maxAngTol && energy < bestEnergy) {
            bestEnergy = energy;
            bestMp = give(candMp);
        }
    }

    if (bestMp.nPatches() == 5) {
        outMp = give(bestMp);
        return true;
    }
    return false;
}

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeModeACoarseLayout()
{
    using namespace internal_splitter;
    SplittingCandidate cand;
    cand.methodUsed = "ModeAPolychord";
    if (m_inputMp.nPatches() == 0) return cand;

    gsMultiPatch<T> mp(m_inputMp);
    gsASG1LocalFitter<T>::unifyPatchOrientations(mp);
    mp.computeTopology(T(1e-4), false, false);

    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    const T minAngTol = static_cast<T>(m_options.getReal("MinCornerAngleDeg"));
    bool hasDegenerateCorners = false;
    bool hasReflexCorners = false;
    bool hasBadJacobian = false;

    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pMinDet = 0, pMaxDet = 0, pMinAbs = 0, pEn = 0;
        evalPatchJacobianAndWinslow(mp.patch(p), 17, pMinDet, pMaxDet, pMinAbs, pEn);
        if (pMinDet <= T(1e-7) && pMaxDet >= -T(1e-7))
            hasBadJacobian = true;

        T angs[4], minA = 0, maxA = 180;
        bool degen = false;
        evalPatchCornerAnglesDeg(mp.patch(p), angs, minA, maxA, degen);
        if (degen || maxA > maxAngTol || minA < minAngTol)
            hasDegenerateCorners = true;
        if (maxA >= T(200))
            hasReflexCorners = true;
    }

    // 1. Single patch with ~180 deg smooth corners (e.g. 1-patch disk/ellipse) -> 5-patch O-grid
    if (mp.nPatches() == 1 && (hasDegenerateCorners || hasBadJacobian) && !hasReflexCorners) {
        gsMultiPatch<T> ogridMp;
        if (splitDegenerateSinglePatchOGrid(mp, ogridMp)) {
            cand.multiPatch = give(ogridMp);
            cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
            cand.distortionEnergy = evaluateDistortionEnergy(cand.multiPatch, &cand.minJacobianDet,
                                                             &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
            cand.methodUsed = "ModeAButterflyOGrid";
            cand.isValid = (cand.minJacobianDet > T(1e-8));
            return cand;
        }
    }

    // 2. If degenerate/reflex/invalid corners remain, extract planar domain and re-split via Mode B
    if (hasDegenerateCorners || hasReflexCorners || hasBadJacobian) {
        gsPlanarDomain<T> extractedDom;
        if (extractPlanarDomainFromMultiPatch(mp, extractedDom)) {
            m_domain = extractedDom;
            m_hasDomain = true;
            if (computeSplitting()) {
                return m_bestCandidate;
            }
        }
    }

    // 3. Otherwise keep existing conforming coarse topology for polychord optimization
    cand.multiPatch = give(mp);
    cand.numPatches = static_cast<int>(cand.multiPatch.nPatches());
    cand.distortionEnergy = evaluateDistortionEnergy(cand.multiPatch, &cand.minJacobianDet,
                                                     &cand.minCornerAngleDeg, &cand.maxCornerAngleDeg);
    cand.isValid = (cand.numPatches > 0 && cand.minJacobianDet > T(1e-8));
    return cand;
}

// ============================================================================
// Automated Winner-Selection Metric, Full 4-Stage Pipeline & Execution
// ============================================================================

template <class T>
gsMultiPatch<T> gsPatchSplitter<T>::compute()
{
    m_degree = std::max(1, m_options.getInt("Degree"));
    int stratOpt = m_options.getInt("Strategy");
    if (stratOpt >= 0 && stratOpt <= 4)
        m_strategy = static_cast<Strategy>(stratOpt);

    // Stage 1: Coarse Conforming Quad Layout (Mode A or Mode B)
    if (m_hasInputMp && m_inputMp.nPatches() > 0 && (!m_hasDomain || m_strategy == Strategy::PolychordOnly)) {
        m_bestCandidate = computeModeACoarseLayout();
    } else if (m_hasDomain && m_domain.numLoops() > 0) {
        computeSplitting();
    } else if (m_hasInputMp && m_inputMp.nPatches() > 0) {
        m_bestCandidate = computeModeACoarseLayout();
    } else {
        m_bestCandidate = SplittingCandidate();
        return gsMultiPatch<T>();
    }

    if (m_bestCandidate.multiPatch.nPatches() == 0)
        return m_bestCandidate.multiPatch;

    gsMultiPatch<T> workMp = m_bestCandidate.multiPatch;

    // Stage 2a: Dual-Sheet / Polychord Strip Reduction (Mode B, aiming for minimal patch count)
    if (m_hasDomain && !m_hasInputMp && workMp.nPatches() > 2) {
        reducePolychordStrips(workMp);
    }

    // Stage 2b: Global Polychord Aspect-Ratio Optimization (zero hanging nodes)
    // Skip re-splitting Mode B multipatches when already single convex quad or >= 6 patches.
    const bool keepModeBCoarse =
        m_hasDomain && !m_hasInputMp &&
        ((m_domain.numHoles() == 0 && m_domain.outer().size() == 4 && workMp.nPatches() == 1) ||
         workMp.nPatches() >= 6);
    if (m_options.getSwitch("EnablePolychordSplitting") && !keepModeBCoarse) {
        optimizePolychordSplits(workMp);
    }

    // Stage 3: Knot Harmonization, Relaxation & Bijective Parameterization
    relaxAndParameterize(workMp, true);

    // Stage 4: ArgyrisSplines Post-Processing & Verification
    formatForArgyris(workMp);

    // Update cached candidate metrics
    m_bestCandidate.multiPatch = workMp;
    m_bestCandidate.numPatches = static_cast<int>(workMp.nPatches());
    m_bestCandidate.distortionEnergy = evaluateDistortionEnergy(
        workMp, &m_bestCandidate.minJacobianDet,
        &m_bestCandidate.minCornerAngleDeg, &m_bestCandidate.maxCornerAngleDeg);
    m_bestCandidate.isValid = (m_bestCandidate.numPatches > 0 && m_bestCandidate.minJacobianDet > T(1e-9));

    return m_bestCandidate.multiPatch;
}

template <class T>
bool gsPatchSplitter<T>::computeSplitting()
{
    if (m_hasInputMp && !m_hasDomain) {
        m_bestCandidate = computeModeACoarseLayout();
        if (m_bestCandidate.multiPatch.nPatches() > 0)
            ensureExactBoundaryConformity(m_bestCandidate.multiPatch);
        return m_bestCandidate.isValid;
    }

    if (m_strategy == Strategy::PolychordOnly && m_hasInputMp) {
        m_bestCandidate = computeModeACoarseLayout();
        if (m_bestCandidate.isValid) {
            optimizePolychordSplits(m_bestCandidate.multiPatch);
            relaxAndParameterize(m_bestCandidate.multiPatch, true);
            m_bestCandidate.numPatches = static_cast<int>(m_bestCandidate.multiPatch.nPatches());
        }
        return m_bestCandidate.isValid;
    }

    if (m_strategy == Strategy::HarmonicCrossField) {
        m_bestCandidate = computeCrossFieldSplitting();
        if (!m_bestCandidate.isValid) {
            m_bestCandidate = computeMedialAxisSplitting();
        }
    } else if (m_strategy == Strategy::MedialAxisSkeleton) {
        m_bestCandidate = computeMedialAxisSplitting();
        if (!m_bestCandidate.isValid) {
            m_bestCandidate = computeTemplateSplitting();
        }
    } else if (m_strategy == Strategy::TemplateDecomposition) {
        m_bestCandidate = computeTemplateSplitting();
    } else {
        SplittingCandidate matCand;
        SplittingCandidate cfCand;
        SplittingCandidate tmplCand;
        // Fast-path: convexish 4-corner hole-free quad -> keep as 1 patch cleanly
        if (m_domain.numHoles() == 0 && m_domain.outer().size() == 4) {
            matCand = computeMedialAxisSplitting();
            if (matCand.isValid && matCand.numPatches == 1) {
                m_bestCandidate = matCand;
                return true;
            }
        }
        cfCand = computeCrossFieldSplitting();
        // Run fast MAT/template candidates without repeating multi-hole triangulation
        if (!matCand.isValid)
            matCand = computeMedialAxisSplitting();
        if (m_domain.numHoles() <= 1)
            tmplCand = computeTemplateSplitting();

        std::vector<SplittingCandidate> validPool;
        if (cfCand.isValid && cfCand.numPatches > 0) validPool.push_back(cfCand);
        if (matCand.isValid && matCand.numPatches > 0) validPool.push_back(matCand);
        if (tmplCand.isValid && tmplCand.numPatches > 0) validPool.push_back(tmplCand);

        if (validPool.empty()) {
            m_bestCandidate = matCand.numPatches > 0 ? matCand : cfCand;
            return false;
        }

        const bool preferSingleQuad = (m_domain.numHoles() == 0 && m_domain.outer().size() == 4);
        auto computeQ = [this](const SplittingCandidate& c) -> T {
            return m_weightPatches * c.numPatches + m_weightDistortion * c.distortionEnergy;
        };

        std::sort(validPool.begin(), validPool.end(),
                  [&](const SplittingCandidate& a, const SplittingCandidate& b) {
                      if (preferSingleQuad && (a.numPatches == 1) != (b.numPatches == 1))
                          return a.numPatches == 1;
                      return computeQ(a) < computeQ(b);
                  });

        m_bestCandidate = validPool.front();
    }

    if (m_bestCandidate.multiPatch.nPatches() > 0) {
        ensureExactBoundaryConformity(m_bestCandidate.multiPatch);
    }

    return m_bestCandidate.isValid;
}

// ============================================================================
// Optimization & Quality Evaluation (Winslow Barrier Functional)
// ============================================================================

template <class T>
T gsPatchSplitter<T>::evaluateDistortionEnergy(const gsMultiPatch<T>& mp, T* minDetOut,
                                               T* minAngleOut, T* maxAngleOut) const
{
    if (mp.nPatches() == 0) {
        if (minDetOut) *minDetOut = T(0);
        if (minAngleOut) *minAngleOut = T(0);
        if (maxAngleOut) *maxAngleOut = T(180);
        return T(1e30);
    }

    T totalEnergy = 0.0;
    T globalMinAbsDet = std::numeric_limits<T>::max();
    T globalMinAng = T(360);
    T globalMaxAng = T(0);
    bool allSingleSign = true;

    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pMinDet = 0, pMaxDet = 0, pMinAbs = 0, pEnergy = 0;
        internal_splitter::evalPatchJacobianAndWinslow(mp.patch(p), 17, pMinDet, pMaxDet, pMinAbs, pEnergy);

        bool pSingleSign = (pMinDet > T(1e-10)) || (pMaxDet < -T(1e-10));
        if (!pSingleSign) {
            allSingleSign = false;
            pEnergy += T(1e6);
        }
        if (pMinAbs < globalMinAbsDet) globalMinAbsDet = pMinAbs;

        T angs[4], pMinA = 0, pMaxA = 180;
        bool degen = false;
        internal_splitter::evalPatchCornerAnglesDeg(mp.patch(p), angs, pMinA, pMaxA, degen);
        globalMinAng = std::min(globalMinAng, pMinA);
        globalMaxAng = std::max(globalMaxAng, pMaxA);

        T angPen = T(0);
        for (int k = 0; k < 4; ++k) {
            T a = angs[k];
            if (a < T(35))  angPen += ((T(35) - a) / T(10)) * ((T(35) - a) / T(10));
            if (a > T(140)) angPen += ((a - T(140)) / T(10)) * ((a - T(140)) / T(10));
        }
        if (degen) angPen += T(1e4);

        totalEnergy += pEnergy + T(5.0) * angPen;
    }

    if (minDetOut) *minDetOut = allSingleSign ? globalMinAbsDet : -globalMinAbsDet;
    if (minAngleOut) *minAngleOut = globalMinAng;
    if (maxAngleOut) *maxAngleOut = globalMaxAng;
    return totalEnergy;
}

template <class T>
bool gsPatchSplitter<T>::optimizeWinslowBarrier(gsMultiPatch<T>& mp, T minDetThreshold) const
{
    if (mp.nPatches() == 0) return false;

    bool allOk = true;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        auto* tp = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;

        T minDet = 0, maxDet = 0, minAbs = 0, curEnergy = 0;
        internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, minDet, maxDet, minAbs, curEnergy);

        // If patch is uniformly negative everywhere, flip parameter orientation without altering geometry
        if (maxDet < -minDetThreshold) {
            tp->reverse(0);
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, minDet, maxDet, minAbs, curEnergy);
        }

        int n_u = tp->basis().size(0);
        int n_v = tp->basis().size(1);
        gsMatrix<T>& coefs = tp->coefs();
        bool valid = (minDet > minDetThreshold) || (maxDet < -minDetThreshold);

        // If Coons interior control points slightly overshoot on deep concave edges, try blending toward
        // corner-bilinear interior CPs while keeping boundary CPs strictly fixed.
        if (!valid && n_u > 2 && n_v > 2) {
            gsMatrix<T> origCoefs = coefs;
            gsPoint2D<T> P00(coefs(0, 0), coefs(0, 1));
            gsPoint2D<T> P10(coefs(n_u - 1, 0), coefs(n_u - 1, 1));
            gsPoint2D<T> P01(coefs((n_v - 1) * n_u, 0), coefs((n_v - 1) * n_u, 1));
            gsPoint2D<T> P11(coefs(n_u - 1 + (n_v - 1) * n_u, 0), coefs(n_u - 1 + (n_v - 1) * n_u, 1));

            const auto& kv_u = tp->basis().knots(0);
            const auto& kv_v = tp->basis().knots(1);
            gsMatrix<T> gr_u = kv_u.greville();
            gsMatrix<T> gr_v = kv_v.greville();

            T bestScore = minDet;
            gsMatrix<T> bestCoefs = origCoefs;

            for (T blend = T(0.2); blend <= T(1.01); blend += T(0.2)) {
                for (int j = 1; j < n_v - 1; ++j) {
                    T eta = gr_v(0, j);
                    for (int i = 1; i < n_u - 1; ++i) {
                        T xi = gr_u(0, i);
                        T bx = (1 - xi) * (1 - eta) * P00.x + xi * (1 - eta) * P10.x +
                               (1 - xi) * eta * P01.x + xi * eta * P11.x;
                        T by = (1 - xi) * (1 - eta) * P00.y + xi * (1 - eta) * P10.y +
                               (1 - xi) * eta * P01.y + xi * eta * P11.y;
                        int idx = i + j * n_u;
                        coefs(idx, 0) = (1 - blend) * origCoefs(idx, 0) + blend * bx;
                        coefs(idx, 1) = (1 - blend) * origCoefs(idx, 1) + blend * by;
                    }
                }
                T md, mxd, ma, en;
                internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, md, mxd, ma, en);
                if (md > bestScore) {
                    bestScore = md;
                    bestCoefs = coefs;
                    curEnergy = en;
                }
            }
            coefs = bestCoefs;
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, minDet, maxDet, minAbs, curEnergy);
            valid = (minDet > minDetThreshold) || (maxDet < -minDetThreshold);
        }

        // Barrier-guarded smoothing / Winslow energy improvement on strictly interior CPs ONLY
        for (int iter = 0; iter < 20 && n_u > 2 && n_v > 2; ++iter) {
            gsMatrix<T> prevCoefs = coefs;
            T prevMinDet = minDet, prevEnergy = curEnergy;
            bool prevValid = valid;

            for (int j = 1; j < n_v - 1; ++j) {
                for (int i = 1; i < n_u - 1; ++i) {
                    int idx   = i + j * n_u;
                    int left  = (i - 1) + j * n_u;
                    int right = (i + 1) + j * n_u;
                    int down  = i + (j - 1) * n_u;
                    int up    = i + (j + 1) * n_u;

                    T lx = T(0.25) * (prevCoefs(left, 0) + prevCoefs(right, 0) + prevCoefs(down, 0) + prevCoefs(up, 0));
                    T ly = T(0.25) * (prevCoefs(left, 1) + prevCoefs(right, 1) + prevCoefs(down, 1) + prevCoefs(up, 1));
                    coefs(idx, 0) = T(0.75) * prevCoefs(idx, 0) + T(0.25) * lx;
                    coefs(idx, 1) = T(0.75) * prevCoefs(idx, 1) + T(0.25) * ly;
                }
            }

            T newMinDet, newMaxDet, newMinAbs, newEnergy;
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 21, newMinDet, newMaxDet, newMinAbs, newEnergy);
            bool newValid = (newMinDet > minDetThreshold) || (newMaxDet < -minDetThreshold);

            bool accept = false;
            if (!prevValid) {
                accept = (newMinDet > prevMinDet + T(1e-12));
            } else {
                accept = newValid && (newEnergy <= prevEnergy * (T(1) + T(1e-6)));
            }

            if (accept) {
                minDet = newMinDet;
                maxDet = newMaxDet;
                curEnergy = newEnergy;
                valid = newValid;
            } else {
                coefs = prevCoefs;
                break;
            }
        }

        if (!valid) allOk = false;
    }

    mp.computeTopology(T(1e-5));
    return allOk;
}

template <class T>
void gsPatchSplitter<T>::ensureExactBoundaryConformity(gsMultiPatch<T>& mp)
{
    if (mp.nPatches() > 0) {
        mp.computeTopology(T(1e-5));
    }
}

// ============================================================================
// Diagnostics & Validation
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::isConforming() const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return false;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        if (mp.patch(p).domainDim() != 2) return false;
    }
    gsMatrix<T> uMat(1, 1), pt1, pt2;
    for (const auto& ifc : mp.interfaces()) {
        auto b1 = mp.patch(ifc.first().patch).boundary(ifc.first().side());
        auto b2 = mp.patch(ifc.second().patch).boundary(ifc.second().side());
        gsMatrix<T> pr1 = b1->parameterRange(), pr2 = b2->parameterRange();
        uMat(0, 0) = T(0.5) * (pr1(0, 0) + pr1(0, 1));
        b1->eval_into(uMat, pt1);
        uMat(0, 0) = T(0.5) * (pr2(0, 0) + pr2(0, 1));
        b2->eval_into(uMat, pt2);
        if (std::hypot(pt1(0, 0) - pt2(0, 0), pt1(1, 0) - pt2(1, 0)) > T(1e-4))
            return false;
    }
    return true;
}

template <class T>
bool gsPatchSplitter<T>::checkJacobians(T* min_det) const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) {
        if (min_det) *min_det = T(0);
        return false;
    }
    T mAbs = std::numeric_limits<T>::max();
    bool ok = true;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        T pAbs = 0;
        if (!internal_splitter::isPatchDetSingleSign(mp.patch(p), 21, T(1e-7), &pAbs))
            ok = false;
        if (pAbs < mAbs) mAbs = pAbs;
    }
    if (min_det) *min_det = ok ? mAbs : -mAbs;
    return ok;
}

template <class T>
bool gsPatchSplitter<T>::checkPatchCornersConvex(T* min_tri_area2) const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return false;
    T worst = std::numeric_limits<T>::max();
    bool ok = true;
    gsMatrix<T> u(2, 1), val;
    const T uv[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        gsPoint2D<T> c[4];
        for (int k = 0; k < 4; ++k) {
            u(0, 0) = uv[k][0]; u(1, 0) = uv[k][1];
            mp.patch(p).eval_into(u, val);
            c[k] = gsPoint2D<T>(val(0, 0), val(1, 0));
        }
        T triA = 0;
        if (!internal_splitter::isStrictlyConvexQuad(c, &triA)) ok = false;
        if (triA < worst) worst = triA;
    }
    if (min_tri_area2) *min_tri_area2 = worst;
    return ok;
}

template <class T>
bool gsPatchSplitter<T>::checkNoSelfIntersections() const
{
    const auto& mp = m_bestCandidate.multiPatch;
    size_t nP = mp.nPatches();
    if (nP == 0) return false;

    struct PatchAABB { T xMin, xMax, yMin, yMax; };
    std::vector<PatchAABB> boxes(nP);
    std::vector<std::vector<gsPoint2D<T>>> polys(nP);

    for (size_t p = 0; p < nP; ++p) {
        polys[p] = internal_splitter::samplePatchBoundaryPoly(mp.patch(p), 20);
        T xMin = std::numeric_limits<T>::max(), xMax = -std::numeric_limits<T>::max();
        T yMin = std::numeric_limits<T>::max(), yMax = -std::numeric_limits<T>::max();
        for (const auto& pt : polys[p]) {
            xMin = std::min(xMin, pt.x); xMax = std::max(xMax, pt.x);
            yMin = std::min(yMin, pt.y); yMax = std::max(yMax, pt.y);
        }
        boxes[p] = {xMin - T(1e-7), xMax + T(1e-7), yMin - T(1e-7), yMax + T(1e-7)};

        int m = static_cast<int>(polys[p].size());
        for (int i = 0; i < m; ++i) {
            for (int j = i + 2; j < m; ++j) {
                if (i == 0 && j == m - 1) continue;
                if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % m],
                                                          polys[p][j], polys[p][(j + 1) % m]))
                    return false;
            }
        }
    }

    gsMatrix<T> u(2, 9), val;
    const T pts[3] = {T(0.25), T(0.5), T(0.75)};
    int col = 0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            u(0, col) = pts[a]; u(1, col) = pts[b];
            ++col;
        }
    }
    std::vector<std::vector<gsPoint2D<T>>> innerPts(nP, std::vector<gsPoint2D<T>>(9));
    for (size_t p = 0; p < nP; ++p) {
        mp.patch(p).eval_into(u, val);
        for (int k = 0; k < 9; ++k) {
            innerPts[p][k] = gsPoint2D<T>(val(0, k), val(1, k));
        }
    }

    for (size_t p = 0; p < nP; ++p) {
        for (size_t q = p + 1; q < nP; ++q) {
            if (boxes[p].xMax < boxes[q].xMin || boxes[p].xMin > boxes[q].xMax ||
                boxes[p].yMax < boxes[q].yMin || boxes[p].yMin > boxes[q].yMax)
                continue;

            int mp1 = static_cast<int>(polys[p].size());
            int mq1 = static_cast<int>(polys[q].size());
            for (int i = 0; i < mp1; ++i) {
                for (int j = 0; j < mq1; ++j) {
                    if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % mp1],
                                                              polys[q][j], polys[q][(j + 1) % mq1]))
                        return false;
                }
            }
            for (int k = 0; k < 9; ++k) {
                if (internal_splitter::pointInPolyStrict(innerPts[p][k], polys[q]))
                    return false;
            }
        }
    }
    return true;
}

template <class T>
bool gsPatchSplitter<T>::validateMultiPatch(
    const gsPlanarDomain<T>& dom,
    const gsMultiPatch<T>& mp,
    T bdrTol,
    std::string* reason,
    T minAngleDeg,
    T maxAngleDeg)
{
    if (mp.nPatches() == 0) {
        if (reason) *reason = "empty multipatch";
        return false;
    }

    const bool hasTopo = (mp.nInterfaces() > 0 || mp.nBoundary() > 0) &&
                         (2 * mp.nInterfaces() + mp.nBoundary() == 4 * mp.nPatches());
    gsMultiPatch<T> mpTopo;
    const gsMultiPatch<T>& topoMp = hasTopo ? mp : mpTopo;
    if (!hasTopo) {
        mpTopo = mp;
        mpTopo.computeTopology(std::max(T(1e-4), bdrTol), false, false);
        if (2 * mpTopo.nInterfaces() + mpTopo.nBoundary() != 4 * mpTopo.nPatches()) {
            if (reason) *reason = "non-conforming topology (hanging nodes or unmatched sides)";
            return false;
        }
    }

    T bdrErr = 0;
    if (!verifyBoundaryMatch(dom, mp, bdrTol, &bdrErr)) {
        if (std::getenv("GISMO_SPLITTER_DEBUG")) {
            std::cerr << "[DBG verifyBoundaryMatch FAIL] bdrErr=" << bdrErr << " tol=" << bdrTol << "\n";
            for (size_t p = 0; p < topoMp.nPatches(); ++p) {
                for (int s = 1; s <= 4; ++s) {
                    boxSide bs = static_cast<boxSide>(s);
                    if (!topoMp.isInterface(patchSide(p, bs))) {
                        auto bc = topoMp.patch(p).boundary(bs);
                        gsMatrix<T> pr = bc->parameterRange(), u(1, 3), val;
                        u << pr(0,0), T(0.5)*(pr(0,0)+pr(0,1)), pr(0,1);
                        bc->eval_into(u, val);
                        std::cerr << "  bnd p=" << p << " s=" << s << ": (" << val(0,0) << "," << val(1,0) << ")->(" << val(0,1) << "," << val(1,1) << ")->(" << val(0,2) << "," << val(1,2) << ")\n";
                    }
                }
            }
        }
        if (reason) *reason = "boundary mismatch";
        return false;
    }

    gsMatrix<T> u(2, 1), val;
    const T uv[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    size_t nP = mp.nPatches();
    std::vector<std::vector<gsPoint2D<T>>> polys(nP);

    for (size_t p = 0; p < nP; ++p) {
        if (!internal_splitter::isPatchDetSingleSign(mp.patch(p), 21, T(1e-7))) {
            if (reason) *reason = "detJ zero or sign flip within patch";
            return false;
        }
        gsPoint2D<T> c[4];
        for (int k = 0; k < 4; ++k) {
            u(0, 0) = uv[k][0]; u(1, 0) = uv[k][1];
            mp.patch(p).eval_into(u, val);
            c[k] = gsPoint2D<T>(val(0, 0), val(1, 0));
        }
        for (int k = 0; k < 4; ++k) {
            int kn = (k + 1) % 4;
            if (std::hypot(c[kn].x - c[k].x, c[kn].y - c[k].y) < T(1e-8)) {
                if (reason) *reason = "non-convex or degenerate quad corners";
                return false;
            }
        }
        T angs[4], minA = 0, maxA = 180;
        bool degen = false;
        internal_splitter::evalPatchCornerAnglesDeg(mp.patch(p), angs, minA, maxA, degen);
        if (degen || minA < minAngleDeg || maxA > maxAngleDeg) {
            if (reason) *reason = "corner tangent angle out of admissible bounds";
            return false;
        }
        polys[p] = internal_splitter::samplePatchBoundaryPoly(mp.patch(p), 20);
        // Only apply straight-chord corner convexity to patches with near-linear boundaries
        // (curved spokes / curved hole arcs can validly curve past straight-chord collinearity).
        T maxEdgeDev = 0;
        {
            const int perSide = static_cast<int>(polys[p].size()) / 4;
            for (int k = 0; k < 4 && perSide > 1; ++k) {
                const gsPoint2D<T>& A = c[k];
                const gsPoint2D<T>& B = c[(k + 1) % 4];
                T len = std::hypot(B.x - A.x, B.y - A.y);
                if (len <= T(1e-12)) continue;
                for (int s = 1; s < perSide; ++s) {
                    const gsPoint2D<T>& q = polys[p][k * perSide + s];
                    T dist = std::abs((B.x - A.x) * (q.y - A.y) - (B.y - A.y) * (q.x - A.x)) / len;
                    maxEdgeDev = std::max(maxEdgeDev, dist);
                }
            }
        }
        if (maxEdgeDev < T(1e-4) && !internal_splitter::isStrictlyConvexQuad(c)) {
            if (reason) *reason = "non-convex or degenerate quad corners";
            return false;
        }
        int m = static_cast<int>(polys[p].size());
        for (int i = 0; i < m; ++i) {
            for (int j = i + 2; j < m; ++j) {
                if (i == 0 && j == m - 1) continue;
                if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % m],
                                                          polys[p][j], polys[p][(j + 1) % m])) {
                    if (reason) *reason = "patch boundary self-intersection";
                    return false;
                }
            }
        }
    }

    struct PatchAABB { T xMin, xMax, yMin, yMax; };
    std::vector<PatchAABB> boxes(nP);
    for (size_t p = 0; p < nP; ++p) {
        T xMin = std::numeric_limits<T>::max(), xMax = -std::numeric_limits<T>::max();
        T yMin = std::numeric_limits<T>::max(), yMax = -std::numeric_limits<T>::max();
        for (const auto& pt : polys[p]) {
            xMin = std::min(xMin, pt.x); xMax = std::max(xMax, pt.x);
            yMin = std::min(yMin, pt.y); yMax = std::max(yMax, pt.y);
        }
        boxes[p] = {xMin - T(1e-7), xMax + T(1e-7), yMin - T(1e-7), yMax + T(1e-7)};
    }

    gsMatrix<T> uInner(2, 9), valInner;
    const T pts[3] = {T(0.25), T(0.5), T(0.75)};
    int colIn = 0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            uInner(0, colIn) = pts[a]; uInner(1, colIn) = pts[b];
            ++colIn;
        }
    }
    std::vector<std::vector<gsPoint2D<T>>> innerPts(nP, std::vector<gsPoint2D<T>>(9));
    for (size_t p = 0; p < nP; ++p) {
        mp.patch(p).eval_into(uInner, valInner);
        for (int k = 0; k < 9; ++k) {
            innerPts[p][k] = gsPoint2D<T>(valInner(0, k), valInner(1, k));
        }
    }

    for (size_t p = 0; p < nP; ++p) {
        for (size_t q = p + 1; q < nP; ++q) {
            if (boxes[p].xMax < boxes[q].xMin || boxes[p].xMin > boxes[q].xMax ||
                boxes[p].yMax < boxes[q].yMin || boxes[p].yMin > boxes[q].yMax)
                continue;

            int mp1 = static_cast<int>(polys[p].size());
            int mq1 = static_cast<int>(polys[q].size());
            for (int i = 0; i < mp1; ++i) {
                for (int j = 0; j < mq1; ++j) {
                    if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % mp1],
                                                              polys[q][j], polys[q][(j + 1) % mq1])) {
                        if (reason) *reason = "cross-patch boundary intersection";
                        return false;
                    }
                }
            }
            for (int k = 0; k < 9; ++k) {
                if (internal_splitter::pointInPolyStrict(innerPts[p][k], polys[q])) {
                    if (reason) *reason = "patch overlap";
                    return false;
                }
            }
        }
    }

    return true;
}

template <class T>
bool gsPatchSplitter<T>::verifyBoundaryMatch(T tol, T* maxDistOut) const
{
    return verifyBoundaryMatch(m_domain, m_bestCandidate.multiPatch, tol, maxDistOut);
}

template <class T>
bool gsPatchSplitter<T>::verifyBoundaryMatch(
    const gsPlanarDomain<T>& dom,
    const gsMultiPatch<T>& mp,
    T tol,
    T* maxDistOut)
{
    if (mp.nPatches() == 0) return false;

    // 1. Collect all original domain curves
    std::vector<const gsGeometry<T>*> origCurves;
    const auto& outer = dom.outer();
    for (int i = 0; i < outer.size(); ++i) {
        origCurves.push_back(&outer.curve(i));
    }
    for (int h = 0; h < dom.numHoles(); ++h) {
        const auto& hole = dom.loop(h + 1);
        for (int j = 0; j < hole.size(); ++j) {
            origCurves.push_back(&hole.curve(j));
        }
    }

    // 2. Collect all MultiPatch exterior boundary curves
    std::vector<typename gsGeometry<T>::uPtr> mpBdrCurves;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        for (int s = 1; s <= 4; ++s) {
            boxSide side = static_cast<boxSide>(s);
            if (!mp.isInterface(patchSide(p, side))) {
                mpBdrCurves.push_back(mp.patch(p).boundary(side));
            }
        }
    }

    // Pre-sample curves once for fast closest-bracket lookup, then refine continuously via golden-section search
    struct SampledCurve {
        const gsGeometry<T>* crv;
        T u0, u1;
        std::vector<T> u, x, y;
        T xMin, xMax, yMin, yMax;
    };
    auto buildSampled = [](const gsGeometry<T>* crv, int nSeg = 40) -> SampledCurve {
        SampledCurve sc;
        sc.crv = crv;
        gsMatrix<T> pr = crv->parameterRange();
        sc.u0 = pr(0, 0); sc.u1 = pr(0, 1);
        sc.u.resize(nSeg + 1);
        sc.x.resize(nSeg + 1);
        sc.y.resize(nSeg + 1);
        gsMatrix<T> uMat(1, nSeg + 1), val;
        for (int i = 0; i <= nSeg; ++i) {
            T ui = sc.u0 + (sc.u1 - sc.u0) * (static_cast<T>(i) / static_cast<T>(nSeg));
            sc.u[i] = ui;
            uMat(0, i) = ui;
        }
        crv->eval_into(uMat, val);
        sc.xMin = std::numeric_limits<T>::max();
        sc.xMax = -std::numeric_limits<T>::max();
        sc.yMin = std::numeric_limits<T>::max();
        sc.yMax = -std::numeric_limits<T>::max();
        for (int i = 0; i <= nSeg; ++i) {
            sc.x[i] = val(0, i);
            sc.y[i] = val(1, i);
            sc.xMin = std::min(sc.xMin, sc.x[i]);
            sc.xMax = std::max(sc.xMax, sc.x[i]);
            sc.yMin = std::min(sc.yMin, sc.y[i]);
            sc.yMax = std::max(sc.yMax, sc.y[i]);
        }
        return sc;
    };

    std::vector<SampledCurve> origSampled, mpSampled;
    origSampled.reserve(origCurves.size());
    for (const auto* c : origCurves) origSampled.push_back(buildSampled(c));
    mpSampled.reserve(mpBdrCurves.size());
    for (const auto& c : mpBdrCurves) mpSampled.push_back(buildSampled(c.get()));

    auto aabbDistSq = [](T px, T py, const SampledCurve& sc) -> T {
        T dx = (px < sc.xMin) ? (sc.xMin - px) : ((px > sc.xMax) ? (px - sc.xMax) : T(0));
        T dy = (py < sc.yMin) ? (sc.yMin - py) : ((py > sc.yMax) ? (py - sc.yMax) : T(0));
        return dx * dx + dy * dy;
    };

    auto pointToCurvesDist = [&](T px, T py, const std::vector<SampledCurve>& targetCurves) -> T {
        T bestPolyD2 = std::numeric_limits<T>::max();
        int bestC = -1;
        int bestSeg = -1;

        for (size_t c = 0; c < targetCurves.size(); ++c) {
            const auto& sc = targetCurves[c];
            if (aabbDistSq(px, py, sc) >= bestPolyD2) continue;

            int nSeg = static_cast<int>(sc.u.size()) - 1;
            for (int i = 0; i < nSeg; ++i) {
                T ax = sc.x[i], ay = sc.y[i];
                T dx = sc.x[i + 1] - ax, dy = sc.y[i + 1] - ay;
                T L2 = dx * dx + dy * dy;
                T t = (L2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((px - ax) * dx + (py - ay) * dy) / L2)) : T(0);
                T qx = ax + t * dx - px, qy = ay + t * dy - py;
                T d2 = qx * qx + qy * qy;
                if (d2 < bestPolyD2) {
                    bestPolyD2 = d2;
                    bestC = static_cast<int>(c);
                    bestSeg = i;
                }
            }
        }

        if (bestC < 0) return std::numeric_limits<T>::max();

        // Refine on best curve via golden-section search
        const auto& sc = targetCurves[bestC];
        int nSeg = static_cast<int>(sc.u.size()) - 1;
        T du = (sc.u1 - sc.u0) / static_cast<T>(nSeg);
        T a = std::max(sc.u0, sc.u[bestSeg] - T(0.5) * du);
        T b = std::min(sc.u1, sc.u[bestSeg + 1] + T(0.5) * du);
        const T phi = T(0.6180339887498949);
        gsMatrix<T> uMat(1, 1), val;
        auto evalD2 = [&](T u) -> T {
            uMat(0, 0) = u;
            sc.crv->eval_into(uMat, val);
            T dx = val(0, 0) - px, dy = val(1, 0) - py;
            return dx * dx + dy * dy;
        };
        T c1 = b - phi * (b - a), c2 = a + phi * (b - a);
        T f1 = evalD2(c1), f2 = evalD2(c2);
        for (int iter = 0; iter < 18; ++iter) {
            if (f1 < f2) {
                b = c2; c2 = c1; f2 = f1;
                c1 = b - phi * (b - a); f1 = evalD2(c1);
            } else {
                a = c1; c1 = c2; f1 = f2;
                c2 = a + phi * (b - a); f2 = evalD2(c2);
            }
        }
        T bestFinal = std::min(bestPolyD2, std::min(f1, f2));
        return std::sqrt(bestFinal);
    };

    T maxDist = 0;

    // Forward check: every original boundary curve point must match MultiPatch boundary
    for (const auto& oSc : origSampled) {
        for (size_t i = 0; i < oSc.u.size(); ++i) {
            T d = pointToCurvesDist(oSc.x[i], oSc.y[i], mpSampled);
            if (d > maxDist) maxDist = d;
        }
    }

    // Backward check: every MultiPatch boundary point must match original domain curves
    for (const auto& bSc : mpSampled) {
        for (size_t i = 0; i < bSc.u.size(); ++i) {
            T d = pointToCurvesDist(bSc.x[i], bSc.y[i], origSampled);
            if (d > maxDist) maxDist = d;
        }
    }

    if (maxDistOut) *maxDistOut = maxDist;
    return (maxDist <= tol);
}

template <class T>
std::vector<int> gsPatchSplitter<T>::getVertexValences() const
{
    return analyzeQuality(m_bestCandidate.multiPatch).intValences;
}

// ============================================================================
// Export & Downstream Interface
// ============================================================================

template <class T>
const gsMultiPatch<T>& gsPatchSplitter<T>::getMultiPatch() const
{
    return m_bestCandidate.multiPatch;
}

template <class T>
gsMultiPatch<T> gsPatchSplitter<T>::getBilinearSkeleton() const
{
    gsMultiPatch<T> skeleton;
    const auto& mp = m_bestCandidate.multiPatch;

    gsKnotVector<T> kv_u(0.0, 1.0, 0, 2, 1, 1);
    gsKnotVector<T> kv_v(0.0, 1.0, 0, 2, 1, 1);

    for (size_t p = 0; p < mp.nPatches(); ++p) {
        const auto* tp = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;
        int nu = tp->basis().size(0);
        int nv = tp->basis().size(1);
        const auto& c = tp->coefs();

        gsMatrix<T> corner_coefs(4, 2);
        corner_coefs.row(0) = c.row(0);
        corner_coefs.row(1) = c.row(nu - 1);
        corner_coefs.row(2) = c.row((nv - 1) * nu);
        corner_coefs.row(3) = c.row(nu - 1 + (nv - 1) * nu);

        skeleton.addPatch(typename gsGeometry<T>::uPtr(
            new gsTensorBSpline<2, T>(kv_u, kv_v, corner_coefs)));
    }

    skeleton.computeTopology();
    return skeleton;
}

template <class T>
void gsPatchSplitter<T>::configureAsG1Generator(gsAsG1GeometryGenerator<T>& generator) const
{
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return;

    int V = static_cast<int>(mp.nPatches());
    generator.setValence(V);
    generator.setDegree(m_degree);

    std::vector<gsVector<T>> corners;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        const auto* tp = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;
        int nu = tp->basis().size(0);
        const auto& c = tp->coefs();
        gsVector<T> pt(2);
        pt << c(nu - 1, 0), c(nu - 1, 1);
        corners.push_back(pt);
    }

    generator.setBoundaryCorners(corners);
    generator.setUseVirtualParallelogram(true);
}

template <class T>
void gsPatchSplitter<T>::exportParaview(const std::string& filename_prefix, unsigned npts) const
{
    if (m_bestCandidate.multiPatch.nPatches() > 0) {
        gsWriteParaview(m_bestCandidate.multiPatch, filename_prefix, npts, true);
    }
}

} // namespace gismo
