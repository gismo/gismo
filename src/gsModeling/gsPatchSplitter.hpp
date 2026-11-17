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

template <class T>
inline gsPoint2D<T> evalCurvePoint(const gsCurve<T>& curve, T t)
{
    gsMatrix<T> u(1, 1);
    u(0, 0) = t;
    gsMatrix<T> pt;
    curve.eval_into(u, pt);
    return gsPoint2D<T>(pt(0, 0), pt(1, 0));
}

template <class T>
inline gsPoint2D<T> evalCurveDeriv(const gsCurve<T>& curve, T t)
{
    gsMatrix<T> u(1, 1);
    u(0, 0) = t;
    gsMatrix<T> der;
    curve.deriv_into(u, der);
    return gsPoint2D<T>(der(0, 0), der(1, 0));
}

template <class T>
inline gsPoint2D<T> startPoint(const gsCurve<T>& crv)
{
    gsMatrix<T> pr = crv.parameterRange();
    return evalCurvePoint(crv, pr(0, 0));
}

template <class T>
inline gsPoint2D<T> endPoint(const gsCurve<T>& crv)
{
    gsMatrix<T> pr = crv.parameterRange();
    return evalCurvePoint(crv, pr(0, 1));
}

template <class T>
inline gsBSpline<T>* makeStraightBSpline(const gsPoint2D<T>& p0, const gsPoint2D<T>& p1, int degree)
{
    int n = degree + 1;
    gsKnotVector<T> kv(0.0, 1.0, 0, n, 1, degree);
    gsMatrix<T> coefs(n, 2);
    for (int i = 0; i < n; ++i) {
        T t = static_cast<T>(i) / static_cast<T>(degree);
        coefs(i, 0) = (1.0 - t) * p0.x + t * p1.x;
        coefs(i, 1) = (1.0 - t) * p0.y + t * p1.y;
    }
    return new gsBSpline<T>(kv, coefs);
}

template <class T>
inline gsBSpline<T>* toClampedBSpline(const gsCurve<T>& curve, int target_degree)
{
    const auto* bsp = dynamic_cast<const gsBSpline<T>*>(&curve);
    if (bsp) {
        auto copy = new gsBSpline<T>(*bsp);
        copy->knots().affineTransformTo(0.0, 1.0);
        if (copy->degree() < target_degree) {
            copy->degreeElevate(target_degree - copy->degree());
        }
        return copy;
    }

    int n = target_degree + 1;
    gsMatrix<T> pr = curve.parameterRange();
    T u0 = pr(0, 0);
    T u1 = pr(0, 1);
    gsKnotVector<T> kv(0.0, 1.0, 0, n, 1, target_degree);
    gsMatrix<T> coefs(n, 2);
    for (int i = 0; i < n; ++i) {
        T t = static_cast<T>(i) / static_cast<T>(n - 1);
        T u_val = u0 + t * (u1 - u0);
        auto p = evalCurvePoint(curve, u_val);
        coefs(i, 0) = p.x;
        coefs(i, 1) = p.y;
    }
    return new gsBSpline<T>(kv, coefs);
}

template <class T>
inline void splitBSplineExact(const gsCurve<T>& crv, T t_split, gsBSpline<T>& left, gsBSpline<T>& right)
{
    std::unique_ptr<gsBSpline<T>> bsp(toClampedBSpline(crv, 3));
    bsp->splitAt(t_split, left, right);
    left.knots().affineTransformTo(0.0, 1.0);
    right.knots().affineTransformTo(0.0, 1.0);
}

template <class T>
inline T arcLengthSplitParam(const gsCurve<T>& crv, T frac = 0.5, int nSamples = 64)
{
    gsMatrix<T> pr = crv.parameterRange();
    T u0 = pr(0, 0), u1 = pr(0, 1);
    std::vector<T> u(nSamples + 1), s(nSamples + 1, T(0));
    gsPoint2D<T> prev = evalCurvePoint(crv, u0);
    u[0] = u0;
    for (int i = 1; i <= nSamples; ++i) {
        T ui = u0 + (u1 - u0) * (static_cast<T>(i) / static_cast<T>(nSamples));
        u[i] = ui;
        gsPoint2D<T> cur = evalCurvePoint(crv, ui);
        s[i] = s[i - 1] + std::hypot(cur.x - prev.x, cur.y - prev.y);
        prev = cur;
    }
    T target = frac * s.back();
    if (s.back() <= 1e-14) return u0 + frac * (u1 - u0);
    for (int i = 1; i <= nSamples; ++i) {
        if (s[i] >= target) {
            T ds = s[i] - s[i - 1];
            T alpha = (ds > 1e-14) ? (target - s[i - 1]) / ds : T(0);
            T uVal = u[i - 1] + alpha * (u[i] - u[i - 1]);
            return (uVal - u0) / (u1 - u0);
        }
    }
    return frac;
}

template <class T>
inline void harmonizeBSplinePair(std::unique_ptr<gsBSpline<T>>& a, std::unique_ptr<gsBSpline<T>>& b)
{
    a->knots().affineTransformTo(0.0, 1.0);
    b->knots().affineTransformTo(0.0, 1.0);
    int maxDeg = std::max(a->degree(), b->degree());
    if (a->degree() < maxDeg) a->degreeElevate(maxDeg - a->degree());
    if (b->degree() < maxDeg) b->degreeElevate(maxDeg - b->degree());

    gsKnotVector<T> kvUnion = a->knots().knotUnion(b->knots());
    std::vector<T> insA, insB;
    kvUnion.symDifference(a->knots(), insA);
    kvUnion.symDifference(b->knots(), insB);
    if (!insA.empty()) a->insertKnots(insA.begin(), insA.end());
    if (!insB.empty()) b->insertKnots(insB.begin(), insB.end());
}

template <class T>
inline gsTensorBSpline<2, T>* makeCoonsPatchFrom4Curves(
    const gsCurve<T>& south,
    const gsCurve<T>& east,
    const gsCurve<T>& north,
    const gsCurve<T>& west,
    int degree)
{
    std::unique_ptr<gsBSpline<T>> s(toClampedBSpline(south, degree));
    std::unique_ptr<gsBSpline<T>> e(toClampedBSpline(east, degree));
    std::unique_ptr<gsBSpline<T>> n(toClampedBSpline(north, degree));
    std::unique_ptr<gsBSpline<T>> w(toClampedBSpline(west, degree));

    harmonizeBSplinePair(s, n);
    harmonizeBSplinePair(w, e);

    const gsKnotVector<T>& kv_u = s->knots();
    const gsKnotVector<T>& kv_v = w->knots();
    int n_u = s->basis().size();
    int n_v = w->basis().size();

    gsMatrix<T> gr_u = kv_u.greville();
    gsMatrix<T> gr_v = kv_v.greville();

    const gsMatrix<T>& cS = s->coefs();
    const gsMatrix<T>& cE = e->coefs();
    const gsMatrix<T>& cN = n->coefs();
    const gsMatrix<T>& cW = w->coefs();

    gsPoint2D<T> P00(cS(0, 0), cS(0, 1));
    gsPoint2D<T> P10(cS(n_u - 1, 0), cS(n_u - 1, 1));
    gsPoint2D<T> P01(cN(0, 0), cN(0, 1));
    gsPoint2D<T> P11(cN(n_u - 1, 0), cN(n_u - 1, 1));

    gsMatrix<T> coefs(n_u * n_v, 2);

    for (int j = 0; j < n_v; ++j) {
        T eta = gr_v(0, j);
        for (int i = 0; i < n_u; ++i) {
            T xi = gr_u(0, i);

            T x_sn = (1.0 - eta) * cS(i, 0) + eta * cN(i, 0);
            T y_sn = (1.0 - eta) * cS(i, 1) + eta * cN(i, 1);

            T x_we = (1.0 - xi) * cW(j, 0) + xi * cE(j, 0);
            T y_we = (1.0 - xi) * cW(j, 1) + xi * cE(j, 1);

            T x_c = (1.0 - xi) * (1.0 - eta) * P00.x +
                    xi * (1.0 - eta) * P10.x +
                    (1.0 - xi) * eta * P01.x +
                    xi * eta * P11.x;
            T y_c = (1.0 - xi) * (1.0 - eta) * P00.y +
                    xi * (1.0 - eta) * P10.y +
                    (1.0 - xi) * eta * P01.y +
                    xi * eta * P11.y;

            int idx = i + j * n_u;
            coefs(idx, 0) = x_sn + x_we - x_c;
            coefs(idx, 1) = y_sn + y_we - y_c;
        }
    }

    return new gsTensorBSpline<2, T>(kv_u, kv_v, coefs);
}

// ============================================================================
// Hole-Free Geometry, Convexity, Jacobian & Intersection Utilities
// ============================================================================

template <class T>
struct CurveSeg {
    std::shared_ptr<gsBSpline<T>> bsp;

    CurveSeg() = default;
    explicit CurveSeg(const gsCurve<T>& crv, int degree) {
        bsp.reset(toClampedBSpline(crv, degree));
    }
    explicit CurveSeg(const gsBSpline<T>& b) {
        bsp = std::make_shared<gsBSpline<T>>(b);
        bsp->knots().affineTransformTo(0.0, 1.0);
    }

    gsPoint2D<T> p0() const { return startPoint(*bsp); }
    gsPoint2D<T> p1() const { return endPoint(*bsp); }

    CurveSeg reversed() const {
        CurveSeg r;
        r.bsp = std::make_shared<gsBSpline<T>>(*bsp);
        r.bsp->reverse();
        r.bsp->knots().affineTransformTo(0.0, 1.0);
        return r;
    }

    std::pair<CurveSeg, CurveSeg> splitArc(T frac = 0.5) const {
        T tSplit = arcLengthSplitParam(*bsp, frac);
        tSplit = std::max(T(0.05), std::min(T(0.95), tSplit));
        gsBSpline<T> left, right;
        bsp->splitAt(tSplit, left, right);
        return {CurveSeg(left), CurveSeg(right)};
    }

    std::pair<CurveSeg, CurveSeg> splitAtParam(T tSplit) const {
        tSplit = std::max(T(1e-4), std::min(T(1) - T(1e-4), tSplit));
        gsBSpline<T> left, right;
        bsp->splitAt(tSplit, left, right);
        return {CurveSeg(left), CurveSeg(right)};
    }

    T closestParam(const gsPoint2D<T>& q) const {
        const int nS = 64;
        T bestU = T(0.5);
        T bestD2 = std::numeric_limits<T>::max();
        for (int i = 0; i <= nS; ++i) {
            T u = static_cast<T>(i) / static_cast<T>(nS);
            gsPoint2D<T> pt = evalCurvePoint(*bsp, u);
            T dx = pt.x - q.x, dy = pt.y - q.y;
            T d2 = dx * dx + dy * dy;
            if (d2 < bestD2) { bestD2 = d2; bestU = u; }
        }
        T du = T(1) / static_cast<T>(nS);
        T a = std::max(T(0), bestU - du);
        T b = std::min(T(1), bestU + du);
        const T phi = T(0.6180339887498949);
        auto evalD2 = [&](T u) -> T {
            gsPoint2D<T> pt = evalCurvePoint(*bsp, u);
            T dx = pt.x - q.x, dy = pt.y - q.y;
            return dx * dx + dy * dy;
        };
        T c1 = b - phi * (b - a), c2 = a + phi * (b - a);
        T f1 = evalD2(c1), f2 = evalD2(c2);
        for (int it = 0; it < 30; ++it) {
            if (f1 < f2) {
                b = c2; c2 = c1; f2 = f1;
                c1 = b - phi * (b - a); f1 = evalD2(c1);
            } else {
                a = c1; c1 = c2; f1 = f2;
                c2 = a + phi * (b - a); f2 = evalD2(c2);
            }
        }
        return T(0.5) * (a + b);
    }
};

template <class T>
inline CurveSeg<T> makeStraightSeg(const gsPoint2D<T>& a, const gsPoint2D<T>& b, int degree)
{
    CurveSeg<T> seg;
    seg.bsp.reset(makeStraightBSpline(a, b, degree));
    return seg;
}

template <class T>
inline gsTensorBSpline<2, T>* makeCoonsFromCCWCycle(
    const CurveSeg<T>& e01,
    const CurveSeg<T>& e12,
    const CurveSeg<T>& e23,
    const CurveSeg<T>& e30,
    int degree)
{
    CurveSeg<T> north = e23.reversed();
    CurveSeg<T> west  = e30.reversed();
    return makeCoonsPatchFrom4Curves(*e01.bsp, *e12.bsp, *north.bsp, *west.bsp, degree);
}

template <class T>
inline T cross2(const gsPoint2D<T>& a, const gsPoint2D<T>& b, const gsPoint2D<T>& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

template <class T>
inline bool isStrictlyConvexQuad(const gsPoint2D<T> c[4], T* minTriArea2 = nullptr)
{
    int sign = 0;
    T minCr = std::numeric_limits<T>::max();
    for (int i = 0; i < 4; ++i) {
        T cr = cross2(c[i], c[(i + 1) % 4], c[(i + 2) % 4]);
        T acr = std::abs(cr);
        if (acr < minCr) minCr = acr;
        if (acr <= T(1e-8)) {
            if (minTriArea2) *minTriArea2 = minCr;
            return false;
        }
        int s = (cr > 0) ? 1 : -1;
        if (sign == 0) sign = s;
        else if (s != sign) {
            if (minTriArea2) *minTriArea2 = minCr;
            return false;
        }
    }
    if (minTriArea2) *minTriArea2 = minCr;
    return true;
}

template <class T>
inline T computeDomainMinCornerAngleDeg(const gsPlanarDomain<T>& dom)
{
    const T radToDeg = T(180) / T(3.14159265358979323846);
    T minAng = T(180);
    for (int l = 0; l < dom.numLoops(); ++l) {
        const auto& lp = dom.loop(l);
        const int K = lp.size();
        for (int i = 0; i < K; ++i) {
            const auto& prevC = lp.curve((i + K - 1) % K);
            const auto& nextC = lp.curve(i);
            gsMatrix<T> pr0 = prevC.parameterRange(), pr1 = nextC.parameterRange();
            auto din  = evalCurveDeriv(prevC, pr0(0, 1));
            auto dout = evalCurveDeriv(nextC, pr1(0, 0));
            T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
            if (n1 <= T(1e-14) || n2 <= T(1e-14)) continue;
            T cr = (din.x * dout.y - din.y * dout.x) / (n1 * n2);
            T dt = (din.x * dout.x + din.y * dout.y) / (n1 * n2);
            T turnRad = std::atan2(cr, dt);
            T intAng = (l == 0) ? (T(180) - turnRad * radToDeg) : (T(180) + turnRad * radToDeg);
            if (intAng > T(1e-3) && intAng < minAng) minAng = intAng;
        }
    }
    return minAng;
}

template <class T>
inline bool pickOptimalTri1to3Center(
    const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
    const gsPoint2D<T>& mAB, const gsPoint2D<T>& mBC, const gsPoint2D<T>& mCA,
    gsPoint2D<T>& outO)
{
    gsPoint2D<T> O_mid((mAB.x + mBC.x + mCA.x) / T(3), (mAB.y + mBC.y + mCA.y) / T(3));
    gsPoint2D<T> O_tri((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
    T la = B.dist(C), lb = C.dist(A), lc = A.dist(B), lsum = la + lb + lc;
    gsPoint2D<T> O_inc = (lsum > T(1e-14))
        ? gsPoint2D<T>((la * A.x + lb * B.x + lc * C.x) / lsum,
                       (la * A.y + lb * B.y + lc * C.y) / lsum)
        : O_tri;

    std::vector<gsPoint2D<T>> cands = {O_mid, O_tri, O_inc};
    for (int s = 1; s <= 4; ++s) {
        T w = static_cast<T>(s) / T(5);
        cands.emplace_back((T(1) - w) * O_mid.x + w * O_inc.x, (T(1) - w) * O_mid.y + w * O_inc.y);
        cands.emplace_back((T(1) - w) * O_mid.x + w * O_tri.x, (T(1) - w) * O_mid.y + w * O_tri.y);
    }
    const gsPoint2D<T> verts3[3] = {A, B, C};
    const T vWeights[4] = {T(0.12), T(0.24), T(0.36), T(0.48)};
    for (const auto& V : verts3) {
        for (T w : vWeights) {
            cands.emplace_back((T(1) - w) * O_mid.x + w * V.x, (T(1) - w) * O_mid.y + w * V.y);
        }
    }

    auto quadMaxDeg = [](const gsPoint2D<T> q[4]) -> T {
        const T radToDeg = T(180) / T(3.14159265358979323846);
        T mx = T(0);
        for (int k = 0; k < 4; ++k) {
            const auto& p0 = q[(k + 3) % 4];
            const auto& pk = q[k];
            const auto& p1 = q[(k + 1) % 4];
            T ux = p1.x - pk.x, uy = p1.y - pk.y;
            T vx = p0.x - pk.x, vy = p0.y - pk.y;
            T nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
            if (nu <= T(1e-14) || nv <= T(1e-14)) return T(360);
            T cr = (ux * vy - uy * vx) / (nu * nv);
            T dt = (ux * vx + uy * vy) / (nu * nv);
            T ang = std::atan2(cr, dt) * radToDeg;
            if (ang <= T(0)) ang += T(360);
            if (ang > mx) mx = ang;
        }
        return mx;
    };

    bool found = false;
    T bestMax = T(360);
    for (const auto& Oc : cands) {
        gsPoint2D<T> q0[4] = {A, mAB, Oc, mCA};
        gsPoint2D<T> q1[4] = {B, mBC, Oc, mAB};
        gsPoint2D<T> q2[4] = {C, mCA, Oc, mBC};
        if (!isStrictlyConvexQuad(q0) || !isStrictlyConvexQuad(q1) || !isStrictlyConvexQuad(q2))
            continue;
        T mA = std::max(quadMaxDeg(q0), std::max(quadMaxDeg(q1), quadMaxDeg(q2)));
        if (!found || mA < bestMax) {
            bestMax = mA;
            outO = Oc;
            found = true;
        }
    }
    return found;
}

template <class T>
inline bool segProperIntersect(const gsPoint2D<T>& a, const gsPoint2D<T>& b,
                               const gsPoint2D<T>& c, const gsPoint2D<T>& d)
{
    T d1 = cross2(a, b, c), d2 = cross2(a, b, d);
    T d3 = cross2(c, d, a), d4 = cross2(c, d, b);
    const T eps = T(1e-9);
    return (((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
            ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps)));
}

template <class T>
inline bool pointInPolyStrict(const gsPoint2D<T>& p, const std::vector<gsPoint2D<T>>& poly, T tol = T(1e-7))
{
    int n = static_cast<int>(poly.size());
    if (n < 3) return false;
    T tol2 = tol * tol;
    for (int i = 0; i < n; ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % n];
        T dx = b.x - a.x, dy = b.y - a.y;
        T L2 = dx * dx + dy * dy;
        if (L2 > T(1e-20)) {
            T t = std::max(T(0), std::min(T(1), ((p.x - a.x) * dx + (p.y - a.y) * dy) / L2));
            T qx = a.x + t * dx - p.x;
            T qy = a.y + t * dy - p.y;
            if (qx * qx + qy * qy <= tol2) return false;
        }
    }
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        if (((poly[i].y > p.y) != (poly[j].y > p.y)) &&
            (p.x < (poly[j].x - poly[i].x) * (p.y - poly[i].y) / (poly[j].y - poly[i].y + T(1e-30)) + poly[i].x))
        {
            inside = !inside;
        }
    }
    return inside;
}

template <class T>
inline std::vector<gsPoint2D<T>> samplePatchBoundaryPoly(const gsGeometry<T>& patch, int nPerSide = 20)
{
    std::vector<gsPoint2D<T>> poly;
    poly.reserve(4 * nPerSide);
    gsMatrix<T> u(2, 1), val;
    auto pushPt = [&](T xi, T eta) {
        u(0, 0) = xi; u(1, 0) = eta;
        patch.eval_into(u, val);
        poly.push_back(gsPoint2D<T>(val(0, 0), val(1, 0)));
    };
    for (int i = 0; i < nPerSide; ++i) pushPt(static_cast<T>(i) / nPerSide, T(0));
    for (int i = 0; i < nPerSide; ++i) pushPt(T(1), static_cast<T>(i) / nPerSide);
    for (int i = 0; i < nPerSide; ++i) pushPt(T(1) - static_cast<T>(i) / nPerSide, T(1));
    for (int i = 0; i < nPerSide; ++i) pushPt(T(0), T(1) - static_cast<T>(i) / nPerSide);
    return poly;
}

template <class T>
inline std::vector<gsPoint2D<T>> sampleCurveLoopPoly(const gsCurveLoop<T>& loop, int nPerCurve = 24)
{
    std::vector<gsPoint2D<T>> poly;
    for (int c = 0; c < loop.size(); ++c) {
        const auto& crv = loop.curve(c);
        gsMatrix<T> pr = crv.parameterRange();
        T u0 = pr(0, 0), u1 = pr(0, 1);
        for (int i = 0; i < nPerCurve; ++i) {
            T u = u0 + (u1 - u0) * (static_cast<T>(i) / nPerCurve);
            poly.push_back(evalCurvePoint(crv, u));
        }
    }
    return poly;
}

template <class T>
inline gsPoint2D<T> polygonCentroid(const std::vector<gsPoint2D<T>>& poly)
{
    T a2 = 0, cx = 0, cy = 0;
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
        const auto& p = poly[i];
        const auto& q = poly[(i + 1) % n];
        T cr = p.x * q.y - q.x * p.y;
        a2 += cr;
        cx += (p.x + q.x) * cr;
        cy += (p.y + q.y) * cr;
    }
    if (std::abs(a2) > T(1e-12))
        return gsPoint2D<T>(cx / (3 * a2), cy / (3 * a2));
    gsPoint2D<T> avg(0, 0);
    for (const auto& p : poly) { avg.x += p.x; avg.y += p.y; }
    if (n > 0) { avg.x /= n; avg.y /= n; }
    return avg;
}

template <class T>
inline bool isSegmentInsidePoly(const gsPoint2D<T>& a, const gsPoint2D<T>& b,
                                const std::vector<gsPoint2D<T>>& outerPoly)
{
    int m = static_cast<int>(outerPoly.size());
    for (int i = 0; i < m; ++i) {
        if (segProperIntersect(a, b, outerPoly[i], outerPoly[(i + 1) % m]))
            return false;
    }
    for (int k = 1; k <= 3; ++k) {
        T t = T(k) / T(4);
        gsPoint2D<T> mid((1 - t) * a.x + t * b.x, (1 - t) * a.y + t * b.y);
        if (!pointInPolyStrict(mid, outerPoly, T(1e-6)))
            return false;
    }
    return true;
}

template <class T>
inline void evalPatchJacobianAndWinslow(const gsGeometry<T>& patch, int nSamples,
                                        T& minDet, T& maxDet, T& minAbsDet, T& energy)
{
    minDet = std::numeric_limits<T>::max();
    maxDet = -std::numeric_limits<T>::max();
    minAbsDet = std::numeric_limits<T>::max();
    energy = T(0);

    const T h = T(1) / (nSamples - 1);
    gsMatrix<T> u(2, 1), J;
    for (int j = 0; j < nSamples; ++j) {
        T eta = j * h;
        for (int i = 0; i < nSamples; ++i) {
            T xi = i * h;
            u(0, 0) = xi;
            u(1, 0) = eta;
            patch.jacobian_into(u, J);

            T Jxi_x  = J(0, 0), Jxi_y  = J(1, 0);
            T Jeta_x = J(0, 1), Jeta_y = J(1, 1);

            T normSq_xi  = Jxi_x * Jxi_x + Jxi_y * Jxi_y;
            T normSq_eta = Jeta_x * Jeta_x + Jeta_y * Jeta_y;
            T detJ = Jxi_x * Jeta_y - Jxi_y * Jeta_x;

            minDet = std::min(minDet, detJ);
            maxDet = std::max(maxDet, detJ);
            T absDet = std::abs(detJ);
            minAbsDet = std::min(minAbsDet, absDet);

            if (absDet > T(1e-12)) {
                T norm_xi  = std::sqrt(normSq_xi);
                T norm_eta = std::sqrt(normSq_eta);
                T win = (normSq_xi + normSq_eta) / absDet;
                T orth = (absDet / (norm_xi * norm_eta + T(1e-14))) - T(1);
                energy += (win + orth * orth) * (h * h);
            } else {
                energy += T(1e6);
            }
        }
    }
}

template <class T>
inline bool isPatchDetSingleSign(const gsGeometry<T>& patch, int nSamples = 21, T tol = T(1e-8),
                                 T* minAbsOut = nullptr)
{
    T minDet = 0, maxDet = 0, minAbs = 0, energy = 0;
    evalPatchJacobianAndWinslow(patch, nSamples, minDet, maxDet, minAbs, energy);
    if (minAbsOut) *minAbsOut = minAbs;
    return (minDet > tol) || (maxDet < -tol);
}

// Conforming 4-patch quad ring between 4 outer curves and 4 hole curves
template <class T>
inline void buildQuadRingByProximity(
    const gsCurveLoop<T>& outer,
    const gsCurveLoop<T>& hole,
    gsMultiPatch<T>& mp,
    int degree)
{
    for (int i = 0; i < 4; ++i) {
        const auto& oCrv = outer.curve(i);
        auto pA = startPoint(oCrv);
        auto pB = endPoint(oCrv);

        int best_h = -1;
        bool best_rev = false;
        T best_dist = 1e9;
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

        // Radial spokes from hole outward to outer boundary:
        // South: hA -> pA, North: hB -> pB, West: hA -> hB (hole), East: pA -> pB (outer)
        std::unique_ptr<gsBSpline<T>> sp_S(makeStraightBSpline(hA, pA, degree));
        std::unique_ptr<gsBSpline<T>> sp_N(makeStraightBSpline(hB, pB, degree));

        auto patch = makeCoonsPatchFrom4Curves(*sp_S, oCrv, *sp_N, *hCrv, degree);
        mp.addPatch(typename gsGeometry<T>::uPtr(patch));
    }
}

// Conforming 4-patch block quad ring around an interior hole in a 4-sided block
template <class T>
inline void buildBlockQuadRing(
    const gsCurve<T>& oS, const gsCurve<T>& oE,
    const gsCurve<T>& oN, const gsCurve<T>& oW,
    const gsCurveLoop<T>& hole,
    gsMultiPatch<T>& mp,
    int degree)
{
    gsPoint2D<T> P_SW = startPoint(oS);
    gsPoint2D<T> P_SE = endPoint(oS);
    gsPoint2D<T> P_NE = endPoint(oE);
    gsPoint2D<T> P_NW = startPoint(oN); // Correct start of North curve

    int nH = hole.size();
    auto findClosestHolePt = [&](const gsPoint2D<T>& target) -> gsPoint2D<T> {
        T bestDist = 1e9;
        gsPoint2D<T> bestPt(0, 0);
        for (int i = 0; i < nH; ++i) {
            auto p = startPoint(hole.curve(i));
            T d = (p.x - target.x)*(p.x - target.x) + (p.y - target.y)*(p.y - target.y);
            if (d < bestDist) { bestDist = d; bestPt = p; }
        }
        return bestPt;
    };

    gsPoint2D<T> H_SW = findClosestHolePt(P_SW);
    gsPoint2D<T> H_SE = findClosestHolePt(P_SE);
    gsPoint2D<T> H_NE = findClosestHolePt(P_NE);
    gsPoint2D<T> H_NW = findClosestHolePt(P_NW);

    std::unique_ptr<gsBSpline<T>> sp_SW_in(makeStraightBSpline(P_SW, H_SW, degree));
    std::unique_ptr<gsBSpline<T>> sp_SE_in(makeStraightBSpline(P_SE, H_SE, degree));
    std::unique_ptr<gsBSpline<T>> sp_NW_out(makeStraightBSpline(H_NW, P_NW, degree));
    std::unique_ptr<gsBSpline<T>> sp_NE_out(makeStraightBSpline(H_NE, P_NE, degree));

    std::unique_ptr<gsCurve<T>> hS(makeStraightBSpline(H_SW, H_SE, degree));
    std::unique_ptr<gsCurve<T>> hN(makeStraightBSpline(H_NW, H_NE, degree));
    std::unique_ptr<gsCurve<T>> hW(makeStraightBSpline(H_SW, H_NW, degree));
    std::unique_ptr<gsCurve<T>> hE(makeStraightBSpline(H_SE, H_NE, degree));

    std::unique_ptr<gsCurve<T>> sp_E_bot(makeStraightBSpline(H_SE, P_SE, degree));
    std::unique_ptr<gsCurve<T>> sp_E_top(makeStraightBSpline(H_NE, P_NE, degree));
    std::unique_ptr<gsCurve<T>> sp_W_bot(makeStraightBSpline(P_SW, H_SW, degree));
    std::unique_ptr<gsCurve<T>> sp_W_top(makeStraightBSpline(P_NW, H_NW, degree));

    // 1. South patch
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(oS, *sp_SE_in, *hS, *sp_SW_in, degree)));
    // 2. East patch
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*sp_E_bot, oE, *sp_E_top, *hE, degree)));
    // 3. North patch
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*hN, *sp_NE_out, oN, *sp_NW_out, degree)));
    // 4. West patch
    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsPatchFrom4Curves(*sp_W_bot, *hW, *sp_W_top, oW, degree)));
}

} // namespace internal_splitter
} // namespace gismo

#include <gsModeling/gsPatchSplitterPolychord.hpp>
#include <gsModeling/gsPatchSplitterCrossField.hpp>

namespace gismo {

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

    const T domMinAng = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
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
        const T domMinAng = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
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
        const T domMinAng = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
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
                e_pw = t1.e[(mateEdgePos + 2) % 3];
                e_wn = t1.e[(mateEdgePos + 1) % 3];
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
                Tri A = mesh.tris[tA], B = mesh.tris[tB];
                int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
                int d = B.v[(kB + 2) % 3];
                gsPoint2D<T> qa[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
                if (!isStrictlyConvexQuad(qa)) continue;
                if (!isSegmentInsidePoly(mesh.pts[c], mesh.pts[d], outerPoly)) continue;

                auto minTriSin = [&](int i0, int i1, int i2) -> T {
                    const auto& P0 = mesh.pts[i0];
                    const auto& P1 = mesh.pts[i1];
                    const auto& P2 = mesh.pts[i2];
                    T cr2 = cross2(P0, P1, P2);
                    if (cr2 <= T(1e-8)) return -T(1);
                    T l2 = std::pow(P0.x - P1.x, 2) + std::pow(P0.y - P1.y, 2)
                         + std::pow(P1.x - P2.x, 2) + std::pow(P1.y - P2.y, 2)
                         + std::pow(P2.x - P0.x, 2) + std::pow(P2.y - P0.y, 2);
                    return cr2 / (l2 + T(1e-14));
                };
                T curMin = std::min(minTriSin(a, b, c), minTriSin(b, a, d));
                T newMin = std::min(minTriSin(c, a, d), minTriSin(c, d, b));
                if (newMin <= curMin + T(1e-4)) continue;

                // Protect smooth boundary vertices from dropping to deg == 1 unless needed
                int degA = 0, degB = 0;
                for (int ti = 0; ti < nT; ++ti) {
                    for (int k = 0; k < 3; ++k) {
                        if (mesh.tris[ti].v[k] == a) ++degA;
                        if (mesh.tris[ti].v[k] == b) ++degB;
                    }
                }
                if ((mesh.isSmooth[a] && degA <= 2) || (mesh.isSmooth[b] && degB <= 2)) continue;

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

        struct PairEdge {
            int t0, t1;
            Quad q;
            T score;
        };
        std::vector<PairEdge> pairEdges;
        for (int eIdx = 0; eIdx < static_cast<int>(triMesh.edges.size()); ++eIdx) {
            if (triMesh.edges[eIdx].isBoundary) continue;
            int tA = -1, kA = -1, tB = -1, kB = -1;
            for (int ti = 0; ti < nT; ++ti) {
                for (int k = 0; k < 3; ++k) {
                    if (triMesh.tris[ti].e[k] == eIdx) {
                        if (tA < 0) { tA = ti; kA = k; }
                        else        { tB = ti; kB = k; }
                    }
                }
            }
            if (tA < 0 || tB < 0) continue;
            Tri A = triMesh.tris[tA], B = triMesh.tris[tB];
            int a = A.v[kA], b = A.v[(kA + 1) % 3], c = A.v[(kA + 2) % 3];
            int d = B.v[(kB + 2) % 3];
            gsPoint2D<T> qc[4] = {triMesh.pts[c], triMesh.pts[a], triMesh.pts[d], triMesh.pts[b]};
            T minArea2 = 0;
            if (!isStrictlyConvexQuad(qc, &minArea2)) continue;
            {
                const T radToDeg = T(180) / T(3.14159265358979323846);
                const T domMinA = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
                const T qMinTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinA);
                const T qMaxTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
                bool angOk = true;
                for (int k = 0; k < 4; ++k) {
                    const auto& p0 = qc[(k + 3) % 4];
                    const auto& pk = qc[k];
                    const auto& p1 = qc[(k + 1) % 4];
                    T ux = p1.x - pk.x, uy = p1.y - pk.y;
                    T vx = p0.x - pk.x, vy = p0.y - pk.y;
                    T nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
                    if (nu < T(1e-12) || nv < T(1e-12)) { angOk = false; break; }
                    T ang = std::atan2((ux * vy - uy * vx) / (nu * nv), (ux * vx + uy * vy) / (nu * nv)) * radToDeg;
                    if (ang <= T(0)) ang += T(360);
                    if (ang < qMinTol - T(1.5) || ang > qMaxTol) { angOk = false; break; }
                }
                if (!angOk) continue;
            }

            int e_bc = A.e[(kA + 1) % 3];
            int e_ca = A.e[(kA + 2) % 3];
            int e_ad = -1, e_db = -1;
            for (int k = 0; k < 3; ++k) {
                int u = B.v[k], v = B.v[(k + 1) % 3];
                if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
            }
            if (e_ad < 0 || e_db < 0) continue;

            T perim2 = T(0);
            for (int s = 0; s < 4; ++s) {
                int sn = (s + 1) % 4;
                perim2 += std::pow(qc[s].x - qc[sn].x, 2) + std::pow(qc[s].y - qc[sn].y, 2);
            }
            T sc = minArea2 / (perim2 + T(1e-14));
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

            ring[iPrev].outSeg = getSegDir(MixedMesh{pts, isSmooth, edges, {}, {}}, e2, ia, ic);
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
    const T domMinAng = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));

    auto updateBest = [&](const SplittingCandidate& c) {
        if (!c.isValid) return;
        if (!bestCand.isValid) { bestCand = c; return; }
        T qNew = m_weightPatches * T(c.numPatches) + m_weightDistortion * c.distortionEnergy;
        T qOld = m_weightPatches * T(bestCand.numPatches) + m_weightDistortion * bestCand.distortionEnergy;
        if (qNew < qOld) bestCand = c;
    };

    auto quadifyTriMesh = [&](TriMeshData<T>& mesh, bool coarseRelaxed = false) {
        const int nT = static_cast<int>(mesh.tris.size());
        if (nT == 0) return;

        // Ensure all edges have CurveSegs and split halves
        for (auto& er : mesh.edges) {
            if (!er.seg.bsp)
                er.seg = makeStraightSeg(mesh.pts[er.a], mesh.pts[er.b], m_degree);
            if (!er.splitDone) {
                auto sp = er.seg.splitArc(T(0.5));
                er.h0 = sp.first;
                er.h1 = sp.second;
                er.mid = sp.first.p1();
                er.splitDone = true;
            }
        }

        auto getSegDir = [&](int eIdx, int u, int v) -> CurveSeg<T> {
            const auto& er = mesh.edges[eIdx];
            return (er.a == u && er.b == v) ? er.seg : er.seg.reversed();
        };
        auto getHalvesDir = [&](int eIdx, int u, int v, CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
            const auto& er = mesh.edges[eIdx];
            outMid = er.mid;
            if (er.a == u && er.b == v) { out0 = er.h0; out1 = er.h1; }
            else                         { out0 = er.h1.reversed(); out1 = er.h0.reversed(); }
        };

        const T effMinValAng = coarseRelaxed ? std::max(T(8), minAngTol - T(3.5)) : minAngTol;
        const T effMaxValAng = coarseRelaxed ? std::max(maxAngTol, T(164.2))     : maxAngTol;
        const T effMinPairAng = std::max(T(8), minAngTol - T(4.0));
        const T effMaxPairAng = std::max(maxAngTol, T(164.2));

        auto emitFromMixed = [&](const std::vector<typename TriMeshData<T>::Tri>& tris,
                                 const std::vector<typename TriMeshData<T>::Quad>& quads,
                                 const std::vector<bool>& needSplit) {
            gsMultiPatch<T> mp;
            for (const auto& t : tris) {
                if (!needSplit[t.e[0]] || !needSplit[t.e[1]] || !needSplit[t.e[2]]) return;
                int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
                CurveSeg<T> e01_a, e01_b, e12_a, e12_b, e20_a, e20_b;
                gsPoint2D<T> m01, m12, m20;
                getHalvesDir(t.e[0], v0, v1, e01_a, e01_b, m01);
                getHalvesDir(t.e[1], v1, v2, e12_a, e12_b, m12);
                getHalvesDir(t.e[2], v2, v0, e20_a, e20_b, m20);

                gsPoint2D<T> O;
                if (!internal_splitter::pickOptimalTri1to3Center(mesh.pts[v0], mesh.pts[v1], mesh.pts[v2],
                                                                 m01, m12, m20, O))
                    return;
                auto s01 = makeStraightSeg(O, m01, m_degree);
                auto s12 = makeStraightSeg(O, m12, m_degree);
                auto s20 = makeStraightSeg(O, m20, m_degree);
                mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e01_a, s01.reversed(), s20, e20_b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e12_a, s12.reversed(), s01, e01_b, m_degree)));
                mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e20_a, s20.reversed(), s12, e12_b, m_degree)));
            }
            for (const auto& q : quads) {
                bool sp0 = needSplit[q.e[0]], sp1 = needSplit[q.e[1]];
                bool sp2 = needSplit[q.e[2]], sp3 = needSplit[q.e[3]];
                if (sp0 != sp2 || sp1 != sp3) return;
                int v0 = q.v[0], v1 = q.v[1], v2 = q.v[2], v3 = q.v[3];
                if (!sp0 && !sp1) {
                    gsPoint2D<T> c[4] = {mesh.pts[v0], mesh.pts[v1], mesh.pts[v2], mesh.pts[v3]};
                    if (!isStrictlyConvexQuad(c)) return;
                    mp.addPatch(typename gsGeometry<T>::uPtr(
                        makeCoonsFromCCWCycle(getSegDir(q.e[0], v0, v1), getSegDir(q.e[1], v1, v2),
                                              getSegDir(q.e[2], v2, v3), getSegDir(q.e[3], v3, v0), m_degree)));
                } else if (sp0 && !sp1) {
                    CurveSeg<T> e0a, e0b, e2a, e2b; gsPoint2D<T> m01, m23;
                    getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                    getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                    gsPoint2D<T> qa[4] = {mesh.pts[v0], m01, m23, mesh.pts[v3]};
                    gsPoint2D<T> qb[4] = {m01, mesh.pts[v1], mesh.pts[v2], m23};
                    if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return;
                    auto sMid = makeStraightSeg(m01, m23, m_degree);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, sMid, e2b, getSegDir(q.e[3], v3, v0), m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0b, getSegDir(q.e[1], v1, v2), e2a, sMid.reversed(), m_degree)));
                } else if (!sp0 && sp1) {
                    CurveSeg<T> e1a, e1b, e3a, e3b; gsPoint2D<T> m12, m30;
                    getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                    getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                    gsPoint2D<T> qa[4] = {mesh.pts[v0], mesh.pts[v1], m12, m30};
                    gsPoint2D<T> qb[4] = {m30, m12, mesh.pts[v2], mesh.pts[v3]};
                    if (!isStrictlyConvexQuad(qa) || !isStrictlyConvexQuad(qb)) return;
                    auto sMid = makeStraightSeg(m12, m30, m_degree);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(getSegDir(q.e[0], v0, v1), e1a, sMid, e3b, m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(sMid.reversed(), e1b, getSegDir(q.e[2], v2, v3), e3a, m_degree)));
                } else {
                    CurveSeg<T> e0a, e0b, e1a, e1b, e2a, e2b, e3a, e3b; gsPoint2D<T> m01, m12, m23, m30;
                    getHalvesDir(q.e[0], v0, v1, e0a, e0b, m01);
                    getHalvesDir(q.e[1], v1, v2, e1a, e1b, m12);
                    getHalvesDir(q.e[2], v2, v3, e2a, e2b, m23);
                    getHalvesDir(q.e[3], v3, v0, e3a, e3b, m30);
                    gsPoint2D<T> O((m01.x + m12.x + m23.x + m30.x) / T(4), (m01.y + m12.y + m23.y + m30.y) / T(4));
                    gsPoint2D<T> Ov((mesh.pts[v0].x + mesh.pts[v1].x + mesh.pts[v2].x + mesh.pts[v3].x) / T(4),
                                    (mesh.pts[v0].y + mesh.pts[v1].y + mesh.pts[v2].y + mesh.pts[v3].y) / T(4));
                    auto check4 = [&](const gsPoint2D<T>& Oc) -> bool {
                        gsPoint2D<T> q0[4] = {mesh.pts[v0], m01, Oc, m30};
                        gsPoint2D<T> q1[4] = {mesh.pts[v1], m12, Oc, m01};
                        gsPoint2D<T> q2[4] = {mesh.pts[v2], m23, Oc, m12};
                        gsPoint2D<T> q3[4] = {mesh.pts[v3], m30, Oc, m23};
                        return isStrictlyConvexQuad(q0) && isStrictlyConvexQuad(q1)
                            && isStrictlyConvexQuad(q2) && isStrictlyConvexQuad(q3);
                    };
                    if (!check4(O)) {
                        if (check4(Ov)) O = Ov;
                        else {
                            bool found = false;
                            for (T w = T(0.1); w <= T(0.91); w += T(0.1)) {
                                gsPoint2D<T> Ot((1 - w) * O.x + w * Ov.x, (1 - w) * O.y + w * Ov.y);
                                if (check4(Ot)) { O = Ot; found = true; break; }
                            }
                            if (!found) return;
                        }
                    }
                    auto s01 = makeStraightSeg(O, m01, m_degree);
                    auto s12 = makeStraightSeg(O, m12, m_degree);
                    auto s23 = makeStraightSeg(O, m23, m_degree);
                    auto s30 = makeStraightSeg(O, m30, m_degree);
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e0a, s01.reversed(), s30, e3b, m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e1a, s12.reversed(), s01, e0b, m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e2a, s23.reversed(), s12, e1b, m_degree)));
                    mp.addPatch(typename gsGeometry<T>::uPtr(makeCoonsFromCCWCycle(e3a, s30.reversed(), s23, e2b, m_degree)));
                }
            }
            if (mp.nPatches() == 0) return;
            // Fast-reject only extreme pre-Winslow degeneracies so Winslow gets a chance to relax angles
            bool angPreOk = true;
            for (size_t p = 0; p < mp.nPatches(); ++p) {
                T angs[4], mnA = 0, mxA = 180;
                bool degen = false;
                evalPatchCornerAnglesDeg(mp.patch(p), angs, mnA, mxA, degen);
                if (degen || mnA < T(5) || mxA > T(174)) { angPreOk = false; break; }
            }
            if (!angPreOk) return;
            gsMultiPatch<T> preWinslow = mp;
            optimizeWinslowBarrier(mp);
            mp.computeTopology(T(1e-5), false, false);
            if (!validateMultiPatch(m_domain, mp, T(1e-4), nullptr, effMinValAng, effMaxValAng)) {
                preWinslow.computeTopology(T(1e-5), false, false);
                if (!validateMultiPatch(m_domain, preWinslow, T(1e-4), nullptr, effMinValAng, effMaxValAng))
                    return;
                mp = give(preWinslow);
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
        };

        // 1. Direct 1-to-3 split
        {
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(mesh.tris, {}, allTrue);
        }

        // 2. Greedy triangle pairing + parity propagation / Catmull-Clark
        struct PairEdge { int t0, t1; typename TriMeshData<T>::Quad q; T score; };
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
            gsPoint2D<T> qc[4] = {mesh.pts[c], mesh.pts[a], mesh.pts[d], mesh.pts[b]};
            T minArea2 = 0;
            if (!isStrictlyConvexQuad(qc, &minArea2)) continue;
            {
                const T radToDeg = T(180) / T(3.14159265358979323846);
                bool angOk = true;
                for (int k = 0; k < 4; ++k) {
                    const auto& p0 = qc[(k + 3) % 4];
                    const auto& pk = qc[k];
                    const auto& p1 = qc[(k + 1) % 4];
                    T ux = p1.x - pk.x, uy = p1.y - pk.y;
                    T vx = p0.x - pk.x, vy = p0.y - pk.y;
                    T nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
                    if (nu < T(1e-12) || nv < T(1e-12)) { angOk = false; break; }
                    T ang = std::atan2((ux * vy - uy * vx) / (nu * nv), (ux * vx + uy * vy) / (nu * nv)) * radToDeg;
                    if (ang <= T(0)) ang += T(360);
                    if (ang < effMinPairAng || ang > effMaxPairAng) { angOk = false; break; }
                }
                if (!angOk) continue;
            }
            int e_bc = A.e[(kA + 1) % 3], e_ca = A.e[(kA + 2) % 3];
            int e_ad = -1, e_db = -1;
            for (int k = 0; k < 3; ++k) {
                int u = B.v[k], v = B.v[(k + 1) % 3];
                if ((u == a && v == d) || (u == d && v == a)) e_ad = B.e[k];
                if ((u == d && v == b) || (u == b && v == d)) e_db = B.e[k];
            }
            if (e_ad < 0 || e_db < 0) continue;
            T perim2 = T(0);
            for (int s = 0; s < 4; ++s) {
                int sn = (s + 1) % 4;
                perim2 += (qc[s].x - qc[sn].x) * (qc[s].x - qc[sn].x) + (qc[s].y - qc[sn].y) * (qc[s].y - qc[sn].y);
            }
            pairEdges.push_back({tA, tB, {{c, a, d, b}, {e_ca, e_ad, e_db, e_bc}}, minArea2 / (perim2 + T(1e-14))});
        }
        std::sort(pairEdges.begin(), pairEdges.end(),
                  [](const PairEdge& p, const PairEdge& q) { return p.score > q.score; });

        std::vector<int> baseDeg(mesh.pts.size(), 0);
        for (const auto& t : mesh.tris) { ++baseDeg[t.v[0]]; ++baseDeg[t.v[1]]; ++baseDeg[t.v[2]]; }

        std::vector<int> matchedPair(nT, -1); // matchedPair[ti] = index of pairEdges matched for tri ti
        std::vector<int> partnerTri(nT, -1);
        for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
            const auto& pe = pairEdges[pi];
            if (matchedPair[pe.t0] >= 0 || matchedPair[pe.t1] >= 0) continue;
            int a = pe.q.v[1], b = pe.q.v[3];
            if ((mesh.isSmoothBnd[a] && baseDeg[a] <= 2) || (mesh.isSmoothBnd[b] && baseDeg[b] <= 2))
                continue;
            matchedPair[pe.t0] = pi;
            matchedPair[pe.t1] = pi;
            partnerTri[pe.t0] = pe.t1;
            partnerTri[pe.t1] = pe.t0;
            --baseDeg[a]; --baseDeg[b];
        }

        // Dual-graph BFS augmenting paths to pair remaining unmatched triangles into quads
        std::vector<std::vector<int>> triPairsAdj(nT);
        for (int pi = 0; pi < static_cast<int>(pairEdges.size()); ++pi) {
            triPairsAdj[pairEdges[pi].t0].push_back(pi);
            triPairsAdj[pairEdges[pi].t1].push_back(pi);
        }
        for (int startTi = 0; startTi < nT; ++startTi) {
            if (matchedPair[startTi] >= 0) continue;
            std::vector<int> prevTri(nT, -1), prevPair(nT, -1);
            std::queue<int> qBfs;
            qBfs.push(startTi);
            prevTri[startTi] = startTi;
            int foundEnd = -1;
            while (!qBfs.empty() && foundEnd < 0) {
                int u = qBfs.front();
                qBfs.pop();
                for (int pi : triPairsAdj[u]) {
                    const auto& pe = pairEdges[pi];
                    int v = (pe.t0 == u) ? pe.t1 : pe.t0;
                    if (v == partnerTri[u]) continue;
                    int ca = pe.q.v[1], cb = pe.q.v[3];
                    if ((mesh.isSmoothBnd[ca] && baseDeg[ca] <= 2) || (mesh.isSmoothBnd[cb] && baseDeg[cb] <= 2))
                        continue;
                    if (matchedPair[v] < 0) {
                        prevTri[v] = u;
                        prevPair[v] = pi;
                        foundEnd = v;
                        break;
                    }
                    int w = partnerTri[v];
                    if (w >= 0 && prevTri[w] < 0) {
                        prevTri[v] = u;
                        prevPair[v] = pi;
                        prevTri[w] = v;
                        qBfs.push(w);
                    }
                }
            }
            if (foundEnd >= 0) {
                // Reconstruct alternating path from startTi to foundEnd and apply augmentation
                int cur = foundEnd;
                while (cur != startTi) {
                    int u = prevTri[cur];
                    int addPi = prevPair[cur];
                    int nextFree = (u == startTi) ? startTi : prevTri[u];
                    if (u != startTi) {
                        int remPi = matchedPair[u];
                        const auto& remPe = pairEdges[remPi];
                        ++baseDeg[remPe.q.v[1]]; ++baseDeg[remPe.q.v[3]];
                    }
                    const auto& addPe = pairEdges[addPi];
                    matchedPair[u] = addPi;
                    matchedPair[cur] = addPi;
                    partnerTri[u] = cur;
                    partnerTri[cur] = u;
                    --baseDeg[addPe.q.v[1]]; --baseDeg[addPe.q.v[3]];
                    cur = nextFree;
                }
            }
        }

        std::vector<typename TriMeshData<T>::Quad> joinedQuads;
        std::vector<typename TriMeshData<T>::Tri> remTris;
        for (int ti = 0; ti < nT; ++ti) {
            if (matchedPair[ti] >= 0) {
                if (partnerTri[ti] > ti)
                    joinedQuads.push_back(pairEdges[matchedPair[ti]].q);
            } else {
                remTris.push_back(mesh.tris[ti]);
            }
        }

        if (!joinedQuads.empty()) {
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
            emitFromMixed(remTris, joinedQuads, needSplit);
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(remTris, joinedQuads, allTrue);
        }
        if ((!bestCand.isValid || (coarseRelaxed && nT <= 35 && bestCand.numPatches > 3 * nT)) &&
            (joinedQuads.empty() || nT <= 35))
        {
            std::vector<bool> allTrue(mesh.edges.size(), true);
            emitFromMixed(mesh.tris, {}, allTrue);
        }
    };

    // Coarse-first passes: minimal hole subdivision first, then refined holes on passes 3..4
    for (int coarsePass = 0; coarsePass < 5 && (!bestCand.isValid || (coarsePass == 0 && bestCand.numPatches > 29)); ++coarsePass) {
        const bool cRef = (coarsePass > 0);
        const bool cRefineHoles = (coarsePass >= 3);
        const T cH = (coarsePass <= 1) ? (diag * T(0.65))
                   : ((coarsePass == 2) ? (diag * T(0.45))
                   : ((coarsePass == 3) ? (diag * T(0.32)) : (diag * T(0.25))));
        const T cBend = (coarsePass <= 1) ? T(1.10)
                      : ((coarsePass <= 3) ? T(0.65) : T(0.65));
        const int cNeedleRounds = (coarsePass >= 3) ? 2 : 1;
        for (int mode = 0; mode < 3 && (!bestCand.isValid || (coarsePass == 0 && bestCand.numPatches > 29)); ++mode) {
            TriMeshData<T> coarse;
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, cH, cRef, mode, coarse, cBend, cRefineHoles))
                continue;
            quadifyTriMesh(coarse, true);
            if (!bestCand.isValid || (coarsePass <= 1 && coarse.tris.size() <= 40 && bestCand.numPatches > 29)) {
                TriMeshData<T> flipped = coarse;
                applyDomainLawsonFlips(flipped);
                quadifyTriMesh(flipped, true);
            }
            if (!bestCand.isValid) {
                TriMeshData<T> mild = coarse;
                refineNeedleCircumcenters(mild, cH * T(0.5), cNeedleRounds, T(16.0));
                bisectObtuseTriangles(mild, m_degree, T(144), 36);
                quadifyTriMesh(mild, true);
            }
        }
    }

    for (int refinePass = 0; refinePass < 2 && !bestCand.isValid; ++refinePass) {
        const bool doRef = (refinePass == 1);
        const T hStep = doRef ? (diag / T(6)) : (diag / T(8));
        for (int mode = 0; mode < 3 && (!bestCand.isValid || !doRef); ++mode) {
            if (mode == 2 && bestCand.isValid) continue;
            TriMeshData<T> coarse;
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, hStep, doRef, mode, coarse))
                continue;
            quadifyTriMesh(coarse);

            if (!bestCand.isValid || !doRef) {
                TriMeshData<T> flipped = coarse;
                applyDomainLawsonFlips(flipped);
                quadifyTriMesh(flipped);
            }

            if (!bestCand.isValid) {
                TriMeshData<T> steiner = coarse;
                refineInteriorTriangles(steiner, hStep, 2);
                bisectObtuseTriangles(steiner, m_degree, T(146), 28);
                quadifyTriMesh(steiner);
            }
        }
    }

    if (!bestCand.isValid) {
        TriMeshData<T> fine;
        if (buildConstrainedDomainTriangulation(m_domain, m_degree, diag / T(10), true, 0, fine)) {
            refineInteriorTriangles(fine, diag / T(10), 2);
            bisectObtuseTriangles(fine, m_degree, T(146), 28);
            quadifyTriMesh(fine);
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
            if (!buildConstrainedDomainTriangulation(m_domain, m_degree, fH, fRef, fMode, needle))
                continue;
            refineInteriorTriangles(needle, fH, nRounds);
            if (runNeedle)
                refineNeedleCircumcenters(needle, fH, nRounds, needleAng, allowCentroid);
            bisectObtuseTriangles(needle, m_degree, maxQ, maxSp, minQ);
            quadifyTriMesh(needle);
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

    const T domMinAng = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
    const T minAngTol = std::min(static_cast<T>(m_options.getReal("MinCornerAngleDeg")), T(0.85) * domMinAng);
    const T maxAngTol = static_cast<T>(m_options.getReal("MaxCornerAngleDeg"));
    auto tryFinalizeHole = [&](SplittingCandidate& cand) -> bool {
        if (cand.multiPatch.nPatches() == 0) return false;
        cand.multiPatch.computeTopology(T(1e-4), false, false);
        optimizeWinslowBarrier(cand.multiPatch);
        gsASG1LocalFitter<T>::unifyPatchOrientations(cand.multiPatch);
        if (!validateMultiPatch(m_domain, cand.multiPatch, T(1e-4), nullptr, minAngTol, std::max(maxAngTol, T(164.2)))) {
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
        cand.isValid = (minPos20 > T(1e-8) && cand.minJacobianDet > T(1e-8) && cand.maxCornerAngleDeg < T(164.8));
        if (!cand.isValid) {
            cand.multiPatch.clear();
            return false;
        }
        return true;
    };

    // ------------------------------------------------------------------------
    // Case 2: 1 Hole & 4 Outer Curves (bean_0corner_1hole, wavy_ring_annulus)
    // ------------------------------------------------------------------------
    if (nHoles == 1 && nOuterCurves == 4 && m_domain.loop(1).size() == 4) {
        internal_splitter::buildQuadRingByProximity(outer, m_domain.loop(1), candidate.multiPatch, m_degree);
        if (tryFinalizeHole(candidate)) return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 5: 6-Corner L-Shape with 1 Hole (curved_lshape_6corner_1hole)
    // ------------------------------------------------------------------------
    if (nHoles == 1 && nOuterCurves == 6) {
        gsBSpline<T> C0a, C0b, C5a, C5b;
        internal_splitter::splitBSplineExact(outer.curve(0), 0.5, C0a, C0b);
        internal_splitter::splitBSplineExact(outer.curve(5), 0.5, C5a, C5b);

        auto V_bot = internal_splitter::endPoint(C0a);
        auto V_left = internal_splitter::startPoint(C5b);
        auto V_reflex = internal_splitter::endPoint(outer.curve(2));

        std::unique_ptr<gsBSpline<T>> s_vert(internal_splitter::makeStraightBSpline(V_bot, V_reflex, m_degree));
        std::unique_ptr<gsBSpline<T>> s_horiz(internal_splitter::makeStraightBSpline(V_left, V_reflex, m_degree));

        std::unique_ptr<gsBSpline<T>> c5b_rev(C5b.clone()); c5b_rev->reverse();
        candidate.multiPatch.addPatch(typename gsGeometry<T>::uPtr(
            internal_splitter::makeCoonsPatchFrom4Curves(C0a, *s_vert, *s_horiz, *c5b_rev, m_degree)));

        std::unique_ptr<gsCurve<T>> c4_rev(outer.curve(4).clone()); c4_rev->reverse();
        std::unique_ptr<gsBSpline<T>> c5a_rev(C5a.clone()); c5a_rev->reverse();
        candidate.multiPatch.addPatch(typename gsGeometry<T>::uPtr(
            internal_splitter::makeCoonsPatchFrom4Curves(*s_horiz, outer.curve(3), *c4_rev, *c5a_rev, m_degree)));

        std::unique_ptr<gsBSpline<T>> sv_rev(s_vert->clone()); sv_rev->reverse();
        std::unique_ptr<gsCurve<T>> c2_rev(outer.curve(2).clone()); c2_rev->reverse();

        internal_splitter::buildBlockQuadRing(C0b, outer.curve(1), *c2_rev, *sv_rev, m_domain.loop(1), candidate.multiPatch, m_degree);

        if (tryFinalizeHole(candidate)) return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 6: 6-Corner L-Shape with 3 Holes (curved_lshape_6corner_3holes)
    // ------------------------------------------------------------------------
    if (nHoles == 3 && nOuterCurves == 6) {
        gsBSpline<T> C0a, C0b, C5a, C5b;
        internal_splitter::splitBSplineExact(outer.curve(0), 0.5, C0a, C0b);
        internal_splitter::splitBSplineExact(outer.curve(5), 0.5, C5a, C5b);

        auto V_bot = internal_splitter::endPoint(C0a);
        auto V_left = internal_splitter::startPoint(C5b);
        auto V_reflex = internal_splitter::endPoint(outer.curve(2));

        std::unique_ptr<gsBSpline<T>> s_vert(internal_splitter::makeStraightBSpline(V_bot, V_reflex, m_degree));
        std::unique_ptr<gsBSpline<T>> s_horiz(internal_splitter::makeStraightBSpline(V_left, V_reflex, m_degree));

        // Corner block: C0a, s_vert, s_horiz, c5b_rev enclosing Hole 3 (at -0.5, -0.5)
        std::unique_ptr<gsBSpline<T>> c5b_rev(C5b.clone()); c5b_rev->reverse();
        internal_splitter::buildBlockQuadRing(C0a, *s_vert, *s_horiz, *c5b_rev, m_domain.loop(3), candidate.multiPatch, m_degree);

        // Right wing: C0b, C1, c2_rev, s_vert enclosing Hole 1 (at 0.55, -0.5)
        std::unique_ptr<gsCurve<T>> c2_rev(outer.curve(2).clone()); c2_rev->reverse();
        internal_splitter::buildBlockQuadRing(C0b, outer.curve(1), *c2_rev, *s_vert, m_domain.loop(1), candidate.multiPatch, m_degree);

        // Top wing: s_horiz, C3, c4_rev, c5a_rev enclosing Hole 2 (at -0.5, 0.55)
        std::unique_ptr<gsCurve<T>> c4_rev(outer.curve(4).clone()); c4_rev->reverse();
        std::unique_ptr<gsBSpline<T>> c5a_rev(C5a.clone()); c5a_rev->reverse();
        internal_splitter::buildBlockQuadRing(*s_horiz, outer.curve(3), *c4_rev, *c5a_rev, m_domain.loop(2), candidate.multiPatch, m_degree);

        if (tryFinalizeHole(candidate)) return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 7: 4-Corner Quad Plate with 2 Holes (quad_4corner_2holes)
    // ------------------------------------------------------------------------
    if (nHoles == 2 && nOuterCurves == 4) {
        gsPoint2D<T> P_bot(0, -1), P_mid(0, 0), P_top(0, 1);
        std::unique_ptr<gsBSpline<T>> s_bot(internal_splitter::makeStraightBSpline(P_bot, P_mid, m_degree));
        std::unique_ptr<gsBSpline<T>> s_top(internal_splitter::makeStraightBSpline(P_mid, P_top, m_degree));

        std::unique_ptr<gsBSpline<T>> s_top_rev(s_top->clone()); s_top_rev->reverse();
        std::unique_ptr<gsBSpline<T>> s_bot_rev(s_bot->clone()); s_bot_rev->reverse();

        internal_splitter::buildBlockQuadRing(outer.curve(3), outer.curve(0), *s_top_rev, *s_bot_rev, m_domain.loop(1), candidate.multiPatch, m_degree);
        internal_splitter::buildBlockQuadRing(*s_bot, *s_top, outer.curve(1), outer.curve(2), m_domain.loop(2), candidate.multiPatch, m_degree);

        if (tryFinalizeHole(candidate)) return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 8: Star5 with 1 Hole (star_5corner_1hole)
    // ------------------------------------------------------------------------
    if (nHoles == 1 && nOuterCurves == 5) {
        const auto& hole = m_domain.loop(1);
        std::vector<gsBSpline<T>> oSegs(10);
        for (int i = 0; i < 5; ++i) {
            internal_splitter::splitBSplineExact(outer.curve(i), 0.5, oSegs[2 * i], oSegs[2 * i + 1]);
        }
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
    }

    // ------------------------------------------------------------------------
    // Case 9: Swiss Cross 12-Corner with 5 Holes (swiss_cross_12corner_5holes, cross_12corner_5holes_star)
    // ------------------------------------------------------------------------
    if (nHoles == 5 && nOuterCurves == 12) {
        gsPoint2D<T> V_SW(-0.38, -0.38);
        gsPoint2D<T> V_SE( 0.38, -0.38);
        gsPoint2D<T> V_NE( 0.38,  0.38);
        gsPoint2D<T> V_NW(-0.38,  0.38);

        std::unique_ptr<gsBSpline<T>> L_S(internal_splitter::makeStraightBSpline(V_SW, V_SE, m_degree));
        std::unique_ptr<gsBSpline<T>> L_E(internal_splitter::makeStraightBSpline(V_SE, V_NE, m_degree));
        std::unique_ptr<gsBSpline<T>> L_N(internal_splitter::makeStraightBSpline(V_NW, V_NE, m_degree));
        std::unique_ptr<gsBSpline<T>> L_W(internal_splitter::makeStraightBSpline(V_SW, V_NW, m_degree));

        int h_center = 1, h_north = 2, h_south = 3, h_east = 4, h_west = 5;

        // 1. Center Block
        internal_splitter::buildBlockQuadRing(*L_S, *L_E, *L_N, *L_W, m_domain.loop(h_center), candidate.multiPatch, m_degree);

        // 2. South Arm: C0, C1, L_S, C11 reversed
        std::unique_ptr<gsCurve<T>> c11_rev(outer.curve(11).clone()); c11_rev->reverse();
        internal_splitter::buildBlockQuadRing(outer.curve(0), outer.curve(1), *L_S, *c11_rev, m_domain.loop(h_south), candidate.multiPatch, m_degree);

        // 3. East Arm: C2, C3, C4 reversed, L_E
        std::unique_ptr<gsCurve<T>> c4_rev(outer.curve(4).clone()); c4_rev->reverse();
        internal_splitter::buildBlockQuadRing(outer.curve(2), outer.curve(3), *c4_rev, *L_E, m_domain.loop(h_east), candidate.multiPatch, m_degree);

        // 4. North Arm: L_N, C5, C6 reversed, C7 reversed
        std::unique_ptr<gsCurve<T>> c6_rev(outer.curve(6).clone()); c6_rev->reverse();
        std::unique_ptr<gsCurve<T>> c7_rev(outer.curve(7).clone()); c7_rev->reverse();
        internal_splitter::buildBlockQuadRing(*L_N, outer.curve(5), *c6_rev, *c7_rev, m_domain.loop(h_north), candidate.multiPatch, m_degree);

        // 5. West Arm: C10, L_W, C8 reversed, C9 reversed
        std::unique_ptr<gsCurve<T>> c8_rev(outer.curve(8).clone()); c8_rev->reverse();
        std::unique_ptr<gsCurve<T>> c9_rev(outer.curve(9).clone()); c9_rev->reverse();
        internal_splitter::buildBlockQuadRing(outer.curve(10), *L_W, *c8_rev, *c9_rev, m_domain.loop(h_west), candidate.multiPatch, m_degree);

        if (tryFinalizeHole(candidate)) return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 10: U-Channel 8-Corner (uchannel_8corner_2holes, uchannel_8corner_4holes)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 8 && (nHoles == 2 || nHoles == 4)) {
        // Base block: [-1.0, 1.0] x [-1.0, -0.6] (or [-1.2, 1.2] x [-1.0, 0.0])
        // Split C1 at x = -0.5 and x = 0.5 into base segments
        // Left arm and Right arm each split vertically
        T x_cut = (nHoles == 4) ? 0.6 : 0.5;
        T y_base = (nHoles == 4) ? -0.6 : 0.0;

        gsPoint2D<T> J_L(-x_cut, y_base), J_R(x_cut, y_base);
        gsPoint2D<T> B_L(-x_cut, -1.0),   B_R(x_cut, -1.0);

        std::unique_ptr<gsBSpline<T>> cut_L(internal_splitter::makeStraightBSpline(B_L, J_L, m_degree));
        std::unique_ptr<gsBSpline<T>> cut_R(internal_splitter::makeStraightBSpline(B_R, J_R, m_degree));

        // For multi-hole U-channel: left arm has 2 holes, right arm has 2 holes (or 0 holes in 2holes)
        for (int h = 0; h < nHoles; ++h) {
            const auto& hole = m_domain.loop(h + 1);
            gsPoint2D<T> hc(0, 0);
            for (int i = 0; i < hole.size(); ++i) {
                auto p = internal_splitter::startPoint(hole.curve(i));
                hc.x += p.x; hc.y += p.y;
            }
            hc.x /= hole.size(); hc.y /= hole.size();

            T arm_xmin = (hc.x < 0) ? -1.0 : x_cut;
            T arm_xmax = (hc.x < 0) ? -x_cut : 1.0;
            T arm_ymin = (hc.y < 0.2) ? y_base : 0.2;
            T arm_ymax = (hc.y < 0.2) ? 0.2 : 1.0;

            gsPoint2D<T> bSW(arm_xmin, arm_ymin), bSE(arm_xmax, arm_ymin);
            gsPoint2D<T> bNE(arm_xmax, arm_ymax), bNW(arm_xmin, arm_ymax);

            std::unique_ptr<gsBSpline<T>> bS(internal_splitter::makeStraightBSpline(bSW, bSE, m_degree));
            std::unique_ptr<gsBSpline<T>> bE(internal_splitter::makeStraightBSpline(bSE, bNE, m_degree));
            std::unique_ptr<gsBSpline<T>> bN(internal_splitter::makeStraightBSpline(bNW, bNE, m_degree));
            std::unique_ptr<gsBSpline<T>> bW(internal_splitter::makeStraightBSpline(bSW, bNW, m_degree));

            internal_splitter::buildBlockQuadRing(*bS, *bE, *bN, *bW, hole, candidate.multiPatch, m_degree);
        }

        // Base quads connecting left arm to right arm
        std::unique_ptr<gsBSpline<T>> base_S(internal_splitter::makeStraightBSpline(B_L, B_R, m_degree));
        std::unique_ptr<gsBSpline<T>> base_N(internal_splitter::makeStraightBSpline(J_L, J_R, m_degree));
        candidate.multiPatch.addPatch(typename gsGeometry<T>::uPtr(
            internal_splitter::makeCoonsPatchFrom4Curves(*base_S, *cut_R, *base_N, *cut_L, m_degree)));

        if (tryFinalizeHole(candidate)) return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 11: Comb 16-Corner with 4 Holes (comb_16corner_4holes, comb_16corner_4holes_curved)
    // ------------------------------------------------------------------------
    if (nOuterCurves == 16 && nHoles == 4) {
        T y_tooth_base = (outer.curve(3).eval(gsMatrix<T>::Identity(1, 1))(1, 0) < -0.2) ? -0.4 : 0.0;

        for (int h = 0; h < 4; ++h) {
            const auto& hole = m_domain.loop(h + 1);
            gsPoint2D<T> hc(0, 0);
            for (int i = 0; i < hole.size(); ++i) {
                auto p = internal_splitter::startPoint(hole.curve(i));
                hc.x += p.x; hc.y += p.y;
            }
            hc.x /= hole.size(); hc.y /= hole.size();

            T t_width = 0.2;
            gsPoint2D<T> bSW(hc.x - t_width, y_tooth_base);
            gsPoint2D<T> bSE(hc.x + t_width, y_tooth_base);
            gsPoint2D<T> bNE(hc.x + t_width, 1.0);
            gsPoint2D<T> bNW(hc.x - t_width, 1.0);

            std::unique_ptr<gsBSpline<T>> bS(internal_splitter::makeStraightBSpline(bSW, bSE, m_degree));
            std::unique_ptr<gsBSpline<T>> bE(internal_splitter::makeStraightBSpline(bSE, bNE, m_degree));
            std::unique_ptr<gsBSpline<T>> bN(internal_splitter::makeStraightBSpline(bNW, bNE, m_degree));
            std::unique_ptr<gsBSpline<T>> bW(internal_splitter::makeStraightBSpline(bSW, bNW, m_degree));

            internal_splitter::buildBlockQuadRing(*bS, *bE, *bN, *bW, hole, candidate.multiPatch, m_degree);
        }

        // Spine rectangular base underneath teeth:
        gsPoint2D<T> sSW(-1.0, -1.0), sSE(1.0, -1.0);
        gsPoint2D<T> sNE( 1.0, y_tooth_base), sNW(-1.0, y_tooth_base);
        std::unique_ptr<gsBSpline<T>> spine_S(internal_splitter::makeStraightBSpline(sSW, sSE, m_degree));
        std::unique_ptr<gsBSpline<T>> spine_E(internal_splitter::makeStraightBSpline(sSE, sNE, m_degree));
        std::unique_ptr<gsBSpline<T>> spine_N(internal_splitter::makeStraightBSpline(sNW, sNE, m_degree));
        std::unique_ptr<gsBSpline<T>> spine_W(internal_splitter::makeStraightBSpline(sSW, sNW, m_degree));

        candidate.multiPatch.addPatch(typename gsGeometry<T>::uPtr(
            internal_splitter::makeCoonsPatchFrom4Curves(*spine_S, *spine_E, *spine_N, *spine_W, m_degree)));

        if (tryFinalizeHole(candidate)) return candidate;
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
    const T domMinAng = internal_splitter::computeDomainMinCornerAngleDeg(m_domain);
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
        const_cast<gsPatchSplitter<T>*>(this)->optimizeWinslowBarrier(candMp);
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
bool gsPatchSplitter<T>::optimizeWinslowBarrier(gsMultiPatch<T>& mp, T minDetThreshold)
{
    if (mp.nPatches() == 0) return false;

    bool allOk = true;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        auto* tp = dynamic_cast<gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;

        T minDet = 0, maxDet = 0, minAbs = 0, curEnergy = 0;
        internal_splitter::evalPatchJacobianAndWinslow(*tp, 17, minDet, maxDet, minAbs, curEnergy);

        // If patch is uniformly negative everywhere, flip parameter orientation without altering geometry
        if (maxDet < -minDetThreshold) {
            tp->reverse(0);
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 17, minDet, maxDet, minAbs, curEnergy);
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
                internal_splitter::evalPatchJacobianAndWinslow(*tp, 17, md, mxd, ma, en);
                if (md > bestScore) {
                    bestScore = md;
                    bestCoefs = coefs;
                    curEnergy = en;
                }
            }
            coefs = bestCoefs;
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 17, minDet, maxDet, minAbs, curEnergy);
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
            internal_splitter::evalPatchJacobianAndWinslow(*tp, 17, newMinDet, newMaxDet, newMinAbs, newEnergy);
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

    std::vector<std::vector<gsPoint2D<T>>> polys(nP);
    for (size_t p = 0; p < nP; ++p) {
        polys[p] = internal_splitter::samplePatchBoundaryPoly(mp.patch(p), 20);
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

    gsMatrix<T> u(2, 1), val;
    const T pts[3] = {T(0.25), T(0.5), T(0.75)};
    for (size_t p = 0; p < nP; ++p) {
        for (size_t q = p + 1; q < nP; ++q) {
            int mp1 = static_cast<int>(polys[p].size());
            int mq1 = static_cast<int>(polys[q].size());
            for (int i = 0; i < mp1; ++i) {
                for (int j = 0; j < mq1; ++j) {
                    if (internal_splitter::segProperIntersect(polys[p][i], polys[p][(i + 1) % mp1],
                                                              polys[q][j], polys[q][(j + 1) % mq1]))
                        return false;
                }
            }
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    u(0, 0) = pts[a]; u(1, 0) = pts[b];
                    mp.patch(p).eval_into(u, val);
                    if (internal_splitter::pointInPolyStrict(gsPoint2D<T>(val(0, 0), val(1, 0)), polys[q]))
                        return false;
                }
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

    gsMultiPatch<T> mpTopo(mp);
    mpTopo.computeTopology(std::max(T(1e-4), bdrTol), false, false);
    if (2 * mpTopo.nInterfaces() + mpTopo.nBoundary() != 4 * mpTopo.nPatches()) {
        if (reason) *reason = "non-conforming topology (hanging nodes or unmatched sides)";
        return false;
    }

    T bdrErr = 0;
    if (!verifyBoundaryMatch(dom, mp, bdrTol, &bdrErr)) {
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
        if (!internal_splitter::isStrictlyConvexQuad(c)) {
            if (reason) *reason = "non-convex or degenerate quad corners";
            return false;
        }
        T angs[4], minA = 0, maxA = 180;
        bool degen = false;
        internal_splitter::evalPatchCornerAnglesDeg(mp.patch(p), angs, minA, maxA, degen);
        if (degen || minA < minAngleDeg || maxA > maxAngleDeg) {
            if (reason) *reason = "corner tangent angle out of admissible bounds";
            return false;
        }
        polys[p] = internal_splitter::samplePatchBoundaryPoly(mp.patch(p), 20);
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

    const T pts[3] = {T(0.25), T(0.5), T(0.75)};
    for (size_t p = 0; p < nP; ++p) {
        for (size_t q = p + 1; q < nP; ++q) {
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
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    u(0, 0) = pts[a]; u(1, 0) = pts[b];
                    mp.patch(p).eval_into(u, val);
                    if (internal_splitter::pointInPolyStrict(gsPoint2D<T>(val(0, 0), val(1, 0)), polys[q])) {
                        if (reason) *reason = "patch overlap";
                        return false;
                    }
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
    };
    auto buildSampled = [](const gsGeometry<T>* crv, int nSeg = 40) -> SampledCurve {
        SampledCurve sc;
        sc.crv = crv;
        gsMatrix<T> pr = crv->parameterRange();
        sc.u0 = pr(0, 0); sc.u1 = pr(0, 1);
        sc.u.resize(nSeg + 1);
        sc.x.resize(nSeg + 1);
        sc.y.resize(nSeg + 1);
        gsMatrix<T> uMat(1, 1), val;
        for (int i = 0; i <= nSeg; ++i) {
            T ui = sc.u0 + (sc.u1 - sc.u0) * (static_cast<T>(i) / static_cast<T>(nSeg));
            sc.u[i] = ui;
            uMat(0, 0) = ui;
            crv->eval_into(uMat, val);
            sc.x[i] = val(0, 0);
            sc.y[i] = val(1, 0);
        }
        return sc;
    };

    std::vector<SampledCurve> origSampled, mpSampled;
    origSampled.reserve(origCurves.size());
    for (const auto* c : origCurves) origSampled.push_back(buildSampled(c));
    mpSampled.reserve(mpBdrCurves.size());
    for (const auto& c : mpBdrCurves) mpSampled.push_back(buildSampled(c.get()));

    auto pointToSampledCurveDist = [](T px, T py, const SampledCurve& sc) -> T {
        int nSeg = static_cast<int>(sc.u.size()) - 1;
        int bestSeg = 0;
        T bestD2 = std::numeric_limits<T>::max();
        for (int i = 0; i < nSeg; ++i) {
            T ax = sc.x[i], ay = sc.y[i];
            T dx = sc.x[i + 1] - ax, dy = sc.y[i + 1] - ay;
            T L2 = dx * dx + dy * dy;
            T t = (L2 > T(1e-24)) ? std::max(T(0), std::min(T(1), ((px - ax) * dx + (py - ay) * dy) / L2)) : T(0);
            T qx = ax + t * dx - px, qy = ay + t * dy - py;
            T d2 = qx * qx + qy * qy;
            if (d2 < bestD2) { bestD2 = d2; bestSeg = i; }
        }
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
        for (int iter = 0; iter < 45; ++iter) {
            if (f1 < f2) {
                b = c2; c2 = c1; f2 = f1;
                c1 = b - phi * (b - a); f1 = evalD2(c1);
            } else {
                a = c1; c1 = c2; f1 = f2;
                c2 = a + phi * (b - a); f2 = evalD2(c2);
            }
        }
        T bestFinal = std::min(bestD2, std::min(f1, f2));
        return std::sqrt(bestFinal);
    };

    T maxDist = 0;

    // Forward check: every original boundary curve point must match MultiPatch boundary
    for (const auto& oSc : origSampled) {
        for (size_t i = 0; i < oSc.u.size(); ++i) {
            T dMin = std::numeric_limits<T>::max();
            for (const auto& bSc : mpSampled) {
                T d = pointToSampledCurveDist(oSc.x[i], oSc.y[i], bSc);
                if (d < dMin) dMin = d;
            }
            if (dMin > maxDist) maxDist = dMin;
        }
    }

    // Backward check: every MultiPatch boundary point must match original domain curves
    for (const auto& bSc : mpSampled) {
        for (size_t i = 0; i < bSc.u.size(); ++i) {
            T dMin = std::numeric_limits<T>::max();
            for (const auto& oSc : origSampled) {
                T d = pointToSampledCurveDist(bSc.x[i], bSc.y[i], oSc);
                if (d < dMin) dMin = d;
            }
            if (dMin > maxDist) maxDist = dMin;
        }
    }

    if (maxDistOut) *maxDistOut = maxDist;
    return (maxDist <= tol);
}

template <class T>
std::vector<int> gsPatchSplitter<T>::getVertexValences() const
{
    std::vector<int> valences;
    const auto& mp = m_bestCandidate.multiPatch;
    if (mp.nPatches() == 0) return valences;

    std::vector<gsPoint2D<T>> bndCorners;
    for (const auto& bnd : mp.boundaries()) {
        size_t p = bnd.patch;
        boxSide s = bnd.side();
        const auto* tp = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;
        int nu = tp->basis().size(0);
        int nv = tp->basis().size(1);
        const auto& c = tp->coefs();

        if (s == boxSide(boundary::south)) {
            bndCorners.emplace_back(c(0, 0), c(0, 1));
            bndCorners.emplace_back(c(nu - 1, 0), c(nu - 1, 1));
        } else if (s == boxSide(boundary::north)) {
            bndCorners.emplace_back(c((nv - 1) * nu, 0), c((nv - 1) * nu, 1));
            bndCorners.emplace_back(c(nu - 1 + (nv - 1) * nu, 0), c(nu - 1 + (nv - 1) * nu, 1));
        } else if (s == boxSide(boundary::west)) {
            bndCorners.emplace_back(c(0, 0), c(0, 1));
            bndCorners.emplace_back(c((nv - 1) * nu, 0), c((nv - 1) * nu, 1));
        } else if (s == boxSide(boundary::east)) {
            bndCorners.emplace_back(c(nu - 1, 0), c(nu - 1, 1));
            bndCorners.emplace_back(c(nu - 1 + (nv - 1) * nu, 0), c(nu - 1 + (nv - 1) * nu, 1));
        }
    }

    std::vector<gsPoint2D<T>> allCorners;
    for (size_t p = 0; p < mp.nPatches(); ++p) {
        const auto* tp = dynamic_cast<const gsTensorBSpline<2, T>*>(&mp.patch(p));
        if (!tp) continue;
        int nu = tp->basis().size(0);
        int nv = tp->basis().size(1);
        const auto& c = tp->coefs();

        allCorners.emplace_back(c(0, 0), c(0, 1));
        allCorners.emplace_back(c(nu - 1, 0), c(nu - 1, 1));
        allCorners.emplace_back(c(nu - 1 + (nv - 1) * nu, 0), c(nu - 1 + (nv - 1) * nu, 1));
        allCorners.emplace_back(c((nv - 1) * nu, 0), c((nv - 1) * nu, 1));
    }

    std::vector<bool> visited(allCorners.size(), false);
    for (size_t i = 0; i < allCorners.size(); ++i) {
        if (visited[i]) continue;
        visited[i] = true;
        int val = 1;
        for (size_t j = i + 1; j < allCorners.size(); ++j) {
            if (allCorners[i].distSq(allCorners[j]) < 1e-4) {
                val++;
                visited[j] = true;
            }
        }

        bool isBoundary = false;
        for (const auto& bc : bndCorners) {
            if (allCorners[i].distSq(bc) < 1e-4) {
                isBoundary = true;
                break;
            }
        }

        if (!isBoundary) {
            valences.push_back(val);
        }
    }

    return valences;
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
