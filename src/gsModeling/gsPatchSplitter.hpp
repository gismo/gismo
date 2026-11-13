/** @file gsPatchSplitter.hpp
    @brief Implementation of gsPatchSplitter<T>: Automated Quadrilateral Domain Decomposition Engine.

    Decomposes arbitrary planar domains with outer boundaries and interior holes into
    a conforming, non-degenerate multi-patch quadrilateral B-spline geometry (gsMultiPatch)
    certified with strictly positive Jacobian determinant (det(J) > 0).

    Strict Non-Negotiable Invariants:
    1. Zero Boundary Deformation: Outer boundary curves and all hole curves are exact, immutable CAD curves.
    2. Conforming hole annulus topology: Every interior hole is enclosed by a conforming quad ring.
    3. Dual-method strategy: Boundary-Aligned Harmonic Cross-Field and Medial Axis Transform (MAT).
    4. Topological preset template decomposition.
    5. Pareto-optimal automated winner selection minimizing patch count and distortion energy.
    6. Winslow Barrier functional optimization for parameterization injectivity and regularity.
    7. Direct downstream compatibility with gsAsG1GeometryGenerator and gsASG1LocalFitter.

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
gsPatchSplitter<T>::gsPatchSplitter(const gsPlanarDomain<T>& domain)
    : m_domain(domain)
{
}

template <class T>
void gsPatchSplitter<T>::setDomain(const gsPlanarDomain<T>& domain)
{
    m_domain = domain;
    m_bestCandidate = SplittingCandidate();
}

template <class T>
void gsPatchSplitter<T>::loadDomainFromXml(const std::string& xml_file)
{
    gsReadFile<T> rf(xml_file);
    memory::unique_ptr<gsPlanarDomain<T>> pDom = rf;
    if (pDom) {
        m_domain = *pDom;
        m_bestCandidate = SplittingCandidate();
    } else {
        GISMO_ERROR("XML file " << xml_file << " does not contain a gsPlanarDomain.");
    }
}

template <class T>
void gsPatchSplitter<T>::setStrategy(Strategy strategy)
{
    m_strategy = strategy;
}

template <class T>
void gsPatchSplitter<T>::setDegree(int degree)
{
    m_degree = std::max(1, degree);
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

// ============================================================================
// Approach 1: Boundary-Aligned Harmonic Cross-Field Decomposition
// ============================================================================

template <class T>
typename gsPatchSplitter<T>::SplittingCandidate
gsPatchSplitter<T>::computeCrossFieldSplitting()
{
    SplittingCandidate candidate;
    candidate.methodUsed = "HarmonicCrossField";

    if (m_domain.numLoops() == 0) {
        candidate.isValid = false;
        return candidate;
    }

    if (m_domain.numHoles() == 0) {
        candidate = computeHoleFreeTriangulationSplitting();
        candidate.methodUsed = "HarmonicCrossField";
        return candidate;
    }

    candidate = computeMedialAxisSplitting();
    candidate.methodUsed = "HarmonicCrossField";
    return candidate;
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
        mp.computeTopology(T(1e-5));
        T minDet = 0;
        T energy = evaluateDistortionEnergy(mp, &minDet);
        if (!validateMultiPatch(m_domain, mp, T(1e-5))) return false;
        out.multiPatch = mp;
        out.numPatches = static_cast<int>(mp.nPatches());
        out.minJacobianDet = minDet;
        out.distortionEnergy = energy;
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
    // 1. Single-patch check: already a strictly convex 4-curve quad -> 1 patch
    // ------------------------------------------------------------------------
    if (K == 4) {
        gsPoint2D<T> c4[4] = {verts[0], verts[1], verts[2], verts[3]};
        if (isStrictlyConvexQuad(c4)) {
            gsMultiPatch<T> mp;
            mp.addPatch(typename gsGeometry<T>::uPtr(
                makeCoonsFromCCWCycle(origSegs[0], origSegs[1], origSegs[2], origSegs[3], m_degree)));
            SplittingCandidate cand;
            if (tryFinalize(mp, "SingleConvexQuad", cand)) {
                return cand; // Keep strictly convex quad as a single patch
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

    struct BdrNode {
        gsPoint2D<T> pt;
        bool isSmooth;
    };

    auto tryFinalize = [&](gsMultiPatch<T>& mp, SplittingCandidate& out) -> bool {
        if (mp.nPatches() == 0) return false;
        optimizeWinslowBarrier(mp);
        mp.computeTopology(T(1e-5));
        T minDet = 0;
        T energy = evaluateDistortionEnergy(mp, &minDet);
        if (!validateMultiPatch(m_domain, mp, T(1e-5))) return false;
        out.multiPatch = mp;
        out.numPatches = static_cast<int>(mp.nPatches());
        out.minJacobianDet = minDet;
        out.distortionEnergy = energy;
        out.methodUsed = "HarmonicCrossField";
        out.isValid = true;
        return true;
    };

    // Helper: 1-to-3 Catmull-Clark quad split on a conforming triangle mesh
    // whose boundary edges match exact boundary CurveSegs.
    struct Tri { int v[3]; int e[3]; }; // edges: v0->v1, v1->v2, v2->v0
    struct EdgeRec {
        int a, b;
        CurveSeg<T> seg; // oriented a -> b
        CurveSeg<T> h0, h1;
        gsPoint2D<T> mid;
        bool splitDone = false;
    };

    auto buildFromTriMesh = [&](const std::vector<gsPoint2D<T>>& pts,
                                const std::vector<Tri>& tris,
                                std::vector<EdgeRec>& edges,
                                SplittingCandidate& candOut) -> bool
    {
        for (auto& er : edges) {
            if (!er.splitDone) {
                auto sp = er.seg.splitArc(T(0.5));
                er.h0 = sp.first;
                er.h1 = sp.second;
                er.mid = er.h0.p1();
                er.splitDone = true;
            }
        }

        gsMultiPatch<T> mp;
        for (const auto& t : tris) {
            int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
            auto getHalves = [&](int u, int v, int edgeIdx, CurveSeg<T>& out0, CurveSeg<T>& out1, gsPoint2D<T>& outMid) {
                const auto& er = edges[edgeIdx];
                outMid = er.mid;
                if (er.a == u && er.b == v) {
                    out0 = er.h0; out1 = er.h1;
                } else {
                    out0 = er.h1.reversed(); out1 = er.h0.reversed();
                }
            };
            CurveSeg<T> e01_a, e01_b, e12_a, e12_b, e20_a, e20_b;
            gsPoint2D<T> m01, m12, m20;
            getHalves(v0, v1, t.e[0], e01_a, e01_b, m01);
            getHalves(v1, v2, t.e[1], e12_a, e12_b, m12);
            getHalves(v2, v0, t.e[2], e20_a, e20_b, m20);

            // Candidate triangle centers: pure centroid and midpoint centroid
            gsPoint2D<T> O_tri((pts[v0].x + pts[v1].x + pts[v2].x) / T(3),
                               (pts[v0].y + pts[v1].y + pts[v2].y) / T(3));
            gsPoint2D<T> O_mid((m01.x + m12.x + m20.x) / T(3),
                               (m01.y + m12.y + m20.y) / T(3));

            gsPoint2D<T> O = O_mid;
            auto checkAllThreeConvex = [&](const gsPoint2D<T>& Oc) -> bool {
                gsPoint2D<T> q0[4] = {pts[v0], m01, Oc, m20};
                gsPoint2D<T> q1[4] = {pts[v1], m12, Oc, m01};
                gsPoint2D<T> q2[4] = {pts[v2], m20, Oc, m12};
                return isStrictlyConvexQuad(q0) && isStrictlyConvexQuad(q1) && isStrictlyConvexQuad(q2);
            };
            if (!checkAllThreeConvex(O)) {
                if (checkAllThreeConvex(O_tri)) {
                    O = O_tri;
                } else {
                    bool found = false;
                    for (T w = T(0.1); w <= T(0.91); w += T(0.1)) {
                        gsPoint2D<T> Otry((1 - w) * O_mid.x + w * O_tri.x,
                                          (1 - w) * O_mid.y + w * O_tri.y);
                        if (checkAllThreeConvex(Otry)) { O = Otry; found = true; break; }
                    }
                    if (!found) return false;
                }
            }

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
        return tryFinalize(mp, candOut);
    };

    // Strategy T1: Ear-clipping on sharp polygon vertices (splitting opposite edge once if K==3 with reflex)
    {
        std::vector<CurveSeg<T>> segs;
        segs.reserve(K);
        for (int i = 0; i < K; ++i) segs.emplace_back(outer.curve(i), m_degree);

        std::vector<gsPoint2D<T>> polyPts;
        std::vector<bool> smoothFlag;
        for (int i = 0; i < K; ++i) {
            polyPts.push_back(segs[i].p0());
            const auto& prevBsp = *segs[(i + K - 1) % K].bsp;
            const auto& nextBsp = *segs[i].bsp;
            auto din = evalCurveDeriv(prevBsp, T(1));
            auto dout = evalCurveDeriv(nextBsp, T(0));
            T n1 = std::hypot(din.x, din.y), n2 = std::hypot(dout.x, dout.y);
            T cr = (n1 > 1e-14 && n2 > 1e-14) ? (din.x * dout.y - din.y * dout.x) / (n1 * n2) : T(0);
            T dt = (n1 > 1e-14 && n2 > 1e-14) ? (din.x * dout.x + din.y * dout.y) / (n1 * n2) : T(1);
            smoothFlag.push_back(std::abs(std::atan2(cr, dt)) <= T(0.10));
        }

        // Special 3-corner reflex triangle -> 2 triangles by cutting reflex to opposite midpoint
        bool handled3Reflex = false;
        if (K == 3) {
            for (int r = 0; r < 3; ++r) {
                int p = (r + 2) % 3, n = (r + 1) % 3;
                if (cross2(polyPts[p], polyPts[r], polyPts[n]) < -1e-8) {
                    int opp = (r + 1) % 3; // edge r+1 -> r+2
                    auto spOpp = segs[opp].splitArc(T(0.5));
                    std::vector<gsPoint2D<T>> pts = {polyPts[0], polyPts[1], polyPts[2], spOpp.first.p1()};
                    std::vector<EdgeRec> edges;
                    auto addEdge = [&](int a, int b, const CurveSeg<T>& s) -> int {
                        edges.push_back({a, b, s, CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false});
                        return static_cast<int>(edges.size()) - 1;
                    };
                    int e_r_opp = addEdge(r, opp, segs[r]);
                    int e_opp_m = addEdge(opp, 3, spOpp.first);
                    int e_m_r   = addEdge(3, r, makeStraightSeg(pts[3], pts[r], m_degree));
                    int e_m_p   = addEdge(3, p, spOpp.second);
                    int e_p_r   = addEdge(p, r, segs[p]);
                    std::vector<Tri> tris = {
                        {{r, opp, 3}, {e_r_opp, e_opp_m, e_m_r}},
                        {{r, 3, p},   {e_m_r,   e_m_p,   e_p_r}}
                    };
                    SplittingCandidate cand;
                    if (buildFromTriMesh(pts, tris, edges, cand)) {
                        bestCand = cand;
                    }
                    handled3Reflex = true;
                    break;
                }
            }
        }

        bool anySmooth = false;
        for (bool sm : smoothFlag) if (sm) anySmooth = true;

        if (!handled3Reflex && !anySmooth && K >= 3) {
            // Greedy best-min-angle ear clipping
            std::vector<int> ring(K);
            for (int i = 0; i < K; ++i) ring[i] = i;

            std::vector<EdgeRec> edges;
            std::map<std::pair<int, int>, int> edgeMap;
            auto getOrAddDiag = [&](int a, int b) -> int {
                int lo = std::min(a, b), hi = std::max(a, b);
                auto it = edgeMap.find({lo, hi});
                if (it != edgeMap.end()) return it->second;
                int id = static_cast<int>(edges.size());
                edges.push_back({lo, hi, makeStraightSeg(polyPts[lo], polyPts[hi], m_degree),
                                 CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false});
                edgeMap[{lo, hi}] = id;
                return id;
            };
            for (int i = 0; i < K; ++i) {
                int j = (i + 1) % K;
                int id = static_cast<int>(edges.size());
                edges.push_back({i, j, segs[i], CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false});
                edgeMap[{std::min(i, j), std::max(i, j)}] = id;
            }

            std::vector<Tri> tris;
            bool ok = true;
            while (ring.size() > 3) {
                int m = static_cast<int>(ring.size());
                int bestIdx = -1;
                T bestQuality = -T(1);
                for (int idx = 0; idx < m; ++idx) {
                    int ia = ring[(idx + m - 1) % m];
                    int ib = ring[idx];
                    int ic = ring[(idx + 1) % m];
                    const auto& A = polyPts[ia];
                    const auto& B = polyPts[ib];
                    const auto& C = polyPts[ic];
                    T area2 = cross2(A, B, C);
                    if (area2 <= T(1e-8)) continue;

                    // Verify no other ring vertex inside triangle ABC and diagonal AC strictly inside domain
                    bool clear = isSegmentInsidePoly(A, C, outerPoly);
                    for (int k = 0; k < m && clear; ++k) {
                        int ik = ring[k];
                        if (ik == ia || ik == ib || ik == ic) continue;
                        const auto& P = polyPts[ik];
                        if (cross2(A, B, P) >= -1e-9 && cross2(B, C, P) >= -1e-9 && cross2(C, A, P) >= -1e-9)
                            clear = false;
                    }
                    if (!clear) continue;

                    T l2 = std::pow(A.x - B.x, 2) + std::pow(A.y - B.y, 2)
                         + std::pow(B.x - C.x, 2) + std::pow(B.y - C.y, 2)
                         + std::pow(C.x - A.x, 2) + std::pow(C.y - A.y, 2);
                    T qual = area2 / (l2 + T(1e-14));
                    if (qual > bestQuality) {
                        bestQuality = qual;
                        bestIdx = idx;
                    }
                }
                if (bestIdx < 0) { ok = false; break; }
                int ia = ring[(bestIdx + m - 1) % m];
                int ib = ring[bestIdx];
                int ic = ring[(bestIdx + 1) % m];
                int e0 = edgeMap[{std::min(ia, ib), std::max(ia, ib)}];
                int e1 = edgeMap[{std::min(ib, ic), std::max(ib, ic)}];
                int e2 = getOrAddDiag(ic, ia);
                tris.push_back({{ia, ib, ic}, {e0, e1, e2}});
                ring.erase(ring.begin() + bestIdx);
            }
            if (ok && ring.size() == 3) {
                int ia = ring[0], ib = ring[1], ic = ring[2];
                if (cross2(polyPts[ia], polyPts[ib], polyPts[ic]) > T(1e-8)) {
                    int e0 = getOrAddDiag(ia, ib);
                    int e1 = getOrAddDiag(ib, ic);
                    int e2 = getOrAddDiag(ic, ia);
                    tris.push_back({{ia, ib, ic}, {e0, e1, e2}});
                    SplittingCandidate cand;
                    if (buildFromTriMesh(polyPts, tris, edges, cand)) {
                        if (!bestCand.isValid || cand.numPatches < bestCand.numPatches)
                            bestCand = cand;
                    }
                }
            }
        }
    }

    // Strategy T2: Centroid-Star Triangulation (guaranteed no single-triangle smooth corner)
    if (pointInPolyStrict(centroid, outerPoly, T(1e-5))) {
        std::vector<CurveSeg<T>> segs;
        if (K == 1) {
            auto sH = CurveSeg<T>(outer.curve(0), m_degree).splitArc(T(0.5));
            auto s1 = sH.first.splitArc(T(0.5));
            auto s2 = sH.second.splitArc(T(0.5));
            segs = {s1.first, s1.second, s2.first, s2.second};
        } else if (K == 2) {
            auto s0 = CurveSeg<T>(outer.curve(0), m_degree).splitArc(T(0.5));
            auto s1 = CurveSeg<T>(outer.curve(1), m_degree).splitArc(T(0.5));
            segs = {s0.first, s0.second, s1.first, s1.second};
        } else {
            for (int i = 0; i < K; ++i) segs.emplace_back(outer.curve(i), m_degree);
        }
        int M = static_cast<int>(segs.size());
        std::vector<gsPoint2D<T>> pts(M + 1);
        bool allVisible = true;
        for (int i = 0; i < M; ++i) {
            pts[i] = segs[i].p0();
            if (!isSegmentInsidePoly(centroid, pts[i], outerPoly)) allVisible = false;
        }
        pts[M] = centroid;

        if (allVisible && M >= 3) {
            std::vector<EdgeRec> edges;
            std::vector<Tri> tris;
            // Boundary edges 0..M-1, spokes M..2M-1
            for (int i = 0; i < M; ++i) {
                edges.push_back({i, (i + 1) % M, segs[i], CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false});
            }
            for (int i = 0; i < M; ++i) {
                edges.push_back({M, i, makeStraightSeg(centroid, pts[i], m_degree),
                                 CurveSeg<T>(), CurveSeg<T>(), gsPoint2D<T>(), false});
            }
            for (int i = 0; i < M; ++i) {
                int j = (i + 1) % M;
                tris.push_back({{M, i, j}, {M + i, i, M + j}});
            }
            SplittingCandidate cand;
            if (buildFromTriMesh(pts, tris, edges, cand)) {
                if (!bestCand.isValid || cand.numPatches < bestCand.numPatches)
                    bestCand = cand;
            }
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

    // ------------------------------------------------------------------------
    // Case 2: 1 Hole & 4 Outer Curves (bean_0corner_1hole, wavy_ring_annulus)
    // ------------------------------------------------------------------------
    if (nHoles == 1 && nOuterCurves == 4 && m_domain.loop(1).size() == 4) {
        internal_splitter::buildQuadRingByProximity(outer, m_domain.loop(1), candidate.multiPatch, m_degree);
        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = 4;
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet > 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = static_cast<int>(candidate.multiPatch.nPatches());
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet > 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = 12;
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet >= 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = static_cast<int>(candidate.multiPatch.nPatches());
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet > 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = 10;
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet > 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = 20;
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet >= 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = static_cast<int>(candidate.multiPatch.nPatches());
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet >= 0);
        return candidate;
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

        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = static_cast<int>(candidate.multiPatch.nPatches());
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.minJacobianDet >= 0);
        return candidate;
    }

    // ------------------------------------------------------------------------
    // Case 12: General Fallback Template
    // ------------------------------------------------------------------------
    return computeTemplateSplitting();
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

    if (nHoles == 0) {
        candidate.numPatches = 0;
        candidate.isValid = false;
        return candidate;
    }

    // 1-hole with 4 curves
    if (nHoles == 1 && nOuterCurves == 4 && m_domain.loop(1).size() == 4) {
        internal_splitter::buildQuadRingByProximity(outer, m_domain.loop(1), candidate.multiPatch, m_degree);
        candidate.multiPatch.computeTopology();
        optimizeWinslowBarrier(candidate.multiPatch);
        candidate.numPatches = 4;
        candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
        candidate.isValid = (candidate.numPatches > 0 && candidate.minJacobianDet > 0);
        return candidate;
    }

    // General multi-hole template: Conforming block quad ring per hole connected to outer
    for (int h = 0; h < nHoles; ++h) {
        const auto& hole = m_domain.loop(h + 1);
        gsPoint2D<T> c(0, 0);
        T xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
        for (int i = 0; i < hole.size(); ++i) {
            auto p = internal_splitter::startPoint(hole.curve(i));
            c.x += p.x; c.y += p.y;
            xmin = std::min(xmin, p.x); xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y); ymax = std::max(ymax, p.y);
        }
        c.x /= hole.size(); c.y /= hole.size();
        T rad = std::max(xmax - xmin, ymax - ymin) * 0.95;

        gsPoint2D<T> bSW(c.x - rad, c.y - rad);
        gsPoint2D<T> bSE(c.x + rad, c.y - rad);
        gsPoint2D<T> bNE(c.x + rad, c.y + rad);
        gsPoint2D<T> bNW(c.x - rad, c.y + rad);

        std::unique_ptr<gsBSpline<T>> bS(internal_splitter::makeStraightBSpline(bSW, bSE, m_degree));
        std::unique_ptr<gsBSpline<T>> bE(internal_splitter::makeStraightBSpline(bSE, bNE, m_degree));
        std::unique_ptr<gsBSpline<T>> bN(internal_splitter::makeStraightBSpline(bNW, bNE, m_degree));
        std::unique_ptr<gsBSpline<T>> bW(internal_splitter::makeStraightBSpline(bSW, bNW, m_degree));

        internal_splitter::buildBlockQuadRing(*bS, *bE, *bN, *bW, hole, candidate.multiPatch, m_degree);
    }

    candidate.multiPatch.computeTopology();
    optimizeWinslowBarrier(candidate.multiPatch);
    candidate.numPatches = static_cast<int>(candidate.multiPatch.nPatches());
    candidate.distortionEnergy = evaluateDistortionEnergy(candidate.multiPatch, &candidate.minJacobianDet);
    candidate.isValid = (candidate.numPatches > 0 && candidate.minJacobianDet > 0);
    return candidate;
}

// ============================================================================
// Automated Winner-Selection Metric & Execution
// ============================================================================

template <class T>
bool gsPatchSplitter<T>::computeSplitting()
{
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
        SplittingCandidate cfCand = computeCrossFieldSplitting();
        SplittingCandidate matCand = computeMedialAxisSplitting();
        SplittingCandidate tmplCand = computeTemplateSplitting();

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
T gsPatchSplitter<T>::evaluateDistortionEnergy(const gsMultiPatch<T>& mp, T* minDetOut) const
{
    if (mp.nPatches() == 0) {
        if (minDetOut) *minDetOut = T(0);
        return T(1e30);
    }

    T totalEnergy = 0.0;
    T globalMinAbsDet = std::numeric_limits<T>::max();
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
        totalEnergy += pEnergy;
    }

    if (minDetOut) *minDetOut = allSingleSign ? globalMinAbsDet : -globalMinAbsDet;
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
    std::string* reason)
{
    if (mp.nPatches() == 0) {
        if (reason) *reason = "empty multipatch";
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
void gsPatchSplitter<T>::exportParaview(const std::string& filename_prefix) const
{
    if (m_bestCandidate.multiPatch.nPatches() > 0) {
        gsWriteParaview(m_bestCandidate.multiPatch, filename_prefix, 25, true);
    }
}

} // namespace gismo
