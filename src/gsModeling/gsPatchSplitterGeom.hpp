/** @file gsPatchSplitterGeom.hpp
    @brief 2D geometric primitives, CAD B-spline curve operations, polygon utilities,
           Winslow/Jacobian patch evaluators, and 1D interval slicing for gsPatchSplitter.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsCoonsPatch.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsKnotVector.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <limits>
#include <utility>
#include <iostream>

namespace gismo {
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
inline T arcLengthSplitParam(const gsCurve<T>& crv, T frac = 0.5, int nSamples = 40)
{
    gsMatrix<T> pr = crv.parameterRange();
    T u0 = pr(0, 0), u1 = pr(0, 1);
    std::vector<T> u(nSamples + 1), s(nSamples + 1, T(0));
    gsMatrix<T> uMat(1, nSamples + 1), valMat;
    for (int i = 0; i <= nSamples; ++i) {
        T ui = u0 + (u1 - u0) * (static_cast<T>(i) / static_cast<T>(nSamples));
        u[i] = ui;
        uMat(0, i) = ui;
    }
    crv.eval_into(uMat, valMat);
    for (int i = 1; i <= nSamples; ++i) {
        T dx = valMat(0, i) - valMat(0, i - 1);
        T dy = valMat(1, i) - valMat(1, i - 1);
        s[i] = s[i - 1] + std::hypot(dx, dy);
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

    gsPoint2D<T> p0() const {
        if (!bsp) return gsPoint2D<T>(0, 0);
        const auto& c = bsp->coefs();
        return gsPoint2D<T>(c(0, 0), c(0, 1));
    }
    gsPoint2D<T> p1() const {
        if (!bsp) return gsPoint2D<T>(0, 0);
        const auto& c = bsp->coefs();
        const auto n = c.rows() - 1;
        return gsPoint2D<T>(c(n, 0), c(n, 1));
    }

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
        if (!bsp) return T(0.5);
        const int nS = 40;
        gsMatrix<T> uMat(1, nS + 1), ptsMat;
        for (int i = 0; i <= nS; ++i) {
            uMat(0, i) = static_cast<T>(i) / static_cast<T>(nS);
        }
        bsp->eval_into(uMat, ptsMat);

        T bestU = T(0.5);
        T bestD2 = std::numeric_limits<T>::max();
        for (int i = 0; i <= nS; ++i) {
            T dx = ptsMat(0, i) - q.x, dy = ptsMat(1, i) - q.y;
            T d2 = dx * dx + dy * dy;
            if (d2 < bestD2) { bestD2 = d2; bestU = uMat(0, i); }
        }
        T du = T(1) / static_cast<T>(nS);
        T a = std::max(T(0), bestU - du);
        T b = std::min(T(1), bestU + du);
        const T phi = T(0.6180339887498949);
        gsMatrix<T> uPt(1, 1), ptVal;
        auto evalD2 = [&](T u) -> T {
            uPt(0, 0) = u;
            bsp->eval_into(uPt, ptVal);
            T dx = ptVal(0, 0) - q.x, dy = ptVal(1, 0) - q.y;
            return dx * dx + dy * dy;
        };
        T c1 = b - phi * (b - a), c2 = a + phi * (b - a);
        T f1 = evalD2(c1), f2 = evalD2(c2);
        for (int it = 0; it < 18; ++it) {
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
    int n = degree + 1;
    gsKnotVector<T> kv(0.0, 1.0, 0, n, 1, degree);
    gsMatrix<T> coefs(n, 2);
    for (int i = 0; i < n; ++i) {
        T t = static_cast<T>(i) / static_cast<T>(degree);
        coefs(i, 0) = (1.0 - t) * a.x + t * b.x;
        coefs(i, 1) = (1.0 - t) * a.y + t * b.y;
    }
    CurveSeg<T> seg;
    seg.bsp = std::make_shared<gsBSpline<T>>(kv, coefs);
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
        // Compute signed area of loop polygon to determine orientation
        T signedArea2 = T(0);
        {
            std::vector<gsPoint2D<T>> poly;
            poly.reserve(K * 8);
            for (int i = 0; i < K; ++i) {
                const auto& c = lp.curve(i);
                gsMatrix<T> pr = c.parameterRange();
                for (int s = 0; s < 8; ++s) {
                    T t = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * T(s) / T(8);
                    poly.push_back(evalCurvePoint(c, t));
                }
            }
            for (size_t i = 0; i < poly.size(); ++i) {
                size_t j = (i + 1) % poly.size();
                signedArea2 += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
            }
        }
        // Outer loop CCW (area > 0) -> domain on left; Hole loop CW (area < 0) -> domain on left
        const bool domainOnLeft = (l == 0) ? (signedArea2 >= T(0)) : (signedArea2 <= T(0));
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
            T intAng = domainOnLeft ? (T(180) - turnRad * radToDeg) : (T(180) + turnRad * radToDeg);
            if (intAng > T(1e-3) && intAng < minAng) minAng = intAng;
        }
    }
    return minAng;
}

template <class T>
inline std::vector<gsPoint2D<T>> generateTri1to3CenterCandidates(
    const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
    const gsPoint2D<T>& mAB, const gsPoint2D<T>& mBC, const gsPoint2D<T>& mCA,
    gsPoint2D<T>* outO_inc = nullptr, T* outLsum = nullptr)
{
    gsPoint2D<T> O_mid((mAB.x + mBC.x + mCA.x) / T(3), (mAB.y + mBC.y + mCA.y) / T(3));
    gsPoint2D<T> O_tri((A.x + B.x + C.x) / T(3), (A.y + B.y + C.y) / T(3));
    T la = B.dist(C), lb = C.dist(A), lc = A.dist(B), lsum = la + lb + lc;
    gsPoint2D<T> O_inc = (lsum > T(1e-14))
        ? gsPoint2D<T>((la * A.x + lb * B.x + lc * C.x) / lsum,
                       (la * A.y + lb * B.y + lc * C.y) / lsum)
        : O_tri;
    if (outO_inc) *outO_inc = O_inc;
    if (outLsum) *outLsum = lsum;

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
    return cands;
}

template <class T>
inline bool pickOptimalTri1to3Center(
    const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
    const gsPoint2D<T>& mAB, const gsPoint2D<T>& mBC, const gsPoint2D<T>& mCA,
    gsPoint2D<T>& outO)
{
    std::vector<gsPoint2D<T>> cands = generateTri1to3CenterCandidates(A, B, C, mAB, mBC, mCA);

    auto quadAngleMarginDeg = [](const gsPoint2D<T> q[4], T& minAngOut, T& maxAngOut) -> T {
        const T radToDeg = T(180) / T(3.14159265358979323846);
        T mn = T(360), mx = T(0);
        for (int k = 0; k < 4; ++k) {
            const auto& p0 = q[(k + 3) % 4];
            const auto& pk = q[k];
            const auto& p1 = q[(k + 1) % 4];
            T ux = p1.x - pk.x, uy = p1.y - pk.y;
            T vx = p0.x - pk.x, vy = p0.y - pk.y;
            T nu = std::hypot(ux, uy), nv = std::hypot(vx, vy);
            if (nu <= T(1e-14) || nv <= T(1e-14)) { minAngOut = T(0); maxAngOut = T(360); return -T(1); }
            T cr = (ux * vy - uy * vx) / (nu * nv);
            T dt = (ux * vx + uy * vy) / (nu * nv);
            T ang = std::atan2(cr, dt) * radToDeg;
            if (ang <= T(0)) ang += T(360);
            mn = std::min(mn, ang);
            mx = std::max(mx, ang);
        }
        minAngOut = mn;
        maxAngOut = mx;
        return std::min(mn, T(180) - mx);
    };

    bool found = false;
    T bestScore = -T(1e9);
    for (const auto& Oc : cands) {
        gsPoint2D<T> q0[4] = {A, mAB, Oc, mCA};
        gsPoint2D<T> q1[4] = {B, mBC, Oc, mAB};
        gsPoint2D<T> q2[4] = {C, mCA, Oc, mBC};
        if (!isStrictlyConvexQuad(q0) || !isStrictlyConvexQuad(q1) || !isStrictlyConvexQuad(q2))
            continue;
        T mn0, mx0, mn1, mx1, mn2, mx2;
        T e0 = quadAngleMarginDeg(q0, mn0, mx0);
        T e1 = quadAngleMarginDeg(q1, mn1, mx1);
        T e2 = quadAngleMarginDeg(q2, mn2, mx2);
        T minMargin = std::min(e0, std::min(e1, e2));
        T maxA = std::max(mx0, std::max(mx1, mx2));
        T minA = std::min(mn0, std::min(mn1, mn2));
        T sc = minMargin + T(0.15) * minA - T(0.15) * maxA;
        if (!found || sc > bestScore) {
            bestScore = sc;
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
    const int total = 4 * nPerSide;
    std::vector<gsPoint2D<T>> poly(total);
    gsMatrix<T> u(2, total), val;
    int col = 0;
    for (int i = 0; i < nPerSide; ++i) { u(0, col) = static_cast<T>(i) / nPerSide; u(1, col) = T(0); ++col; }
    for (int i = 0; i < nPerSide; ++i) { u(0, col) = T(1); u(1, col) = static_cast<T>(i) / nPerSide; ++col; }
    for (int i = 0; i < nPerSide; ++i) { u(0, col) = T(1) - static_cast<T>(i) / nPerSide; u(1, col) = T(1); ++col; }
    for (int i = 0; i < nPerSide; ++i) { u(0, col) = T(0); u(1, col) = T(1) - static_cast<T>(i) / nPerSide; ++col; }
    patch.eval_into(u, val);
    for (int i = 0; i < total; ++i) {
        poly[i] = gsPoint2D<T>(val(0, i), val(1, i));
    }
    return poly;
}

template <class T>
inline std::vector<gsPoint2D<T>> sampleCurveLoopPoly(const gsCurveLoop<T>& loop, int nPerCurve = 24)
{
    std::vector<gsPoint2D<T>> poly;
    poly.reserve(loop.size() * nPerCurve);
    gsMatrix<T> uMat(1, nPerCurve), valMat;
    for (int c = 0; c < loop.size(); ++c) {
        const auto& crv = loop.curve(c);
        gsMatrix<T> pr = crv.parameterRange();
        T u0 = pr(0, 0), u1 = pr(0, 1);
        for (int i = 0; i < nPerCurve; ++i) {
            uMat(0, i) = u0 + (u1 - u0) * (static_cast<T>(i) / static_cast<T>(nPerCurve));
        }
        crv.eval_into(uMat, valMat);
        for (int i = 0; i < nPerCurve; ++i) {
            poly.push_back(gsPoint2D<T>(valMat(0, i), valMat(1, i)));
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
inline const gsMatrix<T>& getSamplingGrid2D(int nSamples)
{
    static thread_local std::map<int, gsMatrix<T>> cache;
    auto it = cache.find(nSamples);
    if (it != cache.end()) return it->second;

    const T h = T(1) / (nSamples - 1);
    const int N = nSamples * nSamples;
    gsMatrix<T>& grid = cache[nSamples];
    grid.resize(2, N);
    for (int j = 0; j < nSamples; ++j) {
        T eta = j * h;
        for (int i = 0; i < nSamples; ++i) {
            int p = j * nSamples + i;
            grid(0, p) = i * h;
            grid(1, p) = eta;
        }
    }
    return grid;
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
    const int N = nSamples * nSamples;
    const gsMatrix<T>& stdGrid = getSamplingGrid2D<T>(nSamples);

    const gsMatrix<T> supp = patch.support();
    const T u0 = supp(0, 0), u1 = supp(0, 1);
    const T v0 = supp(1, 0), v1 = supp(1, 1);
    const T du = u1 - u0, dv = v1 - v0;

    gsMatrix<T> scaledGrid;
    const gsMatrix<T>* uGridPtr = &stdGrid;
    const bool isUnit = (std::abs(u0) < T(1e-9) && std::abs(u1 - T(1)) < T(1e-9) &&
                         std::abs(v0) < T(1e-9) && std::abs(v1 - T(1)) < T(1e-9));
    if (!isUnit) {
        scaledGrid.resize(2, N);
        for (int p = 0; p < N; ++p) {
            scaledGrid(0, p) = u0 + stdGrid(0, p) * du;
            scaledGrid(1, p) = v0 + stdGrid(1, p) * dv;
        }
        uGridPtr = &scaledGrid;
    }

    gsMatrix<T> J;
    patch.jacobian_into(*uGridPtr, J);

    const T hSq = h * h;
    for (int p = 0; p < N; ++p) {
        int col0 = 2 * p;
        int col1 = col0 + 1;

        T Jxi_x  = J(0, col0), Jxi_y  = J(1, col0);
        T Jeta_x = J(0, col1), Jeta_y = J(1, col1);

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
            energy += (win + orth * orth) * hSq;
        } else {
            energy += T(1e6);
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

// Corner Tangent Angle & Aspect Ratio Diagnostics
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
            if (ang > T(175.0) || ang < T(5.0))
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

// Concatenate consecutive oriented CurveSeg<T> segments into a single exact clamped gsBSpline<T>.
template <class T>
inline std::shared_ptr<gsBSpline<T>> joinCurveSegsExact(
    const std::vector<CurveSeg<T>>& segs, int degree)
{
    if (segs.empty())
        return nullptr;
    if (segs.size() == 1) {
        auto copy = std::make_shared<gsBSpline<T>>(*segs[0].bsp);
        copy->knots().affineTransformTo(T(0), T(1));
        if (copy->degree() < degree)
            copy->degreeElevate(degree - copy->degree());
        return copy;
    }

    int d = degree;
    for (const auto& s : segs)
        d = std::max(d, s.bsp->degree());

    const int K = static_cast<int>(segs.size());
    std::vector<std::unique_ptr<gsBSpline<T>>> pieces(K);
    std::vector<T> lengths(K, T(0));
    T totalLen = T(0);
    for (int i = 0; i < K; ++i) {
        pieces[i].reset(new gsBSpline<T>(*segs[i].bsp));
        pieces[i]->knots().affineTransformTo(T(0), T(1));
        if (pieces[i]->degree() < d)
            pieces[i]->degreeElevate(d - pieces[i]->degree());
        T L = T(0);
        gsPoint2D<T> prevPt = evalCurvePoint(*pieces[i], T(0));
        for (int s = 1; s <= 24; ++s) {
            gsPoint2D<T> curPt = evalCurvePoint(*pieces[i], static_cast<T>(s) / T(24));
            L += std::hypot(curPt.x - prevPt.x, curPt.y - prevPt.y);
            prevPt = curPt;
        }
        lengths[i] = L;
        totalLen += L;
    }
    const T minFrac = T(1e-3);
    if (totalLen <= T(1e-14)) {
        for (int i = 0; i < K; ++i) lengths[i] = T(1);
        totalLen = static_cast<T>(K);
    } else {
        T adjTotal = T(0);
        for (int i = 0; i < K; ++i) {
            lengths[i] = std::max(lengths[i], minFrac * totalLen);
            adjTotal += lengths[i];
        }
        totalLen = adjTotal;
    }

    std::vector<T> bp(K + 1, T(0));
    for (int i = 0; i < K; ++i)
        bp[i + 1] = bp[i] + lengths[i] / totalLen;
    bp.back() = T(1);

    std::vector<T> allKnots;
    allKnots.reserve(K * (d + 4) + d + 2);
    for (int m = 0; m < d + 1; ++m)
        allKnots.push_back(T(0));

    int totalCoefs = 1;
    for (int i = 0; i < K; ++i) {
        const auto& kv = pieces[i]->knots();
        const T s0 = bp[i], s1 = bp[i + 1];
        const index_t kvSz = static_cast<index_t>(kv.size());
        for (index_t j = d + 1; j < kvSz - (d + 1); ++j) {
            T u = s0 + (s1 - s0) * kv[j];
            if (u > allKnots.back() + T(1e-14) && u < s1 - T(1e-14))
                allKnots.push_back(u);
        }
        if (i + 1 < K) {
            for (int m = 0; m < d; ++m)
                allKnots.push_back(s1);
        }
        totalCoefs += static_cast<int>(pieces[i]->coefs().rows()) - 1;
    }
    for (int m = 0; m < d + 1; ++m)
        allKnots.push_back(T(1));

    // Match knot vector length precisely if near-duplicate interior knots were skipped
    const int expectedKnots = totalCoefs + d + 1;
    if (static_cast<int>(allKnots.size()) != expectedKnots) {
        // Fallback: resample pieces identically via raw clamped concatenation without skipping
        allKnots.clear();
        for (int m = 0; m < d + 1; ++m) allKnots.push_back(T(0));
        for (int i = 0; i < K; ++i) {
            const auto& kv = pieces[i]->knots();
            const T s0 = bp[i], s1 = bp[i + 1];
            const index_t kvSz = static_cast<index_t>(kv.size());
            for (index_t j = d + 1; j < kvSz - (d + 1); ++j) {
                T u = std::max(s0 + T(1e-12), std::min(s1 - T(1e-12), s0 + (s1 - s0) * kv[j]));
                allKnots.push_back(u);
            }
            if (i + 1 < K) {
                for (int m = 0; m < d; ++m) allKnots.push_back(s1);
            }
        }
        for (int m = 0; m < d + 1; ++m) allKnots.push_back(T(1));
    }

    gsMatrix<T> coefs(totalCoefs, 2);
    coefs.row(0) = pieces[0]->coefs().row(0);
    int rowIdx = 1;
    for (int i = 0; i < K; ++i) {
        const gsMatrix<T>& pc = pieces[i]->coefs();
        const int nr = static_cast<int>(pc.rows());
        for (int r = 1; r < nr - 1; ++r)
            coefs.row(rowIdx++) = pc.row(r);
        if (i + 1 < K) {
            const gsMatrix<T>& nc = pieces[i + 1]->coefs();
            coefs.row(rowIdx++) = T(0.5) * (pc.row(nr - 1) + nc.row(0));
        } else {
            coefs.row(rowIdx++) = pc.row(nr - 1);
        }
    }

    gsKnotVector<T> combinedKv(static_cast<short_t>(d), allKnots.begin(), allKnots.end());
    return std::make_shared<gsBSpline<T>>(combinedKv, coefs);
}

// Convert a hole gsCurveLoop<T> into an ordered CCW sequence of CurveSeg<T>.
template <class T>
inline std::vector<CurveSeg<T>> orientLoopSegments(const gsCurveLoop<T>& loop, int degree, bool wantCCW = true)
{
    const int K = loop.size();
    std::vector<CurveSeg<T>> segs;
    segs.reserve(K);
    for (int i = 0; i < K; ++i)
        segs.emplace_back(loop.curve(i), degree);

    T area2 = T(0);
    std::vector<gsPoint2D<T>> poly;
    poly.reserve(K * 8);
    for (int i = 0; i < K; ++i) {
        for (int s = 0; s < 8; ++s) {
            T u = static_cast<T>(s) / T(8);
            poly.push_back(evalCurvePoint(*segs[i].bsp, u));
        }
    }
    for (size_t i = 0; i < poly.size(); ++i) {
        size_t j = (i + 1) % poly.size();
        area2 += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    if ((wantCCW && area2 < T(0)) || (!wantCCW && area2 > T(0))) {
        std::reverse(segs.begin(), segs.end());
        for (auto& s : segs) s = s.reversed();
    }
    return segs;
}

template <class T>
inline bool isHoleSpokeClearOfHole(
    const gsPoint2D<T>& holePt, const gsPoint2D<T>& outerPt,
    const std::vector<gsPoint2D<T>>& holePoly)
{
    const int N = static_cast<int>(holePoly.size());
    if (N < 3) return true;
    // Check midpoint + sample points are strictly outside hole interior
    for (int s = 1; s <= 7; ++s) {
        T t = static_cast<T>(s) / T(8);
        gsPoint2D<T> q((T(1) - t) * holePt.x + t * outerPt.x,
                       (T(1) - t) * holePt.y + t * outerPt.y);
        if (pointInPolyStrict(q, holePoly, T(1e-7)))
            return false;
    }
    for (int i = 0; i < N; ++i) {
        if (segProperIntersect(holePt, outerPt, holePoly[i], holePoly[(i + 1) % N]))
            return false;
    }
    return true;
}

// Compute sharp corner turn angles (degrees) at start of each CCW hole segment.
template <class T>
inline std::vector<T> computeLoopVertexTurnDegs(const std::vector<CurveSeg<T>>& segs)
{
    const int K = static_cast<int>(segs.size());
    const T radToDeg = T(180) / T(3.14159265358979323846);
    std::vector<T> turns(K, T(0));
    for (int i = 0; i < K; ++i) {
        const auto& prevBsp = *segs[(i + K - 1) % K].bsp;
        const auto& nextBsp = *segs[i].bsp;
        auto d0 = evalCurveDeriv(prevBsp, T(1));
        auto d1 = evalCurveDeriv(nextBsp, T(0));
        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
        if (n0 > T(1e-14) && n1 > T(1e-14)) {
            T cr = (d0.x * d1.y - d0.y * d1.x) / (n0 * n1);
            T dt = (d0.x * d1.x + d0.y * d1.y) / (n0 * n1);
            turns[i] = std::abs(std::atan2(cr, dt)) * radToDeg;
        }
    }
    return turns;
}

// Evaluate quality of 4 candidate Coons patches inside a block ring
template <class T>
inline bool evaluateRingPatches(
    const std::vector<std::unique_ptr<gsGeometry<T>>>& patches,
    const std::vector<gsPoint2D<T>>& holePoly,
    T& scoreOut, T& minDetOut, T& maxAngOut,
    std::string* failReason = nullptr)
{
    minDetOut = T(1e9);
    maxAngOut = T(0);
    T minAng = T(180);
    T totalWin = T(0);
    T maxAR = T(1);
    for (size_t pi = 0; pi < patches.size(); ++pi) {
        const auto& p = patches[pi];
        T mnD = 0, mxD = 0, mnA = 0, en = 0;
        evalPatchJacobianAndWinslow(*p, 21, mnD, mxD, mnA, en);
        // Orientation must be single-sign positive (matching validateMultiPatch 21x21 tol=1e-7)
        if (mnD <= T(2e-7)) {
            if (failReason) *failReason = "p=" + std::to_string(pi) + " mnD=" + std::to_string(mnD) + " mxD=" + std::to_string(mxD);
            return false;
        }
        minDetOut = std::min(minDetOut, mnD);
        totalWin += en;

        // Corner angles and chord aspect ratio of patch
        const T radToDeg = T(180) / T(3.14159265358979323846);
        const gsPoint2D<T> corners[4] = {{0,0}, {1,0}, {1,1}, {0,1}};
        gsPoint2D<T> cPos[4];
        gsMatrix<T> pt(2, 1), jac, cVal;
        for (int c = 0; c < 4; ++c) {
            pt(0, 0) = corners[c].x; pt(1, 0) = corners[c].y;
            p->eval_into(pt, cVal);
            cPos[c] = gsPoint2D<T>(cVal(0, 0), cVal(1, 0));
            p->jacobian_into(pt, jac);
            T xu = jac(0, 0), yu = jac(1, 0), xv = jac(0, 1), yv = jac(1, 1);
            // Adjust signs to interior directions along parameter boundaries
            T sx = (corners[c].x > T(0.5)) ? T(-1) : T(1);
            T sy = (corners[c].y > T(0.5)) ? T(-1) : T(1);
            T ax = sx * xu, ay = sx * yu, bx = sy * xv, by = sy * yv;
            T na = std::hypot(ax, ay), nb = std::hypot(bx, by);
            if (na <= T(1e-12) || nb <= T(1e-12)) {
                if (failReason) *failReason = "p=" + std::to_string(pi) + " c=" + std::to_string(c) + " degenEdge";
                return false;
            }
            T cosA = std::max(T(-1), std::min(T(1), (ax * bx + ay * by) / (na * nb)));
            T ang = std::acos(cosA) * radToDeg;
            minAng = std::min(minAng, ang);
            maxAngOut = std::max(maxAngOut, ang);
        }
        T minEdge = T(1e18), maxEdge = T(0);
        for (int e = 0; e < 4; ++e) {
            T el = std::hypot(cPos[(e + 1) % 4].x - cPos[e].x, cPos[(e + 1) % 4].y - cPos[e].y);
            minEdge = std::min(minEdge, el);
            maxEdge = std::max(maxEdge, el);
        }
        if (minEdge <= T(1e-12)) {
            if (failReason) *failReason = "p=" + std::to_string(pi) + " degenChord";
            return false;
        }
        maxAR = std::max(maxAR, maxEdge / minEdge);

        // Check interior points of patch do not intrude into hole
        gsMatrix<T> uv(2, 1), val;
        for (int r = 1; r <= 3; ++r) {
            for (int c = 1; c <= 3; ++c) {
                uv(0, 0) = T(r) / T(4);
                uv(1, 0) = T(c) / T(4);
                p->eval_into(uv, val);
                if (pointInPolyStrict(gsPoint2D<T>(val(0, 0), val(1, 0)), holePoly, T(1e-7))) {
                    if (failReason) *failReason = "p=" + std::to_string(pi) + " interiorInHole";
                    return false;
                }
            }
        }
    }
    if (minAng < T(4.5) || maxAngOut > T(164.5)) {
        if (failReason) *failReason = "ang=[" + std::to_string(minAng) + "," + std::to_string(maxAngOut) + "]";
        return false;
    }
    if (patches.size() > 4 && maxAR > T(30)) {
        if (failReason) *failReason = "maxAR=" + std::to_string(maxAR);
        return false;
    }
    scoreOut = totalWin + T(200) / std::max(minAng, T(1)) + T(200) / std::max(T(180) - maxAngOut, T(1))
             + (minAng < T(12) ? T(1e4) : T(0)) + (minAng < T(8) ? T(1e5) : T(0))
             + (maxAngOut > T(158) ? T(1e4) : T(0)) + (maxAngOut > T(162) ? T(1e5) : T(0))
             + T(5) * std::max(T(0), maxAR - T(12)) + (maxAR > T(20) ? T(5e3) : T(0));
    return true;
}

template <class T>
inline std::vector<CurveSeg<T>> splitSegMulti(const CurveSeg<T>& seg, const std::vector<T>& ts)
{
    std::vector<CurveSeg<T>> out;
    if (ts.empty()) {
        out.push_back(seg);
        return out;
    }
    CurveSeg<T> rem = seg;
    T prevT = T(0);
    for (size_t i = 0; i < ts.size(); ++i) {
        T t = ts[i];
        T localT = (T(1) - prevT > T(1e-14)) ? (t - prevT) / (T(1) - prevT) : T(0.5);
        localT = std::max(T(1e-4), std::min(T(1) - T(1e-4), localT));
        auto sp = rem.splitAtParam(localT);
        out.push_back(sp.first);
        rem = sp.second;
        prevT = t;
    }
    out.push_back(rem);
    return out;
}

template <class T>
struct SeparatorGrid1D {
    static std::vector<T> findClearCuts(const std::vector<std::pair<T, T>>& intervals, T minGap) {
        if (intervals.empty()) return {};
        auto ivals = intervals;
        std::sort(ivals.begin(), ivals.end());
        std::vector<T> cuts;
        T curMax = ivals[0].second;
        for (size_t i = 1; i < ivals.size(); ++i) {
            if (ivals[i].first > curMax + minGap) {
                cuts.push_back(T(0.5) * (curMax + ivals[i].first));
                curMax = ivals[i].second;
            } else {
                curMax = std::max(curMax, ivals[i].second);
            }
        }
        return cuts;
    }

    static int gridMaxItems(const std::vector<T>& xCuts, const std::vector<T>& yCuts,
                            const std::vector<gsPoint2D<T>>& centers) {
        int nx = static_cast<int>(xCuts.size()) + 1;
        int ny = static_cast<int>(yCuts.size()) + 1;
        std::vector<int> counts(nx * ny, 0);
        int mx = 0;
        for (const auto& c : centers) {
            int ix = static_cast<int>(std::upper_bound(xCuts.begin(), xCuts.end(), c.x) - xCuts.begin());
            int iy = static_cast<int>(std::upper_bound(yCuts.begin(), yCuts.end(), c.y) - yCuts.begin());
            mx = std::max(mx, ++counts[iy * nx + ix]);
        }
        return mx;
    }

    static void pruneRedundantCuts(std::vector<T>& xCuts, std::vector<T>& yCuts,
                                   const std::vector<std::pair<T, T>>& xIvals,
                                   const std::vector<std::pair<T, T>>& yIvals,
                                   const std::vector<gsPoint2D<T>>& centers) {
        if (gridMaxItems(xCuts, {}, centers) <= 1 && !xCuts.empty()) {
            yCuts.clear();
        } else if (gridMaxItems({}, yCuts, centers) <= 1 && !yCuts.empty()) {
            xCuts.clear();
        } else if (gridMaxItems(xCuts, yCuts, centers) <= 1 && (xCuts.size() + yCuts.size() > 1)) {
            auto cutMargin = [&](T val, const std::vector<std::pair<T, T>>& ivals) -> T {
                T mn = T(1e18);
                for (const auto& iv : ivals) {
                    mn = std::min(mn, std::min(std::abs(val - iv.first), std::abs(val - iv.second)));
                }
                return mn;
            };
            while (xCuts.size() + yCuts.size() > 1) {
                bool bestIsX = true;
                int bestIdx = -1;
                T smallestM = T(1e18);
                for (size_t i = 0; i < xCuts.size(); ++i) {
                    std::vector<T> xc = xCuts;
                    xc.erase(xc.begin() + static_cast<long>(i));
                    if ((!xc.empty() || !yCuts.empty()) && gridMaxItems(xc, yCuts, centers) <= 1) {
                        T m = cutMargin(xCuts[i], xIvals);
                        if (m < smallestM) { smallestM = m; bestIsX = true; bestIdx = static_cast<int>(i); }
                    }
                }
                for (size_t i = 0; i < yCuts.size(); ++i) {
                    std::vector<T> yc = yCuts;
                    yc.erase(yc.begin() + static_cast<long>(i));
                    if ((!xCuts.empty() || !yc.empty()) && gridMaxItems(xCuts, yc, centers) <= 1) {
                        T m = cutMargin(yCuts[i], yIvals);
                        if (m < smallestM) { smallestM = m; bestIsX = false; bestIdx = static_cast<int>(i); }
                    }
                }
                if (bestIdx < 0) break;
                if (bestIsX) xCuts.erase(xCuts.begin() + bestIdx);
                else         yCuts.erase(yCuts.begin() + bestIdx);
            }
        }
    }
};

// Evaluate exact corner angles (deg) and min(theta, 180-theta) from 4 CCW boundary CurveSegs.
template <class T>
inline void eval4SegCornerAnglesDeg(
    const CurveSeg<T>& c0, const CurveSeg<T>& c1,
    const CurveSeg<T>& c2, const CurveSeg<T>& c3,
    T& minAngOut, T& maxAngOut, T& effMinOut)
{
    const CurveSeg<T>* cs[4] = {&c0, &c1, &c2, &c3};
    const T radToDeg = T(180) / T(3.14159265358979323846);
    minAngOut = T(360);
    maxAngOut = T(0);
    for (int k = 0; k < 4; ++k) {
        const auto& cOut = *cs[k];
        const auto& cIn  = *cs[(k + 3) % 4];
        gsPoint2D<T> d0 = evalCurveDeriv(*cOut.bsp, T(0));
        gsPoint2D<T> d1 = evalCurveDeriv(*cIn.bsp,  T(1));
        T n0 = std::hypot(d0.x, d0.y), n1 = std::hypot(d1.x, d1.y);
        if (n0 <= T(1e-14) || n1 <= T(1e-14)) {
            minAngOut = T(0); maxAngOut = T(360); effMinOut = T(0); return;
        }
        gsPoint2D<T> tOut(d0.x / n0, d0.y / n0);
        gsPoint2D<T> tBack(-d1.x / n1, -d1.y / n1);
        T cr = tOut.x * tBack.y - tOut.y * tBack.x;
        T dt = tOut.x * tBack.x + tOut.y * tBack.y;
        T ang = std::atan2(cr, dt) * radToDeg;
        if (ang <= T(0)) ang += T(360);
        minAngOut = std::min(minAngOut, ang);
        maxAngOut = std::max(maxAngOut, ang);
    }
    effMinOut = std::min(minAngOut, T(180) - maxAngOut);
}

// Curved-tangent-aware 1-to-3 triangle center selection maximizing min(theta, 180-theta)
template <class T>
inline bool pickOptimalTri1to3CenterCurved(
    const gsPoint2D<T>& A, const gsPoint2D<T>& B, const gsPoint2D<T>& C,
    const CurveSeg<T>& e01a, const CurveSeg<T>& e01b,
    const CurveSeg<T>& e12a, const CurveSeg<T>& e12b,
    const CurveSeg<T>& e20a, const CurveSeg<T>& e20b,
    bool isBnd01, bool isBnd12, bool isBnd20,
    gsPoint2D<T>& outO)
{
    const gsPoint2D<T> m01 = e01a.p1(), m12 = e12a.p1(), m20 = e20a.p1();
    if (!isBnd01 && !isBnd12 && !isBnd20)
        return pickOptimalTri1to3Center(A, B, C, m01, m12, m20, outO);

    gsPoint2D<T> O_inc;
    T lsum = 0;
    std::vector<gsPoint2D<T>> cands = generateTri1to3CenterCandidates(A, B, C, m01, m12, m20, &O_inc, &lsum);

    auto unitTan = [&](const CurveSeg<T>& seg, T u, bool flip) -> gsPoint2D<T> {
        gsPoint2D<T> d = evalCurveDeriv(*seg.bsp, u);
        T nm = std::hypot(d.x, d.y);
        if (nm <= T(1e-14)) return gsPoint2D<T>(0, 0);
        T sgn = flip ? T(-1) : T(1);
        return gsPoint2D<T>(sgn * d.x / nm, sgn * d.y / nm);
    };

    const gsPoint2D<T> tA0 = unitTan(e01a, T(0), false), tA1 = unitTan(e20b, T(1), true);
    const gsPoint2D<T> tB0 = unitTan(e12a, T(0), false), tB1 = unitTan(e01b, T(1), true);
    const gsPoint2D<T> tC0 = unitTan(e20a, T(0), false), tC1 = unitTan(e12b, T(1), true);
    const gsPoint2D<T> tM01_in = unitTan(e01a, T(1), true), tM01_out = unitTan(e01b, T(0), false);
    const gsPoint2D<T> tM12_in = unitTan(e12a, T(1), true), tM12_out = unitTan(e12b, T(0), false);
    const gsPoint2D<T> tM20_in = unitTan(e20a, T(1), true), tM20_out = unitTan(e20b, T(0), false);
    T avgH = lsum / T(3);
    auto addNormalShoots = [&](const gsPoint2D<T>& mPt, const gsPoint2D<T>& tOut) {
        if (std::hypot(tOut.x, tOut.y) <= T(1e-12)) return;
        gsPoint2D<T> nrm(-tOut.y, tOut.x);
        for (T alpha : {T(0.18), T(0.30), T(0.45)}) {
            gsPoint2D<T> pN(mPt.x + alpha * avgH * nrm.x, mPt.y + alpha * avgH * nrm.y);
            cands.push_back(pN);
            cands.emplace_back(T(0.5) * (pN.x + O_inc.x), T(0.5) * (pN.y + O_inc.y));
        }
    };
    if (isBnd01) addNormalShoots(m01, tM01_out);
    if (isBnd12) addNormalShoots(m12, tM12_out);
    if (isBnd20) addNormalShoots(m20, tM20_out);

    const T radToDeg = T(180) / T(3.14159265358979323846);
    auto angVec = [&](const gsPoint2D<T>& u, const gsPoint2D<T>& v) -> T {
        T nu = std::hypot(u.x, u.y), nv = std::hypot(v.x, v.y);
        if (nu <= T(1e-14) || nv <= T(1e-14)) return T(0);
        T a = std::atan2(u.x * v.y - u.y * v.x, u.x * v.x + u.y * v.y) * radToDeg;
        if (a <= T(0)) a += T(360);
        return a;
    };
    auto vecDir = [&](const gsPoint2D<T>& P, const gsPoint2D<T>& Q) -> gsPoint2D<T> {
        T dx = Q.x - P.x, dy = Q.y - P.y, nm = std::hypot(dx, dy);
        return (nm > T(1e-14)) ? gsPoint2D<T>(dx / nm, dy / nm) : gsPoint2D<T>(0, 0);
    };

    bool found = false;
    T bestScore = -T(1e9);
    for (const auto& Oc : cands) {
        gsPoint2D<T> q0[4] = {A, m01, Oc, m20};
        gsPoint2D<T> q1[4] = {B, m12, Oc, m01};
        gsPoint2D<T> q2[4] = {C, m20, Oc, m12};
        if (!isStrictlyConvexQuad(q0) || !isStrictlyConvexQuad(q1) || !isStrictlyConvexQuad(q2))
            continue;
        gsPoint2D<T> s01 = vecDir(Oc, m01), ms01(-s01.x, -s01.y);
        gsPoint2D<T> s12 = vecDir(Oc, m12), ms12(-s12.x, -s12.y);
        gsPoint2D<T> s20 = vecDir(Oc, m20), ms20(-s20.x, -s20.y);
        T aVals[12] = {
            angVec(tA0, tA1), angVec(ms01, tM01_in),  angVec(s20, s01), angVec(tM20_out, ms20),
            angVec(tB0, tB1), angVec(ms12, tM12_in),  angVec(s01, s12), angVec(tM01_out, ms01),
            angVec(tC0, tC1), angVec(ms20, tM20_in),  angVec(s12, s20), angVec(tM12_out, ms12)
        };
        T mnA = T(360), mxA = T(0);
        for (T av : aVals) {
            if (av <= T(1e-6) || av >= T(359)) { mnA = T(0); mxA = T(360); break; }
            mnA = std::min(mnA, av);
            mxA = std::max(mxA, av);
        }
        T effMin = std::min(mnA, T(180) - mxA);
        T sc = effMin + T(0.15) * mnA - T(0.15) * mxA;
        if (!found || sc > bestScore) {
            bestScore = sc;
            outO = Oc;
            found = true;
        }
    }
    return found;
}


} // namespace internal_splitter
} // namespace gismo
