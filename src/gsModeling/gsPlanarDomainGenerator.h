/** @file gsPlanarDomainGenerator.h
    @brief Compact 2D Planar Boundary Domain & Global-Distortion Generator in G+Smo.

    Generates planar B-spline boundary loops, multi-hole domains (gsPlanarDomain),
    and non-bilinear multi-patch geometries (gsMultiPatch) with global diffeomorphic
    mapping distortions (shear-waves, swirl, bend-taper, radial pinch, B-spline FFD
    cages, and custom gsFunction / gsFunctionExpr warps) that preserve spline B-Rep
    compatibility with gsPatchSplitter and AS-G1 reparametrization.

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gismo.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsMultiPatch.h>
#include <gsIO/gsFileManager.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsNurbsCreator.h>
#include <gsNurbs/gsTensorBSpline.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace gismo {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// gsPoint2D: Lightweight 2D Euclidean Vector Helper
// ============================================================================

template <class T = real_t>
struct gsPoint2D {
    T x = 0;
    T y = 0;

    gsPoint2D() = default;
    gsPoint2D(T x_, T y_) : x(x_), y(y_) {}

    gsPoint2D operator+(const gsPoint2D& o) const { return {x + o.x, y + o.y}; }
    gsPoint2D operator-(const gsPoint2D& o) const { return {x - o.x, y - o.y}; }
    gsPoint2D operator*(T s) const { return {x * s, y * s}; }
    gsPoint2D operator/(T s) const { return {x / s, y / s}; }
    gsPoint2D operator-() const { return {-x, -y}; }

    bool operator==(const gsPoint2D& o) const {
        return std::abs(x - o.x) < T(1e-11) && std::abs(y - o.y) < T(1e-11);
    }
    bool operator!=(const gsPoint2D& o) const { return !(*this == o); }

    T dot(const gsPoint2D& o) const { return x * o.x + y * o.y; }
    T cross(const gsPoint2D& o) const { return x * o.y - y * o.x; }
    T normSq() const { return x * x + y * y; }
    T norm() const { return std::hypot(x, y); }
    T dist(const gsPoint2D& o) const { return (*this - o).norm(); }
    T distSq(const gsPoint2D& o) const { return (*this - o).normSq(); }

    gsPoint2D normalized() const {
        T n = norm();
        return n > T(1e-14) ? *this / n : gsPoint2D(0, 0);
    }
    gsPoint2D normalCCW() const { return {-y, x}; }
    gsPoint2D normalCW() const { return {y, -x}; }
    gsPoint2D rotated(T a) const {
        T c = std::cos(a), s = std::sin(a);
        return {x * c - y * s, x * s + y * c};
    }
};

template <class T>
inline gsPoint2D<T> operator*(T s, const gsPoint2D<T>& p) { return p * s; }

// ============================================================================
// gsTriMesh: Triangulated 2D Mesh Representation
// ============================================================================

template <class T = real_t>
struct gsTriMesh {
    struct Triangle {
        int v0 = 0, v1 = 0, v2 = 0;
        Triangle() = default;
        Triangle(int a, int b, int c) : v0(a), v1(b), v2(c) {}
    };

    std::vector<gsPoint2D<T>> vertices;
    std::vector<Triangle> elements;

    T area() const {
        T total = 0;
        for (const auto& t : elements) {
            const auto &a = vertices[t.v0], &b = vertices[t.v1], &c = vertices[t.v2];
            total += T(0.5) * std::abs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
        }
        return total;
    }
};

// ============================================================================
// Flat-Array Earcut Triangulator & Polygon Geometry Helpers
// ============================================================================

namespace internal_earcut {

template <class T>
inline T signedArea(const std::vector<gsPoint2D<T>>& pts) {
    T acc = 0;
    const int n = static_cast<int>(pts.size());
    for (int i = 0; i < n; ++i)
        acc += pts[i].cross(pts[(i + 1) % n]);
    return T(0.5) * acc;
}

template <class T>
inline bool pointInPolygon(const gsPoint2D<T>& pt, const std::vector<gsPoint2D<T>>& poly) {
    bool inside = false;
    const int n = static_cast<int>(poly.size());
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const auto &a = poly[i], &b = poly[j];
        if (((a.y > pt.y) != (b.y > pt.y)) &&
            (pt.x < (b.x - a.x) * (pt.y - a.y) / (b.y - a.y + T(1e-30)) + a.x))
            inside = !inside;
    }
    return inside;
}

template <class T>
inline T distPointToSegment(const gsPoint2D<T>& p, const gsPoint2D<T>& a, const gsPoint2D<T>& b) {
    gsPoint2D<T> ab = b - a;
    T len2 = ab.normSq();
    if (len2 < T(1e-24)) return p.dist(a);
    T t = std::max(T(0), std::min(T(1), (p - a).dot(ab) / len2));
    return p.dist(a + ab * t);
}

template <class T>
inline T distToPolygon(const gsPoint2D<T>& p, const std::vector<gsPoint2D<T>>& poly) {
    T best = T(1e30);
    const int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i)
        best = std::min(best, distPointToSegment(p, poly[i], poly[(i + 1) % n]));
    return best;
}

template <class T>
inline bool segmentsIntersect(const gsPoint2D<T>& p1, const gsPoint2D<T>& p2,
                              const gsPoint2D<T>& p3, const gsPoint2D<T>& p4) {
    if (p1.distSq(p3) < T(1e-12) || p1.distSq(p4) < T(1e-12) ||
        p2.distSq(p3) < T(1e-12) || p2.distSq(p4) < T(1e-12))
        return false;
    auto ccw = [](const gsPoint2D<T>& a, const gsPoint2D<T>& b, const gsPoint2D<T>& c) {
        return (c.y - a.y) * (b.x - a.x) > (b.y - a.y) * (c.x - a.x);
    };
    return (ccw(p1, p3, p4) != ccw(p2, p3, p4)) && (ccw(p1, p2, p3) != ccw(p1, p2, p4));
}

template <class T>
struct Node {
    int v;
    T x, y;
    int prev, next;
    bool steiner = false;
};

template <class T>
inline T orient(const Node<T>& p, const Node<T>& q, const Node<T>& r) {
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
}

template <class T>
inline bool inTri(T ax, T ay, T bx, T by, T cx, T cy, T px, T py) {
    return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 &&
           (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 &&
           (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
}

template <class T>
inline int filterRing(std::vector<Node<T>>& ns, int start, int end = -1) {
    if (start < 0) return -1;
    if (end < 0) end = start;
    int p = start;
    bool again = true;
    while (again || p != end) {
        again = false;
        const auto &cur = ns[p], &nxt = ns[cur.next], &prv = ns[cur.prev];
        bool dup = std::abs(cur.x - nxt.x) < T(1e-11) && std::abs(cur.y - nxt.y) < T(1e-11);
        if (!cur.steiner && (dup || std::abs(orient(prv, cur, nxt)) < T(1e-14))) {
            ns[cur.next].prev = cur.prev;
            ns[cur.prev].next = cur.next;
            p = end = cur.prev;
            if (p == ns[p].next) return -1;
            again = true;
        } else {
            p = cur.next;
        }
    }
    return end;
}

template <class T>
inline bool locallyInside(const std::vector<Node<T>>& ns, int a, int b) {
    const auto &na = ns[a], &nb = ns[b], &prv = ns[na.prev], &nxt = ns[na.next];
    return orient(prv, na, nxt) < 0
               ? (orient(na, nb, nxt) >= 0 && orient(na, prv, nb) >= 0)
               : (orient(na, nb, prv) < 0 || orient(na, nxt, nb) < 0);
}

template <class T>
inline int splitRing(std::vector<Node<T>>& ns, int a, int b) {
    int a2 = static_cast<int>(ns.size());
    ns.push_back({ns[a].v, ns[a].x, ns[a].y, 0, 0, false});
    int b2 = static_cast<int>(ns.size());
    ns.push_back({ns[b].v, ns[b].x, ns[b].y, 0, 0, false});

    int an = ns[a].next, bp = ns[b].prev;
    ns[a].next = b;   ns[b].prev = a;
    ns[a2].next = an; ns[an].prev = a2;
    ns[b2].next = a2; ns[a2].prev = b2;
    ns[bp].next = b2; ns[b2].prev = bp;
    return b2;
}

template <class T>
inline bool isEar(const std::vector<Node<T>>& ns, int ear) {
    const auto &b = ns[ear], &a = ns[b.prev], &c = ns[b.next];
    if (orient(a, b, c) >= 0) return false;
    int p = ns[b.next].next;
    while (p != b.prev) {
        const auto& np = ns[p];
        if (inTri(a.x, a.y, b.x, b.y, c.x, c.y, np.x, np.y) &&
            orient(ns[np.prev], np, ns[np.next]) >= 0)
            return false;
        p = np.next;
    }
    return true;
}

template <class T>
inline void earcutRing(std::vector<Node<T>>& ns, int ear,
                       std::vector<typename gsTriMesh<T>::Triangle>& tris, int pass = 0) {
    if (ear < 0) return;
    int stop = ear, curr = ear;
    while (ns[curr].prev != ns[curr].next) {
        int prv = ns[curr].prev, nxt = ns[curr].next;
        if (isEar(ns, curr)) {
            tris.emplace_back(ns[prv].v, ns[curr].v, ns[nxt].v);
            ns[nxt].prev = prv;
            ns[prv].next = nxt;
            curr = stop = ns[nxt].next;
            continue;
        }
        curr = nxt;
        if (curr == stop) {
            if (pass == 0) earcutRing(ns, filterRing(ns, curr), tris, 1);
            break;
        }
    }
}

} // namespace internal_earcut

// ============================================================================
// gsPlanarWarp: Bijective 2D Diffeomorphisms (det DPhi > 0 Guaranteed)
// ============================================================================

/**
 * @brief Smooth global 2D diffeomorphism Phi : R^2 -> R^2 for domain warping.
 *
 * Uses compositions of volume-preserving nonlinear shears (det J = 1), radial
 * swirls (det J = 1), parabolic bends (det J = 1), and bounded smooth tapers
 * (det J > 0.25) so warped spline loops never fold or self-intersect.
 */
template <class T = real_t>
class gsPlanarWarp : public gsFunction<T> {
public:
    enum class Kind { None, Wave, Swirl, Bend, Pinch, Shear, Weird };

    gsPlanarWarp(const std::string& type = "none", T strength = T(0.75),
                 T frequency = T(2.6), unsigned seed = 42,
                 gsPoint2D<T> center = {0, 0}, T radius = T(1.0))
        : m_strength(strength), m_freq(frequency), m_center(center),
          m_scale(std::max(radius, T(1e-3)))
    {
        std::string t = type;
        std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        std::mt19937 rng(seed);

        if (t.empty() || t == "none" || std::abs(strength) < T(1e-12)) {
            m_kind = Kind::None;
        } else if (t == "wave" || t == "sin" || t == "ripple") {
            m_kind = Kind::Wave;
        } else if (t == "swirl" || t == "vortex" || t == "twist" || t == "spiral") {
            m_kind = Kind::Swirl;
        } else if (t == "bend" || t == "banana" || t == "taper" || t == "arch") {
            m_kind = Kind::Bend;
        } else if (t == "pinch" || t == "barrel" || t == "bulge" ||
                   t == "hourglass" || t == "peanut" || t == "dumbbell") {
            m_kind = Kind::Pinch;
        } else if (t == "shear" || t == "skew") {
            m_kind = Kind::Shear;
        } else {
            // "weird", "random", or any composite stress-test mode
            m_kind = Kind::Weird;
            std::uniform_real_distribution<T> uAng(T(0), T(2.0 * M_PI));
            std::uniform_real_distribution<T> uSign(T(-1.0), T(1.0));
            for (int k = 0; k < 3; ++k) {
                m_angles[k] = uAng(rng);
                m_phases[k] = uAng(rng);
                m_amps[k]   = (T(0.55) + T(0.45) * std::abs(uSign(rng))) * (uSign(rng) >= 0 ? T(1) : T(-1));
            }
            m_vortexCtr  = {T(0.18) * uSign(rng), T(0.18) * uSign(rng)};
            m_pinchPh[0] = uAng(rng);
            m_pinchPh[1] = uAng(rng);
            m_swirl = (T(1.1) + T(0.7) * std::abs(uSign(rng))) * (uSign(rng) >= 0 ? T(1) : T(-1));
            m_taper = T(0.42) * uSign(rng);
            m_rot   = T(0.45) * uSign(rng);
        }
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 2; }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override {
        GISMO_ASSERT(u.rows() >= 2, "gsPlanarWarp expects 2xN input points.");
        const index_t N = u.cols();
        result.resize(2, N);
        if (m_kind == Kind::None) {
            result = u.topRows(2);
            return;
        }

        const T L = m_scale;
        const T s = m_strength;
        const T w = m_freq;

        for (index_t j = 0; j < N; ++j) {
            // Normalized coordinates relative to domain center & scale
            T x = (u(0, j) - m_center.x) / L;
            T y = (u(1, j) - m_center.y) / L;

            switch (m_kind) {
            case Kind::Wave: {
                // Composed horizontal + vertical nonlinear shears: det J == 1 everywhere
                x += T(0.42) * s * std::sin(w * y + T(0.35));
                y += T(0.42) * s * std::sin(w * x - T(0.55));
                x += T(0.18) * s * std::cos(T(1.75) * w * y + T(1.1));
                break;
            }
            case Kind::Swirl: {
                // Steep off-center vortex twist + nonlinear shear: det J == 1 everywhere
                T dx = x - T(0.12), dy = y + T(0.08);
                T r2 = dx * dx + dy * dy;
                T th = T(2.35) * s * std::exp(-T(1.45) * r2);
                T c = std::cos(th), sn = std::sin(th);
                x = T(0.12) + c * dx - sn * dy;
                y = T(-0.08) + sn * dx + c * dy;
                y += T(0.22) * s * std::sin(w * x);
                break;
            }
            case Kind::Bend: {
                // Skew + deep parabolic/cubic arch + exponential keystone flare: det J > 0.63 everywhere
                x += T(0.28) * s * y;
                y += s * (T(0.88) * (x * x - T(0.35)) + T(0.18) * x * x * x);
                x *= std::exp(T(0.45) * std::tanh(s) * std::tanh(T(1.1) * y));
                break;
            }
            case Kind::Pinch: {
                // Exponential hourglass/peanut necking + bulging + sinusoidal shear: det J > 0 everywhere
                T amp = T(0.58) * std::tanh(T(1.1) * s);
                y *= std::exp(-amp * std::cos(T(2.2) * x));
                x *= std::exp(T(0.32) * std::tanh(T(1.1) * s) * std::sin(T(2.0) * y));
                y += T(0.26) * s * std::sin(w * x);
                break;
            }
            case Kind::Shear: {
                // Strong cubic + sinusoidal S-shear: det J == 1 everywhere
                x += s * (T(0.55) * y + T(0.28) * y * y * y);
                y += T(0.38) * s * std::sin(w * x);
                break;
            }
            case Kind::Weird: {
                // 3 oblique shear-waves + exponential necking + off-center vortex + taper + arch
                for (int k = 0; k < 3; ++k) {
                    T ca = std::cos(m_angles[k]), sa = std::sin(m_angles[k]);
                    T perp = -sa * x + ca * y;
                    T disp = T(0.32) * s * m_amps[k] * std::sin(w * (T(0.85) + T(0.4) * k) * perp + m_phases[k]);
                    x += disp * ca;
                    y += disp * sa;
                }
                y *= std::exp(T(0.36) * std::tanh(s) * std::cos(T(2.1) * x + m_pinchPh[0]));
                x *= std::exp(T(0.28) * std::tanh(s) * std::sin(T(1.9) * y + m_pinchPh[1]));
                T vx = x - m_vortexCtr.x, vy = y - m_vortexCtr.y;
                T r2 = vx * vx + vy * vy;
                T th = m_swirl * s * std::exp(-T(1.35) * r2) + m_rot * std::tanh(s);
                T c = std::cos(th), sn = std::sin(th);
                T nx = m_vortexCtr.x + c * vx - sn * vy;
                T ny = m_vortexCtr.y + sn * vx + c * vy;
                x = nx * std::exp(m_taper * std::tanh(s) * std::tanh(T(1.0) * ny));
                y = ny + T(0.38) * s * (nx * nx - T(0.3));
                break;
            }
            case Kind::None:
                break;
            }

            result(0, j) = m_center.x + L * x;
            result(1, j) = m_center.y + L * y;
        }
    }

private:
    Kind m_kind = Kind::None;
    T m_strength = T(0.75);
    T m_freq = T(2.6);
    gsPoint2D<T> m_center{0, 0};
    T m_scale = T(1.0);
    T m_angles[3] = {0, 1, 2};
    T m_phases[3] = {0, 0, 0};
    T m_amps[3]   = {1, 1, 1};
    gsPoint2D<T> m_vortexCtr{0, 0};
    T m_pinchPh[2] = {0, 0};
    T m_swirl = T(0.5), m_taper = T(0.2), m_rot = T(0.15);
};

// ============================================================================
// gsPlanarDomainGenerator: Core Domain, Distortion & B-Rep Generator
// ============================================================================

template <class T = real_t>
class gsPlanarDomainGenerator {
public:
    using Point2D = gsPoint2D<T>;
    using TriMesh = gsTriMesh<T>;

    struct DistortionParams {
        std::string type = "none"; // none, wave, swirl, bend, pinch, shear, weird, ffd
        T strength       = T(0.75);
        T frequency      = T(2.6);
        int refine       = 3;      // uniform B-spline refinements before Greville interpolation
        unsigned seed    = 42;

        bool active() const {
            return !type.empty() && type != "none" && std::abs(strength) > T(1e-12);
        }
    };

    struct Parameters {
        int num_corners        = 4;
        bool convex            = true;
        bool star_shaped       = true;
        int num_holes          = 0;
        std::string hole_shape = "mixed";
        T hole_radius_ratio    = T(0.18);
        int degree             = 3;
        bool curved            = false;
        T curvature            = T(0.3);
        int edge_waves         = 0;     // >0 adds S-shaped / multi-inflection boundary waves
        T scale                = T(1.0);
        T randomness           = T(0.0);
        unsigned seed          = 42;
        DistortionParams distortion;
    };

    // ------------------------------------------------------------------------
    // B-Spline Curve & Loop Constructors
    // ------------------------------------------------------------------------

    static gsBSpline<T>* makeClampedBSpline(const std::vector<Point2D>& cps, int degree = 3) {
        const int n = static_cast<int>(cps.size());
        int p = std::min(degree, n - 1);
        GISMO_ASSERT(p >= 1, "At least 2 control points are required.");
        gsKnotVector<T> kv(0.0, 1.0, n - (p + 1), p + 1, 1, p);
        gsMatrix<T> coefs(n, 2);
        for (int i = 0; i < n; ++i) {
            coefs(i, 0) = cps[i].x;
            coefs(i, 1) = cps[i].y;
        }
        auto* bsp = new gsBSpline<T>(kv, coefs);
        if (p < degree) bsp->degreeElevate(degree - p);
        return bsp;
    }

    static gsCurveLoop<T>* makeCurveLoop(const std::vector<std::vector<Point2D>>& curves_cps,
                                         int degree = 3, bool is_outer = true) {
        auto* loop = new gsCurveLoop<T>();
        for (const auto& cps : curves_cps)
            loop->insertCurve(makeClampedBSpline(cps, degree));
        ensureOrientation(loop, is_outer);
        return loop;
    }

    static std::vector<Point2D> sampleCurve(const gsCurve<T>& curve, int num_samples = 30) {
        num_samples = std::max(2, num_samples);
        gsMatrix<T> pr = curve.parameterRange();
        gsMatrix<T> u(1, num_samples), eval_pts;
        for (int i = 0; i < num_samples; ++i)
            u(0, i) = pr(0, 0) + (pr(0, 1) - pr(0, 0)) * (T(i) / T(num_samples - 1));
        curve.eval_into(u, eval_pts);

        std::vector<Point2D> pts(num_samples);
        for (int i = 0; i < num_samples; ++i)
            pts[i] = {eval_pts(0, i), eval_pts(1, i)};
        return pts;
    }

    static std::vector<Point2D> sampleLoop(const gsCurveLoop<T>& loop, int num_samples_per_curve = 30) {
        std::vector<Point2D> poly;
        const int nc = loop.numCurves();
        poly.reserve(nc * num_samples_per_curve);
        for (int c = 0; c < nc; ++c) {
            auto pts = sampleCurve(loop.curve(c), num_samples_per_curve);
            int take = static_cast<int>(pts.size()) - 1;
            for (int i = 0; i < take; ++i) poly.push_back(pts[i]);
        }
        return poly;
    }

    static T computeSignedArea(const gsCurveLoop<T>& loop, int samples_per_curve = 30) {
        return internal_earcut::signedArea(sampleLoop(loop, samples_per_curve));
    }

    static void ensureOrientation(gsCurveLoop<T>* loop, bool should_be_ccw, int samples_per_curve = 30) {
        if ((computeSignedArea(*loop, samples_per_curve) > 0) != should_be_ccw)
            loop->reverse();
    }

    /**
     * @brief Extract C^0 corner vertices where consecutive curves meet at a non-smooth angle.
     * @param sharp_only If true (default), filters out G^1-smooth joints (turn angle < 5 deg).
     */
    static std::vector<Point2D> extractCorners(const gsPlanarDomain<T>& domain,
                                               bool sharp_only = true,
                                               T smooth_tol_deg = T(5.0)) {
        std::vector<Point2D> corners;
        const auto& outer = domain.outer();
        const int nc = outer.numCurves();
        if (nc <= 1) return corners;

        const T cosTol = std::cos(smooth_tol_deg * T(M_PI / 180.0));
        gsMatrix<T> uEnd(1, 1), uStart(1, 1), d0, d1, pEnd;
        for (int i = 0; i < nc; ++i) {
            const auto& c0 = outer.curve(i);
            const auto& c1 = outer.curve((i + 1) % nc);
            gsMatrix<T> pr0 = c0.parameterRange(), pr1 = c1.parameterRange();
            uEnd(0, 0)   = pr0(0, 1);
            uStart(0, 0) = pr1(0, 0);

            if (sharp_only) {
                c0.deriv_into(uEnd, d0);
                c1.deriv_into(uStart, d1);
                Point2D t0(d0(0, 0), d0(1, 0)), t1(d1(0, 0), d1(1, 0));
                if (t0.normalized().dot(t1.normalized()) > cosTol)
                    continue; // G^1 smooth joint
            }

            if (c0.coefs().rows() > 0) {
                const index_t last = c0.coefs().rows() - 1;
                corners.emplace_back(c0.coefs()(last, 0), c0.coefs()(last, 1));
            } else {
                c0.eval_into(uEnd, pEnd);
                corners.emplace_back(pEnd(0, 0), pEnd(1, 0));
            }
        }
        return corners;
    }

    // ------------------------------------------------------------------------
    // Global Diffeomorphic Spline Distortion (B-Rep & MultiPatch)
    // ------------------------------------------------------------------------

    /**
     * @brief Build a bijective 2D B-spline Free-Form Deformation (FFD) cage over a bounding box.
     * Automatically scales down control point perturbations until min det(J) >= 0.12.
     */
    static typename gsTensorBSpline<2, T>::uPtr makeFFDCage(const gsMatrix<T>& bbox,
                                                            T strength = T(0.75),
                                                            unsigned seed = 42,
                                                            int gridN = 6,
                                                            int degree = 3) {
        gridN = std::max(gridN, degree + 1);
        const T x0 = bbox(0, 0), x1 = bbox(0, 1), y0 = bbox(1, 0), y1 = bbox(1, 1);
        const T wx = std::max(x1 - x0, T(1e-3)), wy = std::max(y1 - y0, T(1e-3));
        const T padX = T(0.15) * wx, padY = T(0.15) * wy;

        gsKnotVector<T> kvx(x0 - padX, x1 + padX, gridN - degree - 1, degree + 1, 1, degree);
        gsKnotVector<T> kvy(y0 - padY, y1 + padY, gridN - degree - 1, degree + 1, 1, degree);
        gsMatrix<T> gx = kvx.greville(), gy = kvy.greville();

        gsMatrix<T> baseCoefs(gridN * gridN, 2), delta(gridN * gridN, 2);
        delta.setZero();
        std::mt19937 rng(seed);
        std::uniform_real_distribution<T> dist(T(-1.0), T(1.0));
        const T hx = (wx + 2 * padX) / T(gridN - 1);
        const T hy = (wy + 2 * padY) / T(gridN - 1);

        // Global affine skew + strong interior FFD control point jitter
        const T skewX = T(0.28) * dist(rng), skewY = T(0.28) * dist(rng);
        const T cx = T(0.5) * (x0 + x1), cy = T(0.5) * (y0 + y1);

        for (int j = 0; j < gridN; ++j) {
            for (int i = 0; i < gridN; ++i) {
                const int idx = i + j * gridN;
                const T px = gx(0, i), py = gy(0, j);
                baseCoefs(idx, 0) = px;
                baseCoefs(idx, 1) = py;
                delta(idx, 0) = skewX * (py - cy);
                delta(idx, 1) = skewY * (px - cx);
                if (i > 0 && i < gridN - 1 && j > 0 && j < gridN - 1) {
                    delta(idx, 0) += T(1.25) * hx * dist(rng);
                    delta(idx, 1) += T(1.25) * hy * dist(rng);
                }
            }
        }

        // Backtrack amplitude until cage Jacobian is strictly positive over the domain bbox
        T alpha = strength;
        gsMatrix<T> samplePts(2, 30 * 30), J;
        const T chkX0 = x0 - T(0.05) * wx, chkX1 = x1 + T(0.05) * wx;
        const T chkY0 = y0 - T(0.05) * wy, chkY1 = y1 + T(0.05) * wy;
        for (int j = 0; j < 30; ++j)
            for (int i = 0; i < 30; ++i) {
                samplePts(0, i + 30 * j) = chkX0 + (chkX1 - chkX0) * (T(i) / T(29));
                samplePts(1, i + 30 * j) = chkY0 + (chkY1 - chkY0) * (T(j) / T(29));
            }

        for (int it = 0; it < 20; ++it) {
            typename gsTensorBSpline<2, T>::uPtr cage(
                new gsTensorBSpline<2, T>(kvx, kvy, baseCoefs + alpha * delta));
            bool ok = true;
            for (index_t k = 0; k < samplePts.cols(); ++k) {
                cage->jacobian_into(samplePts.col(k), J);
                if (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0) < T(0.12)) {
                    ok = false;
                    break;
                }
            }
            if (ok) return cage;
            alpha *= T(0.82);
        }
        return typename gsTensorBSpline<2, T>::uPtr(new gsTensorBSpline<2, T>(kvx, kvy, baseCoefs));
    }

    /**
     * @brief Push a curve through a 2D mapping Phi via knot refinement + Greville interpolation.
     * Exact endpoints u=0 and u=1 are Greville points, preserving exact C^0 loop watertightness.
     */
    static gsBSpline<T>* distortCurve(const gsCurve<T>& crv, const gsFunction<T>& phi,
                                      int num_refine = 2, int target_degree = 3) {
        std::unique_ptr<gsGeometry<T>> cloned(crv.clone().release());
        auto* rawBsp = dynamic_cast<gsBSpline<T>*>(cloned.get());
        std::unique_ptr<gsBSpline<T>> bsp;
        if (rawBsp) {
            cloned.release();
            bsp.reset(rawBsp);
            bsp->knots().affineTransformTo(0.0, 1.0);
            if (bsp->degree() < target_degree)
                bsp->degreeElevate(target_degree - bsp->degree());
        } else {
            auto pts = sampleCurve(crv, std::max(4, target_degree + 1));
            bsp.reset(makeClampedBSpline(pts, target_degree));
        }

        if (num_refine > 0) bsp->uniformRefine((1 << std::min(num_refine, 6)) - 1);

        gsMatrix<T> anchors, pts, warped;
        bsp->basis().anchors_into(anchors);
        bsp->eval_into(anchors, pts);
        phi.eval_into(pts, warped);

        auto geom = bsp->basis().interpolateAtAnchors(warped);
        return static_cast<gsBSpline<T>*>(geom.release());
    }

    /**
     * @brief Apply an arbitrary 2D gsFunction (e.g. gsPlanarWarp, gsFunctionExpr, FFD cage)
     * to every loop and hole of a gsPlanarDomain, returning a new spline B-Rep domain.
     */
    static gsPlanarDomain<T>* distortDomain(const gsPlanarDomain<T>& domain,
                                            const gsFunction<T>& phi,
                                            int num_refine = 2,
                                            int target_degree = 3) {
        auto warpLoop = [&](const gsCurveLoop<T>& src, bool is_outer) {
            auto* dst = new gsCurveLoop<T>();
            for (int c = 0; c < src.numCurves(); ++c) {
                int deg = std::max(target_degree, static_cast<int>(src.curve(c).degree()));
                dst->insertCurve(distortCurve(src.curve(c), phi, num_refine, deg));
            }
            ensureOrientation(dst, is_outer);
            return dst;
        };

        auto* out = new gsPlanarDomain<T>(warpLoop(domain.outer(), true));
        for (int h = 1; h <= domain.numHoles(); ++h)
            out->insertHole(warpLoop(domain.loop(h), false));
        return out;
    }

    /**
     * @brief Apply a named DistortionParams preset (wave, swirl, bend, pinch, shear, weird, ffd)
     * to a gsPlanarDomain.
     */
    static gsPlanarDomain<T>* distortDomain(const gsPlanarDomain<T>& domain,
                                            const DistortionParams& dp,
                                            int target_degree = 3) {
        if (!dp.active()) return new gsPlanarDomain<T>(domain);

        const int effRefine = dp.refine + (std::abs(dp.strength) >= T(0.8) ? 1 : 0);
        gsMatrix<T> bb = domain.boundingBox();
        std::string t = dp.type;
        std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        if (t == "ffd" || t == "cage") {
            auto cage = makeFFDCage(bb, dp.strength, dp.seed, 6, 3);
            return distortDomain(domain, *cage, effRefine, target_degree);
        }

        Point2D ctr(T(0.5) * (bb(0, 0) + bb(0, 1)), T(0.5) * (bb(1, 0) + bb(1, 1)));
        T rad = T(0.5) * std::max(bb(0, 1) - bb(0, 0), bb(1, 1) - bb(1, 0));
        gsPlanarWarp<T> warp(dp.type, dp.strength, dp.frequency, dp.seed, ctr, rad);
        return distortDomain(domain, warp, effRefine, target_degree);
    }

    /**
     * @brief Apply an arbitrary 2D gsFunction warp to all patches of a gsMultiPatch,
     * preserving exact C^0 interface conformity across shared Greville anchors.
     */
    static gsMultiPatch<T>* distortMultiPatch(const gsMultiPatch<T>& mp,
                                              const gsFunction<T>& phi,
                                              int num_refine = 1,
                                              int target_degree = 3) {
        auto* out = new gsMultiPatch<T>();
        for (size_t p = 0; p < mp.nPatches(); ++p) {
            std::unique_ptr<gsGeometry<T>> g(mp.patch(p).clone().release());
            if (g->degree(0) < target_degree) g->degreeElevate(target_degree - g->degree(0), 0);
            if (g->degree(1) < target_degree) g->degreeElevate(target_degree - g->degree(1), 1);
            if (num_refine > 0) g->uniformRefine((1 << std::min(num_refine, 5)) - 1);

            gsMatrix<T> anchors, pts, warped;
            g->basis().anchors_into(anchors);
            g->eval_into(anchors, pts);
            phi.eval_into(pts, warped);
            out->addPatch(g->basis().interpolateAtAnchors(warped));
        }
        out->computeTopology();
        return out;
    }

    static gsMultiPatch<T>* distortMultiPatch(const gsMultiPatch<T>& mp,
                                              const DistortionParams& dp,
                                              int target_degree = 3) {
        if (!dp.active()) return new gsMultiPatch<T>(mp);

        gsMatrix<T> bb;
        mp.boundingBox(bb);
        std::string t = dp.type;
        std::transform(t.begin(), t.end(), t.begin(), ::tolower);
        if (t == "ffd" || t == "cage") {
            auto cage = makeFFDCage(bb, dp.strength, dp.seed, 6, 3);
            return distortMultiPatch(mp, *cage, std::max(2, dp.refine - 1), target_degree);
        }

        Point2D ctr(T(0.5) * (bb(0, 0) + bb(0, 1)), T(0.5) * (bb(1, 0) + bb(1, 1)));
        T rad = T(0.5) * std::max(bb(0, 1) - bb(0, 0), bb(1, 1) - bb(1, 0));
        gsPlanarWarp<T> warp(dp.type, dp.strength, dp.frequency, dp.seed, ctr, rad);
        return distortMultiPatch(mp, warp, std::max(2, dp.refine - 1), target_degree);
    }

    // ------------------------------------------------------------------------
    // Constrained Earcut Triangulation (with Rightmost-First Hole Bridging)
    // ------------------------------------------------------------------------

    static TriMesh triangulate(const gsPlanarDomain<T>& domain, int num_samples_per_curve = 25) {
        TriMesh mesh;
        auto outer_poly = sampleLoop(domain.outer(), num_samples_per_curve);
        if (outer_poly.size() < 3) return mesh;
        if (internal_earcut::signedArea(outer_poly) < 0)
            std::reverse(outer_poly.begin(), outer_poly.end());

        const int num_holes = domain.numHoles();
        std::vector<std::vector<Point2D>> holes_poly(num_holes);
        for (int h = 0; h < num_holes; ++h) {
            holes_poly[h] = sampleLoop(domain.loop(h + 1), num_samples_per_curve);
            if (internal_earcut::signedArea(holes_poly[h]) > 0)
                std::reverse(holes_poly[h].begin(), holes_poly[h].end());
        }

        std::vector<internal_earcut::Node<T>> nodes;
        nodes.reserve(outer_poly.size() + num_holes * 36);

        auto addRing = [&](const std::vector<Point2D>& poly) -> int {
            const int n = static_cast<int>(poly.size());
            if (n < 3) return -1;
            const int vOff = static_cast<int>(mesh.vertices.size());
            const int nOff = static_cast<int>(nodes.size());
            mesh.vertices.insert(mesh.vertices.end(), poly.begin(), poly.end());
            for (int i = 0; i < n; ++i)
                nodes.push_back({vOff + i, poly[i].x, poly[i].y,
                                 nOff + (i + n - 1) % n, nOff + (i + 1) % n, false});
            return internal_earcut::filterRing(nodes, nOff);
        };

        int outer_head = addRing(outer_poly);
        if (outer_head < 0) return mesh;

        struct HoleRing { int head; int maxNode; T maxX; };
        std::vector<HoleRing> hRings;
        for (const auto& hp : holes_poly) {
            int hd = addRing(hp);
            if (hd < 0) continue;
            int cur = hd, best = hd;
            do {
                if (nodes[cur].x > nodes[best].x) best = cur;
                cur = nodes[cur].next;
            } while (cur != hd);
            hRings.push_back({hd, best, nodes[best].x});
        }

        // Process holes from rightmost to leftmost so bridges never cross unprocessed holes
        std::sort(hRings.begin(), hRings.end(),
                  [](const HoleRing& a, const HoleRing& b) { return a.maxX > b.maxX; });

        std::vector<std::pair<Point2D, Point2D>> bndEdges;
        auto collectEdges = [&](const std::vector<Point2D>& poly) {
            const int n = static_cast<int>(poly.size());
            for (int i = 0; i < n; ++i) bndEdges.emplace_back(poly[i], poly[(i + 1) % n]);
        };
        collectEdges(outer_poly);
        for (const auto& hp : holes_poly) collectEdges(hp);

        for (const auto& hr : hRings) {
            // Find mutually visible bridge from rightmost hole vertex (or fallback hole vertex)
            int bestP = -1, bestH = -1;
            T bestDist2 = T(1e30);

            int hCand = hr.maxNode;
            for (int pass = 0; pass < 2 && bestP < 0; ++pass) {
                int hCur = hCand;
                do {
                    Point2D ph(nodes[hCur].x, nodes[hCur].y);
                    int pCur = outer_head;
                    do {
                        Point2D pp(nodes[pCur].x, nodes[pCur].y);
                        if (pass == 0 && pp.x < ph.x - T(1e-9)) {
                            pCur = nodes[pCur].next;
                            continue;
                        }
                        T d2 = ph.distSq(pp);
                        if (d2 < bestDist2 &&
                            internal_earcut::locallyInside(nodes, pCur, hCur) &&
                            internal_earcut::locallyInside(nodes, hCur, pCur)) {
                            bool hit = false;
                            for (const auto& e : bndEdges) {
                                if (internal_earcut::segmentsIntersect(pp, ph, e.first, e.second)) {
                                    hit = true;
                                    break;
                                }
                            }
                            if (!hit) {
                                bool insideDom = true;
                                for (T frac : {T(0.25), T(0.5), T(0.75)}) {
                                    Point2D q = pp + (ph - pp) * frac;
                                    if (!internal_earcut::pointInPolygon(q, outer_poly)) {
                                        insideDom = false;
                                        break;
                                    }
                                    for (const auto& hp : holes_poly) {
                                        if (internal_earcut::pointInPolygon(q, hp)) {
                                            insideDom = false;
                                            break;
                                        }
                                    }
                                    if (!insideDom) break;
                                }
                                if (insideDom) { bestDist2 = d2; bestP = pCur; bestH = hCur; }
                            }
                        }
                        pCur = nodes[pCur].next;
                    } while (pCur != outer_head);
                    if (pass == 0) break; // Pass 0 tests rightmost vertex only (O(N)); Pass 1 scans all if needed
                    hCur = nodes[hCur].next;
                } while (hCur != hCand);
            }

            if (bestP >= 0 && bestH >= 0) {
                bndEdges.emplace_back(Point2D(nodes[bestP].x, nodes[bestP].y),
                                      Point2D(nodes[bestH].x, nodes[bestH].y));
                outer_head = internal_earcut::splitRing(nodes, bestP, bestH);
                outer_head = internal_earcut::filterRing(nodes, outer_head);
            }
        }

        std::vector<typename TriMesh::Triangle> raw;
        internal_earcut::earcutRing(nodes, outer_head, raw);

        mesh.elements.reserve(raw.size());
        for (const auto& tri : raw) {
            Point2D c = (mesh.vertices[tri.v0] + mesh.vertices[tri.v1] + mesh.vertices[tri.v2]) / T(3.0);
            if (!internal_earcut::pointInPolygon(c, outer_poly)) continue;
            bool inHole = false;
            for (const auto& hp : holes_poly)
                if (internal_earcut::pointInPolygon(c, hp)) { inHole = true; break; }
            if (!inHole) mesh.elements.push_back(tri);
        }
        return mesh;
    }

    // ------------------------------------------------------------------------
    // Main Generator Entry Points
    // ------------------------------------------------------------------------

    static gsPlanarDomain<T>* generate(const Parameters& params) {
        std::mt19937 rng(params.seed);
        gsCurveLoop<T>* outer = nullptr;
        const int K = params.num_corners;

        if (K <= 0) {
            outer = generate0Corner(params.convex, params.star_shaped, params.degree,
                                    params.curvature, params.scale, params.randomness, rng);
        } else if (K == 1) {
            outer = generate1Corner(params.convex, params.degree, params.scale);
        } else if (K == 2) {
            outer = generate2Corner(params.convex, params.star_shaped, params.degree, params.scale);
        } else {
            outer = generateNCorner(K, params.convex, params.star_shaped, params.degree,
                                    params.curved, params.curvature, params.edge_waves,
                                    params.scale, params.randomness, rng);
        }

        ensureOrientation(outer, true);
        std::unique_ptr<gsPlanarDomain<T>> domain(new gsPlanarDomain<T>(outer));

        if (params.num_holes > 0) {
            auto holes = generateHoles(*outer, params.num_holes, params.hole_shape,
                                       params.hole_radius_ratio, params.degree,
                                       params.scale, params.randomness, rng);
            for (auto* h : holes) domain->insertHole(h);
        }

        if (params.distortion.active())
            return distortDomain(*domain, params.distortion, params.degree);

        return domain.release();
    }

    // ------------------------------------------------------------------------
    // Canonical Outer Loop Builders
    // ------------------------------------------------------------------------

    static gsCurveLoop<T>* generate0Corner(bool convex, bool star_shaped, int degree,
                                           T curvature, T scale, T randomness, std::mt19937& rng) {
        if (!star_shaped) return generateSmoothHorseshoe(scale, degree);

        if (convex && randomness < T(1e-4)) {
            return makeEllipseLoop({0, 0}, scale * T(1.3), scale * T(0.9), degree, true);
        }

        // Smooth radial star/bean/amoeba curve
        std::uniform_real_distribution<T> phDist(0, T(2.0 * M_PI));
        T ph1 = phDist(rng), ph2 = phDist(rng);
        auto r_func = [=](T th) {
            T r = convex ? T(1.0) : (T(1.0) + T(0.35) * std::cos(th) - (T(0.2) + T(0.1) * curvature) * std::sin(T(2.0) * th));
            if (randomness > T(1e-4))
                r += randomness * (T(0.18) * std::sin(T(3.0) * th + ph1) + T(0.12) * std::cos(T(5.0) * th + ph2));
            return scale * std::max(T(0.25), r);
        };
        return fitRadialClosedCurve(r_func, randomness > T(1e-4) ? 8 : 4, degree);
    }

    static gsCurveLoop<T>* generate1Corner(bool convex, int degree, T scale) {
        T L = scale * T(1.2), H = scale * T(0.85);
        std::vector<std::vector<Point2D>> curves;
        if (convex) {
            curves = {{{L, 0}, {L * T(0.4), H * T(1.1)}, {-L * T(0.7), H * T(0.8)}, {-L, 0}},
                      {{-L, 0}, {-L * T(0.7), -H * T(0.8)}, {L * T(0.4), -H * T(1.1)}, {L, 0}}};
        } else {
            curves = {{{0, -L}, {-L * T(0.9), -L * T(0.2)}, {-L * T(0.9), H * T(0.9)}, {0, H * T(0.4)}},
                      {{0, H * T(0.4)}, {L * T(0.9), H * T(0.9)}, {L * T(0.9), -L * T(0.2)}, {0, -L}}};
        }
        return makeCurveLoop(curves, degree, true);
    }

    static gsCurveLoop<T>* generate2Corner(bool convex, bool star_shaped, int degree, T scale) {
        T L = scale * T(1.2), H = scale * T(0.8);
        if (convex) {
            return makeCurveLoop({{{-L, 0}, {-L * T(0.4), H}, {L * T(0.4), H}, {L, 0}},
                                  {{L, 0}, {L * T(0.4), -H}, {-L * T(0.4), -H}, {-L, 0}}},
                                 degree, true);
        }
        T Hin = star_shaped ? H * T(0.45) : H * T(0.85);
        return makeCurveLoop({{{-L, 0}, {-L * T(0.5), H * T(1.3)}, {L * T(0.5), H * T(1.3)}, {L, 0}},
                              {{L, 0}, {L * T(0.4), Hin}, {-L * T(0.4), Hin}, {-L, 0}}},
                             degree, true);
    }

    static gsCurveLoop<T>* generateNCorner(int K, bool convex, bool star_shaped,
                                           int degree, bool curved, T curvature, int edge_waves,
                                           T scale, T randomness, std::mt19937& rng) {
        if (!star_shaped)
            return generateNonStarDomain(K, degree, curved, curvature, edge_waves, scale, randomness, rng);

        std::vector<T> angles(K);
        T d_theta = T(2.0 * M_PI) / T(K);
        std::uniform_real_distribution<T> aJit(T(-0.35), T(0.35)), rJit(T(-0.22), T(0.22));
        for (int i = 0; i < K; ++i)
            angles[i] = i * d_theta + (randomness > T(1e-4) ? aJit(rng) * d_theta * randomness : T(0));
        std::sort(angles.begin(), angles.end());

        std::vector<Point2D> corners;
        corners.reserve(K);
        for (int i = 0; i < K; ++i) {
            T r = convex ? scale
                         : ((i % 2 == 1) ? scale * (T(0.45) - T(0.2) * curvature)
                                         : scale * (T(1.1) + T(0.2) * curvature));
            if (randomness > T(1e-4)) r *= (T(1.0) + rJit(rng) * randomness);
            corners.emplace_back(r * std::cos(angles[i]), r * std::sin(angles[i]));
        }
        return buildCurvesFromCorners(corners, degree, curved, curvature, true, edge_waves);
    }

    static gsCurveLoop<T>* generateNonStarDomain(int K, int degree, bool curved, T curvature,
                                                 int edge_waves, T s, T randomness, std::mt19937& rng) {
        std::vector<Point2D> c;
        if (K == 6) {
            c = {{-s, -s}, {s, -s}, {s, 0}, {0, 0}, {0, s}, {-s, s}};
        } else if (K == 8) {
            c = {{-s * T(1.2), s}, {-s * T(1.2), -s}, {s * T(1.2), -s}, {s * T(1.2), s},
                 {s * T(0.5), s},  {s * T(0.5), 0},   {-s * T(0.5), 0}, {-s * T(0.5), s}};
        } else if (K == 12) {
            T a = s * T(0.4), b = s;
            c = {{-a, b}, {a, b}, {a, a}, {b, a}, {b, -a}, {a, -a},
                 {a, -b}, {-a, -b}, {-a, -a}, {-b, -a}, {-b, a}, {-a, a}};
        } else {
            int teeth = std::max(1, (K - 4) / 4);
            T x0 = -s * T(1.2), x1 = s * T(1.2), yb = -s * T(0.8), yt = s * T(0.8), yr = 0;
            T w = (x1 - x0) / T(2 * teeth - 1);
            c = {{x0, yt}, {x0, yb}, {x1, yb}, {x1, yt}};
            for (int i = 0; i < teeth - 1; ++i) {
                T xr = x1 - (2 * i + 1) * w, xl = xr - w;
                c.insert(c.end(), {{xr, yt}, {xr, yr}, {xl, yr}, {xl, yt}});
            }
            // If K == 16 -> teeth=3? Wait: 4 + 4*(teeth) = K -> teeth = (K-4)/4 -> 3 valleys = 4 fingers!
            // With teeth = (K-4)/4 valleys: w = (x1 - x0) / (2*teeth + 1) gives exact K = 4 + 4*teeth corners.
            if (static_cast<int>(c.size()) != K) {
                c = {{x0, yt}, {x0, yb}, {x1, yb}, {x1, yt}};
                w = (x1 - x0) / T(2 * teeth + 1);
                for (int i = 0; i < teeth; ++i) {
                    T xr = x1 - (2 * i + 1) * w, xl = xr - w;
                    c.insert(c.end(), {{xr, yt}, {xr, yr}, {xl, yr}, {xl, yt}});
                }
            }
        }

        if (randomness > T(1e-4)) {
            std::uniform_real_distribution<T> jit(T(-0.12), T(0.12));
            for (auto& pt : c) {
                pt.x += s * randomness * jit(rng);
                pt.y += s * randomness * jit(rng);
            }
        }
        return buildCurvesFromCorners(c, degree, curved, curvature, true, edge_waves);
    }

    static gsCurveLoop<T>* buildCurvesFromCorners(const std::vector<Point2D>& corners,
                                                  int degree, bool curved, T curvature,
                                                  bool is_outer = true, int edge_waves = 0) {
        const int nc = static_cast<int>(corners.size());
        std::vector<std::vector<Point2D>> curves_cps;
        curves_cps.reserve(nc);

        for (int i = 0; i < nc; ++i) {
            Point2D p0 = corners[i], p1 = corners[(i + 1) % nc];
            Point2D d = p1 - p0;
            T len = d.norm();
            Point2D nrm = (is_outer ? d.normalCCW() : d.normalCW()).normalized();
            int nCtrl = std::max(degree + 1, edge_waves > 0 ? 2 * edge_waves + degree + 1 : degree + 1);

            std::vector<Point2D> cps;
            cps.reserve(nCtrl);
            for (int j = 0; j < nCtrl; ++j) {
                T t = T(j) / T(nCtrl - 1);
                T offset = 0;
                if (len > T(1e-8) && j > 0 && j < nCtrl - 1) {
                    if (edge_waves > 0)
                        offset = len * curvature * T(0.22) * std::sin(T(edge_waves * M_PI) * t);
                    else if (curved)
                        offset = len * curvature * T(0.35) * ((nCtrl == 4) ? T(1.0) : T(4.0) * t * (T(1.0) - t));
                }
                cps.push_back(p0 + d * t + nrm * offset);
            }
            curves_cps.push_back(std::move(cps));
        }
        return makeCurveLoop(curves_cps, degree, is_outer);
    }

    // ------------------------------------------------------------------------
    // Unified Medial-Clearance Hole Placement & Single Hole Builder
    // ------------------------------------------------------------------------

    static std::vector<gsCurveLoop<T>*> generateHoles(const gsCurveLoop<T>& outer_loop,
                                                      int num_holes,
                                                      const std::string& shape,
                                                      T radius_ratio,
                                                      int degree,
                                                      T scale,
                                                      T randomness,
                                                      std::mt19937& rng) {
        std::vector<gsCurveLoop<T>*> holes;
        if (num_holes <= 0) return holes;

        auto outer_poly = sampleLoop(outer_loop, 32);
        T x0 = T(1e30), x1 = T(-1e30), y0 = T(1e30), y1 = T(-1e30);
        for (const auto& p : outer_poly) {
            x0 = std::min(x0, p.x); x1 = std::max(x1, p.x);
            y0 = std::min(y0, p.y); y1 = std::max(y1, p.y);
        }

        // Sample interior clearance field on an odd grid (includes exact symmetry axes & origin)
        struct Cand { Point2D pt; T dOuter; };
        std::vector<Cand> cands;
        const int G = 33;
        std::uniform_real_distribution<T> jit(T(-0.35), T(0.35));
        const T dx = (x1 - x0) / T(G + 1), dy = (y1 - y0) / T(G + 1);

        auto tryAdd = [&](Point2D p) {
            if (!internal_earcut::pointInPolygon(p, outer_poly)) return;
            T d = internal_earcut::distToPolygon(p, outer_poly);
            if (d > scale * T(0.03)) cands.push_back({p, d});
        };
        tryAdd({T(0.5) * (x0 + x1), T(0.5) * (y0 + y1)});
        for (int j = 1; j <= G; ++j) {
            for (int i = 1; i <= G; ++i) {
                Point2D p(x0 + i * dx, y0 + j * dy);
                if (randomness > T(1e-4)) {
                    p.x += randomness * dx * jit(rng);
                    p.y += randomness * dy * jit(rng);
                }
                tryAdd(p);
            }
        }
        if (cands.empty()) return holes;

        // Greedy farthest-point / medial-axis selection maximizing clearance to boundary and prior holes
        std::vector<Point2D> chosen;
        std::vector<T> chosenOuterDist;
        for (int h = 0; h < num_holes; ++h) {
            int bestIdx = -1;
            T bestScore = T(-1);
            for (size_t i = 0; i < cands.size(); ++i) {
                T minInter = T(1e30);
                for (const auto& c : chosen)
                    minInter = std::min(minInter, T(0.5) * cands[i].pt.dist(c));
                // Slight tie-break preference toward larger boundary clearance
                T score = std::min(cands[i].dOuter, minInter) + T(0.05) * cands[i].dOuter;
                if (minInter < scale * T(0.04)) continue;
                if (score > bestScore) { bestScore = score; bestIdx = static_cast<int>(i); }
            }
            if (bestIdx < 0) break;
            chosen.push_back(cands[bestIdx].pt);
            chosenOuterDist.push_back(cands[bestIdx].dOuter);
        }

        static const std::vector<std::string> pool = {
            "circle", "square", "triangle", "diamond", "hexagon", "star", "cross"
        };
        const T r_target = scale * radius_ratio;
        for (size_t i = 0; i < chosen.size(); ++i) {
            T safe_r = std::min(r_target, T(0.62) * chosenOuterDist[i]);
            for (size_t j = 0; j < chosen.size(); ++j)
                if (i != j) safe_r = std::min(safe_r, T(0.38) * chosen[i].dist(chosen[j]));
            if (safe_r < scale * T(0.02)) continue;

            std::string sh = (shape == "mixed") ? pool[i % pool.size()] : shape;
            T rot = (randomness > T(1e-4)) ? T(2.0 * M_PI) * jit(rng) : T(0);
            holes.push_back(makeSingleHole(chosen[i], safe_r, sh, degree, rot));
        }
        return holes;
    }

    static gsCurveLoop<T>* makeSingleHole(const Point2D& center, T radius,
                                          const std::string& shape, int degree,
                                          T rotation = T(0)) {
        std::string s = shape;
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);

        if (s == "circle" || s == "ellipse") {
            T rx = radius, ry = (s == "ellipse") ? radius * T(0.65) : radius;
            return makeEllipseLoop(center, rx, ry, degree, false, rotation);
        }

        std::vector<Point2D> pts;
        auto regPoly = [&](int K, T rScale, T aspectX, T th0) {
            for (int i = 0; i < K; ++i) {
                T th = th0 - i * T(2.0 * M_PI) / T(K); // CW
                Point2D loc(rScale * radius * aspectX * std::cos(th), rScale * radius * std::sin(th));
                pts.push_back(center + loc.rotated(rotation));
            }
        };

        if (s == "square")        regPoly(4, T(0.75) * std::sqrt(T(2.0)), T(1.0), T(0.75 * M_PI));
        else if (s == "triangle") regPoly(3, T(1.0), T(1.0), T(0.5 * M_PI));
        else if (s == "diamond")  regPoly(4, T(1.0), T(0.75), T(0.5 * M_PI));
        else if (s == "hexagon")  regPoly(6, T(1.0), T(1.0), T(0.5 * M_PI));
        else if (s == "star") {
            for (int i = 0; i < 10; ++i) {
                T th = T(0.5 * M_PI) - i * T(M_PI / 5.0);
                T r = (i % 2 == 0) ? radius : radius * T(0.48);
                pts.push_back(center + Point2D(r * std::cos(th), r * std::sin(th)).rotated(rotation));
            }
        } else if (s == "cross") {
            T a = radius * T(0.35), b = radius * T(0.9);
            std::vector<Point2D> raw = {{-a, b}, {a, b}, {a, a}, {b, a}, {b, -a}, {a, -a},
                                        {a, -b}, {-a, -b}, {-a, -a}, {-b, -a}, {-b, a}, {-a, a}};
            for (const auto& q : raw) pts.push_back(center + q.rotated(rotation));
        } else {
            return makeEllipseLoop(center, radius, radius, degree, false, rotation);
        }

        return buildCurvesFromCorners(pts, degree, false, T(0), false);
    }

    // ------------------------------------------------------------------------
    // Table-Driven Built-in Presets (with Optional Global Distortion)
    // ------------------------------------------------------------------------

    static gsPlanarDomain<T>* preset(const std::string& name, T scale = T(1.0),
                                     int degree = 3,
                                     const DistortionParams& dp = DistortionParams()) {
        std::string n = name;
        std::transform(n.begin(), n.end(), n.begin(), ::tolower);
        std::replace(n.begin(), n.end(), '-', '_');

        if (n == "wavy_ring" || n == "pinched_annulus") {
            auto r_func = [scale, n](T th) {
                return (n == "pinched_annulus")
                           ? scale * (T(1.05) + T(0.28) * std::cos(T(2.0) * th) + T(0.12) * std::sin(T(3.0) * th))
                           : scale * (T(1.1) + T(0.18) * std::sin(T(4.0) * th));
            };
            std::unique_ptr<gsPlanarDomain<T>> dom(
                new gsPlanarDomain<T>(fitRadialClosedCurve(r_func, 6, degree)));
            dom->insertHole(makeSingleHole({0, 0}, scale * T(0.42),
                                           n == "pinched_annulus" ? "ellipse" : "circle", degree, T(0.4)));
            return dp.active() ? distortDomain(*dom, dp, degree) : dom.release();
        }

        struct Spec {
            const char* key;
            int corners; bool convex; bool star; bool curved; T curv;
            int holes; const char* hShape; int waves; T rand; const char* defWarp;
        };
        static const Spec table[] = {
            // Smooth 0-corner
            {"circle",               0, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"disk",                 0, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"ellipse",              0, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"bean",                 0, 0, 1, 0, 0.30, 0, "circle", 0, 0.00, "none"},
            {"kidney",               0, 0, 1, 0, 0.30, 0, "circle", 0, 0.00, "none"},
            {"amoeba",               0, 0, 1, 1, 0.35, 0, "mixed",  0, 0.85, "wave"},
            {"horseshoe_smooth",     0, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"c_shape_smooth",       0, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"annulus",              0, 1, 1, 0, 0.00, 1, "circle", 0, 0.00, "none"},
            {"donut",                0, 1, 1, 0, 0.00, 1, "circle", 0, 0.00, "none"},
            // 1 & 2 corners
            {"teardrop",             1, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"heart",                1, 0, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"lens",                 2, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"football",             2, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"crescent",             2, 0, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"crescent_3holes",      2, 0, 1, 0, 0.00, 3, "mixed",  0, 0.00, "none"},
            {"crescent_2corner_3holes", 2, 0, 1, 0, 0.00, 3, "mixed", 0, 0.00, "none"},
            // Polygons & curved polygons
            {"triangle",             3, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"tri",                  3, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"curved_triangle",      3, 1, 1, 1, 0.35, 0, "circle", 0, 0.00, "none"},
            {"quad",                 4, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"square",               4, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"rectangle",            4, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"curved_quad",          4, 1, 1, 1, 0.30, 0, "circle", 0, 0.00, "none"},
            {"bulged_quad",          4, 1, 1, 1, 0.45, 0, "circle", 0, 0.00, "none"},
            {"indented_quad",        4, 0, 1, 1,-0.30, 0, "circle", 0, 0.00, "none"},
            {"wavy_quad",            4, 1, 1, 1, 0.35, 0, "circle", 2, 0.25, "none"},
            {"pentagon",             5, 1, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"curved_pentagon",      5, 1, 1, 1, 0.25, 0, "circle", 0, 0.00, "none"},
            {"random_polygon",       7, 0, 1, 1, 0.25, 2, "mixed",  1, 0.70, "weird"},
            // Stars
            {"star",                10, 0, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"star5",               10, 0, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"twisted_star",        10, 0, 1, 1, 0.20, 1, "circle", 0, 0.25, "swirl"},
            {"starburst",           12, 0, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"star12",              12, 0, 1, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"curved_starburst",    12, 0, 1, 1, 0.30, 0, "circle", 0, 0.00, "none"},
            {"starburst_3holes",    12, 0, 1, 0, 0.00, 3, "mixed",  0, 0.00, "none"},
            {"starburst_12corner_3holes", 12, 0, 1, 0, 0.00, 3, "mixed", 0, 0.00, "none"},
            // Non-star shapes (L-shape, U-channel, Cross, Comb)
            {"lshape",               6, 0, 0, 0, 0.00, 0, "square", 0, 0.00, "none"},
            {"l_shape",              6, 0, 0, 0, 0.00, 0, "square", 0, 0.00, "none"},
            {"curved_lshape",        6, 0, 0, 1, 0.25, 0, "square", 0, 0.00, "none"},
            {"curved_l_shape",       6, 0, 0, 1, 0.25, 0, "square", 0, 0.00, "none"},
            {"warped_lshape",        6, 0, 0, 1, 0.20, 2, "mixed",  1, 0.25, "weird"},
            {"u_channel",            8, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"uchannel",             8, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"curved_uchannel",      8, 0, 0, 1, 0.20, 0, "circle", 0, 0.00, "none"},
            {"swiss_cross",         12, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"cross",               12, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"bulged_cross",        12, 0, 0, 1, 0.25, 0, "circle", 0, 0.00, "none"},
            {"comb",                16, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"comb_16corner",       16, 0, 0, 0, 0.00, 0, "circle", 0, 0.00, "none"},
            {"curved_comb",         16, 0, 0, 1, 0.15, 0, "circle", 0, 0.00, "none"},
            // Multi-hole benchmarks
            {"plate_with_hole",      4, 1, 1, 0, 0.00, 1, "circle", 0, 0.00, "none"},
            {"plate_with_2holes",    4, 1, 1, 0, 0.00, 2, "circle", 0, 0.00, "none"},
            {"plate_5holes",         4, 1, 1, 0, 0.00, 5, "mixed",  0, 0.00, "none"},
            {"plate_with_5holes",    4, 1, 1, 0, 0.00, 5, "mixed",  0, 0.00, "none"},
            {"plate_8holes",         4, 1, 1, 0, 0.00, 8, "mixed",  0, 0.00, "none"},
            {"plate_with_8holes",    4, 1, 1, 0, 0.00, 8, "mixed",  0, 0.00, "none"},
            {"plate_4corner_8holes_weird", 4, 1, 1, 0, 0.00, 8, "mixed", 0, 0.00, "none"},
            {"lshape_with_hole",     6, 0, 0, 0, 0.00, 1, "square", 0, 0.00, "none"},
            {"curved_lshape_1hole",  6, 0, 0, 1, 0.20, 1, "square", 0, 0.00, "none"},
            {"curved_lshape_6corner_1hole", 6, 0, 0, 1, 0.20, 1, "square", 0, 0.00, "none"},
            {"lshape_3holes",        6, 0, 0, 0, 0.20, 3, "mixed",  0, 0.00, "none"},
            {"curved_lshape_3holes", 6, 0, 0, 1, 0.20, 3, "mixed",  0, 0.00, "none"},
            {"uchannel_2holes",      8, 0, 0, 0, 0.00, 2, "circle", 0, 0.00, "none"},
            {"uchannel_8corner_2holes", 8, 0, 0, 0, 0.00, 2, "circle", 0, 0.00, "none"},
            {"uchannel_4holes",      8, 0, 0, 0, 0.20, 4, "mixed",  0, 0.00, "none"},
            {"curved_uchannel_4holes", 8, 0, 0, 1, 0.20, 4, "mixed", 0, 0.00, "none"},
            {"cross_5holes",        12, 0, 0, 0, 0.20, 5, "mixed",  0, 0.00, "none"},
            {"bulged_cross_5holes", 12, 0, 0, 1, 0.20, 5, "mixed",  0, 0.00, "none"},
            {"comb_4holes",         16, 0, 0, 0, 0.20, 4, "mixed",  0, 0.00, "none"},
            {"comb_16corner_4holes",16, 0, 0, 0, 0.20, 4, "mixed",  0, 0.00, "none"},
            {"curved_comb_4holes",  16, 0, 0, 1, 0.20, 4, "mixed",  0, 0.00, "none"},
            {"comb_16corner_4holes_curved", 16, 0, 0, 1, 0.20, 4, "mixed", 0, 0.00, "none"}
        };

        for (const auto& s : table) {
            if (n == s.key) {
                Parameters p;
                p.num_corners = s.corners;
                p.convex      = s.convex;
                p.star_shaped = s.star;
                p.curved      = s.curved;
                p.curvature   = s.curv;
                p.num_holes   = s.holes;
                p.hole_shape  = s.hShape;
                p.edge_waves  = s.waves;
                p.randomness  = s.rand;
                p.scale       = scale;
                p.degree      = degree;
                if (dp.active()) {
                    p.distortion = dp;
                } else if (std::string(s.defWarp) != "none") {
                    p.distortion.type     = s.defWarp;
                    p.distortion.strength = T(0.75);
                }
                return generate(p);
            }
        }
        GISMO_ERROR("Unknown preset name: " << name);
        return nullptr;
    }

    // ------------------------------------------------------------------------
    // Non-Bilinear Multi-Patch Geometries (via Base Topology + Global Warp)
    // ------------------------------------------------------------------------

    static gsMultiPatch<T>* generateNonBilinearMultiPatch(const std::string& name, int degree = 3,
                                                          const DistortionParams& dp = DistortionParams()) {
        std::string n = name;
        std::transform(n.begin(), n.end(), n.begin(), ::tolower);
        std::replace(n.begin(), n.end(), '-', '_');

        std::unique_ptr<gsMultiPatch<T>> mp(new gsMultiPatch<T>());
        DistortionParams effDp = dp;

        if (n == "curved_2patch" || n == "two_patch_curved" || n == "2patch_curved") {
            *mp = gsNurbsCreator<T>::BSplineSquareGrid(2, 1, T(1.0), T(-1.0), T(-0.5));
            if (!effDp.active()) { effDp.type = "wave"; effDp.strength = T(0.80); }
        } else if (n == "grid_4patch" || n == "warped_grid") {
            *mp = gsNurbsCreator<T>::BSplineSquareGrid(2, 2, T(1.0), T(-1.0), T(-1.0));
            if (!effDp.active()) { effDp.type = "weird"; effDp.strength = T(0.80); }
        } else if (n == "curved_lshape_3patch" || n == "lshape_3patch") {
            mp->addPatch(gsNurbsCreator<T>::BSplineRectangle(T(-1), T(-1), T(0), T(0)));
            mp->addPatch(gsNurbsCreator<T>::BSplineRectangle(T(0),  T(-1), T(1), T(0)));
            mp->addPatch(gsNurbsCreator<T>::BSplineRectangle(T(-1), T(0),  T(0), T(1)));
            mp->computeTopology();
            if (!effDp.active()) { effDp.type = "bend"; effDp.strength = T(0.75); }
        } else if (n == "curved_ring_4patch" || n == "ring_4patch") {
            int nu = degree + 1;
            gsKnotVector<T> kv(0.0, 1.0, 0, nu, 1, degree);
            for (int q = 0; q < 4; ++q) {
                T th0 = q * T(0.5 * M_PI), th1 = (q + 1) * T(0.5 * M_PI);
                gsMatrix<T> cq(nu * nu, 2);
                for (int j = 0; j < nu; ++j) {
                    T r = T(0.5) + (T(j) / T(nu - 1)) * T(0.5);
                    for (int i = 0; i < nu; ++i) {
                        T th = th0 + (T(i) / T(nu - 1)) * (th1 - th0);
                        cq(i + j * nu, 0) = r * std::cos(th);
                        cq(i + j * nu, 1) = r * std::sin(th);
                    }
                }
                mp->addPatch(typename gsGeometry<T>::uPtr(new gsTensorBSpline<2, T>(kv, kv, cq)));
            }
            mp->computeTopology();
            if (!effDp.active()) return mp.release();
        } else {
            GISMO_ERROR("Unknown multi-patch preset: " << name);
        }

        return distortMultiPatch(*mp, effDp, degree);
    }

    // ------------------------------------------------------------------------
    // Export Suite (XML, ParaView VTP/VTU/PVD, Mathematica w/ Knots & Holes)
    // ------------------------------------------------------------------------

    static void saveXml(const gsPlanarDomain<T>& domain, const std::string& filename) {
        gsFileData<T> fd;
        fd << domain;
        fd.save(filename);
    }

    static void saveMultiPatchXml(const gsMultiPatch<T>& mp, const std::string& filename) {
        gsFileData<T> fd;
        fd << mp;
        fd.save(filename);
    }

    static TriMesh saveParaView(const gsPlanarDomain<T>& domain,
                                const std::string& base_path,
                                int num_samples_per_curve = 30) {
        std::string dir  = gsFileManager::getPath(base_path);
        std::string stem = gsFileManager::getBasename(base_path);
        if (stem.empty()) stem = base_path;
        std::string full = dir.empty() ? stem : (dir + "/" + stem);

        TriMesh mesh = triangulate(domain, num_samples_per_curve);

        // 1. Surface VTP + UnstructuredGrid VTU
        std::vector<std::vector<int>> triPolys;
        triPolys.reserve(mesh.elements.size());
        for (const auto& t : mesh.elements) triPolys.push_back({t.v0, t.v1, t.v2});
        writeVtp(full + "_surface.vtp", mesh.vertices, {}, triPolys, false);
        writeVtu(mesh, full + ".vtu");

        // 2. Boundary curves VTP
        std::vector<Point2D> crvPts;
        std::vector<std::vector<int>> crvLines;
        for (int l = 0; l < domain.numLoops(); ++l) {
            const auto& lp = domain.loop(l);
            for (int c = 0; c < lp.numCurves(); ++c) {
                auto pts = sampleCurve(lp.curve(c), num_samples_per_curve);
                std::vector<int> idx(pts.size());
                for (size_t i = 0; i < pts.size(); ++i) {
                    idx[i] = static_cast<int>(crvPts.size());
                    crvPts.push_back(pts[i]);
                }
                crvLines.push_back(std::move(idx));
            }
        }
        writeVtp(full + "_curves.vtp", crvPts, crvLines, {}, false);

        // 3. Sharp corners VTP
        auto corners = extractCorners(domain);
        writeVtp(full + "_corners.vtp", corners, {}, {}, true);

        // 4. Control net VTP
        std::vector<Point2D> cnetPts;
        std::vector<std::vector<int>> cnetLines;
        for (int l = 0; l < domain.numLoops(); ++l) {
            const auto& lp = domain.loop(l);
            for (int c = 0; c < lp.numCurves(); ++c) {
                const auto& cf = lp.curve(c).coefs();
                std::vector<int> idx(cf.rows());
                for (index_t i = 0; i < cf.rows(); ++i) {
                    idx[i] = static_cast<int>(cnetPts.size());
                    cnetPts.emplace_back(cf(i, 0), cf(i, 1));
                }
                if (idx.size() >= 2) cnetLines.push_back(std::move(idx));
            }
        }
        writeVtp(full + "_cnet.vtp", cnetPts, cnetLines, {}, false);

        // 5. PVD Collection
        std::ofstream pvd(full + ".pvd");
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n"
            << "    <DataSet timestep=\"0\" group=\"surface\" part=\"0\" file=\"" << stem << "_surface.vtp\"/>\n"
            << "    <DataSet timestep=\"0\" group=\"boundary_curves\" part=\"1\" file=\"" << stem << "_curves.vtp\"/>\n"
            << "    <DataSet timestep=\"0\" group=\"corners\" part=\"2\" file=\"" << stem << "_corners.vtp\"/>\n"
            << "    <DataSet timestep=\"0\" group=\"control_net\" part=\"3\" file=\"" << stem << "_cnet.vtp\"/>\n"
            << "  </Collection>\n</VTKFile>\n";
        return mesh;
    }

    static void saveMathematica(const gsPlanarDomain<T>& domain, const std::string& filename) {
        std::ofstream m(filename);
        m << "(* Mathematica B-Spline visualization for G+Smo Planar Domain *)\n";
        auto writeLoop = [&](const gsCurveLoop<T>& lp) {
            m << "{\n";
            for (int c = 0; c < lp.numCurves(); ++c) {
                const auto& crv = lp.curve(c);
                m << "  BSplineCurve[{";
                for (index_t i = 0; i < crv.coefs().rows(); ++i)
                    m << "{" << crv.coefs()(i, 0) << "," << crv.coefs()(i, 1) << "}"
                      << (i + 1 < crv.coefs().rows() ? "," : "");
                m << "}, SplineDegree -> " << crv.degree();
                if (const auto* bsp = dynamic_cast<const gsBSpline<T>*>(&crv)) {
                    m << ", SplineKnots -> {";
                    const auto& kv = bsp->knots();
                    for (size_t k = 0; k < kv.size(); ++k)
                        m << kv[k] << (k + 1 < kv.size() ? "," : "");
                    m << "}";
                }
                m << "]" << (c + 1 < lp.numCurves() ? ",\n" : "\n");
            }
            m << "}";
        };

        m << "outerCurves = "; writeLoop(domain.outer()); m << ";\n";
        m << "holeLoops = {";
        for (int h = 1; h <= domain.numHoles(); ++h) {
            writeLoop(domain.loop(h));
            if (h < domain.numHoles()) m << ", ";
        }
        m << "};\n";

        auto corners = extractCorners(domain);
        m << "corners = {";
        for (size_t i = 0; i < corners.size(); ++i)
            m << "{" << corners[i].x << "," << corners[i].y << "}" << (i + 1 < corners.size() ? "," : "");
        m << "};\n";
        m << "Graphics[{{Thick, Blue, outerCurves}, {Thick, Orange, holeLoops}, "
          << "{Red, PointSize[0.02], Point[corners]}}, Axes -> True, AspectRatio -> Automatic]\n";
    }

private:
    static gsCurveLoop<T>* makeEllipseLoop(Point2D ctr, T rx, T ry, int degree,
                                           bool is_outer, T rot = T(0)) {
        const T kx = rx * T(4.0 / 3.0) * (std::sqrt(T(2.0)) - T(1.0));
        const T ky = ry * T(4.0 / 3.0) * (std::sqrt(T(2.0)) - T(1.0));
        std::vector<std::vector<Point2D>> quads = {
            {{rx, 0},  {rx, ky},   {kx, ry},   {0, ry}},
            {{0, ry},  {-kx, ry},  {-rx, ky},  {-rx, 0}},
            {{-rx, 0}, {-rx, -ky}, {-kx, -ry}, {0, -ry}},
            {{0, -ry}, {kx, -ry},  {rx, -ky},  {rx, 0}}
        };
        for (auto& q : quads)
            for (auto& p : q) p = ctr + p.rotated(rot);
        return makeCurveLoop(quads, degree, is_outer);
    }

    static gsCurveLoop<T>* fitRadialClosedCurve(std::function<T(T)> r_func, int segs, int degree) {
        std::vector<std::vector<Point2D>> quads;
        quads.reserve(segs);
        const T d_th = T(2.0 * M_PI) / T(segs), eps = T(1e-4);
        for (int s = 0; s < segs; ++s) {
            T th0 = s * d_th, th1 = (s + 1) * d_th;
            T r0 = r_func(th0), r1 = r_func(th1);
            Point2D p0(r0 * std::cos(th0), r0 * std::sin(th0));
            Point2D p3(r1 * std::cos(th1), r1 * std::sin(th1));
            T dr0 = (r_func(th0 + eps) - r_func(th0 - eps)) / (T(2.0) * eps);
            T dr1 = (r_func(th1 + eps) - r_func(th1 - eps)) / (T(2.0) * eps);
            Point2D t0(dr0 * std::cos(th0) - r0 * std::sin(th0), dr0 * std::sin(th0) + r0 * std::cos(th0));
            Point2D t1(dr1 * std::cos(th1) - r1 * std::sin(th1), dr1 * std::sin(th1) + r1 * std::cos(th1));
            T chord = p0.dist(p3);
            quads.push_back({p0, p0 + t0.normalized() * (chord / T(3.0)),
                             p3 - t1.normalized() * (chord / T(3.0)), p3});
        }
        return makeCurveLoop(quads, degree, true);
    }

    static gsCurveLoop<T>* generateSmoothHorseshoe(T scale, int degree) {
        T ro = scale * T(1.1), ri = scale * T(0.45), rc = T(0.5) * (ro - ri);
        T ko = ro * T(4.0 / 3.0), ki = ri * T(4.0 / 3.0), kc = rc * T(4.0 / 3.0);
        return makeCurveLoop({{{ro, 0},  {ro, ko},   {-ro, ko},  {-ro, 0}},
                              {{-ro, 0}, {-ro, -kc}, {-ri, -kc}, {-ri, 0}},
                              {{-ri, 0}, {-ri, ki},  {ri, ki},   {ri, 0}},
                              {{ri, 0},  {ri, -kc},  {ro, -kc},  {ro, 0}}},
                             degree, true);
    }

    static void writeVtp(const std::string& filename,
                         const std::vector<Point2D>& pts,
                         const std::vector<std::vector<int>>& lines,
                         const std::vector<std::vector<int>>& polys,
                         bool as_verts) {
        std::ofstream f(filename);
        f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <PolyData>\n"
          << "    <Piece NumberOfPoints=\"" << pts.size() << "\" NumberOfVerts=\"" << (as_verts ? pts.size() : 0)
          << "\" NumberOfLines=\"" << lines.size() << "\" NumberOfPolys=\"" << polys.size() << "\">\n"
          << "      <Points><DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& p : pts) f << "        " << p.x << " " << p.y << " 0\n";
        f << "      </DataArray></Points>\n";

        auto writeCells = [&](const char* tag, const std::vector<std::vector<int>>& cells) {
            if (cells.empty()) return;
            f << "      <" << tag << ">\n        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
            for (const auto& c : cells) {
                f << "          ";
                for (int id : c) f << id << " ";
                f << "\n";
            }
            f << "        </DataArray>\n        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
            int off = 0;
            for (const auto& c : cells) { off += static_cast<int>(c.size()); f << off << " "; }
            f << "\n        </DataArray>\n      </" << tag << ">\n";
        };

        if (as_verts && !pts.empty()) {
            std::vector<std::vector<int>> v(pts.size());
            for (size_t i = 0; i < pts.size(); ++i) v[i] = {static_cast<int>(i)};
            writeCells("Verts", v);
        }
        writeCells("Lines", lines);
        writeCells("Polys", polys);
        f << "    </Piece>\n  </PolyData>\n</VTKFile>\n";
    }

    static void writeVtu(const TriMesh& mesh, const std::string& filename) {
        std::ofstream f(filename);
        f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n"
          << "    <Piece NumberOfPoints=\"" << mesh.vertices.size() << "\" NumberOfCells=\"" << mesh.elements.size() << "\">\n"
          << "      <Points><DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& v : mesh.vertices) f << "        " << v.x << " " << v.y << " 0\n";
        f << "      </DataArray></Points>\n      <Cells>\n"
          << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const auto& t : mesh.elements) f << "          " << t.v0 << " " << t.v1 << " " << t.v2 << "\n";
        f << "        </DataArray>\n        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
        for (size_t i = 1; i <= mesh.elements.size(); ++i) f << (i * 3) << " ";
        f << "\n        </DataArray>\n        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
        for (size_t i = 0; i < mesh.elements.size(); ++i) f << "5 ";
        f << "\n        </DataArray>\n      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    }
};

} // namespace gismo
