#include <gismo.h>
#include <gsModeling/gsPlanarDomainGenerator.h>
#include <gsModeling/gsPatchSplitter.h>
#include <gsModeling/gsPatchSplitter.hpp>
#include <gsIO/gsFileData.h>
#include <gsIO/gsFileManager.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace gismo;

using Pt2 = gsPoint2D<real_t>;

static std::vector<Pt2> straightControlPoints(const Pt2& a, const Pt2& b, int degree = 3)
{
    std::vector<Pt2> cps;
    cps.reserve(degree + 1);
    for (int i = 0; i <= degree; ++i)
    {
        real_t t = static_cast<real_t>(i) / static_cast<real_t>(degree);
        cps.push_back(a + (b - a) * t);
    }
    return cps;
}

static gsPlanarDomain<real_t>* makePolygonDomain(const std::vector<Pt2>& pts, int degree = 3)
{
    std::vector<std::vector<Pt2>> curves;
    int n = static_cast<int>(pts.size());
    curves.reserve(n);
    for (int i = 0; i < n; ++i)
        curves.push_back(straightControlPoints(pts[i], pts[(i + 1) % n], degree));

    auto* loop = gsPlanarDomainGenerator<real_t>::makeCurveLoop(curves, degree, true);
    return new gsPlanarDomain<real_t>(loop);
}

// 1. The "Dart" / Boomerang (Non-Convex Quadrilateral, K=4)
static gsPlanarDomain<real_t>* makeDartDomain(int degree = 3)
{
    std::vector<Pt2> pts = {
        {0.0, 0.0},
        {2.0, 1.0},
        {1.0, 0.2},
        {2.0, -1.0}
    };
    return makePolygonDomain(pts, degree);
}

// 2. The Serpentine S-Duct (K=4 sharp corners, 2 straight inlet/outlet, 2 cubic S-curves)
// Interpolates (1,0)->(2,1.5)->(0,3)->(1,4.5) on right wall and (0,4.5)->(-1,3)->(1,1.5)->(0,0) on left wall
// at t = 0, 1/3, 2/3, 1 so the channel bends past x=1 and x=0 with zero valid straight diagonals.
static gsPlanarDomain<real_t>* makeSerpentineSDuctDomain(int degree = 3)
{
    std::vector<std::vector<Pt2>> curves;
    // Inlet: (0,0) -> (1,0)
    curves.push_back(straightControlPoints(Pt2(0.0, 0.0), Pt2(1.0, 0.0), degree));
    // Right wall cubic Bezier interpolating (1,0), (2,1.5), (0,3), (1,4.5)
    curves.push_back({
        Pt2(1.0, 0.0),
        Pt2(5.5, 1.5),
        Pt2(-3.5, 3.0),
        Pt2(1.0, 4.5)
    });
    // Outlet: (1,4.5) -> (0,4.5)
    curves.push_back(straightControlPoints(Pt2(1.0, 4.5), Pt2(0.0, 4.5), degree));
    // Left wall cubic Bezier interpolating (0,4.5), (-1,3), (1,1.5), (0,0)
    curves.push_back({
        Pt2(0.0, 4.5),
        Pt2(-4.5, 3.0),
        Pt2(4.5, 1.5),
        Pt2(0.0, 0.0)
    });

    auto* loop = gsPlanarDomainGenerator<real_t>::makeCurveLoop(curves, degree, true);
    return new gsPlanarDomain<real_t>(loop);
}

// 3. The Filleted U-Bracket (4 sharp G0 corners at open ends, smooth G1 inner & outer fillets)
static gsPlanarDomain<real_t>* makeFilletedUBracketDomain(int degree = 3)
{
    const real_t x_out = 1.2;
    const real_t x_in  = 0.5;
    const real_t y_top = 1.4;
    const real_t y_arc = 0.2;
    const real_t r_in  = 0.5;
    const real_t y_bot = -1.0;
    const real_t r_out = 0.4;
    const real_t kappa = 4.0 / 3.0 * (std::sqrt(2.0) - 1.0);
    const real_t k_in  = r_in * kappa;
    const real_t k_out = r_out * kappa;

    std::vector<std::vector<Pt2>> curves = {
        // 0. Left outer vertical wall
        straightControlPoints(Pt2(-x_out, y_top), Pt2(-x_out, y_bot + r_out), degree),
        // 1. Outer bottom-left G1 fillet
        {
            Pt2(-x_out, y_bot + r_out),
            Pt2(-x_out, y_bot + r_out - k_out),
            Pt2(-x_out + r_out - k_out, y_bot),
            Pt2(-x_out + r_out, y_bot)
        },
        // 2. Outer bottom horizontal edge
        straightControlPoints(Pt2(-x_out + r_out, y_bot), Pt2(x_out - r_out, y_bot), degree),
        // 3. Outer bottom-right G1 fillet
        {
            Pt2(x_out - r_out, y_bot),
            Pt2(x_out - r_out + k_out, y_bot),
            Pt2(x_out, y_bot + r_out - k_out),
            Pt2(x_out, y_bot + r_out)
        },
        // 4. Right outer vertical wall
        straightControlPoints(Pt2(x_out, y_bot + r_out), Pt2(x_out, y_top), degree),
        // 5. Right prong top (sharp G0)
        straightControlPoints(Pt2(x_out, y_top), Pt2(x_in, y_top), degree),
        // 6. Right inner vertical wall
        straightControlPoints(Pt2(x_in, y_top), Pt2(x_in, y_arc), degree),
        // 7. Inner right G1 semicircular half
        {
            Pt2(x_in, y_arc),
            Pt2(x_in, y_arc - k_in),
            Pt2(k_in, y_arc - r_in),
            Pt2(0.0, y_arc - r_in)
        },
        // 8. Inner left G1 semicircular half
        {
            Pt2(0.0, y_arc - r_in),
            Pt2(-k_in, y_arc - r_in),
            Pt2(-x_in, y_arc - k_in),
            Pt2(-x_in, y_arc)
        },
        // 9. Left inner vertical wall
        straightControlPoints(Pt2(-x_in, y_arc), Pt2(-x_in, y_top), degree),
        // 10. Left prong top (sharp G0)
        straightControlPoints(Pt2(-x_in, y_top), Pt2(-x_out, y_top), degree)
    };

    auto* loop = gsPlanarDomainGenerator<real_t>::makeCurveLoop(curves, degree, true);
    return new gsPlanarDomain<real_t>(loop);
}

// 4. The Heptagon Trap (K=7 regular heptagon -> CentroidStarFan creates valence nu=7)
static gsPlanarDomain<real_t>* makeHeptagonDomain(int degree = 3)
{
    std::vector<Pt2> pts;
    pts.reserve(7);
    for (int i = 0; i < 7; ++i)
    {
        real_t theta = 2.0 * EIGEN_PI * static_cast<real_t>(i) / 7.0;
        pts.emplace_back(std::cos(theta), std::sin(theta));
    }
    return makePolygonDomain(pts, degree);
}

// 5. The "Hourglass Pinch" (K=4 with deep inward quadratic arcs; waist width = 0.2 at x=1)
static gsPlanarDomain<real_t>* makeHourglassDomain(int degree = 3)
{
    // Degree-elevate exact quadratic Bezier curves:
    // Left quadratic:  (0,1) -> (1.8,0) -> (0,-1) passes through (0.9, 0) at t=0.5
    // Right quadratic: (2,-1) -> (0.2,0) -> (2,1) passes through (1.1, 0) at t=0.5
    std::vector<std::vector<Pt2>> curves = {
        // Bottom edge: (0,-1) -> (2,-1)
        straightControlPoints(Pt2(0.0, -1.0), Pt2(2.0, -1.0), degree),
        // Right inward quadratic arc (degree-elevated to cubic): (2,-1) -> (1.1,0) -> (2,1)
        {
            Pt2(2.0, -1.0),
            Pt2(0.8, -1.0 / 3.0),
            Pt2(0.8,  1.0 / 3.0),
            Pt2(2.0,  1.0)
        },
        // Top edge: (2,1) -> (0,1)
        straightControlPoints(Pt2(2.0, 1.0), Pt2(0.0, 1.0), degree),
        // Left inward quadratic arc (degree-elevated to cubic): (0,1) -> (0.9,0) -> (0,-1)
        {
            Pt2(0.0,  1.0),
            Pt2(1.2,  1.0 / 3.0),
            Pt2(1.2, -1.0 / 3.0),
            Pt2(0.0, -1.0)
        }
    };

    auto* loop = gsPlanarDomainGenerator<real_t>::makeCurveLoop(curves, degree, true);
    return new gsPlanarDomain<real_t>(loop);
}

// 6. The 3-Way "Y-Junction" Branching Pipe (K=6 convex corners, 3 smooth reflex saddle curves)
static gsPlanarDomain<real_t>* makeYJunctionDomain(int degree = 3)
{
    const real_t L = 2.0;
    const real_t w = 0.45;

    auto rotPoint = [](const Pt2& p, real_t angle) -> Pt2 {
        real_t c = std::cos(angle), s = std::sin(angle);
        return Pt2(p.x * c - p.y * s, p.x * s + p.y * c);
    };

    // Base saddle control points for sector 0 (from arm 0 top corner to arm 1 bottom corner):
    // B0 = (L, w) -> P1 = (0.60, 0.08) -> P2 (reflection of P1 across 60 deg) -> A1 (rot 120 of (L,-w))
    Pt2 B0(L, w);
    Pt2 P1(0.60, 0.08);
    // Reflect (x,y) across angle 60 deg: (x*cos(120)+y*sin(120), x*sin(120)-y*cos(120))
    const real_t c120 = -0.5;
    const real_t s120 = std::sqrt(3.0) * 0.5;
    Pt2 P2(P1.x * c120 + P1.y * s120, P1.x * s120 - P1.y * c120);
    Pt2 A1 = rotPoint(Pt2(L, -w), 2.0 * EIGEN_PI / 3.0);

    std::vector<std::vector<Pt2>> curves;
    for (int k = 0; k < 3; ++k)
    {
        real_t th = 2.0 * EIGEN_PI * static_cast<real_t>(k) / 3.0;
        Pt2 Ak = rotPoint(Pt2(L, -w), th);
        Pt2 Bk = rotPoint(Pt2(L, w), th);
        // Opening k
        curves.push_back(straightControlPoints(Ak, Bk, degree));
        // Reflex saddle k -> k+1
        curves.push_back({
            rotPoint(B0, th),
            rotPoint(P1, th),
            rotPoint(P2, th),
            rotPoint(A1, th)
        });
    }

    auto* loop = gsPlanarDomainGenerator<real_t>::makeCurveLoop(curves, degree, true);
    return new gsPlanarDomain<real_t>(loop);
}

// 7. The Asymmetric Staircase / "Tetris Step" (K=8, Non-Star, Non-Comb)
static gsPlanarDomain<real_t>* makeStaircaseDomain(int degree = 3)
{
    std::vector<Pt2> pts = {
        {0.0, 0.0},
        {3.0, 0.0},
        {3.0, 1.0},
        {2.0, 1.0},
        {2.0, 2.0},
        {1.0, 2.0},
        {1.0, 3.0},
        {0.0, 3.0}
    };
    return makePolygonDomain(pts, degree);
}

static bool segProperIntersect(const Pt2& a, const Pt2& b, const Pt2& c, const Pt2& d)
{
    auto cross = [](const Pt2& p, const Pt2& q, const Pt2& r) {
        return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
    };
    real_t d1 = cross(a, b, c), d2 = cross(a, b, d);
    real_t d3 = cross(c, d, a), d4 = cross(c, d, b);
    const real_t eps = 1e-9;
    return (((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
            ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps)));
}

static bool validateDomainGeometry(const gsPlanarDomain<real_t>& dom,
                                   int& outSharp, int& outSmooth, int& outReflex,
                                   real_t& outArea, real_t& outMaxGap)
{
    const auto& outer = dom.outer();
    int K = outer.size();
    outSharp = 0;
    outSmooth = 0;
    outReflex = 0;
    outMaxGap = 0.0;

    for (int i = 0; i < K; ++i)
    {
        const auto& cPrev = outer.curve((i + K - 1) % K);
        const auto& cNext = outer.curve(i);
        gsMatrix<real_t> u1(1, 1), u0(1, 1), p1, p0, d1, d0;
        u1(0, 0) = cPrev.parameterRange()(0, 1);
        u0(0, 0) = cNext.parameterRange()(0, 0);
        cPrev.eval_into(u1, p1);
        cNext.eval_into(u0, p0);
        real_t gap = std::hypot(p1(0, 0) - p0(0, 0), p1(1, 0) - p0(1, 0));
        outMaxGap = std::max(outMaxGap, gap);

        cPrev.deriv_into(u1, d1);
        cNext.deriv_into(u0, d0);
        real_t n1 = std::hypot(d1(0, 0), d1(1, 0));
        real_t n2 = std::hypot(d0(0, 0), d0(1, 0));
        if (n1 > 1e-14 && n2 > 1e-14)
        {
            real_t cr = (d1(0, 0) * d0(1, 0) - d1(1, 0) * d0(0, 0)) / (n1 * n2);
            real_t dt = (d1(0, 0) * d0(0, 0) + d1(1, 0) * d0(1, 0)) / (n1 * n2);
            real_t ang = std::atan2(cr, dt);
            if (std::abs(ang) > 0.10)
            {
                ++outSharp;
                if (ang < -0.10) ++outReflex;
            }
            else
            {
                ++outSmooth;
            }
        }
    }

    auto poly = gsPlanarDomainGenerator<real_t>::sampleLoop(outer, 40);
    outArea = gsPlanarDomainGenerator<real_t>::computeSignedArea(outer, 40);

    bool selfInt = false;
    int m = static_cast<int>(poly.size());
    for (int i = 0; i < m && !selfInt; ++i)
    {
        for (int j = i + 2; j < m; ++j)
        {
            if (i == 0 && j == m - 1) continue;
            if (segProperIntersect(poly[i], poly[(i + 1) % m], poly[j], poly[(j + 1) % m]))
            {
                selfInt = true;
                break;
            }
        }
    }

    return (dom.numHoles() == 0) && (outMaxGap < 1e-10) && (outArea > 0.0) && !selfInt;
}

int main()
{
    const std::string outDir = "filedata/domain2d/noHolesCounterExamples";
    gsFileManager::mkdir(outDir);

    struct Spec {
        std::string name;
        gsPlanarDomain<real_t>* (*factory)(int);
    };

    std::vector<Spec> specs = {
        {"01_dart_4corner_nonconvex",        makeDartDomain},
        {"02_serpentine_s_duct",             makeSerpentineSDuctDomain},
        {"03_filleted_u_bracket",            makeFilletedUBracketDomain},
        {"04_heptagon_7corner_convex",       makeHeptagonDomain},
        {"05_hourglass_pinch",               makeHourglassDomain},
        {"06_y_junction_3way",               makeYJunctionDomain},
        {"07_asymmetric_staircase_8corner",  makeStaircaseDomain}
    };

    std::cout << std::left
              << std::setw(34) << "Geometry"
              << " | K  | Shp/Smt/Ref | Area     | Gap      | Valid | Splitter Probe\n";
    std::cout << std::string(110, '-') << "\n";

    int validCount = 0;
    for (const auto& sp : specs)
    {
        std::unique_ptr<gsPlanarDomain<real_t>> dom(sp.factory(3));
        int nSharp = 0, nSmooth = 0, nReflex = 0;
        real_t area = 0.0, maxGap = 0.0;
        bool valid = validateDomainGeometry(*dom, nSharp, nSmooth, nReflex, area, maxGap);
        if (valid) ++validCount;

        const std::string prefix = outDir + "/" + sp.name;
        gsPlanarDomainGenerator<real_t>::saveXml(*dom, prefix + ".xml");
        gsPlanarDomainGenerator<real_t>::saveParaView(*dom, prefix, 40);
        gsPlanarDomainGenerator<real_t>::saveMathematica(*dom, prefix + ".m");

        // Probe current gsPatchSplitter to observe counter-example behavior
        gsPatchSplitter<real_t> splitter(*dom);
        splitter.setStrategy(gsPatchSplitter<real_t>::Strategy::AutomaticBest);
        splitter.setDegree(3);
        bool splitOk = splitter.computeSplitting();
        std::string splitInfo = splitOk
            ? (splitter.bestCandidate().methodUsed + " (P=" +
               std::to_string(splitter.getMultiPatch().nPatches()) +
               ", minDet=" + std::to_string(splitter.bestCandidate().minJacobianDet) + ")")
            : "FAILED";

        std::cout << std::setw(34) << sp.name
                  << " | " << std::setw(2) << dom->outer().size()
                  << " | " << nSharp << "/" << nSmooth << "/" << nReflex << "       "
                  << " | " << std::fixed << std::setprecision(4) << std::setw(8) << area
                  << " | " << std::scientific << std::setprecision(1) << std::setw(8) << maxGap
                  << " | " << (valid ? "PASS " : "FAIL ")
                  << " | " << splitInfo << "\n";
    }

    std::cout << std::string(110, '-') << "\n";
    std::cout << "Generated " << validCount << " / " << specs.size()
              << " valid domains in " << outDir << "\n";
    return (validCount == static_cast<int>(specs.size())) ? 0 : 1;
}
