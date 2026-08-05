/** @file as_g1_biharmonic_bubble_multigrid_example.cpp
 *
 *  @brief AS-G1 biharmonic multigrid solver with a manufactured tensor-product
 *         B-spline "bubble" solution.
 *
 *  Because the bubble and its gradient vanish continuously on the domain
 *  boundary, the homogeneous biharmonic boundary conditions (u = 0 and
 *  grad(u).n = 0) hold EXACTLY, so the AS-G1 numerical solution converges with
 *  optimal L2/H1/H2 rates. The linear system is solved with the generic
 *  gsMultiGridOp-preconditioned CG provided by gsAsG1BiharmonicMG
 *  (gsMultiGrid/gsAsG1BiharmonicMG.h), reusing the shared AS-G1 DOF coupling
 *  and gluing-data machinery instead of a hand-rolled V-cycle.
 *
 *  This file is part of the G+Smo library.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Author(s): F. Hasanova, S. Takacs
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gismo.h>

#include <gsMultiGrid/gsAsG1BiharmonicMG.h>

using namespace gismo;
using namespace gismo::expr;

// ISO-like timestamp string YYYYMMDD_HHMMSS for the default output folder.
static std::string getTimestampString()
{
    std::time_t now_c = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());
    struct tm *parts = std::localtime(&now_c);
    std::ostringstream oss;
    oss << (1900 + parts->tm_year)
        << std::setfill('0') << std::setw(2) << (1 + parts->tm_mon)
        << std::setfill('0') << std::setw(2) << parts->tm_mday << "_"
        << std::setfill('0') << std::setw(2) << parts->tm_hour
        << std::setfill('0') << std::setw(2) << parts->tm_min
        << std::setfill('0') << std::setw(2) << parts->tm_sec;
    return oss.str();
}

// Physical Laplacian of a wrapped scalar function (value only), used to build
// the biharmonic energy-form right-hand side l(v) = (Delta u_ex, Delta v).
template <typename T>
class gsLaplacianFunction : public gsFunction<T>
{
    const gsFunction<T>* m_f;
public:
    explicit gsLaplacianFunction(const gsFunction<T>& f) : m_f(&f) {}
    GISMO_CLONE_FUNCTION(gsLaplacianFunction)
    short_t domainDim() const override { return m_f->domainDim(); }
    short_t targetDim() const override { return 1; }
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> d2;
        m_f->deriv2_into(u, d2);
        const short_t d = m_f->domainDim();
        result = d2.topRows(d).colwise().sum();
    }
};

// Largest axis-aligned rectangle inscribed in the multi-patch domain. The
// external boundary sides are obtained directly from the box topology
// (mp.bBegin()/bEnd()); they are sampled into a polygon, and a rectangle is
// grown inside it with a grid + ray-casting point-in-polygon test.
template <typename T>
void computeInscribedRectangle(const gsMultiPatch<T>& mp,
                               T& rect_xmin, T& rect_xmax,
                               T& rect_ymin, T& rect_ymax)
{
    struct BdSeg { std::vector<std::pair<T,T>> pts; };
    std::vector<BdSeg> segs;
    const int K = 40;
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        BdSeg seg;
        const int si = it->side().m_index;
        const short_t dir  = static_cast<short_t>((si - 1) / 2);
        const short_t tang = 1 - dir;
        const T fixedVal = ((si - 1) % 2) ? T(1) : T(0);
        for (int k = 0; k <= K; ++k)
        {
            gsMatrix<T> uv(2, 1);
            uv(dir, 0) = fixedVal;
            uv(tang, 0) = T(k) / T(K);
            gsMatrix<T> xy;
            mp.patch(it->patch).eval_into(uv, xy);
            seg.pts.push_back({xy(0, 0), xy(1, 0)});
        }
        segs.push_back(seg);
    }

    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    if (segs.empty())
    {
        rect_xmin = bbox(0, 0); rect_xmax = bbox(0, 1);
        rect_ymin = bbox(1, 0); rect_ymax = bbox(1, 1);
        return;
    }

    auto dsq = [](std::pair<T,T> a, std::pair<T,T> b) -> T {
        return (a.first-b.first)*(a.first-b.first) + (a.second-b.second)*(a.second-b.second);
    };

    std::vector<std::pair<T,T>> poly;
    std::vector<bool> used(segs.size(), false);
    used[0] = true;
    for (auto& p : segs[0].pts) poly.push_back(p);
    for (size_t step = 1; step < segs.size(); ++step)
    {
        auto last = poly.back();
        int best = -1; T bestD = T(1e30); bool rev = false;
        for (size_t k = 0; k < segs.size(); ++k)
        {
            if (used[k]) continue;
            T d0 = dsq(last, segs[k].pts.front());
            T d1 = dsq(last, segs[k].pts.back());
            if (d0 < bestD) { bestD = d0; best = (int)k; rev = false; }
            if (d1 < bestD) { bestD = d1; best = (int)k; rev = true;  }
        }
        if (best < 0) break;
        used[best] = true;
        auto& pts = segs[best].pts;
        if (rev) for (int k = (int)pts.size()-2; k >= 0; --k) poly.push_back(pts[k]);
        else     for (size_t k = 1; k < pts.size(); ++k)      poly.push_back(pts[k]);
    }

    const T xmin = bbox(0,0), xmax = bbox(0,1), ymin = bbox(1,0), ymax = bbox(1,1);
    const T Lx = xmax - xmin, Ly = ymax - ymin;
    const index_t N = 80;
    const T dx = Lx / N, dy = Ly / N;
    auto inPoly = [&](T px, T py) -> bool {
        bool in = false;
        for (size_t i = 0, j = poly.size()-1; i < poly.size(); j = i++) {
            T yi = poly[i].second, yj = poly[j].second;
            T xi = poly[i].first,  xj = poly[j].first;
            if (((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi))
                in = !in;
        }
        return in;
    };
    std::vector<std::vector<bool>> grid(N+1, std::vector<bool>(N+1, false));
    for (index_t i = 0; i <= N; ++i)
        for (index_t j = 0; j <= N; ++j)
            grid[i][j] = inPoly(xmin + i*dx, ymin + j*dy);

    index_t il = 0, ir = N, jb = 0, jt = N;
    bool chg = true;
    while (chg && il < ir && jb < jt) {
        chg = false;
        for (index_t j = jb; j <= jt; ++j) if (!grid[il][j]) { il++; chg = true; break; }
        for (index_t j = jb; j <= jt; ++j) if (!grid[ir][j]) { ir--; chg = true; break; }
        for (index_t i = il; i <= ir; ++i) if (!grid[i][jb]) { jb++; chg = true; break; }
        for (index_t i = il; i <= ir; ++i) if (!grid[i][jt]) { jt--; chg = true; break; }
    }
    bool ok = false;
    while (!ok && il + 1 < ir && jb + 1 < jt) {
        ok = true;
        for (index_t i = il; i <= ir && ok; ++i)
            for (index_t j = jb; j <= jt && ok; ++j)
                if (!grid[i][j]) {
                    ok = false;
                    index_t dl = i - il, dr = ir - i, db = j - jb, dt = jt - j;
                    index_t m = std::min({dl, dr, db, dt});
                    if      (m == dl) il++;
                    else if (m == dr) ir--;
                    else if (m == db) jb++;
                    else              jt--;
                }
    }
    if (il + 2 < ir) { il++; ir--; }
    if (jb + 2 < jt) { jb++; jt--; }
    rect_xmin = xmin + il*dx; rect_xmax = xmin + ir*dx;
    rect_ymin = ymin + jb*dy; rect_ymax = ymin + jt*dy;
}

// Clamped degree-p knot vector on [lo,hi] containing p+2 simple bump knots.
template <typename T>
gsKnotVector<T> makeBubbleKV(T lo, T hi, T a, T b, index_t p)
{
    std::vector<T> kn;
    kn.reserve(3 * p + 4);
    for (index_t i = 0; i <= p; ++i) kn.push_back(lo);
    for (index_t i = 0; i <= p + 1; ++i)
        kn.push_back(a * T(p + 1 - i) / T(p + 1) + b * T(i) / T(p + 1));
    for (index_t i = 0; i <= p; ++i) kn.push_back(hi);
    return gsKnotVector<T>(static_cast<short_t>(p), kn.begin(), kn.end());
}

int main(int argc, char *argv[])
{
    using T = real_t;

    std::string geometry("domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
    std::string outDir("");
    index_t degree = 3;
    index_t minRefinements = 2;
    index_t maxRefinements = 4;
    index_t numSmooth = 3;
    real_t margin = 0.15;
    bool wcycle = false;

    gsCmdLine cmd("AS-G1 Multigrid Biharmonic Solver with Manufactured Bubble Solution.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
    cmd.addInt("d", "degree", "Target discretization degree (minimum 3).", degree);
    cmd.addInt("m", "minRefinements", "Minimum uniform refinement level (coarsest level, min 2).", minRefinements);
    cmd.addInt("r", "maxRefinements", "Maximum uniform refinement level (finest level).", maxRefinements);
    cmd.addInt("s", "smooth", "Number of pre-/post-smoothing steps.", numSmooth);
    cmd.addReal("", "margin", "Fractional margin between bubble rectangle and bounding box (0..0.5).", margin);
    cmd.addSwitch("wcycle", "Use W-cycles (numCycles = 2) instead of V-cycles.", wcycle);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    if (degree < 3) degree = 3;
    if (minRefinements < 2) minRefinements = 2;
    if (maxRefinements < minRefinements) maxRefinements = minRefinements;
    if (numSmooth < 1) numSmooth = 1;
    if (margin <= 0 || margin >= 0.5) margin = 0.15;

    if (outDir.empty())
        outDir = "experiments/as_g1_bubble_multigrid_" + getTimestampString();
    if (outDir.back() != '/') outDir += '/';
    gsFileManager::mkdir(outDir);

    gsInfo << "\n=========================================================================================================\n";
    gsInfo << "  AS-G1 Geometric Multigrid Biharmonic Solver with Manufactured B-Spline Bubble Solution\n";
    gsInfo << "  Solver:          gsMultiGridOp (Galerkin, symmetric Gauss-Seidel, "
           << numSmooth << " smoothing steps) + CG, " << (wcycle ? "W-cycle" : "V-cycle") << "\n";
    gsInfo << "  Geometry File:   " << geometry << "\n";
    gsInfo << "  Refinement:      r_min = " << minRefinements << " to r_max = " << maxRefinements << "\n";
    gsInfo << "  Output Dir:      " << outDir << "\n";
    gsInfo << "=========================================================================================================\n\n";

    gsMultiPatch<T>::uPtr basePtr = gsReadFile<>(geometry);
    if (!basePtr) {
        gsInfo << "Error: Cannot read geometry " << geometry << "\n";
        return -1;
    }
    gsMultiPatch<T> base = *basePtr;
    base.computeTopology();

    // Bubble support: inscribed rectangle (geometry-aware) and bounding box.
    gsMatrix<T> bbox;
    base.boundingBox(bbox);
    const T x_min = bbox(0,0), x_max = bbox(0,1);
    const T y_min = bbox(1,0), y_max = bbox(1,1);
    const T Lx = (x_max - x_min > 1e-12) ? (x_max - x_min) : T(1);
    const T Ly = (y_max - y_min > 1e-12) ? (y_max - y_min) : T(1);

    T ix0, ix1, iy0, iy1;
    computeInscribedRectangle(base, ix0, ix1, iy0, iy1);
    const T iLx = ix1 - ix0, iLy = iy1 - iy0;
    const T a_x = ix0 + margin * iLx, b_x = ix1 - margin * iLx;
    const T a_y = iy0 + margin * iLy, b_y = iy1 - margin * iLy;
    const T pad_x = 0.05 * Lx, pad_y = 0.05 * Ly;
    const index_t pB = degree;

    // Tensor-product B-spline bubble (single central bump embedded in a clamped
    // basis spanning the padded bounding box, so it is safe to evaluate anywhere).
    gsKnotVector<T> kvx = makeBubbleKV(x_min - pad_x, x_max + pad_x, a_x, b_x, pB);
    gsKnotVector<T> kvy = makeBubbleKV(y_min - pad_y, y_max + pad_y, a_y, b_y, pB);
    gsTensorBSplineBasis<2, T> bubbleBasis(kvx, kvy);
    const index_t size0 = kvx.size() - pB - 1;
    const index_t bumpIdx = (pB + 1) + (pB + 1) * size0;
    gsMatrix<T> bubbleCoefs(bubbleBasis.size(), 1);
    bubbleCoefs.setZero();
    bubbleCoefs(bumpIdx, 0) = 1.0;
    gsTensorBSpline<2, T> bubble(bubbleBasis, give(bubbleCoefs));
    gsLaplacianFunction<T> bubbleLap(bubble);

    // Table header.
    gsInfo << std::setw(5)  << "r"
           << std::setw(12) << "h_max"
           << std::setw(10) << "N_free"
           << std::setw(14) << "L2 Error" << std::setw(9) << "L2 Rate"
           << std::setw(14) << "H1 Error" << std::setw(9) << "H1 Rate"
           << std::setw(14) << "H2 Error" << std::setw(9) << "H2 Rate"
           << std::setw(10) << "MG Iters"
           << std::setw(14) << "Solve Time"
           << "\n";
    gsInfo << std::string(118, '-') << "\n" << std::flush;

    T prev_h = 0, prev_l2 = 0, prev_h1 = 0, prev_h2 = 0;

    for (index_t maxRef = minRefinements; maxRef <= maxRefinements; ++maxRef) {
        try {
            gsAsG1BiharmonicMG<T> solver(base, degree, minRefinements, maxRef, numSmooth);
            solver.setNumSmooth(numSmooth, numSmooth);
            if (wcycle) solver.setNumCycles(2);

            const gsMultiPatch<T> &mp = solver.fineMultiPatch();
            gsMultiBasis<T> dbasis(mp);
            const T h_max = std::pow(0.5, maxRef);

            gsPiecewiseFunction<T> bubbleField, bubbleLapField;
            for (size_t i = 0; i < mp.nPatches(); ++i) {
                bubbleField.addPiece(bubble);
                bubbleLapField.addPiece(bubbleLap);
            }

            // Energy-form RHS l(v) = (Delta u_ex, Delta v); homogeneous clamped BC.
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());
            auto G_map = A.getMap(mp);
            auto u_space = A.getSpace(dbasis);
            auto lap_ex = A.getCoeff(bubbleLapField, G_map);
            A.initSystem();
            A.assemble(ilapl(u_space, G_map) * lap_ex * meas(G_map));
            gsMatrix<T> b_free = solver.transformFine().transpose() * A.rhs();

            // Solve via multigrid-preconditioned CG.
            gsMatrix<T> c_free;
            auto t0_solve = std::chrono::high_resolution_clock::now();
            index_t mgIters = solver.solve(b_free, c_free, 200, T(1e-8));
            auto t1_solve = std::chrono::high_resolution_clock::now();
            double solveTime_ms = std::chrono::duration<double, std::milli>(t1_solve - t0_solve).count();

            // Reconstruct multi-patch field and evaluate error against the bubble.
            gsMultiPatch<T> solField = solver.reconstruct(c_free);
            gsExprEvaluator<T> ev;
            ev.setIntegrationDomain(dbasis.domain());
            auto G_map_ev = ev.getMap(mp);
            auto u_exact_ev = ev.getVariable(bubbleField, G_map_ev);
            auto u_sol_ev = ev.getVariable(solField);

            const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
            const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev)).sqNorm() * meas(G_map_ev)));
            const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev)).sqNorm() * meas(G_map_ev)));

            T l2_rate = 0, h1_rate = 0, h2_rate = 0;
            if (maxRef > minRefinements && prev_h > 0) {
                const T q = std::log(prev_h / h_max);
                l2_rate = std::log(prev_l2 / l2err) / q;
                h1_rate = std::log(prev_h1 / h1err) / q;
                h2_rate = std::log(prev_h2 / h2err) / q;
            }

            gsInfo << std::setw(5)  << maxRef
                   << std::setw(12) << std::scientific << std::setprecision(4) << h_max
                   << std::setw(10) << solver.freeDofs()
                   << std::setw(14) << std::scientific << std::setprecision(4) << l2err;
            if (maxRef > minRefinements) gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << l2_rate; else gsInfo << std::setw(9) << "-";
            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
            if (maxRef > minRefinements) gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h1_rate; else gsInfo << std::setw(9) << "-";
            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
            if (maxRef > minRefinements) gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h2_rate; else gsInfo << std::setw(9) << "-";
            gsInfo << std::setw(10) << mgIters
                   << std::setw(14) << std::fixed << std::setprecision(2) << solveTime_ms << " ms\n" << std::flush;

            prev_h = h_max; prev_l2 = l2err; prev_h1 = h1err; prev_h2 = h2err;

            // Export VTK plots to the output folder.
            gsField<T> solF(mp, solField);
            gsWriteParaview<>(solF, outDir + "as_g1_bubble_multigrid_sol_r" + std::to_string(maxRef), 1000);
            gsField<T> exactF(mp, bubbleField, false);
            gsWriteParaview<>(exactF, outDir + "as_g1_bubble_multigrid_exact_r" + std::to_string(maxRef), 1000);

        } catch (const std::exception &e) {
            gsInfo << "ERROR on level " << maxRef << ": " << e.what() << "\n" << std::flush;
        }
    }

    gsInfo << std::string(118, '-') << "\n";
    gsInfo << "VTK plot files exported to directory: " << outDir << "\n";
    gsInfo << "=========================================================================================================\n";

    return 0;
}
