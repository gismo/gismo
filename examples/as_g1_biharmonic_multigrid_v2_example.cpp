/** @file as_g1_biharmonic_multigrid_v2_example.cpp
 *
 *  @brief Analysis-Suitable G1 (AS-G1) biharmonic multigrid solver, version 2.
 *
 *  This example rewrites as_g1_biharmonic_multigrid_example.cpp. The earlier
 *  version (i) used a manufactured solution sin(a*pi*x)*cos(a*pi*y) that does
 *  NOT satisfy the homogeneous clamped biharmonic boundary conditions it
 *  imposes -- so the discrete error never converged -- and (ii) relied on a
 *  hand-rolled V-cycle whose Jacobi/Schwarz smoother did not cope with the
 *  ill-conditioned fourth-order operator, so the PCG iteration stagnated at the
 *  iteration cap.
 *
 *  The rewrite:
 *    - Uses a tensor-product B-spline "bubble" manufactured solution that
 *      vanishes together with its derivatives on the domain boundary, so the
 *      homogeneous biharmonic BCs hold exactly and the discrete error converges
 *      at the optimal rates.
 *    - Builds the AS-G1 DOF coupling with makeMapperForArgyrisBasis()
 *      (gsModeling/gsAsG1Basis.hpp) and the gluing data with computeGluingData()
 *      (gsModeling/gsAsG1Domain.hpp), instead of duplicating the interface
 *      matching inline.
 *    - Solves with the generic gsMultiGridOp (gsMultiGrid/gsMultiGrid.h) as a
 *      preconditioner for conjugate gradients, with Galerkin coarsening,
 *      symmetric Gauss-Seidel smoothing and a configurable (larger) number of
 *      smoothing steps.
 *
 *  For reference it also runs the original hand-rolled gsAsG1Multigrid solver on
 *  the same operator, so the iteration counts and timings can be compared.
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
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <gismo.h>

#include <gsMultiGrid/gsAsG1BiharmonicMG.h>
#include <gsMultiGrid/gsAsG1Multigrid.h> // original implementation (for comparison)

using namespace gismo;
using namespace gismo::expr;

// ---------------------------------------------------------------------------
// Physical Laplacian of a wrapped scalar function (value only), used to build
// the biharmonic energy-form right-hand side l(v) = (Delta u_ex, Delta v).
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Largest axis-aligned rectangle inscribed in the multi-patch domain.
// (Same geometry-aware algorithm as as_g1_biharmonic_bubble_example.cpp.)
// ---------------------------------------------------------------------------
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

// Original hand-rolled preconditioned CG (verbatim from the old example) for
// the reference comparison.
template<typename T>
index_t solvePCG_orig(const gsSparseMatrix<T>& A, const gsAsG1Multigrid<T>& mg,
                      const gsMatrix<T>& b, gsMatrix<T>& x,
                      index_t maxIter, T tol, index_t finestLevel)
{
    const index_t n = A.rows();
    x = gsMatrix<T>::Zero(n, 1);
    gsMatrix<T> r = b - A * x;
    T r0 = r.norm();
    if (r0 < 1e-15) return 0;
    gsMatrix<T> z = gsMatrix<T>::Zero(n, 1);
    mg.vCycle(finestLevel, z, r, 2, 2);
    gsMatrix<T> p = z;
    T rz = (r.transpose() * z)(0, 0);
    for (index_t iter = 1; iter <= maxIter; ++iter)
    {
        gsMatrix<T> Ap = A * p;
        T pAp = (p.transpose() * Ap)(0, 0);
        if (std::abs(pAp) < 1e-18) return iter;
        T alpha = rz / pAp;
        x += alpha * p;
        r -= alpha * Ap;
        if (r.norm() / r0 < tol) return iter;
        gsMatrix<T> zn = gsMatrix<T>::Zero(n, 1);
        mg.vCycle(finestLevel, zn, r, 2, 2);
        T rzn = (r.transpose() * zn)(0, 0);
        p = zn + (rzn / rz) * p;
        rz = rzn;
    }
    return maxIter;
}

int main(int argc, char* argv[])
{
    using T = real_t;

    std::string dir("filedata/domain2d/2patch/multipatch/");
    std::string file("");
    index_t degree = 3;
    index_t minRef = 2;
    index_t maxRef = 4;
    index_t numSmooth = 3;
    real_t margin = 0.15;
    bool noCompare = false;
    std::string outDir("");
    bool wcycle = false;

    gsCmdLine cmd("AS-G1 biharmonic multigrid solver (v2): "
                  "shared DOF mapper + generic gsMultiGridOp preconditioned CG.");
    cmd.addString("d", "dir", "Directory of multipatch geometry XML files.", dir);
    cmd.addString("f", "file", "Solve a single geometry XML file instead of the whole directory.", file);
    cmd.addInt("p", "degree", "Discretization degree (minimum 3).", degree);
    cmd.addInt("m", "minRef", "Coarsest refinement level (>= 2).", minRef);
    cmd.addInt("r", "maxRef", "Finest refinement level.", maxRef);
    cmd.addInt("s", "smooth", "Number of pre-/post-smoothing steps.", numSmooth);
    cmd.addReal("g", "margin", "Bubble rectangle margin (0..0.5).", margin);
    cmd.addString("o", "outDir", "Output folder for PVD/VTK plots (empty = no plotting).", outDir);
    cmd.addSwitch("wcycle", "Use W-cycles (numCycles = 2) instead of V-cycles.", wcycle);
    cmd.addSwitch("noCompare", "Skip the reference run of the original solver.", noCompare);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (degree < 3)  degree = 3;
    if (minRef < 2)  minRef = 2;
    if (maxRef < minRef) maxRef = minRef;
    if (numSmooth < 1) numSmooth = 1;
    if (margin <= 0 || margin >= 0.5) margin = 0.15;
    if (dir.back() != '/') dir += '/';
    if (!outDir.empty())
    {
        if (outDir.back() != '/') outDir += '/';
        gsFileManager::mkdir(outDir);
    }

    // Geometry list.
    std::vector<std::string> geoms;
    if (!file.empty())
        geoms.push_back(file);
    else
    {
        const char* names[] = {
            "2_patch_rectangle.xml",
            "2_patch_rectangle_non_bilinear.xml",
            "two_patch_mathematica.xml",
            "two_patch_mathematica_non_bilinear.xml",
            "three_patch_belt.xml",
            "three_patch_belt_non_bilinear.xml",
            "4_patch_3_valence_not_bilinear.xml",
            "5_patch_5_valence_not_bilinear.xml",
            "square_6_patch.xml",
            "square_6_patch_non_bilinear.xml",
            "6_patch_4_valence_each_not_bilinear.xml",
            "6_patch_6_valence_not_bilinear.xml",
            "experiment_multivalence_not_bilinear.xml",
            "weirdo_multivalence.xml",
            "weirdo_multivalence_non_bilinear.xml"
        };
        for (const char* n : names)
        {
            std::string full = dir + n;
            if (gsFileManager::fileExists(full)) geoms.push_back(full);
        }
    }

    gsInfo << "\n=====================================================================================================\n";
    gsInfo << "  AS-G1 Biharmonic Multigrid Solver (v2)\n";
    gsInfo << "  DOF coupling : makeMapperForArgyrisBasis + computeGluingData\n";
    gsInfo << "  Solver       : gsMultiGridOp (Galerkin, symmetric Gauss-Seidel, "
           << numSmooth << " smoothing steps) + CG\n";
    gsInfo << "  Cycle type   : " << (wcycle ? "W-cycle" : "V-cycle") << "\n";
    gsInfo << "  Levels       : r = " << minRef << " .. " << maxRef
           << "   Degree = " << degree << "\n";
    gsInfo << "  Geometries   : " << geoms.size() << "\n";
    if (!outDir.empty()) gsInfo << "  Plot folder  : " << outDir << "\n";
    gsInfo << "=====================================================================================================\n";

    index_t nConverged = 0, nTotal = 0;
    struct Summary { std::string name; index_t iters; T relres; T h2; bool ok; };
    std::vector<Summary> summaries;

    for (const std::string& geom : geoms)
    {
        std::string shortName = geom;
        size_t slash = shortName.find_last_of('/');
        if (slash != std::string::npos) shortName = shortName.substr(slash + 1);

        gsInfo << "\n-----------------------------------------------------------------------------------------------------\n";
        gsInfo << "Geometry: " << shortName << "\n";

        gsMultiPatch<T>::uPtr basePtr = gsReadFile<>(geom);
        if (!basePtr)
        {
            gsInfo << "  Cannot read geometry, skipping.\n";
            continue;
        }
        gsMultiPatch<T> base = *basePtr;
        base.computeTopology();

        // Bubble manufactured solution (fixed across refinement levels).
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

        gsPiecewiseFunction<T> bubbleField, bubbleLapField;
        for (size_t i = 0; i < base.nPatches(); ++i)
        {
            bubbleField.addPiece(bubble);
            bubbleLapField.addPiece(bubbleLap);
        }

        // Convergence table (NEW solver).
        gsInfo << std::setw(5)  << "r"
               << std::setw(10) << "N_free"
               << std::setw(13) << "L2 Error"  << std::setw(8) << "rate"
               << std::setw(13) << "H1 Error"  << std::setw(8) << "rate"
               << std::setw(13) << "H2 Error"  << std::setw(8) << "rate"
               << std::setw(8)  << "MG-it"
               << std::setw(12) << "rel.res."
               << std::setw(12) << "solve[ms]"
               << "\n";
        gsInfo << std::string(102, '-') << "\n" << std::flush;

        T ph = 0, pl2 = 0, ph1 = 0, ph2 = 0;
        Summary lastSum{shortName, 0, 0, 0, false};

        for (index_t r = minRef; r <= maxRef; ++r)
        {
            try
            {
                gsAsG1BiharmonicMG<T> solver(base, degree, minRef, r, numSmooth);
                solver.setNumSmooth(numSmooth, numSmooth);
                if (wcycle) solver.setNumCycles(2);

                const gsMultiPatch<T>& mp = solver.fineMultiPatch();
                gsMultiBasis<T> dbasis(mp);

                // Energy-form RHS: (Delta u_ex, Delta v).
                gsExprAssembler<T> A(1, 1);
                A.setIntegrationDomain(dbasis.domain());
                auto G = A.getMap(mp);
                auto u = A.getSpace(dbasis);
                auto lap_ex = A.getCoeff(bubbleLapField, G);
                A.initSystem();
                A.assemble(ilapl(u, G) * lap_ex * meas(G));
                gsMatrix<T> F_free = solver.transformFine().transpose() * A.rhs();

                // Solve with multigrid-preconditioned CG.
                gsMatrix<T> c_free;
                T relres = 0;
                auto t0 = std::chrono::high_resolution_clock::now();
                index_t iters = solver.solve(F_free, c_free, 200, T(1e-8), &relres);
                auto t1 = std::chrono::high_resolution_clock::now();
                double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

                // Error evaluation.
                gsMultiPatch<T> sol = solver.reconstruct(c_free);
                gsExprEvaluator<T> ev;
                ev.setIntegrationDomain(dbasis.domain());
                auto Gev = ev.getMap(mp);
                auto uex = ev.getVariable(bubbleField, Gev);
                auto ush = ev.getVariable(sol);
                const T l2 = std::sqrt(ev.integral((ush - uex).sqNorm() * meas(Gev)));
                const T h1 = std::sqrt(ev.integral((igrad(ush, Gev) - igrad(uex)).sqNorm() * meas(Gev)));
                const T h2 = std::sqrt(ev.integral((ihess(ush, Gev) - ihess(uex)).sqNorm() * meas(Gev)));

                const T h = std::pow(0.5, r);
                T rl2 = 0, rh1 = 0, rh2 = 0;
                if (r > minRef && ph > 0)
                {
                    const T q = std::log(ph / h);
                    rl2 = std::log(pl2 / l2) / q;
                    rh1 = std::log(ph1 / h1) / q;
                    rh2 = std::log(ph2 / h2) / q;
                }

                gsInfo << std::setw(5) << r
                       << std::setw(10) << solver.freeDofs()
                       << std::setw(13) << std::scientific << std::setprecision(3) << l2;
                if (r > minRef) gsInfo << std::setw(8) << std::fixed << std::setprecision(2) << rl2; else gsInfo << std::setw(8) << "-";
                gsInfo << std::setw(13) << std::scientific << std::setprecision(3) << h1;
                if (r > minRef) gsInfo << std::setw(8) << std::fixed << std::setprecision(2) << rh1; else gsInfo << std::setw(8) << "-";
                gsInfo << std::setw(13) << std::scientific << std::setprecision(3) << h2;
                if (r > minRef) gsInfo << std::setw(8) << std::fixed << std::setprecision(2) << rh2; else gsInfo << std::setw(8) << "-";
                gsInfo << std::setw(8) << iters
                       << std::setw(12) << std::scientific << std::setprecision(2) << relres
                       << std::setw(12) << std::fixed << std::setprecision(1) << ms
                       << "\n" << std::flush;

                ph = h; pl2 = l2; ph1 = h1; ph2 = h2;
                lastSum = Summary{shortName, iters, relres, h2, (relres < 1e-8)};

                // Optional PVD/VTK export of numerical and exact solutions.
                if (!outDir.empty())
                {
                    std::string stem = shortName;
                    size_t dot = stem.find_last_of('.');
                    if (dot != std::string::npos) stem = stem.substr(0, dot);
                    const std::string base = outDir + "v2_" + stem + "_r" + std::to_string(r);
                    gsField<T> solPlot(mp, sol);
                    gsWriteParaview<>(solPlot, base + "_sol", 1000);
                    gsField<T> exPlot(mp, bubbleField, false);
                    gsWriteParaview<>(exPlot, base + "_exact", 1000);
                }
            }
            catch (const std::exception& e)
            {
                gsInfo << std::setw(5) << r << "   ERROR: " << e.what() << "\n" << std::flush;
            }
        }

        ++nTotal;
        if (lastSum.ok) ++nConverged;
        summaries.push_back(lastSum);

        // Reference comparison against the original hand-rolled solver on the
        // same finest operator (RHS = K * ones, tol 1e-8).
        if (!noCompare)
        {
            try
            {
                gsAsG1BiharmonicMG<T> newSolver(base, degree, minRef, maxRef, numSmooth);
                newSolver.setNumSmooth(numSmooth, numSmooth);
                if (wcycle) newSolver.setNumCycles(2);
                const gsSparseMatrix<T>& Knew = newSolver.stiffnessMatrix();
                gsMatrix<T> xref = gsMatrix<T>::Ones(Knew.rows(), 1);
                gsMatrix<T> bnew = Knew * xref, xnew;
                T relNew = 0;
                auto tn0 = std::chrono::high_resolution_clock::now();
                index_t itNew = newSolver.solve(bnew, xnew, 200, T(1e-8), &relNew);
                auto tn1 = std::chrono::high_resolution_clock::now();
                double msNew = std::chrono::duration<double, std::milli>(tn1 - tn0).count();

                gsAsG1Multigrid<T> oldSolver(base, degree, minRef, maxRef);
                index_t oldFinest = oldSolver.numLevels() - 1;
                const gsSparseMatrix<T>& Kold = oldSolver.stiffnessMatrix(oldFinest);
                gsMatrix<T> xrefO = gsMatrix<T>::Ones(Kold.rows(), 1);
                gsMatrix<T> bold = Kold * xrefO, xold;
                auto to0 = std::chrono::high_resolution_clock::now();
                index_t itOld = solvePCG_orig(Kold, oldSolver, bold, xold, 200, T(1e-8), oldFinest);
                auto to1 = std::chrono::high_resolution_clock::now();
                double msOld = std::chrono::duration<double, std::milli>(to1 - to0).count();
                T relOld = (bold - Kold * xold).norm() / bold.norm();

                gsInfo << "  [solver comparison @ r=" << maxRef << ", tol 1e-8, rhs = K*1]\n";
                gsInfo << "     original (gsAsG1Multigrid V-cycle PCG): "
                       << std::setw(4) << itOld << " it, "
                       << std::fixed << std::setprecision(1) << msOld << " ms, rel.res "
                       << std::scientific << std::setprecision(2) << relOld
                       << (relOld < 1e-8 ? "" : "  (NOT converged)") << "\n";
                gsInfo << "     new      (gsMultiGridOp CG, " << numSmooth << " smooth): "
                       << std::setw(4) << itNew << " it, "
                       << std::fixed << std::setprecision(1) << msNew << " ms, rel.res "
                       << std::scientific << std::setprecision(2) << relNew
                       << (relNew < 1e-8 ? "" : "  (NOT converged)") << "\n" << std::flush;
            }
            catch (const std::exception& e)
            {
                gsInfo << "  [solver comparison skipped: " << e.what() << "]\n" << std::flush;
            }
        }
    }

    gsInfo << "\n=====================================================================================================\n";
    gsInfo << "  SUMMARY (finest level r = " << maxRef << ")\n";
    gsInfo << std::string(101, '-') << "\n";
    gsInfo << std::setw(45) << "geometry" << std::setw(10) << "MG-it"
           << std::setw(14) << "rel.res." << std::setw(14) << "H2 error"
           << std::setw(14) << "converged" << "\n";
    gsInfo << std::string(101, '-') << "\n";
    for (const Summary& s : summaries)
        gsInfo << std::setw(45) << s.name << std::setw(10) << s.iters
               << std::setw(14) << std::scientific << std::setprecision(2) << s.relres
               << std::setw(14) << std::scientific << std::setprecision(2) << s.h2
               << std::setw(14) << (s.ok ? "yes" : "NO") << "\n";
    gsInfo << std::string(101, '-') << "\n";
    gsInfo << "  Converged: " << nConverged << " / " << nTotal << " geometries.\n";
    gsInfo << "=====================================================================================================\n";

    return (nConverged == nTotal) ? 0 : 1;
}
