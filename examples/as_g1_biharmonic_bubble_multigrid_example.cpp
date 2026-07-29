/** @file as_g1_biharmonic_bubble_multigrid_example.cpp
 *
 *  @brief AS-G1 Biharmonic Multigrid Solver with Manufactured Tensor-Product B-Spline Bubble Solution.
 *
 *  Because the manufactured bubble solution u_bubble and its gradient vanish continuously
 *  on the domain boundary, homogeneous biharmonic boundary conditions (u = 0 and grad(u).n = 0)
 *  are satisfied EXACTLY on the domain boundary.
 *  This guarantees that the computed AS-G1 numerical solution converges asymptotically
 *  to the exact manufactured solution with optimal L2, H1, and H2 convergence rates.
 *
 *  This file is part of the G+Smo library.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Author(s): Antigravity AI, F. Hasanova, S. Takacs
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <gismo.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <chrono>
#include <ctime>
#include <vector>

#include <gsMultiGrid/gsAsG1Multigrid.h>

using namespace gismo;
using namespace gismo::expr;

std::string getBubbleTimestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
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

template <typename T>
class gsLaplacianFunction : public gsFunction<T> {
  const gsFunction<T> *m_f;

public:
  explicit gsLaplacianFunction(const gsFunction<T> &f) : m_f(&f) {}

  GISMO_CLONE_FUNCTION(gsLaplacianFunction)

  short_t domainDim() const override { return m_f->domainDim(); }
  short_t targetDim() const override { return 1; }

  void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const override {
    gsMatrix<T> d2;
    m_f->deriv2_into(u, d2);
    const short_t d = m_f->domainDim();
    result = d2.topRows(d).colwise().sum();
  }
};

template <typename T>
void computeInscribedRectangle(const gsMultiPatch<T>& mp,
                                T& rect_xmin, T& rect_xmax,
                                T& rect_ymin, T& rect_ymax)
{
    std::set<std::pair<size_t, index_t>> ifcSides;
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        ifcSides.insert({it->first().patch, it->first().side()});
        ifcSides.insert({it->second().patch, it->second().side()});
    }

    struct BdSeg { std::vector<std::pair<T,T>> pts; };
    std::vector<BdSeg> segs;
    const int K = 40;

    for (size_t pi = 0; pi < mp.nPatches(); ++pi)
    {
        for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
        {
            if (ifcSides.count({pi, side.m_index})) continue;

            BdSeg seg;
            const int si       = side.m_index;
            const short_t dir  = static_cast<short_t>((si - 1) / 2);
            const short_t tang = 1 - dir;
            const T fixedVal   = ((si - 1) % 2) ? T(1) : T(0);

            for (int k = 0; k <= K; ++k)
            {
                gsMatrix<T> uv(2, 1);
                uv(dir,  0) = fixedVal;
                uv(tang, 0) = T(k) / T(K);
                gsMatrix<T> xy;
                mp.patch(pi).eval_into(uv, xy);
                seg.pts.push_back({xy(0, 0), xy(1, 0)});
            }
            segs.push_back(seg);
        }
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
        int best = -1;
        T bestD = T(1e30);
        bool rev = false;
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
        if (rev)
            for (int k = (int)pts.size()-2; k >= 0; --k) poly.push_back(pts[k]);
        else
            for (size_t k = 1; k < pts.size(); ++k) poly.push_back(pts[k]);
    }

    const T xmin = bbox(0, 0), xmax = bbox(0, 1);
    const T ymin = bbox(1, 0), ymax = bbox(1, 1);
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
                    index_t dl = i - il, dr = ir - i;
                    index_t db = j - jb, dt = jt - j;
                    index_t m = std::min({dl, dr, db, dt});
                    if      (m == dl) il++;
                    else if (m == dr) ir--;
                    else if (m == db) jb++;
                    else              jt--;
                }
    }

    if (il + 2 < ir) { il++; ir--; }
    if (jb + 2 < jt) { jb++; jt--; }

    rect_xmin = xmin + il * dx; rect_xmax = xmin + ir * dx;
    rect_ymin = ymin + jb * dy; rect_ymax = ymin + jt * dy;
}

template<typename T>
index_t solvePCG(const gsSparseMatrix<T>& A,
                 const gsAsG1Multigrid<T>& mg,
                 const gsMatrix<T>& b,
                 gsMatrix<T>& x,
                 index_t maxIter = 100,
                 T tol = 1e-8,
                 index_t finestLevel = 1)
{
    const index_t n = A.rows();
    x = gsMatrix<T>::Zero(n, 1);

    gsMatrix<T> r = b - A * x;
    T r0_norm = r.norm();
    if (r0_norm < 1e-15) return 0;

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

        T r_norm = r.norm();
        T relRes = r_norm / r0_norm;

        if (relRes < tol)
        {
            return iter;
        }

        gsMatrix<T> z_new = gsMatrix<T>::Zero(n, 1);
        mg.vCycle(finestLevel, z_new, r, 2, 2);

        T rz_new = (r.transpose() * z_new)(0, 0);
        T beta = rz_new / rz;

        p = z_new + beta * p;
        rz = rz_new;
    }

    return maxIter;
}

int main(int argc, char *argv[]) {
    using T = real_t;

    std::string geometry("domain2d/2patch/multipatch/weirdo_multivalence_non_bilinear.xml");
    std::string outDir("");
    index_t degree = 3;
    index_t minRefinements = 2;
    index_t maxRefinements = 4;
    real_t margin = 0.15;

    gsCmdLine cmd("AS-G1 Multigrid Biharmonic Solver with Manufactured Bubble Solution.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
    cmd.addInt("d", "degree", "Target discretization degree (minimum 3).", degree);
    cmd.addInt("m", "minRefinements", "Minimum uniform refinement level (coarsest level, min 2).", minRefinements);
    cmd.addInt("r", "maxRefinements", "Maximum uniform refinement level (finest level).", maxRefinements);
    cmd.addReal("", "margin", "Fractional margin between bubble rectangle and bounding box (0..0.5).", margin);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    if (outDir.empty()) {
        outDir = "experiments/as_g1_bubble_multigrid_" + getBubbleTimestamp();
    }
    if (outDir.back() != '/') outDir += '/';
    gsFileManager::mkdir(outDir);

    if (degree < 3) degree = 3;
    if (margin <= 0 || margin >= 0.5) margin = 0.15;

    gsInfo << "\n=========================================================================================================\n";
    gsInfo << "  AS-G1 Geometric Multigrid Biharmonic Solver with Manufactured B-Spline Bubble Solution\n";
    gsInfo << "  Geometry File:   " << geometry << "\n";
    gsInfo << "  Refinement:      r_min = " << minRefinements << " to r_max = " << maxRefinements << "\n";
    gsInfo << "  Output Dir:      " << outDir << "\n";
    gsInfo << "=========================================================================================================\n\n";

    // 1. Compute inscribed rectangle for bubble support (geometry-aware)
    T inscribed_xmin, inscribed_xmax, inscribed_ymin, inscribed_ymax;
    {
        gsMultiPatch<T>::uPtr tmpMp = gsReadFile<>(geometry);
        if (!tmpMp) {
            gsInfo << "Error: Cannot read geometry " << geometry << "\n";
            return -1;
        }
        tmpMp->computeTopology();
        computeInscribedRectangle(*tmpMp, inscribed_xmin, inscribed_xmax, inscribed_ymin, inscribed_ymax);
    }

    // Table Header
    gsInfo << std::setw(5)  << "r"
           << std::setw(12) << "h_max"
           << std::setw(10) << "N_free"
           << std::setw(14) << "L2 Error" << std::setw(9) << "L2 Rate"
           << std::setw(14) << "H1 Error" << std::setw(9) << "H1 Rate"
           << std::setw(14) << "H2 Error" << std::setw(9) << "H2 Rate"
           << std::setw(10) << "PCG Iters"
           << std::setw(14) << "Solve Time"
           << "\n";
    gsInfo << std::string(118, '-') << "\n" << std::flush;

    T prev_h = 0;
    T prev_l2 = 0;
    T prev_h1 = 0;
    T prev_h2 = 0;

    for (index_t maxRef = minRefinements; maxRef <= maxRefinements; ++maxRef) {
        try {
            gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
            if (!mpPtr) {
                gsInfo << "Error: Cannot read geometry " << geometry << "\n";
                return -1;
            }

            gsMultiPatch<T> mpFine = *mpPtr;
            mpFine.computeTopology();
            if (mpFine.patch(0).basis().degree(0) < degree) {
                mpFine.degreeElevate(degree - mpFine.patch(0).basis().degree(0));
            }
            const short_t deg = mpFine.patch(0).basis().degree(0);
            const index_t mult = std::max<index_t>(deg - 1, 1);
            for (index_t i = 0; i < maxRef; ++i) mpFine.uniformRefine(1, mult);

            const T h_max = std::pow(0.5, maxRef);

            // Bounding box of geometry
            gsMatrix<T> bbox;
            mpFine.boundingBox(bbox);
            T x_min = bbox(0, 0), x_max = bbox(0, 1);
            T y_min = bbox(1, 0), y_max = bbox(1, 1);
            T Lx = (x_max - x_min > 1e-12) ? (x_max - x_min) : T(1);
            T Ly = (y_max - y_min > 1e-12) ? (y_max - y_min) : T(1);

            // Construct manufactured tensor-product B-spline bubble
            const index_t p_bubble = deg;
            const T insc_Lx = inscribed_xmax - inscribed_xmin;
            const T insc_Ly = inscribed_ymax - inscribed_ymin;
            const T a_x = inscribed_xmin + margin * insc_Lx;
            const T b_x = inscribed_xmax - margin * insc_Lx;
            const T a_y = inscribed_ymin + margin * insc_Ly;
            const T b_y = inscribed_ymax - margin * insc_Ly;

            const T pad_x = 0.05 * Lx;
            const T pad_y = 0.05 * Ly;

            auto makeBubbleKV = [](T lo, T hi, T a, T b, index_t p) {
                std::vector<T> kn;
                kn.reserve(3 * p + 4);
                for (index_t i = 0; i <= p; ++i) kn.push_back(lo);
                for (index_t i = 0; i <= p + 1; ++i) kn.push_back(a * T(p + 1 - i) / T(p + 1) + b * T(i) / T(p + 1));
                for (index_t i = 0; i <= p; ++i) kn.push_back(hi);
                return gsKnotVector<T>(static_cast<short_t>(p), kn.begin(), kn.end());
            };

            gsKnotVector<T> kv_x = makeBubbleKV(x_min - pad_x, x_max + pad_x, a_x, b_x, p_bubble);
            gsKnotVector<T> kv_y = makeBubbleKV(y_min - pad_y, y_max + pad_y, a_y, b_y, p_bubble);
            gsTensorBSplineBasis<2, T> bubbleBasis(kv_x, kv_y);

            const index_t size0 = kv_x.size() - p_bubble - 1;
            const index_t bumpIdx = (p_bubble + 1) + (p_bubble + 1) * size0;

            gsMatrix<T> bubbleCoefs(bubbleBasis.size(), 1);
            bubbleCoefs.setZero();
            bubbleCoefs(bumpIdx, 0) = 1.0;
            gsTensorBSpline<2, T> bubble(bubbleBasis, give(bubbleCoefs));

            gsPiecewiseFunction<T> bubbleField;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) bubbleField.addPiece(bubble);

            gsLaplacianFunction<T> bubbleLap(bubble);
            gsPiecewiseFunction<T> bubbleLapField;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) bubbleLapField.addPiece(bubbleLap);

            // Construct AS-G1 Multigrid Hierarchy
            auto t0_setup = std::chrono::high_resolution_clock::now();
            gsAsG1Multigrid<T> mgSolver(*mpPtr, degree, minRefinements, maxRef);
            auto t1_setup = std::chrono::high_resolution_clock::now();

            index_t finestLevel = mgSolver.numLevels() - 1;
            const gsSparseMatrix<T>& K_fine = mgSolver.stiffnessMatrix(finestLevel);
            index_t nFree = mgSolver.freeDofs(finestLevel);

            // Assemble Disjoint Load Vector F_disjoint = (Delta B_a, Delta u_bubble)
            gsMultiBasis<T> dbasis(mpFine);
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());

            auto G_map = A.getMap(mpFine);
            auto u_space = A.getSpace(dbasis);
            auto lap_ex = A.getCoeff(bubbleLapField, G_map);

            A.initSystem();
            A.assemble(ilapl(u_space, G_map) * lap_ex * meas(G_map));
            const gsMatrix<T>& F_disjoint = A.rhs();

            // Project RHS to free AS-G1 DOFs: b_free = S_fine^T * F_disjoint
            gsMatrix<T> gd = computeGluingData(mpFine, T(1e-8), 0);
            std::vector<gsArgyrisEmbedding<T>> argBasis;
            gsVector<index_t> patchDofSizes(mpFine.nPatches());
            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                argBasis.push_back(deriveArgyrisBasisEmbedding(
                    dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mpFine.patch(i).basis()),
                    gsMatrix<T>(gd.row(i)), mpFine.patch(i)));
                patchDofSizes[i] = argBasis[i].matrix.cols();
            }

            gsConstantFunction<T> zero;
            gsBoundaryConditions<T> bc;
            for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
                bc.add(it->patch, it->side(), "ValuesAndDerivatives", zero);
            gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis, bc);

            index_t nDisjoint = 0;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) nDisjoint += argBasis[i].matrix.rows();

            gsSparseEntries<T> tFreeEntries;
            index_t rowOffset = 0;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                const gsSparseMatrix<T> &Ai = argBasis[i].matrix;
                const index_t patchSize = mapper.patchSize(i);
                for (index_t j = 0; j < patchSize; ++j) {
                    if (mapper.is_boundary(j, i)) continue;
                    index_t gIdx = mapper.index(j, i);
                    for (typename gsSparseMatrix<T>::InnerIterator it(Ai, j); it; ++it) {
                        tFreeEntries.add(rowOffset + it.row(), gIdx, it.value());
                    }
                }
                rowOffset += Ai.rows();
            }

            gsSparseMatrix<T> S_fine(nDisjoint, nFree);
            S_fine.setFrom(tFreeEntries);

            gsMatrix<T> b_free = S_fine.transpose() * F_disjoint;

            // Solve via Multigrid PCG
            gsMatrix<T> u_pcg(nFree, 1);
            auto t0_solve = std::chrono::high_resolution_clock::now();
            index_t pcgIters = solvePCG(K_fine, mgSolver, b_free, u_pcg, 100, 1e-8, finestLevel);
            auto t1_solve = std::chrono::high_resolution_clock::now();
            double solveTime_ms = std::chrono::duration<double, std::milli>(t1_solve - t0_solve).count();

            // Reconstruct multi-patch field
            gsMatrix<T> c_disjoint = S_fine * u_pcg;
            gsMultiPatch<T> solField;
            index_t cOffset = 0;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                const index_t sz = argBasis[i].matrix.rows();
                gsMatrix<T> ci = c_disjoint.block(cOffset, 0, sz, 1);
                cOffset += sz;
                const gsTensorBSplineBasis<2, T> &tb = dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mpFine.patch(i).basis());
                solField.addPatch(tb.makeGeometry(give(ci)));
            }

            // Error Evaluation against manufactured u_bubble
            gsExprEvaluator<T> ev;
            ev.setIntegrationDomain(dbasis.domain());
            auto G_map_ev = ev.getMap(mpFine);

            auto u_exact_ev = ev.getVariable(bubbleField, G_map_ev);
            auto u_sol_ev = ev.getVariable(solField);

            const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
            const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev)).sqNorm() * meas(G_map_ev)));
            const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev)).sqNorm() * meas(G_map_ev)));

            // Convergence Rates
            T l2_rate = 0, h1_rate = 0, h2_rate = 0;
            if (maxRef > minRefinements && prev_h > 0) {
                const T h_ratio = prev_h / h_max;
                l2_rate = std::log(prev_l2 / l2err) / std::log(h_ratio);
                h1_rate = std::log(prev_h1 / h1err) / std::log(h_ratio);
                h2_rate = std::log(prev_h2 / h2err) / std::log(h_ratio);
            }

            // Print Row
            gsInfo << std::setw(5)  << maxRef
                   << std::setw(12) << std::scientific << std::setprecision(4) << h_max
                   << std::setw(10) << nFree
                   << std::setw(14) << std::scientific << std::setprecision(4) << l2err;
            if (maxRef > minRefinements)
                gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << l2_rate;
            else
                gsInfo << std::setw(9) << "-";

            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
            if (maxRef > minRefinements)
                gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h1_rate;
            else
                gsInfo << std::setw(9) << "-";

            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
            if (maxRef > minRefinements)
                gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h2_rate;
            else
                gsInfo << std::setw(9) << "-";

            gsInfo << std::setw(10) << pcgIters
                   << std::setw(14) << std::fixed << std::setprecision(2) << solveTime_ms << " ms\n" << std::flush;

            prev_h = h_max;
            prev_l2 = l2err;
            prev_h1 = h1err;
            prev_h2 = h2err;

            // Export VTK plots to timestamped output folder
            std::string solName = outDir + "as_g1_bubble_multigrid_sol_r" + std::to_string(maxRef);
            std::string exactName = outDir + "as_g1_bubble_multigrid_exact_r" + std::to_string(maxRef);

            gsField<T> solF(mpFine, solField);
            gsWriteParaview<>(solF, solName, 1000);

            gsField<T> exactF(mpFine, bubbleField, false);
            gsWriteParaview<>(exactF, exactName, 1000);

        } catch (const std::exception &e) {
            gsInfo << "ERROR on level " << maxRef << ": " << e.what() << "\n" << std::flush;
        }
    }

    gsInfo << std::string(118, '-') << "\n";
    gsInfo << "VTK plot files exported to directory: " << outDir << "\n";
    gsInfo << "=========================================================================================================\n";

    return 0;
}
