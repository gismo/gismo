/** @file as_g1_biharmonic_bubble_example.cpp

    @brief Biharmonic solver with a manufactured tensor-product B-spline
           "bubble" solution over Analysis-Suitable G1 (AS-G1) multi-patch
           geometries.

    The manufactured solution is a single tensor-product B-spline basis
    function of degree p supported on a rectangle that is automatically
    fitted inside the domain (the bounding box shrunk by a margin). Because
    the bubble and all its derivatives vanish outside that rectangle, the
    homogeneous biharmonic boundary conditions (u = 0 and grad(u).n = 0) are
    satisfied exactly on the domain boundary.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <gismo.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <vector>

#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

using namespace gismo;
using namespace gismo::expr;

// Scalar function returning the (physical) Laplacian of a wrapped scalar
// function. Used to build the biharmonic energy-form right-hand side
//     l(v) = (Delta u_ex, Delta v)
// with only the *value* of Delta u_ex, which is what the expression
// assembler can consume from a coefficient function.
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
    // For a scalar function of d variables, deriv2 stores the pure second
    // derivatives (f_xx, f_yy, ...) in the first d rows.
    const short_t d = m_f->domainDim();
    result = d2.topRows(d).colwise().sum();
  }
};

// Compute the largest axis-aligned rectangle inscribed in the multi-patch
// domain Omega.  The algorithm:
//   1. Samples external boundary curves to build a polygon approximating dOmega.
//   2. Uses a grid + ray-casting point-in-polygon test.
//   3. Greedily shrinks from the bounding box until the rectangle is fully
//      contained within the domain.
template <typename T>
void computeInscribedRectangle(const gsMultiPatch<T>& mp,
                                T& rect_xmin, T& rect_xmax,
                                T& rect_ymin, T& rect_ymax)
{

    // ---- 2. Sample external boundary curves ----
    struct BdSeg { std::vector<std::pair<T,T>> pts; };
    std::vector<BdSeg> segs;
    const int K = 40; // samples per boundary side

    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        const boxSide side = it->side();
        const index_t pi = it->patch;

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

    // Fallback: use bounding box if no external boundary found
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    if (segs.empty())
    {
        rect_xmin = bbox(0, 0); rect_xmax = bbox(0, 1);
        rect_ymin = bbox(1, 0); rect_ymax = bbox(1, 1);
        return;
    }

    // ---- 3. Order segments into a closed polygon ----
    auto dsq = [](std::pair<T,T> a, std::pair<T,T> b) -> T {
        return (a.first-b.first)*(a.first-b.first) +
               (a.second-b.second)*(a.second-b.second);
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
            for (int k = (int)pts.size()-2; k >= 0; --k)
                poly.push_back(pts[k]);
        else
            for (size_t k = 1; k < pts.size(); ++k)
                poly.push_back(pts[k]);
    }

    // ---- 4. Build containment grid (ray-casting even-odd rule) ----
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
            if (((yi > py) != (yj > py)) &&
                (px < (xj - xi) * (py - yi) / (yj - yi) + xi))
                in = !in;
        }
        return in;
    };

    std::vector<std::vector<bool>> grid(N+1, std::vector<bool>(N+1, false));
    for (index_t i = 0; i <= N; ++i)
        for (index_t j = 0; j <= N; ++j)
            grid[i][j] = inPoly(xmin + i*dx, ymin + j*dy);

    // ---- 5. Find inscribed rectangle by greedy shrinking ----
    index_t il = 0, ir = N, jb = 0, jt = N;

    // Shrink boundary rows/columns that contain outside cells
    bool chg = true;
    while (chg && il < ir && jb < jt) {
        chg = false;
        for (index_t j = jb; j <= jt; ++j)
            if (!grid[il][j]) { il++; chg = true; break; }
        for (index_t j = jb; j <= jt; ++j)
            if (!grid[ir][j]) { ir--; chg = true; break; }
        for (index_t i = il; i <= ir; ++i)
            if (!grid[i][jb]) { jb++; chg = true; break; }
        for (index_t i = il; i <= ir; ++i)
            if (!grid[i][jt]) { jt--; chg = true; break; }
    }

    // Verify entire interior and shrink further if holes remain
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

    // Safety: shrink by one grid cell to account for polygon approximation
    if (il + 2 < ir) { il++; ir--; }
    if (jb + 2 < jt) { jb++; jt--; }

    rect_xmin = xmin + il * dx;
    rect_xmax = xmin + ir * dx;
    rect_ymin = ymin + jb * dy;
    rect_ymax = ymin + jt * dy;
}

int main(int argc, char *argv[]) {
  using T = real_t;

  std::string geometry("domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
  std::string outDir("");
  index_t degree = 3;
  index_t numGaussPerSpan = 0;
  index_t maxRefinements = 4;
  index_t bubbleDegree = 0; // 0 = use the discretization degree
  real_t margin = 0.15;     // fractional margin between rectangle and bounding box
  bool plot = false;

  gsCmdLine cmd("AS-G1 Biharmonic Equation Solver with a Manufactured Tensor-Product B-spline Bubble Solution.");
  cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
  cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
  cmd.addInt("d", "degree", "Target discretization degree (minimum 3).", degree);
  cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).", numGaussPerSpan);
  cmd.addInt("r", "refinements", "Maximum uniform refinement levels to test.", maxRefinements);
  cmd.addInt("b", "bubbleDegree", "Degree p of the tensor-product B-spline bubble (0 = use discretization degree).", bubbleDegree);
  cmd.addReal("m", "margin", "Fractional margin between the bubble rectangle and the domain bounding box (0..0.5).", margin);
  cmd.addSwitch("plot", "Export exact manufactured solution and numerical solution to VTK files.", plot);

  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  std::string prefix;
  if (!outDir.empty()) {
    if (outDir.back() != '/')
      outDir += '/';
    gsFileManager::mkdir(outDir);
    prefix = outDir;
    gsInfo << "Output directory: " << outDir << "\n";
  }

  if (degree < 3)
    degree = 3;

  if (margin <= 0 || margin >= 0.5)
    margin = 0.15;

  // ----- Compute inscribed rectangle (geometry-aware, done once) -----
  T inscribed_xmin, inscribed_xmax, inscribed_ymin, inscribed_ymax;
  {
    gsMultiPatch<T>::uPtr tmpMp = gsReadFile<>(geometry);
    if (!tmpMp) {
      gsInfo << "Error: Cannot read geometry " << geometry << "\n";
      return -1;
    }
    tmpMp->computeTopology();
    gsInfo << "Computing domain-inscribed rectangle for bubble support...\n" << std::flush;
    computeInscribedRectangle(*tmpMp, inscribed_xmin, inscribed_xmax,
                               inscribed_ymin, inscribed_ymax);
    T rLx = inscribed_xmax - inscribed_xmin;
    T rLy = inscribed_ymax - inscribed_ymin;
    gsInfo << "  Inscribed rect (before margin): ["
           << inscribed_xmin << ", " << inscribed_xmax
           << "] x [" << inscribed_ymin << ", " << inscribed_ymax
           << "]  (size " << rLx << " x " << rLy << ")\n" << std::flush;
  }

  gsInfo << "\n======================================================================\n";
  gsInfo << "AS-G1 Biharmonic Solver (Delta^2 u = f) with Manufactured Bubble Solution\n";
  gsInfo << "Manufactured Solution: Tensor-Product B-spline Bubble on an interior rectangle\n";
  gsInfo << "Boundary Condition: u = 0 and grad(u).n = 0 (naturally zero on boundary)\n";
  gsInfo << "======================================================================\n\n";

  // Table header
  gsInfo << std::setw(5) << "r"
         << std::setw(12) << "h_max"
         << std::setw(10) << "N_free"
         << std::setw(14) << "L2 Error"
         << std::setw(10) << "L2 Rate"
         << std::setw(14) << "H1 Error"
         << std::setw(10) << "H1 Rate"
         << std::setw(14) << "H2 Error"
         << std::setw(10) << "H2 Rate"
         << "\n";
  gsInfo << std::string(99, '-') << "\n" << std::flush;

  T prev_h = 0;
  T prev_l2 = 0;
  T prev_h1 = 0;
  T prev_h2 = 0;

  for (index_t ref = 0; ref <= maxRefinements; ++ref) {
    try {
      // ---- Read geometry ----
      gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
      if (!mpPtr) {
        gsInfo << "Error: Cannot read geometry " << geometry << "\n";
        return -1;
      }
      gsMultiPatch<T> &mp = *mpPtr;
      mp.computeTopology();

      // Bounding box of geometry for the interior bubble rectangle
      gsMatrix<T> bbox;
      mp.boundingBox(bbox);
      T x_min = bbox(0, 0);
      T x_max = bbox(0, 1);
      T y_min = bbox(1, 0);
      T y_max = bbox(1, 1);

      T Lx = (x_max - x_min > 1e-12) ? (x_max - x_min) : T(1);
      T Ly = (y_max - y_min > 1e-12) ? (y_max - y_min) : T(1);

      // ---- Degree Elevation ----
      const short_t inputDeg = mp.patch(0).basis().degree(0);
      if (inputDeg < degree) {
        const short_t elev = degree - inputDeg;
        mp.degreeElevate(elev);
      }

      // ---- Refinement ----
      const short_t deg = mp.patch(0).basis().degree(0);
      const index_t mult = std::max<index_t>(deg - 1, 1);
      for (index_t i = 0; i < ref; ++i)
        mp.uniformRefine(1, mult);

      const T h_max = std::pow(0.5, ref);

      // ---- Compute Gluing Data ----
      gsMatrix<T> gd = computeGluingData(mp, T(1e-8), numGaussPerSpan);

      // ---- Build per-patch interface embeddings ----
      std::vector<gsArgyrisEmbedding<T>> argBasis;
      gsVector<index_t> patchDofSizes(mp.nPatches());
      for (size_t i = 0; i < mp.nPatches(); ++i) {
        argBasis.push_back(deriveArgyrisBasisEmbedding(
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis()),
            gsMatrix<T>(gd.row(i)), mp.patch(i)));
        patchDofSizes[i] = argBasis[i].matrix.cols();
      }

      // ---- Boundary conditions ----
      gsConstantFunction<T> zero;
      gsBoundaryConditions<T> bc;
      for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
          bc.add(it->patch, it->side(), "ValuesAndDerivatives", zero);

      // ---- Set up gsDofMapper ----
      gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis, bc);

      const index_t nFree = mapper.freeSize();
      const index_t nBnd = mapper.boundarySize();

      // ---- Build global transformation matrices T_free and T_global ----
      index_t nDisjointBSpline = 0;
      for (size_t i = 0; i < mp.nPatches(); ++i)
        nDisjointBSpline += argBasis[i].matrix.rows();

      gsSparseEntries<T> tFreeEntries, tGlobalEntries;
      index_t rowOffset = 0;

      for (size_t i = 0; i < mp.nPatches(); ++i) {
        const gsSparseMatrix<T> &Ai = argBasis[i].matrix; // nBSpline_i x nArg_i
        const index_t patchSize = mapper.patchSize(i);

        for (index_t j = 0; j < patchSize; ++j) {
          bool isBnd = mapper.is_boundary(j, i);
          index_t gIdx = isBnd ? mapper.bindex(j, i) : mapper.index(j, i);

          for (typename gsSparseMatrix<T>::InnerIterator it(Ai, j); it; ++it) {
            const index_t bsplineRow = rowOffset + it.row();
            const T val = it.value();

            if (isBnd) {
              tGlobalEntries.add(bsplineRow, nFree + gIdx, val);
            } else {
              tFreeEntries.add(bsplineRow, gIdx, val);
              tGlobalEntries.add(bsplineRow, gIdx, val);
            }
          }
        }
        rowOffset += Ai.rows();
      }

      gsSparseMatrix<T> T_free(nDisjointBSpline, nFree);
      T_free.setFrom(tFreeEntries);

      gsSparseMatrix<T> T_global(nDisjointBSpline, nFree + nBnd);
      T_global.setFrom(tGlobalEntries);

      gsMultiBasis<T> dbasis(mp);

      // ----------------------------------------------------------------
      // Manufactured solution: tensor-product B-spline "bubble".
      //
      // We fit a rectangle [a_x,b_x] x [a_y,b_y] inside the domain by
      // shrinking the bounding box by the requested margin, so the
      // rectangle covers the majority of the domain while staying clear of
      // the boundary.
      //
      // On this rectangle the bubble is a SINGLE tensor-product B-spline
      // basis function of degree p_bubble. Following the requested
      // construction, its (uniform) knots on [a,b] are
      //     xi_i = a*(p+1-i)/(p+1) + b*i/(p+1),  i = 0..p+1
      // in each direction, i.e. p+2 simple knots. This single basis
      // function is C^{p-1} and vanishes together with its derivatives up
      // to order p-1 at the rectangle edges, so extending it by zero gives
      // a globally C^{p-1} bubble. Consequently u = 0 and grad(u).n = 0
      // hold exactly on the domain boundary (for p >= 2).
      //
      // To be able to evaluate the bubble at every quadrature point of the
      // domain, we embed those p+2 bump knots into a clamped knot vector
      // that spans the whole (slightly padded) bounding box and activate
      // only the central bump coefficient. The resulting function equals
      // the desired bump on [a,b] and is identically zero elsewhere, so it
      // is safe to evaluate anywhere inside the bounding box.
      // ----------------------------------------------------------------
      const index_t p_bubble = (bubbleDegree > 0) ? bubbleDegree : deg;

      const T insc_Lx = inscribed_xmax - inscribed_xmin;
      const T insc_Ly = inscribed_ymax - inscribed_ymin;
      const T a_x = inscribed_xmin + margin * insc_Lx;
      const T b_x = inscribed_xmax - margin * insc_Lx;
      const T a_y = inscribed_ymin + margin * insc_Ly;
      const T b_y = inscribed_ymax - margin * insc_Ly;

      // Pad the outer clamp slightly beyond the bounding box so that
      // quadrature points lying exactly on the boundary never fall outside
      // the knot domain (which would be an out-of-range B-spline eval).
      const T pad_x = 0.05 * Lx;
      const T pad_y = 0.05 * Ly;

      // Clamped knot vector on [lo, hi] of degree p that contains the p+2
      // uniform bump knots xi_0=a..xi_{p+1}=b as simple interior knots. The
      // basis function of index p+1 has support exactly [a, b] and is the
      // requested degree-p bubble.
      auto makeBubbleKV = [](T lo, T hi, T a, T b, index_t p) {
        std::vector<T> kn;
        kn.reserve(3 * p + 4);
        for (index_t i = 0; i <= p; ++i) // p+1 clamped copies of lo
          kn.push_back(lo);
        for (index_t i = 0; i <= p + 1; ++i) // xi_0=a .. xi_{p+1}=b (simple)
          kn.push_back(a * T(p + 1 - i) / T(p + 1) + b * T(i) / T(p + 1));
        for (index_t i = 0; i <= p; ++i) // p+1 clamped copies of hi
          kn.push_back(hi);
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

      // Replicate the bubble as one piece per patch so it can be used as a
      // physical-coordinate coefficient over the whole multi-patch domain.
      gsPiecewiseFunction<T> bubbleField;
      for (size_t i = 0; i < mp.nPatches(); ++i)
        bubbleField.addPiece(bubble);

      // Laplacian of the bubble (as a value coefficient) for the RHS.
      gsLaplacianFunction<T> bubbleLap(bubble);
      gsPiecewiseFunction<T> bubbleLapField;
      for (size_t i = 0; i < mp.nPatches(); ++i)
        bubbleLapField.addPiece(bubbleLap);

      if (ref == 0) {
        gsInfo << "Bubble rectangle: [" << a_x << ", " << b_x << "] x [" << a_y
               << ", " << b_y << "]  (degree " << p_bubble
               << ", single bump embedded in a " << bubbleBasis.size()
               << "-function basis)\n\n"
               << std::flush;
      }

      // ---- Assemble disjoint biharmonic matrix and energy-form RHS ----
      // Because the bubble vanishes together with its gradient on the
      // boundary, the biharmonic weak form reduces to
      //     a(u, v) = (Delta^2 u_ex, v) = (Delta u_ex, Delta v),
      // so the right-hand side only requires second derivatives of the
      // manufactured solution.
      gsExprAssembler<T> A(1, 1);
      A.setIntegrationDomain(dbasis.domain());

      auto G_map = A.getMap(mp);
      auto u_space = A.getSpace(dbasis);
      auto lap_ex = A.getCoeff(bubbleLapField, G_map);

      A.initSystem();
      A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                 ilapl(u_space, G_map) * lap_ex * meas(G_map));

      const gsSparseMatrix<T> &K_disjoint = A.matrix();
      const gsMatrix<T> &F_disjoint = A.rhs();

      gsSparseMatrix<T> K_free = T_free.transpose() * K_disjoint * T_free;
      gsMatrix<T> F_free = T_free.transpose() * F_disjoint;

      gsSparseSolver<T>::LU solver(K_free);
      gsMatrix<T> c_free = solver.solve(F_free);

      gsMatrix<T> c_global(nFree + nBnd, 1);
      c_global.topRows(nFree) = c_free;
      c_global.bottomRows(nBnd).setZero();

      gsMatrix<T> c_disjoint = T_global * c_global;

      gsMultiPatch<T> sol;
      index_t offset = 0;
      for (size_t i = 0; i < mp.nPatches(); ++i) {
        const index_t sz = argBasis[i].matrix.rows();
        gsMatrix<T> ci = c_disjoint.block(offset, 0, sz, 1);
        offset += sz;
        const gsTensorBSplineBasis<2, T> &tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis());
        sol.addPatch(tb.makeGeometry(give(ci)));
      }

      // ---- Error Evaluation ----
      gsExprEvaluator<T> ev;
      ev.setIntegrationDomain(dbasis.domain());
      auto G_map_ev = ev.getMap(mp);

      auto u_exact_ev = ev.getVariable(bubbleField, G_map_ev);
      auto u_sol_ev = ev.getVariable(sol);

      const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
      const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev)).sqNorm() * meas(G_map_ev)));
      const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev)).sqNorm() * meas(G_map_ev)));

      // Rates
      T l2_rate = 0, h1_rate = 0, h2_rate = 0;
      if (ref > 2 && prev_h > 0) {
        const T h_ratio = prev_h / h_max;
        l2_rate = std::log(prev_l2 / l2err) / std::log(h_ratio);
        h1_rate = std::log(prev_h1 / h1err) / std::log(h_ratio);
        h2_rate = std::log(prev_h2 / h2err) / std::log(h_ratio);
      }

      gsInfo << std::setw(5) << ref
             << std::setw(12) << std::scientific << std::setprecision(4) << h_max
             << std::setw(10) << nFree
             << std::setw(14) << std::scientific << std::setprecision(4) << l2err;
      if (ref > 2)
        gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << l2_rate;
      else
        gsInfo << std::setw(10) << "-";

      gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
      if (ref > 2)
        gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << h1_rate;
      else
        gsInfo << std::setw(10) << "-";

      gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
      if (ref > 2)
        gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << h2_rate;
      else
        gsInfo << std::setw(10) << "-";

      gsInfo << "\n" << std::flush;

      prev_h = h_max;
      prev_l2 = l2err;
      prev_h1 = h1err;
      prev_h2 = h2err;

      if (plot || !outDir.empty()) {
        std::string solName = prefix + "as_g1_biharmonic_bubble_sol_r" + std::to_string(ref);
        std::string exactName = prefix + "as_g1_biharmonic_bubble_exact_r" + std::to_string(ref);

        gsField<T> solField(mp, sol);
        gsWriteParaview<>(solField, solName, 1000);

        gsField<T> exactField(mp, bubbleField, false);
        gsWriteParaview<>(exactField, exactName, 1000);

        gsInfo << " -> Exported VTK plots: " << solName << " & " << exactName << "\n";
      }
    } catch (const std::exception &e) {
      gsInfo << "ERROR: " << e.what() << "\n" << std::flush;
    }
  }

  gsInfo << std::string(99, '-') << "\n";
  gsInfo << "Done.\n";
  return 0;
}
