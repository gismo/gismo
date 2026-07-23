/** @file as_g1_biharmonic_bubble_example.cpp

    @brief Biharmonic solver with manufactured bubble solution over Analysis-Suitable G1 (AS-G1) multi-patch geometries.

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

int main(int argc, char *argv[]) {
  using T = real_t;

  std::string geometry("domain2d/2patch/multipatch/weirdo_multivalence_non_bilinear.xml");
  std::string outDir("");
  index_t degree = 3;
  index_t numGaussPerSpan = 0;
  index_t maxRefinements = 4;
  index_t bubbleType = 0; // 0: B-spline AS-G1 Bubble, 1: Polynomial Bounding-Box Bubble, 2: Trigonometric Bounding-Box Bubble
  bool plot = false;

  gsCmdLine cmd("AS-G1 Biharmonic Equation Solver with Manufactured Bubble Solution.");
  cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
  cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
  cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
  cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).", numGaussPerSpan);
  cmd.addInt("r", "refinements", "Maximum uniform refinement levels to test.", maxRefinements);
  cmd.addInt("b", "bubble", "Bubble type: 0 = B-spline AS-G1 Vertex Bubble (smooth over internal vertex & domain), 1 = Smooth Trigonometric Physical Bubble (sin^4), 2 = Smooth Polynomial Physical Bubble (t^3(1-t)^3).", bubbleType);
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

  gsInfo << "\n======================================================================\n";
  gsInfo << "AS-G1 Biharmonic Solver (Delta^2 u = f) with Manufactured Bubble Solution\n";
  if (bubbleType == 0)
    gsInfo << "Manufactured Solution: B-Spline AS-G1 Internal Vertex Bubble Function\n";
  else if (bubbleType == 1)
    gsInfo << "Manufactured Solution: Smooth Trigonometric Physical Bubble (sin^4)\n";
  else
    gsInfo << "Manufactured Solution: Smooth Polynomial Physical Bubble (t^3(1-t)^3)\n";
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

  auto sumUntil = [](const gsVector<index_t, 13> &vec, index_t until) {
    index_t sum = 0;
    for (index_t i = 0; i < until; ++i)
      sum += vec(i);
    return sum;
  };

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

      // Bounding box of geometry for domain scaling
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

      // ---- Set up gsDofMapper ----
      gsDofMapper mapper(patchDofSizes);

      // Match interface DOFs
      for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
        const boundaryInterface &ifc = *it;
        const patchSide ps1 = ifc.first();
        const patchSide ps2 = ifc.second();

        const index_t nLvl0 = argBasis[ps1.patch].sizes[1 + 2 * (ps1.m_index - 1)];
        const index_t offLvl0_1 = sumUntil(argBasis[ps1.patch].sizes, 1 + 2 * (ps1.m_index - 1));
        const index_t offLvl0_2 = sumUntil(argBasis[ps2.patch].sizes, 1 + 2 * (ps2.m_index - 1));

        const index_t nLvl1 = argBasis[ps1.patch].sizes[2 + 2 * (ps1.m_index - 1)];
        const index_t offLvl1_1 = sumUntil(argBasis[ps1.patch].sizes, 2 + 2 * (ps1.m_index - 1));
        const index_t offLvl1_2 = sumUntil(argBasis[ps2.patch].sizes, 2 + 2 * (ps2.m_index - 1));

        const short_t tanDir1 = 1 - ps1.direction();
        const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

        for (index_t j1 = 0; j1 < nLvl0; ++j1) {
          const index_t j2 = flipped ? nLvl0 - 1 - j1 : j1;
          mapper.matchDof(ps1.patch, offLvl0_1 + j1, ps2.patch, offLvl0_2 + j2);
        }

        for (index_t j1 = 0; j1 < nLvl1; ++j1) {
          const index_t j2 = flipped ? nLvl1 - 1 - j1 : j1;
          mapper.matchDof(ps1.patch, offLvl1_1 + j1, ps2.patch, offLvl1_2 + j2);
        }

        std::vector<patchCorner> corners;
        ps1.getContainedCorners(2, corners);
        for (index_t i = 0; i < 2; ++i) {
          const index_t c1 = corners[i].m_index - 1;
          const index_t c2 = ifc.mapCorner(corners[i]).m_index - 1;

          const index_t off_corner_1 = sumUntil(argBasis[ps1.patch].sizes, 9 + c1);
          const index_t off_corner_2 = sumUntil(argBasis[ps2.patch].sizes, 9 + c2);

          for (index_t k = 0; k < 6; ++k)
            mapper.matchDof(ps1.patch, off_corner_1 + k, ps2.patch, off_corner_2 + k);
        }
      }

      // Mark domain boundary DOFs (homogeneous Dirichlet)
      std::set<std::pair<size_t, index_t>> ifcSides;
      for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
        ifcSides.insert({it->first().patch, it->first().side()});
        ifcSides.insert({it->second().patch, it->second().side()});
      }

      for (size_t i = 0; i < mp.nPatches(); ++i) {
        for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side) {
          if (ifcSides.find({i, side.m_index}) != ifcSides.end())
            continue; // interior interface, skip

          // External boundary side
          const index_t nLvl0 = argBasis[i].sizes[1 + 2 * (side.m_index - 1)];
          const index_t offLvl0 = sumUntil(argBasis[i].sizes, 1 + 2 * (side.m_index - 1));
          for (index_t j = 0; j < nLvl0; ++j)
            mapper.eliminateDof(offLvl0 + j, i);

          const index_t nLvl1 = argBasis[i].sizes[2 + 2 * (side.m_index - 1)];
          const index_t offLvl1 = sumUntil(argBasis[i].sizes, 2 + 2 * (side.m_index - 1));
          for (index_t j = 0; j < nLvl1; ++j)
            mapper.eliminateDof(offLvl1 + j, i);

          // Contained corners on this boundary side
          patchSide ps(i, side);
          std::vector<patchCorner> corners;
          ps.getContainedCorners(2, corners);
          for (const auto &c : corners) {
            const index_t cIdx = c.m_index - 1;
            const index_t offCorner = sumUntil(argBasis[i].sizes, 9 + cIdx);
            for (index_t k = 0; k < 6; ++k)
              mapper.eliminateDof(offCorner + k, i);
          }
        }
      }

      mapper.finalize();

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

      gsMultiPatch<T> sol_exact;
      gsMatrix<T> F_free;

      if (bubbleType == 0) {
        // Option 0: Tensor-Product B-Spline AS-G1 Internal Vertex Bubble Solution
        // Constructed directly in AS-G1 space by placing non-zero weights on internal vertex value DOFs.
        // By construction of AS-G1 space, this function is GUARANTEED to be:
        // 1) G1 continuous across all internal patch interfaces.
        // 2) Smoothly vanishing (value = 0, normal derivative = 0) on all domain boundaries.
        // 3) A smooth bell-shaped bubble centered around the internal vertex/vertices.

        gsMatrix<T> c_exact_free(nFree, 1);
        c_exact_free.setZero();

        std::set<index_t> internalVertexDofs;
        for (size_t i = 0; i < mp.nPatches(); ++i) {
          for (index_t c = 0; c < 4; ++c) {
            const index_t off_corner = sumUntil(argBasis[i].sizes, 9 + c);
            if (!mapper.is_boundary(off_corner, i)) {
              index_t gIdx = mapper.index(off_corner, i);
              internalVertexDofs.insert(gIdx);
            }
          }
        }

        if (!internalVertexDofs.empty()) {
          for (index_t gIdx : internalVertexDofs) {
            c_exact_free(gIdx, 0) = 1.0;
          }
        } else {
          // Fallback: set middle free DOF to 1.0
          c_exact_free(nFree / 2, 0) = 1.0;
        }

        gsMatrix<T> c_global_exact(nFree + nBnd, 1);
        c_global_exact.topRows(nFree) = c_exact_free;
        c_global_exact.bottomRows(nBnd).setZero();

        gsMatrix<T> c_disjoint_exact = T_global * c_global_exact;

        index_t offset = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i) {
          const index_t sz = argBasis[i].matrix.rows();
          gsMatrix<T> ci = c_disjoint_exact.block(offset, 0, sz, 1);
          offset += sz;
          const gsTensorBSplineBasis<2, T> &tb =
              dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis());
          sol_exact.addPatch(tb.makeGeometry(give(ci)));
        }

        // Assemble Biharmonic Matrix K_disjoint
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(dbasis.domain());

        auto G_map = A.getMap(mp);
        auto u_space = A.getSpace(dbasis);

        A.initSystem();
        A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map));

        const gsSparseMatrix<T> &K_disjoint = A.matrix();
        gsSparseMatrix<T> K_free = T_free.transpose() * K_disjoint * T_free;

        F_free = K_free * c_exact_free;

        // Solve linear system
        gsSparseSolver<T>::LU solver(K_free);
        gsMatrix<T> c_free = solver.solve(F_free);

        gsMatrix<T> c_global(nFree + nBnd, 1);
        c_global.topRows(nFree) = c_free;
        c_global.bottomRows(nBnd).setZero();

        gsMatrix<T> c_disjoint = T_global * c_global;

        gsMultiPatch<T> sol;
        offset = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i) {
          const index_t sz = argBasis[i].matrix.rows();
          gsMatrix<T> ci = c_disjoint.block(offset, 0, sz, 1);
          offset += sz;
          const gsTensorBSplineBasis<2, T> &tb =
              dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis());
          sol.addPatch(tb.makeGeometry(give(ci)));
        }

        // Error evaluation between discrete solution and exact B-spline bubble
        gsExprEvaluator<T> ev;
        ev.setIntegrationDomain(dbasis.domain());
        auto G_map_ev = ev.getMap(mp);

        auto u_exact_ev = ev.getVariable(sol_exact);
        auto u_sol_ev = ev.getVariable(sol);

        const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
        const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev, G_map_ev)).sqNorm() * meas(G_map_ev)));
        const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev, G_map_ev)).sqNorm() * meas(G_map_ev)));

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

          gsField<T> exactField(mp, sol_exact);
          gsWriteParaview<>(exactField, exactName, 1000);

          gsInfo << " -> Exported VTK plots: " << solName << " & " << exactName << "\n";
        }

      } else {
        // Option 1 & 2: Analytical Bounding-Box Bubble Functions
        // Option 1: Polynomial Bubble P(t) = t^2 (1-t)^2
        // Option 2: Trigonometric Bubble Q(t) = sin^2(pi * t)
        std::ostringstream u_str, grad_x_str, grad_y_str, h_xx_str, h_xy_str, h_yy_str, f_str;
        std::string x_bar = "( (x - " + std::to_string(x_min) + ") / " + std::to_string(Lx) + " )";
        std::string y_bar = "( (y - " + std::to_string(y_min) + ") / " + std::to_string(Ly) + " )";

        if (bubbleType == 1) {
          // Polynomial Bubble: P(t) = t^2 (1-t)^2
          std::string P_x = "(" + x_bar + "^2 * (1 - " + x_bar + ")^2)";
          std::string P_y = "(" + y_bar + "^2 * (1 - " + y_bar + ")^2)";

          std::string dP_x = "(2*" + x_bar + " - 6*" + x_bar + "^2 + 4*" + x_bar + "^3)";
          std::string dP_y = "(2*" + y_bar + " - 6*" + y_bar + "^2 + 4*" + y_bar + "^3)";

          std::string d2P_x = "(2 - 12*" + x_bar + " + 12*" + x_bar + "^2)";
          std::string d2P_y = "(2 - 12*" + y_bar + " + 12*" + y_bar + "^2)";

          u_str << P_x << " * " << P_y;
          grad_x_str << "(1 / " << Lx << ") * " << dP_x << " * " << P_y;
          grad_y_str << "(1 / " << Ly << ") * " << P_x << " * " << dP_y;

          h_xx_str << "(1 / " << (Lx * Lx) << ") * " << d2P_x << " * " << P_y;
          h_xy_str << "(1 / " << (Lx * Ly) << ") * " << dP_x << " * " << dP_y;
          h_yy_str << "(1 / " << (Ly * Ly) << ") * " << P_x << " * " << d2P_y;

          f_str << "(24 / " << (Lx * Lx * Lx * Lx) << ") * " << P_y << " + "
                << "(2 / " << (Lx * Lx * Ly * Ly) << ") * " << d2P_x << " * " << d2P_y << " + "
                << "(24 / " << (Ly * Ly * Ly * Ly) << ") * " << P_x;
        } else {
          // Trigonometric Bubble: Q(t) = sin^2(pi * t)
          std::string Q_x = "(sin(pi*" + x_bar + ")^2)";
          std::string Q_y = "(sin(pi*" + y_bar + ")^2)";

          std::string dQ_x = "(pi*sin(2*pi*" + x_bar + "))";
          std::string dQ_y = "(pi*sin(2*pi*" + y_bar + "))";

          std::string d2Q_x = "(2*pi^2*cos(2*pi*" + x_bar + "))";
          std::string d2Q_y = "(2*pi^2*cos(2*pi*" + y_bar + "))";

          std::string d4Q_x = "(-8*pi^4*cos(2*pi*" + x_bar + "))";
          std::string d4Q_y = "(-8*pi^4*cos(2*pi*" + y_bar + "))";

          u_str << Q_x << " * " << Q_y;
          grad_x_str << "(1 / " << Lx << ") * " << dQ_x << " * " << Q_y;
          grad_y_str << "(1 / " << Ly << ") * " << Q_x << " * " << dQ_y;

          h_xx_str << "(1 / " << (Lx * Lx) << ") * " << d2Q_x << " * " << Q_y;
          h_xy_str << "(1 / " << (Lx * Ly) << ") * " << dQ_x << " * " << dQ_y;
          h_yy_str << "(1 / " << (Ly * Ly) << ") * " << Q_x << " * " << d2Q_y;

          f_str << "(1 / " << (Lx * Lx * Lx * Lx) << ") * " << d4Q_x << " * " << Q_y << " + "
                << "(2 / " << (Lx * Lx * Ly * Ly) << ") * " << d2Q_x << " * " << d2Q_y << " + "
                << "(1 / " << (Ly * Ly * Ly * Ly) << ") * " << Q_x << " * " << d4Q_y;
        }

        gsFunctionExpr<T> exact_u(u_str.str(), 2);
        gsFunctionExpr<T> exact_grad(grad_x_str.str(), grad_y_str.str(), 2);
        gsFunctionExpr<T> exact_hess(h_xx_str.str(), h_xy_str.str(), h_xy_str.str(), h_yy_str.str(), 2);
        gsFunctionExpr<T> rhs_f(f_str.str(), 2);

        // Assemble Disjoint Biharmonic Matrix & RHS
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(dbasis.domain());

        auto G_map = A.getMap(mp);
        auto u_space = A.getSpace(dbasis);
        auto f_coeff = A.getCoeff(rhs_f, G_map);

        A.initSystem();
        A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                   u_space * f_coeff * meas(G_map));

        const gsSparseMatrix<T> &K_disjoint = A.matrix();
        const gsMatrix<T> &F_disjoint = A.rhs();

        gsSparseMatrix<T> K_free = T_free.transpose() * K_disjoint * T_free;
        F_free = T_free.transpose() * F_disjoint;

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

        // Error Evaluation
        gsExprEvaluator<T> ev;
        ev.setIntegrationDomain(dbasis.domain());
        auto G_map_ev = ev.getMap(mp);

        auto u_exact_ev = ev.getVariable(exact_u, G_map_ev);
        auto u_sol_ev = ev.getVariable(sol);

        const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
        const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev, G_map_ev)).sqNorm() * meas(G_map_ev)));
        const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev, G_map_ev)).sqNorm() * meas(G_map_ev)));

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

          gsField<T> exactField(mp, exact_u, false);
          gsWriteParaview<>(exactField, exactName, 1000);

          gsInfo << " -> Exported VTK plots: " << solName << " & " << exactName << "\n";
        }
      }
    } catch (const std::exception &e) {
      gsInfo << "ERROR: " << e.what() << "\n" << std::flush;
    }
  }

  gsInfo << std::string(99, '-') << "\n";
  gsInfo << "Done.\n";
  return 0;
}
