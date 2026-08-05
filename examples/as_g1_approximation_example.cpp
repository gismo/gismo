/** @file as_g1_approximation_example.cpp

    @brief L2 function approximation and convergence test over Analysis-Suitable G1 (AS-G1) multi-patch geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <algorithm>
#include <cmath>
#include <gismo.h>
#include <iomanip>
#include <iostream>
#include <vector>

#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

using namespace gismo;

int main(int argc, char *argv[]) {
  using T = real_t;

  std::string geometry("domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
  std::string outDir("");
  index_t degree = 3;
  index_t numGaussPerSpan = 0;
  index_t maxRefinements = 5;
  bool plot = false;

  gsCmdLine cmd("AS-G1 Function Approximation and Convergence Experiment.");
  cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
  cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
  cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
  cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).", numGaussPerSpan);
  cmd.addInt("r", "refinements", "Maximum uniform refinement levels to test.", maxRefinements);
  cmd.addSwitch("plot", "Export target function and numerical solution to VTK files.", plot);

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
  gsInfo << "Target Exact Function: f(x, y) = sin(pi * x) * cos(pi * y)\n";
  gsInfo << "======================================================================\n\n";

  // Table header
  gsInfo << std::setw(5) << "r"
         << std::setw(12) << "h_max"
         << std::setw(12) << "N_global"
         << std::setw(16) << "L2 Error"
         << std::setw(12) << "L2 Rate"
         << std::setw(16) << "H1 Error"
         << std::setw(12) << "H1 Rate"
         << "\n";
  gsInfo << std::string(85, '-') << "\n" << std::flush;

  T prev_h = 0;
  T prev_l2 = 0;
  T prev_h1 = 0;

  for (index_t ref = 2; ref <= maxRefinements; ++ref) {
    try {
    // ---- Read geometry ----
    gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) {
      gsInfo << "Error: Cannot read geometry " << geometry << "\n";
      return -1;
    }
    gsMultiPatch<T> &mp = *mpPtr;
    mp.computeTopology();

    // ---- Degree Elevation ----
    const short_t inputDeg = mp.patch(0).basis().degree(0);
    if (inputDeg < degree) {
      const short_t elev = degree - inputDeg;
      mp.degreeElevate(elev);
    }

    gsInfo << "Refinement level " << ref << " ... " << std::flush;

    // ---- Refinement ----
    const short_t deg = mp.patch(0).basis().degree(0);
    const index_t mult = std::max<index_t>(deg - 1, 1);
    for (index_t i = 0; i < ref; ++i)
      mp.uniformRefine(1, mult);

    // Estimate h_max (characteristic mesh size on parameter domain / refined elements)
    const T h_max = std::pow(0.5, ref);

    // ---- Compute Gluing Data ----
    gsMatrix<T> gd = computeGluingData(mp, T(1e-8), numGaussPerSpan);

    // ---- Build per-patch Argyris embeddings ----
    std::vector<gsArgyrisEmbedding<T>> argBasis;
    for (size_t i = 0; i < mp.nPatches(); ++i)
      argBasis.push_back(deriveArgyrisBasisEmbedding(
          dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis()),
          gsMatrix<T>(gd.row(i)), mp.patch(i)));

    // ---- AS-G1 DOF coupling (no boundary conditions for L2 approximation) ----
    gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis);
    const index_t nGlobal = mapper.freeSize();

    // ---- Stacked transformation matrix (disjoint B-spline DOFs <- global DOFs) ----
    index_t nDisjoint = 0;
    for (size_t i = 0; i < mp.nPatches(); ++i)
      nDisjoint += argBasis[i].matrix.rows();

    gsSparseEntries<T> tEntries;
    index_t rowOffset = 0;
    for (size_t i = 0; i < mp.nPatches(); ++i) {
      const gsSparseMatrix<T> &Ai = argBasis[i].matrix;
      for (index_t j = 0; j < mapper.patchSize(i); ++j) {
        const index_t gIdx = mapper.index(j, i);
        for (typename gsSparseMatrix<T>::InnerIterator it(Ai, j); it; ++it)
          tEntries.add(rowOffset + it.row(), gIdx, it.value());
      }
      rowOffset += Ai.rows();
    }
    gsSparseMatrix<T> T_global(nDisjoint, nGlobal);
    T_global.setFrom(tEntries);

    // ---- Assemble Patch Mass Matrix & Load Vector ----
    gsMultiBasis<T> dbasis(mp);
    gsExprAssembler<T> A(1, 1);
    A.setIntegrationDomain(dbasis.domain());

    auto G_map = A.getMap(mp);
    auto u_space = A.getSpace(dbasis);

    gsFunctionExpr<T> exact_f("sin(pi*x)*cos(pi*y)", 2);
    auto f_coeff = A.getCoeff(exact_f, G_map);

    A.initSystem();
    A.assemble(u_space * u_space.tr() * meas(G_map), u_space * f_coeff * meas(G_map));

    const gsSparseMatrix<T> &M_disjoint = A.matrix();
    const gsMatrix<T> &F_disjoint = A.rhs();

    // ---- Global System Assembly & Solve ----
    gsSparseMatrix<T> M_global = T_global.transpose() * M_disjoint * T_global;
    gsMatrix<T> F_global = T_global.transpose() * F_disjoint;

    gsSparseSolver<T>::LU solver(M_global);
    gsMatrix<T> c_global = solver.solve(F_global);

    // ---- Reconstruct Multi-patch Solution Field ----
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
    auto f_exact_ev = ev.getVariable(exact_f, G_map_ev);
    auto u_sol_ev = ev.getVariable(sol);

    const T l2err = std::sqrt(ev.integral((u_sol_ev - f_exact_ev).sqNorm() * meas(G_map_ev)));

    gsFunctionExpr<T> exact_grad("pi*cos(pi*x)*cos(pi*y)", "-pi*sin(pi*x)*sin(pi*y)", 2);
    auto grad_exact_ev = ev.getVariable(exact_grad, G_map_ev);

    const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() * meas(G_map_ev)));

    // Rates
    T l2_rate = 0;
    T h1_rate = 0;
    if (ref > 2 && prev_h > 0) {
      const T h_ratio = prev_h / h_max;
      l2_rate = std::log(prev_l2 / l2err) / std::log(h_ratio);
      h1_rate = std::log(prev_h1 / h1err) / std::log(h_ratio);
    }

    gsInfo << std::setw(5) << ref
           << std::setw(12) << std::scientific << std::setprecision(4) << h_max
           << std::setw(12) << nGlobal
           << std::setw(16) << std::scientific << std::setprecision(6) << l2err;
    if (ref > 2)
      gsInfo << std::setw(12) << std::fixed << std::setprecision(2) << l2_rate;
    else
      gsInfo << std::setw(12) << "-";

    gsInfo << std::setw(16) << std::scientific << std::setprecision(6) << h1err;
    if (ref > 2)
      gsInfo << std::setw(12) << std::fixed << std::setprecision(2) << h1_rate;
    else
      gsInfo << std::setw(12) << "-";

    gsInfo << "\n" << std::flush;

    prev_h = h_max;
    prev_l2 = l2err;
    prev_h1 = h1err;

    // Optional VTK output for original target function and approximated solution
    if (plot || !outDir.empty()) {
      std::string solName = prefix + "as_g1_approx_sol_r" + std::to_string(ref);
      std::string exactName = prefix + "as_g1_approx_exact_r" + std::to_string(ref);

      gsField<T> solField(mp, sol);
      gsWriteParaview<>(solField, solName, 1000);

      gsField<T> exactField(mp, exact_f, false);
      gsWriteParaview<>(exactField, exactName, 1000);

      gsInfo << " -> Exported VTK plots: " << solName << " & " << exactName << "\n";
    }
    } catch (const std::exception &e) {
      gsInfo << "ERROR: " << e.what() << "\n" << std::flush;
    }
  }

  gsInfo << std::string(85, '-') << "\n";
  gsInfo << "Done.\n";
  return 0;
}
