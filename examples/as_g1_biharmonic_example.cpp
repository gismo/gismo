/** @file as_g1_biharmonic_example.cpp

    @brief Biharmonic solver and convergence test over Analysis-Suitable G1
   (AS-G1) multi-patch geometries.

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
using namespace gismo::expr;

int main(int argc, char *argv[]) {
  using T = real_t;

  std::string geometry(
      "domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
  std::string outDir("");
  index_t degree = 3;
  index_t numGaussPerSpan = 0;
  index_t minRefinements = 2;
  index_t maxRefinements = 4;
  index_t freqA = 1;
  bool plot = false;

  gsCmdLine cmd("AS-G1 Biharmonic Equation Solver and Convergence Experiment.");
  cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
  cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
  cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
  cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).",
             numGaussPerSpan);
  cmd.addInt("m", "minRefinements",
             "Minimum uniform refinement level to test (minimum 2 for AS-G1).",
             minRefinements);
  cmd.addInt("r", "refinements", "Maximum uniform refinement levels to test.",
             maxRefinements);
  cmd.addInt("a", "frequency",
             "Frequency integer factor 'a' in sin(a*pi*x)*cos(a*pi*y).",
             freqA);
  cmd.addSwitch("plot",
                "Export target function and numerical solution to VTK files.",
                plot);

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

  if (minRefinements < 2)
    minRefinements = 2;

  if (freqA < 1)
    freqA = 1;

  // ---- Manifactured Solution ----
  const std::string s_a = std::to_string(freqA);
  const std::string u_expr = "sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
  const std::string grad_x_expr =       s_a + "*pi*cos(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
  const std::string grad_y_expr = "-" + s_a + "*pi*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
  const std::string rhs_expr = std::to_string(4 * freqA * freqA * freqA * freqA) +
      "*pi^4*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
  const std::string hess_xx = "-" + std::to_string(freqA * freqA) + "*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
  const std::string hess_xy = "-" + std::to_string(freqA * freqA) + "*pi^2*cos(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
  const std::string hess_yy = "-" + std::to_string(freqA * freqA) + "*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
  gsFunctionExpr<T> exact_u(u_expr, 2);
  gsFunctionExpr<T> exact_grad(grad_x_expr, grad_y_expr, 2);
  gsFunctionExpr<T> exact_hess(hess_xx, hess_xy, hess_xy, hess_yy, 2);

  gsInfo << "\n======================================================================\n";
  gsInfo << "AS-G1 Biharmonic Solver (Delta^2 u = f)\n";
  gsInfo << "Target Exact Function: u(x, y) = " << u_expr << "\n";
  gsInfo << "======================================================================\n\n";

  // ---- Table header ----
  gsInfo << std::setw(5) << "r" << std::setw(12) << "h_max" << std::setw(10)
         << "N_free" << std::setw(14) << "L2 Error" << std::setw(10)
         << "L2 Rate" << std::setw(14) << "H1 Error" << std::setw(10)
         << "H1 Rate" << std::setw(14) << "H2 Error" << std::setw(10)
         << "H2 Rate"
         << "\n";
  gsInfo << std::string(99, '-') << "\n" << std::flush;

  T prev_h = 0;
  T prev_l2 = 0;
  T prev_h1 = 0;
  T prev_h2 = 0;

  for (index_t ref = minRefinements; ref <= maxRefinements; ++ref) {
      // ---- Read geometry ----
      gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
      if (!mpPtr)
      {
          gsInfo << "Error: Cannot read geometry " << geometry << "\n";
          return -1;
      }
      gsMultiPatch<T> &mp = *mpPtr;
      mp.computeTopology();

      // ---- Degree Elevation ----
      const short_t inputDeg = mp.patch(0).basis().degree(0); //TODO: fixme
      if (inputDeg < degree)
      {
          const short_t elev = degree - inputDeg;
          mp.degreeElevate(elev);
      }

      // ---- Get multi basis ----
      gsMultiBasis<T> dbasis(mp);

      // ---- Refinement ----
      const index_t mult = 2;
      for (index_t i = 0; i < ref; ++i)
        dbasis.uniformRefine(1, mult);

      const T h_max = std::pow(0.5, ref); //TODO: fixme



      // ---- Define boundary conditions ----
      gsConstantFunction<T> zero;
      gsBoundaryConditions<T> bc;
      for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
          bc.add(it->patch, it->side(), "ValuesAndDerivatives", zero);

      // ---- Compute Gluing Data ----
      gsMatrix<T> gd = computeGluingData(mp, T(1e-8), numGaussPerSpan);

      // ---- Get boundary corners and asociated normals ----
      std::vector<std::vector<patchCorner>> vertices = getBoundaryVertices(mp, bc, "ValuesAndDerivatives");
      std::vector<gsVector<T>> normalsAtVertices;
      std::vector<cornerCondition> cc;
      gsMatrix<T> normalsForPatches(mp.nPatches(), 2*4);
      normalsForPatches.setZero();
      for (size_t i=0; i<vertices.size(); ++i)
      {
          gsVector<T> normal = getOuterNormalDerivative(mp, vertices[i]);
          for (size_t j=0; j<vertices[i].size(); ++j)
              normalsForPatches.block(vertices[i][j].patch, 2*(vertices[i][j].m_index-1), 1, 2) = normal.transpose();
          cc.push_back(cornerCondition{vertices[i][0], normal.norm()<1e-6 ? cornerConditionType::all : cornerConditionType::valuesNormals});
          normalsAtVertices.push_back(give(normal));
      }
      // Impose (1,0) as normal vector where we do not have a joint normal vector
      for (size_t i=0; i<mp.nPatches(); ++i)
          for (index_t j=0; j<4; ++j)
              if (normalsForPatches.block(i,2*j,1,2).norm() < 1e-6)
              {
                  normalsForPatches(i,2*j)=1;
                  normalsForPatches(i,2*j+1)=0;
              }

      // ---- Build per-patch interface embeddings ----
      std::vector<gsArgyrisEmbedding<T>> argBasis;
      gsVector<index_t> patchDofSizes(mp.nPatches());
      gsBlockSparseMatrix<T> argBasisGlobalCreator(mp.nPatches(), mp.nPatches());
      for (size_t i = 0; i < mp.nPatches(); ++i)
      {
        argBasis.push_back(deriveArgyrisBasisEmbedding(
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[i]),
            gsMatrix<T>(gd.row(i)),
            gsMatrix<T>(normalsForPatches.row(i)),
            mp.patch(i)
        ));
        patchDofSizes[i] = argBasis[i].matrix.cols();
        argBasisGlobalCreator.set(i,i,argBasis[i].matrix);
      }
      gsSparseMatrix<T> argBasisGlobal = argBasisGlobalCreator;

      // ---- Set up gsDofMapper ----
      gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis, bc, cc);

      // ---- Build global transformation matrices T_free and T_bnd ----
      gsSparseMatrix<T> T_global = argBasisGlobal * asEmbeddingMatrix<T>(mapper.size(), mapper.asVector()).transpose();
      const index_t nFree = mapper.freeSize(), nBnd = mapper.boundarySize();
      gsSparseMatrix<T> T_free = T_global.leftCols(nFree);
      gsSparseMatrix<T> T_bnd = T_global.rightCols(nBnd);

      // ---- Evaluate Dirichlet boundary values vector g_bnd ----
      // TODO: Derive inhomogenous boundary data (edges+vertex) properly
      gsMatrix<T> sol_bnd;
      {
          gsExprAssembler<T> A(1, 1);
          A.setIntegrationDomain(dbasis.domain());
          auto G_map = A.getMap(mp);
          auto u_space = A.getSpace(dbasis);
          auto u_coeff = A.getCoeff(exact_u, G_map);

          A.initSystem();
          A.assemble(u_space * u_space.tr() * meas(G_map),
                     u_space * u_coeff * meas(G_map));

          gsSparseMatrix<T> M_global = T_global.transpose() * A.matrix() * T_global;
          gsMatrix<T> F_global = T_global.transpose() * A.rhs();

          gsMatrix<T> sol_global;
          makeSparseCholeskySolver(M_global)->apply(F_global, sol_global);
          sol_bnd = sol_global.bottomRows(nBnd);
      }

      // ---- Assemble Discontinuous Biharmonic Matrix & RHS ----
      gsSparseMatrix<T> K_discont;
      gsMatrix<T> F_discont;
      {
          gsExprAssembler<T> A(1, 1);
          A.setIntegrationDomain(dbasis.domain());

          auto G_map = A.getMap(mp);
          auto u_space = A.getSpace(dbasis);

          gsFunctionExpr<T> rhs_f(rhs_expr, 2);
          auto f_coeff = A.getCoeff(rhs_f, G_map);

          A.initSystem();
          A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                    u_space * f_coeff * meas(G_map));

          K_discont = give(A.matrix());
          F_discont = give(A.rhs());
      }

      // ---- Global Linear System Assembly & Solve ----
      gsSparseMatrix<T> K_free = T_free.transpose() * K_discont * T_free;
      gsMatrix<T> F_free = T_free.transpose() * (F_discont - K_discont * (T_bnd * sol_bnd));
      gsMatrix<T> sol_free;
      makeSparseCholeskySolver(K_free)->apply(F_free, sol_free);

      // Reconstruct full global vector sol_discont
      gsMatrix<T> sol_discont = T_free * sol_free + T_bnd * sol_bnd;

      // Reconstruct multi-patch solution field
      gsMultiPatch<T> sol;
      index_t offset = 0;
      for (size_t i = 0; i < mp.nPatches(); ++i)
      {
          const index_t sz = argBasis[i].matrix.rows();
          gsMatrix<T> ci = sol_discont.block(offset, 0, sz, 1);
          offset += sz;
          const gsTensorBSplineBasis<2, T> &tb = dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[i]);
          sol.addPatch(tb.makeGeometry(give(ci)));
      }

      // ---- Error Evaluation ----
      gsExprEvaluator<T> ev;
      ev.setIntegrationDomain(dbasis.domain());
      auto G_map_ev = ev.getMap(mp);

      auto u_exact_ev = ev.getVariable(exact_u, G_map_ev);
      auto grad_exact_ev = ev.getVariable(exact_grad, G_map_ev);
      auto hess_exact_ev = reshape(ev.getVariable(exact_hess, G_map_ev), 2, 2);

      auto u_sol_ev = ev.getVariable(sol);

      const T l2err = std::sqrt(
          ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
      const T h1err = std::sqrt(ev.integral(
          (igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() *
          meas(G_map_ev)));
      const T h2err = std::sqrt(
          ev.integral((ihess(u_sol_ev, G_map_ev) - hess_exact_ev).sqNorm() *
                      meas(G_map_ev)));

      // Rates
      T l2_rate = 0, h1_rate = 0, h2_rate = 0;
      if (ref > 2 && prev_h > 0) {
        const T h_ratio = prev_h / h_max;
        l2_rate = std::log(prev_l2 / l2err) / std::log(h_ratio);
        h1_rate = std::log(prev_h1 / h1err) / std::log(h_ratio);
        h2_rate = std::log(prev_h2 / h2err) / std::log(h_ratio);
      }

      gsInfo << std::setw(5) << ref << std::setw(12) << std::scientific
             << std::setprecision(4) << h_max << std::setw(10) << nFree
             << std::setw(14) << std::scientific << std::setprecision(4)
             << l2err;
      if (ref > 2)
        gsInfo << std::setw(10) << std::fixed << std::setprecision(2)
               << l2_rate;
      else
        gsInfo << std::setw(10) << "-";

      gsInfo << std::setw(14) << std::scientific << std::setprecision(4)
             << h1err;
      if (ref > 2)
        gsInfo << std::setw(10) << std::fixed << std::setprecision(2)
               << h1_rate;
      else
        gsInfo << std::setw(10) << "-";

      gsInfo << std::setw(14) << std::scientific << std::setprecision(4)
             << h2err;
      if (ref > 2)
        gsInfo << std::setw(10) << std::fixed << std::setprecision(2)
               << h2_rate;
      else
        gsInfo << std::setw(10) << "-";

      gsInfo << "\n" << std::flush;

      prev_h = h_max;
      prev_l2 = l2err;
      prev_h1 = h1err;
      prev_h2 = h2err;

      // ---- Optional VTK output for solution and original exact function----
      if (plot || !outDir.empty())
      {
        std::string solName =
            prefix + "as_g1_biharmonic_sol_r" + std::to_string(ref);
        std::string exactName =
            prefix + "as_g1_biharmonic_exact_r" + std::to_string(ref);

        gsField<T> solField(mp, sol);
        gsWriteParaview<>(solField, solName, 1000);

        gsField<T> exactField(mp, exact_u, false);
        gsWriteParaview<>(exactField, exactName, 1000);

        gsInfo << " -> Exported VTK plots: " << solName << " & " << exactName
               << "\n";
      }
  }

  gsInfo << std::string(99, '-') << "\n";
  gsInfo << "Done.\n";
  return 0;
}
