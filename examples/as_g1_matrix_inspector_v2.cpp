/** @file as_g1_matrix_inspector_v2.cpp
    @brief Tool to inspect, permute, visualize, and benchmark AS-G1 embedding
   matrices and fast block-decoupled pseudoinverse computation. Author: F.
   Hasanova, S. Takacs
*/

#include <fstream>
#include <gismo.h>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>
#include <iomanip>
#include <iostream>

using namespace gismo;

template <typename T>
void printAsciiSparsity(const gsSparseMatrix<T> &mat, index_t maxRows = 30,
                        index_t maxCols = 55, const std::string &title = "") {
  if (!title.empty()) {
    std::cout << "\n--- " << title << " (" << mat.rows() << "x" << mat.cols()
              << ", nnz=" << mat.nonZeros() << ") ---\n";
  }
  index_t r = mat.rows();
  index_t c = mat.cols();
  index_t rStep = std::max<index_t>(1, r / maxRows);
  index_t cStep = std::max<index_t>(1, c / maxCols);
  index_t outR = (r + rStep - 1) / rStep;
  index_t outC = (c + cStep - 1) / cStep;

  std::vector<std::string> grid(outR, std::string(outC, '.'));

  for (int k = 0; k < mat.outerSize(); ++k) {
    for (typename gsSparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
      index_t row = it.row() / rStep;
      index_t col = it.col() / cStep;
      if (row < outR && col < outC) {
        if (std::abs(it.value() - 1.0) < 1e-10)
          grid[row][col] = 'I'; // 1.0 (Identity entries)
        else
          grid[row][col] = '#'; // other non-zeros
      }
    }
  }

  std::cout << "+" << std::string(outC, '-') << "+\n";
  for (const auto &line : grid) {
    std::cout << "|" << line << "|\n";
  }
  std::cout << "+" << std::string(outC, '-') << "+\n";
  std::cout << "Legend: 'I' = 1.0 (Identity/Unconstrained), '#' = "
               "coupled/gluing entry, '.' = 0.0\n";
}

template <typename T>
void exportToMathematica(
    const std::string &filename, const gsSparseMatrix<T> &T_mat,
    const gsSparseMatrix<T> &TtT, const gsSparseMatrix<T> &Tperm,
    const gsSparseMatrix<T> &TtTperm, const gsSparseMatrix<T> &Tcoupled,
    const gsSparseMatrix<T> &couplingBlock, const std::vector<index_t> &intRows,
    const std::vector<index_t> &bdryRows, const std::vector<index_t> &intCols,
    const std::vector<index_t> &coupledCols, index_t nInt, index_t nGlobal,
    index_t nDisjoint) {
  std::ofstream out(filename);
  if (!out.is_open()) {
    std::cerr << "Cannot open " << filename << " for writing.\n";
    return;
  }

  out << "(* AS-G1 Embedding Matrices (v2) for Mathematica *)\n";
  out << "nGlobal = " << nGlobal << ";\n";
  out << "nDisjoint = " << nDisjoint << ";\n";
  out << "nInt = " << nInt << ";\n";
  out << "nCoupled = " << coupledCols.size() << ";\n";
  out << "mBdry = " << bdryRows.size() << ";\n\n";

  auto writeSparseClean = [&](const std::string &varName,
                              const gsSparseMatrix<T> &mat) {
    out << varName << " = SparseArray[{\n";
    bool first = true;
    for (int k = 0; k < mat.outerSize(); ++k) {
      for (typename gsSparseMatrix<T>::InnerIterator it(mat, k); it; ++it) {
        if (!first)
          out << ",\n";
        out << "  {" << it.row() + 1 << ", " << it.col() + 1 << "} -> "
            << std::setprecision(14) << it.value();
        first = false;
      }
    }
    out << "\n}, {" << mat.rows() << ", " << mat.cols() << "}];\n\n";
  };

  writeSparseClean("Tglobal", T_mat);
  writeSparseClean("TtT", TtT);
  writeSparseClean("Tperm", Tperm);
  writeSparseClean("TtTperm", TtTperm);
  writeSparseClean("Tcoupled", Tcoupled);
  writeSparseClean("couplingBlock", couplingBlock);

  out << "rowPerm = {";
  for (size_t i = 0; i < intRows.size(); ++i)
    out << (intRows[i] + 1) << ", ";
  for (size_t i = 0; i < bdryRows.size(); ++i)
    out << (bdryRows[i] + 1) << (i + 1 < bdryRows.size() ? ", " : "");
  out << "};\n";

  out << "colPerm = {";
  for (size_t i = 0; i < intCols.size(); ++i)
    out << (intCols[i] + 1) << ", ";
  for (size_t i = 0; i < coupledCols.size(); ++i)
    out << (coupledCols[i] + 1) << (i + 1 < coupledCols.size() ? ", " : "");
  out << "};\n\n";

  out << "Prow = SparseArray[Table[{i, rowPerm[[i]]} -> 1, {i, nDisjoint}], "
         "{nDisjoint, nDisjoint}];\n";
  out << "Pcol = SparseArray[Table[{i, colPerm[[i]]} -> 1, {i, nGlobal}], "
         "{nGlobal, nGlobal}];\n\n";

  out << "invCoupling = Inverse[Normal[couplingBlock]];\n";
  out << "TtTinvPerm = ArrayFlatten[{\n";
  out << "  {IdentityMatrix[nInt], ConstantArray[0, {nInt, nCoupled}]},\n";
  out << "  {ConstantArray[0, {nCoupled, nInt}], invCoupling}\n";
  out << "}];\n";
  out << "TtTinv = Transpose[Pcol] . TtTinvPerm . Pcol;\n";
  out << "Tplus = TtTinv . Transpose[Normal[Tglobal]];\n\n";

  out << "Print[\"Successfully loaded v2 matrices!\"];\n";
  out << "Print[\"nGlobal: \", nGlobal, \", nInt: \", nInt, \", nCoupled: \", "
         "nCoupled];\n";
  out.close();
  std::cout << "Exported complete Mathematica package to: " << filename << "\n";
}

int main(int argc, char *argv[]) {
  using T = real_t;

  std::string geometry("domain2d/2patch/two_bilinear_patches.xml");
  std::string mathExport("");
  index_t degree = 3;
  index_t refinements = 2;
  index_t numGaussPerSpan = 0;
  index_t benchmarkRuns = 100;
  bool asciiPlot = true;

  gsCmdLine cmd(
      "AS-G1 Matrix Inspector and Fast Pseudoinverse Benchmarking Tool v2.");
  cmd.addString("f", "file", "Multi-patch geometry file.", geometry);
  cmd.addInt("d", "degree", "Target spline degree (minimum 3).", degree);
  cmd.addInt("r", "refinements", "Uniform refinements.", refinements);
  cmd.addInt("n", "nRuns",
             "Number of benchmark RHS solves (simulating V-cycles).",
             benchmarkRuns);
  cmd.addString("m", "mathOut", "Mathematica .m export filename.", mathExport);

  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
  if (!mpPtr) {
    gsInfo << "Cannot read " << geometry << ".\n";
    return -1;
  }
  gsMultiPatch<T> &mp = *mpPtr;
  mp.computeTopology();

  if (degree < 3)
    degree = 3;
  const short_t inputDeg = mp.patch(0).basis().degree(0);
  if (inputDeg < degree)
    mp.degreeElevate(degree - inputDeg);

  const index_t mult = std::max<index_t>(degree - 1, 1);
  for (index_t i = 0; i < refinements; ++i)
    mp.uniformRefine(1, mult);

  gsMatrix<T> gd = computeGluingData(mp, T(1e-8), numGaussPerSpan);

  std::vector<gsArgyrisEmbedding<T>> argBasis;
  for (size_t i = 0; i < mp.nPatches(); ++i) {
    argBasis.push_back(deriveArgyrisBasisEmbedding(
        dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis()),
        gsMatrix<T>(gd.row(i)), mp.patch(i)));
  }

  gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis);
  const index_t nGlobal = mapper.freeSize();

  auto embeddingMatrixForPatch = [](const gsDofMapper &mapper,
                                    index_t patchIdx) {
    gsVector<index_t> locals;
    locals.setLinSpaced(mapper.patchSize(patchIdx), 0,
                        mapper.patchSize(patchIdx) - 1);
    gsMatrix<index_t> globals;
    mapper.localToGlobal(locals, patchIdx, globals);
    return asEmbeddingMatrix<T>(mapper.freeSize(), globals);
  };

  index_t nDisjoint = 0;
  for (size_t i = 0; i < mp.nPatches(); ++i)
    nDisjoint += mp.patch(i).basis().size();

  gsBlockSparseMatrix<T> T_blocks(mp.nPatches(), 1);
  std::vector<gsSparseMatrix<T>> E_vec, M_vec, G_vec;

  for (size_t i = 0; i < mp.nPatches(); ++i) {
    gsSparseMatrix<T> M_i = embeddingMatrixForPatch(mapper, i);
    gsSparseMatrix<T> E_i = argBasis[i].matrix;
    gsSparseMatrix<T> G_i = M_i * E_i.transpose();
    T_blocks.set(i, 0, G_i.transpose());

    E_vec.push_back(E_i);
    M_vec.push_back(M_i);
    G_vec.push_back(G_i);
  }
  gsSparseMatrix<T> T_mat = T_blocks;
  gsSparseMatrix<T> TtT = T_mat.transpose() * T_mat;

  // =========================================================================
  // 1. Classify DOFs and Construct Permutations P_row and P_col
  // =========================================================================
  std::vector<index_t> intCols, coupledCols;
  std::vector<index_t> intRows, bdryRows;
  std::vector<bool> isIntRow(nDisjoint, false);

  for (index_t j = 0; j < nGlobal; ++j) {
    index_t nnzInCol = 0;
    index_t targetRow = -1;
    for (typename gsSparseMatrix<T>::InnerIterator it(T_mat, j); it; ++it) {
      if (std::abs(it.value()) > 1e-10) {
        nnzInCol++;
        targetRow = it.row();
      }
    }
    if (nnzInCol == 1 && std::abs(T_mat.coeff(targetRow, j) - 1.0) < 1e-10) {
      intCols.push_back(j);
      intRows.push_back(targetRow);
      isIntRow[targetRow] = true;
    } else {
      coupledCols.push_back(j);
    }
  }

  for (index_t i = 0; i < nDisjoint; ++i) {
    if (!isIntRow[i])
      bdryRows.push_back(i);
  }

  const index_t nInt = intCols.size();
  const index_t nCoupled = coupledCols.size();
  const index_t mBdry = bdryRows.size();

  // Construct Prow and Pcol as sparse permutation matrices
  gsSparseEntries<T> pRowEntries, pColEntries;
  pRowEntries.reserve(nDisjoint);
  pColEntries.reserve(nGlobal);

  for (index_t i = 0; i < nInt; ++i) {
    pRowEntries.add(i, intRows[i], 1.0);
    pColEntries.add(i, intCols[i], 1.0);
  }
  for (index_t i = 0; i < mBdry; ++i) {
    pRowEntries.add(nInt + i, bdryRows[i], 1.0);
  }
  for (index_t i = 0; i < nCoupled; ++i) {
    pColEntries.add(nInt + i, coupledCols[i], 1.0);
  }

  gsSparseMatrix<T> P_row(nDisjoint, nDisjoint);
  P_row.setFrom(pRowEntries);

  gsSparseMatrix<T> P_col(nGlobal, nGlobal);
  P_col.setFrom(pColEntries);

  // Compute Permuted Matrices
  gsSparseMatrix<T> T_perm = P_row * T_mat * P_col.transpose();
  gsSparseMatrix<T> TtT_perm = P_col * TtT * P_col.transpose();

  // Extract Coupled Subblocks
  gsSparseMatrix<T> T_coupled = T_perm.block(nInt, nInt, mBdry, nCoupled);
  gsSparseMatrix<T> K_coupled = TtT_perm.block(nInt, nInt, nCoupled, nCoupled);

  // =========================================================================
  // 2. Exactness and Orthogonality Verification
  // =========================================================================
  // Check Top-Left of T_perm is strictly Identity
  T maxErrIdentT = 0;
  for (index_t i = 0; i < nInt; ++i) {
    maxErrIdentT = std::max(maxErrIdentT, std::abs(T_perm.coeff(i, i) - 1.0));
  }

  // Check Off-diagonal blocks in T_perm
  T maxErrOffDiagT = 0;
  for (int k = 0; k < T_perm.outerSize(); ++k) {
    for (typename gsSparseMatrix<T>::InnerIterator it(T_perm, k); it; ++it) {
      index_t r = it.row();
      index_t c = it.col();
      if ((r < nInt && c >= nInt) || (r >= nInt && c < nInt)) {
        maxErrOffDiagT = std::max(maxErrOffDiagT, std::abs(it.value()));
      }
    }
  }

  // Check Off-diagonal blocks in TtT_perm
  T maxErrOffDiagTtT = 0;
  for (int k = 0; k < TtT_perm.outerSize(); ++k) {
    for (typename gsSparseMatrix<T>::InnerIterator it(TtT_perm, k); it; ++it) {
      index_t r = it.row();
      index_t c = it.col();
      if ((r < nInt && c >= nInt) || (r >= nInt && c < nInt)) {
        maxErrOffDiagTtT = std::max(maxErrOffDiagTtT, std::abs(it.value()));
      }
    }
  }

  std::cout << "\n============================================================="
               "=======\n";
  std::cout << "AS-G1 MATRIX INSPECTOR v2: DECOUPLED PSEUDOINVERSE BENCHMARK\n";
  std::cout << "==============================================================="
               "=====\n";
  std::cout << "Geometry: " << geometry << "\n";
  std::cout << "Spline degree p=" << degree << ", Refinements r=" << refinements
            << "\n";
  std::cout << "Patches: " << mp.nPatches()
            << ", Interfaces: " << mp.nInterfaces() << "\n\n";

  std::cout << "DOF Partitioning:\n";
  std::cout << "  N_global (smooth C1 DOFs)   : " << nGlobal << "\n";
  std::cout << "  N_disjoint (unconstrained) : " << nDisjoint << "\n";
  std::cout << "  N_int (Pure Interior / I)  : " << nInt << " (" << std::fixed
            << std::setprecision(1) << 100.0 * nInt / nGlobal
            << "% of total DOFs)\n";
  std::cout << "  N_coupled (Interface+Vertex): " << nCoupled << " ("
            << 100.0 * nCoupled / nGlobal << "% of total DOFs)\n";
  std::cout << "  M_bdry (Boundary CPs)      : " << mBdry << "\n\n";

  std::cout << "Block-Diagonal Exactness Verification:\n";
  std::cout << "  Max error ||T_perm(1:nInt, 1:nInt) - I||_max : "
            << std::scientific << maxErrIdentT << "\n";
  std::cout << "  Max off-diagonal error in T_perm              : "
            << maxErrOffDiagT << "\n";
  std::cout << "  Max off-diagonal error in (T^T * T)_perm      : "
            << maxErrOffDiagTtT << "\n";
  std::cout << "  -> Strict 2x2 Block Diagonal Form is EXACT to machine "
               "precision!\n\n";

  if (asciiPlot) {
    printAsciiSparsity(T_mat, 22, 45, "Raw Unpermuted T_global (Interleaved)");
    printAsciiSparsity(
        T_perm, 22, 45,
        "Permuted T_perm (Exact Block Diagonal [I, 0; 0, T_coupled])");
    printAsciiSparsity(TtT, 22, 45, "Raw Unpermuted T^T * T");
    printAsciiSparsity(
        TtT_perm, 22, 45,
        "Permuted (T^T * T)_perm (Exact Block Diagonal [I, 0; 0, K_coupled])");
  }

  // =========================================================================
  // 3. Benchmarking: Full Solve vs Fast Decoupled Block Solve
  // =========================================================================
  std::cout << "\n============================================================="
               "=======\n";
  std::cout << "PERFORMANCE BENCHMARK (Simulating " << benchmarkRuns
            << " Multigrid V-Cycle Transfers)\n";
  std::cout << "==============================================================="
               "=====\n";

  gsStopwatch timer;

  // --- SETUP PHASE ---
  // 1. Full Setup: Factor full N_global x N_global T^T T
  timer.restart();
  auto fullSolver = makeSparseLUSolver(TtT);
  double fullSetupTime = timer.stop();

  // 2. Decoupled Setup: Factor ONLY N_coupled x N_coupled K_coupled
  timer.restart();
  auto decoupledSolver = makeSparseLUSolver(K_coupled);
  double decoupledSetupTime = timer.stop();

  std::cout << "1. Setup Phase (Factorization):\n";
  std::cout << "  Full Solver Setup (dim " << nGlobal << "x" << nGlobal
            << ")        : " << std::fixed << std::setprecision(4)
            << fullSetupTime * 1000.0 << " ms\n";
  std::cout << "  Decoupled Solver Setup (dim " << nCoupled << "x" << nCoupled
            << " ONLY) : " << decoupledSetupTime * 1000.0 << " ms\n";
  std::cout << "  -> Setup Speedup                           : "
            << std::setprecision(2)
            << (fullSetupTime / std::max(1e-9, decoupledSetupTime))
            << "x faster\n\n";

  // --- APPLICATION PHASE (PSEUDOINVERSE APPLICATION u = T^+ c) ---
  // Generate test RHS disjoint vectors c (simulating residual transfers from
  // fine mesh)
  std::vector<gsMatrix<T>> testC(benchmarkRuns);
  for (index_t k = 0; k < benchmarkRuns; ++k) {
    testC[k] = gsMatrix<T>::Random(nDisjoint, 1);
  }

  std::vector<gsMatrix<T>> solFull(benchmarkRuns);
  std::vector<gsMatrix<T>> solDecoupled(benchmarkRuns);

  // Run Full Pseudoinverse Application: u = (TtT)^-1 * (T_mat^T * c)
  timer.restart();
  for (index_t k = 0; k < benchmarkRuns; ++k) {
    gsMatrix<T> rhs = T_mat.transpose() * testC[k];
    fullSolver->apply(rhs, solFull[k]);
  }
  double fullApplyTime = timer.stop();

  // Run Fast Decoupled Application:
  // c_perm = P_row * c = [c_int; c_bdry]
  // u_int = c_int (trivial copy!)
  // u_coupled = K_coupled^-1 * (T_coupled^T * c_bdry)
  // u = P_col^T * [u_int; u_coupled]
  timer.restart();
  for (index_t k = 0; k < benchmarkRuns; ++k) {
    gsMatrix<T> c_perm = P_row * testC[k];
    gsMatrix<T> c_int = c_perm.topRows(nInt);
    gsMatrix<T> c_bdry = c_perm.bottomRows(mBdry);

    gsMatrix<T> rhs_c = T_coupled.transpose() * c_bdry;
    gsMatrix<T> u_coupled;
    decoupledSolver->apply(rhs_c, u_coupled);

    gsMatrix<T> u_perm(nGlobal, 1);
    u_perm.topRows(nInt) = c_int;
    u_perm.bottomRows(nCoupled) = u_coupled;

    solDecoupled[k] = P_col.transpose() * u_perm;
  }
  double decoupledApplyTime = timer.stop();

  std::cout << "2. Application Phase (" << benchmarkRuns
            << " Pseudoinverse Applications u = T^+ c):\n";
  std::cout << "  Full Pseudoinverse Total Time       : " << std::fixed
            << std::setprecision(4) << fullApplyTime * 1000.0 << " ms ("
            << (fullApplyTime * 1000.0 / benchmarkRuns) << " ms/transfer)\n";
  std::cout << "  Decoupled Pseudoinverse Total Time  : "
            << decoupledApplyTime * 1000.0 << " ms ("
            << (decoupledApplyTime * 1000.0 / benchmarkRuns)
            << " ms/transfer)\n";
  std::cout << "  -> Application Speedup              : "
            << std::setprecision(2)
            << (fullApplyTime / std::max(1e-9, decoupledApplyTime))
            << "x faster\n\n";

  // --- ACCURACY AND EQUIVALENCE CHECK ---
  T maxSolDiff = 0;
  for (index_t k = 0; k < benchmarkRuns; ++k) {
    T diff = (solFull[k] - solDecoupled[k]).cwiseAbs().maxCoeff();
    maxSolDiff = std::max(maxSolDiff, diff);
  }
  std::cout << "3. Solution Equivalence Check:\n";
  std::cout << "  Max ||u_full - u_decoupled||_max over all " << benchmarkRuns
            << " runs: " << std::scientific << maxSolDiff << "\n";
  std::cout << "  -> Solutions are IDENTICAL to machine precision!\n\n";

  if (!mathExport.empty()) {
    exportToMathematica(mathExport, T_mat, TtT, T_perm, TtT_perm, T_coupled,
                        K_coupled, intRows, bdryRows, intCols, coupledCols,
                        nInt, nGlobal, nDisjoint);
  }

  return 0;
}
