/** @file as_g1_two_patch_basis.cpp

    @brief Build a C1 conforming basis across two patches over AS-G1 geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <algorithm>
#include <gismo.h>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

using namespace gismo;

// ====================================================================
// main – two-patch AS-G1 basis
// ====================================================================

int main(int argc, char *argv[]) {
  using T = real_t;

  std::string geometry("domain2d/2patch/two_bilinear_patches.xml");
  std::string outDir("");
  index_t numGaussPerSpan = 0;
  index_t refinements = 2;
  index_t plot = -1;

  gsCmdLine cmd("AS-G1 conforming two-patch basis.");
  cmd.addString("f", "file", "Multi-patch geometry file.", geometry);
  cmd.addString("o", "outDir",
                "Output directory for pvd/vts files (created if needed).",
                outDir);
  cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).",
             numGaussPerSpan);
  cmd.addInt("r", "refinements", "Uniform refinements before proceeding.",
             refinements);
  cmd.addInt("p", "plot",
             "Basis function to plot: -1=none, -2=ALL, 0..N=single.", plot);

  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  // ---- Create output directory if specified ----
  std::string prefix;
  if (!outDir.empty()) {
    // Ensure trailing slash
    if (outDir.back() != '/')
      outDir += '/';
    // Create the directory (works on Linux/macOS)
    gsFileManager::mkdir(outDir);
    prefix = outDir;
    gsInfo << "Output directory: " << outDir << "\n";
  }

  // ---- Read geometry ----
  gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
  if (!mpPtr) {
    gsInfo << "Cannot read " << geometry << ".\n";
    return -1;
  }
  gsMultiPatch<T> &mp = *mpPtr;
  mp.computeTopology();

  GISMO_ENSURE(mp.nPatches() == 2,
               "This example is designed for exactly two patches.");
  GISMO_ENSURE(mp.nInterfaces() == 1,
               "Expected exactly one interface between the two patches.");

  gsInfo << "Patches: " << mp.nPatches() << "  Interfaces: " << mp.nInterfaces()
         << "\n";

  // ---- Ensure minimum degree 3 ----
  const short_t inputDeg = mp.patch(0).basis().degree(0);
  if (inputDeg < 3) {
    const short_t elev = 3 - inputDeg;
    mp.degreeElevate(elev);
    gsInfo << "Input degree " << inputDeg << " → elevated by " << elev
           << " to degree 3\n";
  } else {
    gsInfo << "Input degree " << inputDeg << " (no elevation needed)\n";
  }

  // ---- Refine ----
  // Insert knots with multiplicity = degree - 1 (gives C^1 continuity),
  // which is the minimum needed for the Argyris construction
  // (sideBasis must have minInteriorMultiplicity > 1).
  if (refinements > 0) {
    const short_t deg = mp.patch(0).basis().degree(0);
    const index_t mult = std::max<index_t>(deg - 1, 1);
    for (index_t i = 0; i < refinements; ++i)
      mp.uniformRefine(1, mult);
    gsInfo << "Refined " << refinements << " time(s), knot mult = " << mult
           << " (degree " << deg << ")\n";
  }

  // ---- Identify the interface ----
  const boundaryInterface &ifc = *mp.iBegin();
  const patchSide ps1 = ifc.first();
  const patchSide ps2 = ifc.second();

  gsInfo << "Interface: patch " << ps1.patch << " side " << ps1.side()
         << " <-> patch " << ps2.patch << " side " << ps2.side() << "\n";

  // ---- Compute gluing data for the interface ----
  gsMatrix<T> gd =
      computeGluingDataForInterface(mp, ifc, T(1e-8), numGaussPerSpan)
          .transpose();

  gsInfo << "\nGluing data:\n"
         << "  Patch " << ps1.patch << " side " << ps1.side()
         << ": alpha=" << gd(0, 0) << "*(1-t)+" << gd(0, 1) << "*t"
         << "  beta=" << gd(0, 2) << "*(1-t)+" << gd(0, 3) << "*t\n"
         << "  Patch " << ps2.patch << " side " << ps2.side()
         << ": alpha=" << gd(0, 4) << "*(1-t)+" << gd(0, 5) << "*t"
         << "  beta=" << gd(0, 6) << "*(1-t)+" << gd(0, 7) << "*t\n\n";

  // ---- Build per-patch interface-side embeddings ----
  const gsTensorBSplineBasis<2, T> &tb1 =
      dynamic_cast<const gsTensorBSplineBasis<2, T> &>(
          mp.patch(ps1.patch).basis());
  const gsTensorBSplineBasis<2, T> &tb2 =
      dynamic_cast<const gsTensorBSplineBasis<2, T> &>(
          mp.patch(ps2.patch).basis());

  gsAsG1Embedding<T> argBasis1 =
      deriveArgyrisBasisEmbedding<T>(tb1, ps1.side(), gd.leftCols(4), T(1e-12));
  gsAsG1Embedding<T> argBasis2 = deriveArgyrisBasisEmbedding<T>(
      tb2, ps2.side(), gd.rightCols(4), T(1e-12));

  const gsSparseMatrix<T> &E1 = argBasis1.matrix;
  const gsSparseMatrix<T> &E2 = argBasis2.matrix;

  gsInfo << "Patch " << ps1.patch << " interface embedding: " << E1.rows()
         << " x " << E1.cols() << "\n";
  gsInfo << "Patch " << ps2.patch << " interface embedding: " << E2.rows()
         << " x " << E2.cols() << "\n";

  // ====================================================================
  // Setup of gsDofMapper
  // ====================================================================

  GISMO_ENSURE(
      argBasis1.sizes[1] == argBasis2.sizes[1],
      "Level 0 (function values) interfaces sizes differ across interface ("
          << argBasis1.sizes[1] << " vs " << argBasis1.sizes[1] << ").");
  GISMO_ENSURE(
      argBasis1.sizes[2] == argBasis2.sizes[2],
      "Level 1 (crossing derivatives) interfaces sizes differ across interface "
      "(" << argBasis1.sizes[2]
          << " vs " << argBasis2.sizes[2] << ").");

  const index_t nInt1 = argBasis1.sizes[0];
  const index_t nInt2 = argBasis2.sizes[0];
  const index_t nLvl0 = argBasis1.sizes[1];
  const index_t nLvl1 = argBasis1.sizes[2];

  gsDofMapper mapper;
  {
    gsVector<index_t> patchDofSizes(2);
    patchDofSizes[0] = argBasis1.matrix.cols();
    patchDofSizes[1] = argBasis2.matrix.cols();
    mapper = gsDofMapper(patchDofSizes);

    // Check interface orientation: if the tangential directions
    // run in opposite directions, we need to reverse the DOF mapping
    // for patch 2's shared columns.
    const short_t tanDir1 = 1 - ps1.direction();
    const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

    gsInfo << "Interface orientation: " << (flipped ? "FLIPPED" : "aligned")
           << "\n";

    for (index_t j1 = 0; j1 < nLvl0; ++j1) {
      // If flipped, DOF j1 on patch 1 corresponds to DOF (nLvl0-1-j1) on
      // patch 2.
      const index_t j2 = flipped ? nLvl0 - 1 - j1 : j1;
      gsInfo << "mapper.matchDof(" << ps1.patch << ", " << nInt1 + j1 << ", "
             << ps2.patch << ", " << nInt2 + j2 << ")\n";
      mapper.matchDof(ps1.patch, nInt1 + j1, ps2.patch, nInt2 + j2);
    }
    for (index_t j1 = 0; j1 < nLvl1; ++j1) {
      // If flipped, DOF j1 on patch 1 corresponds to DOF (nLvl1-1-j1) on
      // patch 2.
      const index_t j2 = flipped ? nLvl1 - 1 - j1 : j1;
      gsInfo << "mapper.matchDof(" << ps1.patch << ", " << nInt1 + nLvl0 + j1
             << ", " << ps2.patch << ", " << nInt2 + nLvl0 + j2 << ")\n";
      mapper.matchDof(ps1.patch, nInt1 + nLvl0 + j1, ps2.patch,
                      nInt2 + nLvl0 + j2);
    }
    mapper.finalize();
  }

  const index_t nGlobal = mapper.freeSize();
  gsInfo << "\nGlobal DOFs: " << nGlobal << " = " << nInt1 << " (int1) + "
         << nInt2 << " (int2) + " << nLvl0
         << " (level0: shared function values) + " << nLvl1
         << " (level1: shared crossing derivatives)\n";

  GISMO_ASSERT(nGlobal == nInt1 + nInt2 + nLvl0 + nLvl1, "Size missmatch.");

  gsSparseMatrix<T> cc = deriveCornerEmbedding(tb1, mp.patch(ps1.patch),
                                               gsMatrix<T>(gd.leftCols(4)));
  gsInfo << "Corner embedding:\n" << cc << "\n";

  // ====================================================================
  // Extract global embedding matrices from gsDofMapper
  // ====================================================================

  auto embeddingMatrixForPatch = [](const gsDofMapper &mapper,
                                    index_t patchIdx) {
    gsVector<index_t> locals;
    locals.setLinSpaced(mapper.patchSize(patchIdx), 0,
                        mapper.patchSize(patchIdx) - 1);
    gsMatrix<index_t> globals;
    mapper.localToGlobal(locals, patchIdx, globals);
    return asEmbeddingMatrix<T>(mapper.freeSize(), globals);
  };

  gsSparseMatrix<T> G1 =
      E1 * embeddingMatrixForPatch(mapper, ps1.patch).transpose();
  gsSparseMatrix<T> G2 =
      E2 * embeddingMatrixForPatch(mapper, ps2.patch).transpose();

  gsInfo << "\nGlobal-to-patch matrices:\n"
         << "  G1: " << G1.rows() << " x " << G1.cols() << "\n"
         << "  G2: " << G2.rows() << " x " << G2.cols() << "\n";

  // ====================================================================
  // Numerical G1 smoothness verification
  // ====================================================================
  {
    const short_t normDir1 = ps1.direction(), tanDir1_ = 1 - normDir1;
    const short_t normDir2 = ps2.direction(), tanDir2_ = 1 - normDir2;
    const bool par1_ = ps1.parameter(), par2_ = ps2.parameter();
    const bool ifcFlipped = !ifc.dirOrientation(ps1, tanDir1_);

    gsMatrix<T> sup1 = mp.patch(ps1.patch).support();
    gsMatrix<T> sup2 = mp.patch(ps2.patch).support();
    const T ifcCoord1 = sup1(normDir1, par1_ ? 1 : 0);
    const T ifcCoord2 = sup2(normDir2, par2_ ? 1 : 0);
    const T t1a = sup1(tanDir1_, 0), t1b = sup1(tanDir1_, 1);
    const T t2a = sup2(tanDir2_, 0), t2b = sup2(tanDir2_, 1);

    const index_t nCheck = 21;
    T maxValErr = 0, maxGradErr = 0;

    T maxErrInt = 0, maxErrTrace = 0, maxErrL2 = 0;

    for (index_t idx = 0; idx < nGlobal; ++idx) {
      gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
      globalVec(idx) = T(1);
      gsVector<T> c1 = G1 * globalVec;
      gsVector<T> c2 = G2 * globalVec;

      auto func1 = tb1.makeGeometry(c1);
      auto func2 = tb2.makeGeometry(c2);

      T thisMaxGrad = 0;

      for (index_t i = 0; i < nCheck; ++i) {
        T s = T(i) / T(nCheck - 1);
        T t1 = t1a + s * (t1b - t1a);
        T s2 = ifcFlipped ? (1.0 - s) : s;
        T t2 = t2a + s2 * (t2b - t2a);

        gsMatrix<T> pt1(2, 1), pt2(2, 1);
        pt1(normDir1, 0) = ifcCoord1;
        pt1(tanDir1_, 0) = t1;
        pt2(normDir2, 0) = ifcCoord2;
        pt2(tanDir2_, 0) = t2;

        // Check C0: values must match
        gsMatrix<T> v1 = func1->eval(pt1);
        gsMatrix<T> v2 = func2->eval(pt2);
        T valErr = std::abs(v1(0, 0) - v2(0, 0));
        maxValErr = std::max(maxValErr, valErr);

        // Check G1: physical gradients must match
        gsMatrix<T> df1, df2, dG1, dG2;
        func1->deriv_into(pt1, df1);
        func2->deriv_into(pt2, df2);
        mp.patch(ps1.patch).deriv_into(pt1, dG1);
        mp.patch(ps2.patch).deriv_into(pt2, dG2);

        gsMatrix<T> J1(2, 2), J2(2, 2);
        J1(0, 0) = dG1(0, 0);
        J1(0, 1) = dG1(1, 0);
        J1(1, 0) = dG1(2, 0);
        J1(1, 1) = dG1(3, 0);
        J2(0, 0) = dG2(0, 0);
        J2(0, 1) = dG2(1, 0);
        J2(1, 0) = dG2(2, 0);
        J2(1, 1) = dG2(3, 0);

        gsVector<T> paramGrad1(2), paramGrad2(2);
        paramGrad1(0) = df1(0, 0);
        paramGrad1(1) = df1(1, 0);
        paramGrad2(0) = df2(0, 0);
        paramGrad2(1) = df2(1, 0);

        gsVector<T> physGrad1 = J1.inverse().transpose() * paramGrad1;
        gsVector<T> physGrad2 = J2.inverse().transpose() * paramGrad2;

        T gradErr = (physGrad1 - physGrad2).norm();
        maxGradErr = std::max(maxGradErr, gradErr);
        thisMaxGrad = std::max(thisMaxGrad, gradErr);
      }

      // Classify DOF
      if (idx < nInt1 + nInt2)
        maxErrInt = std::max(maxErrInt, thisMaxGrad);
      else if (idx < nInt1 + nInt2 + nLvl0)
        maxErrTrace = std::max(maxErrTrace, thisMaxGrad);
      else
        maxErrL2 = std::max(maxErrL2, thisMaxGrad);
    }

    gsInfo << "\n=== G1 Smoothness Check (physical gradient) ===\n"
           << "  Max C0 (value) error:      " << maxValErr << "\n"
           << "  Max grad (physical) error:  " << maxGradErr << "\n"
           << "    Interior DOFs:    " << maxErrInt << "\n"
           << "    Trace DOFs:       " << maxErrTrace << "\n"
           << "    D-deriv DOFs:     " << maxErrL2 << "\n";
    if (maxValErr < 1e-8 && maxGradErr < 1e-3)
      gsInfo << "  STATUS: PASS\n";
    else
      gsInfo << "  STATUS: FAIL\n";
  }

  // ====================================================================
  // Plot global basis functions on both patches
  // ====================================================================

  // Build a name and category tag for each global DOF
  // Global layout: [p0_interior | p1_interior | ifc_trace | ifc_dderiv]
  auto basisName = [&](index_t idx) -> std::string {
    if (idx < nInt1)
      return "p" + std::to_string(ps1.patch) + "_int_" + std::to_string(idx);
    else if (idx < nInt1 + nInt2)
      return "p" + std::to_string(ps2.patch) + "_int_" +
             std::to_string(idx - nInt1);
    else if (idx < nInt1 + nInt2 + nLvl0)
      return "ifc_trace_" + std::to_string(idx - nInt1 - nInt2);
    else
      return "ifc_dderiv_" + std::to_string(idx - nInt1 - nInt2 - nLvl0);
  };

  // Determine which basis functions to plot
  std::vector<index_t> toPlot;
  if (plot == -2) {
    toPlot.resize(nGlobal);
    std::iota(toPlot.begin(), toPlot.end(), 0);
    gsInfo << "\nPlotting all " << nGlobal << " basis functions ...\n";
  } else if (plot >= 0 && plot < nGlobal) {
    toPlot.push_back(plot);
    gsInfo << "\nPlotting global basis function " << plot << " ("
           << basisName(plot) << ") ...\n";
  } else if (plot >= nGlobal) {
    gsInfo << "\nBasis function index " << plot << " out of range [0, "
           << nGlobal - 1 << "].\n";
  }

  for (const index_t idx : toPlot) {
    gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
    globalVec(idx) = T(1);

    gsVector<T> c1 = G1 * globalVec;
    gsVector<T> c2 = G2 * globalVec;

    gsMatrix<T> coefs1 = c1;
    gsMatrix<T> coefs2 = c2;

    gsMultiPatch<T> sol;
    sol.addPatch(tb1.makeGeometry(give(coefs1)));
    sol.addPatch(tb2.makeGeometry(give(coefs2)));

    gsField<T> field(mp, sol);

    std::string fname = prefix + "as_g1_" + basisName(idx);
    gsWriteParaview<>(field, fname, 1000);
    gsInfo << "  [" << idx << "] " << fname << ".pvd\n";
  }

  if (!toPlot.empty()) {
    gsInfo << "\nFile naming convention:\n"
           << "  p<N>_int_<K>   = patch N interior DOF K\n"
           << "  ifc_trace_<K>  = shared interface trace DOF K\n"
           << "  ifc_dderiv_<K> = shared interface d-derivative DOF K\n";
  }

  // ====================================================================
  // Print summary information about the global basis
  // ====================================================================
  gsInfo << "\n=== Summary ===\n";
  gsInfo << "  Tensor basis size per patch: " << tb1.size() << "\n";
  gsInfo << "  Per-patch embedding cols:    " << E1.cols() << " (int=" << nInt1
         << " nLvl0=" << nLvl0 << " nLvl1=" << nLvl1 << ")\n";
  gsInfo << "  Global DOFs:                 " << nGlobal << "\n";
  gsInfo << "    Patch " << ps1.patch << " interior: " << nInt1 << "\n";
  gsInfo << "    Patch " << ps2.patch << " interior: " << nInt2 << "\n";
  gsInfo << "    Shared interface trace:     " << nLvl0 << "\n";
  gsInfo << "    Shared interface d-deriv:   " << nLvl1 << "\n";
  gsInfo << "\nTo plot basis function k:\n"
         << "  ./bin/as_g1_two_patch_basis -f <file> -r <ref> -p k\n"
         << "To plot ALL basis functions:\n"
         << "  ./bin/as_g1_two_patch_basis -f <file> -r <ref> -p -2\n";

  gsInfo << "\nDone.\n";
  return 0;
}
