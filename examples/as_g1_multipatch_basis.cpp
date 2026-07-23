/** @file as_g1_multipatch_basis.cpp

    @brief Build a C1 conforming basis across multiple patches over AS-G1
   geometry.

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

  std::string geometry(
      "domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
  std::string outDir("");
  index_t degree = 3;
  index_t numGaussPerSpan = 0;
  index_t refinements = 2;
  index_t plot = -1;

  gsCmdLine cmd("AS-G1 conforming multipatch basis.");
  cmd.addString("f", "file", "Multi-patch geometry file.", geometry);
  cmd.addString("o", "outDir",
                "Output directory for pvd/vts files (created if needed).",
                outDir);
  cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
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

  gsInfo << "Patches: " << mp.nPatches() << "  Interfaces: " << mp.nInterfaces()
         << "\n";

  // ---- Ensure minimum degree ----
  if (degree < 3) degree = 3;
  const short_t inputDeg = mp.patch(0).basis().degree(0);
  if (inputDeg < degree) {
    const short_t elev = degree - inputDeg;
    mp.degreeElevate(elev);
    gsInfo << "Input degree " << inputDeg << " → elevated by " << elev
           << " to degree " << degree << "\n";
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
  gsMatrix<T> gd = computeGluingData(mp, T(1e-8), numGaussPerSpan);

  gsInfo << "\nGluing data:\n" << gd << "\n\n";

  // ---- Build per-patch interface-side embeddings ----
  std::vector<gsArgyrisEmbedding<T>> argBasis;
  for (size_t i = 0; i < mp.nPatches(); ++i) {
    argBasis.push_back(deriveArgyrisBasisEmbedding(
        dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mp.patch(i).basis()),
        gsMatrix<T>(gd.row(i)), mp.patch(i)));
  }

  for (size_t i = 0; i < mp.nPatches(); ++i) {
    gsInfo << "Patch " << i << " interface embedding:\n"
           << argBasis[i].matrix << "\n\n";
  }

  // ====================================================================
  // Setup of gsDofMapper
  // ====================================================================
<<<<<<< HEAD

  /*const index_t nInt1 = argBasis1.sizes[0];
  const index_t nInt2 = argBasis2.sizes[0];
  const index_t nLvl0 = argBasis1.sizes[1];
  const index_t nLvl1 = argBasis1.sizes[2];*/

  auto sumUntil = [](const gsVector<index_t, 13> &vec, index_t until) {
    index_t sum = 0;
    for (index_t i = 0; i < until; ++i)
      sum += vec(i);
    return sum;
  };

  gsDofMapper mapper(patchDofSizes);
  for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
    const boundaryInterface &ifc = *it;
    const patchSide ps1 = ifc.first();
    const patchSide ps2 = ifc.second();

    const index_t nLvl0 = argBasis[ps1.patch].sizes[1 + 2 * (ps1.m_index - 1)];
    const index_t offLvl0_1 =
        sumUntil(argBasis[ps1.patch].sizes, 1 + 2 * (ps1.m_index - 1));
    const index_t offLvl0_2 =
        sumUntil(argBasis[ps2.patch].sizes, 1 + 2 * (ps2.m_index - 1));

    const index_t nLvl1 = argBasis[ps1.patch].sizes[2 + 2 * (ps1.m_index - 1)];
    const index_t offLvl1_1 =
        sumUntil(argBasis[ps1.patch].sizes, 2 + 2 * (ps1.m_index - 1));
    const index_t offLvl1_2 =
        sumUntil(argBasis[ps2.patch].sizes, 2 + 2 * (ps2.m_index - 1));

    GISMO_ASSERT(nLvl0 == argBasis[ps2.patch].sizes[1 + 2 * (ps2.m_index - 1)],
                 "Dimension missmatch.");
    GISMO_ASSERT(nLvl1 == argBasis[ps2.patch].sizes[2 + 2 * (ps2.m_index - 1)],
                 "Dimension missmatch.");

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
      gsInfo << "mapper.matchDof(" << ps1.patch << ", " << offLvl0_1 + j1
             << ", " << ps2.patch << ", " << offLvl0_2 + j2 << ")\n";
      mapper.matchDof(ps1.patch, offLvl0_1 + j1, ps2.patch, offLvl0_2 + j2);
    }
    for (index_t j1 = 0; j1 < nLvl1; ++j1) {
      // If flipped, DOF j1 on patch 1 corresponds to DOF (nLvl1-1-j1) on
      // patch 2.
      const index_t j2 = flipped ? nLvl1 - 1 - j1 : j1;
      gsInfo << "mapper.matchDof(" << ps1.patch << ", " << offLvl1_1 + j1
             << ", " << ps2.patch << ", " << offLvl1_2 + j2 << ")\n";
      mapper.matchDof(ps1.patch, offLvl1_1 + j1, ps2.patch, offLvl1_2 + j2);
    }

    // Also match the corners
    std::vector<patchCorner> corners;
    ps1.getContainedCorners(2, corners);
    GISMO_ASSERT(corners.size() == 2, "Unexpected number of corners");
    for (index_t i = 0; i < 2; ++i) {
      const index_t c1 = corners[i].m_index - 1;
      const index_t c2 = ifc.mapCorner(corners[i]).m_index - 1;

      const index_t off_corner_1 = sumUntil(argBasis[ps1.patch].sizes, 9 + c1);
      const index_t off_corner_2 = sumUntil(argBasis[ps2.patch].sizes, 9 + c2);

      for (index_t i = 0; i < 6; ++i)
        mapper.matchDof(ps1.patch, off_corner_1 + i, ps2.patch,
                        off_corner_2 + i);
    }
  }
  mapper.finalize();
=======
  gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis);
>>>>>>> b20f4a27035e18cb3256a8084589e752abd35372

  const index_t nGlobal = mapper.freeSize();
  gsInfo << "\nGlobal DOFs: " << nGlobal << "\n";

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

  std::vector<gsSparseMatrix<T>> G;
  for (size_t i = 0; i < mp.nPatches(); ++i) {
    gsInfo << "Patch " << i << " embedding matrix:\n"
           << argBasis[i].matrix.rows() << "x" << argBasis[i].matrix.cols()
           << "\nembeddingMatrixForPatch: "
           << embeddingMatrixForPatch(mapper, i).rows() << "x"
           << embeddingMatrixForPatch(mapper, i).cols() << "\n\n";
    G.push_back(embeddingMatrixForPatch(mapper, i) *
                argBasis[i].matrix.transpose());
  }

  // ====================================================================
  // Numerical G1 smoothness verification
  // ====================================================================
  for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
    const boundaryInterface &ifc = *it;
    const patchSide ps1 = ifc.first();
    const patchSide ps2 = ifc.second();
    const gsSparseMatrix<T> &G1 = G[ps1.patch].transpose();
    const gsSparseMatrix<T> &G2 = G[ps2.patch].transpose();
    const gsTensorBSplineBasis<2, T> &tb1 =
        dynamic_cast<const gsTensorBSplineBasis<2, T> &>(
            mp.patch(ps1.patch).basis());
    const gsTensorBSplineBasis<2, T> &tb2 =
        dynamic_cast<const gsTensorBSplineBasis<2, T> &>(
            mp.patch(ps2.patch).basis());

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
    }

    gsInfo << "\n=== G1 Smoothness Check (physical gradient) ===\n"
           << "  Max C0 (value) error:      " << maxValErr << "\n"
           << "  Max grad (physical) error:  " << maxGradErr << "\n";
    if (maxValErr < 1e-8 && maxGradErr < 1e-3)
      gsInfo << "  STATUS: PASS\n";
    else
      gsInfo << "  STATUS: FAIL\n";
  }

  // ====================================================================
  // Plot global basis functions on both patches
  // ====================================================================

  // Determine which basis functions to plot
    std::vector<index_t> toPlot;
    if (plot == -2)
    {
        toPlot.resize(nGlobal);
        std::iota(toPlot.begin(), toPlot.end(), 0);
        gsInfo << "\nPlotting all " << nGlobal << " basis functions ...\n";
    }
    else if (plot >= 0 && plot < nGlobal)
    {
        toPlot.push_back(plot);
        gsInfo << "\nPlotting global basis function " << plot << " ...\n";
    }
    else if (plot >= nGlobal)
    {
        gsInfo << "\nBasis function index " << plot
               << " out of range [0, " << nGlobal - 1 << "].\n";
    }

    for (const index_t idx : toPlot)
    {
        gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
        globalVec(idx) = T(1);

        gsMultiPatch<T> sol;
        for (std::size_t k=0; k<mp.nPatches(); ++k)
            sol.addPatch(mp.patch(k).basis().makeGeometry(give(G[k].transpose() * globalVec)));

        gsField<T> field(mp, sol);

        std::string fname = outDir + "as_g1_" + std::to_string(idx);
        gsWriteParaview<>(field, fname, 1000);
        gsInfo << "  [" << idx << "] " << fname << ".pvd\n";
    }

    // ====================================================================
    // Print summary information about the global basis
    // ====================================================================
    gsInfo << "\n=== Summary ===\n";
    gsInfo << "\nTo plot basis function k:\n"
          << "  ./bin/as_g1_multipatch_basis -f <file> -r <ref> -p k\n"
          << "To plot ALL basis functions:\n"
          << "  ./bin/as_g1_multipatch_basis -f <file> -r <ref> -p -2\n";

  gsInfo << "\nDone.\n";
  return 0;
}
