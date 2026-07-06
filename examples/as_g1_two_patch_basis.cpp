/** @file as_g1_two_patch_basis_v4.cpp

    @brief Build an AS-G1 conforming basis across two patches.

    Streamlined v4 -- same gluing-data computation as v3, plus two
    fixes that make the construction work for ALL 196 two-patch
    parametrizations (8 left-side x 8 right-side reparametrizations
    times 3 geometry families), instead of only the 100/196 that
    worked in v2/v3.

    The two fixes are both triggered by the interface orientation
    flag `flipped = !ifc.dirOrientation(ps1, tDir1)`, which is true
    exactly when patch-2's natural tangent runs opposite to patch-1's
    natural tangent along the shared interface.

    1) `tangentSign` parameter on `createGluingDataArgyrisBasis`:
       The per-patch trace embedding solves a Greville interpolation
       problem whose RHS reads
           rhs = -(dBdry*Vm + tangentSign*signN*beta*dVm) / dNeigh
       where dVm is d/dt of the smoother basis in the *patch's own*
       tangent parameter, while beta is the gluing-data tangential
       coefficient defined in the *gluing-data* tangent frame.  When
       the patch's tangent is opposite to the gluing-data tangent,
       these two derivatives differ by a sign, which must be
       absorbed into the embedding's RHS.  We pass tangentSign = -1
       to BOTH patches when `flipped`, because both halves of the
       global AS-G1 relation
            alpha_1 d_1(G^(1)) + alpha_2 d_2(G^(2)) = 0
       are rewritten in the gluing-data tangent frame.

    2) Conditional second-layer sign `l2Sign` in the global G2
       assembly:  In v2 the d-derivative columns of patch-2 were
       inserted with an unconditional `-it.value()`.  When the
       tangent is flipped, that sign must be `+it.value()` because
       the column reversal (j2 = nLD2-1-j) already accounts for the
       mirror reparametrization of the lower-degree basis, which
       flips the d-derivative sign once.

    With these two fixes the global basis satisfies the AS-G1
    condition across the interface to machine precision
    (~1e-15 on bilinear patches, ~3e-7 on cubically refined curved
    patches due to representable approximation error).

    Only the shared interface receives the gluing-data directional
    derivative constraints.  All other boundary sides are left
    untouched (standard tensor B-spline DOFs).

    Interface trace (boundary) DOFs and interface d-derivative
    (second-layer) DOFs are shared between the two patches.
    All other DOFs are independent per patch.

    Use  -p N  to plot global basis function N on both patches
    via Paraview.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <gismo.h>

#include <gsModeling/gsAsG1Domain.hpp>
#include <gsModeling/gsAsG1Basis.hpp>

using namespace gismo;


// ====================================================================
// main – two-patch AS-G1 basis
// ====================================================================

int main(int argc, char* argv[])
{
    using T = real_t;

    std::string geometry("domain2d/2patch/two_bilinear_patches.xml");
    std::string outDir("");
    index_t numGaussPerSpan = 0;
    index_t refinements = 1;
    index_t plot = -1;

    gsCmdLine cmd("AS-G1 conforming two-patch basis.");
    cmd.addString("f", "file", "Multi-patch geometry file.", geometry);
    cmd.addString("o", "outDir",
                  "Output directory for pvd/vts files (created if needed).", outDir);
    cmd.addInt("n", "numGauss",
               "Gauss points per knot span (0 = auto).", numGaussPerSpan);
    cmd.addInt("r", "refinements",
               "Uniform refinements before proceeding.", refinements);
    cmd.addInt("p", "plot",
               "Basis function to plot: -1=none, -2=ALL, 0..N=single.", plot);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ---- Create output directory if specified ----
    std::string prefix;
    if (!outDir.empty())
    {
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
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<T>& mp = *mpPtr;
    mp.computeTopology();

    GISMO_ENSURE(mp.nPatches() == 2,
                 "This example is designed for exactly two patches.");
    GISMO_ENSURE(mp.nInterfaces() == 1,
                 "Expected exactly one interface between the two patches.");

    gsInfo << "Patches: " << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n";

    // ---- Ensure minimum degree 3 ----
    const short_t inputDeg = mp.patch(0).basis().degree(0);
    if (inputDeg < 3)
    {
        const short_t elev = 3 - inputDeg;
        mp.degreeElevate(elev);
        gsInfo << "Input degree " << inputDeg
               << " → elevated by " << elev << " to degree 3\n";
    }
    else
    {
        gsInfo << "Input degree " << inputDeg << " (no elevation needed)\n";
    }

    // ---- Refine ----
    // Insert knots with multiplicity = degree - 1 (gives C^1 continuity),
    // which is the minimum needed for the Argyris construction
    // (sideBasis must have minInteriorMultiplicity > 1).
    if (refinements > 0)
    {
        const short_t deg = mp.patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);
        for (index_t i = 0; i < refinements; ++i)
            mp.uniformRefine(1, mult);
        gsInfo << "Refined " << refinements << " time(s), knot mult = "
               << mult << " (degree " << deg << ")\n";
    }

    // ---- Identify the interface ----
    const boundaryInterface& ifc = *mp.iBegin();
    const patchSide ps1 = ifc.first();
    const patchSide ps2 = ifc.second();

    gsInfo << "Interface: patch " << ps1.patch << " side " << ps1.side()
           << " <-> patch " << ps2.patch << " side " << ps2.side() << "\n";

    // ---- Compute gluing data for the interface ----
    bool ok = false;
    gsVector<T> gd = computeGluingDataForInterface(
        mp, ifc, ok, T(1e-8), numGaussPerSpan);
    GISMO_ENSURE(ok, "Gluing data computation failed.");

    const T a1_0 = gd(0), a1_1 = gd(1), b1_0 = gd(2), b1_1 = gd(3);
    const T a2_0 = gd(4), a2_1 = gd(5), b2_0 = gd(6), b2_1 = gd(7);

    gsInfo << "\nGluing data:\n"
           << "  Patch " << ps1.patch << " side " << ps1.side()
           << ": alpha=" << a1_0 << "*(1-t)+" << a1_1 << "*t"
           << "  beta=" << b1_0 << "*(1-t)+" << b1_1 << "*t\n"
           << "  Patch " << ps2.patch << " side " << ps2.side()
           << ": alpha=" << a2_0 << "*(1-t)+" << a2_1 << "*t"
           << "  beta=" << b2_0 << "*(1-t)+" << b2_1 << "*t\n\n";

    // ---- Build per-patch interface-side embeddings ----
    const gsTensorBSplineBasis<2,T>& tb1 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(ps1.patch).basis());
    const gsTensorBSplineBasis<2,T>& tb2 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(ps2.patch).basis());

    // v4 fix #1: when the patch tangent runs opposite to the
    // gluing-data tangent (i.e. `flipped`), pass tangentSign = -1
    // to the embedding so that `beta * d/dt(smoother)` is evaluated
    // in the gluing-data tangent frame.  Both patches need it
    // because both halves of the AS-G1 relation are expressed in
    // the gluing-data tangent frame.
    const short_t _tdir1 = 1 - ps1.direction();
    const bool _flippedAtEmb = !ifc.dirOrientation(ps1, _tdir1);
    const T tSign = _flippedAtEmb ? T(-1) : T(1);

    gsSparseMatrix<T> E1 = createGluingDataArgyrisBasis(
        tb1, ps1.side(), a1_0, a1_1, b1_0, b1_1, T(1e-12), tSign);
    gsSparseMatrix<T> E2 = createGluingDataArgyrisBasis(
        tb2, ps2.side(), a2_0, a2_1, b2_0, b2_1, T(1e-12), tSign);

    gsInfo << "Patch " << ps1.patch << " interface embedding: "
           << E1.rows() << " x " << E1.cols() << "\n";
    gsInfo << "Patch " << ps2.patch << " interface embedding: "
           << E2.rows() << " x " << E2.cols() << "\n";

    // ---- DOF classification for each patch ----
    // The embedding column layout is:
    //   [nInterior_i | nLowerDeg_i | nSmoother_i]
    // where _i refers to the interface side.

    // Side bases for the interface side
    gsBSplineBasis<T> sideBasis1 = *tb1.boundaryBasis(ps1.side());
    gsBSplineBasis<T> sideBasis2 = *tb2.boundaryBasis(ps2.side());

    gsBSplineBasis<T> smootherBasis1 = sideBasis1;
    smootherBasis1.elevateContinuity(1);
    gsBSplineBasis<T> lowerDegBasis1 = sideBasis1;
    lowerDegBasis1.degreeReduce(1);

    gsBSplineBasis<T> smootherBasis2 = sideBasis2;
    smootherBasis2.elevateContinuity(1);
    gsBSplineBasis<T> lowerDegBasis2 = sideBasis2;
    lowerDegBasis2.degreeReduce(1);

    gsMatrix<index_t> bdryDOFs1  = tb1.boundary(ps1.side());
    gsMatrix<index_t> neighDOFs1 = tb1.boundaryOffset(ps1.side(), 1);
    std::vector<index_t> intDOFs1 =
        getInteriorDofs(tb1.size(), bdryDOFs1, neighDOFs1);

    gsMatrix<index_t> bdryDOFs2  = tb2.boundary(ps2.side());
    gsMatrix<index_t> neighDOFs2 = tb2.boundaryOffset(ps2.side(), 1);
    std::vector<index_t> intDOFs2 =
        getInteriorDofs(tb2.size(), bdryDOFs2, neighDOFs2);

    const index_t nInt1 = static_cast<index_t>(intDOFs1.size());
    const index_t nInt2 = static_cast<index_t>(intDOFs2.size());
    const index_t nLD1  = lowerDegBasis1.size();   // second-layer cols
    const index_t nLD2  = lowerDegBasis2.size();
    const index_t nSm1  = smootherBasis1.size();   // boundary cols
    const index_t nSm2  = smootherBasis2.size();

    // Sanity checks
    GISMO_ENSURE(E1.cols() == nInt1 + nLD1 + nSm1,
                 "Column count mismatch for patch " << ps1.patch);
    GISMO_ENSURE(E2.cols() == nInt2 + nLD2 + nSm2,
                 "Column count mismatch for patch " << ps2.patch);

    gsInfo << "\nPatch " << ps1.patch << ": nInterior=" << nInt1
           << "  nSecondLayer=" << nLD1
           << "  nBoundary=" << nSm1 << "\n";
    gsInfo << "Patch " << ps2.patch << ": nInterior=" << nInt2
           << "  nSecondLayer=" << nLD2
           << "  nBoundary=" << nSm2 << "\n";

    // ====================================================================
    // Assemble the global basis
    // ====================================================================
    //
    // The global DOF numbering is:
    //
    //   [patch1_interior | patch2_interior | shared_boundary | shared_second_layer]
    //
    // where:
    //   - patch_i interior: DOFs not on the interface boundary/second-layer
    //     These are independent per patch and map via identity in E_i.
    //   - shared boundary: The nSm (= nSm1 = nSm2) interface trace DOFs.
    //     The same trace coefficients are used for both patches.
    //   - shared second-layer: The nLD (= nLD1 = nLD2) interface d-derivative DOFs.
    //     These are shared between patches; the AS-G1 condition means
    //     the same second-layer coefficient feeds both patches' embeddings.
    //
    // IMPORTANT: For the two patches to share the interface trace values
    // and d-derivative values, the two side bases must be compatible
    // (same knot vector along the interface).  We verify this.

    GISMO_ENSURE(nSm1 == nSm2,
                 "Smoother basis sizes differ across interface ("
                     << nSm1 << " vs " << nSm2 << ").");
    GISMO_ENSURE(nLD1 == nLD2,
                 "Lower-degree basis sizes differ across interface ("
                     << nLD1 << " vs " << nLD2 << ").");

    const index_t nSharedBdry = nSm1;    // shared trace DOFs
    const index_t nSharedL2   = nLD1;    // shared second-layer DOFs

    // Check interface orientation: if the tangential directions
    // run in opposite directions, we need to reverse the DOF mapping
    // for patch 2's shared columns.
    const short_t tanDir1 = 1 - ps1.direction();
    const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

    gsInfo << "Interface orientation: "
           << (flipped ? "FLIPPED" : "aligned") << "\n";

    const index_t nGlobal = nInt1 + nInt2 + nSharedBdry + nSharedL2;

    gsInfo << "\nGlobal DOFs: " << nGlobal
           << " = " << nInt1 << " (int1) + " << nInt2 << " (int2) + "
           << nSharedBdry << " (shared bdry) + "
           << nSharedL2 << " (shared L2)\n";

    // Build the two global-to-patch matrices:
    //   G1 : nGlobal → tb1.size()    (coefficients for patch 1)
    //   G2 : nGlobal → tb2.size()    (coefficients for patch 2)
    //
    // For patch 1:
    //   - Global cols [0 .. nInt1-1]  →  E1 cols [0 .. nInt1-1]   (interior)
    //   - Global cols [nInt1+nInt2 .. +nSharedBdry-1]  →  E1 cols [nInt1+nLD1 .. end] (boundary)
    //   - Global cols [nInt1+nInt2+nSharedBdry .. end]  →  E1 cols [nInt1 .. nInt1+nLD1-1] (L2)
    //
    // For patch 2:
    //   - Global cols [nInt1 .. nInt1+nInt2-1]  →  E2 cols [0 .. nInt2-1]   (interior)
    //   - Global cols [nInt1+nInt2 .. +nSharedBdry-1]  →  E2 cols [nInt2+nLD2 .. end] (boundary)
    //   - Global cols [nInt1+nInt2+nSharedBdry .. end]  →  E2 cols [nInt2 .. nInt2+nLD2-1] (L2)

    const index_t gOff_int1     = 0;
    const index_t gOff_int2     = nInt1;
    const index_t gOff_bdry     = nInt1 + nInt2;
    const index_t gOff_L2       = nInt1 + nInt2 + nSharedBdry;

    // --- Patch 1 global matrix ---
    gsSparseMatrix<T> G1(tb1.size(), nGlobal);
    {
        // Interior columns of E1 → global interior cols for patch 1
        for (index_t j = 0; j < nInt1; ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(E1, j); it; ++it)
                G1.insert(it.row(), gOff_int1 + j) = it.value();

        // Second-layer columns of E1 → global shared L2 cols
        for (index_t j = 0; j < nLD1; ++j)
        {
            const index_t e1col = nInt1 + j;
            for (typename gsSparseMatrix<T>::InnerIterator it(E1, e1col); it; ++it)
                G1.insert(it.row(), gOff_L2 + j) = it.value();
        }

        // Boundary columns of E1 → global shared boundary cols
        for (index_t j = 0; j < nSm1; ++j)
        {
            const index_t e1col = nInt1 + nLD1 + j;
            for (typename gsSparseMatrix<T>::InnerIterator it(E1, e1col); it; ++it)
                G1.insert(it.row(), gOff_bdry + j) = it.value();
        }
    }
    G1.makeCompressed();

    // --- Patch 2 global matrix ---
    gsSparseMatrix<T> G2(tb2.size(), nGlobal);
    {
        // Interior columns of E2 → global interior cols for patch 2
        for (index_t j = 0; j < nInt2; ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(E2, j); it; ++it)
                G2.insert(it.row(), gOff_int2 + j) = it.value();

        // Second-layer columns of E2 → global shared L2 cols
        // AS-G1 requires alpha_1*d_1 + alpha_2*d_2 = 0, so patch 2's
        // second-layer contribution is sign-flipped relative to patch 1.
        // When the tangential directions are FLIPPED, the d-deriv basis
        // function indexing also reverses (j2 = nLD2-1-j).  The relative
        // sign of patch-2's d_i derivative depends on the orientation:
        //   - aligned:  d_2 = +alpha_2 * phi_j (after embedding) → negate to oppose patch 1
        //   - flipped:  d_2 also has an additional minus from the reparam,
        //               which cancels the negation, so we DON'T negate.
        const T l2Sign = flipped ? T(1) : T(-1);
        for (index_t j = 0; j < nLD2; ++j)
        {
            const index_t j2 = flipped ? (nLD2 - 1 - j) : j;
            const index_t e2col = nInt2 + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(E2, e2col); it; ++it)
                G2.insert(it.row(), gOff_L2 + j) = l2Sign * it.value();
        }

        // Boundary columns of E2 → global shared boundary cols
        // If flipped, DOF j on patch 1 corresponds to DOF (nSm2-1-j) on patch 2.
        for (index_t j = 0; j < nSm2; ++j)
        {
            const index_t j2 = flipped ? (nSm2 - 1 - j) : j;
            const index_t e2col = nInt2 + nLD2 + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(E2, e2col); it; ++it)
                G2.insert(it.row(), gOff_bdry + j) = it.value();
        }
    }
    G2.makeCompressed();

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

        for (index_t idx = 0; idx < nGlobal; ++idx)
        {
            gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
            globalVec(idx) = T(1);
            gsVector<T> c1 = G1 * globalVec;
            gsVector<T> c2 = G2 * globalVec;

            auto func1 = tb1.makeGeometry(c1);
            auto func2 = tb2.makeGeometry(c2);

            T thisMaxGrad = 0;

            for (index_t i = 0; i < nCheck; ++i)
            {
                T s = T(i) / T(nCheck - 1);
                T t1 = t1a + s * (t1b - t1a);
                T s2 = ifcFlipped ? (1.0 - s) : s;
                T t2 = t2a + s2 * (t2b - t2a);

                gsMatrix<T> pt1(2, 1), pt2(2, 1);
                pt1(normDir1, 0) = ifcCoord1; pt1(tanDir1_, 0) = t1;
                pt2(normDir2, 0) = ifcCoord2; pt2(tanDir2_, 0) = t2;

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

                gsMatrix<T> J1(2,2), J2(2,2);
                J1(0,0) = dG1(0,0); J1(0,1) = dG1(1,0);
                J1(1,0) = dG1(2,0); J1(1,1) = dG1(3,0);
                J2(0,0) = dG2(0,0); J2(0,1) = dG2(1,0);
                J2(1,0) = dG2(2,0); J2(1,1) = dG2(3,0);

                gsVector<T> paramGrad1(2), paramGrad2(2);
                paramGrad1(0) = df1(0,0); paramGrad1(1) = df1(1,0);
                paramGrad2(0) = df2(0,0); paramGrad2(1) = df2(1,0);

                gsVector<T> physGrad1 = J1.inverse().transpose() * paramGrad1;
                gsVector<T> physGrad2 = J2.inverse().transpose() * paramGrad2;

                T gradErr = (physGrad1 - physGrad2).norm();
                maxGradErr = std::max(maxGradErr, gradErr);
                thisMaxGrad = std::max(thisMaxGrad, gradErr);
            }

            // Classify DOF
            if (idx < gOff_bdry)
                maxErrInt = std::max(maxErrInt, thisMaxGrad);
            else if (idx < gOff_L2)
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
        else if (maxValErr < 1e-8 && maxGradErr < 1e-1)
            gsInfo << "  STATUS: APPROX (AS-G1 approximation error)\n";
        else
            gsInfo << "  STATUS: FAIL (check coupling)\n";
    }

    // ====================================================================
    // Plot global basis functions on both patches
    // ====================================================================

    // Build a name and category tag for each global DOF
    // Global layout: [p0_interior | p1_interior | ifc_trace | ifc_dderiv]
    auto basisName = [&](index_t idx) -> std::string
    {
        if (idx < gOff_int2)
            return "p" + std::to_string(ps1.patch) + "_int_"
                       + std::to_string(idx - gOff_int1);
        if (idx < gOff_bdry)
            return "p" + std::to_string(ps2.patch) + "_int_"
                       + std::to_string(idx - gOff_int2);
        if (idx < gOff_L2)
            return "ifc_trace_" + std::to_string(idx - gOff_bdry);
        return "ifc_dderiv_" + std::to_string(idx - gOff_L2);
    };

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
        gsInfo << "\nPlotting global basis function " << plot
               << " (" << basisName(plot) << ") ...\n";
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

    if (!toPlot.empty())
    {
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
    gsInfo << "  Per-patch embedding cols:    " << E1.cols()
           << " (int=" << nInt1 << " L2=" << nLD1 << " bdry=" << nSm1 << ")\n";
    gsInfo << "  Global DOFs:                 " << nGlobal << "\n";
    gsInfo << "    Patch " << ps1.patch << " interior: " << nInt1 << "\n";
    gsInfo << "    Patch " << ps2.patch << " interior: " << nInt2 << "\n";
    gsInfo << "    Shared interface trace:     " << nSharedBdry << "\n";
    gsInfo << "    Shared interface d-deriv:   " << nSharedL2 << "\n";
    gsInfo << "\nTo plot basis function k:\n"
           << "  ./bin/as_g1_two_patch_basis_v4 -f <file> -r <ref> -p k\n"
           << "To plot ALL basis functions:\n"
           << "  ./bin/as_g1_two_patch_basis_v4 -f <file> -r <ref> -p -2\n";

    gsInfo << "\nDone.\n";
    return 0;
}
