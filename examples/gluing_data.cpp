/** @file gluing_data_v4.cpp

    @brief Compute AS-G1 gluing data for multi-patch geometries (v4).

    ## v4 changes relative to v3

    The gluing-data **computation** is mathematically unchanged from v3:
    we solve exactly the same linear least-squares system for
    `(alpha_k, beta_k)`.  v4 only differs in:

      * Extended documentation explaining the tangent-frame convention
        and what consumers (Argyris embeddings) must do when the
        interface orientation flag `flipped` is true.
      * The per-interface output now prints the recommended
        `tangentSign` value that consumers should pass to their
        `createGluingDataArgyrisBasis(..., tangentSign)` call.

    ## What this computes

    For every interface of a 2D multi-patch geometry, finds linear
    functions alpha_k(t), beta_k(t) (k=1,2) along the shared edge that
    satisfy the AS-G1 conditions

        alpha_1(t) * D2(t)  +  alpha_2(t) * D1(t)  =  0,           (G0)
        beta_1(t)  * D2(t)  -  beta_2(t)  * D1(t)  =  D3(t),       (G1)

    with the normalisation  integral(alpha_1 + alpha_2) dt = 1,
    where

        D1 = s1 * det( dG1/dn , dG1/dt )
        D2 = s2 * det( dG2/dn , dG2/dt )
        D3 = s1*s2 * det( dG1/dn , dG2/dn )

    and s_k = (-1)^{par_k} corrects for the side parameter.

    ## Tangent-frame convention (important!)

    The parameter `t` runs along the **patch-1 natural tangent** mapped
    to [0,1].  Patch 2's samples are aligned to this parameter by
    reversing patch 2's natural tangent whenever
    `flipped = !ifc.dirOrientation(ps1, tDir1)` is true.  Therefore:

      * `(alpha_1, beta_1)` are linear coefficients in patch 1's own
        natural tangent (always the same as the gluing-data tangent).
      * `(alpha_2, beta_2)` are linear coefficients in the **gluing-data
        tangent**, not in patch 2's natural tangent.  When `flipped`
        is true, patch 2's natural tangent runs opposite to the
        gluing-data tangent.

    A consumer that builds the Argyris-type trace embedding
    `createGluingDataArgyrisBasis(tb, side, a0, a1, b0, b1,
                                  eps, tangentSign)` must therefore
    pass `tangentSign = flipped ? -1 : +1` for **both** patches when
    the gluing data was computed in the gluing-data frame as above.

    ## Return convention
    Per interface, 8 numbers in patch-1 tangential parametrisation:
        [ a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1 ]
    where _0/_1 are values at t=0 / t=1.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <gismo.h>

#include <gsModeling/gsAsG1Domain.hpp>


using namespace gismo;



// ====================================================================
// Entry point
// ====================================================================
int main(int argc, char* argv[])
{
    using T = real_t;

    std::string geometry("domain2d/2patch/two_bicubic_patches.xml");
    index_t numGaussPerSpan = 0;
    real_t eps = 1e-8;

    gsCmdLine cmd("Compute AS-G1 gluing data (v4 — adds tangentSign hint).");
    cmd.addString("f", "file",     "Input multi-patch file.",        geometry);
    cmd.addInt   ("n", "numGauss", "Gauss points per span (0=auto).", numGaussPerSpan);
    cmd.addReal  ("e", "eps",      "Residual tolerance.",             eps);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<>& mp = *mpPtr;
    mp.computeTopology();

    gsInfo << "Patches: "    << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n";

    gsMatrix<T> M = computeGluingData<T>(mp, T(eps), numGaussPerSpan);

    static const char* sideNames[4] = {"west", "east", "south", "north"};
    gsInfo << "\n=== Gluing data matrix (" << M.rows() << " x " << M.cols() << ") ===\n";
    for (index_t p = 0; p < M.rows(); ++p) {
        gsInfo << "Patch " << p << ":\n";
        for (index_t s = 0; s < 4; ++s)
            gsInfo << "  " << std::setw(5) << sideNames[s]
                   << ":  alpha=" << M(p,4*s+0) << "*(1-t)+" << M(p,4*s+1) << "*t"
                   << "   beta="  << M(p,4*s+2) << "*(1-t)+" << M(p,4*s+3) << "*t\n";
    }
    return 0;
}
