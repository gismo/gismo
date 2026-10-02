/** @file immersed_tetmesh_clip_quadrature_example.cpp

    @brief Stand-alone 3D cut-cell quadrature for a tetrahedral-mesh domain
    against a uniform background grid, by exact tet-cell clipping (no
    Gauss-Green identity, no level set). This is the direct 3D analogue of
    the mesh-clipping path of immersed_gauss_green_quadrature_example.cpp
    (optional/gsOpenCascade).

    The clipping machinery -- the ASCII gmsh MSH 4.1 reader, the exact
    six-split tet clip, the Sutherland-Hodgman boundary triangle clip, the
    collapsed-Gauss volume/boundary rules, `tetClipQuadrature` and the
    unclipped-mesh oracles -- lives in gsTetMeshClip.h (namespace
    gsTetClip), whose file header documents the theory, the exact-predicate
    conventions, the degree argument, the reader restriction and the traps.
    This driver checks that quadrature (T1-T9, Tx) against those oracles and
    against analytic cube moments.

    Algoim baseline (`--algoim`, default off, INFO only, never affects
    ALL CHECKS PASS/FAIL or the exit code):
    - builds a `gsSurfMesh<real_t>` from the driver's own derived,
      outward-oriented boundary triangles (`M.bdrTri`), in the mesh's own
      coordinates (no unit-box normalization -- the per-cell comparison
      against the exact clip rule below is only meaningful if both rules
      see the same geometry in the same coordinates), and wraps it in a
      `gsMeshSignedDist` over the background box;
    - runs `gsAlgoimAdaptiveRule` -- the rule behind quadrature rule 15 --
      with the option set rule 15 forwards (`src/gsAssembler/gsQuadrature.h`,
      `quA=1.0`, `quB=1` from the assembler defaults, `maxDepth` from
      `--maxDepth`, default 0, everything else at
      `gsAlgoimAdaptiveRule::defaultOptions()`) on every background cell of
      the same `Grid3` used by the clip rule, once for the volume
      (`{phi<0}`) and once for the surface (`{phi==0}`);
    - prints A0-A5 INFO lines: level-set build/sign-check, total volume,
      volume moments a+b+c<=2p against the unclipped-mesh oracle of T2, the worst 5
      per-cell |V_clip-V_algoim|/h^3, Algoim surface area against the exact
      mesh boundary area, and the timings of both methods.

    Trap: when a boundary triangle's plane exactly coincides with a
    background grid plane (`--case cube` at the default `--n0`), phi is
    identically zero on a whole cell face there. Algoim's height-function
    construction (called from `gsAlgoimAdaptiveRule::mapTo`) can fail to
    terminate on that
    degeneracy (see the split-plane-on-the-interface guard,
    `optional/gsAlgoim/gsAlgoimAdaptiveRule.h:452-457`), exhausting memory.
    `algoimBaselineUsable` detects the axis-aligned-face case and skips the
    Algoim run instead of hanging; `--case cube --n0 3` at `-r 0` or `-r 1`
    does not trigger it (h = 2/3 or 1/3: the cube faces at +-0.5 are not grid
    lines) and runs Algoim normally; from `-r 2` on they are grid lines and it
    is skipped. The detector is specific to that one coincidence, not to
    every way `mapTo` can exhaust memory, so wrap every `--algoim`
    invocation under a memory ulimit regardless, e.g.
    `( ulimit -v 8000000; ./immersed_tetmesh_clip_quadrature_example
    --case rotcube -r 1 --algoim )`, with a `timeout` besides.

    `--plot` (default off, independent of `--algoim`) writes the clip
    volume/boundary quadrature nodes, and (with `--algoim`) the Algoim
    volume/surface nodes, as ParaView point sets in the current working
    directory. Trap: at the default p=2 every clip sub-tet carries 392
    nodes, so a rotcube run at `-r 1` (33 260 pieces) already produces about
    13 million points and a VTP of several hundred MB. Nothing here
    subsamples; use a small `-k`/`-r` (e.g. `-k 1 -r 0`) for a plotting run.

    Example command lines:
      ./immersed_tetmesh_clip_quadrature_example
      ./immersed_tetmesh_clip_quadrature_example --case cube -r 2
      ./immersed_tetmesh_clip_quadrature_example --case sphere
      ./immersed_tetmesh_clip_quadrature_example --case rotcube --fit
      ./immersed_tetmesh_clip_quadrature_example --case none -f volumes/tetmesh_sphere.msh
      ( ulimit -v 8000000; ./immersed_tetmesh_clip_quadrature_example --case rotcube -r 1 --algoim )
      ./immersed_tetmesh_clip_quadrature_example --case rotcube -k 1 -r 0 --plot

    Reference:
    - M.J. Powell, T. Abel, "An exact general remeshing scheme applied to
      physically conservative voxelization", J. Comput. Phys. 297 (2015)
      340-356.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include <gsDomain/gsMeshLevelSet.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "gsTetMeshClip.h"

using namespace gismo;
using namespace gsTetClip;

namespace {

//----------------------------------------------------------------------------
// Analytic moments for the cube cases.
//----------------------------------------------------------------------------

real_t analyticCubeAlignedMoment(int a, int b, int c, real_t lo, real_t hi)
{
    auto term = [](int k, real_t lo_, real_t hi_)
    { return (std::pow(hi_,k+1) - std::pow(lo_,k+1)) / (real_t)(k+1); };
    return term(a,lo,hi) * term(b,lo,hi) * term(c,lo,hi);
}

/// int_{[-s/2,s/2]^3} (c+R xi)^(a,b,c) dxi by tensor Gauss, p+1 points per
/// direction: the integrand's total degree in xi is a+b+c <= 2p <= 2p+1, so
/// the rule is exact.
real_t analyticCubeRotatedMoment(int a, int b, int c, real_t s, const Vec3 & cvec,
                                 const real_t R[3][3], index_t p)
{
    const Rule1D Rg = gauss01(p+1);
    const index_t m = p+1;
    std::vector<real_t> xi(m), wxi(m);
    for (index_t k = 0; k != m; ++k) { xi[k] = s*(Rg.u[k]-0.5); wxi[k] = s*Rg.w[k]; }

    KahanSum total;
    for (index_t i = 0; i != m; ++i)
    for (index_t j = 0; j != m; ++j)
    for (index_t k = 0; k != m; ++k)
    {
        const real_t x1 = xi[i], x2 = xi[j], x3 = xi[k];
        const real_t x = cvec[0] + R[0][0]*x1 + R[0][1]*x2 + R[0][2]*x3;
        const real_t y = cvec[1] + R[1][0]*x1 + R[1][1]*x2 + R[1][2]*x3;
        const real_t z = cvec[2] + R[2][0]*x1 + R[2][1]*x2 + R[2][2]*x3;
        const real_t wgt = wxi[i]*wxi[j]*wxi[k];
        total.add(wgt*std::pow(x,a)*std::pow(y,b)*std::pow(z,c));
    }
    return total.value();
}

//----------------------------------------------------------------------------
// Test-case constants.
//----------------------------------------------------------------------------

const real_t ROTCUBE_SIDE = 0.8;
const real_t ROTCUBE_THETA_Z_DEG = 17.0;
const real_t ROTCUBE_THETA_X_DEG = 29.0;
const Vec3   ROTCUBE_TRANSLATE = {0.03, -0.02, 0.01};

const real_t CUBE_HALF = 0.5;

const real_t SPHERE_RADIUS = 0.55;

void rotcubeMatrix(real_t R[3][3])
{
    const real_t tz = ROTCUBE_THETA_Z_DEG * (real_t)EIGEN_PI / 180.0;
    const real_t tx = ROTCUBE_THETA_X_DEG * (real_t)EIGEN_PI / 180.0;
    const real_t cz = math::cos(tz), sz = math::sin(tz);
    const real_t cx = math::cos(tx), sx = math::sin(tx);
    real_t Rz[3][3] = { {cz,-sz,0}, {sz,cz,0}, {0,0,1} };
    real_t Rx[3][3] = { {1,0,0}, {0,cx,-sx}, {0,sx,cx} };
    for (int i = 0; i != 3; ++i)
        for (int j = 0; j != 3; ++j)
        {
            real_t s = 0;
            for (int kk = 0; kk != 3; ++kk) s += Rx[i][kk]*Rz[kk][j];
            R[i][j] = s;
        }
}

//----------------------------------------------------------------------------
// Checks (T1-T9 plus the high-degree check Tx).
//----------------------------------------------------------------------------

bool checkT1(const TetClipStats & stats, real_t L)
{
    const real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t worst = 0; bool pass = true;
    for (size_t t = 0; t != stats.tetVol.size(); ++t)
    {
        const real_t diff = math::abs(stats.tetClippedVol[t] - stats.tetVol[t]);
        const real_t bound = 1e-14*stats.tetVol[t] + 16*eps*L*stats.tetDiam[t]*stats.tetDiam[t];
        if (diff > bound) pass = false;
        worst = math::max(worst, diff/stats.tetVol[t]);
    }
    gsInfo << "T1 per-tet partition: " << fmtSci(worst) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT2(const TetClipQuadrature & Q, const TetMesh & M, index_t p, real_t V_mesh,
            std::vector<real_t> & oracleMoments, const std::vector<std::array<int,3> > & moments)
{
    real_t worst = 0;
    oracleMoments.resize(moments.size());
    for (size_t k = 0; k != moments.size(); ++k)
    {
        const std::array<int,3> & m = moments[k];
        oracleMoments[k] = tetVolumeOracleMoment(M, m[0],m[1],m[2], p);
        worst = math::max(worst, scaledErr(volMoment(Q,m[0],m[1],m[2]), oracleMoments[k], V_mesh));
    }
    const bool pass = worst <= 1e-13;
    gsInfo << "T2 volume moments vs unclipped oracle: " << fmtSci(worst) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT3(const TetClipQuadrature & Q, const std::string & caseName, index_t p,
            const std::vector<std::array<int,3> > & moments, const std::vector<real_t> & oracleMoments,
            real_t s, const Vec3 & cvec, const real_t R[3][3])
{
    if ("cube" != caseName && "rotcube" != caseName)
    {
        gsInfo << "T3 volume moments vs analytic: skipped (case=" << caseName << ")\n";
        return true;
    }
    const bool aligned = ("cube" == caseName);
    const real_t V_exact = aligned ? 1.0 : s*s*s;

    real_t worst = 0, worstOracle = 0;
    for (size_t k = 0; k != moments.size(); ++k)
    {
        const std::array<int,3> & m = moments[k];
        const real_t analytic = aligned ? analyticCubeAlignedMoment(m[0],m[1],m[2], -CUBE_HALF, CUBE_HALF)
                                        : analyticCubeRotatedMoment(m[0],m[1],m[2], s, cvec, R, p);
        worst = math::max(worst, scaledErr(volMoment(Q,m[0],m[1],m[2]), analytic, V_exact));
        worstOracle = math::max(worstOracle, scaledErr(oracleMoments[k], analytic, V_exact));
    }
    const bool pass = worst <= 1e-13;
    gsInfo << "T3 volume moments vs analytic: " << fmtSci(worst) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    gsInfo << "T3 info: oracle vs analytic scaled diff=" << fmtSci(worstOracle) << "  INFO\n";
    return pass;
}

bool checkT4(const TetClipQuadrature & Q, real_t L)
{
    const real_t eps = std::numeric_limits<real_t>::epsilon();
    const real_t hh = Q.grid.h*Q.grid.h;
    const real_t box = Q.grid.h*hh;
    const real_t tolAbs = 16*eps*L*hh;
    real_t cutMin = std::numeric_limits<real_t>::infinity();
    real_t cutMax = -std::numeric_limits<real_t>::infinity();
    bool pass = true;
    for (size_t id = 0; id != Q.vol.size(); ++id)
    {
        const real_t vol = cellVolume(Q.vol[id]);
        if (vol < -(1e-14*box+tolAbs) || vol > box*(1.0+1e-14)+tolAbs) pass = false;
        if (Cut == Q.status[id])
        {
            cutMin = math::min(cutMin, vol/box);
            cutMax = math::max(cutMax, vol/box);
        }
    }
    gsInfo << "T4 cell bounds: cut-cell min " << fmtSci(cutMin) << " max " << fmtSci(cutMax)
          << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT5(const TetClipQuadrature & Q, real_t V_mesh)
{
    KahanSum sumVol;
    for (const CellRule3 & r : Q.vol) sumVol.add(cellVolume(r));
    const real_t err = math::abs(sumVol.value() - V_mesh) / V_mesh;
    const bool pass = err <= 1e-13;
    gsInfo << "T5 per-cell sum: " << fmtSci(err) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT6(const TetClipQuadrature & Q, real_t P_mesh)
{
    const real_t clippedArea = bdrMoment(Q,0,0,0,-1);
    const real_t err = scaledErr(clippedArea, P_mesh, P_mesh);
    const bool pass = err <= 1e-13;
    gsInfo << "T6 boundary area: " << fmtSci(err) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT7(const TetClipQuadrature & Q)
{
    real_t worst = 0;
    for (int d = 0; d != 3; ++d) worst = math::max(worst, math::abs(bdrMoment(Q,0,0,0,d)));
    const bool pass = worst <= 1e-13;
    gsInfo << "T7 oint n dS: " << fmtSci(worst) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT8(const TetClipQuadrature & Q)
{
    const real_t V_Q = volMoment(Q,0,0,0);
    const real_t err = math::abs(bdrFluxXdotN(Q) - 3.0*V_Q) / math::max((real_t)1.0, V_Q);
    const bool pass = err <= 1e-13;
    gsInfo << "T8 oint x.n dS = 3V: " << fmtSci(err) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

bool checkT9(const TetMesh & M, const Grid3 & grid, const std::string & caseName,
            real_t s, const Vec3 & cvec, const real_t R[3][3])
{
    if ("rotcube" == caseName)
    {
        real_t worst = 0;
        for (int cx = -1; cx <= 1; cx += 2)
        for (int cy = -1; cy <= 1; cy += 2)
        for (int cz = -1; cz <= 1; cz += 2)
        {
            const Vec3 xi = { cx*0.5*s, cy*0.5*s, cz*0.5*s };
            const Vec3 target = { cvec[0] + R[0][0]*xi[0]+R[0][1]*xi[1]+R[0][2]*xi[2],
                                  cvec[1] + R[1][0]*xi[0]+R[1][1]*xi[1]+R[1][2]*xi[2],
                                  cvec[2] + R[2][0]*xi[0]+R[2][1]*xi[1]+R[2][2]*xi[2] };
            real_t best = std::numeric_limits<real_t>::infinity();
            for (const Vec3 & P : M.P)
                best = math::min(best, math::sqrt(dot3(sub3(P,target), sub3(P,target))));
            worst = math::max(worst, best/s);
        }
        const bool pass = worst <= 1e-14;
        gsInfo << "T9 geometry check (rotcube corners): " << fmtSci(worst) << "  " << (pass ? "PASS" : "FAIL") << "\n";
        return pass;
    }
    if ("cube" == caseName)
    {
        std::vector<real_t> X, Y, Z;
        gridLines(grid, X, Y, Z);
        std::set<index_t> nodeIdx;
        for (const std::array<index_t,4> & tri : M.bdrTri)
            for (int v = 0; v != 3; ++v) nodeIdx.insert(tri[v]);
        index_t nOn = 0;
        for (index_t idx : nodeIdx)
        {
            const Vec3 & P = M.P[idx];
            bool on = false;
            for (index_t i = 1; i != grid.n && !on; ++i) if (P[0] == X[i]) on = true;
            for (index_t j = 1; j != grid.n && !on; ++j) if (P[1] == Y[j]) on = true;
            for (index_t k = 1; k != grid.n && !on; ++k) if (P[2] == Z[k]) on = true;
            if (on) ++nOn;
        }
        const index_t nB = (index_t)nodeIdx.size();
        const bool pass = (nOn == nB) && (nB > 0);
        gsInfo << "T9 geometry check (cube on-knot nodes): N_on=" << nOn << " N_b=" << nB
              << "  " << (pass ? "PASS" : "FAIL") << "\n";
        return pass;
    }
    gsInfo << "T9 geometry check: skipped (case=" << caseName << ")\n";
    return true;
}

/// High-degree check whose exponents (a=b=c=2p, total degree 6p) sit
/// outside the T2/T3 moment set, whose total degree never exceeds 2p:
/// exercises the rule degree beyond what the moment set reaches. On the
/// shipped meshes the volume half detects a 4p volume rule (scaled error
/// about 5e-12 on rotcube/cube, about 1.03e-12 on the sphere), while the
/// boundary half detects a 2p boundary rule but not a 4p one: the clipped
/// boundary pieces are small enough that the 4p-vs-6p error falls below
/// 1e-12.
///
/// Scaled by the moment's OWN magnitude, never by V_mesh/P_mesh: this
/// specific moment is several orders of magnitude smaller than V_mesh
/// (e.g. about 2.8e7x on the shipped rotcube mesh), so a V_mesh-scaled
/// relative error would mask a real under-integration by that same factor,
/// making the check unable to fail regardless of tolerance. The volume
/// moment is always positive (the integrand is a product of even powers),
/// so its own oracle value is a safe scale. The boundary moment carries the
/// outward-normal component n_x, which is exactly 0 by symmetry on a
/// background-box-centred aligned cube (the +x and -x faces integrate the
/// same even-power moment with opposite-signed normals over the same
/// [y,z] domain); scaling by the PURE moment (no normal factor, always
/// positive) rather than by that possibly-vanishing quantity keeps the
/// scaled error well defined on every case, exactly as scaledErr's own
/// refScale argument is designed for.
bool checkTx(const TetClipQuadrature & Q, const TetMesh & M, index_t p)
{
    const int e = 2*p;

    const real_t volOracle = tetVolumeOracleMoment(M, e,e,e, p);
    const real_t volErr = scaledErr(volMoment(Q,e,e,e), volOracle, math::abs(volOracle));
    const bool passVol = volErr <= 1e-12;
    gsInfo << "Tx volume moment x^" << e << " y^" << e << " z^" << e << ": " << fmtSci(volErr)
          << "  " << (passVol ? "PASS" : "FAIL") << "\n";

    const real_t bdrOracleN    = triBoundaryOracleMoment(M, e,e,e, 0, p);
    const real_t bdrOraclePure = triBoundaryOracleMoment(M, e,e,e, -1, p);
    const real_t bdrErr = scaledErr(bdrMoment(Q,e,e,e,0), bdrOracleN, bdrOraclePure);
    const bool passBdr = bdrErr <= 1e-12;
    gsInfo << "Tx boundary moment oint x^" << e << " y^" << e << " z^" << e << " nx dS: " << fmtSci(bdrErr)
          << "  " << (passBdr ? "PASS" : "FAIL") << "\n";

    return passVol && passBdr;
}

//----------------------------------------------------------------------------
// Algoim (rule 15) baseline: an independent adaptive-octree quadrature on
// the same background grid, compared cell-by-cell against the exact clip
// rule. INFO only -- see the file header ("Algoim baseline") for scope.
//----------------------------------------------------------------------------

/// False iff some boundary triangle of \a M has all three vertices lying
/// exactly (to 1e-12) on the same background grid plane of \a grid, i.e.
/// phi is identically zero on a whole cell face there. Algoim's
/// height-function construction (called from gsAlgoimAdaptiveRule::mapTo())
/// can fail to terminate on that degeneracy and exhaust
/// memory -- see the split-plane-on-the-interface guard,
/// optional/gsAlgoim/gsAlgoimAdaptiveRule.h:452-457, which shifts a
/// subdivision plane off an interface that passes exactly through a box
/// centre for the same reason. See the 2D analogue algoimBaselineUsable in
/// optional/gsOpenCascade/examples/immersed_gauss_green_quadrature_example.cpp.
bool algoimBaselineUsable(const TetMesh & M, const Grid3 & grid)
{
    std::vector<real_t> X, Y, Z;
    gridLines(grid, X, Y, Z);
    const std::vector<real_t> * L[3] = {&X, &Y, &Z};

    for (const std::array<index_t,4> & tri : M.bdrTri)
        for (int d = 0; d != 3; ++d)
        {
            const std::vector<real_t> & Ld = *L[d];
            for (index_t i = 0; i <= grid.n; ++i)
            {
                bool onPlane = true;
                for (int v = 0; v != 3 && onPlane; ++v)
                    onPlane = (math::abs(M.P[tri[v]][d] - Ld[i]) <= 1e-12);
                if (onPlane) return false;
            }
        }
    return true;
}

/// Builds the gsAlgoimAdaptiveRule (rule 15) baseline on every cell of
/// \a grid and prints the A0-A5 INFO comparison against the exact clip
/// quadrature \a Q. On return, \a algoimVolNodes/Weights and
/// \a algoimSrfNodes/Weights hold the concatenated Algoim volume and
/// surface-cut nodes over all cells, for --plot.
void runAlgoimBaseline(const TetMesh & M, const Grid3 & grid, index_t p, index_t maxDepth,
                       const TetClipQuadrature & Q, real_t V_mesh, real_t P_mesh, real_t clipElapsed,
                       const std::vector<std::array<int,3> > & moments,
                       const std::vector<real_t> & oracleMoments,
                       gsMatrix<real_t> & algoimVolNodes, gsVector<real_t> & algoimVolWeights,
                       gsMatrix<real_t> & algoimSrfNodes, gsVector<real_t> & algoimSrfWeights)
{
    // The option set quadrature rule 15 forwards to gsAlgoimAdaptiveRule:
    // quA=1.0, quB=1 (the gsExprAssembler option defaults, gsExprAssembler.h:1049-1050);
    // maxDepth stays at its gsAlgoimAdaptiveRule default (0) unless overridden by --maxDepth.
    const real_t quA = 1.0;
    const index_t quB = 1;

    // --- A0: signed-distance level set from the driver's own outward-oriented
    // boundary triangles, in the mesh's own coordinates (no unit-box
    // normalization: the per-cell comparison below needs both rules to see
    // the same geometry in the same coordinates). Only nodes referenced by
    // M.bdrTri are added -- M.P also holds interior tet nodes, which would
    // otherwise become isolated vertices.
    gsStopwatch lsClk;
    gsSurfMesh<real_t> surf;
    std::vector<gsSurfMesh<real_t>::Vertex> vmap(M.P.size());
    for (const std::array<index_t,4> & tri : M.bdrTri)
    {
        for (int v = 0; v != 3; ++v)
        {
            const index_t idx = tri[v];
            if (!vmap[idx].is_valid())
            {
                gsSurfMesh<real_t>::Point pt;
                pt << M.P[idx][0], M.P[idx][1], M.P[idx][2];
                vmap[idx] = surf.add_vertex(pt);
            }
        }
        const gsSurfMesh<real_t>::Face f = surf.add_triangle(vmap[tri[0]], vmap[tri[1]], vmap[tri[2]]);
        GISMO_ENSURE(f.is_valid(), "runAlgoimBaseline: add_triangle rejected a boundary face "
                    "(would create a complex edge) -- the derived boundary mesh is not manifold.");
    }

    gsMatrix<real_t> bbox(3,2);
    bbox(0,0) = grid.x0;                 bbox(0,1) = grid.x0 + grid.n*grid.h;
    bbox(1,0) = grid.y0;                 bbox(1,1) = grid.y0 + grid.n*grid.h;
    bbox(2,0) = grid.z0;                 bbox(2,1) = grid.z0 + grid.n*grid.h;
    gsMeshSignedDist<real_t> phi(surf, bbox);   // BVH built here
    const real_t lsTime = lsClk.stop();

    // Sign check: the mesh bbox is strictly inside the background box
    // (GISMO_ENSURE'd in main), so both facts below are guaranteed; a
    // failure here is a level-set construction bug, not a data issue.
    const index_t T = (index_t)M.tet.size();
    for (index_t ti : {index_t(0), T/2, T-1})
    {
        Vec3 c = {0,0,0};
        for (int v = 0; v != 4; ++v) c = add3(c, scale3(0.25, M.P[M.tet[ti][v]]));
        gsMatrix<real_t> pt(3,1); pt << c[0], c[1], c[2];
        gsMatrix<real_t> val; phi.eval_into(pt, val);
        GISMO_ENSURE(val(0,0) < 0, "runAlgoimBaseline: sign check failed at tet " << ti
                    << " centroid (phi=" << val(0,0) << " >= 0).");
    }
    for (int cx = 0; cx != 2; ++cx)
    for (int cy = 0; cy != 2; ++cy)
    for (int cz = 0; cz != 2; ++cz)
    {
        gsMatrix<real_t> pt(3,1);
        pt << (cx ? bbox(0,1) : bbox(0,0)), (cy ? bbox(1,1) : bbox(1,0)), (cz ? bbox(2,1) : bbox(2,0));
        gsMatrix<real_t> val; phi.eval_into(pt, val);
        GISMO_ENSURE(val(0,0) > 0, "runAlgoimBaseline: sign check failed at a background-box "
                    "corner (phi=" << val(0,0) << " <= 0).");
    }

    gsInfo << "A0 level set: bdr vertices=" << surf.n_vertices() << " faces=" << surf.n_faces()
          << " sign check PASS build time=" << fmtSci(lsTime) << "s  INFO\n";

    // --- Algoim rule set (quA/quB declared at the top of this function).
    // At maxDepth=0 the adaptive rule never classifies a box, going straight
    // to the leaf branch (gsAlgoimAdaptiveRule.h:433-441); everything but
    // maxDepth stays at gsAlgoimAdaptiveRule::defaultOptions()
    // (indicator="integralChange", indicatorTol=1e-2, nFallback=0 -> p+1,
    // LipschitzConstant=1.0, splitShift=0.05). No custom box classifier is
    // installed: deeper than depth 0 the built-in midpoint+Lipschitz test is
    // rigorous here since phi is a true signed distance, hence 1-Lipschitz.
    gsOptionList o = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
    o.setInt ("dim", -1);              // volume {phi<0}
    o.setReal("quA", quA);
    o.setInt ("quB", quB);
    o.setInt ("maxDepth", maxDepth);
    gsAlgoimAdaptiveRule<real_t> volRule(phi, (short_t)p, o);
    gsOptionList osrf = o; osrf.setInt("dim", 3);   // surface {phi==0}
    gsAlgoimAdaptiveRule<real_t> srfRule(phi, (short_t)p, osrf);

    std::vector<real_t> X, Y, Z;
    gridLines(grid, X, Y, Z);
    const index_t n = grid.n;
    const real_t h3 = grid.h*grid.h*grid.h;
    const index_t nCells = n*n*n;

    std::vector<gsMatrix<real_t> > cellVolPts((size_t)nCells);
    std::vector<gsVector<real_t> > cellVolWts((size_t)nCells);
    std::vector<gsMatrix<real_t> > cellSrfPts((size_t)nCells);
    std::vector<gsVector<real_t> > cellSrfWts((size_t)nCells);

    gsStopwatch volClk;
    for (index_t k = 0; k != n; ++k)
    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const index_t id = i + n*(j + n*k);
        gsVector<real_t> lower(3), upper(3);
        lower << X[i], Y[j], Z[k];
        upper << X[i+1], Y[j+1], Z[k+1];
        volRule.mapTo(lower, upper, cellVolPts[id], cellVolWts[id]);
    }
    const real_t volTime = volClk.stop();

    gsStopwatch srfClk;
    for (index_t k = 0; k != n; ++k)
    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const index_t id = i + n*(j + n*k);
        gsVector<real_t> lower(3), upper(3);
        lower << X[i], Y[j], Z[k];
        upper << X[i+1], Y[j+1], Z[k+1];
        gsMatrix<real_t> interior; gsVector<real_t> interiorW;
        srfRule.mapToSeparated(lower, upper, interior, interiorW, cellSrfPts[id], cellSrfWts[id]);
    }
    const real_t srfTime = srfClk.stop();

    // --- A1: total volume ---
    KahanSum clipVolTotal, algoimVolTotal;
    index_t totalVolNodes = 0;
    for (const CellRule3 & r : Q.vol)
        for (index_t k = 0; k != r.weights.size(); ++k) clipVolTotal.add(r.weights[k]);
    for (index_t id = 0; id != nCells; ++id)
    {
        for (index_t k = 0; k != cellVolWts[id].size(); ++k) algoimVolTotal.add(cellVolWts[id][k]);
        totalVolNodes += cellVolPts[id].cols();
    }
    const real_t V_clip = clipVolTotal.value(), V_algoim = algoimVolTotal.value();
    gsInfo << "A1 volume: V_mesh=" << fmtSci(V_mesh) << " clip=" << fmtSci(V_clip)
          << " (scaled err " << fmtSci(scaledErr(V_clip, V_mesh, V_mesh)) << ") Algoim="
          << fmtSci(V_algoim) << " (scaled err " << fmtSci(scaledErr(V_algoim, V_mesh, V_mesh))
          << ")  INFO\n";

    // --- A2: volume moments a+b+c<=2p vs the unclipped-mesh oracle (T2) ---
    real_t clipMomErr = 0, algoimMomErr = 0;
    for (size_t k = 0; k != moments.size(); ++k)
    {
        const std::array<int,3> & m = moments[k];
        clipMomErr = math::max(clipMomErr,
            scaledErr(volMoment(Q, m[0], m[1], m[2]), oracleMoments[k], V_mesh));

        KahanSum algMom;
        for (index_t id = 0; id != nCells; ++id)
        {
            const gsMatrix<real_t> & nd = cellVolPts[id];
            const gsVector<real_t> & wt = cellVolWts[id];
            for (index_t kk = 0; kk != wt.size(); ++kk)
                algMom.add(wt[kk]*std::pow(nd(0,kk),m[0])*std::pow(nd(1,kk),m[1])*std::pow(nd(2,kk),m[2]));
        }
        algoimMomErr = math::max(algoimMomErr, scaledErr(algMom.value(), oracleMoments[k], V_mesh));
    }
    gsInfo << "A2 moments a+b+c<=2p vs oracle: clip max scaled err=" << fmtSci(clipMomErr)
          << " Algoim max scaled err=" << fmtSci(algoimMomErr) << "  INFO\n";

    // --- A3: per-cell |V_clip-V_algoim|/h^3, worst 5 ---
    struct WorstCell { index_t i,j,k; real_t vClip, vAlgoim, diff; };
    std::vector<WorstCell> worst; worst.reserve(nCells);
    real_t maxDiff = 0;
    for (index_t k = 0; k != n; ++k)
    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const index_t id = i + n*(j + n*k);
        const real_t vClip = cellVolume(Q.vol[id]);
        const real_t vAlgoim = cellVolWts[id].size() ? cellVolWts[id].sum() : 0.0;
        const real_t diff = math::abs(vClip - vAlgoim);
        maxDiff = math::max(maxDiff, diff/h3);
        worst.push_back(WorstCell{i,j,k, vClip, vAlgoim, diff});
    }
    std::sort(worst.begin(), worst.end(),
             [](const WorstCell & a, const WorstCell & b) { return a.diff > b.diff; });

    gsInfo << "A3 per-cell max |V_clip-V_algoim|/h^3=" << fmtSci(maxDiff) << "  INFO\n";
    gsInfo << "A3 worst cells:\n";
    static const char * statusName[] = {"Cut", "Full"};   // Empty=-1 handled below
    for (index_t r = 0; r != 5 && r != (index_t)worst.size(); ++r)
    {
        const WorstCell & w = worst[r];
        const index_t id = w.i + n*(w.j + n*w.k);
        const real_t cx = 0.5*(X[w.i]+X[w.i+1]), cy = 0.5*(Y[w.j]+Y[w.j+1]), cz = 0.5*(Z[w.k]+Z[w.k+1]);
        const int st = Q.status[id];
        const std::string stName = (Empty == st) ? "Empty" : statusName[st];
        gsInfo << "  (" << w.i << "," << w.j << "," << w.k << ") centre=(" << cx << "," << cy << "," << cz
              << ") status=" << stName << " clip/h^3=" << fmtSci(w.vClip/h3)
              << " Algoim/h^3=" << fmtSci(w.vAlgoim/h3) << " diff/h^3=" << fmtSci(w.diff/h3) << "\n";
    }

    // --- A4: surface area ---
    KahanSum srfTotal;
    index_t totalSrfNodes = 0;
    for (index_t id = 0; id != nCells; ++id)
    {
        for (index_t k = 0; k != cellSrfWts[id].size(); ++k) srfTotal.add(cellSrfWts[id][k]);
        totalSrfNodes += cellSrfPts[id].cols();
    }
    const real_t areaAlgoim = srfTotal.value();
    gsInfo << "A4 surface: P_mesh=" << fmtSci(P_mesh) << " Algoim area=" << fmtSci(areaAlgoim)
          << " (scaled err " << fmtSci(scaledErr(areaAlgoim, P_mesh, P_mesh))
          << ") nZeroSurfaceLeaves=" << srfRule.stats().nZeroSurfaceLeaves
          << " unresolvedSurfaceMeasure=" << fmtSci(srfRule.stats().unresolvedSurfaceMeasure) << "  INFO\n";

    // --- A5: timings ---
    gsInfo << "A5 timing: clip=" << fmtSci(clipElapsed) << "s levelset build=" << fmtSci(lsTime)
          << "s Algoim volume=" << fmtSci(volTime) << "s (nodes=" << totalVolNodes
          << ", nFallbackLeaves=" << volRule.stats().nFallbackLeaves
          << ", nSubBoxes=" << volRule.stats().nSubBoxes << ") Algoim surface=" << fmtSci(srfTime)
          << "s (nodes=" << totalSrfNodes << ")  INFO\n";

    // Concatenate for --plot.
    index_t volTotal = 0, srfTotalN = 0;
    for (const gsVector<real_t> & w : cellVolWts) volTotal += w.size();
    for (const gsVector<real_t> & w : cellSrfWts) srfTotalN += w.size();
    algoimVolNodes.resize(3, volTotal); algoimVolWeights.resize(volTotal);
    algoimSrfNodes.resize(3, srfTotalN); algoimSrfWeights.resize(srfTotalN);
    index_t colV = 0, colS = 0;
    for (index_t id = 0; id != nCells; ++id)
    {
        const index_t mV = cellVolWts[id].size();
        if (mV) { algoimVolNodes.block(0,colV,3,mV) = cellVolPts[id]; algoimVolWeights.segment(colV,mV) = cellVolWts[id]; colV += mV; }
        const index_t mS = cellSrfWts[id].size();
        if (mS) { algoimSrfNodes.block(0,colS,3,mS) = cellSrfPts[id]; algoimSrfWeights.segment(colS,mS) = cellSrfWts[id]; colS += mS; }
    }
}

/// Writes the clip volume/boundary quadrature nodes, and (if the Algoim
/// pointers are non-null) the Algoim volume/surface nodes, as ParaView point
/// sets in the current working directory, weight carried as the point value.
/// See the file header's Size trap: at the default p=2 every clip sub-tet
/// carries 392 nodes.
void writePlot(const TetClipQuadrature & Q,
              const gsMatrix<real_t> * algoimVolNodes, const gsVector<real_t> * algoimVolWeights,
              const gsMatrix<real_t> * algoimSrfNodes, const gsVector<real_t> * algoimSrfWeights)
{
    index_t totalVol = 0, totalBdr = 0;
    for (const CellRule3 & r : Q.vol) totalVol += r.weights.size();
    for (const CellBdrRule3 & r : Q.bdr) totalBdr += r.weights.size();

    gsMatrix<real_t> Xv(1,totalVol), Yv(1,totalVol), Zv(1,totalVol), Vv(1,totalVol);
    index_t col = 0;
    for (const CellRule3 & r : Q.vol)
        for (index_t k = 0; k != r.weights.size(); ++k)
        { Xv(0,col)=r.nodes(0,k); Yv(0,col)=r.nodes(1,k); Zv(0,col)=r.nodes(2,k); Vv(0,col)=r.weights[k]; ++col; }
    gsWriteParaviewPoints(Xv, Yv, Zv, Vv, "tetclip_volume_nodes", 12);

    gsMatrix<real_t> Xb(1,totalBdr), Yb(1,totalBdr), Zb(1,totalBdr), Vb(1,totalBdr);
    col = 0;
    for (const CellBdrRule3 & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
        { Xb(0,col)=r.nodes(0,k); Yb(0,col)=r.nodes(1,k); Zb(0,col)=r.nodes(2,k); Vb(0,col)=r.weights[k]; ++col; }
    gsWriteParaviewPoints(Xb, Yb, Zb, Vb, "tetclip_boundary_nodes", 12);

    std::string written = "tetclip_volume_nodes.vtp tetclip_boundary_nodes.vtp";

    if (algoimVolNodes)
    {
        const index_t m = algoimVolWeights->size();
        gsMatrix<real_t> Xa(1,m), Ya(1,m), Za(1,m), Va(1,m);
        for (index_t k = 0; k != m; ++k)
        { Xa(0,k)=(*algoimVolNodes)(0,k); Ya(0,k)=(*algoimVolNodes)(1,k); Za(0,k)=(*algoimVolNodes)(2,k); Va(0,k)=(*algoimVolWeights)[k]; }
        gsWriteParaviewPoints(Xa, Ya, Za, Va, "algoim_volume_nodes", 12);
        written += " algoim_volume_nodes.vtp";
    }
    if (algoimSrfNodes)
    {
        const index_t m = algoimSrfWeights->size();
        gsMatrix<real_t> Xa(1,m), Ya(1,m), Za(1,m), Va(1,m);
        for (index_t k = 0; k != m; ++k)
        { Xa(0,k)=(*algoimSrfNodes)(0,k); Ya(0,k)=(*algoimSrfNodes)(1,k); Za(0,k)=(*algoimSrfNodes)(2,k); Va(0,k)=(*algoimSrfWeights)[k]; }
        gsWriteParaviewPoints(Xa, Ya, Za, Va, "algoim_surface_nodes", 12);
        written += " algoim_surface_nodes.vtp";
    }
    gsInfo << "ParaView output written: " << written << "\n";
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    index_t p = 2;
    index_t r = 1;
    index_t n0 = 4;
    std::string caseName = "rotcube";
    std::string fileName;
    bool fit = false;
    bool algoim = false;
    bool plot = false;
    index_t maxDepth = 0;

    gsCmdLine cmd("3D cut-cell quadrature from a tetrahedral mesh by exact tet-cell clipping "
                 "against a uniform background grid, verified against unclipped-mesh oracles "
                 "and (cube cases) analytic moments.");
    cmd.addInt   ("k", "degree", "Collapsed-Gauss rule degree parameter p", p);
    cmd.addInt   ("r", "refine", "Refinement level: n0*2^r cells per direction", r);
    cmd.addInt   ("",  "n0",     "Cells per direction at r = 0", n0);
    cmd.addString("",  "case",   "Test case: rotcube | cube | sphere | none", caseName);
    cmd.addString("f", "file",   "ASCII gmsh MSH 4.1 tet-mesh file (overrides --case's default)", fileName);
    cmd.addSwitch("fit", "Background box = cube of side 1.1*Lmax centred at the mesh bbox centre", fit);
    cmd.addSwitch("algoim", "Run the gsAlgoimAdaptiveRule (rule 15) baseline on the same grid and "
                  "print an A0-A5 INFO comparison against the exact clip rule; never affects "
                  "ALL CHECKS PASS or the exit code. Skipped with a message if a boundary "
                  "triangle lies exactly on a background grid plane (--case cube).", algoim);
    cmd.addInt   ("",  "maxDepth", "Algoim adaptive subdivision depth (--algoim only); "
                  "0 = plain rule-15 defaults", maxDepth);
    cmd.addSwitch("plot", "Write the clip (and, with --algoim, Algoim) quadrature nodes as "
                  "ParaView point sets in the current directory. At the default p=2 every clip "
                  "sub-tet carries 392 nodes; use a small -k/-r (e.g. -k 1 -r 0) here.", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (p < 1)  { gsWarn << "-k/--degree must be >= 1\n"; return EXIT_FAILURE; }
    if (r < 0)  { gsWarn << "-r/--refine must be >= 0\n"; return EXIT_FAILURE; }
    if (n0 < 1) { gsWarn << "--n0 must be >= 1\n"; return EXIT_FAILURE; }
    if ("rotcube" != caseName && "cube" != caseName && "sphere" != caseName && "none" != caseName)
    { gsWarn << "--case must be one of rotcube|cube|sphere|none\n"; return EXIT_FAILURE; }
    if ("none" == caseName && fileName.empty())
    { gsWarn << "--case none requires -f <file>\n"; return EXIT_FAILURE; }
    if (fit && "cube" == caseName)
    { gsWarn << "--fit is rejected together with --case cube\n"; return EXIT_FAILURE; }
    if (maxDepth < 0) { gsWarn << "--maxDepth must be >= 0\n"; return EXIT_FAILURE; }

    std::string defaultFile;
    if ("rotcube" == caseName) defaultFile = "volumes/tetmesh_cube_rotated.msh";
    else if ("cube" == caseName) defaultFile = "volumes/tetmesh_cube_aligned.msh";
    else if ("sphere" == caseName) defaultFile = "volumes/tetmesh_sphere.msh";
    const std::string requestedFile = fileName.empty() ? defaultFile : fileName;
    const std::string resolvedFile = gsFileManager::find(requestedFile);
    if (resolvedFile.empty())
    { gsWarn << "-f " << requestedFile << " not found\n"; return EXIT_FAILURE; }

    const TetMesh M = readMsh41(resolvedFile);

    const index_t n = n0 * (index_t(1) << r);
    Grid3 grid;
    if (fit)
    {
        real_t Lmax = 0;
        for (int d = 0; d != 3; ++d) Lmax = math::max(Lmax, M.hi[d]-M.lo[d]);
        const real_t side = 1.1*Lmax;
        grid.x0 = 0.5*(M.lo[0]+M.hi[0]) - 0.5*side;
        grid.y0 = 0.5*(M.lo[1]+M.hi[1]) - 0.5*side;
        grid.z0 = 0.5*(M.lo[2]+M.hi[2]) - 0.5*side;
        grid.h = side/(real_t)n;
        grid.n = n;
    }
    else
    {
        grid.x0 = grid.y0 = grid.z0 = -1.0;
        grid.h = 2.0/(real_t)n;
        grid.n = n;
    }

    GISMO_ENSURE(M.lo[0] > grid.x0 && M.hi[0] < grid.x0 + grid.n*grid.h &&
                M.lo[1] > grid.y0 && M.hi[1] < grid.y0 + grid.n*grid.h &&
                M.lo[2] > grid.z0 && M.hi[2] < grid.z0 + grid.n*grid.h,
                "main: mesh bbox [" << M.lo[0] << "," << M.hi[0] << "] x ["
                  << M.lo[1] << "," << M.hi[1] << "] x [" << M.lo[2] << "," << M.hi[2]
                  << "] does not lie strictly inside the background box ["
                  << grid.x0 << "," << grid.x0+grid.n*grid.h << "] x ["
                  << grid.y0 << "," << grid.y0+grid.n*grid.h << "] x ["
                  << grid.z0 << "," << grid.z0+grid.n*grid.h << "]");

    gsStopwatch clk;
    TetClipStats stats;
    const TetClipQuadrature Q = tetClipQuadrature(M, grid, p, stats);
    const real_t elapsed = clk.stop();

    index_t Vb = 0, Fb = (index_t)M.bdrTri.size();
    std::set<index_t> bNodes;
    std::set<std::pair<index_t,index_t> > bEdges;
    for (const std::array<index_t,4> & tri : M.bdrTri)
    {
        for (int v = 0; v != 3; ++v) bNodes.insert(tri[v]);
        for (int e = 0; e != 3; ++e)
        {
            index_t u = tri[e], w = tri[(e+1)%3];
            if (u > w) std::swap(u,w);
            bEdges.insert(std::make_pair(u,w));
        }
    }
    Vb = (index_t)bNodes.size();
    const index_t Eb = (index_t)bEdges.size();
    const index_t chi = Vb - Eb + Fb;

    gsInfo << "mesh: file=" << resolvedFile << " nodes=" << M.nNodesFile << " tets=" << M.tet.size()
          << " boundary triangles=" << M.bdrTri.size() << " chi=Vb-Eb+Fb=" << chi << "  INFO\n";

    const index_t volDeg = volDegree(p), bdrDeg = bdrDegree(p);
    const index_t mu = ceilHalf(volDeg+3), mv = ceilHalf(volDeg+2), mw = ceilHalf(volDeg+1);
    const index_t mub = ceilHalf(bdrDeg+2), mvb = ceilHalf(bdrDeg+1);
    gsInfo << "stats: pieces=" << stats.nPieces << " bdr pieces=" << stats.nBdrPieces
          << " cut=" << stats.nCut << " full=" << stats.nFull << " empty=" << stats.nEmpty
          << " volDeg=" << volDeg << " bdrDeg=" << bdrDeg
          << " vol nodes/piece=" << mu*mv*mw << " bdr nodes/piece=" << mub*mvb
          << " time=" << fmtSci(elapsed) << "s  INFO\n";

    const real_t L = boxCornerAbsMax(grid);
    const real_t V_mesh = meshVolumeExact(M);
    const real_t P_mesh = unclippedBoundaryArea(M);
    const std::vector<std::array<int,3> > moments = momentSet3(p);

    real_t s = 0; Vec3 cvec = {0,0,0}; real_t R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    if ("rotcube" == caseName) { s = ROTCUBE_SIDE; cvec = ROTCUBE_TRANSLATE; rotcubeMatrix(R); }

    bool allPass = true;
    allPass = checkT1(stats, L) && allPass;
    std::vector<real_t> oracleMoments;
    allPass = checkT2(Q, M, p, V_mesh, oracleMoments, moments) && allPass;
    allPass = checkT3(Q, caseName, p, moments, oracleMoments, s, cvec, R) && allPass;
    allPass = checkT4(Q, L) && allPass;
    allPass = checkT5(Q, V_mesh) && allPass;
    allPass = checkT6(Q, P_mesh) && allPass;
    allPass = checkT7(Q) && allPass;
    allPass = checkT8(Q) && allPass;
    allPass = checkT9(M, grid, caseName, s, cvec, R) && allPass;
    allPass = checkTx(Q, M, p) && allPass;

    if ("sphere" == caseName)
    {
        const real_t exactVol = (4.0/3.0)*(real_t)EIGEN_PI*SPHERE_RADIUS*SPHERE_RADIUS*SPHERE_RADIUS;
        const real_t diff = scaledErr(V_mesh, exactVol, exactVol);
        gsInfo << "sphere info: mesh volume=" << fmtSci(V_mesh) << " 4/3 pi r^3=" << fmtSci(exactVol)
              << " scaled diff=" << fmtSci(diff) << "  INFO\n";
    }

    gsMatrix<real_t> algoimVolNodes, algoimSrfNodes;
    gsVector<real_t> algoimVolWeights, algoimSrfWeights;
    bool haveAlgoim = false;

    if (algoim)
    {
        if (!algoimBaselineUsable(M, grid))
        {
            gsWarn << "Algoim baseline skipped: a boundary triangle lies exactly on a background "
                      "grid plane (phi is identically zero on a whole cell face there), which "
                      "makes gsAlgoimAdaptiveRule::mapTo exhaust memory. See the file header.\n";
            gsInfo << "Algoim baseline skipped\n";
        }
        else
        {
            runAlgoimBaseline(M, grid, p, maxDepth, Q, V_mesh, P_mesh, elapsed, moments, oracleMoments,
                              algoimVolNodes, algoimVolWeights, algoimSrfNodes, algoimSrfWeights);
            haveAlgoim = true;
        }
    }

    if (plot)
        writePlot(Q, haveAlgoim ? &algoimVolNodes : nullptr, haveAlgoim ? &algoimVolWeights : nullptr,
                  haveAlgoim ? &algoimSrfNodes : nullptr, haveAlgoim ? &algoimSrfWeights : nullptr);

    if (allPass) gsInfo << "ALL CHECKS PASS\n";
    else         gsInfo << "SOME CHECKS FAILED\n";

    return allPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
