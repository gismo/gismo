/** @file immersed_tetmesh_poisson_example.cpp

    @brief Driver for the tet-mesh immersed Poisson solver built on the
    streamed lookup-rule quadrature machinery of gsImmersedLookupRule.h
    (namespace gsTetClip). Dispatches on `--study`:

      - `check`: verifies (a) that ClipStreamer's per-cell
        volume/boundary rules, wrapped by VolLookupRule/BdrLookupRule, are
        bitwise identical to gsTetMeshClip.h's reference
        tetClipQuadrature(), for both the on-the-fly streamer and a table
        built from it; (b) that QuadratureScope/VolumeQuadratureScope/
        BoundaryQuadratureScope save and restore assembler/evaluator state
        correctly, that the partition-of-unity rhs sum and the evaluator's
        volume integral under the scopes both match the reference clipped
        volume to 1e-12 relative, that a boundary assembly with
        BdrNormalField registered reproduces Int u |n| dS = Int u dS to
        1e-12 relative, and that a ghost-penalty matrix assembled
        through the scopes is bitwise identical to one assembled without
        any custom quadrature installed; (c) that BdrNormalField throws on
        any point set other than the exact nodes a BdrLookupRule produced;
        (d) that gsTetClipSignDomain.h's TetClipSignDomain reproduces
        ClipStreamer's per-cell status as single-cell kd-tree leaf signs
        with zero mismatches, correct element/ghost/skeleton-face counts,
        the boundary-piece invariant, a firing multi-cell guard and an
        exact identity geometry, for sphere/rotcube at r=0..3.
      - `volume`: the volume/area/flux/moment GATE on the streamed clip
        (`--mode clip`), Tchakaloff-compressed (`--mode tchakaloff`),
        moment-fitted (`--mode momrule`) and Algoim (`--mode algoim`,
        sphere-only, also on `--geo affine|bubble` via the pulled-back level
        set phi o G -- see the Curved backgrounds paragraph below) cell
        rules over the same uniform background grid.
        Prints one `GATE` row per check (`PASS`/`FAIL`/`REPORT`, or
        `DEFERRED` for `--mode tchakaloff`, see below), one `GATEDIAG`
        diagnostics row per (case, mode, r), and a final
        `GATE SUMMARY required=<N> pass=<N> fail=<N>` line; the exit
        code is nonzero iff any required row FAILs. `fluxmom_global`/
        `fluxmom_cell` (the vector-moment, divergence-theorem check the
        Nitsche flux term needs -- scalar-only `moments_vol`/`moments_bdr`
        cannot see a wrong served normal) are Required for EVERY mode
        including momrule, so a `--mode all`/`--mode momrule` run's exit
        code depends on them too: `--mode momrule` FAILs both rows on both
        `--case sphere` and `--case rotcube` (see `poisson` below for why
        that does not also block the solve there).
        Tchakaloff's rows are computed and printed exactly like the other
        modes' (value/ref/relerr/tol, failCells, per-cell INFO lines), but
        every otherwise-Required row prints `wouldPass=0|1 DEFERRED`
        instead of `PASS`/`FAIL` and is excluded from `GATE SUMMARY`
        (`required=0` for a tchakaloff-only run): its NNLS compressor's KKT
        stopping test bounds only the dual gradient, not the primal
        residual (see the TODO on `nnlsLawsonHansonImpl` in
        `gsTchakaloffRule.h`), which can leave a cut cell's residual well
        above `resTol` even though the test itself is satisfied. One
        `GATE-DEFERRED mode=tchakaloff reason=nnls-kkt-stall` line prints
        per (case, r).
        For `--geo affine|bubble` and `--mode clip`, the rows come from
        `gatePhysical` instead of `gateTetMode`/`tetModeTotals` (see the
        Curved backgrounds paragraph below): `pb_finite`/`pb_residual`/
        `pb_inside`/`pb_orientation`/`detJ_min` are always Required;
        `volume`/`area`/`flux` are Required for `--geo affine` (exact up
        to rounding there) and Report for `--geo bubble`, where they
        differ from V_mesh/A_mesh/3V_mesh by a genuine GEOMETRIC quantity
        (Omega_s = G(the mesh polyhedron) is bounded by curved images of
        the mesh's straight edges, not by the straight mesh itself)
        rather than by quadrature error. `fluxmom_global` is Required for
        `--geo affine` and, for `--geo bubble`, Required at p >=
        `BUBBLE_FLUXMOM_REQUIRED_MIN_P` (the measured value sits at the
        rounding floor there, per the measured sweep documented at
        `BUBBLE_FLUXMOM_TOL_REQUIRED`) and Report below that p (its own
        reference rule is not exact under a curved G, see `gatePhysical`
        and `BUBBLE_FLUXMOM_TOL`). No `moments_*`/
        `fluxmom_cell`/`residual_*`/`minweight_*`/`run` rows exist on this
        path (clip only, no compressor).
      - `poisson`: for each (case, mode, r) with r = 0..rMax, builds the
        mode's quadrature tables (R1), runs the `--study volume` gate (R2)
        on those SAME objects and prints its row (R3); only if the gate
        passed, solves the 3D immersed Poisson problem. Exception: for
        `--mode momrule`, a FAIL confined to `fluxmom_global`/`fluxmom_cell`
        does not withhold the solve (`GATE-NOTE mode=momrule fluxmom=FAIL
        solve=report` prints instead of `POISSON-SKIP ... reason=gate-FAIL`)
        -- momrule's served normal is only a pseudonormal at polyhedral
        edges/vertices, a known, already-documented gap (see the Per-mode
        list below), not a reason to skip the solve on every case. Any
        OTHER Required row failing (for momrule or any other mode) still
        skips it. For `--geo affine|bubble` the gate is `gatePhysical`
        instead (see `volume` above and the Curved backgrounds paragraph
        below): the momrule exception does not apply there, since momrule
        itself is rejected outright for a curved G.

          -Delta u = f in Omega,  u = g on d(Omega),

        with symmetric Nitsche on the immersed boundary and a ghost penalty
        on the faces of cut cells, on a degree-p Cartesian B-spline
        background over [-1,1]^3 (dofs outside Omega eliminated). Omega is
        the tet-mesh polyhedron for the tet modes (clip/tchakaloff/momrule)
        and the analytic ball `(x-0.03)^2+(y+0.02)^2+(z-0.01)^2=0.3025` for
        `algoim`. Manufactured solution (smooth on R^3, non-polynomial,
        non-symmetric):

          u(x,y,z) = sin(pi*x/2) * cos(pi*y/3) * exp(z/2),
          f(x,y,z) = -Delta u = (13*pi^2/36 - 1/4) * u  (~ 3.3141 * u),
          grad u   = exp(z/2) * ( (pi/2) cos(pi*x/2)cos(pi*y/3),
                                  -(pi/3) sin(pi*x/2)sin(pi*y/3),
                                   (1/2) sin(pi*x/2)cos(pi*y/3) ),
          g = u on d(Omega).

        Weak form (find u_h with, for every test v_h):

          (grad u_h, grad v_h)_Omega - <grad u_h . n, v_h>_Gamma
            - <u_h, grad v_h . n>_Gamma + (gamma/h) <u_h, v_h>_Gamma
            + g_h(u_h, v_h)
          = (f, v_h)_Omega - <g, grad v_h . n>_Gamma + (gamma/h) <g, v_h>_Gamma,

        gamma = 6(p+1)^2 by default (`--gamma`), Gamma the immersed
        boundary, n its stored outward unit normal (gsTetClip::BdrNormalField).
        On the default background (`--geo identity`, G the identity, J = I)
        n is already physical and no Nanson factor is needed. On a curved
        background (`--geo affine|bubble`, see below) n and the boundary
        weight are made physical by gsTetClip::PullbackBdrSource, via
        Nanson's formula n dA = det(J) J^{-T} N dA0, before assembly; every
        volume/boundary/error expression in the solve uses G (the boundary
        normal is served already physical, uncomposed, at the rule's own
        parametric nodes). The ghost penalty
        g_h(u,v) = gamma_g h^(2p-1) sum_F Int_F [d_n^p u][d_n^p v] dS over
        the faces with at least one cut neighbour (`--ghost`, on by
        default; gamma_g = 10^-(p+1) by default, `--ghostCoef`) is, by
        default, assembled through a SECOND, identity map Gid
        (gsTetClip::identityBoxGeometry on the same grid), not through G:
        d_n^p is the parametric p-th
        normal-derivative jump of the background space across a face, and
        assembling it through G instead would only scale it by the smooth
        G-weight the two share across the face (G is C^p there), never
        change which faces see a nonzero jump. This split (Gid vs. G for
        the ghost term) is a driver-local design choice, not a method drawn
        from the literature; `--ghostMap G` (debug) routes it through G
        instead, for comparison. Only order k = p is assembled: the
        background space is C^(p-1), so every lower-order jump is
        identically zero in exact arithmetic, and assembling it would
        inject nothing but the deterministic (faceShift*h)^2
        face-evaluation-offset floor as spurious signal (faceShift = 1e-6,
        `poisson2_ghost_penalty_example.cpp` file header). Both terms match
        `poisson2_nitsche_immersed_example.cpp`/
        `poisson2_ghost_penalty_example.cpp` term by term, but drop their
        `surfMeas` Nanson factor (a known bug with a physical normal on a
        curved geometry map, harmless only because it is identically 1
        there): the boundary weights served by gsTetClip::BdrLookupRule are
        already PHYSICAL surface measure (via PullbackBdrSource when G is
        curved), so no such factor belongs here. The Nitsche penalty's h
        is, for a curved G, h * min_i sigma_min(J(x_i)) sampled on an 11^3
        tensor grid (odd, so it includes the bubble's face centres, where
        |grad b| = 2 is attained); for `identity` it is h itself. Here
        sigma_min(J) <= 1 for the bubble at any eps (J = I + rank one) and 0.912
        for the default affine map, so hNitsche <= h and the penalty is at
        least that of the parametric h. The sampled minimum
        (gsTetClip::minSingularValueJ, a necessary, not sufficient,
        non-degeneracy check) can only over-estimate the true box-wide
        minimum, so hNitsche is never smaller, and the penalty never
        larger, than with the exact sigma_min. This is a sampled heuristic,
        not a coercivity proof.

        Errors (L2 and the H1-SEMINORM, not L2+seminorm) are measured
        against the REFERENCE rule of the mode's own geometry:
        for the tet modes, the streamer's own uncompressed clip rule at a
        finer Full-cell Gauss order (p+3 vs the solve's p+1; on Cut cells
        this is identically the solve rule for `--mode clip`); for
        `algoim`, a serial pass of a SEPARATE `gsAlgoimAdaptiveRule`
        (quA/maxDepth one notch above the solve rule's, LipschitzConstant
        raised to 4.0 -- see the reference-rule code for why the default
        1.0 is not a valid global Lipschitz bound of
        phi = |x-c|^2 - R^2 -- on `--geo identity`; a rigorous per-`--geo`
        bound on |grad(phi o G)| on `--geo affine|bubble`, see
        `pulledBackSphereLipschitz`). `--mode algoim` additionally classifies its
        own background-cell status from the assembled `gsImplicitTrimmedDomain`'s
        OWN Lobatto leaf signs (not from any Algoim-side status), because
        `allElements()`/`beginBdr()` iterate that classification, not
        Algoim's; the resulting per-cell agreement/disagreement counts print
        as one `ALGOIM-CONSISTENCY` line per r (report-only).

        Pre-asymptotic note: at r = 0 (n = 4, h = 0.5) the whole sphere
        (R = 0.55) sits in a 3x3x3-ish block of cut cells, with slivers near
        x ~ -0.52 and x ~ 0.58 -- no cell is Full. That regime is
        pre-asymptotic by construction, which is why the EoC bar below only
        ever looks at the LAST refinement pair.

        Per mode:
          `clip`        solved at every r = 0..rMax (M1's own reference
                        rule); EoC bar (p = 2, last pair): L2 in
                        [2.7,3.3], H1-seminorm in [1.7,2.3].
          `tchakaloff`  DEFERRED (its `--study volume` gate rows are
                        DEFERRED, never PASS): every r prints
                        `POISSON-SKIP ... reason=gate-deferred`, no solve,
                        `EOC ... REPORT`.
          `momrule`, `algoim`  solved and reported, no EoC bar: the served
                        boundary normal is only a pseudonormal at polyhedral
                        edges for momrule, and the Algoim adaptive rule's
                        box-classification bound is not rigorous.
        A `BUDGET` line (clip only, geometric extrapolation of the r=1->2
        growth factor) predicts the r=3 cost right after the r=2 row.
        `--study poisson` exits nonzero iff a tet-mode gate row FAILed, any
        `zeroRows>0`, any non-finite solve, or an `EOC ... FAIL` line was
        printed; REPORT/INFO lines never fail it. Exception: for
        `--mode momrule`, a gate FAIL confined to `fluxmom_global`/
        `fluxmom_cell` alone does NOT fail the run (see the `poisson` bullet
        above) -- the solve still runs and is judged on its own rows.

      Curved backgrounds (`--geo identity|affine|bubble`, default
      `identity`): the background map G: [-1,1]^3 -> R^3 (gsImmersedPullback.h,
      namespace gsTetClip) the whole clip pipeline runs against. `--eps`
      (default 0.3, |eps| < 0.5) is the bubble amplitude, used only for
      `--geo bubble`. `--mesh` overrides the case's own mesh file (e.g. for
      a gmsh-refined mesh); `--case` still selects the case's analytic
      references (V_exact/A_exact, the sphere level set). For a non-identity
      G, the mesh is read in PHYSICAL coordinates and pulled back once per
      (case, r) to a PARAMETRIC one (gsTetClip::pullBack) BEFORE any
      ClipStreamer or assembly runs, so that the clip pipeline still clips
      against the same parametric box [-1,1]^3; a `PULLBACK ... PASS|FAIL`
      row reports the pull-back's own acceptance check (max residual vs. a
      mesh-scaled tolerance, no failed/outside/inverted vertex), and a FAIL
      aborts the study with `PULLBACK-ABORT` before anything is assembled.
      Curved G is supported for `--mode clip` and, `--case sphere` only,
      `--mode algoim` (below) -- momrule/tchakaloff still assume an identity
      background and are rejected up front for `--geo != identity`. Clip's
      own volume/area/flux/moment gate runs in physical space there
      (`gatePhysical`, before any assembly):
      pull-back diagnostics (`pb_finite`/`pb_residual`/`pb_inside`/
      `pb_orientation`) and `detJ_min` (det(J) > 0 at every served volume
      node) are always Required; `volume`/`area`/`flux` are Required for
      `--geo affine` (the map is affine, so every check is exact to
      rounding) and Report for `--geo bubble`: they differ from
      V_mesh/A_mesh/3V_mesh there by a genuine GEOMETRIC quantity
      (Omega_s = G(the mesh polyhedron) is bounded by curved images of
      the mesh's straight edges). `fluxmom_global` is Required for
      `--geo affine` and, for `--geo bubble`, Required at p >=
      `BUBBLE_FLUXMOM_REQUIRED_MIN_P` (per the measured sweep documented
      at `BUBBLE_FLUXMOM_TOL_REQUIRED`, the measured value sits at the
      rounding floor there) and Report below that p, where its own
      reference rule is not exact once G is curved (`BUBBLE_FLUXMOM_TOL`).
      Every `GATE` row printed on this path carries ` geo=<g>`
      right after `mode=clip`; the identity path never does.
      `--mode algoim` on `--geo affine|bubble` integrates the PULLED-BACK
      level set phi-hat(u) = phi(G(u)) directly in PARAMETER space (no tet
      mesh, no pull-back of a mesh vertex): its surface weights and normals
      are made physical by the same `gsTetClip::PullbackBdrSource` Nanson
      map clip uses, and its gate gains two Required rows there, `detJ_min`
      and `normal_dev` (deviation of the served normal from the analytic
      (x-c)/|x-c|), both absent on `--geo identity`.
      `--ghostMap Gid|G` (default `Gid`) and `--dumpSystem <prefix>`
      (default off) are debug-only options: the former swaps which map the
      ghost penalty's `dnk` uses (see the ghost-penalty paragraph above),
      the latter dumps the assembled system (CSC arrays + rhs, raw binary)
      to `<prefix>_<case>_<mode>_r<r>.bin` for `--study poisson`, for
      byte-level comparison across a run. `--ghostMap G` aborts from `dnk`'s
      own ENSURE (`dnk_expr.h`) for any non-axis-aligned G, i.e.
      `--geo affine|bubble`: `dnk` requires a diagonal Jacobian, which, of
      the `--geo` choices, only `identity` provides.

    Example command lines:
      ./immersed_tetmesh_poisson_example
      ./immersed_tetmesh_poisson_example --study check --case sphere -k 2 -r 1
      ./immersed_tetmesh_poisson_example --case rotcube --n0 4 -r 1
      ./immersed_tetmesh_poisson_example --study volume --case all --mode clip -r 3
      ./immersed_tetmesh_poisson_example --study volume --case all --mode tchakaloff -r 3
      ./immersed_tetmesh_poisson_example --study volume --case all --mode momrule -r 3
      ( ulimit -v 8000000; ./immersed_tetmesh_poisson_example --study volume \
        --case sphere --mode algoim -r 3 )
      ./immersed_tetmesh_poisson_example --study poisson --case sphere --mode clip -r 3
      ./immersed_tetmesh_poisson_example --study poisson --case rotcube --mode momrule -r 3
      ( ulimit -v 8000000; ./immersed_tetmesh_poisson_example --study poisson \
        --case sphere --mode algoim -r 3 )
      ./immersed_tetmesh_poisson_example --study poisson --geo affine --mode clip \
        --case sphere -r 2
      ./immersed_tetmesh_poisson_example --study poisson --geo bubble --eps 0.3 \
        --mode clip --case rotcube -r 2
      ./immersed_tetmesh_poisson_example --study poisson --case sphere --mode clip \
        -r 1 --ghostMap G --dumpSystem /tmp/dump

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include "gsImmersedLookupRule.h"
#include "gsTetClipSignDomain.h"
#include "gsImmersedPullback.h"
#include "gsTchakaloffRule.h"
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include <gsDomain/gsMeshLevelSet.h>
#include <gsAssembler/gsMomentRule.h>

#include <algorithm>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace gismo;

namespace {

//----------------------------------------------------------------------------
// CLI configuration.
//----------------------------------------------------------------------------

struct Config
{
    std::string study    = "check";
    std::string caseName = "all";
    std::string mode     = "all";
    index_t p    = 2;
    index_t rMax = 1;
    index_t n0   = 4;
    real_t  gamma      = -1;    // Nitsche penalty; <=0 -> 6*(p+1)^2  (--study poisson)
    index_t ghostOn    = 1;     // 1 = ghost penalty on, 0 = off      (--study poisson)
    real_t  ghostCoef  = -1;    // ghost penalty gamma_g; <0 -> 10^-(p+1) (--study poisson)
    index_t kappaDense = 6000;  // dense-eigensolver dof threshold    (--study poisson)
    index_t kappaMaxIt = 500;   // power/inverse iteration cap        (--study poisson)
    std::string geo          = "identity"; // identity|affine|bubble background map
    real_t      eps           = 0.3;   // bubble amplitude, |eps|<0.5 (--geo bubble)
    std::string meshOverride;          // --mesh: overrides the case's mesh file
    std::string ghostMap      = "Gid"; // debug: Gid|G, the map inside the ghost dnk
    std::string dumpSystem;            // debug (--study poisson): raw system dump prefix
};

//----------------------------------------------------------------------------
// Small helpers.
//----------------------------------------------------------------------------

/// Bitwise (memcmp) equality of two gsMatrix<real_t>: same rows/cols, then
/// identical bytes over the whole (column-major, contiguous) buffer. Guards
/// size 0 (never calls memcmp on an empty buffer).
bool bitwiseEqual(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
{
    if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
    if (0 == A.size()) return true;
    return 0 == std::memcmp(A.data(), B.data(), (size_t)A.size()*sizeof(real_t));
}

/// Bitwise (memcmp) equality of two gsVector<real_t>.
bool bitwiseEqual(const gsVector<real_t> & a, const gsVector<real_t> & b)
{
    if (a.size() != b.size()) return false;
    if (0 == a.size()) return true;
    return 0 == std::memcmp(a.data(), b.data(), (size_t)a.size()*sizeof(real_t));
}

/// Counts, per node and per axis, occurrences of a coordinate landing
/// below the cell's own lower face or at-or-above its own upper face (the
/// half-open box [xL,xR) x [yB,yT) x [zN,zF) the background basis treats
/// as "belonging" to this element). INFO only: a rule node sitting exactly
/// on or past a knot line is a legitimate degeneracy the basis may
/// evaluate in a neighbouring knot span, not an error.
void countOutsideBox(const gsMatrix<real_t> & nodes, real_t xL, real_t xR, real_t yB, real_t yT,
                     real_t zN, real_t zF, long & below, long & atOrAbove)
{
    for (index_t c = 0; c != nodes.cols(); ++c)
    {
        const real_t x = nodes(0,c), y = nodes(1,c), z = nodes(2,c);
        if (x < xL) ++below; else if (x >= xR) ++atOrAbove;
        if (y < yB) ++below; else if (y >= yT) ++atOrAbove;
        if (z < zN) ++below; else if (z >= zF) ++atOrAbove;
    }
}

/// Case name -> ASCII gmsh MSH 4.1 file, gsFileManager-resolved by the
/// caller.
std::string meshFile(const std::string & caseName)
{
    if ("sphere"  == caseName) return "volumes/tetmesh_sphere.msh";
    if ("rotcube" == caseName) return "volumes/tetmesh_cube_rotated.msh";
    GISMO_ERROR("meshFile: unknown case '" << caseName << "'");
}

/// `--geo` string -> gsTetClip::BgMapKind. GISMO_ERRORs on anything else
/// (main() validates the string first, so this is only ever reached with a
/// value already known to be one of the three).
gsTetClip::BgMapKind geoKind(const std::string & geo)
{
    if ("identity" == geo) return gsTetClip::BgMapKind::Identity;
    if ("affine"   == geo) return gsTetClip::BgMapKind::Affine;
    if ("bubble"   == geo) return gsTetClip::BgMapKind::Bubble;
    GISMO_ERROR("geoKind: unknown --geo '" << geo << "'");
}

/// Builds the background space on grid \a g from gsTetClipSignDomain.h's/
/// gsImmersedPullback.h's shared helpers: \a mp is the background map
/// chosen by \a kind -- gsTetClip::identityBoxGeometry(g) (Identity, exactly
/// affine degree-1, physical == parametric everywhere on the box, not only
/// at grid lines), gsTetClip::affineBoxGeometry(g,A,b) with the default
/// (A,b) (gsTetClip::defaultAffineMap, Affine), or
/// gsTetClip::bubbleBoxGeometry(g,eps,c) with the default direction
/// (gsTetClip::defaultBubbleDirection, Bubble; \a eps unused otherwise).
/// \a mb wraps gsTetClip::backgroundBasis(g,p), the degree-p, n-elements-
/// per-direction basis whose breaks are bitwise \a g's own grid lines. The
/// two are independently built, so they need not agree on element count:
/// the assembler evaluates \a mp's geometry map at quadrature points
/// generated from \a mb's own elements, and spline evaluation only
/// requires the two bases to share a parameter range, which both take from
/// the same \a g.
void makeBackground(const gsTetClip::Grid3 & g, index_t p, gsTetClip::BgMapKind kind, real_t eps,
                    gsMultiPatch<real_t> & mp, gsMultiBasis<real_t> & mb)
{
    switch (kind)
    {
    case gsTetClip::BgMapKind::Identity:
        mp = gsTetClip::identityBoxGeometry(g);
        break;
    case gsTetClip::BgMapKind::Affine:
    {
        gsMatrix<real_t> A; gsVector<real_t> b;
        gsTetClip::defaultAffineMap(A, b);
        mp = gsTetClip::affineBoxGeometry(g, A, b);
        break;
    }
    case gsTetClip::BgMapKind::Bubble:
        mp = gsTetClip::bubbleBoxGeometry(g, eps, gsTetClip::defaultBubbleDirection());
        break;
    }
    mb = gsMultiBasis<real_t>(gsTetClip::backgroundBasis(g, p));
}

//----------------------------------------------------------------------------
// --study check, part (a): streaming/table bitwise equality.
//----------------------------------------------------------------------------

struct CheckACounters
{
    index_t cut = 0, full = 0, empty = 0;
    index_t statusMismatch = 0, volMismatch = 0, bdrMismatch = 0,
            normalMismatch = 0, tableMismatch = 0;
    real_t V_cut = 0;
    long belowLower = 0, atOrAboveUpper = 0;
};

/// Runs check (a): builds the reference tetClipQuadrature() \a Q and a
/// ClipStreamer \a S on the same (mesh,grid,p), and compares every cell's
/// status/volume-rule/boundary-rule/normal-field output of the streamer-
/// backed VolLookupRule/BdrLookupRule/BdrNormalField against \a Q, then
/// builds a VolCellTable/BdrCellTable from \a Q (AFTER those comparisons,
/// moving \a Q's per-cell rules in with give() to avoid holding two full
/// copies of the clip output at once) and compares the table-backed rules
/// against the streamer-backed ones. Returns \a S to the caller (kept
/// alive for check (b)/(c) at the last refinement level); \a Q and the
/// tables are local and freed on return.
CheckACounters checkA(const memory::shared_ptr<const gsTetClip::TetMesh> & M,
                      const gsTetClip::Grid3 & grid, index_t p,
                      memory::shared_ptr<gsTetClip::ClipStreamer> & Sout)
{
    CheckACounters c;

    gsTetClip::TetClipStats stats;
    gsTetClip::TetClipQuadrature Q = gsTetClip::tetClipQuadrature(*M, grid, p, stats);
    memory::shared_ptr<gsTetClip::ClipStreamer> S =
        memory::make_shared(new gsTetClip::ClipStreamer(M, grid, p));

    c.cut = S->numCut(); c.full = S->numFull(); c.empty = S->numEmpty();

    const gsVector<index_t> nG = gsVector<index_t>::Constant(3, p+1);
    gsTetClip::VolLookupRule volR(S->index(), S, nG);
    gsTetClip::BdrLookupRule bdrR(S->index(), S);
    gsTetClip::BdrNormalField nf(S->index(), S);

    gsTetClip::KahanSum vcut;
    const index_t n = grid.n;
    const std::vector<real_t> & X = S->index()->X;
    const std::vector<real_t> & Y = S->index()->Y;
    const std::vector<real_t> & Z = S->index()->Z;

    for (index_t k = 0; k != n; ++k)
    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const size_t id = S->index()->id(i,j,k);
        gsVector<real_t> lower(3), upper(3);
        lower << X[i], Y[j], Z[k];
        upper << X[i+1], Y[j+1], Z[k+1];

        const int qStatus = Q.status[id];
        const int sStatus = S->index()->status[id];
        if (sStatus != qStatus) ++c.statusMismatch;

        gsMatrix<real_t> vn; gsVector<real_t> vw;
        volR.mapTo(lower, upper, vn, vw);

        if (gsTetClip::Cut == qStatus)
        {
            if (0 == vn.cols() || !bitwiseEqual(vn, Q.vol[id].nodes) || !bitwiseEqual(vw, Q.vol[id].weights))
                ++c.volMismatch;
        }
        else if (gsTetClip::Full == qStatus)
        {
            gsMatrix<real_t> gn; gsVector<real_t> gw;
            gsGaussRule<real_t>(nG).mapTo(lower, upper, gn, gw);
            if (!bitwiseEqual(vn, gn) || !bitwiseEqual(vw, gw)) ++c.volMismatch;
        }
        else // Empty
        {
            if (0 != vn.cols()) ++c.volMismatch;
        }

        if (gsTetClip::Cut == sStatus)
        {
            long b = 0, a = 0;
            countOutsideBox(vn, X[i],X[i+1], Y[j],Y[j+1], Z[k],Z[k+1], b, a);
            c.belowLower += b; c.atOrAboveUpper += a;
            for (index_t kk = 0; kk != vw.size(); ++kk) vcut.add(vw[kk]);
        }

        gsMatrix<real_t> bn; gsVector<real_t> bw;
        bdrR.mapTo(lower, upper, bn, bw);
        if (!bitwiseEqual(bn, Q.bdr[id].nodes) || !bitwiseEqual(bw, Q.bdr[id].weights))
            ++c.bdrMismatch;
        {
            long b = 0, a = 0;
            countOutsideBox(bn, X[i],X[i+1], Y[j],Y[j+1], Z[k],Z[k+1], b, a);
            c.belowLower += b; c.atOrAboveUpper += a;
        }

        if (Q.bdr[id].nodes.cols() > 0)
        {
            gsMatrix<real_t> res;
            bool threw = false;
            try { nf.eval_into(Q.bdr[id].nodes, res); }
            catch (const std::runtime_error &) { threw = true; }
            if (threw || !bitwiseEqual(res, Q.bdr[id].normals)) ++c.normalMismatch;
        }
    }
    c.V_cut = vcut.value();

    // Table path: fill AFTER the Q-vs-streamer comparisons above, moving
    // Q's per-cell rules out rather than copying them (Q.vol: Cut cells
    // only; Q.bdr: every cell).
    const size_t nCells = (size_t)n*(size_t)n*(size_t)n;
    memory::shared_ptr<gsTetClip::VolCellTable> volT = memory::make_shared(new gsTetClip::VolCellTable());
    memory::shared_ptr<gsTetClip::BdrCellTable> bdrT = memory::make_shared(new gsTetClip::BdrCellTable());
    volT->cell.resize(nCells);
    bdrT->cell.resize(nCells);
    for (size_t id = 0; id != nCells; ++id)
    {
        if (gsTetClip::Cut == Q.status[id]) volT->cell[id] = give(Q.vol[id]);
        bdrT->cell[id] = give(Q.bdr[id]);
    }

    gsTetClip::VolLookupRule volRT(S->index(), volT, nG);
    gsTetClip::BdrLookupRule bdrRT(S->index(), bdrT);

    for (index_t k = 0; k != n; ++k)
    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        gsVector<real_t> lower(3), upper(3);
        lower << X[i], Y[j], Z[k];
        upper << X[i+1], Y[j+1], Z[k+1];

        gsMatrix<real_t> vn1, vn2; gsVector<real_t> vw1, vw2;
        volR.mapTo(lower, upper, vn1, vw1);
        volRT.mapTo(lower, upper, vn2, vw2);
        if (!bitwiseEqual(vn1,vn2) || !bitwiseEqual(vw1,vw2)) ++c.tableMismatch;

        gsMatrix<real_t> bn1, bn2; gsVector<real_t> bw1, bw2;
        bdrR.mapTo(lower, upper, bn1, bw1);
        bdrRT.mapTo(lower, upper, bn2, bw2);
        if (!bitwiseEqual(bn1,bn2) || !bitwiseEqual(bw1,bw2)) ++c.tableMismatch;
    }

    Sout = S;
    return c;
}

//----------------------------------------------------------------------------
// --study check, part (b): scopes and ghost bitwise equality.
//----------------------------------------------------------------------------

struct ScopeChecks
{
    bool stateAssembler     = false; // check 1
    bool partitionOfUnity   = false; // check 2
    bool evaluatorIntegral  = false; // check 3
    bool quDimPresentBranch = false; // check 4
    bool nestedThrow        = false; // check 5
    bool boundaryField      = false; // check 6

    bool allPass() const
    {
        return stateAssembler && partitionOfUnity && evaluatorIntegral &&
               quDimPresentBranch && nestedThrow && boundaryField;
    }
};

/// Runs one ghost-penalty assembly (pattern: poisson2_ghost_penalty_example.cpp,
/// solveOne()). With \a withScopes == false, no quadrature factory is ever
/// installed (plain option-driven Gauss on every background element, the
/// SAME pattern path as the scoped volume step, which is why the two
/// resulting matrices are bitwise comparable). With \a withScopes == true,
/// every checkable behaviour of QuadratureScope/VolumeQuadratureScope/
/// BoundaryQuadratureScope/BdrNormalField is probed and folded into
/// \a checks (non-null only on the withScopes run). \a nf must already be
/// registered on the SAME \a S as \a volFactory/\a bdrFactory build from.
gsSparseMatrix<real_t> ghostRun(bool withScopes, const gsMultiPatch<real_t> & mp,
                                gsMultiBasis<real_t> & mb, index_t p,
                                memory::shared_ptr<gsTetClip::ClipStreamer> S,
                                memory::shared_ptr<gsImplicitTrimmedDomain<3,real_t> > trDom,
                                gsTetClip::BdrNormalField & nf, real_t V_ref, ScopeChecks * checks)
{
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::space       space;

    gsExprAssembler<real_t> A(1,1);
    geometryMap G = A.getMap(mp);
    space u = A.getSpace(mb);
    gsBoundaryConditions<real_t> bcNone;
    u.setup(bcNone, dirichlet::none, 0);
    A.setIntegrationElements(mb);
    A.initSystem();

    const gsVector<index_t> nG = gsVector<index_t>::Constant(3, p+1);

    // VOLUME step, both runs: same pattern path (every background element,
    // full grid), differing only in which quadrature produced the points.
    if (withScopes)
    {
        gsExprAssembler<real_t>::QuadratureFactory volFactory =
            gsTetClip::makeVolLookupFactory(S->index(), S, nG);
        {
            gsTetClip::VolumeQuadratureScope<gsExprAssembler<real_t> > s(A, volFactory);
            if (checks)
                checks->stateAssembler = A.hasCustomQuadrature() && (-1 == A.options().getInt("quDim"));
            A.assemble(u * meas(G));
        }
        if (checks)
            checks->stateAssembler = checks->stateAssembler &&
                !A.hasCustomQuadrature() && !gsTetClip::hasIntOption(A.options(), "quDim");
    }
    else
    {
        A.assemble(u * meas(G));
    }

    if (withScopes)
    {
        // Partition of unity: the B-splines sum to 1 and every dof is
        // free, so Sum_i rhs_i = Int 1 dV under the rule. This validates
        // keying, weights and thread-safety of the rules under `#pragma
        // omp parallel`, NOT per-entry scatter into the rhs vector.
        if (checks)
        {
            gsTetClip::KahanSum rhsSum;
            const gsMatrix<real_t> & rhs = A.rhs();
            for (index_t cc = 0; cc != rhs.cols(); ++cc)
                for (index_t r = 0; r != rhs.rows(); ++r)
                    rhsSum.add(rhs(r,cc));
            const real_t rel = math::abs(rhsSum.value() - V_ref) / math::abs(V_ref);
            checks->partitionOfUnity = (rel <= 1e-12);
        }

        // Evaluator: shares A's exprData (hence the current integration
        // domain) but NOT its installed factory, so it needs its own
        // scope.
        {
            gsExprEvaluator<real_t> ev(A);
            gsExprAssembler<real_t>::QuadratureFactory volFactory2 =
                gsTetClip::makeVolLookupFactory(S->index(), S, nG);
            real_t integralVal = 0;
            {
                gsTetClip::VolumeQuadratureScope<gsExprEvaluator<real_t> > se(ev, volFactory2);
                integralVal = ev.integral(meas(G));
            }
            if (checks)
            {
                const real_t rel = math::abs(integralVal - V_ref) / math::abs(V_ref);
                checks->evaluatorIntegral = (rel <= 1e-12) &&
                    !ev.hasCustomQuadrature() && !gsTetClip::hasIntOption(ev.options(), "quDim");
            }
        }

        // quDim-present branch: a pre-existing int option must be restored
        // exactly, both for the volume (-1) and boundary (2) scope.
        {
            A.options().addInt("quDim", "test", 7);
            bool ok = true;
            {
                gsExprAssembler<real_t>::QuadratureFactory vf =
                    gsTetClip::makeVolLookupFactory(S->index(), S, nG);
                gsTetClip::VolumeQuadratureScope<gsExprAssembler<real_t> > sv(A, vf);
            }
            ok = ok && (7 == A.options().getInt("quDim")) && !A.hasCustomQuadrature();
            {
                gsExprAssembler<real_t>::QuadratureFactory bf =
                    gsTetClip::makeBdrLookupFactory(S->index(), S);
                bool quDimWas2 = false;
                {
                    gsTetClip::BoundaryQuadratureScope<gsExprAssembler<real_t> > sb(A, bf);
                    quDimWas2 = (2 == A.options().getInt("quDim"));
                }
                ok = ok && quDimWas2 && (7 == A.options().getInt("quDim")) && !A.hasCustomQuadrature();
            }
            A.options().remove("quDim");
            if (checks) checks->quDimPresentBranch = ok;
        }

        // Nested-scope throw: entering a second scope while one is still
        // installed must throw, and the first scope must still unwind
        // correctly (normal C++ stack unwinding of the try block).
        {
            bool threw = false;
            try
            {
                gsExprAssembler<real_t>::QuadratureFactory vf =
                    gsTetClip::makeVolLookupFactory(S->index(), S, nG);
                gsTetClip::VolumeQuadratureScope<gsExprAssembler<real_t> > s1(A, vf);
                gsExprAssembler<real_t>::QuadratureFactory bf =
                    gsTetClip::makeBdrLookupFactory(S->index(), S);
                gsTetClip::BoundaryQuadratureScope<gsExprAssembler<real_t> > s2(A, bf);
            }
            catch (const std::runtime_error &) { threw = true; }
            if (checks)
                checks->nestedThrow = threw && !A.hasCustomQuadrature() && !gsTetClip::hasIntOption(A.options(), "quDim");
        }
    }

    A.setIntegrationDomain(trDom);

    if (withScopes)
    {
        std::vector<patchSide> bdrImm(1, patchSide(0, boundary::none));
        auto nfC = A.getCoeff(nf); // NO geometry argument: evaluated at the rule's own nodes
        const gsMatrix<real_t> r0 = A.rhs();
        bool ok = true;
        {
            gsExprAssembler<real_t>::QuadratureFactory bf =
                gsTetClip::makeBdrLookupFactory(S->index(), S);
            gsTetClip::BoundaryQuadratureScope<gsExprAssembler<real_t> > sb(A, bf);
            try
            {
                A.assembleBdr(bdrImm, u);
                const gsMatrix<real_t> r1 = A.rhs();
                A.assembleBdr(bdrImm, u * nfC.norm());
                const gsMatrix<real_t> r2 = A.rhs();

                const gsMatrix<real_t> d1 = r1 - r0;
                gsTetClip::KahanSum d1sum;
                for (index_t r = 0; r != d1.rows(); ++r) d1sum.add(d1(r,0));
                const gsMatrix<real_t> resid = r2 - r1 - d1;
                ok = (d1sum.value() > 0) && (resid.norm() <= 1e-12*d1.norm());
            }
            catch (const std::runtime_error &) { ok = false; }
        }
        if (checks)
            checks->boundaryField = ok && !A.hasCustomQuadrature() && !gsTetClip::hasIntOption(A.options(), "quDim");
    }

    auto dJ = dnk(u.jump(), G.left(), p);
    A.computePatternGhost(dJ * dJ.tr());
    A.assembleGhost(dJ * dJ.tr());

    gsSparseMatrix<real_t> Mg = A.matrix();
    Mg.makeCompressed();
    return Mg;
}

//----------------------------------------------------------------------------
// --study check, part (c): normal-field throw behaviour.
//----------------------------------------------------------------------------

/// Picks the Cut cell with the most boundary columns, checks the positive
/// case (exact nodes -> exact stored normals, no throw) and four negative
/// cases that must each throw std::runtime_error. Prints CHECK(c).
bool checkC(const std::string & caseName, index_t n, memory::shared_ptr<gsTetClip::ClipStreamer> S)
{
    const gsTetClip::CellIndex & idx = *S->index();
    const index_t nCells1D = idx.grid.n;
    const size_t nCells = (size_t)nCells1D*(size_t)nCells1D*(size_t)nCells1D;

    size_t bestId = nCells; index_t bestCols = -1;
    gsMatrix<real_t> bestNodes, bestNormals;
    for (size_t id = 0; id != nCells; ++id)
    {
        if (gsTetClip::Cut != idx.status[id]) continue;
        gsMatrix<real_t> nn, nrm; gsVector<real_t> ww;
        S->bdrRule(id, nn, ww, nrm);
        if ((index_t)nn.cols() > bestCols)
        { bestCols = (index_t)nn.cols(); bestId = id; bestNodes = nn; bestNormals = nrm; }
    }
    GISMO_ENSURE(nCells != bestId && bestCols >= 2,
                "checkC: no Cut cell with >= 2 boundary nodes found for --case " << caseName);

    gsTetClip::BdrNormalField nf(S->index(), S);

    gsInfo << "note: CHECK(c)'s negative cases below are expected to print stderr 'Error' lines "
             "from GISMO_ENSURE/GISMO_ERROR; this is normal, not a failure.  INFO\n";

    gsMatrix<real_t> res;
    bool positiveOk = false;
    try
    {
        nf.eval_into(bestNodes, res);
        positiveOk = bitwiseEqual(res, bestNormals);
    }
    catch (const std::runtime_error &) { positiveOk = false; }

    index_t throws = 0;

    { // 1: nodes with one coordinate one ULP off.
        gsMatrix<real_t> P = bestNodes;
        P(0,0) = std::nextafter(P(0,0), std::numeric_limits<real_t>::infinity());
        try { gsMatrix<real_t> r; nf.eval_into(P, r); } catch (const std::runtime_error &) { ++throws; }
    }
    { // 2: one fewer column.
        gsMatrix<real_t> P = bestNodes.leftCols(bestNodes.cols()-1);
        try { gsMatrix<real_t> r; nf.eval_into(P, r); } catch (const std::runtime_error &) { ++throws; }
    }
    { // 3: columns 0 and 1 swapped (order).
        gsMatrix<real_t> P = bestNodes;
        P.col(0).swap(P.col(1));
        try { gsMatrix<real_t> r; nf.eval_into(P, r); } catch (const std::runtime_error &) { ++throws; }
    }
    { // 4: a single point at the centre of cell id 0 (corner (-1,-1,-1)),
      // far from either mesh, which therefore has no boundary rule.
        gsMatrix<real_t> P(3,1);
        P(0,0) = 0.5*(idx.X[0]+idx.X[1]);
        P(1,0) = 0.5*(idx.Y[0]+idx.Y[1]);
        P(2,0) = 0.5*(idx.Z[0]+idx.Z[1]);
        try { gsMatrix<real_t> r; nf.eval_into(P, r); } catch (const std::runtime_error &) { ++throws; }
    }

    const bool pass = positiveOk && (4 == throws);
    gsInfo << "CHECK(c) case=" << caseName << " n=" << n << " positive=" << (positiveOk ? "ok" : "FAIL")
          << " throws=" << throws << "/4  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

//----------------------------------------------------------------------------
// --study check, part (d): TetClipSignDomain classification.
//----------------------------------------------------------------------------

/// The trimmed-domain sign of a gsTetClip::CellStatus value
/// (gsTetClipSignDomain.h's file-header table), duplicated here rather than
/// called through TetClipSignDomain: the checks below need a reference sign
/// array BEFORE (check A/D/F) and INDEPENDENTLY of (check C) constructing
/// the domain under test, so the mapping cannot be read off the class
/// itself without making every check trust the very code it is meant to
/// verify.
short_t toDomainSign(int status)
{
    switch (status)
    {
    case gsTetClip::Full:  return -1;
    case gsTetClip::Cut:   return  0;
    case gsTetClip::Empty: return  1;
    default: GISMO_ERROR("toDomainSign: status=" << status << " is not a gsTetClip::CellStatus value.");
    }
}

/// Walks one sign family of \a dom (SignOp/\a familySign matched pairs:
/// InteriorSign/-1, BoundarySign/0, ExteriorSign/1), recovers each visited
/// element's background cell id from its midpoint (the SAME colOf lookup
/// every gsTetClip rule uses), and checks it against the reference status
/// array \a idx.status both via the domain's OWN cached sign (it.sign())
/// and via a fresh toDomainSign() conversion -- the two must always agree
/// with each other and with \a familySign, since it.sign() is exactly what
/// selected this element into this family in the first place. Every
/// visited cell id increments \a visited[id], so the caller can check the
/// three families tile the grid exactly (each id visited once).
template<typename SignOp>
void walkFamilySign(gsTetClip::TetClipSignDomain & dom, const gsTetClip::CellIndex & idx,
                    short_t familySign, std::vector<index_t> & visited, index_t & mismatches)
{
    gsDomain<real_t>::iterator e = dom.end<SignOp>();
    for (gsDomain<real_t>::iterator it = dom.begin<SignOp>(); it < e; ++it)
    {
        const gsVector<real_t> mid = 0.5*(it.lowerCorner() + it.upperCorner());
        const index_t i = gsTetClip::colOf(mid[0], idx.X);
        const index_t j = gsTetClip::colOf(mid[1], idx.Y);
        const index_t k = gsTetClip::colOf(mid[2], idx.Z);
        const size_t id = idx.id(i,j,k);
        if (toDomainSign(idx.status[id]) != familySign || it.sign() != familySign) ++mismatches;
        ++visited[(index_t)id];
    }
}

/// Per-(case,r) outcome of every classification sub-check (A-H: reference-
/// bitwise equality, leaf structure, sign-family tiling, element/ghost/
/// skeleton counts, the boundary-piece invariant, the multi-cell guard, and
/// the identity geometry). Check A only runs at r <= 1 (aRun); when it does
/// not run, aPass stays true so it never gates allPass().
struct SignDomainChecks
{
    bool aRun = false, aPass = true; index_t aMismatches = 0;
    bool bPass = true;
    bool cPass = true; index_t cMismatches = 0;
    index_t dFull = 0, dCut = 0, dEmpty = 0; bool dPass = true;
    bool ePass = true; index_t eBad = 0;
    size_t fGhostRef = 0, fSkelRef = 0, fGhostDom = 0, fSkelDom = 0; bool fPass = true;
    bool gPass = true;
    bool hPass = true; real_t hMaxErr = 0;

    bool allPass() const
    { return aPass && bPass && cPass && dPass && ePass && fPass && gPass && hPass; }
};

/// Runs checks A-H (see SignDomainChecks) for one (mesh, grid, p)
/// at refinement level \a r. Builds a gsTetClip::ClipStreamer \a S as the
/// status source: ClipStreamer's own constructor already runs one full clip
/// to classify every cell, at the SAME asymptotic cost as
/// gsTetClip::tetClipQuadrature but without storing an n^3 array of
/// per-cell rules (see gsImmersedLookupRule.h's own doxygen on
/// ClipStreamer). Check A additionally builds the full reference
/// gsTetClip::tetClipQuadrature() in an inner scope, freed before
/// returning, and only at \a r <= 1: at n=8, tetClipQuadrature keeps every
/// clipped node of every cell live at once, several GB for the sphere case.
SignDomainChecks checkSignDomainOne(const memory::shared_ptr<const gsTetClip::TetMesh> & Mshared,
                                    const std::string & caseName, const gsTetClip::Grid3 & grid,
                                    index_t p, index_t r)
{
    SignDomainChecks c;
    const index_t n = grid.n;
    const size_t nCells = (size_t)n*(size_t)n*(size_t)n;

    gsTetClip::ClipStreamer S(Mshared, grid, p);
    const gsTetClip::CellIndex & idx = *S.index();

    // Check A: bitwise against the reference tetClipQuadrature, r <= 1 only.
    // Piece presence (nVolPieces[id]>0, nBdrPieces[id]>0) is read from
    // S.volRule()/S.bdrRule()'s own returned column count, not derived from
    // idx.status: deriving it from status would make volPresentOk trivially
    // true whenever statusOk already holds (both status decisions are
    // "touched at all" tests), which is exactly the kind of check that
    // cannot fail and so proves nothing on its own.
    if (r <= 1)
    {
        c.aRun = true;
        gsTetClip::TetClipStats stats;
        gsTetClip::TetClipQuadrature Q = gsTetClip::tetClipQuadrature(*Mshared, grid, p, stats);
        index_t listed = 0;
        for (size_t id = 0; id != nCells; ++id)
        {
            const bool statusOk = (Q.status[id] == idx.status[id]);

            gsMatrix<real_t> svn; gsVector<real_t> svw;
            S.volRule(id, svn, svw);
            const bool qVolPresent = Q.vol[id].weights.size() > 0;
            const bool sVolPresent = svw.size() > 0;
            const bool volPresentOk = (qVolPresent == sVolPresent);
            gsTetClip::CellRule3 tmp; tmp.weights = svw;
            const real_t sVol = gsTetClip::cellVolume(tmp);
            const real_t qVol = gsTetClip::cellVolume(Q.vol[id]);
            const bool volOk = (qVol == sVol);

            gsMatrix<real_t> sbn; gsVector<real_t> sbw; gsMatrix<real_t> sbnrm;
            S.bdrRule(id, sbn, sbw, sbnrm);
            const bool qBdrPresent = Q.bdr[id].weights.size() > 0;
            const bool sBdrPresent = sbw.size() > 0;
            const bool bdrPresentOk = (qBdrPresent == sBdrPresent);

            if (!(statusOk && volPresentOk && bdrPresentOk && volOk))
            {
                ++c.aMismatches;
                if (listed < 10)
                {
                    gsInfo << "  CHECK(d.A) mismatch case=" << caseName << " n=" << n << " id=" << id
                          << " Qstatus=" << Q.status[id] << " Sstatus=" << idx.status[id]
                          << " qVolPresent=" << qVolPresent << " sVolPresent=" << sVolPresent
                          << " qBdrPresent=" << qBdrPresent << " sBdrPresent=" << sBdrPresent
                          << " qVol=" << gsTetClip::fmtSci(qVol) << " sVol=" << gsTetClip::fmtSci(sVol) << "\n";
                    ++listed;
                }
            }
        }
        c.aPass = (0 == c.aMismatches);
    }

    // The domain under test.
    gsTetClip::TetClipSignDomain dom(grid, idx.status, p);

    // Check B: construction and leaf structure.
    c.bPass = (1 == dom.numLevels());
    {
        auto it = dom.tree().beginLeafIterator();
        size_t nLeaves = 0;
        bool leavesOk = true;
        while (it.good())
        {
            leavesOk = leavesOk && (0 == it.data().level()) &&
                (1 == (it.data().upperCorner() - it.data().lowerCorner()).prod());
            ++nLeaves;
            it.next();
        }
        c.bPass = c.bPass && leavesOk && (nLeaves == nCells);
    }
    for (short_t j = 0; j != 3; ++j)
    {
        const std::vector<real_t> & br  = dom.breaks(0, j);
        const std::vector<real_t> & ref = (0==j) ? idx.X : (1==j) ? idx.Y : idx.Z;
        bool same = (br.size() == ref.size());
        for (size_t kk = 0; same && kk != br.size(); ++kk) same = (br[kk] == ref[kk]);
        c.bPass = c.bPass && same;
    }

    // Check C: zero sign mismatches, families tile the grid exactly.
    {
        std::vector<index_t> visited(nCells, 0);
        index_t mismatches = 0;
        walkFamilySign<InteriorSign>(dom, idx, -1, visited, mismatches);
        walkFamilySign<BoundarySign>(dom, idx,  0, visited, mismatches);
        walkFamilySign<ExteriorSign>(dom, idx,  1, visited, mismatches);
        bool tiled = true;
        for (size_t id = 0; id != nCells; ++id) tiled = tiled && (1 == visited[id]);
        c.cMismatches = mismatches;
        c.cPass = (0 == mismatches) && tiled;
    }

    // Check D: counts.
    {
        index_t full = 0, cut = 0, empty = 0;
        for (size_t id = 0; id != nCells; ++id)
        {
            if      (gsTetClip::Full  == idx.status[id]) ++full;
            else if (gsTetClip::Cut   == idx.status[id]) ++cut;
            else                                         ++empty;
        }
        c.dFull = full; c.dCut = cut; c.dEmpty = empty;
        c.dPass = (dom.numElements<InteriorSign>() == (size_t)full) &&
                  (dom.numElements<BoundarySign>() == (size_t)cut) &&
                  (dom.numElements<ExteriorSign>() == (size_t)empty);
    }

    // Check E: boundary-piece invariant. A boundary piece in a Full/Empty
    // cell would be silently dropped: beginBdr() ignores its side argument
    // and walks only sign-0 leaves (gsTrimmedDomain.h:288-291), so any
    // Nitsche integral built on this domain would never see it.
    {
        index_t bad = 0, listed = 0;
        const real_t box = grid.h*grid.h*grid.h;
        for (size_t id = 0; id != nCells; ++id)
        {
            gsMatrix<real_t> bn; gsVector<real_t> bw; gsMatrix<real_t> bnrm;
            S.bdrRule(id, bn, bw, bnrm);
            const bool bdrPresent = bw.size() > 0;

            gsMatrix<real_t> vn; gsVector<real_t> vw;
            S.volRule(id, vn, vw);
            const bool volPresent = vw.size() > 0;

            const bool cut = (gsTetClip::Cut == idx.status[id]);
            const bool ok = (!bdrPresent || cut) && (!cut || volPresent);
            if (!ok)
            {
                ++bad;
                if (listed < 10)
                {
                    gsTetClip::CellRule3 tmp; tmp.weights = vw;
                    const real_t vol = gsTetClip::cellVolume(tmp);
                    gsInfo << "  CHECK(d.E) violation case=" << caseName << " n=" << n << " id=" << id
                          << " status=" << idx.status[id] << " bdrPresent=" << bdrPresent
                          << " volPresent=" << volPresent
                          << " |vol-h^3|/h^3=" << gsTetClip::fmtSci(math::abs(vol-box)/box) << "\n";
                    ++listed;
                }
            }
        }
        c.eBad = bad;
        c.ePass = (0 == bad);
    }

    // Check F: ghost/skeleton counts, recomputed directly from
    // toDomainSign(idx.status) over every axis-neighbour interior face
    // (SkeletonFace/GhostFace predicates, gsTrimmedDomain.h:154-169), and
    // compared to dom.numGhostFaces()/numSkeletonFaces() -- called
    // single-threaded here, before anything parallel could touch the
    // lazily built, non-thread-safe face sign grid (gsTrimmedDomain.h:
    // 237-240).
    {
        size_t ghostRef = 0, skelRef = 0;
        for (index_t k = 0; k != n; ++k)
        for (index_t j = 0; j != n; ++j)
        for (index_t i = 0; i != n; ++i)
        {
            const short_t sHere = toDomainSign(idx.status[idx.id(i,j,k)]);
            if (i+1 < n)
            {
                const short_t sNext = toDomainSign(idx.status[idx.id(i+1,j,k)]);
                if (sHere<=0 && sNext<=0) { ++skelRef; if (0==sHere || 0==sNext) ++ghostRef; }
            }
            if (j+1 < n)
            {
                const short_t sNext = toDomainSign(idx.status[idx.id(i,j+1,k)]);
                if (sHere<=0 && sNext<=0) { ++skelRef; if (0==sHere || 0==sNext) ++ghostRef; }
            }
            if (k+1 < n)
            {
                const short_t sNext = toDomainSign(idx.status[idx.id(i,j,k+1)]);
                if (sHere<=0 && sNext<=0) { ++skelRef; if (0==sHere || 0==sNext) ++ghostRef; }
            }
        }
        c.fGhostRef = ghostRef; c.fSkelRef = skelRef;
        c.fGhostDom = dom.numGhostFaces(); c.fSkelDom = dom.numSkeletonFaces();
        c.fPass = (ghostRef == c.fGhostDom) && (skelRef == c.fSkelDom);
    }

    // Check G: the multi-cell guard fires (points from cell (0,0,0) and
    // cell (n-1,n-1,n-1), a direct call outside any OpenMP region), and the
    // single-cell positive control (the 8 corners of cell 0) returns the
    // constant cell sign.
    {
        gsMatrix<real_t> P(3,2);
        P(0,0) = 0.5*(idx.X[0]+idx.X[1]);     P(1,0) = 0.5*(idx.Y[0]+idx.Y[1]);     P(2,0) = 0.5*(idx.Z[0]+idx.Z[1]);
        P(0,1) = 0.5*(idx.X[n-1]+idx.X[n]);   P(1,1) = 0.5*(idx.Y[n-1]+idx.Y[n]);   P(2,1) = 0.5*(idx.Z[n-1]+idx.Z[n]);
        bool threw = false;
        try { dom.sign(P); } catch (const std::runtime_error &) { threw = true; }
        gsInfo << "  CHECK(d.G) note: the GISMO_ENSURE 'Error' line just above (if present) is "
                 "the expected guard firing, not a failure.  INFO\n";

        gsMatrix<real_t> corners(3,8);
        index_t col = 0;
        for (int dz = 0; dz != 2; ++dz)
        for (int dy = 0; dy != 2; ++dy)
        for (int dx = 0; dx != 2; ++dx)
        { corners(0,col) = idx.X[dx]; corners(1,col) = idx.Y[dy]; corners(2,col) = idx.Z[dz]; ++col; }
        bool posOk = true;
        try
        {
            const gsVector<short_t> s = dom.sign(corners);
            const short_t expect = toDomainSign(idx.status[0]);
            for (index_t kk = 0; kk != s.size(); ++kk) posOk = posOk && (s[kk] == expect);
        }
        catch (const std::runtime_error &) { posOk = false; }

        c.gPass = threw && posOk;
    }

    // Check H: the identity geometry, evaluated at every grid vertex.
    {
        gsMultiPatch<real_t> G = gsTetClip::identityBoxGeometry(grid);
        gsMatrix<real_t> V(3, (index_t)((n+1)*(n+1)*(n+1)));
        index_t col = 0;
        for (index_t kk = 0; kk <= n; ++kk)
        for (index_t jj = 0; jj <= n; ++jj)
        for (index_t ii = 0; ii <= n; ++ii)
        { V(0,col) = idx.X[ii]; V(1,col) = idx.Y[jj]; V(2,col) = idx.Z[kk]; ++col; }
        gsMatrix<real_t> Ev;
        G.patch(0).eval_into(V, Ev);
        const real_t tol = 64*std::numeric_limits<real_t>::epsilon()*gsTetClip::boxCornerAbsMax(grid);
        c.hMaxErr = (Ev - V).array().abs().maxCoeff();
        c.hPass = (c.hMaxErr <= tol);
    }

    return c;
}

/// Runs checkSignDomainOne() for case in {sphere, rotcube} x r in {0,1,2,3},
/// INDEPENDENTLY of the driver's own --case/-r options (this section always
/// covers all 8 (case,r) pairs). Prints one summary line per pair and
/// returns whether every one of them passed.
bool checkSignDomain(index_t n0, index_t p)
{
    bool allPass = true;
    const std::vector<std::string> cases = { "sphere", "rotcube" };

    for (const std::string & caseName : cases)
    {
        const std::string resolvedFile = gsFileManager::find(meshFile(caseName));
        GISMO_ENSURE(!resolvedFile.empty(), "checkSignDomain: mesh file for --case " << caseName
                    << " not found.");
        const memory::shared_ptr<const gsTetClip::TetMesh> Mshared =
            memory::make_shared(new gsTetClip::TetMesh(gsTetClip::readMsh41(resolvedFile)));

        for (index_t r = 0; r <= 3; ++r)
        {
            const index_t n = n0 << r;
            gsTetClip::Grid3 grid;
            grid.x0 = grid.y0 = grid.z0 = -1.0;
            grid.h  = 2.0/(real_t)n;
            grid.n  = n;
            GISMO_ENSURE(Mshared->lo[0] > grid.x0 && Mshared->hi[0] < grid.x0 + grid.n*grid.h &&
                        Mshared->lo[1] > grid.y0 && Mshared->hi[1] < grid.y0 + grid.n*grid.h &&
                        Mshared->lo[2] > grid.z0 && Mshared->hi[2] < grid.z0 + grid.n*grid.h,
                        "checkSignDomain: mesh bbox does not lie strictly inside the background "
                        "box for --case " << caseName << " at n=" << n);

            gsStopwatch sw;
            const SignDomainChecks c = checkSignDomainOne(Mshared, caseName, grid, p, r);
            const real_t elapsed = sw.stop();

            const bool pass = c.allPass();
            gsInfo << "sign-domain case=" << caseName << " r=" << r << " n=" << n
                  << " full=" << c.dFull << " cut=" << c.dCut << " empty=" << c.dEmpty
                  << " mismatches=" << c.cMismatches
                  << " ghost=" << c.fGhostDom << "/" << c.fGhostRef
                  << " skel=" << c.fSkelDom << "/" << c.fSkelRef
                  << " A=" << (!c.aRun ? "skip" : (c.aPass ? "PASS" : "FAIL"))
                  << " B=" << (c.bPass ? "PASS" : "FAIL") << " C=" << (c.cPass ? "PASS" : "FAIL")
                  << " D=" << (c.dPass ? "PASS" : "FAIL") << " E=" << (c.ePass ? "PASS" : "FAIL")
                  << " F=" << (c.fPass ? "PASS" : "FAIL") << " G=" << (c.gPass ? "PASS" : "FAIL")
                  << " H=" << (c.hPass ? "PASS" : "FAIL")
                  << " time=" << gsTetClip::fmtSci(elapsed) << "s  " << (pass ? "PASS" : "FAIL") << "\n";
            allPass = allPass && pass;
        }
    }

    return allPass;
}

//----------------------------------------------------------------------------
// --study check driver.
//----------------------------------------------------------------------------

bool runCheckStudy(const Config & cfg)
{
    std::vector<std::string> cases;
    if ("all" == cfg.caseName) { cases.push_back("sphere"); cases.push_back("rotcube"); }
    else cases.push_back(cfg.caseName);

    bool allPass = true;

    for (const std::string & caseName : cases)
    {
        const std::string resolvedFile = gsFileManager::find(meshFile(caseName));
        GISMO_ENSURE(!resolvedFile.empty(), "runCheckStudy: mesh file for --case " << caseName
                    << " not found.");
        const memory::shared_ptr<const gsTetClip::TetMesh> M =
            memory::make_shared(new gsTetClip::TetMesh(gsTetClip::readMsh41(resolvedFile)));

        for (index_t r = 0; r <= cfg.rMax; ++r)
        {
            const index_t n = cfg.n0 << r;
            gsTetClip::Grid3 grid;
            grid.x0 = grid.y0 = grid.z0 = -1.0;
            grid.h = 2.0/(real_t)n;
            grid.n = n;
            GISMO_ENSURE(M->lo[0] > grid.x0 && M->hi[0] < grid.x0 + grid.n*grid.h &&
                        M->lo[1] > grid.y0 && M->hi[1] < grid.y0 + grid.n*grid.h &&
                        M->lo[2] > grid.z0 && M->hi[2] < grid.z0 + grid.n*grid.h,
                        "runCheckStudy: mesh bbox does not lie strictly inside the background "
                        "box for --case " << caseName << " at n=" << n);

            memory::shared_ptr<gsTetClip::ClipStreamer> S;
            const CheckACounters ca = checkA(M, grid, cfg.p, S);

            const bool aPass = (0 == ca.statusMismatch && 0 == ca.volMismatch && 0 == ca.bdrMismatch &&
                               0 == ca.normalMismatch && 0 == ca.tableMismatch && ca.cut > 0);
            gsInfo << "CHECK(a) case=" << caseName << " n=" << n << " cut=" << ca.cut
                  << " full=" << ca.full << " empty=" << ca.empty
                  << " statusMismatch=" << ca.statusMismatch << " volMismatch=" << ca.volMismatch
                  << " bdrMismatch=" << ca.bdrMismatch << " normalMismatch=" << ca.normalMismatch
                  << " tableMismatch=" << ca.tableMismatch << "  " << (aPass ? "PASS" : "FAIL") << "\n";
            gsInfo << "INFO case=" << caseName << " n=" << n << " V_cut=" << gsTetClip::fmtSci(ca.V_cut)
                  << " nodesOutsideBoxBelowLower=" << ca.belowLower
                  << " nodesOutsideBoxAtOrAboveUpper=" << ca.atOrAboveUpper << "  INFO\n";
            allPass = allPass && aPass;

            if (r == cfg.rMax)
            {
                gsMultiPatch<real_t> mp; gsMultiBasis<real_t> mb;
                makeBackground(grid, cfg.p, gsTetClip::BgMapKind::Identity, 0.0, mp, mb);
                gsTensorBSplineBasis<3,real_t> * tbs =
                    dynamic_cast<gsTensorBSplineBasis<3,real_t> *>(&mb.basis(0));
                GISMO_ENSURE(tbs, "runCheckStudy: background basis is not a tensor B-spline basis.");

                gsFunctionExpr<real_t> phi("(x-0.03)^2+(y+0.02)^2+(z-0.01)^2-0.3025", 3);
                memory::shared_ptr<gsImplicitTrimmedDomain<3,real_t> > trDom =
                    memory::make_shared(new gsImplicitTrimmedDomain<3,real_t>(phi, *tbs));
                GISMO_ENSURE(1 == trDom->numLevels(), "runCheckStudy: gsImplicitTrimmedDomain "
                            "requires a single kd-tree level.");
                // Warm the (lazily built, non-thread-safe) face sign grid
                // once, single-threaded, before any OpenMP region touches it.
                const index_t nGhost = (index_t)trDom->numGhostFaces();
                const index_t nCutA  = (index_t)trDom->numElementsBdr(boundary::none);
                gsInfo << "INFO case=" << caseName << " n=" << n << " nGhost=" << nGhost
                      << " nCutA=" << nCutA << "  INFO\n";

                gsTetClip::BdrNormalField nf(S->index(), S);
                const real_t V_ref = ca.V_cut + (real_t)ca.full*grid.h*grid.h*grid.h;

                ScopeChecks checks;
                const gsSparseMatrix<real_t> M2 = ghostRun(true,  mp, mb, cfg.p, S, trDom, nf, V_ref, &checks);
                const gsSparseMatrix<real_t> M1 = ghostRun(false, mp, mb, cfg.p, S, trDom, nf, V_ref, nullptr);

                bool sameShape = (M1.rows() == M2.rows() && M1.cols() == M2.cols() &&
                                  M1.nonZeros() == M2.nonZeros());
                bool bitwiseGhost = sameShape;
                if (sameShape && M1.nonZeros() > 0)
                {
                    bitwiseGhost = bitwiseGhost &&
                        0 == std::memcmp(M1.outerIndexPtr(), M2.outerIndexPtr(),
                                         (size_t)(M1.outerSize()+1)*sizeof(index_t)) &&
                        0 == std::memcmp(M1.innerIndexPtr(), M2.innerIndexPtr(),
                                         (size_t)M1.nonZeros()*sizeof(index_t)) &&
                        0 == std::memcmp(M1.valuePtr(), M2.valuePtr(),
                                         (size_t)M1.nonZeros()*sizeof(real_t));
                }
                bool hasNonzeroValue = false;
                for (index_t kk = 0; kk != M1.nonZeros() && !hasNonzeroValue; ++kk)
                    if (0.0 != M1.valuePtr()[kk]) hasNonzeroValue = true;

                const bool bPass = checks.allPass() && (nGhost > 0) && (M1.nonZeros() > 0) &&
                                  hasNonzeroValue && bitwiseGhost;
                gsInfo << "CHECK(b) case=" << caseName << " n=" << n << " ghostFaces=" << nGhost
                      << " nnz=" << M1.nonZeros() << " ghost matrix bitwise equal  "
                      << (bPass ? "PASS" : "FAIL") << "\n";
                allPass = allPass && bPass;

                const bool cPass = checkC(caseName, n, S);
                allPass = allPass && cPass;
            }
        }
    }

    // TetClipSignDomain classification: sphere/rotcube x r=0..3, independent
    // of --case/-r (see checkSignDomain()'s own doxygen).
    const bool signDomainPass = checkSignDomain(cfg.n0, cfg.p);
    allPass = allPass && signDomainPass;

    // gsImmersedPullback.h self-test: pull-back, Nanson, minSingularValueJ,
    // over both shipped meshes and every background map. --geo does not
    // otherwise change --study check.
    const bool pullbackPass = gsTetClip::pullbackSelfTest();
    allPass = allPass && pullbackPass;

    gsInfo << (allPass ? "ALL CHECKS PASS\n" : "SOME CHECKS FAILED\n");
    return allPass;
}

//----------------------------------------------------------------------------
// --study volume: background-grid and mesh-loading helpers.
//----------------------------------------------------------------------------

/// The fixed box [-1,1]^3's grid at refinement level \a r: n = \a n0 << \a r
/// cells per direction. No mesh check (see makeGrid()): used where the grid
/// is needed before the mesh exists in the grid's own, parametric,
/// coordinates, i.e. before the pull-back.
gsTetClip::Grid3 makeBoxGrid(index_t n0, index_t r)
{
    const index_t n = n0 << r;
    gsTetClip::Grid3 grid;
    grid.x0 = grid.y0 = grid.z0 = -1.0;
    grid.h  = 2.0/(real_t)n;
    grid.n  = n;
    return grid;
}

/// Background grid for case \a caseName's mesh \a M at refinement level
/// \a r: makeBoxGrid(n0,r), then GISMO_ENSUREs \a M's bounding box lies
/// strictly inside it. \a who names the caller in the error message. \a M
/// may be in physical or parametric coordinates -- the ENSURE only compares
/// it against the grid it is itself expressed in.
gsTetClip::Grid3 makeGrid(const gsTetClip::TetMesh & M, index_t n0, index_t r, const std::string & who)
{
    const gsTetClip::Grid3 grid = makeBoxGrid(n0, r);
    GISMO_ENSURE(M.lo[0] > grid.x0 && M.hi[0] < grid.x0 + grid.n*grid.h &&
                M.lo[1] > grid.y0 && M.hi[1] < grid.y0 + grid.n*grid.h &&
                M.lo[2] > grid.z0 && M.hi[2] < grid.z0 + grid.n*grid.h,
                who << ": mesh bbox does not lie strictly inside the background box at n=" << (n0 << r));
    return grid;
}

/// Loads case \a caseName's PHYSICAL mesh: \a cfg.meshOverride if
/// non-empty, else the case's own file (meshFile()), resolved by
/// gsFileManager::find and read by gsTetClip::readMsh41.
gsTetClip::PhysTetMesh loadPhysMesh(const Config & cfg, const std::string & caseName)
{
    const std::string file = cfg.meshOverride.empty() ? meshFile(caseName) : cfg.meshOverride;
    const std::string resolvedFile = gsFileManager::find(file);
    GISMO_ENSURE(!resolvedFile.empty(), "loadPhysMesh: mesh file for --case " << caseName
                << " not found.");
    return gsTetClip::PhysTetMesh{ gsTetClip::readMsh41(resolvedFile) };
}

//----------------------------------------------------------------------------
// --study volume: test-case constants.
//----------------------------------------------------------------------------

const real_t SPHERE_CENTER_X =  0.03;
const real_t SPHERE_CENTER_Y = -0.02;
const real_t SPHERE_CENTER_Z =  0.01;
const real_t SPHERE_RADIUS   =  0.55;
const real_t ROTCUBE_SIDE    =  0.8;

/// `fluxmom_global` under `--geo bubble` measures TRUNCATION error (see
/// \ref gatePhysical's own doxygen: neither the 3p+1 Full-cell reference
/// rule nor the clip volume/boundary rules are exact once G is curved),
/// with two parts visible in the sweep below. (a) Cut-cell/boundary
/// truncation of the clip rules at tet-piece scale: the only part at
/// r=0, where no cell is Full; at p=1 it is ~1e-15/1e-14 on the L=0
/// meshes and falls to ~3e-17/7e-17 on the finer L=1 meshes. (b)
/// Full-cell truncation of the 3p+1 tensor Gauss rule, present from r=1
/// on: at p=1 it dominates (~1e-12/5e-13 at r=1, ~6e-15 to 9e-15 at r=2)
/// and is nearly mesh-independent (L=0 vs L=1 agree to 2-3 significant
/// digits at r=1), since the Full cells are nearly the same set for both
/// meshes. By p=2 both parts are below the ~1e-17 rounding floor (every
/// L=0 sphere row at p=2/p=3 except p=2/r=2 prints 3.134969e-17). The
/// decision rule below therefore fires on the p=1 contrast between r=0
/// (no Full cells) and r>=1, on floor noise at p>=2, and on the p=1-vs-
/// p>=2 gap -- not on an error that grows under refinement. Sweep command:
///
///   cd /home/hverhelst/Code/gismo_immersed
///   SW=/tmp/sweep; mkdir -p $SW
///   for P in 1 2 3; do for C in sphere rotcube; do
///     OMP_NUM_THREADS=2 stdbuf -oL -eL timeout 900 \
///       build/bin/immersed_tetmesh_poisson_example --study volume \
///       --geo bubble --mode clip --case $C -k $P -r 2 2>&1 \
///       | tee $SW/sweep_${C}_L0_p${P}.txt
///   done; done
///   # L=1: same loop with --mesh <build/tetmesh_refined/*_L1.msh>, the
///   # once-gmsh-refined sphere/rotcube meshes (filedata/volumes/refine_tetmesh.sh)
///
/// case=sphere  geo=bubble  eps=0.3  fluxmom_global relerr (case x L x p x r)
///   L=0 p=1: r=0 1.936731e-15  r=1 1.182447e-12  r=2 9.434489e-15
///   L=0 p=2: r=0 3.134969e-17  r=1 3.134969e-17  r=2 1.402001e-17
///   L=0 p=3: r=0 3.134969e-17  r=1 3.134969e-17  r=2 3.134969e-17
///   L=1 p=1: r=0 2.922918e-17  r=1 1.183332e-12  r=2 9.052505e-15
///   L=1 p=2: r=0 2.394455e-17  r=1 2.394455e-17  r=2 3.135077e-17
///   L=1 p=3: not run (memory) -- the sweep process was killed by an
///            idle-memory-pressure reaper before any row printed, and was
///            not restarted.
/// case=rotcube
///   L=0 p=1: r=0 1.140206e-14  r=1 5.175777e-13  r=2 5.766169e-15
///   L=0 p=2: r=0 3.090603e-17  r=1 3.090603e-17  r=2 6.910799e-18
///   L=0 p=3: r=0 2.360487e-17  r=1 3.090603e-17  r=2 3.090603e-17
///   L=1 p=1: r=0 6.528357e-17  r=1 5.164703e-13  r=2 6.070591e-15
///   L=1 p=2: r=0 1.227810e-17  r=1 1.055700e-17  r=2 1.055700e-17
///   L=1 p=3: not run (memory), same reason as sphere L=1 p=3 above.
///
/// Decision rule (applied literally to the 10 series that ran): GROWING
/// -- relerr(r=2) > relerr(r=0) in the same (case,L,p) series -- holds for
/// sphere/L0/p1, sphere/L1/p1, sphere/L1/p2 and rotcube/L0/p3,
/// rotcube/L1/p1. SPREAD = max/min over every entry that ran =
/// 1.183332e-12 / 6.910799e-18 ~ 1.71e5 > 100. Both conditions alone
/// already trigger the same branch of the rule for a SINGLE tolerance
/// spanning every p: REQUIRED = false, TOL = 10 x max = 1.183332e-11
/// (`BUBBLE_FLUXMOM_TOL`), giving REPORT at every p.
///
/// The row is instead tiered by p (`BUBBLE_FLUXMOM_REQUIRED_MIN_P`,
/// `BUBBLE_FLUXMOM_TOL_REQUIRED`), since the single-tolerance rule above
/// conflates two different regimes. At p=1 the Full-cell 3p+1 rule
/// carries genuine, not rounding-level, truncation (part (b) above):
/// ~1.18e-12/5.18e-13 at r=1, ~9e-15/6e-15 at r=2 -- any tolerance placed
/// there is a measured bound on that truncation, not a derived one, so
/// p=1 stays Report at `BUBBLE_FLUXMOM_TOL`. At p>=2 every one of the 10
/// series' p=2/p=3 entries sits at the rounding floor -- the largest is
/// 3.135077e-17 (sphere/L1/p2/r=2), the rest 6.9e-18 to 3.1e-17 -- about
/// 3e3 (three orders of magnitude) below `BUBBLE_FLUXMOM_TOL_REQUIRED` =
/// 1e-13, which in turn sits about ten orders of magnitude below this
/// same `fluxmom_global` row under a wrong served normal
/// (`--study volume --case all --mode all -r 1`, `--geo identity`
/// default: `--mode momrule`'s `fluxmom_global` FAILs at relerr
/// 2.017442e-03 to 2.950480e-02 across sphere/rotcube, r=0/1; momrule
/// itself is rejected outright for curved G, see the header): tight
/// enough to catch a quadrature regression, loose enough never to fire
/// on rounding. So the row is Required at `BUBBLE_FLUXMOM_TOL_REQUIRED`
/// for p >= `BUBBLE_FLUXMOM_REQUIRED_MIN_P` (= 2), and stays Report at
/// `BUBBLE_FLUXMOM_TOL` for p < `BUBBLE_FLUXMOM_REQUIRED_MIN_P`.
const real_t  BUBBLE_FLUXMOM_TOL               = 1.183332e-11;
const index_t BUBBLE_FLUXMOM_REQUIRED_MIN_P    = 2;
const real_t  BUBBLE_FLUXMOM_TOL_REQUIRED      = 1e-13;

real_t sphereVolumeExact() { return (4.0/3.0)*EIGEN_PI*SPHERE_RADIUS*SPHERE_RADIUS*SPHERE_RADIUS; }
real_t sphereAreaExact()   { return 4.0*EIGEN_PI*SPHERE_RADIUS*SPHERE_RADIUS; }
real_t rotcubeVolumeExact() { return ROTCUBE_SIDE*ROTCUBE_SIDE*ROTCUBE_SIDE; }
real_t rotcubeAreaExact()   { return 6.0*ROTCUBE_SIDE*ROTCUBE_SIDE; }

/// First line of a (possibly multi-line) exception message.
std::string firstLine(const std::string & s)
{
    const size_t nl = s.find('\n');
    return (std::string::npos == nl) ? s : s.substr(0, nl);
}

//----------------------------------------------------------------------------
// Cell-rule construction (shared with --study poisson).
//----------------------------------------------------------------------------

/// gsQuadRule returning ONE cell's precomputed points unchanged, whatever box
/// it is handed; the adapter that lets gsMomentRule (which owns and calls a
/// gsQuadRule) compress a rule that ClipStreamer already produced, without
/// re-clipping. Holds non-owning pointers: valid only while the pointed-to
/// matrices live (one loop iteration).
class FixedPointsRule : public gsQuadRule<real_t>
{
public:
    FixedPointsRule(const gsMatrix<real_t> * nodes, const gsVector<real_t> * weights)
    : m_nodes(nodes), m_weights(weights) {}

    using gsQuadRule<real_t>::mapTo;
    void mapTo(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
              gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    {
        GISMO_UNUSED(lower); GISMO_UNUSED(upper);
        nodes = *m_nodes;
        weights = *m_weights;
    }

private:
    const gsMatrix<real_t> * m_nodes;
    const gsVector<real_t> * m_weights;
};

/// Per-cell outcome of one cell-rule pass. "In" = the uncompressed clip rule
/// of the cell, "Out" = the rule the mode SERVES for that cell (== In for
/// clip mode).
struct CellLog
{
    index_t nVolIn = 0, nVolOut = 0, nBdrIn = 0, nBdrOut = 0;
    real_t  volSum = 0, bdrSum = 0, fluxSum = 0;          // Kahan within the cell, of the served rule
    real_t  minWVol = std::numeric_limits<real_t>::infinity(),
            minWBdr = std::numeric_limits<real_t>::infinity(); // over the served weights (+inf if none)
    index_t negWVol = 0, negWBdr = 0;
    real_t  momErrVol = 0, momErrBdr = 0;                 // per-cell Q_2p rel. moment error (compressed modes)
    real_t  fluxErrBdr = 0;                               // per-cell E_K, gateTetMode's fluxmom_cell row
    bool    hasFluxErrBdr = false;                        // fluxErrBdr meaningful (clip boundary area > 0)
    real_t  resVol = 0, resBdr = 0;                       // tchakaloff: max over levelResidual
    bool    okVol = true, okBdr = true;                   // tchakaloff: result.ok
    index_t rankVol = -1, rankBdr = -1, levelsVol = 0, levelsBdr = 0;
    bool    passThroughVol = false, passThroughBdr = false; // momrule: stats().nPassThroughElements > 0
    real_t  compressSec = 0;                              // time spent in the compressor for this cell
    std::string error;                                    // non-empty iff an exception was caught
};

/// Outcome of one `buildCellRules` pass. `vol`/`bdr`, when non-null, hold
/// the SERVED rule per cell id: `vol->cell[id]` is a `CellRule3` (nodes
/// 3xm, weights m); `bdr->cell[id]` is a `CellBdrRule3` (nodes 3xm,
/// weights m, unit outward normals 3xm). All PHYSICAL coordinates. Both
/// are null in clip mode, where the served rule is the streamer's own
/// uncompressed clip rule (nothing extra is stored).
struct CellRulePass
{
    std::string mode;                                      // "clip" | "tchakaloff" | "momrule"
    memory::shared_ptr<gsTetClip::VolCellTable> vol;       // null in clip mode; else cell.size()==n^3
    memory::shared_ptr<gsTetClip::BdrCellTable> bdr;       // null in clip mode; else cell.size()==n^3
    std::vector<size_t> work;                              // processed ids, ascending
    std::vector<CellLog> log;                              // size n^3 (non-work entries default)
    real_t seconds = 0;                                    // wall-clock of the pass
};

/// Q_2p moments of a discrete measure in the tensor Legendre basis
/// orthonormal on [\a lower,\a upper] (gsTetClip::legendreOrthonormal,
/// deg = 2p, K = (2p+1)^3). \a normals == nullptr: K entries, index
/// k = kx + (deg+1)*(ky + (deg+1)*kz). \a normals != nullptr: 4K entries,
/// block b*K + k weighted by 1, n_x, n_y, n_z (b = 0..3), the stacked-
/// Tchakaloff layout. One KahanSum per entry. O(N K) time, O(K) memory:
/// NEVER forms the N x K Vandermonde (400k x 125 doubles = 400 MB/thread).
void cellMoments(const gsMatrix<real_t> & nodes, const gsVector<real_t> & weights,
                 const gsMatrix<real_t> * normals, const gsVector<real_t> & lower,
                 const gsVector<real_t> & upper, index_t p, gsVector<real_t> & m)
{
    const index_t deg = 2*p, n1 = deg+1, K = n1*n1*n1;
    const index_t blocks = (nullptr == normals) ? 1 : 4;
    std::vector<gsTetClip::KahanSum> acc((size_t)(blocks*K));

    gsVector<real_t> vx(n1), vy(n1), vz(n1);
    for (index_t c = 0; c != nodes.cols(); ++c)
    {
        gsTetClip::legendreOrthonormal(nodes(0,c), lower(0), upper(0), deg, vx);
        gsTetClip::legendreOrthonormal(nodes(1,c), lower(1), upper(1), deg, vy);
        gsTetClip::legendreOrthonormal(nodes(2,c), lower(2), upper(2), deg, vz);

        for (index_t kz = 0; kz != n1; ++kz)
        for (index_t ky = 0; ky != n1; ++ky)
        for (index_t kx = 0; kx != n1; ++kx)
        {
            const index_t k = kx + n1*(ky + n1*kz);
            const real_t base = weights[c]*vx[kx]*vy[ky]*vz[kz];
            acc[(size_t)k].add(base);
            if (nullptr != normals)
            {
                acc[(size_t)(K+k)].add(base*(*normals)(0,c));
                acc[(size_t)(2*K+k)].add(base*(*normals)(1,c));
                acc[(size_t)(3*K+k)].add(base*(*normals)(2,c));
            }
        }
    }

    m.resize(blocks*K);
    for (index_t k = 0; k != blocks*K; ++k) m[k] = acc[(size_t)k].value();
}

/// max_k |mOut_k - mIn_k| / max_k |mIn_k| (0 if both are all-zero; +inf if
/// only \a mIn is all-zero).
real_t momentRelErr(const gsVector<real_t> & mIn, const gsVector<real_t> & mOut)
{
    GISMO_ENSURE(mIn.size() == mOut.size(), "momentRelErr: size mismatch.");
    real_t maxIn = 0, maxDiff = 0;
    for (index_t k = 0; k != mIn.size(); ++k)
    {
        maxIn   = math::max(maxIn,   math::abs(mIn[k]));
        maxDiff = math::max(maxDiff, math::abs(mIn[k]-mOut[k]));
    }
    if (0 == maxIn)
        return (0 == maxDiff) ? (real_t)0 : std::numeric_limits<real_t>::infinity();
    return maxDiff/maxIn;
}

/// Value and d/dx of gsTetClip::legendreOrthonormal's basis at \a x: \a vals
/// as legendreOrthonormal itself; \a derivs(k) = d/dx [sqrt((2k+1)/(b-a)) P_k(t)],
/// t = (2x-a-b)/(b-a), via the companion three-term recurrence for the
/// UNNORMALIZED Legendre derivative, P'_{n+1}(t) = P'_{n-1}(t) + (2n+1) P_n(t)
/// (P'_0 = 0; standard identity (2n+1) P_n = P'_{n+1} - P'_{n-1}), chain-ruled
/// by dt/dx = 2/(b-a) and scaled by the same normalization constant as \a vals.
void legendreOrthonormalDeriv(real_t x, real_t a, real_t b, index_t deg,
                              gsVector<real_t> & vals, gsVector<real_t> & derivs)
{
    GISMO_ASSERT(deg >= 0, "legendreOrthonormalDeriv: deg must be >= 0");
    GISMO_ASSERT(b > a, "legendreOrthonormalDeriv: invalid interval [a,b]");

    const real_t t = (2*x - a - b)/(b - a);
    const real_t dtdx = 2.0/(b - a);
    vals.resize(deg+1); derivs.resize(deg+1);

    real_t pnm1 = 0.0, pn = 1.0;      // P_{-1}, P_0
    real_t dpnm1 = 0.0, dpn = 0.0;    // P'_{-1}, P'_0
    vals(0) = pn; derivs(0) = dpn;
    for (index_t j = 0; j != deg; ++j)
    {
        const real_t pnm2 = pnm1; pnm1 = pn;
        pn = (static_cast<real_t>(2*j+1)*t*pnm1 - static_cast<real_t>(j)*pnm2)
           / static_cast<real_t>(j+1);
        const real_t dpnm2 = dpnm1; dpnm1 = dpn;
        dpn = dpnm2 + static_cast<real_t>(2*j+1)*pnm1;
        vals(j+1) = pn; derivs(j+1) = dpn;
    }

    for (index_t k = 0; k <= deg; ++k)
    {
        const real_t c = std::sqrt(static_cast<real_t>(2*k+1)/(b-a));
        vals(k) *= c; derivs(k) *= c*dtdx;
    }
}

/// d(Q_k)/dx_i moments of a volume quadrature rule (\a nodes, \a weights)
/// against the SAME tensor Legendre basis \ref cellMoments uses (deg = 2p,
/// K = (2p+1)^3, index k = kx + (deg+1)*(ky + (deg+1)*kz)), block layout
/// i*K+k for i = 0,1,2 (dQ_k/dx, dQ_k/dy, dQ_k/dz). Feeds gateTetMode's
/// fluxmom_global row's divergence-theorem right-hand side. One KahanSum
/// per entry, same O(N K) discipline as \ref cellMoments.
void volDerivMoments(const gsMatrix<real_t> & nodes, const gsVector<real_t> & weights,
                     const gsVector<real_t> & lower, const gsVector<real_t> & upper,
                     index_t p, gsVector<real_t> & m)
{
    const index_t deg = 2*p, n1 = deg+1, K = n1*n1*n1;
    std::vector<gsTetClip::KahanSum> acc((size_t)(3*K));

    gsVector<real_t> vx(n1), vy(n1), vz(n1), dvx(n1), dvy(n1), dvz(n1);
    for (index_t c = 0; c != nodes.cols(); ++c)
    {
        legendreOrthonormalDeriv(nodes(0,c), lower(0), upper(0), deg, vx, dvx);
        legendreOrthonormalDeriv(nodes(1,c), lower(1), upper(1), deg, vy, dvy);
        legendreOrthonormalDeriv(nodes(2,c), lower(2), upper(2), deg, vz, dvz);

        for (index_t kz = 0; kz != n1; ++kz)
        for (index_t ky = 0; ky != n1; ++ky)
        for (index_t kx = 0; kx != n1; ++kx)
        {
            const index_t k = kx + n1*(ky + n1*kz);
            const real_t w = weights[c];
            acc[(size_t)(0*K+k)].add(w*dvx[kx]*vy[ky]*vz[kz]);
            acc[(size_t)(1*K+k)].add(w*vx[kx]*dvy[ky]*vz[kz]);
            acc[(size_t)(2*K+k)].add(w*vx[kx]*vy[ky]*dvz[kz]);
        }
    }

    m.resize(3*K);
    for (index_t k = 0; k != 3*K; ++k) m[k] = acc[(size_t)k].value();
}

/// factor(k) = c(kx,dx)*c(ky,dy)*c(kz,dz), c(j,d) = sqrt((2j+1)/d): the
/// scale gsTetClip::legendreOrthonormal applies to the RAW Legendre product
/// P_kx(u)P_ky(v)P_kz(w) (sup norm 1 on the mapped box, since |P_j|<=1 on
/// [-1,1]) to produce the ORTHONORMAL basis \ref cellMoments and
/// \ref volDerivMoments actually evaluate. Dividing one of their moment
/// entries by factor(k) recovers the raw-Legendre moment: ||q||_inf = 1 on
/// the mapped box, so |int_Gamma q n_i ds| <= area there, independent of
/// box size. The orthonormal scale itself grows like width^(-3/2) (three
/// factors of sqrt((2k+1)/width)), so a fixed absolute (rounding-level)
/// mismatch in an orthonormal moment inflates without bound as the box
/// shrinks under refinement -- exactly the failure mode this factor
/// removes. The moments_vol/moments_bdr rows (momentRelErr) do not use this factor: they compare orthonormal
/// moments relative to the cell's own largest moment, which is already scale-free.
void legendreRawFactor(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
                       index_t deg, gsVector<real_t> & factor)
{
    const index_t n1 = deg+1, K = n1*n1*n1;
    gsVector<real_t> cx(n1), cy(n1), cz(n1);
    for (index_t j = 0; j <= deg; ++j)
    {
        cx[j] = math::sqrt(static_cast<real_t>(2*j+1)/(upper(0)-lower(0)));
        cy[j] = math::sqrt(static_cast<real_t>(2*j+1)/(upper(1)-lower(1)));
        cz[j] = math::sqrt(static_cast<real_t>(2*j+1)/(upper(2)-lower(2)));
    }
    factor.resize(K);
    for (index_t kz = 0; kz != n1; ++kz)
    for (index_t ky = 0; ky != n1; ++ky)
    for (index_t kx = 0; kx != n1; ++kx)
        factor[kx + n1*(ky + n1*kz)] = cx[kx]*cy[ky]*cz[kz];
}

/// One omp-parallel pass over the work cells (Cut cells, or cells with a
/// non-empty triangle bucket). Mode clip: nothing stored. tchakaloff/
/// momrule: the compressed rules are stored in the tables (Cut-cell volume,
/// every work cell's boundary). \a phiH is required for momrule (boundary
/// normals), ignored otherwise. \a momentCheck computes CellLog::momErr*
/// and CellLog::fluxErrBdr/hasFluxErrBdr (compressed modes only; the latter
/// is gateTetMode's per-cell fluxmom_cell input, a with-normals moment
/// comparison independent of momErrBdr's own mode-dependent normal
/// weighting). Every exception inside the parallel region is
/// caught into the cell's own CellLog::error: one escaping the `#pragma omp
/// parallel` region would terminate the process (GISMO_ERROR/ENSURE throw
/// std::runtime_error).
///
/// Thread-safe: each loop iteration writes only `log[id]`, `vol->cell[id]`
/// and `bdr->cell[id]`, which are distinct, pre-sized elements -- no two
/// iterations touch the same memory. Complexity per work cell: one clip of
/// the cell (the streamer's own cost) plus the compressor (tchakaloff/
/// momrule only), plus O(N_cell*K) for \a momentCheck's moment comparison
/// (K = (2p+1)^3, N_cell the cell's input node count).
CellRulePass buildCellRules(const gsTetClip::ClipStreamer & S, const std::string & mode, index_t p,
                            const gsMeshSignedDist<real_t> * phiH, bool momentCheck)
{
    GISMO_ENSURE("clip" == mode || "tchakaloff" == mode || "momrule" == mode,
                "buildCellRules: mode must be clip|tchakaloff|momrule, got '" << mode << "'.");
    GISMO_ENSURE("momrule" != mode || nullptr != phiH,
                "buildCellRules: momrule mode requires a mesh level set.");

    CellRulePass P;
    P.mode = mode;
    const bool compressed = ("clip" != mode);

    memory::shared_ptr<const gsTetClip::CellIndex> idx = S.index();
    const index_t n = idx->grid.n;
    const size_t N3 = (size_t)n*(size_t)n*(size_t)n;

    P.log.assign(N3, CellLog());
    if (compressed)
    {
        P.vol = memory::make_shared(new gsTetClip::VolCellTable());
        P.bdr = memory::make_shared(new gsTetClip::BdrCellTable());
        P.vol->cell.resize(N3);
        P.bdr->cell.resize(N3);
    }

    for (size_t id = 0; id != N3; ++id)
        if (gsTetClip::Cut == idx->status[id] || !S.triBuckets()[id].empty())
            P.work.push_back(id);

    gsStopwatch sw;
    const index_t W = (index_t)P.work.size();
    #pragma omp parallel for schedule(dynamic,1)
    for (index_t w = 0; w < W; ++w)
    {
        const size_t id = P.work[(size_t)w];
        CellLog & L = P.log[id];
        try
        {
            index_t i,j,k; idx->ijk(id, i,j,k);
            gsVector<real_t> lower(3), upper(3);
            lower << idx->X[i], idx->Y[j], idx->Z[k];
            upper << idx->X[i+1], idx->Y[j+1], idx->Z[k+1];

            // --- Volume, Cut cells only. ---
            if (gsTetClip::Cut == idx->status[id])
            {
                gsMatrix<real_t> nd; gsVector<real_t> wt;
                S.volRule(id, nd, wt);
                L.nVolIn = nd.cols();

                gsMatrix<real_t> sn; gsVector<real_t> sw_;
                if ("clip" == mode)
                {
                    sn = nd; sw_ = wt;
                }
                else if ("tchakaloff" == mode)
                {
                    gsStopwatch csw;
                    const gsTetClip::TchakaloffResult R =
                        gsTetClip::tchakaloffCompress(nd, wt, lower, upper, p);
                    L.compressSec += csw.stop();

                    sn.resize(3, (index_t)R.indices.size());
                    for (size_t c = 0; c != R.indices.size(); ++c)
                        sn.col((index_t)c) = nd.col(R.indices[c]);
                    sw_ = R.weights;

                    L.okVol   = R.ok;
                    L.resVol  = R.levelResidual.empty() ? (real_t)0
                              : *std::max_element(R.levelResidual.begin(), R.levelResidual.end());
                    L.rankVol   = R.rank;
                    L.levelsVol = R.levels;
                }
                else // momrule
                {
                    gsStopwatch csw;
                    typename gsMomentRule<real_t>::uPtr mr = gsMomentRule<real_t>::make(
                        gsQuadRule<real_t>::uPtr(new FixedPointsRule(&nd, &wt)),
                        gsMomentRule<real_t>::exactnessOrder((short_t)p, 3));
                    mr->mapTo(lower, upper, sn, sw_);
                    L.compressSec += csw.stop();
                    L.passThroughVol = mr->stats().nPassThroughElements > 0;
                }

                if (momentCheck && "clip" != mode)
                {
                    gsVector<real_t> mIn, mOut;
                    cellMoments(nd, wt, nullptr, lower, upper, p, mIn);
                    cellMoments(sn, sw_, nullptr, lower, upper, p, mOut);
                    L.momErrVol = momentRelErr(mIn, mOut);
                }

                L.nVolOut = sn.cols();
                gsTetClip::KahanSum vsum;
                for (index_t c = 0; c != sw_.size(); ++c)
                {
                    vsum.add(sw_[c]);
                    if (sw_[c] < L.minWVol) L.minWVol = sw_[c];
                    if (sw_[c] < 0) ++L.negWVol;
                }
                L.volSum = vsum.value();

                if (compressed)
                    P.vol->cell[id] = gsTetClip::CellRule3{ give(sn), give(sw_) };
            }

            // --- Boundary, every work id. ---
            gsMatrix<real_t> bn, bnrm; gsVector<real_t> bw;
            S.bdrRule(id, bn, bw, bnrm);
            L.nBdrIn = bn.cols();

            if (bn.cols() > 0)
            {
                gsMatrix<real_t> sn, snrm; gsVector<real_t> sw_;
                if ("clip" == mode)
                {
                    sn = bn; sw_ = bw; snrm = bnrm;
                }
                else if ("tchakaloff" == mode)
                {
                    gsStopwatch csw;
                    const gsTetClip::TchakaloffResult R =
                        gsTetClip::tchakaloffCompressBoundary(bn, bw, bnrm, lower, upper, p);
                    L.compressSec += csw.stop();

                    sn.resize(3, (index_t)R.indices.size());
                    snrm.resize(3, (index_t)R.indices.size());
                    for (size_t c = 0; c != R.indices.size(); ++c)
                    {
                        sn.col((index_t)c)   = bn.col(R.indices[c]);
                        snrm.col((index_t)c) = bnrm.col(R.indices[c]);
                    }
                    sw_ = R.weights;

                    L.okBdr   = R.ok;
                    L.resBdr  = R.levelResidual.empty() ? (real_t)0
                              : *std::max_element(R.levelResidual.begin(), R.levelResidual.end());
                    L.rankBdr   = R.rank;
                    L.levelsBdr = R.levels;
                }
                else // momrule
                {
                    gsStopwatch csw;
                    typename gsMomentRule<real_t>::uPtr mr = gsMomentRule<real_t>::make(
                        gsQuadRule<real_t>::uPtr(new FixedPointsRule(&bn, &bw)),
                        gsMomentRule<real_t>::exactnessOrder((short_t)p, 3));
                    mr->mapTo(lower, upper, sn, sw_);
                    L.compressSec += csw.stop();
                    L.passThroughBdr = mr->stats().nPassThroughElements > 0;

                    // The Muller-Kummer-Oberlack convention: the unit gradient
                    // of the mesh signed distance, evaluated at the (generally
                    // off-surface) served nodes; on pass-through cells the
                    // nodes are exactly on the surface and this is the outward
                    // pseudonormal.
                    gsMatrix<real_t> g;
                    phiH->deriv_into(sn, g);
                    snrm.resize(3, sn.cols());
                    for (index_t c = 0; c != sn.cols(); ++c)
                    {
                        const real_t nrm = g.col(c).norm();
                        GISMO_ENSURE(nrm > 0, "buildCellRules: zero-norm level-set gradient at a "
                                    "served momrule boundary node.");
                        snrm.col(c) = g.col(c)/nrm;
                    }
                }

                if (momentCheck && "clip" != mode)
                {
                    // With-normals Q_2p moments of the clip reference (bn/bw/bnrm)
                    // against the mode's served rule (sn/sw_/snrm): feeds
                    // fluxErrBdr below for every mode, and IS momErrBdr's own
                    // comparison for tchakaloff (whose moments_bdr row already
                    // weights by the served normal).
                    gsVector<real_t> mInN, mOutN;
                    cellMoments(bn, bw, &bnrm, lower, upper, p, mInN);
                    cellMoments(sn, sw_, &snrm, lower, upper, p, mOutN);

                    if ("tchakaloff" == mode)
                    {
                        L.momErrBdr = momentRelErr(mInN, mOutN);
                    }
                    else // momrule: K scalar moments, no normal weighting --
                         // this is exactly the gap fluxErrBdr below closes,
                         // since momrule's served normal (the mesh level-set
                         // gradient at an off-surface node, wrong near
                         // polyhedral edges/vertices) never enters mIn/mOut.
                    {
                        gsVector<real_t> mIn, mOut;
                        cellMoments(bn, bw, nullptr, lower, upper, p, mIn);
                        cellMoments(sn, sw_, nullptr, lower, upper, p, mOut);
                        L.momErrBdr = momentRelErr(mIn, mOut);
                    }

                    // fluxmom_cell (gateTetMode): E_K = max_{q,i} of the
                    // with-normals block's mismatch, RAW-Legendre-normalized
                    // (legendreRawFactor -- ||q||_inf = 1 on this cell's own
                    // box; mInN/mOutN themselves stay orthonormal, the basis
                    // momErrBdr above compares), divided by the clip
                    // reference's own boundary area A_K = sum(bw) (not
                    // sw_'s: A_K must be independent of the mode under test,
                    // exactly as gateTetMode's fluxmom_global uses the clip
                    // volume rule for its own reference side).
                    gsTetClip::KahanSum areaClip;
                    for (index_t c = 0; c != bw.size(); ++c) areaClip.add(bw[c]);
                    const real_t A_K = areaClip.value();
                    if (0.0 != A_K)
                    {
                        const index_t deg = 2*p, n1 = deg+1, K = n1*n1*n1;
                        gsVector<real_t> factor;
                        legendreRawFactor(lower, upper, deg, factor);
                        real_t eK = 0;
                        for (index_t ii = 0; ii != 3; ++ii)
                            for (index_t kk = 0; kk != K; ++kk)
                                eK = math::max(eK, math::abs(mOutN[(ii+1)*K+kk]-mInN[(ii+1)*K+kk])
                                               /(factor[kk]*A_K));
                        L.fluxErrBdr = eK;
                        L.hasFluxErrBdr = true;
                    }
                }

                L.nBdrOut = sn.cols();
                gsTetClip::KahanSum bsum, fsum;
                for (index_t c = 0; c != sw_.size(); ++c)
                {
                    bsum.add(sw_[c]);
                    fsum.add(sw_[c]*sn.col(c).dot(snrm.col(c)));
                    if (sw_[c] < L.minWBdr) L.minWBdr = sw_[c];
                    if (sw_[c] < 0) ++L.negWBdr;
                }
                L.bdrSum  = bsum.value();
                L.fluxSum = fsum.value();

                if (compressed)
                    P.bdr->cell[id] = gsTetClip::CellBdrRule3{ give(sn), give(sw_), give(snrm) };
            }
        }
        catch (const std::exception & e)
        {
            L.error = e.what();
        }
    }
    P.seconds = sw.stop();

    return P;
}

//----------------------------------------------------------------------------
// Gate evaluation.
//----------------------------------------------------------------------------

/// One printed `GATE` row: a scalar check (`value` vs `ref`, `relerr` vs
/// `tol`) or a per-cell check whose `failCells` lists the ascending ids
/// that fail it (only the first 20 are ever printed). `positive` marks a
/// row whose `tol` prints as the literal `positive` (PASS iff every served
/// weight is strictly greater than 0).
struct GateRow
{
    std::string check;
    real_t value = 0, ref = 0, relerr = 0, tol = 0;
    enum Kind { Required, Report } kind = Report;
    bool positive = false;            // tol printed as "positive": PASS iff every weight > 0
    bool pass = true;                 // meaningful for Required rows only
    std::vector<size_t> failCells;    // ascending cell ids (printed: first 20)
};
/// All `GateRow`s of one (case, mode, r) gate, plus the AND of every
/// `Required` row's `pass`, computed by \ref gateTetMode / \ref gateAlgoim /
/// \ref gatePhysical exactly as for every mode (see \ref printGateRow / \ref runVolumeStudy
/// for tchakaloff's separate DEFERRED printing/counting policy).
struct GateResult
{
    std::string caseName, mode; index_t r = 0, n = 0;
    std::string geo;                  // non-identity --geo only; empty on the identity path
    std::vector<GateRow> rows;
    bool requiredPass = true;         // AND over the Required rows
};

/// Volume, area and flux totals of one gate pass.
struct Totals { real_t V = 0, A = 0, F = 0; };

/// Two-level-Kahan V_q/A_q/F_q of \a P over the streamer \a S's grid: a
/// per-cell KahanSum, then a KahanSum over cells in INCREASING id, so the
/// totals are independent of the OMP thread count and schedule. clip:
/// Cut cells reuse log[id].volSum/bdrSum/fluxSum, Full cells get a fresh
/// tensor-Gauss rule of order (p+1,p+1,p+1) (the SAME rule VolLookupRule
/// serves them), Empty cells contribute 0. Compressed modes: a serial pass
/// through VolLookupRule/the boundary table -- the SAME lookup path the
/// Poisson study assembles with.
Totals tetModeTotals(const CellRulePass & P, const gsTetClip::ClipStreamer & S, index_t p)
{
    const gsTetClip::CellIndex & idx = *S.index();
    const index_t n = idx.grid.n;
    const size_t N3 = (size_t)n*(size_t)n*(size_t)n;
    const gsVector<index_t> nG = gsVector<index_t>::Constant(3, p+1);

    gsTetClip::KahanSum totV, totA, totF;

    memory::shared_ptr<gsTetClip::VolLookupRule> volLookup;
    if ("clip" != P.mode)
        volLookup = memory::make_shared(new gsTetClip::VolLookupRule(S.index(), P.vol, nG));

    for (size_t id = 0; id != N3; ++id)
    {
        index_t i,j,k; idx.ijk(id, i,j,k);
        gsVector<real_t> lower(3), upper(3);
        lower << idx.X[i], idx.Y[j], idx.Z[k];
        upper << idx.X[i+1], idx.Y[j+1], idx.Z[k+1];

        if ("clip" == P.mode)
        {
            if (gsTetClip::Full == idx.status[id])
            {
                gsMatrix<real_t> nd; gsVector<real_t> w;
                gsGaussRule<real_t>(nG).mapTo(lower, upper, nd, w);
                gsTetClip::KahanSum cellV;
                for (index_t c = 0; c != w.size(); ++c) cellV.add(w[c]);
                totV.add(cellV.value());
            }
            else if (gsTetClip::Cut == idx.status[id])
                totV.add(P.log[id].volSum);

            totA.add(P.log[id].bdrSum);
            totF.add(P.log[id].fluxSum);
        }
        else
        {
            gsMatrix<real_t> nd; gsVector<real_t> w;
            volLookup->mapTo(lower, upper, nd, w);
            gsTetClip::KahanSum cellV;
            for (index_t c = 0; c != w.size(); ++c) cellV.add(w[c]);
            totV.add(cellV.value());

            gsMatrix<real_t> bn, bnrm; gsVector<real_t> bw;
            P.bdr->bdrRule(id, bn, bw, bnrm);
            gsTetClip::KahanSum cellA, cellF;
            for (index_t c = 0; c != bw.size(); ++c)
            {
                cellA.add(bw[c]);
                cellF.add(bw[c]*bn.col(c).dot(bnrm.col(c)));
            }
            totA.add(cellA.value());
            totF.add(cellF.value());
        }
    }

    Totals t; t.V = totV.value(); t.A = totA.value(); t.F = totF.value();
    return t;
}

/// Gate of a tet mode (clip/tchakaloff/momrule) on the pass \a P built from
/// \a S. Rows: `volume`, `area`, `flux`, `volume_exact`, `area_exact`
/// (Report always), `fluxmom_global` (every mode), `moments_vol`/
/// `moments_bdr`/`fluxmom_cell` (compressed modes only), `residual_vol`/
/// `residual_bdr` (tchakaloff only), `minweight_vol`/`minweight_bdr`
/// (Required/positive for tchakaloff, Report otherwise), `run`. Tolerances:
/// 1e-13 (clip) / 1e-12 (compressed) for volume/area, same for flux except
/// it is Report for momrule; moments/fluxmom 1e-12; residuals 1e-13; `run`
/// requires 0 errored cells.
///
/// Flux identity: for the mesh's closed polyhedral boundary, the
/// divergence theorem gives `\oint x.n dS = \int div(x) dV = 3 V_mesh`,
/// which is why the `flux` row's reference is `3*V_mesh` rather than a
/// quantity computed from the quadrature itself.
///
/// `fluxmom_global`/`fluxmom_cell`: the SCALAR moments `moments_vol`/
/// `moments_bdr` check (Q_2p Legendre moments, no normal) are exactly what
/// gsMomentRule reproduces by construction, so momrule cannot fail them;
/// what a wrong served normal breaks is the Nitsche flux term
/// `int_Gamma (grad u . n) v ds`, which needs the VECTOR moments
/// `int_Gamma q n_i ds`. `fluxmom_global` is the Q_2p generalization of the
/// `flux` row's divergence-theorem identity `int_Gamma q n_i dS =
/// int_Omega d(q)/dx_i dV`: the left side (B) is built from the mode's own
/// SERVED boundary rule/normals (what the Nitsche assembly actually uses),
/// the right side (V) from the clip volume rule -- a REFERENCE rule here,
/// independent of the mode under test. With the identity background map
/// used here it integrates polynomials on the polyhedron to degree 6p;
/// under a non-identity map it is only a reference, not an exact rule.
/// `fluxmom_cell` is the same identity's per-Cut-cell counterpart, each
/// mode's served rule against the clip reference on the SAME cell (Q
/// mapped to that cell's own box, `CellLog::fluxErrBdr`, computed in
/// \ref buildCellRules) -- it isolates which cells carry the bad normal
/// instead of only the domain-wide sum.
///
/// \c GateRow::kind and \c GateRow::pass are computed identically for
/// every mode, including tchakaloff: the DEFERRED downgrade of its
/// Required rows (they print but do not count toward `GATE SUMMARY`) is a
/// printing/counting policy applied by the caller (\ref printGateRow,
/// \ref runVolumeStudy), not a change to this function's verdicts. Both fluxmom
/// rows count as Required for every mode including momrule; only
/// \ref runPoissonStudy's solve/skip decision treats a momrule-only
/// fluxmom failure specially (see its own comment).
GateResult gateTetMode(const CellRulePass & P, const gsTetClip::ClipStreamer & S, index_t p,
                       real_t V_mesh, real_t A_mesh, real_t V_exact, real_t A_exact,
                       const std::string & caseName, index_t r)
{
    const gsTetClip::CellIndex & idx = *S.index();
    const size_t N3 = P.log.size();
    const std::string & mode = P.mode;

    GateResult G;
    G.caseName = caseName; G.mode = mode; G.r = r; G.n = idx.grid.n;

    const Totals tot = tetModeTotals(P, S, p);

    const real_t tolTight = 1e-13, tolLoose = 1e-12;
    const real_t volTol  = ("clip" == mode) ? tolTight : tolLoose;
    const real_t areaTol = volTol;
    const real_t fluxTol = ("clip" == mode) ? tolTight : tolLoose;
    const bool   fluxRequired = ("momrule" != mode);

    auto scalarRow = [&](const std::string & name, real_t value, real_t ref, real_t tol,
                         GateRow::Kind kind)
    {
        GateRow row;
        row.check = name; row.value = value; row.ref = ref; row.tol = tol; row.kind = kind;
        row.relerr = math::abs(value-ref)/math::abs(ref);
        row.pass = (row.relerr <= tol);
        G.rows.push_back(row);
        if (GateRow::Required == kind) G.requiredPass = G.requiredPass && row.pass;
    };

    scalarRow("volume", tot.V, V_mesh, volTol, GateRow::Required);
    scalarRow("area",   tot.A, A_mesh, areaTol, GateRow::Required);
    scalarRow("flux",   tot.F, 3.0*V_mesh, fluxTol,
             fluxRequired ? GateRow::Required : GateRow::Report);
    scalarRow("volume_exact", tot.V, V_exact, 0, GateRow::Report);
    scalarRow("area_exact",   tot.A, A_exact, 0, GateRow::Report);

    // fluxmom_global: max_{q,i} |B_{q,i}-V_{q,i}| / A_mesh, q the tensor
    // Legendre Q_2p basis mapped to the fixed background box [-1,1]^3 (so
    // ||q||_inf = 1). B is the mode's own served boundary rule/normals
    // (bdrSrc); V is the clip reference volume rule (S), independent of
    // mode -- see the function doxygen above for why this row exists.
    {
        const index_t deg = 2*p, n1 = deg+1, K = n1*n1*n1;
        gsVector<real_t> lo(3), hi(3); lo << -1,-1,-1; hi << 1,1,1;
        const gsVector<index_t> nG = gsVector<index_t>::Constant(3, p+1);

        const gsTetClip::BdrCellSource * bdrSrc = ("clip" == mode)
            ? static_cast<const gsTetClip::BdrCellSource*>(&S)
            : static_cast<const gsTetClip::BdrCellSource*>(P.bdr.get());

        std::vector<gsTetClip::KahanSum> accB((size_t)(4*K)), accV((size_t)(3*K));

        for (size_t id = 0; id != N3; ++id)
        {
            index_t i,j,k; idx.ijk(id, i,j,k);
            gsVector<real_t> lower(3), upper(3);
            lower << idx.X[i], idx.Y[j], idx.Z[k];
            upper << idx.X[i+1], idx.Y[j+1], idx.Z[k+1];

            gsMatrix<real_t> bn, bnrm; gsVector<real_t> bw;
            bdrSrc->bdrRule(id, bn, bw, bnrm);
            if (bn.cols() > 0)
            {
                gsVector<real_t> mCell;
                cellMoments(bn, bw, &bnrm, lo, hi, p, mCell);
                for (index_t e = 0; e != 4*K; ++e) accB[(size_t)e].add(mCell[e]);
            }

            gsMatrix<real_t> nd; gsVector<real_t> w;
            if (gsTetClip::Cut == idx.status[id])
                S.volRule(id, nd, w);
            else if (gsTetClip::Full == idx.status[id])
                gsGaussRule<real_t>(nG).mapTo(lower, upper, nd, w);
            if (nd.cols() > 0)
            {
                gsVector<real_t> mCell;
                volDerivMoments(nd, w, lo, hi, p, mCell);
                for (index_t e = 0; e != 3*K; ++e) accV[(size_t)e].add(mCell[e]);
            }
        }

        gsVector<real_t> factor;
        legendreRawFactor(lo, hi, deg, factor);

        real_t worst = 0;
        for (index_t ii = 0; ii != 3; ++ii)
            for (index_t k = 0; k != K; ++k)
            {
                const real_t B = accB[(size_t)((ii+1)*K+k)].value();
                const real_t V = accV[(size_t)(ii*K+k)].value();
                worst = math::max(worst, math::abs(B-V)/(factor[k]*A_mesh));
            }

        GateRow row; row.check = "fluxmom_global"; row.value = worst; row.ref = 0;
        row.tol = tolLoose; row.kind = GateRow::Required; row.relerr = worst;
        row.pass = (worst <= tolLoose);
        G.rows.push_back(row);
        G.requiredPass = G.requiredPass && row.pass;
    }

    if ("clip" != mode)
    {
        // moments_vol: worst-case relative moment error over Cut cells.
        {
            real_t worst = 0; std::vector<size_t> fail;
            for (size_t id = 0; id != N3; ++id)
                if (gsTetClip::Cut == idx.status[id])
                {
                    const real_t e = P.log[id].momErrVol;
                    worst = math::max(worst, e);
                    if (!(e <= tolLoose)) fail.push_back(id);
                }
            GateRow row; row.check = "moments_vol"; row.value = worst; row.ref = 0;
            row.tol = tolLoose; row.kind = GateRow::Required; row.relerr = worst;
            row.failCells = fail; row.pass = (worst <= tolLoose) && fail.empty();
            G.rows.push_back(row);
            G.requiredPass = G.requiredPass && row.pass;
        }
        // moments_bdr: same, over cells with a boundary rule.
        {
            real_t worst = 0; std::vector<size_t> fail;
            for (size_t id = 0; id != N3; ++id)
                if (P.log[id].nBdrIn > 0)
                {
                    const real_t e = P.log[id].momErrBdr;
                    worst = math::max(worst, e);
                    if (!(e <= tolLoose)) fail.push_back(id);
                }
            GateRow row; row.check = "moments_bdr"; row.value = worst; row.ref = 0;
            row.tol = tolLoose; row.kind = GateRow::Required; row.relerr = worst;
            row.failCells = fail; row.pass = (worst <= tolLoose) && fail.empty();
            G.rows.push_back(row);
            G.requiredPass = G.requiredPass && row.pass;
        }
        // fluxmom_cell: worst-case E_K (CellLog::fluxErrBdr, computed in
        // buildCellRules) over Cut cells with a nonzero clip boundary area;
        // cells with A_K == 0 exactly (CellLog::hasFluxErrBdr false) are
        // skipped, not counted as pass or fail.
        {
            real_t worst = 0; std::vector<size_t> fail;
            for (size_t id = 0; id != N3; ++id)
                if (gsTetClip::Cut == idx.status[id] && P.log[id].hasFluxErrBdr)
                {
                    const real_t e = P.log[id].fluxErrBdr;
                    worst = math::max(worst, e);
                    if (!(e <= tolLoose)) fail.push_back(id);
                }
            GateRow row; row.check = "fluxmom_cell"; row.value = worst; row.ref = 0;
            row.tol = tolLoose; row.kind = GateRow::Required; row.relerr = worst;
            row.failCells = fail; row.pass = (worst <= tolLoose) && fail.empty();
            G.rows.push_back(row);
            G.requiredPass = G.requiredPass && row.pass;
        }
    }

    if ("tchakaloff" == mode)
    {
        // residual_vol
        {
            real_t worst = 0; std::vector<size_t> fail;
            for (size_t id = 0; id != N3; ++id)
                if (gsTetClip::Cut == idx.status[id])
                {
                    const CellLog & L = P.log[id];
                    worst = math::max(worst, L.resVol);
                    if (!L.okVol || !(L.resVol < tolTight)) fail.push_back(id);
                }
            GateRow row; row.check = "residual_vol"; row.value = worst; row.ref = 0;
            row.tol = tolTight; row.kind = GateRow::Required; row.relerr = worst;
            row.failCells = fail; row.pass = fail.empty();
            G.rows.push_back(row);
            G.requiredPass = G.requiredPass && row.pass;
        }
        // residual_bdr
        {
            real_t worst = 0; std::vector<size_t> fail;
            for (size_t id = 0; id != N3; ++id)
                if (P.log[id].nBdrIn > 0)
                {
                    const CellLog & L = P.log[id];
                    worst = math::max(worst, L.resBdr);
                    if (!L.okBdr || !(L.resBdr < tolTight)) fail.push_back(id);
                }
            GateRow row; row.check = "residual_bdr"; row.value = worst; row.ref = 0;
            row.tol = tolTight; row.kind = GateRow::Required; row.relerr = worst;
            row.failCells = fail; row.pass = fail.empty();
            G.rows.push_back(row);
            G.requiredPass = G.requiredPass && row.pass;
        }
    }

    // minweight_vol / minweight_bdr: Report except tchakaloff (Req
    // positive). The positivity failCells are read straight from the
    // served tables (`!(w > 0)` per cell, so NaN fails and 0 fails)
    // rather than from CellLog::negW*, which only counts weights < 0
    // and would silently pass a NaN weight or leave failCells empty for
    // a served weight of exactly 0. A cell whose pass threw before
    // storing its rule leaves an empty (0-column) table entry, so it
    // never adds to failVol/failBdr here -- it is caught by the run row
    // instead.
    {
        real_t minVol = std::numeric_limits<real_t>::infinity();
        real_t minBdr = std::numeric_limits<real_t>::infinity();
        for (size_t w = 0; w != P.work.size(); ++w)
        {
            const CellLog & L = P.log[P.work[w]];
            minVol = math::min(minVol, L.minWVol);
            minBdr = math::min(minBdr, L.minWBdr);
        }
        const bool positiveMode = ("tchakaloff" == mode);

        std::vector<size_t> failVol, failBdr;
        if (positiveMode)
        {
            auto anyNonPositive = [](const gsVector<real_t> & wts) -> bool
            {
                for (index_t c = 0; c != wts.size(); ++c)
                    if (!(wts[c] > 0)) return true;
                return false;
            };
            for (size_t w = 0; w != P.work.size(); ++w)
            {
                const size_t id = P.work[w];
                if (gsTetClip::Cut == idx.status[id] && anyNonPositive(P.vol->cell[id].weights))
                    failVol.push_back(id);
                if (anyNonPositive(P.bdr->cell[id].weights))
                    failBdr.push_back(id);
            }
        }

        GateRow rv; rv.check = "minweight_vol"; rv.value = minVol; rv.ref = 0; rv.relerr = 0;
        rv.kind = positiveMode ? GateRow::Required : GateRow::Report; rv.positive = positiveMode;
        if (positiveMode) { rv.failCells = failVol; rv.pass = failVol.empty(); }
        G.rows.push_back(rv);
        if (positiveMode) G.requiredPass = G.requiredPass && rv.pass;

        GateRow rb; rb.check = "minweight_bdr"; rb.value = minBdr; rb.ref = 0; rb.relerr = 0;
        rb.kind = positiveMode ? GateRow::Required : GateRow::Report; rb.positive = positiveMode;
        if (positiveMode) { rb.failCells = failBdr; rb.pass = failBdr.empty(); }
        G.rows.push_back(rb);
        if (positiveMode) G.requiredPass = G.requiredPass && rb.pass;
    }

    // run: cells whose pass threw.
    {
        std::vector<size_t> fail;
        for (size_t id = 0; id != N3; ++id)
            if (!P.log[id].error.empty()) fail.push_back(id);
        GateRow row; row.check = "run"; row.value = (real_t)fail.size(); row.ref = 0;
        row.tol = 0; row.kind = GateRow::Required; row.relerr = row.value;
        row.failCells = fail; row.pass = fail.empty();
        G.rows.push_back(row);
        G.requiredPass = G.requiredPass && row.pass;
    }

    return G;
}

/// Max columns per Gmap evaluation in gatePhysical; bounds that call's basis-evaluation temporaries.
static const index_t GATE_MAP_CHUNK = 8192;

/// Evaluates \a Gmap on \a nodes in consecutive column blocks of at most
/// \ref GATE_MAP_CHUNK, writing into the pre-sized outputs \a detJ/\a X.
/// \a detJ != nullptr: chunked `gsGeometry::compute`, one `NEED_VALUE|
/// NEED_DERIV` call per block when \a X != nullptr too (so det(J) and the
/// physical points come from the SAME evaluation) or `NEED_DERIV` alone
/// otherwise; `(*detJ)[c]` the determinant of the fixed-size
/// `gsMatrix<real_t,3,3>` \a fd.jacobian(c) (cofactor path,
/// `gsFuncData.h:382-387`: column j = dG/du_j), `X` from `fd.values[0]`
/// (`gsGeometry.hpp:595,604`: numCo x numPt = 3 x len). \a detJ ==
/// nullptr: \a X must be non-null; chunked `gsGeometry::eval_into`, the
/// sequential per-point accumulation (`gsBasis.hpp:38-53,151-181`) \ref
/// gatePhysical's boundary points use, since their Jacobian is produced
/// inside \ref gsTetClip::PullbackBdrSource and is unreachable from here
/// -- chunking this path changes nothing but memory, never the
/// arithmetic. Each block is copied into a local `gsMatrix<real_t>`
/// first: both `compute` and `eval_into` take a `const gsMatrix<T>&`,
/// and a `middleCols` expression is not one.
static void gateMapChunked(const gsGeometry<real_t> & Gmap, const gsMatrix<real_t> & nodes,
                           std::vector<real_t> * detJ, gsMatrix<real_t> * X)
{
    const index_t N = nodes.cols();
    if (nullptr != detJ) detJ->resize((size_t)N);
    if (nullptr != X) X->resize(3, N);

    if (nullptr != detJ)
    {
        for (index_t c0 = 0; c0 < N; c0 += GATE_MAP_CHUNK)
        {
            const index_t len = math::min(GATE_MAP_CHUNK, N-c0);
            const gsMatrix<real_t> block = nodes.middleCols(c0, len);
            gsFuncData<real_t> fd((nullptr != X) ? (NEED_DERIV | NEED_VALUE) : NEED_DERIV);
            Gmap.compute(block, fd);
            for (index_t k = 0; k != len; ++k)
            {
                const gsMatrix<real_t,3,3> J = fd.jacobian(k);
                (*detJ)[(size_t)(c0+k)] = J.determinant();
            }
            if (nullptr != X) X->middleCols(c0, len) = fd.values[0];
        }
    }
    else if (nullptr != X)
    {
        for (index_t c0 = 0; c0 < N; c0 += GATE_MAP_CHUNK)
        {
            const index_t len = math::min(GATE_MAP_CHUNK, N-c0);
            const gsMatrix<real_t> block = nodes.middleCols(c0, len);
            gsMatrix<real_t> Xblk;
            Gmap.eval_into(block, Xblk);
            X->middleCols(c0, len) = Xblk;
        }
    }
}

/// Per-cell partial results of \ref gatePhysical's parallel loop: written
/// once per cell id by that loop, then reduced serially, in increasing
/// id, into the shared accumulators. `hasVol`/`hasRef`/`hasBdr` gate
/// which fields are meaningful; the moment vectors stay empty otherwise,
/// so a fully-populated table costs O(N3) scalars plus O(K) per
/// CONTRIBUTING cell (K = (2p+1)^3), not O(N3*K).
struct GatePhysicalCellPartial
{
    bool hasVol = false;
    real_t cellV = 0;
    real_t detMin = std::numeric_limits<real_t>::infinity();
    bool bad = false;

    bool hasRef = false;
    gsVector<real_t> mV;

    bool hasBdr = false;
    real_t cellA = 0, cellF = 0;
    gsVector<real_t> mB;

    real_t excess = 0;

    std::exception_ptr err;
};

/// Physical-space gate for `--mode clip` on a curved background
/// (`--geo affine|bubble`): the analogue of \ref gateTetMode for a case
/// where `CellLog::volSum/bdrSum/fluxSum` (built by \ref buildCellRules
/// straight from `S.volRule`/`S.bdrRule`) are PARAMETRIC and therefore
/// wrong here -- they lack det(J), the Nanson factor and G(u).n. This
/// function recomputes every quantity itself: the cells run in parallel
/// (`#pragma omp parallel for schedule(dynamic)`) into per-cell partial
/// storage (\ref GatePhysicalCellPartial), followed by a SERIAL reduction
/// over cell ids 0..n^3-1 in increasing order, so the result is
/// independent of the OMP thread count and schedule -- the same
/// two-level Kahan discipline (a per-cell KahanSum, then a KahanSum over
/// cells) \ref tetModeTotals uses. Every GISMO_ENSURE in this function
/// runs OUTSIDE the parallel region -- the two entry checks (\a kind,
/// \a st.nVerts) before it, the G(box) check after it -- and so stays
/// catchable by the caller: thrown from inside an omp-parallel region it
/// would call std::terminate instead. An exception raised inside the
/// parallel loop itself is
/// caught per cell into that cell's own `std::exception_ptr`, and the one
/// from the lowest failing id is rethrown, serially, once the loop has
/// finished. Must run BEFORE assembly.
///
/// Pull-back rows (`pb_finite/pb_residual/pb_inside/pb_orientation`):
/// \a st's own per-vertex diagnostics from \ref gsTetClip::pullBack,
/// Required regardless of \a kind -- a bad pull-back poisons every later
/// row. `pb_residual`'s tolerance is against \a st.maxResidual/scale, not
/// an absolute number, so it stays meaningful across the two mesh scales
/// (`--mesh`); a NaN or infinite residual divides through to a NaN/inf
/// relerr, which fails the `<=` test. `detJ_min` uses the SERVED volume
/// rule's own nodes (the ones the solve/gate rule actually visits), one
/// evaluation of det(J) per node via \ref gateMapChunked's chunked
/// evaluation (one `NEED_VALUE|NEED_DERIV` `gsGeometry::compute` call when
/// the physical points are also needed there, `NEED_DERIV` alone
/// otherwise), det(J) from a fixed-size `gsMatrix<real_t,3,3>` (the
/// cofactor determinant path, `fd.jacobian(c)`, `gsFuncData.h:382-387`)
/// rather than the dynamic-size LU path.
///
/// Boundary weights and normals come from \a physBdr
/// (`gsTetClip::PullbackBdrSource`), Nanson's formula n dA = det(J)
/// J^{-T} N dA0 applied to \a S's own parametric boundary rule -- see
/// `gsImmersedPullback.h`'s own doxygen for the cofactor form it uses.
///
/// `fluxmom_global` is the divergence-theorem identity on Omega_s =
/// G(parametric polyhedron): for the tensor Legendre basis q mapped to
/// [lo_phys,hi_phys] (deg = 2p, K = (2p+1)^3),
///   oint_{dOmega_s} q n_i dS = int_{Omega_s} dq/dx_i dV,
/// the left side (`accB`) built from \a physBdr's own served rule (what
/// the Nitsche flux term actually uses), the right side (`accV`) from a
/// SEPARATE reference volume rule, independent of the served one.
///
/// Full-cell reference rule order 3p+1, not the served p+1: sufficient,
/// and derived here, for `--geo identity|affine`, where G is affine
/// (G_j(u) = sum_i A_ji u_i + b_j, degree 1 in every u_i). For such G, a
/// tensor Legendre basis function q_k(x,y,z) = L_kx(x) L_ky(y) L_kz(z)
/// (kx,ky,kz <= 2p) composed with G has, in u_i alone, degree
/// <= kx+ky+kz <= 6p: each of L_kx(G_1(u)), L_ky(G_2(u)), L_kz(G_3(u))
/// contributes at most kx, ky, kz respectively to the u_i-degree, since
/// G_1, G_2, G_3 are each LINEAR in u_i. Differentiating w.r.t. one u_i
/// (chain-ruled through G, itself affine, so dG/du adds no extra degree)
/// costs one degree: d/dx_i [q(G(u))] has per-direction degree <= 6p-1 in
/// u. Tensor Gauss with m points per direction is exact to per-direction
/// degree 2m-1, so m = 3p+1 is exact (2m-1 = 6p+1 >= 6p-1); the served
/// p+1 rule is exact only for G = I, where the same argument gives degree
/// <= 2p. Cut cells need no separate reference rule under identity/
/// affine: the clip volume rule is exact to TOTAL degree 6p (`volDegree`)
/// >= the total degree of d/dx_i[q o G] there, and the clip boundary rule
/// (`bdrDegree`, also 6p) is exact against total deg(q o G) <= 6p with
/// cof(J) constant -- so for `--geo affine` every fluxmom term is exact
/// and the row is a rounding-level check (observed ~1e-16 to ~1e-17 for
/// both cases at p = 1..3).
///
/// The bubble map is NOT affine: G_j(u) = u_j + eps*b(u)*chat_j has
/// degree 2, not 1, in every u_i (from b(u) = prod_i(1-u_i^2)), so the
/// argument above does not apply -- L_kx(G_1(u)) alone already reaches
/// per-direction degree up to 2*kx, and neither the 3p+1 Full-cell rule
/// nor the clip volume/boundary rule (whose own exactness is stated in
/// the PARAMETRIC polyhedron's degree, not in q o G's degree once G is
/// curved) is exact for it. `fluxmom_global` therefore measures genuine,
/// if often numerically small for the tested eps, TRUNCATION error under
/// `--geo bubble`, split by p instead of covered by one tolerance: its
/// per-p tier and tolerance (`BUBBLE_FLUXMOM_REQUIRED_MIN_P`,
/// `BUBBLE_FLUXMOM_TOL`, `BUBBLE_FLUXMOM_TOL_REQUIRED`) come from a
/// measured sweep instead of being derived, documented at the constants
/// themselves.
///
/// Physical box for q. \a lo_phys/\a hi_phys are the componentwise min/max
/// of \a Gmap evaluated at \a Gmap.support()'s 8 corners. This is exactly
/// G(box) in two cases that cover both maps this driver serves: an affine
/// G, whose extreme points are always its corners' images; and any G that
/// fixes the box boundary pointwise (the bubble does, by construction --
/// its displacement is proportional to prod(1-x_i^2), which vanishes on
/// every face of [-1,1]^3), so G(box) = box exactly regardless of the
/// interior deformation. `||q||_inf = 1` (the normalisation
/// `legendreRawFactor` divides through by) holds only inside that box, so
/// every point q is evaluated at (the fluxmom reference-rule volume nodes
/// and the boundary nodes) is checked against it: the largest excess over
/// the whole pass must be at rounding level, GISMO_ENSUREd against
/// 1e-12*max(1,||hi_phys-lo_phys||_inf) after the loop.
///
/// Complexity: one extra volume clip (`S.volRule`) and one extra boundary
/// clip (`physBdr.bdrRule`) per cell -- this function does not reuse any
/// rule already served elsewhere -- plus a (3p+1)^3-point Gauss rule on
/// every Full cell; O(N*K) time and O(K) memory for the moment
/// accumulators (K = (2p+1)^3), matching \ref cellMoments/
/// \ref volDerivMoments's own complexity note. The geometry-map
/// evaluation temporaries (\ref gateMapChunked) are bounded by
/// O(GATE_MAP_CHUNK) columns per call, independent of a cell's own point
/// count; the per-cell partial storage (\ref GatePhysicalCellPartial) is
/// O(K) per contributing cell and O(1) for every other cell.
///
/// Driver-level diagnostic: parallel over cells (OpenMP); must not be called from inside a parallel region.
GateResult gatePhysical(const gsTetClip::PhysTetMesh & phys, const gsTetClip::PullbackStats & st,
                        const gsGeometry<real_t> & Gmap, gsTetClip::BgMapKind kind,
                        const gsTetClip::ClipStreamer & S,
                        const gsTetClip::BdrCellSource & physBdr,
                        index_t p, real_t V_exact, real_t A_exact,
                        const std::string & caseName, const std::string & geoName, index_t r)
{
    GISMO_ENSURE(gsTetClip::BgMapKind::Identity != kind,
                "gatePhysical: called with an identity background map (use gateTetMode instead).");
    GISMO_ENSURE(st.nVerts == (index_t)phys.m.P.size(),
                "gatePhysical: PullbackStats.nVerts does not match the physical mesh's own vertex count.");

    const gsTetClip::CellIndex & idx = *S.index();
    const size_t N3 = (size_t)idx.grid.n*(size_t)idx.grid.n*(size_t)idx.grid.n;

    GateResult G;
    G.caseName = caseName; G.mode = "clip"; G.geo = geoName; G.r = r; G.n = idx.grid.n;

    const real_t V_mesh = gsTetClip::meshVolumeExact(phys);
    const real_t A_mesh = gsTetClip::unclippedBoundaryArea(phys);

    real_t scale = 1.0;
    for (const gsTetClip::Vec3 & x : phys.m.P)
        scale = math::max(scale, math::max(std::abs(x[0]), math::max(std::abs(x[1]), std::abs(x[2]))));

    // Physical box for q: G(box)'s corners (see the doxygen above for why
    // the corners suffice for both --geo affine and --geo bubble).
    const gsMatrix<real_t> parBox = Gmap.support();   // 3x2: row d, col 0/1 = min/max
    gsMatrix<real_t> corners(3,8);
    {
        index_t c = 0;
        for (int cz = 0; cz != 2; ++cz)
        for (int cy = 0; cy != 2; ++cy)
        for (int cx = 0; cx != 2; ++cx)
        {
            corners(0,c) = parBox(0,cx); corners(1,c) = parBox(1,cy); corners(2,c) = parBox(2,cz);
            ++c;
        }
    }
    gsMatrix<real_t> cornersPhys;
    Gmap.eval_into(corners, cornersPhys);
    gsVector<real_t> lo_phys(3), hi_phys(3);
    for (index_t d = 0; d != 3; ++d)
    {
        lo_phys[d] = cornersPhys(d,0); hi_phys[d] = cornersPhys(d,0);
        for (index_t c = 1; c != 8; ++c)
        {
            lo_phys[d] = math::min(lo_phys[d], cornersPhys(d,c));
            hi_phys[d] = math::max(hi_phys[d], cornersPhys(d,c));
        }
    }

    const index_t deg = 2*p, n1 = deg+1, K = n1*n1*n1;
    std::vector<gsTetClip::KahanSum> accV((size_t)(3*K)), accB((size_t)(4*K));
    gsTetClip::KahanSum totV, totA, totF;
    real_t detMin = std::numeric_limits<real_t>::infinity();
    std::vector<size_t> detFail;
    real_t maxExcess = 0;

    // fullRuleP1 is the SAME order the solve/assembly rule serves Full
    // cells (solvePoissonOnDomain's own makeVolLookupFactory(idx, volSrc,
    // Constant(3, p+1)) call).
    const gsGaussRule<real_t> fullRuleP1(gsVector<index_t>::Constant(3, p+1));
    const gsGaussRule<real_t> fullRule3p1(gsVector<index_t>::Constant(3, 3*p+1));

    std::vector<GatePhysicalCellPartial> part(N3);

    const long N3L = (long)N3;
    #pragma omp parallel for schedule(dynamic)
    for (long idL = 0; idL < N3L; ++idL)
    {
        const size_t id = (size_t)idL;
        GatePhysicalCellPartial & C = part[id];
        try
        {
            index_t i,j,k; idx.ijk(id, i,j,k);
            gsVector<real_t> lower(3), upper(3);
            lower << idx.X[i], idx.Y[j], idx.Z[k];
            upper << idx.X[i+1], idx.Y[j+1], idx.Z[k+1];

            const bool isCut  = (gsTetClip::Cut  == idx.status[id]);
            const bool isFull = (gsTetClip::Full == idx.status[id]);

            auto trackExcess = [&](const gsMatrix<real_t> & Xm)
            {
                for (index_t c = 0; c != Xm.cols(); ++c)
                    for (index_t d = 0; d != 3; ++d)
                    {
                        C.excess = math::max(C.excess, lo_phys[d]-Xm(d,c));
                        C.excess = math::max(C.excess, Xm(d,c)-hi_phys[d]);
                    }
            };

            // 1. Served volume rule -> volume total and detJ_min/detFail.
            gsMatrix<real_t> nd; gsVector<real_t> w;
            std::vector<real_t> detJ; gsMatrix<real_t> X;
            if (isCut)       S.volRule(id, nd, w);
            else if (isFull) fullRuleP1.mapTo(lower, upper, nd, w);

            if (nd.cols() > 0)
            {
                if (isCut) gateMapChunked(Gmap, nd, &detJ, &X);
                else       gateMapChunked(Gmap, nd, &detJ, nullptr);

                gsTetClip::KahanSum cellV;
                bool cellBad = false;
                real_t detMinCell = std::numeric_limits<real_t>::infinity();
                for (index_t c = 0; c != nd.cols(); ++c)
                {
                    const real_t dJ = detJ[(size_t)c];
                    detMinCell = math::min(detMinCell, dJ);
                    if (!(dJ > 0)) cellBad = true;
                    cellV.add(w[c]*math::abs(dJ));
                }
                C.hasVol = true;
                C.detMin = detMinCell;
                C.bad = cellBad;
                C.cellV = cellV.value();
            }

            // 2. fluxmom V side: a reference rule, NOT the served one (see
            // the doxygen above for the 3p+1 exactness argument). Cut cells
            // reuse the same clip rule and det(J)/X already computed in
            // step 1, by pointer; the clip rule is already exact to total
            // degree 6p, which dominates.
            const gsMatrix<real_t> * ndRefP = nullptr; const gsVector<real_t> * wRefP = nullptr;
            const std::vector<real_t> * detJrefP = nullptr; const gsMatrix<real_t> * XrefP = nullptr;
            gsMatrix<real_t> ndRefFull; gsVector<real_t> wRefFull;
            std::vector<real_t> detJrefFull; gsMatrix<real_t> XrefFull;
            if (isCut && nd.cols() > 0)
            {
                ndRefP = &nd; wRefP = &w; detJrefP = &detJ; XrefP = &X;
            }
            else if (isFull)
            {
                fullRule3p1.mapTo(lower, upper, ndRefFull, wRefFull);
                if (ndRefFull.cols() > 0)
                {
                    gateMapChunked(Gmap, ndRefFull, &detJrefFull, &XrefFull);
                    ndRefP = &ndRefFull; wRefP = &wRefFull;
                    detJrefP = &detJrefFull; XrefP = &XrefFull;
                }
            }
            if (nullptr != ndRefP && ndRefP->cols() > 0)
            {
                trackExcess(*XrefP);
                gsVector<real_t> wJ(ndRefP->cols());
                for (index_t c = 0; c != ndRefP->cols(); ++c)
                    wJ[c] = (*wRefP)[c]*math::abs((*detJrefP)[(size_t)c]);
                volDerivMoments(*XrefP, wJ, lo_phys, hi_phys, p, C.mV);
                C.hasRef = true;
            }

            // 3. Boundary: physBdr already serves PHYSICAL weights/normals
            // (Nanson); nodes stay parametric, so Xb = Gmap(bn) is needed
            // for the moment evaluation and for the flux dot product.
            // Called for every id regardless of status, as gateTetMode's
            // own fluxmom_global block does: a boundary piece can live in
            // a cell that is not Cut.
            gsMatrix<real_t> bn, bnrm; gsVector<real_t> bw;
            physBdr.bdrRule(id, bn, bw, bnrm);
            if (bn.cols() > 0)
            {
                gsMatrix<real_t> Xb;
                gateMapChunked(Gmap, bn, nullptr, &Xb);
                trackExcess(Xb);
                gsTetClip::KahanSum cellA, cellF;
                for (index_t c = 0; c != bn.cols(); ++c)
                {
                    cellA.add(bw[c]);
                    cellF.add(bw[c]*Xb.col(c).dot(bnrm.col(c)));
                }
                C.hasBdr = true;
                C.cellA = cellA.value();
                C.cellF = cellF.value();
                cellMoments(Xb, bw, &bnrm, lo_phys, hi_phys, p, C.mB);
            }
        }
        catch (...)
        {
            C.err = std::current_exception();
        }
    }

    // Serial reduction, in increasing cell id: bitwise independent of the
    // thread count and schedule. First the lowest-id exception, if any --
    // an exception escaping the omp-parallel region above would call
    // std::terminate instead of being catchable here.
    for (size_t id = 0; id != N3; ++id)
        if (part[id].err) std::rethrow_exception(part[id].err);

    for (size_t id = 0; id != N3; ++id)
    {
        GatePhysicalCellPartial & C = part[id];
        if (C.hasVol)
        {
            detMin = math::min(detMin, C.detMin);
            if (C.bad) detFail.push_back(id);
            totV.add(C.cellV);
        }
        if (C.hasRef)
            for (index_t e = 0; e != 3*K; ++e) accV[(size_t)e].add(C.mV[e]);
        if (C.hasBdr)
        {
            totA.add(C.cellA);
            totF.add(C.cellF);
            for (index_t e = 0; e != 4*K; ++e) accB[(size_t)e].add(C.mB[e]);
        }
        maxExcess = math::max(maxExcess, C.excess);
        C.mV.resize(0); C.mB.resize(0);
    }

    real_t widthInf = 0;
    for (index_t d = 0; d != 3; ++d) widthInf = math::max(widthInf, hi_phys[d]-lo_phys[d]);
    GISMO_ENSURE(maxExcess <= 1e-12*math::max((real_t)1, widthInf),
                "gatePhysical: a physical evaluation point sits outside G(box) by " << maxExcess
                << ", past the exactness-domain tolerance -- lo_phys/hi_phys is not G(box).");

    const real_t V = totV.value(), A = totA.value(), F = totF.value();

    auto push = [&](GateRow row)
    {
        G.rows.push_back(row);
        if (GateRow::Required == row.kind) G.requiredPass = G.requiredPass && row.pass;
    };
    auto scalarRow = [&](const std::string & name, real_t value, real_t ref, real_t tol,
                         GateRow::Kind kind)
    {
        GateRow row; row.check = name; row.value = value; row.ref = ref; row.tol = tol; row.kind = kind;
        row.relerr = math::abs(value-ref)/math::abs(ref);
        row.pass = (row.relerr <= tol);
        push(row);
    };

    {
        GateRow row; row.check = "pb_finite"; row.value = (real_t)st.nFailed; row.ref = 0; row.tol = 0;
        row.kind = GateRow::Required; row.relerr = row.value; row.pass = (0 == st.nFailed);
        push(row);
    }
    {
        GateRow row; row.check = "pb_residual"; row.value = st.maxResidual; row.ref = 0; row.tol = 1e-12;
        row.kind = GateRow::Required; row.relerr = st.maxResidual/scale;
        row.pass = (row.relerr <= 1e-12);
        push(row);
    }
    {
        GateRow row; row.check = "pb_inside"; row.value = (real_t)st.nOutside; row.ref = 0; row.tol = 0;
        row.kind = GateRow::Required; row.relerr = row.value; row.pass = (0 == st.nOutside);
        push(row);
    }
    {
        GateRow row; row.check = "pb_orientation"; row.value = (real_t)st.nInverted; row.ref = 0; row.tol = 0;
        row.kind = GateRow::Required; row.relerr = row.value;
        row.pass = (0 == st.nInverted && st.minParamTetVol > 0);
        push(row);
    }
    {
        GateRow row; row.check = "detJ_min"; row.value = detMin; row.ref = 0; row.relerr = 0; row.tol = 0;
        row.kind = GateRow::Required; row.positive = true; row.failCells = detFail;
        row.pass = detFail.empty() && (detMin > 0);
        push(row);
    }

    const GateRow::Kind physKind = (gsTetClip::BgMapKind::Affine == kind) ? GateRow::Required : GateRow::Report;
    scalarRow("volume", V, V_mesh, 1e-13, physKind);
    scalarRow("area",   A, A_mesh, 1e-13, physKind);
    scalarRow("flux",   F, 3.0*V_mesh, 1e-13, physKind);
    scalarRow("volume_exact", V, V_exact, 0, GateRow::Report);
    scalarRow("area_exact",   A, A_exact, 0, GateRow::Report);

    {
        gsVector<real_t> factor;
        legendreRawFactor(lo_phys, hi_phys, deg, factor);
        real_t worst = 0;
        for (index_t ii = 0; ii != 3; ++ii)
            for (index_t kk = 0; kk != K; ++kk)
            {
                const real_t B = accB[(size_t)((ii+1)*K+kk)].value();
                const real_t Vd = accV[(size_t)(ii*K+kk)].value();
                // `!(e <= worst)` rather than math::max: a NaN entry must
                // propagate into `worst` and fail the row, not be dropped.
                const real_t e = math::abs(B-Vd)/(factor[kk]*A);
                if (!(e <= worst)) worst = e;
            }

        const bool isAffine = (gsTetClip::BgMapKind::Affine == kind);
        const bool bubbleRequired = (p >= BUBBLE_FLUXMOM_REQUIRED_MIN_P);
        const real_t tol = isAffine ? (real_t)1e-12
                          : (bubbleRequired ? BUBBLE_FLUXMOM_TOL_REQUIRED : BUBBLE_FLUXMOM_TOL);
        const GateRow::Kind fmKind = isAffine ? GateRow::Required
                                    : (bubbleRequired ? GateRow::Required : GateRow::Report);
        GateRow row; row.check = "fluxmom_global"; row.value = worst; row.ref = 0;
        row.tol = tol; row.kind = fmKind; row.relerr = worst; row.pass = (worst <= tol);
        push(row);
    }

    gsInfo << "GATEINFO case=" << caseName << " geo=" << geoName << " r=" << r
          << " pullback nVerts=" << st.nVerts << " nFailed=" << st.nFailed
          << " nOutside=" << st.nOutside << " nInverted=" << st.nInverted
          << " maxResidual=" << gsTetClip::fmtSci(st.maxResidual)
          << " scale=" << gsTetClip::fmtSci(scale)
          << " minParamTetVol=" << gsTetClip::fmtSci(st.minParamTetVol)
          << " physBox=[" << gsTetClip::fmtSci(lo_phys[0]) << "," << gsTetClip::fmtSci(hi_phys[0]) << "]x["
          << gsTetClip::fmtSci(lo_phys[1]) << "," << gsTetClip::fmtSci(hi_phys[1]) << "]x["
          << gsTetClip::fmtSci(lo_phys[2]) << "," << gsTetClip::fmtSci(hi_phys[2]) << "]  INFO\n";

    return G;
}

//----------------------------------------------------------------------------
// Algoim mode (analytic sphere).
//----------------------------------------------------------------------------

/// The pulled-back sphere level set phi-hat(u) = phi(G(u)) = |G(u)-c|^2 - R^2
/// (c/R = SPHERE_CENTER_X/Y/Z/SPHERE_RADIUS), composing the physical sphere
/// level set with a background map G so that Algoim's box quadrature
/// (gsAlgoimAdaptiveRule, gsImplicitTrimmedDomain -- both call only
/// eval_into/deriv_into, never hessian_into) runs directly in PARAMETRIC
/// space. \a m_G is stored by reference, NOT owned: it must outlive every
/// use of this function, including any Algoim rule or domain built over it.
/// No mutable state: \ref gsImplicitTrimmedDomain evaluates \ref eval_into
/// concurrently. Neither \ref eval_into nor \ref deriv_into may
/// GISMO_ENSURE -- an omp-parallel caller could evaluate either, and a
/// thrown exception escaping such a region calls std::terminate instead of
/// being catchable (see \ref gatePhysical's own doxygen).
class PulledBackSphere : public gsFunction<real_t>
{
public:
    explicit PulledBackSphere(const gsGeometry<real_t> & G) : m_G(G) {}

    virtual short_t domainDim() const { return 3; }
    virtual short_t targetDim() const { return 1; }

    /// phi-hat(u) = |G(u)-c|^2 - R^2.
    virtual void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const
    {
        gsMatrix<real_t> x;
        m_G.eval_into(u, x);
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const real_t dx = x(0,k)-SPHERE_CENTER_X, dy = x(1,k)-SPHERE_CENTER_Y,
                        dz = x(2,k)-SPHERE_CENTER_Z;
            result(0,k) = dx*dx+dy*dy+dz*dz-SPHERE_RADIUS*SPHERE_RADIUS;
        }
    }

    /// Chain rule: grad(phi-hat)(u) = J(u)^T * grad(phi)(G(u))
    /// = J(u)^T * 2*(G(u)-c), J = dG/du (gsFuncData.h:382-387: column j of
    /// jacobian(k) is dG/du_j).
    virtual void deriv_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const
    {
        gsMatrix<real_t> x;
        m_G.eval_into(u, x);
        gsFuncData<real_t> fd(NEED_DERIV);
        m_G.compute(u, fd);
        result.resize(3, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            gsVector<real_t> g(3);
            g << 2.0*(x(0,k)-SPHERE_CENTER_X), 2.0*(x(1,k)-SPHERE_CENTER_Y), 2.0*(x(2,k)-SPHERE_CENTER_Z);
            result.col(k) = fd.jacobian(k).transpose()*g;
        }
    }

    GISMO_CLONE_FUNCTION(PulledBackSphere)

private:
    const gsGeometry<real_t> & m_G;
};

/// Rigorous global bound L-hat >= max_u |grad(phi-hat)(u)| over G's support
/// box, for the reference rule's box classifier ONLY (\ref
/// gsAlgoimAdaptiveRule's LipschitzConstant is read solely by classify(),
/// which runs only when the recursion subdivides, i.e. only for
/// maxDepth > 0 -- both solve rules in \ref runAlgoimPass use maxDepth=0
/// and are unaffected by it). From
/// |grad(phi-hat)(u)| = |J(u)^T * 2(G(u)-c)| <= sigma_max(J(u)) * 2|G(u)-c|:
///  - \a kind == Affine (J = A constant): the |Au+b-c| factor's box maximum
///    is attained at a corner (convex function), so max_corner|G(corner)-c|
///    is the EXACT box maximum of |G(u)-c|; the sigma_max(A) factor is only
///    a submultiplicative bound (|A^T v| <= sigma_max(A)|v|), tight only
///    when the maximizing corner's G(corner)-c aligns with A's top
///    right-singular vector. L-hat is therefore rigorous but, in general,
///    not the exact maximum of |grad(phi-hat)| itself.
///  - \a kind == Bubble (G(u) = u + eps*b(u)*chat, b(u) = prod_i(1-u_i^2),
///    |chat| = 1): sigma_max(I + eps*chat*grad(b)^T) <= 1 + |eps|*max|grad(b)|
///    = 1 + 2|eps| (|grad_i(b)| <= 2 on the box, attained at a face centre --
///    see the file header's Nitsche-h paragraph), and
///    |G(u)-c| <= |u-c| + |eps|*b(u) <= |u-c| + |eps| since 0 <= b <= 1;
///    |u-c|'s box maximum is at a corner (convex), giving
///    L-hat = 2*(1+2|eps|)*(max_corner|corner-c| + |eps|). Both factors are
///    rigorous bounds, not exact.
/// The identity path keeps LipschitzConstant = 4.0 (a valid global bound of
/// phi itself) untouched; this function is never called there.
real_t pulledBackSphereLipschitz(const gsGeometry<real_t> & G, gsTetClip::BgMapKind kind, real_t eps)
{
    GISMO_ENSURE(gsTetClip::BgMapKind::Identity != kind,
                "pulledBackSphereLipschitz: not meaningful for --geo identity.");

    const gsMatrix<real_t> supp = G.support();   // 3x2: row d, col 0/1 = min/max
    gsMatrix<real_t> corners(3,8);
    {
        index_t c = 0;
        for (int cz = 0; cz != 2; ++cz)
        for (int cy = 0; cy != 2; ++cy)
        for (int cx = 0; cx != 2; ++cx)
        {
            corners(0,c) = supp(0,cx); corners(1,c) = supp(1,cy); corners(2,c) = supp(2,cz);
            ++c;
        }
    }

    if (gsTetClip::BgMapKind::Affine == kind)
    {
        gsMatrix<real_t> centre(3,1);
        centre.col(0) = 0.5*(supp.col(0)+supp.col(1));
        gsFuncData<real_t> fd(NEED_DERIV);
        G.compute(centre, fd);
        const gsMatrix<real_t,3,3> A = fd.jacobian(0);
        gsMatrix<real_t>::JacobiSVD svdA(A);
        const real_t sigmaMax = svdA.singularValues().maxCoeff();

        gsMatrix<real_t> X;
        G.eval_into(corners, X);
        real_t maxDist = 0;
        for (index_t k = 0; k != 8; ++k)
        {
            const real_t dx = X(0,k)-SPHERE_CENTER_X, dy = X(1,k)-SPHERE_CENTER_Y,
                        dz = X(2,k)-SPHERE_CENTER_Z;
            maxDist = math::max(maxDist, math::sqrt(dx*dx+dy*dy+dz*dz));
        }
        return 2.0*sigmaMax*maxDist;
    }
    else   // Bubble
    {
        real_t maxCornerDist = 0;
        for (index_t k = 0; k != 8; ++k)
        {
            const real_t dx = corners(0,k)-SPHERE_CENTER_X, dy = corners(1,k)-SPHERE_CENTER_Y,
                        dz = corners(2,k)-SPHERE_CENTER_Z;
            maxCornerDist = math::max(maxCornerDist, math::sqrt(dx*dx+dy*dy+dz*dz));
        }
        return 2.0*(1.0+2.0*math::abs(eps))*(maxCornerDist+math::abs(eps));
    }
}

/// Serial 21^3-tensor-grid (corners included) sample of |grad(phi-hat)|
/// against the rigorous bound \a Lhat from \ref pulledBackSphereLipschitz,
/// run once per (geo, r) BEFORE assembly (a GISMO_ENSURE thrown from inside
/// the omp-parallel element loops would std::terminate instead of being
/// catchable). Prints one ALGOIM-LIPSCHITZ INFO line.
void checkPulledBackSphereLipschitz(const PulledBackSphere & phiHat, const gsGeometry<real_t> & G,
                                    const std::string & geo, index_t r, real_t Lhat)
{
    const gsMatrix<real_t> supp = G.support();
    const index_t nPerDir = 21;
    gsMatrix<real_t> pts(3, nPerDir*nPerDir*nPerDir);
    index_t col = 0;
    for (index_t k = 0; k != nPerDir; ++k)
    for (index_t j = 0; j != nPerDir; ++j)
    for (index_t i = 0; i != nPerDir; ++i)
    {
        pts(0,col) = supp(0,0) + (supp(0,1)-supp(0,0))*(real_t)i/(real_t)(nPerDir-1);
        pts(1,col) = supp(1,0) + (supp(1,1)-supp(1,0))*(real_t)j/(real_t)(nPerDir-1);
        pts(2,col) = supp(2,0) + (supp(2,1)-supp(2,0))*(real_t)k/(real_t)(nPerDir-1);
        ++col;
    }
    gsMatrix<real_t> g;
    phiHat.deriv_into(pts, g);
    real_t maxSampled = 0;
    for (index_t c = 0; c != g.cols(); ++c)
        maxSampled = math::max(maxSampled, g.col(c).norm());

    GISMO_ENSURE(maxSampled <= Lhat, "checkPulledBackSphereLipschitz: sampled |grad(phi-hat)|="
                << maxSampled << " exceeds the bound L=" << Lhat << " at geo=" << geo << " r=" << r);
    gsInfo << "ALGOIM-LIPSCHITZ geo=" << geo << " r=" << r << " L=" << gsTetClip::fmtSci(Lhat)
          << " sampledMax=" << gsTetClip::fmtSci(maxSampled) << "  INFO\n";
}

/// gsAlgoimAdaptiveRule option set of the old runAlgoimBaseline: \a dim,
/// \a quA, \a quB, \a maxDepth over defaultOptions().
gsOptionList algoimOptions(short_t dim, real_t quA = 1.0, index_t quB = 1, index_t maxDepth = 0)
{
    gsOptionList o = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
    o.setInt("dim", dim);
    o.setReal("quA", quA);
    o.setInt("quB", quB);
    o.setInt("maxDepth", maxDepth);
    return o;
}

/// Two-level-Kahan V/A/F totals (see \ref tetModeTotals for why: a
/// per-cell KahanSum, then a KahanSum over cells in increasing id) of one
/// serial Algoim pass over every background cell, plus the minimum weight,
/// point and negative-weight counts over all cells, and per-cell errors.
/// `volStats`/`srfStats` are `gsAlgoimAdaptiveRule::Stats` accumulated
/// across all cells for the volume and surface rules respectively.
/// \a pulledBack, \a minDetJ and \a maxNormalDev are set only on the
/// pulled-back path (\ref runAlgoimPass called with a non-null \c G): the
/// minimum det(J) served over every volume AND surface node, and the
/// largest deviation of the served (Nanson) normal from the analytic
/// (x-c)/|x-c|, which coincide in exact arithmetic since
/// J^-T*nu-hat prop-to J^-T*J^T*grad(phi) = grad(phi).
struct AlgoimPass
{
    real_t V = 0, A = 0, F = 0,
           minWVol = std::numeric_limits<real_t>::infinity(),
           minWBdr = std::numeric_limits<real_t>::infinity();
    index_t volPts = 0, bdrPts = 0, negWVol = 0, negWBdr = 0;
    std::vector<size_t> errorCells;
    std::vector<std::string> errors;
    gsAlgoimAdaptiveRule<real_t>::Stats volStats, srfStats;
    real_t seconds = 0;
    bool pulledBack = false;
    real_t minDetJ = std::numeric_limits<real_t>::infinity(), maxNormalDev = 0;
};

/// Serial loop over ALL n^3 cells (as runAlgoimBaseline): volume rule
/// dim=-1 via mapTo, surface rule dim=3 via mapToSeparated (the cut part is
/// the surface rule). Per-cell Kahan, then id order. gsAlgoimAdaptiveRule is
/// not thread-safe, hence the serial loop. Every cell's calls are wrapped in
/// try/catch: an exception is recorded against the cell id and \a phi is NOT
/// reclassified.
///
/// \a G == nullptr (default): the IDENTITY path -- \a phi is the physical
/// sphere level set, normals are the analytic (x-c)/|x-c| computed here
/// (Algoim itself returns none), and every total is already physical.
///
/// \a G != nullptr: the PULLED-BACK path. \a phi is phi-hat = phi o G on
/// the SAME G (a \ref PulledBackSphere), so Algoim runs in parametric space.
/// The volume table keeps the PARAMETRIC (nodes, weights); the physical
/// volume total is Sigma w*det(J) instead, with \a minDetJ tracked over
/// every served volume node. The boundary table stores the PARAMETRIC
/// normal nu-hat = grad(phi-hat)/|grad(phi-hat)| at the cut nodes (from one
/// \c phi.deriv_into(cut, g) call per cell) -- it is NEVER Nanson-mapped in
/// place, so callers (the poisson study, \ref runVolumeStudy) get back a
/// table still usable as a parametric \c BdrCellSource. The area/flux
/// totals and \a minDetJ/\a maxNormalDev are nonetheless computed here, in
/// PHYSICAL units, by wrapping that same table in a temporary
/// \ref gsTetClip::PullbackBdrSource and re-looping over it in id order
/// (the same two-level Kahan discipline) -- Nanson is therefore evaluated
/// once for this pass's own totals and, independently, once more whenever
/// the caller later wraps the returned table for assembly; the STORED table
/// itself is never overwritten, so no quantity is Nanson-mapped twice.
///
/// \a volTabOut / \a bdrTabOut (both optional, default null): when given, a
/// COPY of every cell's Algoim rule is also stored there (\c cell.size() ==
/// n^3), keyed by the SAME flat id this loop already uses, so the poisson
/// study can back a VolLookupRule/BdrLookupRule with the tables instead of
/// re-running Algoim; \a runVolumeStudy passes neither for the identity
/// path (unaffected). On the pulled-back path a null \a bdrTabOut is backed
/// by a LOCAL table instead (needed regardless, to compute the boundary
/// totals via PullbackBdrSource).
AlgoimPass runAlgoimPass(const gsTetClip::Grid3 & grid, index_t p, const gsFunction<real_t> & phi,
                         gsTetClip::VolCellTable * volTabOut = nullptr,
                         gsTetClip::BdrCellTable * bdrTabOut = nullptr,
                         const gsGeometry<real_t> * G = nullptr)
{
    AlgoimPass A;
    A.pulledBack = (nullptr != G);
    gsAlgoimAdaptiveRule<real_t> volRule(phi, (short_t)p, algoimOptions(-1));
    gsAlgoimAdaptiveRule<real_t> srfRule(phi, (short_t)p, algoimOptions(3));

    std::vector<real_t> X, Y, Z;
    gsTetClip::gridLines(grid, X, Y, Z);
    const index_t n = grid.n;
    const size_t nCells = (size_t)n*(size_t)n*(size_t)n;

    gsTetClip::BdrCellTable localBdrTab;
    gsTetClip::BdrCellTable * bdrTabUse = bdrTabOut;
    if (nullptr != G && nullptr == bdrTabUse) bdrTabUse = &localBdrTab;

    if (nullptr != volTabOut) volTabOut->cell.resize(nCells);
    if (nullptr != bdrTabUse) bdrTabUse->cell.resize(nCells);

    gsStopwatch sw;
    gsTetClip::KahanSum totV, totA, totF;
    real_t minDetJ = std::numeric_limits<real_t>::infinity();

    for (index_t k = 0; k != n; ++k)
    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const size_t id = (size_t)i + (size_t)n*((size_t)j + (size_t)n*(size_t)k);
        try
        {
            gsVector<real_t> lower(3), upper(3);
            lower << X[i], Y[j], Z[k];
            upper << X[i+1], Y[j+1], Z[k+1];

            gsMatrix<real_t> vn; gsVector<real_t> vw;
            volRule.mapTo(lower, upper, vn, vw);
            A.volPts += vn.cols();
            gsTetClip::KahanSum cellV;
            if (nullptr == G)
            {
                for (index_t c = 0; c != vw.size(); ++c)
                {
                    cellV.add(vw[c]);
                    if (vw[c] < A.minWVol) A.minWVol = vw[c];
                    if (vw[c] < 0) ++A.negWVol;
                }
            }
            else
            {
                gsFuncData<real_t> fd(NEED_DERIV);
                G->compute(vn, fd);
                for (index_t c = 0; c != vw.size(); ++c)
                {
                    const real_t dJ = fd.jacobian(c).determinant();
                    minDetJ = math::min(minDetJ, dJ);
                    const real_t wPhys = vw[c]*dJ;
                    cellV.add(wPhys);
                    if (wPhys < A.minWVol) A.minWVol = wPhys;
                    if (wPhys < 0) ++A.negWVol;
                }
            }
            totV.add(cellV.value());
            if (nullptr != volTabOut) volTabOut->cell[id] = gsTetClip::CellRule3{ vn, vw };

            gsMatrix<real_t> interior, cut; gsVector<real_t> interiorW, cutW;
            srfRule.mapToSeparated(lower, upper, interior, interiorW, cut, cutW);
            A.bdrPts += cut.cols();

            if (nullptr == G)
            {
                gsMatrix<real_t> cutNormals(3, cut.cols());
                gsTetClip::KahanSum cellA, cellF;
                for (index_t c = 0; c != cutW.size(); ++c)
                {
                    cellA.add(cutW[c]);
                    const real_t nx = cut(0,c)-SPHERE_CENTER_X, ny = cut(1,c)-SPHERE_CENTER_Y,
                                 nz = cut(2,c)-SPHERE_CENTER_Z;
                    const real_t nrm = math::sqrt(nx*nx+ny*ny+nz*nz);
                    GISMO_ENSURE(nrm > 0, "runAlgoimPass: zero-radius surface node.");
                    cellF.add(cutW[c]*(cut(0,c)*nx + cut(1,c)*ny + cut(2,c)*nz)/nrm);
                    if (cutW[c] < A.minWBdr) A.minWBdr = cutW[c];
                    if (cutW[c] < 0) ++A.negWBdr;
                    if (nullptr != bdrTabOut)
                    { cutNormals(0,c) = nx/nrm; cutNormals(1,c) = ny/nrm; cutNormals(2,c) = nz/nrm; }
                }
                totA.add(cellA.value());
                totF.add(cellF.value());
                if (nullptr != bdrTabOut)
                    bdrTabOut->cell[id] = gsTetClip::CellBdrRule3{ cut, cutW, cutNormals };
            }
            else
            {
                // Parametric normal nu-hat = grad(phi-hat)/|grad(phi-hat)|:
                // the table stores it unchanged, PARAMETRIC -- Nanson is
                // applied only when a PullbackBdrSource later wraps this
                // table, never by writing into the table itself.
                // minWBdr/negWBdr are tracked below, over the Nanson-mapped
                // (physical) weight, not over the parametric cutW here.
                gsMatrix<real_t> g;
                phi.deriv_into(cut, g);
                gsMatrix<real_t> cutNormals(3, cut.cols());
                for (index_t c = 0; c != cutW.size(); ++c)
                {
                    const real_t gn = g.col(c).norm();
                    GISMO_ENSURE(gn > 0, "runAlgoimPass: zero gradient at a pulled-back surface node.");
                    cutNormals.col(c) = g.col(c)/gn;
                }
                bdrTabUse->cell[id] = gsTetClip::CellBdrRule3{ cut, cutW, cutNormals };
            }
        }
        catch (const std::exception & e)
        {
            A.errorCells.push_back(id);
            A.errors.push_back(e.what());
        }
    }
    A.seconds = sw.stop();
    A.volStats = volRule.stats();
    A.srfStats = srfRule.stats();

    if (nullptr == G)
    {
        A.V = totV.value(); A.A = totA.value(); A.F = totF.value();
    }
    else
    {
        // Wrap the stored PARAMETRIC table in a PullbackBdrSource and
        // re-loop in id order, the SAME two-level Kahan discipline used
        // above, now over PHYSICAL weights/normals: minWBdr/negWBdr are
        // tracked here (over w, physical), not in the loop above (cutW,
        // parametric).
        gsTetClip::PullbackBdrSource pb(*bdrTabUse, *G);
        real_t maxNormalDev = 0;
        for (size_t id = 0; id != nCells; ++id)
        {
            gsMatrix<real_t> nodes, nrm; gsVector<real_t> w;
            pb.bdrRule(id, nodes, w, nrm);
            if (0 == nodes.cols()) continue;

            gsMatrix<real_t> x;
            G->eval_into(nodes, x);
            gsFuncData<real_t> fd(NEED_DERIV);
            G->compute(nodes, fd);

            gsTetClip::KahanSum cellA, cellF;
            for (index_t c = 0; c != nodes.cols(); ++c)
            {
                minDetJ = math::min(minDetJ, fd.jacobian(c).determinant());
                if (w[c] < A.minWBdr) A.minWBdr = w[c];
                if (w[c] < 0) ++A.negWBdr;
                cellA.add(w[c]);
                cellF.add(w[c]*x.col(c).dot(nrm.col(c)));

                const real_t nx = x(0,c)-SPHERE_CENTER_X, ny = x(1,c)-SPHERE_CENTER_Y,
                            nz = x(2,c)-SPHERE_CENTER_Z;
                const real_t nrmA = math::sqrt(nx*nx+ny*ny+nz*nz);
                gsVector<real_t> analytic(3); analytic << nx/nrmA, ny/nrmA, nz/nrmA;
                maxNormalDev = math::max(maxNormalDev, (nrm.col(c)-analytic).norm());
            }
            totA.add(cellA.value());
            totF.add(cellF.value());
        }
        A.V = totV.value(); A.A = totA.value(); A.F = totF.value();
        A.minDetJ = minDetJ;
        A.maxNormalDev = maxNormalDev;
    }

    return A;
}

/// Gate of the Algoim baseline. Report-only except `run` (and, on the
/// pulled-back path, `detJ_min`/`normal_dev`): this rule runs on EVERY cell
/// (like runAlgoimBaseline), unrelated to the clip status or to a later
/// trimmed-domain classification, so its volume/area/flux rows are
/// informational, not a correctness gate -- \a P.V/A/F already carry
/// PHYSICAL totals on both paths (\ref runAlgoimPass), compared against the
/// same analytic sphere either way. `detJ_min`/`normal_dev` exist only when
/// \a P.pulledBack: det(J) > 0 at every served node, and the served (Nanson)
/// normal agrees with the analytic one to rounding.
GateResult gateAlgoim(const AlgoimPass & P, index_t p, const std::string & caseName, index_t r, index_t n)
{
    GISMO_UNUSED(p);
    GateResult G; G.caseName = caseName; G.mode = "algoim"; G.r = r; G.n = n;

    const real_t V_exact = sphereVolumeExact(), A_exact = sphereAreaExact(), F_exact = 3.0*V_exact;

    auto push = [&](GateRow row)
    {
        if (GateRow::Required == row.kind) G.requiredPass = G.requiredPass && row.pass;
        G.rows.push_back(row);
    };
    auto reportRow = [&](const std::string & name, real_t value, real_t ref)
    {
        GateRow row; row.check = name; row.value = value; row.ref = ref;
        row.relerr = math::abs(value-ref)/math::abs(ref); row.kind = GateRow::Report;
        push(row);
    };
    reportRow("volume", P.V, V_exact);
    reportRow("area",   P.A, A_exact);
    reportRow("flux",   P.F, F_exact);

    GateRow mv; mv.check = "minweight_vol"; mv.value = P.minWVol; mv.ref = 0; mv.relerr = 0;
    mv.kind = GateRow::Report; push(mv);
    GateRow mb; mb.check = "minweight_bdr"; mb.value = P.minWBdr; mb.ref = 0; mb.relerr = 0;
    mb.kind = GateRow::Report; push(mb);

    if (P.pulledBack)
    {
        GateRow dj; dj.check = "detJ_min"; dj.value = P.minDetJ; dj.ref = 0; dj.relerr = 0; dj.tol = 0;
        dj.kind = GateRow::Required; dj.positive = true; dj.pass = (P.minDetJ > 0);
        push(dj);

        GateRow nd; nd.check = "normal_dev"; nd.value = P.maxNormalDev; nd.ref = 0; nd.tol = 1e-12;
        nd.kind = GateRow::Required; nd.relerr = P.maxNormalDev; nd.pass = (P.maxNormalDev <= 1e-12);
        push(nd);
    }

    GateRow run; run.check = "run"; run.value = (real_t)P.errorCells.size(); run.ref = 0;
    run.tol = 0; run.kind = GateRow::Required; run.relerr = run.value;
    run.failCells = P.errorCells; run.pass = P.errorCells.empty();
    push(run);

    return G;
}

//----------------------------------------------------------------------------
// Printing (--study volume only).
//----------------------------------------------------------------------------

/// Prints one `GATE` row. \a deferred (tchakaloff only, see \ref runVolumeStudy
/// for why its Required rows are not gated on) leaves \a row's value/ref/
/// relerr/tol/failCells computation untouched -- it only replaces the
/// PASS/FAIL status token of an otherwise-Required row with `DEFERRED`
/// and prints the underlying verdict as `wouldPass=0|1`, so the row stays
/// fully auditable without being counted into `GATE SUMMARY`. Report rows
/// are unaffected by \a deferred.
void printGateRow(const GateResult & G, const GateRow & row, bool deferred = false)
{
    gsInfo << "GATE case=" << G.caseName << " mode=" << G.mode;
    if (!G.geo.empty()) gsInfo << " geo=" << G.geo;
    gsInfo << " r=" << G.r << " n=" << G.n
          << " check=" << row.check << " value=" << gsTetClip::fmtSci(row.value)
          << " ref=" << gsTetClip::fmtSci(row.ref) << " relerr=" << gsTetClip::fmtSci(row.relerr)
          << " tol=" << (GateRow::Report == row.kind ? "report" :
                        (row.positive ? "positive" : gsTetClip::fmtSci(row.tol)));
    if (!row.failCells.empty())
    {
        gsInfo << " failCells=" << row.failCells.size() << ":";
        const size_t nShow = std::min((size_t)20, row.failCells.size());
        for (size_t c = 0; c != nShow; ++c)
        {
            if (0 != c) gsInfo << ",";
            gsInfo << row.failCells[c];
        }
    }
    const bool isDeferredRow = deferred && (GateRow::Required == row.kind);
    if (isDeferredRow) gsInfo << " wouldPass=" << (row.pass ? 1 : 0);
    gsInfo << " " << (GateRow::Report == row.kind ? "REPORT" :
                      (isDeferredRow ? "DEFERRED" : (row.pass ? "PASS" : "FAIL"))) << "\n";
}

/// No-relax policy for a failing moments_*/residual_* row: one INFO line
/// per failing cell (id, nIn, nOut, momErr, rank, levels), so a BLOCKED
/// report can quote the exact numbers without re-running anything.
void printMomentResidualFailures(const GateResult & G, const CellRulePass & P, const GateRow & row)
{
    if (row.failCells.empty()) return;
    if ("moments_vol" != row.check && "moments_bdr" != row.check &&
        "residual_vol" != row.check && "residual_bdr" != row.check) return;
    const bool isVol = (row.check.find("vol") != std::string::npos);
    for (size_t id : row.failCells)
    {
        const CellLog & L = P.log[id];
        gsInfo << "INFO case=" << G.caseName << " mode=" << G.mode << " r=" << G.r
              << " check=" << row.check << " id=" << id
              << " nIn="  << (isVol ? L.nVolIn  : L.nBdrIn)
              << " nOut=" << (isVol ? L.nVolOut : L.nBdrOut)
              << " momErr=" << gsTetClip::fmtSci(isVol ? L.momErrVol : L.momErrBdr)
              << " rank="   << (isVol ? L.rankVol   : L.rankBdr)
              << " levels=" << (isVol ? L.levelsVol : L.levelsBdr) << "  INFO\n";
    }
}

/// One INFO line per failing `fluxmom_cell` cell (every one, not just the
/// first 20 the `GATE` row itself prints): id and E_K
/// (`CellLog::fluxErrBdr`), so the per-cell distribution of the flux-moment
/// error can be ranked from a single run.
void printFluxmomCellFailures(const GateResult & G, const CellRulePass & P, const GateRow & row)
{
    if ("fluxmom_cell" != row.check) return;
    for (size_t id : row.failCells)
        gsInfo << "INFO case=" << G.caseName << " mode=" << G.mode << " r=" << G.r
              << " check=fluxmom_cell id=" << id
              << " E_K=" << gsTetClip::fmtSci(P.log[id].fluxErrBdr) << "  INFO\n";
}

/// One INFO line per `run`-failure cell (every one, not just the first 5
/// GATEERR lines below): the input point counts CellLog already recorded
/// before the exception, so a cell whose own compressor call threw still
/// leaves a diagnosable trace of how large its input was.
void printRunFailureInfo(const GateResult & G, const CellRulePass & P, const GateRow & row)
{
    if ("run" != row.check) return;
    for (size_t id : row.failCells)
    {
        const CellLog & L = P.log[id];
        gsInfo << "INFO case=" << G.caseName << " mode=" << G.mode << " r=" << G.r
              << " check=run id=" << id << " nVolIn=" << L.nVolIn << " nBdrIn=" << L.nBdrIn
              << "  INFO\n";
    }
}

/// Up to 5 GATEERR lines (via gsWarn) for a run failure's failing cells.
void printRunErrors(const std::string & caseName, const std::string & mode, index_t r,
                    const std::vector<size_t> & failCells, const std::vector<std::string> & msgs)
{
    const size_t nShow = std::min((size_t)5, failCells.size());
    for (size_t s = 0; s != nShow; ++s)
        gsWarn << "GATEERR case=" << caseName << " mode=" << mode << " r=" << r
              << " cell=" << failCells[s] << " what=" << firstLine(msgs[s]) << "\n";
}

/// GATEDIAG row for a tet mode: totals in 64-bit (volIn/volOut/bdrIn/bdrOut
/// reach O(10^8) at r=3, and index_t may be 32-bit). Fields that do not
/// apply to \a G.mode (rank/twoLevel: tchakaloff only; passThrough: momrule
/// only) print "na".
void printGateDiagTet(const GateResult & G, const CellRulePass & P, const gsTetClip::ClipStreamer & S,
                      real_t gateEvalSeconds)
{
    const gsTetClip::CellIndex & idx = *S.index();
    long long volIn=0, volOut=0, bdrIn=0, bdrOut=0, negWVol=0, negWBdr=0;
    index_t rankVolMin = std::numeric_limits<index_t>::max(), rankVolMax = -1;
    index_t rankBdrMin = std::numeric_limits<index_t>::max(), rankBdrMax = -1;
    index_t twoLevelVol=0, twoLevelBdr=0, passThroughVol=0, passThroughBdr=0;
    gsTetClip::KahanSum compressCpu;

    for (size_t id = 0; id != P.log.size(); ++id)
    {
        const CellLog & L = P.log[id];
        compressCpu.add(L.compressSec);
        if (gsTetClip::Cut == idx.status[id])
        {
            volIn += L.nVolIn; volOut += L.nVolOut; negWVol += L.negWVol;
            if ("tchakaloff" == G.mode && L.error.empty())
            {
                rankVolMin = math::min(rankVolMin, L.rankVol);
                rankVolMax = math::max(rankVolMax, L.rankVol);
                if (2 == L.levelsVol) ++twoLevelVol;
            }
            if ("momrule" == G.mode && L.passThroughVol) ++passThroughVol;
        }
    }
    for (size_t w = 0; w != P.work.size(); ++w)
    {
        const CellLog & L = P.log[P.work[w]];
        bdrIn += L.nBdrIn; bdrOut += L.nBdrOut; negWBdr += L.negWBdr;
        if (L.nBdrIn > 0)
        {
            if ("tchakaloff" == G.mode && L.error.empty())
            {
                rankBdrMin = math::min(rankBdrMin, L.rankBdr);
                rankBdrMax = math::max(rankBdrMax, L.rankBdr);
                if (2 == L.levelsBdr) ++twoLevelBdr;
            }
            if ("momrule" == G.mode && L.passThroughBdr) ++passThroughBdr;
        }
    }

    auto ratioOf = [](long long in, long long out) -> std::string
    { return (0 == out) ? "na" : gsTetClip::fmtSci((real_t)in/(real_t)out); };
    auto rankStr = [&](index_t lo, index_t hi) -> std::string
    {
        if ("tchakaloff" != G.mode || hi < 0) return "na/na";
        std::ostringstream oss; oss << lo << "/" << hi; return oss.str();
    };
    auto countOrNa = [&](const std::string & wantMode, index_t v) -> std::string
    {
        if (G.mode != wantMode) return "na";
        std::ostringstream oss; oss << v; return oss.str();
    };

    gsInfo << "GATEDIAG case=" << G.caseName << " mode=" << G.mode << " r=" << G.r << " n=" << G.n
          << " cut=" << S.numCut() << " full=" << S.numFull() << " empty=" << S.numEmpty()
          << " work=" << P.work.size()
          << " volIn=" << volIn << " volOut=" << volOut << " bdrIn=" << bdrIn << " bdrOut=" << bdrOut
          << " ratioVol=" << ratioOf(volIn,volOut) << " ratioBdr=" << ratioOf(bdrIn,bdrOut)
          << " negWVol=" << negWVol << " negWBdr=" << negWBdr
          << " rankVol=" << rankStr(rankVolMin, rankVolMax) << " rankBdr=" << rankStr(rankBdrMin, rankBdrMax)
          << " twoLevelVol=" << countOrNa("tchakaloff", twoLevelVol)
          << " twoLevelBdr=" << countOrNa("tchakaloff", twoLevelBdr)
          << " passThroughVol=" << countOrNa("momrule", passThroughVol)
          << " passThroughBdr=" << countOrNa("momrule", passThroughBdr)
          << " compressCpu=" << gsTetClip::fmtSci(compressCpu.value()) << "s"
          << " passTime=" << gsTetClip::fmtSci(P.seconds) << "s"
          << " gateTime=" << gsTetClip::fmtSci(P.seconds + gateEvalSeconds) << "s  INFO\n";
}

/// `GATEDIAG` row for the Algoim mode: point counts, negative-weight
/// counts and the raw `gsAlgoimAdaptiveRule::Stats` accumulated across
/// every background cell of \a P.
void printGateDiagAlgoim(const AlgoimPass & P, const std::string & caseName, index_t r, index_t n)
{
    gsInfo << "GATEDIAG case=" << caseName << " mode=algoim r=" << r << " n=" << n
          << " cells=" << (long long)n*(long long)n*(long long)n
          << " volPts=" << P.volPts << " bdrPts=" << P.bdrPts
          << " negWVol=" << P.negWVol << " negWBdr=" << P.negWBdr
          << " nSubBoxes=" << P.volStats.nSubBoxes << " nFallbackLeaves=" << P.volStats.nFallbackLeaves
          << " nZeroSurfaceLeaves=" << P.srfStats.nZeroSurfaceLeaves
          << " unresolvedSurfaceMeasure=" << gsTetClip::fmtSci(P.srfStats.unresolvedSurfaceMeasure)
          << " passTime=" << gsTetClip::fmtSci(P.seconds) << "s  INFO\n";
}

//----------------------------------------------------------------------------
// --study volume driver.
//----------------------------------------------------------------------------

/// Runs the volume/area/flux/moment gate over `--case`/`--mode`/`-r`. Builds
/// each (case, r)'s ClipStreamer once and gates clip/tchakaloff/momrule in
/// turn on it (each CellRulePass goes out of scope before the next mode
/// starts, so at most one compressed table set is alive at a time), then
/// the sphere-only Algoim baseline. Returns whether every counted Required
/// row PASSed over the whole run; `main` maps this to the process exit
/// code. Tchakaloff's would-be-Required rows print `wouldPass=0|1 DEFERRED`
/// (see \ref printGateRow) and are excluded from that count and from
/// `GATE SUMMARY`, because its NNLS compressor's KKT stopping test does
/// not bound the primal residual (the TODO on `nnlsLawsonHansonImpl` in
/// `gsTchakaloffRule.h`): a tchakaloff-only run therefore always reports
/// `required=0` and exits successfully regardless of `wouldPass`.
bool runVolumeStudy(const Config & cfg)
{
    const gsTetClip::BgMapKind kind = geoKind(cfg.geo);

    std::vector<std::string> cases;
    if ("all" == cfg.caseName) { cases.push_back("sphere"); cases.push_back("rotcube"); }
    else cases.push_back(cfg.caseName);

    std::vector<std::string> modes;
    if ("all" == cfg.mode)
    { modes = { "clip", "tchakaloff", "momrule", "algoim" }; }
    else modes.push_back(cfg.mode);

    const bool wantAlgoim  = std::find(modes.begin(), modes.end(), "algoim")    != modes.end();
    const bool wantMomrule = std::find(modes.begin(), modes.end(), "momrule")   != modes.end();
    std::vector<std::string> tetModes;
    for (const std::string & m : modes) if ("algoim" != m) tetModes.push_back(m);

    index_t reqTotal = 0, passTotal = 0;

    for (const std::string & caseName : cases)
    {
        const gsTetClip::PhysTetMesh phys = loadPhysMesh(cfg, caseName);

        const real_t V_mesh = gsTetClip::meshVolumeExact(phys);
        const real_t A_mesh = gsTetClip::unclippedBoundaryArea(phys);
        const real_t V_exact = ("sphere" == caseName) ? sphereVolumeExact()  : rotcubeVolumeExact();
        const real_t A_exact = ("sphere" == caseName) ? sphereAreaExact()    : rotcubeAreaExact();

        gsInfo << "GATEINFO case=" << caseName << " V_mesh=" << gsTetClip::fmtSci(V_mesh)
              << " A_mesh=" << gsTetClip::fmtSci(A_mesh) << " V_exact=" << gsTetClip::fmtSci(V_exact)
              << " A_exact=" << gsTetClip::fmtSci(A_exact) << "  INFO\n";

        if ("rotcube" == caseName && wantAlgoim)
            gsInfo << "GATEINFO case=rotcube mode=algoim skipped: no analytic level set  INFO\n";

        memory::shared_ptr<gsMeshSignedDist<real_t> > phiH;
        if (wantMomrule && !tetModes.empty())
        {
            gsSurfMesh<real_t> surf;
            std::vector<gsSurfMesh<real_t>::Vertex> vmap(phys.m.P.size());
            for (const std::array<index_t,4> & tri : phys.m.bdrTri)
            {
                for (int v = 0; v != 3; ++v)
                {
                    const index_t vi = tri[v];
                    if (!vmap[vi].is_valid())
                    {
                        gsSurfMesh<real_t>::Point pt; pt << phys.m.P[vi][0], phys.m.P[vi][1], phys.m.P[vi][2];
                        vmap[vi] = surf.add_vertex(pt);
                    }
                }
                const gsSurfMesh<real_t>::Face f = surf.add_triangle(vmap[tri[0]], vmap[tri[1]], vmap[tri[2]]);
                GISMO_ENSURE(f.is_valid(), "runVolumeStudy: add_triangle rejected a boundary face "
                            "(would create a complex edge).");
            }
            gsMatrix<real_t> bbox(3,2);
            bbox(0,0) = -1.0; bbox(0,1) = 1.0;
            bbox(1,0) = -1.0; bbox(1,1) = 1.0;
            bbox(2,0) = -1.0; bbox(2,1) = 1.0;
            phiH = memory::make_shared(new gsMeshSignedDist<real_t>(surf, bbox));
        }

        real_t pbScale = 1.0;
        if (gsTetClip::BgMapKind::Identity != kind)
            for (const gsTetClip::Vec3 & x : phys.m.P)
                pbScale = math::max(pbScale, math::max(std::abs(x[0]),
                          math::max(std::abs(x[1]), std::abs(x[2]))));

        for (index_t r = 0; r <= cfg.rMax; ++r)
        {
            const gsTetClip::Grid3 grid = makeBoxGrid(cfg.n0, r);
            const index_t n = grid.n;

            // Pull-back (gsImmersedPullback.h): runs BEFORE any ClipStreamer
            // or assembly, even when tetModes is empty (--mode algoim on a
            // curved background is rejected by main(), but the grid/mesh
            // ENSURE below must still run for identity --mode algoim runs).
            gsMultiPatch<real_t> mp; gsMultiBasis<real_t> mb;
            makeBackground(grid, cfg.p, kind, cfg.eps, mp, mb);
            const gsGeometry<real_t> & bgMap = mp.patch(0);

            gsTetClip::PullbackStats pbStats;
            const gsTetClip::ParamTetMesh param =
                gsTetClip::pullBack(phys, bgMap, kind, grid, pbStats);

            if (gsTetClip::BgMapKind::Identity != kind)
            {
                const real_t tol = 1e-12*pbScale;
                const bool pbPass = (0 == pbStats.nFailed && 0 == pbStats.nOutside &&
                                    0 == pbStats.nInverted && pbStats.maxResidual <= tol);
                gsInfo << "PULLBACK case=" << caseName << " geo=" << cfg.geo << " r=" << r
                      << " n=" << n << " nVerts=" << pbStats.nVerts << " nFailed=" << pbStats.nFailed
                      << " nOutside=" << pbStats.nOutside << " nInverted=" << pbStats.nInverted
                      << " maxResidual=" << gsTetClip::fmtSci(pbStats.maxResidual)
                      << " tol=" << gsTetClip::fmtSci(tol)
                      << " minParamTetVol=" << gsTetClip::fmtSci(pbStats.minParamTetVol)
                      << " " << (pbPass ? "PASS" : "FAIL") << "\n";
                if (!pbPass)
                {
                    gsWarn << "PULLBACK-ABORT case=" << caseName << " geo=" << cfg.geo << " r=" << r
                          << ": pull-back failed its acceptance check\n";
                    return false;
                }
            }

            makeGrid(param.m, cfg.n0, r, "runVolumeStudy");   // strict-inside ENSURE, parametric

            const memory::shared_ptr<gsTetClip::TetMesh> Mclip =
                memory::make_shared(new gsTetClip::TetMesh(param.m));

            if (!tetModes.empty())
            {
                gsStopwatch sswatch;
                const memory::shared_ptr<gsTetClip::ClipStreamer> S =
                    memory::make_shared(new gsTetClip::ClipStreamer(Mclip, grid, cfg.p));
                const real_t streamerTime = sswatch.stop();

                gsInfo << "GATEINFO case=" << caseName << " r=" << r << " n=" << n
                      << " streamerTime=" << gsTetClip::fmtSci(streamerTime) << "s"
                      << " cut=" << S->numCut() << " full=" << S->numFull()
                      << " empty=" << S->numEmpty() << "  INFO\n";

                if (gsTetClip::BgMapKind::Identity == kind)
                {
                    for (const std::string & mode : tetModes)
                    {
                        const CellRulePass P = buildCellRules(*S, mode, cfg.p, phiH.get(), true);

                        gsStopwatch gsw;
                        const GateResult G =
                            gateTetMode(P, *S, cfg.p, V_mesh, A_mesh, V_exact, A_exact, caseName, r);
                        const real_t gateEvalSeconds = gsw.stop();

                        // The NNLS compressor's KKT stopping test does not bound the
                        // primal residual (gsTchakaloffRule.h, nnlsLawsonHansonImpl),
                        // so tchakaloff's rows are not treated as pass/fail gates:
                        // they print DEFERRED and never enter reqTotal/passTotal.
                        const bool deferred = ("tchakaloff" == mode);

                        std::vector<size_t> runFail; std::vector<std::string> runMsgs;
                        for (const GateRow & row : G.rows)
                        {
                            printGateRow(G, row, deferred);
                            printMomentResidualFailures(G, P, row);
                            printFluxmomCellFailures(G, P, row);
                            printRunFailureInfo(G, P, row);
                            if (!deferred)
                            {
                                reqTotal  += (GateRow::Required == row.kind) ? 1 : 0;
                                passTotal += (GateRow::Required == row.kind && row.pass) ? 1 : 0;
                            }
                            if ("run" == row.check)
                                for (size_t id : row.failCells) { runFail.push_back(id); runMsgs.push_back(P.log[id].error); }
                        }
                        printRunErrors(caseName, mode, r, runFail, runMsgs);
                        printGateDiagTet(G, P, *S, gateEvalSeconds);
                        if (deferred)
                            gsInfo << "GATE-DEFERRED mode=tchakaloff reason=nnls-kkt-stall\n";
                    }
                }
                else
                {
                    // Curved background (clip only; main() rejects the other
                    // modes): the physical-space gate recomputes volume/
                    // area/flux/fluxmom_global directly from S and bgMap
                    // (gatePhysical), since CellLog's totals are parametric
                    // here. bdrSrc is the Nanson-mapped boundary source of
                    // this (case, r).
                    const gsTetClip::PullbackBdrSource bdrSrc(*S, bgMap);
                    const GateResult G = gatePhysical(phys, pbStats, bgMap, kind, *S, bdrSrc,
                                                      cfg.p, V_exact, A_exact, caseName, cfg.geo, r);
                    for (const GateRow & row : G.rows)
                    {
                        printGateRow(G, row);
                        reqTotal  += (GateRow::Required == row.kind) ? 1 : 0;
                        passTotal += (GateRow::Required == row.kind && row.pass) ? 1 : 0;
                    }
                }
            }

            if (wantAlgoim && "sphere" == caseName)
            {
                gsFunctionExpr<real_t> phi("(x-0.03)^2+(y+0.02)^2+(z-0.01)^2-0.3025", 3);
                const bool pulledV = !gsTetClip::isIdentityMap(bgMap);
                PulledBackSphere phiHatV(bgMap);
                const gsFunction<real_t> & phiUseV = pulledV
                    ? static_cast<const gsFunction<real_t> &>(phiHatV)
                    : static_cast<const gsFunction<real_t> &>(phi);
                const AlgoimPass AP = runAlgoimPass(grid, cfg.p, phiUseV, nullptr, nullptr,
                                                    pulledV ? &bgMap : nullptr);
                const GateResult G = gateAlgoim(AP, cfg.p, caseName, r, n);
                for (const GateRow & row : G.rows)
                {
                    printGateRow(G, row);
                    reqTotal  += (GateRow::Required == row.kind) ? 1 : 0;
                    passTotal += (GateRow::Required == row.kind && row.pass) ? 1 : 0;
                }
                printRunErrors(caseName, "algoim", r, AP.errorCells, AP.errors);
                printGateDiagAlgoim(AP, caseName, r, n);
            }
        }
    }

    gsInfo << "GATE SUMMARY required=" << reqTotal << " pass=" << passTotal
          << " fail=" << (reqTotal-passTotal) << "\n";
    return (reqTotal == passTotal);
}

//----------------------------------------------------------------------------
// --study poisson: 3D immersed Poisson (symmetric Nitsche + ghost penalty).
// See the file header for the weak form, the manufactured solution, the
// ghost-penalty scaling and why only jump order k = p is assembled.
//----------------------------------------------------------------------------

/// Quadrature-cost/positivity diagnostics printed on every `POISSON` row:
/// \a nqVol is the SOLVE rule's total volume node count (served Cut-cell
/// nodes summed over every Cut cell, plus #Full cells * (p+1)^3 -- exactly
/// what \ref gsTetClip::VolLookupRule serves); \a nqBdr is its total
/// boundary node count; \a minw is the minimum served VOLUME weight over
/// Cut cells (+inf if there are none).
struct PoissonStats
{
    long long nqVol = 0, nqBdr = 0;
    real_t minw = std::numeric_limits<real_t>::infinity();
};

/// \ref PoissonStats of a tet-mode pass \a P (clip/tchakaloff/momrule) over
/// streamer \a S: \a P.log already carries the served per-cell node counts
/// and minimum weight from \ref buildCellRules's own loop, for every mode
/// including clip, so this is a serial O(n^3) read of already-materialized
/// data -- no re-streaming, no re-compression.
PoissonStats statsFromTetPass(const CellRulePass & P, const gsTetClip::ClipStreamer & S, index_t p)
{
    const gsTetClip::CellIndex & idx = *S.index();
    PoissonStats st;
    long long nv = 0, nb = 0;
    for (size_t id = 0; id != P.log.size(); ++id)
        if (gsTetClip::Cut == idx.status[id])
        {
            nv += P.log[id].nVolOut;
            st.minw = math::min(st.minw, P.log[id].minWVol);
        }
    nv += (long long)S.numFull() * (long long)(p+1)*(long long)(p+1)*(long long)(p+1);
    for (size_t w = 0; w != P.work.size(); ++w) nb += P.log[P.work[w]].nBdrOut;
    st.nqVol = nv; st.nqBdr = nb;
    return st;
}

/// Builds the boundary-mesh signed distance \a momrule needs for its served
/// normals. Independent of \ref runVolumeStudy's own inline construction:
/// the two build sites live under different lifetime scopes (one
/// gsMeshSignedDist per (case, r) there; one per (case, mode) here).
memory::shared_ptr<gsMeshSignedDist<real_t> > buildMeshSignedDist(const gsTetClip::TetMesh & M)
{
    gsSurfMesh<real_t> surf;
    std::vector<gsSurfMesh<real_t>::Vertex> vmap(M.P.size());
    for (const std::array<index_t,4> & tri : M.bdrTri)
    {
        for (int v = 0; v != 3; ++v)
        {
            const index_t vi = tri[v];
            if (!vmap[vi].is_valid())
            {
                gsSurfMesh<real_t>::Point pt; pt << M.P[vi][0], M.P[vi][1], M.P[vi][2];
                vmap[vi] = surf.add_vertex(pt);
            }
        }
        const gsSurfMesh<real_t>::Face f = surf.add_triangle(vmap[tri[0]], vmap[tri[1]], vmap[tri[2]]);
        GISMO_ENSURE(f.is_valid(), "buildMeshSignedDist: add_triangle rejected a boundary face "
                    "(would create a complex edge).");
    }
    gsMatrix<real_t> bbox(3,2);
    bbox(0,0) = -1.0; bbox(0,1) = 1.0;
    bbox(1,0) = -1.0; bbox(1,1) = 1.0;
    bbox(2,0) = -1.0; bbox(2,1) = 1.0;
    return memory::make_shared(new gsMeshSignedDist<real_t>(surf, bbox));
}

/// Resets the process's peak-RSS high-water mark (Linux >= 4.0's
/// /proc/self/clear_refs, mode 5), so a later \ref peakRssKB reads the peak
/// of only the work since this call. A silent no-op if the pseudo-file is
/// unwritable (a sandboxed or non-Linux run).
void resetPeakRss()
{
    std::ofstream f("/proc/self/clear_refs");
    if (f.is_open()) f << "5\n";
}

/// Peak resident set size (VmHWM, /proc/self/status) since the last \ref
/// resetPeakRss call, in KiB; -1 if unavailable.
long long peakRssKB()
{
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
        if (0 == line.compare(0, 6, "VmHWM:"))
        {
            std::istringstream iss(line.substr(6));
            long long kb = -1; iss >> kb; return kb;
        }
    return -1;
}

/// Plain-scalar outcome of one (case, mode, r) assemble-and-solve, Steps
/// 4-9 of the file header's pipeline. \a t_asm covers dof elimination
/// through the end of assembly; \a t_solve the LU factorization and the
/// primal solve. \a t_err the two reference-rule error integrals.
/// \a zeroRows > 0 or !\a finite means Steps 8-9 were skipped (errors/
/// conditioning are left at their defaults); \a kappaMethod is empty in
/// that case too.
struct PoissonRunResult
{
    index_t ndof = 0, zeroRows = 0;
    bool finite = false;
    real_t L2  = std::numeric_limits<real_t>::quiet_NaN();
    real_t H1s = std::numeric_limits<real_t>::quiet_NaN();
    real_t kappa = 0; bool kappaIndef = true;
    bool kappaConverged = false;                 // dense: always true; iter: both loops converged
    std::string kappaMethod;                    // "dense" | "iter" | "" (not computed)
    long long kItPower = -1, kItInverse = -1;
    real_t lmin = 0, lmax = 0, asym = 0;
    real_t t_asm = 0, t_solve = 0, t_err = 0;
};

/// Assembles and solves the symmetric-Nitsche + ghost-penalty immersed
/// Poisson problem on the background space (\a mp, \a mb) restricted to
/// \a dom (a TetClipSignDomain for the tet modes, a
/// gsImplicitTrimmedDomain<3,real_t> for algoim -- both a
/// gsTrimmedDomain<3,real_t>), using \a idx / \a volSrc / \a bdrSrc for the
/// SOLVE rule (gsTetClip::VolLookupRule/BdrLookupRule/BdrNormalField,
/// nGaussFull = (p+1)^3, exact on Q_{2p+1}) and \a refFactory for the two
/// reference-rule error integrals (the mode's own, finer, reference
/// geometry). See the file header for the weak form and the ghost-penalty
/// scaling. Conditioning uses a dense SelfAdjointEigenSolver at N <=
/// \a kappaDense (always converged), else power/inverse iteration, each
/// capped at \a kappaMaxIt steps, from the deterministic start
/// v_i = 1 + 1e-3*(i mod 7), reusing the primal solve's own LU
/// factorization for the inverse iteration. \ref PoissonRunResult::kappaConverged
/// is true iff BOTH loops exited through their own convergence test rather
/// than the cap; when false, at least one loop hit the cap, so lmax (a
/// Rayleigh quotient, <= lambda_max) and/or lmin (>= lambda_min for a
/// symmetric positive-definite K) may be far from the extreme eigenvalue,
/// and kappa = lmax/lmin is then an underestimate.
///
/// \a mpId is a SECOND background map -- always exactly the identity on the
/// same grid (gsTetClip::identityBoxGeometry), even when \a mp is too --
/// through which the ghost penalty's dnk is assembled by default
/// (\a ghostUsesG false); \a ghostUsesG true routes it through \a mp
/// instead (debug, `--ghostMap G`). \a hNitsche is the Nitsche penalty's
/// own h (identity: h itself; curved: h scaled by a sampled min singular
/// value of \a mp's Jacobian, computed by the caller -- see the file
/// header). The ghost term itself always keeps the parametric \a h. When
/// \a dumpFile is non-empty, the assembled system (CSC arrays + rhs) is
/// written to it right after assembly, before the solve.
PoissonRunResult solvePoissonOnDomain(
    const gsMultiPatch<real_t> & mp, gsMultiBasis<real_t> & mb, const gsMultiPatch<real_t> & mpId,
    memory::shared_ptr<gsTrimmedDomain<3,real_t> > dom,
    memory::shared_ptr<const gsTetClip::CellIndex> idx,
    memory::shared_ptr<const gsTetClip::VolCellSource> volSrc,
    memory::shared_ptr<const gsTetClip::BdrCellSource> bdrSrc,
    gsExprAssembler<real_t>::QuadratureFactory refFactory,
    index_t p, real_t h, real_t hNitsche, real_t gammaEff, real_t gtEff, bool ghostOn,
    bool ghostUsesG, index_t kappaDense, index_t kappaMaxIt,
    const gsFunctionExpr<real_t> & u_exact, const gsFunctionExpr<real_t> & f_rhs,
    const std::string & dumpFile)
{
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::space       space;
    typedef gsExprAssembler<real_t>::solution    solution;

    GISMO_ENSURE(&mpId != &mp && gsTetClip::isIdentityMap(mpId.patch(0)),
                "solvePoissonOnDomain: mpId must be a distinct, exactly identity background map.");

    PoissonRunResult R;

    gsExprAssembler<real_t> A(1,1);
    geometryMap G = A.getMap(mp);
    space u = A.getSpace(mb);
    auto ff  = A.getCoeff(f_rhs, G);
    auto g_D = A.getCoeff(u_exact, G);
    gsTetClip::BdrNormalField nf(idx, bdrSrc);
    auto nfC = A.getCoeff(nf);                          // NO geometry: nf must see the rule's own nodes
    std::vector<patchSide> bdrImm(1, patchSide(0, boundary::none));
    const geometryMap Gg = ghostUsesG ? G : A.getMap(mpId);
    auto dJ = dnk(u.jump(), Gg.left(), p);

    gsStopwatch swAsm;

    // Exterior-dof elimination: a dof is kept iff it is active on the
    // midpoint of some Interior or Cut element. Equivalent to
    // gsDofMapperCreator's createMapper on its non-conforming, BC-free path
    // (gsDofMapperCreator.hpp), built directly through the one public
    // gsDofMapper constructor, gsDofMapper(gsVector<index_t> patchDofSizes, nComp).
    {
        gsVector<index_t> patchDofSizes(mb.nPieces());
        for (index_t k = 0; k != patchDofSizes.size(); ++k)
            patchDofSizes[k] = mb.basis(k).size();
        gsDofMapper mapper(patchDofSizes, 1);
        std::vector<bool> keep(mb.basis(0).size(), false);
        gsMatrix<real_t> centre(3,1);
        gsMatrix<index_t> act;

        gsDomain<real_t>::iterator eInt = dom->end<InteriorSign>();
        for (gsDomain<real_t>::iterator it = dom->beginInterior(); it < eInt; ++it)
        {
            centre.col(0) = 0.5*(it.lowerCorner()+it.upperCorner());
            mb.basis(0).active_into(centre, act);
            for (index_t i = 0; i != act.rows(); ++i) keep[act(i,0)] = true;
        }
        gsDomain<real_t>::iterator eCut = dom->endBdr(boundary::none);
        for (gsDomain<real_t>::iterator it = dom->beginBdr(boundary::none); it < eCut; ++it)
        {
            centre.col(0) = 0.5*(it.lowerCorner()+it.upperCorner());
            mb.basis(0).active_into(centre, act);
            for (index_t i = 0; i != act.rows(); ++i) keep[act(i,0)] = true;
        }
        for (index_t i = 0; i != (index_t)keep.size(); ++i)
            if (!keep[i]) mapper.eliminateDof(i, 0);
        mapper.finalize();
        u.setupMapper(mapper);
        const_cast<expr::gsFeSpace<real_t>&>(u).fixedPart().setZero(mapper.boundarySize(), 1);
    }

    // Assembly. Order and scoping are mandatory: an installed factory is
    // used for EVERY parametric direction of a ghost/skeleton face, and
    // assembleBdr's fixedDirection for the immersed boundary (0) collides
    // with an x-direction ghost face, so a volume/boundary factory must
    // never be installed while computePattern[Ghost]/assembleGhost run.
    A.setIntegrationElements(mb);
    A.initSystem();
    A.computePattern(igrad(u) * igrad(u).tr());                    // full mesh, no factory
    A.setIntegrationDomain(dom);
    if (ghostOn) A.computePatternGhost(dJ * dJ.tr());               // no factory
    {
        gsTetClip::VolumeQuadratureScope<gsExprAssembler<real_t> > s(
            A, gsTetClip::makeVolLookupFactory(idx, volSrc, gsVector<index_t>::Constant(3, p+1)));
        A.assemble(igrad(u,G) * igrad(u,G).tr() * meas(G),  u * ff * meas(G));
    }
    {
        gsTetClip::BoundaryQuadratureScope<gsExprAssembler<real_t> > s(
            A, gsTetClip::makeBdrLookupFactory(idx, bdrSrc));
        A.assembleBdr(bdrImm,
            - (igrad(u,G) * nfC) * u.tr()
            - u * (igrad(u,G) * nfC).tr()
            + (gammaEff / hNitsche) * u * u.tr() );
        A.assembleBdr(bdrImm,
            - (igrad(u,G) * nfC) * g_D
            + (gammaEff / hNitsche) * u * g_D );
    }
    GISMO_ENSURE(!A.hasCustomQuadrature(), "solvePoissonOnDomain: a quadrature factory leaked "
                "past its scope.");
    if (ghostOn)
        A.assembleGhost(gtEff * math::pow(h, (real_t)(2*p-1)) * dJ * dJ.tr());   // LAST, unscoped

    R.t_asm = swAsm.stop();
    R.ndof = A.numDofs();

    // Raw system dump (debug, --dumpSystem): CSC arrays + rhs, right after
    // assembly, before the solve -- so two runs that differ only in
    // --ghostMap can be compared byte-for-byte.
    if (!dumpFile.empty())
    {
        gsSparseMatrix<real_t> K = A.matrix();
        K.makeCompressed();
        std::ofstream ofs(dumpFile.c_str(), std::ios::binary);
        GISMO_ENSURE(ofs, "solvePoissonOnDomain: failed to open dump file '" << dumpFile << "'");
        const long long dumpRows = K.rows(), dumpCols = K.cols(), dumpNnz = K.nonZeros();
        ofs.write(reinterpret_cast<const char*>(&dumpRows), sizeof(dumpRows));
        ofs.write(reinterpret_cast<const char*>(&dumpCols), sizeof(dumpCols));
        ofs.write(reinterpret_cast<const char*>(&dumpNnz),  sizeof(dumpNnz));
        ofs.write(reinterpret_cast<const char*>(K.outerIndexPtr()),
                  (size_t)(K.outerSize()+1)*sizeof(index_t));
        ofs.write(reinterpret_cast<const char*>(K.innerIndexPtr()), (size_t)dumpNnz*sizeof(index_t));
        ofs.write(reinterpret_cast<const char*>(K.valuePtr()),      (size_t)dumpNnz*sizeof(real_t));
        const gsMatrix<real_t> & rhs = A.rhs();
        ofs.write(reinterpret_cast<const char*>(rhs.data()), (size_t)rhs.rows()*sizeof(real_t));
        gsInfo << "SYSTEM-DUMP file=" << dumpFile << " rows=" << dumpRows << " nnz=" << dumpNnz << "\n";
    }

    // Zero-row check -- counted, not ENSUREd: a failing row is a
    // reportable outcome here, not a programming error.
    {
        const gsSparseMatrix<real_t> & M = A.matrix();
        gsVector<real_t> rowSq = gsVector<real_t>::Zero(M.rows());
        for (index_t c = 0; c != M.outerSize(); ++c)
            for (gsSparseMatrix<real_t>::iterator it(M, c); it; ++it)
                rowSq[it.row()] += it.value()*it.value();
        for (index_t i = 0; i != rowSq.rows(); ++i)
            if (0.0 == rowSq[i]) ++R.zeroRows;
    }
    if (R.zeroRows > 0) return R;

    // Direct LU solve (gsEigen::SparseLU/COLAMD): a non-finite result is
    // diagnosable (succeed()/allFinite()), where an iterative solver
    // hitting its cap would instead return a capped, garbage iterate; the
    // same factorization also serves the inverse iteration below.
    gsStopwatch swSolve;
    gsSparseSolver<real_t>::LU solver;
    solver.compute(A.matrix());
    gsMatrix<real_t> solVector;                 // NAMED: getSolution keeps a pointer to it
    solVector = solver.solve(A.rhs());
    R.finite = solver.succeed() && solVector.allFinite();
    R.t_solve = swSolve.stop();
    if (!R.finite) return R;

    solution u_sol = A.getSolution(u, solVector);

    // Errors, reference rule of the mode's own geometry. No factory:
    // ev.integral would silently fall back to full-box Gauss on cut cells.
    {
        gsStopwatch swErr;
        gsExprEvaluator<real_t> ev(A);                 // shares A's exprData/domain; NOT its factory
        auto u_ex = ev.getVariable(u_exact, G);
        real_t l2sq, h1ssq;
        {
            gsTetClip::VolumeQuadratureScope<gsExprEvaluator<real_t> > s(ev, refFactory);
            GISMO_ENSURE(ev.hasCustomQuadrature(), "solvePoissonOnDomain: reference-rule scope "
                        "failed to install its factory.");
            l2sq  = ev.integral((u_ex - u_sol).sqNorm() * meas(G));
            h1ssq = ev.integral((igrad(u_ex) - igrad(u_sol,G)).sqNorm() * meas(G));   // SEMINORM only
        }
        R.L2  = math::sqrt(math::max(l2sq,  (real_t)0));
        R.H1s = math::sqrt(math::max(h1ssq, (real_t)0));
        R.t_err = swErr.stop();
    }

    // Conditioning estimate.
    {
        const gsSparseMatrix<real_t> & K = A.matrix();
        const index_t N = A.numDofs();

        const gsSparseMatrix<real_t> diff = K - gsSparseMatrix<real_t>(K.transpose());
        real_t maxDiff = 0, maxK = 0;
        for (index_t c = 0; c != diff.outerSize(); ++c)
            for (gsSparseMatrix<real_t>::iterator it(diff,c); it; ++it)
                maxDiff = math::max(maxDiff, math::abs(it.value()));
        for (index_t c = 0; c != K.outerSize(); ++c)
            for (gsSparseMatrix<real_t>::iterator it(K,c); it; ++it)
                maxK = math::max(maxK, math::abs(it.value()));
        R.asym = (0 == maxK) ? (real_t)0 : maxDiff/maxK;

        if (N <= kappaDense)
        {
            const gsMatrix<real_t> dense = K.toDense();
            gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> es(dense, gsEigen::EigenvaluesOnly);
            GISMO_ENSURE(gsEigen::Success == es.info(), "solvePoissonOnDomain: dense eigensolver "
                        "failed.");
            R.lmin = es.eigenvalues().minCoeff();
            R.lmax = es.eigenvalues().maxCoeff();
            R.kappaMethod = "dense";
            R.kappaConverged = true;
        }
        else
        {
            gsVector<real_t> v0(N);
            for (index_t i = 0; i != N; ++i) v0[i] = 1.0 + 1e-3*(real_t)(i % 7);
            v0 /= v0.norm();

            gsVector<real_t> v = v0;
            real_t lambda = 0; long long itPower = 0; bool powerConverged = false;
            for (; itPower < kappaMaxIt; )
            {
                const gsVector<real_t> w = K*v;
                const real_t nrm = w.norm();
                GISMO_ENSURE(nrm > 0, "solvePoissonOnDomain: power iteration hit a zero vector.");
                v = w/nrm;
                const real_t lambdaNew = v.dot(K*v);
                const bool converged = (itPower > 0)
                                      && (math::abs(lambdaNew-lambda) <= 1e-8*math::abs(lambdaNew));
                lambda = lambdaNew;
                ++itPower;
                if (converged) { powerConverged = true; break; }
            }
            R.lmax = lambda;

            v = v0;
            real_t lambdaMin = 0; long long itInverse = 0; bool inverseConverged = false;
            for (; itInverse < kappaMaxIt; )
            {
                const gsVector<real_t> w = solver.solve(v);
                const real_t nrm = w.norm();
                GISMO_ENSURE(nrm > 0 && w.allFinite(), "solvePoissonOnDomain: inverse iteration "
                            "produced a non-finite or zero vector.");
                v = w/nrm;
                const real_t rq = v.dot(K*v)/v.dot(v);              // Rayleigh quotient: sign survives
                const bool converged = (itInverse > 0)
                                      && (math::abs(rq-lambdaMin) <= 1e-8*math::abs(rq));
                lambdaMin = rq;
                ++itInverse;
                if (converged) { inverseConverged = true; break; }
            }
            R.lmin = lambdaMin;
            R.kappaMethod = "iter";
            R.kItPower = itPower; R.kItInverse = itInverse;
            R.kappaConverged = powerConverged && inverseConverged;
        }
        R.kappaIndef = !(R.lmin > 0);
        R.kappa = R.kappaIndef ? (real_t)0 : R.lmax/R.lmin;
    }

    return R;
}

/// One `POISSON`/`POISSON-SKIP` data row's bookkeeping, kept per (case,
/// mode) across the r loop so eocL2/eocH1 can read the PREVIOUS solved r's
/// errors (log2(e_{r-1}/e_r)); a skipped or failed r leaves \a solved
/// false, which is what prints `-` for the eoc of the FOLLOWING row too.
struct PoissonHist
{
    bool solved = false;
    index_t n = 0, ndof = 0;
    real_t L2 = 0, H1s = 0;
    real_t kappa = 0; bool kappaIndef = true, kappaConverged = false;
};

std::string fmtFixed2(real_t v)
{
    std::ostringstream os; os << std::fixed << std::setprecision(2) << v; return os.str();
}

/// Prints one `POISSON` data row. \a prev is the previous r's
/// \ref PoissonHist for this (case, mode) (r = 0 passes a default-
/// constructed, unsolved one), used only for the eocL2/eocH1 `-` fallback.
/// \a statusSuffix is appended verbatim (e.g. `" FAIL"` for a zero-row or
/// non-finite outcome, where L2/H1s/eoc/kappa are all left at their
/// not-computed defaults); empty for a normal, fully solved row.
void printPoissonRow(const std::string & caseName, const std::string & mode, index_t r, index_t n,
                     const PoissonRunResult & R, index_t nCut, index_t nGhost,
                     const PoissonStats & stats, const PoissonHist & prev,
                     real_t t_tab, real_t t_wall, const std::string & statusSuffix = "")
{
    const bool haveEoc = prev.solved && R.finite && (0==R.zeroRows);
    const std::string eocL2 = haveEoc ? fmtFixed2(math::log(prev.L2/R.L2)/math::log((real_t)2))  : "-";
    const std::string eocH1 = haveEoc ? fmtFixed2(math::log(prev.H1s/R.H1s)/math::log((real_t)2)) : "-";
    const std::string kIt = ("iter" == R.kappaMethod)
        ? (std::to_string(R.kItPower) + "/" + std::to_string(R.kItInverse)) : "-";

    gsInfo << "POISSON case=" << caseName << " mode=" << mode << " r=" << r << " n=" << n
          << " ndof=" << R.ndof << " nCut=" << nCut << " nGhost=" << nGhost
          << " L2=" << gsTetClip::fmtSci(R.L2) << " H1s=" << gsTetClip::fmtSci(R.H1s)
          << " eocL2=" << eocL2 << " eocH1=" << eocH1
          << " nq_vol=" << stats.nqVol << " nq_bdr=" << stats.nqBdr
          << " minw=" << gsTetClip::fmtSci(stats.minw)
          << " kappa=" << (R.kappaIndef ? std::string("(indef)") : gsTetClip::fmtSci(R.kappa))
          << " kappaMethod=" << (R.kappaMethod.empty() ? "-" : R.kappaMethod)
          << " kIt=" << kIt
          << " kConv=" << (R.kappaMethod.empty() ? "-" : (R.kappaConverged ? "1" : "0"))
          << " lmin=" << gsTetClip::fmtSci(R.lmin) << " lmax=" << gsTetClip::fmtSci(R.lmax)
          << " asym=" << gsTetClip::fmtSci(R.asym)
          << " zeroRows=" << R.zeroRows << " finite=" << (R.finite ? 1 : 0)
          << " t_tab=" << gsTetClip::fmtSci(t_tab) << "s"
          << " t_asm=" << gsTetClip::fmtSci(R.t_asm) << "s"
          << " t_solve=" << gsTetClip::fmtSci(R.t_solve) << "s"
          << " t_err=" << gsTetClip::fmtSci(R.t_err) << "s"
          << " t_wall=" << gsTetClip::fmtSci(t_wall) << "s" << statusSuffix << "\n";
}

/// Prints the r/n/ndof/L2/eocL2/H1s/eocH1/kappa/eocK table of every SOLVED
/// row of \a hist for one (case, mode), immediately before its `EOC`
/// verdict. eocK_r = log(kappa_r/kappa_{r-1}) / log(h_{r-1}/h_r), which
/// collapses to log2(kappa_r/kappa_{r-1}) under this table's own uniform-
/// bisection convention (the same one eocL2/eocH1 already assume); an
/// h^-2 conditioning growth therefore shows as eocK ~= 2. Left as `-` for
/// r=0, an indefinite kappa at either endpoint, or either endpoint's
/// iterative eigenvalue estimate not having converged (kappa would then be
/// an underestimate, making the ratio meaningless).
void printPoissonTable(const std::string & caseName, const std::string & mode,
                       const std::vector<PoissonHist> & hist)
{
    gsInfo << "\n=== Poisson convergence  case=" << caseName << "  mode=" << mode << " ===\n";
    gsInfo << std::right << std::setw(4) << "r" << std::setw(6) << "n" << std::setw(8) << "ndof"
          << std::setw(13) << "L2" << std::setw(8) << "eocL2"
          << std::setw(13) << "H1s" << std::setw(8) << "eocH1"
          << std::setw(13) << "kappa" << std::setw(8) << "eocK" << "\n";
    bool havePrev = false; PoissonHist prev;
    for (size_t r = 0; r != hist.size(); ++r)
    {
        const PoissonHist & h = hist[r];
        if (!h.solved) continue;
        const std::string eocL2 = (havePrev && prev.solved)
            ? fmtFixed2(math::log(prev.L2/h.L2)/math::log((real_t)2)) : "-";
        const std::string eocH1 = (havePrev && prev.solved)
            ? fmtFixed2(math::log(prev.H1s/h.H1s)/math::log((real_t)2)) : "-";
        const std::string kappaStr = h.kappaIndef ? std::string("(indef)") : gsTetClip::fmtSci(h.kappa,3);
        const bool haveEocK = havePrev && prev.solved
            && !prev.kappaIndef && !h.kappaIndef && prev.kappaConverged && h.kappaConverged;
        const std::string eocK = haveEocK
            ? fmtFixed2(math::log(h.kappa/prev.kappa)/math::log((real_t)2)) : "-";
        gsInfo << std::setw(4) << r << std::setw(6) << h.n << std::setw(8) << h.ndof
              << std::setw(13) << gsTetClip::fmtSci(h.L2,3) << std::setw(8) << eocL2
              << std::setw(13) << gsTetClip::fmtSci(h.H1s,3) << std::setw(8) << eocH1
              << std::setw(13) << kappaStr << std::setw(8) << eocK << "\n";
        prev = h; havePrev = true;
    }
}

/// Prints the `EOC` verdict line for one (case, mode), and
/// returns whether it PASSed (FAIL is the only verdict this function can
/// fail the study with; REPORT never does). \a cap is the last-pair cap
/// (3 for every mode reached through this function; tchakaloff never solves
/// and never calls this). \a hist is indexed by r, size rMax+1.
bool printEocVerdict(const std::string & caseName, const std::string & mode, index_t p,
                     index_t rMax, index_t cap, const std::vector<PoissonHist> & hist)
{
    const index_t lo = cap-1, hi = cap;
    const bool reached = (rMax >= cap) && (lo >= 0) && (hi < (index_t)hist.size());
    const bool pairOk = reached && hist[lo].solved && hist[hi].solved;

    std::string L2last = "-", H1last = "-";
    real_t eocL2 = 0, eocH1 = 0;
    if (pairOk)
    {
        eocL2 = math::log(hist[lo].L2/hist[hi].L2)/math::log((real_t)2);
        eocH1 = math::log(hist[lo].H1s/hist[hi].H1s)/math::log((real_t)2);
        L2last = fmtFixed2(eocL2); H1last = fmtFixed2(eocH1);
    }

    std::string verdict = "REPORT";
    bool failsStudy = false;
    if ("clip" == mode && 2 == p)
    {
        if (!reached)
        {
            verdict = "REPORT";
        }
        else if (!pairOk)
        {
            verdict = "FAIL"; failsStudy = true;
        }
        else
        {
            const bool ok = (eocL2 >= 2.7 && eocL2 <= 3.3) && (eocH1 >= 1.7 && eocH1 <= 2.3);
            verdict = ok ? "PASS" : "FAIL";
            failsStudy = !ok;
        }
    }
    else if ("momrule" == mode || "algoim" == mode)
    {
        verdict = "REPORT";
    }
    else // p != 2 (any mode)
    {
        verdict = "REPORT";
    }

    gsInfo << "EOC case=" << caseName << " mode=" << mode << " p=" << p
          << " lastPair=" << lo << "-" << hi
          << " L2last=" << L2last << " H1last=" << H1last << " " << verdict << "\n";
    return !failsStudy;
}

/// Runs the 3D immersed Poisson study over `--case`/`--mode`/`-r`: for each
/// (case, mode, r), builds the mode's quadrature (R1, reusing \ref
/// buildCellRules / \ref runAlgoimPass unchanged from `--study volume`),
/// gates it on those SAME objects (R2, \ref gateTetMode / \ref gateAlgoim)
/// and prints the row (R3, \ref printGateRow); only on a passing gate
/// solves the Poisson problem (\ref solvePoissonOnDomain) and prints one
/// `POISSON` row. See the file header for the per-mode behaviour, the
/// pre-asymptotic note and the EoC bar; tchakaloff never solves (its
/// `--study volume` gate rows are DEFERRED, never PASS).
bool runPoissonStudy(const Config & cfg)
{
    const bool ghostOn = (1 == cfg.ghostOn);
    const real_t gammaEff = (cfg.gamma <= 0) ? 6.0*(real_t)(cfg.p+1)*(real_t)(cfg.p+1) : cfg.gamma;
    const real_t gtEff    = (cfg.ghostCoef < 0) ? math::pow(10.0, -(real_t)(cfg.p+1)) : cfg.ghostCoef;
    const index_t cap = 3;
    const gsTetClip::BgMapKind kind = geoKind(cfg.geo);
    const bool ghostUsesG = ("G" == cfg.ghostMap);

    gsInfo << "POISSON-STUDY p=" << cfg.p << " n0=" << cfg.n0 << " rMax=" << cfg.rMax
          << " gamma=" << gammaEff << " ghost=" << (ghostOn?1:0) << " ghostCoef=" << gtEff
          << " u=sin(pi*x/2)*cos(pi*y/3)*exp(z/2)\n";
    if (gsTetClip::BgMapKind::Identity != kind)
        gsInfo << "POISSON-GEO geo=" << cfg.geo << " eps=" << cfg.eps
              << " ghostMap=" << cfg.ghostMap << " sigmaNPerDir=11\n";

    gsFunctionExpr<real_t> u_exact("sin(pi*x/2)*cos(pi*y/3)*exp(z/2)", 3);
    gsFunctionExpr<real_t> f_rhs("(13*pi^2/36-1/4)*sin(pi*x/2)*cos(pi*y/3)*exp(z/2)", 3);
    gsFunctionExpr<real_t> phiSphere("(x-0.03)^2+(y+0.02)^2+(z-0.01)^2-0.3025", 3);

    std::vector<std::string> cases;
    if ("all" == cfg.caseName) { cases.push_back("sphere"); cases.push_back("rotcube"); }
    else cases.push_back(cfg.caseName);

    std::vector<std::string> modes;
    if ("all" == cfg.mode) modes = { "clip", "tchakaloff", "momrule", "algoim" };
    else modes.push_back(cfg.mode);

    bool studyOk = true;

    for (const std::string & caseName : cases)
    {
        const gsTetClip::PhysTetMesh phys = loadPhysMesh(cfg, caseName);
        const real_t V_mesh = gsTetClip::meshVolumeExact(phys);
        const real_t A_mesh = gsTetClip::unclippedBoundaryArea(phys);
        const real_t V_exact = ("sphere" == caseName) ? sphereVolumeExact()  : rotcubeVolumeExact();
        const real_t A_exact = ("sphere" == caseName) ? sphereAreaExact()    : rotcubeAreaExact();

        real_t pbScale = 1.0;
        if (gsTetClip::BgMapKind::Identity != kind)
            for (const gsTetClip::Vec3 & x : phys.m.P)
                pbScale = math::max(pbScale, math::max(std::abs(x[0]),
                          math::max(std::abs(x[1]), std::abs(x[2]))));

        for (const std::string & mode : modes)
        {
            if ("algoim" == mode && "rotcube" == caseName)
            {
                gsInfo << "POISSON-SKIP case=rotcube mode=algoim reason=sphere-only  INFO\n";
                continue;
            }

            if ("tchakaloff" == mode)
            {
                for (index_t r = 0; r <= cfg.rMax; ++r)
                {
                    const gsTetClip::Grid3 grid = makeGrid(phys.m, cfg.n0, r, "runPoissonStudy");
                    gsInfo << "POISSON-SKIP case=" << caseName << " mode=tchakaloff r=" << r
                          << " n=" << grid.n << " reason=gate-deferred\n";
                }
                gsInfo << "EOC case=" << caseName << " mode=tchakaloff p=" << cfg.p
                      << " lastPair=" << (cap-1) << "-" << cap << " L2last=- H1last=- REPORT\n";
                continue;
            }

            memory::shared_ptr<gsMeshSignedDist<real_t> > phiH;
            if ("momrule" == mode) phiH = buildMeshSignedDist(phys.m);

            std::vector<PoissonHist> hist(cfg.rMax+1);
            std::vector<real_t> tAsmHist(cfg.rMax+1, 0), tErrHist(cfg.rMax+1, 0), tWallHist(cfg.rMax+1, 0);
            std::vector<long long> nqVolHist(cfg.rMax+1, 0);

            for (index_t r = 0; r <= cfg.rMax; ++r)
            {
                gsStopwatch swWall;
                resetPeakRss();
                const gsTetClip::Grid3 grid = makeBoxGrid(cfg.n0, r);
                const index_t n = grid.n;
                const real_t h = grid.h;

                index_t nCut = 0, nGhost = 0;
                PoissonStats stats;
                real_t t_tab = 0;
                bool gatePass = true;

                gsMultiPatch<real_t> mp; gsMultiBasis<real_t> mb;
                makeBackground(grid, cfg.p, kind, cfg.eps, mp, mb);
                const gsGeometry<real_t> & bgMap = mp.patch(0);
                gsTensorBSplineBasis<3,real_t> * tbs =
                    dynamic_cast<gsTensorBSplineBasis<3,real_t>*>(&mb.basis(0));
                GISMO_ENSURE(tbs, "runPoissonStudy: background basis is not a tensor B-spline basis.");

                const std::string dumpFile = cfg.dumpSystem.empty() ? std::string()
                    : cfg.dumpSystem + "_" + caseName + "_" + mode + "_r" + std::to_string(r) + ".bin";
                const gsMultiPatch<real_t> mpId = gsTetClip::identityBoxGeometry(grid);

                PoissonRunResult R;

                if ("algoim" != mode)
                {
                    // --- Pull-back (gsImmersedPullback.h), always before
                    // any ClipStreamer or assembly. ---
                    gsTetClip::PullbackStats pbStats;
                    const gsTetClip::ParamTetMesh param =
                        gsTetClip::pullBack(phys, bgMap, kind, grid, pbStats);

                    if (gsTetClip::BgMapKind::Identity != kind)
                    {
                        const real_t tol = 1e-12*pbScale;
                        const bool pbPass = (0 == pbStats.nFailed && 0 == pbStats.nOutside &&
                                            0 == pbStats.nInverted && pbStats.maxResidual <= tol);
                        gsInfo << "PULLBACK case=" << caseName << " geo=" << cfg.geo << " r=" << r
                              << " n=" << n << " nVerts=" << pbStats.nVerts << " nFailed=" << pbStats.nFailed
                              << " nOutside=" << pbStats.nOutside << " nInverted=" << pbStats.nInverted
                              << " maxResidual=" << gsTetClip::fmtSci(pbStats.maxResidual)
                              << " tol=" << gsTetClip::fmtSci(tol)
                              << " minParamTetVol=" << gsTetClip::fmtSci(pbStats.minParamTetVol)
                              << " " << (pbPass ? "PASS" : "FAIL") << "\n";
                        if (!pbPass)
                        {
                            gsWarn << "PULLBACK-ABORT case=" << caseName << " geo=" << cfg.geo
                                  << " r=" << r << ": pull-back failed its acceptance check\n";
                            return false;
                        }
                    }

                    makeGrid(param.m, cfg.n0, r, "runPoissonStudy");   // strict-inside ENSURE, parametric

                    const memory::shared_ptr<gsTetClip::TetMesh> Mclip =
                        memory::make_shared(new gsTetClip::TetMesh(param.m));

                    // --- Quadrature build + gate: tet modes. ---
                    gsStopwatch swTab;
                    const memory::shared_ptr<gsTetClip::ClipStreamer> S =
                        memory::make_shared(new gsTetClip::ClipStreamer(Mclip, grid, cfg.p));

                    // Nanson-mapped boundary source (non-identity G only);
                    // also the Nitsche bdrSrc below.
                    memory::shared_ptr<const gsTetClip::BdrCellSource> pullbackBdrSrc;
                    if (gsTetClip::BgMapKind::Identity != kind)
                        pullbackBdrSrc = memory::make_shared(new gsTetClip::PullbackBdrSource(*S, bgMap));

                    CellRulePass P;
                    bool solveGatePass;
                    if (gsTetClip::BgMapKind::Identity == kind)
                    {
                        P = buildCellRules(*S, mode, cfg.p, phiH.get(), true);
                        const GateResult G =
                            gateTetMode(P, *S, cfg.p, V_mesh, A_mesh, V_exact, A_exact, caseName, r);
                        for (const GateRow & row : G.rows) printGateRow(G, row, false);
                        t_tab = swTab.stop();

                        // fluxmom_global/fluxmom_cell are Required for every mode
                        // (gateTetMode), including momrule: they count fully
                        // toward `--study volume`'s GATE SUMMARY. Here, though,
                        // momrule solves iff every OTHER Required row passed --
                        // its served normal is only a pseudonormal at polyhedral
                        // edges/vertices (buildCellRules), so a fluxmom-only
                        // failure is the expected, already-documented gap
                        // (file header), not a reason to withhold the solve.
                        solveGatePass = G.requiredPass;
                        if (!solveGatePass && "momrule" == mode)
                        {
                            bool onlyFluxmomFail = true;
                            for (const GateRow & row : G.rows)
                                if (GateRow::Required == row.kind && !row.pass &&
                                    "fluxmom_global" != row.check && "fluxmom_cell" != row.check)
                                    onlyFluxmomFail = false;
                            if (onlyFluxmomFail)
                            {
                                solveGatePass = true;
                                gsInfo << "GATE-NOTE mode=momrule fluxmom=FAIL solve=report\n";
                            }
                        }

                        if (!solveGatePass)
                        {
                            gsInfo << "POISSON-SKIP case=" << caseName << " mode=" << mode << " r=" << r
                                  << " n=" << n << " reason=gate-FAIL\n";
                            gatePass = false;
                        }
                    }
                    else
                    {
                        // Curved background: only --mode clip reaches here
                        // (main() rejects the rest); buildCellRules still runs
                        // (statsFromTetPass below needs P.log). The physical-
                        // space gate recomputes volume/area/flux/
                        // fluxmom_global directly from S and bgMap
                        // (gatePhysical), since CellLog's totals are
                        // parametric here; the P1 momrule fluxmom exception
                        // does not apply -- momrule is rejected outright for
                        // a curved G.
                        P = buildCellRules(*S, "clip", cfg.p, nullptr, true);
                        const GateResult G = gatePhysical(phys, pbStats, bgMap, kind, *S,
                                                          *pullbackBdrSrc, cfg.p, V_exact, A_exact,
                                                          caseName, cfg.geo, r);
                        for (const GateRow & row : G.rows) printGateRow(G, row, false);
                        t_tab = swTab.stop();
                        solveGatePass = G.requiredPass;
                        if (!solveGatePass)
                        {
                            gsInfo << "POISSON-SKIP case=" << caseName << " mode=" << mode << " r=" << r
                                  << " n=" << n << " reason=gate-FAIL\n";
                            gatePass = false;
                        }
                    }

                    if (solveGatePass)
                    {
                        stats = statsFromTetPass(P, *S, cfg.p);

                        memory::shared_ptr<const gsTetClip::CellIndex> idx = S->index();
                        GISMO_ENSURE(idx.get() == S->index().get(), "runPoissonStudy: rule/domain "
                                    "CellIndex identity broken.");

                        memory::shared_ptr<gsTetClip::TetClipSignDomain> tdom =
                            memory::make_shared(new gsTetClip::TetClipSignDomain(grid, idx->status,
                                                                                 (short_t)cfg.p));
                        tdom->numGhostFaces(); tdom->numElementsBdr(boundary::none);  // warm, single-threaded
                        nGhost = (index_t)tdom->numGhostFaces();
                        nCut   = (index_t)tdom->numElementsBdr(boundary::none);

                        memory::shared_ptr<const gsTetClip::VolCellSource> volSrc;
                        memory::shared_ptr<const gsTetClip::BdrCellSource> bdrSrc;
                        if (gsTetClip::BgMapKind::Identity == kind)
                        {
                            if ("clip" == mode) { volSrc = S; bdrSrc = S; }
                            else                { volSrc = P.vol; bdrSrc = P.bdr; }
                        }
                        else
                        {
                            volSrc = S;
                            bdrSrc = pullbackBdrSrc;
                        }

                        gsExprAssembler<real_t>::QuadratureFactory refFactory =
                            gsTetClip::makeVolLookupFactory(idx,
                                memory::shared_ptr<const gsTetClip::VolCellSource>(S),
                                gsVector<index_t>::Constant(3, cfg.p+3));

                        real_t hNitsche = h;
                        if (gsTetClip::BgMapKind::Identity != kind)
                        {
                            const real_t sigmaMinJ = gsTetClip::minSingularValueJ(bgMap, 11);
                            hNitsche = h * sigmaMinJ;
                            gsInfo << "NITSCHE-H case=" << caseName << " r=" << r
                                  << " h=" << gsTetClip::fmtSci(h) << " sigmaMinJ=" << gsTetClip::fmtSci(sigmaMinJ)
                                  << " hNitsche=" << gsTetClip::fmtSci(hNitsche) << "\n";
                        }

                        memory::shared_ptr<gsTrimmedDomain<3,real_t> > domBase = tdom;
                        R = solvePoissonOnDomain(mp, mb, mpId, domBase, idx, volSrc, bdrSrc, refFactory,
                                                 cfg.p, h, hNitsche, gammaEff, gtEff, ghostOn, ghostUsesG,
                                                 cfg.kappaDense, cfg.kappaMaxIt, u_exact, f_rhs, dumpFile);
                    }
                }
                else
                {
                    // Pulled-back predicate from the ACTUAL background map,
                    // not cfg.geo, so identity cannot drift.
                    const bool pulled = !gsTetClip::isIdentityMap(bgMap);

                    if (!pulled)
                        makeGrid(phys.m, cfg.n0, r, "runPoissonStudy");   // strict-inside ENSURE, identity only

                    // phi-hat = phi o G on the pulled-back path, else the
                    // physical phiSphere.
                    memory::unique_ptr<PulledBackSphere> phiHatPtr;
                    if (pulled) phiHatPtr = memory::make_unique(new PulledBackSphere(bgMap));
                    const gsFunction<real_t> & phiUse = pulled
                        ? static_cast<const gsFunction<real_t> &>(*phiHatPtr)
                        : static_cast<const gsFunction<real_t> &>(phiSphere);

                    real_t Lhat = 4.0;   // identity default (LipschitzConstant of the reference rule)
                    if (pulled)
                    {
                        Lhat = pulledBackSphereLipschitz(bgMap, kind, cfg.eps);
                        checkPulledBackSphereLipschitz(*phiHatPtr, bgMap, cfg.geo, r, Lhat);
                    }

                    // --- Quadrature build + gate: algoim. ---
                    gsStopwatch swTab;
                    memory::shared_ptr<gsTetClip::VolCellTable> volTab =
                        memory::make_shared(new gsTetClip::VolCellTable());
                    memory::shared_ptr<gsTetClip::BdrCellTable> bdrTab =
                        memory::make_shared(new gsTetClip::BdrCellTable());
                    const AlgoimPass AP = runAlgoimPass(grid, cfg.p, phiUse, volTab.get(), bdrTab.get(),
                                                        pulled ? &bgMap : nullptr);
                    const GateResult G = gateAlgoim(AP, cfg.p, caseName, r, n);
                    for (const GateRow & row : G.rows) printGateRow(G, row, false);
                    t_tab = swTab.stop();
                    // Algoim's gate rows are report-only: they never skip the solve.

                    // algoim-index CellIndex from the domain's own Lobatto
                    // leaf signs (the assembler/evaluator iterate THIS
                    // classification, not Algoim's own).
                    memory::shared_ptr<gsImplicitTrimmedDomain<3,real_t> > adom =
                        memory::make_shared(new gsImplicitTrimmedDomain<3,real_t>(phiUse, *tbs));
                    GISMO_ENSURE(1 == adom->numLevels(), "runPoissonStudy: gsImplicitTrimmedDomain "
                                "produced more than one kd-tree level.");
                    adom->numGhostFaces(); adom->numElementsBdr(boundary::none);  // warm, single-threaded
                    nGhost = (index_t)adom->numGhostFaces();
                    nCut   = (index_t)adom->numElementsBdr(boundary::none);

                    memory::shared_ptr<gsTetClip::CellIndex> algoimIdx =
                        memory::make_shared(new gsTetClip::CellIndex());
                    algoimIdx->grid = grid;
                    gsTetClip::gridLines(grid, algoimIdx->X, algoimIdx->Y, algoimIdx->Z);
                    const size_t N3 = (size_t)n*(size_t)n*(size_t)n;
                    algoimIdx->status.assign(N3, -2);   // sentinel: must be overwritten exactly once

                    auto markCell = [&](const gsVector<real_t> & lo_, const gsVector<real_t> & hi_,
                                        int status)
                    {
                        const gsVector<real_t> mid = 0.5*(lo_+hi_);
                        const size_t id = algoimIdx->cellOfPoint(mid[0], mid[1], mid[2]);
                        GISMO_ENSURE(-2 == algoimIdx->status[id], "runPoissonStudy: algoim cell id "
                                    << id << " visited twice while classifying domain leaf signs.");
                        algoimIdx->status[id] = status;
                    };
                    {
                        gsDomain<real_t>::iterator e = adom->end<InteriorSign>();
                        for (gsDomain<real_t>::iterator it = adom->beginInterior(); it < e; ++it)
                            markCell(it.lowerCorner(), it.upperCorner(), (int)gsTetClip::Full);
                    }
                    {
                        gsDomain<real_t>::iterator e = adom->end<BoundarySign>();
                        for (gsDomain<real_t>::iterator it = adom->begin<BoundarySign>(); it < e; ++it)
                            markCell(it.lowerCorner(), it.upperCorner(), (int)gsTetClip::Cut);
                    }
                    {
                        gsDomain<real_t>::iterator e = adom->end<ExteriorSign>();
                        for (gsDomain<real_t>::iterator it = adom->begin<ExteriorSign>(); it < e; ++it)
                            markCell(it.lowerCorner(), it.upperCorner(), (int)gsTetClip::Empty);
                    }
                    for (size_t id = 0; id != N3; ++id)
                        GISMO_ENSURE(-2 != algoimIdx->status[id], "runPoissonStudy: algoim cell id "
                                    << id << " was never visited while classifying domain leaf signs.");

                    // ALGOIM-CONSISTENCY: a=lostBdr, b=lostVol, c=fullNotBox, d=cutNoVol
                    // (see the file header for their definitions).
                    long long lostBdr=0, lostVol=0, fullNotBox=0, cutNoVol=0;
                    const real_t h3 = h*h*h;
                    for (size_t id = 0; id != N3; ++id)
                    {
                        const bool domCut   = (gsTetClip::Cut   == algoimIdx->status[id]);
                        const bool domFull  = (gsTetClip::Full  == algoimIdx->status[id]);
                        const bool domEmpty = (gsTetClip::Empty == algoimIdx->status[id]);

                        if (bdrTab->cell[id].nodes.cols() > 0 && !domCut) ++lostBdr;
                        if (volTab->cell[id].nodes.cols() > 0 && domEmpty) ++lostVol;
                        if (domFull)
                        {
                            gsTetClip::KahanSum ks;
                            for (index_t c = 0; c != volTab->cell[id].weights.size(); ++c)
                                ks.add(volTab->cell[id].weights[c]);
                            if (math::abs(ks.value()-h3) > 1e-12*h3) ++fullNotBox;
                        }
                        if (domCut && 0 == volTab->cell[id].nodes.cols()) ++cutNoVol;
                    }
                    gsInfo << "ALGOIM-CONSISTENCY case=sphere r=" << r << " lostBdr=" << lostBdr
                          << " lostVol=" << lostVol << " fullNotBox=" << fullNotBox
                          << " cutNoVol=" << cutNoVol << "  INFO\n";

                    // Stats: nq_vol/minw over algoim-index Cut cells (what the solve
                    // integrates); nq_bdr sums every Algoim surface node over ALL n^3
                    // cells (not restricted to algoim-index Cut), so it equals what the
                    // solve integrates only when ALGOIM-CONSISTENCY's lostBdr is 0.
                    {
                        long long nv = 0, nb = 0, numFullA = 0;
                        real_t minw = std::numeric_limits<real_t>::infinity();
                        for (size_t id = 0; id != N3; ++id)
                        {
                            if (gsTetClip::Cut == algoimIdx->status[id])
                            {
                                nv += volTab->cell[id].nodes.cols();
                                for (index_t c = 0; c != volTab->cell[id].weights.size(); ++c)
                                    minw = math::min(minw, volTab->cell[id].weights[c]);
                            }
                            else if (gsTetClip::Full == algoimIdx->status[id]) ++numFullA;
                            nb += bdrTab->cell[id].nodes.cols();
                        }
                        nv += numFullA*(long long)(cfg.p+1)*(long long)(cfg.p+1)*(long long)(cfg.p+1);
                        stats.nqVol = nv; stats.nqBdr = nb; stats.minw = minw;
                    }

                    // Reference table for the error integrals: algoim-index Cut cells only.
                    memory::shared_ptr<gsTetClip::VolCellTable> refTab =
                        memory::make_shared(new gsTetClip::VolCellTable());
                    refTab->cell.resize(N3);
                    {
                        gsOptionList oRef = algoimOptions(-1, /*quA*/2.0, /*quB*/1, /*maxDepth*/1);
                        // phi = |x-c|^2 - R^2 is NOT 1-Lipschitz: |grad phi| = 2|x-c| <=
                        // 2*(sqrt(3)+|c|) ~= 3.54 on [-1,1]^3. The default LipschitzConstant=1.0
                        // would under-bound the maxDepth=1 box classifier; 4.0 is a valid global
                        // bound on --geo identity. On a curved background, Lhat (computed above by
                        // pulledBackSphereLipschitz) bounds |grad(phi-hat)| instead.
                        oRef.setReal("LipschitzConstant", Lhat);
                        gsAlgoimAdaptiveRule<real_t> refRule(phiUse, (short_t)cfg.p, oRef);
                        std::vector<real_t> X,Y,Z; gsTetClip::gridLines(grid, X, Y, Z);
                        for (size_t id = 0; id != N3; ++id)
                        {
                            if (gsTetClip::Cut != algoimIdx->status[id]) continue;
                            index_t i,j,k; algoimIdx->ijk(id, i,j,k);
                            gsVector<real_t> lower(3), upper(3);
                            lower << X[i],Y[j],Z[k]; upper << X[i+1],Y[j+1],Z[k+1];
                            refRule.mapTo(lower, upper, refTab->cell[id].nodes, refTab->cell[id].weights);
                        }
                    }
                    gsExprAssembler<real_t>::QuadratureFactory refFactory =
                        gsTetClip::makeVolLookupFactory(algoimIdx, refTab,
                                                        gsVector<index_t>::Constant(3, cfg.p+3));

                    memory::shared_ptr<const gsTetClip::CellIndex> algoimIdxConst = algoimIdx;
                    // volTab stays PARAMETRIC (meas(G) is in every volume expression);
                    // bdrTab is wrapped in a PullbackBdrSource for assembly on the
                    // pulled-back path (Nanson), so it is applied exactly once here,
                    // independently of runAlgoimPass's own totals pass above.
                    memory::shared_ptr<const gsTetClip::VolCellSource> volSrc = volTab;
                    memory::shared_ptr<const gsTetClip::BdrCellSource> bdrSrc = pulled
                        ? memory::shared_ptr<const gsTetClip::BdrCellSource>(
                              memory::make_shared(new gsTetClip::PullbackBdrSource(*bdrTab, bgMap)))
                        : memory::shared_ptr<const gsTetClip::BdrCellSource>(bdrTab);
                    memory::shared_ptr<gsTrimmedDomain<3,real_t> > domBase = adom;

                    real_t hNitsche = h;
                    if (pulled)
                    {
                        const real_t sigmaMinJ = gsTetClip::minSingularValueJ(bgMap, 11);
                        hNitsche = h * sigmaMinJ;
                        gsInfo << "NITSCHE-H case=" << caseName << " r=" << r
                              << " h=" << gsTetClip::fmtSci(h) << " sigmaMinJ=" << gsTetClip::fmtSci(sigmaMinJ)
                              << " hNitsche=" << gsTetClip::fmtSci(hNitsche) << "\n";
                    }

                    R = solvePoissonOnDomain(mp, mb, mpId, domBase, algoimIdxConst, volSrc, bdrSrc, refFactory,
                                             cfg.p, h, hNitsche, gammaEff, gtEff, ghostOn, ghostUsesG,
                                             cfg.kappaDense, cfg.kappaMaxIt, u_exact, f_rhs, dumpFile);
                }

                const real_t t_wall = swWall.stop();
                const long long vmhwm = peakRssKB();

                if (!gatePass)
                {
                    studyOk = false;
                    hist[r].solved = false;
                    gsInfo << "POISSON-RSS case=" << caseName << " mode=" << mode << " r=" << r
                          << " vmhwm_kB=" << vmhwm << "\n";
                    continue;
                }

                const PoissonHist prev = (r > 0) ? hist[r-1] : PoissonHist();

                if (R.zeroRows > 0 || !R.finite)
                {
                    printPoissonRow(caseName, mode, r, n, R, nCut, nGhost, stats, prev,
                                    t_tab, t_wall, " FAIL");
                    studyOk = false; hist[r].solved = false;
                    gsInfo << "POISSON-RSS case=" << caseName << " mode=" << mode << " r=" << r
                          << " vmhwm_kB=" << vmhwm << "\n";
                    continue;
                }

                printPoissonRow(caseName, mode, r, n, R, nCut, nGhost, stats, prev, t_tab, t_wall);
                if ("iter" == R.kappaMethod && !R.kappaConverged)
                    gsWarn << "KAPPA-UNCONVERGED case=" << caseName << " mode=" << mode << " r=" << r
                          << " kIt=" << R.kItPower << "/" << R.kItInverse << "\n";
                gsInfo << "POISSON-RSS case=" << caseName << " mode=" << mode << " r=" << r
                      << " vmhwm_kB=" << vmhwm << "\n";

                hist[r].solved = true; hist[r].n = n; hist[r].ndof = R.ndof;
                hist[r].L2 = R.L2; hist[r].H1s = R.H1s;
                hist[r].kappa = R.kappa; hist[r].kappaIndef = R.kappaIndef;
                hist[r].kappaConverged = R.kappaConverged;
                tAsmHist[r] = R.t_asm; tErrHist[r] = R.t_err; tWallHist[r] = t_wall;
                nqVolHist[r] = stats.nqVol;

                // Clip budget line: right after r=2, before r=3 starts, so it
                // is a genuine PREDICTOR of the r=3 cost (not a postmortem).
                if ("clip" == mode && 2 == r && hist[1].solved && hist[2].solved
                    && tAsmHist[1] > 0 && tErrHist[1] > 0 && tWallHist[1] > 0)
                {
                    // Geometric extrapolation of the r=1 -> 2 growth factor, one step further.
                    const real_t t_asm_est  = tAsmHist[2]*tAsmHist[2]/tAsmHist[1];
                    const real_t t_err_est  = tErrHist[2]*tErrHist[2]/tErrHist[1];
                    const real_t t_wall_est = tWallHist[2]*tWallHist[2]/tWallHist[1];
                    gsInfo << "BUDGET case=" << caseName << " mode=clip r=3 n=" << (cfg.n0*8)
                          << " t_asm_est=" << gsTetClip::fmtSci(t_asm_est) << "s"
                          << " t_err_est=" << gsTetClip::fmtSci(t_err_est) << "s"
                          << " t_wall_est=" << gsTetClip::fmtSci(t_wall_est) << "s"
                          << " method=geometric nqvol1=" << nqVolHist[1] << " nqvol2=" << nqVolHist[2]
                          << "\n";
                }
            }

            printPoissonTable(caseName, mode, hist);
            const bool eocPass = printEocVerdict(caseName, mode, cfg.p, cfg.rMax, cap, hist);
            studyOk = studyOk && eocPass;
        }
    }

    gsInfo << (studyOk ? "POISSON STUDY PASS" : "POISSON STUDY FAIL") << "\n";
    return studyOk;
}

} // anonymous namespace

//----------------------------------------------------------------------------
// main
//----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    std::string study    = "check";
    std::string caseName = "all";
    std::string mode     = "all";
    index_t p    = 2;
    index_t rMax = 1;
    index_t n0   = 4;
    real_t  gamma      = -1;
    index_t ghostOn    = 1;
    real_t  ghostCoef  = -1;
    index_t kappaDense = 6000;
    index_t kappaMaxIt = 500;
    std::string geo          = "identity";
    real_t      eps           = 0.3;
    std::string meshOverride;
    std::string ghostMap      = "Gid";
    std::string dumpSystem;

    gsCmdLine cmd("Streamed lookup-rule quadrature checks, the volume/area/flux/moment gate, and "
                 "the 3D immersed Poisson solver (symmetric Nitsche + ghost penalty) for the "
                 "tet-mesh immersed quadrature machinery (gsImmersedLookupRule.h): "
                 "streaming/table bitwise equality against the reference clip "
                 "(gsTetMeshClip.h), RAII quadrature-scope save/restore plus a bitwise "
                 "ghost-penalty-matrix equality check, boundary normal-field throw behaviour "
                 "(--study check); a per-mode quadrature gate over the clip/Tchakaloff/"
                 "moment-fitting/Algoim cell rules (--study volume); and, per (case, mode, r), "
                 "that same gate followed by a Poisson solve with an EoC table (--study poisson).");
    cmd.addString("", "study", "Study to run: check | volume | poisson", study);
    cmd.addString("", "case",  "Test case: sphere | rotcube | all (= both, sphere first)", caseName);
    cmd.addString("", "mode", "Quadrature mode for --study volume|poisson: all|clip|tchakaloff|"
                  "momrule|algoim (algoim is sphere-only; ignored by --study check)", mode);
    cmd.addInt   ("k", "degree", "Background-space degree / clip-rule parameter p", p);
    cmd.addInt   ("r", "refine", "The study runs r = 0..rMax, n = n0*2^r cells per direction", rMax);
    cmd.addInt   ("",  "n0",     "Cells per direction at r = 0", n0);
    cmd.addReal  ("",  "gamma",     "Nitsche penalty coefficient (<=0: auto 6*(p+1)^2) "
                  "(--study poisson)", gamma);
    cmd.addInt   ("",  "ghost",     "Ghost penalty: 1 on, 0 off (--study poisson)", ghostOn);
    cmd.addReal  ("",  "ghostCoef", "Ghost penalty coefficient gamma_g (<0: auto 10^-(p+1)) "
                  "(--study poisson)", ghostCoef);
    cmd.addInt   ("",  "kappaDense", "Dense-eigensolver dof threshold for the conditioning "
                  "estimate (--study poisson)", kappaDense);
    cmd.addInt   ("",  "kappaMaxIt", "Power/inverse iteration cap for the iterative conditioning "
                  "estimate (--study poisson)", kappaMaxIt);
    cmd.addString("",  "geo", "Background map: identity | affine | bubble", geo);
    cmd.addReal  ("",  "eps", "Bubble amplitude, |eps| < 0.5 (--geo bubble)", eps);
    cmd.addString("",  "mesh", "Overrides the case's mesh file (e.g. a gmsh-refined mesh); "
                  "--case still selects the analytic references (--study volume|poisson, "
                  "--case sphere|rotcube only)", meshOverride);
    cmd.addString("",  "ghostMap", "(debug) map inside the ghost penalty's dnk: Gid | G", ghostMap);
    cmd.addString("",  "dumpSystem", "(debug, --study poisson) file prefix for a raw binary dump "
                  "of the assembled system; empty writes nothing", dumpSystem);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (p < 1)    { gsWarn << "-k/--degree must be >= 1\n"; return EXIT_FAILURE; }
    if (rMax < 0) { gsWarn << "-r/--refine must be >= 0\n"; return EXIT_FAILURE; }
    if (n0 < 1)   { gsWarn << "--n0 must be >= 1\n"; return EXIT_FAILURE; }
    if ("sphere" != caseName && "rotcube" != caseName && "all" != caseName)
    { gsWarn << "--case must be one of sphere|rotcube|all\n"; return EXIT_FAILURE; }
    if ("all" != mode && "clip" != mode && "tchakaloff" != mode && "momrule" != mode && "algoim" != mode)
    { gsWarn << "--mode must be one of all|clip|tchakaloff|momrule|algoim\n"; return EXIT_FAILURE; }
    if ("volume" == study && "rotcube" == caseName && "algoim" == mode)
    {
        gsWarn << "--study " << study << " --case rotcube --mode algoim: rotcube has no analytic "
                 "level set; use --case all or --case sphere\n";
        return EXIT_FAILURE;
    }
    if ("poisson" == study)
    {
        if (0 != ghostOn && 1 != ghostOn) { gsWarn << "--ghost must be 0 or 1\n"; return EXIT_FAILURE; }
        if (kappaDense < 0) { gsWarn << "--kappaDense must be >= 0\n"; return EXIT_FAILURE; }
        if (kappaMaxIt < 1) { gsWarn << "--kappaMaxIt must be >= 1\n"; return EXIT_FAILURE; }
    }
    if ("identity" != geo && "affine" != geo && "bubble" != geo)
    { gsWarn << "--geo must be one of identity|affine|bubble\n"; return EXIT_FAILURE; }
    if ("bubble" == geo && !(std::isfinite(eps) && math::abs(eps) < 0.5))
    { gsWarn << "--eps must be finite with |eps| < 0.5 for --geo bubble\n"; return EXIT_FAILURE; }
    if ("Gid" != ghostMap && "G" != ghostMap)
    { gsWarn << "--ghostMap must be Gid or G\n"; return EXIT_FAILURE; }
    if (!meshOverride.empty())
    {
        if ("volume" != study && "poisson" != study)
        { gsWarn << "--mesh requires --study volume or poisson\n"; return EXIT_FAILURE; }
        if ("sphere" != caseName && "rotcube" != caseName)
        { gsWarn << "--mesh requires --case sphere or rotcube\n"; return EXIT_FAILURE; }
        if (gsFileManager::find(meshOverride).empty())
        { gsWarn << "--mesh file '" << meshOverride << "' not found\n"; return EXIT_FAILURE; }
    }
    if (("volume" == study || "poisson" == study) && "identity" != geo)
    {
        if ("momrule" == mode || "tchakaloff" == mode)
        {
            gsWarn << "--mode " << mode << " is only supported for --geo identity\n";
            return EXIT_FAILURE;
        }
        if ("all" == mode)
        {
            gsWarn << "--mode all is not supported for --geo != identity; use --mode clip\n";
            return EXIT_FAILURE;
        }
    }
    if (!dumpSystem.empty() && "poisson" != study)
    { gsWarn << "--dumpSystem requires --study poisson\n"; return EXIT_FAILURE; }

    Config cfg;
    cfg.study = study; cfg.caseName = caseName; cfg.mode = mode; cfg.p = p; cfg.rMax = rMax; cfg.n0 = n0;
    cfg.gamma = gamma; cfg.ghostOn = ghostOn; cfg.ghostCoef = ghostCoef; cfg.kappaDense = kappaDense;
    cfg.kappaMaxIt = kappaMaxIt;
    cfg.geo = geo; cfg.eps = eps; cfg.meshOverride = meshOverride; cfg.ghostMap = ghostMap;
    cfg.dumpSystem = dumpSystem;

    bool ok = false;
    try
    {
        if ("check" == study)
        {
            ok = runCheckStudy(cfg);
        }
        else if ("volume" == study)
        {
            ok = runVolumeStudy(cfg);
        }
        else if ("poisson" == study)
        {
            ok = runPoissonStudy(cfg);
        }
        else
        {
            gsWarn << "--study must be check|volume|poisson\n";
            return EXIT_FAILURE;
        }
    }
    catch (const std::exception & e)
    {
        gsWarn << "ABORT study=" << study << " what=" << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
