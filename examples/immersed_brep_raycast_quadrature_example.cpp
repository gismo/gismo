/** @file immersed_brep_raycast_quadrature_example.cpp

    @brief Geometry loaders, exact oracles, a self-test battery and a 3D
    background-grid quadrature rule built from the ray-cast primitives of
    gsBRepRayCast.h: a uniform grid over Geometry::bg, refined by an octree
    to a per-cell rectangular leaf certified for a three-level (ray cast +
    bisected splits + plain Gauss) rule, with an exact wall term for
    axis-aligned planar boundary elements and a sign-sampling fallback for
    leaves that cannot be certified. No Algoim baseline and no plotting.

    The whole driver is SERIAL. ::detail::fitNcoef (gsBRepRayCast.h) keeps
    function-local static std::map caches, and the call chain
    ::castLine -> ::splitElement -> ::detail::fillElementFields -> fitNcoef
    can insert into them from anywhere the octree recursion reaches; running
    the recursion (or the leaf rule's repeated ::castLine calls) across
    threads would race on those caches. gsBRepSignedDist below does not
    itself touch those caches (its footpoint solve uses only gsGeometry
    evaluation and a thread_local scratch vector), but it is kept serial
    as well, with no `#pragma omp` anywhere in this file, so that the
    driver contains no OpenMP at all.

    Four closed, untrimmed spline BReps are supported (--geom):
    - cube: a BSplineCube of full side s = 0.6, rotated by
      (phiz,phiy,phix) = --rot (default 0.3,0.5,0.7) about the origin and
      translated to an off-lattice centre, its 6 boundary faces as a
      6-patch gsMultiPatch;
    - cubeAligned: an axis-aligned BSplineCube of side 0.5 centred at
      (0.5,0.5,0.5), so every face lies exactly on the plane x_k in
      {0.25,0.75} and bg = [0,1]^3 exactly;
    - sphere: a NurbsSphere of radius 0.4 at an off-lattice, off-seam centre
      (single rational patch, 8 Bezier elements, one per octant);
    - duck: filedata/breps/3D/duck_BRep.xml, a 6-patch polynomial
      biquadratic shell (original coordinates, no rescaling).

    Per-patch orientation signs come from fixOrientation (below). Its
    geometric side-pairing only matches sides ACROSS patches, so a
    single-patch BRep (the sphere) is assigned sgn = {+1} directly, ahead of
    that pairing; every BRep, single- or multi-patch, then goes through the
    same global closure/volume validation.

    The volume-moment oracle (::momentOracle) evaluates
    Integral_Omega x^a y^b z^c dV = Sum_elements Integral (S)
    F(S).(sgn.S_u x S_v)_x du dv, F = x^(a+1)/(a+1).y^b.z^c, by the
    divergence theorem, using brc::normalNumerator = w^4.sgn.(S_u x S_v) so
    that no division by w ever occurs; see its doxygen for the exact
    Gauss-node count. Cube and sphere moments have closed-form references
    (::cubeMoment, ::sphereMoment) used only by --selftest.

    --selftest runs a fixed battery over all four geometries (fixed seed
    12345 for std::mt19937) and prints one PASS/FAIL line per check; see
    ::selftestAll for the exact list and order. Every line's "max" is the
    worst observed (error / that check's own tolerance) ratio, so a line is
    PASS iff max <= 1 and the printed "tol" is always 1 -- this normalizes
    checks whose natural tolerance is per-sample or per-element (e.g.
    1e-12*Nscale(element)) onto one comparable scale.

    Example command lines:
      ./immersed_brep_raycast_quadrature_example --selftest
      ./immersed_brep_raycast_quadrature_example --geom duck
      ./immersed_brep_raycast_quadrature_example --geom cube --rot 0.1 --rot 0.2 --rot 0.3

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include "gsBRepRayCast.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

namespace {

std::string fmtSci(real_t v, int prec = 6)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

} // anonymous namespace

// =============================================================================
//  Orientation fix-up. A single-patch BRep (the sphere) is handled ahead of
//  the geometric side-pairing below, which only matches sides ACROSS
//  patches and therefore cannot itself resolve a one-patch shell -- see the
//  file doxygen for why the sphere needs it.
// =============================================================================

/// Fill \a x = S(uv) and the Jacobian columns \a xu, \a xv at a single point.
static void patchFrame(const gsGeometry<real_t> & g, const gsVector<real_t> & uv,
                       gsVector<real_t,3> & x,
                       gsVector<real_t,3> & xu, gsVector<real_t,3> & xv)
{
    gsMatrix<real_t> u(2,1), val, der;
    u.col(0) = uv;
    g.eval_into (u, val);   // 3 x 1
    g.deriv_into(u, der);   // 6 x 1, row 2*k+j = d x_k / d u_j
    for (short_t k = 0; k != 3; ++k)
    {
        x [k] = val(k,0);
        xu[k] = der(2*k + 0, 0);
        xv[k] = der(2*k + 1, 0);
    }
}

/// Induced boundary traversal tangent of \a side, for the orientation given by
/// n = x_u x x_v.  For the unit square with n = +z the boundary runs
/// counter-clockwise seen from +z: south +u, east +v, north -u, west -v.
static gsVector<real_t,3> edgeTangent(const gsGeometry<real_t> & g,
                                      const boxSide & side,
                                      const gsVector<real_t> & uv)
{
    gsVector<real_t,3> x, xu, xv;
    patchFrame(g, uv, x, xu, xv);
    const short_t d   = side.direction();     // fixed parametric direction
    const bool    par = side.parameter();     // false -> 0, true -> 1
    const real_t  sgn = (par ? (real_t)1 : (real_t)-1) * (d == 0 ? (real_t)1 : (real_t)-1);
    return sgn * (d == 0 ? xv : xu);          // derivative along the free direction
}

/// One side of one patch.
struct EdgeRef { index_t patch; boxSide side; };

/// Sample \a side of patch \a g at \a nSample points; \a pars is 2 x nSample
/// (patch parameters), \a pts is 3 x nSample (physical points).
static void edgeSamples(const gsGeometry<real_t> & g, const boxSide & side,
                        index_t nSample, gsMatrix<real_t> & pars, gsMatrix<real_t> & pts)
{
    const gsMatrix<real_t> supp = g.support();
    const short_t d    = side.direction();
    const short_t free = 1 - d;
    const real_t  fix  = side.parameter() ? supp(d,1) : supp(d,0);

    pars.resize(2, nSample);
    for (index_t i = 0; i != nSample; ++i)
    {
        const real_t t = (real_t)i / (real_t)(nSample - 1);
        pars(d,    i) = fix;
        pars(free, i) = supp(free,0) + t * (supp(free,1) - supp(free,0));
    }
    g.eval_into(pars, pts);
}

/// Symmetric (sampled) Hausdorff distance between two point sets.
static real_t hausdorff(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
{
    real_t h = 0;
    for (index_t i = 0; i != A.cols(); ++i)
    {
        real_t m = std::numeric_limits<real_t>::max();
        for (index_t j = 0; j != B.cols(); ++j)
            m = std::min(m, (B.col(j) - A.col(i)).squaredNorm());
        h = std::max(h, m);
    }
    for (index_t j = 0; j != B.cols(); ++j)
    {
        real_t m = std::numeric_limits<real_t>::max();
        for (index_t i = 0; i != A.cols(); ++i)
            m = std::min(m, (B.col(j) - A.col(i)).squaredNorm());
        h = std::max(h, m);
    }
    return math::sqrt(h);
}

/// Parameter of the point on \a side of patch \a g whose image is closest to
/// \a target, found by dense sampling of the edge.
static gsVector<real_t> matchOnEdge(const gsGeometry<real_t> & g,
                                    const boxSide & side,
                                    const gsVector<real_t,3> & target,
                                    index_t nSample = 400)
{
    gsMatrix<real_t> pars, pts;
    edgeSamples(g, side, nSample, pars, pts);

    index_t best = 0;
    real_t  bestD = std::numeric_limits<real_t>::max();
    for (index_t i = 0; i != nSample; ++i)
    {
        const real_t dd = (pts.col(i) - target).squaredNorm();
        if (dd < bestD) { bestD = dd; best = i; }
    }
    return pars.col(best);
}

/// Oriented surface integrals of the BREP:
///   \a fluxN  = oint n dS   (must vanish for a closed, consistently oriented shell)
///   \a volume = 1/3 oint x.n dS
///   \a area   = oint |n| dS
static void brepIntegrals(const gsMultiPatch<real_t> & brep,
                          const std::vector<real_t> & sgn,
                          gsVector<real_t,3> & fluxN, real_t & volume, real_t & area)
{
    fluxN.setZero(); volume = 0; area = 0;

    const index_t nSub = 4;
    gsVector<index_t> nnodes(2);
    nnodes << 10, 10;
    gsGaussRule<real_t> rule(nnodes);

    gsMatrix<real_t> pts, vals, ders;
    gsVector<real_t> wts;

    for (size_t p = 0; p != brep.nPatches(); ++p)
    {
        const gsGeometry<real_t> & g = brep.patch(p);
        for (auto & elem : g.basis().domain()->allElements())
        {
            const gsVector<real_t> elo = elem.lowerCorner();
            const gsVector<real_t> ehi = elem.upperCorner();
            for (index_t si = 0; si != nSub; ++si)
            for (index_t sj = 0; sj != nSub; ++sj)
            {
            gsVector<real_t> slo(2), shi(2);
            slo[0] = elo[0] + (ehi[0]-elo[0]) * (real_t)si     / (real_t)nSub;
            shi[0] = elo[0] + (ehi[0]-elo[0]) * (real_t)(si+1) / (real_t)nSub;
            slo[1] = elo[1] + (ehi[1]-elo[1]) * (real_t)sj     / (real_t)nSub;
            shi[1] = elo[1] + (ehi[1]-elo[1]) * (real_t)(sj+1) / (real_t)nSub;

            rule.mapTo(slo, shi, pts, wts);
            g.eval_into (pts, vals);
            g.deriv_into(pts, ders);

            for (index_t i = 0; i != pts.cols(); ++i)
            {
                gsVector<real_t,3> x, xu, xv;
                for (short_t k = 0; k != 3; ++k)
                {
                    x [k] = vals(k,i);
                    xu[k] = ders(2*k + 0, i);
                    xv[k] = ders(2*k + 1, i);
                }
                const gsVector<real_t,3> n = sgn[p] * xu.cross(xv);
                fluxN  += wts[i] * n;
                volume += wts[i] * x.dot(n) / (real_t)3;
                area   += wts[i] * n.norm();
            }
            }
        }
    }
}

/// Assign a per-patch orientation multiplier so that all patch normals point
/// consistently outward. A single-patch BRep (the sphere) needs no pairing:
/// the geometric side-pairing below only matches sides ACROSS patches (a
/// same-patch side is skipped), so it cannot itself resolve a one-patch
/// shell. Everything from the pairing onward is otherwise unchanged.
static std::vector<real_t> fixOrientation(const gsMultiPatch<real_t> & brep,
                                          real_t & refVolume, real_t & refArea,
                                          bool verbose = true)
{
    const size_t nP = brep.nPatches();
    std::vector<real_t> sgn(nP, 0.0);        // 0 = not yet visited

    if (nP == 1)
    {
        sgn[0] = 1.0;
    }
    else
    {
        // ---- pair up the patch sides GEOMETRICALLY ----------------------------
        const index_t nS = 40;
        std::vector<EdgeRef>          edges;
        std::vector<gsMatrix<real_t> > epts;
        for (size_t p = 0; p != nP; ++p)
            for (short_t s = 1; s <= 4; ++s)
            {
                EdgeRef er; er.patch = (index_t)p; er.side = boxSide(s);
                gsMatrix<real_t> pars, pts;
                edgeSamples(brep.patch(p), er.side, nS, pars, pts);
                edges.push_back(er);
                epts.push_back(pts);
            }

        const real_t diag = (brep.patch(0).coefs().colwise().maxCoeff()
                           - brep.patch(0).coefs().colwise().minCoeff()).norm();
        const real_t tol  = 1e-6 * std::max(diag, (real_t)1);

        const size_t nE = edges.size();
        std::vector<index_t> partner(nE, -1);
        for (size_t a = 0; a != nE; ++a)
        {
            real_t  bestH = std::numeric_limits<real_t>::max();
            index_t bestB = -1;
            for (size_t b = 0; b != nE; ++b)
            {
                if (edges[b].patch == edges[a].patch) continue;
                const real_t h = hausdorff(epts[a], epts[b]);
                if (h < bestH) { bestH = h; bestB = (index_t)b; }
            }
            if (bestH < tol) partner[a] = bestB;
        }
        for (size_t a = 0; a != nE; ++a)
        {
            GISMO_ENSURE(partner[a] >= 0,
                         "Side " << edges[a].side.index() << " of patch " << edges[a].patch
                         << " has no matching side on any other patch: the BREP is not "
                         "watertight, so the enclosed volume is undefined.");
            GISMO_ENSURE(partner[partner[a]] == (index_t)a,
                         "Side matching is not mutual for side " << edges[a].side.index()
                         << " of patch " << edges[a].patch << ".");
        }

        // ---- BFS over the geometric adjacency graph ---------------------------
        sgn[0] = 1.0;
        std::vector<size_t> queue(1, 0);
        for (size_t qi = 0; qi != queue.size(); ++qi)
        {
            const size_t p = queue[qi];
            for (size_t a = 0; a != nE; ++a)
            {
                if ((size_t)edges[a].patch != p) continue;
                const EdgeRef & self  = edges[a];
                const EdgeRef & other = edges[partner[a]];
                const size_t q = (size_t)other.patch;
                if (sgn[q] != 0.0) continue;     // already fixed

                const gsGeometry<real_t> & gS = brep.patch(p);
                const gsGeometry<real_t> & gO = brep.patch(q);

                gsMatrix<real_t> parsS, ptsS;
                edgeSamples(gS, self.side, 7, parsS, ptsS);

                real_t vote = 0;
                for (index_t i = 1; i + 1 < parsS.cols(); ++i)   // skip the corners
                {
                    const gsVector<real_t>   uvS = parsS.col(i);
                    const gsVector<real_t,3> xS  = ptsS.col(i);
                    const gsVector<real_t>   uvO = matchOnEdge(gO, other.side, xS);
                    const gsVector<real_t,3> tS  = edgeTangent(gS, self.side,  uvS);
                    const gsVector<real_t,3> tO  = edgeTangent(gO, other.side, uvO);

                    const real_t den = tS.norm() * tO.norm();
                    if (den > 0) vote += tS.dot(tO) / den;
                }

                sgn[q] = (vote < 0) ? sgn[p] : -sgn[p];
                queue.push_back(q);
            }
        }

        for (size_t p = 0; p != nP; ++p)
            GISMO_ENSURE(sgn[p] != 0.0,
                         "Patch " << p << " is not connected to patch 0: the BREP is "
                         "not a single connected shell.");
    }

    // ---- global validation -------------------------------------------------
    gsVector<real_t,3> fluxN;
    brepIntegrals(brep, sgn, fluxN, refVolume, refArea);

    if (refVolume < 0)                       // consistent but inward: flip all
    {
        for (size_t p = 0; p != nP; ++p) sgn[p] = -sgn[p];
        brepIntegrals(brep, sgn, fluxN, refVolume, refArea);
    }

    if (verbose)
    {
        gsInfo << "  patch orientation multipliers :";
        for (size_t p = 0; p != nP; ++p) gsInfo << (sgn[p] > 0 ? " +" : " -");
        gsInfo << "\n";
        gsInfo << "  |oint n dS| (closure)         : " << fluxN.norm() << "\n";
        gsInfo << "  surface area                  : "
               << std::setprecision(8) << refArea << "\n";
        gsInfo << "  volume  1/3 oint x.n dS       : "
               << std::setprecision(8) << refVolume << "\n";
    }

    GISMO_ENSURE(fluxN.norm() < 1e-8 * std::max(refArea, (real_t)1),
                 "BREP orientation fix failed: |oint n dS| = " << fluxN.norm()
                 << " (surface area " << refArea << "). The input is either not "
                 "watertight or the per-patch orientation could not be resolved.");
    GISMO_ENSURE(refVolume > 0, "Non-positive enclosed volume: " << refVolume);

    return sgn;
}

// =============================================================================
//  Geometry loaders
// =============================================================================

/// A loaded, oriented BRep plus the analytic data (rotation \a R, centre
/// \a c, cube side \a s or sphere radius \a r) needed by --selftest's
/// closed-form references. \a bg is the 3x2 background box (AABB of the
/// control net, padded by 5% of its own largest edge, EXCEPT for
/// cubeAligned, which fixes bg = [0,1]^3 so that every face lies on an
/// integer-fraction plane); \a L is its largest edge.
struct Geometry
{
    std::string name;
    gsMultiPatch<real_t> mp;
    std::vector<real_t> sgn;
    gsMatrix<real_t> bg;
    real_t L;
    gsMatrix<real_t> R;
    gsVector<real_t> c;
    real_t s, r;
};

static Geometry loadGeometry(const std::string & geom, const std::vector<real_t> & rot)
{
    Geometry g;
    g.name = geom;
    g.s = 0; g.r = 0;

    if (geom == "cube")
    {
        const real_t s = (real_t)0.6;
        gsVector<real_t> c(3); c << 0.5+0.0131, 0.5-0.0217, 0.5+0.0073;
        std::vector<real_t> rr = rot;
        if (rr.empty()) { rr.push_back(0.3); rr.push_back(0.5); rr.push_back(0.7); }
        GISMO_ENSURE(rr.size() == 3, "--rot needs exactly 3 values (phiz phiy phix), got " << rr.size() << ".");

        gsNurbsCreator<real_t>::TensorBSpline3Ptr cube =
            gsNurbsCreator<real_t>::BSplineCube(s, (real_t)-0.5, (real_t)-0.5, (real_t)-0.5);
        gsNurbsCreator<real_t>::rotate3D(*cube, rr[0], rr[1], rr[2]);
        cube->translate(c);

        for (short_t side = 1; side <= 6; ++side)
            g.mp.addPatch(cube->boundary(boxSide(side)));

        gsMatrix<real_t> Rz(3,3), Ry(3,3), Rx(3,3);
        Rz << math::cos(rr[0]), -math::sin(rr[0]), 0,
              math::sin(rr[0]),  math::cos(rr[0]), 0,
              0,                 0,                1;
        Ry << math::cos(rr[1]),  0, math::sin(rr[1]),
              0,                 1, 0,
             -math::sin(rr[1]),  0, math::cos(rr[1]);
        Rx << 1, 0,                 0,
              0, math::cos(rr[2]), -math::sin(rr[2]),
              0, math::sin(rr[2]),  math::cos(rr[2]);

        g.R = Rz*Ry*Rx;
        g.c = c;
        g.s = s;
    }
    else if (geom == "cubeAligned")
    {
        gsNurbsCreator<real_t>::TensorBSpline3Ptr cube =
            gsNurbsCreator<real_t>::BSplineCube((real_t)0.5, (real_t)0, (real_t)0, (real_t)0);
        for (short_t side = 1; side <= 6; ++side)
            g.mp.addPatch(cube->boundary(boxSide(side)));

        g.R = gsMatrix<real_t>::Identity(3,3);
        g.c.resize(3); g.c << (real_t)0.5, (real_t)0.5, (real_t)0.5;
        g.s = (real_t)0.5;
    }
    else if (geom == "sphere")
    {
        const real_t r = (real_t)0.4;
        gsVector<real_t> c(3); c << 0.5+0.0123, 0.5-0.0071, 0.5+0.0049;
        g.mp.addPatch(gsNurbsCreator<real_t>::NurbsSphere(r, c[0], c[1], c[2]));
        g.c = c;
        g.r = r;
    }
    else if (geom == "duck")
    {
        gsReadFile<real_t>("breps/3D/duck_BRep.xml", g.mp);
        GISMO_ENSURE(g.mp.nPatches() > 0, "Could not read a gsMultiPatch from 'breps/3D/duck_BRep.xml'.");
    }
    else
        GISMO_ENSURE(false, "Unknown --geom '" << geom << "' (expected cube, cubeAligned, sphere, duck).");

    real_t refVolume = 0, refArea = 0;
    g.sgn = fixOrientation(g.mp, refVolume, refArea, false);

    if (geom == "cubeAligned")
    {
        g.bg.resize(3,2);
        g.bg.col(0).setZero();
        g.bg.col(1).setOnes();
    }
    else
    {
        gsMatrix<real_t> box(3,2);
        box.col(0) = g.mp.patch(0).coefs().colwise().minCoeff().transpose();
        box.col(1) = g.mp.patch(0).coefs().colwise().maxCoeff().transpose();
        for (size_t p = 1; p != g.mp.nPatches(); ++p)
        {
            box.col(0) = box.col(0).cwiseMin(g.mp.patch(p).coefs().colwise().minCoeff().transpose());
            box.col(1) = box.col(1).cwiseMax(g.mp.patch(p).coefs().colwise().maxCoeff().transpose());
        }
        const real_t Lmax = (box.col(1)-box.col(0)).maxCoeff();
        const real_t pad  = (real_t)0.05*Lmax;
        box.col(0).array() -= pad;
        box.col(1).array() += pad;
        g.bg = box;
    }
    g.L = (g.bg.col(1)-g.bg.col(0)).maxCoeff();

    return g;
}

// =============================================================================
//  Bezier-element box: one knot-span element of one patch, together with the
//  AABB of its active control points (a convex-hull bound on the surface).
//  Used only to PRE-classify a background box via boxClassBRep, never as the
//  cut test itself (see the octree driver section below for why:
//  BezBox::overlaps is CLOSED, which would double-count a face lying
//  exactly on a box plane).
// =============================================================================
struct BezBox
{
    index_t patch;
    real_t  plo[2], pup[2];   // parameter box on the patch
    real_t  lo[3],  up[3];    // control-point AABB in physical space

    /// Squared distance from \a x to this box (0 inside).  This is a *lower
    /// bound* on the squared distance from x to the surface over this element,
    /// and is the only quantity that may be used for branch-and-bound pruning.
    real_t sqDistLowerBound(const real_t * x) const
    {
        real_t d2 = 0;
        for (short_t k = 0; k != 3; ++k)
        {
            const real_t e = (x[k] < lo[k]) ? lo[k] - x[k]
                           : (x[k] > up[k]) ? x[k] - up[k] : (real_t)0;
            d2 += e * e;
        }
        return d2;
    }

    bool overlaps(const gsVector<real_t> & clo, const gsVector<real_t> & chi) const
    {
        return up[0] >= clo[0] && lo[0] <= chi[0]
            && up[1] >= clo[1] && lo[1] <= chi[1]
            && up[2] >= clo[2] && lo[2] <= chi[2];
    }
};

/// Build the table of Bezier-element boxes for a whole BREP.
static std::vector<BezBox> makeBezBoxes(const gsMultiPatch<real_t> & brep)
{
    std::vector<BezBox> boxes;
    gsMatrix<real_t>  centre(2,1);
    gsMatrix<index_t> act;

    for (size_t p = 0; p != brep.nPatches(); ++p)
    {
        const gsGeometry<real_t> & g = brep.patch(p);
        const gsBasis<real_t>    & b = g.basis();
        const gsMatrix<real_t>   & C = g.coefs();

        for (auto & elem : b.domain()->allElements())
        {
            const gsVector<real_t> lo = elem.lowerCorner();
            const gsVector<real_t> hi = elem.upperCorner();
            centre.col(0) = (real_t)0.5 * (lo + hi);

            // Active (= non-vanishing) basis functions on this element; the
            // surface over the element lies in the hull of their coefficients.
            b.active_into(centre, act);

            BezBox bb;
            bb.patch = (index_t)p;
            for (short_t d = 0; d != 2; ++d) { bb.plo[d] = lo[d]; bb.pup[d] = hi[d]; }
            for (short_t k = 0; k != 3; ++k)
            {
                bb.lo[k] =  std::numeric_limits<real_t>::max();
                bb.up[k] = -std::numeric_limits<real_t>::max();
            }
            for (index_t i = 0; i != act.rows(); ++i)
                for (short_t k = 0; k != 3; ++k)
                {
                    const real_t c = C(act(i,0), k);
                    bb.lo[k] = std::min(bb.lo[k], c);
                    bb.up[k] = std::max(bb.up[k], c);
                }
            boxes.push_back(bb);
        }
    }
    return boxes;
}

// =============================================================================
//  gsBRepSignedDist<T> : signed distance to a spline BREP
//
//  Convention (matches gsCutCellRule and the other immersed examples):
//    phi < 0  -> inside      phi > 0  -> outside
//
//  eval_into/deriv_into evaluate one query point at a time, with no
//  `#pragma omp`: the class itself is thread-safe per query (it never
//  calls into gsBRepRayCast.h), but this driver is serial throughout
//  (see the file doxygen).
// =============================================================================
template<class T>
class gsBRepSignedDist : public gsFunction<T>
{
public:
    GISMO_CLONE_FUNCTION(gsBRepSignedDist)

    gsBRepSignedDist(const gsMultiPatch<T> & brep,
                     std::vector<BezBox>     boxes,
                     std::vector<T>          sgn,
                     const gsMatrix<T> &     bbox,
                     index_t maxIter = 24, T tol = (T)1e-12)
    : m_brep(&brep), m_boxes(give(boxes)), m_sgn(give(sgn)), m_bbox(bbox),
      m_maxIter(maxIter), m_tol(tol)
    {
        GISMO_ENSURE(!m_boxes.empty(),
                     "gsBRepSignedDist needs at least one Bezier-element box.");
    }

    short_t     domainDim() const override { return 3; }
    short_t     targetDim() const override { return 1; }
    gsMatrix<T> support()   const override { return m_bbox; }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k < u.cols(); ++k)
        {
            T phi; gsVector<T,3> grad;
            evalOne(u.col(k), phi, grad, false);
            result(0,k) = phi;
        }
    }

    /// Analytic gradient: grad(phi) = sign * (x - c)/|x - c|, free once the
    /// footpoint c is known.  (gsMeshSignedDist falls back to central
    /// differences here, which costs 6 extra footpoint solves per gradient.)
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(3, u.cols());
        for (index_t k = 0; k < u.cols(); ++k)
        {
            T phi; gsVector<T,3> grad;
            evalOne(u.col(k), phi, grad, true);
            result.col(k) = grad;
        }
    }

private:
    /// Footpoint solve + sign for a single query point.
    void evalOne(const gsVector<T> & xq, T & phi, gsVector<T,3> & grad,
                 bool wantGrad) const
    {
        const T x[3] = { xq[0], xq[1], xq[2] };

        // --- branch and bound over the Bezier-element boxes -----------------
        // One O(nBoxes) pass for the lower bounds (no sort: this runs on every
        // single phi evaluation, and algoim asks for a great many), then solve
        // the nearest box to prime the bound and only visit boxes that can
        // still beat it.
        const size_t nB = m_boxes.size();
        static thread_local std::vector<T> lb;       // reused across calls
        lb.resize(nB);

        size_t iMin = 0;
        T lbMin = std::numeric_limits<T>::max();
        for (size_t b = 0; b != nB; ++b)
        {
            lb[b] = m_boxes[b].sqDistLowerBound(x);
            if (lb[b] < lbMin) { lbMin = lb[b]; iMin = b; }
        }

        gsVector<T,3> bestC, bestN, c, n;
        T bestD2 = footpoint(m_boxes[iMin], xq, bestC, bestN);

        for (size_t b = 0; b != nB; ++b)
        {
            // Prune ONLY on the control-hull lower bound.  Pruning on anything
            // derived from a box-clamped Newton result would be unsound: a
            // neighbouring element's constrained minimum can undercut the true
            // global distance and would skip the box holding the real footpoint.
            if (b == iMin || lb[b] >= bestD2) continue;

            const T d2 = footpoint(m_boxes[b], xq, c, n);
            if (d2 < bestD2) { bestD2 = d2; bestC = c; bestN = n; }
        }

        const T dist = math::sqrt(bestD2);
        const gsVector<T,3> w = xq - bestC;
        const T s = (w.dot(bestN) < 0) ? (T)-1 : (T)1;

        phi = s * dist;
        if (wantGrad)
            grad = (dist > (T)1e-14) ? (s / dist) * w
                                     : bestN.normalized();   // on the surface
    }

    /// Gauss-Newton minimisation of 1/2 |S(u) - x|^2 over the parameter box of
    /// \a bb, seeded at its centre and projected back into the box each step.
    /// Returns the squared distance; \a c is the footpoint and \a n the
    /// oriented normal there.
    T footpoint(const BezBox & bb, const gsVector<T> & xq,
                gsVector<T,3> & c, gsVector<T,3> & n) const
    {
        const gsGeometry<T> & g = m_brep->patch(bb.patch);

        gsVector<T> uv(2);
        uv << (T)0.5*(bb.plo[0] + bb.pup[0]), (T)0.5*(bb.plo[1] + bb.pup[1]);

        gsVector<T,3> x, xu, xv, r;
        gsMatrix<T,2,2> A;
        gsVector<T,2> rhs, du;

        patchFrame(g, uv, x, xu, xv);
        r = x - xq;
        T f = r.squaredNorm();

        T mu = (T)1e-8;                       // Levenberg damping
        for (index_t it = 0; it != m_maxIter; ++it)
        {
            A(0,0) = xu.dot(xu); A(0,1) = xu.dot(xv);
            A(1,0) = A(0,1);     A(1,1) = xv.dot(xv);
            A(0,0) += mu; A(1,1) += mu;
            rhs << -xu.dot(r), -xv.dot(r);

            const T det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
            if (math::abs(det) < (T)1e-30) break;
            du[0] = ( A(1,1)*rhs[0] - A(0,1)*rhs[1]) / det;
            du[1] = (-A(1,0)*rhs[0] + A(0,0)*rhs[1]) / det;

            // Projected step with simple backtracking.
            bool improved = false;
            T step = (T)1;
            for (index_t ls = 0; ls != 8; ++ls, step *= (T)0.5)
            {
                gsVector<T> uvT(2);
                for (short_t d = 0; d != 2; ++d)
                    uvT[d] = std::min(std::max(uv[d] + step*du[d],
                                               (T)bb.plo[d]), (T)bb.pup[d]);
                if ((uvT - uv).norm() < m_tol) break;

                gsVector<T,3> xT, xuT, xvT;
                patchFrame(g, uvT, xT, xuT, xvT);
                const T fT = (xT - xq).squaredNorm();
                if (fT < f)
                {
                    uv = uvT; x = xT; xu = xuT; xv = xvT;
                    r = x - xq; f = fT; improved = true;
                    mu = std::max(mu * (T)0.5, (T)1e-12);
                    break;
                }
            }
            if (!improved) { mu *= (T)10; if (mu > (T)1e6) break; }
        }

        c = x;
        n = m_sgn[bb.patch] * xu.cross(xv);
        return f;
    }

    const gsMultiPatch<T> * m_brep;
    std::vector<BezBox>     m_boxes;
    std::vector<T>          m_sgn;
    gsMatrix<T>             m_bbox;
    index_t                 m_maxIter;
    T                       m_tol;
};

// =============================================================================
//  Cell classification (counterpart of boxClassSAT in the mesh example)
//    -1 = fully inside, 0 = cut, +1 = fully outside
//
//  A NECESSARY pre-classifier only (see its call site in the octree driver):
//  a control-hull "hit" only means "possibly cut" (a face lying exactly on a
//  box plane also hits the neighbouring box), so the actual cut test used to
//  drive the octree is the exact, half-open BezElement::lo/up overlap
//  (::overlapR3 below), never this function's BezBox::overlaps.
// =============================================================================
static int boxClassBRep(const std::vector<BezBox> & boxes,
                        const std::vector<int>    & cIdx,
                        const gsFunction<real_t>  & phi,
                        const gsVector<real_t>    & lo,
                        const gsVector<real_t>    & hi)
{
    gsMatrix<real_t> ctr(3,1), val;
    ctr.col(0) = (real_t)0.5 * (lo + hi);
    phi.eval_into(ctr, val);
    const real_t p = val(0,0);

    bool hit = false;
    for (size_t i = 0; i != cIdx.size(); ++i)
        if (boxes[cIdx[i]].overlaps(lo, hi)) { hit = true; break; }

    if (!hit) return (p < 0) ? -1 : 1;

    // A control-hull overlap is a *necessary* condition only, so "hit" alone
    // means "possibly cut".  A true signed distance is 1-Lipschitz, hence any
    // point of the cell is within the half-diagonal of the centre: if
    // |phi(centre)| exceeds that radius the whole cell is strictly on one side.
    // (Half-diagonal, not half-edge -- a half-edge radius makes "confidently
    // inside" wrong near the surface, which is a silent volume error.)
    const real_t radius = (real_t)0.5 * (hi - lo).norm();
    if (p < -radius) return -1;
    if (p >  radius) return  1;
    return 0;
}

// =============================================================================
//  Exact volume-moment oracle (divergence theorem on the Bezier decomposition)
// =============================================================================

/**
    @brief Integral_Omega x^a y^b z^c dV via the divergence theorem on the
    Bezier elements \a els:

        Integral_Omega x^a y^b z^c dV
            = Sum_elements Integral_[uv box] F(S(u,v)) . (sgn.(S_u x S_v))_x du dv,
              F(x,y,z) = x^(a+1)/(a+1) . y^b . z^c,

    since div(F,0,0) = x^a y^b z^c. sgn.(S_u x S_v) is read off as
    brc::normalNumerator(hom,sgn,uv)/w^4, so the integrand is evaluated
    without ever dividing by w except by that single w^4 (the homogeneous
    net's own weight component, never zero).

    Gauss node count per direction, n = a+b+c, p,q the element's bidegree:
    - polynomial elements: m_u = ((n+3)p+1)/2, m_v = ((n+3)q+1)/2 (integer
      division). This is exact because F(S(u,v)) has bidegree ((n+1)p,(n+1)q)
      and (S_u x S_v)_x has bidegree (2p-1,2q-1), so the integrand has
      bidegree ((n+3)p-1,(n+3)q-1), exactly integrated by a Gauss rule of
      ceil(((n+3)p-1+1)/2) = ((n+3)p+1)/2 (integer division) nodes.
    - rational elements: the same count plus \a extraRational per
      direction, since the integrand is then a genuine rational function
      (no closed bidegree bound); the analytic comparison in --selftest is
      what decides whether the count is adequate -- raise it on failure,
      never loosen a tolerance.
*/
static real_t momentOracle(const std::vector<brc::BezElement<real_t> > & els,
                           int a, int b, int c, int extraRational = 20)
{
    const int n = a+b+c;
    real_t total = 0;

    for (std::vector<brc::BezElement<real_t> >::const_iterator el = els.begin(); el != els.end(); ++el)
    {
        const index_t p = el->hom.basis().degree(0), q = el->hom.basis().degree(1);
        index_t mu = ((n+3)*p+1)/2;
        index_t mv = ((n+3)*q+1)/2;
        if (el->rational) { mu += extraRational; mv += extraRational; }

        gsVector<index_t> nnodes(2); nnodes << mu, mv;
        gsGaussRule<real_t> rule(nnodes);

        gsMatrix<real_t> lo(2,1), hi(2,1);
        lo << el->uv[0][0], el->uv[1][0];
        hi << el->uv[0][1], el->uv[1][1];

        gsMatrix<real_t> pts; gsVector<real_t> wts;
        rule.mapTo(lo, hi, pts, wts);

        gsMatrix<real_t> N;
        brc::normalNumerator(el->hom, el->sgn, pts, N); // 3 x nPts = sgn.w^4.(S_u x S_v)

        gsMatrix<real_t> val;
        el->hom.eval_into(pts, val); // 4 x nPts

        for (index_t i = 0; i != pts.cols(); ++i)
        {
            const real_t w = val(3,i);
            const real_t x = val(0,i)/w, y = val(1,i)/w, z = val(2,i)/w;
            const real_t F = std::pow(x, a+1)/(real_t)(a+1) * std::pow(y, (real_t)b) * std::pow(z, (real_t)c);
            total += wts[i] * F * N(0,i) / (w*w*w*w);
        }
    }
    return total;
}

// =============================================================================
//  Analytic references (self-test only)
// =============================================================================

/// Exact Integral x^a y^b z^c over Omega = c + R.[-s/2,s/2]^3, via an
/// n/2+1-node-per-direction 3D Gauss rule on [-s/2,s/2]^3 mapped by
/// x = c + R.eta (a rotation, so dV = d(eta) exactly).
static real_t cubeMoment(const gsVector<real_t> & c, const gsMatrix<real_t> & R, real_t s,
                         int a, int b, int cExp)
{
    const int n = a+b+cExp;
    const index_t nodes = n/2 + 1;
    gsVector<index_t> nn(3); nn << nodes, nodes, nodes;
    gsGaussRule<real_t> rule(nn);

    gsVector<real_t> lo(3), hi(3);
    lo << -s/2, -s/2, -s/2;
    hi <<  s/2,  s/2,  s/2;
    gsMatrix<real_t> pts; gsVector<real_t> wts;
    rule.mapTo(lo, hi, pts, wts);

    real_t total = 0;
    for (index_t i = 0; i != pts.cols(); ++i)
    {
        const gsVector<real_t> x = c + R*pts.col(i);
        total += wts[i] * std::pow(x[0],a) * std::pow(x[1],b) * std::pow(x[2],cExp);
    }
    return total;
}

/// Binomial coefficient C(n,k) for small nonnegative integers, exact in
/// double precision (each step is an integer division).
static real_t binomCoeff(int n, int k)
{
    real_t v = 1;
    for (int t = 1; t <= k; ++t) v = v*(real_t)(n-k+t)/(real_t)t;
    return v;
}

/// Centred moment M(i,j,l) = Integral X^i Y^j Z^l dV over a ball of radius
/// \a r at the origin: 0 unless i,j,l are all even, otherwise
/// M = 2.Gamma((i+1)/2).Gamma((j+1)/2).Gamma((l+1)/2)/Gamma((i+j+l+3)/2)
///     . r^(i+j+l+3)/(i+j+l+3).
static real_t ballCenteredMoment(int i, int j, int l, real_t r)
{
    if (i%2 != 0 || j%2 != 0 || l%2 != 0) return (real_t)0;
    const real_t num = 2.0 * std::tgamma((i+1)/2.0) * std::tgamma((j+1)/2.0) * std::tgamma((l+1)/2.0);
    const real_t den = std::tgamma((i+j+l+3)/2.0);
    const int nTot = i+j+l+3;
    return num/den * std::pow(r, (real_t)nTot) / (real_t)nTot;
}

/// Integral (c_x+X)^a (c_y+Y)^b (c_z+Z)^c dV over a ball of radius \a r
/// centred at \a c, by binomial expansion about the origin-centred moments.
static real_t sphereMoment(const gsVector<real_t> & c, real_t r, int a, int b, int cExp)
{
    real_t total = 0;
    for (int i = 0; i <= a; ++i)
    for (int j = 0; j <= b; ++j)
    for (int l = 0; l <= cExp; ++l)
    {
        const real_t coef = binomCoeff(a,i)*binomCoeff(b,j)*binomCoeff(cExp,l);
        const real_t powc = std::pow(c[0],a-i)*std::pow(c[1],b-j)*std::pow(c[2],cExp-l);
        total += coef*powc*ballCenteredMoment(i,j,l,r);
    }
    return total;
}

// =============================================================================
//  --selftest
// =============================================================================

/// Running state of the self-test battery: the failure count decides the
/// process exit code, the RNG is shared (and re-seeded once, fixed) across
/// every check so the whole run is reproducible.
struct SelfTestState { int nFail; std::mt19937 rng; };

/// Prints one "PASS <label>  max=<r> tol=1" (or FAIL) line and updates
/// \a st.nFail. \a worstRatio is the worst observed (error / that check's
/// own tolerance) ratio over every sample the check took: the check passes
/// iff worstRatio <= 1, which is what "tol=1" always means here (see the
/// file doxygen for why every check reports on this common scale).
static void report(SelfTestState & st, const std::string & label, real_t worstRatio,
                   const std::string & extra = "")
{
    const bool pass = worstRatio <= (real_t)1.0;
    gsInfo << (pass ? "PASS " : "FAIL ") << label
           << "  max=" << fmtSci(worstRatio) << " tol=" << fmtSci((real_t)1.0);
    if (!extra.empty()) gsInfo << "  " << extra;
    gsInfo << "\n";
    if (!pass) ++st.nFail;
}

static void checkSphereStructure(SelfTestState & st, const Geometry & sph)
{
    const gsTensorNurbs<2,real_t> * nurbs = dynamic_cast<const gsTensorNurbs<2,real_t>*>(&sph.mp.patch(0));
    GISMO_ENSURE(nurbs, "sphere-structure: patch 0 is not a gsTensorNurbs<2,real_t>.");

    bool ok = true;
    real_t maxErr = 0;
    const real_t scale = sph.c.cwiseAbs().maxCoeff() + sph.r;

    ok = ok && (nurbs->basis().source().degree(0) == 2) && (nurbs->basis().source().degree(1) == 2);

    // Exact knot vectors, per gsNurbsCreator<T>::NurbsSphere: u = {0,0,0,.5,.5,1,1,1},
    // v = {0,0,0,.25,.25,.5,.5,.75,.75,1,1,1} (both degree 2). Knot values are exact
    // dyadic rationals in double precision, so this is a bitwise comparison, not a
    // tolerance check.
    const gsKnotVector<real_t> expectedU((real_t)0,(real_t)1,1,3,2);
    const gsKnotVector<real_t> expectedV((real_t)0,(real_t)1,3,3,2);
    ok = ok && (nurbs->basis().source().knots(0) == expectedU);
    ok = ok && (nurbs->basis().source().knots(1) == expectedV);

    const gsMatrix<real_t> & C = nurbs->coefs();
    const gsMatrix<real_t> & W = nurbs->weights();
    ok = ok && (C.rows() == 45);

    for (index_t i = 0; i != 5 && ok; ++i)
    {
        const index_t r0 = i+5*0, r8 = i+5*8;
        if (!(C.row(r0) == C.row(r8)) || W(r0,0) != W(r8,0)) ok = false;
    }

    gsVector<real_t,3> north; north << sph.c[0], sph.c[1], sph.c[2]+sph.r;
    gsVector<real_t,3> south; south << sph.c[0], sph.c[1], sph.c[2]-sph.r;
    for (index_t i = 0; i != 5; ++i)
    {
        const index_t r2 = i+5*2, r6 = i+5*6;
        const real_t e2 = (C.row(r2).transpose()-north).norm();
        const real_t e6 = (C.row(r6).transpose()-south).norm();
        maxErr = std::max(maxErr, std::max(e2,e6)/((real_t)1e-13*scale));
        if (e2 > (real_t)1e-13*scale || e6 > (real_t)1e-13*scale) ok = false;
    }

    for (int t = 0; t != 20; ++t)
    {
        const real_t v = (real_t)0.25*(real_t)t/(real_t)19;
        gsMatrix<real_t> uv1(2,1); uv1 << 0, v;
        gsMatrix<real_t> uv2(2,1); uv2 << 0, (real_t)0.5-v;
        gsMatrix<real_t> p1, p2;
        nurbs->eval_into(uv1,p1); nurbs->eval_into(uv2,p2);
        const real_t e = (p1-p2).norm();
        maxErr = std::max(maxErr, e/((real_t)1e-14*scale));
        if (e > (real_t)1e-14*scale) ok = false;
    }

    std::vector<brc::BezElement<real_t> > els = brc::extractElements(sph.mp, sph.sgn);
    ok = ok && (els.size() == 8);
    for (std::vector<brc::BezElement<real_t> >::const_iterator e = els.begin(); e != els.end(); ++e)
        ok = ok && e->rational;

    report(st, "sphere-structure", ok ? maxErr : (real_t)1e300);
}

static void checkSphereClosure(SelfTestState & st, const Geometry & sph, std::mt19937 & rng)
{
    real_t worst = 0;
    std::uniform_real_distribution<real_t> u01(0.0, 1.0);
    for (int i = 0; i != 1000; ++i)
    {
        gsMatrix<real_t> uv(2,1); uv << u01(rng), u01(rng);
        gsMatrix<real_t> S;
        sph.mp.patch(0).eval_into(uv, S);
        const real_t d = (S.col(0)-sph.c).norm();
        const real_t e = std::abs(d-sph.r);
        worst = std::max(worst, e/((real_t)1e-14*sph.r));
    }

    gsVector<real_t,3> fluxN; real_t vol, area;
    brepIntegrals(sph.mp, sph.sgn, fluxN, vol, area);
    worst = std::max(worst, fluxN.norm()/((real_t)1e-12*area));

    report(st, "sphere-closure", worst);
}

static void checkNormalNumerator(SelfTestState & st, const std::string & label, const Geometry & geo,
                                 const std::vector<brc::BezElement<real_t> > & els, std::mt19937 & rng)
{
    real_t worst = 0;
    std::uniform_real_distribution<real_t> u01(0.0, 1.0);
    for (std::vector<brc::BezElement<real_t> >::const_iterator el = els.begin(); el != els.end(); ++el)
    {
        gsTensorBSplineBasis<2,real_t> Bref(gsKnotVector<real_t>(0,1,0,el->dU+1),
                                            gsKnotVector<real_t>(0,1,0,el->dV+1));
        gsTensorBSpline<2,real_t> Nspline(Bref, el->Ncoef);

        for (int s = 0; s != 20; ++s)
        {
            const real_t ru = u01(rng), rv = u01(rng);
            const real_t uu = el->uv[0][0] + ru*(el->uv[0][1]-el->uv[0][0]);
            const real_t vv = el->uv[1][0] + rv*(el->uv[1][1]-el->uv[1][0]);

            gsMatrix<real_t> refuv(2,1); refuv << ru, rv;
            gsMatrix<real_t> Nb;
            Nspline.eval_into(refuv, Nb);

            gsMatrix<real_t> uv(2,1); uv << uu, vv;
            gsMatrix<real_t> der;
            geo.mp.patch(el->patch).deriv_into(uv, der);
            gsVector<real_t,3> Su, Sv;
            for (short_t k = 0; k != 3; ++k) { Su[k] = der(2*k,0); Sv[k] = der(2*k+1,0); }

            gsMatrix<real_t> hval;
            el->hom.eval_into(uv, hval);
            const real_t w = hval(3,0);
            const gsVector<real_t,3> Nexact = el->sgn * std::pow(w,4) * Su.cross(Sv);

            const real_t e = (Nb.col(0)-Nexact).cwiseAbs().maxCoeff();
            const real_t tol = (real_t)1e-12*el->Nscale;
            worst = std::max(worst, e/tol);
        }
    }
    report(st, label, worst);
}

static void checkSignBounds(SelfTestState & st, const std::string & label,
                            const std::vector<brc::BezElement<real_t> > & els, std::mt19937 & rng)
{
    real_t worst = 0;
    std::uniform_real_distribution<real_t> u01(0.0, 1.0);
    int weak[3] = {0,0,0}, strict[3] = {0,0,0}, uncert[3] = {0,0,0};

    for (std::vector<brc::BezElement<real_t> >::const_iterator el = els.begin(); el != els.end(); ++el)
    {
        for (short_t k = 0; k != 3; ++k)
        {
            if      (el->certStrict[k] != 0) ++strict[k];
            else if (el->certWeak[k]   != 0) ++weak[k];
            else                             ++uncert[k];
        }

        for (int s = 0; s != 20; ++s)
        {
            const real_t ru = u01(rng), rv = u01(rng);
            const real_t uu = el->uv[0][0] + ru*(el->uv[0][1]-el->uv[0][0]);
            const real_t vv = el->uv[1][0] + rv*(el->uv[1][1]-el->uv[1][0]);
            gsMatrix<real_t> uv(2,1); uv << uu, vv;
            gsMatrix<real_t> N;
            brc::normalNumerator(el->hom, el->sgn, uv, N);

            for (short_t k = 0; k != 3; ++k)
            {
                const real_t tol = (real_t)1e-13*el->Nscale;
                const real_t lo = el->Nmin[k]-tol, hi = el->Nmax[k]+tol;
                const real_t v = N(k,0);
                real_t violation = 0;
                if (v < lo) violation = lo-v;
                if (v > hi) violation = std::max(violation, v-hi);
                worst = std::max(worst, violation/std::max(tol, (real_t)1e-300));
            }
        }
    }

    std::ostringstream extra;
    extra << "weak=(" << weak[0] << "," << weak[1] << "," << weak[2] << ")"
          << " strict=(" << strict[0] << "," << strict[1] << "," << strict[2] << ")"
          << " uncertified=(" << uncert[0] << "," << uncert[1] << "," << uncert[2] << ")";
    report(st, label, worst, extra.str());
}

static void checkHitsAndWinding(SelfTestState & st, const std::string & gname, const Geometry & geo,
                                const std::vector<brc::BezElement<real_t> > & els,
                                index_t nlines, std::mt19937 & rng)
{
    real_t worstResidual = 0;
    bool allCastOK = true;
    bool windingAllOK = true;

    for (short_t k = 0; k != 3; ++k)
    {
        const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);
        std::uniform_real_distribution<real_t> di(geo.bg(i,0), geo.bg(i,1));
        std::uniform_real_distribution<real_t> dj(geo.bg(j,0), geo.bg(j,1));

        for (index_t n = 0; n < nlines; ++n)
        {
            const real_t yi = di(rng), yj = dj(rng);
            std::vector<brc::Hit<real_t> > hits;
            const bool castOK = brc::castLine(els, k, yi, yj, hits, geo.L);
            if (!castOK) allCastOK = false;

            for (std::vector<brc::Hit<real_t> >::const_iterator h = hits.begin(); h != hits.end(); ++h)
            {
                gsMatrix<real_t> uv(2,1); uv << h->u, h->v;
                gsMatrix<real_t> val;
                geo.mp.patch(h->patch).eval_into(uv, val);
                const real_t res = std::max(std::abs(val(i,0)-yi), std::abs(val(j,0)-yj));
                worstResidual = std::max(worstResidual, res/((real_t)1e-13*geo.L));
            }

            if (!brc::windingOK(hits, k)) windingAllOK = false;
            if (hits.size() % 2 != 0) windingAllOK = false;
        }
    }

    report(st, gname+"-hits-on-surface", allCastOK ? worstResidual : (real_t)1e300);
    report(st, gname+"-winding", (allCastOK && windingAllOK) ? (real_t)0 : (real_t)1e300);
}

static void cubeSlabRoots(const Geometry & geo, short_t k, real_t yi, real_t yj, std::vector<real_t> & roots)
{
    roots.clear();
    const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);

    gsVector<real_t,3> P0; P0.setZero();
    P0[i] = yi; P0[j] = yj; P0[k] = 0;
    const gsVector<real_t,3> a = geo.R.transpose()*(P0-geo.c);
    gsVector<real_t,3> ek; ek.setZero(); ek[k] = 1;
    const gsVector<real_t,3> d = geo.R.transpose()*ek;

    real_t tEntry = -std::numeric_limits<real_t>::max();
    real_t tExit  =  std::numeric_limits<real_t>::max();
    bool valid = true;
    for (short_t m = 0; m != 3; ++m)
    {
        const real_t lo = -geo.s/2, hi = geo.s/2;
        if (d[m] == 0)
        {
            if (a[m] < lo || a[m] > hi) { valid = false; break; }
        }
        else
        {
            real_t t0 = (lo-a[m])/d[m];
            real_t t1 = (hi-a[m])/d[m];
            if (t0 > t1) std::swap(t0,t1);
            tEntry = std::max(tEntry,t0);
            tExit  = std::min(tExit,t1);
        }
    }
    if (valid && tEntry <= tExit)
    {
        roots.push_back(tEntry);
        roots.push_back(tExit);
    }
}

static void sphereRoots(const Geometry & geo, short_t k, real_t yi, real_t yj, std::vector<real_t> & roots)
{
    roots.clear();
    const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);
    const real_t d2 = geo.r*geo.r - (yi-geo.c[i])*(yi-geo.c[i]) - (yj-geo.c[j])*(yj-geo.c[j]);
    if (d2 < 0) return;
    const real_t d = std::sqrt(d2);
    roots.push_back(geo.c[k]-d);
    roots.push_back(geo.c[k]+d);
}

static void checkAnalyticRoots(SelfTestState & st, const std::string & gname, const Geometry & geo,
                               const std::vector<brc::BezElement<real_t> > & els,
                               index_t nlines, std::mt19937 & rng)
{
    real_t worst = 0;
    bool ok = true;
    index_t nExcluded = 0;

    for (short_t k = 0; k != 3; ++k)
    {
        const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);
        std::uniform_real_distribution<real_t> di(geo.bg(i,0), geo.bg(i,1));
        std::uniform_real_distribution<real_t> dj(geo.bg(j,0), geo.bg(j,1));

        for (index_t n = 0; n < nlines; ++n)
        {
            const real_t yi = di(rng), yj = dj(rng);
            std::vector<brc::Hit<real_t> > hits;
            brc::castLine(els, k, yi, yj, hits, geo.L);

            std::vector<real_t> roots;
            if (gname == "sphere") sphereRoots(geo, k, yi, yj, roots);
            else                   cubeSlabRoots(geo, k, yi, yj, roots);
            std::sort(roots.begin(), roots.end());

            bool nearTangent = false;
            for (std::vector<brc::Hit<real_t> >::const_iterator h = hits.begin(); h != hits.end(); ++h)
                if (std::abs(h->n[k]) < (real_t)1e-2) nearTangent = true;
            if (nearTangent) { ++nExcluded; continue; }

            if (hits.size() != roots.size()) { ok = false; continue; }
            for (size_t m = 0; m != hits.size(); ++m)
                worst = std::max(worst, std::abs(hits[m].t-roots[m])/((real_t)1e-13*geo.L));
        }
    }

    std::ostringstream extra; extra << "excluded=" << nExcluded;
    report(st, gname+"-analytic-roots", ok ? worst : (real_t)1e300, extra.str());
}

static void checkMoments(SelfTestState & st, const std::string & gname, const Geometry & geo,
                         const std::vector<brc::BezElement<real_t> > & els)
{
    real_t worst = 0;
    const real_t rho = geo.bg.cwiseAbs().maxCoeff();
    const real_t V = (gname == "sphere") ? ((real_t)4.0/(real_t)3.0)*EIGEN_PI*geo.r*geo.r*geo.r
                                          : geo.s*geo.s*geo.s;

    for (int a = 0; a <= 4; ++a)
    for (int b = 0; a+b <= 4; ++b)
    for (int c = 0; a+b+c <= 4; ++c)
    {
        const real_t oracle = momentOracle(els, a, b, c);
        const real_t ref = (gname == "sphere") ? sphereMoment(geo.c, geo.r, a, b, c)
                                                : cubeMoment(geo.c, geo.R, geo.s, a, b, c);
        const int n = a+b+c;
        const real_t err = std::abs(oracle-ref)/(V*std::pow(rho,(real_t)n));
        worst = std::max(worst, err/(real_t)1e-14);
    }
    report(st, gname+"-moments", worst);
}

static void checkDuckVolume(SelfTestState & st, const Geometry & duck, const std::vector<brc::BezElement<real_t> > & els)
{
    const real_t V = momentOracle(els, 0, 0, 0);
    const real_t ratioAnchor = std::abs(V-(real_t)1.1905049)/(real_t)5e-8;

    gsVector<real_t,3> fluxN; real_t refVolume, refArea;
    brepIntegrals(duck.mp, duck.sgn, fluxN, refVolume, refArea);
    const real_t ratioRel = (std::abs(V-refVolume)/std::abs(refVolume))/(real_t)1e-13;

    std::ostringstream extra;
    extra << "V=" << std::setprecision(10) << V;
    report(st, "duck-volume", std::max(ratioAnchor,ratioRel), extra.str());
}

static void checkDuckOrientation(SelfTestState & st, const Geometry & duck)
{
    const int expected[6] = { 1,-1,-1, 1, 1,-1 };
    bool ok = (duck.sgn.size() == 6);
    for (size_t p = 0; p != duck.sgn.size() && ok; ++p)
        if ((duck.sgn[p] > 0 ? 1 : -1) != expected[p]) ok = false;
    report(st, "duck-orientation", ok ? (real_t)0 : (real_t)1e300);
}

static void checkDuckMerge(SelfTestState & st, const Geometry & duck,
                          const std::vector<brc::BezElement<real_t> > & els, std::mt19937 & rng)
{
    real_t worst = 0;
    bool ok = true;
    std::uniform_real_distribution<real_t> u01(0.0, 1.0);
    const index_t nLines = 20;
    const index_t nP = (index_t)duck.mp.nPatches();

    for (index_t n = 0; n < nLines; ++n)
    {
        const index_t p = n % nP;
        const gsGeometry<real_t> & patch = duck.mp.patch(p);
        const gsTensorBSpline<2,real_t> * tb = dynamic_cast<const gsTensorBSpline<2,real_t>*>(&patch);
        GISMO_ENSURE(tb, "duck-merge: expected a gsTensorBSpline<2,real_t> patch.");

        const gsKnotVector<real_t> & ku = tb->basis().knots(0);
        gsKnotVector<real_t>::uiterator itBegin = ku.ubegin()+1;
        gsKnotVector<real_t>::uiterator itEnd   = ku.uend()-1;
        index_t nInterior = 0;
        for (gsKnotVector<real_t>::uiterator it = itBegin; it != itEnd; ++it) ++nInterior;
        GISMO_ENSURE(nInterior > 0, "duck-merge: patch " << p << " has no interior u-knot.");

        const index_t knotIdx = (n/nP) % nInterior;
        gsKnotVector<real_t>::uiterator uit = itBegin;
        for (index_t t = 0; t < knotIdx; ++t) ++uit;
        const real_t ustar = *uit;

        const gsMatrix<real_t> supp = patch.support();
        const real_t v = supp(1,0) + u01(rng)*(supp(1,1)-supp(1,0));

        gsMatrix<real_t> uv(2,1); uv << ustar, v;
        gsMatrix<real_t> xval, der;
        patch.eval_into(uv, xval);
        patch.deriv_into(uv, der);
        gsVector<real_t,3> Su, Sv;
        for (short_t kk = 0; kk != 3; ++kk) { Su[kk] = der(2*kk,0); Sv[kk] = der(2*kk+1,0); }
        const gsVector<real_t,3> nrm = duck.sgn[p] * Su.cross(Sv);

        short_t k = 0; real_t best = std::abs(nrm[0]);
        for (short_t kk = 1; kk != 3; ++kk)
            if (std::abs(nrm[kk]) > best) { best = std::abs(nrm[kk]); k = kk; }

        const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);
        const real_t yi = xval(i,0), yj = xval(j,0);

        std::vector<brc::Hit<real_t> > hits;
        const bool castOK = brc::castLine(els, k, yi, yj, hits, duck.L);
        if (!castOK) { ok = false; continue; }

        index_t nClose = 0;
        real_t bestDist = std::numeric_limits<real_t>::max();
        for (std::vector<brc::Hit<real_t> >::const_iterator h = hits.begin(); h != hits.end(); ++h)
        {
            gsMatrix<real_t> huv(2,1); huv << h->u, h->v;
            gsMatrix<real_t> hx;
            duck.mp.patch(h->patch).eval_into(huv, hx);
            const real_t dist = (hx.col(0)-xval.col(0)).norm();
            bestDist = std::min(bestDist, dist);
            if (dist <= (real_t)1e-12*duck.L) ++nClose;
        }
        if (nClose != 1) ok = false;
        worst = std::max(worst, bestDist/((real_t)1e-12*duck.L));

        if (!brc::windingOK(hits, k)) ok = false;
    }

    report(st, "duck-merge", ok ? worst : (real_t)1e300);
}

static void checkCubeAlignedR3(SelfTestState & st, const Geometry & ca, const std::vector<brc::BezElement<real_t> > & els)
{
    bool ok = true;
    real_t worst = 0;

    {
        std::vector<brc::Hit<real_t> > hits;
        const bool castOK = brc::castLine(els, (short_t)2, (real_t)0.4, (real_t)0.6, hits, ca.L);
        ok = ok && castOK && hits.size() == 2;
        if (ok)
        {
            ok = ok && (hits[0].t == (real_t)0.25) && (hits[1].t == (real_t)0.75);
            ok = ok && (hits[0].n[2] == (real_t)-1) && (hits[1].n[2] == (real_t)1);
        }

        int windingAtA; std::vector<brc::Hit<real_t> > in;
        brc::clipHits(hits, (real_t)0.5, (real_t)0.75, in, windingAtA);
        ok = ok && in.size() == 1 && in[0].t == (real_t)0.75 && windingAtA == 1;

        brc::clipHits(hits, (real_t)0.75, (real_t)1.0, in, windingAtA);
        ok = ok && in.empty() && windingAtA == 0;

        brc::clipHits(hits, (real_t)0.0, (real_t)0.25, in, windingAtA);
        ok = ok && in.size() == 1 && in[0].t == (real_t)0.25 && windingAtA == 0;
    }

    {
        std::vector<brc::Hit<real_t> > hits;
        const bool castOK = brc::castLine(els, (short_t)0, (real_t)0.6, (real_t)0.25, hits, ca.L);
        ok = ok && castOK && hits.size() == 2;
        if (ok && hits.size() == 2)
            ok = ok && (hits[0].t == (real_t)0.25) && (hits[1].t == (real_t)0.75);
    }

    {
        std::vector<brc::EdgeCurve<real_t> > edges = brc::elementEdges(els);
        ok = ok && edges.size() == 24;

        index_t total = 0;
        for (std::vector<brc::EdgeCurve<real_t> >::const_iterator e = edges.begin(); e != edges.end(); ++e)
        {
            std::vector<brc::PlaneCrossing<real_t> > out;
            brc::curvePlaneCrossings(*e, (short_t)2, (real_t)0.5, out);
            for (std::vector<brc::PlaneCrossing<real_t> >::const_iterator cc = out.begin(); cc != out.end(); ++cc)
                worst = std::max(worst, std::abs(cc->x[2]-(real_t)0.5)/(real_t)1e-15);
            total += (index_t)out.size();
        }
        ok = ok && total == 8;

        total = 0;
        for (std::vector<brc::EdgeCurve<real_t> >::const_iterator e = edges.begin(); e != edges.end(); ++e)
        {
            std::vector<brc::PlaneCrossing<real_t> > out;
            brc::curvePlaneCrossings(*e, (short_t)2, (real_t)0.25, out);
            const gsMatrix<real_t> & CC = e->hom.coefs();
            const real_t zFirst = CC(0,2)/CC(0,3);
            const real_t zLast  = CC(CC.rows()-1,2)/CC(CC.rows()-1,3);
            const int expectedDir = (zLast > zFirst) ? 1 : ((zLast < zFirst) ? -1 : 0);
            for (std::vector<brc::PlaneCrossing<real_t> >::const_iterator cc = out.begin(); cc != out.end(); ++cc)
            {
                worst = std::max(worst, std::abs(cc->x[2]-(real_t)0.25)/(real_t)1e-15);
                if (cc->dir != expectedDir) ok = false;
            }
            total += (index_t)out.size();
        }
        ok = ok && total == 8;
    }

    report(st, "cubeAligned-R3", ok ? worst : (real_t)1e300);
}

static void checkEdgeCrossings(SelfTestState & st, const std::string & gname, const Geometry & geo,
                               const std::vector<brc::BezElement<real_t> > & els, std::mt19937 & rng)
{
    real_t worst = 0;
    bool ok = true;
    std::vector<brc::EdgeCurve<real_t> > edges = brc::elementEdges(els);

    std::uniform_int_distribution<int> dk(0,2);
    for (int trial = 0; trial != 50; ++trial)
    {
        const short_t k = (short_t)dk(rng);
        std::uniform_real_distribution<real_t> dX(geo.bg(k,0), geo.bg(k,1));
        const real_t X = dX(rng);

        for (std::vector<brc::EdgeCurve<real_t> >::const_iterator e = edges.begin(); e != edges.end(); ++e)
        {
            std::vector<brc::PlaneCrossing<real_t> > cr;
            brc::curvePlaneCrossings(*e, k, X, cr);
            for (std::vector<brc::PlaneCrossing<real_t> >::const_iterator cc = cr.begin(); cc != cr.end(); ++cc)
                worst = std::max(worst, std::abs((*cc).x[k]-X)/((real_t)1e-13*geo.L));

            const index_t nS = 1000;
            index_t flips = 0;
            bool prev = false;
            for (index_t s = 0; s <= nS; ++s)
            {
                const real_t tt = e->hom.domainStart()
                                 + (e->hom.domainEnd()-e->hom.domainStart())*(real_t)s/(real_t)nS;
                gsMatrix<real_t> pt(1,1); pt(0,0) = tt;
                gsMatrix<real_t> val;
                e->hom.eval_into(pt, val);
                const bool above = (val(k,0)-X*val(3,0)) > 0;
                if (s > 0 && above != prev) ++flips;
                prev = above;
            }
            if ((index_t)cr.size() != flips) ok = false;
        }
    }
    report(st, gname+"-edge-crossings", ok ? worst : (real_t)1e300);
}

static int selftestAll(index_t nlines, const std::vector<real_t> & rot)
{
    SelfTestState st; st.nFail = 0; st.rng.seed(12345u);

    Geometry cube      = loadGeometry("cube", rot);
    Geometry cubeAlign = loadGeometry("cubeAligned", std::vector<real_t>());
    Geometry sphere    = loadGeometry("sphere", std::vector<real_t>());
    Geometry duck      = loadGeometry("duck", std::vector<real_t>());

    checkSphereStructure(st, sphere);
    checkSphereClosure(st, sphere, st.rng);

    std::vector<brc::BezElement<real_t> > elsCube      = brc::extractElements(cube.mp, cube.sgn);
    std::vector<brc::BezElement<real_t> > elsCubeAlign = brc::extractElements(cubeAlign.mp, cubeAlign.sgn);
    std::vector<brc::BezElement<real_t> > elsSphere    = brc::extractElements(sphere.mp, sphere.sgn);
    std::vector<brc::BezElement<real_t> > elsDuck      = brc::extractElements(duck.mp, duck.sgn);

    struct Entry
    {
        std::string name;
        const Geometry * g;
        const std::vector<brc::BezElement<real_t> > * els;
        bool hasAnalytic;
    };
    Entry entries[4] = {
        { "cube",        &cube,      &elsCube,      true  },
        { "cubeAligned", &cubeAlign, &elsCubeAlign, true  },
        { "sphere",      &sphere,    &elsSphere,    true  },
        { "duck",        &duck,      &elsDuck,      false },
    };

    for (int e = 0; e != 4; ++e)
    {
        checkNormalNumerator(st, entries[e].name+"-normal-numerator", *entries[e].g, *entries[e].els, st.rng);
        checkSignBounds(st, entries[e].name+"-sign-bounds", *entries[e].els, st.rng);
        checkHitsAndWinding(st, entries[e].name, *entries[e].g, *entries[e].els, nlines, st.rng);
        if (entries[e].hasAnalytic)
            checkAnalyticRoots(st, entries[e].name, *entries[e].g, *entries[e].els, nlines, st.rng);
        if (entries[e].name == "duck")
            checkDuckVolume(st, duck, elsDuck);
        else
            checkMoments(st, entries[e].name, *entries[e].g, *entries[e].els);
    }

    checkDuckOrientation(st, duck);
    checkDuckMerge(st, duck, elsDuck, st.rng);
    checkCubeAlignedR3(st, cubeAlign, elsCubeAlign);
    checkEdgeCrossings(st, "cube", cube, elsCube, st.rng);
    checkEdgeCrossings(st, "duck", duck, elsDuck, st.rng);

    if (st.nFail == 0) { gsInfo << "SELFTEST: ALL PASS\n"; return EXIT_SUCCESS; }
    gsInfo << "SELFTEST: " << st.nFail << " FAILED\n";
    return EXIT_FAILURE;
}

// =============================================================================
//  Octree + certification + the three-level leaf rule + the wall surface term
//  + the sign-sampling fallback: turns the ray-cast primitives above into a
//  background-grid quadrature over Geometry::bg.
// =============================================================================

/// The exact, half-open element/box overlap test on ::BezElement::lo/up
/// (never ::BezBox::overlaps, whose CLOSED test would double-count a face
/// lying exactly on a box plane): \a el overlaps [\a lo,\a hi] iff
/// el.up[d] > lo[d] && el.lo[d] <= hi[d] for every d.
static bool overlapR3(const brc::BezElement<real_t> & el, const real_t lo[3], const real_t hi[3])
{
    for (short_t d = 0; d != 3; ++d)
        if (!(el.up[d] > lo[d] && el.lo[d] <= hi[d])) return false;
    return true;
}

/// A Bezier element together with the number of ::brc::splitElement
/// generations separating it from its top-level element (0 at extraction).
struct SubEl
{
    brc::BezElement<real_t> el;
    int depth;
};

/// The 6 ordered (k,j) leaf-direction pairs, in the fixed tie-break order
/// used by ::certifyBox.
static const short_t g_pairOrder[6][2] = { {0,1}, {0,2}, {1,2}, {1,0}, {2,0}, {2,1} };

/// True iff \a el is an axis-aligned planar "wall" in direction \a k: one of
/// its two coordinates transverse to k is degenerate (lo == up, so the
/// element's own AABB is a face, never overlapping any e_k ray -- brc's
/// half-open AABB test at castLineOneElement, gsBRepRayCast.h:760, skips it
/// unconditionally) AND the element is bilinear (4 control points), so its
/// image equals that degenerate rectangle exactly (BSplineCube is trilinear,
/// hence every cubeAligned face element is bidegree (1,1)). A degenerate
/// element that is NOT bilinear is simply not a wall; k stays inadmissible
/// for it rather than triggering an ENSURE.
static bool isWallIn(const brc::BezElement<real_t> & el, short_t k)
{
    const short_t t0 = (short_t)((k+1)%3), t1 = (short_t)((k+2)%3);
    if (!(el.lo[t0] == el.up[t0] || el.lo[t1] == el.up[t1])) return false;
    return el.hom.basis().degree(0) == 1 && el.hom.basis().degree(1) == 1;
}

/// True iff the Bernstein coefficients of N_j on \a el are not of mixed
/// sign, up to the header's own certificate margin delta = 1e-10*Nscale
/// (identically zero is allowed).
static bool monoJ(const brc::BezElement<real_t> & el, short_t j)
{
    const real_t delta = (real_t)1e-10 * el.Nscale;
    return el.Nmin[j] >= -delta || el.Nmax[j] <= delta;
}

/// Pair (k,j) is ok for leaf \a s iff s is a wall in k (castLine skips it for
/// every e_k ray, so no certificate is needed there), or s carries a weak
/// certificate in k AND is monotone in j.
static bool pairOkForLeaf(const SubEl & s, short_t k, short_t j)
{
    if (isWallIn(s.el, k)) return true;
    return s.el.certWeak[k] != 0 && monoJ(s.el, j);
}

/// True iff \a s is "individually ok" for at least one of the 6 ordered pairs.
static bool leafOk(const SubEl & s)
{
    for (int t = 0; t != 6; ++t)
        if (pairOkForLeaf(s, g_pairOrder[t][0], g_pairOrder[t][1])) return true;
    return false;
}

/// Refines \a subs in place: repeatedly replaces any s with s.depth <
/// certSplit whose own AABB extent exceeds the box's, or that is not
/// individually ok for any pair, by its ::brc::splitElement children that
/// still overlap [\a lo,\a hi] (::overlapR3) -- until a fixed point.
static void refineSubs(std::vector<SubEl> & subs, const real_t lo[3], const real_t hi[3], index_t certSplit)
{
    const real_t extB = std::max(hi[0]-lo[0], std::max(hi[1]-lo[1], hi[2]-lo[2]));
    bool changed = true;
    while (changed)
    {
        changed = false;
        std::vector<SubEl> next;
        next.reserve(subs.size());
        for (std::vector<SubEl>::const_iterator s = subs.begin(); s != subs.end(); ++s)
        {
            const real_t extS = std::max(s->el.up[0]-s->el.lo[0],
                                 std::max(s->el.up[1]-s->el.lo[1], s->el.up[2]-s->el.lo[2]));
            const bool needSplit = (s->depth < (int)certSplit) && (extS > extB || !leafOk(*s));
            if (needSplit)
            {
                const std::vector<brc::BezElement<real_t> > kids = brc::splitElement(s->el);
                for (std::vector<brc::BezElement<real_t> >::const_iterator c = kids.begin(); c != kids.end(); ++c)
                    if (overlapR3(*c, lo, hi))
                    { SubEl ns; ns.el = *c; ns.depth = s->depth+1; next.push_back(ns); }
                changed = true;
            }
            else next.push_back(*s);
        }
        subs.swap(next);
    }
}

/// score(s,d) = |mean over rows of Ncoef(.,d)| / Nscale.
static real_t scoreOne(const SubEl & s, short_t d)
{
    return std::abs(s.el.Ncoef.col(d).mean()) / s.el.Nscale;
}

/// S_d for the box, direction \a k fixing which leaves count as "visible"
/// (not a wall in k): the min of ::scoreOne over those leaves, or 1 if none.
static real_t scoreDir(const std::vector<SubEl> & subs, short_t k, short_t d)
{
    real_t best = std::numeric_limits<real_t>::max();
    bool any = false;
    for (std::vector<SubEl>::const_iterator s = subs.begin(); s != subs.end(); ++s)
    {
        if (isWallIn(s->el, k)) continue;
        any = true;
        best = std::min(best, scoreOne(*s, d));
    }
    return any ? best : (real_t)1;
}

/// Picks the leaf direction pair (k,j): among the 6 ordered pairs admissible
/// for every leaf in \a subs (::pairOkForLeaf), maximises (S_k,S_j)
/// lexicographically (::scoreDir), ties broken by ::g_pairOrder. Returns
/// false if no pair is admissible (the box is uncertifiable).
static bool certifyBox(const std::vector<SubEl> & subs, short_t & kOut, short_t & jOut)
{
    bool found = false;
    real_t bestSk = 0, bestSj = 0;
    for (int t = 0; t != 6; ++t)
    {
        const short_t k = g_pairOrder[t][0], j = g_pairOrder[t][1];
        bool admissible = true;
        for (std::vector<SubEl>::const_iterator s = subs.begin(); s != subs.end(); ++s)
            if (!pairOkForLeaf(*s, k, j)) { admissible = false; break; }
        if (!admissible) continue;

        const real_t Sk = scoreDir(subs, k, k);
        const real_t Sj = scoreDir(subs, k, j);
        if (!found || Sk > bestSk || (Sk == bestSk && Sj > bestSj))
        { found = true; bestSk = Sk; bestSj = Sj; kOut = k; jOut = j; }
    }
    return found;
}

/// The moment set {x^a y^b z^c : a+b+c <= 2p}, built once per run.
struct MomentSet
{
    index_t p;
    std::vector<std::array<int,3> > triples;

    void build(index_t p_)
    {
        p = p_;
        triples.clear();
        const int maxDeg = (int)(2*p);
        for (int a = 0; a <= maxDeg; ++a)
        for (int b = 0; a+b <= maxDeg; ++b)
        for (int c = 0; a+b+c <= maxDeg; ++c)
        {
            std::array<int,3> tr; tr[0] = a; tr[1] = b; tr[2] = c;
            triples.push_back(tr);
        }
    }
};

/// Running quadrature state accumulated over one octree pass.
struct Accumulators
{
    std::vector<real_t> moments;

    real_t totalVol;
    index_t nVolNodes;
    index_t nNegWeight;

    real_t A, XN;
    gsVector<real_t,3> F;
    index_t nSurfNodes, nWallNodes;
    real_t lostSurface;
    real_t minAbsNk;

    real_t Afb;             // fallback-leaf surface area (sign-sampled)
    index_t nFbSurfNodes;   // fallback-leaf surface sample count
    real_t maxStopBracket;  // worst (xHi-xLo) at which a bisection was stopped by a castLine false

    index_t cellsIn, cellsOut, cellsCut;
    index_t certifiedLeaves[3];
    index_t uncutSubBoxes;
    index_t fallbackUncertifiable, fallbackCastLineFalse, fallbackFlipBudget;

    index_t castLineCalls;
    index_t castFalseNode, castFalseSample, castFalseBisect;

    std::vector<real_t>  cellVol;         // domain (inside) volume per cell
    std::vector<real_t>  cellBoxTileVol;  // sum of terminal-box geometric volumes per cell
    std::vector<index_t> cellVolNodes;    // volume-node count per cell
    std::vector<index_t> cellNTerminal;   // terminal-box count per cell
    std::vector<bool>    isCut;           // true iff this top-level cell classified cut (depth 0)

    Accumulators()
    : totalVol(0), nVolNodes(0), nNegWeight(0), A(0), XN(0),
      nSurfNodes(0), nWallNodes(0), lostSurface(0),
      minAbsNk(std::numeric_limits<real_t>::max()),
      Afb(0), nFbSurfNodes(0), maxStopBracket(0),
      cellsIn(0), cellsOut(0), cellsCut(0), uncutSubBoxes(0),
      fallbackUncertifiable(0), fallbackCastLineFalse(0), fallbackFlipBudget(0),
      castLineCalls(0), castFalseNode(0), castFalseSample(0), castFalseBisect(0)
    {
        F.setZero();
        certifiedLeaves[0] = certifiedLeaves[1] = certifiedLeaves[2] = 0;
    }

    void initCells(index_t N3)
    {
        cellVol.assign(N3, (real_t)0);
        cellBoxTileVol.assign(N3, (real_t)0);
        cellVolNodes.assign(N3, (index_t)0);
        cellNTerminal.assign(N3, (index_t)0);
        isCut.assign(N3, false);
    }

    index_t fallbackTotal() const
    { return fallbackUncertifiable + fallbackCastLineFalse + fallbackFlipBudget; }
    index_t castFalseTotal() const
    { return castFalseNode + castFalseSample + castFalseBisect; }
};

/// Adds the contribution of one quadrature node (\a x,\a y,\a z, weight \a w)
/// to every moment in \a ms, using powers of x,y,z built once per node.
static void addMomentNode(Accumulators & acc, const MomentSet & ms, real_t x, real_t y, real_t z, real_t w)
{
    const int maxDeg = (int)(2*ms.p);
    std::vector<real_t> xp(maxDeg+1), yp(maxDeg+1), zp(maxDeg+1);
    xp[0] = yp[0] = zp[0] = (real_t)1;
    for (int d = 1; d <= maxDeg; ++d)
    { xp[d] = xp[d-1]*x; yp[d] = yp[d-1]*y; zp[d] = zp[d-1]*z; }

    for (size_t t = 0; t != ms.triples.size(); ++t)
    {
        const std::array<int,3> & tr = ms.triples[t];
        acc.moments[t] += w * xp[tr[0]] * yp[tr[1]] * zp[tr[2]];
    }
}

/// Shared, read-only run configuration plus the mutable ::Accumulators,
/// threaded through the whole octree recursion.
struct Ctx
{
    const std::vector<brc::BezElement<real_t> > & els;
    const std::vector<BezBox> & bezBoxes;
    const std::vector<int> & allBoxIdx;
    const gsBRepSignedDist<real_t> & phi;
    gsVector<index_t> nq;    // nq[0]=level 1 (k), nq[1]=level 2 (j), nq[2]=level 3 (i, outer)
    index_t nsamp;
    index_t certSplit;
    index_t maxDepth;
    real_t L;
    const MomentSet & ms;
    Accumulators acc;
    index_t surfSub;      // ::fallbackSurface's uv sub-box refinement level, 2^surfSub per direction
    bool forceFallback;   // diagnostic: every cut leaf takes the fallback path (::processBox)
};

/// A single ::brc::castLine call at the point with p[j]=xJ, p[i]=yI (\a k,
/// \a j, \a i the chosen leaf directions). {(k+1)%3,(k+2)%3} is the same SET
/// as {i,j}, but not necessarily in that ORDER: which of i,j equals
/// (k+1)%3 depends on k. Building the full 3-point \a p and reading back
/// p[(k+1)%3], p[(k+2)%3] by actual coordinate index is what stays correct
/// regardless of that order, whereas passing (xJ,yI) or (yI,xJ) positionally
/// would silently swap the two transverse coordinates for half of the k
/// values. Counts the call in \a ctx regardless of outcome.
static bool castAt(Ctx & ctx, short_t k, short_t j, short_t i, real_t yI, real_t xJ,
                   std::vector<brc::Hit<real_t> > & hits)
{
    real_t p[3];
    p[j] = xJ;
    p[i] = yI;
    ++ctx.acc.castLineCalls;
    return brc::castLine(ctx.els, k, p[(k+1)%3], p[(k+2)%3], hits, ctx.L);
}

/// The level-2 signature at a sample x_j: the running winding count at z_B
/// (::brc::windingAt, NOT a parity) plus the ordered list of element indices
/// hit strictly inside (z_B,z_T].
struct Sig { int cnt; std::vector<index_t> elems; };

static bool sigEqual(const Sig & a, const Sig & b)
{ return a.cnt == b.cnt && a.elems == b.elems; }

static bool computeSig(Ctx & ctx, short_t k, short_t j, short_t i,
                       real_t yI, real_t xJ, real_t zB, real_t zT, Sig & sig)
{
    std::vector<brc::Hit<real_t> > hits;
    if (!castAt(ctx, k, j, i, yI, xJ, hits)) return false;

    sig.cnt = brc::windingAt(hits, k, zB);
    sig.elems.clear();
    for (std::vector<brc::Hit<real_t> >::const_iterator h = hits.begin(); h != hits.end(); ++h)
        if (zB < h->t && h->t <= zT) sig.elems.push_back(h->elem);
    return true;
}

/// Shared state of one level-2 line's bisection budget (4000 signature
/// evaluations, across every ::findFlips call triggered by that line).
struct BisectCtx
{
    index_t budget, used;
    bool failed, budgetExceeded;
    std::vector<real_t> emitted;
};

/// Locates the split points between a signature change at (\a xLo,\a sLo)
/// and (\a xHi,\a sHi) by predicate-flip bisection (never a line-plane
/// solve): recurses on genuine three-way splits, iterates on a two-way
/// split, and emits \a xHi once the interval can no longer be halved. A
/// ::brc::castLine false at the bisection midpoint (both \a xLo and \a xHi
/// already resolved) stops only that bisection: the midpoint is emitted as a
/// split point in place of the unreachable flip, \a ctx.acc.castFalseBisect
/// counts the event, and \a ctx.acc.maxStopBracket tracks the worst bracket
/// width this ever happened at -- a sliver of that width is then integrated
/// by the wrong side's smooth piece, costing at most O(bracket width) times
/// the integrand jump there. \a bctx.failed is not set by this branch: it
/// stays reserved for the flip budget (::BisectCtx), set only a few lines
/// above when \a bctx.used reaches \a bctx.budget.
static void findFlips(Ctx & ctx, short_t k, short_t j, short_t i, real_t yI, real_t zB, real_t zT,
                      real_t xLo, Sig sLo, real_t xHi, Sig sHi, BisectCtx & bctx)
{
    if (bctx.failed) return;
    for (;;)
    {
        const real_t mid = xLo + (real_t)0.5*(xHi-xLo);
        if (!(xLo < mid && mid < xHi)) { bctx.emitted.push_back(xHi); return; }

        if (bctx.used >= bctx.budget) { bctx.failed = true; bctx.budgetExceeded = true; return; }
        ++bctx.used;

        Sig sMid;
        if (!computeSig(ctx, k, j, i, yI, mid, zB, zT, sMid))
        {
            ++ctx.acc.castFalseBisect;
            ctx.acc.maxStopBracket = std::max(ctx.acc.maxStopBracket, xHi - xLo);
            bctx.emitted.push_back(mid);
            return;
        }

        if      (sigEqual(sMid, sLo)) { xLo = mid; }
        else if (sigEqual(sMid, sHi)) { xHi = mid; }
        else
        {
            findFlips(ctx, k, j, i, yI, zB, zT, xLo, sLo, mid, sMid, bctx);
            if (bctx.failed) return;
            findFlips(ctx, k, j, i, yI, zB, zT, mid, sMid, xHi, sHi, bctx);
            return;
        }
    }
}

/// Why ::evalLeaf did not commit a leaf.
enum class LeafFail { None, CastLineFalse, FlipBudget };

/// Staged output of one leaf, committed only on success (::LeafFail::None).
struct LeafBuffer
{
    std::vector<std::array<real_t,3> > volPts;
    std::vector<real_t>                volWts;
    std::vector<std::array<real_t,3> > surfPts;
    std::vector<std::array<real_t,3> > surfN;
    std::vector<real_t>                surfWts;
    real_t minAbsNk;
    LeafBuffer() : minAbsNk(std::numeric_limits<real_t>::max()) {}
};

/// The three-level certified leaf rule: a plain Gauss outer level along e_i
/// with no split points, a split-point Gauss level along e_j (split points
/// located by predicate-flip bisection), and an exact ray-cast level along
/// e_k providing both the inside sub-intervals (::brc::insideIntervals) for
/// the volume and every hit for the surface. The outer integrand, as a
/// function of y_i, is only piecewise smooth: it has kinks (and jumps, at
/// faces normal to e_i) wherever the level-2 split structure changes, and
/// plain Gauss across those is not exact. Between split points level 2 is
/// exact only for integrands polynomial in x_j (flat faces), so on curved
/// surfaces neither of the two outer levels reaches machine precision.
/// A ::brc::castLine false at a level-2 sample or a level-1 node fails the
/// leaf immediately (never read as "no hit"); a false at a ::findFlips
/// bisection midpoint does not (see ::findFlips) -- every site counts its
/// events by kind in \a ctx.acc.
static LeafFail evalLeaf(Ctx & ctx, const real_t lo[3], const real_t hi[3],
                         short_t k, short_t j, short_t i, LeafBuffer & buf)
{
    const real_t zB = lo[k], zT = hi[k];

    gsGaussRule<real_t> ruleI(ctx.nq[2]);
    gsMatrix<real_t> nodesI; gsVector<real_t> wtsI;
    ruleI.mapTo(lo[i], hi[i], nodesI, wtsI);

    for (index_t oi = 0; oi != nodesI.cols(); ++oi)
    {
        const real_t yI = nodesI(0, oi);
        const real_t Wi = wtsI[oi];

        // ---- Level 2: samples, signature and bisection -> split points ---
        const index_t M = ctx.nsamp;
        std::vector<real_t> xs(M+1);
        std::vector<Sig> sigs(M+1);
        xs[0] = lo[j]; xs[M] = hi[j];
        for (index_t m = 1; m < M; ++m)
            xs[m] = lo[j] + (hi[j]-lo[j]) * (real_t)m / (real_t)M;

        for (index_t m = 0; m <= M; ++m)
            if (!computeSig(ctx, k, j, i, yI, xs[m], zB, zT, sigs[m]))
            { ++ctx.acc.castFalseSample; return LeafFail::CastLineFalse; }

        BisectCtx bctx; bctx.budget = 4000; bctx.used = 0; bctx.failed = false; bctx.budgetExceeded = false;
        std::vector<real_t> splitPts;
        splitPts.push_back(lo[j]);
        for (index_t m = 0; m < M; ++m)
        {
            if (bctx.failed) break;
            if (sigEqual(sigs[m], sigs[m+1])) continue;

            bctx.emitted.clear();
            findFlips(ctx, k, j, i, yI, zB, zT, xs[m], sigs[m], xs[m+1], sigs[m+1], bctx);
            for (size_t e = 0; e != bctx.emitted.size(); ++e)
                splitPts.push_back(bctx.emitted[e]);
        }
        if (bctx.failed)
            return bctx.budgetExceeded ? LeafFail::FlipBudget : LeafFail::CastLineFalse;

        splitPts.push_back(hi[j]);
        std::sort(splitPts.begin(), splitPts.end());

        // ---- Level 2 Gauss on each remaining sub-interval, then level 1 ---
        for (size_t q = 0; q+1 < splitPts.size(); ++q)
        {
            const real_t sLo = splitPts[q], sHi = splitPts[q+1];
            if (!(sHi > sLo)) continue;

            gsGaussRule<real_t> ruleJ(ctx.nq[1]);
            gsMatrix<real_t> nodesJ; gsVector<real_t> wtsJ;
            ruleJ.mapTo(sLo, sHi, nodesJ, wtsJ);

            for (index_t oj = 0; oj != nodesJ.cols(); ++oj)
            {
                const real_t xJ = nodesJ(0, oj);
                const real_t Wj = wtsJ[oj];

                std::vector<brc::Hit<real_t> > hits;
                if (!castAt(ctx, k, j, i, yI, xJ, hits))
                { ++ctx.acc.castFalseNode; return LeafFail::CastLineFalse; }

                std::vector<std::pair<real_t,real_t> > iv;
                brc::insideIntervals(hits, k, zB, zT, iv);
                for (size_t r = 0; r != iv.size(); ++r)
                {
                    gsGaussRule<real_t> ruleK(ctx.nq[0]);
                    gsMatrix<real_t> nodesK; gsVector<real_t> wtsK;
                    ruleK.mapTo(iv[r].first, iv[r].second, nodesK, wtsK);
                    for (index_t ok = 0; ok != nodesK.cols(); ++ok)
                    {
                        std::array<real_t,3> pt;
                        pt[k] = nodesK(0,ok); pt[j] = xJ; pt[i] = yI;
                        buf.volPts.push_back(pt);
                        buf.volWts.push_back(Wi*Wj*wtsK[ok]);
                    }
                }

                for (std::vector<brc::Hit<real_t> >::const_iterator h = hits.begin(); h != hits.end(); ++h)
                {
                    if (!(zB < h->t && h->t <= zT)) continue;
                    std::array<real_t,3> pt, nn;
                    pt[k] = h->t; pt[j] = xJ; pt[i] = yI;
                    for (short_t d = 0; d != 3; ++d) nn[d] = h->n[d];
                    buf.surfPts.push_back(pt);
                    buf.surfN.push_back(nn);
                    buf.surfWts.push_back(Wi*Wj/std::abs(h->n[k]));
                    buf.minAbsNk = std::min(buf.minAbsNk, std::abs(h->n[k]));
                }
            }
        }
    }
    return LeafFail::None;
}

/// Adds a terminal box's own geometric volume to its cell's tiling sum
/// (::Accumulators::cellBoxTileVol), independent of in/out/leaf status --
/// this is the quantity the "leaf-tiling" check partitions [lo,hi] with.
static void markTerminal(Ctx & ctx, const real_t lo[3], const real_t hi[3], index_t cellId)
{
    const real_t boxVol = (hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2]);
    ctx.acc.cellBoxTileVol[cellId] += boxVol;
    ++ctx.acc.cellNTerminal[cellId];
}

/// An uncut box (::boxClassBRep cls != 0): outside boxes (\a cls > 0)
/// contribute nothing to the domain integral, so no quadrature is evaluated
/// at all there; inside boxes (\a cls < 0) get the full tensor Gauss rule
/// nq[0] x nq[1] x nq[2], exact for any polynomial the certified leaves
/// elsewhere in the same cell would also integrate exactly.
static void addUncutBox(Ctx & ctx, const real_t lo[3], const real_t hi[3], index_t cellId, int cls)
{
    markTerminal(ctx, lo, hi, cellId);
    if (cls > 0) return;

    gsVector<real_t> loV(3), hiV(3);
    loV << lo[0], lo[1], lo[2];
    hiV << hi[0], hi[1], hi[2];
    gsGaussRule<real_t> rule(ctx.nq);
    gsMatrix<real_t> pts; gsVector<real_t> wts;
    rule.mapTo(loV, hiV, pts, wts);

    real_t vol = 0;
    for (index_t c = 0; c != pts.cols(); ++c)
    {
        const real_t w = wts[c];
        addMomentNode(ctx.acc, ctx.ms, pts(0,c), pts(1,c), pts(2,c), w);
        vol += w;
        if (w < 0) ++ctx.acc.nNegWeight;
    }
    ctx.acc.totalVol += vol;
    ctx.acc.cellVol[cellId] += vol;
    ctx.acc.cellVolNodes[cellId] += (index_t)pts.cols();
    ctx.acc.nVolNodes += (index_t)pts.cols();
}

/// The exact rectangle surface term of one wall leaf \a el of the chosen
/// direction \a k: intersects the wall's own AABB rectangle with [lo,hi] in
/// the two coordinates transverse to its degenerate coordinate \a d, then
/// integrates the constant normal there with a tensor Gauss rule
/// nq[1] x nq[2] on the two in-plane coordinates e1 < e2 (ascending).
static void addWallContribution(Ctx & ctx, const brc::BezElement<real_t> & el, short_t k,
                                const real_t lo[3], const real_t hi[3])
{
    const short_t t0 = (short_t)((k+1)%3), t1 = (short_t)((k+2)%3);
    short_t d = (el.lo[t0] == el.up[t0]) ? t0 : t1;
    const real_t c = el.lo[d];

    short_t e[2]; int ei = 0;
    for (short_t dd = 0; dd != 3; ++dd) if (dd != d) e[ei++] = dd;

    const real_t alpha0 = std::max(el.lo[e[0]], lo[e[0]]);
    const real_t beta0  = std::min(el.up[e[0]], hi[e[0]]);
    const real_t alpha1 = std::max(el.lo[e[1]], lo[e[1]]);
    const real_t beta1  = std::min(el.up[e[1]], hi[e[1]]);
    if (beta0 <= alpha0 || beta1 <= alpha1) return;

    gsGaussRule<real_t> r0(ctx.nq[1]);
    gsMatrix<real_t> nodes0; gsVector<real_t> wts0;
    r0.mapTo(alpha0, beta0, nodes0, wts0);
    gsGaussRule<real_t> r1(ctx.nq[2]);
    gsMatrix<real_t> nodes1; gsVector<real_t> wts1;
    r1.mapTo(alpha1, beta1, nodes1, wts1);

    gsMatrix<real_t> uvC(2,1);
    uvC(0,0) = (real_t)0.5*(el.uv[0][0]+el.uv[0][1]);
    uvC(1,0) = (real_t)0.5*(el.uv[1][0]+el.uv[1][1]);
    gsMatrix<real_t> N;
    brc::normalNumerator(el.hom, el.sgn, uvC, N);
    const gsVector<real_t,3> nrm = N.col(0)/N.col(0).norm();

    for (index_t a = 0; a != nodes0.cols(); ++a)
    for (index_t b = 0; b != nodes1.cols(); ++b)
    {
        real_t pt[3];
        pt[d] = c; pt[e[0]] = nodes0(0,a); pt[e[1]] = nodes1(0,b);
        const real_t w = wts0[a]*wts1[b];
        ctx.acc.A += w;
        for (short_t dd = 0; dd != 3; ++dd) ctx.acc.F[dd] += w*nrm[dd];
        ctx.acc.XN += w*(pt[0]*nrm[0] + pt[1]*nrm[1] + pt[2]*nrm[2]);
        ++ctx.acc.nWallNodes;
    }
}

/// Commits a successful ::evalLeaf buffer: volume/surface accumulators, the
/// leaf's own tiling contribution, and the wall term of every wall leaf of
/// the chosen \a k among \a subs.
static void commitLeaf(Ctx & ctx, const LeafBuffer & buf, const std::vector<SubEl> & subs,
                       short_t k, const real_t lo[3], const real_t hi[3], index_t cellId)
{
    markTerminal(ctx, lo, hi, cellId);

    real_t leafVol = 0;
    for (size_t q = 0; q != buf.volWts.size(); ++q)
    {
        const real_t w = buf.volWts[q];
        const std::array<real_t,3> & pt = buf.volPts[q];
        addMomentNode(ctx.acc, ctx.ms, pt[0], pt[1], pt[2], w);
        leafVol += w;
        if (w < 0) ++ctx.acc.nNegWeight;
    }
    ctx.acc.totalVol += leafVol;
    ctx.acc.cellVol[cellId] += leafVol;
    ctx.acc.cellVolNodes[cellId] += (index_t)buf.volWts.size();
    ctx.acc.nVolNodes += (index_t)buf.volWts.size();

    for (size_t q = 0; q != buf.surfWts.size(); ++q)
    {
        const real_t w = buf.surfWts[q];
        const std::array<real_t,3> & nn = buf.surfN[q];
        const std::array<real_t,3> & pt = buf.surfPts[q];
        ctx.acc.A += w;
        for (short_t d = 0; d != 3; ++d) ctx.acc.F[d] += w*nn[d];
        ctx.acc.XN += w*(pt[0]*nn[0] + pt[1]*nn[1] + pt[2]*nn[2]);
        ++ctx.acc.nSurfNodes;
    }
    ctx.acc.minAbsNk = std::min(ctx.acc.minAbsNk, buf.minAbsNk);

    for (std::vector<SubEl>::const_iterator s = subs.begin(); s != subs.end(); ++s)
        if (isWallIn(s->el, k))
            addWallContribution(ctx, s->el, k, lo, hi);
}

/// A running Kahan-compensated sum (Kahan 1965): tracks the low-order bits
/// rounding drops on each \a add so that a long run of same-magnitude,
/// same-sign terms accumulates with error independent of the term count,
/// instead of growing linearly with it as a plain running sum does.
struct KahanSum
{
    real_t sum, c;
    KahanSum() : sum(0), c(0) {}
    void add(real_t x)
    {
        const real_t y = x - c;
        const real_t t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
};

/// Surface sign-sampling for a fallback leaf B = (\a lo,\a hi]: a composite
/// Gauss rule on the 2^surfSub x 2^surfSub uniform uv sub-boxes of every
/// \a subs element (::Ctx::surfSub), kept iff the physical image lies in B
/// under the exact half-open membership predicate.
///
/// Why a composite rule on uv sub-boxes, not ::brc::splitElement. Splitting
/// \a surfSub times and taking Gauss on each child is mathematically
/// identical (uniform uv midpoint splits, the same rule this integrates
/// piecewise), but ::brc::splitElement's ::detail::fillElementFields also
/// refits the Bernstein normal certificate (a linear solve) for every child,
/// which this rule does not need: only the geometry (::brc::normalNumerator
/// evaluated directly on the sub-box's own Gauss nodes) is used here.
///
/// Why the clamp into \a s.el's own AABB is required. A face with a constant
/// coordinate c (e.g. an axis-aligned wall) evaluates to c +- 1 ulp (the
/// header's exact-constant rule, gsBRepRayCast.h). ::overlapR3 assigned the
/// element to its owning box using that UNEVALUATED AABB; without the clamp,
/// the evaluated point can land on the wrong side of the plane and be kept
/// by no terminal box at all. Clamping moves the point by at most rounding
/// (it lies in the control-point AABB by the convex-hull property); the
/// clamp itself uses \a s.el's own AABB -- the one ::overlapR3 filtered
/// \a subs with -- never [\a lo,\a hi], or a point at a shared box boundary
/// could be clamped away from the side that actually owns it. The
/// membership test that follows is then the plain half-open predicate
/// against [\a lo,\a hi], the current leaf.
///
/// Partition argument. Terminal boxes tile the background box \a bg with
/// the half-open membership predicate, so every surface point belongs to
/// exactly one terminal box: a certified leaf's ray rule integrates exactly
/// Gamma cap B (hits kept iff zB < t <= zT, the wall term sharing
/// ::overlapR3's ownership), and a fallback leaf's \a subs cover Gamma cap B
/// because every piece with an image point in B has an AABB overlapping B.
/// So this rule neither double-counts nor drops surface area as a set; its
/// only error is quadrature error, concentrated in pieces straddling the
/// leaf boundary, where the membership indicator is itself sampled rather
/// than integrated exactly.
///
/// A, XN and each component of F are accumulated here in Kahan-compensated
/// running sums, folded into \a ctx.acc once at the end of the leaf: plain
/// running sums of the O(10^2 - 10^6) same-sign per-node weights this rule
/// produces (a full uv sub-box refinement of every sub-element, versus the
/// O(10) nodes ::addWallContribution ever summed per wall) drift with the
/// node count rather than its square root, since the rounding here is not
/// sign-random -- empirically 1e-12 to 7e-11 absolute at surfSub=2 on the
/// exact cubeAligned case (n0*2^r = 4, 8, 16) without compensation, measured
/// against a long-double running sum of the same terms; with compensation
/// both agree to machine precision.
static void fallbackSurface(Ctx & ctx, const std::vector<SubEl> & subs,
                            const real_t lo[3], const real_t hi[3])
{
    const index_t ns = index_t(1) << ctx.surfSub;
    gsVector<index_t> nn(2); nn << ctx.nq[0], ctx.nq[0];
    gsGaussRule<real_t> rule(nn);

    KahanSum kA, kXN, kF[3];
    index_t nKept = 0;

    for (std::vector<SubEl>::const_iterator s = subs.begin(); s != subs.end(); ++s)
    {
        const brc::BezElement<real_t> & el = s->el;
        const real_t u0 = el.uv[0][0], du = (el.uv[0][1]-u0)/(real_t)ns;
        const real_t v0 = el.uv[1][0], dv = (el.uv[1][1]-v0)/(real_t)ns;

        for (index_t a = 0; a != ns; ++a)
        for (index_t b = 0; b != ns; ++b)
        {
            gsMatrix<real_t> uvLo(2,1), uvHi(2,1);
            uvLo << u0 + a*du,     v0 + b*dv;
            uvHi << u0 + (a+1)*du, v0 + (b+1)*dv;
            gsMatrix<real_t> uv; gsVector<real_t> gw;
            rule.mapTo(uvLo, uvHi, uv, gw);

            gsMatrix<real_t> val;
            el.hom.eval_into(uv, val);
            gsMatrix<real_t> N;
            brc::normalNumerator(el.hom, el.sgn, uv, N);

            for (index_t c = 0; c != uv.cols(); ++c)
            {
                const real_t w = val(3,c);
                real_t x[3];
                for (short_t d = 0; d != 3; ++d)
                    x[d] = std::min(std::max(val(d,c)/w, el.lo[d]), el.up[d]);

                bool keep = true;
                for (short_t d = 0; d != 3; ++d)
                    if (!(lo[d] < x[d] && x[d] <= hi[d])) { keep = false; break; }
                if (!keep) continue;

                const real_t nrm = N.col(c).norm();
                if (nrm == 0) continue;

                const real_t weight = gw[c]*nrm/(w*w*w*w);
                const gsVector<real_t,3> nvec = N.col(c)/nrm;

                kA.add(weight);
                for (short_t d = 0; d != 3; ++d) kF[d].add(weight*nvec[d]);
                kXN.add(weight*(x[0]*nvec[0] + x[1]*nvec[1] + x[2]*nvec[2]));
                ++nKept;
            }
        }
    }

    ctx.acc.A += kA.sum;
    for (short_t d = 0; d != 3; ++d) ctx.acc.F[d] += kF[d].sum;
    ctx.acc.XN += kXN.sum;
    ctx.acc.Afb += kA.sum;
    ctx.acc.nFbSurfNodes += nKept;
}

/// A fallback leaf (uncertifiable box, or a certified box whose ::evalLeaf
/// failed, at \a maxDepth): sign-samples a plain tensor Gauss rule nq[0]^3
/// on [lo,hi] for the volume, keeping the nodes with phi < 0 at their own
/// Gauss weight, and integrates its own part of the surface Gamma cap B by
/// ::fallbackSurface on \a subs, the leaf's refined, ::overlapR3-filtered
/// sub-elements. Never calls ::addWallContribution: ::fallbackSurface's
/// uv-sub-box sampling already covers walls (degenerate flat elements) the
/// same way it covers curved pieces.
static void fallbackLeaf(Ctx & ctx, const real_t lo[3], const real_t hi[3],
                         const std::vector<SubEl> & subs, index_t cellId,
                         bool uncertifiable, LeafFail reason)
{
    markTerminal(ctx, lo, hi, cellId);

    if      (uncertifiable)                     ++ctx.acc.fallbackUncertifiable;
    else if (reason == LeafFail::FlipBudget)     ++ctx.acc.fallbackFlipBudget;
    else                                         ++ctx.acc.fallbackCastLineFalse;

    gsVector<real_t> loV(3), hiV(3);
    loV << lo[0], lo[1], lo[2];
    hiV << hi[0], hi[1], hi[2];
    gsVector<index_t> nqF(3); nqF << ctx.nq[0], ctx.nq[0], ctx.nq[0];
    gsGaussRule<real_t> rule(nqF);
    gsMatrix<real_t> pts; gsVector<real_t> wts;
    rule.mapTo(loV, hiV, pts, wts);

    gsMatrix<real_t> phiVals;
    ctx.phi.eval_into(pts, phiVals);

    real_t vol = 0; index_t nKept = 0;
    for (index_t c = 0; c != pts.cols(); ++c)
    {
        if (phiVals(0,c) >= 0) continue;
        const real_t w = wts[c];
        addMomentNode(ctx.acc, ctx.ms, pts(0,c), pts(1,c), pts(2,c), w);
        vol += w; ++nKept;
        if (w < 0) ++ctx.acc.nNegWeight;
    }
    ctx.acc.totalVol += vol;
    ctx.acc.cellVol[cellId] += vol;
    ctx.acc.cellVolNodes[cellId] += nKept;
    ctx.acc.nVolNodes += nKept;

    fallbackSurface(ctx, subs, lo, hi);
}

/// The octree driver, recursive in \a depth over one background cell
/// \a cellId: pre-classify via ::boxClassBRep (necessary condition only),
/// filter [\a lo,\a hi]'s own overlapping sub-elements out of \a parentSubs
/// and refine them (::refineSubs), certify a leaf direction pair
/// (::certifyBox) and run ::evalLeaf on success, else split into 8 children
/// sharing one midpoint (\a depth < maxDepth) or fall back
/// (::fallbackLeaf) at \a maxDepth. \a ctx.forceFallback (diagnostic, off by
/// default) skips ::certifyBox and ::evalLeaf entirely, forcing every cut
/// leaf down the same recurse-to-maxDepth-then-fallback path an
/// uncertifiable box would take -- the only way to exercise
/// ::fallbackSurface sharply on a geometry whose faces lie exactly on the
/// background grid's own planes, where the gate would otherwise certify
/// every leaf and the fallback path would never run.
static void processBox(Ctx & ctx, const real_t lo[3], const real_t hi[3], index_t depth,
                       const std::vector<SubEl> & parentSubs, index_t cellId)
{
    gsVector<real_t> loV(3), hiV(3);
    loV << lo[0], lo[1], lo[2];
    hiV << hi[0], hi[1], hi[2];
    const int cls = boxClassBRep(ctx.bezBoxes, ctx.allBoxIdx, ctx.phi, loV, hiV);

    if (cls != 0)
    {
        addUncutBox(ctx, lo, hi, cellId, cls);
        if (depth == 0) { if (cls < 0) ++ctx.acc.cellsIn; else ++ctx.acc.cellsOut; }
        else ++ctx.acc.uncutSubBoxes;
        return;
    }
    if (depth == 0) { ++ctx.acc.cellsCut; ctx.acc.isCut[cellId] = true; }

    std::vector<SubEl> subs;
    subs.reserve(parentSubs.size());
    for (std::vector<SubEl>::const_iterator s = parentSubs.begin(); s != parentSubs.end(); ++s)
        if (overlapR3(s->el, lo, hi)) subs.push_back(*s);

    refineSubs(subs, lo, hi, ctx.certSplit);

    if (subs.empty())
    {
        gsMatrix<real_t> ctr(3,1);
        ctr(0,0) = (real_t)0.5*(lo[0]+hi[0]);
        ctr(1,0) = (real_t)0.5*(lo[1]+hi[1]);
        ctr(2,0) = (real_t)0.5*(lo[2]+hi[2]);
        gsMatrix<real_t> val;
        ctx.phi.eval_into(ctr, val);
        addUncutBox(ctx, lo, hi, cellId, (val(0,0) < 0) ? -1 : 1);
        ++ctx.acc.uncutSubBoxes;
        return;
    }

    short_t k = 0, j = 0;
    bool certified = false;
    LeafFail fail = LeafFail::CastLineFalse;

    if (!ctx.forceFallback)
    {
        certified = certifyBox(subs, k, j);
        if (certified)
        {
            LeafBuffer buf;
            const short_t i = (short_t)(3 - k - j);
            fail = evalLeaf(ctx, lo, hi, k, j, i, buf);
            if (fail == LeafFail::None)
            {
                commitLeaf(ctx, buf, subs, k, lo, hi, cellId);
                ++ctx.acc.certifiedLeaves[k];
                return;
            }
        }
    }

    if (depth < ctx.maxDepth)
    {
        const real_t mid[3] = { (real_t)0.5*(lo[0]+hi[0]),
                                 (real_t)0.5*(lo[1]+hi[1]),
                                 (real_t)0.5*(lo[2]+hi[2]) };
        for (int cx = 0; cx != 2; ++cx)
        for (int cy = 0; cy != 2; ++cy)
        for (int cz = 0; cz != 2; ++cz)
        {
            const real_t clo[3] = { cx ? mid[0] : lo[0], cy ? mid[1] : lo[1], cz ? mid[2] : lo[2] };
            const real_t chi[3] = { cx ? hi[0] : mid[0], cy ? hi[1] : mid[1], cz ? hi[2] : mid[2] };
            processBox(ctx, clo, chi, depth+1, subs, cellId);
        }
    }
    else
        fallbackLeaf(ctx, lo, hi, subs, cellId, !certified, fail);
}

/// The scalar results of one octree pass, used both for the verbose checks
/// and for a single --ladder line.
struct RunResult
{
    real_t volErr, momErr, fluxErr, xnErr, areaErr;
    index_t fallbackCount, castFalseTotal;
    real_t time;
    bool allPass;
    index_t fbUncert, fbCast, fbBudget, cfNode, cfSample, cfBisect;
    real_t maxBracketRel, Afb;
};

static void printCheckLine(index_t & nFail, const std::string & label, bool pass, bool info,
                           const std::string & rest)
{
    gsInfo << (info ? "INFO" : (pass ? "PASS" : "FAIL")) << " " << label << " " << rest << "\n";
    if (!info && !pass) ++nFail;
}

/// Runs one full octree pass over Geometry::bg for \a geomName, prints the
/// header block and check lines when \a verbose, and returns the scalar
/// results ::RunResult (used unconditionally, including by --ladder).
static RunResult runOnce(const std::string & geomName, const std::vector<real_t> & rot,
                         index_t n0, index_t r, index_t maxDepth, index_t p,
                         gsVector<index_t> nq, index_t nsamp, index_t certSplit,
                         index_t surfSub, bool forceFallback, bool verbose)
{
    Geometry geo = loadGeometry(geomName, rot);
    std::vector<brc::BezElement<real_t> > els = brc::extractElements(geo.mp, geo.sgn);
    std::vector<BezBox> bezBoxes = makeBezBoxes(geo.mp);
    std::vector<int> allBoxIdx(bezBoxes.size());
    for (size_t t = 0; t != allBoxIdx.size(); ++t) allBoxIdx[t] = (int)t;
    gsBRepSignedDist<real_t> phi(geo.mp, bezBoxes, geo.sgn, geo.bg);

    MomentSet ms; ms.build(p);

    Ctx ctx = { els, bezBoxes, allBoxIdx, phi, nq, nsamp, certSplit, maxDepth, geo.L, ms, Accumulators(),
                surfSub, forceFallback };
    ctx.acc.moments.assign(ms.triples.size(), (real_t)0);

    const index_t N  = n0 * (index_t(1) << r);
    const index_t N3 = N*N*N;
    ctx.acc.initCells(N3);

    std::vector<real_t> Pd[3];
    for (short_t d = 0; d != 3; ++d)
    {
        Pd[d].resize(N+1);
        Pd[d][0] = geo.bg(d,0);
        Pd[d][N] = geo.bg(d,1);
        for (index_t m = 1; m < N; ++m)
            Pd[d][m] = geo.bg(d,0) + (geo.bg(d,1)-geo.bg(d,0)) * (real_t)m / (real_t)N;
    }

    struct CellBox { real_t lo[3], hi[3]; };
    std::vector<CellBox> cellBoxes(N3);

    gsStopwatch sw; sw.restart();
    for (index_t a = 0; a != N; ++a)
    for (index_t b = 0; b != N; ++b)
    for (index_t c = 0; c != N; ++c)
    {
        const real_t lo[3] = { Pd[0][a],   Pd[1][b],   Pd[2][c]   };
        const real_t hi[3] = { Pd[0][a+1], Pd[1][b+1], Pd[2][c+1] };
        const index_t cellId = a + N*(b + N*c);
        for (short_t d = 0; d != 3; ++d) { cellBoxes[cellId].lo[d] = lo[d]; cellBoxes[cellId].hi[d] = hi[d]; }

        std::vector<SubEl> topSubs;
        topSubs.reserve(els.size());
        for (size_t e = 0; e != els.size(); ++e)
        { SubEl s; s.el = els[e]; s.depth = 0; topSubs.push_back(s); }

        processBox(ctx, lo, hi, 0, topSubs, cellId);
    }
    const real_t elapsed = (real_t)sw.stop();

    const bool aligned = (geomName == "cubeAligned") && (((n0 << r) % 4) == 0);
    const real_t eps = std::numeric_limits<real_t>::epsilon();

    // ---- check 1: cell-sum -------------------------------------------------
    real_t cellVolSum = 0;
    for (index_t cId = 0; cId != N3; ++cId) cellVolSum += ctx.acc.cellVol[cId];
    const real_t cellSumTol = (real_t)ctx.acc.nVolNodes * eps * ctx.acc.totalVol;
    const bool cellSumPass = std::abs(cellVolSum - ctx.acc.totalVol) <= cellSumTol;

    // ---- check 2: cell-bounds ----------------------------------------------
    bool cellBoundsPass = (ctx.acc.nNegWeight == 0);
    real_t minRatio =  std::numeric_limits<real_t>::max();
    real_t maxRatio = -std::numeric_limits<real_t>::max();
    for (index_t cId = 0; cId != N3; ++cId)
    {
        const real_t Bc = (cellBoxes[cId].hi[0]-cellBoxes[cId].lo[0])
                         * (cellBoxes[cId].hi[1]-cellBoxes[cId].lo[1])
                         * (cellBoxes[cId].hi[2]-cellBoxes[cId].lo[2]);
        const real_t Vc = ctx.acc.cellVol[cId];
        const index_t nc = ctx.acc.cellVolNodes[cId];
        const real_t tol = Bc*(real_t)1e-14 + (real_t)nc*eps*Bc;
        if (!(Vc >= 0 && Vc <= Bc + tol)) cellBoundsPass = false;
        const real_t ratio = (Bc > 0) ? Vc/Bc : (real_t)0;
        minRatio = std::min(minRatio, ratio);
        maxRatio = std::max(maxRatio, ratio);
    }

    // ---- check 3: leaf-tiling -----------------------------------------------
    bool leafTilingPass = true;
    real_t worstTiling = 0;
    for (index_t cId = 0; cId != N3; ++cId)
    {
        if (!ctx.acc.isCut[cId]) continue;
        const real_t Bc = (cellBoxes[cId].hi[0]-cellBoxes[cId].lo[0])
                         * (cellBoxes[cId].hi[1]-cellBoxes[cId].lo[1])
                         * (cellBoxes[cId].hi[2]-cellBoxes[cId].lo[2]);
        const real_t diff = std::abs(ctx.acc.cellBoxTileVol[cId] - Bc);
        const real_t tol = (real_t)ctx.acc.cellNTerminal[cId] * (real_t)4 * eps * Bc;
        worstTiling = std::max(worstTiling, diff);
        if (diff > tol) leafTilingPass = false;
    }

    // ---- check 4: fallback ---------------------------------------------------
    const index_t fallbackCount = ctx.acc.fallbackTotal();
    const bool fallbackPass = (fallbackCount == 0);

    // ---- check 5: moments ------------------------------------------------
    const real_t rho = geo.bg.cwiseAbs().maxCoeff();
    const real_t V0 = momentOracle(els, 0, 0, 0);
    real_t maxMomErr = 0, volErr = 0;
    for (size_t t = 0; t != ms.triples.size(); ++t)
    {
        const std::array<int,3> & tr = ms.triples[t];
        const real_t ref = momentOracle(els, tr[0], tr[1], tr[2]);
        const int n = tr[0]+tr[1]+tr[2];
        const real_t err = std::abs(ctx.acc.moments[t] - ref) / (V0 * std::pow(rho, (real_t)n));
        maxMomErr = std::max(maxMomErr, err);
        if (tr[0] == 0 && tr[1] == 0 && tr[2] == 0) volErr = err;
    }
    const bool momPass = (maxMomErr <= (real_t)1e-12);

    // ---- checks 6-8: surface ----------------------------------------------
    gsVector<real_t,3> fluxNRef; real_t volRef, areaRef;
    brepIntegrals(geo.mp, geo.sgn, fluxNRef, volRef, areaRef);
    const real_t fluxErr = ctx.acc.F.norm() / areaRef;
    const real_t xnErr   = std::abs(ctx.acc.XN - (real_t)3*volRef) / ((real_t)3*volRef);
    const real_t areaErr = std::abs(ctx.acc.A - areaRef) / areaRef;
    const bool fluxPass = (fluxErr <= (real_t)1e-12);
    const bool xnPass   = (xnErr   <= (real_t)1e-12);
    const bool areaPass = (areaErr <= (real_t)1e-12);

    index_t nFail = 0;
    if (verbose)
    {
        gsInfo << "geometry '" << geomName << "': " << geo.mp.nPatches() << " patches, "
               << els.size() << " Bezier elements\n";
        gsInfo << "orientation signs:";
        for (size_t pp = 0; pp != geo.sgn.size(); ++pp) gsInfo << (geo.sgn[pp] > 0 ? " +" : " -");
        gsInfo << "\n";
        gsInfo << "|oint n dS| = " << fmtSci(fluxNRef.norm()) << "\n";
        gsInfo << "area        = " << std::setprecision(10) << areaRef << "\n";
        gsInfo << "V (oracle)  = " << std::setprecision(10) << V0 << "\n";
        gsInfo << "background box (rows x,y,z; cols lo,hi):\n" << geo.bg << "\n";

        int weak[3] = {0,0,0}, strict[3] = {0,0,0}, uncert[3] = {0,0,0};
        for (std::vector<brc::BezElement<real_t> >::const_iterator el = els.begin(); el != els.end(); ++el)
            for (short_t k = 0; k != 3; ++k)
            {
                if      (el->certStrict[k] != 0) ++strict[k];
                else if (el->certWeak[k]   != 0) ++weak[k];
                else                              ++uncert[k];
            }
        gsInfo << "weak/strict/uncertified elements per k:\n";
        for (short_t k = 0; k != 3; ++k)
            gsInfo << "  k=" << k << ": weak=" << weak[k] << " strict=" << strict[k]
                   << " uncertified=" << uncert[k] << "\n";

        gsInfo << "grid: n0=" << n0 << " r=" << r << " N=" << N << " maxDepth=" << maxDepth
               << " p=" << p << " nq=(" << nq[0] << "," << nq[1] << "," << nq[2] << ")"
               << " nsamp=" << nsamp << " certSplit=" << certSplit
               << " aligned=" << (aligned ? 1 : 0);
        if (ctx.forceFallback) gsInfo << " forceFallback=1";
        gsInfo << "\n";
        gsInfo << "cells: total=" << N3 << " in=" << ctx.acc.cellsIn
               << " out=" << ctx.acc.cellsOut << " cut=" << ctx.acc.cellsCut << "\n";
        gsInfo << "leaves: certified=" << (ctx.acc.certifiedLeaves[0]+ctx.acc.certifiedLeaves[1]+ctx.acc.certifiedLeaves[2])
               << " (k=0:" << ctx.acc.certifiedLeaves[0] << " k=1:" << ctx.acc.certifiedLeaves[1]
               << " k=2:" << ctx.acc.certifiedLeaves[2] << ")"
               << " uncut=" << ctx.acc.uncutSubBoxes << " fallback=" << fallbackCount
               << " (uncertifiable=" << ctx.acc.fallbackUncertifiable
               << " castLineFalse=" << ctx.acc.fallbackCastLineFalse
               << " flipBudget=" << ctx.acc.fallbackFlipBudget << ")\n";
        gsInfo << "castLine: calls=" << ctx.acc.castLineCalls
               << " false(node/sample/bisect)=" << ctx.acc.castFalseNode << "/"
               << ctx.acc.castFalseSample << "/" << ctx.acc.castFalseBisect << "\n";
        gsInfo << "nodes: volume=" << ctx.acc.nVolNodes << " surface=" << ctx.acc.nSurfNodes
               << " wall=" << ctx.acc.nWallNodes << "  lostSurface=" << fmtSci(ctx.acc.lostSurface)
               << "  min|n_k|=" << fmtSci(ctx.acc.minAbsNk) << "  time=" << fmtSci(elapsed) << "s\n";
        gsInfo << "bisect-stop: count=" << ctx.acc.castFalseBisect
               << " maxBracket=" << fmtSci(ctx.acc.maxStopBracket)
               << " maxBracket/L=" << fmtSci(ctx.acc.maxStopBracket/geo.L) << "\n";
        gsInfo << "fallback-surface: nodes=" << ctx.acc.nFbSurfNodes
               << " A=" << fmtSci(ctx.acc.Afb) << " surfSub=" << surfSub << "\n";

        std::ostringstream r1; r1 << "err=" << fmtSci(std::abs(cellVolSum-ctx.acc.totalVol)) << " tol=" << fmtSci(cellSumTol);
        printCheckLine(nFail, "cell-sum", cellSumPass, false, r1.str());

        std::ostringstream r2; r2 << "minRatio=" << fmtSci(minRatio) << " maxRatio=" << fmtSci(maxRatio)
                                  << " negWeights=" << ctx.acc.nNegWeight;
        printCheckLine(nFail, "cell-bounds", cellBoundsPass, false, r2.str());

        std::ostringstream r3; r3 << "worst=" << fmtSci(worstTiling);
        printCheckLine(nFail, "leaf-tiling", leafTilingPass, false, r3.str());

        std::ostringstream r4; r4 << "count=" << fallbackCount;
        printCheckLine(nFail, "fallback", fallbackPass, !aligned || ctx.forceFallback, r4.str());

        std::ostringstream r5; r5 << "maxScaledErr=" << fmtSci(maxMomErr)
                                  << " (" << ms.triples.size() << " moments, a+b+c<=" << (2*p) << ")";
        printCheckLine(nFail, "moments", momPass, !aligned || ctx.forceFallback, r5.str());

        std::ostringstream r6; r6 << "err=" << fmtSci(fluxErr);
        printCheckLine(nFail, "surface-flux", fluxPass, !aligned, r6.str());

        std::ostringstream r7; r7 << "err=" << fmtSci(xnErr);
        printCheckLine(nFail, "surface-xn", xnPass, !aligned, r7.str());

        std::ostringstream r8; r8 << "err=" << fmtSci(areaErr) << " A=" << std::setprecision(10) << ctx.acc.A
                                  << " Aref=" << std::setprecision(10) << areaRef;
        printCheckLine(nFail, "surface-area", areaPass, !aligned, r8.str());

        if (nFail == 0) gsInfo << "CHECKS: ALL PASS\n";
        else            gsInfo << "CHECKS: " << nFail << " FAILED\n";
    }

    const bool volReqOk = ctx.forceFallback || (fallbackPass && momPass);

    RunResult res;
    res.volErr = volErr;
    res.momErr = maxMomErr;
    res.fluxErr = fluxErr;
    res.xnErr = xnErr;
    res.areaErr = areaErr;
    res.fallbackCount = fallbackCount;
    res.castFalseTotal = ctx.acc.castFalseTotal();
    res.time = elapsed;
    res.allPass = cellSumPass && cellBoundsPass && leafTilingPass
                && (!aligned || (volReqOk && fluxPass && xnPass && areaPass));
    res.fbUncert = ctx.acc.fallbackUncertifiable;
    res.fbCast   = ctx.acc.fallbackCastLineFalse;
    res.fbBudget = ctx.acc.fallbackFlipBudget;
    res.cfNode   = ctx.acc.castFalseNode;
    res.cfSample = ctx.acc.castFalseSample;
    res.cfBisect = ctx.acc.castFalseBisect;
    res.maxBracketRel = ctx.acc.maxStopBracket / geo.L;
    res.Afb = ctx.acc.Afb;
    return res;
}

/// --ladder: an nq-ladder (nq=1..6 at the given maxDepth) then a
/// maxDepth-ladder (0..maxDepth at the given \a nqGiven), one LADDER line
/// each; INFO only, the exit code is unaffected.
static void runLadders(const std::string & geomName, const std::vector<real_t> & rot,
                       index_t n0, index_t r, index_t maxDepth, index_t p,
                       index_t nsamp, index_t certSplit, index_t surfSub, bool forceFallback,
                       const gsVector<index_t> & nqGiven)
{
    for (index_t n = 1; n <= 6; ++n)
    {
        gsVector<index_t> nq(3); nq << n, n, n;
        const RunResult res = runOnce(geomName, rot, n0, r, maxDepth, p, nq, nsamp, certSplit,
                                      surfSub, forceFallback, false);
        gsInfo << "LADDER nq=" << n << " maxDepth=" << maxDepth
               << " vol=" << fmtSci(res.volErr) << " mom=" << fmtSci(res.momErr)
               << " flux=" << fmtSci(res.fluxErr) << " xn=" << fmtSci(res.xnErr)
               << " area=" << fmtSci(res.areaErr) << " fallback=" << res.fallbackCount
               << " castFalse=" << res.castFalseTotal << " time=" << fmtSci(res.time)
               << " fb(u/c/b)=" << res.fbUncert << "/" << res.fbCast << "/" << res.fbBudget
               << " cf(n/s/b)=" << res.cfNode << "/" << res.cfSample << "/" << res.cfBisect
               << " maxBr=" << fmtSci(res.maxBracketRel) << " Afb=" << fmtSci(res.Afb) << "\n";
    }
    for (index_t d = 0; d <= maxDepth; ++d)
    {
        const RunResult res = runOnce(geomName, rot, n0, r, d, p, nqGiven, nsamp, certSplit,
                                      surfSub, forceFallback, false);
        gsInfo << "LADDER nq=" << nqGiven[0] << " maxDepth=" << d
               << " vol=" << fmtSci(res.volErr) << " mom=" << fmtSci(res.momErr)
               << " flux=" << fmtSci(res.fluxErr) << " xn=" << fmtSci(res.xnErr)
               << " area=" << fmtSci(res.areaErr) << " fallback=" << res.fallbackCount
               << " castFalse=" << res.castFalseTotal << " time=" << fmtSci(res.time)
               << " fb(u/c/b)=" << res.fbUncert << "/" << res.fbCast << "/" << res.fbBudget
               << " cf(n/s/b)=" << res.cfNode << "/" << res.cfSample << "/" << res.cfBisect
               << " maxBr=" << fmtSci(res.maxBracketRel) << " Afb=" << fmtSci(res.Afb) << "\n";
    }
}

// =============================================================================
int main(int argc, char * argv[])
{
    std::string geomName = "cube";
    std::vector<real_t> rot;
    index_t n0 = 4;
    index_t r = 0;
    index_t maxDepth = 3;
    index_t p = 2;
    std::vector<index_t> nqIn;
    index_t nsamp = 8;
    index_t certSplit = 8;
    index_t nlines = 200;
    index_t surfSub = 2;
    bool selftest = false;
    bool ladder = false;
    bool forceFallback = false;

    gsCmdLine cmd("Background-grid quadrature on a closed spline BRep via ray casting "
                 "(gsBRepRayCast.h): a uniform grid over Geometry::bg, an octree refining "
                 "every cut cell, a three-level certified leaf rule (ray-cast hits + "
                 "bisected splits + plain Gauss), an exact wall term for axis-aligned "
                 "planar boundary elements, and a sign-sampling fallback that integrates "
                 "both the sign-sampled volume and, by uv-sub-box Gauss sampling, its own "
                 "part of the surface, with a --selftest battery over the ray-cast "
                 "primitives.");
    cmd.addString("", "geom", "Geometry: cube, cubeAligned, sphere, duck", geomName);
    cmd.addMultiReal("", "rot", "Cube rotation angles phiz phiy phix (repeat 3 times)", rot);
    cmd.addInt("", "n0", "Background grid: n0*2^r cells per direction over Geometry::bg", n0);
    cmd.addInt("r", "refine", "Uniform grid refinement steps", r);
    cmd.addInt("", "maxDepth", "Maximum octree depth inside a cut cell", maxDepth);
    cmd.addInt("p", "degree", "Moment set a+b+c <= 2p; default nq = p+1", p);
    cmd.addMultiInt("", "nq", "Gauss node counts (0, 1 or 3 values -> all levels, or levels 1,2,3; default p+1)", nqIn);
    cmd.addInt("", "nsamp", "Level-2 sample count M along e_j", nsamp);
    cmd.addInt("", "certSplit", "Maximum sub-element refinement depth for certification", certSplit);
    cmd.addInt("", "nlines", "Random lines per direction in --selftest", nlines);
    cmd.addInt("", "surfSub", "Fallback-leaf surface sampling: 2^surfSub uv sub-boxes per "
              "sub-element per direction (subs that exhausted certSplit before reaching "
              "leaf scale may still exceed the leaf extent)", surfSub);
    cmd.addSwitch("selftest", "Run the fixed self-test battery over all 4 geometries", selftest);
    cmd.addSwitch("ladder", "Print an nq- and a maxDepth-ladder after the checks", ladder);
    cmd.addSwitch("forceFallback", "Diagnostic: force every cut leaf down the fallback path "
                 "(skips certifyBox/evalLeaf), to exercise the surface sign-sampling rule "
                 "sharply", forceFallback);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (selftest)
        return selftestAll(nlines, rot);

    gsVector<index_t> nq(3);
    if      (nqIn.empty())     { nq[0] = nq[1] = nq[2] = p+1; }
    else if (nqIn.size() == 1) { nq[0] = nq[1] = nq[2] = nqIn[0]; }
    else if (nqIn.size() == 3) { nq[0] = nqIn[0]; nq[1] = nqIn[1]; nq[2] = nqIn[2]; }
    else GISMO_ENSURE(false, "--nq needs 0, 1 or 3 values, got " << nqIn.size() << ".");

    const RunResult res = runOnce(geomName, rot, n0, r, maxDepth, p, nq, nsamp, certSplit,
                                  surfSub, forceFallback, true);

    if (ladder)
        runLadders(geomName, rot, n0, r, maxDepth, p, nsamp, certSplit, surfSub, forceFallback, nq);

    return res.allPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
