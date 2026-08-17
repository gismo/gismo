/** @file immersed_brep_quadrature_example.cpp

    @brief Immersed (cut-cell) quadrature points generated directly from a
           spline BREP -- a watertight gsMultiPatch of untrimmed spline
           surfaces embedded in R^3.

    This is the spline-BREP counterpart of
    immersed_fcm_octree_adaptive_volume_example.cpp (triangle meshes, SAT
    classification) and of example_cutcell_levelset_occ.cpp (OpenCascade, 2D).
    The geometric input here is neither a mesh nor an implicit function: it is
    a boundary representation made of tensor-product spline patches, e.g.
    filedata/breps/3D/duck_BRep.xml (6 biquadratic patches, 12 interfaces, no
    free boundary -> a closed 2-manifold shell).

    Pipeline
    --------
      1. Load the BREP as a gsMultiPatch<real_t> (parDim 2, geoDim 3) and
         verify it is watertight (nBoundary() == 0).

      2. Build a table of *Bezier-element boxes*: for every knot-span element
         of every patch, the axis-aligned bounding box of that element's
         active control points.  By the convex-hull property this bounds the
         surface over the element.  It is the spline analogue of the
         per-triangle AABBs used by the mesh version, and there are far fewer
         of them (~630 for the duck), so a flat array with a linear-scan
         overlap filter is enough -- no AABB tree.

      3. Fix the patch orientations.  Parametric interface matching does NOT
         imply consistent normal orientation, and the duck is a concrete
         counter-example: as loaded, |oint n dS| = 3.09 (must be 0 for a
         closed consistently-oriented shell) and the divergence-theorem volume
         is not translation invariant.  We recover a per-patch sign s_p by BFS
         over the shell, comparing the *induced boundary traversal direction*
         along each shared edge (opposite traversal <=> consistent
         orientation).  The tangent test is used rather than a normal test
         because tangents along a shared edge are parallel by construction,
         whereas normals across a C0 crease need not be.  The result is
         validated globally with |oint n dS| ~ 0 and V > 0.  The geometry is
         never modified; s_p is only a multiplier on x_u x x_v.

         The shared edges are found GEOMETRICALLY, not from the <MultiPatch>
         block: in duck_BRep.xml the stored patch adjacency is right but the
         stored side indices are not -- none of the 12 declared interfaces
         joins two coincident edges (Hausdorff distance ~1 rather than ~0).
         Matching sides by geometry is both correct for that file and
         independent of whether a topology block is present at all.  For the
         duck the 12 true pairs match to ~1e-14 while the next-best candidate
         is 0.63 away, so the pairing is unambiguous.

      4. gsBRepSignedDist<T>: a gsFunction giving the signed distance to the
         BREP, phi < 0 inside (the convention of gsCutCellRule and of the
         other immersed examples).  Unsigned distance comes from a
         branch-and-bound over the Bezier boxes plus a box-projected
         Gauss-Newton footpoint solve; the sign comes from the oriented normal
         at the footpoint.  The gradient is analytic.

      5. Volume quadrature on the cut cells comes from gsAlgoimAdaptiveRule,
         exactly as in the mesh-based example.  Boundary quadrature does NOT
         go through gsAlgoimGenericRule with dim = D: for a BREP the surface is
         the input geometry, so Gauss points are placed on the patches directly
         (see brepSurfaceQuadrature below for why, and for what the level-set
         route does instead).

    Verification anchors (duck_BRep.xml, original coordinates)
    ---------------------------------------------------------
      volume  V    = 1.1905049      ( = 1/3 oint x.n dS,   polynomial => exact )
      surface area = 7.3954796     ( = oint |x_u x x_v|,   needs subdivision )
    Both are printed by step 3 and are computed independently of the
    quadrature, so they validate the whole pipeline.

      6. The volume rule of step 5 can optionally be swapped for the
         moment-fitting rule (--momentFit): the adaptive point cloud of every
         cut cell is then compressed onto the fixed tensor Gauss grid of that
         cell, so the number of quadrature points per cut cell no longer
         depends on the octree depth.  The surface quadrature of step 5 is
         never compressed and is untouched by the flag.
         Caveat: the moment-fitting rule exposes no BoxClassifier hook, so
         inside it the octree children are classified by the midpoint +
         Lipschitz test on phi instead of the Bezier-box test used on the plain
         adaptive path (the background cells handed to the rule are still
         selected by boxClassBRep).

    Usage:
      ./bin/immersed_brep_quadrature_example -r 4 --quadDepth 2 -o out_duck --plot
      ./bin/immersed_brep_quadrature_example -r 4 --quadDepth 2 --momentFit

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
//  Bezier-element box: one knot-span element of one patch, together with the
//  AABB of its active control points (a convex-hull bound on the surface).
//  Replaces the per-triangle AABB array of the mesh-based example.
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

// -----------------------------------------------------------------------------
//  Small geometry helpers on a 2->3 spline patch
// -----------------------------------------------------------------------------

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
//  Orientation fix-up
// =============================================================================

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
/// \a target, found by dense sampling of the edge.  The two sides of a matched
/// pair are the same curve, so the match is exact up to the sampling
/// resolution -- which is then only used to evaluate a tangent, a smooth
/// quantity along the edge.
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
/// with the per-patch orientation multipliers \a sgn.
static void brepIntegrals(const gsMultiPatch<real_t> & brep,
                          const std::vector<real_t> & sgn,
                          gsVector<real_t,3> & fluxN, real_t & volume, real_t & area)
{
    fluxN.setZero(); volume = 0; area = 0;

    // x.(x_u x x_v) and (x_u x x_v) are polynomial (degree <= 5 per direction
    // for these biquadratic patches), so Gauss alone is exact for the volume
    // and for the closure test.  |x_u x x_v| is a square root and is NOT, so the
    // area needs element subdivision as well: with 6 nodes on whole elements it
    // is only accurate to ~2e-4, which would otherwise be misread as an error in
    // the surface quadrature.
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
/// consistently outward.  Returns the multipliers; throws (GISMO_ENSURE) if the
/// result does not close.
static std::vector<real_t> fixOrientation(const gsMultiPatch<real_t> & brep,
                                          real_t & refVolume, real_t & refArea,
                                          bool verbose = true)
{
    const size_t nP = brep.nPatches();
    std::vector<real_t> sgn(nP, 0.0);        // 0 = not yet visited

    // ---- pair up the patch sides GEOMETRICALLY ----------------------------
    // The <MultiPatch> topology is deliberately not used here.  In
    // filedata/breps/3D/duck_BRep.xml the recorded patch *adjacency* is right
    // but the recorded *side indices* are not: none of the 12 declared
    // interfaces joins two coincident edges (Hausdorff distance ~1 instead of
    // ~0).  Matching the sides by geometry is both correct for that file and
    // independent of whether a topology block is present at all.
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

            // Vote over a few interior samples of the shared edge.  Consistent
            // orientation <=> the two faces traverse the edge in *opposite*
            // directions, i.e. tS . tO < 0.  Tangents along a shared edge are
            // parallel by construction, so this test stays well conditioned
            // even across a C0 crease -- unlike comparing the two normals.
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

            // vote < 0 -> opposite traversal -> same orientation as p.
            sgn[q] = (vote < 0) ? sgn[p] : -sgn[p];
            queue.push_back(q);
        }
    }

    for (size_t p = 0; p != nP; ++p)
        GISMO_ENSURE(sgn[p] != 0.0,
                     "Patch " << p << " is not connected to patch 0: the BREP is "
                     "not a single connected shell.");

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
        gsInfo << "  sides paired geometrically    : " << nE/2 << "  (watertight)\n";
        gsInfo << "  patch orientation multipliers :";
        for (size_t p = 0; p != nP; ++p) gsInfo << (sgn[p] > 0 ? " +" : " -");
        gsInfo << "\n";
        gsInfo << "  |oint n dS| (closure)         : " << fluxN.norm() << "\n";
        gsInfo << "  surface area                  : "
               << std::setprecision(8) << refArea << "\n";
        gsInfo << "  volume  1/3 oint x.n dS       : "
               << std::setprecision(8) << refVolume << "\n";
    }

    // A closed, consistently oriented shell has oint n dS = 0 exactly.  This
    // assert is what catches a mis-voted flip (which happens where the tangent
    // test is ill-conditioned), so it must fail loudly rather than be printed.
    GISMO_ENSURE(fluxN.norm() < 1e-8 * std::max(refArea, (real_t)1),
                 "BREP orientation fix failed: |oint n dS| = " << fluxN.norm()
                 << " (surface area " << refArea << "). The input is either not "
                 "watertight or the per-patch orientation could not be resolved.");
    GISMO_ENSURE(refVolume > 0, "Non-positive enclosed volume: " << refVolume);

    return sgn;
}

// =============================================================================
//  gsBRepSignedDist<T> : signed distance to a spline BREP
//
//  Convention (matches gsCutCellRule and the other immersed examples):
//    phi < 0  -> inside      phi > 0  -> outside
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
#       pragma omp parallel for if(u.cols() > 32)
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
#       pragma omp parallel for if(u.cols() > 32)
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
//  Direct surface quadrature on the BREP.
//
//  For an immersed BREP the boundary does not have to be recovered from the
//  level set at all: it *is* the input geometry, exactly and in parametrised
//  form.  Placing Gauss points on the Bezier elements of the patches therefore
//  integrates over the whole boundary to quadrature accuracy (the total weight
//  reproduces the surface area to ~1e-11 here) with exact unit normals.
//
//  Each Bezier element is first split into nSub x nSub parameter sub-boxes so
//  that no sub-box image is much larger than a background cell; every resulting
//  point then belongs to a well-defined cut cell, which is what an assembler
//  needs for Nitsche terms.
//
//  LIMITATION: only the *total* is that accurate.  The per-cell split is
//  O(sub-element size): a Gauss point near a cell face is assigned to whichever
//  cell contains it, while its weight represents a patch of surface that may
//  stick slightly outside that cell.  nSub bounds this but does not remove it.
//  Users needing an exact per-cell split must clip the sub-boxes against the
//  cell faces in parameter space.
//
//  This replaces gsAlgoimGenericRule with dim = D for the surface part.  That
//  rule is exact for the smooth polynomial level sets it is tested with, but it
//  is driven by the interval-arithmetic bound of gsAlgoimFunctionWrapper, whose
//  remainder eps is estimated from finite-difference second derivatives
//  (gsAlgoimFunctionWrapper.h:91-110).  A signed distance function is not
//  twice differentiable at the medial axis, which cut cells of a thin feature
//  do contain, so that estimate diverges and the surface weights blow up --
//  measured here as an area of 6.0, 12.7, 17.0 and then inf over refinement
//  levels 0..3, while the volume from the same level set converged to 5e-4.
// =============================================================================
static void brepSurfaceQuadrature(const gsMultiPatch<real_t> & brep,
                                  const std::vector<real_t> & sgn,
                                  real_t cellSize, index_t nGauss,
                                  gsMatrix<real_t> & pts, gsVector<real_t> & wts,
                                  gsMatrix<real_t> & nrm)
{
    std::vector<real_t> P, N, W;

    gsVector<index_t> nn(2); nn << nGauss, nGauss;
    gsGaussRule<real_t> gr(nn);
    gsMatrix<real_t> qp, vals, ders;
    gsVector<real_t> qw;

    for (size_t p = 0; p != brep.nPatches(); ++p)
    {
        const gsGeometry<real_t> & g = brep.patch(p);
        const gsBasis<real_t>    & b = g.basis();
        const gsMatrix<real_t>   & C = g.coefs();
        gsMatrix<real_t>  centre(2,1);
        gsMatrix<index_t> act;

        for (auto & elem : b.domain()->allElements())
        {
            const gsVector<real_t> elo = elem.lowerCorner();
            const gsVector<real_t> ehi = elem.upperCorner();

            // Physical size of this element, from the control-point hull.
            centre.col(0) = (real_t)0.5 * (elo + ehi);
            b.active_into(centre, act);
            gsVector<real_t,3> clo, chi;
            clo.setConstant( std::numeric_limits<real_t>::max());
            chi.setConstant(-std::numeric_limits<real_t>::max());
            for (index_t i = 0; i != act.rows(); ++i)
                for (short_t k = 0; k != 3; ++k)
                {
                    clo[k] = std::min(clo[k], C(act(i,0), k));
                    chi[k] = std::max(chi[k], C(act(i,0), k));
                }

            // Split so that a sub-box image stays well below one cell.
            const real_t diam = (chi - clo).norm();
            index_t nSub = (index_t)std::ceil(2 * diam / std::max(cellSize, (real_t)1e-12));
            nSub = std::min(std::max(nSub, (index_t)1), (index_t)32);

            for (index_t si = 0; si != nSub; ++si)
            for (index_t sj = 0; sj != nSub; ++sj)
            {
                gsVector<real_t> slo(2), shi(2);
                slo[0] = elo[0] + (ehi[0]-elo[0]) * (real_t)si     / (real_t)nSub;
                shi[0] = elo[0] + (ehi[0]-elo[0]) * (real_t)(si+1) / (real_t)nSub;
                slo[1] = elo[1] + (ehi[1]-elo[1]) * (real_t)sj     / (real_t)nSub;
                shi[1] = elo[1] + (ehi[1]-elo[1]) * (real_t)(sj+1) / (real_t)nSub;

                gr.mapTo(slo, shi, qp, qw);
                g.eval_into (qp, vals);
                g.deriv_into(qp, ders);

                for (index_t i = 0; i != qp.cols(); ++i)
                {
                    gsVector<real_t,3> xu, xv;
                    for (short_t k = 0; k != 3; ++k)
                    { xu[k] = ders(2*k,i); xv[k] = ders(2*k+1,i); }
                    gsVector<real_t,3> n = sgn[p] * xu.cross(xv);
                    const real_t nl = n.norm();
                    if (nl <= 0) continue;                  // degenerate point
                    n /= nl;
                    for (short_t k = 0; k != 3; ++k)
                    { P.push_back(vals(k,i)); N.push_back(n[k]); }
                    W.push_back(qw[i] * nl);                // dS = |x_u x x_v| du dv
                }
            }
        }
    }

    const index_t M = (index_t)W.size();
    pts.resize(3,M); nrm.resize(3,M); wts.resize(M);
    for (index_t i = 0; i != M; ++i)
    {
        for (short_t k = 0; k != 3; ++k)
        { pts(k,i) = P[3*i+k]; nrm(k,i) = N[3*i+k]; }
        wts[i] = W[i];
    }
}

// Append physical points + scalar values as a 4 x N point cloud.
static void appendCloud(gsMatrix<real_t> & cloud, const gsMatrix<real_t> & phys,
                        const gsVector<real_t> & scalar)
{
    if (phys.cols() == 0) return;
    const index_t c = cloud.cols();
    cloud.conservativeResize(4, c + phys.cols());
    cloud.block(0, c, 3, phys.cols()) = phys;
    cloud.row(3).segment(c, phys.cols()) = scalar.transpose();
}

// =============================================================================
int main(int argc, char * argv[])
{
    std::string filename  = "breps/3D/duck_BRep.xml";
    std::string out       = "output_brep_quadrature";
    index_t     numRefine = 3;
    real_t      fill      = 0.9;
    index_t     quadDepth = 0;
    std::string indicator = "integralChange";
    real_t      indicatorTol = 1e-2;
    bool        plot      = false;
    bool        noSurface = false;
    index_t     selfTest  = 0;
    bool        momentFit = false;
    real_t      alpha     = (real_t)0.0;

    gsCmdLine cmd("Immersed cut-cell quadrature from a spline BREP (gsMultiPatch shell).");
    cmd.addString("f", "file",          "BREP file (gsMultiPatch, parDim 2, geoDim 3)", filename);
    cmd.addString("o", "output",        "Output folder for ParaView files",   out);
    cmd.addInt   ("r", "uniformRefine", "Number of uniform refinement steps",  numRefine);
    cmd.addReal  ("",  "fill",          "Fill fraction of [0,1]^3 (0..1)",     fill);
    cmd.addInt   ("",  "quadDepth",     "Adaptive octree depth for cut cells", quadDepth);
    cmd.addString("",  "indicator",     "Adaptive indicator: uniform, fallback or "
                                        "integralChange",                      indicator);
    cmd.addReal  ("",  "indicatorTol",  "integralChange acceptance tolerance", indicatorTol);
    cmd.addSwitch("noSurface", "Skip the surface (boundary) quadrature",       noSurface);
    cmd.addInt   ("",  "selfTest",      "Check phi against gsMultiPatch::closestPointTo "
                                        "at this many random points (0 = off)", selfTest);
    cmd.addSwitch("plot", "Create ParaView output",                            plot);
    cmd.addSwitch("momentFit", "Use moment-fitting compression for the volume "
                               "rule (surface rule unaffected)",               momentFit);
    cmd.addReal  ("",  "alpha",         "Fictitious-domain weight of the outside "
                                        "material in the moment-fitting rule "
                                        "(0 = drop it; only 0 is supported)",  alpha);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Capability regression of the move to the core gsMomentRule: the FCM
    // blend w <- (1-alpha)w + alpha*w^G needs the analytic full-cell Gauss
    // weight w^G, which the geometry-free compressor does not have.
    GISMO_ENSURE(alpha == (real_t)0,
                 "--alpha (fictitious-domain blend) is not available with the "
                 "core gsMomentRule; use --alpha 0.");

    // -------------------------------------------------------------------------
    // 1. Load the spline BREP
    // -------------------------------------------------------------------------
    gsMultiPatch<real_t> brep;
    gsReadFile<real_t>(filename, brep);
    GISMO_ENSURE(brep.nPatches() > 0, "Could not read a gsMultiPatch from '" << filename << "'.");
    GISMO_ENSURE(brep.domainDim() == 2 && brep.targetDim() == 3,
                 "Expected a surface BREP (parDim 2, geoDim 3), got parDim "
                 << brep.domainDim() << ", geoDim " << brep.targetDim() << ".");
    gsInfo << "Spline BREP '" << filename << "':\n";
    gsInfo << "  patches / interfaces (in file): " << brep.nPatches()
           << " / " << brep.nInterfaces() << "\n";
    if (brep.nBoundary() != 0)
        gsWarn << "  The stored topology reports " << brep.nBoundary()
               << " free boundary side(s). Watertightness is re-established "
                  "geometrically below; if that fails the file is genuinely open.\n";

    // -------------------------------------------------------------------------
    // 2. Orientation fix-up + reference volume/area, on the ORIGINAL coordinates
    // -------------------------------------------------------------------------
    real_t refVolume = 0, refArea = 0;
    std::vector<real_t> sgn = fixOrientation(brep, refVolume, refArea);

    // -------------------------------------------------------------------------
    // 3. Rescale into [0,1]^3 so the background box is the unit cube.
    //    NOTE two different factors: volumes scale with scale^3, areas with
    //    scale^2.  Reference values above are in the ORIGINAL coordinates.
    // -------------------------------------------------------------------------
    gsMatrix<real_t> obox = brep.patch(0).coefs().colwise().minCoeff().transpose();
    obox.conservativeResize(3,2);
    obox.col(1) = brep.patch(0).coefs().colwise().maxCoeff().transpose();
    for (size_t p = 0; p != brep.nPatches(); ++p)
    {
        obox.col(0) = obox.col(0).cwiseMin(brep.patch(p).coefs().colwise().minCoeff().transpose());
        obox.col(1) = obox.col(1).cwiseMax(brep.patch(p).coefs().colwise().maxCoeff().transpose());
    }
    const gsVector<real_t> centre  = (real_t)0.5 * (obox.col(0) + obox.col(1));
    const real_t           extent  = (obox.col(1) - obox.col(0)).maxCoeff();
    const real_t           scale   = fill / extent;
    const real_t           invS3   = (real_t)1 / (scale*scale*scale);
    const real_t           invS2   = (real_t)1 / (scale*scale);

    for (size_t p = 0; p != brep.nPatches(); ++p)
    {
        gsMatrix<real_t> & C = brep.patch(p).coefs();
        C = ((C.rowwise() - centre.transpose()) * scale).array() + (real_t)0.5;
    }
    gsInfo << "  rescaled into [0,1]^3, scale  : " << scale << "\n";

    // -------------------------------------------------------------------------
    // 4. Bezier-element boxes + the level set
    // -------------------------------------------------------------------------
    std::vector<BezBox> boxes = makeBezBoxes(brep);
    gsInfo << "  Bezier-element boxes          : " << boxes.size() << "\n";

    gsMatrix<real_t> bbox(3,2); bbox.col(0).setZero(); bbox.col(1).setOnes();
    gsBRepSignedDist<real_t> phi(brep, boxes, sgn, bbox);

    // ---- sign sanity, before any quadrature runs ----------------------------
    {
        gsMatrix<real_t> corners(3,8), cval;
        index_t c = 0;
        for (int i = 0; i != 2; ++i)
        for (int j = 0; j != 2; ++j)
        for (int k = 0; k != 2; ++k, ++c)
            corners.col(c) << (real_t)i, (real_t)j, (real_t)k;
        phi.eval_into(corners, cval);
        gsInfo << "  phi at [0,1]^3 corners        : min "
               << cval.minCoeff() << " (must be > 0)\n";
        GISMO_ENSURE(cval.minCoeff() > 0,
                     "phi is negative at a corner of the background box: the sign "
                     "convention or the patch orientation is wrong.");
    }

    // ---- optional independent check of magnitude AND sign -------------------
    // Two references, neither of which shares any code with gsBRepSignedDist:
    //   magnitude : distance to a dense Gauss-point sampling of the surface
    //               (an upper bound, tight to the sample spacing);
    //   sign      : the generalised winding number
    //               w(x) = 1/(4 pi) oint (y-x).n / |y-x|^3 dS,
    //               ~1 inside and ~0 outside, evaluated by the same samples.
    // gsMultiPatch::closestPointTo is deliberately NOT used: its per-patch
    // Newton frequently fails to converge here (residuals ~1e-1), so it is the
    // less trustworthy of the two.  It also returns ||x-y||/sqrt(2) rather than
    // the Euclidean distance, since gsSquaredDistance is 1/2||.||^2
    // (gsGeometry.hpp:339-357; note the sqrt(2*..) correction at line 373).
    if (selfTest > 0)
    {
        gsMatrix<real_t> sPts, sNrm;      // 3 x M samples and unit normals
        gsVector<real_t> sW;              // area weights
        {
            std::vector<real_t> px, py, pz, nx, ny, nz, wv;
            gsVector<index_t> nn(2); nn << 6, 6;
            gsGaussRule<real_t> gr(nn);
            gsMatrix<real_t> qp, vals, ders;  gsVector<real_t> qw;
            for (size_t p = 0; p != brep.nPatches(); ++p)
            {
                const gsGeometry<real_t> & g = brep.patch(p);
                for (auto & elem : g.basis().domain()->allElements())
                {
                    gr.mapTo(elem.lowerCorner(), elem.upperCorner(), qp, qw);
                    g.eval_into(qp, vals);  g.deriv_into(qp, ders);
                    for (index_t i = 0; i != qp.cols(); ++i)
                    {
                        gsVector<real_t,3> xu, xv;
                        for (short_t k = 0; k != 3; ++k)
                        { xu[k] = ders(2*k,i); xv[k] = ders(2*k+1,i); }
                        gsVector<real_t,3> n = sgn[p] * xu.cross(xv);
                        const real_t nl = n.norm();
                        if (nl <= 0) continue;
                        px.push_back(vals(0,i)); py.push_back(vals(1,i)); pz.push_back(vals(2,i));
                        n /= nl;
                        nx.push_back(n[0]); ny.push_back(n[1]); nz.push_back(n[2]);
                        wv.push_back(qw[i] * nl);
                    }
                }
            }
            const index_t M = (index_t)px.size();
            sPts.resize(3,M); sNrm.resize(3,M); sW.resize(M);
            for (index_t i = 0; i != M; ++i)
            {
                sPts(0,i)=px[i]; sPts(1,i)=py[i]; sPts(2,i)=pz[i];
                sNrm(0,i)=nx[i]; sNrm(1,i)=ny[i]; sNrm(2,i)=nz[i];
                sW[i]=wv[i];
            }
        }
        gsInfo << "\nSelf-test (" << selfTest << " points, "
               << sPts.cols() << " surface samples):\n";

        real_t  maxAbsErr = 0, maxGradErr = 0, maxGradNorm = 0, minGradNorm = 1e30;
        index_t nSignBad = 0, nSignChecked = 0;
        const real_t pi4 = (real_t)(4.0*3.14159265358979323846);

        for (index_t t = 0; t != selfTest; ++t)
        {
            gsVector<real_t> x(3);
            for (short_t k = 0; k != 3; ++k)
                x[k] = (real_t)((t*7919 + k*104729 + 13) % 10007) / (real_t)10007;

            gsMatrix<real_t> xm(3,1), pv; xm.col(0) = x;
            phi.eval_into(xm, pv);

            real_t dRef = std::numeric_limits<real_t>::max(), w = 0;
            for (index_t i = 0; i != sPts.cols(); ++i)
            {
                const gsVector<real_t,3> d = sPts.col(i) - x;
                const real_t r = d.norm();
                dRef = std::min(dRef, r);
                w += sW[i] * d.dot(sNrm.col(i)) / (r*r*r);
            }
            w /= pi4;

            // |phi| must not exceed the sampled distance (which over-estimates).
            maxAbsErr = std::max(maxAbsErr, math::abs(pv(0,0)) - dRef);

            // The winding integral is only trustworthy away from the surface.
            if (dRef > 0.02)
            {
                ++nSignChecked;
                const bool insideRef = (w > 0.5);
                const bool insidePhi = (pv(0,0) < 0);
                if (insideRef != insidePhi) ++nSignBad;
            }

            // Analytic gradient vs central differences.  |grad phi| must be 1
            // for a true distance function; algoim's surface rule divides by a
            // gradient component, so an inconsistent gradient shows up there
            // long before it shows up in the volume.
            gsMatrix<real_t> gA;
            phi.deriv_into(xm, gA);
            gsVector<real_t,3> gFD;
            const real_t h = 1e-6;
            for (short_t k = 0; k != 3; ++k)
            {
                gsMatrix<real_t> xp = xm, xmm = xm, vp, vm;
                xp(k,0) += h; xmm(k,0) -= h;
                phi.eval_into(xp, vp); phi.eval_into(xmm, vm);
                gFD[k] = (vp(0,0) - vm(0,0)) / (2*h);
            }
            maxGradErr  = std::max(maxGradErr, (gA.col(0) - gFD).norm());
            maxGradNorm = std::max(maxGradNorm, (real_t)gA.col(0).norm());
            minGradNorm = std::min(minGradNorm, (real_t)gA.col(0).norm());
        }
        gsInfo << "  max ( |phi| - sampledDist )   : " << maxAbsErr
               << "   (must be <= ~sample spacing)\n";
        gsInfo << "  sign mismatches vs winding no.: " << nSignBad
               << " / " << nSignChecked << "\n";
        gsInfo << "  max |gradAnalytic - gradFD|   : " << maxGradErr << "\n";
        gsInfo << "  |grad phi| range              : [" << minGradNorm
               << ", " << maxGradNorm << "]   (must be ~1)\n";
    }

    // -------------------------------------------------------------------------
    // 5. Background box and output setup
    // -------------------------------------------------------------------------
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsMultiBasis<> dbasis(mp, true);
    gsFileManager::mkdir(out);

    gsInfo << "\nReference (original coordinates):  volume " << std::setprecision(8)
           << refVolume << ",  area " << refArea << "\n";
    gsInfo << "Adaptive octree depth : " << quadDepth << "\n\n";
    gsInfo << std::setw(7)  << "refine" << std::setw(9)  << "cells"
           << std::setw(16) << "vol(phys)"  << std::setw(12) << "volErr"
           << std::setw(16) << "area(phys)" << std::setw(12) << "areaErr"
           << std::setw(11) << "volPts" << std::setw(11) << "cutCells"
           << std::setw(11) << "srfPts" << std::setw(9) << "stray" << "\n";

    std::unique_ptr<gsParaviewCollection> colInterior, colCut, colSurf, colAll, colBg;
    if (plot)
    {
        colInterior.reset(new gsParaviewCollection(out + "/points_interior/interior"));
        colCut     .reset(new gsParaviewCollection(out + "/points_cutcells/cutcells"));
        colSurf    .reset(new gsParaviewCollection(out + "/points_surface/surface"));
        colAll     .reset(new gsParaviewCollection(out + "/points_all/all"));
        colBg      .reset(new gsParaviewCollection(out + "/background/background"));
    }

    std::vector<index_t> cellsHist, quadHist, bdryHist, subBoxHist,
                         cutQuadHist, intSubHist;
    std::vector<real_t>  volPhyHist, areaPhyHist;

    // -------------------------------------------------------------------------
    // 6. Refinement loop
    // -------------------------------------------------------------------------
    for (index_t r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();
        gsTensorBSplineBasis<3,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<3,real_t>*>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        gsOptionList adaptiveOptions = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
        adaptiveOptions.setInt   ("maxDepth",         quadDepth);
        adaptiveOptions.setInt   ("nFallback",        4);
        adaptiveOptions.setString("indicator",        indicator);
        adaptiveOptions.setReal  ("indicatorTol",     indicatorTol);
        // A true signed distance is 1-Lipschitz, so the default 1.0 is exact here.
        gsAlgoimAdaptiveRule<real_t> adaptiveRule(phi, *tbsPtr, adaptiveOptions);

        // Optional replacement of the VOLUME rule only: the same adaptive rule
        // runs underneath (identical maxDepth/nFallback/indicator/indicatorTol),
        // but the cut-cell point cloud is compressed onto the tensor Gauss grid
        // of the cell.  Rebuilt every refinement, so the stats printed below
        // refer to this level only.  The surface quadrature below comes from
        // the BREP directly and never passes through this rule.
        // Note: total mass is preserved identically by the partition of unity
        // of the Lagrange basis, so reproducing the volume is a wiring check,
        // not an accuracy check.
        // The same boxClassBRep closure is installed on the OWNED adaptive rule
        // below (the compressor carries no classification of its own), so its
        // subdivision classifies sub-boxes exactly as the adaptive path does
        // and the two agree on the underlying point count.
        //
        // The compressor OWNS the adaptive rule (the gsQuadRule hierarchy has
        // no clone(), so ownership transfer is the only non-slicing option).
        // A non-const raw observer is kept to install the per-cell classifier;
        // it is re-assigned together with the rule at every refinement.
        gsMomentRule<real_t>::uPtr     momentFitRule;
        gsAlgoimAdaptiveRule<real_t> * mfUnderlying = nullptr;
        if (momentFit)
        {
            // Output points per direction, exactly as the former moment-fitting
            // rule computed them: n = round(quA*maxDegree) + quB, clamped to 1.
            const long nRaw = std::lround(adaptiveOptions.askReal("quA", 1.0)
                                          * static_cast<double>(tbsPtr->maxDegree()))
                            + static_cast<long>(adaptiveOptions.askInt("quB", 1));
            const index_t mfOrder1d = nRaw > 0 ? static_cast<index_t>(nRaw) : 1;

            gsAlgoimAdaptiveRule<real_t> * u =
                new gsAlgoimAdaptiveRule<real_t>(phi, *tbsPtr, adaptiveOptions);
            mfUnderlying  = u;
            momentFitRule = gsMomentRule<real_t>::make(
                                gsQuadRule<real_t>::uPtr(u),
                                gsVector<index_t>::Constant(3, mfOrder1d));
        }

        const std::vector<real_t> bx = tbsPtr->knots(0).breaks();
        const std::vector<real_t> by = tbsPtr->knots(1).breaks();
        const std::vector<real_t> bz = tbsPtr->knots(2).breaks();

        // ---- classify every background element ------------------------------
        real_t volInterior = 0;
        gsMatrix<real_t> intCtr(3,0);
        gsVector<real_t> intVol;
        std::vector<gsVector<real_t> >  cutLo, cutHi;
        std::vector<std::vector<int> >  cutCand;

        for (size_t i = 0; i+1 < bx.size(); ++i)
        for (size_t j = 0; j+1 < by.size(); ++j)
        for (size_t k = 0; k+1 < bz.size(); ++k)
        {
            gsVector<real_t> lo(3), hi(3);
            lo << bx[i], by[j], bz[k];
            hi << bx[i+1], by[j+1], bz[k+1];

            // AABB pre-filter: only Bezier elements whose hull box meets the cell.
            std::vector<int> cIdx;
            for (size_t b = 0; b != boxes.size(); ++b)
                if (boxes[b].overlaps(lo, hi)) cIdx.push_back((int)b);

            const int cls = boxClassBRep(boxes, cIdx, phi, lo, hi);
            if (cls > 0) continue;                     // outside
            if (cls < 0)                               // inside: exact box volume
            {
                const real_t v = (hi - lo).prod();
                volInterior += v;
                if (plot)
                {
                    const index_t c = intCtr.cols();
                    intCtr.conservativeResize(3, c + 1);
                    intCtr.col(c) = (real_t)0.5 * (lo + hi);
                    intVol.conservativeResize(c + 1);
                    intVol[c] = v;
                }
                continue;
            }
            cutLo.push_back(lo);
            cutHi.push_back(hi);
            cutCand.push_back(give(cIdx));
        }
        const int nCut = (int)cutLo.size();

        // ---- quadrature on the cut cells ------------------------------------
        gsMatrix<real_t> quadIn(4,0), quadCut(4,0), quadSurf(4,0), quadAll(4,0);
        real_t  volCut = 0, areaSum = 0;
        index_t nQuad = 0, nBdry = 0, nInteriorSub = 0;

        // The moment-fitting rule keeps mutable scratch state (see its header),
        // so it must not be shared between threads: the parallel region is
        // disabled in that mode.  (This build has GISMO_WITH_OPENMP=OFF, so the
        // clause is a guard, not something exercised here.)
#       pragma omp parallel for schedule(dynamic) if(!momentFit) \
                reduction(+:volCut,nQuad,nBdry,nInteriorSub)
        for (int i = 0; i < nCut; ++i)
        {
            gsMatrix<real_t> insCtr(3,0), cutP(3,0);
            gsVector<real_t> insVol, cutW;

            // The geometry-aware classifier of this cell, shared by both
            // paths: cutCand[i] is complete for every sub-box of cell i (a
            // Bezier box meeting a sub-box meets the cell too, so it is in
            // cIdx), which is what makes it valid to hand the same closure to
            // the adaptive subdivision.
            // i is captured by value because setClassifier() stores a copy of
            // the closure that outlives the loop -- a reference would dangle.
            auto cellClass = [&, i](const gsVector<real_t> & childLo,
                                    const gsVector<real_t> & childHi) -> int
            {
                return boxClassBRep(boxes, cutCand[i], phi, childLo, childHi);
            };

            if (momentFit)
            {
                // One compressed cloud per cell; there is no separated
                // interior output, so everything is counted as cut-cell
                // points (volPts == cutRulePts, octreePts == 0 in this mode).
                // The classifier is installed per cell (it captures i) on the
                // OWNED adaptive rule -- the compressor has no classification
                // of its own -- which is safe only because the parallel region
                // above is disabled in this mode.
                mfUnderlying->setClassifier(cellClass);
                momentFitRule->mapTo(cutLo[i], cutHi[i], cutP, cutW);
            }
            else
            {
                adaptiveRule.mapToSeparated(
                    cutLo[i], cutHi[i], insCtr, insVol, cutP, cutW, cellClass);
            }

            volCut       += insVol.sum() + cutW.sum();
            nQuad        += insCtr.cols() + cutP.cols();
            nBdry        += cutP.cols();
            nInteriorSub += insCtr.cols();

            if (plot)
            {
                gsMatrix<real_t> physIn, physCut;
                mp.patch(0).eval_into(insCtr, physIn);
                mp.patch(0).eval_into(cutP,   physCut);
#               pragma omp critical
                {
                    appendCloud(quadIn,   physIn,   insVol);
                    appendCloud(quadCut,  physCut,  cutW);
                    appendCloud(quadAll,  physIn,   insVol);
                    appendCloud(quadAll,  physCut,  cutW);
                }
            }
        }

        // ---- surface quadrature, taken directly from the BREP ---------------
        index_t nSurfPts = 0, nSurfStray = 0;
        if (!noSurface)
        {
            const real_t h = (real_t)1 / (real_t)tbsPtr->component(0).numElements();
            gsMatrix<real_t> sP, sN;
            gsVector<real_t> sW;
            brepSurfaceQuadrature(brep, sgn, h, 5, sP, sW, sN);
            areaSum  = sW.sum();
            nSurfPts = sP.cols();

            // Cross-check: every boundary point must fall in a cell the
            // classifier called "cut".  A stray point means the classification
            // missed a piece of the interface, which would silently lose volume.
            const index_t nc = tbsPtr->component(0).numElements();
            std::set<index_t> cutSet;
            for (int c = 0; c != nCut; ++c)
            {
                const index_t ci = (index_t)std::floor(0.5*(cutLo[c][0]+cutHi[c][0]) * nc);
                const index_t cj = (index_t)std::floor(0.5*(cutLo[c][1]+cutHi[c][1]) * nc);
                const index_t ck = (index_t)std::floor(0.5*(cutLo[c][2]+cutHi[c][2]) * nc);
                cutSet.insert((ci*nc + cj)*nc + ck);
            }
            for (index_t i = 0; i != sP.cols(); ++i)
            {
                const index_t ci = std::min(std::max((index_t)std::floor(sP(0,i)*nc), (index_t)0), nc-1);
                const index_t cj = std::min(std::max((index_t)std::floor(sP(1,i)*nc), (index_t)0), nc-1);
                const index_t ck = std::min(std::max((index_t)std::floor(sP(2,i)*nc), (index_t)0), nc-1);
                if (cutSet.find((ci*nc + cj)*nc + ck) == cutSet.end()) ++nSurfStray;
            }
            if (nSurfStray > 0)
                gsWarn << "  r=" << r << ": " << nSurfStray << " of " << sP.cols()
                       << " boundary points lie outside every cell classified as cut.\n";

            if (plot) appendCloud(quadSurf, sP, sW);
        }

        const real_t  volParam = volInterior + volCut;
        const real_t  volPhys  = volParam * invS3;
        const real_t  areaPhys = areaSum   * invS2;
        const index_t cellsPerDim = tbsPtr->component(0).numElements();

        gsInfo << std::setw(7)  << r << std::setw(9) << cellsPerDim
               << std::setw(16) << std::fixed << std::setprecision(6) << volPhys
               << std::setw(12) << std::scientific << std::setprecision(2)
               << math::abs(volPhys - refVolume)
               << std::setw(16) << std::fixed << std::setprecision(6) << areaPhys
               << std::setw(12) << std::scientific << std::setprecision(2)
               << math::abs(areaPhys - refArea)
               << std::setw(11) << nQuad << std::setw(11) << nCut
               << std::setw(11) << nSurfPts << std::setw(9) << nSurfStray << "\n";

        // -- Moment-fitting compression statistics (this refinement only: the
        //    rule is rebuilt every r).  outQPs/minW/negW describe the
        //    compressed cells only; cells whose underlying cloud already fits
        //    in the output grid are passed through unchanged (passThrough) and
        //    contribute passQPs, not outQPs.  The honest compression ratio is
        //    underQPs/emitted, with emitted = outQPs + passQPs.
        if (momentFit)
        {
            const gsMomentRule<real_t>::Stats & mfs = momentFitRule->stats();
            gsInfo << "momfit: outQPs=" << mfs.nOutputQPs
                   << " passQPs=" << mfs.nPassThroughQPs
                   << " emitted=" << (mfs.nOutputQPs + mfs.nPassThroughQPs)
                   << " underQPs=" << mfs.nUnderlyingQPs
                   << " passThrough=" << mfs.nPassThroughElements
                   << " minW=" << std::scientific << std::setprecision(3) << mfs.minWeight
                   << " negW=" << mfs.nNegativeWeights << "\n";
        }

        cellsHist  .push_back(cellsPerDim);
        quadHist   .push_back(nQuad);
        bdryHist   .push_back(nSurfPts);
        subBoxHist .push_back(nCut);
        cutQuadHist.push_back(nBdry);
        intSubHist .push_back(nInteriorSub);
        volPhyHist .push_back(volPhys);
        areaPhyHist.push_back(areaPhys);

        // ---- ParaView -------------------------------------------------------
        if (plot)
        {
            const std::string rs = std::to_string(r);
            gsMatrix<real_t> physInt;
            mp.patch(0).eval_into(intCtr, physInt);
            appendCloud(quadIn,  physInt, intVol);
            appendCloud(quadAll, physInt, intVol);

            if (quadIn.cols())
            { gsWriteParaviewPoints(quadIn,   out+"/points_interior/interior_r"+rs);
              colInterior->addPart("interior_r"+rs+".vtp", r); }
            if (quadCut.cols())
            { gsWriteParaviewPoints(quadCut,  out+"/points_cutcells/cutcells_r"+rs);
              colCut->addPart("cutcells_r"+rs+".vtp", r); }
            if (quadSurf.cols())
            { gsWriteParaviewPoints(quadSurf, out+"/points_surface/surface_r"+rs);
              colSurf->addPart("surface_r"+rs+".vtp", r); }
            if (quadAll.cols())
            { gsWriteParaviewPoints(quadAll,  out+"/points_all/all_r"+rs);
              colAll->addPart("all_r"+rs+".vtp", r); }

            gsMesh<real_t> bgMesh(dbasis.basis(0));
            gsWriteParaview(bgMesh, out + "/background/background_r" + rs);
            colBg->addPart("background_r"+rs+".vtp", r);
        }
    }

    // Full-precision values of the finest level: the table above is rounded to
    // 6 digits, which cannot show whether two runs (e.g. with and without
    // --momentFit) agree to round-off.
    if (!volPhyHist.empty())
        gsInfo << "\nfinal (r=" << numRefine << ")"
               << (momentFit ? " [momentFit]" : "")
               << ": vol(phys)=" << std::scientific << std::setprecision(16)
               << volPhyHist.back()
               << "  area(phys)=" << areaPhyHist.back() << "\n";

    // -------------------------------------------------------------------------
    // 7. Results file
    // -------------------------------------------------------------------------
    gsFileManager::mkdir(out + "/results");
    {
        std::ofstream txt((out + "/results/volume.txt").c_str());
        txt << "# Immersed quadrature from the spline BREP '" << filename << "'\n";
        txt << "# patches: " << brep.nPatches() << ", interfaces: " << brep.nInterfaces()
            << ", Bezier-element boxes: " << boxes.size() << "\n";
        txt << "# adaptive octree depth (quadDepth): " << quadDepth << "\n";
        txt << "# rescaling factor (param/phys): " << scale
            << "   (volume ~ scale^3, area ~ scale^2)\n";
        txt << "# reference volume (1/3 oint x.n dS): " << std::setprecision(10) << refVolume << "\n";
        txt << "# reference area   (oint |n| dS)    : " << std::setprecision(10) << refArea   << "\n\n";
        txt << std::setw(8) << "refine" << std::setw(10) << "cells/dim"
            << std::setw(20) << "vol(phys)"  << std::setw(14) << "volErr"
            << std::setw(20) << "area(phys)" << std::setw(14) << "areaErr"
            << std::setw(12) << "volPts" << std::setw(12) << "srfPts"
            << std::setw(12) << "cutCells" << std::setw(12) << "cutRulePts"
            << std::setw(12) << "octreePts" << "\n";
        for (size_t k = 0; k != volPhyHist.size(); ++k)
            txt << std::setw(8) << k << std::setw(10) << cellsHist[k]
                << std::setw(20) << std::fixed << std::setprecision(8) << volPhyHist[k]
                << std::setw(14) << std::scientific << std::setprecision(3)
                << std::abs(volPhyHist[k] - refVolume)
                << std::setw(20) << std::fixed << std::setprecision(8) << areaPhyHist[k]
                << std::setw(14) << std::scientific << std::setprecision(3)
                << std::abs(areaPhyHist[k] - refArea)
                << std::setw(12) << quadHist[k] << std::setw(12) << bdryHist[k]
                << std::setw(12) << subBoxHist[k] << std::setw(12) << cutQuadHist[k]
                << std::setw(12) << intSubHist[k] << "\n";
    }
    gsInfo << "\nResults written to: " << out << "/results/volume.txt\n";

    // -------------------------------------------------------------------------
    // 8. Geometry + level-set output
    // -------------------------------------------------------------------------
    if (plot)
    {
        gsWriteParaview(brep, out + "/brep", 1000);
        gsField<real_t> lsField(mp, phi, true);
        gsWriteParaview(lsField, out + "/levelset", 8000);

        colInterior->save(); colCut->save(); colAll->save(); colBg->save();
        if (!noSurface) colSurf->save();
        gsInfo << "ParaView files written to: " << out << "\n";
    }

    return EXIT_SUCCESS;
}
