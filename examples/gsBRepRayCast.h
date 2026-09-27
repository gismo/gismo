/** @file gsBRepRayCast.h

    @brief Ray-casting primitives on a closed, untrimmed spline BRep
    (gsMultiPatch<T>, parDim 2, geoDim 3, B-spline or NURBS patches).

    Everything here works on the exact Bezier-element decomposition of the
    BRep and on the exact tensor-Bernstein coefficients of the (unnormalized)
    surface normal, never on a facetted approximation. Four primitives are
    provided:

    - P1 (::extractElements, ::splitElement, ::BezElement): the Bezier
      elements of the BRep, each carrying its homogeneous control net, its uv
      box in the original patch parameters, its control-point AABB, and its
      normal certificate (P2).
    - P2 (::normalNumerator): the normal numerator N = w^4 (S_u x S_v),
      certified sign-wise via its own Bernstein coefficients.
    - P3 (::hitElement, ::castLine, ::windingOK, ::clipHits): the exact
      axis-aligned ray/element intersection and the full line cast through a
      set of elements, with winding-number verification.
    - P4 (::elementEdges, ::curvePlaneCrossings): element boundary
      iso-curves and their crossings with an axis-aligned plane.

    Convention shared by P2 and P3: for a line along e_k, the transverse
    coordinate indices are i = (k+1)%3, j = (k+2)%3 (cyclic), so that the 2x2
    Jacobian [[dS_i/du, dS_i/dv],[dS_j/du, dS_j/dv]] has determinant exactly
    (S_u x S_v)_k, and the ray parameter t is the coordinate x_k itself.

    The inclusion and side decisions below are exact predicates with no
    tolerance. Root acceptance (residual <= 1e-13*L), the duplicate-hit
    merge (1e-12*span) and the certificate margin (1e-10*Nscale) are
    tolerance-based:
    - a line with transverse coordinates (y_i, y_j) meets an element iff
      lo_i <= y_i < up_i && lo_j <= y_j < up_j (half-open, so a line lying ON
      a coplanar face's plane skips that face);
    - a hit belongs to a clip range (a, b] iff a < t && t <= b, so t == a
      goes to the lower box;
    - the plane-side predicate is above(x) := x_k > X (x_k == X is not
      above).

    Exact-constant rule: evaluating a B-spline component whose control
    values are all bitwise equal returns that value to +-1 ulp, not exactly.
    So whenever every projected control coordinate of an element or curve is
    bitwise equal to one value c, this file uses c directly and never
    evaluates: P3's hit parameter t is c when component k is constant on the
    element, and P4's plane-crossing count is read off the (unevaluated)
    homogeneous control coefficients, so a curve constant in x_k never
    reports a crossing. Curve endpoints are likewise taken from the first and
    last control point (endpoint interpolation), never evaluated, so that
    neighbouring curves share bitwise-identical corners.

    Hand-off to whatever consumes these primitives next (octree/grid
    traversal, a quadrature rule): the driver's background box (Geometry::bg;
    its largest edge is the length scale L passed to ::castLine and
    ::hitElement) is the authoritative bound on where a caller may ask for a
    ray -- for the
    axis-aligned cube geometry that bound coincides exactly with element
    faces only when it is fixed to [0,1]^3 AND any background lattice
    resolution n0 (refined by r uniform steps) satisfies n0*2^r == 0 (mod 4);
    outside that relation cube faces no longer land on lattice planes.
    ::BezElement exposes both certWeak and certStrict because a consumer
    that only trusts certStrict will find every one of the sphere's 8
    elements uncertified in every direction (N === 0 identically on the
    pole/meridian edges), even though ::castLine works correctly on them
    under certWeak alone.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gismo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

namespace gismo {
namespace brc {

/**
    @brief One Bezier element of a BRep patch, on the 4-column homogeneous
    control net X = (w*S_x, w*S_y, w*S_z, w).

    \a hom is built by ::extractElements (via splitAtMult(1) on the whole
    patch cast to a 4-column projective gsTensorBSpline<2,T>) or by
    ::splitElement (via uniformSplit), so its knot values are the ORIGINAL
    patch parameters and \a uv == hom.support().

    \a lo, \a up bound the projected control points x = X/w (valid since
    w > 0 is enforced at extraction), a convex-hull AABB of the element's
    image.

    \a Ncoef, \a Nmin/\a Nmax/\a Nscale and \a certWeak/\a certStrict are the
    Bernstein certificate of sgn*N described at ::normalNumerator.
*/
template<class T>
struct BezElement
{
    index_t patch;                  ///< patch id in the input gsMultiPatch
    T       sgn;                    ///< outward normal = sgn * S_u x S_v
    bool    rational;                ///< true if any weight of the source patch != 1
    gsTensorBSpline<2,T> hom;       ///< Bezier piece, 4 homogeneous coef columns, ORIGINAL knot values
    T uv[2][2];                     ///< uv[d][0..1] = hom.support()(d,0..1)
    T lo[3], up[3];                 ///< AABB of projected control points x = X/w
    index_t dU, dV;                 ///< bidegree of N (see ::normalNumerator)
    gsMatrix<T> Ncoef;              ///< (dU+1)(dV+1) x 3 Bernstein coefs of sgn*N over the element's uv box, u fastest
    T Nmin[3], Nmax[3], Nscale;     ///< per-component min/max Bernstein coef; Nscale = max|coef| over all 3 components
    int certWeak[3], certStrict[3]; ///< +1 / -1 / 0 per component, see ::normalNumerator
};

/**
    @brief One side of a Bezier element, cut out of an element's Bezier
    piece: the isoparametric curve where direction \a dirFixed equals \a par,
    on the same 4-column homogeneous net as the element (so \a hom is
    exactly (w*S)(par, .) or (w*S)(., par), not the projected curve).
*/
template<class T>
struct EdgeCurve
{
    index_t patch;    ///< patch id of the element this edge came from
    short_t dirFixed; ///< the fixed parametric direction (0 or 1)
    T       par;      ///< the value at which that direction is fixed
    gsBSpline<T> hom; ///< the homogeneous edge curve
};

/// One crossing of an ::EdgeCurve with the axis-aligned plane x_k = X.
template<class T>
struct PlaneCrossing
{
    T s;              ///< parameter, in the curve's own domain, of the crossing
    gsVector<T,3> x;  ///< physical point of the crossing
    int dir;          ///< +1 if x_k goes from <= X to > X as s increases, -1 otherwise
};

/// One intersection of a ray along e_k with a ::BezElement's surface, found by ::hitElement.
template<class T>
struct Hit
{
    index_t patch;      ///< patch id of the hit element
    T u, v;              ///< parameters of the hit point, in the original patch parameters
    T t;                 ///< the ray parameter, i.e. the coordinate x_k of the hit point
    gsVector<T,3> n;     ///< unit outward normal sgn*N/|N| at (u,v)
    bool onBoundary;     ///< true if (u,v) lies within 1e-12*span of the element's uv box edge
    index_t elem;        ///< index into the element vector passed to ::castLine; -1 when produced by ::hitElement directly
};

// Forward declaration: detail::fitNcoef below calls normalNumerator, whose
// definition sits after extractElements/splitElement for readability (it is
// their sibling primitive P2). Argument-dependent lookup alone would not
// find it, since gsTensorBSpline's associated namespace is gismo, not
// gismo::brc, so ordinary unqualified lookup needs this declaration in scope.
template<class T>
void normalNumerator(const gsTensorBSpline<2,T> & hom, T sgn, const gsMatrix<T> & uv, gsMatrix<T> & N);

namespace detail {

/// Returns the 4-column projective gsTensorBSpline<2,T> (w*S_x, w*S_y, w*S_z, w)
/// of a 2->3 B-spline or NURBS patch, and whether any of its weights differs from 1.
template<class T>
gsTensorBSpline<2,T> homogenize(const gsGeometry<T> & patch, bool & rational)
{
    if (const gsTensorNurbs<2,T> * nurbs = dynamic_cast<const gsTensorNurbs<2,T>*>(&patch))
    {
        const gsMatrix<T> & w = nurbs->weights();
        rational = (w.array() != (T)1).any();
        return gsTensorBSpline<2,T>(nurbs->basis().source(),
                                    nurbs->basis().projectiveCoefs(nurbs->coefs()));
    }
    const gsTensorBSpline<2,T> * bs = dynamic_cast<const gsTensorBSpline<2,T>*>(&patch);
    GISMO_ENSURE(bs, "gsBRepRayCast: patch is neither a gsTensorBSpline<2,T> nor a gsTensorNurbs<2,T>.");
    rational = false;
    gsMatrix<T> C(bs->coefs().rows(), 4);
    C.leftCols(3) = bs->coefs();
    C.col(3).setOnes();
    return gsTensorBSpline<2,T>(bs->basis(), give(C));
}

/// Bernstein coefficients of \a sgn*N over the Bezier basis of bidegree
/// (\a dU,\a dV) on the reference square [0,1]^2, mapped affinely to
/// [\a lo0,\a hi0] x [\a lo1,\a hi1]. The reference collocation matrix is
/// invariant under that affine map, so it is factorized once per (dU,dV)
/// and cached across all elements sharing that bidegree.
template<class T>
void fitNcoef(index_t dU, index_t dV, T lo0, T hi0, T lo1, T hi1,
              const gsTensorBSpline<2,T> & hom, T sgn, gsMatrix<T> & Ncoef)
{
    typedef decltype(std::declval<gsMatrix<T>&>().fullPivLu()) LU_t;
    static std::map<std::pair<index_t,index_t>, gsMatrix<T> > s_anchors;
    static std::map<std::pair<index_t,index_t>, LU_t>         s_lu;

    const std::pair<index_t,index_t> key(dU, dV);
    typename std::map<std::pair<index_t,index_t>, gsMatrix<T> >::iterator itA = s_anchors.find(key);
    if (itA == s_anchors.end())
    {
        gsTensorBSplineBasis<2,T> Bref(gsKnotVector<T>(0,1,0,dU+1), gsKnotVector<T>(0,1,0,dV+1));
        gsMatrix<T> A = Bref.anchors();
        gsMatrix<T> C = Bref.collocationMatrix(A);
        itA = s_anchors.insert(std::make_pair(key, give(A))).first;
        s_lu.insert(std::make_pair(key, C.fullPivLu()));
    }

    const gsMatrix<T> & A = itA->second;
    gsMatrix<T> uv(2, A.cols());
    uv.row(0) = (hi0-lo0)*A.row(0); uv.row(0).array() += lo0;
    uv.row(1) = (hi1-lo1)*A.row(1); uv.row(1).array() += lo1;

    gsMatrix<T> N;
    normalNumerator(hom, sgn, uv, N);          // 3 x nPts
    Ncoef = s_lu.at(key).solve(N.transpose()); // nPts x 3
}

/// Fills every derived field of \a el (uv box, AABB, bidegree, Bernstein
/// certificate) from its already-set \a patch, \a sgn, \a rational, \a hom.
template<class T>
void fillElementFields(BezElement<T> & el)
{
    const gsMatrix<T> supp = el.hom.support();
    for (short_t d = 0; d != 2; ++d) { el.uv[d][0] = supp(d,0); el.uv[d][1] = supp(d,1); }

    const gsMatrix<T> & C = el.hom.coefs();
    for (short_t k = 0; k != 3; ++k)
    {
        el.lo[k] =  std::numeric_limits<T>::max();
        el.up[k] = -std::numeric_limits<T>::max();
    }
    for (index_t r = 0; r != C.rows(); ++r)
    {
        const T w = C(r,3);
        GISMO_ENSURE(w > 0, "gsBRepRayCast: non-positive weight encountered in a Bezier element.");
        for (short_t k = 0; k != 3; ++k)
        {
            const T x = C(r,k)/w;
            el.lo[k] = std::min(el.lo[k], x);
            el.up[k] = std::max(el.up[k], x);
        }
    }

    const index_t p = el.hom.basis().degree(0), q = el.hom.basis().degree(1);
    el.dU = el.rational ? 4*p-1 : 2*p-1;
    el.dV = el.rational ? 4*q-1 : 2*q-1;

    fitNcoef(el.dU, el.dV, el.uv[0][0], el.uv[0][1], el.uv[1][0], el.uv[1][1],
             el.hom, el.sgn, el.Ncoef);

    for (short_t k = 0; k != 3; ++k)
    {
        el.Nmin[k] = el.Ncoef.col(k).minCoeff();
        el.Nmax[k] = el.Ncoef.col(k).maxCoeff();
    }
    el.Nscale = el.Ncoef.cwiseAbs().maxCoeff();

    const T delta = (T)1e-10 * el.Nscale;
    for (short_t k = 0; k != 3; ++k)
    {
        el.certWeak[k]   = (el.Nmin[k] >= -delta && el.Nmax[k] >  delta) ?  1
                         : (el.Nmax[k] <=  delta && el.Nmin[k] < -delta) ? -1 : 0;
        el.certStrict[k] = (el.Nmin[k] >  delta) ?  1
                         : (el.Nmax[k] < -delta) ? -1 : 0;
    }
}

/// True iff every control point of \a C, projected by its 4th (weight) column,
/// equals its first row's projection bitwise; \a c receives that value.
template<class T>
bool constantCoord(const gsMatrix<T> & C, short_t k, T & c)
{
    c = C(0,k)/C(0,3);
    for (index_t r = 1; r < C.rows(); ++r)
        if (C(r,k)/C(r,3) != c)
            return false;
    return true;
}

} // namespace detail

/**
    @brief Bezier elements of a spline BRep, per-patch, from splitAtMult(1)
    on the 4-column projective gsTensorBSpline<2,T> (never gsMultiPatch::extractBezier,
    which rebuilds every element on [0,1] and loses both the patch id and the uv box).

    \a sgn[p] must already carry the outward-orientation multiplier of patch p.
*/
template<class T>
std::vector<BezElement<T> > extractElements(const gsMultiPatch<T> & mp, const std::vector<T> & sgn)
{
    GISMO_ENSURE(sgn.size() == mp.nPatches(), "extractElements: sgn.size() must equal mp.nPatches().");

    std::vector<BezElement<T> > elements;
    for (size_t p = 0; p != mp.nPatches(); ++p)
    {
        bool rational = false;
        const gsTensorBSpline<2,T> hom = detail::homogenize(mp.patch(p), rational);

        std::vector<gsGeometry<T>*> pieces = hom.splitAtMult(1);
        for (gsGeometry<T> * piece : pieces)
        {
            BezElement<T> el;
            el.patch    = (index_t)p;
            el.sgn      = sgn[p];
            el.rational = rational;
            el.hom      = *static_cast<gsTensorBSpline<2,T>*>(piece);
            detail::fillElementFields(el);
            elements.push_back(give(el));
        }
        freeAll(pieces);
    }
    return elements;
}

/**
    @brief Splits \a el into its 4 uniformSplit children (knot-range midpoint
    in each of u,v), recomputing the Bernstein certificate (P2) for each.
    Children keep \a el's \a patch and \a sgn; their uv box is their own
    support(). Used by ::hitElement (Newton stall/non-convergence) and by
    ::castLine (uncertified elements).
*/
template<class T>
std::vector<BezElement<T> > splitElement(const BezElement<T> & el)
{
    std::vector<gsGeometry<T>*> pieces = el.hom.uniformSplit(-1);

    std::vector<BezElement<T> > children;
    children.reserve(pieces.size());
    for (gsGeometry<T> * piece : pieces)
    {
        BezElement<T> child;
        child.patch    = el.patch;
        child.sgn      = el.sgn;
        child.rational = el.rational;
        child.hom      = *static_cast<gsTensorBSpline<2,T>*>(piece);
        detail::fillElementFields(child);
        children.push_back(give(child));
    }
    freeAll(pieces);
    return children;
}

/**
    @brief The normal numerator N := sgn * w^4 * (S_u x S_v), evaluated exactly
    from the homogeneous net without ever dividing by w.

    Writing the homogeneous coordinates as X = (w*S_x, w*S_y, w*S_z) and the
    weight w (the hom patch's 4th component), the affine derivatives are
    S_u = (X_u*w - X*w_u)/w^2 =: A/w^2 and S_v =: B/w^2, so
    N = sgn*(A x B) = sgn*w^4*(S_u x S_v); for w === 1 this reduces to
    sgn*(S_u x S_v). \a N is 3 x nPts, one column per column of \a uv.

    Bidegree of N (used to size ::BezElement::Ncoef), with p,q the element's
    degrees in u,v:
    - polynomial elements (w === 1, so w_u = w_v = 0): A = X_u, B = X_v have
      bidegree (p-1,q) and (p,q-1) (one derivative each removes one degree
      in that direction); a component of A x B is a difference of products
      of one A- and one B-component, bidegree (p-1)+p, q+(q-1) =
      (2p-1, 2q-1).
    - rational elements: w itself has bidegree (p,q), so X_u*w and X*w_u
      (bidegree (p-1,q)+(p,q) and (p,q)+(p-1,q)) both give A bidegree
      (2p-1, 2q); symmetrically B = X_v*w - X*w_v has bidegree (2p, 2q-1).
      A component of A x B is then (2p-1)+(2p), (2q)+(2q-1) = (4p-1, 4q-1).

    Certificate (in ::BezElement, with the relative margin delta = 1e-10*Nscale,
    Nscale = max|Bernstein coef| over all 3 components):
    - certWeak[k] = +1 iff Nmin[k] >= -delta && Nmax[k] > delta (all Bernstein
      coefs of N_k are non-negative up to the margin, and not uniformly ~0);
      -1 symmetrically; else 0.
    - certStrict[k] = +1 iff Nmin[k] > delta; -1 iff Nmax[k] < -delta; else 0.

    The property relied on: a Bernstein polynomial with all coefficients >= 0,
    not all 0, is > 0 throughout the OPEN element (convex-hull property), so
    certStrict[k] != 0 makes the e_k-projection (u,v) -> (S_i,S_j) a local
    diffeomorphism there. certWeak is
    the same test relaxed by the margin, needed because a polynomial that is
    nonnegative in the open element can still have individual Bernstein
    coefficients slightly negative (the basis is not the monomial basis): the
    sphere's 8 elements are never certStrict in any k (N === 0 identically on
    the pole/meridian edges forces at least one Bernstein coefficient to be
    exactly 0), so ::hitElement is only ever dispatched there under certWeak.
    ::windingOK is the safety net for the injectivity assumption this
    certificate makes on small Bezier elements.
*/
template<class T>
void normalNumerator(const gsTensorBSpline<2,T> & hom, T sgn, const gsMatrix<T> & uv, gsMatrix<T> & N)
{
    gsMatrix<T> val, der;
    hom.eval_into (uv, val); // 4 x n
    hom.deriv_into(uv, der); // 8 x n, row 2*k+j = d X_k / d u_j

    const index_t n = uv.cols();
    N.resize(3, n);
    for (index_t c = 0; c != n; ++c)
    {
        const T w  = val(3,c);
        const T wu = der(2*3+0,c);
        const T wv = der(2*3+1,c);

        gsVector<T,3> A, B;
        for (short_t k = 0; k != 3; ++k)
        {
            const T Xk  = val(k,c);
            const T Xku = der(2*k+0,c);
            const T Xkv = der(2*k+1,c);
            A[k] = Xku*w - Xk*wu;
            B[k] = Xkv*w - Xk*wv;
        }
        N.col(c) = sgn * A.cross(B);
    }
}

namespace detail {

/// Safeguarded-Newton solve of S_i(u,v) = y_i, S_j(u,v) = y_j inside \a el's
/// uv box, starting at the box center and clamping every iterate to the box.
///
/// The STOPPING rule (step < 1e-15*span, a singular Jacobian, or 50
/// iterations) is deliberately independent of the residual: Newton's
/// quadratic convergence means the step size, not the residual, is the
/// tight indicator of "no further progress possible in double precision",
/// and stopping as soon as the residual first drops below the 1e-13*L
/// ACCEPT threshold (checked only once, below, after stopping) would give
/// up one or more Newton steps early -- exactly the steps that would have
/// squared an O(1e-13*L) residual down to O(1e-26*L), i.e. machine
/// precision. Every derived quantity computed from (u,v) afterwards (in
/// particular the hit parameter t = S_k(u,v), a value Newton never
/// constrains directly) inherits (u,v)'s own error, so under-converging
/// here would silently violate the caller's stated 1e-13*L accuracy on t.
///
/// Returns 1 (accepted: final residual <= 1e-13*L), 0 (stalled: stopped
/// with the residual still above tolerance, either via a tiny step or a
/// singular Jacobian, or clamped at the iteration budget) or -1 (ran out of
/// 50 iterations while still taking unclamped, non-tiny steps: genuine
/// non-convergence, not a stall).
template<class T>
int newtonSolve(const BezElement<T> & el, short_t i, short_t j, T yi, T yj, T L, T & u, T & v)
{
    const T lo0 = el.uv[0][0], hi0 = el.uv[0][1];
    const T lo1 = el.uv[1][0], hi1 = el.uv[1][1];
    const T maxSpan = std::max(hi0-lo0, hi1-lo1);

    u = (T)0.5*(lo0+hi0);
    v = (T)0.5*(lo1+hi1);
    bool lastClamped = false, tinyStep = false, singular = false;

    int it = 0;
    for (; it != 50; ++it)
    {
        gsMatrix<T> pt(2,1); pt(0,0) = u; pt(1,0) = v;
        gsMatrix<T> val, der;
        el.hom.eval_into (pt, val);
        el.hom.deriv_into(pt, der);

        const T w = val(3,0), wu = der(6,0), wv = der(7,0);
        const T Gi = val(i,0) - yi*w, Gj = val(j,0) - yj*w;

        const T J00 = der(2*i,0)-yi*wu, J01 = der(2*i+1,0)-yi*wv;
        const T J10 = der(2*j,0)-yj*wu, J11 = der(2*j+1,0)-yj*wv;
        const T det = J00*J11 - J01*J10;
        if (det == 0) { singular = true; break; }

        const T du = ( J11*(-Gi) - J01*(-Gj))/det;
        const T dv = (-J10*(-Gi) + J00*(-Gj))/det;

        T un = u+du, vn = v+dv;
        lastClamped = false;
        if (un < lo0) { un = lo0; lastClamped = true; } if (un > hi0) { un = hi0; lastClamped = true; }
        if (vn < lo1) { vn = lo1; lastClamped = true; } if (vn > hi1) { vn = hi1; lastClamped = true; }

        const T step = std::max(std::abs(un-u), std::abs(vn-v));
        u = un; v = vn;

        if (step < (T)1e-15*maxSpan) { tinyStep = true; break; }
    }
    const bool ranOutOfIters = (it == 50);

    gsMatrix<T> pt(2,1); pt(0,0) = u; pt(1,0) = v;
    gsMatrix<T> val;
    el.hom.eval_into(pt, val);
    const T w = val(3,0);
    // Accept on the projected residual max(|S_i-y_i|,|S_j-y_j|); the
    // homogeneous residual G = X - y*w equals w times it and, where w < 1,
    // would accept a larger projected error.
    const T res = std::max(std::abs(val(i,0)/w-yi), std::abs(val(j,0)/w-yj));
    if (res <= (T)1e-13*L)
        return 1;

    const bool stalled = tinyStep || singular || (ranOutOfIters && lastClamped);
    return stalled ? 0 : -1;
}

/// True iff (\a qi,\a qj) lies in (or on) the convex hull of the 2D points
/// \a pts (2 x n), via a monotone-chain hull followed by a half-plane test
/// against every hull edge. Ties favour "inside" (a point exactly on an
/// edge counts as inside), which is the safe direction for the caller
/// below: it would rather over-include a target than silently exclude one
/// that a genuine root could still reach.
template<class T>
bool pointInConvexHull(const gsMatrix<T> & pts, T qi, T qj)
{
    const index_t n = pts.cols();
    if (n == 0) return false;
    if (n == 1) return qi == pts(0,0) && qj == pts(1,0);

    std::vector<index_t> idx(n);
    for (index_t t = 0; t != n; ++t) idx[t] = t;
    std::sort(idx.begin(), idx.end(), [&pts](index_t a, index_t b)
    {
        if (pts(0,a) != pts(0,b)) return pts(0,a) < pts(0,b);
        return pts(1,a) < pts(1,b);
    });

    // Signed area x2 of (o,a,b); > 0 iff o->a->b turns left (CCW).
    const auto cross = [&pts](index_t o, index_t a, index_t b) -> T
    {
        return (pts(0,a)-pts(0,o))*(pts(1,b)-pts(1,o))
             - (pts(1,a)-pts(1,o))*(pts(0,b)-pts(0,o));
    };

    std::vector<index_t> hull;
    for (index_t t = 0; t != n; ++t)
    {
        const index_t p = idx[t];
        while (hull.size() >= 2 && cross(hull[hull.size()-2], hull[hull.size()-1], p) <= 0)
            hull.pop_back();
        hull.push_back(p);
    }
    const size_t lowerSize = hull.size()+1;
    for (index_t t = n-2; t >= 0; --t)
    {
        const index_t p = idx[t];
        while (hull.size() >= lowerSize && cross(hull[hull.size()-2], hull[hull.size()-1], p) <= 0)
            hull.pop_back();
        hull.push_back(p);
    }
    hull.pop_back(); // duplicate of the first point

    if (hull.size() < 3)
    {
        T lo0 = pts(0,hull[0]), hi0 = lo0, lo1 = pts(1,hull[0]), hi1 = lo1;
        for (size_t t = 0; t != hull.size(); ++t)
        {
            lo0 = std::min(lo0, pts(0,hull[t])); hi0 = std::max(hi0, pts(0,hull[t]));
            lo1 = std::min(lo1, pts(1,hull[t])); hi1 = std::max(hi1, pts(1,hull[t]));
        }
        return qi >= lo0 && qi <= hi0 && qj >= lo1 && qj <= hi1;
    }

    const size_t m = hull.size();
    for (size_t t = 0; t != m; ++t)
    {
        const index_t a = hull[t], b = hull[(t+1)%m];
        const T cr = (pts(0,b)-pts(0,a))*(qj-pts(1,a)) - (pts(1,b)-pts(1,a))*(qi-pts(0,a));
        if (cr < 0) return false;
    }
    return true;
}

/// Recursive body of ::hitElement: on a Newton stall/non-convergence,
/// subdivides via ::splitElement and recurses (depth <= maxSplit,
/// certification inherited); see ::hitElement for the full return contract,
/// including the certStrict-scoped rule that decides 0 vs. -1 at
/// \a maxSplit. A -1 from any child propagates immediately;
/// otherwise only the first child 1 found is kept, which drops a second,
/// distinct root inside the same certified element (a non-injective
/// e_k-projection) rather than reporting it -- ::windingOK is the safety
/// net for exactly that assumption, not this function.
template<class T>
int hitElementImpl(const BezElement<T> & el, short_t k, T yi, T yj, Hit<T> & hit, T L, int depth, int maxSplit)
{
    const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);
    if (!(el.lo[i] <= yi && yi < el.up[i] && el.lo[j] <= yj && yj < el.up[j]))
        return 0;

    GISMO_ASSERT(el.certWeak[k] != 0, "hitElement: called on an element uncertified in direction k.");

    T u, v;
    const int code = newtonSolve(el, i, j, yi, yj, L, u, v);

    if (code == 1)
    {
        const T tol = (T)1e-12*std::max(el.uv[0][1]-el.uv[0][0], el.uv[1][1]-el.uv[1][0]);
        u = std::min(std::max(u, el.uv[0][0]), el.uv[0][1]);
        v = std::min(std::max(v, el.uv[1][0]), el.uv[1][1]);
        const bool onBoundary = (u-el.uv[0][0] <= tol) || (el.uv[0][1]-u <= tol)
                              || (v-el.uv[1][0] <= tol) || (el.uv[1][1]-v <= tol);

        gsMatrix<T> uv(2,1); uv(0,0) = u; uv(1,0) = v;
        gsMatrix<T> N;
        normalNumerator(el.hom, el.sgn, uv, N);
        const T nrm = N.col(0).norm();
        if (nrm == 0)
            return -1; // pole: N vanishes at the footpoint

        hit.patch      = el.patch;
        hit.u          = u;
        hit.v          = v;
        hit.n          = N.col(0)/nrm;
        hit.onBoundary = onBoundary;
        hit.elem       = -1;

        T c;
        if (constantCoord(el.hom.coefs(), k, c))
            hit.t = c;
        else
        {
            gsMatrix<T> val;
            el.hom.eval_into(uv, val);
            hit.t = val(k,0)/val(3,0);
        }
        return 1;
    }

    if (depth >= maxSplit)
    {
        // On a leaf strictly certified in direction k, N_k is bounded away
        // from 0 (BezElement's certStrict), so the transverse Jacobian is
        // nonsingular everywhere on the leaf and the e_k-projection is a
        // local diffeomorphism there. A stall or non-convergence on such a
        // leaf is ASSUMED to mean there is no root (local invertibility does
        // not guarantee Newton finds one that exists); newtonSolve's label is
        // used as-is and windingOK is the safety net for a missed hit.
        //
        // On a leaf that is only weakly certified, N_k can
        // vanish somewhere on the leaf -- in particular along a coordinate
        // singularity such as a sphere pole, where the AABB test (necessary,
        // never sufficient) never shrinks under repeated subdivision, so
        // Newton can be asked to solve a box that provably has no root and
        // one that has a root but is numerically hard to converge in look
        // identical to newtonSolve. The convex hull of the leaf's own
        // (projected) control net tells them apart, sound up to
        // floating-point rounding in its half-plane test: a rational Bezier
        // patch's image lies inside the hull of its own control points, so
        // a target outside it has no root here; a target inside it might
        // still have one, so the failure must propagate (-1) rather than be
        // silently read as "no hit" (0).
        //
        // The hull test is coarser than a nonsingular Jacobian, though: near
        // an ORDINARY curved element edge (not a singularity) it can still
        // include a target that the neighbouring element already resolved
        // correctly, which is why it is confined to leaves that are not
        // already known-nonsingular via certStrict.
        if (el.certStrict[k] != 0)
            return code;

        gsMatrix<T> proj(2, el.hom.coefs().rows());
        for (index_t r = 0; r != el.hom.coefs().rows(); ++r)
        {
            proj(0,r) = el.hom.coefs()(r,i)/el.hom.coefs()(r,3);
            proj(1,r) = el.hom.coefs()(r,j)/el.hom.coefs()(r,3);
        }
        return pointInConvexHull(proj, yi, yj) ? -1 : 0;
    }

    const std::vector<BezElement<T> > kids = splitElement(el);
    bool found = false;
    Hit<T> best;
    for (typename std::vector<BezElement<T> >::const_iterator c = kids.begin(); c != kids.end(); ++c)
    {
        Hit<T> h;
        const int r = hitElementImpl(*c, k, yi, yj, h, L, depth+1, maxSplit);
        if (r < 0) return -1;
        if (r == 1 && !found) { best = h; found = true; }
    }
    if (found) { hit = best; return 1; }
    return 0;
}

} // namespace detail

/**
    @brief Ray/element intersection along e_k, i.e. the exact root of
    S_i(u,v) = y_i, S_j(u,v) = y_j inside \a el (i,j the transverse indices
    of k, see the file doxygen). Only meaningful on an element already
    certified in direction k (::BezElement::certWeak[k] != 0), whose AABB
    meets the line -- ::castLine enforces both.

    Returns 1 on a hit: a converged (u,v) with a nonzero normal. Every other
    return is 0 or -1, from one of these sources (at any recursion depth up
    to \a maxSplit, not only at the leaf):
    - 0 if the element's AABB does not meet the line at all;
    - -1 if a converged (u,v) has a zero normal (a pole);
    - short of a hit, ::detail::newtonSolve reports a stall (tiny step,
      singular Jacobian, or ran out of iterations while clamped) or a
      non-convergence (ran out of iterations while still taking large,
      non-tiny, unclamped steps). Below \a maxSplit this subdivides
      (::splitElement) and recurses regardless of which of the two it was;
      the four children's results are combined by propagating any -1
      immediately and otherwise keeping the first 1 found (0 if none).
    - AT \a maxSplit itself, the rule depends on the leaf's own
      ::BezElement::certStrict[k]:
      - strictly certified (certStrict[k] != 0, so N_k is bounded away from
        0 and the e_k-projection is a local diffeomorphism on the whole
        leaf): newtonSolve's own stall/non-convergence label is taken as-is
        (stall -> 0, non-convergence -> -1). This ASSUMES that a failure to
        converge on such a leaf means there is no root on it; local
        invertibility does not guarantee that Newton finds an existing root,
        so a missed hit is possible in principle and ::windingOK is the
        safety net that detects it.
      - not strictly certified (in particular along a coordinate
        singularity such as a sphere pole, where N_k can vanish and the
        AABB test -- necessary, never sufficient -- never shrinks under
        subdivision): the label is replaced by a test against the convex
        hull of the leaf's own projected control net (0 iff (\a yi,\a yj)
        lies outside it, else -1), sound up to floating-point rounding in
        its half-plane test. This is coarser than the certStrict rule
        (near an ordinary curved edge the hull of a non-strict leaf can
        still contain a target a neighbouring element already resolved),
        which is exactly why it is confined to leaves that are not already
        known-nonsingular.

    \a L is the length scale (max edge of the background box), used to make
    the residual tolerance 1e-13*L scale-invariant.
*/
template<class T>
int hitElement(const BezElement<T> & el, short_t k, T yi, T yj, Hit<T> & hit, T L, int maxSplit = 8)
{
    return detail::hitElementImpl(el, k, yi, yj, hit, L, 0, maxSplit);
}

/// Winding-number check along e_k: hits are grouped by exactly equal t (so a
/// merged multi-element crossing counts once), each group adds its signed
/// entries/exits (n_k < 0 enters, n_k > 0 exits) at once, and after every
/// group the running count must lie in {0,1} and must end at 0. Any n_k == 0
/// fails immediately (the ray is then tangent to the surface, which the
/// certificate in ::normalNumerator should have already excluded upstream).
template<class T>
bool windingOK(const std::vector<Hit<T> > & hits, short_t k)
{
    int count = 0;
    size_t idx = 0;
    while (idx != hits.size())
    {
        const T t0 = hits[idx].t;
        int delta = 0;
        size_t j = idx;
        for (; j != hits.size() && hits[j].t == t0; ++j)
        {
            if (hits[j].n[k] == 0) return false;
            delta += (hits[j].n[k] < 0) ? 1 : -1;
        }
        count += delta;
        if (count != 0 && count != 1) return false;
        idx = j;
    }
    return count == 0;
}

/// Running winding count along e_k just after every hit with t <= \a a, i.e.
/// after every group of bitwise-equal t up to and including \a a: each hit
/// contributes +1 (n[k] < 0, an entry) or -1 (n[k] > 0, an exit), applied in
/// t order (\a hits must already be sorted by t, as ::castLine produces).
/// Precondition: ::windingOK(hits,k). Under that precondition this equals
/// ::clipHits' windingAtA (nBefore % 2) numerically -- every hit contributes
/// an odd number, so a sum of nBefore of them is congruent to nBefore mod 2,
/// and ::windingOK's own invariant already keeps that sum in {0,1} at every
/// group boundary. The two are computed separately because ::insideIntervals
/// needs this same signed, per-hit computation again, one group at a time,
/// to detect the walk's later 0->1/1->0 transitions -- a single parity at
/// one point (::clipHits' use) does not give that.
template<class T>
int windingAt(const std::vector<Hit<T> > & hits, short_t k, T a)
{
    int count = 0;
    for (typename std::vector<Hit<T> >::const_iterator h = hits.begin(); h != hits.end() && h->t <= a; ++h)
        count += (h->n[k] < 0) ? 1 : -1;
    GISMO_ENSURE(count == 0 || count == 1,
                 "windingAt: running winding count " << count << " outside {0,1}; "
                 "hits must satisfy windingOK.");
    return count;
}

/// The inside sub-intervals of (\a a, \a b] along e_k, i.e. the t-ranges
/// where the running winding count (::windingAt) equals 1, clipped to that
/// half-open range. Hits are walked in groups of bitwise-equal t (the same
/// grouping as ::windingOK), each group's signed entries/exits applied at
/// once: a 0->1 transition opens an interval at that t, a 1->0 transition
/// closes it there. A zero-length interval (the open and close t bitwise
/// equal) is skipped, never tolerance-merged with its neighbour.
template<class T>
void insideIntervals(const std::vector<Hit<T> > & hits, short_t k, T a, T b,
                     std::vector<std::pair<T,T> > & iv)
{
    iv.clear();
    int cnt = windingAt(hits, k, a);
    T start = a;

    size_t idx = 0;
    while (idx != hits.size() && hits[idx].t <= a) ++idx;

    while (idx != hits.size() && hits[idx].t <= b)
    {
        const T t0 = hits[idx].t;
        int delta = 0;
        size_t j = idx;
        for (; j != hits.size() && hits[j].t == t0; ++j)
            delta += (hits[j].n[k] < 0) ? 1 : -1;

        const int prev = cnt;
        cnt += delta;
        if (prev == 0 && cnt == 1) start = t0;
        else if (prev == 1 && cnt == 0 && t0 > start) iv.push_back(std::make_pair(start, t0));
        idx = j;
    }

    if (cnt == 1 && b > start)
        iv.push_back(std::make_pair(start, b));
}

namespace detail {

/// Recursive body of ::castLine for a single top-level element: subdivides
/// uncertified elements (::splitElement, depth <= maxSplit) and otherwise
/// dispatches to ::hitElement, appending any hit to \a hits. Returns false
/// (without aborting the sweep over sibling/child elements) on an
/// uncertified piece surviving to maxSplit or a failing ::hitElement call.
template<class T>
bool castLineOneElement(const BezElement<T> & el, short_t k, T yi, T yj,
                        std::vector<Hit<T> > & hits, T L, int depth, int maxSplit)
{
    const short_t i = (short_t)((k+1)%3), j = (short_t)((k+2)%3);
    if (!(el.lo[i] <= yi && yi < el.up[i] && el.lo[j] <= yj && yj < el.up[j]))
        return true;

    if (el.certWeak[k] == 0)
    {
        if (depth >= maxSplit)
            return false;
        const std::vector<BezElement<T> > kids = splitElement(el);
        bool ok = true;
        for (typename std::vector<BezElement<T> >::const_iterator c = kids.begin(); c != kids.end(); ++c)
            if (!castLineOneElement(*c, k, yi, yj, hits, L, depth+1, maxSplit))
                ok = false;
        return ok;
    }

    Hit<T> hit;
    const int r = hitElement(el, k, yi, yj, hit, L);
    if (r < 0) return false;
    if (r == 1) hits.push_back(hit);
    return true;
}

} // namespace detail

/**
    @brief Full ray cast along e_k through every element of \a els whose AABB
    meets the line at transverse coordinates (\a yi,\a yj): elements
    uncertified in direction k are subdivided (::splitElement) and recursed
    into (depth <= \a maxSplit, inherited certification, since subdivision
    coefficients are convex combinations); certified elements go through
    ::hitElement.

    Every appended ::Hit records, in \a elem, the index into \a els of the
    top-level element ::castLineOneElement was called on (never one of its
    ::splitElement children) -- this is the index a caller must use to look
    a hit back up in the same \a els it passed in, e.g. to test whether a hit
    lies within a clip range's own element list.

    \a hits is sorted by t, then same-patch pairs within 1e-12*spanU(patch)
    in u and 1e-12*spanV(patch) in v are merged (checking every pair of the
    same patch, never across patches) -- this removes the duplicate hit two
    neighbouring elements report on a shared edge.

    Returns false if any uncertified element remains at \a maxSplit, if any
    ::hitElement call fails (returns -1, including a zero normal), or if the
    final deduplicated hit list fails ::windingOK; otherwise returns
    ::windingOK's result.
*/
template<class T>
bool castLine(const std::vector<BezElement<T> > & els, short_t k, T yi, T yj,
              std::vector<Hit<T> > & hits, T L, int maxSplit = 16)
{
    hits.clear();
    bool ok = true;
    for (size_t e = 0; e != els.size(); ++e)
    {
        const size_t before = hits.size();
        if (!detail::castLineOneElement(els[e], k, yi, yj, hits, L, 0, maxSplit))
            ok = false;
        for (size_t q = before; q != hits.size(); ++q)
            hits[q].elem = (index_t)e;
    }
    if (!ok) return false;

    std::sort(hits.begin(), hits.end(),
              [](const Hit<T> & a, const Hit<T> & b) { return a.t < b.t; });

    std::map<index_t, std::array<T,4> > span; // patch -> (minU,maxU,minV,maxV)
    for (typename std::vector<BezElement<T> >::const_iterator el = els.begin(); el != els.end(); ++el)
    {
        typename std::map<index_t, std::array<T,4> >::iterator it = span.find(el->patch);
        if (it == span.end())
            span[el->patch] = {el->uv[0][0], el->uv[0][1], el->uv[1][0], el->uv[1][1]};
        else
        {
            it->second[0] = std::min(it->second[0], el->uv[0][0]);
            it->second[1] = std::max(it->second[1], el->uv[0][1]);
            it->second[2] = std::min(it->second[2], el->uv[1][0]);
            it->second[3] = std::max(it->second[3], el->uv[1][1]);
        }
    }

    std::vector<bool> removed(hits.size(), false);
    for (size_t a = 0; a < hits.size(); ++a)
    {
        if (removed[a]) continue;
        for (size_t b = a+1; b < hits.size(); ++b)
        {
            if (removed[b] || hits[a].patch != hits[b].patch) continue;
            const std::array<T,4> & sp = span.at(hits[a].patch);
            const T spanU = sp[1]-sp[0], spanV = sp[3]-sp[2];
            if (std::abs(hits[a].u-hits[b].u) <= (T)1e-12*spanU &&
                std::abs(hits[a].v-hits[b].v) <= (T)1e-12*spanV)
                removed[b] = true;
        }
    }
    std::vector<Hit<T> > merged;
    merged.reserve(hits.size());
    for (size_t a = 0; a < hits.size(); ++a)
        if (!removed[a]) merged.push_back(hits[a]);
    hits.swap(merged);

    return windingOK(hits, k);
}

/**
    @brief Splits \a hits (already sorted by t, as produced by ::castLine)
    into the sub-range with a < t <= b (the "(a,b]" clip convention) and the
    winding count \a windingAtA just before \a a.

    Every ::Hit contributes exactly +-1 to the running winding count
    (::windingOK), so once that count is known to lie in {0,1} after each
    group of equal t (which ::windingOK verifies), it equals, at any point
    along the line, the PARITY of the number of hits seen so far -- \a
    windingAtA is computed that way, without needing the ray direction k or
    the per-hit normal.
*/
template<class T>
void clipHits(const std::vector<Hit<T> > & hits, T a, T b, std::vector<Hit<T> > & in, int & windingAtA)
{
    in.clear();
    index_t nBefore = 0;
    for (typename std::vector<Hit<T> >::const_iterator h = hits.begin(); h != hits.end(); ++h)
    {
        if (h->t <= a) ++nBefore;
        if (a < h->t && h->t <= b) in.push_back(*h);
    }
    windingAtA = nBefore % 2;
}

/**
    @brief The boundary iso-curves of every element in \a els, on the same
    4-column homogeneous net as the element (via gsTensorBSpline::slice, so
    endpoints are exact, not evaluated).

    Every element contributes its lower-u (dirFixed=0, par=uv[0][0]) and
    lower-v (dirFixed=1, par=uv[1][0]) edges. An upper edge (par=uv[d][1]) is
    added only where that value equals the patch's domain end in direction d
    (the maximum uv[d][1] over all elements of that patch), so each edge of
    the patch's own boundary/interior-knot-line structure appears exactly
    once from its owning element, by exact knot equality -- not once per
    adjacent element pair. There is deliberately no cross-patch
    deduplication: two BRep patches sharing a physical edge each contribute
    their own copy, and any resulting duplicate crossing only ever produces a
    zero-length interval downstream, never a spurious one.

    Degenerate curves (every projected control point bitwise equal, as on
    the sphere's pole/meridian edges) are dropped.
*/
template<class T>
std::vector<EdgeCurve<T> > elementEdges(const std::vector<BezElement<T> > & els)
{
    std::map<index_t, T> hiU, hiV;
    for (typename std::vector<BezElement<T> >::const_iterator el = els.begin(); el != els.end(); ++el)
    {
        typename std::map<index_t,T>::iterator itU = hiU.find(el->patch);
        if (itU == hiU.end()) hiU[el->patch] = el->uv[0][1];
        else itU->second = std::max(itU->second, el->uv[0][1]);
        typename std::map<index_t,T>::iterator itV = hiV.find(el->patch);
        if (itV == hiV.end()) hiV[el->patch] = el->uv[1][1];
        else itV->second = std::max(itV->second, el->uv[1][1]);
    }

    std::vector<EdgeCurve<T> > edges;
    for (typename std::vector<BezElement<T> >::const_iterator el = els.begin(); el != els.end(); ++el)
    {
        std::vector<std::pair<short_t,T> > cands;
        cands.push_back(std::make_pair((short_t)0, el->uv[0][0]));
        cands.push_back(std::make_pair((short_t)1, el->uv[1][0]));
        if (el->uv[0][1] == hiU[el->patch]) cands.push_back(std::make_pair((short_t)0, el->uv[0][1]));
        if (el->uv[1][1] == hiV[el->patch]) cands.push_back(std::make_pair((short_t)1, el->uv[1][1]));

        for (size_t c = 0; c != cands.size(); ++c)
        {
            gsBSpline<T> curve;
            el->hom.slice(cands[c].first, cands[c].second, curve);

            const gsMatrix<T> & CC = curve.coefs();
            T c0;
            bool degenerate = detail::constantCoord(CC, 0, c0);
            if (degenerate)
            {
                T ctmp;
                degenerate = detail::constantCoord(CC, 1, ctmp) && detail::constantCoord(CC, 2, ctmp);
            }
            if (degenerate) continue;

            EdgeCurve<T> ec;
            ec.patch    = el->patch;
            ec.dirFixed = cands[c].first;
            ec.par      = cands[c].second;
            ec.hom      = curve;
            edges.push_back(give(ec));
        }
    }
    return edges;
}

namespace detail {

/// Recursive body of ::curvePlaneCrossings on one curve segment: the
/// exact-constant rule, the sign-change count that isolates 0/1 crossing
/// per piece, the >=2-change subdivision (::gsBSpline::splitAt, depth <=
/// 40) and the 1-change bisection are all described at ::curvePlaneCrossings.
template<class T>
void crossingsOnCurve(const gsBSpline<T> & curve, short_t k, T X, int depth, std::vector<PlaneCrossing<T> > & out)
{
    const gsMatrix<T> & C = curve.coefs();
    const index_t n = C.rows();

    // Exact-constant rule, applied first: if x_k is bitwise the same value c
    // on every control point, the curve can only ever be entirely on one
    // side of the plane, by construction -- never a crossing, regardless of
    // c. This is checked on the PROJECTED coordinate (not on g_m below),
    // because for a rational curve g_m = X_{k,m} - X*w_m can round
    // differently across m even when X_{k,m}/w_m is exactly the same c for
    // every m: c was obtained by a division, and c*w_m need not round-trip
    // back to the stored X_{k,m} bit-for-bit when w_m varies across m, which
    // can flip an individual g_m's sign spuriously.
    T c;
    if (constantCoord(C, k, c)) return;

    std::vector<bool> above(n);
    for (index_t m = 0; m != n; ++m)
        above[m] = (C(m,k) - X*C(m,3)) > 0;

    index_t changes = 0;
    for (index_t m = 1; m != n; ++m)
        if (above[m] != above[m-1]) ++changes;

    if (changes == 0) return;

    if (changes >= 2)
    {
        GISMO_ENSURE(depth < 40, "curvePlaneCrossings: exceeded the recursion depth of 40.");
        const T mid = (T)0.5*(curve.domainStart() + curve.domainEnd());
        gsBSpline<T> left, right;
        curve.splitAt(mid, left, right);
        crossingsOnCurve(left,  k, X, depth+1, out);
        crossingsOnCurve(right, k, X, depth+1, out);
        return;
    }

    // Exactly one strict sign change over the control-point sequence: by
    // Descartes' rule for Bernstein polynomials, the curve's g(s) = x_k(s)-X
    // has at most one root in the open parameter interval, counted with
    // multiplicity, and end predicates differ (the change count is odd).
    T lo = curve.domainStart(), hi = curve.domainEnd();
    bool bLo = above[0], bHi = above[n-1];
    gsVector<T,3> xLo, xHi;
    for (short_t kk = 0; kk != 3; ++kk)
    {
        xLo[kk] = C(0,kk)/C(0,3);
        xHi[kk] = C(n-1,kk)/C(n-1,3);
    }

    for (int it = 0; it != 200; ++it)
    {
        const T mid = (T)0.5*(lo+hi);
        if (mid == lo || mid == hi) break;

        gsMatrix<T> pt(1,1); pt(0,0) = mid;
        gsMatrix<T> val;
        curve.eval_into(pt, val);
        const T w = val(3,0);
        const bool bMid = val(k,0)/w > X; // predicate "above" on the projected point
        gsVector<T,3> xMid;
        for (short_t kk = 0; kk != 3; ++kk) xMid[kk] = val(kk,0)/w;

        if (bMid == bLo) { lo = mid; bLo = bMid; xLo = xMid; }
        else              { hi = mid; bHi = bMid; xHi = xMid; }
    }

    PlaneCrossing<T> pc;
    pc.s   = hi;
    pc.x   = xHi;
    pc.dir = bHi ? 1 : -1;
    out.push_back(pc);
}

} // namespace detail

/**
    @brief Crossings of \a c with the axis-aligned plane x_k = \a X, found by
    predicate-flip bisection on the "above" predicate (x_k > X) of projected
    points, after isolating them via the strict sign changes of
    g_m = X_{k,m} - X*w_m over the curve's homogeneous control points (see
    ::detail::crossingsOnCurve for the exact-Descartes argument). Zero
    coefficients count as "not above" -- matching x_k == X not being
    above -- so a curve that TOUCHES the plane at a control point without
    crossing it still registers the flip there (e.g. a vertical edge from
    z=0.25 to z=0.75 against the plane z=0.25 has g=(0,+), one change, and
    R3's half-open (a,b] clip convention needs that crossing at its lower
    end). The exact-constant rule is checked explicitly, before the g_m sign
    count: if x_k is bitwise the same value on every control point, there is
    no crossing, full stop -- this is NOT the same as "every g_m has the same
    sign", since for a rational curve g_m = X_{k,m} - X*w_m can round
    differently across m (a division-then-multiplication round-trip through
    w_m) even when the projected coordinate is bitwise constant, which could
    otherwise manufacture a spurious change.

    Each isolated piece with more than one change is bisected in the
    parameter (gsBSpline::splitAt, which never rescales), inheriting the
    curve's original parameter values, down to a recursion depth of 40.
*/
template<class T>
void curvePlaneCrossings(const EdgeCurve<T> & c, short_t k, T X, std::vector<PlaneCrossing<T> > & out)
{
    out.clear();
    detail::crossingsOnCurve(c.hom, k, X, 0, out);
}

} // namespace brc
} // namespace gismo
