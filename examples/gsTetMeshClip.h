/** @file gsTetMeshClip.h

    @brief Exact cut-cell quadrature of a tetrahedral-mesh domain against a
    uniform Cartesian background grid, by tet-cell clipping; ASCII gmsh MSH
    4.1 reader; unclipped-mesh oracles. This is the direct 3D analogue of the
    mesh-clipping path of immersed_gauss_green_quadrature_example.cpp
    (optional/gsOpenCascade).

    Theory (Powell & Abel, "An exact general remeshing scheme applied to
    physically conservative voxelization", J. Comput. Phys. 297 (2015)
    340-356): a tetrahedral cell is intersected with a Cartesian background
    cell by six successive half-space splits, x >= xL, x <= xR, y >= yB,
    y <= yT, z >= zN, z <= zF (closed, non-strict half-spaces). Every
    surviving convex piece of a tet-halfspace intersection is itself a tet
    or a triangular prism (a prism splits into 3 tets); after six splits the
    original tet has become a set of sub-tets that exactly partitions
    tet ^ cell, so each sub-tet integrated with an EXACT tet rule sums to
    the exact cell volume contribution -- no approximation beyond floating
    point rounding.

    Predicate and interpolation conventions (exact, no tolerance):
    - a vertex is classified "in" a lower half-space x >= c (resp. upper,
      x <= c) by the closed, non-strict comparison; a vertex exactly ON the
      plane counts as "in" for both a lower and an upper split at the same
      c, so a mixed edge always has exactly one "in" and one "out" endpoint;
    - the intersection of a mixed edge (A,B) with the plane coordinate[axis]
      = c sets that coordinate EXACTLY to c and interpolates the other two
      from the endpoint with the smaller coordinate along axis (so the
      result does not depend on the edge's direction, exactly as
      intersectX/intersectY in the 2D driver); additionally, if either
      endpoint already lies exactly on the plane, that endpoint is returned
      verbatim rather than recomputed through the affine formula, which is
      not guaranteed to reproduce it bit-for-bit. On a mixed edge only the
      "in" endpoint can carry that value (the "out" endpoint fails the
      non-strict test, hence lies strictly on the other side), so this is
      unambiguous;
    - a piece with an exactly zero signed volume is dropped by an exact
      `0.0 == D` test, never a tolerance. A duplicated vertex at P0, or
      P2 == P3, gives D == 0 exactly (D is then a triple product with a
      literally repeated column). A duplicated vertex at P1 -- P1 == P2 or
      P1 == P3 -- computes D as e1.(e1 x e3) for some edge e3, which is
      mathematically zero but rounds to O(eps) rather than landing on the
      exact bit pattern 0.0; such a piece is kept, with a near-zero |D|
      contributing a near-zero weight, and its rounding-level residual is
      absorbed by T1's tolAbs term rather than by this test.

    Vertex classification and split cases, tet (a,b,c,d), I(p,q) the
    intersection of edge (p,q) with the current plane:
    - 4 in: keep the tet unchanged;
    - 0 in: drop;
    - 1 in (a in): one tet (a, I(a,b), I(a,c), I(a,d));
    - 3 in (a,b,c in; d out): prism (a,b,c | I(a,d),I(b,d),I(c,d));
    - 2 in (a,b in; c,d out): prism (a,I(a,c),I(a,d) | b,I(b,c),I(b,d));
    - a prism (A,B,C | A',B',C') is always convex here (a tet cut by one
      half-space), so it splits into 3 tets without further case analysis:
      (A,B,C,A'), (B,C,A',B'), (C,A',B',C').
    Applying the 6 splits in the fixed order above to a working SET of
    pieces (starting from {tet}) grows that set by at most 3x per split;
    each output piece not yet known to be degenerate is re-classified by
    the next split from scratch.

    Collapsed-Gauss (Duffy / Stroud conical) tet rule. For a tet
    (P0,P1,P2,P3) and (u,v,w) in [0,1]^3,

      x(u,v,w) = P0 + u(P1-P0) + uv(P2-P1) + uvw(P3-P2)
      |J|      = |D| u^2 v,   D = det[P1-P0, P2-P0, P3-P0]
      weight   = wu wv ww u^2 v |D|                          (int u^2 v = 1/6)

    A monomial of total degree Dg in x has u-degree <= Dg, v-degree <= Dg,
    w-degree <= Dg; multiplying by the Jacobian u^2 v raises those to
    u <= Dg+2, v <= Dg+1, w <= Dg. Gauss-Legendre with m points is exact to
    degree 2m-1, hence m_u = ceil((Dg+3)/2), m_v = ceil((Dg+2)/2),
    m_w = ceil((Dg+1)/2). tetClipQuadrature uses Dg = volDeg = volDegree(p) = 6p:
    the (p+1)^3-node tensor-Gauss mass matrix of a degree-p background space
    has entries containing x^{2p} y^{2p} z^{2p}, total degree 6p, so 6p is
    the 3D mass-matrix parity degree. That gives
    (m_u,m_v,m_w) = (3p+2, 3p+1, 3p+1), i.e. 392 nodes per sub-tet at p=2.
    Using Dg = 4p instead is a one-line change of the single function
    `volDegree` below, nothing else: it would give
    (m_u,m_v,m_w) = (2p+2, 2p+1, 2p+1).

    Collapsed-Gauss triangle rule for the boundary. For a triangle
    (Q0,Q1,Q2) and (u,v) in [0,1]^2,

      x(u,v) = Q0 + u(Q1-Q0) + uv(Q2-Q1)
      weight = wu wv u |(Q1-Q0) x (Q2-Q0)|                   (int u = 1/2)

    By the same Jacobian argument, u-degree <= Dg+1, v-degree <= Dg, so
    m_u = ceil((Dg+2)/2), m_v = ceil((Dg+1)/2). The boundary rule uses
    Dg = bdrDeg = bdrDegree(p) = 6p, the same target degree as volDeg:
    boundary integrands built from products of two degree-p tensor
    functions contain x^{2p} y^{2p} z^{2p}, which restricted to a planar
    boundary triangle is a polynomial of total degree 6p in (u,v). This
    gives m_u = m_v = 3p+1.

    Reader restriction: ASCII gmsh MSH format 4.1 only ($MeshFormat line
    "4.1 0 8"), 4-node tetrahedra (element type 4) only; every other
    element type (surface triangles gmsh writes for the volume's own
    boundary mesh, edges, points) is skipped by line count, never parsed.
    Boundary triangles are derived from the tets themselves (a face shared
    by exactly one tet), never read from the file, and oriented outward
    using the owning tet's fourth vertex.

    Complexity: reading is O(N) in file size; the boundary derivation sorts
    4T face keys, O(T log T); the clipping loop is
    O(#tets x (n + overlapped cells x 6 splits x <=3 pieces)), n = cells
    per direction (colOf locates the cell range by a linear scan of the
    grid lines); the boundary loop is
    O(#boundary triangles x (n + overlapped cells x 6 clips)).

    Traps documented, not handled:
    - a boundary piece on a grid line lands in the LOWER cell along every
      axis (colOf below), even when that cell is Empty -- this is a
      legitimate degeneracy (the mesh face lies exactly on a knot plane),
      not a bug, and T9 of immersed_tetmesh_clip_quadrature_example
      (--case cube) exists specifically to prove it occurs on the shipped
      test data;
    - a Full cell (its clipped volume equals its box within the T4
      tolerance) KEEPS its clipped rule rather than being replaced by a
      tensor-Gauss rule: the clipped rule is already exact to volDeg, the
      Full label is purely informational, and swapping in a tensor rule
      would let a sliver below the classification tolerance silently
      change the integral;
    - the mesh bounding box must lie strictly inside the background box
      (GISMO_ENSURE); there is no support for a mesh touching or exceeding
      the background box;
    - quadrature never post-processes mesh coordinates: no snapping, no
      hand-editing of a .msh file. A geometry/grid combination that
      produces a genuine face-on-knot-plane degeneracy (--case cube) is
      exercised deliberately, not avoided.

    Reference:
    - M.J. Powell, T. Abel, "An exact general remeshing scheme applied to
      physically conservative voxelization", J. Comput. Phys. 297 (2015)
      340-356.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gismo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

/// Tet-mesh / Cartesian-grid exact clipping quadrature (see the file header).
namespace gsTetClip
{
using namespace gismo;

/// Formats \a v in std::scientific notation with \a prec digits after the decimal point (prec+1 significant digits).
inline std::string fmtSci(real_t v, int prec = 6)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

/// Neumaier's variant of Kahan compensated summation: see the 2D driver
/// (immersed_gauss_green_quadrature_example.cpp) for the rounding-error
/// argument. Used for every sum a check below compares.
class KahanSum
{
public:
    void add(real_t x)
    {
        const real_t t = m_sum + x;
        m_c += (math::abs(m_sum) >= math::abs(x)) ? ((m_sum-t)+x) : ((x-t)+m_sum);
        m_sum = t;
    }
    real_t value() const { return m_sum + m_c; }
private:
    real_t m_sum = 0, m_c = 0;
};

/// e = |computed-exact| / max(|exact|, refScale); refScale is the oracle
/// volume (volume checks) or area (boundary checks), so a moment that
/// vanishes exactly still has a well-defined scaled error.
inline real_t scaledErr(real_t computed, real_t exact, real_t refScale)
{ return math::abs(computed - exact) / math::max(math::abs(exact), refScale); }

/// ceil(x/2) for a non-negative integer x.
inline index_t ceilHalf(index_t x) { return (x+1)/2; }

/// Collapsed-tet volume-rule exactness degree, the single edit point for the
/// 4p/6p choice discussed in the file header (change the body only, e.g. to
/// `return 4*p;`).
inline index_t volDegree(index_t p) { return 6*p; }
/// Collapsed-triangle boundary-rule exactness degree (same parity argument
/// as volDegree, see the file header).
inline index_t bdrDegree(index_t p) { return 6*p; }

/// Exponent triples (a,b,c) with a+b+c <= 2p, the moment set T2/T3 check
/// against the volume oracle and the analytic cube moments.
inline std::vector<std::array<int,3> > momentSet3(index_t p)
{
    std::vector<std::array<int,3> > s;
    for (int a = 0; a <= 2*p; ++a)
        for (int b = 0; b <= 2*p-a; ++b)
            for (int c = 0; c <= 2*p-a-b; ++c)
                s.push_back({a,b,c});
    return s;
}

//----------------------------------------------------------------------------
// Small 3-vector helpers on std::array<real_t,3> (kept free-function, not an
// operator-overloaded class, to match the amount of arithmetic actually
// needed here without hiding the collapsed-map formulas behind operators).
//----------------------------------------------------------------------------

typedef std::array<real_t,3> Vec3;

inline Vec3 sub3(const Vec3 & A, const Vec3 & B)
{ return {A[0]-B[0], A[1]-B[1], A[2]-B[2]}; }
inline Vec3 add3(const Vec3 & A, const Vec3 & B)
{ return {A[0]+B[0], A[1]+B[1], A[2]+B[2]}; }
inline Vec3 scale3(real_t s, const Vec3 & A)
{ return {s*A[0], s*A[1], s*A[2]}; }
inline real_t dot3(const Vec3 & A, const Vec3 & B)
{ return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; }
inline Vec3 cross3(const Vec3 & A, const Vec3 & B)
{ return { A[1]*B[2]-A[2]*B[1], A[2]*B[0]-A[0]*B[2], A[0]*B[1]-A[1]*B[0] }; }
/// Signed volume x6 of the tet spanned by the origin and A,B,C (used with
/// A,B,C already relative to a common vertex, e.g. P1-P0,P2-P0,P3-P0).
inline real_t det3(const Vec3 & A, const Vec3 & B, const Vec3 & C)
{ return dot3(A, cross3(B,C)); }

//----------------------------------------------------------------------------
// Data structures.
//----------------------------------------------------------------------------

/// A tetrahedral volume mesh with its derived outward-oriented boundary
/// triangles, read by readMsh41().
struct TetMesh
{
    std::vector<Vec3>                    P;         // node coordinates (compact index)
    std::vector<std::array<index_t,4> >  tet;        // compact node indices
    std::vector<std::array<index_t,4> >  bdrTri;     // a,b,c (outward-oriented), opp = owner tet's 4th vertex
    std::vector<index_t>                 bdrOwner;   // owner tet of each boundary triangle
    index_t nNodesFile = 0;                          // node count reported by $Nodes
    real_t lo[3], hi[3];                             // bbox of tet vertices
};

/// Uniform Cartesian background grid, n cells per direction; cell id =
/// i + n*(j + n*k).
struct Grid3      { real_t x0, y0, z0, h; index_t n; };            // cell id = i + n*(j + n*k)
/// Per-cell volume quadrature rule: 3 x m nodes (physical coordinates), m
/// weights.
struct CellRule3  { gsMatrix<real_t> nodes; gsVector<real_t> weights; };                       // 3 x m, m
/// Per-cell boundary quadrature rule: 3 x m nodes, m weights, 3 x m outward
/// unit normals, all in physical coordinates.
struct CellBdrRule3 { gsMatrix<real_t> nodes; gsVector<real_t> weights; gsMatrix<real_t> normals; }; // 3 x m, m, 3 x m

/// Classification of a background cell against the tet mesh.
enum CellStatus { Empty = -1, Cut = 0, Full = 1 };

/// The full clip quadrature: one CellRule3/CellBdrRule3 per background cell
/// (size n^3 each), plus the per-cell CellStatus.
struct TetClipQuadrature
{
    Grid3 grid;
    std::vector<int>         status;   // CellStatus per cell id
    std::vector<CellRule3>   vol;      // size n^3
    std::vector<CellBdrRule3> bdr;     // size n^3
};

/// Per-run diagnostics of tetClipQuadrature(): piece/cell counts and, per
/// tet, its clipped-vs-exact volume and diameter (used by T1).
struct TetClipStats
{
    index_t nTets = 0, nBdrTris = 0, nPieces = 0, nBdrPieces = 0, nCut = 0, nFull = 0, nEmpty = 0;
    std::vector<real_t> tetClippedVol, tetVol, tetDiam;   // per tet
};

/// Column index of a point along one axis, given that axis's own precomputed
/// grid-line positions \a X (size n+1, X[i] = x0+i*h): the number of interior
/// lines X[1..n-1] with x > X[i] -- never a floor or a tolerance, so a point
/// exactly ON a line always belongs to the LOWER cell (i-1). \a X must be the
/// SAME array every other predicate in this file compares against (computed
/// once per grid), not a fresh x0+i*h at each call site, since two
/// independently rounded sums can disagree by an ULP and silently break the
/// exact-partition guarantees the clipping and classification rely on. Used
/// for all three axes (x, y and z alike).
inline index_t colOf(real_t x, const std::vector<real_t> & X)
{
    index_t c = 0;
    for (size_t i = 1; i + 1 < X.size(); ++i)
        if (x > X[i]) ++c;
    return c;
}

/// Fills the n+1 grid-line positions along each axis, X[i]/Y[i]/Z[i] =
/// x0/y0/z0 + i*h.
inline void gridLines(const Grid3 & grid, std::vector<real_t> & X, std::vector<real_t> & Y, std::vector<real_t> & Z)
{
    X.resize(grid.n+1); Y.resize(grid.n+1); Z.resize(grid.n+1);
    for (index_t i = 0; i <= grid.n; ++i)
    { X[i] = grid.x0 + i*grid.h; Y[i] = grid.y0 + i*grid.h; Z[i] = grid.z0 + i*grid.h; }
}

/// Largest absolute background-box coordinate, the "L" of every rounding
/// bound below (the classification tolerance, T1, T4): an ULP at that
/// magnitude is the floor no exact-arithmetic argument can beat once a
/// clipped point's coordinate is computed by interpolation.
inline real_t boxCornerAbsMax(const Grid3 & grid)
{
    const index_t n = grid.n;
    return math::max(math::max(math::abs(grid.x0), math::abs(grid.x0 + n*grid.h)),
           math::max(math::max(math::abs(grid.y0), math::abs(grid.y0 + n*grid.h)),
                     math::max(math::abs(grid.z0), math::abs(grid.z0 + n*grid.h))));
}

/// Kahan-summed total weight of a CellRule3, its cell volume.
inline real_t cellVolume(const CellRule3 & r)
{
    KahanSum s;
    for (index_t k = 0; k != r.weights.size(); ++k) s.add(r.weights[k]);
    return s.value();
}

//----------------------------------------------------------------------------
// Exact-at-plane intersection (shared by the tet split and the boundary
// polygon clip).
//----------------------------------------------------------------------------

/// Intersection of segment (A,B) with the axis-aligned plane
/// coordinate[axis] = c, EXACT at that coordinate (set directly to \a c, not
/// recomputed through the affine formula) and interpolating the other two
/// coordinates from the endpoint with the smaller coordinate along \a axis,
/// so the result does not depend on the segment's direction -- the 3D
/// analogue of intersectX/intersectY in the 2D driver. If either endpoint
/// already lies exactly on the plane it is returned verbatim: on a mixed
/// edge (exactly one endpoint "in" a closed half-space, the other strictly
/// "out") only the "in" endpoint can equal \a c, so this is unambiguous.
/// Without the early return, an "in" endpoint that happens to be the one
/// with the LARGER coordinate along \a axis makes the general formula's own
/// t evaluate to exactly 1 (a shared subexpression divided by itself), and
/// S + 1.0*(E-S) is not guaranteed to round back to E bit-for-bit, even
/// though mathematically it must equal E; the early return sidesteps that
/// roundoff entirely. Two cells sharing the plane get bitwise-identical
/// points only when both clip the same edge (both then interpolate from the
/// same pair of endpoints); once an earlier split has already shortened the
/// edge differently in each cell, the two computed points agree to
/// rounding only -- an O(eps) gap or overlap that T1's tolAbs term absorbs.
inline Vec3 intersectPlane(int axis, real_t c, const Vec3 & A, const Vec3 & B)
{
    if (A[axis] == c) return A;
    if (B[axis] == c) return B;
    const bool AisS = (A[axis] <= B[axis]);
    const Vec3 & S = AisS ? A : B;
    const Vec3 & E = AisS ? B : A;
    const real_t t = (c - S[axis]) / (E[axis] - S[axis]);
    Vec3 R;
    for (int d = 0; d != 3; ++d)
        R[d] = (d == axis) ? c : S[d] + t*(E[d]-S[d]);
    return R;
}

//----------------------------------------------------------------------------
// Tet x half-space split (exact vertex classification, Section "Vertex
// classification and split cases" above).
//----------------------------------------------------------------------------

typedef std::array<Vec3,4> Tet4;

/// Splits tet \a T by the closed half-space coordinate[axis] >= c (\a lowerSide
/// true, a cell's lower face) or coordinate[axis] <= c (false) and returns the
/// 0, 1 or 3 sub-tets covering T ^ half-space (split cases in the file header).
inline std::vector<Tet4> splitTetByHalfspace(const Tet4 & T, int axis, real_t c, bool lowerSide)
{
    bool in[4];
    for (int i = 0; i != 4; ++i)
        in[i] = lowerSide ? (T[i][axis] >= c) : (T[i][axis] <= c);
    const int nIn = (int)in[0] + (int)in[1] + (int)in[2] + (int)in[3];

    std::vector<Tet4> out;
    if (4 == nIn) { out.push_back(T); return out; }
    if (0 == nIn) return out;

    auto I = [&](int i, int j) { return intersectPlane(axis, c, T[i], T[j]); };

    if (1 == nIn)
    {
        int a = -1;
        for (int i = 0; i != 4; ++i) if (in[i]) a = i;
        int o[3]; int k = 0;
        for (int i = 0; i != 4; ++i) if (i != a) o[k++] = i;
        out.push_back({ T[a], I(a,o[0]), I(a,o[1]), I(a,o[2]) });
        return out;
    }
    if (3 == nIn)
    {
        int d = -1;
        for (int i = 0; i != 4; ++i) if (!in[i]) d = i;
        int o[3]; int k = 0;
        for (int i = 0; i != 4; ++i) if (i != d) o[k++] = i;
        const int a = o[0], b = o[1], cc = o[2];
        const Vec3 Ap = T[a], Bp = T[b], Cp = T[cc];
        const Vec3 Aq = I(a,d), Bq = I(b,d), Cq = I(cc,d);
        out.push_back({Ap, Bp, Cp, Aq});
        out.push_back({Bp, Cp, Aq, Bq});
        out.push_back({Cp, Aq, Bq, Cq});
        return out;
    }
    // nIn == 2
    int inIdx[2], outIdx[2]; int ki = 0, ko = 0;
    for (int i = 0; i != 4; ++i) { if (in[i]) inIdx[ki++] = i; else outIdx[ko++] = i; }
    const int a = inIdx[0], b = inIdx[1], cc = outIdx[0], dd = outIdx[1];
    const Vec3 Ap = T[a], Bp = T[b];
    const Vec3 Aq1 = I(a,cc), Aq2 = I(a,dd);
    const Vec3 Bq1 = I(b,cc), Bq2 = I(b,dd);
    out.push_back({Ap, Aq1, Aq2, Bp});
    out.push_back({Aq1, Aq2, Bp, Bq1});
    out.push_back({Aq2, Bp, Bq1, Bq2});
    return out;
}

/// Applies the 6 fixed splits (x>=xL, x<=xR, y>=yB, y<=yT, z>=zN, z<=zF) in
/// order to the working set of pieces, starting from {T}.
inline std::vector<Tet4> clipTetToCell(const Tet4 & T, real_t xL, real_t xR, real_t yB, real_t yT,
                                real_t zN, real_t zF)
{
    std::vector<Tet4> pieces(1, T);
    auto applySplit = [&pieces](int axis, real_t c, bool lowerSide)
    {
        std::vector<Tet4> next;
        next.reserve(pieces.size()*3);
        for (const Tet4 & piece : pieces)
        {
            const std::vector<Tet4> sub = splitTetByHalfspace(piece, axis, c, lowerSide);
            next.insert(next.end(), sub.begin(), sub.end());
        }
        pieces.swap(next);
    };
    applySplit(0, xL, true);  applySplit(0, xR, false);
    applySplit(1, yB, true);  applySplit(1, yT, false);
    applySplit(2, zN, true);  applySplit(2, zF, false);
    return pieces;
}

//----------------------------------------------------------------------------
// Boundary triangle x box clip (Sutherland-Hodgman, 3D points, axis/plane/
// side parameters instead of std::functions).
//----------------------------------------------------------------------------

typedef std::vector<Vec3> Poly3;

/// One Sutherland-Hodgman half-space clip of the (possibly non-triangular)
/// polygon \a in against coordinate[axis] {>=,<=} c.
inline Poly3 clipHalfSpace(const Poly3 & in, int axis, real_t c, bool lowerSide)
{
    const size_t m = in.size();
    Poly3 out;
    if (0 == m) return out;
    out.reserve(m+1);
    for (size_t k = 0; k != m; ++k)
    {
        const Vec3 & S = in[k];
        const Vec3 & E = in[(k+1)%m];
        const bool sIn = lowerSide ? (S[axis] >= c) : (S[axis] <= c);
        const bool eIn = lowerSide ? (E[axis] >= c) : (E[axis] <= c);
        if (sIn) out.push_back(S);
        if (sIn != eIn) out.push_back(intersectPlane(axis, c, S, E));
    }
    return out;
}

/// Clips the triangle (A,B,C) against the cell box by the 6 fixed
/// half-space clips, Sutherland-Hodgman.
inline Poly3 clipTriangleToBox3(const Vec3 & A, const Vec3 & B, const Vec3 & C,
                         real_t xL, real_t xR, real_t yB, real_t yT, real_t zN, real_t zF)
{
    Poly3 poly = {A, B, C};
    poly = clipHalfSpace(poly, 0, xL, true);
    poly = clipHalfSpace(poly, 0, xR, false);
    poly = clipHalfSpace(poly, 1, yB, true);
    poly = clipHalfSpace(poly, 1, yT, false);
    poly = clipHalfSpace(poly, 2, zN, true);
    poly = clipHalfSpace(poly, 2, zF, false);
    return poly;
}

//----------------------------------------------------------------------------
// 1D Gauss nodes mapped to [0,1] (the map used by every collapsed rule
// below): u = 0.5*(1+z), w = 0.5*wz.
//----------------------------------------------------------------------------

/// A 1D quadrature rule on [0,1], nodes \a u and weights \a w.
struct Rule1D { std::vector<real_t> u, w; };

/// Gauss-Legendre rule with \a m points, mapped from [-1,1] to [0,1].
inline Rule1D gauss01(index_t m)
{
    gsGaussRule<real_t> g(m);
    const gsMatrix<real_t> & z = g.referenceNodes();
    const gsVector<real_t> & wz = g.referenceWeights();
    Rule1D r; r.u.resize(m); r.w.resize(m);
    for (index_t k = 0; k != m; ++k) { r.u[k] = 0.5*(1.0+z(0,k)); r.w[k] = 0.5*wz[k]; }
    return r;
}

/// Appends the collapsed-Gauss tet rule for the already-known-nondegenerate
/// tet (P0,P1,P2,P3), |D| = absD, to the per-cell node/weight vectors.
inline void addCollapsedTetNodes(const Vec3 & P0, const Vec3 & P1, const Vec3 & P2, const Vec3 & P3,
                          real_t absD, const Rule1D & Ru, const Rule1D & Rv, const Rule1D & Rw,
                          std::vector<real_t> & vx, std::vector<real_t> & vy, std::vector<real_t> & vz,
                          std::vector<real_t> & vw)
{
    const Vec3 e1 = sub3(P1,P0), e2 = sub3(P2,P1), e3 = sub3(P3,P2);
    for (size_t iu = 0; iu != Ru.u.size(); ++iu)
    {
        const real_t u = Ru.u[iu];
        const Vec3 Pu = add3(P0, scale3(u, e1));
        for (size_t iv = 0; iv != Rv.u.size(); ++iv)
        {
            const real_t v = Rv.u[iv];
            const Vec3 Puv = add3(Pu, scale3(u*v, e2));
            const real_t wuv = Ru.w[iu]*Rv.w[iv]*u*u*v*absD;
            for (size_t iw = 0; iw != Rw.u.size(); ++iw)
            {
                const real_t w = Rw.u[iw];
                const Vec3 Pt = add3(Puv, scale3(u*v*w, e3));
                vx.push_back(Pt[0]); vy.push_back(Pt[1]); vz.push_back(Pt[2]);
                vw.push_back(wuv*Rw.w[iw]);
            }
        }
    }
}

/// Appends the collapsed-Gauss triangle rule for the already-known-
/// nondegenerate fan triangle (Q0,Q1,Q2), \a mag = |(Q1-Q0)x(Q2-Q0)|, to the
/// per-cell boundary node/weight/normal vectors. \a normal is the WHOLE
/// boundary triangle's own outward unit normal, constant across every fan
/// piece of it.
inline void addCollapsedTriNodes(const Vec3 & Q0, const Vec3 & Q1, const Vec3 & Q2, real_t mag,
                          const Rule1D & Ru, const Rule1D & Rv, const Vec3 & normal,
                          std::vector<real_t> & bx, std::vector<real_t> & by, std::vector<real_t> & bz,
                          std::vector<real_t> & bw, std::vector<real_t> & bnx, std::vector<real_t> & bny,
                          std::vector<real_t> & bnz)
{
    const Vec3 e1 = sub3(Q1,Q0), e2 = sub3(Q2,Q1);
    for (size_t iu = 0; iu != Ru.u.size(); ++iu)
    {
        const real_t u = Ru.u[iu];
        const Vec3 Pu = add3(Q0, scale3(u, e1));
        for (size_t iv = 0; iv != Rv.u.size(); ++iv)
        {
            const real_t v = Rv.u[iv];
            const Vec3 Pt = add3(Pu, scale3(u*v, e2));
            bx.push_back(Pt[0]); by.push_back(Pt[1]); bz.push_back(Pt[2]);
            bw.push_back(Ru.w[iu]*Rv.w[iv]*u*mag);
            bnx.push_back(normal[0]); bny.push_back(normal[1]); bnz.push_back(normal[2]);
        }
    }
}

//----------------------------------------------------------------------------
// ASCII gmsh MSH 4.1 reader + boundary derivation.
//----------------------------------------------------------------------------

/// Reads an ASCII gmsh MSH 4.1 file: 4-node tetrahedra (element type 4)
/// only, every other element block skipped by line count. Derives the
/// boundary triangles from the tets (a face owned by exactly one tet,
/// O(T log T) by sorting the 4T face keys), oriented outward using the
/// owner tet's fourth vertex.
inline TetMesh readMsh41(const std::string & file)
{
    std::ifstream in(file.c_str());
    GISMO_ENSURE(in.good(), "readMsh41: cannot open " << file);

    std::string line;
    while (std::getline(in, line) && line != "$MeshFormat") { }
    GISMO_ENSURE(!in.eof(), "readMsh41: no $MeshFormat section in " << file);
    {
        std::string verLine;
        std::getline(in, verLine);
        std::istringstream iss(verLine);
        real_t version = 0; int filetype = -1, datasize = 0;
        iss >> version >> filetype >> datasize;
        GISMO_ENSURE(version == 4.1, "readMsh41: " << file << " has MeshFormat version "
                    << version << " (need 4.1); regenerate with -format msh41 and Mesh.Binary = 0");
        GISMO_ENSURE(filetype == 0, "readMsh41: " << file << " is a binary MSH file; "
                    "regenerate with -format msh41 and Mesh.Binary = 0");
        std::getline(in, line); // $EndMeshFormat
    }

    TetMesh M;
    std::vector<index_t> tag2idx;

    while (std::getline(in, line))
    {
        if ("$Nodes" == line)
        {
            std::string hdr; std::getline(in, hdr);
            std::istringstream hiss(hdr);
            long numEntityBlocks = 0, numNodes = 0, minTag = 0, maxTag = 0;
            hiss >> numEntityBlocks >> numNodes >> minTag >> maxTag;
            tag2idx.assign((size_t)maxTag+1, (index_t)-1);
            M.P.assign((size_t)numNodes, Vec3());
            M.nNodesFile = (index_t)numNodes;

            index_t nextIdx = 0;
            for (long b = 0; b != numEntityBlocks; ++b)
            {
                std::string blkHdr; std::getline(in, blkHdr);
                std::istringstream biss(blkHdr);
                long dim = 0, tag = 0, param = 0, cnt = 0;
                biss >> dim >> tag >> param >> cnt;
                GISMO_ENSURE(0 == param, "readMsh41: parametric node block in " << file
                            << " is not supported");
                std::vector<long> tags((size_t)cnt);
                for (long kk = 0; kk != cnt; ++kk)
                { std::string tl; std::getline(in, tl); tags[kk] = std::stol(tl); }
                for (long kk = 0; kk != cnt; ++kk)
                {
                    std::string cl; std::getline(in, cl);
                    std::istringstream ciss(cl);
                    real_t x = 0, y = 0, z = 0; ciss >> x >> y >> z;
                    M.P[nextIdx] = {x,y,z};
                    tag2idx[tags[kk]] = nextIdx;
                    ++nextIdx;
                }
            }
            std::getline(in, line); // $EndNodes
            continue;
        }
        if ("$Elements" == line)
        {
            std::string hdr; std::getline(in, hdr);
            std::istringstream hiss(hdr);
            long numEntityBlocks = 0, numElements = 0, minTag = 0, maxTag = 0;
            hiss >> numEntityBlocks >> numElements >> minTag >> maxTag;
            for (long b = 0; b != numEntityBlocks; ++b)
            {
                std::string blkHdr; std::getline(in, blkHdr);
                std::istringstream biss(blkHdr);
                long dim = 0, tag = 0, etype = 0, cnt = 0;
                biss >> dim >> tag >> etype >> cnt;
                if (4 != etype)
                {
                    for (long kk = 0; kk != cnt; ++kk) std::getline(in, line);
                    continue;
                }
                for (long kk = 0; kk != cnt; ++kk)
                {
                    std::string el; std::getline(in, el);
                    std::istringstream eiss(el);
                    long etag = 0, n0 = 0, n1 = 0, n2 = 0, n3 = 0;
                    eiss >> etag >> n0 >> n1 >> n2 >> n3;
                    GISMO_ENSURE(n0 < (long)tag2idx.size() && n1 < (long)tag2idx.size() &&
                                n2 < (long)tag2idx.size() && n3 < (long)tag2idx.size() &&
                                tag2idx[n0] >= 0 && tag2idx[n1] >= 0 && tag2idx[n2] >= 0 && tag2idx[n3] >= 0,
                                "readMsh41: element " << etag << " in " << file << " references an unknown node tag");
                    M.tet.push_back({tag2idx[n0], tag2idx[n1], tag2idx[n2], tag2idx[n3]});
                }
            }
            std::getline(in, line); // $EndElements
            continue;
        }
        if (!line.empty() && '$' == line[0] && line.substr(0,4) != "$End")
        {
            const std::string endTag = "$End" + line.substr(1);
            while (std::getline(in, line) && line != endTag) { }
        }
    }

    GISMO_ENSURE(!M.tet.empty(), "readMsh41: no type-4 tetrahedra found in " << file);

    for (size_t t = 0; t != M.tet.size(); ++t)
    {
        const std::array<index_t,4> & tt = M.tet[t];
        const real_t D = det3(sub3(M.P[tt[1]],M.P[tt[0]]), sub3(M.P[tt[2]],M.P[tt[0]]), sub3(M.P[tt[3]],M.P[tt[0]]));
        GISMO_ENSURE(0.0 != D, "readMsh41: tet " << t << " in " << file << " is flat (zero volume)");
    }

    for (int d = 0; d != 3; ++d) { M.lo[d] = std::numeric_limits<real_t>::infinity(); M.hi[d] = -M.lo[d]; }
    for (const std::array<index_t,4> & tt : M.tet)
        for (int v = 0; v != 4; ++v)
            for (int d = 0; d != 3; ++d)
            { M.lo[d] = math::min(M.lo[d], M.P[tt[v]][d]); M.hi[d] = math::max(M.hi[d], M.P[tt[v]][d]); }

    // Boundary derivation: face omitting local vertex i, key = sorted node
    // indices, plus (owner tet, opp = tet's local vertex i).
    struct FaceRec { std::array<index_t,3> key, tri; index_t owner, opp; };
    std::vector<FaceRec> faces;
    faces.reserve(M.tet.size()*4);
    for (size_t t = 0; t != M.tet.size(); ++t)
    {
        const std::array<index_t,4> & tt = M.tet[t];
        for (int i = 0; i != 4; ++i)
        {
            std::array<index_t,3> f; int k = 0;
            for (int j = 0; j != 4; ++j) if (j != i) f[k++] = tt[j];
            std::array<index_t,3> key = f;
            std::sort(key.begin(), key.end());
            faces.push_back({key, f, (index_t)t, tt[i]});
        }
    }
    std::sort(faces.begin(), faces.end(),
             [](const FaceRec & a, const FaceRec & b) { return a.key < b.key; });

    size_t idx = 0;
    while (idx != faces.size())
    {
        size_t j = idx;
        while (j != faces.size() && faces[j].key == faces[idx].key) ++j;
        const size_t mult = j - idx;
        GISMO_ENSURE(mult <= 2, "readMsh41: non-manifold face in " << file);
        if (1 == mult)
        {
            const FaceRec & fr = faces[idx];
            index_t a = fr.tri[0], b = fr.tri[1], c = fr.tri[2];
            const Vec3 & Pa = M.P[a], & Pb = M.P[b], & Pc = M.P[c], & Po = M.P[fr.opp];
            const real_t orient = dot3(cross3(sub3(Pb,Pa), sub3(Pc,Pa)), sub3(Po,Pa));
            GISMO_ENSURE(0.0 != orient, "readMsh41: degenerate boundary-face orientation in " << file);
            if (orient > 0.0) std::swap(b,c);
            M.bdrTri.push_back({a,b,c,fr.opp});
            M.bdrOwner.push_back(fr.owner);
        }
        idx = j;
    }

    return M;
}

//----------------------------------------------------------------------------
// Volume + boundary clip quadrature.
//----------------------------------------------------------------------------

/// Exact cell-by-cell tet-clip quadrature of the tet mesh \a M against the
/// uniform background grid \a grid. Volume rule: collapsed-Gauss tet rule,
/// degree volDeg = volDegree(p) (6p by default), on every sub-tet surviving
/// the 6 half-space splits of every (tet, overlapped cell) pair -- each
/// sub-tet carries m_u*m_v*m_w = (3p+2)(3p+1)(3p+1) nodes (392 at p=2).
/// Boundary rule: Sutherland-Hodgman clip of every boundary triangle
/// against every overlapped cell box, fan-triangulated, integrated with the
/// collapsed-Gauss triangle rule, degree bdrDeg = bdrDegree(p) (6p by
/// default), carrying the owner triangle's own outward unit normal -- each
/// fan triangle carries m_u*m_v = (3p+1)^2 nodes (49 at p=2).
///
/// Complexity: O(#tets x (n + overlapped cells per tet x 6 splits x <=3
/// pieces)) for the volume loop, plus O(#boundary triangles x (n +
/// overlapped cells per triangle x 6 clips)) for the boundary loop, n =
/// cells per direction: colOf locates the cell range by a linear scan of
/// the grid lines, O(n) per tet/triangle; each split/clip is O(1) (a tet
/// split produces at most 3 pieces, a triangle clipped by 6 half-spaces has
/// at most 9 vertices). See the file header for the split cases and the
/// degree argument.
inline TetClipQuadrature tetClipQuadrature(const TetMesh & M, const Grid3 & grid, index_t p,
                                    TetClipStats & stats)
{
    stats = TetClipStats();
    const index_t n = grid.n;
    std::vector<real_t> X, Y, Z;
    gridLines(grid, X, Y, Z);

    const index_t volDeg = volDegree(p);
    const Rule1D Ru = gauss01(ceilHalf(volDeg+3));
    const Rule1D Rv = gauss01(ceilHalf(volDeg+2));
    const Rule1D Rw = gauss01(ceilHalf(volDeg+1));

    const index_t bdrDeg = bdrDegree(p);
    const Rule1D Rbu = gauss01(ceilHalf(bdrDeg+2));
    const Rule1D Rbv = gauss01(ceilHalf(bdrDeg+1));

    const size_t nCells = (size_t)n*n*n;
    std::vector<std::vector<real_t> > vx(nCells), vy(nCells), vz(nCells), vw(nCells);
    std::vector<std::vector<real_t> > bx(nCells), by(nCells), bz(nCells), bw(nCells),
                                       bnx(nCells), bny(nCells), bnz(nCells);
    std::vector<char> touched(nCells, 0);

    stats.nTets = (index_t)M.tet.size();
    stats.nBdrTris = (index_t)M.bdrTri.size();
    stats.tetClippedVol.assign(M.tet.size(), 0.0);
    stats.tetVol.assign(M.tet.size(), 0.0);
    stats.tetDiam.assign(M.tet.size(), 0.0);

    for (size_t t = 0; t != M.tet.size(); ++t)
    {
        const std::array<index_t,4> & tt = M.tet[t];
        const Tet4 T = { M.P[tt[0]], M.P[tt[1]], M.P[tt[2]], M.P[tt[3]] };
        const real_t D = det3(sub3(T[1],T[0]), sub3(T[2],T[0]), sub3(T[3],T[0]));
        stats.tetVol[t] = math::abs(D)/6.0;

        real_t diam = 0;
        for (int i = 0; i != 4; ++i)
            for (int j = i+1; j != 4; ++j)
                diam = math::max(diam, math::sqrt(dot3(sub3(T[j],T[i]), sub3(T[j],T[i]))));
        stats.tetDiam[t] = diam;

        real_t xmin=T[0][0], xmax=T[0][0], ymin=T[0][1], ymax=T[0][1], zmin=T[0][2], zmax=T[0][2];
        for (int v = 1; v != 4; ++v)
        {
            xmin = math::min(xmin, T[v][0]); xmax = math::max(xmax, T[v][0]);
            ymin = math::min(ymin, T[v][1]); ymax = math::max(ymax, T[v][1]);
            zmin = math::min(zmin, T[v][2]); zmax = math::max(zmax, T[v][2]);
        }
        const index_t i0 = colOf(xmin,X), i1 = colOf(xmax,X);
        const index_t j0 = colOf(ymin,Y), j1 = colOf(ymax,Y);
        const index_t k0 = colOf(zmin,Z), k1 = colOf(zmax,Z);

        KahanSum tetClipped;
        for (index_t k = k0; k <= k1; ++k)
        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
        {
            const std::vector<Tet4> pieces = clipTetToCell(T, X[i],X[i+1], Y[j],Y[j+1], Z[k],Z[k+1]);
            const size_t id = (size_t)i + (size_t)n*((size_t)j + (size_t)n*(size_t)k);
            for (const Tet4 & piece : pieces)
            {
                const real_t Dp = det3(sub3(piece[1],piece[0]), sub3(piece[2],piece[0]), sub3(piece[3],piece[0]));
                if (0.0 == Dp) continue;
                const real_t absDp = math::abs(Dp);
                ++stats.nPieces;
                touched[id] = 1;
                tetClipped.add(absDp/6.0);
                addCollapsedTetNodes(piece[0],piece[1],piece[2],piece[3], absDp, Ru,Rv,Rw,
                                    vx[id],vy[id],vz[id],vw[id]);
            }
        }
        stats.tetClippedVol[t] = tetClipped.value();
    }

    for (size_t bt = 0; bt != M.bdrTri.size(); ++bt)
    {
        const std::array<index_t,4> & bTri = M.bdrTri[bt];
        const Vec3 & A = M.P[bTri[0]], & B = M.P[bTri[1]], & C = M.P[bTri[2]];
        const Vec3 crAB = cross3(sub3(B,A), sub3(C,A));
        const real_t triMag = math::sqrt(dot3(crAB,crAB));
        const Vec3 normal = scale3(1.0/triMag, crAB);

        const real_t xmin = math::min(A[0], math::min(B[0],C[0])), xmax = math::max(A[0], math::max(B[0],C[0]));
        const real_t ymin = math::min(A[1], math::min(B[1],C[1])), ymax = math::max(A[1], math::max(B[1],C[1]));
        const real_t zmin = math::min(A[2], math::min(B[2],C[2])), zmax = math::max(A[2], math::max(B[2],C[2]));
        const index_t i0 = colOf(xmin,X), i1 = colOf(xmax,X);
        const index_t j0 = colOf(ymin,Y), j1 = colOf(ymax,Y);
        const index_t k0 = colOf(zmin,Z), k1 = colOf(zmax,Z);

        for (index_t k = k0; k <= k1; ++k)
        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
        {
            const Poly3 poly = clipTriangleToBox3(A,B,C, X[i],X[i+1], Y[j],Y[j+1], Z[k],Z[k+1]);
            if (poly.size() < 3) continue;
            const size_t id = (size_t)i + (size_t)n*((size_t)j + (size_t)n*(size_t)k);
            for (size_t f = 1; f+1 < poly.size(); ++f)
            {
                const Vec3 & Q0 = poly[0], & Q1 = poly[f], & Q2 = poly[f+1];
                const Vec3 cr = cross3(sub3(Q1,Q0), sub3(Q2,Q0));
                if (0.0 == cr[0] && 0.0 == cr[1] && 0.0 == cr[2]) continue;
                const real_t mag = math::sqrt(dot3(cr,cr));
                ++stats.nBdrPieces;
                addCollapsedTriNodes(Q0,Q1,Q2, mag, Rbu,Rbv, normal,
                                    bx[id],by[id],bz[id],bw[id], bnx[id],bny[id],bnz[id]);
            }
        }
    }

    const real_t L = boxCornerAbsMax(grid);
    const real_t eps = std::numeric_limits<real_t>::epsilon();
    const real_t hh = grid.h*grid.h;
    const real_t tolAbs = 16*eps*L*hh;
    const real_t box = grid.h*hh;

    TetClipQuadrature Q;
    Q.grid = grid;
    Q.status.assign(nCells, (int)Empty);
    Q.vol.assign(nCells, CellRule3());
    Q.bdr.assign(nCells, CellBdrRule3());

    for (size_t id = 0; id != nCells; ++id)
    {
        CellRule3 vr;
        vr.nodes.resize(3, (index_t)vx[id].size());
        vr.weights.resize((index_t)vx[id].size());
        for (size_t kk = 0; kk != vx[id].size(); ++kk)
        {
            vr.nodes(0,(index_t)kk) = vx[id][kk]; vr.nodes(1,(index_t)kk) = vy[id][kk]; vr.nodes(2,(index_t)kk) = vz[id][kk];
            vr.weights[(index_t)kk] = vw[id][kk];
        }
        Q.vol[id] = vr;

        CellBdrRule3 br;
        br.nodes.resize(3, (index_t)bx[id].size());
        br.weights.resize((index_t)bx[id].size());
        br.normals.resize(3, (index_t)bx[id].size());
        for (size_t kk = 0; kk != bx[id].size(); ++kk)
        {
            br.nodes(0,(index_t)kk) = bx[id][kk]; br.nodes(1,(index_t)kk) = by[id][kk]; br.nodes(2,(index_t)kk) = bz[id][kk];
            br.weights[(index_t)kk] = bw[id][kk];
            br.normals(0,(index_t)kk) = bnx[id][kk]; br.normals(1,(index_t)kk) = bny[id][kk]; br.normals(2,(index_t)kk) = bnz[id][kk];
        }
        Q.bdr[id] = br;

        if (!touched[id]) { Q.status[id] = Empty; ++stats.nEmpty; continue; }

        const real_t vol = cellVolume(Q.vol[id]);
        if (vol >= box-(1e-14*box+tolAbs) && vol <= box*(1.0+1e-14)+tolAbs)
        { Q.status[id] = Full; ++stats.nFull; }
        else
        { Q.status[id] = Cut; ++stats.nCut; }
    }

    return Q;
}

//----------------------------------------------------------------------------
// Independent-code-path oracles (reversed vertex order, one extra Gauss
// point per direction), unclipped mesh integrals.
//----------------------------------------------------------------------------

/// Exact whole-mesh volume, Kahan-summed per-tet |det|/6.
inline real_t meshVolumeExact(const TetMesh & M)
{
    KahanSum s;
    for (const std::array<index_t,4> & tt : M.tet)
        s.add(math::abs(det3(sub3(M.P[tt[1]],M.P[tt[0]]), sub3(M.P[tt[2]],M.P[tt[0]]), sub3(M.P[tt[3]],M.P[tt[0]])))/6.0);
    return s.value();
}

/// Exact whole-mesh boundary area, Kahan-summed per-triangle 0.5*|cross|.
inline real_t unclippedBoundaryArea(const TetMesh & M)
{
    KahanSum s;
    for (const std::array<index_t,4> & tri : M.bdrTri)
    {
        const Vec3 cr = cross3(sub3(M.P[tri[1]],M.P[tri[0]]), sub3(M.P[tri[2]],M.P[tri[0]]));
        s.add(0.5*math::sqrt(dot3(cr,cr)));
    }
    return s.value();
}

/// Volume moment int x^a y^b z^c dV over the WHOLE (unclipped) mesh, by the
/// same collapsed-tet map used by tetClipQuadrature but with the vertex
/// order REVERSED (P3,P2,P1,P0) and one extra Gauss point in each direction
/// -- an independent code path from the clipped rule, used by T2/T3 to
/// validate it, and exact to the same volDeg = 6p as the clipped rule.
inline real_t tetVolumeOracleMoment(const TetMesh & M, int a, int b, int c, index_t p)
{
    const index_t volDeg = volDegree(p);
    const Rule1D Ru = gauss01(ceilHalf(volDeg+3)+1);
    const Rule1D Rv = gauss01(ceilHalf(volDeg+2)+1);
    const Rule1D Rw = gauss01(ceilHalf(volDeg+1)+1);

    KahanSum total;
    for (const std::array<index_t,4> & tt : M.tet)
    {
        const Vec3 & R0 = M.P[tt[3]], & R1 = M.P[tt[2]], & R2 = M.P[tt[1]], & R3 = M.P[tt[0]];
        const real_t absD = math::abs(det3(sub3(R1,R0), sub3(R2,R0), sub3(R3,R0)));
        const Vec3 e1 = sub3(R1,R0), e2 = sub3(R2,R1), e3 = sub3(R3,R2);
        for (size_t iu = 0; iu != Ru.u.size(); ++iu)
        {
            const real_t u = Ru.u[iu];
            const Vec3 Pu = add3(R0, scale3(u,e1));
            for (size_t iv = 0; iv != Rv.u.size(); ++iv)
            {
                const real_t v = Rv.u[iv];
                const Vec3 Puv = add3(Pu, scale3(u*v,e2));
                const real_t wuv = Ru.w[iu]*Rv.w[iv]*u*u*v*absD;
                for (size_t iw = 0; iw != Rw.u.size(); ++iw)
                {
                    const real_t w = Rw.u[iw];
                    const Vec3 Pt = add3(Puv, scale3(u*v*w,e3));
                    total.add(wuv*Rw.w[iw]*std::pow(Pt[0],a)*std::pow(Pt[1],b)*std::pow(Pt[2],c));
                }
            }
        }
    }
    return total.value();
}

/// Boundary moment oint x^a y^b z^c n_comp dS (comp<0: no normal factor)
/// over the WHOLE (unclipped) boundary mesh, by the collapsed-triangle map
/// with the vertex order REVERSED and one extra Gauss point in each
/// direction, exact to the same bdrDeg = 6p as the clipped rule. \a normal
/// (the correctly outward-oriented one) is used for the n_comp factor
/// regardless of the reversed traversal, which only changes how the
/// surface points themselves are generated.
inline real_t triBoundaryOracleMoment(const TetMesh & M, int a, int b, int c, int comp, index_t p)
{
    const index_t bdrDeg = bdrDegree(p);
    const Rule1D Ru = gauss01(ceilHalf(bdrDeg+2)+1);
    const Rule1D Rv = gauss01(ceilHalf(bdrDeg+1)+1);

    KahanSum total;
    for (const std::array<index_t,4> & tri : M.bdrTri)
    {
        const Vec3 & A = M.P[tri[0]], & B = M.P[tri[1]], & C = M.P[tri[2]];
        const Vec3 crAB = cross3(sub3(B,A), sub3(C,A));
        const real_t mag = math::sqrt(dot3(crAB,crAB));
        const Vec3 normal = scale3(1.0/mag, crAB);

        const Vec3 & Q0 = C, & Q1 = B, & Q2 = A;   // reversed traversal
        const Vec3 e1 = sub3(Q1,Q0), e2 = sub3(Q2,Q1);
        for (size_t iu = 0; iu != Ru.u.size(); ++iu)
        {
            const real_t u = Ru.u[iu];
            const Vec3 Pu = add3(Q0, scale3(u,e1));
            for (size_t iv = 0; iv != Rv.u.size(); ++iv)
            {
                const real_t v = Rv.u[iv];
                const Vec3 Pt = add3(Pu, scale3(u*v,e2));
                real_t m = std::pow(Pt[0],a)*std::pow(Pt[1],b)*std::pow(Pt[2],c);
                if (comp >= 0) m *= normal[comp];
                total.add(Ru.w[iu]*Rv.w[iv]*u*mag*m);
            }
        }
    }
    return total.value();
}

/// Clipped-quadrature volume moment int x^a y^b z^c dV, summed over every
/// cell's CellRule3.
inline real_t volMoment(const TetClipQuadrature & Q, int a, int b, int c)
{
    KahanSum total;
    for (const CellRule3 & r : Q.vol)
        for (index_t k = 0; k != r.weights.size(); ++k)
            total.add(r.weights[k]*std::pow(r.nodes(0,k),a)*std::pow(r.nodes(1,k),b)*std::pow(r.nodes(2,k),c));
    return total.value();
}

/// Clipped-quadrature boundary moment oint x^a y^b z^c n_comp dS, summed
/// over every cell's CellBdrRule3 (comp < 0: no normal factor, else times
/// the normal's \a comp component).
inline real_t bdrMoment(const TetClipQuadrature & Q, int a, int b, int c, int comp)
{
    KahanSum total;
    for (const CellBdrRule3 & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
        {
            real_t m = std::pow(r.nodes(0,k),a)*std::pow(r.nodes(1,k),b)*std::pow(r.nodes(2,k),c);
            if (comp >= 0) m *= r.normals(comp,k);
            total.add(r.weights[k]*m);
        }
    return total.value();
}

/// Clipped-quadrature flux oint x.n dS, the divergence-theorem check
/// against 3*volume (T8).
inline real_t bdrFluxXdotN(const TetClipQuadrature & Q)
{
    KahanSum total;
    for (const CellBdrRule3 & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
            total.add(r.weights[k]*(r.nodes(0,k)*r.normals(0,k) + r.nodes(1,k)*r.normals(1,k)
                                    + r.nodes(2,k)*r.normals(2,k)));
    return total.value();
}

} // namespace gsTetClip
