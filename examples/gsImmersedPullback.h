/** @file gsImmersedPullback.h

    @brief Curved single-patch background maps G: [-1,1]^3 -> R^3 for the P1
    tet-mesh clip pipeline (gsTetMeshClip.h / gsImmersedLookupRule.h), and the
    machinery that lets that pipeline run against one: pull-back of a
    physical tet mesh to parameter space, a boundary-rule decorator that
    turns a parametric quadrature rule into a physical one via Nanson's
    formula, and a self-test.

    Design. Every rule in gsImmersedLookupRule.h (ClipStreamer's volRule()/
    bdrRule(), and everything built on top) clips and integrates a TetMesh
    directly, in whatever coordinates that mesh happens to carry -- it has no
    notion of a background geometry map. So running the SAME machinery on a
    curved G needs two things it does not provide on its own:
    - a tet mesh whose vertices are already G^{-1} of the physical ones, so
      that ClipStreamer's clip-and-integrate is performed in PARAMETER space
      against the SAME axis-aligned grid [-1,1]^3 it already assumes
      (pullBack(), producing a ParamTetMesh -- a distinct type from
      PhysTetMesh so the two coordinate systems cannot be silently mixed at
      a call site; see the two structs below);
    - a way to turn the resulting PARAMETRIC boundary rule (nodes, weights,
      unit normals, all still living in [-1,1]^3) into a PHYSICAL one, since
      a Nitsche or flux integral needs physical measure and physical normals
      (PullbackBdrSource, Nanson's formula below).

    Two background maps are provided beyond the pre-existing affine identity
    (gsTetClipSignDomain.h's identityBoxGeometry): affineBoxGeometry (an
    arbitrary, still affine, map -- exercises pullBack's closed-form initial
    guess and PullbackBdrSource's Nanson transform under a CONSTANT Jacobian)
    and bubbleBoxGeometry (a genuinely curved map, built in closed form like
    identityBoxGeometry -- identity control points plus one moved centre
    coefficient, no linear solve -- exercises both under a Jacobian that
    VARIES over the box).

    Nanson's formula. For a diffeomorphism G with Jacobian J = dG, a physical
    outward area element transforms from a parametric one as
        n dA = det(J) J^{-T} N dA0.
    Writing J's columns as a,b,c and cof(J) = det(J) J^{-T} for its
    cofactor matrix, cof(J) N = N0 (b x c) + N1 (c x a) + N2 (a x b) for
    N = (N0,N1,N2) -- this is the adjugate/cofactor identity for a 3x3
    matrix, and PullbackBdrSource::bdrRule uses exactly this form so it never
    has to form J^{-1} (a determinant-only, cross-product-only computation
    that stays well-conditioned even where J is close to singular in one
    direction but not all three).

    Newton pull-back and its two traps (see pullBack() below for the full
    per-vertex algorithm): gsFunction::newtonRaphson's return code is NOT a
    reliable "did this converge to a true root" signal once withSupport
    clamping is in play -- it reports success whenever the clamped update
    step shrinks below the accuracy target, including when that step landed
    on the box boundary because the true root lies outside [-1,1]^3
    entirely. pullBack() therefore treats the return code only as a coarse
    nFailed filter and computes an INDEPENDENT residual -- G evaluated at the
    returned parameter value, compared against the original physical point
    -- as the real correctness gate every acceptance check in this file's
    self-test relies on.

    Bernstein/Weyl facts used by bubbleBoxGeometry (proved in its own
    doxygen below): on [-1,1], 1-x^2 = 2 B_1(t) with t=(x+1)/2 the degree-2
    Bernstein parametrisation, so the degree-(2,2,2) bubble
    b(x,y,z) = (1-x^2)(1-y^2)(1-z^2) is 8 times the SINGLE degree-2 tensor
    Bernstein basis function centred at the box's own centre -- which is
    also the middle Greville anchor of the corresponding B-spline basis. So
    G = identity + eps*b*that-one-control-point-only is exact and needs no
    linear solve. Its Jacobian J = I + eps*chat*grad(b)^T is a RANK-ONE
    perturbation of the identity, so Weyl's inequality gives
    sigma_min(J) >= 1 - |eps|*max|grad(b)| directly, with max|grad(b)| = 2
    attained at face centres -- the bound the self-test's sampled
    minSingularValueJ check is measured against.

    Complexity: pullBack() is O(nVerts * Newton iterations) for the
    per-vertex inversion (OpenMP-parallel, dynamic schedule -- each iteration
    touches only its own output slot) plus O(nTets) for the serial
    orientation/volume tally; PullbackBdrSource::bdrRule(id,...) is O(m) in
    the number of nodes its wrapped parametric rule returns for that cell;
    minSingularValueJ is O(nPerDir^3) geometry evaluations plus one 3x3 SVD
    each.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include "gsImmersedLookupRule.h"
#include "gsTetClipSignDomain.h"

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

namespace gsTetClip
{
using namespace gismo;

//----------------------------------------------------------------------------
// Physical / parametric mesh types and their Phys-only oracles.
//----------------------------------------------------------------------------

/// Which closed-form background map a ParamTetMesh (and the geometry object
/// passed alongside it) was built from. Identity needs no Newton at all
/// (pullBack() special-cases it to a bitwise copy); Affine gives pullBack()
/// an exact closed-form initial guess; Bubble falls back to the physical
/// point itself as the initial guess (the map is close enough to the
/// identity on its intended eps range that this converges).
enum class BgMapKind { Identity, Affine, Bubble };

/// A tet mesh in PHYSICAL coordinates, exactly as read by readMsh41(). Wraps
/// TetMesh only to give it a type distinct from ParamTetMesh, so a call site
/// cannot pass a parametric mesh where a physical one (or vice versa) is
/// expected without a compile error -- see the deleted ParamTetMesh
/// overloads of meshVolumeExact/unclippedBoundaryArea below.
struct PhysTetMesh { TetMesh m; };

/// A tet mesh in PARAMETRIC coordinates ([-1,1]^3), produced by pullBack().
/// Deliberately carries no conversion to/from TetMesh or PhysTetMesh: every
/// consumer (ClipStreamer, PullbackBdrSource) that needs the raw TetMesh
/// reads .m explicitly, so the physical/parametric distinction is visible at
/// every call site rather than erased by an implicit conversion.
struct ParamTetMesh { TetMesh m; };

/// Exact whole-mesh volume of a PHYSICAL tet mesh; forwards to
/// meshVolumeExact(const TetMesh&) (gsTetMeshClip.h).
inline real_t meshVolumeExact(const PhysTetMesh & M) { return meshVolumeExact(M.m); }
/// Exact whole-mesh boundary area of a PHYSICAL tet mesh; forwards to
/// unclippedBoundaryArea(const TetMesh&) (gsTetMeshClip.h).
inline real_t unclippedBoundaryArea(const PhysTetMesh & M) { return unclippedBoundaryArea(M.m); }
/// Deleted: these oracles integrate physical volume/area, which a
/// PARAMETRIC mesh's coordinates do not carry (pull the mesh back through
/// its background map's own Jacobian first, or query the corresponding
/// PhysTetMesh instead). Kept as an overload, not simply omitted, so misuse
/// is a compile error rather than a silent wrong-coordinate integral.
real_t meshVolumeExact(const ParamTetMesh &) = delete;
real_t unclippedBoundaryArea(const ParamTetMesh &) = delete;

//----------------------------------------------------------------------------
// Background geometry maps.
//----------------------------------------------------------------------------

/// An exactly affine map G(x) = A*x + b on the box of \a g: same
/// construction as gsTetClipSignDomain.h's identityBoxGeometry (degree 1,
/// one element per direction, control points placed directly, no linear
/// solve) but with the identity's own anchors first transformed by
/// (A, b) -- exact by the same linear-precision argument, since A*x+b is
/// itself affine.
///
/// \param A 3x3, ENSUREd to have positive determinant (an orientation-
///           reversing background map is out of scope for every consumer
///           in this file -- pullBack's orientation tally, in particular,
///           assumes G preserves tet handedness away from inversion).
/// \param b 3, ENSUREd to have size 3.
inline gsMultiPatch<real_t> affineBoxGeometry(const Grid3 & g, const gsMatrix<real_t> & A,
                                              const gsVector<real_t> & b)
{
    GISMO_ENSURE(3 == A.rows() && 3 == A.cols(),
                "affineBoxGeometry: A must be 3x3, got " << A.rows() << "x" << A.cols());
    GISMO_ENSURE(A.determinant() > 0,
                "affineBoxGeometry: A.determinant() must be positive, got " << A.determinant());
    GISMO_ENSURE(3 == b.size(), "affineBoxGeometry: b.size() must be 3, got " << b.size());

    std::vector<real_t> X, Y, Z;
    gridLines(g, X, Y, Z);
    const std::vector<real_t> Lx = { X.front(), X.back() };
    const std::vector<real_t> Ly = { Y.front(), Y.back() };
    const std::vector<real_t> Lz = { Z.front(), Z.back() };
    gsTensorBSplineBasis<3,real_t> gb(gridKnots(Lx,1), gridKnots(Ly,1), gridKnots(Lz,1));
    gsMatrix<real_t> coefs = ((A * gb.anchors()).colwise() + b).transpose();
    return gsMultiPatch<real_t>(*gb.makeGeometry(give(coefs)));
}

/// A curved map G(x) = x + eps*b(x)*chat on the box of \a g, with
/// chat = c/|c| and the degree-(2,2,2) bubble b(x,y,z) = (1-x^2)(1-y^2)(1-z^2).
///
/// Construction: closed form, not a linear solve. On [-1,1] with
/// t=(x+1)/2 in [0,1], 1-x^2 = (1-x)(1+x) = 2(1-t)*2t = 4t(1-t) = 2*B_1(t),
/// where B_1(t) = 2t(1-t) is the degree-2 Bernstein basis function centred
/// at t=1/2 -- so b = (2B_1(t_x))(2B_1(t_y))(2B_1(t_z)) =
/// 8*B_1(t_x)B_1(t_y)B_1(t_z), a tensor-Bernstein monomial that is
/// non-zero only at the SINGLE degree-2
/// tensor control point sitting at the box centre. Linear precision places
/// the identity's own control points at the basis' Greville anchors
/// {-1,0,1} per direction (27 anchors total, tensor order
/// i + 3*(j + 3*k)); the centre anchor -- index 13 = 1 + 3*(1 + 3*1) -- is
/// therefore the ONLY coefficient the bubble term touches, and it moves by
/// exactly 8*eps*chat. This reproduces the same map an
/// interpolate-at-anchors linear solve would (up to rounding), bitwise
/// deterministically and without ever forming a system matrix.
///
/// Its Jacobian is J = I + eps*chat*grad(b)^T, a RANK-ONE perturbation of
/// the identity (grad(b) is a single vector field, not a matrix). Its
/// determinant follows from the matrix determinant lemma,
/// det(J) = 1 + eps*chat.grad(b) >= 1 - 2|eps|. Its smallest singular
/// value follows from Weyl's perturbation inequality for singular values,
/// sigma_min(J) >= sigma_min(I) - ||eps*chat*grad(b)^T||_2
/// = 1 - |eps|*|grad(b)| >= 1 - 2|eps|, since the 2-norm of the rank-one
/// term is exactly |eps|*|grad(b)| (|chat| = 1) and max|grad(b)| = 2 is
/// attained at the six face centres of the box (at the face centre
/// (+-1,0,0), grad(b) = (-+2,0,0); likewise for the other faces). G fixes
/// every point of the box's own boundary pointwise: on any face, one of
/// the three factors of b vanishes identically.
///
/// \param eps ENSUREd |eps| < 0.5 (keeps the Weyl bound above 0, though the
///            true sigma_min may still be smaller for eps close to that
///            edge -- see minSingularValueJ for a direct sampled check).
/// \param c   ENSUREd |c| > 0; normalised internally to chat = c/|c|.
inline gsMultiPatch<real_t> bubbleBoxGeometry(const Grid3 & g, real_t eps, const Vec3 & c)
{
    const real_t cnorm = std::sqrt(dot3(c,c));
    GISMO_ENSURE(cnorm > 0, "bubbleBoxGeometry: |c| must be positive");
    GISMO_ENSURE(std::abs(eps) < 0.5, "bubbleBoxGeometry: |eps| must be < 0.5, got " << eps);

    std::vector<real_t> X, Y, Z;
    gridLines(g, X, Y, Z);
    GISMO_ENSURE(std::abs(X.front()+1) <= 1e-14 && std::abs(X.back()-1) <= 1e-14 &&
                std::abs(Y.front()+1) <= 1e-14 && std::abs(Y.back()-1) <= 1e-14 &&
                std::abs(Z.front()+1) <= 1e-14 && std::abs(Z.back()-1) <= 1e-14,
                "bubbleBoxGeometry: the grid box must be [-1,1]^3 to within 1e-14, got "
                "[" << X.front() << "," << X.back() << "] x [" << Y.front() << "," << Y.back()
                << "] x [" << Z.front() << "," << Z.back() << "]");

    const std::vector<real_t> Lx = { X.front(), X.back() };
    const std::vector<real_t> Ly = { Y.front(), Y.back() };
    const std::vector<real_t> Lz = { Z.front(), Z.back() };
    gsTensorBSplineBasis<3,real_t> gb(gridKnots(Lx,2), gridKnots(Ly,2), gridKnots(Lz,2));

    const gsMatrix<real_t> anchors = gb.anchors(); // 3 x 27
    GISMO_ENSURE(27 == anchors.cols(),
                "bubbleBoxGeometry: expected 27 anchors, got " << anchors.cols());
    GISMO_ENSURE(anchors.col(13).cwiseAbs().maxCoeff() <= 1e-14,
                "bubbleBoxGeometry: anchors.col(13) is not within 1e-14 of the origin: "
                << anchors.col(13).transpose());

    const Vec3 chat = scale3(1.0/cnorm, c);
    gsMatrix<real_t> coefs = anchors.transpose(); // 27 x 3
    coefs(13,0) += 8*eps*chat[0];
    coefs(13,1) += 8*eps*chat[1];
    coefs(13,2) += 8*eps*chat[2];
    return gsMultiPatch<real_t>(*gb.makeGeometry(give(coefs)));
}

/// True iff \a G's control points equal its own basis' Greville anchors
/// (transposed to N x d) to within \a tol -- i.e. G reproduces the
/// identity function exactly, the same linear-precision fact
/// identityBoxGeometry/affineBoxGeometry construct FROM. Returns false
/// immediately on any shape mismatch (never throws on that). \a tol = 0
/// is exact equality, which is what identityBoxGeometry itself produces.
inline bool isIdentityMap(const gsGeometry<real_t> & G, real_t tol = 0)
{
    const gsMatrix<real_t> anchorsT = G.basis().anchors().transpose();
    if (G.coefs().rows() != anchorsT.rows() || G.coefs().cols() != anchorsT.cols())
        return false;
    return (G.coefs() - anchorsT).cwiseAbs().maxCoeff() <= tol;
}

//----------------------------------------------------------------------------
// Pull-back of a physical tet mesh to parameter space.
//----------------------------------------------------------------------------

/// Diagnostics of one pullBack() call. Every count is over the mesh's
/// nVerts vertices except nInverted, which is over its tets.
struct PullbackStats { index_t nVerts=0, nFailed=0, nOutside=0, nInverted=0; real_t maxResidual=0, minParamTetVol=0; };

/// Pulls \a M back through the background map \a G: returns a ParamTetMesh
/// whose connectivity (tet, bdrTri, bdrOwner, nNodesFile) is copied
/// unchanged from \a M.m and whose vertices are G^{-1}(x) for every
/// physical vertex x, so that gsImmersedLookupRule.h's clip machinery can
/// run against it directly on the SAME [-1,1]^3 grid \a g uses elsewhere.
///
/// Identity is special-cased to a bitwise copy (no Newton, no rounding):
/// \a st is filled as every other case would (residual 0, the orientation
/// pass below still runs, since it is O(nTets) and catches a caller error
/// in \a M as readily as in a curved case).
///
/// Otherwise, per vertex i with physical position x_i:
/// 1. an initial guess u0 -- for Bubble, x_i itself (the map stays close to
///    the identity over its intended |eps| range); for Affine,
///    u0 = c + J^{-1}(x_i - G(c)) with c the centre of G's own support and
///    J its Jacobian there (both computed ONCE, before the per-vertex loop
///    -- exact for a genuinely affine G, up to rounding, since J is then
///    constant);
/// 2. u0 clamped into G.support() (gsFunction::newtonRaphson's withSupport
///    path GISMO_ASSERTs the INITIAL point already lies in the support);
/// 3. G.newtonRaphson(x_i, u, true, 1e-13, 250), run under
///    `#pragma omp parallel for schedule(dynamic)` over vertices -- each
///    iteration owns only its own P[i] and code[i], so there is no race;
/// 4. an INDEPENDENT residual r_i = |G(u_i) - x_i|, from one batched
///    G.eval_into() over every new vertex at once, because
///    newtonRaphson's return code alone cannot be trusted (see the file
///    header) -- it reports success whenever the CLAMPED update step is
///    small, even when that means the point converged onto the box
///    boundary because its true pre-image lies outside [-1,1]^3 entirely;
/// 5. a serial tally into \a st: nFailed counts a -1 return code or any
///    non-finite parameter component; maxResidual is the max r_i (+inf if
///    any r_i is non-finite); nOutside counts vertices not STRICTLY inside
///    the open box (matching the driver's own makeGrid strict-interior
///    convention); nInverted counts tets whose parametric signed volume
///    D_param has flipped sign against the physical D_phys, or is exactly
///    zero (readMsh41 does not orient tets, so the comparison is always
///    against each tet's OWN physical sign, never a fixed convention);
///    minParamTetVol is min over tets of sign(D_phys)*D_param/6, the same
///    orientation normalisation, so a positive value here means every tet
///    kept a consistent, non-degenerate orientation under the pull-back.
///
/// nInverted == 0 also guarantees the copied bdrTri stay outward in
/// parameter space: readMsh41's own outward orientation test (the boundary
/// derivation block of readMsh41, gsTetMeshClip.h) is itself a signed tet
/// volume built from the SAME owner-tet vertices this tally already checks,
/// so a tet that kept its orientation keeps its boundary faces outward too.
///
/// Never throws on bad stats -- callers gate on \a st themselves.
/// Complexity: O(nVerts * Newton iterations) for the parallel loop plus
/// O(nTets) for the serial tally.
inline ParamTetMesh pullBack(const PhysTetMesh & M, const gsGeometry<real_t> & G, BgMapKind kind,
                             const Grid3 & g, PullbackStats & st)
{
    GISMO_ENSURE(3 == G.targetDim() && 3 == G.domainDim(),
                "pullBack: G must be a 3->3 map, got domainDim=" << G.domainDim()
                << " targetDim=" << G.targetDim());

    st = PullbackStats();
    ParamTetMesh out{ M.m };
    const index_t nP = (index_t)out.m.P.size();
    st.nVerts = nP;

    std::vector<real_t> residual((size_t)nP, 0.0);
    std::vector<int> code((size_t)nP, 0);

    if (BgMapKind::Identity != kind)
    {
        const gsMatrix<real_t> supp = G.support(); // 3x2: [lo hi]

        gsVector<real_t> c(3), Gc(3);
        gsMatrix<real_t,3,3> Jinv;
        if (BgMapKind::Affine == kind)
        {
            c = 0.5*(supp.col(0) + supp.col(1));
            gsFuncData<real_t> fdc(NEED_VALUE | NEED_DERIV);
            G.compute(c, fdc);
            const gsMatrix<real_t,3,3> Jc = fdc.jacobian(0);
            Gc = fdc.values[0].col(0);
            Jinv = Jc.inverse();
        }

        gsMatrix<real_t> Pnew(3, nP);
        const long nPl = (long)nP;
        #pragma omp parallel for schedule(dynamic)
        for (long ii = 0; ii < nPl; ++ii)
        {
            const size_t i = (size_t)ii;
            const Vec3 & x = M.m.P[i];
            gsVector<real_t> xVec(3); xVec << x[0], x[1], x[2];

            gsVector<real_t> u(3);
            if (BgMapKind::Affine == kind) u = c + Jinv*(xVec - Gc);
            else                           u = xVec;
            u = u.cwiseMax(supp.col(0)).cwiseMin(supp.col(1));

            code[i] = G.newtonRaphson(xVec, u, true, (real_t)1e-13, 250);
            Pnew.col(ii) = u;
            out.m.P[i] = { u[0], u[1], u[2] };
        }

        gsMatrix<real_t> vals(3, nP);
        G.eval_into(Pnew, vals);
        for (index_t i = 0; i < nP; ++i)
        {
            const Vec3 & x = M.m.P[(size_t)i];
            const real_t dx = vals(0,i)-x[0], dy = vals(1,i)-x[1], dz = vals(2,i)-x[2];
            residual[(size_t)i] = std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    }

    std::vector<real_t> X, Y, Z;
    gridLines(g, X, Y, Z);

    st.nFailed = 0;
    st.maxResidual = 0;
    st.nOutside = 0;
    for (index_t i = 0; i < nP; ++i)
    {
        const Vec3 & u = out.m.P[(size_t)i];
        const bool finiteU = std::isfinite(u[0]) && std::isfinite(u[1]) && std::isfinite(u[2]);
        if (-1 == code[(size_t)i] || !finiteU) ++st.nFailed;

        const real_t r = residual[(size_t)i];
        if (!std::isfinite(r)) st.maxResidual = std::numeric_limits<real_t>::infinity();
        else st.maxResidual = math::max(st.maxResidual, r);

        if (!(u[0] > X.front() && u[0] < X.back() &&
              u[1] > Y.front() && u[1] < Y.back() &&
              u[2] > Z.front() && u[2] < Z.back()))
            ++st.nOutside;
    }

    st.nInverted = 0;
    st.minParamTetVol = std::numeric_limits<real_t>::infinity();
    for (const std::array<index_t,4> & tt : M.m.tet)
    {
        const Vec3 & Px0=M.m.P[tt[0]], & Px1=M.m.P[tt[1]], & Px2=M.m.P[tt[2]], & Px3=M.m.P[tt[3]];
        const real_t Dphys = det3(sub3(Px1,Px0), sub3(Px2,Px0), sub3(Px3,Px0));

        const Vec3 & Pu0=out.m.P[tt[0]], & Pu1=out.m.P[tt[1]], & Pu2=out.m.P[tt[2]], & Pu3=out.m.P[tt[3]];
        const real_t Dparam = det3(sub3(Pu1,Pu0), sub3(Pu2,Pu0), sub3(Pu3,Pu0));

        if (0.0 == Dparam || (Dparam > 0) != (Dphys > 0)) ++st.nInverted;
        const real_t signPhys = (Dphys > 0) ? (real_t)1 : (real_t)-1;
        st.minParamTetVol = math::min(st.minParamTetVol, signPhys*Dparam/6.0);
    }

    for (int d = 0; d != 3; ++d)
    { out.m.lo[d] = std::numeric_limits<real_t>::infinity(); out.m.hi[d] = -out.m.lo[d]; }
    for (const std::array<index_t,4> & tt : out.m.tet)
        for (int v = 0; v != 4; ++v)
            for (int d = 0; d != 3; ++d)
            {
                out.m.lo[d] = math::min(out.m.lo[d], out.m.P[tt[v]][d]);
                out.m.hi[d] = math::max(out.m.hi[d], out.m.P[tt[v]][d]);
            }

    return out;
}

//----------------------------------------------------------------------------
// Nanson boundary decorator.
//----------------------------------------------------------------------------

/// Wraps a PARAMETRIC BdrCellSource (typically a ClipStreamer built on a
/// ParamTetMesh) and turns its output into a PHYSICAL one under the
/// background map \a G, via Nanson's formula n dA = det(J) J^{-T} N dA0
/// (see the file header). Nodes are left UNCHANGED -- they stay parametric,
/// because BdrNormalField (gsImmersedLookupRule.h) matches normals back to
/// quadrature nodes bitwise, and only weights/normals carry the physical
/// transform here.
///
/// const and thread-safe: every call reads only \a paramSrc and \a G, both
/// stored by reference and both read-only from here, so concurrent calls
/// (e.g. from an OpenMP loop over cell id) never race. Both referenced
/// objects must outlive this one.
///
/// Uses the cofactor form cof(J)*that-normal = det(J)*J^{-T}*that-normal
/// rather than forming J^{-1} directly: cof(J)*N for J's columns (a,b,c)
/// and N=(N0,N1,N2) is N0*(b x c) + N1*(c x a) + N2*(a x b), so the whole
/// transform is a handful of cross products and one dot product, with no
/// linear solve and no risk of amplifying an ill-conditioned inverse.
class PullbackBdrSource : public BdrCellSource
{
public:
    /// \a paramSrc and \a G must outlive this object; neither is copied.
    PullbackBdrSource(const BdrCellSource & paramSrc, const gsGeometry<real_t> & G)
    : m_paramSrc(paramSrc), m_G(G) {}

    void bdrRule(size_t id, gsMatrix<real_t> & nodes, gsVector<real_t> & weights,
                gsMatrix<real_t> & normals) const override
    {
        m_paramSrc.bdrRule(id, nodes, weights, normals);
        const index_t m = nodes.cols();
        if (0 == m) return;

        gsFuncData<real_t> fd(NEED_DERIV);
        m_G.compute(nodes, fd);

        for (index_t k = 0; k != m; ++k)
        {
            const gsMatrix<real_t,3,3> J = fd.jacobian(k);
            const Vec3 a = { J(0,0), J(1,0), J(2,0) };
            const Vec3 b = { J(0,1), J(1,1), J(2,1) };
            const Vec3 c = { J(0,2), J(1,2), J(2,2) };
            const Vec3 nhat = { normals(0,k), normals(1,k), normals(2,k) };

            const Vec3 mvec = add3(add3(scale3(nhat[0], cross3(b,c)),
                                        scale3(nhat[1], cross3(c,a))),
                                   scale3(nhat[2], cross3(a,b)));
            const real_t mmag = std::sqrt(dot3(mvec,mvec));

            weights[k] *= mmag;

            const real_t detJ = det3(a,b,c);
            const real_t s = (detJ >= 0.0) ? (real_t)1 : (real_t)-1;
            normals(0,k) = s*mvec[0]/mmag;
            normals(1,k) = s*mvec[1]/mmag;
            normals(2,k) = s*mvec[2]/mmag;
        }
    }

private:
    const BdrCellSource & m_paramSrc;
    const gsGeometry<real_t> & m_G;
};

//----------------------------------------------------------------------------
// Sampled singular-value bound.
//----------------------------------------------------------------------------

/// Minimum, over an \a nPerDir^3 tensor grid on G.support() (including its
/// corners), of the smallest singular value of G's Jacobian. Since it is a
/// SAMPLED minimum, it can only OVER-estimate the true minimum over the
/// whole box -- it is a necessary, not sufficient, non-degeneracy check.
inline real_t minSingularValueJ(const gsGeometry<real_t> & G, index_t nPerDir)
{
    GISMO_ENSURE(nPerDir >= 2, "minSingularValueJ: nPerDir must be >= 2, got " << nPerDir);
    const gsMatrix<real_t> supp = G.support(); // 3x2
    const index_t N = nPerDir*nPerDir*nPerDir;

    gsMatrix<real_t> pts(3, N);
    index_t col = 0;
    for (index_t k = 0; k != nPerDir; ++k)
    for (index_t j = 0; j != nPerDir; ++j)
    for (index_t i = 0; i != nPerDir; ++i)
    {
        pts(0,col) = supp(0,0) + (supp(0,1)-supp(0,0))*i/(real_t)(nPerDir-1);
        pts(1,col) = supp(1,0) + (supp(1,1)-supp(1,0))*j/(real_t)(nPerDir-1);
        pts(2,col) = supp(2,0) + (supp(2,1)-supp(2,0))*k/(real_t)(nPerDir-1);
        ++col;
    }

    gsFuncData<real_t> fd(NEED_DERIV);
    G.compute(pts, fd);

    real_t sMin = std::numeric_limits<real_t>::infinity();
    for (index_t k = 0; k != N; ++k)
    {
        const gsMatrix<real_t,3,3> J = fd.jacobian(k);
        gsMatrix<real_t>::JacobiSVD svd(J);
        sMin = math::min(sMin, svd.singularValues().minCoeff());
    }
    return sMin;
}

//----------------------------------------------------------------------------
// Test-map constants and self-test.
//----------------------------------------------------------------------------

/// The affine test map of the P2 ladder: A = R_z(20 deg) * R_x(15 deg) * S,
/// S = [[1,0.15,0],[0,1,0.1],[0,0,1]], b = (0.02,-0.01,0.03), with
/// R_z(theta) = [[cos,-sin,0],[sin,cos,0],[0,0,1]] and
/// R_x(phi) = [[1,0,0],[0,cos,-sin],[0,sin,cos]]. A.determinant() > 0
/// (a rotation composed with an upper-triangular unit-diagonal shear), so
/// it satisfies affineBoxGeometry's own ENSURE.
inline void defaultAffineMap(gsMatrix<real_t> & A, gsVector<real_t> & b)
{
    const real_t thetaZ = 20.0*std::acos(-1.0)/180.0;
    const real_t phiX   = 15.0*std::acos(-1.0)/180.0;

    gsMatrix<real_t> Rz(3,3), Rx(3,3), S(3,3);
    Rz << std::cos(thetaZ), -std::sin(thetaZ), 0,
          std::sin(thetaZ),  std::cos(thetaZ), 0,
          0,                 0,                1;
    Rx << 1, 0,               0,
          0, std::cos(phiX), -std::sin(phiX),
          0, std::sin(phiX),  std::cos(phiX);
    S  << 1,    0.15, 0,
          0,    1,    0.1,
          0,    0,    1;

    A = Rz*Rx*S;
    b.resize(3);
    b << 0.02, -0.01, 0.03;
}

/// Default bubble direction c = (1, 0.7, -0.5), NOT normalised --
/// bubbleBoxGeometry normalises it to chat internally.
inline Vec3 defaultBubbleDirection() { return { 1.0, 0.7, -0.5 }; }

/// Self-check of every entity in this header, over both shipped test
/// meshes (sphere, rotcube). Prints one
/// `PULLBACK-CHECK <case> <map> n=<n> <quantity> value=<v> tol=<t> PASS|FAIL`
/// line per check and one
/// `PULLBACK-TIME <case> <block> seconds=<t>` line per timed block
/// (identity, affine, bubble, bubble-divergence) -- the tet-mesh Poisson
/// driver's (`immersed_tetmesh_poisson_example`) `--study check` calls
/// this on every run, so its cost must stay visible. Returns
/// false (after printing a FAIL line, never throwing) if either mesh file
/// cannot be resolved by gsFileManager::find, and false iff any check
/// failed.
inline bool pullbackSelfTest()
{
    bool allPass = true;

    auto reportCheck = [&](const std::string & caseName, const std::string & mapName, index_t n,
                           const std::string & quantity, real_t value, real_t tol, bool pass)
    {
        gsInfo << "PULLBACK-CHECK " << caseName << " " << mapName << " n=" << n << " "
               << quantity << " value=" << fmtSci(value,6) << " tol=" << fmtSci(tol,3)
               << " " << (pass ? "PASS" : "FAIL") << "\n";
        if (!pass) allPass = false;
    };

    auto makeG = [](index_t n)
    {
        Grid3 g; g.x0 = g.y0 = g.z0 = -1.0; g.h = 2.0/(real_t)n; g.n = n; return g;
    };
    auto strictlyInside = [](const TetMesh & m, const Grid3 & g)
    {
        return m.lo[0] > g.x0 && m.hi[0] < g.x0 + g.n*g.h &&
               m.lo[1] > g.y0 && m.hi[1] < g.y0 + g.n*g.h &&
               m.lo[2] > g.z0 && m.hi[2] < g.z0 + g.n*g.h;
    };

    gsMatrix<real_t> A; gsVector<real_t> b;
    defaultAffineMap(A, b);
    const Vec3 cdir = defaultBubbleDirection();
    const real_t eps = 0.3;

    const char * cases[2] = { "sphere", "rotcube" };
    const char * files[2] = { "volumes/tetmesh_sphere.msh", "volumes/tetmesh_cube_rotated.msh" };

    for (int ci = 0; ci != 2; ++ci)
    {
        const std::string caseName = cases[ci];
        const std::string resolved = gsFileManager::find(files[ci]);
        if (resolved.empty())
        {
            reportCheck(caseName, "-", 0, "mesh-file-found", 0.0, 0.0, false);
            continue;
        }
        const PhysTetMesh phys{ readMsh41(resolved) };

        real_t maxAbsX = 1.0;
        for (const Vec3 & x : phys.m.P)
            maxAbsX = math::max(maxAbsX, math::max(std::abs(x[0]), math::max(std::abs(x[1]), std::abs(x[2]))));

        //--------------------------------------------------------------
        // Block: identity, n = 4.
        //--------------------------------------------------------------
        gsStopwatch swId;
        {
            const index_t n = 4;
            const Grid3 g = makeG(n);
            const gsMultiPatch<real_t> mp = identityBoxGeometry(g);
            const gsGeometry<real_t> & G = mp.patch(0);

            PullbackStats st;
            const ParamTetMesh param = pullBack(phys, G, BgMapKind::Identity, g, st);

            const bool bitwiseP = (phys.m.P.size() == param.m.P.size()) &&
                0 == std::memcmp(phys.m.P.data(), param.m.P.data(), phys.m.P.size()*sizeof(Vec3));
            reportCheck(caseName, "identity", n, "P-bitwise-equal", bitwiseP?1.0:0.0, 0.0, bitwiseP);

            const bool tetEq = (phys.m.tet == param.m.tet);
            reportCheck(caseName, "identity", n, "tet-equal", tetEq?1.0:0.0, 0.0, tetEq);
            const bool bdrTriEq = (phys.m.bdrTri == param.m.bdrTri);
            reportCheck(caseName, "identity", n, "bdrTri-equal", bdrTriEq?1.0:0.0, 0.0, bdrTriEq);
            const bool bdrOwnerEq = (phys.m.bdrOwner == param.m.bdrOwner);
            reportCheck(caseName, "identity", n, "bdrOwner-equal", bdrOwnerEq?1.0:0.0, 0.0, bdrOwnerEq);

            reportCheck(caseName, "identity", n, "maxResidual", st.maxResidual, 0.0, 0.0==st.maxResidual);
            reportCheck(caseName, "identity", n, "nFailed", (real_t)st.nFailed, 0.0, 0==st.nFailed);
            reportCheck(caseName, "identity", n, "nOutside", (real_t)st.nOutside, 0.0, 0==st.nOutside);
            reportCheck(caseName, "identity", n, "nInverted", (real_t)st.nInverted, 0.0, 0==st.nInverted);

            const bool idMap = isIdentityMap(G);
            reportCheck(caseName, "identity", n, "isIdentityMap", idMap?1.0:0.0, 0.0, idMap);
            const real_t sMinId = minSingularValueJ(G, 5);
            reportCheck(caseName, "identity", n, "minSingularValueJ", sMinId, 1e-15,
                       std::abs(sMinId-1.0) <= 1e-15);
        }
        gsInfo << "PULLBACK-TIME " << caseName << " identity seconds=" << fmtSci(swId.stop(),3) << "\n";

        //--------------------------------------------------------------
        // Block: affine, n = 4 (plus Nanson at n = 4 and n = 8).
        //--------------------------------------------------------------
        gsStopwatch swAff;
        {
            const index_t n = 4;
            const Grid3 g = makeG(n);
            const gsMultiPatch<real_t> mp = affineBoxGeometry(g, A, b);
            const gsGeometry<real_t> & G = mp.patch(0);

            PullbackStats st;
            const ParamTetMesh param = pullBack(phys, G, BgMapKind::Affine, g, st);

            reportCheck(caseName, "affine", n, "nFailed", (real_t)st.nFailed, 0.0, 0==st.nFailed);
            reportCheck(caseName, "affine", n, "nOutside", (real_t)st.nOutside, 0.0, 0==st.nOutside);
            reportCheck(caseName, "affine", n, "nInverted", (real_t)st.nInverted, 0.0, 0==st.nInverted);
            const real_t resTol = 1e-12*maxAbsX;
            reportCheck(caseName, "affine", n, "maxResidual", st.maxResidual, resTol, st.maxResidual<=resTol);

            const bool idMap = isIdentityMap(G);
            reportCheck(caseName, "affine", n, "isIdentityMap", idMap?1.0:0.0, 0.0, !idMap);

            const real_t sMin = minSingularValueJ(G, 5);
            gsMatrix<real_t>::JacobiSVD svdA(A);
            const real_t sMinA = svdA.singularValues().minCoeff();
            const real_t diffSv = std::abs(sMin - sMinA);
            reportCheck(caseName, "affine", n, "minSingularValueJ-diff", diffSv, 1e-14, diffSv<=1e-14);

            for (index_t nn : { (index_t)4, (index_t)8 })
            {
                const Grid3 gN = makeG(nn);
                const gsMultiPatch<real_t> mpN = affineBoxGeometry(gN, A, b);
                const gsGeometry<real_t> & GN = mpN.patch(0);

                PullbackStats stN;
                const ParamTetMesh paramN = pullBack(phys, GN, BgMapKind::Affine, gN, stN);

                if (!strictlyInside(paramN.m, gN))
                {
                    reportCheck(caseName, "affine-nanson", nn, "param-bbox-inside-box", 0.0, 0.0, false);
                    continue;
                }

                const memory::shared_ptr<TetMesh> paramMesh = memory::make_shared(new TetMesh(paramN.m));
                ClipStreamer streamer(paramMesh, gN, 1);
                PullbackBdrSource nanson(streamer, GN);

                KahanSum areaSum, fluxSum;
                real_t maxNormDev = 0;
                const size_t nCells = (size_t)gN.n*(size_t)gN.n*(size_t)gN.n;
                for (size_t id = 0; id != nCells; ++id)
                {
                    gsMatrix<real_t> nodes, normals; gsVector<real_t> weights;
                    nanson.bdrRule(id, nodes, weights, normals);
                    const index_t m = nodes.cols();
                    if (0 == m) continue;
                    gsMatrix<real_t> gx(3, m);
                    GN.eval_into(nodes, gx);
                    for (index_t kk = 0; kk != m; ++kk)
                    {
                        areaSum.add(weights[kk]);
                        const real_t flux = gx(0,kk)*normals(0,kk) + gx(1,kk)*normals(1,kk) + gx(2,kk)*normals(2,kk);
                        fluxSum.add(weights[kk]*flux);
                        const real_t nnorm = std::sqrt(normals(0,kk)*normals(0,kk) +
                                                       normals(1,kk)*normals(1,kk) + normals(2,kk)*normals(2,kk));
                        maxNormDev = math::max(maxNormDev, std::abs(nnorm-1.0));
                    }
                }

                const real_t areaExact = unclippedBoundaryArea(phys);
                const real_t areaErr = scaledErr(areaSum.value(), areaExact, areaExact);
                reportCheck(caseName, "affine-nanson", nn, "area-relerr", areaErr, 1e-13, areaErr<=1e-13);

                const real_t volExact3 = 3.0*meshVolumeExact(phys);
                const real_t fluxErr = scaledErr(fluxSum.value(), volExact3, std::abs(volExact3));
                reportCheck(caseName, "affine-nanson", nn, "flux-relerr", fluxErr, 1e-13, fluxErr<=1e-13);

                reportCheck(caseName, "affine-nanson", nn, "normals-unit-maxdev", maxNormDev, 1e-14,
                           maxNormDev <= 1e-14);
            }
        }
        gsInfo << "PULLBACK-TIME " << caseName << " affine seconds=" << fmtSci(swAff.stop(),3) << "\n";

        //--------------------------------------------------------------
        // Block: bubble (eps = 0.3, defaultBubbleDirection()), n = 4.
        //--------------------------------------------------------------
        gsStopwatch swBub;
        gsMultiPatch<real_t> mpBub;
        {
            const Grid3 g4 = makeG(4);
            mpBub = bubbleBoxGeometry(g4, eps, cdir);
            const gsGeometry<real_t> & G = mpBub.patch(0);

            const real_t cnorm = std::sqrt(dot3(cdir,cdir));
            const Vec3 chat = scale3(1.0/cnorm, cdir);
            const index_t nS = 11;
            gsMatrix<real_t> pts(3, nS*nS*nS);
            index_t col = 0;
            for (index_t k = 0; k != nS; ++k)
            for (index_t j = 0; j != nS; ++j)
            for (index_t i = 0; i != nS; ++i)
            {
                pts(0,col) = -1.0 + 2.0*i/(real_t)(nS-1);
                pts(1,col) = -1.0 + 2.0*j/(real_t)(nS-1);
                pts(2,col) = -1.0 + 2.0*k/(real_t)(nS-1);
                ++col;
            }
            gsMatrix<real_t> vals(3, pts.cols());
            G.eval_into(pts, vals);

            real_t maxAbsErr = 0.0;
            for (index_t c = 0; c != pts.cols(); ++c)
            {
                const real_t x=pts(0,c), y=pts(1,c), z=pts(2,c);
                const real_t bub = (1.0-x*x)*(1.0-y*y)*(1.0-z*z);
                maxAbsErr = math::max(maxAbsErr, std::abs(vals(0,c) - (x + eps*bub*chat[0])));
                maxAbsErr = math::max(maxAbsErr, std::abs(vals(1,c) - (y + eps*bub*chat[1])));
                maxAbsErr = math::max(maxAbsErr, std::abs(vals(2,c) - (z + eps*bub*chat[2])));
            }
            reportCheck(caseName, "bubble", nS, "analytic-maxabs", maxAbsErr, 1e-14, maxAbsErr<=1e-14);

            const bool idMap = isIdentityMap(G);
            reportCheck(caseName, "bubble", 4, "isIdentityMap", idMap?1.0:0.0, 0.0, !idMap);
            const real_t sMinBub = minSingularValueJ(G, 21);
            reportCheck(caseName, "bubble", 21, "minSingularValueJ", sMinBub, 0.4-1e-12, sMinBub >= 0.4-1e-12);

            PullbackStats stBub;
            const ParamTetMesh paramBub = pullBack(phys, G, BgMapKind::Bubble, g4, stBub);

            reportCheck(caseName, "bubble", 4, "nFailed", (real_t)stBub.nFailed, 0.0, 0==stBub.nFailed);
            reportCheck(caseName, "bubble", 4, "nOutside", (real_t)stBub.nOutside, 0.0, 0==stBub.nOutside);
            reportCheck(caseName, "bubble", 4, "nInverted", (real_t)stBub.nInverted, 0.0, 0==stBub.nInverted);

            bool allFinite = true;
            for (const Vec3 & u : paramBub.m.P)
                if (!(std::isfinite(u[0]) && std::isfinite(u[1]) && std::isfinite(u[2]))) allFinite = false;
            reportCheck(caseName, "bubble", 4, "u-finite", allFinite?1.0:0.0, 0.0, allFinite);

            const real_t resTolBub = 1e-12*maxAbsX;
            reportCheck(caseName, "bubble", 4, "maxResidual", stBub.maxResidual, resTolBub,
                       stBub.maxResidual<=resTolBub);
            reportCheck(caseName, "bubble", 4, "minParamTetVol", stBub.minParamTetVol, 0.0,
                       stBub.minParamTetVol > 0.0);
        }
        gsInfo << "PULLBACK-TIME " << caseName << " bubble seconds=" << fmtSci(swBub.stop(),3) << "\n";

        //--------------------------------------------------------------
        // Block: bubble-divergence -- Nanson with varying J, n = 4, p = 2.
        //--------------------------------------------------------------
        gsStopwatch swBd;
        {
            const index_t n = 4;
            const Grid3 g = makeG(n);
            const gsGeometry<real_t> & G = mpBub.patch(0);

            PullbackStats st2;
            const ParamTetMesh param2 = pullBack(phys, G, BgMapKind::Bubble, g, st2);

            if (!strictlyInside(param2.m, g))
            {
                reportCheck(caseName, "bubble-divergence", n, "param-bbox-inside-box", 0.0, 0.0, false);
            }
            else
            {
                const memory::shared_ptr<TetMesh> paramMesh = memory::make_shared(new TetMesh(param2.m));
                ClipStreamer streamer(paramMesh, g, 2);
                PullbackBdrSource nanson(streamer, G);

                KahanSum areaSum, fluxSum, volSum;
                const size_t nCells = (size_t)g.n*(size_t)g.n*(size_t)g.n;
                for (size_t id = 0; id != nCells; ++id)
                {
                    gsMatrix<real_t> nodes, normals; gsVector<real_t> weights;
                    nanson.bdrRule(id, nodes, weights, normals);
                    const index_t m = nodes.cols();
                    if (m > 0)
                    {
                        gsMatrix<real_t> gx(3, m);
                        G.eval_into(nodes, gx);
                        for (index_t kk = 0; kk != m; ++kk)
                        {
                            areaSum.add(weights[kk]);
                            const real_t flux = gx(0,kk)*normals(0,kk) + gx(1,kk)*normals(1,kk) + gx(2,kk)*normals(2,kk);
                            fluxSum.add(weights[kk]*flux);
                        }
                    }

                    gsMatrix<real_t> vnodes; gsVector<real_t> vweights;
                    streamer.volRule(id, vnodes, vweights);
                    const index_t mv = vnodes.cols();
                    if (mv > 0)
                    {
                        gsFuncData<real_t> fd(NEED_DERIV);
                        G.compute(vnodes, fd);
                        for (index_t kk = 0; kk != mv; ++kk)
                        {
                            const gsMatrix<real_t,3,3> J = fd.jacobian(kk);
                            volSum.add(vweights[kk]*J.determinant());
                        }
                    }
                }

                const real_t lhs = fluxSum.value();
                const real_t rhs = 3.0*volSum.value();
                const real_t relErr = scaledErr(lhs, rhs, std::abs(rhs));
                reportCheck(caseName, "bubble-divergence", n, "flux-vs-divJ-relerr", relErr, 1e-12, relErr<=1e-12);

                gsInfo << "REPORT " << caseName << " bubble-nanson-area sum=" << fmtSci(areaSum.value(),10)
                       << " unclippedBoundaryArea=" << fmtSci(unclippedBoundaryArea(phys),10) << "\n";
            }
        }
        gsInfo << "PULLBACK-TIME " << caseName << " bubble-divergence seconds=" << fmtSci(swBd.stop(),3) << "\n";
    }

    return allPass;
}

} // namespace gsTetClip
