/** @file immersed_stokes_cylinder_example.cpp

    @brief Geometry layer of an immersed Stokes cylinder-in-channel benchmark
    (Schaefer-Turek 2D-1): the channel [xiIn,xiOut] x [xiBot,xiTop] = L x H,
    L = 2.2, H = 0.41, with a cylinder of radius R = 0.05 centred at
    c = (0.2, 0.2), the whole configuration rotated by an angle theta about c
    and immersed in a Cartesian B-spline background.

    Geometry is expressed through a rotated frame xi = R(-theta)(x - c) and a
    single scalar level set built from three "normalized factor" half-spaces
    (the two channel walls in each direction, f(xi1) and g(xi2), and, when a
    cylinder is present, the disk indicator h(xi1,xi2)) combined by the
    R-function conjunction a ^ b = a + b + sqrt(a^2+b^2), applied twice
    (f^g, then (f^g)^h). This R-function level set is C^infinity away from
    the zero set of BOTH conjoined arguments simultaneously: the derivative
    ratios a/s, b/s (s = sqrt(a^2+b^2)) have a genuine, unavoidable
    discontinuity exactly at the four channel corners (where f = g = 0
    simultaneously) -- the level set itself is still continuous there (the
    corner is a single point of the zero set, phi = 0 exactly, from either
    side), but its gradient is not unique. LevelSet::deriv_into resolves the
    resulting 0/0 by setting a/s = b/s = 0 there: any bounded value is a
    legitimate weak sub-gradient at an isolated non-differentiable point,
    and this choice keeps every downstream gsMatrix finite (see
    levelSetSelfCheck, item (c)).

    Background lattice: h_r = 0.41/(8*2^r), a Cartesian tensor mesh anchored
    at c (so c is always exactly a lattice vertex, independent of r and
    theta), padded by one cell on each side of the rotated channel's bounding
    box. Distances from c to the three geometric features that can align
    with a lattice line are 0.2/h = 160*2^r/41, 0.21/h = 168*2^r/41 and
    2.0/h = 1600*2^r/41; since 41 is odd none of these is ever a dyadic
    rational, so no wall coincides exactly with a lattice line or with a
    dyadic adaptive-subdivision face at theta in {0, 90} -- which matters
    because the Algoim surface rule hangs on exact tangency/coincidence
    (see immersed_stokes_skeleton_example.cpp for the same argument at a
    different box size).

    Quadrature: the default volume/interface rule is
    gsQuadrature::AlgoimAdaptiveRule (15), a direct adaptive dispatch
    (src/gsAssembler/gsQuadrature.h) which needs three knobs raised
    above their library defaults for this geometry. |grad phi| reaches about
    2*|grad h| = 4|xi|/R^2 near the cylinder and about |grad g| = H/0.205^2
    near the walls, both far above the rule's built-in box-classification
    Lipschitz constant of 1 (gsAlgoimAdaptiveRule::classify): an
    under-estimated Lipschitz constant misclassifies a genuinely cut child
    box as uncut and silently drops area. --lipschitz therefore defaults to
    200 and is checked at run time against 1.5 times the largest observed
    |grad phi| on the cut cells (maxGradOnCutCells); a violation aborts
    rather than reporting numbers nobody can trust. --maxDepth defaults to 6
    (the library default 0 is plain, non-adaptive Algoim, which drops
    corner and 90-degree-arc surface branches) and --indicatorTol defaults
    to 1e-9 (the library default 1e-2 only bounds the per-cut-cell error at
    about 1e-3 with O(1/h) cut cells along a wall).

    Exact surface rule for the four straight walls: a benchmark measurement
    on an axis-aligned square found that the adaptive rule recovers a
    re-entrant-corner branch only at FIRST order in maxDepth (perimeter
    error 0.867 at depth 0 down to 5.56e-3 at depth 6), so the four channel
    corners would cost the wall length O(h * 2^-maxDepth) of accuracy no
    matter how deep the adaptive tree goes. ImmersedSurfaceRule below
    replaces the adaptive rule for the piece-length (surface) integrals with
    an EXACT rule on the four rotated wall segments (Gauss-Legendre on the
    part of each segment clipped to the current element box, via a
    Liang-Barsky clip) and, when a cylinder is present, a second adaptive
    Algoim rule confined to a cylinder-only level set h(x) = (R^2 -
    |x-c|^2)/R^2. It is installed as a gsExprAssembler/gsExprEvaluator
    QuadratureFactory only for the duration of the surface integrals
    (SurfaceQuadratureScope) and must never be left installed across a
    volume call or a ghost/skeleton face loop, both of which consult the
    same factory (gsExprAssembler::assembleGhost/assembleSkeleton) and would
    silently receive a rule that is only valid on {phi == 0}.

    Non-owned level-set lifetime: gsImplicitTrimmedDomain stores its level
    set via make_shared_not_owned (src/gsDomain/gsImplicitTrimmedDomain.h),
    i.e. it does NOT keep the function alive. Every LevelSet/CylinderLevelSet
    object used to build a trimmed domain in this file is therefore a named
    local variable whose scope outlives that domain and every quadrature
    rule built from it.

    Example command lines:
      ./immersed_stokes_cylinder_example --study geometry --r0 0 --r1 4 --theta 0
      ./immersed_stokes_cylinder_example --study geometry --r0 0 --r1 4 --theta 17
      ./immersed_stokes_cylinder_example --study geometry --case channel --theta 17 --r0 0 --r1 2
      ./immersed_stokes_cylinder_example --plot -r 1 --theta 17

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>

#include <cmath>
#include <iomanip>
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

/// One rigid configuration of the Schaefer-Turek 2D-1 channel: the frame map
/// xi = R(-theta)(x - c), the R-function level set built from it, and the
/// piece classification/normal used by the boundary integrals. Kept small
/// and copyable (by value) so every gsFunction wrapper below can just store
/// one.
struct ChannelGeometry
{
    enum Piece { Inflow = 0, Outflow = 1, Bottom = 2, Top = 3, Cylinder = 4 };

    real_t theta;   // radians
    real_t ct, st;
    real_t cx, cy;
    bool   withCylinder;

    const real_t R     = 0.05;
    const real_t xiIn  = -0.2, xiOut = 2.0, xiBot = -0.2, xiTop = 0.21;

    ChannelGeometry(real_t thetaDeg, bool withCyl)
    : theta(thetaDeg * (real_t)EIGEN_PI / (real_t)180.0),
      ct(math::cos(theta)), st(math::sin(theta)),
      cx(0.2), cy(0.2), withCylinder(withCyl)
    { }

    void toXi(real_t x, real_t y, real_t & xi1, real_t & xi2) const
    {
        const real_t dx = x - cx, dy = y - cy;
        xi1 =  ct*dx + st*dy;
        xi2 = -st*dx + ct*dy;
    }

    void fromXi(real_t xi1, real_t xi2, real_t & x, real_t & y) const
    {
        x = cx + ct*xi1 - st*xi2;
        y = cy + st*xi1 + ct*xi2;
    }

    /// Rotates a xi-frame vector (gradient or normal component pair) into
    /// physical (x,y) components: e_xi1 = (ct,st), e_xi2 = (-st,ct).
    void rotateToPhysical(real_t g1, real_t g2, real_t & gx, real_t & gy) const
    { gx = ct*g1 - st*g2; gy = st*g1 + ct*g2; }

    /// R-conjunction a^b = a + b + sqrt(a^2+b^2) and its chain rule through
    /// two coordinates x1, x2 (da1 = da/dx1 etc.). At s == 0 exactly (the
    /// channel corners) both weights a/s, b/s are set to 0: phi is not
    /// differentiable there, and any bounded value is a legitimate choice.
    static void conjunction(real_t a, real_t b, real_t da1, real_t da2,
                            real_t db1, real_t db2,
                            real_t & val, real_t & dval1, real_t & dval2)
    {
        const real_t s  = math::sqrt(a*a + b*b);
        const real_t wa = (0.0 == s) ? 0.0 : a / s;
        const real_t wb = (0.0 == s) ? 0.0 : b / s;
        val   = a + b + s;
        dval1 = (1.0 + wa) * da1 + (1.0 + wb) * db1;
        dval2 = (1.0 + wa) * da2 + (1.0 + wb) * db2;
    }

    void evalPhi(real_t x, real_t y, real_t & phi) const
    {
        real_t xi1, xi2; toXi(x, y, xi1, xi2);
        const real_t f = (xi1 + 0.2) * (xi1 - 2.0) / 1.21;
        const real_t g = (xi2 + 0.2) * (xi2 - 0.21) / 0.042025;
        const real_t s1 = math::sqrt(f*f + g*g);
        const real_t phi1 = f + g + s1;
        if (!withCylinder) { phi = phi1; return; }
        const real_t h = (R*R - (xi1*xi1 + xi2*xi2)) / (R*R);
        const real_t s2 = math::sqrt(phi1*phi1 + h*h);
        phi = phi1 + h + s2;
    }

    void evalPhiGrad(real_t x, real_t y, real_t & phi, real_t & gx, real_t & gy) const
    {
        real_t xi1, xi2; toXi(x, y, xi1, xi2);
        const real_t f   = (xi1 + 0.2) * (xi1 - 2.0) / 1.21;
        const real_t df1 = (2.0*xi1 - 1.8) / 1.21;
        const real_t g   = (xi2 + 0.2) * (xi2 - 0.21) / 0.042025;
        const real_t dg2 = (2.0*xi2 - 0.01) / 0.042025;

        real_t phi1, dphi1_1, dphi1_2;
        conjunction(f, g, df1, 0.0, 0.0, dg2, phi1, dphi1_1, dphi1_2);

        real_t g1, g2;
        if (!withCylinder)
        {
            phi = phi1; g1 = dphi1_1; g2 = dphi1_2;
        }
        else
        {
            const real_t h   = (R*R - (xi1*xi1 + xi2*xi2)) / (R*R);
            const real_t dh1 = -2.0*xi1 / (R*R);
            const real_t dh2 = -2.0*xi2 / (R*R);
            conjunction(phi1, h, dphi1_1, dphi1_2, dh1, dh2, phi, g1, g2);
        }
        rotateToPhysical(g1, g2, gx, gy);
    }

    Piece classify(real_t x, real_t y) const
    {
        real_t xi1, xi2; toXi(x, y, xi1, xi2);
        const real_t f = (xi1 + 0.2) * (xi1 - 2.0) / 1.21;
        const real_t g = (xi2 + 0.2) * (xi2 - 0.21) / 0.042025;
        const real_t af = math::abs(f), ag = math::abs(g);

        if (withCylinder)
        {
            const real_t h  = (R*R - (xi1*xi1 + xi2*xi2)) / (R*R);
            const real_t ah = math::abs(h);
            if (af <= ag && af <= ah) return (xi1 < 0.9) ? Inflow : Outflow;
            if (ag <= ah)             return (xi2 < 0.005) ? Bottom : Top;
            return Cylinder;
        }
        if (af <= ag) return (xi1 < 0.9) ? Inflow : Outflow;
        return (xi2 < 0.005) ? Bottom : Top;
    }

    void outwardNormal(real_t x, real_t y, real_t & nx, real_t & ny) const
    {
        real_t xi1, xi2; toXi(x, y, xi1, xi2);
        real_t n1 = 0, n2 = 0;
        switch (classify(x, y))
        {
        case Inflow:  n1 = -1; n2 =  0; break;
        case Outflow: n1 =  1; n2 =  0; break;
        case Bottom:  n1 =  0; n2 = -1; break;
        case Top:     n1 =  0; n2 =  1; break;
        case Cylinder:
        {
            const real_t r = math::sqrt(xi1*xi1 + xi2*xi2);
            n1 = (r > 0.0) ? -xi1/r : 0.0;
            n2 = (r > 0.0) ? -xi2/r : 0.0;
            break;
        }
        }
        rotateToPhysical(n1, n2, nx, ny);
    }

    /// Physical corners of the rotated rectangle, ccw starting at
    /// (xiIn,xiBot): columns 0..3 = BL, BR, TR, TL. Consecutive columns
    /// (mod 4) bound, in order, the Bottom, Outflow, Top and Inflow walls.
    gsMatrix<real_t> corners() const
    {
        gsMatrix<real_t> C(2,4);
        real_t x, y;
        fromXi(xiIn,  xiBot, x, y); C(0,0) = x; C(1,0) = y;
        fromXi(xiOut, xiBot, x, y); C(0,1) = x; C(1,1) = y;
        fromXi(xiOut, xiTop, x, y); C(0,2) = x; C(1,2) = y;
        fromXi(xiIn,  xiTop, x, y); C(0,3) = x; C(1,3) = y;
        return C;
    }

    real_t exactArea()     const { return withCylinder ? (2.2*0.41 - M_PI*R*R) : 2.2*0.41; }
    real_t exactInflow()   const { return 0.41; }
    real_t exactOutflow()  const { return 0.41; }
    real_t exactBottom()   const { return 2.2; }
    real_t exactTop()      const { return 2.2; }
    real_t exactCylinder() const { return 2.0*M_PI*R; }
};

/// The channel/cylinder level set phi, evaluated at physical points (the
/// background geometry map is the identity). deriv_into is load-bearing:
/// the Algoim wrapper calls it for its interval (Taylor) bounds and for the
/// surface normal (optional/gsAlgoim/gsAlgoimFunctionWrapper.h).
class LevelSet : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(LevelSet)

    LevelSet(const ChannelGeometry & geom, const gsMatrix<real_t> & box)
    : m_geom(geom), m_box(box) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }
    gsMatrix<real_t> support() const override { return m_box; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            m_geom.evalPhi(u(0,k), u(1,k), result(0,k));
    }

    /// Row 0 = d(phi)/dx, row 1 = d(phi)/dy (gsFunction.h:122-164 layout for
    /// a scalar target: targetDim*domainDim = 2 rows, one column per point).
    void deriv_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(2, u.cols());
        real_t phi;
        for (index_t k = 0; k != u.cols(); ++k)
            m_geom.evalPhiGrad(u(0,k), u(1,k), phi, result(0,k), result(1,k));
    }

private:
    ChannelGeometry  m_geom;
    gsMatrix<real_t> m_box;
};

/// Cylinder-only level set h(x) = (R^2 - |x-c|^2)/R^2, negative in the
/// fluid (outside the disk). Used solely to build the ImmersedSurfaceRule's
/// cylinder sub-rule on a rotation-invariant geometry (no xi frame needed).
class CylinderLevelSet : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(CylinderLevelSet)

    CylinderLevelSet(real_t cx, real_t cy, real_t R, const gsMatrix<real_t> & box)
    : m_cx(cx), m_cy(cy), m_R(R), m_box(box) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }
    gsMatrix<real_t> support() const override { return m_box; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const real_t dx = u(0,k)-m_cx, dy = u(1,k)-m_cy;
            result(0,k) = (m_R*m_R - (dx*dx + dy*dy)) / (m_R*m_R);
        }
    }

    void deriv_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(2, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            result(0,k) = -2.0*(u(0,k)-m_cx) / (m_R*m_R);
            result(1,k) = -2.0*(u(1,k)-m_cy) / (m_R*m_R);
        }
    }

private:
    real_t m_cx, m_cy, m_R;
    gsMatrix<real_t> m_box;
};

/// Unit outward normal of classify(x), fed to getCoeff as n_imm by an
/// immersed Nitsche form. Only used at Space == 0 (a data coefficient), so
/// no deriv_into.
class PieceNormal : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(PieceNormal)
    PieceNormal(const ChannelGeometry & geom) : m_geom(geom) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 2; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(2, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            m_geom.outwardNormal(u(0,k), u(1,k), result(0,k), result(1,k));
    }

private:
    ChannelGeometry m_geom;
};

/// chi_D = 1 on Inflow, Bottom, Top and Cylinder (the Dirichlet part of the
/// immersed boundary); 0 on Outflow.
class DirichletIndicator : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(DirichletIndicator)
    DirichletIndicator(const ChannelGeometry & geom) : m_geom(geom) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            result(0,k) = (ChannelGeometry::Outflow == m_geom.classify(u(0,k),u(1,k))) ? 0.0 : 1.0;
    }

private:
    ChannelGeometry m_geom;
};

/// 1 where classify(x) == piece, else 0. Used to isolate one boundary
/// piece's length under a single integralBdr call over the whole immersed
/// boundary.
class PieceIndicator : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(PieceIndicator)
    PieceIndicator(const ChannelGeometry & geom, ChannelGeometry::Piece piece)
    : m_geom(geom), m_piece(piece) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            result(0,k) = (m_piece == m_geom.classify(u(0,k),u(1,k))) ? 1.0 : 0.0;
    }

private:
    ChannelGeometry m_geom;
    ChannelGeometry::Piece m_piece;
};

/// (real_t)classify(x), for diagnostics and Paraview plotting.
class PieceId : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(PieceId)
    PieceId(const ChannelGeometry & geom) : m_geom(geom) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            result(0,k) = (real_t)m_geom.classify(u(0,k),u(1,k));
    }

private:
    ChannelGeometry m_geom;
};

/// The Cartesian background lattice, anchored so that c is always exactly a
/// lattice vertex, and padded by one cell on each side of the rotated
/// channel's axis-aligned bounding box.
struct Background
{
    gsMultiPatch<real_t> mp;
    gsMultiBasis<>        dbasis;
    real_t h;
    index_t nx, ny;
    gsMatrix<real_t> box;   // 2x2: column 0 = lower, column 1 = upper
};

Background makeBackground(const ChannelGeometry & geom, index_t k, index_t r)
{
    const real_t h = 0.41 / (8.0 * (real_t)(index_t(1) << r));

    const gsMatrix<real_t> C = geom.corners();
    const real_t xmin = C.row(0).minCoeff(), xmax = C.row(0).maxCoeff();
    const real_t ymin = C.row(1).minCoeff(), ymax = C.row(1).maxCoeff();

    const index_t i0 = (index_t)std::floor((xmin - geom.cx)/h) - 1;
    const index_t i1 = (index_t)std::ceil ((xmax - geom.cx)/h) + 1;
    const index_t j0 = (index_t)std::floor((ymin - geom.cy)/h) - 1;
    const index_t j1 = (index_t)std::ceil ((ymax - geom.cy)/h) + 1;

    const index_t nx = i1 - i0, ny = j1 - j0;
    const real_t x0 = geom.cx + (real_t)i0*h, x1 = x0 + (real_t)nx*h;
    const real_t y0 = geom.cy + (real_t)j0*h, y1 = y0 + (real_t)ny*h;

    gsKnotVector<real_t> kvx(x0, x1, nx - 1, k + 1);
    gsKnotVector<real_t> kvy(y0, y1, ny - 1, k + 1);
    gsTensorBSplineBasis<2,real_t> tbs(kvx, kvy);

    Background bg;
    bg.dbasis = gsMultiBasis<>(tbs);
    bg.mp = gsMultiPatch<real_t>(*gsNurbsCreator<real_t>::BSplineRectangleWithPara(x0, y0, x1, y1));
    bg.h = h; bg.nx = nx; bg.ny = ny;
    bg.box.resize(2,2);
    bg.box(0,0) = x0; bg.box(1,0) = y0; bg.box(0,1) = x1; bg.box(1,1) = y1;
    return bg;
}

/// Registers the immersed-quadrature options on a fresh gsOptionList
/// (gsExprAssembler and gsExprEvaluator both qualify): sets quRule, adds
/// quDim (volume by default, -1) plus the three adaptive knobs this
/// geometry needs above their library defaults.
void setQuadratureOptions(gsOptionList & opt, index_t quRule, index_t maxDepth,
                          real_t indicatorTol, real_t lipschitz)
{
    opt.setInt ("quRule", quRule);
    opt.addInt ("quDim", "Surface (phi==0) quadrature selector", -1);
    opt.addInt ("maxDepth", "Adaptive quadrature maximum subdivision depth", maxDepth);
    opt.addReal("indicatorTol", "integralChange acceptance tolerance, relative to the box measure",
                indicatorTol);
    opt.addReal("LipschitzConstant", "Lipschitz constant for adaptive box classification", lipschitz);
}

/// Four checks certifying the level set before any quadrature is trusted:
/// (a) sign at three reference points, (b) deriv_into against central
/// differences (h=1e-6) at six interior points, (c) deriv_into stays finite
/// at the four exact channel corners (the NaN-guard case), (d) phi is zero
/// to 1e-12 at eight points on the pieces, with classify() agreeing.
bool levelSetSelfCheck(const ChannelGeometry & geom)
{
    gsMatrix<real_t> hugeBox(2,2);
    hugeBox << -1e3, 1e3, -1e3, 1e3;
    LevelSet phi(geom, hugeBox);

    real_t x, y, val;

    // (a) sign at three reference points.
    geom.fromXi(1.0, 0.0, x, y); geom.evalPhi(x, y, val);
    if (!(val < 0.0))
    { gsWarn << "levelSetSelfCheck (a) FAILED: phi(c+R(1,0)) = " << val << " is not < 0\n"; return false; }

    geom.fromXi(3.0, 0.0, x, y); geom.evalPhi(x, y, val);
    if (!(val > 0.0))
    { gsWarn << "levelSetSelfCheck (a) FAILED: phi(c+R(3,0)) = " << val << " is not > 0\n"; return false; }

    if (geom.withCylinder)
    {
        geom.evalPhi(geom.cx, geom.cy, val);
        if (!(val > 0.0))
        { gsWarn << "levelSetSelfCheck (a) FAILED: phi(c) = " << val << " is not > 0\n"; return false; }
    }

    // (b) deriv_into against central differences, step 1e-6, at six
    // interior, non-corner points.
    const real_t xiPts[6][2] = {
        {0.5, 0.1}, {1.5, -0.1}, {0.08, 0.0}, {-0.1, 0.15}, {1.9, 0.2}, {0.0, -0.07}
    };
    const real_t fdStep = 1e-6;
    real_t worstB = 0;
    for (index_t i = 0; i != 6; ++i)
    {
        geom.fromXi(xiPts[i][0], xiPts[i][1], x, y);
        gsMatrix<real_t> pt(2,1), grad;
        pt << x, y;
        phi.deriv_into(pt, grad);

        real_t vxp, vxm, vyp, vym;
        geom.evalPhi(x+fdStep, y, vxp); geom.evalPhi(x-fdStep, y, vxm);
        geom.evalPhi(x, y+fdStep, vyp); geom.evalPhi(x, y-fdStep, vym);
        const real_t fdx = (vxp-vxm)/(2.0*fdStep), fdy = (vyp-vym)/(2.0*fdStep);

        const real_t denom = math::max(math::sqrt(fdx*fdx+fdy*fdy), (real_t)1.0);
        const real_t rel = math::sqrt(math::pow(grad(0,0)-fdx,2) + math::pow(grad(1,0)-fdy,2)) / denom;
        worstB = math::max(worstB, rel);
    }
    if (worstB > 1e-6)
    { gsWarn << "levelSetSelfCheck (b) FAILED: worst relative deriv_into mismatch = "
             << fmtSci(worstB) << " > 1e-6\n"; return false; }

    // (c) deriv_into stays finite at the four exact rotated channel
    // corners: f == 0 and g == 0 simultaneously there (s == 0 in the
    // conjunction), exercising the NaN guard.
    const gsMatrix<real_t> C = geom.corners();
    for (index_t i = 0; i != 4; ++i)
    {
        gsMatrix<real_t> pt(2,1), grad;
        pt << C(0,i), C(1,i);
        phi.deriv_into(pt, grad);
        if (!grad.allFinite())
        { gsWarn << "levelSetSelfCheck (c) FAILED: deriv_into not finite at corner " << i << "\n"; return false; }
    }

    // (d) phi == 0 to 1e-12 at eight points on the pieces, and classify()
    // agrees at each.
    real_t worstD = 0;
    auto checkPiece = [&](real_t xi1, real_t xi2, ChannelGeometry::Piece expect) -> bool
    {
        real_t px, py;
        geom.fromXi(xi1, xi2, px, py);
        real_t v; geom.evalPhi(px, py, v);
        worstD = math::max(worstD, math::abs(v));
        if (math::abs(v) > 1e-12)
        { gsWarn << "levelSetSelfCheck (d) FAILED: |phi| = " << fmtSci(math::abs(v))
                 << " > 1e-12 at piece " << (int)expect << "\n"; return false; }
        if (geom.classify(px, py) != expect)
        { gsWarn << "levelSetSelfCheck (d) FAILED: classify() mismatch at piece " << (int)expect << "\n";
          return false; }
        return true;
    };
    const real_t mid1 = 0.5*(geom.xiIn + geom.xiOut), mid2 = 0.5*(geom.xiBot + geom.xiTop);
    if (!checkPiece(geom.xiIn,  mid2, ChannelGeometry::Inflow))  return false;
    if (!checkPiece(geom.xiOut, mid2, ChannelGeometry::Outflow)) return false;
    if (!checkPiece(mid1, geom.xiBot, ChannelGeometry::Bottom))  return false;
    if (!checkPiece(mid1, geom.xiTop, ChannelGeometry::Top))     return false;
    if (geom.withCylinder)
    {
        for (index_t a = 0; a != 4; ++a)
        {
            const real_t alpha = (real_t)a * 0.5 * M_PI;
            if (!checkPiece(geom.R*math::cos(alpha), geom.R*math::sin(alpha), ChannelGeometry::Cylinder))
                return false;
        }
    }

    gsInfo << "levelSetSelfCheck theta=" << (geom.theta*180.0/M_PI) << " withCylinder=" << geom.withCylinder
           << ": worst deriv mismatch=" << fmtSci(worstB) << ", worst |phi| on pieces=" << fmtSci(worstD)
           << ": OK\n";
    return true;
}

/// Maximum |grad phi| over 5x5 Lobatto samples of every cut cell of \a dom,
/// the same sampling gsTrimmedDomain::_classifyLeaf uses internally
/// (src/gsDomain/gsTrimmedDomain.h). Governs the Lipschitz-constant check:
/// gsAlgoimAdaptiveRule::classify treats a sub-box as uncut whenever
/// |phi(mid)| > LipschitzConstant * 0.5*diag, so an under-estimated
/// constant silently drops area.
real_t maxGradOnCutCells(const LevelSet & phi, const gsImplicitTrimmedDomain<2,real_t> & dom)
{
    gsVector<index_t> numNodes(2); numNodes.setConstant(5);
    gsLobattoRule<real_t> QR(numNodes);
    gsMatrix<real_t> pts, grad; gsVector<real_t> wts;

    real_t maxGrad = 0;
    gsDomain<real_t>::iterator eCut = dom.endBdr(boundary::none);
    for (gsDomain<real_t>::iterator it = dom.beginBdr(boundary::none); it < eCut; ++it)
    {
        QR.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
        phi.deriv_into(pts, grad);
        for (index_t c = 0; c != grad.cols(); ++c)
            maxGrad = math::max(maxGrad, math::sqrt(grad(0,c)*grad(0,c) + grad(1,c)*grad(1,c)));
    }
    return maxGrad;
}

/// Liang-Barsky clip of the parametric segment P(t) = P0 + t*D, t in [0,1],
/// against the axis-aligned box [lower,upper]. Returns false when the
/// (possibly empty, possibly degenerate) intersection has zero length.
bool clipSegmentToBox(const gsVector<real_t,2> & P0, const gsVector<real_t,2> & D,
                      const gsVector<real_t> & lower, const gsVector<real_t> & upper,
                      real_t & t0, real_t & t1)
{
    real_t tE = 0.0, tL = 1.0;
    const real_t p[4] = { -D[0], D[0], -D[1], D[1] };
    const real_t q[4] = { P0[0]-lower[0], upper[0]-P0[0], P0[1]-lower[1], upper[1]-P0[1] };
    for (index_t i = 0; i != 4; ++i)
    {
        if (0.0 == p[i])
        {
            if (q[i] < 0.0) return false;   // parallel and outside this edge
        }
        else
        {
            const real_t r = q[i] / p[i];
            if (p[i] < 0.0) { if (r > tL) return false; if (r > tE) tE = r; }
            else            { if (r < tE) return false; if (r < tL) tL = r; }
        }
    }
    if (!(tE < tL)) return false;
    t0 = tE; t1 = tL;
    return true;
}

/// Exact surface rule for the immersed boundary, installed only for the
/// piece-length (surface) integrals (see the file header). Concurrency
/// invariant: the QuadratureFactory (makeSurfaceQuadratureFactory) returns a
/// FRESH ImmersedSurfaceRule on every call, and gsExprAssembler/
/// gsExprEvaluator call the factory once per OpenMP worker (or once per
/// boundary loop, for the sequential gsExprEvaluator::computeBdr_impl this
/// file actually drives), never sharing one rule instance across threads.
/// The only state shared BETWEEN instances is the read-only cylinder level
/// set/domain (via shared_ptr, built once per refinement level -- see
/// runGeometryStudy). A single ImmersedSurfaceRule instance is NOT safe for
/// concurrent mapTo() calls: m_cylRule wraps a gsAlgoimAdaptiveDirectRule
/// whose internal Stats counters are mutable, non-atomic state.
class ImmersedSurfaceRule : public gsQuadRule<real_t>
{
public:
    /// \a cylDomain is nullptr for --case channel; otherwise it must already
    /// be classified (built once per refinement level, shared across every
    /// factory call of that level -- see makeSurfaceQuadratureFactory).
    /// Re-classifying the whole background lattice inside every one of the
    /// (#pieces + 1) factory calls per level would multiply the cost of
    /// gsImplicitTrimmedDomain's construction by that many, for no benefit:
    /// the classification does not depend on which piece is being
    /// integrated.
    ImmersedSurfaceRule(const ChannelGeometry & geom, index_t k,
                        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > cylDomain,
                        const gsOptionList & options)
    : m_cx(geom.cx), m_cy(geom.cy), m_R(geom.R), m_nG(2*k+2), m_corners(geom.corners()),
      m_cylDomain(give(cylDomain))
    {
        if (m_cylDomain)
        {
            gsOptionList cylOpts = options;   // carries quRule, quDim==2 and the adaptive knobs
            gsVector<short_t> degs(2); degs << k, k;
            m_cylRule = gsQuadrature::getPtr(*m_cylDomain, cylOpts, 0, degs);
        }
    }

    using gsQuadRule<real_t>::mapTo;
    void mapTo(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
              gsMatrix<real_t> & nodes, gsVector<real_t> & weights) const override
    {
        std::vector<gsMatrix<real_t>> nodeParts;
        std::vector<gsVector<real_t>> wtParts;

        for (index_t e = 0; e != 4; ++e)
        {
            const gsVector<real_t,2> P0 = m_corners.col(e);
            const gsVector<real_t,2> P1 = m_corners.col((e+1) % 4);
            const gsVector<real_t,2> D  = P1 - P0;

            real_t t0, t1;
            if (!clipSegmentToBox(P0, D, lower, upper, t0, t1)) continue;

            // nG = 2k+2: a rotated line intersected with a product of two
            // degree-k tensor B-splines is a polynomial of degree <= 4k in
            // the line parameter, and Gauss with n points is exact to
            // degree 2n-1 >= 4k+1.
            gsGaussRule<real_t> gr(m_nG);
            gsMatrix<real_t> tNodes; gsVector<real_t> tWeights;
            gr.mapTo(t0, t1, tNodes, tWeights);

            gsMatrix<real_t> pts(2, tNodes.cols());
            for (index_t i = 0; i != tNodes.cols(); ++i)
                pts.col(i) = P0 + tNodes(0,i) * D;
            gsVector<real_t> wts = tWeights * D.norm();

            nodeParts.push_back(give(pts));
            wtParts.push_back(give(wts));
        }

        if (m_cylRule)
        {
            const real_t diag = (upper - lower).norm();
            const gsVector<real_t,2> mid = 0.5*(lower + upper);
            const real_t dist = math::sqrt(math::pow(mid[0]-m_cx,2) + math::pow(mid[1]-m_cy,2));
            if (dist >= m_R - diag && dist <= m_R + diag)
            {
                gsMatrix<real_t> cylNodes; gsVector<real_t> cylWts;
                m_cylRule->mapTo(lower, upper, cylNodes, cylWts);
                if (cylNodes.cols() > 0)
                {
                    nodeParts.push_back(give(cylNodes));
                    wtParts.push_back(give(cylWts));
                }
            }
        }

        index_t total = 0;
        for (const auto & p : nodeParts) total += p.cols();
        nodes.resize(2, total);
        weights.resize(total);
        index_t col = 0;
        for (size_t i = 0; i != nodeParts.size(); ++i)
        {
            nodes.middleCols(col, nodeParts[i].cols()) = nodeParts[i];
            weights.segment(col, wtParts[i].size()) = wtParts[i];
            col += nodeParts[i].cols();
        }
    }

private:
    real_t m_cx, m_cy, m_R;
    index_t m_nG;
    gsMatrix<real_t> m_corners;
    memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > m_cylDomain;
    gsQuadRule<real_t>::uPtr m_cylRule;
};

/// Builds the QuadratureFactory that returns a fresh ImmersedSurfaceRule.
/// \a cylDomain (nullptr for --case channel) is built ONCE per refinement
/// level by the caller and shared by every factory call of that level.
gsExprEvaluator<>::QuadratureFactory
makeSurfaceQuadratureFactory(const ChannelGeometry & geom, index_t k,
                             memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > cylDomain)
{
    return [geom, k, cylDomain]
           (const gsDomain<real_t> &, const gsBasis<real_t> *,
            const gsOptionList & options, index_t, short_t,
            const gsVector<short_t> &) -> gsQuadRule<real_t>::uPtr
    {
        return gsQuadRule<real_t>::uPtr(new ImmersedSurfaceRule(geom, k, cylDomain, options));
    };
}

/// RAII scope installing a custom QuadratureFactory on an assembler or
/// evaluator, valid ONLY for surface (quDim==2) calls made within its
/// lifetime; restores plain option-driven volume quadrature on destruction.
/// Never let this scope's lifetime span a volume integral or a
/// ghost/skeleton face loop: both consult the very same factory
/// (gsExprAssembler::assembleGhost/assembleSkeleton,
/// src/gsAssembler/gsExprAssembler.h) and would silently receive a rule
/// that is only valid on {phi == 0}.
template<class ExprObj>
class SurfaceQuadratureScope
{
public:
    SurfaceQuadratureScope(ExprObj & obj, typename ExprObj::QuadratureFactory factory)
    : m_obj(obj)
    {
        m_obj.setQuadratureFactory(give(factory));
        m_obj.options().setInt("quDim", 2);
    }
    ~SurfaceQuadratureScope()
    {
        m_obj.options().setInt("quDim", -1);
        m_obj.clearQuadratureFactory();
    }
    SurfaceQuadratureScope(const SurfaceQuadratureScope&) = delete;
    SurfaceQuadratureScope& operator=(const SurfaceQuadratureScope&) = delete;
private:
    ExprObj & m_obj;
};

void runGeometryStudy(const std::string & caseName, real_t thetaDeg, index_t k,
                     index_t r0, index_t r1, index_t quRule, index_t maxDepth,
                     real_t indicatorTol, real_t lipschitz)
{
    const bool withCylinder = ("cylinder" == caseName);
    ChannelGeometry geom(thetaDeg, withCylinder);

    gsInfo << "\n=== geometry  theta=" << thetaDeg << " k=" << k << " case=" << caseName
           << "  quRule=" << quRule << " maxDepth=" << maxDepth
           << " indicatorTol=" << fmtSci(indicatorTol) << " lipschitz=" << lipschitz << " ===\n";

    gsInfo << std::right << std::setw(3) << "r" << std::setw(11) << "h"
           << std::setw(10) << "nx x ny" << std::setw(7) << "nCut"
           << std::setw(14) << "areaErr" << std::setw(14) << "inErr" << std::setw(14) << "outErr"
           << std::setw(14) << "botErr" << std::setw(14) << "topErr";
    if (withCylinder) gsInfo << std::setw(14) << "cylErr";
    gsInfo << std::setw(14) << "sumErr" << std::setw(14) << "etaMin"
           << std::setw(14) << "consist" << std::setw(10) << "maxGrad" << std::setw(9) << "time(s)\n";

    for (index_t r = r0; r <= r1; ++r)
    {
        gsStopwatch clk;

        Background bg = makeBackground(geom, k, r);
        LevelSet phi(geom, bg.box);

        gsTensorBSplineBasis<2,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&bg.dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");
        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > tr_domain =
            memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(phi, *tbsPtr));
        GISMO_ENSURE(1 == tr_domain->numLevels(), "gsTrimmedDomain requires a single kd-tree level.");
        const index_t nCut = static_cast<index_t>(tr_domain->numElementsBdr(boundary::none));

        const real_t maxGrad = maxGradOnCutCells(phi, *tr_domain);
        GISMO_ENSURE(lipschitz >= 1.5*maxGrad,
                    "--lipschitz " << lipschitz << " is below 1.5*maxGrad = " << 1.5*maxGrad
                    << ": raise --lipschitz, or gsAlgoimAdaptiveRule::classify may misclassify a "
                       "genuinely cut sub-box as uncut and silently drop area.");

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(bg.dbasis);
        ev.setIntegrationDomain(tr_domain);
        setQuadratureOptions(ev.options(), quRule, maxDepth, indicatorTol, lipschitz);

        auto G = ev.getMap(bg.mp);
        const real_t area = ev.integral(meas(G));

        // eta_min: manual per-cut-cell loop with the SAME volume rule ev
        // uses (empty degs so quadratureDegrees() -- which finds none in
        // meas(G) -- is reproduced exactly).
        const gsVector<short_t> degs;
        gsQuadRule<real_t>::uPtr qr = gsQuadrature::getPtr(*tr_domain, ev.options(), -1, degs);
        gsMatrix<real_t> qnodes; gsVector<real_t> qwts;
        real_t etaMin = 1.0, cutArea = 0.0;
        gsVector<real_t> etaMinCorner;
        gsDomain<real_t>::iterator eCut = tr_domain->endBdr(boundary::none);
        for (gsDomain<real_t>::iterator it = tr_domain->beginBdr(boundary::none); it < eCut; ++it)
        {
            qr->mapTo(it.lowerCorner(), it.upperCorner(), qnodes, qwts);
            const real_t cellArea = qwts.sum();
            cutArea += cellArea;
            const real_t eta = cellArea / (bg.h*bg.h);
            if (eta < etaMin) { etaMin = eta; etaMinCorner = it.lowerCorner(); }
        }
        index_t nInterior = 0;
        gsDomain<real_t>::iterator eInt = tr_domain->end<InteriorSign>();
        for (gsDomain<real_t>::iterator it = tr_domain->beginInterior(); it < eInt; ++it)
            ++nInterior;
        const real_t consistArea = cutArea + (real_t)nInterior * bg.h*bg.h;
        const real_t consistRel  = math::abs(consistArea - area) / math::max(area, (real_t)1e-30);

        // Surface (piece-length) integrals: exact walls + adaptive cylinder,
        // installed only for the calls inside this block.
        std::vector<patchSide> bdr_immersed{patchSide(0, boundary::none)};
        PieceNormal pieceNormal(geom);

        std::vector<ChannelGeometry::Piece> pieces = {
            ChannelGeometry::Inflow, ChannelGeometry::Outflow,
            ChannelGeometry::Bottom, ChannelGeometry::Top
        };
        if (withCylinder) pieces.push_back(ChannelGeometry::Cylinder);

        std::vector<PieceIndicator> indicators;
        for (auto p : pieces) indicators.emplace_back(geom, p);

        // Cylinder-only level set/domain, built ONCE for this r and shared
        // by every surface-quadrature factory call below (see
        // ImmersedSurfaceRule and makeSurfaceQuadratureFactory): the
        // non-owned level-set lifetime rule (file header) applies here too,
        // so cylLevelSet must outlive cylDomain and every rule built from it.
        CylinderLevelSet cylLevelSet(geom.cx, geom.cy, geom.R, bg.box);
        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > cylDomain;
        if (withCylinder)
        {
            cylDomain = memory::make_shared(
                new gsImplicitTrimmedDomain<2,real_t>(cylLevelSet, *tbsPtr));
            GISMO_ENSURE(1 == cylDomain->numLevels(),
                        "gsTrimmedDomain requires a single kd-tree level.");
        }

        real_t total = 0, sumLen = 0;
        std::vector<real_t> lens(pieces.size(), 0.0);
        {
            SurfaceQuadratureScope<gsExprEvaluator<> > scope(
                ev, makeSurfaceQuadratureFactory(geom, k, cylDomain));

            auto n_imm    = ev.getVariable(pieceNormal, G);
            auto surfMeas = meas(G) * (jac(G).inv().tr() * n_imm).norm();

            for (size_t i = 0; i != pieces.size(); ++i)
            {
                auto ind = ev.getVariable(indicators[i], G);
                lens[i] = ev.integralBdr(ind * surfMeas, bdr_immersed);
                sumLen += lens[i];
            }
            total = ev.integralBdr(surfMeas, bdr_immersed);
        }

        const real_t areaErr = math::abs(area - geom.exactArea()) / geom.exactArea();
        const real_t exactWall[4] = { geom.exactInflow(), geom.exactOutflow(),
                                      geom.exactBottom(), geom.exactTop() };
        real_t pieceErr[5] = {0,0,0,0,0};
        for (size_t i = 0; i != pieces.size(); ++i)
        {
            const real_t exact = (pieces[i] == ChannelGeometry::Cylinder)
                                ? geom.exactCylinder() : exactWall[pieces[i]];
            pieceErr[i] = math::abs(lens[i] - exact) / exact;
        }
        const real_t exactTotal = geom.exactInflow() + geom.exactOutflow()
                                + geom.exactBottom() + geom.exactTop()
                                + (withCylinder ? geom.exactCylinder() : 0.0);
        const real_t sumErr = math::abs(sumLen - exactTotal) / exactTotal;
        const real_t sumVsTotal = math::abs(sumLen - total) / math::max(total, (real_t)1e-30);

        const real_t elapsed = clk.stop();

        gsInfo << std::setw(3) << r << std::setw(11) << fmtSci(bg.h,3)
               << std::setw(4) << bg.nx << " x" << std::setw(5) << bg.ny
               << std::setw(7) << nCut << std::setw(14) << fmtSci(areaErr);
        for (size_t i = 0; i != 4; ++i) gsInfo << std::setw(14) << fmtSci(pieceErr[i]);
        if (withCylinder) gsInfo << std::setw(14) << fmtSci(pieceErr[4]);
        gsInfo << std::setw(14) << fmtSci(sumErr) << std::setw(14) << fmtSci(etaMin)
               << std::setw(14) << fmtSci(consistRel) << std::setw(10) << fmtSci(maxGrad,2)
               << std::setw(9) << fmtSci(elapsed,2) << "\n";

        gsInfo << "    sum-vs-total rel diff = " << fmtSci(sumVsTotal)
               << " (sum=" << fmtSci(sumLen) << ", total=" << fmtSci(total) << ")"
               << ", eta_min at (" << etaMinCorner[0] << ", " << etaMinCorner[1] << ")\n";
    }
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    index_t degree       = 2;
    index_t refine       = 1;
    index_t r0           = 0;
    index_t r1           = 3;
    real_t  theta        = 0;
    std::string study    = "geometry";
    std::string caseName = "cylinder";
    index_t quRule       = gsQuadrature::AlgoimAdaptiveRule;
    index_t maxDepth     = 6;
    real_t  indicatorTol = 1e-9;
    real_t  lipschitz    = 200;
    bool    plot         = false;
    std::string outFolder = "output_immersed_stokes_cylinder";

    gsCmdLine cmd("Geometry layer of an immersed Stokes cylinder-in-channel benchmark "
                 "(Schaefer-Turek 2D-1), rotated by theta about the cylinder centre and "
                 "immersed in a Cartesian B-spline background.");
    cmd.addInt   ("k", "degree",       "Spline degree of the background basis", degree);
    cmd.addInt   ("r", "refine",       "Refinement level for single-level runs and --plot", refine);
    cmd.addInt   ("",  "r0",           "Study range: first refinement level", r0);
    cmd.addInt   ("",  "r1",           "Study range: last refinement level", r1);
    cmd.addReal  ("",  "theta",        "Rotation angle in degrees (not reduced modulo 360)", theta);
    cmd.addString("",  "study",        "Study: geometry (the only value currently implemented)", study);
    cmd.addString("",  "case",         "cylinder or channel", caseName);
    cmd.addInt   ("q", "quRule",       "Immersed rule: 11=CutCell, 12=Algoim, 13=Octree, "
                                       "15=AlgoimAdaptive (14=MomentFitting is rejected)", quRule);
    cmd.addInt   ("",  "maxDepth",     "Adaptive quadrature maximum subdivision depth", maxDepth);
    cmd.addReal  ("",  "indicatorTol", "integralChange acceptance tolerance", indicatorTol);
    cmd.addReal  ("",  "lipschitz",    "Lipschitz constant for adaptive box classification", lipschitz);
    cmd.addSwitch("plot", "Write Paraview output at refinement level -r", plot);
    cmd.addString("o", "output",       "Output folder", outFolder);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if ("geometry" != study)
    { gsWarn << "--study must be geometry (no other study is currently implemented)\n"; return EXIT_FAILURE; }

    if ("cylinder" != caseName && "channel" != caseName)
    { gsWarn << "--case must be cylinder or channel\n"; return EXIT_FAILURE; }

    if (degree < 1)
    { gsWarn << "-k/--degree must be >= 1\n"; return EXIT_FAILURE; }

    if (r0 < 0 || r1 < r0)
    { gsWarn << "--r0/--r1 must satisfy 0 <= r0 <= r1\n"; return EXIT_FAILURE; }

    if (refine < 0)
    { gsWarn << "-r/--refine must be >= 0\n"; return EXIT_FAILURE; }

    if (quRule != gsQuadrature::CutCellRule && quRule != gsQuadrature::AlgoimRule
        && quRule != gsQuadrature::OctreeRule && quRule != gsQuadrature::AlgoimAdaptiveRule
        && quRule != gsQuadrature::MomentFittingRule)
    { gsWarn << "-q/--quRule must be an immersed rule (11, 12, 13 or 15)\n"; return EXIT_FAILURE; }
    if (gsQuadrature::MomentFittingRule == quRule)
    {
        // gsQuadrature::makeMomentFittingPtr GISMO_ENSUREs quDim < 0: moment
        // fitting compresses points onto a volume tensor grid, which is
        // meaningless on the lower-dimensional zero level set, and this
        // driver's piece-length integrals unconditionally assemble with
        // quDim == 2 (via ImmersedSurfaceRule, itself independent of "-q"
        // except for the cylinder sub-rule).
        gsWarn << "-q/--quRule 14 (MomentFittingRule) cannot be used: it refuses surface "
                  "quadrature (quDim >= 0), and the geometry study's surface quadrature for the "
                  "piece lengths assembles with quDim == 2. Use 11, 12, 13 or 15.\n";
        return EXIT_FAILURE;
    }

    // Startup self-check on the geometry this run actually uses (so
    // --theta 360 certifies phi at theta=360, not a stand-in for theta=0).
    {
        ChannelGeometry g(theta, "cylinder" == caseName);
        if (!levelSetSelfCheck(g))
        { gsWarn << "levelSetSelfCheck FAILED for theta=" << theta << " case=" << caseName << "\n";
          return EXIT_FAILURE; }
    }

    std::string outPath = outFolder;
    if (outPath.empty()) outPath = "output_immersed_stokes_cylinder";
    const bool isAbsolutePath = (!outPath.empty() && outPath[0] == '/');
    if (!isAbsolutePath) outPath = gsFileManager::getCurrentPath() + "/" + outPath;
    const std::string out = gsFileManager::getCanonicRepresentation(outPath);
    if (plot) gsFileManager::mkdir(out);

    runGeometryStudy(caseName, theta, degree, r0, r1, quRule, maxDepth, indicatorTol, lipschitz);

    if (plot)
    {
        const bool withCylinder = ("cylinder" == caseName);
        ChannelGeometry geom(theta, withCylinder);
        Background bg = makeBackground(geom, degree, refine);
        LevelSet phi(geom, bg.box);
        PieceId pieceId(geom);

        gsWriteParaview(phi, bg.box, out + "/levelset", 40000);
        gsWriteParaview(pieceId, bg.box, out + "/pieceid", 40000);
        gsMesh<> mesh(bg.dbasis.basis(0));
        gsWriteParaview(mesh, out + "/background_mesh");
        gsInfo << "Paraview output written to " << out << "\n";
    }

    return EXIT_SUCCESS;
}
