/** @file stokes_cylinder_bodyfitted_example.cpp

    @brief Body-fitted Taylor-Hood Stokes reference for the Schaefer-Turek
    2D-1 benchmark ("flow around a cylinder").

    Solves steady Stokes (no convection) on the 5-patch NURBS Schaefer-Turek
    geometry,

        -div sigma(u,p) = 0,   div u = 0   in Omega = (0,2.2) x (0,0.41) \ B_R(c),
        sigma(u,p) = 2*mu*grad^s(u) - p*I,   mu = 1e-3,

    c = (0.2, 0.2), R = 0.05, with isogeometric Taylor-Hood elements
    (pressure degree k_ref, velocity degree k_ref+1, same regularity) and
    strong Dirichlet conditions:

      inflow  x = 0        : u = (4*Um*y*(H-y)/H^2, 0), Um = 0.3, H = 0.41
      walls   y = 0, 0.41   : u = 0
      cylinder              : u = 0
      outflow x = 2.2        : natural (traction-free); no pressure pin --
                                the outflow removes the constant-pressure
                                kernel.

    Boundary sides (G+Smo convention west=1,east=2,south=3,north=4; verified
    against filedata/planar/flow_around_cylinder.xml <boundary> block):

      role                | (patch, side)
      --------------------|-----------------------------------------
      inflow   x = 0       | (0, north)
      top wall y = 0.41    | (1, north), (4, south)
      bottom wall y = 0    | (3, north), (4, north)
      cylinder             | (0, south), (1, south), (2, south), (3, south)
      outflow  x = 2.2      | (4, east)

    Patch 4 is parametrized left-handed (det J < 0): its v=0 (south) row is
    the physical line y = 0.41 (top wall) and its v=1 (north) row is y = 0
    (bottom wall). meas(G) and nv(G) absorb the sign automatically
    (nv(G) flips by sgn(det J), gsCore/gsFunction.hpp); the area and
    cylinder-length startup checks below would fail if this were wrong.

    Geometry correction ("snap"). The XML coefficients are truncated to 5
    digits: the arc corners are stored as 0.16464/0.23536 (exact
    0.2 -+ R/sqrt(2) = 0.164644660940673 / 0.235355339059327) and the arc
    midpoints as 0.12929/0.27071 (exact 0.2 -+ R*sqrt(2) =
    0.129289321881345 / 0.270710678118655); only the NURBS weight
    0.707106781186548 is exact. Left uncorrected this puts a floor of about
    7e-6 (1.4e-4 relative) under every comparison against an immersed
    solver that uses the exact circle, a floor self-convergence alone can
    never reveal. Each arc control point P is therefore replaced by
    c + rad*(P-c)/|P-c| (rad = R at the arc ends, R*sqrt(2) at the
    midpoint) -- a projection that depends on the stored point only through
    its direction, so corners shared between neighbouring patches are
    recomputed identically from both sides and the interfaces stay
    conforming.

    Quantities of interest. c_D, c_L follow the Schaefer-Turek scaling
    c_{D,L} = F_{1,2} / (rho*Ubar^2*R), rho = 1, Ubar = 0.2, R = 0.05, so
    c_{D,L} = F_{1,2}/0.002. The drag/lift force is computed two ways:

    1. Reaction functional (volume form). With ell_i = -e_i*chi(|x-c|) a
       smooth cutoff (chi = 1 for r <= 0.075, 0 for r >= 0.15, a quintic
       blend in between, so grad(chi) is supported only in that annulus)
       and the discrete residual R(u,p;w) = 2*mu*grad^s(u):grad^s(w) -
       p*div(w) tested with w = ell_i:

         R(u,p;ell_i) = -integral_Omega e_i . sigma_h(u,p) . grad(chi) dx = F_i,

       using grad(ell_i) = -e_i (x) grad(chi) and div(ell_i) = -e_i.grad(chi),
       so 2*mu*grad^s(u):grad^s(ell_i) - p*div(ell_i)
        = -e_i . [mu*(grad(u)+grad(u)^T) - p*I] . grad(chi) = -e_i.sigma_h.grad(chi).
       Since sigma_h is exactly what -div(sigma_h) balances in the discrete
       weak form, this reproduces the boundary reaction on the cylinder
       without differentiating a discrete traction there.
    2. Direct surface integral: F_i = -integral_(cylinder) e_i . sigma_h . n_Omega,
       n_Omega = the outward normal of Omega (into the cylinder). The two
       are algebraically identical at the continuous level (integration by
       parts, div(sigma) = 0 in Omega, chi = 1 on the cylinder and
       grad(chi) = 0 near every other boundary) and are printed side by
       side as a discrete cross-check.

    Delta p = p(0.15, 0.2) - p(0.25, 0.2), the pressure at the upstream and
    downstream arc midpoints of the cylinder (both ON the boundary).

    Reference-file format (-o/--output; for read-back by an immersed solver).
    The finest solved level is written with gsFileData labels:

      "geometry" -- the snapped 5-patch NURBS multipatch, parameter domain
                    [0,1]^2 per patch.
      "velocity" -- 5-patch multipatch, targetDim 2; patch i lives on the
                    parameter domain of "geometry" patch i, i.e. the
                    physical value at x is velocity.patch(i).eval(xi) with
                    geometry.patch(i)(xi) = x. Includes the Dirichlet dofs
                    (gsFeSolution::extract fills them from the fixed part).
      "pressure" -- 5-patch multipatch, targetDim 1, same parametrization.
      "qoi"      -- 6x1 gsMatrix: [c_D_reac, c_L_reac, dp, c_D_surf,
                    c_L_surf, k_ref].

    A comment records k_ref, r, mu, Um and ndofs; a round-trip read-back
    checks patch counts and reproduces Delta p to 1e-12.

    Example command lines:
      ./stokes_cylinder_bodyfitted_example -k 3 -r 2
      ./stokes_cylinder_bodyfitted_example -k 3 -r 4 --study self
      ./stokes_cylinder_bodyfitted_example -k 3 -r 2 -o ref.xml

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>

#include <cmath>
#include <iomanip>
#include <limits>
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

std::string fmtFix(real_t v, int prec = 2)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec) << v;
    return os.str();
}

/// Gradient of the radial cutoff
///   chi(r) = 1                                  (r <= rIn)
///          = 0                                  (r >= rOut)
///          = 1 - (6*t^5 - 15*t^4 + 10*t^3),  t = (r-rIn)/(rOut-rIn)   (rIn<r<rOut)
/// used by the reaction functional: grad(chi)(x) = chi'(r)*(x-c)/r,
/// chi'(r) = -30*t^2*(1-t)^2/(rOut-rIn) on the annulus, 0 elsewhere. Only
/// eval_into is needed by the reaction functional (Ju is taken from the
/// discrete solution, not from this function); deriv_into is intentionally
/// omitted.
class gsCutoffGrad : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(gsCutoffGrad)

    gsCutoffGrad(gsVector<real_t,2> c, real_t rIn, real_t rOut)
    : m_c(c), m_rIn(rIn), m_rOut(rOut)
    { GISMO_ENSURE(m_rOut > m_rIn, "gsCutoffGrad: rOut must exceed rIn."); }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 2; }

    // The cutoff is a single function of the physical coordinate, defined
    // identically on every patch of the multipatch it is composed with
    // (gsFunction::piece() otherwise asserts k==0, since a plain gsFunction
    // is assumed single-subdomain; gsFunctionExpr/gsConstantFunction use the
    // same override for the same reason).
    const gsCutoffGrad & piece(const index_t) const override { return *this; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(2, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
        {
            const gsVector<real_t,2> d = u.col(k) - m_c;
            const real_t r = d.norm();
            if (r <= m_rIn || r >= m_rOut) { result.col(k).setZero(); continue; }
            const real_t t = (r - m_rIn) / (m_rOut - m_rIn);
            const real_t chiPrime = -30.0 * t * t * (1.0 - t) * (1.0 - t) / (m_rOut - m_rIn);
            result.col(k) = (chiPrime / r) * d;
        }
    }

    std::ostream & print(std::ostream & os) const override
    { return os << "gsCutoffGrad(c=(" << m_c[0] << "," << m_c[1] << "), rIn=" << m_rIn << ", rOut=" << m_rOut << ")"; }

private:
    gsVector<real_t,2> m_c;
    real_t m_rIn, m_rOut;
};

/// Compares gsCutoffGrad::eval_into against a central finite difference
/// (step 1e-6) of chi at points sampled across the annulus, and checks
/// that the gradient is exactly zero just outside it (r = 0.07 < rIn and
/// r = 0.16 > rOut). Guards against a sign or exponent typo in the closed
/// form silently degrading the reaction QoI instead of failing visibly.
bool checkCutoffGrad(const gsCutoffGrad & gchi, const gsVector<real_t,2> & c,
                      real_t rIn, real_t rOut)
{
    auto chi = [&](real_t x, real_t y) -> real_t
    {
        const real_t r = std::sqrt((x - c[0]) * (x - c[0]) + (y - c[1]) * (y - c[1]));
        if (r <= rIn)  return 1.0;
        if (r >= rOut) return 0.0;
        const real_t t = (r - rIn) / (rOut - rIn);
        return 1.0 - (6.0 * t * t * t * t * t - 15.0 * t * t * t * t + 10.0 * t * t * t);
    };

    const real_t h = 1e-6;
    real_t worst = 0.0;
    const real_t radii[4] = { 0.08, 0.10, 0.12, 0.14 };
    for (real_t rr : radii)
        for (index_t a = 0; a != 8; ++a)
        {
            const real_t theta = (real_t)a / 8.0 * 2.0 * EIGEN_PI;
            const real_t x = c[0] + rr * std::cos(theta), y = c[1] + rr * std::sin(theta);
            gsMatrix<real_t> pt(2,1); pt << x, y;
            gsMatrix<real_t> val; gchi.eval_into(pt, val);
            const real_t fdx = (chi(x + h, y) - chi(x - h, y)) / (2.0 * h);
            const real_t fdy = (chi(x, y + h) - chi(x, y - h)) / (2.0 * h);
            const real_t relx = std::abs(val(0,0) - fdx) / std::max(std::abs(fdx), 1.0);
            const real_t rely = std::abs(val(1,0) - fdy) / std::max(std::abs(fdy), 1.0);
            worst = std::max(worst, std::max(relx, rely));
        }

    gsMatrix<real_t> ptIn(2,1), ptOut(2,1), valIn, valOut;
    ptIn  << c[0] + 0.07, c[1];
    ptOut << c[0] + 0.16, c[1];
    gchi.eval_into(ptIn, valIn);
    gchi.eval_into(ptOut, valOut);
    const bool zerosOk = valIn.isZero(0.0) && valOut.isZero(0.0);

    gsInfo << "gsCutoffGrad self-check: worst relative FD error = " << fmtSci(worst)
           << ", exactly zero outside the annulus: " << (zerosOk ? "yes" : "no") << "\n";
    GISMO_ENSURE(worst <= 1e-6, "gsCutoffGrad gradient disagrees with a finite difference: "
                 << worst << " > 1e-6.");
    GISMO_ENSURE(zerosOk, "gsCutoffGrad is not exactly zero at r=0.07 or r=0.16.");
    return true;
}

} // anonymous namespace

/// Result of one Taylor-Hood Stokes solve: dof counts, the two drag/lift
/// pairs, Delta p, sanity quantities, and -- when requested -- the
/// extracted velocity/pressure multipatches for the -o reference file.
struct QoIResult
{
    index_t ndofs = 0;
    real_t cD_reac = 0, cL_reac = 0, cD_surf = 0, cL_surf = 0, dp = 0;
    real_t divInt = 0, divL2 = 0, fluxIn = 0, fluxOut = 0;
    index_t pid0 = -1, pid1 = -1;
    gsVector<real_t> preim0, preim1;
    gsMultiPatch<real_t> velocity, pressure;
};

/// Assembles and solves the Taylor-Hood Stokes system on basisV/basisP
/// (pressure degree kref, velocity kref+1, r uniform refinements of mp),
/// then evaluates the QoIs described in the file header. When
/// \a printChecks is set, also prints and GISMO_ENSUREs the domain-area and
/// cylinder-length identities on the refined velocity basis built here --
/// run once, on the first solve of a study, since the checks are
/// refinement-independent up to quadrature accuracy.
QoIResult solveStokes(const gsMultiPatch<real_t> & mp, index_t kref, index_t r, real_t quA,
                       const gsCutoffGrad & cutoffGrad, bool extract, bool printChecks)
{
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::space       space;
    typedef gsExprAssembler<real_t>::solution    solution;

    GISMO_ENSURE(kref >= 2, "k_ref = " << kref << " < 2 would need degreeReduce on a "
                 "rational pressure basis (setDegree calls degreeReduce below the stored "
                 "degree 2); use --kref >= 2.");

    const real_t mu = 1e-3, Um = 0.3, H = 0.41;
    const gsVector<real_t,2> cylCenter = gsVector<real_t,2>::vec(0.2, 0.2);
    const real_t cylRadius = 0.05;
    const real_t rho = 1.0, Ubar = 0.2;
    const real_t scale = rho * Ubar * Ubar * cylRadius;   // = 0.002

    gsMultiBasis<real_t> basisP(mp);
    basisP.setDegree((short_t)kref);
    for (index_t i = 0; i != r; ++i) basisP.uniformRefine();
    gsMultiBasis<real_t> basisV = basisP;
    basisV.degreeElevate(1);

    std::ostringstream inflowExpr;
    inflowExpr << "4*" << Um << "*y*(" << H << "-y)/" << H << "^2";
    gsFunctionExpr<real_t> inflow(inflowExpr.str(), "0", 2);
    gsConstantFunction<real_t> zero2(0., 0., 2);

    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::north, condition_type::dirichlet, inflow, 0, false, -1);
    bc.addCondition(1, boundary::north, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(4, boundary::south, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(3, boundary::north, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(4, boundary::north, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(0, boundary::south, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(1, boundary::south, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(2, boundary::south, condition_type::dirichlet, zero2,  0, false, -1);
    bc.addCondition(3, boundary::south, condition_type::dirichlet, zero2,  0, false, -1);
    bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(2, 2);
    A.options().setReal("quA", quA);

    geometryMap G = A.getMap(mp);
    space v = A.getSpace(basisV, 2, 0);
    space p = A.getSpace(basisP, 1, 1);
    v.setup(bc, dirichlet::l2Projection, 0);   // C0 velocity across interfaces
    p.setup(bc, dirichlet::l2Projection, 0);   // no conditions reach unknown 1; C0 pressure
    A.setIntegrationElements(basisV);
    A.initSystem();

    gsMatrix<real_t> solVector;   // NAMED: getSolution keeps a pointer to it
    solution v_sol = A.getSolution(v, solVector);
    solution p_sol = A.getSolution(p, solVector);

    auto Jv = ijac(v, G);

    gsStopwatch time;
    A.assemble(
          mu * ((Jv.cwisetr() + Jv) % Jv.tr()) * meas(G)
        , -idiv(v, G) * p.tr() * meas(G)
        , -p * idiv(v, G).tr() * meas(G)
    );
    const real_t tAssemble = time.stop();

    time.restart();
    gsSparseSolver<real_t>::LU solver;
    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());
    const real_t tSolve = time.stop();
    GISMO_ENSURE(solVector.allFinite(), "Solve produced a non-finite solution vector.");

    QoIResult res;
    res.ndofs = A.numDofs();
    gsInfo << "  [kref=" << kref << " r=" << r << "] ndofs=" << res.ndofs
           << "  assembly=" << fmtFix(tAssemble,3) << "s  solve=" << fmtFix(tSolve,3) << "s\n";

    gsExprEvaluator<real_t> ev(A);
    ev.options().setReal("quA", quA);

    std::vector<patchSide> cylSides = { patchSide(0, boundary::south), patchSide(1, boundary::south),
                                         patchSide(2, boundary::south), patchSide(3, boundary::south) };

    if (printChecks)
    {
        gsInfo << "Assembler options (quA is used identically by the evaluator below, so the "
                  "div-u/flux identity in the sanity block is exact up to the assembler's own "
                  "quadrature, not merely to the tolerance below):\n" << A.options() << "\n";

        const real_t area = ev.integral(meas(G));
        const real_t areaExact = 0.902 - 0.0025 * EIGEN_PI;
        const real_t areaRelErr = std::abs(area - areaExact) / areaExact;
        const real_t length = ev.integralBdr(nv(G).norm(), cylSides);
        const real_t lengthExact = 0.1 * EIGEN_PI;
        const real_t lengthRelErr = std::abs(length - lengthExact) / lengthExact;

        gsInfo << "Startup check -- domain area = " << area << " (exact " << areaExact
               << ", rel.err. " << fmtSci(areaRelErr) << ")\n";
        gsInfo << "Startup check -- cylinder length = " << length << " (exact " << lengthExact
               << ", rel.err. " << fmtSci(lengthRelErr) << ")\n";
        GISMO_ENSURE(areaRelErr <= 1e-8, "Domain area relative error " << areaRelErr << " > 1e-8.");
        GISMO_ENSURE(lengthRelErr <= 1e-8, "Cylinder length relative error " << lengthRelErr << " > 1e-8.");
    }

    // Sanity: discrete continuity and inflow/outflow mass flux.
    res.divInt = ev.integral(idiv(v_sol, G) * meas(G));
    res.divL2  = std::sqrt(ev.integral(idiv(v_sol, G).sqNorm() * meas(G)));

    std::vector<patchSide> inflowSide  = { patchSide(0, boundary::north) };
    std::vector<patchSide> outflowSide = { patchSide(4, boundary::east) };
    res.fluxIn  = ev.integralBdr(v_sol.tr() * nv(G), inflowSide);
    res.fluxOut = ev.integralBdr(v_sol.tr() * nv(G), outflowSide);

    // Drag/lift: reaction functional and direct surface integral.
    gsConstantFunction<real_t> e1(1., 0., 2), e2(0., 1., 2);
    auto gchi = ev.getVariable(cutoffGrad, G);
    auto E1   = ev.getVariable(e1, G);
    auto E2   = ev.getVariable(e2, G);
    auto Ju   = ijac(v_sol, G);   // 2x2 gradient of the vector solution

    const real_t F1 = -ev.integral( E1.tr() * ( mu * (Ju + Ju.tr()) * gchi - p_sol.val() * gchi ) * meas(G) );
    const real_t F2 = -ev.integral( E2.tr() * ( mu * (Ju + Ju.tr()) * gchi - p_sol.val() * gchi ) * meas(G) );

    const real_t F1s = -ev.integralBdr( E1.tr() * ( mu * (Ju + Ju.tr()) * nv(G) - p_sol.val() * nv(G) ), cylSides );
    const real_t F2s = -ev.integralBdr( E2.tr() * ( mu * (Ju + Ju.tr()) * nv(G) - p_sol.val() * nv(G) ), cylSides );

    GISMO_ENSURE(F1 > 0, "Drag (reaction functional) is not positive: F1 = " << F1);

    res.cD_reac = F1  / scale;  res.cL_reac = F2  / scale;
    res.cD_surf = F1s / scale;  res.cL_surf = F2s / scale;

    // Delta p at the upstream/downstream arc midpoints, both ON the boundary.
    gsMatrix<real_t> dpPts(2, 2);
    dpPts.col(0) << 0.15, 0.2;
    dpPts.col(1) << 0.25, 0.2;
    gsVector<index_t> pids;
    gsMatrix<real_t> preim;
    mp.locatePoints(dpPts, pids, preim, (real_t)1e-12);
    for (index_t i = 0; i != 2; ++i)
    {
        if (-1 == pids[i])
        {
            pids(i) = (0 == i) ? 0 : 2;
            preim(0, i) = 0.5; preim(1, i) = 0.0;
            gsInfo << "  Delta-p point " << i << ": locatePoints returned -1 (boundary Newton "
                      "landed just outside the parameter box); using the known fallback "
                      "preimage (patch " << pids[i] << ", (0.5,0)).\n";
        }
        gsMatrix<real_t> chk;
        mp.patch(pids[i]).eval_into(preim.col(i), chk);
        const real_t err = (chk.col(0) - dpPts.col(i)).norm();
        GISMO_ENSURE(err <= 1e-10, "Delta-p preimage " << i << " does not map back to the "
                     "target point: " << err << " > 1e-10.");
        gsInfo << "  Delta-p point " << i << ": patch " << pids[i] << ", preimage = ("
               << preim(0,i) << ", " << preim(1,i) << ")\n";
    }
    gsVector<real_t> pt0 = preim.col(0), pt1 = preim.col(1);
    const real_t pA = ev.eval(p_sol, pt0, pids[0])(0, 0);
    const real_t pB = ev.eval(p_sol, pt1, pids[1])(0, 0);
    res.dp = pA - pB;
    res.pid0 = pids[0]; res.pid1 = pids[1];
    res.preim0 = pt0;   res.preim1 = pt1;

    if (extract)
    {
        v_sol.extract(res.velocity);
        p_sol.extract(res.pressure);
    }

    return res;
}

/// Writes the reference XML documented in the file header and immediately
/// reads it back, checking patch counts and that Delta p at the (pid,
/// preim) pair recorded in \a res reproduces res.dp to 1e-12 -- the
/// same read-back a consumer of the file performs.
void writeReference(const gsMultiPatch<real_t> & mp, const QoIResult & res,
                     index_t kref, index_t r, const std::string & fn)
{
    gsMatrix<real_t> qoi(6, 1);
    qoi << res.cD_reac, res.cL_reac, res.dp, res.cD_surf, res.cL_surf, (real_t)kref;

    gsFileData<real_t> fd;
    fd.addWithLabel(mp,           "geometry");
    fd.addWithLabel(res.velocity, "velocity");
    fd.addWithLabel(res.pressure, "pressure");
    fd.addWithLabel(qoi,          "qoi");
    std::ostringstream note;
    note << "stokes_cylinder_bodyfitted_example reference: kref=" << kref << " r=" << r
         << " mu=1e-3 Um=0.3 ndofs=" << res.ndofs;
    fd.addComment(note.str());
    fd.save(fn);
    gsInfo << "Wrote reference file '" << fn << "' (kref=" << kref << ", r=" << r << ").\n";

    gsFileData<real_t> fd2(fn);
    gsMultiPatch<real_t> geomChk, velChk, presChk;
    fd2.getLabel("geometry", geomChk);
    fd2.getLabel("velocity", velChk);
    fd2.getLabel("pressure", presChk);
    GISMO_ENSURE(5 == geomChk.nPatches() && 5 == velChk.nPatches() && 5 == presChk.nPatches(),
                 "Round-trip patch count mismatch: geometry=" << geomChk.nPatches()
                 << " velocity=" << velChk.nPatches() << " pressure=" << presChk.nPatches());

    gsMatrix<real_t> pA, pB;
    presChk.patch(res.pid0).eval_into(res.preim0, pA);
    presChk.patch(res.pid1).eval_into(res.preim1, pB);
    const real_t dpChk = pA(0, 0) - pB(0, 0);
    const real_t err = std::abs(dpChk - res.dp);
    gsInfo << "Round-trip check: dp reproduced to " << fmtSci(err) << " (dp=" << res.dp << ").\n";
    GISMO_ENSURE(err <= 1e-12, "Round-trip dp mismatch: " << err << " > 1e-12.");
}

int main(int argc, char * argv[])
{
    std::string fn = "planar/flow_around_cylinder.xml";
    index_t kref = 3;
    index_t numRefine = 2;
    std::string study = "single";
    std::string output = "";
    real_t quA = 2.0;

    gsCmdLine cmd("Body-fitted Taylor-Hood Stokes reference on the Schaefer-Turek cylinder geometry.");
    cmd.addString("f", "file",   "Input geometry file (gsFileData, relative to filedata/)", fn);
    cmd.addInt   ("k", "kref",   "Pressure degree (velocity degree = kref+1); >= 2", kref);
    cmd.addInt   ("r", "refine", "Uniform refinement level (self-convergence upper level R)", numRefine);
    cmd.addString("",  "study",  "single | self", study);
    cmd.addString("o", "output", "Write geometry+velocity+pressure+QoIs to this XML file", output);
    cmd.addReal  ("",  "quA",    "Evaluator/assembler quadrature parameter quA*deg+quB", quA);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE("single" == study || "self" == study, "--study must be 'single' or 'self', got '" << study << "'.");

    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>(fn, mp);
    GISMO_ENSURE(5 == mp.nPatches() && 5 == mp.nInterfaces() && 10 == mp.nBoundary(),
                 "Unexpected topology: nPatches=" << mp.nPatches() << " nInterfaces=" << mp.nInterfaces()
                 << " nBoundary=" << mp.nBoundary() << " (expected 5/5/10).");

    const gsVector<real_t,2> cylCenter = gsVector<real_t,2>::vec(0.2, 0.2);
    const real_t cylRadius = 0.05;

    // Snap the four arc patches (0-3) onto the exact circle; see file header.
    for (index_t p = 0; p != 4; ++p)
    {
        gsMatrix<real_t> & coefs = mp.patch(p).coefs();
        for (index_t i = 0; i != 3; ++i)
        {
            const gsVector<real_t,2> P = coefs.row(i).transpose();
            const gsVector<real_t,2> d = (P - cylCenter).normalized();
            const real_t rad = (1 == i) ? std::sqrt(2.0) * cylRadius : cylRadius;
            coefs.row(i) = (cylCenter + rad * d).transpose();
        }
    }

    // Machine check: the snapped arcs lie on the exact circle to 1e-13.
    {
        gsMatrix<real_t> ptsUV(2, 201);
        for (index_t j = 0; j != 201; ++j) { ptsUV(0,j) = (real_t)j / 200.0; ptsUV(1,j) = 0.0; }
        real_t worstCircle = 0.0;
        for (index_t p = 0; p != 4; ++p)
        {
            gsMatrix<real_t> vals;
            mp.patch(p).eval_into(ptsUV, vals);
            for (index_t j = 0; j != 201; ++j)
                worstCircle = std::max(worstCircle, std::abs((vals.col(j) - cylCenter).norm() - cylRadius));
        }
        gsInfo << "Circle deviation (max over patches 0-3, 201 pts each): " << fmtSci(worstCircle) << "\n";
        GISMO_ENSURE(worstCircle <= 1e-13, "Snapped arcs deviate from the exact circle by "
                     << worstCircle << " > 1e-13.");
    }

    const real_t rIn = 0.075, rOut = 0.15;
    gsCutoffGrad cutoffGrad(cylCenter, rIn, rOut);
    checkCutoffGrad(cutoffGrad, cylCenter, rIn, rOut);

    if ("single" == study)
    {
        QoIResult res = solveStokes(mp, kref, numRefine, quA, cutoffGrad, !output.empty(), /*printChecks*/ true);

        gsInfo << "\n=== single run  kref=" << kref << "  r=" << numRefine << "  ndofs=" << res.ndofs << " ===\n";
        gsInfo << "c_D (reaction)  = " << fmtSci(res.cD_reac,10) << "\n";
        gsInfo << "c_D (surface)   = " << fmtSci(res.cD_surf,10) << "\n";
        gsInfo << "c_L (reaction)  = " << fmtSci(res.cL_reac,10) << "\n";
        gsInfo << "c_L (surface)   = " << fmtSci(res.cL_surf,10) << "\n";
        gsInfo << "Delta p         = " << fmtSci(res.dp,10) << "\n";
        gsInfo << "integral(div u) = " << fmtSci(res.divInt) << ",  ||div u||_L2 = " << fmtSci(res.divL2) << "\n";
        gsInfo << "flux_in = " << fmtSci(res.fluxIn,10) << " (|flux_in-(-0.082)| = " << fmtSci(std::abs(res.fluxIn+0.082))
               << "),  flux_out = " << fmtSci(res.fluxOut,10)
               << ",  flux_in+flux_out = " << fmtSci(res.fluxIn + res.fluxOut) << "\n";

        if (!output.empty()) writeReference(mp, res, kref, numRefine, output);
    }
    else // self
    {
        gsInfo << "\n=== self-convergence study  kref=" << kref << "  r=0.." << numRefine << " ===\n";
        gsInfo << std::right
               << std::setw(4)  << "r"     << std::setw(8)  << "ndofs"
               << std::setw(15) << "cD_reac" << std::setw(15) << "cD_surf"
               << std::setw(15) << "cL_reac" << std::setw(15) << "cL_surf"
               << std::setw(15) << "dp"
               << std::setw(13) << "|dcD|" << std::setw(13) << "|dcL|" << std::setw(13) << "|ddp|"
               << std::setw(9)  << "ratio" << std::setw(15) << "|cD_r-cD_s|" << "\n";

        real_t prevCD = 0, prevCL = 0, prevDp = 0;
        std::vector<real_t> diffsCD;
        QoIResult last;
        for (index_t rr = 0; rr <= numRefine; ++rr)
        {
            QoIResult res = solveStokes(mp, kref, rr, quA, cutoffGrad,
                                         !output.empty() && rr == numRefine, rr == 0);
            last = res;

            const bool havePrev = (rr > 0);
            const real_t dCD = havePrev ? std::abs(res.cD_reac - prevCD) : std::numeric_limits<real_t>::quiet_NaN();
            const real_t dCL = havePrev ? std::abs(res.cL_reac - prevCL) : std::numeric_limits<real_t>::quiet_NaN();
            const real_t dDp = havePrev ? std::abs(res.dp      - prevDp) : std::numeric_limits<real_t>::quiet_NaN();
            if (havePrev) diffsCD.push_back(dCD);
            const real_t ratio = (diffsCD.size() >= 2)
                ? diffsCD[diffsCD.size()-2] / diffsCD.back()
                : std::numeric_limits<real_t>::quiet_NaN();
            const real_t crossCheck = std::abs(res.cD_reac - res.cD_surf);

            gsInfo << std::setw(4) << rr << std::setw(8) << res.ndofs
                   << std::setw(15) << fmtSci(res.cD_reac) << std::setw(15) << fmtSci(res.cD_surf)
                   << std::setw(15) << fmtSci(res.cL_reac) << std::setw(15) << fmtSci(res.cL_surf)
                   << std::setw(15) << fmtSci(res.dp)
                   << std::setw(13) << fmtSci(dCD,2) << std::setw(13) << fmtSci(dCL,2) << std::setw(13) << fmtSci(dDp,2)
                   << std::setw(9)  << fmtFix(ratio) << std::setw(15) << fmtSci(crossCheck) << "\n";

            prevCD = res.cD_reac; prevCL = res.cL_reac; prevDp = res.dp;
        }

        if (!output.empty()) writeReference(mp, last, kref, numRefine, output);
    }

    return 0;
}
