/** @file rh-adaptive_fitting_example.cpp

    @brief rh-adaptive fitting of a parametrized point cloud with a composed
    geometry G = S ∘ σ.

    The point cloud (uv, xyz) is fitted by a THB-spline geometry S whose
    parametrization is composed with a square-domain map σ (gsSquareDomain):
    the data parameters uv live in the σ-domain Ω̃ = [0,1]^2 and the fitted
    map is G(uv) = S(σ(uv)). Three building blocks can be combined in any
    order through the --schedule string:

      F  (Fit)      penalized least-squares fit of the S coefficients with
                    the parameters fixed at ξ_i = σ(uv_i) — equivalent to
                    collocation of the composed basis at uv_i.
      R  (r-adapt)  composite spline relocation of σ via
                    gsAdaptiveParametrization<ValueBased> + HLBFGS, driven by
                    a monitor spline least-squares fitted to the (normalized)
                    pointwise fitting errors over the S-domain Ω̂ (parametric
                    monitor, ω = 1/sqrt(1+θ f)).
      D  (direct)   direct minimization of the true LS fitting error
                    Σ_i ||S(σ(uv_i;α)) - x_i||² over the σ controls α
                    (gsOptProblem + HLBFGS, analytic gradient, fold barrier
                    on det J_σ) — variable projection when alternated with F.
      H  (h-adapt)  THB refinement of the cells (in Ω̂) containing data
                    points whose error exceeds the threshold (or the
                    (1-refPercent) error percentile), with cell extension —
                    same marking as gsHFitting, but applied at ξ_i = σ(uv_i).

    The schedule string is repeated up to -i times, e.g.
      --schedule FRH   :  [Fit → r-adapt → h-refine] per cycle
      --schedule FHR   :  [Fit → h-refine → r-adapt] per cycle
      --schedule FRFRH :  two Fit/r rounds before each h-refinement
      --schedule F     :  plain fitting, no adaptivity

    Initialization: --sigmaInit radial replaces the identity start of σ by a
    prescribed C² radial map that EXPANDS the parametrization on the ring
    r = --sigmaR0 around the domain centre and is the identity for
    r >= --sigmaTaper — a hand-made "answer" for a feature known to sit on
    that ring, useful to separate "the optimizer cannot find a good σ" from
    "no good σ exists in this space" (see gsRadialSigma below).

    Expected input (same as fitting_example.cpp): an XML file containing
      Matrix id 0 : 2 x N parameter values (rescaled here to [0,1]^2)
      Matrix id 1 : 3 x N (or 2 x N) point coordinates

    NOTE on cost: the r-adaptivity energy is integrated on the finest tensor
    level of the THB basis (gsAdaptiveParametrization requires a tensor
    integration basis and internally raises the degree to p_S * p_σ), so R
    steps become more expensive as h-refinement proceeds.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
// Direct (D) step objective over the sigma controls; shared with the other
// composed-fitting drivers.
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsHLBFGS/gsHLBFGS.h>

#include <fstream>
#include <iomanip>

using namespace gismo;

namespace {

// Thin-plate smoothing matrix, adapted from gsFitting<T>::applySmoothing
// (src/gsModeling/gsFitting.hpp): adds lambda * int (sum of squared second
// derivatives) to A.
void applySmoothing(const gsBasis<real_t> & basis, real_t lambda,
                    gsSparseMatrix<real_t> & A)
{
    const short_t dim    = basis.domainDim();
    const short_t stride = dim * (dim + 1) / 2;

    gsVector<index_t> numNodes(dim);
    for (short_t i = 0; i != dim; ++i)
        numNodes[i] = basis.degree(i);
    gsGaussRule<real_t> QuRule(numNodes);

    gsMatrix<real_t> quNodes, der2, localA;
    gsVector<real_t> quWeights;
    gsMatrix<index_t> actives;

    typename gsBasis<real_t>::domainIter domIt    = basis.domain()->beginAll();
    typename gsBasis<real_t>::domainIter domItEnd = basis.domain()->endAll();
    for (; domIt < domItEnd; ++domIt)
    {
        QuRule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights);
        basis.deriv2_into(quNodes, der2);
        basis.active_into(domIt.centerPoint(), actives);
        const index_t numActive = actives.rows();
        localA.setZero(numActive, numActive);

        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            const real_t weight = quWeights[k] * lambda;
            for (index_t i = 0; i != numActive; ++i)
                for (index_t j = 0; j != numActive; ++j)
                {
                    real_t localAij = 0;
                    for (short_t s = 0; s < stride; s++)
                    {
                        // pure second derivatives once, mixed ones twice
                        const real_t dij = der2(i * stride + s, k) *
                                           der2(j * stride + s, k);
                        localAij += (s < dim) ? dij : 2 * dij;
                    }
                    localA(i, j) += weight * localAij;
                }
        }

        for (index_t i = 0; i != numActive; ++i)
            for (index_t j = 0; j != numActive; ++j)
                A.coeffRef(actives(i, 0), actives(j, 0)) += localA(i, j);
    }
}

// Penalized least-squares fit of values (d x N, at params 2 x N in the basis
// domain) returning the coefficient matrix (basis.size() x d).
gsMatrix<real_t> leastSquaresFit(const gsBasis<real_t> & basis,
                                 const gsMatrix<real_t> & params,
                                 const gsMatrix<real_t> & values,
                                 real_t lambda)
{
    gsSparseMatrix<real_t> C = basis.collocationMatrix(params);
    gsSparseMatrix<real_t> A = C.transpose() * C;
    if (lambda > 0)
        applySmoothing(basis, lambda, A);
    // Tiny Tikhonov shift guards basis functions without data support
    for (index_t i = 0; i != A.rows(); ++i)
        A.coeffRef(i, i) += 1e-12;
    A.makeCompressed();
    gsMatrix<real_t> b = C.transpose() * values.transpose();

    gsSparseSolver<real_t>::SimplicialLDLT solver(A);
    GISMO_ENSURE(solver.info() == gsEigen::Success,
                 "Least-squares factorization failed");
    return solver.solve(b);
}

bool isCellAlreadyInserted(const gsVector<index_t, 2> & a_cell,
                           const std::vector<index_t> & cells)
{
    for (size_t i = 0; i != cells.size(); i += a_cell.rows())
    {
        index_t commonEntries = 0;
        for (index_t col = 0; col != a_cell.rows(); col++)
            if (cells[i + col] == a_cell[col])
                commonEntries++;
        if (commonEntries == a_cell.rows())
            return true;
    }
    return false;
}

// Appends the refinement box around one marked parameter (adapted from
// gsHFitting<2,T>::appendBox); param lives in the S-domain Ω̂.
void appendBox(std::vector<index_t> & boxes,
               std::vector<index_t> & cells,
               const gsTHBSplineBasis<2> & basis,
               const gsVector<real_t> & param,
               const std::vector<unsigned> & ext)
{
    const int maxLvl = basis.maxLevel();
    const gsTensorBSplineBasis<2> & tBasis = *(basis.getBases()[maxLvl]);

    gsVector<index_t, 2> a_cell;
    for (short_t dim = 0; dim != 2; dim++)
    {
        const gsKnotVector<real_t> & kv = tBasis.component(dim).knots();
        a_cell(dim) = kv.uFind(param(dim)).uIndex();
    }

    if (isCellAlreadyInserted(a_cell, cells))
        return;
    for (index_t col = 0; col != a_cell.rows(); col++)
        cells.push_back(a_cell[col]);

    gsVector<index_t, 2> a_cell_upp = a_cell + gsVector<index_t, 2>::Ones();
    const int cell_lvl = basis.tree().query3(a_cell, a_cell_upp, maxLvl) + 1;

    gsVector<index_t> box(5);
    box[0] = cell_lvl;
    for (short_t dim = 0; dim != 2; dim++)
    {
        const unsigned numBreaks = basis.numBreaks(cell_lvl, dim) - 1;
        unsigned lowIndex = 0;
        if (cell_lvl < maxLvl)
            lowIndex = (a_cell(dim) >> (maxLvl - cell_lvl));
        else
            lowIndex = (a_cell(dim) << (cell_lvl - maxLvl));

        box[1 + dim] = ((lowIndex > ext[dim]) ? (lowIndex - ext[dim]) : 0);
        box[3 + dim] = ((lowIndex + ext[dim] + 1 < numBreaks) ?
                        (lowIndex + ext[dim] + 1) : numBreaks);
    }
    for (index_t col = 0; col != box.rows(); col++)
        boxes.push_back(box[col]);
}

// Refinement boxes for all parameters whose error exceeds the threshold
// (adapted from gsHFitting<2,T>::getBoxes).
std::vector<index_t> getBoxes(const gsTHBSplineBasis<2> & basis,
                              const std::vector<real_t> & errors,
                              const gsMatrix<real_t> & params,
                              real_t threshold,
                              const std::vector<unsigned> & ext)
{
    std::vector<index_t> cells, boxes;
    for (size_t i = 0; i != errors.size(); i++)
        if (threshold <= errors[i])
            appendBox(boxes, cells, basis, params.col(i), ext);
    return boxes;
}

// Error value such that refPercent of the points lie above it (adapted from
// gsHFitting<2,T>::setRefineThreshold).
real_t refineThreshold(const std::vector<real_t> & errors, real_t refPercent)
{
    std::vector<real_t> errorsCopy = errors;
    const size_t i = cast<real_t, size_t>(errorsCopy.size() * (1.0 - refPercent));
    typename std::vector<real_t>::iterator pos = errorsCopy.begin() + i;
    std::nth_element(errorsCopy.begin(), pos, errorsCopy.end());
    return *pos;
}

// Distance (in the S-domain) to the step curve of the lshape_step dataset,
// i.e. the boundary of the cut quadrant [0.5,1]^2 inside (0,1)^2.
real_t distToStep(real_t x, real_t y)
{
    const real_t a = 0.5 - x, b = 0.5 - y;
    if (a <= 0 && b <= 0)                      // inside the cut quadrant
        return std::min(-a, -b);
    return math::sqrt(std::max(a, 0.0) * std::max(a, 0.0) +
                      std::max(b, 0.0) * std::max(b, 0.0));
}

// Analytic band monitor f(xi) = exp(-dist(xi, step)^2 / w^2) for the
// lshape_step dataset. Bypasses the LS monitor fit entirely: the monitor
// amplitude is exactly 1 and the band width w is a free parameter, giving a
// clean test of the shear-stiffness hypothesis (H1: the Winslow integrand
// caps the det J contrast for thin bands regardless of monitor amplitude).
class gsAnalyticBandMonitor : public gsFunction<real_t>
{
public:
    gsAnalyticBandMonitor(real_t width, bool invert)
    : m_w(width), m_invert(invert) {}

    GISMO_CLONE_FUNCTION(gsAnalyticBandMonitor)

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    gsMatrix<real_t> support() const override
    {
        gsMatrix<real_t> res(2, 2);
        res << 0, 1, 0, 1;
        return res;
    }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t i = 0; i != u.cols(); i++)
        {
            const real_t d = distToStep(u(0, i), u(1, i));
            const real_t f = math::exp(-(d * d) / (m_w * m_w));
            result(0, i) = m_invert ? 1.0 - f : f;
        }
    }

private:
    real_t m_w;
    bool m_invert;
};

// Minimum of det J_sigma, sampled on a 7x7 point grid PER ELEMENT of sigma's
// knot mesh (one fixed global grid undersamples fine sigma meshes, exactly as
// the D-step barrier did before it moved to a Gauss rule on the same mesh).
// Shared by the D step and by the --sigmaInit radial startup check.
real_t minDetJSigma(const gsSquareDomain<real_t> & domain)
{
    const gsBasis<real_t> & sb = domain.domain().basis();
    gsVector<unsigned> np(2);
    np << 7, 7;
    gsMatrix<real_t> pts, dsig;
    real_t mn = std::numeric_limits<real_t>::max();
    for (auto & elem : sb.domain()->allElements())
    {
        pts = gsPointGrid(elem.lowerCorner(), elem.upperCorner(), np);
        domain.deriv_into(pts, dsig);
        for (index_t q = 0; q != pts.cols(); q++)
            mn = std::min(mn, dsig(0, q) * dsig(3, q) - dsig(1, q) * dsig(2, q));
    }
    return mn;
}

// Prescribed radial grading of the sigma map (--sigmaInit radial):
//
//   sigma_a(xi) = c0 + (g(r)/r) (xi - c0),   r = |xi - c0|,  c0 = (0.5, 0.5)
//
// The radial profile g is built from its DERIVATIVE as a zero-mean pair of C^2
// bumps B(u) = (1-u^2)^3 (support |u| < 1), integrated analytically:
//
//   g'(r) = 1 + (M-1) B((r-R0)/w) - ccomp B((r-rc)/Wc),  rc = Wc = taper/2,
//   ccomp = (M-1) w / Wc,
//   g(r)  = r + (M-1) w  IB((r-R0)/w) - ccomp Wc IB((r-rc)/Wc),
//   IB(u) = u - u^3 + (3/5) u^5 - (1/7) u^7   (antiderivative of B, IB(1)=16/35)
//
// DIRECTION (do not invert): the grade M > 1 makes sigma EXPAND on the ring
// r = R0, det J_sigma = g'(r) g(r)/r ~ 3 there for the defaults. The basis of S
// is evaluated at xi = sigma(uv), so a data feature sitting on the ring is
// mapped onto a THICKER annulus of the S-domain and therefore catches more of
// S's uniform knot lines. This mirrors the r-adaptive optimum measured on the
// lshape run (det J_sigma = 2.42 on the step, 0.54 away from it).
//
// The compensation bump returns the stretched length over [0, taper], so
// g(taper) = taper exactly and sigma_a is the identity for r >= taper; with
// taper < 0.5 that covers the whole boundary of [0,1]^2. g is monotone (no
// fold) iff ccomp < 1, i.e. (M-1) w < taper/2, and then min g' = 1 - ccomp,
// attained at r = rc.
class gsRadialSigma : public gsFunction<real_t>
{
public:
    gsRadialSigma(real_t R0, real_t taper, real_t grade, real_t width)
    : m_R0(R0), m_taper(taper), m_M(grade), m_w(width) {}

    GISMO_CLONE_FUNCTION(gsRadialSigma)

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 2; }

    gsMatrix<real_t> support() const override
    {
        gsMatrix<real_t> res(2, 2);
        res << 0, 1, 0, 1;
        return res;
    }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(2, u.cols());
        for (index_t i = 0; i != u.cols(); i++)
        {
            const real_t dx = u(0, i) - 0.5, dy = u(1, i) - 0.5;
            const real_t r  = math::sqrt(dx * dx + dy * dy);
            // identity outside the taper disk; the centre is a removable
            // singularity of g(r)/r (g'(0) = 1, g(r) - r = O(r^4))
            const real_t s = (r < 1e-12 || r >= m_taper) ? 1.0 : profile(r) / r;
            result(0, i) = 0.5 + s * dx;
            result(1, i) = 0.5 + s * dy;
        }
    }

    /// Radial profile g(r); g(0) = 0 and g(taper) = taper by construction
    real_t profile(real_t r) const
    {
        const real_t Wc    = 0.5 * m_taper;
        const real_t ccomp = (m_M - 1.0) * m_w / Wc;
        return r + (m_M - 1.0) * m_w * bumpInt((r - m_R0) / m_w)
                 - ccomp * Wc * bumpInt((r - Wc) / Wc);
    }

private:
    /// IB(clamp(u,-1,1)) + 16/35, i.e. int_{-1}^{u} B: runs from 0 to 32/35
    static real_t bumpInt(real_t u)
    {
        const real_t t  = (u < -1.0) ? -1.0 : ((u > 1.0) ? 1.0 : u);
        const real_t t2 = t * t;
        return t * (1.0 - t2 * (1.0 - t2 * (0.6 - t2 / 7.0))) + 16.0 / 35.0;
    }

    real_t m_R0, m_taper, m_M, m_w;
};

} // anonymous namespace

int main(int argc, char *argv[])
{
    // Data & fitting options (mirrors fitting_example.cpp where sensible)
    std::string fn = "fitting/deepdrawingC.xml";
    index_t deg_x = 2, deg_y = 2;
    index_t numKnots = 2;      // interior knots of the initial S basis
    index_t numURef = 0;       // initial uniform refinements of S
    real_t lambda = 1e-6;      // fitting smoothing weight
    // Adaptivity schedule
    std::string schedule = "FRH";
    index_t maxIter = 3;       // number of schedule cycles
    real_t tolerance = 1e-2;   // pointwise error tolerance (stop criterion)
    real_t threshold = -1;     // marking threshold (-1: use refPercent)
    real_t refPercent = 0.1;   // fraction of points to mark when threshold<0
    index_t extension = 2;     // cell extension of marked boxes
    // Sigma map
    index_t sigmaDeg = 2;
    index_t sigmaRef = 3;      // sigma mesh = level-0 analysis basis refined this often
    bool noSlide = false;      // fix boundary controls of sigma
    // Prescribed initial sigma (see gsRadialSigma above)
    std::string sigmaInit = "identity"; // identity | radial
    real_t sigmaR0 = 0.30;     // radius of the expansion ring
    real_t sigmaTaper = 0.45;  // sigma is the identity for r >= sigmaTaper
    real_t sigmaGrade = 4.0;   // radial stretch at the ring (>1: EXPANSION)
    real_t sigmaWidth = 0.05;  // half-width of the ring bump
    // Monitor (error spline) and relocation
    real_t theta = 1.0;        // monitor smoothing θ in ω=1/sqrt(1+θf)
    real_t penalty = 1e-2;     // fold penalty of the relocation energy
    index_t errKnots = 5;      // interior knots of the error-spline basis
    index_t errDeg = 2;        // degree of the error-spline basis
    real_t errLambda = 1e-4;   // smoothing weight of the error-spline fit
    bool invertMonitor = false; // monitor on 1-f: expand (not contract) sigma at high error
    bool unitDomain = false;    // relocate against a unit square instead of S: the
                                // Winslow term of the fitted S otherwise dominates
                                // the monitor (S is far from isometric)
    real_t analyticW = -1;      // >0: analytic band monitor of this width around
                                // the lshape_step curve instead of the LS-fitted
                                // error spline (H1 shear-stiffness test)
    // Direct sigma optimization (D step). The barrier is averaged over the
    // grid, so a localized near-fold contributes ~ mu * nFolded/M * v^2:
    // mu must be large (>=1e3) to dominate LS gains of order 1e-2 (mu=1
    // provably lets sigma fold and the step gets rejected).
    real_t dirMu = 1000.0;     // weight of the fold/box penalty of the D step
    real_t dirEps = 5e-2;      // det J_sigma floor of the fold barrier
    index_t dirBarrierMode = 0; // 0: Gauss quadrature on sigma's knot mesh
                               // (resolution delegated to the basis, scales
                               // with sigmaRef); 1: Bezier-coefficient
                               // barrier (no quadrature, quB ignored)
    index_t dirQuB = 8;        // extra Gauss nodes per direction of the
                               // sigma-mesh barrier rule (nodes = deg + quB).
                               // ~10 nodes per element per direction are
                               // needed or HLBFGS folds sigma between the
                               // nodes (validated at p_sigma = 2: quB = 4
                               // folds, quB = 8 is clean and reproduces the
                               // hand-tuned uniform-grid results)
    index_t dirSkip = 0;       // skip D steps in the first dirSkip cycles: adapting
                               // sigma against an S that resolves none of the
                               // features distorts the map and poisons later cycles
    bool dirCheckGrad = false; // FD-validate the D-step gradient once
    // Directional tensor refinement (T step)
    index_t refineDir = -1;    // 0 or 1: refine that parametric direction only;
                               // -1: both directions (isotropic)
    index_t refineKnots = 1;   // knots inserted per span by the T step
    index_t optIter = 100;     // HLBFGS iterations per R step
    index_t optVerbose = 1;    // HLBFGS verbosity
    real_t optTol = 1e-6;      // HLBFGS minimal gradient length
    // Output
    std::string output = "rh_fitting_output";
    index_t nSamples = 1000;
    bool plot = false;
    bool save = false;

    gsCmdLine cmd("rh-adaptive fitting of a parametrized point cloud with a "
                  "composed THB geometry G = S(sigma(u)). The --schedule "
                  "string over {F,R,D,H,T} (Fit, r-adapt, direct sigma opt, "
                  "h-refine, directional tensor refine) is repeated up to -i "
                  "times; e.g. FRH, FDH, FRFRH, FT, F. T and H are mutually "
                  "exclusive (T rebuilds S from the plain tensor basis).");
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addInt("x", "deg_x", "Degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "Degree in y direction", deg_y);
    cmd.addInt("n", "iknots", "Interior knots of the initial basis (each direction)", numKnots);
    cmd.addInt("r", "urefine", "Initial uniform refinement steps", numURef);
    cmd.addReal("s", "lambda", "Fitting smoothing coefficient", lambda);
    cmd.addString("", "schedule", "Cycle string over {F,R,H,D,T}", schedule);
    cmd.addInt("i", "iter", "Maximum number of schedule cycles", maxIter);
    cmd.addReal("e", "tolerance", "Pointwise error tolerance (stop criterion)", tolerance);
    cmd.addReal("t", "threshold", "Marking threshold (-1: refPercent percentile)", threshold);
    cmd.addReal("p", "refPercent", "Fraction of points to mark for refinement", refPercent);
    cmd.addInt("q", "extension", "Extension size of marked cells", extension);
    cmd.addInt("E", "sigmaDeg", "Degree of the sigma map", sigmaDeg);
    cmd.addInt("R", "sigmaRef", "Uniform refinements of sigma's mesh w.r.t. the "
               "LEVEL-0 analysis basis (2^sigmaRef elements per direction, so "
               "sigma is always dyadically nested with every analysis level)",
               sigmaRef);
    cmd.addSwitch("noslide", "Fix the boundary control points of sigma (no sliding)", noSlide);
    cmd.addString("", "sigmaInit", "Initial sigma map: 'identity' or 'radial' "
                  "(prescribed C^2 radial grading, expanding at r = sigmaR0)", sigmaInit);
    cmd.addReal("", "sigmaR0", "Radius of the ring on which the radial init expands", sigmaR0);
    cmd.addReal("", "sigmaTaper", "Radius beyond which the radial init is the identity "
                "(must be < 0.5)", sigmaTaper);
    cmd.addReal("", "sigmaGrade", "Radial stretch g'(R0) of the radial init "
                "(>1 expands sigma at the ring)", sigmaGrade);
    cmd.addReal("", "sigmaWidth", "Half-width of the radial-init ring bump", sigmaWidth);
    cmd.addReal("T", "theta", "Monitor smoothing parameter theta", theta);
    cmd.addReal("P", "penalty", "Fold penalty of the relocation energy", penalty);
    cmd.addInt("", "errKnots", "Interior knots of the error-spline basis", errKnots);
    cmd.addInt("", "errDeg", "Degree of the error-spline basis", errDeg);
    cmd.addReal("", "errLambda", "Smoothing weight of the error-spline fit", errLambda);
    cmd.addSwitch("invertMonitor", "Drive the monitor with 1-f so that sigma expands "
                  "(det J large) where the error is large", invertMonitor);
    cmd.addSwitch("unitDomain", "Relocation energy on a unit-square template instead "
                  "of the fitted S (monitor-only relocation)", unitDomain);
    cmd.addReal("", "analyticW", "Width of an analytic Gaussian band monitor around "
                "the lshape_step curve (>0 replaces the LS error spline; "
                "lshape_step.xml only)", analyticW);
    cmd.addReal("", "dirMu", "Fold/box penalty weight of the direct (D) step", dirMu);
    cmd.addReal("", "dirEps", "det J_sigma floor of the D-step fold barrier", dirEps);
    cmd.addInt("", "dirBarrierMode", "D-step fold barrier: 0 = sampled Gauss rule "
               "on sigma's mesh, 1 = Bezier-coefficient barrier", dirBarrierMode);
    cmd.addInt("", "dirQuB", "Extra Gauss nodes per direction of the sigma-mesh "
               "barrier rule", dirQuB);
    cmd.addInt("", "dirSkip", "Skip D steps in the first dirSkip cycles", dirSkip);
    cmd.addSwitch("dirCheckGrad", "Finite-difference check of the D-step gradient "
                  "(first D step only)", dirCheckGrad);
    cmd.addInt("", "refineDir", "T-step refinement direction (0, 1, or -1 for both)", refineDir);
    cmd.addInt("", "refineKnots", "Knots inserted per span by the T step", refineKnots);
    cmd.addInt("", "optIter", "HLBFGS iterations per R step", optIter);
    cmd.addInt("", "optVerbose", "HLBFGS verbosity", optVerbose);
    cmd.addReal("", "optTol", "HLBFGS minimal gradient length", optTol);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addInt("", "nSamples", "Number of samples for ParaView output", nSamples);
    cmd.addSwitch("plot", "Create ParaView visualization files", plot);
    cmd.addSwitch("save", "Save result in XML format", save);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(deg_x > 0 && deg_y > 0, "Degrees must be positive");
    GISMO_ENSURE(extension >= 0, "Extension must be non-negative");
    GISMO_ENSURE(tolerance >= 0, "Tolerance must be non-negative");
    GISMO_ENSURE(dirBarrierMode == 0 || dirBarrierMode == 1,
                 "--dirBarrierMode must be 0 (sampled) or 1 (coefficient), got "
                 << dirBarrierMode);
    std::string schedUp;   // upper-cased copy: the schedule is dispatched
                           // case-insensitively, so validate the same way
    for (char c : schedule)
    {
        const char C = std::toupper(c);
        GISMO_ENSURE(C == 'F' || C == 'R' || C == 'H' || C == 'D' || C == 'T',
                     "Invalid schedule character '" << c << "' (use F, R, H, D, T)");
        schedUp.push_back(C);
    }
    GISMO_ENSURE(!(schedUp.find('T') != std::string::npos &&
                   schedUp.find('H') != std::string::npos),
        "T rebuilds S from the tensor basis and would discard THB levels; "
        "T and H are mutually exclusive in a schedule.");
    // gsTensorBasis::uniformRefine indexes m_bases[dir] unchecked: an
    // out-of-range direction is silent UB in a release build.
    GISMO_ENSURE(refineDir >= -1 && refineDir <= 1,
                 "--refineDir must be 0, 1 or -1 (both directions)");
    GISMO_ENSURE(sigmaInit == "identity" || sigmaInit == "radial",
                 "--sigmaInit must be 'identity' or 'radial' (got '" << sigmaInit << "')");
    if (sigmaInit == "radial")
    {
        // A pure radial map would push the corners (r up to 0.707) out of
        // [0,1]^2; keeping the taper inside the inscribed circle also leaves
        // every boundary control (r >= 0.5) at the identity.
        GISMO_ENSURE(sigmaTaper > 0 && sigmaTaper < 0.5,
                     "--sigmaTaper must lie in (0, 0.5): sigma is the identity for "
                     "r >= taper and r reaches 0.5 on the boundary of [0,1]^2");
        // Support of the ring bump must stay inside [0, taper]: otherwise
        // g(0) != 0 (blow-up of g(r)/r at the centre) or g(taper) != taper
        // (jump at the taper ring).
        GISMO_ENSURE(sigmaWidth > 0 && sigmaWidth <= sigmaR0 &&
                     sigmaR0 + sigmaWidth <= sigmaTaper,
                     "the ring bump support [sigmaR0-sigmaWidth, sigmaR0+sigmaWidth] "
                     "must lie inside [0, sigmaTaper]");
        GISMO_ENSURE(sigmaGrade > 0, "--sigmaGrade must be positive");
        // Monotonicity of g: min g' = 1 - (grade-1)*width/(taper/2) > 0
        GISMO_ENSURE((sigmaGrade - 1.0) * sigmaWidth < 0.5 * sigmaTaper,
                     "sigmaGrade too large: need (grade-1)*width < taper/2");
    }
    if (threshold > 0 && threshold > tolerance)
    {
        gsInfo << "Refinement threshold is over tolerance, setting it to the tolerance.\n";
        threshold = tolerance;
    }

    gsFileManager::mkdir(output);
    output += gsFileManager::getNativePathSeparator();

    //! [Read data]
    // id 0: 2 x N parameter values; id 1: d x N point coordinates
    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv);
    fd_in.getId<gsMatrix<> >(1, xyz);
    GISMO_ENSURE(uv.cols() == xyz.cols() && uv.rows() == 2, "Wrong input");
    const index_t N = uv.cols();

    // Rescale the parameters to the sigma domain [0,1]^2
    for (index_t r = 0; r != 2; r++)
    {
        const real_t mn = uv.row(r).minCoeff(), mx = uv.row(r).maxCoeff();
        GISMO_ENSURE(mx > mn, "Degenerate parameter range in direction " << r);
        uv.row(r) = (uv.row(r).array() - mn) / (mx - mn);
    }
    //! [Read data]

    //! [Create bases and sigma]
    // THB basis for S on [0,1]^2
    gsKnotVector<> ku(0, 1, numKnots, deg_x + 1);
    gsKnotVector<> kv(0, 1, numKnots, deg_y + 1);
    gsTensorBSplineBasis<2> tbasis(ku, kv);
    tbasis.uniformRefine((1 << numURef) - 1);
    gsTHBSplineBasis<2> thb(tbasis);

    // Sigma map (identity initially) on [0,1]^2
    // Sigma's mesh is a refinement LEVEL of the level-0 analysis basis, not a
    // free knot count.  makeIntegrationBasis (the paper's super mesh) unions the
    // sigma and analysis knot lines; a non-nested pair costs roughly their SUM
    // of elements, a nested pair only the finer of the two -- cf. Sec.
    // "Numerical integration scheme": "if either the space associated with sigma
    // or the analysis space is nested within the other, no additional knot lines
    // are introduced, leading to a more efficient quadrature process".  Every
    // analysis level is dyadic, so tying sigma to the same ladder makes a
    // non-nested sigma inexpressible.
    GISMO_ENSURE(sigmaRef >= 0, "sigmaRef must be non-negative");
    const index_t sigmaKnots = (1 << sigmaRef) - 1;
    gsKnotVector<> ks(0, 1, sigmaKnots, sigmaDeg + 1);
    gsTensorBSplineBasis<2> sbasis(ks, ks);
    gsSquareDomain<real_t> domain(sbasis);
    domain.options().addSwitch("Slide", "Boundary controls slide along the boundary", !noSlide);
    domain.applyOptions();

    // min det(J_sigma): set by the radial init below and by every R/D step
    real_t minDetJ = 1.0;

    // --sigmaInit radial: replace the identity start of sigma by the analytic
    // map gsRadialSigma, interpolated at the anchors of sbasis. Only the FREE
    // controls are overwritten, so the eliminated boundary components keep
    // their identity values and sigma remains a map of [0,1]^2 onto itself.
    if (sigmaInit == "radial")
    {
        // Resolution warning: interpolation at the anchors can only reproduce
        // g' if the ring bump spans several sigma elements. Below that the
        // interpolant overshoots and folds (measured on circle_band with the
        // defaults: min det J = -0.16 at 15 sigma knots (--sigmaRef 4),
        // +0.013 at 25, 0.13 at
        // 39, 0.22 at 79 against the analytic 1-ccomp times g/r = 0.232).
        const real_t hSigma = 1.0 / (sigmaKnots + 1);
        if (2 * sigmaWidth < 3 * hSigma)
            gsWarn << "Radial sigma init: the ring bump (width " << 2 * sigmaWidth
                   << ") spans only " << 2 * sigmaWidth / hSigma << " sigma elements; "
                      "g' will be aliased and sigma may fold. Raise --sigmaRef "
                      "(>= " << math::ceil(math::log(math::ceil(3.0 / (2 * sigmaWidth))) / math::log(2.0))
                   << " here) or --sigmaWidth.\n";

        gsRadialSigma sigmaA(sigmaR0, sigmaTaper, sigmaGrade, sigmaWidth);
        gsMatrix<> anch = sbasis.anchors(), vals;
        sigmaA.eval_into(anch, vals);                     // 2 x sbasis.size()
        gsGeometry<>::uPtr sgeo = sbasis.interpolateAtAnchors(vals);

        const gsDofMapper & mp = domain.mapper();
        gsVector<> ctrl = domain.getControls();
        for (index_t k = 0; k != sgeo->coefs().rows(); k++)
            for (index_t dcomp = 0; dcomp != 2; dcomp++)
                if (mp.is_free(k, 0, dcomp))
                    ctrl[mp.index(k, 0, dcomp)] = sgeo->coefs()(k, dcomp);
        domain.setControls(ctrl);

        // The realized sigma is NOT sgeo (the eliminated boundary components
        // stay at the identity), so the fold check must run on the domain.
        // The interpolant also only resolves the ring bump if sigma's mesh is
        // fine enough: with 2w below the sigma element size, g' is aliased and
        // the realized det J drops well below the analytic min 1 - ccomp.
        minDetJ = minDetJSigma(domain);
        gsInfo << "min det J_sigma (radial init) = " << minDetJ << "\n";
        GISMO_ENSURE(minDetJ > 0,
                     "Radial sigma init is folded (min det J_sigma = " << minDetJ
                     << "): lower --sigmaGrade, widen --sigmaWidth, or refine "
                     "sigma with --sigmaRef");
    }

    // Basis for the error monitor spline on the S-domain
    gsKnotVector<> ke(0, 1, errKnots, errDeg + 1);
    gsTensorBSplineBasis<2> ebasis(ke, ke);
    //! [Create bases and sigma]

    std::vector<unsigned> ext(2, static_cast<unsigned>(extension));

    gsInfo << "Fitting " << N << " samples.\n";
    gsInfo << "----------------\n";
    gsInfo << "Schedule           : " << schedule << " x " << maxIter << "\n";
    gsInfo << "S basis            : deg (" << deg_x << "," << deg_y << "), "
           << thb.size() << " DoFs\n";
    gsInfo << "Sigma basis        : deg " << sigmaDeg << ", "
           << domain.nControls() << " free controls"
           << (noSlide ? " (fixed boundary)" : " (sliding boundary)") << "\n";
    gsInfo << "Error tolerance    : " << tolerance << "\n";
    if (threshold >= 0)
        gsInfo << "Ref. threshold     : " << threshold << "\n";
    else
        gsInfo << "Cell refinement    : " << 100 * refPercent << "%\n";
    gsInfo << "Smoothing lambda   : " << lambda << "\n";
    gsInfo << "Monitor theta      : " << theta << "\n";
    gsInfo << "----------------\n";

    // ------------------------------------------------------------------
    // State shared by the F/R/H steps
    // ------------------------------------------------------------------
    gsGeometry<>::uPtr geom;             // current S
    gsGeometry<>::uPtr errFun;           // last monitor spline (for plotting)
    gsMatrix<> xi;                       // ξ_i = σ(uv_i), 2 x N
    std::vector<real_t> errors(N, 0.0);  // pointwise errors ||G(uv_i)-xyz_i||
    real_t minErr = 0, maxErr = 0, rmse = 0;
    // (minDetJ is declared with the sigma map above: the radial init sets it
    //  before any step runs)

    auto updateParams = [&]()
    {
        domain.eval_into(uv, xi);
        // Guard against numerical over/undershoot of σ
        xi = xi.cwiseMax(0.0).cwiseMin(1.0);
    };

    auto computeErrors = [&]()
    {
        gsMatrix<> vals;
        geom->eval_into(xi, vals);
        minErr = std::numeric_limits<real_t>::max();
        maxErr = rmse = 0;
        for (index_t i = 0; i != N; i++)
        {
            errors[i] = (vals.col(i) - xyz.col(i)).norm();
            minErr = std::min(minErr, errors[i]);
            maxErr = std::max(maxErr, errors[i]);
            rmse += errors[i] * errors[i];
        }
        rmse = math::sqrt(rmse / N);
    };

    // F: penalized least-squares fit of S at ξ = σ(uv)
    auto fitStep = [&]()
    {
        updateParams();
        gsMatrix<> coefs = leastSquaresFit(thb, xi, xyz, lambda);
        geom = thb.makeGeometry(give(coefs));
        computeErrors();
    };

    // R: relocate σ, driven by an error-monitor spline over the S-domain
    // (or by the analytic band monitor when --analyticW is set)
    auto rStep = [&]()
    {
        gsAnalyticBandMonitor bandMon(analyticW, invertMonitor);
        const gsFunction<real_t> * monitor = nullptr;
        if (analyticW > 0)
            monitor = &bandMon;
        else
        {
            // Monitor: normalized pointwise errors fitted on the coarse basis
            gsMatrix<> e(1, N);
            const real_t scale = (maxErr > 0) ? maxErr : 1.0;
            // Standard orientation contracts sigma (det J small) where f is large;
            // for fitting that spreads the knot preimages apart at the feature.
            // Inverted (f~ = 1-f) sigma expands there: knot preimages cluster at
            // the high-error region, raising the effective resolution of S(sigma).
            for (index_t i = 0; i != N; i++)
                e(0, i) = invertMonitor ? 1.0 - errors[i] / scale
                                        : errors[i] / scale;
            gsMatrix<> ecoefs = leastSquaresFit(ebasis, xi, e, errLambda);
            // Clamp the coefficients: the B-spline convex hull then guarantees
            // f in [0,1], keeping the monitor 1/sqrt(1+theta*f) well-defined
            // (least-squares overshoot in data-poor cells can otherwise make
            // 1+theta*f negative and NaN the relocation energy).
            ecoefs = ecoefs.cwiseMax(0.0).cwiseMin(1.0);
            errFun = ebasis.makeGeometry(give(ecoefs));
            monitor = errFun.get();
        }

        gsHLBFGS<real_t> optimizer;
        optimizer.options().setInt("MaxIterations", optIter);
        optimizer.options().setInt("Verbose", optVerbose);
        optimizer.options().setReal("MinGradLen", optTol);

        // gsAdaptiveParametrization needs a tensor integration basis: use the
        // finest tensor level of the THB basis (never coarser than S).
        const gsTensorBSplineBasis<2> & itensor = thb.tensorLevel(thb.maxLevel());
        // The relocation energy is a monitor-weighted Winslow energy of the
        // composed map: with the fitted S its (unweighted) Winslow term measures
        // S's own parametric distortion and can dominate the monitor entirely.
        // --unitDomain substitutes an isometric template so that only the
        // monitor drives sigma.
        gsGeometry<>::uPtr unitSq;
        if (unitDomain)
            unitSq = gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0);
        gsAdaptiveParametrization<real_t, MonitorMode::ValueBased>
            reloc(domain, unitDomain ? *unitSq : *geom, *monitor, itensor, optimizer, true);
        reloc.options().setReal("Penalty", penalty);
        reloc.options().setReal("Smoothing", theta);

        const gsVector<real_t> uPrev = domain.getControls();
        reloc.solve();
        minDetJ = reloc.computeMinJacobian();

        // Reject a diverged or folded relocation step
        const gsVector<real_t> uNew = domain.getControls();
        if (!uNew.allFinite() || minDetJ <= 0)
        {
            gsWarn << "Relocation step rejected (min det Js = " << minDetJ
                   << "); sigma restored.\n";
            domain.setControls(uPrev);
            minDetJ = reloc.computeMinJacobian();
        }

        // σ moved: the data points now sample S elsewhere
        updateParams();
        computeErrors();

        // Diagnostic: does det J_sigma follow the monitor? corr(detJ, f) at
        // the data points separates "relocation ignores the monitor" from
        // "monitor spline misrepresents the error field" (corr(detJ, err)).
        {
            gsMatrix<> dsigma, fvals;
            domain.deriv_into(uv, dsigma);
            monitor->eval_into(xi, fvals);
            auto corr = [&](const std::vector<real_t> & a, const std::vector<real_t> & b)
            {
                real_t ma = 0, mb = 0;
                for (index_t i = 0; i != N; i++) { ma += a[i]; mb += b[i]; }
                ma /= N; mb /= N;
                real_t sab = 0, saa = 0, sbb = 0;
                for (index_t i = 0; i != N; i++)
                {
                    sab += (a[i] - ma) * (b[i] - mb);
                    saa += (a[i] - ma) * (a[i] - ma);
                    sbb += (b[i] - mb) * (b[i] - mb);
                }
                return (saa > 0 && sbb > 0) ? sab / math::sqrt(saa * sbb) : 0.0;
            };
            std::vector<real_t> detJ(N), f(N);
            for (index_t i = 0; i != N; i++)
            {
                detJ[i] = dsigma(0, i) * dsigma(3, i) - dsigma(1, i) * dsigma(2, i);
                f[i] = fvals(0, i);
            }
            gsInfo << "    [R diag] corr(detJ, monitor f) = " << corr(detJ, f)
                   << ", corr(detJ, error) = " << corr(detJ, errors) << "\n";

            // H1 contrast measure: mean det J of the points mapped into the
            // band (dist <= w) vs the far quartile (largest dist). The
            // shear-stiffness hypothesis predicts the ratio stays ~1.2 for
            // thin bands and grows once w spans several sigma elements.
            if (analyticW > 0)
            {
                std::vector<real_t> dist(N);
                for (index_t i = 0; i != N; i++)
                    dist[i] = distToStep(xi(0, i), xi(1, i));
                std::vector<real_t> sorted = dist;
                const size_t q3 = (3 * static_cast<size_t>(N)) / 4;
                std::nth_element(sorted.begin(), sorted.begin() + q3, sorted.end());
                const real_t farThr = sorted[q3];
                real_t mBand = 0, mFar = 0, mn = detJ[0], mx = detJ[0];
                index_t nBand = 0, nFar = 0;
                for (index_t i = 0; i != N; i++)
                {
                    mn = std::min(mn, detJ[i]);
                    mx = std::max(mx, detJ[i]);
                    if (dist[i] <= analyticW) { mBand += detJ[i]; nBand++; }
                    if (dist[i] >= farThr)    { mFar  += detJ[i]; nFar++;  }
                }
                if (nBand > 0) mBand /= nBand;
                if (nFar > 0)  mFar  /= nFar;
                gsInfo << "    [H1] w = " << analyticW
                       << " | band mean detJ = " << mBand << " (" << nBand << " pts)"
                       << " | far mean detJ = " << mFar << " (" << nFar << " pts)"
                       << " | contrast = " << ((mFar > 0) ? mBand / mFar : 0.0)
                       << " | detJ range = [" << mn << ", " << mx << "]\n";
            }
        }
    };

    // D: direct minimization of the true LS fitting error over the sigma
    // controls (variable projection with F). Unlike R, this cannot increase
    // the LS error beyond line-search tolerance: the optimizer descends on
    // the reported objective itself.
    bool dBarrierPrinted = false;
    auto dStep = [&]()
    {
        const gsFoldBarrierMode dirMode = dirBarrierMode == 0
            ? gsFoldBarrierMode::Sampled : gsFoldBarrierMode::Coefficient;
        gsOptFit<real_t> obj(domain, *geom, uv, xyz, dirMu, dirEps,
                             dirMode, dirQuB);
        if (!dBarrierPrinted)
        {
            gsInfo << "    [D] barrier quadrature: " << obj.nBarrierPts()
                   << " points ("
                   << (dirBarrierMode == 1 ? "Bezier-coefficient barrier"
                                           : "Gauss on the sigma knot mesh")
                   << ")\n";
            dBarrierPrinted = true;
        }
        if (dirCheckGrad)
        {
            gsInfo << "    [D diag] max rel FD gradient error = "
                   << gsCheckSigmaGradient(obj, domain) << "\n";
            dirCheckGrad = false;
        }
        gsHLBFGS<real_t> optimizer;
        optimizer.options().setInt("MaxIterations", optIter);
        optimizer.options().setInt("Verbose", optVerbose);
        optimizer.options().setReal("MinGradLen", optTol);
        optimizer.setProblem(&obj);

        const gsVector<real_t> uPrev = domain.getControls();
        const real_t rmsePrev = rmse;
        optimizer.solve(uPrev);
        gsVector<real_t> uNew = optimizer.currentDesign();
        domain.setControls(uNew);

        // Fold check on a per-element grid of sigma's knot mesh (a fixed
        // global grid undersamples fine sigma meshes, like the barrier did)
        minDetJ = minDetJSigma(domain);

        // Reject a diverged or folded step (same policy as R)
        if (!uNew.allFinite() || minDetJ <= 0)
        {
            gsWarn << "Direct step rejected (min det Js = " << minDetJ
                   << "); sigma restored.\n";
            domain.setControls(uPrev);
            minDetJ = minDetJSigma(domain);
        }
        updateParams();
        computeErrors();
        if (rmse > rmsePrev && rmsePrev > 0)
            gsInfo << "    [D diag] rmse increased (" << rmsePrev << " -> " << rmse
                   << "): fold/box penalty active or line-search tolerance.\n";
    };

    // H: THB refinement of over-threshold cells (marking in the S-domain).
    // Refining basis AND geometry keeps G unchanged (and the errors valid)
    // until the next fit.
    auto hStep = [&]() -> bool
    {
        const real_t thr = (threshold >= 0) ? threshold
                                            : refineThreshold(errors, refPercent);
        std::vector<index_t> boxes = getBoxes(thb, errors, xi, thr, ext);
        if (boxes.empty())
            return false;
        thb.refineElements(boxes);
        geom->refineElements(boxes);
        return true;
    };

    // T: directional refinement of the plain tensor basis. NEVER route this
    // through the THB basis: gsHTensorBasis::uniformRefine(dir) calls a
    // direction-agnostic m_tree.multiplyByTwo() and its dir==-1 assert is
    // compiled out under NDEBUG (silent wrong answer). tbasis is a stale
    // copy once thb.refineElements ran, which is why T excludes H.
    // Unlike H, T does NOT carry the geometry along (the basis is rebuilt,
    // not refined inside the hierarchy): geom and the errors stay those of
    // the pre-T fit until the next F, so on a T row of the CSV dofs_S is
    // post-refinement while the error columns are the previous fit's.
    // Consequence for schedule design: put an F between T and R (e.g. FTFR,
    // not FTR), otherwise the R step relocates sigma against the pre-T
    // geometry and error monitor and gains nothing from T's new DoFs.
    // (the `dirty` flag is set by the caller in the dispatch switch below,
    // where it is in scope, exactly as for the F/R/D/H steps)
    auto tStep = [&]()
    {
        tbasis.uniformRefine(refineKnots, 1, static_cast<short_t>(refineDir));
        thb = gsTHBSplineBasis<2>(tbasis);
    };

    // ------------------------------------------------------------------
    // Adaptive loop
    // ------------------------------------------------------------------
    std::ofstream csv(output + "convergence.csv");
    csv << "cycle,step,dofs_S,dofs_sigma,minErr,maxErr,rmse,pctBelowTol,minDetJsigma,time\n";

    gsStopwatch timer;
    bool converged = false;
    bool dirty = false;   // true when sigma or the basis changed after the last fit
    for (index_t it = 0; it != maxIter && !converged; ++it)
    {
        gsInfo << "----------------\n";
        gsInfo << "Cycle " << it << " [" << schedule << "]\n";
        for (char c : schedule)
        {
            const char C = std::toupper(c);
            timer.restart();
            bool refined = true;
            switch (C)
            {
            case 'F':
                fitStep();
                dirty = false;
                break;
            case 'R':
                if (!geom) fitStep(); // R needs errors: bootstrap fit
                rStep();
                dirty = true;
                break;
            case 'D':
                if (it < dirSkip) continue; // sigma adaptation deferred
                if (!geom) fitStep(); // D needs S: bootstrap fit
                dStep();
                dirty = true;
                break;
            case 'H':
                if (!geom) fitStep(); // H needs errors: bootstrap fit
                refined = hStep();
                dirty = dirty || refined;
                break;
            case 'T':
                tStep();              // no bootstrap fit: T needs no errors
                dirty = true;
                break;
            }
            const real_t time = timer.stop();

            const real_t pctBelow =
                100.0 * std::count_if(errors.begin(), errors.end(),
                                      [&](real_t e) { return e < tolerance; }) / N;

            gsInfo << "  " << C
                   << " | S DoFs: " << std::setw(6) << thb.size()
                   << " | sigma DoFs: " << std::setw(4) << domain.nControls()
                   << " | err (min/max): " << std::scientific << std::setprecision(3)
                   << minErr << " / " << maxErr
                   << " | rmse: " << rmse
                   << " | < tol: " << std::fixed << std::setprecision(1) << pctBelow << "%"
                   << " | min det Js: " << std::scientific << std::setprecision(2) << minDetJ
                   << " | " << std::fixed << std::setprecision(2) << time << "s\n"
                   << std::defaultfloat;

            csv << it << "," << C << "," << thb.size() << "," << domain.nControls()
                << "," << minErr << "," << maxErr << "," << rmse << "," << pctBelow
                << "," << minDetJ << "," << time << "\n";

            if (C == 'H' && !refined)
                gsInfo << "    (no cells marked for refinement)\n";

            if (geom && maxErr < tolerance)
            {
                gsInfo << "Error tolerance achieved in cycle " << it << ".\n";
                converged = true;
                break;
            }
        }
    }
    // The schedule ended on an R/H step: conclude with a fit so that the
    // reported result reflects the final basis and sigma.
    if (dirty && !converged)
    {
        timer.restart();
        fitStep();
        const real_t time = timer.stop();
        const real_t pctBelow =
            100.0 * std::count_if(errors.begin(), errors.end(),
                                  [&](real_t e) { return e < tolerance; }) / N;
        gsInfo << "  F (final)"
               << " | S DoFs: " << std::setw(6) << thb.size()
               << " | err (min/max): " << std::scientific << std::setprecision(3)
               << minErr << " / " << maxErr
               << " | rmse: " << rmse
               << " | < tol: " << std::fixed << std::setprecision(1) << pctBelow << "%\n"
               << std::defaultfloat;
        csv << maxIter << ",F," << thb.size() << "," << domain.nControls()
            << "," << minErr << "," << maxErr << "," << rmse << "," << pctBelow
            << "," << minDetJ << "," << time << "\n";
    }
    csv.close();
    gsInfo << "----------------\n";
    gsInfo << "Finished: " << (converged ? "tolerance reached" : "iteration budget spent")
           << ". Final S DoFs: " << thb.size()
           << ", max error: " << maxErr << ", rmse: " << rmse << "\n";

    // ------------------------------------------------------------------
    // Monitor-orientation diagnostic: correlation of det J_sigma with the
    // pointwise error. Inverted monitor should give corr > 0 (sigma expands
    // at high error -> knot preimages cluster there); standard gives corr < 0.
    // ------------------------------------------------------------------
    if (geom)
    {
        gsMatrix<> dsigma;
        domain.deriv_into(uv, dsigma); // 4 x N, column: [dx/du dx/dv dy/du dy/dv]
        std::ofstream djcsv(output + "detj_vs_error.csv");
        djcsv << "u,v,detJsigma,error\n";
        real_t mJ = 0, mE = 0;
        std::vector<real_t> detJ(N);
        for (index_t i = 0; i != N; i++)
        {
            detJ[i] = dsigma(0, i) * dsigma(3, i) - dsigma(1, i) * dsigma(2, i);
            djcsv << uv(0, i) << "," << uv(1, i) << ","
                  << detJ[i] << "," << errors[i] << "\n";
            mJ += detJ[i]; mE += errors[i];
        }
        djcsv.close();
        mJ /= N; mE /= N;
        real_t sJE = 0, sJJ = 0, sEE = 0;
        for (index_t i = 0; i != N; i++)
        {
            sJE += (detJ[i] - mJ) * (errors[i] - mE);
            sJJ += (detJ[i] - mJ) * (detJ[i] - mJ);
            sEE += (errors[i] - mE) * (errors[i] - mE);
        }
        // The summary blocks above leave the stream in fixed/setprecision(1),
        // which would round the correlation to a single digit ("0.1"): reset
        // to default float format with 4 significant digits.
        gsInfo.unsetf(std::ios_base::floatfield);
        gsInfo << std::setprecision(4);
        if (sJJ > 0 && sEE > 0)
            gsInfo << "corr(det J_sigma, error) = "
                   << sJE / math::sqrt(sJJ * sEE)
                   << "  (monitor " << (invertMonitor ? "inverted" : "standard")
                   << ")\n";
    }

    // ------------------------------------------------------------------
    // Output
    // ------------------------------------------------------------------
    if (plot && geom)
    {
        gsComposedGeometry<real_t> cgeom(domain.domain(), *geom);
        gsMultiPatch<> mp;
        mp.addPatch(cgeom);
        gsWriteParaview(mp, output + "cgeom", nSamples, true);
        gsWriteParaview(*geom, output + "geom", nSamples, true);
        gsWriteParaview(domain.domain(), output + "domain", nSamples, true, true);
        if (errFun)
            gsWriteParaview(*errFun, output + "errorSpline", nSamples);
        gsWriteParaviewPoints(xyz, output + "points");
        gsInfo << "ParaView files written to " << output << "\n";
    }

    if (save && geom)
    {
        gsFileData<> fd;
        fd << *geom;             // fitted S
        fd << domain.domain();   // sigma map
        fd.dump(output + "rh_fitting_out");
        gsInfo << "Result written to " << output << "rh_fitting_out.xml\n";
    }
    else if (!save)
        gsInfo << "Done. Re-run with --save to export the result as XML.\n";

    return 0;
}
