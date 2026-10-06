/** @file momfit_operator_probe_example.cpp

    @brief Do the negative weights of solve-free moment fitting damage the
           assembled operator?  A measurement, not a claim.

    Every verification of the moment-fitting rule so far integrates SCALARS
    (area, volume, polynomials).  Mass is preserved identically by the rule
    (partition of unity of the nodal Lagrange basis), so those checks say
    nothing about what Stats::nNegativeWeights does to a MATRIX.  This driver
    assembles the immersed mass matrix

        M_ij = \int_{\Omega} phi_i phi_j ,   \Omega = { phi_LS < 0 },

    twice on the same background mesh and the same elements -- once with
    gsAlgoimAdaptiveRule (the accurate reference) and once with the core
    gsMomentRule wrapped around an identically configured gsAlgoimAdaptiveRule
    -- and reports the spectrum of both.  Compression fires whenever the
    underlying point cloud exceeds the n^d output grid (the pass-through
    guard), which at every configuration below it does on every cut element.

    Geometry: a disk immersed in [0,1]^2, level set

        phi_LS(x,y) = sqrt((x-cx)^2 + (y-cy)^2) - R          ( < 0 inside )

    This is the EXACT signed distance function of the disk, hence 1-Lipschitz,
    which is what makes the midpoint box classification of both rules
    (|phi(mid)| > L*halfdiag  =>  uncut, with the default LipschitzConstant=1)
    sound.  Note that the often-used quadratic form (x-cx)^2+(y-cy)^2-R^2 is
    NOT a distance function (|grad| = 2r), so it would need L >= 2*(R+diag).

    QUADRATURE ORDER (important).  Moment fitting is exact for integrands that
    its nodal grid interpolates, i.e. tensor degree <= n-1 per variable, where
    n = quA*p + quB points per direction.  The mass integrand phi_i*phi_j has
    tensor degree 2p, so a meaningful probe needs

        n >= 2p + 1        (default here: quA = 2, quB = 1  =>  n = 2p+1)

    At quA = 1 (n = p+1, the setting of quadrature_circle_area_example) one
    measures UNDER-INTEGRATION, not compression.  The driver enforces the
    inequality and prints n.

    DOF RETENTION.  With alpha = 0 the outside material is dropped, so basis
    functions supported entirely outside \Omega receive no quadrature point at
    all and their row/column of M is EXACTLY zero.  Their eigenvalue is a
    structural zero whose sign is pure roundoff, so they are excluded: DoF i is
    retained iff M_adaptive(i,i) > 0.  The criterion is taken from the
    reference rule, so both matrices are compared on the identical index set;
    the retained count is printed.  Small-but-positive diagonals are KEPT --
    they are the badly-cut DoFs, i.e. exactly the ones under investigation.

    Reported per refinement and per rule:
      sum(w)      total quadrature mass (must equal the disk area; the two
                  rules must agree to ~1e-12 -- a partition-of-unity identity,
                  NOT evidence of accuracy)
      lambda_min  smallest eigenvalue of the retained block (dense solver)
      lambda_max  largest eigenvalue
      cond        lambda_max/lambda_min (meaningless, flagged, if lambda_min<=0)
      lam_minJ    smallest eigenvalue of the Jacobi-scaled block (see below);
      condJ       and its condition number.  The SPD decision is taken here,
                  because the raw lambda_min of an immersed mass matrix sits at
                  the roundoff level of lambda_max and its sign is noise.
      SPD         lam_minJ > 0 ?  (inertia of the scaled block = inertia of M)
      min diag    min_i M(i,i) over the retained set.  Independent of the
                  eigensolver: min diag < 0 alone proves indefiniteness, since
                  e_i^T M e_i = M_ii.
      nNeg/minW   Stats::nNegativeWeights and Stats::minWeight of the
                  moment-fitting rule (reset per refinement).

    Examples:
      ./bin/momfit_operator_probe_example
      ./bin/momfit_operator_probe_example -r 4 -e 3 --quadDepth 3
      ./bin/momfit_operator_probe_example --cx 0.4871 --cy 0.5133 --radius 0.3771

    MEASURED (2026-08-03, disk in [0,1]^2, r = 0..3, indicator = uniform).
    Nine configurations; every one produced negative weights except the last:

      A  defaults (p=2, n=5=2p+1, depth 2, centred)      SPD at all r
      B  off-centre/non-dyadic (cx .4871, cy .5133,
         R .3771 -> generic sliver cuts)                 SPD at all r
      C  quA=1 (n=3 < 2p+1, UNDER-integrated)            INDEFINITE at all r
      D  p=3 (n=7=2p+1)                                  marginal flag at r=2
      D1 p=3, quadDepth=1                                marginal flag at r=2
      D2 p=3, quadDepth=3                                no flag at all
      D3 p=3, off-centre as B                            no flag at all
      E  quadDepth=0 (crudest underlying rule)           SPD at all r
      F  alpha=0.1                                       SPD, and no negative
                                                         weight at all
                                                         (NOTE: --alpha is no
                                                         longer available, see
                                                         the guard in main())

    Reading of it: negative WEIGHTS alone do not make the operator indefinite.
    With 312 negative weights at r=3, configuration A reproduces the reference
    spectrum to every printed digit.  What does break positivity is dropping
    below the exactness order (C, the quA=1 default of the other quadrature
    examples): there a diagonal entry itself goes negative (-1.3e-06 at r=3
    against a reference +9.7e-17), at every refinement, resolved by both tests.

    D is NOT a counter-example at exactness order.  The flagged DoF has
    reference diagonal 1.195e-25 against lam_max 1.35e-02 -- diagRange 3.0e+22,
    i.e. eps*diagRange = 6.7e+06 -- so it sits far below the numerical
    resolution of M, and its sign is not a property of the rule: the same DoF
    gives minDiag = -9.0e-26 at depth 2, -6.8e-26 at depth 1 and +8.4e-26 at
    depth 3, and does not reappear when the disk is moved off the mesh
    symmetry.  The driver therefore prints it as "no*"/MARGINAL and keeps it
    out of the verdict.  Such DoFs are the ones an immersed method must remove
    or stabilise anyway (the REFERENCE rule's own cond is 1e+13..1e+16 there),
    and they are the case non-negative moment fitting (arXiv:2604.15921) would
    exclude by construction.  Blending (F) removes negative weights altogether
    at the price of integrating a different measure.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#include <gsAlgoim/gsAlgoimAdaptiveRule.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

// ---------------------------------------------------------------------------
// Assembly of the immersed mass matrix by a direct element loop.
//
// The embedding map is the identity unit square, so the parametric and the
// physical frame coincide (|det J| = 1) and the element quadrature weights are
// physical weights.  Complexity: O(nElements * nQP * (p+1)^{2d}).
//
// \a elemMass receives, per visited element (in iteration order, identical for
// every rule since the element loop is the rule-independent domain loop), the
// sum of the emitted weights.  That per-element trace is what localises a
// mass mismatch between two rules to the element that classified differently.
// ---------------------------------------------------------------------------
template <class Rule>
void assembleMass(const gsBasis<real_t>              & basis,
                  const gsImplicitTrimmedDomain<2,real_t> & domain,
                  const Rule                        & rule,
                  gsSparseMatrix<real_t>            & M,
                  std::vector<real_t>               & elemMass,
                  real_t                            & minWeight)
{
    const index_t nDofs = basis.size();

    gsSparseEntries<real_t> entries;
    gsMatrix<real_t>  pts, vals;
    gsVector<real_t>  wts;
    gsMatrix<index_t> act;

    elemMass.clear();
    minWeight = std::numeric_limits<real_t>::max();

    for (auto it = domain.beginAll(); it != domain.endAll(); ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);

        real_t mass = 0;
        for (index_t k = 0; k != wts.size(); ++k)
        {
            mass += wts[k];
            minWeight = math::min(minWeight, wts[k]);
        }
        elemMass.push_back(mass);

        if (0 == pts.cols())
            continue;                       // outside element, alpha = 0

        basis.active_into(pts, act);
        basis.eval_into  (pts, vals);

        for (index_t k = 0; k != pts.cols(); ++k)
            for (index_t i = 0; i != act.rows(); ++i)
                for (index_t j = 0; j != act.rows(); ++j)
                    entries.add(act(i,k), act(j,k),
                                wts[k] * vals(i,k) * vals(j,k));
    }

    M.resize(nDofs, nDofs);
    M.setFrom(entries);          // duplicate triplets are summed
}

// ---------------------------------------------------------------------------
// Spectrum of the retained block of \a M (rows/cols listed in \a keep), raw
// and Jacobi-scaled.
//
// Why the scaled numbers are the ones the verdict rests on: on a badly cut
// mesh the raw diagonal of an immersed mass matrix spans many orders of
// magnitude (a DoF supported on a sliver has M_ii ~ 1e-18 here), so
// lambda_min lands at the roundoff level of lambda_max and its SIGN is noise
// -- for either rule.  The congruence transform
//
//     \hat M = D^{-1/2} M D^{-1/2},   D = diag(M_adaptive) > 0
//
// preserves inertia exactly (Sylvester's law), so sign(lambda_min(\hat M))
// answers the definiteness question, but is computed on a matrix with unit
// diagonal where it is numerically meaningful.  The SAME D (from the
// reference rule) scales both matrices, so the two columns stay comparable.
// ---------------------------------------------------------------------------
struct Spectrum
{
    real_t lambdaMin  = 0;   // raw
    real_t lambdaMax  = 0;
    real_t minDiag    = 0;
    real_t lambdaMinS = 0;   // Jacobi-scaled (inertia-preserving)
    real_t lambdaMaxS = 0;
};

Spectrum spectrumOf(const gsSparseMatrix<real_t> & M,
                    const std::vector<index_t>  & keep,
                    const gsVector<real_t>      & scale)  // sqrt(diag(M_adaptive))
{
    const index_t m = static_cast<index_t>(keep.size());
    const gsMatrix<real_t> dense = M.toDense();

    gsMatrix<real_t> S(m, m), Sh(m, m);
    for (index_t a = 0; a != m; ++a)
        for (index_t b = 0; b != m; ++b)
        {
            S (a,b) = dense(keep[a], keep[b]);
            Sh(a,b) = S(a,b) / (scale[a] * scale[b]);
        }

    Spectrum s;
    s.minDiag = std::numeric_limits<real_t>::max();
    for (index_t a = 0; a != m; ++a)
        s.minDiag = math::min(s.minDiag, S(a,a));

    gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> es(S);
    GISMO_ENSURE(gsEigen::Success == es.info(), "Eigen solver failed.");
    s.lambdaMin = es.eigenvalues().minCoeff();
    s.lambdaMax = es.eigenvalues().maxCoeff();

    gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> esh(Sh);
    GISMO_ENSURE(gsEigen::Success == esh.info(), "Eigen solver failed (scaled).");
    s.lambdaMinS = esh.eigenvalues().minCoeff();
    s.lambdaMaxS = esh.eigenvalues().maxCoeff();
    return s;
}

// Pretty-print helpers -------------------------------------------------------
static std::string fmtSci(real_t v, int prec = 6)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

static std::string fmtCond(real_t lmin, real_t lmax)
{
    if (lmin <= 0)
        return "(indef)";
    return fmtSci(lmax / lmin, 2);
}

int main(int argc, char *argv[])
{
    index_t numRefine = 3;      // levels r = 0..numRefine, 2^(r+1) elements/dir
    index_t degree    = 2;
    index_t quadDepth = 2;      // maxDepth of the underlying adaptive rule
    real_t  radius    = 0.4;
    real_t  alpha     = 0.0;
    real_t  cx        = 0.5;
    real_t  cy        = 0.5;
    real_t  quA       = 2.0;    // n = quA*p + quB points/direction, need n >= 2p+1
    index_t quB       = 1;
    std::string indicator = "uniform";

    gsCmdLine cmd("Probe: do moment-fitted negative weights make the immersed "
                  "mass matrix indefinite?");
    cmd.addInt ("r", "refine",    "Number of refinement levels (mesh 2^(r+1) per dir)", numRefine);
    cmd.addInt ("e", "degree",    "Degree of the background tensor B-spline basis", degree);
    cmd.addInt ("d", "quadDepth", "Adaptive subdivision depth (maxDepth) of the "
                                  "underlying cut-cell rule", quadDepth);
    cmd.addReal("R", "radius",    "Disk radius", radius);
    cmd.addReal("a", "alpha",     "Fictitious-domain weight of the outside material "
                                  "(only 0 is supported)", alpha);
    cmd.addReal("x", "cx",        "Disk centre, x", cx);
    cmd.addReal("y", "cy",        "Disk centre, y", cy);
    cmd.addReal("A", "quA",       "Output nodes per direction: n = quA*degree + quB", quA);
    cmd.addInt ("B", "quB",       "Output nodes per direction: n = quA*degree + quB", quB);
    cmd.addString("i", "indicator", "Adaptive indicator: uniform | fallback | integralChange",
                  indicator);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Capability regression of the move to the core gsMomentRule: the FCM
    // blend w <- (1-alpha)w + alpha*w^G needs the analytic full-cell Gauss
    // weight w^G, which the geometry-free compressor does not have.  (Its own
    // optional alpha is the integrand's material field, a different object.)
    GISMO_ENSURE(alpha == (real_t)0,
                 "--alpha (fictitious-domain blend) is not available with the "
                 "core gsMomentRule; use --alpha 0.");

    const index_t n1d = math::max((index_t)1,
        (index_t)std::lround(quA * (double)degree) + quB);

    gsInfo << "=== Moment-fitting operator probe: mass-matrix positivity ===\n\n";
    gsInfo << "Level set     : sqrt((x-" << cx << ")^2+(y-" << cy << ")^2) - " << radius
           << "   (exact SDF, 1-Lipschitz)\n";
    gsInfo << "Background    : [0,1]^2, tensor B-spline, degree " << degree
           << ", meshes 2^(r+1) elements/dir, r = 0.." << numRefine << "\n";
    gsInfo << "Output grid   : n = quA*p+quB = " << quA << "*" << degree << "+" << quB
           << " = " << n1d << " points/direction"
           << "  (mass integrand has tensor degree 2p = " << 2*degree
           << ", exactness needs n >= " << 2*degree+1 << ")\n";
    if (n1d < 2*degree + 1)
        gsInfo << "  *** WARNING: n < 2p+1 -- the moment-fitted mass matrix is "
                  "UNDER-INTEGRATED by construction;\n"
                  "      what follows then measures interpolation error, not "
                  "compression. ***\n";
    gsInfo << "Underlying    : gsAlgoimAdaptiveRule, indicator = " << indicator
           << ", maxDepth = " << quadDepth << "\n";
    gsInfo << "Moment fitting: compression fires whenever the underlying cloud "
              "exceeds n^d (pass-through guard), alpha = " << alpha << "\n";
    gsInfo << "Exact area    : " << std::setprecision(12) << EIGEN_PI*radius*radius
           << " (disk fully inside [0,1]^2 assumed)\n";
    gsInfo << "DoF retention : i kept iff M_adaptive(i,i) > 0 (alpha = 0 leaves the\n"
              "                rows of purely-outside DoFs exactly zero; their\n"
              "                eigenvalue sign would be roundoff).  Same index set\n"
              "                for both rules.\n\n";

    std::vector<index_t> hardLevels;         // resolved indefiniteness
    std::vector<index_t> marginalLevels;     // flagged only below the resolution of M
    std::vector<index_t> negWeightLevels;    // refinements where nNegativeWeights > 0
    real_t worstMassGap  = 0;                // max over levels of |sumW_ad - sumW_mf|
    real_t worstElemGap  = 0;                // max over levels/elements of the same, per element
    real_t worstAdMinWeight = std::numeric_limits<real_t>::max();  // min weight of the reference rule
    real_t worstAdLamMin    = 0;             // most negative RAW lam_min of a reference row
    real_t worstAdLamMax    = 0;             // its lam_max, for the eps*lam_max scale

    gsInfo << "  Columns: lam_min/lam_max/cond are RAW; lam_minJ/condJ are for the\n"
              "  Jacobi-scaled block D^-1/2 M D^-1/2, D = diag(M_adaptive) (same inertia).\n"
              "  minDiag is min_i M_ii over the retained set.\n"
              "  SPD:  yes = definite;  NO = indefinite and RESOLVED, i.e. either\n"
              "        minDiag < -eps*lam_max (a negative diagonal at a magnitude M can\n"
              "        represent) or lam_minJ < -10*eps*diagRange (a scaled eigenvalue\n"
              "        below what the Jacobi amplification of roundoff could fabricate);\n"
              "        no* = flagged, but only through DoFs whose entries lie AT OR BELOW\n"
              "        the numerical resolution of M -- see diagRange on the sub-line.\n"
              "  diagRange = max_i D_ii / min_i D_ii of the reference diagonal.  When\n"
              "  eps*diagRange >~ 1 the retained set contains numerically-null DoFs, the\n"
              "  Jacobi scaling amplifies assembly error by that factor, and the sign of\n"
              "  lam_minJ near those DoFs carries no information.\n"
              "  The minWeight column is filled for BOTH rules: the reference rule emits\n"
              "  only non-negative weights, which makes M_adaptive PSD by construction --\n"
              "  so any '(indef)' in a reference row is measured roundoff, not a defect,\n"
              "  and it calibrates the noise floor of the raw lam_min column.\n"
              "  LIMITATION: in the eps*diagRange >> 1 regime the scaled test cannot fire\n"
              "  (its threshold exceeds any eigenvalue of a unit-diagonal block), so\n"
              "  detection there rests on minDiag alone; an indefiniteness carried purely\n"
              "  by off-diagonal structure with a positive diagonal would be reported as\n"
              "  no*.  None of the configurations in the MEASURED block above is masked\n"
              "  this way, but the probe cannot rule it out in general.\n\n";
    gsInfo << "  r  mesh    dofs      rule      sum(w)          lam_min    lam_max    "
              "cond      lam_minJ   condJ     SPD  minDiag    nNeg  minWeight\n";
    gsInfo << "  --------------------------------------------------------------------"
              "-----------------------------------------------------------------\n";

    for (index_t r = 0; r <= numRefine; ++r)
    {
        const index_t nEl = (index_t)1 << (r+1);      // elements per direction

        gsKnotVector<real_t> kv(0.0, 1.0, nEl-1, degree+1);
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        const std::string phiStr = "sqrt((x-" + std::to_string(cx) + ")^2+(y-"
                                 + std::to_string(cy) + ")^2)-" + std::to_string(radius);
        gsFunctionExpr<real_t> phi(phiStr, 2);

        gsImplicitTrimmedDomain<2,real_t> domain(phi, basis);

        // Same options for both rules, so the underlying cut-cell rule of the
        // compressor is literally the same object type and configuration as the
        // reference rule.
        gsOptionList opt = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
        opt.setReal  ("quA", quA);
        opt.setInt   ("quB", quB);
        opt.setInt   ("maxDepth", quadDepth);
        opt.setString("indicator", indicator);

        gsAlgoimAdaptiveRule<real_t> adaptive(phi, basis, opt);

        // The compressor OWNS its adaptive rule (the gsQuadRule hierarchy has
        // no clone(), so ownership transfer is the only non-slicing option).
        gsMomentRule<real_t> momfit(
            gsQuadRule<real_t>::uPtr(new gsAlgoimAdaptiveRule<real_t>(phi, basis, opt)),
            gsVector<index_t>::Constant(2, n1d));

        gsSparseMatrix<real_t> Mad, Mmf;
        std::vector<real_t> elemAd, elemMf;

        real_t minWAd = 0, minWMf = 0;   // smallest weight EMITTED by each rule

        gsStopwatch clock;
        assembleMass(basis, domain, adaptive, Mad, elemAd, minWAd);
        const double tAd = clock.stop();

        clock.restart();
        assembleMass(basis, domain, momfit, Mmf, elemMf, minWMf);
        const double tMf = clock.stop();

        // Noise-floor calibration.  If the reference rule emitted no negative
        // weight, M_adaptive = sum_k w_k phi phi^T is PSD BY CONSTRUCTION, and
        // so is every principal submatrix of it.  Any negative lam_min printed
        // in an adaptive row is then, by definition, numerical error -- which
        // calibrates how much of a raw lam_min sign is meaningful at all.
        worstAdMinWeight = math::min(worstAdMinWeight, minWAd);

        // --- DoF retention, from the reference rule -------------------------
        const gsMatrix<real_t> denseAd = Mad.toDense();
        std::vector<index_t> keep;
        for (index_t i = 0; i != Mad.rows(); ++i)
            if (denseAd(i,i) > 0)
                keep.push_back(i);
        GISMO_ENSURE(!keep.empty(), "No DoF touches the physical domain.");

        // Jacobi scaling from the reference diagonal, shared by both matrices.
        gsVector<real_t> scale(keep.size());
        real_t dMin = std::numeric_limits<real_t>::max(), dMax = 0;
        for (size_t a = 0; a != keep.size(); ++a)
        {
            const real_t d = denseAd(keep[a], keep[a]);
            scale[a] = math::sqrt(d);
            dMin = math::min(dMin, d);
            dMax = math::max(dMax, d);
        }
        const real_t diagRange = dMax / dMin;

        clock.restart();
        const Spectrum sAd = spectrumOf(Mad, keep, scale);
        const Spectrum sMf = spectrumOf(Mmf, keep, scale);
        const double tEig = clock.stop();

        // --- quadrature mass ------------------------------------------------
        real_t sumAd = 0, sumMf = 0, elemGap = 0;
        index_t nElemGap = 0;
        GISMO_ENSURE(elemAd.size() == elemMf.size(),
                     "The two rules visited a different number of elements.");
        for (size_t e = 0; e != elemAd.size(); ++e)
        {
            sumAd += elemAd[e];
            sumMf += elemMf[e];
            const real_t g = math::abs(elemAd[e] - elemMf[e]);
            elemGap = math::max(elemGap, g);
            if (g > 1e-14) ++nElemGap;
        }
        worstMassGap = math::max(worstMassGap, math::abs(sumAd - sumMf));
        worstElemGap = math::max(worstElemGap, elemGap);

        const typename gsMomentRule<real_t>::Stats & st = momfit.stats();
        const bool compressed = (st.minWeight < std::numeric_limits<real_t>::max());
        if (st.nNegativeWeights > 0) negWeightLevels.push_back(r);

        // Definiteness verdict per rule.  A flag counts as RESOLVED only if it
        // survives the two resolution limits of this measurement: eps*lam_max
        // for the raw diagonal, and the Jacobi amplification eps*diagRange for
        // the scaled eigenvalue.  Below those, the sign says nothing -- see the
        // legend above.  The two tests are deliberately different quantities,
        // so a resolved flag never rests on one number alone.
        const real_t eps = std::numeric_limits<real_t>::epsilon();
        auto classify = [&](const Spectrum & s, bool & flagged) -> const char *
        {
            flagged = (s.minDiag <= 0) || (s.lambdaMinS <= 0);
            const bool resolved = (s.minDiag  < -eps * s.lambdaMax)
                               || (s.lambdaMinS < -10 * eps * diagRange);
            if (!flagged)  return "yes";
            return resolved ? "NO" : "no*";
        };
        bool flagAd = false, flagMf = false;
        const char * spdAd = classify(sAd, flagAd);
        const char * spdMf = classify(sMf, flagMf);
        if (sAd.lambdaMin < worstAdLamMin)
        {
            worstAdLamMin = sAd.lambdaMin;
            worstAdLamMax = sAd.lambdaMax;
        }
        if (flagMf)
        {
            if (0 == std::string("NO").compare(spdMf)) hardLevels.push_back(r);
            else                                       marginalLevels.push_back(r);
        }

        std::ostringstream mesh, dofs;
        mesh << nEl << "x" << nEl;
        dofs << keep.size() << "/" << Mad.rows();

        gsInfo << "  " << std::setw(1) << r
               << "  " << std::left << std::setw(7) << mesh.str()
               << " "  << std::setw(9) << dofs.str() << std::right
               << " adaptive "
               << " " << std::setw(15) << fmtSci(sumAd, 8)
               << " " << std::setw(10) << fmtSci(sAd.lambdaMin, 3)
               << " " << std::setw(10) << fmtSci(sAd.lambdaMax, 3)
               << " " << std::setw(9)  << fmtCond(sAd.lambdaMin, sAd.lambdaMax)
               << " " << std::setw(10) << fmtSci(sAd.lambdaMinS, 3)
               << " " << std::setw(9)  << fmtCond(sAd.lambdaMinS, sAd.lambdaMaxS)
               << " " << std::setw(4)  << spdAd
               << " " << std::setw(10) << fmtSci(sAd.minDiag, 3)
               << "     -  " << fmtSci(minWAd, 3) << "\n";

        gsInfo << "  " << std::setw(1) << " "
               << "  " << std::left << std::setw(7) << " "
               << " "  << std::setw(9) << " " << std::right
               << " momfit   "
               << " " << std::setw(15) << fmtSci(sumMf, 8)
               << " " << std::setw(10) << fmtSci(sMf.lambdaMin, 3)
               << " " << std::setw(10) << fmtSci(sMf.lambdaMax, 3)
               << " " << std::setw(9)  << fmtCond(sMf.lambdaMin, sMf.lambdaMax)
               << " " << std::setw(10) << fmtSci(sMf.lambdaMinS, 3)
               << " " << std::setw(9)  << fmtCond(sMf.lambdaMinS, sMf.lambdaMaxS)
               << " " << std::setw(4)  << spdMf
               << " " << std::setw(10) << fmtSci(sMf.minDiag, 3)
               << "  " << std::setw(4) << st.nNegativeWeights
               << "  " << (compressed ? fmtSci(st.minWeight, 3) : std::string("n/a"))
               << "\n";

        // "elems" counts EVERY element the rule visited (this probe loops
        // domain.beginAll()), not only the cut ones: it is larger than the old
        // "cut elems" figure by the inside+outside element count, by design --
        // the compressor carries no classification and cannot distinguish them.
        // "compressed" is nElements - nPassThroughElements; an element whose
        // wrapped rule returned nothing is counted as compressed-with-zero-
        // output, which the "output QPs" field already discriminates.
        gsInfo << "        elems " << st.nElements
               << ", compressed " << (st.nElements - st.nPassThroughElements)
               << ", output QPs " << st.nOutputQPs
               << " vs underlying " << st.nUnderlyingQPs
               << " | diagRange " << fmtSci(diagRange, 2)
               << " (eps*diagRange " << fmtSci(eps*diagRange, 2) << ")"
               << " | mass gap: total " << fmtSci(math::abs(sumAd-sumMf), 2)
               << ", max/elem " << fmtSci(elemGap, 2)
               << " (" << nElemGap << " elems > 1e-14)"
               << " | t_assemble " << std::fixed << std::setprecision(3)
               << tAd << "s / " << tMf << "s, t_eig " << tEig << "s\n";
    }

    // ------------------------------------------------------------------ verdict
    gsInfo << "\n";
    gsInfo << "Mass (partition-of-unity) check: max |sum(w)_adaptive - sum(w)_momfit| "
              "over the sweep = " << fmtSci(worstMassGap, 3)
           << " (per element: " << fmtSci(worstElemGap, 3) << ")\n";
    if (alpha > 0)
        gsInfo << "  alpha = " << alpha << " > 0: the two rules integrate DIFFERENT "
                  "measures (moment fitting\n"
                  "  adds alpha times the outside material), so this gap is expected to "
                  "be O(alpha)\n"
                  "  and is NOT a defect.  The identity only holds at alpha = 0.\n\n";
    else
        gsInfo << "  At alpha = 0 this identity is exact by construction; it is a sanity "
                  "check on the\n"
                  "  element classification, NOT evidence of accuracy.\n\n";

    gsInfo << "Noise floor: the reference rule's smallest emitted weight over the sweep is "
           << fmtSci(worstAdMinWeight, 3) << ".\n";
    if (worstAdMinWeight >= 0)
    {
        gsInfo << "  It is non-negative, so M_adaptive = sum_k w_k phi phi^T is PSD by "
                  "construction,\n"
                  "  and so is every principal submatrix of it.\n";
        if (worstAdLamMin < 0)
            gsInfo << "  Hence the most negative RAW lam_min printed in a reference row, "
                   << fmtSci(worstAdLamMin, 3) << ",\n"
                      "  is by definition pure numerical error -- and it sits right at "
                      "eps*lam_max = "
                   << fmtSci(std::numeric_limits<real_t>::epsilon() * worstAdLamMax, 3)
                   << ".\n  That is the measured scale at which a raw lam_min sign stops "
                      "meaning anything,\n  for EITHER rule; it is why the verdict uses "
                      "the resolution tests above.\n\n";
        else
            gsInfo << "  No reference row went negative in this sweep, so this run gives "
                      "no direct\n  calibration of the raw noise floor (other "
                      "configurations do).\n\n";
    }
    else
        gsInfo << "  It is NEGATIVE -- unexpected for the reference rule; the PSD-by-"
                  "construction\n  calibration above does not apply to this run.\n\n";

    std::ostringstream cfg;
    cfg << "p=" << degree << ", n=" << n1d << ", maxDepth=" << quadDepth
        << ", indicator=" << indicator << ", alpha=" << alpha
        << ", R=" << radius << ", centre=(" << cx << "," << cy << "), r=0.."
        << numRefine;

    if (negWeightLevels.empty())
    {
        gsInfo << "NOTE: no negative moment-fitted weight occurred in this sweep -- the\n"
                  "      probe did NOT exercise the phenomenon it exists to measure.\n"
                  "      Move the disk off the dyadic mesh symmetry (--cx/--cy/--radius)\n"
                  "      or lower the quadrature order.\n";
    }
    else
    {
        gsInfo << "Negative moment-fitted weights occurred at refinement(s): ";
        for (size_t i = 0; i != negWeightLevels.size(); ++i)
            gsInfo << negWeightLevels[i] << (i+1 < negWeightLevels.size() ? ", " : "");
        gsInfo << "\n";
    }

    auto listOf = [](const std::vector<index_t> & v)
    {
        std::ostringstream os;
        for (size_t i = 0; i != v.size(); ++i)
            os << v[i] << (i+1 < v.size() ? ", " : "");
        return os.str();
    };

    if (!marginalLevels.empty())
        gsInfo << "MARGINAL (no*): at refinement(s) " << listOf(marginalLevels)
               << " the moment-fitted block was flagged only through DoFs whose\n"
                  "  entries lie at or below the numerical resolution of M (see "
                  "diagRange).  Such a\n"
                  "  flag is NOT stable: it can flip with the depth of the underlying "
                  "rule.  It is\n"
                  "  reported, not counted as a verdict.\n";

    if (hardLevels.empty())
    {
        gsInfo << "VERDICT: moment fitting produced NO resolved indefinite mass matrix at "
                  "any refinement of this sweep (" << cfg.str() << "): every negative "
                  "eigenvalue/diagonal seen stayed below the numerical resolution of the "
                  "matrix, despite the negative weights reported above.\n";
    }
    else
    {
        gsInfo << "VERDICT: moment fitting produced an INDEFINITE mass matrix at "
                  "refinement(s) " << listOf(hardLevels)
               << " of this sweep (" << cfg.str() << "): a negative diagonal entry and/or "
                  "a negative scaled eigenvalue beyond the resolution limits there.\n";
    }

    return EXIT_SUCCESS;
}
