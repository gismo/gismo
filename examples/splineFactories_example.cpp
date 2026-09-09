/** @file splineFactories_example.cpp

    @brief Compares spline factory methods (toBezier, squared, cubed, grad, lapl)
    against L2 projection, quasi-interpolation, and interpolateAtAnchors for a
    2D scalar tensor B-spline.

    For each derived quantity the factory result serves as the reference.
    All three alternative methods project / interpolate the same derived function
    onto the factory basis and should reproduce the factory coefficients to
    floating-point precision.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include <gismo.h>

using namespace gismo;

// =============================================================================
// Function wrappers: pointwise evaluation of derived quantities from a 2D
// scalar B-spline c. Required because gsExprAssembler::getCoeff calls clone().
// =============================================================================

template <class T>
class SquaredSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
public:
    explicit SquaredSplineFunction(const gsTensorBSpline<2,T>& c) : m_c(c) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        m_c.eval_into(u, result);                         // 1 × npts
        result = result.cwiseProduct(result);
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(SquaredSplineFunction)
};

template <class T>
class CubedSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
public:
    explicit CubedSplineFunction(const gsTensorBSpline<2,T>& c) : m_c(c) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> val;
        m_c.eval_into(u, val);
        result = val.cwiseProduct(val.cwiseProduct(val));  // c³ element-wise
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(CubedSplineFunction)
};

// ∂c/∂x_{dir} wrapper; dir ∈ {0,1}
template <class T>
class GradSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
    short_t m_dir;
public:
    GradSplineFunction(const gsTensorBSpline<2,T>& c, short_t dir)
        : m_c(c), m_dir(dir) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> deriv;
        m_c.deriv_into(u, deriv);   // (domainDim*targetDim) × npts = 2 × npts
        result = deriv.row(m_dir);  // 1 × npts
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(GradSplineFunction)
};

// Δc = ∂²c/∂x² + ∂²c/∂y² (second-order Voigt rows 0 and 1)
template <class T>
class LaplSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
public:
    explicit LaplSplineFunction(const gsTensorBSpline<2,T>& c) : m_c(c) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> d2;
        m_c.deriv2_into(u, d2);     // 3 × npts for scalar 2D: [∂²/∂x², ∂²/∂y², ∂²/∂x∂y]
        result = d2.row(0) + d2.row(1);
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(LaplSplineFunction)
};

// =============================================================================
// Helpers
// =============================================================================

static gsTensorBSpline<2,real_t> makeIdentityGeo()
{
    gsKnotVector<real_t> kv(0.0, 1.0, 0, 2); // [0 0 1 1]
    gsMatrix<real_t> coefs(4, 2);
    coefs << 0,0,  1,0,  0,1,  1,1;
    return gsTensorBSpline<2,real_t>(kv, kv, coefs);
}

template <class T>
void printRow(const std::string& qty, const std::string& method,
              const gsMatrix<T>& ref, const gsMatrix<T>& other, double secs)
{
    gsMatrix<T> diff = ref - other;
    gsInfo << "  " << std::left
           << std::setw(10) << qty
           << std::setw(22) << method
           << std::scientific << std::setprecision(3)
           << std::setw(14) << diff.norm()
           << std::setw(14) << diff.template lpNorm<gsEigen::Infinity>()
           << secs << "\n";
}

template <class T>
void printRef(const std::string& qty, double secs)
{
    gsInfo << "  " << std::left
           << std::setw(10) << qty
           << std::setw(22) << "Factory (ref)"
           << std::setw(14) << "(exact)"
           << std::setw(14) << "(exact)"
           << secs << "\n";
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[])
{
    index_t degree   = 3;
    index_t nRefine  = 2;
    index_t numReps  = 5;

    gsCmdLine cmd("Compare B-spline factory methods (squared, cubed, grad, lapl) "
                  "against L2 projection, quasi-interpolation, and interpolateAtAnchors.");
    cmd.addInt("p","degree",  "B-spline degree",              degree);
    cmd.addInt("r","refine",  "Number of uniform refinements", nRefine);
    cmd.addInt("n","nreps",   "Repetitions for timing",        numReps);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // -------------------------------------------------------------------------
    // 1.  Build a 2D scalar B-spline c approximating sin(πx)·sin(πy)
    // -------------------------------------------------------------------------
    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);  // 1 interior knot at 0.5
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);

    // Interpolate f at Greville points → spline c
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);                       // 1 × n_bfns
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsInfo << "Spline c: degree=" << degree
           << ", dofs=" << basis.size()
           << " (" << basis.size(0) << "×" << basis.size(1) << ")\n\n";

    // -------------------------------------------------------------------------
    // 2.  Identity geometry map (bilinear parametrisation of [0,1]²)
    // -------------------------------------------------------------------------
    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    // -------------------------------------------------------------------------
    // 3.  Stopwatch and output setup
    // -------------------------------------------------------------------------
    gsStopwatch timer;

    gsInfo << std::left
           << "  " << std::setw(10) << "Quantity"
           << std::setw(22) << "Method"
           << std::setw(14) << "L2 coef err"
           << std::setw(14) << "Linf coef err"
           << "Time (s)\n";
    gsInfo << std::string(80, '-') << "\n";

    // =========================================================================
    // toBezier
    // =========================================================================
    {
        timer.restart();
        for (index_t k = 0; k < numReps; ++k) { auto r = c.toBezier(); (void)r; }
        double tFact = timer.stop() / numReps;
        auto bz = c.toBezier();

        // Functional verification: sample 100 points, check max |c(x)-bz(x)|
        gsMatrix<real_t> pts(2, 100);
        pts.setRandom();
        pts = (pts.array() * 0.5 + 0.5).matrix();  // map (-1,1) → (0,1)
        gsMatrix<real_t> vc, vbz;
        c.eval_into(pts, vc);
        bz.eval_into(pts, vbz);
        real_t funcErr = (vc - vbz).norm();

        gsInfo << "  " << std::left
               << std::setw(10) << "toBezier"
               << std::setw(22) << "Factory (ref)"
               << std::setw(14) << "(func check)"
               << std::scientific << std::setprecision(3)
               << std::setw(14) << funcErr
               << tFact << "  (func err at 100 pts)\n";
    }

    // =========================================================================
    // squared  c² → basis degree 2p, reduced continuity
    // =========================================================================
    {
        // --- Factory (reference) ---
        timer.restart();
        for (index_t k = 0; k < numReps; ++k) { auto r = c.squared(); (void)r; }
        double tFact = timer.stop() / numReps;
        auto c2 = c.squared();
        printRef<real_t>("squared", tFact);

        const gsTensorBSplineBasis<2,real_t>& bSq = c2.basis();
        const gsMatrix<real_t>& refCoefs = c2.coefs();

        SquaredSplineFunction<real_t> sq(c);

        // --- L2 projection ---
        gsMatrix<real_t> coefs_l2;
        timer.restart();
        gsL2Projection<real_t>::project(bSq, idGeo, sq, coefs_l2);
        double tL2 = timer.stop();
        printRow<real_t>("squared", "L2 projection", refCoefs, coefs_l2, tL2);

        // --- Quasi-interpolation ---
        gsMatrix<real_t> coefs_qi;
        timer.restart();
        gsQuasiInterpolate<real_t>::localIntpl(bSq, sq, coefs_qi);
        double tQI = timer.stop();
        printRow<real_t>("squared", "Quasi-interp", refCoefs, coefs_qi, tQI);

        // --- interpolateAtAnchors ---
        timer.restart();
        gsMatrix<real_t> anch = bSq.anchors();     // 2 × n_bfns
        gsMatrix<real_t> vals;
        sq.eval_into(anch, vals);                   // 1 × n_bfns
        auto geom_ia = bSq.interpolateAtAnchors(vals);
        double tIA = timer.stop();
        printRow<real_t>("squared", "InterpAtAnchors", refCoefs, geom_ia->coefs(), tIA);
    }

    // =========================================================================
    // cubed  c³ → basis degree 3p, reduced continuity
    // =========================================================================
    {
        timer.restart();
        for (index_t k = 0; k < numReps; ++k) { auto r = c.cubed(); (void)r; }
        double tFact = timer.stop() / numReps;
        auto c3 = c.cubed();
        printRef<real_t>("cubed", tFact);

        const gsTensorBSplineBasis<2,real_t>& bCu = c3.basis();
        const gsMatrix<real_t>& refCoefs = c3.coefs();

        CubedSplineFunction<real_t> cu(c);

        gsMatrix<real_t> coefs_l2;
        timer.restart();
        gsL2Projection<real_t>::project(bCu, idGeo, cu, coefs_l2);
        printRow<real_t>("cubed", "L2 projection", refCoefs, coefs_l2, timer.stop());

        gsMatrix<real_t> coefs_qi;
        timer.restart();
        gsQuasiInterpolate<real_t>::localIntpl(bCu, cu, coefs_qi);
        printRow<real_t>("cubed", "Quasi-interp", refCoefs, coefs_qi, timer.stop());

        timer.restart();
        gsMatrix<real_t> anch = bCu.anchors();
        gsMatrix<real_t> vals;
        cu.eval_into(anch, vals);
        auto geom_ia = bCu.interpolateAtAnchors(vals);
        printRow<real_t>("cubed", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    }

    // =========================================================================
    // grad x  ∂c/∂x → degree (p-1, p)
    // =========================================================================
    {
        timer.restart();
        for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(0); (void)r; }
        double tFact = timer.stop() / numReps;
        auto cdx = c.grad(0);
        printRef<real_t>("grad_x", tFact);

        const gsTensorBSplineBasis<2,real_t>& bGx = cdx.basis();
        const gsMatrix<real_t>& refCoefs = cdx.coefs();

        GradSplineFunction<real_t> gx(c, 0);

        gsMatrix<real_t> coefs_l2;
        timer.restart();
        gsL2Projection<real_t>::project(bGx, idGeo, gx, coefs_l2);
        printRow<real_t>("grad_x", "L2 projection", refCoefs, coefs_l2, timer.stop());

        gsMatrix<real_t> coefs_qi;
        timer.restart();
        gsQuasiInterpolate<real_t>::localIntpl(bGx, gx, coefs_qi);
        printRow<real_t>("grad_x", "Quasi-interp", refCoefs, coefs_qi, timer.stop());

        timer.restart();
        gsMatrix<real_t> anch = bGx.anchors();
        gsMatrix<real_t> vals;
        gx.eval_into(anch, vals);
        auto geom_ia = bGx.interpolateAtAnchors(vals);
        printRow<real_t>("grad_x", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    }

    // =========================================================================
    // grad y  ∂c/∂y → degree (p, p-1)
    // =========================================================================
    {
        timer.restart();
        for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(1); (void)r; }
        double tFact = timer.stop() / numReps;
        auto cdy = c.grad(1);
        printRef<real_t>("grad_y", tFact);

        const gsTensorBSplineBasis<2,real_t>& bGy = cdy.basis();
        const gsMatrix<real_t>& refCoefs = cdy.coefs();

        GradSplineFunction<real_t> gy(c, 1);

        gsMatrix<real_t> coefs_l2;
        timer.restart();
        gsL2Projection<real_t>::project(bGy, idGeo, gy, coefs_l2);
        printRow<real_t>("grad_y", "L2 projection", refCoefs, coefs_l2, timer.stop());

        gsMatrix<real_t> coefs_qi;
        timer.restart();
        gsQuasiInterpolate<real_t>::localIntpl(bGy, gy, coefs_qi);
        printRow<real_t>("grad_y", "Quasi-interp", refCoefs, coefs_qi, timer.stop());

        timer.restart();
        gsMatrix<real_t> anch = bGy.anchors();
        gsMatrix<real_t> vals;
        gy.eval_into(anch, vals);
        auto geom_ia = bGy.interpolateAtAnchors(vals);
        printRow<real_t>("grad_y", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    }

    // =========================================================================
    // lapl  Δc = ∂²c/∂x² + ∂²c/∂y² → degree (p-2, p-2) in each direction
    // =========================================================================
    if (degree >= 2)
    {
        timer.restart();
        for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
        double tFact = timer.stop() / numReps;
        auto cl = c.lapl();
        printRef<real_t>("lapl", tFact);

        const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
        const gsMatrix<real_t>& refCoefs = cl.coefs();

        LaplSplineFunction<real_t> lp(c);

        gsMatrix<real_t> coefs_l2;
        timer.restart();
        gsL2Projection<real_t>::project(bL, idGeo, lp, coefs_l2);
        printRow<real_t>("lapl", "L2 projection", refCoefs, coefs_l2, timer.stop());

        gsMatrix<real_t> coefs_qi;
        timer.restart();
        gsQuasiInterpolate<real_t>::localIntpl(bL, lp, coefs_qi);
        printRow<real_t>("lapl", "Quasi-interp", refCoefs, coefs_qi, timer.stop());

        timer.restart();
        gsMatrix<real_t> anch = bL.anchors();
        gsMatrix<real_t> vals;
        lp.eval_into(anch, vals);
        auto geom_ia = bL.interpolateAtAnchors(vals);
        printRow<real_t>("lapl", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    }
    else
    {
        gsInfo << "  (lapl skipped: requires degree >= 2)\n";
    }

    gsInfo << std::string(80, '-') << "\n";
    gsInfo << "Note: all three comparison methods should reproduce factory coefficients\n"
           << "to near floating-point precision, since each derived quantity is exactly\n"
           << "representable in the factory-produced basis.\n";

    return EXIT_SUCCESS;
}
