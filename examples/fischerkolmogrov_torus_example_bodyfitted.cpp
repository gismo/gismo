/** @file fischerkolmogrov_torus_example_bodyfitted.cpp

    @brief Fisher-Kolmogorov tumour growth on a body-fitted SOLID TORUS
           (genus 1), as the reference for the immersed/body-fitted comparison.

    The solver (generalized-alpha time stepping, Newton loop, tanh tumour seed,
    biomarkers, csv output) is exactly the one of
    fischerkolmogrov_prostate_example_bodyfitted.cpp. Only the geometry differs:
    instead of reading a patch from the xml file, an EXACT solid torus is built
    as a single NURBS patch from gismo's own primitives,

        torus = gsNurbsCreator::NurbsDisk  (cross-section, minor radius r)
                revolved along
                gsNurbsCreator::NurbsCircle (toroidal direction, major radius R)

    i.e. with the disk control points (a_ij, b_ij), weights w_ij and the unit
    circle control points (cx_k, cy_k), weights v_k:

        P_ijk = c + ( (R + r a_ij) cx_k , (R + r a_ij) cy_k , r b_ij ),
        weight  w_ij * v_k.

    Parametric directions:  u, v : disk cross-section,  w : toroidal angle.
    The toroidal seam (w = 0 and w = 1 are the same surface) is glued as a
    patch interface, so the discrete solution is continuous across it.
    The disk patch has four corner points with a degenerate (180 deg) Jacobian;
    they are isolated points and never quadrature points.

    Default torus: centre (0.5,0.5,0.5), R = 0.3, r = 0.15, so it spans
    [0.05,0.95]^2 x [0.35,0.65] inside [0,1]^3. Its largest extent is 0.9, which
    is also what gsNormalizeToUnitBox (fill 0.9) gives in the immersed example,
    so a triangulation of this torus is left unchanged by that normalization.
    Exact volume: 2 pi^2 R r^2 (printed for checking V_P).

    UNITS. All model input is PHYSICAL: lengths in mm, time in days, D in
    mm^2/day, rho in 1/day. The torus lives in [0,1]^3 and one box unit is
    L mm (-L, default 50: outer diameter 45 mm, tube radius 7.5 mm, prostate
    sized). The conversion to box units (D/L^2, R/L, ...) is done once, in
    fkToBox(); the csv reports volumes in mm^3. Before the time loop the
    program prints the dimensionless numbers that decide whether the tumour
    grows or is diluted by diffusion (rho*T, D/R^2 vs rho, front width vs h).

    The bounding-box centre is in the HOLE of the torus, so the default tumour
    centre is the point on the core circle at +x. The seed must fit in the
    tube: R_seed + 2 sigma < r*L. The seed peak is -a (default 0.5) so that the
    core can get denser over time (a saturated seed, a = 1, cannot).

    Literature-order prostate values, two years between MRI scans
    (D = 1e-3 mm^2/day, rho = 5e-3 1/day, T = 146 x 5 = 730 days):

    ./build/bin/fischerkolmogrov_torus_example_bodyfitted \
        -D 1e-3 -p 5e-3 -t 5 -N 146 -R 5 -s 1 -a 0.5 -r 3 -L 50 \
        -o 00_output_FK_torus --plot

    For a slightly displaced intial condition (not on the torus' axis)
    ./build/bin/fischerkolmogrov_torus_example_bodyfitted \
        -D 1e-3 -p 5e-3 -t 5 -N 146 -R 3 -s 0.75 -a 0.5 -r 4 -L 50 \
        --xc 7.88 --yc 25 --zc 27.12 \
        -o 00_output_FK_torus --plot

    rm -rf 00_output_FK_torus
    ./build/bin/fischerkolmogrov_torus_example_bodyfitted \
        -D 1e-3 -p 5e-3 -t 5 -N 146 -R 3 -s 0.75 -a 0.5 \
        -r 4 --refTor 1 --plotPts 200000 -L 50 \
        --xc 7.88 --yc 25 --zc 27.12 \
        -o 00_output_FK_torus --plot

    The xml file (-f) is only used for rho_inf and the Newton tolerance; its
    model parameters (id 3, in box units) are ignored on purpose.

    Author(s): L. Venta-Viñuela
*/
#include <gismo.h>

using namespace gismo;

// =============================================================================
// LINEAR SOLVER (as in adaptive_phase-field_fracture_example_DC_*.cpp)
//   "cgdiag": CG with a diagonal (Jacobi) preconditioner (default).
//   "cgic":   CG with an incomplete-Cholesky preconditioner.
//   "direct": Eigen SimplicialLDLT (Pardiso LDLT if G+Smo has Pardiso). Exact,
//             but its time and memory grow superlinearly with the DoFs in 3D.
// The FK tangent alpha_m M + alpha_f gamma dt (D K - rho(1-2N) M) is mass
// dominated for small dt; the iteration count (LinIts in computational.csv)
// shows how well conditioned it is (see also --cond).
// =============================================================================
template <class Solver, class T>
index_t fkSolveCG(Solver & solver, const char * name,
                  const gsSparseMatrix<T> & K, const gsMatrix<T> & b, gsMatrix<T> & x,
                  index_t maxIt, T tol)
{
    solver.setMaxIterations(maxIt);
    solver.setTolerance(tol);
    solver.compute(K);
    GISMO_ENSURE(solver.info()==gsEigen::Success,
                 name<<": the preconditioner setup failed.");
    x = solver.solve(b);
    GISMO_ENSURE(solver.info()==gsEigen::Success,
                 name<<" did not converge: "<<solver.iterations()<<"/"<<maxIt<<" iteration(s), "
                 "error = "<<solver.error()<<" > "<<tol<<". Raise --cgMaxIt or --cgTol, or use "
                 "--solver cgic / direct.");
    return solver.iterations();
}

/// Solves K x = b with \a linSolver; returns the number of iterations (1 for direct).
template <class T>
index_t fkSolveLinear(const std::string & linSolver,
                      const gsSparseMatrix<T> & K, const gsMatrix<T> & b, gsMatrix<T> & x,
                      index_t maxIt, T tol)
{
    typedef gsEigen::SparseMatrix<T,0,index_t> EigenMat;
    if (linSolver=="cgic")
    {
        gsEigen::ConjugateGradient<EigenMat, gsEigen::Lower|gsEigen::Upper,
                                   gsEigen::IncompleteCholesky<T,gsEigen::Lower,
                                            gsEigen::AMDOrdering<index_t> > > solver;
        return fkSolveCG(solver,"CG (incomplete Cholesky)",K,b,x,maxIt,tol);
    }
    else if (linSolver=="cgdiag")
    {
        gsEigen::ConjugateGradient<EigenMat, gsEigen::Lower|gsEigen::Upper,
                                   gsEigen::DiagonalPreconditioner<T> > solver;
        return fkSolveCG(solver,"CG (diagonal)",K,b,x,maxIt,tol);
    }
    GISMO_ENSURE(linSolver=="direct","Unknown linear solver \""<<linSolver<<"\" (cgdiag|cgic|direct)");
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
    typename gsSparseSolver<T>::SimplicialLDLT solver;
#endif
    solver.compute(K);
    GISMO_ENSURE(solver.info()==gsEigen::Success,"Factorisation of the FK matrix failed.");
    x = solver.solve(b);
    return 1;
}

// =============================================================================
// CONDITION NUMBER (--cond): spectral condition number kappa = lmax/lmin of a
// symmetric positive definite matrix, estimated with the power iteration (lmax)
// and the inverse power iteration with an LDLT factorization (lmin). Done for
// the system matrix K and for the Jacobi-scaled S = D^-1/2 K D^-1/2, which is
// the matrix CG with the diagonal preconditioner effectively sees.
// =============================================================================
template <class T>
struct fkCond { T lmax, lmin, kappa; index_t itMax, itMin; };

template <class T>
fkCond<T> fkConditionNumber(const gsSparseMatrix<T> & K, index_t maxIt = 1000, T tol = 1e-8)
{
    fkCond<T> c;
    const index_t n = K.rows();
    gsMatrix<T> x = gsMatrix<T>::Ones(n,1), y;
    for (index_t i = 0; i != n; ++i) x(i,0) += 0.1 * math::sin(1.0 + i); // generic start
    x /= x.norm();
    // power iteration: largest eigenvalue
    T lam = 0, lamOld = 0;
    for (c.itMax = 1; c.itMax <= maxIt; ++c.itMax)
    {
        y = K * x;
        lam = x.col(0).dot(y.col(0)); // Rayleigh quotient
        x = y / y.norm();
        if (c.itMax > 1 && math::abs(lam - lamOld) < tol * math::abs(lam)) break;
        lamOld = lam;
    }
    c.lmax = lam;
    // inverse power iteration: smallest eigenvalue
    typename gsSparseSolver<T>::SimplicialLDLT ldlt;
    ldlt.compute(K);
    GISMO_ENSURE(ldlt.info()==gsEigen::Success, "--cond: factorization failed.");
    x = gsMatrix<T>::Ones(n,1);
    for (index_t i = 0; i != n; ++i) x(i,0) += 0.1 * math::cos(2.0 + i);
    x /= x.norm();
    lamOld = 0;
    for (c.itMin = 1; c.itMin <= maxIt; ++c.itMin)
    {
        y = ldlt.solve(x);
        lam = x.col(0).dot(y.col(0)); // = 1/lmin estimate
        x = y / y.norm();
        if (c.itMin > 1 && math::abs(lam - lamOld) < tol * math::abs(lam)) break;
        lamOld = lam;
    }
    c.lmin  = 1.0 / lam;
    c.kappa = c.lmax / c.lmin;
    return c;
}

/// Jacobi-scaled copy D^-1/2 K D^-1/2 (unit diagonal)
template <class T>
gsSparseMatrix<T> fkJacobiScaled(const gsSparseMatrix<T> & K)
{
    gsVector<T> d = K.diagonal().cwiseAbs().cwiseSqrt().cwiseInverse();
    gsSparseMatrix<T> S = K;
    for (index_t k = 0; k < S.outerSize(); ++k)
        for (typename gsSparseMatrix<T>::InnerIterator it(S, k); it; ++it)
            it.valueRef() *= d[it.row()] * d[it.col()];
    return S;
}

/// Writes the --cond results to the log and returns them as csv lines
template <class T>
std::string fkReportCond(const gsSparseMatrix<T> & K)
{
    const fkCond<T> c = fkConditionNumber(K);
    const fkCond<T> s = fkConditionNumber(fkJacobiScaled(K));
    gsInfo << "Condition number of the system matrix (first Newton matrix, n = " << K.rows() << "):\n"
           << "  K             : lmax = " << c.lmax << ", lmin = " << c.lmin << ", kappa = " << c.kappa
           << "  (" << c.itMax << "/" << c.itMin << " power/inverse iterations)\n"
           << "  D^-1/2 K D^-1/2: lmax = " << s.lmax << ", lmin = " << s.lmin << ", kappa = " << s.kappa
           << "  (" << s.itMax << "/" << s.itMin << " power/inverse iterations)\n";
    std::ostringstream os;
    os.precision(10);
    os << "cond_K," << c.kappa << "\ncond_K_lmax," << c.lmax << "\ncond_K_lmin," << c.lmin << "\n"
       << "cond_jacobi," << s.kappa << "\ncond_jacobi_lmax," << s.lmax << "\ncond_jacobi_lmin," << s.lmin << "\n";
    return os.str();
}

// =============================================================================
// GEOMETRY: exact single-patch NURBS solid torus
// =============================================================================
struct torusParams
{
    real_t R  = 0.30; // major radius (centre of the tube to the torus axis)
    real_t r  = 0.15; // minor radius (tube radius)
    real_t cx = 0.5, cy = 0.5, cz = 0.5; // centre; the axis is parallel to z
};

// periodic = true : no seam interface; the solution basis is made periodic in
//                   the toroidal direction w instead (C1 across w = 0 = 1).
// periodic = false: the seam faces w = 0 and w = 1 are glued as a C0 interface.
template <class T>
gsMultiPatch<T> makeTorus(const torusParams & tp, bool periodic = true)
{
    GISMO_ENSURE(tp.r > 0 && tp.R > tp.r, "Need R > r > 0 for a torus.");

    // Cross-section: exact unit disk (3x3 biquadratic NURBS), and the exact
    // unit circle (9 control points, quadratic) for the toroidal direction.
    typename gsNurbsCreator<T>::TensorNurbs2Ptr disk = gsNurbsCreator<T>::NurbsDisk(1, 0, 0);
    typename gsNurbsCreator<T>::NurbsPtr        circ = gsNurbsCreator<T>::NurbsCircle(1, 0, 0);

    const gsMatrix<T> & Cd = disk->coefs();   // nd x 2
    const gsMatrix<T> & Wd = disk->weights(); // nd x 1
    const gsMatrix<T> & Cc = circ->coefs();   // nc x 2
    const gsMatrix<T> & Wc = circ->weights(); // nc x 1
    const index_t nd = Cd.rows(), nc = Cc.rows();

    // Tensor ordering: the first direction runs fastest, so the disk index
    // (itself u-fastest) runs inside the toroidal index k.
    gsMatrix<T> coefs(nd * nc, 3), wts(nd * nc, 1);
    for (index_t k = 0; k != nc; ++k)
        for (index_t i = 0; i != nd; ++i)
        {
            const index_t g   = i + nd * k;
            const T       rho = tp.R + tp.r * Cd(i, 0); // distance to the axis
            coefs(g, 0) = tp.cx + rho * Cc(k, 0);
            coefs(g, 1) = tp.cy + rho * Cc(k, 1);
            coefs(g, 2) = tp.cz + tp.r * Cd(i, 1);
            wts(g, 0)   = Wd(i, 0) * Wc(k, 0);
        }

    gsTensorNurbs<3, T> torus(disk->knots(0), disk->knots(1), circ->knots(),
                              give(coefs), give(wts));

    gsMultiPatch<T> mp;
    mp.addPatch(torus);
    if (!periodic)
        mp.addInterface(0, boundary::front, 0, boundary::back); // toroidal seam w=0 <-> w=1
    // The 4 disk sides (west/east/south/north) are the torus surface. In the
    // periodic case front/back (w = 0, 1) are also listed as boundaries, but no
    // condition is imposed anywhere (natural zero flux) and the periodic basis
    // has no functions that live only on them, so they carry no constraint.
    mp.addAutoBoundaries();
    return mp;
}

// =============================================================================
// FISHER-KOLMOGOROV MODEL PARAMETERS (user-settable, for literature calibration)
// =============================================================================
// ALL USER INPUT IS IN PHYSICAL UNITS: lengths in mm, time in days.
// The geometry (and the solver) live in the box [0,1]^3; one box unit is L mm.
// The conversion to box units is done in ONE place, fkToBox(), so a literature
// value (e.g. D = 1e-3 mm^2/day) can be typed in as it is.
//   D_box     = D / L^2      [box^2/day]
//   R_box     = R / L,  sigma_box = sigma / L,  centre_box = centre / L
//   rho, dt   unchanged      [1/day], [day]
// Defaults: prostate cancer values of the order reported in the literature
// (D ~ 1e-3 mm^2/day, rho ~ 5e-3 1/day), with T = dt*Nsteps ~ 1-2 years.
template <class T>
struct fkParams
{
    T D     = 1.0e-3; // tumor cell diffusivity          [mm^2/day]
    T rho   = 5.0e-3; // net proliferation rate          [1/day]
    T R     = 5.0;    // radius of the initial tumor     [mm]
    T sigma = 1.0;    // interface width of the IC       [mm]
    T amp   = 0.5;    // peak density of the IC, N0 <= amp <= 1 [-]
    T Nth   = 0.15;   // threshold defining Omega_T      [-]
    T L     = 50.0;   // length of one box unit          [mm]
    // Centre of the initial tumor [mm]. NaN: a point on the torus core circle.
    T xc = std::numeric_limits<T>::quiet_NaN();
    T yc = std::numeric_limits<T>::quiet_NaN();
    T zc = std::numeric_limits<T>::quiet_NaN();
};

// Physical (mm, day) -> solver (box, day) units.
template <class T>
fkParams<T> fkToBox(const fkParams<T> & p)
{
    fkParams<T> b = p;
    b.D     = p.D / (p.L * p.L);
    b.R     = p.R / p.L;
    b.sigma = p.sigma / p.L;
    b.xc    = p.xc / p.L;
    b.yc    = p.yc / p.L;
    b.zc    = p.zc / p.L;
    b.L     = 1;
    return b;
}

// Smooth spherical tumor seed:  N(x,0) = amp*0.5*(1 - tanh((||x - x_tumor|| - R)/sigma))
// (expects BOX units, i.e. the output of fkToBox)
template <short_t dim, class T>
gsFunctionExpr<T> tumorSeed(const fkParams<T> & par, const gsMultiPatch<T> & mp)
{
    // Default seed location: centre of the bounding box of the domain
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    gsVector<T,3> c;
    c.setZero();
    for (short_t d = 0; d != dim; d++)
        c[d] = 0.5 * (bbox(d,0) + bbox(d,1));

    if (!math::isnan(par.xc))            c[0] = par.xc;
    if (!math::isnan(par.yc))            c[1] = par.yc;
    if (dim == 3 && !math::isnan(par.zc)) c[2] = par.zc;

    std::ostringstream expr;
    expr.precision(16);
    expr << par.amp << "*0.5*(1-tanh((sqrt((x-(" << c[0] << "))^2+(y-(" << c[1] << "))^2";
    if (dim == 3)
        expr << "+(z-(" << c[2] << "))^2";
    expr << ")-(" << par.R << "))/(" << par.sigma << ")))";

    gsInfo << "Initial condition: N(x,0) = " << expr.str() << "\n";
    return gsFunctionExpr<T>(expr.str(), dim);
}

template <class T, class Expr>
T computeInt(const gsMultiBasis<T>  & integrationBasis,
                const gsFunctionSet<T> & geometryMap,
                const Expr & sourceFunction1,
                const gsOptionList     & options)
{
  
    // Create an assembler
    gsExprAssembler<T> A(1,1);
    A.options().update(options,gsOptionList::addIfUnknown); 

    // Set the integration elements
    A.setIntegrationElements(integrationBasis);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    gsExprEvaluator<T> ev(A);
    ev.options().update(A.options(),gsOptionList::addIfUnknown);
    ev.options().setSwitch("SameElement",false); // add warning to remind that im deactivating sameelemnt        

    return ev.integral(sourceFunction1 * meas(G));

}

// =============================================================================
// SOLVER ROUTINE
// =============================================================================
template <short_t dim, class T>
void solve(gsMultiPatch<T> & mp,
           gsBoundaryConditions<T> & bc,
           gsOptionList & CHopt,
           gsOptionList & TIMEopt,
           gsOptionList & Aopt,
           real_t & dt,
           index_t & maxSteps,
           index_t & numRefine,
           const fkParams<T> & parPhys,
           std::string out,
           bool plot,
           index_t refTor  = 0,      // extra refinements in the toroidal direction w only
           index_t plotPts = 50000,  // ParaView sample points per frame
           real_t volExact = 0,      // exact |Omega| in mm^3, for the csv
           bool periodic = true,     // C1 (periodic) instead of C0 (glued) at the seam
           bool plotMesh = false,    // write the computational mesh to ParaView
           const std::string & linSolver = "cgdiag", // cgdiag | cgic | direct
           index_t cgMaxIt = 20000,
           real_t  cgTol   = 1e-8,   // relative to ||rhs||
           bool    cond    = false)  // estimate the condition number of the first matrix
{
    // Everything below runs in BOX units (lengths / L); the csv converts back.
    const fkParams<T> par = fkToBox(parPhys);
    const real_t L3 = parPhys.L * parPhys.L * parPhys.L; // box volume -> mm^3
    gsStopwatch swTotal; // wall time of the whole solve (setup + time loop)
    // -------------------------------------------------------------------------
    // 1. Basis & Problem Setup (Uniform Tensor Product Basis, No Adaptivity)
    // -------------------------------------------------------------------------
    // Solution space: quadratic B-splines (the numerator of the NURBS geometry,
    // i.e. the same knots; the geometry itself stays the exact NURBS torus).
    // Target: C1 everywhere, as the C1 quadratic NURBS space of Lorenzo et al.
    // (2024).
    gsMultiBasis<> dbasis(mp, true);
    gsTensorBSplineBasis<3,T> * tb =
        dynamic_cast<gsTensorBSplineBasis<3,T>*>(&dbasis.basis(0));
    GISMO_ENSURE(tb, "Expected a tensor B-spline solution basis.");
    // (a) NurbsCircle has double knots at 0.25, 0.5, 0.75, so the inherited basis
    //     is only C0 at 90, 180 and 270 deg. Drop one copy, before refining while
    //     these are the only interior knots. The disk directions u, v have no
    //     interior knots, so refinement gives single knots (C1) there.
    tb->knots(2).reduceMultiplicity(1);
    for (index_t r = 0; r != numRefine; ++r)
        for (size_t p = 0; p != dbasis.nBases(); ++p)
        {
            dbasis.basis(p).uniformRefine();
        }
    // Toroidal elements are ~2.4x longer than the cross-section ones: refine w only
    for (index_t r = 0; r != refTor; ++r)
        for (size_t p = 0; p != dbasis.nBases(); ++p)
            dbasis.basis(p).uniformRefine(1, 1, 2);
    // (b) The seam w = 0 = 1. Periodic B-splines: after refinement the w-knots
    //     are uniform, h = 1/n. Replace the clamped knot vector by the UNCLAMPED
    //     uniform one  {-p h, ..., -h, 0, h, ..., 1, 1+h, ..., 1+p h}  (same
    //     domain [0,1], same n elements, n+p functions). Its first p functions
    //     (which stick out below w = 0) and its last p (which stick out above
    //     w = 1) are the two halves of the same periodic B-splines, so their
    //     DoFs are identified in the mapper below (periodicPairs). The result is
    //     the periodic quadratic space: C1 across the seam, n functions in w.
    //     (gismo's own setPeriodic() is not used: on this clamped knot vector it
    //      shrinks the domain and loses elements.)
    std::vector<std::pair<index_t,index_t> > periodicPairs;
    if (periodic)
    {
        const short_t pw = tb->degree(2);
        const index_t nEl = tb->knots(2).numElements();
        const T h = 1.0 / nEl;
        // uniform spacing, single interior knots, clamped ends (p+1 fold)
        const std::vector<T> uk = tb->knots(2).unique();
        for (size_t k = 1; k < uk.size(); ++k)
            GISMO_ENSURE(math::abs(uk[k] - uk[k-1] - h) < 1e-12,
                         "Periodic w-basis needs uniform knots.");
        GISMO_ENSURE(static_cast<index_t>(tb->knots(2).size()) == nEl + 2*pw + 1,
                     "Periodic w-basis needs single interior knots.");
        std::vector<T> kw;
        for (index_t i = -pw; i <= nEl + pw; ++i)
            kw.push_back(i * h);
        tb->knots(2) = gsKnotVector<T>(kw, pw);
        const index_t nu = tb->component(0).size(), nv = tb->component(1).size(),
                      nw = tb->component(2).size(); // = nEl + pw
        GISMO_ENSURE(nw == nEl + pw && tb->knots(2).numElements() == nEl,
                     "Unexpected periodic w-basis.");
        for (index_t j = 0; j != pw; ++j)
            for (index_t b = 0; b != nv; ++b)
                for (index_t a = 0; a != nu; ++a)
                    periodicPairs.push_back(std::make_pair(tb->index(a, b, j),
                                                           tb->index(a, b, nEl + j)));
        gsInfo << "Toroidal basis: periodic (C1 at the seam), " << nEl
               << " functions in w (" << nw << " before identifying the "
               << pw << " wrapped ones).\n";
    }
    else
        gsInfo << "Toroidal basis: glued seam interface (C0 at the seam).\n";
    gsInfo << "Elements (u x v x w): " << tb->component(0).numElements() << " x "
           << tb->component(1).numElements() << " x " << tb->component(2).numElements()
           << " = " << tb->numElements() << "\n";

    // Physical element size [mm]: per element h_e = its longest edge (chord
    // between the mapped corners). Also the longest cross-section (u,v) and
    // toroidal (w) edge separately. Same definition as h in the immersed code
    // (cube edge), so the two can be put on one axis.
    real_t hMax = 0, hMin = std::numeric_limits<real_t>::max(), hSum = 0;
    real_t hCrossMax = 0, hTorMax = 0;
    {
        typename gsBasis<T>::domainIter domIt    = tb->domain()->beginAll();
        typename gsBasis<T>::domainIter domItEnd = tb->domain()->endAll();
        gsMatrix<T> c(3, 8), x;
        index_t nE = 0;
        for (; domIt < domItEnd; ++domIt)
        {
            const gsVector<T> lo = domIt.lowerCorner(), up = domIt.upperCorner();
            for (index_t k = 0; k != 8; ++k)
                for (index_t d = 0; d != 3; ++d)
                    c(d, k) = ((k >> d) & 1) ? up[d] : lo[d];
            mp.patch(0).eval_into(c, x);
            real_t he = 0;
            for (index_t k = 0; k != 8; ++k)
                for (index_t d = 0; d != 3; ++d)
                    if (!((k >> d) & 1))
                    {
                        const real_t len = (x.col(k) - x.col(k | (1 << d))).norm() * parPhys.L;
                        he = math::max(he, len);
                        if (d == 2) hTorMax   = math::max(hTorMax, len);
                        else        hCrossMax = math::max(hCrossMax, len);
                    }
            hMax = math::max(hMax, he);
            hMin = math::min(hMin, he);
            hSum += he;
            ++nE;
        }
        hSum /= nE; // mean
    }
    gsInfo << "Element size (longest edge per element): h_max = " << hMax << " mm, h_min = "
           << hMin << " mm, h_mean = " << hSum << " mm; longest cross-section edge "
           << hCrossMax << " mm, longest toroidal edge " << hTorMax << " mm\n";

    // Computational mesh for ParaView: the element boundaries of the refined
    // solution basis, mapped through the exact geometry (curved edges).
    //   mesh_volume : all element edges (3D wireframe)
    //   mesh_surface: element edges on the torus surface (4 disk sides)
    if (plotMesh)
    {
        gsFileManager::mkdir(out);
        // Exact torus SURFACE (its 4 disk sides), finely sampled so it renders
        // smooth: 300 x 300 points per side. (Writing the volume patch instead
        // would also draw the seam disk and the degenerate disk corners, and
        // with few points it looks faceted.)
        {
            gsMultiPatch<T> torusSurf;
            const boxSide gsides[4] = {boundary::west, boundary::east,
                                       boundary::south, boundary::north};
            for (index_t s = 0; s != 4; ++s)
                torusSurf.addPatch(mp.patch(0).boundary(gsides[s]));
            gsWriteParaview(torusSurf, out + "/torus_geometry", 90000);
        }
        writeSingleCompMesh(*tb, mp.patch(0), out + "/mesh_volume", 2); // ~130 MB at -r 4 with 4
        gsParaviewCollection surf(out + "/mesh_surface");
        const boxSide sides[4] = {boundary::west, boundary::east,
                                  boundary::south, boundary::north};
        for (index_t s = 0; s != 4; ++s)
        {
            typename gsBasis<T>::uPtr bb = tb->boundaryBasis(sides[s]);
            typename gsGeometry<T>::uPtr bg = mp.patch(0).boundary(sides[s]);
            const std::string fn = out + "/mesh_surface_" + util::to_string(s);
            writeSingleCompMesh(*bb, *bg, fn, 6);
            surf.addPart(gsFileManager::getFilename(fn) + ".vtp");
        }
        surf.save();
        gsInfo << "Mesh written: " << out << "/mesh_volume.vtp, "
               << out << "/mesh_surface.pvd, " << out << "/torus_geometry.pvd\n";
    }

    real_t hmax    = dbasis.basis(0).getMaxCellLength();
    real_t lambda  = CHopt.askReal("lambda", 1.0 / (32.0 * pow(EIGEN_PI, 2)));
    real_t penalty = 1e4 * lambda;

    // -------------------------------------------------------------------------
    // 2. Time Integration Setup (Generalized-Alpha Method)
    // -------------------------------------------------------------------------
    real_t rho_inf = TIMEopt.askReal("rho_inf", 0.5);
    real_t alpha_m = 0.5 * (3.0 - rho_inf) / (1.0 + rho_inf);
    real_t alpha_f = 1.0 / (1.0 + rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;
    
    // Parameters (fixed model constants, set from the command line)
    real_t D   = par.D;   // tumor cell diffusivity
    real_t rho = par.rho; // net proliferation rate
    gsInfo << "Solver (box) units: D = " << D << " box^2/day, rho = " << rho
           << " 1/day, R = " << par.R << ", sigma = " << par.sigma << ", amp = " << par.amp
           << ", N_th = " << par.Nth << "\n";
    // unknowns = basis functions minus the wrapped copies merged at the periodic seam
    gsInfo << "Discretization: numRefine = " << numRefine << ", DoFs (unknowns) = "
           << dbasis.size() - (index_t)periodicPairs.size() << " (" << dbasis.size()
           << " basis functions, " << periodicPairs.size() << " merged at the seam)"
           << ", dt = " << dt << ", steps = " << maxSteps << "\n";

    // Initialization biomarkers
    real_t V_T, V_N, V_P, N_dash, N_T, A_p;
    real_t int_N, int_H; // \int_Omega N dOmega  and  \int_Omega H(N-N_th) dOmega

    // Output (identical columns in the immersed example):
    //  biomarkers.csv    one row per time step, t = 0 included. Volumes in mm^3,
    //                    time in days. V_T = int_H; growth = \int_Omega rho N(1-N)
    //                    [mm^3/day]; dintN_dt = (int_N^n - int_N^{n-1})/dt.
    //                    Mass balance (zero flux): dintN_dt ~ growth.
    //  computational.csv one row per time step: Newton iterations, timings.
    //  computational_summary.csv  quantity,value: DoFs, quadrature points,
    //                    IC projection error, setup/total times, means.
    std::ofstream csvFile, compFile;
    gsFileManager::mkdir(out);
    csvFile.open(out+"/biomarkers.csv");
    csvFile << std::setprecision(12); // V_P_relerr needs more than the default 6 digits
    // V_exact = 2 pi^2 R r^2 L^3 [mm^3], V_P_relerr = |V_P - V_exact| / V_exact
    csvFile << "TimeStep,Time,NumDOFs,int_N,int_H,V_P,N_dash,N_T,A_p,N_max,N_min,"
               "V_exact,V_P_relerr,V_T,V_N,growth,dintN_dt\n";
    compFile.open(out+"/computational.csv");
    compFile << std::setprecision(8)
             << "TimeStep,Time,NewtonIts,LinIts,ResidualNorm,StepSolve_s,StepTotal_s,Wall_s\n";
    
    real_t N_hat_th = par.Nth;  // Cell density threshold defining Omega_T (tumor)
    real_t k_steep = 200.0;  // Steepness factor for Heaviside step function


    gsConstantFunction<> alpha_f_func(alpha_f, dim);
    gsConstantFunction<> alpha_m_func(alpha_m, dim);

    // -------------------------------------------------------------------------
    // 3. Assembler & Expression Setup
    // -------------------------------------------------------------------------
    gsExprAssembler<> A(1, 1);
    A.options().setSwitch("SameElement", Aopt.askSwitch("SameElement", true));
    A.setIntegrationElements(dbasis);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(mp);
    space       w = A.getSpace(dbasis);

    gsMatrix<> Nnew, dNnew, Nold, dNold;
    gsMatrix<> Nalpha, dNalpha, dNupdate, Q;

    gsMultiPatch<> mp_nold, mp_dnold;

    auto nold  = A.getCoeff(mp_nold);
    auto dnold = A.getCoeff(mp_dnold);

    solution nnew_sol  = A.getSolution(w, Nnew);
    solution dnnew_sol = A.getSolution(w, dNnew);

    auto af = A.getCoeff(alpha_f_func);
    auto am = A.getCoeff(alpha_m_func);

    // Interpolated state variables
    auto N    = nold  + af * (nnew_sol - nold); // N is the tumor cell density
    auto dN    = dnold + am * (dnnew_sol - dnold);
    auto gradN = igrad(nold, G) + af * (igrad(nnew_sol, G) - igrad(nold, G));

    // Domain residual of  dN/dt = div(D grad N) + rho*N*(1-N):
    //   \int w*dN + \int D*grad(w).grad(N) - \int rho*w*N*(1-N) = 0
    auto residual = w * dN +
                    D * igrad(w, G) * gradN.tr() +
                    (-rho) * w * N.val() * (1.0 - N.val());

    w.setup(bc, dirichlet::l2Projection, 0);
    if (!periodicPairs.empty())
    {   // periodic seam: identify the wrapped DoFs (see (b) above)
        gsDofMapper pmap(dbasis);
        for (size_t k = 0; k != periodicPairs.size(); ++k)
            pmap.matchDof(0, periodicPairs[k].first, 0, periodicPairs[k].second);
        pmap.finalize();
        w.setupMapper(give(pmap));
    }
    A.initSystem();

    gsSparseMatrix<> K, K_const, K_nitsche;
    gsInfo << "Linear solver: " << linSolver
           << (linSolver == "direct" ? "" : " (CG tol " + util::to_string(cgTol) + ")") << "\n";

    // -------------------------------------------------------------------------
    // 4. Initial Condition Setup (smooth hyperbolic-tangent tumor seed)
    // -------------------------------------------------------------------------
    gsFunctionExpr<T> N0 = tumorSeed<dim,T>(par, mp);

    // L2-projection of the seed onto the solution space, with the seed
    // evaluated at PHYSICAL points: getCoeff(N0, G) composes N0 with the
    // geometry map, so the tanh ball is round in x,y,z and independent of the
    // parametrization. (gsL2Projection::project uses getCoeff(N0) WITHOUT G,
    // i.e. it evaluates N0 at the parametric point (u,v,w) -- on any non-
    // identity map that gives a seed deformed with the geometry.)
    // Projecting in the assembler's own space also respects the glued seam.
    gsMatrix<T> N0coefs;
    real_t icRelErr = 0, tIC = 0;
    {
        gsStopwatch swIC;
        auto f0 = A.getCoeff(N0, G);
        A.initSystem();
        A.assemble(w * w.tr() * meas(G), w * f0 * meas(G));
        // mass matrix: SPD and well conditioned, same linear solver as the steps
        fkSolveLinear<T>(linSolver, A.matrix(), A.rhs(), Nold, cgMaxIt, cgTol);
        tIC = swIC.stop();

        solution n0_sol = A.getSolution(w, Nold);
        n0_sol.extract(N0coefs, 0);

        // relative L2(Omega) error ||N_h(0) - N0|| / ||N0|| (unit-free, so the
        // same number in both codes)
        gsExprEvaluator<T> ev0(A);
        const real_t errIC = math::sqrt(ev0.integral((n0_sol - f0).sqNorm() * meas(G)));
        const real_t nrmIC = math::sqrt(ev0.integral(f0.sqNorm() * meas(G)));
        icRelErr = errIC / nrmIC;
        gsInfo << "IC L2-projection error (physical space): " << errIC
               << ", relative " << icRelErr << "\n";
        A.initSystem(); // clear matrix/rhs before the constant terms below
    }

    mp_nold.addPatch(dbasis.basis(0).makeGeometry(N0coefs));
    mp_dnold.addPatch(dbasis.basis(0).makeGeometry(gsMatrix<T>::Zero(dbasis.basis(0).size(),1)));

    dNold.setZero(A.numDofs(), 1);

    // -------------------------------------------------------------------------
    // 5. Pre-assemble Constant Terms Matrix
    // -------------------------------------------------------------------------
    A.assemble(meas(G) * (w * w.tr() * alpha_m +
               (alpha_f * gamma * dt) * (D * igrad(w, G) * igrad(w, G).tr())));
    K_const = A.giveMatrix();

    // -------------------------------------------------------------------------
    // 6. ParaView output of the cell density N (one time step per frame)
    // -------------------------------------------------------------------------
    // N and the geometry are sampled on the SAME tensor grid of the geometry's
    // parameter domain [0,1]^3 and written as one structured grid per frame.
    // (gsParaviewCollection with an expression evaluator samples the field on
    //  the bounding box of the SOLUTION basis instead. With the periodic
    //  (unclamped) w-knots that box is [-2h, 1+2h] in w, so the field array
    //  gets a different number of points than the geometry grid and ParaView
    //  rejects it: "The data array in the element may be too short".)
    const gsMatrix<T> parBox = mp.patch(0).support(); // [0,1]^3
    gsVector<unsigned> np(3);
    np.setConstant(static_cast<unsigned>(math::ceil(math::pow((T)plotPts, (T)1/3))));
    const gsMatrix<T> plotGrid = gsPointGrid<T>(parBox.col(0), parBox.col(1), np);
    const gsMatrix<T> plotXYZ  = mp.patch(0).eval(plotGrid); // static geometry
    gsParaviewCollection collection(out + "/N");
    if (plot) gsFileManager::mkdir(out + "/N_pvd");
    index_t frame = 0;
    auto writeFrame = [&](T time)
    {
        const gsMatrix<T> vals = mp_nold.patch(0).eval(plotGrid);
        const std::string rel = "N_pvd/N_" + util::to_string(frame++) + ".vts";
        std::ofstream f((out + "/" + rel).c_str());
        f << std::setprecision(8);
        f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n"
          << "<StructuredGrid WholeExtent=\"0 " << np[0]-1 << " 0 " << np[1]-1 << " 0 " << np[2]-1 << "\">\n"
          << "<Piece Extent=\"0 " << np[0]-1 << " 0 " << np[1]-1 << " 0 " << np[2]-1 << "\">\n"
          << "<PointData Scalars=\"N\">\n<DataArray type=\"Float64\" Name=\"N\" format=\"ascii\">\n";
        for (index_t j = 0; j != vals.cols(); ++j) f << vals(0,j) << " ";
        f << "\n</DataArray>\n</PointData>\n<Points>\n"
          << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (index_t j = 0; j != plotXYZ.cols(); ++j)
            f << plotXYZ(0,j) << " " << plotXYZ(1,j) << " " << plotXYZ(2,j) << " ";
        f << "\n</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
        collection.addPart(rel, time);
    };
    // Max/min of N, sampled on a parametric grid of the patch (diagnostic only:
    // tells whether the tumour core gets denser or is being diluted). Also on
    // the geometry's [0,1]^3, not on the (wider) periodic basis support.
    const gsMatrix<T> sampleGrid = gsPointGrid<T>(parBox, 40000);
    auto rangeN = [&](T & Nmax, T & Nmin)
    {
        gsMatrix<T> vals = mp_nold.patch(0).eval(sampleGrid);
        Nmax = vals.maxCoeff();
        Nmin = vals.minCoeff();
    };
    T Nmax, Nmin;
    rangeN(Nmax, Nmin);
    gsInfo << "Initial condition: max N = " << Nmax << ", min N = " << Nmin << "\n";

    const T icNmax = Nmax, icNmin = Nmin;

    if (plot)
        writeFrame(0); // initial condition; frame i of the collection = row i of the csv file

    // Biomarkers of the current solution Nnew (row n of biomarkers.csv)
    real_t growth = 0, int_N_prev = 0;
    auto writeBio = [&](index_t n, real_t time)
    {
        // Indicator expression for tumor domain Omega_T (where N >= N_th) | N_th is 0.15
        auto tumor_indicator = 1.0 / (1.0 + (-(k_steep) * (nnew_sol.val() - N_hat_th)).exp());

        // ==============================================================
        // Quantities of interest of the forward simulation
        int_N = computeInt(dbasis, mp, nnew_sol.val(), A.options());  // \int_Omega N dOmega
        int_H = computeInt(dbasis, mp, tumor_indicator, A.options()); // \int_Omega H(N-N_th) dOmega
        gsInfo << "  int_Omega N dOmega          = " << int_N << "\n"
               << "  int_Omega H(N-N_th) dOmega  = " << int_H << "\n";

        // Compute biomarkers
        V_T = int_H; // Tumor volume
        auto VN_field = nnew_sol.val() * tumor_indicator; // needed tumor_indicator because we integrate over Omega_T
        V_N = computeInt(dbasis, mp, VN_field, A.options()); // total tumor cell volume

        V_P = computeInt(dbasis, mp, 1.0, A.options()); // prostate volume

        N_dash = V_N / V_T; // mean normalized tumor cell density
        N_T = N_dash * V_T / V_P; // total tumor index
        auto A_p_field = rho * nnew_sol.val() * (1.0 - nnew_sol.val()) * tumor_indicator; //we integrate over Omega_T
        A_p = computeInt(dbasis,mp,A_p_field,A.options()); // mean proliferation activity of the tumor
        // Mass balance: with zero flux, d/dt \int_Omega N = \int_Omega rho N (1-N)
        auto growth_field = rho * nnew_sol.val() * (1.0 - nnew_sol.val()); // over ALL of Omega
        growth = computeInt(dbasis, mp, growth_field, A.options());
        // ==============================================================

        // Write in a csv file
        rangeN(Nmax, Nmin);
        gsInfo << "  V_T = " << int_H*L3 << " mm^3, max N = " << Nmax << ", min N = " << Nmin << "\n";

        csvFile << n << "," << time << "," << A.numDofs() << ","
                << int_N*L3 << "," << int_H*L3 << ","
                << V_P*L3 << "," << N_dash << "," << N_T << "," << A_p*L3 << ","
                << Nmax << "," << Nmin << ","
                << volExact << "," << math::abs(V_P*L3 - volExact) / volExact << ","
                << V_T*L3 << "," << V_N*L3 << "," << growth*L3 << ",";
        if (n > 0) csvFile << (int_N - int_N_prev)*L3 / dt;
        csvFile << "\n";
        csvFile.flush();
        int_N_prev = int_N;
    };
    Nnew = Nold;
    writeBio(0, 0.0);
    const real_t V_P_relerr = math::abs(V_P*L3 - volExact) / volExact;

    // Quadrature points: (Gauss rule per element) x (elements)
    const index_t nQuad = gsQuadrature::get(dbasis.basis(0), A.options()).numNodes()
                          * dbasis.basis(0).numElements();
    const real_t tSetup = swTotal.stop();

    // -------------------------------------------------------------------------
    // 7. Time Integration & Newton-Raphson Solving Loop
    // -------------------------------------------------------------------------
    index_t maxIt = 50;
    real_t tol    = TIMEopt.askReal("tol", 1e-4);
    index_t itsTotal = 0, linItsTotal = 0;
    std::string condCsv; // filled by --cond
    real_t  tSolveTotal = 0;
    gsStopwatch swLoop;

    for (index_t step = 0; step < maxSteps; ++step)
    {
        gsInfo << "Time Step " << step + 1 << " / " << maxSteps << "\n";
        gsStopwatch swStep;

        Nnew  = Nold;
        dNnew = (gamma - 1.0) / gamma * dNold;

        index_t its = maxIt, linIts = 0; // Newton its, linear-solver its (summed over Newton)
        real_t  resNorm = 0;
        for (index_t it = 0; it < maxIt; ++it)
        {
            Nalpha  = Nold  + alpha_f * (Nnew  - Nold);
            dNalpha = dNold + alpha_m * (dNnew - dNold);

            // Domain Residual Assembly
            A.initMatrix();
            A.clearRhs();
            A.assemble(residual * meas(G));
            Q = A.rhs();
            resNorm = Q.norm();

            if (it > 0 && resNorm < tol)
            {
                gsInfo << "  Converged in " << it << " iterations.\n";
                its = it;
                break;
            }

            // Tangent Matrix Assembly
            A.initMatrix();
            A.assemble(meas(G) * (-alpha_f * gamma * dt) *
                                 (rho * w * (1.0 - 2.0 * N.val()) * w.tr()));

            K = A.giveMatrix();
            K += K_const;

            // Condition number of the first system matrix (--cond)
            if (cond && step == 0 && it == 0)
                condCsv = fkReportCond(K);

            // Solve and Update
            const index_t li = fkSolveLinear<T>(linSolver, K, -Q, dNupdate, cgMaxIt, cgTol);
            linIts      += li;
            linItsTotal += li;

            dNnew += dNupdate;
            Nnew.noalias() += (gamma * dt) * dNupdate;
        }
        if (its == maxIt)
            gsWarn << "Newton did not converge in " << maxIt << " iterations (|Q| = "
                   << resNorm << ").\n";
        const real_t tSolve = swStep.stop(); // Newton loop only

        // Advance step
        nnew_sol.extract(mp_nold);
        dnnew_sol.extract(mp_dnold);
        Nold  = Nnew;
        dNold = dNnew;

        if (plot)
            writeFrame((step+1)*dt);

        writeBio(step+1, (step+1)*dt);

        itsTotal    += its;
        tSolveTotal += tSolve;
        compFile << step+1 << "," << (step+1)*dt << "," << its << "," << linIts << "," << resNorm << ","
                 << tSolve << "," << swStep.stop() << "," << swTotal.stop() << "\n";
        compFile.flush();
    }
    csvFile.close();
    compFile.close();
    if (plot) collection.save();

    const real_t tLoop = swLoop.stop(), tTotal = swTotal.stop();
    std::ofstream sumFile((out + "/computational_summary.csv").c_str());
    sumFile << std::setprecision(10) << "quantity,value\n"
            << "method,bodyfitted\n"
            << "numRefine," << numRefine << "\nrefTor," << refTor << "\n"
            << "degree," << tb->degree(0) << "\n"
            << "seam," << (periodic ? "periodic_C1" : "glued_C0") << "\n"
            << "elements," << dbasis.basis(0).numElements() << "\n"
            << "h_max_mm," << hMax << "\nh_min_mm," << hMin << "\nh_mean_mm," << hSum << "\n"
            << "h_cross_max_mm," << hCrossMax << "\nh_tor_max_mm," << hTorMax << "\n"
            << "DoFs," << A.numDofs() << "\n"
            << "DoFs_active," << A.numDofs() << "\n"
            << "quad_points_interior," << nQuad << "\nquad_points_cut,0\n"
            << "quad_points_total," << nQuad << "\n"
            << "time_quad_build_s,0\n"
            << "time_IC_projection_s," << tIC << "\n"
            << "IC_L2_relerr," << icRelErr << "\n"
            << "IC_Nmax," << icNmax << "\nIC_Nmin," << icNmin << "\n"
            << "V_P_relerr," << V_P_relerr << "\n"
            << "dt," << dt << "\nsteps," << maxSteps << "\n"
            << "newton_its_total," << itsTotal << "\n"
            << "newton_its_mean," << (maxSteps ? (real_t)itsTotal / maxSteps : 0) << "\n"
            << "lin_solver," << linSolver << "\ncg_tol," << cgTol << "\n"
            << "lin_its_total," << linItsTotal << "\n"
            << "lin_its_per_newton_mean," << (itsTotal ? (real_t)linItsTotal / itsTotal : 0) << "\n"
            << condCsv
            << "time_setup_s," << tSetup << "\n"
            << "time_solve_total_s," << tSolveTotal << "\n"
            << "time_solve_per_step_mean_s," << (maxSteps ? tSolveTotal / maxSteps : 0) << "\n"
            << "time_loop_s," << tLoop << "\n"
            << "time_per_step_mean_s," << (maxSteps ? tLoop / maxSteps : 0) << "\n"
            << "time_total_s," << tTotal << "\n";
    gsInfo << "Times: setup " << tSetup << " s, time loop " << tLoop << " s ("
           << (maxSteps ? tLoop / maxSteps : 0) << " s/step, Newton solve "
           << (maxSteps ? tSolveTotal / maxSteps : 0) << " s/step, "
           << (maxSteps ? (real_t)itsTotal / maxSteps : 0) << " its/step), total "
           << tTotal << " s\n";

}

// =============================================================================
// MAIN FUNCTION
// =============================================================================
int main(int argc, char *argv[])
{
    std::string fn("pde/fisher_kolmogorov_cube_bvp.xml");
    std::string out("output");

    // Pre-scan the command line for the setup file, so that the option lists it
    // contains can serve as the defaults of all the other arguments below.
    // Precedence: built-in default  <  xml file  <  command line.
    for (int i = 1; i+1 != argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-f" || arg == "--file") fn = argv[i+1];
    }

    // Read input problem settings (the geometry of the xml file is NOT used)
    gsFileData<> fd(fn);
    gsInfo << "Loaded configuration file: " << fd.lastPath() << "\n";

    gsOptionList CHopt, TIMEopt, Aopt;
    fd.getId(3, CHopt);   // Physics options
    fd.getId(4, TIMEopt); // Integrator options
    fd.getId(5, Aopt);    // Assembler options

    // Model parameters: PHYSICAL units (mm, day), taken from the command line
    // only. The labels of xml id 3 are deliberately NOT read: in the cube xml
    // they are in box units (e.g. radius = 0.15), and reading them would
    // silently override the mm/day defaults below.
    fkParams<real_t> par;

    // Time integration (xml id 4: rho_inf, tol) and discretization.
    // Default horizon T = dt * Nsteps = 730 days (2 years between MRI scans).
    real_t  dt        = 5.0;
    index_t maxSteps  = 146;
    index_t numRefine = 3;
    bool    plot      = false;

    // Torus geometry
    torusParams tp;
    tp.R  = Aopt.askReal("torusR",  tp.R );
    tp.r  = Aopt.askReal("torusr",  tp.r );
    tp.cx = Aopt.askReal("torusCx", tp.cx);
    tp.cy = Aopt.askReal("torusCy", tp.cy);
    tp.cz = Aopt.askReal("torusCz", tp.cz);

    // Command-line interface
    gsCmdLine cmd("Fisher-Kolmogorov Generalized-Alpha Solver on a body-fitted solid torus (forward problem)");
    cmd.addReal("t", "dt", "Time step size [day]", dt);
    cmd.addInt("N", "Nsteps", "Number of time steps", maxSteps);
    cmd.addInt("r", "numRefine", "Number of uniform refinement steps", numRefine);
    index_t refTor = 0, plotPts = 50000;
    cmd.addInt("", "refTor", "Extra refinement steps in the toroidal direction only", refTor);
    cmd.addInt("", "plotPts", "ParaView sample points per frame (~plotPts^(1/3) per direction)", plotPts);
    cmd.addString("f", "file", "Input XML setup file", fn);
    cmd.addString( "o", "output", "Output directory", out);
    cmd.addSwitch("plot", "Write ParaView files of N for every time step", plot);
    bool C0seam = false, plotMesh = false;
    cmd.addSwitch("C0seam", "Glue the toroidal seam as a C0 interface instead of a periodic (C1) basis", C0seam);
    cmd.addSwitch("plotMesh", "Write the computational mesh (element edges) to ParaView", plotMesh);
    std::string linSolver = "cgdiag";
    index_t cgMaxIt = 20000;
    real_t  cgTol   = 1e-8;
    bool    cond    = false;
    cmd.addString("", "solver", "Linear solver: cgdiag (default) | cgic | direct", linSolver);
    cmd.addInt   ("", "cgMaxIt", "Max CG iterations per linear solve", cgMaxIt);
    cmd.addReal  ("", "cgTol", "CG tolerance, relative to ||rhs||", cgTol);
    cmd.addSwitch("cond", "Estimate the condition number of the first system matrix (K and Jacobi-scaled)", cond);
    cmd.addReal("D", "diffusivity", "Tumor cell diffusivity D [mm^2/day]", par.D);
    cmd.addReal("p", "rho", "Net proliferation rate rho [1/day]", par.rho);
    cmd.addReal("R", "radius", "Radius R of the initial tumor [mm]", par.R);
    cmd.addReal("s", "sigma", "Interface width sigma of the initial tumor [mm]", par.sigma);
    cmd.addReal("a", "amp", "Peak density of the initial tumor (0..1]", par.amp);
    cmd.addReal("L", "length", "Physical length of one box unit [mm]", par.L);
    cmd.addReal("", "Nth", "Cell density threshold N_th defining Omega_T", par.Nth);
    cmd.addReal("", "xc", "x-coordinate of the tumor centre [mm] (default: on the torus core circle)", par.xc);
    cmd.addReal("", "yc", "y-coordinate of the tumor centre [mm] (default: on the torus core circle)", par.yc);
    cmd.addReal("", "zc", "z-coordinate of the tumor centre [mm] (default: on the torus core circle)", par.zc);
    cmd.addReal("", "torusR", "Major radius of the torus [box units]", tp.R);
    cmd.addReal("", "torusr", "Minor (tube) radius of the torus [box units]", tp.r);
    cmd.addReal("", "torusCx", "x-coordinate of the torus centre [box units]", tp.cx);
    cmd.addReal("", "torusCy", "y-coordinate of the torus centre [box units]", tp.cy);
    cmd.addReal("", "torusCz", "z-coordinate of the torus centre [box units]", tp.cz);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // -------------------------------------------------------------------------
    // Geometry: exact single-patch NURBS solid torus
    // -------------------------------------------------------------------------
    gsMultiPatch<> mp = makeTorus<real_t>(tp, !C0seam);
    const real_t L = par.L;
    const real_t volExact = 2.0 * EIGEN_PI * EIGEN_PI * tp.R * tp.r * tp.r * L * L * L;
    gsInfo << "Torus (box units): centre (" << tp.cx << "," << tp.cy << "," << tp.cz
           << "), R = " << tp.R << ", r = " << tp.r << "\n"
           << "Torus (1 box unit = " << L << " mm): R = " << tp.R * L << " mm, r = "
           << tp.r * L << " mm, outer diameter " << 2 * (tp.R + tp.r) * L
           << " mm, exact volume 2 pi^2 R r^2 = " << volExact << " mm^3\n";
    gsMatrix<> bb;
    mp.boundingBox(bb);
    gsInfo << "Torus bounding box: [" << bb(0,0) << "," << bb(0,1) << "] x ["
           << bb(1,0) << "," << bb(1,1) << "] x [" << bb(2,0) << "," << bb(2,1) << "]\n";
    if (bb.minCoeff() < 0 || bb.maxCoeff() > 1)
        gsWarn << "The torus does not fit inside [0,1]^3.\n";

    // The bounding-box centre is in the hole: default tumour centre on the core circle
    if (math::isnan(par.xc) && math::isnan(par.yc) && math::isnan(par.zc))
    {
        par.xc = (tp.cx + tp.R) * L;
        par.yc = tp.cy * L;
        par.zc = tp.cz * L;
    }
    if (par.R + 2 * par.sigma > tp.r * L)
        gsWarn << "Tumour seed R + 2 sigma = " << par.R + 2 * par.sigma
               << " mm exceeds the tube radius r = " << tp.r * L
               << " mm: the seed is cut by the torus surface.\n";
    GISMO_ENSURE(par.amp > 0 && par.amp <= 1, "The seed amplitude must lie in (0,1].");

    // -------------------------------------------------------------------------
    // Regime diagnostics, all in physical units. They tell in advance whether
    // the tumour will GROW (and get denser) or be DILUTED by diffusion over
    // the simulated horizon T. For Fisher-Kolmogorov:
    //   rho*T          : logistic growth; N grows ~exp(rho*T) while N << 1.
    //                    rho*T < 1  => hardly any growth in the horizon.
    //   D/R^2 vs rho   : diffusive loss of a seed of radius R vs proliferation.
    //                    D/R^2 > rho => the seed spreads out faster than it
    //                    grows: max N and V_T DECREASE (diffusion-dominated).
    //   c = 2 sqrt(D rho)   : speed of the invasion front once established.
    //   l = sqrt(D/rho)     : front width scale; the tanh front spans ~4 l
    //                    and has to be resolved by the mesh (h <~ 2 l).
    // -------------------------------------------------------------------------
    {
        const real_t T     = dt * maxSteps;
        const real_t c     = 2 * math::sqrt(par.D * par.rho);
        const real_t ell   = math::sqrt(par.D / par.rho);
        const real_t nEl   = math::pow(2.0, (real_t)numRefine);
        const real_t hCross= 2 * tp.r * L / nEl;                          // across the tube
        const real_t hTor  = 2 * EIGEN_PI * (tp.R + tp.r) * L / (4 * nEl * math::pow(2.0, (real_t)refTor)); // along the outer rim
        gsInfo << "Regime over T = " << T << " days:\n"
               << "  rho*T          = " << par.rho * T << "  (growth factor exp(rho T) = "
               << math::exp(par.rho * T) << ")\n"
               << "  D/R^2          = " << par.D / (par.R * par.R) << " 1/day  vs  rho = "
               << par.rho << " 1/day\n"
               << "  front speed c  = " << c << " mm/day  (c*T = " << c * T << " mm)\n"
               << "  diffusion len. = sqrt(D T) = " << math::sqrt(par.D * T) << " mm\n"
               << "  front width    ~ 4 sqrt(D/rho) = " << 4 * ell << " mm,  mesh h ~ "
               << hCross << " mm (cross-section), " << hTor << " mm (toroidal)\n";
        if (par.rho * T < 1)
            gsWarn << "rho*T < 1: almost no proliferation within the horizon.\n";
        if (par.D / (par.R * par.R) > par.rho)
            gsWarn << "D/R^2 > rho: diffusion-dominated, the tumour will be diluted "
                      "(max N and V_T decrease) instead of growing.\n";
        if (math::max(hCross, hTor) > 2 * 4 * ell)
            gsWarn << "The mesh is much coarser than the front width: refine (-r).\n";
    }

    // No boundary conditions: natural zero-flux on the whole torus surface
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsFileManager::mkdir(out);
    // --plot also writes the geometry and the computational mesh (they are
    // static, written once); --plotMesh writes them without the N frames.
    solve<3, real_t>(mp, bc, CHopt, TIMEopt, Aopt, dt, maxSteps, numRefine, par, out, plot, refTor, plotPts, volExact,
                     !C0seam, plot || plotMesh, linSolver, cgMaxIt, cgTol, cond);

    gsInfo << "Exact torus volume (compare with V_P in the csv): " << volExact << " mm^3\n";
    return EXIT_SUCCESS;
}
