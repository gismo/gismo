/** @file fischerkolmogrov_torus_example_immersed.cpp

    @brief Fisher-Kolmogorov tumour growth on an IMMERSED solid torus, the
           counterpart of fischerkolmogrov_torus_example_bodyfitted.cpp.

    Same model, same inputs, same outputs as the body-fitted torus example
    (generalized-alpha, Newton, tanh seed, biomarkers, csv columns, physical
    units mm/day). The ONLY difference is how the domain is represented and
    integrated, exactly as in fischerkolmogrov_prostate_example_immersed.cpp:

      - the domain is Omega = {phi < 0} inside the background box [0,1]^3,
        with phi the EXACT signed distance of the torus (no triangle mesh):
            phi(x) = sqrt( (sqrt((x-cx)^2+(y-cy)^2) - R)^2 + (z-cz)^2 ) - r,
        so the immersed domain is the very same torus as the exact NURBS
        patch of the body-fitted example (same R, r, centre, defaults);
      - the background box carries a plain tensor B-spline basis (degree
        1 + e, default quadratic like the body-fitted one);
      - interior cells: tensor Gauss rule; cut cells: Algoim octree rule;
        outside cells: dropped. Zero flux on {phi = 0} is natural (no terms).

    PARAMETRIC vs PHYSICAL SPACE
      Body-fitted: a NURBS map G : [0,1]^3 (u,v,w) -> torus. Everything is
        pulled back: grad = J^{-T} grad_u, dx = |det J| du, and functions of
        x (the seed) must be evaluated at G(u) (getCoeff(f, G)).
      Immersed: the background map is the IDENTITY of [0,1]^3, so parameter
        space == physical space (in box units) -- no pull-back, meas = 1,
        basis derivatives are already physical, the seed is evaluated
        directly at x. The geometry enters ONLY through phi, i.e. through
        which quadrature points exist.
      Both: 1 box unit = L mm (-L); all inputs are mm/day and converted once
        in fkToBox(); csv volumes are mm^3; ParaView coordinates are box
        units in BOTH examples, so their outputs overlay directly.

    COMPARING THE TWO
      Run both with the same physics flags; compare output_data.csv
      (V_P vs the exact 2 pi^2 R r^2 L^3, V_T(t), int_N(t), N_max(t)) and
      overlay the N fields in ParaView. Use matched element sizes (L=50):
        body-fitted -r k      : h ~ 15/2^k mm across the tube,
                    --refTor t : h ~ 141/(4 2^(k+t)) mm along the outer rim
        immersed    -r m      : h = 50/2^m mm (isotropic)
      e.g.  body-fitted -r 3 --refTor 1  (1.9 / 2.2 mm)  vs  immersed -r 5 (1.6 mm)
            body-fitted -r 4 --refTor 1  (0.9 / 1.1 mm)  vs  immersed -r 6 (0.8 mm)
      and check that both converge to the same curves under refinement.

    ./build/bin/fischerkolmogrov_torus_example_immersed \
        -D 1e-3 -p 5e-3 -t 5 -N 146 -R 3 -s 0.75 -a 0.5 -r 5 -L 50 \
        -o 00_output_FK_torus_immersed --plot

    Author(s): L. Venta-Viñuela
*/
#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>

#include <iomanip>
#include <map>
#include <vector>

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
// FISHER-KOLMOGOROV MODEL PARAMETERS -- PHYSICAL UNITS (identical to the
// body-fitted torus example): lengths in mm, time in days.
// =============================================================================
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
    // Centre of the initial tumor [mm]. NaN: on the core circle at 45 deg.
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

// =============================================================================
// GEOMETRY: exact signed distance of the solid torus (box units)
// =============================================================================
struct torusParams
{
    real_t R  = 0.30; // major radius (centre of the tube to the torus axis)
    real_t r  = 0.15; // minor radius (tube radius)
    real_t cx = 0.5, cy = 0.5, cz = 0.5; // centre; the axis is parallel to z
};

/// phi(x) = sqrt( (sqrt((x-cx)^2+(y-cy)^2) - R)^2 + (z-cz)^2 ) - r
/// Negative inside the solid torus, |grad phi| = 1. Its gradient is undefined
/// on the core circle (deep inside) and on the axis (in the hole); there it is
/// set to 0 -- neither set is near {phi = 0}, where Algoim needs it.
template <class T>
class gsTorusLevelSet : public gsFunction<T>
{
public:
    GISMO_CLONE_FUNCTION(gsTorusLevelSet)

    explicit gsTorusLevelSet(const torusParams & tp) : m_tp(tp)
    {
        m_bbox.resize(3, 2);
        m_bbox.col(0).setZero();
        m_bbox.col(1).setOnes();
    }

    short_t     domainDim() const override { return 3; }
    short_t     targetDim() const override { return 1; }
    gsMatrix<T> support()   const override { return m_bbox; }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t i = 0; i != u.cols(); ++i)
        {
            const T dx = u(0,i) - m_tp.cx, dy = u(1,i) - m_tp.cy, dz = u(2,i) - m_tp.cz;
            const T q  = math::sqrt(dx*dx + dy*dy) - m_tp.R;
            result(0,i) = math::sqrt(q*q + dz*dz) - m_tp.r;
        }
    }

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        const T eps = 1e-14;
        result.resize(3, u.cols());
        for (index_t i = 0; i != u.cols(); ++i)
        {
            const T dx = u(0,i) - m_tp.cx, dy = u(1,i) - m_tp.cy, dz = u(2,i) - m_tp.cz;
            const T rh = math::sqrt(dx*dx + dy*dy);
            const T q  = rh - m_tp.R;
            const T d  = math::sqrt(q*q + dz*dz);
            if (d < eps) { result.col(i).setZero(); continue; }
            const T fr = (rh < eps) ? (T)0 : q / (d * rh);
            result(0,i) = fr * dx;
            result(1,i) = fr * dy;
            result(2,i) = dz / d;
        }
    }

    std::ostream & print(std::ostream & os) const override
    { os << "Signed distance of a solid torus (R=" << m_tp.R << ", r=" << m_tp.r << ")"; return os; }

private:
    torusParams m_tp;
    gsMatrix<T> m_bbox;
};

// =============================================================================
// FINITE-CELL (IMMERSED) DISCRETIZATION SETTINGS
// =============================================================================
struct fcmParams
{
    index_t     degElevate   = 1;     // background basis degree = 1 + this
    index_t     quadDepth    = 2;     // octree depth of the cut-cell quadrature
    std::string indicator    = "integralChange"; // uniform|fallback|integralChange
    real_t      indicatorTol = 1e-2;  // integralChange acceptance tolerance
    index_t     plotRes      = 64;    // sample points per direction of the N field
    bool        plotPoints   = false; // also dump N at the quadrature points
    bool        plotMesh     = true;  // also write the background grid
    // Output coordinates: box units, like the body-fitted torus example, so
    // that the two results overlay in ParaView (x_out = toPhys*x_box + shift).
    real_t      toPhys       = 1.0;
    real_t      shift[3]     = {0, 0, 0};
};

// Map 3 x n box points (in place) to the output coordinates.
template <class T>
void fkToPhys(const fcmParams & fcm, gsMatrix<T> & pts)
{
    for (index_t d = 0; d != 3; ++d)
        pts.row(d) = (pts.row(d).array() * fcm.toPhys + fcm.shift[d]).matrix();
}

// Smooth spherical tumor seed (BOX units, i.e. the output of fkToBox):
//   N(x,0) = amp*0.5*(1 - tanh((||x - x_tumor|| - R)/sigma))
// The background map is the identity, so x here is the physical point.
template <short_t dim, class T>
gsFunctionExpr<T> tumorSeed(const fkParams<T> & par)
{
    std::ostringstream expr;
    expr.precision(16);
    expr << par.amp << "*0.5*(1-tanh((sqrt((x-(" << par.xc << "))^2+(y-(" << par.yc
         << "))^2+(z-(" << par.zc << "))^2)-(" << par.R << "))/(" << par.sigma << ")))";
    gsInfo << "Initial condition: N(x,0) = " << expr.str() << "\n";
    return gsFunctionExpr<T>(expr.str(), dim);
}

// =============================================================================
// IMMERSED QUADRATURE BLOCK
// =============================================================================
// This replaces gsExprAssembler's integration elements. One block holds the
// quadrature of one background cell that meets Omega = {phi < 0}, together with
// everything the assembly needs at those points. The background grid is static
// during the whole simulation, so the blocks -- basis values and derivatives
// included -- are built once and reused by every Newton iteration of every time
// step; only the coefficient vectors change.
template <class T>
struct fkQuadBlock
{
    gsMatrix<T>          pts;   // dim x nq       nodes (parametric == physical)
    gsVector<T>          wts;   // nq             weights
    gsMatrix<index_t>    act;   // na x nq        global index of active function i at point q
    gsMatrix<index_t>    loc;   // na x nq        slot of act(i,q) inside \a glob
    std::vector<index_t> glob;  // nu             local -> global map of the block
    gsMatrix<T>          bv;    // na x nq        basis values
    gsMatrix<T>          bd;    // (na*dim) x nq  basis derivatives
};

/// Evaluate the basis at (\a pts, \a wts) and cache it in \a b, together with a
/// contiguous local numbering of the functions active anywhere in the block.
/// The local numbering is needed because a quadrature point sitting exactly on
/// a cell face can be attributed to the neighbouring knot span, so the active
/// set is not guaranteed to be identical for every point of the block.
template <class T>
void fkBuildBlock(const gsBasis<T> & basis,
                  const gsMatrix<T> & pts, const gsVector<T> & wts,
                  fkQuadBlock<T> & b)
{
    b.pts = pts;
    b.wts = wts;
    basis.active_into(pts, b.act);
    basis.eval_into  (pts, b.bv);
    basis.deriv_into (pts, b.bd);

    b.glob.clear();
    b.loc.resize(b.act.rows(), b.act.cols());
    std::map<index_t,index_t> slot;
    for (index_t q = 0; q != b.act.cols(); ++q)
        for (index_t i = 0; i != b.act.rows(); ++i)
        {
            const std::pair<std::map<index_t,index_t>::iterator,bool> ins =
                slot.insert(std::make_pair(b.act(i,q), (index_t)b.glob.size()));
            if (ins.second) b.glob.push_back(b.act(i,q));
            b.loc(i,q) = ins.first->second;
        }
}

// =============================================================================
// VOLUMETRIC FIELD OUTPUT
// =============================================================================
// Writes one VTK XML structured grid (.vts) holding BOTH the cell density N and
// the level set phi, sampled on a uniform grid of the background box. Carrying
// phi in the same file is what makes the output usable for an immersed result:
// a single ParaView "Clip -> Scalars: phi -> Value 0 -> Invert" shows N on the
// prostate only, instead of the whole background box. (gsWriteParaview has no
// overload that puts two fields in one file, hence the hand-written .vts.)
//
// N is MASKED to 0 outside {phi<0}. Outside the immersed domain the spline is
// not a solution of anything: those coefficients are only held in place by the
// tiny background regularization, so they drift to large arbitrary values
// (O(10) after a few steps) which would otherwise dominate the ParaView colour
// range and hide the tumour. Masking costs nothing -- the equation is only
// solved, and the biomarkers only integrated, on {phi<0}.
template <class T>
void fkWriteField(const std::string & fn,
                  const gsVector<index_t> & res,  // sample points per direction
                  const gsMatrix<T> & gridPts,    // 3 x np, x fastest then y then z
                  const gsMatrix<T> & Nvals,      // 1 x np
                  const gsMatrix<T> & phiVals)    // 1 x np
{
    const index_t np = gridPts.cols();
    std::ofstream f(fn.c_str());
    f << std::setprecision(7);
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <StructuredGrid WholeExtent=\"0 " << res[0]-1 << " 0 " << res[1]-1
      << " 0 " << res[2]-1 << "\">\n"
      << "    <Piece Extent=\"0 " << res[0]-1 << " 0 " << res[1]-1
      << " 0 " << res[2]-1 << "\">\n"
      << "      <PointData Scalars=\"N\">\n"
      << "        <DataArray type=\"Float32\" Name=\"N\" format=\"ascii\">\n";
    for (index_t i = 0; i != np; ++i)
        f << Nvals(0,i) << (0 == (i+1) % 12 ? "\n" : " ");
    f << "\n        </DataArray>\n"
      << "        <DataArray type=\"Float32\" Name=\"phi\" format=\"ascii\">\n";
    for (index_t i = 0; i != np; ++i)
        f << phiVals(0,i) << (0 == (i+1) % 12 ? "\n" : " ");
    f << "\n        </DataArray>\n"
      << "      </PointData>\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (index_t i = 0; i != np; ++i)
        f << gridPts(0,i) << " " << gridPts(1,i) << " " << gridPts(2,i) << "\n";
    f << "        </DataArray>\n"
      << "      </Points>\n"
      << "    </Piece>\n"
      << "  </StructuredGrid>\n"
      << "</VTKFile>\n";
    f.close();
}

// =============================================================================
// SOLVER ROUTINE
// =============================================================================
template <short_t dim, class T>
void solve(gsMultiPatch<T> & mp,
           const gsFunction<T> & phi,
           gsBoundaryConditions<T> & bc,
           gsOptionList & CHopt,
           gsOptionList & TIMEopt,
           gsOptionList & Aopt,
           real_t & dt,
           index_t & maxSteps,
           index_t & numRefine,
           const fkParams<T> & parPhys,
           const fcmParams & fcm,
           std::string out,
           bool plot,
           real_t volExact = 0,  // exact |Omega| in mm^3 (0: unknown)
           const std::string & linSolver = "cgdiag", // cgdiag | cgic | direct
           index_t cgMaxIt = 20000,
           real_t  cgTol   = 1e-8,   // relative to ||rhs||
           bool    cond    = false)  // estimate the condition number of the first matrix
{
    // Everything below runs in BOX units (lengths / L); the csv converts back.
    const fkParams<T> par = fkToBox(parPhys);
    const real_t L3 = parPhys.L * parPhys.L * parPhys.L; // box volume -> mm^3
    gsStopwatch swTotal; // wall time of the whole solve (setup + time loop)
    // The immersed boundary {phi==0} carries the natural zero-flux condition,
    // which contributes no term to the weak form: the xml boundary conditions
    // (which describe the boundary of the background box) are not used.
    GISMO_UNUSED(bc);

    // -------------------------------------------------------------------------
    // 1. Basis & Problem Setup (Uniform Tensor Product Basis, No Adaptivity)
    //    The patch is now the BACKGROUND box [0,1]^dim, not a body-fitted patch.
    // -------------------------------------------------------------------------
    gsMultiBasis<T> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + fcm.degElevate);
    for (index_t r = 0; r != numRefine; ++r)
        for (size_t p = 0; p != dbasis.nBases(); ++p)
        {
            dbasis.basis(p).uniformRefine();
        }

    gsTensorBSplineBasis<dim,T> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<dim,T> *>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline background basis.");
    const index_t deg   = dbasis.maxCwiseDegree();
    const index_t nDofs = tbsPtr->size();

    real_t hmax    = dbasis.basis(0).getMaxCellLength();
    real_t lambda  = CHopt.askReal("lambda", 1.0 / (32.0 * pow(EIGEN_PI, 2)));
    real_t penalty = 1e4 * lambda;
    GISMO_UNUSED(hmax); GISMO_UNUSED(penalty);

    // -------------------------------------------------------------------------
    // 1b. IMMERSED QUADRATURE (this is the only part that differs from the
    //     body-fitted example). Classify the background cells and collect the
    //     quadrature of Omega = {phi < 0} once: interior cells get a tensor
    //     Gauss rule, cut cells the adaptive Algoim octree rule.
    // -------------------------------------------------------------------------
    const index_t nFallback = 4; // uniform samples/direction when Algoim returns 0 nodes

    gsImplicitTrimmedDomain<dim,T> tr_domain(phi, *tbsPtr);

    gsGaussRule<T> gauss(gsVector<index_t>::Constant(dim, deg + 1));
    gsOptionList qopt = gsAlgoimAdaptiveRule<T>::defaultOptions();
    qopt.setInt   ("maxDepth",     fcm.quadDepth);
    qopt.setInt   ("nFallback",    nFallback);
    qopt.setString("indicator",    fcm.indicator);
    qopt.setReal  ("indicatorTol", fcm.indicatorTol);
    gsAlgoimAdaptiveRule<T> volumeRule(phi, *tbsPtr, qopt);

    std::vector<fkQuadBlock<T> > blocks; // quadrature of Omega, one entry per cell
    gsMatrix<T> qpts;
    gsVector<T> qwts;
    index_t nQuadInterior = 0, nQuadCut = 0, nInteriorCells = 0, nCutCells = 0;

    gsStopwatch swQuad;
    for (auto it = tr_domain.beginInterior(); it != tr_domain.template end<InteriorSign>(); ++it)
    {
        gauss.mapTo(it.lowerCorner(), it.upperCorner(), qpts, qwts);
        ++nInteriorCells;
        if (qpts.cols() == 0) continue;
        blocks.push_back(fkQuadBlock<T>());
        fkBuildBlock(*tbsPtr, qpts, qwts, blocks.back());
        nQuadInterior += qpts.cols();
    }
    for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
    {
        volumeRule.mapTo(it.lowerCorner(), it.upperCorner(), qpts, qwts);
        ++nCutCells;
        if (qpts.cols() == 0) continue;
        blocks.push_back(fkQuadBlock<T>());
        fkBuildBlock(*tbsPtr, qpts, qwts, blocks.back());
        nQuadCut += qpts.cols();
    }
    const double tQuad = swQuad.stop();

    // Measure of the immersed domain: the plain sum of the weights. This is the
    // volume the biomarkers below are normalized with, so it is worth printing.
    T volOmega = 0;
    for (size_t e = 0; e != blocks.size(); ++e)
        volOmega += blocks[e].wts.sum();

    // The blocks cache the basis values/derivatives at every quadrature point,
    // which is what makes the repeated Newton assemblies cheap -- but it is also
    // what the memory of a run is spent on, so report it.
    size_t cacheBytes = 0;
    for (size_t e = 0; e != blocks.size(); ++e)
        cacheBytes += sizeof(T)*(blocks[e].bv.size() + blocks[e].bd.size()
                                 + blocks[e].pts.size() + blocks[e].wts.size())
                    + sizeof(index_t)*(blocks[e].act.size() + blocks[e].loc.size());

    // DoFs that actually carry the solution: those whose support meets Omega.
    // The others are frozen by the regularization; nDofs counts the whole box.
    index_t nActive = 0;
    std::vector<index_t> activeIdx(nDofs, -1); // global dof -> active index (-1: inactive)
    {
        std::vector<bool> touched(nDofs, false);
        for (size_t e = 0; e != blocks.size(); ++e)
            for (size_t i = 0; i != blocks[e].glob.size(); ++i)
                touched[blocks[e].glob[i]] = true;
        for (index_t i = 0; i != nDofs; ++i)
            if (touched[i]) activeIdx[i] = nActive++;
    }

    gsInfo << "Immersed quadrature: " << nInteriorCells << " interior cells ("
           << nQuadInterior << " pts) + " << nCutCells << " cut cells ("
           << nQuadCut << " pts), total " << (nQuadInterior+nQuadCut)
           << " points, built in " << tQuad << " s ("
           << (cacheBytes >> 20) << " MB cached).\n"
           << "  quadDepth = " << fcm.quadDepth << ", indicator = " << fcm.indicator
           << ", |Omega| = " << volOmega * L3 << " mm^3\n"
           << "Active DoFs (support meets Omega): " << nActive << " of " << nDofs << "\n";
    if (volExact > 0)
    {   // geometric integration error of the cut-cell quadrature
        const real_t volH = volOmega * L3;
        gsInfo << std::setprecision(15)
               << "Volume: |Omega_h| = " << volH << " mm^3, exact = " << volExact
               << " mm^3\n" << std::setprecision(4) << std::scientific
               << "Volume: relative error |V_h - V|/V = "
               << math::abs(volH - volExact) / volExact << ", signed (V_h - V)/V = "
               << (volH - volExact) / volExact << "\n"
               << std::defaultfloat << std::setprecision(6);
    }

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
    gsInfo << "Discretization: numRefine = " << numRefine << ", degree = " << deg
           << ", DoFs = " << nDofs << ", dt = " << dt << ", steps = " << maxSteps << "\n";

    // Initialization biomarkers
    real_t V_T, V_N, V_P, N_dash, N_T, A_p;
    real_t int_N, int_H; // \int_Omega N dOmega  and  \int_Omega H(N-N_th) dOmega

    gsFileManager::mkdir(out);
    // Output: same three files and columns as the body-fitted torus example
    // (biomarkers.csv, computational.csv, computational_summary.csv).
    std::ofstream csvFile, compFile;
    csvFile.open(out+"/biomarkers.csv");
    csvFile << std::setprecision(12); // V_P_relerr needs more than the default 6 digits
    // V_exact = 2 pi^2 R r^2 L^3 [mm^3], V_P_relerr = |V_P - V_exact| / V_exact
    // growth = \int_Omega rho N(1-N) [mm^3/day], dintN_dt = (int_N^n - int_N^{n-1})/dt
    csvFile << "TimeStep,Time,NumDOFs,int_N,int_H,V_P,N_dash,N_T,A_p,N_max,N_min,"
               "V_exact,V_P_relerr,V_T,V_N,growth,dintN_dt\n";
    compFile.open(out+"/computational.csv");
    compFile << std::setprecision(8)
             << "TimeStep,Time,NewtonIts,LinIts,ResidualNorm,StepSolve_s,StepTotal_s,Wall_s\n";

    real_t N_hat_th = par.Nth;  // Cell density threshold defining Omega_T (tumor)
    real_t k_steep = 200.0;  // Steepness factor for Heaviside step function

    // -------------------------------------------------------------------------
    // 3. Assembly on the immersed quadrature points
    //    (this replaces the gsExprAssembler expressions of the body-fitted
    //     example; the background map is the identity on [0,1]^dim, so
    //     meas(G) == 1 and the parametric derivatives are the physical ones)
    // -------------------------------------------------------------------------

    // Scatter a local dense block (indexed by fkQuadBlock::glob) into triplets.
    auto scatter = [](const fkQuadBlock<T> & b, const gsMatrix<T> & Ke,
                      gsSparseEntries<T> & trip)
    {
        const index_t nu = (index_t)b.glob.size();
        for (index_t i = 0; i != nu; ++i)
            for (index_t j = 0; j != nu; ++j)
                if ((T)0 != Ke(i,j))
                    trip.add(b.glob[i], b.glob[j], Ke(i,j));
    };

    // Constant part of the tangent, scaled by \a coeff:
    //    alpha_m * w*w^T  +  alpha_f*gamma*dt * D * grad(w).grad(w)^T
    // \a coeff is 1 on Omega and the tiny regularization weight on the
    // background cells.
    auto addConstant = [&](const fkQuadBlock<T> & b, T coeff, gsSparseEntries<T> & trip)
    {
        const index_t na = b.act.rows(), nu = (index_t)b.glob.size();
        gsMatrix<T> Ke(nu,nu);
        Ke.setZero();
        for (index_t q = 0; q != b.pts.cols(); ++q)
        {
            const T w = b.wts[q] * coeff;
            for (index_t i = 0; i != na; ++i)
            {
                const index_t li = b.loc(i,q);
                for (index_t j = 0; j != na; ++j)
                {
                    T gg = 0;
                    for (short_t d = 0; d != dim; ++d)
                        gg += b.bd(i*dim+d,q) * b.bd(j*dim+d,q);
                    Ke(li, b.loc(j,q)) += w * ( alpha_m * b.bv(i,q) * b.bv(j,q)
                                              + (alpha_f*gamma*dt) * D * gg );
                }
            }
        }
        scatter(b, Ke, trip);
    };

    // Domain residual of  dN/dt = div(D grad N) + rho*N*(1-N):
    //   \int w*dN + \int D*grad(w).grad(N) - \int rho*w*N*(1-N) = 0,
    // integrated over the immersed quadrature points of Omega. When \a trip is
    // non-null, the reaction part of the tangent,
    //   (-alpha_f*gamma*dt) * rho * w * (1 - 2N) * w^T,
    // is accumulated in the same sweep (it needs the very same N at the very
    // same points).
    auto assembleFK = [&](const gsMatrix<T> & Na, const gsMatrix<T> & dNa,
                          gsMatrix<T> & Q, gsSparseEntries<T> * trip)
    {
        Q.setZero(nDofs, 1);
        gsVector<T> gN(dim);
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const index_t na = b.act.rows(), nu = (index_t)b.glob.size();
            gsVector<T> Fe(nu);
            Fe.setZero();
            gsMatrix<T> Ke;
            if (trip) Ke.setZero(nu,nu);

            for (index_t q = 0; q != b.pts.cols(); ++q)
            {
                // N, grad(N) and dN/dt at the quadrature point
                T Nv = 0, dNv = 0;
                gN.setZero();
                for (index_t i = 0; i != na; ++i)
                {
                    const index_t gi = b.act(i,q);
                    Nv  += b.bv(i,q) * Na (gi,0);
                    dNv += b.bv(i,q) * dNa(gi,0);
                    for (short_t d = 0; d != dim; ++d)
                        gN[d] += b.bd(i*dim+d,q) * Na(gi,0);
                }

                const T w = b.wts[q];
                for (index_t i = 0; i != na; ++i)
                {
                    T gg = 0;
                    for (short_t d = 0; d != dim; ++d)
                        gg += b.bd(i*dim+d,q) * gN[d];
                    Fe[b.loc(i,q)] += w * ( b.bv(i,q) * dNv + D * gg
                                          - rho * b.bv(i,q) * Nv * (1.0 - Nv) );
                }

                if (trip)
                {
                    const T cf = w * (-alpha_f*gamma*dt) * rho * (1.0 - 2.0*Nv);
                    for (index_t i = 0; i != na; ++i)
                    {
                        const index_t li = b.loc(i,q);
                        for (index_t j = 0; j != na; ++j)
                            Ke(li, b.loc(j,q)) += cf * b.bv(i,q) * b.bv(j,q);
                    }
                }
            }

            for (index_t i = 0; i != nu; ++i)
                Q(b.glob[i],0) += Fe[i];
            if (trip) scatter(b, Ke, *trip);
        }
    };

    // Quantities of interest, all of them integrals over Omega of a function of
    // N alone -- this is what computeInt() did with the expression evaluator on
    // the body-fitted patch, now evaluated on the immersed quadrature points.
    auto biomarkers = [&](const gsMatrix<T> & coefs,
                          real_t & o_int_N, real_t & o_int_H, real_t & o_V_N,
                          real_t & o_V_P, real_t & o_A_p, real_t & o_growth)
    {
        o_int_N = o_int_H = o_V_N = o_V_P = o_A_p = o_growth = 0;
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const index_t na = b.act.rows();
            for (index_t q = 0; q != b.pts.cols(); ++q)
            {
                T Nv = 0;
                for (index_t i = 0; i != na; ++i)
                    Nv += b.bv(i,q) * coefs(b.act(i,q),0);
                // Indicator of the tumor domain Omega_T (where N >= N_th)
                const T H = 1.0 / (1.0 + math::exp(-k_steep*(Nv - N_hat_th)));
                const T w = b.wts[q];
                o_int_N += w * Nv;              // \int_Omega N
                o_int_H += w * H;               // \int_Omega H(N-N_th) = V_T
                o_V_N   += w * Nv * H;          // total tumor cell volume
                o_V_P   += w;                   // prostate volume
                o_A_p   += w * rho * Nv * (1.0 - Nv) * H; // proliferation activity
                o_growth+= w * rho * Nv * (1.0 - Nv);     // over ALL of Omega (mass balance)
            }
        }
    };

    // -------------------------------------------------------------------------
    // 4. Initial Condition Setup (smooth hyperbolic-tangent tumor seed)
    // -------------------------------------------------------------------------
    gsFunctionExpr<T> N0 = tumorSeed<dim,T>(par);

    // L2-projection of the seed onto the background basis. The background map
    // is the identity, so gsL2Projection's parametric evaluation of N0 IS the
    // physical one here (unlike on the body-fitted torus).
    // There are no eliminated degrees of freedom here (the immersed boundary
    // condition is natural), so the coefficients ARE the unknowns.
    gsMatrix<T> N0coefs;
    gsStopwatch swIC;
    gsInfo << "IC L2-projection error: "
           << gsL2Projection<T>::project(dbasis, mp, N0, N0coefs) << "\n";
    const real_t tIC = swIC.stop();
    // Relative L2(Omega) error ||N_h(0) - N0|| / ||N0|| on the immersed
    // quadrature (unit-free, so the same number as in the body-fitted code;
    // the projection above itself is over the whole background box).
    real_t icRelErr = 0;
    {
        real_t e2 = 0, n2 = 0;
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const gsMatrix<T> f = N0.eval(b.pts);
            for (index_t q = 0; q != b.pts.cols(); ++q)
            {
                T Nv = 0;
                for (index_t i = 0; i != b.act.rows(); ++i)
                    Nv += b.bv(i,q) * N0coefs(b.act(i,q),0);
                e2 += b.wts[q] * (Nv - f(0,q)) * (Nv - f(0,q));
                n2 += b.wts[q] * f(0,q) * f(0,q);
            }
        }
        icRelErr = math::sqrt(e2 / n2);
        gsInfo << "IC L2(Omega) relative projection error: " << icRelErr << "\n";
    }

    gsMatrix<T> Nnew, dNnew, Nold, dNold;
    gsMatrix<T> Nalpha, dNalpha, dNupdate, Q;

    Nold = N0coefs;
    dNold.setZero(nDofs, 1);

    gsInfo << "Linear solver: " << linSolver
           << (linSolver == "direct" ? "" : " (CG tol " + util::to_string(cgTol) + ")") << "\n";
    gsSparseMatrix<T> K, K_const;

    // -------------------------------------------------------------------------
    // 5. Pre-assemble Constant Terms Matrix
    // -------------------------------------------------------------------------
    {
        gsSparseEntries<T> trip;
        for (size_t e = 0; e != blocks.size(); ++e)
            addConstant(blocks[e], (T)1, trip);

        // Tiny full-background regularization: keeps the tangent nonsingular for
        // the degrees of freedom whose support never touches {phi<0}. Exactly
        // as in the immersed Poisson examples, it is added to the MATRIX only,
        // so it contributes nothing to the residual and those dofs simply stay
        // at their initial value.
        fkQuadBlock<T> bg;
        for (auto it = tbsPtr->domain()->beginAll(); it != tbsPtr->domain()->endAll(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), qpts, qwts);
            fkBuildBlock(*tbsPtr, qpts, qwts, bg);
            addConstant(bg, (T)1e-10, trip);
        }

        K_const.resize(nDofs, nDofs);
        K_const.setFrom(trip);
    }

    // -------------------------------------------------------------------------
    // 6. ParaView output of the cell density N (one time step per frame).
    //    Primary output: N as a volumetric FIELD on a uniform sample grid of the
    //    background box, written together with the level set phi in the same
    //    .vts file, so that one ParaView "Clip -> Scalars: phi -> Value 0"
    //    restricts it to the prostate. Optionally (--plotPoints) the same N
    //    sampled at the immersed quadrature points, which is the honest picture
    //    of where the weak form is actually integrated.
    // -------------------------------------------------------------------------
    std::unique_ptr<gsParaviewCollection> colField, colPts, colMesh;
    gsMatrix<T> cloud, gridPts, gridPhys, phiGrid, Nvals;
    gsVector<index_t> res(3);

    if (plot)
    {
        gsFileManager::mkdir(out + "/N");
        colField.reset(new gsParaviewCollection(out + "/N/N"));

        // Uniform sample grid of [0,1]^3, x fastest then y then z (VTK order).
        res.setConstant(math::max(fcm.plotRes, (index_t)2));
        gridPts.resize(3, res.prod());
        index_t c = 0;
        for (index_t k = 0; k != res[2]; ++k)
        for (index_t j = 0; j != res[1]; ++j)
        for (index_t i = 0; i != res[0]; ++i, ++c)
        {
            gridPts(0,c) = (T)i / (T)(res[0]-1);
            gridPts(1,c) = (T)j / (T)(res[1]-1);
            gridPts(2,c) = (T)k / (T)(res[2]-1);
        }
        // phi is static, so it is evaluated once and reused by every frame
        // (each evaluation is a BVH query against the triangle mesh).
        gsStopwatch swPhi;
        phi.eval_into(gridPts, phiGrid);
        phiGrid *= fcm.toPhys;               // distance in mesh units
        gridPhys = gridPts;                  // written coordinates: mesh units
        fkToPhys(fcm, gridPhys);
        gsInfo << "Plot grid: " << res[0] << "^3 = " << res.prod()
               << " points, level set evaluated in " << swPhi.stop() << " s.\n";
    }
    if (plot && fcm.plotMesh)
    {
        // The background grid the immersed quadrature lives on. The Poisson
        // examples rewrite it at every refinement level; here it is static for
        // the whole simulation, so it is written once as the single part of
        // its own collection.
        gsFileManager::mkdir(out + "/background");
        colMesh.reset(new gsParaviewCollection(out + "/background/background"));
        gsMesh<T> bgMesh(dbasis.basis(0));
        for (size_t v = 0; v != bgMesh.numVertices(); ++v)
            for (index_t d = 0; d != 3; ++d)
                (*bgMesh.vertices()[v])[d] = (*bgMesh.vertices()[v])[d] * fcm.toPhys
                                           + fcm.shift[d];
        gsWriteParaview(bgMesh, out + "/background/grid");
        colMesh->addPart("grid.vtp", 0);
    }
    if (plot && fcm.plotPoints)
    {
        cloud.resize(dim+1, nQuadInterior + nQuadCut);
        gsFileManager::mkdir(out + "/Npts");
        colPts.reset(new gsParaviewCollection(out + "/Npts/Npts"));
    }

    // N sampled on the background grid -> one .vts (N and phi) per time step.
    auto writeField = [&](const gsMatrix<T> & coefs, index_t frame, real_t time)
    {
        typename gsGeometry<T>::uPtr geo = dbasis.basis(0).makeGeometry(coefs);
        geo->eval_into(gridPts, Nvals);
        for (index_t i = 0; i != Nvals.cols(); ++i) // see fkWriteField
            if (phiGrid(0,i) >= 0) Nvals(0,i) = 0;
        const std::string fn = "N_" + util::to_string(frame);
        fkWriteField<T>(out + "/N/" + fn + ".vts", res, gridPhys, Nvals, phiGrid);
        colField->addPart(fn + ".vts", time);
    };

    // N sampled at the immersed quadrature points (4 x N: xyz + value).
    auto writeCloud = [&](const gsMatrix<T> & coefs, index_t frame, real_t time)
    {
        index_t c = 0;
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const index_t na = b.act.rows();
            for (index_t q = 0; q != b.pts.cols(); ++q, ++c)
            {
                T Nv = 0;
                for (index_t i = 0; i != na; ++i)
                    Nv += b.bv(i,q) * coefs(b.act(i,q),0);
                cloud.col(c).head(dim) = b.pts.col(q);
                cloud(dim,c) = Nv;
            }
        }
        gsMatrix<T> xyz = cloud.topRows(dim);
        fkToPhys(fcm, xyz);
        cloud.topRows(dim) = xyz;
        const std::string fn = "Npts_" + util::to_string(frame);
        gsWriteParaviewPoints(cloud, out + "/Npts/" + fn);
        colPts->addPart(fn + ".vtp", time);
    };

    // frame i of the collection = row i-1 of the csv file
    auto writeOut = [&](const gsMatrix<T> & coefs, index_t frame, real_t time)
    {
        if (!plot) return;
        writeField(coefs, frame, time);
        if (fcm.plotPoints) writeCloud(coefs, frame, time);
    };

    writeOut(Nold, 0, 0.0); // initial condition

    // DEBUG: max/min of N over the immersed quadrature points (i.e. inside the
    // prostate). The coefficients themselves are not used, since dofs whose
    // support misses {phi<0} stay frozen at their initial value.
    auto rangeN = [&](const gsMatrix<T> & coefs, T & Nmax, T & Nmin)
    {
        Nmax = -std::numeric_limits<T>::max();
        Nmin =  std::numeric_limits<T>::max();
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            for (index_t q = 0; q != b.pts.cols(); ++q)
            {
                T Nv = 0;
                for (index_t i = 0; i != b.act.rows(); ++i)
                    Nv += b.bv(i,q) * coefs(b.act(i,q),0);
                Nmax = math::max(Nmax, Nv);
                Nmin = math::min(Nmin, Nv);
            }
        }
    };
    T Nmax, Nmin;
    rangeN(Nold, Nmax, Nmin);
    gsInfo << "Initial condition: max N = " << Nmax << ", min N = " << Nmin << "\n";
    const T icNmax = Nmax, icNmin = Nmin;

    // Biomarkers of the coefficients Nold (row n of biomarkers.csv)
    real_t growth = 0, int_N_prev = 0;
    auto writeBio = [&](index_t n, real_t time)
    {
        // Quantities of interest of the forward simulation, integrated over the
        // immersed domain with the very same quadrature as the weak form.
        biomarkers(Nold, int_N, int_H, V_N, V_P, A_p, growth);
        gsInfo << "  int_Omega N dOmega          = " << int_N << "\n"
               << "  int_Omega H(N-N_th) dOmega  = " << int_H << "\n";
        rangeN(Nold, Nmax, Nmin);
        gsInfo << "  max N = " << Nmax << ", min N = " << Nmin << "\n";

        V_T    = int_H;         // Tumor volume
        N_dash = V_N / V_T;     // mean normalized tumor cell density
        N_T    = N_dash * V_T / V_P; // total tumor index

        // Integrals are computed in box units; the csv reports mm^3.
        // N_dash and N_T are ratios and need no conversion.
        const real_t vol = L3;
        csvFile << n << "," << time << "," << nActive << ","
                << int_N*vol << "," << int_H*vol << ","
                << V_P*vol << "," << N_dash << "," << N_T << "," << A_p*vol << ","
                << Nmax << "," << Nmin << ","
                << volExact << "," << math::abs(V_P*vol - volExact) / volExact << ","
                << V_T*vol << "," << V_N*vol << "," << growth*vol << ",";
        if (n > 0) csvFile << (int_N - int_N_prev)*vol / dt;
        csvFile << "\n";
        csvFile.flush();
        int_N_prev = int_N;
    };
    writeBio(0, 0.0);
    const real_t V_P_relerr = math::abs(V_P*L3 - volExact) / volExact;
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

            // Domain Residual + reaction Tangent Assembly, on the immersed
            // quadrature points (both come out of the same sweep because they
            // need N at the same points).
            gsSparseEntries<T> trip;
            assembleFK(Nalpha, dNalpha, Q, &trip);

            resNorm = Q.norm();
            if (it > 0 && resNorm < tol)
            {
                gsInfo << "  Converged in " << it << " iterations.\n";
                its = it;
                break;
            }

            K.resize(nDofs, nDofs);
            K.setFrom(trip);
            K += K_const;

            // Condition number of the first system matrix (--cond)
            if (cond && step == 0 && it == 0)
            {
                // The full matrix also contains the inactive background dofs,
                // held only by the 1e-10 regularization: its kappa measures that
                // regularization. The ACTIVE block (dofs whose support meets
                // Omega) is the one that carries the cut-cell conditioning.
                gsSparseEntries<T> ta;
                for (index_t k = 0; k < K.outerSize(); ++k)
                    for (typename gsSparseMatrix<T>::InnerIterator itK(K, k); itK; ++itK)
                        if (activeIdx[itK.row()] >= 0 && activeIdx[itK.col()] >= 0)
                            ta.add(activeIdx[itK.row()], activeIdx[itK.col()], itK.value());
                gsSparseMatrix<T> Ka(nActive, nActive);
                Ka.setFrom(ta);
                gsInfo << "[full matrix incl. inactive background dofs]\n";
                condCsv = fkReportCond(K);
                gsInfo << "[active dofs only]\n";
                std::string ca = fkReportCond(Ka);
                // prefix the active-block keys: cond_K -> cond_active_K, ...
                for (size_t pos = 0; (pos = ca.find("cond_", pos)) != std::string::npos; pos += 12)
                    ca.replace(pos, 5, "cond_active_");
                condCsv += ca;
            }

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
        Nold  = Nnew;
        dNold = dNnew;

        writeOut(Nold, step+1, (step+1)*dt);

        writeBio(step+1, (step+1)*dt);

        itsTotal    += its;
        tSolveTotal += tSolve;
        compFile << step+1 << "," << (step+1)*dt << "," << its << "," << linIts << "," << resNorm << ","
                 << tSolve << "," << swStep.stop() << "," << swTotal.stop() << "\n";
        compFile.flush();
    }
    csvFile.close();
    compFile.close();
    if (plot)
    {
        colField->save();
        if (fcm.plotPoints) colPts->save();
        if (fcm.plotMesh) colMesh->save();
    }

    const real_t tLoop = swLoop.stop(), tTotal = swTotal.stop();
    std::ofstream sumFile((out + "/computational_summary.csv").c_str());
    sumFile << std::setprecision(10) << "quantity,value\n"
            << "method,immersed\n"
            << "numRefine," << numRefine << "\nquadDepth," << fcm.quadDepth << "\n"
            << "degree," << deg << "\n"
            << "elements," << nInteriorCells + nCutCells << "\n"
            << "elements_interior," << nInteriorCells << "\nelements_cut," << nCutCells << "\n"
            // uniform background cubes: h = edge = L / 2^r (cut cells have the
            // same cell size; only the part inside Omega is smaller)
            << "h_max_mm," << parPhys.L / math::pow(2.0, (real_t)numRefine) << "\n"
            << "h_min_mm," << parPhys.L / math::pow(2.0, (real_t)numRefine) << "\n"
            << "h_mean_mm," << parPhys.L / math::pow(2.0, (real_t)numRefine) << "\n"
            << "DoFs," << nDofs << "\n"
            << "DoFs_active," << nActive << "\n"
            << "quad_points_interior," << nQuadInterior << "\n"
            << "quad_points_cut," << nQuadCut << "\n"
            << "quad_points_total," << nQuadInterior + nQuadCut << "\n"
            << "quad_cache_MB," << (cacheBytes >> 20) << "\n"
            << "time_quad_build_s," << tQuad << "\n"
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
    gsInfo << "Times: quadrature " << tQuad << " s, setup " << tSetup << " s, time loop "
           << tLoop << " s (" << (maxSteps ? tLoop / maxSteps : 0) << " s/step, Newton solve "
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

    for (int i = 1; i+1 != argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-f" || arg == "--file") fn = argv[i+1];
    }

    // The xml file is only used for rho_inf and the Newton tolerance (id 4);
    // its model parameters (id 3, box units) are ignored on purpose, exactly
    // as in the body-fitted torus example.
    gsFileData<> fd(fn);
    gsInfo << "Loaded configuration file: " << fd.lastPath() << "\n";
    gsOptionList CHopt, TIMEopt, Aopt;
    fd.getId(3, CHopt);
    fd.getId(4, TIMEopt);
    fd.getId(5, Aopt);

    // Physical model parameters (mm, day) and time stepping: SAME defaults
    // as fischerkolmogrov_torus_example_bodyfitted.cpp.
    fkParams<real_t> par;
    real_t  dt        = 5.0;
    index_t maxSteps  = 146;
    index_t numRefine = 5;   // background box: h = L/2^5 = 1.6 mm for L = 50
    bool    plot      = false;

    torusParams tp;
    fcmParams   fcm;
    bool noPlotMesh = !fcm.plotMesh;

    gsCmdLine cmd("Fisher-Kolmogorov Generalized-Alpha Solver on an IMMERSED solid torus (forward problem)");
    cmd.addReal("t", "dt", "Time step size [day]", dt);
    cmd.addInt("N", "Nsteps", "Number of time steps", maxSteps);
    cmd.addInt("r", "numRefine", "Number of uniform refinement steps of the background box", numRefine);
    cmd.addString("f", "file", "Input XML setup file (rho_inf, tol)", fn);
    cmd.addString( "o", "output", "Output directory", out);
    cmd.addSwitch("plot", "Write ParaView files of N for every time step", plot);
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
    cmd.addReal("", "xc", "x-coordinate of the tumor centre [mm] (default: on the core circle at 45 deg)", par.xc);
    cmd.addReal("", "yc", "y-coordinate of the tumor centre [mm] (default: on the core circle at 45 deg)", par.yc);
    cmd.addReal("", "zc", "z-coordinate of the tumor centre [mm] (default: on the core circle at 45 deg)", par.zc);
    cmd.addReal("", "torusR", "Major radius of the torus [box units]", tp.R);
    cmd.addReal("", "torusr", "Minor (tube) radius of the torus [box units]", tp.r);
    cmd.addReal("", "torusCx", "x-coordinate of the torus centre [box units]", tp.cx);
    cmd.addReal("", "torusCy", "y-coordinate of the torus centre [box units]", tp.cy);
    cmd.addReal("", "torusCz", "z-coordinate of the torus centre [box units]", tp.cz);
    // Immersed / finite-cell quadrature
    cmd.addInt("e", "degreeElevation", "Degree elevation of the (linear) background basis; 1 = quadratic", fcm.degElevate);
    cmd.addInt("", "quadDepth", "Octree subdivision depth of the cut-cell quadrature", fcm.quadDepth);
    cmd.addString("", "indicator", "Adaptive indicator: uniform, fallback or integralChange", fcm.indicator);
    cmd.addReal("", "indicatorTol", "integralChange acceptance tolerance", fcm.indicatorTol);
    cmd.addInt("", "plotRes", "Sample points per direction of the plotted N field", fcm.plotRes);
    cmd.addSwitch("plotPoints", "Additionally dump N at the immersed quadrature points", fcm.plotPoints);
    cmd.addSwitch("noPlotMesh", "Do not write the background grid", noPlotMesh);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    fcm.plotMesh = !noPlotMesh;

    // -------------------------------------------------------------------------
    // Geometry: exact torus level set in the background box [0,1]^3
    // -------------------------------------------------------------------------
    GISMO_ENSURE(tp.r > 0 && tp.R > tp.r, "Need R > r > 0 for a torus.");
    const real_t L = par.L;
    const real_t volExact = 2.0 * EIGEN_PI * EIGEN_PI * tp.R * tp.r * tp.r * L * L * L;
    gsInfo << "Torus (box units): centre (" << tp.cx << "," << tp.cy << "," << tp.cz
           << "), R = " << tp.R << ", r = " << tp.r << "\n"
           << "Torus (1 box unit = " << L << " mm): R = " << tp.R * L << " mm, r = "
           << tp.r * L << " mm, outer diameter " << 2 * (tp.R + tp.r) * L
           << " mm, exact volume 2 pi^2 R r^2 = " << volExact << " mm^3\n";
    if (tp.cx - tp.R - tp.r < 0 || tp.cx + tp.R + tp.r > 1 ||
        tp.cy - tp.R - tp.r < 0 || tp.cy + tp.R + tp.r > 1 ||
        tp.cz - tp.r < 0 || tp.cz + tp.r > 1)
        gsWarn << "The torus does not fit inside [0,1]^3.\n";

    gsTorusLevelSet<real_t> impl_fun(tp);

    // Default tumour centre: on the core circle at 45 deg (same as body-fitted)
    if (math::isnan(par.xc) && math::isnan(par.yc) && math::isnan(par.zc))
    {
        par.xc = (tp.cx + tp.R * math::sqrt(0.5)) * L;
        par.yc = (tp.cy + tp.R * math::sqrt(0.5)) * L;
        par.zc = tp.cz * L;
    }
    GISMO_ENSURE(!math::isnan(par.xc) && !math::isnan(par.yc) && !math::isnan(par.zc),
                 "Give all three of --xc --yc --zc, or none.");
    if (par.R + 2 * par.sigma > tp.r * L)
        gsWarn << "Tumour seed R + 2 sigma = " << par.R + 2 * par.sigma
               << " mm exceeds the tube radius r = " << tp.r * L
               << " mm: the seed is cut by the torus surface.\n";
    GISMO_ENSURE(par.amp > 0 && par.amp <= 1, "The seed amplitude must lie in (0,1].");

    // Regime diagnostics (identical to the body-fitted example, h isotropic)
    {
        const real_t T   = dt * maxSteps;
        const real_t c   = 2 * math::sqrt(par.D * par.rho);
        const real_t ell = math::sqrt(par.D / par.rho);
        const real_t h   = L / math::pow(2.0, (real_t)numRefine);
        gsInfo << "Regime over T = " << T << " days:\n"
               << "  rho*T          = " << par.rho * T << "  (growth factor exp(rho T) = "
               << math::exp(par.rho * T) << ")\n"
               << "  D/R^2          = " << par.D / (par.R * par.R) << " 1/day  vs  rho = "
               << par.rho << " 1/day\n"
               << "  front speed c  = " << c << " mm/day  (c*T = " << c * T << " mm)\n"
               << "  diffusion len. = sqrt(D T) = " << math::sqrt(par.D * T) << " mm\n"
               << "  front width    ~ 4 sqrt(D/rho) = " << 4 * ell << " mm,  background h = "
               << h << " mm\n";
        if (par.rho * T < 1)
            gsWarn << "rho*T < 1: almost no proliferation within the horizon.\n";
        if (par.D / (par.R * par.R) > par.rho)
            gsWarn << "D/R^2 > rho: diffusion-dominated, the tumour will be diluted "
                      "(max N and V_T decrease) instead of growing.\n";
        if (h > 2 * 4 * ell)
            gsWarn << "The mesh is much coarser than the front width: refine (-r).\n";
    }

    // Background box [0,1]^3: the identity map, parameter space == physical
    // space (box units). The torus enters only through impl_fun.
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    solve<3, real_t>(mp, impl_fun, bc, CHopt, TIMEopt, Aopt,
                     dt, maxSteps, numRefine, par, fcm, out, plot, volExact,
                     linSolver, cgMaxIt, cgTol, cond);

    gsInfo << "Exact torus volume (compare with V_P in the csv): " << volExact << " mm^3\n";
    return EXIT_SUCCESS;
}
