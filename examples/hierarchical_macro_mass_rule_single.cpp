/** @file hierarchical_macro_mass_rule_single.cpp

    @brief Nonlinear moment-fitting quadrature rule for one atomic cell of a
    hierarchical THB-spline domain.

    For a selected gsHTree leaf (= one atomic cell = one macro region):
      1. Builds a local target basis T_r = S^{2p}_{p-1}(D_r) (degree-2p, C^{p-1}
         tensor B-spline over the macro's knot lines).
      2. Computes exact target moments b_a = integral_{D_r} Phi_a using the
         B-spline integral formula (no quadrature needed for moments).
      3. Initialises n_q = ceil(rank/3) + safety quadrature points, one per
         element sub-grid, with log-weight parametrisation (always positive).
      4. Solves the nonlinear moment system F_a = sum_q exp(s_q) Phi_a(xi_q,eta_q) - b_a = 0
         via Levenberg-Marquardt, using element-local coordinates to enforce
         containment without projection.
      5. Validates the result against all mass product moments B_i B_j on D_r.

    Usage:
        ./hierarchical_macro_mass_rule_single [-c N] [-s S] [-i maxIter] [-t tol]

    Options:
        -c N     Atomic cell index to solve (default: smallest cell with ~8 elements)
        -s S     Safety margin: extra points above ceil(rank/3) (default: 8)
        -i N     Max LM iterations (default: 200)
        -t tol   Convergence tolerance on ||F|| (default: 1e-10)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <unordered_map>
#include <cmath>

using namespace gismo;

// ============================================================
// Shared infrastructure (same as hierarchical_macro_mass_system.cpp)
// ============================================================

struct ElementBox
{
    gsVector<real_t> lower;
    gsVector<real_t> upper;
};

struct AtomicCell
{
    index_t            level;
    gsVector<index_t>  lowerIndex;
    gsVector<index_t>  upperIndex;
    gsVector<real_t>   lowerParam;
    gsVector<real_t>   upperParam;
    std::vector<ElementBox> elements;
};

template<short_t d, class T>
std::vector<AtomicCell> collectAtomicCells(const gsHTensorBasis<d,T> & thb)
{
    GISMO_ASSERT( !thb.manualLevels(),
        "collectAtomicCells: manualLevels bases not supported here." );

    std::vector<AtomicCell> cells;
    auto leaf = thb.tree().beginLeafIterator();

    for (; leaf.good(); leaf.next())
    {
        AtomicCell cell;
        cell.level      = leaf.level();
        cell.lowerIndex = leaf.lowerCorner();
        cell.upperIndex = leaf.upperCorner();
        cell.lowerParam.resize(d);
        cell.upperParam.resize(d);

        for (short_t dd = 0; dd < d; ++dd)
        {
            const gsKnotVector<T>& kv =
                thb.tensorLevel(cell.level).component(dd).knots();
            cell.lowerParam[dd] = kv( cell.lowerIndex[dd] );
            cell.upperParam[dd] = kv( cell.upperIndex[dd] );
        }

        // Enumerate all knot-span elements inside this leaf
        gsVector<index_t> alpha(d);
        alpha = cell.lowerIndex;

        bool more = true;
        while (more)
        {
            ElementBox elem;
            elem.lower.resize(d);
            elem.upper.resize(d);
            for (short_t dd = 0; dd < d; ++dd)
            {
                const gsKnotVector<T>& kv =
                    thb.tensorLevel(cell.level).component(dd).knots();
                elem.lower[dd] = kv( alpha[dd]     );
                elem.upper[dd] = kv( alpha[dd] + 1 );
            }
            cell.elements.push_back(elem);

            more = false;
            for (short_t dd = d - 1; dd >= 0; --dd)
            {
                if (alpha[dd] < cell.upperIndex[dd] - 1)
                { ++alpha[dd]; more = true; break; }
                alpha[dd] = cell.lowerIndex[dd];
            }
        }
        cells.push_back(cell);
    }
    return cells;
}

// ============================================================
// Target basis: S^{2p}_{p-1}(D_r) over the macro's breakpoints
// ============================================================

/**
 * Build a degree-2p, C^{p-1} tensor B-spline basis over the atomic cell.
 * Internal knot multiplicity = p+1 = (2p) - (p-1).
 * 1D dimension = (p+1)*n + p.
 *
 * Uses gsKnotVector(uKnots, degree, regularity) constructor which sets
 *   boundary multiplicity = degree+1
 *   internal multiplicity = degree - regularity
 */
template<class T>
gsTensorBSplineBasis<2,T>
buildTargetBasis(const AtomicCell & cell,
                 const gsHTensorBasis<2,T> & thb)
{
    const int p = (int)thb.degree(0);  // assume isotropic
    const int q = 2 * p;               // target degree

    gsKnotVector<T> kvs[2];
    for (int dd = 0; dd < 2; ++dd)
    {
        const gsKnotVector<T>& kv_lev =
            thb.tensorLevel(cell.level).component(dd).knots();

        // Collect unique breakpoints in parametric space
        std::vector<T> uKnots;
        for (index_t i = cell.lowerIndex[dd]; i <= cell.upperIndex[dd]; ++i)
            uKnots.push_back( kv_lev( i ) );  // unique-knot access: kv(i)

        // gsKnotVector(uKnots, degree, regularity):
        //   boundary mult = q+1, internal mult = q - regularity = p+1
        kvs[dd] = gsKnotVector<T>(uKnots, q, p - 1);
    }
    return gsTensorBSplineBasis<2,T>(kvs[0], kvs[1]);
}

// ============================================================
// Exact target moments via B-spline integral formula
// ============================================================

/**
 * Compute b_k = integral_{D_r} Phi_k  using the exact formula
 *   integral N_{i,q} = (tau[i+q+1] - tau[i]) / (q+1)
 * where tau is the REPEATED knot vector.
 *
 * Tensor product:  b[i + j*n0] = integral_x[i] * integral_y[j].
 */
template<class T>
gsVector<T>
computeTargetMoments(const gsTensorBSplineBasis<2,T> & target)
{
    const index_t q  = target.degree(0);  // = 2p
    const index_t n0 = target.size(0);    // 1D count in x
    const index_t n1 = target.size(1);    // 1D count in y

    // Repeated knot access: kv[i]  (not kv(i) which is unique-knot)
    const gsKnotVector<T>& kv0 = target.component(0).knots();
    const gsKnotVector<T>& kv1 = target.component(1).knots();

    gsVector<T> m0(n0), m1(n1);
    for (index_t i = 0; i < n0; ++i)
        m0[i] = ( kv0[i + q + 1] - kv0[i] ) / T(q + 1);
    for (index_t j = 0; j < n1; ++j)
        m1[j] = ( kv1[j + q + 1] - kv1[j] ) / T(q + 1);

    // Ordering: index(i,j) = i + j*n0  (x fastest, per gsTensorBSplineBasis::index)
    gsVector<T> b(n0 * n1);
    for (index_t j = 0; j < n1; ++j)
        for (index_t i = 0; i < n0; ++i)
            b[i + j * n0] = m0[i] * m1[j];

    return b;
}

// ============================================================
// MacroRule: quadrature points in element-local coordinates
// ============================================================

/**
 * Quadrature rule parametrised for the LM solver.
 *
 * Each point q belongs to a fixed element (elemOf[q]).  The unknowns are:
 *   xi[q]   ∈ (eps, 1-eps)  — element-local x coordinate
 *   eta[q]  ∈ (eps, 1-eps)  — element-local y coordinate
 *   logw[q] ∈ R             — log(weight), ensures w > 0
 *
 * Physical positions and weights are derived via applyState.
 *
 * State vector z layout: z[3q]=xi, z[3q+1]=eta, z[3q+2]=logw.
 */
struct MacroRule
{
    index_t                nq;
    std::vector<index_t>   elemOf;
    gsVector<real_t>       xi;
    gsVector<real_t>       eta;
    gsVector<real_t>       logw;
    gsMatrix<real_t>       pts;    // 2 × nq, derived
    gsVector<real_t>       wts;    // nq, derived (= exp(logw))
};

// Rebuild pts and wts from (xi, eta, logw) with optional clamping
static const real_t kEps = real_t(1e-7);

template<class T>
void applyState(MacroRule & rule, const gsVector<T> & z,
                const AtomicCell & cell)
{
    const index_t nq = rule.nq;
    for (index_t q = 0; q < nq; ++q)
    {
        rule.xi[q]   = std::max(kEps, std::min(T(1) - kEps, z[3*q    ]));
        rule.eta[q]  = std::max(kEps, std::min(T(1) - kEps, z[3*q + 1]));
        rule.logw[q] = z[3*q + 2];

        const index_t e  = rule.elemOf[q];
        const T wx = cell.elements[e].upper[0] - cell.elements[e].lower[0];
        const T wy = cell.elements[e].upper[1] - cell.elements[e].lower[1];
        rule.pts(0, q) = cell.elements[e].lower[0] + rule.xi[q]  * wx;
        rule.pts(1, q) = cell.elements[e].lower[1] + rule.eta[q] * wy;
        rule.wts[q]    = std::exp(rule.logw[q]);
    }
}

template<class T>
gsVector<T> extractState(const MacroRule & rule)
{
    gsVector<T> z(3 * rule.nq);
    for (index_t q = 0; q < rule.nq; ++q)
    {
        z[3*q    ] = rule.xi[q];
        z[3*q + 1] = rule.eta[q];
        z[3*q + 2] = rule.logw[q];
    }
    return z;
}

/**
 * Distribute nq points among n_e elements.
 * Each element gets k = floor(nq/ne) or floor(nq/ne)+1 points,
 * placed at Gauss-Legendre node positions (better initial residual than
 * uniform grid) on a nx × ny tensor grid in (eps, 1-eps)^2.
 * Initial weight for element e: area(e) / k_e.
 */
template<class T>
MacroRule initMacroRule(const AtomicCell & cell, index_t nq)
{
    MacroRule rule;
    rule.nq = nq;

    const index_t ne = (index_t)cell.elements.size();
    rule.elemOf.resize(nq);
    rule.xi.resize(nq);
    rule.eta.resize(nq);
    rule.logw.resize(nq);
    rule.pts.resize(2, nq);
    rule.wts.resize(nq);

    // Distribute points among elements
    std::vector<index_t> ptsPerElem(ne, nq / ne);
    for (index_t e = 0; e < nq % ne; ++e)
        ++ptsPerElem[e];

    // Cache 1D GL nodes in [0,1] by point count (avoid repeated computation)
    std::map<index_t, std::vector<T>> glCache;
    auto glNodes1D = [&](index_t n) -> const std::vector<T>&
    {
        auto it = glCache.find(n);
        if (it != glCache.end()) return it->second;
        gsGaussRule<T> rule1d(n);
        gsMatrix<T> pts1d;
        gsVector<T> wts1d;
        rule1d.mapTo(T(0), T(1), pts1d, wts1d);  // nodes in [0,1]
        std::vector<T>& nodes = glCache[n];
        nodes.resize(n);
        for (index_t i = 0; i < n; ++i) nodes[i] = pts1d(0, i);
        return nodes;
    };

    index_t qidx = 0;
    for (index_t e = 0; e < ne; ++e)
    {
        const index_t k = ptsPerElem[e];
        if (k == 0) continue;

        const T area_e =
            (cell.elements[e].upper[0] - cell.elements[e].lower[0]) *
            (cell.elements[e].upper[1] - cell.elements[e].lower[1]);
        const T w_e = area_e / T(k);

        // Tensor grid: nx × ny, take first k points; use GL nodes for positions
        const index_t nx = (index_t)std::ceil(std::sqrt((double)k));
        const index_t ny = (k + nx - 1) / nx;

        const std::vector<T>& gx = glNodes1D(nx);
        const std::vector<T>& gy = glNodes1D(ny);

        index_t placed = 0;
        for (index_t jy = 0; jy < ny && placed < k; ++jy)
        {
            for (index_t jx = 0; jx < nx && placed < k; ++jx, ++placed, ++qidx)
            {
                rule.elemOf[qidx] = e;
                // GL position in [eps, 1-eps]
                rule.xi[qidx]  = kEps + (T(1) - T(2)*kEps) * gx[jx];
                rule.eta[qidx] = kEps + (T(1) - T(2)*kEps) * gy[jy];
                rule.logw[qidx] = std::log(w_e);
            }
        }
    }
    GISMO_ASSERT(qidx == nq, "Point distribution error");

    // Fill pts and wts from (xi, eta, logw)
    const gsVector<T> z = extractState<T>(rule);
    applyState<T>(rule, z, cell);
    return rule;
}

// ============================================================
// Residual and Jacobian
// ============================================================

/**
 * Compute F and J for the moment system:
 *
 *   F[k] = sum_q exp(logw_q) * Phi_k(x_q, y_q) - b[k]
 *
 * State variables: z = [xi_0, eta_0, logw_0, xi_1, ...]
 *
 *   dF[k]/d xi_q   = w_q * (d_x Phi_k)(x_q,y_q) * wx_e
 *   dF[k]/d eta_q  = w_q * (d_y Phi_k)(x_q,y_q) * wy_e
 *   dF[k]/d logw_q = w_q * Phi_k(x_q,y_q)
 *
 * Row index in J corresponds to the global target basis function index
 * returned by active_into (scatter assignment).
 *
 * Derivative layout for gsBasis::deriv_into (d=2):
 *   derivs(2*k,   q) = d/dx of active function k at point q
 *   derivs(2*k+1, q) = d/dy of active function k at point q
 */
template<class T>
void evalResJac(const MacroRule                  & rule,
                const gsTensorBSplineBasis<2,T>  & target,
                const AtomicCell                 & cell,
                const gsVector<T>                & b,
                gsVector<T>                      & F,
                gsSparseMatrix<T>                & J,
                double* t_basis = nullptr,   // out: time in basis eval (s)
                double* t_assemble = nullptr) // out: time in assembly (s)
{
    const index_t nq = rule.nq;
    const index_t nt = b.size();

    F.setZero(nt);
    gsSparseEntries<T> triplets;

    // Evaluate target basis at all quadrature points
    gsMatrix<index_t> acts;
    gsMatrix<T>       vals, derivs;

    gsStopwatch sw;
    sw.restart();
    target.active_into(rule.pts, acts);   // (na × nq)
    target.eval_into  (rule.pts, vals);   // (na × nq)
    target.deriv_into (rule.pts, derivs); // (2*na × nq)
    if (t_basis) *t_basis += sw.stop();

    const index_t na = acts.rows();

    sw.restart();
    for (index_t q = 0; q < nq; ++q)
    {
        const T wq = rule.wts[q];
        const index_t e  = rule.elemOf[q];
        const T wx = cell.elements[e].upper[0] - cell.elements[e].lower[0];
        const T wy = cell.elements[e].upper[1] - cell.elements[e].lower[1];

        for (index_t k = 0; k < na; ++k)
        {
            const index_t row = acts(k, q);       // global target basis index
            const T phi    = vals(k, q);
            const T dphi_x = derivs(2*k,     q);
            const T dphi_y = derivs(2*k + 1, q);

            F[row] += wq * phi;

            // Chain rule: d/d xi  = d/dx * wx,  d/d logw = w * phi
            triplets.add(row, 3*q    , wq * dphi_x * wx);
            triplets.add(row, 3*q + 1, wq * dphi_y * wy);
            triplets.add(row, 3*q + 2, wq * phi);
            // J(row, 3*q    ) = wq * dphi_x * wx;
            // J(row, 3*q + 1) = wq * dphi_y * wy;
            // J(row, 3*q + 2) = wq * phi;
        }
    }
    F -= b;
    if (t_assemble) *t_assemble += sw.stop();
    J.clear();
    J.resize(nt, 3*nq);
    J.setFrom(triplets);
}

// ============================================================
// Levenberg-Marquardt solver (SVD-based)
// ============================================================

/**
 * LM solver for the nonlinear moment system using dense normal equations.
 *
 * At every accepted step we recompute J and form the normal equations
 *
 *   (J^T J + lambda I) dz = -J^T F
 *
 * solved via Cholesky (LDLT).  This is O(nz^3) per step vs. O(nt^2 * nz)
 * for thin-SVD and is ~100x faster in practice because nz (= 3*nq = 414)
 * is modest.
 *
 * A lambda floor of  max(diag(J^T J)) * 1e-12  prevents lambda from
 * collapsing to 1e-20 (which made the original normal-equations solver
 * stall).  This is equivalent to the sqrt(eps)*s_max^2 floor used in the
 * SVD formulation.
 *
 * @return true if ||F|| < tol within maxIter iterations
 */
template<class T>
bool solveMacroRule(MacroRule                       & rule,
                    const gsTensorBSplineBasis<2,T> & target,
                    const AtomicCell                & cell,
                    const gsVector<T>               & b,
                    index_t maxIter, T tol)
{
    const index_t nq = rule.nq;
    const index_t nz = 3 * nq;

    gsVector<T> F;
    gsSparseMatrix<T> J;

    // Timing accumulators (seconds)
    double t_basis = 0, t_assemble = 0, t_jtj = 0, t_solve = 0, t_apply = 0;
    index_t n_evalResJac = 0;

    gsStopwatch sw;
    sw.restart(); evalResJac<T>(rule, target, cell, b, F, J, &t_basis, &t_assemble); ++n_evalResJac;

    T res0 = F.norm();
    T res  = res0;

    gsInfo << "  Initial ||F|| = " << std::scientific << res0 << "\n";

    // Form J^T J and J^T F once (updated on every accepted step)
    sw.restart();
    gsSparseMatrix<T> JtJ = J.transpose() * J;   // (nz × nz), symmetric
    gsVector<T> JtF = J.transpose() * F;
    t_jtj += sw.stop();

    // Initial lambda: Marquardt heuristic scaled by largest diagonal entry
    T diagMax = JtJ.diagonal().maxCoeff();
    T lambda  = T(1e-3) * diagMax;
    if (lambda < T(1e-15)) lambda = T(1e-3);

    gsSparseMatrix<T> Id(nz, nz);
    Id.setIdentity();
    gsVector<T> z = extractState<T>(rule);

    for (index_t iter = 0; iter < maxIter; ++iter)
    {
        // Lambda floor: prevents collapse into numerical null space of J^T J
        // (J is underdetermined: nz > nt, so J^T J has null-space dimension nz-nt)
        const T lambdaFloor = diagMax * T(1e-12);
        lambda = std::max(lambda, lambdaFloor);

        // Solve (J^T J + lambda I) dz = -J^T F  via Cholesky LDLT
        sw.restart();
        typename gsSparseSolver<T>::CGDiagonal solver(JtJ + lambda * Id);
        gsVector<T> dz = solver.solve(-JtF);
        t_solve += sw.stop();

        // Trial state
        gsVector<T> z_trial = z + dz;

        sw.restart();
        MacroRule trial = rule;
        applyState<T>(trial, z_trial, cell);
        t_apply += sw.stop();

        gsVector<T> F_trial;
        gsSparseMatrix<T> J_trial;
        evalResJac<T>(trial, target, cell, b, F_trial, J_trial, &t_basis, &t_assemble);
        ++n_evalResJac;
        T res_trial = F_trial.norm();

        if (res_trial < res)
        {
            // Accept: update state and re-form J^T J
            rule = trial;
            z    = z_trial;
            for (index_t q = 0; q < nq; ++q)
            {
                z[3*q    ] = rule.xi[q];
                z[3*q + 1] = rule.eta[q];
            }
            F   = F_trial;
            J   = J_trial;
            res = res_trial;

            sw.restart();
            JtJ    = J.transpose() * J;
            JtF    = J.transpose() * F;
            diagMax = JtJ.diagonal().maxCoeff();
            t_jtj += sw.stop();

            lambda = std::max(lambdaFloor, lambda * T(0.1));

            if (res < tol)
            {
                gsInfo << "  Converged at iter " << iter + 1
                       << "  ||F|| = " << std::scientific << res << "\n";
                break;
            }
        }
        else
        {
            // Reject: increase lambda, J unchanged
            lambda *= T(10.0);
            if (lambda > T(1e30))
            {
                gsInfo << "  Diverged (lambda too large) at iter " << iter + 1 << "\n";
                break;
            }
        }

        if (iter < 5 || (iter + 1) % 20 == 0)
            gsInfo << "  iter " << std::setw(5) << iter + 1
                   << "  ||F|| = " << std::scientific << std::setprecision(4) << res
                   << "  lambda = " << std::scientific << lambda << "\n";
    }

    // ---- Timing report ----
    const double t_total = t_basis + t_assemble + t_jtj + t_solve + t_apply;
    gsInfo << "\n  === Timing breakdown (nq=" << nq
           << ", nt=" << b.size() << ", nz=" << nz << ") ===\n";
    gsInfo << std::fixed << std::setprecision(3);
    gsInfo << "  evalResJac calls : " << n_evalResJac << "\n";
    gsInfo << "  basis eval       : " << t_basis    << " s  ("
           << std::setprecision(1) << 100.*t_basis   /t_total << "%)\n";
    gsInfo << "  J assembly       : " << std::setprecision(3) << t_assemble
           << " s  (" << std::setprecision(1) << 100.*t_assemble/t_total << "%)\n";
    gsInfo << "  J^T J form       : " << std::setprecision(3) << t_jtj
           << " s  (" << std::setprecision(1) << 100.*t_jtj    /t_total << "%)\n";
    gsInfo << "  LDLT solve       : " << std::setprecision(3) << t_solve
           << " s  (" << std::setprecision(1) << 100.*t_solve  /t_total << "%)\n";
    gsInfo << "  applyState       : " << std::setprecision(3) << t_apply
           << " s  (" << std::setprecision(1) << 100.*t_apply  /t_total << "%)\n";
    gsInfo << "  accounted total  : " << std::setprecision(3) << t_total << " s\n\n";

    if (res < tol) return true;
    gsInfo << "  Max iter reached.  ||F|| = " << std::scientific << res << "\n";
    return false;
}

// ============================================================
// Validation against THB mass products
// ============================================================

/**
 * Validates the macro rule against all mass product moments B_i * B_j on D_r.
 *
 * Uses tensor-level-l evaluation (equivalent to THB at points inside the
 * atomic cell, where truncation corrections vanish) with the SAME_ELEMENT
 * optimisation for the reference Gauss integration.
 *
 * Comparison:
 *   reference: integral_{D_r} B_i B_j   (high-order Gauss)
 *   rule:      sum_q w_q B_i(x_q) B_j(x_q)
 */
template<class T>
void validateRule(const MacroRule          & rule,
                  const AtomicCell         & cell,
                  const gsHTensorBasis<2,T>& thb,
                  index_t overQuad = 4)
{
    // --- 1. Reference moments via high-order Gauss ---
    struct PairHash {
        size_t operator()(std::pair<index_t,index_t> p) const noexcept {
            return std::hash<size_t>()(
                (size_t(uint32_t(p.first)) << 32) | size_t(uint32_t(p.second)));
        }
    };
    std::unordered_map<std::pair<index_t,index_t>, index_t, PairHash> rowOfPair;
    rowOfPair.reserve(1024);
    std::vector<T> refMoments;

    gsOptionList opts;
    opts.addInt   ("quRule","Quadrature rule",gsQuadrature::GaussLegendre);
    opts.addReal  ("quA",   "Points",         T(1));
    opts.addInt   ("quB",   "Points",         (int)overQuad);
    opts.addSwitch("overInt","Over-int.",     false);
    typename gsQuadRule<T>::uPtr gauss =
        gsQuadrature::getPtr(thb.tensorLevel(cell.level), opts);

    const gsTensorBSplineBasis<2,T>& tensor =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(
            thb.tensorLevel(cell.level));

    gsFuncData<T> fd;
    fd.flags = NEED_VALUE | NEED_ACTIVE | SAME_ELEMENT;
    gsMatrix<T> nodes; gsVector<T> weights;

    for (const ElementBox & e : cell.elements)
    {
        gauss->mapTo(e.lower, e.upper, nodes, weights);
        tensor.compute(nodes, fd);

        const index_t na = fd.allActives().rows();
        const index_t nq = (index_t)weights.size();

        // Cache pair rows (determined once from actives at col 0)
        std::vector<index_t> rowCache;
        rowCache.reserve((na*(na+1))/2);
        for (index_t a = 0; a < na; ++a)
        {
            const index_t bi = fd.allActives()(a, 0);
            for (index_t c = a; c < na; ++c)
            {
                const index_t bj = fd.allActives()(c, 0);
                const index_t ii = std::min(bi, bj), jj = std::max(bi, bj);
                auto key = std::make_pair(ii, jj);
                auto it  = rowOfPair.find(key);
                index_t row;
                if (it == rowOfPair.end())
                {
                    row = (index_t)rowOfPair.size();
                    rowOfPair[key] = row;
                    refMoments.push_back(T(0));
                }
                else row = it->second;
                rowCache.push_back(row);
            }
        }

        for (index_t q = 0; q < nq; ++q)
        {
            const T w = weights[q];
            index_t k = 0;
            for (index_t a = 0; a < na; ++a)
                for (index_t c = a; c < na; ++c, ++k)
                    refMoments[rowCache[k]] +=
                        w * fd.allValues()(a,q) * fd.allValues()(c,q);
        }
    }

    // --- 2. Rule moments at the quadrature points ---
    const index_t numPairs = (index_t)rowOfPair.size();
    std::vector<T> ruleMoments(numPairs, T(0));

    // Evaluate tensor basis at all rule points (no SAME_ELEMENT: pts span elements)
    gsFuncData<T> fd2;
    fd2.flags = NEED_VALUE | NEED_ACTIVE;
    tensor.compute(rule.pts, fd2);

    const index_t na2 = fd2.allActives().rows();

    for (index_t q = 0; q < rule.nq; ++q)
    {
        const T wq = rule.wts[q];
        for (index_t a = 0; a < na2; ++a)
        {
            const index_t bi = fd2.allActives()(a, q);
            for (index_t c = a; c < na2; ++c)
            {
                const index_t bj = fd2.allActives()(c, q);
                const index_t ii = std::min(bi, bj), jj = std::max(bi, bj);
                auto it = rowOfPair.find({ii, jj});
                if (it != rowOfPair.end())
                    ruleMoments[it->second] +=
                        wq * fd2.allValues()(a, q) * fd2.allValues()(c, q);
            }
        }
    }

    // --- 3. Error report ---
    T maxAbsErr = T(0), sumSqErr = T(0);
    for (index_t r = 0; r < numPairs; ++r)
    {
        const T err = std::abs(ruleMoments[r] - refMoments[r]);
        maxAbsErr = std::max(maxAbsErr, err);
        sumSqErr += err * err;
    }

    const T sumW = rule.wts.sum();
    const T cellArea =
        (cell.upperParam[0] - cell.lowerParam[0]) *
        (cell.upperParam[1] - cell.lowerParam[1]);
    const T minW = rule.wts.minCoeff();
    const T maxW = rule.wts.maxCoeff();

    gsInfo << "\n  Validation (THB mass products B_i B_j on D_r):\n";
    gsInfo << "    Pairs:            " << numPairs << "\n";
    gsInfo << "    Max |err|:        " << std::scientific << maxAbsErr << "\n";
    gsInfo << "    RMS |err|:        " << std::scientific
           << std::sqrt(sumSqErr / T(numPairs)) << "\n";
    gsInfo << "    Sum of weights:   " << std::fixed << std::setprecision(8)
           << sumW << "  (cell area = " << cellArea << ")\n";
    gsInfo << "    Min/max weight:   " << minW << " / " << maxW << "\n";
    gsInfo << "    Negative weights: 0  (log-weight parametrisation)\n";
}

// ============================================================
// main
// ============================================================

int main(int argc, char * argv[])
{
    index_t cellId  = -1;   // -1 = auto-select smallest ~8-element cell
    index_t safety  = 8;    // extra points above ceil(rank/3)
    index_t maxIter = 200;
    real_t  tol     = 1e-10;

    gsCmdLine cmd("Nonlinear moment-fitting quadrature rule for one atomic cell.");
    cmd.addInt ("c", "cellId",  "Atomic cell to solve (-1 = auto)",    cellId);
    cmd.addInt ("s", "safety",  "Extra points above ceil(rank/3)",      safety);
    cmd.addInt ("i", "maxIter", "Max LM iterations",                    maxIter);
    cmd.addReal("t", "tol",     "Convergence tolerance on ||F||",       tol);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ------------------------------------------------------------------
    // Build THB-spline basis with L-shaped local refinement
    // (same as hierarchical_macro_mass_system.cpp)
    // ------------------------------------------------------------------
    gsKnotVector<> kv1(0, 1.0, 0, 3, 1);
    kv1.uniformRefine();
    kv1.uniformRefine();
    gsKnotVector<> kv2 = kv1;

    gsTensorBSplineBasis<2, real_t> tbsb(kv1, kv2);
    gsTHBSplineBasis<2, real_t>     thb(tbsb);

    {
        thb.uniformRefine();
        for (index_t lvl = 1; lvl <= 5; ++lvl)
        {
            const index_t N = 4 * (1 << lvl);
            std::vector<index_t> boxes;
            boxes.insert(boxes.end(), {lvl, 0, 0, N, 4});
            boxes.insert(boxes.end(), {lvl, 0, 0, 4, N});
            thb.refineElements(boxes);
        }
    }

    const index_t deg = thb.degree(0);
    gsInfo << "Basis: degree=" << deg
           << ", maxLevel=" << thb.maxLevel()
           << ", DOF=" << thb.size() << "\n\n";

    // ------------------------------------------------------------------
    // Collect atomic cells
    // ------------------------------------------------------------------
    const std::vector<AtomicCell> cells = collectAtomicCells<2, real_t>(thb);
    gsInfo << "Atomic cells: " << cells.size() << "\n";

    // ------------------------------------------------------------------
    // Select target cell
    // ------------------------------------------------------------------
    if (cellId < 0)
    {
        // Auto-select: prefer smallest cell (~8 elements)
        index_t bestScore = std::numeric_limits<index_t>::max();
        for (index_t r = 0; r < (index_t)cells.size(); ++r)
        {
            const index_t ne = (index_t)cells[r].elements.size();
            const index_t score = std::abs(ne - 8);
            if (score < bestScore) { bestScore = score; cellId = r; }
        }
        gsInfo << "Auto-selected cell " << cellId
               << " (elements=" << cells[cellId].elements.size() << ")\n\n";
    }
    else
    {
        GISMO_ENSURE(cellId < (index_t)cells.size(),
            "Cell index " << cellId << " out of range (0.." << cells.size()-1 << ")");
    }

    const AtomicCell & cell = cells[cellId];
    const index_t ne = (index_t)cell.elements.size();
    const index_t p  = deg;
    const index_t nx = cell.upperIndex[0] - cell.lowerIndex[0];
    const index_t ny = cell.upperIndex[1] - cell.lowerIndex[1];
    const index_t rank = ((p+1)*nx + p) * ((p+1)*ny + p);  // corrected bound
    const index_t nq   = (rank + 2) / 3 + safety;          // ceil(rank/3) + s
    const index_t nGauss = ne * (p+1) * (p+1);

    gsInfo << "Selected cell:  id=" << cellId
           << "  level=" << cell.level
           << "  elements=" << ne
           << "  [" << cell.lowerParam.transpose()
           << " , " << cell.upperParam.transpose() << "]\n";
    gsInfo << "  rank bound (corrected): " << rank << "\n";
    gsInfo << "  n_q = ceil(rank/3) + safety = "
           << (rank+2)/3 << " + " << safety << " = " << nq << "\n";
    gsInfo << "  Gauss pts (standard):   " << nGauss << "\n";
    gsInfo << "  Ratio n_q / Gauss:      "
           << std::fixed << std::setprecision(3)
           << (real_t)nq / (real_t)nGauss << "\n\n";

    if (nq > 500)
        gsInfo << "WARNING: n_q=" << nq << " is large. Consider a smaller cell "
               << "(-c N with small N) for a first test.\n\n";

    // ------------------------------------------------------------------
    // Build target basis and compute exact moments
    // ------------------------------------------------------------------
    gsInfo << "Building target basis S^{" << 2*p << "}_{" << (p-1)
           << "} on cell " << cellId << " ...\n";
    const gsTensorBSplineBasis<2,real_t> target =
        buildTargetBasis<real_t>(cell, thb);

    const index_t n_target = target.size();
    gsInfo << "  Target basis size: " << n_target
           << "  (= " << target.size(0) << " x " << target.size(1) << ")\n";
    GISMO_ENSURE(n_target == rank,
        "Target basis size " << n_target << " != rank bound " << rank);

    gsInfo << "Computing exact target moments ...\n";
    const gsVector<real_t> b = computeTargetMoments<real_t>(target);

    // Sanity: sum of all moments = cell area
    const real_t cellArea =
        (cell.upperParam[0] - cell.lowerParam[0]) *
        (cell.upperParam[1] - cell.lowerParam[1]);
    gsInfo << "  Sum of moments: " << std::fixed << std::setprecision(8)
           << b.sum() << "  (cell area = " << cellArea << ")\n\n";

    // ------------------------------------------------------------------
    // Initialise quadrature rule
    // ------------------------------------------------------------------
    gsInfo << "Initialising " << nq << " points across " << ne << " elements ...\n";
    MacroRule rule = initMacroRule<real_t>(cell, nq);

    // Compute initial residual (informative)
    {
        gsVector<real_t> F0; gsSparseMatrix<real_t> J0;
        evalResJac<real_t>(rule, target, cell, b, F0, J0);
        gsInfo << "  Initial ||F|| = " << std::scientific << F0.norm() << "\n\n";
    }

    // ------------------------------------------------------------------
    // Levenberg-Marquardt solve
    // ------------------------------------------------------------------
    gsInfo << "Solving LM (maxIter=" << maxIter
           << ", tol=" << std::scientific << tol << ") ...\n";
    const bool converged =
        solveMacroRule<real_t>(rule, target, cell, b, maxIter, real_t(tol));

    gsInfo << "\nSolution:\n";
    gsInfo << "  Converged:      " << (converged ? "YES" : "NO") << "\n";
    gsInfo << "  n_q:            " << nq << "\n";
    gsInfo << "  Gauss pts:      " << nGauss << "\n";
    gsInfo << "  Ratio:          " << std::fixed << std::setprecision(3)
           << (real_t)nq / (real_t)nGauss << "\n";
    gsInfo << "  Sum of weights: " << std::fixed << std::setprecision(8)
           << rule.wts.sum() << "  (cell area = " << cellArea << ")\n";
    gsInfo << "  Min/max weight: " << rule.wts.minCoeff()
           << " / " << rule.wts.maxCoeff() << "\n";

    // ------------------------------------------------------------------
    // Validate against THB mass products
    // ------------------------------------------------------------------
    gsInfo << "\nValidating against THB mass product moments ...\n";
    validateRule<real_t>(rule, cell, thb);

    return converged ? 0 : 1;
}
