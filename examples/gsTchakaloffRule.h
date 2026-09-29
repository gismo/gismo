/** @file gsTchakaloffRule.h

    @brief Tchakaloff-style compression of a discrete positive measure on a
    single box cell, by nonnegative least squares on a reduced orthogonal
    moment system.

    Given nodes x_i, weights w_i >= 0 (optionally unit normals n_i) that
    represent a measure mu on a box [lower, upper] subset R^3, this file
    produces a SUBSET of the input nodes, with strictly positive weights,
    that reproduces every moment of Q_2p := polynomials of degree <= 2p in
    each direction (K = (2p+1)^3 functions):

    - volume variant:   int q dmu           for every q in Q_2p        (K constraints)
    - boundary variant: int q dmu, int q n_x dmu, int q n_y dmu, int q n_z dmu
                         for every q in Q_2p                          (4K constraints)

    The boundary variant preserves int q n dS exactly because the compressed
    rule keeps a subset of the input nodes together with their exact
    (unmodified) normals.

    Theory (Tchakaloff / Caratheodory compression). Let M (N x m) be the
    matrix whose columns are the m constraint functionals evaluated at the N
    input nodes (m = K for volume, m = 4K for boundary, built from an
    orthonormal tensor-Legendre basis so M is well scaled regardless of the
    cell's position or size). The exact constraint on the compressed weights
    w is M^T w = M^T w_in. Since w_in itself is feasible (w_in >= 0 and
    M^T w_in = M^T w_in trivially), the nonnegative least-squares problem
        min_{w>=0} || M^T w - M^T w_in ||_2
    has optimal value 0, and any BASIC optimal solution of an NNLS problem
    with m rows has at most m nonzero entries (Caratheodory's theorem via the
    Lawson-Hanson active-set method: every passive-set matrix at every outer
    iteration has full column rank, so the passive set never exceeds
    rank(M) <= m). This is why compressing N >> K nodes down to at most
    rank(M) <= m nonzero weights is always possible in principle, and why
    Lawson-Hanson NNLS (not a generic LP or interior-point solver) is the
    right tool: it is built to exploit exactly this zero-residual,
    small-support structure.

    M itself is rank-deficient whenever the input measure is degenerate in a
    way the polynomial basis cannot resolve. An axis-aligned facet z = const
    with constant normal n = e_z is the sharpest example: every z-direction
    Legendre mode of order > 0 is evaluated at the same z on every node, so it
    is a scalar multiple of the order-0 mode there, and n_z == 1 makes the
    "n_z * V" column block bitwise identical to the plain "V" block --
    collapsing the boundary system's rank to (2p+1)^2 (the in-plane degrees of
    freedom only), far below both K = (2p+1)^3 and m = 4K; see this file's
    self-test cases F-axis (rank collapses exactly to (2p+1)^2) and F-tilted (a
    non-axis-aligned facet, rank still < 4K but not to a closed form).
    The reduction below never forms M^T M or attempts to invert anything
    rank-deficient: it factors M itself (N x m, tall) via a column-pivoted
    QR, keeps only the leading full-rank block of size r = rank(M), and
    hands Lawson-Hanson an r x N system with ORTHONORMAL rows (condition
    number exactly 1, independent of how M was scaled or how degenerate it
    is). See gsTetClip::detail::reduceChunk for the exact linear algebra and
    ::nnlsLawsonHanson for the active-set method.

    References:
    - C. L. Lawson and R. J. Hanson, *Solving Least Squares Problems*,
      Prentice-Hall, 1974, ch. 23 (algorithm NNLS).
    - A. Sommariva and M. Vianello, "Compression of multivariate discrete
      measures and applications", Numer. Funct. Anal. Optim. 36 (2015)
      1198-1223.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gismo.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

namespace gsTetClip {
using namespace gismo;

/// vals(k) = sqrt((2k+1)/(b-a)) * P_k(t),  t = (2x-a-b)/(b-a),  k = 0..deg
/// (orthonormal on [a,b]: \int_a^b vals_k vals_l dx = delta_kl)
inline void legendreOrthonormal(real_t x, real_t a, real_t b, index_t deg, gsVector<real_t>& vals)
{
    GISMO_ASSERT(deg >= 0, "deg must be >= 0");
    GISMO_ASSERT(b > a, "invalid interval [a,b]");

    const real_t t = (2 * x - a - b) / (b - a);
    vals.resize(deg + 1);

    // Three-term Legendre recurrence P_{j+1} = ((2j+1) t P_j - j P_{j-1}) / (j+1),
    // P_0 = 1, P_{-1} = 0 (same recurrence as gsGaussRule.hpp, but keeping
    // every intermediate degree instead of only the last one).
    real_t pnm1 = 0.0, pn = 1.0;
    vals(0) = pn;
    for (index_t j = 0; j < deg; ++j)
    {
        const real_t pnm2 = pnm1;
        pnm1 = pn;
        pn = (static_cast<real_t>(2 * j + 1) * t * pnm1 - static_cast<real_t>(j) * pnm2)
             / static_cast<real_t>(j + 1);
        vals(j + 1) = pn;
    }

    const real_t invLen = 1.0 / (b - a);
    for (index_t k = 0; k <= deg; ++k)
        vals(k) *= std::sqrt(static_cast<real_t>(2 * k + 1) * invLen);
}

// Internal helpers, not part of the public API listed in this file's brief.
// Kept in a nested namespace because gsTetClip is shared with
// gsTetMeshClip.h, which defines generic-sounding names of its own (e.g. a
// Vec3 cross3); unqualified helpers here would collide when a caller
// includes both headers. Opened early (rather than once, below
// ::legendreVandermonde) only for ::ensureAllocatable, which
// ::legendreVandermonde must call before its own resize().
namespace detail {

/// Guards an Eigen (rows x cols) allocation against the 32-bit index_t
/// product silently wrapping on this build (index_t is int, and gsEigen's
/// own Index is the same type via EIGEN_DEFAULT_DENSE_INDEX_TYPE, see
/// build/gsCore/gsConfig.h). Computes the product in 64-bit and fails
/// loudly before an overflowed size reaches Eigen's resize() -- a
/// GISMO_ENSURE, not a GISMO_ASSERT, because asserts are compiled out
/// under -DNDEBUG and this guard must fire in a Release build.
inline void ensureAllocatable(std::int64_t rows, std::int64_t cols, const char* where)
{
    const std::int64_t prod = rows * cols;
    GISMO_ENSURE(prod <= (std::int64_t)std::numeric_limits<index_t>::max(),
                 where << ": allocation " << rows << " x " << cols << " = " << prod
                       << " exceeds the 32-bit Eigen index bound "
                       << std::numeric_limits<index_t>::max());
}

} // namespace detail

/// V (N x K), K=(deg+1)^3, column index  k = kx + (deg+1)*(ky + (deg+1)*kz),
/// V(i,k) = Lx_kx(x_i) Ly_ky(y_i) Lz_kz(z_i), with L orthonormal on [lower_d, upper_d].
inline void legendreVandermonde(const gsMatrix<real_t>& nodes, const gsVector<real_t>& lower,
                                 const gsVector<real_t>& upper, index_t deg, gsMatrix<real_t>& V)
{
    GISMO_ASSERT(nodes.rows() == 3, "nodes must be 3 x N");
    const index_t N = nodes.cols();
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    detail::ensureAllocatable(N, K, "legendreVandermonde");
    V.resize(N, K);

    gsVector<real_t> vx, vy, vz;
    for (index_t i = 0; i < N; ++i)
    {
        legendreOrthonormal(nodes(0, i), lower(0), upper(0), deg, vx);
        legendreOrthonormal(nodes(1, i), lower(1), upper(1), deg, vy);
        legendreOrthonormal(nodes(2, i), lower(2), upper(2), deg, vz);
        for (index_t kz = 0; kz < n1; ++kz)
            for (index_t ky = 0; ky < n1; ++ky)
                for (index_t kx = 0; kx < n1; ++kx)
                    V(i, kx + n1 * (ky + n1 * kz)) = vx(kx) * vy(ky) * vz(kz);
    }
}

struct NnlsResult
{
    gsVector<real_t> w;        // length N, w >= 0
    index_t iterations = 0;    // outer (column-adding) iterations
    real_t  residual   = 0;    // ||A w - c||_2 / ||c||_2
    bool    converged  = false;// stopped by the residual or KKT test, not by the cap
};

// Internal helpers, not part of the public API listed in this file's brief.
// Kept in a nested namespace because gsTetClip is shared with
// gsTetMeshClip.h, which defines generic-sounding names of its own (e.g. a
// Vec3 cross3); unqualified helpers here would collide when a caller
// includes both headers.
namespace detail {

/// Dense column-major matrix type of the bundled Eigen (gsEigen), used
/// throughout this file's linear algebra (QR factorizations and their
/// Ref-bound in-place variants). Kept in gsTetClip::detail, not in the
/// shared gsTetClip namespace, since gsTetMeshClip.h also lives there.
typedef gsEigen::Matrix<real_t, gsEigen::Dynamic, gsEigen::Dynamic> DMat;

/// Why a stop test in ::nnlsLawsonHansonImpl ended the outer loop. Zero
/// (c == 0) returns before the loop. Residual: ||c - A w|| <= 1e-14*cn.
/// ZEmpty: every column is passive. Kkt: no column of Z has a reduced
/// gradient above 1e-15*cn (the Lawson-Hanson optimality certificate,
/// up to that threshold). Full: |P| == r. Exhausted: positive-gradient
/// candidates remain but every one was rejected by the rank or
/// feasibility test -- NOT an optimality certificate. Cap: maxIter hit.
/// All but Cap set NnlsResult::converged = true, so a caller holding
/// only NnlsResult cannot tell them apart; \c NnlsDiagnostics exists
/// to separate them, and to record, for a call whose returned residual is above
/// tolerance despite a "converged" stop, whether that is because the
/// incremental iterate was not yet the passive-set least-squares
/// minimiser (\c gPnormStop away from 0) or because the KKT scan itself,
/// restricted to Z, genuinely found no positive-gradient candidate above
/// the fixed \c 1e-15*cn threshold even though the residual has not
/// reached \c resTol.
enum class NnlsStop { Zero, Residual, ZEmpty, Kkt, Full, Exhausted, Cap };

/// Per-call NNLS diagnostics, filled only when a non-NULL sink is passed to
/// ::nnlsLawsonHansonImpl. Every field not documented inline is exactly the
/// quantity its name states, evaluated for the RETURNED w unless suffixed
/// \c Stop (evaluated instead at the outer-loop stop test that ended the
/// call, i.e. for the incremental iterate before the terminal re-solve).
/// \c condAP, the orthogonality/factor/gradient checks and the fresh
/// solve's rejected residual are diagnostic-only: computing them is gated
/// by `if (diag)` at every call site, so passing NULL costs nothing and
/// changes no floating-point operation the non-NULL path also executes.
struct NnlsDiagnostics
{
    NnlsStop stop = NnlsStop::Cap;
    index_t r = 0, N = 0, k = 0, iterations = 0, deletions = 0;
    real_t resInc = 0;        // ||c - A_P z_inc||/cn, the incremental iterate, before the terminal re-solve
    real_t resFresh = -1;     // ||c - A_P z_fresh||/cn, fresh ColPivHouseholderQR solve (computed even if rejected)
    bool   freshAccepted = false;
    index_t freshNonPos = 0;  // # entries of z_fresh that are <= 0
    real_t freshMin = 0;      // min entry of z_fresh
    real_t res = 0;           // the returned residual (== NnlsResult::residual)
    real_t orthErr = 0;       // max |(Q^T Q - I)_ij| at exit
    real_t factErr = 0;       // ||Q[:,0:k] R_k - A_P||_F / ||A_P||_F at exit
    real_t dErr = 0;          // ||d - Q^T c|| / cn at exit
    real_t gPnorm = 0;        // ||A_P^T (c - A w)|| / cn for the returned w
    real_t maxGZ = 0;         // max_{j in Z} g_j / cn for the returned w (<= 0 means KKT holds on Z)
    real_t gPnormStop = 0;    // same as gPnorm, but for the incremental iterate at the stop test that ended the loop
    real_t maxGZStop = 0;     // same as maxGZ, at that stop test (what the KKT scan actually saw)
    index_t rejNu = 0, rejDiag = 0, rejZ = 0; // rejections by cause, LAST outer iteration only
    index_t rejNuTotal = 0, rejDiagTotal = 0, rejZTotal = 0;
    real_t rho = 0, minRii = 0; // at exit: max colNorm over P, min |R_jj| over j<k
    real_t condAP = 0;        // sigma_max/sigma_min of A_P (JacobiSVD of the r x k matrix), at exit
    index_t refreshes = 0;    // reserved: factor refreshes performed; always 0 (no refresh mechanism exists)
};

/// One ::reduceChunk call's diagnostics (\c level 1 or 2, \c chunk its
/// index within that level, -1 for the level-2 union) plus that call's own
/// ::NnlsDiagnostics. Appended, one per call, to a ::ReduceDiagnostics by
/// ::reduceDriver in call order.
struct ChunkDiag { index_t level = 0, chunk = 0, n = 0, rank = 0; real_t wMin = 0, wMax = 0; NnlsDiagnostics nnls; };

/// One ::reduceDriver call's full diagnostic trail, one ::ChunkDiag per
/// ::reduceChunk call it made (level-1 chunks in call order, then the
/// level-2 union call when present).
struct ReduceDiagnostics { std::vector<ChunkDiag> calls; };

/// Length of the leading run of (non-increasing, by construction of a
/// column-pivoted QR) diagonal magnitudes |R_ii| that exceed
/// rankTol*|R_00|. This is the rank definition used throughout this file;
/// Eigen's own ColPivHouseholderQR::rank() is deliberately not used (it
/// scales its threshold by the maximal pivot found anywhere on the
/// diagonal, which coincides with |R_00| only when the diagonal is
/// monotone -- true here because of the column pivoting, but the leading
/// run is the form the linear algebra in ::reduceChunk actually needs: Q1
/// must be the first r columns of Q, so r must be the length of a
/// PREFIX, not a count of large entries anywhere on the diagonal).
///
/// Takes a Ref<const DMat> rather than a DMat so that it binds without a
/// copy to qr.matrixQR() when qr factors an in-place Ref<DMat> (the
/// overload used by ::reduceChunk) as well as to a plain DMat or gsMatrix
/// (every other call site), which would otherwise materialise a full
/// n x m temporary on every call.
inline index_t leadingRunRank(gsEigen::Ref<const DMat> R, real_t rankTol)
{
    const index_t m = std::min(R.rows(), R.cols());
    if (m == 0)
        return 0;
    const real_t r00 = std::abs(R(0, 0));
    if (r00 == 0.0)
        return 0;
    index_t r = 0;
    while (r < m && std::abs(R(r, r)) > rankTol * r00)
        ++r;
    return r;
}

/// Builds the constraint matrix M (N x m) whose columns are the moment
/// functionals: M = V for the volume variant (normals == NULL, m = K), or
/// M = [V | diag(n_x)V | diag(n_y)V | diag(n_z)V] for the boundary variant
/// (m = 4K); block d*K+k holds n_d(x_i)*V(i,k). Shared by ::reduceChunk and
/// by the self-test's Legendre moment check.
inline void buildConstraintMatrix(const gsMatrix<real_t>& nodes, const gsMatrix<real_t>* normals,
                                   const gsVector<real_t>& lower, const gsVector<real_t>& upper,
                                   index_t deg, gsMatrix<real_t>& M)
{
    gsMatrix<real_t> V;
    legendreVandermonde(nodes, lower, upper, deg, V);
    if (normals == NULL)
    {
        M = V;
        return;
    }

    GISMO_ASSERT(normals->rows() == 3 && normals->cols() == nodes.cols(), "normals must be 3 x N");
    const index_t n = nodes.cols();
    const index_t K = V.cols();
    ensureAllocatable(n, 4 * (std::int64_t)K, "buildConstraintMatrix");
    M.resize(n, 4 * K);
    M.leftCols(K) = V;
    for (index_t d = 0; d < 3; ++d)
    {
        const gsVector<real_t> nd = normals->row(d).transpose();
        M.middleCols((d + 1) * K, K) = V.array().colwise() * nd.array();
    }
}

/// Same constraint matrix as ::buildConstraintMatrix (bitwise identical
/// output: same three legendreOrthonormal calls per node, same
/// left-to-right product vx(kx)*vy(ky)*vz(kz), same column layout), built
/// row by row directly into \a M with no N x K Vandermonde temporary. Used
/// where M must be factored in place (::reduceChunk). ::buildConstraintMatrix
/// is a separate implementation (it does not delegate here), so the
/// self-test's and the benchmark's moment checks run on a code path
/// independent of the one this file's reduction actually factors.
inline void fillConstraintMatrix(const gsMatrix<real_t>& nodes, const gsMatrix<real_t>* normals,
                                  const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t deg,
                                  DMat& M)
{
    GISMO_ASSERT(nodes.rows() == 3, "nodes must be 3 x N");
    const index_t n = nodes.cols();
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    const index_t m = (normals != NULL) ? 4 * K : K;
    if (normals != NULL)
        GISMO_ASSERT(normals->rows() == 3 && normals->cols() == n, "normals must be 3 x N");
    ensureAllocatable(n, m, "fillConstraintMatrix");
    M.resize(n, m);

    gsVector<real_t> vx, vy, vz;
    for (index_t i = 0; i < n; ++i)
    {
        legendreOrthonormal(nodes(0, i), lower(0), upper(0), deg, vx);
        legendreOrthonormal(nodes(1, i), lower(1), upper(1), deg, vy);
        legendreOrthonormal(nodes(2, i), lower(2), upper(2), deg, vz);
        for (index_t kz = 0; kz < n1; ++kz)
            for (index_t ky = 0; ky < n1; ++ky)
                for (index_t kx = 0; kx < n1; ++kx)
                {
                    const index_t k = kx + n1 * (ky + n1 * kz);
                    const real_t v = vx(kx) * vy(ky) * vz(kz);
                    M(i, k) = v;
                    if (normals != NULL)
                        for (index_t d = 0; d < 3; ++d)
                            M(i, (d + 1) * K + k) = v * (*normals)(d, i);
                }
    }
}

/// Lawson-Hanson active-set NNLS: min ||A w - c||_2 s.t. w >= 0.  A is r x N.
/// Records per-call diagnostics into \a diag when non-NULL (see
/// ::NnlsDiagnostics); every diagnostic computation is gated by
/// `if (diag)`, so a NULL sink executes exactly the floating-point
/// operations the undiagnosed algorithm always did.
///
/// Algorithm NNLS (Lawson & Hanson 1974, ch. 23): maintains a passive set P
/// (candidate nonzero weights) and its complement Z, and at every outer
/// iteration moves the column of Z with the largest gradient component
/// A^T(c - A w) into P, resolves the unconstrained least-squares problem
/// restricted to P, and -- if that solution has any nonpositive entry --
/// backtracks w along the segment towards the new solution until the first
/// coordinate would go negative, evicting it (and any other coordinate that
/// rounds to <= 0) back to Z. The KKT stopping test (no column of Z has a
/// positive gradient) certifies a global optimum because the objective is
/// convex.
///
/// TODO: the KKT test above bounds only the DUAL
/// gradient (\c NnlsDiagnostics::maxGZStop), not the primal residual, and the two are related
/// quadratically rather than linearly. At a passive-set solution with a
/// feasible input w_in >= 0 (A w_in = c exactly -- e.g. a chunk's own
/// unreduced weights, which are always such a point of that chunk's own
/// reduced system by construction), the dual vector satisfies
/// g_Z^T w_in,Z = ||r||^2, where r is the current iterate's residual: the
/// gradient is quadratic in the residual, not linear in it. An ABSOLUTE
/// dual-vector tolerance tau can therefore stop Lawson-Hanson with ||r|| of
/// order sqrt(tau*||w_in||_1) -- a tolerance tight enough for one candidate
/// pool's column-norm scale is too loose for another. This has been
/// observed at about 4e-9 on small-mass level-1 chunks of real cut cells,
/// while the same nodes solved as one un-chunked pool converge to about
/// 1e-15. Candidate remedies, none applied here: a residual-relative or
/// feasibility-certified stopping test; a smaller candidate pool kept to a
/// single level (no chunking); or replacing NNLS with an LP (simplex) solve
/// or an explicit recombination scheme.
///
/// Passive-set solve: an incrementally updated QR of A_P, not a
/// refactorization (Lawson & Hanson 1974, ch. 24; Golub & Van Loan,
/// *Matrix Computations*, 4th ed., §6.5). Q is r x r orthogonal, R's
/// leading k x k block (k = |P|) is upper triangular, and d = Q^T c is kept
/// in step by construction: Q [R_k; 0] = A_P and every update is one of
///  - a Householder reflector, when a column is appended: it acts only on
///    Q's and d's rows/columns k..r-1, so it leaves Q[:,0:k], R_k and hence
///    A_P themselves untouched -- a rejected tentative append (rank test or
///    z(k) <= 0) therefore costs nothing to undo, it needs no rollback;
///  - a sequence of Givens rotations chasing the Hessenberg bulge left to
///    right, when a column is deleted.
/// Each update costs O(r^2) (a Householder append) or O((k-j)*(r+k)) (a
/// deletion at position j), instead of refactoring A_P from scratch at
/// O(r*k^2) per outer iteration. On every converged exit with a nonempty P,
/// one terminal ColPivHouseholderQR solve of A_P from scratch (O(r*k^2))
/// replaces the returned weights with its result whenever every entry is
/// strictly positive, which removes whatever drift the updates
/// accumulated; the updated factors only ever steered the active-set
/// decisions along the way. Cost per call: O(iters*(r*N + r^2)) plus
/// O(r*k^2) for the terminal solve; refactoring A_P from scratch on every
/// outer iteration would cost O(iters*(r*N + r*k^2)).
///
/// Stability. Every update is an orthogonal transformation applied in
/// floating point, and each is backward stable, so after s updates the
/// computed R is the exact R of a matrix within O(s*eps)*||A_P|| of A_P,
/// and Q is orthogonal to O(s*eps); s is the number of outer iterations
/// plus deletions for that solve, a small multiple of r, so the worst-case
/// growth is linear in the number of updates (typical growth is closer to
/// sqrt(s)). Errors accumulate additively here and are never amplified by
/// a condition number, which is why this updates A_P itself rather than
/// the normal equations: A has orthonormal rows, which bounds ||A_P|| <= 1
/// from above only, and sigma_min(A_P) is unbounded below except by the
/// anti-cycling guard (nnlsRankTol = 1e-11, i.e. kappa(A_P) up to about
/// 1e11 is accepted). An updated Cholesky factor of A_P^T A_P would square
/// this to kappa^2 ~ 1e22, far past 1/eps, and could lose every digit or
/// break down on problems this guard accepts; updating A_P directly keeps
/// the least-squares error scaling with kappa(A_P), not its square.
inline NnlsResult nnlsLawsonHansonImpl(const gsMatrix<real_t>& A, const gsVector<real_t>& c, index_t maxIter,
                                        NnlsDiagnostics* diag)
{
    typedef gsEigen::Matrix<real_t, gsEigen::Dynamic, gsEigen::Dynamic> DMat;
    const real_t nnlsRankTol = 1e-11;

    const index_t r = A.rows();
    const index_t N = A.cols();
    NnlsResult out;
    out.w = gsVector<real_t>::Zero(N);

    const real_t cn = c.norm();
    if (cn == 0.0)
    {
        out.converged = true;
        if (diag)
        {
            *diag = NnlsDiagnostics();
            diag->stop = NnlsStop::Zero;
            diag->r = r; diag->N = N;
        }
        return out;
    }

    gsVector<real_t> w = gsVector<real_t>::Zero(N);
    std::vector<bool> inP(N, false);
    std::vector<index_t> P;
    P.reserve(r);

    // Incremental-QR state for the passive-set least-squares problem (see
    // the doxygen above): Q [R.topLeftCorner(k,k); 0] = A_P, d = Q^T c,
    // AP's leading k columns hold A's passive-set columns verbatim (so the
    // top-of-loop residual is a contiguous r x k gemv instead of an r x N
    // one), and colNorm(j) = ||AP.col(j)||, used as the pivoted-QR
    // reference scale in the rank test below.
    DMat Q = DMat::Identity(r, r);
    DMat R(r, r);
    DMat AP(r, r);
    gsVector<real_t> d = c;
    gsVector<real_t> colNorm(r);
    gsVector<real_t> work(r);
    index_t k = 0;
    index_t iterations = 0;
    index_t deletions = 0;
    index_t rejNuTotal = 0, rejDiagTotal = 0, rejZTotal = 0;
    index_t rejNu = 0, rejDiag = 0, rejZ = 0; // this outer iteration's counts

    // Removes column at position j (0 <= j < k) from the k-column factor:
    // shifts columns j+1..k-1 of R/AP/colNorm and the matching P entries
    // one to the left (alias-safe, ascending), which leaves R upper
    // Hessenberg in columns j..k-2, then chases that Hessenberg bulge to
    // zero with k-1-j Givens rotations (Golub & Van Loan §6.5), applied
    // to R and d on the left (as G^T) and to Q on the right (as G), so
    // that Q[R_k;0] = A_P and d = Q^T c keep holding for the shrunken
    // k-1 column factor.
    auto deleteColumn = [&](index_t j)
    {
        const index_t kOld = k;
        for (index_t s = j + 1; s < kOld; ++s)
        {
            R.col(s - 1) = R.col(s);
            AP.col(s - 1) = AP.col(s);
            colNorm(s - 1) = colNorm(s);
        }
        P.erase(P.begin() + j);

        for (index_t i = j; i < kOld - 1; ++i)
        {
            gsEigen::JacobiRotation<real_t> G;
            G.makeGivens(R(i, i), R(i + 1, i), &R(i, i));
            R(i + 1, i) = 0.0;
            const index_t width = kOld - 2 - i;
            if (width > 0)
                R.middleCols(i + 1, width).applyOnTheLeft(i, i + 1, G.adjoint());
            Q.applyOnTheRight(i, i + 1, G);
            d.applyOnTheLeft(i, i + 1, G.adjoint());
        }
        --k;
        ++deletions;
    };

    // Fills *diag with the state at an outer-loop stop test (the "Stop"
    // suffixed fields use the INCREMENTAL iterate w, i.e. before the
    // terminal re-solve below runs; the un-suffixed exit fields are filled
    // separately, after that re-solve, once the returned w is final). A
    // no-op when diag is NULL, so every call site below costs nothing on
    // the default path.
    auto fillDiag = [&](NnlsStop stopReason)
    {
        if (!diag)
            return;
        diag->stop = stopReason;
        diag->r = r; diag->N = N; diag->k = k;
        diag->iterations = iterations;
        diag->deletions = deletions;
        diag->rejNu = rejNu; diag->rejDiag = rejDiag; diag->rejZ = rejZ;
        diag->rejNuTotal = rejNuTotal + rejNu;
        diag->rejDiagTotal = rejDiagTotal + rejDiag;
        diag->rejZTotal = rejZTotal + rejZ;

        gsVector<real_t> wPd(k);
        for (index_t j = 0; j < k; ++j)
            wPd(j) = w(P[j]);
        const gsVector<real_t> residStop = (k > 0) ? gsVector<real_t>(c - AP.leftCols(k) * wPd) : c;
        diag->resInc = residStop.norm() / cn;

        const gsVector<real_t> gStop = A.transpose() * residStop;
        real_t gPsum = 0.0;
        for (index_t j = 0; j < k; ++j)
            gPsum += gStop(P[j]) * gStop(P[j]);
        diag->gPnormStop = std::sqrt(gPsum) / cn;
        real_t maxGZs = 0.0;
        bool anyZ = false;
        for (index_t j = 0; j < N; ++j)
            if (!inP[j])
            {
                if (!anyZ || gStop(j) > maxGZs) maxGZs = gStop(j);
                anyZ = true;
            }
        diag->maxGZStop = anyZ ? (maxGZs / cn) : 0.0;

        real_t rhoNow = 0.0;
        for (index_t j = 0; j < k; ++j) rhoNow = std::max(rhoNow, colNorm(j));
        diag->rho = rhoNow;
        real_t minRiiNow = (k > 0) ? std::abs(R(0, 0)) : 0.0;
        for (index_t j = 0; j < k; ++j) minRiiNow = std::min(minRiiNow, std::abs(R(j, j)));
        diag->minRii = minRiiNow;
    };

    while (iterations < maxIter)
    {
        gsVector<real_t> wP(k);
        for (index_t j = 0; j < k; ++j)
            wP(j) = w(P[j]);
        const gsVector<real_t> resid = (k > 0) ? gsVector<real_t>(c - AP.leftCols(k) * wP) : c;
        const real_t residNorm = resid.norm();
        const gsVector<real_t> g = A.transpose() * resid;

        bool zEmpty = (static_cast<index_t>(P.size()) == N);
        real_t maxG = 1e-15 * cn;
        for (index_t j = 0; j < N; ++j)
            if (!inP[j] && g(j) > maxG)
                maxG = g(j);
        const bool anyCandidate = (maxG > 1e-15 * cn);

        if (residNorm <= 1e-14 * cn)
        {
            out.converged = true;
            fillDiag(NnlsStop::Residual);
            break;
        }
        if (zEmpty)
        {
            out.converged = true;
            fillDiag(NnlsStop::ZEmpty);
            break;
        }
        if (!anyCandidate)
        {
            out.converged = true;
            fillDiag(NnlsStop::Kkt);
            break;
        }
        if (k == r)
        {
            // A_{P union t} would need r+1 columns in R^r for any t: no
            // column can be appended. This avoids the O(|Z|*N) scan that
            // rejecting every remaining candidate in Z one by one would cost.
            out.converged = true;
            fillDiag(NnlsStop::Full);
            break;
        }

        // Reset this outer iteration's rejection tally, accumulating the
        // PREVIOUS iteration's into the running totals first.
        rejNuTotal += rejNu; rejDiagTotal += rejDiag; rejZTotal += rejZ;
        rejNu = 0; rejDiag = 0; rejZ = 0;

        // Grow the passive set P by one column: the candidate in Z with the
        // largest gradient component, retried with that column excluded
        // whenever the anti-cycling guard below rejects it.
        std::vector<bool> excludedNow(N, false);
        bool advanced = false;
        bool exhausted = false;
        gsVector<real_t> zP;

        while (!advanced && !exhausted)
        {
            index_t t = -1;
            real_t tg = 1e-15 * cn;
            for (index_t j = 0; j < N; ++j)
                if (!inP[j] && !excludedNow[j] && g(j) > tg)
                {
                    tg = g(j);
                    t = j;
                }
            if (t < 0)
            {
                exhausted = true;
                break;
            }

            ++iterations;

            const gsVector<real_t> at = A.col(t);
            const real_t atNorm = at.norm();
            gsVector<real_t> u = Q.transpose() * at;
            const real_t nu = u.tail(r - k).norm();
            real_t rho = atNorm;
            for (index_t j = 0; j < k; ++j)
                rho = std::max(rho, colNorm(j));

            // Rank test before touching Q, R or d: the tentative diagonal
            // would be |R_kk| = nu = dist(a_t, span(A_P)) (L&H ch. 23's own
            // independence test), and rho is what the pivoted |R_00| of
            // A_{P union t} would be (pivoting always selects the largest
            // column norm first). The second clause guards against rho
            // having grown since P's own columns were accepted.
            const bool nuOk = (nu > nnlsRankTol * rho);
            bool fullRank = nuOk;
            if (fullRank)
                for (index_t j = 0; j < k; ++j)
                    if (!(std::abs(R(j, j)) > nnlsRankTol * rho))
                    {
                        fullRank = false;
                        break;
                    }

            if (!fullRank)
            {
                // Q, R, d, AP and colNorm are all untouched: nothing to undo.
                // Distinguishes the two clauses above only for diagnostics
                // (see NnlsDiagnostics::rejNu/rejDiag): the candidate's own
                // independence test (nu) versus the passive-diagonal guard,
                // which tests P and so rejects the candidate for a reason
                // that has nothing to do with it.
                if (!nuOk) ++rejNu; else ++rejDiag;
                excludedNow[t] = true;
                continue;
            }

            real_t tau, beta;
            u.tail(r - k).makeHouseholderInPlace(tau, beta);
            if (r - k > 1)
            {
                Q.rightCols(r - k).applyHouseholderOnTheRight(u.tail(r - k - 1), tau, work.data());
                d.tail(r - k).applyHouseholderOnTheLeft(u.tail(r - k - 1), tau, work.data());
            }
            R.col(k).head(k) = u.head(k);
            R(k, k) = beta;

            const gsVector<real_t> z =
                R.topLeftCorner(k + 1, k + 1).template triangularView<gsEigen::Upper>().solve(d.head(k + 1));

            if (z(k) <= 0.0)
            {
                // The reflector above acted only on columns/rows k..r-1 of
                // Q and d: Q[:,0:k], R_k and AP.leftCols(k) are exactly as
                // before this tentative insertion (see the doxygen above),
                // so rejecting it costs nothing beyond not incrementing k.
                ++rejZ;
                excludedNow[t] = true;
                continue;
            }

            AP.col(k) = at;
            colNorm(k) = atNorm;
            P.push_back(t);
            inP[t] = true;
            ++k;

            zP = z;
            advanced = true;
        }

        if (exhausted)
        {
            out.converged = true;
            fillDiag(NnlsStop::Exhausted);
            break;
        }
        if (iterations >= maxIter)
        {
            fillDiag(NnlsStop::Cap);
            break;
        }

        // Restore feasibility: while the unconstrained solution z on the
        // current P has a nonpositive entry, move w towards z as far as the
        // first coordinate that would cross zero allows, evict it (and any
        // coordinate that rounds to <= 0) back to Z, and re-solve on the
        // shrunken P. Capped at |P|+1 passes as a safety net -- each pass
        // strictly shrinks P, so it always terminates well before the cap.
        const index_t capPasses = k + 1;
        for (index_t pass = 0; pass < capPasses; ++pass)
        {
            bool anyNeg = false;
            bool first = true;
            real_t alpha = 0.0;
            index_t qPos = -1;
            for (index_t j = 0; j < k; ++j)
            {
                if (zP(j) <= 0.0)
                {
                    anyNeg = true;
                    const real_t wj = w(P[j]);
                    const real_t denom = wj - zP(j);
                    const real_t a = (denom > 0.0) ? wj / denom : 0.0;
                    if (first || a < alpha)
                    {
                        alpha = a;
                        qPos = j;
                        first = false;
                    }
                }
            }
            if (!anyNeg)
                break;

            for (index_t j = 0; j < k; ++j)
            {
                const index_t idx = P[j];
                w(idx) = w(idx) + alpha * (zP(j) - w(idx));
            }
            w(P[qPos]) = 0.0;

            // Positions to delete, ascending; deleteColumn is applied in
            // descending order below so an earlier deletion never
            // invalidates a later one's position.
            std::vector<index_t> toRemove;
            for (index_t j = 0; j < k; ++j)
            {
                const index_t idx = P[j];
                if (j == qPos || w(idx) <= 0.0)
                {
                    w(idx) = 0.0;
                    inP[idx] = false;
                    toRemove.push_back(j);
                }
            }
            for (std::size_t ri = toRemove.size(); ri > 0; --ri)
                deleteColumn(toRemove[ri - 1]);

            if (k == 0)
            {
                zP.resize(0);
                break;
            }

            zP = R.topLeftCorner(k, k).template triangularView<gsEigen::Upper>().solve(d.head(k));
        }

        // Accept the restored solution as the new iterate and start the next
        // outer iteration with a clean exclusion set.
        for (index_t j = 0; j < k; ++j)
            w(P[j]) = zP(j);
    }

    // Terminal re-solve (only on a converged exit with a nonempty P; the
    // maxIter cap is left untouched): a from-scratch column-pivoted QR
    // solve of A_P removes whatever drift the incremental updates
    // accumulated. Guarded by strict positivity because the updated
    // iterate, not an infeasible from-scratch solve, must be returned when
    // rounding has pushed some coefficient to exactly the feasibility
    // boundary.
    if (out.converged && k > 0)
    {
        gsEigen::ColPivHouseholderQR<DMat> qrFinal(AP.leftCols(k));
        const gsVector<real_t> zFinal = qrFinal.solve(c);
        bool allPositive = true;
        for (index_t j = 0; j < k; ++j)
            if (!(zFinal(j) > 0.0))
            {
                allPositive = false;
                break;
            }
        if (diag)
        {
            index_t nonPos = 0;
            real_t minZ = zFinal.size() ? zFinal(0) : 0.0;
            for (index_t j = 0; j < k; ++j)
            {
                if (!(zFinal(j) > 0.0)) ++nonPos;
                minZ = std::min(minZ, zFinal(j));
            }
            diag->freshAccepted = allPositive;
            diag->freshNonPos = nonPos;
            diag->freshMin = minZ;
            const gsVector<real_t> residFresh = c - AP.leftCols(k) * zFinal;
            diag->resFresh = residFresh.norm() / cn;
        }
        if (allPositive)
            for (index_t j = 0; j < k; ++j)
                w(P[j]) = zFinal(j);
    }

    out.w = w;
    out.iterations = iterations;
    const gsVector<real_t> finalResid = c - A * w;
    out.residual = finalResid.norm() / cn;

    // Diagnostics at exit, for the RETURNED w/k/P (post terminal re-solve,
    // if one ran): factor-consistency and KKT quantities the caller cannot
    // otherwise see, since NnlsResult only carries the scalar residual.
    if (diag)
    {
        diag->res = out.residual;
        diag->iterations = iterations;
        diag->deletions = deletions;

        const gsVector<real_t> gFinal = A.transpose() * finalResid;
        real_t gPsum = 0.0;
        for (index_t j = 0; j < k; ++j)
            gPsum += gFinal(P[j]) * gFinal(P[j]);
        diag->gPnorm = std::sqrt(gPsum) / cn;
        real_t maxGZf = 0.0;
        bool anyZf = false;
        for (index_t j = 0; j < N; ++j)
            if (!inP[j])
            {
                if (!anyZf || gFinal(j) > maxGZf) maxGZf = gFinal(j);
                anyZf = true;
            }
        diag->maxGZ = anyZf ? (maxGZf / cn) : 0.0;

        // Orthogonality / factorization / rotated-rhs consistency of the
        // incrementally updated Q, R_k, d at exit, measured against A_P and c
        // directly: how far the Givens/Householder updates have drifted from
        // an exact QR of the final passive set.
        const DMat QtQ = Q.transpose() * Q;
        real_t orthErr = 0.0;
        for (index_t ii = 0; ii < r; ++ii)
            for (index_t jj = 0; jj < r; ++jj)
                orthErr = std::max(orthErr, std::abs(QtQ(ii, jj) - (ii == jj ? real_t(1) : real_t(0))));
        diag->orthErr = orthErr;

        if (k > 0)
        {
            // R's storage is never zeroed and deleteColumn's left-shift
            // leaves stale entries strictly below the diagonal (harmless
            // to the algorithm itself, which only ever reads R through an
            // Upper triangularView -- the two triangular solves above).
            // Reading the raw dense block here instead would report those
            // stale entries as factor error that was never there.
            const DMat APk = AP.leftCols(k);
            const DMat QRk = Q.leftCols(k) * DMat(R.topLeftCorner(k, k).template triangularView<gsEigen::Upper>());
            const real_t apNorm = APk.norm();
            diag->factErr = (apNorm > 0.0) ? (QRk - APk).norm() / apNorm : (QRk - APk).norm();
        }
        else
        {
            diag->factErr = 0.0;
        }

        const gsVector<real_t> dCheck = Q.transpose() * c;
        diag->dErr = (d - dCheck).norm() / cn;

        diag->condAP = 0.0;
        if (k > 0)
        {
            gsEigen::JacobiSVD<DMat> svd(AP.leftCols(k));
            const gsVector<real_t>& sv = svd.singularValues();
            const real_t smax = sv(0), smin = sv(sv.size() - 1);
            diag->condAP = (smin > 0.0) ? (smax / smin) : std::numeric_limits<real_t>::infinity();
        }
    }

    return out;
}

} // namespace detail

/// Lawson-Hanson active-set NNLS: min ||A w - c||_2 s.t. w >= 0. A is r x N.
/// See detail::nnlsLawsonHansonImpl for the algorithm, its stability
/// analysis and the theory references; this is the diagnostics-free public
/// entry point (forwards with a NULL sink).
inline NnlsResult nnlsLawsonHanson(const gsMatrix<real_t>& A, const gsVector<real_t>& c, index_t maxIter)
{
    return detail::nnlsLawsonHansonImpl(A, c, maxIter, NULL);
}

struct TchakaloffOptions
{
    /// Leading-run cutoff |R_ii| > rankTol*|R_00| (see ::leadingRunRank).
    ///
    /// The reduced system A = Q1^T has orthonormal rows by construction, so
    /// retaining more directions (a smaller rankTol) never worsens the
    /// conditioning the NNLS solve sees: A A^T = I_r regardless of r. The
    /// input weights satisfy every constraint exactly (A w_in = c), so the
    /// retained system stays feasible with w >= 0 at any r. Householder QR's
    /// own rounding floor on M is O(eps*sqrt(N)*||M||), on the order of 5e-15
    /// for the chunk sizes this file targets, so a 1e-13 cutoff keeps about a
    /// factor of 20 margin above that floor. A discarded direction of size
    /// rankTol*|R_00| contributes up to about rankTol*sqrt(N/r)*max|m| to the
    /// moment error (see the discarded-R22 bound in ::reduceChunk's doxygen
    /// below), so rankTol must sit well below the moment tolerance a caller
    /// checks against.
    real_t  rankTol   = 1e-13;
    real_t  resTol    = 1e-13; // required relative residual per level
    index_t chunkSize = 0;     // B; <=0 -> 4 * (#constraint columns): 4K volume, 16K boundary
    index_t maxIterFactor = 3; // NNLS outer cap = maxIterFactor * N
};

struct TchakaloffResult
{
    std::vector<index_t> indices;       // strictly increasing, into the INPUT columns
    gsVector<real_t>     weights;       // weights(j) > 0 belongs to input node indices[j]
    std::vector<real_t>  levelResidual; // one entry per level; level 1 = max over its chunks
    std::vector<index_t> levelRank;     // one entry per level; level 1 = MIN over its chunks
    index_t rank       = 0;             // r of the final level's system
    index_t iterations = 0;             // total NNLS outer iterations over all chunks and levels
    index_t levels     = 0;             // 1 or 2
    bool    ok         = false;         // all levelResidual < resTol, all NNLS converged, all weights > 0
};

namespace detail {

struct ChunkReduceResult
{
    std::vector<index_t> keptLocal; // indices into 0..n-1 of this chunk
    gsVector<real_t>     keptWeights;
    index_t rank       = 0;
    real_t  residual   = 0;
    index_t iterations = 0;
    bool    converged  = false;
};

/// Reduces one chunk (its own nodes/weights[/normals]) to a nonnegative
/// subset reproducing its moments.
///
/// Linear algebra (M is N x m, N = chunk size, m = K or 4K):
///  1. Factor M itself (not M^T): M*Pi = Q*R, Pi a column permutation, Q
///     N x N orthogonal, R N x m upper-trapezoidal with non-increasing
///     |R_ii| (guaranteed by column pivoting).
///  2. r = ::leadingRunRank(R, rankTol). Q1 = first r columns of Q; every
///     column of the discarded R22 block has norm <= |R_rr| <= rankTol*|R_00|.
///  3. Reduced system: A = Q1^T (r x N, orthonormal rows => condition number
///     exactly 1), c = Q1^T * w_in. NNLS then finds w >= 0 with A w = c to
///     within resTol, which by construction is a near-exact zero of the
///     original moment residual M^T w - M^T w_in (the discarded part scales
///     with rankTol, not with resTol).
///
/// Memory, in doubles (n = chunk size, m = constraint columns, r = rank):
/// M is filled directly into the storage that an in-place ColPivHouseholderQR
/// factors (Ref<DMat> binds to a DMat with no copy), and both are released
/// before A is formed, so M, its factorization and A are never all three
/// alive at once. The peaks are n*(m+r) while Q1 is being formed, 2*n*r <=
/// n*(m+r) while A is the transpose of Q1, and n*r + O(r^2+n) during the
/// NNLS solve -- n*(m+r) overall, half of what keeping M's factorization,
/// Q1 and A alive together would cost.
inline ChunkReduceResult reduceChunk(const gsMatrix<real_t>& nodesChunk, const gsVector<real_t>& wChunk,
                                      const gsMatrix<real_t>* normChunk, const gsVector<real_t>& lower,
                                      const gsVector<real_t>& upper, index_t deg, const TchakaloffOptions& opt,
                                      ReduceDiagnostics* diag = NULL)
{
    const index_t n = nodesChunk.cols();

    DMat M;
    fillConstraintMatrix(nodesChunk, normChunk, lower, upper, deg, M);

    index_t r;
    DMat Q1;
    {
        gsEigen::ColPivHouseholderQR<gsEigen::Ref<DMat> > qr(M);
        r = leadingRunRank(qr.matrixQR(), opt.rankTol);
        Q1 = qr.householderQ() * DMat::Identity(n, r); // moved, not copied
    } // qr's factorization is gone once this block ends
    M.resize(0, 0); // release the n x m storage

    const gsVector<real_t> c = Q1.transpose() * wChunk;
    gsMatrix<real_t> A = Q1.transpose(); // r x n
    Q1.resize(0, 0); // release the n x r storage

    // 64-bit product, clamped to index_t's range: at the sizes this file
    // targets maxIterFactor*n fits comfortably in 32 bits (maxIterFactor is
    // a small constant), but the multiplication itself must not be done in
    // index_t when n approaches 2^31/maxIterFactor.
    const std::int64_t maxIter64 = (std::int64_t)opt.maxIterFactor * (std::int64_t)n;
    const index_t maxIter =
        static_cast<index_t>(std::min<std::int64_t>(maxIter64, std::numeric_limits<index_t>::max()));
    NnlsDiagnostics nnlsDiag;
    const NnlsResult nr = nnlsLawsonHansonImpl(A, c, maxIter, diag ? &nnlsDiag : NULL);

    ChunkReduceResult cr;
    cr.rank = r;
    cr.residual = nr.residual;
    cr.iterations = nr.iterations;
    cr.converged = nr.converged;
    for (index_t j = 0; j < n; ++j)
        if (nr.w(j) > 0.0)
            cr.keptLocal.push_back(j);
    cr.keptWeights.resize(static_cast<index_t>(cr.keptLocal.size()));
    for (std::size_t j = 0; j < cr.keptLocal.size(); ++j)
        cr.keptWeights(static_cast<index_t>(j)) = nr.w(cr.keptLocal[j]);

    // level/chunk are left at their ChunkDiag default (0); the caller
    // (::reduceDriver, the only one with that context) stamps the entry
    // just appended once this call returns.
    if (diag)
    {
        ChunkDiag cd;
        cd.n = n;
        cd.rank = r;
        cd.wMin = (wChunk.size() > 0) ? wChunk.minCoeff() : 0.0;
        cd.wMax = (wChunk.size() > 0) ? wChunk.maxCoeff() : 0.0;
        cd.nnls = nnlsDiag;
        diag->calls.push_back(cd);
    }
    return cr;
}

/// Shared volume/boundary reduction driver: chunked level-1 reduction (M is
/// rebuilt per chunk, inside ::reduceChunk, and never formed for the full
/// input), followed by an optional level-2 reduction of the union of
/// level-1 survivors. There are never more than 2 levels.
///
/// Cost: per chunk, QR is O(n*m^2); forming Q1 is O(n*m*r); NNLS is
/// O(iters*(r*n + r^2)) plus O(r*k^2) for its terminal re-solve (see
/// ::nnlsLawsonHanson's doxygen).
///
/// Memory. Each ::reduceChunk call peaks at 8*n*(m+r) bytes (its own
/// doxygen derives this; n is that call's own node count). Level 1 has
/// n <= ceil(N/nc) <= B with nc = ceil(N/B) chunks. Level 2, when present,
/// reduces the union of level-1 survivors, whose size N_u = sum_i r_i is
/// bounded by nc*m (each chunk keeps at most its own rank r_i <= m
/// survivors). The peak over the whole call is therefore
/// 8*max(ceil(N/nc), N_u)*(m+r) bytes + O(N) -- **not** O(B*m + N): for
/// N large enough that N_u = nc*m dominates ceil(N/nc) (N >~ B^2/r), the
/// bound is O(N*r*(m+r)/B), i.e. O(N*m/2) for the volume variant
/// (B = 4m, r = m). Worked example, volume, p = 3, N = 2*10^5: B = 4*343 =
/// 1372, nc = ceil(N/B) = 146, N_u = 146*343 = 50078, so the peak is about
/// 8*50078*686 bytes ~= 262 MiB.
inline TchakaloffResult reduceDriver(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                                      const gsMatrix<real_t>* normals, const gsVector<real_t>& lower,
                                      const gsVector<real_t>& upper, index_t p, const TchakaloffOptions& opt,
                                      ReduceDiagnostics* diag = NULL)
{
    TchakaloffResult result;
    const bool isBoundary = (normals != NULL);
    const index_t N = nodes.cols();

    GISMO_ENSURE(nodes.rows() == 3, "nodes must be 3 x N");
    GISMO_ENSURE(weights.size() == N, "weights size mismatch");
    for (index_t i = 0; i < N; ++i)
        GISMO_ENSURE(weights(i) >= 0.0, "weights must be >= 0");
    if (isBoundary)
        GISMO_ENSURE(normals->rows() == 3 && normals->cols() == N, "normals must be 3 x N");
    GISMO_ENSURE(lower.size() == 3 && upper.size() == 3, "lower/upper must have size 3");
    for (index_t d = 0; d < 3; ++d)
        GISMO_ENSURE(lower(d) < upper(d), "lower must be < upper componentwise");
    GISMO_ENSURE(p >= 1, "p must be >= 1");

    const real_t sumW = weights.sum();
    if (N == 0 || sumW == 0.0)
    {
        result.ok = true;
        result.levels = 0;
        return result;
    }

    const index_t deg = 2 * p;
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    const index_t m = isBoundary ? 4 * K : K;
    const index_t B = (opt.chunkSize > 0) ? opt.chunkSize : 4 * m;

    // 64-bit chunk-boundary arithmetic: (ci*N) can exceed 2^31-1 once
    // nChunks-1 grows past about sqrt(2^31*B) (N ~ 1.04e6 at vol p=2), within
    // this file's real-data range. 64-bit so ci*N cannot wrap; the quotient
    // is the same formula as the 32-bit one wherever that did not overflow,
    // so every non-overflowing chunk boundary is bitwise identical (the
    // V-random/S-cap-chunk self-test cases depend on this).
    const std::int64_t nChunks64 = (N > B) ? ((std::int64_t)N + B - 1) / B : 1;
    const index_t nChunks = static_cast<index_t>(nChunks64); // <= N, fits index_t

    std::vector<index_t> keptIdx;
    std::vector<real_t>  keptWVec;
    real_t  levelRes1  = 0.0;
    index_t levelRank1 = -1;
    index_t totalIters = 0;
    bool allConverged  = true;

    for (index_t ci = 0; ci < nChunks; ++ci)
    {
        const std::int64_t lo64 = (std::int64_t)ci * N / nChunks64;
        const std::int64_t hi64 = (std::int64_t)(ci + 1) * N / nChunks64;
        const index_t lo = static_cast<index_t>(lo64); // <= N, fits index_t
        const index_t hi = static_cast<index_t>(hi64); // <= N, fits index_t
        const index_t n = hi - lo;

        const gsMatrix<real_t> nodesChunk = nodes.middleCols(lo, n);
        const gsVector<real_t> wChunk = weights.segment(lo, n);
        gsMatrix<real_t> normChunk;
        if (isBoundary)
            normChunk = normals->middleCols(lo, n);

        const ChunkReduceResult cr =
            reduceChunk(nodesChunk, wChunk, isBoundary ? &normChunk : NULL, lower, upper, deg, opt, diag);
        if (diag)
        {
            diag->calls.back().level = 1;
            diag->calls.back().chunk = ci;
        }

        levelRes1 = std::max(levelRes1, cr.residual);
        levelRank1 = (levelRank1 < 0) ? cr.rank : std::min(levelRank1, cr.rank);
        totalIters += cr.iterations;
        allConverged = allConverged && cr.converged;

        for (std::size_t j = 0; j < cr.keptLocal.size(); ++j)
        {
            keptIdx.push_back(lo + cr.keptLocal[j]);
            keptWVec.push_back(cr.keptWeights(static_cast<index_t>(j)));
        }
    }
    result.levelResidual.push_back(levelRes1);
    result.levelRank.push_back(levelRank1);

    if (nChunks > 1)
    {
        const index_t Nu = static_cast<index_t>(keptIdx.size());
        gsMatrix<real_t> unionNodes(3, Nu);
        gsVector<real_t> unionW(Nu);
        gsMatrix<real_t> unionNorm;
        if (isBoundary)
            unionNorm.resize(3, Nu);
        for (index_t j = 0; j < Nu; ++j)
        {
            unionNodes.col(j) = nodes.col(keptIdx[j]);
            unionW(j) = keptWVec[j];
            if (isBoundary)
                unionNorm.col(j) = normals->col(keptIdx[j]);
        }

        const ChunkReduceResult cr2 =
            reduceChunk(unionNodes, unionW, isBoundary ? &unionNorm : NULL, lower, upper, deg, opt, diag);
        if (diag)
        {
            diag->calls.back().level = 2;
            diag->calls.back().chunk = -1;
        }

        result.levelResidual.push_back(cr2.residual);
        result.levelRank.push_back(cr2.rank);
        totalIters += cr2.iterations;
        allConverged = allConverged && cr2.converged;

        std::vector<index_t> finalIdx;
        std::vector<real_t>  finalW;
        for (std::size_t j = 0; j < cr2.keptLocal.size(); ++j)
        {
            finalIdx.push_back(keptIdx[cr2.keptLocal[j]]);
            finalW.push_back(cr2.keptWeights(static_cast<index_t>(j)));
        }
        keptIdx.swap(finalIdx);
        keptWVec.swap(finalW);
        result.levels = 2;
        result.rank = cr2.rank;
    }
    else
    {
        result.levels = 1;
        result.rank = levelRank1;
    }

    std::vector<index_t> order(keptIdx.size());
    for (std::size_t j = 0; j < order.size(); ++j)
        order[j] = static_cast<index_t>(j);
    std::sort(order.begin(), order.end(),
              [&keptIdx](index_t a, index_t b) { return keptIdx[a] < keptIdx[b]; });

    result.indices.resize(keptIdx.size());
    result.weights.resize(static_cast<index_t>(keptIdx.size()));
    for (std::size_t j = 0; j < order.size(); ++j)
    {
        result.indices[j] = keptIdx[order[j]];
        result.weights(static_cast<index_t>(j)) = keptWVec[order[j]];
    }

    result.iterations = totalIters;

    bool allPositive = true;
    for (index_t j = 0; j < result.weights.size(); ++j)
        if (!(result.weights(j) > 0.0))
            allPositive = false;
    bool allResOk = true;
    for (std::size_t j = 0; j < result.levelResidual.size(); ++j)
        if (!(result.levelResidual[j] < opt.resTol))
            allResOk = false;

    result.ok = allResOk && allConverged && allPositive;
    return result;
}

} // namespace detail

/// Volume Tchakaloff compression: returns a subset of \a nodes, with
/// strictly positive weights, that reproduces int q dmu for every q in
/// Q_2p (K = (2p+1)^3 tensor-Legendre polynomials of degree <= 2p in each
/// direction), where mu = sum_i weights(i) delta_{nodes(i)} is the input
/// measure. The guarantee holds only when the returned \c ok is true, and
/// then up to the per-level NNLS residual opt.resTol and the rank
/// truncation at opt.rankTol (see TchakaloffOptions::rankTol).
///
/// \a nodes is 3 x N in the physical coordinates of the cell [\a lower,
/// \a upper]; \a weights has length N with every entry >= 0; \a lower and
/// \a upper have length 3 with lower < upper componentwise; \a p >= 1.
/// GISMO_ENSURE checks all of the above on entry. N == 0 or
/// weights.sum() == 0 is not an error: the result is empty with
/// ok = true, levels = 0. \a opt collects rankTol, resTol, chunkSize and
/// maxIterFactor (TchakaloffOptions).
///
/// In the returned TchakaloffResult, \c indices is strictly increasing
/// into the input columns 0..N-1, weights(j) > 0 belongs to input node
/// indices[j], and \c ok holds iff every level residual is < opt.resTol,
/// every NNLS solve converged, and every output weight is > 0 (see
/// detail::reduceDriver for the chunking and the at-most-two-level
/// composition that produce this result).
/// When \c ok is false the returned rule carries no moment guarantee: it
/// may be empty (e.g. an NNLS solve stopped by the maxIterFactor cap
/// before selecting any node) or only approximate; callers must check
/// \c ok.
///
/// Cost per reduction, with N the number of nodes in that reduction (a
/// level-1 chunk, or the level-2 union) and m the number of constraint
/// columns (m = K here): QR is O(N m^2); forming Q1 is O(N m r); NNLS is
/// O(iters*(r N + r^2)) plus O(r k^2) for its terminal re-solve (k = final
/// |P| <= r; see ::nnlsLawsonHanson's doxygen for the incremental-QR
/// update scheme), r = rank(M) <= m. See ::reduceDriver's doxygen for the
/// memory bound across the (at most two) chunked levels this cost is paid
/// over.
///
/// References: C. L. Lawson and R. J. Hanson, *Solving Least Squares
/// Problems*, Prentice-Hall, 1974, ch. 23-24 (algorithm NNLS and updating a
/// QR decomposition); G. H. Golub and C. F. Van Loan, *Matrix
/// Computations*, 4th ed., Johns Hopkins, 2013, §6.5 (updating matrix
/// factorizations); A. Sommariva and M. Vianello, "Compression of
/// multivariate discrete measures and applications", Numer. Funct. Anal.
/// Optim. 36 (2015) 1198-1223 (see this file's @file block for the full
/// theory).
inline TchakaloffResult tchakaloffCompress(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                                            const gsVector<real_t>& lower, const gsVector<real_t>& upper,
                                            index_t p, const TchakaloffOptions& opt = TchakaloffOptions())
{
    return detail::reduceDriver(nodes, weights, NULL, lower, upper, p, opt);
}

/// Boundary Tchakaloff compression: as tchakaloffCompress, and under the
/// same conditions (\c ok true, up to opt.resTol and opt.rankTol) also
/// reproduces int q n_x dmu, int q n_y dmu and int q n_z dmu for every
/// q in Q_2p (4K constraints total, m = 4K), which preserves int q n dS
/// exactly because the retained nodes keep their exact input \a normals
/// unmodified.
///
/// \a normals is the boundary-specific argument: 3 x N unit vectors, one
/// per node. See ::tchakaloffCompress for \a nodes, \a weights, \a lower,
/// \a upper, \a p, \a opt, the input checks, the N == 0 / zero-weight
/// case, the output contract and the cost bound (with m = 4K here).
inline TchakaloffResult tchakaloffCompressBoundary(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                                                    const gsMatrix<real_t>& normals, const gsVector<real_t>& lower,
                                                    const gsVector<real_t>& upper, index_t p,
                                                    const TchakaloffOptions& opt = TchakaloffOptions())
{
    return detail::reduceDriver(nodes, weights, &normals, lower, upper, p, opt);
}

// ---------------------------------------------------------------------
// Self-test.
// ---------------------------------------------------------------------

namespace detail {

inline index_t ipow(index_t base, index_t exp)
{
    index_t res = 1;
    for (index_t i = 0; i < exp; ++i)
        res *= base;
    return res;
}

/// Tensor Gauss quadrature of \a nPts points per direction, on each of
/// nSub^d equal sub-boxes of [lower, upper] (d = lower.size()). Used to
/// build redundant (but exact for the target polynomial degree) synthetic
/// measures for the self-test.
inline void tensorGaussSubdivided(const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t nSub,
                                   index_t nPts, gsMatrix<real_t>& nodes, gsVector<real_t>& weights)
{
    const index_t d = lower.size();
    gsVector<index_t> numNodes(d);
    numNodes.setConstant(nPts);
    gsGaussRule<real_t> g(numNodes);

    const index_t nSubCells = ipow(nSub, d);
    const index_t nPerCell = ipow(nPts, d);
    nodes.resize(d, nSubCells * nPerCell);
    weights.resize(nSubCells * nPerCell);

    index_t col = 0;
    for (index_t c = 0; c < nSubCells; ++c)
    {
        gsVector<real_t> subLo(d), subHi(d);
        index_t rem = c;
        for (index_t dd = 0; dd < d; ++dd)
        {
            const index_t ii = rem % nSub;
            rem /= nSub;
            const real_t hh = (upper(dd) - lower(dd)) / static_cast<real_t>(nSub);
            subLo(dd) = lower(dd) + static_cast<real_t>(ii) * hh;
            subHi(dd) = subLo(dd) + hh;
        }
        gsMatrix<real_t> subNodes;
        gsVector<real_t> subWeights;
        g.mapTo(subLo, subHi, subNodes, subWeights);
        nodes.middleCols(col, subNodes.cols()) = subNodes;
        weights.segment(col, subWeights.size()) = subWeights;
        col += subNodes.cols();
    }
}

inline void buildVGauss(const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t p,
                         gsMatrix<real_t>& nodes, gsVector<real_t>& weights)
{
    tensorGaussSubdivided(lower, upper, 2, 2 * p + 2, nodes, weights);
}

inline void buildVRandom(const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t K, real_t vol,
                          std::mt19937& rng, gsMatrix<real_t>& nodes, gsVector<real_t>& weights)
{
    const index_t N = 3 * (4 * K) + 17;
    nodes.resize(3, N);
    weights.resize(N);
    std::uniform_real_distribution<real_t> ux(lower(0), upper(0));
    std::uniform_real_distribution<real_t> uy(lower(1), upper(1));
    std::uniform_real_distribution<real_t> uz(lower(2), upper(2));
    std::uniform_real_distribution<real_t> uw(0.5, 1.5);
    for (index_t i = 0; i < N; ++i)
    {
        nodes(0, i) = ux(rng);
        nodes(1, i) = uy(rng);
        nodes(2, i) = uz(rng);
        weights(i) = uw(rng) * vol / static_cast<real_t>(N);
    }
}

inline void buildFAxis(const gsVector<real_t>& mid, real_t h, index_t p, gsMatrix<real_t>& nodes,
                        gsVector<real_t>& weights, gsMatrix<real_t>& normals)
{
    gsVector<real_t> lo2(2), hi2(2);
    lo2(0) = mid(0) - h / 4.0; lo2(1) = mid(1) - h / 4.0;
    hi2(0) = mid(0) + h / 4.0; hi2(1) = mid(1) + h / 4.0;

    gsMatrix<real_t> nodes2D;
    tensorGaussSubdivided(lo2, hi2, 2, 2 * p + 2, nodes2D, weights);

    const index_t N = nodes2D.cols();
    nodes.resize(3, N);
    normals.resize(3, N);
    for (index_t i = 0; i < N; ++i)
    {
        nodes(0, i) = nodes2D(0, i);
        nodes(1, i) = nodes2D(1, i);
        nodes(2, i) = mid(2);
        normals(0, i) = 0.0;
        normals(1, i) = 0.0;
        normals(2, i) = 1.0;
    }
}

/// Explicit 3-vector cross product. gsVector<real_t> is dynamic-size, and
/// Eigen's MatrixBase::cross() is only enabled for vectors whose SIZE IS
/// FIXED AT COMPILE TIME to 3 (a GISMO_ASSERT in debug, a silent
/// precondition in release) -- a runtime size of 3 does not satisfy it, so
/// it must not be used on gsVector here.
inline gsVector<real_t> cross3(const gsVector<real_t>& a, const gsVector<real_t>& b)
{
    gsVector<real_t> c(3);
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

inline void buildFTilted(const gsVector<real_t>& mid, real_t h, index_t p, gsMatrix<real_t>& nodes,
                          gsVector<real_t>& weights, gsMatrix<real_t>& normals)
{
    gsVector<real_t> n(3);
    n(0) = 1.0; n(1) = 2.0; n(2) = 3.0;
    n /= std::sqrt(14.0);
    gsVector<real_t> ez(3);
    ez(0) = 0.0; ez(1) = 0.0; ez(2) = 1.0;
    gsVector<real_t> t1 = cross3(n, ez);
    t1.normalize();
    gsVector<real_t> t2 = cross3(n, t1);

    gsVector<real_t> lo2(2), hi2(2);
    lo2(0) = -h / 4.0; lo2(1) = -h / 4.0;
    hi2(0) =  h / 4.0; hi2(1) =  h / 4.0;

    gsMatrix<real_t> st;
    tensorGaussSubdivided(lo2, hi2, 2, 2 * p + 2, st, weights);

    const index_t N = st.cols();
    nodes.resize(3, N);
    normals.resize(3, N);
    for (index_t i = 0; i < N; ++i)
    {
        nodes.col(i) = mid + st(0, i) * t1 + st(1, i) * t2;
        normals.col(i) = n;
    }
}

inline void buildSCap(const gsVector<real_t>& mid, index_t p, gsMatrix<real_t>& nodes,
                       gsVector<real_t>& weights, gsMatrix<real_t>& normals)
{
    const real_t R = 0.55;
    const real_t theta0 = 1.0, phi0 = 0.7;
    gsVector<real_t> d(3);
    d(0) = std::sin(theta0) * std::cos(phi0);
    d(1) = std::sin(theta0) * std::sin(phi0);
    d(2) = std::cos(theta0);
    const gsVector<real_t> cs = mid - R * d;

    gsVector<real_t> lo2(2), hi2(2);
    lo2(0) = theta0 - 0.15; lo2(1) = phi0 - 0.15;
    hi2(0) = theta0 + 0.15; hi2(1) = phi0 + 0.15;

    gsMatrix<real_t> tp;
    gsVector<real_t> wtp;
    tensorGaussSubdivided(lo2, hi2, 3, 2 * p + 4, tp, wtp);

    const index_t N = tp.cols();
    nodes.resize(3, N);
    weights.resize(N);
    normals.resize(3, N);
    for (index_t i = 0; i < N; ++i)
    {
        const real_t th = tp(0, i);
        const real_t ph = tp(1, i);
        const real_t sn = std::sin(th);
        gsVector<real_t> dir(3);
        dir(0) = sn * std::cos(ph);
        dir(1) = sn * std::sin(ph);
        dir(2) = std::cos(th);
        nodes.col(i) = cs + R * dir;
        normals.col(i) = dir;
        weights(i) = R * R * sn * wtp(i);
    }
}

/// Moments of the scaled monomials mu_abc(x) = ((x-mid)/h)^a ((y-mid)/h)^b
/// ((z-mid)/h)^c, 0 <= a,b,c <= deg, stacked in the same block layout as
/// ::buildConstraintMatrix (plain moments, then times n_x, n_y, n_z when
/// normals != NULL). Deliberately independent of ::legendreVandermonde /
/// ::legendreOrthonormal, so it can catch a broken recurrence or
/// normalisation that a check built from the same code could not.
inline void monomialMoments(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                             const gsMatrix<real_t>* normals, const gsVector<real_t>& mid, real_t h,
                             index_t deg, gsVector<real_t>& mom)
{
    const index_t n1 = deg + 1;
    const index_t K = n1 * n1 * n1;
    const index_t nBlocks = (normals != NULL) ? 4 : 1;
    mom = gsVector<real_t>::Zero(nBlocks * K);

    const index_t N = nodes.cols();
    gsVector<real_t> mx(n1), my(n1), mz(n1);
    for (index_t i = 0; i < N; ++i)
    {
        const real_t sx = (nodes(0, i) - mid(0)) / h;
        const real_t sy = (nodes(1, i) - mid(1)) / h;
        const real_t sz = (nodes(2, i) - mid(2)) / h;
        mx(0) = 1.0; my(0) = 1.0; mz(0) = 1.0;
        for (index_t k = 1; k < n1; ++k)
        {
            mx(k) = mx(k - 1) * sx;
            my(k) = my(k - 1) * sy;
            mz(k) = mz(k - 1) * sz;
        }
        for (index_t kz = 0; kz < n1; ++kz)
            for (index_t ky = 0; ky < n1; ++ky)
                for (index_t kx = 0; kx < n1; ++kx)
                {
                    const real_t basis = mx(kx) * my(ky) * mz(kz);
                    const index_t kidx = kx + n1 * (ky + n1 * kz);
                    mom(kidx) += weights(i) * basis;
                    if (normals != NULL)
                    {
                        mom(K + kidx)     += weights(i) * (*normals)(0, i) * basis;
                        mom(2 * K + kidx) += weights(i) * (*normals)(1, i) * basis;
                        mom(3 * K + kidx) += weights(i) * (*normals)(2, i) * basis;
                    }
                }
    }
}

/// Gathers the compressed nodes/weights[/normals] at \a res.indices and
/// runs checks 1-5 shared by every self-test case (result validity, weight
/// positivity, index ordering, the count-vs-rank bound, the Legendre and
/// independent-monomial moment checks). Every check below is evaluated
/// unconditionally (not short-circuited on an earlier failure) so the
/// accumulated \a failReason names every check that actually failed, except
/// that the moment-check gather is skipped -- not merely reported -- when
/// \a res.indices contains an out-of-range entry, since gathering through
/// such an index would read out of bounds. Case-specific analytic checks
/// are applied by the caller on top of this.
inline bool genericChecks(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                           const gsMatrix<real_t>* normals, const gsVector<real_t>& lower,
                           const gsVector<real_t>& upper, const gsVector<real_t>& mid, real_t h, index_t p,
                           const TchakaloffResult& res, real_t& maxMomErr, real_t& maxMonErr,
                           std::string& failReason)
{
    bool pass = true;
    std::string reason;
    // Accumulates every failing check into one semicolon-separated reason
    // string, instead of the first check's failReason hiding the rest.
    auto fail = [&pass, &reason](const char* msg)
    {
        pass = false;
        if (!reason.empty())
            reason += "; ";
        reason += msg;
    };

    const index_t N = nodes.cols();
    const index_t deg = 2 * p;

    if (!res.ok)
        fail("res.ok == false");
    for (std::size_t j = 0; j < res.levelResidual.size(); ++j)
        if (!(res.levelResidual[j] < 1e-13))
            fail("levelResidual >= 1e-13");
    for (index_t j = 0; j < res.weights.size(); ++j)
        if (!(res.weights(j) > 0.0))
            fail("non-positive output weight");

    bool indicesValid = true;
    for (std::size_t j = 0; j < res.indices.size(); ++j)
    {
        if (res.indices[j] < 0 || res.indices[j] >= N)
        {
            fail("index out of [0,N)");
            indicesValid = false;
        }
        if (j > 0 && !(res.indices[j] > res.indices[j - 1]))
            fail("indices not strictly increasing");
    }

    const index_t count = static_cast<index_t>(res.indices.size());
    if (!(count <= res.rank))
        fail("count > rank");

    maxMomErr = 0.0;
    maxMonErr = 0.0;
    if (indicesValid)
    {
        gsMatrix<real_t> nodesOut(3, count);
        gsMatrix<real_t> normOut;
        if (normals != NULL)
            normOut.resize(3, count);
        for (index_t j = 0; j < count; ++j)
        {
            nodesOut.col(j) = nodes.col(res.indices[j]);
            if (normals != NULL)
                normOut.col(j) = normals->col(res.indices[j]);
        }

        gsMatrix<real_t> Min, Mout;
        buildConstraintMatrix(nodes, normals, lower, upper, deg, Min);
        buildConstraintMatrix(nodesOut, normals != NULL ? &normOut : NULL, lower, upper, deg, Mout);
        const gsVector<real_t> momIn = Min.transpose() * weights;
        const gsVector<real_t> momOut = Mout.transpose() * res.weights;
        const real_t maxIn = momIn.array().abs().maxCoeff();
        const real_t maxDiff = (momOut - momIn).array().abs().maxCoeff();
        maxMomErr = (maxIn > 0.0) ? maxDiff / maxIn : maxDiff;
        if (!(maxMomErr <= 1e-12))
            fail("Legendre moment error > 1e-12");

        gsVector<real_t> monIn, monOut;
        monomialMoments(nodes, weights, normals, mid, h, deg, monIn);
        monomialMoments(nodesOut, res.weights, normals != NULL ? &normOut : NULL, mid, h, deg, monOut);
        const real_t maxMonIn = monIn.array().abs().maxCoeff();
        const real_t maxMonDiff = (monOut - monIn).array().abs().maxCoeff();
        maxMonErr = (maxMonIn > 0.0) ? maxMonDiff / maxMonIn : maxMonDiff;
        if (!(maxMonErr <= 1e-12))
            fail("monomial moment error > 1e-12");
    }
    else
    {
        fail("moment checks skipped: out-of-range index");
    }

    failReason = reason;
    return pass;
}

inline void printCaseLine(const std::string& id, index_t p, index_t N, const TchakaloffResult& res,
                           real_t maxMomErr, real_t maxMonErr, double timeSec, bool pass)
{
    gsInfo << (pass ? "PASS" : "FAIL") << "  " << id << "  p=" << p << "  N=" << N
           << "  levels=" << res.levels << "  r=" << res.rank << "  levelRank=";
    for (std::size_t i = 0; i < res.levelRank.size(); ++i)
    {
        if (i) gsInfo << ",";
        gsInfo << res.levelRank[i];
    }
    gsInfo << "  count=" << res.indices.size() << "  res=";
    for (std::size_t i = 0; i < res.levelResidual.size(); ++i)
    {
        if (i) gsInfo << ",";
        gsInfo << res.levelResidual[i];
    }
    gsInfo << "  maxMomErr=" << maxMomErr << "  maxMonErr=" << maxMonErr << "  iters=" << res.iterations
           << "  time=" << timeSec << "\n";
}

} // namespace detail

/// Runs every case below for p = pmin..pmax, prints one PASS/FAIL line per case; true iff all PASS.
inline bool tchakaloffSelfTest(index_t pmin = 1, index_t pmax = 3, bool verbose = true)
{
    std::mt19937 rng(12345);

    gsVector<real_t> lower(3), upper(3), mid(3);
    lower(0) = 0.10; lower(1) = -0.20; lower(2) = 0.30;
    const real_t h = 0.25;
    for (index_t d = 0; d < 3; ++d)
    {
        upper(d) = lower(d) + h;
        mid(d) = lower(d) + 0.5 * h;
    }
    const real_t vol = h * h * h;

    index_t nCases = 0, nPass = 0;

    for (index_t p = pmin; p <= pmax; ++p)
    {
        const index_t deg = 2 * p;
        const index_t K = (deg + 1) * (deg + 1) * (deg + 1);

        // --- V-gauss ---
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes;
            gsVector<real_t> weights;
            detail::buildVGauss(lower, upper, p, nodes, weights);
            TchakaloffOptions opt;
            opt.chunkSize = nodes.cols();
            const TchakaloffResult res = tchakaloffCompress(nodes, weights, lower, upper, p, opt);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            bool pass = detail::genericChecks(nodes, weights, NULL, lower, upper, mid, h, p, res, maxMomErr,
                                               maxMonErr, reason);
            auto fail = [&pass, &reason](const std::string& msg)
            {
                pass = false;
                if (!reason.empty()) reason += "; ";
                reason += msg;
            };

            if (res.rank != K)
                fail("r != K");

            const index_t count = static_cast<index_t>(res.indices.size());
            bool indicesValid = true;
            for (index_t j = 0; j < count; ++j)
                if (res.indices[j] < 0 || res.indices[j] >= nodes.cols())
                    indicesValid = false;
            if (indicesValid)
            {
                gsMatrix<real_t> nodesOut(3, count);
                for (index_t j = 0; j < count; ++j)
                    nodesOut.col(j) = nodes.col(res.indices[j]);
                gsMatrix<real_t> Vout;
                legendreVandermonde(nodesOut, lower, upper, deg, Vout);
                const gsVector<real_t> momOut = Vout.transpose() * res.weights;
                const real_t target0 = std::sqrt(vol);
                if (!(std::abs(momOut(0) - target0) <= 1e-12 * target0))
                    fail("m_0 != sqrt(vol)");
                bool higherModesZero = true;
                for (index_t k = 1; k < momOut.size(); ++k)
                    if (!(std::abs(momOut(k)) <= 1e-12 * target0))
                        higherModesZero = false;
                if (!higherModesZero)
                    fail("m_k != 0 for k>0");
            }
            else
            {
                fail("analytic moment check skipped: out-of-range index");
            }
            const real_t sumWOut = res.weights.sum();
            if (!(std::abs(sumWOut - vol) <= 1e-12 * vol))
                fail("sum(w_out) != vol");

            detail::printCaseLine("V-gauss", p, nodes.cols(), res, maxMomErr, maxMonErr, t, pass);
            if (!pass && verbose) gsInfo << "  FAILED CHECK: " << reason << "\n";
            ++nCases; nPass += pass ? 1 : 0;
        }

        // --- V-random ---
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes;
            gsVector<real_t> weights;
            detail::buildVRandom(lower, upper, K, vol, rng, nodes, weights);
            TchakaloffOptions opt;
            const TchakaloffResult res = tchakaloffCompress(nodes, weights, lower, upper, p, opt);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            bool pass = detail::genericChecks(nodes, weights, NULL, lower, upper, mid, h, p, res, maxMomErr,
                                               maxMonErr, reason);
            auto fail = [&pass, &reason](const std::string& msg)
            {
                pass = false;
                if (!reason.empty()) reason += "; ";
                reason += msg;
            };

            if (res.levels != 2)
                fail("levels != 2");
            for (std::size_t j = 0; j < res.levelRank.size(); ++j)
                if (res.levelRank[j] != K)
                    fail("levelRank != K");
            const real_t sumWIn = weights.sum();
            const real_t sumWOut = res.weights.sum();
            if (!(std::abs(sumWOut - sumWIn) <= 1e-12 * sumWIn))
                fail("sum(w_out) != sum(w_in)");

            detail::printCaseLine("V-random", p, nodes.cols(), res, maxMomErr, maxMonErr, t, pass);
            if (!pass && verbose) gsInfo << "  FAILED CHECK: " << reason << "\n";
            ++nCases; nPass += pass ? 1 : 0;
        }

        // --- F-axis ---
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes, normals;
            gsVector<real_t> weights;
            detail::buildFAxis(mid, h, p, nodes, weights, normals);
            TchakaloffOptions opt;
            const TchakaloffResult res = tchakaloffCompressBoundary(nodes, weights, normals, lower, upper, p, opt);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            bool pass = detail::genericChecks(nodes, weights, &normals, lower, upper, mid, h, p, res, maxMomErr,
                                               maxMonErr, reason);
            auto fail = [&pass, &reason](const std::string& msg)
            {
                pass = false;
                if (!reason.empty()) reason += "; ";
                reason += msg;
            };

            const index_t expectR = (2 * p + 1) * (2 * p + 1);
            if (res.rank != expectR)
                fail("r != (2p+1)^2");
            const real_t area = (h / 2.0) * (h / 2.0);
            const real_t sumWOut = res.weights.sum();
            if (!(std::abs(sumWOut - area) <= 1e-12 * area))
                fail("sum(w_out) != (h/2)^2");

            detail::printCaseLine("F-axis", p, nodes.cols(), res, maxMomErr, maxMonErr, t, pass);
            if (!pass && verbose) gsInfo << "  FAILED CHECK: " << reason << "\n";
            ++nCases; nPass += pass ? 1 : 0;
        }

        // --- F-tilted ---
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes, normals;
            gsVector<real_t> weights;
            detail::buildFTilted(mid, h, p, nodes, weights, normals);
            TchakaloffOptions opt;
            const TchakaloffResult res = tchakaloffCompressBoundary(nodes, weights, normals, lower, upper, p, opt);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            bool pass = detail::genericChecks(nodes, weights, &normals, lower, upper, mid, h, p, res, maxMomErr,
                                               maxMonErr, reason);
            auto fail = [&pass, &reason](const std::string& msg)
            {
                pass = false;
                if (!reason.empty()) reason += "; ";
                reason += msg;
            };

            if (!(res.rank < 4 * K))
                fail("r not < 4K");
            const real_t area = (h / 2.0) * (h / 2.0);
            const real_t sumWOut = res.weights.sum();
            if (!(std::abs(sumWOut - area) <= 1e-12 * area))
                fail("sum(w_out) != (h/2)^2");

            detail::printCaseLine("F-tilted", p, nodes.cols(), res, maxMomErr, maxMonErr, t, pass);
            if (!pass && verbose)
                gsInfo << "  FAILED CHECK: " << reason << " (r=" << res.rank << ")\n";
            ++nCases; nPass += pass ? 1 : 0;
        }

        // --- S-cap ---
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes, normals;
            gsVector<real_t> weights;
            detail::buildSCap(mid, p, nodes, weights, normals);

            bool insideCell = true;
            for (index_t i = 0; i < nodes.cols() && insideCell; ++i)
                for (index_t d = 0; d < 3; ++d)
                    if (nodes(d, i) < lower(d) || nodes(d, i) > upper(d))
                        insideCell = false;

            TchakaloffOptions opt;
            const TchakaloffResult res = tchakaloffCompressBoundary(nodes, weights, normals, lower, upper, p, opt);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            bool pass = detail::genericChecks(nodes, weights, &normals, lower, upper, mid, h, p, res, maxMomErr,
                                               maxMonErr, reason);
            auto fail = [&pass, &reason](const std::string& msg)
            {
                pass = false;
                if (!reason.empty()) reason += "; ";
                reason += msg;
            };
            if (!insideCell)
                fail("S-cap nodes fall outside the cell");

            const real_t sumWIn = weights.sum();
            const real_t sumWOut = res.weights.sum();
            if (!(std::abs(sumWOut - sumWIn) <= 1e-12 * sumWIn))
                fail("sum(w_out) != sum(w_in)");

            detail::printCaseLine("S-cap", p, nodes.cols(), res, maxMomErr, maxMonErr, t, pass);
            if (!pass && verbose) gsInfo << "  FAILED CHECK: " << reason << "\n";
            ++nCases; nPass += pass ? 1 : 0;
        }

        // --- S-cap-chunk ---
        {
            gsStopwatch sw;
            gsMatrix<real_t> nodes, normals;
            gsVector<real_t> weights;
            detail::buildSCap(mid, p, nodes, weights, normals);

            TchakaloffOptions opt;
            opt.chunkSize = (nodes.cols() + 2) / 3; // ceil(N/3)
            const TchakaloffResult res = tchakaloffCompressBoundary(nodes, weights, normals, lower, upper, p, opt);
            const double t = sw.stop();

            real_t maxMomErr = 0, maxMonErr = 0;
            std::string reason;
            bool pass = detail::genericChecks(nodes, weights, &normals, lower, upper, mid, h, p, res, maxMomErr,
                                               maxMonErr, reason);
            auto fail = [&pass, &reason](const std::string& msg)
            {
                pass = false;
                if (!reason.empty()) reason += "; ";
                reason += msg;
            };

            if (res.levels != 2)
                fail("levels != 2");
            const real_t sumWIn = weights.sum();
            const real_t sumWOut = res.weights.sum();
            if (!(std::abs(sumWOut - sumWIn) <= 1e-12 * sumWIn))
                fail("sum(w_out) != sum(w_in)");

            detail::printCaseLine("S-cap-chunk", p, nodes.cols(), res, maxMomErr, maxMonErr, t, pass);
            if (!pass && verbose) gsInfo << "  FAILED CHECK: " << reason << "\n";
            ++nCases; nPass += pass ? 1 : 0;
        }
    }

    if (nPass == nCases)
        gsInfo << "ALL " << nCases << " CASES PASS\n";
    else
        gsInfo << (nCases - nPass) << " OF " << nCases << " CASES FAILED\n";

    return nPass == nCases;
}

} // namespace gsTetClip
