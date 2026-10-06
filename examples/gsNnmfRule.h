/** @file gsNnmfRule.h

    @brief NNMF compression of the discrete positive measure of one cell
    (clip nodes and weights, optionally with unit normals) onto a small
    deterministic pool of its own nodes, by nonnegative least squares on a
    reduced orthogonal moment system with EXTERNAL target moments.

    Given nodes x_i, weights w_i >= 0 (optionally unit normals n_i) that
    represent a measure mu on a box [lower, upper] subset R^3, this file
    returns a SUBSET of the input nodes, with strictly positive weights, that
    reproduces every moment of Q_2p (polynomials of degree <= 2p in each
    direction, K = (2p+1)^3 functions) in the orthonormal tensor-Legendre
    basis of the box, to the relative maximum error NnmfOptions::momTol:

    - volume variant (::nnmfCompress):          int q dmu for q in Q_2p, m = K constraints
    - boundary variant (::nnmfCompressBoundary): int q dmu and int q n_d dmu,
      d = x,y,z, for q in Q_2p, m = 4K constraints

    The boundary variant keeps exact (unmodified) normals of the kept nodes.

    Algorithm.
    1. The exact moments mIn (m entries) of the full input are streamed with
       compensated summation, O(N m) time and O(m) memory; the N x m
       Vandermonde-type matrix of the full input is only formed in the
       full-set fallback round of step 5.
    2. A deterministic pool of input nodes is taken (no randomness, no
       hashing; the result is bitwise reproducible), in one of two modes:
       - NnmfOptions::blockSize == 0 (global stride): nPool = min(N,
         poolFactor*K) nodes at the indices floor(j*N/nPool), j = 0..nPool-1;
       - NnmfOptions::blockSize = B > 0 (per piece): the input is piece-major
         with B nodes per piece, N = nPieces*B, and piece b draws
         s_b = min(B, ceil(T*W_b/W)) quantile levels (none if W_b == 0),
         T = poolFactor*K initially, W_b the weight sum of the piece and W
         the total. The nodes are the quantiles of the piece's weights at
         golden-ratio-spread levels, see detail::nnmfBlockPoolIndices; a
         node hit by several levels is listed once, so the pool can be much
         smaller than T and need not grow when T grows. Growth multiplies
         T; T >= N, or a pool that happens to be the whole input, is the
         full set.
    3. Let M (nPool x m) be the constraint matrix of the pool and
       M Pi = Q R its column-pivoted QR, r = leadingRunRank(R) its numerical
       rank, Q1 the first r columns of Q, R11 the leading r x r block of R.
       The weights w of the pool must satisfy M^T w = mIn, i.e.
       R^T Q^T w = Pi^T mIn. The first r rows read
           R11^T (Q1^T w) = b,   b = (Pi^T mIn)(0..r-1),
       so c = R11^{-T} b is the target of the reduced system
           min_{w >= 0} || Q1^T w - c ||_2
       (orthonormal rows, condition number 1), solved by Lawson-Hanson NNLS.
       The remaining rows of the identity hold if and only if mIn lies in
       range(M^T), which is checked in step 4 and not assumed.
    4. Acceptance on the FULL moments: the kept nodes (w_j > 0) are streamed
       through the same moment routine and e = max|mOut - mIn| / max|mIn|
       must not exceed momTol. The NNLS convergence flag and reduced
       residual are not consulted.
    5. If the round is rejected and the pool is smaller than N, the pool is
       grown to min(N, growthFactor*nPool) (per piece: T to
       min(N, growthFactor*T)) and step 3 repeats. A round whose
       pool is the whole input (reached by growth or because N <=
       poolFactor*K) is the full-set fallback and the last round: its target
       is Q1^T w_in directly (feasible by construction, as in the Tchakaloff
       reduction of gsTchakaloffRule.h). If it is still rejected the result
       is returned with ok = false; nothing is thrown for non-acceptance.
    6. The result lists the kept input indices (ascending), their weights and
       the diagnostics of the final round. With NnmfOptions::recordRounds the
       result also carries one NnmfRound record per round (pool size, rank,
       NNLS iterations, convergence and residual, moment error, acceptance).

    Memory. The peak per round is about 8*nPool*(m+r) bytes (M, then Q1 and
    A = Q1^T), as in gsTetClip::detail::reduceChunk. The full-set fallback
    is the expensive round, with N in place of nPool.

    Thread safety. The file uses no OpenMP and no mutable static or global
    state; concurrent calls on distinct arguments share nothing, so it can be
    called from inside a caller's parallel loop, and the results do not
    depend on the thread count.

    References:
    - M. Garhuom and A. Duester, Comput. Mech. 70 (2022) 1059-1081.
    - R. Legrain, Comput. Math. Appl. 99 (2021) 270-291.
    - C. L. Lawson and R. J. Hanson, *Solving Least Squares Problems*,
      Prentice-Hall, 1974, ch. 23 (algorithm NNLS).

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gismo.h>

#include "gsTchakaloffRule.h"
#include "gsTetMeshClip.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace gsTetClip {
using namespace gismo;

/// Options of ::nnmfCompress and ::nnmfCompressBoundary.
struct NnmfOptions
{
    real_t  rankTol       = 1e-13; ///< leadingRunRank cutoff |R_ii| > rankTol*|R_00|
    real_t  momTol        = 1e-12; ///< acceptance: relative max moment error <= momTol
    index_t poolFactor    = 3;     ///< L: initial pool = min(N, L*K), K = (2p+1)^3
    index_t growthFactor  = 2;     ///< pool multiplier after a rejected round
    index_t maxIterFactor = 3;     ///< NNLS outer-iteration cap = maxIterFactor * nPool
    index_t blockSize     = 0;     ///< points per input piece (input is piece-major); 0 = global stride pool
    bool    recordRounds  = false; ///< fill NnmfResult::roundLog
};

/// Diagnostics of one round of ::nnmfCompress / ::nnmfCompressBoundary.
struct NnmfRound
{
    index_t poolSize       = 0;     ///< pool size of the round
    index_t rank           = 0;     ///< numerical rank of the pool's constraint matrix
    index_t kept           = 0;     ///< number of nodes with positive NNLS weight
    index_t nnlsIterations = 0;     ///< NNLS outer iterations of the round
    bool    nnlsConverged  = false; ///< NNLS stopped before its iteration cap
    real_t  nnlsResidual   = 0;     ///< relative residual reported by NNLS
    real_t  momErr         = 0;     ///< relative max moment error on the full moments
    bool    allPos         = false; ///< all kept weights are positive
    bool    accepted       = false; ///< the round met momTol
    bool    full           = false; ///< the pool was the whole input
};

/// Result of ::nnmfCompress and ::nnmfCompressBoundary.
struct NnmfResult
{
    std::vector<index_t> indices;  ///< kept input columns, strictly increasing
    gsVector<real_t>     weights;  ///< weights(j) > 0 belongs to indices[j]
    index_t rank          = 0;     ///< numerical rank r of the final round
    real_t  momErr        = 0;     ///< relative max moment error of the final round
    index_t poolSize      = 0;     ///< pool size of the final round
    index_t rounds        = 0;     ///< number of rounds run (0 iff the input is empty)
    bool    fallbackFull  = false; ///< the final round used the whole input as pool
    bool    ok            = false; ///< the final round met momTol
    index_t nnlsIterations = 0;    ///< total NNLS outer iterations over all rounds
    std::vector<NnmfRound> roundLog; ///< per-round records, only if NnmfOptions::recordRounds
};

namespace detail {

/// Moments of the discrete measure (nodes, weights[, normals]) in the tensor
/// Legendre basis orthonormal on [lower, upper], degree 2p per direction,
/// K = (2p+1)^3. normals == NULL: K entries, k = kx + n1*(ky + n1*kz).
/// Otherwise 4K entries, block b*K + k weighted by 1, n_x, n_y, n_z. One
/// compensated sum per entry; O(N K) time, O(K) memory, no N x K matrix.
inline void nnmfStreamMoments(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                              const gsMatrix<real_t>* normals, const gsVector<real_t>& lower,
                              const gsVector<real_t>& upper, index_t p, gsVector<real_t>& m)
{
    const index_t deg = 2 * p, n1 = deg + 1, K = n1 * n1 * n1;
    const index_t blocks = (NULL == normals) ? 1 : 4;
    std::vector<KahanSum> acc((size_t)(blocks * K));

    gsVector<real_t> vx(n1), vy(n1), vz(n1);
    for (index_t c = 0; c != nodes.cols(); ++c)
    {
        legendreOrthonormal(nodes(0, c), lower(0), upper(0), deg, vx);
        legendreOrthonormal(nodes(1, c), lower(1), upper(1), deg, vy);
        legendreOrthonormal(nodes(2, c), lower(2), upper(2), deg, vz);

        for (index_t kz = 0; kz != n1; ++kz)
        for (index_t ky = 0; ky != n1; ++ky)
        for (index_t kx = 0; kx != n1; ++kx)
        {
            const index_t k = kx + n1 * (ky + n1 * kz);
            const real_t base = weights[c] * vx[kx] * vy[ky] * vz[kz];
            acc[(size_t)k].add(base);
            if (NULL != normals)
            {
                acc[(size_t)(K + k)].add(base * (*normals)(0, c));
                acc[(size_t)(2 * K + k)].add(base * (*normals)(1, c));
                acc[(size_t)(3 * K + k)].add(base * (*normals)(2, c));
            }
        }
    }

    m.resize(blocks * K);
    for (index_t k = 0; k != blocks * K; ++k)
        m[k] = acc[(size_t)k].value();
}

/// max_k |mOut_k - mIn_k| / max_k |mIn_k|; 0 if both are all-zero, +inf if
/// only mIn is all-zero. A NaN entry yields a NaN or infinite result, never
/// a small one.
inline real_t nnmfMomentRelErr(const gsVector<real_t>& mIn, const gsVector<real_t>& mOut)
{
    GISMO_ENSURE(mIn.size() == mOut.size(), "nnmfMomentRelErr: size mismatch.");
    real_t maxIn = 0, maxDiff = 0;
    bool nan = false;
    for (index_t k = 0; k != mIn.size(); ++k)
    {
        const real_t d = math::abs(mIn[k] - mOut[k]);
        if (d != d || mIn[k] != mIn[k])
            nan = true;
        maxIn   = math::max(maxIn,   math::abs(mIn[k]));
        maxDiff = math::max(maxDiff, d);
    }
    if (nan)
        return std::numeric_limits<real_t>::quiet_NaN();
    if (0 == maxIn)
        return (0 == maxDiff) ? (real_t)0 : std::numeric_limits<real_t>::infinity();
    return maxDiff / maxIn;
}

/// Strided pool indices idx_j = floor(j*N/nPool), j = 0..nPool-1 (64-bit
/// product); strictly increasing because nPool <= N.
inline void nnmfPoolIndices(index_t N, index_t nPool, std::vector<index_t>& idx)
{
    idx.resize((size_t)nPool);
    for (index_t j = 0; j != nPool; ++j)
        idx[(size_t)j] = (index_t)(((std::int64_t)j * (std::int64_t)N) / (std::int64_t)nPool);
}

/// Per-piece pool indices for piece-major input of N = nPieces*blockSize
/// points, concatenated in piece order, hence strictly increasing. The pool
/// follows the measure: with W_b the weight sum of piece b and W the total,
/// piece b contributes s_b = min(blockSize, ceil(target*W_b/W)) nodes (none
/// if W_b == 0), chosen by the inverse CDF of the piece's weights at the
/// golden-ratio-spread levels W_b*frac(a_b + j/s_b), j = 0..s_b-1, with
/// a_b = frac(b*0.6180339887498949). Heavy nodes are hit repeatedly and
/// listed once, light nodes near the collapsed vertices are rarely drawn, so
/// the pool is a deterministic quantile sample of the input measure.
/// Plain double arithmetic, no RNG: bitwise independent of the thread count.
inline void nnmfBlockPoolIndices(const gsVector<real_t>& weights, index_t blockSize,
                                 std::int64_t target, std::vector<index_t>& idx)
{
    const index_t N = weights.size();
    const index_t nPieces = N / blockSize;
    std::vector<double> Wb((size_t)nPieces);
    double W = 0;
    for (index_t b = 0; b != nPieces; ++b)
    {
        double acc = 0;
        for (index_t o = 0; o != blockSize; ++o)
            acc += weights(b * blockSize + o);
        Wb[(size_t)b] = acc;
        W += acc;
    }
    idx.clear();
    std::vector<index_t> off;
    std::vector<double> lev;
    for (index_t b = 0; b != nPieces; ++b)
    {
        const double wb = Wb[(size_t)b];
        if (!(wb > 0))
            continue;
        const index_t base = b * blockSize;
        const index_t s = static_cast<index_t>(std::min<double>(
            (double)blockSize, std::max(1.0, std::ceil((double)target * wb / W))));
        lev.clear();
        double a = (double)b * 0.6180339887498949;
        a -= std::floor(a);
        for (index_t j = 0; j != s; ++j)
        {
            double u = a + (double)j / (double)s;
            u -= std::floor(u);
            lev.push_back(wb * u);
        }
        std::sort(lev.begin(), lev.end());
        off.clear();
        double cum = weights(base);
        index_t o = 0;
        for (std::vector<double>::const_iterator l = lev.begin(); l != lev.end(); ++l)
        {
            while (o + 1 < blockSize && cum <= *l)
                cum += weights(base + ++o);
            off.push_back(o);
        }
        std::sort(off.begin(), off.end());
        const std::vector<index_t>::iterator e = std::unique(off.begin(), off.end());
        for (std::vector<index_t>::const_iterator it = off.begin(); it != e; ++it)
            idx.push_back(base + *it);
    }
}

/// One round on the pool given by \a idx. When \a full is false the target
/// is the external moment vector \a mIn through the y-solve of the file
/// header (step 3); when true (idx is the whole input) the target is
/// Q1^T w_in. Fills the kept pool positions in \a kept (ascending, into idx),
/// their weights, and the rank.
inline void nnmfSolvePool(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                          const gsMatrix<real_t>* normals, const gsVector<real_t>& lower,
                          const gsVector<real_t>& upper, index_t deg, const std::vector<index_t>& idx,
                          bool full, const gsVector<real_t>& mIn, const NnmfOptions& opt,
                          std::vector<index_t>& kept, gsVector<real_t>& keptW, index_t& rank,
                          index_t& nnlsIter, bool& nnlsConverged, real_t& nnlsResidual)
{
    const index_t n = (index_t)idx.size();
    gsMatrix<real_t> pn(3, n), pnor;
    gsVector<real_t> pw(n);
    if (NULL != normals)
        pnor.resize(3, n);
    for (index_t j = 0; j != n; ++j)
    {
        pn.col(j) = nodes.col(idx[(size_t)j]);
        pw(j) = weights(idx[(size_t)j]);
        if (NULL != normals)
            pnor.col(j) = normals->col(idx[(size_t)j]);
    }

    DMat M;
    fillConstraintMatrix(pn, (NULL != normals) ? &pnor : NULL, lower, upper, deg, M);
    pn.resize(0, 0);
    pnor.resize(0, 0);

    index_t r;
    DMat Q1;
    gsVector<real_t> c;
    {
        gsEigen::ColPivHouseholderQR<gsEigen::Ref<DMat> > qr(M);
        r = leadingRunRank(qr.matrixQR(), opt.rankTol);
        if (r > 0)
        {
            Q1 = qr.householderQ() * DMat::Identity(n, r);
            if (full)
                c = Q1.transpose() * pw;
            else
            {
                gsVector<real_t> b(r);
                for (index_t i = 0; i < r; ++i)
                    b(i) = mIn(qr.colsPermutation().indices()(i));
                c = qr.matrixQR().topLeftCorner(r, r).transpose()
                        .triangularView<gsEigen::Lower>().solve(b);
            }
        }
    }
    M.resize(0, 0);
    rank = r;
    kept.clear();
    keptW.resize(0);
    nnlsConverged = false;
    nnlsResidual  = 0;
    if (r == 0)
        return;

    gsMatrix<real_t> A = Q1.transpose();
    Q1.resize(0, 0);
    const std::int64_t maxIter64 = (std::int64_t)opt.maxIterFactor * (std::int64_t)n;
    const index_t maxIter = static_cast<index_t>(
        std::min<std::int64_t>(maxIter64, std::numeric_limits<index_t>::max()));
    const NnlsResult nr = nnlsLawsonHanson(A, c, maxIter);
    nnlsIter += nr.iterations;
    nnlsConverged = nr.converged;
    nnlsResidual  = nr.residual;

    index_t cnt = 0;
    for (index_t j = 0; j != n; ++j)
        if (nr.w(j) > 0.0)
            ++cnt;
    keptW.resize(cnt);
    kept.reserve((size_t)cnt);
    cnt = 0;
    for (index_t j = 0; j != n; ++j)
        if (nr.w(j) > 0.0)
        {
            kept.push_back(j);
            keptW(cnt++) = nr.w(j);
        }
}

/// Shared implementation of ::nnmfCompress (normals == NULL) and
/// ::nnmfCompressBoundary.
inline void nnmfDriver(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                       const gsMatrix<real_t>* normals, const gsVector<real_t>& lower,
                       const gsVector<real_t>& upper, index_t p, const NnmfOptions& opt,
                       NnmfResult& res)
{
    res = NnmfResult();

    GISMO_ENSURE(nodes.rows() == 3, "nnmfCompress: nodes must be 3 x N.");
    GISMO_ENSURE(weights.size() == nodes.cols(), "nnmfCompress: weights size must equal the node count.");
    if (NULL != normals)
        GISMO_ENSURE(normals->rows() == 3 && normals->cols() == nodes.cols(),
                     "nnmfCompress: normals must be 3 x N.");
    GISMO_ENSURE(lower.size() == 3 && upper.size() == 3, "nnmfCompress: lower/upper must have size 3.");
    for (index_t d = 0; d < 3; ++d)
        GISMO_ENSURE(upper(d) > lower(d), "nnmfCompress: invalid box.");
    GISMO_ENSURE(p >= 0, "nnmfCompress: p must be >= 0.");
    GISMO_ENSURE(opt.poolFactor >= 1 && opt.growthFactor >= 2,
                 "nnmfCompress: poolFactor must be >= 1 and growthFactor >= 2.");

    const index_t N = nodes.cols();
    for (index_t i = 0; i < N; ++i)
        GISMO_ENSURE(std::isfinite(weights(i)) && weights(i) >= 0.0,
                     "nnmfCompress: weights must be finite and nonnegative.");

    if (N == 0)
    {
        res.ok = true;
        return;
    }

    const index_t deg = 2 * p, n1 = deg + 1, K = n1 * n1 * n1;

    gsVector<real_t> mIn;
    nnmfStreamMoments(nodes, weights, normals, lower, upper, p, mIn);

    const bool blocked = (opt.blockSize > 0);
    std::int64_t target = (std::int64_t)opt.poolFactor * (std::int64_t)K;
    index_t nPool = static_cast<index_t>(
        std::min<std::int64_t>((std::int64_t)N, (std::int64_t)opt.poolFactor * (std::int64_t)K));
    if (blocked)
    {
        GISMO_ENSURE(N % opt.blockSize == 0, "nnmfCompress: N must be a multiple of blockSize.");
    }

    std::vector<index_t> idx, kept;
    gsVector<real_t> keptW;
    for (;;)
    {
        bool full;
        if (blocked)
        {
            full = (target >= (std::int64_t)N);
            if (!full)
            {
                nnmfBlockPoolIndices(weights, opt.blockSize, target, idx);
                full = ((index_t)idx.size() == N);
            }
            if (full)
            {
                idx.resize((size_t)N);
                for (index_t i = 0; i != N; ++i)
                    idx[(size_t)i] = i;
            }
            nPool = (index_t)idx.size();
        }
        else
        {
            full = (nPool == N);
            nnmfPoolIndices(N, nPool, idx);
        }
        index_t rank = 0, roundIter = res.nnlsIterations;
        bool nnlsConv = false;
        real_t nnlsRes = 0;
        nnmfSolvePool(nodes, weights, normals, lower, upper, deg, idx, full, mIn, opt,
                      kept, keptW, rank, res.nnlsIterations, nnlsConv, nnlsRes);
        roundIter = res.nnlsIterations - roundIter;
        ++res.rounds;

        const index_t nk = (index_t)kept.size();
        std::vector<index_t> keptIdx((size_t)nk);
        gsMatrix<real_t> kn(3, nk), knor;
        if (NULL != normals)
            knor.resize(3, nk);
        gsVector<real_t> kw(nk);
        bool allPos = true;
        for (index_t j = 0; j != nk; ++j)
        {
            const index_t g = idx[(size_t)kept[(size_t)j]];
            keptIdx[(size_t)j] = g;
            kn.col(j) = nodes.col(g);
            if (NULL != normals)
                knor.col(j) = normals->col(g);
            kw(j) = keptW(j);
            allPos = allPos && (kw(j) > 0.0);
        }
        gsVector<real_t> mOut;
        nnmfStreamMoments(kn, kw, (NULL != normals) ? &knor : NULL, lower, upper, p, mOut);
        const real_t e = nnmfMomentRelErr(mIn, mOut);
        const bool accept = allPos && (e <= opt.momTol);

        if (opt.recordRounds)
        {
            NnmfRound rd;
            rd.poolSize = nPool; rd.rank = rank; rd.kept = nk;
            rd.nnlsIterations = roundIter; rd.nnlsConverged = nnlsConv;
            rd.nnlsResidual = nnlsRes; rd.momErr = e; rd.allPos = allPos;
            rd.accepted = accept; rd.full = full;
            res.roundLog.push_back(rd);
        }

        if (accept || full)
        {
            res.indices      = keptIdx;
            res.weights      = kw;
            res.rank         = rank;
            res.momErr       = e;
            res.poolSize     = nPool;
            res.fallbackFull = full;
            res.ok           = accept;
            return;
        }

        if (blocked)
            target = std::min<std::int64_t>((std::int64_t)N, (std::int64_t)opt.growthFactor * target);
        else
            nPool = static_cast<index_t>(
                std::min<std::int64_t>((std::int64_t)N, (std::int64_t)opt.growthFactor * (std::int64_t)nPool));
    }
}

} // namespace detail

/// Compresses the volume measure (\a nodes 3 x N, \a weights >= 0) of one
/// box cell [\a lower, \a upper] onto a subset of its nodes with strictly
/// positive weights that reproduces all K = (2p+1)^3 Q_2p moments in the
/// orthonormal tensor-Legendre basis (see the file header for the
/// algorithm). \a res.ok tells whether the relative max moment error is
/// <= opt.momTol; non-acceptance never throws. N == 0 gives ok = true,
/// an empty rule and rounds = 0. Input violations (shapes, box, p, option
/// ranges, non-finite or negative weights) are GISMO_ENSURE failures.
/// Peak memory per round about 8*nPool*(K+rank) bytes; the full-set
/// fallback uses N for nPool. The candidate pool is the global stride pool
/// (opt.blockSize == 0) or the per-piece pool (opt.blockSize > 0, N must be a
/// multiple of it, input piece-major); see the file header, step 2.
/// opt.recordRounds fills res.roundLog. Thread-safe, bitwise deterministic.
inline void nnmfCompress(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                         const gsVector<real_t>& lower, const gsVector<real_t>& upper, index_t p,
                         const NnmfOptions& opt, NnmfResult& res)
{
    detail::nnmfDriver(nodes, weights, NULL, lower, upper, p, opt, res);
}

/// Boundary variant of ::nnmfCompress: \a normals (3 x N, unit normals of
/// the input nodes) are carried by the kept nodes, and the 4K stacked moments
/// [1, n_x, n_y, n_z] x Q_2p are reproduced. Peak memory per round about
/// 8*nPool*(4K+rank) bytes. Same contract otherwise, including the two pool
/// modes and opt.recordRounds.
inline void nnmfCompressBoundary(const gsMatrix<real_t>& nodes, const gsVector<real_t>& weights,
                                 const gsMatrix<real_t>& normals, const gsVector<real_t>& lower,
                                 const gsVector<real_t>& upper, index_t p,
                                 const NnmfOptions& opt, NnmfResult& res)
{
    detail::nnmfDriver(nodes, weights, &normals, lower, upper, p, opt, res);
}

} // namespace gsTetClip
