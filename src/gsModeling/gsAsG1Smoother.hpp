/** @file gsAsG1Smoother.hpp

    @brief Specialized O(N) smoothers for AS-G1 multi-patch biharmonic multigrid solvers.

    Provides four A-symmetric smoother families sharing unified subspace/banded primitives:
      1. gsAsG1SkeletonSmootherOp  : Point Gauss-Seidel + Coupled Ring-1 skeleton Cholesky.
      2. gsAsG1BlockSmootherOp     : Multiplicative block Gauss-Seidel (alternating interior
                                     line strips + Coupled Ring-1 skeleton).
      3. gsAsG1FastDiagSmootherOp  : Patch-interior Fast Diagonalization + Ring-1 skeleton
                                     (with optional Point Gauss-Seidel hybrid sweep).
      4. gsAsG1ScmsSmootherOp      : O(p N_h) h- and p-robust Subspace-Corrected Mass Smoother
                                     (SCMS) + overlapping skeleton Schwarz smoother, supporting
                                     multiplicative, additive, analytical, and Sogn-Takacs
                                     hybrid (Point GS + SCMS + skeleton) modes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsSimplePreconditioners.h>

namespace gismo
{

// ============================================================================
// Section 1: Shared Internal Linear-Algebra, Banded, and Patch Primitives
// ============================================================================
namespace internal
{

/// Fast O(1)-on-nonzero check to detect identically zero iterates (e.g. coarse-grid pre-smooth).
template <class T>
inline bool isExactlyZero(const gsMatrix<T> &x)
{
    const T *ptr = x.data();
    const index_t n = x.size();
    for (index_t i = 0; i < n; ++i)
        if (ptr[i] != T(0))
            return false;
    return true;
}

/// Extract the principal submatrix A(dofs, dofs) from a symmetric ColMajor sparse matrix A.
/// Uses symmetry (column `row` of ColMajor A equals row `row` of A) to avoid storing a
/// duplicate RowMajor copy of A. `globToLoc` must be size >= A.rows() and initialized to -1;
/// it is restored to -1 on exit.
template <class T>
gsSparseMatrix<T> extractSymmetricSubmatrix(
    const gsSparseMatrix<T> &A,
    const std::vector<index_t> &dofs,
    std::vector<index_t> &globToLoc,
    index_t estNnzPerRow = 25)
{
    const index_t k = static_cast<index_t>(dofs.size());
    for (index_t i = 0; i < k; ++i)
        if (dofs[i] >= 0)
            globToLoc[dofs[i]] = i;

    gsSparseEntries<T> se;
    se.reserve(k * estNnzPerRow);
    for (index_t i = 0; i < k; ++i)
    {
        const index_t row = dofs[i];
        if (row < 0)
            continue;
        for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
        {
            const index_t j = globToLoc[it.row()];
            if (j >= 0)
                se.add(i, j, it.value());
        }
    }

    for (index_t idx : dofs)
        if (idx >= 0)
            globToLoc[idx] = -1;

    gsSparseMatrix<T> A_sub(k, k);
    A_sub.setFrom(se);
    return A_sub;
}

/// Exact sparse Cholesky correction block on an arbitrary subset of free DOFs.
/// Shared by Skeleton, Block, FastDiag, and SCMS smoothers.
template <class T>
class SubspaceCholeskyBlock
{
public:
    typedef gsEigen::SimplicialLDLT<gsSparseMatrix<T>> EigenLDLT;

    SubspaceCholeskyBlock() = default;

    bool init(const gsSparseMatrix<T> &A,
              const std::vector<index_t> &rawDofs,
              std::vector<index_t> &globToLoc,
              index_t estNnzPerRow = 24)
    {
        const index_t nFree = A.rows();
        m_dofs.clear();
        m_dofs.reserve(rawDofs.size());
        for (index_t idx : rawDofs)
        {
            if (idx >= 0 && idx < nFree && globToLoc[idx] < 0)
            {
                globToLoc[idx] = static_cast<index_t>(m_dofs.size());
                m_dofs.push_back(idx);
            }
        }
        const index_t k = static_cast<index_t>(m_dofs.size());
        if (k == 0)
            return false;

        gsSparseEntries<T> se;
        se.reserve(k * estNnzPerRow);
        for (index_t i = 0; i < k; ++i)
        {
            const index_t row = m_dofs[i];
            for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
            {
                const index_t j = globToLoc[it.row()];
                if (j >= 0)
                    se.add(i, j, it.value());
            }
        }
        for (index_t idx : m_dofs)
            globToLoc[idx] = -1;

        gsSparseMatrix<T> A_sub(k, k);
        A_sub.setFrom(se);
        m_solver.reset(new EigenLDLT());
        m_solver->compute(A_sub);
        return true;
    }

    bool empty() const { return !m_solver || m_dofs.empty(); }
    index_t size() const { return static_cast<index_t>(m_dofs.size()); }
    const std::vector<index_t> &dofs() const { return m_dofs; }

    /// Multiplicative subspace correction: x_S += A_SS^{-1} (rhs - A x)_S.
    void relax(const gsSparseMatrix<T> &A,
               const gsMatrix<T> &rhs,
               gsMatrix<T> &x,
               bool xIsZero = false) const
    {
        if (empty())
            return;

        const index_t k = static_cast<index_t>(m_dofs.size());
        const index_t ncols = rhs.cols();
        m_rSub.resize(k, ncols);

        if (ncols == 1)
        {
            const T *rhsPtr = rhs.data();
            T *xPtr = x.data();
            T *rPtr = m_rSub.data();
            if (xIsZero)
            {
                for (index_t i = 0; i < k; ++i)
                    rPtr[i] = rhsPtr[m_dofs[i]];
            }
            else
            {
                for (index_t i = 0; i < k; ++i)
                {
                    const index_t row = m_dofs[i];
                    T acc = rhsPtr[row];
                    for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
                        acc -= it.value() * xPtr[it.row()];
                    rPtr[i] = acc;
                }
            }
            m_eSub = m_solver->solve(m_rSub);
            const T *ePtr = m_eSub.data();
            for (index_t i = 0; i < k; ++i)
                xPtr[m_dofs[i]] += ePtr[i];
        }
        else
        {
            for (index_t i = 0; i < k; ++i)
            {
                const index_t row = m_dofs[i];
                m_rSub.row(i) = rhs.row(row);
                if (!xIsZero)
                    for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
                        m_rSub.row(i) -= it.value() * x.row(it.row());
            }
            m_eSub = m_solver->solve(m_rSub);
            for (index_t i = 0; i < k; ++i)
                x.row(m_dofs[i]) += m_eSub.row(i);
        }
    }

    /// Additive subspace correction from a precomputed global residual r: x_S += weight * A_SS^{-1} r_S.
    void applyAdditive(const gsMatrix<T> &r, gsMatrix<T> &x, T weight) const
    {
        if (empty())
            return;

        const index_t k = static_cast<index_t>(m_dofs.size());
        const index_t ncols = r.cols();
        m_rSub.resize(k, ncols);

        if (ncols == 1)
        {
            const T *rGlob = r.data();
            T *rPtr = m_rSub.data();
            for (index_t i = 0; i < k; ++i)
                rPtr[i] = rGlob[m_dofs[i]];

            m_eSub = m_solver->solve(m_rSub);
            T *xPtr = x.data();
            const T *ePtr = m_eSub.data();
            for (index_t i = 0; i < k; ++i)
                xPtr[m_dofs[i]] += weight * ePtr[i];
        }
        else
        {
            for (index_t i = 0; i < k; ++i)
                m_rSub.row(i) = r.row(m_dofs[i]);
            m_eSub = m_solver->solve(m_rSub);
            for (index_t i = 0; i < k; ++i)
                x.row(m_dofs[i]) += weight * m_eSub.row(i);
        }
    }

private:
    std::vector<index_t> m_dofs;
    std::unique_ptr<EigenLDLT> m_solver;
    mutable gsMatrix<T> m_rSub;
    mutable gsMatrix<T> m_eSub;
};

/// Narrow-banded SPD Cholesky solver (A = L L^T) stored in flat column-major band layout.
/// Avoids AMD ordering, indirect sparse indexing, heap allocations, and dense transposes:
///   - solveLeftInPlace(X)  computes X <- (L L^T)^{-1} X in O(bw * rows * cols).
///   - solveRightInPlace(X) computes X <- X (L L^T)^{-1} in O(bw * rows * cols)
///     directly along contiguous columns of column-major X (zero transposes).
template <class T>
class BandedSPDCholesky
{
public:
    BandedSPDCholesky() : m_n(0), m_bw(0) {}

    void compute(const gsSparseMatrix<T> &A)
    {
        m_n = A.rows();
        if (m_n == 0)
        {
            m_bw = 0;
            m_L.resize(0, 0);
            m_invDiag.resize(0);
            return;
        }

        m_bw = 0;
        for (index_t col = 0; col < m_n; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(A, col); it; ++it)
                if (it.row() >= col)
                    m_bw = std::max(m_bw, it.row() - col);

        const index_t stride = m_bw + 1;
        m_L.setZero(stride, m_n);
        m_invDiag.resize(m_n);

        for (index_t col = 0; col < m_n; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(A, col); it; ++it)
                if (it.row() >= col)
                    m_L(it.row() - col, col) = it.value();

        // Right-looking banded Cholesky (LAPACK dpbtrf)
        T *Ldata = m_L.data();
        for (index_t j = 0; j < m_n; ++j)
        {
            T *colJ = Ldata + j * stride;
            const T diag = (colJ[0] > T(1e-30)) ? colJ[0] : T(1e-30);
            const T ljj = std::sqrt(diag);
            const T invLjj = T(1) / ljj;
            colJ[0] = ljj;
            m_invDiag(j) = invLjj;

            const index_t len = std::min(m_bw, m_n - 1 - j);
            for (index_t k = 1; k <= len; ++k)
                colJ[k] *= invLjj;

            for (index_t c = 1; c <= len; ++c)
            {
                const T ljc = colJ[c];
                if (ljc == T(0))
                    continue;
                T *colTarget = Ldata + (j + c) * stride;
                const index_t rem = len - c;
                for (index_t d = 0; d <= rem; ++d)
                    colTarget[d] -= ljc * colJ[c + d];
            }
        }
    }

    index_t size() const { return m_n; }
    index_t bandwidth() const { return m_bw; }

    /// In-place left solve: X <- (L L^T)^{-1} X for column-major X (m_n x ncols).
    void solveLeftInPlace(gsMatrix<T> &X) const
    {
        if (m_n == 0 || X.cols() == 0)
            return;

        const index_t stride = m_bw + 1;
        const T *Ldata = m_L.data();
        const T *invD = m_invDiag.data();
        const index_t ncols = X.cols();

        for (index_t col = 0; col < ncols; ++col)
        {
            T *x = X.data() + col * m_n;
            // Forward solve L y = b
            for (index_t j = 0; j < m_n; ++j)
            {
                const T yj = x[j] * invD[j];
                x[j] = yj;
                const index_t len = std::min(m_bw, m_n - 1 - j);
                const T *lcol = Ldata + j * stride;
                for (index_t k = 1; k <= len; ++k)
                    x[j + k] -= lcol[k] * yj;
            }
            // Backward solve L^T x = y
            for (index_t j = m_n - 1; j >= 0; --j)
            {
                T sum = x[j];
                const index_t len = std::min(m_bw, m_n - 1 - j);
                const T *lcol = Ldata + j * stride;
                for (index_t k = 1; k <= len; ++k)
                    sum -= lcol[k] * x[j + k];
                x[j] = sum * invD[j];
            }
        }
    }

    /// In-place left solve for a vector.
    void solveLeftInPlace(gsVector<T> &x) const
    {
        if (m_n == 0)
            return;

        const index_t stride = m_bw + 1;
        const T *Ldata = m_L.data();
        const T *invD = m_invDiag.data();
        T *xp = x.data();

        for (index_t j = 0; j < m_n; ++j)
        {
            const T yj = xp[j] * invD[j];
            xp[j] = yj;
            const index_t len = std::min(m_bw, m_n - 1 - j);
            const T *lcol = Ldata + j * stride;
            for (index_t k = 1; k <= len; ++k)
                xp[j + k] -= lcol[k] * yj;
        }
        for (index_t j = m_n - 1; j >= 0; --j)
        {
            T sum = xp[j];
            const index_t len = std::min(m_bw, m_n - 1 - j);
            const T *lcol = Ldata + j * stride;
            for (index_t k = 1; k <= len; ++k)
                sum -= lcol[k] * xp[j + k];
            xp[j] = sum * invD[j];
        }
    }

    /// In-place right solve: X <- X (L L^T)^{-1} = (X L^{-T}) L^{-1} for column-major X (nrows x m_n).
    /// Operates purely via contiguous column AXPYs (zero transposes, zero allocations).
    void solveRightInPlace(gsMatrix<T> &X) const
    {
        if (m_n == 0 || X.rows() == 0)
            return;

        const index_t stride = m_bw + 1;
        const T *Ldata = m_L.data();
        const T *invD = m_invDiag.data();

        // Step 1: Y L^T = B  (forward sweep over columns i = 0 .. m_n - 1)
        for (index_t i = 0; i < m_n; ++i)
        {
            X.col(i) *= invD[i];
            const index_t len = std::min(m_bw, m_n - 1 - i);
            const T *lcol = Ldata + i * stride;
            for (index_t k = 1; k <= len; ++k)
                X.col(i + k).noalias() -= lcol[k] * X.col(i);
        }

        // Step 2: X L = Y  (backward sweep over columns j = m_n - 1 .. 0)
        for (index_t j = m_n - 1; j >= 0; --j)
        {
            const index_t len = std::min(m_bw, m_n - 1 - j);
            const T *lcol = Ldata + j * stride;
            for (index_t k = 1; k <= len; ++k)
                X.col(j).noalias() -= lcol[k] * X.col(j + k);
            X.col(j) *= invD[j];
        }
    }

private:
    index_t m_n;
    index_t m_bw;
    gsMatrix<T> m_L;
    gsVector<T> m_invDiag;
};

/// Structured representation of the localized 1D nullspace basis Q_0 (m x nD0) of S^{D,0}.
/// Eliminates sparse-dense matrix products during SCMS smoothing:
///   - Identity  (p = 3, r_ev = 0)  : Q_0 = I_m (O(1) pass-through).
///   - BlockDiag (m >= 2 * (p - 1)) : diag(V_L, I_{nMid}, V_R) with tiny dense corner GEMMs.
///   - Dense     (coarse fallback)  : small dense GEMM (m < 2*(p - 1)).
template <class T>
class StructuredNullspaceBasis
{
public:
    enum class Mode { Identity, BlockDiag, Dense };

    StructuredNullspaceBasis()
        : m_mode(Mode::Identity), m_m(0), m_nD0(0), m_kb(0), m_kNull(0), m_nMid(0)
    {}

    void setIdentity(index_t m)
    {
        m_mode = Mode::Identity;
        m_m = m;
        m_nD0 = m;
        m_kb = m_kNull = m_nMid = 0;
        m_VL.resize(0, 0);
        m_VR.resize(0, 0);
        m_dense.resize(0, 0);
    }

    void setBlockDiag(index_t m, index_t kb, const gsMatrix<T> &VL, const gsMatrix<T> &VR)
    {
        m_mode = Mode::BlockDiag;
        m_m = m;
        m_kb = kb;
        m_kNull = VL.cols();
        m_nMid = m - 2 * kb;
        m_nD0 = 2 * m_kNull + m_nMid;
        m_VL = VL;
        m_VR = VR;
        m_dense.resize(0, 0);
    }

    void setDense(const gsMatrix<T> &Q0_dense)
    {
        m_mode = Mode::Dense;
        m_m = Q0_dense.rows();
        m_nD0 = Q0_dense.cols();
        m_kb = m_kNull = m_nMid = 0;
        m_VL.resize(0, 0);
        m_VR.resize(0, 0);
        m_dense = Q0_dense;
    }

    bool isIdentity() const { return m_mode == Mode::Identity; }
    index_t rows() const { return m_m; }
    index_t cols() const { return m_nD0; }

    /// Build sparse representation once during setup for M_{0,0} = Q_0^T M_0 Q_0.
    gsSparseMatrix<T> toSparseMatrix() const
    {
        gsSparseEntries<T> se;
        if (m_mode == Mode::Identity)
        {
            se.reserve(m_m);
            for (index_t i = 0; i < m_m; ++i)
                se.add(i, i, T(1));
        }
        else if (m_mode == Mode::BlockDiag)
        {
            se.reserve(2 * m_kb * m_kNull + m_nMid);
            for (index_t c = 0; c < m_kNull; ++c)
                for (index_t r = 0; r < m_kb; ++r)
                    if (std::abs(m_VL(r, c)) > T(1e-14))
                        se.add(r, c, m_VL(r, c));
            for (index_t i = 0; i < m_nMid; ++i)
                se.add(m_kb + i, m_kNull + i, T(1));
            for (index_t c = 0; c < m_kNull; ++c)
                for (index_t r = 0; r < m_kb; ++r)
                    if (std::abs(m_VR(r, c)) > T(1e-14))
                        se.add(m_m - m_kb + r, m_kNull + m_nMid + c, m_VR(r, c));
        }
        else
        {
            return m_dense.sparseView(T(1), T(1e-14));
        }
        gsSparseMatrix<T> Q(m_m, m_nD0);
        Q.setFrom(se);
        return Q;
    }

    /// out (nD0 x c) = Q_0^T * R (m x c)
    void applyTransposeLeft(const gsMatrix<T> &R, gsMatrix<T> &out) const
    {
        if (m_mode == Mode::Identity)
        {
            out = R;
        }
        else if (m_mode == Mode::BlockDiag)
        {
            out.resize(m_nD0, R.cols());
            if (m_kNull > 0)
            {
                out.topRows(m_kNull).noalias() = m_VL.transpose() * R.topRows(m_kb);
                out.bottomRows(m_kNull).noalias() = m_VR.transpose() * R.bottomRows(m_kb);
            }
            if (m_nMid > 0)
                out.middleRows(m_kNull, m_nMid) = R.middleRows(m_kb, m_nMid);
        }
        else
        {
            out.resize(m_nD0, R.cols());
            out.noalias() = m_dense.transpose() * R;
        }
    }

    /// out (m x c) = Q_0 * W (nD0 x c)
    void applyLeft(const gsMatrix<T> &W, gsMatrix<T> &out) const
    {
        if (m_mode == Mode::Identity)
        {
            out = W;
        }
        else if (m_mode == Mode::BlockDiag)
        {
            out.resize(m_m, W.cols());
            if (m_kNull > 0)
            {
                out.topRows(m_kb).noalias() = m_VL * W.topRows(m_kNull);
                out.bottomRows(m_kb).noalias() = m_VR * W.bottomRows(m_kNull);
            }
            else
            {
                out.topRows(m_kb).setZero();
                out.bottomRows(m_kb).setZero();
            }
            if (m_nMid > 0)
                out.middleRows(m_kb, m_nMid) = W.middleRows(m_kNull, m_nMid);
        }
        else
        {
            out.resize(m_m, W.cols());
            out.noalias() = m_dense * W;
        }
    }

    /// out (r x nD0) = M (r x m) * Q_0
    void applyRight(const gsMatrix<T> &M, gsMatrix<T> &out) const
    {
        if (m_mode == Mode::Identity)
        {
            out = M;
        }
        else if (m_mode == Mode::BlockDiag)
        {
            out.resize(M.rows(), m_nD0);
            if (m_kNull > 0)
            {
                out.leftCols(m_kNull).noalias() = M.leftCols(m_kb) * m_VL;
                out.rightCols(m_kNull).noalias() = M.rightCols(m_kb) * m_VR;
            }
            if (m_nMid > 0)
                out.middleCols(m_kNull, m_nMid) = M.middleCols(m_kb, m_nMid);
        }
        else
        {
            out.resize(M.rows(), m_nD0);
            out.noalias() = M * m_dense;
        }
    }

    /// out (r x m) = U (r x nD0) * Q_0^T
    void applyTransposeRight(const gsMatrix<T> &U, gsMatrix<T> &out) const
    {
        if (m_mode == Mode::Identity)
        {
            out = U;
        }
        else if (m_mode == Mode::BlockDiag)
        {
            out.resize(U.rows(), m_m);
            if (m_kNull > 0)
            {
                out.leftCols(m_kb).noalias() = U.leftCols(m_kNull) * m_VL.transpose();
                out.rightCols(m_kb).noalias() = U.rightCols(m_kNull) * m_VR.transpose();
            }
            else
            {
                out.leftCols(m_kb).setZero();
                out.rightCols(m_kb).setZero();
            }
            if (m_nMid > 0)
                out.middleCols(m_kb, m_nMid) = U.middleCols(m_kNull, m_nMid);
        }
        else
        {
            out.resize(U.rows(), m_m);
            out.noalias() = U * m_dense.transpose();
        }
    }

private:
    Mode m_mode;
    index_t m_m, m_nD0, m_kb, m_kNull, m_nMid;
    gsMatrix<T> m_VL, m_VR;
    gsMatrix<T> m_dense;
};

/// Descriptor for a patch's clamped interior (m1 x m2) tensor B-spline DOF grid.
struct PatchInteriorGrid
{
    index_t m1 = 0;
    index_t m2 = 0;
    bool hasFree = false;
    bool allFree = false;
    std::vector<index_t> dofs; // column-major: dofs[l1 + m1 * l2], -1 if constrained
};

/// Build a patch's clamped interior DOF grid in O(m1 * m2) without set_difference overhead.
template <class T>
PatchInteriorGrid extractPatchInteriorGrid(
    const gsTensorBSplineBasis<2, T> &tb,
    const gsArgyrisEmbedding<T> &emb,
    const gsDofMapper &mapper,
    size_t patchIdx)
{
    PatchInteriorGrid grid;
    const index_t nInt = emb.sizes[0];
    const index_t n1 = tb.component(0).size();
    const index_t n2 = tb.component(1).size();
    const index_t m1 = n1 - 4;
    const index_t m2 = n2 - 4;
    if (nInt <= 0 || m1 <= 0 || m2 <= 0)
        return grid;

    grid.m1 = m1;
    grid.m2 = m2;
    const index_t nPatch = m1 * m2;
    grid.dofs.assign(nPatch, -1);

    index_t nFreeCount = 0;
    if (nInt == nPatch)
    {
        // Standard AS-G1 layout: interior local DOF loc = l1 + m1 * l2 maps to (l1 + 2, l2 + 2)
        for (index_t loc = 0; loc < nPatch; ++loc)
        {
            if (mapper.is_free(loc, patchIdx))
            {
                grid.dofs[loc] = mapper.index(loc, patchIdx);
                ++nFreeCount;
            }
        }
    }
    else
    {
        const std::vector<index_t> intDofs = getInteriorDofs(tb);
        for (index_t loc = 0; loc < nInt; ++loc)
        {
            if (!mapper.is_free(loc, patchIdx))
                continue;
            const index_t tensorIdx = intDofs[loc];
            const index_t l1 = (tensorIdx % n1) - 2;
            const index_t l2 = (tensorIdx / n1) - 2;
            if (l1 >= 0 && l1 < m1 && l2 >= 0 && l2 < m2)
            {
                grid.dofs[l1 + m1 * l2] = mapper.index(loc, patchIdx);
                ++nFreeCount;
            }
        }
    }

    grid.hasFree = (nFreeCount > 0);
    grid.allFree = (nFreeCount == nPatch);
    return grid;
}

/// Cache key identifying a univariate B-spline basis (degree + knot sequence).
template <class T>
struct UnivariateBasisKey
{
    short_t degree;
    std::vector<T> knots;

    static UnivariateBasisKey fromBasis(const gsBSplineBasis<T> &b1d)
    {
        UnivariateBasisKey key;
        key.degree = b1d.degree();
        key.knots.assign(b1d.knots().begin(), b1d.knots().end());
        return key;
    }

    bool operator==(const UnivariateBasisKey &other) const
    {
        if (degree != other.degree || knots.size() != other.knots.size())
            return false;
        for (size_t i = 0; i < knots.size(); ++i)
            if (std::abs(knots[i] - other.knots[i]) > T(1e-14))
                return false;
        return true;
    }
};

/// Assemble 1D clamped interior mass M0 and 2nd-derivative biharmonic B0 (m x m, m = n - 4)
/// together with the minimum knot span hMin.
template <class T>
void assemble1DClampedMassAndBiharmonic(
    const gsBSplineBasis<T> &b1d,
    gsSparseMatrix<T> &M0,
    gsSparseMatrix<T> &B0,
    T &hMin)
{
    const index_t n = b1d.size();
    const index_t m = n - 4;
    const short_t deg = b1d.degree();

    const std::vector<T> uknots = b1d.knots().unique();
    hMin = (uknots.size() >= 2) ? (uknots[1] - uknots[0]) : T(1);
    for (size_t k = 1; k + 1 < uknots.size(); ++k)
    {
        const T span = uknots[k + 1] - uknots[k];
        if (span > T(0) && span < hMin)
            hMin = span;
    }
    if (hMin <= T(0))
        hMin = T(1);

    gsSparseEntries<T> seM, seB;
    seM.reserve(m * (2 * deg + 1));
    seB.reserve(m * (2 * deg + 1));

    gsGaussRule<T> quRule(deg + 1);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;
    std::vector<gsMatrix<T>> ev;
    gsMatrix<index_t> actives;

    for (size_t k = 0; k + 1 < uknots.size(); ++k)
    {
        quRule.mapTo(uknots[k], uknots[k + 1], quNodes, quWeights);
        b1d.evalAllDers_into(quNodes, 2, ev);
        b1d.active_into(quNodes.col(0), actives);

        const index_t nAct = actives.rows();
        const index_t nQ = quWeights.rows();
        for (index_t r = 0; r < nAct; ++r)
        {
            const index_t ir = actives(r, 0) - 2;
            if (ir < 0 || ir >= m)
                continue;
            for (index_t c = 0; c < nAct; ++c)
            {
                const index_t ic = actives(c, 0) - 2;
                if (ic < 0 || ic >= m)
                    continue;
                T mVal = T(0), bVal = T(0);
                for (index_t q = 0; q < nQ; ++q)
                {
                    const T w = quWeights(q, 0);
                    mVal += w * ev[0](r, q) * ev[0](c, q);
                    bVal += w * ev[2](r, q) * ev[2](c, q);
                }
                seM.add(ir, ic, mVal);
                seB.add(ir, ic, bVal);
            }
        }
    }

    M0.resize(m, m);
    B0.resize(m, m);
    M0.setFrom(seM);
    B0.setFrom(seB);
}

/// Summary of patch physical pullback metric sampled on an nSample x nSample midpoint grid.
template <class T>
struct PatchMetric
{
    T c1 = T(1);
    T c2 = T(1);
    T detMean = T(1);
    T cBarG = T(1);
};

/// Vectorized batch evaluation of patch Jacobian pullback metric coefficients.
template <class T>
PatchMetric<T> samplePatchMetric(const gsGeometry<T> &geo,
                                 bool useMetricWeights,
                                 int nSample = 4)
{
    const gsMatrix<T> pr = geo.parameterRange();
    const index_t nPts = nSample * nSample;
    gsMatrix<T> pts(2, nPts);
    index_t idx = 0;
    for (int sy = 0; sy < nSample; ++sy)
    {
        const T v = pr(1, 0) + (T(sy) + T(0.5)) / T(nSample) * (pr(1, 1) - pr(1, 0));
        for (int sx = 0; sx < nSample; ++sx)
        {
            pts(0, idx) = pr(0, 0) + (T(sx) + T(0.5)) / T(nSample) * (pr(0, 1) - pr(0, 0));
            pts(1, idx) = v;
            ++idx;
        }
    }

    gsMatrix<T> Jall;
    geo.jacobian_into(pts, Jall); // 2 x (2 * nPts), block k at cols [2*k, 2*k+1]

    T sumC1 = T(0), sumC2 = T(0), sumDet = T(0), maxBound = T(0);
    for (index_t k = 0; k < nPts; ++k)
    {
        const T j00 = Jall(0, 2 * k);
        const T j01 = Jall(0, 2 * k + 1);
        const T j10 = Jall(1, 2 * k);
        const T j11 = Jall(1, 2 * k + 1);

        const T g00 = j00 * j00 + j10 * j10;
        const T g01 = j00 * j01 + j10 * j11;
        const T g11 = j01 * j01 + j11 * j11;
        const T detG = std::max(g00 * g11 - g01 * g01, T(1e-30));
        const T absDetJ = std::sqrt(detG);
        const T invDetG = T(1) / detG;
        const T ginv11 = g11 * invDetG;
        const T ginv12 = -g01 * invDetG;
        const T ginv22 = g00 * invDetG;

        sumC1 += absDetJ * ginv11 * ginv11;
        sumC2 += absDetJ * ginv22 * ginv22;
        sumDet += absDetJ;
        const T frobSq = ginv11 * ginv11 + T(2) * ginv12 * ginv12 + ginv22 * ginv22;
        maxBound = std::max(maxBound, absDetJ * frobSq);
    }

    PatchMetric<T> pm;
    const T invN = T(1) / T(nPts);
    pm.detMean = sumDet * invN;
    pm.cBarG = std::max(maxBound, T(1e-12));
    if (useMetricWeights)
    {
        pm.c1 = std::max(sumC1 * invN, T(1e-12));
        pm.c2 = std::max(sumC2 * invN, T(1e-12));
    }
    else
    {
        pm.c1 = T(1);
        pm.c2 = T(1);
    }
    return pm;
}

/// Calibrate damping factor `safetyFactor / lambda_max(L^{-1} A_pp)` via single-SpMV
/// A_pp-inner-product power iteration with early relative convergence stopping.
template <class T, class ApplyOpFn>
T calibratePatchOperatorScale(const gsSparseMatrix<T> &A_pp,
                              const PatchInteriorGrid &grid,
                              ApplyOpFn &&applyOp,
                              T safetyFactor,
                              int maxIter = 15,
                              T relTol = T(1e-3))
{
    const index_t m1 = grid.m1;
    const index_t m2 = grid.m2;
    const index_t nPatch = m1 * m2;

    gsVector<T> v(nPatch);
    for (index_t i = 0; i < nPatch; ++i)
        v(i) = (grid.dofs[i] >= 0)
            ? (std::sin(T(i + 1) * T(1.37)) + T(0.5) * std::cos(T(i + 1) * T(0.23)))
            : T(0);

    gsVector<T> Av = A_pp * v;
    T vAv = v.dot(Av);
    if (vAv <= T(1e-30))
        return T(1);

    const T invInit = T(1) / std::sqrt(vAv);
    v *= invInit;
    Av *= invInit;

    T lamMax = T(1);
    T prevLam = T(0);
    gsMatrix<T> R_mat(m1, m2), E_mat(m1, m2);
    gsVector<T> z(nPatch), Az(nPatch);

    for (int iter = 0; iter < maxIter; ++iter)
    {
        std::copy(Av.data(), Av.data() + nPatch, R_mat.data());
        applyOp(R_mat, E_mat);
        std::copy(E_mat.data(), E_mat.data() + nPatch, z.data());
        if (!grid.allFree)
        {
            for (index_t i = 0; i < nPatch; ++i)
                if (grid.dofs[i] < 0)
                    z(i) = T(0);
        }

        const T zAv = z.dot(Av); // since v^T A v == 1, Rayleigh quotient is simply z^T (A v)
        if (zAv > T(0))
            lamMax = zAv;

        if (iter >= 4 && prevLam > T(0) && std::abs(lamMax - prevLam) <= relTol * lamMax)
            break;
        prevLam = lamMax;

        // Single SpMV per iteration: reuse Az as next Av after A-normalization
        Az.noalias() = A_pp * z;
        const T zAz = z.dot(Az);
        if (zAz <= T(1e-30))
            break;
        const T invNorm = T(1) / std::sqrt(zAz);
        v.noalias() = invNorm * z;
        Av.noalias() = invNorm * Az;
    }

    return (lamMax > T(1e-12)) ? (safetyFactor / lamMax) : T(1);
}

/// Gather a patch interior's residual (rhs - A * x) into column-major R_mat (m1 x m2).
template <class T>
void gatherPatchResidual(const gsSparseMatrix<T> &A,
                         const PatchInteriorGrid &grid,
                         const gsMatrix<T> &rhs,
                         const gsMatrix<T> &x,
                         index_t col,
                         bool xIsZero,
                         gsMatrix<T> &R_mat)
{
    const index_t nPatch = grid.m1 * grid.m2;
    const index_t nRows = rhs.rows();
    const T *rhsCol = rhs.data() + col * nRows;
    const T *xCol = x.data() + col * nRows;
    T *rPtr = R_mat.data();
    const index_t *dofs = grid.dofs.data();

    if (grid.allFree)
    {
        if (xIsZero)
        {
            for (index_t i = 0; i < nPatch; ++i)
                rPtr[i] = rhsCol[dofs[i]];
        }
        else
        {
            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = dofs[i];
                T acc = rhsCol[row];
                for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
                    acc -= it.value() * xCol[it.row()];
                rPtr[i] = acc;
            }
        }
    }
    else
    {
        for (index_t i = 0; i < nPatch; ++i)
        {
            const index_t row = dofs[i];
            if (row < 0)
            {
                rPtr[i] = T(0);
                continue;
            }
            T acc = rhsCol[row];
            if (!xIsZero)
                for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
                    acc -= it.value() * xCol[it.row()];
            rPtr[i] = acc;
        }
    }
}

/// Scatter a weighted patch correction x[dof] += weight * E_mat[i].
template <class T>
void scatterPatchCorrection(const PatchInteriorGrid &grid,
                            const gsMatrix<T> &E_mat,
                            T weight,
                            index_t col,
                            gsMatrix<T> &x)
{
    const index_t nPatch = grid.m1 * grid.m2;
    T *xCol = x.data() + col * x.rows();
    const T *ePtr = E_mat.data();
    const index_t *dofs = grid.dofs.data();

    if (grid.allFree)
    {
        if (weight == T(1))
        {
            for (index_t i = 0; i < nPatch; ++i)
                xCol[dofs[i]] += ePtr[i];
        }
        else
        {
            for (index_t i = 0; i < nPatch; ++i)
                xCol[dofs[i]] += weight * ePtr[i];
        }
    }
    else
    {
        for (index_t i = 0; i < nPatch; ++i)
        {
            const index_t row = dofs[i];
            if (row >= 0)
                xCol[row] += weight * ePtr[i];
        }
    }
}

/// Check whether distinct patch interior DOF sets have zero mutual coupling in A.
/// Guaranteed true on all AS-G1 Galerkin levels; enables safe OpenMP patch parallelism.
template <class T>
bool arePatchInteriorsDecoupled(const gsSparseMatrix<T> &A,
                                const std::vector<const PatchInteriorGrid *> &grids)
{
    if (grids.size() <= 1)
        return true;

    std::vector<int> owner(A.rows(), -1);
    for (size_t p = 0; p < grids.size(); ++p)
        for (index_t dof : grids[p]->dofs)
            if (dof >= 0)
                owner[dof] = static_cast<int>(p);

    for (size_t p = 0; p < grids.size(); ++p)
    {
        const int pid = static_cast<int>(p);
        for (index_t row : grids[p]->dofs)
        {
            if (row < 0)
                continue;
            for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
            {
                const int colOwner = owner[it.row()];
                if (colOwner >= 0 && colOwner != pid)
                    return false;
            }
        }
    }
    return true;
}

} // namespace internal

// ============================================================================
// Section 2: Public DOF & Block Collectors
// ============================================================================

/**
 * @brief Collect global free DOF indices on the AS-G1 interface/vertex skeleton
 *        augmented with an O(p)-layer interior B-spline boundary strip (Section 5.1).
 *
 * @param r_p  Total B-spline strip index width per patch side:
 *             - r_p < 0 (default): uses r_p = p + 1 (minimal p-robust strip of Remark 5,
 *               including all p - 1 interior B-spline layers touching dOmega_hat;
 *               equals 4 for p = 3, matching collectCoupledRing1Dofs).
 *             - r_p == 2: non-overlapping skeleton V_h^{skel} (Remark 5).
 *             - r_p >= 3: custom strip width (e.g., (4*c0 + 3)*p from Theorem 4).
 */
template <class T>
std::vector<index_t> collectCoupledStripDofs(
    const gsMultiPatch<T> &mp,
    const gsMultiBasis<T> &dbasis,
    const std::vector<gsArgyrisEmbedding<T>> &argBasis,
    const gsDofMapper &mapper,
    index_t r_p = -1,
    bool includeOuterBoundaryRing = true)
{
    const size_t nPatches = mp.nPatches();
    const index_t nFree = mapper.freeSize();
    std::vector<bool> inSkeleton(nFree, false);

    // Precompute active interface sides once in O(N_ifc) instead of per-patch scanning
    std::vector<bool> sideActive(nPatches * 5, includeOuterBoundaryRing);
    if (!includeOuterBoundaryRing)
    {
        for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
        {
            sideActive[static_cast<size_t>(it->first().patch)  * 5 + it->first().m_index]  = true;
            sideActive[static_cast<size_t>(it->second().patch) * 5 + it->second().m_index] = true;
        }
    }

    for (size_t p = 0; p < nPatches; ++p)
    {
        const index_t nInt = argBasis[p].sizes[0];
        const index_t nLocal = argBasis[p].matrix.cols();
        const gsTensorBSplineBasis<2, T> &tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
        const index_t n1 = tb.component(0).size();
        const index_t n2 = tb.component(1).size();
        const index_t m1 = n1 - 4;
        const index_t m2 = n2 - 4;
        const index_t r1 = (r_p < 0) ? (tb.degree(0) + 1) : std::max<index_t>(2, r_p);
        const index_t r2 = (r_p < 0) ? (tb.degree(1) + 1) : std::max<index_t>(2, r_p);

        const bool actW = sideActive[p * 5 + 1];
        const bool actE = sideActive[p * 5 + 2];
        const bool actS = sideActive[p * 5 + 3];
        const bool actN = sideActive[p * 5 + 4];

        // Interior DOFs within the strip [2..r_d-1] or [n_d-r_d..n_d-3]
        if (nInt > 0 && (r1 > 2 || r2 > 2))
        {
            if (m1 > 0 && m2 > 0 && nInt == m1 * m2)
            {
                for (index_t loc = 0; loc < nInt; ++loc)
                {
                    if (!mapper.is_free(loc, p))
                        continue;
                    const index_t i1 = (loc % m1) + 2;
                    const index_t i2 = (loc / m1) + 2;
                    if (((i1 < r1) && actW) || ((i1 >= n1 - r1) && actE) ||
                        ((i2 < r2) && actS) || ((i2 >= n2 - r2) && actN))
                        inSkeleton[mapper.index(loc, p)] = true;
                }
            }
            else
            {
                const std::vector<index_t> intDofs = getInteriorDofs(tb);
                for (index_t loc = 0; loc < nInt; ++loc)
                {
                    if (!mapper.is_free(loc, p))
                        continue;
                    const index_t tensorIdx = intDofs[loc];
                    const index_t i1 = tensorIdx % n1;
                    const index_t i2 = tensorIdx / n1;
                    if (((i1 < r1) && actW) || ((i1 >= n1 - r1) && actE) ||
                        ((i2 < r2) && actS) || ((i2 >= n2 - r2) && actN))
                        inSkeleton[mapper.index(loc, p)] = true;
                }
            }
        }

        // All non-interior free edge and corner Argyris DOFs
        for (index_t loc = nInt; loc < nLocal; ++loc)
            if (mapper.is_free(loc, p))
                inSkeleton[mapper.index(loc, p)] = true;
    }

    std::vector<index_t> result;
    result.reserve(nFree / 4);
    for (index_t i = 0; i < nFree; ++i)
        if (inSkeleton[i])
            result.push_back(i);
    return result;
}

/**
 * @brief Backwards-compatible wrapper for fixed 2-layer interior strip (r_p = 4).
 */
template <class T>
std::vector<index_t> collectCoupledRing1Dofs(
    const gsMultiPatch<T> &mp,
    const gsMultiBasis<T> &dbasis,
    const std::vector<gsArgyrisEmbedding<T>> &argBasis,
    const gsDofMapper &mapper,
    bool includeOuterBoundaryRing = true)
{
    return collectCoupledStripDofs(mp, dbasis, argBasis, mapper, 4, includeOuterBoundaryRing);
}

/**
 * @brief Build overlapping/paired 1D patch-interior strip blocks (u-lines and v-lines)
 *        followed by the Coupled Ring-1 skeleton block for alternating line-block smoothing.
 */
template <class T>
std::vector<std::vector<index_t>> collectInteriorLineAndRing1Blocks(
    const gsMultiPatch<T> &mp,
    const gsMultiBasis<T> &dbasis,
    const std::vector<gsArgyrisEmbedding<T>> &argBasis,
    const gsDofMapper &mapper,
    index_t stripWidth = 2,
    index_t stripStep = 2,
    bool includeOuterBoundaryRing = true)
{
    std::vector<std::vector<index_t>> uBlocks;
    std::vector<std::vector<index_t>> vBlocks;

    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        const gsTensorBSplineBasis<2, T> &tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
        const internal::PatchInteriorGrid grid =
            internal::extractPatchInteriorGrid(tb, argBasis[p], mapper, p);
        if (!grid.hasFree)
            continue;

        const index_t m1 = grid.m1;
        const index_t m2 = grid.m2;

        const index_t w2 = std::min(stripWidth, m2);
        index_t lastJ0 = -1;
        for (index_t j = 0; j < m2; j += stripStep)
        {
            const index_t j0 = std::min(j, m2 - w2);
            if (j0 == lastJ0)
                continue;
            lastJ0 = j0;
            std::vector<index_t> blk;
            blk.reserve(m1 * w2);
            for (index_t jj = j0; jj < j0 + w2; ++jj)
                for (index_t ii = 0; ii < m1; ++ii)
                {
                    const index_t g = grid.dofs[ii + m1 * jj];
                    if (g >= 0)
                        blk.push_back(g);
                }
            if (!blk.empty())
                uBlocks.push_back(std::move(blk));
        }

        const index_t w1 = std::min(stripWidth, m1);
        index_t lastI0 = -1;
        for (index_t i = 0; i < m1; i += stripStep)
        {
            const index_t i0 = std::min(i, m1 - w1);
            if (i0 == lastI0)
                continue;
            lastI0 = i0;
            std::vector<index_t> blk;
            blk.reserve(w1 * m2);
            for (index_t ii = i0; ii < i0 + w1; ++ii)
                for (index_t jj = 0; jj < m2; ++jj)
                {
                    const index_t g = grid.dofs[ii + m1 * jj];
                    if (g >= 0)
                        blk.push_back(g);
                }
            if (!blk.empty())
                vBlocks.push_back(std::move(blk));
        }
    }

    std::vector<std::vector<index_t>> blocks;
    blocks.reserve(uBlocks.size() + vBlocks.size() + 1);
    for (auto &b : uBlocks) blocks.push_back(std::move(b));
    for (auto &b : vBlocks) blocks.push_back(std::move(b));

    std::vector<index_t> ring1 = collectCoupledRing1Dofs(
        mp, dbasis, argBasis, mapper, includeOuterBoundaryRing);
    if (!ring1.empty())
        blocks.push_back(std::move(ring1));

    return blocks;
}

// ============================================================================
// Section 3: Public Smoother Operators
// ============================================================================

/**
 * @brief Multiplicative Point Gauss-Seidel + Coupled Skeleton Ring-1 Subspace Smoother.
 *
 *  - step(rhs, x):  forward point Gauss-Seidel -> exact skeleton Cholesky correction.
 *  - stepT(rhs, x): exact skeleton Cholesky correction -> backward point Gauss-Seidel.
 */
template <class T>
class gsAsG1SkeletonSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::shared_ptr<gsAsG1SkeletonSmootherOp> Ptr;
    typedef memory::unique_ptr<gsAsG1SkeletonSmootherOp> uPtr;
    typedef gsPreconditionerOp<T> Base;

    gsAsG1SkeletonSmootherOp(const gsSparseMatrix<T> &mat,
                             const std::vector<index_t> &skeletonDofs)
        : m_mat(mat)
    {
        m_gs = makeGaussSeidelOp(m_mat);
        std::vector<index_t> globToLoc(m_mat.rows(), -1);
        m_skeleton.init(m_mat, skeletonDofs, globToLoc, 24);
    }

    static uPtr make(const gsSparseMatrix<T> &mat,
                     const std::vector<index_t> &skeletonDofs)
    {
        return uPtr(new gsAsG1SkeletonSmootherOp(mat, skeletonDofs));
    }

    void relaxSkeleton(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        m_skeleton.relax(m_mat, rhs, x);
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        m_gs->step(rhs, x);
        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        m_skeleton.relax(m_mat, rhs, x, internal::isExactlyZero(x));
        m_gs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    typename gsPreconditionerOp<T>::uPtr m_gs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
};

/**
 * @brief Multiplicative Block Gauss-Seidel Smoother with Exact Block Cholesky Solves.
 *
 *  - step(rhs, x):  relaxes blocks 0 .. K-1 in forward order.
 *  - stepT(rhs, x): relaxes blocks K-1 .. 0 in reverse order (exact A-adjoint).
 */
template <class T>
class gsAsG1BlockSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::shared_ptr<gsAsG1BlockSmootherOp> Ptr;
    typedef memory::unique_ptr<gsAsG1BlockSmootherOp> uPtr;
    typedef gsPreconditionerOp<T> Base;

    gsAsG1BlockSmootherOp(const gsSparseMatrix<T> &mat,
                          const std::vector<std::vector<index_t>> &blocks)
        : m_mat(mat)
    {
        std::vector<index_t> globToLoc(m_mat.rows(), -1);
        m_blocks.reserve(blocks.size());
        for (const auto &rawDofs : blocks)
        {
            internal::SubspaceCholeskyBlock<T> blk;
            if (blk.init(m_mat, rawDofs, globToLoc, 16))
                m_blocks.push_back(std::move(blk));
        }
    }

    static uPtr make(const gsSparseMatrix<T> &mat,
                     const std::vector<std::vector<index_t>> &blocks)
    {
        return uPtr(new gsAsG1BlockSmootherOp(mat, blocks));
    }

    void relaxBlock(size_t b, const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        m_blocks[b].relax(m_mat, rhs, x);
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        const bool xZero0 = internal::isExactlyZero(x);
        for (size_t b = 0; b < m_blocks.size(); ++b)
            m_blocks[b].relax(m_mat, rhs, x, xZero0 && b == 0);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        const bool xZero0 = internal::isExactlyZero(x);
        for (size_t b = m_blocks.size(); b-- > 0;)
            m_blocks[b].relax(m_mat, rhs, x, xZero0 && b + 1 == m_blocks.size());
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    std::vector<internal::SubspaceCholeskyBlock<T>> m_blocks;
};

/**
 * @brief Patch-Interior Fast Diagonalization (FD) + Coupled Ring-1 Skeleton Smoother
 *        for p-robust and anisotropic biharmonic multigrid smoothing.
 *
 *  - step(rhs, x):  [fwd Point GS ->] Patch FD sweeps -> Ring-1 Cholesky
 *  - stepT(rhs, x): Ring-1 Cholesky -> Patch FD sweeps [-> bwd Point GS]
 */
template <class T>
class gsAsG1FastDiagSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::shared_ptr<gsAsG1FastDiagSmootherOp> Ptr;
    typedef memory::unique_ptr<gsAsG1FastDiagSmootherOp> uPtr;
    typedef gsPreconditionerOp<T> Base;

    gsAsG1FastDiagSmootherOp(const gsSparseMatrix<T> &mat,
                             const gsMultiPatch<T> &mp,
                             const gsMultiBasis<T> &dbasis,
                             const std::vector<gsArgyrisEmbedding<T>> &argBasis,
                             const gsDofMapper &mapper,
                             T beta = T(0),
                             bool usePointGS = true,
                             bool includeOuterBoundaryRing = true)
        : m_mat(mat), m_usePointGS(usePointGS)
    {
        if (m_usePointGS)
            m_gs = makeGaussSeidelOp(m_mat);

        std::vector<index_t> globToLoc(m_mat.rows(), -1);

        // 1. Setup Coupled Ring-1 skeleton exact Cholesky block
        const std::vector<index_t> skelDofs = collectCoupledRing1Dofs(
            mp, dbasis, argBasis, mapper, includeOuterBoundaryRing);
        m_skeleton.init(m_mat, skelDofs, globToLoc, 24);

        // 2. Setup per-patch Fast Diagonalization operators (caching 1D eigenproblems by knot vector)
        std::vector<std::pair<internal::UnivariateBasisKey<T>, std::shared_ptr<UnivariateFD>>> fdCache;
        auto getUnivariateFD = [&](const gsBSplineBasis<T> &b1d) -> std::shared_ptr<UnivariateFD> {
            const auto key = internal::UnivariateBasisKey<T>::fromBasis(b1d);
            for (const auto &entry : fdCache)
                if (entry.first == key)
                    return entry.second;

            auto u = std::make_shared<UnivariateFD>();
            gsSparseMatrix<T> M_sp, B_sp;
            T hMin = T(1);
            internal::assemble1DClampedMassAndBiharmonic(b1d, M_sp, B_sp, hMin);
            const gsMatrix<T> M_int = M_sp.toDense();
            const gsMatrix<T> B_int = B_sp.toDense();
            typename gsMatrix<T>::GenSelfAdjEigenSolver ges(B_int, M_int, gsEigen::ComputeEigenvectors);
            u->Q = ges.eigenvectors();
            u->lam = ges.eigenvalues().cwiseMax(T(0));
            fdCache.emplace_back(key, u);
            return u;
        };

        m_patches.reserve(mp.nPatches());
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsTensorBSplineBasis<2, T> &tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
            PatchFD fd;
            fd.grid = internal::extractPatchInteriorGrid(tb, argBasis[p], mapper, p);
            if (!fd.grid.hasFree)
                continue;

            fd.u1 = getUnivariateFD(tb.component(0));
            fd.u2 = getUnivariateFD(tb.component(1));

            const internal::PatchMetric<T> pm = internal::samplePatchMetric(mp.patch(p), true, 3);
            const T c0 = beta * pm.detMean;

            fd.invDiag.resize(fd.grid.m1, fd.grid.m2);
            for (index_t b = 0; b < fd.grid.m2; ++b)
            {
                const T s2 = std::sqrt(pm.c2 * fd.u2->lam(b));
                for (index_t a = 0; a < fd.grid.m1; ++a)
                {
                    const T s1 = std::sqrt(pm.c1 * fd.u1->lam(a));
                    const T denom = (s1 + s2) * (s1 + s2) + c0;
                    fd.invDiag(a, b) = (denom > T(0)) ? (T(1) / denom) : T(0);
                }
            }

            const gsSparseMatrix<T> A_pp =
                internal::extractSymmetricSubmatrix(m_mat, fd.grid.dofs, globToLoc, 25);

            Workspace ws;
            const T scale = internal::calibratePatchOperatorScale(
                A_pp, fd.grid,
                [&](const gsMatrix<T> &R, gsMatrix<T> &E) { applyPatchFD(fd, R, E, ws); },
                T(0.90), /*maxIter=*/12);
            fd.invDiag *= scale;

            m_patches.push_back(std::move(fd));
        }

        std::vector<const internal::PatchInteriorGrid *> gridPtrs;
        gridPtrs.reserve(m_patches.size());
        for (const auto &fd : m_patches)
            gridPtrs.push_back(&fd.grid);
        m_patchesDecoupled = internal::arePatchInteriorsDecoupled(m_mat, gridPtrs);

#ifdef _OPENMP
        m_workspaces.resize(static_cast<size_t>(omp_get_max_threads()));
#else
        m_workspaces.resize(1);
#endif
    }

    static uPtr make(const gsSparseMatrix<T> &mat,
                     const gsMultiPatch<T> &mp,
                     const gsMultiBasis<T> &dbasis,
                     const std::vector<gsArgyrisEmbedding<T>> &argBasis,
                     const gsDofMapper &mapper,
                     T beta = T(0),
                     bool usePointGS = true,
                     bool includeOuterBoundaryRing = true)
    {
        return uPtr(new gsAsG1FastDiagSmootherOp(
            mat, mp, dbasis, argBasis, mapper, beta, usePointGS, includeOuterBoundaryRing));
    }

    void relaxPatch(size_t pIdx, const gsMatrix<T> &rhs, gsMatrix<T> &x, bool xIsZero = false) const
    {
#ifdef _OPENMP
        Workspace &ws = m_workspaces[static_cast<size_t>(omp_get_thread_num())];
#else
        Workspace &ws = m_workspaces[0];
#endif
        const PatchFD &fd = m_patches[pIdx];
        const index_t ncols = rhs.cols();
        ws.R_mat.resize(fd.grid.m1, fd.grid.m2);
        ws.E_mat.resize(fd.grid.m1, fd.grid.m2);

        for (index_t col = 0; col < ncols; ++col)
        {
            internal::gatherPatchResidual(m_mat, fd.grid, rhs, x, col, xIsZero, ws.R_mat);
            applyPatchFD(fd, ws.R_mat, ws.E_mat, ws);
            internal::scatterPatchCorrection(fd.grid, ws.E_mat, T(1), col, x);
        }
    }

    void relaxSkeleton(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        m_skeleton.relax(m_mat, rhs, x);
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_usePointGS)
            m_gs->step(rhs, x);

        const bool xZero = internal::isExactlyZero(x);
#ifdef _OPENMP
        if (m_patchesDecoupled && m_patches.size() > 1)
        {
#pragma omp parallel for schedule(static)
            for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
                relaxPatch(static_cast<size_t>(p), rhs, x, xZero);
        }
        else
#endif
        {
            for (size_t p = 0; p < m_patches.size(); ++p)
                relaxPatch(p, rhs, x, xZero && (m_patchesDecoupled || p == 0));
        }

        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        m_skeleton.relax(m_mat, rhs, x, internal::isExactlyZero(x));

        const bool xZero = internal::isExactlyZero(x);
#ifdef _OPENMP
        if (m_patchesDecoupled && m_patches.size() > 1)
        {
#pragma omp parallel for schedule(static)
            for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
                relaxPatch(static_cast<size_t>(p), rhs, x, xZero);
        }
        else
#endif
        {
            for (size_t p = m_patches.size(); p-- > 0;)
                relaxPatch(p, rhs, x, xZero && (m_patchesDecoupled || p + 1 == m_patches.size()));
        }

        if (m_usePointGS)
            m_gs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    struct UnivariateFD
    {
        gsMatrix<T> Q;
        gsVector<T> lam;
    };

    struct PatchFD
    {
        internal::PatchInteriorGrid grid;
        std::shared_ptr<UnivariateFD> u1, u2;
        gsMatrix<T> invDiag;
    };

    struct Workspace
    {
        gsMatrix<T> R_mat, T_mat, Z_mat, E_mat;
    };

    static void applyPatchFD(const PatchFD &fd,
                             const gsMatrix<T> &R_mat,
                             gsMatrix<T> &E_mat,
                             Workspace &ws)
    {
        ws.T_mat.noalias() = fd.u1->Q.transpose() * R_mat;
        ws.Z_mat.noalias() = ws.T_mat * fd.u2->Q;
        ws.Z_mat.array() *= fd.invDiag.array();
        ws.T_mat.noalias() = fd.u1->Q * ws.Z_mat;
        E_mat.noalias() = ws.T_mat * fd.u2->Q.transpose();
    }

    gsSparseMatrix<T> m_mat;
    bool m_usePointGS;
    bool m_patchesDecoupled = true;
    typename gsPreconditionerOp<T>::uPtr m_gs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
    std::vector<PatchFD> m_patches;
    mutable std::vector<Workspace> m_workspaces;
};

/**
 * @brief O(p N_h) h- and p-robust Subspace-Corrected Mass Smoother (SCMS) + Overlapping
 *        Skeleton Schwarz Smoother for multi-patch AS-G1 biharmonic discretizations.
 *
 * On each patch interior S^0_{p,p-2,h}(Omega_hat), splits each 1D spline space into:
 *   - S^{D,0} (vanishing even boundary derivatives u^{(2j)} = 0 for 2j <= p - 2),
 *     inverted in O(p m) via StructuredNullspaceBasis Q_0 and zero-transpose
 *     BandedSPDCholesky of M_{0,0} = Q_0^T M_0 Q_0.
 *   - S^{D,1} = (S^{D,0})^{\perp_{W_0}} (dim <= p - 2), simultaneously diagonalizing
 *     B_0 and M_0.
 *
 * Coupled with an exact sparse Cholesky solve on the enlarged skeleton \widetilde{V}_h^{skel}:
 *  - Multiplicative mode (additive = false, default):
 *      step(rhs, x):  [fwd Point GS ->] patch-interior SCMS -> skeleton Cholesky
 *      stepT(rhs, x): skeleton Cholesky -> patch-interior SCMS [-> bwd Point GS]
 *  - Hybrid mode (usePointGS = true, Sogn & Takacs 2019 Section 5.3):
 *      wraps Patch SCMS + Skeleton in symmetric Point Gauss-Seidel sweeps to eliminate
 *      geometry-distortion degradation on strongly curved/stretched multi-patch domains.
 *  - Additive mode (additive = true, Theorem 5):
 *      x <- x + 0.5 * ( sum_i scale_i I_i (L_i^{int})^{-1} I_i^T r + I_skel A_skel^{-1} I_skel^T r )
 */
template <class T>
class gsAsG1ScmsSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::shared_ptr<gsAsG1ScmsSmootherOp> Ptr;
    typedef memory::unique_ptr<gsAsG1ScmsSmootherOp> uPtr;
    typedef gsPreconditionerOp<T> Base;

    gsAsG1ScmsSmootherOp(const gsSparseMatrix<T> &mat,
                         const gsMultiPatch<T> &mp,
                         const gsMultiBasis<T> &dbasis,
                         const std::vector<gsArgyrisEmbedding<T>> &argBasis,
                         const gsDofMapper &mapper,
                         T beta = T(0),
                         index_t r_p = -1,
                         bool additive = false,
                         bool useMetricWeights = true,
                         bool calibrateByPowerMethod = true,
                         bool includeOuterBoundaryRing = true,
                         bool usePointGS = false)
        : m_mat(mat), m_additive(additive), m_usePointGS(usePointGS)
    {
        if (m_usePointGS)
            m_gs = makeGaussSeidelOp(m_mat);

        std::vector<index_t> globToLoc(m_mat.rows(), -1);

        // 1. Setup enlarged skeleton Cholesky solver on \widetilde{V}_h^{skel}
        const std::vector<index_t> skelDofs = collectCoupledStripDofs(
            mp, dbasis, argBasis, mapper, r_p, includeOuterBoundaryRing);
        m_skeleton.init(m_mat, skelDofs, globToLoc, 32);

        // 2. Setup O(p N_h^{(i)}) Kronecker SCMS on each patch interior, caching 1D bases
        //    (W_0 = c_d (B_0 + sigma0Hat * M_0) factors out c_d completely!)
        std::vector<std::pair<internal::UnivariateBasisKey<T>, std::shared_ptr<UnivariateSCMS>>> uCache;
        auto getUnivariateSCMS = [&](const gsBSplineBasis<T> &b1d) -> std::shared_ptr<UnivariateSCMS> {
            const auto key = internal::UnivariateBasisKey<T>::fromBasis(b1d);
            for (const auto &entry : uCache)
                if (entry.first == key)
                    return entry.second;

            auto u = std::make_shared<UnivariateSCMS>();
            buildUnivariateSCMS(b1d, /*useSharpSigma0=*/calibrateByPowerMethod, *u);
            uCache.emplace_back(key, u);
            return u;
        };

        m_patches.reserve(mp.nPatches());
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsTensorBSplineBasis<2, T> &tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
            PatchSCMS scms;
            scms.grid = internal::extractPatchInteriorGrid(tb, argBasis[p], mapper, p);
            if (!scms.grid.hasFree)
                continue;

            const internal::PatchMetric<T> pm =
                internal::samplePatchMetric(mp.patch(p), useMetricWeights, 4);
            scms.c0 = beta * pm.detMean;

            scms.u1 = getUnivariateSCMS(tb.component(0));
            scms.u2 = getUnivariateSCMS(tb.component(1));

            const T sigma1 = pm.c1 * scms.u1->sigma0Hat;
            const T sigma2 = pm.c2 * scms.u2->sigma0Hat;
            scms.invSigma00 = T(1) / (sigma1 + sigma2 + scms.c0);

            if (scms.u1->nD1 > 0)
            {
                scms.invDiag10.resize(scms.u1->nD1);
                for (index_t a = 0; a < scms.u1->nD1; ++a)
                {
                    const T denom = pm.c1 * scms.u1->b11Hat(a)
                                  + (sigma2 + scms.c0) * scms.u1->mu11(a);
                    scms.invDiag10(a) = (denom > T(0)) ? (T(1) / denom) : T(0);
                }
            }
            if (scms.u2->nD1 > 0)
            {
                scms.invDiag01.resize(scms.u2->nD1);
                for (index_t b = 0; b < scms.u2->nD1; ++b)
                {
                    const T denom = pm.c2 * scms.u2->b11Hat(b)
                                  + (sigma1 + scms.c0) * scms.u2->mu11(b);
                    scms.invDiag01(b) = (denom > T(0)) ? (T(1) / denom) : T(0);
                }
            }
            if (scms.u1->nD1 > 0 && scms.u2->nD1 > 0)
            {
                scms.invDiag11.resize(scms.u1->nD1, scms.u2->nD1);
                for (index_t b = 0; b < scms.u2->nD1; ++b)
                    for (index_t a = 0; a < scms.u1->nD1; ++a)
                    {
                        const T denom = pm.c1 * scms.u1->b11Hat(a) * scms.u2->mu11(b)
                                      + scms.u1->mu11(a) * pm.c2 * scms.u2->b11Hat(b)
                                      + scms.c0 * scms.u1->mu11(a) * scms.u2->mu11(b);
                        scms.invDiag11(a, b) = (denom > T(0)) ? (T(1) / denom) : T(0);
                    }
            }

            if (calibrateByPowerMethod)
            {
                const gsSparseMatrix<T> A_pp =
                    internal::extractSymmetricSubmatrix(m_mat, scms.grid.dofs, globToLoc, 25);
                Workspace ws;
                scms.scale = internal::calibratePatchOperatorScale(
                    A_pp, scms.grid,
                    [&](const gsMatrix<T> &R, gsMatrix<T> &E) { applyUnscaledPatchSCMS(scms, R, E, ws); },
                    T(0.92), /*maxIter=*/15);
            }
            else
            {
                // Analytical bound (8 * c_bar_G)^{-1} from Theorem 5 / Eq. (55)
                scms.scale = T(1) / (T(8) * pm.cBarG);
            }

            m_patches.push_back(std::move(scms));
        }

        std::vector<const internal::PatchInteriorGrid *> gridPtrs;
        gridPtrs.reserve(m_patches.size());
        for (const auto &scms : m_patches)
            gridPtrs.push_back(&scms.grid);
        m_patchesDecoupled = internal::arePatchInteriorsDecoupled(m_mat, gridPtrs);

#ifdef _OPENMP
        m_workspaces.resize(static_cast<size_t>(omp_get_max_threads()));
#else
        m_workspaces.resize(1);
#endif
    }

    static uPtr make(const gsSparseMatrix<T> &mat,
                     const gsMultiPatch<T> &mp,
                     const gsMultiBasis<T> &dbasis,
                     const std::vector<gsArgyrisEmbedding<T>> &argBasis,
                     const gsDofMapper &mapper,
                     T beta = T(0),
                     index_t r_p = -1,
                     bool additive = false,
                     bool useMetricWeights = true,
                     bool calibrateByPowerMethod = true,
                     bool includeOuterBoundaryRing = true,
                     bool usePointGS = false)
    {
        return uPtr(new gsAsG1ScmsSmootherOp(
            mat, mp, dbasis, argBasis, mapper, beta, r_p,
            additive, useMetricWeights, calibrateByPowerMethod,
            includeOuterBoundaryRing, usePointGS));
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_additive)
        {
            if (m_usePointGS)
                m_gs->step(rhs, x);
            stepAdditive(rhs, x);
            if (m_usePointGS)
                m_gs->stepT(rhs, x);
            return;
        }

        if (m_usePointGS)
            m_gs->step(rhs, x);

        const bool xZero = internal::isExactlyZero(x);
#ifdef _OPENMP
        if (m_patchesDecoupled && m_patches.size() > 1)
        {
#pragma omp parallel for schedule(static)
            for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
                relaxPatch(static_cast<size_t>(p), rhs, x, xZero);
        }
        else
#endif
        {
            for (size_t p = 0; p < m_patches.size(); ++p)
                relaxPatch(p, rhs, x, xZero && (m_patchesDecoupled || p == 0));
        }

        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_additive)
        {
            step(rhs, x);
            return;
        }

        m_skeleton.relax(m_mat, rhs, x, internal::isExactlyZero(x));

        const bool xZero = internal::isExactlyZero(x);
#ifdef _OPENMP
        if (m_patchesDecoupled && m_patches.size() > 1)
        {
#pragma omp parallel for schedule(static)
            for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
                relaxPatch(static_cast<size_t>(p), rhs, x, xZero);
        }
        else
#endif
        {
            for (size_t p = m_patches.size(); p-- > 0;)
                relaxPatch(p, rhs, x, xZero && (m_patchesDecoupled || p + 1 == m_patches.size()));
        }

        if (m_usePointGS)
            m_gs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    struct UnivariateSCMS
    {
        index_t m = 0;
        index_t nD0 = 0;
        index_t nD1 = 0;
        T sigma0Hat = T(1);                       // unit-metric inverse constant (physical sigma0 = c_d * sigma0Hat)
        internal::StructuredNullspaceBasis<T> Q0; // m x nD0 localized basis of S^{D,0}
        internal::BandedSPDCholesky<T> cholM00;   // narrow-banded Cholesky of M_{0,0} = Q_0^T M_0 Q_0
        gsMatrix<T> Q1;                           // m x nD1 basis of S^{D,1} (diagonalizing B_0, M_0)
        gsVector<T> b11Hat;                       // unit-metric q_a^T B_0 q_a
        gsVector<T> mu11;                         // q_a^T M_0 q_a
    };

    struct PatchSCMS
    {
        internal::PatchInteriorGrid grid;
        std::shared_ptr<UnivariateSCMS> u1, u2;
        T c0 = T(0);
        T invSigma00 = T(1);
        gsVector<T> invDiag10;
        gsVector<T> invDiag01;
        gsMatrix<T> invDiag11;
        T scale = T(1);
    };

    struct Workspace
    {
        gsMatrix<T> R_mat, E_mat;
        gsMatrix<T> T0, T1;
        gsMatrix<T> W00, W10, W01, W11;
        gsMatrix<T> U0, U1;
    };

    /// Estimate lambda_max(M_{0,0}^{-1} B_{0,0}) in O(p * nD0) using banded Cholesky power iteration.
    static T estimate1DInverseConstant(const internal::BandedSPDCholesky<T> &cholM00,
                                       const gsSparseMatrix<T> &B00,
                                       T fallbackSigma)
    {
        const index_t n = cholM00.size();
        if (n <= 0)
            return fallbackSigma;

        gsVector<T> v(n);
        for (index_t i = 0; i < n; ++i)
            v(i) = std::sin(T(i + 1) * T(1.73)) + T(0.5) * std::cos(T(i + 1) * T(0.37));
        const T vn = v.norm();
        if (vn <= T(0))
            return fallbackSigma;
        v /= vn;

        T lamMax = T(0);
        T prevLam = T(0);
        gsVector<T> Bv(n), w(n);
        for (int iter = 0; iter < 20; ++iter)
        {
            Bv.noalias() = B00 * v;
            const T vBv = v.dot(Bv);
            if (vBv <= T(1e-30))
                break;

            w = Bv;
            cholM00.solveLeftInPlace(w); // w = M_{0,0}^{-1} (B_{0,0} v)
            const T wBv = w.dot(Bv);
            const T lam = (wBv > T(0)) ? (wBv / vBv) : T(0);
            if (lam > T(0))
                lamMax = lam;

            const T wn = w.norm();
            if (wn <= T(1e-30))
                break;
            v = w / wn;

            if (iter >= 5 && prevLam > T(0) && std::abs(lamMax - prevLam) <= T(1e-4) * lamMax)
                break;
            prevLam = lamMax;
        }
        return (lamMax > T(1e-12)) ? (T(1.02) * lamMax) : fallbackSigma;
    }

    static void buildUnivariateSCMS(const gsBSplineBasis<T> &b1d,
                                    bool useSharpSigma0,
                                    UnivariateSCMS &out)
    {
        const short_t deg = b1d.degree();
        const index_t n = b1d.size();
        const index_t m = n - 4;
        out.m = m;

        gsSparseMatrix<T> M0, B0;
        T hMin = T(1);
        internal::assemble1DClampedMassAndBiharmonic(b1d, M0, B0, hMin);

        // Lemma 4 analytical bound on S^{D,0}_{p,p-2,h}
        const T analyticSigma0 = T(2520) / (hMin * hMin * hMin * hMin);

        const index_t r_ev = (deg - 2) / 2; // even derivatives 2j <= p - 2 per end
        if (r_ev <= 0)
        {
            out.nD0 = m;
            out.nD1 = 0;
            out.Q0.setIdentity(m);
            out.cholM00.compute(M0);
            out.sigma0Hat = useSharpSigma0
                ? estimate1DInverseConstant(out.cholM00, B0, analyticSigma0)
                : analyticSigma0;
            return;
        }

        // Evaluate boundary even derivatives on hMin-scaled knots for high-p stability
        const std::vector<T> uknots = b1d.knots().unique();
        const T uMin = uknots.front();
        const T uMaxScaled = (uknots.back() - uMin) / hMin;

        std::vector<T> scaledKnots;
        scaledKnots.reserve(b1d.knots().size());
        for (typename gsKnotVector<T>::const_iterator it = b1d.knots().begin();
             it != b1d.knots().end(); ++it)
            scaledKnots.push_back((*it - uMin) / hMin);

        gsKnotVector<T> kvScaled(deg, scaledKnots.begin(), scaledKnots.end());
        gsBSplineBasis<T> bScaled(kvScaled);

        gsMatrix<T> endPts(1, 2);
        endPts(0, 0) = T(0);
        endPts(0, 1) = uMaxScaled;
        std::vector<gsMatrix<T>> evEnd;
        gsMatrix<index_t> actEnd;
        bScaled.evalAllDers_into(endPts, static_cast<int>(2 * r_ev), evEnd);
        bScaled.active_into(endPts, actEnd);

        const index_t kb = deg - 1;
        gsMatrix<T> C_full;

        if (m >= 2 * kb)
        {
            gsMatrix<T> C_L(r_ev, kb), C_R(r_ev, kb);
            C_L.setZero();
            C_R.setZero();

            for (index_t j = 1; j <= r_ev; ++j)
            {
                for (index_t r = 0; r < actEnd.rows(); ++r)
                {
                    const index_t irL = actEnd(r, 0) - 2;
                    if (irL >= 0 && irL < kb)
                        C_L(j - 1, irL) = evEnd[2 * j](r, 0);

                    const index_t irR = actEnd(r, 1) - 2;
                    if (irR >= m - kb && irR < m)
                        C_R(j - 1, irR - (m - kb)) = evEnd[2 * j](r, 1);
                }
                const T nL = C_L.row(j - 1).norm();
                if (nL > T(0)) C_L.row(j - 1) /= nL;
                const T nR = C_R.row(j - 1).norm();
                if (nR > T(0)) C_R.row(j - 1) /= nR;
            }

            const gsMatrix<T> V_L = C_L.jacobiSvd(gsEigen::ComputeFullV).matrixV();
            const gsMatrix<T> V_R = C_R.jacobiSvd(gsEigen::ComputeFullV).matrixV();
            const index_t kNull = kb - r_ev;
            out.Q0.setBlockDiag(m, kb, V_L.rightCols(kNull), V_R.rightCols(kNull));
            out.nD0 = out.Q0.cols();

            C_full.setZero(2 * r_ev, m);
            C_full.topLeftCorner(r_ev, kb)     = C_L;
            C_full.bottomRightCorner(r_ev, kb) = C_R;
        }
        else
        {
            // Coarse-mesh fallback (m < 2*(p - 1))
            gsMatrix<T> C_raw(2 * r_ev, m);
            C_raw.setZero();
            index_t nRows = 0;
            for (index_t j = 1; j <= r_ev; ++j)
            {
                for (int side = 0; side < 2; ++side)
                {
                    gsVector<T> rowVec(m);
                    rowVec.setZero();
                    for (index_t r = 0; r < actEnd.rows(); ++r)
                    {
                        const index_t ir = actEnd(r, side) - 2;
                        if (ir >= 0 && ir < m)
                            rowVec(ir) = evEnd[2 * j](r, side);
                    }
                    const T rn = rowVec.norm();
                    if (rn > T(1e-14))
                        C_raw.row(nRows++) = (rowVec / rn).transpose();
                }
            }
            C_raw.conservativeResize(nRows, m);

            auto svd = C_raw.jacobiSvd(gsEigen::ComputeFullU | gsEigen::ComputeFullV);
            const gsVector<T> sv = svd.singularValues();
            index_t rk = 0;
            if (sv.size() > 0 && sv(0) > T(0))
                for (index_t i = 0; i < sv.size(); ++i)
                    if (sv(i) > T(1e-10) * sv(0))
                        ++rk;

            out.nD0 = m - rk;
            out.Q0.setDense(svd.matrixV().rightCols(out.nD0));
            C_full = svd.matrixU().leftCols(rk).transpose() * C_raw;
        }

        if (out.nD0 > 0)
        {
            const gsSparseMatrix<T> Q0_sp = out.Q0.toSparseMatrix();
            const gsSparseMatrix<T> M00 = Q0_sp.transpose() * M0 * Q0_sp;
            out.cholM00.compute(M00);

            if (useSharpSigma0)
            {
                const gsSparseMatrix<T> B00 = Q0_sp.transpose() * B0 * Q0_sp;
                out.sigma0Hat = estimate1DInverseConstant(out.cholM00, B00, analyticSigma0);
            }
            else
            {
                out.sigma0Hat = analyticSigma0;
            }
        }
        else
        {
            out.sigma0Hat = analyticSigma0;
        }

        // Construct S^{D,1} = (S^{D,0})^{\perp_{W0Hat}} = span(W0Hat^{-1} C_full^T)
        if (C_full.rows() == 0)
        {
            out.nD1 = 0;
            return;
        }

        const gsSparseMatrix<T> W0Hat = B0 + out.sigma0Hat * M0;
        internal::BandedSPDCholesky<T> cholW0;
        cholW0.compute(W0Hat);

        gsMatrix<T> Q1_tilde = C_full.transpose();
        cholW0.solveLeftInPlace(Q1_tilde);

        // Normalize columns in W0Hat-norm and W0Hat-orthonormalize via small symmetric eigensolve
        gsMatrix<T> WQ = W0Hat * Q1_tilde;
        for (index_t k = 0; k < Q1_tilde.cols(); ++k)
        {
            const T wNormSq = Q1_tilde.col(k).dot(WQ.col(k));
            if (wNormSq > T(0))
            {
                const T invNorm = T(1) / std::sqrt(wNormSq);
                Q1_tilde.col(k) *= invNorm;
                WQ.col(k) *= invNorm;
            }
        }

        const gsMatrix<T> G_W = Q1_tilde.transpose() * WQ;
        gsEigen::SelfAdjointEigenSolver<typename gsMatrix<T>::Base> esW(G_W, gsEigen::ComputeEigenvectors);
        const gsVector<T> ew = esW.eigenvalues();
        const T ewMax = (ew.size() > 0) ? ew.maxCoeff() : T(0);

        std::vector<index_t> keep;
        for (index_t k = 0; k < ew.size(); ++k)
            if (ew(k) > T(1e-11) * ewMax)
                keep.push_back(k);

        out.nD1 = static_cast<index_t>(keep.size());
        if (out.nD1 == 0)
            return;

        gsMatrix<T> Q1_hat(m, out.nD1);
        for (index_t idx = 0; idx < out.nD1; ++idx)
        {
            const index_t k = keep[idx];
            Q1_hat.col(idx).noalias() = (Q1_tilde * esW.eigenvectors().col(k)) / std::sqrt(ew(k));
        }

        // Simultaneously diagonalize M_{1,1} and B_{1,1} inside the W0Hat-orthonormal basis
        const gsMatrix<T> M11_hat = Q1_hat.transpose() * (M0 * Q1_hat);
        gsEigen::SelfAdjointEigenSolver<typename gsMatrix<T>::Base> esM(M11_hat, gsEigen::ComputeEigenvectors);
        out.Q1.noalias() = Q1_hat * esM.eigenvectors();

        const gsMatrix<T> BQ1 = B0 * out.Q1;
        const gsMatrix<T> MQ1 = M0 * out.Q1;
        out.b11Hat.resize(out.nD1);
        out.mu11.resize(out.nD1);
        for (index_t a = 0; a < out.nD1; ++a)
        {
            out.b11Hat(a) = std::max(T(0), out.Q1.col(a).dot(BQ1.col(a)));
            out.mu11(a)   = std::max(T(1e-30), out.Q1.col(a).dot(MQ1.col(a)));
        }
    }

    /// Factored 4-block Kronecker SCMS apply with zero heap allocations and zero transposes.
    static void applyUnscaledPatchSCMS(const PatchSCMS &scms,
                                       const gsMatrix<T> &R_mat,
                                       gsMatrix<T> &E_mat,
                                       Workspace &ws)
    {
        const UnivariateSCMS &u1 = *scms.u1;
        const UnivariateSCMS &u2 = *scms.u2;

        // Ultra-fast p = 3 path (S^{D,1} = {0}, Q_0 = I): pure banded Cholesky left + right solves
        if (u1.nD1 == 0 && u2.nD1 == 0 && u1.Q0.isIdentity() && u2.Q0.isIdentity())
        {
            E_mat = R_mat;
            u1.cholM00.solveLeftInPlace(E_mat);
            u2.cholM00.solveRightInPlace(E_mat);
            E_mat *= scms.invSigma00;
            return;
        }

        // 1. Left projections (each computed ONCE instead of twice)
        if (u1.nD0 > 0)
            u1.Q0.applyTransposeLeft(R_mat, ws.T0);
        if (u1.nD1 > 0)
        {
            ws.T1.resize(u1.nD1, u2.m);
            ws.T1.noalias() = u1.Q1.transpose() * R_mat;
        }

        // 2. Subspace (*, 0): combine Blocks (0,0) and (1,0) before single right Q_0^{(2)T}
        if (u2.nD0 > 0)
        {
            if (u1.nD0 > 0)
            {
                u2.Q0.applyRight(ws.T0, ws.W00);
                u1.cholM00.solveLeftInPlace(ws.W00);
                u2.cholM00.solveRightInPlace(ws.W00);
                ws.W00 *= scms.invSigma00;
                u1.Q0.applyLeft(ws.W00, ws.U0);
            }
            else
            {
                ws.U0.setZero(u1.m, u2.nD0);
            }

            if (u1.nD1 > 0)
            {
                u2.Q0.applyRight(ws.T1, ws.W10);
                for (index_t a = 0; a < u1.nD1; ++a)
                    ws.W10.row(a) *= scms.invDiag10(a);
                u2.cholM00.solveRightInPlace(ws.W10);
                ws.U0.noalias() += u1.Q1 * ws.W10;
            }

            u2.Q0.applyTransposeRight(ws.U0, E_mat);
        }
        else
        {
            E_mat.setZero(u1.m, u2.m);
        }

        // 3. Subspace (*, 1): combine Blocks (0,1) and (1,1) before single right Q_1^{(2)T}
        if (u2.nD1 > 0)
        {
            if (u1.nD0 > 0)
            {
                ws.W01.resize(u1.nD0, u2.nD1);
                ws.W01.noalias() = ws.T0 * u2.Q1;
                u1.cholM00.solveLeftInPlace(ws.W01);
                for (index_t b = 0; b < u2.nD1; ++b)
                    ws.W01.col(b) *= scms.invDiag01(b);
                u1.Q0.applyLeft(ws.W01, ws.U1);
            }
            else
            {
                ws.U1.setZero(u1.m, u2.nD1);
            }

            if (u1.nD1 > 0)
            {
                ws.W11.resize(u1.nD1, u2.nD1);
                ws.W11.noalias() = ws.T1 * u2.Q1;
                ws.W11.array() *= scms.invDiag11.array();
                ws.U1.noalias() += u1.Q1 * ws.W11;
            }

            E_mat.noalias() += ws.U1 * u2.Q1.transpose();
        }
    }

    void relaxPatch(size_t pIdx, const gsMatrix<T> &rhs, gsMatrix<T> &x, bool xIsZero = false) const
    {
#ifdef _OPENMP
        Workspace &ws = m_workspaces[static_cast<size_t>(omp_get_thread_num())];
#else
        Workspace &ws = m_workspaces[0];
#endif
        const PatchSCMS &scms = m_patches[pIdx];
        const index_t ncols = rhs.cols();
        ws.R_mat.resize(scms.grid.m1, scms.grid.m2);
        ws.E_mat.resize(scms.grid.m1, scms.grid.m2);

        for (index_t col = 0; col < ncols; ++col)
        {
            internal::gatherPatchResidual(m_mat, scms.grid, rhs, x, col, xIsZero, ws.R_mat);
            applyUnscaledPatchSCMS(scms, ws.R_mat, ws.E_mat, ws);
            internal::scatterPatchCorrection(scms.grid, ws.E_mat, scms.scale, col, x);
        }
    }

    void stepAdditive(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        const bool xZero = internal::isExactlyZero(x);
        const gsMatrix<T> *rPtr = &rhs;
        if (!xZero)
        {
            m_resBuf.resize(rhs.rows(), rhs.cols());
            m_resBuf.noalias() = rhs - m_mat * x;
            rPtr = &m_resBuf;
        }
        const gsMatrix<T> &r = *rPtr;
        const index_t ncols = rhs.cols();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (m_patches.size() > 1)
#endif
        for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
        {
#ifdef _OPENMP
            Workspace &ws = m_workspaces[static_cast<size_t>(omp_get_thread_num())];
#else
            Workspace &ws = m_workspaces[0];
#endif
            const PatchSCMS &scms = m_patches[static_cast<size_t>(p)];
            const index_t nPatch = scms.grid.m1 * scms.grid.m2;
            ws.R_mat.resize(scms.grid.m1, scms.grid.m2);
            ws.E_mat.resize(scms.grid.m1, scms.grid.m2);
            const T w = T(0.5) * scms.scale;

            for (index_t col = 0; col < ncols; ++col)
            {
                const T *rCol = r.data() + col * r.rows();
                T *rmPtr = ws.R_mat.data();
                const index_t *dofs = scms.grid.dofs.data();
                if (scms.grid.allFree)
                {
                    for (index_t i = 0; i < nPatch; ++i)
                        rmPtr[i] = rCol[dofs[i]];
                }
                else
                {
                    for (index_t i = 0; i < nPatch; ++i)
                    {
                        const index_t row = dofs[i];
                        rmPtr[i] = (row >= 0) ? rCol[row] : T(0);
                    }
                }

                applyUnscaledPatchSCMS(scms, ws.R_mat, ws.E_mat, ws);
                internal::scatterPatchCorrection(scms.grid, ws.E_mat, w, col, x);
            }
        }

        m_skeleton.applyAdditive(r, x, T(0.5));
    }

    gsSparseMatrix<T> m_mat;
    bool m_additive;
    bool m_usePointGS;
    bool m_patchesDecoupled = true;
    typename gsPreconditionerOp<T>::uPtr m_gs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
    std::vector<PatchSCMS> m_patches;
    mutable std::vector<Workspace> m_workspaces;
    mutable gsMatrix<T> m_resBuf;
};

} // namespace gismo
