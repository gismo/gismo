/** @file gsAsG1Smoother.hpp

    @brief Specialized O(N) smoothers for AS-G1 multi-patch biharmonic multigrid solvers.

    Provides four A-symmetric smoother families sharing unified subspace/banded primitives:
      1. gsAsG1SkeletonSmootherOp  : Point Gauss-Seidel + Coupled Ring-1 skeleton Cholesky.
      2. gsAsG1BlockSmootherOp     : Multiplicative block Gauss-Seidel (alternating interior
                                     line strips + Coupled Ring-1 skeleton).
      3. gsAsG1FastDiagSmootherOp  : Patch-interior Fast Diagonalization + Ring-1 skeleton
                                     (with optional Point Gauss-Seidel hybrid sweep).
      4. gsAsG1ScmsSmootherOp      : O(p N_h) h- and p-robust Subspace-Corrected Mass Smoother
                                     (SCMS) based on the scaled-norm L2-orthogonal odd-derivative
                                     S_0 splitting (single pure Kronecker mass block on S_0 x S_0
                                     + exact sparse or L2-orthogonal skeleton Cholesky), with
                                     fast branchless Point Gauss-Seidel hybrid mode.

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
        if (m_solver->info() != gsEigen::Success)
        {
            gsWarn << "SubspaceCholeskyBlock: SimplicialLDLT factorization failed (size=" << k
                   << "); disabling block.\n";
            m_solver.reset();
            m_dofs.clear();
            return false;
        }
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
        m_eSub.resize(k, ncols);

        if (ncols == 1)
        {
            const T *rhsPtr = rhs.data();
            T *xPtr = x.data();
            T *rPtr = m_rSub.data();
            const index_t *dofsPtr = m_dofs.data();
            if (xIsZero)
            {
                for (index_t i = 0; i < k; ++i)
                    rPtr[i] = rhsPtr[dofsPtr[i]];
            }
            else if (A.isCompressed())
            {
                typedef typename gsSparseMatrix<T>::StorageIndex StorageIndex;
                const StorageIndex *outer = A.outerIndexPtr();
                const StorageIndex *inner = A.innerIndexPtr();
                const T *val = A.valuePtr();
                for (index_t i = 0; i < k; ++i)
                {
                    const index_t row = dofsPtr[i];
                    T acc = rhsPtr[row];
                    const StorageIndex end = outer[row + 1];
                    for (StorageIndex p = outer[row]; p < end; ++p)
                        acc -= val[p] * xPtr[inner[p]];
                    rPtr[i] = acc;
                }
            }
            else
            {
                for (index_t i = 0; i < k; ++i)
                {
                    const index_t row = dofsPtr[i];
                    T acc = rhsPtr[row];
                    for (typename gsSparseMatrix<T>::InnerIterator it(A, row); it; ++it)
                        acc -= it.value() * xPtr[it.row()];
                    rPtr[i] = acc;
                }
            }
            m_eSub = m_solver->solve(m_rSub);
            const T *ePtr = m_eSub.data();
            for (index_t i = 0; i < k; ++i)
                xPtr[dofsPtr[i]] += ePtr[i];
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
        m_eSub.resize(k, ncols);

        if (ncols == 1)
        {
            const T *rGlob = r.data();
            T *rPtr = m_rSub.data();
            const index_t *dofsPtr = m_dofs.data();
            for (index_t i = 0; i < k; ++i)
                rPtr[i] = rGlob[dofsPtr[i]];

            m_eSub = m_solver->solve(m_rSub);
            T *xPtr = x.data();
            const T *ePtr = m_eSub.data();
            for (index_t i = 0; i < k; ++i)
                xPtr[dofsPtr[i]] += weight * ePtr[i];
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

/// Fast, multi-RHS, reference-free Point Gauss-Seidel smoother for symmetric compressed
/// ColMajor sparse matrices. Caches reciprocal diagonals and diagonal positions so that:
///   - Multi-column right-hand sides are handled natively (no single-RHS restriction).
///   - No internal reference to A is stored (safe across moves/copies of owning smoothers).
///   - Inner loops are branchless on raw CSC pointers, with a triangular half-sweep fast
///     path whenever the initial iterate is identically zero.
template <class T>
class PointGaussSeidelSmoother
{
public:
    typedef typename gsSparseMatrix<T>::StorageIndex StorageIndex;

    PointGaussSeidelSmoother() : m_sortedTri(false) {}

    void init(const gsSparseMatrix<T> &A)
    {
        GISMO_ASSERT(A.rows() == A.cols(), "PointGaussSeidelSmoother requires a square matrix.");
        GISMO_ASSERT(A.isCompressed(), "PointGaussSeidelSmoother requires a compressed matrix.");

        const index_t n = A.rows();
        m_invDiag.resize(n);
        m_diagPos.assign(static_cast<size_t>(n), StorageIndex(-1));
        m_sortedTri = (n > 0);

        const StorageIndex *outer = A.outerIndexPtr();
        const StorageIndex *inner = A.innerIndexPtr();
        const T *val = A.valuePtr();

        for (index_t i = 0; i < n; ++i)
        {
            const StorageIndex start = outer[i];
            const StorageIndex end   = outer[i + 1];
            StorageIndex prev = StorageIndex(-1);
            T d = T(0);
            bool foundDiag = false;
            for (StorageIndex k = start; k < end; ++k)
            {
                const StorageIndex j = inner[k];
                if (j <= prev)
                    m_sortedTri = false;
                prev = j;
                if (static_cast<index_t>(j) == i)
                {
                    d = val[k];
                    m_diagPos[static_cast<size_t>(i)] = k;
                    foundDiag = true;
                }
            }
            if (!foundDiag)
                m_sortedTri = false;

            if (std::abs(d) > T(1e-30))
                m_invDiag(i) = T(1) / d;
            else
            {
                m_invDiag(i) = T(0);
                gsWarn << "PointGaussSeidelSmoother: zero or near-zero diagonal at row " << i << ".\n";
            }
        }
    }

    /// Forward Point Gauss-Seidel sweep (exact on symmetric ColMajor compressed A).
    void stepForward(const gsSparseMatrix<T> &A,
                     const gsMatrix<T> &rhs,
                     gsMatrix<T> &x,
                     bool xIsZero = false) const
    {
        const index_t n = A.rows();
        const index_t ncols = rhs.cols();
        if (n == 0 || ncols == 0)
            return;

        const StorageIndex *outer = A.outerIndexPtr();
        const StorageIndex *inner = A.innerIndexPtr();
        const T *val = A.valuePtr();
        const T *invD = m_invDiag.data();

        for (index_t c = 0; c < ncols; ++c)
        {
            const T *f = rhs.data() + c * n;
            T *xp = x.data() + c * n;
            if (xIsZero && m_sortedTri)
            {
                const StorageIndex *dpos = m_diagPos.data();
                for (index_t i = 0; i < n; ++i)
                {
                    T sum = T(0);
                    const StorageIndex end = dpos[i];
                    for (StorageIndex k = outer[i]; k < end; ++k)
                        sum += val[k] * xp[inner[k]];
                    xp[i] = (f[i] - sum) * invD[i];
                }
            }
            else
            {
                for (index_t i = 0; i < n; ++i)
                {
                    T sum = T(0);
                    const StorageIndex end = outer[i + 1];
                    for (StorageIndex k = outer[i]; k < end; ++k)
                        sum += val[k] * xp[inner[k]];
                    xp[i] += (f[i] - sum) * invD[i];
                }
            }
        }
    }

    /// Backward Point Gauss-Seidel sweep (exact A-adjoint of stepForward on symmetric A).
    void stepBackward(const gsSparseMatrix<T> &A,
                      const gsMatrix<T> &rhs,
                      gsMatrix<T> &x,
                      bool xIsZero = false) const
    {
        const index_t n = A.rows();
        const index_t ncols = rhs.cols();
        if (n == 0 || ncols == 0)
            return;

        const StorageIndex *outer = A.outerIndexPtr();
        const StorageIndex *inner = A.innerIndexPtr();
        const T *val = A.valuePtr();
        const T *invD = m_invDiag.data();

        for (index_t c = 0; c < ncols; ++c)
        {
            const T *f = rhs.data() + c * n;
            T *xp = x.data() + c * n;
            if (xIsZero && m_sortedTri)
            {
                const StorageIndex *dpos = m_diagPos.data();
                for (index_t i = n - 1; i >= 0; --i)
                {
                    T sum = T(0);
                    const StorageIndex end = outer[i + 1];
                    for (StorageIndex k = dpos[i] + 1; k < end; ++k)
                        sum += val[k] * xp[inner[k]];
                    xp[i] = (f[i] - sum) * invD[i];
                }
            }
            else
            {
                for (index_t i = n - 1; i >= 0; --i)
                {
                    T sum = T(0);
                    const StorageIndex end = outer[i + 1];
                    for (StorageIndex k = outer[i]; k < end; ++k)
                        sum += val[k] * xp[inner[k]];
                    xp[i] += (f[i] - sum) * invD[i];
                }
            }
        }
    }

private:
    gsVector<T> m_invDiag;
    std::vector<StorageIndex> m_diagPos;
    bool m_sortedTri;
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

/// Assemble 1D clamped interior mass M0 (m x m, m = n - 4) and optionally the full
/// n x n 1D mass matrix Mfull (when buildFullMass is true), plus minimum knot span hMin.
/// Avoids computing 2nd derivatives when only mass matrices are needed (new S_0 SCMS).
template <class T>
void assemble1DMassMatrices(
    const gsBSplineBasis<T> &b1d,
    gsSparseMatrix<T> &M0,
    T &hMin,
    bool buildFullMass = false,
    gsSparseMatrix<T> *Mfull = nullptr)
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

    gsSparseEntries<T> seM0, seMfull;
    seM0.reserve(std::max<index_t>(0, m) * (2 * deg + 1));
    if (buildFullMass && Mfull)
        seMfull.reserve(n * (2 * deg + 1));

    gsGaussRule<T> quRule(deg + 1);
    gsMatrix<T> quNodes, ev;
    gsVector<T> quWeights;
    gsMatrix<index_t> actives;

    for (size_t k = 0; k + 1 < uknots.size(); ++k)
    {
        quRule.mapTo(uknots[k], uknots[k + 1], quNodes, quWeights);
        b1d.eval_into(quNodes, ev);
        b1d.active_into(quNodes.col(0), actives);

        const index_t nAct = actives.rows();
        const index_t nQ = quWeights.rows();
        for (index_t r = 0; r < nAct; ++r)
        {
            const index_t ar = actives(r, 0);
            const index_t ir = ar - 2;
            const bool rInt = (ir >= 0 && ir < m);
            for (index_t c = 0; c < nAct; ++c)
            {
                const index_t ac = actives(c, 0);
                const index_t ic = ac - 2;
                const bool cInt = (ic >= 0 && ic < m);
                if (!rInt && !cInt && !(buildFullMass && Mfull))
                    continue;
                if ((!rInt || !cInt) && !(buildFullMass && Mfull))
                    continue;

                T mVal = T(0);
                for (index_t q = 0; q < nQ; ++q)
                    mVal += quWeights(q, 0) * ev(r, q) * ev(c, q);

                if (rInt && cInt)
                    seM0.add(ir, ic, mVal);
                if (buildFullMass && Mfull)
                    seMfull.add(ar, ac, mVal);
            }
        }
    }

    M0.resize(std::max<index_t>(0, m), std::max<index_t>(0, m));
    M0.setFrom(seM0);
    if (buildFullMass && Mfull)
    {
        Mfull->resize(n, n);
        Mfull->setFrom(seMfull);
    }
}

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
        else if (A.isCompressed())
        {
            typedef typename gsSparseMatrix<T>::StorageIndex StorageIndex;
            const StorageIndex *outer = A.outerIndexPtr();
            const StorageIndex *inner = A.innerIndexPtr();
            const T *val = A.valuePtr();
            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = dofs[i];
                T acc = rhsCol[row];
                const StorageIndex end = outer[row + 1];
                for (StorageIndex k = outer[row]; k < end; ++k)
                    acc -= val[k] * xCol[inner[k]];
                rPtr[i] = acc;
            }
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
 *        augmented with a boundary strip of interior B-splines per patch side.
 *
 * @param r_p  Total B-spline strip index width per patch side:
 *             - r_p == -1 (default): exact boundary-supported B-spline strip
 *               r_d = max(4, p_d) covering {B_0, ..., B_{p_d - 1}} (Section 5.3:
 *               {B_{p_d}, ..., B_{n_d - p_d - 1}} belong identically to S_0;
 *               equals 4 for p = 3, 4 and p for p >= 5).
 *             - r_p == -2: legacy p + 1 boundary strip.
 *             - r_p == -3: minimal odd-derivative constraint strip max(4, 2*floor((p-3)/2) + 2).
 *             - r_p == 2:  non-overlapping skeleton V_h^{skel}.
 *             - r_p >= 3:  explicit custom strip width.
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
    auto resolveStrip = [](short_t deg, index_t rp) -> index_t {
        if (rp >= 2)
            return rp;
        if (rp == -2)
            return static_cast<index_t>(deg) + 1;
        if (rp == -3)
        {
            const index_t ncTilde = std::max<index_t>(0, (static_cast<index_t>(deg) - 3) / 2);
            return std::max<index_t>(4, 2 * ncTilde + 2);
        }
        return std::max<index_t>(4, static_cast<index_t>(deg));
    };

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
        const index_t r1 = resolveStrip(tb.degree(0), r_p);
        const index_t r2 = resolveStrip(tb.degree(1), r_p);

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
        m_mat.makeCompressed();
        m_pgs.init(m_mat);
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
        m_pgs.stepForward(m_mat, rhs, x, internal::isExactlyZero(x));
        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        const bool xZero0 = internal::isExactlyZero(x);
        m_skeleton.relax(m_mat, rhs, x, xZero0);
        m_pgs.stepBackward(m_mat, rhs, x, xZero0 && m_skeleton.empty());
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    internal::PointGaussSeidelSmoother<T> m_pgs;
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
        m_mat.makeCompressed();
        if (m_usePointGS)
            m_pgs.init(m_mat);

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
        const bool xZero0 = internal::isExactlyZero(x);
        if (m_usePointGS)
            m_pgs.stepForward(m_mat, rhs, x, xZero0);

        const bool xZero = !m_usePointGS && xZero0;
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
        const bool xZero0 = internal::isExactlyZero(x);
        m_skeleton.relax(m_mat, rhs, x, xZero0);

        const bool xZero = xZero0 && m_skeleton.empty();
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
            m_pgs.stepBackward(m_mat, rhs, x, false);
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
    internal::PointGaussSeidelSmoother<T> m_pgs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
    std::vector<PatchFD> m_patches;
    mutable std::vector<Workspace> m_workspaces;
};

/**
 * @brief O(p N_h) h- and p-robust Subspace-Corrected Mass Smoother (SCMS) + Skeleton
 *        Schwarz Smoother for multi-patch AS-G1 biharmonic discretizations (Section 5).
 *
 * Based on the scaled-norm L2-orthogonal odd-derivative splitting:
 *   On each 1D spline space S = S_{p,p-2,h}(0,1), define the constrained subspace
 *     S_0 = { s in S : s^{(r)}(0) = s^{(r)}(1) = 0  for all r in R_p },
 *     R_p = {0} U { 2*l + 1 : l >= 0, 2*l + 1 <= p - 2 }  (|R_p| = floor((p+1)/2)),
 *   which satisfies S_0 subset S^0 = S cap H_0^2(0,1) and the p-robust inverse inequality
 *   |s_0|_{H^2}^2 <= 2520 h^{-4} ||s_0||_{L^2}^2 (Lemma 3).
 *
 *   - For p = 3 and p = 4: S_0 = S^0 (P_0 = I, zero boundary constraints/SVDs).
 *   - For p >= 5: only n_c_tilde = floor((p - 3)/2) odd derivatives r = 3, 5, ..., 2*n_c_tilde + 1
 *     are enforced on k_b = 2*n_c_tilde interior B-splines per endpoint via tiny SVDs.
 *   - On each patch interior V_{h,00}^{(i)} = S_0 x S_0, a SINGLE pure Kronecker mass block
 *       (L_{00}^{(i)})^{-1} = scale_i * (P_0^{(1)} M_{00,1}^{-1} P_0^{(1)T}
 *                                       (x) P_0^{(2)} M_{00,2}^{-1} P_0^{(2)T})
 *     is applied in O(p N_h^{(i)}) via StructuredNullspaceBasis and zero-transpose
 *     BandedSPDCholesky -- eliminating all Q_1, W_0, eigensolve, and 4-block overhead.
 *
 * Coupled with an exact skeleton Cholesky solver:
 *  - Practical sparse local skeleton (useL2OrthSkeleton = false, default; Remark 6(2)):
 *      exact sparse Cholesky on the AS-G1 interface/vertex skeleton plus a minimal
 *      O(p/2)-layer interior strip covering the odd-derivative constraint support.
 *  - Exact L2(Omega_hat)-orthogonal skeleton (useL2OrthSkeleton = true; Theorem 5 / Eq. (60)):
 *      exact non-overlapping direct sum V_h^{AS-G1} = (sum_i V_{h,00}^{(i)}) (+) V_h^{skel,perp}
 *      with implicit O(p N_h) application of I_skel^perp and (I_skel^perp)^T.
 *  - Multiplicative mode (additive = false, default):
 *      step(rhs, x):  [fwd Point GS ->] patch-interior S_00 mass -> skeleton Cholesky
 *      stepT(rhs, x): skeleton Cholesky -> patch-interior S_00 mass [-> bwd Point GS]
 *  - Hybrid mode (usePointGS = true, Sogn & Takacs 2019 Section 5.3):
 *      wraps Patch S_00 mass + Skeleton in fast branchless Point Gauss-Seidel sweeps to
 *      eliminate geometry-distortion / anisotropy / C^1-knot degradation.
 *  - Additive mode (additive = true, Theorem 5 / Eq. (61)):
 *      x <- x + 0.5 * ( sum_i scale_i I_{00}^{(i)} (M_{00}^{(i)})^{-1} (I_{00}^{(i)})^T r
 *                       + I_skel A_skel^{-1} I_skel^T r ).
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
                         bool usePointGS = false,
                         bool useL2OrthSkeleton = false,
                         T multSafetyFactor = T(1.50))
        : m_mat(mat),
          m_additive(additive),
          m_usePointGS(usePointGS),
          m_useL2OrthSkeleton(useL2OrthSkeleton)
    {
        m_mat.makeCompressed();

        if (m_usePointGS)
            m_pgs.init(m_mat);

        std::vector<index_t> globToLoc(m_mat.rows(), -1);

        // 1. Setup practical sparse skeleton Cholesky solver (unless exact L2-orth is requested)
        if (!m_useL2OrthSkeleton)
        {
            const std::vector<index_t> skelDofs = collectCoupledStripDofs(
                mp, dbasis, argBasis, mapper, r_p, includeOuterBoundaryRing);
            m_skeleton.init(m_mat, skelDofs, globToLoc, 32);
        }

        // 2. Setup O(p N_h^{(i)}) pure Kronecker mass smoother on each patch interior S_0 x S_0,
        //    caching unique 1D bases across patches and directions.
        std::vector<std::pair<internal::UnivariateBasisKey<T>, std::shared_ptr<UnivariateS0>>> uCache;
        auto getUnivariateS0 = [&](const gsBSplineBasis<T> &b1d) -> std::shared_ptr<UnivariateS0> {
            const auto key = internal::UnivariateBasisKey<T>::fromBasis(b1d);
            for (const auto &entry : uCache)
                if (entry.first == key)
                    return entry.second;

            auto u = std::make_shared<UnivariateS0>();
            buildUnivariateS0(b1d, /*needFullMass=*/m_useL2OrthSkeleton, *u);
            uCache.emplace_back(key, u);
            return u;
        };

        m_patches.reserve(mp.nPatches());
        m_totalInteriorDofs = 0;
        const T safetyFactor = m_additive ? T(0.92) : multSafetyFactor;

        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsTensorBSplineBasis<2, T> &tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
            PatchS00 scms;
            scms.patchIdx = p;
            scms.grid = internal::extractPatchInteriorGrid(tb, argBasis[p], mapper, p);
            if (!scms.grid.hasFree)
                continue;

            scms.u1 = getUnivariateS0(tb.component(0));
            scms.u2 = getUnivariateS0(tb.component(1));
            m_totalInteriorDofs += scms.grid.m1 * scms.grid.m2;

            if (scms.u1->n0 > 0 && scms.u2->n0 > 0)
            {
                if (calibrateByPowerMethod)
                {
                    const gsSparseMatrix<T> A_pp =
                        internal::extractSymmetricSubmatrix(m_mat, scms.grid.dofs, globToLoc, 25);
                    Workspace ws;
                    scms.scale = internal::calibratePatchOperatorScale(
                        A_pp, scms.grid,
                        [&](const gsMatrix<T> &R, gsMatrix<T> &E) { applyUnscaledPatchS00(scms, R, E, ws); },
                        safetyFactor, /*maxIter=*/20);
                }
                else
                {
                    // Exact analytical bound from Section 5.3 / Eq. (58):
                    //   L_{00}^{(i)} = (2 * cBarG * 2520 * (h1^{-4} + h2^{-4}) + beta * detMean) * (M0 x M0)
                    //   (= (4 * cBarG * 2520 * h^{-4} + beta * detMean) * (M0 x M0) when h1 == h2 == h)
                    const internal::PatchMetric<T> pm =
                        internal::samplePatchMetric(mp.patch(p), useMetricWeights, 4);
                    const T sigma12 = T(2) * (scms.u1->analyticSigma + scms.u2->analyticSigma);
                    const T denom = pm.cBarG * sigma12 + beta * pm.detMean;
                    scms.scale = (denom > T(0)) ? (T(1) / denom) : T(1);
                }
            }
            else
            {
                scms.scale = T(0);
            }

            m_patches.push_back(std::move(scms));
        }

        // 3. Optional exact L2(Omega_hat)-orthogonal skeleton (Theorem 5 / Eq. (60))
        if (m_useL2OrthSkeleton)
            initL2OrthSkeleton(mp, argBasis, mapper, globToLoc);

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
                     bool usePointGS = false,
                     bool useL2OrthSkeleton = false,
                     T multSafetyFactor = T(1.50))
    {
        return uPtr(new gsAsG1ScmsSmootherOp(
            mat, mp, dbasis, argBasis, mapper, beta, r_p,
            additive, useMetricWeights, calibrateByPowerMethod,
            includeOuterBoundaryRing, usePointGS, useL2OrthSkeleton,
            multSafetyFactor));
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_additive)
        {
            if (m_usePointGS)
                m_pgs.stepForward(m_mat, rhs, x, internal::isExactlyZero(x));
            stepAdditive(rhs, x);
            if (m_usePointGS)
                m_pgs.stepBackward(m_mat, rhs, x, false);
            return;
        }

        const bool xZero0 = internal::isExactlyZero(x);
        if (m_usePointGS)
            m_pgs.stepForward(m_mat, rhs, x, xZero0);

        const bool xZeroPatch = !m_usePointGS && xZero0;
        relaxAllPatches(rhs, x, xZeroPatch, /*reverse=*/false);

        if (m_useL2OrthSkeleton)
            relaxL2OrthSkeleton(rhs, x, false);
        else
            m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_additive)
        {
            step(rhs, x);
            return;
        }

        const bool xZero0 = internal::isExactlyZero(x);
        bool skelEmpty = false;
        if (m_useL2OrthSkeleton)
        {
            skelEmpty = !m_l2SkelSolver || (m_l2Nskel == 0);
            relaxL2OrthSkeleton(rhs, x, xZero0);
        }
        else
        {
            skelEmpty = m_skeleton.empty();
            m_skeleton.relax(m_mat, rhs, x, xZero0);
        }

        const bool xZeroPatch = xZero0 && skelEmpty;
        relaxAllPatches(rhs, x, xZeroPatch, /*reverse=*/true);

        if (m_usePointGS)
            m_pgs.stepBackward(m_mat, rhs, x, false);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    typedef gsEigen::SimplicialLDLT<gsSparseMatrix<T>> EigenLDLT;

    struct UnivariateS0
    {
        index_t n = 0;                          // full 1D B-spline dimension
        index_t m = 0;                          // clamped interior dimension n - 4
        index_t n0 = 0;                         // dim(S_0)
        index_t nPerp = 0;                      // dim(S_perp) = m - n0
        T hMin = T(1);
        T analyticSigma = T(2520);              // 2520 / hMin^4 (Lemma 3)
        internal::StructuredNullspaceBasis<T> P0; // m x n0 localized basis of S_0
        internal::BandedSPDCholesky<T> cholM00;   // narrow-banded Cholesky of M_{00} = P_0^T M_0 P_0
        // Populated only when useL2OrthSkeleton == true:
        gsSparseMatrix<T> Pperp;                // m x nPerp orthogonal complement of P_0 in S^0
        gsSparseMatrix<T> Mcross;               // n0 x n = P_0^T M_full[2:n-2, 0:n]
        gsSparseMatrix<T> McrossT;              // n x n0 = Mcross^T (compressed ColMajor)
    };

    struct PatchS00
    {
        size_t patchIdx = 0;
        internal::PatchInteriorGrid grid;
        std::shared_ptr<UnivariateS0> u1, u2;
        T scale = T(1);
    };

    struct PatchL2OrthData
    {
        size_t scmsIdx = 0;
        gsSparseMatrix<T> E_free;               // (n1 * n2) x nFreeLoc
        std::vector<index_t> freeGlobDofs;      // size nFreeLoc
    };

    struct Workspace
    {
        gsMatrix<T> R_mat, E_mat;
        gsMatrix<T> T0, W00, U0;
    };

    static void buildUnivariateS0(const gsBSplineBasis<T> &b1d,
                                  bool needFullMass,
                                  UnivariateS0 &out)
    {
        const short_t deg = b1d.degree();
        const index_t n = b1d.size();
        const index_t m = n - 4;
        out.n = n;
        out.m = std::max<index_t>(0, m);

        gsSparseMatrix<T> M0, Mfull;
        T hMin = T(1);
        internal::assemble1DMassMatrices(b1d, M0, hMin, needFullMass, needFullMass ? &Mfull : nullptr);
        out.hMin = hMin;
        out.analyticSigma = T(2520) / (hMin * hMin * hMin * hMin);

        if (m <= 0)
        {
            out.n0 = 0;
            out.nPerp = 0;
            return;
        }

        // Helper to extract M_strip = Mfull[2:n-2, 0:n] (m x n) when exact L2-orth is requested
        auto buildMcross = [&](const gsSparseMatrix<T> &P0_sp) {
            if (!needFullMass || out.n0 <= 0)
            {
                out.Mcross.resize(0, n);
                out.McrossT.resize(n, 0);
                return;
            }
            gsSparseEntries<T> seStrip;
            seStrip.reserve(Mfull.nonZeros());
            for (index_t col = 0; col < n; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(Mfull, col); it; ++it)
                {
                    const index_t ir = it.row() - 2;
                    if (ir >= 0 && ir < m)
                        seStrip.add(ir, col, it.value());
                }
            gsSparseMatrix<T> Mstrip(m, n);
            Mstrip.setFrom(seStrip);
            out.Mcross = P0_sp.transpose() * Mstrip;
            out.Mcross.makeCompressed();
            out.McrossT = out.Mcross.transpose();
            out.McrossT.makeCompressed();
        };

        // Odd boundary derivatives inside S^0: r in {3, 5, ..., 2*ncTilde + 1} <= p - 2
        const index_t ncTilde = std::max<index_t>(0, (static_cast<index_t>(deg) - 3) / 2);
        if (ncTilde == 0)
        {
            // p = 3 and p = 4: R_p = {0, 1}, so S_0 = S^0 with zero extra constraints!
            out.n0 = m;
            out.nPerp = 0;
            out.P0.setIdentity(m);
            out.cholM00.compute(M0);
            if (needFullMass)
            {
                out.Pperp.resize(m, 0);
                buildMcross(out.P0.toSparseMatrix());
            }
            return;
        }

        // Evaluate odd derivatives r = 3, 5, ..., rMax on hMin-scaled knots for high-p stability
        const index_t rMax = 2 * ncTilde + 1;
        const index_t kb = rMax - 1; // exact support width in S^0: {B_2, ..., B_{rMax}}

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
        bScaled.evalAllDers_into(endPts, static_cast<int>(rMax), evEnd);
        bScaled.active_into(endPts, actEnd);

        if (m >= 2 * kb)
        {
            gsMatrix<T> C_L(ncTilde, kb), C_R(ncTilde, kb);
            C_L.setZero();
            C_R.setZero();

            for (index_t l = 0; l < ncTilde; ++l)
            {
                const index_t der = 2 * l + 3; // 3, 5, ..., rMax
                for (index_t r = 0; r < actEnd.rows(); ++r)
                {
                    const index_t irL = actEnd(r, 0) - 2;
                    if (irL >= 0 && irL < kb)
                        C_L(l, irL) = evEnd[der](r, 0);

                    const index_t irR = actEnd(r, 1) - 2;
                    if (irR >= m - kb && irR < m)
                        C_R(l, irR - (m - kb)) = evEnd[der](r, 1);
                }
                const T nL = C_L.row(l).norm();
                if (nL > T(0)) C_L.row(l) /= nL;
                const T nR = C_R.row(l).norm();
                if (nR > T(0)) C_R.row(l) /= nR;
            }

            const gsMatrix<T> V_L = C_L.jacobiSvd(gsEigen::ComputeFullV).matrixV();
            const gsMatrix<T> V_R = C_R.jacobiSvd(gsEigen::ComputeFullV).matrixV();
            const index_t kNull = kb - ncTilde; // == ncTilde
            out.P0.setBlockDiag(m, kb, V_L.rightCols(kNull), V_R.rightCols(kNull));
            out.n0 = out.P0.cols();
            out.nPerp = 2 * ncTilde;

            if (needFullMass)
            {
                gsSparseEntries<T> sePerp;
                sePerp.reserve(2 * kb * ncTilde);
                for (index_t c = 0; c < ncTilde; ++c)
                {
                    for (index_t r = 0; r < kb; ++r)
                    {
                        if (std::abs(V_L(r, c)) > T(1e-14))
                            sePerp.add(r, c, V_L(r, c));
                        if (std::abs(V_R(r, c)) > T(1e-14))
                            sePerp.add(m - kb + r, ncTilde + c, V_R(r, c));
                    }
                }
                out.Pperp.resize(m, out.nPerp);
                out.Pperp.setFrom(sePerp);
                out.Pperp.makeCompressed();
            }
        }
        else
        {
            // Coarse-mesh fallback (m < 2 * kb)
            gsMatrix<T> C_raw(2 * ncTilde, m);
            C_raw.setZero();
            index_t nRows = 0;
            for (index_t l = 0; l < ncTilde; ++l)
            {
                const index_t der = 2 * l + 3;
                for (int side = 0; side < 2; ++side)
                {
                    gsVector<T> rowVec(m);
                    rowVec.setZero();
                    for (index_t r = 0; r < actEnd.rows(); ++r)
                    {
                        const index_t ir = actEnd(r, side) - 2;
                        if (ir >= 0 && ir < m)
                            rowVec(ir) = evEnd[der](r, side);
                    }
                    const T rn = rowVec.norm();
                    if (rn > T(1e-14))
                        C_raw.row(nRows++) = (rowVec / rn).transpose();
                }
            }
            C_raw.conservativeResize(nRows, m);

            auto svd = C_raw.jacobiSvd(gsEigen::ComputeFullV);
            const gsVector<T> sv = svd.singularValues();
            index_t rk = 0;
            if (sv.size() > 0 && sv(0) > T(0))
                for (index_t i = 0; i < sv.size(); ++i)
                    if (sv(i) > T(1e-10) * sv(0))
                        ++rk;

            out.n0 = m - rk;
            out.nPerp = rk;
            out.P0.setDense(svd.matrixV().rightCols(out.n0));
            if (needFullMass)
            {
                out.Pperp = svd.matrixV().leftCols(rk).sparseView(T(1), T(1e-14));
                out.Pperp.makeCompressed();
            }
        }

        if (out.n0 > 0)
        {
            const gsSparseMatrix<T> P0_sp = out.P0.toSparseMatrix();
            const gsSparseMatrix<T> M00 = P0_sp.transpose() * M0 * P0_sp;
            out.cholM00.compute(M00);
            if (needFullMass)
                buildMcross(P0_sp);
        }
        else if (needFullMass)
        {
            out.Mcross.resize(0, n);
            out.McrossT.resize(n, 0);
        }
    }

    /// Single-block pure Kronecker mass apply on S_0 x S_0 (zero heap allocations, zero transposes).
    static void applyUnscaledPatchS00(const PatchS00 &scms,
                                      const gsMatrix<T> &R_mat,
                                      gsMatrix<T> &E_mat,
                                      Workspace &ws)
    {
        const UnivariateS0 &u1 = *scms.u1;
        const UnivariateS0 &u2 = *scms.u2;

        if (u1.n0 <= 0 || u2.n0 <= 0)
        {
            E_mat.setZero(u1.m, u2.m);
            return;
        }

        // Ultra-fast p = 3, 4 path (P_0 = I): pure banded Cholesky left + right solves
        if (u1.P0.isIdentity() && u2.P0.isIdentity())
        {
            E_mat = R_mat;
            u1.cholM00.solveLeftInPlace(E_mat);
            u2.cholM00.solveRightInPlace(E_mat);
            return;
        }

        // General p >= 5 path: single S_0 x S_0 Kronecker mass block
        u1.P0.applyTransposeLeft(R_mat, ws.T0);
        u2.P0.applyRight(ws.T0, ws.W00);
        u1.cholM00.solveLeftInPlace(ws.W00);
        u2.cholM00.solveRightInPlace(ws.W00);
        u1.P0.applyLeft(ws.W00, ws.U0);
        u2.P0.applyTransposeRight(ws.U0, E_mat);
    }

    void relaxPatch(size_t pIdx, const gsMatrix<T> &rhs, gsMatrix<T> &x, bool xIsZero = false) const
    {
#ifdef _OPENMP
        Workspace &ws = m_workspaces[static_cast<size_t>(omp_get_thread_num())];
#else
        Workspace &ws = m_workspaces[0];
#endif
        const PatchS00 &scms = m_patches[pIdx];
        if (scms.scale == T(0))
            return;

        const index_t ncols = rhs.cols();
        ws.R_mat.resize(scms.grid.m1, scms.grid.m2);
        ws.E_mat.resize(scms.grid.m1, scms.grid.m2);

        for (index_t col = 0; col < ncols; ++col)
        {
            internal::gatherPatchResidual(m_mat, scms.grid, rhs, x, col, xIsZero, ws.R_mat);
            applyUnscaledPatchS00(scms, ws.R_mat, ws.E_mat, ws);
            internal::scatterPatchCorrection(scms.grid, ws.E_mat, scms.scale, col, x);
        }
    }

    void relaxAllPatches(const gsMatrix<T> &rhs, gsMatrix<T> &x, bool xZero, bool reverse) const
    {
        const size_t nP = m_patches.size();
        if (nP == 0)
            return;
#ifdef _OPENMP
        if (m_patchesDecoupled && nP > 1 && m_totalInteriorDofs >= 4000)
        {
            const int nThreads = std::min<int>(omp_get_max_threads(), static_cast<int>(nP));
#pragma omp parallel for schedule(static) num_threads(nThreads)
            for (index_t p = 0; p < static_cast<index_t>(nP); ++p)
                relaxPatch(static_cast<size_t>(p), rhs, x, xZero);
            return;
        }
#endif
        if (!reverse)
        {
            for (size_t p = 0; p < nP; ++p)
                relaxPatch(p, rhs, x, xZero && (m_patchesDecoupled || p == 0));
        }
        else
        {
            for (size_t p = nP; p-- > 0;)
                relaxPatch(p, rhs, x, xZero && (m_patchesDecoupled || p + 1 == nP));
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
        const size_t nP = m_patches.size();

#ifdef _OPENMP
        const bool useOmp = (nP > 1 && m_totalInteriorDofs >= 4000);
        const int nThreads = useOmp ? std::min<int>(omp_get_max_threads(), static_cast<int>(nP)) : 1;
#pragma omp parallel for schedule(static) if (useOmp) num_threads(nThreads)
#endif
        for (index_t p = 0; p < static_cast<index_t>(nP); ++p)
        {
#ifdef _OPENMP
            Workspace &ws = m_workspaces[static_cast<size_t>(omp_get_thread_num())];
#else
            Workspace &ws = m_workspaces[0];
#endif
            const PatchS00 &scms = m_patches[static_cast<size_t>(p)];
            if (scms.scale == T(0))
                continue;

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

                applyUnscaledPatchS00(scms, ws.R_mat, ws.E_mat, ws);
                internal::scatterPatchCorrection(scms.grid, ws.E_mat, w, col, x);
            }
        }

        if (m_useL2OrthSkeleton)
            applyL2OrthSkeletonAdditive(r, x, T(0.5));
        else
            m_skeleton.applyAdditive(r, x, T(0.5));
    }

    // ------------------------------------------------------------------------
    // Exact L2(Omega_hat)-Orthogonal Skeleton (Theorem 5 / Eq. (60))
    // ------------------------------------------------------------------------
    void initL2OrthSkeleton(const gsMultiPatch<T> &mp,
                            const std::vector<gsArgyrisEmbedding<T>> &argBasis,
                            const gsDofMapper &mapper,
                            std::vector<index_t> &globToLoc)
    {
        (void)mp;
        const index_t nFree = m_mat.rows();
        std::vector<bool> isPatchInt(static_cast<size_t>(nFree), false);
        for (const auto &scms : m_patches)
            for (index_t dof : scms.grid.dofs)
                if (dof >= 0)
                    isPatchInt[static_cast<size_t>(dof)] = true;

        std::vector<index_t> evColOfDof(static_cast<size_t>(nFree), -1);
        index_t nEV = 0;
        gsSparseEntries<T> seIloc;
        seIloc.reserve(nFree * 2);
        for (index_t g = 0; g < nFree; ++g)
        {
            if (!isPatchInt[static_cast<size_t>(g)])
            {
                evColOfDof[static_cast<size_t>(g)] = nEV;
                seIloc.add(g, nEV, T(1));
                ++nEV;
            }
        }

        index_t nextCol = nEV;
        std::vector<std::vector<index_t>> patchSkelCols(m_patches.size());
        m_l2Patches.resize(m_patches.size());

        for (size_t idx = 0; idx < m_patches.size(); ++idx)
        {
            const PatchS00 &scms = m_patches[idx];
            const size_t p = scms.patchIdx;
            const UnivariateS0 &u1 = *scms.u1;
            const UnivariateS0 &u2 = *scms.u2;
            const index_t m1 = scms.grid.m1;
            const index_t nInt = argBasis[p].sizes[0];
            const index_t nLocal = argBasis[p].matrix.cols();
            const index_t nFull = argBasis[p].matrix.rows();

            // Store free-DOF columns of patch embedding E^{(p)} for implicit O(p N_h) apply
            PatchL2OrthData &l2p = m_l2Patches[idx];
            l2p.scmsIdx = idx;
            gsSparseEntries<T> seE;
            seE.reserve(argBasis[p].matrix.nonZeros());
            std::vector<index_t> seenEV;
            for (index_t loc = 0; loc < nLocal; ++loc)
            {
                if (!mapper.is_free(loc, p))
                    continue;
                const index_t g = mapper.index(loc, p);
                const index_t colIdx = static_cast<index_t>(l2p.freeGlobDofs.size());
                l2p.freeGlobDofs.push_back(g);
                for (typename gsSparseMatrix<T>::InnerIterator it(argBasis[p].matrix, loc); it; ++it)
                    seE.add(it.row(), colIdx, it.value());

                if (loc >= nInt && evColOfDof[static_cast<size_t>(g)] >= 0)
                    seenEV.push_back(evColOfDof[static_cast<size_t>(g)]);
            }
            l2p.E_free.resize(nFull, static_cast<index_t>(l2p.freeGlobDofs.size()));
            l2p.E_free.setFrom(seE);
            l2p.E_free.makeCompressed();

            std::sort(seenEV.begin(), seenEV.end());
            seenEV.erase(std::unique(seenEV.begin(), seenEV.end()), seenEV.end());
            patchSkelCols[idx] = std::move(seenEV);

            // Add patch-local interior orthogonal complement columns (0,perp), (perp,0), (perp,perp)
            if (u1.nPerp > 0 || u2.nPerp > 0)
            {
                const gsSparseMatrix<T> P0_1 = u1.P0.toSparseMatrix();
                const gsSparseMatrix<T> P0_2 = u2.P0.toSparseMatrix();

                auto addTensorCols = [&](const gsSparseMatrix<T> &A1, const gsSparseMatrix<T> &A2) {
                    if (A1.cols() == 0 || A2.cols() == 0)
                        return;
                    for (index_t c2 = 0; c2 < A2.cols(); ++c2)
                    {
                        for (index_t c1 = 0; c1 < A1.cols(); ++c1)
                        {
                            bool hasNnz = false;
                            const index_t col = nextCol;
                            for (typename gsSparseMatrix<T>::InnerIterator it2(A2, c2); it2; ++it2)
                            {
                                for (typename gsSparseMatrix<T>::InnerIterator it1(A1, c1); it1; ++it1)
                                {
                                    const index_t g = scms.grid.dofs[it1.row() + m1 * it2.row()];
                                    if (g >= 0)
                                    {
                                        seIloc.add(g, col, it1.value() * it2.value());
                                        hasNnz = true;
                                    }
                                }
                            }
                            if (hasNnz)
                            {
                                patchSkelCols[idx].push_back(col);
                                ++nextCol;
                            }
                        }
                    }
                };

                addTensorCols(P0_1, u2.Pperp);
                addTensorCols(u1.Pperp, P0_2);
                addTensorCols(u1.Pperp, u2.Pperp);
            }
        }

        m_l2Nskel = nextCol;
        if (m_l2Nskel == 0)
            return;

        m_Iloc.resize(nFree, m_l2Nskel);
        m_Iloc.setFrom(seIloc);
        m_Iloc.makeCompressed();
        m_IlocT = m_Iloc.transpose();
        m_IlocT.makeCompressed();

        // Exact Galerkin assembly of A_skel^perp = (I_skel^perp)^T A_h I_skel^perp
        const gsSparseMatrix<T> AI_loc = m_mat * m_Iloc;
        const gsSparseMatrix<T> A_loc = m_IlocT * AI_loc;

        gsSparseEntries<T> seSkel;
        seSkel.reserve(A_loc.nonZeros() * 4);
        for (index_t c = 0; c < A_loc.cols(); ++c)
            for (typename gsSparseMatrix<T>::InnerIterator it(A_loc, c); it; ++it)
                seSkel.add(it.row(), c, it.value());

        for (size_t idx = 0; idx < m_patches.size(); ++idx)
        {
            const PatchS00 &scms = m_patches[idx];
            const UnivariateS0 &u1 = *scms.u1;
            const UnivariateS0 &u2 = *scms.u2;
            const std::vector<index_t> &skCols = patchSkelCols[idx];
            const index_t kp = static_cast<index_t>(skCols.size());
            if (kp == 0 || u1.n0 <= 0 || u2.n0 <= 0)
                continue;

            const PatchL2OrthData &l2p = m_l2Patches[idx];
            const index_t m1 = scms.grid.m1;
            const index_t m2 = scms.grid.m2;
            const index_t nPatch = m1 * m2;

            // Map global DOF -> local column(s) of E_free on this patch
            std::vector<std::pair<index_t, index_t>> dofCols;
            dofCols.reserve(l2p.freeGlobDofs.size());
            for (size_t lc = 0; lc < l2p.freeGlobDofs.size(); ++lc)
                dofCols.emplace_back(l2p.freeGlobDofs[lc], static_cast<index_t>(lc));
            std::sort(dofCols.begin(), dofCols.end());

            // Build C_p (nPatch x kp): patch interior L2-projection correction -Q_00(I_loc[:, skCols])
            gsMatrix<T> C_p(nPatch, kp);
            gsMatrix<T> V_hat(u1.n, u2.n), tmp(u1.n0, u2.n), Y(u1.n0, u2.n0), U0(m1, u2.n0), W(m1, m2);
            for (index_t j = 0; j < kp; ++j)
            {
                const index_t s = skCols[j];
                V_hat.setZero();
                T *vPtr = V_hat.data();
                for (typename gsSparseMatrix<T>::InnerIterator it(m_Iloc, s); it; ++it)
                {
                    const index_t g = it.row();
                    const T coeff = it.value();
                    auto rangeIt = std::lower_bound(
                        dofCols.begin(), dofCols.end(), std::make_pair(g, index_t(-1)));
                    for (; rangeIt != dofCols.end() && rangeIt->first == g; ++rangeIt)
                    {
                        const index_t lc = rangeIt->second;
                        for (typename gsSparseMatrix<T>::InnerIterator eit(l2p.E_free, lc); eit; ++eit)
                            vPtr[eit.row()] += coeff * eit.value();
                    }
                }

                tmp.noalias() = u1.Mcross * V_hat;
                Y.noalias() = tmp * u2.McrossT;
                u1.cholM00.solveLeftInPlace(Y);
                u2.cholM00.solveRightInPlace(Y);
                u1.P0.applyLeft(Y, U0);
                u2.P0.applyTransposeRight(U0, W);

                T *cCol = C_p.data() + j * nPatch;
                const T *wPtr = W.data();
                for (index_t i = 0; i < nPatch; ++i)
                    cCol[i] = (scms.grid.dofs[i] >= 0) ? (-wPtr[i]) : T(0);
            }

            // Extract B_p = AI_loc[grid.dofs, skCols] (nPatch x kp)
            for (index_t i = 0; i < nPatch; ++i)
                if (scms.grid.dofs[i] >= 0)
                    globToLoc[scms.grid.dofs[i]] = i;

            gsMatrix<T> B_p(nPatch, kp);
            B_p.setZero();
            for (index_t j = 0; j < kp; ++j)
            {
                const index_t s = skCols[j];
                for (typename gsSparseMatrix<T>::InnerIterator it(AI_loc, s); it; ++it)
                {
                    const index_t locRow = globToLoc[it.row()];
                    if (locRow >= 0)
                        B_p(locRow, j) = it.value();
                }
            }
            for (index_t dof : scms.grid.dofs)
                if (dof >= 0)
                    globToLoc[dof] = -1;

            const gsSparseMatrix<T> A_pp =
                internal::extractSymmetricSubmatrix(m_mat, scms.grid.dofs, globToLoc, 25);
            const gsMatrix<T> AC_p = A_pp * C_p;
            const gsMatrix<T> B_t_C = B_p.transpose() * C_p;
            gsMatrix<T> Delta_p = B_t_C + B_t_C.transpose();
            Delta_p.noalias() += C_p.transpose() * AC_p;

            for (index_t c = 0; c < kp; ++c)
                for (index_t r = 0; r < kp; ++r)
                {
                    const T val = Delta_p(r, c);
                    if (std::abs(val) > T(1e-16))
                        seSkel.add(skCols[r], skCols[c], val);
                }
        }

        gsSparseMatrix<T> A_skel_perp(m_l2Nskel, m_l2Nskel);
        A_skel_perp.setFrom(seSkel);
        m_l2SkelSolver.reset(new EigenLDLT());
        m_l2SkelSolver->compute(A_skel_perp);
        if (m_l2SkelSolver->info() != gsEigen::Success)
        {
            gsWarn << "gsAsG1ScmsSmootherOp: L2-orthogonal skeleton Cholesky failed (size="
                   << m_l2Nskel << "); disabling skeleton block.\n";
            m_l2SkelSolver.reset();
        }
    }

    /// Implicit O(p N_h) application of y = (I_skel^perp)^T * r.
    void applyL2OrthPerpT(const gsMatrix<T> &r, gsMatrix<T> &ySkel) const
    {
        const index_t nFree = r.rows();
        const index_t ncols = r.cols();
        m_l2GlobBuf = r;

        for (const PatchL2OrthData &l2p : m_l2Patches)
        {
            const PatchS00 &scms = m_patches[l2p.scmsIdx];
            const UnivariateS0 &u1 = *scms.u1;
            const UnivariateS0 &u2 = *scms.u2;
            if (u1.n0 <= 0 || u2.n0 <= 0)
                continue;

            const index_t m1 = scms.grid.m1;
            const index_t m2 = scms.grid.m2;
            const index_t nPatch = m1 * m2;
            gsMatrix<T> R_int(m1, m2), T0, Z, tmp(u1.n, u2.n0), W_hat(u1.n, u2.n);

            for (index_t col = 0; col < ncols; ++col)
            {
                const T *rCol = r.data() + col * nFree;
                T *rmPtr = R_int.data();
                for (index_t i = 0; i < nPatch; ++i)
                {
                    const index_t g = scms.grid.dofs[i];
                    rmPtr[i] = (g >= 0) ? rCol[g] : T(0);
                }

                u1.P0.applyTransposeLeft(R_int, T0);
                u2.P0.applyRight(T0, Z);
                u1.cholM00.solveLeftInPlace(Z);
                u2.cholM00.solveRightInPlace(Z);
                tmp.noalias() = u1.McrossT * Z;
                W_hat.noalias() = tmp * u2.Mcross;

                const T *wPtr = W_hat.data();
                T *gCol = m_l2GlobBuf.data() + col * nFree;
                const index_t nFreeLoc = static_cast<index_t>(l2p.freeGlobDofs.size());
                for (index_t lc = 0; lc < nFreeLoc; ++lc)
                {
                    T acc = T(0);
                    for (typename gsSparseMatrix<T>::InnerIterator it(l2p.E_free, lc); it; ++it)
                        acc += it.value() * wPtr[it.row()];
                    gCol[l2p.freeGlobDofs[lc]] -= acc;
                }
            }
        }

        ySkel.resize(m_l2Nskel, ncols);
        ySkel.noalias() = m_IlocT * m_l2GlobBuf;
    }

    /// Implicit O(p N_h) update: x += weight * I_skel^perp * eSkel.
    void applyL2OrthPerpAdd(const gsMatrix<T> &eSkel, gsMatrix<T> &x, T weight) const
    {
        const index_t nFree = x.rows();
        const index_t ncols = eSkel.cols();
        m_l2GlobBuf.resize(nFree, ncols);
        m_l2GlobBuf.noalias() = m_Iloc * eSkel;
        x.noalias() += weight * m_l2GlobBuf;

        for (const PatchL2OrthData &l2p : m_l2Patches)
        {
            const PatchS00 &scms = m_patches[l2p.scmsIdx];
            const UnivariateS0 &u1 = *scms.u1;
            const UnivariateS0 &u2 = *scms.u2;
            if (u1.n0 <= 0 || u2.n0 <= 0)
                continue;

            const index_t m1 = scms.grid.m1;
            const index_t m2 = scms.grid.m2;
            gsMatrix<T> V_hat(u1.n, u2.n), tmp(u1.n0, u2.n), Y(u1.n0, u2.n0), U0(m1, u2.n0), W(m1, m2);

            for (index_t col = 0; col < ncols; ++col)
            {
                V_hat.setZero();
                T *vPtr = V_hat.data();
                const T *uCol = m_l2GlobBuf.data() + col * nFree;
                const index_t nFreeLoc = static_cast<index_t>(l2p.freeGlobDofs.size());
                for (index_t lc = 0; lc < nFreeLoc; ++lc)
                {
                    const T coeff = uCol[l2p.freeGlobDofs[lc]];
                    if (coeff == T(0))
                        continue;
                    for (typename gsSparseMatrix<T>::InnerIterator it(l2p.E_free, lc); it; ++it)
                        vPtr[it.row()] += coeff * it.value();
                }

                tmp.noalias() = u1.Mcross * V_hat;
                Y.noalias() = tmp * u2.McrossT;
                u1.cholM00.solveLeftInPlace(Y);
                u2.cholM00.solveRightInPlace(Y);
                u1.P0.applyLeft(Y, U0);
                u2.P0.applyTransposeRight(U0, W);

                internal::scatterPatchCorrection(scms.grid, W, -weight, col, x);
            }
        }
    }

    void relaxL2OrthSkeleton(const gsMatrix<T> &rhs, gsMatrix<T> &x, bool xIsZero) const
    {
        if (!m_l2SkelSolver || m_l2Nskel == 0)
            return;
        const gsMatrix<T> *rPtr = &rhs;
        if (!xIsZero)
        {
            m_resBuf.resize(rhs.rows(), rhs.cols());
            m_resBuf.noalias() = rhs - m_mat * x;
            rPtr = &m_resBuf;
        }
        applyL2OrthPerpT(*rPtr, m_l2Rsub);
        m_l2Esub = m_l2SkelSolver->solve(m_l2Rsub);
        applyL2OrthPerpAdd(m_l2Esub, x, T(1));
    }

    void applyL2OrthSkeletonAdditive(const gsMatrix<T> &r, gsMatrix<T> &x, T weight) const
    {
        if (!m_l2SkelSolver || m_l2Nskel == 0)
            return;
        applyL2OrthPerpT(r, m_l2Rsub);
        m_l2Esub = m_l2SkelSolver->solve(m_l2Rsub);
        applyL2OrthPerpAdd(m_l2Esub, x, weight);
    }

    gsSparseMatrix<T> m_mat;
    bool m_additive;
    bool m_usePointGS;
    bool m_useL2OrthSkeleton;
    bool m_patchesDecoupled = true;
    index_t m_totalInteriorDofs = 0;
    internal::PointGaussSeidelSmoother<T> m_pgs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
    std::vector<PatchS00> m_patches;
    mutable std::vector<Workspace> m_workspaces;
    mutable gsMatrix<T> m_resBuf;

    // Exact L2-orthogonal skeleton state (populated only when m_useL2OrthSkeleton == true)
    index_t m_l2Nskel = 0;
    gsSparseMatrix<T> m_Iloc, m_IlocT;
    std::vector<PatchL2OrthData> m_l2Patches;
    std::unique_ptr<EigenLDLT> m_l2SkelSolver;
    mutable gsMatrix<T> m_l2GlobBuf, m_l2Rsub, m_l2Esub;
};

} // namespace gismo
