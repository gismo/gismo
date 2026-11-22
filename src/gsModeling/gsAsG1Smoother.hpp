/** @file gsAsG1Smoother.hpp

    @brief Specialized O(N) smoothers for AS-G1 multi-patch biharmonic multigrid solvers.

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
#include <vector>

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsSimplePreconditioners.h>
#include <gsModeling/gsAsG1Basis.hpp>

namespace gismo
{

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
    const index_t nFree = mapper.freeSize();
    std::vector<bool> inSkeleton(nFree, false);

    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        const auto &sz = argBasis[p].sizes;
        const index_t nLocal = argBasis[p].matrix.cols();
        const gsTensorBSplineBasis<2, T> &tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
        const index_t n1 = tb.component(0).size();
        const index_t n2 = tb.component(1).size();
        const index_t r1 = (r_p < 0) ? (tb.degree(0) + 1) : std::max<index_t>(2, r_p);
        const index_t r2 = (r_p < 0) ? (tb.degree(1) + 1) : std::max<index_t>(2, r_p);
        const std::vector<index_t> intDofs = getInteriorDofs(tb);

        bool sideActive[5] = {false, includeOuterBoundaryRing, includeOuterBoundaryRing,
                              includeOuterBoundaryRing, includeOuterBoundaryRing};
        if (!includeOuterBoundaryRing)
        {
            for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
            {
                if (it->first().patch  == static_cast<index_t>(p)) sideActive[it->first().m_index]  = true;
                if (it->second().patch == static_cast<index_t>(p)) sideActive[it->second().m_index] = true;
            }
        }

        for (index_t loc = 0; loc < nLocal; ++loc)
        {
            if (!mapper.is_free(loc, p))
                continue;
            const index_t glob = mapper.index(loc, p);
            if (loc < sz[0])
            {
                const index_t tensorIdx = intDofs[loc];
                const index_t i1 = tensorIdx % n1;
                const index_t i2 = tensorIdx / n1;
                const bool rW = (i1 < r1),       rE = (i1 >= n1 - r1);
                const bool rS = (i2 < r2),       rN = (i2 >= n2 - r2);
                if ((rW && sideActive[1]) || (rE && sideActive[2]) ||
                    (rS && sideActive[3]) || (rN && sideActive[4]))
                    inSkeleton[glob] = true;
            }
            else
            {
                inSkeleton[glob] = true;
            }
        }
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
 * @brief Multiplicative Point Gauss-Seidel + Coupled Skeleton Ring-1 Subspace Smoother.
 *
 *  - step(rhs, x):  applies a forward point Gauss-Seidel sweep followed by an
 *                   exact sparse Cholesky correction on the O(sqrt(N)) skeleton
 *                   block (free edge/vertex DOFs + first 2 interior B-spline layers).
 *  - stepT(rhs, x): applies the skeleton Cholesky correction followed by a
 *                   reverse point Gauss-Seidel sweep, making the pre/post
 *                   smoothing pair in gsMultiGridOp strictly A-symmetric.
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
        : m_mat(mat), m_skeletonDofs(skeletonDofs)
    {
        m_gs = makeGaussSeidelOp(m_mat);

        const index_t nSub = static_cast<index_t>(m_skeletonDofs.size());
        if (nSub > 0)
        {
            gsSparseEntries<T> se;
            se.reserve(nSub);
            for (index_t i = 0; i < nSub; ++i)
                se.add(m_skeletonDofs[i], i, T(1));
            gsSparseMatrix<T> P(m_mat.rows(), nSub);
            P.setFrom(se);

            // Pre-extract m_subRows = P^T * A (nSub x N) so residual restriction costs O(sqrt(N))
            m_subRows = P.transpose() * m_mat;
            gsSparseMatrix<T> A_sub = m_subRows * P;
            m_subSolver = makeSparseCholeskySolver(A_sub);
        }
    }

    static uPtr make(const gsSparseMatrix<T> &mat,
                     const std::vector<index_t> &skeletonDofs)
    {
        return uPtr(new gsAsG1SkeletonSmootherOp(mat, skeletonDofs));
    }

    void relaxSkeleton(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        if (!m_subSolver)
            return;

        const index_t nSub = static_cast<index_t>(m_skeletonDofs.size());
        gsMatrix<T> r_sub = -(m_subRows * x);
        for (index_t i = 0; i < nSub; ++i)
            r_sub(i, 0) += rhs(m_skeletonDofs[i], 0);

        gsMatrix<T> e_sub;
        m_subSolver->apply(r_sub, e_sub);

        for (index_t i = 0; i < nSub; ++i)
            x(m_skeletonDofs[i], 0) += e_sub(i, 0);
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        m_gs->step(rhs, x);
        relaxSkeleton(rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        relaxSkeleton(rhs, x);
        m_gs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override
    {
        return makeMatrixOp(m_mat);
    }

    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    std::vector<index_t> m_skeletonDofs;
    gsSparseMatrix<T> m_subRows;
    typename gsPreconditionerOp<T>::uPtr m_gs;
    typename gsLinearOperator<T>::uPtr m_subSolver;
};

/**
 * @brief Build overlapping/paired 1D patch-interior strip blocks (u-lines and v-lines)
 *        followed by the Coupled Ring-1 skeleton block for alternating line-block smoothing.
 *
 * Grouping two adjacent interior B-spline layers per strip handles both high aspect-ratio
 * patch anisotropy and C^1 double-knot B-spline pairs, while the trailing Ring-1 block
 * resolves cross-interface/boundary bending locking.
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
        const index_t nInt = argBasis[p].sizes[0];
        if (nInt <= 0)
            continue;

        const gsTensorBSplineBasis<2, T> &tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
        const index_t n1 = tb.component(0).size();
        const index_t n2 = tb.component(1).size();
        const index_t m1 = n1 - 4;
        const index_t m2 = n2 - 4;
        if (m1 <= 0 || m2 <= 0)
            continue;

        const std::vector<index_t> intDofs = getInteriorDofs(tb);
        std::vector<index_t> intGrid(m1 * m2, -1);
        for (index_t loc = 0; loc < nInt; ++loc)
        {
            if (!mapper.is_free(loc, p))
                continue;
            const index_t tensorIdx = intDofs[loc];
            const index_t i1 = tensorIdx % n1;
            const index_t i2 = tensorIdx / n1;
            const index_t l1 = i1 - 2;
            const index_t l2 = i2 - 2;
            if (l1 >= 0 && l1 < m1 && l2 >= 0 && l2 < m2)
                intGrid[l1 + m1 * l2] = mapper.index(loc, p);
        }

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
                    const index_t g = intGrid[ii + m1 * jj];
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
                    const index_t g = intGrid[ii + m1 * jj];
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

/**
 * @brief Multiplicative Block Gauss-Seidel Smoother with Exact Block Cholesky Solves.
 *
 *  - step(rhs, x):  relaxes blocks 0 .. K-1 in forward order.
 *  - stepT(rhs, x): relaxes blocks K-1 .. 0 in reverse order, making pre/post
 *                   smoothing in gsMultiGridOp strictly A-symmetric.
 */
template <class T>
class gsAsG1BlockSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::shared_ptr<gsAsG1BlockSmootherOp> Ptr;
    typedef memory::unique_ptr<gsAsG1BlockSmootherOp> uPtr;
    typedef gsPreconditionerOp<T> Base;
    typedef typename gsSparseSolver<T>::SimplicialLDLT SubSolver;

    gsAsG1BlockSmootherOp(const gsSparseMatrix<T> &mat,
                          const std::vector<std::vector<index_t>> &blocks)
        : m_mat(mat), m_matRow(mat)
    {
        const index_t nFree = m_mat.rows();
        std::vector<index_t> globToLoc(nFree, -1);

        m_blocks.reserve(blocks.size());
        m_solvers.reserve(blocks.size());

        for (size_t b = 0; b < blocks.size(); ++b)
        {
            if (blocks[b].empty())
                continue;

            std::vector<index_t> dofs;
            dofs.reserve(blocks[b].size());
            for (index_t idx : blocks[b])
            {
                if (idx >= 0 && idx < nFree && globToLoc[idx] < 0)
                {
                    globToLoc[idx] = static_cast<index_t>(dofs.size());
                    dofs.push_back(idx);
                }
            }

            const index_t k = static_cast<index_t>(dofs.size());
            if (k == 0)
                continue;

            gsSparseEntries<T> se;
            se.reserve(k * 16);
            for (index_t i = 0; i < k; ++i)
            {
                const index_t row = dofs[i];
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                {
                    const index_t j = globToLoc[it.col()];
                    if (j >= 0)
                        se.add(i, j, it.value());
                }
            }

            for (index_t idx : dofs)
                globToLoc[idx] = -1;

            gsSparseMatrix<T> A_sub(k, k);
            A_sub.setFrom(se);

            m_solvers.emplace_back(new SubSolver());
            m_solvers.back()->compute(A_sub);
            m_blocks.push_back(std::move(dofs));
        }
    }

    static uPtr make(const gsSparseMatrix<T> &mat,
                     const std::vector<std::vector<index_t>> &blocks)
    {
        return uPtr(new gsAsG1BlockSmootherOp(mat, blocks));
    }

    void relaxBlock(size_t b, const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        const std::vector<index_t> &dofs = m_blocks[b];
        const index_t k = static_cast<index_t>(dofs.size());
        const index_t ncols = rhs.cols();

        gsMatrix<T> r_sub(k, ncols);
        if (ncols == 1)
        {
            for (index_t i = 0; i < k; ++i)
            {
                const index_t row = dofs[i];
                T acc = rhs(row, 0);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    acc -= it.value() * x(it.col(), 0);
                r_sub(i, 0) = acc;
            }
        }
        else
        {
            for (index_t i = 0; i < k; ++i)
            {
                const index_t row = dofs[i];
                r_sub.row(i) = rhs.row(row);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    r_sub.row(i) -= it.value() * x.row(it.col());
            }
        }

        const gsMatrix<T> e_sub = m_solvers[b]->solve(r_sub);
        for (index_t i = 0; i < k; ++i)
            x.row(dofs[i]) += e_sub.row(i);
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        for (size_t b = 0; b < m_blocks.size(); ++b)
            relaxBlock(b, rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        for (size_t b = m_blocks.size(); b-- > 0;)
            relaxBlock(b, rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override
    {
        return makeMatrixOp(m_mat);
    }

    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    gsSparseMatrix<T, RowMajor> m_matRow;
    std::vector<std::vector<index_t>> m_blocks;
    std::vector<std::unique_ptr<SubSolver>> m_solvers;
};

/**
 * @brief Patch-Interior Fast Diagonalization (FD) + Coupled Ring-1 Skeleton Smoother
 *        for p-robust and anisotropic biharmonic multigrid smoothing.
 *
 * On each patch interior (m1 x m2 tensor B-spline subspace [2..n-3]^2), solves the
 * 1D clamped biharmonic-mass generalized eigenproblems (B_d, M_d) in u and v to form
 * an O(m1 m2 (m1 + m2)) Kronecker FD solver calibrated by the patch physical metric
 * and power-iteration Rayleigh damping. Combined symmetrically with Point GS and
 * Coupled Ring-1 skeleton Cholesky:
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
    typedef typename gsSparseSolver<T>::SimplicialLDLT SubSolver;

    gsAsG1FastDiagSmootherOp(const gsSparseMatrix<T> &mat,
                             const gsMultiPatch<T> &mp,
                             const gsMultiBasis<T> &dbasis,
                             const std::vector<gsArgyrisEmbedding<T>> &argBasis,
                             const gsDofMapper &mapper,
                             T beta = T(0),
                             bool usePointGS = true,
                             bool includeOuterBoundaryRing = true)
        : m_mat(mat), m_matRow(mat), m_usePointGS(usePointGS)
    {
        if (m_usePointGS)
            m_gs = makeGaussSeidelOp(m_mat);

        const index_t nFree = m_mat.rows();
        std::vector<index_t> globToLoc(nFree, -1);

        // 1. Setup Coupled Ring-1 skeleton exact Cholesky block
        m_skeletonDofs = collectCoupledRing1Dofs(
            mp, dbasis, argBasis, mapper, includeOuterBoundaryRing);
        if (!m_skeletonDofs.empty())
        {
            const index_t nSkel = static_cast<index_t>(m_skeletonDofs.size());
            for (index_t i = 0; i < nSkel; ++i)
                globToLoc[m_skeletonDofs[i]] = i;

            gsSparseEntries<T> se;
            se.reserve(nSkel * 24);
            for (index_t i = 0; i < nSkel; ++i)
            {
                const index_t row = m_skeletonDofs[i];
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                {
                    const index_t j = globToLoc[it.col()];
                    if (j >= 0)
                        se.add(i, j, it.value());
                }
            }
            for (index_t idx : m_skeletonDofs)
                globToLoc[idx] = -1;

            gsSparseMatrix<T> A_skel(nSkel, nSkel);
            A_skel.setFrom(se);
            m_skelSolver.reset(new SubSolver());
            m_skelSolver->compute(A_skel);
        }

        // 2. Setup per-patch Fast Diagonalization operators on interior tensor subspaces
        m_patches.reserve(mp.nPatches());
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const index_t nInt = argBasis[p].sizes[0];
            if (nInt <= 0)
                continue;

            const gsTensorBSplineBasis<2, T> &tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
            const auto &b1 = tb.component(0);
            const auto &b2 = tb.component(1);
            const index_t n1 = b1.size();
            const index_t n2 = b2.size();
            const index_t m1 = n1 - 4;
            const index_t m2 = n2 - 4;
            if (m1 <= 0 || m2 <= 0)
                continue;

            PatchFD fd;
            fd.m1 = m1;
            fd.m2 = m2;
            fd.dofs.assign(m1 * m2, -1);

            const std::vector<index_t> intDofs = getInteriorDofs(tb);
            bool hasFree = false;
            for (index_t loc = 0; loc < nInt; ++loc)
            {
                if (!mapper.is_free(loc, p))
                    continue;
                const index_t tensorIdx = intDofs[loc];
                const index_t l1 = (tensorIdx % n1) - 2;
                const index_t l2 = (tensorIdx / n1) - 2;
                if (l1 >= 0 && l1 < m1 && l2 >= 0 && l2 < m2)
                {
                    fd.dofs[l1 + m1 * l2] = mapper.index(loc, p);
                    hasFree = true;
                }
            }
            if (!hasFree)
                continue;

            // Assemble 1D interior mass and biharmonic (2nd-derivative) matrices
            gsMatrix<T> M1, B1, M2, B2;
            assemble1DInteriorMassAndBiharmonic(b1, M1, B1);
            assemble1DInteriorMassAndBiharmonic(b2, M2, B2);

            typename gsMatrix<T>::GenSelfAdjEigenSolver ges1(B1, M1, gsEigen::ComputeEigenvectors);
            typename gsMatrix<T>::GenSelfAdjEigenSolver ges2(B2, M2, gsEigen::ComputeEigenvectors);
            fd.Q1 = ges1.eigenvectors();
            fd.Q2 = ges2.eigenvectors();
            const gsVector<T> lam1 = ges1.eigenvalues().cwiseMax(T(0));
            const gsVector<T> lam2 = ges2.eigenvalues().cwiseMax(T(0));

            // Estimate average physical metric coefficients on patch p
            T c1 = T(0), c2 = T(0), detMean = T(0);
            const gsMatrix<T> pr = mp.patch(p).parameterRange();
            const int nSample = 3;
            gsMatrix<T> pt(2, 1);
            for (int sy = 0; sy < nSample; ++sy)
            {
                pt(1, 0) = pr(1, 0) + (T(sy) + T(0.5)) / T(nSample) * (pr(1, 1) - pr(1, 0));
                for (int sx = 0; sx < nSample; ++sx)
                {
                    pt(0, 0) = pr(0, 0) + (T(sx) + T(0.5)) / T(nSample) * (pr(0, 1) - pr(0, 0));
                    const gsMatrix<T> J = mp.patch(p).jacobian(pt);
                    const gsMatrix<T> g = J.transpose() * J;
                    const T detG = std::max(g.determinant(), T(1e-30));
                    const T absDetJ = std::sqrt(detG);
                    const T ginv11 = g(1, 1) / detG;
                    const T ginv22 = g(0, 0) / detG;
                    c1 += absDetJ * ginv11 * ginv11;
                    c2 += absDetJ * ginv22 * ginv22;
                    detMean += absDetJ;
                }
            }
            c1 /= T(nSample * nSample);
            c2 /= T(nSample * nSample);
            detMean /= T(nSample * nSample);
            const T c0 = beta * detMean;

            fd.invDiag.resize(m1, m2);
            for (index_t b = 0; b < m2; ++b)
            {
                const T s2 = std::sqrt(c2 * lam2(b));
                for (index_t a = 0; a < m1; ++a)
                {
                    const T s1 = std::sqrt(c1 * lam1(a));
                    const T denom = (s1 + s2) * (s1 + s2) + c0;
                    fd.invDiag(a, b) = (denom > T(0)) ? (T(1) / denom) : T(0);
                }
            }

            // Extract sparse patch-interior diagonal block A_pp and calibrate damping via power method
            const index_t nPatch = m1 * m2;
            for (index_t i = 0; i < nPatch; ++i)
                if (fd.dofs[i] >= 0)
                    globToLoc[fd.dofs[i]] = i;

            gsSparseEntries<T> se;
            se.reserve(nPatch * 25);
            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = fd.dofs[i];
                if (row < 0)
                    continue;
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                {
                    const index_t j = globToLoc[it.col()];
                    if (j >= 0)
                        se.add(i, j, it.value());
                }
            }
            for (index_t idx : fd.dofs)
                if (idx >= 0)
                    globToLoc[idx] = -1;

            gsSparseMatrix<T> A_pp(nPatch, nPatch);
            A_pp.setFrom(se);

            // 12 steps of power iteration in the A_pp-inner product for lambda_max(A_fd^{-1} A_pp)
            gsVector<T> v(nPatch);
            for (index_t i = 0; i < nPatch; ++i)
                v(i) = (fd.dofs[i] >= 0)
                    ? (std::sin(T(i + 1) * T(1.37)) + T(0.5) * std::cos(T(i + 1) * T(0.23)))
                    : T(0);
            T vNorm = v.norm();
            if (vNorm > T(0))
                v /= vNorm;

            T lamMax = T(1);
            gsMatrix<T> R_mat(m1, m2), Z_mat(m1, m2), E_mat(m1, m2);
            gsVector<T> z(nPatch);
            for (int iter = 0; iter < 12; ++iter)
            {
                const gsVector<T> Av = A_pp * v;
                const T vAv = v.dot(Av);
                if (vAv <= T(1e-30))
                    break;

                std::copy(Av.data(), Av.data() + nPatch, R_mat.data());
                Z_mat.noalias() = fd.Q1.transpose() * R_mat * fd.Q2;
                Z_mat.array() *= fd.invDiag.array();
                E_mat.noalias() = fd.Q1 * Z_mat * fd.Q2.transpose();
                std::copy(E_mat.data(), E_mat.data() + nPatch, z.data());
                for (index_t i = 0; i < nPatch; ++i)
                    if (fd.dofs[i] < 0)
                        z(i) = T(0);

                const T zAv = z.dot(Av);
                if (zAv > T(0))
                    lamMax = zAv / vAv;

                const gsVector<T> Az = A_pp * z;
                const T zAz = z.dot(Az);
                if (zAz <= T(1e-30))
                    break;
                v = z / std::sqrt(zAz);
            }

            if (lamMax > T(1e-12))
                fd.invDiag *= (T(0.90) / lamMax);

            m_patches.push_back(std::move(fd));
        }
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

    void relaxPatch(size_t pIdx, const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        const PatchFD &fd = m_patches[pIdx];
        const index_t m1 = fd.m1;
        const index_t m2 = fd.m2;
        const index_t nPatch = m1 * m2;
        const index_t ncols = rhs.cols();

        gsMatrix<T> R_mat(m1, m2), Z_mat(m1, m2), E_mat(m1, m2);
        for (index_t col = 0; col < ncols; ++col)
        {
            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = fd.dofs[i];
                if (row < 0)
                {
                    R_mat.data()[i] = T(0);
                    continue;
                }
                T acc = rhs(row, col);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    acc -= it.value() * x(it.col(), col);
                R_mat.data()[i] = acc;
            }

            Z_mat.noalias() = fd.Q1.transpose() * R_mat * fd.Q2;
            Z_mat.array() *= fd.invDiag.array();
            E_mat.noalias() = fd.Q1 * Z_mat * fd.Q2.transpose();

            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = fd.dofs[i];
                if (row >= 0)
                    x(row, col) += E_mat.data()[i];
            }
        }
    }

    void relaxSkeleton(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        if (!m_skelSolver)
            return;

        const index_t nSub = static_cast<index_t>(m_skeletonDofs.size());
        const index_t ncols = rhs.cols();
        gsMatrix<T> r_sub(nSub, ncols);
        if (ncols == 1)
        {
            for (index_t i = 0; i < nSub; ++i)
            {
                const index_t row = m_skeletonDofs[i];
                T acc = rhs(row, 0);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    acc -= it.value() * x(it.col(), 0);
                r_sub(i, 0) = acc;
            }
        }
        else
        {
            for (index_t i = 0; i < nSub; ++i)
            {
                const index_t row = m_skeletonDofs[i];
                r_sub.row(i) = rhs.row(row);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    r_sub.row(i) -= it.value() * x.row(it.col());
            }
        }

        const gsMatrix<T> e_sub = m_skelSolver->solve(r_sub);
        for (index_t i = 0; i < nSub; ++i)
            x.row(m_skeletonDofs[i]) += e_sub.row(i);
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_usePointGS)
            m_gs->step(rhs, x);
        for (size_t p = 0; p < m_patches.size(); ++p)
            relaxPatch(p, rhs, x);
        relaxSkeleton(rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        relaxSkeleton(rhs, x);
        for (size_t p = m_patches.size(); p-- > 0;)
            relaxPatch(p, rhs, x);
        if (m_usePointGS)
            m_gs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override
    {
        return makeMatrixOp(m_mat);
    }

    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    struct PatchFD
    {
        index_t m1, m2;
        std::vector<index_t> dofs;
        gsMatrix<T> Q1;
        gsMatrix<T> Q2;
        gsMatrix<T> invDiag;
    };

    static void assemble1DInteriorMassAndBiharmonic(
        const gsBSplineBasis<T> &b1d,
        gsMatrix<T> &M_int,
        gsMatrix<T> &B_int)
    {
        const index_t n = b1d.size();
        const index_t m = n - 4;
        M_int.setZero(m, m);
        B_int.setZero(m, m);

        std::vector<T> uknots = b1d.knots().unique();
        gsGaussRule<T> quRule(b1d.degree() + 1);
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
                    M_int(ir, ic) += mVal;
                    B_int(ir, ic) += bVal;
                }
            }
        }
    }

    gsSparseMatrix<T> m_mat;
    gsSparseMatrix<T, RowMajor> m_matRow;
    bool m_usePointGS;
    typename gsPreconditionerOp<T>::uPtr m_gs;
    std::vector<index_t> m_skeletonDofs;
    std::unique_ptr<SubSolver> m_skelSolver;
    std::vector<PatchFD> m_patches;
};

/**
 * @brief O(N_h) h- and p-robust Subspace-Corrected Mass Smoother (SCMS) + Overlapping
 *        Skeleton Schwarz Smoother for multi-patch AS-G1 biharmonic discretizations (Section 5).
 *
 * On each patch interior S^0_{p,p-2,h}(Omega_hat), decomposes the 1D spline spaces into
 * the generalized Dirichlet subspace S^{D,0} (where |u|_{H^2}^2 <= 2520 h^{-4} ||u||_{L^2}^2,
 * solved in O(p m) via a localized banded basis Q_0 and sparse banded Cholesky of M_{0,0})
 * and its (B_0 + sigma_0 M_0)-orthogonal boundary complement S^{D,1} of dimension <= p - 2
 * (simultaneously diagonalized in O(p^3)). Coupled with an exact sparse Cholesky solve on
 * the enlarged skeleton space \widetilde{V}_h^{skel} of strip width r_p:
 *
 *  - Multiplicative mode (additive = false, default, Remark 5):
 *      step(rhs, x):  parallel patch-interior SCMS -> skeleton Cholesky
 *      stepT(rhs, x): skeleton Cholesky -> parallel patch-interior SCMS  (exact A-adjoint)
 *  - Additive mode (additive = true, Eq. (57) / Theorem 5):
 *      x <- x + 0.5 * ( sum_i I_i (L_i^{int})^{-1} I_i^T r + I_skel A_skel^{-1} I_skel^T r )
 */
template <class T>
class gsAsG1ScmsSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::shared_ptr<gsAsG1ScmsSmootherOp> Ptr;
    typedef memory::unique_ptr<gsAsG1ScmsSmootherOp> uPtr;
    typedef gsPreconditionerOp<T> Base;
    typedef typename gsSparseSolver<T>::SimplicialLDLT SubSolver;

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
                         bool includeOuterBoundaryRing = true)
        : m_mat(mat), m_matRow(mat), m_additive(additive)
    {
        const index_t nFree = m_mat.rows();
        std::vector<index_t> globToLoc(nFree, -1);

        // 1. Setup enlarged skeleton Cholesky solver on \widetilde{V}_h^{skel}
        m_skeletonDofs = collectCoupledStripDofs(
            mp, dbasis, argBasis, mapper, r_p, includeOuterBoundaryRing);
        if (!m_skeletonDofs.empty())
        {
            const index_t nSkel = static_cast<index_t>(m_skeletonDofs.size());
            for (index_t i = 0; i < nSkel; ++i)
                globToLoc[m_skeletonDofs[i]] = i;

            gsSparseEntries<T> se;
            se.reserve(nSkel * 32);
            for (index_t i = 0; i < nSkel; ++i)
            {
                const index_t row = m_skeletonDofs[i];
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                {
                    const index_t j = globToLoc[it.col()];
                    if (j >= 0)
                        se.add(i, j, it.value());
                }
            }
            for (index_t idx : m_skeletonDofs)
                globToLoc[idx] = -1;

            gsSparseMatrix<T> A_skel(nSkel, nSkel);
            A_skel.setFrom(se);
            m_skelSolver.reset(new SubSolver());
            m_skelSolver->compute(A_skel);
        }

        // 2. Setup O(p N_h^{(i)}) Kronecker SCMS on each patch interior S^0_{p,p-2,h}
        m_patches.reserve(mp.nPatches());
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const index_t nInt = argBasis[p].sizes[0];
            if (nInt <= 0)
                continue;

            const gsTensorBSplineBasis<2, T> &tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[p]);
            const auto &b1 = tb.component(0);
            const auto &b2 = tb.component(1);
            const index_t n1 = b1.size();
            const index_t n2 = b2.size();
            const index_t m1 = n1 - 4;
            const index_t m2 = n2 - 4;
            if (m1 <= 0 || m2 <= 0)
                continue;

            PatchSCMS scms;
            scms.m1 = m1;
            scms.m2 = m2;
            scms.dofs.assign(m1 * m2, -1);

            const std::vector<index_t> intDofs = getInteriorDofs(tb);
            bool hasFree = false;
            for (index_t loc = 0; loc < nInt; ++loc)
            {
                if (!mapper.is_free(loc, p))
                    continue;
                const index_t tensorIdx = intDofs[loc];
                const index_t l1 = (tensorIdx % n1) - 2;
                const index_t l2 = (tensorIdx / n1) - 2;
                if (l1 >= 0 && l1 < m1 && l2 >= 0 && l2 < m2)
                {
                    scms.dofs[l1 + m1 * l2] = mapper.index(loc, p);
                    hasFree = true;
                }
            }
            if (!hasFree)
                continue;

            // Sample patch physical geometry metric (and analytical c_bar_G bound)
            T c1 = T(1), c2 = T(1), detMean = T(1), cBarG = T(1);
            samplePatchMetric(mp.patch(p), useMetricWeights, c1, c2, detMean, cBarG);
            scms.c0 = beta * detMean;

            buildUnivariateSCMS(b1, c1, scms.u1);
            buildUnivariateSCMS(b2, c2, scms.u2);

            // Precompute small 1D/2D diagonal scaling tables for L_00, L_10, L_01, L_11
            scms.invSigma00 = T(1) / (scms.u1.sigma0 + scms.u2.sigma0 + scms.c0);

            if (scms.u1.nD1 > 0)
            {
                scms.invDiag10.resize(scms.u1.nD1);
                for (index_t a = 0; a < scms.u1.nD1; ++a)
                {
                    const T denom = scms.u1.b11(a) + (scms.u2.sigma0 + scms.c0) * scms.u1.mu11(a);
                    scms.invDiag10(a) = (denom > T(0)) ? (T(1) / denom) : T(0);
                }
            }
            if (scms.u2.nD1 > 0)
            {
                scms.invDiag01.resize(scms.u2.nD1);
                for (index_t b = 0; b < scms.u2.nD1; ++b)
                {
                    const T denom = scms.u2.b11(b) + (scms.u1.sigma0 + scms.c0) * scms.u2.mu11(b);
                    scms.invDiag01(b) = (denom > T(0)) ? (T(1) / denom) : T(0);
                }
            }
            if (scms.u1.nD1 > 0 && scms.u2.nD1 > 0)
            {
                scms.invDiag11.resize(scms.u1.nD1, scms.u2.nD1);
                for (index_t b = 0; b < scms.u2.nD1; ++b)
                    for (index_t a = 0; a < scms.u1.nD1; ++a)
                    {
                        const T denom = scms.u1.b11(a) * scms.u2.mu11(b)
                                      + scms.u1.mu11(a) * scms.u2.b11(b)
                                      + scms.c0 * scms.u1.mu11(a) * scms.u2.mu11(b);
                        scms.invDiag11(a, b) = (denom > T(0)) ? (T(1) / denom) : T(0);
                    }
            }

            scms.scale = T(1);
            if (calibrateByPowerMethod)
            {
                // Extract sparse patch-interior stiffness block A_pp and calibrate lambda_max
                const index_t nPatch = m1 * m2;
                for (index_t i = 0; i < nPatch; ++i)
                    if (scms.dofs[i] >= 0)
                        globToLoc[scms.dofs[i]] = i;

                gsSparseEntries<T> se;
                se.reserve(nPatch * 25);
                for (index_t i = 0; i < nPatch; ++i)
                {
                    const index_t row = scms.dofs[i];
                    if (row < 0)
                        continue;
                    for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    {
                        const index_t j = globToLoc[it.col()];
                        if (j >= 0)
                            se.add(i, j, it.value());
                    }
                }
                for (index_t idx : scms.dofs)
                    if (idx >= 0)
                        globToLoc[idx] = -1;

                gsSparseMatrix<T> A_pp(nPatch, nPatch);
                A_pp.setFrom(se);

                gsVector<T> v(nPatch);
                for (index_t i = 0; i < nPatch; ++i)
                    v(i) = (scms.dofs[i] >= 0)
                        ? (std::sin(T(i + 1) * T(1.37)) + T(0.5) * std::cos(T(i + 1) * T(0.23)))
                        : T(0);
                const T vNorm = v.norm();
                if (vNorm > T(0))
                    v /= vNorm;

                T lamMax = T(1);
                gsMatrix<T> R_mat(m1, m2), E_mat(m1, m2);
                gsVector<T> z(nPatch);
                for (int iter = 0; iter < 15; ++iter)
                {
                    const gsVector<T> Av = A_pp * v;
                    const T vAv = v.dot(Av);
                    if (vAv <= T(1e-30))
                        break;

                    std::copy(Av.data(), Av.data() + nPatch, R_mat.data());
                    applyUnscaledPatchSCMS(scms, R_mat, E_mat);
                    std::copy(E_mat.data(), E_mat.data() + nPatch, z.data());
                    for (index_t i = 0; i < nPatch; ++i)
                        if (scms.dofs[i] < 0)
                            z(i) = T(0);

                    const T zAv = z.dot(Av);
                    if (zAv > T(0))
                        lamMax = zAv / vAv;

                    const gsVector<T> Az = A_pp * z;
                    const T zAz = z.dot(Az);
                    if (zAz <= T(1e-30))
                        break;
                    v = z / std::sqrt(zAz);
                }
                scms.scale = (lamMax > T(1e-12)) ? (T(0.92) / lamMax) : T(1);
            }
            else
            {
                // Analytical bound (8 * c_bar_G)^{-1} from Eq. (55) of Section 5.2
                scms.scale = T(1) / (T(8) * cBarG);
            }

            m_patches.push_back(std::move(scms));
        }
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
                     bool includeOuterBoundaryRing = true)
    {
        return uPtr(new gsAsG1ScmsSmootherOp(
            mat, mp, dbasis, argBasis, mapper, beta, r_p,
            additive, useMetricWeights, calibrateByPowerMethod, includeOuterBoundaryRing));
    }

    void step(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_additive)
        {
            stepAdditive(rhs, x);
            return;
        }
        for (size_t p = 0; p < m_patches.size(); ++p)
            relaxPatch(p, rhs, x);
        relaxSkeleton(rhs, x);
    }

    void stepT(const gsMatrix<T> &rhs, gsMatrix<T> &x) const override
    {
        if (m_additive)
        {
            stepAdditive(rhs, x);
            return;
        }
        relaxSkeleton(rhs, x);
        for (size_t p = m_patches.size(); p-- > 0;)
            relaxPatch(p, rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override
    {
        return makeMatrixOp(m_mat);
    }

    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    struct UnivariateSCMS
    {
        index_t m = 0;
        index_t nD0 = 0;
        index_t nD1 = 0;
        T sigma0 = T(1);
        gsSparseMatrix<T> Q0;               // m x nD0 sparse block-diagonal basis of S^{D,0}
        std::unique_ptr<SubSolver> solverM00; // banded Cholesky of M_{0,0} = Q_0^T M_0 Q_0
        gsMatrix<T> Q1;                     // m x nD1 basis of S^{D,1} (simultaneously diagonalizing B_0, M_0)
        gsVector<T> b11;                    // diagonal entries c_d * q_a^T B_0 q_a > 0
        gsVector<T> mu11;                   // diagonal entries q_a^T M_0 q_a > 0
    };

    struct PatchSCMS
    {
        index_t m1 = 0, m2 = 0;
        std::vector<index_t> dofs;
        UnivariateSCMS u1, u2;
        T c0 = T(0);
        T invSigma00 = T(1);
        gsVector<T> invDiag10;
        gsVector<T> invDiag01;
        gsMatrix<T> invDiag11;
        T scale = T(1);
    };

    static void samplePatchMetric(const gsGeometry<T> &geo,
                                  bool useMetricWeights,
                                  T &c1, T &c2, T &detMean, T &cBarG)
    {
        const gsMatrix<T> pr = geo.parameterRange();
        const int nSample = 4;
        T sumC1 = T(0), sumC2 = T(0), sumDet = T(0), maxBound = T(0);
        gsMatrix<T> pt(2, 1);
        for (int sy = 0; sy < nSample; ++sy)
        {
            pt(1, 0) = pr(1, 0) + (T(sy) + T(0.5)) / T(nSample) * (pr(1, 1) - pr(1, 0));
            for (int sx = 0; sx < nSample; ++sx)
            {
                pt(0, 0) = pr(0, 0) + (T(sx) + T(0.5)) / T(nSample) * (pr(0, 1) - pr(0, 0));
                const gsMatrix<T> J = geo.jacobian(pt);
                const gsMatrix<T> g = J.transpose() * J;
                const T detG = std::max(g.determinant(), T(1e-30));
                const T absDetJ = std::sqrt(detG);
                const T ginv11 = g(1, 1) / detG;
                const T ginv12 = -g(0, 1) / detG;
                const T ginv22 = g(0, 0) / detG;
                sumC1  += absDetJ * ginv11 * ginv11;
                sumC2  += absDetJ * ginv22 * ginv22;
                sumDet += absDetJ;
                const T frobSq = ginv11 * ginv11 + T(2) * ginv12 * ginv12 + ginv22 * ginv22;
                maxBound = std::max(maxBound, absDetJ * frobSq);
            }
        }
        const T invN = T(1) / T(nSample * nSample);
        detMean = sumDet * invN;
        cBarG   = std::max(maxBound, T(1e-12));
        if (useMetricWeights)
        {
            c1 = std::max(sumC1 * invN, T(1e-12));
            c2 = std::max(sumC2 * invN, T(1e-12));
        }
        else
        {
            c1 = T(1);
            c2 = T(1);
        }
    }

    static void assemble1DSparseInteriorMatrices(
        const gsBSplineBasis<T> &b1d,
        gsSparseMatrix<T> &M0,
        gsSparseMatrix<T> &B0,
        T &h_d)
    {
        const index_t n = b1d.size();
        const index_t m = n - 4;
        const short_t deg = b1d.degree();

        std::vector<T> uknots = b1d.knots().unique();
        const index_t nElem = std::max<index_t>(1, static_cast<index_t>(uknots.size()) - 1);
        h_d = (uknots.back() - uknots.front()) / T(nElem);

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
            const index_t nQ   = quWeights.rows();
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

    static void buildUnivariateSCMS(const gsBSplineBasis<T> &b1d,
                                    T c_d,
                                    UnivariateSCMS &out)
    {
        const short_t deg = b1d.degree();
        const index_t n   = b1d.size();
        const index_t m   = n - 4;
        out.m = m;

        gsSparseMatrix<T> M0, B0;
        T h_d = T(1);
        assemble1DSparseInteriorMatrices(b1d, M0, B0, h_d);

        // Lemma 4: C_inv = 2520 on S^{D,0}_{p,p-2,h}
        out.sigma0 = T(2520) * c_d / (h_d * h_d * h_d * h_d);
        const gsSparseMatrix<T> W0 = c_d * B0 + out.sigma0 * M0;

        const index_t r_ev = (deg - 2) / 2; // number of even derivatives 2j <= p - 2 per end
        if (r_ev <= 0)
        {
            out.nD0 = m;
            out.nD1 = 0;
            gsSparseEntries<T> seI;
            seI.reserve(m);
            for (index_t i = 0; i < m; ++i)
                seI.add(i, i, T(1));
            out.Q0.resize(m, m);
            out.Q0.setFrom(seI);

            out.solverM00.reset(new SubSolver());
            out.solverM00->compute(M0);
            return;
        }

        // Evaluate boundary even derivatives on unit-spaced knots [0, nElem] for high-p stability
        std::vector<T> uknots = b1d.knots().unique();
        const index_t nElem = std::max<index_t>(1, static_cast<index_t>(uknots.size()) - 1);
        const T uMin = uknots.front();

        std::vector<T> scaledKnots;
        scaledKnots.reserve(b1d.knots().size());
        for (typename gsKnotVector<T>::const_iterator it = b1d.knots().begin();
             it != b1d.knots().end(); ++it)
            scaledKnots.push_back((*it - uMin) / h_d);

        gsKnotVector<T> kvScaled(deg, scaledKnots.begin(), scaledKnots.end());
        gsBSplineBasis<T> bScaled(kvScaled);

        gsMatrix<T> endPts(1, 2);
        endPts(0, 0) = T(0);
        endPts(0, 1) = T(nElem);
        std::vector<gsMatrix<T>> evEnd;
        gsMatrix<index_t> actEnd;
        bScaled.evalAllDers_into(endPts, static_cast<int>(2 * r_ev), evEnd);
        bScaled.active_into(endPts, actEnd);

        const index_t kb = deg - 1; // at most p - 1 interior B-splines have nonzero derivs at an end
        gsMatrix<T> C_full;

        if (m >= 2 * kb)
        {
            // Disjoint left and right boundary blocks -> sparse block-diagonal Q_0
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
            const index_t nMid  = m - 2 * kb;
            out.nD0 = 2 * kNull + nMid;

            gsSparseEntries<T> seQ0;
            seQ0.reserve(2 * kb * kNull + nMid);
            for (index_t c = 0; c < kNull; ++c)
                for (index_t r = 0; r < kb; ++r)
                    if (std::abs(V_L(r, r_ev + c)) > T(1e-14))
                        seQ0.add(r, c, V_L(r, r_ev + c));

            for (index_t i = 0; i < nMid; ++i)
                seQ0.add(kb + i, kNull + i, T(1));

            for (index_t c = 0; c < kNull; ++c)
                for (index_t r = 0; r < kb; ++r)
                    if (std::abs(V_R(r, r_ev + c)) > T(1e-14))
                        seQ0.add(m - kb + r, kNull + nMid + c, V_R(r, r_ev + c));

            out.Q0.resize(m, out.nD0);
            out.Q0.setFrom(seQ0);

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
            const gsMatrix<T> Q0_dense = svd.matrixV().rightCols(out.nD0);
            out.Q0 = Q0_dense.sparseView(T(1), T(1e-14));
            C_full = svd.matrixU().leftCols(rk).transpose() * C_raw;
        }

        if (out.nD0 > 0)
        {
            const gsSparseMatrix<T> M00 = out.Q0.transpose() * M0 * out.Q0;
            out.solverM00.reset(new SubSolver());
            out.solverM00->compute(M00);
        }

        // Construct S^{D,1} = (S^{D,0})^{\perp_{W_0}} = span(W_0^{-1} C_full^T)
        if (C_full.rows() == 0)
        {
            out.nD1 = 0;
            return;
        }

        SubSolver solverW0;
        solverW0.compute(W0);
        const gsMatrix<T> C_T = C_full.transpose();
        gsMatrix<T> Q1_tilde = solverW0.solve(C_T);

        // Normalize columns in W_0-norm and W_0-orthonormalize via symmetric eigensolve
        for (index_t k = 0; k < Q1_tilde.cols(); ++k)
        {
            const T wNormSq = Q1_tilde.col(k).dot(W0 * Q1_tilde.col(k));
            if (wNormSq > T(0))
                Q1_tilde.col(k) /= std::sqrt(wNormSq);
        }

        const gsMatrix<T> G_W = Q1_tilde.transpose() * (W0 * Q1_tilde);
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
            Q1_hat.col(idx) = (Q1_tilde * esW.eigenvectors().col(k)) / std::sqrt(ew(k));
        }

        // Simultaneously diagonalize M_{1,1} and c_d B_{1,1} inside the W_0-orthonormal basis
        const gsMatrix<T> M11_hat = Q1_hat.transpose() * (M0 * Q1_hat);
        gsEigen::SelfAdjointEigenSolver<typename gsMatrix<T>::Base> esM(M11_hat, gsEigen::ComputeEigenvectors);
        out.Q1 = Q1_hat * esM.eigenvectors();

        out.b11.resize(out.nD1);
        out.mu11.resize(out.nD1);
        for (index_t a = 0; a < out.nD1; ++a)
        {
            out.b11(a)  = std::max(T(0), c_d * out.Q1.col(a).dot(B0 * out.Q1.col(a)));
            out.mu11(a) = std::max(T(1e-30), out.Q1.col(a).dot(M0 * out.Q1.col(a)));
        }
    }

    static void applyUnscaledPatchSCMS(const PatchSCMS &scms,
                                       const gsMatrix<T> &R_mat,
                                       gsMatrix<T> &E_mat)
    {
        E_mat.setZero(scms.m1, scms.m2);
        const UnivariateSCMS &u1 = scms.u1;
        const UnivariateSCMS &u2 = scms.u2;

        // Block (0, 0): ((sigma0_1 + sigma0_2 + c0) M_{0,0}^{(1)} \otimes M_{0,0}^{(2)})^{-1}
        if (u1.nD0 > 0 && u2.nD0 > 0)
        {
            const gsMatrix<T> tmp1 = u1.Q0.transpose() * R_mat;
            const gsMatrix<T> Y00  = tmp1 * u2.Q0;
            const gsMatrix<T> Z00  = u1.solverM00->solve(Y00);
            const gsMatrix<T> Z00T = Z00.transpose();
            const gsMatrix<T> W00  = u2.solverM00->solve(Z00T).transpose();
            const gsMatrix<T> tmp2 = u1.Q0 * W00;
            E_mat.noalias() += scms.invSigma00 * (tmp2 * u2.Q0.transpose());
        }

        // Block (1, 0): (diag(b11 + (sigma0_2 + c0)*mu11) \otimes M_{0,0}^{(2)})^{-1}
        if (u1.nD1 > 0 && u2.nD0 > 0)
        {
            const gsMatrix<T> tmp1 = u1.Q1.transpose() * R_mat;
            gsMatrix<T> Y10        = tmp1 * u2.Q0;
            Y10 = scms.invDiag10.asDiagonal() * Y10;
            const gsMatrix<T> Y10T = Y10.transpose();
            const gsMatrix<T> W10  = u2.solverM00->solve(Y10T).transpose();
            const gsMatrix<T> tmp2 = u1.Q1 * W10;
            E_mat.noalias() += tmp2 * u2.Q0.transpose();
        }

        // Block (0, 1): (M_{0,0}^{(1)} \otimes diag(b22 + (sigma0_1 + c0)*mu22))^{-1}
        if (u1.nD0 > 0 && u2.nD1 > 0)
        {
            const gsMatrix<T> tmp1 = u1.Q0.transpose() * R_mat;
            const gsMatrix<T> Y01  = tmp1 * u2.Q1;
            gsMatrix<T> W01        = u1.solverM00->solve(Y01);
            W01 = W01 * scms.invDiag01.asDiagonal();
            const gsMatrix<T> tmp2 = u1.Q0 * W01;
            E_mat.noalias() += tmp2 * u2.Q1.transpose();
        }

        // Block (1, 1): diagonalized (c1 B_{1,1}^{(1)} \otimes M_{1,1}^{(2)} + M_{1,1}^{(1)} \otimes c2 B_{1,1}^{(2)})^{-1}
        if (u1.nD1 > 0 && u2.nD1 > 0)
        {
            gsMatrix<T> W11 = (u1.Q1.transpose() * R_mat) * u2.Q1;
            W11.array() *= scms.invDiag11.array();
            E_mat.noalias() += (u1.Q1 * W11) * u2.Q1.transpose();
        }
    }

    void relaxPatch(size_t pIdx, const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        const PatchSCMS &scms = m_patches[pIdx];
        const index_t m1 = scms.m1;
        const index_t m2 = scms.m2;
        const index_t nPatch = m1 * m2;
        const index_t ncols = rhs.cols();

        gsMatrix<T> R_mat(m1, m2), E_mat(m1, m2);
        for (index_t col = 0; col < ncols; ++col)
        {
            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = scms.dofs[i];
                if (row < 0)
                {
                    R_mat.data()[i] = T(0);
                    continue;
                }
                T acc = rhs(row, col);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    acc -= it.value() * x(it.col(), col);
                R_mat.data()[i] = acc;
            }

            applyUnscaledPatchSCMS(scms, R_mat, E_mat);
            for (index_t i = 0; i < nPatch; ++i)
            {
                const index_t row = scms.dofs[i];
                if (row >= 0)
                    x(row, col) += scms.scale * E_mat.data()[i];
            }
        }
    }

    void relaxSkeleton(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        if (!m_skelSolver)
            return;

        const index_t nSub  = static_cast<index_t>(m_skeletonDofs.size());
        const index_t ncols = rhs.cols();
        gsMatrix<T> r_sub(nSub, ncols);
        if (ncols == 1)
        {
            for (index_t i = 0; i < nSub; ++i)
            {
                const index_t row = m_skeletonDofs[i];
                T acc = rhs(row, 0);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    acc -= it.value() * x(it.col(), 0);
                r_sub(i, 0) = acc;
            }
        }
        else
        {
            for (index_t i = 0; i < nSub; ++i)
            {
                const index_t row = m_skeletonDofs[i];
                r_sub.row(i) = rhs.row(row);
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_matRow, row); it; ++it)
                    r_sub.row(i) -= it.value() * x.row(it.col());
            }
        }

        const gsMatrix<T> e_sub = m_skelSolver->solve(r_sub);
        for (index_t i = 0; i < nSub; ++i)
            x.row(m_skeletonDofs[i]) += e_sub.row(i);
    }

    void stepAdditive(const gsMatrix<T> &rhs, gsMatrix<T> &x) const
    {
        const gsMatrix<T> r = rhs - m_mat * x;
        const index_t ncols = rhs.cols();

        for (size_t pIdx = 0; pIdx < m_patches.size(); ++pIdx)
        {
            const PatchSCMS &scms = m_patches[pIdx];
            const index_t nPatch = scms.m1 * scms.m2;
            gsMatrix<T> R_mat(scms.m1, scms.m2), E_mat(scms.m1, scms.m2);
            for (index_t col = 0; col < ncols; ++col)
            {
                for (index_t i = 0; i < nPatch; ++i)
                {
                    const index_t row = scms.dofs[i];
                    R_mat.data()[i] = (row >= 0) ? r(row, col) : T(0);
                }
                applyUnscaledPatchSCMS(scms, R_mat, E_mat);
                const T w = T(0.5) * scms.scale;
                for (index_t i = 0; i < nPatch; ++i)
                {
                    const index_t row = scms.dofs[i];
                    if (row >= 0)
                        x(row, col) += w * E_mat.data()[i];
                }
            }
        }

        if (m_skelSolver)
        {
            const index_t nSub = static_cast<index_t>(m_skeletonDofs.size());
            gsMatrix<T> r_sub(nSub, ncols);
            for (index_t i = 0; i < nSub; ++i)
                r_sub.row(i) = r.row(m_skeletonDofs[i]);

            const gsMatrix<T> e_sub = m_skelSolver->solve(r_sub);
            for (index_t i = 0; i < nSub; ++i)
                x.row(m_skeletonDofs[i]) += T(0.5) * e_sub.row(i);
        }
    }

    gsSparseMatrix<T> m_mat;
    gsSparseMatrix<T, RowMajor> m_matRow;
    bool m_additive;
    std::vector<index_t> m_skeletonDofs;
    std::unique_ptr<SubSolver> m_skelSolver;
    std::vector<PatchSCMS> m_patches;
};

} // namespace gismo
