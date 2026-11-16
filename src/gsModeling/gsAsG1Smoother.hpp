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
 *        augmented with the adjacent two-layer interior B-spline ring (Ring 1).
 *
 * Including the first two interior B-spline layers (i_n in {2,3} or {n-4,n-3})
 * alongside free edge-trace, edge-crossing, and corner-jet DOFs resolves the
 * O(h^{-4}) normal-curvature (d_nn u) bending locking across interfaces and
 * boundaries while keeping the subspace dimension O(sqrt(N)).
 */
template <class T>
std::vector<index_t> collectCoupledRing1Dofs(
    const gsMultiPatch<T> &mp,
    const gsMultiBasis<T> &dbasis,
    const std::vector<gsArgyrisEmbedding<T>> &argBasis,
    const gsDofMapper &mapper,
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
                const bool rW = (i1 <= 3), rE = (i1 >= n1 - 4);
                const bool rS = (i2 <= 3), rN = (i2 >= n2 - 4);
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

} // namespace gismo
