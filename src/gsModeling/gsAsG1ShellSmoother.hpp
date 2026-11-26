/** @file gsAsG1ShellSmoother.hpp

    @brief 3-component multi-patch and geometric multigrid smoothers & domain
           decomposition preconditioners for Kirchhoff-Love shells on AS-G1
           surfaces (Hasanova, S. Takacs, T. Takacs 2026).

    Implements:
      - gsNodalBlockGaussSeidelOp      ("nodal_gs")   : 3x3 nodal block GS removing A_3 tilt locking
      - gsAsG1ShellSkeletonSmootherOp  ("ring1")      : Nodal GS + 3-component Ring-1 skeleton Cholesky
      - gsAsG1ShellBlockSmootherOp     ("line_ring1") : 3-component u-strip / v-strip / Ring-1 LDLT
      - gsAsG1ShellScmsSmootherOp      ("scms" / "scms_hybrid"):
            O(p N_h) h- and p-robust tangent-normal frame SCMS (2nd-order membrane + 4th-order bending)
            + r_p = p + 1 overlapping 3-component skeleton Schwarz (matching Lemma 10 O(ph) support)
      - gsAsG1ShellPatchSchwarzPrecOp  ("schwarz")    : Single-level parallel Patch-Interior + Skeleton PCG

    This file is part of the G+Smoother / G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsModeling/gsAsG1Smoother.hpp>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiPatch.h>
#include <gsSolver/gsPreconditioner.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace gismo {

namespace shell_internal {

/// Group free global DOFs across components c in {0, 1, 2} that share the same
/// scalar AS-G1 basis node (patch, localArgyrisCol).
template <class T>
std::vector<std::vector<index_t>> collectNodalDofGroups(
    const std::vector<gsArgyrisEmbedding<T>>& argBasis,
    const gsDofMapper& mapper,
    index_t nComp = 3)
{
    const index_t nFree = mapper.freeSize();
    std::vector<bool> visited(nFree, false);
    std::vector<std::vector<index_t>> groups;
    groups.reserve(static_cast<size_t>(nFree / std::max<index_t>(nComp, 1) + 1));

    for (size_t p = 0; p < argBasis.size(); ++p)
    {
        const index_t nLoc = argBasis[p].matrix.cols();
        for (index_t loc = 0; loc < nLoc; ++loc)
        {
            // Check if any component at this scalar node is free and unvisited
            bool anyUnvisited = false;
            std::vector<index_t> grp;
            grp.reserve(nComp);
            for (index_t c = 0; c < nComp; ++c)
            {
                if (mapper.is_free(loc, static_cast<index_t>(p), c))
                {
                    const index_t g = mapper.index(loc, static_cast<index_t>(p), c);
                    if (g >= 0 && g < nFree)
                    {
                        if (!visited[g])
                            anyUnvisited = true;
                        grp.push_back(g);
                    }
                }
            }
            if (anyUnvisited && !grp.empty())
            {
                std::vector<index_t> fresh;
                fresh.reserve(grp.size());
                for (index_t g : grp)
                {
                    if (!visited[g])
                    {
                        visited[g] = true;
                        fresh.push_back(g);
                    }
                }
                if (!fresh.empty())
                    groups.push_back(std::move(fresh));
            }
        }
    }
    return groups;
}

/// Collect 3-component free DOFs on the AS-G1 skeleton + r_p interior boundary rings,
/// ordered node-by-node across components (u_x, u_y, u_z) for minimal Cholesky fill-in.
template <class T>
std::vector<index_t> collectShellCoupledStripDofs(
    const gsMultiPatch<T>& mp,
    const gsMultiBasis<T>& dbasis,
    const std::vector<gsArgyrisEmbedding<T>>& argBasis,
    const gsDofMapper& mapper,
    index_t r_p = 1,
    bool includeOuterBoundaryRing = true,
    index_t nComp = 3)
{
    const index_t nFree = mapper.freeSize();
    std::vector<bool> inSet(nFree, false);
    std::vector<index_t> dofs;

    auto addNode = [&](index_t patch, index_t localCol) {
        for (index_t c = 0; c < nComp; ++c)
        {
            if (mapper.is_free(localCol, patch, c))
            {
                const index_t g = mapper.index(localCol, patch, c);
                if (g >= 0 && g < nFree && !inSet[g])
                {
                    inSet[g] = true;
                    dofs.push_back(g);
                }
            }
        }
    };

    // Determine which sides of each patch are interfaces
    const size_t nP = mp.nPatches();
    std::vector<std::array<bool, 4>> isInterfaceSide(nP, {false, false, false, false});
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        isInterfaceSide[it->first().patch][it->first().m_index - 1] = true;
        isInterfaceSide[it->second().patch][it->second().m_index - 1] = true;
    }

    for (size_t p = 0; p < nP; ++p)
    {
        const gsVector<index_t, 13>& sz = argBasis[p].sizes;
        const index_t nInt = sz[0];
        const index_t nCols = argBasis[p].matrix.cols();

        // 1. All non-interior Argyris DOFs (edges level-0, level-1, and corners)
        for (index_t col = nInt; col < nCols; ++col)
            addNode(static_cast<index_t>(p), col);

        // 2. Interior boundary rings within distance r_p of interfaces (or all sides)
        //    (r_p == 0 disables rings; r_p < 0 selects the default p + 1 overlapping ring width)
        if (r_p == 0 || nInt == 0)
            continue;

        const gsTensorBSplineBasis<2, T>& tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T>&>(dbasis[p]);
        const index_t nu = tb.component(0).size();
        const index_t nv = tb.component(1).size();
        const index_t mu = nu - 4;
        const index_t mv = nv - 4;
        if (mu <= 0 || mv <= 0 || mu * mv != nInt)
            continue;

        const index_t effRp = (r_p < 0) ? (std::max(tb.degree(0), tb.degree(1)) + 1) : r_p;
        const bool west  = includeOuterBoundaryRing || isInterfaceSide[p][0];
        const bool east  = includeOuterBoundaryRing || isInterfaceSide[p][1];
        const bool south = includeOuterBoundaryRing || isInterfaceSide[p][2];
        const bool north = includeOuterBoundaryRing || isInterfaceSide[p][3];

        for (index_t j = 0; j < mv; ++j)
        {
            for (index_t i = 0; i < mu; ++i)
            {
                const bool nearSide =
                    (west  && i < effRp) ||
                    (east  && i >= mu - effRp) ||
                    (south && j < effRp) ||
                    (north && j >= mv - effRp);
                if (nearSide)
                    addNode(static_cast<index_t>(p), i + j * mu);
            }
        }
    }
    return dofs;
}

/// Sample mean tangent-normal orthonormal frame R_p = [e1, e2, n_mean] in SO(3)
/// and mean metric factors (detMean, invA11Mean, invA22Mean) on a 3D surface patch.
template <class T>
struct ShellPatchGeometryInfo
{
    gsMatrix<T, 3, 3> R; // columns: e1, e2, nMean (R^T R = I)
    T detMean = T(1);
    T invA11Mean = T(1);
    T invA22Mean = T(1);
    T invA11SqMean = T(1);
    T invA22SqMean = T(1);
    T curvFrobSqMean = T(0); // mean ||A^{-1} B||_F^2 (membrane-curvature coupling on normal mode)
};

template <class T>
ShellPatchGeometryInfo<T> sampleShellPatchGeometry(const gsGeometry<T>& geo, index_t nSamples = 4)
{
    ShellPatchGeometryInfo<T> info;
    const gsMatrix<T> supp = geo.support();
    gsMatrix<T> pts(2, nSamples * nSamples);
    for (index_t j = 0; j < nSamples; ++j)
        for (index_t i = 0; i < nSamples; ++i)
        {
            const T tu = (T(i) + T(0.5)) / T(nSamples);
            const T tv = (T(j) + T(0.5)) / T(nSamples);
            pts(0, i + j * nSamples) = supp(0, 0) + tu * (supp(0, 1) - supp(0, 0));
            pts(1, i + j * nSamples) = supp(1, 0) + tv * (supp(1, 1) - supp(1, 0));
        }

    gsMatrix<T> d1 = geo.deriv(pts);
    gsMatrix<T> d2 = geo.deriv2(pts);
    gsVector<T, 3> avgN = gsVector<T, 3>::Zero();
    gsVector<T, 3> avgA1 = gsVector<T, 3>::Zero();
    T sumDet = 0, sumA11 = 0, sumA22 = 0, sumA11Sq = 0, sumA22Sq = 0, sumCurv2 = 0;
    index_t valid = 0;

    for (index_t q = 0; q < pts.cols(); ++q)
    {
        gsVector<T, 3> A1, A2, A11, A22, A12;
        for (int c = 0; c < 3; ++c)
        {
            A1(c)  = d1(2 * c + 0, q);
            A2(c)  = d1(2 * c + 1, q);
            A11(c) = d2(3 * c + 0, q);
            A22(c) = d2(3 * c + 1, q);
            A12(c) = d2(3 * c + 2, q);
        }
        const T a11 = A1.squaredNorm();
        const T a22 = A2.squaredNorm();
        const T a12 = A1.dot(A2);
        const T detA = a11 * a22 - a12 * a12;
        if (detA <= T(1e-24)) continue;

        const T dS = std::sqrt(detA);
        const gsVector<T, 3> n = A1.cross(A2) / dS;
        avgN += n;
        avgA1 += A1 / std::max(A1.norm(), T(1e-15));
        const T inv11 = a22 / detA;
        const T inv22 = a11 / detA;
        const T inv12 = -a12 / detA;

        // Second fundamental form B_{alpha beta} = A_{alpha,beta} . n and ||A^{-1} B||_F^2
        const T b11 = A11.dot(n), b22 = A22.dot(n), b12 = A12.dot(n);
        const T m00 = inv11 * b11 + inv12 * b12;
        const T m01 = inv11 * b12 + inv12 * b22;
        const T m10 = inv12 * b11 + inv22 * b12;
        const T m11 = inv12 * b12 + inv22 * b22;
        const T curvFrob2 = m00 * m00 + T(2) * m01 * m10 + m11 * m11;

        sumDet   += dS;
        sumA11   += inv11;
        sumA22   += inv22;
        sumA11Sq += inv11 * inv11;
        sumA22Sq += inv22 * inv22;
        sumCurv2 += curvFrob2;
        ++valid;
    }

    if (valid > 0)
    {
        info.detMean        = sumDet / T(valid);
        info.invA11Mean     = sumA11 / T(valid);
        info.invA22Mean     = sumA22 / T(valid);
        info.invA11SqMean   = sumA11Sq / T(valid);
        info.invA22SqMean   = sumA22Sq / T(valid);
        info.curvFrobSqMean = sumCurv2 / T(valid);
    }

    gsVector<T, 3> nMean;
    if (avgN.norm() > T(1e-12))
        nMean = avgN.normalized();
    else
        nMean << T(0), T(0), T(1);
    gsVector<T, 3> e1 = avgA1 - avgA1.dot(nMean) * nMean;
    if (e1.norm() < T(1e-8))
    {
        gsVector<T, 3> trial;
        if (std::abs(nMean(0)) < T(0.9))
            trial << T(1), T(0), T(0);
        else
            trial << T(0), T(1), T(0);
        e1 = trial - trial.dot(nMean) * nMean;
    }
    e1.normalize();
    gsVector<T, 3> e2 = nMean.cross(e1).normalized();
    info.R.col(0) = e1;
    info.R.col(1) = e2;
    info.R.col(2) = nMean;
    return info;
}

/// Assemble 1D interior mass M0 and 2nd-order stiffness K0 on S^0_{p,k,Z} (indices 2..n-3).
template <class T>
void assembleUnivariateInteriorSecondOrder(const gsBSplineBasis<T>& b1d,
                                           gsSparseMatrix<T>& M0,
                                           gsSparseMatrix<T>& K0)
{
    const index_t n = b1d.size();
    const index_t m = n - 4;
    M0.resize(m, m);
    K0.resize(m, m);
    if (m <= 0) return;

    gsGaussRule<T> rule(b1d.degree() + 2);
    gsMatrix<T> qNodes;
    gsVector<T> qWeights;
    rule.mapToAll(b1d.knots().breaks(), qNodes, qWeights);

    std::vector<gsMatrix<T>> ders;
    b1d.evalAllDers_into(qNodes, 1, ders);
    gsMatrix<index_t> act = b1d.active(qNodes);

    gsSparseEntries<T> eM, eK;
    const index_t nQ = qNodes.cols();
    eM.reserve(nQ * (b1d.degree() + 1) * (b1d.degree() + 1));
    eK.reserve(nQ * (b1d.degree() + 1) * (b1d.degree() + 1));

    for (index_t q = 0; q < nQ; ++q)
    {
        const T w = qWeights(q);
        for (index_t a = 0; a < act.rows(); ++a)
        {
            const index_t ia = act(a, q) - 2;
            if (ia < 0 || ia >= m) continue;
            const T va = ders[0](a, q), da = ders[1](a, q);
            for (index_t b = 0; b < act.rows(); ++b)
            {
                const index_t ib = act(b, q) - 2;
                if (ib < 0 || ib >= m) continue;
                eM.add(ia, ib, w * va * ders[0](b, q));
                eK.add(ia, ib, w * da * ders[1](b, q));
            }
        }
    }
    M0.setFrom(eM);
    K0.setFrom(eK);
}

/// Assemble 1D interior mass M0 and 4th-order biharmonic stiffness B0 on S^0_{p,k,Z} (indices 2..n-3).
template <class T>
void assembleUnivariateInteriorFourthOrder(const gsBSplineBasis<T>& b1d,
                                           gsSparseMatrix<T>& M0,
                                           gsSparseMatrix<T>& B0)
{
    const index_t n = b1d.size();
    const index_t m = n - 4;
    M0.resize(m, m);
    B0.resize(m, m);
    if (m <= 0) return;

    gsGaussRule<T> rule(b1d.degree() + 2);
    gsMatrix<T> qNodes;
    gsVector<T> qWeights;
    rule.mapToAll(b1d.knots().breaks(), qNodes, qWeights);

    std::vector<gsMatrix<T>> ders;
    b1d.evalAllDers_into(qNodes, 2, ders);
    gsMatrix<index_t> act = b1d.active(qNodes);

    gsSparseEntries<T> eM, eB;
    const index_t nQ = qNodes.cols();
    eM.reserve(nQ * (b1d.degree() + 1) * (b1d.degree() + 1));
    eB.reserve(nQ * (b1d.degree() + 1) * (b1d.degree() + 1));

    for (index_t q = 0; q < nQ; ++q)
    {
        const T w = qWeights(q);
        for (index_t a = 0; a < act.rows(); ++a)
        {
            const index_t ia = act(a, q) - 2;
            if (ia < 0 || ia >= m) continue;
            const T va = ders[0](a, q), dda = ders[2](a, q);
            for (index_t b = 0; b < act.rows(); ++b)
            {
                const index_t ib = act(b, q) - 2;
                if (ib < 0 || ib >= m) continue;
                eM.add(ia, ib, w * va * ders[0](b, q));
                eB.add(ia, ib, w * dda * ders[2](b, q));
            }
        }
    }
    M0.setFrom(eM);
    B0.setFrom(eB);
}

/// Estimate lambda_max(M0^{-1} K0) via banded Cholesky power iteration.
template <class T>
T estimate1DMaxEig(const internal::BandedSPDCholesky<T>& cholM,
                   const gsSparseMatrix<T>& K)
{
    const index_t n = cholM.size();
    if (n <= 0) return T(1);
    gsVector<T> v(n);
    for (index_t i = 0; i < n; ++i)
        v(i) = std::sin(T(i + 1) * T(1.37)) + T(0.5) * std::cos(T(i + 1) * T(0.73));
    v.normalize();
    T lamMax = T(1), prev = T(0);
    gsVector<T> Kv(n), w(n);
    for (int it = 0; it < 20; ++it)
    {
        Kv.noalias() = K * v;
        const T vKv = v.dot(Kv);
        if (vKv <= T(1e-30)) break;
        w = Kv;
        cholM.solveLeftInPlace(w);
        const T wKv = w.dot(Kv);
        const T lam = (wKv > T(0)) ? (wKv / vKv) : T(0);
        if (lam > T(0)) lamMax = lam;
        const T wn = w.norm();
        if (wn <= T(1e-30)) break;
        v = w / wn;
        if (it >= 5 && prev > T(0) && std::abs(lamMax - prev) <= T(1e-4) * lamMax)
            break;
        prev = lamMax;
    }
    return T(1.02) * lamMax;
}

} // namespace shell_internal

// ====================================================================
// 1. 3x3 Nodal Block Gauss-Seidel Smoother (A-symmetric step / stepT)
// ====================================================================
template <class T>
class gsNodalBlockGaussSeidelOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::unique_ptr<gsNodalBlockGaussSeidelOp> uPtr;

    gsNodalBlockGaussSeidelOp(const gsSparseMatrix<T>& mat,
                              const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                              const gsDofMapper& mapper,
                              index_t nComp = 3)
        : m_mat(mat)
    {
        const auto groups = shell_internal::collectNodalDofGroups(argBasis, mapper, nComp);
        m_blocks.reserve(groups.size());

        for (const auto& grp : groups)
        {
            NodalBlock blk;
            blk.sz = static_cast<int>(grp.size());
            for (int a = 0; a < blk.sz; ++a)
                blk.dofs[a] = grp[a];

            // Extract small dense symmetric diagonal block (1x1, 2x2, or 3x3) and invert
            gsMatrix<T> D(blk.sz, blk.sz);
            D.setZero();
            for (int a = 0; a < blk.sz; ++a)
            {
                const index_t row = blk.dofs[a];
                for (typename gsSparseMatrix<T>::InnerIterator it(m_mat, row); it; ++it)
                {
                    for (int b = 0; b < blk.sz; ++b)
                        if (it.row() == blk.dofs[b])
                            D(a, b) = it.value();
                }
            }
            // Regularize tiny/zero diagonal safely
            for (int a = 0; a < blk.sz; ++a)
                if (std::abs(D(a, a)) < T(1e-30))
                    D(a, a) = T(1);

            gsMatrix<T> Dinv = D.ldlt().solve(gsMatrix<T>::Identity(blk.sz, blk.sz));
            for (int a = 0; a < blk.sz; ++a)
                for (int b = 0; b < blk.sz; ++b)
                    blk.invD[a * 3 + b] = Dinv(a, b);
            m_blocks.push_back(blk);
        }
    }

    static uPtr make(const gsSparseMatrix<T>& mat,
                     const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                     const gsDofMapper& mapper,
                     index_t nComp = 3)
    {
        return uPtr(new gsNodalBlockGaussSeidelOp(mat, argBasis, mapper, nComp));
    }

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        const index_t ncols = rhs.cols();
        for (size_t b = 0; b < m_blocks.size(); ++b)
            relaxBlock(m_blocks[b], rhs, x, ncols);
    }

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        const index_t ncols = rhs.cols();
        for (size_t b = m_blocks.size(); b-- > 0;)
            relaxBlock(m_blocks[b], rhs, x, ncols);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    struct NodalBlock
    {
        int sz = 0;
        index_t dofs[3] = {-1, -1, -1};
        T invD[9] = {0};
    };

    void relaxBlock(const NodalBlock& blk, const gsMatrix<T>& rhs, gsMatrix<T>& x, index_t ncols) const
    {
        for (index_t col = 0; col < ncols; ++col)
        {
            const T* rPtr = rhs.data() + col * rhs.rows();
            T* xPtr = x.data() + col * x.rows();
            T res[3] = {0, 0, 0};
            for (int a = 0; a < blk.sz; ++a)
            {
                const index_t row = blk.dofs[a];
                T acc = rPtr[row];
                for (typename gsSparseMatrix<T>::InnerIterator it(m_mat, row); it; ++it)
                    acc -= it.value() * xPtr[it.row()];
                res[a] = acc;
            }
            for (int a = 0; a < blk.sz; ++a)
            {
                T upd = T(0);
                for (int c = 0; c < blk.sz; ++c)
                    upd += blk.invD[a * 3 + c] * res[c];
                xPtr[blk.dofs[a]] += upd;
            }
        }
    }

    gsSparseMatrix<T> m_mat;
    std::vector<NodalBlock> m_blocks;
};

// ====================================================================
// 2. Ring-1 Skeleton + 3x3 Nodal GS Smoother ("ring1")
// ====================================================================
template <class T>
class gsAsG1ShellSkeletonSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::unique_ptr<gsAsG1ShellSkeletonSmootherOp> uPtr;

    gsAsG1ShellSkeletonSmootherOp(const gsSparseMatrix<T>& mat,
                                  const gsMultiPatch<T>& mp,
                                  const gsMultiBasis<T>& dbasis,
                                  const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                                  const gsDofMapper& mapper,
                                  index_t r_p = 1)
        : m_mat(mat)
    {
        m_nodalGs = gsNodalBlockGaussSeidelOp<T>::make(m_mat, argBasis, mapper, 3);
        std::vector<index_t> globToLoc(m_mat.rows(), -1);
        const std::vector<index_t> skelDofs =
            shell_internal::collectShellCoupledStripDofs(mp, dbasis, argBasis, mapper, r_p, true, 3);
        m_skeleton.init(m_mat, skelDofs, globToLoc, 48);
    }

    static uPtr make(const gsSparseMatrix<T>& mat,
                     const gsMultiPatch<T>& mp,
                     const gsMultiBasis<T>& dbasis,
                     const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                     const gsDofMapper& mapper,
                     index_t r_p = 1)
    {
        return uPtr(new gsAsG1ShellSkeletonSmootherOp(mat, mp, dbasis, argBasis, mapper, r_p));
    }

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        m_nodalGs->step(rhs, x);
        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        m_skeleton.relax(m_mat, rhs, x, internal::isExactlyZero(x));
        m_nodalGs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    typename gsNodalBlockGaussSeidelOp<T>::uPtr m_nodalGs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
};

// ====================================================================
// 3. 3-Component Line + Ring-1 Multiplicative Block Smoother ("line_ring1")
// ====================================================================
template <class T>
class gsAsG1ShellBlockSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::unique_ptr<gsAsG1ShellBlockSmootherOp> uPtr;

    gsAsG1ShellBlockSmootherOp(const gsSparseMatrix<T>& mat,
                               const gsMultiPatch<T>& mp,
                               const gsMultiBasis<T>& dbasis,
                               const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                               const gsDofMapper& mapper,
                               index_t r_p = 1)
        : m_mat(mat)
    {
        std::vector<index_t> globToLoc(m_mat.rows(), -1);
        const index_t nFree = m_mat.rows();

        // 1. Interior u-lines and v-lines (node-interleaved across c = 0, 1, 2)
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(dbasis[p]);
            const index_t mu = tb.component(0).size() - 4;
            const index_t mv = tb.component(1).size() - 4;
            if (mu <= 0 || mv <= 0 || mu * mv != argBasis[p].sizes[0])
                continue;

            // u-lines (fixed j)
            for (index_t j = 0; j < mv; ++j)
            {
                std::vector<index_t> line;
                line.reserve(3 * mu);
                for (index_t i = 0; i < mu; ++i)
                    for (index_t c = 0; c < 3; ++c)
                        if (mapper.is_free(i + j * mu, static_cast<index_t>(p), c))
                        {
                            const index_t g = mapper.index(i + j * mu, static_cast<index_t>(p), c);
                            if (g >= 0 && g < nFree) line.push_back(g);
                        }
                internal::SubspaceCholeskyBlock<T> blk;
                if (blk.init(m_mat, line, globToLoc, 36))
                    m_blocks.push_back(std::move(blk));
            }

            // v-lines (fixed i)
            for (index_t i = 0; i < mu; ++i)
            {
                std::vector<index_t> line;
                line.reserve(3 * mv);
                for (index_t j = 0; j < mv; ++j)
                    for (index_t c = 0; c < 3; ++c)
                        if (mapper.is_free(i + j * mu, static_cast<index_t>(p), c))
                        {
                            const index_t g = mapper.index(i + j * mu, static_cast<index_t>(p), c);
                            if (g >= 0 && g < nFree) line.push_back(g);
                        }
                internal::SubspaceCholeskyBlock<T> blk;
                if (blk.init(m_mat, line, globToLoc, 36))
                    m_blocks.push_back(std::move(blk));
            }
        }

        // 2. Coupled 3-component Ring-1 / Strip skeleton block
        const std::vector<index_t> skelDofs =
            shell_internal::collectShellCoupledStripDofs(mp, dbasis, argBasis, mapper, r_p, true, 3);
        internal::SubspaceCholeskyBlock<T> skelBlk;
        if (skelBlk.init(m_mat, skelDofs, globToLoc, 64))
            m_blocks.push_back(std::move(skelBlk));
    }

    static uPtr make(const gsSparseMatrix<T>& mat,
                     const gsMultiPatch<T>& mp,
                     const gsMultiBasis<T>& dbasis,
                     const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                     const gsDofMapper& mapper,
                     index_t r_p = 1)
    {
        return uPtr(new gsAsG1ShellBlockSmootherOp(mat, mp, dbasis, argBasis, mapper, r_p));
    }

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        for (size_t i = 0; i < m_blocks.size(); ++i)
            m_blocks[i].relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        for (size_t i = m_blocks.size(); i-- > 0;)
            m_blocks[i].relax(m_mat, rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    std::vector<internal::SubspaceCholeskyBlock<T>> m_blocks;
};

// ====================================================================
// 4. O(p N_h) h- and p-Robust Tangent-Normal SCMS + r_p = p+1 Skeleton
//    Smoother ("scms" / "scms_hybrid", Hasanova-Takacs-Takacs 2026)
// ====================================================================
template <class T>
class gsAsG1ShellScmsSmootherOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::unique_ptr<gsAsG1ShellScmsSmootherOp> uPtr;

    gsAsG1ShellScmsSmootherOp(const gsSparseMatrix<T>& mat,
                              const gsMultiPatch<T>& mp,
                              const gsMultiBasis<T>& dbasis,
                              const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                              const gsDofMapper& mapper,
                              T E, T nu, T thickness,
                              index_t r_p = -1,
                              bool useNodalGS = true)
        : m_mat(mat), m_useNodalGS(useNodalGS)
    {
        if (m_useNodalGS)
            m_nodalGs = gsNodalBlockGaussSeidelOp<T>::make(m_mat, argBasis, mapper, 3);

        // Default overlapping strip width r_p = p + 1 matching Lemma 10 O(ph) Hermite support
        index_t maxDeg = 3;
        for (size_t p = 0; p < mp.nPatches(); ++p)
            maxDeg = std::max<index_t>(maxDeg, std::max(dbasis[p].degree(0), dbasis[p].degree(1)));
        const index_t effRp = (r_p < 0) ? (maxDeg + 1) : r_p;

        std::vector<index_t> globToLoc(m_mat.rows(), -1);
        const std::vector<index_t> skelDofs =
            shell_internal::collectShellCoupledStripDofs(mp, dbasis, argBasis, mapper, effRp, true, 3);
        m_skeleton.init(m_mat, skelDofs, globToLoc, 64);

        const T memFac  = E * thickness / (T(1) - nu * nu);
        const T bendFac = E * thickness * thickness * thickness / (T(12) * (T(1) - nu * nu));

        // Build per-patch O(p m_1 m_2) tangent-normal rotated SCMS blocks
        m_patches.reserve(mp.nPatches());
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(dbasis[p]);
            const index_t mu = tb.component(0).size() - 4;
            const index_t mv = tb.component(1).size() - 4;
            if (mu <= 0 || mv <= 0 || mu * mv != argBasis[p].sizes[0])
                continue;

            PatchShellSCMS ps;
            ps.mu = mu;
            ps.mv = mv;
            ps.nInt = mu * mv;
            ps.dofs.resize(3 * ps.nInt, -1);
            bool anyFree = false;
            for (index_t k = 0; k < ps.nInt; ++k)
            {
                for (index_t c = 0; c < 3; ++c)
                {
                    if (mapper.is_free(k, static_cast<index_t>(p), c))
                    {
                        ps.dofs[c * ps.nInt + k] = mapper.index(k, static_cast<index_t>(p), c);
                        anyFree = true;
                    }
                }
            }
            if (!anyFree) continue;

            const auto gInfo = shell_internal::sampleShellPatchGeometry(mp.patch(p), 4);
            ps.R = gInfo.R;

            // 1D banded mass Cholesky + 2nd-order membrane & 4th-order bending SCMS symbols
            build1DMembraneBending(tb.component(0), ps.u1);
            build1DMembraneBending(tb.component(1), ps.u2);

            const T cMem1  = memFac  * gInfo.detMean * gInfo.invA11Mean;
            const T cMem2  = memFac  * gInfo.detMean * gInfo.invA22Mean;
            const T cBend1 = bendFac * gInfo.detMean * gInfo.invA11SqMean;
            const T cBend2 = bendFac * gInfo.detMean * gInfo.invA22SqMean;
            // Curvature membrane stiffness c_curv * M_u \otimes M_v acting on the normal mode u_n
            const T cCurv  = memFac  * gInfo.detMean * gInfo.curvFrobSqMean;

            ps.invMemSig  = T(1) / std::max(cMem1 * ps.u1.sigMem + cMem2 * ps.u2.sigMem, T(1e-30));
            ps.invBendSig = T(1) / std::max(cBend1 * ps.u1.sigBend + cBend2 * ps.u2.sigBend + cCurv, T(1e-30));

            // Calibrate global patch damping factor scale so lambda_max(M_p^{-1} A_pp) <= 0.95
            gsSparseMatrix<T> A_pp = internal::extractSymmetricSubmatrix(m_mat, ps.dofs, globToLoc, 48);
            ps.scale = calibratePatchScale(ps, A_pp);
            m_patches.push_back(std::move(ps));
        }
    }

    static uPtr make(const gsSparseMatrix<T>& mat,
                     const gsMultiPatch<T>& mp,
                     const gsMultiBasis<T>& dbasis,
                     const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                     const gsDofMapper& mapper,
                     T E, T nu, T thickness,
                     index_t r_p = -1,
                     bool useNodalGS = true)
    {
        return uPtr(new gsAsG1ShellScmsSmootherOp(
            mat, mp, dbasis, argBasis, mapper, E, nu, thickness, r_p, useNodalGS));
    }

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        if (m_useNodalGS)
            m_nodalGs->step(rhs, x);

        const bool xZero = internal::isExactlyZero(x);
        #pragma omp parallel for schedule(static) if(m_patches.size() > 1)
        for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
            relaxPatch(m_patches[p], rhs, x, xZero);

        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        m_skeleton.relax(m_mat, rhs, x, internal::isExactlyZero(x));

        const bool xZero = internal::isExactlyZero(x);
        #pragma omp parallel for schedule(static) if(m_patches.size() > 1)
        for (index_t p = 0; p < static_cast<index_t>(m_patches.size()); ++p)
            relaxPatch(m_patches[p], rhs, x, xZero);

        if (m_useNodalGS)
            m_nodalGs->stepT(rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    struct Univariate1D
    {
        internal::BandedSPDCholesky<T> cholM;
        T sigMem = T(1);
        T sigBend = T(1);
    };

    struct PatchShellSCMS
    {
        index_t mu = 0, mv = 0, nInt = 0;
        std::vector<index_t> dofs; // 3 * nInt (comp * nInt + k)
        gsMatrix<T, 3, 3> R;
        Univariate1D u1, u2;
        T invMemSig = T(1);
        T invBendSig = T(1);
        T scale = T(1);
    };

    static void build1DMembraneBending(const gsBSplineBasis<T>& b1d, Univariate1D& out)
    {
        gsSparseMatrix<T> M0, K0, M0b, B0;
        shell_internal::assembleUnivariateInteriorSecondOrder(b1d, M0, K0);
        shell_internal::assembleUnivariateInteriorFourthOrder(b1d, M0b, B0);
        out.cholM.compute(M0);
        out.sigMem  = shell_internal::estimate1DMaxEig(out.cholM, K0);
        out.sigBend = shell_internal::estimate1DMaxEig(out.cholM, B0);
    }

    /// Apply unscaled rotated patch-interior SCMS: y = M_p^{-1} r (O(p * mu * mv) cost).
    static void applyPatchPrecond(const PatchShellSCMS& ps, const gsVector<T>& rLoc, gsVector<T>& yLoc)
    {
        const index_t mu = ps.mu, mv = ps.mv, n = ps.nInt;
        yLoc.resize(3 * n);

        // Rotate residual into patch mean tangent-normal frame: r_rot = R^T * r
        gsMatrix<T> C0(mu, mv), C1(mu, mv), C2(mu, mv);
        for (index_t k = 0; k < n; ++k)
        {
            const index_t i = k % mu, j = k / mu;
            const T rx = rLoc(0 * n + k);
            const T ry = rLoc(1 * n + k);
            const T rz = rLoc(2 * n + k);
            C0(i, j) = ps.R(0, 0) * rx + ps.R(1, 0) * ry + ps.R(2, 0) * rz;
            C1(i, j) = ps.R(0, 1) * rx + ps.R(1, 1) * ry + ps.R(2, 1) * rz;
            C2(i, j) = ps.R(0, 2) * rx + ps.R(1, 2) * ry + ps.R(2, 2) * rz;
        }

        // Apply O(p * mu * mv) banded tensor-product mass inverse with membrane/bending scaling
        ps.u1.cholM.solveLeftInPlace(C0);  ps.u2.cholM.solveRightInPlace(C0); C0 *= ps.invMemSig;
        ps.u1.cholM.solveLeftInPlace(C1);  ps.u2.cholM.solveRightInPlace(C1); C1 *= ps.invMemSig;
        ps.u1.cholM.solveLeftInPlace(C2);  ps.u2.cholM.solveRightInPlace(C2); C2 *= ps.invBendSig;

        // Rotate back to Cartesian coordinates: y = R * y_rot
        for (index_t k = 0; k < n; ++k)
        {
            const index_t i = k % mu, j = k / mu;
            const T y0 = C0(i, j), y1 = C1(i, j), y2 = C2(i, j);
            yLoc(0 * n + k) = ps.R(0, 0) * y0 + ps.R(0, 1) * y1 + ps.R(0, 2) * y2;
            yLoc(1 * n + k) = ps.R(1, 0) * y0 + ps.R(1, 1) * y1 + ps.R(1, 2) * y2;
            yLoc(2 * n + k) = ps.R(2, 0) * y0 + ps.R(2, 1) * y1 + ps.R(2, 2) * y2;
        }
    }

    static T calibratePatchScale(const PatchShellSCMS& ps, const gsSparseMatrix<T>& A_pp)
    {
        const index_t m = A_pp.rows();
        if (m == 0) return T(1);
        gsVector<T> v(m), Av(m), w(m);
        for (index_t i = 0; i < m; ++i)
            v(i) = std::sin(T(i + 1) * T(1.13)) + T(0.3) * std::cos(T(i + 1) * T(0.47));
        v.normalize();
        T lamMax = T(1), prev = T(0);
        for (int it = 0; it < 20; ++it)
        {
            Av.noalias() = A_pp * v;
            const T vAv = v.dot(Av);
            if (vAv <= T(1e-30)) break;
            applyPatchPrecond(ps, Av, w);
            const T wAv = w.dot(Av);
            const T lam = (wAv > T(0)) ? (wAv / vAv) : T(0);
            if (lam > T(0)) lamMax = lam;
            const T wn = w.norm();
            if (wn <= T(1e-30)) break;
            v = w / wn;
            if (it >= 5 && prev > T(0) && std::abs(lamMax - prev) <= T(1e-3) * lamMax)
                break;
            prev = lamMax;
        }
        return (lamMax > T(1e-12)) ? (T(0.95) / lamMax) : T(1);
    }

    void relaxPatch(const PatchShellSCMS& ps, const gsMatrix<T>& rhs, gsMatrix<T>& x, bool xZero) const
    {
        const index_t m = static_cast<index_t>(ps.dofs.size());
        gsVector<T> rLoc(m), yLoc(m);
        for (index_t col = 0; col < rhs.cols(); ++col)
        {
            const T* rPtr = rhs.data() + col * rhs.rows();
            T* xPtr = x.data() + col * x.rows();
            for (index_t k = 0; k < m; ++k)
            {
                const index_t row = ps.dofs[k];
                if (row < 0) { rLoc(k) = T(0); continue; }
                T acc = rPtr[row];
                if (!xZero)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_mat, row); it; ++it)
                        acc -= it.value() * xPtr[it.row()];
                rLoc(k) = acc;
            }
            applyPatchPrecond(ps, rLoc, yLoc);
            for (index_t k = 0; k < m; ++k)
                if (ps.dofs[k] >= 0)
                    xPtr[ps.dofs[k]] += ps.scale * yLoc(k);
        }
    }

    gsSparseMatrix<T> m_mat;
    bool m_useNodalGS;
    typename gsNodalBlockGaussSeidelOp<T>::uPtr m_nodalGs;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
    std::vector<PatchShellSCMS> m_patches;
};

// ====================================================================
// 5. Single-Level Multi-Patch Overlapping Patch-Interior + Skeleton
//    Domain Decomposition Preconditioner ("schwarz", SPD for PCG)
// ====================================================================
template <class T>
class gsAsG1ShellPatchSchwarzPrecOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:
    typedef memory::unique_ptr<gsAsG1ShellPatchSchwarzPrecOp> uPtr;

    gsAsG1ShellPatchSchwarzPrecOp(const gsSparseMatrix<T>& mat,
                                  const gsMultiPatch<T>& mp,
                                  const gsMultiBasis<T>& dbasis,
                                  const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                                  const gsDofMapper& mapper,
                                  index_t r_p = -1)
        : m_mat(mat)
    {
        const index_t nFree = m_mat.rows();
        std::vector<index_t> globToLoc(nFree, -1);

        // 1. Decoupled patch-interior 3-component exact sparse Cholesky blocks (parallelizable!)
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const index_t nInt = argBasis[p].sizes[0];
            if (nInt <= 0) continue;
            std::vector<index_t> pDofs;
            pDofs.reserve(3 * nInt);
            for (index_t k = 0; k < nInt; ++k)
                for (index_t c = 0; c < 3; ++c)
                    if (mapper.is_free(k, static_cast<index_t>(p), c))
                    {
                        const index_t g = mapper.index(k, static_cast<index_t>(p), c);
                        if (g >= 0 && g < nFree) pDofs.push_back(g);
                    }
            internal::SubspaceCholeskyBlock<T> blk;
            if (blk.init(m_mat, pDofs, globToLoc, 48))
                m_patchInteriors.push_back(std::move(blk));
        }

        // 2. Overlapping 3-component skeleton block (r_p = p + 1 by default)
        index_t maxDeg = 3;
        for (size_t p = 0; p < mp.nPatches(); ++p)
            maxDeg = std::max<index_t>(maxDeg, std::max(dbasis[p].degree(0), dbasis[p].degree(1)));
        const index_t effRp = (r_p < 0) ? (maxDeg + 1) : r_p;

        const std::vector<index_t> skelDofs =
            shell_internal::collectShellCoupledStripDofs(mp, dbasis, argBasis, mapper, effRp, true, 3);
        m_skeleton.init(m_mat, skelDofs, globToLoc, 64);
    }

    static uPtr make(const gsSparseMatrix<T>& mat,
                     const gsMultiPatch<T>& mp,
                     const gsMultiBasis<T>& dbasis,
                     const std::vector<gsArgyrisEmbedding<T>>& argBasis,
                     const gsDofMapper& mapper,
                     index_t r_p = -1)
    {
        return uPtr(new gsAsG1ShellPatchSchwarzPrecOp(mat, mp, dbasis, argBasis, mapper, r_p));
    }

    // Symmetric multiplicative preconditioner: interiors -> skeleton -> interiors (exact SPD!)
    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const override
    {
        x.setZero(input.rows(), input.cols());
        step(input, x);
        stepT(input, x);
    }

    void step(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        const bool xZero = internal::isExactlyZero(x);
        #pragma omp parallel for schedule(static) if(m_patchInteriors.size() > 1)
        for (index_t p = 0; p < static_cast<index_t>(m_patchInteriors.size()); ++p)
            m_patchInteriors[p].relax(m_mat, rhs, x, xZero);
        m_skeleton.relax(m_mat, rhs, x);
    }

    void stepT(const gsMatrix<T>& rhs, gsMatrix<T>& x) const override
    {
        m_skeleton.relax(m_mat, rhs, x, internal::isExactlyZero(x));
        #pragma omp parallel for schedule(static) if(m_patchInteriors.size() > 1)
        for (index_t p = 0; p < static_cast<index_t>(m_patchInteriors.size()); ++p)
            m_patchInteriors[p].relax(m_mat, rhs, x);
    }

    typename gsLinearOperator<T>::Ptr underlyingOp() const override { return makeMatrixOp(m_mat); }
    index_t rows() const override { return m_mat.rows(); }
    index_t cols() const override { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    std::vector<internal::SubspaceCholeskyBlock<T>> m_patchInteriors;
    internal::SubspaceCholeskyBlock<T> m_skeleton;
};

} // namespace gismo
