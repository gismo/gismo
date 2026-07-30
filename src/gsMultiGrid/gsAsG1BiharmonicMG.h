/** @file gsAsG1BiharmonicMG.h
 *
 *  @brief Geometric multigrid solver for the biharmonic equation on
 *         Analysis-Suitable G1 (AS-G1) multi-patch spaces.
 *
 *  In contrast to the earlier hand-rolled solver (gsAsG1Multigrid), this
 *  implementation
 *    - builds the AS-G1 degree-of-freedom coupling with the shared helper
 *      \ref makeMapperForArgyrisBasis (from gsModeling/gsAsG1Basis.hpp) and the
 *      gluing data from \ref computeGluingData (gsModeling/gsAsG1Domain.hpp),
 *      instead of re-implementing the interface matching inline, and
 *    - assembles the level hierarchy into the generic \ref gsMultiGridOp from
 *      gsMultiGrid/gsMultiGrid.h (with Galerkin coarsening, symmetric
 *      Gauss-Seidel smoothing and a configurable number of smoothing steps),
 *      used as a preconditioner for a conjugate-gradient iteration, instead of
 *      a bespoke V-cycle.
 *
 *  This file is part of the G+Smo library.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Author(s): F. Hasanova, S. Takacs
 */

#pragma once

#include <algorithm>
#include <memory>
#include <set>
#include <vector>

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsDofMapper.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsPde/gsBoundaryConditions.h>

#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

#include <gsMultiGrid/gsMultiGrid.h>
#include <gsSolver/gsSimplePreconditioners.h>
#include <gsSolver/gsConjugateGradient.h>

namespace gismo
{

/** @brief AS-G1 biharmonic geometric multigrid solver.
 *
 *  Builds a hierarchy of uniformly refined AS-G1 spaces (levels
 *  @a minRefinements .. @a maxRefinements) and wraps it into a
 *  \ref gsMultiGridOp preconditioner.  The finest-level stiffness matrix and
 *  the exact (space-nested) inter-level prolongations are provided explicitly;
 *  the coarse-grid operators are obtained by Galerkin projection inside
 *  \ref gsMultiGridOp.
 *
 *  @ingroup Solver
 */
template<typename T>
class gsAsG1BiharmonicMG
{
public:
    typedef gsSparseMatrix<T, RowMajor> SpTransfer;

    /// Constructor: sets up the whole multigrid hierarchy.
    gsAsG1BiharmonicMG(const gsMultiPatch<T>& geometry,
                       index_t degree         = 3,
                       index_t minRefinements = 2,
                       index_t maxRefinements = 4,
                       index_t numGaussPerSpan = 0)
    : m_degree(degree < 3 ? 3 : degree),
      m_minRef(minRefinements < 0 ? 0 : minRefinements),
      m_maxRef(maxRefinements < m_minRef ? m_minRef : maxRefinements),
      m_numGauss(numGaussPerSpan),
      m_numLevels(m_maxRef - m_minRef + 1)
    {
        setup(geometry);
    }

    /// Number of grid levels in the hierarchy.
    index_t numLevels() const { return m_numLevels; }

    /// Number of free (interior, coupled) DOFs on the finest level.
    index_t freeDofs() const { return m_Tfree.back().cols(); }

    /// Finest-level AS-G1 biharmonic stiffness matrix.
    const gsSparseMatrix<T>& stiffnessMatrix() const { return m_Kfine; }

    /// Finest-level transformation matrix (disjoint B-spline DOFs <- free DOFs).
    const gsSparseMatrix<T>& transformFine() const { return m_Tfree.back(); }

    /// Finest-level multi-patch (degree-elevated and refined).
    const gsMultiPatch<T>& fineMultiPatch() const { return m_mp.back(); }

    /// The underlying multigrid preconditioner.
    gsMultiGridOp<T>& multigrid() const { return *m_mg; }

    /** @brief Solve @a K x = @a b with CG preconditioned by the multigrid V-cycle.
     *
     *  @param[in]  b          Right-hand side in finest free-DOF space.
     *  @param[out] x          Solution vector.
     *  @param[in]  maxIter    Maximum number of CG iterations.
     *  @param[in]  tol        Relative residual tolerance.
     *  @param[out] finalError Achieved relative residual (optional).
     *  @return                Number of CG iterations performed.
     */
    index_t solve(const gsMatrix<T>& b, gsMatrix<T>& x,
                  index_t maxIter = 200, T tol = T(1e-8),
                  T* finalError = nullptr) const
    {
        gsConjugateGradient<T> cg(m_Kfine, m_mgPtr);
        cg.setMaxIterations(maxIter);
        cg.setTolerance(tol);
        x.setZero(m_Kfine.rows(), b.cols());
        cg.solve(b, x);
        if (finalError) *finalError = cg.error();
        return cg.iterations();
    }

    /// Number of eliminated (boundary) DOFs on the finest level.
    index_t boundaryDofs() const { return m_TbndFine.cols(); }

    /// Finest-level boundary transformation (disjoint B-spline DOFs <- boundary DOFs).
    const gsSparseMatrix<T>& transformBoundaryFine() const { return m_TbndFine; }

    /// Finest-level interior/boundary coupling block K_fb = T_free^T K T_bnd.
    const gsSparseMatrix<T>& couplingBlock() const { return m_Kfb; }

    /** @brief Boundary DOF values realising an inhomogeneous clamped condition.
     *
     *  Returns the boundary-DOF coefficients of the L2 projection of @a uex
     *  onto the finest AS-G1 space. These carry the trace and normal-derivative
     *  data of @a uex and can be prescribed as inhomogeneous boundary values.
     */
    gsMatrix<T> projectBoundaryValues(const gsFunction<T>& uex) const
    {
        using namespace gismo::expr;
        const gsMultiPatch<T>& mp = m_mp.back();
        gsPiecewiseFunction<T> uf;
        for (size_t i = 0; i < mp.nPatches(); ++i) uf.addPiece(uex);

        gsMultiBasis<T> dbasis(mp);
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(dbasis.domain());
        auto G = A.getMap(mp);
        auto u = A.getSpace(dbasis);
        auto f = A.getCoeff(uf, G);
        A.initSystem();
        A.assemble(u * u.tr() * meas(G), u * f * meas(G));

        const gsSparseMatrix<T> Tfull = hconcat(m_Tfree.back(), m_TbndFine);
        gsSparseMatrix<T> Mfull = Tfull.transpose() * A.matrix() * Tfull;
        gsMatrix<T> rhs = Tfull.transpose() * A.rhs();
        typename gsSparseSolver<T>::SimplicialLDLT solver;
        solver.compute(Mfull);
        gsMatrix<T> clift = solver.solve(rhs);
        return clift.bottomRows(m_TbndFine.cols());
    }

    /// Apply the boundary lifting to a free-DOF right-hand side:
    /// F_free  ->  F_free - K_fb * cBnd.
    gsMatrix<T> liftedRhs(const gsMatrix<T>& Ffree, const gsMatrix<T>& cBnd) const
    { return Ffree - m_Kfb * cBnd; }

    /// Reconstruct a multi-patch field from a finest-level free-DOF vector
    /// (boundary DOFs are homogeneous, i.e. zero).
    gsMultiPatch<T> reconstruct(const gsMatrix<T>& cFree) const
    { return reconstructFromDisjoint(m_Tfree.back() * cFree); }

    /// Reconstruct a multi-patch field from free and prescribed boundary DOFs.
    gsMultiPatch<T> reconstruct(const gsMatrix<T>& cFree, const gsMatrix<T>& cBnd) const
    { return reconstructFromDisjoint(m_Tfree.back() * cFree + m_TbndFine * cBnd); }

    /// Set the number of pre-/post-smoothing steps (helps for the badly
    /// conditioned biharmonic operator).
    void setNumSmooth(index_t pre, index_t post)
    {
        m_mg->setNumPreSmooth(pre);
        m_mg->setNumPostSmooth(post);
    }

    /// Set the number of coarse-grid cycles (1 = V-cycle, 2 = W-cycle).
    void setNumCycles(index_t n) { if (m_numLevels > 1) m_mg->setNumCycles(n); }

private:

    // Reconstruct a multi-patch field from a disjoint (per-patch) coefficient vector.
    gsMultiPatch<T> reconstructFromDisjoint(const gsMatrix<T>& cDisjoint) const
    {
        const gsMultiPatch<T>& mp = m_mp.back();
        gsMultiPatch<T> field;
        index_t offset = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            const index_t sz = m_patchDisjointRows.back()[i];
            gsMatrix<T> ci = cDisjoint.block(offset, 0, sz, 1);
            offset += sz;
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(mp.patch(i).basis());
            field.addPatch(tb.makeGeometry(give(ci)));
        }
        return field;
    }

    // Horizontal concatenation of two sparse matrices with equal row count.
    static gsSparseMatrix<T> hconcat(const gsSparseMatrix<T>& A, const gsSparseMatrix<T>& B)
    {
        GISMO_ASSERT(A.rows() == B.rows(), "hconcat: row size mismatch.");
        gsSparseEntries<T> e;
        for (index_t k = 0; k < A.outerSize(); ++k)
            for (typename gsSparseMatrix<T>::InnerIterator it(A, k); it; ++it)
                e.add(it.row(), it.col(), it.value());
        for (index_t k = 0; k < B.outerSize(); ++k)
            for (typename gsSparseMatrix<T>::InnerIterator it(B, k); it; ++it)
                e.add(it.row(), A.cols() + it.col(), it.value());
        gsSparseMatrix<T> C(A.rows(), A.cols() + B.cols());
        C.setFrom(e);
        return C;
    }

    // Build the coupled AS-G1 space of one level. Returns the transformation for
    // the free (interior) DOFs; if @a Tbnd is non-null it also receives the
    // transformation for the eliminated boundary DOFs.
    gsSparseMatrix<T> buildLevelTransform(const gsMultiPatch<T>& mp,
                                          const gsMatrix<T>& gd,
                                          std::vector<index_t>& patchRows,
                                          gsSparseMatrix<T>* Tbnd = nullptr) const
    {
        // Per-patch Argyris embeddings.
        std::vector<gsArgyrisEmbedding<T>> argBasis;
        argBasis.reserve(mp.nPatches());
        patchRows.assign(mp.nPatches(), 0);
        index_t nDisjoint = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            argBasis.push_back(deriveArgyrisBasisEmbedding(
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(mp.patch(i).basis()),
                gsMatrix<T>(gd.row(i)),
                const_cast<gsGeometry<T>&>(mp.patch(i))));
            patchRows[i] = argBasis[i].matrix.rows();
            nDisjoint += patchRows[i];
        }

        // Mark all external boundary sides as clamped (u and grad(u).n given).
        std::set<std::pair<size_t, index_t>> ifcSides;
        for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
        {
            ifcSides.insert({it->first().patch,  it->first().side()});
            ifcSides.insert({it->second().patch, it->second().side()});
        }
        gsBoundaryConditions<T> bc;
        typename gsBoundaryConditions<T>::function_ptr nullFun;
        for (size_t i = 0; i < mp.nPatches(); ++i)
            for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
                if (ifcSides.find({i, side.m_index}) == ifcSides.end())
                    bc.add(i, side, "ValuesAndDerivatives", nullFun);

        // Interface coupling + boundary elimination via the shared helper.
        gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis, bc);
        const index_t nFree = mapper.freeSize();
        const index_t nBnd  = mapper.boundarySize();

        gsSparseEntries<T> ef, eb;
        index_t rowOffset = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            const gsSparseMatrix<T>& Ai = argBasis[i].matrix;
            const index_t patchSize = mapper.patchSize(i);
            for (index_t j = 0; j < patchSize; ++j)
            {
                const bool isB = mapper.is_boundary(j, i);
                const index_t g = isB ? mapper.bindex(j, i) : mapper.index(j, i);
                for (typename gsSparseMatrix<T>::InnerIterator it(Ai, j); it; ++it)
                {
                    if (isB) eb.add(rowOffset + it.row(), g, it.value());
                    else     ef.add(rowOffset + it.row(), g, it.value());
                }
            }
            rowOffset += Ai.rows();
        }

        gsSparseMatrix<T> Tfree(nDisjoint, nFree);
        Tfree.setFrom(ef);
        if (Tbnd)
        {
            Tbnd->resize(nDisjoint, nBnd);
            Tbnd->setFrom(eb);
        }
        return Tfree;
    }

    // Assemble the disjoint (block-diagonal, per patch) biharmonic stiffness
    // matrix for one level.
    gsSparseMatrix<T> assembleDisjointStiffness(const gsMultiPatch<T>& mp) const
    {
        using namespace gismo::expr;
        gsMultiBasis<T> dbasis(mp);
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(dbasis.domain());
        auto G = A.getMap(mp);
        auto u = A.getSpace(dbasis);
        A.initSystem();
        A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G));
        return A.matrix();
    }

    // Per-patch B-spline refinement transfer (fine <- coarse), assembled into a
    // block-diagonal disjoint transfer.
    gsSparseMatrix<T> buildRefinement(const gsMultiPatch<T>& coarse,
                                      index_t mult) const
    {
        gsSparseEntries<T> entries;
        index_t rowOff = 0, colOff = 0;
        for (size_t p = 0; p < coarse.nPatches(); ++p)
        {
            gsTensorBSplineBasis<2, T> basis =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(coarse.patch(p).basis());
            gsSparseMatrix<T, RowMajor> Rp;
            basis.uniformRefine_withTransfer(Rp, 1, mult);
            for (int k = 0; k < Rp.outerSize(); ++k)
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(Rp, k); it; ++it)
                    entries.add(rowOff + it.row(), colOff + it.col(), it.value());
            rowOff += Rp.rows();
            colOff += Rp.cols();
        }
        gsSparseMatrix<T> R(rowOff, colOff);
        R.setFrom(entries);
        return R;
    }

    void setup(const gsMultiPatch<T>& geometry)
    {
        m_mp.resize(m_numLevels);
        m_Tfree.resize(m_numLevels);
        m_patchDisjointRows.resize(m_numLevels);

        gsMultiPatch<T> base = geometry;
        base.computeTopology();
        const gsMatrix<T> gd = computeGluingData(base, T(1e-8), m_numGauss);

        // Build each level (coarse -> fine).
        std::vector<std::shared_ptr<typename gsSparseSolver<T>::SimplicialLDLT>> stsSolvers(m_numLevels);
        for (index_t l = 0; l < m_numLevels; ++l)
        {
            gsMultiPatch<T> mp = base;
            const short_t inputDeg = mp.patch(0).basis().degree(0);
            if (inputDeg < m_degree) mp.degreeElevate(m_degree - inputDeg);
            const short_t deg  = mp.patch(0).basis().degree(0);
            const index_t mult = std::max<index_t>(deg - 1, 1);
            const index_t ref  = m_minRef + l;
            for (index_t r = 0; r < ref; ++r) mp.uniformRefine(1, mult);

            m_mp[l] = mp;
            m_Tfree[l] = buildLevelTransform(mp, gd, m_patchDisjointRows[l],
                                             (l == m_numLevels - 1) ? &m_TbndFine : nullptr);

            // Normal-equation solver for the exact prolongation.
            gsSparseMatrix<T> StS = m_Tfree[l].transpose() * m_Tfree[l];
            stsSolvers[l] = std::make_shared<typename gsSparseSolver<T>::SimplicialLDLT>();
            stsSolvers[l]->compute(StS);
        }

        // Finest-level stiffness matrix and interior/boundary coupling block.
        {
            gsSparseMatrix<T> Kdis = assembleDisjointStiffness(m_mp.back());
            m_Kfine = m_Tfree.back().transpose() * Kdis * m_Tfree.back();
            m_Kfine.makeCompressed();
            m_Kfb = m_Tfree.back().transpose() * Kdis * m_TbndFine;
            m_Kfb.makeCompressed();
        }

        // Inter-level prolongations P_l : level l -> level l+1 (fine <- coarse).
        std::vector<SpTransfer> transfers(m_numLevels - 1);
        for (index_t l = 0; l + 1 < m_numLevels; ++l)
        {
            const short_t deg  = m_mp[l].patch(0).basis().degree(0);
            const index_t mult = std::max<index_t>(deg - 1, 1);
            gsSparseMatrix<T> R = buildRefinement(m_mp[l], mult);

            // M = T_fine^T * R * T_coarse ; P = (T_fine^T T_fine)^{-1} M.
            gsSparseMatrix<T> M = m_Tfree[l + 1].transpose() * (R * m_Tfree[l]);
            gsMatrix<T> Pdense = stsSolvers[l + 1]->solve(gsMatrix<T>(M));
            transfers[l] = SpTransfer(Pdense.sparseView(T(1), T(1e-10)));
        }

        // Assemble the multigrid preconditioner (Galerkin coarsening).
        m_mg = gsMultiGridOp<T>::make(m_Kfine, transfers);
        m_mg->setNumPreSmooth(3);
        m_mg->setNumPostSmooth(3);
        m_mg->setSymmSmooth(true);
        if (m_numLevels > 1) m_mg->setNumCycles(1);

        // Symmetric Gauss-Seidel smoothers on every level.
        for (index_t l = 0; l < m_numLevels; ++l)
        {
            std::shared_ptr<gsSparseMatrix<T>> Kl =
                std::make_shared<gsSparseMatrix<T>>(m_mg->matrix(l));
            m_mg->setSmoother(l, makeSymmetricGaussSeidelOp(Kl));
        }

        m_mgPtr = typename gsLinearOperator<T>::Ptr(m_mg);
    }

    index_t m_degree, m_minRef, m_maxRef, m_numGauss, m_numLevels;

    std::vector<gsMultiPatch<T>>    m_mp;                 // per-level geometry
    std::vector<gsSparseMatrix<T>>  m_Tfree;             // per-level transform
    std::vector<std::vector<index_t>> m_patchDisjointRows;
    gsSparseMatrix<T>               m_Kfine;             // finest stiffness
    gsSparseMatrix<T>               m_TbndFine;          // finest boundary transform
    gsSparseMatrix<T>               m_Kfb;               // interior-boundary coupling

    typename gsMultiGridOp<T>::Ptr  m_mg;
    typename gsLinearOperator<T>::Ptr m_mgPtr;
};

} // namespace gismo
