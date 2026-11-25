/** @file gsAsG1ShellSolver.hpp

    @brief Multi-level geometric multigrid & single-level multi-patch domain
           decomposition solver hierarchy for Kirchhoff-Love shells on AS-G1
           surfaces (Hasanova, S. Takacs, T. Takacs 2026).

    Supports:
      - Arbitrary degree p >= 3 and interior smoothness 1 <= k <= p - 2
        with strictly nested dyadic refinement uniformRefine(1, p - k),
      - Non-isoparametric polynomial spaces S_{p,k,Z} over polynomial OR
        exact rational AS-G1 NURBS surfaces (Paper Section 3.1),
      - Appendix A normalized gluing data & shared 3D vertex tangent frames,
      - Component-block exact sparse prolongation P_{l-1 -> l},
      - Direct (SimplicialLDLT), Single-level Schwarz-PCG, V/W-cycle MG, and MG-PCG.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsModeling/gsAsG1Domain.hpp>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1KLShellAssembler.hpp>
#include <gsModeling/gsAsG1ShellSmoother.hpp>
#include <gsMultiGrid/gsMultiGrid.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsSolver/gsSimplePreconditioners.h>
#include <gsUtils/gsStopwatch.h>

#include <cmath>
#include <string>
#include <vector>

namespace gismo {

template <typename T>
struct gsAsG1ShellLevel
{
    gsMultiBasis<T> dbasis;
    std::vector<gsArgyrisEmbedding<T>> argBasis;         // derivative-scaled embeddings
    std::vector<gsArgyrisEmbedding<T>> unscaledArgBasis; // populated on finest level for exact projection
    std::vector<gsVector<T>> colScales;
    gsDofMapper mapper;
    gsSparseMatrix<T> T_free;
    gsSparseMatrix<T> T_bnd;
    T hLevel = T(1);
};

template <typename T>
struct gsAsG1ShellSolveStats
{
    index_t nFree = 0;
    index_t nBnd = 0;
    index_t iterations = 0;
    T relResidual = T(0);
    double setupTimeSec = 0.0;
    double solveTimeSec = 0.0;
    bool converged = false;
};

/**
 * @brief Nested multi-level AS-G1 Kirchhoff-Love shell solver hierarchy.
 */
template <typename T>
class gsAsG1ShellHierarchy
{
public:
    gsAsG1ShellHierarchy(const gsMultiPatch<T>& mp,
                         const gsBoundaryConditions<T>& bc,
                         short_t degreeP,
                         short_t regularityK,
                         index_t minRef,
                         index_t maxRef,
                         bool normalizeGluing = true,
                         T gluingEps = T(1e-5))
        : m_mp(mp), m_bc(bc), m_p(degreeP),
          m_k(regularityK < 0 ? static_cast<short_t>(degreeP - 2) : regularityK),
          m_minRef(minRef), m_maxRef(maxRef)
    {
        GISMO_ENSURE(m_mp.targetDim() == 3, "gsAsG1ShellHierarchy expects a 3D surface (targetDim == 3).");
        GISMO_ENSURE(m_p >= 3 && m_k >= 1 && m_k <= m_p - 2,
                     "Require p >= 3 and 1 <= k <= p - 2.");
        GISMO_ENSURE(m_maxRef >= m_minRef, "Require maxRef >= minRef.");

        // 1. Compute 3D surface gluing data (with optional Appendix A normalization)
        m_gluing = computeGluingData(m_mp, gluingEps);
        if (normalizeGluing)
        {
            auto sideCol = [](const patchSide& ps) -> index_t {
                return 4 * (static_cast<index_t>(ps.side()) - 1);
            };
            for (auto it = m_mp.iBegin(); it != m_mp.iEnd(); ++it)
            {
                const index_t p1 = it->first().patch, p2 = it->second().patch;
                const index_t c1 = sideCol(it->first()), c2 = sideCol(it->second());
                gsVector<T> gd(8);
                gd.head(4) = m_gluing.row(p1).segment(c1, 4).transpose();
                gd.tail(4) = m_gluing.row(p2).segment(c2, 4).transpose();
                normalizeGluingDataAppendixA(gd);
                m_gluing.row(p1).segment(c1, 4) = gd.head(4).transpose();
                m_gluing.row(p2).segment(c2, 4) = gd.tail(4).transpose();
            }
        }

        // 2. Compute shared 3D vertex tangent frames (stride = 10) and corner conditions per component
        gsMatrix<index_t> cornerSmooth;
        m_cornerFrames = computeSurfaceCornerFrames(m_mp, &cornerSmooth, T(1e-3));
        m_ccPerComp = buildCornerConditionsPerComponent(m_mp, m_bc, cornerSmooth, m_cornerFrames);

        // 3. Build nested (p, k) levels l = 0 .. L (L = maxRef - minRef)
        const index_t nLevels = m_maxRef - m_minRef + 1;
        m_levels.resize(nLevels);
        const short_t mult = m_p - m_k;

        for (index_t l = 0; l < nLevels; ++l)
        {
            if (l == 0)
            {
                m_levels[0].dbasis = buildAsG1PolynomialMultiBasis(m_mp, m_p, m_k, m_minRef);
            }
            else
            {
                m_levels[l].dbasis = m_levels[l - 1].dbasis;
                m_levels[l].dbasis.uniformRefine(1, mult);
            }
            buildLevelData(l, l == nLevels - 1);
        }

        // 4. Compute exact component-block sparse intergrid prolongations P_{l -> l+1}
        m_transfers.resize(nLevels > 1 ? nLevels - 1 : 0);
        for (index_t l = 0; l + 1 < nLevels; ++l)
            m_transfers[l] = computeVectorProlongation(m_levels[l], m_levels[l + 1], mult);
    }

    index_t numLevels() const { return static_cast<index_t>(m_levels.size()); }
    const gsAsG1ShellLevel<T>& level(index_t l) const { return m_levels[l]; }
    const gsAsG1ShellLevel<T>& finest() const { return m_levels.back(); }
    const gsMatrix<T>& gluingData() const { return m_gluing; }
    const gsMatrix<T>& cornerFrames() const { return m_cornerFrames; }
    short_t degree() const { return m_p; }
    short_t regularity() const { return m_k; }

    /**
     * @brief Solve K_free * u_free = F_free using one of:
     *        - "direct"     : SimplicialLDLT sparse direct solver
     *        - "schwarz_cg" : Single-level multi-patch Patch-Interior + Skeleton DD-PCG
     *        - "mg"         : Stationary V-cycle / W-cycle geometric multigrid
     *        - "cg"         : Multigrid-Preconditioned Conjugate Gradient (MG-PCG)
     *
     * Smoother choices for "mg" / "cg":
     *        "nodal_gs", "ring1", "line_ring1", "scms", "scms_hybrid"
     */
    gsAsG1ShellSolveStats<T> solve(
        const gsSparseMatrix<T>& K_free,
        const gsVector<T>& F_free,
        const gsKLShellMaterial<T>& mat,
        gsVector<T>& u_free,
        const std::string& solverType = "cg",
        const std::string& smootherType = "line_ring1",
        index_t numCycles = 1,
        index_t preSmooth = 1,
        index_t postSmooth = 1,
        index_t maxIter = 200,
        T tol = T(1e-10),
        index_t scmsRp = -1) const
    {
        GISMO_ENSURE(solverType == "direct" || solverType == "schwarz_cg" ||
                     solverType == "mg" || solverType == "cg",
                     "Unknown shell solver type: " << solverType
                     << " (expected \"direct\", \"schwarz_cg\", \"mg\", or \"cg\").");

        gsAsG1ShellSolveStats<T> stats;
        stats.nFree = finest().mapper.freeSize();
        stats.nBnd  = finest().mapper.boundarySize();

        gsStopwatch timer;
        const T fNorm = F_free.norm();
        u_free.setZero(K_free.rows());

        if (solverType == "direct" || (m_levels.size() == 1 && (solverType == "mg" || solverType == "cg")))
        {
            timer.restart();
            gsEigen::SimplicialLDLT<gsSparseMatrix<T>> ldlt;
            ldlt.compute(K_free);
            stats.setupTimeSec = timer.stop();

            timer.restart();
            u_free = ldlt.solve(F_free);
            stats.solveTimeSec = timer.stop();
            stats.iterations = 1;
            stats.relResidual = (fNorm > T(0)) ? (F_free - K_free * u_free).norm() / fNorm : T(0);
            stats.converged = (ldlt.info() == gsEigen::Success) && (stats.relResidual < T(1e-6));
            return stats;
        }

        if (solverType == "schwarz_cg")
        {
            timer.restart();
            auto prec = gsAsG1ShellPatchSchwarzPrecOp<T>::make(
                K_free, m_mp, finest().dbasis, finest().argBasis, finest().mapper, scmsRp);
            stats.setupTimeSec = timer.stop();

            timer.restart();
            gsConjugateGradient<T> cg(K_free, typename gsLinearOperator<T>::Ptr(prec.release()));
            cg.setTolerance(tol);
            cg.setMaxIterations(maxIter);
            gsMatrix<T> xMat = gsMatrix<T>::Zero(K_free.rows(), 1);
            gsMatrix<T> bMat = F_free;
            cg.solve(bMat, xMat);
            u_free = xMat.col(0);
            stats.solveTimeSec = timer.stop();
            stats.iterations = cg.iterations();
            stats.relResidual = (fNorm > T(0)) ? (F_free - K_free * u_free).norm() / fNorm : T(0);
            stats.converged = (stats.relResidual <= tol * T(10));
            return stats;
        }

        // Multi-level Geometric Multigrid (V/W-cycle or MG-PCG)
        timer.restart();
        typename gsMultiGridOp<T>::uPtr mg = gsMultiGridOp<T>::make(K_free, m_transfers);
        mg->setNumCycles(numCycles);
        mg->setNumPreSmooth(preSmooth);
        mg->setNumPostSmooth(postSmooth);
        mg->setCoarseSolver(makeSparseCholeskySolver(mg->matrix(0)));

        const index_t nL = static_cast<index_t>(m_levels.size());
        for (index_t l = 1; l < nL; ++l)
        {
            typename gsMultiGridOp<T>::PrecondPtr sm;
            if (smootherType == "nodal_gs")
                sm = typename gsMultiGridOp<T>::PrecondPtr(
                    gsNodalBlockGaussSeidelOp<T>::make(
                        mg->matrix(l), m_levels[l].argBasis, m_levels[l].mapper, 3).release());
            else if (smootherType == "ring1")
                sm = typename gsMultiGridOp<T>::PrecondPtr(
                    gsAsG1ShellSkeletonSmootherOp<T>::make(
                        mg->matrix(l), m_mp, m_levels[l].dbasis, m_levels[l].argBasis,
                        m_levels[l].mapper, 1).release());
            else if (smootherType == "line_ring1")
                sm = typename gsMultiGridOp<T>::PrecondPtr(
                    gsAsG1ShellBlockSmootherOp<T>::make(
                        mg->matrix(l), m_mp, m_levels[l].dbasis, m_levels[l].argBasis,
                        m_levels[l].mapper, 1).release());
            else if (smootherType == "scms")
                sm = typename gsMultiGridOp<T>::PrecondPtr(
                    gsAsG1ShellScmsSmootherOp<T>::make(
                        mg->matrix(l), m_mp, m_levels[l].dbasis, m_levels[l].argBasis,
                        m_levels[l].mapper, mat.E, mat.nu, mat.thickness, scmsRp, false).release());
            else if (smootherType == "scms_hybrid")
                sm = typename gsMultiGridOp<T>::PrecondPtr(
                    gsAsG1ShellScmsSmootherOp<T>::make(
                        mg->matrix(l), m_mp, m_levels[l].dbasis, m_levels[l].argBasis,
                        m_levels[l].mapper, mat.E, mat.nu, mat.thickness, scmsRp, true).release());
            else
                GISMO_ERROR("Unknown shell smoother type: " << smootherType);
            mg->setSmoother(l, sm);
        }
        stats.setupTimeSec = timer.stop();

        timer.restart();
        gsMatrix<T> xMat = gsMatrix<T>::Zero(K_free.rows(), 1);
        gsMatrix<T> bMat = F_free;
        if (solverType == "mg")
        {
            index_t it = 0;
            for (it = 1; it <= maxIter; ++it)
            {
                mg->step(bMat, xMat);
                const T rel = (fNorm > T(0)) ? (F_free - K_free * xMat.col(0)).norm() / fNorm : T(0);
                if (rel <= tol)
                    break;
            }
            stats.iterations = std::min(it, maxIter);
        }
        else if (solverType == "cg")
        {
            gsConjugateGradient<T> cg(K_free, typename gsLinearOperator<T>::Ptr(mg.release()));
            cg.setTolerance(tol);
            cg.setMaxIterations(maxIter);
            cg.solve(bMat, xMat);
            stats.iterations = cg.iterations();
        }
        else
        {
            GISMO_ERROR("Unknown shell solver type: " << solverType);
        }
        u_free = xMat.col(0);
        stats.solveTimeSec = timer.stop();
        stats.relResidual = (fNorm > T(0)) ? (F_free - K_free * u_free).norm() / fNorm : T(0);
        stats.converged = (stats.relResidual <= tol * T(10));
        return stats;
    }

private:
    static std::vector<cornerConditions> buildCornerConditionsPerComponent(
        const gsMultiPatch<T>& mp,
        const gsBoundaryConditions<T>& bc,
        const gsMatrix<index_t>& cornerSmooth,
        const gsMatrix<T>& cornerFrames)
    {
        std::vector<cornerConditions> ccPerComp(3);
        const auto& vertices = mp.topology().vertices();
        const index_t stride = cornerFrames.cols() / 4;

        for (index_t comp = 0; comp < 3; ++comp)
        {
            auto hasSideLabel = [&](const patchSide& ps, const std::string& label) -> bool {
                for (auto it = bc.begin(label); it != bc.end(label); ++it)
                    if (it->unknown() == comp && it->ps == ps)
                        return true;
                return false;
            };
            auto sideHasV = [&](const patchSide& ps) -> bool {
                return hasSideLabel(ps, "ValuesAndDerivatives")
                    || hasSideLabel(ps, "Clamped")
                    || hasSideLabel(ps, "Values")
                    || hasSideLabel(ps, "Dirichlet");
            };
            auto sideHasD = [&](const patchSide& ps) -> bool {
                return hasSideLabel(ps, "ValuesAndDerivatives")
                    || hasSideLabel(ps, "Clamped")
                    || hasSideLabel(ps, "Derivatives");
            };

            for (size_t v = 0; v < vertices.size(); ++v)
            {
                const auto& grp = vertices[v];
                if (grp.empty()) continue;

                // Enumerate incident boundary sides in the exact same order as
                // getSurfaceBoundaryCornerData so bSides[0] aligns with frame dir1
                // and bSides[1] aligns with frame dir2.
                std::vector<patchSide> bSides;
                for (size_t i = 0; i < grp.size(); ++i)
                {
                    std::vector<patchSide> sides;
                    grp[i].getContainingSides(2, sides);
                    for (size_t j = 0; j < sides.size(); ++j)
                        if (mp.isBoundary(sides[j]))
                            bSides.push_back(sides[j]);
                }
                if (bSides.empty()) continue;

                const index_t p0 = grp[0].patch;
                const index_t c0 = grp[0].m_index - 1;
                const bool isSmooth = (cornerSmooth(p0, c0) != 0);

                cornerConditionType ctype = cornerConditionType::none;
                if (bSides.size() == 2)
                {
                    const bool v0 = sideHasV(bSides[0]), d0 = sideHasD(bSides[0]);
                    const bool v1 = sideHasV(bSides[1]), d1 = sideHasD(bSides[1]);

                    if (isSmooth)
                    {
                        // Case B: smooth boundary vertex, frame [d1 = n_co, d2 = t]
                        // Clamped (V+D) fixes u, d_n, d_t, d_tt + k*d_n, d_nt and MUST leave d_nn (j=3) free.
                        if ((v0 || v1) && (d0 || d1))
                            ctype = cornerConditionType::valuesNormals;
                        else if (v0 || v1)
                            ctype = cornerConditionType::valuesSmooth;
                        else if (d0 || d1)
                            ctype = cornerConditionType::derivatives;
                    }
                    else
                    {
                        // Case A: sharp boundary corner, frame [d1 = t_0 (bSides[0]), d2 = t_1 (bSides[1])]
                        gsVector<T, 3> t0, t1;
                        t0 << cornerFrames(p0, stride * c0 + 0),
                              cornerFrames(p0, stride * c0 + 1),
                              cornerFrames(p0, stride * c0 + 2);
                        t1 << cornerFrames(p0, stride * c0 + 3),
                              cornerFrames(p0, stride * c0 + 4),
                              cornerFrames(p0, stride * c0 + 5);
                        const bool isOrthogonal = (std::abs(t0.dot(t1)) < T(1e-3));

                        if (v0 && v1)
                        {
                            // Both edges constrain values; any normal-slope constraint eliminates all 6 jet DOFs.
                            ctype = (d0 || d1) ? cornerConditionType::all
                                               : cornerConditionType::valuesSharp;
                        }
                        else if (v0 && d0)
                        {
                            // Side 0 clamped, side 1 free (or orthogonal symmetry): leave d^2/ds_1^2 (j=4) free.
                            ctype = (!d1 || isOrthogonal) ? cornerConditionType::clampedSide1
                                                          : cornerConditionType::all;
                        }
                        else if (v1 && d1)
                        {
                            // Side 1 clamped, side 0 free (or orthogonal symmetry): leave d^2/ds_0^2 (j=3) free.
                            ctype = (!d0 || isOrthogonal) ? cornerConditionType::clampedSide2
                                                          : cornerConditionType::all;
                        }
                        else if (v0)
                        {
                            // Side 0 values only, side 1 free or symmetry derivative
                            ctype = d1 ? cornerConditionType::valuesSide1DerSide2
                                       : cornerConditionType::valuesSide1;
                        }
                        else if (v1)
                        {
                            // Side 1 values only, side 0 free or symmetry derivative
                            ctype = d0 ? cornerConditionType::derSide1ValuesSide2
                                       : cornerConditionType::valuesSide2;
                        }
                        else if (d0 && d1)
                        {
                            ctype = cornerConditionType::derivativesSharp;
                        }
                        else if (d0)
                        {
                            ctype = cornerConditionType::derivativesSide1;
                        }
                        else if (d1)
                        {
                            ctype = cornerConditionType::derivativesSide2;
                        }
                    }
                }
                else
                {
                    // Fallback for non-manifold boundary vertices
                    index_t nVal = 0, nDer = 0;
                    for (size_t k = 0; k < bSides.size(); ++k)
                    {
                        if (sideHasV(bSides[k])) ++nVal;
                        if (sideHasD(bSides[k])) ++nDer;
                    }
                    if (nVal > 0 && nDer > 0)
                        ctype = cornerConditionType::all;
                    else if (nVal >= 2)
                        ctype = cornerConditionType::valuesSharp;
                    else if (nVal == 1)
                        ctype = cornerConditionType::valuesSmooth;
                    else if (nDer > 0)
                        ctype = cornerConditionType::derivatives;
                }

                if (ctype != cornerConditionType::none)
                    ccPerComp[comp].emplace_back(grp[0], ctype);
            }
        }
        return ccPerComp;
    }

    void buildLevelData(index_t l, bool isFinest)
    {
        gsAsG1ShellLevel<T>& ld = m_levels[l];
        const size_t nP = m_mp.nPatches();
        ld.argBasis.resize(nP);

        // Compute representative parametric mesh size hLevel as maximum knot span across all patches & directions
        T hMax = T(0);
        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tbp =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(ld.dbasis[p]);
            for (short_t d = 0; d < 2; ++d)
            {
                const std::vector<T> brks = tbp.component(d).knots().breaks();
                for (size_t e = 0; e + 1 < brks.size(); ++e)
                    hMax = std::max(hMax, brks[e + 1] - brks[e]);
            }
        }
        ld.hLevel = (hMax > T(0)) ? hMax : T(1);

        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(ld.dbasis[p]);
            gsGeometry<T>& geoRef = const_cast<gsGeometry<T>&>(m_mp.patch(p));
            ld.argBasis[p] = deriveArgyrisBasisEmbedding(
                tb, gsMatrix<T>(m_gluing.row(p)), gsMatrix<T>(m_cornerFrames.row(p)), geoRef);
        }

        if (isFinest)
            ld.unscaledArgBasis = ld.argBasis;

        ld.colScales = scaleArgyrisDerivativeBlocks(ld.argBasis, ld.hLevel);
        ld.mapper = makeVectorMapperForArgyrisBasis(m_mp, ld.argBasis, 3, m_bc, m_ccPerComp);
        buildVectorArgyrisTransferMatrices(ld.argBasis, ld.mapper, 3, ld.T_free, ld.T_bnd);
    }

    /// Compute exact component-block sparse prolongation P_{c -> f} for all 3 components.
    gsSparseMatrix<T, RowMajor> computeVectorProlongation(
        const gsAsG1ShellLevel<T>& coarse,
        const gsAsG1ShellLevel<T>& fine,
        short_t mult) const
    {
        const size_t nP = m_mp.nPatches();
        gsBlockSparseMatrix<T> R_block(nP, nP);
        for (size_t p = 0; p < nP; ++p)
        {
            gsSparseMatrix<T, RowMajor> localTr;
            const gsTensorBSplineBasis<2, T>& cBasis =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(coarse.dbasis[p]);
            gsTensorBSplineBasis<2, T> copyBasis(cBasis);
            copyBasis.uniformRefine_withTransfer(localTr, 1, mult);
            gsSparseMatrix<T> localTrCol = localTr;
            R_block.set(p, p, localTrCol);
        }
        const gsSparseMatrix<T> R_scalar = R_block;
        const index_t nDiscontCoarse = R_scalar.cols();
        const index_t nDiscontFine   = R_scalar.rows();

        const index_t nFreeFine   = fine.mapper.freeSize();
        const index_t nFreeCoarse = coarse.mapper.freeSize();

        // Extract per-component slices of T_free (since gsDofMapper orders free DOFs contiguously per component)
        gsSparseEntries<T> allEntries;
        index_t fColOff = 0, cColOff = 0;

        for (index_t comp = 0; comp < 3; ++comp)
        {
            const index_t nFc = fine.mapper.freeSize(comp);
            const index_t nCc = coarse.mapper.freeSize(comp);
            if (nFc == 0 || nCc == 0)
            {
                fColOff += nFc;
                cColOff += nCc;
                continue;
            }

            const gsSparseMatrix<T> Tf_c = fine.T_free.block(comp * nDiscontFine, fColOff, nDiscontFine, nFc);
            const gsSparseMatrix<T> Tc_c = coarse.T_free.block(comp * nDiscontCoarse, cColOff, nDiscontCoarse, nCc);

            const gsSparseMatrix<T> M_gram = Tf_c.transpose() * Tf_c;
            const gsSparseMatrix<T> RHS_mat = Tf_c.transpose() * (R_scalar * Tc_c);

            gsEigen::SimplicialLDLT<gsSparseMatrix<T>> chol;
            chol.compute(M_gram);

            const index_t batchSize = 128;
            const index_t nBatches = (nCc + batchSize - 1) / batchSize;
            std::vector<gsSparseEntries<T>> batchEnt(nBatches);

            #pragma omp parallel for schedule(dynamic)
            for (index_t b = 0; b < nBatches; ++b)
            {
                const index_t bStart = b * batchSize;
                const index_t bCols = std::min(batchSize, nCc - bStart);
                gsMatrix<T> rhsB = gsMatrix<T>(RHS_mat.middleCols(bStart, bCols));
                gsMatrix<T> P_b = chol.solve(rhsB);
                for (index_t j = 0; j < bCols; ++j)
                {
                    for (index_t i = 0; i < nFc; ++i)
                    {
                        const T val = P_b(i, j);
                        if (std::abs(val) > T(1e-12))
                            batchEnt[b].add(fColOff + i, cColOff + bStart + j, val);
                    }
                }
            }
            for (index_t b = 0; b < nBatches; ++b)
                for (const auto& t : batchEnt[b])
                    allEntries.add(t.row(), t.col(), t.value());

            fColOff += nFc;
            cColOff += nCc;
        }

        gsSparseMatrix<T> P_col(nFreeFine, nFreeCoarse);
        P_col.setFrom(allEntries);
        gsSparseMatrix<T, RowMajor> P_row = P_col;
        return P_row;
    }

    gsMultiPatch<T> m_mp;
    gsBoundaryConditions<T> m_bc;
    short_t m_p;
    short_t m_k;
    index_t m_minRef;
    index_t m_maxRef;
    gsMatrix<T> m_gluing;
    gsMatrix<T> m_cornerFrames;
    std::vector<cornerConditions> m_ccPerComp;
    std::vector<gsAsG1ShellLevel<T>> m_levels;
    std::vector<gsSparseMatrix<T, RowMajor>> m_transfers;
};

} // namespace gismo
