/** @file gsAsG1Multigrid.hpp
 *
 *  @brief Implementation of Analysis-Suitable G1 (AS-G1) Multi-Patch Multigrid Solver.
 *
 *  This file is part of the G+Smo library.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Author(s): Antigravity AI, F. Hasanova, S. Takacs
 */

#pragma once

#include <gsMultiGrid/gsAsG1Multigrid.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

namespace gismo
{

template<typename T>
gsAsG1Multigrid<T>::gsAsG1Multigrid(const gsMultiPatch<T>& geometry,
                                     index_t degree,
                                     index_t minRefinements,
                                     index_t maxRefinements,
                                     index_t numGaussPerSpan)
    : m_baseGeometry(geometry),
      m_degree(degree < 3 ? 3 : degree),
      m_minRefinements(minRefinements < 2 ? 2 : minRefinements),
      m_maxRefinements(maxRefinements < minRefinements ? minRefinements : maxRefinements),
      m_numGaussPerSpan(numGaussPerSpan),
      m_numLevels(m_maxRefinements - m_minRefinements + 1)
{
    setupHierarchy();
}

template<typename T>
void gsAsG1Multigrid<T>::setupHierarchy()
{
    using namespace gismo::expr;

    m_mpLevels.resize(m_numLevels);
    m_S.resize(m_numLevels);
    m_R.resize(m_numLevels - 1);
    m_K.resize(m_numLevels);
    m_diagK.resize(m_numLevels);
    m_ifcDofs.resize(m_numLevels);
    m_StSSolvers.resize(m_numLevels);
    m_ifcBlockSolvers.resize(m_numLevels);
    m_coarseSolvers.resize(m_numLevels);
    m_nFree.resize(m_numLevels);
    m_nDisjoint.resize(m_numLevels);

    m_baseGeometry.computeTopology();
    gsMatrix<T> globalGd = computeGluingData(m_baseGeometry, T(1e-8), m_numGaussPerSpan);

    for (index_t l = 0; l < m_numLevels; ++l)
    {
        index_t ref = m_minRefinements + l;
        m_mpLevels[l] = m_baseGeometry;

        const short_t inputDeg = m_mpLevels[l].patch(0).basis().degree(0);
        if (inputDeg < m_degree)
        {
            m_mpLevels[l].degreeElevate(m_degree - inputDeg);
        }

        const short_t deg = m_mpLevels[l].patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);
        for (index_t r = 0; r < ref; ++r)
        {
            m_mpLevels[l].uniformRefine(1, mult);
        }

        std::vector<gsArgyrisEmbedding<T>> argBasis;
        gsVector<index_t> patchDofSizes(m_mpLevels[l].nPatches());
        for (size_t i = 0; i < m_mpLevels[l].nPatches(); ++i)
        {
            argBasis.push_back(deriveArgyrisBasisEmbedding(
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(m_mpLevels[l].patch(i).basis()),
                gsMatrix<T>(globalGd.row(i)), m_mpLevels[l].patch(i)));
            patchDofSizes[i] = argBasis[i].matrix.cols();
        }

        gsDofMapper mapper(patchDofSizes);

        for (auto it = m_mpLevels[l].iBegin(); it != m_mpLevels[l].iEnd(); ++it)
        {
            const boundaryInterface &ifc = *it;
            const patchSide ps1 = ifc.first();
            const patchSide ps2 = ifc.second();

            const index_t nLvl0 = argBasis[ps1.patch].sizes[1 + 2 * (ps1.m_index - 1)];
            const index_t offLvl0_1 = sumUntil(argBasis[ps1.patch].sizes, 1 + 2 * (ps1.m_index - 1));
            const index_t offLvl0_2 = sumUntil(argBasis[ps2.patch].sizes, 1 + 2 * (ps2.m_index - 1));

            const index_t nLvl1 = argBasis[ps1.patch].sizes[2 + 2 * (ps1.m_index - 1)];
            const index_t offLvl1_1 = sumUntil(argBasis[ps1.patch].sizes, 2 + 2 * (ps1.m_index - 1));
            const index_t offLvl1_2 = sumUntil(argBasis[ps2.patch].sizes, 2 + 2 * (ps2.m_index - 1));

            const short_t tanDir1 = 1 - ps1.direction();
            const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

            for (index_t j1 = 0; j1 < nLvl0; ++j1)
            {
                const index_t j2 = flipped ? nLvl0 - 1 - j1 : j1;
                mapper.matchDof(ps1.patch, offLvl0_1 + j1, ps2.patch, offLvl0_2 + j2);
            }
            for (index_t j1 = 0; j1 < nLvl1; ++j1)
            {
                const index_t j2 = flipped ? nLvl1 - 1 - j1 : j1;
                mapper.matchDof(ps1.patch, offLvl1_1 + j1, ps2.patch, offLvl1_2 + j2);
            }

            std::vector<patchCorner> corners;
            ps1.getContainedCorners(2, corners);
            for (index_t i = 0; i < 2; ++i)
            {
                const index_t c1 = corners[i].m_index - 1;
                const index_t c2 = ifc.mapCorner(corners[i]).m_index - 1;

                const index_t off_corner_1 = sumUntil(argBasis[ps1.patch].sizes, 9 + c1);
                const index_t off_corner_2 = sumUntil(argBasis[ps2.patch].sizes, 9 + c2);

                for (index_t k = 0; k < 6; ++k)
                    mapper.matchDof(ps1.patch, off_corner_1 + k, ps2.patch, off_corner_2 + k);
            }
        }

        std::set<std::pair<size_t, index_t>> ifcSides;
        for (auto it = m_mpLevels[l].iBegin(); it != m_mpLevels[l].iEnd(); ++it)
        {
            ifcSides.insert({it->first().patch, it->first().side()});
            ifcSides.insert({it->second().patch, it->second().side()});
        }

        for (size_t i = 0; i < m_mpLevels[l].nPatches(); ++i)
        {
            for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side)
            {
                if (ifcSides.find({i, side.m_index}) != ifcSides.end())
                    continue;

                const index_t nLvl0 = argBasis[i].sizes[1 + 2 * (side.m_index - 1)];
                const index_t offLvl0 = sumUntil(argBasis[i].sizes, 1 + 2 * (side.m_index - 1));
                for (index_t j = 0; j < nLvl0; ++j)
                    mapper.eliminateDof(offLvl0 + j, i);

                const index_t nLvl1 = argBasis[i].sizes[2 + 2 * (side.m_index - 1)];
                const index_t offLvl1 = sumUntil(argBasis[i].sizes, 2 + 2 * (side.m_index - 1));
                for (index_t j = 0; j < nLvl1; ++j)
                    mapper.eliminateDof(offLvl1 + j, i);

                patchSide ps(i, side);
                std::vector<patchCorner> corners;
                ps.getContainedCorners(2, corners);
                for (const auto &c : corners)
                {
                    const index_t cIdx = c.m_index - 1;
                    const index_t offCorner = sumUntil(argBasis[i].sizes, 9 + cIdx);
                    for (index_t k = 0; k < 6; ++k)
                        mapper.eliminateDof(offCorner + k, i);
                }
            }
        }

        mapper.finalize();

        m_nFree[l] = mapper.freeSize();
        m_nDisjoint[l] = 0;
        for (size_t i = 0; i < m_mpLevels[l].nPatches(); ++i)
            m_nDisjoint[l] += argBasis[i].matrix.rows();

        gsSparseEntries<T> tFreeEntries;
        index_t rowOffset = 0;

        for (size_t i = 0; i < m_mpLevels[l].nPatches(); ++i)
        {
            const gsSparseMatrix<T> &Ai = argBasis[i].matrix;
            const index_t patchSize = mapper.patchSize(i);

            for (index_t j = 0; j < patchSize; ++j)
            {
                if (mapper.is_boundary(j, i)) continue;
                index_t gIdx = mapper.index(j, i);

                for (typename gsSparseMatrix<T>::InnerIterator it(Ai, j); it; ++it)
                {
                    tFreeEntries.add(rowOffset + it.row(), gIdx, it.value());
                }
            }
            rowOffset += Ai.rows();
        }

        m_S[l].resize(m_nDisjoint[l], m_nFree[l]);
        m_S[l].setFrom(tFreeEntries);

        gsSparseMatrix<T> StS = m_S[l].transpose() * m_S[l];
        m_StSSolvers[l] = std::make_shared<typename gsSparseSolver<T>::SimplicialLDLT>();
        m_StSSolvers[l]->compute(StS);

        // Classify Interface DOFs (columns of StS with non-zeros > 1 or diagonal != 1.0)
        for (index_t col = 0; col < StS.cols(); ++col)
        {
            index_t colNnz = 0;
            T diagVal = 0;
            for (typename gsSparseMatrix<T>::InnerIterator it(StS, col); it; ++it)
            {
                colNnz++;
                if (it.row() == col) diagVal = it.value();
            }
            if (colNnz > 1 || std::abs(diagVal - 1.0) > 1e-10)
            {
                m_ifcDofs[l].push_back(col);
            }
        }

        // Assemble Disjoint Biharmonic Stiffness Matrix K_disjoint
        gsMultiBasis<T> dbasis(m_mpLevels[l]);
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(dbasis.domain());

        auto G_map = A.getMap(m_mpLevels[l]);
        auto u_space = A.getSpace(dbasis);

        A.initSystem();
        A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map));

        const gsSparseMatrix<T> &K_disjoint = A.matrix();
        m_K[l] = m_S[l].transpose() * K_disjoint * m_S[l];

        m_diagK[l].resize(m_nFree[l]);
        for (index_t j = 0; j < m_nFree[l]; ++j)
        {
            m_diagK[l][j] = m_K[l].coeff(j, j);
            if (std::abs(m_diagK[l][j]) < 1e-14) m_diagK[l][j] = 1.0;
        }

        // Extract and factorize interface submatrix K_ifc
        index_t nIfc = m_ifcDofs[l].size();
        if (nIfc > 0)
        {
            gsSparseEntries<T> ifcEntries;
            for (index_t i = 0; i < nIfc; ++i)
            {
                index_t gi = m_ifcDofs[l][i];
                for (typename gsSparseMatrix<T>::InnerIterator it(m_K[l], gi); it; ++it)
                {
                    index_t gj = it.row();
                    auto locIt = std::find(m_ifcDofs[l].begin(), m_ifcDofs[l].end(), gj);
                    if (locIt != m_ifcDofs[l].end())
                    {
                        index_t j = std::distance(m_ifcDofs[l].begin(), locIt);
                        ifcEntries.add(j, i, it.value());
                    }
                }
            }
            gsSparseMatrix<T> K_ifc(nIfc, nIfc);
            K_ifc.setFrom(ifcEntries);

            m_ifcBlockSolvers[l] = std::make_shared<typename gsSparseSolver<T>::SimplicialLDLT>();
            m_ifcBlockSolvers[l]->compute(K_ifc);
        }

        if (l == 0)
        {
            m_coarseSolvers[0] = std::make_shared<typename gsSparseSolver<T>::SimplicialLDLT>();
            m_coarseSolvers[0]->compute(m_K[0]);
        }
    }

    for (index_t l = 1; l < m_numLevels; ++l)
    {
        const short_t deg = m_mpLevels[l - 1].patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);

        gsSparseEntries<T> rEntries;
        index_t rRowOffset = 0;
        index_t rColOffset = 0;

        for (size_t p = 0; p < m_mpLevels[l - 1].nPatches(); ++p)
        {
            gsTensorBSplineBasis<2, T> coarseBasis = dynamic_cast<const gsTensorBSplineBasis<2, T> &>(m_mpLevels[l - 1].patch(p).basis());
            gsSparseMatrix<T, RowMajor> R_patch;
            coarseBasis.uniformRefine_withTransfer(R_patch, 1, mult);

            for (int k = 0; k < R_patch.outerSize(); ++k)
            {
                for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(R_patch, k); it; ++it)
                {
                    rEntries.add(rRowOffset + it.row(), rColOffset + it.col(), it.value());
                }
            }
            rRowOffset += R_patch.rows();
            rColOffset += R_patch.cols();
        }

        m_R[l - 1].resize(m_nDisjoint[l], m_nDisjoint[l - 1]);
        m_R[l - 1].setFrom(rEntries);
    }
}

template<typename T>
gsMatrix<T> gsAsG1Multigrid<T>::prolongate(index_t l, const gsMatrix<T>& coarseVector) const
{
    gsMatrix<T> v_coarse = m_S[l - 1] * coarseVector;
    gsMatrix<T> v_fine_disjoint = m_R[l - 1] * v_coarse;
    gsMatrix<T> rhs_fine = m_S[l].transpose() * v_fine_disjoint;
    return m_StSSolvers[l]->solve(rhs_fine);
}

template<typename T>
gsMatrix<T> gsAsG1Multigrid<T>::restrict(index_t l, const gsMatrix<T>& fineResidual) const
{
    gsMatrix<T> w_fine = m_StSSolvers[l]->solve(fineResidual);
    gsMatrix<T> z_fine = m_S[l] * w_fine;
    gsMatrix<T> z_coarse = m_R[l - 1].transpose() * z_fine;
    return m_S[l - 1].transpose() * z_coarse;
}

template<typename T>
void gsAsG1Multigrid<T>::smooth(index_t l, gsMatrix<T>& x, const gsMatrix<T>& b, index_t numSteps, T omega) const
{
    const auto& K_l = m_K[l];
    const auto& diagK_l = m_diagK[l];
    const index_t n = m_nFree[l];
    const index_t nIfc = m_ifcDofs[l].size();

    gsMatrix<T> res(n, 1);
    for (index_t step = 0; step < numSteps; ++step)
    {
        res = b - K_l * x;

        // Step 1: Damped Jacobi for interior DOFs
        for (index_t i = 0; i < n; ++i)
        {
            x(i, 0) += omega * res(i, 0) / diagK_l[i];
        }

        // Step 2: Exact solve on Interface Block (Overlapping Schwarz coupling)
        if (nIfc > 0 && m_ifcBlockSolvers[l])
        {
            gsMatrix<T> rIfc(nIfc, 1);
            res = b - K_l * x;
            for (index_t i = 0; i < nIfc; ++i)
            {
                rIfc(i, 0) = res(m_ifcDofs[l][i], 0);
            }

            gsMatrix<T> dxIfc = m_ifcBlockSolvers[l]->solve(rIfc);

            for (index_t i = 0; i < nIfc; ++i)
            {
                x(m_ifcDofs[l][i], 0) += dxIfc(i, 0);
            }
        }
    }
}

template<typename T>
void gsAsG1Multigrid<T>::vCycle(index_t l, gsMatrix<T>& x, const gsMatrix<T>& b, index_t preSmooth, index_t postSmooth) const
{
    if (l == 0)
    {
        x = m_coarseSolvers[0]->solve(b);
        return;
    }

    smooth(l, x, b, preSmooth);

    gsMatrix<T> res = b - m_K[l] * x;
    gsMatrix<T> coarseRes = restrict(l, res);

    gsMatrix<T> coarseCorrection = gsMatrix<T>::Zero(m_nFree[l - 1], 1);
    vCycle(l - 1, coarseCorrection, coarseRes, preSmooth, postSmooth);

    gsMatrix<T> fineCorrection = prolongate(l, coarseCorrection);
    x += fineCorrection;

    smooth(l, x, b, postSmooth);
}

template<typename T>
bool gsAsG1Multigrid<T>::solve(gsMatrix<T>& x, const gsMatrix<T>& b, index_t maxCycles, T tol, index_t preSmooth, index_t postSmooth) const
{
    const index_t L = m_numLevels - 1;
    x = gsMatrix<T>::Zero(m_nFree[L], 1);

    T initialResNorm = (b - m_K[L] * x).norm();
    if (initialResNorm < 1e-15) return true;

    gsInfo << std::setw(8) << "Cycle" << std::setw(18) << "Residual Norm" << std::setw(18) << "Rel Residual" << "\n";
    gsInfo << std::string(44, '-') << "\n";

    for (index_t cycle = 1; cycle <= maxCycles; ++cycle)
    {
        vCycle(L, x, b, preSmooth, postSmooth);
        T resNorm = (b - m_K[L] * x).norm();
        T relRes = resNorm / initialResNorm;

        gsInfo << std::setw(8) << cycle 
               << std::setw(18) << std::scientific << std::setprecision(6) << resNorm 
               << std::setw(18) << std::scientific << std::setprecision(6) << relRes << "\n";

        if (relRes < tol)
        {
            gsInfo << "AS-G1 Multigrid converged in " << cycle << " V-cycles.\n";
            return true;
        }
    }

    gsInfo << "Warning: AS-G1 Multigrid did not reach tolerance in " << maxCycles << " cycles.\n";
    return false;
}

} // namespace gismo
