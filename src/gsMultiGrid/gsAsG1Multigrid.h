/** @file gsAsG1Multigrid.h
 *
 *  @brief Analysis-Suitable G1 (AS-G1) Multi-Patch Multigrid Solver for Biharmonic Equations.
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

#include <vector>
#include <memory>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <gsCore/gsForwardDeclarations.h>
#include <gsSolver/gsPreconditioner.h>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

namespace gismo
{

/** @brief Multigrid Solver for Analysis-Suitable G1 (AS-G1) Multi-Patch Spaces
 *
 *  This class implements an exact space-nesting geometric multigrid solver
 *  specifically designed for AS-G1 discretization of 4th-order biharmonic equations.
 *  Inter-level transfer operators (prolongation and restriction) rely on a sparse
 *  normal equations projection (S_l^T S_l)^{-1} S_l^T R_{l-1->l} S_{l-1}.
 *  Interface smoothing is handled via an Overlapping Block-Schwarz Smoother.
 *
 *  @ingroup Solver
 */
template<typename T>
class gsAsG1Multigrid
{
public:
    typedef memory::shared_ptr<gsAsG1Multigrid> Ptr;
    typedef memory::unique_ptr<gsAsG1Multigrid> uPtr;

    /// Constructor
    gsAsG1Multigrid(const gsMultiPatch<T>& geometry,
                    index_t degree = 3,
                    index_t minRefinements = 2,
                    index_t maxRefinements = 4,
                    index_t numGaussPerSpan = 0);

    /// Set up hierarchy, transformation matrices S_l, refinement matrices R_{l-1->l}, and stiffness matrices K_l
    void setupHierarchy();

    /// Perform inter-level prolongation P_{l-1 -> l} u_{l-1}
    gsMatrix<T> prolongate(index_t l, const gsMatrix<T>& coarseVector) const;

    /// Perform inter-level restriction R_{l -> l-1} r_l = P_{l-1 -> l}^T r_l
    gsMatrix<T> restrict(index_t l, const gsMatrix<T>& fineResidual) const;

    /// Apply Overlapping Block-Schwarz Smoother on level l
    void smooth(index_t l, gsMatrix<T>& x, const gsMatrix<T>& b, index_t numSteps = 2, T omega = 0.5) const;

    /// Execute single Multigrid V-cycle on level l
    void vCycle(index_t l, gsMatrix<T>& x, const gsMatrix<T>& b, index_t preSmooth = 2, index_t postSmooth = 2) const;

    /// Solve system K_L x = b on the finest level using V-cycles
    bool solve(gsMatrix<T>& x, const gsMatrix<T>& b, index_t maxCycles = 50, T tol = 1e-8, index_t preSmooth = 2, index_t postSmooth = 2) const;

    /// Get number of levels
    index_t numLevels() const { return m_numLevels; }

    /// Get stiffness matrix on level l
    const gsSparseMatrix<T>& stiffnessMatrix(index_t l) const { return m_K[l]; }

    /// Get free DOFs count on level l
    index_t freeDofs(index_t l) const { return m_nFree[l]; }

private:
    gsMultiPatch<T> m_baseGeometry;
    index_t m_degree;
    index_t m_minRefinements;
    index_t m_maxRefinements;
    index_t m_numGaussPerSpan;
    index_t m_numLevels;

    std::vector<gsMultiPatch<T>> m_mpLevels;
    std::vector<gsSparseMatrix<T>> m_S;               // Transformation matrices S_l
    std::vector<gsSparseMatrix<T>> m_R;               // B-spline refinement matrices R_{l-1->l}
    std::vector<gsSparseMatrix<T>> m_K;               // Biharmonic stiffness matrices K_l
    std::vector<gsVector<T>> m_diagK;                 // Diagonal of stiffness matrices K_l
    std::vector<std::vector<index_t>> m_ifcDofs;      // Interface DOFs per level
    std::vector<std::shared_ptr<typename gsSparseSolver<T>::SimplicialLDLT>> m_StSSolvers;       // Solvers for (S_l^T S_l)^{-1}
    std::vector<std::shared_ptr<typename gsSparseSolver<T>::SimplicialLDLT>> m_ifcBlockSolvers; // Solvers for K_interface
    std::vector<std::shared_ptr<typename gsSparseSolver<T>::SimplicialLDLT>> m_coarseSolvers;   // Exact solver for coarsest level K_0
    std::vector<index_t> m_nFree;
    std::vector<index_t> m_nDisjoint;
};

} // namespace gismo

#include <gsMultiGrid/gsAsG1Multigrid.hpp>
