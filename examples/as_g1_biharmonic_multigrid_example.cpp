/** @file as_g1_biharmonic_multigrid_example.cpp

    @brief Multigrid solver for the first biharmonic problem
           (\beta u + \Delta^2 u = f, u = g_1, \partial_n u = g_2 on \partial\Omega)
           over Analysis-Suitable G1 (AS-G1) multi-patch geometries with
           Symmetric Gauss-Seidel smoother.

    Supports inhomogeneous manufactured solutions as well as homogeneous (zero)
    boundary conditions (u = 0, \partial_n u = 0) using a compactly supported
    tensor-product B-spline bubble, a smooth trigonometric bubble, or a
    box-conforming trigonometric solution, with exact L2, H1, and H2 error
    evaluation and optional PVD/VTK export.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gismo.h>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>
#include <gsModeling/gsAsG1Smoother.hpp>
#include <gsMultiGrid/gsMultiGrid.h>
#include <gsUtils/gsStopwatch.h>

using namespace gismo;
using namespace gismo::expr;

static void mkdirRecursive(const std::string &path)
{
    if (path.empty())
        return;
    std::string cur;
    for (size_t i = 0; i < path.size(); ++i)
    {
        cur.push_back(path[i]);
        if (path[i] == '/' && cur.size() > 1)
            gsFileManager::mkdir(cur);
    }
    if (cur.back() != '/')
        gsFileManager::mkdir(cur);
}

static std::string sanitizeStem(const std::string &path)
{
    std::string base = gsFileManager::getBasename(path);
    for (char &c : base)
        if (c == ' ' || c == '/')
            c = '_';
    return base;
}

// Struct to hold assembled level data
template <class T>
struct LevelData
{
    index_t ref;
    gsMultiBasis<T> dbasis;
    std::vector<gsArgyrisEmbedding<T>> argBasis;
    gsDofMapper mapper;
    gsSparseMatrix<T> T_free;
    gsSparseMatrix<T> T_bnd;
    gsSparseMatrix<T> K_free;
    gsMatrix<T> F_free;
    gsMatrix<T> sol_bnd;
    std::vector<index_t> skeletonDofs;
    std::vector<std::vector<index_t>> smootherBlocks;
};

// Setup geometry and Argyris C1 embedding on a single grid level (without PDE assembly)
template <class T>
LevelData<T> setupLevelGeometry(
    const gsMultiPatch<T> &mp,
    const gsMatrix<T> &gd,
    const gsMatrix<T> &normalsForPatches,
    const gsBoundaryConditions<T> &bc,
    const std::vector<cornerCondition> &cc,
    index_t ref,
    bool isFinest,
    const std::string &smootherMode = "")
{
    LevelData<T> ld;
    ld.ref = ref;
    ld.dbasis = gsMultiBasis<T>(mp);

    const index_t mult = 2;
    for (index_t i = 0; i < ref; ++i)
        ld.dbasis.uniformRefine(1, mult);

    // Uniform parametric cell scale for this level (shared across all patches so interface/vertex matching stays exact)
    T hLevel = ld.dbasis[0].getMinCellLength();
    for (size_t i = 1; i < mp.nPatches(); ++i)
        hLevel = std::min(hLevel, ld.dbasis[i].getMinCellLength());
    if (hLevel <= T(0))
        hLevel = std::pow(T(0.5), ref);

    gsBlockSparseMatrix<T> argBasisGlobalCreator(mp.nPatches(), mp.nPatches());
    for (size_t i = 0; i < mp.nPatches(); ++i)
    {
        ld.argBasis.push_back(deriveArgyrisBasisEmbedding(
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(ld.dbasis[i]),
            gsMatrix<T>(gd.row(i)),
            gsMatrix<T>(normalsForPatches.row(i)),
            mp.patch(i)));

        // Rescale derivative DOFs uniformly (1st deriv by 1/h, 2nd deriv by 1/h^2) to balance T^T T and K_free
        const auto &sz = ld.argBasis[i].sizes;
        gsMatrix<T> scaleVec = gsMatrix<T>::Ones(ld.argBasis[i].matrix.cols(), 1);
        index_t off = sz[0];
        for (int s = 0; s < 4; ++s)
        {
            off += sz[1 + 2 * s];
            for (index_t k = 0; k < sz[2 + 2 * s]; ++k)
                scaleVec(off + k, 0) = T(1) / hLevel;
            off += sz[2 + 2 * s];
        }
        for (int c = 0; c < 4; ++c)
        {
            scaleVec(off + 1, 0) = T(1) / hLevel;
            scaleVec(off + 2, 0) = T(1) / hLevel;
            scaleVec(off + 3, 0) = T(1) / (hLevel * hLevel);
            scaleVec(off + 4, 0) = T(1) / (hLevel * hLevel);
            scaleVec(off + 5, 0) = T(1) / (hLevel * hLevel);
            off += 6;
        }
        ld.argBasis[i].matrix = ld.argBasis[i].matrix * scaleVec.asDiagonal();
        argBasisGlobalCreator.set(i, i, ld.argBasis[i].matrix);
    }
    gsSparseMatrix<T> argBasisGlobal = argBasisGlobalCreator;

    ld.mapper = makeMapperForArgyrisBasis(mp, ld.argBasis, bc, cc);
    gsSparseMatrix<T> T_global = argBasisGlobal * asEmbeddingMatrix<T>(ld.mapper.size(), ld.mapper.asVector()).transpose();

    const index_t nFree = ld.mapper.freeSize();
    ld.T_free = T_global.leftCols(nFree);

    if (isFinest)
    {
        const index_t nBnd = ld.mapper.boundarySize();
        ld.T_bnd = T_global.rightCols(nBnd);
    }

    if (smootherMode == "ring1")
        ld.skeletonDofs = collectCoupledRing1Dofs(mp, ld.dbasis, ld.argBasis, ld.mapper);
    else if (smootherMode == "line_ring1")
        ld.smootherBlocks = collectInteriorLineAndRing1Blocks(mp, ld.dbasis, ld.argBasis, ld.mapper);

    return ld;
}

// Assemble first biharmonic PDE system and boundary projection ONCE on the finest grid
template <class T>
void assembleFinestSystem(
    const gsMultiPatch<T> &mp,
    const gsBoundaryConditions<T> &bc,
    LevelData<T> &ld,
    const gsFunctionExpr<T> *exact_u,
    const gsFunctionExpr<T> *exact_grad,
    const gsFunctionSet<T> *rhs_f,
    const gsBubbleData<T> *bubbleData,
    T beta,
    bool homogeneousBC)
{
    const index_t nBnd = ld.mapper.boundarySize();
    const bool useBubble = (bubbleData != nullptr);

    // Boundary Dirichlet values projection
    if (!homogeneousBC && nBnd > 0)
    {
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(ld.dbasis.domain());
        auto G_map = A.getMap(mp);
        auto u_space = A.getSpace(ld.dbasis);
        auto u_coeff = A.getCoeff(*exact_u, G_map);
        auto u_grad_coeff = A.getCoeff(*exact_grad, G_map);

        A.initSystem();
        A.assembleBdr(bc.get("ValuesAndDerivatives"),
                      u_space * u_space.tr() * meas(G_map) + (igrad(u_space, G_map) * unv(G_map)) * (igrad(u_space, G_map) * unv(G_map)).tr() * meas(G_map),
                      u_space * u_coeff * meas(G_map) + (igrad(u_space, G_map) * unv(G_map)) * (u_grad_coeff.tr() * unv(G_map)) * meas(G_map));

        gsSparseMatrix<T> M_bnd = ld.T_bnd.transpose() * A.matrix() * ld.T_bnd;
        gsMatrix<T> F_bnd = ld.T_bnd.transpose() * A.rhs();
        makeSparseCholeskySolver(M_bnd)->apply(F_bnd, ld.sol_bnd);
    }
    else
    {
        ld.sol_bnd.setZero(nBnd, 1);
    }

    // Discontinuous Biharmonic matrix & RHS
    gsSparseMatrix<T> K_discont;
    gsMatrix<T> F_discont;
    {
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(ld.dbasis.domain());
        if (useBubble)
            A.options().setReal("quA", 2.0);
        auto G_map = A.getMap(mp);
        auto u_space = A.getSpace(ld.dbasis);

        A.initSystem();
        if (useBubble)
        {
            auto u_ex_coeff   = A.getCoeff(bubbleData->field, G_map);
            auto lap_ex_coeff = A.getCoeff(bubbleData->lapField, G_map);
            if (beta > T(0))
                A.assemble((beta * (u_space * u_space.tr()) + ilapl(u_space, G_map) * ilapl(u_space, G_map).tr()) * meas(G_map),
                           (beta * (u_space * u_ex_coeff) + ilapl(u_space, G_map) * lap_ex_coeff) * meas(G_map));
            else
                A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                           ilapl(u_space, G_map) * lap_ex_coeff * meas(G_map));
        }
        else
        {
            auto f_coeff = A.getCoeff(*rhs_f, G_map);
            if (beta > T(0))
                A.assemble((beta * (u_space * u_space.tr()) + ilapl(u_space, G_map) * ilapl(u_space, G_map).tr()) * meas(G_map),
                           u_space * f_coeff * meas(G_map));
            else
                A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                           u_space * f_coeff * meas(G_map));
        }

        K_discont = give(A.matrix());
        F_discont = give(A.rhs());
    }

    ld.K_free = ld.T_free.transpose() * K_discont * ld.T_free;
    if (!homogeneousBC && nBnd > 0)
        ld.F_free = ld.T_free.transpose() * (F_discont - K_discont * (ld.T_bnd * ld.sol_bnd));
    else
        ld.F_free = ld.T_free.transpose() * F_discont;
}

// Compute intergrid prolongation matrix P_{coarse -> fine}
template <class T>
gsSparseMatrix<T, RowMajor> computeIntergridProlongation(
    const LevelData<T> &coarse,
    const LevelData<T> &fine,
    size_t nPatches)
{
    gsBlockSparseMatrix<T> R_block(nPatches, nPatches);
    for (size_t p = 0; p < nPatches; ++p)
    {
        gsSparseMatrix<T, RowMajor> localTr;
        const gsTensorBSplineBasis<2, T> &cBasis =
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(coarse.dbasis[p]);
        gsTensorBSplineBasis<2, T> cBasisCopy(cBasis);
        cBasisCopy.uniformRefine_withTransfer(localTr, 1, 2);
        gsSparseMatrix<T> localTrCol = localTr;
        R_block.set(p, p, localTrCol);
    }
    gsSparseMatrix<T> R_disjoint = R_block;

    // P_{c->f} = (T_free_fine^T * T_free_fine)^{-1} * (T_free_fine^T * R_disjoint * T_free_coarse)
    gsSparseMatrix<T> M_gram = fine.T_free.transpose() * fine.T_free;
    gsSparseMatrix<T> RHS_mat = fine.T_free.transpose() * (R_disjoint * coarse.T_free);

    auto chol = makeSparseCholeskySolver(M_gram);
    const index_t nFine = fine.T_free.cols();
    const index_t nCoarse = coarse.T_free.cols();
    const index_t batchSize = 256;
    const index_t nBatches = (nCoarse + batchSize - 1) / batchSize;

    std::vector<gsSparseEntries<T>> thread_entries(nBatches);

#pragma omp parallel for schedule(dynamic)
    for (index_t b = 0; b < nBatches; ++b)
    {
        const index_t bStart = b * batchSize;
        const index_t bCols = std::min(batchSize, nCoarse - bStart);
        gsMatrix<T> RHS_batch = RHS_mat.middleCols(bStart, bCols).toDense();
        gsMatrix<T> P_batch;
        chol->apply(RHS_batch, P_batch);
        for (index_t c = 0; c < bCols; ++c)
        {
            const index_t j = bStart + c;
            for (index_t i = 0; i < nFine; ++i)
            {
                const T val = P_batch(i, c);
                if (std::abs(val) > 1e-12)
                    thread_entries[b].add(i, j, val);
            }
        }
    }

    gsSparseEntries<T> entries;
    size_t total_entries = 0;
    for (index_t b = 0; b < nBatches; ++b)
        total_entries += thread_entries[b].size();
    entries.reserve(total_entries);
    for (index_t b = 0; b < nBatches; ++b)
        entries.insert(entries.end(), thread_entries[b].begin(), thread_entries[b].end());

    gsSparseMatrix<T> P_sparse(nFine, nCoarse);
    P_sparse.setFromTriplets(entries.begin(), entries.end());
    gsSparseMatrix<T, RowMajor> P_rowMajor = P_sparse;
    return P_rowMajor;
}

// Run multigrid experiment on a given geometry
template <class T>
bool runDomainExperiment(
    const std::string &geomPath,
    index_t degree,
    index_t bubbleDeg,
    index_t minRef,
    index_t maxRef,
    index_t numCycles,
    index_t numPreSmooth,
    index_t numPostSmooth,
    const std::string &solverType,
    const std::string &smootherType,
    index_t scmsRp,
    T tol,
    index_t maxIt,
    bool homogeneousBC,
    const std::string &exactType,
    T freqA,
    T beta,
    T forceVal,
    T eps,
    T margin,
    bool verifyDirect,
    bool plot,
    const std::string &outDir,
    const std::string &plotSubDir,
    const std::string &txtOut)
{
    const bool useBubble     = homogeneousBC && (exactType == "bubble" || exactType == "sin_bubble");
    const bool useConstForce = homogeneousBC && (exactType == "const");

    const std::string smootherLabel =
        (smootherType == "scms") ? "SCMS Multiplicative + Overlapping Skeleton (Remark 5)" :
        (smootherType == "scms_add") ? "SCMS Additive Schwarz + Overlapping Skeleton (Theorem 5, isotropic)" :
        (smootherType == "scms_add_metric") ? "SCMS Additive Schwarz + Overlapping Skeleton (metric-weighted)" :
        (smootherType == "scms_add_analytic") ? "SCMS Additive Schwarz (Theorem 5, analytic 1/(8*c_bar_G))" :
        (smootherType == "scms_nonov") ? "SCMS Multiplicative + Non-Overlapping Skeleton (r_p=2)" :
        (smootherType == "fd_ring1") ? "Point GS + Patch Interior FD + Coupled Ring-1" :
        (smootherType == "line_ring1") ? "Alternating Line Strips + Coupled Ring-1" :
        (smootherType == "ring1") ? "Point GS + Coupled Ring-1 Skeleton" :
        (smootherType == "sgs" || (smootherType == "auto" && solverType != "cg")) ? "Symmetric Gauss-Seidel" :
        "Point Gauss-Seidel (fwd pre, bwd post)";

    gsInfo << "\n" << std::string(80, '=') << "\n";
    gsInfo << "Domain: " << geomPath << "\n";
    gsInfo << "PDE: " << beta << "*u + Delta^2 u = f, first biharmonic (u, d_n u on boundary)\n";
    gsInfo << "Degree p = " << degree << " | Levels: r = " << minRef << " -> " << maxRef
           << " | Smoother: " << smootherLabel << " (" << numPreSmooth << " pre, " << numPostSmooth << " post)\n";
    gsInfo << "Solver: " << (solverType == "cg" ? "PCG (Multigrid Preconditioner)" : "Direct Multigrid Iteration")
           << " | Cycles: " << (numCycles == 1 ? "V-cycle" : "W-cycle") << " | Tol = " << tol << "\n";
    gsInfo << "Boundary condition: " << (homogeneousBC ? "Homogeneous zero (u=0, d_n u=0)" : "Manufactured solution projection")
           << " | Exact: " << (homogeneousBC ? exactType : "sin*cos") << "\n";
    gsInfo << std::string(80, '=') << "\n";

    typename gsMultiPatch<T>::uPtr mpPtr;
    try
    {
        mpPtr = gsReadFile<>(geomPath);
    }
    catch (const std::exception &e)
    {
        gsInfo << "Skipping geometry (read failed: " << e.what() << ")\n";
        return false;
    }
    if (!mpPtr)
    {
        gsInfo << "Error: Cannot load geometry " << geomPath << "\n";
        return false;
    }
    gsMultiPatch<T> &mp = *mpPtr;
    mp.computeTopology();

    // Degree elevation if necessary
    const short_t minDeg = gsMultiBasis<T>(mp).minCwiseDegree();
    if (minDeg < degree)
    {
        const short_t elev = degree - minDeg;
        mp.degreeElevate(elev);
    }

    // Gluing data & boundary corners
    gsMatrix<T> gd;
    try
    {
        gd = computeGluingData(mp, eps, 0);
    }
    catch (const std::exception &e)
    {
        gsInfo << "Skipping geometry (gluing data failed: " << e.what() << ")\n";
        return false;
    }

    gsConstantFunction<T> zero;
    gsBoundaryConditions<T> bc;
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
        bc.add(it->patch, it->side(), "ValuesAndDerivatives", zero);

    std::vector<std::vector<patchCorner>> vertices = getBoundaryVertices(mp, bc, "ValuesAndDerivatives");
    std::vector<cornerCondition> cc;
    gsMatrix<T> normalsForPatches(mp.nPatches(), 2 * 4);
    normalsForPatches.setZero();
    try
    {
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            gsVector<T> normal = getOuterNormalDerivative(mp, vertices[i]);
            for (size_t j = 0; j < vertices[i].size(); ++j)
                normalsForPatches.block(vertices[i][j].patch, 2 * (vertices[i][j].m_index - 1), 1, 2) = normal.transpose();
            cc.push_back(cornerCondition{vertices[i][0], normal.norm() < T(1e-6) ? cornerConditionType::all : cornerConditionType::valuesNormals});
        }
    }
    catch (const std::exception &e)
    {
        gsInfo << "Skipping geometry (corner classification failed: " << e.what() << ")\n";
        return false;
    }
    for (size_t i = 0; i < mp.nPatches(); ++i)
        for (index_t j = 0; j < 4; ++j)
            if (normalsForPatches.block(i, 2 * j, 1, 2).norm() < T(1e-6))
            {
                normalsForPatches(i, 2 * j) = 1;
                normalsForPatches(i, 2 * j + 1) = 0;
            }

    // Manufactured / source solution
    gsBubbleData<T> bubbleData;
    std::unique_ptr<gsFunctionExpr<T>> exact_u_ptr, exact_grad_ptr, exact_hess_ptr, rhs_f_ptr;
    std::unique_ptr<gsConstantFunction<T>> const_rhs_ptr;

    if (useBubble)
    {
        const index_t bDeg = (bubbleDeg > 0) ? bubbleDeg : (degree + 2);
        if (exactType == "sin_bubble")
            bubbleData = makeTrigBubble(mp, std::max<index_t>(bDeg, 4), margin);
        else
            bubbleData = makeBSplineBubble(mp, bDeg, margin);

        gsInfo << "  Bubble box      : [" << bubbleData.a_x << ", " << bubbleData.b_x
               << "] x [" << bubbleData.a_y << ", " << bubbleData.b_y << "]\n";
    }
    else if (useConstForce)
    {
        const_rhs_ptr.reset(new gsConstantFunction<T>(forceVal, 2));
        gsInfo << "  Constant force  : f = " << forceVal << " (no closed-form exact solution)\n";
    }
    else if (homogeneousBC)
    {
        gsBoxExactExprs ex = makeBoxHomogeneousExpressions(mp, false, freqA, beta);
        exact_u_ptr.reset(new gsFunctionExpr<T>(ex.u_expr, 2));
        exact_grad_ptr.reset(new gsFunctionExpr<T>(ex.grad_x_expr, ex.grad_y_expr, 2));
        exact_hess_ptr.reset(new gsFunctionExpr<T>(ex.hess_xx, ex.hess_xy, ex.hess_xy, ex.hess_yy, 2));
        rhs_f_ptr.reset(new gsFunctionExpr<T>(ex.rhs_expr, 2));
    }
    else
    {
        const std::string s_a = std::to_string(freqA);
        const std::string s_b = std::to_string(beta);
        const std::string u_expr      = "sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string grad_x_expr =        s_a + "*pi*cos(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string grad_y_expr = "-" +  s_a + "*pi*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        const std::string rhs_expr    = "(" + s_b + " + 4*(" + s_a + ")^4*pi^4)*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string hess_xx     = "-" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string hess_xy     = "-" + s_a + "^2*pi^2*cos(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        const std::string hess_yy     = "-" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";

        exact_u_ptr.reset(new gsFunctionExpr<T>(u_expr, 2));
        exact_grad_ptr.reset(new gsFunctionExpr<T>(grad_x_expr, grad_y_expr, 2));
        exact_hess_ptr.reset(new gsFunctionExpr<T>(hess_xx, hess_xy, hess_xy, hess_yy, 2));
        rhs_f_ptr.reset(new gsFunctionExpr<T>(rhs_expr, 2));
    }

    // Set up grid hierarchy from minRef to maxRef
    const std::string levelSmootherMode = (smootherType == "ring1" || smootherType == "line_ring1") ? smootherType : std::string();
    const index_t nLevels = maxRef - minRef + 1;
    const index_t finest = nLevels - 1;
    std::vector<LevelData<T>> levels(nLevels);
    gsInfo << "Setting up grid hierarchy (" << nLevels << " levels)...\n";
    double totalSetupTime = 0;
    try
    {
        for (index_t l = 0; l < nLevels; ++l)
        {
            const index_t ref = minRef + l;
            const bool isFinest = (l == finest);
            gsStopwatch timer;
            timer.restart();
            levels[l] = setupLevelGeometry(mp, gd, normalsForPatches, bc, cc, ref, isFinest, l > 0 ? levelSmootherMode : std::string());
            if (isFinest)
            {
                const gsFunctionSet<T> *rhsPtr = useConstForce
                    ? static_cast<const gsFunctionSet<T> *>(const_rhs_ptr.get())
                    : static_cast<const gsFunctionSet<T> *>(rhs_f_ptr.get());
                assembleFinestSystem(
                    mp, bc, levels[l],
                    exact_u_ptr.get(), exact_grad_ptr.get(), rhsPtr,
                    useBubble ? &bubbleData : nullptr,
                    beta, homogeneousBC);
            }
            const double tLevel = timer.stop();
            totalSetupTime += tLevel;
            gsInfo << "  Level " << l << " (ref=" << ref << "): " << levels[l].mapper.freeSize() << " DOFs ("
                   << (isFinest ? "Geometry + PDE Assembly: " : "Geometry only: ") << std::fixed << std::setprecision(4) << tLevel << " s)\n";
        }
    }
    catch (const std::exception &e)
    {
        gsInfo << "Skipping geometry (setup/assembly failed: " << e.what() << ")\n";
        return false;
    }

    // Compute intergrid transfer matrices
    double totalTransferTime = 0;
    std::vector<gsSparseMatrix<T, RowMajor>> transferMatrices(nLevels - 1);
    for (index_t l = 0; l < nLevels - 1; ++l)
    {
        gsStopwatch timer;
        timer.restart();
        transferMatrices[l] = computeIntergridProlongation(levels[l], levels[l + 1], mp.nPatches());
        const double tTr = timer.stop();
        totalTransferTime += tTr;
        gsInfo << "  Transfer P[" << l << "->" << l + 1 << "]: " << transferMatrices[l].rows()
               << " x " << transferMatrices[l].cols() << " (" << std::fixed << std::setprecision(4) << tTr << " s)\n";
    }

    // Set up gsMultiGridOp
    typename gsMultiGridOp<T>::uPtr mg = gsMultiGridOp<T>::make(levels[finest].K_free, transferMatrices);
    mg->setNumCycles(numCycles);
    mg->setNumPreSmooth(numPreSmooth);
    mg->setNumPostSmooth(numPostSmooth);
    mg->setCoarseSolver(makeSparseCholeskySolver(mg->matrix(0)));

    gsStopwatch smTimer;
    smTimer.restart();
    for (index_t l = 1; l < nLevels; ++l)
    {
        typename gsMultiGridOp<T>::PrecondPtr sm;
        if (smootherType == "scms")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1ScmsSmootherOp<T>::make(
                    mg->matrix(l), mp, levels[l].dbasis, levels[l].argBasis, levels[l].mapper,
                    beta, scmsRp, /*additive=*/false, /*useMetricWeights=*/true, /*calibrateByPowerMethod=*/true).release());
        else if (smootherType == "scms_add")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1ScmsSmootherOp<T>::make(
                    mg->matrix(l), mp, levels[l].dbasis, levels[l].argBasis, levels[l].mapper,
                    beta, scmsRp, /*additive=*/true, /*useMetricWeights=*/false, /*calibrateByPowerMethod=*/true).release());
        else if (smootherType == "scms_add_metric")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1ScmsSmootherOp<T>::make(
                    mg->matrix(l), mp, levels[l].dbasis, levels[l].argBasis, levels[l].mapper,
                    beta, scmsRp, /*additive=*/true, /*useMetricWeights=*/true, /*calibrateByPowerMethod=*/true).release());
        else if (smootherType == "scms_add_analytic")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1ScmsSmootherOp<T>::make(
                    mg->matrix(l), mp, levels[l].dbasis, levels[l].argBasis, levels[l].mapper,
                    beta, scmsRp, /*additive=*/true, /*useMetricWeights=*/false, /*calibrateByPowerMethod=*/false).release());
        else if (smootherType == "scms_nonov")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1ScmsSmootherOp<T>::make(
                    mg->matrix(l), mp, levels[l].dbasis, levels[l].argBasis, levels[l].mapper,
                    beta, /*r_p=*/2, /*additive=*/false, /*useMetricWeights=*/true, /*calibrateByPowerMethod=*/true).release());
        else if (smootherType == "fd_ring1")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1FastDiagSmootherOp<T>::make(
                    mg->matrix(l), mp, levels[l].dbasis, levels[l].argBasis, levels[l].mapper, beta).release());
        else if (smootherType == "line_ring1")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1BlockSmootherOp<T>::make(mg->matrix(l), levels[l].smootherBlocks).release());
        else if (smootherType == "ring1")
            sm = typename gsMultiGridOp<T>::PrecondPtr(
                gsAsG1SkeletonSmootherOp<T>::make(mg->matrix(l), levels[l].skeletonDofs).release());
        else if (smootherType == "sgs" || (smootherType == "auto" && solverType != "cg"))
            sm = typename gsMultiGridOp<T>::PrecondPtr(makeSymmetricGaussSeidelOp(mg->matrix(l)).release());
        else
            sm = typename gsMultiGridOp<T>::PrecondPtr(makeGaussSeidelOp(mg->matrix(l)).release());
        mg->setSmoother(l, sm);
    }
    const double smootherSetupTime = smTimer.stop();

    // Solve system on finest level
    const gsSparseMatrix<T> &A_fine = levels[finest].K_free;
    const gsMatrix<T> &rhs_fine = levels[finest].F_free;
    gsMatrix<T> sol_free;
    sol_free.setZero(A_fine.rows(), 1);

    gsStopwatch solveTimer;
    solveTimer.restart();

    index_t numIterations = 0;
    T finalRelRes = 1.0;
    gsMatrix<T> errHist;

    if (solverType == "cg")
    {
        gsOptionList opt;
        opt.addReal("Tolerance", "Tolerance", tol);
        opt.addInt("MaxIterations", "MaxIterations", maxIt);
        gsConjugateGradient<T> cg(A_fine, typename gsMultiGridOp<T>::Ptr(mg.release()));
        cg.setOptions(opt);
        cg.solveDetailed(rhs_fine, sol_free, errHist);
        numIterations = errHist.rows() - 1;
        finalRelRes = errHist(numIterations, 0);
    }
    else
    {
        gsMatrix<T> r = rhs_fine - A_fine * sol_free;
        const T r0_norm = r.norm();
        finalRelRes = 1.0;
        errHist.resize(maxIt + 1, 1);
        errHist(0, 0) = 1.0;

        for (index_t it = 1; it <= maxIt; ++it)
        {
            mg->step(rhs_fine, sol_free);
            r = rhs_fine - A_fine * sol_free;
            finalRelRes = r.norm() / (r0_norm > 1e-15 ? r0_norm : 1.0);
            errHist(it, 0) = finalRelRes;
            numIterations = it;
            if (finalRelRes < tol)
                break;
        }
    }
    const double solveTime = solveTimer.stop();

    const T avgRate = (numIterations > 0 && finalRelRes > T(0)) ? std::pow(finalRelRes, T(1) / T(numIterations)) : T(0);
    gsInfo << "\n--- Multigrid Convergence Results ---\n";
    gsInfo << "  Iterations      : " << numIterations << "\n";
    gsInfo << "  Avg. Conv. Rate : " << std::fixed << std::setprecision(4) << avgRate << "\n";
    gsInfo << "  Final Rel. Res. : " << std::scientific << std::setprecision(4) << finalRelRes << "\n";
    gsInfo << "  Solve Time      : " << std::fixed << std::setprecision(4) << solveTime << " s\n";
    gsInfo << "  Smoother Setup  : " << std::fixed << std::setprecision(4) << smootherSetupTime << " s\n";
    gsInfo << "  Setup Time      : " << std::fixed << std::setprecision(4) << totalSetupTime << " s\n";
    gsInfo << "  Transfer Time   : " << std::fixed << std::setprecision(4) << totalTransferTime << " s\n";
    gsInfo << "  Total Time      : " << std::fixed << std::setprecision(4) << (totalSetupTime + totalTransferTime + smootherSetupTime + solveTime) << " s\n";

    if (verifyDirect)
    {
        gsMatrix<T> sol_direct;
        gsStopwatch dirTimer;
        dirTimer.restart();
        makeSparseCholeskySolver(A_fine)->apply(rhs_fine, sol_direct);
        const double dirTime = dirTimer.stop();
        const T diffNorm = (sol_free - sol_direct).norm() / std::max(T(1e-15), sol_direct.norm());
        gsInfo << "  Direct Cholesky : " << std::fixed << std::setprecision(4) << dirTime
               << " s | Rel. Diff vs MG = " << std::scientific << std::setprecision(4) << diffNorm << "\n";
    }

    // Reconstruct full solution and compute exact errors (both inhomogeneous and homogeneous BCs)
    gsMatrix<T> sol_discont = levels[finest].T_free * sol_free;
    if (!homogeneousBC && levels[finest].mapper.boundarySize() > 0)
        sol_discont += levels[finest].T_bnd * levels[finest].sol_bnd;

    gsMultiPatch<T> sol;
    index_t offset = 0;
    for (size_t i = 0; i < mp.nPatches(); ++i)
    {
        const index_t sz = levels[finest].argBasis[i].matrix.rows();
        gsMatrix<T> ci = sol_discont.block(offset, 0, sz, 1);
        offset += sz;
        const gsTensorBSplineBasis<2, T> &tb = dynamic_cast<const gsTensorBSplineBasis<2, T> &>(levels[finest].dbasis[i]);
        sol.addPatch(tb.makeGeometry(give(ci)));
    }

    gsExprEvaluator<T> ev;
    ev.setIntegrationDomain(levels[finest].dbasis.domain());
    if (useBubble)
        ev.options().setReal("quA", 2.0);
    auto G_map_ev = ev.getMap(mp);
    auto u_sol_ev = ev.getVariable(sol);

    T l2err = 0.0, h1err = 0.0, h2err = 0.0;
    if (useBubble)
    {
        auto u_exact_ev = ev.getVariable(bubbleData.field, G_map_ev);
        l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
        h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev)).sqNorm() * meas(G_map_ev)));
        h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev)).sqNorm() * meas(G_map_ev)));
        gsInfo << "  L2 Error        : " << std::scientific << std::setprecision(4) << l2err << "\n";
        gsInfo << "  H1 Error        : " << std::scientific << std::setprecision(4) << h1err << "\n";
        gsInfo << "  H2 Error        : " << std::scientific << std::setprecision(4) << h2err << "\n";
    }
    else if (useConstForce)
    {
        const T l2norm = std::sqrt(ev.integral(u_sol_ev.sqNorm() * meas(G_map_ev)));
        const T h1norm = std::sqrt(ev.integral(igrad(u_sol_ev, G_map_ev).sqNorm() * meas(G_map_ev)));
        const T h2norm = std::sqrt(ev.integral(ihess(u_sol_ev, G_map_ev).sqNorm() * meas(G_map_ev)));
        gsInfo << "  ||u||_L2        : " << std::scientific << std::setprecision(4) << l2norm << "\n";
        gsInfo << "  |u|_H1          : " << std::scientific << std::setprecision(4) << h1norm << "\n";
        gsInfo << "  |u|_H2          : " << std::scientific << std::setprecision(4) << h2norm << "\n";
    }
    else
    {
        auto u_exact_ev    = ev.getVariable(*exact_u_ptr, G_map_ev);
        auto grad_exact_ev = ev.getVariable(*exact_grad_ptr, G_map_ev);
        auto hess_exact_ev = reshape(ev.getVariable(*exact_hess_ptr, G_map_ev), 2, 2);

        l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
        h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() * meas(G_map_ev)));
        h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - hess_exact_ev).sqNorm() * meas(G_map_ev)));
        gsInfo << "  L2 Error        : " << std::scientific << std::setprecision(4) << l2err << "\n";
        gsInfo << "  H1 Error        : " << std::scientific << std::setprecision(4) << h1err << "\n";
        gsInfo << "  H2 Error        : " << std::scientific << std::setprecision(4) << h2err << "\n";
    }

    if (plot && !outDir.empty())
    {
        std::string prefix = outDir;
        if (prefix.back() != '/')
            prefix += '/';
        if (!plotSubDir.empty())
        {
            prefix += plotSubDir;
            if (prefix.back() != '/')
                prefix += '/';
        }
        mkdirRecursive(prefix);

        const std::string solName   = prefix + "sol";
        const std::string exactName = prefix + "exact";

        gsField<T> solField(mp, sol);
        gsWriteParaview<>(solField, solName, 1000);

        if (useBubble)
        {
            gsField<T> exactField(mp, bubbleData.field, false);
            gsWriteParaview<>(exactField, exactName, 1000);
            gsInfo << "  Exported PVD    : " << solName << ".pvd & " << exactName << ".pvd\n";
        }
        else if (!useConstForce)
        {
            gsField<T> exactField(mp, *exact_u_ptr, false);
            gsWriteParaview<>(exactField, exactName, 1000);
            gsInfo << "  Exported PVD    : " << solName << ".pvd & " << exactName << ".pvd\n";
        }
        else
        {
            gsInfo << "  Exported PVD    : " << solName << ".pvd\n";
        }
    }

    const bool converged = (finalRelRes <= tol);
    gsInfo << "  Status          : " << (converged ? "CONVERGED" : "NOT CONVERGED") << "\n";

    if (!txtOut.empty())
    {
        const bool exists = gsFileManager::fileExists(txtOut);
        std::ofstream outfile;
        outfile.open(txtOut.c_str(), std::ios_base::app);

        if (!exists)
            outfile << "geometry\t"
                       "patches\t"
                       "degree\t"
                       "beta\t"
                       "minRef\t"
                       "maxRef\t"
                       "nLevels\t"
                       "nDOFs\t"
                       "cycles\t"
                       "preSmooth\t"
                       "postSmooth\t"
                       "solver\t"
                       "exact\t"
                       "tol\t"
                       "iter\t"
                       "relRes\t"
                       "setupTime\t"
                       "transferTime\t"
                       "solveTime\t"
                       "totalTime\t"
                       "l2Error\t"
                       "h1Error\t"
                       "h2Error\t"
                       "status\n";

        outfile << geomPath << "\t"
                << mp.nPatches() << "\t"
                << degree << "\t"
                << beta << "\t"
                << minRef << "\t"
                << maxRef << "\t"
                << nLevels << "\t"
                << levels[finest].mapper.freeSize() << "\t"
                << numCycles << "\t"
                << numPreSmooth << "\t"
                << numPostSmooth << "\t"
                << solverType << "\t"
                << (homogeneousBC ? exactType : "inhomo") << "\t"
                << tol << "\t"
                << numIterations << "\t"
                << finalRelRes << "\t"
                << totalSetupTime << "\t"
                << totalTransferTime << "\t"
                << solveTime << "\t"
                << (totalSetupTime + totalTransferTime + solveTime) << "\t"
                << l2err << "\t"
                << h1err << "\t"
                << h2err << "\t"
                << (converged ? "CONVERGED" : "FAILED") << "\n";
    }

    return converged;
}

int main(int argc, char *argv[])
{
    using T = real_t;

    std::string geometry("domain2d/2patch/two_bilinear_patches.xml");
    std::string outDir("");
    std::string txtOut("");
    std::string solverType("cg");
    std::string smootherType("auto");
    std::string suite("classic");
    std::string exactType("bubble");
    index_t degree = 3;
    index_t bubbleDeg = 0;
    index_t minRef = 2;
    index_t maxRef = 4;
    index_t numCycles = 1;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t maxIt = 300;
    index_t scmsRp = -1;
    T tol = 1e-8;
    T freqA = 1.0;
    T beta = 0.0;
    T forceVal = 1.0;
    T eps = 1e-3;
    T margin = 0.15;
    bool runAll = false;
    bool homogeneousBC = false;
    bool constForce = false;
    bool verifyDirect = false;
    bool plot = false;

    gsCmdLine cmd("AS-G1 Biharmonic Multigrid Solver with Symmetric Gauss-Seidel Smoother.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for Paraview VTK/PVD files.", outDir);
    cmd.addString("t", "txtOut", "Filename for tab-separated text results output.", txtOut);
    cmd.addString("s", "solver", "Iterative solver: 'cg' (PCG) or 'mg' (Direct Multigrid Iteration).", solverType);
    cmd.addString("", "smoother", "Smoother: 'auto', 'gs', 'sgs', 'ring1', 'line_ring1', 'fd_ring1', 'scms', 'scms_add', 'scms_add_metric', 'scms_add_analytic', or 'scms_nonov'.", smootherType);
    cmd.addInt("", "scmsRp", "SCMS skeleton strip width r_p (-1 = p+1, 2 = non-overlapping, >=3 = custom).", scmsRp);
    cmd.addString("", "suite", "Benchmark suite when --all is set: classic, generator, noholes, holes, all.", suite);
    cmd.addString("", "exact", "Homogeneous exact solution (--homo): 'bubble' (B-spline), 'sin_bubble', 'box', or 'const'.", exactType);
    cmd.addInt("d", "degree", "Spline degree (minimum 3 for C1 AS-G1).", degree);
    cmd.addInt("", "bubbleDeg", "Degree/power of bubble exact solution (0 = degree + 2).", bubbleDeg);
    cmd.addInt("m", "minRef", "Coarsest grid refinement level (minimum 2).", minRef);
    cmd.addInt("r", "maxRef", "Finest grid refinement level.", maxRef);
    cmd.addInt("c", "cycles", "Multigrid cycle type: 1 = V-cycle, 2 = W-cycle.", numCycles);
    cmd.addInt("", "pre", "Number of pre-smoothing steps.", numPreSmooth);
    cmd.addInt("", "post", "Number of post-smoothing steps.", numPostSmooth);
    cmd.addInt("", "maxIt", "Maximum solver iterations.", maxIt);
    cmd.addReal("e", "tol", "Relative residual tolerance.", tol);
    cmd.addReal("a", "frequency", "Frequency factor in manufactured solution.", freqA);
    cmd.addReal("b", "beta", "Reaction parameter beta >= 0 in (beta u + Delta^2 u = f).", beta);
    cmd.addReal("", "force", "Constant RHS force value when --constForce or --exact const is used.", forceVal);
    cmd.addReal("", "eps", "Tolerance for AS-G1 gluing data.", eps);
    cmd.addReal("", "margin", "Fractional margin for bubble support inside inscribed rectangle (0..0.5).", margin);
    cmd.addSwitch("homo", "Use homogeneous zero Dirichlet boundary conditions.", homogeneousBC);
    cmd.addSwitch("constForce", "Use homogeneous zero BCs with constant force f = forceVal (no exact solution).", constForce);
    cmd.addSwitch("verify", "Verify iterative solution against direct Cholesky solver.", verifyDirect);
    cmd.addSwitch("all", "Run test suite across multi-patch benchmark domains.", runAll);
    cmd.addSwitch("plot", "Export Paraview visualization.", plot);

    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int rv)
    {
        return rv;
    }

    if (constForce)
    {
        homogeneousBC = true;
        exactType = "const";
    }

    if (degree < 3)
        degree = 3;
    if (minRef < 2)
        minRef = 2;
    if (maxRef < minRef)
        maxRef = minRef;
    if (margin <= T(0) || margin >= T(0.5))
        margin = T(0.15);

    if (!runAll)
    {
        runDomainExperiment<T>(
            geometry, degree, bubbleDeg, minRef, maxRef,
            numCycles, numPreSmooth, numPostSmooth,
            solverType, smootherType, scmsRp, tol, maxIt, homogeneousBC, exactType, freqA, beta, forceVal, eps, margin,
            verifyDirect, plot, outDir, sanitizeStem(geometry), txtOut);
    }
    else
    {
        struct SuiteEntry { std::string tag; std::string path; };
        std::vector<SuiteEntry> benchmarkFiles;

        auto addList = [&](const std::string &tag, const std::vector<std::string> &list) {
            for (const auto &f : list)
                benchmarkFiles.push_back({tag, f});
        };

        auto addClassic = [&]() {
            addList("classic", {
                "domain2d/2patch/two_bilinear_patches.xml",
                "domain2d/2patch/2_patch_rectangle.xml",
                "domain2d/2patch/2_patch_rectangle_non_bilinear.xml",
                "domain2d/2patch/three_patch_belt.xml",
                "domain2d/2patch/three_patch_belt_non_bilinear.xml",
                "domain2d/2patch/4_patch_3_valence_not_bilinear.xml",
                "domain2d/2patch/5_patch_5_valence_not_bilinear.xml",
                "domain2d/2patch/6_patch_4_valence_each_not_bilinear.xml",
                "domain2d/2patch/6_patch_6_valence_not_bilinear.xml",
                "domain2d/2patch/square_6_patch.xml",
                "domain2d/2patch/square_6_patch_non_bilinear.xml",
                "domain2d/2patch/two_patch_mathematica.xml",
                "domain2d/2patch/two_patch_mathematica_non_bilinear.xml",
                "domain2d/2patch/weirdo_multivalence.xml",
                "domain2d/2patch/weirdo_multivalence_non_bilinear.xml",
                "domain2d/2patch/experiment_multivalence_not_bilinear.xml",
                "domain2d/2patch/convex_3_valence.xml",
                "domain2d/2patch/convex_3_valence_non_bilinear.xml",
                "domain2d/2patch/convex_5_valence.xml",
                "domain2d/2patch/convex_5_valence_non_bilinear.xml",
                "domain2d/2patch/convex_6_valence.xml",
                "domain2d/2patch/convex_6_valence_non_bilinear.xml",
                "domain2d/2patch/blob_geometry.xml",
                "domain2d/2patch/blob_geometry_non_bilinear.xml",
                "domain2d/2patch/circle/asg1_circle_5patch.xml",
                "domain2d/2patch/circle/asg1_general_5patch_v2.xml",
                "domain2d/2patch/circle/asg1_general_5patch_v4.xml"
            });
        };

        auto addGenerator = [&]() {
            addList("generator", {
                "domain2d/as_g1_generator/2_patch_rectangle_asg1.xml",
                "domain2d/as_g1_generator/two_squares_asg1.xml",
                "domain2d/as_g1_generator/ldomain_asg1.xml",
                "domain2d/as_g1_generator/lake_asg1.xml",
                "domain2d/as_g1_generator/as_g1_tri_p3.xml",
                "domain2d/as_g1_generator/as_g1_pentagon_p3.xml",
                "domain2d/as_g1_generator/as_g1_pentagon_p3_r1.xml",
                "domain2d/as_g1_generator/as_g1_hex_p3.xml",
                "domain2d/as_g1_generator/as_g1_v3_p3.xml",
                "domain2d/as_g1_generator/as_g1_v4_p3.xml",
                "domain2d/as_g1_generator/as_g1_v5_p3_r1.xml",
                "domain2d/as_g1_generator/as_g1_v6_p3.xml",
                "domain2d/as_g1_generator/as_g1_v7_p3.xml",
                "domain2d/as_g1_generator/as_g1_patch3_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch4_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch5_polygon_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch5_star_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch5_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch6_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch7_sym.xml",
                "domain2d/as_g1_generator/as_g1_patch3_random.xml",
                "domain2d/as_g1_generator/as_g1_patch4_random.xml",
                "domain2d/as_g1_generator/as_g1_patch5_random.xml",
                "domain2d/as_g1_generator/as_g1_patch6_random.xml",
                "domain2d/as_g1_generator/as_g1_patch7_random.xml",
                "domain2d/as_g1_generator/as_g1_pentagon_random_s101.xml",
                "domain2d/as_g1_generator/as_g1_hex_random_s777.xml",
                "domain2d/as_g1_generator/as_g1_star_v5_p3.xml",
                "domain2d/as_g1_generator/as_g1_ellipse_v5_p3.xml",
                "domain2d/as_g1_generator/yeti_asg1.xml",
                "domain2d/lshape_9patch_non_bilinear_trig_p3_asg1.xml",
                "domain2d/ring_7patch_non_bilinear_swirl_p3_asg1.xml"
            });
        };

        auto addNoHoles = [&]() {
            addList("noholes", {
                "domain2d/as_g1_no_holes_reparametrized/benchmark_6/quad_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/benchmark_6/curved_quad_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/benchmark_6/lshape_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/benchmark_6/curved_lshape_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/benchmark_6/circle_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/comb_no_holes/comb_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/comb_no_holes/curved_comb_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/bulged_quad_non_bilinear_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/teardrop_1corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/lens_2corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/crescent_2corner_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/curved_triangle_3corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/concave_triangle_3corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/star_triangle_3corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/curved_quad_4corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/concave_quad_4corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/star_quad_4corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/curved_pentagon_5corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/concave_pentagon_5corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/curved_hexagon_6corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/star_hexagon_6corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/curved_heptagon_7corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/star_heptagon_7corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/heart_1corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/1to7corners/heart_2corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/01_dart_4corner_nonconvex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/02_serpentine_s_duct_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/03_filleted_u_bracket_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/04_heptagon_7corner_convex_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/05_hourglass_pinch_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/06_y_junction_3way_asg1.xml",
                "domain2d/as_g1_no_holes_reparametrized/noHolesCounterExamples/07_asymmetric_staircase_8corner_asg1.xml"
            });
        };

        auto addHoles = [&]() {
            addList("holes", {
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/bean_0corner_1hole_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_ring_4patch_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/wavy_ring_annulus_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_lshape_6corner_1hole_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/star_5corner_1hole_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/crescent_2corner_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_2patch_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_lshape_3patch_asg1.xml"
            });
        };

        if (suite == "classic")
            addClassic();
        else if (suite == "generator")
            addGenerator();
        else if (suite == "noholes")
            addNoHoles();
        else if (suite == "holes")
            addHoles();
        else
        {
            addClassic();
            addGenerator();
            addNoHoles();
            addHoles();
        }

        gsInfo << "\n======================================================================\n";
        gsInfo << "RUNNING FIRST BIHARMONIC MULTIGRID SUITE (" << suite << ")\n";
        gsInfo << "======================================================================\n";

        index_t numPassed = 0;
        index_t numFailed = 0;

        for (const auto &entry : benchmarkFiles)
        {
            const std::string &file = entry.path;
            if (!gsFileManager::fileExists(file) && !gsFileManager::fileExists("filedata/" + file))
            {
                gsInfo << "Skipping non-existent: " << file << "\n";
                continue;
            }
            const std::string subDir = entry.tag + "/" + sanitizeStem(file);
            bool ok = runDomainExperiment<T>(
                file, degree, bubbleDeg, minRef, maxRef,
                numCycles, numPreSmooth, numPostSmooth,
                solverType, smootherType, scmsRp, tol, maxIt, homogeneousBC, exactType, freqA, beta, forceVal, eps, margin,
                verifyDirect, plot, outDir, subDir, txtOut);
            if (ok)
                numPassed++;
            else
                numFailed++;
        }

        gsInfo << "\n======================================================================\n";
        gsInfo << "BENCHMARK SUMMARY: " << numPassed << " PASSED, " << numFailed << " FAILED\n";
        gsInfo << "======================================================================\n";
    }

    return 0;
}
