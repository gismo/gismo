/** @file as_g1_biharmonic2_multigrid_example.cpp

    @brief Multigrid solver for the second biharmonic problem
           (\beta u + \Delta^2 u = f, u = g_1, \Delta u = g_2 on \partial\Omega)
           over Analysis-Suitable G1 (AS-G1) multi-patch geometries with
           Symmetric Gauss-Seidel smoother.

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
#include <string>
#include <vector>

#include <gismo.h>
#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>
#include <gsMultiGrid/gsMultiGrid.h>
#include <gsUtils/gsStopwatch.h>

using namespace gismo;
using namespace gismo::expr;

// Struct to hold assembled level data
template <class T>
struct LevelData2
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
};

// Setup geometry and Argyris AS-G1 embedding on a single grid level (without PDE assembly)
template <class T>
LevelData2<T> setupLevelGeometry2(
    const gsMultiPatch<T> &mp,
    const gsMatrix<T> &gd,
    const gsMatrix<T> &framesForPatches,
    const gsBoundaryConditions<T> &bc,
    const std::vector<cornerCondition> &cc,
    index_t ref,
    bool isFinest)
{
    LevelData2<T> ld;
    ld.ref = ref;
    ld.dbasis = gsMultiBasis<T>(mp);

    const index_t mult = 2;
    for (index_t i = 0; i < ref; ++i)
        ld.dbasis.uniformRefine(1, mult);

    gsBlockSparseMatrix<T> argBasisGlobalCreator(mp.nPatches(), mp.nPatches());
    for (size_t i = 0; i < mp.nPatches(); ++i)
    {
        ld.argBasis.push_back(deriveArgyrisBasisEmbedding(
            dynamic_cast<const gsTensorBSplineBasis<2, T> &>(ld.dbasis[i]),
            gsMatrix<T>(gd.row(i)),
            gsMatrix<T>(framesForPatches.row(i)),
            mp.patch(i)));
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

    return ld;
}

// Assemble second biharmonic PDE system and boundary projection ONCE on the finest grid
template <class T>
void assembleFinestSystem2(
    const gsMultiPatch<T> &mp,
    const gsBoundaryConditions<T> &bc,
    LevelData2<T> &ld,
    const gsFunctionExpr<T> &exact_u,
    const gsFunctionExpr<T> &exact_laplace,
    const gsFunctionExpr<T> &rhs_f,
    T beta,
    bool homogeneousBC)
{
    const index_t nBnd = ld.mapper.boundarySize();

    // Essential boundary Dirichlet projection for u = g1
    if (!homogeneousBC && nBnd > 0)
    {
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(ld.dbasis.domain());
        auto G_map   = A.getMap(mp);
        auto u_space = A.getSpace(ld.dbasis);
        auto u_coeff = A.getCoeff(exact_u, G_map);

        A.initSystem();
        A.assembleBdr(bc.get("Values"),
                      u_space * u_space.tr() * meas(G_map),
                      u_space * u_coeff * meas(G_map));

        gsSparseMatrix<T> M_bnd = ld.T_bnd.transpose() * A.matrix() * ld.T_bnd;
        gsMatrix<T>       F_bnd = ld.T_bnd.transpose() * A.rhs();
        makeSparseCholeskySolver(M_bnd)->apply(F_bnd, ld.sol_bnd);
    }
    else
    {
        ld.sol_bnd.setZero(nBnd, 1);
    }

    // Discontinuous Second Biharmonic matrix & RHS (including natural BC Delta u = g2)
    gsSparseMatrix<T> K_discont;
    gsMatrix<T> F_discont;
    {
        gsExprAssembler<T> A(1, 1);
        A.setIntegrationDomain(ld.dbasis.domain());
        auto G_map    = A.getMap(mp);
        auto u_space  = A.getSpace(ld.dbasis);
        auto f_coeff  = A.getCoeff(rhs_f, G_map);
        auto g2_coeff = A.getCoeff(exact_laplace, G_map);

        A.initSystem();
        if (beta > T(0))
            A.assemble((beta * (u_space * u_space.tr()) + ilapl(u_space, G_map) * ilapl(u_space, G_map).tr()) * meas(G_map),
                       u_space * f_coeff * meas(G_map));
        else
            A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                       u_space * f_coeff * meas(G_map));

        if (!homogeneousBC)
            A.assembleBdr(bc.get("Values"), (igrad(u_space, G_map) * unv(G_map)) * g2_coeff.tr() * meas(G_map));

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
gsSparseMatrix<T, RowMajor> computeIntergridProlongation2(
    const LevelData2<T> &coarse,
    const LevelData2<T> &fine,
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
    gsSparseMatrix<T> M_gram  = fine.T_free.transpose() * fine.T_free;
    gsSparseMatrix<T> RHS_mat = fine.T_free.transpose() * (R_disjoint * coarse.T_free);

    auto chol = makeSparseCholeskySolver(M_gram);
    const index_t nFine     = fine.T_free.cols();
    const index_t nCoarse   = coarse.T_free.cols();
    const index_t batchSize = 256;
    const index_t nBatches  = (nCoarse + batchSize - 1) / batchSize;

    std::vector<gsSparseEntries<T>> thread_entries(nBatches);

#pragma omp parallel for schedule(dynamic)
    for (index_t b = 0; b < nBatches; ++b)
    {
        const index_t bStart = b * batchSize;
        const index_t bCols  = std::min(batchSize, nCoarse - bStart);
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
bool runDomainExperiment2(
    const std::string &geomPath,
    index_t degree,
    index_t minRef,
    index_t maxRef,
    index_t numCycles,
    index_t numPreSmooth,
    index_t numPostSmooth,
    const std::string &solverType,
    T tol,
    index_t maxIt,
    bool homogeneousBC,
    T freqA,
    T beta,
    T eps,
    bool verifyDirect,
    bool plot,
    const std::string &outDir,
    const std::string &txtOut)
{
    gsInfo << "\n" << std::string(80, '=') << "\n";
    gsInfo << "Domain: " << geomPath << "\n";
    gsInfo << "PDE: " << beta << "*u + Delta^2 u = f, u = g1, Delta u = g2 on boundary\n";
    gsInfo << "Degree p = " << degree << " | Levels: r = " << minRef << " -> " << maxRef
           << " | Smoother: Symmetric Gauss-Seidel (" << numPreSmooth << " pre, " << numPostSmooth << " post)\n";
    gsInfo << "Solver: " << (solverType == "cg" ? "PCG (Multigrid Preconditioner)" : "Direct Multigrid Iteration")
           << " | Cycles: " << (numCycles == 1 ? "V-cycle" : "W-cycle") << " | Tol = " << tol << "\n";
    gsInfo << "Boundary condition: " << (homogeneousBC ? "Homogeneous (u=0, Delta u=0)" : "Manufactured solution projection (u=g1, Delta u=g2)") << "\n";
    gsInfo << std::string(80, '=') << "\n";

    typename gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geomPath);
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

    // Manufactured / source solution
    const std::string s_a = std::to_string(freqA);
    const std::string s_b = std::to_string(beta);
    std::string u_expr, grad_x_expr, grad_y_expr, hess_xx, hess_xy, hess_yy, lapl_expr, rhs_expr;
    if (homogeneousBC)
    {
        u_expr      = "sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        grad_x_expr =        s_a + "*pi*cos(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        grad_y_expr =        s_a + "*pi*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        hess_xx     = "-" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        hess_xy     =       s_a + "^2*pi^2*cos(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        hess_yy     = "-" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        lapl_expr   = "0";
        rhs_expr    = "(" + s_b + " + 4*(" + s_a + ")^4*pi^4)*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
    }
    else
    {
        u_expr      = "sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        grad_x_expr =        s_a + "*pi*cos(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        grad_y_expr = "-" + s_a + "*pi*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        hess_xx     = "-" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        hess_xy     = "-" + s_a + "^2*pi^2*cos(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        hess_yy     = "-" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        lapl_expr   = "-2*(" + s_a + ")^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        rhs_expr    = "(" + s_b + " + 4*(" + s_a + ")^4*pi^4)*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
    }

    gsFunctionExpr<T> exact_u(u_expr, 2);
    gsFunctionExpr<T> exact_grad(grad_x_expr, grad_y_expr, 2);
    gsFunctionExpr<T> exact_hess(hess_xx, hess_xy, hess_xy, hess_yy, 2);
    gsFunctionExpr<T> exact_laplace(lapl_expr, 2);
    gsFunctionExpr<T> rhs_f(rhs_expr, 2);

    // Boundary conditions ("Values" only)
    gsConstantFunction<T> zero;
    gsBoundaryConditions<T> bc;
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
        bc.add(it->patch, it->side(), "Values", zero);

    // Gluing data & boundary corners (Case A sharp vs Case B smooth)
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

    std::vector<std::vector<patchCorner>> vertices = getBoundaryVertices(mp, bc, "Values");
    std::vector<cornerCondition> cc;
    gsMatrix<T> framesForPatches(mp.nPatches(), 5 * 4);
    framesForPatches.setZero();
    for (size_t i = 0; i < mp.nPatches(); ++i)
        for (index_t c = 0; c < 4; ++c)
        {
            framesForPatches(i, 5 * c + 0) = T(1);
            framesForPatches(i, 5 * c + 3) = T(-1);
        }

    index_t nSmoothVertices = 0, nSharpCorners = 0;
    try
    {
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            gsBoundaryCornerData<T> bcd = getBoundaryCornerData(mp, vertices[i]);
            for (size_t j = 0; j < vertices[i].size(); ++j)
            {
                const index_t p = vertices[i][j].patch;
                const index_t c = vertices[i][j].m_index - 1;
                framesForPatches(p, 5 * c + 0) = bcd.dir1(0);
                framesForPatches(p, 5 * c + 1) = bcd.dir1(1);
                framesForPatches(p, 5 * c + 2) = bcd.dir2(0);
                framesForPatches(p, 5 * c + 3) = bcd.dir2(1);
                framesForPatches(p, 5 * c + 4) = bcd.curvature;
            }
            cc.push_back(cornerCondition{
                vertices[i][0],
                bcd.isSmooth ? cornerConditionType::valuesSmooth : cornerConditionType::valuesSharp
            });
            if (bcd.isSmooth)
                ++nSmoothVertices;
            else
                ++nSharpCorners;
        }
    }
    catch (const std::exception &e)
    {
        gsInfo << "Skipping geometry (corner classification failed: " << e.what() << ")\n";
        return false;
    }
    gsInfo << "Patches: " << mp.nPatches() << " | Boundary vertices: " << vertices.size()
           << " (sharp=" << nSharpCorners << ", smooth=" << nSmoothVertices << ")\n";

    // Set up grid hierarchy from minRef to maxRef
    const index_t nLevels = maxRef - minRef + 1;
    const index_t finest  = nLevels - 1;
    std::vector<LevelData2<T>> levels(nLevels);
    gsInfo << "Setting up grid hierarchy (" << nLevels << " levels)...\n";
    double totalSetupTime = 0;
    for (index_t l = 0; l < nLevels; ++l)
    {
        const index_t ref = minRef + l;
        const bool isFinest = (l == finest);
        gsStopwatch timer;
        timer.restart();
        levels[l] = setupLevelGeometry2(mp, gd, framesForPatches, bc, cc, ref, isFinest);
        if (isFinest)
            assembleFinestSystem2(mp, bc, levels[l], exact_u, exact_laplace, rhs_f, beta, homogeneousBC);
        const double tLevel = timer.stop();
        totalSetupTime += tLevel;
        gsInfo << "  Level " << l << " (ref=" << ref << "): " << levels[l].mapper.freeSize() << " DOFs ("
               << (isFinest ? "Geometry + PDE Assembly: " : "Geometry only: ") << std::fixed << std::setprecision(4) << tLevel << " s)\n";
    }

    // Compute intergrid transfer matrices
    double totalTransferTime = 0;
    std::vector<gsSparseMatrix<T, RowMajor>> transferMatrices(nLevels - 1);
    for (index_t l = 0; l < nLevels - 1; ++l)
    {
        gsStopwatch timer;
        timer.restart();
        transferMatrices[l] = computeIntergridProlongation2(levels[l], levels[l + 1], mp.nPatches());
        const double tTr = timer.stop();
        totalTransferTime += tTr;
        gsInfo << "  Transfer P[" << l << "->" << l + 1 << "]: " << transferMatrices[l].rows()
               << " x " << transferMatrices[l].cols() << " (" << std::fixed << std::setprecision(4) << tTr << " s)\n";
    }

    // Set up gsMultiGridOp with Galerkin coarse-grid operators
    typename gsMultiGridOp<T>::uPtr mg = gsMultiGridOp<T>::make(levels[finest].K_free, transferMatrices);
    mg->setNumCycles(numCycles);
    mg->setNumPreSmooth(numPreSmooth);
    mg->setNumPostSmooth(numPostSmooth);

    // Level 0 uses exact direct solver by default in gsMultiGridOp; levels 1..L use Symmetric Gauss-Seidel smoother.
    // (Setting Gauss-Seidel on level 0 too matches as_g1_biharmonic_multigrid_example.cpp; wait: let's check both!
    // In as_g1_biharmonic_multigrid_example.cpp lines 367-371, setSmoother is called for l = 0..nLevels-1,
    // because gsMultiGridOp::step applies m_coarseSolver at l=0 unless overridden, wait: setSmoother(0, sm) actually sets smoother[0] while level 0 uses base solver! Calling it on 0..nLevels-1 matches 1st biharmonic identically.)
    for (index_t l = 0; l < nLevels; ++l)
    {
        typename gsMultiGridOp<T>::PrecondPtr sm(makeGaussSeidelOp(mg->matrix(l)).release());
        mg->setSmoother(l, sm);
    }

    // Solve system on finest level
    const gsSparseMatrix<T> &A_fine   = levels[finest].K_free;
    const gsMatrix<T>       &rhs_fine = levels[finest].F_free;
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
    else // Direct Multigrid Iteration
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

    gsInfo << "\n--- Multigrid Convergence Results ---\n";
    gsInfo << "  Iterations      : " << numIterations << "\n";
    gsInfo << "  Final Rel. Res. : " << std::scientific << std::setprecision(4) << finalRelRes << "\n";
    gsInfo << "  Solve Time      : " << std::fixed << std::setprecision(4) << solveTime << " s\n";
    gsInfo << "  Setup Time      : " << std::fixed << std::setprecision(4) << totalSetupTime << " s\n";
    gsInfo << "  Transfer Time   : " << std::fixed << std::setprecision(4) << totalTransferTime << " s\n";
    gsInfo << "  Total Time      : " << std::fixed << std::setprecision(4) << (totalSetupTime + totalTransferTime + solveTime) << " s\n";

    // Optional direct solver verification
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

    // Error evaluation against exact solution
    T l2err = 0.0, h1err = 0.0, h2err = 0.0;
    if (!homogeneousBC)
    {
        gsMatrix<T> sol_discont = levels[finest].T_free * sol_free + levels[finest].T_bnd * levels[finest].sol_bnd;
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
        auto G_map_ev      = ev.getMap(mp);
        auto u_exact_ev    = ev.getVariable(exact_u, G_map_ev);
        auto grad_exact_ev = ev.getVariable(exact_grad, G_map_ev);
        auto hess_exact_ev = reshape(ev.getVariable(exact_hess, G_map_ev), 2, 2);
        auto u_sol_ev      = ev.getVariable(sol);

        l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
        h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() * meas(G_map_ev)));
        h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - hess_exact_ev).sqNorm() * meas(G_map_ev)));

        gsInfo << "  L2 Error        : " << std::scientific << std::setprecision(4) << l2err << "\n";
        gsInfo << "  H1 Error        : " << std::scientific << std::setprecision(4) << h1err << "\n";
        gsInfo << "  H2 Error        : " << std::scientific << std::setprecision(4) << h2err << "\n";

        if (plot && !outDir.empty())
        {
            std::string prefix = outDir;
            if (prefix.back() != '/')
                prefix += '/';
            gsFileManager::mkdir(prefix);
            std::string solName = prefix + "as_g1_biharmonic2_mg_sol";
            gsField<T> solField(mp, sol);
            gsWriteParaview<>(solField, solName, 1000);
            gsInfo << "  Exported Paraview plot: " << solName << "\n";
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
            outfile << "geometry\tpatches\tsharpCorners\tsmoothVertices\tdegree\tbeta\tminRef\tmaxRef\tnLevels\tnDOFs\t"
                       "cycles\tpreSmooth\tpostSmooth\tsolver\ttol\titer\trelRes\t"
                       "setupTime\ttransferTime\tsolveTime\ttotalTime\tl2Error\th1Error\th2Error\tstatus\n";

        outfile << geomPath << "\t"
                << mp.nPatches() << "\t"
                << nSharpCorners << "\t"
                << nSmoothVertices << "\t"
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
    std::string suite("classic");
    index_t degree = 3;
    index_t minRef = 2;
    index_t maxRef = 4;
    index_t numCycles = 1;
    index_t numPreSmooth = 1;
    index_t numPostSmooth = 1;
    index_t maxIt = 300;
    T tol = 1e-8;
    T freqA = 1.0;
    T beta = 0.0;
    T eps = 1e-3;
    bool runAll = false;
    bool homogeneousBC = false;
    bool verifyDirect = false;
    bool plot = false;

    gsCmdLine cmd("AS-G1 Second Biharmonic Multigrid Solver with Symmetric Gauss-Seidel Smoother.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for Paraview VTK files.", outDir);
    cmd.addString("t", "txtOut", "Filename for tab-separated text results output.", txtOut);
    cmd.addString("s", "solver", "Iterative solver: 'cg' (PCG) or 'mg' (Direct Multigrid Iteration).", solverType);
    cmd.addString("", "suite", "Benchmark suite when --all is set: classic, generator, noholes, holes, all.", suite);
    cmd.addInt("d", "degree", "Spline degree (minimum 3 for C1 AS-G1).", degree);
    cmd.addInt("m", "minRef", "Coarsest grid refinement level (minimum 2).", minRef);
    cmd.addInt("r", "maxRef", "Finest grid refinement level.", maxRef);
    cmd.addInt("c", "cycles", "Multigrid cycle type: 1 = V-cycle, 2 = W-cycle.", numCycles);
    cmd.addInt("", "pre", "Number of pre-smoothing steps.", numPreSmooth);
    cmd.addInt("", "post", "Number of post-smoothing steps.", numPostSmooth);
    cmd.addInt("", "maxIt", "Maximum solver iterations.", maxIt);
    cmd.addReal("e", "tol", "Relative residual tolerance.", tol);
    cmd.addReal("a", "frequency", "Frequency factor in manufactured solution.", freqA);
    cmd.addReal("b", "beta", "Reaction parameter beta >= 0 in (beta u + Delta^2 u = f).", beta);
    cmd.addReal("", "eps", "Tolerance for AS-G1 gluing data.", eps);
    cmd.addSwitch("homo", "Use homogeneous boundary conditions (u=0, Delta u=0).", homogeneousBC);
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

    if (degree < 3)
        degree = 3;
    if (minRef < 2)
        minRef = 2;
    if (maxRef < minRef)
        maxRef = minRef;

    if (!runAll)
    {
        runDomainExperiment2<T>(
            geometry, degree, minRef, maxRef,
            numCycles, numPreSmooth, numPostSmooth,
            solverType, tol, maxIt, homogeneousBC, freqA, beta, eps,
            verifyDirect, plot, outDir, txtOut);
    }
    else
    {
        std::vector<std::string> benchmarkFiles;

        auto addClassic = [&]() {
            const std::vector<std::string> list = {
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
            };
            benchmarkFiles.insert(benchmarkFiles.end(), list.begin(), list.end());
        };

        auto addGenerator = [&]() {
            const std::vector<std::string> list = {
                "domain2d/as_g1_generator/2_patch_rectangle_asg1.xml",
                "domain2d/as_g1_generator/two_squares_asg1.xml",
                "domain2d/as_g1_generator/ldomain_asg1.xml",
                "domain2d/as_g1_generator/lake_asg1.xml",
                "domain2d/as_g1_generator/as_g1_tri_p3.xml",
                "domain2d/as_g1_generator/as_g1_pentagon_p3.xml",
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
                "domain2d/as_g1_generator/as_g1_pentagon_random_s101.xml",
                "domain2d/as_g1_generator/as_g1_hex_random_s777.xml",
                "domain2d/as_g1_generator/as_g1_star_v5_p3.xml",
                "domain2d/as_g1_generator/as_g1_ellipse_v5_p3.xml",
                "domain2d/as_g1_generator/yeti_asg1.xml",
                "domain2d/lshape_9patch_non_bilinear_trig_p3_asg1.xml",
                "domain2d/ring_7patch_non_bilinear_swirl_p3_asg1.xml"
            };
            benchmarkFiles.insert(benchmarkFiles.end(), list.begin(), list.end());
        };

        auto addNoHoles = [&]() {
            const std::vector<std::string> list = {
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
            };
            benchmarkFiles.insert(benchmarkFiles.end(), list.begin(), list.end());
        };

        auto addHoles = [&]() {
            const std::vector<std::string> list = {
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/bean_0corner_1hole_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_ring_4patch_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/wavy_ring_annulus_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_lshape_6corner_1hole_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/star_5corner_1hole_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/crescent_2corner_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_2patch_asg1.xml",
                "domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/curved_lshape_3patch_asg1.xml"
            };
            benchmarkFiles.insert(benchmarkFiles.end(), list.begin(), list.end());
        };

        if (suite == "classic")
            addClassic();
        else if (suite == "generator")
            addGenerator();
        else if (suite == "noholes")
            addNoHoles();
        else if (suite == "holes")
            addHoles();
        else // "all"
        {
            addClassic();
            addGenerator();
            addNoHoles();
            addHoles();
        }

        gsInfo << "\n======================================================================\n";
        gsInfo << "RUNNING SECOND BIHARMONIC MULTIGRID SUITE (" << suite << ")\n";
        gsInfo << "======================================================================\n";

        index_t numPassed = 0;
        index_t numFailed = 0;

        for (const auto &file : benchmarkFiles)
        {
            if (!gsFileManager::fileExists(file) && !gsFileManager::fileExists("filedata/" + file))
            {
                gsInfo << "Skipping non-existent: " << file << "\n";
                continue;
            }
            bool ok = runDomainExperiment2<T>(
                file, degree, minRef, maxRef,
                numCycles, numPreSmooth, numPostSmooth,
                solverType, tol, maxIt, homogeneousBC, freqA, beta, eps,
                verifyDirect, false, "", txtOut);
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
