/** @file as_g1_biharmonic2_example.cpp

    @brief Direct solver and convergence test for the second biharmonic problem
           (\beta u + \Delta^2 u = f, u = g_1, \Delta u = g_2 on \partial\Omega)
           over Analysis-Suitable G1 (AS-G1) multi-patch geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <algorithm>
#include <cmath>
#include <gismo.h>
#include <iomanip>
#include <iostream>
#include <vector>

#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

using namespace gismo;
using namespace gismo::expr;

int main(int argc, char *argv[])
{
    using T = real_t;

    std::string geometry(
        "domain2d/2patch/two_bilinear_patches.xml");
    std::string outDir("");
    index_t degree = 3;
    index_t numGaussPerSpan = 0;
    index_t minRefinements = 2;
    index_t maxRefinements = 4;
    T freqA = 1;
    T beta = 0;
    T eps = 1e-3;
    bool homogeneousBC = false;
    bool plot = false;

    gsCmdLine cmd("AS-G1 Second Biharmonic Equation Solver and Convergence Experiment.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
    cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
    cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).", numGaussPerSpan);
    cmd.addInt("m", "minRefinements", "Minimum uniform refinement level to test (minimum 2 for AS-G1).", minRefinements);
    cmd.addInt("r", "refinements", "Maximum uniform refinement levels to test.", maxRefinements);
    cmd.addReal("a", "frequency", "Frequency factor 'a' in manufactured solution.", freqA);
    cmd.addReal("b", "beta", "Reaction parameter beta >= 0 in (beta u + Delta^2 u = f).", beta);
    cmd.addReal("e", "eps", "Tolerance for AS-G1 gluing data.", eps);
    cmd.addSwitch("homo", "Use homogeneous boundary conditions (u = 0, Delta u = 0).", homogeneousBC);
    cmd.addSwitch("plot", "Export target function and numerical solution to VTK files.", plot);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    std::string prefix;
    if (!outDir.empty())
    {
        if (outDir.back() != '/')
            outDir += '/';
        gsFileManager::mkdir(outDir);
        prefix = outDir;
        gsInfo << "Output directory: " << outDir << "\n";
    }

    if (degree < 3)
        degree = 3;

    if (minRefinements < 0)
        minRefinements = 0;

    // ---- Read geometry ----
    gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr)
    {
        gsInfo << "Error: Cannot read geometry " << geometry << "\n";
        return -1;
    }
    gsMultiPatch<T> &mp = *mpPtr;
    mp.computeTopology();

    // ---- Degree Elevation ----
    {
        const short_t minDeg = gsMultiBasis<T>(mp).minCwiseDegree();
        const short_t maxDeg = gsMultiBasis<T>(mp).maxCwiseDegree();
        if (minDeg == maxDeg)
            gsInfo << "Geometry has degree " << minDeg << "\n";
        else
            gsInfo << "Geometry has degree between " << minDeg << " and " << maxDeg << "\n";

        if (minDeg < degree)
        {
            const short_t elev = degree - minDeg;
            gsInfo << "Elevate degree by " << elev << " steps.\n";
            mp.degreeElevate(elev);
        }
    }

    // ---- Manufactured Solution ----
    const std::string s_a = std::to_string(freqA);
    const std::string s_b = std::to_string(beta);
    std::string u_expr, grad_x_expr, grad_y_expr, hess_xx, hess_xy, hess_yy, lapl_expr, rhs_expr;
    if (homogeneousBC)
    {
        // Note: exact only on [0,1]^2; tests homogeneous BCs u = 0, Delta u = 0 on general domains
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

    gsInfo << "\n======================================================================\n";
    gsInfo << "AS-G1 Second Biharmonic Solver (" << beta << "*u + Delta^2 u = f, u = g1, Delta u = g2)\n";
    gsInfo << "Domain: " << geometry << " (" << mp.nPatches() << " patches)\n";
    gsInfo << "Target Exact Function: u(x, y) = " << u_expr << "\n";
    gsInfo << "Boundary Mode: " << (homogeneousBC ? "Homogeneous (u=0, Delta u=0)" : "Manufactured (u=g1, Delta u=g2)") << "\n";
    gsInfo << "======================================================================\n\n";

    // ---- Define boundary conditions ("Values" only; normal derivatives stay free) ----
    gsConstantFunction<T> zero;
    gsBoundaryConditions<T> bc;
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
        bc.add(it->patch, it->side(), "Values", zero);

    // ---- Compute Gluing Data ----
    gsMatrix<T> gd = computeGluingData(mp, eps, numGaussPerSpan);

    // ---- Get boundary corners and local frames (Case A: sharp, Case B: smooth) ----
    std::vector<std::vector<patchCorner>> vertices = getBoundaryVertices(mp, bc, "Values");
    std::vector<cornerCondition> cc;
    gsMatrix<T> framesForPatches(mp.nPatches(), 5 * 4);
    framesForPatches.setZero();

    // Default interior frame: d1 = (1, 0), d2 = (0, -1), curvature = 0
    for (size_t i = 0; i < mp.nPatches(); ++i)
        for (index_t c = 0; c < 4; ++c)
        {
            framesForPatches(i, 5 * c + 0) = T(1);
            framesForPatches(i, 5 * c + 3) = T(-1);
        }

    index_t nSmoothVertices = 0, nSharpCorners = 0;
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
    gsInfo << "Boundary vertices: " << vertices.size()
           << " (Case A sharp corners: " << nSharpCorners
           << ", Case B smooth vertices: " << nSmoothVertices << ")\n\n";

    // ---- Table header ----
    gsInfo << std::setw(5) << "r" << std::setw(12) << "h_max" << std::setw(10)
           << "N_free" << std::setw(14) << "L2 Error" << std::setw(10)
           << "L2 Rate" << std::setw(14) << "H1 Error" << std::setw(10)
           << "H1 Rate" << std::setw(14) << "H2 Error" << std::setw(10)
           << "H2 Rate" << "\n";
    gsInfo << std::string(99, '-') << "\n" << std::flush;

    T prev_h = 0;
    T prev_l2 = 0;
    T prev_h1 = 0;
    T prev_h2 = 0;

    for (index_t ref = minRefinements; ref <= maxRefinements; ++ref)
    {
        // ---- Get multi-basis and refine ----
        gsMultiBasis<T> dbasis(mp);
        const index_t mult = 2;
        for (index_t i = 0; i < ref; ++i)
            dbasis.uniformRefine(1, mult);

        T h_max = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
            h_max = std::max(h_max, dbasis[i].getMaxCellLength());
        if (h_max <= 0)
            h_max = std::pow(T(0.5), ref);

        // ---- Build per-patch Argyris AS-G1 embeddings ----
        std::vector<gsArgyrisEmbedding<T>> argBasis;
        gsBlockSparseMatrix<T> argBasisGlobalCreator(mp.nPatches(), mp.nPatches());
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            argBasis.push_back(deriveArgyrisBasisEmbedding(
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[i]),
                gsMatrix<T>(gd.row(i)),
                gsMatrix<T>(framesForPatches.row(i)),
                mp.patch(i)
            ));
            argBasisGlobalCreator.set(i, i, argBasis[i].matrix);
        }
        gsSparseMatrix<T> argBasisGlobal = argBasisGlobalCreator;

        // ---- Set up gsDofMapper ----
        gsDofMapper mapper = makeMapperForArgyrisBasis(mp, argBasis, bc, cc);

        // ---- Build global transformation matrices T_free and T_bnd ----
        gsSparseMatrix<T> T_global = argBasisGlobal * asEmbeddingMatrix<T>(mapper.size(), mapper.asVector()).transpose();

        const index_t nFree = mapper.freeSize(), nBnd = mapper.boundarySize();
        gsSparseMatrix<T> T_free = T_global.leftCols(nFree);
        gsSparseMatrix<T> T_bnd  = T_global.rightCols(nBnd);

        // ---- Evaluate Dirichlet boundary values u = g1 via L^2(dOmega) projection ----
        gsMatrix<T> sol_bnd;
        if (!homogeneousBC && nBnd > 0)
        {
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());

            auto G_map   = A.getMap(mp);
            auto u_space = A.getSpace(dbasis);
            auto u_coeff = A.getCoeff(exact_u, G_map);

            A.initSystem();
            A.assembleBdr(bc.get("Values"),
                          u_space * u_space.tr() * meas(G_map),
                          u_space * u_coeff * meas(G_map));

            gsSparseMatrix<T> M_bnd = T_bnd.transpose() * A.matrix() * T_bnd;
            gsMatrix<T>       F_bnd = T_bnd.transpose() * A.rhs();
            makeSparseCholeskySolver(M_bnd)->apply(F_bnd, sol_bnd);
        }
        else
        {
            sol_bnd.setZero(nBnd, 1);
        }

        // ---- Assemble Discontinuous Second Biharmonic Matrix & RHS ----
        gsSparseMatrix<T> K_discont;
        gsMatrix<T> F_discont;
        {
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());

            auto G_map    = A.getMap(mp);
            auto u_space  = A.getSpace(dbasis);
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

        // ---- Global Linear System Assembly & Solve ----
        gsSparseMatrix<T> K_free = T_free.transpose() * K_discont * T_free;
        gsMatrix<T> F_free;
        if (!homogeneousBC && nBnd > 0)
            F_free = T_free.transpose() * (F_discont - K_discont * (T_bnd * sol_bnd));
        else
            F_free = T_free.transpose() * F_discont;
        gsMatrix<T> sol_free;
        makeSparseCholeskySolver(K_free)->apply(F_free, sol_free);

        // Reconstruct full global vector sol_discont
        gsMatrix<T> sol_discont = T_free * sol_free + T_bnd * sol_bnd;

        // Reconstruct multi-patch solution field
        gsMultiPatch<T> sol;
        index_t offset = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            const index_t sz = argBasis[i].matrix.rows();
            gsMatrix<T> ci = sol_discont.block(offset, 0, sz, 1);
            offset += sz;
            const gsTensorBSplineBasis<2, T> &tb = dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[i]);
            sol.addPatch(tb.makeGeometry(give(ci)));
        }

        // ---- Error Evaluation ----
        gsExprEvaluator<T> ev;
        ev.setIntegrationDomain(dbasis.domain());
        auto G_map_ev = ev.getMap(mp);

        auto u_exact_ev    = ev.getVariable(exact_u, G_map_ev);
        auto grad_exact_ev = ev.getVariable(exact_grad, G_map_ev);
        auto hess_exact_ev = reshape(ev.getVariable(exact_hess, G_map_ev), 2, 2);

        auto u_sol_ev = ev.getVariable(sol);

        const T l2err = std::sqrt(
            ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
        const T h1err = std::sqrt(
            ev.integral((igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() * meas(G_map_ev)));
        const T h2err = std::sqrt(
            ev.integral((ihess(u_sol_ev, G_map_ev) - hess_exact_ev).sqNorm() * meas(G_map_ev)));

        // Rates
        T l2_rate = 0, h1_rate = 0, h2_rate = 0;
        if (ref > minRefinements && prev_h > 0)
        {
            const T h_ratio = prev_h / h_max;
            l2_rate = std::log(prev_l2 / l2err) / std::log(h_ratio);
            h1_rate = std::log(prev_h1 / h1err) / std::log(h_ratio);
            h2_rate = std::log(prev_h2 / h2err) / std::log(h_ratio);
        }

        gsInfo << std::setw(5) << ref << std::setw(12) << std::scientific
               << std::setprecision(4) << h_max << std::setw(10) << nFree
               << std::setw(14) << std::scientific << std::setprecision(4)
               << l2err;
        if (ref > minRefinements)
            gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << l2_rate;
        else
            gsInfo << std::setw(10) << "-";

        gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
        if (ref > minRefinements)
            gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << h1_rate;
        else
            gsInfo << std::setw(10) << "-";

        gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
        if (ref > minRefinements)
            gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << h2_rate;
        else
            gsInfo << std::setw(10) << "-";

        gsInfo << "\n" << std::flush;

        prev_h  = h_max;
        prev_l2 = l2err;
        prev_h1 = h1err;
        prev_h2 = h2err;

        // ---- Optional VTK output ----
        if (plot || !outDir.empty())
        {
            std::string solName   = prefix + "as_g1_biharmonic2_sol_r"   + std::to_string(ref);
            std::string exactName = prefix + "as_g1_biharmonic2_exact_r" + std::to_string(ref);

            gsField<T> solField(mp, sol);
            gsWriteParaview<>(solField, solName, 1000);

            gsField<T> exactField(mp, exact_u, false);
            gsWriteParaview<>(exactField, exactName, 1000);

            gsInfo << " -> Exported VTK plots: " << solName << " & " << exactName << "\n";
        }
    }

    gsInfo << std::string(99, '-') << "\n";
    gsInfo << "Done.\n";
    return 0;
}
