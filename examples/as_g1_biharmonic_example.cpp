/** @file as_g1_biharmonic_example.cpp

    @brief Direct solver and convergence test for the first (and optionally second)
           biharmonic problem over Analysis-Suitable G1 (AS-G1) multi-patch geometries.

    Supports inhomogeneous manufactured solutions as well as homogeneous (zero)
    boundary conditions using a compactly supported tensor-product B-spline bubble,
    a smooth trigonometric bubble, or a box-conforming trigonometric solution,
    with exact L2, H1, and H2 error evaluation.

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
#include <memory>
#include <string>
#include <vector>

#include <gsModeling/gsAsG1Basis.hpp>
#include <gsModeling/gsAsG1Domain.hpp>

using namespace gismo;
using namespace gismo::expr;

int main(int argc, char *argv[])
{
    using T = real_t;

    std::string geometry("domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
    std::string outDir("");
    std::string exactType("bubble");
    index_t degree = 3;
    index_t bubbleDeg = 0;
    index_t numGaussPerSpan = 0;
    index_t minRefinements = 2;
    index_t maxRefinements = 4;
    T freqA = 1;
    T beta = 0;
    T forceVal = 1;
    T eps = 1e-3;
    T margin = 0.15;
    bool homogeneousBC = false;
    bool constForce = false;
    bool secondBiharmonic = false;
    bool plot = false;

    gsCmdLine cmd("AS-G1 Biharmonic Equation Solver and Convergence Experiment.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
    cmd.addString("", "exact", "Homogeneous exact solution (--homo): 'bubble' (B-spline), 'sin_bubble', 'box', or 'const'.", exactType);
    cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
    cmd.addInt("", "bubbleDeg", "Degree/power of bubble exact solution (0 = degree + 2).", bubbleDeg);
    cmd.addInt("n", "numGauss", "Gauss points per knot span (0 = auto).", numGaussPerSpan);
    cmd.addInt("m", "minRefinements", "Minimum uniform refinement level to test.", minRefinements);
    cmd.addInt("r", "refinements", "Maximum uniform refinement levels to test.", maxRefinements);
    cmd.addReal("a", "frequency", "Frequency factor 'a' in manufactured trigonometric solution.", freqA);
    cmd.addReal("b", "beta", "Reaction parameter beta >= 0 in (beta u + Delta^2 u = f).", beta);
    cmd.addReal("", "force", "Constant RHS force value when --constForce or --exact const is used.", forceVal);
    cmd.addReal("e", "eps", "Tolerance for AS-G1 gluing data.", eps);
    cmd.addReal("", "margin", "Fractional margin for bubble support inside inscribed rectangle (0..0.5).", margin);
    cmd.addSwitch("homo", "Use homogeneous zero boundary conditions.", homogeneousBC);
    cmd.addSwitch("constForce", "Use homogeneous zero BCs with constant force f = forceVal (no exact solution).", constForce);
    cmd.addSwitch("second", "Solve second biharmonic problem (u=g1, Delta u=g2) instead of first (u=g1, d_n u=g2).", secondBiharmonic);
    cmd.addSwitch("plot", "Export target exact function and numerical solution to VTK/PVD files.", plot);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    if (constForce)
    {
        homogeneousBC = true;
        exactType = "const";
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
    if (maxRefinements < minRefinements)
        maxRefinements = minRefinements;
    if (margin <= T(0) || margin >= T(0.5))
        margin = T(0.15);

    // ---- Read geometry ----
    typename gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
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

    // ---- Manufactured Exact Solution / Constant Force ----
    const bool useBubble     = homogeneousBC && (exactType == "bubble" || exactType == "sin_bubble");
    const bool useConstForce = homogeneousBC && (exactType == "const");
    gsBubbleData<T> bubbleData;
    std::unique_ptr<gsFunctionExpr<T>> exact_u_ptr, exact_grad_ptr, exact_hess_ptr, exact_lap_ptr, rhs_f_ptr;
    std::unique_ptr<gsConstantFunction<T>> const_rhs_ptr;

    if (useBubble)
    {
        const index_t bDeg = (bubbleDeg > 0) ? bubbleDeg : (degree + 2);
        if (exactType == "sin_bubble")
            bubbleData = makeTrigBubble(mp, std::max<index_t>(bDeg, 4), margin);
        else
            bubbleData = makeBSplineBubble(mp, bDeg, margin);
    }
    else if (useConstForce)
    {
        const_rhs_ptr.reset(new gsConstantFunction<T>(forceVal, 2));
    }
    else if (homogeneousBC)
    {
        gsBoxExactExprs ex = makeBoxHomogeneousExpressions(mp, secondBiharmonic, freqA, beta);
        exact_u_ptr.reset(new gsFunctionExpr<T>(ex.u_expr, 2));
        exact_grad_ptr.reset(new gsFunctionExpr<T>(ex.grad_x_expr, ex.grad_y_expr, 2));
        exact_hess_ptr.reset(new gsFunctionExpr<T>(ex.hess_xx, ex.hess_xy, ex.hess_xy, ex.hess_yy, 2));
        exact_lap_ptr.reset(new gsFunctionExpr<T>(ex.lapl_expr, 2));
        rhs_f_ptr.reset(new gsFunctionExpr<T>(ex.rhs_expr, 2));
    }
    else
    {
        const std::string s_a = std::to_string(freqA);
        const std::string s_b = std::to_string(beta);
        const std::string u_expr      = "sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string grad_x_expr =        s_a + "*pi*cos(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string grad_y_expr = "-" +  s_a + "*pi*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        const std::string hess_xx     = "-" +  s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string hess_xy     = "-" +  s_a + "^2*pi^2*cos(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
        const std::string hess_yy     = "-" +  s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string lapl_expr   = "-2*(" + s_a + ")^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
        const std::string rhs_expr    = "(" + s_b + " + 4*(" + s_a + ")^4*pi^4)*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";

        exact_u_ptr.reset(new gsFunctionExpr<T>(u_expr, 2));
        exact_grad_ptr.reset(new gsFunctionExpr<T>(grad_x_expr, grad_y_expr, 2));
        exact_hess_ptr.reset(new gsFunctionExpr<T>(hess_xx, hess_xy, hess_xy, hess_yy, 2));
        exact_lap_ptr.reset(new gsFunctionExpr<T>(lapl_expr, 2));
        rhs_f_ptr.reset(new gsFunctionExpr<T>(rhs_expr, 2));
    }

    gsInfo << "\n======================================================================\n";
    gsInfo << "AS-G1 " << (secondBiharmonic ? "Second" : "First")
           << " Biharmonic Solver (" << beta << "*u + Delta^2 u = f)\n";
    gsInfo << "Domain: " << geometry << " (" << mp.nPatches() << " patches)\n";
    gsInfo << "Boundary Mode: "
           << (homogeneousBC
                   ? (secondBiharmonic ? "Homogeneous (u = 0, Delta u = 0)"
                                       : "Homogeneous (u = 0, d_n u = 0)")
                   : (secondBiharmonic ? "Manufactured (u = g1, Delta u = g2)"
                                       : "Manufactured (u = g1, d_n u = g2)"))
           << "\n";
    if (useBubble)
    {
        const index_t bDeg = (bubbleDeg > 0) ? bubbleDeg : (degree + 2);
        gsInfo << "Exact Solution: " << (exactType == "sin_bubble" ? "Trig bubble (q=" : "B-spline bubble (deg=")
               << bDeg << ") on [" << bubbleData.a_x << ", " << bubbleData.b_x << "] x ["
               << bubbleData.a_y << ", " << bubbleData.b_y << "]\n";
    }
    else if (useConstForce)
        gsInfo << "Right-Hand Side: Constant force f = " << forceVal << " (no closed-form exact solution)\n";
    else if (homogeneousBC)
        gsInfo << "Exact Solution: Box trigonometric (" << (secondBiharmonic ? "sin*sin" : "sin^2*sin^2") << ")\n";
    else
        gsInfo << "Exact Solution: sin(" << freqA << "*pi*x)*cos(" << freqA << "*pi*y)\n";
    gsInfo << "======================================================================\n\n";

    // ---- Define boundary conditions ----
    gsConstantFunction<T> zero;
    gsBoundaryConditions<T> bc;
    const std::string bcLabel = secondBiharmonic ? "Values" : "ValuesAndDerivatives";
    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
        bc.add(it->patch, it->side(), bcLabel, zero);

    // ---- Compute Gluing Data ----
    gsMatrix<T> gd = computeGluingData(mp, eps, numGaussPerSpan);

    // ---- Get boundary corners and associated normals / frames ----
    std::vector<std::vector<patchCorner>> vertices = getBoundaryVertices(mp, bc, bcLabel);
    std::vector<cornerCondition> cc;
    gsMatrix<T> cornerDataForPatches;

    if (secondBiharmonic)
    {
        cornerDataForPatches.setZero(mp.nPatches(), 5 * 4);
        for (size_t i = 0; i < mp.nPatches(); ++i)
            for (index_t c = 0; c < 4; ++c)
            {
                cornerDataForPatches(i, 5 * c + 0) = T(1);
                cornerDataForPatches(i, 5 * c + 3) = T(-1);
            }
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            gsBoundaryCornerData<T> bcd = getBoundaryCornerData(mp, vertices[i]);
            for (size_t j = 0; j < vertices[i].size(); ++j)
            {
                const index_t p = vertices[i][j].patch;
                const index_t c = vertices[i][j].m_index - 1;
                cornerDataForPatches(p, 5 * c + 0) = bcd.dir1(0);
                cornerDataForPatches(p, 5 * c + 1) = bcd.dir1(1);
                cornerDataForPatches(p, 5 * c + 2) = bcd.dir2(0);
                cornerDataForPatches(p, 5 * c + 3) = bcd.dir2(1);
                cornerDataForPatches(p, 5 * c + 4) = bcd.curvature;
            }
            cc.push_back(cornerCondition{
                vertices[i][0],
                bcd.isSmooth ? cornerConditionType::valuesSmooth : cornerConditionType::valuesSharp
            });
        }
    }
    else
    {
        cornerDataForPatches.setZero(mp.nPatches(), 2 * 4);
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            gsVector<T> normal = getOuterNormalDerivative(mp, vertices[i]);
            for (size_t j = 0; j < vertices[i].size(); ++j)
                cornerDataForPatches.block(vertices[i][j].patch, 2 * (vertices[i][j].m_index - 1), 1, 2) = normal.transpose();
            cc.push_back(cornerCondition{
                vertices[i][0],
                normal.norm() < T(1e-6) ? cornerConditionType::all : cornerConditionType::valuesNormals
            });
        }
        for (size_t i = 0; i < mp.nPatches(); ++i)
            for (index_t j = 0; j < 4; ++j)
                if (cornerDataForPatches.block(i, 2 * j, 1, 2).norm() < T(1e-6))
                {
                    cornerDataForPatches(i, 2 * j) = T(1);
                    cornerDataForPatches(i, 2 * j + 1) = T(0);
                }
    }

    // ---- Table header ----
    if (useConstForce)
    {
        gsInfo << std::setw(5) << "r" << std::setw(12) << "h_max" << std::setw(10)
               << "N_free" << std::setw(14) << "||u||_L2" << std::setw(10)
               << "-" << std::setw(14) << "|u|_H1" << std::setw(10)
               << "-" << std::setw(14) << "|u|_H2" << std::setw(10)
               << "-" << "\n";
    }
    else
    {
        gsInfo << std::setw(5) << "r" << std::setw(12) << "h_max" << std::setw(10)
               << "N_free" << std::setw(14) << "L2 Error" << std::setw(10)
               << "L2 Rate" << std::setw(14) << "H1 Error" << std::setw(10)
               << "H1 Rate" << std::setw(14) << "H2 Error" << std::setw(10)
               << "H2 Rate"
               << "\n";
    }
    gsInfo << std::string(99, '-') << "\n" << std::flush;

    T prev_h = 0;
    T prev_l2 = 0;
    T prev_h1 = 0;
    T prev_h2 = 0;

    for (index_t ref = minRefinements; ref <= maxRefinements; ++ref)
    {
        // ---- Get multi basis & refine ----
        gsMultiBasis<T> dbasis(mp);
        const index_t mult = 2;
        for (index_t i = 0; i < ref; ++i)
            dbasis.uniformRefine(1, mult);

        T h_max = 0;
        for (size_t i = 0; i < mp.nPatches(); ++i)
            h_max = std::max(h_max, dbasis[i].getMaxCellLength());
        if (h_max <= T(0))
            h_max = std::pow(T(0.5), ref);

        // ---- Build per-patch interface embeddings ----
        std::vector<gsArgyrisEmbedding<T>> argBasis;
        gsBlockSparseMatrix<T> argBasisGlobalCreator(mp.nPatches(), mp.nPatches());
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            argBasis.push_back(deriveArgyrisBasisEmbedding(
                dynamic_cast<const gsTensorBSplineBasis<2, T> &>(dbasis[i]),
                gsMatrix<T>(gd.row(i)),
                gsMatrix<T>(cornerDataForPatches.row(i)),
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

        // ---- Evaluate Dirichlet boundary values vector g_bnd ----
        gsMatrix<T> sol_bnd;
        if (!homogeneousBC && nBnd > 0)
        {
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());

            auto G_map   = A.getMap(mp);
            auto u_space = A.getSpace(dbasis);
            auto u_coeff = A.getCoeff(*exact_u_ptr, G_map);
            auto u_grad_coeff = A.getCoeff(*exact_grad_ptr, G_map);

            A.initSystem();
            if (secondBiharmonic)
            {
                A.assembleBdr(bc.get("Values"),
                              u_space * u_space.tr() * meas(G_map),
                              u_space * u_coeff * meas(G_map));
            }
            else
            {
                A.assembleBdr(bc.get("ValuesAndDerivatives"),
                              u_space * u_space.tr() * meas(G_map) + (igrad(u_space, G_map) * unv(G_map)) * (igrad(u_space, G_map) * unv(G_map)).tr() * meas(G_map),
                              u_space * u_coeff * meas(G_map) + (igrad(u_space, G_map) * unv(G_map)) * (u_grad_coeff.tr() * unv(G_map)) * meas(G_map));
            }

            gsSparseMatrix<T> M_bnd = T_bnd.transpose() * A.matrix() * T_bnd;
            gsMatrix<T>       F_bnd = T_bnd.transpose() * A.rhs();
            makeSparseCholeskySolver(M_bnd)->apply(F_bnd, sol_bnd);
        }
        else
        {
            sol_bnd.setZero(nBnd, 1);
        }

        // ---- Assemble Discontinuous Biharmonic Matrix & RHS ----
        gsSparseMatrix<T> K_discont;
        gsMatrix<T> F_discont;
        {
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());
            if (useBubble)
                A.options().setReal("quA", 2.0);

            auto G_map   = A.getMap(mp);
            auto u_space = A.getSpace(dbasis);

            A.initSystem();
            if (useBubble)
            {
                auto u_ex_coeff   = A.getCoeff(bubbleData.field, G_map);
                auto lap_ex_coeff = A.getCoeff(bubbleData.lapField, G_map);
                if (beta > T(0))
                    A.assemble((beta * (u_space * u_space.tr()) + ilapl(u_space, G_map) * ilapl(u_space, G_map).tr()) * meas(G_map),
                               (beta * (u_space * u_ex_coeff) + ilapl(u_space, G_map) * lap_ex_coeff) * meas(G_map));
                else
                    A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                               ilapl(u_space, G_map) * lap_ex_coeff * meas(G_map));
            }
            else if (useConstForce)
            {
                auto f_coeff = A.getCoeff(*const_rhs_ptr, G_map);
                if (beta > T(0))
                    A.assemble((beta * (u_space * u_space.tr()) + ilapl(u_space, G_map) * ilapl(u_space, G_map).tr()) * meas(G_map),
                               u_space * f_coeff * meas(G_map));
                else
                    A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                               u_space * f_coeff * meas(G_map));
            }
            else
            {
                auto f_coeff = A.getCoeff(*rhs_f_ptr, G_map);
                if (beta > T(0))
                    A.assemble((beta * (u_space * u_space.tr()) + ilapl(u_space, G_map) * ilapl(u_space, G_map).tr()) * meas(G_map),
                               u_space * f_coeff * meas(G_map));
                else
                    A.assemble(ilapl(u_space, G_map) * ilapl(u_space, G_map).tr() * meas(G_map),
                               u_space * f_coeff * meas(G_map));

                if (secondBiharmonic && !homogeneousBC)
                {
                    auto g2_coeff = A.getCoeff(*exact_lap_ptr, G_map);
                    A.assembleBdr(bc.get("Values"), (igrad(u_space, G_map) * unv(G_map)) * g2_coeff.tr() * meas(G_map));
                }
            }

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
        gsMatrix<T> sol_discont = T_free * sol_free;
        if (!homogeneousBC && nBnd > 0)
            sol_discont += T_bnd * sol_bnd;

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
        if (useBubble)
            ev.options().setReal("quA", 2.0);
        auto G_map_ev = ev.getMap(mp);
        auto u_sol_ev = ev.getVariable(sol);

        T l2err = 0, h1err = 0, h2err = 0;
        if (useBubble)
        {
            auto u_exact_ev = ev.getVariable(bubbleData.field, G_map_ev);
            l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
            h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev)).sqNorm() * meas(G_map_ev)));
            h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev)).sqNorm() * meas(G_map_ev)));
        }
        else if (useConstForce)
        {
            l2err = std::sqrt(ev.integral(u_sol_ev.sqNorm() * meas(G_map_ev)));
            h1err = std::sqrt(ev.integral(igrad(u_sol_ev, G_map_ev).sqNorm() * meas(G_map_ev)));
            h2err = std::sqrt(ev.integral(ihess(u_sol_ev, G_map_ev).sqNorm() * meas(G_map_ev)));
        }
        else
        {
            auto u_exact_ev    = ev.getVariable(*exact_u_ptr, G_map_ev);
            auto grad_exact_ev = ev.getVariable(*exact_grad_ptr, G_map_ev);
            auto hess_exact_ev = reshape(ev.getVariable(*exact_hess_ptr, G_map_ev), 2, 2);

            l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
            h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() * meas(G_map_ev)));
            h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - hess_exact_ev).sqNorm() * meas(G_map_ev)));
        }

        // Rates
        T l2_rate = 0, h1_rate = 0, h2_rate = 0;
        if (!useConstForce && ref > minRefinements && prev_h > 0)
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
        if (!useConstForce && ref > minRefinements)
            gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << l2_rate;
        else
            gsInfo << std::setw(10) << "-";

        gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
        if (!useConstForce && ref > minRefinements)
            gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << h1_rate;
        else
            gsInfo << std::setw(10) << "-";

        gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
        if (!useConstForce && ref > minRefinements)
            gsInfo << std::setw(10) << std::fixed << std::setprecision(2) << h2_rate;
        else
            gsInfo << std::setw(10) << "-";

        gsInfo << "\n" << std::flush;

        prev_h  = h_max;
        prev_l2 = l2err;
        prev_h1 = h1err;
        prev_h2 = h2err;

        // ---- Optional VTK output for solution and exact function ----
        if (plot || !outDir.empty())
        {
            const std::string tag = secondBiharmonic ? "as_g1_biharmonic2" : "as_g1_biharmonic";
            std::string solName   = prefix + tag + "_sol_r"   + std::to_string(ref);
            std::string exactName = prefix + tag + "_exact_r" + std::to_string(ref);

            gsField<T> solField(mp, sol);
            gsWriteParaview<>(solField, solName, 1000);

            if (useBubble)
            {
                gsField<T> exactField(mp, bubbleData.field, false);
                gsWriteParaview<>(exactField, exactName, 1000);
                gsInfo << " -> Exported VTK plots: " << solName << ".pvd & " << exactName << ".pvd\n";
            }
            else if (!useConstForce)
            {
                gsField<T> exactField(mp, *exact_u_ptr, false);
                gsWriteParaview<>(exactField, exactName, 1000);
                gsInfo << " -> Exported VTK plots: " << solName << ".pvd & " << exactName << ".pvd\n";
            }
            else
            {
                gsInfo << " -> Exported VTK plot: " << solName << ".pvd\n";
            }
        }
    }

    gsInfo << std::string(99, '-') << "\n";
    gsInfo << "Done.\n";
    return 0;
}
