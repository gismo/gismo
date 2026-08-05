/** @file as_g1_biharmonic_multigrid_example.cpp
 *
 *  @brief Analysis-Suitable G1 (AS-G1) geometric multigrid solver for the
 *         biharmonic equation, with convergence analysis across refinement
 *         levels and timestamped VTK plotting.
 *
 *  The manufactured solution is u(x,y) = sin(a*pi*x)*cos(a*pi*y). Since this
 *  does not vanish on the domain boundary, INHOMOGENEOUS clamped boundary
 *  conditions (trace u and normal derivative grad(u).n from the exact solution)
 *  are prescribed via an L2 projection onto the AS-G1 space; the interior
 *  problem is solved with the generic gsMultiGridOp-preconditioned CG provided
 *  by gsAsG1BiharmonicMG (gsMultiGrid/gsAsG1BiharmonicMG.h), which reuses the
 *  shared AS-G1 DOF coupling and gluing-data machinery instead of a
 *  hand-rolled V-cycle.
 *
 *  This file is part of the G+Smo library.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Author(s): F. Hasanova, S. Takacs
 */

#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <gismo.h>

#include <gsMultiGrid/gsAsG1BiharmonicMG.h>

using namespace gismo;
using namespace gismo::expr;

// ISO-like timestamp string YYYYMMDD_HHMMSS for the default output folder.
static std::string getTimestampString()
{
    std::time_t now_c = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::now());
    struct tm *parts = std::localtime(&now_c);
    std::ostringstream oss;
    oss << (1900 + parts->tm_year)
        << std::setfill('0') << std::setw(2) << (1 + parts->tm_mon)
        << std::setfill('0') << std::setw(2) << parts->tm_mday << "_"
        << std::setfill('0') << std::setw(2) << parts->tm_hour
        << std::setfill('0') << std::setw(2) << parts->tm_min
        << std::setfill('0') << std::setw(2) << parts->tm_sec;
    return oss.str();
}

int main(int argc, char *argv[])
{
    using T = real_t;

    std::string geometry("domain2d/2patch/weirdo_multivalence_non_bilinear.xml");
    std::string outDir("");
    index_t degree = 3;
    index_t minRefinements = 2;
    index_t maxRefinements = 4;
    index_t numSmooth = 3;
    index_t freqA = 1;
    bool wcycle = false;

    gsCmdLine cmd("AS-G1 Multigrid Biharmonic Convergence & Plotting Tool.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
    cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
    cmd.addInt("m", "minRefinements", "Minimum uniform refinement level (coarsest level, min 2).", minRefinements);
    cmd.addInt("r", "maxRefinements", "Maximum uniform refinement level (finest level).", maxRefinements);
    cmd.addInt("s", "smooth", "Number of pre-/post-smoothing steps.", numSmooth);
    cmd.addInt("a", "frequency", "Frequency integer factor 'a' in sin(a*pi*x)*cos(a*pi*y).", freqA);
    cmd.addSwitch("wcycle", "Use W-cycles (numCycles = 2) instead of V-cycles.", wcycle);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    if (degree < 3) degree = 3;
    if (minRefinements < 2) minRefinements = 2;
    if (maxRefinements < minRefinements) maxRefinements = minRefinements;
    if (numSmooth < 1) numSmooth = 1;
    if (freqA < 1) freqA = 1;

    if (outDir.empty())
        outDir = "experiments/as_g1_multigrid_" + getTimestampString();
    if (outDir.back() != '/') outDir += '/';
    gsFileManager::mkdir(outDir);

    gsInfo << "\n=========================================================================================================\n";
    gsInfo << "  AS-G1 Geometric Multigrid Biharmonic Convergence & Benchmark Suite\n";
    gsInfo << "  Target Function: u(x, y) = sin(" << freqA << "*pi*x) * cos(" << freqA << "*pi*y)\n";
    gsInfo << "  Boundary:        inhomogeneous clamped (u and grad(u).n from exact solution)\n";
    gsInfo << "  Solver:          gsMultiGridOp (Galerkin, symmetric Gauss-Seidel, "
           << numSmooth << " smoothing steps) + CG, " << (wcycle ? "W-cycle" : "V-cycle") << "\n";
    gsInfo << "  Geometry File:   " << geometry << "\n";
    gsInfo << "  Refinement:      r_min = " << minRefinements << " to r_max = " << maxRefinements << "\n";
    gsInfo << "  Output Dir:      " << outDir << "\n";
    gsInfo << "=========================================================================================================\n\n";

    // Manufactured solution and its (physical) Laplacian.
    const std::string s_a = std::to_string(freqA);
    gsFunctionExpr<T> exact_u("sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)", 2);
    gsFunctionExpr<T> lap_u("-2*" + s_a + "^2*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)", 2);

    gsMultiPatch<T>::uPtr basePtr = gsReadFile<>(geometry);
    if (!basePtr) {
        gsInfo << "Error: Cannot read geometry " << geometry << "\n";
        return -1;
    }
    gsMultiPatch<T> base = *basePtr;
    base.computeTopology();

    // Table header.
    gsInfo << std::setw(5)  << "r"
           << std::setw(12) << "h_max"
           << std::setw(10) << "N_free"
           << std::setw(14) << "L2 Error" << std::setw(9) << "L2 Rate"
           << std::setw(14) << "H1 Error" << std::setw(9) << "H1 Rate"
           << std::setw(14) << "H2 Error" << std::setw(9) << "H2 Rate"
           << std::setw(10) << "MG Iters"
           << std::setw(14) << "Solve Time"
           << "\n";
    gsInfo << std::string(118, '-') << "\n" << std::flush;

    T prev_h = 0, prev_l2 = 0, prev_h1 = 0, prev_h2 = 0;

    for (index_t maxRef = minRefinements; maxRef <= maxRefinements; ++maxRef) {
        try {
            gsAsG1BiharmonicMG<T> solver(base, degree, minRefinements, maxRef, numSmooth);
            solver.setNumSmooth(numSmooth, numSmooth);
            if (wcycle) solver.setNumCycles(2);

            const gsMultiPatch<T> &mp = solver.fineMultiPatch();
            gsMultiBasis<T> dbasis(mp);
            const T h_max = std::pow(0.5, maxRef);

            gsPiecewiseFunction<T> uField, lapField;
            for (size_t i = 0; i < mp.nPatches(); ++i) {
                uField.addPiece(exact_u);
                lapField.addPiece(lap_u);
            }

            // Prescribe inhomogeneous boundary values from the exact solution.
            gsMatrix<T> cBnd = solver.projectBoundaryValues(exact_u);

            // Energy-form RHS on interior test functions: (Delta u_ex, Delta v).
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());
            auto G_map = A.getMap(mp);
            auto u_space = A.getSpace(dbasis);
            auto lap_ex = A.getCoeff(lapField, G_map);
            A.initSystem();
            A.assemble(ilapl(u_space, G_map) * lap_ex * meas(G_map));
            gsMatrix<T> F_free = solver.transformFine().transpose() * A.rhs();
            gsMatrix<T> F_lifted = solver.liftedRhs(F_free, cBnd);

            // Solve interior system via multigrid-preconditioned CG.
            gsMatrix<T> c_free;
            auto t0_solve = std::chrono::high_resolution_clock::now();
            index_t mgIters = solver.solve(F_lifted, c_free, 200, T(1e-8));
            auto t1_solve = std::chrono::high_resolution_clock::now();
            double solveTime_ms = std::chrono::duration<double, std::milli>(t1_solve - t0_solve).count();

            // Reconstruct field (interior + prescribed boundary) and evaluate error.
            gsMultiPatch<T> solField = solver.reconstruct(c_free, cBnd);
            gsExprEvaluator<T> ev;
            ev.setIntegrationDomain(dbasis.domain());
            auto G_map_ev = ev.getMap(mp);
            auto u_exact_ev = ev.getVariable(uField, G_map_ev);
            auto u_sol_ev = ev.getVariable(solField);

            const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
            const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - igrad(u_exact_ev)).sqNorm() * meas(G_map_ev)));
            const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - ihess(u_exact_ev)).sqNorm() * meas(G_map_ev)));

            T l2_rate = 0, h1_rate = 0, h2_rate = 0;
            if (maxRef > minRefinements && prev_h > 0) {
                const T q = std::log(prev_h / h_max);
                l2_rate = std::log(prev_l2 / l2err) / q;
                h1_rate = std::log(prev_h1 / h1err) / q;
                h2_rate = std::log(prev_h2 / h2err) / q;
            }

            gsInfo << std::setw(5)  << maxRef
                   << std::setw(12) << std::scientific << std::setprecision(4) << h_max
                   << std::setw(10) << solver.freeDofs()
                   << std::setw(14) << std::scientific << std::setprecision(4) << l2err;
            if (maxRef > minRefinements) gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << l2_rate; else gsInfo << std::setw(9) << "-";
            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
            if (maxRef > minRefinements) gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h1_rate; else gsInfo << std::setw(9) << "-";
            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
            if (maxRef > minRefinements) gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h2_rate; else gsInfo << std::setw(9) << "-";
            gsInfo << std::setw(10) << mgIters
                   << std::setw(14) << std::fixed << std::setprecision(2) << solveTime_ms << " ms\n" << std::flush;

            prev_h = h_max; prev_l2 = l2err; prev_h1 = h1err; prev_h2 = h2err;

            // Export VTK plots to the output folder.
            gsField<T> solF(mp, solField);
            gsWriteParaview<>(solF, outDir + "as_g1_multigrid_sol_r" + std::to_string(maxRef), 1000);
            gsField<T> exactF(mp, uField, false);
            gsWriteParaview<>(exactF, outDir + "as_g1_multigrid_exact_r" + std::to_string(maxRef), 1000);

        } catch (const std::exception &e) {
            gsInfo << "ERROR on level " << maxRef << ": " << e.what() << "\n" << std::flush;
        }
    }

    gsInfo << std::string(118, '-') << "\n";
    gsInfo << "VTK plot files exported to directory: " << outDir << "\n";
    gsInfo << "=========================================================================================================\n";

    return 0;
}
