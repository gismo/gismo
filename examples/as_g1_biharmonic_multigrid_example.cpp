/** @file as_g1_biharmonic_multigrid_example.cpp
 *
 *  @brief Example demonstrating Analysis-Suitable G1 (AS-G1) Multigrid Solver for Biharmonic Equations,
 *         including convergence analysis across refinement levels and timestamped VTK plotting.
 *
 *  This file is part of the G+Smo library.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Author(s): Antigravity AI, F. Hasanova, S. Takacs
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <sstream>

#include <gismo.h>
#include <gsMultiGrid/gsAsG1Multigrid.h>

using namespace gismo;
using namespace gismo::expr;

// Helper to generate ISO-like timestamp string: YYYYMMDD_HHMMSS
std::string getTimestampString() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
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

template<typename T>
index_t solvePCG(const gsSparseMatrix<T>& A,
                 const gsAsG1Multigrid<T>& mg,
                 const gsMatrix<T>& b,
                 gsMatrix<T>& x,
                 index_t maxIter = 100,
                 T tol = 1e-8,
                 index_t finestLevel = 1)
{
    const index_t n = A.rows();
    x = gsMatrix<T>::Zero(n, 1);

    gsMatrix<T> r = b - A * x;
    T r0_norm = r.norm();
    if (r0_norm < 1e-15) return 0;

    gsMatrix<T> z = gsMatrix<T>::Zero(n, 1);
    mg.vCycle(finestLevel, z, r, 2, 2);

    gsMatrix<T> p = z;
    T rz = (r.transpose() * z)(0, 0);

    for (index_t iter = 1; iter <= maxIter; ++iter)
    {
        gsMatrix<T> Ap = A * p;
        T pAp = (p.transpose() * Ap)(0, 0);
        if (std::abs(pAp) < 1e-18) return iter;

        T alpha = rz / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        T r_norm = r.norm();
        T relRes = r_norm / r0_norm;

        if (relRes < tol)
        {
            return iter;
        }

        gsMatrix<T> z_new = gsMatrix<T>::Zero(n, 1);
        mg.vCycle(finestLevel, z_new, r, 2, 2);

        T rz_new = (r.transpose() * z_new)(0, 0);
        T beta = rz_new / rz;

        p = z_new + beta * p;
        rz = rz_new;
    }

    return maxIter;
}

int main(int argc, char *argv[]) {
    using T = real_t;

    std::string geometry("domain2d/2patch/multipatch/weirdo_multivalence_non_bilinear.xml");
    std::string outDir("");
    index_t degree = 3;
    index_t minRefinements = 2;
    index_t maxRefinements = 4;
    index_t freqA = 1;
    bool plot = false;

    gsCmdLine cmd("AS-G1 Multigrid Biharmonic Convergence & Plotting Tool.");
    cmd.addString("f", "file", "Multi-patch geometry XML file.", geometry);
    cmd.addString("o", "outDir", "Output directory for VTK/PVD files.", outDir);
    cmd.addInt("d", "degree", "Target degree (minimum 3).", degree);
    cmd.addInt("m", "minRefinements", "Minimum uniform refinement level (coarsest level, min 2).", minRefinements);
    cmd.addInt("r", "maxRefinements", "Maximum uniform refinement level (finest level).", maxRefinements);
    cmd.addInt("a", "frequency", "Frequency integer factor 'a' in sin(a*pi*x)*cos(a*pi*y).", freqA);
    cmd.addSwitch("plot", "Export target function and numerical solution to VTK files.", plot);

    try {
        cmd.getValues(argc, argv);
    } catch (int rv) {
        return rv;
    }

    // Set up output directory inside experiments/
    if (outDir.empty()) {
        outDir = "experiments/as_g1_multigrid_" + getTimestampString();
    }
    if (outDir.back() != '/') outDir += '/';
    gsFileManager::mkdir(outDir);

    gsInfo << "\n=========================================================================================================\n";
    gsInfo << "  AS-G1 Geometric Multigrid Biharmonic Convergence & Benchmark Suite\n";
    gsInfo << "  Target Function: u(x, y) = sin(" << freqA << "*pi*x) * cos(" << freqA << "*pi*y)\n";
    gsInfo << "  Geometry File:   " << geometry << "\n";
    gsInfo << "  Refinement:      r_min = " << minRefinements << " to r_max = " << maxRefinements << "\n";
    gsInfo << "  Output Dir:      " << outDir << "\n";
    gsInfo << "=========================================================================================================\n\n";

    // Setup Target Exact Functions
    const std::string s_a = std::to_string(freqA);
    const std::string u_expr = "sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
    const std::string grad_x_expr = s_a + "*pi*cos(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
    const std::string grad_y_expr = "-" + s_a + "*pi*sin(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
    const std::string rhs_expr = std::to_string(4 * freqA * freqA * freqA * freqA) + "*pi^4*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";

    const std::string hess_xx = "-" + std::to_string(freqA * freqA) + "*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";
    const std::string hess_xy = "-" + std::to_string(freqA * freqA) + "*pi^2*cos(" + s_a + "*pi*x)*sin(" + s_a + "*pi*y)";
    const std::string hess_yy = "-" + std::to_string(freqA * freqA) + "*pi^2*sin(" + s_a + "*pi*x)*cos(" + s_a + "*pi*y)";

    gsFunctionExpr<T> exact_u(u_expr, 2);
    gsFunctionExpr<T> exact_grad(grad_x_expr, grad_y_expr, 2);
    gsFunctionExpr<T> exact_hess(hess_xx, hess_xy, hess_xy, hess_yy, 2);
    gsFunctionExpr<T> rhs_f(rhs_expr, 2);

    // Table Header
    gsInfo << std::setw(5)  << "r"
           << std::setw(12) << "h_max"
           << std::setw(10) << "N_free"
           << std::setw(14) << "L2 Error" << std::setw(9) << "L2 Rate"
           << std::setw(14) << "H1 Error" << std::setw(9) << "H1 Rate"
           << std::setw(14) << "H2 Error" << std::setw(9) << "H2 Rate"
           << std::setw(10) << "PCG Iters"
           << std::setw(14) << "Solve Time"
           << "\n";
    gsInfo << std::string(118, '-') << "\n" << std::flush;

    T prev_h = 0;
    T prev_l2 = 0;
    T prev_h1 = 0;
    T prev_h2 = 0;

    for (index_t maxRef = minRefinements; maxRef <= maxRefinements; ++maxRef) {
        try {
            gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
            if (!mpPtr) {
                gsInfo << "Error: Cannot read geometry " << geometry << "\n";
                return -1;
            }

            // Construct Multigrid solver for level 2 -> maxRef
            auto t0_setup = std::chrono::high_resolution_clock::now();
            gsAsG1Multigrid<T> mgSolver(*mpPtr, degree, minRefinements, maxRef);
            auto t1_setup = std::chrono::high_resolution_clock::now();
            double setupTime_ms = std::chrono::duration<double, std::milli>(t1_setup - t0_setup).count();

            index_t finestLevel = mgSolver.numLevels() - 1;
            const gsSparseMatrix<T>& K_fine = mgSolver.stiffnessMatrix(finestLevel);
            index_t nFree = mgSolver.freeDofs(finestLevel);
            const T h_max = std::pow(0.5, maxRef);

            // Re-get finest geometry map & bases
            gsMultiPatch<T> mpFine = *mpPtr;
            mpFine.computeTopology();
            if (mpFine.patch(0).basis().degree(0) < degree) {
                mpFine.degreeElevate(degree - mpFine.patch(0).basis().degree(0));
            }
            const short_t deg = mpFine.patch(0).basis().degree(0);
            const index_t mult = std::max<index_t>(deg - 1, 1);
            for (index_t i = 0; i < maxRef; ++i) mpFine.uniformRefine(1, mult);

            // Assemble RHS vector b_free
            gsMultiBasis<T> dbasis(mpFine);
            gsExprAssembler<T> A(1, 1);
            A.setIntegrationDomain(dbasis.domain());
            auto G_map = A.getMap(mpFine);
            auto u_space = A.getSpace(dbasis);
            auto f_coeff = A.getCoeff(rhs_f, G_map);

            A.initSystem();
            A.assemble(u_space * f_coeff * meas(G_map));
            const gsMatrix<T>& F_disjoint = A.rhs();

            // Project RHS to AS-G1 free DOFs
            gsMatrix<T> gd = computeGluingData(mpFine, T(1e-8), 0);
            std::vector<gsArgyrisEmbedding<T>> argBasis;
            gsVector<index_t> patchDofSizes(mpFine.nPatches());
            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                argBasis.push_back(deriveArgyrisBasisEmbedding(
                    dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mpFine.patch(i).basis()),
                    gsMatrix<T>(gd.row(i)), mpFine.patch(i)));
                patchDofSizes[i] = argBasis[i].matrix.cols();
            }

            gsDofMapper mapper(patchDofSizes);
            for (auto it = mpFine.iBegin(); it != mpFine.iEnd(); ++it) {
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

                for (index_t j1 = 0; j1 < nLvl0; ++j1) {
                    const index_t j2 = flipped ? nLvl0 - 1 - j1 : j1;
                    mapper.matchDof(ps1.patch, offLvl0_1 + j1, ps2.patch, offLvl0_2 + j2);
                }
                for (index_t j1 = 0; j1 < nLvl1; ++j1) {
                    const index_t j2 = flipped ? nLvl1 - 1 - j1 : j1;
                    mapper.matchDof(ps1.patch, offLvl1_1 + j1, ps2.patch, offLvl1_2 + j2);
                }

                std::vector<patchCorner> corners;
                ps1.getContainedCorners(2, corners);
                for (index_t i = 0; i < 2; ++i) {
                    const index_t c1 = corners[i].m_index - 1;
                    const index_t c2 = ifc.mapCorner(corners[i]).m_index - 1;
                    const index_t off_corner_1 = sumUntil(argBasis[ps1.patch].sizes, 9 + c1);
                    const index_t off_corner_2 = sumUntil(argBasis[ps2.patch].sizes, 9 + c2);
                    for (index_t k = 0; k < 6; ++k)
                        mapper.matchDof(ps1.patch, off_corner_1 + k, ps2.patch, off_corner_2 + k);
                }
            }

            std::set<std::pair<size_t, index_t>> ifcSides;
            for (auto it = mpFine.iBegin(); it != mpFine.iEnd(); ++it) {
                ifcSides.insert({it->first().patch, it->first().side()});
                ifcSides.insert({it->second().patch, it->second().side()});
            }

            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                for (boxSide side = boxSide::getFirst(2); side != boxSide::getEnd(2); ++side) {
                    if (ifcSides.find({i, side.m_index}) != ifcSides.end()) continue;
                    const index_t nLvl0 = argBasis[i].sizes[1 + 2 * (side.m_index - 1)];
                    const index_t offLvl0 = sumUntil(argBasis[i].sizes, 1 + 2 * (side.m_index - 1));
                    for (index_t j = 0; j < nLvl0; ++j) mapper.eliminateDof(offLvl0 + j, i);
                    const index_t nLvl1 = argBasis[i].sizes[2 + 2 * (side.m_index - 1)];
                    const index_t offLvl1 = sumUntil(argBasis[i].sizes, 2 + 2 * (side.m_index - 1));
                    for (index_t j = 0; j < nLvl1; ++j) mapper.eliminateDof(offLvl1 + j, i);

                    patchSide ps(i, side);
                    std::vector<patchCorner> corners;
                    ps.getContainedCorners(2, corners);
                    for (const auto &c : corners) {
                        const index_t cIdx = c.m_index - 1;
                        const index_t offCorner = sumUntil(argBasis[i].sizes, 9 + cIdx);
                        for (index_t k = 0; k < 6; ++k) mapper.eliminateDof(offCorner + k, i);
                    }
                }
            }
            mapper.finalize();

            index_t nDisjoint = 0;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) nDisjoint += argBasis[i].matrix.rows();

            gsSparseEntries<T> tFreeEntries;
            index_t rowOffset = 0;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                const gsSparseMatrix<T> &Ai = argBasis[i].matrix;
                const index_t patchSize = mapper.patchSize(i);
                for (index_t j = 0; j < patchSize; ++j) {
                    if (mapper.is_boundary(j, i)) continue;
                    index_t gIdx = mapper.index(j, i);
                    for (typename gsSparseMatrix<T>::InnerIterator it(Ai, j); it; ++it) {
                        tFreeEntries.add(rowOffset + it.row(), gIdx, it.value());
                    }
                }
                rowOffset += Ai.rows();
            }

            gsSparseMatrix<T> S_fine(nDisjoint, nFree);
            S_fine.setFrom(tFreeEntries);

            gsMatrix<T> b_free = S_fine.transpose() * F_disjoint;

            // Solve via Multigrid PCG
            gsMatrix<T> u_pcg(nFree, 1);
            auto t0_solve = std::chrono::high_resolution_clock::now();
            index_t pcgIters = solvePCG(K_fine, mgSolver, b_free, u_pcg, 100, 1e-8, finestLevel);
            auto t1_solve = std::chrono::high_resolution_clock::now();
            double solveTime_ms = std::chrono::duration<double, std::milli>(t1_solve - t0_solve).count();

            // Reconstruct multi-patch field for error analysis & plotting
            gsMatrix<T> c_disjoint = S_fine * u_pcg;
            gsMultiPatch<T> solField;
            index_t cOffset = 0;
            for (size_t i = 0; i < mpFine.nPatches(); ++i) {
                const index_t sz = argBasis[i].matrix.rows();
                gsMatrix<T> ci = c_disjoint.block(cOffset, 0, sz, 1);
                cOffset += sz;
                const gsTensorBSplineBasis<2, T> &tb = dynamic_cast<const gsTensorBSplineBasis<2, T> &>(mpFine.patch(i).basis());
                solField.addPatch(tb.makeGeometry(give(ci)));
            }

            // Error Evaluation
            gsExprEvaluator<T> ev;
            ev.setIntegrationDomain(dbasis.domain());
            auto G_map_ev = ev.getMap(mpFine);

            auto u_exact_ev = ev.getVariable(exact_u, G_map_ev);
            auto grad_exact_ev = ev.getVariable(exact_grad, G_map_ev);
            auto hess_exact_ev = reshape(ev.getVariable(exact_hess, G_map_ev), 2, 2);
            auto u_sol_ev = ev.getVariable(solField);

            const T l2err = std::sqrt(ev.integral((u_sol_ev - u_exact_ev).sqNorm() * meas(G_map_ev)));
            const T h1err = std::sqrt(ev.integral((igrad(u_sol_ev, G_map_ev) - grad_exact_ev.tr()).sqNorm() * meas(G_map_ev)));
            const T h2err = std::sqrt(ev.integral((ihess(u_sol_ev, G_map_ev) - hess_exact_ev).sqNorm() * meas(G_map_ev)));

            // Convergence Rates
            T l2_rate = 0, h1_rate = 0, h2_rate = 0;
            if (maxRef > minRefinements && prev_h > 0) {
                const T h_ratio = prev_h / h_max;
                l2_rate = std::log(prev_l2 / l2err) / std::log(h_ratio);
                h1_rate = std::log(prev_h1 / h1err) / std::log(h_ratio);
                h2_rate = std::log(prev_h2 / h2err) / std::log(h_ratio);
            }

            // Print Row
            gsInfo << std::setw(5)  << maxRef
                   << std::setw(12) << std::scientific << std::setprecision(4) << h_max
                   << std::setw(10) << nFree
                   << std::setw(14) << std::scientific << std::setprecision(4) << l2err;
            if (maxRef > minRefinements)
                gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << l2_rate;
            else
                gsInfo << std::setw(9) << "-";

            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h1err;
            if (maxRef > minRefinements)
                gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h1_rate;
            else
                gsInfo << std::setw(9) << "-";

            gsInfo << std::setw(14) << std::scientific << std::setprecision(4) << h2err;
            if (maxRef > minRefinements)
                gsInfo << std::setw(9) << std::fixed << std::setprecision(2) << h2_rate;
            else
                gsInfo << std::setw(9) << "-";

            gsInfo << std::setw(10) << pcgIters
                   << std::setw(14) << std::fixed << std::setprecision(2) << solveTime_ms << " ms\n" << std::flush;

            prev_h = h_max;
            prev_l2 = l2err;
            prev_h1 = h1err;
            prev_h2 = h2err;

            // Export VTK plots to timestamped output folder
            std::string solName = outDir + "as_g1_multigrid_sol_r" + std::to_string(maxRef);
            std::string exactName = outDir + "as_g1_multigrid_exact_r" + std::to_string(maxRef);

            gsField<T> solF(mpFine, solField);
            gsWriteParaview<>(solF, solName, 1000);

            gsField<T> exactF(mpFine, exact_u, false);
            gsWriteParaview<>(exactF, exactName, 1000);

        } catch (const std::exception &e) {
            gsInfo << "ERROR on level " << maxRef << ": " << e.what() << "\n" << std::flush;
        }
    }

    gsInfo << std::string(118, '-') << "\n";
    gsInfo << "VTK plot files exported to directory: " << outDir << "\n";
    gsInfo << "=========================================================================================================\n";

    return 0;
}
