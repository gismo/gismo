/** @file as_g1_biharmonic_multigrid_sincos_example.cpp
 *
 *  @brief AS-G1 biharmonic multigrid solver for the manufactured solution
 *         u(x,y) = sin(a*pi*x) * cos(a*pi*y) with INHOMOGENEOUS clamped
 *         boundary conditions.
 *
 *  Unlike as_g1_biharmonic_multigrid_v2_example.cpp (which uses a bubble that
 *  vanishes on the boundary and hence homogeneous BCs), this example prescribes
 *  the boundary values from the manufactured solution itself: on the domain
 *  boundary both the trace u and the normal derivative grad(u).n are set to the
 *  values of sin(a*pi*x)*cos(a*pi*y). The boundary DOFs are obtained by an L2
 *  projection of the exact solution onto the finest AS-G1 space (carrying the
 *  trace and normal-derivative data); the interior problem is then solved with
 *  the gsMultiGridOp-preconditioned CG from gsAsG1BiharmonicMG.
 *
 *  The multigrid cycle type (V- or W-cycle) is selectable and the numerical and
 *  exact solutions can be exported as PVD/VTK to a chosen folder.
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
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <gismo.h>

#include <gsMultiGrid/gsAsG1BiharmonicMG.h>

using namespace gismo;
using namespace gismo::expr;

int main(int argc, char* argv[])
{
    using T = real_t;

    std::string dir("filedata/domain2d/2patch/multipatch/");
    std::string file("");
    std::string outDir("");
    index_t degree = 3;
    index_t minRef = 2;
    index_t maxRef = 4;
    index_t numSmooth = 3;
    index_t freqA = 1;
    bool wcycle = false;

    gsCmdLine cmd("AS-G1 biharmonic multigrid solver for sin(a*pi*x)*cos(a*pi*y) "
                  "with inhomogeneous clamped boundary conditions.");
    cmd.addString("d", "dir", "Directory of multipatch geometry XML files.", dir);
    cmd.addString("f", "file", "Solve a single geometry XML file instead of the whole directory.", file);
    cmd.addString("o", "outDir", "Output folder for PVD/VTK plots (empty = no plotting).", outDir);
    cmd.addInt("p", "degree", "Discretization degree (minimum 3).", degree);
    cmd.addInt("m", "minRef", "Coarsest refinement level (>= 2).", minRef);
    cmd.addInt("r", "maxRef", "Finest refinement level.", maxRef);
    cmd.addInt("s", "smooth", "Number of pre-/post-smoothing steps.", numSmooth);
    cmd.addInt("a", "frequency", "Frequency factor a in sin(a*pi*x)*cos(a*pi*y).", freqA);
    cmd.addSwitch("wcycle", "Use W-cycles (numCycles = 2) instead of V-cycles.", wcycle);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (degree < 3)  degree = 3;
    if (minRef < 2)  minRef = 2;
    if (maxRef < minRef) maxRef = minRef;
    if (numSmooth < 1) numSmooth = 1;
    if (freqA < 1) freqA = 1;
    if (dir.back() != '/') dir += '/';
    if (!outDir.empty())
    {
        if (outDir.back() != '/') outDir += '/';
        gsFileManager::mkdir(outDir);
    }

    // Manufactured solution and its Laplacian (physical coordinates).
    const std::string sa = std::to_string(freqA);
    const std::string u_expr   = "sin(" + sa + "*pi*x)*cos(" + sa + "*pi*y)";
    const std::string lap_expr = "-2*" + sa + "^2*pi^2*sin(" + sa + "*pi*x)*cos(" + sa + "*pi*y)";
    gsFunctionExpr<T> exact_u(u_expr, 2);
    gsFunctionExpr<T> lap_u(lap_expr, 2);

    // Geometry list.
    std::vector<std::string> geoms;
    if (!file.empty())
        geoms.push_back(file);
    else
    {
        const char* names[] = {
            "2_patch_rectangle.xml",
            "2_patch_rectangle_non_bilinear.xml",
            "two_patch_mathematica.xml",
            "two_patch_mathematica_non_bilinear.xml",
            "three_patch_belt.xml",
            "three_patch_belt_non_bilinear.xml",
            "4_patch_3_valence_not_bilinear.xml",
            "5_patch_5_valence_not_bilinear.xml",
            "square_6_patch.xml",
            "square_6_patch_non_bilinear.xml",
            "6_patch_4_valence_each_not_bilinear.xml",
            "6_patch_6_valence_not_bilinear.xml",
            "experiment_multivalence_not_bilinear.xml",
            "weirdo_multivalence.xml",
            "weirdo_multivalence_non_bilinear.xml"
        };
        for (const char* n : names)
        {
            std::string full = dir + n;
            if (gsFileManager::fileExists(full)) geoms.push_back(full);
        }
    }

    gsInfo << "\n=====================================================================================================\n";
    gsInfo << "  AS-G1 Biharmonic Multigrid Solver -- sin(" << freqA << "*pi*x)*cos(" << freqA << "*pi*y)\n";
    gsInfo << "  Boundary     : inhomogeneous clamped (u and grad(u).n from exact solution)\n";
    gsInfo << "  Solver       : gsMultiGridOp (Galerkin, symmetric Gauss-Seidel, "
           << numSmooth << " smoothing steps) + CG\n";
    gsInfo << "  Cycle type   : " << (wcycle ? "W-cycle" : "V-cycle") << "\n";
    gsInfo << "  Levels       : r = " << minRef << " .. " << maxRef
           << "   Degree = " << degree << "\n";
    gsInfo << "  Geometries   : " << geoms.size() << "\n";
    if (!outDir.empty()) gsInfo << "  Plot folder  : " << outDir << "\n";
    gsInfo << "=====================================================================================================\n";

    index_t nConverged = 0, nTotal = 0;
    struct Summary { std::string name; index_t iters; T relres; T h2; bool ok; };
    std::vector<Summary> summaries;

    for (const std::string& geom : geoms)
    {
        std::string shortName = geom;
        size_t slash = shortName.find_last_of('/');
        if (slash != std::string::npos) shortName = shortName.substr(slash + 1);
        std::string stem = shortName;
        size_t dot = stem.find_last_of('.');
        if (dot != std::string::npos) stem = stem.substr(0, dot);

        gsInfo << "\n-----------------------------------------------------------------------------------------------------\n";
        gsInfo << "Geometry: " << shortName << "\n";

        gsMultiPatch<T>::uPtr basePtr = gsReadFile<>(geom);
        if (!basePtr)
        {
            gsInfo << "  Cannot read geometry, skipping.\n";
            continue;
        }
        gsMultiPatch<T> base = *basePtr;
        base.computeTopology();

        gsInfo << std::setw(5)  << "r"
               << std::setw(10) << "N_free"
               << std::setw(10) << "N_bnd"
               << std::setw(13) << "L2 Error"  << std::setw(8) << "rate"
               << std::setw(13) << "H1 Error"  << std::setw(8) << "rate"
               << std::setw(13) << "H2 Error"  << std::setw(8) << "rate"
               << std::setw(8)  << "MG-it"
               << std::setw(12) << "rel.res."
               << std::setw(12) << "solve[ms]"
               << "\n";
        gsInfo << std::string(112, '-') << "\n" << std::flush;

        T ph = 0, pl2 = 0, ph1 = 0, ph2 = 0;
        Summary lastSum{shortName, 0, 0, 0, false};

        for (index_t r = minRef; r <= maxRef; ++r)
        {
            try
            {
                gsAsG1BiharmonicMG<T> solver(base, degree, minRef, r, numSmooth);
                solver.setNumSmooth(numSmooth, numSmooth);
                if (wcycle) solver.setNumCycles(2);

                const gsMultiPatch<T>& mp = solver.fineMultiPatch();
                gsMultiBasis<T> dbasis(mp);

                gsPiecewiseFunction<T> uField, lapField;
                for (size_t i = 0; i < mp.nPatches(); ++i)
                {
                    uField.addPiece(exact_u);
                    lapField.addPiece(lap_u);
                }

                // Prescribe inhomogeneous boundary values from the exact solution.
                gsMatrix<T> cBnd = solver.projectBoundaryValues(exact_u);

                // Energy-form RHS on interior test functions: (Delta u_ex, Delta v).
                gsExprAssembler<T> A(1, 1);
                A.setIntegrationDomain(dbasis.domain());
                auto G = A.getMap(mp);
                auto u = A.getSpace(dbasis);
                auto lap_ex = A.getCoeff(lapField, G);
                A.initSystem();
                A.assemble(ilapl(u, G) * lap_ex * meas(G));
                gsMatrix<T> F_free = solver.transformFine().transpose() * A.rhs();

                // Move the prescribed boundary contribution to the right-hand side.
                gsMatrix<T> F_lifted = solver.liftedRhs(F_free, cBnd);

                // Solve interior system with multigrid-preconditioned CG.
                gsMatrix<T> c_free;
                T relres = 0;
                auto t0 = std::chrono::high_resolution_clock::now();
                index_t iters = solver.solve(F_lifted, c_free, 200, T(1e-8), &relres);
                auto t1 = std::chrono::high_resolution_clock::now();
                double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

                // Reconstruct field (interior + prescribed boundary) and evaluate error.
                gsMultiPatch<T> sol = solver.reconstruct(c_free, cBnd);
                gsExprEvaluator<T> ev;
                ev.setIntegrationDomain(dbasis.domain());
                auto Gev = ev.getMap(mp);
                auto uex = ev.getVariable(uField, Gev);
                auto ush = ev.getVariable(sol);
                const T l2 = std::sqrt(ev.integral((ush - uex).sqNorm() * meas(Gev)));
                const T h1 = std::sqrt(ev.integral((igrad(ush, Gev) - igrad(uex)).sqNorm() * meas(Gev)));
                const T h2 = std::sqrt(ev.integral((ihess(ush, Gev) - ihess(uex)).sqNorm() * meas(Gev)));

                const T h = std::pow(0.5, r);
                T rl2 = 0, rh1 = 0, rh2 = 0;
                if (r > minRef && ph > 0)
                {
                    const T q = std::log(ph / h);
                    rl2 = std::log(pl2 / l2) / q;
                    rh1 = std::log(ph1 / h1) / q;
                    rh2 = std::log(ph2 / h2) / q;
                }

                gsInfo << std::setw(5) << r
                       << std::setw(10) << solver.freeDofs()
                       << std::setw(10) << solver.boundaryDofs()
                       << std::setw(13) << std::scientific << std::setprecision(3) << l2;
                if (r > minRef) gsInfo << std::setw(8) << std::fixed << std::setprecision(2) << rl2; else gsInfo << std::setw(8) << "-";
                gsInfo << std::setw(13) << std::scientific << std::setprecision(3) << h1;
                if (r > minRef) gsInfo << std::setw(8) << std::fixed << std::setprecision(2) << rh1; else gsInfo << std::setw(8) << "-";
                gsInfo << std::setw(13) << std::scientific << std::setprecision(3) << h2;
                if (r > minRef) gsInfo << std::setw(8) << std::fixed << std::setprecision(2) << rh2; else gsInfo << std::setw(8) << "-";
                gsInfo << std::setw(8) << iters
                       << std::setw(12) << std::scientific << std::setprecision(2) << relres
                       << std::setw(12) << std::fixed << std::setprecision(1) << ms
                       << "\n" << std::flush;

                ph = h; pl2 = l2; ph1 = h1; ph2 = h2;
                lastSum = Summary{shortName, iters, relres, h2, (relres < 1e-8)};

                // Optional PVD/VTK export of numerical and exact solutions.
                if (!outDir.empty())
                {
                    const std::string bname = outDir + "sincos_" + stem + "_r" + std::to_string(r);
                    gsField<T> solPlot(mp, sol);
                    gsWriteParaview<>(solPlot, bname + "_sol", 1000);
                    gsField<T> exPlot(mp, uField, false);
                    gsWriteParaview<>(exPlot, bname + "_exact", 1000);
                }
            }
            catch (const std::exception& e)
            {
                gsInfo << std::setw(5) << r << "   ERROR: " << e.what() << "\n" << std::flush;
            }
        }

        ++nTotal;
        if (lastSum.ok) ++nConverged;
        summaries.push_back(lastSum);
    }

    gsInfo << "\n=====================================================================================================\n";
    gsInfo << "  SUMMARY (finest level r = " << maxRef << ", " << (wcycle ? "W-cycle" : "V-cycle") << ")\n";
    gsInfo << std::string(101, '-') << "\n";
    gsInfo << std::setw(45) << "geometry" << std::setw(10) << "MG-it"
           << std::setw(14) << "rel.res." << std::setw(14) << "H2 error"
           << std::setw(14) << "converged" << "\n";
    gsInfo << std::string(101, '-') << "\n";
    for (const Summary& s : summaries)
        gsInfo << std::setw(45) << s.name << std::setw(10) << s.iters
               << std::setw(14) << std::scientific << std::setprecision(2) << s.relres
               << std::setw(14) << std::scientific << std::setprecision(2) << s.h2
               << std::setw(14) << (s.ok ? "yes" : "NO") << "\n";
    gsInfo << std::string(101, '-') << "\n";
    gsInfo << "  Converged: " << nConverged << " / " << nTotal << " geometries.\n";
    gsInfo << "=====================================================================================================\n";

    return (nConverged == nTotal) ? 0 : 1;
}
