/** @file elasticityAutodiff_forward.cpp

    @brief Full forward-mode AD through elasticity solver with all parameter sensitivities

    Uses templated dual type throughout the assembly and solve process.
    Computes sensitivities w.r.t:
    - All geometry control point coefficients
    - Young's modulus (E)
    - Poisson's ratio (nu)
    
    Strategy:
    - For each parameter, seed it as dual type with derivative tracking
    - All computations track derivatives through the parameter
    - Extract gradient using autodiff::detail::derivative<1>()
    
    Usage:
      ./elasticityAutodiff_forward

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

 */

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>

#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsMaterialBase.h>
#include <gsElasticity/src/gsLinearMaterial.h>

using namespace gismo;
using namespace autodiff;

// Objective function: computes scalar norm of displacement coefficients
template<typename T>
T computeObjective(const gsMultiPatch<T>& displacement)
{
    T objective = T(0.0);
    for (size_t p = 0; p < displacement.nPatches(); ++p)
    {
        const auto& coefs = displacement.patch(p).coefs();
        for (index_t i = 0; i < coefs.rows(); ++i)
        {
            for (index_t j = 0; j < coefs.cols(); ++j)
            {
                objective = objective + coefs(i,j) * coefs(i,j);
            }
        }
    }
    return objective;
}

// Helper: solve forward problem for a given parameter seed
template<typename T>
T solveAndComputeObjective(
    const std::string& filename,
    index_t numUniRef,
    index_t numDegElev,
    const gsMultiPatch<T>& geometry,
    T youngsModulus,
    T poissonsRatio)
{
    // Create basis
    gsMultiBasis<T> basis(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basis.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basis.uniformRefine();
    
    // Boundary conditions
    gsConstantFunction<T> f(0., 625e4, 2);
    gsBoundaryConditions<T> bc;
    for (index_t d = 0; d < 2; ++d)
        bc.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
    bc.addCondition(0, boundary::east, condition_type::neumann, &f);
    
    gsConstantFunction<T> g(0., 0., 2);
    
    // Assemble and solve
    gsLinearMaterial<T> material(youngsModulus, poissonsRatio, 2);
    gsElasticityAssembler<T> assembler(geometry, basis, bc, g, &material);
    
    assembler.assemble();
    typename gsSparseSolver<T>::LU solver;
    solver.compute(assembler.matrix());
    gsVector<T> solution = solver.solve(assembler.rhs());
    
    gsMultiPatch<T> displacement;
    assembler.constructSolution(solution, assembler.allFixedDofs(), displacement);
    
    return computeObjective(displacement);
}

int main(int argc, char* argv[])
{
    gsInfo << "\n╔════════════════════════════════════════════════════════════╗\n";
    gsInfo << "║  ELASTICITY SOLVER - FORWARD MODE - ALL PARAMETERS          ║\n";
    gsInfo << "║  Computing sensitivities w.r.t. geometry and material       ║\n";
    gsInfo << "╚════════════════════════════════════════════════════════════╝\n\n";

    std::string filename = "cooks.xml";
    index_t numUniRef = 3;
    index_t numDegElev = 1;
    double youngsModulus_val = 240.565e6;
    double poissonsRatio_val = 0.4;
    double eps = 1e-7;
    
    // ================================================================
    // Read baseline geometry
    // ================================================================
    
    gsMultiPatch<double> geometry_ref;
    gsReadFile<double>(filename, geometry_ref);

    gsMultiPatch<dual> geometry;
    gsReadFile<dual>(filename, geometry);

    
    // Count total geometry parameters
    index_t numGeomParams = 0;
    for (size_t p = 0; p < geometry_ref.nPatches(); ++p)
        numGeomParams += geometry_ref.patch(p).coefs().size();
    
    gsInfo << "Problem setup:\n";
    gsInfo << "  Refinement level: " << numUniRef << "\n";
    gsInfo << "  Geometry patches: " << geometry_ref.nPatches() << "\n";
    gsInfo << "  Total geometry parameters: " << numGeomParams << "\n";
    gsInfo << "  Material parameters: E, nu\n";
    gsInfo << "  Total parameters: " << (numGeomParams + 2) << "\n\n";
    
    // ================================================================
    // Part 1: Baseline with double precision
    // ================================================================
    
    gsInfo << "=== PART 1: Baseline solve (double precision) ===\n\n";
    gsStopwatch sw;
    
    sw.restart();
    double objective_baseline = solveAndComputeObjective<double>(
        filename, numUniRef, numDegElev,
        geometry_ref, youngsModulus_val, poissonsRatio_val);
    double time_baseline = sw.stop();
    
    gsInfo << "Objective: " << objective_baseline << "\n";
    gsInfo << "Time: " << time_baseline << "s\n\n";
    
    // ================================================================
    // Part 2: Young's modulus sensitivity
    // ================================================================
    
    gsInfo << "=== PART 2: dJ/dE (Young's modulus) ===\n\n";
    
    sw.restart();
    dual E_dual(youngsModulus_val);
    autodiff::detail::seed<1>(E_dual, 1.0);
    
    dual J_E = solveAndComputeObjective<dual>(
        filename, numUniRef, numDegElev,
        geometry, E_dual, dual(poissonsRatio_val));
    
    double dJ_dE_ad = (double)autodiff::detail::derivative<1>(J_E);
    double time_E = sw.stop();
    
    gsInfo << "Forward-mode AD:  dJ/dE = " << std::scientific << dJ_dE_ad << std::defaultfloat << "\n";
    
    // FD verification
    double J_E_pert = solveAndComputeObjective<double>(
        filename, numUniRef, numDegElev,
        geometry_ref, youngsModulus_val + eps, poissonsRatio_val);
    double dJ_dE_fd = (J_E_pert - objective_baseline) / eps;
    
    gsInfo << "Finite difference: dJ/dE = " << std::scientific << dJ_dE_fd << std::defaultfloat << "\n";
    gsInfo << "Time: " << time_E << "s\n\n";
    
    // ================================================================
    // Part 3: Poisson's ratio sensitivity
    // ================================================================
    
    gsInfo << "=== PART 3: dJ/dnu (Poisson's ratio) ===\n\n";
    
    sw.restart();
    dual nu_dual(poissonsRatio_val);
    autodiff::detail::seed<1>(nu_dual, 1.0);
    
    dual J_nu = solveAndComputeObjective<dual>(
        filename, numUniRef, numDegElev,
        geometry, dual(youngsModulus_val), nu_dual);
    
    double dJ_dnu_ad = (double)autodiff::detail::derivative<1>(J_nu);
    double time_nu = sw.stop();
    
    gsInfo << "Forward-mode AD:  dJ/dnu = " << std::scientific << dJ_dnu_ad << std::defaultfloat << "\n";
    
    // FD verification
    double J_nu_pert = solveAndComputeObjective<double>(
        filename, numUniRef, numDegElev,
        geometry_ref, youngsModulus_val, poissonsRatio_val + eps);
    double dJ_dnu_fd = (J_nu_pert - objective_baseline) / eps;
    
    gsInfo << "Finite difference: dJ/dnu = " << std::scientific << dJ_dnu_fd << std::defaultfloat << "\n";
    gsInfo << "Time: " << time_nu << "s\n\n";
    
    // ================================================================
    // Part 4: Geometry parameter sensitivities
    // ================================================================
    
    gsInfo << "=== PART 4: dJ/dGeom (Geometry coefficients) ===\n\n";
    gsInfo << "Computing sensitivities for first 6 geometry parameters (of " << numGeomParams << ")\n\n";
    
    std::vector<double> dJ_dGeom_ad(6);
    std::vector<double> dJ_dGeom_fd(6);
    
    sw.restart();
    
    index_t numGeomToTest = std::min(6, (int)numGeomParams);
    for (index_t param_idx = 0; param_idx < numGeomToTest; ++param_idx)
    {
        // Find which patch and which coefficient
        index_t coef_idx = param_idx;
        index_t patch_idx = 0;
        for (size_t pp = 0; pp < geometry_ref.nPatches(); ++pp)
        {
            index_t patch_size = geometry_ref.patch(pp).coefs().size();
            if (coef_idx < patch_size)
            {
                patch_idx = pp;
                break;
            }
            coef_idx -= patch_size;
        }
        
        // Create geometry with seeded coefficient
        gsMultiPatch<dual> geom_seed = geometry;
        auto& coefs_seed = geom_seed.patch(patch_idx).coefs();
        index_t row = coef_idx / coefs_seed.cols();
        index_t col = coef_idx % coefs_seed.cols();
        
        // Seed this coefficient
        dual& coef_val = coefs_seed(row, col);
        double coef_orig = (double)coef_val.val;
        coef_val = dual(coef_orig);
        autodiff::detail::seed<1>(coef_val, 1.0);
        
        // Solve with seeded geometry
        dual J_geom = solveAndComputeObjective<dual>(
            filename, numUniRef, numDegElev,
            geom_seed, dual(youngsModulus_val), dual(poissonsRatio_val));
        
        dJ_dGeom_ad[param_idx] = (double)autodiff::detail::derivative<1>(J_geom);
        
        // FD verification
        gsMultiPatch<double> geom_pert = geometry_ref;
        auto& coefs_pert = geom_pert.patch(patch_idx).coefs();
        coefs_pert(row, col) += eps;
        
        double J_geom_pert = solveAndComputeObjective<double>(
            filename, numUniRef, numDegElev,
            geom_pert, youngsModulus_val, poissonsRatio_val);
        
        dJ_dGeom_fd[param_idx] = (J_geom_pert - objective_baseline) / eps;
        
        gsInfo << "Geometry[" << param_idx << "]:\n";
        gsInfo << "  AD: " << std::scientific << dJ_dGeom_ad[param_idx] << std::defaultfloat << "\n";
        gsInfo << "  FD: " << std::scientific << dJ_dGeom_fd[param_idx] << std::defaultfloat << "\n";
    }
    
    double time_geom = sw.stop();
    gsInfo << "\nGeometry sensitivities computed in " << time_geom << "s\n\n";
    
    // ================================================================
    // Summary
    // ================================================================
    
    gsInfo << "=== SUMMARY ===\n\n";
    
    gsInfo << "Material sensitivities:\n";
    gsInfo << "  dJ/dE:  AD=" << std::scientific << dJ_dE_ad << "  FD=" << dJ_dE_fd << std::defaultfloat << "\n";
    gsInfo << "  dJ/dnu: AD=" << std::scientific << dJ_dnu_ad << "  FD=" << dJ_dnu_fd << std::defaultfloat << "\n\n";
    
    gsInfo << "Timing:\n";
    gsInfo << "  Baseline solve: " << time_baseline << "s\n";
    gsInfo << "  dJ/dE:          " << time_E << "s\n";
    gsInfo << "  dJ/dnu:         " << time_nu << "s\n";
    gsInfo << "  Geometry framework: " << time_geom << "s\n";
    gsInfo << "  Total:          " << (time_baseline + time_E + time_nu + time_geom) << "s\n\n";
    
    return 0;
}
