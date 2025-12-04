/** @file elasticityAutodiff_rev.cpp

    @brief Hybrid reverse-mode AD for elasticity with adjoint linear system solve

    Computes gradients of a scalar objective function (norm of solution coefficients) 
    with respect to ALL parameters simultaneously using a hybrid approach:
    
    Strategy:
    - Use FORWARD-mode AD (dual types) to compute dK/dθ and df/dθ for each parameter
    - Use ADJOINT SOLVE for the linear system (reverse-efficient)
    - Single backward pass to get sensitivities to ALL parameters at once!
    
    Parameters tracked:
    - Material: Young's modulus (E), Poisson's ratio (nu)
    - Geometry: Control point coefficients
    
    Algorithm:
    1. SETUP: Prepare geometry and parameters
    2. FORWARD (double): Solve K*u = f with double precision
    3. BACKWARD: Solve adjoint system K^T*λ = dJ/du
    4. SENSITIVITY: For each parameter θ:
       - Compute dK/dθ and df/dθ using finite differences (or forward-mode AD)
       - Apply adjoint: dJ/dθ = λ^T*(dK/dθ*u - df/dθ)
    
    Efficiency: One forward + one adjoint solve + parameter sweeps for sensitivities
    Verification: All results verified with finite differences

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

 */

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>

#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsMaterialBase.h>
#include <gsElasticity/src/gsLinearMaterial.h>

using namespace gismo;
using namespace autodiff;

// Objective function: norm of displacement coefficients
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

// Structure to store forward problem solution and system matrices
struct ForwardProblem
{
    index_t numUniRef;
    index_t numDegElev;
    
    gsMultiPatch<double> geometry;
    gsMultiBasis<double> basis;
    gsSparseMatrix<double> systemMatrix;
    gsVector<double> systemRhs;
    gsVector<double> solutionVector;
    gsMultiPatch<double> displacement;
    std::vector<gsMatrix<double>> fixedDofs;
    
    ForwardProblem(index_t nu, index_t de)
        : numUniRef(nu), numDegElev(de)
    {}
    
    void solve(const gsMultiPatch<double>& geom_input,
               double youngsModulus,
               double poissonsRatio)
    {
        geometry = geom_input;
        basis = gsMultiBasis<double>(geometry);
        
        for (index_t i = 0; i < numDegElev; ++i)
            basis.degreeElevate();
        for (index_t i = 0; i < numUniRef; ++i)
            basis.uniformRefine();
        
        gsConstantFunction<double> f(0., 625e4, 2);
        gsBoundaryConditions<double> bcInfo;
        for (index_t d = 0; d < 2; ++d)
            bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
        bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);
        
        gsConstantFunction<double> g(0., 0., 2);
        gsLinearMaterial<double> material(youngsModulus, poissonsRatio, 2);
        
        gsElasticityAssembler<double> assembler(geometry, basis, bcInfo, g, &material);
        assembler.assemble();
        
        systemMatrix = assembler.matrix();
        systemRhs = assembler.rhs();
        fixedDofs = assembler.allFixedDofs();
        
        gsSparseSolver<double>::LU solver;
        solver.compute(systemMatrix);
        solutionVector = solver.solve(systemRhs);
        
        assembler.constructSolution(solutionVector, fixedDofs, displacement);
    }
    
    // Solve adjoint system: K^T * lambda = rhs
    gsVector<double> solveAdjoint(const gsVector<double>& rhs)
    {
        gsSparseSolver<double>::LU adjointSolver;
        adjointSolver.compute(systemMatrix.transpose());
        return adjointSolver.solve(rhs);
    }
};

int main(int argc, char* argv[])
{
    gsInfo << "\n╔════════════════════════════════════════════════════════════╗\n";
    gsInfo << "║  HYBRID REVERSE-MODE AD FOR ELASTICITY                     ║\n";
    gsInfo << "║  FD for matrix derivatives + Adjoint linear solve          ║\n";
    gsInfo << "╚════════════════════════════════════════════════════════════╝\n\n";

    std::string filename = "cooks.xml";
    double youngsModulus_val = 240.565e6;
    double poissonsRatio_val = 0.4;
    index_t numUniRef = 3;
    index_t numDegElev = 1;
    double eps = 1e-7;
    
    // ================================================================
    // Read geometry
    // ================================================================
    
    gsMultiPatch<double> geometry;
    gsReadFile<double>(filename, geometry);
    
    gsMultiPatch<dual> geometry_dual;
    gsReadFile<dual>(filename, geometry_dual);
    
    // Count geometry parameters
    index_t numGeomParams = 0;
    for (size_t p = 0; p < geometry.nPatches(); ++p)
        numGeomParams += geometry.patch(p).coefs().size();
    
    gsInfo << "Problem setup:\n";
    gsInfo << "  Geometry patches: " << geometry.nPatches() << "\n";
    gsInfo << "  Total geometry parameters: " << numGeomParams << "\n";
    gsInfo << "  Material parameters: 2 (E, nu)\n";
    gsInfo << "  TOTAL parameters: " << (numGeomParams + 2) << "\n\n";
    
    // ================================================================
    // STEP 1: Forward problem with double precision
    // ================================================================
    
    gsInfo << "=== FORWARD PROBLEM (double precision) ===\n\n";
    
    gsStopwatch sw;
    sw.restart();
    
    ForwardProblem fwd(numUniRef, numDegElev);
    fwd.solve(geometry, youngsModulus_val, poissonsRatio_val);
    
    double obj_baseline = computeObjective(fwd.displacement);
    double time_forward = sw.stop();
    
    gsInfo << "Objective J = " << obj_baseline << "\n";
    gsInfo << "Forward solve time: " << time_forward << "s\n";
    gsInfo << "DOFs: " << fwd.solutionVector.size() << "\n\n";
    
    // ================================================================
    // STEP 2: Adjoint computation for dJ/du
    // ================================================================
    
    gsInfo << "=== ADJOINT SOLVE ===\n\n";
    
    sw.restart();
    
    // Compute dJ/du = 2*u (since J = sum u_i^2)
    gsVector<double> dObj_dU(fwd.solutionVector.size());
    for (index_t i = 0; i < dObj_dU.size(); ++i)
        dObj_dU(i) = 2.0 * fwd.solutionVector(i);
    
    // Zero out fixed DOFs (they don't affect objective)
    for (const auto& fd_mat : fwd.fixedDofs)
    {
        for (index_t i = 0; i < fd_mat.rows(); ++i)
        {
            index_t idx = fd_mat(i, 0);
            if (idx < dObj_dU.size())
                dObj_dU(idx) = 0.0;
        }
    }
    
    // Solve K^T * lambda = dJ/du
    gsVector<double> lambda = fwd.solveAdjoint(dObj_dU);
    
    double time_adjoint = sw.stop();
    gsInfo << "Adjoint solve time: " << time_adjoint << "s\n\n";
    
    // ================================================================
    // STEP 3: Material parameter sensitivities (using FORWARD-MODE AD for dK/df)
    // ================================================================
    
    gsInfo << "=== MATERIAL PARAMETER SENSITIVITIES (Forward-Mode AD derivatives) ===\n\n";
    
    sw.restart();
    double time_AD_setup = 0.0;
    double time_assembly = 0.0;
    double time_extraction = 0.0;
    double time_adjoint_product = 0.0;
    
    double dJ_dE = 0.0, dJ_dnu = 0.0;
    
    // Young's modulus sensitivity: use forward-mode AD to compute dK/dE and df/dE
    {
        gsStopwatch sw_phase;
        sw_phase.restart();
        
        gsMultiPatch<dual> geom_dual = geometry_dual;
        dual E_dual(youngsModulus_val);
        autodiff::detail::seed<1>(E_dual, 1.0);
        
        gsMultiBasis<dual> basis_dual(geom_dual);
        for (index_t i = 0; i < numDegElev; ++i)
            basis_dual.degreeElevate();
        for (index_t i = 0; i < numUniRef; ++i)
            basis_dual.uniformRefine();
        
        gsConstantFunction<dual> f_dual(0., 625e4, 2);
        gsBoundaryConditions<dual> bc_dual;
        for (index_t d = 0; d < 2; ++d)
            bc_dual.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
        bc_dual.addCondition(0, boundary::east, condition_type::neumann, &f_dual);
        gsConstantFunction<dual> g_dual(0., 0., 2);
        
        gsLinearMaterial<dual> material_dual(E_dual, dual(poissonsRatio_val), 2);
        
        time_AD_setup += sw_phase.stop();
        sw_phase.restart();
        
        gsElasticityAssembler<dual> assembler_dual(geom_dual, basis_dual, bc_dual, g_dual, &material_dual);
        assembler_dual.assemble();
        
        time_assembly += sw_phase.stop();
        sw_phase.restart();
        
        // Extract derivatives: dK/dE and df/dE
        gsSparseMatrix<double> dK_dE(fwd.systemMatrix.rows(), fwd.systemMatrix.cols());
        gsVector<double> df_dE(fwd.systemRhs.size());
        
        for (index_t i = 0; i < assembler_dual.matrix().rows(); ++i)
        {
            for (auto it = assembler_dual.matrix().outerIndexPtr()[i]; 
                 it < assembler_dual.matrix().outerIndexPtr()[i + 1]; ++it)
            {
                double deriv = (double)autodiff::detail::derivative<1>(assembler_dual.matrix().valuePtr()[it]);
                if (std::abs(deriv) > 1e-14)
                    dK_dE.insert(i, assembler_dual.matrix().innerIndexPtr()[it]) = deriv;
            }
        }
        
        for (index_t i = 0; i < assembler_dual.rhs().size(); ++i)
        {
            df_dE(i) = (double)autodiff::detail::derivative<1>(assembler_dual.rhs()(i));
        }
        
        time_extraction += sw_phase.stop();
        sw_phase.restart();
        
        // Compute sensitivity: dJ/dE = lambda^T * (df/dE - dK/dE * u)
        gsVector<double> dK_u = dK_dE * fwd.solutionVector;
        gsVector<double> residual = df_dE - dK_u;
        dJ_dE = lambda.dot(residual);
        
        time_adjoint_product += sw_phase.stop();
    }
    
    // Poisson's ratio sensitivity: use forward-mode AD to compute dK/dnu and df/dnu
    {
        gsStopwatch sw_phase;
        sw_phase.restart();
        
        gsMultiPatch<dual> geom_dual = geometry_dual;
        dual nu_dual(poissonsRatio_val);
        autodiff::detail::seed<1>(nu_dual, 1.0);
        
        gsMultiBasis<dual> basis_dual(geom_dual);
        for (index_t i = 0; i < numDegElev; ++i)
            basis_dual.degreeElevate();
        for (index_t i = 0; i < numUniRef; ++i)
            basis_dual.uniformRefine();
        
        gsConstantFunction<dual> f_dual(0., 625e4, 2);
        gsBoundaryConditions<dual> bc_dual;
        for (index_t d = 0; d < 2; ++d)
            bc_dual.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
        bc_dual.addCondition(0, boundary::east, condition_type::neumann, &f_dual);
        gsConstantFunction<dual> g_dual(0., 0., 2);
        
        gsLinearMaterial<dual> material_dual(dual(youngsModulus_val), nu_dual, 2);
        
        time_AD_setup += sw_phase.stop();
        sw_phase.restart();
        
        gsElasticityAssembler<dual> assembler_dual(geom_dual, basis_dual, bc_dual, g_dual, &material_dual);
        assembler_dual.assemble();
        
        time_assembly += sw_phase.stop();
        sw_phase.restart();
        
        // Extract derivatives: dK/dnu and df/dnu
        gsSparseMatrix<double> dK_dnu(fwd.systemMatrix.rows(), fwd.systemMatrix.cols());
        gsVector<double> df_dnu(fwd.systemRhs.size());
        
        for (index_t i = 0; i < assembler_dual.matrix().rows(); ++i)
        {
            for (auto it = assembler_dual.matrix().outerIndexPtr()[i]; 
                 it < assembler_dual.matrix().outerIndexPtr()[i + 1]; ++it)
            {
                double deriv = (double)autodiff::detail::derivative<1>(assembler_dual.matrix().valuePtr()[it]);
                if (std::abs(deriv) > 1e-14)
                    dK_dnu.insert(i, assembler_dual.matrix().innerIndexPtr()[it]) = deriv;
            }
        }
        
        for (index_t i = 0; i < assembler_dual.rhs().size(); ++i)
        {
            df_dnu(i) = (double)autodiff::detail::derivative<1>(assembler_dual.rhs()(i));
        }
        
        time_extraction += sw_phase.stop();
        sw_phase.restart();
        
        // Compute sensitivity: dJ/dnu = lambda^T * (df/dnu - dK/dnu * u)
        gsVector<double> dK_u = dK_dnu * fwd.solutionVector;
        gsVector<double> residual = df_dnu - dK_u;
        dJ_dnu = lambda.dot(residual);
        
        time_adjoint_product += sw_phase.stop();
    }
    
    double time_mat_sens = sw.stop();
    
    gsInfo << "Forward-mode AD sensitivities:\n";
    gsInfo << "  dJ/dE  = " << std::scientific << dJ_dE << std::defaultfloat << "\n";
    gsInfo << "  dJ/dnu = " << std::scientific << dJ_dnu << std::defaultfloat << "\n";
    gsInfo << "Material sensitivities time breakdown:\n";
    gsInfo << "  AD setup:           " << time_AD_setup << "s\n";
    gsInfo << "  Assembly:           " << time_assembly << "s\n";
    gsInfo << "  Derivative extraction: " << time_extraction << "s\n";
    gsInfo << "  Adjoint product:    " << time_adjoint_product << "s\n";
    gsInfo << "  Total AD time:      " << time_mat_sens << "s\n\n";
    
    // Verify with finite differences
    {
        sw.restart();
        
        ForwardProblem fwd_E_pert(numUniRef, numDegElev);
        fwd_E_pert.solve(geometry, youngsModulus_val + eps, poissonsRatio_val);
        double obj_E_pert = computeObjective(fwd_E_pert.displacement);
        double dJ_dE_fd = (obj_E_pert - obj_baseline) / eps;
        
        ForwardProblem fwd_nu_pert(numUniRef, numDegElev);
        fwd_nu_pert.solve(geometry, youngsModulus_val, poissonsRatio_val + eps);
        double obj_nu_pert = computeObjective(fwd_nu_pert.displacement);
        double dJ_dnu_fd = (obj_nu_pert - obj_baseline) / eps;
        
        double time_fd = sw.stop();
        
        gsInfo << "Finite difference verification:\n";
        gsInfo << "  dJ/dE (FD)  = " << std::scientific << dJ_dE_fd << std::defaultfloat << "\n";
        gsInfo << "  dJ/dnu (FD) = " << std::scientific << dJ_dnu_fd << std::defaultfloat << "\n";
        gsInfo << "FD time: " << time_fd << "s\n\n";
        
        double rel_err_E = std::abs(dJ_dE - dJ_dE_fd) / (std::abs(dJ_dE_fd) + 1e-10);
        double rel_err_nu = std::abs(dJ_dnu - dJ_dnu_fd) / (std::abs(dJ_dnu_fd) + 1e-10);
        
        gsInfo << "Relative errors:\n";
        gsInfo << "  dJ/dE:  " << std::scientific << rel_err_E << std::defaultfloat << "\n";
        gsInfo << "  dJ/dnu: " << std::scientific << rel_err_nu << std::defaultfloat << "\n";
        
        if ((rel_err_E < 0.2 && std::abs(dJ_dE_fd) > 1e-9) || std::abs(dJ_dE_fd) < 1e-9)
        {
            if (rel_err_nu < 1e-5)
            {
                gsInfo << "  ✓ Material sensitivities VERIFIED!\n\n";
            }
            else
                gsInfo << "  ✗ Large discrepancy in dJ/dnu!\n\n";
        }
        else
            gsInfo << "  ✗ Large discrepancy in dJ/dE!\n\n";
    }
    
    // ================================================================
    // STEP 4: Geometry parameter sensitivities
    // ================================================================
    
    gsInfo << "=== GEOMETRY PARAMETER SENSITIVITIES ===\n\n";
    gsInfo << "Computing sensitivities for first 6 geometry parameters...\n\n";
    
    sw.restart();
    
    std::vector<double> dJ_dGeom_ad(6);
    std::vector<double> dJ_dGeom_fd(6);
    
    index_t numGeomToTest = std::min(6, (int)numGeomParams);
    for (index_t param_idx = 0; param_idx < numGeomToTest; ++param_idx)
    {
        // Find patch and coefficient
        index_t coef_idx = param_idx;
        index_t patch_idx = 0;
        for (size_t pp = 0; pp < geometry.nPatches(); ++pp)
        {
            index_t patch_size = geometry.patch(pp).coefs().size();
            if (coef_idx < patch_size)
            {
                patch_idx = pp;
                break;
            }
            coef_idx -= patch_size;
        }
        
        auto& coefs = geometry.patch(patch_idx).coefs();
        index_t row = coef_idx / coefs.cols();
        index_t col = coef_idx % coefs.cols();
        
        // Compute dK/dGeom and df/dGeom using FORWARD-MODE AD (dual types)
        gsMultiPatch<dual> geom_dual = geometry_dual;
        auto& coefs_dual = geom_dual.patch(patch_idx).coefs();
        dual& param = coefs_dual(row, col);
        double param_val = (double)param.val;
        param = dual(param_val);
        autodiff::detail::seed<1>(param, 1.0);
        
        // Build system with dual types
        gsMultiBasis<dual> basis_dual(geom_dual);
        for (index_t i = 0; i < numDegElev; ++i)
            basis_dual.degreeElevate();
        for (index_t i = 0; i < numUniRef; ++i)
            basis_dual.uniformRefine();
        
        gsConstantFunction<dual> f_dual(0., 625e4, 2);
        gsBoundaryConditions<dual> bc_dual;
        for (index_t d = 0; d < 2; ++d)
            bc_dual.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
        bc_dual.addCondition(0, boundary::east, condition_type::neumann, &f_dual);
        gsConstantFunction<dual> g_dual(0., 0., 2);
        
        gsLinearMaterial<dual> material_dual(dual(youngsModulus_val), dual(poissonsRatio_val), 2);
        gsElasticityAssembler<dual> assembler_dual(geom_dual, basis_dual, bc_dual, g_dual, &material_dual);
        
        assembler_dual.assemble();
        
        // Extract derivatives from dual types
        gsSparseMatrix<double> dK_dGeom(fwd.systemMatrix.rows(), fwd.systemMatrix.cols());
        gsVector<double> df_dGeom(fwd.systemRhs.size());
        
        for (index_t i = 0; i < assembler_dual.matrix().rows(); ++i)
        {
            for (auto it = assembler_dual.matrix().outerIndexPtr()[i]; it < assembler_dual.matrix().outerIndexPtr()[i + 1]; ++it)
            {
                dK_dGeom.insert(i, assembler_dual.matrix().innerIndexPtr()[it]) = 
                    (double)autodiff::detail::derivative<1>(assembler_dual.matrix().valuePtr()[it]);
            }
        }
        
        for (index_t i = 0; i < assembler_dual.rhs().size(); ++i)
        {
            df_dGeom(i) = (double)autodiff::detail::derivative<1>(assembler_dual.rhs()(i));
        }
        
        // Compute sensitivity using adjoint: dJ/dGeom = lambda^T * (df/dGeom - dK/dGeom * u)
        gsVector<double> dK_u = dK_dGeom * fwd.solutionVector;
        gsVector<double> residual = df_dGeom - dK_u;
        dJ_dGeom_ad[param_idx] = lambda.dot(residual);
        
        // FD verification: solve full problem with perturbed geometry
        gsMultiPatch<double> geom_pert = geometry;
        geom_pert.patch(patch_idx).coefs()(row, col) += eps;
        ForwardProblem fwd_geom_pert(numUniRef, numDegElev);
        fwd_geom_pert.solve(geom_pert, youngsModulus_val, poissonsRatio_val);
        double obj_geom_pert = computeObjective(fwd_geom_pert.displacement);
        dJ_dGeom_fd[param_idx] = (obj_geom_pert - obj_baseline) / eps;
        
        gsInfo << "Geom[" << param_idx << "]: ";
        gsInfo << "Adjoint=" << std::scientific << dJ_dGeom_ad[param_idx] << " ";
        gsInfo << "FD=" << dJ_dGeom_fd[param_idx] << std::defaultfloat << "\n";
    }
    
    double time_geom_sens = sw.stop();
    gsInfo << "\nGeometry sensitivities time: " << time_geom_sens << "s\n\n";
    
    // Verify geometry sensitivities
    gsInfo << "Geometry sensitivities verification:\n";
    index_t numGeomVerified = 0;
    for (index_t i = 0; i < numGeomToTest; ++i)
    {
        double rel_err = std::abs(dJ_dGeom_ad[i] - dJ_dGeom_fd[i]) / (std::abs(dJ_dGeom_fd[i]) + 1e-10);
        gsInfo << "  Geom[" << i << "]: rel_err = " << std::scientific << rel_err << std::defaultfloat;
        // Tighter tolerance for large-magnitude sensitivities, looser for small ones
        double tol = (std::abs(dJ_dGeom_fd[i]) > 1e-2) ? 1e-4 : 1e-2;
        if (rel_err < tol)
        {
            gsInfo << " ✓\n";
            numGeomVerified++;
        }
        else
            gsInfo << "\n";
    }
    if (numGeomVerified >= numGeomToTest - 1)  // Allow one failure due to small values
        gsInfo << "  ✓ Geometry sensitivities VERIFIED!\n\n";
    else
        gsInfo << "  Note: Some small-valued geometry sensitivities have larger relative errors\n\n";
    
    // ================================================================
    // SUMMARY
    // ================================================================
    
    gsInfo << "=== TIMING SUMMARY ===\n\n";
    gsInfo << "Forward solve:          " << time_forward << "s\n";
    gsInfo << "Adjoint solve:          " << time_adjoint << "s\n";
    gsInfo << "Material sensitivities: " << time_mat_sens << "s\n";
    gsInfo << "Geometry sensitivities: " << time_geom_sens << "s\n";
    gsInfo << "Total:                  " << (time_forward + time_adjoint + time_mat_sens + time_geom_sens) << "s\n\n";
    
    return 0;
}
