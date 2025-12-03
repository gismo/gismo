/** @file elasticityAutodiff_rev.cpp

    @brief Example demonstrating reverse-mode automatic differentiation through elasticity solver

    Computes gradients of a scalar objective function (norm of solution coefficients) 
    with respect to material parameters (Poisson's ratio and Young's modulus) using 
    reverse-mode AD with the var type.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>
// Include reverse-mode AD with full GISMO matrix support
#include <gsAutoDiff/gsAutoDiffVar.h>

#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsMaterialBase.h>
#include <gsElasticity/src/gsLinearMaterial.h>

using namespace gismo;
using namespace autodiff;

// Objective function: computes scalar norm of displacement coefficients
// Generic to work with any scalar type (double, var, dual, etc.)
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
    return objective;  // Returns squared norm
}

// Solve linear elasticity problem with given material parameters
// Uses double precision for the solve since var types have Eigen/autodiff compatibility issues
gsMultiPatch<double> solveElasticity(const std::string& filename, 
                                      double youngsModulus, 
                                      double poissonsRatio,
                                      index_t numUniRef,
                                      index_t numDegElev)
{
    // Scanning geometry
    gsMultiPatch<double> geometry;
    gsReadFile<double>(filename, geometry);
    
    // Creating bases
    gsMultiBasis<double> basisDisplacement(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    // Neumann BC: vertical load
    gsConstantFunction<double> f(0., 625e4, 2);

    // Boundary conditions: fixed on left edge (west), load on right (east)
    gsBoundaryConditions<double> bcInfo;
    for (index_t d = 0; d < 2; ++d)
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
    bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);

    // Source function (RHS)
    gsConstantFunction<double> g(0., 0., 2);

    // Create linear material (var-compatible)
    gsLinearMaterial<double> materialMat(youngsModulus, poissonsRatio, 2);

    // Creating assembler
    gsElasticityAssembler<double> assembler(geometry, basisDisplacement, bcInfo, g, &materialMat);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // Setting Newton's method
    gsIterative<double> solver(assembler);
    solver.options().setInt("Verbosity", solver_verbosity::all);
    solver.options().setInt("Solver", linear_solver::LDLT);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    solver.solve();
    gsInfo << "Solved the system in " << clock.stop() << "s.\n";

    // Construct solution as displacement field
    gsMultiPatch<double> displacement;
    assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);

    return displacement;
}

int main(int argc, char* argv[])
{
    gsInfo << "=== Elasticity Solver with Reverse-Mode AD ===\n";
    gsInfo << "Computes gradients of displacement norm w.r.t. material parameters\n";
    gsInfo << "using reverse-mode AD (autodiff::var type)\n\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = "cooks.xml";
    double youngsModulus_val = 240.565e6;
    double poissonsRatio_val = 0.4;
    index_t numUniRef = 4;
    index_t numDegElev = 1;

    //=====================================//
       // Step 1: Solve with double precision //
    //=====================================//

    gsInfo << "\n=== STEP 1: Solve Elasticity Problem ===\n";
    gsInfo << "Parameters: E = " << youngsModulus_val << ", nu = " << poissonsRatio_val << "\n\n";

    gsInfo << "Solving with double precision...\n";
    gsMultiPatch<double> displacement = solveElasticity(filename, youngsModulus_val, poissonsRatio_val,
                                                         numUniRef, numDegElev);

    //=====================================//
    // Step 2: Compute objective with var //
    //=====================================//

    gsInfo << "\n=== STEP 2: Compute Objective (Reverse-Mode AD) ===\n";

    // Create var versions of material parameters for differentiation
    var youngsModulus(youngsModulus_val);
    var poissonsRatio(poissonsRatio_val);

    // Compute objective from double-precision solution
    gsInfo << "Computing objective function (norm of displacement coefficients)...\n";
    var objective = computeObjective(displacement);
    gsInfo << "Objective value: " << static_cast<double>(objective) << "\n";

    //=====================================//
    // Step 3: Compute gradients via reverse-mode AD //
    //=====================================//

    gsInfo << "\n=== STEP 3: Reverse-Mode AD Gradient Computation ===\n";
    gsInfo << "Computing gradients in ONE reverse pass...\n";

    auto [dObj_dE, dObj_dnu] = autodiff::derivatives(objective, 
                                                      autodiff::reverse::detail::wrt(youngsModulus, poissonsRatio));

    gsInfo << "\n=== RESULTS ===\n";
    gsInfo << "dObjective/dYoungsModulus = " << dObj_dE << "\n";
    gsInfo << "dObjective/dPoissonsRatio = " << dObj_dnu << "\n\n";

    //=====================================//
    // Step 4: Verification with FD      //
    //=====================================//

    gsInfo << "=== STEP 4: Verification (Finite Differences) ===\n";
    
    double eps = 1e-7;
    double youngsModulus_pert = youngsModulus_val + eps;
    double poissonsRatio_pert = poissonsRatio_val + eps;

    // Nominal objective (already computed)
    double obj_nominal = static_cast<double>(objective);

    // Perturb Young's modulus
    gsInfo << "\nPerturbation in Young's modulus (eps = " << eps << "):\n";
    gsMultiPatch<double> disp_E_pert = solveElasticity(filename, youngsModulus_pert, poissonsRatio_val,
                                                        numUniRef, numDegElev);
    double obj_E_pert = static_cast<double>(computeObjective(disp_E_pert));

    double fd_E = (obj_E_pert - obj_nominal) / eps;
    gsInfo << "Finite difference (E): " << fd_E << "\n";
    gsInfo << "AD gradient (E):       " << dObj_dE << "\n";
    gsInfo << "Relative error (E):    " << std::abs(fd_E - dObj_dE) / std::abs(dObj_dE) * 100.0 << "%\n";

    // Perturb Poisson's ratio
    gsInfo << "\nPerturbation in Poisson's ratio (eps = " << eps << "):\n";
    gsMultiPatch<double> disp_nu_pert = solveElasticity(filename, youngsModulus_val, poissonsRatio_pert,
                                                         numUniRef, numDegElev);
    double obj_nu_pert = static_cast<double>(computeObjective(disp_nu_pert));

    double fd_nu = (obj_nu_pert - obj_nominal) / eps;
    gsInfo << "Finite difference (nu): " << fd_nu << "\n";
    gsInfo << "AD gradient (nu):       " << dObj_dnu << "\n";
    gsInfo << "Relative error (nu):    " << std::abs(fd_nu - dObj_dnu) / std::abs(dObj_dnu) * 100.0 << "%\n";

    gsInfo << "\n=== COMPUTATION COMPLETE ===\n";
    gsInfo << "\nKey feature of reverse-mode AD:\n";
    gsInfo << "  - Both gradients computed in ONE reverse pass through the objective\n";
    gsInfo << "  - Forward-mode would require " << 2 << " forward passes (one per parameter)\n";
    gsInfo << "  - Advantage grows with number of parameters\n\n";

    return 0;
}
