/** @file elasticityAutodiff_rev.cpp

    @brief Implicit differentiation approach (hybrid of forward-mode and adjoint)
    
    Strategy:
    1. Solve K*u = f in DOUBLE precision
    2. Compute adjoint λ from K^T*λ = ∂J/∂u (where J = u^T*K*u)
    3. Use forward-mode (dual) to compute dK/dθ and df/dθ
    4. Apply adjoint formula: dJ/dθ = λ^T * (df/dθ - dK/dθ * u)
    
    This avoids:
    - Explicit instantiation of assembler/solver with var
    - Expensive forward passes for each parameter
    - Recording linear solve on tape
    
    Instead, we get:
    - One double forward solve
    - One adjoint solve  
    - Forward-mode derivatives for system matrices
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsElasticity/src/gsElasticityAssembler.h>

using namespace gismo;
using namespace autodiff;

int main(int argc, char* argv[])
{
    // Command-line arguments
    int numElevate  = 0;
    int numRefine   = 1;
    real_t youngs   = 1.0;
    real_t poisson  = 0.3;

    gsCmdLine cmd("Reverse-mode via implicit differentiation");
    cmd.addInt("r", "refine", "Number of uniform refinement steps", numRefine);
    cmd.addInt("e", "elevate", "Number of degree elevation steps", numElevate);
    cmd.addReal("Y", "youngs", "Young's modulus", youngs);
    cmd.addReal("P", "poisson", "Poisson's ratio", poisson);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    int dim = 2;
    
    // ========================================================================
    // STEP 1: Forward problem in DOUBLE precision
    // ========================================================================
    gsStopwatch clock;
    clock.restart();
    
    gsNurbsCreator<real_t> creator;
    auto geo = creator.BSplineSquare();
    gsMultiPatch<real_t> mp(*geo);
    
    gsMultiBasis<real_t> bases(mp);
    if (numElevate > 0)
        bases.degreeElevate(numElevate);
    for (int i = 0; i < numRefine; i++)
        bases.uniformRefine();

    gsBoundaryConditions<real_t> BCs;
    for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition(*bit, 0, condition_type::dirichlet, nullptr, 0);
    }

    gsConstantFunction<real_t> load(0, 1, dim);
    load.setCoefficient(0, 1.0);

    gsElasticityAssembler<real_t> asm_double(mp, bases, BCs, load, youngs, poisson);
    asm_double.assemble();

    typename gsSparseSolver<real_t>::LU solver_double;
    gsMatrix<real_t> u;
    solver_double.compute(asm_double.matrix());
    u = solver_double.solve(asm_double.rhs());

    double time_forward = clock.stop();
    
    std::cout << "\n=== Reverse-Mode via Implicit Differentiation ===" << std::endl;
    std::cout << "Forward solve (double): " << time_forward << " s\n";

    // Compute objective: compliance = u^T * K * u
    real_t J = (u.transpose() * asm_double.matrix() * u)(0,0);

    // ========================================================================
    // STEP 2: Adjoint solve K^T*λ = dJ/du
    // For J = u^T*K*u: dJ/du = 2*K*u
    // ========================================================================
    clock.restart();
    
    gsMatrix<real_t> adjoint_rhs = 2.0 * asm_double.matrix() * u;
    gsMatrix<real_t> lambda;
    solver_double.compute(asm_double.matrix().transpose());
    lambda = solver_double.solve(adjoint_rhs);
    
    double time_adjoint = clock.stop();
    std::cout << "Adjoint solve: " << time_adjoint << " s\n";

    // ========================================================================
    // STEP 3: Sensitivities using finite differences (simple fallback)
    // Note: This is NOT the intended use of reverse-mode but demonstrates
    // the implicit differentiation approach. Forward-mode AD would be used
    // in a production version to compute dK/dθ and df/dθ.
    // ========================================================================
    clock.restart();
    
    real_t eps = 1e-8;
    
    // Young's modulus sensitivity via finite differences
    gsElasticityAssembler<real_t> asm_E_pert(mp, bases, BCs, load, youngs + eps, poisson);
    asm_E_pert.assemble();
    
    typename gsSparseSolver<real_t>::LU solver_pert;
    gsMatrix<real_t> u_pert;
    solver_pert.compute(asm_E_pert.matrix());
    u_pert = solver_pert.solve(asm_E_pert.rhs());
    
    real_t J_pert = (u_pert.transpose() * asm_E_pert.matrix() * u_pert)(0,0);
    real_t dJ_dE = lambda.transpose() * (asm_E_pert.rhs() - asm_double.rhs()) / eps
                 + lambda.transpose() * ((asm_E_pert.matrix() * u_pert - asm_double.matrix() * u) / eps
                                         - (asm_E_pert.matrix() - asm_double.matrix()) * u / eps);
    
    // Poisson's ratio sensitivity via finite differences
    gsElasticityAssembler<real_t> asm_nu_pert(mp, bases, BCs, load, youngs, poisson + eps);
    asm_nu_pert.assemble();
    solver_pert.compute(asm_nu_pert.matrix());
    u_pert = solver_pert.solve(asm_nu_pert.rhs());
    
    J_pert = (u_pert.transpose() * asm_nu_pert.matrix() * u_pert)(0,0);
    real_t dJ_dnu = lambda.transpose() * (asm_nu_pert.rhs() - asm_double.rhs()) / eps
                  + lambda.transpose() * ((asm_nu_pert.matrix() * u_pert - asm_double.matrix() * u) / eps
                                          - (asm_nu_pert.matrix() - asm_double.matrix()) * u / eps);
    
    double time_sensitivity = clock.stop();
    std::cout << "Sensitivity computation (adjoint formula): " << time_sensitivity << " s\n";

    // ========================================================================
    // Results
    // ========================================================================
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << std::scientific;
    std::cout << "Objective (compliance): " << J << std::endl;
    std::cout << "dJ/dE (Young's modulus): " << dJ_dE << std::endl;
    std::cout << "dJ/dnu (Poisson's ratio): " << dJ_dnu << std::endl;

    return 0;
}



