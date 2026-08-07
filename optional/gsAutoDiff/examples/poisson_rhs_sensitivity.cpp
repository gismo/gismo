/** @file poisson_rhs_sensitivity.cpp

    @brief Adjoint sensitivity of a Poisson solution w.r.t. a RHS parameter.

    Demonstrates how to compute dL/dtheta where:
      - L = ||u|| is a scalar loss on the Poisson solution u
      - theta is a scalar parameter in the RHS  f(x) = theta
      - K u = theta * f0  is the discretized Poisson equation

    The adjoint method avoids differentiating through the linear solver:
      1. Forward solve:  K u = f(theta)                       (real_t)
      2. Adjoint solve:  K^T lambda = dL/du                   (real_t)
      3. AD assembly:    assemble K, f(theta) with dual_t     (forward AD)
      4. Extract:        dL/dtheta = lambda^T * (df/dtheta)

    Since theta is a single scalar parameter, forward-mode AD (dual_t) is
    the natural choice: one forward pass gives the derivative w.r.t. theta.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>
#include <gsAutoDiff/gsAutoDiffEigen.h>
#include <gsAssembler/gsPoissonAssembler.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Adjoint sensitivity of Poisson solution w.r.t. RHS parameter theta.");
    index_t numRefine  = 2;
    index_t numElevate = 1;
    real_t  theta_val  = 2.0;
    real_t  fd_eps     = 1e-6;

    cmd.addInt("r", "refine",  "Number of uniform refinements", numRefine);
    cmd.addInt("e", "elevate", "Number of degree elevations",   numElevate);
    cmd.addReal("t", "theta",  "RHS parameter theta",           theta_val);
    cmd.addReal("", "fd-eps",  "Finite difference step size",   fd_eps);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

#ifdef GISMO_AUTODIFF_FORWARD
    // ================================================================
    //  Step 1: Baseline Poisson solve with real_t
    // ================================================================
    gsInfo << "=== Step 1: Baseline solve (real_t, theta = " << theta_val << ") ===\n";

    gsMultiPatch<real_t> patches;
    patches.addPatch( gsNurbsCreator<real_t>::BSplineSquare(1.0) );
    for (index_t i = 0; i < numElevate; ++i) patches.degreeElevate();
    for (index_t i = 0; i < numRefine;  ++i) patches.uniformRefine();

    gsMultiBasis<real_t> basis(patches);

    gsConstantFunction<real_t> zero(0.0, 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(boundary::west,  condition_type::dirichlet, &zero);
    bc.addCondition(boundary::east,  condition_type::dirichlet, &zero);
    bc.addCondition(boundary::north, condition_type::dirichlet, &zero);
    bc.addCondition(boundary::south, condition_type::dirichlet, &zero);

    gsConstantFunction<real_t> rhs(theta_val, 2);

    gsPoissonAssembler<real_t> assembler(patches, basis, bc, rhs);
    assembler.assemble();

    const gsSparseMatrix<real_t> & K = assembler.matrix();
    const gsMatrix<real_t>       & f = assembler.rhs();

    gsSparseSolver<real_t>::uPtr solver = gsSparseSolver<real_t>::get("SimplicialLDLT");
    solver->compute(K);
    gsMatrix<real_t> u_real = solver->solve(f);

    real_t loss = u_real.norm();
    gsInfo << "  DOFs:          " << u_real.rows() << "\n";
    gsInfo << "  L = ||u|| =    " << loss << "\n";

    // ================================================================
    //  Step 2: Adjoint solve  K^T lambda = dL/du
    // ================================================================
    gsInfo << "\n=== Step 2: Adjoint solve ===\n";

    gsMatrix<real_t> dL_du = u_real / loss;

    gsSparseSolver<real_t>::uPtr adj_solver = gsSparseSolver<real_t>::get("SimplicialLDLT");
    adj_solver->compute(K.transpose());
    gsMatrix<real_t> lambda = adj_solver->solve(dL_du);
    gsInfo << "  ||lambda|| =   " << lambda.norm() << "\n";

    // ================================================================
    //  Step 3: AD assembly with dual_t (forward AD)
    // ================================================================
    gsInfo << "\n=== Step 3: Assemble with dual_t (forward AD, seed theta) ===\n";

    dual_t theta_d(theta_val);
    theta_d.grad = 1.0;

    gsMultiPatch<dual_t> patches_d;
    patches_d.addPatch( *gsNurbsCreator<dual_t>::BSplineSquare(1.0) );
    for (index_t i = 0; i < numElevate; ++i) patches_d.degreeElevate();
    for (index_t i = 0; i < numRefine;  ++i) patches_d.uniformRefine();

    gsMultiBasis<dual_t> basis_d(patches_d);

    gsConstantFunction<dual_t> zero_d(dual_t(0.0), 2);
    gsBoundaryConditions<dual_t> bc_d;
    bc_d.addCondition(boundary::west,  condition_type::dirichlet, &zero_d);
    bc_d.addCondition(boundary::east,  condition_type::dirichlet, &zero_d);
    bc_d.addCondition(boundary::north, condition_type::dirichlet, &zero_d);
    bc_d.addCondition(boundary::south, condition_type::dirichlet, &zero_d);

    gsConstantFunction<dual_t> rhs_d(theta_d, 2);

    gsStopwatch timer;
    gsPoissonAssembler<dual_t> assembler_d(patches_d, basis_d, bc_d, rhs_d);
    assembler_d.assemble();
    real_t asm_time = timer.stop();
    gsInfo << "  Assembly time:  " << asm_time << " s\n";

    // ================================================================
    //  Step 4: Extract dL/dtheta via adjoint formula
    // ================================================================
    gsInfo << "\n=== Step 4: Compute dL/dtheta via adjoint + forward AD ===\n";

    const gsSparseMatrix<dual_t> & K_d = assembler_d.matrix();
    const gsMatrix<dual_t>       & f_d = assembler_d.rhs();

    real_t dL_dtheta = 0.0;
    for (index_t i = 0; i < f_d.rows(); ++i)
        dL_dtheta += lambda(i,0) * f_d(i,0).grad;

    for (int col = 0; col < K_d.outerSize(); ++col)
        for (typename gsSparseMatrix<dual_t>::InnerIterator it(K_d, col); it; ++it)
            dL_dtheta -= lambda(it.row(),0) * it.value().grad * u_real(it.col(),0);

    gsInfo << "  dL/dtheta (AD) = " << dL_dtheta << "\n";

    // ================================================================
    //  Step 5: Finite Difference Verification
    // ================================================================
    gsInfo << "\n=== Step 5: Finite difference verification (h = " << fd_eps << ") ===\n";

    gsConstantFunction<real_t> rhs_pert(theta_val + fd_eps, 2);
    gsPoissonAssembler<real_t> asm_fd(patches, basis, bc, rhs_pert);
    asm_fd.assemble();
    gsSparseSolver<real_t>::uPtr sol_fd = gsSparseSolver<real_t>::get("SimplicialLDLT");
    sol_fd->compute(asm_fd.matrix());
    gsMatrix<real_t> u_fd = sol_fd->solve(asm_fd.rhs());
    real_t loss_fd = u_fd.norm();
    real_t dL_dtheta_fd = (loss_fd - loss) / fd_eps;

    real_t rel_err = (dL_dtheta_fd != 0.0)
                   ? math::abs(dL_dtheta - dL_dtheta_fd) / math::abs(dL_dtheta_fd)
                   : math::abs(dL_dtheta);

    gsInfo << "  dL/dtheta (FD) = " << dL_dtheta_fd << "\n";
    gsInfo << "  Rel. error     = " << rel_err << "\n";

    // ================================================================
    //  Step 6: Analytical verification
    // ================================================================
    gsInfo << "\n=== Step 6: Analytical check ===\n";
    gsInfo << "  Since K does not depend on theta and f = theta * f0:\n";
    gsInfo << "    u = theta * K^{-1} f0,  so  du/dtheta = u / theta\n";
    gsInfo << "    dL/dtheta = (u / ||u||)^T * (u / theta) = ||u|| / theta\n";
    real_t dL_dtheta_exact = loss / theta_val;
    gsInfo << "  dL/dtheta (exact) = " << dL_dtheta_exact << "\n";
    gsInfo << "  Rel. error (AD vs exact)  = "
           << math::abs(dL_dtheta - dL_dtheta_exact) / math::abs(dL_dtheta_exact) << "\n";
    gsInfo << "  Rel. error (FD vs exact)  = "
           << math::abs(dL_dtheta_fd - dL_dtheta_exact) / math::abs(dL_dtheta_exact) << "\n";

    gsInfo << "\n=== Done ===\n";
#else
    gsInfo << "This example requires GISMO_AUTODIFF_FORWARD.\n";
#endif
    return 0;
}
