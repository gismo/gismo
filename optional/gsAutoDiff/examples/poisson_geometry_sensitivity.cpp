/** @file poisson_geometry_sensitivity.cpp

    @brief Adjoint sensitivity of a Poisson solution w.r.t. geometry control points.

    Demonstrates how to compute dL/dp where:
      - L = ||u|| is a scalar loss on the Poisson solution u
      - p are geometry control points
      - K(p) u = f  is the discretized Poisson equation

    The adjoint method avoids differentiating through the linear solver:
      1. Forward solve:  K u = f              (real_t)
      2. Adjoint solve:  K^T lambda = dL/du   (real_t)
      3. AD assembly:    assemble K(p), f(p)   (var_t)
      4. Backpropagate:  R = lambda^T (f - K u), then gradient(R, p)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>
#include <gsAutoDiff/gsAutoDiffEigen.h>
#include <gsAssembler/gsPoissonAssembler.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Adjoint sensitivity of Poisson solution w.r.t. geometry control points.");
    index_t numRefine  = 1;
    index_t numElevate = 1;
    real_t  fd_eps     = 1e-6;

    cmd.addInt("r", "refine",  "Number of uniform refinements", numRefine);
    cmd.addInt("e", "elevate", "Number of degree elevations",   numElevate);
    cmd.addReal("", "fd-eps",  "Finite difference step size",   fd_eps);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

#ifdef GISMO_AUTODIFF_BACKWARD
    // ================================================================
    //  Step 1: Baseline Poisson solve with real_t
    // ================================================================
    gsInfo << "=== Step 1: Baseline solve (real_t) ===\n";

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

    gsConstantFunction<real_t> rhs(1.0, 2);

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
    //  Step 3: AD assembly with var_t
    // ================================================================
    gsInfo << "\n=== Step 3: Assemble with var_t (reverse AD) ===\n";

    gsMultiPatch<var_t> patches_v;
    patches_v.addPatch( *gsNurbsCreator<var_t>::BSplineSquare(1.0) );
    for (index_t i = 0; i < numElevate; ++i) patches_v.degreeElevate();
    for (index_t i = 0; i < numRefine;  ++i) patches_v.uniformRefine();

    gsMatrix<var_t> & coefs_v = patches_v.patch(0).coefs();
    const index_t nCoefs = coefs_v.rows() * coefs_v.cols();

    gsVector<var_t> params(nCoefs);
    for (index_t k = 0; k < nCoefs; ++k)
    {
        params(k) = coefs_v.data()[k];
        coefs_v.data()[k] = params(k);
    }

    gsMultiBasis<var_t> basis_v(patches_v);
    gsConstantFunction<var_t> zero_v(var_t(0.0), 2);
    gsBoundaryConditions<var_t> bc_v;
    bc_v.addCondition(boundary::west,  condition_type::dirichlet, &zero_v);
    bc_v.addCondition(boundary::east,  condition_type::dirichlet, &zero_v);
    bc_v.addCondition(boundary::north, condition_type::dirichlet, &zero_v);
    bc_v.addCondition(boundary::south, condition_type::dirichlet, &zero_v);

    gsConstantFunction<var_t> rhs_v(var_t(1.0), 2);

    gsStopwatch timer;
    gsPoissonAssembler<var_t> assembler_v(patches_v, basis_v, bc_v, rhs_v);
    assembler_v.assemble();
    real_t asm_time = timer.stop();
    gsInfo << "  Assembly time:  " << asm_time << " s\n";

    // ================================================================
    //  Step 4: Backpropagate  R = lambda^T (f_v - K_v u)
    // ================================================================
    gsInfo << "\n=== Step 4: Backpropagate via reverse AD ===\n";

    const gsSparseMatrix<var_t> & K_v = assembler_v.matrix();
    const gsMatrix<var_t>       & f_v = assembler_v.rhs();

    var_t R(0.0);
    for (index_t i = 0; i < f_v.rows(); ++i)
        R += var_t(lambda(i,0)) * f_v(i,0);

    for (int col = 0; col < K_v.outerSize(); ++col)
        for (typename gsSparseMatrix<var_t>::InnerIterator it(K_v, col); it; ++it)
            R -= var_t(lambda(it.row(),0)) * it.value() * var_t(u_real(it.col(),0));

    timer.restart();
    gsVector<real_t> dL_dp = autodiff::reverse::detail::gradient(R, params);
    real_t rev_time = timer.stop();
    gsInfo << "  Backprop time:  " << rev_time << " s\n";

    gsInfo << "\n  dL/dp (first " << math::min(index_t(10), nCoefs) << " of " << nCoefs << "):\n  ";
    gsInfo << dL_dp.head(math::min(index_t(10), nCoefs)).transpose() << "\n";
    gsInfo << "  ||dL/dp|| =    " << dL_dp.norm() << "\n";

    // ================================================================
    //  Step 5: Finite Difference Verification
    // ================================================================
    gsInfo << "\n=== Step 5: Finite difference verification (h = " << fd_eps << ") ===\n";

    const index_t nCheck = math::min(index_t(6), nCoefs);
    gsVector<real_t> dL_dp_fd(nCheck);
    for (index_t k = 0; k < nCheck; ++k)
    {
        patches.patch(0).coefs().data()[k] += fd_eps;
        gsMultiBasis<real_t> basis_fd(patches);
        gsPoissonAssembler<real_t> asm_fd(patches, basis_fd, bc, rhs);
        asm_fd.assemble();
        gsSparseSolver<real_t>::uPtr sol_fd = gsSparseSolver<real_t>::get("SimplicialLDLT");
        sol_fd->compute(asm_fd.matrix());
        gsMatrix<real_t> u_fd = sol_fd->solve(asm_fd.rhs());
        real_t loss_fd = u_fd.norm();
        dL_dp_fd(k) = (loss_fd - loss) / fd_eps;
        patches.patch(0).coefs().data()[k] -= fd_eps;
    }

    gsInfo << "  Coeff  |   AD         |   FD         |   Rel. Error\n";
    gsInfo << "  -------+--" << std::string(12,'-') << "-+--" << std::string(12,'-') << "-+--" << std::string(12,'-') << "\n";
    for (index_t k = 0; k < nCheck; ++k)
    {
        real_t rel = (dL_dp_fd(k) != 0.0)
                   ? math::abs(dL_dp(k) - dL_dp_fd(k)) / math::abs(dL_dp_fd(k))
                   : math::abs(dL_dp(k));
        gsInfo << "  " << std::setw(5) << k
               << "  | " << std::setw(12) << dL_dp(k)
               << " | " << std::setw(12) << dL_dp_fd(k)
               << " | " << std::setw(12) << rel << "\n";
    }

    gsInfo << "\n=== Done ===\n";
#else
    gsInfo << "This example requires GISMO_AUTODIFF_BACKWARD.\n";
#endif
    return 0;
}
