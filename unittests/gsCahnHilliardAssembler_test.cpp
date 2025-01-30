/** @file gsCahnHilliardAssembler_test.cpp

    @brief Tests for gsCahnHilliardAssembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#include "gismo_unittest.h"
#include <gsAssembler/gsCahnHilliardAssembler.h>

SUITE(gsCahnHilliardAssembler_test)
{
    TEST(FullSimulation)
    {
        const real_t dt = 1e-3;
        const index_t maxSteps = 100;
        const index_t numRefine = 4;
        const index_t numElevate = 1;

        // Geometry
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineSquare());

        // Initial condition
        gsFunctionExpr<> source("0.1 * cos(2*pi*x) * cos(2*pi*y)",2);

        // Boundary conditions
        gsBoundaryConditions<> bc;
        bc.addCondition(boundary::west,condition_type::weak_clamped,0);
        bc.addCondition(boundary::east,condition_type::weak_clamped,0);
        bc.addCondition(boundary::south,condition_type::weak_clamped,0);
        bc.addCondition(boundary::north,condition_type::weak_clamped,0);
        bc.setGeoMap(mp);

        // Parameters for CH
        // real_t theta = 1.5;
        real_t lambda = 1.5e-3;
        // real_t M0 = 0.005;
        real_t penalty = 1e4*lambda;

        // Time integrator options
        real_t rho_inf = 0.5;
        real_t tol     = 1e-4;
        // real_t t_rho   = 0.9;
        real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
        real_t alpha_f = 1 / (1+rho_inf);
        real_t gamma   = 0.5 + alpha_m - alpha_f;

        // Basis
        gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

        // Elevate and p-refine the basis to order p + numElevate
        // where p is the highest degree in the bases
        dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

        // h-refine each basis
        for (int r =0; r < numRefine; ++r)
            dbasis.uniformRefine();

        gsCahnHilliardAssembler<real_t> assembler(mp, dbasis, bc);
        assembler.options().setReal("Lambda",lambda);
        assembler.options().setReal("Penalty",penalty);
        assembler.options().setSwitch("AssembleWeakBCs",false);
        // assembler.options().setReal("M0",M0);
        assembler.initialize();

        gsSparseMatrix<> K_nitsche, M;
        if (bc.get("Weak Clamped").size()!=0 && !assembler.options().getSwitch("AssembleWeakBCs"))
        {
            assembler.assembleNitscheMatrix();
            assembler.matrix_into(K_nitsche); // .matrix_into() moves the matrix A into K_nitsche (avoids having two matrices A and K_nitsche)
        }
        assembler.assembleMassMatrix();
        assembler.matrix_into(M);

        gsSparseSolver<>::CGDiagonal solver;

        gsMatrix<> Q, Qnitsche;
        gsSparseMatrix<> K;

        // Solution vector and solution variable
        gsMatrix<> Cnew, Calpha, Cold;
        gsMatrix<> dCnew,dCalpha,dCold, dCupdate;
        gsMultiPatch<> cnew, dcnew;

        // // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),source,tmp);
        // real_t error = gsL2Projection<real_t>::project(dbasis,source,mp,tmp);  // 3rd arg has to be multipatch

        gsMultiPatch<> cold;
        cold.addPatch(*dbasis.basis(0).makeGeometry(give(tmp)));
        assembler.constructSolution(cold,Cold);


        Calpha = Cold;
        dCold.setZero(assembler.numDofs(),1);

        real_t Q0norm = 1, Qnorm = 10;
        real_t dt_old = dt;
        // real_t t_err = 1;
        index_t lmax = 1;
        std::vector<gsMatrix<>> Csols(2);

        real_t tmp_alpha_m = 1;
        real_t tmp_alpha_f = 1;
        real_t tmp_gamma   = 1;

        real_t time = 0;
        bool converged = false;
        index_t maxIt = 100;

        assembler.initialize(); // Initialize the system (outside the loops)

        for (index_t step = 0; step!=maxSteps; step++)
        {
            for (index_t dt_it = 0; dt_it != lmax; dt_it++)
            {
                gsTestInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
                tmp_alpha_m = tmp_alpha_f = tmp_gamma = 1;

                for (index_t k = 0; k!=1; k++)
                {
                    converged = false;
                    std::string method = (k==0) ? "Backward Euler " : "Generalized Alpha ";
                    // ==================================================================================
                    // Predictor
                    Cnew = Cold;
                    dCnew = (tmp_gamma-1)/tmp_gamma * dCold;

                    Q0norm = 1;
                    Qnorm = 10;

                    for (index_t it = 0; it!= maxIt; it++)
                    {
                        Calpha.noalias()  = Cold  + tmp_alpha_f * ( Cnew  - Cold );
                        dCalpha.noalias() = dCold + tmp_alpha_m * ( dCnew - dCold);

                        assembler.constructSolution(Calpha,  cnew);
                        assembler.constructSolution(dCalpha,dcnew);
                        assembler.assembleResidual(cnew, dcnew);
                        assembler.rhs_into(Q);

                        if (bc.get("Weak Clamped").size()!=0 && !assembler.options().getSwitch("AssembleWeakBCs"))
                        {
                            assembler.assembleNitscheVector(cnew,dcnew);
                            assembler.rhs_into(Qnitsche);
                            Q.noalias() += Qnitsche; // add the residual term from Nitche (using the matrix )
                        }

                        if (it == 0) Q0norm = Q.norm();
                        else         Qnorm = Q.norm();

                        if (it>0 && Qnorm/Q0norm < tol)
                        {
                            gsTestInfo<<"\t\t"<<method<<"converged in "<<it<<" iterations\n";
                            converged = true;
                            break;
                        }
                        else if (it==maxIt-1)
                        {
                            gsTestInfo<<"\t\t"<<method<<"did not converge!\n";
                            converged = false;
                            break;
                        }


                        // Assembly of the tangent stiffness matrix (K_m and K_f simultaneously) %%
                        assembler.assembleJacobian(cnew, dcnew);
                        assembler.matrix_into(K);
                        K *= (tmp_alpha_f * tmp_gamma * dt);
                        K += tmp_alpha_m * M;

                        if (bc.get("Weak Clamped").size()!=0 && !assembler.options().getSwitch("AssembleWeakBCs"))
                            K += (tmp_alpha_f * tmp_gamma * dt) * K_nitsche; // add the Nitsche term to the stiffness matrix

                        solver.compute(K);
                        dCupdate = solver.solve(-Q);

                        dCnew += dCupdate;
                        Cnew.noalias() += (tmp_gamma*dt)*dCupdate;
                    }
                    if (!converged)
                        break;

                    // %% Switch to generalized-alpha parameters (k=1)
                    tmp_alpha_m = alpha_m;
                    tmp_alpha_f = alpha_f;
                    tmp_gamma = gamma;
                }
            }

            time += dt_old;
            Cold = Cnew;
            dCold = dCnew;
        }

        gsExprEvaluator<> ev;
        ev.setIntegrationElements(dbasis);
        auto c = ev.getVariable(cnew);

        gsMatrix<> pts;
        // In the following points, the solution should be close to 1:
        // (0,0), (0,1), (1,0), (1,1), (0.5,0.5)
        pts.resize(2,5);
        pts << 0, 0, 1, 1, 0.5,
               0, 1, 0, 1, 0.5;

        for (index_t i = 0; i < pts.cols(); i++)
            CHECK_CLOSE(ev.eval(c,pts.col(i)).value(),1.0,1e-3);

        // In the following points, the solution should be close to -1:
        // (0.5,0), (0,0.5), (0.5,1), (1,0.5)
        pts.resize(2,4);
        pts << 0.5, 0, 0.5, 1,
               0, 0.5, 1, 0.5;

        for (index_t i = 0; i < pts.cols(); i++)
            CHECK_CLOSE(ev.eval(c,pts.col(i)).value(),-1,1e-3);

        // The integral of c should be close to 0
        CHECK_CLOSE(math::abs(ev.integral(c)),0.,1e-5);
    }
}
