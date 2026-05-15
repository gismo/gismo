/** @file gsCahnHilliardAssembler_example.cpp

    @brief Tutorial on how to the gsCahnHilliardAssembler class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
               L. Venta Viñuela (UniPv)
*/

//! [Include namespace]
#include <gismo.h>
#include <gsAssembler/gsCahnHilliardAssembler.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    real_t dt = 1e-3;
    index_t maxSteps = 10;

    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    index_t verbose = 1;
    bool random = false;

    /* TIME INTEGRATION OPTIONS */
    // Generalized-alpha method parameters

    // time stepping options
    /* NONLINEAR SOLVER OPTIONS */
    index_t maxIt = 50;

    std::string fn("pde/cahn_hilliard_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addReal( "t", "dt","dt parameter",dt); // -t () or --dt ()
    cmd.addInt ( "N", "Nsteps", "Number of time steps",  maxSteps );
    cmd.addInt ( "p", "PlotMod", "Modulo for plotting",  plotmod );
    cmd.addInt ( "v", "verbose", "Verbosity level",  verbose );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("random", "Random initial condition of the CH problem", random);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> source;
    fd.getId(1, source); // id=1: initial condition function
    gsInfo<<"Initial condition function "<< source << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList CHopt;
    fd.getId(3, CHopt); // id=3: reference solution

    real_t theta    = CHopt.askReal("theta",1.5);
    real_t lambda   = CHopt.askReal("lambda",1/(32*pow(EIGEN_PI,2)));
    real_t M0       = CHopt.askReal("M0",0.005);
    real_t penalty  = 1e4*lambda;

    gsOptionList TIMEopt;
    fd.getId(4, TIMEopt); // id=4: time integrator options

    real_t rho_inf = TIMEopt.askReal("rho_inf",0.5);
    real_t tol = TIMEopt.askReal("tol",1e-4);
    real_t t_rho = TIMEopt.askReal("t_rho",0.9);
    real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    real_t alpha_f = 1 / (1+rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;

    gsOptionList Aopt;
    fd.getId(5, Aopt); // id=5: assembler options
    //! [Read input file]

    //! [Refinement]
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

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::LU solver;
#endif


    gsMatrix<> Q, Qnitsche;
    gsSparseMatrix<> K;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";
    gsInfo<<"Initial condition.."<<"\n";

    // Solution vector and solution variable
    gsMatrix<> Cnew, Calpha, Cold;
    gsMatrix<> dCnew,dCalpha,dCold, dCupdate;
    gsMultiPatch<> cnew;
    if (random)
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Random initial condition %%%%%%%%%%%%%%%%%%%%%%%%
        gsMatrix<> tmp = gsMatrix<>::Random(assembler.numDofs(),1);
        Cold = tmp.array()*CHopt.askReal("ampl",0.005); //random uniform variable in [-0.05,0.05]
        Cold.array() += CHopt.askReal("mean",0.0); // 0.45
    }
    else
    {
        // // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),source,tmp);
        // real_t error = gsL2Projection<real_t>::project(dbasis,source,mp,tmp);  // 3rd arg has to be multipatch
        // if (verbose>0) gsInfo << "L2 projection error "<<error<<"\n";

        gsMultiPatch<> cold;
        cold.addPatch(*dbasis.basis(0).makeGeometry(give(tmp)));
        assembler.constructSolution(cold,Cold);
    }

    Calpha = Cold;
    dCold.setZero(assembler.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;


    gsExprEvaluator<> ev;
    ev.setIntegrationElements(dbasis);
    auto c = ev.getVariable(cnew);
    gsParaviewCollection collection("ParaviewOutput/solution",&ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);

    real_t dt_old = dt;
    real_t t_err = 1;
    index_t lmax = 1;
    std::vector<gsMatrix<>> Csols(2);

    real_t tmp_alpha_m = 1;
    real_t tmp_alpha_f = 1;
    real_t tmp_gamma   = 1;

    real_t time = 0;
    bool converged = false;

    assembler.initialize(); // Initialize the system (outside the loops)

    for (index_t step = 0; step!=maxSteps; step++)
    {
        for (index_t dt_it = 0; dt_it != lmax; dt_it++)
        {
            if (verbose>0) gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
            tmp_alpha_m = tmp_alpha_f = tmp_gamma = 1;

            for (index_t k = 0; k!=2; k++)
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
                    assembler.assembleResidual(cnew);
                    assembler.rhs_into(Q);
                    Q += M * dCalpha; // add mass matrix contribution: M*dC

                    if (bc.get("Weak Clamped").size()!=0 && !assembler.options().getSwitch("AssembleWeakBCs"))
                    {
                        assembler.assembleNitscheVector(cnew);
                        assembler.rhs_into(Qnitsche);
                        Q.noalias() += Qnitsche; // add the residual term from Nitche (using the matrix )
                    }

                    if (it == 0) Q0norm = Q.norm();
                    else         Qnorm = Q.norm();

                    if (verbose==2) gsInfo<<"\t\tNR iter   "<<it<<": res = "<<Qnorm/Q0norm<<"\n";

                    if (it>0 && Qnorm/Q0norm < tol)
                    {
                        if (verbose>0) gsInfo<<"\t\t"<<method<<"converged in "<<it<<" iterations\n";
                        converged = true;
                        break;
                    }
                    else if (it==maxIt-1)
                    {
                        if (verbose>0) gsInfo<<"\t\t"<<method<<"did not converge!\n";
                        converged = false;
                        break;
                    }


                    // Assembly of the tangent stiffness matrix (K_m and K_f simultaneously) %%
                    assembler.assembleJacobian(cnew);
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

                // %% For time step adaptivity %%
                // Csols[k] = Cnew; // k=0: BE, k=1: alpha
            }

            // %%%%%%%%%% For time step adaptivity %%%%%%%%%%
            // if (converged)
            // {
            //     t_err = (Csols[0] - Csols[1]).norm() / (Csols[1]).norm();
            //     dt_old = dt;
            //     // dt *= t_rho * math::sqrt(TOL / t_err);
            //     if (t_err < TOL)
            //         break;
            // }
            // else
            // {
            //     dt_old = dt;
            //     // dt *= t_rho;
            // }
        }

        time += dt_old;
        Cold = Cnew;
        dCold = dCnew;

        //! [Export visualization in ParaView]
        if (plot && step % plotmod==0)
        {
            assembler.constructSolution(Cnew,  cnew);
            collection.newTimeStep(&mp);
            collection.addField(c,"numerical solution");
            collection.saveTimeStep();
        }
    }

    if (plot)
    {
        collection.save();
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return EXIT_SUCCESS;

}// end main
