/** @file gsCahnHilliardAssembler_multipatch_example.cpp

    @brief Tutorial on how to the gsCahnHilliardAssembler class on multi-patches

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

    index_t verbose = 1;

    /* TIME INTEGRATION OPTIONS */
    // Generalized-alpha method parameters

    // time stepping options
    /* NONLINEAR SOLVER OPTIONS */
    index_t maxIt = 50;

    std::string fn = "pde/cahn_hilliard_multipatch_bvp.xml";
    std::string fn_g;

    gsCmdLine cmd("Tutorial on solving a Cahn-Hilliard problem on multi-patches.");
    cmd.addReal( "t", "dt","dt parameter",dt); // -t () or --dt ()
    cmd.addInt ( "N", "Nsteps", "Number of time steps",  maxSteps );
    cmd.addInt ( "p", "PlotMod", "Modulo for plotting",  plotmod );
    cmd.addInt ( "v", "verbose", "Verbosity level",  verbose );
    cmd.addString( "f", "file", "Input XML file for the options", fn );
    cmd.addString( "g", "geo",  "Input XML file for the geometry (if not in the other file)", fn_g);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ENSURE(fn.empty()==false,"Options file not provided");

    gsMultiPatch<> mp;
    gsMultiBasis<> mb;
    gsSparseMatrix<> cf;
    gsMappedBasis<2, real_t> mbasis;

    // Get options
    gsFileData<> fd;
    if (!fn_g.empty())
        fd.read(fn_g);
    else
        fd.read(fn);

    fd.getFirst(mp);
    gsInfo<<"Read geometry from file "<<fd.lastPath()<<"\n";
    gsInfo<<"* patches:\n"<<mp<<"\n";
    fd.getFirst(mb);
    gsInfo<<"* bases:\n"<<mb<<"\n";
    fd.getFirst(cf);
    gsInfo<<"* mapper "<<cf.rows()<<" x "<<cf.cols()<<"\n";
    mbasis.init(mb, cf);

    fd.clear();
    if (!fn_g.empty())
        fd.read(fn);

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

    // Set boundary conditions:
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    // Set all BCs to homogeneous clamped
    for ( gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        bc.addCondition( *bit, condition_type::clamped,0);

    // Define the assembler
    gsCahnHilliardAssembler<real_t> assembler(mp, mb, bc);
    assembler.setSpaceBasis(mbasis);
    assembler.options().setReal("Lambda",lambda);
    assembler.options().setReal("Penalty",penalty);
    assembler.options().setInt("Continuity",0);
    assembler.options().setSwitch("AssembleWeakBCs",false);
    // assembler.options().setReal("M0",M0);
    assembler.initialize();

    // gsSparseMatrix<> K_nitsche;
    // if (bc.get("Weak Clamped").size()!=0 && !assembler.options().getSwitch("AssembleWeakBCs"))
    // {
    //     assembler.assembleNitscheMatrix();
    //     assembler.matrix_into(K_nitsche); // .matrix_into() moves the matrix A into K_nitsche (avoids having two matrices A and K_nitsche)
    // }
    gsSparseMatrix<> M;
    assembler.assembleMassMatrix();
    assembler.matrix_into(M);

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
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
    gsMappedSpline<2,real_t> cnew, dcnew;

    // %%%%%%%%%%%%%%%%%%%%%%%% Random initial condition %%%%%%%%%%%%%%%%%%%%%%%%
    gsMatrix<> tmp = gsMatrix<>::Random(assembler.numDofs(),1);
    Cold = tmp.array()*CHopt.askReal("ampl",0.005); //random uniform variable in [-0.05,0.05]
    Cold.array() += CHopt.askReal("mean",0.0); // 0.45

    Calpha = Cold;
    dCold.setZero(assembler.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;

    // USE A MESH FOR PLOTTING
    gsSurfMesh mesh = mp.toMesh();
    auto pid = mesh.get_vertex_property<index_t>("v:patch");
    auto aid = mesh.get_vertex_property<index_t>("v:anchor");
    auto field = mesh.add_vertex_property<real_t>("v:field");

    // Time step settings
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
                    assembler.assembleJacobian(cnew, dcnew);
                    assembler.matrix_into(K);
                    K *= (tmp_alpha_f * tmp_gamma * dt);
                    K += tmp_alpha_m * M;

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

            gsExprEvaluator<> ev;
            auto c = ev.getVariable(cnew);

            gsVector<> pt;
            for (auto vit = mesh.vertices_begin(); vit < mesh.vertices_end(); ++vit)
            {
                index_t k = pid[*vit];
                pt = mp.patch(k).basis().anchor(aid[*vit]);
                field[*vit] = ev.eval( c , pt, k ).value();//any expression
            }
            std::string fileName = "ParaviewOutput/solution"+util::to_string(step);
            gsWriteParaview(mesh,fileName, {"v:field"} );
        }
    }

    if (!plot)
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return EXIT_SUCCESS;

}// end main
