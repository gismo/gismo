/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
               L. Venta Viñuela (UniPv)
    
    
    Run a simple Cahn-Hilliard example with an analytical initial condition "0.1 * cos(2*pi*x) * cos(2*pi*y)" (Strong enforcement) (Gomez et al., 2014)
    ./bin/cahn-hilliard_example --plot -N 80 --plot
    
    Run a simple Cahn-Hilliard example with an analytical initial condition "0.1 * cos(2*pi*x) * cos(2*pi*y)" (Nitsche) (Bracco et al., 2023)
    ./bin/cahn-hilliard_example --plot -N 80 --nitsche --plot

    Run a simple Cahn-Hilliard example with a random normal initial concentration distribution of mean 0.0 until (almost) equilibrium (Nitsche)
    ./bin/cahn-hilliard_example --plot -N 1000 --nitsche --initial --plot
    
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    real_t dt = 1e-3;
    index_t maxSteps = 10;

    bool output = true;

    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    bool last = false;
    real_t mean = 0.0;

    bool random = false;

    std::string fn("pde/cahn_hilliard_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addReal( "t", "dt","dt parameter",dt); // -t () or --dt ()
    cmd.addInt ( "N", "Nsteps", "Number of time steps",  maxSteps );
    cmd.addInt ( "p", "PlotMod", "Modulo for plotting",  plotmod );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
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
    real_t penalty  = CHopt.askReal("penalty",1e4*lambda);

    gsOptionList TIMEopt;
    fd.getId(4, TIMEopt); // id=4: time integrator options

    gsOptionList Aopt;
    fd.getId(5, Aopt); // id=5: assembler options
    //! [Read input file]

    // Determine maximum mesh size
    real_t hmax = 0;
    for (size_t p=0; p!=mp.nPatches(); p++)
        hmax = math::max(hmax, mp.basis(p).getMaxCellLength());

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        dbasis.uniformRefine();

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    // geometryMap G = A.getMap(surface);
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space w = A.getSpace(dbasis);

    // basis.init(dbasis, cf);

    // Solution vector and solution variable
    gsMatrix<> Cnew, Calpha, Cold;
    gsMatrix<> dCnew,dCalpha,dCold, dCupdate;

    solution c = A.getSolution(w, Calpha); // C
    solution dc = A.getSolution(w, dCalpha); // \dot{C}

    gsSparseMatrix<> K_nitsche; // empty variable

    // Assemble the Nitsche BC on the sides with Neumann condition
    A.initSystem();
    A.assembleBdr(bc.get("Neumann"), - lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                  penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                  lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
    K_nitsche = A.giveMatrix(); // .giveMatrix() moves the matrix A into K_nitche (avoids having two matrices A and K_nitsche)


    // auto mu_c = 1.0 / (2.0*theta) * (c / (1.0-c).val()).log() + 1 - 2*c;
    // auto dmu_c= 1.0 / (2.0*theta) * igrad(c,G) / (c - c*c).val() - 2.0 * igrad(c,G);

     // auto mu_c = pow(c,3).val() - c.val();
    // auto dmu_c= igrad(c,G) * (- 1.0 + 3.0 * (c*c).val());

    // auto mu_c= -c.val() * (1.0 - (c*c).val());
    // auto dmu_c = -igrad(c,G) * (1.0 - (c*c).val()) - c.val() * (1.0 - 2.0 * c.val() * igrad(c,G));
    // auto M_c  = M0 * c * (1.0-c.val());
    // auto dM_c = M0 * igrad(c,G) - 2.0 * M0 * igrad(c,G);

    // Derivatives of the double well potential (Gomez et al., 2008)
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val();
    auto dM_c = 0.0 * igrad(c,G);

    // auto residual = w*dc + M_c.val()*igrad(w,G)*dmu_c.tr() +
    //                 lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + M_c.val() * ilapl(w,G)*lambda*ilapl(c,G);
    // auto residual = w*dc + M_c.val()*igrad(w,G)*dmu_c.tr() +
    //                     lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + M_c.val() * ilapl(w,G)*lambda*ilapl(c,G);
    // // auto residual = w*dc + // M
    //                 igrad(w,G)  * (- 1.0 + 3.0 * (c*c).val()) * igrad(c,G).tr() + // F_bar
    //                 //igrad(w,G)  * (-1) * igrad(c,G).tr() + // F_bar
    //                 // lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + // term gradient mobility!
    //                 ilapl(w,G)*lambda*ilapl(c,G).val(); // K_laplacian

    auto residual = w*dc + // M
                    M_c.val() * igrad(w,G)  * dmu_c * igrad(c,G).tr() + // F_bar
                    M_c.val() * ilapl(w,G)*lambda*ilapl(c,G).val(); // K_laplacian
                    // lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + // term gradient mobility!

    //! [Problem setup]

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solver;
#   else
    gsSparseSolver<>::LU solver;
#endif

    // Generalized-alpha method parameters
    real_t rho_inf = TIMEopt.askReal("rho_inf",0.5);
    real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    real_t alpha_f = 1 / (1+rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;
    // time stepping options
    index_t maxIt = 50;

    gsMatrix<> Q, Q1, Q2;
    gsSparseMatrix<> K, K_m, K_f;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";

    // Setup the space (compute Dirichlet BCs)
    w.setup(bc, dirichlet::l2Projection, 0);

    gsInfo<<"Initial condition.."<<"\n";

    if (random)
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Random initial condition %%%%%%%%%%%%%%%%%%%%%%%%
        gsMatrix<> tmp = gsMatrix<>::Random(A.numDofs(),1);
        Cold = tmp.array()*CHopt.askReal("ampl"); //random uniform variable in [-0.05,0.05]
        Cold.array() += CHopt.askReal("mean"); // 0.45
    }
    else
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        Cold.setZero(A.numDofs(),1);
        real_t error = gsL2Projection<real_t>::projectFunction(dbasis,source,mp,tmp);  // 3rd arg has to be multipatch
        // gsInfo << "L2 projection error "<<error<<"\n";
        for (index_t i = 0; i < dbasis.basis(0).size(); i++)
            if (w.mapper().is_free(i))
                Cold(w.mapper().index(i),0) = tmp(i,0);
    }
    
    Calpha = Cold;
    dCold.setZero(A.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = 1e-4;

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 1000);

    real_t dt_old = dt;
    real_t t_rho = TIMEopt.askReal("t_rho",0.9);
    real_t t_err = 1;
    index_t lmax = 1;
    real_t TOL = 1e-3;
    std::vector<gsMatrix<>> Csols(2);

    real_t tmp_alpha_m = 1;
    real_t tmp_alpha_f = 1;
    real_t tmp_gamma   = 1;

    real_t time = 0;
    bool converged = false;

    A.initSystem(); // Initialize the system (outside the loops)

    for (index_t step = 0; step!=maxSteps; step++)
    {
        for (index_t dt_it = 0; dt_it != lmax; dt_it++)
        {
            gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
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
                    A.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
                    Calpha.noalias()  = Cold  + tmp_alpha_f * ( Cnew  - Cold );
                    dCalpha.noalias() = dCold + tmp_alpha_m * ( dCnew - dCold);
                    
                    A.assemble(residual * meas(G));
                    Q = A.rhs();

                    if (bc.get("Neumann").size()!=0)
                    {      
                        Q.noalias() += K_nitsche * Calpha; // add the residual term from Nitche (using the matrix )
                        // Old code lines:
                        // A.assembleBdr(bc.get("Neumann"), - igrad(w,G) * nv(G) * lambda * ilapl(c,G).val()); // consistency term
                        // A.assembleBdr(bc.get("Neumann"),  (igrad(w,G) * nv(G).normalized()) * hmax * penalty * (igrad(c,G) * nv(G)) ); // penalty term
                        // A.assembleBdr(bc.get("Neumann"), - lambda * ilapl(w,G) * igrad(c,G) * nv(G)); // symmetry term
                    }

                    if (it == 0) Q0norm = Q.norm();
                    else         Qnorm = Q.norm();

                    gsInfo<<"\t\tNR iter   "<<it<<": res = "<<Qnorm/Q0norm<<"\n";
                    
                    if (it>0 && Qnorm/Q0norm < tol)
                    {
                        gsInfo<<"\t\t"<<method<<"converged in "<<it<<" iterations\n";
                        converged = true;
                        break;
                    }
                    else if (it==maxIt-1)
                    {
                        gsInfo<<"\t\t"<<method<<"did not converge!\n";
                        converged = false;
                        break;
                    }

                    // // gsInfo<<"Assembly K_m\n";
                    // A.assembleJacobian( residual * meas(G), dc );
                    // K_m = tmp_alpha_m * A.matrix(); 

                    // // gsInfo<<"Assembly K_f\n";
                    // A.assembleJacobian( residual * meas(G), c );
                    // K_f = tmp_alpha_f * tmp_gamma * dt * A.matrix();

                    // A.assemble(M_c.val() * igrad(w,G).tr() * igrad(w,G) * (- 1.0 + 3.0 * (c*c).val()) +
                    //             lambda*ilapl(w,G).tr()*igrad(w,G)*dM_c.tr() +
                    //             M_c.val() * ilapl(w,G)*lambda*ilapl(w,G).tr());

                    // A.assemble(meas(G) * (igrad(w,G) * (- 1.0 + 3.0 * (c*c).val())* igrad(w,G).tr()  + // K_f1
                    //                     igrad(w,G) * ((6.0 * c.val()) * igrad(c,G).tr() * w.tr()) + // K_f2
                    //                     // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility
                    //                     lambda * ilapl(w,G) * ilapl(w,G).tr())); // K_laplacian

                    //%% Assembly of the tangent stiffness matrix (K_m and K_f simultaneously) %%
                    A.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
                    A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                                        (tmp_alpha_f * tmp_gamma * dt)* (dmu_c *igrad(w,G) * igrad(w,G).tr() + // K_f1
                                        ddmu_c * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                                        lambda * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian                    
                                        // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility
                    
                    K = A.matrix(); 

                    if (bc.get("Neumann").size()!=0)
                        K += (tmp_alpha_f * tmp_gamma * dt) * K_nitsche; // add the Nitsche term to the stiffness matrix

#ifdef GISMO_WITH_SUPERLU
                    if (0==k)
                        solver.analyzePattern(K);
                    solver.factorize(K);
#else
                    solver.compute(K);
#endif
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
            Calpha = Cnew;
            // collection.newTimeStep(&mp);
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
