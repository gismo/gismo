/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    real_t theta = 1.5;
    real_t lambda = 1;
    real_t L0 = 1;
    real_t M0 = 1;
    real_t dt = 1e-3;
    index_t maxSteps = 10;

    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    bool last = false;
    std::string fn("pde/poisson2d_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addReal( "t", "dt","dt parameter",dt);
    cmd.addReal( "T", "theta","Theta parameter",theta);
    cmd.addReal( "l", "lambda","lambda parameter",lambda);
    cmd.addReal( "L", "L0","L0 parameter",L0);
    cmd.addReal( "M", "M0","M0 parameter",M0);
    cmd.addInt ( "N", "Nsteps", "Number of time steps",  maxSteps );
    cmd.addInt ( "m", "PlotMod", "Modulo for plotting",  plotmod );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> data(fn);
    gsMultiPatch<> mp;
    //mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    gsMultiBasis<> mb;
    gsSparseMatrix<> cf;
    gsMappedBasis<2, real_t> mbasis;

    data.getFirst(mp);
    data.getFirst(mb);
    data.getFirst(cf);
    

    gsBoundaryConditions<> bc;
    // TODO
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    //gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    //dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    //for (int r =0; r < numRefine; ++r)
    //   dbasis.uniformRefine();

    //for (size_t p = 0; p!=dbasis.nBases(); p++)
        //gsInfo<<"Basis "<<p<<": "<<dbasis.basis(p).maxDegree()<<"\n";

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(mb);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);
     
    // Set the discretization space
    space w = A.getSpace(mbasis);

    mbasis.init(mb, cf);

    // Solution vector and solution variable
    gsMatrix<> Cnew, Calpha, Cold;
    gsMatrix<> dCnew,dCalpha,dCold, dCupdate;

    solution c = A.getSolution(w, Calpha); // C
    solution dc = A.getSolution(w, dCalpha); // \dot{C}

    // real_t N2   = L0*L0/lambda;
    real_t N2   = 41.7313;
    auto mu_c = 1.0 / (2.0*theta) * (c / (1.0-c).val()).log() + 1 - 2*c;
    auto dmu_c= 1.0 / (2.0*theta) * igrad(c,G) / (c - c*c).val() - 2.0 * igrad(c,G);

    auto M_c  = M0 * c * (1.0-c.val());
    auto dM_c = M0 * igrad(c,G) - 2.0 * M0 * igrad(c,G);

    auto residual = w*dc + N2*M_c.val()*igrad(w,G)*dmu_c.tr() +
                    ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + M_c.val() * ilapl(w,G)*ilapl(c,G);

    //! [Problem setup]

    // Define linear solver
    gsSparseSolver<>::CGDiagonal solver;

    // TIME INTEGRATION
    // constants
    real_t rho_inf = 0.5;
    real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    real_t alpha_f = 1 / (1+rho_inf);
    real_t delta   = 0.5 + alpha_m - alpha_f;
    // time stepping options
    index_t maxIt = 200;

    gsMatrix<> Q;
    gsSparseMatrix<> K, K_m, K_f;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    // Setup the space (compute Dirichlet BCs)
    w.setup(bc, dirichlet::l2Projection, 0);

    // Initialize the system
    gsMatrix<> tmp = gsMatrix<>::Random(A.numDofs(),1);
    Cold = tmp.array()*0.1/2; //random uniform variable in [-0.05,0.05]
    Cold.array() += 0.45;
    dCold.setZero(A.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = 1e-4;

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 16);

    real_t t_rho = 0.9;
    real_t t_err = 1;
    index_t lmax = 10;
    real_t TOL = 1e-3;
    std::vector<gsMatrix<>> Csols(2);

    real_t tmp_alpha_m = 1;
    real_t tmp_alpha_f = 1;
    real_t tmp_delta   = 1;

    real_t time = 0;
    bool converged = false;
    for (index_t step = 0; step!=maxSteps; step++)
    {
        for (index_t dt_it = 0; dt_it != lmax; dt_it++)
        {
            gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
            tmp_alpha_m = tmp_alpha_f = tmp_delta = 1;
            for (index_t k = 0; k!=2; k++)
            {
                converged = false;
                std::string method = (k==0) ? "Backward Euler " : "Generalized Alpha ";
                // ==================================================================================
                // Predictor
                Cnew = Cold;
                dCnew = (tmp_delta-1)/tmp_delta * dCold;

                Q0norm = 1;
                Qnorm = 10;

                for (index_t it = 0; it!= maxIt; it++)
                {
                    A.initSystem();
                    A.initVector(1);
                    Calpha  = Cold  + tmp_alpha_f * ( Cnew  - Cold );
                    dCalpha = dCold + tmp_alpha_m * ( dCnew - dCold);
                    A.assemble(residual * meas(G));
                    Q = A.rhs();

                    if (it == 0) Q0norm = Q.norm();
                    else         Qnorm = Q.norm();

                    // gsInfo<<"\t\tIteration "<<it<<": Qnorm = "<<Qnorm<<"; Q0norm = "<<Q0norm<<"; Qnorm/Q0norm = "<<Qnorm/Q0norm<<"\n";
                    if (Qnorm/Q0norm < tol)
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

                    A.assembleJacobian( residual * meas(G), dc );
                    K_m = tmp_alpha_m * A.matrix();

                    A.assembleJacobian( residual * meas(G), c );
                    K_f = tmp_alpha_f * tmp_delta * dt * A.matrix();

                    K = K_m + K_f;

                    solver.compute(K);
                    dCupdate = solver.solve(-Q);


                    dCnew += dCupdate;
                    Cnew  += tmp_delta*dt*dCupdate;
                }
                if (!converged)
                    break;
                // ==================================================================================
                tmp_alpha_m = alpha_m;
                tmp_alpha_f = alpha_f;
                tmp_delta = delta;

                Csols[k] = Cnew; // k=0: BE, k=1: alpha
            }

            if (converged)
            {
                t_err = (Csols[0] - Csols[1]).norm() / (Csols[1]).norm();
                dt *= t_rho * math::sqrt(TOL / t_err);
                if (t_err < TOL)
                    break;
            }
            else
                dt *= t_rho;
        }

        time += dt;
        Cold = Cnew;
        dCold = dCnew;

        //! [Export visualization in ParaView]
        if (plot && step % plotmod==0)
        {
            // Calpha = Cnew;
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




    // solver.compute( A.matrix() );
    // solVector = solver.solve(A.rhs());

    return EXIT_SUCCESS;

}// end main
