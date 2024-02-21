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
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    gsBoundaryConditions<> bc;
    // TODO
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    for (int r =0; r < numRefine; ++r)
        dbasis.uniformRefine();

    for (size_t p = 0; p!=dbasis.nBases(); p++)
        gsInfo<<"Basis "<<p<<": "<<dbasis.basis(p).maxDegree()<<"\n";

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
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space w = A.getSpace(dbasis);

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
    index_t maxIt = 10;

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
    A.initSystem();
    gsMatrix<> tmp = gsMatrix<>::Random(A.numDofs(),1);
    Cold = tmp.array()*0.1/2;
    Cold.array() += 0.5;
    dCold.setZero(A.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = 1e-4;

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 16);

    for (index_t step = 0; step!=maxSteps; step++)
    {
        // Predictor
        Cnew = Cold;
        dCnew = (delta-1)/delta * dCold;

        Q0norm = 1;
        Qnorm = 10;
        gsInfo<<"Time "<<step*dt<<" ("<<step<<"/"<<maxSteps<<"):\n";

        for (index_t it = 0; it!= maxIt; it++)
        {
            Calpha  = Cold  + alpha_f * ( Cnew  - Cold );
            dCalpha = dCold + alpha_m * ( dCnew - dCold);
            A.assemble(residual * meas(G));
            Q = A.rhs();

            if (it == 0) Q0norm = Q.norm();
            else         Qnorm = Q.norm();

            gsInfo<<"\tIteration "<<it<<": Qnorm = "<<Qnorm<<"; Q0norm = "<<Q0norm<<"; Qnorm/Q0norm = "<<Qnorm/Q0norm<<"\n";
            if (Qnorm/Q0norm < tol)
                break;

            A.assembleJacobian( residual * meas(G), dc );
            K_m = alpha_m * A.matrix();

            A.assembleJacobian( residual * meas(G), c );
            K_f = alpha_f * delta * dt * A.matrix();

            K = K_m + K_f;

            solver.compute(K);
            dCupdate = solver.solve(-Q);


            dCnew += dCupdate;
            Cnew  += delta*dt*dCupdate;
        }

        Cold = Cnew;
        dCold = dCnew;

        //! [Export visualization in ParaView]
        if (plot)
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
