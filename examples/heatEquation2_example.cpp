/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t steps = 10;
    bool last{false}, export_b64{false};
    std::string fn("pde/heat2d_square_ibvp1.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "N", "steps", "Number of time steps", steps );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addSwitch(
        "plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        dbasis.uniformRefine();

    // Set heat conduction coefficient
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);

    // Time integration coefficient (0.0 = explicit, 1.0 = implicit)
    real_t theta = 1.0;

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);
    u.setup(bc, dirichlet::homogeneous, 0);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    // Assemble mass matrix
    A.initSystem();
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();
    gsSparseMatrix<> K;

    // Enforce Neumann conditions to right-hand side
    auto g_Neumann = A.getBdrFunction(G);
    // Set initial condition
    solVector.setZero(A.numDofs(),1);

    // Plot initial solution
    gsParaviewCollection collection("solution");
    if (plot)
    {
        std::string fileName = "solution_" + util::to_string(0);
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview( u_sol   , G, fileName);
        for (size_t p=0; p!=mp.nPatches(); p++)
        {
          fileName = "solution_" + util::to_string(0) + std::to_string(p);
          collection.addTimestep(fileName,0.,".vts");
        }
    }

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;


    // Time step
    real_t time = 0.0;
    real_t maxTime = 0.1;
    real_t dt = maxTime / steps;

    // Assemble F
    gsMatrix<> F, Fpp;
    gsFunctionExpr<> * f_ptr = dynamic_cast<gsFunctionExpr<>*>(&f);
    f_ptr->set_t(0.0);
    A.assembleBdr(bc.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
    A.assemble( u * ff * meas(G) );
    F = A.rhs();

    // Assemble the stiffness matrix K; assumes constant BCs
    u.setup(bc, dirichlet::l2Projection, 0); // NOTE:
    A.initSystem();
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) );
    K = A.matrix();

    gsMatrix<> RHS;
    for (index_t i = 1; i < steps; ++i)
    {
        time = i*dt;
        f_ptr->set_t(time);

        // Assemble the RHS (can be time dependent)
        A.assembleBdr(bc.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
        A.assemble( u * ff * meas(G) );
        Fpp = A.rhs();
        RHS = theta*dt*Fpp + (1.0-theta)*dt*F + (M-dt*(1.0-theta)*K)*solVector;

        gsDebugVar(RHS.norm());

        // Solve the linear system
        solVector = solver.compute(M + dt*theta*K).solve(RHS);

        if (plot)
        {
            std::string fileName = "solution_" + util::to_string(i);
            ev.options().setSwitch("plot.elements", true);
            ev.options().setInt("plot.npts", 1000);
            ev.writeParaview( u_sol   , G, fileName);
            for (size_t p=0; p!=mp.nPatches(); p++)
            {
                fileName = "solution_" + util::to_string(i) + std::to_string(p);
                collection.addTimestep(fileName,time,".vts");
            }
        }

        F = Fpp;
    }

    if (plot)
    {
        collection.save();
    }

    return  EXIT_SUCCESS;
}
