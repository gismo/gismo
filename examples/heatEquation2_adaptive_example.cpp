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
    index_t numElevate = 1;
    index_t steps = 100;
    index_t plotmod = 1;

    real_t maxTime = 1;
    bool last{false}, export_b64{false};
    std::string fn("pde/heat2d_square_laser_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "N", "steps", "Number of time steps", steps );
    cmd.addReal("T", "Tmax", "", maxTime);
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addInt( "p","plotmod","",plotmod);
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

    gsFunctionExpr<> is;
    fd.getId(3, is); // id=3: initial solution

    gsFunctionExpr<> ms;
    fd.getId(4, ms);

    gsOptionList Aopt;
    fd.getId(5, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis;//true: poly-splines (not NURBS)
    // Cast all patches of the mp object to THB splines
    gsTHBSplineBasis<2,real_t> thb;
    for (size_t k=0; k!=mp.nPatches(); ++k)
    {
        gsTensorBSplineBasis<2,real_t> *b = dynamic_cast< gsTensorBSplineBasis<2,real_t> * > (&mp.basis(k));
        thb = gsTHBSplineBasis<2,real_t>(*b);
        dbasis.addBasis(&thb);
    }

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (size_t k=0; k!=dbasis.nBases(); ++k)
    {
        gsTHBSplineBasis<2,real_t> *b = dynamic_cast< gsTHBSplineBasis<2,real_t> * > (&dbasis.basis(k));
        gsMatrix<> bbox = b->support();
        for (index_t i = 0; i< numRefine; ++i)
            b->refineElements(b->asElements(bbox));
    }

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
    gsMatrix<> Unew;
    gsMultiPatch<> uold;
    solution u_new = A.getSolution(u, Unew);
    auto u_old = A.getCoeff(uold);

    // Assemble mass matrix
    A.initSystem();
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();
    gsSparseMatrix<> K;

    // Enforce Neumann conditions to right-hand side
    auto g_Neumann = A.getBdrFunction(G);
    // Set initial condition
    Unew.setZero(A.numDofs(),1);
    u_new.extract(uold);

    // Plot initial solution
    gsParaviewCollection collection("solution", &ev);
    if (plot)
    {
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints",1e3);
        collection.newTimeStep(&mp);

        collection.addField(u_new,"numerical solution");
        collection.addField(ff,"source");
        // collection.addField((ff-beta*dthetadt+alpha*d2thetadX).sqNorm(),"residual");

        collection.saveTimeStep();
    }

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;


    // Time step
    real_t time = 0.0;
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

    // auto t = 3;
    // auto X = {0,1,2}; // [:]

    // grad_expr(_expr, EIGEN EQUIVALENT OF : )
    // DOES NOT WORK
    // A.assemble( grad(u,t)*u * meas(G) + grad(u,X)*jac(G).ginv() * grad(u,X)*jac(G).ginv().tr() * meas(G) );

    gsAdaptiveMeshing<2,real_t> mesher(dbasis);
    mesher.options().setInt("RefineRule",3);
    mesher.options().setInt("CoarsenRule",3);
    mesher.options().setSwitch("Admissible",true);
    mesher.options().setReal("RefineParam",0.2);
    mesher.options().setReal("CoarsenParam",0.001);
    mesher.options().setInt("MaxLevel",8);
    mesher.getOptions();
    gsHBoxContainer<2,real_t> refine,coarsen;

    gsMatrix<> RHS;
    for (index_t i = 1; i < steps; ++i)
    {
        time = i*dt;
        // // Assemble the OLD RHS (can be time dependent)
        // u.setup(bc, dirichlet::l2Projection, 0); // NOTE:
        // A.initSystem();
        // A.assembleBdr(bc.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
        // A.assemble( u * ff * meas(G) );
        // F = A.rhs();

        // Advance the time
        f_ptr->set_t(time);

        // Assemble mass matrix
        u.setup(bc, dirichlet::homogeneous, 0);
        A.initSystem();
        A.assemble( u * u.tr() * meas(G));
        M = A.matrix();

        // Assemble the stiffness matrix K; assumes constant BCs
        u.setup(bc, dirichlet::l2Projection, 0); // NOTE:
        A.initSystem();
        A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) );
        K = A.matrix();

        // Assemble the RHS (can be time dependent)
        A.assembleBdr(bc.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
        A.assemble( u * ff * meas(G) );
        Fpp = A.rhs();


        // Compose the RHS
        // gsDebugVar(Fpp.size());
        // gsDebugVar(F.size());
        // gsDebugVar(M.rows());
        // gsDebugVar(K.rows());
        // gsDebugVar(Unew.rows());
        RHS = theta*dt*Fpp;// + (1.0-theta)*dt*F + (M-dt*(1.0-theta)*K)*Unew;

        // Solve the linear system
        Unew = solver.compute(M + dt*theta*K).solve(RHS);

        real_t error = ev.integralElWise( (ff - (u_new-u_old)/dt + ilapl(u_new,G)).sqNorm() * meas(G) );
        std::vector<real_t> errors = ev.elementwise();
        real_t energy_int = ev.integral((ilapl(u_new,G))*meas(G));
        real_t energy_tot = ev.integral(((u_new-u_old)/dt + ilapl(u_new,G))*meas(G));
        gsInfo<<"Current error = "<<error<<"; internal energy = "<<energy_int<<"; total energy = "<<energy_tot<<"\n";

        if (plot && i % plotmod==0)
        {
            collection.newTimeStep(&mp);

            collection.addField(u_new,"numerical solution");
            collection.addField(u_old,"old numerical solution");
            collection.addField(ff,"source");
            // collection.addField((ff-(u_new-u_old)/dt+ilapl(u_new,G)).sqNorm(),"residual");
            collection.addField((ff-(u_new-u_old)/dt+ilapl(u_new,G)).sqNorm(),"residual");
            collection.addField(((u_new-u_old)/dt+ilapl(u_new,G)).sqNorm(),"lhs");

            gsElementErrorPlotter<real_t> plotter(dbasis.basis(0),errors);
            auto err = ev.getVariable(plotter);
            // collection.addField(err,"element errors");

            collection.saveTimeStep();


            const gsField<> elemError_eh( mp.patch(0), plotter, true );
            gsWriteParaview(elemError_eh,"error_elem_ref" + std::to_string(i));
        }

        // F = Fpp;
        u_new.extract(uold);

        mesher.markRef_into(errors,refine);
        mesher.markCrs_into(errors,refine,coarsen);
        mesher.refine(refine);
        mesher.unrefine(coarsen);
        mesher.rebuild();

        // Uold = Unew;
    }

    if (plot)
    {
        collection.save();
    }

    return  EXIT_SUCCESS;
}
