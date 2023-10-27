/** @file heat-equation-coupling.cpp

    @brief Heat equation participant for a double coupled heat equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEFunction.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    short_t side = 0;
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("s","side", "Patchside of interface", side);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]

    gsMultiPatch<> patches;
    if (side==0) //left
        patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0));
    else if (side==1) //right
        patches.addPatch(gsNurbsCreator<>::BSplineRectangle(1.0,0.0,2.0,1.0));
    else
        GISMO_ERROR("Side unknown");

    real_t alpha = 3;
    real_t beta  = 1.3;
    real_t time  = 0;
    real_t k_temp = 1;

    // Set external heat-flux to zero
    gsConstantFunction<> f(beta-2-2*alpha,2);
    gsFunctionExpr<> u_ex("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);

    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();
    numRefine = 0;

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

// ----------------------------------------------------------------------------------------------
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(bases);

    gsExprEvaluator<> ev(A);

// ----------------------------------------------------------------------------------------------
    boxSide couplingSide;
    if (side==0) //left
        couplingSide = boxSide(2);
    else if (side==1) //right
        couplingSide = boxSide(1);
    else
        GISMO_ERROR("Side unknown");

    patchSide couplingInterface(0,couplingSide);
    typename gsBasis<real_t>::domainIter domIt = bases.basis(couplingInterface.patch).makeDomainIterator(couplingInterface.side());
    index_t rows = patches.targetDim();
    gsMatrix<> nodes;
    // Start iteration over elements
    gsVector<> tmp;
    index_t k=0;

    gsOptionList quadOptions = A.options();
    // NEED A DIFFERENT RULE FOR dirichlet::interpolation --> see: gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());
    /*
        quadOptions.addInt("quRule","Quadrature rule [1:GaussLegendre,2:GaussLobatto]",1);
        quadOptions.addReal("quA","Number of quadrature points: quA*deg + quB",1.0);
        quadOptions.addInt("quB","Number of quadrature points: quA*deg + quB",1);
    */

    index_t quadSize = 0;
    typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
    for (; domIt->good(); domIt->next(), k++ )
    {
        QuRule = gsQuadrature::getPtr(bases.basis(couplingInterface.patch), quadOptions,couplingInterface.side().direction());
        quadSize+=QuRule->numNodes();
    }
    gsMatrix<> uv(rows,quadSize);
    gsMatrix<> xy(rows,quadSize);

    index_t offset = 0;

    for (domIt->reset(); domIt->good(); domIt->next(), k++ )
    {
        QuRule = gsQuadrature::getPtr(bases.basis(couplingInterface.patch), quadOptions,couplingInterface.side().direction());
        // Map the Quadrature rule to the element
        QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                       nodes, tmp);
        uv.block(0,offset,rows,QuRule->numNodes()) = nodes;

        gsMatrix<> tmp2;
        patches.patch(couplingInterface.patch).eval_into(nodes,tmp2);
        xy.block(0,offset,rows,QuRule->numNodes()) = patches.patch(couplingInterface.patch).eval(nodes);
        offset += QuRule->numNodes();
    }

    std::string membername;
    if (side==0) //left
        membername = "Dirichlet";
    else if (side==1) //right
        membername = "Neumann";
    else
        GISMO_ERROR("Side unknown");

    gsPreCICE<real_t> interface(membername, precice_config);
    std::string meshname = membername + "-Mesh";
    interface.addMesh(meshname,xy);
    real_t precice_dt = interface.initialize();

    index_t meshID = interface.getMeshID(meshname);
    index_t tempID = interface.getDataID("Temperature",meshID);
    index_t fluxID = interface.getDataID("Heat-Flux",meshID);

// ----------------------------------------------------------------------------------------------

    gsBoundaryConditions<> bcInfo;
    gsPreCICEFunction<real_t> g_CD(&interface,meshID,(side==0 ? tempID : fluxID),patches);
    gsPreCICEFunction<real_t> g_CN(&interface,meshID,(side==0 ? tempID : fluxID),patches);
    gsFunction<> * g_C = &u_ex;

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u_ex, 0, false, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &u_ex, 0, false, 0);
    if (side==0)
    {
        bcInfo.addCondition(0, couplingSide,condition_type::dirichlet  , g_C, 0, false, 0);
        bcInfo.addCondition(0, couplingSide.opposite(),condition_type::dirichlet, &u_ex, 0, false, 0);
    }
    else
    {
        bcInfo.addCondition(0, couplingSide,condition_type::neumann  , &g_CN);
        bcInfo.addCondition(0, couplingSide.opposite(),condition_type::dirichlet, &u_ex, 0, false, 0);
    }

    bcInfo.setGeoMap(patches);
    // gsDebugVar(bcInfo);

// ----------------------------------------------------------------------------------------------

    // Generate system matrix and load vector
    // gsInfo << "Assembling mass and stiffness...\n";

    // Set the geometry map
    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(bases);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector, solVector_ini;
    solution u_sol = A.getSolution(u, solVector);
    solution u_ini = A.getSolution(u, solVector);

    u.setup(bcInfo, dirichlet::homogeneous, 0);
    A.initSystem();
    gsDebugVar(A.numDofs());
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();


    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    real_t dt = 0.01;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    auto uex = A.getCoeff(u_ex, G);
    // RHS of the projection
    u.setup(bcInfo, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( u * u.tr() * meas(G), u * uex * meas(G) );
    solver.compute(A.matrix());
    solVector_ini = solVector = solver.solve(A.rhs());

        gsMatrix<> result(1,uv.cols()), tmp2;
    // Write initial data
    if (interface.isActionRequired(interface.actionWriteInitialData()))
    {
        for (index_t k=0; k!=uv.cols(); k++)
        {
            if (side==0)
            {
                gsWarn<<"Write the flux here!!!\n";
                tmp2 = ev.eval( - jac(u_sol) * nv(G).normalized(),uv.col(k));
                result(0,k) = tmp2.at(0);
            }
            else
            {
                tmp2 = ev.eval(u_sol,uv.col(k));
                result(0,k) = tmp2.at(0);
            }
        }
        // gsDebugVar(result);
        interface.writeBlockScalarData(meshID,side==0 ? fluxID : tempID,xy,result);
        interface.markActionFulfilled(interface.actionWriteInitialData());
    }
    interface.initialize_data();

    // interface.readBlockScalarData(meshID,side==0 ? tempID : fluxID,xy,result);
    // gsDebugVar(result);

    // Initialize the RHS for assembly
    if (side==0)
        g_C = &g_CD;
    gsDebugVar(bcInfo);
    A.initSystem();
    A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * uex * meas(G) );
    gsSparseMatrix<> K = A.matrix();

    // Assemble the RHS
    gsVector<> F = dt*A.rhs() + M*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;
    gsMatrix<> solVector_checkpoint = solVector;

    gsParaviewCollection collection("solution");
    gsParaviewCollection exact_collection("exact_solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t t_checkpoint = 0;
    if (plot)
    {
        std::string fileName = "solution_" + util::to_string(timestep);
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview( u_sol   , G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = "solution_" + util::to_string(timestep) + std::to_string(p);
          collection.addTimestep(fileName,time,".vts");
        }

        // fileName = "exact_solution_" + util::to_string(timestep);
        // ev.writeParaview( uex   , G, fileName);
        // for (size_t p=0; p!=patches.nPatches(); p++)
        // {
        //   fileName = "exact_solution_" + util::to_string(timestep) + std::to_string(p);
        //   exact_collection.addTimestep(fileName,time,".vts");
        // }

        ev.writeParaview( u_sol   , G, "initial_solution");

    }

    // time += precice_dt;
    while (interface.isCouplingOngoing())
    {
        // read temperature from interface
        // if (interface.isReadDataAvailable())
        // {
            u_ex = gsFunctionExpr<>("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);
            u.setup(bcInfo, dirichlet::l2Projection, 0); // NOTE:

            A.initSystem();
            A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
            K = A.matrix();
            auto g_Neumann = A.getBdrFunction(G);
            A.assembleBdr(bcInfo.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
            F = dt*A.rhs() + M*solVector;

            // gsMatrix<> result;
            // interface.readBlockScalarData(meshID,side==0 ? tempID : fluxID,xy,result);
            // gsDebugVar(result);
        // }

        // save checkpoint
        if (interface.isActionRequired(interface.actionWriteIterationCheckpoint()))
        {
            gsDebugVar("Write checkpoint");
            F_checkpoint = F0;
            t_checkpoint = time;
            timestep_checkpoint = timestep;
            solVector_checkpoint = solVector;
            interface.markActionFulfilled(interface.actionWriteIterationCheckpoint());
        }

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);

        // solve gismo timestep
        gsInfo << "Solving timestep " << timestep*dt << "...";
        solVector = solver.compute(M + dt*K).solve(F);
        gsInfo<<"Finished\n";
        // write heat fluxes to interface
        if (interface.isWriteDataRequired(dt))
        {
            gsMatrix<> result(1,uv.cols()), tmp;
            for (index_t k=0; k!=uv.cols(); k++)
            {
                // gsDebugVar(ev.eval(nv(G),uv.col(k)));
                // tmp = ev.eval(k_temp * igrad(u_sol,G),uv.col(k));
                // Only exchange y component
                // result(0,k) = -tmp.at(1);
                //
                //
                if (side==0)
                {
                    gsWarn<<"Write the flux here!!!\n";
                    tmp = ev.eval(jac(u_sol) * nv(G).normalized(),uv.col(k));
                    result(0,k) = tmp.at(0);
                }
                else
                {
                    tmp = ev.eval(u_sol,uv.col(k));
                    result(0,k) = tmp.at(0);
                }
            }
            // TODO
            interface.writeBlockScalarData(meshID,side==0 ? fluxID : tempID,xy,result);
        }

        // do the coupling
        precice_dt = interface.advance(dt);

        // advance variables
        time += dt;
        timestep += 1;
        F0 = F;

        if (interface.isActionRequired(interface.actionReadIterationCheckpoint()))
        {
            gsDebugVar("Read checkpoint");
            F0 = F_checkpoint;
            time = t_checkpoint;
            timestep = timestep_checkpoint;
            solVector = solVector_checkpoint;
            interface.markActionFulfilled(interface.actionReadIterationCheckpoint());
        }
        else
        {
            if (timestep % plotmod==0 && plot)
            {
                std::string fileName = "solution_" + util::to_string(timestep);
                ev.options().setSwitch("plot.elements", true);
                ev.options().setInt("plot.npts", 1000);
                ev.writeParaview( u_sol   , G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = "solution_" + util::to_string(timestep) + std::to_string(p);
                  collection.addTimestep(fileName,time,".vts");
                }

                // fileName = "exact_solution_" + util::to_string(timestep);
                // ev.writeParaview( uex   , G, fileName);
                // for (size_t p=0; p!=patches.nPatches(); p++)
                // {
                //   fileName = "exact_solution_" + util::to_string(timestep) + std::to_string(p);
                //   exact_collection.addTimestep(fileName,time,".vts");
                // }
            }
        }
    }

    if (plot)
    {
        collection.save();
        exact_collection.save();
    }


    return  EXIT_SUCCESS;
}
