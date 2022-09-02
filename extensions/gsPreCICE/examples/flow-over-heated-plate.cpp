/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

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
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    //! [Read input file]

    // Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,-0.25,1.0,0.0));

    // Set external heat-flux to zero
    gsFunctionExpr<> f("0",2);

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    // Set degree
    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();
    numRefine = 0;

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";
    
    // Set heat conduction coefficient
    real_t k_temp = 100;
    // Set bottom wall temp
    real_t u_wall = 310;

// ----------------------------------------------------------------------------------------------
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(bases);

    gsExprEvaluator<> ev(A);

// ----------------------------------------------------------------------------------------------
    // Set the interface for the precice coupling
    patchSide couplingInterface(0,boundary::north);

    // Get a domain iterator on the coupling interface
    typename gsBasis<real_t>::domainIter domIt = bases.basis(couplingInterface.patch).makeDomainIterator(couplingInterface.side());

    // Set the dimension of the points
    gsMatrix<> nodes;
    // Start iteration over elements
    gsVector<> tmp;

    gsOptionList quadOptions = A.options();

    // First obtain the size of all quadrature points
    index_t quadSize = 0;
    typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
    for (; domIt->good(); domIt->next() )
    {
        QuRule = gsQuadrature::getPtr(bases.basis(couplingInterface.patch), quadOptions,couplingInterface.side().direction());
        quadSize+=QuRule->numNodes();
    }

    // Initialize parametric coordinates
    gsMatrix<> uv(patches.domainDim(),quadSize);
    // Initialize physical coordinates
    gsMatrix<> xy(patches.targetDim(),quadSize);

    // Grab all quadrature points
    index_t offset = 0;
    for (domIt->reset(); domIt->good(); domIt->next())
    {
        QuRule = gsQuadrature::getPtr(bases.basis(couplingInterface.patch), quadOptions,couplingInterface.side().direction());
        // Map the Quadrature rule to the element
        QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                       nodes, tmp);
        uv.block(0,offset,patches.domainDim(),QuRule->numNodes()) = nodes;

        gsMatrix<> tmp2;
        patches.patch(couplingInterface.patch).eval_into(nodes,tmp2);
        xy.block(0,offset,patches.targetDim(),QuRule->numNodes()) = patches.patch(couplingInterface.patch).eval(nodes);
        offset += QuRule->numNodes();
    }

    // Define precice interface
    gsPreCICE<real_t> interface("Solid", precice_config);
    interface.addMesh("Solid-Mesh",xy);
    real_t precice_dt = interface.initialize();

    index_t meshID = interface.getMeshID("Solid-Mesh");
    index_t tempID = interface.getDataID("Temperature",meshID);
    index_t fluxID = interface.getDataID("Heat-Flux",meshID);

// ----------------------------------------------------------------------------------------------

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(u_wall,2);
    // Homogeneous Neumann
    gsConstantFunction<> g_N(0.0,2);
    // Coupling side
    gsPreCICEFunction<real_t> g_C(&interface,meshID,tempID,patches);
    // Add all BCs
    // Isolated Neumann sides
    bcInfo.addCondition(0, boundary::east,  condition_type::neumann  , &g_N, 0, false, 0);
    bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_N, 0, false, 0);
    // Bottom side (prescribed temp)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0, false, 0);
    // Coupling interface
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g_C, 0, false, 0);
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------

    // Time integration coefficient (0.0 = explicit, 1.0 = implicit)
    real_t theta = 1.0;

    // Set the geometry map
    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(bases);
    u.setup(bcInfo, dirichlet::homogeneous, 0);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    // Assemble mass matrix
    A.initSystem();
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();

    // Assemble stiffness matrix (NOTE: also adds the dirichlet BCs inside the matrixs)
    A.initSystem();
    A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G) );
    gsSparseMatrix<> K = A.matrix();

    // Enforce Neumann conditions to right-hand side
    auto g_Neumann = A.getBdrFunction(G);
    A.assembleBdr(bcInfo.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    // Time step
    real_t dt = 0.01;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    gsConstantFunction<> uwall_fun(u_wall,2);
    auto uwall = A.getCoeff(uwall_fun, G);
    // RHS of the projection
    A.assemble( u * uwall * meas(G) );
    solver.compute(M);
    solVector = solver.solve(A.rhs());

    // Initialize the RHS for assembly
    A.initVector();
    u.setup(bcInfo, dirichlet::l2Projection, 0);

    // Assemble the RHS
    gsVector<> F = dt*A.rhs() + (M-dt*(1.0-theta)*K)*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;

    // Define the solution collection for Paraview
    gsParaviewCollection collection("solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t time = 0;

    // Plot initial solution
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
    }

    // Time integration loop
    while (interface.isCouplingOngoing())
    {
        // read temperature from interface
        if (interface.isReadDataAvailable())
        {
            u.setup(bcInfo, dirichlet::l2Projection, 0); // NOTE:
            A.initSystem();
            A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G) );
            K = A.matrix();

            // Then assemble
            auto g_Neumann = A.getBdrFunction(G);
            A.assembleBdr(bcInfo.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
            A.assemble( u * ff * meas(G) );
            F = theta*dt*A.rhs() + (1.0-theta)*dt*F0 + (M-dt*(1.0-theta)*K)*solVector;
        }

        // save checkpoint
        if (interface.isActionRequired(interface.actionWriteIterationCheckpoint()))
        {
            F_checkpoint = F0;
            timestep_checkpoint = timestep;
            interface.markActionFulfilled(interface.actionWriteIterationCheckpoint());
        }

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);

        // solve gismo timestep
        gsInfo << "Solving timestep " << timestep*dt << "...";
        solVector = solver.compute(M + dt*theta*K).solve(F);
        gsInfo<<"Finished\n";
        // write heat fluxes to interface
        if (interface.isWriteDataRequired(dt))
        {
            gsMatrix<> result(1,uv.cols()), tmp;
            for (index_t k=0; k!=uv.cols(); k++)
            {
                // gsDebugVar(ev.eval(nv(G),uv.col(k)));
                tmp = ev.eval(k_temp * igrad(u_sol,G),uv.col(k));
                // Only exchange y component
                result(0,k) = -tmp.at(1);
            }
            interface.writeBlockScalarData(meshID,fluxID,xy,result);
        }

        // do the coupling
        precice_dt = interface.advance(dt);

        // advance variables
        timestep += 1;
        F0 = F;

        if (interface.isActionRequired(interface.actionReadIterationCheckpoint()))
        {
            F0 = F_checkpoint;
            timestep = timestep_checkpoint;
            interface.markActionFulfilled(interface.actionReadIterationCheckpoint());
        }
        else
        {
            time += dt;
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
            }
        }


    }

    if (plot)
    {
        collection.save();
    }


    return  EXIT_SUCCESS;
}
