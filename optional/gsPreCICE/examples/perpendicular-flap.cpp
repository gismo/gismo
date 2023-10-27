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
// #include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 0;
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
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(-0.05,0.0,0.05,1.0));

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    // Set degree
    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();
    numRefine = 0;

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";
    
    real_t rho = 3000;
    real_t E = 40000;
    real_t nu = 0.3;
    real_t mu = E / (2.0 * (1.0 + nu));
    real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

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
    std::vector<patchSide> couplingInterfaces(3);
    couplingInterfaces[0] = patchSide(0,boundary::east);
    couplingInterfaces[1] = patchSide(0,boundary::north);
    couplingInterfaces[2] = patchSide(0,boundary::west);


    // Set the dimension of the points
    gsMatrix<> nodes;
    // Start iteration over elements
    gsVector<> tmp;

    gsOptionList quadOptions = A.options();

    index_t quadSize = 0;

    for (std::vector<patchSide>::const_iterator it = couplingInterfaces.begin(); it!=couplingInterfaces.end(); it++)
    {
        // Get a domain iterator on the coupling interface
        typename gsBasis<real_t>::domainIter domIt = bases.basis(it->patch).makeDomainIterator(it->side());

        // First obtain the size of all quadrature points
        typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
        for (; domIt->good(); domIt->next() )
        {
            QuRule = gsQuadrature::getPtr(bases.basis(it->patch), quadOptions,it->side().direction());
            quadSize+=QuRule->numNodes();
        }
    }

    // Initialize parametric coordinates
    gsMatrix<> uv(patches.domainDim(),quadSize);
    // Initialize physical coordinates
    gsMatrix<> xy(patches.targetDim(),quadSize);

    // Grab all quadrature points
    index_t offset = 0;

    for (std::vector<patchSide>::const_iterator it = couplingInterfaces.begin(); it!=couplingInterfaces.end(); it++)
    {
        // Get a domain iterator on the coupling interface
        typename gsBasis<real_t>::domainIter domIt = bases.basis(it->patch).makeDomainIterator(it->side());
        typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
        for (domIt->reset(); domIt->good(); domIt->next())
        {
            QuRule = gsQuadrature::getPtr(bases.basis(it->patch), quadOptions,it->side().direction());
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           nodes, tmp);
            uv.block(0,offset,patches.domainDim(),QuRule->numNodes()) = nodes;

            gsMatrix<> tmp2;
            patches.patch(it->patch).eval_into(nodes,tmp2);
            xy.block(0,offset,patches.targetDim(),QuRule->numNodes()) = patches.patch(it->patch).eval(nodes);
            offset += QuRule->numNodes();
        }
    }

    // Define precice interface
    gsPreCICE<real_t> interface("Solid", precice_config);
    interface.addMesh("Solid-Mesh",xy);
    real_t precice_dt = interface.initialize();

    index_t meshID = interface.getMeshID("Solid-Mesh");
    index_t dispID = interface.getDataID("Displacement",meshID);
    index_t forceID = interface.getDataID("Force",meshID);

// ----------------------------------------------------------------------------------------------

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Coupling side
    gsPreCICEVectorFunction<real_t> g_C(&interface,meshID,forceID,patches,patches.geoDim());
    // Add all BCs
    // Coupling interface
    bcInfo.addCondition(0, boundary::north,  condition_type::neumann , &g_C);
    bcInfo.addCondition(0, boundary::east,  condition_type::neumann  , &g_C);
    bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_C);
    // Bottom side (prescribed temp)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0, false, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0, false, 1);
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------

    // // source function, rhs
    // gsConstantFunction<> g(0.,0.,2);

    // // source function, rhs
    // gsConstantFunction<> gravity(0.,0.,2);

    // // creating mass assembler
    // gsMassAssembler<real_t> massAssembler(patches,bases,bcInfo,gravity);
    // massAssembler.options().setReal("Density",rho);
    // massAssembler.assemble();
    // gsDebugVar(massAssembler.matrix().toDense());

    // gsElasticityAssembler<real_t> assembler(patches,bases,bcInfo,g);
    // assembler.options().setReal("YoungsModulus",E);
    // assembler.options().setReal("PoissonsRatio",nu);
    // assembler.assemble();
    // gsDebugVar(assembler.matrix().toDense());

    // Time integration coefficient (0.0 = explicit, 1.0 = implicit)
    real_t theta = 1.0;

    // Set the geometry map
    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(bases,patches.geoDim());

    // Set the solution
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    // Assemble mass matrix
    u.setup(bcInfo, dirichlet::homogeneous, 0);
    A.initSystem();
    A.assemble( rho * u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();

    // Assemble stiffness matrix (NOTE: also adds the dirichlet BCs inside the matrixs)
    u.setup(bcInfo, dirichlet::l2Projection, 0);
    A.initSystem();

    auto v = u;
    A.assemble(
                0.5 * mu * meas(G) *
                ( ijac(v,G) + ijac(v,G).cwisetr() )
                %
                ( ijac(u,G) +  ijac(u,G).cwisetr()).tr()   // 18x18
                +
                lambda * meas(G) *
                ( ijac( v, G).trace()  *  ijac( u, G).trace().tr() )
                );

    gsDebugVar(A.matrix().toDense());

    auto g_Neumann = A.getBdrFunction(G);
    // Assemble Neumann BC term ( to .rhs())
    A.assembleBdr(bcInfo.get("Neumann"), v * g_Neumann * meas(G)  );

    gsSparseMatrix<> K = A.matrix();

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    // Time step
    real_t dt = 0.1;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    solVector.setZero(A.numDofs(),1);

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

    auto g_C_tmp = A.getCoeff(g_C,G);

    // Time integration loop
    while (interface.isCouplingOngoing())
    {
        // read temperature from interface
        if (interface.isReadDataAvailable())
        {
            A.initSystem();
            // Then assemble
            // auto g_Neumann = A.getBdrFunction(G);
            auto g_Neumann = g_C_tmp;
            // gsDebugVar("Assembling neumann ...");
            A.assembleBdr(bcInfo.get("Neumann"), v * g_Neumann * nv(G).norm()  );
            // A.assembleBdr(bcInfo.get("Neumann"), v * nv(G)  );
            // gsDebugVar(A.rhs().transpose());
            // gsDebugVar("Assembling neumann finished");

            // gsDebugVar(uv.col(13));
            // gsDebugVar(ev.eval(g_Neumann,uv.col(13)));
            // gsDebugVar(ev.eval(meas(G),uv.col(13)));
            // gsDebugVar(ev.eval(g_Neumann * meas(G),uv.col(13)));

            // gsDebugVar(ev.eval(v * g_Neumann * meas(G),uv.col(13)));

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
            gsMatrix<> result(patches.geoDim(),uv.cols());
            for (index_t k=0; k!=uv.cols(); k++)
            {
                // gsDebugVar(ev.eval(nv(G),uv.col(k)));
                result.col(k) = ev.eval(u_sol,uv.col(k));
            }
            interface.writeBlockVectorData(meshID,dispID,xy,result);
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
