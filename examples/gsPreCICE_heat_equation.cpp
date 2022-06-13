/** @file heatEquation2_example.cpp

    @brief Solves the heat equation using time-stepping

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore, A. Mantzaflaris
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEFunction.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    bool last = false;
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Testing the heat equation.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    //! [Read input file]

    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(-0.25,0.0,0.0,1.0));

    gsFunctionExpr<> f("0",2);

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
    patchSide couplingInterface(0,boundary::north);
    typename gsBasis<real_t>::domainIter domIt = bases.basis(couplingInterface.patch).makeDomainIterator(couplingInterface.side());
    index_t rows = patches.targetDim();
    gsMatrix<> nodes;
    // Start iteration over elements
    gsVector<> tmp;
    index_t k=0;

    index_t quadSize = 0;
    typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
    for (; domIt->good(); domIt->next(), k++ )
    {
        QuRule = gsQuadrature::getPtr(bases.basis(couplingInterface.patch), A.options());
        quadSize+=QuRule->numNodes();
    }
    gsMatrix<> uv(rows,quadSize);
    gsMatrix<> xy(rows,quadSize);

    index_t offset = 0;

    for (domIt->reset(); domIt->good(); domIt->next(), k++ )
    {
        QuRule = gsQuadrature::getPtr(bases.basis(couplingInterface.patch), A.options());
        // Map the Quadrature rule to the element
        QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                       nodes, tmp);
        uv.block(0,offset,rows,QuRule->numNodes()) = nodes;

        gsMatrix<> tmp2;
        patches.patch(couplingInterface.patch).eval_into(nodes,tmp2);
        gsDebugVar(nodes);
        gsDebugVar(tmp2);
        xy.block(0,offset,rows,QuRule->numNodes()) = patches.patch(couplingInterface.patch).eval(nodes);
        offset += QuRule->numNodes();
    }

    gsDebugVar(uv);
    gsDebugVar(xy);

    gsPreCICE<real_t> interface("Solid", precice_config);
    interface.addMesh("Solid-Mesh",xy);
    real_t precice_dt = interface.initialize();

    index_t meshID = interface.getMeshID("Solid-Mesh");
    index_t tempID = interface.getDataID("Temperature",meshID);
    index_t fluxID = interface.getDataID("Heat-Flux",meshID);

    gsDebugVar(meshID);
    gsDebugVar(fluxID);
    gsDebugVar(tempID);

    gsMatrix<> result;
    // interface.readBlockScalarData(meshID,tempID,xy,result);
    // gsDebugVar(result);


    gsPreCICEFunction<real_t> g_N(&interface,meshID,tempID,patches);

// ----------------------------------------------------------------------------------------------

    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> g_D("0",2);
    bcInfo.addCondition(0, boundary::north,  condition_type::neumann  , &g_N);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D);
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------

    real_t theta = 0.5;

    // Generate system matrix and load vector
    // gsInfo << "Assembling mass and stiffness...\n";
    
    // Set the geometry map
    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(bases);

    u.setup(bcInfo, dirichlet::interpolation, 0);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);


    A.initSystem();
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) );
    gsSparseMatrix<> K = A.matrix();

    A.initSystem();
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();

    // Enforce Neumann conditions to right-hand side
    auto g_Neumann = A.getBdrFunction(G);
    A.assembleBdr(bcInfo.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    real_t dt = 0.01;

    A.assemble( u * ff * meas(G) );
    solVector.resize(A.numDofs(),1);
    solVector.setZero();
    gsVector<> F = dt*A.rhs() + (M-dt*(1.0-theta)*K)*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;
    gsParaviewCollection collection("solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    while (interface.isCouplingOngoing())
    {
        // read temperature from interface
        if (interface.isReadDataAvailable())
        {
            // Then assemble
            A.assemble( u * ff * meas(G) );

            // temperature_values = interface.read_block_scalar_data(temperature_id, vertex_ids)
            // temperature_function = coupling_sample.asfunction(temperature_values)

            // sqr = coupling_sample.integral((ns.u - temperature_function)**2)
            // cons = nutils.solver.optimize('lhs', sqr, droptol=1e-15, constrain=cons0)
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
        gsInfo << "Solving timestep " << timestep*dt << ".\n";
        F = dt*A.rhs() + (M-dt*(1.0-theta)*K)*solVector;
        solVector = solver.compute(M + dt*theta*K).solve(F);

        // write heat fluxes to interface
        if (interface.isWriteDataRequired(dt))
        {
            // TODO
            // interface.writeBlockScalarData()
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
            // WRITE
            /*
                if (plot)
                {
                    std::string fileName = "solution_" + util::to_string(i);
                    ev.options().setSwitch("plot.elements", true);
                    ev.options().setInt("plot.npts", 1000);
                    ev.writeParaview( u_sol   , G, fileName);
                    for (size_t p=0; p!=patches.nPatches(); p++)
                    {
                      fileName = "solution_" + util::to_string(i) + std::to_string(p);
                      collection.addTimestep(fileName,i,".vts");
                    }
                }
            */
        }


    }

    // for ( int i = 1; i<=numSteps; ++i) // for all timesteps
    // {
    //     A.assemble( u * ff * meas(G) );
    //     F = dt*A.rhs() + (M-dt*(1.0-theta)*K)*solVector;
    //     // Compute the system for the timestep i
    //     gsInfo << "Solving timestep " << i*dt << ".\n";
    //     solVector = solver.compute(M + dt*theta*K).solve(F);

    //     if (plot)
    //     {
    //         std::string fileName = "solution_" + util::to_string(i);
    //         ev.options().setSwitch("plot.elements", true);
    //         ev.options().setInt("plot.npts", 1000);
    //         ev.writeParaview( u_sol   , G, fileName);
    //         for (size_t p=0; p!=patches.nPatches(); p++)
    //         {
    //           fileName = "solution_" + util::to_string(i) + std::to_string(p);
    //           collection.addTimestep(fileName,i,".vts");
    //         }
    //     }
    // }

    if (plot)
    {
        collection.save();
        gsFileManager::open("solution.pvd");
    }


    return  EXIT_SUCCESS;
}
