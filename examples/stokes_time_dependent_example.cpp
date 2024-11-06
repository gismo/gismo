#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
        // field IDs
        constexpr index_t PRESSURE_ID = 0;
        constexpr index_t VELOCITY_ID = 1;
        // field dimensions
        constexpr index_t PRESSURE_DIM = 1;
        real_t dt = 0.01; // Time step size
        index_t n_steps = 10; // Number of time steps

        // Setup values for timing
        double setup_time(0), assembly_time_ls(0), solving_time_ls(0),
          plotting_time(0);
        gsStopwatch timer;
        timer.restart();

        //! [Parse command line]
        bool plot = true;
        bool last = false;
        bool compute_error{false};
        std::string file_name("pde/stokes_quadCircle_intf.xml");

        gsCmdLine cmd("Solving a two-dimensional Stokes problem.");
        cmd.addString( "f", "file", "Input XML file", file_name );
        cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
        cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
        cmd.addSwitch("compute-error",
                    "Evaluate the error with respect to the analytical solution "
                    "(evaluation with default options and default file required)",
                    compute_error);

        
        // Material constants
        // kinematic viscosity
        real_t viscosity{1e-6};
        
        // Mesh options
        index_t numElevate = 0;
        cmd.addInt("e", "degreeElevate", "Number of uniform degree elevations",
                 numElevate);
        index_t numRefine = 0;
        cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement loops",
                 numRefine);

        //! [Parse command line]
        try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

        //! [Read input file]

        gsFileData<> file_data(file_name);
        gsInfo << "Loaded file "<< file_data.lastPath() <<"\n";

        gsMultiPatch<> multi_patch;
        file_data.getId(0, multi_patch); // id=0: Multipatch domain

        gsBoundaryConditions<> boundary_conditions;
        file_data.getId(2, boundary_conditions); // id=2: boundary conditions
        boundary_conditions.setGeoMap(multi_patch);
        gsInfo<<"Boundary conditions:\n"<< boundary_conditions <<"\n";

        //! Define function bases for velocity and pressure
        gsMultiBasis<> function_basis_velocity(multi_patch, true);//true: poly-splines (not NURBS)
        gsMultiBasis<> function_basis_pressure(multi_patch, true);//true: poly-splines (not NURBS)

        // Elevate the degree and increase the degree one 
        // additional time for the velocity to obtain Taylor-Hood elements
        for ( size_t i = 0; i < function_basis_velocity.nBases(); ++ i )
        {
            function_basis_velocity[i].setDegreePreservingMultiplicity(function_basis_velocity.maxCwiseDegree() + numElevate + 1);
            function_basis_pressure[i].setDegreePreservingMultiplicity(function_basis_pressure.maxCwiseDegree() + numElevate);
            function_basis_velocity[i].reduceContinuity(1);
        }

        const int geometric_dimension = multi_patch.geoDim();

        // h-refine each basis
        for (int r =0; r < numRefine; ++r) {
            function_basis_velocity.uniformRefine();
            function_basis_pressure.uniformRefine();
        }

        // Output information
        gsInfo << "Summary Velocity:" << std::endl
               << "Patches:" << multi_patch.nPatches() 
               << ", \nMin-degree: " << function_basis_velocity.minCwiseDegree()
               << ", \nMax-degree: " << function_basis_velocity.maxCwiseDegree()
               << std::endl << std::endl;
        gsInfo << "Summary Pressure:" << std::endl
               << "Patches:" << multi_patch.nPatches() 
               << ", \nMin-degree: " << function_basis_pressure.minCwiseDegree()
               << ", \nMax-degree: " << function_basis_pressure.maxCwiseDegree()
               << std::endl << std::endl;
        
    #ifdef _OPENMP
        gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
    #endif

        //! [Problem setup]
        gsExprAssembler<> expression_assembler(2,2);

        gsInfo<<"Active options:\n"<< expression_assembler.options() <<"\n";

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        // Elements used for numerical integration
        expression_assembler.setIntegrationElements(function_basis_velocity);
        gsExprEvaluator<> ev(expression_assembler);

        // Set the geometry map
        geometryMap geometry_map = expression_assembler.getMap(multi_patch);

        // Set the discretization spaces
        space velocity_trial_space = expression_assembler.getSpace(
            function_basis_velocity, geometric_dimension, VELOCITY_ID
            );
        space pressure_trial_space = expression_assembler.getSpace(
            function_basis_pressure, PRESSURE_DIM, PRESSURE_ID
            );

    // Solution vector and solution variables
    gsMatrix<> solution_vector;
    solution pressure_sol_expression = expression_assembler.getSolution(pressure_trial_space, solution_vector);
    solution velocity_sol_expression = expression_assembler.getSolution(velocity_trial_space, solution_vector);

    // Set Dirichlet BCs for velocity
    velocity_trial_space.setup(boundary_conditions, dirichlet::l2Projection, 0);

    // Initialize the system
    expression_assembler.initSystem();
    setup_time += timer.stop();

    gsInfo << "Number of degrees of freedom : " << expression_assembler.numDofs() << std::endl;
    gsInfo << "Number of blocks in the system matrix : " << expression_assembler.numBlocks() << std::endl;

    // Define the matrices for the Stokes system
    auto phys_jacobian = ijac(velocity_trial_space, geometry_map);
    auto mass_matrix = velocity_trial_space * velocity_trial_space.tr() * meas(geometry_map);
    auto viscosity_matrix = viscosity * (phys_jacobian.cwisetr() % phys_jacobian.tr()) * meas(geometry_map);
    auto divergence_matrix = -idiv(velocity_trial_space, geometry_map) * pressure_trial_space.tr() * meas(geometry_map);
    auto divergence_transpose_matrix = pressure_trial_space * idiv(velocity_trial_space, geometry_map).tr() * meas(geometry_map);

    //-----------------------------------The following terms do not change during assembly----------------------------------------
    gsSparseMatrix<> M, system_matrix;
    // A, B, Bt;
    expression_assembler.assemble(mass_matrix);
    M = expression_assembler.matrix();

    expression_assembler.initMatrix();

    expression_assembler.assemble(mass_matrix/dt, viscosity_matrix,divergence_matrix,divergence_transpose_matrix);

    // A = expression_assembler.matrix();

    // expression_assembler.initMatrix();
    // expression_assembler.assemble(divergence_matrix);
    // B = expression_assembler.matrix();

    // expression_assembler.initMatrix();
    // expression_assembler.assemble(divergence_transpose_matrix);
    system_matrix = expression_assembler.matrix();
    gsDebugVar(system_matrix);

    const auto& force_vector = expression_assembler.rhs(); //stationary force vector

    gsDebugVar(force_vector);
    //-----------------------------------The above terms do not change during assembly----------------------------------------

    // Set up system matrix for the Stokes problem

    // Time-stepping loop
    gsParaviewCollection collection("ParaviewOutput/solution");
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 16);

    // Initialize solution_vector with the correct dimensions
    solution_vector.setZero(expression_assembler.numDofs(), 1); 


    // Time-stepping loop
    for (index_t t = 0; t < n_steps; ++t)
    {
        gsInfo << "Time step " << t + 1 << "/" << n_steps << std::endl;

        // Update right-hand side with the previous solution inside of the time loop
        // Before this line
        auto rhs = M * solution_vector / dt + force_vector;

        gsSparseSolver<>::BiCGSTABILUT solver;
        solver.compute(system_matrix);
        gsDebugVar(system_matrix);
        solution_vector = solver.solve(rhs);
        gsDebugVar(solution_vector);



        gsExprEvaluator<> expression_evaluator(expression_assembler);
        gsInfo << "Solution vector size: " << solution_vector.size() << "\n";
        gsInfo << "Expected DoFs: " << expression_assembler.numDofs() << "\n";
        gsInfo << "Expected RHS DoFs: " << solution_vector.size() << "\n";
        // gsInfo << "Expression evaluator size:" <<expression_evaluator<< "\n";

        // Visualization and output
        if (plot)
        {
            gsInfo << "\nStarting the paraview export ..." << std::flush;
            // gsParaviewCollection collection("ParaviewOutput/solution",
            //                                 &expression_evaluator);
            // collection.options().setSwitch("plotElements", true);
            // collection.options().setInt("plotElements.resolution", 16);
            // collection.newTimeStep(&multi_patch);
            // collection.addField(pressure_sol_expression, "pressure");
            // collection.addField(velocity_sol_expression, "velocity");
            collection.saveTimeStep();
        }
    }

    // Final output
    if (plot) collection.save();

    gsInfo << "\nTotal time: " << setup_time + assembly_time_ls + solving_time_ls + plotting_time << std::endl;

    return EXIT_SUCCESS;
}