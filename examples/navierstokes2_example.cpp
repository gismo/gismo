/** @file stokes2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Stokes equation in 2D

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): V. Travnikova
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    // field IDs
    constexpr index_t PRESSURE_ID = 0;
    constexpr index_t VELOCITY_ID = 1;
    // field dimensions
    constexpr index_t PRESSURE_DIM = 1;
    // number of solution and test spaces
    //constexpr index_t NUM_TRIAL = 2;
    //constexpr index_t NUM_TEST = 2;

    // Setup values for timing
    double setup_time(0), assembly_time_ls(0), solving_time_ls(0),
      plotting_time(0);
    gsStopwatch timer;
    timer.restart();

    //! [Parse command line]
    bool plot = false;
    bool last = false;
    bool compute_error{false};
    std::string file_name("pde/stokes_quadCircle.xml");

    gsCmdLine cmd("Tutorial on solving a two-dimensional Stokes problem.");
    cmd.addString( "f", "file", "Input XML file", file_name );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("compute-error",
                "Evaluate the error with respect to the analytical solution "
                "(evaluation with default options and default file required)",
                compute_error);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // Material constants
    real_t viscosity{1e-3};
    //real_t density{1e3};
    //cmd.addReal("v", "visc", "Viscosity", viscosity);
    

    // Mesh options
    index_t numElevate = 0;
    //cmd.addInt("e", "degreeElevate", "Number of uniform degree elevations",
    //         numElevate);
    index_t numRefine = 0;
    //cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement loops",
    //         numRefine);

    //! [Parse command line]

    //! [Read input file]

    gsFileData<> file_data(file_name);
    gsInfo << "Loaded file "<< file_data.lastPath() <<"\n";

    gsMultiPatch<> multi_patch;
    file_data.getId(0, multi_patch); // id=0: Multipatch domain

    //gsFunctionExpr<> source_function_expression;
    //file_data.getId(1, source_function_expression); // id=1: source function
    //gsInfo<<"Source function "<< source_function_expression << "\n";

    gsBoundaryConditions<> boundary_conditions;
    file_data.getId(2, boundary_conditions); // id=2: boundary conditions
    boundary_conditions.setGeoMap(multi_patch);
    gsInfo<<"Boundary conditions:\n"<< boundary_conditions <<"\n";

    //gsFunctionExpr<> ref_sol_expression;
    //file_data.getId(3, ref_sol_expression); // id=3: reference solution

    gsOptionList assembler_options;
    file_data.getId(4, assembler_options); // id=4: assembler options

    //! Define function bases for velocity and pressure
    gsMultiBasis<> function_basis_velocity(multi_patch, true);//true: poly-splines (not NURBS)
    gsMultiBasis<> function_basis_pressure(multi_patch, true);//true: poly-splines (not NURBS)

    // Elevate the degree and increase the degree one 
    // additional time for the velocity to obtain Taylor-Hood elements
    function_basis_velocity.setDegree( function_basis_velocity.maxCwiseDegree() + numElevate + 1);
    function_basis_pressure.setDegree( function_basis_pressure.maxCwiseDegree() + numElevate);

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
    expression_assembler.setOptions(assembler_options);

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
    

    // Set the source term
    //auto source_function = expression_assembler.getCoeff(source_function_expression, geometry_map);

    // Solution vector and solution variables
    gsMatrix<> solution_vector;
    solution pressure_sol_expression = 
        expression_assembler.getSolution(pressure_trial_space, solution_vector);
    solution velocity_sol_expression = 
        expression_assembler.getSolution(velocity_trial_space, solution_vector);

    // Set Dirichlet BCs for velocity
    velocity_trial_space.setup(boundary_conditions, dirichlet::l2Projection, 0);
    // Initialize interfaces for pressure field?

    // Initialize the system
    expression_assembler.initSystem();
    setup_time += timer.stop();

    gsInfo << "Number of degrees of freedom : " << expression_assembler.numDofs()
         << std::endl;
    gsInfo << "Number of blocks in the system matrix : "
         << expression_assembler.numBlocks() << std::endl;
    
    
    // Assembly
    gsInfo << "Starting assembly of linear system ..." << std::flush;
    timer.restart();
    
    // Compute the system matrix and right-hand side
    // Compute the system matrix and right-hand side
    auto phys_jacobian = ijac(velocity_trial_space, geometry_map);
    auto bilin_conti = pressure_trial_space * idiv(velocity_trial_space, geometry_map).tr() * meas(geometry_map);
    auto bilin_press = -idiv(velocity_trial_space, geometry_map) * pressure_trial_space.tr() * meas(geometry_map);
    auto bilin_mu_1 = viscosity * (phys_jacobian.cwisetr() % phys_jacobian.tr()) *
                    meas(geometry_map);
    auto bilin_mu_2 =
      viscosity * (phys_jacobian % phys_jacobian.tr()) * meas(geometry_map);
    
    auto trilin_conv = velocity_trial_space.cwisetr() * (velocity_trial_space * idiv(velocity_trial_space, geometry_map)) * meas(geometry_map);

    expression_assembler.assemble(bilin_conti, bilin_press, bilin_mu_1, bilin_mu_2, trilin_conv);

    assembly_time_ls += timer.stop();
    gsInfo << "\t\tFinished" << std::endl;

    assembly_time_ls += timer.stop();

    // gsDebugVar(A.matrix().toDense());
    // gsDebugVar(A.rhs().transpose()   );

    gsInfo << "\t\tFinished" << std::endl;

  ///////////////////
  // Linear Solver //
  ///////////////////

  gsInfo << "Solving the linear system of equations ..." << std::flush;
  timer.restart();

  const auto& system_matrix = expression_assembler.matrix();
  const auto& rhs_vector = expression_assembler.rhs();

  // Initialize linear solver
  gsSparseSolver<>::CGDiagonal solver;
  
  solver.compute(system_matrix);
  solution_vector = solver.solve(rhs_vector);

  solving_time_ls += timer.stop();
  gsInfo << "\tFinished" << std::endl;

  // Export and visualization
  gsExprEvaluator<> expression_evaluator(expression_assembler);
  gsInfo << "\nStarting the paraview export ..." << std::flush;
  timer.restart();

  gsParaviewCollection collection("ParaviewOutput/solution",
                                    &expression_evaluator);
  collection.options().setSwitch("plotElements", true);
  collection.options().setInt("plotElements.resolution", 16);
  collection.newTimeStep(&multi_patch);
  collection.addField(pressure_sol_expression, "pressure");
  collection.addField(velocity_sol_expression, "velocity");
  collection.saveTimeStep();
  collection.save();

  plotting_time += timer.stop();
  gsInfo << "\tFinished" << std::endl;
  
  gsInfo << "Starting the xml export ..." << std::flush;

  // Export pressure
  gsMatrix<> solution_pressure;
  gsFileData<> output_pressure;
    pressure_sol_expression.extractFull(
        solution_pressure);  // patch-wise solution with BCs
    output_pressure << solution_pressure;
    output_pressure.save("pressure_field.xml");

    // Export velocity
    gsMatrix<> solution_velocity;
    gsFileData<> output;
    velocity_sol_expression.extractFull(
        solution_velocity);  // patch-wise solution with BCs
    output << solution_velocity;
    output.save("velocity_field.xml");

// User output infor timings
  gsInfo << "\n\nTotal time: "
         << setup_time + assembly_time_ls + solving_time_ls + plotting_time
         << std::endl;
  gsInfo << "                       Setup: " << setup_time << std::endl;
  gsInfo << "      Assembly Linear System: " << assembly_time_ls << std::endl;
  gsInfo << "       Solving Linear System: " << solving_time_ls << std::endl;
  gsInfo << "                    Plotting: " << plotting_time << std::endl
         << std::flush;
        

    return EXIT_SUCCESS;

}// end main
