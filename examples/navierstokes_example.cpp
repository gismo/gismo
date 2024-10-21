/** @file navierstokes2_example.cpp

    @brief Using the expression assembler to solve the Navier-Stokes equations in 2D

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
    constexpr index_t NUM_TRIAL = 2;
    constexpr index_t NUM_TEST = 2;

    // Setup values for timing
    double setup_time(0), assembly_time_ls(0), linear_solving_time_ls(0), nonlinear_solving_time_ls(0),
      plotting_time(0);
    gsStopwatch timer;
    timer.restart();

    //! [Parse command line]
    bool plot = false;
    bool last = false;
    bool compute_error{false};
    std::string file_name("pde/navierstokes_coarse_cavity.xml");

    gsCmdLine cmd("Solving a two-dimensional Navier-Stokes problem.");
    cmd.addString( "f", "file", "Input XML file", file_name );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // Material constants
    real_t viscosity{10.0};
    //real_t density{1e3};
    
    // Mesh options
    index_t numElevate = 0;
    cmd.addInt("e", "degreeElevate", "Number of uniform degree elevations",
             numElevate);
    index_t numRefine = 0;
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement loops",
             numRefine);

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
    

    //! [Problem setup]
    gsExprAssembler<> expression_assembler(NUM_TEST, NUM_TRIAL);

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
    solution pressure_sol_expression = 
        expression_assembler.getSolution(pressure_trial_space, solution_vector);
    solution velocity_sol_expression = 
        expression_assembler.getSolution(velocity_trial_space, solution_vector);


    // Set Dirichlet BCs for velocity
    //velocity_trial_space.setup(boundary_conditions, dirichlet::l2Projection, 0);
    velocity_trial_space.setup(boundary_conditions, dirichlet::interpolation, 0);

    solution_vector = gsMatrix<>(expression_assembler.numDofs(), 1);
  
    gsInfo<<"Active options:\n"<< expression_assembler.options() <<"\n";

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
    auto phys_jacobian = ijac(velocity_trial_space, geometry_map);
    auto bilin_conti = pressure_trial_space * idiv(velocity_trial_space, geometry_map).tr() * meas(geometry_map);
    auto bilin_press = -idiv(velocity_trial_space, geometry_map) * pressure_trial_space.tr() * meas(geometry_map);
    auto bilin_mu_1 = viscosity * (phys_jacobian.cwisetr() % phys_jacobian.tr()) * meas(geometry_map);
    //auto bilin_mu_2 =
    //  viscosity * (phys_jacobian % phys_jacobian.tr()) * meas(geometry_map);

    // TEST: define point to evaluate expressions at for debugging
    gsVector<> point(2);
    point << 0.8, 0.1468;
    //gsInfo << "\n eval point physical coordinates" << "\n" << ev.eval(geometry_map, point) << std::endl;
  
    //

    // Assemble the linear (Stokes) system 
    expression_assembler.assemble(bilin_conti, bilin_press, bilin_mu_1);

    assembly_time_ls += timer.stop();
    gsInfo << "\t\tFinished" << std::endl;

  ///////////////////
  // Linear Solver //
  ///////////////////

  gsInfo << "Solving the linear system of equations ..." << std::flush;
  timer.restart();

  //const auto& system_matrix = expression_assembler.matrix();
  //const auto& rhs_vector = expression_assembler.rhs();

  // Initialize linear solver
  gsSparseSolver<>::BiCGSTABILUT solver;
  //gsMatrix<> ker = expression_assembler.matrix().toDense().fullPivLu().kernel();
  //gsInfo << ker << '\n';
  
  // Compute system matrix
  solver.compute(expression_assembler.matrix());

  // Initial solution
  solution_vector = solver.solve(expression_assembler.rhs());

  gsInfo << "\n Stokes velocity \n" << ev.eval(velocity_sol_expression, point) << std::endl;
  velocity_sol_expression.printDetail(gsInfo);

  linear_solving_time_ls += timer.stop();
  gsInfo << "\tLinear solving done" << std::endl;

  timer.restart();

  //TEST: export the initial solution of Stokes
  gsExprEvaluator<> expression_evaluator_linear(expression_assembler);
  gsInfo << "\nStarting the paraview export of the Stokes solution used as initial solution for the nonlinear solver ..." << std::flush;

  gsParaviewCollection collection_linear("ParaviewOutput/solution_stokes",
                                   &expression_evaluator_linear);
  collection_linear.options().setSwitch("plotElements", true);
  collection_linear.options().setInt("plotElements.resolution", 16);
  collection_linear.newTimeStep(&multi_patch);
  collection_linear.addField(pressure_sol_expression, "pressure");
  collection_linear.addField(velocity_sol_expression, "velocity");
  collection_linear.saveTimeStep();
  collection_linear.save();

  // TODO: parse maxIter from command line
  int maxIter {5};
  int l;

  gsInfo << "Entering the nonlinear solver loop ..." << std::flush;

  // ! Newton method
  // expression_assembler.options().setInt("DirichletStrategy", 0); //switch off elimination
  // auto phys_jacobian_sol = ijac(velocity_sol_expression, geometry_map);
  //auto bilin_conti_sol = pressure_trial_space * idiv(velocity_sol_expression, geometry_map).tr() * meas(geometry_map);
  //auto bilin_press_sol = -idiv(velocity_trial_space, geometry_map) * pressure_sol_expression.tr() * meas(geometry_map);
  //auto bilin_mu_1_sol = viscosity * (phys_jacobian_sol.cwisetr() % phys_jacobian.tr()) * meas(geometry_map);
  //auto trilin_conv_1 = velocity_sol_expression.tr() * phys_jacobian * velocity_trial_space.tr() * meas(geometry_map);
  //auto trilin_conv_2 = velocity_trial_space * phys_jacobian_sol.tr() * velocity_sol_expression * meas(geometry_map);

  //auto trilin_rhs = velocity_trial_space.tr() * phys_jacobian * velocity_trial_space.tr() * meas(geometry_map);
  //auto residual = trilin_conv_1 + bilin_mu_1 + trilin_conv_2;// + bilin_press;// + bilin_conti;

  //gsInfo << "\n residual \n" << ev.eval(residual, point) << std::endl;
  //residual.printDetail(gsInfo);

  //expression_assembler.clearMatrix();
  //expression_assembler.clearRhs();

  //gsMatrix<> solvec0;
  // for (l=0; l<maxIter; ++l)
  // {
    //solvec0 = solution_vector;
    //timer.restart();
  //  expression_assembler.clearMatrix();
  //  expression_assembler.clearRhs();
  //  expression_assembler.assembleJacobian(residual, velocity_sol_expression);
    //expression_assembler.assemble(residual);

 // }
  
  
  //gsInfo << expression_assembler.matrix().toDense()<<"\n";
  //gsInfo << expression_assembler.rhs().transpose()<<"\n";
  //gsInfo << solution_vector.transpose()<<" \n";


  //! Picard method
  gsMatrix <> initial_solution;
  for (l=0; l<maxIter; ++l)
  {
    gsInfo<<" ----------------------------------------------------- \n";
    initial_solution = solution_vector;
    timer.restart();
    //velocity_trial_space.setup(boundary_conditions, dirichlet::l2Projection, 0);
    
   expression_assembler.clearMatrix();
   expression_assembler.clearRhs();
  
   auto phys_jacobian = ijac(velocity_trial_space, geometry_map);
   auto bilin_conti = pressure_trial_space * idiv(velocity_trial_space, geometry_map).tr() * meas(geometry_map);
   auto bilin_press = -idiv(velocity_trial_space, geometry_map) * pressure_trial_space.tr() * meas(geometry_map);
   auto bilin_mu_1 = viscosity * (phys_jacobian.cwisetr() % phys_jacobian.tr()) * meas(geometry_map);
   //auto bilin_mu_2 = viscosity * (phys_jacobian % phys_jacobian.tr()) * meas(geometry_map);

   //auto trilin_conv_1 = velocity_sol_expression.tr().val() * phys_jacobian * velocity_trial_space.tr() * meas(geometry_map);

  //gsInfo<<solution_vector.norm()<<"\n";

   //original trilinear term:
   auto trilin_conv_1 = velocity_sol_expression.tr() * phys_jacobian * velocity_trial_space.tr() * meas(geometry_map);

expression_assembler.initSystem();
   expression_assembler.assemble(bilin_conti, bilin_press, bilin_mu_1, trilin_conv_1);
  //gsInfo << expression_assembler.matrix().toDense()<<"\n";
  //gsInfo << expression_assembler.rhs().transpose()<<"\n";
  //gsInfo << solution_vector.transpose()<<" \n";

   solver.compute(expression_assembler.matrix());
   solution_vector = solver.solve(expression_assembler.rhs());
    gsInfo<<solution_vector.norm()<<"\n";
  
   gsInfo << "\n velocity" << l << "\n" << ev.eval(velocity_sol_expression, point) << std::endl;
   velocity_sol_expression.printDetail(gsInfo);

   gsInfo << "\n bilinear form viscosity\n"; bilin_mu_1.printDetail(gsInfo);
   gsInfo << "\n bilinear form pressure\n"; bilin_press.printDetail(gsInfo);
   gsInfo << "\n bilinear form continuity\n"; bilin_conti.printDetail(gsInfo);
   
   gsInfo << "\n trilinear form " << l << "\n";
   trilin_conv_1.printDetail(gsInfo);

   gsInfo << "\n trilinear form " << l << "\n" << ev.eval(trilin_conv_1, point).norm() << std::endl;
   trilin_conv_1.printDetail(gsInfo);

   gsInfo << "Iteration " << l << ": norm = " <<(solution_vector-initial_solution).norm() << std::endl;

    if ((solution_vector-initial_solution).norm()<1e-12) break;
  }

  nonlinear_solving_time_ls += timer.stop();
  gsInfo << "\tNonlinear solving done" << std::endl;
  timer.restart();

  // Export and visualization
  gsExprEvaluator<> expression_evaluator(expression_assembler);
/*
  gsInfo << expression_evaluator.eval(velocity_sol_expression.tr(), point) << '\n';
  gsInfo << expression_evaluator.eval(phys_jacobian, point) << '\n';
  gsInfo << expression_evaluator.eval(velocity_trial_space.tr(), point) << '\n';

  gsInfo << "\n product" <<"\n";
   gsInfo << expression_evaluator.eval(velocity_sol_expression.tr()*phys_jacobian, point) << '\n';

    gsInfo << "\n final_product" <<"\n";
   gsInfo << expression_evaluator.eval(velocity_sol_expression.tr()*phys_jacobian* velocity_trial_space.tr(), point) << '\n';

    gsInfo << (velocity_sol_expression.tr()*phys_jacobian* velocity_trial_space.tr()).Space << '\n';
*/

  gsInfo << "\nStarting the paraview export ..." << std::flush;
  //timer.restart();

  gsParaviewCollection collection("ParaviewOutput/solution_NS",
                                    &expression_evaluator);
  collection.options().setSwitch("plotElements", true);
  collection.options().setInt("plotElements.resolution", 16);
  collection.newTimeStep(&multi_patch);
  collection.addField(pressure_sol_expression, "pressure_NS");
  collection.addField(velocity_sol_expression, "velocity_NS");
  collection.saveTimeStep();
  collection.save();

  plotting_time += timer.stop();
  gsInfo << "\tFinished" << std::endl;
  
  gsInfo << "Starting the xml export ..." << std::flush;

  // Export pressure
  //gsMatrix<> solution_pressure;
  //gsFileData<> output_pressure;
    //pressure_sol_expression.extractFull(
      //  solution_pressure);  // patch-wise solution with BCs
    //output_pressure << solution_pressure;
   // output_pressure.save("pressure_field.xml");

    // Export velocity
    //gsMatrix<> solution_velocity;
    //gsFileData<> output;
    //velocity_sol_expression.extractFull(
     //   solution_velocity);  // patch-wise solution with BCs
    //output << solution_velocity;
    //output.save("velocity_field.xml");

// User output infor timings
  gsInfo << "\n\nTotal time: "
         << setup_time + assembly_time_ls + linear_solving_time_ls + plotting_time
         << std::endl;
  gsInfo << "                       Setup: " << setup_time << std::endl;
  gsInfo << "      Assembly Linear System: " << assembly_time_ls << std::endl;
  gsInfo << "       Solving Linear System: " << linear_solving_time_ls << std::endl;
  gsInfo << "                    Plotting: " << plotting_time << std::endl
         << std::flush;
        

    return EXIT_SUCCESS;

}// end main
