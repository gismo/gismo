/** @file linear_elasticity_example.cpp

  @brief Expression assembler to solve the a linear elasticity problem

  Based on the poisson_example

  The test case is based on the proposed manufactured solution found here:
  https://web.mit.edu/16.20/homepage/4_ElasticityBVP/ElasticityBVP_files/module_4_with_solutions.pdf
  page: 82

  The solution must be degree 2 minimum to produce accurate results, use
  degree elevation and refinement to modify accuracy

  This file is part of the G+Smo library.

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[]) {
  //! [Parse command line]
  bool plot = false;
  index_t n_h_refinements = 0;
  index_t n_deg_elevations = 0;
  bool only_last = false;
  bool compute_error{false};
  index_t sample_rate{9};

  std::string file_name("pde/linear_elasticity_example.xml");

  gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
  cmd.addInt("e", "degreeElevation",
             "Number of degree elevation steps to perform before solving (0: "
             "equalize degree in all directions)",
             n_deg_elevations);
  cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement loops",
             n_h_refinements);
  cmd.addString("f", "file", "Input XML file", file_name);
  cmd.addSwitch("only_last",
                "Solve solely for the only_last level of h-refinement",
                only_last);
  cmd.addSwitch("plot",
                "Create a ParaView visualization file with the solution", plot);
  cmd.addInt("q", "sample-rate", "Samples per spline in paraview export",
             sample_rate);
  cmd.addSwitch("compute-error",
                "Evaluate the error with respect to the analytical solution "
                "(evaluation with default options and default file required)",
                compute_error);

  // Error is computed for an analytical solution with the following
  // Lame constants
  real_t lame_lambda{80000.0}, lame_mu{80000.0};
  cmd.addReal("L", "firstLame", "First Lame constant, material parameter",
              lame_lambda);
  cmd.addReal("M", "secondLame", "Second Lame constant, material parameter",
              lame_mu);

  //! [Parse command line]
  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  //! [Read input file]
  gsFileData<> file_data(file_name);
  gsInfo << "Loaded file " << file_data.lastPath() << "\n";

  gsMultiPatch<> multi_patch;
  file_data.getId(0, multi_patch);  // id=0: Multipatch domain

  gsFunctionExpr<> source_function_expression;
  file_data.getId(1, source_function_expression);  // id=1: source function
  gsInfo << "Source function " << source_function_expression << "\n";

  gsBoundaryConditions<> boundary_conditions;
  file_data.getId(2, boundary_conditions);  // id=2: boundary conditions
  boundary_conditions.setGeoMap(multi_patch);
  gsInfo << "Boundary conditions:\n" << boundary_conditions << "\n";

  gsFunctionExpr<> reference_solution_expr;
  file_data.getId(3, reference_solution_expr);  // id=3: reference solution

  gsOptionList assembler_options = gsAssembler<>::defaultOptions();
  if (file_data.hasId(4)) {
    file_data.getId(4, assembler_options);  // id=4: assembler options
  }
  //! [Read input file]

  //! [Refinement]
  gsMultiBasis<> function_basis(multi_patch, true);
  const int solution_field_dimension = multi_patch.geoDim();

  // Elevate and p-refine the basis
  function_basis.degreeElevate(n_deg_elevations);

  // h-refine each basis
  for (int r = 0; r < n_h_refinements; ++r) {
    function_basis.uniformRefine();
  }
  //! [Refinement]

  //! [Problem setup]
  gsExprAssembler<> expression_assembler(1, 1);
  if (file_data.hasId(4)) {
    expression_assembler.setOptions(assembler_options);
  }

  typedef gsExprAssembler<>::geometryMap GeometryMap;
  typedef gsExprAssembler<>::variable Variable;
  typedef gsExprAssembler<>::space Space;
  typedef gsExprAssembler<>::solution Solution;

  // Elements used for numerical integration
  expression_assembler.setIntegrationElements(function_basis);
  gsExprEvaluator<> ev(expression_assembler);

  // Set the geometry map
  GeometryMap geometry_map = expression_assembler.getMap(multi_patch);

  // Set the discretization space
  Space u_trial =
      expression_assembler.getSpace(function_basis, solution_field_dimension);

  // Set the source term
  auto source_function =
      expression_assembler.getCoeff(source_function_expression, geometry_map);

  // Recover manufactured solution
  auto reference_solution =
      ev.getVariable(reference_solution_expr, geometry_map);

  // Solution vector and solution variable
  gsMatrix<> solution_vector;
  Solution u_solution_expression =
      expression_assembler.getSolution(u_trial, solution_vector);
  //! [Problem setup]

  // ![User output]
  gsInfo << "Problem Overview:\t"
         << "\nGeometric Dim:\t" << solution_field_dimension << "\nPatches:\t"
         << multi_patch.nPatches() << "\nMinimum degree:\t"
         << function_basis.minCwiseDegree() << "\nMaximum degree:\t"
         << function_basis.maxCwiseDegree() << "\n\n";
#ifdef _OPENMP
  gsInfo << "Available threads: " << omp_get_max_threads() << "\n";
#endif
  gsInfo << "Active assembly options:\n"
         << expression_assembler.options() << std::endl;
  // ![User output]

  //! [Assembly]
  gsSparseSolver<>::CGDiagonal linear_solver;

  real_t l2_err{}, h1_err{}, setup_time{}, ma_time{}, slv_time{}, err_time{};
  gsStopwatch timer;
  function_basis.uniformRefine();

  // Apply dirichlet boundary conditions
  u_trial.setup(boundary_conditions, dirichlet::l2Projection, 0);

  // Initialize the system
  expression_assembler.initSystem();
  setup_time += timer.stop();

  // User output
  gsInfo << "Number of degrees of freedom:\t" << expression_assembler.numDofs()
         << std::endl;
  gsInfo << "\nAssembling linear system ... " << std::flush;
  timer.restart();

  // Compute the system matrix and right-hand side
  auto phys_jacobian = ijac(u_trial, geometry_map);

  // Bilinear form
  auto bilin_lambda = lame_lambda * idiv(u_trial, geometry_map) *
                      idiv(u_trial, geometry_map).tr() * meas(geometry_map);
  auto bilin_mu_ =
      lame_mu *
      ((phys_jacobian.cwisetr() + phys_jacobian) % phys_jacobian.tr()) *
      meas(geometry_map);
  auto bilin_combined = (bilin_lambda + bilin_mu_);

  // Linear Form
  auto linear_form = u_trial * source_function * meas(geometry_map);

  expression_assembler.assemble(bilin_combined,  // matrix
                                linear_form      // rhs vector
  );

  // Compute the Neumann terms (will not be used in default file)
  auto neumann_boundary_function =
      expression_assembler.getBdrFunction(geometry_map);

  // Neumann conditions
  expression_assembler.assembleBdr(
      boundary_conditions.get("Neumann"),
      u_trial * neumann_boundary_function * nv(geometry_map).norm());
  //! [Assembly]

  //! [Solver]
  ma_time += timer.stop();
  gsInfo << "Done" << std::endl;
  gsInfo << "Solving Linear System ...... " << std::flush;
  timer.restart();

  // Solve linear system
  linear_solver.compute(expression_assembler.matrix());
  solution_vector = linear_solver.solve(expression_assembler.rhs());
  //! [Solver]

  slv_time += timer.stop();
  gsInfo << "Done" << std::endl;

  //! [Evaluate Error]
  if (compute_error) {
    gsInfo << "Evaluating solution ........ " << std::flush;
    timer.restart();
    l2_err = math::sqrt(
        ev.integral((reference_solution - u_solution_expression).sqNorm() *
                    meas(geometry_map)));

    h1_err = l2_err +
             math::sqrt(ev.integral((igrad(reference_solution) -
                                     igrad(u_solution_expression, geometry_map))
                                        .sqNorm() *
                                    meas(geometry_map)));
    err_time += timer.stop();
    gsInfo << "Done" << std::endl;
  }
  //! [Evaluate Error]

  //! [Summarize User Output]
  gsInfo << "\n\nTotal time: " << setup_time + ma_time + slv_time + err_time
         << "\n";
  gsInfo << "     Setup: " << setup_time << "\n";
  gsInfo << "  Assembly: " << ma_time << "\n";
  gsInfo << "   Solving: " << slv_time << "\n";
  gsInfo << "     Norms: " << err_time << "\n";

  if (compute_error) {
    gsInfo << "\nL2 error: " << std::scientific << std::setprecision(3)
           << l2_err << "\nH1 error: " << std::scientific << h1_err << "\n\n";
  }
  //! [Summarize User Output]

  //! [Export visualization in ParaView]
  if (plot) {
    gsInfo << "Plotting in Paraview ... ";

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", sample_rate);
    collection.newTimeStep(&multi_patch);
    collection.addField(u_solution_expression, "numerical solution");
    if (compute_error) {
      collection.addField(reference_solution - u_solution_expression, "error");
    }
    collection.saveTimeStep();
    collection.save();

    gsFileManager::open("ParaviewOutput/solution.pvd");
    gsInfo << "Done" << std::endl;
  } else {
    gsInfo << "No output created, re-run with --plot to get a ParaView "
              "file containing the solution.\n";
  }
  //! [Export visualization in ParaView]

  return EXIT_SUCCESS;

}  // end main
