/** @file precision_example.cpp

    @brief Example illustrating the effect of quantization and rounding errors in the absence of sufficient precision.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): N. Kohl
*/

#include <gismo.h>

using namespace gismo;

struct DataIn
{
    index_t refine = 5;
    index_t degree = 5;
    bool plot = false;
};

struct DataOut
{
    real_t errorL2;
    real_t errorH1;

    int dofs;
    int solverInfo;
};

///
/// The example showcases the effect of quantization and rounding errors. In a nutshell: the required precision to
/// ensure discretization-error-accurate solutions grows with the condition number of the system matrix.
///
/// In practice, for a fixed problem, the precision requirements therefore grow with the refinement level. The higher
/// the order of the underlying discretization, the faster the precision requirements grow. This test sets up a simple
/// grid study for a Poisson problem on a square domain to demonstrate the issue.
///
/// To play around with the precision, you can build this example with the MPFR library (supported directly through
/// G+SMO: just set
///
///     -DGISMO_COEFF_TYPE=mpfr::mpreal
///
/// via cmake).
///
/// See e.g.,
///
///     Tamstorf, R., Benzaken, J., & McCormick, S. F. (2021).
///     Discretization-error-accurate mixed-precision multigrid solvers.
///     SIAM Journal on Scientific Computing, 43(5), S420-S447.
///
/// for details.
///
DataOut test(DataIn dataIn)
{
    DataOut dataOut;

    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue;

    gsFunctionExpr<> source("0", 2);
    gsFunctionExpr<> solVal("sinh(x) * sin(y)", 2);
    gsFunctionExpr<> sol1der("cosh(x) * sin(y)", "sinh(x) * cos(y)", 2);

    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der);

    gsMultiPatch<> geo(*gsNurbsCreator<>::BSplineSquare());
    gsMultiBasis<> basis(geo);

    for (int i = 0; i < dataIn.degree - 1; ++i)
        basis.degreeElevate();
    for (int i = 0; i < dataIn.refine; ++i)
        basis.uniformRefine();

    // Setting up boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(boundary::west, condition_type::dirichlet, &solution);
    bcInfo.addCondition(boundary::east, condition_type::dirichlet, &solution);
    bcInfo.addCondition(boundary::north, condition_type::dirichlet, &solution);
    bcInfo.addCondition(boundary::south, condition_type::dirichlet, &solution);

    // Assemble system.
    gsPoissonAssembler<real_t> poissonAssembler(geo, basis, bcInfo, source, dirStrategy, intStrategy);
    poissonAssembler.assemble();
    dataOut.dofs = poissonAssembler.numDofs();

    // Setup linear solver.
    gsSparseSolver<real_t>::CGDiagonal solver;
    // Set iterations to more than default.
    solver.setMaxIterations(dataOut.dofs * dataOut.dofs);
    solver.compute(poissonAssembler.matrix());
    gsMatrix<> solVector = solver.solve(poissonAssembler.rhs());
    // Warning the user later if the solver did not converge.
    dataOut.solverInfo = solver.info();

    // Reconstruct solution.
    gsMultiPatch<> mpsol;
    poissonAssembler.constructSolution(solVector, mpsol);
    gsField<> solField(poissonAssembler.patches(), mpsol);

    real_t errorH1Semi = solField.distanceH1(solution, false);
    real_t errorL2 = solField.distanceL2(solution, false);
    real_t errorH1 = math::sqrt(errorH1Semi * errorH1Semi + errorL2 * errorL2);

    dataOut.errorL2 = errorL2;
    dataOut.errorH1 = errorH1;

    // Plot solution in paraview.
    if (dataIn.plot)
    {
        auto baseName = "poisson_r" + std::to_string(dataIn.refine) + "_p" + std::to_string(dataIn.degree);
        gsWriteParaview<>(solField, baseName + "_computed", 5000);
        const gsField<> exact(geo, solution, false);
        gsWriteParaview<>(exact, baseName + "_exact", 5000);
    }

    return dataOut;
}

int main(int argc, char *argv[])
{
    DataIn dataIn;

    index_t defaultDegree = 5;
    dataIn.degree = defaultDegree;
    index_t maxRefine = 7;

    // Setting default to IEEE 754 "double"
    index_t mantissa = 53;

    gsCmdLine cmd("Example for solving the Poisson problem.");
    cmd.addInt("r", "refine", "Number of refinement steps", maxRefine);
    cmd.addInt("p", "degree", "Polynomial degree", dataIn.degree);
    cmd.addInt("m", "mantissa", "Mantissa bits", mantissa);
    cmd.addSwitch("plot", "Plot result in ParaView format", dataIn.plot);

    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int rv)
    {
        return rv;
    }

#ifdef gsMpfr_ENABLED
    mpfr::mpreal::set_default_prec(mantissa);
#else
    gsInfo << "Since MPFR is not enabled, the mantissa width cannot be selected as a parameter. Defaulting to what is "
              "set via real_t.\n";
    mantissa = std::numeric_limits<real_t>::digits;
#endif

    gsInfo << "Parameters:\n";
    gsInfo << " - mantissa bits:  " << mantissa << "\n";
    gsInfo << " - max refinement: " << maxRefine << "\n";
    gsInfo << " - degree:         " << dataIn.degree << "\n";

    gsInfo << "\n";

    real_t expectedRateL2 = math::pow(2, -(dataIn.degree + 1));
    real_t expectedRateH1 = math::pow(2, -dataIn.degree);

    gsInfo << "Expected rates: \n";
    gsInfo << " - L2: " << std::fixed << expectedRateL2 << "\n";
    gsInfo << " - H1: " << std::fixed << expectedRateH1 << "\n";

    gsInfo << "\n";

    gsInfo << "Starting grid study. Observe the rates and compare to the expected values.\n";
    gsInfo << "If the floating point precision is insufficient, they will blow up (and get _worse_ with more "
              "refinement). This should happen with the default values.\n";
    gsInfo
        << "The tipping point can be reached 'quicker' with higher order discretizations or a lower precision floating "
           "point format.\n";
    gsInfo << "\n";

    gsInfo << "refinements || L2 error     | rate     || H1 error     | rate     \n";

    real_t L2prev = 0;
    real_t H1prev = 0;

    for (int r = 1; r <= maxRefine; r++)
    {
        dataIn.refine = r;

        auto dataOut = test(dataIn);

        real_t L2rate = 0;
        real_t H1rate = 0;

        if (r > 1)
        {
            L2rate = dataOut.errorL2 / L2prev;
            H1rate = dataOut.errorH1 / H1prev;
        }

        gsInfo << std::setw(11) << r << " || " << std::scientific << dataOut.errorL2 << " | " << std::fixed << L2rate
               << " || " << std::scientific << dataOut.errorH1 << " | " << std::fixed << H1rate;

        if (L2rate > 2 * expectedRateL2 || H1rate > 2 * expectedRateH1)
        {
            gsInfo << " || at least one rate > 2 * expected\n";
        }
        else
        {
            gsInfo << "\n";
        }

        L2prev = dataOut.errorL2;
        H1prev = dataOut.errorH1;

        if (dataOut.solverInfo != gsEigen::ComputationInfo::Success)
        {
            gsInfo << "Solver did not converge...\n";
        }
    }
}
