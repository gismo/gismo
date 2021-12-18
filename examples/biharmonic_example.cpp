/** @file biharmonic_example.cpp

    @brief A Biharmonic example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

# include <gismo.h>
# include <gsAssembler/gsBiharmonicAssembler.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 5;
    index_t numDegree = 1;
    bool plot = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("r", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "degree", "Polynomial degree", numDegree);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue;

    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                              "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                              "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                              " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    gsMultiPatch<> geo( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    gsMultiBasis<> basis(geo);

    //p-refine to get equal polynomial degree s,t directions (for Annulus)
    basis.degreeElevate(1,0);

    for (int i = 0; i < numDegree; ++i)
        basis.degreeElevate();
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();


    //Setting up oundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &solution);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &solution);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &solution);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &solution);
    //Neumann condition of second kind
    gsBoundaryConditions<> bcInfo2;
    bcInfo2.addCondition( boundary::west,  condition_type::neumann, &laplace);
    bcInfo2.addCondition( boundary::east,  condition_type::neumann, &laplace);
    bcInfo2.addCondition( boundary::north, condition_type::neumann, &laplace);
    bcInfo2.addCondition( boundary::south, condition_type::neumann, &laplace);

    //Initilize solver
    gsBiharmonicAssembler<real_t> BiharmonicAssembler( geo,basis,bcInfo,bcInfo2,source,
                                                       dirStrategy, intStrategy);

    gsInfo<<"Assembling..." << "\n";
    BiharmonicAssembler.assemble();

    gsInfo<<"Solving with direct solver, "<< BiharmonicAssembler.numDofs()<< " DoFs..."<< "\n";
    gsSparseSolver<real_t>::LU solver;
    solver.analyzePattern(BiharmonicAssembler.matrix() );
    solver.factorize(BiharmonicAssembler.matrix());
    gsMatrix<> solVector= solver.solve(BiharmonicAssembler.rhs());

    //Reconstruct solution
    gsMultiPatch<> mpsol;
    BiharmonicAssembler.constructSolution(solVector, mpsol);
    gsField<> solField(BiharmonicAssembler.patches(), mpsol);

    //Contruct the H2 norm, part by part.
    real_t errorH2Semi = solField.distanceH2(solution, false);
    real_t errorH1Semi = solField.distanceH1(solution, false);
    real_t errorL2 = solField.distanceL2(solution, false);
    real_t errorH1 = math::sqrt(errorH1Semi*errorH1Semi + errorL2*errorL2);
    real_t errorH2 = math::sqrt(errorH2Semi*errorH2Semi + errorH1Semi*errorH1Semi + errorL2*errorL2);

    gsInfo << "The L2 error of the solution is : " << errorL2 << "\n";
    gsInfo << "The H1 error of the solution is : " << errorH1 << "\n";
    gsInfo << "The H2 error of the solution is : " << errorH2 << "\n";

    // Plot solution in paraview
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in ParaView...\n";
        gsWriteParaview<>(solField, "Biharmonic2d", 5000);
        const gsField<> exact( geo, solution, false );
        gsWriteParaview<>( exact, "Biharmonic2d_exact", 5000);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return  0;
}
