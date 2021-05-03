/** @file biharmonic_test.cpp

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
    index_t degree = 2;
    std::string kstr = "1";
    bool plot = false;
    bool L2ProjBC = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("r", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "degree", "Spline degree", degree);
    cmd.addString("k", "oscillation", "oscillation degree of manufactured solution", kstr);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "L2BC", "Use L2 projection to calculate Dirichlet BC instead of interpolation", L2ProjBC );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    
    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue; //Not used as we have single patch

    gsFunctionExpr<> solVal("sin("+kstr+"*pi*x) *exp(-"+kstr+"*pi*y)",2);

    gsFunctionExpr<> laplace ("0.0",2);
    
    
    gsFunctionExpr<> source  ("0.0",2);
    
    gsFunctionExpr<>sol1der (kstr+"*pi*cos("+kstr+"*pi*x) *exp(-"+kstr+"*pi*y)",
                             "-"+kstr+"*pi*sin("+kstr+"*pi*x) *exp(-"+kstr+"*pi*y)",2);
    
    gsFunctionExpr<>sol2der ("-"+kstr+"*"+kstr+"*pi*pi*sin("+kstr+"*pi*x) *exp(-"+kstr+"*pi*y)",
                             kstr+"*"+kstr+"*pi*pi*sin("+kstr+"*pi*x) *exp(-"+kstr+"*pi*y)",
                             "-"+kstr+"*"+kstr+"*pi*pi*cos("+kstr+"*pi*x) *exp(-"+kstr+"*pi*y)",2);

    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);


    
    gsMultiPatch<> geo = gsNurbsCreator<>::BSplineSquareGrid(1, 1, 1.0);
    gsMultiBasis<> basis(geo);
    
    
    
    for ( size_t i = 0; i < basis.nBases(); ++ i )
        basis[i].setDegreePreservingMultiplicity(degree);
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();


    //Setting up oundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &solution);
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &solution);
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
    if(L2ProjBC)
        BiharmonicAssembler.options().setInt("DirichletValues", dirichlet::l2Projection);
    
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
    const double PI = 3.14159265358979323846;
    real_t omega= std::stod(kstr)*PI;
    real_t normL2 = math::sqrt((1-math::exp(-2*omega))/(4*omega));
    gsInfo << "The L2 norm of the solution is : " << normL2 << "\n";
    gsInfo << "The relative L2 error of the solution is : " << errorL2/normL2 << "\n";
    
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
