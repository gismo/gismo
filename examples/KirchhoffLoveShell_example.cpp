/** @file KirchhoffLoveShell_example.cpp

    @brief A Kirchhoff-Love example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

# include <gismo.h>
# include <gsAssembler/gsKirchhoffLoveShellAssembler.h>
# include <gsG1Basis/gsG1AuxiliaryMultiplePatches.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 5;
    index_t numDegree = 1;
    bool plot = false;

    gsCmdLine cmd("Example for solving the Kirchhoff-Love problem.");
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

    gsFileData<> fileSrc("KirchhoffLoveGeo/square_curvedInterface.xml");
    gsInfo << "Loaded file " << fileSrc.lastPath() << "\n";

    gsMultiPatch<> geo;
    gsInfo << "Geometry taken correctly \n";
    fileSrc.getId(2, geo);
    geo.computeTopology();
    gsInfo << "Geometry computed correctly\n";
    gsMultiBasis<> basis(geo);

//    p-refine to get equal polynomial degree s,t directions
//    basis.degreeElevate(1,0);

//    for (int i = 0; i < numDegree; ++i)
//        basis.degreeElevate();
//    for (int i = 0; i < numRefine; ++i)
//        basis.uniformRefine();


    for (const boundaryInterface &  item : geo.interfaces() )
    {

        gsG1AuxiliaryMultiplePatches a(geo, item.first().patch, item.second().patch);
        gsMultiPatch<> test(a.reparametrizeG1Interface());
        gsMultiBasis<> testb(test);
        gsInfo << "New: " << testb << "\n";
    }




//    gsWriteParaview(newgeom1, "Geometry", 1000);

//    // Write file .xml of the new geometry
//    gsFileData<> fd;
//    fd << geo;
//    // output is a string. The extention .xml is added automatically
//    fd.save("newGeo");


    //Setting up oundary conditions
    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfo2;
    for (gsMultiPatch<>::const_biterator
             bit = geo.bBegin(); bit != geo.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &solution );
        bcInfo2.addCondition( *bit,  condition_type::neumann, &laplace);
    }




    //Initilize solver
    gsKirchhoffLoveShellAssembler<real_t> KirchhoffLoveShellAssembler( geo,basis,bcInfo,bcInfo2,source,
                                                       dirStrategy, intStrategy);

    gsInfo<<"Assembling..." << "\n";
    KirchhoffLoveShellAssembler.assemble();

    gsInfo<<"Solving with direct solver, "<< KirchhoffLoveShellAssembler.numDofs()<< " DoFs..."<< "\n";
    gsSparseSolver<real_t>::LU solver;
    solver.analyzePattern(KirchhoffLoveShellAssembler.matrix() );
    solver.factorize(KirchhoffLoveShellAssembler.matrix());
    gsMatrix<> solVector= solver.solve(KirchhoffLoveShellAssembler.rhs());

    //Reconstruct solution
    gsMultiPatch<> mpsol;
    KirchhoffLoveShellAssembler.constructSolution(solVector, mpsol);
    gsField<> solField(KirchhoffLoveShellAssembler.patches(), mpsol);

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
        gsWriteParaview<>(solField, "KirchhoffLoveShell2d", 5000);
        const gsField<> exact( geo, solution, false );
        gsWriteParaview<>( exact, "KirchhoffLoveShell2d_exact", 5000);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return  0;
}
