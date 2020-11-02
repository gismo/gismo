/** @file biharmonic_example.cpp

    @brief A Biharmonic example with different possible boundary conditions on a single Patch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn & P. Weinmüller
*/

# include <gismo.h>
# include <gsAssembler/gsBiharmonicAssembler.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 5;
    index_t numDegree = 1;

    index_t geometry = 0;

    bool plot = false;
    bool neumann = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("r", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "degree", "Polynomial degree", numDegree);
    cmd.addInt("g", "geometry", "Geometry type", geometry);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "neumann", "Neumann boundary", neumann);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue;

/*
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<> sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                              "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<> sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                              "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                              " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)",2);


*/
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",3);
    gsFunctionExpr<> sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                              "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",
                              "0",3);
    gsFunctionExpr<> sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                              "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                              " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)",
                              "0",
                              "0",
                              "0",
                              3);
/*
    //gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*(x/cos(45)))*cos(4*pi*y) - cos(4*pi*(x/cos(45))) - cos(4*pi*y))",3);
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(-cos(4*pi*x/cos(45))*1/(cos(45)*cos(45)*cos(45)*cos(45)) + cos(4*pi*y) * "
                              "(-1 + cos(4*pi*x/cos(45)) * (1+1/(cos(45)*cos(45))) * (1+1/(cos(45)*cos(45))) ))",3);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    gsFunctionExpr<> solVal("(cos(4*pi*(x/cos(45))) - 1) * (cos(4*pi*y) - 1) ",3);
    gsFunctionExpr<> sol1der ("-4*pi*1/cos(45) * (cos(4*pi*y) - 1)*sin(4*pi*(x/cos(45)))",
                              "-4*pi*(cos(4*pi*(x/cos(45))) - 1)*sin(4*pi*y)",
                              "0",3);
    gsFunctionExpr<> sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                              "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                              " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)",
                              "0",
                              "0",
                              "0",
                              3);

    gsFunctionExpr<> source  ("0",3);
    gsFunctionExpr<> laplace ("0",3);
    gsFunctionExpr<> solVal("(x/cos(45)-1)*(y-1)",3);
    gsFunctionExpr<> sol1der ("y-1",
                              "x/cos(45)-1",
                              "0",3);
    gsFunctionExpr<> sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                              "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                              " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)",
                              "0",
                              "0",
                              "0",
                              3);

    gsFunctionExpr<> source  ("0",3);
    gsFunctionExpr<> laplace ("0",3);
    gsFunctionExpr<> solVal("z",3);
    gsFunctionExpr<>sol1der ("0",
                            "0",
                             "1",3);
    gsFunctionExpr<>sol2der ("0",
                             "0",
                             "0",
                             "0",
                             "0",
                             "0", 3);
    */

    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    //gsMultiPatch<> geo( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //gsMultiPatch<> geo( *gsNurbsCreator<>::BSplineSquare(1,0.25,0.25) );

    // ======= Geometry =========
    std::string string_geo;
    switch(geometry)
    {
        case 0: // planar
            string_geo = "KirchhoffLoveGeo/square_singlePatch.xml";
            break;
        case 1: // surface
            string_geo = "KirchhoffLoveGeo/squareSurface3d.xml";
            break;

        default:
            gsInfo << "No geometry is used! \n";
            break;
    }

    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> geo;
    fd.getId(0, geo); // id=0: Multipatch domain
    geo.computeTopology();


    gsMultiBasis<> basis(geo);

    //p-refine to get equal polynomial degree s,t directions (for Annulus)
    //basis.degreeElevate(1,0);
    for (int i = 0; i < numDegree; ++i)
        basis.degreeElevate();
    //for (int i = 0; i < numRefine; ++i)
    //    basis.uniformRefine();

    basis.uniformRefine();

    gsInfo << "Degree: " << basis.maxCwiseDegree() << "\n";

    //Setting up boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &solVal);
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &solVal);
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &solVal);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &solVal);

    // Second boundary condition
    gsBoundaryConditions<> bcInfo2;
    if (!neumann)
    {
        //Laplace condition of second kind
        bcInfo2.addCondition( boundary::west,  condition_type::laplace, &laplace);
        bcInfo2.addCondition( boundary::east,  condition_type::laplace, &laplace);
        bcInfo2.addCondition( boundary::north, condition_type::laplace, &laplace);
        bcInfo2.addCondition( boundary::south, condition_type::laplace, &laplace);
    }
    else
    {
        //Neumann condition of second kind
        bcInfo2.addCondition( boundary::west,  condition_type::neumann, &sol1der);
        bcInfo2.addCondition( boundary::east,  condition_type::neumann, &sol1der);
        bcInfo2.addCondition( boundary::north, condition_type::neumann, &sol1der);
        bcInfo2.addCondition( boundary::south, condition_type::neumann, &sol1der);
    }



    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
             "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        basis.uniformRefine();

        //Initilize solver
        gsBiharmonicAssembler<real_t> BiharmonicAssembler(geo, basis, bcInfo, bcInfo2, source,
                                                          dirStrategy, intStrategy);

        //gsInfo << "Assembling..." << "\n";
        gsInfo<< BiharmonicAssembler.numDofs() <<std::flush;

        BiharmonicAssembler.assemble();

        gsInfo<< "." <<std::flush;// Assemblying done

        //gsInfo << "Solving with direct solver, " << BiharmonicAssembler.numDofs() << " DoFs..." << "\n";
//        gsSparseSolver<real_t>::LU solver;
//        solver.analyzePattern(BiharmonicAssembler.matrix());
//        solver.factorize(BiharmonicAssembler.matrix());
//        gsMatrix<> solVector = solver.solve(BiharmonicAssembler.rhs());



        gsSparseSolver<real_t>::CGDiagonal solver;
        solver.analyzePattern(BiharmonicAssembler.matrix());
        solver.factorize(BiharmonicAssembler.matrix());
        gsMatrix<> solVector = solver.solve(BiharmonicAssembler.rhs());

        gsInfo << "rhs: " << BiharmonicAssembler.rhs() << "\n";
        gsInfo << "Matrix: " << BiharmonicAssembler.matrix().toDense() << "\n";
        gsInfo << "sol: " << solVector << "\n";

        gsInfo<< "." <<std::flush; // Linear solving done


        //Reconstruct solution
        gsMultiPatch<> mpsol;
        BiharmonicAssembler.constructSolution(solVector, mpsol);
        gsField<> solField(BiharmonicAssembler.patches(), mpsol);

        //Contruct the H2 norm, part by part.
        real_t errorH2Semi, errorH1Semi;

#pragma omp parallel for
        for (index_t e = 0; e < 3; ++e)
        {
            if(geo.dimensions().first+1 == geo.dimensions().second) // TODO
            {
                if (e == 2)
                    l2err[r] = solField.distanceL2(solution, false);

                h1err.setOnes();
                h2err.setOnes();
            }
            else
            {
                if (e == 0)
                    errorH2Semi = solField.distanceH2(solution, false);
                else if (e == 1)
                    errorH1Semi = solField.distanceH1(solution, false);
                else if (e == 2)
                    l2err[r] = solField.distanceL2(solution, false);
            }
        }

        if(geo.dimensions().first == geo.dimensions().second)
        {
            h1err[r] = math::sqrt(errorH1Semi * errorH1Semi + l2err[r] * l2err[r]);
            h2err[r] = math::sqrt(errorH2Semi * errorH2Semi + errorH1Semi * errorH1Semi + l2err[r] * l2err[r]);

        }

        gsInfo<< ". " <<std::flush; // Error computations done

        // Plot solution in paraview
        if (plot && r == numRefine)
        {
            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in ParaView...\n";
            gsWriteParaview<>(solField, "Biharmonic2d", 5000);
            const gsField<> exact( geo, solVal, false );
            gsWriteParaview<>( exact, "Biharmonic2d_exact", 5000);
        }
        else if (r == numRefine)
            gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                      "file containing the solution.\n";
    }

    //! [Error and convergence rates]
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<h2err.transpose()<<"\n";

    if (numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << ( l2err.head(numRefine).array() /
                  l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              <<( h2err.head(numRefine).array() /
                  h2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]



    return  0;
}
