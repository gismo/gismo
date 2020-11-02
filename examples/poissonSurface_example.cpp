/** @file poisson_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
# include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 2;

    // Number for p-refinement of the computational (trial/test) basis.
    int degree     = 2;

    // Geometry
    index_t geometry = 0;

    //! [Parse command line]
    bool plot = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt("r", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "degree", "Polynomial degree", degree);
    cmd.addInt("g", "geometry", "Geometry type", geometry);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    index_t ddim = 2;
    if (geometry == 1)
        ddim = 3;

    //! [Function data]
    // Define source function
    //gsFunctionExpr<> f("0",3);
    gsFunctionExpr<> f("16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",ddim);
    //gsFunctionExpr<> f("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",ddim);
    // For homogeneous term, we can use this (last argument is the dimension of the domain)
    //gsConstantFunction<> f(0.0, 0.0, 2);

    // Define exact solution (optional)
    gsFunctionExpr<> g("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",ddim);
    //gsFunctionExpr<> g("1",3);

    // Print out source function and solution
    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<"\n\n";
    //! [Function data]

    //! [Geometry data]
    // Define Geometry, must be a gsMultiPatch object
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
    //! [Geometry data]

    // For single patch unit square of quadratic elements use (Note:
    // you need to update the bounadry conditions section for this to
    // work properly!) :
    // patches = gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquare(2));

    // Geometry can also be read from file (if gsMultiPatch):
    // gsReadFile<>("planar/lshape_p2.xml", patches);

    // Define Boundary conditions. Note that if one boundary is
    // "free", eg. if no condition is defined, then it is a natural
    // boundary (zero Neumann condition)
    // Also, remember that a pure Neumann problem has no unique
    // solution, thereforer implies a singular matrix. In this case
    // a corner DoF can be fixed to a given value to obtain a unique solution.
    // (example: bcInfo.addCornerValue(boundary::southwest, value, patch);)

    //! [Boundary conditions]
    gsBoundaryConditions<> bcInfo;
    //Alternatively: You can automatically create Dirichlet boundary
    //conditions using one function (the exact solution) for all
    //boundaries like this:

    for (gsMultiPatch<>::const_biterator
             bit = geo.bBegin(); bit != geo.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }
    //! [Boundary conditions]


    //! [Refinement]
    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( geo );

    // k-refinement (set degree)
    for ( size_t i = 0; i < refine_bases.nBases(); ++ i )
        refine_bases[i].setDegreePreservingMultiplicity(degree);

    //refine_bases.uniformRefine();
    gsInfo << "Degree of the geometry: " << refine_bases.maxCwiseDegree() << "\n";

    //! [Refinement]
    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
             "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        refine_bases.uniformRefine();

        //! [Assemble]
        gsPoissonAssembler<real_t> assembler(geo, refine_bases, bcInfo, f,
            //dirichlet::elimination, iFace::glue);
                                             dirichlet::elimination, iFace::glue);

        // Generate system matrix and load vector
        gsInfo << "Assembling...\n";
        assembler.assemble();
        gsInfo << "Have assembled a system (matrix and load vector) with "
               << assembler.numDofs() << " dofs.\n";
        //! [Assemble]

        //! [Solve]
        // Initialize the conjugate gradient solver
        gsInfo << "Solving...\n";
        gsSparseSolver<>::CGDiagonal solver(assembler.matrix());
        gsMatrix<> solVector = solver.solve(assembler.rhs());
        gsInfo << "Solved the system with CG solver.\n";
        //! [Solve]

        //! [Construct solution]
        // Construct the solution as a scalar field
        gsMultiPatch<> mpsol;
        assembler.constructSolution(solVector, mpsol);
        gsField<> sol(assembler.patches(), mpsol);
        //! [Construct solution]

        l2err[r] = sol.distanceL2(g, false);
        h1err.setOnes();

        if (plot && r == numRefine)
        {
            //! [Plot in Paraview]
            // Write approximate and exact solution to paraview files
            gsInfo << "Plotting in Paraview...\n";
            gsWriteParaview<>(sol, "poisson2d", 1000);
            const gsField<> exact(assembler.patches(), g, false);
            gsWriteParaview<>(exact, "poisson2d_exact", 1000);
            //! [Plot in Paraview]
        }
    }

    //! [Error and convergence rates]
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (numRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << ( l2err.head(numRefine).array() /
                  l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

}// end main
