#/** @file biharmonic_argyris_example.cpp

    @brief Example using the Argyris space for solving biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <omp.h>


//! [Include namespace]
#include <gismo.h>

# include <gsAssembler/gsG1BiharmonicAssembler.h>
# include <gsArgyris/gsC1Argyris.h>
# include <gsArgyris/gsErrorAnalysis/gsArgyrisNorms.h>

# include <gsArgyris/gsC1ArgyrisIO.h>

# include <gsMSplines/gsMappedBasis.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool last       = false;
    index_t discrete_p = 0; // Degree of the geometry + discrete_p = polynomial degree of the discrete space
    index_t discrete_r = 1; // Regularity of the geometry + discrete_r = regularity of the discrete space

    index_t numRefine = 0; // Number of loops = refinement steps

    index_t geometry = 0;
    std::string input;

    bool isogeometric = false;
    bool exactGD = false;
    bool neumann = false;

    bool info = false;
    bool latex = false;

    bool twoPatch = false;

    gsCmdLine cmd("Solving biharmonic equation with Argyris space.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );

    // For the discrete space
    cmd.addInt( "p", "degreeElevate",
                "Number of degree elevation steps to perform before solving", discrete_p );
    cmd.addInt( "r", "regularity",
                "Number of increase Continuity steps to perform before solving",  discrete_r );

    // For computing the convergence rates
    cmd.addInt( "l", "loop",
                "Number of Uniform h-refinement steps to perform before solving",  numRefine );

    // Which method one want to use
    cmd.addSwitch( "exactGD", "To compute the gluing data exact", exactGD );
    cmd.addSwitch( "isogeometric", "Project the basis in isogeometric concept", isogeometric );
    cmd.addSwitch("neumann","Compute the biharmonic with neumann bdy",neumann);

    cmd.addSwitch("twoPatch","Two Patch",twoPatch);

    // Output features
    cmd.addSwitch("latex","Print the rate and error latex-ready",latex);
    cmd.addSwitch("plot", "Plotting the results!",plot);
    cmd.addSwitch("last", "Last case only!",last);
    cmd.addSwitch( "info", "Print information", info );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Exact solution]

    gsFunctionExpr<> source("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<> sol1der("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<> sol2der("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
/*
    gsFunctionExpr<> source("0",2);
    gsFunctionExpr<> laplace("0",2);
    gsFunctionExpr<> solVal("x",2);
    gsFunctionExpr<> sol1der("1",
                             "0",2);
    gsFunctionExpr<> sol2der("0",
                             "0",
                             "0", 2);
    */
    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    //! [Exact solution]

    //! [Problem setup]
    gsFileData<> fd;
    gsMultiPatch<> mp;
    gsMultiBasis<> mb;

    // For input/output stuff
    gsC1ArgyrisIO c1ArgyrisIO;

    gsOptionList optionList;
    optionList.addInt( "degreeElevate", "Number of degree elevation steps to perform before solving", discrete_p );
    optionList.addInt( "regularity", "Number of increase Continuity steps to perform before solving",  discrete_r );
    optionList.addSwitch("isogeometric", "Project the basis in isogeometric concept", isogeometric);
    optionList.addSwitch("exactGD", "To compute the gluing data exact", exactGD);
    optionList.addSwitch("info", "Print information", info);
    optionList.addSwitch("plot", "Plot", plot);

    optionList.addSwitch("twoPatch", "Two Patch", twoPatch);
    //! [Problem setup]

    //! [Read geometry]
    std::string string_geo;
    if (input.empty())
    {
        switch(geometry)
        {
            case 0:
                string_geo = "planar/twoPatches/two_squares_linear.xml";
                break;
            case 1:
                string_geo = "planar/twoPatches/two_squares_linear_with_inner_knot.xml";
                break;
            case 2:
                string_geo = "planar/twoPatches/square_cubic_with_inner_knot.xml";
                break;
            case 3:
                string_geo = "planar/twoPatches/two_squares_linear_partial_matching.xml";
                break;
            case 4:
                string_geo = "planar/twoPatches/two_squares_linear_non_matching.xml";
                break;
            case 5:
                string_geo = "planar/twoPatches/two_squares_linear_with_inner_knot2.xml";
                break;
            case 6:
                string_geo = "planar/twoPatches/square_cubic_with_inner_knot2.xml";
                break;

            case 10:
                string_geo = "planar/twoPatches/square_curved.xml";
                break;
            case 11:
                string_geo = "planar/twoPatches/funny_example.xml";
                break;

            case 100:
                string_geo = "planar/multiPatches/four_squares_linear.xml";
                break;

            case -8 ... -1:
                string_geo = "planar/twoPatches/benchmark/square_linear_" + std::to_string(-geometry) + ".xml";
                break;

            default:
                gsInfo << "No geometry is used! \n";
                break;
        }
    }

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();
    //! [Read geometry]

    //! [Check the input data]
    c1ArgyrisIO.checkInput(mp, optionList);
    //! [Check the input data]

    //! [Boundary condition]
    gsBoundaryConditions<> bcInfo, bcInfo2;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bcInfo.addCondition(*bit, condition_type::dirichlet, &solVal);
        if (!neumann)
            bcInfo2.addCondition(*bit, condition_type::laplace, &laplace);
        else
            bcInfo2.addCondition(*bit,  condition_type::neumann, &sol1der);
    }
    //! [Boundary condition]

    //! [Initialise the discrete space]
    // TODO
    mb = gsMultiBasis<>(mp);
    gsC1Argyris<2, real_t> c1Argyris(mp, optionList);
    gsMappedBasis<2,real_t> mappedBasis;
    /*
    std::vector<gsBasis<real_t> *> test = c1Argyris.getBases();
    gsMultiBasis<> mb_temp(test, mp.topology());
    gsInfo << "TEST: " << mb_temp.nBases() << "\n";
    gsInfo << "TEST 2: " << mb_temp[0].size() << "\n";
    */
    if (last)
    {
        // h-refine
        for (int l =0; l < numRefine; ++l)
             c1Argyris.uniformRefine();

        numRefine = 0;
    }
    //! [Initialise the discrete space]

    gsSparseSolver<>::CGDiagonal solver;

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), h2err(numRefine+1), dofs(numRefine+1), meshSize(numRefine+1);
    gsMatrix<> jumperr(numRefine+1, mp.nInterfaces());

    gsInfo<< "(dot1=got_argyris_space, dot2=assembled, dot3=solved, dot4=got_error)\n";
    for( index_t l = 0; l<=numRefine; ++l)
    {
        gsSparseMatrix<> sparseMatrix_argyris;
        gsMultiBasis<> mb_argyris;

        c1Argyris.init();
        c1Argyris.createArgyrisSpace(); // Slow TODO
        if (plot) {
            gsInfo << "Plot start \n";
            c1Argyris.writeParaviewSinglePatch(0, "inner");
            c1Argyris.writeParaviewSinglePatch(0, "edge");
            c1Argyris.writeParaviewSinglePatch(0, "vertex");

            c1Argyris.writeParaviewSinglePatch(1, "inner");
            c1Argyris.writeParaviewSinglePatch(1, "edge");
            c1Argyris.writeParaviewSinglePatch(1, "vertex");
            gsInfo << "Plot end \n";
        }
        c1Argyris.getMultiBasis(mb_argyris);
        sparseMatrix_argyris = c1Argyris.getSystem();
        mappedBasis.init(mb_argyris, sparseMatrix_argyris.transpose());
        gsInfo<< "." <<std::flush;// Construction of Argyris space done

        gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(mp, mappedBasis, bcInfo, bcInfo2, source, twoPatch);
        g1BiharmonicAssembler.assemble();
        gsInfo<< "." <<std::flush;// Assemblying done

        gsSparseSolver<real_t>::CGDiagonal solver;
        solver.compute(g1BiharmonicAssembler.matrix());
        gsMatrix<real_t> solVector= solver.solve(g1BiharmonicAssembler.rhs());

        gsMatrix<real_t> solFull;
        g1BiharmonicAssembler.constructSolution(solVector, solFull);
        sparseMatrix_argyris = solFull.asDiagonal() * sparseMatrix_argyris;
        c1Argyris.setSystem(sparseMatrix_argyris);
        gsInfo<< "." <<std::flush;// Linear solving done

        mappedBasis.init(mb_argyris, sparseMatrix_argyris.transpose());
        gsArgyrisNorms<real_t> argyrisNorms(mp, mappedBasis, solution);
        argyrisNorms.compute();

        // Collecting data
        meshSize[l] = c1Argyris.getMinMeshSize();
        dofs[l] = g1BiharmonicAssembler.numDofs();

        l2err[l] = argyrisNorms.valueL2();
        h1err[l] = math::sqrt(argyrisNorms.valueH1() * argyrisNorms.valueH1() + l2err[l] * l2err[l]);
        h2err[l] = math::sqrt(argyrisNorms.valueH2() * argyrisNorms.valueH2() +
                argyrisNorms.valueH1() * argyrisNorms.valueH1() + l2err[l] * l2err[l]);
        gsInfo<< ". " <<std::flush;// Error computations done

        // TODO Refine spaces
        c1Argyris.uniformRefine();
    }
    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo << "\nDofs: " << dofs.transpose() << "\n";
    gsInfo << "Mesh-size: " << meshSize.transpose() << "\n";
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

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        c1Argyris.plotParaview("G1Biharmonic",10000);
        mp.uniformRefine(3);
        gsWriteParaview(mp,"Geometry_init",1000,true);
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
