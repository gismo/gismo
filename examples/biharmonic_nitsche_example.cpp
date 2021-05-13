#/** @file biharmonic_nitsche_example.cpp

    @brief Example using the Nitsche method for solving biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <omp.h>

//! [Include namespace]
#include <gismo.h>

# include <gsAssembler/gsBiharmonicNitscheAssembler.h>

// Maybe shift or do sth else
# include <gsArgyris/gsErrorAnalysis/gsC1NitscheNorms.h>


using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool last       = false;
    index_t discrete_p = 3; // Polynomial degree of the discrete space
    index_t discrete_r = 1; // Regularity of the discrete space

    index_t numRefine = 0; // Number of loops = refinement steps

    index_t geometry = 0;
    std::string input;

    bool isogeometric = false;
    bool exactGD = false;
    bool neumann = false;

    bool info = false;
    bool latex = false;
    bool csv = false;
    bool mesh = false;
    bool csv_sol = false;

    bool interpolation = false;

    bool C1Vertex = false;
    bool twoPatch = false;
    bool simplified = false;

    real_t mu = 100;

    gsCmdLine cmd("Solving biharmonic equation with Nitsche method.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );

    // For the discrete space
    cmd.addInt( "p", "discreteDegree",
                "Polynomial degree of the discrete space", discrete_p );
    cmd.addInt( "r", "discreteRegularity",
                "Regularity of the discrete space",  discrete_r );

    // For computing the convergence rates
    cmd.addInt( "l", "loop",
                "Number of Uniform h-refinement steps to perform before solving",  numRefine );

    // Which method one want to use
    cmd.addSwitch( "exactGD", "To compute the gluing data exact", exactGD );
    cmd.addSwitch( "isogeometric", "Project the basis in isogeometric concept", isogeometric );
    cmd.addSwitch("neumann","Compute the biharmonic with neumann bdy",neumann);

    cmd.addSwitch( "interpolation", "Interpolate the basis functions", interpolation );

    cmd.addSwitch("C1Vertex","Using C^1 vertex space",C1Vertex);
    cmd.addSwitch("twoPatch","Two Patch case",twoPatch);
    cmd.addSwitch("simplified","Simplified Argyris space",simplified);

    // Output features
    cmd.addSwitch("latex","Print the rate and error latex-ready",latex);
    cmd.addSwitch("plot", "Plotting the results!",plot);
    cmd.addSwitch("last", "Last case only!",last);
    cmd.addSwitch( "info", "Print information", info );
    cmd.addSwitch( "csv", "Save the output to a csv file", csv );
    cmd.addSwitch( "mesh", "Save the mesh to a csv file", mesh );
    cmd.addSwitch( "solution", "Save the solution to a csv file", csv_sol );

    cmd.addReal("m" , "mu", "Mu for Nitsche", mu);
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

    gsOptionList optionList = cmd.getOptionList();
    gsInfo << optionList << "\n";
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

            case 1000 ... 1999:
                string_geo = "planar/geometries/g" + std::to_string(geometry) + ".xml";
                break;

            case -8 ... -1:
                string_geo = "planar/twoPatches/benchmark/square_linear_" + std::to_string(-geometry) + ".xml";
                break;
            case -18 ... -11:
                string_geo = "planar/twoPatches/benchmark/square_linear_bc_" + std::to_string(-geometry-10) + ".xml";
                break;

            default:
                gsInfo << "No geometry is used! \n";
                break;
        }
    }
    else
    {
        string_geo = input;
    }

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();
    //! [Read geometry]

    //! [Check the input data]

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
    mb = gsMultiBasis<>(mp);
    index_t numDegree = discrete_p - mb.maxCwiseDegree();
    for (int i = 0; i < numDegree; ++i)
        mb.degreeElevate();

    if (last)
    {
        // h-refine
        for (int l =0; l < numRefine; ++l)
            mb.uniformRefine();

        numRefine = 0;
    }
    else
        mb.uniformRefine();
    //! [Initialise the discrete space]

    gsInfo << "Basis: " << mb.basis(0) << "\n";

    gsSparseSolver<>::CGDiagonal solver;

    //! [Solver loop]
    gsMatrix<> jumperr(numRefine+1, mp.nInterfaces()), jumperrRate(numRefine+1, mp.nInterfaces());
    gsMatrix<> normerr(numRefine+1, 8); // 8 == 3 Norms + 3 Rates + h + dofs
    gsMatrix<> time_mat(numRefine+1, 4);
    normerr.setZero(); jumperrRate.setZero();

    gsMultiPatch<> mpsol;
    gsInfo<< "(dot1=none, dot2=assembled, dot3=solved, dot4=got_error)\n";
    gsStopwatch time;
    for( index_t l = 0; l<=numRefine; ++l)
    {
        gsInfo<<"--------------------------------------------------------------\n";

        time.restart();
        gsBiharmonicNitscheAssembler<real_t> biharmonicNitscheAssembler(mp, mb, bcInfo, bcInfo2, source, optionList);
        gsInfo<<"\tDegrees of freedom:\t"<< biharmonicNitscheAssembler.numDofs() <<"\n";
        biharmonicNitscheAssembler.assemble();
        gsInfo<< "." <<std::flush;// Assemblying done

        time_mat(l, 1) = time.stop();
        gsInfo<<"\tSystem assembly:\t"<<time_mat(l, 1)<<"\t[s]\n";

        time.restart();
        gsSparseSolver<real_t>::LU solver;
        solver.compute(biharmonicNitscheAssembler.matrix());
        gsMatrix<real_t> solVector= solver.solve(biharmonicNitscheAssembler.rhs());
        time_mat(l, 2) = time.stop();
        gsInfo<<"\tSolving system:\t\t"<<time_mat(l, 2)<<"\t[s]\n";

        time.restart();
        biharmonicNitscheAssembler.constructSolution(solVector, mpsol);
        gsInfo<< "." <<std::flush;// Linear solving done

        // TODO Error
        gsC1NitscheNorms<real_t> c1NitscheNorms(mp, mpsol, solution);
        c1NitscheNorms.compute();
        real_t errorH2Semi = c1NitscheNorms.valueH2();
        real_t errorH1Semi = c1NitscheNorms.valueH1();
        normerr(l,2) = c1NitscheNorms.valueL2();
        normerr(l,4) = math::sqrt(errorH1Semi*errorH1Semi + normerr(l,2)*normerr(l,2));
        normerr(l,6) = math::sqrt(errorH2Semi*errorH2Semi + errorH1Semi*errorH1Semi + normerr(l,2)*normerr(l,2));

        time_mat(l, 3) = time.stop();
        gsInfo<<"\tError computations:\t"<<time_mat(l, 3)<<"\t[s]\n";

        // Collecting data
        /*
        normerr(l,0) = mb.basis(0).getMinCellLength();
        normerr(l,1) = g1BiharmonicAssembler.numDofs();

        normerr(l,2) = argyrisNorms.valueL2();
        normerr(l,4) = math::sqrt(argyrisNorms.valueH1() * argyrisNorms.valueH1() + normerr(l,2) * normerr(l,2));
        normerr(l,6) = math::sqrt(argyrisNorms.valueH2() * argyrisNorms.valueH2() +
                                  argyrisNorms.valueH1() * argyrisNorms.valueH1() + normerr(l,2) * normerr(l,2));
        jumperr.row(l) = c1ArgyrisJumpNorm.value();
         */
        gsInfo<< ". " <<std::flush;// Error computations done

        // TODO Refine spaces
        mb.uniformRefine();
    }
    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo << "\nDofs: " << normerr.col(1).transpose() << "\n";
    gsInfo << "Mesh-size: " << normerr.col(0).transpose() << "\n";
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<normerr.col(2).transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<normerr.col(4).transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<normerr.col(6).transpose()<<"\n";
    gsInfo<< "Jump error: "<<std::scientific<<jumperr.transpose()<<"\n";

    if (numRefine>0)
    {

        normerr.block(1,3, numRefine, 1) = ( normerr.col(2).head(numRefine).array() /
                                             normerr.col(2).tail(numRefine).array() ).log() / std::log(2.0); // L2
        normerr.block(1, 5, numRefine, 1)  = ( normerr.col(4).head(numRefine).array() /
                                               normerr.col(4).tail(numRefine).array() ).log() / std::log(2.0); // H1
        normerr.block(1, 7, numRefine, 1)  = ( normerr.col(6).head(numRefine).array() /
                                               normerr.col(6).tail(numRefine).array() ).log() / std::log(2.0); // H2

        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << normerr.col(3).transpose() <<"\n";
        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              << normerr.col(5).transpose() <<"\n";
        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              << normerr.col(7).transpose() <<"\n";

    }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsField<> solField(mp, mpsol);
        gsWriteParaview<>(solField, "BiharmonicNitsche", 5000);
        mp.uniformRefine(3);
        gsWriteParaview(mp,"Geometry_init",1000,true);
    }
    //! [Export visualization in ParaView]

    //! [Save results to csv file]
    //! [Save results to csv file]

    //! [Save mesh to csv file]
    //! [Save mesh to csv file]

    //! [Save solution to csv file]
    //! [Save solution to csv file]

    return EXIT_SUCCESS;

}// end main
