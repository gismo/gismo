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
# include <gsArgyris/gsErrorAnalysis/gsC1ArgyrisJumpNorm.h>

# include <gsArgyris/gsC1ArgyrisIO.h>

# include <gsMSplines/gsMappedBasis.h>

#include <boost/filesystem.hpp>

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

    bool twoPatch = false;

    gsCmdLine cmd("Solving biharmonic equation with Argyris space.");
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

    cmd.addSwitch("twoPatch","Two Patch",twoPatch);

    // Output features
    cmd.addSwitch("latex","Print the rate and error latex-ready",latex);
    cmd.addSwitch("plot", "Plotting the results!",plot);
    cmd.addSwitch("last", "Last case only!",last);
    cmd.addSwitch( "info", "Print information", info );
    cmd.addSwitch( "csv", "Save the output to a csv file", csv );
    cmd.addSwitch( "mesh", "Save the mesh to a csv file", mesh );
    cmd.addSwitch( "solution", "Save the solution to a csv file", csv_sol );
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
    gsMatrix<> jumperr(numRefine+1, mp.nInterfaces()), jumperrRate(numRefine+1, mp.nInterfaces());
    gsMatrix<> normerr(numRefine+1, 8); // 8 == 3 Norms + 3 Rates + h + dofs
    gsMatrix<> time_mat(numRefine+1, 4);
    normerr.setZero(); jumperrRate.setZero();

    gsInfo<< "(dot1=got_argyris_space, dot2=assembled, dot3=solved, dot4=got_error)\n";
    gsStopwatch time;
    for( index_t l = 0; l<=numRefine; ++l)
    {
        gsInfo<<"--------------------------------------------------------------\n";
        time.restart();
        gsSparseMatrix<> sparseMatrix_argyris;
        gsMultiBasis<> mb_argyris;

        c1Argyris.init();
        c1Argyris.createArgyrisSpace(); // Slow TODO
        /*
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
        */
        c1Argyris.getMultiBasis(mb_argyris);
        sparseMatrix_argyris = c1Argyris.getSystem();
        mappedBasis.init(mb_argyris, sparseMatrix_argyris.transpose());
        gsInfo<< "." <<std::flush;// Construction of Argyris space done
        time_mat(l, 0) = time.stop();
        gsInfo<<"\tAssembly of mapping:\t"<< time_mat(l, 0) <<"\t[s]\n";

        time.restart();
        gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(mp, mappedBasis, bcInfo, bcInfo2, source, twoPatch);
        gsInfo<<"\tDegrees of freedom:\t"<< g1BiharmonicAssembler.numDofs() <<"\n";
        g1BiharmonicAssembler.assemble();
        gsInfo<< "." <<std::flush;// Assemblying done
        time_mat(l, 1) = time.stop();
        gsInfo<<"\tSystem assembly:\t"<<time_mat(l, 1)<<"\t[s]\n";

        time.restart();
        gsSparseSolver<real_t>::CGDiagonal solver;
        solver.compute(g1BiharmonicAssembler.matrix());
        gsMatrix<real_t> solVector= solver.solve(g1BiharmonicAssembler.rhs());
        time_mat(l, 2) = time.stop();
        gsInfo<<"\tSolving system:\t\t"<<time_mat(l, 2)<<"\t[s]\n";

        time.restart();
        gsMatrix<real_t> solFull;
        g1BiharmonicAssembler.constructSolution(solVector, solFull);
        sparseMatrix_argyris = solFull.asDiagonal() * sparseMatrix_argyris;
        c1Argyris.setSystem(sparseMatrix_argyris);
        gsInfo<< "." <<std::flush;// Linear solving done

        mappedBasis.init(mb_argyris, sparseMatrix_argyris.transpose());
        gsArgyrisNorms<real_t> argyrisNorms(mp, mappedBasis, solution);
        argyrisNorms.compute();
        gsC1ArgyrisJumpNorm<real_t> c1ArgyrisJumpNorm(mp, mappedBasis, solution);
        c1ArgyrisJumpNorm.compute();
        time_mat(l, 3) = time.stop();
        gsInfo<<"\tError computations:\t"<<time_mat(l, 3)<<"\t[s]\n";

        // Collecting data
        normerr(l,0) = c1Argyris.getMinMeshSize();
        normerr(l,1) = g1BiharmonicAssembler.numDofs();

        normerr(l,2) = argyrisNorms.valueL2();
        normerr(l,4) = math::sqrt(argyrisNorms.valueH1() * argyrisNorms.valueH1() + normerr(l,2) * normerr(l,2));
        normerr(l,6) = math::sqrt(argyrisNorms.valueH2() * argyrisNorms.valueH2() +
                argyrisNorms.valueH1() * argyrisNorms.valueH1() + normerr(l,2) * normerr(l,2));
        jumperr.row(l) = c1ArgyrisJumpNorm.value();
        gsInfo<< ". " <<std::flush;// Error computations done

        // TODO Refine spaces
        c1Argyris.uniformRefine();
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

        for (size_t numInt = 0; numInt < mp.interfaces().size(); numInt++ )
        {
            gsVector<> singleInterr = jumperr.col(numInt);
            jumperrRate.block(1, numInt, numRefine, 1) = ( singleInterr.head(numRefine).array() /
                                        singleInterr.tail(numRefine).array() ).log() / std::log(2.0);
            gsInfo<<   "EoC (Interface " + util::to_string(numInt) + "): "<< std::fixed<<std::setprecision(2)
              << jumperrRate.col(numInt).transpose() <<"\n";
        }

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

    //! [Save results to csv file]
    if (csv)
    {
        std::vector<std::string> colNames, colNamesJump, colNamesTime;
        colNames.push_back("h");
        colNames.push_back("dofs");
        colNames.push_back("L2");
        colNames.push_back("Rate");
        colNames.push_back("H1");
        colNames.push_back("Rate");
        colNames.push_back("H2");
        colNames.push_back("Rate");

        for (size_t numInt = 0; numInt < mp.interfaces().size(); numInt++ )
        {
            colNamesJump.push_back("Interface" + util::to_string(numInt) );
            colNamesJump.push_back("Rate");
        }

        colNamesTime.push_back("Mapping");
        colNamesTime.push_back("System");
        colNamesTime.push_back("Solving");
        colNamesTime.push_back("Error");

        std::string name = "";
        for(index_t i = 1; i < argc; i++)
            if ((std::string) argv[i] != "--csv")
                name += (std::string) argv[i];

        std::string fullname = "";
        for(index_t i = 0; i < argc; i++)
            fullname += (std::string) argv[i] + " ";

        std::string path = "../../gismo_results/results/g" + std::to_string(geometry);

        std::string command = "mkdir " + path;
        system(command.c_str());

        path += "/" + name + ".csv";

        std::ofstream file(path);
        c1ArgyrisIO.writeLineString(file, "Command", fullname);
        c1ArgyrisIO.writeBlockMatrix(file, "Error", normerr, colNames, true);
        c1ArgyrisIO.writeBlockMatrix(file, "Jump", jumperr, jumperrRate, colNamesJump);
        c1ArgyrisIO.writeBlockMatrix(file, "Time", time_mat, colNamesTime);
        file.close();
    }
    //! [Save results to csv file]

    //! [Save mesh to csv file]
    if (mesh) {
        unsigned resolution = 100;

        std::string path = "../../gismo_results/results/g" + std::to_string(geometry);

        std::string command = "mkdir " + path;
        system(command.c_str());

        std::string name = "linedata";
        std::string name2 = "points";
        std::string name3 = "points_sorted_";

        std::ofstream file_linedata(path + "/" + name + ".csv");
        std::ofstream file_points(path + "/" + name2 + ".csv");


        size_t offset = 0;
        for (size_t np = 0; np < mp.nPatches(); np++)
        {
            std::ofstream file_points_sorted(path + "/" + name3 + util::to_string(np) + ".csv");

            gsBasis<real_t> &basis = mb.basis(np);
            gsGeometry<real_t> &Geo = mp.patch(np);

            //basis.uniformRefine();

            gsMesh<real_t> sl(basis, resolution);
            Geo.evaluateMesh(sl);

            for (typename std::vector<gsVertex<real_t> *>::const_iterator it = sl.vertices().begin();
                 it != sl.vertices().end(); ++it) {
                const gsVertex<real_t> &vertex = **it;
                file_points << vertex[0] << " ";
                file_points << vertex[1] << "\n";
                //gsInfo << vertex[0] << " "; // 3D
                //gsInfo << vertex[1] << "\n"; // 3D
            }

            for (typename std::vector<gsEdge<real_t> >::const_iterator it = sl.edges().begin();
                 it != sl.edges().end(); ++it) {
                file_linedata << std::fixed << std::setprecision(4) << it->source->getId() + offset << " " << it->target->getId() + offset << "\n";
            }

            offset += sl.numVertices();

            gsMatrix<> points;
            points.setZero(2,resolution);
            gsVector<> point_temp;
            point_temp.setLinSpaced(resolution, 0, 1);

            // v = 0
            points.row(0) = point_temp;
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            // u = 1
            points.setOnes();
            points.row(1) = point_temp;
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            // v = 1
            points.setOnes();
            points.row(0) = point_temp.reverse();
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            // u = 0
            points.setZero();
            points.row(1) = point_temp.reverse();
            file_points_sorted << Geo.eval(points).transpose() << "\n";

            file_points_sorted.close();
        }

        file_linedata.close();
        file_points.close();


    }
    //! [Save mesh to csv file]

    //! [Save solution to csv file]
    if (csv_sol) {
        //unsigned resolution = 100/mp.nPatches();
        unsigned resolution = 50;

        std::string path = "../../gismo_results/results/g" + std::to_string(geometry);

        std::string command = "mkdir " + path;
        system(command.c_str());

        //std::string name = "solution";
        //std::ofstream file_points(path + "/" + name + ".csv");

        for (size_t np = 0; np < mp.nPatches(); np++)
        {
            std::string name = "solution" + std::to_string(np);
            std::ofstream file_points(path + "/" + name + ".csv");

            gsGeometry<real_t> &Geo = mp.patch(np);

            gsMatrix<real_t> ab = Geo.support();
            gsVector<real_t> a = ab.col(0);
            gsVector<real_t> b = ab.col(1);

            gsVector<unsigned> numpoints = uniformSampleCount(a, b, resolution*resolution);
            gsMatrix<real_t> pts = gsPointGrid(a, b, numpoints);

            gsMatrix<real_t> eval_geo = Geo.eval(pts);//pts
            gsMatrix<real_t> eval_field = solVal.eval(eval_geo);

            for (index_t it = 0; it < eval_geo.cols(); ++it) {
                file_points << eval_geo(0, it) << " ";
                file_points << eval_geo(1, it) << " ";
                file_points << eval_field(0, it) << "\n";
            }

            file_points.close();
        }
        //file_points.close();
    }
    //! [Save solution to csv file]

    return EXIT_SUCCESS;

}// end main
