#/** @file biharmonic_multiPatch_example.cpp

    @brief Example using the Argyris space or the Nitsche's method
           for solving biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <omp.h>


//! [Include namespace]
#include <gismo.h>

#include <gsAssembler/gsBiharmonicMethods.h>

#include <gsC1Basis/gsC1ArgyrisIO.h>

#include <sstream>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}



using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool last       = false;
    index_t discrete_p = 3; // Polynomial degree of the discrete space
    index_t discrete_r = 1; // Regularity of the discrete space

    // If no gluing data is chosen, then gluingData_p = discrete_p - 1
    // and gluingData_r = gluingData_p - 1
    index_t gluingData_p = -1; // Polynomial degree of the discrete space
    index_t gluingData_r = -1; // Regularity of the discrete space

    index_t numRefine = 0; // Number of loops = refinement steps

    index_t geometry = 0;
    std::string input;

    bool neumann = false;

    bool info = false;
    bool latex = false;
    bool csv = false;
    bool mesh = false;
    bool csv_sol = false;

    std::string method = "argyris";

    real_t mu = 10; // TODO replace to general

    gsCmdLine cmd("Solving biharmonic equation with Argyris space.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "g", "geometry", "Which geometry",  geometry );

    // For the discrete space
    cmd.addInt( "p", "discreteDegree",
                "Polynomial degree of the discrete space", discrete_p );
    cmd.addInt( "r", "discreteRegularity",
                "Regularity of the discrete space",  discrete_r );

    // For the gluing data space
    cmd.addInt( "P", "gluingDataDegree",
                "Polynomial degree of the gluing data space", gluingData_p );
    cmd.addInt( "R", "gluingDataRegularity",
                "Regularity of the gluing data space",  gluingData_r );

    // For computing the convergence rates
    cmd.addInt( "l", "loop",
                "Number of Uniform h-refinement steps to perform before solving",  numRefine );

    // Which method one want to use
    cmd.addSwitch("neumann","Compute the biharmonic with neumann bdy",neumann);


    // How to compute the biharmonic equation
    cmd.addString("m", "method", "Which method do you want to use: nitsche or argyris", method);

    // Output features
    cmd.addSwitch("latex","Print the rate and error latex-ready",latex);
    cmd.addSwitch("plot", "Plotting the results!",plot);
    cmd.addSwitch("last", "Last case only!",last);
    cmd.addSwitch( "info", "Print information", info );
    cmd.addSwitch( "csv", "Save the output to a csv file", csv );
    cmd.addSwitch( "mesh", "Save the mesh to a csv file", mesh );
    cmd.addSwitch( "solution", "Save the solution to a csv file", csv_sol );

    cmd.addReal("y" , "mu", "Mu for Nitsche", mu); // TODO more general
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
    mb = gsMultiBasis<>(mp);
    // p-refine
    for (size_t np = 0; np < mp.nPatches(); ++np)
        mb.basis(np).setDegree(optionList.getInt("discreteDegree"));

    if (last)
    {
        // h-refine
        for (int l =0; l < numRefine; ++l)
            mb.uniformRefine();

        numRefine = 0;
    }
    else
    {
        index_t p = mb.minCwiseDegree();
        index_t r = optionList.getInt("discreteRegularity");
        mb.uniformRefine(1, p-r);
    }
    //! [Initialise the discrete space]

    //! [Solver loop]
    gsMatrix<> jumperr(numRefine+1, mp.nInterfaces()), jumperrRate(numRefine+1, mp.nInterfaces()),
               nitschePen(numRefine+1, mp.nInterfaces());
    gsMatrix<> normerr(numRefine+1, 10); // 10 == 4 Norms + 4 Rates + h + dofs
    gsMatrix<> time_mat(numRefine+1, 4);
    normerr.setZero(); jumperrRate.setZero(); nitschePen.setZero();

    gsInfo<< "(dot1=got_argyris_space, dot2=assembled, dot3=solved, dot4=got_error)\n";
    gsStopwatch time;
    for( index_t l = 0; l<=numRefine; ++l) {

        gsBiharmonic<real_t> *biharmonic;
        if (method == "argyris")
            biharmonic = new gsBiharmonicArgyris<real_t>(mp, mb, optionList);
        else if (method == "nitsche")
            biharmonic = new gsBiharmonicNitsche<real_t>(mp, mb, optionList);
        else
            biharmonic = new gsBiharmonicArgyris<real_t>(mp, mb, optionList); // Only for avoiding the warning

        gsInfo << "--------------------------------------------------------------\n";

        time.restart();
        biharmonic->init();
        gsInfo << "." << std::flush;// Construction of Argyris space done
        time_mat(l, 0) = time.stop();
        gsInfo << "\tAssembly of mapping:\t" << time_mat(l, 0) << "\t[s]\n";

        time.restart();
        biharmonic->assemble(bcInfo, bcInfo2, source);
        gsInfo<<"\tDegrees of freedom:\t"<< biharmonic->numDofs() <<"\n";
        gsInfo << "." << std::flush;// Assemblying done
        time_mat(l, 1) = time.stop();
        gsInfo << "\tSystem assembly:\t" << time_mat(l, 1) << "\t[s]\n";

        time.restart();
        gsSparseSolver<real_t>::CGDiagonal solver;
        solver.compute(biharmonic->matrix());
        gsMatrix<real_t> solVector = solver.solve(biharmonic->rhs());
        time_mat(l, 2) = time.stop();
        gsInfo << "\tSolving system:\t\t" << time_mat(l, 2) << "\t[s]\n";
        gsInfo << "\tNumber of non-zeros:\t" << biharmonic->matrix().nonZeros() << "\n";
        gsInfo << "." << std::flush;// Linear solving done

        biharmonic->constructSolution(solVector); // Construct Solution

        time.restart();
        biharmonic->error(solution);
        time_mat(l, 3) = time.stop();
        gsInfo << "\tError computations:\t" << time_mat(l, 3) << "\t[s]\n";
        gsInfo << ". " << std::flush;// Error computations done

        // Collecting data
        normerr(l,0) = mb.basis(0).getMinCellLength(); // TODO general
        normerr(l,1) = biharmonic->numDofs();

        normerr(l,2) = biharmonic->valueL2();
        normerr(l,4) = math::sqrt(biharmonic->valueH1() * biharmonic->valueH1() + normerr(l,2) * normerr(l,2));
        normerr(l,6) = math::sqrt(biharmonic->valueH2() * biharmonic->valueH2() +
                                          biharmonic->valueH1() * biharmonic->valueH1() + normerr(l,2) * normerr(l,2));
        normerr(l,8) = biharmonic->valueJumpSum();
        jumperr.row(l) = biharmonic->valueJump();
        if (method == "nitsche")
            nitschePen.row(l) = biharmonic->valuePenalty();
        gsInfo<< ". " <<std::flush;// Error computations done

        // Refine the discrete space
        index_t p = mb.minCwiseDegree();
        index_t r = optionList.getInt("discreteRegularity");
        mb.uniformRefine(1, p-r);
    }
    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo << "\nDofs: " << normerr.col(1).transpose() << "\n";
    gsInfo << "Mesh-size: " << normerr.col(0).transpose() << "\n";
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<normerr.col(2).transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<normerr.col(4).transpose()<<"\n";
    gsInfo<< "H2 error: "<<std::scientific<<normerr.col(6).transpose()<<"\n";
    gsInfo<< "Jump error: "<<std::scientific<<normerr.col(8).transpose()<<"\n";
    gsInfo<< "Jump error single: \n"<<std::scientific<<jumperr.transpose()<<"\n";

    if (method == "nitsche")
        gsInfo<< "Penalty value: \n"<<std::scientific<<nitschePen.transpose()<<"\n";

    if (numRefine>0)
    {

        normerr.block(1,3, numRefine, 1) = ( normerr.col(2).head(numRefine).array() /
                                             normerr.col(2).tail(numRefine).array() ).log() / std::log(2.0); // L2
        normerr.block(1, 5, numRefine, 1)  = ( normerr.col(4).head(numRefine).array() /
                                               normerr.col(4).tail(numRefine).array() ).log() / std::log(2.0); // H1
        normerr.block(1, 7, numRefine, 1)  = ( normerr.col(6).head(numRefine).array() /
                                               normerr.col(6).tail(numRefine).array() ).log() / std::log(2.0); // H2
        normerr.block(1, 9, numRefine, 1)  = ( normerr.col(8).head(numRefine).array() /
                                               normerr.col(8).tail(numRefine).array() ).log() / std::log(2.0); // H2

        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              << normerr.col(3).transpose() <<"\n";
        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              << normerr.col(5).transpose() <<"\n";
        gsInfo<<   "EoC (H2): "<< std::fixed<<std::setprecision(2)
              << normerr.col(7).transpose() <<"\n";
        gsInfo<<   "EoC (Jump): "<< std::fixed<<std::setprecision(2)
              << normerr.col(9).transpose() <<"\n";

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
        //c1Argyris.plotParaview("G1Biharmonic",10000);
        mp.uniformRefine(3);
        gsWriteParaview(mp,"Geometry_init",1000,true);
    }
    //! [Export visualization in ParaView]

    //! [Save results to csv file]
    if (csv)
    {
        std::vector<std::string> colNames, colNamesJump, colNamesTime, colNamesInter;
        colNames.push_back("h");
        colNames.push_back("dofs");
        colNames.push_back("L2");
        colNames.push_back("Rate");
        colNames.push_back("H1");
        colNames.push_back("Rate");
        colNames.push_back("H2");
        colNames.push_back("Rate");
        colNames.push_back("Jump");
        colNames.push_back("Rate");

        for (size_t numInt = 0; numInt < mp.interfaces().size(); numInt++ )
        {
            colNamesJump.push_back("Interface" + util::to_string(numInt) );
            colNamesInter.push_back("Interface" + util::to_string(numInt) );
            colNamesJump.push_back("Rate");
        }

        colNamesTime.push_back("Mapping");
        colNamesTime.push_back("System");
        colNamesTime.push_back("Solving");
        colNamesTime.push_back("Error");

        std::string fullname = "";
        for(index_t i = 0; i < argc; i++)
            fullname += (std::string) argv[i] + " ";

        std::string smallNumber = "";
        if (mu < 1)
            while(mu < 1)
            {
                mu *= 10;
                smallNumber += "0";
            }

        std::string name = "-g" + std::to_string(geometry)
                           + "--" + (method == "argyris" ? "argyris" : "nitsche") +
                           "-p" + std::to_string(discrete_p) +
                           "-r" + std::to_string(discrete_r) +
                           "-l" + std::to_string(numRefine)
                           + (method == "nitsche" ? "-y" + smallNumber + to_string_with_precision(mu,0) : "");

        std::string path = "../../gismo_results/results/g" + std::to_string(geometry) + "/nitsche";

        std::string command = "mkdir " + path;
        system(command.c_str());

        path += "/" + name + ".csv";

        std::ofstream file(path);
        c1ArgyrisIO.writeLineString(file, "Command", fullname);
        c1ArgyrisIO.writeBlockMatrix(file, "Error", normerr, colNames, true);
        c1ArgyrisIO.writeBlockMatrix(file, "Jump", jumperr, jumperrRate, colNamesJump);
        c1ArgyrisIO.writeBlockMatrix(file, "Time", time_mat, colNamesTime);

        if (method == "nitsche")
            c1ArgyrisIO.writeBlockMatrix(file, "Penalty", nitschePen, colNamesInter);

        file.close();
    }
    //! [Save results to csv file]

    //! [Save mesh to csv file]
    if (mesh) {
        std::string path = "../../gismo_results/results/g" + std::to_string(geometry);

        std::string command = "mkdir " + path;
        system(command.c_str());

        c1ArgyrisIO.saveMesh(path, mp, mb, 30);
    }
    //! [Save mesh to csv file]

    //! [Save solution to csv file]
    if (csv_sol) {
        //unsigned resolution = 100/mp.nPatches();
        std::string path = "../../gismo_results/results/g" + std::to_string(geometry);

        std::string command = "mkdir " + path;
        system(command.c_str());

        c1ArgyrisIO.saveSolution(path, mp, mb, solVal, 50);
    }
    //! [Save solution to csv file]

    if (latex)
    {
        for (index_t i = 0; i < normerr.rows(); i++)
        {
            printf("%-5f & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n",
                   normerr(i,0), normerr(i,2), normerr(i,3), normerr(i,4), normerr(i,5),
                   normerr(i,6), normerr(i,7), normerr(i,8), normerr(i,9));
        }
    }


    return EXIT_SUCCESS;

}// end main
