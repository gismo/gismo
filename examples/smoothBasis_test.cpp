/** @file gsCompositeBasis_test.h
    @brief File testing the gsCompositeBasis class.
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsSmoothPatches/gsDPatch.h>
#include <gsArgyris/gsC1Argyris.h>

using namespace gismo;

void writeLineString(std::ofstream & file, std::string command, std::string name )
{
    file<<"# " + command + "\n";
    file<<command;
    file<<"\n";
}

void writeBlockMatrix(std::ofstream & file, std::string command, gsMatrix<> matrix, std::vector<std::string> colName , bool rate = false)
{
    if(matrix.cols() != colName.size())
        gsInfo << "Something in the csv file will be wrong\n";

    file<<"# Start " + command + "\n";
    // Colname
    for (std::vector<std::string>::const_iterator it = colName.begin(); it != colName.end(); it++)
    {
        if (it != std::prev(colName.end()))
            file<<*it<<',';
        else
            file<<*it;
    }
    file<<"\n";

    // Results
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
            if (rate) {
                if (j + 1 == matrix.cols())
                    file << std::fixed << std::setprecision(2) << matrix(i, j); // last is rate
                else if (j%2 == 1 && j > 2)
                    file << std::fixed << std::setprecision(2) << matrix(i, j) << ',';
                else if (j == 0)
                    file << std::fixed << std::setprecision(5) << matrix(i, j) << ',';
                else if (j == 1)
                    file << std::fixed << std::setprecision(0) << matrix(i, j) << ',';
                else
                    file << std::scientific << std::setprecision(5) << matrix(i, j) << ',';
            }
            else {
                if (j + 1 == matrix.cols())
                    file << std::scientific << std::setprecision(5) << matrix(i, j);
                else
                    file << std::scientific << std::setprecision(5) << matrix(i, j) << ',';
            }
        }
        file<<'\n';
    }
    file<<"# End " + command + "\n";
}


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot       = false;
    bool last       = false;
    bool csv        = false;
    index_t numPoints = 100;


    index_t discrete_p = 0; // Degree of the geometry + discrete_p = polynomial degree of the discrete space
    index_t discrete_r = 1; // Regularity of the geometry + discrete_r = regularity of the discrete space

    index_t numRefine = 0; // Number of loops = refinement steps

    index_t smoothing = 0;

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

    cmd.addInt( "s", "smoothing",
                "Which smooth basis functions?", smoothing );

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
    cmd.addSwitch( "csv", "Save results to csv!", csv );
    cmd.addInt( "n", "npts", "Number of plotting points in file",  numPoints );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsOptionList optionList;
    optionList.addInt( "degreeElevate", "Number of degree elevation steps to perform before solving", discrete_p );
    optionList.addInt( "regularity", "Number of increase Continuity steps to perform before solving",  discrete_r );
    optionList.addSwitch("isogeometric", "Project the basis in isogeometric concept", isogeometric);
    optionList.addSwitch("info", "Print information", info);
    optionList.addSwitch("plot", "Plot", plot);

    optionList.addSwitch("twoPatch", "Two Patch", twoPatch);


    gsMultiPatch<> mp;

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
            case -18 ... -11:
                string_geo = "planar/twoPatches/benchmark/square_linear_bc_" + std::to_string(-geometry-10) + ".xml";
                break;

            default:
                gsInfo << "No geometry is used! \n";
                break;
        }

        gsInfo << "Filedata: " << string_geo << "\n";
        gsReadFile<>(string_geo, mp);
        mp.clearTopology();
        mp.computeTopology();
        //! [Read geometry]
    }
    else
    {
        gsReadFile<>(input, mp);
        mp.clearTopology();
        mp.computeTopology();
    }

    gsMultiPatch<> geom = mp;

    // p-refine
    if (smoothing == 0)
    {
        if (discrete_p!=0)
            mp.degreeElevate(discrete_p);

        mp.uniformRefine(3,2);
        if (last)
        {
            // h-refine
            for (int r =0; r < numRefine; ++r)
                mp.uniformRefine();

            numRefine = 0;
        }
    }

    // For smoothing 1
    gsC1Argyris<2, real_t> c1Argyris(mp, optionList);


    gsWriteParaview(mp,"mp",1000,true,false);

    gsMultiBasis<> mb(mp);
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(mb);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    space u = A.getSpace(bb2);

    // Get the map and function
    geometryMap G = A.getMap(geom);
    gsFunctionExpr<> ff;
    ff = gsFunctionExpr<>("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2);
    variable f = A.getCoeff(ff, G);

    // Solution vector and solution variable
    gsMatrix<> sVector;
    solution s = A.getSolution(u, sVector);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), linferr(numRefine+1),
        b2err(numRefine+1), b1err(numRefine+1), binferr(numRefine+1), ddofs(numRefine+1), meshSize(numRefine+1);
    gsVector<index_t> CGiter(numRefine+1);
    gsMatrix<real_t> time_mat(numRefine+1, 4);
    time_mat.setZero();

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;
    gsStopwatch time;
    for( index_t r = 0; r<=numRefine; ++r)
    {
        gsInfo<<"--------------------------------------------------------------\n";
        time.restart();
        if (smoothing == 0)
        {

            gsDPatch<2,real_t> dpatch(mp);
            dpatch.matrix_into(global2local);
            global2local = global2local.transpose();
            geom = dpatch.exportToPatches();
            mb = dpatch.localBasis();
        }
        else if (smoothing == 1)
        {
            c1Argyris.init();
            c1Argyris.createArgyrisSpace(); // Slow TODO
            c1Argyris.getMultiBasis(mb);
            global2local = c1Argyris.getSystem();
            global2local = global2local.transpose();
        }

        time_mat(r, 0) = time.stop();
        gsInfo<<"\tAssembly of mapping:\t"<<time_mat(r, 0)<<"\t[s]\n";

        bb2.init(mb,global2local);
        // gsMappedSpline<2,real_t> mspline(bb2,coefs);
        // geom = mspline.exportToPatches();


        //u.setup(bc, dirichlet::interpolation, 0);

        time.restart();
        // Initialize the system
        A.initSystem();

        ddofs.at(r) = A.numDofs();
        gsInfo<<"\tDegrees of freedom:\t"<< A.numDofs() <<"\n";

        // Compute the system matrix and right-hand side
        A.assemble(  u * u.tr(), u * f);

        // Poisson
//        A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G), - u * ilapl(f) * meas(G) );

        // Biharmonic
        //A.assemble(  ilapl(u,G) * ilapl(u,G).tr() * meas(G), - u * ilapl(f) * meas(G) );

        time_mat(r, 1) = time.stop();
        gsInfo<<"\tSystem assembly:\t"<<time_mat(r, 1)<<"\t[s]\n";

        time.restart();

        solver.compute( A.matrix() );
        sVector = solver.solve(A.rhs());
        CGiter[r] = solver.iterations();

        time_mat(r, 2) = time.stop();
        gsInfo<<"\tSolving system:\t\t"<<time_mat(r, 2)<<"\t[s]\n";
        time.restart();

        l2err[r]= math::sqrt( ev.integral( (f - s).sqNorm()*meas(G) ) / ev.integral(f.sqNorm()*meas(G)) );

        h1err[r]= l2err[r] + math::sqrt(ev.integral( ( igrad(f) - grad(s)*jac(G).inv() ).sqNorm()*meas(G) )/ev.integral( igrad(f).sqNorm()*meas(G) ) );

        linferr[r] = ev.max( f-s ) / ev.max(f);

        time_mat(r, 3) = time.stop();
        gsInfo<<"\tError computations:\t"<<time_mat(r, 3)<<"\t[s]\n"; // This takes longer for the D-patch, probably because there are a lot of points being evaluated, all containing the linear combinations of the MSplines

        meshSize[r] = mb.basis(0).getMaxCellLength();
        if (smoothing == 0)
        {
            mp.uniformRefine(1,2);
            mb = gsMultiBasis<>(mp);
        }
        else if (smoothing == 1)
            c1Argyris.uniformRefine();

    }
    //! [Solver loop]
    gsInfo << time_mat << "\n";

    gsInfo<< "\n\nCG it.: "<<CGiter.transpose()<<"\n";

    //! [Error and convergence rates]
    gsInfo<<"\n* Error\n";
    gsInfo<< "H1    "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "L2    "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "Linf  "<<std::scientific<<linferr.transpose()<<"\n";

    gsMatrix<> rate_l2, rate_h1, rate_linf;

    rate_h1 = ( h1err.head(numRefine).array() /
              h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) ;
    rate_l2 = ( l2err.head(numRefine).array() /
               l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) ;
    rate_linf = ( linferr.head(numRefine).array() /
              linferr.tail(numRefine).array() ).log().transpose() / std::log(2.0) ;

    gsInfo<<"\n* EoC\n";
    gsInfo<< "H1c   0 "<< std::fixed<<std::setprecision(2)
          << rate_h1 <<"\n";
    gsInfo<< "L2c   0 " << std::fixed<<std::setprecision(2)
          << rate_l2 <<"\n";

    gsInfo<<   "Linfc 0 "<< std::fixed<<std::setprecision(2)
          << rate_linf <<"\n";

    //! [Error and convergence rates]

    //! [Export to CSV]
    if(csv)
    {
        std::vector<std::string> colNames, colNamesTime;
        colNames.push_back("h");
        colNames.push_back("dofs");
        colNames.push_back("Linfty");
        colNames.push_back("Rate");
        colNames.push_back("L2");
        colNames.push_back("Rate");
        colNames.push_back("H1");
        colNames.push_back("Rate");

        colNamesTime.push_back("Mapping");
        colNamesTime.push_back("System");
        colNamesTime.push_back("Solving");
        colNamesTime.push_back("Error");


        gsMatrix<> results_csv(linferr.rows(), colNames.size());
        results_csv.setZero();
        results_csv.block(0,0,meshSize.rows(),1) = meshSize; //h
        results_csv.block(0,1,ddofs.rows(),1) = ddofs; //ddofs

        results_csv.col(2) = linferr.col(0);
        results_csv.block(1,3,rate_linf.cols(),1) = rate_linf.transpose();
        results_csv.col(4) = l2err.col(0);
        results_csv.block(1,5,rate_linf.cols(),1) = rate_l2.transpose();
        results_csv.col(6) = h1err.col(0);
        results_csv.block(1,7,rate_linf.cols(),1) = rate_h1.transpose();

        std::ofstream file("test2.csv");
        writeLineString(file, "Command", "-g 0 -k 4 -l 3");
        writeBlockMatrix(file, "Error", results_csv, colNames, true);
        writeBlockMatrix(file, "Time", time_mat, colNamesTime);
        file.close();

        /*
        std::ofstream file_points("points.txt");
        gsMatrix<real_t> ab = mp.patch(0).parameterRange();
        gsVector<real_t> a = ab.col(0);
        gsVector<real_t> b = ab.col(1);
        gsVector<unsigned> np = uniformSampleCount(a,b, numPoints );
        gsMatrix<real_t> pts = gsPointGrid(a,b,np) ;
        gsMatrix<real_t> sol = ff.eval(pts);

        index_t i = 0;
        for(int col = 0; col < np[1]; col++){
            for(int row = 0; row < np[0]; row++, i++) {
                file_points << std::scientific <<
                    pts(0, i) << " " << pts(1, i) << " " << sol(0, i) << "\n";
            }
            file_points << "\n";
        }
        file_points.close();
        */

    }
    //! [Export to CSV]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt   ("plot.npts"    , 500);
        ev.writeParaview( s, G, "solution");
        //ev.writeParaview( grad(s), G, "solution_grad");
        ev.writeParaview( f, G, "solution_ex");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        ev.writeParaview( (f-s), G, "error_pointwise");
        gsWriteParaview( geom, "geom",100,true);
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}