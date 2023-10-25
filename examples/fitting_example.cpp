/** @file fitting_example.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save     = false;
    index_t numURef   = 3;
    index_t iter      = 2;
    index_t deg_x     = 3;
    index_t deg_y     = 3;
    index_t maxPcIter = 0;
    bool correctBoundary = false;
    index_t sepIndex  = -1;
    real_t lambda = 0;
    real_t threshold = 1e-02;
    real_t tolerance = 1e-02;
    index_t extension = 2;
    real_t refPercent = 0.1;
    std::string fn = "";
    std::string outname = "fitting_out";

    std::vector<index_t> modevec;
    std::vector<index_t> corners;

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("b", "b_correction", "apply parater correction also on the boundary points", correctBoundary);
    cmd.addInt("c", "parcor", "Steps of parameter correction", maxPcIter);
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    cmd.addReal("t", "threshold", "error threshold (special valule -1)", threshold);
    cmd.addReal("p", "refPercent", "percentage of points to refine in each iteration", refPercent);
    cmd.addInt("q", "extension", "extension size", extension);
    cmd.addInt("r", "urefine", "initial uniform refinement steps", numURef);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addInt("n", "interiors", "number of interior points belonging to the input point cloud", sepIndex);
    cmd.addMultiInt("m", "modes", "Modes to select", modevec);
    cmd.addMultiInt("g", "corners", "point cloud corners", corners);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (deg_x < 1)
    { gsInfo << "Degree x must be positive.\n";  return 0;}
    if (deg_y < 1)
    { gsInfo << "Degree y must be positive.\n"; return 0;}
    if (extension < 0)
    { gsInfo << "Extension must be non negative.\n"; return 0;}

    if ( tolerance < 0 )
    {
        gsInfo << "Error tolerance cannot be negative, setting it to default value.\n";
        tolerance = 1e-02;
    }

    if (threshold > 0 && threshold > tolerance )
    {
        gsInfo << "Refinement threshold is over tolerance, setting it the same as tolerance.\n";
        threshold = tolerance;
    }


    gsMatrix<> uv, xyz;
    if (fn.size() > 0)
    {
    // -x 2 -y 2 -i 7 -s 2.5e-5 --save -c 1 -n 8137 -e 1e-3 -r 0
    // ! [Read data]
    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    gsFileData<> fd_in(fn);
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    }
    else
    {
      gsFunctionExpr<> source("1/(1.5*exp(sqrt((10*(x-0.5)-3)^2+(10*(y-0.5)-3)^2))) + 1/(1.5*exp(sqrt((10*(x-0.5)+3)^2+(10*(y-0.5)+3)^2))) + 1/(1.5*exp(sqrt((10*(x-0.5))^2+(10*(y-0.5))^2)))",2);

      gsVector<> lower(2), upper(2);
      lower(0) = 0; lower(1) = 0;
      upper(0) = 1; upper(1) = 1;
      uv = uniformPointGrid(lower, upper, 10000);
      gsInfo << uv.rows() << " x " << uv.cols() << "\n";

      // source.eval_into(uv, points);
      gsMatrix<> fval = source.eval(uv);
      xyz.resize(3, fval.cols());
      xyz << uv.row(0), uv.row(1), fval.row(0);
    }

    gsInfo << uv.rows() << " x " << uv.cols() << "\n";
    gsInfo << xyz.rows() << " x " << xyz.cols() << "\n";
    // gsInfo << points << "\n";

    // gsWriteParaviewPoints(xyz, "peaks");
    //! [Read data]

    // gsWriteParaviewPoints(uv, "parameters");
    // gsWriteParaviewPoints(xyz, "points");

    // This is for outputing an XML file, if requested
    gsFileData<> fd, pout;

    // Check if matrix sizes are OK
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    if (sepIndex > xyz.cols() || sepIndex < 0)
    {
      gsInfo << "Apply " << maxPcIter << " parameter correction step to the whole point cloud.\n";
      sepIndex = xyz.cols();
    }
    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    gsInfo << "Parameter domain:\n";
    gsInfo << "u_min : " << u_min << "\n";
    gsInfo << "u_max : " << u_max << "\n";
    gsInfo << "v_min : " << v_min << "\n";
    gsInfo << "v_max : " << v_max << "\n";
    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, 0, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, 0, deg_y+1 ) ;

    // Create a tensor-basis nad apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create Initial hierarchical basis

    gsTHBSplineBasis<2>  THB ( T_tbasis ) ;
    //gsHBSplineBasis<2>  THB ( T_tbasis ) ;

    // Specify extension size in u and v cells
    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);

    // Create hierarchical refinement object
    gsHFitting<2, real_t> ref( uv, xyz, THB, refPercent, ext, lambda);

    std::vector<real_t> errors;
    std::vector<real_t> errors2;
    real_t sum_of_2errors;

    // Print settings summary
    gsInfo<<"Fitting "<< xyz.cols() <<" samples.\n";
    gsInfo<<"----------------\n";
    gsInfo<<"Cell extension     : "<< ext[0]<<" "<<ext[1]<<".\n";
    if ( threshold >= 0.0 )
        gsInfo<<"Ref. threshold     : "<< threshold<<".\n";
    else
        gsInfo<<"Cell refinement    : "<< 100*refPercent<<"%%.\n";
    gsInfo<<"Error tolerance    : "<< tolerance<<".\n";
    gsInfo<<"Smoothing parameter: "<< lambda<<".\n";

    gsStopwatch time;

    std::ofstream file_out;
    file_out.open("LS_results.csv");
    file_out << "m, deg, pen, dofs, it, min, max, mse, rmse\n";

    for(int i = 0; i <= iter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo<<"Iteration "<<i<<".."<<"\n";

        time.restart();
        //ref.nextIteration(tolerance, threshold, maxPcIter);
        if (correctBoundary){
          gsInfo << "Correct also the boundary points\n";
          ref.nextIterationSepBoundary(tolerance, threshold, maxPcIter, sepIndex);
        }
        else{
          if(sepIndex < 0){
            gsInfo << "Apply "<< maxPcIter << "steps of parameter correction to the whole pointcloud.\n";
            ref.nextIteration(tolerance, threshold, maxPcIter);
          }
          else{
            gsInfo << "No parameter correction on boundary points.\n";
          }
          ref.nextIterationFixedBoundary(tolerance, threshold, maxPcIter, sepIndex);
        }
        // if(i == 0){
        //   //ref.nextIteration(tolerance, threshold, maxPcIter);
        //   ref.nextIterationSepBoundary(tolerance, threshold, maxPcIter, sepIndex);
        //   //ref.nextIterationFixedBoundary(tolerance, threshold, maxPcIter, sepIndex);
        // }
        // else{
        //   ref.nextIterationSepBoundary(tolerance, threshold, 1, sepIndex);
        //   //ref.nextIteration(tolerance, threshold, maxPcIter);
        //   //ref.nextIterationFixedBoundary(tolerance, threshold, 0, sepIndex);
        // }
        time.stop();

        gsMesh<> mesh(ref.result()->basis());
        gsWriteParaview(mesh, "LS" + internal::to_string(i) + "_iter_mesh");
        gsWriteParaview(*ref.result(), "LS" + internal::to_string(i) + "geo", 100000, false, true);
        //gsWriteParaview(*ref.result(), internal::to_string(i+1) + "_iter_geo", 100000, false, true);

        errors = ref.pointWiseErrors();
        gsMatrix<> colorPoints(xyz.rows() + 1, xyz.cols());
        gsMatrix<> errorsMat = ref.result()->pointWiseErrors(uv, xyz);
        // for(index_t err = 0; err < errors.size(); err++){
        //     errorsMat(0,err) = errors[err];
        // }
        colorPoints.row(0) = xyz.row(0);
        colorPoints.row(1) = xyz.row(1);
        colorPoints.row(2) = xyz.row(2);
        colorPoints.row(3) = errorsMat;


        gsWriteParaviewPoints(colorPoints, "LS" + internal::to_string(i) + "color");


        ref.get_Error(errors2, 0);
        sum_of_2errors = std::accumulate(errors2.begin(), errors2.end(), 0.0);

        std::vector<real_t> sol_min_max_mse = ref.result()->MinMaxMseErrors(uv, xyz);

        gsInfo<<"Fitting time: "<< time <<"\n";
        gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
        gsInfo<<"DOFs         : "<< ref.result()->basis().size() <<"\n";
        gsInfo<<"Min distance : "<< sol_min_max_mse[0] << std::scientific << "\n";//ref.minPointError() <<"\n";
        std::cout << "Max distance : "<< sol_min_max_mse[1] << std::scientific << "\n";//<< ref.maxPointError() << std::scientific <<"\n";
        std::cout << "MSE    error : "<< sol_min_max_mse[2] << std::scientific << "\n";//<< sum_of_2errors/errors2.size() << std::scientific <<"\n";
        std::cout << "rMSE    error : "<< math::sqrt(sol_min_max_mse[2]) << std::scientific << "\n";//<< sum_of_2errors/errors2.size() << std::scientific <<"\n";
        gsInfo<<"Points below tolerance: "<< 100.0 * ref.numPointsBelow(tolerance)/errors.size()<<"%.\n";

        //"m, deg, pen, dofs, it, min, max, mse\n";
        file_out << xyz.cols() << "," << deg_x << "," << lambda << std::scientific << ","
                 << ref.result()->basis().size() << "," << i << ","
                 << sol_min_max_mse[0] << std::scientific << ","
                 << sol_min_max_mse[1] << std::scientific << ","
                 << sol_min_max_mse[2] << std::scientific << ","
                 << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
                 << sum_of_2errors/errors2.size() << std::scientific << "\n";

        if ( ref.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<<i<<" iterations.\n";
            break;
        }
    }

    gsInfo<<"----------------\n";

    if ( save )
    {
        gsMatrix<> fitting_out_parameters = ref.returnParamValues();

        gsInfo << "Plot parametric field on the geometry.\n";
        index_t numSamples(100000);
        bool plot_mesh = true;
        std::string pname("param");
        gsGeometry<>::uPtr geo = ref.result()->clone();
        gsParamField<real_t> pfld(*geo);
        // gsWriteParaview(*geo, pfld, pname + "_geo", numSamples);
        // gsWriteParaview(*geo, pname + "_msh", numSamples, plot_mesh, false);



        // gsInfo << "Plot the max approximation error of the geometry.\n";
        //
        // gsMatrix<> fieldCoefs = ref.getFieldMaxError(fitting_out_parameters, errors);
        // gsInfo << "Error field coefs:\n" << fieldCoefs << "\n";

        gsInfo<<"Done. Writing solution to file fitting_out.xml\n";

        gsInfo << "Solution:\n" << *ref.result() << "\n";
        fd << *ref.result() ;
        pout << fitting_out_parameters;
        pout << xyz;

        // gsWriteParaviewPoints(fitting_out_parameters, "fitting_out_parameters");
        for(index_t idx = 0; idx < modevec.size(); idx++)
        {
          std::string idxparamname = std::to_string(modevec[idx]) + "_parameter";
          gsMatrix<> print_parameter(2,1);
          print_parameter << fitting_out_parameters.col(modevec[idx]);
          // gsWriteParaviewPoints(print_parameter, idxparamname);
          std::string idxpointname = std::to_string(modevec[idx]) + "_point";
          gsMatrix<> print_point(3,1); print_point << xyz.col(modevec[idx]);
          // gsWriteParaviewPoints(print_point, idxpointname);
        }
        fd.dump("fitting_out");
        pout.dump("data_out");
    }
    else
        gsInfo << "Done. No output created, re-run with --save to get a xml "
                  "file containing the solution.\n";

    return 0;
}
