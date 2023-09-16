/** @file fitting_example.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, D. Mokris
*/

#include <gismo.h>


#ifdef gsParasolid_ENABLED
#include <gsParasolid/gsWriteParasolid.h>
#endif

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save     = false;
    index_t numURef   = 3;
    index_t iter      = 2;
    index_t deg_x     = 2;
    index_t deg_y     = 2;
    index_t u_numKnots = 0;
    index_t v_numKnots = 0;
    index_t maxPcIter = 1;
    real_t lambda = 1e-07;
    real_t tolerance = 1e-02;
    real_t threshold = tolerance;
    index_t extension = 2;
    real_t refPercent = 0.1;
    std::string fn = "fitting/deepdrawingC.xml";

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addInt("c", "parcor", "Steps of parameter correction", maxPcIter);
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addInt("m", "knots_x", "number of interior knots in u-direction", u_numKnots);
    cmd.addInt("n", "knots_y", "number of interior knots in v-direction", v_numKnots);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    cmd.addReal("t", "threshold", "error threshold (special valule -1)", threshold);
    cmd.addReal("p", "refPercent", "percentage of points to refine in each iteration", refPercent);
    cmd.addInt("q", "extension", "extension size", extension);
    cmd.addInt("r", "urefine", "initial uniform refinement steps", numURef);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);
    cmd.addString("d", "data", "Input sample data", fn);

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

    //! [Read data]
    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data]

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(xyz, "points");

    ///////////////////////////// START /////////////////////////////


    std::ofstream file_results;
    file_results.open("r"+internal::to_string(numURef)+"_fiting_example_results.csv");
    file_results << "step, dofs, min, max, mse, percentage\n";

    index_t step = maxPcIter;
    // for (index_t step=0; step <= maxPcIter; step ++)
    {
    std::ofstream max_results;
    std::ofstream mse_results;
    max_results.open("r"+internal::to_string(numURef)+"pc"+internal::to_string(step)+"_max_results.csv");
    max_results << "dofs, max\n";

    mse_results.open("r"+internal::to_string(numURef)+"pc"+internal::to_string(step)+"_mse_results.csv");
    mse_results << "dofs, mse\n";

    // This is for outputing an XML file, if requested
    gsFileData<> fd;

    // Check if matrix sizes are OK
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, u_numKnots, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, v_numKnots, deg_y+1 ) ;

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

    const std::vector<real_t> & errors = ref.pointWiseErrors();
    std::vector<real_t> errors2;
    real_t sum_of_errors2;

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

    for(int i = 0; i <= iter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo<<"Iteration "<<i<<".."<<"\n";

        time.restart();
        //ref.nextIteration(tolerance, threshold, maxPcIter);
        ref.nextIteration(tolerance, threshold, step);
        time.stop();

        gsMesh<> mesh(ref.result()->basis());
        gsMatrix<> uv_fitting = ref.returnParamValues() ;
        gsWriteParaview(mesh, internal::to_string(i+1) + "_iter_mesh_pc" + internal::to_string(step));
        gsWriteParaview(*ref.result(), internal::to_string(i+1) + "_iter_geo_pc" + internal::to_string(step), 100000, false, true);
        gsWriteParaviewPoints(uv_fitting, internal::to_string(i+1) + "_iter_fitting_parameters_pc" + internal::to_string(step));


        // compute mean squared error
        ref.get_Error(errors2, 0);
        sum_of_errors2 = std::accumulate(errors2.begin(), errors2.end(), 0.0);

        gsInfo<<"Fitting time: "<< time <<"\n";


        index_t dofs = ref.result()->basis().size();
        real_t minPointError = ref.minPointError();
        real_t maxPointError = ref.maxPointError();
        real_t mseError = sum_of_errors2/errors2.size();
        real_t percentagePoint = 100.0 * ref.numPointsBelow(tolerance)/errors.size();

        gsMatrix<real_t> matErrors(4,errors.size());
        for(size_t el = 0; el < errors.size(); el++)
        {
            matErrors(0,el) = uv_fitting(0,el);
            matErrors(1,el) = uv_fitting(1,el);
            matErrors(2,el) = 0;
            matErrors(3,el) = errors[el];
        }

        gsWriteParaviewPoints(matErrors, internal::to_string(i+1) + "colors_iter_fitting_parameters_pc" + internal::to_string(step));

        // This is checks the two ways of computing MSE for consistency.
        // TODO: remove!
        const std::vector<real_t> ptErrs = ref.pointWiseErrors();
        real_t mse_Dominik(0);
        for(auto it = ptErrs.begin(); it != ptErrs.end(); ++it)
            mse_Dominik += *it * *it;
        mse_Dominik /= real_t(ptErrs.size());

        gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
        gsInfo<<"Parameter correction steps: "<< step << "\n";
        gsInfo    << "DOFs         : "<< dofs <<"\n";
        // gsInfo<<"Min distance : "<< ref.minPointError() <<" / ";
        // gsInfo<<"Max distance : "<< ref.maxPointError() <<"\n";
        std::cout << "Min distance : "<< minPointError        << std::scientific <<"\n";
        std::cout << "Max distance : "<< maxPointError        << std::scientific <<"\n";
        std::cout << "MSE (Dominik): "<< mse_Dominik          << std::scientific <<"\n";
        std::cout << "MSE (Sofia)  : "<< mseError             << std::scientific <<"\n";
        std::cout << "        RMSE : "<< math::sqrt(mseError) << std::scientific <<"\n";
        gsInfo<<"Points below tolerance: "<< percentagePoint  <<"%.\n";


        //file_results << std::to_string(step) << "," << std::to_string(dofs) << "," << std::ostringstream(std::to_string(minPointError)).str() << "," << std::ostringstream(std::to_string(maxPointError)).str() << "," << std::ostringstream(std::to_string(mseError)).str() << "," << std::ostringstream(std::to_string(percentagePoint)).str() << "\n";
        file_results << step << ","
					 << dofs << ","
					 << minPointError << std::scientific << ","
					 << maxPointError << std::scientific << ","
					 << mseError << std::scientific << ","
					 << percentagePoint << "\n";
        max_results <<  std::to_string(dofs) << "," << std::ostringstream(std::to_string(maxPointError)).str() << "\n";
        mse_results <<  std::to_string(dofs) << "," << std::ostringstream(std::to_string(mseError)).str() << "\n";



        if ( ref.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<< i <<" iterations.\n";
            break;
        }
    } // iterative loop

    max_results.close();
    mse_results.close();

    gsInfo<<"----------------\n";

    if ( save )
    {
        gsInfo<<"Done. Writing solution to file fitting_out_pc" << step << ".xml\n";
        fd << *ref.result() ;

        fd.dump("fitting_out_pc"+internal::to_string(step));

#ifdef gsParasolid_ENABLED
        gsTHBSpline<2>* result = static_cast<gsTHBSpline<2>*>(ref.result());
		extensions::gsWriteParasolid<real_t>(*result, "result_assembly");
		gsTensorBSpline<2> TP_spline;
        result->convertToBSpline(TP_spline);
        extensions::gsWritePK_SHEET(TP_spline, "result_tp");
#endif
    }
    else
        gsInfo << "Done. No output created, re-run with --save to get a xml "
                  "file containing the solution.\n";
    } // parameter correction steps
    file_results.close();
    return 0;
}
