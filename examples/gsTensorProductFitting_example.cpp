/** @file fitting_data_driven.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>
#include<gsHSplines/gsHDDFitting.h>
// #include <gsHSplines/gsHBox.h>
// #include <gsHSplines/gsHBoxContainer.h>
// #include <gsHSplines/gsHBoxUtils.h>
//#include <gsModeling/gsParametrization.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save         = false;
    bool manualRefT   = false;
    bool constantRefT = false;
    bool CellAvgRef   = false;
    index_t maxURef   = 0;
    index_t iter      = 0;
    index_t deg_x     = 2;
    index_t deg_y     = 2;
    index_t maxPcIter = 0;
    real_t lambda     = 1e-09;
    real_t threshold  = 1e-02;
    real_t tolerance  = 1e-02;
    index_t extension = 2;
    real_t refPercent = 0.;
    std::string fn = "../../testdata/ShipHull/shipHullPts_scale_m55.xml";

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("manual", "Set refinement threshold manually", manualRefT);
    cmd.addSwitch("constant", "Set constant refinement threshold", constantRefT);
    cmd.addSwitch("avg", "Set refinement scheme.\nIf true: average cell error;\notherwise: max cell error;", CellAvgRef);
    cmd.addInt("c", "parcor", "Steps of parameter correction", maxPcIter);
    cmd.addInt("i", "iter", "number of iterations", iter);
    cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
    cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    cmd.addReal("t", "threshold", "error threshold (special valule -1)", threshold);
    cmd.addReal("p", "refPercent", "percentage of points to refine in each iteration", refPercent);
    cmd.addInt("q", "extension", "extension size", extension);
    cmd.addInt("r", "urefine", "initial uniform refinement steps", maxURef);
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

    // Specify extension size in u and v cells
    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);




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
    gsWriteParaviewPoints(uv, "input_parameters");
    gsWriteParaviewPoints(xyz, "input_points");


    // Print settings summary
    gsInfo<<"----------------------------------------------------------------\n";
    gsInfo<<"Fitting "<< xyz.cols() <<" samples.\n";
    gsInfo<<"----------------------------------------------------------------\n";
    gsInfo<<"Cell extension     : "<< ext[0]<<" "<<ext[1]<<".\n";
    if ( threshold >= 0.0 )
        gsInfo<<"Ref. threshold     : "<< threshold<<".\n";
    else
        gsInfo<<"Cell refinement    : "<< 100*refPercent<<"%%.\n";
    gsInfo<<"Error tolerance    : "<< tolerance<<".\n";
    gsInfo<<"Smoothing parameter: "<< lambda<<".\n";

    // Check if matrix sizes are OK
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
                   u_max = uv.row(0).maxCoeff(),
                   v_min = uv.row(1).minCoeff(),
                   v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, 0, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, 0, deg_y+1 ) ;

    for (index_t numURef = 0; numURef < maxURef + 1; ++numURef)
    {
        std::string outname = "deg" + internal::to_string(deg_x) + "_ref" + internal::to_string(numURef);

        if (maxURef < 1){
          lambda = 0.;
        }
        // Create a tensor-basis and apply compute initial polynomial approximation
        gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
        T_tbasis.uniformRefine( (1<<numURef)-1 );
        // Create Initial hierarchical basis
        gsTHBSplineBasis<2>  THB ( T_tbasis ) ;
        gsMesh<> TPmesh(THB);
        gsWriteParaview(TPmesh, outname + "TPmesh");


        // Create polynomial refinement object
        gsHFitting<2, real_t> polymodel( uv, xyz, THB, 0., ext, lambda);
        // gsHFitting<2, real_t> ref( uv, xyz, THB, refPercent, ext, lambda);

        // point-wise polynomial error.

        const std::vector<real_t> & errors = polymodel.pointWiseErrors();
        std::vector<real_t> errors2;
        real_t sum_of_2errors;

        gsStopwatch time;
        gsInfo<<"----------------------------------------------------------------\n";
        gsInfo<<"Polynomial fitting:\n";
        time.restart();
        polymodel.computeApproximation(maxPcIter);
        time.stop();
        gsInfo<<"Polynomial fitting time: "<< time <<"\n";


        polymodel.get_Error(errors2, 0);
        sum_of_2errors = std::accumulate(errors2.begin(), errors2.end(), 0.0);

        real_t mse = 0.;
        const std::vector<real_t> & Perrors = polymodel.pointWiseErrors();
        for(index_t el = 0; el != Perrors.size(); el++){
            mse += std::pow(Perrors[el], 2);
        }
        mse = mse / Perrors.size();


        gsInfo<<"Polynomial basis:\n"<< polymodel.result()->basis() <<"\n";
        gsInfo<<"DOFs            : "<< polymodel.result()->basis().size() <<"\n";
        gsInfo<<"Min distance    : "<< polymodel.minPointError() <<"\n";
        gsInfo<<"Max distance    : "<< polymodel.maxPointError() <<"\n";
        gsInfo<<"MSE    error    : "<< mse/3 <<"\n";
        gsInfo<<"Points below tolerance: "<< 100.0 * polymodel.numPointsBelow(tolerance)/errors.size()<<"%.\n";
        gsInfo<<"----------------------------------------------------------------\n";


        // gsMesh<> plotmesh;
        // polymodel.result()->toMesh(plotmesh, uv.cols());
        // gsMatrix<> params(2, plotmesh.numVertices());
        // for(size_t i=0; i<plotmesh.numVertices(); i++)
        // {
        //     size_t index = plotmesh.unsorted(i);
        //     params.col(i) = getParameterPoint(index).transpose();
        // }
        // gsWriteParaview(plotmesh, "paramsField_ref" + internal::to_string(numURef), params);
        index_t numSamples(100);
        bool plot_mesh = true;
        std::string pname("param");
        // gsInfo << "*polymodel.result():\n" << *polymodel.result() << "\n";
        gsGeometry<>::uPtr geo = polymodel.result()->clone();
        gsParamField<real_t> pfld(*geo);
        gsWriteParaview(*geo, pfld, pname + "_geo", numSamples);
        gsWriteParaview(*geo, pname + "_msh", numSamples, plot_mesh, false);


        gsVector<unsigned> numPtsVec(2);
        numPtsVec<<50,50;
        gsMatrix<> ab;
        ab = polymodel.result()->support();
        gsVector<> a = ab.col(0);
        gsVector<> b = ab.col(1);
        gsMatrix<> gridParams = gsPointGrid(a,b, numPtsVec);
        gsMatrix<> gridEval;
        polymodel.result()->eval_into(gridParams, gridEval);
        gsFileData<> ptsOut;
        ptsOut << gridEval;
        ptsOut << gridParams;
        ptsOut.dump(outname + "dhd");

        if ( save )
        {
            gsWriteParaview(*polymodel.result(), outname, numSamples, true);
            // This is for outputing an XML file, if requested
            gsFileData<> fd;
            gsInfo<<"Done. Writing solution to file .xml\n";
            fd << *polymodel.result();

            fd.dump(outname);
        }
        else
            gsInfo << "No output created, re-run with --save to get a xml "
                      "file containing the solution.\n";
    } // iteration on number of uniform refinement

    return 0;
}
