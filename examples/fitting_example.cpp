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
    bool save     = false;  // save
    index_t numURef   = 0;  // r
    index_t numKnots   = 0; // n
    index_t nx   = -1; // a
    index_t ny   = -1; // b
    index_t iter      = 2;  // i
    index_t deg_x     = 2;  // x
    index_t deg_y     = 2;  // y
    index_t maxPcIter = 0;  // c
    real_t lambda = 1e-06;  // s
    real_t threshold = 1e-02; // t
    real_t tolerance = 1e-02; // e
    index_t extension = 2;  // q
    real_t refPercent = 0.1;// p
    std::string fn = "fitting/deepdrawingC.xml"; // d

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
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    cmd.addReal("t", "threshold", "error threshold (special valule -1)", threshold);
    cmd.addReal("p", "refPercent", "percentage of points to refine in each iteration", refPercent);
    cmd.addInt("q", "extension", "extension size", extension);
    cmd.addInt("r", "urefine", "initial uniform refinement steps", numURef);
    cmd.addInt("n", "iknots", "number of interior knots in each direction", numKnots);
    cmd.addInt("a", "uknots", "number of interior knots in u-direction", nx);
    cmd.addInt("b", "vknots", "number of interior knots in v-direction", ny);
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

    gsWriteParaviewPoints(uv, "uv");
    gsWriteParaviewPoints(xyz, "xyz");

    // This is for outputing an XML file, if requested
    gsFileData<> fd;

    // Check if matrix sizes are OK
    // GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  // "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    if( nx < 0)
      nx = numKnots;
    if( ny < 0)
      ny = numKnots;
    gsKnotVector<> u_knots (u_min, u_max, nx, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, ny, deg_y+1 ) ;

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
        ref.nextIteration(tolerance, threshold, maxPcIter);
        time.stop();
        gsInfo<<"Fitting time: "<< time <<"\n";

        gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
        gsInfo<<"Min distance : "<< ref.minPointError() <<" / ";
        gsInfo<<"Max distance : "<< ref.maxPointError() <<"\n";
        gsInfo<<"Points below tolerance: "<< 100.0 * ref.numPointsBelow(tolerance)/errors.size()<<"%.\n";

        if ( ref.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<<i<<" iterations.\n";
            break;
        }
    }

    gsInfo<<"----------------\n";

    if ( save )
    {
        gsInfo<<"Done. Writing solution to file fitting_out.xml\n";
        fd << *ref.result() ;

        fd.dump("fitting_out");
    }
    else
        gsInfo << "Done. No output created, re-run with --save to get a xml "
                  "file containing the solution.\n";

    return 0;
}
