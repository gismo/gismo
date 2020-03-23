/** @file g1fitting_example.cpp

    @brief Demonstrates adaptive fitting of g1 basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 4;

    // Options with default values
    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    //! [Read data]
    // Surface fitting
    // Expected input is a file with matrices with:
    // id 0:  u,v   -- parametric coordinates, size 2 x N
    // id 1:  x,y,z -- corresponding mapped values, size 3 x N
    index_t N = 1000;
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);

    gsMatrix<> uv, xyz;


    //! [Read data]

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
        ref.nextIteration(tolerance, threshold);
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
