/** @file thbSplineBasis_example.cpp

    @brief Tutorial on gsTHBSplineBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

// Look also in basis_example and bSplineBasis_example.

//! [Include namespace]
#include <string>
#include <gismo.h>
//! [Include namespace]

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // --------------- construction of a THB-spline basis ---------------

    //! [constBasis]
    // Set up and construct a knot vector...
    real_t a = 0;                   // starting knot
    real_t b = 1;                   // ending knot
    unsigned interior = 3;          // number of interior knots
    int      degree = 2;
    unsigned multEnd = degree + 1;  // multiplicity at the two end knots
    gsKnotVector<> kv(a, b, interior, multEnd);

    // ...a 2D-tensor-B-spline basis with this knot vector...
    gsTensorBSplineBasis<2,real_t> tens( kv, kv );

    // ...and a 2D-THB-spline basis out of the tensor-B-spline basis.
    gsTHBSplineBasis<2,real_t> thb( tens );
    //! [constBasis]


    gsInfo << "basis before refinement:\n" << thb << std::endl;

    // Export the initial basis to paraview files
    gsWriteParaview(thb, "thb0_init" );

    //! [refViaStdVec]
    std::vector<index_t> box;
    box.push_back( 1 );
    box.push_back( 2 );
    box.push_back( 0 );
    box.push_back( 8 );
    box.push_back( 4 );

    thb.refineElements(box);
    //! [refViaStdVec]

    // Export the refined basis to paraview files
    gsWriteParaview(thb, "thb1_refined" );
    gsInfo << "after refinement," << std::endl;

    //! [stdOpsCout]
    gsInfo << "this basis is:\n" << thb << std::endl;
    //! [stdOpsCout]

    // --------------- "standard" evaluations ---------------
    //! [stdOpsStd]
    gsMatrix<real_t> u(2,3);
    u(0,0) = 0.95;
    u(1,0) = 0.05;

    u(0,1) = 0.95;
    u(1,1) = 0.3;

    u(0,2) = 0.6;
    u(1,2) = 0.9;

    gsMatrix<index_t> resActives;
    gsMatrix<real_t>   resEvals;

    thb.active_into( u, resActives);
    thb.eval_into(   u, resEvals);

    gsInfo << "active functions: \n" << resActives << std::endl;
    gsInfo << "their values:     \n" << resEvals   << std::endl;
    //! [stdOpsStd]

    gsInfo << std::endl;

    // --------------- index-computations ---------------

    //! [indexTransfForw]
    std::vector<unsigned> tmpFlatIndices;
    std::vector<int>      tmpLevels;

    gsInfo << "transform indices\n";
    gsInfo << "global/hier.index  ->  flat tensor index (on level)" << std::endl;
    for( unsigned i = 27; i <= 35; i++)
    {
        // print computed indices/levels
        gsInfo << i;
        gsInfo << "  ->  ";
        gsInfo << thb.flatTensorIndexOf(i);
        gsInfo << "  ( " << thb.levelOf(i) << " )" << std::endl;

        // store indices/levels for reverse transformation later
        tmpFlatIndices.push_back( thb.flatTensorIndexOf(i) );
        tmpLevels.push_back(  thb.levelOf(i) );
    }
    //! [indexTransfForw]


    //! [indexTransfBack]
    gsInfo << std::endl;
    gsInfo << "reverse index transformation\n";
    gsInfo << "flat tensor index (on level)  ->  global/hier.index" << std::endl;
    for( unsigned i = 0; i < tmpLevels.size(); i++ )
    {
        // print global/hierarchical indices
        gsInfo << tmpFlatIndices[i] << "  ( " << tmpLevels[i] << " )";
        gsInfo << "  ->  ";
        gsInfo << thb.flatTensorIndexToHierachicalIndex( tmpFlatIndices[i], tmpLevels[i] ) << std::endl;
    }
    //! [indexTransfBack]

    gsInfo << std::endl;

    // --------------- some gsHTensorBasis-specific functions ---------------
    //! [stdOpsHTens]
    gsVector<index_t> resLevels;
    gsMatrix<index_t> resLowerCorner;

    thb.getLevelUniqueSpanAtPoints(u, resLevels, resLowerCorner);

    gsInfo << "levels:        " << std::endl << resLevels.transpose() << std::endl;
    gsInfo << "lower corners: " << std::endl << resLowerCorner        << std::endl;
    //! [stdOpsHTens]

    gsInfo << std::endl;

    // print the underlying tree
    //! [stdOpsHTensTree]
    gsMatrix<index_t> resUpperCorner;

    thb.tree().getBoxes( resLowerCorner, resUpperCorner, resLevels);

    gsInfo << "levels:        " << std::endl << resLevels      << std::endl;
    gsInfo << "lower corners: " << std::endl << resLowerCorner << std::endl;
    gsInfo << "upper corners: " << std::endl << resUpperCorner << std::endl;
    //! [stdOpsHTensTree]




    // --------------- 2nd local refinement ---------------
    gsInfo << std::endl << std::endl;

    //! [refViaStdVec2]
    box.clear();
    box.push_back( 2 );
    box.push_back( 2 );
    box.push_back( 4 );
    box.push_back( 6 );
    box.push_back( 10 );

    thb.refineElements(box);

    gsInfo << "after 2nd refinement, this basis is:\n" << thb << std::endl;
    //! [refViaStdVec2]

    gsWriteParaview(thb, "thb2_refined" );

    // --------------- plot basis after 1 refinement ---------------
    //! [Plot in Paraview]
    if( plot )
    {
        // Run paraview
        gsFileManager::open("thb1_refined.pvd");
    }
    //! [Plot in Paraview]
    else
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    }
    return EXIT_SUCCESS;
}
