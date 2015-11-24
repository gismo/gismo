/** @file tutorialTHBSplineBasis.cpp

    @brief Tutorial on gsTHBSplineBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

// Look also in tuturialBasis and tutorialBSplineBasis.

//! [Include namespace]
#include <string>
#include <gismo.h>
//! [Include namespace]

using namespace gismo;
using namespace std;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = true;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    const bool ok = cmd.getValues(argc,argv);
    if (!ok) { gsWarn << "Error during parsing the command line!\n"; return 0;}
    //! [Parse command line]

    // --------------- construction of a THB-spline basis ---------------

    //! [constBasis]
    // Set up and construct a knot vector...
    real_t a = 0;                   // starting knot
    real_t b = 1;                   // ending knot
    int degree = 2;
    unsigned interior = 3;          // number of interior knots
    unsigned multEnd = degree + 1;  // multiplicity at the two end knots
    gsKnotVector<> kv(a, b, interior, multEnd);

    // ...a 2D-tensor-B-spline basis with this knot vector...
    gsTensorBSplineBasis<2,real_t,gsKnotVector<> > tens( kv, kv );

    // ...and a 2D-THB-spline basis out of the tensor-B-spline basis.
    gsTHBSplineBasis<2,real_t> thb( tens );
    //! [constBasis]


    cout << "basis before refinement:\n" << thb << endl;

    // Export the initial basis to paraview files
    gsWriteParaview(thb, "thb0_init" );

    //! [refViaStdVec]
    std::vector<unsigned> box;
    box.push_back( 1 );
    box.push_back( 2 );
    box.push_back( 0 );
    box.push_back( 8 );
    box.push_back( 4 );

    thb.refineElements(box);
    //! [refViaStdVec]

    // Export the refined basis to paraview files
    gsWriteParaview(thb, "thb1_refined" );
    cout << "after refinement," << endl;

    //! [stdOpsCout]
    std::cout << "this basis is:\n" << thb << std::endl;
    //! [stdOpsCout]

    // --------------- "standard" evaluations ---------------
    //! [stdOpsStd]
    gsMatrix<real_t> u(2,2);
    u(0,0) = 0.99;
    u(1,0) = 0.01;

    u(0,1) = 0.6;
    u(1,1) = 0.9;

    gsMatrix<unsigned> resActives;
    gsMatrix<real_t>   resEvals;

    thb.active_into( u,resActives);
    thb.eval_into(   u,resEvals);

    std::cout << "active functions: " << std::endl << resActives << std::endl;
    std::cout << "their values:     " << std::endl << resEvals   << std::endl;
    //! [stdOpsStd]

    std::cout << std::endl;

    // --------------- some gsHTensorBasis-specific functions ---------------
    //! [stdOpsHTens]
    gsVector<unsigned> resLevels;
    gsMatrix<unsigned> resLowerCorner;

    thb.getLevelUniqueSpanAtPoints(u, resLevels, resLowerCorner);
    std::cout << "levels:        " << std::endl << resLevels.transpose() << std::endl;
    std::cout << "lower corners: " << std::endl << resLowerCorner        << std::endl;
    //! [stdOpsHTens]

    std::cout << std::endl;

    // print the underlying tree
    //! [stdOpsHTensTree]
    gsMatrix<unsigned> resUpperCorner;
    thb.tree().getBoxes( resLowerCorner, resUpperCorner, resLevels);
    std::cout << "levels:        " << std::endl << resLevels      << std::endl;
    std::cout << "lower corners: " << std::endl << resLowerCorner << std::endl;
    std::cout << "upper corners: " << std::endl << resUpperCorner << std::endl;
    //! [stdOpsHTensTree]


    // --------------- 2nd local refinement ---------------
    std::cout << std::endl << std::endl;

    //! [refViaStdVec2]
    box.clear();
    box.push_back( 2 );
    box.push_back( 2 );
    box.push_back( 4 );
    box.push_back( 6 );
    box.push_back( 10 );

    thb.refineElements(box);

    std::cout << "after 2nd refinement, this basis is:\n" << thb << std::endl;
    //! [refViaStdVec2]

    gsWriteParaview(thb, "thb2_refined" );

    //! [Plot in Paraview]
    if( true ) // (plot)
    {
        // Run paraview
        return system("paraview thb1_refined.pvd &");
    }
    //! [Plot in Paraview]
    else
    {
        gsInfo<<"Quitting.. No output created, re-run with --plot to get a ParaView "
                "file containing Plotting image data.\n";
        return 0;
    }
}
