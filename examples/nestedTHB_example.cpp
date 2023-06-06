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
    int      degree = 3;
    unsigned multEnd = degree + 1;  // multiplicity at the two end knots
    gsKnotVector<> kv(a, b, interior, multEnd);

    // ...a 2D-tensor-B-spline basis with this knot vector...
    gsTensorBSplineBasis<2,real_t> tens( kv, kv )   ;

    // ...and a 2D-THB-spline basis out of the tensor-B-spline basis.
    gsTHBSplineBasis<2,real_t> thb( tens , true);
    gsTHBSplineBasis<2,real_t> thb2( tens , false);
    gsWriteParaview(tens, "thb_level0" );
    //! [constBasis]


    gsInfo << "basis before refinement:\n" << thb << std::endl;

    gsKnotVector<> kv2(a, b, interior, multEnd);//degree-1

    // tens.uniformRefine(1,degree-1);
    // thb.addLevel(tens);

    // ...a 2D-tensor-B-spline basis with this knot vector...
    gsTensorBSplineBasis<2,real_t> tens2( kv2, kv2 );
    tens2.insertKnot(0.375,0);
    tens2.insertKnot(0.375,1);
    tens2.insertKnot(0.875,0);
    tens2.insertKnot(0.875,1);
    thb.addLevel(tens2);
    gsWriteParaview(tens2, "thb_level1" );

    thb.printBases();

    // Export the initial basis to paraview files
    gsWriteParaview(thb, "thb0_init" );


    // ONLY FOD DYADIC
    std::vector<index_t> box;
    box.push_back( 1 );
    box.push_back( 0 );
    box.push_back( 0 );
    box.push_back( 3 );
    box.push_back( 3 );

    thb.refineElements(box); //data is tailored for dyadic refinement.

    // gsMatrix<> rbox(2,2);
    // rbox.col(0)<<0.0,0.0;
    // rbox.col(1)<<0.49,0.49;
    // // rbox<< 0,0, .5, .5 ;
    // thb.refine(rbox); //data is tailored for dyadic refinement.
    // thb2.refine(rbox); //data is tailored for dyadic refinement.
        
    gsInfo << "basis after refinement:\n" << thb << std::endl;
    gsInfo << "uniform basis after refinement:\n" << thb2 << std::endl;

    thb.tree().printLeaves();
    thb.printSpaces();
    //thb.printCharMatrix();

    // Export the refined basis to paraview files
    gsWriteParaview(thb, "thb_refined_first" );
    gsWriteParaview(thb2, "thb2_refined_first" );
    gsInfo << "after refinement," << std::endl;

    gsVector<unsigned> np(2); np<<100,100;
    gsVector<> A(2); A<<0,0;
    gsVector<> B(2); B<<1,1;
    gsMatrix<> grid = gsPointGrid<>(A,B,np);
    gsMatrix<> res;
    thb2.eval_into(grid,res);
    gsVector<> sums = res.colwise().sum();

    if ((sums.array()<1-1e-12 && sums.array()>1+1e-12).count()==0)
        gsInfo<<"The basis has the partition of unity property\n";

    return EXIT_SUCCESS;
}
