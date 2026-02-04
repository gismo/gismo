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
    unsigned interior = 6;          // number of interior knots
    int      degree = 2;
    unsigned multEnd = degree + 1;  // multiplicity at the two end knots
    gsKnotVector<> kv(a, b, interior, multEnd);
    
    // Create 1D B-spline and THB-spline basis
    gsBSplineBasis<real_t> bspline(kv);
    gsTHBSplineBasis<1,real_t> thb(bspline);

    gsInfo << "THB spline basis size: " << thb.size() << "\n";

    // gsInfo << "basis before refinement:\n" << thb << "\n";

    // Export the initial basis to paraview files
    // gsWriteParaview(thb, "basisthb/thb0_init" );

    // // Check support of basis functions
    // for (index_t i = 0; i < thb.size(); ++i)
    // {
    //     gsInfo << "Basis function " << i << " supports: "
    //         << thb.support(i) << "\n";

    //     // Get coefficients of this basis function
    //     gsVector<> coefs = thb.getCoefs(i);
    //     gsInfo << "  Coefficients: " << coefs.transpose() << "\n";
    // }
    gsTHBSplineBasis<1,real_t> *thb_crs = dynamic_cast<gsTHBSplineBasis<1,real_t> *>(&thb);


    std::vector<index_t> box, box2;
    box.push_back(1); // level to refine
    box.push_back(4); // start parametric coordinate
    box.push_back(10); // end parametric coordinate

    thb.refineElements(box);

    
    // Export the refined basis to paraview files
    // gsWriteParaview(thb, "basisthb/thb_refined_first" );
    
    // gsInfo << "after refinement," << std::endl;
    gsInfo << "THB spline basis size: " << thb.size() << "\n";

    gsFunctionExpr<real_t> myFunction("exp(-100*((x-0.5)^2))", 1);
    
    gsMultiBasis<> dbasis_crs;
    gsMatrix<> C_coarse_QI;
    dbasis_crs.addBasis(thb_crs->clone());
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_crs.basis(0), myFunction,C_coarse_QI); 
    gsTHBSpline<1,real_t> thb_coar_exp(*thb_crs,C_coarse_QI); 
    gsVector<> pt1(1,1);
    gsVector<> pt2(1,1);
    gsVector<> pt3(1,1);

    pt1<< 0.365187;
    pt2<< 0.392857;
    pt3<< 0.420527;

    gsInfo << "thb_coar_exp at pt " << pt1.transpose() << ": " << thb_coar_exp.eval(pt1) << "\n";
    gsInfo << "thb_coar_exp at pt " << pt2.transpose() << ": " << thb_coar_exp.eval(pt2) << "\n";
    gsInfo << "thb_coar_exp at pt " << pt3.transpose() << ": " << thb_coar_exp.eval(pt3) << "\n";


    gsMatrix<> C_fine_QI;

    gsMultiBasis<> dbasis_fine;

    gsTHBSplineBasis<1,real_t> *thb_graded = dynamic_cast<gsTHBSplineBasis<1,real_t> *>(&thb);
    dbasis_fine.addBasis(thb_graded->clone());
    gsDebugVar(dbasis_fine.basis(0));
    // gsWriteParaview(*thb_graded,"basispaper/basisfine",1000); // this is the fine mesh
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), myFunction,C_fine_QI); 

    gsInfo<< "size " << dbasis_fine.basis(0).size() << "\n" ;

    // gsInfo<<"fine QI: \n"<<C_fine_QI<<"\n";

    gsGeometry<>::uPtr fine_QI = dbasis_fine.basis(0).makeGeometry((C_fine_QI));
    gsTHBSpline<1,real_t> thb_exp(*thb_graded,C_fine_QI); 

    gsWriteParaview(thb_exp,"function_fine", 200); 

    

    // // // Check support of basis functions
    // for (index_t i = 0; i < thb.size(); ++i)
    // {
    //     gsInfo << "Basis function " << i << " supports: "
    //         << thb.support(i) << "\n";

    //     // Get coefficients of this basis function
    //     gsVector<> coefs = thb.getCoefs(i);
    //     gsInfo << "  Coefficients: " << coefs.transpose() << "\n";
    // }

    // Loop over the refined THB basis
gsVector<> C = C_fine_QI;  // coefficients used to build thb_exp


// Anchors (Greville points) in parametric domain
gsMatrix<> greville = thb_exp.basis().anchors();

// Loop through each basis function
for (index_t i = 0; i < thb_exp.basis().size(); ++i)
{
    // Parametric support of the i-th basis function
    gsMatrix<> supp = thb_exp.basis().support(i);

    // Greville point in parametric domain
    gsVector<> pos = greville.col(i);

    // Corresponding coefficient value
    real_t c = C[i];

    gsInfo << "Basis " << i
           << "\n Level: " << thb_exp.basis().levelOf(i)
           << " Support interval: [" << supp(0,0) << ", " << supp(0,1) << "]"
           << " Position: " << pos.transpose()
           << " Coefficient: " << c
           << "\n";
}



    // for (index_t i = 0; i < thb.size(); ++i)
    // {
    //     if (thb.isTruncated(i))
    //         gsInfo<<"Coeff " << i <<": "<< thb.getCoefs(i).transpose() << "\n";
    // }

    // gsInfo << "size" << thb.size()  << "\n";

    

    // box2.push_back(1); // level to refine
    // box2.push_back(0); // start parametric coordinate
    // box2.push_back(14); // end parametric coordinate



    // thb.refineElements(box2);
    // // Export the refined basis to paraview files
    // gsWriteParaview(thb, "basisthb/thb_level1_" );
    // gsInfo << "after refinement2," << std::endl;


//     // Define element interval
// real_t a = 5.0 / 14.0;
// real_t b = 3.0 / 7.0;

// // Center of the element
// real_t center = 0.5 * (a + b);
// gsInfo << "Center of the element [" << a << ", " << b << "]: " << center << "\n";

// // Find active basis functions on this element (pseudo code, depends on your library)
// gsMatrix<> active = thb.active(center); // 'thb' is your basis object

// gsInfo << "Active basis functions on element [" << a << ", " << b << "]:\n" << active.transpose() << "\n";


// Define interval endpoints
// real_t aa = 5.0 / 14.0;
// real_t bb = 3.0 / 7.0;

// // Calculate center (scalar)
// real_t center_val = 0.5 * (aa + bb);

// // Wrap it into a gsMatrix (1x1 matrix)
// gsMatrix<real_t> center(1, 1);
// center << center_val;  // assign scalar value into 1x1 matrix

// gsInfo << "Center of the element [" << aa << ", " << bb << "]: " << center(0,0) << "\n";

// // Now call your method with gsMatrix
// gsMatrix<index_t> active = thb.active(center);
// gsInfo << "Active basis functions indices on element [" << aa << ", " << bb << "]:\n"
//        << active.transpose() << "\n";


    // //! [stdOpsCout]
    // gsInfo << "this basis is:\n" << thb << std::endl;
    // //! [stdOpsCout]

    // // --------------- "standard" evaluations ---------------
    // //! [stdOpsStd]
    // gsMatrix<real_t> u(2,3);
    // u(0,0) = 0.95;
    // u(1,0) = 0.05;

    // u(0,1) = 0.95;
    // u(1,1) = 0.3;

    // u(0,2) = 0.6;
    // u(1,2) = 0.9;

    // gsMatrix<index_t> resActives;
    // gsMatrix<real_t>   resEvals;

    // thb.active_into( u, resActives);
    // thb.eval_into(   u, resEvals);

    // gsInfo << "active functions: \n" << resActives << std::endl;
    // gsInfo << "their values:     \n" << resEvals   << std::endl;
    // //! [stdOpsStd]

    // gsInfo << std::endl;

    // // --------------- index-computations ---------------

    // //! [indexTransfForw]
    // std::vector<unsigned> tmpFlatIndices;
    // std::vector<int>      tmpLevels;

    // gsInfo << "transform indices\n";
    // gsInfo << "global/hier.index  ->  flat tensor index (on level)" << std::endl;
    // for( unsigned i = 27; i <= 35; i++)
    // {
    //     // print computed indices/levels
    //     gsInfo << i;
    //     gsInfo << "  ->  ";
    //     gsInfo << thb.flatTensorIndexOf(i);
    //     gsInfo << "  ( " << thb.levelOf(i) << " )" << std::endl;

    //     // store indices/levels for reverse transformation later
    //     tmpFlatIndices.push_back( thb.flatTensorIndexOf(i) );
    //     tmpLevels.push_back(  thb.levelOf(i) );
    // }
    // //! [indexTransfForw]


    // //! [indexTransfBack]
    // gsInfo << std::endl;
    // gsInfo << "reverse index transformation\n";
    // gsInfo << "flat tensor index (on level)  ->  global/hier.index" << std::endl;
    // for( unsigned i = 0; i < tmpLevels.size(); i++ )
    // {
    //     // print global/hierarchical indices
    //     gsInfo << tmpFlatIndices[i] << "  ( " << tmpLevels[i] << " )";
    //     gsInfo << "  ->  ";
    //     gsInfo << thb.flatTensorIndexToHierachicalIndex( tmpFlatIndices[i], tmpLevels[i] ) << std::endl;
    // }
    // //! [indexTransfBack]

    // gsInfo << std::endl;

    // // --------------- some gsHTensorBasis-specific functions ---------------
    // //! [stdOpsHTens]
    // gsVector<index_t> resLevels;
    // gsMatrix<index_t> resLowerCorner;

    // thb.getLevelUniqueSpanAtPoints(u, resLevels, resLowerCorner);

    // gsInfo << "levels:        " << std::endl << resLevels.transpose() << std::endl;
    // gsInfo << "lower corners: " << std::endl << resLowerCorner        << std::endl;
    // //! [stdOpsHTens]

    // gsInfo << std::endl;

    // // print the underlying tree
    // //! [stdOpsHTensTree]
    // gsMatrix<index_t> resUpperCorner;

    // thb.tree().getBoxes( resLowerCorner, resUpperCorner, resLevels);

    // gsInfo << "levels:        " << std::endl << resLevels      << std::endl;
    // gsInfo << "lower corners: " << std::endl << resLowerCorner << std::endl;
    // gsInfo << "upper corners: " << std::endl << resUpperCorner << std::endl;
    // //! [stdOpsHTensTree]




    // // --------------- 2nd local refinement ---------------
    // gsInfo << std::endl << std::endl;

    // //! [refViaStdVec2]
    // box.clear();
    // box.push_back( 2 );
    // box.push_back( 2 );
    // box.push_back( 4 );
    // box.push_back( 6 );
    // box.push_back( 10 );

    // thb.refineElements(box);

    // gsInfo << "after 2nd refinement, this basis is:\n" << thb << std::endl;
    // //! [refViaStdVec2]

    // gsWriteParaview(thb, "thb_refined_second" );

    // boxSide side(1);
    // gsMatrix<index_t> result = thb.boundaryOffset(1,0);
    // gsInfo<<"Basis indices along side "<<side<<": "<<result.transpose()<<"\n";

    // for (index_t k=1; k<5; k++)
    // {
    //     boxCorner corner(k);
    //     gsInfo<<"Basis function index in corner "<<corner<<": "<<thb.functionAtCorner(corner)<<"\n";
    // }

    // // --------------- plot basis after 1 refinement ---------------
    // //! [Plot in Paraview]
    // if( plot )
    // {
    //     // Run paraview
    //     gsFileManager::open("thb1_refined.pvd");
    // }
    // //! [Plot in Paraview]
    // else
    // {
    //     gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
    //               "file containing the solution.\n";
    // }
    return EXIT_SUCCESS;
}
