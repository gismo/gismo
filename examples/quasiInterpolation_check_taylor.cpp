/** @file quasiinterpolation_check.cpp

    @brief Different Quasi-Interpolation Schemes.

    Different quasi-interpolation-methods of a function.
    Some of the implemented methods are based on
    "Tom Lyche and Knut Morken. Spline methods draft. 2011"
    http://www.uio.no/studier/emner/matnat/ifi/INF-MAT5340/v11/undervisningsmateriale/

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): 
**/

#include <gismo.h>
#include <gsSolver/gsSolverUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{    
    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2); // Function to interpolate 
    gsFunctionExpr<real_t> mySinus1D("sin(x)",1); // Function to interpolate 

    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>:: BSplineSquare()); // create a multipatch with multibasis!???????

    // Define basis parameters
    int deg = 2;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd = deg+1; // multiplicity at the two end knots

    gsKnotVector<> kv(a, b, interior, multEnd);
    gsBSplineBasis<> bas(kv);

    gsTensorBSplineBasis<2> base2(kv,kv);
    gsTHBSplineBasis<2,real_t> thb;
    //thb = gsTHBSplineBasis<2,real_t>(base2); // create basis with knot vectors
    gsMultiBasis<> dbasis_fine;
    // dbasis_fine.addBasis(thb.clone());

    // dbasis_fine.addBasis(base2);



    gsMatrix<> coefs, coefs_2d;

    gsQuasiInterpolate<real_t>::Taylor(bas, mySinus1D, deg, coefs);
    gsQuasiInterpolate<real_t>::Taylor2D(base2, mySinus, deg, coefs_2d);

    gsInfo << coefs_2d;

    gsTensorBSpline<2,real_t> spline_taylor(base2,coefs_2d); 

    gsMatrix<> coefs_sine;
    gsL2Projection<real_t>::projectFunction(base2,mySinus,mp,coefs_sine);
    gsTensorBSpline<2,real_t> sine_func(base2,coefs_sine); 

    gsWriteParaview(spline_taylor,"approximation_Taylor"); // interpolated function mySinus
    gsWriteParaview(sine_func,"approximation_l2"); // interpolated function mySinus

}// end main
