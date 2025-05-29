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
    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2);
    gsFunctionExpr<> myPolyQuad("-5*x^2 + 3*x*y + y^2 + 22*x + y - 4",2);

    int p = 2;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd = p+1; // multiplicity at the two end knots

    gsKnotVector<> kv(a, b, interior, multEnd);

    gsTensorBSplineBasis<2> bas(kv,kv);
    // gsAdaptiveMeshing<2,real_t> mesher; //????

    int numRef = 5;
    bool passed = true;

    gsTHBSplineBasis<2,real_t> thb, thb1, thb2, thb3, thb_fine;

    gsMultiPatch<> mp;
 
    thb = gsTHBSplineBasis<2,real_t>(bas); // coarse basis

    gsMatrix<> coefs, C_coarse, C_fine, C_coarse2; // for projection

    // Get COARSE coefficients    
    std::vector<real_t> error_list, h_list;
    gsQuasiInterpolate<real_t>::localIntpl(thb, mySinus, coefs); //???? no se si tiene sentido?
    gsMatrix<> C_coarse_initial = coefs;
    gsTHBSpline<2,real_t> C_coarse_geom(thb,C_coarse_initial);

    // Refine basis
    thb1 = gsTHBSplineBasis<2,real_t>(bas); // (after refinement) fine basis
    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 0.5,0.5;
    thb1.refine( refBoxes );
    
    // thb.unrefine()
    gsTHBSpline<2,real_t> C_fine_geom = C_coarse_geom;
    C_fine_geom.refineElements({1,0,0,2,2});

    gsQuasiInterpolate<real_t>::localIntpl(thb,C_fine_geom,C_coarse2);

    gsDebugVar(C_coarse2-C_coarse_initial);

}// end main
