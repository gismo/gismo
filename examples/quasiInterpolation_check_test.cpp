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
    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>:: BSplineSquare()); // create a multipatch with multibasis!???????
 
    
    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2); // Function to interpolate 

    // Define basis parameters
    int p = 2;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd = p+1; // multiplicity at the two end knots

    gsKnotVector<> kv(a, b, interior, multEnd);
    gsTensorBSplineBasis<2> bas(kv,kv);

    gsTHBSplineBasis<2,real_t> thb, thb_fine;
    thb = gsTHBSplineBasis<2,real_t>(bas); // create basis with knot vectors

    gsMatrix<> coefs, C_coarse, C_fine, C_coarse_final; // for projection

    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 0.5,0.5;
    thb.refine( refBoxes );
    gsMatrix<> refBoxes2(2,2);
    refBoxes2.col(0) << 0,0;
    refBoxes2.col(1) << 0.25,0.25;
    thb.refine( refBoxes2 );
    
    // ===== Get the coefficients of the mySinus function with local interpolation =====
    gsQuasiInterpolate<real_t>::localIntpl(thb, mySinus, coefs); 
    gsMatrix<> C_coarse_initial = coefs;
    // Check for coefficients of the coarse basis
    // gsDebugVar(C_coarse_initial.size()); // gives 22, dofs of the coarse basis

    // ======= Construct the gsTHBSpline with the coarse coefficients of mySinus =======
    gsTHBSpline<2,real_t> thb_coarse_basis(thb,C_coarse_initial); 
    gsWriteParaview(thb_coarse_basis,"funcioninter"); // interpolated function mySinus
    
    // ======= Refine basis and get fine coefficients =======
    gsTHBSpline<2,real_t> thb_fine_basis = thb_coarse_basis;
    thb_fine_basis.refineElements({2,0,0,8,8}); // Refine the spline (basis and the coefficients)
    // gsDebugVar(thb_fine_basis.coefs().size()); // gives 100, dofs of the fine basis
    gsTHBSplineBasis<2,real_t> thb_fine_plot = thb_fine_basis.basis(); // take basis for plotß
    // ======================================================

    // ======= Quasi interpolation II (get the coarse coefficients) =======
    gsQuasiInterpolate<real_t>::localIntpl(thb,thb_fine_basis,C_coarse_final); // I get the coarse coefficients

    // ========= For plotting =========
    gsMultiBasis<> dbasis_fine, dbasis_coarse;
    dbasis_fine.addBasis(thb_fine_plot.clone()); // why?
    //gsDebugVar(thb_fine_basis);
    dbasis_coarse.addBasis(thb.clone());

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?
    ///  ================================
    
    gsGeometry<>::uPtr Ccoarse_final_geom = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_final));
    //gsGeometry<>::uPtr Ccoarse_final_geom = dbasis_fine.basis(0).makeGeometry(give(C_coarse_final));
    gsGeometry<>::uPtr Cfine_geom = dbasis_fine.basis(0).makeGeometry(give(thb_fine_basis.coefs()));
    
    auto cfine = ev.getVariable(*Cfine_geom,G);
    auto ccoarse2 = ev.getVariable(*Ccoarse_final_geom,G);

    A.setIntegrationElements(dbasis_coarse); 
    gsParaviewCollection error_collection("Error_Analysis_QI", &ev);
    error_collection.options().setSwitch("plotElements", true);
    error_collection.options().setInt("plotElements.resolution", 4);
    error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);

    error_collection.newTimeStep(&mp); // i need an updated mp!
    error_collection.addField((ccoarse2-cfine).sqNorm(),"error QI");
    error_collection.saveTimeStep();
    error_collection.save();

}// end main
