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

    gsMatrix<> coefs, C_coarse, C_fine, C_coarse_final, C_fine_coefs, coefs_in, coefs_final; // for projection

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
    coefs_in =  coefs;
    // Check for coefficients of the coarse basis
    // gsDebugVar(C_coarse_initial.size()); // gives 22, dofs of the coarse basis
    
    // ====== Plotting error from QI ======
    gsMultiBasis<> dbasis_qi;
    dbasis_qi.addBasis(thb.clone()); // why?

    gsGeometry<>::uPtr sinus_qi = dbasis_qi.basis(0).makeGeometry(give(coefs));
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_qi);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?

    auto c_sinus_initial = ev.getVariable(*sinus_qi,G); // pointer to the QI
    auto cfunction = ev.getVariable(mySinus,G);

    gsParaviewCollection error_collection("Error_Analysis_QI", &ev);
    error_collection.options().setSwitch("plotElements", true);
    error_collection.options().setInt("plotElements.resolution", 4);
    error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 10000);
    error_collection.options().setInt("precision", 30); // 1e-18

    error_collection.newTimeStep(&mp); // i need an updated mp!
    error_collection.addField((cfunction-c_sinus_initial).sqNorm(),"error mySinus");
    // ====================================

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
    coefs_final = C_coarse_final;

    // gsDebugVar((C_coarse_final-C_coarse_initial).norm());
    // gsDebugVar(C_coarse_final-C_coarse_initial);

    // ========= For plotting =========
    gsMultiBasis<> dbasis_fine, dbasis_coarse;
    dbasis_fine.addBasis(thb_fine_plot.clone()); // why?
    //gsDebugVar(thb_fine_basis);
    dbasis_coarse.addBasis(thb.clone());


    // ====== !!!!Error de QI a fine mesh!!!! ======
    // sinus_qi is the pointer to the initial coarse coefficients

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0),*sinus_qi,C_fine_coefs); // I get the coarse coefficients 
    gsGeometry<>::uPtr Cfine_geom = dbasis_fine.basis(0).makeGeometry(give(C_fine_coefs));
    A.setIntegrationElements(dbasis_fine);
    auto c_fine = ev.getVariable(*Cfine_geom,G);
    auto cfunction_2 = ev.getVariable(mySinus,G);
    error_collection.addField((c_sinus_initial-c_fine).sqNorm(),"error_refinement");

    // ======== Plot error from coarsening =======
    gsGeometry<>::uPtr Ccoarse2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_final)); // C_coarse_final    
    // Also works
    // gsGeometry<>::uPtr Ccoarse2 = thb.makeGeometry(give(C_coarse_final)); // C_coarse_final    
    A.setIntegrationElements(dbasis_coarse); // no lo sé, tengo que ver como hacer?????
    auto c_coarse_2 = ev.getVariable(*Ccoarse2,G);
    auto c_fine2= ev.getVariable(*Cfine_geom,G);
    error_collection.addField((c_fine-c_coarse_2).sqNorm(),"error_coarsening");
    error_collection.saveTimeStep();
    error_collection.save();


    // gsExprAssembler<> A(1,1);
    // A.setIntegrationElements(dbasis_fine);
    // gsExprEvaluator<> ev(A);
    // auto G = ev.getMap(mp); // is this correct?
    ///  ================================
    
    // gsGeometry<>::uPtr Ccoarse_final_geom = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_final));
    // //gsGeometry<>::uPtr Ccoarse_final_geom = dbasis_fine.basis(0).makeGeometry(give(C_coarse_final));
    // gsGeometry<>::uPtr Cfine_geom = dbasis_fine.basis(0).makeGeometry(give(thb_fine_basis.coefs()));
    
    // auto cfine = ev.getVariable(*Cfine_geom,G);
    // auto ccoarse2 = ev.getVariable(*Ccoarse_final_geom,G);

    // A.setIntegrationElements(dbasis_coarse); 
    // gsParaviewCollection error_collection("Error_Analysis_QI", &ev);
    // error_collection.options().setSwitch("plotElements", true);
    // error_collection.options().setInt("plotElements.resolution", 4);
    // error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);

    // error_collection.newTimeStep(&mp); // i need an updated mp!
    // error_collection.addField((ccoarse2-cfine).sqNorm(),"error QI");
    // error_collection.saveTimeStep();
    // error_collection.save();

}// end main
