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

    gsMatrix<> coefs, C_coarse, C_fine, C_coarse_final, C_fine_coefs, coefs_in; // for projection
    gsMatrix<> coefs_L2, C_fine_coefs_L2, C_coarse_final_L2, C_fine_coefs_L2_local, C_coarse_final_L2_local, coefs_L2_local; // for projection

    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 0.5,0.5;
    thb.refine( refBoxes );
    gsMatrix<> refBoxes2(2,2);
    refBoxes2.col(0) << 0,0;
    refBoxes2.col(1) << 0.25,0.25;
    thb.refine( refBoxes2 );
    gsMatrix<> refBoxes3(2,2);
    refBoxes3.col(0) << 0.5,0.5;
    refBoxes3.col(1) << 1.0,1.0;
    thb.refine( refBoxes3 );
    gsMatrix<> refBoxes4(2,2);
    refBoxes4.col(0) << 0.75,0.5;
    refBoxes4.col(1) << 1.0,0.75;
    thb.refine( refBoxes4 );
    gsMatrix<> refBoxes5(2,2);
    refBoxes5.col(0) << 0.875,0.5;
    refBoxes5.col(1) << 1.,0.625;
    thb.refine( refBoxes5 );
    // gsMatrix<> refBoxes3(2,2);
    // refBoxes3.col(0) << 0.5,0.0;
    // refBoxes3.col(1) << 1.0,0.5;
    // thb.refine( refBoxes3 );
    // gsMatrix<> refBoxes4(2,2);
    // refBoxes4.col(0) << 0.75,0.0;
    // refBoxes4.col(1) << 1.0,0.25;
    // thb.refine( refBoxes4 );
    // gsMatrix<> refBoxes5(2,2);
    // refBoxes5.col(0) << 0.875,0.0;
    // refBoxes5.col(1) << 1.,0.125;
    // thb.refine( refBoxes5 );
    
    // gsMatrix<> refBoxes4(2,2);
    // refBoxes4.col(0) << 0.5,0.0;
    // refBoxes4.col(1) << 1.0,0.5;
    // thb.refine( refBoxes4 );
    // gsMatrix<> refBoxes5(2,2);
    // refBoxes5.col(0) << 0.75,0.0;
    // refBoxes5.col(1) << 1.,0.25;
    // thb.refine( refBoxes5 );


    // ===== Get the coefficients of the mySinus function with local interpolation =====
    gsQuasiInterpolate<real_t>::localIntpl(thb, mySinus, coefs); 
    gsL2Projection<real_t>::projectFunction(thb,mySinus,mp,coefs_L2);
    gsQuasiInterpolate<real_t>::localL2(thb,mySinus,coefs_L2_local);

    gsMatrix<> C_coarse_initial = coefs;
    gsMatrix<> C_coarse_L2_initial = coefs_L2;

    coefs_in =  coefs;
    // Check for coefficients of the coarse basis
    // gsDebugVar(C_coarse_initial.size()); // gives 22, dofs of the coarse basis
    
    // ====== Plotting error from QI ======
    gsMultiBasis<> dbasis_qi;
    dbasis_qi.addBasis(thb.clone()); // why?

    gsGeometry<>::uPtr sinus_qi = dbasis_qi.basis(0).makeGeometry(give(coefs));
    gsGeometry<>::uPtr sinus_L2 = dbasis_qi.basis(0).makeGeometry(give(coefs_L2));
    gsGeometry<>::uPtr sinus_L2_local = dbasis_qi.basis(0).makeGeometry(give(coefs_L2_local));

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_qi);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?
    
    auto c_sinus_L2 = ev.getVariable(*sinus_L2,G); // pointer to the L2
    auto c_sinus_qi = ev.getVariable(*sinus_qi,G); // pointer to the QI
    auto c_sinus_L2_local = ev.getVariable(*sinus_L2_local,G); // pointer to the L2 local
    auto cfunction = ev.getVariable(mySinus,G);

    gsParaviewCollection error_collection("Error_Analysis_QI", &ev);
    error_collection.options().setSwitch("plotElements", true);
    error_collection.options().setInt("plotElements.resolution", 4);
    error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 10000);
    error_collection.options().setInt("precision", 32); // 1e-32

    error_collection.newTimeStep(&mp); // i need an updated mp!
    error_collection.addField((cfunction-c_sinus_qi).sqNorm(),"error mySinus QI");
    error_collection.addField((cfunction-c_sinus_L2).sqNorm(),"error mySinus L2");
    error_collection.addField((cfunction-c_sinus_L2_local).sqNorm(),"error mySinus L2 local");

    // ====================================

    // ======= Construct the gsTHBSpline with the coarse coefficients of mySinus =======
    gsTHBSpline<2,real_t> thb_coarse_basis(thb,C_coarse_initial); 
    gsWriteParaview(thb_coarse_basis,"funcioninter"); // interpolated function mySinus
    
    // ======= Refine basis and get fine coefficients =======
    gsTHBSpline<2,real_t> thb_fine_basis = thb_coarse_basis;
    thb_fine_basis.refineElements({3,0,0,16,16}); // Refine the spline (basis and the coefficients)
    //thb_fine_basis.refineElements({3,8,0,8,8}); // Refine the spline (basis and the coefficients)
    //thb_fine_basis.refineElements({3,12,8,16,12}); // Refine the spline (basis and the coefficients)
    //thb_fine_basis.refineElements({3,12,8,12,16}); // Refine the spline (basis and the coefficients)
    gsTHBSplineBasis<2,real_t> thb_fine_plot = thb_fine_basis.basis(); // take basis for plot
    // ======================================================
   
    // ========= For plotting =========
    gsMultiBasis<> dbasis_fine, dbasis_coarse;
    dbasis_fine.addBasis(thb_fine_plot.clone()); // why?
    dbasis_coarse.addBasis(thb.clone());

    gsMatrix<> coefs_trial = thb_fine_basis.coefs();
    gsGeometry<>::uPtr geom_fine_coefs = dbasis_fine.basis(0).makeGeometry((coefs_trial));

    // // ======= Quasi interpolation II (get the coarse coefficients) =======
    //gsQuasiInterpolate<real_t>::localIntpl(thb,thb_fine_basis,C_coarse_final); // I get the coarse coefficients 
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0),*geom_fine_coefs,C_coarse_final); // I get the coarse coefficients 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),dbasis_coarse,*geom_fine_coefs,mp,C_coarse_final_L2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_coarse.basis(0),*geom_fine_coefs,C_coarse_final_L2_local);

    // ====== !!!!Error de QI a fine mesh!!!! ======
    // sinus_qi is the pointer to the initial coarse coefficients
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0),*sinus_qi,C_fine_coefs); // I get the fine coefficients 
    gsL2Projection<real_t>::projectFunction(dbasis_fine,*sinus_L2,mp,C_fine_coefs_L2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0),*sinus_L2_local,C_fine_coefs_L2_local);

    gsGeometry<>::uPtr Cfine_geom = dbasis_fine.basis(0).makeGeometry(give(C_fine_coefs));
    gsGeometry<>::uPtr Cfine_geom_L2 = dbasis_fine.basis(0).makeGeometry(give(C_fine_coefs_L2));
    gsGeometry<>::uPtr Cfine_geom_L2_local = dbasis_fine.basis(0).makeGeometry(give(C_fine_coefs_L2_local));

    A.setIntegrationElements(dbasis_fine);
    
    auto c_fine = ev.getVariable(*Cfine_geom,G);
    auto c_fine_L2 = ev.getVariable(*Cfine_geom_L2,G);
    auto c_fine_L2_local = ev.getVariable(*Cfine_geom_L2_local,G);

    auto cfunction_2 = ev.getVariable(mySinus,G);
    error_collection.addField((c_sinus_qi-c_fine).sqNorm(),"error_refinement QI");
    error_collection.addField((c_sinus_L2-c_fine_L2).sqNorm(),"error_refinement L2");
    error_collection.addField((c_sinus_L2_local-c_fine_L2_local).sqNorm(),"error_refinement L2 local");


    gsTHBSpline<2,real_t> thb_coarse_basis2(thb,C_coarse_final); 
    gsTHBSpline<2,real_t> thb_coarse_basis2_L2(thb,C_coarse_final_L2); 

    gsWriteParaview(thb_coarse_basis2,"funcioninter_final"); // interpolated function mySinus
    gsWriteParaview(thb_coarse_basis2_L2,"funcioninter_final_l2"); // interpolated function mySinus

    // // ======== Plot error from coarsening =======
    gsGeometry<>::uPtr Ccoarse2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_final)); // coarse QI  
    gsGeometry<>::uPtr Ccoarse2_L2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_final_L2)); // coarse L2  
    gsGeometry<>::uPtr Ccoarse2_L2_local = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_final_L2_local)); // coarse L2 local

    // Also works
    // gsGeometry<>::uPtr Ccoarse2 = thb.makeGeometry(give(C_coarse_final)); // C_coarse_final    
    A.setIntegrationElements(dbasis_coarse); 
    auto c_coarse_2 = ev.getVariable(*Ccoarse2,G);
    // auto c_fine2= ev.getVariable(*Cfine_geom,G);
    auto c_coarse_2_L2 = ev.getVariable(*Ccoarse2_L2,G);
    auto c_coarse_2_L2_local = ev.getVariable(*Ccoarse2_L2_local,G);


    error_collection.addField((c_fine-c_coarse_2_L2).sqNorm(),"error_coarsening L2");
    error_collection.addField((c_fine-c_coarse_2).sqNorm(),"error_coarsening QI");
    error_collection.addField((c_fine-c_coarse_2_L2_local).sqNorm(),"error_coarsening L2 local");
    
    error_collection.addField((c_coarse_2_L2-c_sinus_qi).sqNorm(),"error_final L2");
    error_collection.addField((c_coarse_2-c_sinus_qi).sqNorm(),"error_final QI");
    error_collection.addField((c_coarse_2_L2_local-c_sinus_qi).sqNorm(),"error_final L2 local");


    error_collection.saveTimeStep();
    error_collection.save();

}// end main
