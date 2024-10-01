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
 
    
    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2);
    //gsFunctionExpr<> myPolyQuad("-5*x^2 + 3*x*y + y^2 + 22*x + y - 4",2);

    int p = 2;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd = p+1; // multiplicity at the two end knots

    gsKnotVector<> kv(a, b, interior, multEnd);

    gsTensorBSplineBasis<2> bas(kv,kv);
    // gsAdaptiveMeshing<2,real_t> mesher; //????
    // gsMultiBasis<> dbasis;
    // dbasis.addBasis(gsTHBSplineBasis<2> (bas(kv,kv)));


    int numRef = 5;
    bool passed = true;

    gsTHBSplineBasis<2,real_t> thb, thb1, thb2, thb3, thb_fine, thb_coarse;
 
    thb = gsTHBSplineBasis<2,real_t>(bas); // coarse basis

    thb_coarse = thb; 

    gsMatrix<> coefs, C_coarse, C_fine, C_coarse_final; // for projection

    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 0.5,0.5;
    thb.refine( refBoxes );
    gsMatrix<> refBoxes2(2,2);
    refBoxes2.col(0) << 0,0;
    refBoxes2.col(1) << 0.25,0.25;
    thb.refine( refBoxes2 );

    // gsDebugVar(thb);
    // gsWriteParaview(thb,"holaprueba");
    
    // Get COARSE coefficients    
    gsQuasiInterpolate<real_t>::localIntpl(thb, mySinus, coefs); 
    gsMatrix<> C_coarse_initial = coefs;
    // construct the thb spline basis with the coarse coefficients
    gsTHBSpline<2,real_t> thb_coarse_basis(thb,C_coarse_initial); 
    gsWriteParaview(thb_coarse_basis,"funcioninter"); // interpolated function  sin cos

    // gsWriteParaview(thb,"QUE"); // coarse basis
    
    // thb.unrefine()
    
    // ======= Refine basis and get fine coefficients=======
    gsTHBSpline<2,real_t> thb_fine_basis = thb_coarse_basis;
    thb_fine_basis.refineElements({2,0,0,8,8}); // how to plot this?


    gsDebugVar(thb_fine_basis.coefs().size());
    
    gsTHBSplineBasis<2,real_t> thb_fine_plot = thb_fine_basis.basis();

    //gsDebugVar(thb_fine_basis.coefs());
    // gsDebugVar(thb_fine_basis);

    // Quasi interpolation -- Fine coefficients!
    gsQuasiInterpolate<real_t>::localIntpl(thb,thb_fine_basis,C_coarse_final);

    //gsDebugVar(C_coarse_final);


    //gsDebugVar(C_coarse_final-C_coarse_initial);
    
    ///  ========= For plotting =========
    gsMultiBasis<> dbasis_fine, dbasis_coarse;
    dbasis_fine.addBasis(thb_fine_plot.clone()); // why?
    //gsDebugVar(thb_fine_basis);
    dbasis_coarse.addBasis(thb.clone());

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?
    ///  ================================
    
    // gsDebugVar(C_coarse_final);
    // gsDebugVar(thb_fine_basis.coefs());

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

    // // gsDebugVar(ev.integralElWise(meas(G) * (ccoarse2)));
    // gsDebugVar(ev.integralElWise(meas(G) * (cfine)));


    // gsDebugVar(C_coarse_final-C_coarse_initial);

    //gsGeometry<>::uPtr C_last = thb_coarse.basis(0).makeGeometry(give(C_coarse_final)); //coarse

   
   
   
   
   
    // auto cfine = ev.getVariable(thb_fine_basis,G);
    // auto ccoarse = ev.getVariable(*C_last,G); //quasi interpolation

    // A.setIntegrationElements(dbasis_coarse); // to do plot in the coarse mesh!
    // gsParaviewCollection error_collection("Error_Analysis", &ev);
    // error_collection.newTimeStep(&mp); // i need an updated mp!
    // error_collection.addField((ccoarse-cfine).sqNorm(),"error QI");

}// end main
