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
    unsigned interior = 0; // number of interior knots
    unsigned multEnd = p+1; // multiplicity at the two end knots

    gsKnotVector<> kv(a, b, interior, multEnd);
    gsTensorBSplineBasis<2> bas(kv,kv);

    gsTHBSplineBasis<2,real_t> thb, thb_fine;
    thb = gsTHBSplineBasis<2,real_t>(bas); // create basis with knot vectors


    gsMultiBasis<> dbasis_fine, dbasis_coarse;

    gsMatrix<> C_coarse, C_fine, C_coarse_l2, C_coarse_l2_local; // for projection
    
    dbasis_coarse.addBasis(thb.clone());

    //gsWriteParaview(thb, "hola");


    // ===== Get the coefficients of the mySinus function with local interpolation =====
    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 1,1;
    thb.refine( refBoxes );
    
    dbasis_fine.addBasis(thb.clone());

    gsMatrix<> C_fine_QI, C_fine_L2, C_fine_L2_local;

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), mySinus, C_fine_QI); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),mySinus,mp,C_fine_L2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0),dbasis_fine.basis(0),mySinus,mp,C_fine_L2_local);

    gsTHBSpline<2,real_t> thb_fine_basis(thb,C_fine_QI); 

    gsGeometry<>::uPtr fine_QI = dbasis_fine.basis(0).makeGeometry(give(C_fine_QI));
    gsGeometry<>::uPtr fine_L2 = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2));
    gsGeometry<>::uPtr fine_L2_local = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2_local));

    gsWriteParaview(thb_fine_basis,"Surface L2 initial QI"); // interpolated function mySinus

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0), *fine_QI, C_coarse); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_coarse, *fine_L2,mp,C_coarse_l2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0), dbasis_coarse.basis(0), *fine_L2_local,mp,C_coarse_l2_local);

    thb.unrefine(refBoxes,0); // coarse mesh!
    gsTHBSpline<2,real_t> thb_fine_plot(thb,C_coarse); 
    gsTHBSpline<2,real_t> thb_fine_plot_l2(thb,C_coarse_l2); 
    gsTHBSpline<2,real_t> thb_fine_plot_l2_local(thb,C_coarse_l2_local); 

    gsWriteParaview(thb_fine_plot,"Surface QI"); // interpolated function mySinus
    gsWriteParaview(thb_fine_plot_l2,"Surface L2"); // interpolated function mySinus
    gsWriteParaview(thb_fine_plot_l2_local,"Surface L2 local"); // interpolated function mySinus

    // ====== Plotting error from QI ======xw
    gsGeometry<>::uPtr sol_coarse = dbasis_coarse.basis(0).makeGeometry(give(C_coarse));
    gsGeometry<>::uPtr sol_coarse_L2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2));
    gsGeometry<>::uPtr sol_coarse_L2_local = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2_local));

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?
    
    auto c_sinus_L2_local = ev.getVariable(*sol_coarse_L2_local,G); // pointer to the LOCAL L2 proj
    auto c_sinus_L2 = ev.getVariable(*sol_coarse_L2,G); // pointer to the L2 proj
    auto c_sinus_qi = ev.getVariable(*sol_coarse,G); // pointer to the QI
    auto cfunction = ev.getVariable(mySinus,G);

    gsParaviewCollection error_col("Error_Simple", &ev);
    error_col.options().setSwitch("plotElements", true);
    error_col.options().setInt("plotElements.resolution", 4);
    error_col.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 10000);
    error_col.options().setInt("precision", 32); // 1e-18

    error_col.newTimeStep(&mp);
    error_col.addField((cfunction-c_sinus_qi).sqNorm(),"error mySinus QI");
    error_col.addField((cfunction-c_sinus_L2).sqNorm(),"error mySinus L2");
    error_col.addField((cfunction-c_sinus_L2_local).sqNorm(),"error mySinus L2 local");

    error_col.saveTimeStep();
    error_col.save();

}// end main
