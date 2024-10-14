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
    

    //gsWriteParaview(thb, "hola");
    //gsDebugVar(thb);

    // ===== Get the coefficients of the mySinus function with local interpolation =====
    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 1,1;
    thb.refine( refBoxes );

    dbasis_coarse.addBasis(thb.clone());


    gsMatrix<> refBoxes1(2,2);
    refBoxes1.col(0) << 0,0;
    refBoxes1.col(1) << 0.5,0.5;
    thb.refine( refBoxes1 );



    // gsMatrix<> refBoxes2(2,2);
    // refBoxes1.col(0) << 0,0;
    // refBoxes1.col(1) << 0.25,0.25;
    // thb.refine( refBoxes2 );
    
    dbasis_fine.addBasis(thb.clone());

    // gsInfo<<"hola\n";
    // gsDebugVar(thb);

    gsMatrix<> C_fine_QI, C_fine_L2, C_fine_L2_local;

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), mySinus, C_fine_QI); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),mySinus,mp,C_fine_L2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0),mySinus,C_fine_L2_local);
    
    gsTHBSpline<2,real_t> thb_fine_basis(thb,C_fine_QI); 
    gsTHBSpline<2,real_t> thb_fine_basis_L2(thb,C_fine_L2); 
    gsTHBSpline<2,real_t> thb_fine_basis_L2_local(thb,C_fine_L2_local); 

    //gsWriteParaview(thb,"hola_fine");

    gsGeometry<>::uPtr fine_QI = dbasis_fine.basis(0).makeGeometry(give(C_fine_QI));
    gsGeometry<>::uPtr fine_L2 = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2));
    gsGeometry<>::uPtr fine_L2_local = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2_local));
    
    gsWriteParaview(thb_fine_basis,"Surface initial QI"); // interpolated function mySinus
    gsWriteParaview(thb_fine_basis_L2,"Surface L2 initial"); // interpolated function mySinus
    gsWriteParaview(thb_fine_basis_L2_local,"Surface L2 local initial"); // interpolated function mySinus

    gsInfo<< "================ FIND  ================\n";

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0), *fine_QI, C_coarse); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_coarse, *fine_L2,mp,C_coarse_l2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_coarse.basis(0), *fine_L2_local,C_coarse_l2_local);
    
    // gsDebugVar(dbasis_coarse.basis(0));
    // gsDebugVar(C_coarse_l2_local-C_coarse);
    // gsDebugVar(C_coarse_l2_local-C_coarse_l2);

    thb.unrefine(refBoxes1,0); // coarse mesh!
    gsTHBSpline<2,real_t> thb_coarse_plot(thb,C_coarse); 
    gsTHBSpline<2,real_t> thb_coarse_plot_l2(thb,C_coarse_l2); 
    gsTHBSpline<2,real_t> thb_coarse_plot_l2_local(thb,C_coarse_l2_local); 

    gsWriteParaview(thb_coarse_plot,"Surface_QI"); // interpolated function mySinus
    gsWriteParaview(thb_coarse_plot_l2,"Surface_L2"); // interpolated function mySinus
    gsWriteParaview(thb_coarse_plot_l2_local,"Surface_L2_local"); // interpolated function mySinus

    // ====== Plotting error from QI ======
    gsGeometry<>::uPtr sol_coarse = dbasis_coarse.basis(0).makeGeometry(give(C_coarse));
    gsGeometry<>::uPtr sol_coarse_L2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2));
    gsGeometry<>::uPtr sol_coarse_L2_local = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2_local));

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_coarse);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?

    auto c_sinus_qi_fine = ev.getVariable(*fine_QI,G); // pointer to the QI
    auto c_sinus_L2_fine = ev.getVariable(*fine_L2,G); // pointer to the L2
    auto c_sinus_L2_local_fine = ev.getVariable(*fine_L2_local,G); // pointer to the L2 local

    auto cfunction = ev.getVariable(mySinus,G);
    auto c_sinus_qi = ev.getVariable(*sol_coarse,G); // pointer to the QI
    auto c_sinus_L2 = ev.getVariable(*sol_coarse_L2,G); // pointer to the L2 proj
    auto c_sinus_L2_local = ev.getVariable(*sol_coarse_L2_local,G); // pointer to the LOCAL L2 proj


    gsParaviewCollection error_col("Error_Simple", &ev);
    error_col.options().setSwitch("plotElements", true);
    error_col.options().setInt("plotElements.resolution", 4);
    error_col.options().setInt("numPoints",(mp.geoDim()==3) ? 20000 : 20000);
    error_col.options().setInt("precision", 40); // 1e-18

    error_col.newTimeStep(&mp);
    error_col.addField((c_sinus_qi_fine-c_sinus_qi).sqNorm(),"error coarsening QI");
    error_col.addField((c_sinus_L2_fine-c_sinus_L2).sqNorm(),"error coarsening L2");
    error_col.addField((c_sinus_L2_local_fine-c_sinus_L2_local).sqNorm(),"error coarsening L2 local");

    error_col.addField((c_sinus_qi_fine-cfunction).sqNorm(),"error function projection QI");
    error_col.addField((c_sinus_L2_fine-cfunction).sqNorm(),"error function projection L2");
    error_col.addField((c_sinus_L2_local_fine-cfunction).sqNorm(),"error function projection L2 local");

    error_col.addField(c_sinus_L2_local_fine-c_sinus_qi_fine,"diff fine local and QI");
    error_col.addField(c_sinus_L2_local-c_sinus_qi,"diff coarse local and QI");
    error_col.addField(c_sinus_L2-c_sinus_qi,"diff coarse L2 and QI");

    error_col.saveTimeStep();
    error_col.save();

}// end main
