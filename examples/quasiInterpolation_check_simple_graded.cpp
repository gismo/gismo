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
    gsMultiPatch<> mpBspline, mp;

    gsTensorBSpline<2,real_t> bspline = *gsNurbsCreator<>::BSplineSquare(1,0,0);
    bspline.degreeElevate(1); // quadratic!!!!! 

    mpBspline.addPatch(bspline);

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.addPatch(thb);
    }


    std::vector<index_t> boxes(5);

    // Initial refinement
    boxes[0] = 1;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = 2;
    boxes[4] = 2;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 2;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = 4;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 3;
    boxes[1] = 2;
    boxes[2] = 0;
    boxes[3] = 8;
    boxes[4] = 6;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 4;
    boxes[1] = 8;
    boxes[2] = 0;
    boxes[3] = 16;
    boxes[4] = 8;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 5;
    boxes[1] = 20;
    boxes[2] = 0;
    boxes[3] = 32;
    boxes[4] = 12;
    mp.patch(0).refineElements(boxes);

    // mp.patch(0).unrefineElements()
    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2); // Function to interpolate 
    gsMultiBasis<> dbasis_fine, dbasis_coarse;


    gsWriteParaview(mp,"init",1,true); //mesh after the refinements

    gsTHBSplineBasis<2,real_t> *thb_graded = dynamic_cast<gsTHBSplineBasis<2,real_t> *>(&mp.basis(0));
    

    dbasis_fine.addBasis(thb_graded->clone());

    std::vector<index_t> boxes_unrefine(5);

    // A VER !!!! NO FUNCIONA LITERALMENTE EL COARSENING... !!!!
    boxes_unrefine[0] = 5;
    boxes_unrefine[1] = 20;
    boxes_unrefine[2] = 0;
    boxes_unrefine[3] = 32;
    boxes_unrefine[4] = 12;
    mp.patch(0).unrefineElements(boxes_unrefine);

    gsWriteParaview(mp,"after",1,true); //mesh after the refinements
    
    gsTHBSplineBasis<2,real_t> *thb_graded_coarse = dynamic_cast<gsTHBSplineBasis<2,real_t> *>(&mp.basis(0));


    dbasis_coarse.addBasis(thb_graded_coarse->clone());

    gsWriteParaview(dbasis_coarse.basis(0),"coarse_thb_graded",1000);
    gsWriteParaview(dbasis_fine.basis(0),"fine_thb_graded",1000);

    gsMatrix<> C_fine_L2, C_fine_QI, C_fine_L2_local;
    gsMatrix<> C_coarse, C_coarse_l2, C_coarse_l2_local;

    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), mySinus, C_fine_QI); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),mySinus,mp,C_fine_L2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0),dbasis_fine.basis(0),mySinus,mp,C_fine_L2_local);
    
    gsGeometry<>::uPtr fine_QI = dbasis_fine.basis(0).makeGeometry(give(C_fine_QI));
    gsGeometry<>::uPtr fine_L2 = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2));
    gsGeometry<>::uPtr fine_L2_local = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2_local));
    
    // gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0), *fine_QI, C_coarse); 
    // gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_coarse, *fine_L2,mp,C_coarse_l2);
    // gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0), dbasis_coarse.basis(0), *fine_L2_local,mp,C_coarse_l2_local);
    
    // gsGeometry<>::uPtr sol_coarse = dbasis_coarse.basis(0).makeGeometry(give(C_coarse));
    // gsGeometry<>::uPtr sol_coarse_L2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2));
    // gsGeometry<>::uPtr sol_coarse_L2_local = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2_local));

    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?

    auto c_sinus_qi_fine = ev.getVariable(*fine_QI,G); // pointer to the QI
    auto c_sinus_L2_fine = ev.getVariable(*fine_L2,G); // pointer to the L2
    auto c_sinus_L2_local_fine = ev.getVariable(*fine_L2_local,G); // pointer to the L2 local
    auto cfunction = ev.getVariable(mySinus,G);

    // auto c_sinus_qi = ev.getVariable(*sol_coarse,G); // pointer to the QI
    // auto c_sinus_L2 = ev.getVariable(*sol_coarse_L2,G); // pointer to the L2 proj
    // auto c_sinus_L2_local = ev.getVariable(*sol_coarse_L2_local,G); // pointer to the LOCAL L2 proj

    gsParaviewCollection error_col("New_error", &ev);
    error_col.options().setSwitch("plotElements", true);
    error_col.options().setInt("plotElements.resolution", 4);
    error_col.options().setInt("numPoints",(mp.geoDim()==3) ? 20000 : 20000);
    error_col.options().setInt("precision", 40); // 1e-18

    error_col.newTimeStep(&mp);
    // error_col.addField((c_sinus_qi_fine-c_sinus_qi).sqNorm(),"error coarsening QI");
    // error_col.addField((c_sinus_L2_fine-c_sinus_L2).sqNorm(),"error coarsening L2");
    // error_col.addField((c_sinus_L2_local_fine-c_sinus_L2_local).sqNorm(),"error coarsening L2 local");
    
    error_col.addField((c_sinus_qi_fine-cfunction).sqNorm(),"error function projection QI");
    error_col.addField((c_sinus_L2_fine-cfunction).sqNorm(),"error function projection L2");
    error_col.addField((c_sinus_L2_local_fine-cfunction).sqNorm(),"error function projection L2 local");

    error_col.saveTimeStep();
    error_col.save();





    // gsDebugVar(mp.basis(0));
    
    /// NO SE PUEDE HACER CON UN RECTANGULO??????????
   
    // gsMultiPatch<> mp;
    // //mp.addPatch(*gsNurbsCreator<>:: BSplineSquare()); // create a multipatch with multibasis!???????
    // // mp.addPatch(*gsNurbsCreator<>:: BSplineRectangle(0,0,1.25,1)); // create a multipatch with multibasis!???????

    // gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2); // Function to interpolate 

    // gsMultiBasis<> dbasis_fine, dbasis_coarse;
    
    // // // Define basis parameters
    // int p = 2;
    // real_t a = 0; // starting knot
    // real_t b = 1; // ending knot
    // unsigned interior = 4; // number of interior knots
    // unsigned multEnd = p+1; // multiplicity at the two end knots

    // gsKnotVector<> kvx(a, b, interior, multEnd);

    // unsigned interior2 = 3; // number of interior knots
    // gsKnotVector<> kvy(a, b, interior2, multEnd);

    // // // Define basis parameters
    // // int p = 2;
    // // real_t a = 0; // starting knot
    // // real_t b = 1; // ending knot
    // // unsigned interior = 4; // number of interior knots
    // // unsigned multEnd = p+1; // multiplicity at the two end knots

    // // gsKnotVector<> kv(a, b, interior, multEnd);
    // gsTensorBSplineBasis<2> bas(kvx,kvy);

    // gsTHBSplineBasis<2,real_t> thb_graded, thb_nongraded;
    // thb_graded = gsTHBSplineBasis<2,real_t>(bas); // create basis with knot vectors
    
    // gsTHBSpline<2,read_t> thb = 
    // mp.addPatch(thb_graded); // create a multipatch with multibasis!???????

    // gsMatrix<> refBoxes(2,2);
    // refBoxes.col(0) << 0.4,0.0;
    // refBoxes.col(1) << 1.0,0.75;
    // thb_graded.refine( refBoxes );

    // gsWriteParaview(thb_graded,"coarse_thb_graded");

    // dbasis_coarse.addBasis(thb_graded.clone());
    
    // // gsMatrix<> refBoxes2(2,2);
    // // refBoxes2.col(0) << 0.7,0.;
    // // refBoxes2.col(1) << 1,0.6;

    // gsMatrix<> refBoxes2(2,2);
    // //refBoxes2.col(0) << 0.7,0.;
    // refBoxes2.col(0) << 0.6,0.;
    // refBoxes2.col(1) << 1,0.5;
    // thb_graded.refine( refBoxes2 );

    // gsWriteParaview(thb_graded,"fine_thb_graded");

    // dbasis_fine.addBasis(thb_graded.clone());

    // gsMatrix<> C_fine_L2, C_fine_QI, C_fine_L2_local;
    // gsMatrix<> C_coarse, C_coarse_l2, C_coarse_l2_local;

    // gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), mySinus, C_fine_QI); 
    // gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),mySinus,mp,C_fine_L2);
    // gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0),dbasis_fine.basis(0),mySinus,mp,C_fine_L2_local);
    
    // gsGeometry<>::uPtr fine_QI = dbasis_fine.basis(0).makeGeometry(give(C_fine_QI));
    // gsGeometry<>::uPtr fine_L2 = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2));
    // gsGeometry<>::uPtr fine_L2_local = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2_local));
    
    // gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0), *fine_QI, C_coarse); 
    // gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0), dbasis_coarse, *fine_L2,mp,C_coarse_l2);
    // gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0), dbasis_coarse.basis(0), *fine_L2_local,mp,C_coarse_l2_local);
    
    // gsGeometry<>::uPtr sol_coarse = dbasis_coarse.basis(0).makeGeometry(give(C_coarse));
    // gsGeometry<>::uPtr sol_coarse_L2 = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2));
    // gsGeometry<>::uPtr sol_coarse_L2_local = dbasis_coarse.basis(0).makeGeometry(give(C_coarse_l2_local));

    // gsExprAssembler<> A(1,1);
    // A.setIntegrationElements(dbasis_fine);
    // gsExprEvaluator<> ev(A);
    // auto G = ev.getMap(mp); // is this correct?

    // auto c_sinus_qi_fine = ev.getVariable(*fine_QI,G); // pointer to the QI
    // auto c_sinus_L2_fine = ev.getVariable(*fine_L2,G); // pointer to the L2
    // auto c_sinus_L2_local_fine = ev.getVariable(*fine_L2_local,G); // pointer to the L2 local

    // auto cfunction = ev.getVariable(mySinus,G);
    // auto c_sinus_qi = ev.getVariable(*sol_coarse,G); // pointer to the QI
    // auto c_sinus_L2 = ev.getVariable(*sol_coarse_L2,G); // pointer to the L2 proj
    // auto c_sinus_L2_local = ev.getVariable(*sol_coarse_L2_local,G); // pointer to the LOCAL L2 proj

    // gsParaviewCollection error_col("New_error", &ev);
    // error_col.options().setSwitch("plotElements", true);
    // error_col.options().setInt("plotElements.resolution", 4);
    // error_col.options().setInt("numPoints",(mp.geoDim()==3) ? 20000 : 20000);
    // error_col.options().setInt("precision", 40); // 1e-18

    // error_col.newTimeStep(&mp);
    // // error_col.addField((c_sinus_qi_fine-c_sinus_qi).sqNorm(),"error coarsening QI");
    // // error_col.addField((c_sinus_L2_fine-c_sinus_L2).sqNorm(),"error coarsening L2");
    // // error_col.addField((c_sinus_L2_local_fine-c_sinus_L2_local).sqNorm(),"error coarsening L2 local");
    
    // // error_col.addField((c_sinus_qi_fine-cfunction).sqNorm(),"error function projection QI");
    // // error_col.addField((c_sinus_L2_fine-cfunction).sqNorm(),"error function projection L2");
    // // error_col.addField((c_sinus_L2_local_fine-cfunction).sqNorm(),"error function projection L2 local");

    // error_col.saveTimeStep();
    // error_col.save();


}// end main
