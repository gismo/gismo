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
    //gsTensorBSpline<2,real_t> bspline = *gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.25);
    gsTensorBSpline<2,real_t> bspline = *gsNurbsCreator<>::BSplineSquare(1,0,0);

    bspline.degreeElevate(1); // quadratic!!!!! 

    int numInitUniformRefine  = 2;  
    for (int i = 0; i < numInitUniformRefine; ++i)
        bspline.uniformRefine();

    mpBspline.addPatch(bspline);
    
    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.addPatch(thb);
    }

    // gsAdaptiveMeshing<2,real_t> mesher(mp);
    // gsHBoxContainer<2,real_t> hboxes;
    // gsHBoxContainer<2,real_t> refine,coarsen;

    std::vector<index_t> boxes(5);

    boxes[0] = 1;
    boxes[1] = 2;
    boxes[2] = 0;
    boxes[3] = 8;
    boxes[4] = 6;
    mp.patch(0).refineElements(boxes);

    // gsHBox<2,real_t> hbox({1,2,0,8,6},mp.basis(0));
    // mp.patch(0).refineElements(hbox.toRefBox());

    boxes[0] = 2;
    boxes[1] = 8;
    boxes[2] = 0;
    boxes[3] = 16;
    boxes[4] = 8;
    mp.patch(0).refineElements(boxes);
    boxes[0] = 3;
    boxes[1] = 20;
    boxes[2] = 0;
    boxes[3] = 32;
    boxes[4] = 12;
    mp.patch(0).refineElements(boxes);

    gsWriteParaview(mp,"init",1,true); //mesh after the refinements
    
    // Pointer to initial refined basis
    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));


    gsMultiBasis<> dbasis_fine, dbasis_coarse;
    // dbasis_fine.addBasis(thb.clone());

    gsTHBSplineBasis<2,real_t> *thb_graded = dynamic_cast<gsTHBSplineBasis<2,real_t> *>(&mp.basis(0));
    dbasis_fine.addBasis(thb_graded->clone());
    gsDebugVar(dbasis_fine.basis(0));
    gsWriteParaview(*thb_graded,"basis_graded",1000); // this is the fine mesh

    // // Approximations
    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2); // Function to interpolate 
    
    // // does not work
    // std::vector<index_t> boxes_unrefine(5);
    // boxes_unrefine[0] = 5;
    // boxes_unrefine[1] = 20;
    // boxes_unrefine[2] = 0;
    // boxes_unrefine[3] = 32;
    // boxes_unrefine[4] = 12;
    // mp.patch(0).unrefineElements(boxes_unrefine);
    // gsWriteParaview(mp,"after_coarsening",1,true); //mesh after the refinements

    gsMatrix<> C_fine_L2, C_fine_QI, C_fine_L2_local;
    gsMatrix<> C_coarse, C_coarse_l2, C_coarse_l2_local;
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), mySinus, C_fine_QI); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),mySinus,mp,C_fine_L2);
    gsQuasiInterpolate<real_t>::localL2(dbasis_fine.basis(0),mySinus,C_fine_L2_local);

    // gsDebugVar(C_fine_L2_local-C_fine_L2); // 1e-6
    // gsDebugVar(C_fine_L2_local-C_fine_QI); // 1e-15
    
    gsGeometry<>::uPtr fine_QI = dbasis_fine.basis(0).makeGeometry(give(C_fine_QI));
    gsGeometry<>::uPtr fine_L2 = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2));
    gsGeometry<>::uPtr fine_L2_local = dbasis_fine.basis(0).makeGeometry(give(C_fine_L2_local));
    
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?
    auto c_sinus_qi_fine = ev.getVariable(*fine_QI); // pointer to the QI
    auto c_sinus_L2_fine = ev.getVariable(*fine_L2); // pointer to the L2
    auto c_sinus_L2_local_fine = ev.getVariable(*fine_L2_local); // pointer to the L2 local
    auto cfunction = ev.getVariable(mySinus);
   
    gsInfo<<"hola\n";

    gsParaviewCollection error_col("New_error", &ev);
    error_col.options().setSwitch("plotElements", true);
    error_col.options().setInt("plotElements.resolution", 4);
    error_col.options().setInt("numPoints",(mp.geoDim()==3) ? 20000 : 20000);
    error_col.options().setInt("precision", 40); // 1e-18
    error_col.newTimeStep(&mp);
    error_col.addField((c_sinus_qi_fine-cfunction).sqNorm(),"error function projection QI");
    error_col.addField((c_sinus_L2_fine-cfunction).sqNorm(),"error function projection L2");
    error_col.addField((c_sinus_L2_local_fine-cfunction).sqNorm(),"error function projection L2 local");
    error_col.saveTimeStep();
    error_col.save();


    // ================ Adaptive meshing ================
    gsAdaptiveMeshing<2,real_t> mesher(mp);
    index_t rule = 3;

    mesher.options().setInt("RefineRule",rule);
    mesher.options().setInt("CoarsenRule",rule);
    mesher.options().setSwitch("Admissible",true); // true for admissible refinement
    if (rule==1)
    {
        mesher.options().setReal("RefineParam",0.3);
        mesher.options().setReal("CoarsenParam",0.1);
    }
    else if (rule==2)
    {
        mesher.options().setReal("RefineParam",0.3);
        mesher.options().setReal("CoarsenParam",0.1);
    }
    else if (rule==3)
    {
        mesher.options().setReal("RefineParam",0.1);
        mesher.options().setReal("CoarsenParam",0.01);
    }

    mesher.getOptions();
    gsHBoxContainer<2,real_t> refine,coarsen;
    // mesher.markRef_into(errors,refine);

    gsVector<index_t,2> low,upp;
    index_t lvl;
    low <<2,5;
    upp <<3,6;
    lvl = 1;
    gsHBox<2>cell(low,upp,lvl,basis);
    
    gsVector<index_t,2> low1,upp1;
    index_t lvl1;
    low1 <<3,3;
    upp1 <<4,4;
    lvl1 = 0;
    gsHBox<2>cell2(low1,upp1,lvl1,basis);

    gsDebugVar(basis);

    // mesher.markRef_into(errors,refine);
    refine.add(cell);
    refine.add(cell2);

    gsWriteParaview(refine,"marked4ref2");

    mesher.refine(refine);
    gsWriteParaview(mp,"end_mp",1,true);

}// end main