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

    gsMatrix<> C_coarse, C_fine; // for projection
    
    dbasis_coarse.addBasis(thb.clone());

    // ===== Get the coefficients of the mySinus function with local interpolation =====
    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 1,1;
    thb.refine( refBoxes );
    
    dbasis_fine.addBasis(thb.clone());

    //gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0), mySinus, C_fine); 
    gsL2Projection<real_t>::projectFunction(dbasis_fine.basis(0),mySinus,mp,C_fine);

    
    gsTHBSpline<2,real_t> thb_fine_basis(thb,C_fine); 

    gsWriteParaview(thb_fine_basis,"funcioninter"); // interpolated function mySinus


    //gsDebugVar(thb_fine_basis.uniformCoarsen());
    // thb_fine_basis.uniformCoarsen(1);
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0), thb_fine_basis, C_coarse); 

    // gsDebugVar(C_fine.size());
    // gsDebugVar(C_coarse.size());


    // ====== Plotting error from QI ======
    // gsMultiBasis<> dbasis_qi;
    // dbasis_qi.addBasis(thb.clone()); // why?

    gsGeometry<>::uPtr sol_coarse = dbasis_coarse.basis(0).makeGeometry(give(C_coarse));
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_fine);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?
    
    auto c_sinus_qi = ev.getVariable(*sol_coarse,G); // pointer to the QI
    auto cfunction = ev.getVariable(mySinus,G);

    gsParaviewCollection error_col("Error_Simple", &ev);
    error_col.options().setSwitch("plotElements", true);
    error_col.options().setInt("plotElements.resolution", 4);
    error_col.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 10000);
    error_col.options().setInt("precision", 32); // 1e-18

    error_col.newTimeStep(&mp);
    error_col.addField((cfunction-c_sinus_qi).sqNorm(),"error mySinus QI");

    error_col.saveTimeStep();
    error_col.save();


}// end main
