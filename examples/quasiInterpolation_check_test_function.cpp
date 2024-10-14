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
 
    
    //gsFunctionExpr<real_t> myfunction("2/(3*exp(sqrt((10*x-3)^2+(10*y-3)^2))) + 2/(3*exp(sqrt((10*x+3)^2+(10*y+3)^2))) + 2/(3*exp(sqrt((10*x)^2+(10*y)^2))) ",2); // Function to interpolate 
    //gsFunctionExpr<real_t> myfunction("2/(3*exp(sqrt((10*(2*x-1)-3)^2+(10*(2*y-1)-3)^2))) + 2/(3*exp(sqrt((10*(2*x-1)+3)^2+(10*(2*y-1)+3)^2))) + 2/(3*exp(sqrt((10*(2*x-1))^2+(10*(2*y-1))^2)))",2);
    //gsFunctionExpr<real_t> myfunction("2/(3*exp(sqrt((20*(x-0.75))^2+(20*(y-0.75))^2))) + 2/(3*exp(sqrt((20*(x-0.25))^2+(20*(y-0.25))^2))) + 2/(3*exp(sqrt((20*(x-0.5))^2+(20*(y-0.5))^2)))",2);
    //gsFunctionExpr<real_t> myfunction("2/(3*exp(sqrt((10*x - 3)^2+(10*y - 3)^2))) + 2/(3*exp(sqrt((10*x + 3)^2+(10*y+3)^2))) + 2/(3*exp(sqrt((10*x)^2+(10*y)^2))) ",2); // Function to interpolate 
    //gsFunctionExpr<real_t> myfunction( "2/(3*exp(sqrt((10*((x + 1)/2) - 3)^2 + (10*((y + 1)/2) - 3)^2))) + 2/(3*exp(sqrt((10*((x + 1)/2) + 3)^2 + (10*((y + 1)/2) + 3)^2))) + 2/(3*exp(sqrt((10*((x + 1)/2))^2 + (10*((y + 1)/2))^2)))",2);
    
    gsFunctionExpr<real_t> myfunction(
        "2/(3*exp(sqrt((10*(2*x - 1) - 3)^2 + (10*(2*y - 1) - 3)^2))) + \
        2/(3*exp(sqrt((10*(2*x - 1) + 3)^2 + (10*(2*y - 1) + 3)^2))) + \
        2/(3*exp(sqrt((10*(2*x - 1))^2 + (10*(2*y - 1))^2)))",
        2);

    // gsFunctionExpr<real_t> myfunction(
    //     "exp(-1/(x*(1-x)))*exp(-1/(y*(1-y)))",
    //     2);
        
    // Define basis parameters
    int p = 2;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 149; // number of interior knots
    unsigned multEnd = p+1; // multiplicity at the two end knots

    gsKnotVector<> kv(a, b, interior, multEnd);
    gsTensorBSplineBasis<2> bas(kv,kv);

    gsTHBSplineBasis<2,real_t> thb, thb_fine;
    thb = gsTHBSplineBasis<2,real_t>(bas); // create basis with knot vectors

    gsMatrix<> coefs, C_coarse, C_fine, C_coarse_final, C_fine_coefs, coefs_in; // for projection
    gsMatrix<> coefs_L2, C_fine_coefs_L2, C_coarse_final_L2, C_fine_coefs_L2_local, C_coarse_final_L2_local, coefs_L2_local; // for projection


    // ===== Get the coefficients of the mySinus function with local interpolation =====
    gsQuasiInterpolate<real_t>::localIntpl(thb, myfunction, coefs); 
    gsL2Projection<real_t>::projectFunction(thb,myfunction,mp,coefs_L2);
    gsQuasiInterpolate<real_t>::localL2(thb,myfunction,coefs_L2_local);
   
    // ====== Plotting error from QI ======
    gsMultiBasis<> dbasis_qi;
    dbasis_qi.addBasis(thb.clone()); // why?

    gsTHBSpline<2,real_t> thb_coarse_basis_L2(thb,coefs_L2); 
    gsTHBSpline<2,real_t> thb_coarse_basis_QI(thb,coefs); 
    gsTHBSpline<2,real_t> thb_coarse_basis_L2_local(thb,coefs_L2_local); 

    // gsWriteParaview(myfunction,"sinus");
    gsWriteParaview(thb_coarse_basis_L2,"function_L2", 100000); 
    gsWriteParaview(thb_coarse_basis_QI,"function_QI",100000); 
    gsWriteParaview(thb_coarse_basis_L2_local,"function_L2_local",100000); 

    // gsWriteParaview(myfunction,"function_original"); 

    gsGeometry<>::uPtr fun_qi = dbasis_qi.basis(0).makeGeometry(give(coefs));
    gsGeometry<>::uPtr fun_L2 = dbasis_qi.basis(0).makeGeometry(give(coefs_L2));
    gsGeometry<>::uPtr fun_L2_local = dbasis_qi.basis(0).makeGeometry(give(coefs_L2_local));

    // gsWriteParaview(*fun_L2_local,"function_L2_local",100000); 


    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis_qi);
    gsExprEvaluator<> ev(A);
    auto G = ev.getMap(mp); // is this correct?

    auto c_fun_L2 = ev.getVariable(*fun_L2,G); // pointer to the L2
    auto c_fun_qi = ev.getVariable(*fun_qi,G); // pointer to the QI
    auto c_fun_L2_local = ev.getVariable(*fun_L2_local,G); // pointer to the L2 local
    auto cfunction = ev.getVariable(myfunction,G);

    gsParaviewCollection error_collection("Error_Analysis_Projection", &ev);
    error_collection.options().setSwitch("plotElements", true);
    error_collection.options().setInt("plotElements.resolution", 4);
    error_collection.options().setInt("numPoints", 100000);
    error_collection.options().setInt("precision", 40); // 1e-32

    error_collection.newTimeStep(&mp); // i need an updated mp!
    error_collection.addField((cfunction-c_fun_qi).sqNorm(),"error function QI");
    error_collection.addField((cfunction-c_fun_L2).sqNorm(),"error function L2");
    error_collection.addField((cfunction-c_fun_L2_local).sqNorm(),"error function L2 local");

    error_collection.addField(cfunction,"funcion");

    error_collection.saveTimeStep();
    error_collection.save();

}// end main
