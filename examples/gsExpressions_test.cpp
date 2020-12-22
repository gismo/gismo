/** @file gsExprIntegral_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsAssembler/gsExprEvaluator.h>

using namespace gismo;

int main(int argc, char *argv[])
{

/*
~ = needs a check but seems correct
? = expression unclear
V = verified correct
X = error/incorrect
    [V] inv_expr
    [V] det_expr
    [V] sqNorm_expr
    [V] norm_expr
    [ ] col_expr
    [ ] _expr
    [V] gsGeometryMap
    [ ] gsFeElement
    [ ] cdiam_expr
    [V] gsFeVariable
    [V] gsFeSpace
    [V] gsFeSolution
    [V] solGrad_expr
    [V] tr_expr
    [ ] temp_expr
    [ ] trace_expr
    [V] adjugate_expr
    [V] reshape_expr
    [ ] replicate_expr
    [ ] flat_expr
    [~] asDiag_expr
    [~] idMat_expr
    [~] sign_expr
    [V] pow_expr
    [?] matrix_by_space_expr
    [?] matrix_by_space_expr_tr
    [ ] value_expr
    [V] grad_expr                           ?? is transposed for space
    [ ] dJacdc_expr
    [?] nabla_expr                          ?? Does not compile
    [?] nabla2_expr                         ?? Does not compile
    [V] onormal_expr
    [V] normalized_expr
    [V] normal_expr
    [ ] tangent_expr
    [V] lapl_expr
    [~] solLapl_expr
    [~] solHess_expr                        !! export format is different from hess_expr
    [X] fform_expr                          ?? Commented
    [ ] jacGinv_expr
    [ ] jacG_expr
    [V] jac_expr
    [V] fjac_expr                           ?? Is transposed
    [V] hess_expr                           ?? Should be 1D
    [ ] dJacG_expr
    [V] meas_expr
    [X] curl_expr                           ?? gives bus error
    [V] mult_expr
    [ ] collapse_expr
    [ ] frprod_expr
    [ ] divide_expr
    [V] add_expr
    [ ] summ_expr
    [V] sub_expr
    [V] symm_expr
    [V] symmetrize_expr
 */

    # define M_PI 3.14159265358979323846
    # define M_R  1.0

    gsMultiPatch<> mp;

    mp.addPatch(gsNurbsCreator<>::NurbsSphere(M_R));
    gsWriteParaview(mp,"mp",1000,1,1);

    mp.computeTopology();
    mp.degreeElevate();
    gsInfo<<mp.patch(0).basis();
    gsMultiBasis<> basis(mp);

    gsWriteParaview(basis.basis(0),"basis");


    //b.basis(0).component(0).uniformRefine();

    gsFunctionExpr<> m_("1","2","3",2);
    gsFunctionExpr<> M_("1","0","0","0","2","0","0","0","3",2);
    gsFunctionExpr<> N_("1","0","1","0","1","0","0","0","1",2);

    // Set the expression assembler
    gsExprAssembler<> A(1,1);
    // Set the parameter mesh as the integration mesh
    A.setIntegrationElements(basis);

    // Set the expression evaluator
    gsExprEvaluator<> ev(A);

    // Define integrant variables
    typedef gsExprAssembler<>::element     element;
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    element     e = ev.getElement();
    geometryMap G = ev.getMap(mp);
    variable    m = ev.getVariable(m_);
    variable    M = ev.getVariable(M_);
    variable    N = ev.getVariable(N_);

    gsMatrix<> result;
    gsMatrix<> exact;
    gsMatrix<> physpoint;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////MATRIX OPERATIONS////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    gsVector<> point(2);
    point.setConstant(0.5);
    gsInfo<< "* Vector:\n";
    ev.eval(m,point);
    gsInfo<< "  Global: "<< ev.allValues().transpose() <<"\n";

    gsInfo<< "* Matrix:\n";
    result = ev.eval(reshape(M,3,3),point);
    gsInfo<< result <<"\n";

    gsInfo<< "* Matrix col(0):\n";
    result = ev.eval(reshape(M,3,3)[0],point);
    gsInfo<<result<<"\n";

    gsInfo<< "* Matrix sign(trace(M - I)):\n"; // - gismo::expr::id(3).temp()
    result = ev.eval(((reshape(M,3,3).inv()-gismo::expr::id(3)).trace()).sgn(),point);
    gsInfo<< result <<"\n";

    gsInfo<< "* Matrix replicate:\n";
    result = ev.eval(replicate(m,1,2),point);
    gsInfo<< result <<"\n";


    /*
        Computes the inverse in two ways: 1) directly with Eigen, 2) with 1/det(M)*adj(M).
        Assessment of:
            - reshape_expr(gsFeVariable)
            - inv_expr(gsFeVariable)
            - det_expr(gsFeVariable)
            - adj_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix inverse:\t";
    result = ev.eval(reshape(M,3,3).inv(),point);
    result -= ev.eval(1.0 / reshape(M,3,3).det() * reshape(M,3,3).adj() ,point);
    gsInfo<<( result.norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
        Computes the norm in two ways: 1) sqnorm, 2) power of the norm.
        Assessment of:
            - reshape_expr(gsFeVariable)
            - sqNorm_expr(gsFeVariable)
            - norm_expr(gsFeVariable)
            - pow_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix sqnorm:\t";
    result = ev.eval(reshape(M,3,3).sqNorm(),point);
    result -= ev.eval(pow(reshape(M,3,3).norm(),2),point);
    gsInfo<<( result.norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
        Symmetrizes directly (symm) and indirectly (N*N^T)
        Assessment of:
            - reshape_expr(gsFeVariable)
            - symm_expr(gsFeVariable)
            - tr_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix symm:\t\t";
    result = ev.eval(reshape(N,3,3).symm(),point);
    result -= ev.eval(reshape(N,3,3)*reshape(N,3,3).tr(),point);
    gsInfo<<( result.norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
        Symmetrizes directly (symmmetrize) and indirectly (N+N^T)
        Assessment of:
            - reshape_expr(gsFeVariable)
            - symmetrize_expr(gsFeVariable)
            - tr_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix symmetrize:\t";
    result = ev.eval(reshape(N,3,3).symmetrize(),point);
    result -= ev.eval((reshape(N,3,3)+reshape(N,3,3).tr()),point);
    gsInfo<<( result.norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////GEOMETRY OPERATIONS//////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*
        Computes the area of the domain and compares to exact (sphere)
        Assessment of:
        - gsExprEvaluator/integral_impl
        - meas_expr
    */
    gsInfo<< "* Domain area (integral):\t";
    real_t num = ev.integral( meas(G) );
    real_t ref = 4*M_PI*M_R*M_R;
    gsInfo<<( std::abs( num-ref ) / ref < 1e-10 ? "passed" : "failed" )<<"\n";

    point.setConstant(0.5);
    physpoint = ev.eval( G,point );
    real_t phi = math::atan2(physpoint(1,0),physpoint(0,0));
    real_t theta = math::atan2(math::sqrt( math::pow(physpoint(1,0),2)+math::pow(physpoint(0,0),2) ), physpoint(2,0));
    // real_t x = M_R * math::sin(theta)*math::cos(phi);
    // real_t y = M_R * math::sin(theta)*math::sin(phi);
    // real_t z = M_R * math::cos(phi);

    gsInfo<< "* Normal vector:\t";
    result = ev.eval( sn(G).normalized(), point );
    exact.resize(3,1);
    exact<<math::cos(theta)*math::sin(phi),math::sin(theta)*math::sin(phi),math::cos(phi);
    exact.normalize();
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    // gsInfo<< "  Result: "<< result <<"\n";
    // gsInfo<< "  Exact:  "<< exact <<"\n";


    mp.clear();
    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(M_R,2*M_R));
    mp.degreeElevate();
    gsWriteParaview(mp,"mp",1000,1,1);
    basis = gsMultiBasis(mp);
    gsWriteParaview(basis.basis(0),"basis");

    /*
        Computes the boundary normal vector and compares to exact
        Assessment of:
        - onormal_expr(gsGeometryMap)
        - normalized_expr
    */
    point<<1.0,0.5;
    physpoint = ev.eval( G,point );
    gsInfo<< "* Outward Normal vector:\n";
    gsVector<real_t,2> resVec = ev.eval( nv(G).normalized(), point );
    gsVector<real_t,2> exVec;
    exVec<<math::sin(phi),math::cos(phi);
    // result = gsVector<real_t,2>();
    phi = math::atan2(physpoint(1,0),physpoint(0,0));
    exact.resize(2,1);
    exact<<math::sin(phi),math::cos(phi);
    gsInfo<<( std::abs( (exVec.transpose()*resVec) - 1 ) < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: "<< ev.allValues().transpose() <<"\n";
    gsInfo<< "  Exact:  "<< exact.transpose()/exact.norm() <<"\n";

    /*
        Computes the fundamental form of a geometry
        Assessment of:
        - fform_expr(gsGeometryMap)
    */
    gsInfo<< "* fform:\n";
    result = ev.eval( fform(G), point );
    gsInfo<< "  Result: \n"<< result <<"\n";
    exact.resize(2,2);
    phi = math::atan2(physpoint(1,0),physpoint(0,0));
    exact<<1.0,
         -M_R*math::cos(phi)*math::sin(phi) + M_R*math::sin(phi)*math::cos(phi),
         -M_R*math::cos(phi)*math::sin(phi) + M_R*math::sin(phi)*math::cos(phi),
         2*M_R*M_R;
    gsInfo<< "  Exact:  \n"<< exact <<"\n";


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////FUNCTION OPERATIONS//////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsFunctionExpr<> a_("x^2","y^2","x*y",2);
    gsFunctionExpr<> b_("x^2+y^2",2);
    gsFunctionExpr<> c_("-y","x","0",2);
    variable    a = ev.getVariable(a_, G);
    variable    b = ev.getVariable(b_, G);
    variable    c = ev.getVariable(c_, G);

    /*
        Computes the gradient of a variable
        Assessment of:
        - grad_expr(gsFeVariable)
    */
    gsInfo<< "* Grad:\t\t";
    result = ev.eval( grad(a), point );
    exact.resize(3,2);
    exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
        Computes the Jacobian of a variable
        Assessment of:
        - fjac_expr(gsFeVariable)
    */
    gsInfo<< "* fjac:\t\t";
    result = ev.eval( fjac(a), point );
    exact.resize(2,3);
    exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
        Computes the Laplacian of a variable
        Assessment of:
        - lapl_expr(gsFeVariable)
    */
    gsInfo<< "* Lapl:\t\t";
    result = ev.eval( lapl(b), point );
    exact.resize(1,1);
    exact<<4;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: "<< result <<"\n";
    gsInfo<< "  Exact:  "<< exact <<"\n";
    gsInfo<< "  Diff:  "<< (result-exact) <<"\n";

    /*
        Computes the Hessian of a variable
        Assessment of:
        - hess_expr(gsFeVariable)
    */
    gsInfo<< "* Hess:\t\t";
    result = ev.eval( hess(b), point );
    exact.resize(2,2);
    exact<<2,0,0,2;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: "<< result <<"\n";
    gsInfo<< "  Exact:  "<< exact <<"\n";
    gsInfo<< "  Diff:  "<< (result-exact) <<"\n";

    /*
        Computes the Jacobian of a (1D) variable
        Assessment of:
        - jac_expr(gsFeVariable)
    */
    gsInfo<< "* jac:\t\t";
    result = ev.eval( jac(b), point );
    exact.resize(1,2);
    exact<<2*physpoint(0,0),2*physpoint(1,0);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
        Computes the curl of a variable
        Assessment of:
        - curl_expr(gsFeVariable)
    */
    // gsInfo<< "* curl:\n";
    // result = ev.eval( curl(a), point );
    // gsInfo<< "  Result: "<< result <<"\n";
    // exact.resize(3,1);
    // exact<<0,0,1;
    // gsInfo<< "  Exact:  "<< exact <<"\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////SPACE OPERATIONS/////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    space u = A.getSpace(basis,1); // to construct solution manually
    space u2 = A.getSpace(basis,2); // for gsFeSolution
    gsMatrix<> coefs = mp.patch(0).coefs();

    /*
        Computes the values on the geometry
        Assessment of:
            - gsFeSpace
            - gsGeometryMap
            - gsFeSolution
    */
    mp.patch(0).eval_into(point,exact);
    exact.transposeInPlace();
    gsInfo<< "* values (space):\t";
    result.transpose() = ev.eval( u, point );
    result *= coefs;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: "<<result<<"\n";
    gsInfo<< "  Exact:  "<< exact <<"\n";

    gsInfo<< "* values (map):\t";
    result.transpose() = ev.eval( G, point );
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: "<<result<<"\n";
    gsInfo<< "  Exact:  "<< exact <<"\n";

    A.initSystem();
    gsMatrix<> solVec;
    solution u_sol = A.getSolution(u2,solVec);
    gsInfo<< "* values (solution):\t";
    solVec = coefs;
    solVec.resize(2*solVec.rows(),1);
    result.transpose() = ev.eval( u_sol, point );
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: "<<result<<"\n";
    gsInfo<< "  Exact:  "<< exact <<"\n";

    /*
        Computes the derivatives on the geometry
        Assessment of:
            - grad_expr(gsFeSpace)
            - jac_expr(gsGeometryMap)
            - solgrad_expr(gsFeSolution)
    */
    mp.patch(0).deriv_into(point,exact);
    exact.resize(2,2);

    gsInfo<< "* deriv (space):\t";
    result.transpose() = ev.eval( grad(u), point );
    result *= coefs;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: \n"<<result<<"\n";
    gsInfo<< "  Exact:  \n"<< exact <<"\n";

    gsInfo<< "* deriv (map) [TRANSPOSED]:\t";
    result = ev.eval( jac(G), point );
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: \n"<<result<<"\n";
    gsInfo<< "  Exact:  \n"<< exact <<"\n";

    gsInfo<< "* deriv (solution) [TRANSPOSED]:\t";
    result = ev.eval( grad(u_sol), point );
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result: \n"<<result<<"\n";
    gsInfo<< "  Exact:  \n"<< exact <<"\n";

    /*
        Computes the second derivatives on the geometry
        Assessment of:
            - hess_expr(gsGeometryMap)
            - solHess_expr(gsFeSolution)
            - solLapl_expr(gsFeSolution)
    */
    mp.patch(0).deriv2_into(point,exact);
    gsDebugVar(exact);

    gsWarn<<"The following expressions have different formats:\n";

    gsMatrix<> result2 = ev.eval( hess(G), point );
    gsInfo<< "  Hess (map): \n"<<result2<<"\n";
    result = ev.eval( shess(u_sol), point );
    gsInfo<< "  Hess (solution): \n"<<result<<"\n";
    result = ev.eval( slapl(u_sol), point );
    gsInfo<< "  Lapl (solution): \n"<<result<<"\n";

    /*
        Computes some expressions for assembly
        Assessment of:
            - transpose_expr(gsFeSolution)
            - transpose_expr(gsFeSpace)
            - solGrad_expr(gsFeSolution)
            - jac_expr(gsFeSpace)
            - mult_expr (type 1 & type 2)
    */

    gsInfo<<"* u grad(u) [sol*jac(space)]:\t";
    exact.transpose() = ev.eval( u_sol.tr() * grad(u_sol), point );
    result = ev.eval( (u_sol.tr() * jac(u2)).tr(), point );
    result *= solVec;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result:"<<result.transpose()<<"\n"; //
    gsInfo<< "  Exact: "<<exact.transpose()<<"\n";

    gsInfo<<"* u grad(u) [space*grad(sol)]:\t";
    result.transpose() = ev.eval( u2 * grad(u_sol), point );
    result *= solVec;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<< "  Result:"<<result.transpose()<<"\n"; //
    gsInfo<< "  Exact: "<<exact.transpose()<<"\n";

    // // COLBLOCKS ERROR
    // result = ev.eval( u2*grad(u_sol) + u_sol * jac(u2) , point );
    // gsInfo<< "  Result (something): \n"<<result<<"\n";

    return EXIT_SUCCESS;
}
