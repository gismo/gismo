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
V = verified with test
X = error/incorrect
    [V] inv_expr
    [V] det_expr
    [V] sqNorm_expr
    [V] norm_expr
    [~] col_expr
    [V] gsGeometryMap
    [ ] gsFeElement
    [ ] cdiam_expr
    [V] gsFeVariable
    [V] gsFeSpace
    [V] gsFeSolution
    [V] solGrad_expr
    [V] tr_expr
    [V] temp_expr
    [V] trace_expr
    [V] adjugate_expr
    [V] reshape_expr
    [~] replicate_expr
    [ ] flat_expr
    [V] asDiag_expr
    [V] idMat_expr
    [V] sign_expr
    [V] pow_expr
    [?] matrix_by_space_expr
    [?] matrix_by_space_expr_tr
    [V] value_expr
    [V] grad_expr                           ?? is transposed for space
    [ ] dJacdc_expr
    [X] nabla_expr                          ?? Does not compile
    [X] nabla2_expr                         ?? Does not compile
    [V] onormal_expr
    [V] normalized_expr
    [V] normal_expr
    [V] tangent_expr
    [V] lapl_expr
    [~] solLapl_expr
    [~] solHess_expr                        !! export format is different from hess_expr
    [X] fform_expr
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

    bool verbose = false;
    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("verbose", "Show result and exact", verbose);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo  <<"To do:\n"
            <<"- add deriv2 expressions to gsExpressions.h and make test case\n"
            <<"\n";

    gsMultiPatch<> mp;

    mp.addPatch(gsNurbsCreator<>::NurbsSphere(M_R));
    mp.computeTopology();
    mp.degreeElevate();
    gsMultiBasis<> basis(mp);



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
    gsVector<> point(2);
    point.setConstant(0.5);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////MATRIX OPERATIONS////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"------------------------------Linear Algebra-----------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    /*
        Takes a column of a matrix
        Assessment of:
         - reshape_expr(gsFeVariable)
         - col_expr

        EVALUATION ONLY
    */
    gsInfo<< "* Matrix col(0):\t";
    result = ev.eval(reshape(M,3,3)[0].temp(),point);
    gsInfo<<"passed \t\tnote: eval only\n";
    if (verbose)
        gsInfo<<"Result:\n"<<result<<"\n";

    /*
        Computes trace(M^-1 - I) with M = diag([1,2,3])
        Assessment of:
         - inv_expr(gsFeVariable)
         - trace_expr(gsFeVariable)
         - val_expr(gsFeVariable)
         - idMat_expr
    */
    gsInfo<< "* Matrix expression:\t"; // - gismo::expr::id(3).temp()
    ev.eval(((reshape(M,3,3).inv()-gismo::expr::id(3)).trace()).val(),point);
    gsInfo<<( std::abs(ev.value() - ( 0.0 - 1./2. - 2./3. ) ) < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<ev.value()<<"\n"
                <<"Exact:\n"<<( 0.0 - 1./2. - 2./3. )<<"\n";

    gsInfo<< "* Matrix expr sign:\t"; // - gismo::expr::id(3).temp()
    ev.eval(((reshape(M,3,3).inv()-gismo::expr::id(3)).trace()).sgn(),point);
    gsInfo<<( std::abs(ev.value() + 1) < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<ev.value()<<"\n"
                <<"Exact:\n"<<-1<<"\n";

    /*
        Computes trace(M^-1 - I) with M = diag([1,2,3])
        Assessment of:
         - inv_expr(gsFeVariable)
         - trace_expr(gsFeVariable)
         - val_expr(gsFeVariable)
         - idMat_expr
    */
    gsInfo<< "* Matrix diag:\t\t"; // - gismo::expr::id(3).temp()
    exact  = ev.eval(reshape(M,3,3),point);
    result = ev.eval(m.asDiag(),point);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    /*
        Replicates a vector
        Assessment of:
         - replicate_expr(gsFeVariable)

        EVALUATION ONLY
    */
    gsInfo<< "* Matrix replicate:\t";
    result = ev.eval(replicate(m,1,2),point);
    gsInfo<<"passed \t\tnote: eval only\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n";
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
    exact  = ev.eval(1.0 / reshape(M,3,3).det() * reshape(M,3,3).adj() ,point);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

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
    exact  = ev.eval(pow(reshape(M,3,3).norm(),2),point);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    /*
        Symmetrizes directly (symm) and indirectly (N*N^T)
        Assessment of:
            - reshape_expr(gsFeVariable)
            - symm_expr(gsFeVariable)
            - tr_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix symm:\t\t";
    result = ev.eval(reshape(N,3,3).symm(),point);
    exact  = ev.eval(reshape(N,3,3)*reshape(N,3,3).tr(),point);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    /*
        Symmetrizes directly (symmmetrize) and indirectly (N+N^T)
        Assessment of:
            - reshape_expr(gsFeVariable)
            - symmetrize_expr(gsFeVariable)
            - tr_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix symmetrize:\t";
    result = ev.eval(reshape(N,3,3).symmetrize(),point);
    exact  = ev.eval((reshape(N,3,3)+reshape(N,3,3).tr()),point);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////GEOMETRY OPERATIONS//////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"--------------------------------Geometry---------------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    /*
        Computes the area of the domain and compares to exact (sphere)
        Assessment of:
        - gsExprEvaluator/integral_impl
        - meas_expr
    */
    /// NOTE: Tolerance is lower!
    gsInfo<< "* Area (integral):\t";
    real_t num = ev.integral( meas(G) );
    real_t ref = 4*M_PI*M_R*M_R;
    gsInfo<<( std::abs( num-ref ) / ref < 1e-4 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<num<<"\n"
                <<"Exact:\n"<<ref<<"\n";

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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";


    mp.clear();
    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(M_R,2*M_R));
    mp.degreeElevate();
    basis = gsMultiBasis(mp);

    /*
        Computes the boundary normal vector and compares to exact
        Assessment of:
        - onormal_expr(gsGeometryMap)
        - normalized_expr
    */
    gsVector<real_t,2> resVec,exVec;
    point<<1.0,0.25;
    physpoint = ev.eval( G,point );
    gsInfo<< "* Plane normal:\t\t";
    resVec = ev.eval( nv(G).normalized(), point );
    phi = math::atan2(physpoint(1,0),physpoint(0,0));
    exVec<<math::cos(phi),math::sin(phi);
    gsInfo<<( std::abs( (exVec.transpose()*resVec) ) - 1  < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: sign might be wrong"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    gsInfo<< "* Plane tangent:\t";
    resVec = ev.eval( tv(G).normalized(), point );
    phi = math::atan2(physpoint(1,0),physpoint(0,0));
    exVec<<-math::sin(phi),math::cos(phi);
    gsInfo<<( std::abs( (exVec.transpose()*resVec) ) - 1  < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: sign might be wrong"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";


    /*
        Computes the fundamental form of a geometry
        Assessment of:
        - fform_expr(gsGeometryMap)
    */
    gsInfo<< "* Fundamental form:\t";
    result = ev.eval( fform(G), point );
    exact.resize(2,2);
    phi = math::atan2(physpoint(1,0),physpoint(0,0));
    exact<<1.0,
         -M_R*math::cos(phi)*math::sin(phi) + M_R*math::sin(phi)*math::cos(phi),
         -M_R*math::cos(phi)*math::sin(phi) + M_R*math::sin(phi)*math::cos(phi),
         2*(2*M_R*2*M_R);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tISSUE KNOWN; expression seems wrong"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////FUNCTION OPERATIONS//////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"------------------------------Function operations------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

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
    gsInfo<< "* Gradient:\t\t";
    result = ev.eval( grad(a), point );
    exact.resize(3,2);
    exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    /*
        Computes the Jacobian of a variable
        Assessment of:
        - fjac_expr(gsFeVariable)
    */
    gsInfo<< "* Function Jacobian:\t";
    result = ev.eval( fjac(a), point );
    exact.resize(2,3);
    exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    /*
        Computes the Jacobian of a (1D) variable
        Assessment of:
        - jac_expr(gsFeVariable)
    */
    gsInfo<< "* 1-D Jacobian:\t\t";
    result = ev.eval( jac(b), point );
    exact.resize(1,2);
    exact<<2*physpoint(0,0),2*physpoint(1,0);
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    /*
        Computes the Laplacian of a variable
        Assessment of:
        - lapl_expr(gsFeVariable)
    */
    gsInfo<< "* Laplacian:\t\t";
    result = ev.eval( lapl(b), point );
    exact.resize(1,1);
    exact<<4;
    gsInfo<<( (result-exact).norm() < 1e-5 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n"
                <<"Diff:\n"<< (result-exact) <<"\n";

    /*
        Computes the Hessian of a variable
        Assessment of:
        - hess_expr(gsFeVariable)
    */
    gsInfo<< "* Hessian:\t\t";
    result = ev.eval( hess(b), point );
    exact.resize(2,2);
    exact<<2,0,0,2;
    gsInfo<<( (result-exact).norm() < 1e-5 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n"
                <<"Diff:\n"<< (result-exact) <<"\n";
    /*
        Computes the curl of a variable
        Assessment of:
        - curl_expr(gsFeVariable)
    */
    // gsInfo<< "* curl:\n";
    // result = ev.eval( curl(a), point );
    // gsInfo<< "  Result:\n"<< result <<"\n";
    // exact.resize(3,1);
    // exact<<0,0,1;
    // gsInfo<< "  Exact:\n"<< exact <<"\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////SPACE OPERATIONS/////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"-------------------------------------------------------------------------"<<"\n";;
    gsInfo<<"------------------------------Space & Solution---------------------------"<<"\n";;
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";;

    mp.clear();
    mp.addPatch(gsNurbsCreator<>::BSplineFatQuarterAnnulus(M_R,2*M_R));
    mp.degreeElevate();
    basis = gsMultiBasis(mp);

    space u = A.getSpace(basis,1); // to construct solution manually
    space u2 = A.getSpace(basis,2); // for gsFeSolution
    gsMatrix<> coefs = mp.patch(0).coefs();

    /*
        Computes the values on the geometry

        Performed using
            - gsFeSolution  (basis*coefs inside gsExprEvaluator)
            - gsFeSpace     (basis*coefs after gsExprEvaluator)

        Assessment of:
            - gsFeSpace
            - gsGeometryMap
            - gsFeSolution
    */
    mp.patch(0).eval_into(point,exact);
    exact.transposeInPlace();
    gsInfo<< "* Values (space):\t";
    result.transpose() = ev.eval( u, point );
    result *= coefs;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    gsInfo<< "* Values (map):\t\t";
    result.transpose() = ev.eval( G, point );
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    A.initSystem();
    gsMatrix<> solVec;
    solution u_sol = A.getSolution(u2,solVec);
    gsInfo<< "* Values (solution):\t";
    solVec = coefs;
    solVec.resize(2*solVec.rows(),1);
    result.transpose() = ev.eval( u_sol, point );
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    /*
        Computes the derivatives on the geometry

        Performed using
            - gsFeSolution  (basis*coefs inside gsExprEvaluator)
            - gsFeSpace     (basis*coefs after gsExprEvaluator)

        Assessment of:
            - grad_expr(gsFeSpace)
            - jac_expr(gsGeometryMap)
            - solgrad_expr(gsFeSolution)
    */
    mp.patch(0).deriv_into(point,exact);
    exact.resize(2,2);

    gsInfo<< "* Deriv  (space):\t";
    result.transpose() = ev.eval( grad(u), point );
    result *= coefs;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    gsInfo<< "* Deriv  (map):\t\t";
    result = ev.eval( jac(G), point );
    gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" );//<<"\n";
    // note: transpose() below should not be needed!
    gsInfo<<"\t\tnote: is transposed!"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    gsInfo<< "* Deriv  (solution):\t";
    result = ev.eval( grad(u_sol), point );
    // note: transpose() below should not be needed!
    gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: is transposed!"<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    /*
        Computes the second derivatives on the geometry

        Assessment of:
            - hess_expr(gsGeometryMap)
            - solHess_expr(gsFeSolution)
            - solLapl_expr(gsFeSolution)
    */
    mp.patch(0).deriv2_into(point,exact);

    gsWarn<<"The following expressions have different formats:\n";

    gsMatrix<> result2 = ev.eval( hess(G), point );
    gsInfo<< "* Hess (map): \n"<<result2<<"\n";
    result = ev.eval( shess(u_sol), point );
    gsInfo<< "* Hess (solution): \n"<<result<<"\n";
    result = ev.eval( slapl(u_sol), point );
    gsInfo<< "* Lapl (solution): \n"<<result<<"\n";

    /*
        Computes some expressions for assembly
        NOTE: the expressions are multiplied by the components such that the
        expression with th solution (exact) is represented.

        Performed using
            - gsFeSolution  (basis*coefs inside gsExprEvaluator)
            - gsFeSpace     (basis*coefs after gsExprEvaluator)

        Assessment of:
            - transpose_expr(gsFeSolution)
            - transpose_expr(gsFeSpace)
            - solGrad_expr(gsFeSolution)
            - jac_expr(gsFeSpace)
            - mult_expr (type 1 & type 2)
    */

    gsInfo<<"* s grad(u):\t\t";
    exact.transpose() = ev.eval( u_sol.tr() * grad(u_sol), point );
    result = ev.eval( (u_sol.tr() * jac(u2)).tr(), point );
    result *= solVec;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact.transpose()<<"\n";

    gsInfo<<"* u grad(s):\t\t";
    result.transpose() = ev.eval( u2 * grad(u_sol), point );
    result *= solVec;
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact.transpose()<<"\n";


    // COLBLOCKS ERROR
    gsInfo<<"* s grad(u) + u grad(s):\t";
    result = ev.eval( u2*grad(u_sol) + u_sol * jac(u2) , point );
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<2*exact.transpose()<<"\n";



    return EXIT_SUCCESS;
}
