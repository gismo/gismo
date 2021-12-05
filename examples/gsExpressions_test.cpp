/** @file gsExpressions_test.cpp

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
  [V] col_expr
  [V] gsGeometryMap
  -->[ ] gsFeElement                         ?? don't know how to test this
  -->[ ] cdiam_expr                          ?? don't know how to test this
  [V] gsFeVariable
  [V] gsFeSpace
  [V] gsFeSolution
  [V] solGrad_expr
  [V] tr_expr
  [V] temp_expr
  [V] trace_expr
  [V] adjugate_expr
  [V] reshape_expr
  [V] replicate_expr
  [V] flat_expr                           note: this is only implemented for 2x2 matrices!
  [V] asDiag_expr
  [V] idMat_expr
  [V] sign_expr
  [V] pow_expr
  -->[???] matrix_by_space_expr           ?? does not give a result.. What should it do?
  -->[???] matrix_by_space_expr_tr        ?? Same as above
  [V] value_expr
  [V] grad_expr                           ?? is transposed for space
  -->[???] dJacdc_expr                       ?? How to check this one?
  -->[X] nabla_expr                          !! Does not compile
  -->[X] nabla2_expr                         !! Does not compile
  [V] onormal_expr
  [V] normalized_expr
  [V] normal_expr
  [V] tangent_expr
  [V] lapl_expr                           !! Tolerance of 1e-5 needed
  -->[V] lapl_expr(gsFeSolution)
  -->[X] fform_expr                          !! expression is wrong?
  [V] jac_expr                            !! there seems to be a problem when using jac(gsFeSpace).cols()
  [V] fjac_expr                           ?? Is transposed
  [V] hess_expr                           !! Tolerance of 1e-5 needed
  [V] hess_expr (gsGeometryMap)
  -->[X] hess_expr (gsFeSolution)         !! does not compile!               !! export format is different from hess_expr(G)
  -->[X] dJacG_expr                       !! does not compile!
  [V] meas_expr
  -->[X] curl_expr                        !! does not compile!
  [V] mult_expr
  -->[ ] collapse_expr                    ?? what does it do?
  -->[X] frprod_expr (v1)                 !! sizes seem to be wrong, maybe also cardinality in trace is wrong
  -->[V] frprod_expr (v2)
  [V] divide_expr
  [V] add_expr
  -->[V] summ_expr                        ?? works, but how to verify best?
  [V] sub_expr
  [V] symm_expr
  [V] symmetrize_expr
*/

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
    gsFunctionExpr<> O_("1","2","3","4",2);
    gsFunctionExpr<> o_("1","4","5",2);

    // Set the expression assembler
    gsExprAssembler<> A(1,1);
    // Set the parameter mesh as the integration mesh
    A.setIntegrationElements(basis);

    // Set the expression evaluator
    gsExprEvaluator<> ev(A);

    // Define integrant variables
    //element     e = ev.getElement();
    auto G = ev.getMap(mp);
    auto m = ev.getVariable(m_);
    auto M = ev.getVariable(M_);
    auto N = ev.getVariable(N_);
    auto O = ev.getVariable(O_);
    auto o = ev.getVariable(o_);

    //auto el= ev.getElement();

    gsMatrix<> result, exact, tmp;
    gsVector<> physpoint, point(2);
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
    */
    gsInfo<< "* Matrix col(0):\t";
    result = ev.eval(reshape(M,3,3)[0].temp(),point);
    exact  = ev.eval(reshape(M,3,3),point);
    if (verbose)
        gsInfo<<"Result:\n"<<result<<"\n"
              <<"Exact:\n"<<exact.col(0)<<"\n";

    gsInfo<<( (result-exact.col(0)).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    // gsDebug<<ev.eval(el.diam().val(),point);

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
    if (verbose)
        gsInfo  <<"Result:\n"<<ev.value()<<"\n"
                <<"Exact:\n"<<( 0.0 - 1./2. - 2./3. )<<"\n";
    gsInfo<<( std::abs(ev.value() - ( 0.0 - 1./2. - 2./3. ) ) < 1e-10 ? "passed" : "failed" )<<"\n";

    gsInfo<< "* Matrix expr sign:\t"; // - gismo::expr::id(3).temp()
    ev.eval(((reshape(M,3,3).inv()-gismo::expr::id(3)).trace()).sgn(),point);
    if (verbose)
        gsInfo  <<"Result:\n"<<ev.value()<<"\n"
                <<"Exact:\n"<<-1<<"\n";
    gsInfo<<( std::abs(ev.value() + 1) < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      XXXX
      Assessment of:
      - XXX
    */
    gsInfo<< "* Matrix diag:\t\t"; // - gismo::expr::id(3).temp()
    exact  = ev.eval(reshape(M,3,3),point);
    result = ev.eval(m.asDiag(),point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Replicates a vector
      Assessment of:
      - replicate_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix replicate:\t";
    result = ev.eval(replicate(m,1,3),point);
    exact = ev.eval(m,point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result.diagonal()<<"\n"
                <<"Exact:\n"<<exact<<"\n";

    gsInfo<<( (result.diagonal()-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Symmetrizes directly (symmetrize) and indirectly (N+N^T)
      Assessment of:
      - reshape_expr(gsFeVariable)
      - symmetrize_expr(gsFeVariable)
      - tr_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix symmetrize:\t";
    result = ev.eval(reshape(N,3,3).symmetrize(),point);
    exact  = ev.eval((reshape(N,3,3)+reshape(N,3,3).tr()),point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    /*
      Flattens a 2x2 matrix [a,b,c,d] to a 3x1 vector [a,b,c+d]
      - flat_expr(gsFeVariable)
    */
    gsInfo<< "* Matrix flatten:\t";
    result = ev.eval(flat(reshape(O,2,2)),point);
    exact  = ev.eval(o,point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

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
    real_t ref = 4*EIGEN_PI*M_R*M_R;
    if (verbose)
        gsInfo  <<"Result:\n"<<num<<"\n"
                <<"Exact:\n"<<ref<<"\n";
    gsInfo<<( std::abs( num-ref ) / ref < 1e-4 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";

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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    mp.clear();
    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(M_R,2*M_R));
    mp.degreeElevate();
    basis = gsMultiBasis<>(mp);

    gsWriteParaview(mp,"mp");

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
    if (verbose)
        gsInfo  <<"Result:\n"<<resVec<<"\n"
                <<"Exact:\n"<<exVec<<"\n";
    gsInfo<<( std::abs( (exVec.transpose()*resVec) ) - 1  < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: sign might be wrong"<<"\n";

    gsInfo<< "* Plane tangent:\t";
    resVec = ev.eval( tv(G).normalized(), point );
    phi = math::atan2(physpoint(1,0),physpoint(0,0));
    exVec<<-math::sin(phi),math::cos(phi);
    if (verbose)
        gsInfo  <<"Result:\n"<<resVec<<"\n"
                <<"Exact:\n"<<exVec<<"\n";
    gsInfo<<( std::abs( (exVec.transpose()*resVec) ) - 1  < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: sign might be wrong"<<"\n";

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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tISSUE KNOWN; expression seems wrong"<<"\n";


    gsInfo<< "* Hess (map): \t";
    result = ev.eval( hess(G), point );
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n";
    gsInfo<<"computed\n";//<<"\n";


    // THIS ONE DOES NOT COMPILE
/*
    gsInfo<< "* dJac (map): \t";
    result = ev.eval( dJac(G), point );
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n";
    gsInfo<<"computed\n";//<<"\n";
*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////FUNCTION OPERATIONS//////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"------------------------------Function operations------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    gsFunctionExpr<> a_("x^2","y^2","x*y",2);// R^2 -> R
    gsFunctionExpr<> b_("x^2+y^2",2);
    gsFunctionExpr<> c_("-y","x","0",2);
    auto a = ev.getVariable(a_, G);
    auto b = ev.getVariable(b_, G);
    //auto c = ev.getVariable(c_, G);

    /*
      Computes the value of a variable
      Assessment of:
      - gsFeVariable
      = avg(gsFeVariable)
    */
    gsInfo<< "* Value:\t\t";
    result = ev.eval( a, point );
    exact  = a_.eval(physpoint);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    /*
    gsInfo<< "* Average:\t\t";
    result.resize(1,1); result.at(0) = ev.integralInterface( avg(a) );
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    */

    /*
      Computes the gradient of a variable
      Assessment of:
      - grad_expr(gsFeVariable)
    */
    gsInfo<< "* Gradient:\t\t";
    result = ev.eval( grad(a), point );
    exact.resize(3,2);
    exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Computes the Jacobian of a variable
      Assessment of:
      - fjac_expr(gsFeVariable)
    */
    gsInfo<< "* Function Jacobian:\t";
    result = ev.eval( jac(a), point );
    //exact.resize(3,2); //done above
    //exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
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
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n"
                <<"Diff:\n"<< (result-exact) <<"\n";
    gsInfo<<( (result-exact).norm() < 1e-5 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";

    /*
      Computes the Hessian of a variable
      Assessment of:
      - hess_expr(gsFeVariable)
    */
    gsInfo<< "* Hessian:\t\t";
    result = ev.eval( hess(b), point );
    exact.resize(2,2);
    exact<<2,0,0,2;
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n"
                <<"Diff:\n"<< (result-exact) <<"\n";
    gsInfo<<( (result-exact).norm() < 1e-5 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";

    // DOES NOT COMPILE
    /*
      Computes the curl of a variable
      Assessment of:
      - curl_expr(gsFeVariable)
    */
    /*
      gsInfo<< "* curl:\n";
      result = ev.eval( curl(a), point );
      gsInfo<< "  Result:\n"<< result <<"\n";
      exact.resize(3,1);
      exact<<0,0,1;
      gsInfo<< "  Exact:\n"<< exact <<"\n";
    */

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////SPACE OPERATIONS/////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"------------------------------Space & Solution---------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    mp.clear();
    mp.addPatch(gsNurbsCreator<>::BSplineFatQuarterAnnulus(M_R,2*M_R));
//    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(M_R,2*M_R));
    mp.degreeElevate();
    basis = gsMultiBasis<>(mp);

    auto u = A.getSpace(basis, 1); // to construct solution manually
    auto u2 = A.getSpace(basis,2); // for gsFeSolution

    u .setup(gsBoundaryConditions<>(), 0, 0);
    u2.setup(gsBoundaryConditions<>(), 0, 0);

    gsMatrix<> coefs = mp.patch(0).coefs();


    gsInfo<<"------------------------------Sandbox---------------------------"<<"\n";

    /*
      Computes the matrix using a space
      Assessment of:
      - matrix_by_space_expr(gsFeVariable,space)
    */
    // WRONG: DOES NOT GIVE RESULT
    gsInfo<< "* Matrix by space:\t";
    result = ev.eval(matrix_by_space(N,u2),point);
    // exact  = ev.eval(o,point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n";
                // <<"Exact:\n"<<exact<<"\n";
    // gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    gsInfo<< "* dJacdc:\t";
    result = ev.eval(dJacdc(a,1),point);
    // exact  = ev.eval(o,point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n";
                // <<"Exact:\n"<<exact<<"\n";
    // gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    // THESE TWO DO NOT COMPILE

    // gsInfo<< "* Nabla:\t";
    // result = ev.eval(nabla(a),point);
    // // exact  = ev.eval(o,point);
    // if (verbose)
    //     gsInfo  <<"Result:\n"<<result<<"\n";
    //             // <<"Exact:\n"<<exact<<"\n";
    // // gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    // gsInfo<< "* Nabla:\t";
    // result = ev.eval(nabla2(m),point);
    // // exact  = ev.eval(o,point);
    // if (verbose)
    //     gsInfo  <<"Result:\n"<<result<<"\n";
    //             // <<"Exact:\n"<<exact<<"\n";
    // // gsInfo<<( (result.transpose()-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    gsInfo<<"------------------------------Sandbox---------------------------"<<"\n";


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
    gsDebugVar(result);
    result *= coefs;
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    gsInfo<< "* Values (map):\t\t";
    result.transpose() = ev.eval( G, point );
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    A.initSystem();
    gsMatrix<> solVec;
    auto u_sol = A.getSolution(u2,solVec);
    gsInfo<< "* Values (solution):\t";
    solVec = coefs;
    solVec.resize(2*solVec.rows(),1);
    result.transpose() = ev.eval( u_sol, point );
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

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
    exact = mp.patch(0).deriv(point);
    exact.resize(2,2);// = transposed Jacobian Matrix

    gsInfo<< "* Gradient  (space):\t";
    result.transpose() = ev.eval( grad(u), point );
    result *= coefs;
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    gsInfo<< "* Jacobian  (map):\t";
    result = ev.eval( jac(G), point );
    // note: transpose() below should not be needed!
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    // Note: deriv(.) returns a column with the gradients
    // transposition is required to obtain Jacobian matrix
    gsInfo<<( (result-exact.transpose()).norm() < 1e-9 ? "passed" : "failed" ) <<"\n";

    gsInfo<< "* Gradient  (solution):\t";
    result = ev.eval( grad(u_sol), point );
    // note: transpose() below should not be needed!
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    // Note: deriv(.) returns columns filled with gradients
    // transposition is required to obtain the gradients in the rows
    gsInfo<<( (result-exact.transpose()).norm() < 1e-10 ? "passed" : "failed" ) <<"\n";


    // DOES NOT WORK (sizes seem to be wrong)
    /*
      Computes the frobenius product between the gradients of a space

      Assessment of:
      - frprod_expr (version 1)
    */
    /*
      gsInfo<< "* frprod:\t";
      result = ev.eval(jac(u2).tr()%jac(u2).tr(),point);
      exact  = ev.eval(jac(u2).tr().trace(),point);
      if (verbose)
          gsInfo  <<"Result:\n"<<result<<"\n"
                  <<"Exact:\n"<<exact<<"\n";
      gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    */
    /*
      Computes the frobenius product between the gradient of a space and an id matrix

      Assessment of:
      - frprod_expr (version 2)
      - idMat_expr
    */
    gsInfo<< "* frprod (space * matrix):\t";
    result = ev.eval(jac(u2).tr()%gismo::expr::id(2),point);
    exact  = ev.eval(jac(u2).tr().trace(),point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact.transpose()<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Scalar Division
      Assessment of:
      - divide_expr (version 1)
    */
    gsInfo<< "* 1/(scalar func):\t";
    result = ev.eval(b*b/b.val(),point);
    exact  = ev.eval(b,point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result(0,0)<<"\n"
                <<"Exact:\n"<<exact(0,0)<<"\n";
    gsInfo<<( (result(0,0)-exact(0,0)) < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Space division
      Assessment of:
      - divide_expr (version 2)
    */
    gsInfo<< "* space/(scalar func):\t";
    result = ev.eval(u/b.val(),point);
    exact  = ev.eval(u,point);
    tmp = ev.eval(b,point);
    exact *= 1/tmp(0,0);
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact.transpose()<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Space division
      Assessment of:
      - divide_expr (version 1)
    */
    gsInfo<< "* space/(scalar func):\t";
    result = ev.eval(1/b.val(),point);
    exact  = ev.eval(b,point);
    if (verbose)
        gsInfo  <<"Result:\n"<<result(0,0)<<"\n"
                <<"Exact:\n"<<1/exact(0,0)<<"\n";
    gsInfo<<( (result(0,0)-1/exact(0,0)) < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Summ
      Assessment of:
      - summ_expr
    */

    // gsInfo<< "* summ:\t";
    // result = ev.eval(summ(u2.tr(),jac(u2)).tr(),point);
    // if (verbose)
    //     gsInfo  <<"Result:\n"<<result<<"\n";
    //             // <<"Exact:\n"<<exact<<"\n";
    // // gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    // gsInfo<<"passed (no exact check)\n";

    /*
      Computes the second derivatives on the geometry

      Assessment of:
      - hess_expr(gsGeometryMap)
      - hess_expr(gsFeSolution)
      - lapl_expr(gsFeSolution)
    */
    mp.patch(0).deriv2_into(point,exact);

    gsWarn<<"The following expressions have different formats:\n";


    result = ev.eval( hess(u_sol), point );
    gsInfo<< "* Hess (solution): \n"<<result<<"\n";

    index_t numDers = mp.domainDim() * (mp.domainDim() + 1) / 2;
    result.resize(mp.targetDim() * numDers,exact.cols());
    gsInfo<< "* Deriv2 (geometryMap): \n"<<exact<<"\n";
    tmp = ev.eval( hess(G), point );
    gsInfo<< "* Hess (map): \n";
    for (index_t k=0; k!=mp.targetDim(); k++)
    {
      result(numDers*k+0,0) = tmp(0,2*k);
      result(numDers*k+1,0) = tmp(1,2*k+1);
      result(numDers*k+2,0) = tmp(1,2*k);
    }
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact.transpose()<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    tmp = exact;
    exact.resize(mp.targetDim(),1);
    result = ev.eval( lapl(u_sol), point );
    for (index_t k=0; k!=mp.targetDim(); k++)
    {
      exact(k,0) = 0;
      for (index_t l=0; l!=mp.domainDim(); l++)
        exact(k,0) += tmp(k*numDers+l,0);
    }
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact.transpose()<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    /*
      Computes some expressions for assembly
      NOTE: the expressions are multiplied by the components such that the
      expression with the solution (exact) is represented.

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
    exact = ev.eval( u_sol.tr() * grad(u_sol).tr(), point );
    //space=2,cb=1 / m1
    auto e1 = (u_sol.tr() * jac(u2).tr()); //note: transposition is blockwise
    result = ev.eval( e1.tr(), point );

    gsDebugVar(jac(u2).rows());
    // DOES NOT WORK!
    // gsDebugVar(jac(u2).cols());

    // This one DOES NOT WORK!! (gives e1.rows()=e1.rows()=2)
    // auto e1 = (jac(u2) * u_sol); //note: transposition is blockwise
    // result = ev.eval( e1, point );

    result *= solVec;
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact.transpose()).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    //--------------------------------------------------------------

    gsInfo<<"* u grad(s):\t\t";
    //space=2,cb=0 / m1

    // This one WORKS
    // auto e2 = ( u2 * grad(u_sol).tr() );
    // result = ev.eval( e2.tr(), point );

    // This one WORKS
    auto e2 = ( grad(u_sol) * u2.tr() );
    result = ev.eval( e2, point );

    result *= solVec;
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact.transpose()).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    gsInfo<<"* s grad(u) + u grad(s):\t";
    //space=1,cb=1
    result = ev.eval( e1.tr() + e2, point );
    result *= solVec;
    if (verbose)
        gsInfo  <<"Result:\n"<<result.transpose()<<"\n"
                <<"Exact:\n"<<2*exact<<"\n";
    gsInfo<<( (result-2*exact.transpose()).norm() < 1e-10 ? "passed" : "failed" )<<"\n";


    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"---------------------------------Assemblers------------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    A.initSystem();
    A.assemble(u2 * u2.tr(),u2 * u_sol);
    gsInfo<<( (A.matrix()*solVec-A.rhs()).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    real_t integral = ev.integral(u_sol.tr() * grad(u_sol) * grad(u_sol).tr() * u_sol);

    A.initSystem();
    auto lhs = u2 * grad(u_sol) * grad(u_sol).tr() * u2.tr();
    auto rhs = u2 * grad(u_sol) * grad(u_sol).tr() * u_sol;
    A.assemble(lhs,rhs);

    result = solVec.transpose() * A.matrix() * solVec;
    GISMO_ENSURE(result.rows()==1,"Result must be scalar.");
    GISMO_ENSURE(result.cols()==1,"Result must be scalar.");
    gsInfo<<( std::abs(result(0,0) - integral) < 1e-10 ? "passed" : "failed" )<<"\n";
    gsInfo<<( std::abs(result(0,0) - integral) < 1e-10 ? "passed" : "failed, value = " + std::to_string(std::abs(result(0,0) - integral)) )<<"\n";
    GISMO_ENSURE(result.rows()==1,"Result must be scalar.");
    GISMO_ENSURE(result.cols()==1,"Result must be scalar.");
    result = A.rhs().transpose() * solVec;
    gsInfo<<( std::abs(result(0,0) - integral) < 1e-10 ? "passed" : "failed, value = " + std::to_string(std::abs(result(0,0) - integral)) )<<"\n";

    A.initSystem();
    auto lhs2 = u_sol.tr() * jac(u2) * jac(u2).tr() * u_sol;
    auto rhs2 = u2 * grad(u_sol) * grad(u_sol).tr() * u_sol;

    // gsDebug<<ev.eval(lhs2,point)<<"\n";
    gsDebugVar(lhs2.rows());
    // THE ONE BELOW GIVES AN ERROR!
    gsDebugVar(lhs2.cols());
    // gsDebug<<ev.eval(rhs2,point)<<"\n";

    // A.assemble(lhs2,rhs2);

    // result = solVec.transpose() * A.matrix() * solVec;
    // GISMO_ENSURE(result.rows()==1,"Result must be scalar.");
    // GISMO_ENSURE(result.cols()==1,"Result must be scalar.");
    // gsInfo<<( std::abs(result(0,0) - integral) < 1e-10 ? "passed" : "failed" )<<"\n";
    // GISMO_ENSURE(result.rows()==1,"Result must be scalar.");
    // GISMO_ENSURE(result.cols()==1,"Result must be scalar.");
    // result = A.rhs().transpose() * solVec;
    // gsInfo<<( std::abs(result(0,0) - integral) < 1e-10 ? "passed" : "failed" )<<"\n";


    return EXIT_SUCCESS;
}
