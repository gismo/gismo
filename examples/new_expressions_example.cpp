/** @file gsExpressions_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsNewExpressions/NewExpressions.h>
#include <gsNewExpressions/ExpressionHelper.h>
using namespace gismo;
using namespace Expr;

template <typename E>
void print(const BaseExpression<E> & expr)
{
    gsInfo << "Expression(order=" << expr.order
       << ", space=" << expr.space
       << ", sizes=";
    if (expr.order>0)
    {
        gsInfo<<"[";
        for (size_t i = 0; i < expr.order; ++i)
        {
            gsInfo << expr.sizes()[i];
            if (i < expr.order - 1)
                gsInfo << ", ";
        }
        gsInfo << "]";
    }
    else
        gsInfo<<"none";
    gsInfo<<")\n";
}


template <typename E>
gsMatrix<typename E::Scalar> eval(ExpressionHelper<typename E::Scalar> & helper,
                                  const E & expr,
                                  const gsMatrix<typename E::Scalar> & pts)
{
    helper.points() = pts;
    expr.parse(helper);
    gsInfo<<expr<<", max derivative = "<<expr.getDerivative()<<"\n";
    helper.precompute();
    return expr.eval(0);
}


int main(int argc, char *argv[])
{

    bool verbose = false;
    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsVector<real_t> βvec = gsVector<real_t>::LinSpaced(3,1.0,3.0);
    gsMatrix<real_t> γmat(3,3);
    γmat << 1.0,2.0,3.0,
            4.0,5.0,6.0,
            7.0,8.0,9.0;

    gsMatrix<real_t> points(2,3);
    points << 0.0, 0.5, 1.0,
              0.0, 0.5, 1.0;

    gsFunctionExpr<> fun("sin(x) + cos(y)",2);

    ExpressionHelper<real_t> helper;
    helper.points() = points;

    auto α = helper.getConstant(2.0);
    auto β = helper.getConstant(βvec);
    auto γ = helper.getConstant(γmat);
    auto δ = helper.getScalarFunction(fun);

    print(α);
    print(β);
    print(γ);
    print(δ);


    gsDebug<<"α + α: "<< eval(helper, α + α, points)<<"\n";
    // gsDebug<<"α + β: "<< eval(helper, α + β, points)<<"\n";
    // gsDebug<<"β + α: "<< eval(helper, β + α, points)<<"\n";
    gsDebug<<"β + β: "<< eval(helper, β + β, points)<<"\n";
    gsDebug<<"γ + γ: "<< eval(helper, γ + γ, points)<<"\n";
    gsDebug<<"δ + δ: "<< eval(helper, δ + δ, points)<<"\n";

    auto β_T = Expr::transpose(β);
    print(β_T);
    // gsDebug<<"β_T: "<< eval(helper, β_T, points)<<"\n";
    // gsDebug<<"β + β_T: "<< eval(helper, β + β_T, points)<<"\n";

    auto Dδ = Expr::grad(δ);
    print(Dδ);
    auto DDδ = Expr::grad(Dδ);
    eval(helper,DDδ,points);
    // gsDebugVar(DDδ.max


    gsDebug<<"α * α: "<< eval(helper, α * α, points)<<"\n";
    gsDebug<<"α * β: "<< eval(helper, α * β, points)<<"\n";
    // gsDebug<<"β * α: "<< eval(helper, β * α, points)<<"\n";
    // gsDebug<<"β * β: "<< eval(helper, β * β, points)<<"\n";
    gsDebug<<"γ * γ: "<< eval(helper, γ * γ, points)<<"\n";
    gsDebug<<"δ * δ: "<< eval(helper, δ * δ, points)<<"\n";

    gsDebug<<"dot(β,β): "<< eval(helper, dot(β,β), points)<<"\n";
    gsDebug<<"cross(β,β): "<< eval(helper, cross(β,β), points)<<"\n";


    // gsDebug<<"Dδ, order = "<<Dδ.order<<", size = "<<Dδ.sizes()[0]<<"\n";

    // ScalarExpression<real_t> scalar_expr(5.0);
    // gsVector<real_t> vector_data = gsVector<real_t>::LinSpaced(3, 1.0, 3.0);
    // VectorExpression<real_t> vector_expr(vector_data);
    // gsMatrix<real_t> matrix_data(2, 2);
    // matrix_data << 1.0, 2.0,
    //                3.0, 4.0;
    // MatrixExpression<real_t> matrix_expr(matrix_data);

    // // Test addition of scalar + scalar
    // auto scalar_plus_scalar = scalar_expr + scalar_expr;
    // gsMatrix<real_t> result_scalar = scalar_plus_scalar.eval(0);
    // std::cout << "Result of scalar + scalar:\n" << result_scalar << std::endl;

    // // // Test addition of scalar + vector
    // // auto scalar_plus_vector = scalar_expr + vector_expr;
    // // gsMatrix<real_t> result_vector = scalar_plus_vector.eval(0);
    // // std::cout << "Result of scalar + vector:\n" << result_vector << std::endl;

    // // // Test addition of vector + scalar
    // // auto vector_plus_scalar = vector_expr + scalar_expr;
    // // gsMatrix<real_t> result_vector2 = vector_plus_scalar.eval(0);
    // // std::cout << "Result of vector + scalar:\n" << result_vector2 << std::endl;

    // // Test addition of two vectors
    // auto vector_plus_vector = vector_expr + vector_expr;
    // gsMatrix<real_t> result_vector3 = vector_plus_vector.eval(0);
    // std::cout << "Result of vector + vector:\n" << result_vector3 << std::endl;

    // // // Test addition of scalar + matrix
    // // auto scalar_plus_matrix = scalar_expr + matrix_expr;
    // // gsMatrix<real_t> result_matrix = scalar_plus_matrix.eval(0);
    // // std::cout << "Result of scalar + matrix:\n" << result_matrix << std::endl;

    // // // Test addition of matrix + scalar
    // // auto matrix_plus_scalar = matrix_expr + scalar_expr;
    // // gsMatrix<real_t> result_matrix2 = matrix_plus_scalar.eval(0);
    // // std::cout << "Result of matrix + scalar:\n" << result_matrix2 << std::endl;

    // gsFunctionExpr<> fun("sin(x) + cos(y)",2);
    // ScalarVariableExpression<real_t> scalar_var_expr(fun); // order 0

    // // // Differetial operators of constant scalars, should all evaluate to null-expressions of the correct order
    // // GradExpression<ScalarExpression<real_t>> grad_expr(scalar_expr); // order 1
    // // DivExpression<GradExpression<ScalarExpression<real_t>>> div_grad_expr(grad_expr); // order 0
    // // GradExpression<GradExpression<ScalarExpression<real_t>>> grad_grad_expr(grad_expr); // order 2

    // // // Differetial operators of variable scalars, should evaluate to differential operators of the correct order
    // // GradExpression<ScalarVariableExpression<real_t>> grad_expr(scalar_var_expr); // order 1
    // // DivExpression<GradExpression<ScalarVariableExpression<real_t>>> div_grad_expr(grad_expr); // order 0
    // // GradExpression<GradExpression<ScalarVariableExpression<real_t>>> grad_grad_expr(grad_expr); // order 2

    // // First-order derivative
    // DerivativeTensor<1,FunctionExpression<0,real_t>> // second order
    // DerivativeTensor<1,FunctionExpression<1,real_t>> // third order

    // // Second-order derivative
    // DerivativeTensor<2,FunctionExpression<0,real_t>> // second order
    // DerivativeTensor<2,FunctionExpression<1,real_t>> // third order

    // //
    // grad(grad(s))
    // div(grad(s)) // alias deriv2(s).sum()
    // curl(f1)

    // // NEW NAMING

    // ExprEvaluator<> ev;
    // ConstantExpression<0,real_t> s = ev.template Constant<0>(0.5);           // Template is not needed, but added for consistency with function (order can be derived from the type)
    // ConstantExpression<1,real_t> v = ev.template Constant<1>(vector_data);   // Template is not needed, but added for consistency with function (order can be derived from the type)
    // ConstantExpression<2,real_t> m = ev.template Constant<2>(matrix_data);   // Template is not needed, but added for consistency with function (order can be derived from the type)
    // FunctionExpression<0,real_t> f = ev.template Function<0>(fun);           // Template is needed because we can only deduce the order from the target dimension of the function, which is runtime
    // FunctionExpression<0,real_t> f0= ev.template Function0(fun.component(0));           // Template is needed because we can only deduce the order from the target dimension of the function, which is runtime
    // FunctionExpression<1,real_t> f1= ev.template Function1(fun);           // Template is needed because we can only deduce the order from the target dimension of the function, which is runtime
    // FunctionExpression<2,real_t> f2= ev.template Function2(fun,n,m); //reshapes          // Template is needed because we can only deduce the order from the target dimension of the function, which is runtime
    // auto hess = grad(grad(s)) = hess(s);
    // auto grad_hess = grad(hess);
    // grad(s*grad(s)) = grad(s) * grad(s) + s * grad(grad(s)) = grad(s) * grad(s) + s * hess(s);
    // ArrayExpression<ConstantExpression<1,real_t>> v_array = v.array(); //
    // auto add_scalar_vector = s+v_array;



    // SpaceExpression<0,Space::Test,real_t> phi= A.ScalarTrialSpace();
    // SpaceExpression<0,Space::Trial,real_t> psi= A.ScalarTestSpace();



    // MultExpression<

    // s * v // always fine in Eigen
    // s / v // not allowed in Eigen, use s / v.array() instead


    // // COMPLETELY AVOID TRANSPOSED VECTORS IN EXPRESSIONS
    // v * w // not allowed because both are column vectors, use v.array() * v.array() instead
    // v.transpose() * w = inner(v,w) = dot(v,w) = v % w    // order 0, also for matrices also for scalar
    // v * w.transpose() = outer(v,w) // order 2, only for vectors
    //                     cross(v,w) // order 1, only for vectors

    //                     times()

    // v.transpose() * m * w = inner(v,m*w)
    // v.transpose() * m.transpose() * w = inner(v,m.transpose()*w)

    // grad(s) // allow
    // inner(Expr::NABLA,s) = grad(s) // shortcut

    // // Transpose expression
    // s.transpose() //  not allows
    // v.transpose() // not allowed
    // m.transpose() // allowed

    // m.flatten() // Eigen::asVector()
    // m.squeeze() //

    // // OPERATORS
    // // * AdditionOperator       (+) -- for everything of same order, and scalar with any order
    // // * SubtractionOperator    (-) -- for everything of same order, and scalar with any order

    // // * ScalarMultiplicationOperator (*) -- scalar * everything, everything * scalar
    // // * MatrixMultiplicationOperator (*) -- mat*mat, mat*vec
    // // * TensorMultiplicationOperator (times) -- allows to multiply tensors of different orders, e.g., vector * matrix, specifying the axis of contraction
    // // * HadamardProductOperator (*) -- vec*vec, mat*mat (automatically activated via .array() method)

    // // * DivisionOperator       (/) -- by scalar, or by vector/matrix with array() method

    // // * InnerProductOperator
    // // * OuterProductOperator
    // // * CrossProductOperator


    // // DIFFERENTIAL OPERATORS
    // // * GradientOperator       (grad) -- grad(s), grad(v), grad(m)
    // // * DivergenceOperator     (div) -- div(v), div(m)
    // // * CurlOperator           (curl) -- curl(v), curl(m)
    // // * HessianOperator        (hess) -- hess(s), hess(v), hess(m)
    // // * LaplacianOperator      (laplacian) -- laplacian(s), laplacian(v), laplacian(m)

    // times(grad(m),0,v,0) // multiply grad(m) with v, contracting the first axis of grad(m) and the first axis of v
    // times(grad(m).slice(0),v,0) // multiply grad(m) with v, contracting the first axis of grad(m) and the first axis of v


    return EXIT_SUCCESS;
}
