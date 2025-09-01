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
    gsInfo << expr<<" (order=" << expr.order
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
    gsInfo<<", deriv=" << expr.getDerivative()
       << ", isConstant=" << (expr.isConstant ? "true" : "false")
       << ")\n";
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

    // bool verbose = false;
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
    auto DDδ = Expr::div(Dδ);
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

    gsInfo << "\n=== Vector Calculus Identity Tests ===\n";

    // Create additional expressions for testing
    gsConstantFunction<> cfun1(5.0,3);
    gsConstantFunction<> cfun2(7.0,3);
    gsConstantFunction<> cfun3(gsVector<>::vec(1.0,2.0,3.0),3);
    gsConstantFunction<> cfun4(gsVector<>::vec(4.0,5.0,6.0),3);
    auto ψ = helper.getScalarFunction(cfun1);  // Another scalar function
    auto φ = helper.getScalarFunction(cfun2);  // Another scalar function
    auto A = helper.getVectorFunction(cfun3);       // Vector field A
    auto B = helper.getVectorFunction(cfun4);       // Vector field B (same as A for testing)

    gsInfo << "\n--- Operator tests ---\n";
    print(ψ);
    print(φ);
    print(A);
    print(B);

    print(grad(ψ));
    print(div(A));
    print(curl(A));
    print(lapl(ψ));

    print(ψ + φ);
    print(ψ - φ);
    print(ψ * φ);
    print(ψ / φ);

    gsInfo << "\n--- Gradient Identities ---\n";
    print(grad(ψ + φ));                      // ∇(ψ + φ) = ∇ψ + ∇φ
    print(grad(ψ * φ));                      // ∇(ψφ) = ψ∇φ + φ∇ψ
    print(grad(ψ / φ));                      // ∇(ψ/φ) = (φ∇ψ - ψ∇φ)/φ²
    print(grad(dot(A, B)));                  // ∇(A·B) = (A·∇)B + (B·∇)A + A×(∇×B) + B×(∇×A)
    print(grad(cross(A, B)));                // ∇(A×B) = (∇A)×B - (∇B)×A

    gsInfo << "\n--- Divergence Identities ---\n";
    print(div(A + B));                       // ∇·(A + B) = ∇·A + ∇·B
    print(div(φ * A));                       // ∇·(φA) = φ∇·A + A·∇φ
    print(div(cross(A, B)));                 // ∇·(A×B) = B·(∇×A) - A·(∇×B)
    print(div(grad(ψ)));                     // ∇·(∇ψ) = ∇²ψ (Laplacian)

    gsInfo << "\n--- Curl Identities ---\n";
    print(curl(A + B));                      // ∇×(A + B) = ∇×A + ∇×B
    print(curl(φ * A));                      // ∇×(φA) = φ(∇×A) + (∇φ)×A
    print(curl(cross(A, B)));                // ∇×(A×B) = A(∇·B) - B(∇·A) + (B·∇)A - (A·∇)B

    gsInfo << "\n--- Laplacian Identities ---\n";
    print(lapl(ψ + φ));                      // ∇²(ψ + φ) = ∇²ψ + ∇²φ
    print(lapl(ψ * φ));                      // ∇²(ψφ) = ψ∇²φ + 2∇ψ·∇φ + φ∇²ψ
    print(lapl(ψ / φ));                      // ∇²(ψ/φ) = (φ∇²ψ - ψ∇²φ - 2∇ψ·∇φ)/φ²

    gsInfo << "\n--- Higher-Order Identities ---\n";
    print(grad(grad(ψ)));                    // ∇(∇ψ) = Hessian matrix
    print(curl(grad(ψ)));                    // ∇×(∇ψ) = 0 (curl of gradient is zero)
    print(div(curl(A)));                     // ∇·(∇×A) = 0 (divergence of curl is zero)
    print(curl(curl(A)));                    // ∇×(∇×A) = ∇(∇·A) - ∇²A
    // print(grad(div(A)));                     // ∇(∇·A)
    // print(lapl(grad(ψ)));                    // ∇²(∇ψ) = ∇(∇²ψ) = gradient of Laplacian

    gsInfo << "\n--- Product Rule Combinations ---\n";
    print(grad(ψ * A));                      // ∇(ψA) = ψ∇A + A⊗∇ψ (tensor product)
    print(div(ψ * A));                       // ∇·(ψA) = ψ∇·A + A·∇ψ
    print(curl(ψ * A));                      // ∇×(ψA) = ψ(∇×A) + (∇ψ)×A
    print(grad(A / φ));                      // ∇(A/φ) = (φ∇A - A⊗∇φ)/φ²

    // gsInfo << "\n--- Triple Product Identities ---\n";
    // print(div(cross(A, cross(B, A))));       // ∇·(A×(B×A))
    // print(curl(cross(A, cross(B, A))));      // ∇×(A×(B×A))
    // print(grad(dot(A, cross(B, A))));        // ∇(A·(B×A))

    return EXIT_SUCCESS;
}
