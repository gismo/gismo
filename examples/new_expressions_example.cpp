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
    gsInfo << expr<<" (order=" << expr.order()
       << ", space=" << expr.space()
       << ", sizes=";
    if (expr.order()>0)
    {
        gsInfo<<"[";
        for (size_t i = 0; i < expr.order(); ++i)
        {
            gsInfo << expr.sizes()[i];
            if (i < expr.order() - 1)
                gsInfo << ", ";
        }
        gsInfo << "]";
    }
    else
        gsInfo<<"none";
    gsInfo<<", deriv=" << expr.getDerivative()
       << ", isConstant=" << (expr.isConstant() ? "true" : "false")
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
    return expr.eval(0)(0,0);
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

    auto α = helper.getConstant(2.0,"α");
    auto β = helper.getConstant(βvec,"β");
    auto γ = helper.getConstant(γmat,"γ");
    auto δ = helper.getScalarFunction(fun,"δ");

    print(α);
    print(β);
    print(γ);
    print(δ);


    gsDebug<<"α + α: "<< eval(helper, α + α, points)<<"\n";
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
    gsInfo << "Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities\n";

    // Create additional expressions for testing
    gsConstantFunction<> cfun1(5.0,3);
    gsConstantFunction<> cfun2(7.0,3);
    gsConstantFunction<> cfun3(gsVector<>::vec(1.0,2.0,3.0),3);
    gsConstantFunction<> cfun4(gsVector<>::vec(4.0,5.0,6.0),3);
    auto ψ = helper.getScalarFunction(cfun1,"ψ");  // Scalar function ψ
    auto φ = helper.getScalarFunction(cfun2,"φ");  // Scalar function φ
    auto A = helper.getVectorFunction(cfun3,"A");  // Vector field A
    auto B = helper.getVectorFunction(cfun4,"B");  // Vector field B

    gsInfo << "\n--- Operator definitions ---\n";
    gsInfo<<"ψ = "; print(ψ);
    gsInfo<<"φ = "; print(φ);
    gsInfo<<"A = "; print(A);
    gsInfo<<"B = "; print(B);

    // ==========================================================================
    // FIRST DERIVATIVE IDENTITIES
    // Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities#First_derivative_identities
    // ==========================================================================

    gsInfo << "\n\n=== FIRST DERIVATIVE IDENTITIES ===\n";
    gsInfo << "Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities#First_derivative_identities\n";

    // Distributive properties
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Distributive_properties
    gsInfo << "\n--- Distributive Properties ---\n";
    gsInfo << "∇(ψ + φ) = ∇ψ + ∇φ\n";
    gsInfo<<"∇(ψ + φ) = "; print(grad(ψ + φ));

    gsInfo << "∇(A + B) = ∇A + ∇B\n";
    gsInfo<<"∇(A + B) = "; print(grad(A + B));

    gsInfo << "∇·(A + B) = ∇·A + ∇·B\n";
    gsInfo<<"∇·(A + B) = "; print(div(A + B));

    gsInfo << "∇×(A + B) = ∇×A + ∇×B\n";
    gsInfo<<"∇×(A + B) = "; print(curl(A + B));

    // Product rule for multiplication by a scalar
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Product_rule_for_multiplication_by_a_scalar
    gsInfo << "\n--- Product Rule for Multiplication by a Scalar ---\n";
    gsInfo << "∇(ψφ) = φ∇ψ + ψ∇φ\n";
    gsInfo<<"∇(ψφ) = "; print(grad(ψ * φ));

    gsInfo << "∇·(ψA) = ψ∇·A + (∇ψ)·A\n";
    gsInfo<<"∇·(ψA) = "; print(div(ψ * A));

    gsInfo << "∇×(ψA) = ψ∇×A + (∇ψ)×A\n";
    gsInfo<<"∇×(ψA) = "; print(curl(ψ * A));

    gsInfo << "∇²(ψφ) = ψ∇²φ + 2∇ψ·∇φ + φ∇²ψ\n";
    gsInfo<<"∇²(ψφ) = "; print(lapl(ψ * φ));

    // Quotient rule for division by a scalar
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Quotient_rule_for_division_by_a_scalar
    gsInfo << "\n--- Quotient Rule for Division by a Scalar ---\n";
    gsInfo << "∇(ψ/φ) = (φ∇ψ - ψ∇φ)/φ²\n";
    gsInfo<<"∇(ψ/φ) = "; print(grad(ψ / φ));

    gsInfo << "∇·(A/φ) = (φ∇·A - ∇φ·A)/φ²\n";
    gsInfo<<"∇·(A/φ) = "; print(div(A / φ));

    gsInfo << "∇×(A/φ) = (φ∇×A - ∇φ×A)/φ²\n";
    gsInfo<<"∇×(A/φ) = "; print(curl(A / φ));

    // Chain rule - would need composed functions, commenting out for now
    gsInfo << "\n--- Chain Rule ---\n";
    gsInfo << "// Chain rule identities require function composition - not easily testable here\n";
    gsInfo << "// ∇(f∘φ) = (f'∘φ)∇φ\n";
    gsInfo << "// ∇(φ∘A) = (∇A)·(∇φ∘A)\n";

    // Dot product rule
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Dot_product_rule
    gsInfo << "\n--- Dot Product Rule ---\n";
    gsInfo << "∇(A·B) = (A·∇)B + (B·∇)A + A×(∇×B) + B×(∇×A)\n";
    gsInfo<<"∇(A·B) = "; print(grad(dot(A, B)));

    // Cross product rule
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Cross_product_rule
    gsInfo << "\n--- Cross Product Rule ---\n";
    gsInfo << "∇·(A×B) = (∇×A)·B - A·(∇×B)\n";
    gsInfo<<"∇·(A×B) = "; print(div(cross(A, B)));

    gsInfo << "∇×(A×B) = A(∇·B) - B(∇·A) + (B·∇)A - (A·∇)B\n";
    gsInfo<<"∇×(A×B) = "; print(curl(cross(A, B)));

    // ==========================================================================
    // SECOND DERIVATIVE IDENTITIES
    // Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities#Second_derivative_identities
    // ==========================================================================

    gsInfo << "\n\n=== SECOND DERIVATIVE IDENTITIES ===\n";
    gsInfo << "Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities#Second_derivative_identities\n";

    // Divergence of curl is zero
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Divergence_of_curl_is_zero
    gsInfo << "\n--- Divergence of Curl is Zero ---\n";
    gsInfo << "∇·(∇×A) = 0\n";
    gsInfo<<"∇·(∇×A) = "; print(div(curl(A)));

    // Divergence of gradient is Laplacian
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Divergence_of_gradient_is_Laplacian
    gsInfo << "\n--- Divergence of Gradient is Laplacian ---\n";
    gsInfo << "Δψ = ∇²ψ = ∇·(∇ψ)\n";
    gsInfo<<"∇·(∇ψ) = "; print(div(grad(ψ)));

    // Curl of gradient is zero
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Curl_of_gradient_is_zero
    gsInfo << "\n--- Curl of Gradient is Zero ---\n";
    gsInfo << "∇×(∇φ) = 0\n";
    gsInfo<<"∇×(∇φ) = "; print(curl(grad(φ)));

    // Curl of curl
    // https://en.wikipedia.org/wiki/Vector_calculus_identities#Curl_of_curl
    gsInfo << "\n--- Curl of Curl ---\n";
    gsInfo << "∇×(∇×A) = ∇(∇·A) - ∇²A\n";
    gsInfo<<"WRONG: ∇×(∇×A) = "; print(curl(curl(A)));

    // ==========================================================================
    // SUMMARY OF IMPORTANT IDENTITIES
    // Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities#Summary_of_important_identities
    // ==========================================================================

    gsInfo << "\n\n=== SUMMARY OF IMPORTANT IDENTITIES ===\n";
    gsInfo << "Reference: https://en.wikipedia.org/wiki/Vector_calculus_identities#Summary_of_important_identities\n";

    // Gradient identities
    gsInfo << "\n--- Gradient ---\n";
    gsInfo << "• ∇(ψ + φ) = ∇ψ + ∇φ\n";
    gsInfo<<"  Result: "; print(grad(ψ + φ));

    gsInfo << "• ∇(ψφ) = φ∇ψ + ψ∇φ\n";
    gsInfo<<"  Result: "; print(grad(ψ * φ));

    gsInfo << "• ∇(ψA) = ∇ψ⊗A + ψ∇A\n";
    gsInfo<<"  Result: "; print(grad(ψ * A));

    gsInfo << "• ∇(A·B) = (A·∇)B + (B·∇)A + A×(∇×B) + B×(∇×A)\n";
    gsInfo<<"  Result: "; print(grad(dot(A, B)));

    // Divergence identities
    gsInfo << "\n--- Divergence ---\n";
    gsInfo << "• ∇·(A + B) = ∇·A + ∇·B\n";
    gsInfo<<"  Result: "; print(div(A + B));

    gsInfo << "• ∇·(ψA) = ψ∇·A + A·∇ψ\n";
    gsInfo<<"  Result: "; print(div(ψ * A));

    gsInfo << "• ∇·(A×B) = (∇×A)·B - (∇×B)·A\n";
    gsInfo<<"  Result: "; print(div(cross(A, B)));

    // Curl identities
    gsInfo << "\n--- Curl ---\n";
    gsInfo << "• ∇×(A + B) = ∇×A + ∇×B\n";
    gsInfo<<"  Result: "; print(curl(A + B));

    gsInfo << "• ∇×(ψA) = ψ(∇×A) + (∇ψ)×A\n";
    gsInfo<<"  Result: "; print(curl(ψ * A));

    gsInfo << "• ∇×(ψ∇φ) = ∇ψ×∇φ\n";
    gsInfo<<"  Result: "; print(curl(ψ * grad(φ)));

    gsInfo << "• ∇×(A×B) = A(∇·B) - B(∇·A) + (B·∇)A - (A·∇)B\n";
    gsInfo<<"  Result: "; print(curl(cross(A, B)));

    // Second derivatives
    gsInfo << "\n--- Second Derivatives ---\n";
    gsInfo << "• ∇·(∇×A) = 0\n";
    gsInfo<<"  Result: "; print(div(curl(A)));

    gsInfo << "• ∇×(∇ψ) = 0\n";
    gsInfo<<"  Result: "; print(curl(grad(ψ)));

    gsInfo << "• ∇·(∇ψ) = ∇²ψ (scalar Laplacian)\n";
    gsInfo<<"  Result: "; print(div(grad(ψ)));

    gsInfo << "• ∇×(∇×A) = ∇(∇·A) - ∇²A (vector Laplacian identity)\n";
    gsInfo<<"  WRONG [needs implementation of grad(div), which is different from div(grad)]: "; print(curl(curl(A)));

    gsInfo << "\n--- Additional Tests ---\n";
    gsInfo<<"∇(ψA) = "; print(grad(ψ * A));                      // ∇(ψA) = ψ∇A + A⊗∇ψ (tensor product)
    gsInfo<<"∇(A/φ) = "; print(grad(A / φ));                      // ∇(A/φ) = (φ∇A - A⊗∇φ)/φ²

    gsInfo << "\n--- Complex Expressions ---\n";
    gsInfo<<"∇·(A×(B×A)) = "; print(div(cross(A, cross(B, A))));       // Triple cross product divergence
    gsInfo<<"∇×(A×(B×A)) = "; print(curl(cross(A, cross(B, A))));      // Triple cross product curl
    gsInfo<<"∇(A·(B×A)) = ";  print(grad(dot(A, cross(B, A))));        // Scalar triple product gradient

    gsInfo << "\n--- TESTING DIRECTIONAL DERIVATIVE (A·∇)B ---\n";

    // Test if we can compute (A·∇)B as dot(A, grad(B))
    gsInfo << "(A·∇)B = ";
    try {
        auto directional_derivative = dot(A, grad(B));
        print(directional_derivative);
    } catch (const std::exception& e) {
        gsInfo << "ERROR: " << e.what() << "\n";
    }

    gsInfo << "\n--- DEBUGGING TEMPLATE MATCHING ---\n";

    // Test what type A and B create when combined with dot
    auto testDot = dot(A, B);
    gsInfo << "Type of dot(A,B): ";
    print(testDot);

    // Test what type grad creates when applied to the dot product
    auto testGradDot = grad(testDot);
    gsInfo << "Type of grad(dot(A,B)): ";
    print(testGradDot);

    // Test curl of curl
    auto testCurl = curl(A);
    gsInfo << "Type of curl(A): ";
    print(testCurl);

    auto testCurlCurl = curl(testCurl);
    gsInfo << "Type of curl(curl(A)): ";
    print(testCurlCurl);


    auto space = helper.getScalarTrialSpace(cfun1,0,"ψ");
    auto space2 = helper.getVectorTrialSpace(cfun3,1,2,"ψ");
    gsDebugVar(space2.space());
    gsDebugVar(space2.order());
    /*
    gsMatrix<> solvector;
    auto sol = helper.getSolution(space,solvector);
    auto sol2 = helper.getSolution(space2,solvector);
    auto dsol= variation(sol,space);
    // auto variation = variation(grad(sol)*space+div(f)*sol,space);

    print(sol);
    print(variation(sol,space));

    print(variation(grad(sol),space));
    print(variation(div(sol2),space));
    print(variation(curl(sol2),space));
    print(variation(lapl(sol),space));


    // print(variation(α+α,space));
    print(variation(sol+α,space));
    // print(variation(sol+sol,space));

    print(variation(sol*α,space));
    print(variation(α*sol,space));
*/
    return EXIT_SUCCESS;
}
