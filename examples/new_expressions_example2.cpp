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

template <typename E>
void printEval(ExpressionHelper<typename E::Scalar> & helper,
                                       const E & expr,
                                       const gsMatrix<typename E::Scalar> & pts)
{
    helper.points() = pts;
    expr.parse(helper);
    gsInfo<<expr<<", max derivative = "<<expr.getDerivative()<<"\n";
    helper.precompute();
    for (index_t k=0; k<pts.cols(); ++k)
    {
        gsInfo<<"At point "<<pts.col(k).transpose()<<":\n";
        ExpressionValue<typename E::Scalar> ev = expr.eval(k);
        for (index_t i=0; i!=ev.rowCardinality(); ++i)
        {
            for (index_t j=0; j!=ev.colCardinality(); ++j)
                gsInfo<<ev(i,j)<<" ";
            gsInfo<<"\n";
        }



    }
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

    gsKnotVector<> kv(0.,1.,3,2);
    gsTensorBSplineBasis<2> tbasis(kv,kv);

    gsInfo<<tbasis.eval(points)<<"\n";

    ExpressionHelper<real_t> helper;
    helper.points() = points;

    auto α = helper.getConstant(2.0,"α");
    auto β = helper.getConstant(βvec,"β");
    auto γ = helper.getConstant(γmat,"γ");
    auto δ = helper.getScalarFunction(fun,"δ");
    auto φ = helper.getScalarTestSpace(tbasis,0,"φ");
    auto η = helper.getScalarTestSpace(tbasis,2,"η");
    auto ψ = helper.getScalarTrialSpace(tbasis,1,"ψ");
    auto Ψ = helper.getVectorTestSpace(tbasis,2,3,"Ψ");
    auto Φ = helper.getVectorTrialSpace(tbasis,2,4,"Φ");

    gsInfo<<"\n--- Expressions ---\n";
    print(α);
    print(β);
    print(γ);
    print(δ);
    print(φ);
    print(ψ);
    print(Φ);
    print(Ψ);
    print(φ+φ);
    print(ψ+ψ);
    print(Φ+Φ);
    print(Ψ+Ψ);
    print(grad(φ));
    print(grad(ψ));
    print(grad(Ψ));
    print(grad(Φ));

    gsInfo<<"\n--- Evaluations ---\n";
    printEval(helper, φ, points);
    printEval(helper, ψ, points);
    printEval(helper, Φ, points);
    printEval(helper, Ψ, points);

    printEval(helper, φ + φ, points);
    printEval(helper, ψ + ψ, points);
    printEval(helper, η + η, points);
    printEval(helper, φ + η, points);
    printEval(helper, Φ + Φ, points);
    printEval(helper, Ψ + Ψ, points);

    
    printEval(helper, φ + ψ, points);
    printEval(helper, Φ + Ψ, points);
    
    printEval(helper, α * φ, points);
    printEval(helper, α * ψ, points);
    printEval(helper, α * Φ, points);
    printEval(helper, α * Ψ, points);

    printEval(helper, grad(φ), points);
    printEval(helper, grad(ψ), points);
    printEval(helper, grad(Ψ), points);
    printEval(helper, grad(Φ), points);

    // printEval(helper, grad(φ)*grad(φ), points

    gsInfo<<"\n--- Done ---\n";

    gsInfo<<"\n--- Old framework ---\n";
    gsExprAssembler<real_t> assembler(1,1);
    gsExprEvaluator<real_t> oldEval(assembler);
    auto φold = assembler.getSpace(tbasis);
    auto Φold = assembler.getSpace(tbasis,2);
    gsInfo<<"φold: at point 0: "<<points.col(0).transpose()<<"\n";
    gsInfo<<oldEval.eval(φold, points.col(0))<<"\n";
    gsInfo<<"point 1: "<<points.col(1).transpose()<<"\n";
    gsInfo<<oldEval.eval(φold, points.col(1))<<"\n";
    gsInfo<<"point 2: "<<points.col(2).transpose()<<"\n";
    gsInfo<<oldEval.eval(φold, points.col(2))<<"\n";
    gsInfo<<"grad(φold) at point 0: "<<points.col(0).transpose()<<"\n";
    gsInfo<<oldEval.eval(expr::grad(φold)*expr::grad(φold).tr(), points.col(0))<<"\n";
    gsInfo<<"grad(φold) at point 1: "<<points.col(1).transpose()<<"\n";
    gsInfo<<oldEval.eval(expr::grad(φold)*expr::grad(φold).tr(), points.col(1))<<"\n";
    gsInfo<<"grad(φold) at point 2: "<<points.col(2).transpose()<<"\n";
    gsInfo<<oldEval.eval(expr::grad(φold)*expr::grad(φold).tr(), points.col(2))<<"\n";

    gsInfo<<"Φold: at point 0: "<<points.col(0).transpose()<<"\n";
    gsInfo<<oldEval.eval(Φold, points.col(0))<<"\n";
    gsInfo<<"point 1: "<<points.col(1).transpose()<<"\n";
    gsInfo<<oldEval.eval(Φold, points.col(1))<<"\n";
    gsInfo<<"point 2: "<<points.col(2).transpose()<<"\n";
    gsInfo<<oldEval.eval(Φold, points.col(2))<<"\n";
    gsInfo<<"jac(Φold) at point 0: "<<points.col(0).transpose()<<"\n";
    gsInfo<<oldEval.eval(expr::jac(Φold), points.col(0))<<"\n";
    gsInfo<<"jac(Φold) at point 1: "<<points.col(1).transpose()<<"\n";
    gsInfo<<oldEval.eval(expr::jac(Φold), points.col(1))<<"\n";
    gsInfo<<"jac(Φold) at point 2: "<<points.col(2).transpose()<<"\n";
    gsInfo<<oldEval.eval(expr::jac(Φold), points.col(2))<<"\n";

    return EXIT_SUCCESS;
}
