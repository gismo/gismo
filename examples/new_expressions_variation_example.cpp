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
        ExpressionResult<typename E::Scalar> ev = expr.eval(k);
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

    ExpressionHelper<real_t> helper;
    helper.points() = points;

    auto α = helper.getConstant(2.0,"α");
    auto β = helper.getConstant(βvec,"β");
    auto γ = helper.getConstant(γmat,"γ");
    auto δ = helper.getScalarFunction(fun,"δ");

    gsKnotVector<> kv(0.,1.,3,2);
    gsTensorBSplineBasis<2> tbasis(kv,kv);
    gsMatrix<> coefs(tbasis.size(),1);
    coefs.setRandom();
    auto φ = helper.getScalarTrialSpace(tbasis,0,"φ");
    auto ψ = helper.getScalarTrialSpace(tbasis,1,"ψ");
    auto u = helper.getSolution(φ,coefs,"u");
    auto v = helper.getSolution(ψ,coefs,"v");

    // TODO: Fix this
    eval(helper, u, points);
    eval(helper, v, points);

    // Get variation of expressions
    // w.r.t φ
    print(variation(α, φ));
    print(variation(β, φ));
    print(variation(γ, φ));
    print(variation(δ, φ));
    print(variation(φ, φ));
    print(variation(ψ, φ));
    print(variation(u, φ));
    print(variation(v, φ));
    // w.r.t ψ
    print(variation(α, ψ));
    print(variation(β, ψ));
    print(variation(γ, ψ));
    print(variation(δ, ψ));
    print(variation(φ, ψ));
    print(variation(ψ, ψ));
    print(variation(u, ψ));
    print(variation(v, ψ));

    print(variation(grad(u), φ));
    print(variation(grad(u), ψ));
    print(variation(grad(v), φ));
    print(variation(grad(v), ψ));
    print(variation(u+v, φ));
    print(variation(u+v, ψ));
    print(variation(u-v, φ));
    print(variation(u-v, ψ));
    print(variation(u*v, φ));
    print(variation(u*v, ψ));
    print(variation(u/v, φ));
    print(variation(u/v, ψ));
    print(variation(grad(u+v), ψ));
    print(variation(div(grad(u+v)), ψ));
    print(variation(grad(u*v), ψ));



    return EXIT_SUCCESS;
}
