#include <gismo.h>
#include <gsNewExpressions/NewExpressions.h>
#include <gsNewExpressions/ExpressionHelper.h>
#include <gsNewExpressions/ExprEvaluator.h>

using namespace gismo;



int main(int argc, char* argv[])
{
    gsCmdLine cmd("new_expressions_evaluator_example.cpp");
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // 1. Setup
    ExpressionHelper<real_t> helper;
    
    // 2. Create expressions
    gsFunctionExpr<> fun("x*y", 2);
    auto f = helper.getScalarFunction(fun, "f");
    auto c = helper.getConstant(2.0, "c");
    auto expr = c * f;

    // 3. Use ExprEvaluator
    ExprEvaluator<real_t> evaluator(helper);
    gsVector<real_t> pt(2);
    pt << 0.5, 0.5;
    auto result_val = evaluator.eval(expr, pt);

    gsInfo << "Evaluating '" << expr << "' at " << pt.transpose() << ": ";
    gsInfo << result_val(0,0) << "\n";
    gsInfo << "Expected: 0.5\n";

    // 4. Test computeIntegral
    gsInfo << "\nComputing integral of 'x*y' over a square domain...\n";
    gsFunctionExpr<> integral_fun("x*y", 2);
    auto integral_expr = helper.getScalarFunction(integral_fun, "integral_expr");
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare()); // A unit square
    
    gsMatrix<real_t> integral_result = evaluator.computeIntegral(integral_expr, *mp.domain());
    
    gsInfo << "Integral of 'x*y' over unit square: " << integral_result(0,0) << "\n";
    gsInfo << "Expected: 0.25 (since integral from 0 to 1 of x dx = 0.5, and y dy = 0.5, product = 0.25)\n";
    
    if (std::abs(integral_result(0,0) - 0.25) < 1e-9) {
        gsInfo << "Integral test PASSED.\n";
    } else {
        gsInfo << "Integral test FAILED.\n";
    }

    gsInfo << "Done.\n";

    return 0;
}
