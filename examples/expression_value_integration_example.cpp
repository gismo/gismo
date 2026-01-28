/** @file expression_value_integration_example.cpp

    @brief Example showing how ExpressionResult integrates with the expression system

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

#include <gsNewExpressions/NewExpressions.h>
#include <gsNewExpressions/ExpressionHelper.h>

using namespace gismo;
using namespace Expr;

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Example: ExpressionResult integration with expression system.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "=== ExpressionResult Integration Example ===\n\n";

    // Setup basis
    gsKnotVector<> kv(0., 1., 2, 2);  // Simple 3-basis setup
    gsTensorBSplineBasis<2> tbasis(kv, kv);
    
    gsInfo << "Basis: " << tbasis.size() << " basis functions\n\n";

    // Create expression helper
    ExpressionHelper<real_t> helper;
    gsMatrix<real_t> points(2, 1);
    points << 0.5, 0.5;
    helper.points() = points;

    // Define test and trial spaces
    auto φ = helper.getScalarTestSpace(tbasis, 0, "φ");
    auto ψ = helper.getScalarTrialSpace(tbasis, 1, "ψ");

    // Example 1: Scalar expression (e.g., constant or function)
    gsInfo << "Example 1: Scalar expression\n";
    auto α = helper.getConstant(2.0, "α");
    gsInfo << "  Expression: " << α << "\n";
    gsInfo << "  Space type: " << α.space() << " (None=0)\n";
    gsInfo << "  Expected cardinality: (1, 1)\n";
    
    ExpressionResult<real_t> scalar_result = makeExpressionResult<real_t>(α.space());
    gsInfo << "  ExpressionResult cardinality: (" << scalar_result.rowCardinality() 
           << ", " << scalar_result.colCardinality() << ")\n\n";

    // Example 2: Test space expression (e.g., φ or ∇φ)
    gsInfo << "Example 2: Test space expression\n";
    gsInfo << "  Expression: " << φ << "\n";
    gsInfo << "  Space type: " << φ.space() << " (Test=1)\n";
    gsInfo << "  Number of test basis functions: " << tbasis.size() << "\n";
    gsInfo << "  Expected cardinality: (" << tbasis.size() << ", 1)\n";
    
    ExpressionResult<real_t> test_result = makeExpressionResult<real_t>(φ.space(), tbasis.size());
    gsInfo << "  ExpressionResult cardinality: (" << test_result.rowCardinality() 
           << ", " << test_result.colCardinality() << ")\n";
    gsInfo << "  Can access test_result(i, 0) for i = 0.." << (tbasis.size()-1) << "\n\n";

    // Example 3: Trial space expression (e.g., ψ or ∇ψ)
    gsInfo << "Example 3: Trial space expression\n";
    gsInfo << "  Expression: " << ψ << "\n";
    gsInfo << "  Space type: " << ψ.space() << " (Trial=2)\n";
    gsInfo << "  Number of trial basis functions: " << tbasis.size() << "\n";
    gsInfo << "  Expected cardinality: (1, " << tbasis.size() << ")\n";
    
    ExpressionResult<real_t> trial_result = makeExpressionResult<real_t>(ψ.space(), 1, tbasis.size());
    gsInfo << "  ExpressionResult cardinality: (" << trial_result.rowCardinality() 
           << ", " << trial_result.colCardinality() << ")\n";
    gsInfo << "  Can access trial_result(0, j) for j = 0.." << (tbasis.size()-1) << "\n\n";

    // Example 4: Bilinear expression (e.g., φ·ψ, φ*ψ, ∇φ·∇ψ)
    gsInfo << "Example 4: Bilinear expression (Test + Trial = Both)\n";
    auto bilinear_expr = φ + ψ;
    gsInfo << "  Expression: " << bilinear_expr << "\n";
    gsInfo << "  Space type: " << bilinear_expr.space() << " (Both=3)\n";
    gsInfo << "  Expected cardinality: (" << tbasis.size() << ", " << tbasis.size() << ")\n";
    
    ExpressionResult<real_t> bilinear_result = makeExpressionResult<real_t>(
        bilinear_expr.space(), tbasis.size(), tbasis.size());
    gsInfo << "  ExpressionResult cardinality: (" << bilinear_result.rowCardinality() 
           << ", " << bilinear_result.colCardinality() << ")\n";
    gsInfo << "  Can access bilinear_result(i, j) for i,j = 0.." << (tbasis.size()-1) << "\n";
    gsInfo << "  Total matrices stored: " << bilinear_result.size() << "\n\n";

    // Example 5: Demonstrating matrix dimensions
    gsInfo << "Example 5: Matrix dimensions for different orders\n";
    
    // Scalar (order 0): returns 1×1 matrix (or n_quad_points×1 for multiple points)
    gsInfo << "  Scalar expression (order 0): matrix is (n_points × 1)\n";
    ExpressionResult<real_t> val0 = makeExpressionResult<real_t>(SpaceType::None, 1, 1, 1, 1);
    gsInfo << "    Example: (" << val0(0,0).rows() << " × " << val0(0,0).cols() << ")\n";
    
    // Vector (order 1): returns d×1 matrix per basis function
    gsInfo << "  Vector expression (order 1): matrix is (d × n_points) where d=dimension\n";
    ExpressionResult<real_t> val1 = makeExpressionResult<real_t>(SpaceType::Test, 3, 1, 2, 1);
    gsInfo << "    Example for 2D: (" << val1(0,0).rows() << " × " << val1(0,0).cols() << ")\n";
    
    // Matrix (order 2): returns d×d matrix per basis function
    gsInfo << "  Matrix expression (order 2): matrix is (d × d)\n";
    ExpressionResult<real_t> val2 = makeExpressionResult<real_t>(SpaceType::Test, 3, 1, 2, 2);
    gsInfo << "    Example for 2D: (" << val2(0,0).rows() << " × " << val2(0,0).cols() << ")\n\n";

    // Example 6: Assembly use case
    gsInfo << "Example 6: Typical assembly pattern\n";
    gsInfo << "  For a bilinear form like ∫ ∇φ·∇ψ dx:\n";
    gsInfo << "  1. Evaluate ∇φ at quadrature points → ExpressionResult with cardinality (N, 1)\n";
    gsInfo << "  2. Evaluate ∇ψ at quadrature points → ExpressionResult with cardinality (1, M)\n";
    gsInfo << "  3. Compute inner product → ExpressionResult with cardinality (N, M)\n";
    gsInfo << "  4. Each entry (i,j) represents contribution to global matrix entry (I[i], J[j])\n\n";

    gsInfo << "=== Example completed successfully! ===\n";

    return EXIT_SUCCESS;
}
