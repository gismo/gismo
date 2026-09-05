/** @file gsEigenIntegration_test.cpp

    @brief Tests for Eigen 5 integration: core gsMatrix/gsVector types,
    Eigen's Tensor module, and Eigen's own AutoDiff umbrella header.

    Eigen's unsupported modules (Tensor, AutoDiff) can be included directly
    in the same translation unit as gismo.h now that the legacy
    '#define Eigen gsEigen' alias has been removed from the codebase.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
*/

#include "gismo_unittest.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/AutoDiff>

using namespace gismo;

SUITE(gsEigenIntegration_test)
{

// -------------------------------------------------------------------------
// Basic Eigen integration via gsMatrix / gsVector
// -------------------------------------------------------------------------

TEST(Matrix_multiply)
{
    gsMatrix<> A(2,2), B(2,2);
    A << 1, 2,
         3, 4;
    B << 5, 6,
         7, 8;
    gsMatrix<> C = A * B;
    CHECK_CLOSE(19.0, C(0,0), 1e-14);
    CHECK_CLOSE(22.0, C(0,1), 1e-14);
    CHECK_CLOSE(43.0, C(1,0), 1e-14);
    CHECK_CLOSE(50.0, C(1,1), 1e-14);
}

TEST(SVD)
{
    // Use Eigen::MatrixXd directly to avoid cross-namespace evaluator issues
    // in Eigen 5 when passing gismo::gsMatrix to JacobiSVD.
    Eigen::MatrixXd A(3,2);
    A << 1, 0,
         0, 2,
         0, 0;
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
    Eigen::VectorXd sv = svd.singularValues();
    CHECK(sv(0) >= sv(1));
    CHECK(sv(1) >= 0.0);
    CHECK_CLOSE(2.0, sv(0), 1e-14);
    CHECK_CLOSE(1.0, sv(1), 1e-14);
    // Reconstruction: U * diag(sv) * V^T == A
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2,2);
    S(0,0) = sv(0); S(1,1) = sv(1);
    Eigen::MatrixXd Ar = svd.matrixU() * S * svd.matrixV().transpose();
    CHECK_CLOSE(0.0, (A - Ar).norm(), 1e-13);
}

// -------------------------------------------------------------------------
// Eigen AutoDiff — first-order (gradient) using AutoDiffScalar
// -------------------------------------------------------------------------

TEST(AutoDiff_gradient)
{
    // f(x,y) = x^2 + 2*y  =>  df/dx = 2x = 6,  df/dy = 2  at (x,y)=(3,1)
    typedef Eigen::Matrix<double, 2, 1> DerivVec;
    typedef Eigen::AutoDiffScalar<DerivVec> ADS;

    ADS x(3.0, DerivVec::Unit(0));  // value=3, seed=[1,0]
    ADS y(1.0, DerivVec::Unit(1));  // value=1, seed=[0,1]

    ADS f = x * x + ADS(2.0) * y;

    CHECK_CLOSE(11.0, f.value(),           1e-14);
    CHECK_CLOSE( 6.0, f.derivatives()(0), 1e-14);  // df/dx
    CHECK_CLOSE( 2.0, f.derivatives()(1), 1e-14);  // df/dy
}

// -------------------------------------------------------------------------
// Eigen AutoDiff — second-order (Hessian) using nested AutoDiffScalar
// -------------------------------------------------------------------------

TEST(AutoDiff_hessian)
{
    // f(x,y) = x^2 * y  at (x,y)=(3,1)
    //   value      = 9
    //   df/dx      = 2xy = 6,  df/dy = x^2 = 9
    //   d^2f/dx^2  = 2y  = 2,  d^2f/dxdy  = 2x = 6,  d^2f/dy^2 = 0

    typedef Eigen::Matrix<double, 2, 1>        InnerDeriv;
    typedef Eigen::AutoDiffScalar<InnerDeriv>  InnerADS;
    typedef Eigen::Matrix<InnerADS, 2, 1>      OuterDeriv;
    typedef Eigen::AutoDiffScalar<OuterDeriv>  ADS2;

    // Construct x: value=3, first-order seed e_0, second-order seeds zero
    ADS2 x, y;
    x.value()              = InnerADS(3.0, InnerDeriv::Unit(0));
    x.derivatives()(0)     = InnerADS(1.0, InnerDeriv::Zero());
    x.derivatives()(1)     = InnerADS(0.0, InnerDeriv::Zero());

    // Construct y: value=1, first-order seed e_1, second-order seeds zero
    y.value()              = InnerADS(1.0, InnerDeriv::Unit(1));
    y.derivatives()(0)     = InnerADS(0.0, InnerDeriv::Zero());
    y.derivatives()(1)     = InnerADS(1.0, InnerDeriv::Zero());

    ADS2 f = x * x * y;

    // value
    CHECK_CLOSE(9.0, f.value().value(),              1e-14);
    // gradient (read from inner derivatives of the outer value)
    CHECK_CLOSE(6.0, f.value().derivatives()(0),    1e-14);  // df/dx
    CHECK_CLOSE(9.0, f.value().derivatives()(1),    1e-14);  // df/dy
    // Hessian row 0: d^2f/dx^2 and d^2f/dxdy
    CHECK_CLOSE(2.0, f.derivatives()(0).derivatives()(0), 1e-14); // H[0,0]
    CHECK_CLOSE(6.0, f.derivatives()(0).derivatives()(1), 1e-14); // H[0,1]
    // Hessian row 1: d^2f/dydx and d^2f/dy^2
    CHECK_CLOSE(6.0, f.derivatives()(1).derivatives()(0), 1e-14); // H[1,0]
    CHECK_CLOSE(0.0, f.derivatives()(1).derivatives()(1), 1e-14); // H[1,1]
}

// -------------------------------------------------------------------------
// Eigen Tensor module — now includable alongside gismo.h since the
// '#define Eigen gsEigen' alias has been removed.
// -------------------------------------------------------------------------

TEST(Tensor_basic)
{
    // Create a rank-3 tensor of shape (2, 3, 4) and zero-initialise it.
    Eigen::Tensor<double, 3> t(2, 3, 4);
    t.setZero();

    // Total element count must be 2*3*4 = 24.
    CHECK_EQUAL(24, static_cast<int>(t.size()));

    // All entries should be exactly 0.
    bool all_zero = true;
    for (int i = 0; i < 2 && all_zero; ++i)
        for (int j = 0; j < 3 && all_zero; ++j)
            for (int k = 0; k < 4 && all_zero; ++k)
                if (t(i, j, k) != 0.0) all_zero = false;
    CHECK(all_zero);

    // A simple element-wise assignment and read-back.
    t(1, 2, 3) = 42.0;
    CHECK_CLOSE(42.0, t(1, 2, 3), 1e-14);
}

// -------------------------------------------------------------------------
// Eigen AutoDiff umbrella header — forward-mode AD using Eigen's own
// AutoDiffScalar type (not gismo's DScalar wrapper).
// -------------------------------------------------------------------------

TEST(AutoDiff_umbrella)
{
    // Derivative of f(x) = x^3 at x = 2 is f'(2) = 3*x^2 = 12.
    typedef Eigen::Matrix<double, 1, 1> DerivVec;
    typedef Eigen::AutoDiffScalar<DerivVec> ADS;

    ADS x;
    x.value() = 2.0;
    x.derivatives() = DerivVec::Ones(); // dx/dx = 1

    ADS f = x * x * x;  // f(x) = x^3

    CHECK_CLOSE( 8.0, f.value(),           1e-14);  // 2^3 = 8
    CHECK_CLOSE(12.0, f.derivatives()(0),  1e-14);  // 3*2^2 = 12

    // Verify linearity: d(f + f)/dx = 2 * df/dx
    ADS g = f + f;
    CHECK_CLOSE(16.0, g.value(),           1e-14);
    CHECK_CLOSE(24.0, g.derivatives()(0),  1e-14);
}

} // SUITE
