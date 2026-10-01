/** @file gsEigenAddons_test.cpp

    @brief Numerical and compile-time coverage for the src/gsEigen/ wrappers
    (KroneckerProduct, MatrixFunctions, IterativeSolvers, SparseExtra): each
    wrapper's module symbols land in namespace gsEigen and compute the
    expected result on G+Smo matrix types.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"
#include <gsEigen/KroneckerProduct.h>
#include <gsEigen/MatrixFunctions.h>
#include <gsEigen/IterativeSolvers.h>
#include <gsEigen/SparseExtra.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <type_traits>
#include <utility>

namespace {
typedef gsMatrix<real_t>                                                  NsDMat;
typedef gsSparseMatrix<real_t>                                            NsSMat;
typedef gsEigen::Matrix<real_t, gsEigen::Dynamic, gsEigen::Dynamic>       NsEDMat;
typedef gsEigen::SparseMatrix<real_t, 0, index_t>                         NsESMat;
}

// Each assert below names a symbol that only exists once its module has been
// included through the gsEigen:: rename: kroneckerProduct/GMRES/MINRES/
// saveMarket/loadMarket are declared solely inside their unsupported module,
// and the two sizeof() checks need the complete MatrixExponentialReturnValue
// / MatrixSquareRootReturnValue types that only MatrixFunctions.h defines
// (MatrixBase::exp()/sqrt() themselves are declared in Core, so decltype on
// them alone would compile even without the module).
static_assert(std::is_same<decltype(gsEigen::kroneckerProduct(std::declval<const NsDMat&>(),
                                                              std::declval<const NsDMat&>()).rows()),
                           gsEigen::Index>::value, "gsEigen::kroneckerProduct (dense) missing");
static_assert(std::is_same<decltype(gsEigen::kroneckerProduct(std::declval<const NsSMat&>(),
                                                              std::declval<const NsSMat&>()).rows()),
                           gsEigen::Index>::value, "gsEigen::kroneckerProduct (sparse) missing");
static_assert(sizeof(gsEigen::MatrixExponentialReturnValue<NsEDMat>) > 0, "gsEigen matrix exp not defined");
static_assert(sizeof(gsEigen::MatrixSquareRootReturnValue<NsEDMat>)  > 0, "gsEigen matrix sqrt not defined");
static_assert(std::is_base_of<gsEigen::IterativeSolverBase<gsEigen::GMRES<NsESMat> >,
                              gsEigen::GMRES<NsESMat> >::value, "gsEigen::GMRES missing");
static_assert(std::is_base_of<gsEigen::IterativeSolverBase<gsEigen::MINRES<NsESMat> >,
                              gsEigen::MINRES<NsESMat> >::value, "gsEigen::MINRES missing");
static_assert(std::is_same<decltype(gsEigen::saveMarket(std::declval<const NsSMat&>(), std::declval<const std::string&>())),
                           bool>::value, "gsEigen::saveMarket missing");
static_assert(std::is_same<decltype(gsEigen::loadMarket(std::declval<NsSMat&>(), std::declval<const std::string&>())),
                           bool>::value, "gsEigen::loadMarket missing");

namespace {

// Shared 2x2 / 2x3 operands and their reference Kronecker product, used by
// both KroneckerProduct tests so the dense and sparse paths check the same
// numbers without either test depending on the other's state.
gsMatrix<real_t> kroneckerA()
{
    gsMatrix<real_t> A(2,2);
    A << 1, -2,
         3,  0;
    return A;
}

gsMatrix<real_t> kroneckerB()
{
    gsMatrix<real_t> B(2,3);
    B << 0, 5, -1,
         6, 0,  7;
    return B;
}

gsMatrix<real_t> kroneckerRef(const gsMatrix<real_t> & A, const gsMatrix<real_t> & B)
{
    gsMatrix<real_t> Kref(4,6);
    for (index_t i = 0; i < 2; ++i)
        for (index_t j = 0; j < 2; ++j)
            Kref.block(2*i, 3*j, 2, 3) = A(i,j) * B;
    return Kref;
}

}

SUITE(gsEigenAddons_test)
{

TEST(KroneckerProduct_dense)
{
    const gsMatrix<real_t> A = kroneckerA();
    const gsMatrix<real_t> B = kroneckerB();
    const gsMatrix<real_t> Kref = kroneckerRef(A, B);

    gsMatrix<real_t> K = gsEigen::kroneckerProduct(A, B);

    CHECK_EQUAL(4, K.rows());
    CHECK_EQUAL(6, K.cols());
    // Every entry of A (x) B is a product of two small integers, which
    // binary64 represents exactly (integers up to 2^53), so the comparison
    // needs no tolerance. B (x) A has the same 4x6 shape, so only comparing
    // entries (not just size) catches a swap of the two factors.
    CHECK_EQUAL(real_t(0), (K - Kref).cwiseAbs().maxCoeff());
}

TEST(KroneckerProduct_sparse)
{
    const gsMatrix<real_t> A = kroneckerA();
    const gsMatrix<real_t> B = kroneckerB();
    const gsMatrix<real_t> Kref = kroneckerRef(A, B);

    const gsSparseMatrix<real_t> As = A.sparseView();
    const gsSparseMatrix<real_t> Bs = B.sparseView();
    gsSparseMatrix<real_t> Ks = gsEigen::kroneckerProduct(As, Bs);

    CHECK_EQUAL(4, Ks.rows());
    CHECK_EQUAL(6, Ks.cols());
    // nnz(A)=3 (the (1,1) entry of A is a structural zero), nnz(B)=4, and the
    // Kronecker product stores exactly nnz(A)*nnz(B)=12 entries: no dense
    // fill-in for the structural zero.
    CHECK_EQUAL(12, Ks.nonZeros());
    CHECK_EQUAL(real_t(0), (gsMatrix<real_t>(Ks.toDense()) - Kref).cwiseAbs().maxCoeff());
}

TEST(MatrixFunctions_exp)
{
    const real_t tol = 1e-12;

    // Diagonal case: exp(diag(d)) = diag(exp(d)) exactly. A coefficient-wise
    // exp() applied to the whole matrix would instead put a spurious 1 in
    // every off-diagonal zero, so this also rejects that wrong formula.
    gsMatrix<real_t> D = gsMatrix<real_t>::Zero(3,3);
    D(0,0) = 0; D(1,1) = 1; D(2,2) = -2;
    gsMatrix<real_t> Dref = gsMatrix<real_t>::Zero(3,3);
    Dref(0,0) = std::exp(0.); Dref(1,1) = std::exp(1.); Dref(2,2) = std::exp(-2.);

    // E5 computes exp() via scaling-and-squaring Pade approximation, backward
    // stable at unit roundoff u ~ 1.1e-16 (matrix_exp_computeUV<...,double>,
    // MatrixExponential.h). Forward error ~ kappa_exp * u, and kappa_exp =
    // ||D||_2 = 2 for this normal D, so tol = 1e-12 (~1e4*u) keeps >= 2 orders
    // of margin; a formula missing a series term errs by O(1), 12 orders
    // above tol.
    gsMatrix<real_t> ED = D.exp();
    CHECK((ED - Dref).cwiseAbs().maxCoeff() <= tol * Dref.cwiseAbs().maxCoeff());

    // Nilpotent case: N^3 = 0, so exp(N) = I + N + N^2/2 exactly, with N^2
    // having the single nonzero entry (0,2) = 1*3 = 3.
    gsMatrix<real_t> N(3,3);
    N << 0, 1, 2,
         0, 0, 3,
         0, 0, 0;
    gsMatrix<real_t> Nref(3,3);
    Nref << 1, 1, 3.5,
            0, 1, 3,
            0, 0, 1;

    gsMatrix<real_t> EN = N.exp();
    CHECK((EN - Nref).cwiseAbs().maxCoeff() <= tol * Nref.cwiseAbs().maxCoeff());
}

TEST(MatrixFunctions_sqrt)
{
    const real_t tol = 1e-12;

    // Symmetric with row sums bounding the spectrum, by Gershgorin, to
    // [1,5]: S is SPD with kappa_2(S) <= 5.
    gsMatrix<real_t> S(3,3);
    S << 4, 1, 0,
         1, 3, 1,
         0, 1, 2;

    gsMatrix<real_t> R = S.sqrt();

    // E5 computes the real square root via a real Schur decomposition plus a
    // quasi-triangular recurrence (MatrixSquareRoot.h). For normal S the
    // residual ||R^2-S||/||S|| is O(n*u) (Higham's alpha(R) = ||R||^2/||S||
    // = 1 here), so tol = 1e-12 keeps ample margin. A coefficient-wise sqrt
    // would fail the first check below by O(1) from the spurious
    // off-diagonal ones it introduces.
    CHECK((R*R - S).cwiseAbs().maxCoeff() <= tol * S.cwiseAbs().maxCoeff());
    CHECK((R - R.transpose()).cwiseAbs().maxCoeff() <= tol * R.cwiseAbs().maxCoeff());
    // The principal square root of an SPD matrix is itself SPD, which rules
    // out the other (non-principal) roots that also satisfy R^2 = S.
    CHECK_EQUAL(gsEigen::Success, R.llt().info());

    gsMatrix<real_t> Ddiag = gsMatrix<real_t>::Zero(3,3);
    Ddiag(0,0) = 4; Ddiag(1,1) = 9; Ddiag(2,2) = 16;
    gsMatrix<real_t> Dref = gsMatrix<real_t>::Zero(3,3);
    Dref(0,0) = 2; Dref(1,1) = 3; Dref(2,2) = 4;

    gsMatrix<real_t> Rd = Ddiag.sqrt();
    CHECK((Rd - Dref).cwiseAbs().maxCoeff() <= tol * Dref.cwiseAbs().maxCoeff());
}

TEST(IterativeSolvers_MINRES)
{
    const index_t n = 20;
    gsSparseMatrix<real_t> A(n,n);
    A.reserve(gsVector<index_t>::Constant(n,3));
    for (index_t i = 0; i < n; ++i)
    {
        A.insert(i,i) = 2;
        if (i > 0)   A.insert(i,i-1) = -1;
        if (i < n-1) A.insert(i,i+1) = -1;
    }
    A.makeCompressed();

    gsVector<real_t> xref(n);
    for (index_t i = 0; i < n; ++i)
        xref(i) = std::sin(static_cast<real_t>(i+1));
    gsVector<real_t> b = A * xref;

    // Default UpLo=Lower, IdentityPreconditioner (MINRES.h).
    gsEigen::MINRES<gsEigen::SparseMatrix<real_t,0,index_t> > s;
    s.setTolerance(1e-10);
    s.setMaxIterations(200);
    s.compute(A);
    gsVector<real_t> x = s.solve(b);

    CHECK_EQUAL(gsEigen::Success, s.info());
    // x0 = 0 is not the solution, so a solver that did nothing would fail
    // this rather than the residual check.
    CHECK(s.iterations() > 0);
    CHECK(s.iterations() <= 200);
    // kappa_2(A) = lambda_max/lambda_min with lambda_k = 2 - 2cos(k*pi/(n+1))
    // for this 1D Laplacian, ~178 for n=20. The recurrence residual MINRES
    // monitors internally drifts from the explicit b-Ax by O(kappa*u) ~ 2e-14
    // relative, so checking the explicit residual at 10x the solver
    // tolerance (1e-9 vs 1e-10) has ample margin.
    CHECK((b - A*x).norm() <= 1e-9 * b.norm());
}

TEST(IterativeSolvers_GMRES)
{
    const index_t n = 20;
    gsSparseMatrix<real_t> A(n,n);
    A.reserve(gsVector<index_t>::Constant(n,3));
    for (index_t i = 0; i < n; ++i)
    {
        A.insert(i,i) = 3;
        if (i > 0)   A.insert(i,i-1) = -1.3;
        if (i < n-1) A.insert(i,i+1) = -0.7;
    }
    A.makeCompressed();

    gsVector<real_t> xref(n);
    for (index_t i = 0; i < n; ++i)
        xref(i) = std::sin(static_cast<real_t>(i+1));
    gsVector<real_t> b = A * xref;

    // Default DiagonalPreconditioner, default restart 30 >= n so this is
    // full (unrestarted) GMRES (GMRES.h). A is strictly diagonally dominant
    // (|3| > 1.3+0.7), hence nonsingular, and genuinely nonsymmetric: it has
    // no closed-form condition number since it is non-normal.
    gsEigen::GMRES<gsEigen::SparseMatrix<real_t,0,index_t> > s;
    s.setTolerance(1e-10);
    s.setMaxIterations(200);
    s.compute(A);
    gsVector<real_t> x = s.solve(b);

    CHECK_EQUAL(gsEigen::Success, s.info());
    CHECK(s.iterations() > 0);
    CHECK(s.iterations() <= 200);
    // The diagonal is constant (3), so Jacobi preconditioning is a uniform
    // scaling and the preconditioned relative residual equals the true one;
    // 1e-9 keeps a 10x margin over the 1e-10 solver tolerance.
    CHECK((b - A*x).norm() <= 1e-9 * b.norm());
}

TEST(SparseExtra_market_roundtrip)
{
    gsSparseMatrix<real_t> M(4,5);
    M.reserve(5);
    M.insert(0,0) = 1.0/3.0;
    M.insert(1,2) = std::sqrt(2.0);
    M.insert(3,1) = 0.1;
    M.insert(2,3) = 1e-300;
    M.insert(3,4) = -6.02214076e23;
    M.makeCompressed();

    const std::string fn = gsFileManager::getTempPath() + "gsEigenAddons_test_market.mtx";
    std::remove(fn.c_str());

    // saveMarket writes scientific notation with digits10+2 = 17 significant
    // digits, binary64's max_digits10 (MarketIO.h), and loadMarket parses
    // with %lg, so the round trip is exact for every entry. sqrt(2) needs
    // all 17 digits to round-trip exactly, 1/3 needs 16 and -6.02214076e23
    // needs 9: these three exercise the written precision. 0.1 and 1e-300
    // would round-trip even at 6 digits; they are present only for
    // exponent-range coverage.
    bool okS = gsEigen::saveMarket(M, fn);
    gsSparseMatrix<real_t> L;
    bool okL = gsEigen::loadMarket(L, fn);

    CHECK(okS);
    CHECK(okL);
    CHECK_EQUAL(4, L.rows());
    CHECK_EQUAL(5, L.cols());
    CHECK_EQUAL(M.nonZeros(), L.nonZeros());
    CHECK_EQUAL(real_t(0), (gsMatrix<real_t>(L.toDense()) - gsMatrix<real_t>(M.toDense())).cwiseAbs().maxCoeff());

    CHECK_EQUAL(0, std::remove(fn.c_str()));
    CHECK(!gsFileManager::fileExists(fn));
}

TEST(wrappers_do_not_leak_rename_macro)
{
#ifdef Eigen
    CHECK(false); // a wrapper left the Eigen->gsEigen renaming macro defined
#else
    CHECK(true);
#endif
}

}
