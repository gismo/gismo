/** @file BiCGSTAB.h

    @brief BiCGSTAB solver gsEigen::gsBiCGSTAB with a relative stopping
    criterion |b - A x| <= tol |b|.

    Modified copy of Eigen 5.0.0 Eigen/src/IterativeLinearSolvers/BiCGSTAB.h.
    The only change in the algorithm is upstream Eigen commit 93ff3888
    (libeigen/eigen!2118): the stopping threshold is scaled by |b|.

    Eigen 5.0.0 and 5.0.1 stop BiCGSTAB on the absolute residual
    |b - A x| <= tol, although tolerance() and error() are documented as, and
    error() returns, the relative residual |b - A x| / |b| (regression of
    Eigen commit 809d266b, which dropped the factor |b| from the former test
    |r|^2 <= tol^2 |b|^2). For |b| < 1 the solver then stops before error()
    meets tolerance() and reports NoConvergence; for |b| > 1 it iterates past
    the requested accuracy, up to maxIterations().

    The class and the iteration function are renamed (gsBiCGSTAB,
    internal::gs_bicgstab) so that they coexist with stock gsEigen::BiCGSTAB,
    which <Eigen/Sparse> always declares.

    TODO: remove this file and point the BiCGSTAB* typedefs in
    gsCore/gsLinearAlgebra.h back to gsEigen::BiCGSTAB once the minimum
    supported Eigen release contains 93ff3888.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Original Eigen authors:
    Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
    Copyright (C) 2012 Désiré Nuentsa-Wakam <desire.nuentsa_wakam@inria.fr>
*/

#pragma once

// gsCore/gsLinearAlgebra.h includes this header and names gsBiCGSTAB in its
// typedefs. When this header is included first, that nested include is
// skipped, so the template must already be declared here.
namespace gsEigen {
template <typename MatrixType_, typename Preconditioner_> class gsBiCGSTAB;
}

#include <gsCore/gsLinearAlgebra.h>

namespace gsEigen {

namespace internal {

/** \internal Low-level bi conjugate gradient stabilized algorithm
 * \param mat The matrix A
 * \param rhs The right hand side vector b
 * \param x On input and initial solution, on output the computed solution.
 * \param precond A preconditioner being able to efficiently solve for an
 *                approximation of Ax=b (regardless of b)
 * \param iters On input the max number of iteration, on output the number of performed iterations.
 * \param tol_error On input the relative tolerance |b - A x| / |b|, on output an estimation of the relative error.
 * \return false in the case of numerical issue, for example a break down of BiCGSTAB.
 */
template <typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
bool gs_bicgstab(const MatrixType& mat, const Rhs& rhs, Dest& x, const Preconditioner& precond, Index& iters,
                 typename Dest::RealScalar& tol_error) {
  typedef typename Dest::RealScalar RealScalar;
  typedef typename Dest::Scalar Scalar;
  typedef Matrix<Scalar, Dynamic, 1> VectorType;
  Index maxIters = iters;

  Index n = mat.cols();
  VectorType r = rhs - mat * x;
  VectorType r0 = r;

  RealScalar r0_norm = r0.stableNorm();
  RealScalar r_norm = r0_norm;
  RealScalar rhs_norm = rhs.stableNorm();
  if (rhs_norm == 0) {
    x.setZero();
    return true;
  }

  // The tolerance is relative: compare |r| against tol |b| (Eigen 93ff3888).
  RealScalar tol = tol_error * rhs_norm;

  Scalar rho(1);
  Scalar alpha(0);
  Scalar w(1);

  VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
  VectorType y(n), z(n);
  VectorType kt(n), ks(n);

  VectorType s(n), t(n);

  RealScalar eps = NumTraits<Scalar>::epsilon();
  Index i = 0;
  Index restarts = 0;

  while (r_norm > tol && i < maxIters) {
    Scalar rho_old = rho;
    rho = r0.dot(r);
    if (numext::abs(rho) / numext::maxi(r0_norm, r_norm) < eps * numext::mini(r0_norm, r_norm)) {
      // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
      // Let's restart with a new r0:
      r = rhs - mat * x;
      r0 = r;
      rho = r.squaredNorm();
      r0_norm = r.stableNorm();
      alpha = Scalar(0);
      w = Scalar(1);
      if (restarts++ == 0) i = 0;
    }
    Scalar beta = (rho / rho_old) * (alpha / w);
    p = r + beta * (p - w * v);

    y = precond.solve(p);

    v.noalias() = mat * y;
    Scalar theta = r0.dot(v);
    // For small angles ∠(r0, v) < eps, random restart.
    RealScalar v_norm = v.stableNorm();
    if (numext::abs(theta) / numext::maxi(r0_norm, v_norm) < eps * numext::mini(r0_norm, v_norm)) {
      r = rhs - mat * x;
      r0.setRandom();
      r0_norm = r0.stableNorm();
      rho = Scalar(1);
      alpha = Scalar(0);
      w = Scalar(1);
      if (restarts++ == 0) i = 0;
      continue;
    }
    alpha = rho / theta;
    s = r - alpha * v;

    z = precond.solve(s);
    t.noalias() = mat * z;

    RealScalar tmp = t.squaredNorm();
    if (tmp > RealScalar(0)) {
      w = t.dot(s) / tmp;
    } else {
      w = Scalar(0);
    }
    x += alpha * y + w * z;
    r = s - w * t;
    r_norm = r.stableNorm();
    ++i;
  }

  tol_error = r_norm / rhs_norm;
  iters = i;
  return true;
}

}  // namespace internal

template <typename MatrixType_, typename Preconditioner_ = DiagonalPreconditioner<typename MatrixType_::Scalar> >
class gsBiCGSTAB;

namespace internal {

template <typename MatrixType_, typename Preconditioner_>
struct traits<gsBiCGSTAB<MatrixType_, Preconditioner_> > {
  typedef MatrixType_ MatrixType;
  typedef Preconditioner_ Preconditioner;
};

}  // namespace internal

/** \brief A bi conjugate gradient stabilized solver for sparse square problems
 *
 * Drop-in replacement for gsEigen::BiCGSTAB with the same interface. The
 * tolerance, set by setTolerance() and NumTraits<Scalar>::epsilon() by
 * default, bounds the relative residual |Ax-b|/|b|, which error() returns.
 *
 * \tparam MatrixType_ the type of the sparse matrix A, can be a dense or a sparse matrix.
 * \tparam Preconditioner_ the type of the preconditioner. Default is DiagonalPreconditioner
 *
 * The maximal number of iterations defaults to twice the size of the problem
 * (see IterativeSolverBase::maxIterations()).
 * By default the iterations start with x=0 as an initial guess of the solution;
 * solveWithGuess() sets the start.
 */
template <typename MatrixType_, typename Preconditioner_>
class gsBiCGSTAB : public IterativeSolverBase<gsBiCGSTAB<MatrixType_, Preconditioner_> > {
  typedef IterativeSolverBase<gsBiCGSTAB> Base;
  using Base::m_error;
  using Base::m_info;
  using Base::m_isInitialized;
  using Base::m_iterations;
  using Base::matrix;

 public:
  typedef MatrixType_ MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::RealScalar RealScalar;
  typedef Preconditioner_ Preconditioner;

 public:
  /** Default constructor. */
  gsBiCGSTAB() : Base() {}

  /** Initialize the solver with matrix \a A for further \c Ax=b solving.
   *
   * \warning this class stores a reference to the matrix A as well as some
   * precomputed values that depend on it. Therefore, if \a A is changed
   * this class becomes invalid. Call compute() to update it with the new
   * matrix A, or modify a copy of A.
   */
  template <typename MatrixDerived>
  explicit gsBiCGSTAB(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

  ~gsBiCGSTAB() {}

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const {
    m_iterations = Base::maxIterations();
    m_error = Base::m_tolerance;

    bool ret = internal::gs_bicgstab(matrix(), b, x, Base::m_preconditioner, m_iterations, m_error);

    m_info = (!ret) ? NumericalIssue : m_error <= Base::m_tolerance ? Success : NoConvergence;
  }
};

}  // namespace gsEigen
