/** @file gsBarrierCore.hpp

    @brief This software facilitates the creation of analysis-suitable
    parameterizations from given boundary representations. Serving as a
    reference implementation, it embodies the methods and concepts detailed
    in Ye Ji's doctoral research. Here, optimization-based (barrier, penalty)
    methods and PDE-based methods are provided. Please refer to the
    implementation for the relevant references.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsConstantFunction.h>

namespace gismo {

/// Garanzha's regularization of the Jacobian determinant D:
/// chi(D,eps) = (D + sqrt(D^2 + eps^2))/2 > 0, with chi -> max(D,0) for
/// eps -> 0. See V. Garanzha, I. Kaporin, Regularization of the barrier
/// variational method, Comput. Math. Math. Phys. 39 (1999), and Garanzha et
/// al., Foldover-free maps in 50 lines of code, ACM TOG 40(4), 2021.
template<class T>
inline T garanzhaChi(const T D, const T eps) {
  const T root = math::sqrt(D * D + eps * eps);
  // For D < 0, D + root suffers from cancellation; use the equivalent form
  // (D + root) = eps^2 / (root - D)
  return (D < 0 ? 0.5 * eps * eps / (root - D) : 0.5 * (D + root));
}

/// Cofactor matrix cof(J) = det(J) J^{-T} of a 2x2 or 3x3 matrix; it is
/// well-defined also for singular J. The rows of cof(J) are the derivatives
/// of det(J) with respect to the rows of J.
template<class MatJ, class MatC>
inline void cofactorMatrix(const MatJ &J, MatC &cof) {
  const index_t dim = J.rows();
  cof.resize(dim, dim);
  if (2 == dim) {
    cof(0, 0) = J(1, 1);
    cof(0, 1) = -J(1, 0);
    cof(1, 0) = -J(0, 1);
    cof(1, 1) = J(0, 0);
  } else {
    GISMO_ASSERT(3 == dim, "Implemented for dimension 2 and 3 only.");
    for (index_t r = 0; r != 3; ++r) {
      // row r is the cross product of the two following rows of J
      const index_t r1 = (r + 1) % 3, r2 = (r + 2) % 3;
      cof(r, 0) = J(r1, 1) * J(r2, 2) - J(r1, 2) * J(r2, 1);
      cof(r, 1) = J(r1, 2) * J(r2, 0) - J(r1, 0) * J(r2, 2);
      cof(r, 2) = J(r1, 0) * J(r2, 1) - J(r1, 1) * J(r2, 0);
    }
  }
}

/// Derivative of the cofactor matrix of J in the direction of a perturbation
/// of row \a b of J by the vector \a g (i.e. of the coordinate \a b of the
/// geometry by a function with parametric gradient \a g).
template<class MatJ, class Vec, class MatC>
inline void cofactorDerivative(const MatJ &J, const index_t b,
                               const Vec &g, MatC &dCof) {
  const index_t dim = J.rows();
  dCof.setZero(dim, dim);
  if (2 == dim) {
    if (0 == b) {
      dCof(1, 0) = -g(1);
      dCof(1, 1) = g(0);
    } else {
      dCof(0, 0) = g(1);
      dCof(0, 1) = -g(0);
    }
  } else {
    const index_t r1 = (b + 1) % 3, r2 = (b + 2) % 3;
    // row r1 = J_{r2} x J_b  ->  J_{r2} x g
    dCof(r1, 0) = J(r2, 1) * g(2) - J(r2, 2) * g(1);
    dCof(r1, 1) = J(r2, 2) * g(0) - J(r2, 0) * g(2);
    dCof(r1, 2) = J(r2, 0) * g(1) - J(r2, 1) * g(0);
    // row r2 = J_b x J_{r1}  ->  g x J_{r1}
    dCof(r2, 0) = g(1) * J(r1, 2) - g(2) * J(r1, 1);
    dCof(r2, 1) = g(2) * J(r1, 0) - g(0) * J(r1, 2);
    dCof(r2, 2) = g(0) * J(r1, 1) - g(1) * J(r1, 0);
  }
}

/// Integrand of a regularized Winslow-type energy at one point, with its
/// gradient and Hessian with respect to the Jacobian matrix J of the map
/// (J vectorized row by row: index c * d + a for d x_c / d xi_a):
///   form 0: f = |J|^2 / chi(D, eps)^(2/d)   (Garanzha's conformal form)
///   form 1: f = |cof J|^2 / chi(D, eps)      (inverse Dirichlet energy:
///           |grad_x xi|^2 det J; in 2D the same as form 0)
/// with D = orientation * det J. Hessian (if \a hessian is given):
///   type 0: Gauss-Newton (f = |r|^2, r = vec(J or cof J) / chi^(p/2)),
///   type 1: exact, type 2: exact, projected onto the positive semidefinite
///   matrices (negative eigenvalues set to zero).
template<short_t d, typename T>
T garanzhaEnergyDerivatives(const gsEigen::Matrix<T, d, d> &J,
                            const T orientation, const T eps,
                            const index_t form, const index_t hessianType,
                            gsEigen::Matrix<T, d * d, 1> *gradient,
                            gsEigen::Matrix<T, d * d, d * d> *hessian) {
  typedef gsEigen::Matrix<T, d, d> MatD;
  typedef gsEigen::Matrix<T, d * d, 1> VecJ;
  typedef gsEigen::Matrix<T, d * d, d * d> MatJ;
  MatD cof, dCof;
  cofactorMatrix(J, cof);
  const T D = orientation * J.row(0).dot(cof.row(0));
  const T root = math::sqrt(D * D + eps * eps);
  const T chi = garanzhaChi(D, eps);
  const T rho = 1 / root;                       // chi' / chi
  const bool inverse = (gsBarrierCore<d, T>::EnergyForm::InverseDirichlet == form && 3 == d);
  const T p = inverse ? T(1) : T(2) / d;
  const T chiP = inverse ? chi : (2 == d ? chi : math::pow(chi, p));
  // E = d vec(cof) / d vec(J), column b * d + e
  MatJ dCofdJ;
  gsEigen::Matrix<T, d, 1> unit;
  for (short_t b = 0; b != d; ++b)
    for (short_t e = 0; e != d; ++e) {
      unit.setZero();
      unit(e) = 1;
      cofactorDerivative(J, b, unit, dCof);
      for (short_t r = 0; r != d; ++r)
        for (short_t c = 0; c != d; ++c) { dCofdJ(r * d + c, b * d + e) = dCof(r, c); }
    }
  VecJ Gvec, gradD;
  for (short_t r = 0; r != d; ++r)
    for (short_t c = 0; c != d; ++c) {
      Gvec(r * d + c) = inverse ? cof(r, c) : J(r, c);
      gradD(r * d + c) = orientation * cof(r, c);
    }
  const T g = Gvec.squaredNorm();
  const T f = g / chiP;
  if (!gradient && !hessian) { return f; }
  // E = d Gvec / d vec(J)
  const VecJ gradG = inverse ? VecJ(2 * dCofdJ.transpose() * Gvec) : VecJ(2 * Gvec);
  if (gradient) { *gradient = (gradG - p * g * rho * gradD) / chiP; }
  if (!hessian) { return f; }
  MatJ EtE;
  if (inverse) { EtE.noalias() = dCofdJ.transpose() * dCofdJ; }
  else { EtE.setIdentity(); }
  const MatJ mixed = gradG * gradD.transpose() + gradD * gradG.transpose();
  if (gsBarrierCore<d, T>::EnergyHessian::GaussNewton == hessianType) {
    *hessian = (2 * EtE - (p * rho / 2) * mixed
                + (p * p * rho * rho * g / 2) * (gradD * gradD.transpose())) / chiP;
    return f;
  }
  // exact: H_g - p rho (grad g grad D^T + sym) + g (p(p+1) rho^2 - p chi''/chi)
  //        grad D grad D^T - p g rho H_D,   chi''/chi = (root - D) / root^3
  MatJ Hg = 2 * EtE;
  if (inverse) {
    // second derivatives of cof (bilinear in J): d^2 cof / dJ_ca dJ_be
    MatD unitJ;
    for (short_t c = 0; c != d; ++c)
      for (short_t a = 0; a != d; ++a) {
        unitJ.setZero();
        unitJ(c, a) = 1;
        for (short_t b = 0; b != d; ++b)
          for (short_t e = 0; e != d; ++e) {
            unit.setZero();
            unit(e) = 1;
            cofactorDerivative(unitJ, b, unit, dCof);
            T sum = 0;
            for (short_t r = 0; r != d; ++r)
              for (short_t q = 0; q != d; ++q) { sum += cof(r, q) * dCof(r, q); }
            Hg(c * d + a, b * d + e) += 2 * sum;
          }
      }
  }
  const T chi2 = (root - D) / (root * root * root);
  *hessian = (Hg - p * rho * mixed
              + g * (p * (p + 1) * rho * rho - p * chi2) * (gradD * gradD.transpose())
              - (p * g * rho * orientation) * dCofdJ) / chiP;
  if (gsBarrierCore<d, T>::EnergyHessian::ProjectedNewton == hessianType) {
    gsEigen::SelfAdjointEigenSolver<MatJ> eig(*hessian);
    const VecJ lambda = eig.eigenvalues().cwiseMax(T(0));
    *hessian = eig.eigenvectors() * lambda.asDiagonal() * eig.eigenvectors().transpose();
  }
  return f;
}

template<short_t d, typename T> class gsRegularizedPDEAssembler;

/// Values and derivatives up to order \a n of the p+1 B-splines of degree
/// \a p that are nonzero on the knot span [U[span], U[span+1]), at \a u
/// (Piegl & Tiller, The NURBS Book, algorithm A2.3). Unlike an evaluation
/// that locates the span of \a u, this gives the functions of the given span
/// also at its end points (one-sided limits from inside the span).
/// ders(k, j): k-th derivative of the j-th function of the span.
template<class T>
void bsplineSpanDerivatives(const gsKnotVector<T> &U, const index_t p,
                            const index_t span, const T u, const index_t n,
                            gsMatrix<T> &ders) {
  ders.setZero(n + 1, p + 1);
  gsMatrix<T> ndu(p + 1, p + 1), a(2, p + 1);
  std::vector<T> left(p + 1), right(p + 1);
  ndu(0, 0) = 1;
  for (index_t j = 1; j <= p; ++j) {
    left[j] = u - U[span + 1 - j];
    right[j] = U[span + j] - u;
    T saved = 0;
    for (index_t r = 0; r < j; ++r) {
      ndu(j, r) = right[r + 1] + left[j - r];
      const T temp = ndu(r, j - 1) / ndu(j, r);
      ndu(r, j) = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu(j, j) = saved;
  }
  for (index_t j = 0; j <= p; ++j) { ders(0, j) = ndu(j, p); }
  const index_t nn = std::min(n, p); // higher derivatives vanish
  for (index_t r = 0; r <= p; ++r) {
    index_t s1 = 0, s2 = 1;
    a(0, 0) = 1;
    for (index_t k = 1; k <= nn; ++k) {
      T d = 0;
      const index_t rk = r - k, pk = p - k;
      if (r >= k) {
        a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
        d = a(s2, 0) * ndu(rk, pk);
      }
      const index_t j1 = (rk >= -1) ? 1 : -rk;
      const index_t j2 = (r - 1 <= pk) ? k - 1 : p - r;
      for (index_t j = j1; j <= j2; ++j) {
        a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
        d += a(s2, j) * ndu(rk + j, pk);
      }
      if (r <= pk) {
        a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
        d += a(s2, k) * ndu(r, pk);
      }
      ders(k, r) = d;
      std::swap(s1, s2);
    }
  }
  T factor = p;
  for (index_t k = 1; k <= nn; ++k) {
    ders.row(k) *= factor;
    factor *= (p - k);
  }
}

/// Coefficients w such that A : H = w . h for a symmetric matrix A and the
/// second derivatives h stored as (11, 22, [33,] 12, [13, 23]).
template<class MatA, class VecW>
inline void packSymmetric(const MatA &A, VecW &w) {
  const index_t dim = A.rows();
  w.resize(dim * (dim + 1) / 2);
  index_t m = dim;
  for (index_t a = 0; a != dim; ++a) {
    w(a) = A(a, a);
    for (index_t b = a + 1; b != dim; ++b) { w(m++) = 2 * A(a, b); }
  }
}

/// Relative tolerance of the tests whether det(J) depends on the free dofs
/// (detIsControllable and the rank test of gsRegularizedPDEAssembler)
const real_t detDependenceTolerance = 1e-10;

/// Tells whether det(J) at a point can be changed by the free control points:
/// \a cof is the cofactor matrix of J, \a derivBasis (dim x N) holds the
/// parametric gradients of the N active basis functions with indices
/// \a actives on patch \a patch. The derivative of det(J) with respect to
/// the coordinate c of the control point i is  cof.row(c) . grad(phi_i).
/// It vanishes for all free control points, e.g., on the edges and corners of
/// the parameter domain, and where the parameterization of a boundary face
/// is degenerate: there det(J) is determined by the boundary alone.
template<class MatC, class Mat>
inline bool detIsControllable(const gsDofMapper &mapper,
                              const gsMatrix<index_t> &actives,
                              const index_t patch,
                              const MatC &cof,
                              const Mat &derivBasis) {
  typedef typename MatC::Scalar T;
  const index_t dim = cof.rows();
  const T tol = T(detDependenceTolerance) * cof.norm()
              * derivBasis.colwise().norm().maxCoeff();
  for (index_t i = 0; i != derivBasis.cols(); ++i) {
    for (index_t c = 0; c != dim; ++c) {
      if (mapper.is_free(actives(i, 0), patch, c)
          && math::abs(cof.row(c).dot(derivBasis.col(i))) > tol) {
        return true;
      }
    }
  }
  return false;
}

/// One-line summary of an Anderson acceleration solve (for Verbose > 0)
template<class Solver>
inline std::string aaSummary(const std::string &stage, const Solver &solver,
                             const double seconds) {
  std::ostringstream os;
  os << "[PDE] " << stage << ": " << solver.iterations() << " iterations, "
     << (solver.converged() ? "converged" : "NOT converged")
     << " (residual " << solver.residualNorm() << "), "
     << solver.preconditionerUpdates() << " preconditioner updates, "
     << solver.restarts() << " restarts, "
     << solver.innerIterations() << " GMRES iterations, "
     << solver.innerFallbacks() << " LU fall-backs, "
     << solver.rejections() << " safeguard rejections; time " << seconds
     << " s (preconditioner " << solver.preconditionerTime()
     << " s, residual " << solver.residualTime()
     << " s, solves " << solver.solveTime() << " s)";
  return os.str();
}

template<short_t d, typename T>
gsOptionList gsBarrierCore<d, T>::defaultOptions() {
  gsOptionList options;

  // Print output 0: no print, 1: summary, 2: iterations and summary
  options.addInt("Verbose", "Print output level", 0);

  // Initialization Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch
  options.addInt("InitialMethod", "Initialization method", 0);

  // Parameterization Method: 0 Barrier patch (default), 1 Penalty patch, 2: PDE patch
  options.addInt("ParamMethod", "Parameterization method", 0);

  // Parameter and stopping criteria for foldover elimination step
  options.addReal("ff_Delta", "Delta for foldover-free optimization", 0.05);
  options.addInt("ff_MaxIterations",
                 "Max iterations for quality improvement",
                 1e4);
  options.addReal("ff_MinGradientLength",
                  "Min gradient length for foldover-free optimization",
                  1e-20);
  options.addReal("ff_MinStepLength",
                  "Min step length for foldover-free optimization",
                  1e-20);

  // Parameters and stopping criteria for quality improvement step
  options.addInt("qi_MaxIterations",
                 "Max iterations for quality improvement",
                 1e4);
  options.addReal("qi_MinGradientLength",
                  "Min gradient length for quality improvement",
                  1e-4);
  options.addReal("qi_MinStepLength",
                  "Min step length for quality improvement",
                  1e-4);

  // Set quadrature rules for objective function and gradient evaluation
  options.addReal("quA",
                  "Quadrature points: quA*deg + quB; For patchRule: Order of target space",
                  1.0);
  options.addInt("quB",
                 "Quadrature points: quA*deg + quB; For patchRule: Regularity of target space",
                 1);
  options.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
  options.addSwitch("overInt", "Apply over-integration?", false);

  // Preconditioner type for AA: 0: NO preconditioning, 1: Full Jacobian preconditioner, 2: Diagonal Jacobian preconditioner, 3: Diagonal Block Jacobian preconditioner
  options.addInt("AAPreconditionType", "Preconditioner type for AA", 0);

  // Update the preconditioner every N_update steps
  options.addInt("N_update", "update preconditioner every N_update steps", 10);

  // window size m for our preconditioned AA solver
  options.addInt("AAwindowsize", "window size for preconditioned AA solver", 5);

  // need improve the parameterization quality for PDE-based method?
  options.addSwitch("needPDEH1", "improve quality by H1 discrezation?", true);

  // options of the (regularized) PDE-based method: group "PDE"
  options.update(defaultPDEOptions().wrapIntoGroup("PDE"),
                 gsOptionList::addIfUnknown);
  return options;
}

template<short_t d, typename T>
gsOptionList gsBarrierCore<d, T>::defaultPDEOptions() {
  // See the page as_parameterization_example for the background of the
  // options and of their default values
  gsOptionList options;
  options.addSwitch("Garanzha", "Garanzha's regularization: det(J) replaced by "
                    "chi(det J, eps) = (det J + sqrt(det J^2 + eps^2)) / 2, "
                    "folded results untangled (off: original method)", true);
  options.addReal("GaranzhaEps", "Smallest regularization parameter eps, "
                  "relative to the mean of |det(J)|", 1e-4);
  options.addInt("LinearSolver", "Linear solvers of the regularized method "
                 "(LinearSolver): 0 automatic (by DirectSolverLimit), 1 direct, "
                 "2 iterative (GMRES / CG)", LinearSolver::Automatic);
  options.addInt("DirectSolverLimit", "Automatic choice of the linear solvers: "
                 "max. number of unknowns for direct solvers", 4000);
  options.addInt("InnerMaxIterations", "Max. number of GMRES iterations per "
                 "solve before falling back to sparse LU", 100);
  options.addInt("Preconditioner", "Preconditioner of GMRES (Preconditioner): "
                 "0 automatic, 1 block Gauss-Seidel by coordinate (sparse LU "
                 "blocks), 2 nodal algebraic multigrid, 3 block Gauss-Seidel by "
                 "coordinate (multigrid blocks)", Preconditioner::Automatic);
  options.addSwitch("FreeInterfaceEnergy", "Free interfaces: chord-scaled "
                    "parameter domains and energy minimization instead of the "
                    "H2 and H1 stages", true);
  options.addSwitch("EnergyFinal", "The energy minimization is the final stage "
                    "(no H1 stage)", false);
  options.addSwitch("EnergyFallback", "Minimize the energy to convergence if "
                    "the H1 stage does not converge", true);
  options.addReal("EnergyTol", "Relative energy decrease at which the final "
                  "energy minimization stops", 1e-11);
  options.addInt("EnergyForm", "Energy of the untangling and of the final "
                 "energy stage (EnergyForm): 0 |J|^2 / chi^(2/d), "
                 "1 |cof J|^2 / chi (the same in 2D)", EnergyForm::Conformal);
  options.addInt("EnergyHessian", "Hessian of the energy minimization "
                 "(EnergyHessian): -1 automatic, 0 Gauss-Newton, 1 exact, "
                 "2 projected Newton", EnergyHessian::Automatic);
  options.addReal("H1Safeguard", "Regularized H1 stage: restart from the best "
                  "iterate when the residual norm exceeds this factor times the "
                  "smallest one (0: off)", 100);
  options.addSwitch("SkipH2", "Experimental: no H2 stage; untangling and H1 "
                    "start from the given interior", false);
  options.addInt("H1Form", "Experimental: regularized H1 stage (H1Form): "
                 "0 original form, 1 weak form of Laplace's equation for the "
                 "inverse map", H1Form::Original);
  options.addInt("ExtraNodes", "Experimental: additional Gauss-Lobatto nodes "
                 "per direction of the regularized stage", 0);
  options.addInt("AdaptiveRounds", "Experimental: max. number of adaptive "
                 "refinements of the quadrature where det(J) <= 0 between the "
                 "nodes (0: off)", 0);
  options.addInt("CheckExtraNodes", "Adaptive refinement (AdaptiveRounds): "
                 "additional Gauss-Lobatto nodes per direction of the check", 4);
  return options;
}

template<short_t d, typename T>
gsOptionList gsBarrierCore<d, T>::withPDEDefaults(const gsOptionList &options) {
  gsOptionList result = defaultPDEOptions().wrapIntoGroup("PDE");
  result.update(options, gsOptionList::addIfUnknown);
  return result;
}

/// Compute the area of the computational domain
template<short_t d, typename T>
T gsBarrierCore<d, T>::computeArea(const gsMultiPatch<T> &mp) {
  return computeAreaInterior(mp);
}

/// Compute the area of a multi-patch representing computational domain
template<short_t d, typename T>
T gsBarrierCore<d, T>::computeAreaInterior(const gsMultiPatch<T> &multiPatch) {

  // Initializing an expression evaluator
  gsExprEvaluator<T> evaluator;

  // Setting integration elements
  evaluator.setIntegrationDomain(multiPatch.domain());

  // Getting the geometry map
  geometryMap geometry = evaluator.getMap(multiPatch);

  // Computing the area by integrating over the determinant of the Jacobian of the geometry map
  T area = evaluator.integral(jac(geometry).det());

  return area;
}

template<short_t d, typename T>
T gsBarrierCore<d, T>::computeAreaBoundary(const gsMultiPatch<T> &)
{
    GISMO_ERROR("NO_IMPLEMENTATION");
}

enum class ParamMethod {
  PDEPatch,
  PenaltyPatch,
  PenaltyPatch2,
  BarrierPatch,
  VariationalHarmonicPatch,
  // Add new methods here...
};

template<short_t d, typename T>
gsMultiPatch<T> gsBarrierCore<d, T>::compute(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             const gsOptionList &options) {
  gsMultiPatch<T> result;

  ParamMethod
      method = static_cast<ParamMethod>(options.askInt("ParamMethod", 0));

  switch (method) {
    case ParamMethod::PenaltyPatch: {
#ifdef gsHLBFGS_ENABLED
      result = computePenaltyPatch(mp, mapper, options);
#else
      GISMO_ERROR("PenaltyPatch not available without gsHLBFGS. Please "
                  "compile with gsHLBGFS enabled");
#endif
      break;
    }
    case ParamMethod::PenaltyPatch2: {
#ifdef gsHLBFGS_ENABLED
      result = computePenaltyPatch2(mp, mapper, options);
#else
      GISMO_ERROR("PenaltyPatch2 not available without gsHLBFGS. Please "
                  "compile with gsHLBGFS enabled");
#endif
      break;
    }
    case ParamMethod::BarrierPatch: {
#ifdef gsHLBFGS_ENABLED
      result = computeBarrierPatch(mp, mapper, options);
#else
      GISMO_ERROR("BarrierPatch not available without gsHLBFGS. Please "
                  "compile with gsHLBGFS enabled");
#endif
      break;
    }
    case ParamMethod::VariationalHarmonicPatch: {
#ifdef gsHLBFGS_ENABLED
      result = computeVHPatch(mp, mapper, options);
#else
      GISMO_ERROR("VHPatch not available without gsHLBFGS. Please "
                  "compile with gsHLBGFS enabled");
#endif
      break;
    }
    case ParamMethod::PDEPatch:
    default: {
      if (method != ParamMethod::PDEPatch) {
        gsWarn << "Invalid ParamMethod value. Defaulting to PDEPatch.\n";
      }
      result = computePDEPatch(mp, mapper, options);
      break;
    }
  }
  return result;
}

#ifdef gsHLBFGS_ENABLED
// Modified Variational Harmonic Method
template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computeVHPatch(const gsMultiPatch<T> &mp,
                                    const gismo::gsDofMapper &mapper,
                                    const gismo::gsOptionList &options) {

  verboseLog("Start variational harmonic parameterization construction...",
             options.askInt("Verbose", 0));

  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);

  gsObjVHPt<d, T> objVHPt(mp, mapper);
  objVHPt.applyOptions(options);

  gsHLBFGS<T> optimizer(&objVHPt);
  setOptimizerOptions<T>(optimizer, options);

  optimizer.solve(initialGuessVector);
  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(optimizer.currentDesign(), mapper, result);

  verboseLog("Finished!", options.askInt("Verbose", 0));

  return result;
}

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePenaltyPatch(const gsMultiPatch<T> &mp,
                                         const gsDofMapper &mapper,
                                         const gsOptionList &options) {

  verboseLog("Start penalty function-based parameterization construction...",
             options.askInt("Verbose", 0));

  // Compute scaledArea and initial guess vector
  T scaledArea = computeAreaInterior(mp);

  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);

  gsObjPenaltyPt<d, T> objPenaltyPt(mp, mapper);
  gsOptionList thisOptions = options;
  thisOptions.addReal("qi_lambda1", "Sets the lambda_1 value for Emips", 1.0);
  thisOptions.addReal("qi_lambda2",
                      "Sets the lambda 2 value for Eunif",
                      1.0 / pow(scaledArea, 2));
  objPenaltyPt.applyOptions(thisOptions);

  gsHLBFGS<T> optimizer(&objPenaltyPt);
  setOptimizerOptions<T>(optimizer, options);

  optimizer.solve(initialGuessVector);
  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(optimizer.currentDesign(), mapper, result);

  verboseLog("Finished!", options.askInt("Verbose", 0));

  return result;
}

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePenaltyPatch2(const gsMultiPatch<T> &mp,
                                          const gsDofMapper &mapper,
                                          const gsOptionList &options) {

  verboseLog("Penalty function-based (2) parameterization construction ...",
             options.askInt("Verbose", 0));

  // Compute scaledArea and initial guess vector
  T scaledArea = computeAreaInterior(mp);

  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);

  gsObjPenaltyPt2<d, T> objPenaltyPt(mp, mapper);
  objPenaltyPt.setEps(options.askReal("penaltyFactor", 0.01) * scaledArea);

  gsOptionList thisOptions = options;
  thisOptions.addReal("qi_lambda1", "Sets the lambda_1 value for Emips", 1.0);
  thisOptions.addReal("qi_lambda2",
                      "Sets the lambda 2 value for Eunif",
                      1.0 / pow(scaledArea, 2));
  objPenaltyPt.applyOptions(thisOptions);

  gsHLBFGS<T> optimizer(&objPenaltyPt);
  setOptimizerOptions<T>(optimizer, options);

  optimizer.solve(initialGuessVector);
  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(optimizer.currentDesign(), mapper, result);

  verboseLog("Finished!", options.askInt("Verbose", 0));

  return result;
}

template<short_t d, typename T>
gsMultiPatch<T>
    gsBarrierCore<d,T>::computeBarrierPatch(const gsMultiPatch<T> &mp,
                                            const gsDofMapper &mapper,
                                            const gsOptionList &options) {

  // Compute scaledArea and initial guess vector
  T scaledArea = computeAreaInterior(mp);

  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);

  // STEP 2: foldover elimination step
  foldoverElimination(mp, mapper, initialGuessVector, scaledArea, options);

  // STEP 3: parameterization quality improvement
  qualityImprovement(mp, mapper, initialGuessVector, scaledArea, options);

  // Update the result with optimized parameterization
  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch(initialGuessVector, mapper, result);

  return result;
}

template<short_t d, typename T>
void gsBarrierCore<d, T>::foldoverElimination(const gsMultiPatch<T> &mp,
                                              const gsDofMapper &mapper,
                                              gsVector<T> &initialGuessVector,
                                              const T &scaledArea,
                                              const gsOptionList &options) {

  verboseLog("Start foldover elimination step...",
             options.askInt("Verbose", 0));
  constexpr T EPSILON = 1e-20;
  constexpr int MAX_ITER = 10;
  gsObjFoldoverFree<d, T> objFoldoverFree(mp, mapper);
  objFoldoverFree.addOptions(options);

  gsHLBFGS<T> optFoldoverFree(&objFoldoverFree);
  optFoldoverFree.options().setInt("MaxIterations",
                                   options.askInt("ff_MaxIterations", 1e4));
  optFoldoverFree.options().setReal("MinGradLen",
                                    options.askReal("ff_MinGradientLength",
                                                    1e-12));
  optFoldoverFree.options().setReal("MinStepLen",
                                    options.askReal("ff_MinStepLength", 1e-12));
  optFoldoverFree.options().setInt("Verbose", options.askInt("Verbose", 0));

  T Efoldover = std::numeric_limits<T>::max();
  for (index_t it = 0; it < MAX_ITER; ++it) {
    T delta = pow(0.1, it) * 5e-2 * scaledArea; // parameter delta
    objFoldoverFree.setDelta(delta);

    optFoldoverFree.solve(initialGuessVector);

    Efoldover = optFoldoverFree.objective();
    initialGuessVector = optFoldoverFree.currentDesign();

    if (Efoldover <= EPSILON) { break; }
  }

  GISMO_ENSURE(Efoldover <= EPSILON,
               "Foldover elimination failed. The foldover-energy value is " +
                   std::to_string(Efoldover) +
                   ". This suggests there may be issues with the input data.");
}

template<short_t d, typename T>
void gsBarrierCore<d, T>::qualityImprovement(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             gsVector<T> &initialGuessVector,
                                             const T &scaledArea,
                                             const gsOptionList &options) {
  verboseLog("Start parameterization quality improvement step...",
             options.askInt("Verbose", 0));
  gsObjQualityImprovePt<d, T> objQualityImprovePt(mp, mapper);

  gsOptionList thisOptions = options;
  thisOptions.addReal("qi_lambda1", "Sets the lambda_1 value for Emips", 1.0);
  thisOptions.addReal("qi_lambda2",
                      "Sets the lambda 2 value for Eunif",
                      1.0 / pow(scaledArea, 2));
  objQualityImprovePt.applyOptions(thisOptions);

  gsHLBFGS<T> optQualityImprovePt(&objQualityImprovePt);
  setOptimizerOptions<T>(optQualityImprovePt, options);

  optQualityImprovePt.solve(initialGuessVector);
  initialGuessVector = optQualityImprovePt.currentDesign();
}

template<short_t d, typename T>
gsObjFoldoverFree<d, T>::gsObjFoldoverFree(const gsMultiPatch<T> &patches,
                                           gsDofMapper mapper)
    :
    m_mp(patches),
    m_mapper(std::move(mapper)),
    m_mb(m_mp) {
  defaultOptions();
  m_assembler.setIntegrationDomain(m_mb.domain());
  m_evaluator = gsExprEvaluator<T>(m_assembler);
}

template<short_t d, typename T>
void gsObjFoldoverFree<d, T>::defaultOptions() {
  // @Ye, make this reasonable default options
  m_options.addReal("ff_Delta", "Sets the delta value", 1e-2);
}

template<short_t d, typename T>
void gsObjFoldoverFree<d, T>::addOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
void gsObjFoldoverFree<d, T>::applyOptions(const gsOptionList &options) {
  m_eps = m_options.getReal("ff_Delta");
  m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
T gsObjFoldoverFree<d, T>::evalObj(const gsAsConstVector<T> &u) const {
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);
  geometryMap G = m_evaluator.getMap(m_mp);

  auto EfoldoverFree = (m_eps - jac(G).det()).ppartval();
  return m_evaluator.integral(EfoldoverFree);
}

template<short_t d, typename T>
void gsObjFoldoverFree<d, T>::gradObj_into(const gsAsConstVector<T> &u,
                                           gsAsVector<T> &result) const {
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);
  geometryMap G = m_assembler.getMap(m_mp);

  // Only call these once if their results don't change.
  space space1 = m_assembler.getSpace(m_mb, d);
  space1.setupMapper(m_mapper);

  // |J|' w.r.t. physical coordinates x and y
  auto derJacDet = frprod2(space1, jac(G).tr().adj());

  gsConstantFunction<T> zeroFunc(gsVector<T>::Zero(d), d);
  auto zeroVar = m_evaluator.getVariable(zeroFunc);

  auto Eder = ternary(m_eps - jac(G).det(), -derJacDet, space1 * zeroVar);

  m_assembler.initSystem();
  m_assembler.assemble(Eder);

  result.resize(m_assembler.rhs().rows());
  std::copy(m_assembler.rhs().data(),
            m_assembler.rhs().data() + m_assembler.rhs().rows(),
            result.data());
}

template<short_t d, typename T>
gsObjQualityImprovePt<d, T>::gsObjQualityImprovePt(
    const gsMultiPatch<T> &patches,
    gsDofMapper mapper)
    :
    m_mp(patches),
    m_mapper(std::move(mapper)),
    m_mb(m_mp) {
    m_assembler.setIntegrationDomain(m_mb.domain());
  m_evaluator = gsExprEvaluator<T>(m_assembler);
}

template<short_t d, typename T>
void gsObjQualityImprovePt<d, T>::defaultOptions() {
  // @Ye, make this reasonable default options
  m_options.addReal("qi_lambda1", "Sets the lambda 1 value", 1.0);
  m_options.addReal("qi_lambda2", "Sets the lambda 2 value", 1.0);
}

template<short_t d, typename T>
void gsObjQualityImprovePt<d, T>::addOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
void gsObjQualityImprovePt<d, T>::applyOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
  m_lambda1 = m_options.getReal("qi_lambda1");
  m_lambda2 = m_options.getReal("qi_lambda2");
  m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
T gsObjQualityImprovePt<d, T>::evalObj(const gsAsConstVector<T> &u) const {
  return evalObj_impl<d>(u);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 2, T>::type
gsObjQualityImprovePt<d, T>::evalObj_impl(const gsAsConstVector<T> &u) const {
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);
  geometryMap G = m_evaluator.getMap(m_mp);

  if (m_evaluator.min(jac(G).det()) < 0) {
    return std::numeric_limits<T>::max();
  } else {
    auto Ewinslow = jac(G).sqNorm() / jac(G).det();
    auto Euniform = pow(jac(G).det(), 2);

    return m_evaluator.integral(m_lambda1 * Ewinslow + m_lambda2 * Euniform);
  }
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 3, T>::type
gsObjQualityImprovePt<d, T>::evalObj_impl(const gsAsConstVector<T> &u) const {
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  geometryMap G = m_evaluator.getMap(m_mp);

  if (m_evaluator.min(jac(G).det()) < 0) {
    return std::numeric_limits<T>::max();
  } else {
    auto Euniform = pow(jac(G).det(), 2);
    auto Ewinslow = 0.5 * (jac(G).sqNorm() * jac(G).inv().sqNorm());

    return m_evaluator.integral(m_lambda1 * Ewinslow + m_lambda2 * Euniform);
  }
}

template<short_t d, typename T>
void gsObjQualityImprovePt<d, T>::gradObj_into(const gsAsConstVector<T> &u,
                                               gsAsVector<T> &result) const {
  gradObj_into_impl<d>(u, result);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 2, T>::type
gsObjQualityImprovePt<d, T>::gradObj_into_impl(const gsAsConstVector<T> &u,
                                               gsAsVector<T> &result) const {
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  geometryMap G = m_assembler.getMap(m_mp);

  space space1 = m_assembler.getSpace(m_mb, d); // 1D space!!
  space1.setupMapper(m_mapper);

  // |J|' w.r.t. physical coordinates x and y
  auto derJacDet = frprod2(space1, jac(G).tr().adj());

  auto Ewinslow = jac(G).sqNorm() / jac(G).det();
  auto derEwinslow = 2.0 / jac(G).det() * (frprod2(space1, jac(G))) -
      Ewinslow.val() / jac(G).det() * derJacDet;
  auto derEuniform = 2.0 * jac(G).det() * derJacDet;

  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());
  return EXIT_SUCCESS;
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 3, T>::type
gsObjQualityImprovePt<d, T>::gradObj_into_impl(const gsAsConstVector<T> &u,
                                               gsAsVector<T> &result) const {
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);
  geometryMap G = m_assembler.getMap(m_mp);

  space space1 = m_assembler.getSpace(m_mb, d); // 1D space!!
  space1.setupMapper(m_mapper);

  //      |J|' w.r.t. physical coordinates x and y
  auto derJacDet = frprod2(space1, jac(G).tr().adj());

  auto derEwinslow = frprod2(space1, (jac(G).inv().sqNorm() * jac(G) - jac(G)
      .sqNorm() * (jac(G).tr() * jac(G) * jac(G).tr()).inv()));

// (gradient of the other choice, |J|^2 / det(J)^(2/3):
//  2 frprod2(u, J) / det(J)^(2/3) - (2/3) |J|^2 / det(J)^(5/3) d det(J))
  auto derEuniform = 2.0 * jac(G).det() * derJacDet;

  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());
  return EXIT_SUCCESS;
}

template<short_t d, typename T>
gsObjVHPt<d, T>::gsObjVHPt(const gsMultiPatch<T> &patches,
                           gsDofMapper mapper)
    :
    m_mp(patches),
    m_mapper(std::move(mapper)),
    m_mb(m_mp) {
  defaultOptions();
  m_assembler.setIntegrationDomain(m_mb.domain());
  m_evaluator = gsExprEvaluator<T>(m_assembler);
}

template<short_t d, typename T>
void gsObjVHPt<d, T>::defaultOptions() {
  // @Ye, make this reasonable default options
  m_options.addReal("qi_lambda1", "Sets the lambda 1 value", 1.0);
  m_options.addReal("qi_lambda2", "Sets the lambda 2 value", 1.0);
}

template<short_t d, typename T>
void gsObjVHPt<d, T>::addOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
void gsObjVHPt<d, T>::applyOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
  m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
T gsObjVHPt<d, T>::evalObj(const gsAsConstVector<T> &u) const {
  return evalObj_impl<d>(u);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 2, T>::type
gsObjVHPt<d, T>::evalObj_impl(const gsAsConstVector<T> &u) const {
  // Convert the free vector to multi-patch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Get the map for the geometry
  geometryMap G = m_evaluator.getMap(m_mp);

  // Get and set up the 1D space
  space space1 = m_assembler.getSpace(m_mb, d);
  space1.setupMapper(m_mapper);

  // Compute the metric matrix
  auto metricMat = jac(G).tr() * jac(G);

  // Compute the Hessian
  auto hessG = hess(G);

  // Calculate the scale factor
  auto scale = metricMat.trace();

  // Calculate the Lx value
  auto Lx = hessG % metricMat.adj() / scale.val();

  // Build the nonlinear system
  auto nonlinearSystem = frprod3(space1, Lx);

  // Initialize and assemble the system
  m_assembler.initSystem();
  m_assembler.assemble(nonlinearSystem);

  // Create the result vector
  gsVector<T> result = gsAsVector<T>(
      const_cast<T *>(m_assembler.rhs().data()),
      m_assembler.rhs().rows());

  return result.squaredNorm();
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 3, T>::type
gsObjVHPt<d, T>::evalObj_impl(const gsAsConstVector<T> &u) const {
  GISMO_NO_IMPLEMENTATION;
}

template<short_t d, typename T>
void gsObjVHPt<d, T>::gradObj2_into(const gsAsConstVector<T> &u,
                                    gsAsVector<T> &result) const {
  gradObj_into_impl<d>(u, result);
}

template<short_t d, typename T>
gsObjPenaltyPt<d, T>::gsObjPenaltyPt(const gsMultiPatch<T> &patches,
                                     gsDofMapper mapper)
    :
    m_mp(patches),
    m_mapper(std::move(mapper)),
    m_mb(m_mp) {
  defaultOptions();
  m_assembler.setIntegrationDomain(m_mb.domain());
  m_evaluator = gsExprEvaluator<T>(m_assembler);
}

template<short_t d, typename T>
void gsObjPenaltyPt<d, T>::defaultOptions() {
  // @Ye, make this reasonable default options
  m_options.addReal("qi_lambda1", "Sets the lambda 1 value", 1.0);
  m_options.addReal("qi_lambda2", "Sets the lambda 2 value", 1.0);
}

template<short_t d, typename T>
void gsObjPenaltyPt<d, T>::addOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
void gsObjPenaltyPt<d, T>::applyOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
  m_lambda1 = m_options.getReal("qi_lambda1");
  m_lambda2 = m_options.getReal("qi_lambda2");
  m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
T gsObjPenaltyPt<d, T>::evalObj(const gsAsConstVector<T> &u) const {
  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Generate map for the geometry
  geometryMap G = m_evaluator.getMap(m_mp);

  // Generate epsilon value
  gsConstantFunction<T> eps1(pow(m_eps, 2.0), d);
  auto eps = m_evaluator.getVariable(eps1);

  // Calculate chi value - penalty function
  auto
      chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));

  // Calculate Ewinslow and Euniform
  auto Ewinslow = jac(G).sqNorm() / pow(chi, (2.0 / static_cast<T>(d)));
  auto Euniform = pow(jac(G).det(), 2.0);

  return m_evaluator.integral(m_lambda1 * Ewinslow + m_lambda2 * Euniform);
}

template<short_t d, typename T>
void gsObjPenaltyPt<d, T>::gradObj_into(const gsAsConstVector<T> &u,
                                        gsAsVector<T> &result) const {
  gradObj_into_impl<d>(u, result);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 2, T>::type
gsObjPenaltyPt<d, T>::gradObj_into_impl(const gsAsConstVector<T> &u,
                                        gsAsVector<T> &result) const {

  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Get map and space
  geometryMap G = m_assembler.getMap(m_mp);
  space space1 = m_assembler.getSpace(m_mb, d);
  space1.setupMapper(m_mapper);

  // Derivative of Jacobian determinant
  auto derJacDet = jac(space1) % jac(G).tr().adj();

  // Define constants for evaluation
  gsConstantFunction<T> eps1(pow(m_eps, 2.0), d);
  gsConstantFunction<T> unit1(1.0, d);
  auto eps = m_evaluator.getVariable(eps1);
  auto unit = m_evaluator.getVariable(unit1);

  // Compute common term, chi and chip
  auto commonTerm = pow(eps.val() + pow(jac(G).det(), 2.0), 0.5);
  auto chi = 0.5 * (jac(G).det() + commonTerm);
  auto chip = 0.5 * (unit.val() + jac(G).det() / commonTerm);

  // Compute Ewinslow
  auto Ewinslow = jac(G).sqNorm() / chi;

  // Compute derivatives of Ewinslow and Euniform
  auto derEwinslow = 1.0 / chi * (2.0 * frprod2(space1, jac(G)) -
      Ewinslow.val() * chip * derJacDet);
  auto derEuniform = 2.0 * jac(G).det() * derJacDet;

  // Assemble the system
  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);

  // Update the result
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());

  return EXIT_SUCCESS;
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 3, T>::type
gsObjPenaltyPt<d, T>::gradObj_into_impl(const gsAsConstVector<T> &u,
                                        gsAsVector<T> &result) const {
  const T twoThirds = 2.0 / 3.0;

  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Get map and space
  geometryMap G = m_assembler.getMap(m_mp);
  space space1 = m_assembler.getSpace(m_mb, d); // 1D space!!
  space1.setupMapper(m_mapper);

  // Derivative of Jacobian determinant
  auto derJacDet = frprod2(space1, jac(G).tr().adj());

  // Define constants for evaluation
  gsConstantFunction<T> eps1(pow(m_eps, 2.0), d);
  gsConstantFunction<T> unit1(1.0, d);
  auto eps = m_evaluator.getVariable(eps1);
  auto unit = m_evaluator.getVariable(unit1);

  // Compute common term, chi and chip
  auto commonTerm = pow(eps.val() + pow(jac(G).det(), 2.0), 0.5);
  auto chi = 0.5 * (jac(G).det() + commonTerm);
  auto chip = 0.5 * (unit.val() + jac(G).det() / commonTerm);

  // Compute Ewinslow
  auto Ewinslow = jac(G).sqNorm() / pow(chi, twoThirds);

  // Compute derivatives of Ewinslow and Euniform
  auto derEwinslow = 2.0 * frprod2(space1, jac(G)) / pow(chi, twoThirds) -
      twoThirds * Ewinslow / chi * chip * derJacDet;
  auto derEuniform = 2.0 * jac(G).det() * derJacDet;

  // Assemble the system
  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);

  // Update the result
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());

  return EXIT_SUCCESS;
}

// gsObjPenaltyPt2: my penalty function
template<short_t d, typename T>
gsObjPenaltyPt2<d, T>::gsObjPenaltyPt2(const gsMultiPatch<T> &patches,
                                       gsDofMapper mapper)
    :
    m_mp(patches),
    m_mapper(std::move(mapper)),
    m_mb(m_mp) {
  defaultOptions();
  m_assembler.setIntegrationDomain(m_mb.domain());
  m_evaluator = gsExprEvaluator<T>(m_assembler);
}

template<short_t d, typename T>
void gsObjPenaltyPt2<d, T>::defaultOptions() {
  // @Ye, make this reasonable default options
  m_options.addReal("qi_lambda1", "Sets the lambda 1 value", 1.0);
  m_options.addReal("qi_lambda2", "Sets the lambda 2 value", 1.0);
}

template<short_t d, typename T>
void gsObjPenaltyPt2<d, T>::addOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
void gsObjPenaltyPt2<d, T>::applyOptions(const gsOptionList &options) {
  m_options.update(options, gsOptionList::addIfUnknown);
  m_lambda1 = m_options.getReal("qi_lambda1");
  m_lambda2 = m_options.getReal("qi_lambda2");
  m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
}

template<short_t d, typename T>
T gsObjPenaltyPt2<d, T>::evalObj(const gsAsConstVector<T> &u) const {
  return evalObj_impl<d>(u);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 2, T>::type
gsObjPenaltyPt2<d, T>::evalObj_impl(const gsAsConstVector<T> &u) const {
  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Get the map G
  geometryMap G = m_evaluator.getMap(m_mp);

  // Constant function for evaluation
  gsConstantFunction<T> epsilon(m_eps, d);
  auto eps = m_evaluator.getVariable(epsilon);

  // Calculation of chi with ternary operation
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exp());
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());

  // Calculation of Ewinslow and Euniform
  auto Ewinslow = jac(G).sqNorm() / chi;
  auto Euniform = pow(jac(G).det(), 2.0);

  // Evaluate the integral and return
  return m_evaluator.integral(m_lambda1 * Ewinslow + m_lambda2 * Euniform);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 3, T>::type
gsObjPenaltyPt2<d, T>::evalObj_impl(const gsAsConstVector<T> &u) const {
  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Get the map G
  geometryMap G = m_evaluator.getMap(m_mp);

  // Define epsilon constant function
  gsConstantFunction<T> epsilon(m_eps, d);
  auto eps = m_evaluator.getVariable(epsilon);

  // Compute chi part and chi
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exp());
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());

  // Compute Ewinslow
  auto Ewinslow = 0.5 * (jac(G).sqNorm() * jac(G).inv().sqNorm()) *
      pow(jac(G).det(), 2.0) / pow(chi, 2.0);

  // Compute Euniform
  auto Euniform = pow(jac(G).det(), 2.0);

  // Compute and return the integral
  return m_evaluator.integral(m_lambda1 * Ewinslow + m_lambda2 * Euniform);
}

template<short_t d, typename T>
void gsObjPenaltyPt2<d, T>::gradObj_into(const gsAsConstVector<T> &u,
                                         gsAsVector<T> &result) const {
  gradObj_into_impl<d>(u, result);
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 2, T>::type
gsObjPenaltyPt2<d, T>::gradObj_into_impl(const gsAsConstVector<T> &u,
                                         gsAsVector<T> &result) const {
  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Get the map G
  geometryMap G = m_assembler.getMap(m_mp);

  // Set up the space mapper
  space spaceMapper = m_assembler.getSpace(m_mb, d);
  spaceMapper.setupMapper(m_mapper);

  // Compute the derivative of |J| with respect to physical coordinates x and y
  auto derJacDet = frprod2(spaceMapper, jac(G).tr().adj());

  // Define constant functions for epsilon and unity
  gsConstantFunction<T> epsilon(m_eps, d);
  gsConstantFunction<T> unity(1.0, d);

  // Get the variables for epsilon and unity
  auto eps = m_evaluator.getVariable(epsilon);
  auto unit = m_evaluator.getVariable(unity);

  // Define chi and its derivative chip
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exp());
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());
  auto chip = ternary(eps.val() - jac(G).det(), chiPPart.val(), unit.val());

  // Define Ewinslow and its derivative
  auto Ewinslow = jac(G).sqNorm() / chi;
  auto derEwinslow = 1.0 / chi * (2.0 * frprod2(spaceMapper, jac(G)) -
      Ewinslow.val() * chip * derJacDet);

  // Define the derivative of Euniform
  // (for the alternative term chi + 1 / chi: (chi' - chi' / chi^2) d det(J))
  auto derEuniform = 2.0 * jac(G).det() * derJacDet;

  // Assemble the system and set the result
  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());

  // Return success
  return EXIT_SUCCESS;
}

template<short_t d, typename T>
template<short_t _d>
typename std::enable_if<_d == 3, T>::type
gsObjPenaltyPt2<d, T>::gradObj_into_impl(const gsAsConstVector<T> &u,
                                         gsAsVector<T> &result) const {
  // Convert the free vector to multipatch
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  // Initialize the map G
  geometryMap G = m_assembler.getMap(m_mp);

  // Initialize and setup the space
  space space1 = m_assembler.getSpace(m_mb, d);
  space1.setupMapper(m_mapper);

  // Compute the derivative of |J| w.r.t. physical coordinates x and y
  auto derJacDet = frprod2(space1, jac(G).tr().adj());

  // Initialize constant functions
  gsConstantFunction<T> epsilonFunction(m_eps, d);
  gsConstantFunction<T> unityFunction(1.0, d);

  // Get the corresponding variables
  auto eps = m_evaluator.getVariable(epsilonFunction);
  auto unit = m_evaluator.getVariable(unityFunction);

  // Compute the chi part
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exp());

  // Ternary operation to compute chi and chip
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());
  auto chip = ternary(eps.val() - jac(G).det(), chiPPart.val(), unit.val());

  // Define additional computations
  auto jacFrobNorm2 = jac(G).sqNorm() * jac(G).inv().sqNorm();
  auto derJacFrob2 = frprod2(space1,
                             (jac(G).inv().sqNorm() * jac(G) - jac(G).sqNorm()
                                 * (jac(G).tr() * jac(G) * jac(G).tr()).inv()));

  auto derEwinslow = derJacFrob2 * pow(jac(G).det(), 2) / pow(chi, 2)
      + jacFrobNorm2 / pow(chi, 2)
          * (jac(G).det() - chip / chi * pow(jac(G).det(), 2)) * derJacDet;
  // (for the alternative term chi + 1 / chi: (chi' - chi' / chi^2) d det(J))
  auto derEuniform = 2.0 * jac(G).det() * derJacDet;

  // Initialize and assemble the system
  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);

  // Compute and return the result
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());
  return EXIT_SUCCESS;
}
#endif

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePDEPatch(const gsMultiPatch<T> &mp,
                                     const gismo::gsDofMapper &mapper,
                                     const gismo::gsOptionList &givenOptions) {
  // Nothing to do without free control points (e.g. a bilinear patch)
  if (mapper.freeSize() == 0) { return mp; }
  // all options of the group PDE are set from here on
  const gsOptionList options = withPDEDefaults(givenOptions);

  const bool regularized = options.getSwitch("PDE.Garanzha");
  // Free interfaces: chord-scaled parameter domains, energy minimization
  if (regularized && options.getSwitch("PDE.FreeInterfaceEnergy")
      && couplesPatches(mp, mapper)) {
    return computePDEPatchFreeInterfaces(mp, mapper, options);
  }

  // Experimental: no H2 stage, untangling and regularized H1 directly from
  // the given interior (e.g. a Coons patch)
  if (regularized && options.getSwitch("PDE.SkipH2")) {
    verboseLog("PDE-based parameterization construction (without H2) ...\n",
               options.askInt("Verbose", 0));
    gsVector<T> x = convertMultiPatchToFreeVector(mp, mapper);
    x = computeRegularizedH1(mp, mapper, x, options);
    gsMultiPatch<T> result = mp;
    convertFreeVectorToMultiPatch<T>(x, mapper, result);
    return result;
  }

  // Volumes, and planar domains with Garanzha's regularization: generic
  // implementation; otherwise the original planar one (PDE.Garanzha = 0)
  if (d != 2 || regularized) {
    return computePDEPatchGeneric(mp, mapper, options);
  }
  return computePDEPatchPlanar(mp, mapper, options);
}

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePDEPatchFreeInterfaces(const gsMultiPatch<T> &mp,
                                                   const gsDofMapper &mapper,
                                                   const gsOptionList &options) {
  verboseLog("PDE-based parameterization construction with free interfaces: "
             "chord-scaled parameter domains, energy minimization ...\n",
             options.askInt("Verbose", 0));
  gsMultiPatch<T> scaled = chordScaled(mp);
  gsOptionList energyOptions = options;
  energyOptions.setSwitch("PDE.EnergyFinal", true);
  gsVector<T> x = convertMultiPatchToFreeVector(scaled, mapper);
  x = computeRegularizedH1(scaled, mapper, x, energyOptions);
  // same control points, original parameter domains
  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(x, mapper, result);
  return result;
}

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePDEPatchPlanar(const gsMultiPatch<T> &mp,
                                           const gsDofMapper &mapper,
                                           const gsOptionList &options) {
  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);


  verboseLog("PDE-based parameterization construction ...\n",
             options.askInt("Verbose", 0));

  gsExprEvaluator<T> evaluator;
  gsExprAssembler<T> assembler;

  gsMultiBasis<T> mb(mp);
  assembler.setIntegrationDomain(mb.domain());

  gsBoundaryConditions<> bc;
  bc.setGeoMap(mp);
  for (auto bit = mp.bBegin(); bit != mp.bEnd(); ++bit) {
    bc.addCondition(*bit, condition_type::dirichlet, nullptr);
  }

  space space1 = assembler.getSpace(mb, d); // 1D space!!
  space1.setup(bc, dirichlet::homogeneous, 0);

  // Function for the Residual
  typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;
  typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;

  gsMultiPatch<T> mpSubstitute = mp;
  Residual_t Residual = [&assembler, &mapper, &mpSubstitute, &space1](
      gsVector<T> const &x) {
    // Convert free vector to MultiPatch
    convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);

    // Calculate geometry map
    geometryMap G = assembler.getMap(mpSubstitute);

    // Compute metric matrix and Hessian of G
    auto metricMat = jac(G).tr() * jac(G);
    auto hessG = hess(G);

    // Compute scale factor
    auto scale = metricMat.trace();

    // Compute Lx
    auto Lx = hessG % metricMat.adj() / scale.val();

    // Compute nonlinear system
    auto nonlinearSystem = frprod3(space1, Lx);

    // Assemble the system
    assembler.initSystem();
    assembler.assemble(nonlinearSystem);

    // Return the assembled system
    return assembler.rhs();
  };

  Jacobian_t Jacobian;
  int preconditionerType = options.askInt("AAPreconditionType", 0);
  switch (preconditionerType) {
    case 1:
      Jacobian = [&assembler, &mapper, &mpSubstitute, &space1](gsVector<T> const &x)
      {
        // diagonal Jacobian matrix as a preconditioner
        convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);
        geometryMap G = assembler.getMap(mpSubstitute);

        auto jacMat = jacScaledLxDiag(space1, G);

        assembler.initSystem();
        assembler.assemble(jacMat);
        return assembler.matrix();
      };
      break;

    case 2:
      Jacobian = [&assembler, &mapper, &mpSubstitute, &space1](gsVector<T> const &x)
      {
        // diagonal-block Jacobian matrix as a preconditioner
        convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);
        geometryMap G = assembler.getMap(mpSubstitute);

        auto jacMat = jacScaledLxDiagBlock(space1, G);

        assembler.initSystem();
        assembler.assemble(jacMat);
        return assembler.matrix();
      };
      break;

    default:
      // analytical Jacobian matrix as a preconditioner
      Jacobian = [&assembler, &mapper, &mpSubstitute, &space1](gsVector<T> const &x)
      {
        convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);
        geometryMap G = assembler.getMap(mpSubstitute);

        auto jacMat = jacScaledLx(space1, G);

        assembler.initSystem();
        assembler.assemble(jacMat);
        return assembler.matrix();
      };
  }

  preAAParam<T> param;
  param.m = options.askInt("AAwindowsize", 5);
  param.usePreconditioning = true;
  param.updatePreconditionerStep = options.askInt("N_update", 5);
  param.epsilon = 1e-5;

  preAApp::AndersonAcceleration<T> AASolver(param);
  gsStopwatch stageTimer;
  gsVector<T> solVector = AASolver.compute(initialGuessVector,Residual, Jacobian);
  verboseLog(aaSummary("H2", AASolver, stageTimer.stop()),
             options.askInt("Verbose", 0));

  if (options.askSwitch("needPDEH1", true)) {
    solVector = planarH1Stage(mp, mapper, solVector, options);
  }

  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(solVector, mapper, result);

  verboseLog("Finished!", options.askInt("Verbose", 0));

  return result;
}

template<short_t d, typename T>
gsVector<T>
gsBarrierCore<d, T>::planarH1Stage(const gsMultiPatch<T> &mp,
                                   const gsDofMapper &mapper,
                                   const gsVector<T> &x,
                                   const gsOptionList &options) {
  gsExprAssembler<T> assembler;
  gsMultiBasis<T> mb(mp);
  assembler.setIntegrationDomain(mb.domain());
  gsBoundaryConditions<> bc;
  bc.setGeoMap(mp);
  for (auto bit = mp.bBegin(); bit != mp.bEnd(); ++bit) {
    bc.addCondition(*bit, condition_type::dirichlet, nullptr);
  }
  space space1 = assembler.getSpace(mb, d);
  space1.setup(bc, dirichlet::homogeneous, 0);
  gsMultiPatch<T> mpSubstitute = mp;

  typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;
  typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;
  Jacobian_t Jacobian;
  preAAParam<T> param;
  param.m = options.askInt("AAwindowsize", 5);
  param.usePreconditioning = true;
  param.updatePreconditionerStep = options.askInt("N_update", 5);
  param.epsilon = 1e-5;

  verboseLog("\nStart parameterization improvement by H1 discretization...",options.askInt("Verbose", 0));
  Residual_t Residual = [&assembler, &mapper, &mpSubstitute, &space1](gsVector<T> const &u)
  {
    // H1 discretization
    convertFreeVectorToMultiPatch<T>(u, mapper, mpSubstitute);
    geometryMap G = assembler.getMap(mpSubstitute);

    auto invJacMat = jac(G).inv();
    auto jacU = jac(space1) * invJacMat;
    auto nonlinearSystem = jacU % invJacMat;

    assembler.initSystem();
    assembler.assemble(nonlinearSystem);
    return assembler.rhs();
  };

  index_t preconditionerType = options.askInt("AAPreconditionType", 0);
  switch (preconditionerType) {
    case 1:
      Jacobian = [&assembler, &mapper, &mpSubstitute, &space1](gsVector<T> const &u)
      {
        // diagonal block Jacobian matrix as a preconditioner
        convertFreeVectorToMultiPatch<T>(u, mapper, mpSubstitute);
        geometryMap G = assembler.getMap(mpSubstitute);

        auto jacMat = jacScaledLxH1DiagBlock(space1, G);

        assembler.initSystem();
        assembler.assemble(jacMat);

        return assembler.matrix();
      };
      break;

    default:
      Jacobian = [&assembler, &mapper, &mpSubstitute,&space1](gsVector<T> const &u)
      {
        // analytical Jacobian matrix as a preconditioner
        convertFreeVectorToMultiPatch<T>(u, mapper, mpSubstitute);
        geometryMap G = assembler.getMap(mpSubstitute);

        auto jacMat = jacScaledLxH1(space1, G);

        assembler.initSystem();
        assembler.assemble(jacMat);

        return assembler.matrix();
      };
  }

  preAApp::AndersonAcceleration<T> AASolver2(param);
  return AASolver2.compute(x, Residual, Jacobian);
}

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePDEPatchGeneric(const gsMultiPatch<T> &mp,
                                           const gsDofMapper &mapper,
                                           const gsOptionList &options) {
  GISMO_ENSURE(d == 2 || d == 3, "The PDE-based method is implemented for "
                                 "dimension 2 and 3 only.");

  const index_t verbose = options.askInt("Verbose", 0);
  verboseLog("PDE-based parameterization construction ...\n", verbose);

  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);

  // H2 stage in its own scope: its assembler, Jacobian and preconditioner
  // are released before the next stage
  gsVector<T> solVector;
  {
    typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;
    typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;

    // Elliptic grid generation equations in non-divergence form (H2),
    // assembled element-wise with the default (Gauss) quadrature of
    // gsExprAssembler, as the expressions scaledEGGRes / scaledEGGJac
    gsRegularizedPDEAssembler<d, T> h2Assembler(
        mp, mapper, gsExprAssembler<T>::defaultOptions(), 1);
    Residual_t Residual = [&](gsVector<T> const &x) {
      gsVector<T> residual;
      h2Assembler.setFreeVector(x);
      h2Assembler.h2Residual(residual);
      return residual;
    };

    // Preconditioner: 0: analytical Jacobian, 1: its diagonal,
    // 2: its diagonal blocks
    index_t preconditionerType = options.askInt("AAPreconditionType", 0);
    if (preconditionerType < 0 || preconditionerType > 2) {
      preconditionerType = 0;
    }
    Jacobian_t Jacobian = [&](gsVector<T> const &x) {
      gsSparseMatrix<T> jacobian;
      h2Assembler.setFreeVector(x);
      h2Assembler.h2Jacobian(preconditionerType, jacobian);
      return jacobian;
    };

    preAAParam<T> param;
    param.m = options.askInt("AAwindowsize", 5);
    param.usePreconditioning = true;
    param.updatePreconditionerStep = options.askInt("N_update", 5);
    param.epsilon = 1e-5;
    if (useIterativeSolver(mapper.freeSize(), options)) {
      param.linearSolver = gmresPreconditioner(options, true);
      param.blocks = componentOfDofs(mp, mapper);
      param.nodes = nodeOfDofs(mp, mapper);
      param.innerMaxIterations = options.getInt("PDE.InnerMaxIterations");
    }

    preAApp::AndersonAcceleration<T> AASolver(param);
    if (verbose > 1) { AASolver.enableIterationInfoPrinting(); }
    else { AASolver.disableIterationInfoPrinting(); }
    gsStopwatch stageTimer;
    solVector = AASolver.compute(initialGuessVector, Residual, Jacobian);
    verboseLog(aaSummary("H2", AASolver, stageTimer.stop()), verbose);
  }
  if (!solVector.allFinite()) { solVector = initialGuessVector; }

  if (options.getSwitch("PDE.Garanzha")) {
    // Robust second stage: untangling + regularized H1 discretization
    solVector = computeRegularizedH1(mp, mapper, solVector, options);
  } else if (options.askSwitch("needPDEH1", true)) {
    gsWarn << "For volumes, the H1 discretization is only available with "
              "Garanzha's regularization (option PDE.Garanzha); the result of "
              "the H2 discretization is returned.\n";
  }

  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(solVector, mapper, result);

  verboseLog("Finished!", verbose);

  return result;
}

template<short_t d, typename T>
std::vector<index_t>
gsBarrierCore<d, T>::componentOfDofs(const gsMultiPatch<T> &mp,
                                     const gsDofMapper &mapper) {
  std::vector<index_t> component(mapper.freeSize(), 0);
  for (size_t ip = 0; ip != mp.nPatches(); ++ip) {
    for (index_t c = 0; c != mp.targetDim(); ++c) {
      for (index_t i = 0; i != index_t(mapper.patchSize(ip, c)); ++i) {
        if (mapper.is_free(i, ip, c)) { component[mapper.index(i, ip, c)] = c; }
      }
    }
  }
  return component;
}

template<short_t d, typename T>
std::vector<index_t>
gsBarrierCore<d, T>::nodeOfDofs(const gsMultiPatch<T> &mp,
                                const gsDofMapper &mapper) {
  // key of a control point: its smallest free unknown
  std::vector<index_t> key(mapper.freeSize(), -1);
  std::vector<index_t> dofs;
  for (size_t ip = 0; ip != mp.nPatches(); ++ip) {
    const index_t numBasis = mapper.patchSize(ip, 0);
    for (index_t i = 0; i != numBasis; ++i) {
      dofs.clear();
      for (index_t c = 0; c != mp.targetDim(); ++c) {
        if (i < index_t(mapper.patchSize(ip, c)) && mapper.is_free(i, ip, c)) {
          dofs.push_back(mapper.index(i, ip, c));
        }
      }
      if (dofs.empty()) { continue; }
      const index_t k = *std::min_element(dofs.begin(), dofs.end());
      for (index_t dof : dofs) {
        key[dof] = (key[dof] < 0) ? k : std::min(key[dof], k);
      }
    }
  }
  // consecutive numbers
  std::vector<index_t> number(mapper.freeSize(), -1), node(mapper.freeSize());
  index_t numNodes = 0;
  for (size_t i = 0; i != key.size(); ++i) {
    const index_t k = key[i] < 0 ? index_t(i) : key[i];
    if (number[k] < 0) { number[k] = numNodes++; }
    node[i] = number[k];
  }
  return node;
}

template<short_t d, typename T>
bool gsBarrierCore<d, T>::couplesPatches(const gsMultiPatch<T> &mp,
                                         const gsDofMapper &mapper) {
  if (mp.nPatches() < 2) { return false; }
  std::vector<index_t> owner(mapper.freeSize(), -1);
  for (size_t ip = 0; ip != mp.nPatches(); ++ip)
    for (index_t c = 0; c != mp.targetDim(); ++c)
      for (index_t i = 0; i != index_t(mapper.patchSize(ip, c)); ++i) {
        if (!mapper.is_free(i, ip, c)) { continue; }
        index_t &o = owner[mapper.index(i, ip, c)];
        if (o >= 0 && o != index_t(ip)) { return true; }
        o = ip;
      }
  return false;
}

template<short_t d, typename T>
gsMultiPatch<T> gsBarrierCore<d, T>::chordScaled(const gsMultiPatch<T> &mp) {
  const size_t np = mp.nPatches();
  // (patch, direction) -> class; union-find
  std::vector<index_t> parent(np * d);
  for (size_t i = 0; i != parent.size(); ++i) { parent[i] = i; }
  std::function<index_t(index_t)> find = [&](index_t i) {
    while (parent[i] != i) { i = parent[i] = parent[parent[i]]; }
    return i;
  };
  // physical length along each direction: mean length of the control
  // polygons of the 2^(d-1) patch edges in that direction
  std::vector<T> length(np * d, T(0));
  std::vector<bool> tensor(np, false);
  for (size_t ip = 0; ip != np; ++ip) {
    const gsTensorBSpline<d, T> *patch =
        dynamic_cast<const gsTensorBSpline<d, T> *>(&mp.patch(ip));
    if (!patch) { continue; }
    tensor[ip] = true;
    const gsTensorBSplineBasis<d, T> &basis = patch->basis();
    gsVector<index_t, d> size;
    for (short_t a = 0; a != d; ++a) { size(a) = basis.size(a); }
    for (short_t a = 0; a != d; ++a) {
      T sum = 0;
      index_t numEdges = 0;
      // corner pattern of the other directions
      for (index_t corner = 0; corner != (1 << (d - 1)); ++corner) {
        gsVector<index_t, d> idx;
        index_t bit = 0;
        for (short_t b = 0; b != d; ++b) {
          if (b == a) { continue; }
          idx(b) = ((corner >> bit++) & 1) ? size(b) - 1 : 0;
        }
        for (index_t i = 0; i + 1 < size(a); ++i) {
          idx(a) = i;
          const index_t i0 = basis.index(idx);
          idx(a) = i + 1;
          sum += (patch->coefs().row(basis.index(idx)) - patch->coefs().row(i0)).norm();
        }
        ++numEdges;
      }
      length[ip * d + a] = sum / numEdges;
    }
  }
  for (auto it = mp.topology().iBegin(); it != mp.topology().iEnd(); ++it) {
    const patchSide &first = it->first(), &second = it->second();
    const gsVector<index_t> dirMap = it->dirMap(first);
    for (short_t a = 0; a != d; ++a) {
      if (a == first.direction()) { continue; }
      const index_t r1 = find(first.patch * d + a);
      const index_t r2 = find(second.patch * d + dirMap(a));
      if (r1 != r2) { parent[r1] = r2; }
    }
  }
  std::vector<T> sum(np * d, T(0));
  std::vector<index_t> count(np * d, 0);
  for (size_t i = 0; i != np * d; ++i) {
    if (!tensor[i / d] || !(length[i] > 0)) { continue; }
    sum[find(i)] += length[i];
    ++count[find(i)];
  }
  // lengths relative to the longest chord (affineTransformTo needs
  // intervals longer than 1e-5)
  T longest = 0;
  for (size_t i = 0; i != np * d; ++i) {
    if (count[i] > 0) { longest = std::max(longest, sum[i] / count[i]); }
  }
  gsMultiPatch<T> result = mp;
  for (size_t ip = 0; ip != np; ++ip) {
    if (!tensor[ip]) { continue; }
    const gsTensorBSpline<d, T> &patch =
        static_cast<const gsTensorBSpline<d, T> &>(mp.patch(ip));
    std::vector<gsKnotVector<T> > knots(d);
    for (short_t a = 0; a != d; ++a) {
      const index_t r = find(ip * d + a);
      knots[a] = patch.basis().knots(a);
      if (count[r] > 0) {
        knots[a].affineTransformTo(0, sum[r] / count[r] / longest);
      }
    }
    gsTensorBSplineBasis<d, T> basis(knots);
    result.patch(ip) = gsTensorBSpline<d, T>(basis, patch.coefs());
  }
  return result;
}

template<short_t d, typename T>
int gsBarrierCore<d, T>::gmresPreconditioner(const gsOptionList &options,
                                             const bool h2Stage) {
  switch (options.getInt("PDE.Preconditioner")) {
    case Preconditioner::BlockGaussSeidelLU:
      return preAAParam<T>::BlockGaussSeidelLU;
    case Preconditioner::NodalAMG:
      return preAAParam<T>::NodalAMG;
    case Preconditioner::BlockGaussSeidelAMG:
      return preAAParam<T>::BlockGaussSeidelAMG;
    default:
      return h2Stage ? preAAParam<T>::BlockGaussSeidelAMG
                     : preAAParam<T>::NodalAMG;
  }
}

template<short_t d, typename T>
bool gsBarrierCore<d, T>::useIterativeSolver(const index_t numDofs,
                                             const gsOptionList &options) {
  switch (options.getInt("PDE.LinearSolver")) {
    case LinearSolver::Direct: return false;
    case LinearSolver::Iterative: return true;
    default: return numDofs > options.getInt("PDE.DirectSolverLimit");
  }
}

/**
 * @brief Element-wise assembly of the quantities used in the regularized
 * second stage of the PDE-based method (planar domains and volumes):
 * Garanzha's regularized Winslow functional with its gradient and
 * Gauss-Newton Hessian, the regularized H1 residual with its Jacobian, and
 * the minimum of det(J) over the quadrature nodes.
 *
 * Away from nodes with degenerate J, the same quantities are available as
 * expressions (garanzhaWinslowEnergy, garanzhaWinslowGrad/Hess,
 * regularizedH1Res/Jac); the expressions skip a node where J is degenerate
 * (detIsControllable), the assembler keeps it if the gradients of the free
 * functions span all directions (see monitored()). Here the element matrices
 * are computed with matrix-matrix products over all quadrature nodes of an
 * element, the sparsity patterns are computed once, and the Hessian is stored
 * as its lower triangle only.
 *
 * Basis functions: for tensor-product B-spline bases, the univariate
 * B-splines of every knot span and their first and second derivatives are
 * tabulated at the univariate quadrature nodes (a few kB in total). The
 * Jacobian matrices of the map at the nodes (energy, minimum of det(J), H1
 * residual) and the projection of the H1 fluxes onto the gradients of the
 * basis functions are computed from these tables by sum factorization; where
 * all derivatives of the active functions are needed (Hessian, Jacobians,
 * H2), they are formed from the tables (same products, in the same order, as
 * gsTensorBasis). For other bases, the derivatives at the quadrature nodes
 * are cached within a memory budget (\a cacheBytes).
 *
 * Cost per element, for degree p and q nodes per direction: the Jacobians
 * of the map at all q^d nodes by sum factorization take
 * O(d^2 (p+1) max(p+1, q)^d) operations instead of O(d^2 (p+1)^d q^d) with
 * the derivatives of all active functions; element matrices (Hessians,
 * Jacobians) take O(d^2 (p+1)^(2d) q^d), as for any element-wise assembly.
 *
 * det(J) is monitored (and the regularized quantities are integrated) at the
 * quadrature nodes at which it depends on the free control points; this is
 * determined once, for the given geometry: on the boundary of the parameter
 * domain it only depends on the fixed control points.
 *
 * Besides the regularized stages, the class assembles the residual and the
 * Jacobian of the (scaled) elliptic grid generation equations in
 * non-divergence form (H2 discretization, the same as scaledEGGRes /
 * scaledEGGJac); for these, use an instance with the quadrature rule of the
 * H2 stage.
 */
template<short_t d, typename T>
class gsRegularizedPDEAssembler {
 public:
  typedef gsBarrierCore<d, T> Core;
  typedef gsEigen::Matrix<T, d, d> MatrixD;

  /// \a orientations: for every patch, +1 or -1, the sign of det(J) of a
  /// bijective parameterization of that patch (patches of a multi-patch may
  /// be oriented differently). \a extraNodes (optional): for every element
  /// (numbered patch by patch in the order of the domain iterators), the
  /// number of additional Gauss-Lobatto nodes per direction to use there
  gsRegularizedPDEAssembler(const gsMultiPatch<T> &mp,
                            const gsDofMapper &mapper,
                            const gsOptionList &quOptions,
                            const std::vector<T> &orientations,
                            const std::vector<index_t> &extraNodes = std::vector<index_t>(),
                            const size_t cacheBytes = size_t(1) << 30)
      : m_mp(mp), m_mb(mp), m_mapper(mapper), m_orientations(orientations),
        m_numDofs(mapper.freeSize()) {
    GISMO_ENSURE(m_orientations.size() == mp.nPatches(),
                 "gsRegularizedPDEAssembler: one orientation per patch expected");
    setup(quOptions, extraNodes, cacheBytes);
  }

  /// Same orientation (+1 or -1) for all patches
  gsRegularizedPDEAssembler(const gsMultiPatch<T> &mp,
                            const gsDofMapper &mapper,
                            const gsOptionList &quOptions,
                            const T orientation,
                            const std::vector<index_t> &extraNodes = std::vector<index_t>(),
                            const size_t cacheBytes = size_t(1) << 30)
      : m_mp(mp), m_mb(mp), m_mapper(mapper),
        m_orientations(mp.nPatches(), orientation),
        m_numDofs(mapper.freeSize()) {
    setup(quOptions, extraNodes, cacheBytes);
  }

 private:
  /// Elements, quadrature nodes, tables and monitored nodes (constructors)
  void setup(const gsOptionList &quOptions, const std::vector<index_t> &extraNodes,
             const size_t cacheBytes) {
    size_t cached = 0;
    gsMatrix<T> ders, coefs;
    MatrixD jacGeom, cof;
    index_t id = 0;
    for (size_t ip = 0; ip != m_mb.nBases(); ++ip) {
      const gsBasis<T> &basis = m_mb.basis(ip);
      const gsTensorBSplineBasis<d, T> *tensorBasis =
          dynamic_cast<const gsTensorBSplineBasis<d, T> *>(&basis);
      // univariate tables of this patch: (direction, first function, nodes)
      std::map<std::array<index_t, 3>, index_t> tableOf;
      // quadrature rules by number of additional nodes
      std::map<index_t, memory::shared_ptr<gsQuadRule<T> > > rules;
      for (auto &domElem : basis.domain()->allElements()) {
        const index_t extra = (size_t(id) < extraNodes.size() ? extraNodes[id] : 0);
        if (!rules.count(extra)) {
          gsOptionList elementOptions = quOptions;
          if (extra > 0) {
            elementOptions.setInt("quRule", gsQuadrature::GaussLobatto);
            elementOptions.setInt("quB", quOptions.getInt("quB") + extra);
          }
          rules[extra] = memory::shared_ptr<gsQuadRule<T> >(
              gsQuadrature::getPtr(basis, elementOptions).release());
        }
        Element e;
        e.id = id++;
        e.patch = ip;
        e.lower = domElem.lowerCorner();
        e.upper = domElem.upperCorner();
        rules[extra]->mapTo(e.lower, e.upper, e.nodes, e.weights);
        basis.active_into(domElem.centerPoint(), e.actives);
        const index_t A = e.actives.rows();
        e.dofs.resize(d * A);
        bool anyFree = false;
        for (index_t c = 0; c != d; ++c) {
          for (index_t i = 0; i != A; ++i) {
            const index_t dof = e.actives(i, 0);
            e.dofs[c * A + i] = m_mapper.is_free(dof, ip, c)
                                ? m_mapper.index(dof, ip, c) : -1;
            anyFree = anyFree || e.dofs[c * A + i] >= 0;
          }
        }
        if (!anyFree) { continue; } // nothing depends on this element

        e.table.fill(-1);
        if (tensorBasis) { setupTables(*tensorBasis, tableOf, e); }
        if (e.table[0] >= 0) {
          tensorBasisValues(e, nullptr, &ders, nullptr);
        } else {
          elementDerivatives(basis, e.nodes, ders);
        }
        localCoefs(e, coefs);
        // det(J) is monitored at the nodes where it depends on free dofs
        e.monitor.assign(e.numNodes(), 0);
        for (index_t k = 0; k != e.numNodes(); ++k) {
          const gsAsConstMatrix<T> dB = ders.reshapeCol(k, d, A);
          jacGeom.noalias() = coefs.transpose() * dB.transpose();
          cofactorMatrix(jacGeom, cof);
          e.monitor[k] = monitored(e, cof, dB);
        }
        if (e.table[0] >= 0) {
          e.nodes.resize(0, 0); // formed from the tables when needed
        } else {
          const size_t bytes = ders.size() * sizeof(T);
          if (cached + bytes <= cacheBytes) {
            e.ders = ders;
            cached += bytes;
          }
        }
        m_elements.push_back(give(e));
      }
    }
    m_numAllElements = id;
  }

 public:
  /// Total number of elements (including those without free dofs)
  index_t numElements() const { return m_numAllElements; }

  /// For every element, the number of nodes of the given (check) rule at
  /// which det(J) depends on the free dofs but orientation * det(J) <= 0
  std::vector<index_t> negativeCheckPoints(const gsOptionList &checkOptions) const {
    std::vector<index_t> result(m_numAllElements, 0);
    gsMatrix<T> nodes, ders, coefs;
    gsVector<T> weights;
    MatrixD jacGeom, cof;
    std::vector<memory::shared_ptr<gsQuadRule<T> > > rules(m_mb.nBases());
    for (const Element &e : m_elements) {
      if (!rules[e.patch]) {
        rules[e.patch] = memory::shared_ptr<gsQuadRule<T> >(
            gsQuadrature::getPtr(m_mb.basis(e.patch), checkOptions).release());
      }
      rules[e.patch]->mapTo(e.lower, e.upper, nodes, weights);
      elementDerivatives(m_mb.basis(e.patch), nodes, ders);
      localCoefs(e, coefs);
      const index_t A = e.actives.rows();
      for (index_t k = 0; k != nodes.cols(); ++k) {
        const gsAsConstMatrix<T> dB = ders.reshapeCol(k, d, A);
        jacGeom.noalias() = coefs.transpose() * dB.transpose();
        cofactorMatrix(jacGeom, cof);
        if (!monitored(e, cof, dB)) { continue; }
        const T det = orientation(e) * jacGeom.row(0).dot(cof.row(0));
        if (!(det > 0)) { ++result[e.id]; }
      }
    }
    return result;
  }

  /// Sets the free control points
  void setFreeVector(const gsVector<T> &x) {
    convertFreeVectorToMultiPatch<T>(x, m_mapper, m_mp);
  }

  /// Minimum of orientation * det(J) over the monitored quadrature nodes
  T minDet() const {
    std::vector<T> partial(numThreads(), std::numeric_limits<T>::max());
    // buffers of a thread
    struct Work {
      gsMatrix<T> coefs, jacobians, work;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &coefs = ws.coefs;
      auto &jacobians = ws.jacobians;
      auto &work = ws.work;
      localCoefs(e, coefs);
      nodeJacobians(e, coefs, jacobians, work);
      T &result = partial[thread];
      for (index_t k = 0; k != e.numNodes(); ++k) {
        if (!e.monitor[k]) { continue; }
        const T det = orientation(e)
            * MatrixD(jacobians.middleCols(k * d, d)).determinant();
        if (det != det) { result = det; return; } // NaN must not go unnoticed
        if (det < result) { result = det; }
      }
    }, false);
    T result = std::numeric_limits<T>::max();
    for (const T value : partial) {
      if (value != value) { return value; }
      result = std::min(result, value);
    }
    return result;
  }

  /// Garanzha's regularized Winslow functional: int |J|^2 / chi^(2/d)
  /// (energy form Core::EnergyForm::Conformal), or int |cof J|^2 / chi
  /// (Core::EnergyForm::InverseDirichlet, see setEnergyForm)
  T energy(const T eps) const {
    std::vector<T> partial(numThreads(), T(0));
    // buffers of a thread
    struct Work {
      gsMatrix<T> coefs, jacobians, work;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &coefs = ws.coefs;
      auto &jacobians = ws.jacobians;
      auto &work = ws.work;
      MatrixD jacGeom, cof;
      T &result = partial[thread];
      localCoefs(e, coefs);
      nodeJacobians(e, coefs, jacobians, work);
      for (index_t k = 0; k != e.numNodes(); ++k) {
        if (!e.monitor[k]) { continue; }
        jacGeom = jacobians.middleCols(k * d, d);
        cofactorMatrix(jacGeom, cof);
        const T det = orientation(e) * jacGeom.row(0).dot(cof.row(0));
        if (Core::EnergyForm::InverseDirichlet == m_energyForm && 3 == d) {
          result += e.weights(k) * cof.squaredNorm() / garanzhaChi(det, eps);
        } else {
          result += e.weights(k) * jacGeom.squaredNorm()
              / chiPower(garanzhaChi(det, eps));
        }
      }
    }, false);
    T result = 0;
    for (const T value : partial) { result += value; }
    return result;
  }

  /// Gradient and (lower triangle of the) Gauss-Newton Hessian of the
  /// regularized Winslow functional
  void gradientHessian(const T eps, gsVector<T> &gradient,
                       gsSparseMatrix<T> &hessian) const {
    gradient.setZero(m_numDofs);
    if (m_patternLower.outer.empty()) { buildPattern(true, m_patternLower); }
    makeMatrix(m_patternLower, hessian);
    // buffers of a thread
    struct Work {
      gsMatrix<T> buffer, coefs, P, Q, X, Y, Dstack, Dweighted, local, laplace;
      gsVector<T> localGradient;
      std::vector<index_t> nodes;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &buffer = ws.buffer;
      auto &coefs = ws.coefs;
      auto &P = ws.P;
      auto &Q = ws.Q;
      auto &X = ws.X;
      auto &Y = ws.Y;
      auto &Dstack = ws.Dstack;
      auto &Dweighted = ws.Dweighted;
      auto &local = ws.local;
      auto &laplace = ws.laplace;
      auto &localGradient = ws.localGradient;
      auto &nodes = ws.nodes;
      MatrixD jacGeom, cof;
      const gsMatrix<T> &ders = derivatives(e, buffer);
      localCoefs(e, coefs);
      const index_t A = e.actives.rows();
      monitoredNodes(e, nodes);
      const index_t nq = nodes.size();
      if (0 == nq) { return; }

      P.resize(d * A, nq);
      Q.resize(d * A, nq);
      Y.resize(d * A, 2 * nq);
      Dstack.resize(d * nq, A);
      Dweighted.resize(d * nq, A);
      localGradient.setZero(d * A);
      for (index_t n = 0; n != nq; ++n) {
        const index_t k = nodes[n];
        const gsAsConstMatrix<T> dB = ders.reshapeCol(k, d, A);
        jacGeom.noalias() = coefs.transpose() * dB.transpose();
        cofactorMatrix(jacGeom, cof);
        const T det = orientation(e) * jacGeom.row(0).dot(cof.row(0));
        const T sqNorm = jacGeom.squaredNorm();
        const T root = math::sqrt(det * det + eps * eps);
        const T chiPow = chiPower(garanzhaChi(det, eps));
        const T rho = 1.0 / root; // chi'/chi
        const T w = e.weights(k);

        // p.row(c) = grad(x_c) . grad(phi), q.row(c) = d(det J)/d(x_c)
        for (index_t c = 0; c != d; ++c) {
          P.col(n).segment(c * A, A).noalias() =
              (jacGeom.row(c) * dB).transpose();
          Q.col(n).segment(c * A, A).noalias() =
              orientation(e) * (cof.row(c) * dB).transpose();
        }
        localGradient.noalias() += (w / chiPow)
            * (2.0 * P.col(n) - (T(2) / d) * rho * sqNorm * Q.col(n));

        // Hessian: c1 (p q^T + q p^T) + c2 q q^T + 2/chi^(2/d) laplace
        const T wp = 2.0 * w / chiPow;
        const T c1 = -(T(1) / d) * rho * wp;
        const T c2 = (T(1) / (d * d)) * rho * rho * sqNorm * wp;
        Y.col(n) = c1 * Q.col(n);
        Y.col(nq + n) = c1 * P.col(n) + c2 * Q.col(n);
        Dstack.middleRows(d * n, d) = dB;
        Dweighted.middleRows(d * n, d) = wp * dB;
      }
      X.resize(d * A, 2 * nq);
      X.leftCols(nq) = P;
      X.rightCols(nq) = Q;
      local.noalias() = X * Y.transpose();
      laplace.noalias() = Dstack.transpose() * Dweighted;
      for (index_t c = 0; c != d; ++c) {
        local.block(c * A, c * A, A, A) += laplace;
      }
      scatterVector(e, localGradient, gradient);
      scatterMatrix(e, local, true, hessian);
    });
  }

  /// Form of the energy: 0: |J|^2 / chi^(2/d), 1: |cof J|^2 / chi (see
  /// garanzhaEnergyDerivatives; the same in 2D)
  void setEnergyForm(const index_t form) { m_energyForm = form; }
  /// Form of the energy (see setEnergyForm)
  index_t energyForm() const { return m_energyForm; }

  /// Gradient and lower triangle of a Hessian of the energy (form set by
  /// setEnergyForm): hessianType 0: Gauss-Newton, 1: exact, 2: exact,
  /// projected onto the positive semidefinite matrices at every node.
  /// (For form 0 with Gauss-Newton, gradientHessian is the faster variant.)
  void energyGradientHessian(const T eps, const index_t hessianType,
                             gsVector<T> &gradient,
                             gsSparseMatrix<T> &hessian) const {
    gradient.setZero(m_numDofs);
    if (m_patternLower.outer.empty()) { buildPattern(true, m_patternLower); }
    makeMatrix(m_patternLower, hessian);
    typedef gsEigen::Matrix<T, d * d, 1> VecJ;
    typedef gsEigen::Matrix<T, d * d, d * d> MatJ;
    // buffers of a thread
    struct Work {
      gsMatrix<T> buffer, coefs, jacobians, work, fluxes, localGradient,
          Dstack, Tcb, local;
      std::vector<MatJ, gsEigen::aligned_allocator<MatJ> > H;
      std::vector<index_t> nodes;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &buffer = ws.buffer;
      auto &coefs = ws.coefs;
      auto &jacobians = ws.jacobians;
      auto &work = ws.work;
      auto &fluxes = ws.fluxes;
      auto &localGradient = ws.localGradient;
      auto &Dstack = ws.Dstack;
      auto &Tcb = ws.Tcb;
      auto &local = ws.local;
      auto &H = ws.H;
      auto &nodes = ws.nodes;
      VecJ gradJ;
      MatrixD jacGeom, N;
      monitoredNodes(e, nodes);
      const index_t nq = nodes.size();
      if (0 == nq) { return; }
      localCoefs(e, coefs);
      nodeJacobians(e, coefs, jacobians, work);
      const gsMatrix<T> &ders = derivatives(e, buffer);
      const index_t A = e.actives.rows();
      fluxes.setZero(d, d * e.numNodes());
      H.resize(nq);
      Dstack.resize(d * nq, A);
      for (index_t n = 0; n != nq; ++n) {
        const index_t k = nodes[n];
        jacGeom = jacobians.middleCols(k * d, d);
        garanzhaEnergyDerivatives<d, T>(jacGeom, orientation(e), eps,
                                        m_energyForm, hessianType, &gradJ, &H[n]);
        H[n] *= e.weights(k);
        for (short_t c = 0; c != d; ++c)
          for (short_t a = 0; a != d; ++a) {
            fluxes(c, k * d + a) = e.weights(k) * gradJ(c * d + a);
          }
        Dstack.middleRows(d * n, d) = ders.reshapeCol(k, d, A);
      }
      testGradients(e, fluxes, localGradient, work);
      gsAsConstVector<T> gradientVector(localGradient.data(), localGradient.size());
      scatterVector(e, gradientVector, gradient);

      // block (c, b) = sum_n dB_n^T N_n dB_n, N_n(a, e') = H_n(c d + a, b d + e')
      local.resize(d * A, d * A);
      Tcb.resize(d * nq, A);
      for (short_t c = 0; c != d; ++c) {
        for (short_t b = 0; b <= c; ++b) {
          for (index_t n = 0; n != nq; ++n) {
            N = H[n].template block<d, d>(c * d, b * d);
            Tcb.middleRows(d * n, d).noalias() = N * Dstack.middleRows(d * n, d);
          }
          local.block(c * A, b * A, A, A).noalias() = Dstack.transpose() * Tcb;
          if (b != c) {
            local.block(b * A, c * A, A, A) = local.block(c * A, b * A, A, A).transpose();
          }
        }
      }
      scatterMatrix(e, local, true, hessian);
    });
  }

  /// Exponent q of chi in the H1 residual, int grad(phi) . adj(g) e_c / chi^q
  void setH1Power(const index_t q) {
    GISMO_ENSURE(q == 1 || q == 2, "The exponent must be 1 or 2.");
    m_h1Power = q;
  }

  /// Residual of the regularized H1 discretization
  void h1Residual(const T eps, gsVector<T> &residual) const {
    residual.setZero(m_numDofs);
    // buffers of a thread
    struct Work {
      gsMatrix<T> coefs, jacobians, fluxes, local, work;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &coefs = ws.coefs;
      auto &jacobians = ws.jacobians;
      auto &fluxes = ws.fluxes;
      auto &local = ws.local;
      auto &work = ws.work;
      MatrixD jacGeom, cof;
      localCoefs(e, coefs);
      nodeJacobians(e, coefs, jacobians, work);
      // fluxes(c, k * d + a) = K(c, a) at node k (zero if not monitored);
      // the residual is  sum_k grad(phi_i) . K e_c
      fluxes.setZero(d, d * e.numNodes());
      for (index_t k = 0; k != e.numNodes(); ++k) {
        if (!e.monitor[k]) { continue; }
        jacGeom = jacobians.middleCols(k * d, d);
        cofactorMatrix(jacGeom, cof);
        const T det = orientation(e) * jacGeom.row(0).dot(cof.row(0));
        const T chi = garanzhaChi(det, eps);
        const T chiPow = (2 == m_h1Power ? chi * chi : chi);
        fluxes.middleCols(k * d, d).noalias() =
            (cof.transpose() * cof) * (e.weights(k) / chiPow);
      }
      testGradients(e, fluxes, local, work);
      gsAsConstVector<T> localVector(local.data(), local.size());
      scatterVector(e, localVector, residual);
    });
  }

  /// Jacobian of the residual of the regularized H1 discretization
  /// (optionally only its diagonal blocks)
  void h1Jacobian(const T eps, const bool diagBlocksOnly,
                  gsSparseMatrix<T> &jacobian) const {
    if (m_patternFull.outer.empty()) { buildPattern(false, m_patternFull); }
    makeMatrix(m_patternFull, jacobian);
    jacobian.coeffs().setZero();

    // dKs[n][b * d + a]: variation of K at node n for the perturbation of
    // the coordinate b by a function with gradient e_a
    // buffers of a thread
    struct Work {
      gsMatrix<T> buffer, coefs, Dstack, Tcb, local;
      std::vector<index_t> nodes;
      std::vector<MatrixD, gsEigen::aligned_allocator<MatrixD> > dKs;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &buffer = ws.buffer;
      auto &coefs = ws.coefs;
      auto &Dstack = ws.Dstack;
      auto &Tcb = ws.Tcb;
      auto &local = ws.local;
      auto &nodes = ws.nodes;
      auto &dKs = ws.dKs;
      MatrixD jacGeom, cof, dCof, K, dK;
      const gsMatrix<T> &ders = derivatives(e, buffer);
      localCoefs(e, coefs);
      const index_t A = e.actives.rows();
      monitoredNodes(e, nodes);
      const index_t nq = nodes.size();
      if (0 == nq) { return; }

      Dstack.resize(d * nq, A);
      dKs.resize(nq * d * d);
      gsEigen::Matrix<T, d, 1> unit;
      for (index_t n = 0; n != nq; ++n) {
        const index_t k = nodes[n];
        const gsAsConstMatrix<T> dB = ders.reshapeCol(k, d, A);
        Dstack.middleRows(d * n, d) = dB;
        jacGeom.noalias() = coefs.transpose() * dB.transpose();
        cofactorMatrix(jacGeom, cof);
        const T det = orientation(e) * jacGeom.row(0).dot(cof.row(0));
        const T root = math::sqrt(det * det + eps * eps);
        const T chi = garanzhaChi(det, eps);
        const T w = 1.0 / (2 == m_h1Power ? chi * chi : chi);
        K.noalias() = (cof.transpose() * cof) * w;
        for (index_t b = 0; b != d; ++b) {
          for (index_t a = 0; a != d; ++a) {
            unit.setZero();
            unit(a) = 1;
            cofactorDerivative(jacGeom, b, unit, dCof);
            // d(1/chi^q) = dInvChiq / chi^q, chi'/chi = 1/sqrt(D^2+eps^2)
            const T dInvChiq = -T(m_h1Power) * orientation(e) * cof(b, a) / root;
            dK.noalias() = dCof.transpose() * cof;
            dK += MatrixD(dK.transpose());
            dKs[n * d * d + b * d + a] =
                e.weights(k) * (w * dK + dInvChiq * K);
          }
        }
      }

      local.setZero(d * A, d * A);
      Tcb.resize(d * nq, A);
      for (index_t c = 0; c != d; ++c) {
        for (index_t b = 0; b != d; ++b) {
          if (diagBlocksOnly && c != b) { continue; }
          // block (c,b) = sum_n dB_n^T N_n dB_n with N_n(e,a) = dK[b][a](c,e)
          for (index_t n = 0; n != nq; ++n) {
            MatrixD N;
            for (index_t a = 0; a != d; ++a) {
              N.col(a) = dKs[n * d * d + b * d + a].row(c).transpose();
            }
            Tcb.middleRows(d * n, d).noalias() =
                N * Dstack.middleRows(d * n, d);
          }
          local.block(c * A, b * A, A, A).noalias() =
              Dstack.transpose() * Tcb;
        }
      }
      scatterMatrix(e, local, false, jacobian);
    });
  }

  /// Residual of the scaled elliptic grid generation equations in
  /// non-divergence form (H2):
  ///   R_{c,i} = int phi_i (adj(g) : Hess(x_c)) / trace(adj(g)) dxi
  void h2Residual(gsVector<T> &residual) const {
    residual.setZero(m_numDofs);
    const index_t nd2 = d * (d + 1) / 2;
    // buffers of a thread
    struct Work {
      gsMatrix<T> values, deriv, deriv2, coefs;
      gsVector<T> local;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &values = ws.values;
      auto &deriv = ws.deriv;
      auto &deriv2 = ws.deriv2;
      auto &coefs = ws.coefs;
      auto &local = ws.local;
      gsEigen::Matrix<T, d * (d + 1) / 2, 1> weights;
      MatrixD jacGeom, cof, adjMetric;
      allBasisValues(e, values, deriv, deriv2);
      localCoefs(e, coefs);
      const index_t A = e.actives.rows();
      local.setZero(d * A);
      for (index_t k = 0; k != e.numNodes(); ++k) {
        jacGeom.noalias() =
            coefs.transpose() * deriv.reshapeCol(k, d, A).transpose();
        cofactorMatrix(jacGeom, cof);
        adjMetric.noalias() = cof.transpose() * cof;
        packSymmetric(adjMetric, weights);
        const T scale = adjMetric.trace();
        const gsAsConstMatrix<T> D2 = deriv2.reshapeCol(k, nd2, A);
        for (index_t c = 0; c != d; ++c) {
          const T Lx = weights.dot(D2 * coefs.col(c)) / scale;
          local.segment(c * A, A) += (e.weights(k) * Lx) * values.col(k);
        }
      }
      scatterVector(e, local, residual);
    });
  }

  /// Jacobian of h2Residual; mode 0: full matrix, 1: its diagonal only,
  /// 2: its diagonal blocks (coordinate components) only
  void h2Jacobian(const index_t mode, gsSparseMatrix<T> &jacobian) const {
    if (m_patternFull.outer.empty()) { buildPattern(false, m_patternFull); }
    makeMatrix(m_patternFull, jacobian);
    jacobian.coeffs().setZero();
    const index_t nd2 = d * (d + 1) / 2;
    typedef gsEigen::Matrix<T, d * (d + 1) / 2, 1> VectorH;
    // dL[c * d + b]: row n = variation of Lx(c) at node n with respect to
    // coordinate b of the active functions
    // buffers of a thread
    struct Work {
      gsMatrix<T> values, deriv, deriv2, coefs, local, weighted;
      std::vector<gsMatrix<T> > dL = std::vector<gsMatrix<T> >(d * d);
      gsEigen::Matrix<T, 1, Dynamic> commonTerm;
    };
    forEachElement<Work>([&](const Element &e, const int thread, Work &ws) {
      GISMO_UNUSED(thread);
      auto &values = ws.values;
      auto &deriv = ws.deriv;
      auto &deriv2 = ws.deriv2;
      auto &coefs = ws.coefs;
      auto &local = ws.local;
      auto &weighted = ws.weighted;
      auto &dL = ws.dL;
      VectorH weights, dWeights, deriv2Geom[d];
      MatrixD jacGeom, cof, dCof, adjMetric, dAdjMetric, value;
      auto &commonTerm = ws.commonTerm;
      gsEigen::Matrix<T, d, 1> unit, Lx;
      allBasisValues(e, values, deriv, deriv2);
      localCoefs(e, coefs);
      const index_t A = e.actives.rows(), nq = e.numNodes();
      for (auto &m : dL) { m.resize(nq, A); }
      weighted.resize(A, nq);
      for (index_t k = 0; k != nq; ++k) {
        const gsAsConstMatrix<T> dB = deriv.reshapeCol(k, d, A);
        const gsAsConstMatrix<T> D2 = deriv2.reshapeCol(k, nd2, A);
        jacGeom.noalias() = coefs.transpose() * dB.transpose();
        cofactorMatrix(jacGeom, cof);
        adjMetric.noalias() = cof.transpose() * cof;
        packSymmetric(adjMetric, weights);
        const T scale = adjMetric.trace();
        for (index_t c = 0; c != d; ++c) {
          deriv2Geom[c].noalias() = D2 * coefs.col(c);
          Lx(c) = weights.dot(deriv2Geom[c]) / scale;
        }
        // adj(g) : Hess(phi_j) / scale
        commonTerm.noalias() = (weights.transpose() * D2) / scale;
        weighted.col(k) = e.weights(k) * values.col(k);
        for (index_t b = 0; b != d; ++b) {
          for (index_t a = 0; a != d; ++a) {
            unit.setZero();
            unit(a) = 1;
            cofactorDerivative(jacGeom, b, unit, dCof);
            dAdjMetric.noalias() = dCof.transpose() * cof;
            dAdjMetric += MatrixD(dAdjMetric.transpose());
            packSymmetric(dAdjMetric, dWeights);
            const T dScale = dAdjMetric.trace();
            for (index_t c = 0; c != d; ++c) {
              value(c, a) = (dWeights.dot(deriv2Geom[c]) - Lx(c) * dScale) / scale;
            }
          }
          for (index_t c = 0; c != d; ++c) {
            if (0 != mode && c != b) { continue; }
            dL[c * d + b].row(k).noalias() = value.row(c) * dB;
            if (c == b) { dL[c * d + b].row(k) += commonTerm; }
          }
        }
      }
      local.setZero(d * A, d * A);
      for (index_t c = 0; c != d; ++c) {
        for (index_t b = 0; b != d; ++b) {
          if (0 != mode && c != b) { continue; }
          if (1 == mode) {
            local.block(b * A, b * A, A, A).diagonal() =
                (weighted.array() * dL[c * d + b].transpose().array())
                    .rowwise().sum();
          } else {
            local.block(c * A, b * A, A, A).noalias() = weighted * dL[c * d + b];
          }
        }
      }
      scatterMatrix(e, local, false, jacobian);
    });
  }

  /// (For tests) number of elements with free dofs, and for element \a i:
  /// its patch, quadrature nodes, active functions and the values, first
  /// and second derivatives of these at the nodes, as used for the assembly
  index_t numActiveElements() const { return m_elements.size(); }
  void elementBasisData(const index_t i, size_t &patch, gsMatrix<T> &nodes,
                        gsMatrix<index_t> &actives, gsMatrix<T> &values,
                        gsMatrix<T> &deriv, gsMatrix<T> &deriv2,
                        bool &tabulated) const {
    const Element &e = m_elements[i];
    patch = e.patch;
    actives = e.actives;
    tabulated = (e.table[0] >= 0);
    allBasisValues(e, values, deriv, deriv2);
    if (!tabulated) { nodes = e.nodes; return; }
    std::vector<gsVector<T> > grid(d);
    for (short_t a = 0; a != d; ++a) {
      grid[a] = m_tables[e.table[a]].nodes;
    }
    gsPointGrid(grid, nodes);
  }

  /// Number of free dofs
  index_t numDofs() const { return m_numDofs; }

 private:
  struct Element {
    index_t id = 0;     // number of the element (patch by patch)
    size_t patch = 0;
    gsVector<T> lower, upper;   // corners of the element
    gsMatrix<T> nodes;          // quadrature nodes (only without tables)
    gsVector<T> weights;
    gsMatrix<index_t> actives;
    std::vector<index_t> dofs;  // global index of (coordinate c, function i), -1: fixed
    std::vector<char> monitor;  // det(J) monitored at node k?
    gsMatrix<T> ders;           // cached derivatives (only without tables)
    std::array<index_t, d> table; // univariate tables per direction, or -1
    index_t numNodes() const { return weights.size(); }
  };

  /// Sparsity pattern (compressed column storage), without values
  struct Pattern {
    std::vector<index_t> outer, inner;
  };

  /// Univariate B-splines of one knot span (degree + 1 functions) and their
  /// first and second derivatives at the univariate quadrature nodes
  struct Table1D {
    index_t degree = 0;
    gsVector<T> nodes;
    std::vector<T> values;  // (order, function, node)
    T &value(index_t order, index_t j, index_t q) {
      return values[(order * (degree + 1) + j) * nodes.size() + q];
    }
    T value(index_t order, index_t j, index_t q) const {
      return values[(order * (degree + 1) + j) * nodes.size() + q];
    }
  };

  /// Orientation (+1 or -1) of the patch of element \a e
  T orientation(const Element &e) const { return m_orientations[e.patch]; }

  static T chiPower(const T chi) {
    return (2 == d ? chi : math::pow(chi, T(2) / d));
  }

  /// First derivatives of the active functions at the quadrature nodes of
  /// \a e, rows i * d + a (function i, direction a), one column per node
  const gsMatrix<T> &derivatives(const Element &e, gsMatrix<T> &buffer) const {
    if (e.table[0] >= 0) {
      tensorBasisValues(e, nullptr, &buffer, nullptr);
      return buffer;
    }
    if (e.ders.size()) { return e.ders; }
    elementDerivatives(m_mb.basis(e.patch), e.nodes, buffer);
    return buffer;
  }

  /// Values (rows i), first derivatives (rows i * d + a) and second
  /// derivatives (rows i * d(d+1)/2 + m, m as in packSymmetric) of the active
  /// functions at the quadrature nodes of \a e
  void allBasisValues(const Element &e, gsMatrix<T> &values,
                      gsMatrix<T> &deriv, gsMatrix<T> &deriv2) const {
    if (e.table[0] >= 0) {
      tensorBasisValues(e, &values, &deriv, &deriv2);
      return;
    }
    // sameElement: the functions of the element (that of the first node,
    // its lower corner) at all nodes, also on its upper knots
    std::vector<gsMatrix<T> > all;
    m_mb.basis(e.patch).evalAllDers_into(e.nodes, 2, all, true);
    values = all[0];
    deriv = all[1];
    deriv2 = all[2];
  }

  /// Univariate tables of the element (tensor-product B-spline basis): the
  /// quadrature rule must be a tensor rule and the active functions must be
  /// numbered lexicographically, otherwise the element keeps the generic
  /// evaluation (e.table[0] < 0).
  void setupTables(const gsTensorBSplineBasis<d, T> &basis,
                   std::map<std::array<index_t, 3>, index_t> &tableOf,
                   Element &e) {
    const gsVector<index_t, d> first = basis.tensorIndex(e.actives(0, 0));
    gsVector<index_t, d> size, numNodes, stride;
    for (short_t a = 0; a != d; ++a) {
      size(a) = basis.component(a).degree() + 1;
    }
    if (size.prod() != e.actives.rows()) { return; }
    // lexicographic numbering of the active functions
    gsVector<index_t, d> v = gsVector<index_t, d>::Zero(), idx;
    index_t r = 0;
    do {
      idx = first + v;
      if (basis.index(idx) != e.actives(r++, 0)) { return; }
    } while (nextLexicographic(v, size));
    // tensor grid of quadrature nodes, first direction fastest
    index_t total = 1;
    for (short_t a = 0; a != d; ++a) {
      stride(a) = total;
      // nodes n * stride(a), n = 0, 1, ...: direction a varies, the others
      // stay at their first coordinate
      index_t n = 1;
      while (n * total < e.numNodes()) {
        bool same = true;
        for (short_t b = 0; b != d; ++b) {
          if (b != a && e.nodes(b, n * total) != e.nodes(b, 0)) { same = false; }
        }
        if (!same) { break; }
        ++n;
      }
      numNodes(a) = n;
      total *= n;
    }
    if (total != e.numNodes()) { return; }
    v.setZero();
    index_t k = 0;
    do {
      for (short_t a = 0; a != d; ++a) {
        if (e.nodes(a, k) != e.nodes(a, v(a) * stride(a))) { return; }
      }
      ++k;
    } while (nextLexicographic(v, numNodes));

    std::array<index_t, d> table;
    gsMatrix<T> ders;
    for (short_t a = 0; a != d; ++a) {
      const std::array<index_t, 3> key = {{index_t(a), first(a), numNodes(a)}};
      auto it = tableOf.find(key);
      if (it == tableOf.end()) {
        const gsBSplineBasis<T> &comp = basis.component(a);
        const index_t p = comp.degree();
        Table1D t;
        t.degree = p;
        t.nodes.resize(numNodes(a));
        t.values.assign(3 * (p + 1) * numNodes(a), T(0));
        for (index_t q = 0; q != numNodes(a); ++q) {
          t.nodes(q) = e.nodes(a, q * stride(a));
          bsplineSpanDerivatives(comp.knots(), p, first(a) + p, t.nodes(q), 2, ders);
          for (index_t o = 0; o != 3; ++o)
            for (index_t j = 0; j <= p; ++j) { t.value(o, j, q) = ders(o, j); }
        }
        it = tableOf.emplace(key, index_t(m_tables.size())).first;
        m_tables.push_back(give(t));
      }
      table[a] = it->second;
    }
    e.table = table;
  }

  /// Values / first / second derivatives (any of them may be skipped) of the
  /// active functions at the nodes of \a e, from the univariate tables; the
  /// products are formed in the same order as in gsTensorBasis
  void tensorBasisValues(const Element &e, gsMatrix<T> *values,
                         gsMatrix<T> *deriv, gsMatrix<T> *deriv2) const {
    const Table1D *tab[d];
    gsVector<index_t, d> size, numNodes;
    for (short_t a = 0; a != d; ++a) {
      tab[a] = &m_tables[e.table[a]];
      size(a) = tab[a]->degree + 1;
      numNodes(a) = tab[a]->nodes.size();
    }
    const index_t A = size.prod(), nq = numNodes.prod();
    const index_t nd2 = d * (d + 1) / 2;
    if (values) { values->resize(A, nq); }
    if (deriv) { deriv->resize(d * A, nq); }
    if (deriv2) { deriv2->resize(nd2 * A, nq); }
    T val[d], der[d], der2[d];
    gsVector<index_t, d> kq = gsVector<index_t, d>::Zero(), v;
    index_t k = 0;
    do {
      v.setZero();
      index_t r = 0;
      do {
        for (short_t a = 0; a != d; ++a) {
          val[a] = tab[a]->value(0, v(a), kq(a));
          der[a] = tab[a]->value(1, v(a), kq(a));
          der2[a] = tab[a]->value(2, v(a), kq(a));
        }
        if (values) {
          T prod = val[0];
          for (short_t i = 1; i < d; ++i) { prod *= val[i]; }
          (*values)(r, k) = prod;
        }
        if (deriv) {
          for (short_t c = 0; c != d; ++c) {
            T prod = der[c];
            for (short_t i = 0; i != d; ++i) {
              if (i != c) { prod *= val[i]; }
            }
            (*deriv)(r * d + c, k) = prod;
          }
        }
        if (deriv2) {
          index_t m = d;
          for (short_t c = 0; c != d; ++c) {
            T prod = der2[c];
            for (short_t i = 0; i != d; ++i) {
              if (i != c) { prod *= val[i]; }
            }
            (*deriv2)(r * nd2 + c, k) = prod;
            for (short_t l = c + 1; l < d; ++l) {
              T mixed = der[c] * der[l];
              for (short_t i = 0; i != d; ++i) {
                if (i != c && i != l) { mixed *= val[i]; }
              }
              (*deriv2)(r * nd2 + m++, k) = mixed;
            }
          }
        }
        ++r;
      } while (nextLexicographic(v, size));
      ++k;
    } while (nextLexicographic(kq, numNodes));
  }


  void localCoefs(const Element &e, gsMatrix<T> &coefs) const {
    const gsMatrix<T> &all = m_mp.patch(e.patch).coefs();
    coefs.resize(e.actives.rows(), d);
    for (index_t i = 0; i != e.actives.rows(); ++i) {
      coefs.row(i) = all.row(e.actives(i, 0));
    }
  }

  /// First derivatives (as deriv_into) of the functions active on the
  /// element whose quadrature nodes are \a nodes (first node: its lower
  /// corner), also at nodes on its upper knots, where deriv_into would take
  /// the functions of the next knot span
  static void elementDerivatives(const gsBasis<T> &basis,
                                 const gsMatrix<T> &nodes, gsMatrix<T> &ders) {
    std::vector<gsMatrix<T> > all;
    basis.evalAllDers_into(nodes, 1, all, true);
    ders.swap(all[1]);
  }

  /// Is det(J) at the node monitored? If the gradients of the free active
  /// functions span all directions (e.g. in the interior of the parameter
  /// domain), det(J) depends on the free dofs unless J is degenerate there,
  /// which must not be overlooked: monitored. Otherwise (on the boundary of
  /// the parameter domain) the dependence is decided by the cofactors
  /// (controllable), which then only involve the fixed control points.
  static bool monitored(const Element &e, const MatrixD &cof,
                        const gsAsConstMatrix<T> &dB) {
    const index_t A = dB.cols();
    MatrixD gram = MatrixD::Zero();
    for (index_t i = 0; i != A; ++i) {
      bool free = false;
      for (index_t c = 0; c != d; ++c) { free = free || e.dofs[c * A + i] >= 0; }
      if (free) { gram.noalias() += dB.col(i) * dB.col(i).transpose(); }
    }
    const gsEigen::SelfAdjointEigenSolver<MatrixD> eig(gram, gsEigen::EigenvaluesOnly);
    const T largest = eig.eigenvalues().maxCoeff();
    if (largest > 0 && eig.eigenvalues().minCoeff() > T(detDependenceTolerance) * largest) {
      return true;
    }
    return controllable(e, cof, dB);
  }

  /// Same criterion as detIsControllable
  static bool controllable(const Element &e, const MatrixD &cof,
                           const gsAsConstMatrix<T> &dB) {
    const index_t A = dB.cols();
    const T tol = T(detDependenceTolerance) * cof.norm()
                * dB.colwise().norm().maxCoeff();
    for (index_t c = 0; c != d; ++c) {
      for (index_t i = 0; i != A; ++i) {
        if (e.dofs[c * A + i] >= 0
            && math::abs(cof.row(c).dot(dB.col(i))) > tol) { return true; }
      }
    }
    return false;
  }

  /// Number of threads of the parallel assembly
  static int numThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  }

  /// Calls body(element, thread) for all elements. With several (OpenMP)
  /// threads, the elements are processed color by color in parallel: the
  /// elements of one color share no free dofs, so that they can add to the
  /// global vectors and matrices concurrently (\a colored = false: any order,
  /// for reductions). With one thread, the elements are processed in their
  /// order (results bit-identical to a serial build).
  /// \a Work: buffers of a thread, constructed on its stack and passed to
  /// every call of \a body as body(element, thread, work).
  template<class Work, class Body>
  void forEachElement(Body &&body, const bool colored = true) const {
#ifdef _OPENMP
    if (omp_get_max_threads() > 1 && m_elements.size() > 1) {
      if (colored && m_colors.empty()) { colorElements(); }
      const index_t numElements = m_elements.size();
#pragma omp parallel
      {
        const int thread = omp_get_thread_num();
        Work work;
        if (colored) {
          for (const std::vector<index_t> &color : m_colors) {
            const index_t size = color.size();
#pragma omp for schedule(dynamic, 4)
            for (index_t i = 0; i < size; ++i) {
              body(m_elements[color[i]], thread, work);
            }
          }
        } else {
#pragma omp for schedule(dynamic, 8)
          for (index_t i = 0; i < numElements; ++i) {
            body(m_elements[i], thread, work);
          }
        }
      }
      return;
    }
#endif
    GISMO_UNUSED(colored);
    Work work;
    for (const Element &e : m_elements) { body(e, 0, work); }
  }

  /// Greedy coloring of the elements: elements with a common free dof get
  /// different colors
  void colorElements() const {
    std::vector<std::vector<index_t> > colorsOfDof(m_numDofs);
    std::vector<char> used;
    m_colors.clear();
    for (size_t i = 0; i != m_elements.size(); ++i) {
      used.assign(m_colors.size() + 1, 0);
      for (index_t g : m_elements[i].dofs) {
        if (g < 0) { continue; }
        for (index_t c : colorsOfDof[g]) { used[c] = 1; }
      }
      index_t color = 0;
      while (used[color]) { ++color; }
      if (color == index_t(m_colors.size())) { m_colors.emplace_back(); }
      m_colors[color].push_back(i);
      for (index_t g : m_elements[i].dofs) {
        if (g >= 0 && (colorsOfDof[g].empty() || colorsOfDof[g].back() != color)) {
          colorsOfDof[g].push_back(color);
        }
      }
    }
  }

  /// Upper bound of the size of the intermediate tensors of the sum
  /// factorization: d * prod_j max(degree_j + 1, nodes_j)
  index_t contractionBound(const Element &e) const {
    index_t bound = d;
    for (short_t j = 0; j != d; ++j) {
      const Table1D &t = m_tables[e.table[j]];
      bound *= std::max<index_t>(t.degree + 1, t.nodes.size());
    }
    return bound;
  }

  /// Quadrature nodes of \a e at which det(J) depends on the free control
  /// points (determined once: on the boundary of the parameter domain this
  /// depends on the fixed control points only)
  static void monitoredNodes(const Element &e, std::vector<index_t> &nodes) {
    nodes.clear();
    for (index_t k = 0; k != e.numNodes(); ++k) {
      if (e.monitor[k]) { nodes.push_back(k); }
    }
  }

  /// Jacobian matrices of the map at the quadrature nodes of \a e:
  /// jacobians(c, k * d + a) = d x_c / d xi_a at node k. For tensor-product
  /// bases by sum factorization: the local control points C(i_1..i_d, c) are
  /// contracted with the univariate tables direction by direction.
  void nodeJacobians(const Element &e, const gsMatrix<T> &coefs,
                     gsMatrix<T> &jacobians, gsMatrix<T> &work) const {
    const index_t nq = e.numNodes();
    jacobians.resize(d, d * nq);
    if (e.table[0] < 0) {
      const gsMatrix<T> &ders = derivatives(e, work);
      const index_t A = e.actives.rows();
      for (index_t k = 0; k != nq; ++k) {
        jacobians.middleCols(k * d, d).noalias() =
            coefs.transpose() * ders.reshapeCol(k, d, A).transpose();
      }
      return;
    }
    typedef gsEigen::Matrix<T, Dynamic, Dynamic> Dense;
    typedef gsEigen::Map<const Dense, 0, gsEigen::OuterStride<> > TableMap;
    // two buffers for the intermediate tensors (no allocation per element)
    const index_t bound = contractionBound(e);
    if (work.size() < 2 * bound) { work.resize(2 * bound, 1); }
    T *buffer[2] = {work.data(), work.data() + bound};
    for (short_t a = 0; a != d; ++a) {
      // X: (s_j x rest) with i_j fastest; after direction j the node index
      // q_j becomes the slowest index: Y(rest, q_j) = sum_i X(i, rest) G(i, q_j)
      const T *X = coefs.data();
      index_t size = coefs.size();
      for (short_t j = 0; j != d; ++j) {
        const Table1D &t = m_tables[e.table[j]];
        const index_t s = t.degree + 1, n = t.nodes.size();
        const index_t rest = size / s;
        // G^T as (n x s) column major: (q, i) -> values[(o * s + i) * n + q]
        const TableMap Gt(&t.values[(j == a ? 1 : 0) * s * n], n, s,
                          gsEigen::OuterStride<>(n));
        gsEigen::Map<Dense> Y(buffer[j % 2], rest, n);
        Y.noalias() = gsEigen::Map<const Dense>(X, s, rest).transpose()
            * Gt.transpose();
        X = Y.data();
        size = rest * n;
      }
      // X: (c, q_1, ..., q_d), c fastest
      const gsEigen::Map<const Dense> flat(X, d, nq);
      for (index_t k = 0; k != nq; ++k) { jacobians.col(k * d + a) = flat.col(k); }
    }
  }

  /// Projection onto the gradients of the active functions:
  /// local(i, c) = sum_k sum_a d phi_i / d xi_a (node k) * fluxes(c, k * d + a)
  /// (the adjoint of nodeJacobians)
  void testGradients(const Element &e, const gsMatrix<T> &fluxes,
                     gsMatrix<T> &local, gsMatrix<T> &work) const {
    const index_t nq = e.numNodes(), A = e.actives.rows();
    local.setZero(A, d);
    if (e.table[0] < 0) {
      const gsMatrix<T> &ders = derivatives(e, work);
      for (index_t k = 0; k != nq; ++k) {
        local.noalias() += ders.reshapeCol(k, d, A).transpose()
            * fluxes.middleCols(k * d, d).transpose();
      }
      return;
    }
    typedef gsEigen::Matrix<T, Dynamic, Dynamic> Dense;
    typedef gsEigen::Map<const Dense, 0, gsEigen::OuterStride<> > TableMap;
    const index_t bound = contractionBound(e);
    if (work.size() < 2 * bound) { work.resize(2 * bound, 1); }
    T *buffer[2] = {work.data(), work.data() + bound};
    for (short_t a = 0; a != d; ++a) {
      // Y: (c, q_1, ..., q_d) = fluxes(c, k * d + a)
      gsEigen::Map<Dense> start(buffer[0], d, nq);
      for (index_t k = 0; k != nq; ++k) { start.col(k) = fluxes.col(k * d + a); }
      const T *Y = buffer[0];
      index_t size = d * nq;
      for (short_t j = d - 1; j >= 0; --j) {
        // Y viewed as (rest x n_j) with q_j slowest; X(i, rest) = sum_q G(i, q) Y(rest, q)
        const Table1D &t = m_tables[e.table[j]];
        const index_t s = t.degree + 1, n = t.nodes.size();
        const index_t rest = size / n;
        const TableMap Gt(&t.values[(j == a ? 1 : 0) * s * n], n, s,
                          gsEigen::OuterStride<>(n));
        gsEigen::Map<Dense> X(buffer[(d - j) % 2], s, rest);
        X.noalias() = Gt.transpose()
            * gsEigen::Map<const Dense>(Y, rest, n).transpose();
        Y = X.data();
        size = s * rest;
      }
      // Y: (i_1, ..., i_d, c), i_1 fastest
      local += gsEigen::Map<const Dense>(Y, A, d);
    }
  }

  /// Adds the element vector \a local (a vector or a map) to \a global
  template<class LocalVector>
  void scatterVector(const Element &e, const LocalVector &local,
                     gsVector<T> &global) const {
    for (size_t l = 0; l != e.dofs.size(); ++l) {
      if (e.dofs[l] >= 0) { global(e.dofs[l]) += local(l); }
    }
  }

  void scatterMatrix(const Element &e, const gsMatrix<T> &local,
                     const bool lowerOnly, gsSparseMatrix<T> &global) const {
    const index_t n = e.dofs.size();
    const index_t *outer = global.outerIndexPtr();
    const index_t *inner = global.innerIndexPtr();
    T *values = global.valuePtr();
    for (index_t jl = 0; jl != n; ++jl) {
      const index_t gj = e.dofs[jl];
      if (gj < 0) { continue; }
      const index_t *first = inner + outer[gj], *last = inner + outer[gj + 1];
      for (index_t il = 0; il != n; ++il) {
        const index_t gi = e.dofs[il];
        if (gi < 0 || (lowerOnly && gi < gj)) { continue; }
        const index_t *pos = std::lower_bound(first, last, gi);
        values[pos - inner] += local(il, jl);
      }
    }
  }

  /// Sparsity pattern (compressed, column major) of the coupling of the
  /// free dofs by the elements (lower triangle only if \a lowerOnly)
  void buildPattern(const bool lowerOnly, Pattern &pattern) const {
    // elements containing each dof
    std::vector<std::vector<index_t> > elementsOf(m_numDofs);
    for (size_t ie = 0; ie != m_elements.size(); ++ie) {
      for (index_t g : m_elements[ie].dofs) {
        if (g >= 0 && (elementsOf[g].empty() || elementsOf[g].back() != index_t(ie))) {
          elementsOf[g].push_back(ie);
        }
      }
    }
    std::vector<index_t> marker(m_numDofs, -1), column;
    std::vector<index_t> outer(m_numDofs + 1, 0), inner;
    for (index_t gj = 0; gj != m_numDofs; ++gj) {
      column.clear();
      for (index_t ie : elementsOf[gj]) {
        for (index_t gi : m_elements[ie].dofs) {
          if (gi < 0 || (lowerOnly && gi < gj) || marker[gi] == gj) { continue; }
          marker[gi] = gj;
          column.push_back(gi);
        }
      }
      std::sort(column.begin(), column.end());
      inner.insert(inner.end(), column.begin(), column.end());
      outer[gj + 1] = inner.size();
    }
    pattern.outer.swap(outer);
    pattern.inner.swap(inner);
  }

  /// A matrix with the given sparsity pattern and zero values
  void makeMatrix(const Pattern &pattern, gsSparseMatrix<T> &matrix) const {
    matrix.resize(m_numDofs, m_numDofs);
    matrix.resizeNonZeros(pattern.inner.size());
    std::copy(pattern.outer.begin(), pattern.outer.end(), matrix.outerIndexPtr());
    std::copy(pattern.inner.begin(), pattern.inner.end(), matrix.innerIndexPtr());
    std::fill(matrix.valuePtr(), matrix.valuePtr() + pattern.inner.size(), T(0));
  }

  gsMultiPatch<T> m_mp;
  gsMultiBasis<T> m_mb;
  const gsDofMapper &m_mapper;
  /// per patch: +1 or -1, the sign of det(J) of a bijective result
  std::vector<T> m_orientations;
  index_t m_numDofs;
  index_t m_h1Power = 2;
  index_t m_energyForm = 0;
  index_t m_numAllElements = 0;
  std::vector<Element> m_elements;
  std::vector<Table1D> m_tables;
  mutable Pattern m_patternLower, m_patternFull;
  mutable std::vector<std::vector<index_t> > m_colors;
};

template<short_t d, typename T>
gsVector<T>
gsBarrierCore<d, T>::computeRegularizedH1(const gsMultiPatch<T> &mp,
                                          const gsDofMapper &mapper,
                                          gsVector<T> x,
                                          const gsOptionList &options) {
  GISMO_ENSURE(d == 2 || d == 3, "The regularized PDE-based method is "
                                 "implemented for dimension 2 and 3 only.");

  const index_t verbose = options.askInt("Verbose", 0);
  const index_t maxRounds = options.getInt("PDE.AdaptiveRounds");
  if (maxRounds <= 0) {
    return regularizedStage(mp, mapper, x, options,
                            std::vector<index_t>(), nullptr);
  }

  // det(J) is checked at the nodes of a denser Gauss-Lobatto rule. Where
  // it is not positive, the rule of the stage is refined in that element and
  // the stage is repeated from the last result.
  const index_t maxExtraNodes = 3;
  std::vector<index_t> extraNodes, negatives;
  gsVector<T> best = x;
  index_t bestNegatives = std::numeric_limits<index_t>::max();
  for (index_t round = 0; ; ++round) {
    const gsVector<T> result = regularizedStage(mp, mapper, x, options,
                                                extraNodes, &negatives);
    index_t numNegative = 0, numElements = 0;
    for (index_t n : negatives) { numNegative += n; numElements += (n > 0); }
    if (verbose > 0) {
      gsInfo << "[PDE] adaptive round " << round << ": " << numNegative
             << " check points with det(J) <= 0 in " << numElements
             << " elements\n";
    }
    if (numNegative < bestNegatives) {
      best = result;
      bestNegatives = numNegative;
    }
    if (0 == numNegative || round == maxRounds) { break; }

    extraNodes.resize(negatives.size(), 0);
    bool refined = false;
    for (size_t e = 0; e != negatives.size(); ++e) {
      if (negatives[e] > 0 && extraNodes[e] < maxExtraNodes) {
        ++extraNodes[e];
        refined = true;
      }
    }
    if (!refined) { break; }
    x = result;
  }
  return best;
}

/**
 * The regularized stage of the PDE-based method for one multi-patch with its
 * dof mapper (see gsBarrierCore::regularizedStage): minimization of
 * Garanzha's regularized Winslow energy (untangling, or to convergence) and
 * the regularized H1 discretization, with the step control that keeps an
 * unfolded parameterization unfolded. The quantities are assembled
 * element-wise (gsRegularizedPDEAssembler) with a Gauss-Lobatto rule, whose
 * nodes include the element boundaries.
 */
template<short_t d, typename T>
class gsRegularizedStage {
 public:
  typedef gsBarrierCore<d, T> Core;
  /// Fixed parameters of the stage
  struct Constants {
    /// every step keeps min det(J) >= stepFraction * (current minimum)
    T stepFraction = T(0.2);
    /// eps = sqrt(epsMin^2 + (epsFactor * min(Dmin, 0))^2), see Garanzha et
    /// al., Foldover-free maps in 50 lines of code, ACM TOG 40(4), 2021
    T epsFactor = T(0.2);
    /// energy minimization to convergence, once unfolded:
    /// eps <= epsShrink * min det(J), so that chi ~ det(J) everywhere
    T epsShrink = T(0.1);
    index_t maxBacktracking = 30;
    T armijo = T(1e-4);
    /// untangling stops once unfolded and the relative energy decrease is
    /// below untangleTol, or if folded and below stagnationTol
    T untangleTol = T(1e-5);
    T stagnationTol = T(1e-12);
    index_t maxUntangleIterations = 200;
    index_t maxFinalIterations = 1000;
    /// forcing term of the inexact Newton steps (CG):
    /// min(maxForcing, sqrt(|g| / |g_0|)), at least minForcing
    T minForcing = T(1e-8), maxForcing = T(0.1);
    index_t maxCGIterations = 500;
    /// tolerance of the H1 stage (residual norm)
    T h1Tolerance = T(1e-5);
  };

  gsRegularizedStage(const gsMultiPatch<T> &mp, const gsDofMapper &mapper,
                     const gsOptionList &options,
                     const std::vector<index_t> &extraNodes)
      : m_mp(mp), m_mapper(mapper), m_options(options),
        m_verbose(options.askInt("Verbose", 0)),
        m_quOptions(lobattoOptions(options.getInt("PDE.ExtraNodes"))),
        m_patchDet(patchDeterminants(mp, mapper, m_quOptions)),
        m_meanDet(meanAbsolute(m_patchDet)),
        m_orientations(signs(m_patchDet)),
        m_epsMin(options.getReal("PDE.GaranzhaEps") * m_meanDet),
        m_fem(mp, mapper, m_quOptions, m_orientations, extraNodes),
        m_iterative(gsBarrierCore<d, T>::useIterativeSolver(m_fem.numDofs(),
                                                            options)) {
    m_fem.setH1Power(Core::H1Form::WeakLaplace == options.getInt("PDE.H1Form")
                         ? 1 : 2);
    index_t form = options.getInt("PDE.EnergyForm");
    if (form != Core::EnergyForm::InverseDirichlet) { form = Core::EnergyForm::Conformal; }
    m_fem.setEnergyForm(form);
  }

  /// False if the computational domain, or one of its patches, has zero
  /// measure
  bool valid() const {
    for (const T patchDet : m_patchDet) {
      if (!(patchDet != 0)) { return false; }
    }
    return m_meanDet > 0;
  }

  /// True if the patches do not all have the same orientation
  bool mixedOrientation() const {
    return std::find(m_orientations.begin(), m_orientations.end(),
                     -m_orientations.front()) != m_orientations.end();
  }

  /// Minimum of orientation * det(J) over the monitored quadrature nodes
  T minDet(const gsVector<T> &x) {
    m_fem.setFreeVector(x);
    return m_fem.minDet();
  }

  /// For every element, the number of nodes of a denser Gauss-Lobatto rule
  /// (PDE.CheckExtraNodes additional nodes) with det(J) <= 0
  std::vector<index_t> negativeCheckPoints(const gsVector<T> &x) {
    m_fem.setFreeVector(x);
    return m_fem.negativeCheckPoints(
        lobattoOptions(m_options.getInt("PDE.CheckExtraNodes")));
  }

  /// Minimization of the regularized energy from \a x (with min det(J)
  /// \a Dmin), by Gauss-Newton or projected Newton steps with Armijo line
  /// search; once unfolded, no step folds again. final = false: untangling
  /// (stops once unfolded and settled); final = true: to convergence
  /// (PDE.EnergyTol), with eps <= epsShrink * min det(J) once unfolded.
  void minimizeEnergy(gsVector<T> &x, T &Dmin, const bool final);

  /// Regularized H1 stage (Anderson acceleration with the Jacobian as
  /// preconditioner and the step control) from the unfolded \a x; returns
  /// true if it converged (then \a x is the result)
  bool solveH1(gsVector<T> &x);

 private:
  /// State of the linear solver of the energy minimization
  struct NewtonSolver {
    typedef gsEigen::SparseMatrix<T, 0, index_t> EigenSparse;
    gsEigen::SimplicialLDLT<EigenSparse, gsEigen::Lower> ldlt;
    gsEigen::ConjugateGradient<EigenSparse, gsEigen::Lower,
        gsEigen::IncompleteCholesky<T, gsEigen::Lower,
                                    gsEigen::AMDOrdering<index_t> > > pcg;
    memory::shared_ptr<preAApp::NodalAMG<T> > amg;
    bool patternAnalyzed = false;
    T gradientNorm0 = -1;
    index_t numCGIterations = 0;
    double timeAssembly = 0, timeSolve = 0;
  };

  /// Gauss-Lobatto rule with 2 + extra nodes per direction beyond the degree
  static gsOptionList lobattoOptions(const index_t extra) {
    gsOptionList result = gsExprAssembler<T>::defaultOptions();
    result.setInt("quRule", gsQuadrature::GaussLobatto);
    result.setReal("quA", 1.0);
    result.setInt("quB", 2 + std::max<index_t>(0, extra));
    return result;
  }

  /// For every patch, the integral of det(J) (its signed area or volume)
  /// divided by the measure of the whole parameter domain. The integral only
  /// depends on the control points on the boundary of the patch; to evaluate
  /// it reliably even if the given interior is far off, the free control
  /// points inside the patch are moved to the centroid of the fixed ones.
  static std::vector<T> patchDeterminants(const gsMultiPatch<T> &mp,
                                          const gsDofMapper &mapper,
                                          const gsOptionList &quOptions);

  /// Mean value of |det(J)|: sum of the absolute values of \a patchDet
  static T meanAbsolute(const std::vector<T> &patchDet) {
    T result = 0;
    for (const T value : patchDet) { result += math::abs(value); }
    return result;
  }

  /// Orientation of every patch: the sign of \a patchDet
  static std::vector<T> signs(const std::vector<T> &patchDet) {
    std::vector<T> result(patchDet.size());
    for (size_t p = 0; p != patchDet.size(); ++p) {
      result[p] = (patchDet[p] < 0 ? T(-1) : T(1));
    }
    return result;
  }

  /// Gradient and Hessian of the energy at the current free vector of the
  /// assembler, and the Newton direction (linearSolve)
  void newtonDirection(const T eps, const index_t hessianType,
                       NewtonSolver &solver, gsVector<T> &gradient,
                       gsSparseMatrix<T> &hessian, gsVector<T> &direction);

  /// Newton direction for the given gradient and (lower triangle of the)
  /// Hessian: sparse LDL^T, or inexact CG preconditioned by NodalAMG
  /// (PDE.Preconditioner = BlockGaussSeidelLU: incomplete Cholesky)
  void linearSolve(NewtonSolver &solver, const gsVector<T> &gradient,
                   const gsSparseMatrix<T> &hessian, gsVector<T> &direction);

  const gsMultiPatch<T> &m_mp;
  const gsDofMapper &m_mapper;
  const gsOptionList &m_options;
  const index_t m_verbose;
  const gsOptionList m_quOptions;
  const Constants m_constants = Constants();
  const std::vector<T> m_patchDet;
  const T m_meanDet;
  const std::vector<T> m_orientations;
  const T m_epsMin;
  gsRegularizedPDEAssembler<d, T> m_fem;
  const bool m_iterative;
};

template<short_t d, typename T>
std::vector<T> gsRegularizedStage<d, T>::patchDeterminants(const gsMultiPatch<T> &mp,
                                                           const gsDofMapper &mapper,
                                                           const gsOptionList &quOptions) {
  gsMultiPatch<T> mpSubstitute = mp;
  gsVector<T> centroid = gsVector<T>::Zero(mp.targetDim());
  gsVector<index_t> numFixed = gsVector<index_t>::Zero(mp.targetDim());
  for (size_t iptch = 0; iptch != mp.nPatches(); ++iptch)
    for (index_t idim = 0; idim != mp.targetDim(); ++idim)
      for (index_t idof = 0; idof != index_t(mapper.patchSize(iptch, idim)); ++idof)
        if (!mapper.is_free(idof, iptch, idim)) {
          centroid(idim) += mp.patch(iptch).coefs()(idof, idim);
          ++numFixed(idim);
        }
  if (numFixed.minCoeff() > 0) {
    for (size_t iptch = 0; iptch != mp.nPatches(); ++iptch) {
      // free control points on the boundary of the patch (free interfaces)
      // are kept: they bound the patch
      const gsMatrix<index_t> boundary = mp.basis(iptch).allBoundary();
      std::vector<bool> onBoundary(mp.patch(iptch).coefs().rows(), false);
      for (index_t i = 0; i != boundary.rows(); ++i) { onBoundary[boundary(i, 0)] = true; }
      for (index_t idim = 0; idim != mp.targetDim(); ++idim)
        for (index_t idof = 0; idof != index_t(mapper.patchSize(iptch, idim)); ++idof)
          if (mapper.is_free(idof, iptch, idim) && !onBoundary[idof])
            mpSubstitute.patch(iptch).coefs()(idof, idim) =
                centroid(idim) / numFixed(idim);
    }
  }

  // measure of the parameter domain
  T parameterMeasure = 0;
  for (size_t iptch = 0; iptch != mp.nPatches(); ++iptch) {
    const gsMatrix<T> box = mp.basis(iptch).support();
    parameterMeasure += (box.col(1) - box.col(0)).prod();
  }

  std::vector<T> result(mp.nPatches());
  for (size_t iptch = 0; iptch != mp.nPatches(); ++iptch) {
    const gsMultiPatch<T> patch(mpSubstitute.patch(iptch));
    gsMultiBasis<T> mb(patch);
    gsExprAssembler<T> assembler;
    assembler.setIntegrationDomain(mb.domain());
    assembler.setOptions(quOptions);
    gsExprEvaluator<T> evaluator(assembler);
    typename gsExprAssembler<T>::geometryMap G = evaluator.getMap(patch);
    result[iptch] = evaluator.integral(jac(G).det()) / parameterMeasure;
  }
  return result;
}

template<short_t d, typename T>
void gsRegularizedStage<d, T>::newtonDirection(const T eps,
                                               const index_t hessianType,
                                               NewtonSolver &solver,
                                               gsVector<T> &gradient,
                                               gsSparseMatrix<T> &hessian,
                                               gsVector<T> &direction) {
  gsStopwatch timer;
  if (Core::EnergyForm::Conformal == m_fem.energyForm()
      && Core::EnergyHessian::GaussNewton == hessianType) {
    m_fem.gradientHessian(eps, gradient, hessian);  // faster variant
  } else {
    m_fem.energyGradientHessian(eps, hessianType, gradient, hessian);
  }
  solver.timeAssembly += timer.stop();
  timer.restart();
  linearSolve(solver, gradient, hessian, direction);
  solver.timeSolve += timer.stop();
}

template<short_t d, typename T>
void gsRegularizedStage<d, T>::linearSolve(NewtonSolver &solver,
                                           const gsVector<T> &gradient,
                                           const gsSparseMatrix<T> &hessian,
                                           gsVector<T> &direction) {
  if (!m_iterative) {
    if (!solver.patternAnalyzed) {
      solver.ldlt.analyzePattern(hessian);
      solver.patternAnalyzed = true;
    }
    solver.ldlt.factorize(hessian);
    direction = -solver.ldlt.solve(gradient);
    return;
  }
  // Inexact Newton: relative tolerance (forcing term). CG starts from zero,
  // which yields a descent direction for any number of CG iterations.
  const T gradientNorm = gradient.norm();
  if (solver.gradientNorm0 < 0) { solver.gradientNorm0 = gradientNorm; }
  const T forcing = std::max<T>(
      m_constants.minForcing,
      std::min<T>(m_constants.maxForcing,
                  math::sqrt(gradientNorm / solver.gradientNorm0)));
  if (Core::Preconditioner::BlockGaussSeidelLU
      == m_options.getInt("PDE.Preconditioner")) {
    // incomplete Cholesky (the preconditioner before NodalAMG)
    const gsVector<T> zero = gsVector<T>::Zero(gradient.size());
    solver.pcg.setTolerance(forcing);
    solver.pcg.compute(hessian);
    direction = -solver.pcg.solveWithGuess(gradient, zero);
    // (Eigen does not count the last iteration if it converged)
    solver.numCGIterations += solver.pcg.iterations()
        + (solver.pcg.iterations() < solver.pcg.maxIterations() ? 1 : 0);
    return;
  }
  if (!solver.amg) {
    solver.amg.reset(new preAApp::NodalAMG<T>());
    solver.amg->setNodes(gsBarrierCore<d, T>::nodeOfDofs(m_mp, m_mapper),
                         gsBarrierCore<d, T>::componentOfDofs(m_mp, m_mapper));
  }
  const gsSparseMatrix<T> full =
      hessian.template selfadjointView<gsEigen::Lower>();
  solver.amg->compute(full);
  gsConjugateGradient<T> cg(full, solver.amg);
  cg.setTolerance(forcing);
  cg.setMaxIterations(m_constants.maxCGIterations);
  gsMatrix<T> solution = gsMatrix<T>::Zero(gradient.size(), 1);
  cg.solve(gradient, solution);
  solver.numCGIterations += cg.iterations();
  direction = -solution.col(0);
}

template<short_t d, typename T>
void gsRegularizedStage<d, T>::minimizeEnergy(gsVector<T> &x, T &Dmin,
                                              const bool final) {
  // Hessian: as set by PDE.EnergyHessian; automatic: projected Newton for the
  // final minimization of volumes, Gauss-Newton otherwise
  index_t hessianType = m_options.getInt("PDE.EnergyHessian");
  if (hessianType < Core::EnergyHessian::GaussNewton
      || hessianType > Core::EnergyHessian::ProjectedNewton) {
    hessianType = (final && 3 == d) ? Core::EnergyHessian::ProjectedNewton
                                    : Core::EnergyHessian::GaussNewton;
  }
  const index_t maxIterations = final ? m_constants.maxFinalIterations
                                      : m_constants.maxUntangleIterations;
  const T tolerance = final ? T(m_options.getReal("PDE.EnergyTol"))
                            : m_constants.untangleTol;
  verboseLog(final ? std::string("\nMinimizing the regularized Winslow "
                                 "energy (final stage)...")
                   : "\nThe starting point is folded (min det(J) / mean det(J) = "
                         + std::to_string(Dmin / m_meanDet)
                         + "). Start untangling...", m_verbose);

  T eps = m_epsMin;
  auto energy = [&](const gsVector<T> &u) {
    m_fem.setFreeVector(u);
    return m_fem.energy(eps);
  };
  NewtonSolver solver;
  gsSparseMatrix<T> hessian;
  gsVector<T> gradient, direction, trial;
  gsStopwatch stageTimer, timer;
  double timeEnergy = 0, timeLineSearch = 0;
  index_t numEnergyEvaluations = 0, numNewtonFallbacks = 0;
  index_t iter = 0;
  for (; iter != maxIterations; ++iter) {
    const T negPart = std::min<T>(Dmin, 0);
    eps = math::sqrt(m_epsMin * m_epsMin + m_constants.epsFactor
                     * m_constants.epsFactor * negPart * negPart);
    // to convergence, once unfolded: chi ~ det(J) everywhere, so that the
    // energy acts as a barrier against degenerate cells (with eps > det(J)
    // it is bounded there)
    if (final && Dmin > 0) { eps = std::min(eps, m_constants.epsShrink * Dmin); }

    timer.restart();
    const T W0 = energy(x);
    ++numEnergyEvaluations;
    timeEnergy += timer.stop();
    m_fem.setFreeVector(x);
    newtonDirection(eps, hessianType, solver, gradient, hessian, direction);
    if (Core::EnergyHessian::GaussNewton != hessianType
        && (!direction.allFinite() || !(gradient.dot(direction) < 0))) {
      // no descent direction from the (projected) Newton step:
      // Gauss-Newton step instead
      ++numNewtonFallbacks;
      newtonDirection(eps, Core::EnergyHessian::GaussNewton, solver, gradient,
                      hessian, direction);
    }

    const T slope = gradient.dot(direction);
    if (!direction.allFinite() || !(slope < 0)) { break; }

    // Armijo line search; once unfolded, never fold again
    timer.restart();
    T step = 1, W1 = W0, D1 = Dmin;
    bool accepted = false;
    for (index_t ib = 0; ib != m_constants.maxBacktracking; ++ib, step *= 0.5) {
      trial = x + step * direction;
      W1 = energy(trial);
      ++numEnergyEvaluations;
      if (!(W1 <= W0 + m_constants.armijo * step * slope)) { continue; }
      D1 = minDet(trial);
      if (Dmin > 0 && !(D1 >= m_constants.stepFraction * Dmin)) { continue; }
      accepted = true;
      break;
    }
    timeLineSearch += timer.stop();
    if (!accepted) { break; }

    const bool wasUnfolded = (Dmin > 0);
    x = trial;
    Dmin = D1;
    if (m_verbose > 1) {
      gsInfo << (final ? "  energy iter " : "  untangling iter ") << iter
             << ": energy = " << W1
             << ", min det(J) / mean det(J) = " << Dmin / m_meanDet << "\n";
    }
    // unfolded and the energy has settled
    if (wasUnfolded && Dmin > 0 && (W0 - W1) <= tolerance * W0) {
      ++iter;
      break;
    }
    // still folded, but no progress any more: give up
    if (!(Dmin > 0) && (W0 - W1) <= m_constants.stagnationTol * W0) {
      ++iter;
      break;
    }
  }
  verboseLog((final ? "Energy minimization finished after " : "Untangling finished after ")
                 + std::to_string(iter)
                 + " iterations, min det(J) / mean det(J) = "
                 + std::to_string(Dmin / m_meanDet), m_verbose);
  if (m_verbose > 0) {
    gsInfo << (final ? "[PDE] energy minimization: " : "[PDE] untangling: ")
           << iter << " iterations, " << numEnergyEvaluations
           << " energy evaluations";
    if (m_iterative) {
      gsInfo << ", " << solver.numCGIterations << " CG iterations";
    }
    gsInfo << "; time " << stageTimer.stop() << " s (assembly "
           << timeEnergy + solver.timeAssembly << " s, linear solver "
           << solver.timeSolve << " s, line search " << timeLineSearch
           << " s); energy form " << m_fem.energyForm() << ", Hessian "
           << hessianType << ", " << numNewtonFallbacks
           << " Gauss-Newton fall-backs\n";
  }
}

template<short_t d, typename T>
bool gsRegularizedStage<d, T>::solveH1(gsVector<T> &x) {
  verboseLog("\nStart parameterization improvement by the regularized H1 "
             "discretization...", m_verbose);
  typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;
  typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;

  Residual_t residual = [&](gsVector<T> const &u) {
    gsVector<T> result;
    m_fem.setFreeVector(u);
    m_fem.h1Residual(m_epsMin, result);
    return result;
  };
  // AAPreconditionType 1: diagonal blocks of the Jacobian, otherwise the full
  // analytical Jacobian
  const bool diagonalBlocks = (1 == m_options.askInt("AAPreconditionType", 0));
  Jacobian_t jacobian = [&](gsVector<T> const &u) {
    gsSparseMatrix<T> result;
    m_fem.setFreeVector(u);
    m_fem.h1Jacobian(m_epsMin, diagonalBlocks, result);
    return result;
  };

  preAAParam<T> param;
  param.m = m_options.askInt("AAwindowsize", 5);
  param.usePreconditioning = true;
  param.updatePreconditionerStep = m_options.askInt("N_update", 5);
  param.epsilon = m_constants.h1Tolerance;
  if (m_iterative) {
    param.linearSolver =
        gsBarrierCore<d, T>::gmresPreconditioner(m_options, false);
    param.blocks = gsBarrierCore<d, T>::componentOfDofs(m_mp, m_mapper);
    param.nodes = gsBarrierCore<d, T>::nodeOfDofs(m_mp, m_mapper);
    param.innerMaxIterations = m_options.getInt("PDE.InnerMaxIterations");
  }
  param.safeguardFactor = m_options.getReal("PDE.H1Safeguard");

  preAApp::AndersonAcceleration<T> solver(param);
  if (m_verbose > 1) { solver.enableIterationInfoPrinting(); }
  else { solver.disableIterationInfoPrinting(); }

  // Backtracking: an unfolded parameterization stays unfolded
  index_t numShortenedSteps = 0;
  double timeFilter = 0;
  // the last accepted iterate and its min det(J): usually the next \a from
  gsVector<T> accepted;
  T acceptedMinDet = 0;
  solver.setStepFilter([&](gsVector<T> const &from, gsVector<T> &to) {
    gsStopwatch filterTimer;
    const bool known = accepted.size() == from.size()
                       && (accepted.array() == from.array()).all();
    const T bound = m_constants.stepFraction * (known ? acceptedMinDet : minDet(from));
    const gsVector<T> direction = to - from;
    T step = 1;
    for (index_t ib = 0; ib != m_constants.maxBacktracking; ++ib, step *= 0.5) {
      to = from + step * direction;
      const T toMinDet = minDet(to);
      if (toMinDet >= bound) {
        if (ib > 0) { ++numShortenedSteps; }
        accepted = to;
        acceptedMinDet = toMinDet;
        timeFilter += filterTimer.stop();
        return true;
      }
    }
    timeFilter += filterTimer.stop();
    return false;
  });

  gsStopwatch stageTimer;
  const gsVector<T> result = solver.compute(x, residual, jacobian);
  if (m_verbose > 0) {
    gsInfo << aaSummary("regularized H1", solver, stageTimer.stop())
           << "; " << numShortenedSteps << " shortened steps, step control "
           << timeFilter << " s\n";
  }
  if (!solver.converged() || !result.allFinite()) {
    verboseLog("The regularized H1 solver did not converge (residual norm "
                   + std::to_string(solver.residualNorm()) + ").", m_verbose);
    return false;
  }
  x = result;
  return true;
}

template<short_t d, typename T>
gsVector<T>
gsBarrierCore<d, T>::regularizedStage(const gsMultiPatch<T> &mp,
                                      const gsDofMapper &mapper,
                                      gsVector<T> x,
                                      const gsOptionList &options,
                                      const std::vector<index_t> &extraNodes,
                                      std::vector<index_t> *negatives) {
  if (negatives) { negatives->clear(); }
  gsRegularizedStage<d, T> stage(mp, mapper, options, extraNodes);
  if (!stage.valid()) {
    gsWarn << "The computational domain has zero measure; the regularized "
              "PDE-based method is skipped.\n";
    return x;
  }
  auto finish = [&](const gsVector<T> &result) {
    if (negatives) { *negatives = stage.negativeCheckPoints(result); }
    return result;
  };

  // 1. Energy minimization: untangling if folded. With PDE.EnergyFinal (set for
  //    free interfaces) it is the last stage: it also runs for an unfolded
  //    start, to convergence, and the H1 stage is skipped.
  const bool energyFinal = options.getSwitch("PDE.EnergyFinal");
  T Dmin = stage.minDet(x);
  if (stage.mixedOrientation()) {
    verboseLog("The patches are oriented differently; det(J) is taken with "
               "the orientation of each patch.", options.askInt("Verbose", 0));
  }
  if (!(Dmin > 0) || energyFinal) { stage.minimizeEnergy(x, Dmin, energyFinal); }
  if (!(Dmin > 0)) {
    gsWarn << "The parameterization could not be untangled: det(J) <= 0 at "
              "some quadrature nodes"
           << (energyFinal ? "" : "; the H1 discretization is skipped") << ".\n";
    return finish(x);
  }
  if (energyFinal || !options.askSwitch("needPDEH1", true)) { return finish(x); }

  // 2. Regularized H1 stage. If it does not converge: back to its start,
  //    which is only unfolded at the quadrature nodes and may be barely so;
  //    the energy is then minimized to convergence from there
  //    (PDE.EnergyFallback).
  if (!stage.solveH1(x)) {
    if (options.getSwitch("PDE.EnergyFallback")) {
      verboseLog("The energy is minimized instead.", options.askInt("Verbose", 0));
      Dmin = stage.minDet(x);
      stage.minimizeEnergy(x, Dmin, true);
    } else {
      verboseLog("Its starting point is returned.", options.askInt("Verbose", 0));
    }
  }
  return finish(x);
}


namespace expr {

template<typename E1, typename E2>
class frprod2_expr;
template<typename E1, typename E2>
class frprod3_expr;
template<class E>
class jacScaledLx_expr;
template<class E>
class jacScaledLxDiag_expr;
template<class E>
class jacScaledLxDiagBlock_expr;
template<class E>
class jacScaledLxH1_expr;
template<class E>
class jacScaledLxH1DiagBlock_expr;

/**

  Products of the test functions of the space u with the vector
  expression M: one entry per active function phi_i and component c,
  phi_i * M_c, i.e. the element vector of the weak form int phi M dxi
  (column vector, d*n x 1). Used for the residual of the H2
  discretization, M = L(x).

 */
template<typename E1, typename E2>
class frprod3_expr : public _expr<frprod3_expr<E1, E2> > {
 public:
  typedef typename E2::Scalar Scalar;
  enum { ScalarValued = 0, Space = E1::Space, ColBlocks = 0 };

 private:
  typename E1::Nested_t _u;
  typename E2::Nested_t _v;

  mutable gsMatrix<Scalar> res, bGrads, b;

 public:

  frprod3_expr(_expr<E1> const &u, _expr<E2> const &v)
      : _u(u), _v(v) {}

  const gsMatrix<Scalar> &eval(const index_t k) const //todo: specialize for nb==1
  {
    auto A = _v.eval(k);
    b = _u.eval(k);

    res.noalias() = b * A;
    return res;
  }

  index_t rows() const { return 1; }
  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    _v.parse(evList);
    _u.parse(evList);
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _v.rowVar(); }

  void print(std::ostream &os) const {
    os << "(";
    _u.print(os);
    os << " % ";
    _v.print(os);
    os << ")";
  }
};

/// Frobenious product (also known as double dot product) operator for expressions
template<typename E1, typename E2>
EIGEN_STRONG_INLINE
frprod3_expr<E1, E2> const frprod3(E1 const &u,
                                   E2 const &M) {
  return frprod3_expr<E1, E2>(u, M);
}

/**
  Expression for Jacobian matrix for PDE-based parameterization construction
*/
template<class E>
class jacScaledLx_expr : public _expr<jacScaledLx_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  jacScaledLx_expr(const E &u, const gsGeometryMap<Scalar> &G) : _u(u), _G(G) {}

  mutable gsMatrix<Scalar> res, derivGeom, deriv2Geom, derivBasis, deriv2Basis;
  mutable gsMatrix<Scalar> dg11dx, dg11dy, dg22dx, dg22dy, dg12dx, dg12dy;
  mutable gsMatrix<Scalar> commonTerm, dLxdx, dLxdy, dLydx, dLydy;

  const gsMatrix<Scalar> &eval(const index_t k) const {
    gsMatrix<Scalar> basis = _u.data().values[0].col(k);

    derivBasis = _u.data().values[1].col(k).transpose();
    deriv2Basis = _u.data().values[2].col(k).transpose();

    derivBasis.blockTransposeInPlace(_u.dim());
    deriv2Basis.blockTransposeInPlace(1 + _u.dim());

    derivGeom = _G.data().values[1].col(k);
    deriv2Geom = _G.data().values[2].col(k);

    Scalar g11 = derivGeom(0) * derivGeom(0) + derivGeom(2) * derivGeom(2);
    Scalar g12 = derivGeom(0) * derivGeom(1) + derivGeom(2) * derivGeom(3);
    Scalar g22 = derivGeom(1) * derivGeom(1) + derivGeom(3) * derivGeom(3);

    Scalar scaleFactor = g11 + g22;

    Scalar Lx =
        (g22 * deriv2Geom(0) + g11 * deriv2Geom(1) - 2.0 * g12 * deriv2Geom(2))
            / scaleFactor;
    Scalar Ly =
        (g22 * deriv2Geom(3) + g11 * deriv2Geom(4) - 2.0 * g12 * deriv2Geom(5))
            / scaleFactor;

    dg11dx.noalias() = 2.0 * derivGeom(0) * derivBasis.row(0);
    dg11dy.noalias() = 2.0 * derivGeom(2) * derivBasis.row(0);
    dg22dx.noalias() = 2.0 * derivGeom(1) * derivBasis.row(1);
    dg22dy.noalias() = 2.0 * derivGeom(3) * derivBasis.row(1);
    dg12dx.noalias() =
        derivGeom(1) * derivBasis.row(0) + derivGeom(0) * derivBasis.row(1);
    dg12dy.noalias() =
        derivGeom(3) * derivBasis.row(0) + derivGeom(2) * derivBasis.row(1);

    commonTerm.noalias() = g22 * deriv2Basis.row(0) + g11 * deriv2Basis.row(1)
        - 2.0 * g12 * deriv2Basis.row(2);
    dLxdx.noalias() = dg22dx * deriv2Geom(0) + dg11dx * deriv2Geom(1)
        - 2.0 * dg12dx * deriv2Geom(2) + commonTerm;
    dLxdy.noalias() = dg22dy * deriv2Geom(0) + dg11dy * deriv2Geom(1)
        - 2.0 * dg12dy * deriv2Geom(2);
    dLydx.noalias() = dg22dx * deriv2Geom(3) + dg11dx * deriv2Geom(4)
        - 2.0 * dg12dx * deriv2Geom(5);
    dLydy.noalias() = dg22dy * deriv2Geom(3) + dg11dy * deriv2Geom(4)
        - 2.0 * dg12dy * deriv2Geom(5) + commonTerm;

    dLxdx = (dLxdx - Lx * (dg11dx + dg22dx)) / scaleFactor;
    dLxdy = (dLxdy - Lx * (dg11dy + dg22dy)) / scaleFactor;
    dLydx = (dLydx - Ly * (dg11dx + dg22dx)) / scaleFactor;
    dLydy = (dLydy - Ly * (dg11dy + dg22dy)) / scaleFactor;

    const index_t A = _u.cardinality() / _u.dim(); // _u.data().actives.rows()
    res.resize(_u.cardinality(), _u.cardinality());
    res.topLeftCorner(A, A).noalias() = basis * dLxdx;
    res.topRightCorner(A, A).noalias() = basis * dLxdy;
    res.bottomLeftCorner(A, A).noalias() = basis * dLydx;
    res.bottomRightCorner(A, A).noalias() = basis * dLydy;

    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_VALUE | NEED_GRAD | NEED_DERIV2;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV | NEED_DERIV2;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  // rowVar / colVar: the spaces of the rows and columns of the assembled
  // matrix; cardinality: the number of local (active) functions of the space
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLx(";
    _u.print(os);
    os << ")";
  }
};

/**
  Expression for Jacobian matrix (diagonal part) for PDE-based parameterization construction
*/
template<class E>
class jacScaledLxDiag_expr : public _expr<jacScaledLxDiag_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  jacScaledLxDiag_expr(const E &u, const gsGeometryMap<Scalar> &G)
      : _u(u), _G(G) {}

  mutable gsMatrix<Scalar> res, derivGeom, deriv2Geom, derivBasis, deriv2Basis;
  mutable gsMatrix<Scalar> dg11dx, dg11dy, dg22dx, dg22dy, dg12dx, dg12dy;
  mutable gsMatrix<Scalar> commonTerm, dLxdx, dLydy;

  const gsMatrix<Scalar> &eval(const index_t k) const {
    derivBasis = _u.data().values[1].col(k).transpose();
    deriv2Basis = _u.data().values[2].col(k).transpose();

    derivBasis.blockTransposeInPlace(_u.dim());
    deriv2Basis.blockTransposeInPlace(1 + _u.dim());

    derivGeom = _G.data().values[1].col(k);
    deriv2Geom = _G.data().values[2].col(k);

    Scalar g11 = derivGeom(0) * derivGeom(0) + derivGeom(2) * derivGeom(2);
    Scalar g12 = derivGeom(0) * derivGeom(1) + derivGeom(2) * derivGeom(3);
    Scalar g22 = derivGeom(1) * derivGeom(1) + derivGeom(3) * derivGeom(3);

    Scalar scaleFactor = g11 + g22;

    Scalar Lx =
        (g22 * deriv2Geom(0) + g11 * deriv2Geom(1) - 2.0 * g12 * deriv2Geom(2))
            / scaleFactor;
    Scalar Ly =
        (g22 * deriv2Geom(3) + g11 * deriv2Geom(4) - 2.0 * g12 * deriv2Geom(5))
            / scaleFactor;

    dg11dx.noalias() = 2.0 * derivGeom(0) * derivBasis.row(0);
    dg11dy.noalias() = 2.0 * derivGeom(2) * derivBasis.row(0);
    dg22dx.noalias() = 2.0 * derivGeom(1) * derivBasis.row(1);
    dg22dy.noalias() = 2.0 * derivGeom(3) * derivBasis.row(1);
    dg12dx.noalias() =
        derivGeom(1) * derivBasis.row(0) + derivGeom(0) * derivBasis.row(1);
    dg12dy.noalias() =
        derivGeom(3) * derivBasis.row(0) + derivGeom(2) * derivBasis.row(1);

    commonTerm.noalias() = g22 * deriv2Basis.row(0) + g11 * deriv2Basis.row(1)
        - 2.0 * g12 * deriv2Basis.row(2);
    dLxdx.noalias() = dg22dx * deriv2Geom(0) + dg11dx * deriv2Geom(1)
        - 2.0 * dg12dx * deriv2Geom(2) + commonTerm;
    dLydy.noalias() = dg22dy * deriv2Geom(3) + dg11dy * deriv2Geom(4)
        - 2.0 * dg12dy * deriv2Geom(5) + commonTerm;

    dLxdx = (dLxdx - Lx * (dg11dx + dg22dx)) / scaleFactor;
    dLydy = (dLydy - Ly * (dg11dy + dg22dy)) / scaleFactor;

    const index_t A = _u.cardinality() / _u.dim(); // _u.data().actives.rows()
    res.resize(_u.cardinality(), _u.cardinality());
    res.setZero();
    res.topLeftCorner(A, A) = (_u.data().values[0].col(k).array()
        * dLxdx.array()).matrix().asDiagonal();
    res.bottomRightCorner(A, A) = (_u.data().values[0].col(k).array()
        * dLydy.array()).matrix().asDiagonal();

    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_VALUE | NEED_GRAD | NEED_DERIV2;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV | NEED_DERIV2;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  // rowVar / colVar: the spaces of the rows and columns of the assembled
  // matrix; cardinality: the number of local (active) functions of the space
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLxDiag(";
    _u.print(os);
    os << ")";
  }
};

/**
  Expression for Jacobian matrix (diagonal block) for PDE-based parameterization construction
*/
template<class E>
class jacScaledLxDiagBlock_expr : public _expr<jacScaledLxDiagBlock_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  jacScaledLxDiagBlock_expr(const E &u, const gsGeometryMap<Scalar> &G)
      : _u(u), _G(G) {}

  mutable gsMatrix<Scalar> res, basis, derivGeom, deriv2Geom, derivBasis,
      deriv2Basis;
  mutable gsMatrix<Scalar> dg11dx, dg11dy, dg22dx, dg22dy, dg12dx, dg12dy;
  mutable gsMatrix<Scalar> commonTerm, dLxdx, dLydy;

  const gsMatrix<Scalar> &eval(const index_t k) const {
    derivBasis = _u.data().values[1].col(k).transpose();
    deriv2Basis = _u.data().values[2].col(k).transpose();

    derivBasis.blockTransposeInPlace(_u.dim());
    deriv2Basis.blockTransposeInPlace(1 + _u.dim());

    derivGeom = _G.data().values[1].col(k);
    deriv2Geom = _G.data().values[2].col(k);

    Scalar g11 = derivGeom(0) * derivGeom(0) + derivGeom(2) * derivGeom(2);
    Scalar g12 = derivGeom(0) * derivGeom(1) + derivGeom(2) * derivGeom(3);
    Scalar g22 = derivGeom(1) * derivGeom(1) + derivGeom(3) * derivGeom(3);
    Scalar scaleFactor = g11 + g22;

    Scalar Lx =
        (g22 * deriv2Geom(0) + g11 * deriv2Geom(1) - 2.0 * g12 * deriv2Geom(2))
            / scaleFactor;
    Scalar Ly =
        (g22 * deriv2Geom(3) + g11 * deriv2Geom(4) - 2.0 * g12 * deriv2Geom(5))
            / scaleFactor;

    dg11dx.noalias() = 2.0 * derivGeom(0) * derivBasis.row(0);
    dg11dy.noalias() = 2.0 * derivGeom(2) * derivBasis.row(0);
    dg22dx.noalias() = 2.0 * derivGeom(1) * derivBasis.row(1);
    dg22dy.noalias() = 2.0 * derivGeom(3) * derivBasis.row(1);
    dg12dx.noalias() =
        derivGeom(1) * derivBasis.row(0) + derivGeom(0) * derivBasis.row(1);
    dg12dy.noalias() =
        derivGeom(3) * derivBasis.row(0) + derivGeom(2) * derivBasis.row(1);

    commonTerm.noalias() = g22 * deriv2Basis.row(0) + g11 * deriv2Basis.row(1)
        - 2.0 * g12 * deriv2Basis.row(2);
    dLxdx.noalias() = dg22dx * deriv2Geom(0) + dg11dx * deriv2Geom(1)
        - 2.0 * dg12dx * deriv2Geom(2) + commonTerm;
    dLydy.noalias() = dg22dy * deriv2Geom(3) + dg11dy * deriv2Geom(4)
        - 2.0 * dg12dy * deriv2Geom(5) + commonTerm;

    dLxdx = (dLxdx - Lx * (dg11dx + dg22dx)) / scaleFactor;
    dLydy = (dLydy - Ly * (dg11dy + dg22dy)) / scaleFactor;

    const index_t A = _u.cardinality() / _u.dim(); // _u.data().actives.rows()
    res.resize(_u.cardinality(), _u.cardinality());
    res.setZero();
    res.topLeftCorner(A, A).noalias() =
        _u.data().values[0].col(k) * dLxdx;
    res.bottomRightCorner(A, A).noalias() =
        _u.data().values[0].col(k) * dLydy;

    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_VALUE | NEED_GRAD | NEED_DERIV2;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV | NEED_DERIV2;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  // rowVar / colVar: the spaces of the rows and columns of the assembled
  // matrix; cardinality: the number of local (active) functions of the space
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLxDiagBlock(";
    _u.print(os);
    os << ")";
  }
};

/**
  Expression for Jacobian matrix (in H1 space) for PDE-based parameterization construction
*/
template<class E>
class jacScaledLxH1_expr : public _expr<jacScaledLxH1_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  jacScaledLxH1_expr(const E &u, const gsGeometryMap<Scalar> &G)
      : _u(u), _G(G) {}

  mutable gsMatrix<Scalar> res, derivBasis;

  const gsMatrix<Scalar> &eval(const index_t k) const {
    gsMatrix<Scalar> jacMat = _G.data().values[1].reshapeCol(k,
                                                             _G.data().dim.first,
                                                             _G.data().dim.second);
    gsMatrix<Scalar> invJacMat = jacMat.inverse();

    derivBasis = _u.data().values[1].col(k).transpose();
    derivBasis.blockTransposeInPlace(_u.dim());
    gsMatrix<Scalar> invJacHatG;
    invJacHatG.noalias() = invJacMat * derivBasis;

    const index_t N = _u.cardinality() / _u.dim(); // _u.data().actives.rows()
    gsMatrix<Scalar> jacdLxdx(N, N);
    gsMatrix<Scalar> jacdLxdy(N, N);
    gsMatrix<Scalar> jacdLydx(N, N);
    gsMatrix<Scalar> jacdLydy(N, N);

    gsMatrix<> temp(2, N);
    for (auto i = 0; i < N; ++i) {
      // for x-direction
      temp.row(0).noalias() = invJacHatG(0, i) * invJacHatG.row(0);
      temp.row(1).noalias() = invJacHatG(1, i) * invJacHatG.row(0);
      gsMatrix<Scalar> dinvJacdx(2, 2);
      dinvJacdx.row(0).noalias() = invJacHatG(0, i) * invJacMat.row(0);
      dinvJacdx.row(1).noalias() = invJacHatG(1, i) * invJacMat.row(0);
      jacdLxdx.col(i).noalias() = -temp.transpose() * invJacMat.col(0)
          - invJacHatG.transpose() * dinvJacdx.col(0);
      jacdLydx.col(i).noalias() = -temp.transpose() * invJacMat.col(1)
          - invJacHatG.transpose() * dinvJacdx.col(1);

      // for y-direction
      temp.row(0).noalias() = invJacHatG(0, i) * invJacHatG.row(1);
      temp.row(1).noalias() = invJacHatG(1, i) * invJacHatG.row(1);
      gsMatrix<Scalar> dinvJacdy(2, 2);
      dinvJacdy.row(0).noalias() = invJacHatG(0, i) * invJacMat.row(1);
      dinvJacdy.row(1).noalias() = invJacHatG(1, i) * invJacMat.row(1);
      jacdLxdy.col(i).noalias() = -temp.transpose() * invJacMat.col(0)
          - invJacHatG.transpose() * dinvJacdy.col(0);
      jacdLydy.col(i).noalias() = -temp.transpose() * invJacMat.col(1)
          - invJacHatG.transpose() * dinvJacdy.col(1);
    }

    res.resize(_u.cardinality(), _u.cardinality());
    res.setZero();
    res.topLeftCorner(N, N).noalias() = jacdLxdx;
    res.topRightCorner(N, N).noalias() = jacdLxdy;
    res.bottomLeftCorner(N, N).noalias() = jacdLydx;
    res.bottomRightCorner(N, N).noalias() = jacdLydy;

    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_GRAD;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  // rowVar / colVar: the spaces of the rows and columns of the assembled
  // matrix; cardinality: the number of local (active) functions of the space
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLx(";
    _u.print(os);
    os << ")";
  }
};

/**
  Expression for Jacobian matrix (in H1 space) for PDE-based parameterization construction
*/
template<class E>
class jacScaledLxH1DiagBlock_expr
    : public _expr<jacScaledLxH1DiagBlock_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  jacScaledLxH1DiagBlock_expr(const E &u, const gsGeometryMap<Scalar> &G)
      : _u(u), _G(G) {}

  mutable gsMatrix<Scalar> res, derivBasis;

  const gsMatrix<Scalar> &eval(const index_t k) const {
    gsMatrix<Scalar> jacMat = _G.data().values[1].reshapeCol(k,
                                                             _G.data().dim.first,
                                                             _G.data().dim.second);
    gsMatrix<Scalar> invJacMat = jacMat.inverse();

    derivBasis = _u.data().values[1].col(k).transpose();
    derivBasis.blockTransposeInPlace(_u.dim());
    gsMatrix<Scalar> invJacHatG = invJacMat * derivBasis;

    const index_t N = _u.cardinality() / _u.dim(); // _u.data().actives.rows()
    // diagonal blocks only: d(L_x)/dx and d(L_y)/dy
    gsMatrix<Scalar> jacdLxdx(N, N);
    gsMatrix<Scalar> jacdLydy(N, N);

    gsMatrix<> temp(2, N);
    for (auto i = 0; i < N; ++i) {
      // for x-direction
      temp.row(0) = invJacHatG(0, i) * invJacHatG.row(0);
      temp.row(1) = invJacHatG(1, i) * invJacHatG.row(0);
      gsMatrix<Scalar> dinvJacdx(2, 2);
      dinvJacdx.row(0) = invJacHatG(0, i) * invJacMat.row(0);
      dinvJacdx.row(1) = invJacHatG(1, i) * invJacMat.row(0);
      jacdLxdx.col(i) = -temp.transpose() * invJacMat.col(0)
          - invJacHatG.transpose() * dinvJacdx.col(0);

      // for y-direction
      temp.row(0) = invJacHatG(0, i) * invJacHatG.row(1);
      temp.row(1) = invJacHatG(1, i) * invJacHatG.row(1);
      gsMatrix<Scalar> dinvJacdy(2, 2);
      dinvJacdy.row(0) = invJacHatG(0, i) * invJacMat.row(1);
      dinvJacdy.row(1) = invJacHatG(1, i) * invJacMat.row(1);
      jacdLydy.col(i) = -temp.transpose() * invJacMat.col(1)
          - invJacHatG.transpose() * dinvJacdy.col(1);
    }

    res.resize(_u.cardinality(), _u.cardinality());
    res.setZero();
    res.topLeftCorner(N, N) = jacdLxdx;
    res.bottomRightCorner(N, N) = jacdLydy;

    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_GRAD;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  // rowVar / colVar: the spaces of the rows and columns of the assembled
  // matrix; cardinality: the number of local (active) functions of the space
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLx(";
    _u.print(os);
    os << ")";
  }
};

/**

  Expression for the product of the jacob of space and matrix A. Return term
  is a column vector (d*n x 1). The result should be the same as the result
  from jac(space) % A, whereas avoids redundant multiply with zero-value components.

 For 2D case:
 frprod2(space,A) = [a_{11}*dN_1/dxi_1, a_{12}*dN_1/dxi_2, a_{11}*dN_2/dxi_1,
 a_{12}*dN_2/dxi_2,...,a_{11}*dN_n/dxi_1, a_{12}*dN_n/dxi_2, a_{21}*dN_1/dxi_1,
 a_{22}*dN_1/dxi_2, a_{21}*dN_2/dxi_1, a_{22}*dN_2/dxi_2,...,a_{21}*dN_n/dxi_1,
 a_{22}*dN_n/dxi_2]
 For 3D case:
 frprod2(space,A) = [a_{11}*dN_1/dxi_1, a_{12}*dN_1/dxi_2, a_{13}*dN_1/dxi_3,...
 a_{11}*dN_2/dxi_1, a_{12}*dN_2/dxi_1, a_{13}*dN_2/dxi_3,...,
 a_{11}*dN_n/dxi_1, a_{12}*dN_n/dxi_1, a_{13}*dN_n/dxi_3,
 a_{21}*dN_1/dxi_1, a_{22}*dN_1/dxi_2, a_{23}*dN_1/dxi_3,...
 a_{21}*dN_2/dxi_1, a_{22}*dN_2/dxi_1, a_{23}*dN_2/dxi_3,...,
 a_{21}*dN_n/dxi_1, a_{22}*dN_n/dxi_1, a_{23}*dN_n/dxi_3,
 a_{31}*dN_1/dxi_1, a_{32}*dN_1/dxi_2, a_{33}*dN_1/dxi_3,...
 a_{31}*dN_2/dxi_1, a_{32}*dN_2/dxi_1, a_{33}*dN_2/dxi_3,...,
 a_{31}*dN_n/dxi_1, a_{32}*dN_n/dxi_1, a_{33}*dN_n/dxi_3,]

NOTE: _u should be a space, _v should NOT be a space (fix with assert)

 */
template<typename E1, typename E2>
class frprod2_expr : public _expr<frprod2_expr<E1, E2> > {
 public:
  typedef typename E2::Scalar Scalar;
  enum { ScalarValued = 0, Space = E1::Space, ColBlocks = 0 };

 private:
  typename E1::Nested_t _u;
  typename E2::Nested_t _v;

  mutable gsMatrix<Scalar> res, bGrads;

 public:

  frprod2_expr(_expr<E1> const &u, _expr<E2> const &v)
      : _u(u), _v(v) {}

  const gsMatrix<Scalar> &eval(const index_t k) const //todo: specialize for nb==1
  {
    auto A = _v.eval(k);
    bGrads = _u.data().values[1].col(k).transpose();
    bGrads.blockTransposeInPlace(_u.dim());
    res.noalias() = A * bGrads;
    res.transposeInPlace();
    res.resize(1, _u.cardinality());
    res.transposeInPlace();
    return res;
  }

  index_t rows() const { return 1; }
  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    _v.parse(evList);
    evList.add(_u);
    _u.data().flags |= NEED_GRAD;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _v.rowVar(); }

  void print(std::ostream &os) const {
    os << "(";
    _u.print(os);
    os << " % ";
    _v.print(os);
    os << ")";
  }
};

/// jacobian matrix of scaled Lx for PDE-based parameterization construction
template<class E>
EIGEN_STRONG_INLINE
jacScaledLx_expr<E> jacScaledLx(const E &u,
                                const gsGeometryMap<typename E::Scalar> &G) {
  return jacScaledLx_expr<E>(u, G);
}

/// diagonal part of jacobian matrix of scaled Lx for PDE-based parameterization construction
template<class E>
EIGEN_STRONG_INLINE
jacScaledLxDiag_expr<E> jacScaledLxDiag(const E &u,
                                        const gsGeometryMap<typename E::Scalar> &G) {
  return jacScaledLxDiag_expr<E>(u, G);
}

/// diagonal block part of jacobian matrix of scaled Lx for PDE-based parameterization construction
template<class E>
EIGEN_STRONG_INLINE
jacScaledLxDiagBlock_expr<E> jacScaledLxDiagBlock(const E &u,
                                                  const gsGeometryMap<typename E::Scalar> &G) {
  return jacScaledLxDiagBlock_expr<E>(u, G);
}

/// jacobian matrix of scaled Lx (in H1 space) for PDE-based parameterization construction
template<class E>
EIGEN_STRONG_INLINE
jacScaledLxH1_expr<E> jacScaledLxH1(const E &u,
                                    const gsGeometryMap<typename E::Scalar> &G) {
  return jacScaledLxH1_expr<E>(u, G);
}

/// jacobian matrix of scaled Lx (in H1 space) for PDE-based parameterization construction
template<class E>
EIGEN_STRONG_INLINE
jacScaledLxH1DiagBlock_expr<E> jacScaledLxH1DiagBlock(const E &u,
                                                      const gsGeometryMap<
                                                          typename E::Scalar> &G) {
  return jacScaledLxH1DiagBlock_expr<E>(u, G);
}

// Frobenious product (also known as double dot product) operator for expressions
template<typename E1, typename E2>
EIGEN_STRONG_INLINE
frprod2_expr<E1, E2> const frprod2(E1 const &u,
                                   E2 const &M) {
  return frprod2_expr<E1, E2>(u, M);
}

/**
  Residual of the (scaled) elliptic grid generation equations in
  non-divergence form (H2 discretization), in two and three dimensions:

     R_{c,i} = phi_i * ( adj(g) : Hess(x_c) ) / trace(adj(g)),   g = J^T J

  For planar domains trace(adj(g)) = trace(g), i.e. this is the residual
  frprod3(u, hess(G) % (jac(G).tr()*jac(G)).adj() / trace).
*/
template<class E>
class scaledEGGRes_expr : public _expr<scaledEGGRes_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;

  // matrices of size dim x dim, dim <= 3, without heap allocation
  typedef gsEigen::Matrix<Scalar, Dynamic, Dynamic, 0, 3, 3> SmallMatrix;
  typedef gsEigen::Matrix<Scalar, Dynamic, 1, 0, 6, 1> SmallVector;

  mutable gsMatrix<Scalar> res;
  mutable SmallMatrix jacGeom, cof, adjMetric;
  mutable SmallVector weights;

 public:
  enum { Space = 1, ScalarValued = 0, ColBlocks = 0 };

  scaledEGGRes_expr(const E &u, const gsGeometryMap<Scalar> &G)
      : _u(u), _G(G) {}

  const gsMatrix<Scalar> &eval(const index_t k) const {
    const index_t dim = _u.dim();
    const index_t A = _u.cardinality() / dim;
    const index_t numDeriv2 = dim * (dim + 1) / 2;

    // jacGeom(c,a) = d x_c / d xi_a
    jacGeom = _G.data().values[1].reshapeCol(k, dim, dim).transpose();
    cofactorMatrix(jacGeom, cof);
    adjMetric.noalias() = cof.transpose() * cof;
    packSymmetric(adjMetric, weights);
    const Scalar scaleFactor = adjMetric.trace();

    const auto deriv2Geom = _G.data().values[2].col(k);
    const auto basis = _u.data().values[0].col(k);

    res.resize(_u.cardinality(), 1);
    for (index_t c = 0; c != dim; ++c) {
      const Scalar Lx =
          weights.dot(deriv2Geom.segment(c * numDeriv2, numDeriv2))
              / scaleFactor;
      res.middleRows(c * A, A).noalias() = Lx * basis;
    }
    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_VALUE;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV | NEED_DERIV2;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return gsNullExpr<Scalar>::get(); }
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "scaledEGGRes(";
    _u.print(os);
    os << ")";
  }
};

/**
  Exact Jacobian matrix of scaledEGGRes_expr with respect to the control
  points, in two and three dimensions.
  mode 0: full matrix, 1: diagonal only, 2: diagonal blocks only
*/
template<class E>
class scaledEGGJac_expr : public _expr<scaledEGGJac_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;
  index_t _mode;

  typedef gsEigen::Matrix<Scalar, Dynamic, Dynamic, 0, 3, 3> SmallMatrix;
  typedef gsEigen::Matrix<Scalar, Dynamic, 1, 0, 6, 1> SmallVector;

  mutable gsMatrix<Scalar> res, derivBasis, commonTerm, dLx;
  mutable SmallMatrix jacGeom, cof, dCof, adjMetric, dAdjMetric, value;
  mutable SmallVector weights, dWeights, Lx;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  scaledEGGJac_expr(const E &u, const gsGeometryMap<Scalar> &G,
                    const index_t mode)
      : _u(u), _G(G), _mode(mode) {}

  const gsMatrix<Scalar> &eval(const index_t k) const {
    const index_t dim = _u.dim();
    const index_t A = _u.cardinality() / dim;
    const index_t numDeriv2 = dim * (dim + 1) / 2;

    jacGeom = _G.data().values[1].reshapeCol(k, dim, dim).transpose();
    cofactorMatrix(jacGeom, cof);
    adjMetric.noalias() = cof.transpose() * cof;
    packSymmetric(adjMetric, weights);
    const Scalar scaleFactor = adjMetric.trace();

    const auto deriv2Geom = _G.data().values[2].col(k);
    const auto basis = _u.data().values[0].col(k);
    derivBasis = _u.data().values[1].reshapeCol(k, dim, A);
    const auto deriv2Basis = _u.data().values[2].reshapeCol(k, numDeriv2, A);

    Lx.resize(dim);
    for (index_t c = 0; c != dim; ++c) {
      Lx(c) = weights.dot(deriv2Geom.segment(c * numDeriv2, numDeriv2))
          / scaleFactor;
    }
    // adj(g) : Hess(phi_j) / scale
    commonTerm.noalias() = (weights.transpose() * deriv2Basis) / scaleFactor;

    res.resize(_u.cardinality(), _u.cardinality());
    if (0 != _mode) { res.setZero(); }

    // The variation of adj(g) = cof^T cof, hence of the residual, is linear
    // in grad(phi_j): it is computed for the unit vectors and combined.
    gsEigen::Matrix<Scalar, 3, 1> unit;
    value.resize(dim, dim);
    for (index_t b = 0; b != dim; ++b) {
      for (index_t a = 0; a != dim; ++a) {
        unit.setZero();
        unit(a) = 1;
        cofactorDerivative(jacGeom, b, unit, dCof);
        dAdjMetric.noalias() = dCof.transpose() * cof;
        dAdjMetric += SmallMatrix(dAdjMetric.transpose());
        packSymmetric(dAdjMetric, dWeights);
        const Scalar dScale = dAdjMetric.trace();

        for (index_t c = 0; c != dim; ++c) {
          value(c, a) = (dWeights.dot(
              deriv2Geom.segment(c * numDeriv2, numDeriv2))
              - Lx(c) * dScale) / scaleFactor;
        }
      }

      // derivative of Lx(c) with respect to the coordinate b of all
      // control points (row vector), times the test functions
      for (index_t c = 0; c != dim; ++c) {
        if (0 != _mode && c != b) { continue; }
        dLx.noalias() = value.row(c) * derivBasis;
        if (c == b) { dLx += commonTerm; }
        if (1 == _mode) {
          res.block(b * A, b * A, A, A).diagonal() =
              basis.cwiseProduct(dLx.transpose());
        } else {
          res.block(c * A, b * A, A, A).noalias() = basis * dLx;
        }
      }
    }
    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_VALUE | NEED_GRAD | NEED_DERIV2;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV | NEED_DERIV2;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "scaledEGGJac(";
    _u.print(os);
    os << ")";
  }
};

/**
  Residual of the H1 discretization of the PDE-based parameterization
  (cf. jacU % invJac with invJac = jac(G).inv(), jacU = jac(u) * invJac),
  where 1/det(J) is replaced by 1/chi(det(J),eps):

     R_{c,i} = grad(phi_i) . ( adj(g) e_c ) / chi^2,
     g = J^T J,  adj(g) = cof(J)^T cof(J)

  It is well-defined for every J. Two and three dimensions.
*/
template<class E>
class regularizedH1Res_expr : public _expr<regularizedH1Res_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;
  Scalar _eps, _orientation;

  typedef gsEigen::Matrix<Scalar, Dynamic, Dynamic, 0, 3, 3> SmallMatrix;

  mutable gsMatrix<Scalar> res, flux, derivBasis;
  mutable SmallMatrix jacGeom, cof, K;

 public:
  enum { Space = 1, ScalarValued = 0, ColBlocks = 0 };

  regularizedH1Res_expr(const E &u, const gsGeometryMap<Scalar> &G,
                        Scalar eps, Scalar orientation)
      : _u(u), _G(G), _eps(eps), _orientation(orientation) {}

  const gsMatrix<Scalar> &eval(const index_t k) const {
    const index_t dim = _u.dim();
    const index_t A = _u.cardinality() / dim;

    // jacGeom(c,a) = d x_c / d xi_a
    jacGeom = _G.data().values[1].reshapeCol(k, dim, dim).transpose();
    cofactorMatrix(jacGeom, cof);
    const Scalar det = _orientation * jacGeom.row(0).dot(cof.row(0));
    const Scalar chi = garanzhaChi(det, _eps);

    derivBasis = _u.data().values[1].reshapeCol(k, dim, A);
    res.resize(_u.cardinality(), 1);
    if (!detIsControllable(_u.mapper(), _u.data().actives, _u.data().patchId,
                           cof, derivBasis)) {
      // det(J) is determined by the boundary: no contribution
      res.setZero();
      return res;
    }

    K.noalias() = (cof.transpose() * cof) / (chi * chi);
    flux.noalias() = K * derivBasis;
    for (index_t c = 0; c != dim; ++c) {
      res.middleRows(c * A, A) = flux.row(c).transpose();
    }
    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_GRAD | NEED_ACTIVE;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return gsNullExpr<Scalar>::get(); }
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "regularizedH1Res(";
    _u.print(os);
    os << ")";
  }
};

/**
  Exact Jacobian matrix of regularizedH1Res_expr with respect to the control
  points (optionally only its diagonal blocks). Two and three dimensions.
*/
template<class E>
class regularizedH1Jac_expr : public _expr<regularizedH1Jac_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;
  Scalar _eps, _orientation;
  bool _diagBlocksOnly;

  typedef gsEigen::Matrix<Scalar, Dynamic, Dynamic, 0, 3, 3> SmallMatrix;

  mutable gsMatrix<Scalar> res, derivBasis, flux[3];
  mutable SmallMatrix jacGeom, cof, dCof, K, dK;

 public:
  enum { Space = 3, ScalarValued = 0, ColBlocks = 0 };

  regularizedH1Jac_expr(const E &u, const gsGeometryMap<Scalar> &G,
                        Scalar eps, Scalar orientation, bool diagBlocksOnly)
      : _u(u), _G(G), _eps(eps), _orientation(orientation),
        _diagBlocksOnly(diagBlocksOnly) {}

  const gsMatrix<Scalar> &eval(const index_t k) const {
    const index_t dim = _u.dim();
    const index_t A = _u.cardinality() / dim;

    jacGeom = _G.data().values[1].reshapeCol(k, dim, dim).transpose();
    cofactorMatrix(jacGeom, cof);
    const Scalar det = _orientation * jacGeom.row(0).dot(cof.row(0));
    const Scalar root = math::sqrt(det * det + _eps * _eps);
    const Scalar chi = garanzhaChi(det, _eps);
    const Scalar w = 1.0 / (chi * chi);

    K.noalias() = (cof.transpose() * cof) * w;
    derivBasis = _u.data().values[1].reshapeCol(k, dim, A);

    res.resize(_u.cardinality(), _u.cardinality());
    if (!detIsControllable(_u.mapper(), _u.data().actives, _u.data().patchId,
                           cof, derivBasis)) {
      res.setZero();
      return res;
    }
    if (_diagBlocksOnly) { res.setZero(); }

    // The variation of K = adj(g)/chi^2 is linear in grad(phi_j): it is
    // computed for the unit vectors and combined.
    gsEigen::Matrix<Scalar, 3, 1> unit;
    for (index_t b = 0; b != dim; ++b) {
      for (index_t a = 0; a != dim; ++a) {
        unit.setZero();
        unit(a) = 1;
        cofactorDerivative(jacGeom, b, unit, dCof);
        // d(1/chi^2) = dInvChi2 / chi^2, using chi'/chi = 1/sqrt(D^2+eps^2)
        // and d(det J) = cof.row(b) . grad(phi_j)
        const Scalar dInvChi2 = -2.0 * _orientation * cof(b, a) / root;

        dK.noalias() = dCof.transpose() * cof;
        dK += SmallMatrix(dK.transpose());
        dK = w * dK + dInvChi2 * K;
        flux[a].noalias() = dK * derivBasis;
      }

      for (index_t c = 0; c != dim; ++c) {
        if (_diagBlocksOnly && c != b) { continue; }
        auto block = res.block(c * A, b * A, A, A);
        block.noalias() = flux[0].row(c).transpose() * derivBasis.row(0);
        for (index_t a = 1; a != dim; ++a) {
          block.noalias() += flux[a].row(c).transpose() * derivBasis.row(a);
        }
      }
    }
    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_GRAD | NEED_ACTIVE;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _u.rowVar(); }
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "regularizedH1Jac(";
    _u.print(os);
    os << ")";
  }
};

/**
  Gradient (isHessian = false) and Gauss-Newton Hessian (isHessian = true) of
  Garanzha's regularized Winslow functional in dimension d = 2, 3,

     W(x) = int |J|_F^2 / chi(det(J),eps)^(2/d),

  with respect to the control points. Writing W = int |r|^2 with
  r = J / chi^(1/d), the Gauss-Newton Hessian 2 (dr)^T (dr) is symmetric
  positive definite for eps > 0.
*/
template<class E, bool isHessian>
class garanzhaWinslow_expr
    : public _expr<garanzhaWinslow_expr<E, isHessian> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;
  Scalar _eps, _orientation;

  typedef gsEigen::Matrix<Scalar, Dynamic, Dynamic, 0, 3, 3> SmallMatrix;

  mutable gsMatrix<Scalar> res, derivBasis, p, q, laplace;
  mutable gsVector<Scalar> pVec, qVec;
  mutable SmallMatrix jacGeom, cof;

 public:
  enum { Space = (isHessian ? 3 : 1), ScalarValued = 0, ColBlocks = 0 };

  garanzhaWinslow_expr(const E &u, const gsGeometryMap<Scalar> &G,
                       Scalar eps, Scalar orientation)
      : _u(u), _G(G), _eps(eps), _orientation(orientation) {}

  const gsMatrix<Scalar> &eval(const index_t k) const {
    const index_t dim = _u.dim();
    const index_t A = _u.cardinality() / dim;

    jacGeom = _G.data().values[1].reshapeCol(k, dim, dim).transpose();
    cofactorMatrix(jacGeom, cof);
    const Scalar det = _orientation * jacGeom.row(0).dot(cof.row(0));
    const Scalar sqNorm = jacGeom.squaredNorm();
    const Scalar root = math::sqrt(det * det + _eps * _eps);
    const Scalar chi = garanzhaChi(det, _eps);
    const Scalar rho = 1.0 / root; // chi'/chi
    const Scalar chiPow = (2 == dim ? chi : math::pow(chi, Scalar(2) / dim));

    derivBasis = _u.data().values[1].reshapeCol(k, dim, A);
    if (!detIsControllable(_u.mapper(), _u.data().actives, _u.data().patchId,
                           cof, derivBasis)) {
      // det(J) is determined by the boundary: the point does not take part
      // in the functional (see garanzhaWinslowEnergy_expr)
      res.setZero(_u.cardinality(), isHessian ? _u.cardinality() : 1);
      return res;
    }
    // p.row(c) = grad(x_c) . grad(phi),  q.row(c) = d(det J)/d(x_c)
    p.noalias() = jacGeom * derivBasis;
    q.noalias() = _orientation * (cof * derivBasis);

    if (!isHessian) {
      res.resize(_u.cardinality(), 1);
      const Scalar c1 = (Scalar(2) / dim) * rho * sqNorm;
      for (index_t c = 0; c != dim; ++c) {
        res.middleRows(c * A, A).noalias() =
            ((2.0 * p.row(c) - c1 * q.row(c)) / chiPow).transpose();
      }
    } else {
      res.resize(_u.cardinality(), _u.cardinality());
      const Scalar w = 2.0 / chiPow;
      const Scalar c1 = -(Scalar(1) / dim) * rho * w;
      const Scalar c2 = (Scalar(1) / (dim * dim)) * rho * rho * sqNorm * w;
      // block (c,b) = c1 (p_c^T q_b + q_c^T p_b) + c2 q_c^T q_b
      //               + delta_cb w grad(phi)^T grad(phi),
      // i.e. a rank-two matrix plus diagonal blocks
      pVec.resize(_u.cardinality());
      qVec.resize(_u.cardinality());
      for (index_t c = 0; c != dim; ++c) {
        pVec.segment(c * A, A) = p.row(c).transpose();
        qVec.segment(c * A, A) = q.row(c).transpose();
      }
      res.noalias() = (c1 * pVec) * qVec.transpose();
      res.noalias() += qVec * (c1 * pVec + c2 * qVec).transpose();
      laplace.noalias() = w * (derivBasis.transpose() * derivBasis);
      for (index_t c = 0; c != dim; ++c) {
        res.block(c * A, c * A, A, A) += laplace;
      }
    }
    return res;
  }

  index_t rows() const { return 1; }

  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_GRAD | NEED_ACTIVE;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV;
  }

  const gsFeSpace<Scalar> &rowVar() const { return _u.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const {
    return isHessian ? _u.rowVar() : gsNullExpr<Scalar>::get();
  }
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << (isHessian ? "garanzhaWinslowHess(" : "garanzhaWinslowGrad(");
    _u.print(os);
    os << ")";
  }
};

/**
  Density of Garanzha's regularized Winslow functional in dimension d = 2, 3:

     |J|_F^2 / chi(det(J),eps)^(2/d)

  Points at which det(J) cannot be changed by the free control points of the
  space \a u (see detIsControllable) do not take part in the functional.
*/
template<class E>
class garanzhaWinslowEnergy_expr
    : public _expr<garanzhaWinslowEnergy_expr<E> > {
 public:
  typedef typename E::Scalar Scalar;

 private:
  typename E::Nested_t _u;
  typename gsGeometryMap<Scalar>::Nested_t _G;
  Scalar _eps, _orientation;

  typedef gsEigen::Matrix<Scalar, Dynamic, Dynamic, 0, 3, 3> SmallMatrix;

  mutable SmallMatrix jacGeom, cof;
  mutable Scalar res;

 public:
  enum { Space = 0, ScalarValued = 1, ColBlocks = 0 };

  garanzhaWinslowEnergy_expr(const E &u, const gsGeometryMap<Scalar> &G,
                             Scalar eps, Scalar orientation)
      : _u(u), _G(G), _eps(eps), _orientation(orientation) {}

  Scalar eval(const index_t k) const {
    const index_t dim = _u.dim();
    const index_t A = _u.data().values[1].rows() / dim;
    jacGeom = _G.data().values[1].reshapeCol(k, dim, dim).transpose();
    cofactorMatrix(jacGeom, cof);
    if (!detIsControllable(_u.mapper(), _u.data().actives, _u.data().patchId,
                           cof, _u.data().values[1].reshapeCol(k, dim, A))) {
      res = 0;
      return res;
    }
    const Scalar det = _orientation * jacGeom.row(0).dot(cof.row(0));
    const Scalar chi = garanzhaChi(det, _eps);
    const Scalar chiPow = (2 == dim ? chi : math::pow(chi, Scalar(2) / dim));
    res = jacGeom.squaredNorm() / chiPow;
    return res;
  }

  index_t rows() const { return 0; }

  index_t cols() const { return 0; }

  void parse(gsExprHelper<Scalar> &evList) const {
    evList.add(_u);
    _u.data().flags |= NEED_GRAD | NEED_ACTIVE;

    evList.add(_G);
    _G.data().flags |= NEED_DERIV;
  }

  const gsFeSpace<Scalar> &rowVar() const { return gsNullExpr<Scalar>::get(); }
  const gsFeSpace<Scalar> &colVar() const { return gsNullExpr<Scalar>::get(); }

  void print(std::ostream &os) const {
    os << "garanzhaWinslowEnergy(";
    _G.print(os);
    os << ")";
  }
};

/// density of Garanzha's regularized Winslow functional
template<class E>
EIGEN_STRONG_INLINE
garanzhaWinslowEnergy_expr<E> garanzhaWinslowEnergy(const E &u,
                                                    const gsGeometryMap<typename E::Scalar> &G,
                                                    const typename E::Scalar eps,
                                                    const typename E::Scalar orientation = 1) {
  return garanzhaWinslowEnergy_expr<E>(u, G, eps, orientation);
}

/// residual of the scaled elliptic grid generation equations (2D and 3D)
template<class E>
EIGEN_STRONG_INLINE
scaledEGGRes_expr<E> scaledEGGRes(const E &u,
                                  const gsGeometryMap<typename E::Scalar> &G) {
  return scaledEGGRes_expr<E>(u, G);
}

/// jacobian matrix of the scaled elliptic grid generation equations (2D and
/// 3D); mode 0: full, 1: diagonal, 2: diagonal blocks
template<class E>
EIGEN_STRONG_INLINE
scaledEGGJac_expr<E> scaledEGGJac(const E &u,
                                  const gsGeometryMap<typename E::Scalar> &G,
                                  const index_t mode = 0) {
  return scaledEGGJac_expr<E>(u, G, mode);
}

/// residual of the H1 discretization with Garanzha's regularization
template<class E>
EIGEN_STRONG_INLINE
regularizedH1Res_expr<E> regularizedH1Res(const E &u,
                                          const gsGeometryMap<typename E::Scalar> &G,
                                          const typename E::Scalar eps,
                                          const typename E::Scalar orientation = 1) {
  return regularizedH1Res_expr<E>(u, G, eps, orientation);
}

/// jacobian matrix of the H1 discretization with Garanzha's regularization
template<class E>
EIGEN_STRONG_INLINE
regularizedH1Jac_expr<E> regularizedH1Jac(const E &u,
                                          const gsGeometryMap<typename E::Scalar> &G,
                                          const typename E::Scalar eps,
                                          const typename E::Scalar orientation = 1,
                                          const bool diagBlocksOnly = false) {
  return regularizedH1Jac_expr<E>(u, G, eps, orientation, diagBlocksOnly);
}

/// gradient of Garanzha's regularized Winslow functional
template<class E>
EIGEN_STRONG_INLINE
garanzhaWinslow_expr<E, false> garanzhaWinslowGrad(const E &u,
                                                   const gsGeometryMap<typename E::Scalar> &G,
                                                   const typename E::Scalar eps,
                                                   const typename E::Scalar orientation = 1) {
  return garanzhaWinslow_expr<E, false>(u, G, eps, orientation);
}

/// Gauss-Newton Hessian of Garanzha's regularized Winslow functional
template<class E>
EIGEN_STRONG_INLINE
garanzhaWinslow_expr<E, true> garanzhaWinslowHess(const E &u,
                                                  const gsGeometryMap<typename E::Scalar> &G,
                                                  const typename E::Scalar eps,
                                                  const typename E::Scalar orientation = 1) {
  return garanzhaWinslow_expr<E, true>(u, G, eps, orientation);
}

/// Ternary expression for conditional evaluation
template <class E0, class E1, class E2>
class ternary_expr : public _expr<ternary_expr<E0, E1, E2>> {
  typename E0::Nested_t _u;
  typename E1::Nested_t _v;
  typename E2::Nested_t _w;

public:
  typedef typename E1::Scalar Scalar;

  explicit ternary_expr(_expr<E0> const &u, _expr<E1> const &v,
                        _expr<E2> const &w)
      : _u(u), _v(v), _w(w) {
    GISMO_ASSERT(E0::ScalarValued, "Condition must be scalar valued");
    GISMO_ASSERT((int)E1::ScalarValued == (int)E2::ScalarValued,
                 "Both v and w must be scalar valued (or not).");
    GISMO_ASSERT((int)E1::ColBlocks == (int)E2::ColBlocks,
                 "Both v and w must be colblocks (or not).");
    GISMO_ASSERT((int)E1::Space == (int)E2::Space,
                 "Both v and w must be space (or not), but E1::Space = "
                     << E1::Space << " and E2::Space = " << E2::Space);
    GISMO_ASSERT(_v.rows() == _w.rows(), "Rows of v and w differ. _v.rows() = "
                                             << _v.rows()
                                             << ", _w.rows() = " << _w.rows());
    GISMO_ASSERT(_v.cols() == _w.cols(),
                 "Columns of v and w differ. _v.cols() = "
                     << _v.cols() << ", _w.cols() = " << _w.cols());
    GISMO_ASSERT(_v.rowVar() == _w.rowVar(), "rowVar of v and w differ.");
    GISMO_ASSERT(_v.colVar() == _w.colVar(), "colVar of v and w differ.");
  }

public:
  enum {
    ScalarValued = E1::ScalarValued,
    ColBlocks = E1::ColBlocks,
    Space = E1::Space
  }; // == E2::Space

  const Temporary_t eval(const index_t k) const {
    return (_u.eval(k) > 0 ? _v.eval(k) : _w.eval(k));
  }

  index_t rows() const { return _v.rows(); }
  index_t cols() const { return _v.cols(); }
  void parse(gsExprHelper<Scalar> &evList) const {
    _u.parse(evList);
    _v.parse(evList);
    _w.parse(evList);
  }

  const gsFeSpace<Scalar> &rowVar() const { return _v.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _v.colVar(); }
};

/// Ternary ternary_expr
template<class E0, class E1, class E2>
EIGEN_STRONG_INLINE ternary_expr<E0, E1, E2> ternary(const E0 &u,
                                                     const E1 &v,
                                                     const E2 &w) {
  return ternary_expr<E0, E1, E2>(u, v, w);
}

} // namespace expr

namespace preAApp {

template<typename Scalar>
void NodalAMG<Scalar>::setNodes(const std::vector<index_t> &nodes, const std::vector<index_t> &components) {
  GISMO_ENSURE(nodes.size() == components.size(),
               "NodalAMG: nodes and components differ in size.");
  const index_t n = nodes.size();
  // order the unknowns node by node
  m_perm.resize(n);
  for (index_t i = 0; i != n; ++i) { m_perm[i] = i; }
  std::stable_sort(m_perm.begin(), m_perm.end(), [&](index_t a, index_t b) {
    return nodes[a] != nodes[b] ? nodes[a] < nodes[b]
                                : components[a] < components[b];
  });
  m_inv.resize(n);
  for (index_t i = 0; i != n; ++i) { m_inv[m_perm[i]] = i; }
  m_nodePtr0.assign(1, 0);
  m_comp0.resize(n);
  for (index_t i = 0; i != n; ++i) {
    m_comp0[i] = components[m_perm[i]];
    if (i > 0 && nodes[m_perm[i]] != nodes[m_perm[i - 1]]) {
      m_nodePtr0.push_back(i);
    }
  }
  if (n > 0) { m_nodePtr0.push_back(n); }
}

template<typename Scalar>
template<typename MatType>
NodalAMG<Scalar> &NodalAMG<Scalar>::compute(const MatType &mat) {
  const index_t n = mat.rows();
  GISMO_ENSURE(index_t(m_perm.size()) == n && mat.cols() == n,
               "NodalAMG: the nodes of the unknowns are not set.");
  m_levels.clear();
  m_levels.emplace_back();
  {
    Level &fine = m_levels.back();
    std::vector<gsEigen::Triplet<Scalar, index_t> > entries;
    entries.reserve(mat.nonZeros());
    for (index_t j = 0; j < mat.outerSize(); ++j)
      for (typename MatType::InnerIterator it(mat, j); it; ++it)
        entries.emplace_back(m_inv[it.row()], m_inv[it.col()], it.value());
    fine.A.resize(n, n);
    fine.A.setFromTriplets(entries.begin(), entries.end());
    fine.nodePtr = m_nodePtr0;
    fine.comp = m_comp0;
  }
  m_info = gsEigen::Success;
  while (true) {
    Level &L = m_levels.back();
    if (!blockDiagonalInverse(L)) { m_info = gsEigen::NumericalIssue; }
    const index_t numNodes = L.nodePtr.size() - 1;
    if (L.A.rows() <= m_settings.coarseSize
        || index_t(m_levels.size()) >= m_settings.maxLevels) { break; }
    std::vector<index_t> agg;
    const index_t numAgg = aggregate(L, agg);
    // stop if the aggregation hardly coarsens
    if (numAgg == 0 || 10 * numAgg > 9 * numNodes) { break; }
    Level C;
    L.P = prolongation(L, agg, numAgg, C);
    L.R = RowMatrix(L.P.transpose());
    const RowMatrix AP = L.A * L.P;
    C.A = RowMatrix(L.R * AP);
    m_levels.push_back(give(C));
  }
  const Level &coarsest = m_levels.back();
  m_coarseDirect = (coarsest.A.rows() <= 4 * m_settings.coarseSize);
  if (m_coarseDirect && coarsest.A.rows() > 0) {
    m_coarseLU.compute(Dense(coarsest.A));
    if (!m_coarseLU.matrixLU().diagonal().allFinite()
        || m_coarseLU.matrixLU().diagonal().cwiseAbs().minCoeff() == Scalar(0)) {
      m_info = gsEigen::NumericalIssue;
    }
  }
  return *this;
}

template<typename Scalar>
Scalar NodalAMG<Scalar>::operatorComplexity() const {
  Scalar sum = 0;
  for (const Level &L : m_levels) { sum += L.A.nonZeros(); }
  return m_levels.empty() ? Scalar(0)
                          : sum / Scalar(m_levels.front().A.nonZeros());
}

template<typename Scalar>
void NodalAMG<Scalar>::apply(const gsMatrix<Scalar> &input, gsMatrix<Scalar> &x) const {
  const index_t n = m_perm.size();
  x.resize(input.rows(), input.cols());
  if (0 == n) { return; }
  Vector b(n), y;
  for (index_t j = 0; j != input.cols(); ++j) {
    for (index_t i = 0; i != n; ++i) { b(i) = input(m_perm[i], j); }
    cycle(0, b, y);
    for (index_t i = 0; i != n; ++i) { x(m_perm[i], j) = y(i); }
  }
}

template<typename Scalar>
bool NodalAMG<Scalar>::blockDiagonalInverse(Level &L) const {
  const index_t numNodes = L.nodePtr.size() - 1;
  L.dinvPtr.assign(numNodes + 1, 0);
  for (index_t I = 0; I != numNodes; ++I) {
    const index_t bs = L.nodePtr[I + 1] - L.nodePtr[I];
    L.dinvPtr[I + 1] = L.dinvPtr[I] + bs * bs;
  }
  L.dinv.assign(L.dinvPtr.back(), Scalar(0));
  bool regular = true;
  for (index_t I = 0; I != numNodes; ++I) {
    const index_t s = L.nodePtr[I], bs = L.nodePtr[I + 1] - s;
    Dense block = Dense::Zero(bs, bs);
    for (index_t i = s; i != s + bs; ++i)
      for (typename RowMatrix::InnerIterator it(L.A, i); it; ++it)
        if (it.col() >= s && it.col() < s + bs) {
          block(i - s, it.col() - s) = it.value();
        }
    gsEigen::FullPivLU<Dense> lu(block);
    Dense inverse;
    if (lu.isInvertible()) {
      inverse = lu.inverse();
    } else {
      inverse = Dense::Zero(bs, bs);
      for (index_t k = 0; k != bs; ++k) {
        if (block(k, k) != Scalar(0)) { inverse(k, k) = 1 / block(k, k); }
        else { regular = false; }
      }
    }
    for (index_t r = 0; r != bs; ++r)
      for (index_t c = 0; c != bs; ++c)
        L.dinv[L.dinvPtr[I] + r * bs + c] = inverse(r, c);
  }
  return regular;
}

template<typename Scalar>
index_t NodalAMG<Scalar>::aggregate(const Level &L, std::vector<index_t> &agg) const {
  const index_t numNodes = L.nodePtr.size() - 1;
  std::vector<index_t> nodeOf(L.A.rows());
  for (index_t I = 0; I != numNodes; ++I)
    for (index_t i = L.nodePtr[I]; i != L.nodePtr[I + 1]; ++i) { nodeOf[i] = I; }
  // Frobenius norms of the node blocks; pass 0: diagonal blocks,
  // pass 1: strong off-diagonal couplings
  std::vector<Scalar> diag(numNodes, 0), acc(numNodes, 0);
  std::vector<char> seen(numNodes, 0);
  std::vector<index_t> touched;
  std::vector<std::vector<std::pair<index_t, Scalar> > > strong(numNodes);
  for (index_t pass = 0; pass != 2; ++pass) {
    for (index_t I = 0; I != numNodes; ++I) {
      touched.clear();
      for (index_t i = L.nodePtr[I]; i != L.nodePtr[I + 1]; ++i)
        for (typename RowMatrix::InnerIterator it(L.A, i); it; ++it) {
          const index_t J = nodeOf[it.col()];
          if (!seen[J]) { seen[J] = 1; touched.push_back(J); }
          acc[J] += it.value() * it.value();
        }
      for (index_t J : touched) {
        const Scalar norm = math::sqrt(acc[J]);
        acc[J] = 0;
        seen[J] = 0;
        if (0 == pass) {
          if (J == I) { diag[I] = norm; }
        } else if (J != I && norm > m_settings.theta * math::sqrt(diag[I] * diag[J])) {
          strong[I].emplace_back(J, norm);
          strong[J].emplace_back(I, norm);
        }
      }
    }
  }
  for (auto &list : strong) { // remove duplicates, keep the max. strength
    std::sort(list.begin(), list.end());
    size_t k = 0;
    for (size_t m = 0; m != list.size(); ++m) {
      if (k > 0 && list[k - 1].first == list[m].first) {
        list[k - 1].second = std::max(list[k - 1].second, list[m].second);
      } else {
        list[k++] = list[m];
      }
    }
    list.resize(k);
  }
  agg.assign(numNodes, -1);
  index_t numAgg = 0;
  // 1. a node and its strong neighbours, if none of them is aggregated yet
  for (index_t I = 0; I != numNodes; ++I) {
    if (agg[I] >= 0) { continue; }
    bool free = true;
    for (const auto &e : strong[I]) {
      if (agg[e.first] >= 0) { free = false; break; }
    }
    if (!free) { continue; }
    agg[I] = numAgg;
    for (const auto &e : strong[I]) { agg[e.first] = numAgg; }
    ++numAgg;
  }
  // 2. the remaining nodes join the most strongly coupled aggregate of 1.
  const std::vector<index_t> first = agg;
  for (index_t I = 0; I != numNodes; ++I) {
    if (agg[I] >= 0) { continue; }
    Scalar best = -1;
    for (const auto &e : strong[I]) {
      if (first[e.first] >= 0 && e.second > best) {
        best = e.second;
        agg[I] = first[e.first];
      }
    }
  }
  // 3. nodes without strong couplings to aggregates
  for (index_t I = 0; I != numNodes; ++I) {
    if (agg[I] >= 0) { continue; }
    agg[I] = numAgg;
    for (const auto &e : strong[I]) {
      if (agg[e.first] < 0) { agg[e.first] = numAgg; }
    }
    ++numAgg;
  }
  return numAgg;
}

template<typename Scalar>
typename NodalAMG<Scalar>::RowMatrix NodalAMG<Scalar>::prolongation(const Level &L, const std::vector<index_t> &agg, const index_t numAgg, Level &C) const {
  const index_t numNodes = L.nodePtr.size() - 1;
  // (component, number of unknowns) of every aggregate
  std::vector<std::vector<std::pair<index_t, index_t> > > comps(numAgg);
  for (index_t I = 0; I != numNodes; ++I)
    for (index_t i = L.nodePtr[I]; i != L.nodePtr[I + 1]; ++i) {
      std::vector<std::pair<index_t, index_t> > &list = comps[agg[I]];
      size_t k = 0;
      while (k != list.size() && list[k].first != L.comp[i]) { ++k; }
      if (k == list.size()) { list.emplace_back(L.comp[i], 1); }
      else { ++list[k].second; }
    }
  C.nodePtr.assign(numAgg + 1, 0);
  C.comp.clear();
  for (index_t a = 0; a != numAgg; ++a) {
    std::sort(comps[a].begin(), comps[a].end());
    C.nodePtr[a + 1] = C.nodePtr[a] + comps[a].size();
    for (const auto &p : comps[a]) { C.comp.push_back(p.first); }
  }
  std::vector<gsEigen::Triplet<Scalar, index_t> > entries;
  entries.reserve(L.A.rows());
  for (index_t I = 0; I != numNodes; ++I)
    for (index_t i = L.nodePtr[I]; i != L.nodePtr[I + 1]; ++i) {
      const index_t a = agg[I];
      for (size_t k = 0; k != comps[a].size(); ++k)
        if (comps[a][k].first == L.comp[i]) {
          entries.emplace_back(i, C.nodePtr[a] + k,
                               1 / math::sqrt(Scalar(comps[a][k].second)));
        }
    }
  RowMatrix P(L.A.rows(), C.nodePtr.back());
  P.setFromTriplets(entries.begin(), entries.end());
  if (!m_settings.smoothProlongation) { return P; }

  // P <- (I - omega D^{-1} A) P, omega = 4 / (3 rho(D^{-1} A))
  entries.clear();
  for (index_t I = 0; I != numNodes; ++I) {
    const index_t s = L.nodePtr[I], bs = L.nodePtr[I + 1] - s;
    for (index_t r = 0; r != bs; ++r)
      for (index_t c = 0; c != bs; ++c)
        entries.emplace_back(s + r, s + c, L.dinv[L.dinvPtr[I] + r * bs + c]);
  }
  RowMatrix Dinv(L.A.rows(), L.A.rows());
  Dinv.setFromTriplets(entries.begin(), entries.end());
  const RowMatrix DA = Dinv * L.A;
  Vector v = Vector::Ones(L.A.rows()), w;
  for (index_t i = 0; i < v.size(); i += 3) { v(i) = Scalar(-0.5); }
  Scalar rho = 1;
  for (index_t k = 0; k != 15; ++k) {
    w = DA * v;
    rho = w.norm() / v.norm();
    if (!(rho > 0)) { return P; }
    v = w / w.norm();
  }
  const RowMatrix smoothed = P - (Scalar(4) / (3 * rho)) * RowMatrix(DA * P);
  return smoothed;
}

template<typename Scalar>
void NodalAMG<Scalar>::smooth(const Level &L, const Vector &b, Vector &x, const bool forward) const {
  const index_t numNodes = L.nodePtr.size() - 1;
  Vector r;
  for (index_t k = 0; k != numNodes; ++k) {
    const index_t I = forward ? k : numNodes - 1 - k;
    const index_t s = L.nodePtr[I], bs = L.nodePtr[I + 1] - s;
    r.resize(bs);
    for (index_t i = s; i != s + bs; ++i) {
      Scalar sum = b(i);
      for (typename RowMatrix::InnerIterator it(L.A, i); it; ++it) {
        sum -= it.value() * x(it.col());
      }
      r(i - s) = sum;
    }
    const Scalar *D = &L.dinv[L.dinvPtr[I]];
    for (index_t p = 0; p != bs; ++p) {
      Scalar dx = 0;
      for (index_t q = 0; q != bs; ++q) { dx += D[p * bs + q] * r(q); }
      x(s + p) += dx;
    }
  }
}

template<typename Scalar>
void NodalAMG<Scalar>::cycle(const size_t l, const Vector &b, Vector &x) const {
  const Level &L = m_levels[l];
  if (l + 1 == m_levels.size()) {
    if (m_coarseDirect) {
      x = m_coarseLU.solve(b);
    } else {  // the coarsening stopped early: smoothing only
      x.setZero(b.size());
      for (index_t k = 0; k != 2 * m_settings.postSweeps; ++k) {
        smooth(L, b, x, true);
        smooth(L, b, x, false);
      }
    }
    return;
  }
  x.setZero(b.size());
  for (index_t k = 0; k != m_settings.preSweeps; ++k) { smooth(L, b, x, true); }
  const Vector coarseRhs = L.R * (b - L.A * x);
  Vector coarseSol;
  cycle(l + 1, coarseRhs, coarseSol);
  x += L.P * coarseSol;
  for (index_t k = 0; k != m_settings.postSweeps; ++k) { smooth(L, b, x, false); }
}

template<typename Scalar>
void BlockGaussSeidel<Scalar>::setBlocks(const std::vector<index_t> &blocks) {
  m_blockOf = blocks;
  const index_t numBlocks = blocks.empty() ? 0
      : *std::max_element(blocks.begin(), blocks.end()) + 1;
  m_members.assign(numBlocks, std::vector<index_t>());
  m_local.resize(blocks.size());
  for (size_t i = 0; i != blocks.size(); ++i) {
    m_local[i] = m_members[blocks[i]].size();
    m_members[blocks[i]].push_back(i);
  }
}

template<typename Scalar>
template<typename MatType>
BlockGaussSeidel<Scalar> &BlockGaussSeidel<Scalar>::compute(const MatType &mat) {
  GISMO_ENSURE(m_blockOf.size() == size_t(mat.rows()) && mat.cols() == mat.rows(),
               "BlockGaussSeidel: the blocks of the unknowns are not set.");
  const size_t numBlocks = m_members.size();
  typedef std::vector<gsEigen::Triplet<Scalar, index_t> > Entries;
  // entries[g][h]: block (g,h) of the lower block triangle
  std::vector<std::vector<Entries> > entries(numBlocks,
                                             std::vector<Entries>(numBlocks));
  for (index_t j = 0; j < mat.outerSize(); ++j) {
    for (typename MatType::InnerIterator it(mat, j); it; ++it) {
      const index_t g = m_blockOf[it.row()], h = m_blockOf[it.col()];
      if (g >= h && it.value() != Scalar(0)) {
        entries[g][h].emplace_back(m_local[it.row()], m_local[it.col()],
                                   it.value());
      }
    }
  }
  m_lower.assign(numBlocks, std::vector<SparseMatrix>(numBlocks));
  for (size_t g = 0; g != numBlocks; ++g) {
    for (size_t h = 0; h < g; ++h) {
      m_lower[g][h].resize(m_members[g].size(), m_members[h].size());
      m_lower[g][h].setFromTriplets(entries[g][h].begin(),
                                    entries[g][h].end());
    }
  }
  m_info = gsEigen::Success;
  m_lu.assign(numBlocks, nullptr);
  m_amg.assign(numBlocks, nullptr);
  for (size_t b = 0; b != numBlocks; ++b) {
    SparseMatrix block(m_members[b].size(), m_members[b].size());
    block.setFromTriplets(entries[b][b].begin(), entries[b][b].end());
    if (m_amgBlocks) {
      std::vector<index_t> nodes(m_members[b].size()), zero(nodes.size(), 0);
      for (size_t i = 0; i != nodes.size(); ++i) { nodes[i] = i; }
      m_amg[b].reset(new NodalAMG<Scalar>());
      m_amg[b]->setNodes(nodes, zero);
      m_amg[b]->compute(block);
      if (m_amg[b]->info() != gsEigen::Success) { m_info = m_amg[b]->info(); }
    } else {
      m_lu[b].reset(new gsEigen::SparseLU<SparseMatrix>());
      m_lu[b]->compute(block);
      if (m_lu[b]->info() != gsEigen::Success) { m_info = m_lu[b]->info(); }
    }
  }
  return *this;
}

template<typename Scalar>
template<typename Rhs>
gsVector<Scalar> BlockGaussSeidel<Scalar>::solve(const gsEigen::MatrixBase<Rhs> &rhs) const {
  gsVector<Scalar> result(rhs.rows());
  gsMatrix<Scalar> blockRhs;
  std::vector<gsMatrix<Scalar> > blockSol(m_members.size());
  for (size_t b = 0; b != m_members.size(); ++b) {
    const std::vector<index_t> &members = m_members[b];
    blockRhs.resize(members.size(), 1);
    for (size_t i = 0; i != members.size(); ++i) { blockRhs(i) = rhs(members[i]); }
    for (size_t h = 0; h < b; ++h) { blockRhs -= m_lower[b][h] * blockSol[h]; }
    if (m_amgBlocks) { m_amg[b]->apply(blockRhs, blockSol[b]); }
    else { blockSol[b] = m_lu[b]->solve(blockRhs); }
    for (size_t i = 0; i != members.size(); ++i) { result(members[i]) = blockSol[b](i); }
  }
  return result;
}

template<typename Scalar>
void BlockGaussSeidel<Scalar>::apply(const gsMatrix<Scalar> &input, gsMatrix<Scalar> &x) const {
  x.resize(input.rows(), input.cols());
  for (index_t j = 0; j != input.cols(); ++j) { x.col(j) = solve(input.col(j)); }
}

} // namespace preAApp

}// namespace gismo
