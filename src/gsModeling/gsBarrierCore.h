/** @file gsBarrierCore.h

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


#include <gsUtils/gsStopwatch.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasis.h>

#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsSolver/gsGMRes.h>


#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif

//#include <gsLBFGSpp/gsLBFGSpp.h>
//#include <gsPreAA/gsPreAA.h>

namespace gismo {
/**
 * \brief Computes a patch parametrization given a set of boundary geometries.
 * Parametrization is not guaranteed to be non-singular. Works for planar surfaces and volumes.
 *
 * \tparam d domain dimension
 * \tparam T Coefficient type
 * \ingroup Modeling
 */

/// Convert the free control points of a multi-patch into a vector
template<typename T = real_t>
gsVector<T> convertMultiPatchToFreeVector(const gsMultiPatch<T> &mp,
                                          const gsDofMapper &mapper) {
  gsVector<T> freeVec(mapper.freeSize());

  for (size_t iptch = 0; iptch != mp.nPatches(); iptch++) {
    for (index_t idim = 0; idim != mp.targetDim(); idim++) {
      for (size_t idof = 0; idof != mapper.patchSize(iptch, idim); idof++)
        // if it is possible to just loop over the free index
        // since this function is called very often during optimization
        if (mapper.is_free(idof, iptch, idim)) {
          index_t idx = mapper.index(idof, iptch, idim);
          freeVec(idx) = mp.patch(iptch).coefs()(idof, idim);
        }
    }
  }
  return freeVec;
}

/// Convert free control points from a vector into a multi-patch
template<typename T = real_t>
void convertFreeVectorToMultiPatch(const gsVector<T> &gsFreeVec,
                                   const gsDofMapper &mapper,
                                   gsMultiPatch<T> &mp) {
  for (size_t iptch = 0; iptch != mp.nPatches(); iptch++) {
    for (index_t idim = 0; idim != mp.targetDim(); idim++) {
      for (size_t idof = 0; idof != mapper.patchSize(iptch, idim); idof++)
        // if it is possible to just loop over the free index
        // since this function is called very often during optimization
        if (mapper.is_free(idof, iptch, idim)) {
          index_t idx = mapper.index(idof, iptch, idim);
          mp.patch(iptch).coefs()(idof, idim) = gsFreeVec(idx);
        }
    }
  }
}

/// helper function to set optimizer options
#ifdef gsHLBFGS_ENABLED
template<typename T>
void setOptimizerOptions(gsHLBFGS<T> &optimizer, const gsOptionList &options) {
  optimizer.options().setInt("MaxIterations",
                             options.askInt("qi_MaxIterations", 1e4));
  optimizer.options().setReal("MinGradLen",
                              options.askReal("qi_MinGradientLength", 1e-4));
  optimizer.options().setReal("MinStepLen",
                              options.askReal("qi_MinStepLength", 1e-4));
  optimizer.options().setInt("Verbose", options.askInt("Verbose", 0));
}
#endif

/// helper function to verbose log
inline void verboseLog(const std::string &message, const index_t &verbose) {
  if (verbose > 0) { gsInfo << message << "\n"; }
}

/// gsBarrierCore
template<short_t d, typename T= real_t>
class gsBarrierCore
{
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  /// construct analysis-suitable parameterization
  static gsMultiPatch<T> compute(const gsMultiPatch<T> &mp,
                                 const gsDofMapper &mapper,
                                 const gsOptionList &options);

  /// Compute the area of the computational domain
  static T computeArea(const gsMultiPatch<T> &mp);

  /// Default options
  static gsOptionList defaultOptions();

  /// Coordinate component of every free unknown (blocks of the block
  /// Gauss-Seidel preconditioner)
  static std::vector<index_t> componentOfDofs(const gsMultiPatch<T> &mp,
                                              const gsDofMapper &mapper);

  /// Control point (node) of every free unknown: the free unknowns of one
  /// control point, also of one shared by several patches, get the same
  /// node number (numbers 0, 1, ...; for the algebraic multigrid method)
  static std::vector<index_t> nodeOfDofs(const gsMultiPatch<T> &mp,
                                         const gsDofMapper &mapper);

 private:
  /// Barrier function-based method
  static gsMultiPatch<T> computeBarrierPatch(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             const gsOptionList &options);

  /// Penalty function-based method (1)
  static gsMultiPatch<T> computePenaltyPatch(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             const gsOptionList &options);

  /// Penalty function-based method (2)
  static gsMultiPatch<T> computePenaltyPatch2(const gsMultiPatch<T> &mp,
                                              const gsDofMapper &mapper,
                                              const gsOptionList &options);

  /// PDE-based methods, including H2 and H1
  static gsMultiPatch<T> computePDEPatch(const gsMultiPatch<T> &mp,
                                         const gsDofMapper &mapper,
                                         const gsOptionList &options);

  /// variational harmonic method
  static gsMultiPatch<T> computeVHPatch(const gsMultiPatch<T> &mp,
                                        const gsDofMapper &mapper,
                                        const gsOptionList &options);

  /// PDE-based method, dimension-independent implementation (used for
  /// volumes): H2 discretization of the elliptic grid generation equations,
  /// followed by computeRegularizedH1()
  static gsMultiPatch<T> computePDEPatchVolume(const gsMultiPatch<T> &mp,
                                               const gsDofMapper &mapper,
                                               const gsOptionList &options);

  /**
   * @brief Robust second stage of the PDE-based method (planar domains and
   * volumes), based on Garanzha's regularization
   * chi(D,eps) = (D + sqrt(D^2 + eps^2))/2 of the Jacobian determinant D.
   *
   * 1. If the given parameterization is folded, it is untangled by minimizing
   *    the regularized Winslow functional int |J|^2 / chi(det J, eps)^(2/d).
   * 2. The H1 discretization is solved with 1/det J replaced by
   *    1/chi(det J, eps), using steps that keep the parameterization unfolded.
   *
   * @param[in]  mp       Multipatch (already scaled to the unit box)
   * @param[in]  mapper   DoF mapper of the free control points
   * @param[in]  x        Free control points to start from (H2 solution)
   * @param[in]  options  Options
   *
   * @return     The free control points of the result
   */
  static gsVector<T> computeRegularizedH1(const gsMultiPatch<T> &mp,
                                          const gsDofMapper &mapper,
                                          gsVector<T> x,
                                          const gsOptionList &options);

  /// One pass of computeRegularizedH1(): untangling (if needed) and the
  /// regularized H1 stage, with \a extraNodes additional Gauss-Lobatto
  /// nodes per direction in the elements (see gsRegularizedPDEAssembler).
  /// If \a negatives is given, it receives the number of check points with
  /// det(J) <= 0 of every element of the result.
  static gsVector<T> regularizedStage(const gsMultiPatch<T> &mp,
                                      const gsDofMapper &mapper,
                                      gsVector<T> x,
                                      const gsOptionList &options,
                                      const std::vector<index_t> &extraNodes,
                                      std::vector<index_t> *negatives);

  /// Value of preAAParam::linearSolver for GMRES in the H2 or H1 stage
  /// (option PDEPreconditioner)
  static int gmresPreconditioner(const gsOptionList &options, bool h2Stage);

  /// Use iterative linear solvers for a problem with this many unknowns?
  static bool useIterativeSolver(const index_t numDofs,
                                 const gsOptionList &options);

 protected:
  /**
   * @brief      Computes the area of a multipatch representing a full domain
   *
   * @param[in]  mp    Multipatch with interior
   *
   * @return     The area.
   */
  static T computeAreaInterior(const gsMultiPatch<T> &mp);

  /**
   * @brief      Computes the area of a multipatch representing boundary curves
   *
   * @param[in]  mp    Boundary curves
   *
   * @return     The area.
   */
  static T computeAreaBoundary(const gsMultiPatch<T> &mp);

 private:
  static void foldoverElimination(const gsMultiPatch<T> &mp,
                                  const gsDofMapper &mapper,
                                  gsVector<T> &initialGuessVector,
                                  const T &scaledArea,
                                  const gsOptionList &options);

  static void qualityImprovement(const gsMultiPatch<T> &mp,
                                 const gsDofMapper &mapper,
                                 gsVector<T> &initialGuessVector,
                                 const T &scaledArea,
                                 const gsOptionList &options);
};

#ifdef gsHLBFGS_ENABLED
template<short_t d, typename T = real_t>
class gsObjFoldoverFree : public gsOptProblem<T> {

 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjFoldoverFree(const gsMultiPatch<T> &patches, gsDofMapper mapper);

  /// Evaluates the objective function at the given point u.
  T evalObj(const gsAsConstVector<T> &u) const override;

  /// Computes the gradient of the objective function at the given point u
  // and stores it in result.
  void gradObj_into(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const override;

  void setDelta(const T &delta) { m_eps = delta; };

  /// Returns the options list for the class instance.
  gsOptionList options() { return m_options; }

  /// Sets the default options for the class instance.
  void defaultOptions();

  /// Adds the given options to the class instance's options.
  void addOptions(const gsOptionList &options);

  /// Applies the current options to the class instance.
  void applyOptions(const gsOptionList &options);

 protected:
  mutable gsMultiPatch<T> m_mp;
  const gsDofMapper m_mapper;
  const gsMultiBasis<T> m_mb;

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;

  gsOptionList m_options;

  T m_eps = T (1e-2);
};

template<short_t d, typename T>
class gsObjQualityImprovePt : public gsOptProblem<T> {
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjQualityImprovePt(const gsMultiPatch<T> &patches, gsDofMapper mapper);

  /// Evaluates the objective function at the given point u.
  T evalObj(const gsAsConstVector<T> &u) const override;

  /// Computes the gradient of the objective function at the given point u
  // and stores it in result.
  void gradObj_into(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const override;

  /// Returns the options list for the class instance.
  gsOptionList options() { return m_options; }

  /// Sets the default options for the class instance.
  void defaultOptions();

  /// Adds the given options to the class instance's options.
  void addOptions(const gsOptionList &options);

  /// Applies the current options to the class instance.
  void applyOptions(const gsOptionList &options);

 private:
  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const {
    GISMO_ERROR("The dimension of target domain should be 2 or 3.");
  }

  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const {
    GISMO_ERROR("The dimension of target domain should be 2 or 3.");
  }

 protected:
  mutable gsMultiPatch<T> m_mp;
  const gsDofMapper m_mapper;
  const gsMultiBasis<T> m_mb;

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;

  gsOptionList m_options;

  T m_lambda1 = 1.0, m_lambda2 = 1.0;
};

template<short_t d, typename T>
class gsObjVHPt : public gsOptProblem<T> {
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjVHPt(const gsMultiPatch<T> &patches,
            gsDofMapper mapper);

  T evalObj(const gsAsConstVector<T> &u) const final;

  void gradObj2_into(const gsAsConstVector<T> &u,
                     gsAsVector<T> &result) const;

  gsOptionList &options() { return m_options; }

  void setEps(T tol) { m_eps = tol; }

  void defaultOptions();

  void addOptions(const gsOptionList &options);

  void applyOptions(const gsOptionList &options);

 private:
  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  evalObj_impl(
      const gsAsConstVector<T> &u) const {GISMO_NO_IMPLEMENTATION; }

  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const {GISMO_NO_IMPLEMENTATION; }

 protected:
  mutable gsMultiPatch<T> m_mp;
  const gsDofMapper m_mapper;
  const gsMultiBasis<T> m_mb;

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;

  gsOptionList m_options;

  T m_lambda1, m_lambda2, m_eps = 4e-4; // 4e-4
};

template<short_t d, typename T>
class gsObjPenaltyPt : public gsOptProblem<T> {
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjPenaltyPt(const gsMultiPatch<T> &patches,
                 gsDofMapper mapper);

  T evalObj(const gsAsConstVector<T> &u) const final;

  void gradObj_into(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const final;

  gsOptionList &options() { return m_options; }

  void setEps(T tol) { m_eps = tol; }

  void defaultOptions();

  void addOptions(const gsOptionList &options);

  void applyOptions(const gsOptionList &options);

 private:

  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const {GISMO_NO_IMPLEMENTATION; }

 protected:
  mutable gsMultiPatch<T> m_mp;
  const gsDofMapper m_mapper;
  const gsMultiBasis<T> m_mb;

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;

  gsOptionList m_options;

  T m_lambda1, m_lambda2, m_eps = 1e-4; // 4e-4
};

template<short_t d, typename T>
class gsObjPenaltyPt2 : public gsOptProblem<T> {
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjPenaltyPt2(const gsMultiPatch<T> &patches,
                  gsDofMapper mapper);

  T evalObj(const gsAsConstVector<T> &u) const final;

  void gradObj_into(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  void setEps(T tol) { m_eps = tol; }

  gsOptionList &options() { return m_options; }

  void defaultOptions();

  void addOptions(const gsOptionList &options);

  void applyOptions(const gsOptionList &options);

 private:
  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  evalObj_impl(const gsAsConstVector<T> &u) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  evalObj_impl(
      const gsAsConstVector<T> &u) const {GISMO_NO_IMPLEMENTATION; }

  template<short_t _d>
  typename std::enable_if<_d == 2, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d == 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const;

  template<short_t _d>
  typename std::enable_if<_d != 2 && _d != 3, T>::type
  gradObj_into_impl(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const {GISMO_NO_IMPLEMENTATION; }

 protected:
  mutable gsMultiPatch<T> m_mp;
  const gsDofMapper m_mapper;
  const gsMultiBasis<T> m_mb;

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;

  gsOptionList m_options;

  T m_lambda1 = 1.0, m_lambda2 = 1.0, m_eps = 1e-3;
};
#endif

///
/// Parameters to control the preconditioned Anderson Acceleration algorithm
///
template <typename Scalar = double>
class preAAParam {
 public:
  ///
  /// Depth: The number of previous iterates used for Andreson Acceleration.
  /// Typically in application \ref m is small, say <= 10.
  /// There is usually little advantage in large \ref m.
  /// The default value is \c 5.
  /// Large values result in excessive computing time.
  ///
  int m;
  ///
  /// Absolute tolerance for convergence test.
  /// This parameter determines the absolute accuracy \f$\epsilon_{abs}\f$
  /// with which the solution is to be found. A minimization terminates when
  /// \f$||F|| < \epsilon_{abs}\f$, where \f$||\cdot||\f$ denotes the
  /// Euclidean (L2) norm. The default value is \c 1e-5.
  ///
  Scalar epsilon;
  ///
  /// Relative tolerance for convergence test.
  /// This parameter determines the relative accuracy \f$\epsilon_{rel}\f$
  /// with which the solution is to be found. A minimization terminates when
  /// \f$||F|| < \max\{\epsilon_{abs}, \epsilon_{rel}||x||\}\f$,
  /// where \f$||\cdot||\f$ denotes the Euclidean (L2) norm.
  /// The default value is \c 1e-5.
  ///
  Scalar epsilon_rel;
  ///
  /// The maximum number of iterations.
  /// The iteration process is terminated when the iteration count exceeds
  /// this parameter. Setting this parameter to zero continues an iteration
  /// process until a convergence or error. The default value is \c 0.
  ///
  int max_iterations;
  ///
  /// The step to update the preconditioner.
  /// update the preconditioner every updatePreconditionerStep step
  ///
  int updatePreconditionerStep;
  ///
  /// Mixing or damping parameter.
  /// It plays exactly the same role as it does for Picard iteration.
  /// One could vary \ref beta as the iteration.
  ///
//  Scalar beta;
  ///
  /// Use preconditioning or not?
  ///
  bool usePreconditioning;
  ///
  /// Linear solver for the preconditioner: 0: sparse LU (default),
  /// 1: GMRES with a block Gauss-Seidel preconditioner (blocks given by
  /// \ref blocks, solved by sparse LU), 2: GMRES with the algebraic
  /// multigrid preconditioner NodalAMG (nodes given by \ref nodes, components
  /// by \ref blocks), 3: GMRES with block Gauss-Seidel whose diagonal blocks
  /// are approximated by NodalAMG. If GMRES does not converge within
  /// \ref innerMaxIterations, sparse LU is used for the rest of the solve.
  ///
  int linearSolver;
  ///
  /// Relative tolerance and max. iterations of the iterative linear solver
  ///
  Scalar innerTolerance;
  int innerMaxIterations;
  ///
  /// A solution of the iterative linear solver is accepted if its true
  /// relative residual is below this bound (gismo's GMRES is left
  /// preconditioned, so its own criterion refers to the preconditioned
  /// residual). The fixed point of the iteration does not depend on the
  /// accuracy of the linear solves, only its convergence speed does.
  ///
  Scalar innerAcceptance;
  ///
  /// Block (e.g. coordinate component) of each unknown, for linearSolver
  /// 1, 2 and 3
  ///
  std::vector<index_t> blocks;
  ///
  /// Node (control point) of each unknown, for linearSolver 2
  ///
  std::vector<index_t> nodes;
  ///
  /// Safeguard against divergence (0: off): if the residual norm exceeds
  /// safeguardFactor times the smallest one so far, the iteration returns to
  /// the iterate with the smallest residual and restarts there with a fresh
  /// preconditioner and a halved (damped) plain step. The damping is relaxed
  /// again whenever the residual reaches a new minimum.
  ///
  Scalar safeguardFactor;

 public:
  ///
  /// Constructor for preAA parameters.
  /// Default values for parameters will be set when the object is created.
  ///
  preAAParam() {
    // clang-format off
    m = 5;
    epsilon = Scalar(1e-5);
    epsilon_rel = Scalar(1e-5);
    max_iterations = 100;
    updatePreconditionerStep = 10;
//    beta = Scalar(1.0);
    usePreconditioning = false;
    linearSolver = 0;
    innerTolerance = Scalar(1e-8);
    innerAcceptance = Scalar(1e-4);
    innerMaxIterations = 100;
    safeguardFactor = 0;
    // clang-format on
  }

  ///
  /// Checking the validity of preconditioned AA parameters.
  /// An `std::invalid_argument` exception will be thrown if some parameter
  /// is invalid.
  ///
  inline void check_param() const 
  {
    GISMO_ENSURE(m >= 0, "'m' must be non-negative");
    GISMO_ENSURE(epsilon >= 0, "'epsilon' must be non-negative");
    GISMO_ENSURE(epsilon_rel >= 0, "'epsilon_rel' must be non-negative");
    GISMO_ENSURE(max_iterations >= 0, "'max_iterations' must be non-negative");
//    GISMO_ENSURE(beta >= 0 && beta <= 1.0, "'beta' must be between 0 and 1");
  }
};

namespace preAApp {

template<typename Func>
inline void measureTime(Func func) {
  auto beforeTime = std::chrono::high_resolution_clock::now();

  func();  // Call the function

  auto afterTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double>(afterTime - beforeTime).count();
  printf("\nTime passed by %.8f sec. \n", duration);
}

/// @brief Algebraic multigrid preconditioner for systems whose unknowns
/// are grouped in nodes (here: the coordinates of one control point).
///
/// Aggregation-based AMG: the nodes are aggregated on the graph of strong
/// couplings between nodes (Frobenius norms of the node blocks), the
/// prolongation interpolates each component by a constant on every aggregate
/// (unsmoothed by default; optionally smoothed by one damped block Jacobi
/// step), the coarse matrices are Galerkin products R A P with R = P^T, and
/// the smoother is the nodal block Gauss-Seidel method (forward before,
/// backward after the coarse grid correction), which inverts the coupling of
/// all components of a node exactly. Unlike a partition of the unknowns by
/// coordinate component, nodal blocks are invariant under rotations of the
/// physical coordinates. One cycle is applied per call; meant as a
/// preconditioner of GMRES. The matrix need not be symmetric.
template<typename Scalar = double>
class NodalAMG : public gsLinearOperator<Scalar> {
 public:
  typedef gsEigen::SparseMatrix<Scalar, gsEigen::RowMajor, index_t> RowMatrix;
  typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, 1> Vector;
  typedef gsEigen::Matrix<Scalar, gsEigen::Dynamic, gsEigen::Dynamic> Dense;

  struct Settings {
    /// Strength threshold: node J is strongly coupled to node I if
    /// |A_IJ| > theta * sqrt(|A_II| |A_JJ|) (Frobenius norms of the blocks)
    Scalar theta = Scalar(0.08);
    /// Block Gauss-Seidel sweeps before / after the coarse grid correction
    int preSweeps = 2, postSweeps = 2;
    /// The coarsening stops at this many unknowns (dense LU below)
    index_t coarseSize = 400;
    int maxLevels = 10;
    /// Smooth the prolongation (smoothed aggregation)
    bool smoothProlongation = false;
  } settings;

  NodalAMG() {}

  /// Sets the node and the component of every unknown (the unknowns of one
  /// node should have different components)
  void setNodes(const std::vector<index_t> &nodes,
                const std::vector<index_t> &components) {
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

  /// Builds the hierarchy for the matrix \a mat
  template<typename MatType>
  NodalAMG &compute(const MatType &mat) {
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
      if (L.A.rows() <= settings.coarseSize
          || int(m_levels.size()) >= settings.maxLevels) { break; }
      std::vector<index_t> agg;
      const index_t numAgg = aggregate(L, agg);
      // stop if the aggregation hardly coarsens
      if (numAgg == 0 || 10 * numAgg > 9 * numNodes) { break; }
      Level C;
      L.P = prolongation(L, agg, numAgg, C);
      L.R = RowMatrix(L.P.transpose());
      const RowMatrix AP = L.A * L.P;
      C.A = RowMatrix(L.R * AP);
      m_levels.push_back(std::move(C));
    }
    const Level &coarsest = m_levels.back();
    m_coarseDirect = (coarsest.A.rows() <= 4 * settings.coarseSize);
    if (m_coarseDirect && coarsest.A.rows() > 0) {
      m_coarseLU.compute(Dense(coarsest.A));
      if (!m_coarseLU.matrixLU().diagonal().allFinite()
          || m_coarseLU.matrixLU().diagonal().cwiseAbs().minCoeff() == Scalar(0)) {
        m_info = gsEigen::NumericalIssue;
      }
    }
    return *this;
  }

  gsEigen::ComputationInfo info() const { return m_info; }

  /// Number of levels of the hierarchy
  index_t numLevels() const { return m_levels.size(); }

  /// Sum of the nonzeros of all levels relative to the finest level
  Scalar operatorComplexity() const {
    Scalar sum = 0;
    for (const Level &L : m_levels) { sum += L.A.nonZeros(); }
    return m_levels.empty() ? Scalar(0)
                            : sum / Scalar(m_levels.front().A.nonZeros());
  }

  /// gsLinearOperator interface: one cycle with zero initial guess
  void apply(const gsMatrix<Scalar> &input, gsMatrix<Scalar> &x) const {
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
  index_t rows() const { return m_perm.size(); }
  index_t cols() const { return m_perm.size(); }

 private:
  struct Level {
    RowMatrix A, P, R;
    /// unknowns of node I: nodePtr[I], ..., nodePtr[I+1]-1
    std::vector<index_t> nodePtr;
    /// component of every unknown
    std::vector<index_t> comp;
    /// inverses of the diagonal node blocks (row-major, packed)
    std::vector<Scalar> dinv;
    std::vector<index_t> dinvPtr;
  };

  /// Inverts the diagonal node blocks (diagonal inverse for singular blocks)
  bool blockDiagonalInverse(Level &L) const {
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

  /// Greedy aggregation of the nodes on the (symmetrized) graph of strong
  /// couplings; returns the number of aggregates
  index_t aggregate(const Level &L, std::vector<index_t> &agg) const {
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
    for (int pass = 0; pass != 2; ++pass) {
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
          } else if (J != I && norm > settings.theta * math::sqrt(diag[I] * diag[J])) {
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

  /// Prolongation: every component is interpolated by a constant on each
  /// aggregate (columns normalized); \a C receives the nodes and components
  /// of the coarse level
  RowMatrix prolongation(const Level &L, const std::vector<index_t> &agg,
                         const index_t numAgg, Level &C) const {
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
    if (!settings.smoothProlongation) { return P; }

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
    for (int k = 0; k != 15; ++k) {
      w = DA * v;
      rho = w.norm() / v.norm();
      if (!(rho > 0)) { return P; }
      v = w / w.norm();
    }
    const RowMatrix smoothed = P - (Scalar(4) / (3 * rho)) * RowMatrix(DA * P);
    return smoothed;
  }

  /// One sweep of the nodal block Gauss-Seidel method
  void smooth(const Level &L, const Vector &b, Vector &x, const bool forward) const {
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

  /// V-cycle on level l with zero initial guess
  void cycle(const size_t l, const Vector &b, Vector &x) const {
    const Level &L = m_levels[l];
    if (l + 1 == m_levels.size()) {
      if (m_coarseDirect) {
        x = m_coarseLU.solve(b);
      } else {  // the coarsening stopped early: smoothing only
        x.setZero(b.size());
        for (int k = 0; k != 2 * settings.postSweeps; ++k) {
          smooth(L, b, x, true);
          smooth(L, b, x, false);
        }
      }
      return;
    }
    x.setZero(b.size());
    for (int k = 0; k != settings.preSweeps; ++k) { smooth(L, b, x, true); }
    const Vector coarseRhs = L.R * (b - L.A * x);
    Vector coarseSol;
    cycle(l + 1, coarseRhs, coarseSol);
    x += L.P * coarseSol;
    for (int k = 0; k != settings.postSweeps; ++k) { smooth(L, b, x, false); }
  }

  std::vector<index_t> m_perm, m_inv, m_nodePtr0, m_comp0;
  std::vector<Level> m_levels;
  gsEigen::PartialPivLU<Dense> m_coarseLU;
  bool m_coarseDirect = true;
  gsEigen::ComputationInfo m_info = gsEigen::Success;
};

/// @brief Block Gauss-Seidel preconditioner: the block lower triangle of a
/// matrix, for a given partition of the unknowns into blocks (e.g.
/// coordinate components), applied by forward block substitution. The
/// diagonal blocks are solved by sparse LU factorizations or, if
/// \a amgBlocks, approximated by one cycle of the algebraic multigrid method
/// NodalAMG (scalar nodes), which is much cheaper to set up for large 3D
/// problems.
template<typename Scalar = double>
class BlockGaussSeidel : public gsLinearOperator<Scalar> {
 public:
  typedef gsEigen::SparseMatrix<Scalar, 0, index_t> SparseMatrix;

  explicit BlockGaussSeidel(const bool amgBlocks = false)
      : m_amgBlocks(amgBlocks) {}

  /// Sets the block of every unknown
  void setBlocks(const std::vector<index_t> &blocks) {
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

  template<typename MatType>
  BlockGaussSeidel &compute(const MatType &mat) {
    GISMO_ASSERT(m_blockOf.size() == size_t(mat.rows()),
                 "The blocks of the unknowns are not set.");
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

  /// Forward block substitution
  template<typename Rhs>
  gsVector<Scalar> solve(const gsEigen::MatrixBase<Rhs> &rhs) const {
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

  gsEigen::ComputationInfo info() const { return m_info; }

  /// gsLinearOperator interface (for gismo's iterative solvers)
  void apply(const gsMatrix<Scalar> &input, gsMatrix<Scalar> &x) const {
    x.resize(input.rows(), input.cols());
    for (index_t j = 0; j != input.cols(); ++j) { x.col(j) = solve(input.col(j)); }
  }
  index_t rows() const { return m_blockOf.size(); }
  index_t cols() const { return m_blockOf.size(); }

 private:
  std::vector<index_t> m_blockOf, m_local;
  std::vector<std::vector<index_t> > m_members;
  std::vector<std::vector<SparseMatrix> > m_lower;
  std::vector<memory::shared_ptr<gsEigen::SparseLU<SparseMatrix> > > m_lu;
  std::vector<memory::shared_ptr<NodalAMG<Scalar> > > m_amg;
  bool m_amgBlocks;
  gsEigen::ComputationInfo m_info = gsEigen::Success;
};

/// @brief Anderson acceleration solver and its (preconditioned) variants
///
/// \ingroup Solver
template<typename Scalar = double>
class AndersonAcceleration {
 public:
  typedef std::function<gsSparseMatrix<Scalar>(gsVector<Scalar> const &)> Jacobian_t;
  typedef std::function<gsVector<Scalar>(gsVector<Scalar> const &)> Residual_t;
  /// Step filter: given the current iterate and a proposed new iterate, it
  /// may shorten the step (by modifying the proposal in place). Returns false
  /// if no admissible step exists.
  typedef std::function<bool(gsVector<Scalar> const &, gsVector<Scalar> &)> StepFilter_t;

  explicit AndersonAcceleration(const preAAParam<Scalar> &param) :
      m_param(param) {
    checkParam();
  }

  /// perform Anderson acceleration iteration
  const gsVector<Scalar> &compute(const gsVector<Scalar> &u0, const Residual_t &F,
                         const Jacobian_t &Jacobian) {
    // print the current solver information
    if (m_printInfo) { printSolverInfo();}

    // initialize the solver and preallocate memory
    initialize(u0);

    // Iteration 0
    updateG(F, Jacobian);

    // check if the initial guess is the solution
    if (m_currResidualNorm < m_param.epsilon) {
      if (m_printInfo)
        printf("You are lucky, the initial guess is exactly the solution.\n\n\n");
      return u0;
    }
    m_bestNorm = m_currResidualNorm;
    m_bestSolution = m_solution;
    if (!acceptStep(m_currentG)) {
      // no admissible step from the initial guess
      m_stalled = true;
      return m_solution;
    }

    // print iteration information
    if (m_printInfo) {
      printf("Iter.      ||F(x)|| \n");
      printIterationInfo();
    }

    // Start iteration
    startIteration(F, Jacobian);

    return m_solution;
  }

  inline void enableIterationInfoPrinting() { m_printInfo = true; }
  inline void disableIterationInfoPrinting() { m_printInfo = false; }

  /// Sets a filter that is applied to every step (e.g. a line search)
  inline void setStepFilter(const StepFilter_t &filter) { m_stepFilter = filter; }

  /// True if the residual norm dropped below the tolerance
  inline bool converged() const { return m_currResidualNorm < m_param.epsilon; }

  /// True if the iteration stopped since the step filter rejected all steps
  inline bool stalled() const { return m_stalled; }

  /// Number of iterations performed
  inline int iterations() const { return m_iter; }

  /// Norm of the last evaluated residual
  inline Scalar residualNorm() const { return m_currResidualNorm; }

  /// Number of preconditioner (Jacobian) updates
  inline int preconditionerUpdates() const { return m_numPreconditionerUpdates; }

  /// Number of restarts after a rejected accelerated step
  inline int restarts() const { return m_numRestarts; }

  /// Wall time (s) spent in assembling and factorizing the preconditioner
  inline double preconditionerTime() const { return m_timePreconditioner; }

  /// Wall time (s) spent in evaluating the residual
  inline double residualTime() const { return m_timeResidual; }

  /// Wall time (s) spent in solving with the preconditioner
  inline double solveTime() const { return m_timeSolve; }

  /// Total number of iterations of the iterative linear solver (GMRES)
  inline int innerIterations() const { return m_numInnerIterations; }

  /// Number of fall-backs from the iterative linear solver to sparse LU
  inline int innerFallbacks() const { return m_numFallbacks; }

  /// Number of iterates rejected by the safeguard
  inline int rejections() const { return m_numRejections; }

 private:
  inline void initialize(const gsVector<Scalar> &u0) {
    m_solution = u0;
    m_dim = m_solution.size();
    m_iter = 0;
    m_sinceRestart = 0;
    m_columnIndex = 0;
    m_stalled = false;
    m_numPreconditionerUpdates = 0;
    m_numRestarts = 0;
    m_numInnerIterations = 0;
    m_numFallbacks = 0;
    m_iterativeFailed = false;
    m_numRejections = 0;
    m_damping = 1;
    m_timePreconditioner = m_timeResidual = m_timeSolve = 0;
    m_blockPreconditioner.reset();
    m_amgPreconditioner.reset();
    if (m_param.linearSolver >= 1 && m_param.linearSolver <= 3) {
      GISMO_ENSURE(m_param.blocks.size() == size_t(m_dim),
                   "preAAParam::blocks must give the block of every unknown.");
    }
    if (1 == m_param.linearSolver || 3 == m_param.linearSolver) {
      m_blockPreconditioner.reset(
          new BlockGaussSeidel<Scalar>(3 == m_param.linearSolver));
      m_blockPreconditioner->setBlocks(m_param.blocks);
    } else if (2 == m_param.linearSolver) {
      GISMO_ENSURE(m_param.nodes.size() == size_t(m_dim),
                   "preAAParam::nodes must give the node of every unknown.");
      m_amgPreconditioner.reset(new NodalAMG<Scalar>());
      m_amgPreconditioner->setNodes(m_param.nodes, m_param.blocks);
    }

    m_currentF.resize(m_dim);
    m_prevdG.resize(m_dim, m_param.m);
    m_prevdF.resize(m_dim, m_param.m);

    m_normalEquationMatrix.resize(m_param.m, m_param.m);
    m_alpha.resize(m_param.m);
    m_scaledF.resize(m_param.m);

    if (m_param.usePreconditioning) { m_preconditioner.resize(m_dim, m_dim); }
  }

  inline void checkParam() {
    m_param.check_param();
  }

  inline void printSolverInfo() {
    printf("\nAnderson Acceleration SOLVER: parameter settings... \n");
    printf("depth                       =     %d\n", m_param.m);
    printf("use preconditioner          =     %d\n",
           m_param.usePreconditioning);
    printf("update preconditioner step  =     %d\n\n",
           m_param.updatePreconditionerStep);
  }

  /// Moves to the proposed iterate, subject to the step filter (if any).
  /// Returns false if the proposal is rejected; the iterate is then unchanged.
  inline bool acceptStep(const gsVector<Scalar> &proposal) {
    if (!m_stepFilter) {
      m_solution = proposal;
      return true;
    }
    gsVector<Scalar> next = proposal;
    if (!next.allFinite() || !m_stepFilter(m_solution, next)) { return false; }
    m_solution = next;
    return true;
  }

  inline void startIteration(const Residual_t &F, const Jacobian_t &Jacobian) {
    m_iter = 1;
    m_sinceRestart = 1;
    if (m_param.m == 0) {
      performPicardIteration(F, Jacobian);
    } else {
      performAAIteration(F, Jacobian);
    }
  }

  inline void performPicardIteration(const Residual_t &F,
                                     const Jacobian_t &Jacobian) {
    while (m_iter < m_param.max_iterations
        && m_currResidualNorm > m_param.epsilon) {
      updateG(F, Jacobian);

      // update the solution
      if (!acceptStep(m_currentG)) { m_stalled = true; break; }

      printIterationInfo();

      ++m_iter;
      ++m_sinceRestart;
      trackIterationInfo();
    }
  }

  inline void performAAIteration(const Residual_t &F,
                                 const Jacobian_t &Jacobian) {
    m_prevdF.col(0) = -m_currentF;
    m_prevdG.col(0) = -m_currentG;

    while (m_iter < m_param.max_iterations &&
        m_currResidualNorm > m_param.epsilon) {
      updateG(F, Jacobian);

      printIterationInfo();

      if (m_param.safeguardFactor > 0) {
        if (!(m_currResidualNorm <= m_param.safeguardFactor * m_bestNorm)) {
          // Divergence: back to the best iterate, restart there with a fresh
          // preconditioner and a damped plain step
          ++m_numRejections;
          m_damping *= 0.5;
          m_solution = m_bestSolution;
          if (m_damping < Scalar(1) / 1024) {
            m_currResidualNorm = m_bestNorm;
            m_stalled = true;
            break;
          }
          ++m_iter;
          m_sinceRestart = 0;
          trackIterationInfo();
          continue;
        }
        if (m_currResidualNorm < m_bestNorm) {
          m_bestNorm = m_currResidualNorm;
          m_bestSolution = m_solution;
          m_damping = std::min<Scalar>(1, 2 * m_damping);
        }
      }

      if (m_sinceRestart == 0) {
        // restart: plain (preconditioned, possibly damped) step with a fresh
        // history
        if (!acceptStep(m_damping < 1 ? gsVector<Scalar>(m_solution + m_damping * m_currentF)
                                      : m_currentG)) { m_stalled = true; break; }
        m_columnIndex = 0;
        m_prevdF.col(0) = -m_currentF;
        m_prevdG.col(0) = -m_currentG;
      } else {
        updatePrevdFAndPrevdG();

        // update alpha and solution (compute normal equation)
        updateAlpha();
        if (!updateSolution()) {
          // The accelerated step was rejected by the step filter: restart
          // from the current iterate with a fresh preconditioner.
          ++m_iter;
          ++m_numRestarts;
          m_sinceRestart = 0;
          trackIterationInfo();
          continue;
        }

        // update the column indices
        updateColumnIndices();
      }

      ++m_iter;
      ++m_sinceRestart;
      trackIterationInfo();
    }
  }

  /// update fixed point function
  inline void updateG(const Residual_t &F, const Jacobian_t &Jacobian) {
    if (m_param.usePreconditioning) {
      updateGWithPreconditioning(F, Jacobian);
    } else {
      updateGWithoutPreconditioning(F);
    }
    m_currentG = m_currentF + m_solution;
  }

  /// update fixed point function with preconditioning
  inline void updateGWithPreconditioning(const Residual_t &F,
                                         const Jacobian_t &Jacobian) {
    gsStopwatch timer;
    if (!(m_sinceRestart % m_param.updatePreconditionerStep)) {
      m_preconditioner = Jacobian(m_solution);
      m_useLU = !(m_blockPreconditioner || m_amgPreconditioner)
                || m_iterativeFailed;
      if (!m_useLU) {
        memory::shared_ptr<gsLinearOperator<Scalar> > preconditioner;
        if (m_blockPreconditioner) {
          m_blockPreconditioner->compute(m_preconditioner);
          m_useLU = (m_blockPreconditioner->info() != gsEigen::Success);
          preconditioner = m_blockPreconditioner;
        } else {
          m_amgPreconditioner->compute(m_preconditioner);
          m_useLU = (m_amgPreconditioner->info() != gsEigen::Success);
          preconditioner = m_amgPreconditioner;
        }
        if (m_useLU) {
          ++m_numFallbacks;
          m_iterativeFailed = true;
        } else {
          m_iterativeSolver.reset(
              new gsGMRes<Scalar>(m_preconditioner, preconditioner));
          m_iterativeSolver->setTolerance(m_param.innerTolerance);
          m_iterativeSolver->setMaxIterations(m_param.innerMaxIterations);
        }
      }
      if (m_useLU) { m_linearSolverPreconditioning.compute(m_preconditioner); }
      ++m_numPreconditionerUpdates;
      m_timePreconditioner += timer.stop();
    }
    if (m_useLU && m_linearSolverPreconditioning.info() != gsEigen::Success) {
      // The preconditioner could not be factorized: no step, and a
      // non-finite residual norm which terminates the iteration.
      m_currentF.setZero(m_dim);
      m_currResidualNorm = std::numeric_limits<Scalar>::quiet_NaN();
      return;
    }
    timer.restart();
    gsVector<Scalar> residual = F(m_solution);
    m_currResidualNorm = residual.norm();
    m_timeResidual += timer.stop();

    timer.restart();
    if (!m_useLU) {
      gsMatrix<Scalar> solution;
      solution.setZero(m_dim, 1);
      m_iterativeSolver->solve(residual, solution);
      m_numInnerIterations += m_iterativeSolver->iterations();
      // gismo's GMRES is left preconditioned: check the true residual
      const Scalar trueResidual = (m_preconditioner * solution - residual).norm();
      m_currentF = -solution.col(0);
      if (!(trueResidual <= m_param.innerAcceptance * residual.norm())
          || !m_currentF.allFinite()) {
        // The iterative solver did not converge: sparse LU from now on
        ++m_numFallbacks;
        m_useLU = true;
        m_iterativeFailed = true;
        gsStopwatch factorizationTimer;
        m_linearSolverPreconditioning.compute(m_preconditioner);
        m_timePreconditioner += factorizationTimer.stop();
        if (m_linearSolverPreconditioning.info() != gsEigen::Success) {
          m_currentF.setZero(m_dim);
          m_currResidualNorm = std::numeric_limits<Scalar>::quiet_NaN();
          return;
        }
      }
    }
    if (m_useLU) {
      m_currentF = -m_linearSolverPreconditioning.solve(residual);
    }
    m_timeSolve += timer.stop();
  }

  /// update fixed point function without preconditioning
  inline void updateGWithoutPreconditioning(const Residual_t &F) {
    m_currentF = F(m_solution);
    m_currResidualNorm = m_currentF.norm();
  }

  inline void printIterationInfo() {
    if (m_printInfo) {
        printf(" %d         %.4e\n", m_iter, cast<Scalar,double>( m_currResidualNorm ));
    }
  }

  inline void updatePrevdFAndPrevdG() {
    m_prevdF.col(m_columnIndex) += m_currentF;
    m_prevdG.col(m_columnIndex) += m_currentG;

    // scale previous dF for better numerical stability
    Scalar scale = std::max(EPSILON, m_prevdF.col(m_columnIndex).norm());
    m_scaledF(m_columnIndex) = scale;
    m_prevdF.col(m_columnIndex) /= scale;
  }

  /// update the coefficients \f$ alpha_i \f$ by solving a Least-Square problem
  inline void updateAlpha() {
    // compute m_mk
    m_mk = std::min(m_param.m, m_sinceRestart);

    if (m_mk == 1) {
      m_alpha(0) = 0;
      Scalar dF_squaredNorm = m_prevdF.col(m_columnIndex).squaredNorm();
      m_normalEquationMatrix(0, 0) = dF_squaredNorm;
      Scalar dF_norm = math::sqrt(dF_squaredNorm);

      // For better numerical stability
      if (dF_norm > EPSILON) {
        // compute alpha = (dF * F) / (dF * dF)
        m_alpha(0) = (m_prevdF.col(m_columnIndex) / dF_norm).dot(
            m_currentF / dF_norm);
      }
    } else {
      // Update the normal equation matrix
      // for the column and row corresponding to the new dF column.
      // note: only one column and one row are needed to be updated.
      gsVector<Scalar> new_inner_prod = (m_prevdF.col(m_columnIndex).transpose()
          * m_prevdF.block(0, 0, m_dim, m_mk)).transpose();
      m_normalEquationMatrix.block(m_columnIndex, 0, 1, m_mk) =
          new_inner_prod.transpose();
      m_normalEquationMatrix.block(0, m_columnIndex, m_mk, 1) =
          new_inner_prod;

      // Solve normal equation: A^{T} A x = A^{T} b
      m_linearSolver.compute(m_normalEquationMatrix.block(0, 0, m_mk, m_mk));
      m_alpha.head(m_mk) = m_linearSolver.solve(
          m_prevdF.block(0, 0, m_dim, m_mk).transpose() * m_currentF);
    }
  }

  /// update the solution; returns false if the step filter rejects the step
  inline bool updateSolution() {
    // Update the current solution (x) using the rescaled alpha
    if (!m_stepFilter) {
      m_solution = m_currentG - m_prevdG.block(0, 0, m_dim, m_mk) *
          ((m_alpha.head(m_mk).array()
              / m_scaledF.head(m_mk).array()).matrix());
      return true;
    }
    const gsVector<Scalar> proposal = m_currentG
        - m_prevdG.block(0, 0, m_dim, m_mk) * ((m_alpha.head(m_mk).array()
            / m_scaledF.head(m_mk).array()).matrix());
    return acceptStep(proposal);
  }

  void updateColumnIndices() {
    m_columnIndex = (m_columnIndex + 1) % m_param.m;
    m_prevdF.col(m_columnIndex) = -m_currentF;
    m_prevdG.col(m_columnIndex) = -m_currentG;
  }

  inline void trackIterationInfo() {
    if (m_trackResidualNorm) { m_residualList.push_back(m_currResidualNorm); }
  }

  // Linear solver for the Least-Square problem
  gsEigen::FullPivLU<gsMatrix<Scalar>> m_linearSolver;
  // Parameters
  const preAAParam<Scalar> m_param;
  int m_dim = -1;
  // Iteration
  int m_iter = 0;
  // Iterations since the last (re)start; equals m_iter if no restart occurred
  int m_sinceRestart = 0;
  int m_mk = 0;
  int m_columnIndex = 0;
  // Residual
  Scalar m_currResidualNorm;
  // Solution
  gsVector<Scalar> m_solution;
  gsVector<Scalar> m_currentF;
  gsVector<Scalar> m_currentG;
  // Anderson extrapolation
  gsMatrix<Scalar> m_prevdG;
  gsMatrix<Scalar> m_prevdF;
  gsMatrix<Scalar> m_normalEquationMatrix;
  gsVector<Scalar> m_alpha;
  gsVector<Scalar> m_scaledF;
  // Preconditioning
  gsSparseMatrix<Scalar> m_preconditioner;
  gsEigen::SparseLU<gsSparseMatrix<Scalar>> m_linearSolverPreconditioning;
  memory::shared_ptr<BlockGaussSeidel<Scalar> > m_blockPreconditioner;
  memory::shared_ptr<NodalAMG<Scalar> > m_amgPreconditioner;
  memory::unique_ptr<gsGMRes<Scalar> > m_iterativeSolver;
  bool m_useLU = true;
  // the iterative solver failed once: sparse LU for the rest of the solve
  bool m_iterativeFailed = false;
  // Optional step filter (line search) and its status
  StepFilter_t m_stepFilter;
  bool m_stalled = false;
  // Statistics
  int m_numPreconditionerUpdates = 0;
  int m_numRestarts = 0;
  int m_numInnerIterations = 0;
  int m_numFallbacks = 0;
  // Safeguard
  int m_numRejections = 0;
  Scalar m_bestNorm = 0, m_damping = 1;
  gsVector<Scalar> m_bestSolution;
  double m_timePreconditioner = 0, m_timeResidual = 0, m_timeSolve = 0;
  // Print information
  bool m_printInfo = true;

  // track residual norm history
  bool m_trackResidualNorm = false;
  std::vector<Scalar> m_residualList;

  const Scalar EPSILON = 1e-14;
};

}

}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBarrierCore.hpp)
#endif
