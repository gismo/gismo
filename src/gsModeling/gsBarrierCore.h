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

#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsAssembler/gsExprAssembler.h>


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
  optimizer.options().setReal("MinGradientLength",
                              options.askReal("qi_MinGradientLength", 1e-4));
  optimizer.options().setReal("MinStepLength",
                              options.askReal("qi_MinStepLength", 1e-4));
  optimizer.options().setInt("Verbose", options.askInt("Verbose", 0));
}
#endif

/// helper function to verbose log
void verboseLog(const std::string &message, const index_t &verbose) {
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
    // clang-format on
  }

  ///
  /// Checking the validity of preconditioned AA parameters.
  /// An `std::invalid_argument` exception will be thrown if some parameter
  /// is invalid.
  ///
  inline void check_param() const {
    if (m < 0)
      throw std::invalid_argument("'m' must be non-negative");
    if (epsilon < 0)
      throw std::invalid_argument("'epsilon' must be non-negative");
    if (epsilon_rel < 0)
      throw std::invalid_argument("'epsilon_rel' must be non-negative");
    if (max_iterations < 0)
      throw std::invalid_argument("'max_iterations' must be non-negative");
//    if (beta < 0 || beta > 1.0)
//      throw std::invalid_argument("'beta' must be between 0 and 1");
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

/// @brief Anderson acceleration solver and its (preconditioned) variants
///
/// \ingroup Solver
template<typename Scalar = double>
class AndersonAcceleration {
 public:
  typedef std::function<gsSparseMatrix<Scalar>(gsVector<Scalar> const &)> Jacobian_t;
  typedef std::function<gsVector<Scalar>(gsVector<Scalar> const &)> Residual_t;

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
      printf("You are lucky, the initial guess is exactly the solution.\n\n\n");
      return u0;
    } else {
      m_solution = m_currentG;
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

 private:
  inline void initialize(const gsVector<Scalar> &u0) {
    m_solution = u0;
    m_dim = m_solution.size();
    m_iter = 0;
    m_columnIndex = 0;

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

  inline void startIteration(const Residual_t &F, const Jacobian_t &Jacobian) {
    m_iter = 1;
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
      m_solution = m_currentG;

      printIterationInfo();

      ++m_iter;
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

      updatePrevdFAndPrevdG();

      // update alpha and solution (compute normal equation)
      updateAlpha();
      updateSolution();

      // update the column indices
      updateColumnIndices();

      ++m_iter;
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
    if (!(m_iter % m_param.updatePreconditionerStep)) {
      m_preconditioner = Jacobian(m_solution);
      m_linearSolverPreconditioning.compute(m_preconditioner);
    }
    gsVector<Scalar> residual = F(m_solution);
    m_currResidualNorm = residual.norm();

    m_currentF = -m_linearSolverPreconditioning.solve(residual);
  }

  /// update fixed point function without preconditioning
  inline void updateGWithoutPreconditioning(const Residual_t &F) {
    m_currentF = F(m_solution);
    m_currResidualNorm = m_currentF.norm();
  }

  inline void printIterationInfo() {
    if (m_printInfo) {
      printf(" %d         %.4e\n", m_iter, static_cast< double >( m_currResidualNorm ));
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
    m_mk = std::min(m_param.m, m_iter);

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

  /// update the solution
  inline void updateSolution() {
    // Update the current solution (x) using the rescaled alpha
    m_solution = m_currentG - m_prevdG.block(0, 0, m_dim, m_mk) *
        ((m_alpha.head(m_mk).array()
            / m_scaledF.head(m_mk).array()).matrix());
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
#endif (edited) 
