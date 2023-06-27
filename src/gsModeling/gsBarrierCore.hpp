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

namespace gismo {

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
  options.addInt("overInt", "Apply over-integration?", 0);

  // Preconditioner type for AA: 0: NO preconditioning, 1: Full Jacobian preconditioner, 2: Diagonal Jacobian preconditioner, 3: Diagonal Block Jacobian preconditioner
  options.addInt("AAPreconditionType", "Preconditioner type for AA", 0);

  // Update the preconditioner every N_update steps
  options.addInt("N_update", "update preconditioner every N_update steps", 10);

  // window size m for our preconditioned AA solver
  options.addInt("AAwindowsize", "window size for preconditioned AA solver", 5);

  // need improve the parameterization quality for PDE-based method?
  options.addSwitch("needPDEH1", "improve quality by H1 discrezation?", true);

  return options;
}

/// Compute the area of the computational domain
template<short_t d, typename T>
T gsBarrierCore<d, T>::computeArea(const gsMultiPatch<T> &mp) {
  return computeAreaInterior(mp);
}

/// Compute the area of a multi-patch representing computational domain
template<short_t d, typename T>
T gsBarrierCore<d, T>::computeAreaInterior(const gsMultiPatch<T> &multiPatch) {
  // Creating a multi-basis from the provided multi-patch
  gsMultiBasis<T> multiBasis(multiPatch);

  // Initializing an expression evaluator
  gsExprEvaluator<T> evaluator;

  // Setting integration elements
  evaluator.setIntegrationElements(multiBasis);

  // Getting the geometry map
  geometryMap geometry = evaluator.getMap(multiPatch);

  // Computing the area by integrating over the determinant of the Jacobian of the geometry map
  T area = evaluator.integral(jac(geometry).det());

  return area;
}

template<short_t d, typename T>
T gsBarrierCore<d, T>::computeAreaBoundary(const gsMultiPatch<T> &mp) {
  GISMO_ERROR("Not implemented");
  return -1;
}

enum class ParamMethod {
  BarrierPatch,
  PenaltyPatch,
  PenaltyPatch2,
  PDEPatch,
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
      result = computePenaltyPatch(mp, mapper, options);
      break;
    }
    case ParamMethod::PenaltyPatch2: {
      result = computePenaltyPatch2(mp, mapper, options);
      break;
    }
    case ParamMethod::PDEPatch: {
      result = computePDEPatch(mp, mapper, options);
      break;
    }
    case ParamMethod::VariationalHarmonicPatch: {
      result = computeVHPatch(mp, mapper, options);
      break;
    }
    case ParamMethod::BarrierPatch:
    default: {
      if (method != ParamMethod::BarrierPatch) {
        gsWarn << "Invalid ParamMethod value. Defaulting to BarrierPatch.\n";
      }
      result = computeBarrierPatch(mp, mapper, options);
      break;
    }
  }
  return result;
}

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
  objPenaltyPt.setEps(options.getReal("ff_Delta") * scaledArea);

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
gsMultiPatch<T> gsBarrierCore<d,
                              T>::computeBarrierPatch(const gsMultiPatch<T> &mp,
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
  optFoldoverFree.options().setReal("MinGradientLength",
                                    options.askReal("ff_MinGradientLength",
                                                    1e-12));
  optFoldoverFree.options().setReal("MinStepLength",
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

  if (Efoldover > EPSILON) {
    throw std::runtime_error(
        "Maximum iterations reached. The foldover-energy value is " +
            std::to_string(Efoldover) +
            ". This suggests there may be issues with the input data.");
  }
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
  m_assembler.setIntegrationElements(m_mb);
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
  static const space space1 = m_assembler.getSpace(m_mb, d);
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
  m_assembler.setIntegrationElements(m_mb);
  m_evaluator = gsExprEvaluator<T>(m_assembler);
//  defaultOptions();
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
//        gsMultiPatch<T> mp = m_mp;
  convertFreeVectorToMultiPatch<T>(u, m_mapper, m_mp);

  geometryMap G = m_evaluator.getMap(m_mp);

  if (m_evaluator.min(jac(G).det()) < 0) {
    return std::numeric_limits<T>::max();
  } else {
    auto Euniform = pow(jac(G).det(), 2);
    auto Ewinslow = 0.5 * (jac(G).sqNorm() * jac(G).inv().sqNorm());

//            // another objective function term - Jiao et al. 2011
//            auto Ewinslow = jac(G).sqNorm() / pow(jac(G).det(), 2.0 / 3.0);

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

//        // gradient of another objective function term
//        auto Ewinslow = jac(G).sqNorm() / pow(jac(G).det(), 2.0 / 3.0);
//        auto derEwinslow =
//                2.0 * frprod2(space1, jac(G)) / pow(jac(G).det(), 2.0 / 3.0) -
//                2.0 * Ewinslow.val() * derJacDet / jac(G).det() / 3.0;
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
  m_assembler.setIntegrationElements(m_mb);
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
  m_assembler.setIntegrationElements(m_mb);
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
//  return evalObj_impl<d>(u);
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
  m_assembler.setIntegrationElements(m_mb);
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
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exponent());
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());

  // Calculation of Ewinslow and Euniform
  auto Ewinslow = jac(G).sqNorm() / chi;
  auto Euniform = chi + 1.0 / chi;

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
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exponent());
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());

  // Compute Ewinslow
  auto Ewinslow = 0.5 * (jac(G).sqNorm() * jac(G).inv().sqNorm()) *
      pow(jac(G).det(), 2.0) / pow(chi, 2.0);

  // Compute Euniform
  auto Euniform = chi + 1.0 / chi;

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
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exponent());
  auto chi =
      ternary(eps.val() - jac(G).det(), chiPPart.val(), jac(G).det().val());
  auto chip = ternary(eps.val() - jac(G).det(), chiPPart.val(), unit.val());

  // Define Ewinslow and its derivative
  auto Ewinslow = jac(G).sqNorm() / chi;
  auto derEwinslow = 1.0 / chi * (2.0 * frprod2(spaceMapper, jac(G)) -
      Ewinslow.val() * chip * derJacDet);

  // Define the derivative of Euniform
  auto derEuniform = (chip - chip / pow(chi, 2)) * derJacDet;

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
  auto chiPPart = eps * ((jac(G).det() - eps.val()).exponent());

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
  auto derEuniform = (chip - chip / pow(chi, 2)) * derJacDet;

  // Initialize and assemble the system
  m_assembler.initSystem();
  m_assembler.assemble(m_lambda1 * derEwinslow + m_lambda2 * derEuniform);

  // Compute and return the result
  result = gsAsVector<T>(const_cast<T *>(m_assembler.rhs().data()),
                         m_assembler.rhs().rows());
  return EXIT_SUCCESS;
}

template<short_t d, typename T>
gsMultiPatch<T>
gsBarrierCore<d, T>::computePDEPatch(const gsMultiPatch<T> &mp,
                                     const gismo::gsDofMapper &mapper,
                                     const gismo::gsOptionList &options) {

  // get initial guess vector
  gsVector<T> initialGuessVector = convertMultiPatchToFreeVector(mp, mapper);

//  std::string solverSetting = "Parameters setting:\n";
//  solverSetting += "\t\t\t Window size = " +
//      std::to_string(options.getInt("AAwindowsize")) + "\n";
//  solverSetting += "\t\t\t Update preconditioner every " +
//      std::to_string(options.getInt("N_update")) + " steps \n";
//  solverSetting += "\t\t\t Preconditioner type: " +
//      std::to_string(options.getInt("AAPreconditionType")) + "\n";
//  solverSetting += "\t\t\t Need improve by H1?  " +
//      std::to_string(options.getSwitch("needPDEH1")) + "\n";
//  verboseLog(solverSetting, options.askInt("Verbose", 0));

  verboseLog("PDE-based parameterization construction ...\n",
             options.askInt("Verbose", 0));

  gsExprEvaluator<T> evaluator;
  gsExprAssembler<T> assembler;

  gsMultiBasis<T> mb(mp);
  assembler.setIntegrationElements(mb);

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
      Jacobian = [&Residual, &assembler, &mapper, &mpSubstitute, &space1](
          gsVector<T> const &x) {
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
      Jacobian = [&Residual, &assembler, &mapper, &mpSubstitute, &space1](
          gsVector<T> const &x) {
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
      Jacobian = [&assembler, &mapper, &mpSubstitute, &space1](
          gsVector<T> const &x) {
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
  gsVector<T> solVector = AASolver.compute(initialGuessVector,
                                           Residual, Jacobian);

//  int m = options.askInt("AAwindowsize", 5);
//  // TODO: use preAApp and remove these resHist etc.
//  AndersonAcceleration<T> solver(m);
//  std::vector<int> iterHist;
//  std::vector<double> resHist, timeHist;
//  gsVector<T> solVector = solver.computePrecond(initialGuessVector,
//                                                Residual, Jacobian, iterHist,
//                                                resHist, timeHist);

  if (options.askSwitch("needPDEH1", true)) {
    verboseLog("\nStart parameterization improvement by H1 discrezation...",
               options.askInt("Verbose", 0));
    Residual = [&assembler, &mapper, &mpSubstitute, &space1](
        gsVector<T> const &x) {
      // H1 discretization
      convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);
      geometryMap G = assembler.getMap(mpSubstitute);

      auto invJacMat = jac(G).inv();
      auto jacU = jac(space1) * invJacMat;
      auto nonlinearSystem = jacU % invJacMat;

      assembler.initSystem();
      assembler.assemble(nonlinearSystem);
      return assembler.rhs();
    };

    preconditionerType = options.askInt("AAPreconditionType", 0);
    switch (preconditionerType) {
      case 1:
        Jacobian = [&assembler, &mapper, &mpSubstitute, &mb, &space1](
            gsVector<T> const &x) {
          // diagonal block Jacobian matrix as a preconditioner
          convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);
          geometryMap G = assembler.getMap(mpSubstitute);

          auto jacMat = jacScaledLxH1DiagBlock(space1, G);

          assembler.initSystem();
          assembler.assemble(jacMat);

          return assembler.matrix();
        };
        break;

      default:
        Jacobian = [&assembler, &mapper, &mpSubstitute,
            &space1](gsVector<T> const &x) {
          // analytical Jacobian matrix as a preconditioner
          convertFreeVectorToMultiPatch<T>(x, mapper, mpSubstitute);
          geometryMap G = assembler.getMap(mpSubstitute);

          auto jacMat = jacScaledLxH1(space1, G);

          assembler.initSystem();
          assembler.assemble(jacMat);

          return assembler.matrix();
        };
    }

    preAApp::AndersonAcceleration<T> AASolver2(param);
    solVector = AASolver2.compute(solVector, Residual, Jacobian);
  }

  gsMultiPatch<T> result = mp;
  convertFreeVectorToMultiPatch<T>(solVector, mapper, result);

  verboseLog("Finished!", options.askInt("Verbose", 0));

  return result;
}

namespace expr {

template<typename E1, typename E2>
class frprod2_expr;
template<typename E1, typename E2>
class frprod3_expr;
template<class E0, class E1, class E2>
class ternary_expr; // ternary expression
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

/*

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
      : _u(u), _v(v) {
    // gsInfo << "expression is space ? "<<E1::Space <<"\n"; _u.print(gsInfo);
    // GISMO_ASSERT(_u.rows() == _v.rows(),
    //              "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
    // GISMO_ASSERT(_u.cols() == _v.cols(),
    //              "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
  }

  const gsMatrix<Scalar> &eval(const index_t k) const //todo: specialize for nb==1
  {
    auto A = _v.eval(k);
    b = _u.eval(k);
//        gsDebugVar(b);

    res.noalias() = b * A;
    return res;
  }

  index_t rows() const { return 1; }
  index_t cols() const { return 1; }

  void parse(gsExprHelper<Scalar> &evList) const {
    _v.parse(evList);
//        evList.add(_u);
//        _u.data().flags |= NEED_GRAD;
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

template<class E0, class E1, class E2>
class ternary_expr : public _expr<ternary_expr<E0, E1, E2> > {
  typename E0::Nested_t _u;
  typename E1::Nested_t _v;
  typename E2::Nested_t _w;
 public:
  typedef typename E1::Scalar Scalar;

  explicit ternary_expr(_expr<E0> const &u,
                        _expr<E1> const &v,
                        _expr<E2> const &w)
      :
      _u(u),
      _v(v),
      _w(w) {
    GISMO_ASSERT(E0::ScalarValued, "Condition must be scalar valued");
    GISMO_ASSERT((int) E1::ScalarValued == (int) E2::ScalarValued,
                 "Both v and w must be scalar valued (or not).");
    GISMO_ASSERT((int) E1::ColBlocks == (int) E2::ColBlocks,
                 "Both v and w must be colblocks (or not).");
    GISMO_ASSERT((int) E1::Space == (int) E2::Space,
                 "Both v and w must be space (or not), but E1::Space = "
                     << E1::Space << " and E2::Space = " << E2::Space);
    GISMO_ASSERT(_v.rows() == _w.rows(),
                 "Rows of v and w differ. _v.rows() = " << _v.rows()
                                                        << ", _w.rows() = "
                                                        << _w.rows());
    GISMO_ASSERT(_v.cols() == _w.cols(),
                 "Columns of v and w differ. _v.cols() = " << _v.cols()
                                                           << ", _w.cols() = "
                                                           << _w.cols());
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
  // { res = eval_impl(_u,_v,_w,k); return  res;}

  index_t rows() const { return _v.rows(); }
  index_t cols() const { return _v.cols(); }
  void parse(gsExprHelper<Scalar> &evList) const {
    _u.parse(evList);
    _v.parse(evList);
    _w.parse(evList);
  }

  const gsFeSpace<Scalar> &rowVar() const { return _v.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _v.colVar(); }

// private:
//     template<class U, class V, class W> static inline
//     typename util::enable_if<U::ScalarValued && V::ScalarValued,AutoReturn_t>::type
//     eval_impl(const U &u, const V & v, const W & w, const index_t k)
//     {
//         gsMatrix<Scalar> res(1,1);
// //        bool test = u.eval(k) > 0;
//         res<<(u.eval(k) > 0 ?  v.eval(k) : w.eval(k));
//         return res;
//     }

//     template<class U, class V, class W> static inline
//     typename util::enable_if<U::ScalarValued && !V::ScalarValued,AutoReturn_t>::type
//     eval_impl(const U &u, const V & v, const W & w, const index_t k)
//     {
//         return u.eval(k) > 0 ? v.eval(k) : w.eval(k);
//     }

//     template<class U, class V, class W> static inline
//     typename util::enable_if<!U::ScalarValued,gsMatrix<Scalar>>::type
//     eval_impl(const U &u, const V & v, const W & w, const index_t k)
//     {
//         GISMO_ERROR("Something went wrong");
//     }
};

/*
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

//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
//        res.setZero();
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
  // TODO: question, what do these parameters mean?
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLx(";
    _u.print(os);
    os << ")";
  }
};

/*
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

//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
  // TODO: question, what do these parameters mean?
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLxDiag(";
    _u.print(os);
    os << ")";
  }
};

/*
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

//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  const gsMatrix<Scalar> &eval(const index_t k) const {
//        basis = _u.data().values[0].col(k);

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
    res.template topLeftCorner(A, A).noalias() =
        _u.data().values[0].col(k) * dLxdx;
    res.template bottomRightCorner(A, A).noalias() =
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
  // TODO: question, what do these parameters mean?
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLxDiagBlock(";
    _u.print(os);
    os << ")";
  }
};

/*
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

//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
    res.template topLeftCorner(N, N).noalias() = jacdLxdx;
    res.template topRightCorner(N, N).noalias() = jacdLxdy;
    res.template bottomLeftCorner(N, N).noalias() = jacdLydx;
    res.template bottomRightCorner(N, N).noalias() = jacdLydy;

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
  // TODO: question, what do these parameters mean?
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLx(";
    _u.print(os);
    os << ")";
  }
};

/*
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

//  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  const gsMatrix<Scalar> &eval(const index_t k) const {
    gsMatrix<Scalar> jacMat = _G.data().values[1].reshapeCol(k,
                                                             _G.data().dim.first,
                                                             _G.data().dim.second);
    gsMatrix<Scalar> invJacMat = jacMat.inverse();

    derivBasis = _u.data().values[1].col(k).transpose();
    derivBasis.blockTransposeInPlace(_u.dim());
    gsMatrix<Scalar> invJacHatG = invJacMat * derivBasis;

    const index_t N = _u.cardinality() / _u.dim(); // _u.data().actives.rows()
    gsMatrix<Scalar> jacdLxdx(N, N);
//        gsMatrix<Scalar> jacdLxdy(N,N);
//        gsMatrix<Scalar> jacdLydx(N,N);
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
//            jacdLydx.col(i) = -temp.transpose()*invJacMat.col(1) - invJacHatG.transpose()*dinvJacdx.col(1);

      // for y-direction
      temp.row(0) = invJacHatG(0, i) * invJacHatG.row(1);
      temp.row(1) = invJacHatG(1, i) * invJacHatG.row(1);
      gsMatrix<Scalar> dinvJacdy(2, 2);
      dinvJacdy.row(0) = invJacHatG(0, i) * invJacMat.row(1);
      dinvJacdy.row(1) = invJacHatG(1, i) * invJacMat.row(1);
//            jacdLxdy.col(i) = -temp.transpose()*invJacMat.col(0) - invJacHatG.transpose()*dinvJacdy.col(0);
      jacdLydy.col(i) = -temp.transpose() * invJacMat.col(1)
          - invJacHatG.transpose() * dinvJacdy.col(1);
    }

    res.resize(_u.cardinality(), _u.cardinality());
    res.setZero();
    res.template topLeftCorner(N, N) = jacdLxdx;
//        res.template topRightCorner(N,N) = jacdLxdy;
//        res.template bottomLeftCorner(N,N) = jacdLydx;
    res.template bottomRightCorner(N, N) = jacdLydy;

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
  // TODO: question, what do these parameters mean?
  index_t cardinality_impl() const { return _u.cardinality_impl(); }

  void print(std::ostream &os) const {
    os << "jacScaledLx(";
    _u.print(os);
    os << ")";
  }
};

/*

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
      : _u(u), _v(v) {
    // gsInfo << "expression is space ? "<<E1::Space <<"\n"; _u.print(gsInfo);
    // GISMO_ASSERT(_u.rows() == _v.rows(),
    //              "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
    // GISMO_ASSERT(_u.cols() == _v.cols(),
    //              "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
  }

  const gsMatrix<Scalar> &eval(const index_t k) const //todo: specialize for nb==1
  {
    auto A = _v.eval(k);
    bGrads = _u.data().values[1].col(k).transpose();
    bGrads.blockTransposeInPlace(_u.dim());
    res.noalias() = A * bGrads;
    res.transposeInPlace();
    res.resize(1, _u.cardinality());
//        res.reshaped(1,_u.cardinality());
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

/// Ternary ternary_expr
template<class E0, class E1, class E2>
EIGEN_STRONG_INLINE
ternary_expr<E0, E1, E2> ternary(const E0 &u,
                                 const E1 &v,
                                 const E2 &w) {
  return ternary_expr<E0, E1, E2>(u, v, w);
}
} // namespace expr

}// namespace gismo
