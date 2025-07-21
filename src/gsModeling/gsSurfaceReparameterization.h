/** @file gsSurfaceReparameterization.h

@brief Provides declaration for SurfaceReparameterization.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): Ye Ji
*/

#pragma once

#include <gsNurbs/gsMobiusMap.h>

#include <gsHLBFGS/gsHLBFGS.h>

namespace gismo
{

template <short_t d = 2, typename T = real_t>
class gsObjFuncSurface : public gsOptProblem<T> {
private:
  using geometryMap = typename gsExprAssembler<T>::geometryMap;
  using space = typename gsExprAssembler<T>::space;
  using solution = typename gsExprAssembler<T>::solution;

  using Base = gsOptProblem<T>;

public:
  explicit gsObjFuncSurface(const gsMultiPatch<T> &patches,
                            const gsMobiusMap<T> &mobiusDomain)
      : m_mp(patches), m_MobiusMap(mobiusDomain), m_lambda1(1.0),
        m_lambda2(1.0), m_area(1)
  {
    defaultOptions();

    gsMatrix<T> bbox;
    m_mp.boundingBox(bbox);
    m_mp.patch(0).translate(-bbox.col(0));
    m_mp.patch(0).scale(1 / (bbox.col(1) - bbox.col(0)).array());

    gsComposedGeometry<T> cgeom(m_MobiusMap, m_mp.patch(0));
    m_evaluator.setIntegrationDomain(cgeom.basis().domain());

    // Set the geometry map
    geometryMap G = m_evaluator.getMap(cgeom);
    m_area = m_evaluator.integral(meas(G));

    // Optimizer properties
    m_numDesignVars = m_MobiusMap.nControls();
    m_curDesign = m_MobiusMap.alpha();
  }

  gsOptionList & options() { return m_options; }

  /// @brief see \ref gsOptProblem.h for more details
  T evalObj(const gsAsConstVector<T> &coefsM) const final
  {
    m_MobiusMap.updateGeom(coefsM);

    gsComposedGeometry<T> cgeom(m_MobiusMap, m_mp.patch(0));
    m_evaluator.setIntegrationDomain(cgeom.basis().domain());

    geometryMap G = m_evaluator.getMap(cgeom);
    auto FFF = jac(G).tr() * jac(G);
    auto m_integration = m_lambda1 * (FFF.trace() / meas(G)).val() +
                         m_lambda2 * pow(FFF.det().val(), 2) / pow(m_area, 2);

    return m_evaluator.integral(m_integration);
  }

  /// @brief see \ref gsOptProblem.h for more details
  void gradObj_into(const gsAsConstVector<T> &u,
                    gsAsVector<T> &result) const override
  {
    const index_t n = u.rows();
    gsMatrix<T> uu = u; // Create a copy
    gsAsVector<T> tmp(uu.data(), n);
    gsAsConstVector<T> ctmp(uu.data(), n);

    const T h = static_cast<T>(1e-6);

    // Central finite difference gradient
    for (index_t i = 0; i < n; ++i) {
      tmp[i] += h;
      const T e1 = this->evalObj(ctmp);
      tmp[i] = u[i] - h;
      const T e2 = this->evalObj(ctmp);
      tmp[i] = u[i];
      result[i] = (e1 - e2) / (2 * h);
    }
  }

  /// @brief Default options
  void defaultOptions()
  {
    m_options.addReal("qi_lambda1", "Sets the lambda 1 value", 1.0);
    m_options.addReal("qi_lambda2", "Sets the lambda 2 value", 1.0);
  }

  /// @brief Adds options to the list
  void addOptions(const gsOptionList &options)
  {
    m_options.update(options, gsOptionList::addIfUnknown);
  }

  void setLambda1(T lambda1)
  {
    if (lambda1 < 0.0 || lambda1 > 1.0) {
      gsWarn << "Lambda1 must be between 0 and 1. Setting to default value 1.0.\n";
      lambda1 = 1.0; // Default value if out of bounds
    }
    m_lambda1 = lambda1;
    m_options.setReal("qi_lambda1", lambda1);
  }

  void setLambda2(T lambda2)
  {
    if (lambda2 < 0.0 || lambda2 > 1.0) {
      gsWarn << "Lambda2 must be between 0 and 1. Setting to default value 1.0.\n";
      lambda2 = 1.0; // Default value if out of bounds
    }
    m_lambda2 = lambda2;
    m_options.setReal("qi_lambda2", lambda2);
  }

  void setLambda(T lambda1, T lambda2)
  {
    setLambda1(lambda1);
    setLambda2(lambda2);
  }

  void setLambdaRatio(T ratio)
  {
    if (ratio < 0.0 || ratio > 1.0) {
      gsWarn << "Ratio must be between 0 and 1. Setting to default value 0.5.\n";
      ratio = 0.5; // Default value if out of bounds
    }
    m_lambda1 = ratio;
    m_lambda2 = 1.0 - ratio;
    m_options.setReal("qi_lambda1", ratio);
    m_options.setReal("qi_lambda2", 1.0 - ratio);
  }

  /// @brief Applies the options
  void applyOptions(const gsOptionList &options)
  {
    m_options.update(options, gsOptionList::addIfUnknown);
    m_lambda1 = m_options.getReal("qi_lambda1");
    m_lambda2 = m_options.getReal("qi_lambda2");
    m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
  }

protected:
  using Base::m_curDesign;
  using Base::m_numDesignVars;

  const gsMultiPatch<T> m_mp;
  const gsDofMapper m_mapper;
  const gsMultiBasis<T> m_mb;

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;

  gsOptionList m_options;

  mutable gsMobiusMap<T> m_MobiusMap;
  T m_lambda1, m_lambda2, m_area;
  gsComposedGeometry<T> m_cgeom;
};

template <class T = real_t>
gsMultiPatch<T> approximateWithBSpline(const gsMultiPatch<T> &mp,
                                       const gsMatrix<T> &coefsMobiusIn)
{

  GISMO_ASSERT(mp.geoDim() == 3, "Only 3D geometry is supported.");
  GISMO_ASSERT(mp.dim() == 2,
               "This function only supports 2D parametric domains.");

  gsMultiPatch<T> result;

  gsMatrix<T, 2, 2> alpha = coefsMobiusIn.reshape(2, 2);
  gsMobiusMap<T> mobiusDomain(alpha);

  for (const auto &patch : mp.patches()) {
    // Generate UV grid points for parameterization
//    gsMatrix<T> uv =
//        gsPointGrid(mp.parameterRange(0), patch->basis().size() * 2);
    gsMatrix<T> uv = mp.basis(0).anchors();

    // Evaluate the Mobius domain mapping
    gsMatrix<T> xieta;
    mobiusDomain.eval_into(uv, xieta);

    // Evaluate geometry
    gsMatrix<T> eval_geo(3, uv.cols());
    patch->eval_into(xieta, eval_geo);

    // Convert the patch's basis into a tensor B-spline basis
    gsTensorBSplineBasis<2, T> &bbasis =
        static_cast<gsTensorBSplineBasis<2, T> &>(patch->basis());

    // Fit surface using gsFitting and adjust parameters
    gsFitting<T> fittingSurface(uv, eval_geo, bbasis);
    fittingSurface.compute();
    fittingSurface.parameterCorrection();

    // Add the fitted patch to the result multipatch
    result.addPatch(*fittingSurface.result());
  }

  return result;
}

/**
 * @brief Class for surface reparameterization using Mobius mapping
 *
 * This class provides a method to reparameterize a surface using Mobius
 * mapping. It takes a multi-patch geometry and an optimizer as input, and
 * generates a reparameterized B-Spline surface.
 *
 * @tparam T The type of the input data (e.g., real_t)
 */
template <typename T> class SurfaceReparameterization
{
public:
  /**
   * @brief Constructor for SurfaceReparameterization
   *
   * @param patches The input multi-patch geometry
   * @param optimizer The optimizer to be used
   */
  // Constructor takes a multi-patch input and alpha matrix for the Mobius
  // domain
  explicit SurfaceReparameterization(const gsMultiPatch<T> &patches,
                                     gsOptimizer<T> &optimizer)
      : m_mp(patches), m_optimizer(optimizer)
  {
    m_mobiusDomain = gsMobiusMap<T>(gsMatrix<T, 2, 2>::Constant(0.5));
  }

  // Run the optimization process and generate the reparameterized B-Spline
  // surface
  gsMultiPatch<T> solve(T lambda1 = 1.0, T lambda2 = 1.0)
  {
    gsObjFuncSurface<2, T> objFuncSurface(m_mp, m_mobiusDomain);
    objFuncSurface.setLambda1(lambda1);
    objFuncSurface.setLambda2(lambda2);

    gsVector<T> initialGuessVector(4);
    initialGuessVector.setConstant(
        0.5); // Adjust this as necessary for better performance

    // Assign the optimizatio problem
    m_optimizer.setProblem(&objFuncSurface);

    // Perform the optimization
    m_optimizer.solve(initialGuessVector);

    return approximateWithBSpline(m_mp, m_optimizer.currentDesign());
  }

private:
  gsMultiPatch<T> m_mp;          // Input multi-patch geometry
  gsMobiusMap<T> m_mobiusDomain; // Mobius domain instance
  gsOptimizer<T> &m_optimizer;   // Optimizer
};

} // namespace gismo
