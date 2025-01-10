/** @file gsSurfaceReparameterization.h

    @brief Provides declaration for SurfaceReparameterization.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
*/

#pragma once

#include <gsNurbs/gsSquareDomain.h>
#include <gsNurbs/gsMobiusDomain.h>

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>

namespace gismo {
  template<short_t d=2, typename T=real_t>
  class gsObjFuncSurface : public gsOptProblem<T> {
   private:
	using geometryMap = typename gsExprAssembler<T>::geometryMap;
	using space = typename gsExprAssembler<T>::space;
	using solution = typename gsExprAssembler<T>::solution;

   public:
	explicit gsObjFuncSurface(const gsMultiPatch<T> &patches, const gsMobiusDomain<2,T> &mobiusDomain) : m_mp(patches), m_MobiusDomain(mobiusDomain) {
	  defaultOptions();

	  gsMatrix<T> bbox;
	  m_mp.boundingBox(bbox);
	  m_mp.patch(0).translate(-bbox.col(0));
	  m_mp.patch(0).scale(1/(bbox.col(1)-bbox.col(0)).array());

	  gsComposedGeometry<T> cgeom(m_MobiusDomain, m_mp.patch(0));
	  gsMultiBasis<T> dbasis(cgeom.basis());
	  m_evaluator.setIntegrationElements(dbasis);

	  // Set the geometry map
	  geometryMap G = m_evaluator.getMap(cgeom);
	  m_area = m_evaluator.integral(meas(G));
	}

	T evalObj(const gsAsConstVector<T> &coefsM) const final;

	void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const override;

	gsOptionList &options() { return m_options; }

	void defaultOptions();

	void addOptions(const gsOptionList &options);

	void applyOptions(const gsOptionList &options);

   protected:
	const gsMultiPatch<T> m_mp;
	const gsDofMapper m_mapper;
	const gsMultiBasis<T> m_mb;

	mutable gsExprEvaluator<T> m_evaluator;
	mutable gsExprAssembler<T> m_assembler;

	gsOptionList m_options;

	T m_lambda1 = 1.0, m_lambda2 = 1.0;
	T m_area = 1;
	gsComposedGeometry<T> m_cgeom;
	mutable gsMobiusDomain<2,T> m_MobiusDomain;
  };

  template<short_t d, typename T>
  T gsObjFuncSurface<d, T>::evalObj(const gsAsConstVector<T>& coefsM) const {
	m_MobiusDomain.updateGeom(coefsM);

	gsComposedGeometry<T> cgeom(m_MobiusDomain, m_mp.patch(0));
	gsMultiBasis<T> dbasis(cgeom.basis());
	m_evaluator.setIntegrationElements(dbasis);

	geometryMap G = m_evaluator.getMap(cgeom);
	auto FFF = jac(G).tr() * jac(G);
	auto m_integration = (FFF.trace() / meas(G)).val() + pow(FFF.det().val(), 2) / pow(m_area, 2);

	return m_evaluator.integral(m_integration);
  }

  template<short_t d, typename T>
  void gsObjFuncSurface<d, T>::gradObj_into(const gsAsConstVector<T>& u, gsAsVector<T>& result) const {
	const index_t n = u.rows();
	gsMatrix<T> uu = u;  // Create a copy
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

  template<short_t d, typename T>
  void gsObjFuncSurface<d, T>::defaultOptions() {
	m_options.addReal("qi_lambda1", "Sets the lambda 1 value", 1.0);
	m_options.addReal("qi_lambda2", "Sets the lambda 2 value", 1.0);
  }

  template<short_t d, typename T>
  void gsObjFuncSurface<d, T>::addOptions(const gsOptionList &options) {
	m_options.update(options, gsOptionList::addIfUnknown);
  }

  template<short_t d, typename T>
  void gsObjFuncSurface<d, T>::applyOptions(const gsOptionList &options) {
	m_options.update(options, gsOptionList::addIfUnknown);
	m_lambda1 = m_options.getReal("qi_lambda1");
	m_lambda2 = m_options.getReal("qi_lambda2");
	m_evaluator.options().update(m_options, gsOptionList::addIfUnknown);
  }

  template<class T = real_t>
  gsMultiPatch<T> convertIntoBSpline(const gsMultiPatch<T>& mp, const gsMatrix<T>& coefsMobiusIn) {
	gsMultiPatch<T> result;

	// Precompute Mobius domain matrix (use structured bindings in C++17)
	gsAsConstVector<T> coefsMobius(coefsMobiusIn.data(), 4);
	gsMatrix<T, 2, 2> alpha;
	alpha << coefsMobius(0), coefsMobius(2), coefsMobius(1), coefsMobius(3);
	gsMobiusDomain<2, T> mobiusDomain(alpha);

	for (const auto& patch : mp.patches()) {
	  // Generate UV grid points for parameterization
	  gsMatrix<T> uv = gsPointGrid(mp.parameterRange(0), patch->basis().size() * 2);

	  // Evaluate the Mobius domain mapping
	  gsMatrix<T> xieta;
	  mobiusDomain.eval_into(uv, xieta);

	  // Evaluate geometry
	  gsMatrix<T> eval_geo(3, uv.cols());
	  patch->eval_into(xieta, eval_geo);

	  // Convert the patch's basis into a tensor B-spline basis
	  auto& bbasis = static_cast<gsTensorBSplineBasis<2, T>&>(patch->basis());

	  // Fit surface using gsFitting and adjust parameters
	  gsFitting<T> fittingSurface(uv, eval_geo, bbasis);
	  fittingSurface.compute();
	  fittingSurface.parameterCorrection();

	  // Add the fitted patch to the result multipatch
	  result.addPatch(*fittingSurface.result());
	}

	return result;
  }


  template <typename T>
  class SurfaceReparameterization {
   public:
	// Constructor takes the multipatch input and alpha matrix for the Mobius domain
	SurfaceReparameterization(const gsMultiPatch<T>& patches)
		: m_mp(patches) {
	  gsMatrix<T, 2, 2> alpha;
	  alpha.setConstant(0.5);
	  m_mobiusDomain = gsMobiusDomain<2, T>(alpha);
	}

	// Run the optimization process and generate the reparameterized B-Spline surface
	gsMultiPatch<T> solve() {
	  gsObjFuncSurface<2, T> objFuncSurface(m_mp, m_mobiusDomain);

	  gsVector<T> initialGuessVector(4);
	  initialGuessVector.setConstant(0.5); // Adjust this as necessary for better performance

	  // Set up the optimizer
	  gsHLBFGS<T> optimizer(&objFuncSurface);
	  optimizer.options().setReal("MinGradientLength", 1e-6);
	  optimizer.options().setReal("MinStepLength", 1e-6);
	  optimizer.options().setInt("MaxIterations", 200);
	  optimizer.options().setInt("Verbose", 0);

	  // Perform the optimization
	  optimizer.solve(initialGuessVector);

	  return convertIntoBSpline(m_mp, optimizer.currentDesign());
	}

   private:
	gsMultiPatch<T> m_mp;                  // Input multi-patch geometry
	gsMobiusDomain<2, T> m_mobiusDomain;   // Mobius domain instance
  };
#endif
} // namespace gismo