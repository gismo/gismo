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

#include <gismo.h>
#include <gsHLBFGS/gsHLBFGS.h>
//#include <gsCore/gsFuncData.h>
//#include <gsAssembler/gsDirichletValues.h>
//#include <gsStructuralAnalysis/gsStaticNewton.h>
//#include <gsLBFGSpp/gsLBFGSpp.h>

using namespace gismo;

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
template<typename T>
void setOptimizerOptions(gsHLBFGS<T>& optimizer, const gsOptionList& options) {
  optimizer.options().setInt("MaxIterations", options.askInt("qi_MaxIterations", 1e4));
  optimizer.options().setReal("MinGradientLength", options.askReal("qi_MinGradientLength", 1e-4));
  optimizer.options().setReal("MinStepLength", options.askReal("qi_MinStepLength", 1e-4));
  optimizer.options().setInt("Verbose", options.askInt("Verbose", 0));
}

/// helper function to verbose log
void verboseLog(const std::string& message, const index_t& verbose) {
  if (verbose > 0) { gsInfo << message << "\n"; }
}

/// gsBarrierCore
template<short_t d, typename T= real_t>
struct gsBarrierCore {
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  static gsMultiPatch<T>
  compute(const gsMultiPatch<T> &mp, const gsDofMapper &mapper,
          const gsOptionList &options);

  static gsMultiPatch<T> computeBarrierPatch(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             const gsOptionList &options);

  static gsMultiPatch<T> computePenaltyPatch(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             const gsOptionList &options);

  static gsMultiPatch<T> computePenaltyPatch2(const gsMultiPatch<T> &mp,
                                              const gsDofMapper &mapper,
                                              const gsOptionList &options);

  // TODO: implement Jochen's PDE-based parameterization method using Newton's iteration
  static gsMultiPatch<T> computePDEPatch(const gsMultiPatch<T> &mp,
                                         const gsDofMapper &mapper,
                                         const gsOptionList &options);

  // Jochen's PDE-based parameterization method using Anderson Acceleration
  static gsMultiPatch<T> computePDEPatchAA(const gsMultiPatch<T> &mp,
                                           const gsDofMapper &mapper,
                                           const gsOptionList &options);

  static gsMultiPatch<T> computePDEPatchAAH1(const gsMultiPatch<T> &mp,
                                             const gsDofMapper &mapper,
                                             const gsOptionList &options);

  /// variational harmonic method
  static gsMultiPatch<T> computeVHPatch(const gsMultiPatch<T> &mp,
                                        const gsDofMapper &mapper,
                                        const gsOptionList &options);

  static gsMultiPatch<T> computeSmoothing(const gsMultiPatch<T> &mp,
                                          const gsDofMapper &mapper,
                                          const gsOptionList &options);

  /// Compute the area of the computational domain
  static T computeArea(const gsMultiPatch<T> &mp);

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

 public:

  /**
   * @brief      Translate and scale a multipatch
   *
   * @param      mp           { parameter_description }
   * @param[in]  translation  { parameter_description }
   * @param[in]  scaling      { parameter_description }
   */
//  static void
//  scaling(gsMultiPatch<T> &mp, const gsVector<T, d> &translation,
//          const gsVector<T, d> &scaling);

  /**
   * @brief      Undo translation and scaling of a multipatch
   *
   * @param      mp           { parameter_description }
   * @param[in]  translation  { parameter_description }
   * @param[in]  scaling      { parameter_description }
   */
//  static void
//  scalingUndo(gsMultiPatch<T> &mp, const gsVector<T, d> &translation,
//              const gsVector<T, d> &scaling);

  static gsOptionList defaultOptions();

 private:
  static void foldoverElimination(const gsMultiPatch<T>& mp,
                                  const gsDofMapper& mapper,
                                  gsVector<T>& initialGuessVector,
                                  const T& scaledArea,
                                  const gsOptionList& options);

  static void qualityImprovement(const gsMultiPatch<T>& mp,
                                 const gsDofMapper& mapper,
                                 gsVector<T>& initialGuessVector,
                                 const T& scaledArea,
                                 const gsOptionList& options);
};

template<short_t d, typename T = real_t>
class gsObjFoldoverFree : public gsOptProblem<T> {

 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjFoldoverFree(const gsMultiPatch<T> &patches, gsDofMapper  mapper);

  /// Evaluates the objective function at the given point u.
  T evalObj(const gsAsConstVector<T> &u) const override;

  /// Computes the gradient of the objective function at the given point u
  // and stores it in result.
  void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const override;

  void setDelta(const T& delta) {m_eps=delta;};

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

  T m_eps = T();  // need to handle later, set m_eps = 0.05*S, T() is the default value of T.
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
  void gradObj_into(const gsAsConstVector<T> &u, gsAsVector<T> &result) const override;

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

  void applyOptions();

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
class gsObjSmoothingPt : public gsOptProblem<T> {
 private:
  typedef typename gsExprAssembler<T>::geometryMap geometryMap;
  typedef typename gsExprAssembler<T>::space space;
  typedef typename gsExprAssembler<T>::solution solution;

 public:
  gsObjSmoothingPt(const gsMultiPatch<T> &patches,
                   gsDofMapper mapper);

  T evalObj(const gsAsConstVector<T> &u) const final;

  void gradObj2_into(const gsAsConstVector<T> &u,
                     gsAsVector<T> &result) const;

  gsOptionList &options() { return m_options; }

  void setEps(T tol) { m_eps = tol; }

  void defaultOptions();

  void addOptions(const gsOptionList &options);

  void applyOptions();

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

}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBarrierCore.hpp)
#endif

#ifdef USE_FLOAT_SCALAR
typedef float Scalar
#else
typedef double Scalar;
#endif

#ifdef EIGEN_DONT_ALIGN
#define EIGEN_ALIGNMENT Eigen::DontAlign
#else
#define EIGEN_ALIGNMENT gsEigen::AutoAlign
#endif

template<int Rows, int Cols, int Options = (gsEigen::ColMajor |
    EIGEN_ALIGNMENT)>
using MatrixT = gsEigen::Matrix<Scalar,
                              Rows,
                              Cols,
                              Options>;  ///< A typedef of the dense matrix of Eigen.
typedef MatrixT<gsEigen::Dynamic, 1>
    VectorX;                    ///< A nd column vector.
typedef MatrixT<gsEigen::Dynamic, gsEigen::Dynamic> MatrixXX;  ///< A n by m
// matrix.

template<typename T>
void
AALoop(const int &iter_, const gsVector<T> &current_F_, gsVector<T> &current_u_,
       gsMatrix<T> &prev_dG_, gsMatrix<T> &prev_dF_, gsMatrix<T> &M_,
       gsVector<T> &theta_, gsVector<T> &dF_scale_, const int &m_,
       const int &dim_, int &col_idx_, const real_t &resNorm_) {
  gsVector<T> G = current_F_ + current_u_;

  if (iter_ == 0) {
    prev_dF_.col(0) = -current_F_;
    prev_dG_.col(0) = -G;
    current_u_ = G;
  } else {
    prev_dF_.col(col_idx_) += current_F_;
    prev_dG_.col(col_idx_) += G;

    Scalar eps = 1e-14;
    Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
    dF_scale_(col_idx_) = scale;
    prev_dF_.col(col_idx_) /= scale;

    int m_k = std::min(m_, iter_);

    if (m_k == 1) {
      theta_(0) = 0;
      Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
      M_(0, 0) = dF_sqrnorm;
      Scalar dF_norm = std::sqrt(dF_sqrnorm);

      if (dF_norm > eps) {
        // compute theta = (dF * F) / (dF * dF)
        theta_(0) = (prev_dF_.col(col_idx_) / dF_norm).dot(
            current_F_ / dF_norm);
      }
    } else {
      // Update the normal equation matrix, for the column and row corresponding to the new dF column
      VectorX new_inner_prod = (prev_dF_.col(col_idx_).transpose()
          * prev_dF_.block(0, 0, dim_,
                           m_k)).transpose();
      M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
      M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

      // Solve normal equation
      gsEigen::CompleteOrthogonalDecomposition<MatrixXX> cod_;
      cod_.compute(M_.block(0, 0, m_k, m_k));
      theta_.head(m_k) = cod_.solve(
          prev_dF_.block(0, 0, dim_, m_k).transpose() *
              current_F_);
    }

    // Use rescaled theata to compute new u
    current_u_ =
        G
            - prev_dG_.block(0, 0, dim_, m_k)
                * ((theta_.head(m_k).array() /
                    dF_scale_.head(m_k).array())
                    .matrix());

    col_idx_ = (col_idx_ + 1) % m_;
    prev_dF_.col(col_idx_) = -current_F_;
    prev_dG_.col(col_idx_) = -G;
  }
}

template<typename T>
void
AALoopSR1(const int &iter_,
          gsVector<T> &current_F_,
          gsVector<T> &current_u_,
          gsMatrix<T> &prev_dG_,
          gsMatrix<T> &prev_dF_,
          gsMatrix<T> &M_,
          gsVector<T> &theta_,
          gsVector<T> &dF_scale_,
          const int &m_,
          const int &dim_,
          int &col_idx_,
          const real_t &resNorm_,
          gsMatrix<T> &H_) {
  gsVector<T> G = current_F_ + current_u_;

  if (iter_ == 0) {
    prev_dF_.col(0) = -current_F_;
    prev_dG_.col(0) = -G;
    current_u_ = G;
  } else {
    prev_dF_.col(col_idx_) += current_F_;
    prev_dG_.col(col_idx_) += G;

    // SR1
//        gsVector<> dx = prev_dG_.col(col_idx_)-prev_dF_.col(col_idx_);
//        gsVector<> commonVec = dx - H_ * prev_dF_.col(col_idx_);
//        H_ += (commonVec * commonVec.transpose()) / (prev_dF_.col(col_idx_).transpose() * commonVec);

//        gsVector<> dx = prev_dG_.col(col_idx_)-prev_dF_.col(col_idx_);
//        real_t commonFactor = dx.transpose() * prev_dF_.col(col_idx_);
////        gsDebugVar(1.0 + prev_dF_.col(col_idx_).transpose()*H_*prev_dF_.col(col_idx_)/commonFactor);
////        gsDebugVar((prev_dG_.col(col_idx_)* prev_dG_.col(col_idx_).transpose()) / commonFactor);
//        real_t num1 = (prev_dF_.col(col_idx_).transpose()*H_*prev_dF_.col(col_idx_)/commonFactor).value();
////        gsDebugVar(num1 * (prev_dG_.col(col_idx_)* prev_dG_.col(col_idx_).transpose()) / commonFactor);
//        H_ += (dx* dx.transpose()) / commonFactor +
//              num1 * (dx* dx.transpose()) / commonFactor -
//              (dx*prev_dF_.col(col_idx_).transpose()*H_+
//               H_*prev_dF_.col(col_idx_)*dx.transpose())/commonFactor;
//        current_F_ = -M.inverse() * current_F_;

    // SR1
    Scalar eps = 1e-14;
    Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
    dF_scale_(col_idx_) = scale;
    prev_dF_.col(col_idx_) /= scale;

    int m_k = std::min(m_, iter_);

    if (m_k == 1) {
      theta_(0) = 0;
      Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
      M_(0, 0) = dF_sqrnorm;
      Scalar dF_norm = std::sqrt(dF_sqrnorm);

      if (dF_norm > eps) {
        // compute theta = (dF * F) / (dF * dF)
        theta_(0) = (prev_dF_.col(col_idx_) / dF_norm).dot(
            current_F_ / dF_norm);
      }
    } else {
      // Update the normal equation matrix, for the column and row corresponding to the new dF column
      VectorX new_inner_prod = (prev_dF_.col(col_idx_).transpose()
          * prev_dF_.block(0, 0, dim_,
                           m_k)).transpose();
      M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
      M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

      // Solve normal equation
      gsEigen::CompleteOrthogonalDecomposition<MatrixXX> cod_;
      cod_.compute(M_.block(0, 0, m_k, m_k));
      theta_.head(m_k) = cod_.solve(
          prev_dF_.block(0, 0, dim_, m_k).transpose() *
              current_F_);
    }

    // Use rescaled theata to compute new u
    current_u_ =
        G
            - prev_dG_.block(0, 0, dim_, m_k)
                * ((theta_.head(m_k).array() /
                    dF_scale_.head(m_k).array())
                    .matrix());

    col_idx_ = (col_idx_ + 1) % m_;
    prev_dF_.col(col_idx_) = -current_F_;
    prev_dG_.col(col_idx_) = -G;
  }
}

template<typename T>
void AALoopAdpt(const int &iter_, const gsVector<T> &current_F_,
                gsVector<T> &current_u_,
                gsMatrix<T> &prev_dG_, gsMatrix<T> &prev_dF_, gsMatrix<T> &M_,
                gsVector<T> &theta_, gsVector<T> &dF_scale_, const int &m_,
                const int &dim_, int &col_idx_, const real_t &resNorm_) {
  gsVector<T> G = current_F_ + current_u_;

  if (iter_ == 0) {
    prev_dF_.col(0) = -current_F_;
    prev_dG_.col(0) = -G;
    current_u_ = G;
  } else {
    prev_dF_.col(col_idx_) += current_F_;
    prev_dG_.col(col_idx_) += G;

    Scalar eps = 1e-14;
    Scalar scale = std::max(eps, prev_dF_.col(col_idx_).norm());
    dF_scale_(col_idx_) = scale;
    prev_dF_.col(col_idx_) /= scale;

    int m_k = std::min(m_, iter_);

    if (m_k == 1) {
      theta_(0) = 0;
      Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
      M_(0, 0) = dF_sqrnorm;
      Scalar dF_norm = std::sqrt(dF_sqrnorm);

      if (dF_norm > eps) {
        // compute theta = (dF * F) / (dF * dF)
        theta_(0) = (prev_dF_.col(col_idx_) / dF_norm).dot(
            current_F_ / dF_norm);
      }
    } else {
      // Update the normal equation matrix, for the column and row corresponding to the new dF column
      VectorX new_inner_prod = (prev_dF_.col(col_idx_).transpose()
          * prev_dF_.block(0, 0, dim_,
                           m_k)).transpose();
      M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
      M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

      // Solve normal equation
      gsEigen::CompleteOrthogonalDecomposition<MatrixXX> cod_;
      cod_.compute(M_.block(0, 0, m_k, m_k));
      theta_.head(m_k) = cod_.solve(
          prev_dF_.block(0, 0, dim_, m_k).transpose() *
              current_F_);
    }

    T gamma = (current_F_ -
        prev_dF_.block(0, 0, dim_, m_k) * theta_.head(m_k)).norm() /
        current_F_.norm();
    T beta = 1.0 - gamma;

    // Use rescaled theata to compute new u
    current_u_ =
        G
            - beta * (prev_dG_.block(0, 0, dim_, m_k)
                * ((theta_.head(m_k).array() /
                    dF_scale_.head(m_k).array())
                    .matrix()));

    col_idx_ = (col_idx_ + 1) % m_;
    prev_dF_.col(col_idx_) = -current_F_;
    prev_dG_.col(col_idx_) = -G;
  }
}

template<typename T>
class AndersonAcceleration {
  typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;
  typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;
 public:
  AndersonAcceleration() = default;

  explicit AndersonAcceleration(int m)
      : m_(m),
        iter_(0),
        col_idx_(0) {
    GISMO_ASSERT(m > 0, "m should be greater than 0");
  }

  AndersonAcceleration(int m, int maxIter, T tol, int updateStep)
      : m_(m),
        maxIter_(maxIter),
        tolerance_(tol),
        updateStep_(updateStep),
        iter_(0),
        col_idx_(0) {
    GISMO_ASSERT(m > 0, "m should be greater than 0");
  }

// 残差向量和Jacobian句柄
  const gsVector<T> &compute(const gsVector<T> &u0,
                             const Residual_t &F,
                             std::vector<int> &iterHist,
                             std::vector<double> &resHist,
                             std::vector<double> &timeHist) {
    current_u_ = u0;
    init();

    double time = 0.0;
    timeHist.push_back(time);
    gsStopwatch timer;
    // Anderson iteration
    while (iter_ < 1000 && resNorm_ > 1e-5) {
      current_F_ = F(current_u_);
      resNorm_ = current_F_.norm();
      gsInfo << "iter: " << iter_ << ", resNorm: " << resNorm_ << "\n";

      AALoop(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
             theta_, dF_scale_, m_, dim_, col_idx_, resNorm_);

      iter_++;
      time += timer.stop();
      // track residual history and timing history
      iterHist.push_back(iter_ - 1);
      resHist.push_back(resNorm_);
      if (iter_ < 1000 && resNorm_ > 1e-5)
        timeHist.push_back(time);

      timer.restart();
    }
    return current_u_;
  }

  const gsVector<T> &computePrecond(const gsVector<T> &u0,
                                    const Residual_t &F,
                                    const Jacobian_t &Jacobian,
                                    std::vector<int> &iterHist,
                                    std::vector<double> &resHist,
                                    std::vector<double> &timeHist) {
    current_u_ = u0;
    init();

    double time = 0.0;
    timeHist.push_back(time);
    gsStopwatch timer;
    T orgresNorm = std::numeric_limits<real_t>::max();

    // Anderson iteration
    while (iter_ < maxIter_ && orgresNorm > tolerance_) {
      if (!(iter_ % updateStep_))
      {
        preconditioner_ = Jacobian(current_u_);
        linearSolver_.compute(preconditioner_);
      }
      gsMatrix<T> orgRes = F(current_u_);
      current_F_ = -linearSolver_.solve(orgRes);

      resNorm_ = current_F_.norm();
      gsInfo << "iter: " << iter_ << ", resNorm: " << resNorm_ << "\n";
      orgresNorm = orgRes.norm();

      AALoop(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
             theta_, dF_scale_, m_, dim_, col_idx_, resNorm_);

      iter_++;

      time += timer.stop();
      // track residual history and timing history
      iterHist.push_back(iter_ - 1);
      resHist.push_back(orgresNorm);
      if (iter_ < 1000 && orgresNorm > 1e-5)
        timeHist.push_back(time);

      timer.restart();
    }

    return current_u_;
  }

  const gsVector<T> &computeSR1(const gsVector<T> &u0, const Residual_t &F) {
    current_u_ = u0;
    init();

//        gsSparseMatrix<T> H(dim_, dim_);
    gsMatrix<T> H(dim_, dim_);
    H.setIdentity();

    // Anderson iteration
    while (iter_ < 1000 && resNorm_ > 1e-5) {
      current_F_ = F(current_u_);
      resNorm_ = current_F_.norm();
      gsInfo << "iter: " << iter_ << ", resNorm: " << resNorm_ << "\n";

//            current_F_ = -H*F(current_u_);
      AALoopSR1(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
                theta_, dF_scale_, m_, dim_, col_idx_, resNorm_, H);

      iter_++;
    }

    return current_u_;
  }

  // m: number of previous iterations used
  // d: dimension of variables
  // u0: initial variable values
  void init() {
    dim_ = current_u_.size();
    current_u_.resize(dim_);
    current_F_.resize(dim_);
    prev_dG_.resize(dim_, m_);
    prev_dF_.resize(dim_, m_);
    preconditioner_.resize(dim_, dim_);

    M_.resize(m_, m_);
    theta_.resize(m_);
    dF_scale_.resize(m_);
  }

 private:
//    Eigen::CompleteOrthogonalDecomposition<MatrixXX> linearSolver_;
  gsSparseSolver<>::LU linearSolver_;
  gsSparseMatrix<T> preconditioner_;
  gsVector<T> current_u_;
  gsVector<T> current_F_;
  gsMatrix<T> prev_dG_;
  gsMatrix<T> prev_dF_;
  gsMatrix<T> M_;        // Normal equations matrix for the computing theta
  gsVector<T> theta_;  // theta value computed from normal equations
  gsVector<T> dF_scale_;        // The scaling factor for each column of prev_dF

  T resNorm_ = std::numeric_limits<real_t>::max();
  int m_ =
      5;        // Number of previous iterates used for Andreson Acceleration
  int dim_ = -1;  // Dimension of variables
  int iter_ = 0;    // Iteration count since initialization
  int col_idx_ = 0;  // Index for history matrix column to store the next value
  int maxIter_ = 1e2; // maximum iteration
  T tolerance_ = 1e-5; // tolerance for convergence test
  int updateStep_ = 10; // update the preconditioner every updateStep step
};

template<typename T>
class AACompositeAA {
  typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;
  typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;

 public:
  AACompositeAA() = default;

  // m - outer loop window size, n - inner loop window size
  AACompositeAA(int m, int n)
      : m_(m),
        n_(n),
        iter_(0),
        col_idx_(0) {
    GISMO_ASSERT(m > 0, "m should be greater than 0");
    GISMO_ASSERT(n > 0, "n should be greater than 0");
  }

  // 残差向量和Jacobian句柄
  const gsVector<T> &compute(const gsVector<T> &u0,
                             const Residual_t &F,
                             std::vector<int> &iterHist,
                             std::vector<double> &resHist,
                             std::vector<double> &timeHist) {
    current_u_ = u0;
    init();

    double time = 0.0;
    timeHist.push_back(time);
    gsStopwatch timer;
    // Anderson iteration
    while (iter_ < 1000 && resNorm_ > 1e-5) {
      current_F_ = F(current_u_);
      resNorm_ = current_F_.norm();
//            gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_ << "\n";

      AALoop(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
             theta_,
             dF_scale_, m_, dim_, col_idx_, resNorm_);
      iter_++;

      time += timer.stop();
      // track residual history and timing history
      iterHist.push_back(iter_ - 1);
      resHist.push_back(resNorm_);
      if (iter_ < 1000 && resNorm_ > 1e-5)
        timeHist.push_back(time);

      timer.restart();

      if (iter_ >= m_) {
        int iterInner = 0;

        int col_idx = 0;
        gsMatrix<T> prev_dG;
        gsMatrix<T> prev_dF;
        gsMatrix<T> M;
        gsVector<T> theta;
        gsVector<T> dF_scale;

        prev_dG.resize(dim_, n_);
        prev_dF.resize(dim_, n_);

        M.resize(n_, n_);
        theta.resize(n_);
        dF_scale.resize(n_);

        while (iterInner <= n_ && resNorm_ > 1e-5) {
          current_F_ = F(current_u_);
          resNorm_ = current_F_.norm();
//                    gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_
//                           << "\n";

          AALoop(iterInner, current_F_, current_u_, prev_dG, prev_dF,
                 M, theta, dF_scale, n_, dim_, col_idx, resNorm_);
          iterInner++;
          iter_++;

          time += timer.stop();
          // track residual history and timing history
          iterHist.push_back(iter_ - 1);
          resHist.push_back(resNorm_);
          if (iter_ < 1000 && resNorm_ > 1e-5)
            timeHist.push_back(time);

          timer.restart();
        }
      }
    }

    return current_u_;
  }

  const gsVector<T> &computePrecond(const gsVector<T> &u0, const Residual_t &F,
                                    const Jacobian_t &Jacobian) {
    current_u_ = u0;
    init();

    // Anderson iteration
    while (iter_ < 1000 && resNorm_ > 1e-5) {
      if (!(iter_ % 10)) {
        preconditioner_ = Jacobian(current_u_);
        linearSolver_.compute(preconditioner_);
      }
      current_F_ = -linearSolver_.solve(F(current_u_));

      resNorm_ = current_F_.norm();
      gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_ << "\n";

      AALoop(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
             theta_,
             dF_scale_, m_, dim_, col_idx_, resNorm_);
      iter_++;

      if (iter_ >= m_) {
        int iterInner = 0;

        int col_idx = 0;
        gsMatrix<T> prev_dG;
        gsMatrix<T> prev_dF;
        gsMatrix<T> M;
        gsVector<T> theta;
        gsVector<T> dF_scale;

        prev_dG.resize(dim_, n_);
        prev_dF.resize(dim_, n_);

        M.resize(n_, n_);
        theta.resize(n_);
        dF_scale.resize(n_);

        while (iterInner <= n_ && resNorm_ > 1e-5) {
          if (!(iter_ % 10)) {
            preconditioner_ = Jacobian(current_u_);
            linearSolver_.compute(preconditioner_);
          }

          current_F_ = -linearSolver_.solve(F(current_u_));
          resNorm_ = current_F_.norm();
          gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_
                 << "\n";

          AALoop(iterInner, current_F_, current_u_, prev_dG, prev_dF,
                 M, theta,
                 dF_scale, n_, dim_, col_idx, resNorm_);
          iterInner++;
          iter_++;
        }
      }
    }
    return current_u_;
  }

  // m: number of previous iterations used
  // d: dimension of variables
  // u0: initial variable values
  void init() {
    dim_ = current_u_.size();
//        current_u_.resize(dim_);
    current_F_.resize(dim_);
    prev_dG_.resize(dim_, m_);
    prev_dF_.resize(dim_, m_);

    M_.resize(m_, m_);
    theta_.resize(m_);
    dF_scale_.resize(m_);
  }

 private:
  gsSparseSolver<>::LU linearSolver_;
  gsSparseMatrix<T> preconditioner_;
  gsVector<T> current_u_;
  gsVector<T> current_F_;
  gsMatrix<T> prev_dG_;
  gsMatrix<T> prev_dF_;
  gsMatrix<T> M_;        // Normal equations matrix for the computing theta
  gsVector<T> theta_;  // theta value computed from normal equations
  gsVector<T> dF_scale_;        // The scaling factor for each column of prev_dF

  real_t resNorm_ = std::numeric_limits<real_t>::max();
  int m_ = 5;        // Number of previous iterates used for AA
  int n_ = 5;     // Widow size of inner loop
  int dim_ = -1;  // Dimension of variables
  int iter_ = 0;    // Iteration count since initialization
  int col_idx_ = 0;  // Index for history matrix column to store the next value
};

template<typename T>
class AAAdaptiveCompositeAA {
  typedef std::function<gsVector<T>(gsVector<T> const &)> Residual_t;
  typedef std::function<gsSparseMatrix<T>(gsVector<T> const &)> Jacobian_t;

 public:
  AAAdaptiveCompositeAA() = default;

  // m - outer loop window size, n - inner loop window size
  AAAdaptiveCompositeAA(int m, int n)
      : m_(m),
        n_(n),
        iter_(0),
        col_idx_(0) {
    GISMO_ASSERT(m > 0, "m should be greater than 0");
    GISMO_ASSERT(n > 0, "n should be greater than 0");
  }

  // 残差向量和Jacobian句柄
  const gsVector<T> &compute(const gsVector<T> &u0,
                             const Residual_t &F,
                             std::vector<int> &iterHist,
                             std::vector<double> &resHist,
                             std::vector<double> &timeHist) {
    current_u_ = u0;
    init();

    double time = 0.0;
    timeHist.push_back(time);
    gsStopwatch timer;
    // Anderson iteration
    while (iter_ < 1000 && resNorm_ > 1e-5) {
      current_F_ = F(current_u_);
      resNorm_ = current_F_.norm();
//            gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_ << "\n";

      AALoopAdpt(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
                 theta_,
                 dF_scale_, m_, dim_, col_idx_, resNorm_);
      iter_++;

      time += timer.stop();
      // track residual history and timing history
      iterHist.push_back(iter_ - 1);
      resHist.push_back(resNorm_);
      if (iter_ < 1000 && resNorm_ > 1e-5)
        timeHist.push_back(time);

      timer.restart();

      if (iter_ >= m_) {
        int iterInner = 0;

        int col_idx = 0;
        gsMatrix<T> prev_dG;
        gsMatrix<T> prev_dF;
        gsMatrix<T> M;
        gsVector<T> theta;
        gsVector<T> dF_scale;

        prev_dG.resize(dim_, n_);
        prev_dF.resize(dim_, n_);

        M.resize(n_, n_);
        theta.resize(n_);
        dF_scale.resize(n_);

        while (iterInner <= n_ && resNorm_ > 1e-5) {
          current_F_ = F(current_u_);
          resNorm_ = current_F_.norm();
//                    gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_
//                           << "\n";

          AALoop(iterInner, current_F_, current_u_, prev_dG, prev_dF,
                 M, theta,
                 dF_scale, n_, dim_, col_idx, resNorm_);
          iterInner++;
          iter_++;

          time += timer.stop();
          // track residual history and timing history
          iterHist.push_back(iter_ - 1);
          resHist.push_back(resNorm_);
          if (iter_ < 1000 && resNorm_ > 1e-5)
            timeHist.push_back(time);

          timer.restart();
        }
      }
    }

    return current_u_;
  }

  const gsVector<T> &computePrecond(const gsVector<T> &u0, const Residual_t &F,
                                    const Jacobian_t &Jacobian) {
    current_u_ = u0;
    init();

    // Anderson iteration
    while (iter_ < 1000 && resNorm_ > 1e-5) {
      if (!(iter_ % 10)) {
        preconditioner_ = Jacobian(current_u_);
        linearSolver_.compute(preconditioner_);
      }
      current_F_ = -linearSolver_.solve(F(current_u_));

      resNorm_ = current_F_.norm();
      gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_ << "\n";

      AALoopAdpt(iter_, current_F_, current_u_, prev_dG_, prev_dF_, M_,
                 theta_, dF_scale_, m_, dim_, col_idx_, resNorm_);
      iter_++;

      if (iter_ >= m_) {
        int iterInner = 0;

        int col_idx = 0;
        gsMatrix<T> prev_dG;
        gsMatrix<T> prev_dF;
        gsMatrix<T> M;
        gsVector<T> theta;
        gsVector<T> dF_scale;

        prev_dG.resize(dim_, n_);
        prev_dF.resize(dim_, n_);

        M.resize(n_, n_);
        theta.resize(n_);
        dF_scale.resize(n_);

        while (iterInner <= n_ && resNorm_ > 1e-5) {
          if (!(iter_ % 10)) {
            preconditioner_ = Jacobian(current_u_);
            linearSolver_.compute(preconditioner_);
          }

          current_F_ = -linearSolver_.solve(F(current_u_));
          resNorm_ = current_F_.norm();
          gsInfo << "iter: " << iter_ << "  resNorm: " << resNorm_
                 << "\n";

          AALoop(iterInner, current_F_, current_u_, prev_dG, prev_dF,
                 M, theta,
                 dF_scale, n_, dim_, col_idx, resNorm_);
          iterInner++;
          iter_++;
        }
      }
    }

    return current_u_;
  }

  // m: number of previous iterations used
  // d: dimension of variables
  // u0: initial variable values
  void init() {
    dim_ = current_u_.size();
    current_u_.resize(dim_);
    current_F_.resize(dim_);
    prev_dG_.resize(dim_, m_);
    prev_dF_.resize(dim_, m_);

    M_.resize(m_, m_);
    theta_.resize(m_);
    dF_scale_.resize(m_);
  }

 private:
  gsSparseSolver<>::LU linearSolver_;
  gsSparseMatrix<T> preconditioner_;
  gsVector<T> current_u_;
  gsVector<T> current_F_;
  gsMatrix<T> prev_dG_;
  gsMatrix<T> prev_dF_;
  gsMatrix<T> M_;        // Normal equations matrix for the computing theta
  gsVector<T> theta_;  // theta value computed from normal equations
  gsVector<T> dF_scale_;        // The scaling factor for each column of prev_dF

  real_t resNorm_ = std::numeric_limits<real_t>::max();
  int m_ =
      5;        // Number of previous iterates used for Andreson Acceleration
  int n_ = 5;     // Widow size of inner loop
  int dim_ = -1;  // Dimension of variables
  int iter_ = 0;    // Iteration count since initialization
  int col_idx_ = 0;  // Index for history matrix column to store the next value
};
