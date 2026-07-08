/** @file gsAsG1GlobalBasis.hpp
    @brief Construct an exact, sparse transformation matrix for the Argyris
   (AS-G1) C1 space.
*/

#pragma once
#include "gsAsG1Basis.hpp"
#include "gsAsG1Domain.hpp"
#include <gismo.h>

namespace gismo {

namespace asg1global {

// ====================================================================
// gsVertexBasisMP -- analytic AS-G1 vertex shape function (one bfID).
//
// Ported from gsUnstructuredSplines (gsApproxC1Utils.h, gsVertexBasis,
// P. Weinmueller & A. Farahat).  Evaluates, in the *standard* frame
// (vertex at (0,0), incident edges along +u/+v), the closed-form
// Kapl-Sangalli-Takacs vertex basis function for a prescribed 6-jet
// column `bfID` of `Phi`, using exact linear alpha/beta from the
// (reparametrized) geometry.
// ====================================================================
template <class T> class gsVertexBasisMP : public gismo::gsFunction<T> {
protected:
  const gsGeometry<T> &m_geo;
  gsBasis<T> &m_basis;
  std::vector<gsBSpline<T>> m_alpha;
  std::vector<gsBSpline<T>> m_beta;
  std::vector<gsBSplineBasis<T>> m_basis_plus;
  std::vector<gsBSplineBasis<T>> m_basis_minus;
  const gsMatrix<T> m_Phi;
  const std::vector<bool> m_kindOfEdge;
  const index_t m_bfID;

public:
  gsVertexBasisMP(const gsGeometry<T> &geo, gsBasis<T> &basis,
                  std::vector<gsBSpline<T>> alpha,
                  std::vector<gsBSpline<T>> beta,
                  std::vector<gsBSplineBasis<T>> basis_plus,
                  std::vector<gsBSplineBasis<T>> basis_minus,
                  const gsMatrix<T> Phi, const std::vector<bool> kindOfEdge,
                  const index_t bfID)
      : m_geo(geo), m_basis(basis), m_alpha(alpha), m_beta(beta),
        m_basis_plus(basis_plus), m_basis_minus(basis_minus), m_Phi(Phi),
        m_kindOfEdge(kindOfEdge), m_bfID(bfID), _vertexBasis_piece(nullptr) {}

  ~gsVertexBasisMP() { delete _vertexBasis_piece; }
  GISMO_CLONE_FUNCTION(gsVertexBasisMP)
  short_t domainDim() const { return 2; }
  short_t targetDim() const { return 1; }

  mutable gsVertexBasisMP<T> *_vertexBasis_piece;
  const gsFunction<T> &piece(const index_t k) const {
    GISMO_UNUSED(k);
    _vertexBasis_piece = new gsVertexBasisMP(*this);
    return *_vertexBasis_piece;
  }

  void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const {
    result.resize(targetDim(), u.cols());
    result.setZero();

    gsMatrix<T> zero;
    zero.setZero(2, 1);

    std::vector<gsMatrix<T>> c_0, c_1;
    std::vector<gsMatrix<T>> c_0_plus, c_1_plus, c_2_plus;
    std::vector<gsMatrix<T>> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
    std::vector<gsMatrix<T>> c_0_minus, c_1_minus;
    for (index_t i = 0; i < 2; i++) {
      gsMatrix<T> b_0, b_1;
      gsMatrix<T> b_0_plus, b_1_plus, b_2_plus;
      gsMatrix<T> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
      gsMatrix<T> b_0_minus, b_1_minus;

      m_basis.component(i).evalSingle_into(0, u.row(i), b_0);
      m_basis.component(i).evalSingle_into(1, u.row(i), b_1);

      m_basis_plus[i].evalSingle_into(0, u.row(i), b_0_plus);
      m_basis_plus[i].evalSingle_into(1, u.row(i), b_1_plus);
      m_basis_plus[i].evalSingle_into(2, u.row(i), b_2_plus);

      m_basis_plus[i].derivSingle_into(0, u.row(i), b_0_plus_deriv);
      m_basis_plus[i].derivSingle_into(1, u.row(i), b_1_plus_deriv);
      m_basis_plus[i].derivSingle_into(2, u.row(i), b_2_plus_deriv);

      m_basis_minus[i].evalSingle_into(0, u.row(i), b_0_minus);
      m_basis_minus[i].evalSingle_into(1, u.row(i), b_1_minus);

      gsMatrix<T> b_1_0, b_1_minus_0;
      m_basis.component(i).derivSingle_into(1, zero.row(i), b_1_0);
      m_basis_minus[i].derivSingle_into(1, zero.row(i), b_1_minus_0);

      T factor_b_1 = 1.0 / b_1_0(0, 0);
      c_0.push_back(b_0 + b_1);
      c_1.push_back(factor_b_1 * b_1);

      T factor_b_1_minus = 1.0 / b_1_minus_0(0, 0);
      c_0_minus.push_back(b_0_minus + b_1_minus);
      c_1_minus.push_back(factor_b_1_minus * b_1_minus);

      gsMatrix<T> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
      m_basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
      m_basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
      m_basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

      T factor_c_1_plus = 1.0 / der_b_1_plus_0(0, 0);
      T factor2_c_1_plus = -der2_b_1_plus_0(0, 0) /
                           (der_b_1_plus_0(0, 0) * der2_b_2_plus_0(0, 0));
      T factor_c_2_plus = 1.0 / der2_b_2_plus_0(0, 0);

      c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
      c_1_plus.push_back(factor_c_1_plus * b_1_plus +
                         factor2_c_1_plus * b_2_plus);
      c_2_plus.push_back(factor_c_2_plus * b_2_plus);

      c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv +
                               b_2_plus_deriv);
      c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv +
                               factor2_c_1_plus * b_2_plus_deriv);
      c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);
    }

    std::vector<gsMatrix<T>> alpha, beta, alpha_0, beta_0, alpha_deriv,
        beta_deriv;
    gsMatrix<T> temp_mat;
    for (index_t i = 0; i < 2; ++i) {
      if (m_kindOfEdge[i]) {
        m_alpha[i].eval_into(u.row(i), temp_mat);
        alpha.push_back(temp_mat);
        m_alpha[i].eval_into(zero.row(0), temp_mat);
        alpha_0.push_back(temp_mat);
        m_alpha[i].deriv_into(zero.row(0), temp_mat);
        alpha_deriv.push_back(temp_mat);
        m_beta[i].eval_into(u.row(i), temp_mat);
        beta.push_back(temp_mat);
        m_beta[i].eval_into(zero.row(0), temp_mat);
        beta_0.push_back(temp_mat);
        m_beta[i].deriv_into(zero.row(0), temp_mat);
        beta_deriv.push_back(temp_mat);
      } else {
        temp_mat.setOnes(1, u.cols());
        alpha.push_back(temp_mat);
        temp_mat.setOnes(1, 1);
        alpha_0.push_back(temp_mat);
        temp_mat.setZero(1, 1);
        alpha_deriv.push_back(temp_mat);
        temp_mat.setZero(1, u.cols());
        beta.push_back(temp_mat);
        temp_mat.setZero(1, 1);
        beta_0.push_back(temp_mat);
        temp_mat.setZero(1, 1);
        beta_deriv.push_back(temp_mat);
      }
    }

    gsMatrix<T> geo_jac = m_geo.jacobian(zero);
    gsMatrix<T> geo_der2 = m_geo.deriv2(zero);

    gsMatrix<T> dd_ik_plus, dd_ik_minus, dd_ik_minus_deriv, dd_ik_plus_deriv;
    dd_ik_minus = -1.0 / (alpha_0[0](0, 0)) *
                  (geo_jac.col(1) + beta_0[0](0, 0) * geo_jac.col(0));
    dd_ik_plus = 1.0 / (alpha_0[1](0, 0)) *
                 (geo_jac.col(0) + beta_0[1](0, 0) * geo_jac.col(1));

    gsMatrix<T> geo_deriv2_12(2, 1), geo_deriv2_11(2, 1), geo_deriv2_22(2, 1);
    geo_deriv2_12.row(0) = geo_der2.row(2);
    geo_deriv2_12.row(1) = geo_der2.row(5);
    geo_deriv2_11.row(0) = geo_der2.row(0);
    geo_deriv2_11.row(1) = geo_der2.row(3);
    geo_deriv2_22.row(0) = geo_der2.row(1);
    geo_deriv2_22.row(1) = geo_der2.row(4);

    gsMatrix<T> alpha_squared_u = alpha_0[0] * alpha_0[0];
    gsMatrix<T> alpha_squared_v = alpha_0[1] * alpha_0[1];

    dd_ik_minus_deriv =
        -1.0 / (alpha_squared_u(0, 0)) *
        ((geo_deriv2_12 + (beta_deriv[0](0, 0) * geo_jac.col(0) +
                           beta_0[0](0, 0) * geo_deriv2_11)) *
             alpha_0[0](0, 0) -
         (geo_jac.col(1) + beta_0[0](0, 0) * geo_jac.col(0)) *
             alpha_deriv[0](0, 0));
    dd_ik_plus_deriv = 1.0 / (alpha_squared_v(0, 0)) *
                       ((geo_deriv2_12 + (beta_deriv[1](0, 0) * geo_jac.col(1) +
                                          beta_0[1](0, 0) * geo_deriv2_22)) *
                            alpha_0[1](0, 0) -
                        (geo_jac.col(0) + beta_0[1](0, 0) * geo_jac.col(1)) *
                            alpha_deriv[1](0, 0));

    std::vector<gsMatrix<T>> d_ik;
    d_ik.push_back(m_Phi.col(0));
    d_ik.push_back(m_Phi.block(1, 0, 2, 6).transpose() * geo_jac.col(0));
    d_ik.push_back(m_Phi.block(1, 0, 2, 6).transpose() * geo_jac.col(1));
    d_ik.push_back(
        (geo_jac(0, 0) * m_Phi.col(3) + geo_jac(1, 0) * m_Phi.col(4)) *
            geo_jac(0, 1) +
        (geo_jac(0, 0) * m_Phi.col(4) + geo_jac(1, 0) * m_Phi.col(5)) *
            geo_jac(1, 1) +
        m_Phi.block(0, 1, 6, 1) * geo_der2.row(2) +
        m_Phi.block(0, 2, 6, 1) * geo_der2.row(5));

    std::vector<gsMatrix<T>> d_ilik_minus, d_ilik_plus;
    d_ilik_minus.push_back(m_Phi.col(0));
    d_ilik_minus.push_back(m_Phi.block(1, 0, 2, 6).transpose() *
                           geo_jac.col(0));
    d_ilik_minus.push_back(
        (geo_jac(0, 0) * m_Phi.col(3) + geo_jac(1, 0) * m_Phi.col(4)) *
            geo_jac(0, 0) +
        (geo_jac(0, 0) * m_Phi.col(4) + geo_jac(1, 0) * m_Phi.col(5)) *
            geo_jac(1, 0) +
        m_Phi.block(0, 1, 6, 1) * geo_der2.row(0) +
        m_Phi.block(0, 2, 6, 1) * geo_der2.row(3));
    d_ilik_minus.push_back(m_Phi.block(1, 0, 2, 6).transpose() * dd_ik_minus);
    d_ilik_minus.push_back(
        (geo_jac(0, 0) * m_Phi.col(3) + geo_jac(1, 0) * m_Phi.col(4)) *
            dd_ik_minus(0, 0) +
        (geo_jac(0, 0) * m_Phi.col(4) + geo_jac(1, 0) * m_Phi.col(5)) *
            dd_ik_minus(1, 0) +
        m_Phi.block(0, 1, 6, 1) * dd_ik_minus_deriv.row(0) +
        m_Phi.block(0, 2, 6, 1) * dd_ik_minus_deriv.row(1));

    d_ilik_plus.push_back(m_Phi.col(0));
    d_ilik_plus.push_back(m_Phi.block(1, 0, 2, 6).transpose() * geo_jac.col(1));
    d_ilik_plus.push_back(
        (geo_jac(0, 1) * m_Phi.col(3) + geo_jac(1, 1) * m_Phi.col(4)) *
            geo_jac(0, 1) +
        (geo_jac(0, 1) * m_Phi.col(4) + geo_jac(1, 1) * m_Phi.col(5)) *
            geo_jac(1, 1) +
        m_Phi.block(0, 1, 6, 1) * geo_der2.row(1) +
        m_Phi.block(0, 2, 6, 1) * geo_der2.row(4));
    d_ilik_plus.push_back(m_Phi.block(1, 0, 2, 6).transpose() * dd_ik_plus);
    d_ilik_plus.push_back(
        (geo_jac(0, 1) * m_Phi.col(3) + geo_jac(1, 1) * m_Phi.col(4)) *
            dd_ik_plus(0, 0) +
        (geo_jac(0, 1) * m_Phi.col(4) + geo_jac(1, 1) * m_Phi.col(5)) *
            dd_ik_plus(1, 0) +
        m_Phi.block(0, 1, 6, 1) * dd_ik_plus_deriv.row(0) +
        m_Phi.block(0, 2, 6, 1) * dd_ik_plus_deriv.row(1));

    result =
        d_ilik_minus.at(0)(m_bfID, 0) *
            (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
             beta[0].cwiseProduct(
                 c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
        d_ilik_minus.at(1)(m_bfID, 0) *
            (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
             beta[0].cwiseProduct(
                 c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
        d_ilik_minus.at(2)(m_bfID, 0) *
            (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
             beta[0].cwiseProduct(
                 c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
        d_ilik_minus.at(3)(m_bfID, 0) *
            alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
        d_ilik_minus.at(4)(m_bfID, 0) *
            alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1)));

    result +=
        d_ilik_plus.at(0)(m_bfID, 0) *
            (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
             beta[1].cwiseProduct(
                 c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
        d_ilik_plus.at(1)(m_bfID, 0) *
            (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
             beta[1].cwiseProduct(
                 c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
        d_ilik_plus.at(2)(m_bfID, 0) *
            (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
             beta[1].cwiseProduct(
                 c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
        d_ilik_plus.at(3)(m_bfID, 0) *
            alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
        d_ilik_plus.at(4)(m_bfID, 0) *
            alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0)));

    result -= d_ik.at(0)(m_bfID, 0) * c_0.at(0).cwiseProduct(c_0.at(1)) +
              d_ik.at(2)(m_bfID, 0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
              d_ik.at(1)(m_bfID, 0) * c_1.at(0).cwiseProduct(c_0.at(1)) +
              d_ik.at(3)(m_bfID, 0) * c_1.at(0).cwiseProduct(c_1.at(1));
  }
};

/// Reparametrize `geo` so the given corner `cIdx` (1..4) maps to the
/// standard frame (vertex at (0,0), incident edges along +u, +v),
/// mirroring gsUnstructuredSplines' gsG1AuxiliaryPatch.  Returns the
/// rotated patch `geoRot` and a permutation `perm` with
/// `perm[rotDof] = origDof`.
template <class T>
void reparamCorner(const gsTensorBSpline<2, T> &geo, index_t cIdx,
                   gsTensorBSpline<2, T> &geoRot, std::vector<index_t> &perm) {
  gsKnotVector<T> kvU = geo.basis().knots(0);
  gsKnotVector<T> kvV = geo.basis().knots(1);
  index_t n0 = geo.basis().size(0), n1 = geo.basis().size(1);
  const index_t N = n0 * n1;
  gsMatrix<T> M(N, 3);
  M.leftCols(2) = geo.coefs();
  for (index_t r = 0; r < N; ++r)
    M(r, 2) = T(r);

  auto antiClock = [&]() {
    gsMatrix<T> Mn(N, 3);
    for (index_t j = 0; j < n0; ++j)
      for (index_t i = 0; i < n1; ++i)
        Mn.row(i + j * n1) = M.row((n0 - 1 - j) + n0 * i);
    M = Mn;
    std::swap(kvU, kvV);
    std::swap(n0, n1);
  };
  auto clock = [&]() {
    gsMatrix<T> Mn(N, 3);
    for (index_t j = n0 - 1; j >= 0; --j)
      for (index_t i = 0; i < n1; ++i)
        Mn.row(i + (n0 - 1 - j) * n1) = M.row((n1 * n0 - 1 - j) - n0 * i);
    M = Mn;
    std::swap(kvU, kvV);
    std::swap(n0, n1);
  };
  auto twice = [&]() {
    gsMatrix<T> Mn(N, 3);
    for (index_t r = 0; r < N; ++r)
      Mn.row(r) = M.row(N - 1 - r);
    M = Mn;
  };
  auto swapAx = [&]() {
    gsMatrix<T> Mn(N, 3);
    for (index_t j = 0; j < n0; ++j)
      for (index_t i = 0; i < n1; ++i)
        Mn.row(i + j * n1) = M.row(j + n0 * i);
    M = Mn;
    std::swap(kvU, kvV);
    std::swap(n0, n1);
  };

  const bool switched = (geo.orientation() == -1);
  if (switched)
    swapAx();
  if (!switched) {
    if (cIdx == 2)
      antiClock();
    else if (cIdx == 3)
      clock();
    else if (cIdx == 4)
      twice();
  } else {
    if (cIdx == 2)
      clock();
    else if (cIdx == 3)
      antiClock();
    else if (cIdx == 4)
      twice();
  }

  gsMatrix<T> coefsRot = M.leftCols(2);
  geoRot = gsTensorBSpline<2, T>(kvU, kvV, coefsRot);
  perm.resize(N);
  for (index_t r = 0; r < N; ++r)
    perm[r] = static_cast<index_t>(std::lround(M(r, 2)));
}

/// 2-D corner index (1..4) at the low/high tangent end of a side.
/// `tdEnd` is 0 for the low tangent end, 1 for the high end.
inline index_t cornerAtSideEnd(const patchSide &s, int tdEnd) {
  const short_t d = s.direction();
  const bool pr = s.parameter();
  index_t upar, vpar;
  if (d == 0) {
    upar = pr ? 1 : 0;
    vpar = tdEnd;
  } else {
    vpar = pr ? 1 : 0;
    upar = tdEnd;
  }
  return 1 + upar + 2 * vpar;
}

} // namespace asg1global

template <typename T> class gsAsG1GlobalBasis {
public:
  gsAsG1GlobalBasis(const gsMultiPatch<T> &mp, const gsMultiBasis<T> &mb)
      : m_mp(mp), m_mb(mb) {
    // Compute edge gluing data natively
    m_M_gd = computeGluingData(m_mp);
    assemble();
  }

  /// Return the transformation matrix: Global_C1_DOFs -> Broken_Patch_DOFs
  gsSparseMatrix<T> getTransformationMatrix() const { return m_M; }

private:
  const gsMultiPatch<T> &m_mp;
  const gsMultiBasis<T> &m_mb;
  gsMatrix<T> m_M_gd;
  gsSparseMatrix<T> m_M;

  struct EdgeInfo {
    int p1;
    boxSide s1;
    int p2;
    boxSide s2; // -1 if boundary
    bool flipped;
  };

  struct Corner {
    int p;
    int c; // 0:(0,0), 1:(1,0), 2:(0,1), 3:(1,1)
  };

  bool getGluing(int p, boxSide s, T &a0, T &a1, T &b0, T &b1,
                 T &tangentSign) const {
    a0 = m_M_gd(p, 4 * (s.index() - 1) + 0);
    a1 = m_M_gd(p, 4 * (s.index() - 1) + 1);
    b0 = m_M_gd(p, 4 * (s.index() - 1) + 2);
    b1 = m_M_gd(p, 4 * (s.index() - 1) + 3);
    tangentSign = 1.0;

    patchSide ps(p, s);
    for (size_t i = 0; i < m_mp.interfaces().size(); ++i) {
      const boundaryInterface &interf = m_mp.interfaces()[i];
      if (interf.second() == ps || interf.first() == ps) {
        // tangentSign must be identical for BOTH sides of the
        // interface: it is a property of the (flipped) edge, not of
        // an individual patch side. It multiplies the beta*dVm term
        // in the trace/smoother (normal-derivative) gamma solve.
        short_t tDir1 = 1 - interf.first().direction();
        bool tangentialFlipped = !interf.dirOrientation(interf.first(), tDir1);
        if (tangentialFlipped)
          tangentSign = -1.0;
        return true;
      }
    }
    // Fallback for boundaries to avoid division by zero (identity mapping)
    a0 = 1.0;
    a1 = 1.0;
    b0 = 0.0;
    b1 = 0.0;
    return false;
  }

  void assemble() {
    int nPatches = m_mp.nPatches();
    std::vector<int> offsets(nPatches + 1, 0);
    for (int p = 0; p < nPatches; ++p)
      offsets[p + 1] = offsets[p] + m_mb.basis(p).size();
    int total_broken = offsets[nPatches];

    // 1. Cache all per-side embeddings
    std::vector<std::vector<gsSparseMatrix<T>>> side_matrices(
        nPatches, std::vector<gsSparseMatrix<T>>(4));
    for (int p = 0; p < nPatches; ++p) {
      for (int s = 1; s <= 4; ++s) {
        boxSide bs(s);
        T a0, a1, b0, b1, tSign;
        getGluing(p, bs, a0, a1, b0, b1, tSign);
        const gsTensorBSplineBasis<2, T> &tb =
            static_cast<const gsTensorBSplineBasis<2, T> &>(m_mb.basis(p));
        side_matrices[p][s - 1] =
            createGluingDataArgyrisBasis(tb, bs, a0, a1, b0, b1, 1e-12, tSign);
      }
    }

    int current_col = 0;
    std::vector<gsEigen::Triplet<T>> triplets;

    // 2. Patch Interiors
    for (int p = 0; p < nPatches; ++p) {
      const gsTensorBSplineBasis<2, T> &tb =
          static_cast<const gsTensorBSplineBasis<2, T> &>(m_mb.basis(p));
      int N1 = tb.size(0);
      int N2 = tb.size(1);
      for (int j = 2; j <= N2 - 3; ++j) {
        for (int i = 2; i <= N1 - 3; ++i) {
          triplets.push_back(
              gsEigen::Triplet<T>(offsets[p] + i + j * N1, current_col++, 1.0));
        }
      }
    }

    // 3. Edge Interiors
    std::vector<EdgeInfo> edges;
    for (size_t i = 0; i < m_mp.interfaces().size(); ++i) {
      const boundaryInterface &interf = m_mp.interfaces()[i];
      edges.push_back({interf.first().patch, interf.first().side(),
                       interf.second().patch, interf.second().side(),
                       !interf.dirOrientation(interf.first(),
                                              1 - interf.first().direction())});
    }
    for (size_t i = 0; i < m_mp.boundaries().size(); ++i) {
      patchSide ps = m_mp.boundaries()[i];
      edges.push_back({ps.patch, ps.side(), -1, boxSide(0), false});
    }

    for (const auto &e : edges) {
      gsSparseMatrix<T> &M1 = side_matrices[e.p1][e.s1.index() - 1];
      const gsTensorBSplineBasis<2, T> &tb1 =
          static_cast<const gsTensorBSplineBasis<2, T> &>(m_mb.basis(e.p1));

      int nInt1 = tb1.size() - tb1.boundary(e.s1).rows() -
                  tb1.boundaryOffset(e.s1, 1).rows();

      gsSparseMatrix<T> *M2_ptr = nullptr;
      int nInt2 = 0;
      if (e.p2 != -1) {
        M2_ptr = &side_matrices[e.p2][e.s2.index() - 1];
        const gsTensorBSplineBasis<2, T> &tb2 =
            static_cast<const gsTensorBSplineBasis<2, T> &>(m_mb.basis(e.p2));
        nInt2 = tb2.size() - tb2.boundary(e.s2).rows() -
                tb2.boundaryOffset(e.s2, 1).rows();
      }

      // In `createGluingDataArgyrisBasis`, columns are ordered: [Interior |
      // LowerDeg | Smoother] We strip out the overlapping endpoints interacting
      // with Vertex elements.
      gsBSplineBasis<T> bS1 = *tb1.boundaryBasis(e.s1);
      bS1.elevateContinuity(1);
      gsBSplineBasis<T> bL1 = *tb1.boundaryBasis(e.s1);
      bL1.degreeReduce(1);
      int nL = bL1.size(), nS = bS1.size();

      for (int k = 2; k <= nL - 3; ++k) {
        for (typename gsSparseMatrix<T>::InnerIterator it(M1, nInt1 + k); it;
             ++it)
          triplets.push_back(gsEigen::Triplet<T>(offsets[e.p1] + it.row(),
                                                 current_col, it.value()));
        if (M2_ptr) {
          int k2 = e.flipped ? (nL - 1 - k) : k;
          // The lower-degree (normal-derivative) block flips sign on
          // patch 2 for non-flipped interfaces, since the two patch
          // normals point in opposite directions across the edge.
          const T l2Sign = e.flipped ? T(1) : T(-1);
          for (typename gsSparseMatrix<T>::InnerIterator it(*M2_ptr,
                                                            nInt2 + k2);
               it; ++it)
            triplets.push_back(gsEigen::Triplet<T>(
                offsets[e.p2] + it.row(), current_col, l2Sign * it.value()));
        }
        current_col++;
      }
      for (int k = 3; k <= nS - 4; ++k) {
        for (typename gsSparseMatrix<T>::InnerIterator it(M1, nInt1 + nL + k);
             it; ++it)
          triplets.push_back(gsEigen::Triplet<T>(offsets[e.p1] + it.row(),
                                                 current_col, it.value()));
        if (M2_ptr) {
          int k2 = e.flipped ? (nS - 1 - k) : k;
          for (typename gsSparseMatrix<T>::InnerIterator it(*M2_ptr,
                                                            nInt2 + nL + k2);
               it; ++it)
            triplets.push_back(gsEigen::Triplet<T>(offsets[e.p2] + it.row(),
                                                   current_col, it.value()));
        }
        current_col++;
      }
    }

    // 4. Vertex functions (Kapl-Sangalli-Takacs closed form, ported
    //    from gsUnstructuredSplines). One block of 6 functions is added
    //    at every vertex where at least two interfaces meet; these are
    //    exactly the near-vertex DOFs that the edge blocks dropped.
    {
      using asg1global::cornerAtSideEnd;
      using asg1global::gsVertexBasisMP;
      using asg1global::reparamCorner;

      // Enumerate unique vertices by coincident corner coordinates.
      std::vector<Corner> all_corners;
      std::vector<gsMatrix<T>> coords;
      for (int p = 0; p < nPatches; ++p) {
        for (int c = 0; c < 4; ++c) {
          gsMatrix<T> uv(2, 1);
          uv(0, 0) = (c == 1 || c == 3) ? 1.0 : 0.0;
          uv(1, 0) = (c == 2 || c == 3) ? 1.0 : 0.0;
          all_corners.push_back({p, c});
          coords.push_back(m_mp.patch(p).eval(uv));
        }
      }

      std::vector<std::vector<Corner>> unique_vertices;
      std::vector<bool> visited(all_corners.size(), false);
      for (size_t i = 0; i < all_corners.size(); ++i) {
        if (visited[i])
          continue;
        std::vector<Corner> group;
        group.push_back(all_corners[i]);
        visited[i] = true;
        for (size_t j = i + 1; j < all_corners.size(); ++j) {
          if (!visited[j] && (coords[i] - coords[j]).norm() < 1e-8) {
            group.push_back(all_corners[j]);
            visited[j] = true;
          }
        }
        unique_vertices.push_back(group);
      }

      for (const auto &v_group : unique_vertices) {
        // Count interfaces meeting at this vertex.
        auto cornerInGroup = [&](int patch, index_t cIdx1) -> bool {
          for (const Corner &cc : v_group)
            if (cc.p == patch && cc.c == cIdx1 - 1)
              return true;
          return false;
        };
        index_t ifCount = 0;
        for (size_t i = 0; i < m_mp.interfaces().size(); ++i) {
          const boundaryInterface &interf = m_mp.interfaces()[i];
          const patchSide ps1 = interf.first(), ps2 = interf.second();
          const bool on1 = cornerInGroup(ps1.patch, cornerAtSideEnd(ps1, 0)) ||
                           cornerInGroup(ps1.patch, cornerAtSideEnd(ps1, 1));
          const bool on2 = cornerInGroup(ps2.patch, cornerAtSideEnd(ps2, 0)) ||
                           cornerInGroup(ps2.patch, cornerAtSideEnd(ps2, 1));
          if (on1 && on2)
            ++ifCount;
        }
        if (ifCount < 2)
          continue; // regular boundary vertex: no block

        const index_t nInc = static_cast<index_t>(v_group.size());

        // Reparametrize each incident patch to the standard frame.
        std::vector<gsTensorBSpline<2, T>> geoRot(nInc);
        std::vector<std::vector<index_t>> perm(nInc);
        std::vector<index_t> patchOf(nInc);
        for (index_t k = 0; k < nInc; ++k) {
          patchOf[k] = v_group[k].p;
          const gsTensorBSpline<2, T> &geo =
              dynamic_cast<const gsTensorBSpline<2, T> &>(
                  m_mp.patch(patchOf[k]));
          reparamCorner<T>(geo, v_group[k].c + 1, geoRot[k], perm[k]);
        }

        // Uniform scaling sigma (same heuristic as gsUnstructuredSplines).
        T pdeg = 0, h_geo = 1;
        for (index_t k = 0; k < nInc; ++k) {
          const gsTensorBSplineBasis<2, T> &b = geoRot[k].basis();
          pdeg = std::max<T>(pdeg, std::max(b.degree(0), b.degree(1)));
        }
        for (index_t k = 0; k < nInc; ++k) {
          const gsTensorBSplineBasis<2, T> &b = geoRot[k].basis();
          for (index_t j = 0; j < 2; ++j)
            h_geo = std::min<T>(h_geo,
                                b.knots(j).at(static_cast<size_t>(pdeg) + 1));
        }
        T sigmaRaw = 0;
        {
          gsMatrix<T> zero;
          zero.setZero(2, 1);
          for (index_t k = 0; k < nInc; ++k) {
            gsMatrix<T> d = geoRot[k].deriv(zero);
            sigmaRaw += d.template lpNorm<gsEigen::Infinity>();
          }
        }
        const T sigmaScale = sigmaRaw * h_geo / (T(nInc) * pdeg);
        const T sigma = (sigmaScale != T(0)) ? T(1) / sigmaScale : T(1);

        gsMatrix<T> Phi(6, 6);
        Phi.setIdentity();
        Phi.col(1) *= sigma;
        Phi.col(2) *= sigma;
        Phi.col(3) *= sigma * sigma;
        Phi.col(4) *= sigma * sigma;
        Phi.col(5) *= sigma * sigma;

        const int v_col = current_col;
        current_col += 6;

        for (index_t k = 0; k < nInc; ++k) {
          const index_t pk = patchOf[k];
          const index_t corner = v_group[k].c;
          gsTensorBSplineBasis<2, T> &tbRot =
              const_cast<gsTensorBSplineBasis<2, T> &>(geoRot[k].basis());

          boxSide w_side =
              (corner == 0 || corner == 3) ? boxSide(1) : boxSide(2);
          gsSparseMatrix<T> &M_West_edge =
              side_matrices[pk][w_side.index() - 1];
          boxSide s_side =
              (corner == 0 || corner == 1) ? boxSide(3) : boxSide(4);
          gsSparseMatrix<T> &M_South_edge =
              side_matrices[pk][s_side.index() - 1];

          // Find side basis functions for West/East edge
          gsBSplineBasis<T> bS_W = *tbRot.boundaryBasis(w_side);
          bS_W.elevateContinuity(1);
          gsBSplineBasis<T> bL_W = *tbRot.boundaryBasis(w_side);
          bL_W.degreeReduce(1);
          int nS_W = bS_W.size(), nL_W = bL_W.size();
          int nInt_W = M_West_edge.cols() - nL_W - nS_W;
          bool w_start = (corner == 0 || corner == 1);
          int w_idxS = w_start ? 0 : nS_W - 3;
          int w_idxL = w_start ? 0 : nL_W - 2;
          std::vector<int> cols_W = {nInt_W + nL_W + w_idxS + 0,
                                     nInt_W + nL_W + w_idxS + 1,
                                     nInt_W + nL_W + w_idxS + 2,
                                     nInt_W + w_idxL + 0, nInt_W + w_idxL + 1};

          // Find side basis functions for South/North edge
          gsBSplineBasis<T> bS_S = *tbRot.boundaryBasis(s_side);
          bS_S.elevateContinuity(1);
          gsBSplineBasis<T> bL_S = *tbRot.boundaryBasis(s_side);
          bL_S.degreeReduce(1);
          int nS_S = bS_S.size(), nL_S = bL_S.size();
          int nInt_S = M_South_edge.cols() - nL_S - nS_S;
          bool s_start = (corner == 0 || corner == 3);
          int s_idxS = s_start ? 0 : nS_S - 3;
          int s_idxL = s_start ? 0 : nL_S - 2;
          std::vector<int> cols_S = {nInt_S + nL_S + s_idxS + 0,
                                     nInt_S + nL_S + s_idxS + 1,
                                     nInt_S + nL_S + s_idxS + 2,
                                     nInt_S + s_idxL + 0, nInt_S + s_idxL + 1};

          // Evaluate Jacobian and Hessian of the reparametrized geometry at
          // (0,0)
          gsMatrix<T> zero_uv = gsMatrix<T>::Zero(2, 1);
          gsMatrix<T> J = geoRot[k].jacobian(zero_uv);
          gsMatrix<T> secDers;
          geoRot[k].deriv2_into(zero_uv, secDers);

          T J00 = J(0, 0), J10 = J(1, 0);
          T J01 = J(0, 1), J11 = J(1, 1);
          T X_uu = secDers(0, 0), X_vv = secDers(1, 0), X_uv = secDers(2, 0);
          T Y_uu = secDers(3, 0), Y_vv = secDers(4, 0), Y_uv = secDers(5, 0);

          // Prepare univariate evaluations of basis functions and their
          // first/second derivatives at 0.
          const gsBSplineBasis<T> &b0 =
              dynamic_cast<const gsBSplineBasis<T> &>(tbRot.component(0));
          const gsBSplineBasis<T> &b1 =
              dynamic_cast<const gsBSplineBasis<T> &>(tbRot.component(1));
          gsMatrix<T> pt = gsMatrix<T>::Zero(1, 1);

          gsMatrix<T> val0_mat, der0_mat, der20_mat;
          b0.eval_into(pt, val0_mat);
          b0.deriv_into(pt, der0_mat);
          b0.deriv2_into(pt, der20_mat);
          gsMatrix<T> val1_mat, der1_mat, der21_mat;
          b1.eval_into(pt, val1_mat);
          b1.deriv_into(pt, der1_mat);
          b1.deriv2_into(pt, der21_mat);

          gsVector<T> eval_val = gsVector<T>::Zero(tbRot.size());
          gsVector<T> eval_u = gsVector<T>::Zero(tbRot.size());
          gsVector<T> eval_v = gsVector<T>::Zero(tbRot.size());
          gsVector<T> eval_uv = gsVector<T>::Zero(tbRot.size());
          gsVector<T> eval_uu = gsVector<T>::Zero(tbRot.size());
          gsVector<T> eval_vv = gsVector<T>::Zero(tbRot.size());

          gsMatrix<index_t> act0 = b0.active(pt);
          gsMatrix<index_t> act1 = b1.active(pt);
          for (int a = 0; a < act0.rows(); ++a) {
            for (int b = 0; b < act1.rows(); ++b) {
              int r = act0(a, 0) + act1(b, 0) * tbRot.size(0);
              eval_val(r) = val0_mat(a, 0) * val1_mat(b, 0);
              eval_u(r) = der0_mat(a, 0) * val1_mat(b, 0);
              eval_v(r) = val0_mat(a, 0) * der1_mat(b, 0);
              eval_uv(r) = der0_mat(a, 0) * der1_mat(b, 0);
              eval_uu(r) = der20_mat(a, 0) * val1_mat(b, 0);
              eval_vv(r) = val0_mat(a, 0) * der21_mat(b, 0);
            }
          }

          gsMatrix<T> M_W = gsMatrix<T>::Zero(5, 5);
          for (int col = 0; col < 5; ++col) {
            for (typename gsSparseMatrix<T>::InnerIterator it(M_West_edge,
                                                              cols_W[col]);
                 it; ++it) {
              M_W(0, col) += it.value() * eval_val(it.row());
              M_W(1, col) += it.value() * eval_v(it.row());  // tangent
              M_W(2, col) += it.value() * eval_vv(it.row()); // second tangent
              M_W(3, col) += it.value() * eval_u(it.row());  // normal
              M_W(4, col) += it.value() * eval_uv(it.row()); // mixed
            }
          }

          gsMatrix<T> M_S = gsMatrix<T>::Zero(5, 5);
          for (int col = 0; col < 5; ++col) {
            for (typename gsSparseMatrix<T>::InnerIterator it(M_South_edge,
                                                              cols_S[col]);
                 it; ++it) {
              M_S(0, col) += it.value() * eval_val(it.row());
              M_S(1, col) += it.value() * eval_u(it.row());  // tangent
              M_S(2, col) += it.value() * eval_uu(it.row()); // second tangent
              M_S(3, col) += it.value() * eval_v(it.row());  // normal
              M_S(4, col) += it.value() * eval_uv(it.row()); // mixed
            }
          }

          std::vector<std::pair<int, int>> patch_indices = {
              {0, 0}, {1, 0}, {0, 1}, {1, 1}};
          gsMatrix<T> M_P(4, 4);
          for (int col = 0; col < 4; ++col) {
            int j1 = patch_indices[col].first;
            int j2 = patch_indices[col].second;
            M_P(0, col) = val0_mat(j1, 0) * val1_mat(j2, 0); // value
            M_P(1, col) = der0_mat(j1, 0) * val1_mat(j2, 0); // deriv u
            M_P(2, col) = val0_mat(j1, 0) * der1_mat(j2, 0); // deriv v
            M_P(3, col) = der0_mat(j1, 0) * der1_mat(j2, 0); // deriv uv
          }

          for (index_t bfID = 0; bfID < 6; ++bfID) {
            T val_phi = (bfID == 0) ? 1.0 : 0.0;
            T d1_phi = (bfID == 1) ? sigma : 0.0;
            T d2_phi = (bfID == 2) ? sigma : 0.0;
            T d11_phi = (bfID == 3) ? sigma * sigma : 0.0;
            T d22_phi = (bfID == 5) ? sigma * sigma : 0.0;
            T d12_phi = (bfID == 4) ? sigma * sigma : 0.0;

            T f_val = val_phi;
            T f_u = J00 * d1_phi + J10 * d2_phi;
            T f_v = J01 * d1_phi + J11 * d2_phi;
            T f_uu = J00 * J00 * d11_phi + 2 * J00 * J10 * d12_phi +
                     J10 * J10 * d22_phi + d1_phi * X_uu + d2_phi * Y_uu;
            T f_vv = J01 * J01 * d11_phi + 2 * J01 * J11 * d12_phi +
                     J11 * J11 * d22_phi + d1_phi * X_vv + d2_phi * Y_vv;
            T f_uv = J00 * J01 * d11_phi + (J00 * J11 + J10 * J01) * d12_phi +
                     J10 * J11 * d22_phi + d1_phi * X_uv + d2_phi * Y_uv;

            gsVector<T> rhs_W(5);
            rhs_W << f_val, f_v, f_vv, f_u, f_uv;
            gsVector<T> coefs_W = M_W.partialPivLu().solve(rhs_W);

            gsVector<T> rhs_S(5);
            rhs_S << f_val, f_u, f_uu, f_v, f_uv;
            gsVector<T> coefs_S = M_S.partialPivLu().solve(rhs_S);

            gsVector<T> rhs_P(4);
            rhs_P << f_val, f_u, f_v, f_uv;
            gsVector<T> coefs_P = M_P.partialPivLu().solve(rhs_P);

            gsVector<T> local_W = gsVector<T>::Zero(tbRot.size());
            for (int col = 0; col < 5; ++col) {
              for (typename gsSparseMatrix<T>::InnerIterator it(M_West_edge,
                                                                cols_W[col]);
                   it; ++it) {
                local_W(it.row()) += coefs_W(col) * it.value();
              }
            }

            gsVector<T> local_S = gsVector<T>::Zero(tbRot.size());
            for (int col = 0; col < 5; ++col) {
              for (typename gsSparseMatrix<T>::InnerIterator it(M_South_edge,
                                                                cols_S[col]);
                   it; ++it) {
                local_S(it.row()) += coefs_S(col) * it.value();
              }
            }

            gsVector<T> local_P = gsVector<T>::Zero(tbRot.size());
            for (int col = 0; col < 4; ++col) {
              int r = act0(patch_indices[col].first, 0) +
                      act1(patch_indices[col].second, 0) * tbRot.size(0);
              local_P(r) += coefs_P(col);
            }

            for (int i = 0; i < tbRot.size(); ++i) {
              T val = local_W(i) + local_S(i) - local_P(i);
              if (std::abs(val) > 1e-12) {
                triplets.push_back(gsEigen::Triplet<T>(offsets[pk] + perm[k][i],
                                                       v_col + bfID, val));
              }
            }
          }
        }
      }
    }
    // 5. Fill global sparse matrix
    m_M.resize(total_broken, current_col);
    m_M.setFromTriplets(triplets.begin(), triplets.end());
  }
};

} // namespace gismo
