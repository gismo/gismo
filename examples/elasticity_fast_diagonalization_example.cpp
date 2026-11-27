/** @file elasticity_fast_diagonalization_example.cpp

  @brief Linear elasticity on one patch, solved by conjugate gradients with
  the Jacobi preconditioner and with the fast diagonalization of
  gsPatchPreconditionersCreator::elasticityFastDiagonalizationOp

  Linear elasticity on one patch with G+Smo's expression assembler, solved by
  conjugate gradients with (a) the Jacobi preconditioner and (b)
  gsPatchPreconditionersCreator::elasticityFastDiagonalizationOp, without a
  geometry (fd), with one scalar per direction from the anchor averages of
  the metric diagonal (fd_geo), with the Gauss-weighted
  means of mu G_kk + (lambda+mu) detJ (dxi_k/dx_c)^2 (fd_kern), and with the
  separable weights the operator builds from a geometry (fd_sep). Prints
  iteration counts and times, and checks the per-direction operator against
  an independently assembled Kronecker sum and the geometry-weighted scalar
  operator against the assembled Laplace matrix of a graded box.

  Usage: elasticity_fast_diagonalization_example -d 3 -g 0 -p 3 -r 4
  prints one line per preconditioner with iterations, solve time and
  condition estimate; -g 1 is the 8:1 box, -g 2 the quarter annulus,
  -n Poisson's ratio.

  This file is part of the G+Smo library.

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Author(s): R. Goswami
*/

#include <gismo.h>
#include <cstdio>
using namespace gismo;
typedef real_t T;

static void run(const char *name, const gsMultiPatch<T> &mp, int p, int r, T nu) {
  const short_t d = mp.geoDim();
  const T E = 1, lam = E * nu / ((1 + nu) * (1 - 2 * nu)), mu = E / (2 * (1 + nu));
  gsMultiBasis<T> mb(mp, true);
  mb.setDegree(p);
  for (int i = 0; i < r; ++i) mb.uniformRefine();
  gsBoundaryConditions<T> bc;
  gsConstantFunction<T> zero(gsVector<T>::Zero(d), d);
  bc.addCondition(0, boundary::west, condition_type::dirichlet, &zero, 0, false, -1);
  bc.setGeoMap(mp);
  gsExprAssembler<T> A(1, 1);
  A.setIntegrationDomain(mb.domain());
  gsExprAssembler<T>::geometryMap G = A.getMap(mp);
  gsExprAssembler<T>::space u = A.getSpace(mb, d);
  u.setup(bc, dirichlet::l2Projection, 0);
  A.initSystem();
  gsStopwatch sw;
  auto J = ijac(u, G);
  A.assemble(lam * idiv(u, G) * idiv(u, G).tr() * meas(G) +
             mu * ((J.cwisetr() + J) % J.tr()) * meas(G));
  const double tasm = sw.stop();
  const gsSparseMatrix<T> &K = A.matrix();
  gsMatrix<T> rhs = gsMatrix<T>::Zero(K.rows(), 1), x;
  // deterministic right-hand side with all frequencies
  for (index_t i = 0; i < rhs.rows(); ++i) rhs(i, 0) = std::sin(0.37 * i) + 0.01 * (i % 7);
  gsOptionList opt = gsAssembler<T>::defaultOptions();
  // All components are clamped on the west side
  gsBoundaryConditions<T> bcs;
  gsConstantFunction<T> zero1(0.0, d);
  bcs.addCondition(boundary::west, condition_type::dirichlet, &zero1);
  static const char *names[] = {"jacobi", "fd", "fd_geo", "fd_kern", "fd_sep"};
  for (int which = 0; which < 5; ++which) {
    sw.restart();
    gsLinearOperator<T>::Ptr P;
    if (which == 0) P = makeJacobiOp(K);
    if (which == 1) P = gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp(mb[0], bc, opt, lam, mu);
    if (which == 2) { // one scalar per direction: anchor averages of the metric diagonal
      gsMatrix<T> pts, jacs;
      mp.patch(0).basis().anchors_into(pts);
      mp.patch(0).jacobian_into(pts, jacs);
      gsVector<T> sc = gsVector<T>::Zero(d);
      for (index_t i = 0; i < pts.cols(); ++i) {
        const gsMatrix<T> Jm = jacs.block(0, i * d, d, d), Ji = Jm.inverse();
        for (int k = 0; k < d; ++k) sc[k] += std::abs(Jm.determinant()) * Ji.row(k).squaredNorm();
      }
      sc /= T(pts.cols());
      gsBlockOp<T>::uPtr B = gsBlockOp<T>::make(d, d);
      for (int c = 0; c < d; ++c) {
        gsVector<T> beta = mu * sc;
        beta[c] += (lam + mu) * sc[c];
        B->addOperator(c, c, gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], bcs, opt, 0, beta));
      }
      P = give(B);
    }
    if (which == 3) { // per-direction scalars: Gauss-weighted means over the patch
      gsMatrix<T> acc = gsMatrix<T>::Zero(d, d), nodes, jacs;
      gsVector<T> wts;
      T vol = 0;
      gsQuadRule<T>::uPtr rule = gsQuadrature::getPtr(mb[0], opt);
      memory::shared_ptr<gsDomain<T>> dom = mb[0].domain();
      gsBasis<T>::domainIter it = dom->beginAll(), end = dom->endAll();
      for (; it < end; ++it) {
        rule->mapTo(it.lowerCorner(), it.upperCorner(), nodes, wts);
        mp.patch(0).jacobian_into(nodes, jacs);
        for (index_t i = 0; i < nodes.cols(); ++i) {
          const gsMatrix<T> Jm = jacs.block(0, i * d, d, d), Ji = Jm.inverse();
          acc += wts[i] * std::abs(Jm.determinant()) * Ji.cwiseAbs2();
          vol += wts[i];
        }
      }
      acc /= vol;
      gsBlockOp<T>::uPtr B = gsBlockOp<T>::make(d, d);
      for (int c = 0; c < d; ++c) {
        gsVector<T> beta(d);
        for (int k = 0; k < d; ++k) beta[k] = mu * acc.row(k).sum() + (lam + mu) * acc(k, c);
        B->addOperator(c, c, gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], bcs, opt, 0, beta));
      }
      P = give(B);
    }
    if (which == 4) P = gsPatchPreconditionersCreator<T>::elasticityFastDiagonalizationOp(mb[0], bc, opt, lam, mu, &mp.patch(0));
    const double tset = sw.stop();
    gsConjugateGradient<T> cg(K, P);
    cg.setTolerance(1e-8);
    cg.setMaxIterations(20000);
    cg.setCalcEigenvalues(true);
    x.setZero(K.rows(), 1);
    sw.restart();
    cg.solve(rhs, x);
    const double tsol = sw.stop();
    std::printf("FDEL geo=%s d=%d p=%d r=%d nu=%.3f dofs=%ld asm=%.3f prec=%s setup=%.4f "
                "it=%d solve=%.4f relres=%.2e cond=%.1f true_relres=%.2e\n",
                name, int(d), p, r, double(nu), long(K.rows()), tasm,
                names[which], tset,
                int(cg.iterations()), tsol, double(cg.error()), double(cg.getConditionNumber()),
                double((rhs - K * x).norm() / rhs.norm()));
    std::fflush(stdout);
  }
}

// The per-direction operator must invert the Kronecker sum built from the
// univariate matrices of the same class.
static int exactness() {
  gsMultiPatch<T> mp(*gsNurbsCreator<T>::BSplineCube(1));
  gsMultiBasis<T> mb(mp, true);
  mb.setDegree(3);
  mb.uniformRefine();
  dynamic_cast<gsTensorBSplineBasis<3, T> &>(mb[0]).component(1).uniformRefine();
  gsBoundaryConditions<T> bc, none;
  gsConstantFunction<T> one(1, 3);
  bc.addCondition(boundary::west, condition_type::dirichlet, &one);
  bc.addCondition(boundary::back, condition_type::dirichlet, &one);
  gsOptionList opt = gsAssembler<T>::defaultOptions();
  gsVector<T> beta(3);
  beta << 0.7, 3.0, 11.0;
  const T alpha = 0.25;
  // Kronecker sum from univariate mass and stiffness matrices
  std::vector<gsSparseMatrix<T>> M(3), S(3);
  for (int k = 0; k < 3; ++k) {
    gsBoundaryConditions<T> b1;
    if (k == 0) b1.addCondition(boundary::west, condition_type::dirichlet, &one);
    if (k == 2) b1.addCondition(boundary::east, condition_type::dirichlet, &one);
    M[k] = gsPatchPreconditionersCreator<T>::massMatrix(mb[0].component(k), b1, opt);
    S[k] = gsPatchPreconditionersCreator<T>::stiffnessMatrix(mb[0].component(k), b1, opt, 0, 1);
  }
  gsSparseMatrix<T> Kx = M[2].kron(M[1]).kron(S[0]), Ky = M[2].kron(S[1]).kron(M[0]),
                    Kz = S[2].kron(M[1]).kron(M[0]), Mm = M[2].kron(M[1]).kron(M[0]);
  gsSparseMatrix<T> Aref = beta[0] * Kx + beta[1] * Ky + beta[2] * Kz + alpha * Mm;
  gsLinearOperator<T>::Ptr P = gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], bc, opt, alpha, beta);
  gsMatrix<T> X = gsMatrix<T>::Random(Aref.rows(), 3), Y, Z = Aref * X;
  P->apply(Z, Y);
  const double err = double((Y - X).norm() / X.norm());
  // scalar interface unchanged
  gsLinearOperator<T>::Ptr P1 = gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], bc, opt, alpha, 2.0);
  gsSparseMatrix<T> A1 = gsPatchPreconditionersCreator<T>::stiffnessMatrix(mb[0], bc, opt, alpha, 2.0);
  Z = A1 * X;
  P1->apply(Z, Y);
  const double err1 = double((Y - X).norm() / X.norm());
  std::printf("FDEL_EXACT n=%ld aniso_err=%.2e scalar_err=%.2e\n", long(Aref.rows()), err, err1);
  return (err < 1e-10 && err1 < 1e-10) ? 0 : 1;
}

// On a box with a different grading per direction the metric is diagonal and
// its entries are products of univariate functions, so the scalar operator
// with a geometry must invert the assembled Laplace matrix. Also prints the
// CG iterations of that operator for the Laplace matrix of the annulus.
static int sep_exactness() {
  int bad = 0;
  for (int d = 2; d <= 3; ++d) {
    gsMultiPatch<T> mp;
    if (d == 2) mp.addPatch(*gsNurbsCreator<T>::BSplineSquare(1));
    else mp.addPatch(*gsNurbsCreator<T>::BSplineCube(1));
    mp.patch(0).degreeElevate(1);
    gsMatrix<T> &c = mp.patch(0).coefs();
    c.col(0) = c.col(0).array().square().matrix();
    c.col(1) = (0.5 * (c.col(1).array() + c.col(1).array().square())).matrix();
    if (d == 3) c.col(2) *= 3;
    mp.computeTopology();
    gsMultiBasis<T> mb(mp, true);
    mb.setDegree(3);
    for (int i = 0; i < (d == 2 ? 3 : 2); ++i) mb.uniformRefine();
    gsBoundaryConditions<T> bc;
    gsConstantFunction<T> zero(0.0, d);
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &zero);
    bc.addCondition(0, boundary::north, condition_type::dirichlet, &zero);
    bc.setGeoMap(mp);
    gsExprAssembler<T> A(1, 1);
    A.setIntegrationDomain(mb.domain());
    gsExprAssembler<T>::geometryMap G = A.getMap(mp);
    gsExprAssembler<T>::space u = A.getSpace(mb);
    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G));
    const gsSparseMatrix<T> &K = A.matrix();
    gsOptionList opt = gsAssembler<T>::defaultOptions();
    gsStopwatch sw;
    gsLinearOperator<T>::Ptr P = gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], mp.patch(0), bc, opt);
    const double tset = sw.stop();
    gsMatrix<T> X = gsMatrix<T>::Random(K.rows(), 2), Y, Z = K * X;
    P->apply(Z, Y);
    const double err = double((Y - X).norm() / X.norm());
    std::printf("FDEL_SEP_EXACT d=%d n=%ld err=%.2e setup=%.4f\n", d, long(K.rows()), err, tset);
    if (!(err < 1e-8)) bad = 1;
  }
  { // Laplace on the quarter annulus
    gsMultiPatch<T> mp(*gsNurbsCreator<T>::BSplineFatQuarterAnnulus(1, 2));
    mp.computeTopology();
    gsMultiBasis<T> mb(mp, true);
    mb.setDegree(3);
    for (int i = 0; i < 4; ++i) mb.uniformRefine();
    gsBoundaryConditions<T> bc;
    gsConstantFunction<T> zero(0.0, 2);
    bc.addCondition(0, boundary::west, condition_type::dirichlet, &zero);
    bc.setGeoMap(mp);
    gsExprAssembler<T> A(1, 1);
    A.setIntegrationDomain(mb.domain());
    gsExprAssembler<T>::geometryMap G = A.getMap(mp);
    gsExprAssembler<T>::space u = A.getSpace(mb);
    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G));
    const gsSparseMatrix<T> &K = A.matrix();
    gsOptionList opt = gsAssembler<T>::defaultOptions();
    gsMatrix<T> rhs(K.rows(), 1), x;
    for (index_t i = 0; i < rhs.rows(); ++i) rhs(i, 0) = std::sin(0.37 * i) + 0.01 * (i % 7);
    for (int which = 0; which < 2; ++which) {
      gsLinearOperator<T>::Ptr P = which == 0
          ? gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], bc, opt)
          : gsPatchPreconditionersCreator<T>::fastDiagonalizationOp(mb[0], mp.patch(0), bc, opt);
      gsConjugateGradient<T> cg(K, P);
      cg.setTolerance(1e-8);
      cg.setMaxIterations(2000);
      cg.setCalcEigenvalues(true);
      x.setZero(K.rows(), 1);
      cg.solve(rhs, x);
      std::printf("FDEL_SEP_LAPLACE geo=annulus prec=%s n=%ld it=%d cond=%.2f\n", which == 0 ? "fd" : "fd_sep",
                  long(K.rows()), int(cg.iterations()), double(cg.getConditionNumber()));
    }
  }
  return bad;
}

int main(int argc, char **argv) {
  int p = 2, r = 3, d = 3, geo = 0;
  T nu = 0.3;
  gsCmdLine cmd("fast diagonalization for elasticity");
  cmd.addInt("p", "degree", "spline degree", p);
  cmd.addInt("r", "refine", "uniform refinements", r);
  cmd.addInt("d", "dim", "2 or 3", d);
  cmd.addInt("g", "geo", "0 unit box, 1 stretched box 8:1(:1), 2 quarter annulus", geo);
  cmd.addReal("n", "nu", "Poisson ratio", nu);
  try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
  if (exactness()) return 1;
  if (sep_exactness()) return 1;
  gsMultiPatch<T> mp;
  const char *name = geo == 0 ? "box" : geo == 1 ? "stretched" : "annulus";
  if (geo == 2) {
    gsTensorBSpline<2, T>::uPtr q = gsNurbsCreator<T>::BSplineFatQuarterAnnulus(1, 2);
    if (d == 2) mp.addPatch(*q);
    else mp.addPatch(*gsNurbsCreator<T>::lift3D(*q, 1));
  } else {
    if (d == 2) mp.addPatch(*gsNurbsCreator<T>::BSplineSquare(1));
    else mp.addPatch(*gsNurbsCreator<T>::BSplineCube(1));
    if (geo == 1) {
      gsMatrix<T> &c = mp.patch(0).coefs();
      c.col(0) *= 8;
    }
  }
  mp.computeTopology();
  run(name, mp, p, r, nu);
  return 0;
}
