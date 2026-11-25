/** @file gsBarrierCore_test.cpp

    @brief Tests of the PDE-based parameterization (gsBarrierCore): energy
    derivatives, element-wise assembly, univariate basis tables, the nodal
    algebraic multigrid preconditioner, free interfaces and a degenerate
    initial guess.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
**/

#include "gismo_unittest.h"

#include <gsAssembler/gsDofMapperCreator.h>
#include <gsModeling/gsBarrierCore.h>
#include <gsModeling/gsBarrierCore.hpp>

#include <random>

namespace
{

typedef gsExprAssembler<real_t>::geometryMap geometryMap;
typedef gsExprAssembler<real_t>::space       space;

/// Tensor-product B-spline patch of degree \a deg with \a nElem elements per
/// direction on [0,1]^d whose control points are a smooth bijective
/// deformation of the Greville abscissae.
template<short_t d>
gsMultiPatch<real_t> warpedPatch(const int deg, const int nElem)
{
    const gsKnotVector<real_t> kv(0, 1, nElem - 1, deg + 1);
    std::vector<gsKnotVector<real_t> > kvs(d, kv);
    gsTensorBSplineBasis<d, real_t> basis(kvs);
    const gsMatrix<real_t> g = basis.anchors();
    gsMatrix<real_t> coefs(g.cols(), d);
    for (index_t k = 0; k != g.cols(); ++k)
    {
        const real_t u = g(0, k), v = g(1, k);
        coefs(k, 0) = u + 0.15 * math::sin(EIGEN_PI * v);
        coefs(k, 1) = v + 0.10 * math::sin(EIGEN_PI * u) + 0.2 * u * v;
        if (d == 3)
        {
            coefs(k, 2) = g(2, k) + 0.1 * math::sin(EIGEN_PI * u);
        }
    }
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsTensorBSpline<d, real_t>(basis, give(coefs)));
    mp.computeTopology();
    return mp;
}

/// Polynomial map of degree 2 whose inverse is harmonic: the conformal map
/// z -> z^2 in the plane of the parameters (a,b), affine in the other one.
gsMultiPatch<real_t> conformalPatch(const short_t d, const int a, const int b)
{
    const gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsMultiPatch<real_t> mp;
    if (d == 2)
    {
        gsTensorBSplineBasis<2, real_t> basis(kv, kv);
        const gsMatrix<real_t> g = basis.anchors();
        gsMatrix<real_t> val(2, g.cols());
        for (index_t k = 0; k != g.cols(); ++k)
        {
            const real_t p = 1 + g(a, k), q = g(b, k);
            val(0, k) = p * p - q * q;
            val(1, k) = 2 * p * q;
        }
        mp.addPatch(basis.interpolateAtAnchors(val));
    }
    else
    {
        gsTensorBSplineBasis<3, real_t> basis(kv, kv, kv);
        const gsMatrix<real_t> g = basis.anchors();
        gsMatrix<real_t> val(3, g.cols());
        const int c = 3 - a - b;
        for (index_t k = 0; k != g.cols(); ++k)
        {
            const real_t p = 1 + g(a, k), q = g(b, k);
            val(0, k) = p * p - q * q;
            val(1, k) = 2 * p * q;
            val(2, k) = 0.7 * g(c, k);
        }
        mp.addPatch(basis.interpolateAtAnchors(val));
    }
    mp.computeTopology();
    return mp;
}

/// Same parameterization on the parameter domain [0,1]^2
gsTensorBSpline<2, real_t> onUnitSquare(const gsTensorBSpline<2, real_t>& g)
{
    gsKnotVector<real_t> ku = g.basis().knots(0), kv = g.basis().knots(1);
    ku.affineTransformTo(0, 1);
    kv.affineTransformTo(0, 1);
    return gsTensorBSpline<2, real_t>(gsTensorBSplineBasis<2, real_t>(ku, kv),
                                      g.coefs());
}

real_t relDiff(const gsVector<real_t>& a, const gsVector<real_t>& b)
{
    return (a - b).norm() / b.norm();
}

/// Test setup: all boundary control points fixed, interfaces matched; an
/// expression assembler with Gauss-Lobatto (regularized stage) or Gauss (H2
/// stage) quadrature for the comparisons with the expressions.
struct Setup
{
    gsMultiPatch<real_t> mp, mpS;
    gsMultiBasis<real_t> mb;
    gsDofMapper mapper;
    gsExprAssembler<real_t> A;
    gsExprEvaluator<real_t> ev;
    gsBoundaryConditions<real_t> bc;
    short_t d;

    Setup(const gsMultiPatch<real_t>& m, const bool lobatto)
    : mp(m), mpS(m), mb(m), d(m.parDim())
    {
        mapper = createMapper(mb, d, false);
        for (gsBoxTopology::const_iiterator it = mb.topology().iBegin();
             it != mb.topology().iEnd(); ++it)
        {
            mb.matchInterface(*it, mapper);
        }
        for (gsMultiPatch<real_t>::const_biterator bit = mp.bBegin();
             bit != mp.bEnd(); ++bit)
        {
            const gsMatrix<index_t> idx = mb.basis(bit->patch).boundary(bit->index());
            for (short_t c = 0; c != d; ++c)
            {
                mapper.markBoundary(bit->patch, idx, c);
            }
            bc.addCondition(*bit, condition_type::dirichlet, nullptr);
        }
        mapper.finalize();
        bc.setGeoMap(mp);
        A.setIntegrationDomain(mb.domain());
        if (lobatto)
        {
            A.options().setInt("quRule", gsQuadrature::GaussLobatto);
            A.options().setInt("quB", 2);
        }
        ev = gsExprEvaluator<real_t>(A);
    }

    gsVector<real_t> freeVector() const
    {
        return convertMultiPatchToFreeVector(mp, mapper);
    }

    template<class F>
    gsVector<real_t> vector(const gsVector<real_t>& x, F make)
    {
        convertFreeVectorToMultiPatch<real_t>(x, mapper, mpS);
        geometryMap G = A.getMap(mpS);
        space u = A.getSpace(mb, d);
        u.setup(bc, dirichlet::homogeneous, 0);
        A.initSystem();
        A.assemble(make(u, G));
        return gsVector<real_t>(A.rhs());
    }

    template<class F>
    gsMatrix<real_t> matrix(const gsVector<real_t>& x, F make)
    {
        convertFreeVectorToMultiPatch<real_t>(x, mapper, mpS);
        geometryMap G = A.getMap(mpS);
        space u = A.getSpace(mb, d);
        u.setup(bc, dirichlet::homogeneous, 0);
        A.initSystem();
        A.assemble(make(u, G));
        return gsMatrix<real_t>(A.matrix());
    }

    real_t energy(const gsVector<real_t>& x, const real_t eps, const real_t orient)
    {
        convertFreeVectorToMultiPatch<real_t>(x, mapper, mpS);
        geometryMap G = ev.getMap(mpS);
        space u = A.getSpace(mb, d);
        u.setup(bc, dirichlet::homogeneous, 0);
        return ev.integral(garanzhaWinslowEnergy(u, G, eps, orient));
    }
};

/// Derivatives of the energy integrands with respect to J (both energy forms,
/// all Hessian types) against central differences
template<short_t d>
void checkEnergyIntegrand(const real_t orient, const bool folded)
{
    typedef gsEigen::Matrix<real_t, d, d> MatD;
    typedef gsEigen::Matrix<real_t, d * d, 1> VecJ;
    typedef gsEigen::Matrix<real_t, d * d, d * d> MatJ;
    const index_t gaussNewton = PDEEnergyHessian::GaussNewton;
    const index_t exact = PDEEnergyHessian::Exact;
    const index_t projected = PDEEnergyHessian::ProjectedNewton;
    const real_t eps = 0.05, h = 1e-6;
    std::srand(17);
    for (int trial = 0; trial != 10; ++trial)
    {
        MatD J = MatD::Identity() + 0.5 * MatD::Random();
        if (folded && trial % 2) { J.row(0) *= -1; }
        for (index_t form = 0; form != 2; ++form)
        {
            VecJ g;
            MatJ H, HGN, HP;
            garanzhaEnergyDerivatives<d, real_t>(J, orient, eps, form, exact, &g, &H);
            garanzhaEnergyDerivatives<d, real_t>(J, orient, eps, form, gaussNewton, nullptr, &HGN);
            garanzhaEnergyDerivatives<d, real_t>(J, orient, eps, form, projected, nullptr, &HP);
            // Gauss-Newton residual: vec(J or cof J) / chi^(p/2)
            const bool inverseDirichlet = (form == PDEEnergyForm::InverseDirichlet && d == 3);
            const real_t p = inverseDirichlet ? real_t(1) : real_t(2) / d;
            MatJ Hfd, Jr;
            real_t errG = 0;
            for (int k = 0; k != d * d; ++k)
            {
                MatD Jp = J, Jm = J;
                Jp(k / d, k % d) += h;
                Jm(k / d, k % d) -= h;
                const real_t fp = garanzhaEnergyDerivatives<d, real_t>(Jp, orient, eps, form, gaussNewton, nullptr, nullptr);
                const real_t fm = garanzhaEnergyDerivatives<d, real_t>(Jm, orient, eps, form, gaussNewton, nullptr, nullptr);
                errG = std::max(errG, math::abs((fp - fm) / (2 * h) - g(k)));
                VecJ gp, gm;
                garanzhaEnergyDerivatives<d, real_t>(Jp, orient, eps, form, gaussNewton, &gp, nullptr);
                garanzhaEnergyDerivatives<d, real_t>(Jm, orient, eps, form, gaussNewton, &gm, nullptr);
                Hfd.col(k) = (gp - gm) / (2 * h);
                VecJ rp, rm;
                for (int s = 0; s != 2; ++s)
                {
                    const MatD& M = s ? Jm : Jp;
                    MatD cof;
                    cofactorMatrix(M, cof);
                    const real_t D = orient * M.row(0).dot(cof.row(0));
                    const MatD G = inverseDirichlet ? cof : M;
                    VecJ r;
                    for (int a = 0; a != d; ++a)
                        for (int b = 0; b != d; ++b) { r(a * d + b) = G(a, b); }
                    (s ? rm : rp) = r / math::pow(garanzhaChi(D, eps), p / 2);
                }
                Jr.col(k) = (rp - rm) / (2 * h);
            }
            CHECK(errG <= 1e-7 * (1 + g.cwiseAbs().maxCoeff()));
            CHECK((Hfd - H).cwiseAbs().maxCoeff() <= 1e-6 * (1 + H.cwiseAbs().maxCoeff()));
            const MatJ GN = 2 * Jr.transpose() * Jr;
            CHECK((GN - HGN).cwiseAbs().maxCoeff() <= 1e-7 * (1 + HGN.cwiseAbs().maxCoeff()));
            gsEigen::SelfAdjointEigenSolver<MatJ> es(HP);
            CHECK(es.eigenvalues().minCoeff()
                  >= -1e-12 * (1 + es.eigenvalues().cwiseAbs().maxCoeff()));
        }
    }
}

/// Element-wise assembly (gsRegularizedPDEAssembler) against the expressions
template<short_t d>
void checkElementAssembly(const gsMultiPatch<real_t>& in, const real_t pert,
                          const real_t orient)
{
    Setup s(in, true);
    std::srand(11);
    gsVector<real_t> x = s.freeVector();
    x += pert * gsVector<real_t>::Random(x.size());
    const real_t eps = 1e-2;
    gsRegularizedPDEAssembler<d, real_t> fem(s.mp, s.mapper, s.A.options(), orient);
    fem.setFreeVector(x);

    const real_t e0 = s.energy(x, eps, orient);
    CHECK(math::abs(fem.energy(eps) - e0) <= 1e-12 * math::abs(e0));

    gsVector<real_t> g;
    gsSparseMatrix<real_t> Hl;
    fem.gradientHessian(eps, g, Hl);
    CHECK(relDiff(g, s.vector(x, [&](space u, geometryMap G)
                  { return garanzhaWinslowGrad(u, G, eps, orient); })) <= 1e-12);
    const gsMatrix<real_t> H0 = s.matrix(x, [&](space u, geometryMap G)
                                         { return garanzhaWinslowHess(u, G, eps, orient); });
    const gsMatrix<real_t> H1(Hl);
    gsMatrix<real_t> Hfull = H1 + gsMatrix<real_t>(H1.transpose());
    Hfull.diagonal() = H1.diagonal();
    CHECK((Hfull - H0).norm() <= 1e-10 * H0.norm());
    CHECK(gsMatrix<real_t>(H1.template triangularView<gsEigen::StrictlyUpper>()).norm() == 0);

    // H2 stage: Gauss rule
    {
        Setup sg(in, false);
        gsRegularizedPDEAssembler<d, real_t> h2(sg.mp, sg.mapper,
                                                gsExprAssembler<real_t>::defaultOptions(), 1);
        h2.setFreeVector(x);
        gsVector<real_t> rh;
        h2.h2Residual(rh);
        CHECK(relDiff(rh, sg.vector(x, [](space u, geometryMap G)
                      { return scaledEGGRes(u, G); })) <= 1e-11);
        gsSparseMatrix<real_t> Jh;
        h2.h2Jacobian(0, Jh);
        const gsMatrix<real_t> J0 = sg.matrix(x, [](space u, geometryMap G)
                                              { return scaledEGGJac(u, G, 0); });
        CHECK((gsMatrix<real_t>(Jh) - J0).norm() <= 1e-11 * J0.norm());
    }

    gsVector<real_t> r;
    fem.h1Residual(eps, r);
    CHECK(relDiff(r, s.vector(x, [&](space u, geometryMap G)
                  { return regularizedH1Res(u, G, eps, orient); })) <= 1e-12);
    gsSparseMatrix<real_t> J;
    fem.h1Jacobian(eps, false, J);
    const gsMatrix<real_t> J0 = s.matrix(x, [&](space u, geometryMap G)
                                         { return regularizedH1Jac(u, G, eps, orient, false); });
    CHECK((gsMatrix<real_t>(J) - J0).norm() <= 1e-12 * J0.norm());

    // H1 Jacobians of both forms against central differences of the residual
    for (index_t q = 1; q <= 2; ++q)
    {
        fem.setH1Power(q);
        fem.setFreeVector(x);
        gsSparseMatrix<real_t> Jq;
        fem.h1Jacobian(eps, false, Jq);
        const gsVector<real_t> v = gsVector<real_t>::Random(x.size());
        const gsVector<real_t> an = Jq * v;
        real_t best = std::numeric_limits<real_t>::max();
        for (real_t h : {1e-4, 1e-5, 1e-6})
        {
            gsVector<real_t> rp, rm;
            fem.setFreeVector(x + h * v);
            fem.h1Residual(eps, rp);
            fem.setFreeVector(x - h * v);
            fem.h1Residual(eps, rm);
            best = std::min(best, relDiff((rp - rm) / (2 * h), an));
        }
        CHECK(best <= 1e-6);
    }
}

/// Univariate basis tables against gismo's evaluation of the basis
template<short_t d>
void checkTables(const gsMultiPatch<real_t>& in, const gsOptionList& quadrature,
                 const std::vector<index_t>& extraNodes = std::vector<index_t>())
{
    Setup s(in, false);
    gsRegularizedPDEAssembler<d, real_t> fem(s.mp, s.mapper, quadrature, 1, extraNodes);
    real_t errV = 0, errD = 0, errD2 = 0;
    index_t tabulated = 0;
    for (index_t i = 0; i != fem.numActiveElements(); ++i)
    {
        size_t patch;
        gsMatrix<real_t> nodes, v, d1, d2;
        gsMatrix<index_t> act;
        bool tab;
        fem.elementBasisData(i, patch, nodes, act, v, d1, d2, tab);
        tabulated += tab;
        const gsBasis<real_t>& basis = s.mb.basis(patch);
        std::vector<gsMatrix<real_t> > ref;
        for (index_t k = 0; k != nodes.cols(); ++k)
        {
            gsMatrix<index_t> a;
            basis.active_into(nodes.col(k), a);
            if (a != act) { continue; } // node on an element boundary
            basis.evalAllDers_into(nodes.col(k), 2, ref);
            errV = std::max(errV, (ref[0] - v.col(k)).cwiseAbs().maxCoeff());
            errD = std::max(errD, (ref[1] - d1.col(k)).cwiseAbs().maxCoeff()
                                  / (1 + ref[1].cwiseAbs().maxCoeff()));
            errD2 = std::max(errD2, (ref[2] - d2.col(k)).cwiseAbs().maxCoeff()
                                    / (1 + ref[2].cwiseAbs().maxCoeff()));
        }
    }
    CHECK_EQUAL(fem.numActiveElements(), tabulated);
    CHECK(errV <= 1e-13);
    CHECK(errD <= 1e-13);
    CHECK(errD2 <= 1e-12);
}

/// Number of grid points (n per direction) with det(J) <= 0
index_t nonPositiveJacobians(const gsGeometry<real_t>& g, const index_t n)
{
    const short_t d = g.parDim();
    const gsMatrix<real_t> box = g.support();
    gsVector<real_t> lo = box.col(0), hi = box.col(1);
    gsVector<unsigned> np(d);
    np.setConstant(n);
    const gsMatrix<real_t> pts = gsPointGrid(lo, hi, np);
    gsMatrix<real_t> der;
    g.deriv_into(pts, der);
    index_t count = 0;
    for (index_t k = 0; k != pts.cols(); ++k)
    {
        count += !(der.reshapeCol(k, d, d).determinant() > 0);
    }
    return count;
}

} // namespace

SUITE(gsBarrierCore_test)
{
    TEST(energyIntegrandDerivatives)
    {
        checkEnergyIntegrand<2>(1, false);
        checkEnergyIntegrand<2>(-1, true);
        checkEnergyIntegrand<3>(1, false);
        checkEnergyIntegrand<3>(1, true);
    }

    TEST(elementAssemblyVsExpressions)
    {
        const gsMultiPatch<real_t> a = warpedPatch<2>(3, 6), b = warpedPatch<3>(2, 3);
        checkElementAssembly<2>(a, 0.01, 1);
        checkElementAssembly<2>(a, 0.3, -1);  // folded
        checkElementAssembly<3>(b, 0.02, 1);
        checkElementAssembly<3>(b, 0.3, 1);   // folded
    }

    TEST(assembledEnergyDerivatives)
    {
        Setup s(warpedPatch<2>(3, 5), true);
        std::srand(19);
        gsVector<real_t> x = s.freeVector();
        x += 0.01 * gsVector<real_t>::Random(x.size());
        gsRegularizedPDEAssembler<2, real_t> fem(s.mp, s.mapper, s.A.options(), 1);
        fem.setFreeVector(x);
        const real_t eps = 1e-2;
        gsVector<real_t> g0, g1;
        gsSparseMatrix<real_t> H0, H1;
        fem.gradientHessian(eps, g0, H0);  // Gauss-Newton fast path
        fem.energyGradientHessian(eps, PDEEnergyHessian::GaussNewton, g1, H1);
        CHECK(relDiff(g1, g0) <= 1e-12);
        CHECK(gsMatrix<real_t>(H1 - H0).norm() <= 1e-12 * gsMatrix<real_t>(H0).norm());
        for (index_t form = 0; form != 2; ++form)
        {
            fem.setEnergyForm(form);
            fem.setFreeVector(x);
            gsVector<real_t> g;
            gsSparseMatrix<real_t> Hl, tmp;
            fem.energyGradientHessian(eps, PDEEnergyHessian::Exact, g, Hl);
            gsMatrix<real_t> Hfull(Hl);
            Hfull += gsMatrix<real_t>(Hl.transpose());
            Hfull.diagonal() /= 2;
            const gsVector<real_t> v = gsVector<real_t>::Random(x.size());
            const gsVector<real_t> Hv = Hfull * v;
            real_t bestG = std::numeric_limits<real_t>::max(), bestH = bestG;
            for (real_t h : {1e-4, 1e-5, 1e-6})
            {
                gsVector<real_t> gp, gm;
                fem.setFreeVector(x + h * v);
                const real_t ep = fem.energy(eps);
                fem.energyGradientHessian(eps, PDEEnergyHessian::GaussNewton, gp, tmp);
                fem.setFreeVector(x - h * v);
                const real_t em = fem.energy(eps);
                fem.energyGradientHessian(eps, PDEEnergyHessian::GaussNewton, gm, tmp);
                bestG = std::min(bestG, math::abs((ep - em) / (2 * h) - g.dot(v))
                                        / math::abs(g.dot(v)));
                bestH = std::min(bestH, relDiff((gp - gm) / (2 * h), Hv));
            }
            CHECK(bestG <= 1e-7);
            CHECK(bestH <= 1e-6);
        }
    }

    TEST(univariateTables)
    {
        gsOptionList lobatto = gsExprAssembler<real_t>::defaultOptions(), gauss = lobatto;
        lobatto.setInt("quRule", gsQuadrature::GaussLobatto);
        lobatto.setInt("quB", 2);
        // non-uniform knots and a C0 knot
        gsMultiPatch<real_t> a = warpedPatch<2>(3, 4);
        gsTensorBSpline<2, real_t>& g = static_cast<gsTensorBSpline<2, real_t>&>(a.patch(0));
        g.insertKnot(0.137, 0);
        g.insertKnot(0.61, 1, 2);
        g.insertKnot(0.42, 0, 3);
        checkTables<2>(a, lobatto);
        checkTables<2>(a, gauss);
        std::vector<index_t> extra(100, 0);
        for (size_t i = 0; i < extra.size(); i += 3) { extra[i] = 2; }
        checkTables<2>(a, lobatto, extra);
        const gsMultiPatch<real_t> b = warpedPatch<3>(2, 3);
        checkTables<3>(b, lobatto);
        checkTables<3>(b, gauss);
    }

    TEST(exactSolutions)
    {
        // maps with harmonic inverse solve the H2 equations and the weak
        // Laplace form of the H1 equations exactly
        for (short_t d = 2; d <= 3; ++d)
        {
            const gsMultiPatch<real_t> mp = conformalPatch(d, 0, 1);
            Setup s(mp, false);
            const gsVector<real_t> x = s.freeVector();
            const gsVector<real_t> rH2 = s.vector(x, [](space u, geometryMap G)
                                                  { return scaledEGGRes(u, G); });
            CHECK(rH2.norm() <= 1e-12);
            gsVector<real_t> r;
            if (d == 2)
            {
                gsRegularizedPDEAssembler<2, real_t> fem(s.mp, s.mapper, s.A.options(), 1);
                fem.setH1Power(1);
                fem.setFreeVector(x);
                fem.h1Residual(1e-12, r);
            }
            else
            {
                gsRegularizedPDEAssembler<3, real_t> fem(s.mp, s.mapper, s.A.options(), 1);
                fem.setH1Power(1);
                fem.setFreeVector(x);
                fem.h1Residual(1e-12, r);
            }
            CHECK(r.norm() <= 1e-12);
        }
    }

    TEST(nodalAMG)
    {
        const short_t d = 2;
        Setup s(warpedPatch<d>(3, 24), true);  // > 400 unknowns: several levels
        const gsVector<real_t> x = s.freeVector();
        gsRegularizedPDEAssembler<d, real_t> fem(s.mp, s.mapper, s.A.options(), 1);
        fem.setFreeVector(x);
        gsSparseMatrix<real_t> A;
        fem.h1Jacobian(1e-2, false, A);
        const index_t n = A.rows();
        const std::vector<index_t> nodes = gsBarrierCore<d, real_t>::nodeOfDofs(s.mp, s.mapper);
        const std::vector<index_t> comps = gsBarrierCore<d, real_t>::componentOfDofs(s.mp, s.mapper);
        const index_t numNodes = *std::max_element(nodes.begin(), nodes.end()) + 1;
        CHECK_EQUAL(n, d * numNodes);

        std::srand(13);
        const gsVector<real_t> b = gsVector<real_t>::Random(n);
        memory::shared_ptr<preAApp::NodalAMG<real_t> > amg(new preAApp::NodalAMG<real_t>());
        amg->setNodes(nodes, comps);
        amg->compute(A);
        CHECK(amg->numLevels() > 1);
        gsGMRes<real_t> gmres(A, amg);
        gmres.setTolerance(1e-10);
        gmres.setMaxIterations(100);
        gsMatrix<real_t> y = gsMatrix<real_t>::Zero(n, 1);
        gmres.solve(b, y);
        CHECK(gmres.iterations() < 60);
        CHECK((A * y - b).norm() <= 1e-8 * b.norm());

        gsMatrix<real_t> z;
        amg->apply(b, z);

        // equivariance under a renumbering of the unknowns
        std::vector<index_t> perm(n);
        for (index_t i = 0; i != n; ++i) { perm[i] = i; }
        std::shuffle(perm.begin(), perm.end(), std::mt19937(5));
        std::vector<index_t> nodesP(n), compsP(n);
        gsVector<real_t> bP(n);
        for (index_t i = 0; i != n; ++i)
        {
            nodesP[perm[i]] = nodes[i];
            compsP[perm[i]] = comps[i];
            bP(perm[i]) = b(i);
        }
        std::vector<gsEigen::Triplet<real_t, index_t> > tr;
        for (index_t j = 0; j < A.outerSize(); ++j)
            for (gsSparseMatrix<real_t>::InnerIterator it(A, j); it; ++it)
            {
                tr.emplace_back(perm[it.row()], perm[it.col()], it.value());
            }
        gsSparseMatrix<real_t> AP(n, n);
        AP.setFromTriplets(tr.begin(), tr.end());
        preAApp::NodalAMG<real_t> amgP;
        amgP.setNodes(nodesP, compsP);
        amgP.compute(AP);
        gsMatrix<real_t> zP;
        amgP.apply(bP, zP);
        real_t err = 0;
        for (index_t i = 0; i != n; ++i)
        {
            err = std::max(err, math::abs(zP(perm[i], 0) - z(i, 0)));
        }
        CHECK(err <= 1e-12 * z.cwiseAbs().maxCoeff());

        // equivariance under a rotation of the physical coordinates:
        // A' = A Q with Q = diag(R, ..., R), then M(A') = Q^T M(A)
        const real_t angle = 0.7;
        gsMatrix<real_t> R(2, 2);
        R << math::cos(angle), -math::sin(angle), math::sin(angle), math::cos(angle);
        std::vector<std::vector<index_t> > dofOf(numNodes, std::vector<index_t>(d));
        for (index_t i = 0; i != n; ++i) { dofOf[nodes[i]][comps[i]] = i; }
        tr.clear();
        for (index_t I = 0; I != numNodes; ++I)
            for (index_t r = 0; r != d; ++r)
                for (index_t k = 0; k != d; ++k)
                {
                    tr.emplace_back(dofOf[I][r], dofOf[I][k], R(r, k));
                }
        gsSparseMatrix<real_t> Q(n, n);
        Q.setFromTriplets(tr.begin(), tr.end());
        const gsSparseMatrix<real_t> AQ = A * Q;
        preAApp::NodalAMG<real_t> amgQ;
        amgQ.setNodes(nodes, comps);
        amgQ.compute(AQ);
        gsMatrix<real_t> zQ;
        amgQ.apply(b, zQ);
        const gsVector<real_t> expected = Q.transpose() * z.col(0);
        CHECK(relDiff(zQ.col(0), expected) <= 1e-10);
    }

    TEST(freeInterfaces)
    {
        // one patch split at u = 1/4; the right part is kept, or rotated by
        // 90 / 180 degrees in its parameter domain; all parameter domains are
        // [0,1]^2. The result must not depend on these choices.
        gsMultiPatch<real_t> full = warpedPatch<2>(3, 8);
        const gsTensorBSpline<2, real_t>& g = static_cast<const gsTensorBSpline<2, real_t>&>(full.patch(0));
        gsTensorBSpline<2, real_t> left, right;
        g.splitAt(0, 0.25, left, right);
        left = onUnitSquare(left);

        gsMatrix<real_t> pts(2, 121), ref;
        for (int j = 0; j != 11; ++j)
            for (int i = 0; i != 11; ++i)
            {
                pts(0, j * 11 + i) = i / 10.0;
                pts(1, j * 11 + i) = j / 10.0;
            }
        for (int variant = 0; variant != 3; ++variant)
        {
            gsTensorBSpline<2, real_t> r = onUnitSquare(right);
            if (variant == 1) { r.reverse(0); r.swapDirections(0, 1); }
            if (variant == 2) { r.reverse(0); r.reverse(1); }
            gsMultiPatch<real_t> mp;
            mp.addPatch(left);
            mp.addPatch(r);
            mp.computeTopology();
            // interior control points moved towards the centroid
            for (size_t p = 0; p != 2; ++p)
            {
                gsMatrix<real_t>& C = mp.patch(p).coefs();
                const gsMatrix<index_t> bnd = mp.basis(p).allBoundary();
                std::vector<bool> onBoundary(C.rows(), false);
                for (index_t i = 0; i != bnd.size(); ++i) { onBoundary[bnd(i, 0)] = true; }
                const gsMatrix<real_t> c = C.colwise().mean();
                for (index_t i = 0; i != C.rows(); ++i)
                {
                    if (!onBoundary[i]) { C.row(i) = 0.5 * (C.row(i) + c); }
                }
            }
            if (variant == 1)
            {
                // the interface has the same parametric length on both sides
                const gsMultiPatch<real_t> scaled = gsBarrierCore<2, real_t>::chordScaled(mp);
                const boundaryInterface& bi = *scaled.topology().iBegin();
                const index_t a = 1 - bi.first().direction();
                const index_t b = bi.dirMap(bi.first())(a);
                const real_t l1 = scaled.basis(bi.first().patch).support()(a, 1);
                const real_t l2 = scaled.basis(bi.second().patch).support()(b, 1);
                CHECK(math::abs(l1 - l2) <= 1e-14 * l1);
            }
            gsDofMapper uncoupled = createMapper(gsMultiBasis<real_t>(mp), 2, false);
            uncoupled.finalize();
            const bool coupled = gsBarrierCore<2, real_t>::couplesPatches(mp, uncoupled);
            CHECK(!coupled);

            gsBarrierPatch<2, real_t> opt(mp, false);
            opt.options().setInt("Verbose", 0);
            opt.options().setReal("PDEEnergyTol", 1e-13);
            opt.compute();
            const gsMultiPatch<real_t>& res = opt.result();
            for (size_t p = 0; p != res.nPatches(); ++p)
            {
                CHECK_EQUAL(0, nonPositiveJacobians(res.patch(p), 30));
            }
            gsMatrix<real_t> xl;
            res.patch(0).eval_into(pts, xl);
            if (variant == 0) { ref = xl; }
            else { CHECK((xl - ref).cwiseAbs().maxCoeff() <= 1e-5); }
        }
    }

    TEST(degenerateInitialGuess)
    {
        // gsBarrierCore::compute directly (no rescaling): square [-1,1]^2,
        // all interior control points exactly at the origin, J = 0 exactly in
        // the elements without boundary control points
        const gsKnotVector<real_t> kv(-1, 1, 10, 4);
        gsTensorBSplineBasis<2, real_t> basis(kv, kv);
        gsMatrix<real_t> C = basis.anchors().transpose();
        const index_t n = basis.size(0);
        for (index_t j = 1; j + 1 < n; ++j)
            for (index_t i = 1; i + 1 < n; ++i) { C.row(j * n + i).setZero(); }
        const gsMultiPatch<real_t> mp(gsTensorBSpline<2, real_t>(basis, C));
        gsMultiBasis<real_t> mb(mp);
        gsDofMapper mapper = createMapper(mb, 2, false);
        for (index_t c = 0; c != 2; ++c) { mapper.markBoundary(0, basis.allBoundary(), c); }
        mapper.finalize();
        gsOptionList opt = gsBarrierCore<2, real_t>::defaultOptions();
        opt.setSwitch("PDESkipH2", true);
        opt.setInt("Verbose", 0);
        const gsMultiPatch<real_t> res = gsBarrierCore<2, real_t>::compute(mp, mapper, opt);
        CHECK_EQUAL(0, nonPositiveJacobians(res.patch(0), 100));
    }
}
