/** @file gsBarrierCore_test.cpp

    @brief Tests of the PDE-based parameterization (gsBarrierCore): energy
    derivatives, element-wise assembly, univariate basis tables, the linear
    solvers and preconditioners, Anderson acceleration, free interfaces and
    degenerate or differently oriented input.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
**/

#include "gismo_unittest.h"

#include <gsAssembler/gsDofMapperCreator.h>
#include <gsModeling/gsBarrierCore.h>
// gsRegularizedPDEAssembler, garanzhaEnergyDerivatives and the bodies of the
// preconditioners are tested directly
#include <gsModeling/gsBarrierCore.hpp>

#include <random>

// Tolerances used below:
// - central differences, the best of the steps 1e-4, 1e-5, 1e-6: truncation
//   error O(h^2) and round-off O(1e-16 / h) leave a relative error of about
//   1e-8 to 1e-9; checked against 1e-6 (derivatives of residuals) or 1e-7;
// - the same quantity computed in two ways (element-wise vs expressions,
//   tables vs evalAllDers): round-off only, 1e-15 to 1e-13; checked against
//   1e-12 (1e-10 for Hessians of folded maps, where large terms cancel);
// - solver results: the stopping tolerances of the solvers involved.

namespace
{

typedef gsExprAssembler<real_t>::geometryMap geometryMap;
typedef gsExprAssembler<real_t>::space       space;

/// Relative difference |a - b| / |b|; infinity if a or b is not finite, so
/// that NaN fails every check
real_t relDiff(const gsVector<real_t>& a, const gsVector<real_t>& b)
{
    if (!a.allFinite() || !b.allFinite())
        return std::numeric_limits<real_t>::infinity();
    return (a - b).norm() / b.norm();
}

/// Maximum of \a err and \a value; infinity if \a value is not finite
/// (std::max(err, NaN) would return err)
real_t maxError(const real_t err, const real_t value)
{
    return math::isfinite(value) ? std::max(err, value)
                                 : std::numeric_limits<real_t>::infinity();
}

/// Tensor-product B-spline patch of degree \a degree with \a numElements
/// elements per direction on [0,1]^d whose control points are a smooth
/// bijective deformation of the Greville abscissae
template<short_t d>
gsMultiPatch<real_t> warpedPatch(const short_t degree, const index_t numElements)
{
    const gsKnotVector<real_t> kv(0, 1, numElements - 1, degree + 1);
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
            coefs(k, 2) = g(2, k) + 0.1 * math::sin(EIGEN_PI * u);
    }
    gsMultiPatch<real_t> mp;
    mp.addPatch(gsTensorBSpline<d, real_t>(basis, give(coefs)));
    mp.computeTopology();
    return mp;
}

/// Polynomial map of degree 2 whose inverse is harmonic: the conformal map
/// z -> z^2 in the plane of the parameters (a,b), affine in the other one
gsMultiPatch<real_t> conformalPatch(const short_t d, const index_t a, const index_t b)
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
        const index_t c = 3 - a - b;
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

/// Moves the control points off the boundary of every patch towards the
/// centroid of the patch's control points by the fraction \a factor
/// (1: all at the centroid)
void collapseInterior(gsMultiPatch<real_t>& mp, const real_t factor)
{
    for (size_t p = 0; p != mp.nPatches(); ++p)
    {
        gsMatrix<real_t>& C = mp.patch(p).coefs();
        const gsMatrix<index_t> boundary = mp.basis(p).allBoundary();
        std::vector<bool> onBoundary(C.rows(), false);
        for (index_t i = 0; i != boundary.size(); ++i)
            onBoundary[boundary(i, 0)] = true;
        const gsMatrix<real_t> centroid = C.colwise().mean();
        for (index_t i = 0; i != C.rows(); ++i)
            if (!onBoundary[i])
                C.row(i) = (1 - factor) * C.row(i) + factor * centroid;
    }
}

/// Number of points of a grid with n points per direction at which det(J)
/// does not have the sign of the patch (the sign of its mean value) or is 0
index_t inconsistentJacobians(const gsGeometry<real_t>& g, const index_t n)
{
    const short_t d = g.parDim();
    const gsMatrix<real_t> box = g.support();
    gsVector<real_t> lo = box.col(0), hi = box.col(1);
    gsVector<unsigned> np(d);
    np.setConstant(n);
    const gsMatrix<real_t> pts = gsPointGrid(lo, hi, np);
    gsMatrix<real_t> der;
    g.deriv_into(pts, der);
    std::vector<real_t> det(pts.cols());
    real_t mean = 0;
    for (index_t k = 0; k != pts.cols(); ++k)
    {
        det[k] = der.reshapeCol(k, d, d).determinant();
        mean += det[k];
    }
    const real_t sign = (mean < 0 ? -1 : 1);
    index_t count = 0;
    for (const real_t D : det)
        count += !(sign * D > 0);
    return count;
}

/// Grid of n x n points on [0,1]^2
gsMatrix<real_t> unitGrid(const index_t n)
{
    gsMatrix<real_t> pts(2, n * n);
    for (index_t j = 0; j != n; ++j)
        for (index_t i = 0; i != n; ++i)
        {
            pts(0, j * n + i) = i / real_t(n - 1);
            pts(1, j * n + i) = j / real_t(n - 1);
        }
    return pts;
}

/// Size of the bounding box of the control points of \a mp
real_t modelSize(const gsMultiPatch<real_t>& mp)
{
    gsMatrix<real_t> box;
    mp.boundingBox(box);
    return (box.col(1) - box.col(0)).maxCoeff();
}

/// Test setup: all boundary control points fixed, interfaces matched; an
/// expression assembler with Gauss-Lobatto (regularized stage) or Gauss (H2
/// stage) quadrature for the comparisons with the expressions
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
                mapper.markBoundary(bit->patch, idx, c);
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

    real_t energy(const gsVector<real_t>& x, const real_t eps, const real_t orientation)
    {
        convertFreeVectorToMultiPatch<real_t>(x, mapper, mpS);
        geometryMap G = ev.getMap(mpS);
        space u = A.getSpace(mb, d);
        u.setup(bc, dirichlet::homogeneous, 0);
        return ev.integral(garanzhaWinslowEnergy(u, G, eps, orientation));
    }
};

/// Derivatives of the energy integrands with respect to J (both energy forms,
/// all Hessian types) against central differences
template<short_t d>
void checkEnergyIntegrand(const real_t orientation, const bool folded)
{
    typedef gsEigen::Matrix<real_t, d, d> MatD;
    typedef gsEigen::Matrix<real_t, d * d, 1> VecJ;
    typedef gsEigen::Matrix<real_t, d * d, d * d> MatJ;
    typedef gsBarrierCore<d, real_t> Core;
    const index_t gaussNewton = Core::EnergyHessian::GaussNewton;
    const index_t exact = Core::EnergyHessian::Exact;
    const index_t projected = Core::EnergyHessian::ProjectedNewton;
    const real_t eps = 0.05, h = 1e-6;
    std::srand(17);
    for (index_t trial = 0; trial != 10; ++trial)
    {
        MatD J = MatD::Identity() + 0.5 * MatD::Random();
        if (folded && trial % 2)
            J.row(0) *= -1;
        for (index_t form = 0; form != 2; ++form)
        {
            VecJ g;
            MatJ H, HGN, HP;
            garanzhaEnergyDerivatives<d, real_t>(J, orientation, eps, form, exact, &g, &H);
            garanzhaEnergyDerivatives<d, real_t>(J, orientation, eps, form, gaussNewton, nullptr, &HGN);
            garanzhaEnergyDerivatives<d, real_t>(J, orientation, eps, form, projected, nullptr, &HP);
            // Gauss-Newton residual: vec(J or cof J) / chi^(p/2)
            const bool inverseDirichlet = (form == Core::EnergyForm::InverseDirichlet && d == 3);
            const real_t p = inverseDirichlet ? real_t(1) : real_t(2) / d;
            MatJ Hfd, Jr;
            real_t errG = 0;
            for (index_t k = 0; k != d * d; ++k)
            {
                MatD Jp = J, Jm = J;
                Jp(k / d, k % d) += h;
                Jm(k / d, k % d) -= h;
                const real_t fp = garanzhaEnergyDerivatives<d, real_t>(Jp, orientation, eps, form, gaussNewton, nullptr, nullptr);
                const real_t fm = garanzhaEnergyDerivatives<d, real_t>(Jm, orientation, eps, form, gaussNewton, nullptr, nullptr);
                errG = maxError(errG, math::abs((fp - fm) / (2 * h) - g(k)));
                VecJ gp, gm;
                garanzhaEnergyDerivatives<d, real_t>(Jp, orientation, eps, form, gaussNewton, &gp, nullptr);
                garanzhaEnergyDerivatives<d, real_t>(Jm, orientation, eps, form, gaussNewton, &gm, nullptr);
                Hfd.col(k) = (gp - gm) / (2 * h);
                VecJ rp, rm;
                for (index_t s = 0; s != 2; ++s)
                {
                    const MatD& M = s ? Jm : Jp;
                    MatD cof;
                    cofactorMatrix(M, cof);
                    const real_t D = orientation * M.row(0).dot(cof.row(0));
                    const MatD G = inverseDirichlet ? cof : M;
                    VecJ r;
                    for (index_t a = 0; a != d; ++a)
                        for (index_t b = 0; b != d; ++b)
                            r(a * d + b) = G(a, b);
                    (s ? rm : rp) = r / math::pow(garanzhaChi(D, eps), p / 2);
                }
                Jr.col(k) = (rp - rm) / (2 * h);
            }
            CHECK(errG <= 1e-7 * (1 + g.cwiseAbs().maxCoeff()));
            CHECK(maxError(0, (Hfd - H).cwiseAbs().maxCoeff())
                  <= 1e-6 * (1 + H.cwiseAbs().maxCoeff()));
            const MatJ GN = 2 * Jr.transpose() * Jr;
            CHECK(maxError(0, (GN - HGN).cwiseAbs().maxCoeff())
                  <= 1e-7 * (1 + HGN.cwiseAbs().maxCoeff()));
            gsEigen::SelfAdjointEigenSolver<MatJ> es(HP);
            CHECK(es.eigenvalues().minCoeff()
                  >= -1e-12 * (1 + es.eigenvalues().cwiseAbs().maxCoeff()));
        }
    }
}

/// Element-wise assembly (gsRegularizedPDEAssembler) against the expressions
template<short_t d>
void checkElementAssembly(const gsMultiPatch<real_t>& in, const real_t perturbation,
                          const real_t orientation)
{
    Setup s(in, true);
    std::srand(11);
    gsVector<real_t> x = s.freeVector();
    x += perturbation * gsVector<real_t>::Random(x.size());
    const real_t eps = 1e-2;
    gsRegularizedPDEAssembler<d, real_t> fem(s.mp, s.mapper, s.A.options(), orientation);
    fem.setFreeVector(x);

    const real_t e0 = s.energy(x, eps, orientation);
    CHECK(math::abs(fem.energy(eps) - e0) <= 1e-12 * math::abs(e0));

    gsVector<real_t> g;
    gsSparseMatrix<real_t> Hl;
    fem.gradientHessian(eps, g, Hl);
    CHECK(relDiff(g, s.vector(x, [&](space u, geometryMap G)
                  { return garanzhaWinslowGrad(u, G, eps, orientation); })) <= 1e-12);
    const gsMatrix<real_t> H0 = s.matrix(x, [&](space u, geometryMap G)
                                         { return garanzhaWinslowHess(u, G, eps, orientation); });
    const gsMatrix<real_t> H1(Hl);
    gsMatrix<real_t> Hfull = H1 + gsMatrix<real_t>(H1.transpose());
    Hfull.diagonal() = H1.diagonal();
    CHECK(maxError(0, (Hfull - H0).norm()) <= 1e-10 * H0.norm());
    // exactly zero by construction: only the lower triangle is stored
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
        CHECK(maxError(0, (gsMatrix<real_t>(Jh) - J0).norm()) <= 1e-11 * J0.norm());
    }

    gsVector<real_t> r;
    fem.h1Residual(eps, r);
    CHECK(relDiff(r, s.vector(x, [&](space u, geometryMap G)
                  { return regularizedH1Res(u, G, eps, orientation); })) <= 1e-12);
    gsSparseMatrix<real_t> J;
    fem.h1Jacobian(eps, false, J);
    const gsMatrix<real_t> J0 = s.matrix(x, [&](space u, geometryMap G)
                                         { return regularizedH1Jac(u, G, eps, orientation, false); });
    CHECK(maxError(0, (gsMatrix<real_t>(J) - J0).norm()) <= 1e-12 * J0.norm());

    // H1 Jacobians of both forms against central differences of the residual
    for (index_t q = 1; q <= 2; ++q)
    {
        fem.setH1Power(q);
        fem.setFreeVector(x);
        gsSparseMatrix<real_t> Jq;
        fem.h1Jacobian(eps, false, Jq);
        const gsVector<real_t> v = gsVector<real_t>::Random(x.size());
        const gsVector<real_t> an = Jq * v;
        real_t best = std::numeric_limits<real_t>::infinity();
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

/// H2 Jacobian (mode 0, exact) against central differences of the H2
/// residual: independent of the helpers shared by the assembler and the
/// expressions. The warped map has a full metric g, so the mixed second
/// derivatives (off-diagonal entries of adj(g)) contribute.
template<short_t d>
void checkH2Jacobian(const gsMultiPatch<real_t>& in)
{
    Setup s(in, false);
    std::srand(23);
    gsVector<real_t> x = s.freeVector();
    x += 0.02 * gsVector<real_t>::Random(x.size());
    gsRegularizedPDEAssembler<d, real_t> h2(s.mp, s.mapper,
                                            gsExprAssembler<real_t>::defaultOptions(), 1);
    h2.setFreeVector(x);
    gsSparseMatrix<real_t> J;
    h2.h2Jacobian(0, J);
    const gsVector<real_t> v = gsVector<real_t>::Random(x.size());
    const gsVector<real_t> an = J * v;
    real_t best = std::numeric_limits<real_t>::infinity();
    for (real_t h : {1e-4, 1e-5, 1e-6})
    {
        gsVector<real_t> rp, rm;
        h2.setFreeVector(x + h * v);
        h2.h2Residual(rp);
        h2.setFreeVector(x - h * v);
        h2.h2Residual(rm);
        best = std::min(best, relDiff((rp - rm) / (2 * h), an));
    }
    CHECK(best <= 1e-6);
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
            if (a != act)
                continue; // node on an element boundary
            basis.evalAllDers_into(nodes.col(k), 2, ref);
            errV = maxError(errV, (ref[0] - v.col(k)).cwiseAbs().maxCoeff());
            errD = maxError(errD, (ref[1] - d1.col(k)).cwiseAbs().maxCoeff()
                                  / (1 + ref[1].cwiseAbs().maxCoeff()));
            errD2 = maxError(errD2, (ref[2] - d2.col(k)).cwiseAbs().maxCoeff()
                                    / (1 + ref[2].cwiseAbs().maxCoeff()));
        }
    }
    CHECK_EQUAL(fem.numActiveElements(), tabulated);
    CHECK(errV <= 1e-13);
    CHECK(errD <= 1e-13);
    CHECK(errD2 <= 1e-12);
}

/// H1 Jacobian of the warped planar patch (fixed boundary) at its control
/// points, with the node and the component of every unknown
struct H1System
{
    gsSparseMatrix<real_t> A;
    std::vector<index_t> nodes, components;

    explicit H1System(const index_t numElements)
    {
        Setup s(warpedPatch<2>(3, numElements), true);
        gsRegularizedPDEAssembler<2, real_t> fem(s.mp, s.mapper, s.A.options(), 1);
        fem.setFreeVector(s.freeVector());
        fem.h1Jacobian(1e-2, false, A);
        nodes = gsBarrierCore<2, real_t>::nodeOfDofs(s.mp, s.mapper);
        components = gsBarrierCore<2, real_t>::componentOfDofs(s.mp, s.mapper);
    }
};

/// gsBarrierPatch result for a single patch (fixed boundary) with the given
/// PDE options (name without the prefix "PDE.", value)
gsMultiPatch<real_t> parameterize(const gsMultiPatch<real_t>& mp,
                                  const std::vector<std::pair<std::string, real_t> >& settings)
{
    gsBarrierPatch<2, real_t> opt(mp, true);
    opt.options().setInt("Verbose", 0);
    for (const auto& s : settings)
    {
        const std::string name = "PDE." + s.first;
        if (name == "PDE.SkipH2" || name == "PDE.EnergyFinal")
            opt.options().setSwitch(name, s.second != 0);
        else if (name == "PDE.H1Safeguard" || name == "PDE.EnergyTol")
            opt.options().setReal(name, s.second);
        else
            opt.options().setInt(name, static_cast<index_t>(s.second));
    }
    opt.compute();
    return opt.result();
}

} // namespace

SUITE(gsBarrierCore_test)
{
    // Gradient and Hessians (Gauss-Newton, exact, projected) of the energy
    // integrands w.r.t. J against central differences, for unfolded and
    // partly folded J and both orientations
    TEST(energyIntegrandDerivatives)
    {
        checkEnergyIntegrand<2>(1, false);
        checkEnergyIntegrand<2>(-1, true);
        checkEnergyIntegrand<3>(1, false);
        checkEnergyIntegrand<3>(1, true);
    }

    // Element-wise assembly against the gsExprAssembler expressions of the
    // same quantities; perturbation 0.3: folded maps
    TEST(elementAssemblyVsExpressions)
    {
        const gsMultiPatch<real_t> a = warpedPatch<2>(3, 6), b = warpedPatch<3>(2, 3);
        checkElementAssembly<2>(a, 0.01, 1);
        checkElementAssembly<2>(a, 0.3, -1);
        checkElementAssembly<3>(b, 0.02, 1);
        checkElementAssembly<3>(b, 0.3, 1);
    }

    // Exact H2 Jacobian against central differences of the H2 residual
    TEST(h2JacobianVsFiniteDifferences)
    {
        checkH2Jacobian<2>(warpedPatch<2>(3, 5));
        checkH2Jacobian<3>(warpedPatch<3>(2, 3));
    }

    // Assembled energy derivatives: the general path (any Hessian type)
    // agrees with the Gauss-Newton fast path; gradient and exact Hessian
    // of both energy forms against central differences of the energy
    TEST(assembledEnergyDerivatives)
    {
        typedef gsBarrierCore<2, real_t> Core;
        Setup s(warpedPatch<2>(3, 5), true);
        std::srand(19);
        gsVector<real_t> x = s.freeVector();
        x += 0.01 * gsVector<real_t>::Random(x.size());
        gsRegularizedPDEAssembler<2, real_t> fem(s.mp, s.mapper, s.A.options(), 1);
        fem.setFreeVector(x);
        const real_t eps = 1e-2;
        gsVector<real_t> g0, g1;
        gsSparseMatrix<real_t> H0, H1;
        fem.gradientHessian(eps, g0, H0);
        fem.energyGradientHessian(eps, Core::EnergyHessian::GaussNewton, g1, H1);
        CHECK(relDiff(g1, g0) <= 1e-12);
        CHECK(maxError(0, gsMatrix<real_t>(H1 - H0).norm()) <= 1e-12 * gsMatrix<real_t>(H0).norm());
        for (index_t form = 0; form != 2; ++form)
        {
            fem.setEnergyForm(form);
            fem.setFreeVector(x);
            gsVector<real_t> g;
            gsSparseMatrix<real_t> Hl, tmp;
            fem.energyGradientHessian(eps, Core::EnergyHessian::Exact, g, Hl);
            gsMatrix<real_t> Hfull(Hl);
            Hfull += gsMatrix<real_t>(Hl.transpose());
            Hfull.diagonal() /= 2;
            const gsVector<real_t> v = gsVector<real_t>::Random(x.size());
            const gsVector<real_t> Hv = Hfull * v;
            real_t bestG = std::numeric_limits<real_t>::infinity(), bestH = bestG;
            for (real_t h : {1e-4, 1e-5, 1e-6})
            {
                gsVector<real_t> gp, gm;
                fem.setFreeVector(x + h * v);
                const real_t ep = fem.energy(eps);
                fem.energyGradientHessian(eps, Core::EnergyHessian::GaussNewton, gp, tmp);
                fem.setFreeVector(x - h * v);
                const real_t em = fem.energy(eps);
                fem.energyGradientHessian(eps, Core::EnergyHessian::GaussNewton, gm, tmp);
                bestG = std::min(bestG, maxError(0, math::abs((ep - em) / (2 * h) - g.dot(v))
                                                    / math::abs(g.dot(v))));
                bestH = std::min(bestH, relDiff((gp - gm) / (2 * h), Hv));
            }
            CHECK(bestG <= 1e-7);
            CHECK(bestH <= 1e-6);
        }
    }

    // Basis values and derivatives from the univariate tables against
    // evalAllDers: non-uniform knots, a C0 knot, Gauss and Gauss-Lobatto
    // rules, extra nodes on some elements, 2D and 3D
    TEST(univariateTables)
    {
        gsOptionList lobatto = gsExprAssembler<real_t>::defaultOptions(), gauss = lobatto;
        lobatto.setInt("quRule", gsQuadrature::GaussLobatto);
        lobatto.setInt("quB", 2);
        gsMultiPatch<real_t> a = warpedPatch<2>(3, 4);
        gsTensorBSpline<2, real_t>& g = static_cast<gsTensorBSpline<2, real_t>&>(a.patch(0));
        g.insertKnot(0.137, 0);
        g.insertKnot(0.61, 1, 2);
        g.insertKnot(0.42, 0, 3);
        checkTables<2>(a, lobatto);
        checkTables<2>(a, gauss);
        std::vector<index_t> extra(100, 0);
        for (size_t i = 0; i < extra.size(); i += 3)
            extra[i] = 2;
        checkTables<2>(a, lobatto, extra);
        const gsMultiPatch<real_t> b = warpedPatch<3>(2, 3);
        checkTables<3>(b, lobatto);
        checkTables<3>(b, gauss);
    }

    // Maps with harmonic inverse (conformal maps) solve the H2 equations and
    // the weak Laplace form of the H1 equations exactly: residuals at
    // round-off level
    TEST(exactSolutions)
    {
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

    // NodalAMG on an H1 Jacobian (1250 unknowns, several levels): GMRES
    // converges in few iterations; the preconditioner is equivariant under
    // a renumbering of the unknowns and under rotations of the coordinates
    TEST(nodalAMG)
    {
        const short_t d = 2;
        const H1System sys(24);
        const gsSparseMatrix<real_t>& A = sys.A;
        const std::vector<index_t>& nodes = sys.nodes;
        const std::vector<index_t>& comps = sys.components;
        const index_t n = A.rows();
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
        // about 20 iterations; a broken hierarchy needs (many) more
        CHECK(gmres.iterations() < 60);
        CHECK((A * y - b).norm() <= 1e-8 * b.norm());

        gsMatrix<real_t> z;
        amg->apply(b, z);

        // renumbering of the unknowns
        std::vector<index_t> perm(n);
        for (index_t i = 0; i != n; ++i)
            perm[i] = i;
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
                tr.emplace_back(perm[it.row()], perm[it.col()], it.value());
        gsSparseMatrix<real_t> AP(n, n);
        AP.setFromTriplets(tr.begin(), tr.end());
        preAApp::NodalAMG<real_t> amgP;
        amgP.setNodes(nodesP, compsP);
        amgP.compute(AP);
        gsMatrix<real_t> zP;
        amgP.apply(bP, zP);
        real_t err = 0;
        for (index_t i = 0; i != n; ++i)
            err = maxError(err, math::abs(zP(perm[i], 0) - z(i, 0)));
        CHECK(err <= 1e-12 * z.cwiseAbs().maxCoeff());

        // rotation of the physical coordinates: A' = A Q with
        // Q = diag(R, ..., R), then M(A') = Q^T M(A)
        const real_t angle = 0.7;
        gsMatrix<real_t> R(2, 2);
        R << math::cos(angle), -math::sin(angle), math::sin(angle), math::cos(angle);
        std::vector<std::vector<index_t> > dofOf(numNodes, std::vector<index_t>(d));
        for (index_t i = 0; i != n; ++i)
            dofOf[nodes[i]][comps[i]] = i;
        tr.clear();
        for (index_t I = 0; I != numNodes; ++I)
            for (index_t r = 0; r != d; ++r)
                for (index_t k = 0; k != d; ++k)
                    tr.emplace_back(dofOf[I][r], dofOf[I][k], R(r, k));
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

    // BlockGaussSeidel (blocks = coordinate components) against a dense
    // solve with the block lower triangle; with AMG blocks below the coarse
    // size, the AMG solves the blocks directly, so both variants agree
    TEST(blockGaussSeidel)
    {
        const H1System sys(8);  // 2 x 81 unknowns
        const gsSparseMatrix<real_t>& A = sys.A;
        const index_t n = A.rows();
        gsMatrix<real_t> lower(A);
        for (index_t i = 0; i != n; ++i)
            for (index_t j = 0; j != n; ++j)
                if (sys.components[j] > sys.components[i])
                    lower(i, j) = 0;
        std::srand(29);
        const gsVector<real_t> b = gsVector<real_t>::Random(n);
        const gsVector<real_t> expected = lower.partialPivLu().solve(b);
        for (const bool amgBlocks : {false, true})
        {
            preAApp::BlockGaussSeidel<real_t> bgs(amgBlocks);
            bgs.setBlocks(sys.components);
            bgs.compute(A);
            CHECK(bgs.info() == gsEigen::Success);
            gsMatrix<real_t> z;
            bgs.apply(b, z);
            // round-off of the direct solves (condition of the blocks ~1e3)
            CHECK(relDiff(z.col(0), expected) <= 1e-10);
        }
        // the blocks must be set before compute()
        preAApp::BlockGaussSeidel<real_t> unset;
        CHECK_THROW(unset.compute(A), std::runtime_error);
    }

    // Anderson acceleration with the Jacobian as preconditioner on the linear
    // problem A x = b: every iterative linear solver converges to the
    // solution of the direct one; with one GMRES iteration per solve, GMRES
    // falls back to sparse LU
    TEST(andersonAccelerationSolvers)
    {
        const H1System sys(16);
        const gsSparseMatrix<real_t>& A = sys.A;
        std::srand(31);
        const gsVector<real_t> b = gsVector<real_t>::Random(A.rows());
        const gsVector<real_t> x0 = gsVector<real_t>::Zero(A.rows());
        preAApp::AndersonAcceleration<real_t>::Residual_t residual =
            [&](const gsVector<real_t>& x) { return gsVector<real_t>(A * x - b); };
        preAApp::AndersonAcceleration<real_t>::Jacobian_t jacobian =
            [&](const gsVector<real_t>&) { return A; };
        gsSparseSolver<real_t>::LU lu(A);
        const gsVector<real_t> exact = lu.solve(b);
        for (index_t solver = preAAParam<real_t>::SparseLU;
             solver <= preAAParam<real_t>::BlockGaussSeidelAMG; ++solver)
        {
            for (const index_t maxInner : {100, 1})
            {
                preAAParam<real_t> param;
                param.usePreconditioning = true;
                param.epsilon = 1e-10;
                param.linearSolver = solver;
                param.blocks = sys.components;
                param.nodes = sys.nodes;
                param.innerMaxIterations = maxInner;
                preAApp::AndersonAcceleration<real_t> aa(param);
                aa.disableIterationInfoPrinting();
                const gsVector<real_t> x = aa.compute(x0, residual, jacobian);
                CHECK(aa.converged());
                // residual 1e-10, condition of A ~1e4
                CHECK(relDiff(x, exact) <= 1e-5);
                if (solver != preAAParam<real_t>::SparseLU)
                {
                    CHECK(aa.innerIterations() > 0);
                    if (1 == maxInner)
                        CHECK(aa.innerFallbacks() > 0);
                }
            }
        }
    }

    // Anderson acceleration: a step filter that rejects every step stops the
    // iteration at the initial guess; with the divergence safeguard, a
    // residual that is NaN away from the initial guess makes the iteration
    // return there and halve the step until it stops (before, a NaN
    // residual ended the loop with that residual)
    TEST(andersonAccelerationSafeguards)
    {
        const index_t n = 5;
        const gsVector<real_t> x0 = gsVector<real_t>::Ones(n), target = gsVector<real_t>::Zero(n);
        gsSparseMatrix<real_t> identity(n, n);
        identity.setIdentity();
        preAApp::AndersonAcceleration<real_t>::Jacobian_t jacobian =
            [&](const gsVector<real_t>&) { return identity; };
        preAAParam<real_t> param;
        param.usePreconditioning = true;
        {
            preAApp::AndersonAcceleration<real_t> aa(param);
            aa.disableIterationInfoPrinting();
            aa.setStepFilter([](const gsVector<real_t>&, gsVector<real_t>&) { return false; });
            const gsVector<real_t> x = aa.compute(
                x0, [&](const gsVector<real_t>& u) { return gsVector<real_t>(u - target); },
                jacobian);
            CHECK(aa.stalled());
            CHECK(!aa.converged());
            CHECK((x - x0).norm() == 0);
        }
        {
            param.safeguardFactor = 100;
            preAApp::AndersonAcceleration<real_t> aa(param);
            aa.disableIterationInfoPrinting();
            const gsVector<real_t> x = aa.compute(
                x0,
                [&](const gsVector<real_t>& u)
                {
                    return (u - x0).norm() == 0
                           ? gsVector<real_t>(u - target)
                           : gsVector<real_t>::Constant(n, std::numeric_limits<real_t>::quiet_NaN());
                },
                jacobian);
            CHECK(aa.stalled());
            CHECK(math::isfinite(aa.residualNorm()));
            CHECK(aa.rejections() >= 10);  // damping 1/2, ..., 1/1024
            CHECK((x - x0).norm() == 0);
        }
    }

    // couplesPatches: true for a mapper with matched interfaces, false for
    // one without
    TEST(couplesPatches)
    {
        gsMultiPatch<real_t> full = warpedPatch<2>(3, 4), mp;
        gsTensorBSpline<2, real_t> left, right;
        static_cast<const gsTensorBSpline<2, real_t>&>(full.patch(0)).splitAt(0, 0.5, left, right);
        mp.addPatch(left);
        mp.addPatch(right);
        mp.computeTopology();
        const Setup coupled(mp, false);
        bool result = gsBarrierCore<2, real_t>::couplesPatches(mp, coupled.mapper);
        CHECK(result);
        gsDofMapper uncoupled = createMapper(gsMultiBasis<real_t>(mp), 2, false);
        uncoupled.finalize();
        result = gsBarrierCore<2, real_t>::couplesPatches(mp, uncoupled);
        CHECK(!result);
    }

    // Free interfaces: one patch split at u = 1/4; the right part is kept,
    // rotated by 90 or 180 degrees, or reversed in one direction (opposite
    // orientation) in its parameter domain; all parameter domains are
    // [0,1]^2. The result must be bijective, with the orientation of each
    // patch, and must not depend on these choices (to the tolerance of the
    // energy minimization, PDE.EnergyTol = 1e-13)
    TEST(freeInterfaces)
    {
        gsMultiPatch<real_t> full = warpedPatch<2>(3, 8);
        const gsTensorBSpline<2, real_t>& g = static_cast<const gsTensorBSpline<2, real_t>&>(full.patch(0));
        gsTensorBSpline<2, real_t> left, right;
        g.splitAt(0, 0.25, left, right);
        left = onUnitSquare(left);

        const gsMatrix<real_t> pts = unitGrid(11);
        gsMatrix<real_t> ref;
        const real_t size = modelSize(full);
        for (index_t variant = 0; variant != 4; ++variant)
        {
            gsTensorBSpline<2, real_t> r = onUnitSquare(right);
            if (variant == 1) { r.reverse(0); r.swapDirections(0, 1); }
            if (variant == 2) { r.reverse(0); r.reverse(1); }
            if (variant == 3) { r.reverse(0); }
            gsMultiPatch<real_t> mp;
            mp.addPatch(left);
            mp.addPatch(r);
            mp.computeTopology();
            collapseInterior(mp, 0.5);
            if (variant == 1)
            {
                // the interface has the same parametric length on both sides
                const gsMultiPatch<real_t> scaled = gsBarrierCore<2, real_t>::chordScaled(mp);
                const boundaryInterface& bi = *scaled.topology().iBegin();
                const index_t a = 1 - bi.first().direction();
                const index_t b = bi.dirMap(bi.first())(a);
                const real_t l1 = scaled.basis(bi.first().patch).support()(a, 1);
                const real_t l2 = scaled.basis(bi.second().patch).support()(b, 1);
                // round-off of the chord lengths
                CHECK(math::abs(l1 - l2) <= 1e-14 * l1);
            }
            gsBarrierPatch<2, real_t> opt(mp, false);
            opt.options().setInt("Verbose", 0);
            opt.options().setReal("PDE.EnergyTol", 1e-13);
            opt.compute();
            const gsMultiPatch<real_t>& res = opt.result();
            for (size_t p = 0; p != res.nPatches(); ++p)
                CHECK_EQUAL(0, inconsistentJacobians(res.patch(p), 30));
            gsMatrix<real_t> xl;
            res.patch(0).eval_into(pts, xl);
            if (variant == 0)
                ref = xl;
            else
                CHECK((xl - ref).cwiseAbs().maxCoeff() <= 1e-5 * size);
        }
    }

    // Full runs with the linear solvers and options that the other tests do
    // not reach (folded start, so that all stages run): direct and iterative
    // solvers (also by the automatic choice) give the same result to the
    // tolerance of the H1 stage; the experimental H1 form, the adaptive
    // quadrature and the energy as the final stage give bijective results
    TEST(linearSolversAndOptions)
    {
        typedef gsBarrierCore<2, real_t> Core;
        gsMultiPatch<real_t> mp = warpedPatch<2>(3, 8);
        collapseInterior(mp, 1);
        const real_t size = modelSize(mp);
        const gsMultiPatch<real_t> direct =
            parameterize(mp, {{"SkipH2", 1}, {"LinearSolver", Core::LinearSolver::Direct}});
        const gsMatrix<real_t> pts = unitGrid(11);
        gsMatrix<real_t> xd, x;
        direct.patch(0).eval_into(pts, xd);
        CHECK_EQUAL(0, inconsistentJacobians(direct.patch(0), 30));
        const std::vector<std::vector<std::pair<std::string, real_t> > > iterative = {
            {{"SkipH2", 1}, {"LinearSolver", Core::LinearSolver::Iterative}},
            {{"SkipH2", 1}, {"DirectSolverLimit", 1}},
            {{"SkipH2", 1}, {"LinearSolver", Core::LinearSolver::Iterative},
             {"Preconditioner", Core::Preconditioner::BlockGaussSeidelLU}},
            {{"LinearSolver", Core::LinearSolver::Iterative}},  // with H2
        };
        for (const auto& settings : iterative)
        {
            const gsMultiPatch<real_t> res = parameterize(mp, settings);
            CHECK_EQUAL(0, inconsistentJacobians(res.patch(0), 30));
            res.patch(0).eval_into(pts, x);
            // H1 tolerance 1e-5 (residual norm)
            CHECK((x - xd).cwiseAbs().maxCoeff() <= 1e-4 * size);
        }
        const std::vector<std::vector<std::pair<std::string, real_t> > > others = {
            {{"SkipH2", 1}, {"H1Form", Core::H1Form::WeakLaplace}},
            {{"SkipH2", 1}, {"AdaptiveRounds", 2}},
            {{"SkipH2", 1}, {"EnergyFinal", 1}},
            {{"SkipH2", 1}, {"EnergyFinal", 1}, {"EnergyHessian", Core::EnergyHessian::ProjectedNewton}},
        };
        for (const auto& settings : others)
            CHECK_EQUAL(0, inconsistentJacobians(parameterize(mp, settings).patch(0), 30));
    }

    // gsBarrierCore::compute directly (no rescaling): square [-1,1]^2, all
    // interior control points exactly at the origin, J = 0 exactly in the
    // elements without boundary control points; the nodes there must still
    // be monitored
    TEST(degenerateInitialGuess)
    {
        const gsKnotVector<real_t> kv(-1, 1, 10, 4);
        gsTensorBSplineBasis<2, real_t> basis(kv, kv);
        gsMatrix<real_t> C = basis.anchors().transpose();
        const index_t n = basis.size(0);
        for (index_t j = 1; j + 1 < n; ++j)
            for (index_t i = 1; i + 1 < n; ++i)
                C.row(j * n + i).setZero();
        const gsMultiPatch<real_t> mp(gsTensorBSpline<2, real_t>(basis, C));
        gsMultiBasis<real_t> mb(mp);
        gsDofMapper mapper = createMapper(mb, 2, false);
        for (index_t c = 0; c != 2; ++c)
            mapper.markBoundary(0, basis.allBoundary(), c);
        mapper.finalize();
        gsOptionList opt = gsBarrierCore<2, real_t>::defaultOptions();
        opt.setSwitch("PDE.SkipH2", true);
        opt.setInt("Verbose", 0);
        const gsMultiPatch<real_t> res = gsBarrierCore<2, real_t>::compute(mp, mapper, opt);
        CHECK_EQUAL(0, inconsistentJacobians(res.patch(0), 100));
    }
}
