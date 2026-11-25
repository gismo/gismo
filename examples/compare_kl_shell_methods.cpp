/** @file compare_kl_shell_methods.cpp

    @brief Head-to-head benchmark, pointwise value-difference analyzer, and
           unified ParaView (.pvd / .vtp) exporter comparing:
             (A) Ours: Conforming AS-G1 (gsAsG1KLShellAssembler + gsAsG1ShellHierarchy)
                       with MG-PCG (line_ring1) and SimplicialLDLT (direct),
             (B) gsKLShell Weak Penalty C0/C1 (gsThinShellAssembler + addWeakC0/addWeakC1),
             (C) gsKLShell + gsC1SurfSpline (Farahat et al. 2023 mapped AS-G1 + direct),
             (D) gsKLShell + gsApproxC1Spline (Verhelst et al. 2024 approx C1 + direct).

           All methods run on the EXACT same gsMultiPatch geometry, same (p, k)
           refinement level, matched Gauss quadrature order, and identical
           source functions and boundary conditions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsModeling/gsAsG1ShellSolver.hpp>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
#endif

#ifdef gsUnstructuredSplines_ENABLED
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace gismo;

namespace {

void mkdirRecursive(const std::string& path)
{
    std::string cur;
    for (size_t i = 0; i < path.size(); ++i)
    {
        cur.push_back(path[i]);
        if (path[i] == '/' && cur.size() > 1)
            ::mkdir(cur.c_str(), 0755);
    }
    if (!cur.empty())
        ::mkdir(cur.c_str(), 0755);
}

void normalizePatchSupports(gsMultiPatch<real_t>& mp)
{
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        if (auto* bsp = dynamic_cast<gsTensorBSpline<2, real_t>*>(&mp.patch(p)))
        {
            bsp->basis().knots(0).affineTransformTo(0.0, 1.0);
            bsp->basis().knots(1).affineTransformTo(0.0, 1.0);
        }
        else if (auto* nurbs = dynamic_cast<gsTensorNurbs<2, real_t>*>(&mp.patch(p)))
        {
            nurbs->basis().source().knots(0).affineTransformTo(0.0, 1.0);
            nurbs->basis().source().knots(1).affineTransformTo(0.0, 1.0);
        }
    }
}

gsMultiPatch<real_t> loadSurface3D(const std::string& filename, bool split2x2)
{
    gsMultiPatch<real_t> mp;
    gsFileData<real_t> fd(filename);
    if (fd.has<gsMultiPatch<real_t>>())
        mp = *fd.getFirst<gsMultiPatch<real_t>>();
    else if (fd.has<gsGeometry<real_t>>())
    {
        auto geoms = fd.getAll<gsGeometry<real_t>>();
        for (size_t i = 0; i < geoms.size(); ++i)
            mp.addPatch(give(*geoms[i]));
    }
    else
        GISMO_ERROR("No geometry found in file: " << filename);

    if (mp.targetDim() == 2)
        mp.embed(3);

    normalizePatchSupports(mp);
    mp.computeTopology();

    if (split2x2)
    {
        mp = mp.uniformSplit();
        normalizePatchSupports(mp);
        mp.computeTopology();
    }
    return mp;
}

bool isAllPolynomialBSpline(const gsMultiPatch<real_t>& mp)
{
    for (size_t p = 0; p < mp.nPatches(); ++p)
        if (!dynamic_cast<const gsTensorBSpline<2, real_t>*>(&mp.patch(p)))
            return false;
    return true;
}

bool hasZeroInteriorKnots(const gsMultiPatch<real_t>& mp)
{
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        const auto* bsp = dynamic_cast<const gsTensorBSpline<2, real_t>*>(&mp.patch(p));
        if (!bsp) return false;
        if (bsp->basis().component(0).knots().unique().size() != 2 ||
            bsp->basis().component(1).knots().unique().size() != 2)
            return false;
    }
    return true;
}

/// Shared 3D physical load function used identically by gsAsG1KLShellAssembler and gsThinShellAssembler.
class SharedShellForceFunction : public gsFunction<real_t>
{
public:
    typedef std::function<gsVector<real_t, 3>(const gsVector<real_t>&)> ForceFn;

    explicit SharedShellForceFunction(ForceFn fn) : m_fn(std::move(fn)) {}
    GISMO_CLONE_FUNCTION(SharedShellForceFunction)

    short_t domainDim() const override { return 3; }
    short_t targetDim() const override { return 3; }
    const gsFunction<real_t>& piece(const index_t) const override { return *this; }

    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override
    {
        const index_t n = u.cols();
        result.resize(3, n);
        for (index_t k = 0; k < n; ++k)
        {
            const gsVector<real_t, 3> val = m_fn(gsVector<real_t>(u.col(k)));
            result.col(k) = val;
        }
    }

    gsVector<real_t, 3> evalPoint(const gsVector<real_t>& x) const { return m_fn(x); }

private:
    ForceFn m_fn;
};

/// Exact analytical displacement function for the planar plate exact-source comparison.
class ExactPlateDisplacementFunction : public gsFunction<real_t>
{
public:
    explicit ExactPlateDisplacementFunction(bool clamped) : m_clamped(clamped) {}
    GISMO_CLONE_FUNCTION(ExactPlateDisplacementFunction)

    short_t domainDim() const override { return 3; }
    short_t targetDim() const override { return 3; }
    const gsFunction<real_t>& piece(const index_t) const override { return *this; }

    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override
    {
        const real_t pi = EIGEN_PI;
        const index_t n = u.cols();
        result.setZero(3, n);
        for (index_t k = 0; k < n; ++k)
        {
            const real_t sx = std::sin(pi * u(0, k));
            const real_t sy = std::sin(pi * u(1, k));
            result(2, k) = m_clamped ? (sx * sx * sy * sy) : (sx * sy);
        }
    }

    void evalJet(index_t comp, const gsVector<real_t>& xP,
                 real_t& val, gsVector<real_t>& g, gsMatrix<real_t>& H) const
    {
        const real_t pi = EIGEN_PI;
        g.setZero(3);
        H.setZero(3, 3);
        if (comp != 2)
        {
            val = 0.0;
            return;
        }
        const real_t x = xP(0), y = xP(1);
        const real_t sx = std::sin(pi * x), cx = std::cos(pi * x);
        const real_t sy = std::sin(pi * y), cy = std::cos(pi * y);
        if (!m_clamped)
        {
            val = sx * sy;
            g(0) = pi * cx * sy;
            g(1) = pi * sx * cy;
            H(0, 0) = -pi * pi * sx * sy;
            H(1, 1) = -pi * pi * sx * sy;
            H(0, 1) = H(1, 0) = pi * pi * cx * cy;
        }
        else
        {
            const real_t s2x = std::sin(2.0 * pi * x), c2x = std::cos(2.0 * pi * x);
            const real_t s2y = std::sin(2.0 * pi * y), c2y = std::cos(2.0 * pi * y);
            const real_t gx = sx * sx, gy = sy * sy;
            const real_t dgx = pi * s2x, dgy = pi * s2y;
            const real_t ddgx = 2.0 * pi * pi * c2x, ddgy = 2.0 * pi * pi * c2y;
            val = gx * gy;
            g(0) = dgx * gy;
            g(1) = gx * dgy;
            H(0, 0) = ddgx * gy;
            H(1, 1) = gx * ddgy;
            H(0, 1) = H(1, 0) = dgx * dgy;
        }
    }

private:
    bool m_clamped;
};

/// Fast O(nPts * (p+1)^2) evaluator for gsMappedBasis solutions (values, 1st & 2nd derivs).
/// Avoids gsMappedBasis::deriv_into's O(nPts * nAct * nGlobal) per-function sparse assembly.
class FastMappedPatch : public gsFunction<real_t>
{
public:
    FastMappedPatch() = default;
    FastMappedPatch(const gsBasis<real_t>& basis, const gsMatrix<real_t>& localCoefs)
        : m_basis(basis.clone().release()), m_coefs(localCoefs) {}
    FastMappedPatch(const FastMappedPatch& other)
        : m_basis(other.m_basis ? other.m_basis->clone().release() : nullptr),
          m_coefs(other.m_coefs) {}
    FastMappedPatch& operator=(const FastMappedPatch& other)
    {
        if (this != &other)
        {
            m_basis.reset(other.m_basis ? other.m_basis->clone().release() : nullptr);
            m_coefs = other.m_coefs;
        }
        return *this;
    }
    FastMappedPatch(FastMappedPatch&&) = default;
    FastMappedPatch& operator=(FastMappedPatch&&) = default;

    GISMO_CLONE_FUNCTION(FastMappedPatch)

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return static_cast<short_t>(m_coefs.cols()); }

    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override
    {
        gsMatrix<index_t> act;
        gsMatrix<real_t> B;
        m_basis->active_into(u, act);
        m_basis->eval_into(u, B);
        const index_t nPts = u.cols(), nAct = act.rows(), dim = m_coefs.cols();
        result.setZero(dim, nPts);
        for (index_t q = 0; q < nPts; ++q)
            for (index_t a = 0; a < nAct; ++a)
            {
                const index_t idx = act(a, q);
                const real_t val = B(a, q);
                for (index_t c = 0; c < dim; ++c)
                    result(c, q) += val * m_coefs(idx, c);
            }
    }

    void deriv_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override
    {
        gsMatrix<index_t> act;
        gsMatrix<real_t> dB;
        m_basis->active_into(u, act);
        m_basis->deriv_into(u, dB);
        const index_t nPts = u.cols(), nAct = act.rows(), dim = m_coefs.cols();
        result.setZero(2 * dim, nPts);
        for (index_t q = 0; q < nPts; ++q)
            for (index_t a = 0; a < nAct; ++a)
            {
                const index_t idx = act(a, q);
                const real_t du = dB(2 * a + 0, q);
                const real_t dv = dB(2 * a + 1, q);
                for (index_t c = 0; c < dim; ++c)
                {
                    result(2 * c + 0, q) += du * m_coefs(idx, c);
                    result(2 * c + 1, q) += dv * m_coefs(idx, c);
                }
            }
    }

    void deriv2_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override
    {
        gsMatrix<index_t> act;
        gsMatrix<real_t> ddB;
        m_basis->active_into(u, act);
        m_basis->deriv2_into(u, ddB);
        const index_t nPts = u.cols(), nAct = act.rows(), dim = m_coefs.cols();
        result.setZero(3 * dim, nPts);
        for (index_t q = 0; q < nPts; ++q)
            for (index_t a = 0; a < nAct; ++a)
            {
                const index_t idx = act(a, q);
                const real_t duu = ddB(3 * a + 0, q);
                const real_t dvv = ddB(3 * a + 1, q);
                const real_t duv = ddB(3 * a + 2, q);
                for (index_t c = 0; c < dim; ++c)
                {
                    result(3 * c + 0, q) += duu * m_coefs(idx, c);
                    result(3 * c + 1, q) += dvv * m_coefs(idx, c);
                    result(3 * c + 2, q) += duv * m_coefs(idx, c);
                }
            }
    }

private:
    std::unique_ptr<gsBasis<real_t>> m_basis;
    gsMatrix<real_t> m_coefs;
};

class FastMappedDisplacement : public gsFunctionSet<real_t>
{
public:
    FastMappedDisplacement() = default;

    void init(const gsMappedBasis<2, real_t>& bb2, const gsMatrix<real_t>& globalCoefs)
    {
        gsMatrix<real_t> localCoefs;
        bb2.global_coef_to_local_coef(globalCoefs, localCoefs);
        m_pieces.clear();
        m_pieces.reserve(bb2.nPatches());
        index_t offset = 0;
        for (size_t p = 0; p < bb2.nPatches(); ++p)
        {
            const index_t sz = bb2.localSize(static_cast<index_t>(p));
            m_pieces.emplace_back(bb2.getBase(static_cast<int>(p)),
                                  localCoefs.middleRows(offset, sz));
            offset += sz;
        }
    }

    GISMO_CLONE_FUNCTION(FastMappedDisplacement)

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return m_pieces.empty() ? 3 : m_pieces[0].targetDim(); }
    index_t nPieces() const override { return static_cast<index_t>(m_pieces.size()); }
    const gsFunction<real_t>& piece(const index_t k) const override { return m_pieces[static_cast<size_t>(k)]; }

    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override { m_pieces[0].eval_into(u, result); }
    void deriv_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override { m_pieces[0].deriv_into(u, result); }
    void deriv2_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const override { m_pieces[0].deriv2_into(u, result); }

private:
    std::vector<FastMappedPatch> m_pieces;
};

struct LogicalBCEntry
{
    patchSide ps;
    short_t comp;
    bool clamped; // false => Values / Dirichlet only; true => Values + Normal Slope (Clamped)
};

struct MethodResult
{
    std::string name;
    bool ok = false;
    std::string note;

    index_t nFree = 0;
    index_t nnz = 0;
    index_t iterations = 0;
    real_t relRes = 0.0;
    double basisTimeSec = 0.0;
    double asmTimeSec = 0.0;
    double solveTimeSec = 0.0;
    double totalTimeSec = 0.0;

    real_t strainEnergy = 0.0; // 0.5 * u^T K u
    real_t dispL2Norm = 0.0;
    real_t dispMaxMag = 0.0;
    real_t probeUz = 0.0;

    real_t c0JumpL2 = 0.0;
    real_t c0JumpMax = 0.0;
    real_t c1JumpL2 = 0.0;
    real_t c1JumpMax = 0.0;
    real_t rotJumpL2 = 0.0;
    real_t rotJumpMax = 0.0;

    // Concrete displacement representation (exactly one is active; both inherit from gsFunctionSet<real_t>)
    bool usesMapped = false;
    gsMultiPatch<real_t> dispMP;
    std::unique_ptr<FastMappedDisplacement> dispMapped;

    const gsFunctionSet<real_t>& dispFunc() const
    {
        if (usesMapped && dispMapped) return *dispMapped;
        return dispMP;
    }
};

struct PairwiseDiff
{
    std::string nameA;
    std::string nameB;
    real_t lInfDiff = 0.0;
    real_t lInfRel = 0.0;
    real_t dxMax = 0.0, dyMax = 0.0, dzMax = 0.0;
    gsVector<real_t, 3> maxDiffPos = gsVector<real_t, 3>::Zero();

    real_t l2Diff = 0.0, l2Rel = 0.0;
    real_t h1Diff = 0.0, h1Rel = 0.0;
    real_t h2Diff = 0.0, h2Rel = 0.0;
    real_t enDiff = 0.0, enRel = 0.0;
};

struct InterfaceJumpSample
{
    gsVector<real_t, 3> pos;
    real_t c0Jump = 0.0;
    real_t c1Jump = 0.0;
    real_t rotJump = 0.0;
};

/// Evaluate surface L2/H1/broken-H2/energy norms and pointwise differences between two
/// gsFunctionSet<real_t> displacement fields defined on the exact same gsMultiPatch mp.
PairwiseDiff computePairwiseDiff(const gsMultiPatch<real_t>& mp,
                                 const gsMultiBasis<real_t>& dbasis,
                                 const gsKLShellMaterial<real_t>& mat,
                                 const gsFunctionSet<real_t>& fA,
                                 const gsFunctionSet<real_t>& fB,
                                 const std::string& nameA,
                                 const std::string& nameB)
{
    PairwiseDiff out;
    out.nameA = nameA;
    out.nameB = nameB;

    const real_t tMem = mat.thickness;
    const real_t tBend = (mat.thickness * mat.thickness * mat.thickness) / 12.0;

    real_t l2E2 = 0.0, l2N2 = 0.0;
    real_t h1E2 = 0.0, h1N2 = 0.0;
    real_t h2E2 = 0.0, h2N2 = 0.0;
    real_t enE2 = 0.0, enN2 = 0.0;
    real_t maxAMag = 0.0;

    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        const auto& tb = dynamic_cast<const gsTensorBSplineBasis<2, real_t>&>(dbasis[p]);
        const gsGeometry<real_t>& geo = mp.patch(p);
        const gsFunctionSet<real_t>& pA = fA.piece(static_cast<index_t>(p));
        const gsFunctionSet<real_t>& pB = fB.piece(static_cast<index_t>(p));

        const std::vector<real_t> brkU = tb.component(0).knots().breaks();
        const std::vector<real_t> brkV = tb.component(1).knots().breaks();
        const gsGaussRule<real_t> ruleU(tb.degree(0) + 3);
        const gsGaussRule<real_t> ruleV(tb.degree(1) + 3);
        const index_t nQu = ruleU.numNodes(), nQv = ruleV.numNodes();
        const index_t nQ = nQu * nQv;

        gsMatrix<real_t> qu, qv, pts(2, nQ);
        gsVector<real_t> wu, wv, wq(nQ);
        gsMatrix<real_t> Xq, gD1, gD2;
        gsMatrix<real_t> vA, d1A, d2A, vB, d1B, d2B;
        gsMatrix<real_t, 3, 3> C;

        for (size_t ev = 0; ev + 1 < brkV.size(); ++ev)
        {
            ruleV.mapTo(brkV[ev], brkV[ev + 1], qv, wv);
            for (size_t eu = 0; eu + 1 < brkU.size(); ++eu)
            {
                ruleU.mapTo(brkU[eu], brkU[eu + 1], qu, wu);
                for (index_t j = 0; j < nQv; ++j)
                    for (index_t i = 0; i < nQu; ++i)
                    {
                        const index_t q = i + j * nQu;
                        pts(0, q) = qu(0, i);
                        pts(1, q) = qv(0, j);
                        wq(q) = wu(i) * wv(j);
                    }

                geo.eval_into(pts, Xq);
                geo.deriv_into(pts, gD1);
                geo.deriv2_into(pts, gD2);

                pA.eval_into(pts, vA);   pA.deriv_into(pts, d1A);   pA.deriv2_into(pts, d2A);
                pB.eval_into(pts, vB);   pB.deriv_into(pts, d1B);   pB.deriv2_into(pts, d2B);

                for (index_t q = 0; q < nQ; ++q)
                {
                    gsSurfaceKinematics3D<real_t> kin;
                    if (!kin.compute(gD1, gD2, q))
                        continue;
                    mat.voigtMatrix(kin.A11inv, kin.A22inv, kin.A12inv, C);

                    gsVector<real_t, 3> uA, uA1, uA2, uA11, uA22, uA12;
                    gsVector<real_t, 3> uB, uB1, uB2, uB11, uB22, uB12;
                    for (int c = 0; c < 3; ++c)
                    {
                        uA(c)   = vA(c, q);
                        uA1(c)  = d1A(2 * c + 0, q);
                        uA2(c)  = d1A(2 * c + 1, q);
                        uA11(c) = d2A(3 * c + 0, q);
                        uA22(c) = d2A(3 * c + 1, q);
                        uA12(c) = d2A(3 * c + 2, q);

                        uB(c)   = vB(c, q);
                        uB1(c)  = d1B(2 * c + 0, q);
                        uB2(c)  = d1B(2 * c + 1, q);
                        uB11(c) = d2B(3 * c + 0, q);
                        uB22(c) = d2B(3 * c + 1, q);
                        uB12(c) = d2B(3 * c + 2, q);
                    }

                    const gsVector<real_t, 3> e   = uA   - uB;
                    const gsVector<real_t, 3> e1  = uA1  - uB1;
                    const gsVector<real_t, 3> e2  = uA2  - uB2;
                    const gsVector<real_t, 3> e11 = uA11 - uB11;
                    const gsVector<real_t, 3> e22 = uA22 - uB22;
                    const gsVector<real_t, 3> e12 = uA12 - uB12;

                    maxAMag = std::max(maxAMag, uA.norm());
                    const real_t eNorm = e.norm();
                    if (eNorm > out.lInfDiff)
                    {
                        out.lInfDiff = eNorm;
                        out.maxDiffPos << Xq(0, q), Xq(1, q), Xq(2, q);
                    }
                    out.dxMax = std::max(out.dxMax, std::abs(e(0)));
                    out.dyMax = std::max(out.dyMax, std::abs(e(1)));
                    out.dzMax = std::max(out.dzMax, std::abs(e(2)));

                    const real_t wDS = wq(q) * kin.dS;
                    l2E2 += wDS * e.squaredNorm();
                    l2N2 += wDS * uA.squaredNorm();

                    auto surfH1Sq = [&](const gsVector<real_t, 3>& d1, const gsVector<real_t, 3>& d2) -> real_t {
                        return kin.A11inv * d1.squaredNorm()
                             + 2.0 * kin.A12inv * d1.dot(d2)
                             + kin.A22inv * d2.squaredNorm();
                    };
                    h1E2 += wDS * surfH1Sq(e1, e2);
                    h1N2 += wDS * surfH1Sq(uA1, uA2);

                    auto surfH2Sq = [&](const gsVector<real_t, 3>& d1, const gsVector<real_t, 3>& d2,
                                        const gsVector<real_t, 3>& d11, const gsVector<real_t, 3>& d22,
                                        const gsVector<real_t, 3>& d12) -> real_t {
                        const gsVector<real_t, 3> c11 = d11 - kin.G1_11 * d1 - kin.G2_11 * d2;
                        const gsVector<real_t, 3> c22 = d22 - kin.G1_22 * d1 - kin.G2_22 * d2;
                        const gsVector<real_t, 3> c12 = d12 - kin.G1_12 * d1 - kin.G2_12 * d2;
                        const real_t a11 = kin.A11inv, a22 = kin.A22inv, a12 = kin.A12inv;
                        real_t sum = 0.0;
                        for (int comp = 0; comp < 3; ++comp)
                        {
                            const real_t h00 = c11(comp), h11 = c22(comp), h01 = c12(comp);
                            const real_t m00 = a11 * h00 + a12 * h01;
                            const real_t m01 = a11 * h01 + a12 * h11;
                            const real_t m10 = a12 * h00 + a22 * h01;
                            const real_t m11 = a12 * h01 + a22 * h11;
                            sum += m00 * m00 + 2.0 * m01 * m10 + m11 * m11;
                        }
                        return sum;
                    };
                    h2E2 += wDS * surfH2Sq(e1, e2, e11, e22, e12);
                    h2N2 += wDS * surfH2Sq(uA1, uA2, uA11, uA22, uA12);

                    auto shellEnDens = [&](const gsVector<real_t, 3>& d1, const gsVector<real_t, 3>& d2,
                                           const gsVector<real_t, 3>& d11, const gsVector<real_t, 3>& d22,
                                           const gsVector<real_t, 3>& d12) -> real_t {
                        gsVector<real_t, 3> epsV, kapV;
                        epsV(0) = kin.A1.dot(d1);
                        epsV(1) = kin.A2.dot(d2);
                        epsV(2) = kin.A1.dot(d2) + kin.A2.dot(d1);
                        kapV(0) = -(d11 - kin.G1_11 * d1 - kin.G2_11 * d2).dot(kin.A3);
                        kapV(1) = -(d22 - kin.G1_22 * d1 - kin.G2_22 * d2).dot(kin.A3);
                        kapV(2) = -2.0 * (d12 - kin.G1_12 * d1 - kin.G2_12 * d2).dot(kin.A3);
                        return tMem * epsV.dot(C * epsV) + tBend * kapV.dot(C * kapV);
                    };
                    enE2 += wDS * shellEnDens(e1, e2, e11, e22, e12);
                    enN2 += wDS * shellEnDens(uA1, uA2, uA11, uA22, uA12);
                }
            }
        }
    }

    out.lInfRel = (maxAMag > 0.0) ? out.lInfDiff / maxAMag : out.lInfDiff;
    out.l2Diff = std::sqrt(std::max(l2E2, 0.0));
    out.l2Rel  = (l2N2 > 0.0) ? std::sqrt(std::max(l2E2, 0.0) / l2N2) : out.l2Diff;
    out.h1Diff = std::sqrt(std::max(h1E2, 0.0));
    out.h1Rel  = (h1N2 > 0.0) ? std::sqrt(std::max(h1E2, 0.0) / h1N2) : out.h1Diff;
    out.h2Diff = std::sqrt(std::max(h2E2, 0.0));
    out.h2Rel  = (h2N2 > 0.0) ? std::sqrt(std::max(h2E2, 0.0) / h2N2) : out.h2Diff;
    out.enDiff = std::sqrt(std::max(enE2, 0.0));
    out.enRel  = (enN2 > 0.0) ? std::sqrt(std::max(enE2, 0.0) / enN2) : out.enDiff;
    return out;
}

/// Evaluate surface L2 norm, max magnitude, and Koiter strain energy of a displacement gsFunctionSet.
void computeFieldStats(const gsMultiPatch<real_t>& mp,
                       const gsMultiBasis<real_t>& dbasis,
                       const gsKLShellMaterial<real_t>& mat,
                       const gsFunctionSet<real_t>& fU,
                       real_t& l2Norm,
                       real_t& maxMag,
                       real_t& strainEnergy)
{
    const real_t tMem = mat.thickness;
    const real_t tBend = (mat.thickness * mat.thickness * mat.thickness) / 12.0;
    real_t l2Sq = 0.0, enSq = 0.0;
    maxMag = 0.0;

    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        const auto& tb = dynamic_cast<const gsTensorBSplineBasis<2, real_t>&>(dbasis[p]);
        const gsGeometry<real_t>& geo = mp.patch(p);
        const gsFunctionSet<real_t>& pU = fU.piece(static_cast<index_t>(p));

        const std::vector<real_t> brkU = tb.component(0).knots().breaks();
        const std::vector<real_t> brkV = tb.component(1).knots().breaks();
        const gsGaussRule<real_t> ruleU(tb.degree(0) + 2);
        const gsGaussRule<real_t> ruleV(tb.degree(1) + 2);
        const index_t nQu = ruleU.numNodes(), nQv = ruleV.numNodes();
        const index_t nQ = nQu * nQv;

        gsMatrix<real_t> qu, qv, pts(2, nQ);
        gsVector<real_t> wu, wv, wq(nQ);
        gsMatrix<real_t> gD1, gD2, vU, d1U, d2U;
        gsMatrix<real_t, 3, 3> C;

        for (size_t ev = 0; ev + 1 < brkV.size(); ++ev)
        {
            ruleV.mapTo(brkV[ev], brkV[ev + 1], qv, wv);
            for (size_t eu = 0; eu + 1 < brkU.size(); ++eu)
            {
                ruleU.mapTo(brkU[eu], brkU[eu + 1], qu, wu);
                for (index_t j = 0; j < nQv; ++j)
                    for (index_t i = 0; i < nQu; ++i)
                    {
                        const index_t q = i + j * nQu;
                        pts(0, q) = qu(0, i);
                        pts(1, q) = qv(0, j);
                        wq(q) = wu(i) * wv(j);
                    }
                geo.deriv_into(pts, gD1);
                geo.deriv2_into(pts, gD2);
                pU.eval_into(pts, vU);
                pU.deriv_into(pts, d1U);
                pU.deriv2_into(pts, d2U);

                for (index_t q = 0; q < nQ; ++q)
                {
                    gsSurfaceKinematics3D<real_t> kin;
                    if (!kin.compute(gD1, gD2, q)) continue;
                    mat.voigtMatrix(kin.A11inv, kin.A22inv, kin.A12inv, C);

                    gsVector<real_t, 3> u, u1, u2, u11, u22, u12;
                    for (int c = 0; c < 3; ++c)
                    {
                        u(c) = vU(c, q);
                        u1(c) = d1U(2 * c + 0, q);   u2(c) = d1U(2 * c + 1, q);
                        u11(c) = d2U(3 * c + 0, q);  u22(c) = d2U(3 * c + 1, q);  u12(c) = d2U(3 * c + 2, q);
                    }
                    maxMag = std::max(maxMag, u.norm());
                    const real_t wDS = wq(q) * kin.dS;
                    l2Sq += wDS * u.squaredNorm();

                    gsVector<real_t, 3> epsV, kapV;
                    epsV(0) = kin.A1.dot(u1);
                    epsV(1) = kin.A2.dot(u2);
                    epsV(2) = kin.A1.dot(u2) + kin.A2.dot(u1);
                    kapV(0) = -(u11 - kin.G1_11 * u1 - kin.G2_11 * u2).dot(kin.A3);
                    kapV(1) = -(u22 - kin.G1_22 * u1 - kin.G2_22 * u2).dot(kin.A3);
                    kapV(2) = -2.0 * (u12 - kin.G1_12 * u1 - kin.G2_12 * u2).dot(kin.A3);
                    enSq += wDS * (tMem * epsV.dot(C * epsV) + tBend * kapV.dot(C * kapV));
                }
            }
        }
    }
    l2Norm = std::sqrt(std::max(l2Sq, 0.0));
    strainEnergy = 0.5 * enSq;
}

/// Measure C0 displacement jumps and C1 normal-derivative / rotation jumps across all interior interfaces.
std::vector<std::vector<InterfaceJumpSample>> computeInterfaceJumps(
    const gsMultiPatch<real_t>& mp,
    const gsMultiBasis<real_t>& dbasis,
    const gsFunctionSet<real_t>& fU,
    real_t& c0L2, real_t& c0Max,
    real_t& c1L2, real_t& c1Max,
    real_t& rotL2, real_t& rotMax)
{
    c0L2 = c0Max = c1L2 = c1Max = rotL2 = rotMax = 0.0;
    real_t c0Sq = 0.0, c1Sq = 0.0, rotSq = 0.0;
    std::vector<std::vector<InterfaceJumpSample>> ifcCurves;

    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface& ifc = *it;
        const patchSide ps1 = ifc.first();
        const patchSide ps2 = ifc.second();
        const short_t nDir1 = ps1.direction(), tDir1 = 1 - nDir1;
        const short_t nDir2 = ps2.direction(), tDir2 = 1 - nDir2;
        const bool flipped = !ifc.dirOrientation(ps1, tDir1);

        const auto& tb1 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>&>(dbasis[ps1.patch]);
        const std::vector<real_t> brks = tb1.component(tDir1).knots().breaks();
        const gsGaussRule<real_t> rule(tb1.degree(tDir1) + 2);
        gsMatrix<real_t> nodes;
        gsVector<real_t> weights;
        rule.mapToAll(brks, nodes, weights);
        const index_t nQ = nodes.cols();

        gsMatrix<real_t> pts1(2, nQ), pts2(2, nQ);
        pts1.row(nDir1).setConstant(ps1.parameter() ? 1.0 : 0.0);
        pts1.row(tDir1) = nodes.row(0);
        pts2.row(nDir2).setConstant(ps2.parameter() ? 1.0 : 0.0);
        for (index_t q = 0; q < nQ; ++q)
            pts2(tDir2, q) = flipped ? (1.0 - nodes(0, q)) : nodes(0, q);

        gsMatrix<real_t> X1, g1D1, g2D1, u1V, u1D1, u2V, u2D1;
        mp.patch(ps1.patch).eval_into(pts1, X1);
        mp.patch(ps1.patch).deriv_into(pts1, g1D1);
        mp.patch(ps2.patch).deriv_into(pts2, g2D1);

        const gsFunctionSet<real_t>& p1 = fU.piece(ps1.patch);
        const gsFunctionSet<real_t>& p2 = fU.piece(ps2.patch);
        p1.eval_into(pts1, u1V);   p1.deriv_into(pts1, u1D1);
        p2.eval_into(pts2, u2V);   p2.deriv_into(pts2, u2D1);

        std::vector<InterfaceJumpSample> curve;
        curve.reserve(static_cast<size_t>(nQ));

        const real_t outSign1 = ps1.parameter() ? 1.0 : -1.0;

        for (index_t q = 0; q < nQ; ++q)
        {
            gsSurfaceKinematics3D<real_t> k1, k2;
            if (!k1.computeFirstOrder(g1D1, q) || !k2.computeFirstOrder(g2D1, q))
                continue;
            const real_t det1 = k1.dS * k1.dS;
            const real_t det2 = k2.dS * k2.dS;
            const gsVector<real_t, 3> a1_1up = (k1.A2.squaredNorm() * k1.A1 - k1.A1.dot(k1.A2) * k1.A2) / det1;
            const gsVector<real_t, 3> a1_2up = (k1.A1.squaredNorm() * k1.A2 - k1.A1.dot(k1.A2) * k1.A1) / det1;
            const gsVector<real_t, 3> a2_1up = (k2.A2.squaredNorm() * k2.A1 - k2.A1.dot(k2.A2) * k2.A2) / det2;
            const gsVector<real_t, 3> a2_2up = (k2.A1.squaredNorm() * k2.A2 - k2.A1.dot(k2.A2) * k2.A1) / det2;

            const gsVector<real_t, 3> tang = (tDir1 == 0) ? k1.A1 : k1.A2;
            const real_t ds = tang.norm();
            const gsVector<real_t, 3> tHat = tang / ds;
            const gsVector<real_t, 3> gOut1 = outSign1 * ((nDir1 == 0) ? k1.A1 : k1.A2);
            const gsVector<real_t, 3> nCo = (gOut1 - gOut1.dot(tHat) * tHat).normalized();

            gsVector<real_t, 3> u1, u2, u1_1, u1_2, u2_1, u2_2;
            for (int c = 0; c < 3; ++c)
            {
                u1(c) = u1V(c, q);
                u2(c) = u2V(c, q);
                u1_1(c) = u1D1(2 * c + 0, q); u1_2(c) = u1D1(2 * c + 1, q);
                u2_1(c) = u2D1(2 * c + 0, q); u2_2(c) = u2D1(2 * c + 1, q);
            }

            // Directional surface derivative along the SAME 3D conormal vector nCo on both patches
            const gsVector<real_t, 3> dnU1 = nCo.dot(a1_1up) * u1_1 + nCo.dot(a1_2up) * u1_2;
            const gsVector<real_t, 3> dnU2 = nCo.dot(a2_1up) * u2_1 + nCo.dot(a2_2up) * u2_2;

            const real_t j0 = (u1 - u2).norm();
            const real_t j1 = (dnU1 - dnU2).norm();
            const real_t jRot = std::abs((dnU1 - dnU2).dot(k1.A3));

            const real_t wds = weights(q) * ds;
            c0Sq  += wds * j0 * j0;
            c1Sq  += wds * j1 * j1;
            rotSq += wds * jRot * jRot;
            c0Max  = std::max(c0Max, j0);
            c1Max  = std::max(c1Max, j1);
            rotMax = std::max(rotMax, jRot);

            InterfaceJumpSample smp;
            smp.pos << X1(0, q), X1(1, q), X1(2, q);
            smp.c0Jump = j0;
            smp.c1Jump = j1;
            smp.rotJump = jRot;
            curve.push_back(smp);
        }
        ifcCurves.push_back(std::move(curve));
    }

    c0L2  = std::sqrt(std::max(c0Sq, 0.0));
    c1L2  = std::sqrt(std::max(c1Sq, 0.0));
    rotL2 = std::sqrt(std::max(rotSq, 0.0));
    return ifcCurves;
}

/// Write ASCII VTK PolyData (.vtp) along all interior interfaces with per-method C0/C1/rotation jump scalars.
void writeInterfaceJumpsVTP(
    const std::string& filename,
    const std::vector<std::string>& methodNames,
    const std::vector<std::vector<std::vector<InterfaceJumpSample>>>& allMethodCurves)
{
    if (allMethodCurves.empty() || allMethodCurves[0].empty())
        return;

    const auto& refCurves = allMethodCurves[0];
    size_t totalPts = 0;
    for (const auto& c : refCurves) totalPts += c.size();
    if (totalPts == 0) return;

    std::ofstream ofs(filename);
    ofs << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <PolyData>\n"
        << "    <Piece NumberOfPoints=\"" << totalPts
        << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << refCurves.size()
        << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
        << "      <PointData>\n";

    for (size_t m = 0; m < methodNames.size(); ++m)
    {
        auto writeArr = [&](const std::string& suffix, int which) {
            ofs << "        <DataArray type=\"Float64\" Name=\"" << methodNames[m] << suffix
                << "\" format=\"ascii\">\n";
            ofs << std::scientific << std::setprecision(6);
            for (const auto& c : allMethodCurves[m])
                for (const auto& smp : c)
                    ofs << (which == 0 ? smp.c0Jump : (which == 1 ? smp.c1Jump : smp.rotJump)) << " ";
            ofs << "\n        </DataArray>\n";
        };
        writeArr("_c0_jump", 0);
        writeArr("_c1_jump", 1);
        writeArr("_rot_jump", 2);
    }

    ofs << "      </PointData>\n"
        << "      <Points>\n"
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n"
        << std::scientific << std::setprecision(8);
    for (const auto& c : refCurves)
        for (const auto& smp : c)
            ofs << smp.pos(0) << " " << smp.pos(1) << " " << smp.pos(2) << " ";
    ofs << "\n        </DataArray>\n"
        << "      </Points>\n"
        << "      <Lines>\n"
        << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t idx = 0;
    for (const auto& c : refCurves)
        for (size_t k = 0; k < c.size(); ++k)
            ofs << (idx++) << " ";
    ofs << "\n        </DataArray>\n"
        << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    size_t off = 0;
    for (const auto& c : refCurves)
    {
        off += c.size();
        ofs << off << " ";
    }
    ofs << "\n        </DataArray>\n"
        << "      </Lines>\n"
        << "    </Piece>\n"
        << "  </PolyData>\n"
        << "</VTKFile>\n";
    gsInfo << "Exported interface jump polylines -> " << filename << "\n";
}

/// Export all methods' displacement fields + pairwise differences onto ONE shared surface grid (.pvd).
void exportUnifiedComparisonPVD(const gsMultiPatch<real_t>& mpConst,
                                const gsMultiBasis<real_t>& dbasis,
                                const std::vector<MethodResult>& results,
                                const gsFunction<real_t>* exactFunc,
                                const std::string& exportPrefix,
                                index_t nPts = 1500)
{
    if (exportPrefix.empty() || results.empty())
        return;

    const size_t slash = exportPrefix.find_last_of("/\\");
    const std::string parentDir = (slash != std::string::npos) ? exportPrefix.substr(0, slash + 1) : "./";
    mkdirRecursive(parentDir + "_vtk");

    gsMultiPatch<real_t> mp = mpConst;
    gsExprEvaluator<real_t> ev;
    ev.setIntegrationDomain(dbasis.domain());
    auto G = ev.getMap(mp);

    gsParaviewCollection col(exportPrefix + "_comparison.pvd", &ev);
    col.options().setInt("numPoints", nPts);
    col.options().setInt("precision", 6);
    col.options().setSwitch("plotElements", true);
    col.options().setInt("plotElements.resolution", -1);
    col.options().setSwitch("plotControlNet", false);
    col.options().setSwitch("makeSubfolder", true);
    col.options().setString("subfolder", "_vtk");

    col.newTimeStep(&mp);

    // Reference method is index 0 (Ours: ASG1_MG or ASG1_Direct)
    const MethodResult& refM = results[0];
    auto uRef = ev.getVariable(refM.dispFunc());

    for (size_t i = 0; i < results.size(); ++i)
    {
        if (!results[i].ok) continue;
        auto uM = ev.getVariable(results[i].dispFunc());
        col.addField(uM, "u_" + results[i].name);
        col.addField(uM.norm(), "u_" + results[i].name + "_mag");

        if (i > 0)
        {
            col.addField(uRef - uM, "diff_" + refM.name + "_vs_" + results[i].name);
            col.addField((uRef - uM).norm(), "diff_" + refM.name + "_vs_" + results[i].name + "_mag");
        }
    }

    if (exactFunc != nullptr)
    {
        auto uEx = ev.getVariable(*exactFunc, G);
        col.addField(uEx, "u_Exact");
        col.addField(uEx.norm(), "u_Exact_mag");
        for (size_t i = 0; i < results.size(); ++i)
        {
            if (!results[i].ok) continue;
            auto uM = ev.getVariable(results[i].dispFunc());
            col.addField((uEx - uM).norm(), "err_" + results[i].name + "_mag");
        }
    }

    col.saveTimeStep();
    col.save();
    gsInfo << "Exported unified comparison PVD   -> " << exportPrefix << "_comparison.pvd\n";
}

void buildBCObjects(const gsMultiPatch<real_t>& mp,
                    const std::vector<LogicalBCEntry>& entries,
                    bool useWeakClampedForKLShell,
                    gsBoundaryConditions<real_t>& bcASG1,
                    gsBoundaryConditions<real_t>& bcKLShell)
{
    static gsConstantFunction<real_t> zero3(0.0, 3);
    static gsConstantFunction<real_t> zero1(0.0, 3);
    bcASG1 = gsBoundaryConditions<real_t>();
    bcKLShell = gsBoundaryConditions<real_t>();

    for (const auto& e : entries)
    {
        bcASG1.add(e.ps.patch, e.ps.side(), e.clamped ? "ValuesAndDerivatives" : "Values", zero3, e.comp);
        bcKLShell.addCondition(e.ps, condition_type::dirichlet, &zero1, 0, false, e.comp);
        if (e.clamped)
        {
            bcKLShell.addCondition(e.ps,
                                   useWeakClampedForKLShell ? condition_type::weak_clamped : condition_type::clamped,
                                   &zero1, 0, false, e.comp);
        }
    }
    bcKLShell.setGeoMap(mp);
}

std::vector<real_t> parsePenaltyList(const std::string& s)
{
    std::vector<real_t> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
        if (!tok.empty())
            out.push_back(std::stod(tok));
    return out;
}

std::string formatPenLabel(real_t alpha)
{
    std::ostringstream oss;
    oss << "Penalty_" << std::scientific << std::setprecision(0) << alpha;
    std::string s = oss.str();
    std::replace(s.begin(), s.end(), '+', 'p');
    std::replace(s.begin(), s.end(), '-', 'm');
    return s;
}

void runHeadToHeadCase(const std::string& caseName,
                       const gsMultiPatch<real_t>& mp,
                       const std::vector<LogicalBCEntry>& bcEntries,
                       const SharedShellForceFunction& forceFunc,
                       real_t normalPressure,
                       const gsKLShellMaterial<real_t>& mat,
                       short_t degreeP,
                       short_t regularityK,
                       index_t minRef,
                       index_t maxRef,
                       const std::vector<real_t>& penaltyScales,
                       bool runMappedC1,
                       bool runApproxC1,
                       bool pinAxialMode,
                       index_t axialComp,
                       const std::function<real_t(const gsFunctionSet<real_t>&)>& evalProbe,
                       const gsFunction<real_t>* exactFunc,
                       const std::string& exportPrefix)
{
    const short_t effK = (regularityK < 0) ? static_cast<short_t>(degreeP - 2) : regularityK;
    const bool allPoly = isAllPolynomialBSpline(mp);
    const bool zeroInnerGeoKnots = hasZeroInteriorKnots(mp);

    gsInfo << "\n==============================================================================================================\n"
           << "HEAD-TO-HEAD KIRCHHOFF-LOVE SHELL COMPARISON: " << caseName << "\n"
           << "  Patches=" << mp.nPatches() << "  Interfaces=" << mp.nInterfaces()
           << "  Boundaries=" << mp.boundaries().size()
           << "  Polynomial=" << (allPoly ? "YES" : "NO(NURBS)")
           << "  Space=S_{p=" << degreeP << ", k=" << effK << "}  levels=" << minRef << ".." << maxRef << "\n"
           << "  Material: E=" << mat.E << ", nu=" << mat.nu << ", t=" << mat.thickness
           << "  quadOrder=p+2\n"
           << "==============================================================================================================\n";

    gsBoundaryConditions<real_t> bcASG1, bcKLShellStrongClamp, bcKLShellWeakClamp;
    buildBCObjects(mp, bcEntries, false, bcASG1, bcKLShellStrongClamp);
    buildBCObjects(mp, bcEntries, true,  bcASG1, bcKLShellWeakClamp);

    std::vector<MethodResult> results;

    // ------------------------------------------------------------------
    // 1. Ours: Conforming AS-G1 (MG-PCG + Direct) on shared (mp, dbasis)
    // ------------------------------------------------------------------
    gsStopwatch sw;
    sw.restart();
    gsAsG1ShellHierarchy<real_t> hier(mp, bcASG1, degreeP, effK, minRef, maxRef, true);
    const double hierTime = sw.stop();
    const auto& fin = hier.finest();
    const gsMultiBasis<real_t> sharedBasis = fin.dbasis;

    gsAsG1KLShellAssembler<real_t> ourAsmbl(mp, sharedBasis, mat, /*quadExtra=*/2);
    gsSparseMatrix<real_t> K_discont, K_free;
    gsVector<real_t> F_discont, F_free, u_bnd;

    sw.restart();
    ourAsmbl.assembleStiffness(K_discont);
    ourAsmbl.assemblePhysicalRHS(
        [&](const gsVector<real_t>& x, gsVector<real_t, 3>& fBody, real_t& pNorm) {
            fBody = forceFunc.evalPoint(x);
            pNorm = normalPressure;
        },
        F_discont);
    ourAsmbl.projectSystem(K_discont, F_discont, fin.T_free, fin.T_bnd, u_bnd, K_free, F_free);

    if (pinAxialMode)
    {
        bool pinned = false;
        for (size_t p = 0; p < mp.nPatches() && !pinned; ++p)
            for (index_t c = 0; c < 4 && !pinned; ++c)
            {
                const index_t cOff = sumUntil(fin.argBasis[p].sizes, 9 + c);
                if (fin.mapper.is_free(cOff, static_cast<index_t>(p), axialComp))
                {
                    const index_t idx = fin.mapper.index(cOff, static_cast<index_t>(p), axialComp);
                    K_free.coeffRef(idx, idx) += mat.E;
                    pinned = true;
                }
            }
        GISMO_ENSURE(pinned, "Failed to pin axial rigid-body mode in ASG1.");
    }
    const double ourAsmTime = sw.stop();

    auto runOurSolver = [&](const std::string& label, const std::string& sType, const std::string& smType) {
        MethodResult res;
        res.name = label;
        res.nFree = K_free.rows();
        res.nnz = K_free.nonZeros();
        res.basisTimeSec = hierTime;
        res.asmTimeSec = ourAsmTime;
        gsVector<real_t> u_free;
        const auto st = hier.solve(K_free, F_free, mat, u_free, sType, smType, 1, 1, 1, 300, 1e-11);
        res.solveTimeSec = st.setupTimeSec + st.solveTimeSec;
        res.totalTimeSec = res.basisTimeSec + res.asmTimeSec + res.solveTimeSec;
        res.iterations = st.iterations;
        res.relRes = st.relResidual;
        res.ok = st.converged;
        res.usesMapped = false;
        gsVector<real_t> c_sol = fin.T_free * u_free;
        res.dispMP = ourAsmbl.buildDisplacementMultiPatch(c_sol);
        computeFieldStats(mp, sharedBasis, mat, res.dispFunc(), res.dispL2Norm, res.dispMaxMag, res.strainEnergy);
        if (evalProbe) res.probeUz = evalProbe(res.dispFunc());
        results.push_back(std::move(res));
    };

    runOurSolver("ASG1_MG",     "cg",     "line_ring1");
    runOurSolver("ASG1_Direct", "direct", "-");

#ifdef gsKLShell_ENABLED
    gsConstantFunction<real_t> tFun(mat.thickness, 3);
    gsConstantFunction<real_t> EFun(mat.E, 3);
    gsConstantFunction<real_t> nuFun(mat.nu, 3);
    gsConstantFunction<real_t> rhoFun(0.0, 3);
    gsConstantFunction<real_t> pressFun(normalPressure, 3);
    std::vector<gsFunctionSet<real_t>*> matParams(2);
    matParams[0] = &EFun;
    matParams[1] = &nuFun;

    // ------------------------------------------------------------------
    // 2. gsKLShell Weak Penalty C0/C1 on the EXACT same (mp, sharedBasis)
    // ------------------------------------------------------------------
    {
        const std::vector<real_t> penList = (mp.nInterfaces() > 0) ? penaltyScales : std::vector<real_t>{1e3};
        for (real_t alphaIfc : penList)
        {
            MethodResult res;
            res.name = (mp.nInterfaces() > 0) ? formatPenLabel(alphaIfc) : "gsKLShell_1P";
            try
            {
                gsMultiPatch<real_t> mpLocal = mp;
                gsMultiBasis<real_t> dbLocal = sharedBasis;
                gsBoundaryConditions<real_t> bcLocal = bcKLShellWeakClamp;
                bcLocal.setGeoMap(mpLocal);

                gsMaterialMatrixLinear<3, real_t> mm(mpLocal, tFun, matParams, rhoFun);
                sw.restart();
                gsThinShellAssembler<3, real_t, true> klAsm(mpLocal, dbLocal, bcLocal, forceFunc, &mm);
                klAsm.options().setReal("ExprAssembler.quA", 1.0);
                klAsm.options().setInt("ExprAssembler.quB", 2);
                klAsm.options().setReal("WeakClamped", 1e3);
                if (normalPressure != 0.0)
                    klAsm.setPressure(pressFun);

                if (mpLocal.nInterfaces() > 0)
                {
                    klAsm.options().setInt("Continuity", -1);
                    klAsm.options().setReal("IfcPenalty", alphaIfc);
                    klAsm.addWeakC0(mpLocal.topology().interfaces());
                    klAsm.addWeakC1(mpLocal.topology().interfaces());
                    klAsm.initInterfaces();
                }
                const double setupT = sw.stop();

                sw.restart();
                klAsm.assemble();
                gsSparseMatrix<real_t> Kkl = klAsm.matrix();
                gsVector<real_t> Fkl = klAsm.rhs();

                if (pinAxialMode)
                {
                    const gsDofMapper klMap = klAsm.getMapper();
                    bool pinned = false;
                    for (size_t k = 0; k < klMap.patchSize(0, axialComp) && !pinned; ++k)
                        if (klMap.is_free(static_cast<index_t>(k), 0, axialComp))
                        {
                            const index_t id = klMap.index(static_cast<index_t>(k), 0, axialComp);
                            Kkl.coeffRef(id, id) += mat.E;
                            pinned = true;
                        }
                    GISMO_ENSURE(pinned, "Failed to pin axial rigid-body mode in gsKLShell.");
                }
                const double asmT = sw.stop();

                res.nFree = Kkl.rows();
                res.nnz = Kkl.nonZeros();
                res.basisTimeSec = setupT;
                res.asmTimeSec = asmT;

                sw.restart();
                gsEigen::SimplicialLDLT<gsSparseMatrix<real_t>> ldlt;
                ldlt.compute(Kkl);
                gsVector<real_t> sol = ldlt.solve(Fkl);
                res.solveTimeSec = sw.stop();
                res.totalTimeSec = res.basisTimeSec + res.asmTimeSec + res.solveTimeSec;
                res.iterations = 1;
                const real_t fN = Fkl.norm();
                res.relRes = (fN > 0.0) ? (Fkl - Kkl * sol).norm() / fN : 0.0;
                res.ok = (ldlt.info() == gsEigen::Success) && std::isfinite(res.relRes);

                if (res.ok)
                {
                    res.usesMapped = false;
                    res.dispMP = klAsm.constructDisplacement(sol);
                    computeFieldStats(mp, sharedBasis, mat, res.dispFunc(),
                                      res.dispL2Norm, res.dispMaxMag, res.strainEnergy);
                    if (evalProbe) res.probeUz = evalProbe(res.dispFunc());
                }
            }
            catch (const std::exception& ex)
            {
                res.ok = false;
                res.note = ex.what();
                gsWarn << "Method " << res.name << " failed: " << ex.what() << "\n";
            }
            results.push_back(std::move(res));
        }
    }

#ifdef gsUnstructuredSplines_ENABLED
    // Helper to refine a polynomial copy of mp to match sharedBasis's exact (degreeP, knot vector)
    auto makeMatchedPolyGeom = [&]() -> gsMultiPatch<real_t> {
        gsMultiPatch<real_t> mpRef = mp;
        for (size_t p = 0; p < mpRef.nPatches(); ++p)
        {
            auto& bsp = dynamic_cast<gsTensorBSpline<2, real_t>&>(mpRef.patch(p));
            const auto& tbTarget = dynamic_cast<const gsTensorBSplineBasis<2, real_t>&>(sharedBasis[p]);
            for (short_t d = 0; d < 2; ++d)
            {
                const short_t dp = tbTarget.degree(d) - bsp.degree(d);
                if (dp > 0) bsp.degreeElevate(dp, d);
                const gsKnotVector<real_t>& kvCur = bsp.basis().component(d).knots();
                const gsKnotVector<real_t>& kvTgt = tbTarget.component(d).knots();
                for (real_t brk : kvTgt.breaks())
                {
                    const index_t mCur = kvCur.multiplicity(brk);
                    const index_t mTgt = kvTgt.multiplicity(brk);
                    if (mTgt > mCur)
                        bsp.insertKnot(brk, d, static_cast<int>(mTgt - mCur));
                }
            }
        }
        mpRef.computeTopology();
        return mpRef;
    };

    // ------------------------------------------------------------------
    // 3. gsKLShell + gsC1SurfSpline (Farahat et al. 2023 mapped AS-G1)
    // ------------------------------------------------------------------
    if (runMappedC1 && allPoly && mp.nInterfaces() > 0)
    {
        MethodResult res;
        res.name = "Mapped_C1";
        try
        {
            sw.restart();
            gsMultiPatch<real_t> mpC1 = makeMatchedPolyGeom();
            gsMultiBasis<real_t> dbC1(mpC1);
            gsC1SurfSpline<2, real_t> smoothC1(mpC1, dbC1);
            smoothC1.init();
            smoothC1.compute();
            gsSparseMatrix<real_t> g2l = smoothC1.getSystem().transpose();
            smoothC1.getMultiBasis(dbC1);
            gsMappedBasis<2, real_t> bb2;
            bb2.init(dbC1, g2l);
            res.basisTimeSec = sw.stop();

            gsBoundaryConditions<real_t> bcLocal = bcKLShellWeakClamp;
            bcLocal.setGeoMap(mpC1);
            gsMaterialMatrixLinear<3, real_t> mm(mpC1, tFun, matParams, rhoFun);

            sw.restart();
            gsThinShellAssembler<3, real_t, true> klAsm(mpC1, dbC1, bcLocal, forceFunc, &mm);
            klAsm.options().setReal("ExprAssembler.quA", 1.0);
            klAsm.options().setInt("ExprAssembler.quB", 2);
            klAsm.options().setInt("Continuity", -1);
            klAsm.options().setReal("WeakClamped", 1e3);
            if (normalPressure != 0.0)
                klAsm.setPressure(pressFun);
            klAsm.setSpaceBasis(bb2);
            klAsm.assemble();
            gsSparseMatrix<real_t> Kkl = klAsm.matrix();
            gsVector<real_t> Fkl = klAsm.rhs();
            res.asmTimeSec = sw.stop();

            res.nFree = Kkl.rows();
            res.nnz = Kkl.nonZeros();

            sw.restart();
            gsEigen::SimplicialLDLT<gsSparseMatrix<real_t>> ldlt;
            ldlt.compute(Kkl);
            gsVector<real_t> sol = ldlt.solve(Fkl);
            res.solveTimeSec = sw.stop();
            res.totalTimeSec = res.basisTimeSec + res.asmTimeSec + res.solveTimeSec;
            res.iterations = 1;
            const real_t fN = Fkl.norm();
            res.relRes = (fN > 0.0) ? (Fkl - Kkl * sol).norm() / fN : 0.0;
            res.ok = (ldlt.info() == gsEigen::Success) && std::isfinite(res.relRes);

            if (res.ok)
            {
                gsMatrix<real_t> solFull = klAsm.fullSolutionVector(sol);
                solFull.resize(solFull.rows() / 3, 3);
                res.usesMapped = true;
                res.dispMapped.reset(new FastMappedDisplacement());
                res.dispMapped->init(bb2, solFull);
                computeFieldStats(mp, sharedBasis, mat, res.dispFunc(),
                                  res.dispL2Norm, res.dispMaxMag, res.strainEnergy);
                if (evalProbe) res.probeUz = evalProbe(res.dispFunc());
            }
        }
        catch (const std::exception& ex)
        {
            res.ok = false;
            res.note = ex.what();
            gsWarn << "Method Mapped_C1 failed: " << ex.what() << "\n";
        }
        results.push_back(std::move(res));
    }

    // ------------------------------------------------------------------
    // 4. gsKLShell + gsApproxC1Spline (Verhelst et al. 2024 approx C1)
    // ------------------------------------------------------------------
    if (runApproxC1 && allPoly && zeroInnerGeoKnots && mp.nInterfaces() > 0)
    {
        MethodResult res;
        res.name = "Approx_C1";
        try
        {
            sw.restart();
            gsMultiPatch<real_t> mpAp = mp;
            gsMultiBasis<real_t> dbAp = sharedBasis;
            gsApproxC1Spline<2, real_t> approxC1(mpAp, dbAp);
            approxC1.options().setSwitch("interpolation", true);
            approxC1.options().setInt("gluingDataDegree", -1);
            approxC1.options().setInt("gluingDataSmoothness", -1);
            gsMappedBasis<2, real_t> bb2;
            approxC1.update(bb2);
            res.basisTimeSec = sw.stop();

            gsBoundaryConditions<real_t> bcLocal = bcKLShellWeakClamp;
            bcLocal.setGeoMap(mpAp);
            gsMaterialMatrixLinear<3, real_t> mm(mpAp, tFun, matParams, rhoFun);

            sw.restart();
            gsThinShellAssembler<3, real_t, true> klAsm(mpAp, dbAp, bcLocal, forceFunc, &mm);
            klAsm.options().setReal("ExprAssembler.quA", 1.0);
            klAsm.options().setInt("ExprAssembler.quB", 2);
            klAsm.options().setInt("Continuity", -1);
            klAsm.options().setReal("WeakClamped", 1e3);
            if (normalPressure != 0.0)
                klAsm.setPressure(pressFun);
            klAsm.setSpaceBasis(bb2);
            klAsm.assemble();
            gsSparseMatrix<real_t> Kkl = klAsm.matrix();
            gsVector<real_t> Fkl = klAsm.rhs();
            res.asmTimeSec = sw.stop();

            res.nFree = Kkl.rows();
            res.nnz = Kkl.nonZeros();

            sw.restart();
            gsEigen::SimplicialLDLT<gsSparseMatrix<real_t>> ldlt;
            ldlt.compute(Kkl);
            gsVector<real_t> sol = ldlt.solve(Fkl);
            res.solveTimeSec = sw.stop();
            res.totalTimeSec = res.basisTimeSec + res.asmTimeSec + res.solveTimeSec;
            res.iterations = 1;
            const real_t fN = Fkl.norm();
            res.relRes = (fN > 0.0) ? (Fkl - Kkl * sol).norm() / fN : 0.0;
            res.ok = (ldlt.info() == gsEigen::Success) && std::isfinite(res.relRes);

            if (res.ok)
            {
                gsMatrix<real_t> solFull = klAsm.fullSolutionVector(sol);
                solFull.resize(solFull.rows() / 3, 3);
                res.usesMapped = true;
                res.dispMapped.reset(new FastMappedDisplacement());
                res.dispMapped->init(bb2, solFull);
                computeFieldStats(mp, sharedBasis, mat, res.dispFunc(),
                                  res.dispL2Norm, res.dispMaxMag, res.strainEnergy);
                if (evalProbe) res.probeUz = evalProbe(res.dispFunc());
            }
        }
        catch (const std::exception& ex)
        {
            res.ok = false;
            res.note = ex.what();
            gsWarn << "Method Approx_C1 failed: " << ex.what() << "\n";
        }
        results.push_back(std::move(res));
    }
#endif // gsUnstructuredSplines_ENABLED
#endif // gsKLShell_ENABLED

    // ------------------------------------------------------------------
    // 5. Compute Interface Jumps & Pairwise Cross-Method Value Differences
    // ------------------------------------------------------------------
    std::vector<std::string> okNames;
    std::vector<std::vector<std::vector<InterfaceJumpSample>>> allCurves;
    for (auto& res : results)
    {
        if (!res.ok) continue;
        auto curves = computeInterfaceJumps(mp, sharedBasis, res.dispFunc(),
                                            res.c0JumpL2, res.c0JumpMax,
                                            res.c1JumpL2, res.c1JumpMax,
                                            res.rotJumpL2, res.rotJumpMax);
        okNames.push_back(res.name);
        allCurves.push_back(std::move(curves));
    }

    // Print Per-Method Summary Table
    gsInfo << "\n--- Method Performance, Energy & Interface Continuity Summary ---\n"
           << std::left
           << std::setw(14) << "Method"
           << std::setw(8)  << "nFree"
           << std::setw(10) << "nnz"
           << std::setw(6)  << "Iter"
           << std::setw(9)  << "Basis(s)"
           << std::setw(9)  << "Asm(s)"
           << std::setw(9)  << "Solve(s)"
           << std::setw(9)  << "Tot(s)"
           << std::setw(12) << "Energy"
           << std::setw(11) << "||u||_max"
           << std::setw(11) << "ProbeUz"
           << std::setw(10) << "[[u]]_L2"
           << std::setw(10) << "[[dn u]]L2"
           << "\n" << std::string(128, '-') << "\n";

    for (const auto& r : results)
    {
        if (!r.ok)
        {
            gsInfo << std::left << std::setw(14) << r.name << "FAILED (" << r.note << ")\n";
            continue;
        }
        gsInfo << std::left
               << std::setw(14) << r.name
               << std::setw(8)  << r.nFree
               << std::setw(10) << r.nnz
               << std::setw(6)  << r.iterations
               << std::fixed << std::setprecision(4)
               << std::setw(9)  << r.basisTimeSec
               << std::setw(9)  << r.asmTimeSec
               << std::setw(9)  << r.solveTimeSec
               << std::setw(9)  << r.totalTimeSec
               << std::scientific << std::setprecision(3)
               << std::setw(12) << r.strainEnergy
               << std::setw(11) << r.dispMaxMag
               << std::fixed << std::setprecision(6)
               << std::setw(11) << r.probeUz
               << std::scientific << std::setprecision(2)
               << std::setw(10) << r.c0JumpL2
               << std::setw(10) << r.c1JumpL2
               << "\n";
    }

    // Print Pairwise Value Differences vs Reference (Ours: ASG1_MG)
    std::vector<PairwiseDiff> diffs;
    if (!results.empty() && results[0].ok)
    {
        gsInfo << "\n--- Pairwise Value Differences vs " << results[0].name << " ---\n"
               << std::left
               << std::setw(16) << "ComparedTo"
               << std::setw(11) << "Linf_abs"
               << std::setw(11) << "Linf_rel"
               << std::setw(11) << "L2_rel"
               << std::setw(11) << "H1_rel"
               << std::setw(11) << "H2_rel"
               << std::setw(11) << "En_rel"
               << std::setw(28) << "MaxDiffPos(x,y,z)"
               << "\n" << std::string(110, '-') << "\n";

        for (size_t i = 1; i < results.size(); ++i)
        {
            if (!results[i].ok) continue;
            PairwiseDiff d = computePairwiseDiff(mp, sharedBasis, mat,
                                                 results[0].dispFunc(), results[i].dispFunc(),
                                                 results[0].name, results[i].name);
            std::ostringstream posStr;
            posStr << "(" << std::fixed << std::setprecision(2)
                   << d.maxDiffPos(0) << "," << d.maxDiffPos(1) << "," << d.maxDiffPos(2) << ")";
            gsInfo << std::left
                   << std::setw(16) << d.nameB
                   << std::scientific << std::setprecision(2)
                   << std::setw(11) << d.lInfDiff
                   << std::setw(11) << d.lInfRel
                   << std::setw(11) << d.l2Rel
                   << std::setw(11) << d.h1Rel
                   << std::setw(11) << d.h2Rel
                   << std::setw(11) << d.enRel
                   << std::setw(28) << posStr.str()
                   << "\n";
            diffs.push_back(d);
        }
    }

    // Export Unified PVD, Interface Jump VTP, and CSV Summary
    if (!exportPrefix.empty())
    {
        const size_t slash = exportPrefix.find_last_of("/\\");
        if (slash != std::string::npos)
            mkdirRecursive(exportPrefix.substr(0, slash + 1));

        exportUnifiedComparisonPVD(mp, sharedBasis, results, exactFunc, exportPrefix);
        if (mp.nInterfaces() > 0)
            writeInterfaceJumpsVTP(exportPrefix + "_interfaces.vtp", okNames, allCurves);

        std::ofstream csv(exportPrefix + "_summary.csv");
        csv << "Method,nFree,nnz,Iter,BasisSec,AsmSec,SolveSec,TotalSec,StrainEnergy,MaxMag,ProbeUz,C0JumpL2,C1JumpL2,RotJumpL2\n";
        for (const auto& r : results)
        {
            if (!r.ok) continue;
            csv << r.name << "," << r.nFree << "," << r.nnz << "," << r.iterations << ","
                << std::setprecision(8) << r.basisTimeSec << "," << r.asmTimeSec << ","
                << r.solveTimeSec << "," << r.totalTimeSec << ","
                << r.strainEnergy << "," << r.dispMaxMag << "," << r.probeUz << ","
                << r.c0JumpL2 << "," << r.c1JumpL2 << "," << r.rotJumpL2 << "\n";
        }
        csv << "\nPairA,PairB,LinfAbs,LinfRel,L2Rel,H1Rel,H2Rel,EnRel,MaxX,MaxY,MaxZ\n";
        for (const auto& d : diffs)
        {
            csv << d.nameA << "," << d.nameB << "," << std::setprecision(8)
                << d.lInfDiff << "," << d.lInfRel << "," << d.l2Rel << ","
                << d.h1Rel << "," << d.h2Rel << "," << d.enRel << ","
                << d.maxDiffPos(0) << "," << d.maxDiffPos(1) << "," << d.maxDiffPos(2) << "\n";
        }
        gsInfo << "Exported comparison CSV summary   -> " << exportPrefix << "_summary.csv\n";
    }
}

} // anonymous namespace

int main(int argc, char* argv[])
{
    std::string geomFile = "asg1_reparametrized_domains/farahat_kl_shells/xml/surfHalfSphere_asg1.xml";
    std::string testMode = "compare";
    std::string bcMode = "simply_supported";
    std::string loadMode = "smooth3d";
    std::string penaltyStr = "1e1,1e3,1e5";
    std::string exportPrefix = "";
    int degreeP = 3;
    int regularityK = -1;
    int minRef = 0;
    int maxRef = 2;
    real_t E = 1e6;
    real_t nu = 0.3;
    real_t thickness = 0.1;
    bool split2x2 = false;
    bool skipMappedC1 = false;
    bool skipApproxC1 = false;

    gsCmdLine cmd("Head-to-Head KL Shell Method Comparator (AS-G1 MG vs gsKLShell Penalty / Mapped C1 / Approx C1).");
    cmd.addString("f", "file", "Input 3D surface XML file.", geomFile);
    cmd.addString("t", "test", "Mode: compare | scordelis_lo | plate_exact.", testMode);
    cmd.addString("b", "bc", "Boundary condition: simply_supported | clamped.", bcMode);
    cmd.addString("l", "load", "Load: smooth3d | gravity | pressure.", loadMode);
    cmd.addString("", "penalties", "Comma-separated gsKLShell IfcPenalty values.", penaltyStr);
    cmd.addString("o", "export", "Output prefix for unified _comparison.pvd, _interfaces.vtp, _summary.csv.", exportPrefix);
    cmd.addInt("p", "degree", "Spline degree p >= 3.", degreeP);
    cmd.addInt("k", "regularity", "Interior regularity 1 <= k <= p-2 (-1 => p-2).", regularityK);
    cmd.addInt("r", "minRef", "Coarsest multigrid level.", minRef);
    cmd.addInt("R", "maxRef", "Finest refinement level.", maxRef);
    cmd.addReal("", "young", "Young's modulus E.", E);
    cmd.addReal("", "nu", "Poisson ratio nu.", nu);
    cmd.addReal("", "thickness", "Shell thickness t.", thickness);
    cmd.addSwitch("split2x2", "Uniformly split each patch 2x2 before comparing.", split2x2);
    cmd.addSwitch("skipMappedC1", "Skip gsC1SurfSpline baseline.", skipMappedC1);
    cmd.addSwitch("skipApproxC1", "Skip gsApproxC1Spline baseline.", skipApproxC1);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    const std::vector<real_t> penalties = parsePenaltyList(penaltyStr);

    if (testMode == "plate_exact")
    {
        // Exact analytical strong-form Kirchhoff plate on 2x2 unit square [0,1]^2 in z = 0
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5);
        mp.embed(3);
        normalizePatchSupports(mp);
        mp.computeTopology();

        const bool isClamped = (bcMode == "clamped");
        gsKLShellMaterial<real_t> mat(E, nu, thickness);
        const real_t D = mat.E * std::pow(mat.thickness, 3) / (12.0 * (1.0 - mat.nu * mat.nu));
        const real_t pi = EIGEN_PI;
        const real_t pi4 = pi * pi * pi * pi;

        std::vector<LogicalBCEntry> bcEntries;
        for (const auto& ps : mp.boundaries())
            for (short_t c = 0; c < 3; ++c)
                bcEntries.push_back(LogicalBCEntry{ps, c, isClamped});

        SharedShellForceFunction forceFn([=](const gsVector<real_t>& x) -> gsVector<real_t, 3> {
            gsVector<real_t, 3> f = gsVector<real_t, 3>::Zero();
            if (!isClamped)
            {
                f(2) = 4.0 * pi4 * D * std::sin(pi * x(0)) * std::sin(pi * x(1));
            }
            else
            {
                const real_t sx = std::sin(pi * x(0)), sy = std::sin(pi * x(1));
                const real_t c2x = std::cos(2.0 * pi * x(0)), c2y = std::cos(2.0 * pi * x(1));
                f(2) = 8.0 * pi4 * D * (c2x * c2y - c2x * sy * sy - sx * sx * c2y);
            }
            return f;
        });

        ExactPlateDisplacementFunction exactFunc(isClamped);
        auto probeCenter = [&](const gsFunctionSet<real_t>& fU) -> real_t {
            for (size_t p = 0; p < mp.nPatches(); ++p)
                for (int k = 0; k < 4; ++k)
                {
                    gsMatrix<real_t> u(2, 1);
                    u << (k % 2), (k / 2);
                    const gsMatrix<real_t> X = mp.patch(p).eval(u);
                    if (std::abs(X(0, 0) - 0.5) < 1e-10 && std::abs(X(1, 0) - 0.5) < 1e-10)
                    {
                        gsMatrix<real_t> v;
                        fU.piece(static_cast<index_t>(p)).eval_into(u, v);
                        return v(2, 0);
                    }
                }
            return 0.0;
        };

        runHeadToHeadCase(isClamped ? "2x2 Clamped Plate (Exact u_z = sin^2(pi x) sin^2(pi y))"
                                    : "2x2 Simply-Supported Plate (Exact u_z = sin(pi x) sin(pi y))",
                          mp, bcEntries, forceFn, 0.0, mat, degreeP, regularityK, minRef, maxRef,
                          penalties, !skipMappedC1, !skipApproxC1, false, 0,
                          probeCenter, &exactFunc, exportPrefix);
        return 0;
    }

    if (testMode == "scordelis_lo")
    {
        const std::string scFile = (geomFile == "asg1_reparametrized_domains/farahat_kl_shells/xml/surfHalfSphere_asg1.xml")
                                 ? "filedata/surfaces/scordelis_lo_roof.xml"
                                 : geomFile;
        gsMultiPatch<real_t> mp = loadSurface3D(scFile, split2x2);
        gsKLShellMaterial<real_t> mat(4.32e8, 0.0, 0.25);

        gsVector<real_t, 3> bbMin = gsVector<real_t, 3>::Constant(1e30);
        gsVector<real_t, 3> bbMax = gsVector<real_t, 3>::Constant(-1e30);
        gsMatrix<real_t> grid(2, 9);
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
            {
                grid(0, i + 3 * j) = i * 0.5;
                grid(1, i + 3 * j) = j * 0.5;
            }
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsMatrix<real_t> X = mp.patch(p).eval(grid);
            for (int c = 0; c < 3; ++c)
            {
                bbMin(c) = std::min(bbMin(c), X.row(c).minCoeff());
                bbMax(c) = std::max(bbMax(c), X.row(c).maxCoeff());
            }
        }
        const gsVector<real_t, 3> bbExt = bbMax - bbMin;
        const real_t scale = std::max(bbExt.maxCoeff(), real_t(1e-12));

        std::vector<LogicalBCEntry> bcEntries;
        std::vector<patchSide> freeSides;
        index_t axialComp = 0;
        for (const auto& ps : mp.boundaries())
        {
            gsMatrix<real_t> pts(2, 3);
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            pts.row(nDir).setConstant(ps.parameter() ? 1.0 : 0.0);
            pts(tDir, 0) = 0.0; pts(tDir, 1) = 0.5; pts(tDir, 2) = 1.0;
            const gsMatrix<real_t> X = mp.patch(ps.patch).eval(pts);
            const real_t dX = X.row(0).maxCoeff() - X.row(0).minCoeff();
            const real_t dY = X.row(1).maxCoeff() - X.row(1).minCoeff();
            const real_t dZ = X.row(2).maxCoeff() - X.row(2).minCoeff();
            if (dZ > 0.05 * bbExt(2) && (dY < 1e-3 * scale || dX < 1e-3 * scale))
            {
                const short_t tComp = (dY < 1e-3 * scale) ? 0 : 1;
                axialComp = 1 - tComp;
                bcEntries.push_back(LogicalBCEntry{ps, tComp, false});
                bcEntries.push_back(LogicalBCEntry{ps, 2, false});
            }
            else
                freeSides.push_back(ps);
        }

        const real_t midAx = 0.5 * (bbMin(axialComp) + bbMax(axialComp));
        auto probeFreeMid = [&](const gsFunctionSet<real_t>& fU) -> real_t {
            real_t minUz = 0.0;
            bool found = false;
            for (const auto& ps : freeSides)
            {
                const short_t nDir = ps.direction(), tDir = 1 - nDir;
                const real_t ss[3] = {0.0, 0.5, 1.0};
                for (int s = 0; s < 3; ++s)
                {
                    gsMatrix<real_t> pt(2, 1);
                    pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                    pt(tDir, 0) = ss[s];
                    const gsMatrix<real_t> X = mp.patch(ps.patch).eval(pt);
                    if (std::abs(X(axialComp, 0) - midAx) <= 1e-3 * scale)
                    {
                        gsMatrix<real_t> v;
                        fU.piece(ps.patch).eval_into(pt, v);
                        if (!found || v(2, 0) < minUz)
                        {
                            minUz = v(2, 0);
                            found = true;
                        }
                    }
                }
            }
            return minUz;
        };

        SharedShellForceFunction gravFn([](const gsVector<real_t>&) -> gsVector<real_t, 3> {
            gsVector<real_t, 3> f; f << 0.0, 0.0, -90.0; return f;
        });

        runHeadToHeadCase(split2x2 ? "Scordelis-Lo Roof (2x2 Split Exact NURBS, ref u_z = -0.3006)"
                                   : "Scordelis-Lo Roof (1-Patch Exact NURBS, ref u_z = -0.3006)",
                          mp, bcEntries, gravFn, 0.0, mat, degreeP, regularityK, minRef, maxRef,
                          penalties, false, false, true, axialComp,
                          probeFreeMid, nullptr, exportPrefix);
        return 0;
    }

    // Default mode: general multi-patch curved 3D shell head-to-head comparison
    gsMultiPatch<real_t> mp = loadSurface3D(geomFile, split2x2);
    gsKLShellMaterial<real_t> mat(E, nu, thickness);
    const bool isClamped = (bcMode == "clamped");

    std::vector<LogicalBCEntry> bcEntries;
    for (const auto& ps : mp.boundaries())
        for (short_t c = 0; c < 3; ++c)
            bcEntries.push_back(LogicalBCEntry{ps, c, isClamped});

    real_t pNorm = 0.0;
    SharedShellForceFunction::ForceFn fLambda;
    if (loadMode == "gravity")
    {
        fLambda = [](const gsVector<real_t>&) -> gsVector<real_t, 3> {
            gsVector<real_t, 3> f; f << 0.0, 0.0, -90.0; return f;
        };
    }
    else if (loadMode == "pressure")
    {
        pNorm = 10.0;
        fLambda = [](const gsVector<real_t>&) -> gsVector<real_t, 3> {
            return gsVector<real_t, 3>::Zero();
        };
    }
    else // "smooth3d"
    {
        fLambda = [](const gsVector<real_t>& x) -> gsVector<real_t, 3> {
            gsVector<real_t, 3> f;
            f << 50.0 * std::sin(x(0)) * std::cos(x(1)),
                 50.0 * std::cos(x(0)) * std::sin(x(2)),
                -90.0 + 40.0 * std::cos(x(0) + x(1));
            return f;
        };
    }
    SharedShellForceFunction forceFn(fLambda);

    auto probePatch0Mid = [&](const gsFunctionSet<real_t>& fU) -> real_t {
        gsMatrix<real_t> pt(2, 1); pt << 0.5, 0.5;
        gsMatrix<real_t> v;
        fU.piece(0).eval_into(pt, v);
        return v(2, 0);
    };

    runHeadToHeadCase(geomFile + " [BC=" + bcMode + ", load=" + loadMode + "]",
                      mp, bcEntries, forceFn, pNorm, mat, degreeP, regularityK, minRef, maxRef,
                      penalties, !skipMappedC1, !skipApproxC1, false, 0,
                      probePatch0Mid, nullptr, exportPrefix);
    return 0;
}
