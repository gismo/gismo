/** @file gsAsG1KLShellAssembler.hpp

    @brief Thread-safe multi-patch Kirchhoff-Love (Koiter) shell assembler over
           non-isoparametric polynomial AS-G1 spaces [V_A^1]^3 on 3D polynomial
           or rational NURBS surfaces.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiPatch.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsIO/gsParaviewCollection.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <queue>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace gismo {

/**
 * @brief Isotropic linear elastic Kirchhoff-Love shell material parameters.
 */
template <typename T>
struct gsKLShellMaterial
{
    T E = T(1e6);
    T nu = T(0.3);
    T thickness = T(0.1);

    gsKLShellMaterial() = default;
    gsKLShellMaterial(T E_, T nu_, T thickness_)
        : E(E_), nu(nu_), thickness(thickness_) {}

    /// Contravariant plane-stress constitutive matrix C in Voigt order (11, 22, 12)
    /// corresponding to strain vector (eps_11, eps_22, 2*eps_12)^T.
    void voigtMatrix(T A11inv, T A22inv, T A12inv, gsMatrix<T, 3, 3>& C) const
    {
        const T fac = E / (T(1) - nu * nu);
        C(0, 0) = fac * A11inv * A11inv;
        C(1, 1) = fac * A22inv * A22inv;
        C(0, 1) = fac * (nu * A11inv * A22inv + (T(1) - nu) * A12inv * A12inv);
        C(1, 0) = C(0, 1);
        C(0, 2) = fac * A11inv * A12inv;
        C(2, 0) = C(0, 2);
        C(1, 2) = fac * A22inv * A12inv;
        C(2, 1) = C(1, 2);
        C(2, 2) = fac * (T(0.5) * (T(1) - nu) * A11inv * A22inv
                       + T(0.5) * (T(1) + nu) * A12inv * A12inv);
    }
};

/**
 * @brief Surface differential geometry kinematics evaluated at a quadrature point.
 */
template <typename T>
struct gsSurfaceKinematics3D
{
    gsVector<T, 3> A1, A2, A11, A22, A12;
    gsVector<T, 3> A3;      // unit surface normal
    gsVector<T, 3> A1up, A2up; // contravariant tangent basis A^1, A^2
    T dS;                   // area measure sqrt(det A)
    T A11inv, A22inv, A12inv;
    // Christoffel symbols Gamma^gamma_{alpha beta} = A_{alpha,beta} . A^gamma
    T G1_11, G2_11;
    T G1_22, G2_22;
    T G1_12, G2_12;
    // Second fundamental form B_{alpha beta} = A_{alpha,beta} . A_3
    T B11, B22, B12;

    /// Lightweight first-order surface kinematics (A1, A2, dS, A3 only).
    bool computeFirstOrder(const gsMatrix<T>& gD1, index_t q)
    {
        for (int c = 0; c < 3; ++c)
        {
            A1(c) = gD1(2 * c + 0, q);
            A2(c) = gD1(2 * c + 1, q);
        }
        const T a11 = A1.dot(A1);
        const T a22 = A2.dot(A2);
        const T a12 = A1.dot(A2);
        const T detA = a11 * a22 - a12 * a12;
        if (!(detA > T(0)) || !std::isfinite(detA))
            return false;

        dS = std::sqrt(detA);
        A3 = A1.cross(A2) / dS;
        return true;
    }

    bool compute(const gsMatrix<T>& gD1, const gsMatrix<T>& gD2, index_t q)
    {
        for (int c = 0; c < 3; ++c)
        {
            A1(c)  = gD1(2 * c + 0, q);
            A2(c)  = gD1(2 * c + 1, q);
            A11(c) = gD2(3 * c + 0, q);
            A22(c) = gD2(3 * c + 1, q);
            A12(c) = gD2(3 * c + 2, q);
        }
        const T a11 = A1.dot(A1);
        const T a22 = A2.dot(A2);
        const T a12 = A1.dot(A2);
        const T detA = a11 * a22 - a12 * a12;
        if (!(detA > T(0)) || !std::isfinite(detA))
            return false;

        dS = std::sqrt(detA);
        const T invDet = T(1) / detA;
        A11inv =  a22 * invDet;
        A22inv =  a11 * invDet;
        A12inv = -a12 * invDet;

        A1up = A11inv * A1 + A12inv * A2;
        A2up = A12inv * A1 + A22inv * A2;

        A3 = A1.cross(A2) / dS;

        G1_11 = A11.dot(A1up); G2_11 = A11.dot(A2up);
        G1_22 = A22.dot(A1up); G2_22 = A22.dot(A2up);
        G1_12 = A12.dot(A1up); G2_12 = A12.dot(A2up);

        B11 = A11.dot(A3);
        B22 = A22.dot(A3);
        B12 = A12.dot(A3);
        return true;
    }
};

/**
 * @brief Exact surface error norms for a 3-component displacement field.
 */
template <typename T>
struct gsKLShellErrors
{
    T l2Err = T(0), l2Norm = T(0);
    T h1Err = T(0), h1Norm = T(0);
    T h2Err = T(0), h2Norm = T(0);
    T energyErr = T(0), energyNorm = T(0);

    T l2Rel()     const { return (l2Norm > T(0)) ? l2Err / l2Norm : l2Err; }
    T h1Rel()     const { return (h1Norm > T(0)) ? h1Err / h1Norm : h1Err; }
    T h2Rel()     const { return (h2Norm > T(0)) ? h2Err / h2Norm : h2Err; }
    T energyRel() const { return (energyNorm > T(0)) ? energyErr / energyNorm : energyErr; }
};

/**
 * @brief Thread-safe disjoint Kirchhoff-Love shell assembler over polynomial
 *        multi-basis S_{p,k,Z} mapped by polynomial or rational NURBS geometry mp.
 */
template <typename T>
class gsAsG1KLShellAssembler
{
public:
    gsAsG1KLShellAssembler(const gsMultiPatch<T>& mp,
                           const gsMultiBasis<T>& dbasis,
                           const gsKLShellMaterial<T>& mat,
                           index_t quadExtra = 2)
        : m_mp(mp), m_dbasis(dbasis), m_mat(mat), m_quadExtra(quadExtra)
    {
        GISMO_ENSURE(m_mp.nPatches() > 0,
                     "gsAsG1KLShellAssembler requires at least one patch.");
        GISMO_ENSURE(m_mp.targetDim() == 3 && m_mp.parDim() == 2,
                     "gsAsG1KLShellAssembler requires a 2D surface in R^3.");
        GISMO_ENSURE(m_dbasis.nBases() == m_mp.nPatches(),
                     "gsAsG1KLShellAssembler: dbasis patch count (" << m_dbasis.nBases()
                     << ") does not match geometry patch count (" << m_mp.nPatches() << ").");
        GISMO_ENSURE(m_mat.E > T(0) && m_mat.thickness > T(0) &&
                     m_mat.nu > T(-1) && m_mat.nu < T(0.5),
                     "gsAsG1KLShellAssembler: invalid material parameters (E=" << m_mat.E
                     << ", nu=" << m_mat.nu << ", thickness=" << m_mat.thickness << ").");
        GISMO_ENSURE(m_quadExtra >= 0,
                     "gsAsG1KLShellAssembler: quadExtra must be non-negative.");

        const size_t nP = m_mp.nPatches();
        m_patchOff.assign(nP + 1, 0);
        m_brkU.resize(nP);
        m_brkV.resize(nP);

        for (size_t p = 0; p < nP; ++p)
        {
            const auto* tb = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&m_dbasis[p]);
            GISMO_ENSURE(tb != nullptr,
                         "gsAsG1KLShellAssembler requires gsTensorBSplineBasis<2, T> on patch " << p << ".");
            const gsMatrix<T> bSupp = tb->support();
            const gsMatrix<T> gSupp = m_mp.patch(p).support();
            const T suppScale = std::max<T>(T(1), gSupp.cwiseAbs().maxCoeff());
            GISMO_ENSURE((bSupp - gSupp).cwiseAbs().maxCoeff() <= T(1e-12) * suppScale,
                         "gsAsG1KLShellAssembler: parametric support mismatch between dbasis and geometry on patch "
                         << p << ".");

            m_patchOff[p + 1] = m_patchOff[p] + tb->size();
            m_brkU[p] = tb->component(0).knots().breaks();
            m_brkV[p] = tb->component(1).knots().breaks();
            const index_t nEu = static_cast<index_t>(m_brkU[p].size()) - 1;
            const index_t nEv = static_cast<index_t>(m_brkV[p].size()) - 1;
            for (index_t ev = 0; ev < nEv; ++ev)
                for (index_t eu = 0; eu < nEu; ++eu)
                    m_elements.push_back(PatchElement{static_cast<index_t>(p), eu, ev});
        }
        m_nDiscontScalar = m_patchOff[nP];

        // Determine consistent per-patch surface normal signs across interfaces via BFS
        // so follower normal pressure p_normal * A3 is globally consistent on orientable surfaces.
        m_patchNormalSign.assign(nP, T(1));
        std::vector<bool> visited(nP, false);
        std::vector<std::vector<std::pair<size_t, T>>> adj(nP);
        for (auto it = m_mp.iBegin(); it != m_mp.iEnd(); ++it)
        {
            const patchSide ps1 = it->first();
            const patchSide ps2 = it->second();
            const size_t p1 = static_cast<size_t>(ps1.patch);
            const size_t p2 = static_cast<size_t>(ps2.patch);
            const gsGeometry<T>& g1 = m_mp.patch(p1);
            const gsGeometry<T>& g2 = m_mp.patch(p2);
            const gsMatrix<T> s1 = g1.support();
            const gsMatrix<T> s2 = g2.support();

            gsMatrix<T> u1(2, 1), u2(2, 1);
            const short_t n1 = ps1.direction(), t1 = 1 - n1;
            const short_t n2 = ps2.direction(), t2 = 1 - n2;
            u1(n1, 0) = ps1.parameter() ? s1(n1, 1) : s1(n1, 0);
            u1(t1, 0) = T(0.5) * (s1(t1, 0) + s1(t1, 1));
            u2(n2, 0) = ps2.parameter() ? s2(n2, 1) : s2(n2, 0);
            u2(t2, 0) = T(0.5) * (s2(t2, 0) + s2(t2, 1));

            gsMatrix<T> d1, d2;
            g1.deriv_into(u1, d1);
            g2.deriv_into(u2, d2);
            gsSurfaceKinematics3D<T> k1, k2;
            if (k1.computeFirstOrder(d1, 0) && k2.computeFirstOrder(d2, 0))
            {
                const T relSign = (k1.A3.dot(k2.A3) < T(0)) ? T(-1) : T(1);
                adj[p1].push_back({p2, relSign});
                adj[p2].push_back({p1, relSign});
            }
        }
        for (size_t seed = 0; seed < nP; ++seed)
        {
            if (visited[seed]) continue;
            std::queue<size_t> q;
            visited[seed] = true;
            m_patchNormalSign[seed] = T(1);
            q.push(seed);
            while (!q.empty())
            {
                const size_t u = q.front();
                q.pop();
                for (const auto& edge : adj[u])
                {
                    if (!visited[edge.first])
                    {
                        visited[edge.first] = true;
                        m_patchNormalSign[edge.first] = m_patchNormalSign[u] * edge.second;
                        q.push(edge.first);
                    }
                }
            }
        }
    }

    index_t scalarDiscontSize() const { return m_nDiscontScalar; }
    index_t vectorDiscontSize() const { return 3 * m_nDiscontScalar; }
    const std::vector<index_t>& patchOffsets() const { return m_patchOff; }
    const std::vector<T>& patchNormalSigns() const { return m_patchNormalSign; }

    /**
     * @brief Assemble the disjoint 3-component Kirchhoff-Love stiffness matrix
     *        K_discont of size (3 * nDiscontScalar) x (3 * nDiscontScalar).
     *
     * Row/column ordering: comp * nDiscontScalar + patchOff[p] + localDof.
     */
    void assembleStiffness(gsSparseMatrix<T>& K_discont) const
    {
        const index_t N = 3 * m_nDiscontScalar;
        const index_t nElems = static_cast<index_t>(m_elements.size());
        const T tMem = m_mat.thickness;
        const T tBend = (m_mat.thickness * m_mat.thickness * m_mat.thickness) / T(12);

        int maxThreads = 1;
#ifdef _OPENMP
        maxThreads = omp_get_max_threads();
#endif
        std::vector<gsSparseMatrix<T>> threadMats(static_cast<size_t>(maxThreads));
        std::atomic<index_t> degenerateCount(0);

        #pragma omp parallel
        {
            int tid = 0;
            int nT = 1;
#ifdef _OPENMP
            tid = omp_get_thread_num();
            nT = omp_get_num_threads();
#endif
            gsSparseEntries<T> locEnt;
            const auto& tb0 = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[0]);
            const index_t nActEst = (tb0.component(0).degree() + 1) * (tb0.component(1).degree() + 1);
            locEnt.reserve(static_cast<size_t>((nElems + nT - 1) / nT) * 9 * nActEst * nActEst);

            gsMatrix<T> qu, qv, pts;
            gsVector<T> wu, wv, wq;
            gsMatrix<T> gD1, gD2;
            std::vector<gsMatrix<T>> ders;
            gsMatrix<index_t> act;
            gsMatrix<T> Ke, Bm, Bb;
            gsMatrix<T, 3, 3> C;
            index_t locDegenerate = 0;

            #pragma omp for schedule(dynamic, 4)
            for (index_t eIdx = 0; eIdx < nElems; ++eIdx)
            {
                const PatchElement& el = m_elements[static_cast<size_t>(eIdx)];
                const size_t p = static_cast<size_t>(el.patch);
                const auto& tb = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
                const gsGeometry<T>& geo = m_mp.patch(p);
                const index_t pOff = m_patchOff[p];

                const gsBSplineBasis<T>& bu = tb.component(0);
                const gsBSplineBasis<T>& bv = tb.component(1);
                const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra);
                const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra);
                const index_t nQu = ruleU.numNodes();
                const index_t nQv = ruleV.numNodes();
                const index_t nQ = nQu * nQv;

                pts.resize(2, nQ);
                wq.resize(nQ);
                ruleU.mapTo(m_brkU[p][el.eu], m_brkU[p][el.eu + 1], qu, wu);
                ruleV.mapTo(m_brkV[p][el.ev], m_brkV[p][el.ev + 1], qv, wv);
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
                tb.evalAllDers_into(pts, 2, ders);
                tb.active_into(pts.col(0), act);
                const index_t nAct = act.rows();
                const index_t nEl = 3 * nAct;
                Ke.setZero(nEl, nEl);
                Bm.resize(3, nEl);
                Bb.resize(3, nEl);

                for (index_t q = 0; q < nQ; ++q)
                {
                    gsSurfaceKinematics3D<T> kin;
                    if (!kin.compute(gD1, gD2, q))
                    {
                        ++locDegenerate;
                        continue;
                    }
                    m_mat.voigtMatrix(kin.A11inv, kin.A22inv, kin.A12inv, C);

                    for (index_t a = 0; a < nAct; ++a)
                    {
                        const T d1  = ders[1](2 * a + 0, q);
                        const T d2  = ders[1](2 * a + 1, q);
                        const T d11 = ders[2](3 * a + 0, q);
                        const T d22 = ders[2](3 * a + 1, q);
                        const T d12 = ders[2](3 * a + 2, q);

                        // Covariant scalar second derivatives N_{a; alpha beta}
                        const T cov11 = -(d11 - kin.G1_11 * d1 - kin.G2_11 * d2);
                        const T cov22 = -(d22 - kin.G1_22 * d1 - kin.G2_22 * d2);
                        const T cov12 = -T(2) * (d12 - kin.G1_12 * d1 - kin.G2_12 * d2);

                        for (index_t c = 0; c < 3; ++c)
                        {
                            const index_t col = c * nAct + a;
                            Bm(0, col) = d1 * kin.A1(c);
                            Bm(1, col) = d2 * kin.A2(c);
                            Bm(2, col) = d2 * kin.A1(c) + d1 * kin.A2(c);

                            Bb(0, col) = cov11 * kin.A3(c);
                            Bb(1, col) = cov22 * kin.A3(c);
                            Bb(2, col) = cov12 * kin.A3(c);
                        }
                    }

                    const T wDS = wq(q) * kin.dS;
                    Ke.noalias() += (wDS * tMem)  * (Bm.transpose() * C * Bm);
                    Ke.noalias() += (wDS * tBend) * (Bb.transpose() * C * Bb);
                }

                for (index_t c1 = 0; c1 < 3; ++c1)
                {
                    const index_t base1 = c1 * m_nDiscontScalar + pOff;
                    for (index_t a = 0; a < nAct; ++a)
                    {
                        const index_t row = base1 + act(a, 0);
                        const index_t rLoc = c1 * nAct + a;
                        for (index_t c2 = 0; c2 < 3; ++c2)
                        {
                            const index_t base2 = c2 * m_nDiscontScalar + pOff;
                            for (index_t b = 0; b < nAct; ++b)
                            {
                                const T val = Ke(rLoc, c2 * nAct + b);
                                if (val != T(0))
                                    locEnt.add(row, base2 + act(b, 0), val);
                            }
                        }
                    }
                }
            }

            if (locDegenerate > 0)
                degenerateCount.fetch_add(locDegenerate, std::memory_order_relaxed);

            threadMats[static_cast<size_t>(tid)].resize(N, N);
            threadMats[static_cast<size_t>(tid)].setFrom(locEnt);
        }

        GISMO_ENSURE(degenerateCount.load() == 0,
                     "gsAsG1KLShellAssembler::assembleStiffness encountered "
                     << degenerateCount.load()
                     << " degenerate/non-finite surface metric quadrature point(s).");

        K_discont.swap(threadMats[0]);
        for (size_t t = 1; t < threadMats.size(); ++t)
            if (threadMats[t].nonZeros() > 0)
                K_discont += threadMats[t];
        K_discont.makeCompressed();
    }

    /**
     * @brief Assemble the exact consistent variational RHS F_discont(v) = a(u_ex, v)
     *        for ANY smooth manufactured 3D displacement field u_ex(x,y,z).
     *
     * @note When compiled with OpenMP, @p evalExactComp must be thread-safe for
     *       concurrent read-only invocation from multiple threads.
     *
     * ExactVecFn signature:
     *   void evalExactComp(index_t comp, const gsVector<T>& xPhys,
     *                      T& val, gsVector<T>& gradPhys, gsMatrix<T>& hessPhys)
     */
    template <typename ExactVecFn>
    void assembleManufacturedRHS(const ExactVecFn& evalExactComp,
                                 gsVector<T>& F_discont) const
    {
        const index_t N = 3 * m_nDiscontScalar;
        const index_t nElems = static_cast<index_t>(m_elements.size());
        F_discont.setZero(N);
        const T tMem = m_mat.thickness;
        const T tBend = (m_mat.thickness * m_mat.thickness * m_mat.thickness) / T(12);
        std::atomic<index_t> degenerateCount(0);

        #pragma omp parallel
        {
            int nT = 1;
#ifdef _OPENMP
            nT = omp_get_num_threads();
#endif
            gsVector<T> fLoc;
            gsVector<T>* fTarget = &F_discont;
            if (nT > 1)
            {
                fLoc.setZero(N);
                fTarget = &fLoc;
            }

            gsMatrix<T> qu, qv, pts;
            gsVector<T> wu, wv, wq;
            gsMatrix<T> Xq, gD1, gD2;
            std::vector<gsMatrix<T>> ders;
            gsMatrix<index_t> act;
            gsMatrix<T, 3, 3> C;
            gsVector<T> gP(3);
            gsMatrix<T> hP(3, 3);
            gsVector<T> Fe;
            index_t locDegenerate = 0;

            #pragma omp for schedule(dynamic, 4)
            for (index_t eIdx = 0; eIdx < nElems; ++eIdx)
            {
                const PatchElement& el = m_elements[static_cast<size_t>(eIdx)];
                const size_t p = static_cast<size_t>(el.patch);
                const auto& tb = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
                const gsGeometry<T>& geo = m_mp.patch(p);
                const index_t pOff = m_patchOff[p];

                const gsBSplineBasis<T>& bu = tb.component(0);
                const gsBSplineBasis<T>& bv = tb.component(1);
                const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra);
                const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra);
                const index_t nQu = ruleU.numNodes();
                const index_t nQv = ruleV.numNodes();
                const index_t nQ = nQu * nQv;

                pts.resize(2, nQ);
                wq.resize(nQ);
                ruleU.mapTo(m_brkU[p][el.eu], m_brkU[p][el.eu + 1], qu, wu);
                ruleV.mapTo(m_brkV[p][el.ev], m_brkV[p][el.ev + 1], qv, wv);
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
                tb.evalAllDers_into(pts, 2, ders);
                tb.active_into(pts.col(0), act);
                const index_t nAct = act.rows();

                Fe.setZero(3 * nAct);

                for (index_t q = 0; q < nQ; ++q)
                {
                    gsSurfaceKinematics3D<T> kin;
                    if (!kin.compute(gD1, gD2, q))
                    {
                        ++locDegenerate;
                        continue;
                    }
                    m_mat.voigtMatrix(kin.A11inv, kin.A22inv, kin.A12inv, C);

                    gsVector<T> xPhys = Xq.col(q);
                    gsVector<T, 3> u1, u2, u11, u22, u12;
                    for (index_t c = 0; c < 3; ++c)
                    {
                        T val = T(0);
                        evalExactComp(c, xPhys, val, gP, hP);
                        u1(c)  = gP.dot(kin.A1);
                        u2(c)  = gP.dot(kin.A2);
                        u11(c) = kin.A1.dot(hP * kin.A1) + gP.dot(kin.A11);
                        u22(c) = kin.A2.dot(hP * kin.A2) + gP.dot(kin.A22);
                        u12(c) = kin.A1.dot(hP * kin.A2) + gP.dot(kin.A12);
                    }

                    // Exact Voigt membrane strain eps(u_ex) and curvature change kappa(u_ex)
                    gsVector<T, 3> epsEx, kapEx;
                    epsEx(0) = kin.A1.dot(u1);
                    epsEx(1) = kin.A2.dot(u2);
                    epsEx(2) = kin.A1.dot(u2) + kin.A2.dot(u1);

                    kapEx(0) = -(u11 - kin.G1_11 * u1 - kin.G2_11 * u2).dot(kin.A3);
                    kapEx(1) = -(u22 - kin.G1_22 * u1 - kin.G2_22 * u2).dot(kin.A3);
                    kapEx(2) = -T(2) * (u12 - kin.G1_12 * u1 - kin.G2_12 * u2).dot(kin.A3);

                    const gsVector<T, 3> Nforce = (wq(q) * kin.dS * tMem)  * (C * epsEx);
                    const gsVector<T, 3> Mmom   = (wq(q) * kin.dS * tBend) * (C * kapEx);

                    for (index_t a = 0; a < nAct; ++a)
                    {
                        const T d1  = ders[1](2 * a + 0, q);
                        const T d2  = ders[1](2 * a + 1, q);
                        const T d11 = ders[2](3 * a + 0, q);
                        const T d22 = ders[2](3 * a + 1, q);
                        const T d12 = ders[2](3 * a + 2, q);

                        const T cov11 = -(d11 - kin.G1_11 * d1 - kin.G2_11 * d2);
                        const T cov22 = -(d22 - kin.G1_22 * d1 - kin.G2_22 * d2);
                        const T cov12 = -T(2) * (d12 - kin.G1_12 * d1 - kin.G2_12 * d2);
                        const T mScalar = cov11 * Mmom(0) + cov22 * Mmom(1) + cov12 * Mmom(2);

                        for (index_t c = 0; c < 3; ++c)
                        {
                            const T bm0 = d1 * kin.A1(c);
                            const T bm1 = d2 * kin.A2(c);
                            const T bm2 = d2 * kin.A1(c) + d1 * kin.A2(c);
                            Fe(c * nAct + a) += bm0 * Nforce(0) + bm1 * Nforce(1) + bm2 * Nforce(2)
                                              + mScalar * kin.A3(c);
                        }
                    }
                }

                for (index_t c = 0; c < 3; ++c)
                {
                    const index_t base = c * m_nDiscontScalar + pOff;
                    for (index_t a = 0; a < nAct; ++a)
                        (*fTarget)(base + act(a, 0)) += Fe(c * nAct + a);
                }
            }

            if (locDegenerate > 0)
                degenerateCount.fetch_add(locDegenerate, std::memory_order_relaxed);

            if (nT > 1)
            {
                #pragma omp critical(gsAsG1KL_mfgRHS)
                F_discont += fLoc;
            }
        }

        GISMO_ENSURE(degenerateCount.load() == 0,
                     "gsAsG1KLShellAssembler::assembleManufacturedRHS encountered "
                     << degenerateCount.load()
                     << " degenerate/non-finite surface metric quadrature point(s).");
    }

    /**
     * @brief Assemble physical surface load RHS:
     *        F(v) = int_S (f_body(X) + p_normal(X) * n) . v dS,
     *        where n = patchNormalSign[p] * A3 is globally oriented across interfaces.
     *
     * @note When compiled with OpenMP, @p evalForceAndPressure must be thread-safe
     *       for concurrent read-only invocation from multiple threads.
     */
    template <typename LoadFn>
    void assemblePhysicalRHS(const LoadFn& evalForceAndPressure,
                             gsVector<T>& F_discont) const
    {
        const index_t N = 3 * m_nDiscontScalar;
        const index_t nElems = static_cast<index_t>(m_elements.size());
        F_discont.setZero(N);
        std::atomic<index_t> degenerateCount(0);

        #pragma omp parallel
        {
            int nT = 1;
#ifdef _OPENMP
            nT = omp_get_num_threads();
#endif
            gsVector<T> fLoc;
            gsVector<T>* fTarget = &F_discont;
            if (nT > 1)
            {
                fLoc.setZero(N);
                fTarget = &fLoc;
            }

            gsMatrix<T> qu, qv, pts;
            gsVector<T> wu, wv, wq;
            gsMatrix<T> Xq, gD1, vals;
            gsMatrix<index_t> act;
            index_t locDegenerate = 0;

            #pragma omp for schedule(dynamic, 4)
            for (index_t eIdx = 0; eIdx < nElems; ++eIdx)
            {
                const PatchElement& el = m_elements[static_cast<size_t>(eIdx)];
                const size_t p = static_cast<size_t>(el.patch);
                const auto& tb = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
                const gsGeometry<T>& geo = m_mp.patch(p);
                const index_t pOff = m_patchOff[p];
                const T nSign = m_patchNormalSign[p];

                const gsBSplineBasis<T>& bu = tb.component(0);
                const gsBSplineBasis<T>& bv = tb.component(1);
                const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra);
                const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra);
                const index_t nQu = ruleU.numNodes(), nQv = ruleV.numNodes();
                const index_t nQ = nQu * nQv;

                pts.resize(2, nQ);
                wq.resize(nQ);
                ruleU.mapTo(m_brkU[p][el.eu], m_brkU[p][el.eu + 1], qu, wu);
                ruleV.mapTo(m_brkV[p][el.ev], m_brkV[p][el.ev + 1], qv, wv);
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
                tb.eval_into(pts, vals);
                tb.active_into(pts.col(0), act);
                const index_t nAct = act.rows();

                for (index_t q = 0; q < nQ; ++q)
                {
                    gsSurfaceKinematics3D<T> kin;
                    if (!kin.computeFirstOrder(gD1, q))
                    {
                        ++locDegenerate;
                        continue;
                    }
                    gsVector<T> xPhys = Xq.col(q);
                    gsVector<T, 3> fBody = gsVector<T, 3>::Zero();
                    T pNorm = T(0);
                    evalForceAndPressure(xPhys, fBody, pNorm);
                    const gsVector<T, 3> load = (wq(q) * kin.dS) * (fBody + (pNorm * nSign) * kin.A3);

                    for (index_t c = 0; c < 3; ++c)
                    {
                        const index_t base = c * m_nDiscontScalar + pOff;
                        for (index_t a = 0; a < nAct; ++a)
                            (*fTarget)(base + act(a, 0)) += vals(a, q) * load(c);
                    }
                }
            }

            if (locDegenerate > 0)
                degenerateCount.fetch_add(locDegenerate, std::memory_order_relaxed);

            if (nT > 1)
            {
                #pragma omp critical(gsAsG1KL_physRHS)
                F_discont += fLoc;
            }
        }

        GISMO_ENSURE(degenerateCount.load() == 0,
                     "gsAsG1KLShellAssembler::assemblePhysicalRHS encountered "
                     << degenerateCount.load()
                     << " degenerate/non-finite surface metric quadrature point(s).");
    }

    /**
     * @brief Project disjoint stiffness and load onto the conforming free AS-G1 space
     *        and incorporate inhomogeneous boundary lifting T_bnd * u_bnd.
     *
     * Pass either u_bnd of size T_bnd.cols(), or an empty u_bnd (size 0) for
     * homogeneous zero boundary conditions. Any other non-zero size mismatch throws.
     */
    static void projectSystem(const gsSparseMatrix<T>& K_discont,
                              const gsVector<T>& F_discont,
                              const gsSparseMatrix<T>& T_free,
                              const gsSparseMatrix<T>& T_bnd,
                              const gsVector<T>& u_bnd,
                              gsSparseMatrix<T>& K_free,
                              gsVector<T>& F_free)
    {
        GISMO_ENSURE(K_discont.rows() == K_discont.cols(),
                     "projectSystem: K_discont must be square.");
        GISMO_ENSURE(T_free.rows() == K_discont.rows(),
                     "projectSystem: T_free row count (" << T_free.rows()
                     << ") does not match K_discont (" << K_discont.rows() << ").");
        GISMO_ENSURE(F_discont.size() == K_discont.rows(),
                     "projectSystem: F_discont size (" << F_discont.size()
                     << ") does not match K_discont (" << K_discont.rows() << ").");

        K_free = T_free.transpose() * K_discont * T_free;
        K_free.makeCompressed();

        if (T_bnd.cols() > 0)
        {
            GISMO_ENSURE(T_bnd.rows() == K_discont.cols(),
                         "projectSystem: T_bnd row count (" << T_bnd.rows()
                         << ") does not match K_discont (" << K_discont.cols() << ").");
            GISMO_ENSURE(u_bnd.size() == 0 || u_bnd.size() == T_bnd.cols(),
                         "projectSystem: u_bnd size (" << u_bnd.size()
                         << ") must be either 0 (homogeneous) or T_bnd.cols() ("
                         << T_bnd.cols() << ").");
            if (u_bnd.size() == T_bnd.cols() && u_bnd.squaredNorm() > T(0))
            {
                const gsVector<T> c_bnd = T_bnd * u_bnd;
                F_free = T_free.transpose() * (F_discont - K_discont * c_bnd);
                return;
            }
        }
        else
        {
            GISMO_ENSURE(u_bnd.size() == 0,
                         "projectSystem: non-empty u_bnd passed when T_bnd has 0 columns.");
        }
        F_free = T_free.transpose() * F_discont;
    }

    /**
     * @brief Compute exact surface L^2, H^1, covariant H^2, and shell energy errors
     *        of a discontinuous B-spline coefficient vector c_discont against u_ex.
     *
     * @note When compiled with OpenMP, @p evalExactComp must be thread-safe.
     */
    template <typename ExactVecFn>
    gsKLShellErrors<T> computeErrors(const gsVector<T>& c_discont,
                                     const ExactVecFn& evalExactComp) const
    {
        GISMO_ENSURE(c_discont.size() == 3 * m_nDiscontScalar,
                     "computeErrors: c_discont size (" << c_discont.size()
                     << ") does not match vectorDiscontSize() (" << (3 * m_nDiscontScalar) << ").");

        const index_t nElems = static_cast<index_t>(m_elements.size());
        const T tMem = m_mat.thickness;
        const T tBend = (m_mat.thickness * m_mat.thickness * m_mat.thickness) / T(12);

        T l2E2 = 0, l2N2 = 0, h1E2 = 0, h1N2 = 0, h2E2 = 0, h2N2 = 0, enE2 = 0, enN2 = 0;
        index_t degenerateCount = 0;

        #pragma omp parallel for reduction(+:l2E2,l2N2,h1E2,h1N2,h2E2,h2N2,enE2,enN2,degenerateCount) schedule(dynamic, 4)
        for (index_t eIdx = 0; eIdx < nElems; ++eIdx)
        {
            const PatchElement& el = m_elements[static_cast<size_t>(eIdx)];
            const size_t p = static_cast<size_t>(el.patch);
            const auto& tb = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
            const gsGeometry<T>& geo = m_mp.patch(p);
            const index_t pOff = m_patchOff[p];

            const gsBSplineBasis<T>& bu = tb.component(0);
            const gsBSplineBasis<T>& bv = tb.component(1);
            const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra + 1);
            const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra + 1);
            const index_t nQu = ruleU.numNodes(), nQv = ruleV.numNodes();
            const index_t nQ = nQu * nQv;

            gsMatrix<T> qu, qv, pts(2, nQ);
            gsVector<T> wu, wv, wq(nQ);
            gsMatrix<T> Xq, gD1, gD2;
            std::vector<gsMatrix<T>> ders;
            gsMatrix<index_t> act;
            gsMatrix<T, 3, 3> C;
            gsVector<T> gP(3);
            gsMatrix<T> hP(3, 3);

            ruleU.mapTo(m_brkU[p][el.eu], m_brkU[p][el.eu + 1], qu, wu);
            ruleV.mapTo(m_brkV[p][el.ev], m_brkV[p][el.ev + 1], qv, wv);
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
            tb.evalAllDers_into(pts, 2, ders);
            tb.active_into(pts.col(0), act);
            const index_t nAct = act.rows();

            for (index_t q = 0; q < nQ; ++q)
            {
                gsSurfaceKinematics3D<T> kin;
                if (!kin.compute(gD1, gD2, q))
                {
                    ++degenerateCount;
                    continue;
                }
                m_mat.voigtMatrix(kin.A11inv, kin.A22inv, kin.A12inv, C);

                gsVector<T> xPhys = Xq.col(q);
                gsVector<T, 3> uEx, uEx1, uEx2, uEx11, uEx22, uEx12;
                gsVector<T, 3> uH = gsVector<T, 3>::Zero();
                gsVector<T, 3> uH1 = gsVector<T, 3>::Zero(), uH2 = gsVector<T, 3>::Zero();
                gsVector<T, 3> uH11 = gsVector<T, 3>::Zero(), uH22 = gsVector<T, 3>::Zero(), uH12 = gsVector<T, 3>::Zero();

                for (index_t c = 0; c < 3; ++c)
                {
                    T val = T(0);
                    evalExactComp(c, xPhys, val, gP, hP);
                    uEx(c)   = val;
                    uEx1(c)  = gP.dot(kin.A1);
                    uEx2(c)  = gP.dot(kin.A2);
                    uEx11(c) = kin.A1.dot(hP * kin.A1) + gP.dot(kin.A11);
                    uEx22(c) = kin.A2.dot(hP * kin.A2) + gP.dot(kin.A22);
                    uEx12(c) = kin.A1.dot(hP * kin.A2) + gP.dot(kin.A12);

                    const index_t base = c * m_nDiscontScalar + pOff;
                    for (index_t a = 0; a < nAct; ++a)
                    {
                        const T coef = c_discont(base + act(a, 0));
                        uH(c)   += coef * ders[0](a, q);
                        uH1(c)  += coef * ders[1](2 * a + 0, q);
                        uH2(c)  += coef * ders[1](2 * a + 1, q);
                        uH11(c) += coef * ders[2](3 * a + 0, q);
                        uH22(c) += coef * ders[2](3 * a + 1, q);
                        uH12(c) += coef * ders[2](3 * a + 2, q);
                    }
                }

                const gsVector<T, 3> e   = uEx   - uH;
                const gsVector<T, 3> e1  = uEx1  - uH1;
                const gsVector<T, 3> e2  = uEx2  - uH2;
                const gsVector<T, 3> e11 = uEx11 - uH11;
                const gsVector<T, 3> e22 = uEx22 - uH22;
                const gsVector<T, 3> e12 = uEx12 - uH12;

                const T wDS = wq(q) * kin.dS;
                l2E2 += wDS * e.squaredNorm();
                l2N2 += wDS * uEx.squaredNorm();

                auto surfH1Sq = [&](const gsVector<T, 3>& d1, const gsVector<T, 3>& d2) -> T {
                    return kin.A11inv * d1.squaredNorm()
                         + T(2) * kin.A12inv * d1.dot(d2)
                         + kin.A22inv * d2.squaredNorm();
                };
                h1E2 += wDS * surfH1Sq(e1, e2);
                h1N2 += wDS * surfH1Sq(uEx1, uEx2);

                // Surface covariant H^2 seminorm: A^{alpha mu} A^{beta nu} (v_{;alpha beta} . v_{;mu nu})
                auto surfH2Sq = [&](const gsVector<T, 3>& d1, const gsVector<T, 3>& d2,
                                    const gsVector<T, 3>& d11, const gsVector<T, 3>& d22,
                                    const gsVector<T, 3>& d12) -> T {
                    const gsVector<T, 3> c11 = d11 - kin.G1_11 * d1 - kin.G2_11 * d2;
                    const gsVector<T, 3> c22 = d22 - kin.G1_22 * d1 - kin.G2_22 * d2;
                    const gsVector<T, 3> c12 = d12 - kin.G1_12 * d1 - kin.G2_12 * d2;
                    const T a11 = kin.A11inv, a22 = kin.A22inv, a12 = kin.A12inv;
                    T sum = T(0);
                    for (int comp = 0; comp < 3; ++comp)
                    {
                        // Frobenius norm of A^{-1/2} H A^{-1/2} = trace(A^{-1} H A^{-1} H)
                        const T h00 = c11(comp), h11 = c22(comp), h01 = c12(comp);
                        const T m00 = a11 * h00 + a12 * h01;
                        const T m01 = a11 * h01 + a12 * h11;
                        const T m10 = a12 * h00 + a22 * h01;
                        const T m11 = a12 * h01 + a22 * h11;
                        sum += m00 * m00 + T(2) * m01 * m10 + m11 * m11;
                    }
                    return sum;
                };
                h2E2 += wDS * surfH2Sq(e1, e2, e11, e22, e12);
                h2N2 += wDS * surfH2Sq(uEx1, uEx2, uEx11, uEx22, uEx12);

                // Shell energy norm
                auto shellEnDens = [&](const gsVector<T, 3>& d1, const gsVector<T, 3>& d2,
                                       const gsVector<T, 3>& d11, const gsVector<T, 3>& d22,
                                       const gsVector<T, 3>& d12) -> T {
                    gsVector<T, 3> epsV, kapV;
                    epsV(0) = kin.A1.dot(d1);
                    epsV(1) = kin.A2.dot(d2);
                    epsV(2) = kin.A1.dot(d2) + kin.A2.dot(d1);
                    kapV(0) = -(d11 - kin.G1_11 * d1 - kin.G2_11 * d2).dot(kin.A3);
                    kapV(1) = -(d22 - kin.G1_22 * d1 - kin.G2_22 * d2).dot(kin.A3);
                    kapV(2) = -T(2) * (d12 - kin.G1_12 * d1 - kin.G2_12 * d2).dot(kin.A3);
                    return tMem * epsV.dot(C * epsV) + tBend * kapV.dot(C * kapV);
                };
                enE2 += wDS * shellEnDens(e1, e2, e11, e22, e12);
                enN2 += wDS * shellEnDens(uEx1, uEx2, uEx11, uEx22, uEx12);
            }
        }

        GISMO_ENSURE(degenerateCount == 0,
                     "gsAsG1KLShellAssembler::computeErrors encountered " << degenerateCount
                     << " degenerate/non-finite surface metric quadrature point(s).");

        gsKLShellErrors<T> out;
        out.l2Err      = std::sqrt(std::max(l2E2, T(0)));
        out.l2Norm     = std::sqrt(std::max(l2N2, T(0)));
        out.h1Err      = std::sqrt(std::max(h1E2, T(0)));
        out.h1Norm     = std::sqrt(std::max(h1N2, T(0)));
        out.h2Err      = std::sqrt(std::max(h2E2, T(0)));
        out.h2Norm     = std::sqrt(std::max(h2N2, T(0)));
        out.energyErr  = std::sqrt(std::max(enE2, T(0)));
        out.energyNorm = std::sqrt(std::max(enN2, T(0)));
        return out;
    }

    /**
     * @brief Evaluate displacement vector u_h at parametric point (u, v) on patch p.
     */
    gsVector<T, 3> evalDisplacement(const gsVector<T>& c_discont,
                                    size_t patchIdx, T u, T v) const
    {
        GISMO_ENSURE(c_discont.size() == 3 * m_nDiscontScalar,
                     "evalDisplacement: c_discont size (" << c_discont.size()
                     << ") does not match vectorDiscontSize() (" << (3 * m_nDiscontScalar) << ").");
        GISMO_ENSURE(patchIdx < m_mp.nPatches(),
                     "evalDisplacement: patchIdx (" << patchIdx
                     << ") out of range [0, " << m_mp.nPatches() << ").");

        const auto& tb = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[patchIdx]);
        gsMatrix<T> pt(2, 1);
        pt << u, v;
        gsMatrix<T> vals = tb.eval(pt);
        gsMatrix<index_t> act = tb.active(pt);
        gsVector<T, 3> disp = gsVector<T, 3>::Zero();
        const index_t pOff = m_patchOff[patchIdx];
        for (index_t c = 0; c < 3; ++c)
        {
            const index_t base = c * m_nDiscontScalar + pOff;
            for (index_t a = 0; a < act.rows(); ++a)
                disp(c) += c_discont(base + act(a, 0)) * vals(a, 0);
        }
        return disp;
    }

    /**
     * @brief Reconstruct the 3-component displacement gsMultiPatch<T> in the
     *        polynomial spline basis dbasis (useful for ParaView export).
     */
    gsMultiPatch<T> buildDisplacementMultiPatch(const gsVector<T>& c_discont) const
    {
        GISMO_ENSURE(c_discont.size() == 3 * m_nDiscontScalar,
                     "buildDisplacementMultiPatch: c_discont size (" << c_discont.size()
                     << ") does not match vectorDiscontSize() (" << (3 * m_nDiscontScalar) << ").");

        gsMultiPatch<T> dispMP;
        for (size_t p = 0; p < m_mp.nPatches(); ++p)
        {
            const auto& tb = static_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
            const index_t nLoc = tb.size();
            const index_t pOff = m_patchOff[p];
            gsMatrix<T> coefs(nLoc, 3);
            for (index_t c = 0; c < 3; ++c)
                coefs.col(c) = c_discont.segment(c * m_nDiscontScalar + pOff, nLoc);
            dispMP.addPatch(typename gsGeometry<T>::uPtr(new gsTensorBSpline<2, T>(tb, coefs)));
        }
        return dispMP;
    }

private:
    struct PatchElement
    {
        index_t patch;
        index_t eu;
        index_t ev;
    };

    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_dbasis;
    gsKLShellMaterial<T> m_mat;
    index_t m_quadExtra;
    std::vector<index_t> m_patchOff;
    index_t m_nDiscontScalar;
    std::vector<std::vector<T>> m_brkU;
    std::vector<std::vector<T>> m_brkV;
    std::vector<PatchElement> m_elements;
    std::vector<T> m_patchNormalSign;
};

} // namespace gismo
