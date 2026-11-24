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

#include <cmath>
#include <vector>

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
        if (detA <= T(1e-28))
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
        GISMO_ENSURE(mp.targetDim() == 3 && mp.parDim() == 2,
                     "gsAsG1KLShellAssembler requires a 2D surface in R^3.");
        const size_t nP = mp.nPatches();
        m_patchOff.assign(nP + 1, 0);
        for (size_t p = 0; p < nP; ++p)
            m_patchOff[p + 1] = m_patchOff[p] + dbasis[p].size();
        m_nDiscontScalar = m_patchOff[nP];
    }

    index_t scalarDiscontSize() const { return m_nDiscontScalar; }
    index_t vectorDiscontSize() const { return 3 * m_nDiscontScalar; }
    const std::vector<index_t>& patchOffsets() const { return m_patchOff; }

    /**
     * @brief Assemble the disjoint 3-component Kirchhoff-Love stiffness matrix
     *        K_discont of size (3 * nDiscontScalar) x (3 * nDiscontScalar).
     *
     * Row/column ordering: comp * nDiscontScalar + patchOff[p] + localDof.
     */
    void assembleStiffness(gsSparseMatrix<T>& K_discont) const
    {
        const size_t nP = m_mp.nPatches();
        const index_t N = 3 * m_nDiscontScalar;
        const T tMem = m_mat.thickness;
        const T tBend = (m_mat.thickness * m_mat.thickness * m_mat.thickness) / T(12);

        std::vector<gsSparseEntries<T>> patchEntries(nP);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
            const gsGeometry<T>& geo = m_mp.patch(p);
            const index_t pOff = m_patchOff[p];

            const gsBSplineBasis<T>& bu = tb.component(0);
            const gsBSplineBasis<T>& bv = tb.component(1);
            const std::vector<T> brkU = bu.knots().breaks();
            const std::vector<T> brkV = bv.knots().breaks();
            const index_t nEu = static_cast<index_t>(brkU.size()) - 1;
            const index_t nEv = static_cast<index_t>(brkV.size()) - 1;

            const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra);
            const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra);
            const index_t nQu = ruleU.numNodes();
            const index_t nQv = ruleV.numNodes();
            const index_t nQ = nQu * nQv;

            gsMatrix<T> qu, qv, pts(2, nQ);
            gsVector<T> wu, wv, wq(nQ);
            gsMatrix<T> gD1, gD2;
            std::vector<gsMatrix<T>> ders;
            gsMatrix<index_t> act;

            gsSparseEntries<T>& locEnt = patchEntries[p];
            const index_t nActEst = (bu.degree() + 1) * (bv.degree() + 1);
            locEnt.reserve(nEu * nEv * 9 * nActEst * nActEst / 2);

            gsMatrix<T> Ke;
            gsMatrix<T, 3, 3> C;

            for (index_t ev = 0; ev < nEv; ++ev)
            {
                ruleV.mapTo(brkV[ev], brkV[ev + 1], qv, wv);
                for (index_t eu = 0; eu < nEu; ++eu)
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
                    tb.evalAllDers_into(pts, 2, ders);
                    tb.active_into(pts.col(0), act);
                    const index_t nAct = act.rows();
                    const index_t nEl = 3 * nAct;
                    Ke.setZero(nEl, nEl);

                    gsMatrix<T> Bm(3, nEl), Bb(3, nEl);

                    for (index_t q = 0; q < nQ; ++q)
                    {
                        gsSurfaceKinematics3D<T> kin;
                        if (!kin.compute(gD1, gD2, q))
                            continue;
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
            }
        }

        gsSparseEntries<T> allEnt;
        size_t totalNnz = 0;
        for (size_t p = 0; p < nP; ++p)
            totalNnz += patchEntries[p].size();
        allEnt.reserve(totalNnz);
        for (size_t p = 0; p < nP; ++p)
            for (const auto& t : patchEntries[p])
                allEnt.add(t.row(), t.col(), t.value());

        K_discont.resize(N, N);
        K_discont.setFrom(allEnt);
    }

    /**
     * @brief Assemble the exact consistent variational RHS F_discont(v) = a(u_ex, v)
     *        for ANY smooth manufactured 3D displacement field u_ex(x,y,z).
     *
     * ExactVecFn signature:
     *   void evalExactComp(index_t comp, const gsVector<T>& xPhys,
     *                      T& val, gsVector<T>& gradPhys, gsMatrix<T>& hessPhys)
     */
    template <typename ExactVecFn>
    void assembleManufacturedRHS(const ExactVecFn& evalExactComp,
                                 gsVector<T>& F_discont) const
    {
        const size_t nP = m_mp.nPatches();
        const index_t N = 3 * m_nDiscontScalar;
        F_discont.setZero(N);
        const T tMem = m_mat.thickness;
        const T tBend = (m_mat.thickness * m_mat.thickness * m_mat.thickness) / T(12);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
            const gsGeometry<T>& geo = m_mp.patch(p);
            const index_t pOff = m_patchOff[p];

            const gsBSplineBasis<T>& bu = tb.component(0);
            const gsBSplineBasis<T>& bv = tb.component(1);
            const std::vector<T> brkU = bu.knots().breaks();
            const std::vector<T> brkV = bv.knots().breaks();
            const index_t nEu = static_cast<index_t>(brkU.size()) - 1;
            const index_t nEv = static_cast<index_t>(brkV.size()) - 1;

            const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra);
            const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra);
            const index_t nQu = ruleU.numNodes();
            const index_t nQv = ruleV.numNodes();
            const index_t nQ = nQu * nQv;

            gsMatrix<T> qu, qv, pts(2, nQ);
            gsVector<T> wu, wv, wq(nQ);
            gsMatrix<T> Xq, gD1, gD2;
            std::vector<gsMatrix<T>> ders;
            gsMatrix<index_t> act;
            gsMatrix<T, 3, 3> C;

            gsVector<T> gP(3);
            gsMatrix<T> hP(3, 3);

            for (index_t ev = 0; ev < nEv; ++ev)
            {
                ruleV.mapTo(brkV[ev], brkV[ev + 1], qv, wv);
                for (index_t eu = 0; eu < nEu; ++eu)
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
                    tb.evalAllDers_into(pts, 2, ders);
                    tb.active_into(pts.col(0), act);
                    const index_t nAct = act.rows();

                    gsVector<T> Fe = gsVector<T>::Zero(3 * nAct);

                    for (index_t q = 0; q < nQ; ++q)
                    {
                        gsSurfaceKinematics3D<T> kin;
                        if (!kin.compute(gD1, gD2, q))
                            continue;
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
                            F_discont(base + act(a, 0)) += Fe(c * nAct + a);
                    }
                }
            }
        }
    }

    /**
     * @brief Assemble physical surface load RHS:
     *        F(v) = int_S (f_body(X) + p_normal(X) * A3) . v dS.
     */
    template <typename LoadFn>
    void assemblePhysicalRHS(const LoadFn& evalForceAndPressure,
                             gsVector<T>& F_discont) const
    {
        const size_t nP = m_mp.nPatches();
        F_discont.setZero(3 * m_nDiscontScalar);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
            const gsGeometry<T>& geo = m_mp.patch(p);
            const index_t pOff = m_patchOff[p];

            const gsBSplineBasis<T>& bu = tb.component(0);
            const gsBSplineBasis<T>& bv = tb.component(1);
            const std::vector<T> brkU = bu.knots().breaks();
            const std::vector<T> brkV = bv.knots().breaks();
            const index_t nEu = static_cast<index_t>(brkU.size()) - 1;
            const index_t nEv = static_cast<index_t>(brkV.size()) - 1;

            const gsGaussRule<T> ruleU(bu.degree() + m_quadExtra);
            const gsGaussRule<T> ruleV(bv.degree() + m_quadExtra);
            const index_t nQu = ruleU.numNodes(), nQv = ruleV.numNodes();
            const index_t nQ = nQu * nQv;

            gsMatrix<T> qu, qv, pts(2, nQ);
            gsVector<T> wu, wv, wq(nQ);
            gsMatrix<T> Xq, gD1, gD2, vals;
            gsMatrix<index_t> act;

            for (index_t ev = 0; ev < nEv; ++ev)
            {
                ruleV.mapTo(brkV[ev], brkV[ev + 1], qv, wv);
                for (index_t eu = 0; eu < nEu; ++eu)
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
                    tb.eval_into(pts, vals);
                    tb.active_into(pts.col(0), act);
                    const index_t nAct = act.rows();

                    for (index_t q = 0; q < nQ; ++q)
                    {
                        gsSurfaceKinematics3D<T> kin;
                        if (!kin.compute(gD1, gD2, q))
                            continue;
                        gsVector<T> xPhys = Xq.col(q);
                        gsVector<T, 3> fBody = gsVector<T, 3>::Zero();
                        T pNorm = T(0);
                        evalForceAndPressure(xPhys, fBody, pNorm);
                        const gsVector<T, 3> load = (wq(q) * kin.dS) * (fBody + pNorm * kin.A3);

                        for (index_t c = 0; c < 3; ++c)
                        {
                            const index_t base = c * m_nDiscontScalar + pOff;
                            for (index_t a = 0; a < nAct; ++a)
                                F_discont(base + act(a, 0)) += vals(a, q) * load(c);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Project disjoint stiffness and load onto the conforming free AS-G1 space
     *        and incorporate inhomogeneous boundary lifting T_bnd * u_bnd.
     */
    static void projectSystem(const gsSparseMatrix<T>& K_discont,
                              const gsVector<T>& F_discont,
                              const gsSparseMatrix<T>& T_free,
                              const gsSparseMatrix<T>& T_bnd,
                              const gsVector<T>& u_bnd,
                              gsSparseMatrix<T>& K_free,
                              gsVector<T>& F_free)
    {
        K_free = T_free.transpose() * K_discont * T_free;
        K_free.makeCompressed();
        if (T_bnd.cols() > 0 && u_bnd.size() == T_bnd.cols())
        {
            const gsVector<T> c_bnd = T_bnd * u_bnd;
            F_free = T_free.transpose() * (F_discont - K_discont * c_bnd);
        }
        else
        {
            F_free = T_free.transpose() * F_discont;
        }
    }

    /**
     * @brief Compute exact surface L^2, H^1, covariant H^2, and shell energy errors
     *        of a discontinuous B-spline coefficient vector c_discont against u_ex.
     */
    template <typename ExactVecFn>
    gsKLShellErrors<T> computeErrors(const gsVector<T>& c_discont,
                                     const ExactVecFn& evalExactComp) const
    {
        const size_t nP = m_mp.nPatches();
        const T tMem = m_mat.thickness;
        const T tBend = (m_mat.thickness * m_mat.thickness * m_mat.thickness) / T(12);

        T l2E2 = 0, l2N2 = 0, h1E2 = 0, h1N2 = 0, h2E2 = 0, h2N2 = 0, enE2 = 0, enN2 = 0;

        #pragma omp parallel for reduction(+:l2E2,l2N2,h1E2,h1N2,h2E2,h2N2,enE2,enN2) schedule(dynamic)
        for (size_t p = 0; p < nP; ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
            const gsGeometry<T>& geo = m_mp.patch(p);
            const index_t pOff = m_patchOff[p];

            const gsBSplineBasis<T>& bu = tb.component(0);
            const gsBSplineBasis<T>& bv = tb.component(1);
            const std::vector<T> brkU = bu.knots().breaks();
            const std::vector<T> brkV = bv.knots().breaks();
            const index_t nEu = static_cast<index_t>(brkU.size()) - 1;
            const index_t nEv = static_cast<index_t>(brkV.size()) - 1;

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

            for (index_t ev = 0; ev < nEv; ++ev)
            {
                ruleV.mapTo(brkV[ev], brkV[ev + 1], qv, wv);
                for (index_t eu = 0; eu < nEu; ++eu)
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
                    tb.evalAllDers_into(pts, 2, ders);
                    tb.active_into(pts.col(0), act);
                    const index_t nAct = act.rows();

                    for (index_t q = 0; q < nQ; ++q)
                    {
                        gsSurfaceKinematics3D<T> kin;
                        if (!kin.compute(gD1, gD2, q))
                            continue;
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
            }
        }

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
        const gsTensorBSplineBasis<2, T>& tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[patchIdx]);
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
        gsMultiPatch<T> dispMP;
        for (size_t p = 0; p < m_mp.nPatches(); ++p)
        {
            const gsTensorBSplineBasis<2, T>& tb =
                dynamic_cast<const gsTensorBSplineBasis<2, T>&>(m_dbasis[p]);
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
    const gsMultiPatch<T>& m_mp;
    const gsMultiBasis<T>& m_dbasis;
    gsKLShellMaterial<T> m_mat;
    index_t m_quadExtra;
    std::vector<index_t> m_patchOff;
    index_t m_nDiscontScalar;
};

} // namespace gismo
