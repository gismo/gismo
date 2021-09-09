/** @file gsG1MultiBasis.h
    @brief Provides G1 Basis.

    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/
#pragma once


namespace gismo
{
template <class T>
class gsG1MultiBasis
{

public:
    gsG1MultiBasis(gsMultiPatch<T> const         & patches,
                   const gsMultiBasis<T>               & bases)
                   : m_mp(patches), m_mb(bases)
    {
        // Create plus minus basis
        basis_pm.clear();
        for (size_t numInt = 0; numInt < m_mp.interfaces().size(); numInt++) // No loop for single patch!!!
        {
            const boundaryInterface &item = m_mp.interfaces()[numInt];

            index_t dir = item.first().m_index < 3 ? 1 : 0; // dir of interface

            // TODO assume that basis 2 is the same space
            gsBSplineBasis<T> basis_int = dynamic_cast<gsBSplineBasis<T>&>(m_mb.basis(item.first().patch).component(dir));

            //gsInfo << "Basis int : " << basis_int.knots().asMatrix() << "\n";

            index_t p = basis_int.degree();
            index_t r = p - basis_int.knots().multiplicityIndex(p+1); // r = p - m

            //gsInfo << "p: " << p << " r: " << r << "\n";

            // first,last,interior,mult_ends,mult_interior
            gsKnotVector<real_t> kv_plus(0, 1, 0, p + 1, p - 1 - r); // p,r+1 //-1 bc r+1
            gsBSplineBasis<T> basis_plus(kv_plus);

            for (size_t i = p + 1; i < basis_int.knots().size() - (p + 1); i += basis_int.knots().multiplicityIndex(i))
                basis_plus.insertKnot(basis_int.knot(i), p - 1 - r);

            basis_pm.push_back(basis_plus);
            //gsInfo << "Basis plus : " << basis_plus.knots().asMatrix() << "\n";

            gsKnotVector<real_t> kv_minus(0, 1, 0, p + 1 - 1, p - 1 - r); // p-1,r //-1 bc p-1
            gsBSplineBasis<T> basis_minus(kv_minus);

            for (size_t i = p + 1; i < basis_int.knots().size() - (p + 1); i += basis_int.knots().multiplicityIndex(i))
                    basis_minus.insertKnot(basis_int.knot(i), p - 1 - r);

            basis_pm.push_back(basis_minus);
            //gsInfo << "Basis minus : " << basis_minus.knots().asMatrix() << "\n";

            // TODO assume that basis 2 is the same space
            basis_geo = dynamic_cast<gsBSplineBasis<T>&>(m_mb.basis(item.first().patch).component(1-dir));
        }
    }

    index_t nPatches() { return m_mp.nPatches(); };
    index_t nBasisFunctions() { return basis_pm.size() > 0 ? basis_pm[0].size() + basis_pm[1].size() : 0; };

    void active_into(const gsMatrix<T> &, gsMatrix<index_t> &, index_t);

    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> >& result, index_t patchIdx);

    void evalAllDers_alpha_S_into(const gsMatrix<T> & points, index_t n, std::vector<gsMatrix<T> >& result, index_t patchIdx, index_t side);
    void evalAllDers_beta_S_into(const gsMatrix<T> & points, index_t n, std::vector<gsMatrix<T> >& result, index_t patchIdx, index_t side);


    void eval_beta_into(const gsMatrix<T> & u, gsMatrix<T>& result);

protected:
    gsMultiPatch<T> const m_mp;
    gsMultiBasis<T> m_mb;

    std::vector<gsBSplineBasis<T>> basis_pm;

    gsBSplineBasis<T> basis_geo;

    std::vector<gsMatrix<index_t>> g1active;
}; // class gsG1MultiBasis

template<class T>
void gsG1MultiBasis<T>::evalAllDers_alpha_S_into(const gsMatrix<T> & points, index_t n,
                                         std::vector<gsMatrix<T> >& result_alpha, index_t patchIdx, index_t side)
{
    result_alpha.resize(n+1);

    GISMO_ASSERT(points.rows() == 2, "Points should be 2D");

    gsMatrix<> uv, ev, ev2;
    if (side==1 || side==3) // West or south
    {
        uv.setZero(2, points.cols());
        uv.row(side == 1 ? 1 : 0) = points.row(side == 1 ? 1 : 0); // West == v, south == u
    }
    else if (side==2 || side==4) // East or north
    {
        uv.setOnes(2, points.cols());
        uv.row(side == 2 ? 1 : 0) = points.row(side == 2 ? 1 : 0); // East == v, north == u
    }

    // ======== Determine bar{alpha^(S)} ========
    const gsGeometry<> & P0 = m_mp.patch(patchIdx);

    result_alpha[0].setZero(1, points.cols());
    P0.deriv_into(uv, ev);
    result_alpha[0] = (patchIdx == 0 ? -1 : 1) * ( ev.row(0).cwiseProduct(ev.row(3)) - ev.row(1).cwiseProduct(ev.row(2)) );

    if (n >= 1)
    {
        result_alpha[1].setZero(1, points.cols());

        std::vector<gsMatrix<>> ders;
        m_mp.patch(patchIdx).basis().evalAllDersFunc_into(uv,m_mp.patch(patchIdx).coefs(),n+1,ders);
        ev = ders[1];
        ev2 = ders[2];

        if (side < 3)
            result_alpha[1] = (patchIdx == 0 ? -1 : 1) * (ev2.row(2).cwiseProduct(ev.row(3)) +
                ev2.row(4).cwiseProduct(ev.row(0)) - ev2.row(1).cwiseProduct(ev.row(2)) -
                ev2.row(5).cwiseProduct(ev.row(1)));
        if (side > 2)
            gsInfo << "not implemented \n";

        if (n >= 2)
        {
            gsMatrix<> ev3;
            ev3 = ders[3];

            if (side < 3)
                result_alpha[2] = (patchIdx == 0 ? -1 : 1) * (-2.0 * ev2.row(5).cwiseProduct(ev2.row(1)) +
                    2.0 * ev2.row(2).cwiseProduct(ev2.row(4)) + ev.row(3).cwiseProduct(ev3.row(2)) -
                    ev.row(1).cwiseProduct(ev3.row(6)) - ev.row(2).cwiseProduct(ev3.row(3)) +
                    ev.row(0).cwiseProduct(ev3.row(7)));
            if (side > 2)
                gsInfo << "not implemented \n";
        }
    }
} // evalAllDers_alpha_S_into


template<class T>
void gsG1MultiBasis<T>::evalAllDers_beta_S_into(const gsMatrix<T> & points, index_t n,
                                                 std::vector<gsMatrix<T> >& result_beta, index_t patchIdx, index_t side)
{
    result_beta.resize(n+1);

    GISMO_ASSERT(points.rows() == 2, "Points should be 2D");

    index_t dir = side < 3 ? 1 : 0;

    gsMatrix<> uv, ev, ev2, ev3;
    if (side==1 || side==3) // West or south
    {
        uv.setZero(2, points.cols());
        uv.row(side == 1 ? 1 : 0) = points.row(side == 1 ? 1 : 0); // West == v, south == u
    }
    else if (side==2 || side==4) // East or north
    {
        uv.setOnes(2, points.cols());
        uv.row(side == 2 ? 1 : 0) = points.row(side == 2 ? 1 : 0); // East == v, north == u
    }

    // ======== Determine bar{alpha^(S)} ========
    const gsGeometry<> & P0 = m_mp.patch(patchIdx);

    const index_t d = m_mp.parDim();
    gsVector<> D0(d);

    result_beta[0].setZero(1, points.cols());
    if (n > 0)
        result_beta[1].setZero(1, points.cols());
    if (n > 1)
        result_beta[2].setZero(1, points.cols());
    for(index_t i = 0; i < uv.cols(); i++)
    {
        P0.jacobian_into(uv.col(i),ev);
        D0 = ev.col(dir);
        real_t D1 = 1/ D0.norm();
        result_beta[0](0,i) = (patchIdx == 0 ? 1 : -1) * D1 * D1 * ev.col(1).transpose() * ev.col(0);

        if (n > 0)
        {
            P0.deriv2_into(uv.col(i), ev2);
            D1 = 1/ D0.squaredNorm();
            real_t D2 = D0.squaredNorm();
            if (dir == 1)
                result_beta[1](0,i) = (patchIdx == 0 ? 1 : -1) * D1 * D1 * (D2*(ev2(2,0)*ev(0,1) + ev2(1,0)*ev(0,0)+
                    ev2(5,0)*ev(1,1) + ev2(4,0)*ev(1,0)) -
                    (ev.col(1).transpose() * ev.col(0))(0,0) * 2.0 * (ev2(1,0)*ev(0,1) + ev2(4,0)*ev(1,1)));
            else if (dir == 0)
                result_beta[1](0,i) = (patchIdx == 0 ? 1 : -1) * D1 * D1 * (D2*(ev2(0,0)*ev(0,1) + ev2(2,0)*ev(0,0)+
                    ev2(3,0)*ev(1,1) + ev2(5,0)*ev(1,0)) -
                    (ev.col(1).transpose() * ev.col(0))(0,0) * 2.0 * (ev2(0,0)*ev(0,0) + ev2(3,0)*ev(1,0)));

            if (n > 1)
            {
                std::vector<gsMatrix<>> ders;
                m_mp.patch(patchIdx).basis().evalAllDersFunc_into(uv.col(i),m_mp.patch(patchIdx).coefs(),3,ders); // TODO before
                ev3 = ders[3];

                if (dir == 1)
                    result_beta[2](0,i) = (patchIdx == 0 ? 1 : -1) * D1 * D1 * D1 * (
                            -4 * D2 * (ev(0,1)*ev2(2,0) + ev(1,1)*ev2(5,0) +
                            ev(0,0)*ev2(1,0) + ev(1,0)*ev2(4,0)) *
                            (ev(0,1)*ev2(1,0) + ev(1,1)*ev2(4,0))
                            +
                            D2 * D2 * (2 * ev2(2,0)*ev2(1,0) + 2 * ev2(5,0)*ev2(4,0) +
                            ev(0,1)*ev3(2,0) + ev(1,1)*ev3(6,0) +
                            ev(0,0)*ev3(3,0) + ev(1,0)*ev3(7,0))
                            +
                            (ev(0,0)*ev(0,1) + ev(1,0)*ev(1,1)) *
                            (8 * (ev(0,1)*ev2(1,0) + ev(1,1)*ev2(4,0)) *
                            (ev(0,1)*ev2(1,0) + ev(1,1)*ev2(4,0)) -
                            2 * (ev(0,1)*ev(0,1) + ev(1,1)*ev(1,1)) *
                            (ev2(1,0)*ev2(1,0) + ev2(4,0)*ev2(4,0) +
                            ev(0,1)*ev3(3,0) + ev(1,1)*ev3(7,0))
                            )
                            );
            }
        }
    }
} // evalAllDers_beta_into


template <class T>
void gsG1MultiBasis<T>::eval_beta_into(const gsMatrix <T> & points, gsMatrix <T> & result)
{
    // TODO just work for 2 patches

    gsMatrix<T> uv0, uv1, ev0, ev1;

    uv0.setOnes(2, points.cols());
    uv0.row(1) = points; // v

    uv1.setZero(2, points.cols());
    uv1.row(1) = points; // v

    const gsGeometry<> &P0 = m_mp.patch(0); // iFace.first().patch = 1
    const gsGeometry<> &P1 = m_mp.patch(1); // iFace.second().patch = 0
    // ======================================

    const index_t d = 2;
    gsMatrix<> D0(d,d);

    // ======== Determine bar{beta} ========
    result.setZero(1, points.cols());
    for (index_t i = 0; i < uv1.cols(); i++) {
        P0.jacobian_into(uv0.col(i), ev0);
        P1.jacobian_into(uv1.col(i), ev1);

        D0.col(0) = ev0.col(0); // (DuFL, *)
        D0.col(1) = ev1.col(0); // (*,DuFR)

        result(0, i) = D0.determinant();
    }
} // eval_beta_into

template <class T>
void gsG1MultiBasis<T>::active_into(const gsMatrix <T> & points, gsMatrix<index_t> & active, index_t patchIdx)
{
    for (size_t numInt = 0; numInt < m_mp.interfaces().size(); numInt++) // No loop for single patch!!!
    {
        const boundaryInterface &item = m_mp.interfaces()[numInt];

        index_t dir = -1, idx_geo = -1;
        if (item.first().patch == patchIdx)
        {
            dir = item.first().m_index < 3 ? 1 : 0; // dir of interface
            idx_geo = item.first().m_index % 2 > 0 ? 1 : basis_geo.size() - 2; // 1 or n-1
        }
        else if (item.second().patch == patchIdx)
        {
            dir = item.second().m_index < 3 ? 1 : 0; // dir of interface
            idx_geo = item.second().m_index % 2 > 0 ? 1 : basis_geo.size() - 2; // 1 or n-1
        }

        gsMatrix<index_t> act_temp, act_plus, act_minus;
        basis_geo.active_into(points.row(1-dir), act_temp);

        for (index_t act_i = 0; act_i < act_temp.rows(); act_i++)
            if (act_temp(act_i,0) == idx_geo) // points is at the interface
            {
                basis_pm[0].active_into(points.row(dir), act_plus);
                basis_pm[1].active_into(points.row(dir), act_minus);
            }

        // For computing the basis function
        g1active.clear();
        g1active.push_back(act_plus);
        g1active.push_back(act_minus);

        // for getting the local dofs
        gsMatrix<index_t>ones(act_plus.rows(), act_plus.cols());
        ones.setOnes();
        active = act_plus + m_mb.basis(patchIdx).size() * ones;
        active.conservativeResize(active.rows() + act_minus.rows(), active.cols());

        ones.setOnes(act_minus.rows(), act_minus.cols());
        active.bottomRows(act_minus.rows()) = act_minus + (basis_pm[0].size() + m_mb.basis(patchIdx).size()) * ones;

    }
} // active_into

template<class T>
void gsG1MultiBasis<T>::evalAllDers_into(const gsMatrix<T> & points, index_t n, std::vector<gsMatrix<T>> & result, index_t patchIdx)
{
/*
    index_t p_size = 5;
    gsVector<> vec;
    vec.setLinSpaced(p_size,0,1);

*/

    if (g1active.size() == 0 )
    {
        gsMatrix<index_t> act_plus, act_minus;
        gsVector<index_t> vec;
        vec.setLinSpaced(basis_pm[0].size(),0,basis_pm[0].size()-1);
        act_plus = vec;
        vec.setLinSpaced(basis_pm[1].size(), 0, basis_pm[1].size()-1);
        act_minus = vec;

        g1active.clear();
        g1active.push_back(act_plus);
        g1active.push_back(act_minus);
    }

    result.resize(n+1);
    result[0].setZero(g1active[0].rows()+g1active[1].rows(), points.cols());
    if (n > 0)
        result[1].setZero(2*(g1active[0].rows()+g1active[1].rows()), points.cols());
    if (n > 1)
        result[2].setZero(3*(g1active[0].rows()+g1active[1].rows()), points.cols());

    const boundaryInterface &item = m_mp.interfaces()[0];

    index_t dir = -1, idx_geo = -1;
    if (item.first().patch == patchIdx)
    {
        dir = item.first().m_index < 3 ? 1 : 0; // dir of interface
        idx_geo = item.first().m_index % 2 > 0 ? 1 : basis_geo.size() - 2; // 1 or n-1
    }
    else if (item.second().patch == patchIdx)
    {
        dir = item.second().m_index < 3 ? 1 : 0; // dir of interface
        idx_geo = item.second().m_index % 2 > 0 ? 1 : basis_geo.size() - 2; // 1 or n-1
    }


    // tau/p
    real_t p = basis_geo.degree();
    real_t tau_1 = basis_geo.knots().at(p + 1); // p + 2 TODO Assume that mesh size is the same!!!

    gsMatrix<T>
            N_0, N_1,
            N_j_minus, N_i_plus,
            der_N_i_plus;

    // For the first derivative
    gsMatrix<T>
            der_N_0, der_N_1,
            der_N_j_minus,
            der2_N_i_plus;

    // For the second derivative
    gsMatrix<T>
            der2_N_0, der2_N_1,
            der2_N_j_minus,
            der3_N_i_plus;


    // Compute gluing data
    index_t side = item.first().patch == patchIdx ? item.first().m_index : item.second().m_index;
    std::vector<gsMatrix<T>> result_alpha, result_beta;

    evalAllDers_alpha_S_into(points, n, result_alpha, patchIdx, side);
    evalAllDers_beta_S_into(points, n, result_beta, patchIdx, side);

    // ======== For modifying beta ========
    gsMatrix<T> zeroOne(2,2);
    zeroOne.setIdentity();

    // For modifying beta
    gsMatrix<T> alpha2, beta2;
    alpha2.setZero(1, zeroOne.cols());
    beta2.setZero(1, zeroOne.cols());

    const gsGeometry<> & PL = m_mp.patch(0); // Only in two Patch case TODO
    const gsGeometry<> & PR = m_mp.patch(1); // Only in two Patch case TODO

    gsMatrix<T> evZeroOne;
    PL.jacobian_into(zeroOne.col(0), evZeroOne);
    alpha2(0,0) = -1 * evZeroOne.determinant(); // alpha^L (0)

    const index_t d = m_mp.parDim();
    gsVector<> D0(d);
    D0 = evZeroOne.col(1); // dir of interface == 1
    real_t D1 = 1/ D0.norm();
    beta2(0,0) = D1 * D1 * evZeroOne.col(1).transpose() * evZeroOne.col(0); // beta^L (0)


    PR.jacobian_into(zeroOne.col(1), evZeroOne);
    alpha2(0,1) = evZeroOne.determinant(); // alpha^R (1)

    D0 = evZeroOne.col(1); // dir of interface == 1
    D1 = 1/ D0.norm();
    beta2(0,1) = -1 * D1 * D1 * evZeroOne.col(1).transpose() * evZeroOne.col(0); // beta^R (1)

    real_t lambdaL = beta2(0,0)/alpha2(0,0);
    real_t lambdaR = beta2(0,1)/alpha2(0,1);

    gsMatrix<> ones;
    ones.setOnes(result_beta[0].rows(), result_beta[0].cols());
    result_beta[0] = result_beta[0] - lambdaL*(ones - points.row(dir)).cwiseProduct(result_alpha[0]) - lambdaR*(points.row(dir)).cwiseProduct(result_alpha[0]);

    if (n > 0)
        result_beta[1] = result_beta[1] - lambdaL*(ones - points.row(dir)).cwiseProduct(result_alpha[1]) - lambdaR*(points.row(dir)).cwiseProduct(result_alpha[1])
            + lambdaL*result_alpha[0] - lambdaR*result_alpha[0];
    if (n > 1)
        result_beta[2] = result_beta[2] - lambdaL*((ones - points.row(dir)).cwiseProduct(result_alpha[2]) - 2.0*result_alpha[1]) -
            lambdaR*((points.row(dir)).cwiseProduct(result_alpha[2]) + 2.0*result_alpha[1]);

    // End compute gluing data



    basis_geo.evalSingle_into(idx_geo == 1 ? 0 : idx_geo + 1, points.row(1-dir),N_0);
    basis_geo.evalSingle_into(idx_geo,points.row(1-dir),N_1);

    basis_geo.derivSingle_into(idx_geo == 1 ? 0 : idx_geo + 1, points.row(1-dir),der_N_0);
    basis_geo.derivSingle_into(idx_geo, points.row(1-dir),der_N_1);

    basis_geo.deriv2Single_into(idx_geo == 1 ? 0 : idx_geo + 1, points.row(1-dir),der2_N_0); // u
    basis_geo.deriv2Single_into(idx_geo, points.row(1-dir),der2_N_1); // u

    // TODO Make matrix multiplication instead of a for loop
    // Compute plus basis functions
    for (index_t i_plus = 0; i_plus < g1active[0].rows(); i_plus++)
    {
        index_t bfID = g1active[0](i_plus,0);

        basis_pm[0].evalSingle_into(bfID, points.row(dir), N_i_plus);
        basis_pm[0].derivSingle_into(bfID, points.row(dir), der_N_i_plus);

        gsMatrix<T> temp = result_beta[0].cwiseProduct(der_N_i_plus);

        result[0].row(i_plus) = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;
        if (n > 0)
        {
            basis_pm[0].deriv2Single_into(bfID, points.row(dir),der2_N_i_plus);

            gsMatrix<T> der_temp = result_beta[1].cwiseProduct(der_N_i_plus) + result_beta[0].cwiseProduct(der2_N_i_plus);

            result[1].row(2*i_plus + 1-dir) = N_i_plus.cwiseProduct(der_N_0 + der_N_1) - temp.cwiseProduct(der_N_1) * tau_1 / p;
            result[1].row(2*i_plus + dir) = der_N_i_plus.cwiseProduct(N_0 + N_1) - der_temp.cwiseProduct(N_1) * tau_1 / p;

            if (n > 1)
            {
                basis_pm[0].evalDerSingle_into(bfID, points.row(dir),3,der3_N_i_plus);

                gsMatrix<T> der2_temp = 2.0*result_beta[1].cwiseProduct(der2_N_i_plus) + result_beta[0].cwiseProduct(der3_N_i_plus) + result_beta[2].cwiseProduct(der_N_i_plus);

                result[2].row(3*i_plus + 1-dir) = N_i_plus.cwiseProduct(der2_N_0 + der2_N_1) - temp.cwiseProduct(der2_N_1) * tau_1 / p;
                result[2].row(3*i_plus + dir) = der2_N_i_plus.cwiseProduct(N_0 + N_1) - der2_temp.cwiseProduct(N_1) * tau_1 / p;
                result[2].row(3*i_plus + 2) = der_N_i_plus.cwiseProduct(der_N_0 + der_N_1) - der_temp.cwiseProduct(der_N_1) * tau_1 / p;
            }
        }

    }
    //  End plus basis

    // Compute minus basis
    for (index_t i_minus = 0; i_minus < g1active[1].rows(); i_minus++)
    {
        index_t bfID = g1active[1](i_minus,0);

        basis_pm[1].evalSingle_into(bfID, points.row(dir), N_j_minus);

        result[0].row(g1active[0].rows() + i_minus) = result_alpha[0].cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p; // (dir == 0 ? -1 : 1)
        if (n > 0)
        {
            basis_pm[1].derivSingle_into(bfID, points.row(dir),der_N_j_minus);

            result[1].row(2*g1active[0].rows() + 2*i_minus + 1-dir) = result_alpha[0].cwiseProduct(N_j_minus.cwiseProduct(der_N_1)) * tau_1 / p;
            result[1].row(2*g1active[0].rows() + 2*i_minus + dir) = (result_alpha[1].cwiseProduct(N_j_minus)+result_alpha[0].cwiseProduct(der_N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
        }
        if (n > 1)
        {
            basis_pm[1].deriv2Single_into(bfID, points.row(dir), der2_N_j_minus);

            result[2].row(3*g1active[0].rows() + 3*i_minus + 1-dir) = result_alpha[0].cwiseProduct(N_j_minus.cwiseProduct(der2_N_1)) * tau_1 / p;
            result[2].row(3*g1active[0].rows() + 3*i_minus + dir) = (2.0*result_alpha[1].cwiseProduct(der_N_j_minus)+result_alpha[0].cwiseProduct(der2_N_j_minus)+result_alpha[2].cwiseProduct(N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
            result[2].row(3*g1active[0].rows() + 3*i_minus + 2) = (result_alpha[1].cwiseProduct(N_j_minus)+result_alpha[0].cwiseProduct(der_N_j_minus)).cwiseProduct(der_N_1) * tau_1 / p;
        }
    }
    // End minus basis
} // deriv2_into


} // namespace gismo

