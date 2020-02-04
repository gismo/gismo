/** @file gsGluingData.h

    @brief Compute the gluing data for the interfaces.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsNurbs/gsBSpline.h>

#include <gsG1Basis/gsGlobalGDAssembler.h>
#include <gsG1Basis/gsLocalGDAssembler.h>



namespace gismo
{

template<class T>
class gsGluingData
{
public:
    gsGluingData()
    { }

    gsGluingData(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> const & mb,
                 boundaryInterface const & iFace,
                 index_t const interior,
                 index_t const reg,
                 index_t const p_t,
                 index_t const r_t,
                 bool const dir,
                 bool const pl)
                 : m_mp(mp), m_mb(mb), m_iFace(iFace), m_p(mb.maxCwiseDegree()), m_k(interior),
                   m_r(reg), p_tilde(p_t), r_tilde(r_t), direct(dir), plot(pl)
    {

        m_gamma = 1.0;

        beta_exact();
        setGlobalGluingData();

        setLocalGluingData();


    }


    void setGlobalGluingData();
    void setLocalGluingData();

    void beta_exact();

    void conditionTest()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        gsMatrix<> temp;
        if (direct)
        {
            gsMatrix<T> a_L, b_L, a_R, b_R;
            gsMatrix<T> points_L, points_R;
            points_L.setOnes(2,1000);
            points_R.setZero(2,1000);
            points_L.row(1) = points;
            points_R.row(1) = points;

            eval_into_alpha_L(points_L,a_L);
            eval_into_alpha_R(points_R,a_R);
            eval_into_beta_L(points_L,b_L);
            eval_into_beta_R(points_R,b_R);

            temp = beta_tilde_L.eval(points);

            //gsInfo << "\nbeta gluing data: \n" << temp.array().abs().maxCoeff() << "\n\n";
        }
        else
        {
            temp = alpha_tilde_L.eval(points).cwiseProduct(beta_tilde_R.eval(points))
                - alpha_tilde_R.eval(points).cwiseProduct(beta_tilde_L.eval(points))
                + beta_bar.eval(points);
        }

        /*
        temp = alpha_j_tilde_L.at(1).eval(points).cwiseProduct(beta_i_tilde_R.at(1).eval(points))
            - alpha_j_tilde_R.at(1).eval(points).cwiseProduct(beta_i_tilde_L.at(1).eval(points))
            + beta_bar.eval(points);
        */
        gsInfo << "\nConditiontest gluing data: \n" << temp.array().abs().maxCoeff() << "\n\n";

    }

    void eval_into_alpha_L(const gsMatrix<T> & points, gsMatrix<T>& result);
    void eval_into_alpha_R(const gsMatrix<T> & points, gsMatrix<T>& result);
    void eval_into_beta_L(const gsMatrix<T> & points, gsMatrix<T>& result);
    void eval_into_beta_R(const gsMatrix<T> & points, gsMatrix<T>& result);

    const gsBSpline<T> get_alpha_tilde_L() const {return alpha_tilde_L; }
    const gsBSpline<T> get_alpha_tilde_R() const {return alpha_tilde_R; }
    const gsBSpline<T> get_beta_tilde_L() const {return beta_tilde_L; }
    const gsBSpline<T> get_beta_tilde_R() const {return beta_tilde_R; }
    const gsBSpline<T> get_beta_bar() const {return beta_bar; }

    const gsBSpline<T> get_beta_i_tilde_L(index_t i) const {return beta_i_tilde_L.at(i); }
    const gsBSpline<T> get_beta_i_tilde_R(index_t i) const {return beta_i_tilde_R.at(i); }
    const gsBSpline<T> get_alpha_j_tilde_L(index_t i) const {return alpha_j_tilde_L.at(i); }
    const gsBSpline<T> get_alpha_j_tilde_R(index_t i) const {return alpha_j_tilde_R.at(i); }


protected:
    // Single Interface and their geometry
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    boundaryInterface m_iFace;

    // Spline spaces (p,k,r)
    index_t m_p, m_k, m_r;
    index_t p_tilde, r_tilde;

    // Gamma
    real_t m_gamma;

    bool direct, plot;

protected:
    // Global Gluing data
    gsBSpline<T> beta_bar;
    gsBSpline<T> alpha_tilde_L, alpha_tilde_R;
    gsBSpline<T> beta_tilde_L, beta_tilde_R;

    // Local Gluing data
    std::vector<gsBSpline<T>> alpha_j_tilde_L, alpha_j_tilde_R;
    std::vector<gsBSpline<T>> beta_i_tilde_L, beta_i_tilde_R;

}; // class gsGluingData


template<class T>
void gsGluingData<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,m_k,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsGlobalGDAssembler<T> globalGdAssembler(bsp_gD,m_mp,m_iFace,m_gamma);
    globalGdAssembler.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver;
    gsVector<> sol_a_L, sol_a_R, sol_b_L, sol_b_R;

    // alpha^S
    solver.compute(globalGdAssembler.matrix_alpha_L());
    sol_a_L = solver.solve(globalGdAssembler.rhs_alpha_L());

    solver.compute(globalGdAssembler.matrix_alpha_R());
    sol_a_R = solver.solve(globalGdAssembler.rhs_alpha_R());

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol_a_L);
    gsBSpline<T> alpha_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    alpha_tilde_L = alpha_tilde_L_2;

    tilde_temp = bsp_gD.makeGeometry(sol_a_R);
    gsBSpline<T> alpha_tilde_R_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    alpha_tilde_R = alpha_tilde_R_2;

    // beta^S
    solver.compute(globalGdAssembler.matrix_beta_L());
    sol_b_L = solver.solve(globalGdAssembler.rhs_beta_L());

    solver.compute(globalGdAssembler.matrix_beta_R());
    sol_b_R = solver.solve(globalGdAssembler.rhs_beta_R());

    tilde_temp = bsp_gD.makeGeometry(sol_b_L);
    gsBSpline<T> beta_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    beta_tilde_L = beta_tilde_L_2;

    tilde_temp = bsp_gD.makeGeometry(sol_b_R);
    gsBSpline<T> beta_tilde_R_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    beta_tilde_R = beta_tilde_R_2;

    //conditionTest();

    if (plot)
    {
        gsWriteParaview(alpha_tilde_L,"alpha_tilde_L",5000);
        gsWriteParaview(alpha_tilde_R,"alpha_tilde_R",5000);

        gsWriteParaview(beta_tilde_L,"beta_tilde_L",5000);
        gsWriteParaview(beta_tilde_R,"beta_tilde_R",5000);
    }
}

template<class T>
void gsGluingData<T>::setLocalGluingData()
{
    gsKnotVector<T> kv_tilde(0,1,m_k,m_mb.basis(0).component(1).maxDegree()+1,m_mb.basis(0).component(1).maxDegree()-1-m_r); // p,r+1 //-1 bc r+1
    gsBSplineBasis<> basis_tilde(kv_tilde);
    index_t n_tilde = basis_tilde.size();

    gsKnotVector<T> kv_bar(0,1,m_k,m_mb.basis(0).component(1).maxDegree()+1-1,m_mb.basis(0).component(1).maxDegree()-1-m_r); // p-1,r //-1 bc p-1
    gsBSplineBasis<> basis_bar(kv_bar);
    index_t n_bar = basis_bar.size();

    // Setting the space for each alpha_tilde, beta_tilde
    alpha_j_tilde_L.resize(n_bar);
    alpha_j_tilde_R.resize(n_bar);
    beta_i_tilde_L.resize(n_tilde);
    beta_i_tilde_R.resize(n_tilde);

    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,m_k,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    for (index_t i = 0; i < n_bar; i++)
    {
        gsMatrix<T> ab = basis_bar.support(i);

        gsKnotVector<T> kv(ab.at(0), ab.at(1),0, p_tilde+1);

        real_t span = bsp_gD.getMaxCellLength();
        real_t temp_knot = ab.at(0) + span;
        while (temp_knot < ab.at(1))
        {
            kv.insert(temp_knot,p_tilde-r_tilde);
            temp_knot += span;
        }
        gsBSplineBasis<T> bsp_geo(kv);

        gsLocalGDAssembler<T> localGdAssembler(bsp_gD,bsp_geo,m_mp, m_iFace, m_gamma,"alpha");
        localGdAssembler.assemble();

        gsSparseSolver<real_t>::CGDiagonal solver;
        gsVector<> sol_L, sol_R;

        // alpha^S
        solver.compute(localGdAssembler.matrix_L());
        sol_L = solver.solve(localGdAssembler.rhs_L());

        solver.compute(localGdAssembler.matrix_R());
        sol_R = solver.solve(localGdAssembler.rhs_R());

        gsGeometry<>::uPtr tilde_temp;
        tilde_temp = bsp_gD.makeGeometry(sol_L);
        gsBSpline<T> a_t_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        alpha_j_tilde_L.at(i) = a_t_L;

        tilde_temp = bsp_gD.makeGeometry(sol_R);
        gsBSpline<T> a_t_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        alpha_j_tilde_R.at(i) = a_t_R;
    }

    for (index_t i = 0; i < n_tilde; i++)
    {
        gsMatrix<T> ab = basis_tilde.support(i);

        gsKnotVector<T> kv(ab.at(0), ab.at(1),0, p_tilde+1);

        real_t span = bsp_gD.getMaxCellLength();
        real_t temp_knot = ab.at(0) + span;
        while (temp_knot < ab.at(1))
        {
            kv.insert(temp_knot,p_tilde-r_tilde);
            temp_knot += span;
        }
        gsBSplineBasis<T> bsp_geo(kv);

        gsLocalGDAssembler<T> localGdAssembler(bsp_gD,bsp_geo,m_mp, m_iFace, m_gamma,"beta");
        localGdAssembler.assemble();

        gsSparseSolver<real_t>::CGDiagonal solver;
        gsVector<> sol_L, sol_R;

        // beta^S
        solver.compute(localGdAssembler.matrix_L());
        sol_L = solver.solve(localGdAssembler.rhs_L());

        solver.compute(localGdAssembler.matrix_R());
        sol_R = solver.solve(localGdAssembler.rhs_R());

        gsGeometry<>::uPtr tilde_temp;
        tilde_temp = bsp_gD.makeGeometry(sol_L);
        gsBSpline<T> b_t_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        beta_i_tilde_L.at(i) = b_t_L;

        tilde_temp = bsp_gD.makeGeometry(sol_R);
        gsBSpline<T> b_t_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        beta_i_tilde_R.at(i) = b_t_R;
    }

    //conditionTest();

} // setLocalGluingData

template <class T>
void gsGluingData<T>::beta_exact()
{
    // first,last,interior,mult_ends,mult_interior,degree
    gsKnotVector<T> kv(0, 1, m_k, 2 * m_p + 1, 2 * m_p - m_r );
    gsBSplineBasis<T> bsp(kv);

    gsMatrix<> greville = bsp.anchors();
    gsMatrix<> uv1, uv2, ev1, ev2;

    const index_t d = m_mp.parDim();
    gsMatrix<> D0(d,d);

    gsGeometry<>::Ptr beta_temp;

    uv1.setZero(2,greville.cols());
    uv1.bottomRows(1) = greville;

    gsAffineFunction<T> ifaceMap(m_mp.getMapForInterface(m_iFace));
    ifaceMap.eval_into(uv1,uv2); // The correspond point on other patch u2

    const gsGeometry<> & P1 = m_mp.patch(m_iFace.first().patch); // iFace.first().patch = 1
    const gsGeometry<> & P2 = m_mp.patch(m_iFace.second().patch); // iFace.second().patch = 0
    // ======================================

    // ======== Determine bar{beta} ========
    for(index_t i = 0; i < uv1.cols(); i++)
    {
        P1.jacobian_into(uv1.col(i),ev1);
        P2.jacobian_into(uv2.col(i),ev2);

        D0.col(0) = ev1.col(0); // (DuFL, *)
        D0.col(1) = ev2.col(0); // (*,DuFR)

        uv1(0,i) = m_gamma * D0.determinant();
    }

    beta_temp = bsp.interpolateData(uv1.topRows(1), uv1.bottomRows(1));
    beta_bar = dynamic_cast<gsBSpline<T> &> (*beta_temp);
}

template<class T>
void gsGluingData<T>::eval_into_alpha_L(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{alpha^(L)} == Patch 0 ========
    const gsGeometry<> & PL = m_mp.patch(0); // iFace.second().patch = 0 = Left

    gsMatrix<T> points_L(2,points.cols());
    points_L.setOnes();
    points_L.row(1) = points.bottomRows(1);

    gsMatrix<T> ev2;
    result.resize(1,points.cols());

    for (index_t i = 0; i < points_L.cols(); i++)
    {
        PL.jacobian_into(points_L.col(i), ev2);
        result(0, i) = - m_gamma * ev2.determinant();
    }

} // eval_into_alpha_L

template<class T>
void gsGluingData<T>::eval_into_alpha_R(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{alpha^(R)} == Patch 1 ========
    const gsGeometry<> & PR = m_mp.patch(1); // iFace.first().patch = 1 = Right

    gsMatrix<T> points_R(2,points.cols());
    points_R.setZero();
    points_R.row(1) = points.bottomRows(1);

    gsMatrix<T> ev1;
    result.resize(1,points.cols());

    for (index_t i = 0; i < points_R.cols(); i++)
    {
        PR.jacobian_into(points_R.col(i), ev1);
        result(0, i) = m_gamma * ev1.determinant(); // erste spalte: alphaL an stelle zweiter spalte
    }

} // eval_into_alpha_R

template<class T>
void gsGluingData<T>::eval_into_beta_L(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{beta}^L ========
    const gsGeometry<> & PL = m_mp.patch(0); // iFace.second().patch = 0 = Left

    gsMatrix<T> points_L(2,points.cols());
    points_L.setOnes();
    points_L.row(1) = points.bottomRows(1);

    gsMatrix<T> ev2, D0;
    result.resize(1,points.cols());

    for(index_t i = 0; i < points_L.cols(); i++)
    {
        PL.jacobian_into(points_L.col(i),ev2);
        D0 = ev2.col(1);
        real_t D1 = 1/ D0.norm();
        result(0,i) = m_gamma * D1 * D1 * ev2.col(0).transpose() * ev2.col(1);
    }

} // eval_into_beta_L

template<class T>
void gsGluingData<T>::eval_into_beta_R(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{beta}^R ========
    const gsGeometry<> & PR = m_mp.patch(1); // iFace.first().patch = 1 = Right

    gsMatrix<T> points_R(2,points.cols());
    points_R.setZero();
    points_R.row(1) = points.bottomRows(1);

    gsMatrix<T> ev1, D0;
    result.resize(1,points.cols());

    for(index_t i = 0; i < points_R.cols(); i++)
    {
        PR.jacobian_into(points_R.col(i),ev1);
        D0 = ev1.col(1);
        real_t D1 = 1/ D0.norm();
        result(0,i) = - m_gamma * D1 * D1 * ev1.col(0).transpose() * ev1.col(1);
    }

} // eval_into_beta_R

} // namespace gismo

