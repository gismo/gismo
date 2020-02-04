/** @file gsGluingData_mp.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsG1Basis/gsGlobalGDAssembler_mp.h>




namespace gismo
{

template<class T>
class gsGluingData_mp
{
public:
    gsGluingData_mp()
    { }

    gsGluingData_mp(gsMultiPatch<T> const & mp,
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

        beta_exact(); // Not exact!!! Interpolation via Greville abscissae // But only necessary for condition test!

        // LoKALE PROJECTION PROGRAMMIEREN
        setGlobalGluingData();

        //conditionTest();
    }

    void beta_exact();
    void setGlobalGluingData();

    boundaryInterface get_iFace() { return  m_iFace;}

    void conditionTest()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        gsMatrix<> temp;

        temp = alpha_tilde_L.eval(points).cwiseProduct(beta_tilde_R.eval(points))
            - alpha_tilde_R.eval(points).cwiseProduct(beta_tilde_L.eval(points))
            + beta_bar.eval(points);

        gsInfo << "\nConditiontest gluing data for patches " << m_iFace.first() << " and " << m_iFace.second() <<
        ":\n" << temp.array().abs().maxCoeff() << "\n\n";

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

}; // class gsGluingData_mp

template <class T>
void gsGluingData_mp<T>::beta_exact()
{
    // first,last,interior,mult_ends,mult_interior,degree
    gsKnotVector<T> kv(0, 1, m_k, 2 * m_p + 1, 2 * m_p - m_r );
    gsBSplineBasis<T> bsp(kv);

    gsMatrix<> greville = bsp.anchors();
    gsMatrix<> uv1, uv2, ev1, ev2;

    const index_t d = m_mp.parDim();
    gsMatrix<> D0(d,d);

    gsGeometry<>::Ptr beta_temp;

    index_t nulleins;
    real_t plusminus;
    if (m_iFace.first().side().index() >= 3) // North and south
    {
        if (m_iFace.first().side().index() == 4) // North
            uv1.setOnes(2,greville.cols());
        if (m_iFace.first().side().index() == 3) // South
            uv1.setZero(2,greville.cols());
        uv1.topRows(1) = greville;
        nulleins = 1;
        plusminus = -1;
    }
    else if (m_iFace.first().side().index() <= 2) // East and west
    {

        if (m_iFace.first().side().index() == 2) // East
            uv1.setOnes(2, greville.cols());
        if (m_iFace.first().side().index() == 1) // West
            uv1.setZero(2, greville.cols());
        uv1.bottomRows(1) = greville;
        nulleins = 0;
        plusminus = 1;
    }
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

        D0.col(0) = ev1.col(nulleins); // (DuFL, *)
        D0.col(1) = ev2.col(nulleins); // (*,DuFR)

        uv1(0,i) = plusminus * m_gamma * D0.determinant();
    }

    beta_temp = bsp.interpolateData(uv1.topRows(1), greville);
    beta_bar = dynamic_cast<gsBSpline<T> &> (*beta_temp);
}

template<class T>
void gsGluingData_mp<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    //std::vector<gsBSplineBasis<T>> bsp_gD;
    gsBSplineBasis<T> bsp_gD(kv);

    gsTensorBSplineBasis<2,real_t> & temp_Tensorbasis_L = dynamic_cast<gsTensorBSplineBasis<2,real_t> &>(m_mb.basis(0));
    gsBSplineBasis<> temp_basis_L = dynamic_cast<gsBSplineBasis<> &>(temp_Tensorbasis_L.component(1));
    gsTensorBSplineBasis<2,real_t> & temp_Tensorbasis_R = dynamic_cast<gsTensorBSplineBasis<2,real_t> &>(m_mb.basis(1));
    gsBSplineBasis<> temp_basis_R = dynamic_cast<gsBSplineBasis<> &>(temp_Tensorbasis_R.component(1));

    index_t max_degree = temp_basis_L.maxDegree(); // The same as right
    for (index_t i = max_degree+1; i < temp_basis_R.knots().size() - (max_degree+1); i = i+(max_degree-m_r))
    {
        //bsp_gD.insertKnot(temp_basis_L.knot(i),p_tilde-r_tilde);
        bsp_gD.insertKnot(temp_basis_R.knot(i),p_tilde-r_tilde);
    }


    gsGlobalGDAssembler_mp<T> globalGdAssembler(bsp_gD,m_mp,m_iFace,m_gamma);
    globalGdAssembler.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver;
    gsVector<> sol_a_L, sol_a_R, sol_b_L, sol_b_R;

    // alpha^S
    solver.compute(globalGdAssembler.matrix_alpha_L());
    sol_a_L = solver.solve(globalGdAssembler.rhs_alpha_L());

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol_a_L);
    gsBSpline<T> alpha_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    alpha_tilde_L = alpha_tilde_L_2;


    // beta^S
    solver.compute(globalGdAssembler.matrix_beta_L());
    sol_b_L = solver.solve(globalGdAssembler.rhs_beta_L());

    tilde_temp = bsp_gD.makeGeometry(sol_b_L);
    gsBSpline<T> beta_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    beta_tilde_L = beta_tilde_L_2;







/*
    gsKnotVector<T> kv_R(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    //std::vector<gsBSplineBasis<T>> bsp_gD;
    gsBSplineBasis<T> bsp_gD_R(kv_R);

    for (index_t i = max_degree+1; i < temp_basis_L.knots().size() - (max_degree+1); i = i+(max_degree-m_r))
    {
        bsp_gD.insertKnot(temp_basis_L.knot(i),p_tilde-r_tilde);
        //bsp_gD_R.insertKnot(temp_basis_R.knot(i),p_tilde-r_tilde);
    }


    gsGlobalGDAssembler_mp<T> globalGdAssembler_R(bsp_gD_R,m_mp,m_iFace,m_gamma);
    globalGdAssembler_R.assemble();
*/
    solver.compute(globalGdAssembler.matrix_alpha_R());
    sol_a_R = solver.solve(globalGdAssembler.rhs_alpha_R());

    solver.compute(globalGdAssembler.matrix_beta_R());
    sol_b_R = solver.solve(globalGdAssembler.rhs_beta_R());

    tilde_temp = bsp_gD.makeGeometry(sol_b_R);
    gsBSpline<T> beta_tilde_R_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    beta_tilde_R = beta_tilde_R_2;

    tilde_temp = bsp_gD.makeGeometry(sol_a_R);
    gsBSpline<T> alpha_tilde_R_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    alpha_tilde_R = alpha_tilde_R_2;

    if (plot)
    {
        gsWriteParaview(alpha_tilde_L,"alpha_tilde_L",5000);
        gsWriteParaview(alpha_tilde_R,"alpha_tilde_R",5000);

        gsWriteParaview(beta_tilde_L,"beta_tilde_L",5000);
        gsWriteParaview(beta_tilde_R,"beta_tilde_R",5000);
    }
}

template<class T>
void gsGluingData_mp<T>::eval_into_alpha_L(const gsMatrix<T> & points, gsMatrix<T> & result)
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
void gsGluingData_mp<T>::eval_into_alpha_R(const gsMatrix<T> & points, gsMatrix<T> & result)
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
void gsGluingData_mp<T>::eval_into_beta_L(const gsMatrix<T> & points, gsMatrix<T> & result)
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
void gsGluingData_mp<T>::eval_into_beta_R(const gsMatrix<T> & points, gsMatrix<T> & result)
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

