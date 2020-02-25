/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsG1Basis/gsGlobalGDAssembler.h>


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
                 index_t patchId_local,
                 bool isBoundary,
                 gsOptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_patchId_local(patchId_local), m_isBoundary(isBoundary), m_optionList(optionList)
    {
        m_gamma = 1.0;

        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");

        m_r = m_optionList.getInt("regularity");

        if (m_optionList.getSwitch("local"))
            gsInfo << "Is not yet implemented \n";
        else
            setGlobalGluingData();

        //beta_exact();

        //conditionTest();
    }

    // Not exact!!! Interpolation via Greville abscissae // But only necessary for condition test!
    void beta_exact();

    // Computed the gluing data globally
    void setGlobalGluingData();

    void conditionTest()
    { /*
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        beta_exact();

        gsMatrix<> temp;

        temp = alpha_tilde_0.eval(points).cwiseProduct(beta_tilde_1.eval(points))
            + alpha_tilde_1.eval(points).cwiseProduct(beta_tilde_0.eval(points))
            - beta_bar.eval(points);

        gsInfo << "Conditiontest gluing data : " << temp.array().abs().maxCoeff() << "\n\n";
*/
    } // Conditiontest

    const gsBSpline<T> get_alpha_tilde_0() const {return alpha_tilde; }
    const gsBSpline<T> get_alpha_tilde_1() const {return alpha_tilde; }
    const gsBSpline<T> get_beta_tilde_0() const {return beta_tilde; }
    const gsBSpline<T> get_beta_tilde_1() const {return beta_tilde; }
    //const gsBSpline<T> get_beta_bar() const {return beta_bar; }

    void eval_into_alpha_0(const gsMatrix<T> & points, gsMatrix<T>& result);
    void eval_into_alpha_1(const gsMatrix<T> & points, gsMatrix<T>& result);
    void eval_into_beta_0(const gsMatrix<T> & points, gsMatrix<T>& result);
    void eval_into_beta_1(const gsMatrix<T> & points, gsMatrix<T>& result);

protected:
    // The geometry for a single interface in the right parametrizations
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    index_t m_patchId_local;
    bool m_isBoundary;
    gsOptionList m_optionList;

    real_t m_gamma;

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    // Regularity of the geometry
    index_t m_r;

protected:
    // Global Gluing data
    //gsBSpline<T> beta_bar;
    gsBSpline<T> alpha_tilde;
    gsBSpline<T> beta_tilde;

}; // class gsGluingData

template <class T>
void gsGluingData<T>::beta_exact()
{ /*
    // Spline space for beta
    index_t m_k, m_p;

    // Assume for now only uniform refinement TODO: Do it for general knotvectors
    gsBasis<T> & basis_0 = m_mb.basis(0).component(1); // patch 0 interface at v-direction
    gsBasis<T> & basis_1 = m_mb.basis(1).component(0); // patch 1 interface at u-direction
    if (basis_0.numElements() >= basis_1.numElements())
        m_k = basis_1.numElements()-1;
    else
        m_k = basis_0.numElements()-1;

    // Need that for beta_exact
    if (basis_0.maxDegree() >= basis_1.maxDegree())
        m_p = basis_1.maxDegree();
    else
        m_p = basis_0.maxDegree();

    // first,last,interior,mult_ends,mult_interior,degree
    gsKnotVector<T> kv(0, 1, m_k, 2 * m_p + 1, 2 * m_p - m_r );
    gsBSplineBasis<T> bsp(kv);

    gsMatrix<> greville = bsp.anchors();
    gsMatrix<> uv0, uv1, ev0, ev1;

    const index_t d = m_mp.parDim();
    gsMatrix<> D0(d,d);

    gsGeometry<>::Ptr beta_temp;

    // Set points for Patch 0
    uv0.setZero(2,greville.cols());
    uv0.bottomRows(1) = greville; // v

    // Set points for Patch 1
    uv1.setZero(2,greville.cols());
    uv1.topRows(1) = greville; // u

    const gsGeometry<> & P0 = m_mp.patch(0);
    const gsGeometry<> & P1 = m_mp.patch(1);
    // ======================================

    // ======== Determine bar{beta} ========
    for(index_t i = 0; i < uv1.cols(); i++)
    {
        P0.jacobian_into(uv0.col(i),ev0);
        P1.jacobian_into(uv1.col(i),ev1);

        D0.col(0) = ev0.col(0); // (DuFL, *)
        D0.col(1) = ev1.col(1); // (*,DuFR)

        uv1(0,i) =  m_gamma * D0.determinant();
    }

    beta_temp = bsp.interpolateData(uv1.topRows(1), greville);
    beta_bar = dynamic_cast<gsBSpline<T> &> (*beta_temp);*/
} // beta_exact

template<class T>
void gsGluingData<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(m_patchId_local).component(m_patchId_local)); // u
    //gsBSplineBasis<> temp_basis_second = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(1).component(1)); // v
/*
    if (temp_basis_first.numElements() >= temp_basis_second.numElements())
    {
        index_t degree = temp_basis_second.maxDegree();
        for (size_t i = degree+1; i < temp_basis_second.knots().size() - (degree+1); i = i+(degree-m_r))
            bsp_gD.insertKnot(temp_basis_second.knot(i),p_tilde-r_tilde);
    }
    else
    {
        index_t degree = temp_basis_first.maxDegree();
        for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i = i+(degree-m_r))
            bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);

    }
*/

    index_t degree = temp_basis_first.maxDegree();
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i = i+(degree-m_r))
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);


    gsGlobalGDAssembler<T> globalGdAssembler(bsp_gD,m_patchId_local,m_mp,m_gamma);
    globalGdAssembler.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver;
    gsVector<> sol_a, sol_b;

    // alpha^S
    solver.compute(globalGdAssembler.matrix_alpha());
    sol_a = solver.solve(globalGdAssembler.rhs_alpha());

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol_a);
    gsBSpline<T> alpha_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    alpha_tilde = alpha_tilde_L_2;


    // beta^S
    solver.compute(globalGdAssembler.matrix_beta());
    sol_b = solver.solve(globalGdAssembler.rhs_beta());

    tilde_temp = bsp_gD.makeGeometry(sol_b);
    gsBSpline<T> beta_tilde_L_2 = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    beta_tilde = beta_tilde_L_2;

    gsWriteParaview(alpha_tilde,"alpha_tilde_L",5000);

    gsWriteParaview(beta_tilde,"beta_tilde_L",5000);

    if (m_optionList.getSwitch("plot"))
    {
        gsWriteParaview(alpha_tilde,"alpha_tilde_L",5000);

        gsWriteParaview(beta_tilde,"beta_tilde_L",5000);
    }
} // setGlobalGluingData


template<class T>
void gsGluingData<T>::eval_into_alpha_0(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{alpha^(L)} == Patch 0 ========
    const gsGeometry<> & PL = m_mp.patch(0); // iFace.second().patch = 0 = Left

    gsMatrix<T> points_L(2,points.cols());
    points_L.setZero();
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
void gsGluingData<T>::eval_into_alpha_1(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{alpha^(R)} == Patch 1 ========
    const gsGeometry<> & PR = m_mp.patch(1); // iFace.first().patch = 1 = Right

    gsMatrix<T> points_R(2,points.cols());
    points_R.setZero();
    points_R.row(1) = points.topRows(1);

    gsMatrix<T> ev1;
    result.resize(1,points.cols());

    for (index_t i = 0; i < points_R.cols(); i++)
    {
        PR.jacobian_into(points_R.col(i), ev1);
        result(0, i) = m_gamma * ev1.determinant(); // erste spalte: alphaL an stelle zweiter spalte
    }

} // eval_into_alpha_R

template<class T>
void gsGluingData<T>::eval_into_beta_0(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{beta}^L ========
    const gsGeometry<> & PL = m_mp.patch(0); // iFace.second().patch = 0 = Left

    gsMatrix<T> points_L(2,points.cols());
    points_L.setZero();
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
void gsGluingData<T>::eval_into_beta_1(const gsMatrix<T> & points, gsMatrix<T> & result)
{
    // ======== Determine bar{beta}^R ========
    const gsGeometry<> & P1 = m_mp.patch(1); // iFace.first().patch = 1 = Right

    gsMatrix<T> points_R(2,points.cols());
    points_R.setZero();
    points_R.row(1) = points.topRows(1);

    gsMatrix<T> ev1, D0;
    result.resize(1,points.cols());

    for(index_t i = 0; i < points_R.cols(); i++)
    {
        P1.jacobian_into(points_R.col(i),ev1);
        D0 = ev1.col(0);
        real_t D1 = 1/ D0.norm();
        points_R(0,i) = m_gamma * D1 * D1 * ev1.col(0).transpose() * ev1.col(1);
    }

} // eval_into_beta_R

} // namespace gismo

