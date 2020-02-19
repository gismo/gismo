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
                 gsOptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        m_gamma = 1.0;

        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");

        m_r = m_optionList.getInt("regularity");

        if (m_optionList.getSwitch("local"))
            gsInfo << "Is not yet implemented \n";
        else
            setGlobalGluingData();

        beta_exact();

        conditionTest();
    }

    // Not exact!!! Interpolation via Greville abscissae // But only necessary for condition test!
    void beta_exact();

    // Computed the gluing data globally
    void setGlobalGluingData();

    void conditionTest()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        beta_exact();

        gsMatrix<> temp;

        temp = alpha_tilde_L.eval(points).cwiseProduct(beta_tilde_R.eval(points))
            + alpha_tilde_R.eval(points).cwiseProduct(beta_tilde_L.eval(points))
            - beta_bar.eval(points);

        gsInfo << "Conditiontest gluing data : " << temp.array().abs().maxCoeff() << "\n\n";

    } // Conditiontest

    const gsBSpline<T> get_alpha_tilde_L() const {return alpha_tilde_L; }
    const gsBSpline<T> get_alpha_tilde_R() const {return alpha_tilde_R; }
    const gsBSpline<T> get_beta_tilde_L() const {return beta_tilde_L; }
    const gsBSpline<T> get_beta_tilde_R() const {return beta_tilde_R; }
    const gsBSpline<T> get_beta_bar() const {return beta_bar; }
protected:
    // The geometry for a single interface in the right parametrizations
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsOptionList m_optionList;

    real_t m_gamma;

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    // Regularity of the geometry
    index_t m_r;

protected:
    // Global Gluing data
    gsBSpline<T> beta_bar;
    gsBSpline<T> alpha_tilde_L, alpha_tilde_R;
    gsBSpline<T> beta_tilde_L, beta_tilde_R;

}; // class gsGluingData

template <class T>
void gsGluingData<T>::beta_exact()
{
    // Spline space for beta
    index_t m_k, m_p;

    // Assume for now only uniform refinement TODO: Do it for general knotvectors
    gsBasis<T> & basis_first = m_mb.basis(0).component(0); // interface at u-direction
    gsBasis<T> & basis_second = m_mb.basis(1).component(1); // interface at v-direction
    if (basis_first.numElements() >= basis_second.numElements())
        m_k = basis_second.numElements()-1;
    else
        m_k = basis_first.numElements()-1;

    // Need that for beta_exact
    if (basis_first.maxDegree() >= basis_second.maxDegree())
        m_p = basis_second.maxDegree();
    else
        m_p = basis_first.maxDegree();

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
    uv0.topRows(1) = greville;

    // Set points for Patch 1
    uv1.setZero(2,greville.cols());
    uv1.bottomRows(1) = greville;


    const gsGeometry<> & P0 = m_mp.patch(0);
    const gsGeometry<> & P1 = m_mp.patch(1);
    // ======================================

    // ======== Determine bar{beta} ========
    for(index_t i = 0; i < uv1.cols(); i++)
    {
        P0.jacobian_into(uv0.col(i),ev0);
        P1.jacobian_into(uv1.col(i),ev1);

        D0.col(1) = ev0.col(1); // (DuFL, *)
        D0.col(0) = ev1.col(0); // (*,DuFR)

        uv1(0,i) =  m_gamma * D0.determinant();
    }

    beta_temp = bsp.interpolateData(uv1.topRows(1), greville);
    beta_bar = dynamic_cast<gsBSpline<T> &> (*beta_temp);
} // beta_exact

template<class T>
void gsGluingData<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(0)); // u
    gsBSplineBasis<> temp_basis_second = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(1).component(1)); // v

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

    gsGlobalGDAssembler<T> globalGdAssembler(bsp_gD,m_mp,m_gamma);
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

    if (m_optionList.getSwitch("plot"))
    {
        gsWriteParaview(alpha_tilde_L,"alpha_tilde_L",5000);
        gsWriteParaview(alpha_tilde_R,"alpha_tilde_R",5000);

        gsWriteParaview(beta_tilde_L,"beta_tilde_L",5000);
        gsWriteParaview(beta_tilde_R,"beta_tilde_R",5000);
    }
} // setGlobalGluingData

} // namespace gismo

