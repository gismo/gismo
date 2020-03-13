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
                 index_t uv,
                 bool isBoundary,
                 gsOptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_uv(uv), m_isBoundary(isBoundary), m_optionList(optionList)
    {
        m_gamma = 1.0;

        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");

        m_r = m_optionList.getInt("regularity");

        if (m_optionList.getSwitch("local"))
            gsInfo << "Is not yet implemented \n";
        else
            setGlobalGluingData();
    }


    // Computed the gluing data globally
    void setGlobalGluingData();

    void beta_exact();

    const gsBSpline<T> get_alpha_tilde() const {return alpha_tilde; }
    const gsBSpline<T> get_beta_tilde() const {return beta_tilde; }

protected:
    // The geometry for a single interface in the right parametrizations
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    index_t m_uv;
    bool m_isBoundary;
    gsOptionList m_optionList;

    real_t m_gamma;

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    // Regularity of the geometry
    index_t m_r;

protected:
    // Global Gluing data
    gsBSpline<T> alpha_tilde;
    gsBSpline<T> beta_tilde;

}; // class gsGluingData


template<class T>
void gsGluingData<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(m_uv)); // u
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

    gsGlobalGDAssembler<T> globalGdAssembler(bsp_gD,m_uv,m_mp,m_gamma,m_isBoundary);
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

} // setGlobalGluingData

} // namespace gismo

