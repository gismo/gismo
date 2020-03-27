/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

# include <gsG1Basis/gsGlobalGDAssembler.h>
# include <gsG1Basis/gsLocalGDAssembler.h>

# include <gsG1Basis/gsGluingData.h>

# include <gsG1Basis/gsG1OptionList.h>

namespace gismo
{

template<class T>
class gsApproxGluingData : public gsGluingData<T>
{
public:
    gsApproxGluingData()
    { }

    gsApproxGluingData(gsMultiPatch<T> const & mp,
                       gsMultiBasis<T> const & mb,
                       index_t uv,
                       bool isBoundary,
                       gsG1OptionList const & optionList)
        : gsGluingData<T>(mp, mb, uv, isBoundary, optionList)
    {
        p_tilde = this->m_optionList.getInt("p_tilde");
        r_tilde = this->m_optionList.getInt("r_tilde");

        m_r = this->m_optionList.getInt("regularity");

    }

    // Computed the gluing data globally
    void setGlobalGluingData();

    // Computed the gluing data locally
    void setLocalGluingData(gsBSplineBasis<> & basis_plus, gsBSplineBasis<> & basis_minus, std::string edgeVertex);

    const gsBSpline<T> get_local_alpha_tilde(index_t i) const {return alpha_minus_tilde[i]; }
    const gsBSpline<T> get_local_beta_tilde(index_t i) const {return beta_plus_tilde[i]; }

protected:

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    // Regularity of the geometry
    index_t m_r;

    std::vector<gsBSpline<T>> alpha_minus_tilde, beta_plus_tilde;

}; // class gsGluingData


template<class T>
void gsApproxGluingData<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(this->m_mb.basis(0).component(this->m_uv)); // u
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


    gsGlobalGDAssembler<T> globalGdAssembler(bsp_gD, this->m_uv, this->m_mp, this->m_gamma, this->m_isBoundary);

    globalGdAssembler.assemble();

    gsSparseSolver<real_t>::CGDiagonal solver;
    gsVector<> sol_a, sol_b;

    // alpha^S
    solver.compute(globalGdAssembler.matrix_alpha());
    sol_a = solver.solve(globalGdAssembler.rhs_alpha());

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol_a);

    gsBSpline<T> alpha_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    this->alpha_tilde = alpha_t;

    // beta^S
    solver.compute(globalGdAssembler.matrix_beta());
    sol_b = solver.solve(globalGdAssembler.rhs_beta());

    tilde_temp = bsp_gD.makeGeometry(sol_b);

    gsBSpline<T> beta_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    this->beta_tilde = beta_t;


} // setGlobalGluingData

template<class T>
void gsApproxGluingData<T>::setLocalGluingData(gsBSplineBasis<> & basis_plus, gsBSplineBasis<> & basis_minus, std::string edgeVertex)
{
    index_t n_plus = basis_plus.size();
    index_t n_minus = basis_minus.size();

    // Setting the space for each alpha_tilde, beta_tilde
    if (edgeVertex == "edge")
    {
        alpha_minus_tilde.resize(n_minus);
        beta_plus_tilde.resize(n_plus);
    }
    else if (edgeVertex == "vertex")
    {
        alpha_minus_tilde.resize(1);
        beta_plus_tilde.resize(1);
    }

    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(this->m_mb.basis(0).component(this->m_uv)); // u

    index_t degree = temp_basis_first.maxDegree();
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i = i+(degree-m_r))
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);

    if (edgeVertex == "edge")
    {
        // Compute alpha_minus
        for (index_t i = 0; i < n_minus; i++)
        {
            gsMatrix<T> ab = basis_minus.support(i);

            gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, p_tilde + 1);

            index_t degree = temp_basis_first.maxDegree();
            for (size_t ii = degree + 1; ii < temp_basis_first.knots().size() - (degree + 1); ii = ii + (degree - m_r))
                if ((temp_basis_first.knot(ii) > ab.at(0)) && (temp_basis_first.knot(ii) < ab.at(1)))
                    kv.insert(temp_basis_first.knot(ii), p_tilde - r_tilde);
            /*
            real_t span = bsp_gD.getMaxCellLength();
            real_t temp_knot = ab.at(0) + span;
            while (temp_knot < ab.at(1))
            {
                kv.insert(temp_knot,p_tilde-r_tilde);
                temp_knot += span;
            }
             */
            gsBSplineBasis<T> bsp_geo(kv);

            // The first basis (bsp_geo) is for the gd, the second for the integral
            gsLocalGDAssembler<T>
                localGdAssembler(bsp_gD, bsp_geo, this->m_uv, this->m_mp, this->m_gamma, this->m_isBoundary, "alpha");
            localGdAssembler.assemble();

            gsSparseSolver<real_t>::CGDiagonal solver;
            gsVector<> sol;

            // alpha^S
            solver.compute(localGdAssembler.matrix());
            sol = solver.solve(localGdAssembler.rhs());

            gsGeometry<>::uPtr tilde_temp;
            tilde_temp = bsp_gD.makeGeometry(sol);
            gsBSpline<T> a_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
            alpha_minus_tilde.at(i) = a_t;
        }
        for (index_t i = 0; i < n_plus; i++)
        {
            gsMatrix<T> ab = basis_plus.support(i);

            gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, p_tilde + 1);

            index_t degree = temp_basis_first.maxDegree();
            for (size_t ii = degree + 1; ii < temp_basis_first.knots().size() - (degree + 1); ii = ii + (degree - m_r))
                if ((temp_basis_first.knot(ii) > ab.at(0)) && (temp_basis_first.knot(ii) < ab.at(1)))
                    kv.insert(temp_basis_first.knot(ii), p_tilde - r_tilde);
            /*
            real_t span = bsp_gD.getMaxCellLength();
            real_t temp_knot = ab.at(0) + span;
            while (temp_knot < ab.at(1))
            {
                kv.insert(temp_knot,p_tilde-r_tilde);
                temp_knot += span;
            }
             */
            gsBSplineBasis<T> bsp_geo(kv);

            // The first basis (bsp_geo) is for the gd, the second for the integral
            gsLocalGDAssembler<T>
                localGdAssembler(bsp_gD, bsp_geo, this->m_uv, this->m_mp, this->m_gamma, this->m_isBoundary, "beta");
            localGdAssembler.assemble();

            gsSparseSolver<real_t>::CGDiagonal solver;
            gsVector<> sol;

            // alpha^S
            solver.compute(localGdAssembler.matrix());
            sol = solver.solve(localGdAssembler.rhs());

            gsGeometry<>::uPtr tilde_temp;
            tilde_temp = bsp_gD.makeGeometry(sol);
            gsBSpline<T> b_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
            beta_plus_tilde.at(i) = b_t;
        }
    } // edge
    else if (edgeVertex == "vertex")
    {
        // ALPHA
        gsMatrix<T> ab = basis_minus.support(0); // FIXED TO SUPP(b_0^-)

        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, p_tilde + 1);

        index_t degree = temp_basis_first.maxDegree();
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i = i + (degree - m_r))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), p_tilde - r_tilde);

        gsBSplineBasis<T> bsp_geo(kv);

        // The first basis (bsp_geo) is for the gd, the second for the integral
        gsLocalGDAssembler<T>
            localGdAssembler(bsp_gD, bsp_geo, this->m_uv, this->m_mp, this->m_gamma, this->m_isBoundary, "alpha");
        localGdAssembler.assemble();

        gsSparseSolver<real_t>::CGDiagonal solver;
        gsVector<> sol;

        // alpha^S
        solver.compute(localGdAssembler.matrix());
        sol = solver.solve(localGdAssembler.rhs());

        gsGeometry<>::uPtr tilde_temp;
        tilde_temp = bsp_gD.makeGeometry(sol);
        gsBSpline<T> a_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        alpha_minus_tilde.at(0) = a_t;

        // BETA
        ab = basis_plus.support(0); // FIXED TO SUPP(b_0^+)

        gsKnotVector<T> kv2(ab.at(0), ab.at(1), 0, p_tilde + 1);

        degree = temp_basis_first.maxDegree();
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i = i + (degree - m_r))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv2.insert(temp_basis_first.knot(i), p_tilde - r_tilde);

        gsBSplineBasis<T> bsp_geo2(kv2);

        // The first basis (bsp_geo) is for the gd, the second for the integral
        gsLocalGDAssembler<T>
            localGdAssembler2(bsp_gD, bsp_geo2, this->m_uv, this->m_mp, this->m_gamma, this->m_isBoundary, "beta");
        localGdAssembler2.assemble();

        // alpha^S
        solver.compute(localGdAssembler2.matrix());
        sol = solver.solve(localGdAssembler2.rhs());

        gsGeometry<>::uPtr tilde_temp2;
        tilde_temp2 = bsp_gD.makeGeometry(sol);
        gsBSpline<T> b_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp2);
        beta_plus_tilde.at(0) = b_t;
    }

} // setLocalGluingData


} // namespace gismo

