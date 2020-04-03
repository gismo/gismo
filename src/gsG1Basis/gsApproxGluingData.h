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
    }

    // Computed the gluing data globally
    void setGlobalGluingData();

    // Computed the gluing data locally
    void setLocalGluingData(gsBSplineBasis<> & basis_plus, gsBSplineBasis<> & basis_minus, std::string edgeVertex);

    const gsBSpline<T> get_local_alpha_tilde(index_t i) const {return alpha_minus_tilde[i]; }
    const gsBSpline<T> get_local_beta_tilde(index_t i) const {return beta_plus_tilde[i]; }

    void plotGluingData(index_t nummGd = 0);

protected:

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    std::vector<gsBSpline<T>> alpha_minus_tilde, beta_plus_tilde;

}; // class gsGluingData

template<class T>
void gsApproxGluingData<T>::plotGluingData(index_t numGd)
{
    if (this->m_optionList.getInt("gluingData") == gluingData::global)
    {
        std::string fileName;
        std::string basename = "ApproxGluingDataGlobal" + util::to_string(numGd);
        gsParaviewCollection collection(basename);

        // Alpha
        fileName = basename + "_alpha_" + util::to_string(0);
        gsWriteParaview(this->alpha_tilde,fileName,1000);
        collection.addPart(fileName,".vtp");
        // Beta
        fileName = basename + "_beta_" + util::to_string(0);
        gsWriteParaview(this->beta_tilde,fileName,1000);
        collection.addPart(fileName,".vtp");

        collection.save();
    }
    else if (this->m_optionList.getInt("gluingData") == gluingData::local)
    {
        std::string fileName;
        std::string basename = "ApproxGluingDataLocal2" + util::to_string(numGd);
        gsParaviewCollection collection(basename);

        for (size_t i = 0; i < alpha_minus_tilde.size(); i++)
        {
            fileName = basename + "_alpha_" + util::to_string(i);
            gsWriteParaview(alpha_minus_tilde[i],fileName,1000);
            collection.addTimestep(fileName,i,".vtp");
        }
        for (size_t i = 0; i < beta_plus_tilde.size(); i++)
        {
            fileName = basename + "_beta_" + util::to_string(i);
            gsWriteParaview(beta_plus_tilde[i],fileName,5000);
            collection.addTimestep(fileName,i,".vtp");
        }
        collection.save();
    }
}


template<class T>
void gsApproxGluingData<T>::setGlobalGluingData()
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(this->m_mb.basis(0).component(this->m_uv));

    index_t degree = temp_basis_first.maxDegree();
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
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
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);

    if (edgeVertex == "edge")
    {
        // Compute alpha_minus
        for (index_t bfID = 0; bfID < n_minus; bfID++)
        {
            gsMatrix<T> ab = basis_minus.support(bfID);

            gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);

            index_t degree = temp_basis_first.maxDegree();
            for (size_t ii = degree + 1; ii < temp_basis_first.knots().size() - (degree + 1); ii += temp_basis_first.knots().multiplicityIndex(ii))
                if ((temp_basis_first.knot(ii) > ab.at(0)) && (temp_basis_first.knot(ii) < ab.at(1)))
                    kv.insert(temp_basis_first.knot(ii), 1);

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
            alpha_minus_tilde.at(bfID) = a_t;
        }
        for (index_t bfID = 0; bfID < n_plus; bfID++)
        {
            gsMatrix<T> ab = basis_plus.support(bfID);

            gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);

            index_t degree = temp_basis_first.maxDegree();
            for (size_t ii = degree + 1; ii < temp_basis_first.knots().size() - (degree + 1); ii += temp_basis_first.knots().multiplicityIndex(ii))
                if ((temp_basis_first.knot(ii) > ab.at(0)) && (temp_basis_first.knot(ii) < ab.at(1)))
                    kv.insert(temp_basis_first.knot(ii), 1);

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
            beta_plus_tilde.at(bfID) = b_t;
        }
    } // edge
    else if (edgeVertex == "vertex")
    {
        // ALPHA
        gsMatrix<T> ab = basis_minus.support(1); // FIXED TO SUPP(b_0^-)

        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);

        index_t degree = temp_basis_first.maxDegree();
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), 1);

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
        ab = basis_plus.support(2); // FIXED TO SUPP(b_0^+)

        gsKnotVector<T> kv2(ab.at(0), ab.at(1), 0, p_tilde + 1);

        degree = temp_basis_first.maxDegree();
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv2.insert(temp_basis_first.knot(i), 1);

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

