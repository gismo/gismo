/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

# include <gsG1Basis/ApproxG1Basis/gsGlobalGDAssembler.h>
# include <gsG1Basis/ApproxG1Basis/gsLocalGDAssembler.h>

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
                       index_t dir,
                       bool isBoundary,
                       gsG1OptionList const & optionList)
        : gsGluingData<T>(mp, mb, dir, isBoundary, optionList)
    {
        p_tilde = this->m_optionList.getInt("p_tilde");
        r_tilde = this->m_optionList.getInt("r_tilde");
    }


    gsApproxGluingData(gsMultiPatch<T> const & mp,
                       gsMultiBasis<T> const & mb,
                       bool isBoundary,
                       gsG1OptionList const & optionList)
        : gsGluingData<T>(mp, mb, isBoundary, optionList)
    {
        p_tilde = this->m_optionList.getInt("p_tilde");
        r_tilde = this->m_optionList.getInt("r_tilde");

        if (this->m_optionList.getInt("gluingData") == gluingData::local)
            gsInfo << "!!!!!! LOCAL GLUING DATA NOT YET IMPLEMENTED!!!!!!" << "\n"; //setLocalGluingData(basis_plus, basis_minus, "edge");
        else if (this->m_optionList.getInt("gluingData") == gluingData::global)
        {
            if (mp.nPatches() == 2)
            {
                setGlobalGluingData(0,1); // Order is important!!!
                setGlobalGluingData(1,0);
            }
            else if (mp.nPatches() == 1)
                setGlobalGluingData(0,1);
            else
                gsInfo << "Gluing data does not work for #patches > 3 \n";
        }
    }

    // Computed the gluing data globally
    void setGlobalGluingData(index_t patchID = 0, index_t uv = 1);

    // Computed the gluing data locally
    void setLocalGluingData(gsBSplineBasis<> & basis_plus, gsBSplineBasis<> & basis_minus, std::string edgeVertex);

    const gsBSpline<T> get_alpha_S_tilde(index_t i) const {return alpha_S_tilde[i]; }
    const gsBSpline<T> get_beta_S_tilde(index_t i) const {return beta_S_tilde[i]; }

    void set_beta_tilde(gsBSpline<T> beta_t) { this->beta_tilde = beta_t; }

    void plotGluingData(index_t nummGd = 0);

    void eval_alpha_into(gsMatrix<T> const points, gsMatrix<> & result)
    {
        result.clear();

        gsMatrix<> uv, ev;
        // alpha^S
        if (this->m_uv==1)
        {
            uv.setZero(2,points.cols());
            uv.bottomRows(1) = points; // v
        }
        else if (this->m_uv==0)
        {
            uv.setZero(2,points.cols());
            uv.topRows(1) = points; // u
        }

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & P0 = this->m_mp.patch(0); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            uv(0, i) = this->m_gamma * ev.determinant();
        }
        if (this->m_isBoundary)
            uv.setOnes();
        result = uv.row(0);
    }

    void eval_beta_into(gsMatrix<T> const points, gsMatrix<> & result)
    {
        result.clear();

        gsMatrix<> uv, ev;
        // beta^S
        if (this->m_uv==1)
        {
            uv.setZero(2,points.cols());
            uv.bottomRows(1) = points; // v
        }
        else if (this->m_uv==0)
        {
            uv.setZero(2,points.cols());
            uv.topRows(1) = points; // u
        }

        const index_t d = this->m_mp.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        const gsGeometry<> & P0 = this->m_mp.patch(0); // iFace.second().patch = 0

        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(this->m_uv);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - this->m_gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);

        }
        if (this->m_isBoundary)
            uv.setZero();
        result = uv.row(0);
    }

    void deriv_alpha_into(gsMatrix<T> const points, gsMatrix<> & result)
    {
        result.clear();

        gsMatrix<> uv, ev, ev2;
        // alpha^S
        if (this->m_uv==1)
        {
            uv.setZero(2,points.cols());
            uv.bottomRows(1) = points; // v
        }
        else if (this->m_uv==0)
        {
            uv.setZero(2,points.cols());
            uv.topRows(1) = points; // u
        }

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & P0 = this->m_mp.patch(0); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            P0.deriv2_into(uv.col(i), ev2);
            if (this->m_uv == 1)
                uv(0, i) = this->m_gamma * (ev2(2,0)*ev(1,1) + ev2(4,0)*ev(0,0) -
                    ev2(1,0)*ev(1,0) - ev2(5,0)*ev(0,1));
            else if (this->m_uv == 0)
                uv(0, i) = this->m_gamma * (ev2(0,0)*ev(1,1) + ev2(5,0)*ev(0,0) -
                    ev2(2,0)*ev(1,0) - ev2(3,0)*ev(0,1));
        }
        if (this->m_isBoundary)
            uv.setZero();
        result = uv.row(0);
    }

    void deriv_beta_into(gsMatrix<T> const points, gsMatrix<> & result)
    {
        result.clear();

        gsMatrix<> uv, ev, ev2;
        // beta^S
        if (this->m_uv==1)
        {
            uv.setZero(2,points.cols());
            uv.bottomRows(1) = points; // v
        }
        else if (this->m_uv==0)
        {
            uv.setZero(2,points.cols());
            uv.topRows(1) = points; // u
        }

        const index_t d = this->m_mp.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        const gsGeometry<> & P0 = this->m_mp.patch(0); // iFace.second().patch = 0

        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            P0.deriv2_into(uv.col(i), ev2);
            D0 = ev.col(this->m_uv);
            real_t D1 = 1/ D0.squaredNorm();
            real_t D2 = D0.squaredNorm();
            if (this->m_uv == 1)
                uv(0,i) = - this->m_gamma * D1 * D1 * (D2*(ev2(2,0)*ev(0,1) + ev2(1,0)*ev(0,0)+
                    ev2(5,0)*ev(1,1) + ev2(4,0)*ev(1,0)) -
                    (ev.col(1).transpose() * ev.col(0))(0,0) * 2.0 * (ev2(1,0)*ev(0,1) + ev2(4,0)*ev(1,1)));
            else if (this->m_uv == 0)
                uv(0,i) = - this->m_gamma * D1 * D1 * (D2*(ev2(0,0)*ev(0,1) + ev2(2,0)*ev(0,0)+
                    ev2(3,0)*ev(1,1) + ev2(5,0)*ev(1,0)) -
                    (ev.col(1).transpose() * ev.col(0))(0,0) * 2.0 * (ev2(0,0)*ev(0,0) + ev2(3,0)*ev(1,0)));

        }
        if (this->m_isBoundary)
            uv.setZero();
        result = uv.row(0);
    }

protected:

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

    std::vector<gsBSpline<T>> alpha_S_tilde, beta_S_tilde;

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
        std::string basename = "ApproxGluingDataLocal" + util::to_string(numGd);
        gsParaviewCollection collection(basename);

        for (size_t i = 0; i < alpha_S_tilde.size(); i++)
        {
            fileName = basename + "_alpha_" + util::to_string(i);
            gsWriteParaview(alpha_S_tilde[i],fileName,1000);
            collection.addTimestep(fileName,i,".vtp");
        }
        for (size_t i = 0; i < beta_S_tilde.size(); i++)
        {
            fileName = basename + "_beta_" + util::to_string(i);
            gsWriteParaview(beta_S_tilde[i],fileName,5000);
            collection.addTimestep(fileName,i,".vtp");
        }
        collection.save();
    }
}


template<class T>
void gsApproxGluingData<T>::setGlobalGluingData(index_t patchID, index_t uv)
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
    gsBSplineBasis<T> bsp_gD(kv);

    gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(this->m_mb.basis(patchID).component(uv));

    index_t degree = temp_basis_first.maxDegree();
    for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
        bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);

    gsGlobalGDAssembler<T> globalGdAssembler(bsp_gD, uv, patchID, this->m_mp, this->m_gamma, this->m_isBoundary, this->m_optionList.getSwitch("twoPatch"));

    globalGdAssembler.assemble();

    gsSparseSolver<real_t>::LU solver;
    gsVector<> sol_a, sol_b;

    // alpha^S
    if (globalGdAssembler.matrix_alpha().rows() != 0)
    {
        solver.compute(globalGdAssembler.matrix_alpha());
        sol_a = solver.solve(globalGdAssembler.rhs_alpha());
    }

    gsGeometry<>::uPtr tilde_temp;
    if (this->m_optionList.getSwitch("twoPatch"))
    {
        gsVector<> sol_a_new(sol_a.rows() + 2);
        sol_a_new.setZero();
        sol_a_new.block(1, 0, sol_a.rows(), 1) = sol_a;
        sol_a_new.at(0) = globalGdAssembler.bdy_alpha()(0, 0);
        sol_a_new.at(sol_a.rows() + 1) = globalGdAssembler.bdy_alpha()(1, 0);


        tilde_temp = bsp_gD.makeGeometry(sol_a_new);
    }
    else
        tilde_temp = bsp_gD.makeGeometry(sol_a);

    gsBSpline<T> alpha_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    alpha_S_tilde.push_back(alpha_t);
/*
    if (patchID == 0)
        gsWriteParaview(alpha_t, "alpha_L", 1000);
    if (patchID == 1)
        gsWriteParaview(alpha_t, "alpha_R", 1000);
*/
    // beta^S
    if (globalGdAssembler.matrix_beta().rows() != 0)
    {
        solver.compute(globalGdAssembler.matrix_beta());
        sol_b = solver.solve(globalGdAssembler.rhs_beta());
    }
    if (this->m_optionList.getSwitch("twoPatch"))
    {
        gsVector<> sol_b_new(sol_b.rows() + 2);
        sol_b_new.setZero();
        sol_b_new.block(1, 0, sol_b.rows(), 1) = sol_b;
        sol_b_new.at(0) = globalGdAssembler.bdy_beta()(0, 0);
        sol_b_new.at(sol_b.rows() + 1) = globalGdAssembler.bdy_beta()(1, 0);
        tilde_temp = bsp_gD.makeGeometry(sol_b_new);
    }
    else
        tilde_temp = bsp_gD.makeGeometry(sol_b);

    gsBSpline<T> beta_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
    beta_S_tilde.push_back(beta_t);

    if (patchID == 0)
        gsWriteParaview(alpha_t, "alpha_L", 1000);
    if (patchID == 1)
        gsWriteParaview(alpha_t, "alpha_R", 1000);

    if (patchID == 0)
        gsWriteParaview(beta_t, "beta_L", 1000);
    if (patchID == 1)
        gsWriteParaview(beta_t, "beta_R", 1000);

} // setGlobalGluingData

template<class T>
void gsApproxGluingData<T>::setLocalGluingData(gsBSplineBasis<> & basis_plus, gsBSplineBasis<> & basis_minus, std::string edgeVertex)
{
    index_t n_plus = basis_plus.size();
    index_t n_minus = basis_minus.size();

    // Setting the space for each alpha_tilde, beta_tilde
    if (edgeVertex == "edge")
    {
        alpha_S_tilde.resize(n_minus);
        beta_S_tilde.resize(n_plus);
    }
    else if (edgeVertex == "vertex")
    {
        alpha_S_tilde.resize(1);
        beta_S_tilde.resize(1);
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

            gsSparseSolver<real_t>::LU solver;
            gsVector<> sol;

            // alpha^S
            solver.compute(localGdAssembler.matrix());
            sol = solver.solve(localGdAssembler.rhs());

            gsMatrix<T> coeffs;
            localGdAssembler.constructSolution(sol,coeffs);

            gsGeometry<>::uPtr tilde_temp;
            tilde_temp = bsp_gD.makeGeometry(coeffs);
            gsBSpline<T> a_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
            alpha_S_tilde.at(bfID) = a_t;
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

            gsSparseSolver<real_t>::LU solver;
            gsVector<> sol;

            // beta^S
            solver.compute(localGdAssembler.matrix());
            sol = solver.solve(localGdAssembler.rhs());

            gsMatrix<T> coeffs;
            localGdAssembler.constructSolution(sol,coeffs);

            gsGeometry<>::uPtr tilde_temp;
            tilde_temp = bsp_gD.makeGeometry(coeffs);
            gsBSpline<T> b_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
            beta_S_tilde.at(bfID) = b_t;
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

        gsSparseSolver<real_t>::LU solver;
        gsVector<> sol;

        // alpha^S
        solver.compute(localGdAssembler.matrix());
        sol = solver.solve(localGdAssembler.rhs());

        gsMatrix<T> coeffs;
        localGdAssembler.constructSolution(sol,coeffs);

        gsGeometry<>::uPtr tilde_temp;
        tilde_temp = bsp_gD.makeGeometry(coeffs);
        gsBSpline<T> a_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp);
        alpha_S_tilde.at(0) = a_t;

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

        gsMatrix<T> coeffs2;
        localGdAssembler2.constructSolution(sol,coeffs2);

        gsGeometry<>::uPtr tilde_temp2;
        tilde_temp2 = bsp_gD.makeGeometry(coeffs2);
        gsBSpline<T> b_t = dynamic_cast<gsBSpline<T> &> (*tilde_temp2);
        beta_S_tilde.at(0) = b_t;
    }

} // setLocalGluingData


} // namespace gismo

