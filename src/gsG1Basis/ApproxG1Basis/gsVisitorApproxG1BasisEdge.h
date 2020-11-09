/** @file gsVisitorG1Basis.h

    @brief Visitor for the G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

# include "gsG1Basis/ApproxG1Basis/gsApproxGluingData.h"
# include <gsG1Basis/gsG1OptionList.h>

namespace gismo
{
template <class T>
class gsVisitorApproxG1BasisEdge
{
public:

    gsVisitorApproxG1BasisEdge()
    {
    }

    void initialize(const gsBasis<T>       & basis, //
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE ;
    }

    // Evaluate on element.
    inline void evaluate(const index_t bfID,
                         std::string typeBf,
                         gsBasis<T>       & basis, //
                         gsBasis<T>       & basis_geo,
                         gsBasis<T>       & basis_plus,
                         gsBasis<T>       & basis_minus,
                         const gsGeometry<T>    & geo, // patch
                         gsMatrix<T>            & quNodes,
                         index_t & uv,
                         gsApproxGluingData<T>  & gluingData,
                         bool & isBoundary,
                         gsG1OptionList g1OptionList,
                         gsBSpline<T> & result_singleEdge)
    {
        h1projection = g1OptionList.getSwitch("h1projection");
        h2projection = g1OptionList.getSwitch("h2projection");

        md.points = quNodes;

/*
#pragma omp critical
        {
            gsBSplineBasis<> basis_bspline = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(0));
            gsInfo << "Basis: " << basis_bspline.knots().asMatrix() << "\n";

            gsBSplineBasis<> basis_bspline22 = dynamic_cast<gsBSplineBasis<real_t> &>(basis.component(1));
            gsInfo << "Basis22: " << basis_bspline22.knots().asMatrix() << "\n";

            gsBSplineBasis<> basis_plus2 = dynamic_cast<gsBSplineBasis<real_t> &>(basis_plus);
            gsInfo << "basis_plus2: " << basis_plus2.knots().asMatrix() << "\n";

            gsBSplineBasis<> basis_minus2 = dynamic_cast<gsBSplineBasis<real_t> &>(basis_minus);
            gsInfo << "basis_minus: " << basis_minus2.knots().asMatrix() << "\n";

            gsBSplineBasis<> basis_geo2 = dynamic_cast<gsBSplineBasis<real_t> &>(basis_geo);
            gsInfo << "basis_geo: " << basis_geo2.knots().asMatrix() << "\n";
        }
*/
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        if (h1projection && !h2projection)
            basis.evalAllDers_into( md.points, 1, basisData);
        else if (!h1projection && h2projection)
            basis.evalAllDers_into( md.points, 2, basisData);
        else
            basis.evalAllDers_into( md.points, 0, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        numActive = actives.rows();

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo);

        real_t p = basis_geo.maxDegree();
        real_t tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> alpha, beta,
            N_0, N_1,
            N_j_minus, N_i_plus,
            der_N_i_plus;

        gsMatrix<T> der_alpha, der_beta,
            der_N_0, der_N_1,
            der_N_j_minus,
            der2_N_i_plus;

        gsMatrix<T> der2_alpha, der2_beta,
            der2_N_0, der2_N_1,
            der2_N_j_minus,
            der3_N_i_plus;

        if (uv == 1) // edge is in v-direction == Patch 0
        {
            if (g1OptionList.getInt("gluingData") == gluingData::global)
            {
                gluingData.get_alpha_S_tilde(1-uv).eval_into(md.points.bottomRows(1),alpha); // v
                gluingData.get_beta_S_tilde(1-uv).eval_into(md.points.bottomRows(1),beta);
            }
            else if (g1OptionList.getInt("gluingData") == gluingData::exact)
            {
                gluingData.eval_alpha_into(md.points.bottomRows(1),alpha); // v
                gluingData.eval_beta_into(md.points.bottomRows(1),beta);
            }
            basis_geo.evalSingle_into(0,md.points.topRows(1),N_0); // u
            basis_geo.evalSingle_into(1,md.points.topRows(1),N_1); // u

            if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
            {
                rhsGrads.setZero(2,md.points.cols());

                basis_geo.derivSingle_into(0,md.points.topRows(1),der_N_0); // u
                basis_geo.derivSingle_into(1,md.points.topRows(1),der_N_1); // u

                if (g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gluingData.get_alpha_S_tilde(1-uv).deriv_into(md.points.bottomRows(1),der_alpha); // v
                    gluingData.get_beta_S_tilde(1-uv).deriv_into(md.points.bottomRows(1),der_beta);
                }
            }
            if (g1OptionList.getSwitch("h2projection"))
            {
                rhs2Der.setZero(3,md.points.cols());

                basis_geo.deriv2Single_into(0,md.points.topRows(1),der2_N_0); // u
                basis_geo.deriv2Single_into(1,md.points.topRows(1),der2_N_1); // u

                if (g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gluingData.get_alpha_S_tilde(1-uv).deriv2_into(md.points.bottomRows(1),der2_alpha); // v
                    gluingData.get_beta_S_tilde(1-uv).deriv2_into(md.points.bottomRows(1),der2_beta);
                }
            }


            // Initialize local matrix/rhs
            if (typeBf == "plus")
            {
                basis_plus.evalSingle_into(bfID,md.points.bottomRows(1),N_i_plus); // v
                basis_plus.derivSingle_into(bfID,md.points.bottomRows(1),der_N_i_plus);

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                    basis_plus.deriv2Single_into(bfID,md.points.bottomRows(1),der2_N_i_plus);
                if (g1OptionList.getSwitch("h2projection"))
                    basis_plus.evalDerSingle_into(bfID,md.points.bottomRows(1),3,der3_N_i_plus);

                if (g1OptionList.getInt("gluingData") == gluingData::local)
                {
                    gsMatrix<> ab = gluingData.get_beta_S_tilde(bfID).support();
                    if ((md.points(1,0) >= ab(0)) && (md.points(1,0) <= ab(1)))
                        gluingData.get_beta_S_tilde(bfID).eval_into(md.points.bottomRows(1),beta);
                    else
                        beta.setZero(1,md.points.cols());
                }

                beta = isBoundary ? beta.setZero() : beta; // For the boundary, only on Patch 0

                //if (bfID == 1 || bfID == basis_plus.size()-2)
                //    beta.setZero();

                gsMatrix<T> temp;
/*                if (g1OptionList.getSwitch("twoPatch") && bfID == 1 && g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gsMatrix<> lambda, null(1,1);
                    null << 0.0;
                    lambda = gluingData.get_beta_tilde().eval(null) * 1/(gluingData.get_alpha_tilde().eval(null)(0, 0));
                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else if (g1OptionList.getSwitch("twoPatch") && bfID == basis_plus.size()-2 && g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gsMatrix<> lambda, one(1,1);
                    one << 1.0;
                    lambda = gluingData.get_beta_tilde().eval(one) * 1/(gluingData.get_alpha_tilde().eval(one)(0, 0));
                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else if (g1OptionList.getSwitch("twoPatch") && bfID == 1 && g1OptionList.getInt("gluingData") == gluingData::exact)
                {
                    gsMatrix<> lambda, null(1,1), alpha_0, beta_0;
                    null << 0.0;
                    gluingData.eval_alpha_into(null,alpha_0); // v
                    gluingData.eval_beta_into(null,beta_0);
                    lambda = beta_0 * 1/(alpha_0(0, 0));
                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else if (g1OptionList.getSwitch("twoPatch") && bfID == basis_plus.size()-2 && g1OptionList.getInt("gluingData") == gluingData::exact)
                {
                    gsMatrix<> lambda, one(1,1), alpha_1, beta_1;
                    one << 1.0;
                    gluingData.eval_alpha_into(one,alpha_1); // v
                    gluingData.eval_beta_into(one,beta_1);
                    lambda = beta_1 * 1/(alpha_1(0, 0));
                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else
*/                    temp = beta.cwiseProduct(der_N_i_plus);


                //gsMatrix<> temp = beta.cwiseProduct(der_N_i_plus);

                rhsVals = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                {
                    gsMatrix<> der_temp;

                    der_temp = der_beta.cwiseProduct(der_N_i_plus) + beta.cwiseProduct(der2_N_i_plus);

                    rhsGrads.row(0) = N_i_plus.cwiseProduct(der_N_0 + der_N_1) - temp.cwiseProduct(der_N_1) * tau_1 / p;
                    rhsGrads.row(1) = der_N_i_plus.cwiseProduct(N_0 + N_1) - der_temp.cwiseProduct(N_1) * tau_1 / p;
                }

                if (g1OptionList.getSwitch("h2projection"))
                {
                    gsMatrix<> der2_temp, der_temp;

                    der_temp = der_beta.cwiseProduct(der_N_i_plus) + beta.cwiseProduct(der2_N_i_plus);
                    der2_temp = 2.0*der_beta.cwiseProduct(der2_N_i_plus) + beta.cwiseProduct(der3_N_i_plus) + der2_beta.cwiseProduct(der_N_i_plus);

                    rhs2Der.row(0) = N_i_plus.cwiseProduct(der2_N_0 + der2_N_1) - temp.cwiseProduct(der2_N_1) * tau_1 / p;
                    rhs2Der.row(1) = der2_N_i_plus.cwiseProduct(N_0 + N_1) - der2_temp.cwiseProduct(N_1) * tau_1 / p;
                    rhs2Der.row(2) = der_N_i_plus.cwiseProduct(der_N_0 + der_N_1) - der_temp.cwiseProduct(der_N_1) * tau_1 / p;
                }

                localMat.setZero(numActive, numActive);
                localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides

            } // n_plus
            else if (typeBf == "minus")
            {
                basis_minus.evalSingle_into(bfID,md.points.bottomRows(1),N_j_minus); // v

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection") )
                    basis_minus.derivSingle_into(bfID,md.points.bottomRows(1),der_N_j_minus); // v

                if (g1OptionList.getSwitch("h2projection") )
                    basis_minus.deriv2Single_into(bfID,md.points.bottomRows(1),der2_N_j_minus); // v

                if (g1OptionList.getInt("gluingData") == gluingData::local)
                {
                    gsMatrix<> ab = gluingData.get_alpha_S_tilde(bfID).support();
                    if ((md.points(1,0) >= ab(0)) && (md.points(1,0) <= ab(1)))
                        gluingData.get_alpha_S_tilde(bfID).eval_into(md.points.bottomRows(1),alpha);
                    else
                        alpha.setZero(1,md.points.cols());
                }

                alpha = isBoundary ? alpha.setOnes() : alpha; // For the boundary, only on Patch 0

                rhsVals = alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;

                //rhsVals.setZero();

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                {
                    rhsGrads.row(0) = alpha.cwiseProduct(N_j_minus.cwiseProduct(der_N_1)) * tau_1 / p;
                    rhsGrads.row(1) = (der_alpha.cwiseProduct(N_j_minus)+alpha.cwiseProduct(der_N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
                }
                if (g1OptionList.getSwitch("h2projection"))
                {
                    rhs2Der.row(0) = alpha.cwiseProduct(N_j_minus.cwiseProduct(der2_N_1)) * tau_1 / p;
                    rhs2Der.row(1) = (2.0*der_alpha.cwiseProduct(der_N_j_minus)+alpha.cwiseProduct(der2_N_j_minus)+der2_alpha.cwiseProduct(N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
                    rhs2Der.row(2) = (der_alpha.cwiseProduct(N_j_minus)+alpha.cwiseProduct(der_N_j_minus)).cwiseProduct(der_N_1) * tau_1 / p;
                }

                localMat.setZero(numActive, numActive);
                localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
            } // n_minus

        } // Patch 0
        else if (uv == 0) // edge is in u-direction
        {
            if (g1OptionList.getInt("gluingData") == gluingData::global)
            {
                gluingData.get_alpha_S_tilde(1-uv).eval_into(md.points.topRows(1),alpha); // u
                gluingData.get_beta_S_tilde(1-uv).eval_into(md.points.topRows(1),beta);
            }
            else if (g1OptionList.getInt("gluingData") == gluingData::exact)
            {
                gluingData.eval_alpha_into(md.points.topRows(1),alpha); // v
                gluingData.eval_beta_into(md.points.topRows(1),beta);
            }
            basis_geo.evalSingle_into(0,md.points.bottomRows(1),N_0); // v
            basis_geo.evalSingle_into(1,md.points.bottomRows(1),N_1); // v

            if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
            {
                rhsGrads.setZero(2,md.points.cols());

                basis_geo.derivSingle_into(0,md.points.bottomRows(1),der_N_0); // u
                basis_geo.derivSingle_into(1,md.points.bottomRows(1),der_N_1); // u

                if (g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gluingData.get_alpha_S_tilde(1-uv).deriv_into(md.points.topRows(1),der_alpha); // v
                    gluingData.get_beta_S_tilde(1-uv).deriv_into(md.points.topRows(1),der_beta);
                }
            }

            if (g1OptionList.getSwitch("h2projection"))
            {
                rhs2Der.setZero(3,md.points.cols());

                basis_geo.deriv2Single_into(0,md.points.bottomRows(1),der2_N_0); // u
                basis_geo.deriv2Single_into(1,md.points.bottomRows(1),der2_N_1); // u

                if (g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gluingData.get_alpha_S_tilde(1-uv).deriv2_into(md.points.topRows(1),der2_alpha); // v
                    gluingData.get_beta_S_tilde(1-uv).deriv2_into(md.points.topRows(1),der2_beta);
                }
            }

            // Initialize local matrix/rhs
            if (typeBf == "plus")
            {
                basis_plus.evalSingle_into(bfID,md.points.topRows(1),N_i_plus); // u
                basis_plus.derivSingle_into(bfID,md.points.topRows(1),der_N_i_plus);

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                    basis_plus.deriv2Single_into(bfID,md.points.topRows(1),der2_N_i_plus);

                if (g1OptionList.getSwitch("h2projection"))
                    basis_plus.evalDerSingle_into(bfID,md.points.topRows(1),3,der3_N_i_plus);

                if (g1OptionList.getInt("gluingData") == gluingData::local)
                {
                    gsMatrix<> ab = gluingData.get_beta_S_tilde(bfID).support();
                    if ((md.points(0,0) >= ab(0)) && (md.points(0,0) <= ab(1)))
                        gluingData.get_beta_S_tilde(bfID).eval_into(md.points.topRows(1),beta);
                    else
                        beta.setZero(1,md.points.cols());
                }



                beta = isBoundary ? beta.setZero() : beta; // For the boundary, only on Patch 0

                //if (bfID == 1 || bfID == basis_plus.size()-2)
                //    beta.setZero();

                gsMatrix<T> temp;
/*                if (g1OptionList.getSwitch("twoPatch") && bfID == 1 && g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gsMatrix<> lambda, null(1,1);
                    null << 0.0;
                    lambda = gluingData.get_beta_tilde().eval(null) * 1/(gluingData.get_alpha_tilde().eval(null)(0, 0));

                    lambda(0,0) = -g1OptionList.getReal("lambda");

                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else if (g1OptionList.getSwitch("twoPatch") && bfID == basis_plus.size()-2 && g1OptionList.getInt("gluingData") == gluingData::global)
                {
                    gsMatrix<> lambda, one(1,1);
                    one << 1.0;
                    lambda = gluingData.get_beta_tilde().eval(one) * 1/(gluingData.get_alpha_tilde().eval(one)(0, 0));

                    lambda(0,0) = -g1OptionList.getReal("lambda2");
                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else if (g1OptionList.getSwitch("twoPatch") && bfID == 1 && g1OptionList.getInt("gluingData") == gluingData::exact)
                {
                    gsMatrix<> lambda, null(1,1), alpha_0, beta_0;
                    null << 0.0;
                    gluingData.eval_alpha_into(null,alpha_0); // v
                    gluingData.eval_beta_into(null,beta_0);
                    lambda = beta_0 * 1/(alpha_0(0, 0));
                    //lambda(0,0) = -g1OptionList.getReal("lambda2");

                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else if (g1OptionList.getSwitch("twoPatch") && bfID == basis_plus.size()-2 && g1OptionList.getInt("gluingData") == gluingData::exact)
                {
                    gsMatrix<> lambda, one(1,1), alpha_1, beta_1;
                    one << 1.0;
                    gluingData.eval_alpha_into(one,alpha_1); // v
                    gluingData.eval_beta_into(one,beta_1);
                    lambda = beta_1 * 1/(alpha_1(0, 0));
                    //lambda(0,0) = -g1OptionList.getReal("lambda2");
                    temp = (beta - lambda * alpha).cwiseProduct(der_N_i_plus);
                }
                else
*/                    temp = beta.cwiseProduct(der_N_i_plus);


                //gsMatrix<T> temp = beta.cwiseProduct(der_N_i_plus);

                //gsInfo << "uv = 0 : " << temp - result_singleEdge.eval(md.points.topRows(1)) << "\n";
                //temp = result_singleEdge.eval(md.points.topRows(1));

                rhsVals = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                {
                    gsMatrix<> der_temp;

                    der_temp = der_beta.cwiseProduct(der_N_i_plus) + beta.cwiseProduct(der2_N_i_plus);

                    rhsGrads.row(1) = N_i_plus.cwiseProduct(der_N_0 + der_N_1) - temp.cwiseProduct(der_N_1) * tau_1 / p;
                    rhsGrads.row(0) = der_N_i_plus.cwiseProduct(N_0 + N_1) - der_temp.cwiseProduct(N_1) * tau_1 / p;
                }
                if (g1OptionList.getSwitch("h2projection"))
                {
                    gsMatrix<> der2_temp, der_temp;

                    der_temp = der_beta.cwiseProduct(der_N_i_plus) + beta.cwiseProduct(der2_N_i_plus);
                    der2_temp = 2.0*der_beta.cwiseProduct(der2_N_i_plus) + beta.cwiseProduct(der3_N_i_plus) + der2_beta.cwiseProduct(der_N_i_plus);

                    rhs2Der.row(1) = N_i_plus.cwiseProduct(der2_N_0 + der2_N_1) - temp.cwiseProduct(der2_N_1) * tau_1 / p;
                    rhs2Der.row(0) = der2_N_i_plus.cwiseProduct(N_0 + N_1) - der2_temp.cwiseProduct(N_1) * tau_1 / p;
                    rhs2Der.row(2) = der_N_i_plus.cwiseProduct(der_N_0 + der_N_1) - der_temp.cwiseProduct(der_N_1) * tau_1 / p;
                }

                localMat.setZero(numActive, numActive);
                localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides

            } // n_tilde
            else if (typeBf == "minus")
            {
                basis_minus.evalSingle_into(bfID,md.points.topRows(1),N_j_minus); // u

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                    basis_minus.derivSingle_into(bfID,md.points.topRows(1),der_N_j_minus); // u

                if (g1OptionList.getSwitch("h2projection"))
                    basis_minus.deriv2Single_into(bfID,md.points.topRows(1),der2_N_j_minus); // u

                if (g1OptionList.getInt("gluingData") == gluingData::local)
                {
                    gsMatrix<> ab = gluingData.get_alpha_S_tilde(bfID).support();
                    if ((md.points(0,0) >= ab(0)) && (md.points(0,0) <= ab(1)))
                        gluingData.get_alpha_S_tilde(bfID).eval_into(md.points.topRows(1),alpha);
                    else
                        alpha.setZero(1,md.points.cols());
                }

                alpha = isBoundary ? alpha.setOnes() : alpha; // For the boundary, only on Patch 0

                rhsVals = - alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;

                //gsInfo << "RHSVALS: " << alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p - alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) << "\n";

                if (g1OptionList.getSwitch("h1projection") || g1OptionList.getSwitch("h2projection"))
                {
                    rhsGrads.row(1) = - alpha.cwiseProduct(N_j_minus.cwiseProduct(der_N_1)) * tau_1 / p;
                    rhsGrads.row(0) = - (der_alpha.cwiseProduct(N_j_minus)+alpha.cwiseProduct(der_N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
                }
                if (g1OptionList.getSwitch("h2projection"))
                {
                    rhs2Der.row(1) = - alpha.cwiseProduct(N_j_minus.cwiseProduct(der2_N_1)) * tau_1 / p;
                    rhs2Der.row(0) = - (2.0*der_alpha.cwiseProduct(der_N_j_minus)+alpha.cwiseProduct(der2_N_j_minus)+der2_alpha.cwiseProduct(N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
                    rhs2Der.row(2) = - (der_alpha.cwiseProduct(N_j_minus)+alpha.cwiseProduct(der_N_j_minus)).cwiseProduct(der_N_1) * tau_1 / p;
                }

                localMat.setZero(numActive, numActive);
                localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
            } // n_bar

        } // Patch 1
    } // evaluate1

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];

        // ( u, v)
        localMat.noalias() =
            bVals * quWeights.asDiagonal() * bVals.transpose();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            if (h1projection || h2projection)
            {
                gsMatrix<T> & bGrads = basisData[1];

                const index_t numGrads = bGrads.rows() / md.dim.first;
                const gsAsConstMatrix<T> grads_k(bGrads.col(k).data(), md.dim.first, numGrads);

                localMat.noalias() += weight * (grads_k.transpose() * grads_k);
                localRhs.noalias() += weight * (grads_k.transpose() * rhsGrads.col(k)) ;
            }

            if (h2projection)
            {
                gsMatrix<T> & bGrads2 = basisData[2];

                const index_t numGrads = bGrads2.rows() / 3; // Only planar
                const gsAsConstMatrix<T> grads_k(bGrads2.col(k).data(), 3, numGrads);

                localMat.noalias() += weight * (grads_k.transpose() * grads_k);
                localRhs.noalias() += weight * (grads_k.transpose() * rhs2Der.col(k)) ;
            }


            // Multiply weight by the geometry measure
            localRhs.noalias() += weight * (bVals.col(k) * rhsVals.col(k).transpose());
            //
        }

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>      & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);
        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs[0], 0, 0);
    }

protected:
    gsMatrix<unsigned> actives;
    std::vector<gsMatrix<T> > basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals, rhsGrads, rhs2Der;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;

    bool h1projection, h2projection;

}; // class gsVisitorG1BasisEdge
} // namespace gismo