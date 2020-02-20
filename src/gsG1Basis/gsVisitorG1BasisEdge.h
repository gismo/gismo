/** @file gsVisitorG1Basis.h

    @brief Visitor for the G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include "gsG1Basis/gsGluingData.h"

namespace gismo
{
template <class T>
class gsVisitorG1BasisEdge
{
public:

    gsVisitorG1BasisEdge()
    {
    }

    void initialize(const gsBasis<T>       & basis, //
                    const gsBasis<T>       & basis_plus,
                    const gsBasis<T>       & basis_minus,
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE ;

        n_plus = basis_plus.size();
        n_minus = basis_minus.size();


        localMat_tilde.resize(n_plus);
        localRhs_tilde.resize(n_plus);

        rhsVals_tilde.resize(n_plus);

        localMat_bar.resize(n_minus);
        localRhs_bar.resize(n_minus);

        rhsVals_bar.resize(n_minus);
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T>       & basis, //
                         gsBasis<T>       & basis_geo,
                         gsBasis<T>       & basis_plus,
                         gsBasis<T>       & basis_minus,
                         const gsGeometry<T>    & geo, // patch
                         gsMatrix<T>            & quNodes,
                         gsGluingData<T>  & gluingData,
                         gsOptionList optionList)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        basis.eval_into(md.points, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        numActive = actives.rows();

        gsMatrix<unsigned> actives_plus, actives_minus;
        basis_plus.active_into(md.points.col(0), actives_plus);
        basis_minus.active_into(md.points.col(0), actives_minus);

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo);

        real_t p = basis_geo.maxDegree();
        real_t tau_1 = bsp_temp.knots().at(p + 2);

        gsMatrix<T> alpha_0, alpha_1, beta_0, beta_1,
                    N_0_patch0, N_1_patch0, N_0_patch1, N_1_patch1,
                    N_j_minus_patch0, N_i_plus_patch0, N_j_minus_patch1, N_i_plus_patch1,
                    der_N_i_plus_patch0, der_N_i_plus_patch1;

        gsMatrix<T> temp_0, temp_1;

        if (optionList.getSwitch("direct"))
        {
            gluingData.eval_into_alpha_0(md.points,alpha_0);
            gluingData.eval_into_alpha_1(md.points,alpha_1);

            gluingData.eval_into_beta_0(md.points,beta_0);
            gluingData.eval_into_beta_1(md.points,beta_1);
        }
        else
        {

            gluingData.get_alpha_tilde_0().eval_into(md.points.bottomRows(1),alpha_0); // v
            gluingData.get_beta_tilde_0().eval_into(md.points.bottomRows(1),beta_0);

            gluingData.get_alpha_tilde_1().eval_into(md.points.topRows(1),alpha_1); // u
            gluingData.get_beta_tilde_1().eval_into(md.points.topRows(1),beta_1);
        }

        basis_geo.evalSingle_into(0,md.points.topRows(1),N_0_patch0); // u
        basis_geo.evalSingle_into(1,md.points.topRows(1),N_1_patch0); // u

        basis_geo.evalSingle_into(0,md.points.bottomRows(1),N_0_patch1); // v
        basis_geo.evalSingle_into(1,md.points.bottomRows(1),N_1_patch1); // v

        // Initialize local matrix/rhs
        for (index_t i = 0; i < n_plus; i++)
        {
            basis_plus.evalSingle_into(i,md.points.bottomRows(1),N_i_plus_patch0); // v
            basis_plus.derivSingle_into(i,md.points.bottomRows(1),der_N_i_plus_patch0);

            basis_plus.evalSingle_into(i,md.points.topRows(1),N_i_plus_patch1); // u
            basis_plus.derivSingle_into(i,md.points.topRows(1),der_N_i_plus_patch1);

            if (optionList.getSwitch("local"))
            {

            }

            temp_0 = beta_0.cwiseProduct(N_1_patch0);
            temp_1 = beta_1.cwiseProduct(N_1_patch1);

            if  ( i == 1 || i == n_plus -2 || i == n_plus -3 || i == 2 )
            {
                gsMatrix<T> lambda_L, lambda_R;

                gsMatrix<T> nulleins(1,1);
                index_t ii = 0;
                if (i == 1 || i == 2)
                {
                    ii = i;
                    nulleins << 0.0;
                }
                if (i == n_plus -3 || i == n_plus -2 )
                {
                    if ( i == n_plus -3)
                        ii = n_minus -3;
                    if ( i == n_plus -2)
                        ii = n_minus -2;
                    nulleins << 1.0;
                }
                if (optionList.getSwitch("direct"))
                {

                }
                else if (optionList.getSwitch("local"))
                {

                }
                else
                {
                    lambda_L = gluingData.get_beta_tilde_0().eval(nulleins) * 1
                        / (gluingData.get_alpha_tilde_0().eval(nulleins)(0, 0));
                    lambda_R = gluingData.get_beta_tilde_1().eval(nulleins) * 1
                        / (gluingData.get_alpha_tilde_1().eval(nulleins)(0, 0));
                }

                if ( i == 1 || i == n_plus -2 ) // MODIFY
                {
                    gsMatrix<T> temp_L_tt = (beta_0 - lambda_L * alpha_0).cwiseProduct(N_1_patch0);
                    gsMatrix<T> temp_R_tt = (beta_1 - lambda_R * alpha_1).cwiseProduct(N_1_patch1); // lambda_L == lambda_R

                    if (geo.id() == 0) // left
                        rhsVals_tilde.at(i) =
                            N_i_plus_patch0.cwiseProduct(N_0_patch0 + N_1_patch0) - temp_L_tt.cwiseProduct(der_N_i_plus_patch0) * tau_1 / p;
                    if (geo.id() == 1) // right
                        rhsVals_tilde.at(i) =
                            N_i_plus_patch1.cwiseProduct(N_0_patch1 + N_1_patch1) - temp_R_tt.cwiseProduct(der_N_i_plus_patch1) * tau_1 / p;
                }
                else
                {
                    if (geo.id() == 0) // left
                        rhsVals_tilde.at(i) =
                            N_i_plus_patch0.cwiseProduct(N_0_patch0 + N_1_patch0) - temp_0.cwiseProduct(der_N_i_plus_patch0) * tau_1 / p;
                    if (geo.id() == 1) // right
                        rhsVals_tilde.at(i) =
                            N_i_plus_patch1.cwiseProduct(N_0_patch1 + N_1_patch1) - temp_1.cwiseProduct(der_N_i_plus_patch1) * tau_1 / p;
                }
            }
            else
            {
                // MINUS für lineare interfaces!!!
                if (geo.id() == 0) // left
                    rhsVals_tilde.at(i) = N_i_plus_patch0.cwiseProduct(N_0_patch0 + N_1_patch0) - temp_0.cwiseProduct(der_N_i_plus_patch0) * tau_1 / p;
                if (geo.id() == 1) // right
                    rhsVals_tilde.at(i) = N_i_plus_patch1.cwiseProduct(N_0_patch1 + N_1_patch1) - temp_1.cwiseProduct(der_N_i_plus_patch1) * tau_1 / p;
            }

            localMat_tilde.at(i).setZero(numActive, numActive);
            localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

            localMat_tilde.at(i).setZero(numActive, numActive);
            localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

        } // n_tilde

        // Initialize local matrix/rhs
        for (index_t i = 0; i < n_minus; i++)
        {

            basis_minus.evalSingle_into(i,md.points.bottomRows(1),N_j_minus_patch0); // v
            basis_minus.evalSingle_into(i,md.points.topRows(1),N_j_minus_patch1); // u

            if (optionList.getSwitch("local"))
            {

            }

            if (geo.id() == 0) // left
                rhsVals_bar.at(i) = alpha_0.cwiseProduct(N_j_minus_patch0.cwiseProduct(N_1_patch0)) * tau_1 / p;
            if (geo.id() == 1) // right
                rhsVals_bar.at(i) = - alpha_1.cwiseProduct(N_j_minus_patch1.cwiseProduct(N_1_patch1)) * tau_1 / p;

            localMat_bar.at(i).setZero(numActive, numActive);
            localRhs_bar.at(i).setZero(numActive, rhsVals_bar.at(i).rows());//multiple right-hand sides

            localMat_bar.at(i).setZero(numActive, numActive);
            localRhs_bar.at(i).setZero(numActive, rhsVals_bar.at(i).rows());//multiple right-hand sides
        } // n_bar

    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;
        for (index_t i = 0; i < n_plus; i++)
        {
            // ( u, v)
            localMat_tilde.at(i).noalias() =
                basisData * quWeights.asDiagonal() *
                    md.measures.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k] * md.measure(k);

                localRhs_tilde.at(i).noalias() += weight * (basisVals.col(k) * rhsVals_tilde.at(i).col(k).transpose());
            }
        }
        for (index_t i = 0; i < n_minus; i++)
        {
            // ( u, v)
            localMat_bar.at(i).noalias() =
                basisData * quWeights.asDiagonal() *
                    md.measures.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k] * md.measure(k);
                localRhs_bar.at(i).noalias() += weight * (basisVals.col(k) * rhsVals_bar.at(i).col(k).transpose());
            }
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              std::vector< gsSparseSystem<T> >     & system_0,
                              std::vector< gsSparseSystem<T> >     & system_1)
    {
        gsMatrix<unsigned> actives_temp;
        for (unsigned i = 0; i < system_0.size(); i++) // n_tilde oder n_bar
        {
            // Map patch-local DoFs to global DoFs
            system_0.at(i).mapColIndices(actives, patchIndex, actives_temp);
            // Add contributions to the system matrix and right-hand side
            system_0.at(i).push(localMat_tilde.at(i), localRhs_tilde.at(i), actives_temp, eliminatedDofs[0], 0, 0);
        }
        for (unsigned i = 0; i < system_1.size(); i++) // n_tilde oder n_bar
        {
            // Map patch-local DoFs to global DoFs
            system_1.at(i).mapColIndices(actives, patchIndex, actives_temp);

            // Add contributions to the system matrix and right-hand side
            system_1.at(i).push(localMat_bar.at(i), localRhs_bar.at(i), actives_temp, eliminatedDofs[1], 0, 0);
        }

    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    std::vector< gsMatrix<T> >  rhsVals_tilde;

    std::vector< gsMatrix<T> >  rhsVals_bar;

protected:
    // Local matrices
    std::vector< gsMatrix<T> > localMat_tilde;
    std::vector< gsMatrix<T> > localRhs_tilde;


    // Local matrices
    std::vector< gsMatrix<T> > localMat_bar;
    std::vector< gsMatrix<T> > localRhs_bar;

    gsMapData<T> md;

    index_t n_plus, n_minus;

}; // class gsVisitorG1BasisEdge
} // namespace gismo