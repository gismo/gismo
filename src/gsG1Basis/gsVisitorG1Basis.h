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
class gsVisitorG1Basis
{
public:

    gsVisitorG1Basis()
    {
    }

    void initialize(const gsBasis<T>       & basis, //
                    const gsBasis<T>       & basis_geo,
                    const gsBasis<T>       & basis_tilde,
                    const gsBasis<T>       & basis_bar,
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE ;

        n_tilde = basis_tilde.size();
        n_bar = basis_bar.size();


        localMat_tilde.resize(n_tilde);
        localRhs_tilde.resize(n_tilde);

        rhsVals_tilde.resize(n_tilde);

        localMat_bar.resize(n_bar);
        localRhs_bar.resize(n_bar);

        rhsVals_bar.resize(n_bar);
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T>       & basis, //
                         gsBasis<T>       & basis_geo,
                         gsBasis<T>       & basis_tilde,
                         gsBasis<T>       & basis_bar,
                         const gsGeometry<T>    & geo, // patch
                         gsMatrix<T>            & quNodes,
                         gsGluingData<T>  & gluingData,
                         const bool & direct,
                         const bool & local)
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

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo);

        real_t p = basis_geo.maxDegree();
        real_t tau_1 = bsp_temp.knots().at(p + 2);

        // Evaluate Gluing data
        gsMatrix<T> a_t_L, a_t_R, N_0_L, N_1_L, N_j_b, N_i_t, b_t_L, b_t_R, der_N_i_t, N_0_R, N_1_R;
        gsMatrix<T> temp_L, temp_R;

        if (direct)
        {
            gluingData.eval_into_alpha_L(md.points,a_t_L);
            gluingData.eval_into_alpha_R(md.points,a_t_R);

            gluingData.eval_into_beta_L(md.points,b_t_L);
            gluingData.eval_into_beta_R(md.points,b_t_R);
        }
        else
        {
            gluingData.get_alpha_tilde_L().eval_into(md.points.bottomRows(1),a_t_L);
            gluingData.get_alpha_tilde_R().eval_into(md.points.bottomRows(1),a_t_R);

            gluingData.get_beta_tilde_L().eval_into(md.points.bottomRows(1),b_t_L);
            gluingData.get_beta_tilde_R().eval_into(md.points.bottomRows(1),b_t_R);
        }

        basis_geo.evalSingle_into(basis_geo.size()-1,md.points.topRows(1),N_0_L); // u
        basis_geo.evalSingle_into(basis_geo.size()-2,md.points.topRows(1),N_1_L); // u

        basis_geo.evalSingle_into(0,md.points.topRows(1),N_0_R); // u
        basis_geo.evalSingle_into(1,md.points.topRows(1),N_1_R); // u

        // Initialize local matrix/rhs
        for (index_t i = 0; i < n_tilde; i++)
        {
            basis_tilde.evalSingle_into(i,md.points.bottomRows(1),N_i_t);
            basis_tilde.derivSingle_into(i,md.points.bottomRows(1),der_N_i_t);

            if (local)
            {
                gluingData.get_beta_i_tilde_L(i).eval_into(md.points.bottomRows(1),b_t_L);
                gluingData.get_beta_i_tilde_R(i).eval_into(md.points.bottomRows(1),b_t_R);
            }

            temp_L = b_t_L.cwiseProduct(N_1_L);
            temp_R = b_t_R.cwiseProduct(N_1_R);

            if  ( i == 1 || i == n_tilde -2 || i == n_tilde -3 || i == 2 )
            {
                gsMatrix<T> lambda_L, lambda_R;

                gsMatrix<T> nulleins(1,1);
                index_t ii = 0;
                if (i == 1 || i == 2)
                {
                    ii = i;
                    nulleins << 0.0;
                }
                if (i == n_tilde -3 || i == n_tilde -2 )
                {
                    if ( i == n_tilde -3)
                        ii = n_bar -3;
                    if ( i == n_tilde -2)
                        ii = n_bar -2;
                    nulleins << 1.0;
                }
                if (direct)
                {
                    gsMatrix<T> a_L, b_L, a_R, b_R;
                    gsMatrix<T> nulleins_L(2,1),  nulleins_R(2,1);
                    nulleins_L.setOnes();
                    nulleins_L.row(1) = nulleins;

                    nulleins_R.setZero();
                    nulleins_R.row(1) = nulleins;

                    gluingData.eval_into_alpha_L(nulleins_L,a_L);
                    gluingData.eval_into_alpha_R(nulleins_R,a_R);

                    gluingData.eval_into_beta_L(nulleins_L,b_L);
                    gluingData.eval_into_beta_R(nulleins_R,b_R);

                    lambda_L = b_L * 1 / (a_L(0, 0));
                    lambda_R = b_R * 1 / (a_R(0, 0));
                }
                else if (local)
                {
                    gluingData.get_alpha_j_tilde_L(ii).eval_into(md.points.bottomRows(1),a_t_L);
                    gluingData.get_alpha_j_tilde_R(ii).eval_into(md.points.bottomRows(1),a_t_R);

                    lambda_L = gluingData.get_beta_i_tilde_L(i).eval(nulleins) * 1
                        / (gluingData.get_alpha_j_tilde_L(ii).eval(nulleins)(0, 0));
                    lambda_R = gluingData.get_beta_i_tilde_R(i).eval(nulleins) * 1
                        / (gluingData.get_alpha_j_tilde_R(ii).eval(nulleins)(0, 0));
                }
                else
                {
                    lambda_L = gluingData.get_beta_tilde_L().eval(nulleins) * 1
                        / (gluingData.get_alpha_tilde_L().eval(nulleins)(0, 0));
                    lambda_R = gluingData.get_beta_tilde_R().eval(nulleins) * 1
                        / (gluingData.get_alpha_tilde_R().eval(nulleins)(0, 0));
                }

                if (i == 1 || i == n_tilde -2 ) // MODIFY
                {
                    gsMatrix<T> temp_L_tt = (b_t_L - lambda_L * a_t_L).cwiseProduct(N_1_L);
                    gsMatrix<T> temp_R_tt = (b_t_R - lambda_L * a_t_R).cwiseProduct(N_1_R); // lambda_L == lambda_R

                    if (geo.id() == 0) // left
                        rhsVals_tilde.at(i) =
                            N_i_t.cwiseProduct(N_0_L + N_1_L) - temp_L_tt.cwiseProduct(der_N_i_t) * tau_1 / p;
                    if (geo.id() == 1) // right
                        rhsVals_tilde.at(i) =
                            N_i_t.cwiseProduct(N_0_R + N_1_R) - temp_R_tt.cwiseProduct(der_N_i_t) * tau_1 / p;
                }
                else
                {
                    if (geo.id() == 0) // left
                        rhsVals_tilde.at(i) =
                            N_i_t.cwiseProduct(N_0_L + N_1_L) - temp_L.cwiseProduct(der_N_i_t) * tau_1 / p;
                    if (geo.id() == 1) // right
                        rhsVals_tilde.at(i) =
                           N_i_t.cwiseProduct(N_0_R + N_1_R) - temp_R.cwiseProduct(der_N_i_t) * tau_1 / p;
                }
            }
            else
            {
                // MINUS für lineare interfaces!!!
                if (geo.id() == 0) // left
                    rhsVals_tilde.at(i) = N_i_t.cwiseProduct(N_0_L + N_1_L) - temp_L.cwiseProduct(der_N_i_t) * tau_1 / p;
                if (geo.id() == 1) // right
                    rhsVals_tilde.at(i) = N_i_t.cwiseProduct(N_0_R + N_1_R) - temp_R.cwiseProduct(der_N_i_t) * tau_1 / p;
            }

            localMat_tilde.at(i).setZero(numActive, numActive);
            localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

            localMat_tilde.at(i).setZero(numActive, numActive);
            localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

        } // n_tilde

        // Initialize local matrix/rhs
        for (index_t i = 0; i < n_bar; i++)
        {
            basis_bar.evalSingle_into(i,md.points.bottomRows(1),N_j_b);

            if (local)
            {
                gluingData.get_alpha_j_tilde_L(i).eval_into(md.points.bottomRows(1),a_t_L);
                gluingData.get_alpha_j_tilde_R(i).eval_into(md.points.bottomRows(1),a_t_R);
            }

            if (geo.id() == 0) // left
                rhsVals_bar.at(i) = a_t_L.cwiseProduct(N_j_b.cwiseProduct(N_1_L)) * tau_1 / p;
            if (geo.id() == 1) // right
                rhsVals_bar.at(i) = a_t_R.cwiseProduct(N_j_b.cwiseProduct(N_1_R)) * tau_1 / p;

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
        for (index_t i = 0; i < n_tilde; i++)
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
        for (index_t i = 0; i < n_bar; i++)
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
                              std::vector< gsSparseSystem<T> >     & system_t,
                              std::vector< gsSparseSystem<T> >     & system_b)
    {
        gsMatrix<unsigned> actives_temp;
        for (unsigned i = 0; i < system_t.size(); i++) // n_tilde oder n_bar
        {
            // Map patch-local DoFs to global DoFs
            system_t.at(i).mapColIndices(actives, patchIndex, actives_temp);
            // Add contributions to the system matrix and right-hand side
            system_t.at(i).push(localMat_tilde.at(i), localRhs_tilde.at(i), actives_temp, eliminatedDofs[0], 0, 0);
        }
        for (unsigned i = 0; i < system_b.size(); i++) // n_tilde oder n_bar
        {
            // Map patch-local DoFs to global DoFs
            system_b.at(i).mapColIndices(actives, patchIndex, actives_temp);

            // Add contributions to the system matrix and right-hand side
            system_b.at(i).push(localMat_bar.at(i), localRhs_bar.at(i), actives_temp, eliminatedDofs[0], 0, 0);
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

    index_t n_tilde, n_bar;

}; // class gsVisitorG1Basis
} // namespace gismo