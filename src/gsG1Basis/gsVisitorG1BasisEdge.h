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
                         bool & isBoundary,
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

        gsMatrix<T> alpha, beta,
                    N_0, N_1,
                    N_j_minus, N_i_plus,
                    der_N_i_plus;

        if (geo.id() == 0) // Patch 0
        {
            if (optionList.getSwitch("direct"))
            {
                gluingData.eval_into_alpha_0(md.points,alpha);
                gluingData.eval_into_beta_0(md.points,beta);
            }
            else
            {
                gluingData.get_alpha_tilde_0().eval_into(md.points.bottomRows(1),alpha); // v
                gluingData.get_beta_tilde_0().eval_into(md.points.bottomRows(1),beta);
            }
            basis_geo.evalSingle_into(0,md.points.topRows(1),N_0); // u
            basis_geo.evalSingle_into(1,md.points.topRows(1),N_1); // u

            // Initialize local matrix/rhs
            for (index_t i = 0; i < n_plus; i++)
            {
                basis_plus.evalSingle_into(i,md.points.bottomRows(1),N_i_plus); // v
                basis_plus.derivSingle_into(i,md.points.bottomRows(1),der_N_i_plus);

                if (optionList.getSwitch("local"))
                {

                }

                beta = isBoundary ? beta.setZero() : beta; // For the boundary, only on Patch 0

                gsMatrix<T> temp = beta.cwiseProduct(N_1);
                rhsVals_tilde.at(i) = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(der_N_i_plus) * tau_1 / p;


                localMat_tilde.at(i).setZero(numActive, numActive);
                localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

                localMat_tilde.at(i).setZero(numActive, numActive);
                localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

            } // n_plus

            // Initialize local matrix/rhs
            for (index_t i = 0; i < n_minus; i++)
            {

                basis_minus.evalSingle_into(i,md.points.bottomRows(1),N_j_minus); // v

                if (optionList.getSwitch("local"))
                {

                }

                alpha = isBoundary ? alpha.setOnes() : alpha; // For the boundary, only on Patch 0

                rhsVals_bar.at(i) = alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));

                localMat_bar.at(i).setZero(numActive, numActive);
                localRhs_bar.at(i).setZero(numActive, rhsVals_bar.at(i).rows());//multiple right-hand sides

                localMat_bar.at(i).setZero(numActive, numActive);
                localRhs_bar.at(i).setZero(numActive, rhsVals_bar.at(i).rows());//multiple right-hand sides
            } // n_minus

        } // Patch 0
        else if (geo.id() == 1) // Patch 1
        {
            if (optionList.getSwitch("direct"))
            {
                gluingData.eval_into_alpha_1(md.points,alpha);
                gluingData.eval_into_beta_1(md.points,beta);
            }
            else
            {
                gluingData.get_alpha_tilde_1().eval_into(md.points.topRows(1),alpha); // u
                gluingData.get_beta_tilde_1().eval_into(md.points.topRows(1),beta);
            }

            basis_geo.evalSingle_into(0,md.points.bottomRows(1),N_0); // v
            basis_geo.evalSingle_into(1,md.points.bottomRows(1),N_1); // v

            // Initialize local matrix/rhs
            for (index_t i = 0; i < n_plus; i++)
            {
                basis_plus.evalSingle_into(i,md.points.topRows(1),N_i_plus); // u
                basis_plus.derivSingle_into(i,md.points.topRows(1),der_N_i_plus);

                if (optionList.getSwitch("local"))
                {

                }

                gsMatrix<> temp = beta.cwiseProduct(N_1);
                rhsVals_tilde.at(i) = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(der_N_i_plus) * tau_1 / p;

                localMat_tilde.at(i).setZero(numActive, numActive);
                localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

                localMat_tilde.at(i).setZero(numActive, numActive);
                localRhs_tilde.at(i).setZero(numActive, rhsVals_tilde.at(i).rows());//multiple right-hand sides

            } // n_tilde

            // Initialize local matrix/rhs
            for (index_t i = 0; i < n_minus; i++)
            {

                basis_minus.evalSingle_into(i,md.points.topRows(1),N_j_minus); // u

                if (optionList.getSwitch("local"))
                {

                }

                rhsVals_bar.at(i) = - alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));

                localMat_bar.at(i).setZero(numActive, numActive);
                localRhs_bar.at(i).setZero(numActive, rhsVals_bar.at(i).rows());//multiple right-hand sides

                localMat_bar.at(i).setZero(numActive, numActive);
                localRhs_bar.at(i).setZero(numActive, rhsVals_bar.at(i).rows());//multiple right-hand sides
            } // n_bar

        } // Patch 1
    } // evaluate

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