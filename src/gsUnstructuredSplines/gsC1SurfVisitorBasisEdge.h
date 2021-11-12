/** @file gsC1SurfEdge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsC1SurfGluingData.h>

namespace gismo
{
    template <class T>
    class gsC1SurfVisitorBasisEdge
    {
    public:

        gsC1SurfVisitorBasisEdge()
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
                             gsC1SurfGluingData<T>  & gluingData,
                             bool & isBoundary)
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

            gsMatrix<T> alpha, beta,
                    N_0, N_1,
                    N_j_minus, N_i_plus,
                    der_N_i_plus;


            if (uv == 1) // edge is in v-direction
            {

                alpha = gluingData.evalAlpha_R(md.points.bottomRows(1));
                beta = gluingData.evalBeta_R(md.points.bottomRows(1));

                basis_geo.evalSingle_into(0,md.points.topRows(1),N_0); // u
                basis_geo.evalSingle_into(1,md.points.topRows(1),N_1); // u

                // Initialize local matrix/rhs
                if (typeBf == "plus")
                {
                    basis_plus.evalSingle_into(bfID,md.points.bottomRows(1),N_i_plus); // v
                    basis_plus.derivSingle_into(bfID,md.points.bottomRows(1),der_N_i_plus);


                    beta = isBoundary ? beta.setZero() : beta; // For the boundary, only on Patch 0

                    gsMatrix<T> temp = beta.cwiseProduct(N_1);
                    rhsVals = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(der_N_i_plus) * tau_1 / p;

                    localMat.setZero(numActive, numActive);
                    localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides

                } // n_plus
                else if (typeBf == "minus")
                {
                    basis_minus.evalSingle_into(bfID,md.points.bottomRows(1),N_j_minus); // v


                    alpha = isBoundary ? alpha.setOnes() : alpha; // For the boundary, only on Patch 0

                    rhsVals = alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));

                    localMat.setZero(numActive, numActive);
                    localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
                } // n_minus

            } // Patch 0
            else if (uv == 0) // edge is in u-direction
            {

                alpha = gluingData.evalAlpha_L(md.points.topRows(1));
                beta = gluingData.evalBeta_L(md.points.topRows(1));

                basis_geo.evalSingle_into(0,md.points.bottomRows(1),N_0); // v
                basis_geo.evalSingle_into(1,md.points.bottomRows(1),N_1); // v

                // Initialize local matrix/rhs
                if (typeBf == "plus")
                {
                    basis_plus.evalSingle_into(bfID,md.points.topRows(1),N_i_plus); // u
                    basis_plus.derivSingle_into(bfID,md.points.topRows(1),der_N_i_plus);

                    beta = isBoundary ? beta.setZero() : beta; // For the boundary, only on Patch 0

                    gsMatrix<> temp = beta.cwiseProduct(N_1);
                    rhsVals = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(der_N_i_plus) * tau_1 / p;

                    localMat.setZero(numActive, numActive);
                    localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides

                } // n_tilde
                else if (typeBf == "minus")
                {
                    basis_minus.evalSingle_into(bfID,md.points.topRows(1),N_j_minus); // u


                    alpha = isBoundary ? alpha.setOnes() : alpha; // For the boundary, only on Patch 0

                    rhsVals = - alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));

                    localMat.setZero(numActive, numActive);
                    localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
                } // n_bar

            } // Patch 1
        } // evaluate1

        inline void assemble(gsDomainIterator<T>    & element,
                             const gsVector<T>      & quWeights)
        {
            gsMatrix<T> & basisVals  = basisData;

            // ( u, v)
            localMat.noalias() =
                    basisData * quWeights.asDiagonal() *
                    md.measures.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k] * md.measure(k);
                localRhs.noalias() += weight * (basisVals.col(k) * rhsVals.col(k).transpose());
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
        gsMatrix<index_t> actives;
        gsMatrix<T> basisData;
        index_t numActive;

    protected:
        // Local values of the right hand side
        gsMatrix<T>  rhsVals;

    protected:
        // Local matrices
        gsMatrix<T> localMat;
        gsMatrix<T> localRhs;

        gsMapData<T> md;

    }; // class gsVisitorG1BasisEdge
} // namespace gismo
