/** @file gsVisitorGluingData.h

    @brief Visitor for the G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

namespace gismo
{


template <class T>
class gsVisitorGlobalGD
{
public:

    gsVisitorGlobalGD()
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
        // md.flags = NEED_MEASURE ;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T>       & basis, //
                         gsMatrix<T>      & quNodes,
                         gsMultiPatch<T> & mp,
                         index_t & gamma)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        basis.eval_into(md.points,basisData);

        numActive = actives.rows();

        // ++++++++++++++++++++++++++++++++
        // Compute alpha^S and beta^S exact
        // ++++++++++++++++++++++++++++++++

        // alpha^S
        gsMatrix<> uv0, uv1, ev0, ev1;

        uv0.setZero(2,md.points.cols());
        uv0.topRows(1) = md.points;

        uv1.setZero(2,md.points.cols());
        uv1.bottomRows(1) = md.points;

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & P0 = mp.patch(0); // iFace.second().patch = 0

        for (index_t i = 0; i < uv0.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i), ev0);
            uv0(0, i) = gamma * ev0.determinant();
        }
        rhsVals_alpha_L = uv0.row(0);

        // ======== Determine bar{alpha^(R)} == Patch 1 ========
        const gsGeometry<> & P1 = mp.patch(1); // iFace.first().patch = 1

        for (index_t i = 0; i < uv1.cols(); i++)
        {
            P1.jacobian_into(uv1.col(i), ev1);
            uv1(0, i) = gamma * ev1.determinant(); // erste spalte: alphaL an stelle zweiter spalte
        }
        rhsVals_alpha_R = uv1.row(0);

        // beta^S
        uv0.setZero(2,md.points.cols());
        uv0.topRows(1) = md.points;

        uv1.setZero(2,md.points.cols());
        uv1.bottomRows(1) = md.points;

        const index_t d = mp.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        for(index_t i = 0; i < uv0.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i),ev0);
            D0 = ev0.col(0);
            real_t D1 = 1/ D0.norm();
            uv0(0,i) = gamma * D1 * D1 * ev0.col(1).transpose() * ev0.col(0);
        }
        rhsVals_beta_L = uv0.row(0);

        // ======== Determine bar{beta}^R ========
        for(index_t i = 0; i < uv1.cols(); i++)
        {
            P1.jacobian_into(uv1.col(i),ev1);
            D0 = ev1.col(1);
            real_t D1 = 1/ D0.norm();
            uv1(0,i) = gamma * D1 * D1 * ev1.col(0).transpose() * ev1.col(1);
        }
        rhsVals_beta_R = uv1.row(0);

        // ++++++++++++++++++++++++++++++++
        // ================================
        // ++++++++++++++++++++++++++++++++

        // Initialize local matrix/rhs
        localMat_L.setZero(numActive, numActive      );
        localRhs_L.setZero(numActive, rhsVals_alpha_L.rows() );//multiple right-hand sides

        localMat_R.setZero(numActive, numActive      );
        localRhs_R.setZero(numActive, rhsVals_alpha_R.rows() );//multiple right-hand sides

        localMat_b_L.setZero(numActive, numActive      );
        localRhs_b_L.setZero(numActive, rhsVals_beta_L.rows() );//multiple right-hand sides

        localMat_b_R.setZero(numActive, numActive      );
        localRhs_b_R.setZero(numActive, rhsVals_beta_R.rows() );//multiple right-hand sides

    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;

        // ( u, v)
        localMat_L.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();

        localMat_R.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose(); // Jacobi == 1

        localMat_b_L.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();

        localMat_b_R.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            localRhs_L.noalias() += weight * (basisVals.col(k) * rhsVals_alpha_L.col(k).transpose());
            localRhs_R.noalias() += weight * (basisVals.col(k) * rhsVals_alpha_R.col(k).transpose());

            localRhs_b_L.noalias() += weight * (basisVals.col(k) * rhsVals_beta_L.col(k).transpose());
            localRhs_b_R.noalias() += weight * (basisVals.col(k) * rhsVals_beta_R.col(k).transpose());
        }

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system_alpha_L,
                              gsSparseSystem<T>     & system_alpha_R,
                              gsSparseSystem<T>     & system_beta_L,
                              gsSparseSystem<T>     & system_beta_R)
    {
        gsMatrix<unsigned> actives_temp;

        // Map patch-local DoFs to global DoFs
        system_alpha_L.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system_alpha_L.push(localMat_L, localRhs_L, actives_temp, eliminatedDofs[0], 0, 0);

        // Map patch-local DoFs to global DoFs
        system_alpha_R.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system_alpha_R.push(localMat_R, localRhs_R, actives_temp, eliminatedDofs[1], 0, 0);

        // Map patch-local DoFs to global DoFs
        system_beta_L.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system_beta_L.push(localMat_b_L, localRhs_b_L, actives_temp, eliminatedDofs[2], 0, 0);

        // Map patch-local DoFs to global DoFs
        system_beta_R.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system_beta_R.push(localMat_b_R, localRhs_b_R, actives_temp, eliminatedDofs[3], 0, 0);


    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals_alpha_L, rhsVals_alpha_R;
    gsMatrix<T>  rhsVals_beta_L, rhsVals_beta_R;

protected:
    // Local matrices
    gsMatrix<T> localMat_L;
    gsMatrix<T>  localRhs_L;

    gsMatrix<T>  localMat_R;
    gsMatrix<T>  localRhs_R;

    gsMatrix<T> localMat_b_L;
    gsMatrix<T>  localRhs_b_L;

    gsMatrix<T>  localMat_b_R;
    gsMatrix<T>  localRhs_b_R;

    gsMapData<T> md;


}; // class gsVisitorGluingData

} // namespace gismo
