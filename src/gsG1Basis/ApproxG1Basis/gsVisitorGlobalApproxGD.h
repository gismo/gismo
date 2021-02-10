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
class gsVisitorGlobalApproxGD
{
public:

    gsVisitorGlobalApproxGD()
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
                         index_t m_uv,
                         gsMultiPatch<T> & mp,
                         index_t m_patchID,
                         real_t & gamma,
                         bool & isBoundary,
                         bool twoPatch)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        //basis.eval_into(md.points,basisData);
        basis.evalAllDers_into( md.points, 1, basisData);

        numActive = actives.rows();

        // ++++++++++++++++++++++++++++++++
        // Compute alpha^S and beta^S exact
        // ++++++++++++++++++++++++++++++++
        real_t D1, lambda0 = 0, lambda1 = 0;
        gsMatrix<> uv, ev, D0;

        gsMatrix<> alpha_S, beta_S;

        if (mp.nPatches() == 2)
        {
            gsMatrix<> zeroOne(2, 2);
            zeroOne.setZero();
            zeroOne(1, 1) = 1.0; // v

            const gsGeometry<> & PR = mp.patch(0); // Right m_uv = 1
            const gsGeometry<> & PL = mp.patch(1); // Left m_uv = 0

            PR.jacobian_into(zeroOne.col(1), ev);
            lambda1 = 1 / ev.determinant(); // alpha_R

            D0 = ev.col(1);
            D1 = 1 / D0.norm();
            lambda1 *= (m_uv == 1 ? -gamma : gamma) * D1 * D1 * ev.col(1).transpose() * ev.col(0);

            PL.jacobian_into(zeroOne.col(0), ev);
            lambda0 = 1 / ev.determinant(); // alpha_L

            D0 = ev.col(0);
            D1 = 1 / D0.norm();
            lambda0 *= (m_uv == 1 ? gamma : -gamma) * D1 * D1 * ev.col(1).transpose() * ev.col(0);
        } // nPatches == 2

        // ======== Determine bar{alpha}^S ========
        uv.setZero(2,md.points.cols());
        uv.row(m_uv) = md.points; // u

        const gsGeometry<> & P0 = mp.patch(m_patchID); // Right
        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();

        }
        if (isBoundary)
            uv.setOnes();
        alpha_S = uv.row(0);

        // ======== Determine bar{beta}^S ========
        uv.setZero(2,md.points.cols());
        uv.row(m_uv) = md.points; // u

        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(m_uv);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);
        }
        if (isBoundary)
            uv.setZero();
        beta_S = uv.row(0);


        // ++++++++++++++++++++++++++++++++
        // ================================
        // ++++++++++++++++++++++++++++++++
        gsMatrix<> ones;
        ones.setOnes(1,md.points.cols());

        rhsVals_alpha = alpha_S;
        if (twoPatch)
            rhsVals_beta = beta_S - lambda0 * (ones - md.points).cwiseProduct(alpha_S) - lambda1 * (md.points).cwiseProduct(alpha_S);
        else
            rhsVals_beta = beta_S;


        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals_alpha.rows() );//multiple right-hand sides


        localMat_b.setZero(numActive, numActive      );
        localRhs_b.setZero(numActive, rhsVals_beta.rows() );//multiple right-hand sides


    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData[0];


        // ( u, v)
        localMat.noalias() =
            basisVals * quWeights.asDiagonal() * basisVals.transpose();

        localMat_b.noalias() =
            basisVals * quWeights.asDiagonal() * basisVals.transpose();


        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            localRhs.noalias() += weight * (basisVals.col(k) * rhsVals_alpha.col(k).transpose());
            localRhs_b.noalias() += weight * (basisVals.col(k) * rhsVals_beta.col(k).transpose());

        }

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system_alpha_L,
                              gsSparseSystem<T>     & system_beta_L)
    {
        gsMatrix<unsigned> actives_temp;

        // Map patch-local DoFs to global DoFs
        system_alpha_L.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system_alpha_L.push(localMat, localRhs, actives_temp, eliminatedDofs[0], 0, 0);

        // Map patch-local DoFs to global DoFs
        system_beta_L.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system_beta_L.push(localMat_b, localRhs_b, actives_temp, eliminatedDofs[1], 0, 0);

    }

protected:
    gsMatrix<unsigned> actives;
    std::vector<gsMatrix<T>> basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals_alpha;
    gsMatrix<T>  rhsVals_beta;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    gsMatrix<T> localMat_b;
    gsMatrix<T>  localRhs_b;

    gsMapData<T> md;


}; // class gsVisitorGluingData

} // namespace gismo
