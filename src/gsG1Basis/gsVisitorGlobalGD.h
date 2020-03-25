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
                         index_t m_uv,
                         gsMultiPatch<T> & mp,
                         index_t & gamma,
                         bool & isBoundary)
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
        gsMatrix<> uv, ev;

        if (m_uv==1)
        {
            uv.setZero(2,md.points.cols());
            uv.bottomRows(1) = md.points; // v
        }
        else if (m_uv==0)
        {
            uv.setZero(2,md.points.cols());
            uv.topRows(1) = md.points; // u
        }


        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & P0 = mp.patch(0); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();
        }
        if (isBoundary)
            uv.setOnes();
        rhsVals_alpha = uv.row(0);


        // beta^S
        if (m_uv==1)
        {
            uv.setZero(2,md.points.cols());
            uv.bottomRows(1) = md.points; // v
        }
        else if (m_uv==0)
        {
            uv.setZero(2,md.points.cols());
            uv.topRows(1) = md.points; // u
        }

        const index_t d = mp.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(m_uv);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);

        }
        if (isBoundary)
            uv.setZero();
        rhsVals_beta = uv.row(0);


        // ++++++++++++++++++++++++++++++++
        // ================================
        // ++++++++++++++++++++++++++++++++

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals_alpha.rows() );//multiple right-hand sides


        localMat_b.setZero(numActive, numActive      );
        localRhs_b.setZero(numActive, rhsVals_beta.rows() );//multiple right-hand sides


    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;

        // ( u, v)
        localMat.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();



        localMat_b.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();


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
        system_beta_L.push(localMat_b, localRhs_b, actives_temp, eliminatedDofs[2], 0, 0);

    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
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
