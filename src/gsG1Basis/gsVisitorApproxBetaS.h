/** @file gsVisitorApproxBetaS.h

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
class gsVisitorApproxBetaS
{
public:

    gsVisitorApproxBetaS()
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
    inline void evaluate(gsBSplineBasis<>  & basis, //
                         gsMatrix<T>      & quNodes,
                         gsMultiPatch<T> & mp,
                         gsG1ASGluingData<T> gd_andrea)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);


        gsMatrix<unsigned > actives2(actives.rows(),1);
        for (index_t i = 0; i < actives2.rows(); i++)
            actives2(i,0) = actives(i,0) + basis.size();

        actives.conservativeResize(2*actives.rows(),actives.cols());
        actives.block(actives.rows()/2,0,actives.rows()/2,1) = actives2;



        // Evaluate basis functions on element
        gsMatrix<> basisData_Value;
        basis.eval_into(md.points,basisData);

        numActive = actives.rows();

        // ++++++++++++++++++++++++++++++++
        // Compute alpha^S and beta^S exact
        // ++++++++++++++++++++++++++++++++

        // alpha^S
        gsMatrix<> uv, ev, alpha_L, alpha_R; // LEFT == u and patch 1

        uv.setZero(2,md.points.cols());
        uv.topRows(1) = md.points; // u

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & PL = mp.patch(1); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            PL.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();

        }
        alpha_L = uv.row(0);

        // RIGHT == v and patch 0
        uv.setZero(2,md.points.cols());
        uv.bottomRows(1) = md.points; // v


        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & PR = mp.patch(0); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            PR.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();

        }
        alpha_R = uv.row(0);

        // BETA
        gsMatrix<> beta;
        gsMatrix<> uv1, uv0, ev1, ev0;

        const index_t d = 2;
        gsMatrix<> D0(d,d);

        gsGeometry<>::Ptr beta_temp;


        uv0.setZero(2,md.points.cols());
        uv0.bottomRows(1) = md.points; // v

        uv1.setZero(2,md.points.cols());
        uv1.topRows(1) = md.points; // u


        const gsGeometry<> & P0 = mp.patch(0); // iFace.first().patch = 1
        const gsGeometry<> & P1 = mp.patch(1); // iFace.second().patch = 0
        // ======================================

        // ======== Determine bar{beta} ========
        for(index_t i = 0; i < uv1.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i),ev0);
            P1.jacobian_into(uv1.col(i),ev1);

            D0.col(1) = ev0.col(0); // (DuFL, *)
            D0.col(0) = ev1.col(1); // (*,DuFR)

            uv0(0,i) = D0.determinant();
        }

        beta = uv0.row(0);

        /*
#pragma omp critical
        {
            gsMatrix<> sol_andrea(7,1);
            sol_andrea << 1.4337782823995119940718723228, 0.8300821382872851650347456598, 0.45277209114169275627759247982, 0.98100616960476005878888372536, 0.56596511177547570436985324704, 0.13205851945361074539775358971, -0.37731006523260224305715837545;

            gsMatrix<> ones(1, md.points.cols());
            ones.setOnes();
            //alpha_R = sol_andrea.row(0) * ( ones - md.points ) + sol_andrea.row(1) * md.points;
            //alpha_L = sol_andrea.row(2) * ( ones - md.points ) + sol_andrea.row(3) * md.points;

            gsMatrix<> beta_andrea = sol_andrea.row(4) * ( md.points.cwiseProduct(md.points) - 2 * md.points + ones ) + 2 * sol_andrea.row(5) * md.points.cwiseProduct( ones - md.points ) + sol_andrea.row(6) * md.points.cwiseProduct(md.points);
            //beta = beta_andrea;

        }
*/
        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        //alpha_L = gd_andrea.evalAlpha_L(md.points);
        //alpha_R = gd_andrea.evalAlpha_R(md.points);
        gsMatrix<> sol = gd_andrea.getSol();
        //beta = sol.row(4) * ( md.points.cwiseProduct(md.points) - 2 * md.points + ones ) + 2 * sol.row(5) * md.points.cwiseProduct( ones - md.points ) + sol.row(6) * md.points.cwiseProduct(md.points);


        real_t lambda = 1/1000000000;

        rhsVals1 = alpha_L.cwiseProduct(beta);
        rhsVals2 = alpha_R.cwiseProduct(beta);





        block11 = alpha_L.cwiseProduct(alpha_L) + lambda * ones;
        block22 = alpha_R.cwiseProduct(alpha_R) + lambda * ones;
        //block11 = alpha_L.cwiseProduct(alpha_L);
        //block22 = alpha_R.cwiseProduct(alpha_R);

        block12 = alpha_L.cwiseProduct(alpha_R);
        block21 = alpha_L.cwiseProduct(alpha_R);

//#pragma omp critical
//        {
//        }
        // ++++++++++++++++++++++++++++++++
        // ================================
        // ++++++++++++++++++++++++++++++++

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals1.rows() );//multiple right-hand sides


    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            localRhs.block(0, 0, actives.rows() / 2, rhsVals1.rows()).noalias() +=
                weight * (basisVals.col(k) * rhsVals1.col(k).transpose());
            localRhs.block(actives.rows() / 2, 0, actives.rows() / 2, rhsVals2.rows()).noalias() +=
                weight * (basisVals.col(k) * rhsVals2.col(k).transpose());


            localMat.block(0, 0, actives.rows() / 2, actives.rows() / 2).noalias() +=
                weight * block11(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(actives.rows() / 2, actives.rows() / 2, actives.rows() / 2, actives.rows() / 2).noalias() +=
                weight * block22(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(0, actives.rows() / 2, actives.rows() / 2, actives.rows() / 2).noalias() +=
                weight * block12(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(actives.rows() / 2, 0, actives.rows() / 2, actives.rows() / 2).noalias() +=
                weight * block21(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        gsMatrix<unsigned> actives_temp;

        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives_temp, eliminatedDofs[0], 0, 0);
    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
    index_t numActive;

    gsMatrix<T> block11, block12, block21, block22;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals1, rhsVals2;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    gsMapData<T> md;


}; // class gsVisitorGluingData

} // namespace gismo
