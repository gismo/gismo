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
                         gsMultiPatch<T> & mp)
    {
        md.points = quNodes;

        const boundaryInterface & item = mp.interfaces()[0];
        index_t uv_L = item.first().index() < 3 ? 1 : 0;
        index_t L = item.first().patch;

        index_t uv_R = item.second().index() < 3 ? 1 : 0;
        index_t R = item.second().patch;

        gsInfo << "L:" << L << ": R : " << R << "\n";
        gsInfo << "L:" << uv_L << ": R : " << uv_R << "\n";

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
        gsMatrix<> uv, ev, alpha_L, alpha_R;

        if (uv_L==1)
        {
            uv.setZero(2,md.points.cols()); // TODO
            uv.bottomRows(1) = md.points; // v
        }
        else if (uv_L==0)
        {
            uv.setZero(2,md.points.cols());
            uv.topRows(1) = md.points; // u
        }

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & PL = mp.patch(L); // iFace.second().patch = 0

        for (index_t i = 0; i < uv.cols(); i++)
        {
            PL.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();

        }
        alpha_L = uv.row(0);

        if (uv_R==1)
        {
            uv.setZero(2,md.points.cols());
            uv.bottomRows(1) = md.points; // v
        }
        else if (uv_R==0)
        {
            uv.setZero(2,md.points.cols());
            uv.topRows(1) = md.points; // u
        }

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & PR = mp.patch(R); // iFace.second().patch = 0

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

        if (uv_R==1)
        {
            uv0.setZero(2,md.points.cols());
            uv0.bottomRows(1) = md.points; // v
        }
        else if (uv_R==0)
        {
            uv0.setZero(2,md.points.cols());
            uv0.topRows(1) = md.points; // u
        }

        if (uv_L==1)
        {
            uv1.setZero(2,md.points.cols());
            uv1.bottomRows(1) = md.points; // v
        }
        else if (uv_L==0)
        {
            uv1.setZero(2,md.points.cols());
            uv1.topRows(1) = md.points; // u
        }

        const gsGeometry<> & P0 = mp.patch(R); // iFace.first().patch = 1
        const gsGeometry<> & P1 = mp.patch(L); // iFace.second().patch = 0
        // ======================================

        // ======== Determine bar{beta} ========
        for(index_t i = 0; i < uv1.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i),ev0);
            P1.jacobian_into(uv1.col(i),ev1);

            D0.col(uv_R) = ev0.col(1-uv_R); // (DuFL, *)
            D0.col(uv_L) = ev1.col(1-uv_L); // (*,DuFR)

            uv0(0,i) = D0.determinant();
        }

        beta = uv0.row(0);

        real_t lambda = 1/100000;

        rhsVals1 = alpha_R.cwiseProduct(beta);
        rhsVals2 = alpha_L.cwiseProduct(beta);

        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        block11 = lambda * ones + alpha_R.cwiseProduct(alpha_R);
        block22 = lambda * ones - alpha_L.cwiseProduct(alpha_L);

        block12 = alpha_L.cwiseProduct(alpha_R);
        block21 = block12;

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

/*
        localMat.block(actives.rows()/2,actives.rows()/2,actives.rows()/2,actives.rows()/2) =
            basisData * quWeights.asDiagonal() * basisData.transpose();

        gsInfo << "actives 2: " << localMat << "\n";


        localMat.block(0,actives.rows()/2,actives.rows()/2,actives.rows()/2).noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();

        localMat.block(actives.rows()/2,0,actives.rows()/2,actives.rows()/2).noalias() =
            basisData * quWeights.asDiagonal() * block21.asDiagonal() * basisData.transpose();
*/

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            localRhs.block(0,0,actives.rows()/2,rhsVals1.rows()).noalias() += weight * (basisVals.col(k) * rhsVals1.col(k).transpose());
            localRhs.block(actives.rows()/2,0,actives.rows()/2,rhsVals1.rows()).noalias() += weight * (basisVals.col(k) * rhsVals2.col(k).transpose());

            localMat.block(0,0,actives.rows()/2,actives.rows()/2).noalias() += weight * block11(0,k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(actives.rows()/2,actives.rows()/2,actives.rows()/2,actives.rows()/2).noalias() += weight * block22(0,k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(0,actives.rows()/2,actives.rows()/2,actives.rows()/2).noalias() += weight * block12(0,k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(actives.rows()/2,0,actives.rows()/2,actives.rows()/2).noalias() += weight * block21(0,k) * (basisVals.col(k) * basisVals.col(k).transpose());
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

    gsVector<T> block11, block12, block21, block22;

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
