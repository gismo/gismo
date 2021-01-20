//
// Created by afarahat on 3/25/20.
//

#pragma once
#include <gismo.h>

namespace gismo
{


template <class T>
class gsG1ASGluingDataVisitorGlobal
{
public:

    gsG1ASGluingDataVisitorGlobal()
    {
    }

    void initialize(const gsBasis<T> & basis,
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
    inline void evaluate(gsMatrix<T>  & quNodes,
                         gsMultiPatch<T> & mp)
    {
        md.points = quNodes;
        actives.setZero(7, 1);
        actives << 0, 1, 2, 3, 4, 5, 6;
        numActive = actives.rows();

        basisData.setZero(numActive * numActive, md.points.cols());
        rhsVals.setZero(numActive, md.points.cols());

        gsGeometry<> & FR = mp.patch(0);
        gsGeometry<> & FL = mp.patch(1);

            gsMatrix<> pointV(FR.parDim(), md.points.cols());
            pointV.setZero();
            pointV.row(1) = md.points;

            gsMatrix<> pointU(FL.parDim(), md.points.cols());
            pointU.setZero();
            pointU.row(0) = md.points;

            gsMatrix<> ones(1, md.points.cols());
            ones.setOnes();

            gsMatrix<> lam = ones / 10000000000;

            gsMatrix<> DuFR(FR.targetDim(), md.points.cols());
            gsMatrix<> DvFR(FR.targetDim(), md.points.cols());
            gsMatrix<> DvFL(FR.targetDim(), md.points.cols());

            gsMatrix<> DuFR_DuFR(1, md.points.cols());
            gsMatrix<> DvFL_DvFL(1, md.points.cols());
            gsMatrix<> DvFR_DvFR(1, md.points.cols());

            gsMatrix<> DvFL_DuFR(1, md.points.cols());
            gsMatrix<> DvFR_DuFR(1, md.points.cols());
            gsMatrix<> DvFL_DvFR(1, md.points.cols());

        for(index_t i = 0; i < md.points.cols(); i++)
        {
            DuFR.col(i) = FR.jacobian(pointV.col(i)).col(0);
            DvFR.col(i) = FR.jacobian(pointV.col(i)).col(1); // Same as DuFL
            DvFL.col(i) = FL.jacobian(pointU.col(i)).col(1);


            // Set scalar product of the jacobian vectors of the geometric mapping
            DuFR_DuFR.col(i) = DuFR.col(i).transpose() * DuFR.col(i);
            DvFL_DvFL.col(i) = DvFL.col(i).transpose() * DvFL.col(i);
            DvFR_DvFR.col(i) = DvFR.col(i).transpose() * DvFR.col(i);

            DvFL_DuFR.col(i) = DvFL.col(i).transpose() * DuFR.col(i);
            DvFR_DuFR.col(i) = DvFR.col(i).transpose() * DuFR.col(i);
            DvFL_DvFR.col(i) = DvFL.col(i).transpose() * DvFR.col(i);
        }



        // Set Matrix 7x7 values
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // alpha_0R


        basisData.row(0) = 2 * ( md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DvFL - lam) - 2 * md.points.cwiseProduct(DvFL_DvFL - lam) + ones.cwiseProduct(DvFL_DvFL - lam));

        basisData.row(1) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFL - lam);
        basisData.row(7) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFL - lam);

        basisData.row(2) = 2 * (md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DuFR) - 2 * md.points.cwiseProduct(DvFL_DuFR) + ones.cwiseProduct(DvFL_DuFR));
        basisData.row(14) = 2 * (md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DuFR) - 2 * md.points.cwiseProduct(DvFL_DuFR) + ones.cwiseProduct(DvFL_DuFR));

        basisData.row(3) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DvFL_DuFR);
        basisData.row(21) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DvFL_DuFR);

        basisData.row(4) = 2 * (ones.cwiseProduct(DvFL_DvFR) - 3 * md.points.cwiseProduct(DvFL_DvFR) + 3 * md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DvFR)
            - md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFL_DvFR));
        basisData.row(28) = 2 * (ones.cwiseProduct(DvFL_DvFR) - 3 * md.points.cwiseProduct(DvFL_DvFR) + 3 * md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DvFR)
            - md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFL_DvFR));

        basisData.row(5) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);
        basisData.row(35) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);

        basisData.row(6) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);
        basisData.row(42) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // alpha_1R

        basisData.row(8) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DvFL - lam);

        basisData.row(9) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DvFL_DuFR);
        basisData.row(15) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DvFL_DuFR);

        basisData.row(10) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DuFR);
        basisData.row(22) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(DvFL_DuFR);

        basisData.row(11) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);
        basisData.row(29) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);

        basisData.row(12) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);
        basisData.row(36) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFL_DvFR);

        basisData.row(13) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFL_DvFR);
        basisData.row(43) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFL_DvFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // alpha_0L

        basisData.row(16) = 2 * (md.points.cwiseProduct(md.points).cwiseProduct(DuFR_DuFR - lam) - 2 * md.points.cwiseProduct(DuFR_DuFR - lam) + ones.cwiseProduct(DuFR_DuFR - lam));

        basisData.row(17) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DuFR_DuFR - lam);
        basisData.row(23) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(DuFR_DuFR - lam);

        basisData.row(18) = 2 * (ones.cwiseProduct(DvFR_DuFR) - 3 * md.points.cwiseProduct(DvFR_DuFR) + 3 * md.points.cwiseProduct(md.points).cwiseProduct(DvFR_DuFR)
            - md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DuFR));
        basisData.row(30) = 2 * (ones.cwiseProduct(DvFR_DuFR) - 3 * md.points.cwiseProduct(DvFR_DuFR) + 3 * md.points.cwiseProduct(md.points).cwiseProduct(DvFR_DuFR)
            - md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DuFR));

        basisData.row(19) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);
        basisData.row(37) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);

        basisData.row(20) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);
        basisData.row(44) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // alpha_1L

        basisData.row(24) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(DuFR_DuFR - lam);

        basisData.row(25) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);
        basisData.row(31) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);

        basisData.row(26) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);
        basisData.row(38) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DuFR);

        basisData.row(27) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DuFR);
        basisData.row(45) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DuFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_0

        basisData.row(32) = 2 * ( ones.cwiseProduct(DvFR_DvFR) - 4 * md.points.cwiseProduct(DvFR_DvFR) + 6 * md.points.cwiseProduct(md.points).cwiseProduct(DvFR_DvFR)
                            - 4 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DvFR)
                            + md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DvFR));

        basisData.row(33) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);
        basisData.row(39) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);

        basisData.row(34) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);
        basisData.row(46) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_1

        basisData.row(40) = 8 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);

        basisData.row(41) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);
        basisData.row(47) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(DvFR_DvFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_2

        basisData.row(48) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(DvFR_DvFR);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Set RHS 7x1 values

        rhsVals.row(0) = -2 * ( lam - md.points.cwiseProduct(lam) );

        rhsVals.row(1) = -2 * md.points.cwiseProduct(lam);

        rhsVals.row(2) = -2 * ( lam - md.points.cwiseProduct(lam) );

        rhsVals.row(3) = -2 * md.points.cwiseProduct(lam);

        rhsVals.row(4) = 0 * ones;

        rhsVals.row(5) = 0 * ones;

        rhsVals.row(6) = 0 * ones;

        // Initialize local matrix/rhs
            localMat.setZero(numActive, numActive);
            localRhs.setZero(numActive, 1 ); //multiple right-hand sides
    }


    // Evaluate on element.
    inline void evaluateBeta(gsMatrix<T>  & quNodes,
                         gsMultiPatch<T> & mp,
                             gsMatrix<> sol)
    {
        md.points = quNodes;
        activesBeta.setZero(4, 1);
        activesBeta << 0, 1, 2, 3;
        numActiveBeta = activesBeta.rows();

        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        gsMatrix<> alpha_R = sol.row(0) * ( ones - md.points ) + sol.row(1) * md.points;
        gsMatrix<> alpha_L = sol.row(2) * ( ones - md.points ) + sol.row(3) * md.points;
        gsMatrix<> beta = sol.row(4) * ( md.points.cwiseProduct(md.points) - 2 * md.points + ones ) + 2 * sol.row(5) * md.points.cwiseProduct( ones - md.points ) + sol.row(6) * md.points.cwiseProduct(md.points);


        gsMatrix<> alpha_R_Squared = alpha_R.cwiseProduct(alpha_R);
        gsMatrix<> alpha_L_Squared = alpha_L.cwiseProduct(alpha_L);

        gsMatrix<> alpha_R_L = alpha_R.cwiseProduct(alpha_L);

        gsMatrix<> lamB = ones / 10000000;


        basisDataBeta.setZero(numActiveBeta * numActiveBeta, md.points.cols());
        rhsValsBeta.setZero(numActiveBeta, md.points.cols());


        // Set Matrix 4x4 values
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_0R

        basisDataBeta.row(0) = 2 * ( md.points.cwiseProduct(md.points).cwiseProduct(alpha_L_Squared - lamB) - 2 * md.points.cwiseProduct(alpha_L_Squared - lamB) + ones.cwiseProduct(alpha_L_Squared - lamB));

        basisDataBeta.row(1) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_L_Squared - lamB);
        basisDataBeta.row(4) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_L_Squared - lamB);

        basisDataBeta.row(2) = 2 * ( md.points.cwiseProduct(md.points).cwiseProduct(alpha_R_L) - 2 * md.points.cwiseProduct(alpha_R_L) + ones.cwiseProduct(alpha_R_L));
        basisDataBeta.row(8) = 2 * ( md.points.cwiseProduct(md.points).cwiseProduct(alpha_R_L) - 2 * md.points.cwiseProduct(alpha_R_L) + ones.cwiseProduct(alpha_R_L));

        basisDataBeta.row(3) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_R_L);
        basisDataBeta.row(12) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_R_L);
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_1R

        basisDataBeta.row(5) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(alpha_L_Squared - lamB);

        basisDataBeta.row(6) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_R_L);
        basisDataBeta.row(9) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_R_L);

        basisDataBeta.row(7) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(alpha_R_L);
        basisDataBeta.row(13) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(alpha_R_L);


// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_0L

        basisDataBeta.row(10) = 2 * ( md.points.cwiseProduct(md.points).cwiseProduct(alpha_R_Squared - lamB) - 2 * md.points.cwiseProduct(alpha_R_Squared - lamB) + ones.cwiseProduct(alpha_R_Squared - lamB));

        basisDataBeta.row(11) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_R_Squared - lamB);
        basisDataBeta.row(14) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(alpha_R_Squared - lamB);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        // beta_1L

        basisDataBeta.row(15) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(alpha_R_Squared - lamB);

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // Set RHS 4x1 values

        rhsValsBeta.row(0) = 2 * beta.cwiseProduct(ones - md.points).cwiseProduct(alpha_L);

        rhsValsBeta.row(1) = 2 * beta.cwiseProduct(md.points).cwiseProduct(alpha_L);

        rhsValsBeta.row(2) = 2 * beta.cwiseProduct(ones - md.points).cwiseProduct(alpha_R);

        rhsValsBeta.row(3) = 2 * beta.cwiseProduct(md.points).cwiseProduct(alpha_R);

        // Initialize local matrix/rhs for beta
        localMatBeta.setZero(numActiveBeta, numActiveBeta);
        localRhsBeta.setZero(numActiveBeta, 1 ); //multiple right-hand sides
    }


    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            for(index_t i = 0; i < numActive; i++)
            {
                localRhs(i, 0) += weight * rhsVals(i,k);

                for(index_t j = 0; j < numActive; j++)
                {
                    localMat(i, j) += weight * basisData(i*7 + j, k);
                }
            }
        }
//        gsInfo << "Matrix : " << localMat << "\n";
    }

    inline void assembleBeta(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            for(index_t i = 0; i < numActiveBeta; i++)
            {
                localRhsBeta(i, 0) += weight * rhsValsBeta(i,k);

                for(index_t j = 0; j < numActiveBeta; j++)
                {
                    localMatBeta(i, j) += weight * basisDataBeta(i*numActiveBeta + j, k);
                }
            }
        }
    }


    inline void localToGlobal(const gsMatrix<T>    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        gsMatrix<unsigned> actives_temp;

        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, 0, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives_temp, eliminatedDofs, 0, 0);

    }


    inline void localToGlobalBeta(const gsMatrix<T>    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        gsMatrix<unsigned> actives_temp;

        // Map patch-local DoFs to global DoFs
        system.mapColIndices(activesBeta, 0, actives_temp);
        // Add contributions to the system matrix and right-hand side
        system.push(localMatBeta, localRhsBeta, actives_temp, eliminatedDofs, 0, 0);

    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<unsigned> activesBeta;
    gsMatrix<T> basisData;
    gsMatrix<T> basisDataBeta;
    index_t numActive;
    index_t numActiveBeta;


protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals;
    gsMatrix<T>  rhsValsBeta;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    // Local matrices
    gsMatrix<T> localMatBeta;
    gsMatrix<T>  localRhsBeta;

    gsMapData<T> md;

}; // class gsVisitorGluingData

} // namespace gismo

