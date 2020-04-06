//
// Created by afarahat on 3/25/20.
//

#pragma once

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
        actives.setZero(6, 1);
        actives << 0, 1, 2, 3, 4, 5;
        numActive = actives.rows();

        basisData.setZero(36, md.points.cols());
        rhsVals.setZero(6, md.points.cols());

        gsGeometry<> & FR = mp.patch(0);
        gsGeometry<> & FL = mp.patch(1);

            gsMatrix<> pointV(FR.parDim(), md.points.cols());
            pointV.setZero();
            pointV.row(1) = md.points;

            gsMatrix<> pointU(FL.parDim(), md.points.cols());
            pointU.setZero();
            pointU.row(0) = md.points;

            DuFR = FR.jacobian(pointV).col(0);
            DvFR = FR.jacobian(pointV).col(1); // Same as DuFL

            DvFL = FL.jacobian(pointU).col(1);

            gsMatrix<> ones(1, md.points.cols());
            ones.setOnes();

            // Set Matrix 6x6 values

            basisData.row(0) = 2 * md.points.cwiseProduct(md.points) * DuFR.norm() * DuFR.norm();

            basisData.row(1) = 2 * md.points.cwiseProduct(ones - md.points) * DvFL.transpose() * DuFR;
            basisData.row(6) = 2 * md.points.cwiseProduct(ones - md.points) * DvFL.transpose() * DuFR;

            basisData.row(2) = 2 * md.points.cwiseProduct(md.points) * DvFL.transpose() * DuFR;
            basisData.row(12) = 2 * md.points.cwiseProduct(md.points) * DvFL.transpose() * DuFR;

            basisData.row(3) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DuFR.transpose() * DvFR;
            basisData.row(18) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DuFR.transpose() * DvFR;

            basisData.row(4) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFR.transpose() * DuFR;
            basisData.row(24) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFR.transpose() * DuFR;

            basisData.row(5) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points) * DvFR.transpose() * DuFR;
            basisData.row(30) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points) * DvFR.transpose() * DuFR;

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            basisData.row(7) = 2 * (md.points.cwiseProduct(md.points) - 2 * md.points + ones) * DvFL.norm() * DvFL.norm();

            basisData.row(8) = 2 * md.points.cwiseProduct(ones - md.points) * DvFL.norm() * DvFL.norm();
            basisData.row(13) = 2 * md.points.cwiseProduct(ones - md.points) * DvFL.norm() * DvFL.norm();

            basisData.row(9) = 2 * (ones - 3 * md.points + 3 * md.points.cwiseProduct(md.points) - md.points.cwiseProduct(md.points).cwiseProduct(md.points)) * DvFL.transpose() * DvFR;
            basisData.row(19) = 2 * (ones - 3 * md.points+ 3 * md.points.cwiseProduct(md.points) - md.points.cwiseProduct(md.points).cwiseProduct(md.points)) * DvFL.transpose() * DvFR;

            basisData.row(10) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;
            basisData.row(25) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;

            basisData.row(11) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;
            basisData.row(31) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            basisData.row(14) = 2 * md.points.cwiseProduct(md.points) * DvFL.norm() * DvFL.norm();

            basisData.row(15) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;
            basisData.row(20) = 2 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;

            basisData.row(16) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;
            basisData.row(26) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFL.transpose() * DvFR;

            basisData.row(17) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points) * DvFL.transpose() * DvFR;
            basisData.row(32) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points) * DvFL.transpose() * DvFR;

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            basisData.row(21) = 2 * ( ones - 4 * md.points + 6 * md.points.cwiseProduct(md.points) - 4 * md.points.cwiseProduct(md.points).cwiseProduct(md.points) +
                                md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(md.points) ) * DvFR.norm() * DvFR.norm();

            basisData.row(22) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();
            basisData.row(27) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();

            basisData.row(23) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();
            basisData.row(33) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            basisData.row(28) = 8 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();

            basisData.row(29) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();
            basisData.row(34) = 4 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFR.norm() * DvFR.norm();

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            basisData.row(35) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(md.points).cwiseProduct(md.points) * DvFR.norm() * DvFR.norm();


            // Set RHS 6x1 values

            rhsVals.row(0) = 2 * md.points.cwiseProduct(ones - md.points) * DuFR.norm() * DuFR.norm();

            rhsVals.row(1) = 2 * (md.points.cwiseProduct(md.points) - 2 * md.points + ones) * DvFL.transpose() * DuFR;

            rhsVals.row(2) = 2 * md.points.cwiseProduct(ones - md.points) * DvFL.transpose() * DuFR;

            rhsVals.row(3) = 2 * (ones - 3 * md.points + 3 * md.points.cwiseProduct(md.points) - md.points.cwiseProduct(md.points).cwiseProduct(md.points)) * DvFR.transpose() * DuFR;

            rhsVals.row(4) = 4 * md.points.cwiseProduct(ones - md.points).cwiseProduct(ones - md.points) * DvFR.transpose() * DuFR;

            rhsVals.row(5) = 2 * md.points.cwiseProduct(md.points).cwiseProduct(ones - md.points) * DvFR.transpose() * DuFR;


            // Initialize local matrix/rhs
            localMat.setZero(numActive, numActive);
            localRhs.setZero(numActive, rhsVals.rows() ); //multiple right-hand sides
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];


            for(index_t i = 0; i < 6; i++)
            {
                localRhs(i, 0) += weight * rhsVals(i,k);
                for(index_t j = i; j < 6; j++)
                {
                    localMat(i, j) += weight * basisData(i*6 + j, k);
                    localMat(j, i) += weight * basisData(i*6 + j, k);
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

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    gsMapData<T> md;

    gsMatrix<> DuFR;
    gsMatrix<> DvFR;

    gsMatrix<> DvFL;
}; // class gsVisitorGluingData

} // namespace gismo

