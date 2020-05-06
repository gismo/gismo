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
class gsVisitorGlobalApproxGD3Lambda
{
public:

    gsVisitorGlobalApproxGD3Lambda()
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
    inline void evaluate(std::vector<gsBSplineBasis<T>>  & basis, //
                         gsMatrix<T>      & quNodes,
                         gsMultiPatch<T> & mp,
                         gsMatrix<> & old_sol,
                         gsG1OptionList optionList)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        std::vector<gsMatrix<unsigned>> actives_temp;

        numActive = 0;
        for (size_t i = 0; i < basis.size(); i++)
        {
            gsMatrix<unsigned> temp_mat;
            basis[i].active_into(md.points.col(0),temp_mat);
            actives_temp.push_back(temp_mat);
            numActive += temp_mat.rows();
            numActive_vec.push_back(numActive);
        }
        numActive += 1; // lambda

        actives.resize(numActive,1);
        index_t ii = 0, jj = 0;
        for (size_t i = 0; i < basis.size(); i++)
        {
            for (index_t j = 0; j < actives_temp[i].rows(); j++)
            {
                actives(ii,0) = actives_temp[i](j,0) + jj;
                ii += 1;
            }
            jj += basis[i].size();
        }
        actives(numActive-1,0) = jj;

        gsVector<> numBasisFunctions(4);
        numBasisFunctions.setZero();
        for (size_t i = 0; i < basis.size(); i++)
            numBasisFunctions[i+1] = numBasisFunctions[i] + basis[i].size();



        // Evaluate basis functions on element
        basis[0].eval_into(md.points,basisData1);
        basis[1].eval_into(md.points,basisData2);
        basis[2].eval_into(md.points,basisData3);

        // ++++++++++++++++++++++++++++++++
        // Create alpha^S (x_old) and beta^S (x_old)
        // ++++++++++++++++++++++++++++++++
        gsGeometry<>::uPtr tilde_temp;

        tilde_temp = basis[0].makeGeometry(old_sol.block(0,0,numBasisFunctions[1],1));
        gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis[1].makeGeometry(old_sol.block(numBasisFunctions[1],0,numBasisFunctions[2] - numBasisFunctions[1],1));
        gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis[2].makeGeometry(old_sol.block(numBasisFunctions[2],0,numBasisFunctions[3] - numBasisFunctions[2],1));
        gsBSpline<T> beta = dynamic_cast<gsBSpline<T> &> (*tilde_temp);


        gsMatrix<> a_L, a_R, b;
        alpha_L.eval_into(md.points,a_L);
        alpha_R.eval_into(md.points,a_R);

        beta.eval_into(md.points,b);
        // ++++++++++++++++++++++++++++++++

        // ++++++++++++++++++++++++++++++++
        // Compute the mapping
        // ++++++++++++++++++++++++++++++++
        gsGeometry<> & FR = mp.patch(0);
        gsGeometry<> & FL = mp.patch(1);

        gsMatrix<> pointV(FR.parDim(), md.points.cols());
        pointV.setZero();
        pointV.row(1) = md.points;

        gsMatrix<> pointU(FL.parDim(), md.points.cols());
        pointU.setZero();
        pointU.row(0) = md.points;

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
        // ++++++++++++++++++++++++++++++++


        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        real_t l = old_sol.bottomRows(1)(0,0);
        gsMatrix<> lam = ones * l;

        // ++++++++++++++++++++++++++++++++
        // Compute rhs
        // ++++++++++++++++++++++++++++++++

        rhsVals1 = - (a_L.cwiseProduct(DuFR_DuFR) + a_R.cwiseProduct(DvFL_DuFR) + b.cwiseProduct(DvFR_DuFR) + lam.cwiseProduct(a_L - ones));
        rhsVals2 = - (a_L.cwiseProduct(DvFL_DuFR) + a_R.cwiseProduct(DvFL_DvFL) + b.cwiseProduct(DvFL_DvFR) + lam.cwiseProduct(a_R - ones));
        rhsVals3 = - (a_L.cwiseProduct(DvFR_DuFR) + a_R.cwiseProduct(DvFL_DvFR) + b.cwiseProduct(DvFR_DvFR));
        rhsVals4 = - ((a_L - ones).cwiseProduct(a_L - ones) + (a_R - ones).cwiseProduct(a_R - ones));
        // ++++++++++++++++++++++++++++++++

        // ++++++++++++++++++++++++++++++++
        // Compute Jacobian matrix
        // ++++++++++++++++++++++++++++++++

        block11 = DuFR_DuFR + lam.cwiseProduct(a_L - ones);
        block12 = DvFL_DuFR;
        block13 = DvFR_DuFR;
        block14 = a_L - ones;

        block21 = DvFL_DuFR;
        block22 = DvFL_DvFL + lam.cwiseProduct(a_R - ones);
        block23 = DvFL_DvFR;
        block24 = a_R - ones;

        block31 = DvFR_DuFR;
        block32 = DvFL_DvFR;
        block33 = DvFR_DvFR;

        block41 = a_L - ones;
        block42 = a_R - ones;

        // ++++++++++++++++++++++++++++++++

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals1.rows() );//multiple right-hand sides

    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals1  = basisData1;
        gsMatrix<T> & basisVals2  = basisData2;
        gsMatrix<T> & basisVals3  = basisData3;

        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            localRhs.block(0, 0, numActive_vec[0], rhsVals1.rows()).noalias() += weight * (basisVals1.col(k) * rhsVals1.col(k).transpose());
            localRhs.block(numActive_vec[0], 0, numActive_vec[1] - numActive_vec[0], rhsVals1.rows()).noalias() += weight * (basisVals2.col(k) * rhsVals2.col(k).transpose());
            localRhs.block(numActive_vec[1], 0, numActive_vec[2] - numActive_vec[1], rhsVals1.rows()).noalias() += weight * (basisVals3.col(k) * rhsVals3.col(k).transpose());
            localRhs.block(numActive_vec[2], 0, 1, rhsVals1.rows()).noalias() += weight * (ones.col(k) * rhsVals4.col(k).transpose());

            // first block row
            localMat.block(0, 0, numActive_vec[0], numActive_vec[0]).noalias() +=
                weight * block11(0, k) * (basisVals1.col(k) * basisVals1.col(k).transpose());

            localMat.block(0, numActive_vec[0], numActive_vec[0], numActive_vec[1] - numActive_vec[0]).noalias() +=
                weight * block12(0, k) * (basisVals1.col(k) * basisVals2.col(k).transpose());

            localMat.block(0, numActive_vec[1], numActive_vec[0], numActive_vec[2] - numActive_vec[1]).noalias() +=
                weight * block13(0, k) * (basisVals1.col(k) * basisVals3.col(k).transpose());

            localMat.block(0, numActive_vec[2], numActive_vec[0], 1).noalias() +=
                weight * block14(0, k) * (basisVals1.col(k) * ones.col(k).transpose());

            // second block row
            localMat.block(numActive_vec[0], 0, numActive_vec[1] - numActive_vec[0], numActive_vec[0]).noalias() +=
                weight * block21(0, k) * (basisVals2.col(k) * basisVals1.col(k).transpose());

            localMat.block(numActive_vec[0], numActive_vec[0], numActive_vec[1] - numActive_vec[0], numActive_vec[1] - numActive_vec[0]).noalias() +=
                weight * block22(0, k) * (basisVals2.col(k) * basisVals2.col(k).transpose());

            localMat.block(numActive_vec[0], numActive_vec[1], numActive_vec[1] - numActive_vec[0], numActive_vec[2] - numActive_vec[1]).noalias() +=
                weight * block23(0, k) * (basisVals2.col(k) * basisVals3.col(k).transpose());

            localMat.block(numActive_vec[0], numActive_vec[2], numActive_vec[1] - numActive_vec[0], 1).noalias() +=
                weight * block24(0, k) * (basisVals2.col(k) * ones.col(k).transpose());

            // third block row
            localMat.block(numActive_vec[1], 0, numActive_vec[2] - numActive_vec[1], numActive_vec[0]).noalias() +=
                weight * block31(0, k) * (basisVals3.col(k) * basisVals1.col(k).transpose());

            localMat.block(numActive_vec[1], numActive_vec[0], numActive_vec[2] - numActive_vec[1], numActive_vec[1] - numActive_vec[0]).noalias() +=
                weight * block32(0, k) * (basisVals3.col(k) * basisVals2.col(k).transpose());

            localMat.block(numActive_vec[1], numActive_vec[1], numActive_vec[2] - numActive_vec[1], numActive_vec[2] - numActive_vec[1]).noalias() +=
                weight * block33(0, k) * (basisVals3.col(k) * basisVals3.col(k).transpose());

            // fourth block row
            localMat.block(numActive_vec[2], 0, 1, numActive_vec[0]).noalias() +=
                weight * block41(0, k) * (ones.col(k) * basisVals1.col(k).transpose());

            localMat.block(numActive_vec[2], numActive_vec[0], 1, numActive_vec[1] - numActive_vec[0]).noalias() +=
                weight * block42(0, k) * (ones.col(k) * basisVals2.col(k).transpose());
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
    gsMatrix<T> basisData1, basisData2, basisData3;
    index_t numActive;
    std::vector<index_t> numActive_vec;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals1, rhsVals2, rhsVals3, rhsVals4;
    gsMatrix<T>  block11, block12, block13, block14,
        block21, block22, block23, block24,
        block31, block32, block33,
        block41, block42;


protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    gsMapData<T> md;


}; // class gsVisitorGluingData

} // namespace gismo
