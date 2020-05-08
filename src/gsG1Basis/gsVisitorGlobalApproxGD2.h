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
class gsVisitorGlobalApproxGD2
{
public:

    gsVisitorGlobalApproxGD2()
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
                         gsMatrix<T> & old_sol,
                         gsG1OptionList g1OptionList)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        gsMatrix<unsigned> actives_temp;
        basis[0].active_into(md.points.col(0), actives_temp);

        actives.resize(actives_temp.rows()*4,1);
        for (index_t j = 0; j < 4; j++)
            for (index_t i = 0; i < actives_temp.rows(); i++)
                actives(i + j*actives_temp.rows(),0) = actives_temp(i,0) + j*basis[0].size();

        // Evaluate basis functions on element
        basis[0].eval_into(md.points,basisData);

        numActive = actives.rows();

        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        gsMatrix<> lam = ones * g1OptionList.getReal("lambda");

        // ++++++++++++++++++++++++++++++++
        // Create alpha^S (x_old) and beta^S (x_old)
        // ++++++++++++++++++++++++++++++++
        gsGeometry<>::uPtr tilde_temp;

        tilde_temp = basis[0].makeGeometry(old_sol.block(0,0,basis[0].size(),1));
        gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis[0].makeGeometry(old_sol.block(basis[0].size(),0,basis[0].size(),1));
        gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis[0].makeGeometry(old_sol.block(2*basis[0].size(),0,basis[0].size(),1));
        gsBSpline<T> beta_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis[0].makeGeometry(old_sol.block(3*basis[0].size(),0,basis[0].size(),1));
        gsBSpline<T> beta_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);


        gsMatrix<> a_L, a_R, b_L, b_R;
        alpha_L.eval_into(md.points,a_L);
        alpha_R.eval_into(md.points,a_R);

        beta_L.eval_into(md.points,b_L);
        beta_R.eval_into(md.points,b_R);
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


        // ++++++++++++++++++++++++++++++++
        // Compute rhs
        // ++++++++++++++++++++++++++++++++
        rhsVals1 = - a_L.cwiseProduct(lam + DuFR_DuFR + 2*b_R.cwiseProduct(DvFR_DuFR) + b_R.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR));
        rhsVals1 -= a_R.cwiseProduct(DvFL_DuFR + b_R.cwiseProduct(DvFL_DvFR) + b_L.cwiseProduct(DvFR_DuFR) + b_L.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR));
        rhsVals1 += lam;

        rhsVals2 = - a_R.cwiseProduct(lam + DvFL_DvFL + 2*b_L.cwiseProduct(DvFL_DvFR) + b_L.cwiseProduct(b_L).cwiseProduct(DvFR_DvFR));
        rhsVals2 -= a_L.cwiseProduct(DvFL_DuFR + b_L.cwiseProduct(DvFR_DuFR) + b_R.cwiseProduct(DvFL_DvFR) + b_L.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR));
        rhsVals2 += lam;

        rhsVals3 = - b_L.cwiseProduct(lam + a_R.cwiseProduct(a_R).cwiseProduct(DvFR_DvFR));
        rhsVals3 -= a_R.cwiseProduct(a_R.cwiseProduct(DvFL_DvFR) + a_L.cwiseProduct(DvFR_DuFR) + a_L.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR));

        rhsVals4 = - b_R.cwiseProduct(lam + a_L.cwiseProduct(a_L).cwiseProduct(DvFR_DvFR));
        rhsVals4 -= a_L.cwiseProduct(a_R.cwiseProduct(DvFL_DvFR) + a_L.cwiseProduct(DvFR_DuFR) + a_R.cwiseProduct(b_L).cwiseProduct(DvFR_DvFR));
        // ++++++++++++++++++++++++++++++++

        // ++++++++++++++++++++++++++++++++
        // Compute Jacobian matrix
        // ++++++++++++++++++++++++++++++++

        block11 = lam + DuFR_DuFR + 2*b_R.cwiseProduct(DvFR_DuFR) + b_R.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR);
        block12 = DvFL_DuFR + b_R.cwiseProduct(DvFL_DvFR) + b_L.cwiseProduct(DvFR_DuFR) + b_L.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR);
        block13 = a_R.cwiseProduct(DvFR_DuFR + b_R.cwiseProduct(DvFR_DvFR));
        block14 = a_L.cwiseProduct(2*DvFR_DuFR + 2*b_R.cwiseProduct(DvFR_DvFR)) + a_R.cwiseProduct(DvFL_DvFR + b_L.cwiseProduct(DvFR_DvFR));

        block21 = DvFL_DuFR + b_L.cwiseProduct(DvFR_DuFR) + b_R.cwiseProduct(DvFL_DvFR) + b_L.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR);
        block22 = lam + DvFL_DvFL + 2*b_L.cwiseProduct(DvFL_DvFR) + b_L.cwiseProduct(b_L).cwiseProduct(DvFR_DvFR);
        block23 = a_R.cwiseProduct(2*DvFL_DvFR + 2*b_L.cwiseProduct(DvFR_DvFR)) + a_L.cwiseProduct(DvFR_DuFR + b_R.cwiseProduct(DvFR_DvFR));
        block24 = a_L.cwiseProduct(DvFL_DvFR + b_L.cwiseProduct(DvFR_DvFR));

        block31 = a_R.cwiseProduct(DvFR_DuFR + DvFR_DvFR.cwiseProduct(b_R));
        block32 = b_L.cwiseProduct(2*a_R.cwiseProduct(DvFR_DvFR)) + 2*a_R.cwiseProduct(DvFL_DvFR) + a_L.cwiseProduct(DvFR_DuFR) + a_L.cwiseProduct(b_R).cwiseProduct(DvFR_DvFR);
        block33 = lam + a_R.cwiseProduct(a_R).cwiseProduct(DvFR_DvFR);
        block34 = a_R.cwiseProduct(a_L.cwiseProduct(DvFR_DvFR));

        block41 = b_R.cwiseProduct(2*a_L.cwiseProduct(DvFR_DvFR)) + 2*a_L.cwiseProduct(DvFR_DuFR) + a_R.cwiseProduct(DvFL_DvFR) + a_R.cwiseProduct(b_L).cwiseProduct(DvFR_DvFR);;
        block42 = a_L.cwiseProduct(DvFL_DvFR + DvFR_DvFR.cwiseProduct(b_L));
        block43 = a_L.cwiseProduct(a_R.cwiseProduct(DvFR_DvFR));
        block44 = lam + a_L.cwiseProduct(a_L).cwiseProduct(DvFR_DvFR);
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

            localRhs.block(0, 0, numActive / 4, rhsVals1.rows()).noalias() += weight * (basisVals.col(k) * rhsVals1.col(k).transpose());
            localRhs.block(numActive / 4, 0, numActive / 4, rhsVals1.rows()).noalias() += weight * (basisVals.col(k) * rhsVals2.col(k).transpose());
            localRhs.block(2*numActive / 4, 0, numActive / 4, rhsVals1.rows()).noalias() += weight * (basisVals.col(k) * rhsVals3.col(k).transpose());
            localRhs.block(3*numActive / 4, 0, numActive / 4, rhsVals1.rows()).noalias() += weight * (basisVals.col(k) * rhsVals4.col(k).transpose());

            // first block row
            localMat.block(0, 0, numActive / 4, numActive / 4).noalias() +=
                weight * block11(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(0, numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block12(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(0, 2 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block13(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(0, 3 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block14(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            // second block row
            localMat.block(numActive / 4, 0, numActive / 4, numActive / 4).noalias() +=
                weight * block21(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(numActive / 4, numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block22(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(numActive / 4, 2 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block23(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(numActive / 4, 3 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block24(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            // third block row
            localMat.block(2*numActive / 4, 0, numActive / 4, numActive / 4).noalias() +=
                weight * block31(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(2*numActive / 4, numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block32(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(2*numActive / 4, 2 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block33(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(2*numActive / 4, 3 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block34(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            // forth block row
            localMat.block(3*numActive / 4, 0, numActive / 4, numActive / 4).noalias() +=
                weight * block41(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(3*numActive / 4, numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block42(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(3*numActive / 4, 2 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block43(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

            localMat.block(3*numActive / 4, 3 * numActive / 4, numActive / 4, numActive / 4).noalias() +=
                weight * block44(0, k) * (basisVals.col(k) * basisVals.col(k).transpose());

        }

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        gsMatrix<unsigned> actives_temp;

        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives_temp);
        //gsInfo << "actives " << actives << " : " << actives_temp << "\n";
        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives_temp, eliminatedDofs[0], 0, 0);

    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals1, rhsVals2, rhsVals3, rhsVals4;
    gsMatrix<T>  block11, block12, block13, block14,
                 block21, block22, block23, block24,
                 block31, block32, block33, block34,
                 block41, block42, block43, block44;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    gsMapData<T> md;


}; // class gsVisitorGluingData

} // namespace gismo
