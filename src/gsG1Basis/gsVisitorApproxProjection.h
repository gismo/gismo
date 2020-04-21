/** @file gsVisitorApproxProjection.h

    @brief Visitor for the local approximate gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once
# include <gsG1Basis/gsApproxGluingData.h>

namespace gismo
{


template <class T>
class gsVisitorApproxProjection
{
public:

    gsVisitorApproxProjection()
    {
    }

    void initialize(gsBSplineBasis<T> & basis_target, //
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis_target.dim() );
        for (int i = 0; i < basis_target.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis_target.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE ;
    }

    // Evaluate on element.
    inline void evaluate(gsBSplineBasis<T> & basis_target, //
                         gsMatrix<T>      & quNodes,
                         gsMultiPatch<T> & mp,
                         gsBSplineBasis<T> & basis_plus,
                         gsBSplineBasis<T> & basis_minus,
                         gsApproxGluingData<T> & gluingData,
                         gsOptionList optionList)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis_target.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        basis_target.eval_into(md.points,basisData);

        numActive = actives.rows();

        // ++++++++++++++++++++++++++++++++
        // Compute alpha^S and beta^S exact
        // ++++++++++++++++++++++++++++++++
        // alpha^S
        gsMatrix<> uv, ev;

        uv.setZero(2,md.points.cols());
        uv.bottomRows(1) = md.points; // v

        const index_t d = mp.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        const gsGeometry<> & P0 = mp.patch(0); // iFace.second().patch = 0

        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(1);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - D1 * D1 * ev.col(1).transpose() * ev.col(0);

        }

        gsMatrix<> der_b_plus, b_plus, b_minus, beta, alpha;
        basis_plus.evalSingle_into(optionList.getInt("basisID"),md.points,b_plus);
        basis_plus.derivSingle_into(optionList.getInt("basisID"),md.points,der_b_plus);
        basis_minus.evalSingle_into(optionList.getInt("basisID"),md.points,b_minus);

        if (optionList.getInt("gluingData") == 0) // global
            gluingData.get_beta_tilde().eval_into(md.points,beta);
        else if (optionList.getInt("gluingData") == 1) // local
            gluingData.get_local_beta_tilde(optionList.getInt("basisID")).eval_into(md.points,beta);

        if (optionList.getInt("gluingData") == 0) // global
            gluingData.get_alpha_tilde().eval_into(md.points,alpha);
        else if (optionList.getInt("gluingData") == 1) // local
            gluingData.get_local_alpha_tilde(optionList.getInt("basisID")).eval_into(md.points,alpha);

        //rhsVals = b_plus + beta.cwiseProduct(der_b_plus);
        rhsVals = b_plus + beta.cwiseProduct(der_b_plus);


        // ++++++++++++++++++++++++++++++++
        // ================================
        // ++++++++++++++++++++++++++++++++

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides

    }

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;

        // ( u, v)
        localMat.noalias() =
            basisData * quWeights.asDiagonal() * basisData.transpose();


        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k];
            localRhs.noalias() += weight * (basisVals.col(k) * rhsVals.col(k).transpose());
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

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T>  localRhs;

    gsMapData<T> md;


}; // class gsVisitorApproxProjection

} // namespace gismo
