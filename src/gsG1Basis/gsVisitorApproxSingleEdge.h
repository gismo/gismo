/** @file gsVisitorApproxSingleEdge.h

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
class gsVisitorApproxSingleEdge
{
public:

    gsVisitorApproxSingleEdge()
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
                         gsApproxGluingData<T> & gluingData,
                         gsG1OptionList optionList,
                         index_t m_uv,
                         index_t bfID)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis_target.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        basis_target.eval_into(md.points,basisData);

        numActive = actives.rows();

        gsMatrix<> der_b_plus, b_plus, beta;
        //basis_plus.evalSingle_into(bfID,md.points,b_plus);
        basis_plus.derivSingle_into(bfID,md.points,der_b_plus);


        if (optionList.getInt("gluingData") == gluingData::global) // global
            gluingData.get_beta_tilde().eval_into(md.points,beta);
        else if (optionList.getInt("gluingData") == gluingData::local) // local
            gluingData.get_local_beta_tilde(bfID).eval_into(md.points,beta);

        rhsVals = beta.cwiseProduct(der_b_plus);
        //rhsVals = beta;


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
