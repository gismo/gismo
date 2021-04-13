/** @file gsApproxArgyrisEdgeBasisVisitor.h

    @brief Visitor for the G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

# include <gsArgyris/gsGluingData/gsApproxGluingData.h>


namespace gismo
{
template <short_t d, class T>
class gsApproxArgyrisEdgeBasisVisitor
{
public:

    gsApproxArgyrisEdgeBasisVisitor()
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
    }

    // Evaluate on element.
    inline void evaluate(const gsGeometry<T>    & geo,
                         gsBasis<T>       & basis,
                         gsBasis<T>       & basis_plus,
                         gsBasis<T>       & basis_minus,
                         gsBasis<T>       & basis_geo,
                         gsApproxGluingData<d, T> approxGluingData,
                         gsMatrix<T>            & quNodes,
                         const index_t bfID,
                         const std::string typeBf,
                         const index_t & dir,
                         bool isboundary,
                         const gsOptionList optionList)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 0, basisData);

        numActive = actives.rows();

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo);

        real_t p = bsp_temp.degree();
        real_t tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        // Initialize local matrix/rhs
        if (typeBf == "plus")
        {
            gsMatrix<T> beta,
                    N_0, N_1,
                    N_i_plus,
                    der_N_i_plus;

            if (!isboundary)
                approxGluingData.betaS(dir).eval_into(quNodes.row(dir),beta); // 1-dir == PatchID
            else
                beta.setZero(1, quNodes.cols());

            basis_geo.evalSingle_into(0,quNodes.row(1-dir),N_0); // u
            basis_geo.evalSingle_into(1,quNodes.row(1-dir),N_1); // u

            basis_plus.evalSingle_into(bfID,quNodes.row(dir),N_i_plus); // v
            basis_plus.derivSingle_into(bfID,quNodes.row(dir),der_N_i_plus);

            gsMatrix<T> temp = beta.cwiseProduct(der_N_i_plus);
            rhsVals = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;

            localMat.setZero(numActive, numActive);
            localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides

        } // n_plus
        else if (typeBf == "minus")
        {

            gsMatrix<T> alpha,
                    N_1,
                    N_j_minus;

            if (!isboundary)
                approxGluingData.alphaS(dir).eval_into(quNodes.row(dir),alpha); // 1-dir == PatchID
            else
                alpha.setOnes(1, quNodes.cols());

            basis_geo.evalSingle_into(1,quNodes.row(1-dir),N_1); // u

            basis_minus.evalSingle_into(bfID,quNodes.row(dir),N_j_minus); // v

            if (!isboundary)
                rhsVals = (dir == 0 ? -1 : 1) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;
            else
                rhsVals = (dir == 0 ? -1 : 1) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));

            localMat.setZero(numActive, numActive);
            localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
        } // n_minus

    } // evaluate

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];

        // ( u, v)
        localMat.noalias() =
            bVals * quWeights.asDiagonal() * bVals.transpose();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            // Multiply weight by the geometry measure
            const T weight = quWeights[k];

            // Multiply weight by the geometry measure
            localRhs.noalias() += weight * (bVals.col(k) * rhsVals.col(k).transpose());
            //
        }

    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              gsSparseSystem<T>      & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);
        // Add contributions to the system matrix and right-hand side
        system.push(localMat, localRhs, actives, eliminatedDofs[0], 0, 0);
    }

protected:
    gsMatrix<index_t> actives;
    std::vector<gsMatrix<T> > basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T>  rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

}; // class gsVisitorG1BasisEdge
} // namespace gismo