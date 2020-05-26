/** @file gsVisitorLocalApproxGD.h

    @brief Visitor for the local approximate gluing data.

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
class gsVisitorLocalApproxGD
{
public:

    gsVisitorLocalApproxGD()
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
        md.flags = NEED_MEASURE ;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T>       & basis, //
                         gsMatrix<T>      & quNodes,
                         gsMultiPatch<T> & mp,
                         index_t & dir,
                         index_t & gamma,
                         std::string & alphaBeta,
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
        gsMatrix<> uv, ev;
        if (alphaBeta == "alpha")
        {
            // alpha^S
            if (dir==1)
            {
                uv.setZero(2,md.points.cols());
                uv.bottomRows(1) = md.points; // v
            }
            else if (dir==0)
            {
                uv.setZero(2,md.points.cols());
                uv.topRows(1) = md.points; // u
            }

            // ======== Determine bar{alpha^(L)} == Patch 0 ========
            const gsGeometry<> & P0 = mp.patch(0); // iFace.second().patch = 0

            for (index_t i = 0; i < uv.cols(); i++)
            {
                P0.jacobian_into(uv.col(i), ev);
                uv(0, i) = gamma * ev.determinant();
            }
            if (isBoundary)
                uv.setOnes();
            rhsVals = uv.row(0);
        }
        else if (alphaBeta == "beta")
        {
            // beta^S
            if (dir==1)
            {
                uv.setZero(2,md.points.cols());
                uv.bottomRows(1) = md.points; // v
            }
            else if (dir==0)
            {
                uv.setZero(2,md.points.cols());
                uv.topRows(1) = md.points; // u
            }

            const index_t d = mp.parDim();
            gsVector<> D0(d);

            // ======== Determine bar{beta}^L ========
            const gsGeometry<> & P0 = mp.patch(0); // iFace.second().patch = 0

            for(index_t i = 0; i < uv.cols(); i++)
            {
                P0.jacobian_into(uv.col(i),ev);
                D0 = ev.col(dir);
                real_t D1 = 1/ D0.norm();
                uv(0,i) = - gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);

            }
            if (isBoundary)
                uv.setZero();
            rhsVals = uv.row(0);
        }
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


}; // class gsVisitorApproxGD

} // namespace gismo
