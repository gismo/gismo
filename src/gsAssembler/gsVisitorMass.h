/** @file gsVisitorMass.h

    @brief Mass visitor for assmbling element mass matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{


/** 
    @brief The visitor computes element mass integrals
*/
template <class T>
class gsVisitorMass
{
public:

    gsVisitorMass(bool _useMapper = true)
    : useMapper(_useMapper)
    { }

    static void initialize(const gsBasis<T> & basis, 
                           gsQuadRule<T> & rule, 
                           unsigned & evFlags )
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_MEASURE;
    }


    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         // todo: add element here for efficiency
                         gsMatrix<T>            & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(quNodes.col(0) , actives);
        const index_t numActive = actives.rows();
 
        // Evaluate basis functions on element
        basis.eval_into(quNodes, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Initialize local matrix/rhs
        localMat.resize(numActive, numActive);
        localMat.setZero();
    }

    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);
        
            localMat.noalias() += weight * ( basisData.col(k) * basisData.col(k).transpose() );
        }
    }
    
    void localToGlobal(const gsDofMapper     & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const int               patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        if ( useMapper ) // Translate local Dofs to global dofs ?
            mapper.localToGlobal(actives, patchIndex, actives);

        const index_t numActive = actives.rows();

        for (index_t i = 0; i < numActive; ++i)
        {
            const int ii = actives(i,0); // N_i
            for (index_t j = 0; j < numActive; ++j)
            {
                const int jj = actives(j,0); // N_j

                // store lower triangular part only
                if ( jj <= ii )
                    sysMatrix.coeffRef(ii, jj) += localMat(i, j); // N_i*N_j
            }
        }
    }

protected:

    bool useMapper;

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<unsigned> actives;

    // Local matrix
    gsMatrix<T> localMat;
};


} // namespace gismo

