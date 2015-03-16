/** @file gsVisitorPoisson.h

    @brief Poisson equation element visitor.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{


template <class T>
class gsVisitorPoisson
{
public:

    /// Constructor with the right hand side function of the Poisson equation
    gsVisitorPoisson(const gsFunction<T> & rhs) : 
    rhs_ptr(&rhs)
    { }

    gsVisitorPoisson() :     rhs_ptr(NULL) { }


    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometryEvaluator<T> & geoEval,
                         gsMatrix<T> const      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);
        numActive = actives.rows();
        
        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 1, basisData);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);// is this generic ??
        
        // Evaluate right-hand side at the geometry points
        rhs_ptr->eval_into( geoEval.values(), rhsVals ); // to do: parametric rhs ?
        
        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive      );
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const typename gsMatrix<T>::Block bVals  = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block bGrads =
            basisData.bottomRows( numActive * element.dim() );

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, physGrad);
            
            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() ) ;
            localMat.noalias() += weight * (physGrad.transpose() * physGrad);
        }
        //gsDebugVar(localRhs.transpose() );
        //gsDebugVar(localMat.asVector().transpose() );
    }
    
    inline void localToGlobal(const gsDofMapper     & mapper,
                              const gsMatrix<T>     & eliminatedDofs,
                              const int patchIndex,
                              gsSparseMatrix<T>     & sysMatrix,
                              gsMatrix<T>           & rhsMatrix )
    {
        //Assert eliminatedDofs.rows() == mapper.boundarySize()

        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        //const int numActive = actives.rows();
        
        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = actives(i);
            if ( mapper.is_free_index(ii) )
            {
                rhsMatrix.row(ii) += localRhs.row(i);
                
                for (index_t j=0; j < numActive; ++j)
                {
                    const int jj = actives(j);
                    if ( mapper.is_free_index(jj) )
                    {
                        // Matrix is symmetric, so we could store only
                        // lower triangular part
                        //if ( jj <= ii ) 
                        sysMatrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                    {
                        rhsMatrix.row(ii).noalias() -= localMat(i, j) * 
                            eliminatedDofs.row( mapper.global_to_bindex(jj) );
                    }
                }
            }
        }
    }

protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;

protected:
    // Basis values
    gsMatrix<T>        basisData;
    gsMatrix<T>        physGrad;
    gsMatrix<unsigned> actives;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
};


} // namespace gismo

