/** @file gsVisitorResidual.h

    @brief Visitor for element-wise residual for the Poisson equation.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Kleiss et al.
*/


#pragma once

namespace gismo
{


template <class T>
class gsVisitorResidual
{
public:

    gsVisitorResidual(const gsFunction<T> & rhs, const gsField<T> & solution)
    : rhs_ptr(&rhs),
      sol_ptr(&solution)
    { }

    static void initialize(const gsBasis<T> & basis, 
                           gsQuadRule<T> & rule, 
                           int & nonZerosPerCol,
                           unsigned & evFlags ) // replace with geoEval ?
    {
        // Allocate sparse matrix
        nonZerosPerCol = 0;
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
        {
            numQuadNodes[i] = basis.degree(i) + 1;
        }
        
        // Setup Quadrature
        rule.setNodes(numQuadNodes);

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE|NEED_GRAD_TRANSFORM; // NEED_2ND_DERIV
    }


    // Evaluate on element.
    void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                  gsGeometryEvaluator<T> & geoEval,
                  gsMatrix<T>            & quNodes, 
                  gsVector<T>            & quWeights)
    {
        // Evaluate the rhs on the quadrature nodes
        rhs_ptr->eval_into(quNodes, rhsVals);

        // Evaluate the Hessian of the discrete solution // (!) TO DO
        sol_ptr->deriv2_into(quNodes, hessVals);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Initialize local error value
        error = 0;

        // Fetch the element index (assumes all quad. nodes are on the same element)
        elIndex = basis.elementIndex( quNodes.col(0) );
    }

    // assemble on element
    void assemble(gsDomainIterator<T>    & element, 
                  gsGeometryEvaluator<T> & geoEval,
                  const index_t k, T weight       )
    {
        // Multiply quadrature weight by the geometry measure
        weight *= geoEval.measure(k);
        
        // Compute Transformation of the Hessian to the physical domain
        geoEval.transformHessian(k, hessVals, phHessVals);

        //rhsVals.col(k) // for multiple rhs
        error += weight * (  rhsVals(0,k) + phHessVals.trace() )*
                          (  rhsVals(0,k) + phHessVals.trace() );
    }
    
    void localToGlobal(const std::vector< gsDofMapper *>  & mappers,
                       const gsMatrix<T>     & eliminatedDofs,
                       const int patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        rhsMatrix(elIndex, 0 ) = error;
    }

private:

    const gsFunction<T> * rhs_ptr;
    const gsField<T>    * sol_ptr;

private:

    // Gradient values
    gsMatrix<T>      rhsVals;
    gsMatrix<T>      hessVals;
    gsMatrix<T>      phHessVals;

    // Element index
    index_t elIndex;
    
    // Local error
    T error;
};


} // namespace gismo

