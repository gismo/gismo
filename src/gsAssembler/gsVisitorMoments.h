/** @file gsVisitorMoments.h

    @brief Element visitor for moment vector.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/** \brief Visitor for the moment vector of a function
 *
 * Assembles the linear term
 * \f[ (f,v)_\Omega \f]
 *
 */

template <class T, bool paramCoef = false>
class gsVisitorMoments
{
public:
    
    /** \brief Constructor for gsVisitorMoments.
     *
     * \param[in] rhs Given right-hand-side function/source term that, for
     */
    /// Constructor with the right hand side function of the Poisson equation
    gsVisitorMoments(const gsFunction<T> & rhs)
    : rhs_ptr(&rhs)
    { }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T>    & rule,
                    unsigned         & evFlags )
    {
        // Setup Quadrature
        const real_t quA = options.getReal("quA");
        const index_t quB = options.getInt ("quB");
        rule = gsGaussRule<T>(basis, quA, quB);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_VALUE;
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
        basis.evalAllDers_into( quNodes, 0, basisData);
        
        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval.evaluateAt(quNodes);// is this generic ??
        
        // Evaluate right-hand side at the geometry points paramCoef
        // specifies whether the right hand side function should be
        // evaluated in parametric(true) or physical (false)
        rhs_ptr->eval_into( (paramCoef ?  quNodes :  geoEval.values() ), rhsVals );
        
        // Initialize local matrix/rhs
        localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand sides
    }
    
    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];

        // localRhs.noalias() = quWeights.array() *  geoEval.measures().array() *  
        //                      bVals * rhsVals.transpose();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);
            
            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() ) ;
        }
        //gsDebugVar(localRhs.transpose() );
        //gsDebugVar(localMat.asVector().transpose() );
    }
    
    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the right-hand side
        system.pushToRhs(localRhs, actives, 0);
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
            }
        }
    }

protected:
    // Right hand side
    const gsFunction<T> * rhs_ptr;

protected:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<unsigned> actives;
    index_t numActive;

protected:
    // Local values of the right hand side
    gsMatrix<T> rhsVals;

protected:
    // Local matrices
    gsMatrix<T> localRhs;
};


} // namespace gismo

