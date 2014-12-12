
#pragma once

namespace gismo
{


template <class T>
class gsVisitorNitsche
{
public:

    gsVisitorNitsche(const gsFunction<T> & dirdata, T _penalty, boxSide s) : 
    dirdata_ptr(&dirdata),penalty(_penalty), side(s)
    { }

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T> & rule,
                    unsigned & evFlags  )
    {
        const int dir = side.direction();
        gsVector<int> numQuadNodes ( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
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

        // Evaluate basis values and derivatives on element
        basis.evalAllDers_into( quNodes, 1, basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Dirichlet data
        dirdata_ptr->eval_into(geoEval.values(), dirData);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, dirdata_ptr->targetDim() );
    }

    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        const index_t numActive = actives.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

        const typename gsMatrix<T>::Block bVals  = basisData.block(0,k,numActive,1);
        const typename gsMatrix<T>::Block bGrads = 
            basisData.bottomRows( numActive * element.dim() );

        // Compute the outer normal vector on the side
        geoEval.outerNormal(k, side, unormal);

        // Multiply quadrature weight by the geometry measure
        const T weight = quWeights[k] *unormal.norm();   

        // Compute the unit normal vector 
        unormal.normalize();
        
        // Compute physical gradients at k as a Dim x NumActive matrix
        geoEval.transformGradients(k, bGrads, pGrads);
        
        // Get penalty parameter
        const T mu = penalty / element.getCellSize();

        // Sum up quadrature point evaluations
        localRhs.noalias() += weight * (( pGrads.transpose() * unormal - mu * bVals )
                                        * dirData.col(k).transpose() );

        localMat.noalias() += weight * ( bVals * unormal.transpose() * pGrads
                           +  (bVals * unormal.transpose() * pGrads).transpose()
                           -  mu * bVals * bVals.transpose() );
        }
    }
    
    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>     & eliminatedDofs,
                       const int patchIndex,
                       gsSparseMatrix<T>     & sysMatrix,
                       gsMatrix<T>           & rhsMatrix )
    {
        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            rhsMatrix.row(jj) -= localRhs.row(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = actives(i);
                if ( jj <= ii ) // assuming symmetric problem (!) probably we must not.
                    sysMatrix( ii, jj ) -= localMat(i,j);
            }
        }

    }

private:
    // Dirichlet function
    const gsFunction<T> * dirdata_ptr;

    // Penalty constant
    T penalty;

    // Side
    boxSide side;

private:
    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<T>      pGrads;
    gsMatrix<unsigned> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> dirData;

    // Local  matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

};


} // namespace gismo
