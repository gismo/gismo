
#pragma once

namespace gismo
{


template <class T>
class gsVisitorNeumann
{
public:

    gsVisitorNeumann(const gsFunction<T> & neudata, boundary::side s) : 
    neudata_ptr(&neudata), side(s)
    { }

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule, 
                    unsigned         & evFlags )
    {
        const int dir = direction(side);
        gsVector<int> numQuadNodes ( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_JACOBIAN;
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

        // Evaluate the Neumann data
        neudata_ptr->eval_into(geoEval.values(), neuData);

        // Initialize local matrix/rhs
        localRhs.setZero(numActive, neudata_ptr->targetDim() );
    }

    inline void assemble(gsDomainIterator<T>    & element, 
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, side, unormal);
            
            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k] * unormal.norm();   
            
            localRhs.noalias() += weight * basisData.col(k) * neuData.col(k).transpose() ;
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

        // Push element contribution to the global load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            // convert local dof index to global dof index
            const unsigned jj = actives(j);
            if (mapper.is_free_index(jj))
                rhsMatrix.row(jj) += localRhs.row(j);
        }
    }

private:

    
    // Neumann function
    const gsFunction<T> * neudata_ptr;
    boundary::side side;

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<unsigned> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> neuData;

    // Local matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

};


} // namespace gismo
