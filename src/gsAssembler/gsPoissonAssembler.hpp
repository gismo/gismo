/** @file gsPoissonAssembler.hpp

    @brief Provides assembler implementation for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/

#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals and load vector


namespace gismo
{

template<class T>
void gsPoissonAssembler<T>::applyOptions()
{
    // Would be efficient not to recreate the mappers every time.
    // For this, a , "bool force" argument can be added, defaulting to false
    // Commented for now

    //if ( m_dirStrategy != options.dirStrategy ||
     //    m_intStrategy != options.intStrategy )
    //{
        const bool conforming = ( m_options.intStrategy == iFace::glue );
        m_dofMappers.resize(1);
        if ( m_options.dirStrategy == dirichlet::elimination )
            m_bases.front().getMapper(conforming, m_bConditions, m_dofMappers.front() );
        else
            m_bases.front().getMapper(conforming, m_dofMappers.front() );
        
        m_dofs = m_dofMappers.front().freeSize();

   // }
}

template<class T>
void gsPoissonAssembler<T>::setOptions(const gsAssemblerOptions& options)
{
    m_options = options;

    applyOptions();
}

template<class T>
void gsPoissonAssembler<T>::assembleNitsche()
{
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); ++it )
    {
        gsVisitorNitsche<T> nitsche(*it->function(), penalty(it->patch()), it->side());
            
        // Note: it->unknown()
        this->apply(nitsche, it->patch(), it->side() );
    }
}

template<class T>
void gsPoissonAssembler<T>::assembleNeumann()
{
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.neumannBegin();
          it != m_bConditions.neumannEnd(); ++it )
    {
        gsVisitorNeumann<T> neumann(*it->function(), it->side());

        // Note: it->unknown()
        this->apply(neumann, it->patch(), it->side() );
    }
}

template<class T>
void gsPoissonAssembler<T>::penalizeDirichlet()
{
    static const T PP = 1e9; // magic number
    
    // Store mapper
    gsDofMapper old = m_dofMappers[0];

    // Recompute dof mapper
    const bool conforming = ( m_options.intStrategy == iFace::conforming );
    m_bases.front().getMapper(conforming, m_bConditions, 0, m_dofMappers[0]);

    // Compute dirichlet values
    computeDirichletDofs();

    // apply BCs
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); ++it )
    {            
        const gsBasis<T> & basis = (m_bases[0])[it->patch()];

        gsMatrix<unsigned> bnd = safe( basis.boundary(it->side() ) ); 
        for (index_t k=0; k!= bnd.size(); ++k)
        {
            const index_t ii = old            .index ( bnd(k) , it->patch() );
            const index_t bb = m_dofMappers[0].bindex( bnd(k) , it->patch() );

            m_matrix(ii,ii) = PP;
            m_rhs.row(ii)   = PP * m_ddof.row(bb);
        }
    }

    // restore mapper
    m_dofMappers[0] = old;
}

template<class T>
void gsPoissonAssembler<T>::assemble()
{

    if ( m_options.dirStrategy == dirichlet::elimination ||
         m_options.dirStrategy == dirichlet::penalize     )
    {
        computeDirichletDofs();
    }

    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
        nonZerosPerCol *= 3 * m_bases.front().maxDegree(i) + 1; // need more for DG !

    m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
        
    // Resize the load vector
    m_rhs.setZero(m_dofs, m_rhsFun->targetDim() );


    // Assemble volume stiffness and load vector integrals
    gsVisitorPoisson<T> poisson(*m_rhsFun);
    for (unsigned np=0; np < m_patches.nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(poisson, np);
    }

    // If requested, force Dirichlet boundary conditions by Nitsche's method
    if ( m_options.dirStrategy == dirichlet::nitsche )
        assembleNitsche();

    // If requested, force Dirichlet boundary conditions by diagonal penalization
    if ( m_options.dirStrategy == dirichlet::penalize )
        penalizeDirichlet();

    // Enforce Neumann boundary conditions
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T>
void gsPoissonAssembler<T>::computeDirichletDofs()
{
    switch (m_options.dirValues)
    {
    case dirichlet::homogeneous:
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros
        m_ddof.setZero( m_dofMappers[0].boundarySize(), m_rhsFun->targetDim() );
        break;
    case dirichlet::interpolation:        
        computeDirichletDofsIntpl();
        break;
    case dirichlet::l2Projection:
        computeDirichletDofsL2Proj();
        break;
    case dirichlet::user:
        // 
        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values.");
    }

    // Corner values
    const gsDofMapper & mapper = m_dofMappers.front();
    for ( typename gsBoundaryConditions<T>::const_citerator
              it = m_bConditions.cornerBegin();
          it != m_bConditions.cornerEnd(); ++it )
    {
        const int i  = m_bases[it->unknown][it->patch].functionAtCorner(it->corner);
        const int ii = mapper.bindex( i , it->patch );
        m_ddof.row(ii).setConstant(it->value);
    }
    
}    
    

// SKleiss: Note that this implementation is not useable for (T)HB-Splines!
//
// 1. Computation of the Dirichlet values explicitly uses gsPointGrid(rr)
// Computing a grid of evaluation points does not make sense for any locally
// refined basis.
// Also, "component(i)" is used.
// I'm afraid this makes sense ONLY FOR TENSOR-PRODUCT bases.
//
// 2. As of now (16.May 2014), the boundaryBasis of (T)HB-spline basis is not
// implemented, as far as I know.
//
// 3. gsInterpolate uses the anchors of the boundary basis.
// With truncated hierarchical B-splines, the use of classical anchors does
// not work, because functions might be truncated to zero at these points.
template<class T> // 
void gsPoissonAssembler<T>::computeDirichletDofsIntpl()
{
    const gsDofMapper & mapper = m_dofMappers.front();
    
    m_ddof.resize( mapper.boundarySize(), m_rhsFun->targetDim() ); //--mrhs
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); ++it )
    {
        const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk])[k];

        // Get dofs on this boundary
        gsMatrix<unsigned> * boundary = basis.boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary->size(); ++i)
            {
                const int ii= mapper.bindex( (*boundary)(i) , k );
                m_ddof.row(ii).setZero();
            }
            delete boundary;
            continue;
        }

        // Get the side information
        int dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( this->patches().parDim() );

        for ( int i=0; i < this->patches().parDim(); ++i)
        {
            if ( i==dir )
            {
                gsVector<T> b(1); 
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back( basis.component(i).anchors()->transpose() );
            }
        }

        GISMO_ASSERT(it->function()->targetDim() == m_rhsFun->targetDim(), 
                     "Given Dirichlet boundary function does not match problem dimension.\n");

        // Compute dirichlet values
        gsMatrix<T> fpts;
        if ( it->parametric() )
            fpts = it->function()->eval( gsPointGrid( rr ) );
        else
            fpts = it->function()->eval( m_patches[it->patch()].eval(  gsPointGrid( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = basis.boundaryBasis(it->side());
        gsGeometry<T> * geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t l=0; l!= boundary->size(); ++l)
        {
            const int ii= mapper.bindex( (*boundary)(l) , it->patch() );

            m_ddof.row(ii) = dVals.row(l);
        }

        delete h;
        delete geo;
        delete boundary;
    }
}

template<class T>
void gsPoissonAssembler<T>::computeDirichletDofsL2Proj()
{
    const gsDofMapper & mapper = m_dofMappers.front();

    m_ddof.resize( mapper.boundarySize(), m_rhsFun->targetDim() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>       globProjRhs;   
    globProjRhs.setZero( mapper.boundarySize(), m_rhsFun->targetDim() );

    // Temporaries
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals;

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_bConditions.dirichletBegin();
          iter != m_bConditions.dirichletEnd(); ++iter )
    {
        const int unk = iter->unknown();
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        typename gsGeometry<T>::Evaluator geoEval( m_patches[patchIdx].evaluator(NEED_MEASURE));

        // Set up quadrature. The number of integration points in the direction
        // that is NOT along the element has to be set to 1.
        gsVector<index_t> numQuNodes( basis.dim() );
        for( int i=0; i < basis.dim(); i++)
            numQuNodes[i] = (basis.degree(i)+1);
        numQuNodes[ iter->side().direction()] = 1;

        gsGaussRule<T> bdQuRule(numQuNodes);

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                          quNodes, quWeights);

            geoEval->evaluateAt( quNodes );

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_patches[patchIdx].eval( quNodes ) );

            basis.eval_into( quNodes, basisVals);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(quNodes.col(0), globIdxAct );
            mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index( globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            // Do the actual assembly:
            for( index_t k=0; k < quNodes.cols(); k++ )
            {
                const T weight_k = quWeights[k] * geoEval->measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const unsigned ii = mapper.global_to_bindex( globIdxAct( i ));

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const unsigned j = eltBdryFcts[j0];
                        const unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j

                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    gsSparseSolver<>::CGDiagonal solver;
    m_ddof = solver.compute( globProjMat ).solve ( globProjRhs );
    
} // computeDirichletDofsL2Proj


template<class T>
gsField<T> *  gsPoissonAssembler<T>::constructSolution(const gsMatrix<T>& solVector) const
//gsField<T> & result ) const
{
    //For IETI, this assertion cannot be passed, because the the global system is never assembled.
    //probably add ana accesor to m_rhs to sets its size accordingly.
    //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    const gsDofMapper & mapper = m_dofMappers.front();

    GISMO_ASSERT(solVector.rows() == mapper.freeSize(), "Something went wrong, solution vector is not OK.");

    std::vector<gsFunction<T> * > sols ;

    const index_t dim = m_rhsFun->targetDim();
    gsMatrix<T> coeffs;
    
    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {    
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[0][p].size();
        coeffs.resize( sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof.row( mapper.bindex(i, p) );
            }
        }
        
        sols.push_back( m_bases[0][p].makeGeometry( give(coeffs) ) );
    }

    //result = gsField<T>(m_patches, sols);
    return new gsField<T>(m_patches, sols);
}


template<class T>
void gsPoissonAssembler<T>::constructSolution(const gsMatrix<T>& solVector, 
                                              gsMultiPatch<T>& result) const
{
    // we might need to get a result even without having the system ..
    //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    const gsDofMapper & mapper = m_dofMappers.front();

    GISMO_ASSERT(solVector.rows() == mapper.freeSize(), "Something went wrong, solution vector is not OK.");

    result.clear(); // result is cleared first
    
    const index_t dim = m_rhsFun->targetDim();
    
    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {    
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[0][p].size();
        gsMatrix<T> coeffs( sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof.row( mapper.bindex(i, p) );
            }
        }
        
        result.addPatch( m_bases[0][p].makeGeometry( give(coeffs) ) );
    }

    // AM: result topology ?
}

}// namespace gismo
